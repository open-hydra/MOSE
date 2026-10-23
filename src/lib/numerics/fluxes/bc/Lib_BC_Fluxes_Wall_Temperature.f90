module MOSE_Lib_BC_Fluxes_Wall_Temperature
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties, Wall_State_Temperature

  implicit none
  public

contains

  !> Wall of prescribed temperature T_wall: the viscous flux of the face. Its ghost values, turbulence variables included,
  !> are those of Fill_Ghost_Cell (Ghost_Wall in MOSE_Lib_Ghost), the one set of ghost cells of the fluxes and of the source
  !> terms (an SRM grain, 502, falls back on this wall when it does not burn: Ghost_SRM_Turbulence takes its wall values).
  subroutine BC_Wall_Temperature ( Im, Jm, Km, Fm, Blk, T_wall, Ovar, w_wall, k_rough, rec, w_grad )
    use MOSE_Config_Types_m, only: obj_rans
    use MOSE_Lib_Diffusive, only: Inverse3
    use MOSE_Lib_Fluid
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: T_wall
    type(MOSE_block_type), intent(inout), target :: Blk
    real(R8), optional, dimension(8), intent(out) :: Ovar
    real(R8), optional, dimension(3), intent(in) :: w_wall
    real(R8), optional, intent(in) :: k_rough   ! sand-grain roughness of the face (absent = smooth)
    ! the record of the face: its wall stencil (diffusive-metric = stencil)
    type(MOSE_bc_type), optional, intent(in), target :: rec
    real(R8), optional, dimension(3,3), intent(in) :: w_grad   ! with w_wall: its gradient, d w_wall(i) / d x(j)
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, q
    real(R8) :: Normal(3), Area, Dist, Mi(3,3), Prim(nprim), rho, Rgas, T, rho_wall, dl
    real(R8), pointer, contiguous :: M(:,:)   ! the metric of the face, read in place
    real(R8) :: Prim_wall(nprim), mil, kl, Gradient(nprim,3), Stress(3), Flux(nprim)
    integer  :: vgrad_i
    real(R8) :: gvec1, gvec2, gvec3, VelGrad(3,3)
    real(R8) :: kr, mit, kt
    integer  :: tb(2)
    logical  :: wst

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    kr = 0d0
    if ( present(k_rough) ) kr = k_rough

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    Dist = 1d-20 ! since the flux is computed at the wall
    M => Blk % M (Im,Jm,Km) % c
    ! diffusive-metric = stencil: the metric of the wall face, with the differences of the wall values along the wall
    ! below (Wall_Face_Stencils in MOSE_Lib_Stencil_Metrics); otherwise the metric of the cell and no tangential difference
    ! but those of the velocity of a wall that moves (below)
    wst = .false.
    if ( present(rec) ) wst = rec % wstencil
    if ( wst ) M => rec % wM

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    ! Initialization
    Gradient = 0d0

    ! Boundary layer: p(wall) = p(cell). Also, non catalytic wall, so R(wall) = R(cell) and ci(wall) = ci(cell), but not roi.
    ! Laminar viscosity/conductivity at interface computed at Twall. The same wall state as the ghost cells of the wall
    ! (Ghost_Wall in MOSE_Lib_Ghost).
    call Wall_State_Temperature ( Prim, rho, Rgas, T_wall, rho_wall, Prim_wall(1:nsc), mil, kl )

    ! Difference in Across-face direction: symmetry; tangential direction: null (below: with the wall stencil the
    ! differences of the wall temperature, on a wall that moves those of its velocity)
    Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3                  ! Velocity
    if ( present(w_wall) ) then
      Gradient(nu:nw,Dir) = ( Prim(nu:nw) - w_wall ) * modfm3  ! Velocity relative to wall
    else
      Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3               ! Standard no-slip (w_wall=0)
    end if
    Gradient(np,Dir) = ( T - T_wall ) * modfm3                  ! Temperature
    if ( wst ) then
      select case (Dir)
        case (1) ; tb = [2, 3]
        case (2) ; tb = [1, 3]
        case default ; tb = [1, 2]
      end select
      Gradient(np,tb(1)) = rec % wdT(1) ; Gradient(np,tb(2)) = rec % wdT(2)
    end if
    ! A wall that moves with a velocity linear in x (a stationary face of a rotating frame): the differences of the wall
    ! velocity along the position differences of the metric, the columns of M^-1 (the stencil of the wall face, or the
    ! edges of the cell), from its gradient
    if ( present(w_grad) ) then
      Mi = Inverse3 ( M )                                       ! columns: the position differences of M
      do q = 1, 3
        if ( q /= Dir ) Gradient(nu:nw,q) = matmul ( w_grad, Mi(:,q) )
      enddo
    end if

    if (model==2)  then
      dl = blk % yn(im,jm,km)
      call RANS_Set_Wall_Values( mil, Prim(nt:nprim) * rho_wall / rho, Prim_Wall(nt:nprim), dl, kr )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3 
    end if

    ! Transformation

    !> Computational -> physical transform.  Explicit loop, not
    !> `Gradient = matmul(Gradient,M)`: Gradient aliases itself there, so the
    !> compiler would build the result in a temporary on every face.
    do vgrad_i = 1, nprim
      gvec1 = Gradient(vgrad_i,1) ; gvec2 = Gradient(vgrad_i,2) ; gvec3 = Gradient(vgrad_i,3)
      Gradient(vgrad_i,1) = gvec1*M(1,1) + gvec2*M(2,1) + gvec3*M(3,1)
      Gradient(vgrad_i,2) = gvec1*M(1,2) + gvec2*M(2,2) + gvec3*M(3,2)
      Gradient(vgrad_i,3) = gvec1*M(1,3) + gvec2*M(2,3) + gvec3*M(3,3)
    end do
    !> Contiguous copy of the velocity-gradient rows: passing the strided
    !> section directly makes each consumer build its own temporary.
    VelGrad(1,1)=Gradient(nu,1); VelGrad(1,2)=Gradient(nu,2); VelGrad(1,3)=Gradient(nu,3)
    VelGrad(2,1)=Gradient(nv,1); VelGrad(2,2)=Gradient(nv,2); VelGrad(2,3)=Gradient(nv,3)
    VelGrad(3,1)=Gradient(nw,1); VelGrad(3,2)=Gradient(nw,2); VelGrad(3,3)=Gradient(nw,3)

    ! Eddy viscosity at the wall. A smooth wall zeroes the turbulence variables
    ! there, so it is zero; a rough wall leaves nit_wall > 0 (SA-rough).
    mit = 0d0
    kt  = 0d0
    if ( model==2 .and. kr > 0d0 ) then
      call Eddy_Viscosity ( mut=mit, rans_variables=Prim_Wall(nt:nprim), mul=mil, rho=rho_wall, &
                            vel_gradient=VelGrad, walldist=0d0, k_rough=kr )
      kt = mit * f_cp ( Prim_wall(1:nsc), T_wall, rho_wall ) / obj_rans%Prt
    end if

    ! Stress vector
    Stress = Stress_Vector ( VelGrad, Normal, mil, mit, Prim(nt:) )

    ! Fluxes
    Flux = 0.0d0 
    Flux(nu:nw) = Stress * Area
    Flux(np) = Area * ( kl + kt ) * dot_product ( Gradient(np,:), Normal )
    ! A wall that moves: the work of the wall stress on the wall velocity, (tau . n) . w_wall, as a face between cells
    ! carries (tau . n) . u of its velocity
    if ( present(w_wall) ) Flux(np) = Flux(np) + Area * dot_product ( Stress, w_wall )

    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=VelGrad, &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=dist )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), into=-modfm2, Tw=T_wall, &
                                                    rhow=rho_wall, mu=mil,                          &
                                                    y=blk%dl(im,jm,km)%c(dir)*0.5d0,                &
                                                    qw=( kl + kt ) * dot_product ( Gradient(np,:), Normal ), &
                                                    exit_array=Ovar)

  end subroutine BC_Wall_Temperature

end module MOSE_Lib_BC_Fluxes_Wall_Temperature