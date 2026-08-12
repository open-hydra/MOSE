module MOSE_Lib_BC_Fluxes_Wall_Temperature
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_Wall_Temperature ( Im, Jm, Km, Fm, Blk, T_wall, Ovar, w_wall )
    use MOSE_Lib_Fluid
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: T_wall
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(out) :: Ovar
    real(R8), optional, dimension(3), intent(in) :: w_wall
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, ig, jg, kg
    real(R8) :: Normal(3), Area, Dist, M(3,3), Prim(nprim), rho, Rgas, T, rho_wall, dl
    real(R8) :: Prim_wall(nprim), mil, kl, Gradient(nprim,3), Stress(3), Flux(nprim)
    integer  :: vgrad_i
    real(R8) :: gvec1, gvec2, gvec3, VelGrad(3,3)

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    Dist = 1d-20 ! since the flux is computed at the wall
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    ! Initialization
    Gradient = 0d0

    ! Boundary layer: p(wall) = p(cell). Also, non catalytic wall, so R(wall) = R(cell) and ci(wall) = ci(cell), but not roi
    rho_wall = Prim(np) / ( Rgas * T_wall )
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho

    ! Laminar viscosity/conductivity at interface computed at Twall
    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )

    ! Difference in Across-face direction: symmetry; tangential direction: null
    Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3                  ! Velocity
    if ( present(w_wall) ) then
      Gradient(nu:nw,Dir) = ( Prim(nu:nw) - w_wall ) * modfm3  ! Velocity relative to wall
    else
      Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3               ! Standard no-slip (w_wall=0)
    end if
    Gradient(np,Dir) = ( T - T_wall ) * modfm3                  ! Temperature

    if (model==2)  then
      dl = blk % yn(im,jm,km)
      call RANS_Set_Wall_Values( mil, Prim_Wall(nt:nprim), dl )
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

    ! Stress vector
    Stress = Stress_Vector ( VelGrad, Normal, mil, 0d0, Prim(nt:) )

    ! Fluxes
    Flux = 0.0d0 
    Flux(nu:nw) = Stress * Area
    Flux(np) = Area * kl * dot_product ( Gradient(np,:), Normal )

    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=VelGrad, &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=dist )
      ! Ghost-cell extrapolation of RANS variables
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), Tw=T_wall,          &
                                                    rhow=rho_wall, mu=mil,                          &
                                                    y=blk%dl(im,jm,km)%c(dir)*0.5d0,                &
                                                    qw=kl * dot_product ( Gradient(np,:), Normal ), &
                                                    exit_array=Ovar)

  end subroutine BC_Wall_Temperature

end module MOSE_Lib_BC_Fluxes_Wall_Temperature