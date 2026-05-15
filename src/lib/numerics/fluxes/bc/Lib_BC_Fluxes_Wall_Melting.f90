module MOSE_Lib_BC_Fluxes_Wall_Melting
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_Wall_Melting ( Im, Jm, Km, Fm, Blk, T_wall, q_rad, Ovar )
    use MOSE_Lib_Fluid
    use MOSE_Lib_GSI
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: T_wall, q_rad
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    ! Local
    integer :: s, modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, Ig, Jg, Kg
    real(R8) :: Normal(3), Area, Dist, dl, M(3,3), Prim(nprim), rho, Rgas, T, Gradient(nprim,3)
    real(R8) :: rho_wall, Prim_Wall(nprim), mil, kl, q_conv, omega(nsc), mdot
    real(R8) :: Stress(3), Blowing(3), Sum_Omegai_Hwi, Flux(nprim)


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

    ! Difference in Across-face direction: symmetry; tangential direction: null
    Gradient(np,Dir) = ( T - T_wall ) * modfm3                  ! Temperature

    Gradient(np,:) = matmul ( Gradient(np,:), M )               ! Temperature gradient
    rho_wall = Prim(np) / ( Rgas * T_wall )                     ! Boundary layer: p(wall) = p(cell), R(wall) = R(cell) and ci(wall) = ci(cell), but not roi
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho             ! Partial densities at wall

    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl ) ! transport

    if (model==2) then ! turbulence variables
      dl = blk % yn(im,jm,km)
      call RANS_Set_Wall_Values( mil, Prim_Wall(nt:nprim), dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
    end if

    q_conv = kl * dot_product( Gradient(np,:), Normal)          ! Convective heat flux

    call GSI_Melting ( q_conv + q_rad, T_wall, omega, mdot )    ! Gas-surface interaction for paraffin mass flux

    Blowing = mdot / rho_wall * Normal                          ! Blowing velocity
    Gradient(nu:nw,Dir) = ( Prim(nu:nw) - Blowing ) * modfm3    ! Velocity gradient due to blowing in direction normal to the face

    Gradient(nu:nw,:) = matmul ( Gradient(nu:nw,:), M )         ! Velocity gradient transformation

    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, 0d0, Prim(nt:) )

    Sum_Omegai_Hwi = 0.d0
    do s = 1, nsc
      Sum_Omegai_Hwi = Sum_Omegai_Hwi + omega(s) * f_tabT( T_wall,s,h_tab )
    end do

    ! Fluxes
    Flux = 0.0d0
    Flux(1:nsc) = - omega * Area
    Flux(nu:nw) = ( Stress - mdot * Blowing ) * Area
    Flux(np) = Area * ( q_conv + q_rad - Sum_Omegai_Hwi - mdot * 0.5d0 * sum( Blowing**2 ) + &
                        dot_product( Stress, Blowing ) ) ! check sign q_rad!!!

    if (model==2) then
      Gradient(nt:nprim,:) = matmul ( Gradient(nt:nprim,:), M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, area=area, &
                                 normal=normal, dist=dist )
      ! Ghost-cell extrapolation of RANS variables
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), qw=q_conv+q_rad, &
                                                    mdot=mdot, y=blk%dl(im,jm,km)%c(dir)*0.5d0, &
                                                    Tw=T_wall, rhow=rho_wall, mu=mil, exit_array=Ovar)

  end subroutine BC_Wall_Melting

end module MOSE_Lib_BC_Fluxes_Wall_Melting