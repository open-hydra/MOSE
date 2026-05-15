module MOSE_Lib_BC_Fluxes_SRM
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_SRM ( Im, Jm, Km, Fm, Blk, aCoeff, n, pRef, Taf, haf, BC_ci, Ovar )
    use MOSE_Lib_Fluid
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic
    use MOSE_Global_m
    implicit none

    integer,  intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: aCoeff, n, pRef, Taf, haf, BC_ci(nsc)
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, Ig, Jg, Kg
    real(R8) :: Normal(3), Area, Dist, dl, M(3,3), Prim(nprim), rho, Rgas, Gradient(nprim,3)
    real(R8) :: rho_wall, Prim_Wall(nprim), mil, kl, mdot, hf, SFlocal, T
    real(R8) :: Stress(3), Blowing(3), Flux(nprim)
    real(R8) :: q_conv, q_rad

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    Dist = Blk % yn(Im,Jm,Km)
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    rho_wall = Prim(np) / (Rgas*Taf)

    mdot = aCoeff*(Prim(np)/pRef)**n * (-modFm2)
    Blowing = (mdot/rho_wall)*Normal ! Blowing velocity

    q_conv = 0.d0
    Stress = 0.d0
    Gradient = 0d0
    if (model>0) then
      Prim_wall(1:nsc) = BC_ci * rho_wall ! Partial densities at wall
      T = Prim(np) / ( rho * Rgas )
      call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, Taf, mil, kl )

      Gradient(np,Dir) = ( T - Taf ) * modfm3
      Gradient(np,:) = matmul ( Gradient(np,:), M )
      q_conv = kl * dot_product(Gradient(np,:), normal)

      Gradient(nu:nw,Dir) = ( Prim(nu:nw) - Blowing ) * modfm3    ! Velocity gradient due to blowing in direction normal to the face
      Gradient(nu:nw,:) = matmul ( Gradient(nu:nw,:), M )         ! Velocity gradient transformation
      Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, 0d0, Prim(nt:) ) 
    endif

    ! Fluxes
    Flux(1:nsc) = - mdot * Area * BC_ci
    Flux(nu:nw) = ( Stress - mdot * Blowing ) * Area
    Flux(np) = Area * (q_conv - mdot*( haf + 0.5d0*sum( Blowing**2 ) ) &
                        + dot_product( Stress, Blowing ) ) ! check sign q_rad!!!
    if (model==2) then
      dl = blk % yn(im,jm,km)
      call RANS_Set_Blowing_Wall ( rho=rho_wall, mil=mil, &
                                   rans_variables=Prim_Wall(nt:nprim), &
                                   tau=stress, mdot=mdot, dist=dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
      Gradient(nt:nprim,:) = matmul ( Gradient(nt:nprim,:), M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=1d-20 )
      ! Ghost-cell extrapolation of RANS variables
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux
    
    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), qw=q_conv, &
                                                    mdot=mdot, y=blk%dl(im,jm,km)%c(dir)*0.5d0, &
                                                    Tw=Taf, rhow=rho_wall, mu=mil, exit_array=Ovar)
    
  end subroutine BC_SRM

end module MOSE_Lib_BC_Fluxes_SRM