module MOSE_Lib_Diffusive
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Diffusive_Flux, Compute_Diffusive_Flux

contains

  subroutine Diffusive_Flux ( normal, area, waldis1, waldis2, Prim1, Prim2, Prim3, Prim4, Prim5, &
                              Prim6, Prim7, Prim8, Prim9, Prim10, M1, M2, Res1, Res2, a, b, c, &
                              Sc, Sct, Prt, soot_enabled )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in)  :: a, b, c
    logical, intent(in)  :: soot_enabled
    real(R8), intent(in) :: Sc, Sct, Prt
    real(R8), intent(in) :: normal(3), area, waldis1, waldis2
    real(R8), intent(in), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4, Prim5, Prim6
    real(R8), intent(in), dimension(nprim) :: Prim7, Prim8, Prim9, Prim10
    real(R8), intent(in), dimension(3,3) :: M1, M2
    real(R8), intent(inout), dimension(nprim) :: Res1, Res2
    ! Local
    real(R8) :: rho1, Rgas, T1, rho2, T2, Gradient(nprim,3), Prim(nprim), M(3,3), waldis, Flux(nprim)

    ! Gradient in the same direction of the face: 1 and 2
    call co_rotot_Rtot ( Prim1(1:nsc), rho1, Rgas )
    T1 = Prim1(np) / ( rho1 * Rgas )
    call co_rotot_Rtot ( Prim2(1:nsc), rho2, Rgas )
    T2 = Prim2(np) / ( rho2 * Rgas )

    Gradient ( 1:nsc, a ) = Prim2 ( 1:nsc ) / rho2 - Prim1 ( 1:nsc ) / rho1 ! concentration gradient
    Gradient ( np, a ) = T2 - T1 ! temperature gradient
    Gradient ( nu:nw, a ) = Prim2 ( nu:nw ) - Prim1 ( nu:nw ) ! velocity gradient
      if (nprim>np) then
      Gradient ( np+1:nprim, a ) = Prim2 ( np+1:nprim ) / rho2 - Prim1 ( np+1:nprim ) / rho1 ! RANS variable gradient
      waldis = 0.5d0 * ( waldis1 + waldis2 ) ! distance to nearest wall
    end if

    ! Gradient in tangential directions: 3-10
    call Tangential_Gradient ( Prim3, Prim4, Prim5, Prim6,  Gradient(:,b) )
    call Tangential_Gradient ( Prim7, Prim8, Prim9, Prim10, Gradient(:,c) )
    
    M = 0.5d0 * ( M1 + M2 )
    Gradient = matmul ( Gradient, M )

    Prim = 0.5d0 * ( Prim1 + Prim2 )
    call Compute_Diffusive_Flux ( Prim, Gradient, area, normal, waldis, Flux, Sc, Sct, Prt, soot_enabled )

    Res1 = Res1 - Flux
    Res2 = Res2 + Flux

  end subroutine Diffusive_Flux

  
  subroutine Tangential_Gradient ( Prim1, Prim2, Prim3, Prim4, Gradient )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4
    real(R8), intent(out), dimension(nprim) :: Gradient
    ! Local
    real(R8) :: rho1, rho2, rho3, rho4, T1, T2, T3, T4, Rgas

    call co_rotot_Rtot ( Prim1(1:nsc), rho1, Rgas )
    T1 = Prim1(np) / ( rho1 * Rgas )
    call co_rotot_Rtot ( Prim2(1:nsc), rho2, Rgas )
    T2 = Prim2(np) / ( rho2 * Rgas )
    call co_rotot_Rtot ( Prim3(1:nsc), rho3, Rgas )
    T3 = Prim3(np) / ( rho3 * Rgas )
    call co_rotot_Rtot ( Prim4(1:nsc), rho4, Rgas )
    T4 = Prim4(np) / ( rho4 * Rgas )

    Gradient ( 1:nsc ) = ( Prim2 ( 1:nsc ) / rho2 - Prim1 ( 1:nsc ) / rho1 + &
                            Prim4 ( 1:nsc ) / rho4 - Prim3 ( 1:nsc ) / rho3 ) * 0.25d0

    Gradient ( np ) = ( T2 - T1 + T4 - T3 ) * 0.25d0

    Gradient ( nu:nw ) = ( Prim2 ( nu:nw ) - Prim1 ( nu:nw ) + Prim4 ( nu:nw ) - Prim3 ( nu:nw ) ) * 0.25d0

    if (nprim>np) then
      Gradient ( np+1:nprim ) = ( Prim2 ( np+1:nprim ) / rho2 - Prim1 ( np+1:nprim ) / rho1 + &
                                  Prim4 ( np+1:nprim ) / rho4 - Prim3 ( np+1:nprim ) / rho3 ) * 0.25d0
    end if

  end subroutine Tangential_Gradient


  subroutine Compute_Diffusive_Flux ( Prim, Gradient, area, normal, waldis, Flux, Sc, Sct, Prt, soot_enabled )
    use MOSE_Global_m
    use MOSE_Lib_Fluid
    use MOSE_Lib_RANS
    use MOSE_Mod_Soot, only: Soot_Diffusive_Flux
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: Prim(nprim), Gradient(nprim,3), area, normal(3), waldis
    real(R8), intent(in)  :: Sc, Sct, Prt
    logical, intent(in)   :: soot_enabled
    real(R8), intent(out) :: Flux(nprim)
    ! Local
    integer :: s, T_i, Tint(2)
    real(R8) :: rho, Rgas, T, Tdiff, cp, mil, kl, mie, kappa
    real(R8) :: Dm(nsc), stress(3), DiffHFlux, DmGradYi

    ! Thermodynamic and transport properties at the interface
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    T_i     = int(T)
    Tdiff   = T - T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1
    cp = f_cp_expr ( Prim(1:nsc), Tint, Tdiff, rho )
    call co_k_mi_lam_Wilke_expr ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )
    
    ! Eddy viscosity
    mie = 0d0
    if (model==2) then
      call Eddy_Viscosity ( mut=mie, rans_variables=Prim(nt:nprim), &
                            mul=mil, rho=rho, vel_gradient=Gradient(nu:nw,:), &
                            walldist=waldis )
    end if

    Dm (1:nsc) = ( mil/Sc + mie/Sct ) / rho ! binary coefficient computed from Schmidt = mi/rho*Dm

    kappa = kl + mie*cp/Prt ! Laminar + turbulent conductivity

    Stress = stress_vector ( Gradient(nu:nw,:), normal, mil, mie, prim(nt:) ) ! Stress tensor in cartesian components

    ! Fluxes computation
    DiffHFlux = 0.0d0
    Flux = 0.0d0

    ! Diffusive mass fluxes
    do s = 1, nsc
      DmGradYi  = Dm(s) * dot_product ( Gradient(s,:), normal )     ! Dm * dyi/dn
      Flux(s) = rho * area * DmGradYi                               ! rho * Dm *dyi/dn * Area
      DiffHFlux = DiffHFlux + DmGradYi * f_tabT_expr ( s, h_tab, Tint, Tdiff )
    enddo

    ! Correcting the mass flux for possible errors (sum must be zero)
    Flux(1:nsc) = Flux(1:nsc) - sum(Flux(1:nsc))*Prim(1:nsc)/rho

    ! Momentum flux
    Flux(nu:nw) = area * Stress

    ! Diffusive enthalpy flux
    Flux(np) = area * ( dot_product ( Stress, Prim(nu:nw) ) + &
                kappa * dot_product ( Gradient(np,:), normal ) + DiffHFlux * rho )

    ! Turbulence variables diffusive flux
    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho, &
                                 area=area, normal=normal, dist=waldis )
    end if

    ! Soot diffusive fluxes
    if (soot_enabled) then
      call Soot_Diffusive_Flux ( flux=Flux, &
                                 prim=Prim, &
                                 T_gradient=Gradient(np,:), &
                                 mu=mil+mie, &
                                 rho=rho, &
                                 T=T, &
                                 area=area, & 
                                 normal=normal )
    end if

  end subroutine Compute_Diffusive_Flux

end module MOSE_Lib_Diffusive