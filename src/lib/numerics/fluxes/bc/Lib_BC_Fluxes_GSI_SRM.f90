module MOSE_Lib_BC_Fluxes_GSI_SRM
  !< Solid rocket motor grain combustion boundary condition (BC 502).
  !<
  !< A solid propellant carries its own oxidizer, so the burning surface is a
  !< deflagration front rather than a gasification front: the rate is set by the
  !< flame standing off the grain, not by the surface energy balance. Both the rate
  !< and the flame temperature are therefore measured properties of the propellant,
  !< supplied per boundary, and no wall temperature iteration is needed.
  !<
  !< The regression rate follows Saint-Robert's (Vieille's) law, evaluated in mass
  !< flux form and normalised by a reference pressure,
  !<
  !<   mdot = a * ( p_cell / p_ref )**n
  !<
  !< so `a` is the propellant mass flux at `p_ref`, in kg/(m^2 s). The grain density
  !< is *not* applied here - a linear burn rate coefficient has to be multiplied by
  !< it before it is handed over.
  !<
  !< Everything downstream of the rate law - wall gas state, blowing velocity,
  !< viscous stress, flux assembly, RANS wall treatment - is the blowing wall shared
  !< with the ablative boundary conditions, in `MOSE_Lib_BC_Fluxes_GSI_Core`.
  !<
  !< Two features distinguish this boundary condition from the ablative ones:
  !<
  !<  * the surface imposes the *combustion product* composition on the wall state,
  !<    rather than inheriting the composition of the adjacent cell. The gas touching
  !<    a burning grain is the flame product mixture, and the wall transport
  !<    properties that set the conductive flux follow from it;
  !<  * the heat reaching the surface plays no part in setting the rate, so `q_conv`
  !<    is passed through to the gas energy flux but never fed back.
  !<
  !< The flame temperature T_af of the record is the total temperature of the products,
  !< that of the gas at rest at the end of the flame (the chamber of an equilibrium
  !< calculation): the products enter with the total enthalpy h(T_af) per unit mass,
  !< that of the propellant, at the static temperature of `SRM_Static_Temperature`. In
  !< mode 2 (an inflow) the grain is the inlet 404 at that static temperature (BC_Fluxes).

  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, f_tabT, h_tab, cp_tab, Ri_tab
  use MOSE_Lib_BC_Fluxes_GSI_Core, only: wall_face_t, Get_Wall_Face, Wall_Gas_State, &
                                         Ablation_Wall_Fluxes, Surface_Is_Ablating, &
                                         RANS_WALL_BLOWING

  implicit none
  private

  public :: BC_SRM, SRM_Wall_Turbulence, SRM_Static_Temperature

contains

  subroutine BC_SRM ( Im, Jm, Km, Fm, Blk, aCoeff, n, pRef, Taf, y_af, Ovar, burning, T_surface )
    !< Apply the grain combustion boundary condition to face `Fm` of cell `(Im,Jm,Km)`.
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: aCoeff   !< burn rate coefficient, as a mass flux at `pRef` (kg/(m^2 s))
    real(R8), intent(in) :: n        !< pressure exponent
    real(R8), intent(in) :: pRef     !< reference pressure of the burn rate law
    real(R8), intent(in) :: Taf      !< adiabatic flame temperature, imposed at the wall
    real(R8), intent(in), dimension(1:nsc) :: y_af !< mass fractions of the combustion products
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    logical,  intent(out), optional :: burning     !< .false. if the surface does not burn
    real(R8), intent(out), optional :: T_surface   !< wall temperature for the caller's fallback
    ! Local
    logical :: burns
    type(wall_face_t) :: face
    real(R8) :: Prim(nprim), rho, Rgas, T_cell
    real(R8) :: rho_wall, Prim_wall(nprim), mu_wall, k_wall
    real(R8) :: grad_T(3), q_conv, omega(nsc), mdot, T_w

    call Get_Wall_Face ( Im, Jm, Km, Fm, Blk, face )

    ! Boundary cell state
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T_cell = Prim(np) / ( rho * Rgas )

    ! Saint-Robert's law, signed by the face orientation the way GSI_Source signs the
    ! ablation models: (-modfm2) is +1 on a low face and -1 on a high one, so the mass
    ! is injected into the domain whichever side of the block the grain sits on.
    mdot  = aCoeff * ( Prim(np) / pRef )**n * (-face % modfm2)
    omega = mdot * y_af

    ! T_af is the total temperature of the products: they leave the wall at the static
    ! temperature T_w of h(T_w) + u_b^2/2 = h(T_af), u_b = mdot/rho_wall, so that the gas
    ! gets back the total enthalpy h(T_af) per unit mass, that of the propellant. T_w, known
    ! up front, gives the wall state and the conductive flux in a single pass. The wall sees
    ! the combustion products, not the chamber gas, hence the prescribed composition.
    T_w = SRM_Static_Temperature ( Taf, abs(mdot), Prim(np), y_af )
    call Wall_Gas_State ( face, Prim, rho, Rgas, T_cell, T_w, &
                          rho_wall, Prim_wall(1:nsc), mu_wall, k_wall, grad_T, q_conv, &
                          y_wall = y_af )

    ! An inhibited surface (a = 0) injects nothing, and a failed state must not be
    ! turned into a flux. Leave the residual alone and let the caller close the
    ! boundary with an impermeable wall at the flame temperature instead (no flow:
    ! static and total temperature are T_af).
    burns = Surface_Is_Ablating ( face, mdot, Taf )
    if ( present(burning)   ) burning   = burns
    if ( present(T_surface) ) T_surface = Taf
    if ( .not. burns ) return

    ! The gas loses q_conv to the surface and gets back the enthalpy of the injected
    ! products at T_w from the thermodynamic tables and their kinetic energy, h(T_af) per
    ! unit mass. No radiative term: the grain is not heated from outside the domain.
    call Ablation_Wall_Fluxes ( Im, Jm, Km, face, Blk, Prim, rho, T_w, rho_wall, Prim_wall, &
                                mu_wall, grad_T, q_gas=q_conv, q_wall=q_conv, &
                                omega=omega, mdot=mdot, &
                                rans_wall_model=RANS_WALL_BLOWING, Ovar=Ovar )

  end subroutine BC_SRM


  subroutine SRM_Wall_Turbulence ( Im, Jm, Km, Fm, Blk, aCoeff, n, pRef, Taf, y_af, rho_wall, rans_wall )
    !< Wall values of the turbulence variables of face `Fm` of cell `(Im,Jm,Km)`, those that the boundary flux of the
    !< grain takes in the same stage, for its ghost cells (Ghost_SRM_Turbulence in `MOSE_Lib_Ghost`): while the surface
    !< burns, the blowing wall of `Ablation_Wall_Fluxes` (the wall state of `BC_SRM`, the stress of the blowing velocity
    !< gradient, `RANS_Set_Blowing_Wall`); when it does not, the impermeable wall at the flame temperature that the caller
    !< of `BC_SRM` closes the face with (`BC_Wall_Temperature`: `RANS_Set_Wall_Values`, smooth).
    use MOSE_Lib_BC_Fluxes, only: Wall_State_Temperature
    use MOSE_Lib_Fluid, only: Stress_Vector
    use MOSE_Lib_RANS, only: RANS_Set_Blowing_Wall, RANS_Set_Wall_Values
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(in) :: Blk
    real(R8), intent(in) :: aCoeff, n, pRef, Taf      !< as in `BC_SRM`
    real(R8), intent(in), dimension(1:nsc) :: y_af   !< mass fractions of the combustion products
    real(R8), intent(out) :: rho_wall                !< density of the wall state
    real(R8), intent(out) :: rans_wall(nrans)        !< turbulence variables at the wall, times rho_wall
    ! Local
    type(wall_face_t) :: face
    real(R8) :: Prim(nprim), rho, Rgas, T_cell, roi_wall(nsc), mu_wall, k_wall, grad_T(3), q_conv, mdot
    real(R8) :: Blowing(3), VelGrad(3,3), Stress(3), Rans_zero(max(nRANS,1))

    call Get_Wall_Face ( Im, Jm, Km, Fm, Blk, face )
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T_cell = Prim(np) / ( rho * Rgas )
    mdot = aCoeff * ( Prim(np) / pRef )**n * (-face % modfm2)
    call Wall_Gas_State ( face, Prim, rho, Rgas, T_cell, SRM_Static_Temperature ( Taf, abs(mdot), Prim(np), y_af ), &
                          rho_wall, roi_wall, mu_wall, k_wall, grad_T, q_conv, y_wall = y_af )

    if ( Surface_Is_Ablating ( face, mdot, Taf ) ) then
      ! the operations of Ablation_Wall_Fluxes
      Blowing = mdot / rho_wall * face % normal
      VelGrad = 0d0
      VelGrad(:,face % dir) = ( Prim(nu:nw) - Blowing ) * face % modfm3
      VelGrad = matmul ( VelGrad, face % M )
      ! the turbulence variables entering the stress (Reynolds stresses, k) are zero at the wall
      Rans_zero = 0d0
      Stress = Stress_Vector ( VelGrad, face % normal, mu_wall, 0d0, Rans_zero )
      call RANS_Set_Blowing_Wall ( rho=rho_wall, mil=mu_wall, rans_variables=rans_wall, &
                                   tau=Stress, mdot=mdot, dist=Blk % yn(Im,Jm,Km) )
    else
      call Wall_State_Temperature ( Prim, rho, Rgas, Taf, rho_wall, roi_wall, mu_wall, k_wall )
      call RANS_Set_Wall_Values ( mu_wall, Prim(nt:nprim) * rho_wall / rho, rans_wall, Blk % yn(Im,Jm,Km), 0d0 )
    end if

  end subroutine SRM_Wall_Turbulence

  pure real(R8) function SRM_Static_Temperature ( T0, g, p, y ) result ( T )
    !< Static temperature of the products of mass fractions `y` that cross a face with the mass flux `g` >= 0
    !< (kg/(m^2 s)) at the pressure `p` and have the total temperature `T0` (the T_af of a grain, 502): the root of
    !< h(T) + u(T)^2/2 = h(T0), u = g R T / p, with h and R of the tables of the gas (f_tabT, Ri_tab), as every state of
    !< MOSE, so that the gas carries the total enthalpy h(T0). The left side grows with T and is not below h(T0) at T0:
    !< Newton from T0 comes down to the root (T0 - T = u^2/(2 cp): 2.5e-3 K at 20 kg/(m^2 s) and 7 MPa for products of
    !< 371 J/(kg K)), in two or three steps; cp of the tables as the slope. g = 0 gives T0.
    implicit none
    real(R8), intent(in) :: T0, g, p, y(nsc)
    ! Local
    integer  :: it, s
    real(R8) :: R, a2, h0, F, dF, dT

    T = T0
    if ( g == 0d0 ) return
    R  = sum ( y * Ri_tab(1:nsc) )
    a2 = ( g * R / p )**2                       ! u^2 = a2 T^2
    h0 = 0d0
    do s = 1, nsc
      h0 = h0 + y(s) * f_tabT ( T0, s, h_tab )
    end do
    do it = 1, 50
      F  = 0.5d0 * a2 * T * T - h0
      dF = a2 * T
      do s = 1, nsc
        F  = F  + y(s) * f_tabT ( T, s, h_tab )
        dF = dF + y(s) * f_tabT ( T, s, cp_tab )
      end do
      dT = F / dF
      T  = T - dT
      if ( abs(dT) <= 1d-12 * T0 ) exit
    end do

  end function SRM_Static_Temperature

end module MOSE_Lib_BC_Fluxes_GSI_SRM
