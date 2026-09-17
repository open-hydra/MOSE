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

  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot
  use MOSE_Lib_BC_Fluxes_GSI_Core, only: wall_face_t, Get_Wall_Face, Wall_Gas_State, &
                                         Ablation_Wall_Fluxes, Surface_Is_Ablating, &
                                         RANS_WALL_BLOWING

  implicit none
  private

  public :: BC_SRM

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
    real(R8) :: grad_T(3), q_conv, omega(nsc), mdot

    call Get_Wall_Face ( Im, Jm, Km, Fm, Blk, face )

    ! Boundary cell state
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T_cell = Prim(np) / ( rho * Rgas )

    ! The wall temperature is the flame temperature, known up front, so the wall state
    ! and the conductive flux follow in a single pass. The wall sees the combustion
    ! products, not the chamber gas, hence the prescribed composition.
    call Wall_Gas_State ( face, Prim, rho, Rgas, T_cell, Taf, &
                          rho_wall, Prim_wall(1:nsc), mu_wall, k_wall, grad_T, q_conv, &
                          y_wall = y_af )

    ! Saint-Robert's law, signed by the face orientation the way GSI_Source signs the
    ! ablation models: (-modfm2) is +1 on a low face and -1 on a high one, so the mass
    ! is injected into the domain whichever side of the block the grain sits on.
    mdot  = aCoeff * ( Prim(np) / pRef )**n * (-face % modfm2)
    omega = mdot * y_af

    ! An inhibited surface (a = 0) injects nothing, and a failed state must not be
    ! turned into a flux. Leave the residual alone and let the caller close the
    ! boundary with an impermeable wall at the flame temperature instead.
    burns = Surface_Is_Ablating ( face, mdot, Taf )
    if ( present(burning)   ) burning   = burns
    if ( present(T_surface) ) T_surface = Taf
    if ( .not. burns ) return

    ! The gas loses q_conv to the surface and gets back the enthalpy of the injected
    ! products, evaluated at the flame temperature from the thermodynamic tables. No
    ! radiative term: the grain is not heated from outside the domain.
    call Ablation_Wall_Fluxes ( Im, Jm, Km, face, Blk, Prim, rho, Taf, rho_wall, Prim_wall, &
                                mu_wall, grad_T, q_gas=q_conv, q_wall=q_conv, &
                                omega=omega, mdot=mdot, &
                                rans_wall_model=RANS_WALL_BLOWING, Ovar=Ovar )

  end subroutine BC_SRM

end module MOSE_Lib_BC_Fluxes_GSI_SRM
