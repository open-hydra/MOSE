module MOSE_Lib_BC_Fluxes_GSI_Melting
  !< Melting wall boundary condition (BC 503).
  !<
  !< The wall temperature is prescribed - it is the melting temperature of the solid -
  !< so, unlike pyrolysis and surface reactions, no iteration is needed: the heat
  !< reaching the surface is known straight away and the melting mass flux follows
  !< from the surface energy balance. The wall state and the flux assembly are shared
  !< with the other ablative boundary conditions, in
  !< `MOSE_Lib_BC_Fluxes_GSI_Core`.

  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Lib_GSI, only: GSI_Melting
  use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot
  use MOSE_Lib_BC_Fluxes_GSI_Core, only: wall_face_t, Get_Wall_Face, Wall_Gas_State, &
                                              Ablation_Wall_Fluxes, RANS_WALL_NO_BLOWING, &
                                              Surface_Is_Ablating

  implicit none
  private

  public :: BC_Wall_Melting

contains

  subroutine BC_Wall_Melting ( Im, Jm, Km, Fm, Blk, cp_wall, T_wall, T_i_wall, Dh_wall, q_rad, y_wall, &
                               Ovar, ablating, T_surface )
    !< Apply the melting wall boundary condition to face `Fm` of cell `(Im,Jm,Km)`.
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: cp_wall    !< specific heat of the solid
    real(R8), intent(in) :: T_wall     !< prescribed wall (melting) temperature
    real(R8), intent(in) :: T_i_wall   !< initial temperature of the solid
    real(R8), intent(in) :: Dh_wall    !< latent heat of melting
    real(R8), intent(in) :: q_rad      !< radiative heat flux at the wall
    real(R8), intent(in), dimension(1:nsc) :: y_wall !< mass fractions of the melted material
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    logical,  intent(out), optional :: ablating    !< .false. if the surface is inert
    real(R8), intent(out), optional :: T_surface   !< wall temperature for the caller's fallback
    ! Local
    logical :: ablates
    type(wall_face_t) :: face
    real(R8) :: Prim(nprim), rho, Rgas, T_cell
    real(R8) :: rho_wall, Prim_wall(nprim), mu_wall, k_wall
    real(R8) :: grad_T(3), q_conv, q_rad_face, omega(nsc), mdot

    call Get_Wall_Face ( Im, Jm, Km, Fm, Blk, face )

    ! Boundary cell state
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T_cell = Prim(np) / ( rho * Rgas )

    ! The wall temperature is known, so the wall state and the conductive flux follow
    ! in a single pass.
    call Wall_Gas_State ( face, Prim, rho, Rgas, T_cell, T_wall, &
                          rho_wall, Prim_wall(1:nsc), mu_wall, k_wall, grad_T, q_conv )

    ! Orient the radiative flux like the conductive one. q_conv already carries the
    ! factor (-modfm2), through modfm3 = -2*modfm2, so the incoming radiation has to
    ! carry it too for their sum to be the heat absorbed by this particular face.
    ! Because mdot follows from that sum, it -- and with it omega -- comes out
    ! correctly signed for the face, the way GSI_Source signs the pyrolysis models.
    ! On a low face (-modfm2) is +1, so only the high faces change.
    q_rad_face = q_rad * (-face % modfm2)

    ! The heat reaching the surface is what melts the solid and sets the mass flux.
    call GSI_Melting ( cp_wall, T_wall, T_i_wall, Dh_wall, q_conv + q_rad_face, y_wall, omega, mdot )

    ! Nothing melts unless the surface is a net receiver of heat: with q_conv + q_rad
    ! negative the energy balance returns a non-positive mdot, which would have the
    ! wall swallowing gas. Leave the residual alone and let the caller close the
    ! boundary with an impermeable wall at the melting temperature instead.
    ablates = Surface_Is_Ablating ( face, mdot, T_wall )
    if ( present(ablating)  ) ablating  = ablates
    if ( present(T_surface) ) T_surface = T_wall
    if ( .not. ablates ) return

    ! The surface absorbs q_conv + q_rad and spends it melting the solid, but only
    ! q_conv comes out of the gas: the radiation arrives from outside the domain.
    call Ablation_Wall_Fluxes ( Im, Jm, Km, face, Blk, Prim, rho, T_wall, rho_wall, Prim_wall, &
                                mu_wall, grad_T, q_gas=q_conv, q_wall=q_conv + q_rad_face, &
                                omega=omega, mdot=mdot, &
                                rans_wall_model=RANS_WALL_NO_BLOWING, Ovar=Ovar )

  end subroutine BC_Wall_Melting

end module MOSE_Lib_BC_Fluxes_GSI_Melting
