module MOSE_Lib_BC_Fluxes_GSI_Char_Pyrolysis
  !< Pyrolysing ablator wall boundary condition (BC 506).
  !<
  !< A charring ablator: its char reacts with the incoming gas, as the wall of BC 505, and the gas of the pyrolysing
  !< virgin material is blown through it, mdot_pg = phi mdot_ch, with the equilibrium composition of a table of
  !< INPUT/<prefix>pyrolysis.dat at the wall pressure and temperature (GSI_Char_Pyrolysis, MOSE_Lib_GSI). The surface
  !< energy balance is that of steady conduction in the solid, the heat the surface absorbs less mdot h_v0, the enthalpy
  !< of the virgin material. The machinery is that of BC 505 (MOSE_Lib_BC_Fluxes_GSI_Core); with phi = 0 and h_v0 = 0
  !< this boundary condition is BC 505 bit for bit.

  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Lib_BC_Fluxes_GSI_Core, only: BC_Wall_Ablation, GSI_KIND_CHAR_PYROLYSIS, gsi_pyrolysis_t

  implicit none
  private

  public :: BC_Wall_Char_Pyrolysis

contains

  subroutine BC_Wall_Char_Pyrolysis ( Im, Jm, Km, Fm, Blk, GSI_id, pyro_table, phi, h_v0, qrad, Ovar, ablating, T_surface )
    !< Apply the pyrolysing ablator wall boundary condition to face `Fm` of cell `(Im,Jm,Km)`.
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    integer, intent(in) :: GSI_id       !< reaction model of the char, see `GSI_REAC_*`
    integer, intent(in) :: pyro_table   !< table of INPUT/<prefix>pyrolysis.dat, read when phi > 0
    real(R8), intent(in) :: phi         !< mdot_pg / mdot_ch
    real(R8), intent(in) :: h_v0        !< enthalpy of the virgin material, on the datum of the gas tables [J/kg]
    real(R8), intent(in) :: qrad        !< radiative heat flux at the wall
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    logical,  intent(out), optional :: ablating    !< .false. if the surface is inert
    real(R8), intent(out), optional :: T_surface   !< wall temperature for the caller's fallback

    call BC_Wall_Ablation ( Im, Jm, Km, Fm, Blk, GSI_KIND_CHAR_PYROLYSIS, GSI_id, qrad, Ovar=Ovar, &
                            ablating=ablating, T_surface=T_surface, pyro=gsi_pyrolysis_t(pyro_table, phi, h_v0) )

  end subroutine BC_Wall_Char_Pyrolysis

end module MOSE_Lib_BC_Fluxes_GSI_Char_Pyrolysis
