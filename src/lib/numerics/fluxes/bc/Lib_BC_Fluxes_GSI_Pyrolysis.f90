module MOSE_Lib_BC_Fluxes_GSI_Pyrolysis
  !< Pyrolysing wall boundary condition (BC 504).
  !<
  !< A solid fuel decomposes under an Arrhenius law in the wall temperature and blows
  !< its pyrolysis products into the domain. The wall temperature is unknown and is
  !< obtained from the surface energy balance; all of that machinery is shared with
  !< the surface reaction boundary condition and lives in
  !< `MOSE_Lib_BC_Fluxes_GSI_Core`.

  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Lib_BC_Fluxes_GSI_Core, only: BC_Wall_Ablation, GSI_KIND_PYROLYSIS

  implicit none
  private

  public :: BC_Wall_Pyrolysis

contains

  subroutine BC_Wall_Pyrolysis ( Im, Jm, Km, Fm, Blk, GSI_id, qrad, y_wall, Ovar, ablating, T_surface )
    !< Apply the pyrolysing wall boundary condition to face `Fm` of cell `(Im,Jm,Km)`.
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    integer, intent(in) :: GSI_id                   !< pyrolysis model, see `GSI_PYRO_*`
    real(R8), intent(in) :: qrad                    !< radiative heat flux at the wall
    real(R8), intent(in), dimension(1:nsc) :: y_wall !< pyrolysis product mass fractions
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    logical,  intent(out), optional :: ablating    !< .false. if the surface is inert
    real(R8), intent(out), optional :: T_surface   !< wall temperature for the caller's fallback

    call BC_Wall_Ablation ( Im, Jm, Km, Fm, Blk, GSI_KIND_PYROLYSIS, GSI_id, qrad, &
                            Ovar=Ovar, y_wall=y_wall, ablating=ablating, T_surface=T_surface )

  end subroutine BC_Wall_Pyrolysis

end module MOSE_Lib_BC_Fluxes_GSI_Pyrolysis
