module MOSE_Lib_BC_Fluxes_GSI_Reactions
  !< Heterogeneous surface reaction wall boundary condition (BC 505).
  !<
  !< The wall material reacts with the incoming gas, consuming some species and
  !< releasing others. Unlike pyrolysis the products follow from the reaction
  !< mechanism itself, so no product composition has to be supplied. The wall
  !< temperature is again the root of the surface energy balance; the machinery is
  !< shared with the pyrolysis boundary condition and lives in
  !< `MOSE_Lib_BC_Fluxes_GSI_Core`.

  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Lib_BC_Fluxes_GSI_Core, only: BC_Wall_Ablation, GSI_KIND_REACTIONS

  implicit none
  private

  public :: BC_Wall_Reactions

contains

  subroutine BC_Wall_Reactions ( Im, Jm, Km, Fm, Blk, GSI_id, qrad, Ovar, ablating, T_surface )
    !< Apply the surface reaction wall boundary condition to face `Fm` of cell `(Im,Jm,Km)`.
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    integer, intent(in) :: GSI_id       !< reaction model, see `GSI_REAC_*`
    real(R8), intent(in) :: qrad        !< radiative heat flux at the wall
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    logical,  intent(out), optional :: ablating    !< .false. if the surface is inert
    real(R8), intent(out), optional :: T_surface   !< wall temperature for the caller's fallback

    call BC_Wall_Ablation ( Im, Jm, Km, Fm, Blk, GSI_KIND_REACTIONS, GSI_id, qrad, Ovar=Ovar, &
                            ablating=ablating, T_surface=T_surface )

  end subroutine BC_Wall_Reactions

end module MOSE_Lib_BC_Fluxes_GSI_Reactions
