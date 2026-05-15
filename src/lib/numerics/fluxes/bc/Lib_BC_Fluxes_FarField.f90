module MOSE_Lib_BC_Fluxes_FarField
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  private

  public :: BC_FarField

contains


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 407: Ambient BC — T0 + p0 (inflow) or back pressure BC_pback (outflow).
  ! Flow direction is determined from the interior velocity field.
  subroutine BC_FarField ( Bm, Im, Jm, Km, Fm, Blk, BC_T0, BC_p0, BC_pback, &
                                 BC_rel_fac, BC_alpha, BC_beta, BC_ci, BC_RANS, BC_mdot )
    use MOSE_Config_Types_m,        only: obj_rot
    use MOSE_Lib_RotatingFrame,     only: RF_Convert_BC_Inflow
    use FLINT_Lib_Thermodynamic,    only: f_Rtot
    use MOSE_Lib_BC_Fluxes_Inflow,  only: Stagn_Inflow, setup_inflow_geometry
    use MOSE_Lib_BC_Fluxes_Outflow, only: Outflow
    implicit none
    type(MOSE_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Km, Fm
    real(R8), intent(in)    :: BC_T0, BC_p0, BC_pback, BC_rel_fac, BC_alpha, BC_beta
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    real(R8), intent(out)   :: BC_mdot
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k
    real(R8) :: Normal(3), Area, t_Vec(3), BC_Sign
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8) :: T0_loc, p0_loc, Mach_loc, alpha_loc, beta_loc, Face_Rgas_loc
    real(R8) :: Flux(nprim)

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Int_k, Normal, Area, t_Vec, BC_Sign,   &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,    &
                                  Bound_Sound, Bound_Gamma, Riem )

    T0_loc    = BC_T0
    p0_loc    = BC_p0
    Mach_loc  = 0.0d0
    alpha_loc = BC_alpha
    beta_loc  = BC_beta

    Flux = 0d0

    if (BC_Sign * Un >= 0d0) then
      ! Outflow: apply back-pressure outlet treatment
      call Outflow( Bound_Prim, BC_Sign, Un, Normal, Area, BC_pback, Flux, BC_mdot )
    else
      ! Inflow
      if (obj_rot%enabled) then
        Face_Rgas_loc = f_Rtot(BC_ci)
        call RF_Convert_BC_Inflow( Blk=Blk, Fm=Fm, Im=Im, Jm=Jm, Km=Km,   &
                                    BC_Mach_abs=0.0d0, gamma=Bound_Gamma,   &
                                    Rgas=Face_Rgas_loc, T0_or_Tstat=T0_loc, &
                                    p0_or_pstat=p0_loc, alpha=alpha_loc,     &
                                    beta=beta_loc, Mach_rel=Mach_loc )
      end if
      call Stagn_Inflow ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem,       &
                                 alpha_loc, beta_loc, BC_ci, BC_RANS, T0_loc, p0_loc, &
                                 BC_rel_fac, BC_Sign, t_Vec, Normal, Area, Flux, BC_mdot )
    end if

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux

  end subroutine BC_FarField

end module MOSE_Lib_BC_Fluxes_FarField