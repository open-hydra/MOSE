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
  ! The gas enters when p0 is above the pressure that the gas of the cell reaches at rest
  ! (Stagn_Inflow), whatever the sign of the cell velocity; otherwise the face is an outlet at
  ! BC_pback (Outflow, or Outlet_Moving_In when the gas of the cell moves in).
  subroutine BC_FarField ( Bm, Im, Jm, Km, Fm, Blk, BC_T0, BC_p0, BC_pback, &
                                 BC_rel_fac, BC_alpha, BC_beta, BC_ci, BC_RANS, BC_mdot )
    use MOSE_Config_Types_m,        only: obj_rot
    use MOSE_Lib_RotatingFrame,     only: RF_Convert_BC_Inflow
    use FLINT_Lib_Thermodynamic,    only: f_Rtot
    use MOSE_Lib_BC_Fluxes_Inflow,  only: Stagn_Inflow, setup_inflow_geometry, Compute_Inflow_Direction
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
    logical  :: enters

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Int_k, Normal, Area, t_Vec, BC_Sign,   &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,    &
                                  Bound_Sound, Bound_Gamma, Riem )

    T0_loc    = BC_T0
    p0_loc    = BC_p0
    Mach_loc  = 0.0d0
    ! The token normal gives the face normal, as for the other inlets; resolved before the
    ! conversion of the rotating frame, which turns the direction it is given
    call Compute_Inflow_Direction ( Normal, BC_alpha, BC_beta, alpha_loc, beta_loc )

    Flux = 0d0

    ! Inflow when the ambient drives gas in (Stagn_Inflow decides, with the total conditions seen
    ! in the rotating frame), otherwise the back-pressure outlet treatment
    if (obj_rot%enabled) then
      Face_Rgas_loc = f_Rtot(BC_ci)
      call RF_Convert_BC_Inflow( Blk=Blk, Fm=Fm, Im=Im, Jm=Jm, Km=Km,   &
                                  BC_Mach_abs=0.0d0, gamma=Bound_Gamma,   &
                                  Rgas=Face_Rgas_loc, T0_or_Tstat=T0_loc, &
                                  p0_or_pstat=p0_loc, alpha=alpha_loc,     &
                                  beta=beta_loc, Mach_rel=Mach_loc )
    end if
    call Stagn_Inflow ( Bound_Prim, Bound_Sound, Bound_Gamma, Riem,           &
                               alpha_loc, beta_loc, BC_ci, BC_RANS, T0_loc, p0_loc, &
                               BC_rel_fac, BC_Sign, t_Vec, Normal, Area, Flux, BC_mdot, enters )
    if ( .not. enters ) then
      if ( BC_Sign * Un >= 0d0 ) then
        ! The gas of the cell is at rest or moves out: the outlet treatment of 406 at BC_pback
        call Outflow( Bound_Prim, BC_Sign, Un, Normal, Area, BC_pback, Flux, BC_mdot )
      else
        ! The gas of the cell moves in, and the ambient does not feed it. Outflow would raise the
        ! exit pressure to the critical pressure of the cell, computed with |Un|, which is above
        ! p_stop from a normal Mach number of about 0.37 (gamma = 1.4), and from normal Mach 1 it
        ! takes the state of the cell: the face would let in the gas that the ambient cannot feed
        call Outlet_Moving_In( Bound_Prim, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, &
                               BC_Sign, Un, Normal, Area, BC_pback, Flux, BC_mdot )
      end if
    end if

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux

  end subroutine BC_FarField


  !─────────────────────────────────────────────────────────────────────────────
  ! Outlet of a 407 face whose cell gas moves in while the ambient does not feed it (p0 <= p_stop).
  ! The face state lies on the invariant Riem = a + g w that reaches the face from inside (w along
  ! the outward normal, g = (gamma-1)/2) and on the isentrope of the cell: at the face pressure p_f,
  ! a_f = a (p_f/p)**((gamma-1)/(2 gamma)) and w_f = (Riem - a_f)/g. So w_f >= 0 exactly when
  ! a_f <= Riem, that is p_f <= p_stop = p (Riem/a)**(2 gamma/(gamma-1)), the pressure at which the
  ! outgoing wave brings the gas of the cell to rest. The face takes p_f = min(BC_pback, p_stop): the
  ! gas leaves at BC_pback, and with BC_pback at or above p_stop the face is at rest at p_stop.
  ! With Riem <= 0 (gas moving in at 2 a/(gamma-1) or faster) every p_f > 0 gives w_f < 0: no face
  ! state lets the gas out, and the limit p_f -> 0, the vacuum that the gas leaves behind, carries no
  ! flux (rho_f -> 0 with w_f bounded); so does a BC_pback <= 0.
  subroutine Outlet_Moving_In ( Bound_Prim, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, &
                                BC_Sign, Un, Normal, Area, BC_pback, Flux, BC_mdot )
    use FLINT_Lib_Thermodynamic, only: f_tabT, h_tab
    implicit none
    real(R8), intent(in)  :: Bound_Prim(nprim), Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8), intent(in)  :: BC_Sign, Un, Normal(3), Area, BC_pback
    real(R8), intent(out) :: Flux(nprim), BC_mdot
    ! Local
    integer  :: s
    real(R8) :: p_stop, Face_Sound, Face_Un, Face_rho, Face_T, Face_Enthalpy, Face_Prim(nprim)

    Flux    = 0d0
    BC_mdot = 0d0
    if ( Riem <= 0d0 .or. BC_pback <= 0d0 ) return

    p_stop = Bound_Prim(np) * ( Riem / Bound_Sound )**( 2d0 * Bound_Gamma / ( Bound_Gamma - 1d0 ) )
    if ( BC_pback < p_stop ) then
      Face_Prim(np) = BC_pback
      Face_Sound    = Bound_Sound * ( BC_pback / Bound_Prim(np) )**( ( Bound_Gamma - 1d0 ) / ( 2d0 * Bound_Gamma ) )
    else
      Face_Prim(np) = p_stop
      Face_Sound    = Riem                      ! the gas at rest at the face: w_f = 0
    end if
    Face_Un  = BC_Sign * ( Riem - Face_Sound ) / ( 0.5d0 * ( Bound_Gamma - 1d0 ) )
    Face_rho = Bound_Gamma * Face_Prim(np) / Face_Sound**2
    Face_Prim(1:nsc) = Bound_Prim(1:nsc) / Bound_rho * Face_rho             ! mass fractions of the cell
    Face_Prim(nu:nw) = Bound_Prim(nu:nw) - Un * Normal + Face_Un * Normal   ! tangential velocity of the cell
    if ( nprim > np ) Face_Prim(np+1:nprim) = Bound_Prim(np+1:nprim) / Bound_rho * Face_rho

    Face_T = Face_Prim(np) / ( Bound_Rgas * Face_rho )
    Face_Enthalpy = 0d0
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT( Face_T, s, h_tab )
      Flux(s) = Face_Prim(s) * Area * Face_Un
    end do
    BC_mdot = Face_rho * Area * Face_Un
    Flux(nu:nw) = Face_Prim(nu:nw) * BC_mdot + Face_Prim(np) * Area * Normal
    Flux(np) = BC_mdot * ( 0.5d0 * sum( Face_Prim(nu:nw)**2 ) + Face_Enthalpy )
    if ( nprim > np ) Flux(np+1:nprim) = Face_Prim(np+1:nprim) / Face_rho * BC_mdot

  end subroutine Outlet_Moving_In

end module MOSE_Lib_BC_Fluxes_FarField