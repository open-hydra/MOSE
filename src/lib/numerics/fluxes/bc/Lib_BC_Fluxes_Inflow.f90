module MOSE_Lib_BC_Fluxes_Inflow
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  private

  public :: BC_Inlet_StagnCond, BC_Inlet_MassFlux_T0, BC_Inlet_MassFlux_T, BC_Inlet_Supersonic_Static
  public :: Stagn_Inflow, setup_inflow_geometry

contains

  !─────────────────────────────────────────────────────────────────────────────
  ! BC 401/402: Subsonic inlet — stagnation total temperature T0 and total
  ! pressure p0.  Flow direction from BC_alpha/BC_beta (or face-normal if specified).
  subroutine BC_Inlet_StagnCond ( Bm, Im, Jm, Km, Fm, Blk, BC_T0, BC_p0, &
                                   BC_rel_fac, BC_alpha, BC_beta, BC_ci, BC_RANS, BC_mdot, error )
    use MOSE_Config_Types_m, only: obj_rot
    use MOSE_Lib_RotatingFrame, only: RF_Convert_BC_Inflow
    use FLINT_Lib_Thermodynamic
    implicit none
    type(MOSE_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Km, Fm
    real(R8), intent(in)    :: BC_T0, BC_p0, BC_rel_fac, BC_alpha, BC_beta
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    real(R8), intent(out)   :: BC_mdot
    integer,  intent(inout) :: error
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k
    real(R8) :: Normal(3), Area, t_Vec(3), BC_Sign, alpha, beta
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8) :: T0_loc, p0_loc, alpha_loc, beta_loc, Face_Rgas_loc, Mach_rel_dummy
    real(R8) :: Flux(nprim)

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Int_k, Normal, Area, t_Vec, BC_Sign,   &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,    &
                                  Bound_Sound, Bound_Gamma, Riem )

    ! Reject if flow is directed outward (pure inflow BC)
    if (BC_Sign * Un >= 0d0) then
      error = 1
      return
    end if

    call Compute_Inflow_Direction (Normal, BC_alpha, BC_beta, alpha, beta)

    T0_loc    = BC_T0
    p0_loc    = BC_p0
    alpha_loc = alpha
    beta_loc  = beta

    if (obj_rot%enabled) then
      Face_Rgas_loc = f_Rtot(BC_ci)
      call RF_Convert_BC_Inflow( Blk=Blk, Fm=Fm, Im=Im, Jm=Jm, Km=Km,     &
                                  BC_Mach_abs=0.0d0, gamma=Bound_Gamma,     &
                                  Rgas=Face_Rgas_loc, T0_or_Tstat=T0_loc,   &
                                  p0_or_pstat=p0_loc, alpha=alpha_loc,       &
                                  beta=beta_loc, Mach_rel=Mach_rel_dummy )
    end if

    Flux = 0d0
    call Stagn_Inflow ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem,      &
                               alpha_loc, beta_loc, BC_ci, BC_RANS, T0_loc, p0_loc, &
                               BC_rel_fac, BC_Sign, t_Vec, Normal, Area, Flux, BC_mdot )

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux

  end subroutine BC_Inlet_StagnCond


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 403: Inlet — total temperature T0 + prescribed mass flux g [kg/m2/s].
  ! Pressure is extrapolated from the interior.
  subroutine BC_Inlet_MassFlux_T0 ( Bm, Im, Jm, Km, Fm, Blk, BC_T0, BC_g, &
                                     BC_rel_fac, BC_alpha, BC_beta, BC_ci, BC_RANS, error )
    use FLINT_Lib_Thermodynamic
    implicit none
    type(MOSE_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Km, Fm
    real(R8), intent(in)    :: BC_T0, BC_g, BC_rel_fac, BC_alpha, BC_beta
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    integer,  intent(inout) :: error
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k, s
    real(R8) :: Normal(3), Area, t_Vec(3), BC_Sign, fmass
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, alpha, beta
    real(R8) :: E1, E2, cp3, Rgas3, T3, rho3, p3, Un3, Ut3, Ub3
    real(R8) :: Face_Prim(nprim), Face_rho, Face_Rgas, Flux(nprim)
    real(R8) :: b_Vec(3), b_Mod, XA, XB, XC, XD, XE, XF, check

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Int_k, Normal, Area, t_Vec, BC_Sign,   &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,    &
                                  Bound_Sound, Bound_Gamma, Riem )

    if (BC_Sign * Un >= 0d0) then
      error = 1
      return
    end if

    call Compute_Inflow_Direction (Normal, BC_alpha, BC_beta, alpha, beta)

    ! Vector product n x t (why?)
    b_Vec(1) = Normal(2)*t_Vec(3) - Normal(3)*t_Vec(2)
    b_Vec(2) = Normal(3)*t_Vec(1) - Normal(1)*t_Vec(3)
    b_Vec(3) = Normal(1)*t_Vec(2) - Normal(2)*t_Vec(1)
    b_Mod = norm2 ( b_Vec )
    b_Vec = b_Vec/b_Mod

    ! Obscure esotheric stuff
    XA = Normal(2) - Normal(1) * tan(alpha)
    XB =  t_Vec(2) -  t_Vec(1) * tan(alpha)
    XC =  b_Vec(2) -  b_Vec(1) * tan(alpha)
    XD = Normal(3) - Normal(1) * tan(beta)
    XE =  t_Vec(3) -  t_Vec(1) * tan(beta)
    XF =  b_Vec(3) -  b_vec(1) * tan(beta)
    check = ( XC * XE - XF * XB )
    E1 = ( XF * XA - XC * XD ) / check ! ??
    E2 = ( XD * XB - XE * XA ) / check ! ??

    ! Gas properties from species composition
    cp3   = 0d0
    Rgas3 = 0d0
    do s = 1, nsc
      cp3   = cp3   + BC_ci(s) * f_tabT(BC_T0, s, cp_tab)
      Rgas3 = Rgas3 + BC_ci(s) * Ri_tab(s)
    end do

    ! T0 is total temperature: recover static T3 iteratively
    !   T3 = T0 - u^2/(2*cp)  with u = g/rho = g*Rgas3*T3/p3
    p3  = Bound_Prim(np)   ! pressure extrapolated from interior
    T3  = (-cp3 + sqrt(cp3**2 + 2d0*cp3*BC_T0*(BC_g*Rgas3/p3)**2)) * (p3 / (BC_g*Rgas3))**2

    rho3 = p3 / (Rgas3 * T3)
    Un3  = BC_g / rho3
    Ut3  = Un3 * E1
    Ub3  = Un3 * E2

    Face_Rgas      = Rgas3
    Face_rho       = rho3
    Face_Prim(1:nsc)  = Face_rho * BC_ci
    Face_Prim(nu:nw)  = Un3 * Normal + Ut3 * t_Vec + Ub3 * b_Vec
    Face_Prim(np)     = p3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    call Compute_Flux_from_Face (Face_Prim, Face_rho, Face_Rgas, Area, Normal, Flux, fmass)
    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux

  end subroutine BC_Inlet_MassFlux_T0


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 404: Inlet — static temperature T + prescribed mass flux g [kg/m2/s].
  ! Pressure is extrapolated from the interior.
  subroutine BC_Inlet_MassFlux_T ( Bm, Im, Jm, Km, Fm, Blk, BC_T, BC_g, &
                                    BC_rel_fac, BC_alpha, BC_beta, BC_ci, BC_RANS, error )
    use FLINT_Lib_Thermodynamic
    implicit none
    type(MOSE_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Km, Fm
    real(R8), intent(in)    :: BC_T, BC_g, BC_rel_fac, BC_alpha, BC_beta
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    integer,  intent(inout) :: error
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k, s
    real(R8) :: Normal(3), Area, t_Vec(3), BC_Sign, fmass
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, alpha, beta
    real(R8) :: E1, E2, Rgas3, rho3, p3, Un3, Ut3, Ub3
    real(R8) :: Face_Prim(nprim), Face_rho, Face_Rgas, Flux(nprim)
    real(R8) :: b_Vec(3), b_Mod, XA, XB, XC, XD, XE, XF, check

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Int_k, Normal, Area, t_Vec, BC_Sign,   &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,    &
                                  Bound_Sound, Bound_Gamma, Riem )

    if (BC_Sign * Un >= 0d0) then
      error = 1
      return
    end if

    call Compute_Inflow_Direction (Normal, BC_alpha, BC_beta, alpha, beta)

    ! Vector product n x t (why?)
    b_Vec(1) = Normal(2)*t_Vec(3) - Normal(3)*t_Vec(2)
    b_Vec(2) = Normal(3)*t_Vec(1) - Normal(1)*t_Vec(3)
    b_Vec(3) = Normal(1)*t_Vec(2) - Normal(2)*t_Vec(1)
    b_Mod = norm2 ( b_Vec )
    b_Vec = b_Vec/b_Mod

    ! Obscure esotheric stuff
    XA = Normal(2) - Normal(1) * tan(alpha)
    XB =  t_Vec(2) -  t_Vec(1) * tan(alpha)
    XC =  b_Vec(2) -  b_Vec(1) * tan(alpha)
    XD = Normal(3) - Normal(1) * tan(beta)
    XE =  t_Vec(3) -  t_Vec(1) * tan(beta)
    XF =  b_Vec(3) -  b_vec(1) * tan(beta)
    check = ( XC * XE - XF * XB )
    E1 = ( XF * XA - XC * XD ) / check ! ??
    E2 = ( XD * XB - XE * XA ) / check ! ??

    Rgas3 = 0d0
    do s = 1, nsc
      Rgas3 = Rgas3 + BC_ci(s) * Ri_tab(s)
    end do

    p3   = Bound_Prim(np)      ! pressure extrapolated from interior
    rho3 = p3 / (Rgas3 * BC_T)
    Un3  = BC_g / rho3
    Ut3  = Un3 * E1
    Ub3  = Un3 * E2

    Face_Rgas         = Rgas3
    Face_rho          = rho3
    Face_Prim(1:nsc)  = Face_rho * BC_ci
    Face_Prim(nu:nw)  = Un3 * Normal + Ut3 * t_Vec + Ub3 * b_Vec
    Face_Prim(np)     = p3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    call Compute_Flux_from_Face (Face_Prim, Face_rho, Face_Rgas, Area, Normal, Flux, fmass)
    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux

  end subroutine BC_Inlet_MassFlux_T


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 406: Supersonic inlet — Mach + static temperature T and static pressure p.
  ! Note: BC_Supersonic_Inflow treats its BC_p0/BC_T0 arguments as static values,
  ! so the BC 406 static inputs map directly.  The full 3-cell Riemann stencil is
  ! preserved exactly as for BC 405.
  subroutine BC_Inlet_Supersonic_Static ( Bm, Im, Jm, Km, Fm, Blk, BC_Mach, BC_T, BC_p, &
                                           BC_rel_fac, BC_alpha, BC_beta, BC_ci, BC_RANS, BC_mdot, error )
    implicit none
    type(MOSE_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Km, Fm
    real(R8), intent(in)    :: BC_Mach, BC_T, BC_p, BC_rel_fac, BC_alpha, BC_beta
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    real(R8), intent(out)   :: BC_mdot
    integer,  intent(inout) :: error
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k
    real(R8) :: Normal(3), Area, t_Vec(3), BC_Sign
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8) :: Flux(nprim)

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Int_k, Normal, Area, t_Vec, BC_Sign,   &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,    &
                                  Bound_Sound, Bound_Gamma, Riem )

    if (BC_Sign * Un >= 0d0) then
      error = 1
      return
    end if

    ! BC_Supersonic_Inflow treats BC_p0/BC_T0 as static p/T — pass through directly.
    ! Blending with the interior and the full Riemann stencil are handled inside.
    Flux = 0d0
    call Supersonic_Inflow ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p, &
                                 BC_T, BC_ci, BC_RANS, BC_rel_fac, BC_alpha, BC_beta,         &
                                 Normal, Area, Flux, BC_mdot )

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux

  end subroutine BC_Inlet_Supersonic_Static


  !! Low-level routines to be called from the high-level BC subroutines above.  These work on primitive arrays directly and are not aware of the block structure.


  !─────────────────────────────────────────────────────────────────────────────
  ! Supersonic_Inflow (low-level, works on primitive arrays directly)
  subroutine Supersonic_Inflow ( Prim2, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p0, BC_T0, BC_ci, & 
                                    BC_RANS, BC_rel_fac, BC_alpha, BC_beta, Normal, Area, Flux, Fmass )
    use MOSE_Lib_Limiters, only : rlimiter
    use MOSE_Mod_Riemann
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in) :: modfm, modfm1, modfm2
    real(R8), intent(in) :: Prim2(nprim), Int_Prim(nprim), BC_Mach, BC_p0, BC_T0, BC_rel_fac
    real(R8), intent(in) :: BC_ci(nsc), BC_RANS(1:), BC_alpha, BC_beta, Normal(3), Area
    real(R8), intent(out) :: Flux(nprim), Fmass
    ! Local
    integer :: s
    real(R8) :: Sup_Prim(nprim), Sup_T, Sup_Rgas, Sup_rho, Sup_Sound, Prim1(nprim), Prim3(nprim)
    real(R8), dimension(nprim) :: Diff32, Diff21, Slope, Face_Prim, Prim_L, Prim_R
    real(R8) :: rho_L, rho_R, Rgas_L, Rgas_R, Sound_L, Sound_R, F_r, F_u, F_v, F_w, F_E
    real(R8) :: su, Sel_L, Sel_R
    real(R8) :: flow_Rgas, flow_Sound, flow_Mach, BC_Mach_local, BC_p0_local
    real(R8) :: alpha_eff, beta_eff

    ! compute flow properties in the interior cell
    flow_Rgas = f_Rtot( Prim2(1:nsc) )
    flow_Sound = f_ss ( Prim2(1:nsc), Prim2(np), sum(Prim2(1:nsc)), flow_Rgas )
    flow_Mach = norm2( Prim2(nu:nw) ) / flow_Sound

    ! blended Mach number at the boundary
    BC_Mach_local = BC_rel_fac*BC_Mach + (1-BC_rel_fac)*flow_Mach
    ! blended p static at the boundary
    BC_p0_local = BC_rel_fac*BC_p0 + (1-BC_rel_fac)*Prim2(np)

    call Compute_Inflow_Direction (Normal, BC_alpha, BC_beta, alpha_eff, beta_eff)

    ! supersonic inflow. BC_Mach enforced at boundary, T0 and p0 are static
    Sup_Prim(np) = BC_p0_local
    Sup_T = BC_T0
    Sup_Rgas = sum ( BC_ci*Ri_Tab )
    Sup_rho = Sup_Prim(np) / ( Sup_Rgas * Sup_T )
    Sup_Prim(1:nsc) = Sup_rho * BC_ci
    Sup_Sound = f_ss ( Sup_Prim(1:nsc), Sup_Prim(np), Sup_rho, Sup_Rgas )
    Sup_Prim(nu) = BC_Mach_local * Sup_Sound * cos(alpha_eff) * cos(beta_eff)
    Sup_Prim(nv) = BC_Mach_local * Sup_Sound * sin(alpha_eff) * cos(beta_eff)
    Sup_Prim(nw) = BC_Mach_local * Sup_Sound * sin(beta_eff)
    if (model==2) Sup_Prim(nt:nprim) = Sup_rho * BC_RANS

    ! building a 3 cell stencil to extrapolate the solution at the interface

    ! stencil cell 2: boundary cell

    ! stencil cell 1 is the ghost cell (odd faces) or the interior cell (even faces)
    Prim1 = modfm*Sup_Prim + modfm1*Int_Prim

    ! stencil cell 3 is the interior cell (odd faces) or the ghost cell (even faces)
    Prim3 = modfm*Int_Prim + modfm1*Sup_Prim
    
    ! The state at the interface from the interior side is reconstructed: (ro,u,p)
    Diff32 = Prim3 - Prim2
    Diff21 = Prim2 - Prim1
    do s = 1, nprim
      Slope(s) = rlimiter ( Diff32(s), Diff21(s) )
    end do
    
    ! Extrapolation at the interface from the boundary cell
    Face_Prim = Prim2 + 0.5d0 * Slope * modfm2

    ! state L of riemann problem is the ghost state for even faces and boundary reconstructed state for odd faces
    Prim_L = modfm*Sup_Prim + modfm1*Face_Prim
    call co_rotot_Rtot ( Prim_L(1:nsc), rho_L, Rgas_L )
    Sound_L = f_ss ( Prim_L(1:nsc), Prim_L(np), rho_L, Rgas_L )

    ! state R of riemann problem is the boundary reconstructed state for odd faces and the ghost state for even faces
    Prim_R = modfm*Face_Prim + modfm1*Sup_Prim
    call co_rotot_Rtot ( Prim_R(1:nsc), rho_R, Rgas_R )
    Sound_R = f_ss ( Prim_R(1:nsc), Prim_R(np), rho_R, Rgas_R )

    call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(np), Sound_L, rho_L, &
                   Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(np), Sound_R, rho_R, &
                   1d0, Normal(1), Normal(2), Normal(3), F_r, F_u, F_v, F_w, F_E)

    su = sign ( 1d0, F_r )
    Sel_L = ( 1d0+su ) / 2d0  ! state left selector
    Sel_R = ( su-1d0 ) / 2d0  ! state right selector

    ! Fluxes
    Fmass = F_r * Area
    Flux = 0.0d0
    Flux(1:nsc) = Fmass * ( Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R )
    Flux(nu) = F_u * Area
    Flux(nv) = F_v * Area
    Flux(nw) = F_w * Area
    Flux(np) = F_E * Area

    if (model==2) &
      Flux(nt:nprim) = Fmass * ( Sel_L*Prim_L(nt:nprim)/rho_L - Sel_R*Prim_R(nt:nprim)/rho_R )

  end subroutine Supersonic_Inflow

  !─────────────────────────────────────────────────────────────────────────────
  ! Subsonic inflow: prescribed stagnation conditions (T0, p0).
  ! Solves the characteristic equations to find the face state compatible with
  ! the interior Riemann invariant and the prescribed total conditions.
  subroutine Stagn_Inflow ( Prim, Un, Sound, Gamma, Riem, alpha, beta, &
                                 BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, &
                                 t_Vec, Normal, Area, Flux, Fmass )
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: Prim(nprim), Un, Sound, Gamma, Riem, alpha, beta
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(1:), BC_T0, BC_p0, BC_rel_fac
    real(R8), intent(in)  :: BC_Sign, t_Vec(3), Normal(3), Area
    real(R8), intent(out) :: Flux(nprim), Fmass
    ! Local
    integer  :: s, iter_discr
    real(R8) :: Face_Rgas, Sound0, Coef_D, Coef_C, YA, YB, YC, Discr, Un_Sign, Un3, Sound3
    real(R8) :: p3, Ut3, Ub3, BC_p0_try
    real(R8) :: Face_T, Face_Sound, Face_Prim(nprim), Face_rho, Face_Enthalpy, Face_Un
    real(R8) :: XA, XB, XC, XD, XE, XF, check, E1, E2, b_Vec(3), b_Mod

    ! Vector product n x t (why?)
    b_Vec(1) = Normal(2)*t_Vec(3) - Normal(3)*t_Vec(2)
    b_Vec(2) = Normal(3)*t_Vec(1) - Normal(1)*t_Vec(3)
    b_Vec(3) = Normal(1)*t_Vec(2) - Normal(2)*t_Vec(1)
    b_Mod = norm2 ( b_Vec )
    b_Vec = b_Vec/b_Mod

    ! Obscure esotheric stuff
    XA = Normal(2) - Normal(1) * tan(alpha)
    XB =  t_Vec(2) -  t_Vec(1) * tan(alpha)
    XC =  b_Vec(2) -  b_Vec(1) * tan(alpha)
    XD = Normal(3) - Normal(1) * tan(beta)
    XE =  t_Vec(3) -  t_Vec(1) * tan(beta)
    XF =  b_Vec(3) -  b_vec(1) * tan(beta)
    check = ( XC * XE - XF * XB )
    E1 = ( XF * XA - XC * XD ) / check ! ??
    E2 = ( XD * XB - XE * XA ) / check ! ??

    discr = -1.d0
    BC_p0_try = BC_p0
    iter_discr = 0
    do while (discr<0)
      iter_discr = iter_discr + 1
      ! Mach in BC input >0: Subsonic inflow with assigned total pressure and temperature (po To)
      Face_Rgas = sum ( BC_ci * Ri_tab )
      Sound0 = sqrt ( Gamma * Face_Rgas * BC_T0 )
      Coef_D = ( 1d0 + E1**2 + E2**2 ) * ( Gamma - 1d0 )/2d0
      Coef_C = ( BC_p0_try / Prim(np) )**( ( Gamma - 1d0 ) / Gamma )
      Coef_C = Coef_C * ( ( Sound / Sound0 )**2 )
      YA = ( ( Gamma - 1d0 ) / 2d0 )**2 + Coef_C * Coef_D
      YB = BC_Sign * Riem * ( ( Gamma - 1d0 ) / 2d0 )
      YC = Riem**2 - Coef_C * Sound0**2
      Discr = YB**2 - YA * YC
      if ( ( discr < 0d0 ) .and. ( abs(discr) < 1d-2 ) ) Discr = 0d0
      if ( ( discr < 0d0 ) .and. ( abs(discr) > 1d-2) ) then
       BC_p0_try = BC_rel_fac*BC_p0_try + (1-BC_rel_fac)*Prim(np)
      end if
      if ( iter_discr > 1000 ) then
        Discr = 0d0
        exit
      end if
    enddo

    Un_Sign = sign ( 1d0, Un )  ! Sign of normal velocity at the boundary interface
    Un3 = ( YB + Un_Sign * sqrt( Discr ) ) / YA
    Sound3 = sqrt ( Sound0**2 - Coef_D * Un3**2 )
    p3 = BC_p0_try * ( Sound3 / Sound0 )**( 2d0 * Gamma / ( Gamma - 1d0 ) )
    Ut3 = Un3 * E1
    Ub3 = Un3 * E2
    Face_Sound = Sound3
    Face_Prim(np) = p3
    Face_rho = Gamma * Face_Prim(np) / Face_Sound**2
    Face_Prim(1:nsc) = BC_ci * Face_rho
    Face_Prim(nu:nw) = Un3 * Normal + Ut3 * t_Vec + Ub3 * b_Vec
    Face_Un = Un3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    Face_Enthalpy = 0d0
    Face_T = Face_Prim(np) / ( Face_Rgas * Face_rho )
    Flux = 0.0d0
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT( Face_T,s,h_tab )
      Flux(s) = Face_Prim(s) * Face_Un * Area
    end do
    Fmass = Face_rho * Face_Un * Area
    Flux(nu:nw) = Fmass * Face_Prim(nu:nw) + Face_Prim(np) * Area * Normal
    Flux(np) = Fmass * ( 0.5d0 * sum( Face_Prim(nu:nw)**2 ) + Face_Enthalpy )
    if (model==2) Flux(nt:nprim) = Face_Prim(nt:nprim) / Face_rho * Fmass 

  end subroutine Stagn_Inflow


  !! Common geometric/metric computations for every inflow BC. Called at the start of each per-type routine.


  pure subroutine Compute_Inflow_Direction(normal, alpha_in, beta_in, alpha_out, beta_out)
    real(R8), intent(in) :: normal(3), alpha_in, beta_in
    real(R8), intent(out) :: alpha_out, beta_out
    real(R8) :: alpha_rad, beta_rad, cos_alpha, cos_beta

    ! Resolve 'normal' direction sentinel (parse_dir_tok returns huge(R8) for token 'normal').
    if (alpha_in >= 0.5_R8 * huge(1.0_R8)) then
      alpha_out = atan(Normal(2) / (Normal(1) + 1d-20))
    else
      alpha_out = alpha_in
    end if
    if (beta_in >= 0.5_R8 * huge(1.0_R8)) then
      beta_out  = atan(Normal(3) / (Normal(1) + 1d-20))
    else
      beta_out  = beta_in
    end if

  end subroutine Compute_Inflow_Direction


  !─────────────────────────────────────────────────────────────────────────────
  ! Setup_Inflow_Geometry: compute all geometric/metric quantities common to
  ! every inflow BC.  Called at the start of each per-type routine.
  pure subroutine Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm,               &
                                     modfm, modfm1, modfm2, modfm3,      &
                                     Int_i, Int_j, Int_k,                &
                                     Normal, Area, t_Vec, BC_Sign,       &
                                     Bound_Prim, Int_Prim,               &
                                     Un, Bound_rho, Bound_Rgas,          &
                                     Bound_Sound, Bound_Gamma, Riem )
    use FLINT_Lib_Thermodynamic
    implicit none
    type(MOSE_block_type), intent(in)  :: Blk
    integer,  intent(in)               :: Im, Jm, Km, Fm
    integer,  intent(out)              :: modfm, modfm1, modfm2, modfm3
    integer,  intent(out)              :: Int_i, Int_j, Int_k
    real(R8), intent(out)              :: Normal(3), Area, t_Vec(3), BC_Sign
    real(R8), intent(out)              :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8), intent(out)              :: Un, Bound_rho, Bound_Rgas
    real(R8), intent(out)              :: Bound_Sound, Bound_Gamma, Riem
    ! Local
    integer  :: Dir, Face_i, Face_j, Face_k
    real(R8) :: t_Mod

    call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )

    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)
    Int_k = Km + guide(Fm,3)

    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    Bound_Prim = Blk % P(:,Im,Jm,Km)
    Int_Prim   = Blk % P(:,Int_i,Int_j,Int_k)

    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area   = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    select case ( Fm )
      case (1,2)
        t_Vec = Blk % node(Face_i,Face_j,Face_k) % c - &
                Blk % node(Face_i,Face_j-1,Face_k-1) % c
      case (3,4)
        t_Vec = Blk % node(Face_i,Face_j,Face_k) % c - &
                Blk % node(Face_i-1,Face_j,Face_k-1) % c
      case (5,6)
        t_Vec = Blk % node(Face_i,Face_j,Face_k) % c - &
                Blk % node(Face_i-1,Face_j-1,Face_k) % c
    end select
    t_Mod = norm2(t_Vec)
    t_Vec = t_Vec / t_Mod

    BC_Sign = real(modfm2)   ! -1 for faces 1,3,5 ; +1 for faces 2,4,6

    Un          = dot_product(Bound_Prim(nu:nw), Normal)
    call co_rotot_Rtot(Bound_Prim(1:nsc), Bound_rho, Bound_Rgas)
    Bound_Sound = f_ss   (Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas)
    Bound_Gamma = f_gamma(Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas)
    Riem        = Bound_Sound + BC_Sign * 0.5d0 * (Bound_Gamma - 1d0) * Un

  end subroutine Setup_Inflow_Geometry


  !─────────────────────────────────────────────────────────────────────────────
  ! Compute_Flux_from_Face: assemble the convective flux vector from a known
  ! face state.  Used by all inflow/outflow routines.
  pure subroutine Compute_Flux_from_Face ( Face_Prim, Face_rho, Face_Rgas, Area, Normal, Flux, Fmass )
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: Face_Prim(nprim), Face_rho, Face_Rgas, Area, Normal(3)
    real(R8), intent(out) :: Flux(nprim), Fmass
    ! Local
    integer  :: s
    real(R8) :: Face_Un, Face_T, Face_Enthalpy

    Face_Un       = dot_product(Face_Prim(nu:nw), Normal)
    Face_T        = Face_Prim(np) / (Face_Rgas * Face_rho)
    Fmass         = Face_rho * Face_Un * Area
    Face_Enthalpy = 0d0
    Flux          = 0d0
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s)/Face_rho * f_tabT(Face_T,s,h_tab)
      Flux(s)       = Face_Prim(s) * Face_Un * Area
    end do
    Flux(nu:nw) = Fmass * Face_Prim(nu:nw) + Face_Prim(np) * Area * Normal
    Flux(np)    = Fmass * (0.5d0 * sum(Face_Prim(nu:nw)**2) + Face_Enthalpy)
    if (nprim > np) Flux(np+1:nprim) = Face_Prim(np+1:nprim) / Face_rho * Fmass

  end subroutine Compute_Flux_from_Face


end module MOSE_Lib_BC_Fluxes_Inflow