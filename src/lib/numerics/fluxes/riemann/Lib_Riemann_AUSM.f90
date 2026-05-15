!  Key References for AUSM-family Schemes

!- AUSM (Advection Upstream Splitting Method)
!
!  Liou, M.-S., and Steffen Jr., C. J. (1993)
!  A New Flux Splitting Scheme. Journal of Computational Physics, 107, 23–39. DOI: 10.1006/jcph.1993.1132
!
!  Summary:
!  Introduced the AUSM concept — splitting inviscid flux into convective (Mach-number-based) and pressure (acoustic) components.
!  Foundation of all subsequent AUSM-family methods.

!- AUSM⁺ (Enhanced AUSM)
!
!  Liou, M.-S. (1996)
!  A Sequel to AUSM: AUSM⁺. Journal of Computational Physics, 129, 364–382. DOI: 10.1006/jcph.1996.0256
!
!  Summary:
!  Refined the Mach number and pressure splitting for smoother behavior near ∣𝑀∣=1, improved stability, and better low-Mach performance.
!  Standard for compressible flow solvers.

!- AUSM⁺-up (All-Speed Variant)
!
!  Liou, M.-S. (2006)
!  A Sequel to AUSM, Part II: AUSM⁺-up for All Speeds. Journal of Computational Physics, 214, 137–170. DOI: 10.1016/j.jcp.2005.09.020
!
!  Summary:
!  Introduced a pressure-based velocity dissipation term (“up”) for all-speed flows; very effective for subsonic/incompressible regimes.

!- AUSM⁺-up2 (All-Speed Variant for Hypersonics)
!
!  Kitamura, K., & Shima, E.
!  Towards shock-stable and accurate hypersonic heating computations: A new pressure flux for AUSM-family schemes.
!
!  Summary:
!  An improved AUSM+-up version to take into account dissipation scaling with Mach number.

module MOSE_Lib_Riemann_AUSM
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Config_Types_m, only: obj_riemann

  implicit none
  private
  public :: riemann_Hanel, riemann_AUSMp
  public :: riemann_AUSMp_up, riemann_AUSMp_up2

contains

  subroutine riemann_Hanel(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,switch,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0, f_gamma, f_Rtot
    implicit none
    ! --- inputs ---
    real(R8), intent(in)  :: dl(nsc), ul, vl, wl, pl, al, dltot
    real(R8), intent(in)  :: dr(nsc), ur, vr, wr, pr, ar, drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: switch
    ! --- outputs ---
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_E
    ! --- local ---
    real(R8) :: unl, unr
    real(R8) :: H_L, H_R
    real(R8) :: uLp, uRm
    real(R8) :: pLp, pRm, phalf

    !====================================================
    ! Normal velocities and densities
    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz

    !====================================================
    ! Thermodynamics
    H_L = H0(pl, dl, sqrt(ul*ul + vl*vl + wl*wl))
    H_R = H0(pr, dr, sqrt(ur*ur + vr*vr + wr*wr))

    !====================================================
    ! Velocity splitting
    if (abs(unl) <= al) then
      uLp = (unl + al)**2/(4.d0*al)
    else
      uLp = 0.5d0*(unl + abs(unl))
    endif

    if (abs(unr) <= ar) then
      uRm = (-(unr - ar)**2)/(4.d0*ar)
    else
      uRm = 0.5d0*(unr - abs(unr))
    endif

    !====================================================
    ! Pressure splitting
    if (abs(unl) <= al) then
      pLp = pl*uLp*(2d0 - unl/al)/al
    else
      pLp = pl*uLp/unl
    endif

    if (abs(unr) <= ar) then
      pRm = pr*uRm*(-2d0 - unr/ar)/ar
    else
      pRm = pr*uRm/unr
    endif

    phalf = pLp + pRm

    !====================================================
    ! Final fluxes
    F_r = uLp*dltot + uRm*drtot

    F_u = uLp*(dltot*ul) + uRm*(drtot*ur) + phalf*nx
    F_v = uLp*(dltot*vl) + uRm*(drtot*vr) + phalf*ny
    F_w = uLp*(dltot*wl) + uRm*(drtot*wr) + phalf*nz

    F_E = uLp*(dltot*H_L) + uRm*(drtot*H_R)

  end subroutine riemann_Hanel


  !> @brief AUSM⁺ Riemann solver for compressible Euler equations.
  !>
  !> Computes interface fluxes using Liou’s AUSM⁺ scheme, which separates
  !> convective (Mach-number-based) and pressure (acoustic) flux components.
  !> It ensures smooth transitions between subsonic and supersonic regimes
  !> and suppresses pressure oscillations near shocks.
  !>
  !> @details
  !> The AUSM⁺ method modifies the original AUSM by introducing improved
  !> Mach number and pressure splitting functions:
  !>   - Quadratic splitting near |M| < 1 for continuity,
  !>   - Pressure-weighted pressure splitting for smoother coupling.
  !> This version provides reliable shock capturing and low dissipation.
  !>
  !> @references
  !>   - Liou, M.-S., and Steffen, C. J. Jr. (1993),
  !>     "A New Flux Splitting Scheme," Journal of Computational Physics, 107, 23–39.
  !>   - Liou, M.-S. (1996),
  !>     "A Sequel to AUSM: AUSM⁺," Journal of Computational Physics, 129, 364–382.
  !>
  !> @param[in]  dl(nsc)   Species densities (left state)
  !> @param[in]  ul,vl,wl  Velocity components (left state)
  !> @param[in]  pl,al     Pressure and sound speed (left state)
  !> @param[in]  dltot     Total density (left state)
  !> @param[in]  dr(nsc)   Species densities (right state)
  !> @param[in]  ur,vr,wr  Velocity components (right state)
  !> @param[in]  pr,ar     Pressure and sound speed (right state)
  !> @param[in]  drtot     Total density (right state)
  !> @param[in]  nx,ny,nz  Unit normal components
  !> @param[out] F_r,F_u,F_v,F_w,F_E Fluxes of mass, momentum, and energy
  !>
  !> @ingroup Lib_RiemannPrivateProcedure
  subroutine riemann_AUSMp(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,switch,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0, f_gamma, f_Rtot
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: switch
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! specific
    real(R8), parameter :: beta = 0.125d0, alpha = 0.1875d0
    real(R8) :: unl, unr 
    real(R8) :: ML, MR
    real(R8) :: MplusL, MminusR, Mhalf
    real(R8) :: pplusL, pminusR, phalf
    real(R8) :: H_L, H_R, gamma_L, gamma_R, Rgasl, Rgasr
    real(R8) :: a_half, mdotm, mdotp
    real(R8) :: astarL2, astarR2, astarL, astarR, atildeL, atildeR

    ! normal velocities
    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz

    ! total enthalpies
    H_L = H0(pl, dl, sqrt(ul**2 + vl**2 + wl**2))
    H_R = H0(pr, dr, sqrt(ur**2 + vr**2 + wr**2))

    ! --- Numerical interface sound speed ---
    ! ! Liou's strict AUSM+
    ! Rgasl = f_Rtot(dl)
    ! Rgasr = f_Rtot(dr)
    ! gamma_L = f_gamma(dl, pl, dltot, Rgasl)
    ! gamma_R = f_gamma(dr, pr, drtot, Rgasr)
    ! ! critical speed squared a*^2 = 2*(gamma-1)*H/(gamma+1)
    ! astarL2 = 2d0*(gamma_L - 1d0)/(gamma_L + 1.d0)*H_L
    ! astarR2 = 2d0*(gamma_R - 1d0)/(gamma_R + 1.d0)*H_R
    ! astarL = sqrt(astarL2)
    ! astarR = sqrt(astarR2)
    ! ! atilde = a*^2 / max(a*, |u|)
    ! atildeL = astarL2 / max(astarL, abs(unl))
    ! atildeR = astarR2 / max(astarR, abs(unr))
    ! ! interface numerical sound speed
    ! a_half = min(atildeL, atildeR)

    ! simplified AUSM+ sound speed (Liou, PROGRESS TOWARDS AN IMPROVED CFD METHOD: AUSM+)
    !a_half = 0.5d0 * (al + ar)
    a_half = sqrt( al * ar )

    ! --- Interface Mach numbers use a_half ---
    ML = unl / a_half
    MR = unr / a_half

    ! --- AUSM+ Mach/pressure splitting polynomials ---
    if (abs(ML) < 1.d0) then
      MplusL = 0.25d0*(ML + 1.d0)**2 + beta*(ML**2 - 1.d0)**2
      pplusL = 0.25d0*(ML + 1.d0)**2 * (2.d0 - ML) + alpha*ML*(ML**2 - 1.d0)**2
    else
      MplusL = 0.5d0*(ML + abs(ML))
      pplusL = 0.5d0*(1.d0 + sign(1.d0, ML))
    endif

    if (abs(MR) < 1.d0) then
      MminusR = -0.25d0*(MR - 1.d0)**2 - beta*(MR**2 - 1.d0)**2
      pminusR = 0.25d0*(MR - 1.d0)**2 * (2.d0 + MR) - alpha*MR*(MR**2 - 1.d0)**2
    else
      MminusR = 0.5d0*(MR - abs(MR))
      pminusR = 0.5d0*(1.d0 - sign(1.d0, MR))
    endif

    ! combined interface quantities
    Mhalf = MplusL + MminusR
    phalf = pplusL*pl + pminusR*pr

    ! mass fluxes
    mdotp = dltot * a_half * 0.5d0 * (Mhalf + abs(Mhalf)) ! max(0d0,Mhalf)
    mdotm = drtot * a_half * 0.5d0 * (Mhalf - abs(Mhalf)) ! min(0d0,Mhalf)

    ! Mass
    F_r = mdotp + mdotm

    ! Momentum
    F_u = mdotp * ul + mdotm * ur + nx * phalf
    F_v = mdotp * vl + mdotm * vr + ny * phalf
    F_w = mdotp * wl + mdotm * wr + nz * phalf

    ! Energy
    F_E = mdotp * H_L + mdotm * H_R

  end subroutine riemann_AUSMp


  !> @brief AUSM⁺-up Riemann solver for compressible Euler equations (SU2-style)
  !>
  !> Computes interface fluxes using Liou’s AUSM⁺ augmented with the AUSM⁺-up
  !> pressure/velocity diffusion terms (pressure-diffusion in mass flux and
  !> velocity-diffusion in pressure flux). This implementation follows the
  !> SU2 numerical realization.
  !>
  !> @details
  !> - Uses Liou's a* / \tilde a numerical interface sound-speed and forms
  !>   interface Mach numbers with aF (the interface sound speed).
  !> - Mach and pressure splitting polynomials use SU2 dynamic alpha/beta
  !>   scaling (fa-based) for smooth blending across |M|=1.
  !> - Adds AUSM⁺-up diffusion terms:
  !>     Mp = -(Kp/fa) * phi(Mbar) * (pR - pL) / (rho_half * aF^2)
  !>     Pu = -Ku*fa*betaL*betaR*2*rho_half*aF*(uR - uL)
  !>   with phi(x) = max(1 - sigma*x^2, 0).
  !> - Mass flux uses mdot = aF*( max(mF,0)*rhoL + min(mF,0)*rhoR )
  !> - Momentum: mdot * u_upwind + p_half * n
  !> - Energy: mdot * H_upwind (H upwind is total enthalpy per unit mass)
  !>
  !> @references
  !>   - Liou, M.-S. (2006), "A sequel to AUSM, Part II: AUSM⁺-up for all speeds", J. Comput. Phys., 214, 137–170.
  !>   - SU2 open-source implementation (AUSM⁺-up numerical realization; used as template).
  !>
  !> @param[in]  dl(nsc)   Species densities (left state)
  !> @param[in]  ul,vl,wl  Velocity components (left state)
  !> @param[in]  pl,al     Pressure and sound speed (left state)
  !> @param[in]  dltot     Total density (left state)
  !> @param[in]  dr(nsc)   Species densities (right state)
  !> @param[in]  ur,vr,wr  Velocity components (right state)
  !> @param[in]  pr,ar     Pressure and sound speed (right state)
  !> @param[in]  drtot     Total density (right state)
  !> @param[in]  nx,ny,nz  Unit normal components
  !> @param[out] F_r,F_u,F_v,F_w,F_E Fluxes of mass, momentum, and energy
  !>
  !> @ingroup Lib_RiemannPrivateProcedure
  subroutine riemann_AUSMp_up(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,switch,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0, f_gamma, f_Rtot
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: switch
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! AUSM+-up parameters
    real(R8), parameter :: Kp = 0.25d0
    real(R8), parameter :: Ku = 0.75d0
    real(R8), parameter :: sigma = 1.0d0
    real(R8), parameter :: beta =  0.125d0
    real(R8) :: unl, unr
    real(R8) :: H_L, H_R
    real(R8) :: gamma_L, gamma_R, Rgasl, Rgasr
    real(R8) :: astarL, astarR, astarL2, astarR2
    real(R8) :: atildeL, atildeR, aF
    real(R8) :: mL, mR, MFsq, Mrefsq, fa, alpha
    real(R8) :: mLP, mRM, betaLP, betaRM
    real(R8) :: rhoF, Mp, Pu
    real(R8) :: mF, mdot
    real(R8) :: upwind_u, upwind_v, upwind_w
    real(R8) :: upwind_H
    real(R8) :: phalf, p1, p2

    ! --- projected normal velocities ---
    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz

    ! total enthalpies
    H_L = H0(pl, dl, sqrt(ul**2 + vl**2 + wl**2))
    H_R = H0(pr, dr, sqrt(ur**2 + vr**2 + wr**2))

    ! --- Liou's numerical interface sound speed ---
    Rgasl = f_Rtot(dl)
    Rgasr = f_Rtot(dr)
    gamma_L = f_gamma(dl, pl, dltot, Rgasl)
    gamma_R = f_gamma(dr, pr, drtot, Rgasr)
    astarL2 = 2d0*(gamma_L - 1d0)/(gamma_L + 1.d0)*H_L
    astarR2 = 2d0*(gamma_R - 1d0)/(gamma_R + 1.d0)*H_R
    astarL = sqrt(astarL2)
    astarR = sqrt(astarR2)

    ! atilde = a*^2 / max(a*, |u|)
    atildeL = astarL2 / max(astarL,  unl)
    atildeR = astarR2 / max(astarR, -unr)

    aF = min(atildeL, atildeR)

    ! --- dimensionless Machs at interface ---
    mL = unl / aF
    mR = unr / aF

    ! --- Mbar related quantities
    MFsq = 0.5d0*(mL*mL + mR*mR)
    Mrefsq = min(1.0d0, max(MFsq, obj_riemann%Minf*obj_riemann%Minf))

    fa = 2.d0*sqrt(Mrefsq) - Mrefsq

    alpha = 0.1875d0 * (-4.d0 + 5.d0 * fa * fa)

    ! --- compute split polynomials ---
    if (abs(mL) <= 1.d0) then
      p1 = 0.25d0 * (mL + 1.d0)**2
      p2 = (mL*mL - 1.d0)**2
      mLP    = p1 + beta * p2
      betaLP = p1 * (2.d0 - mL) + alpha * mL * p2
    else
      mLP    = 0.5d0 * (mL + abs(mL))
      betaLP = mLP/mL
    endif

    if (abs(mR) <= 1.d0) then
      p1 = 0.25d0 * (mR - 1.d0)**2
      p2 = (mR*mR - 1.d0)**2
      mRM    = -p1 - beta * p2
      betaRM = p1 * (2.d0 + mR) - alpha * mR * p2
    else
      mRM    = 0.5d0 * (mR - abs(mR))
      betaRM = mRM/mR
    endif

    ! --- pressure and mass diffusion terms (SU2 forms) ---
    rhoF = 0.5d0 * (dltot + drtot)

    Mp = - (Kp / fa) * max(0.d0, 1.d0 - sigma * MFsq) * (pr - pl) / (rhoF * aF * aF)

    Pu = - Ku * fa * betaLP * betaRM * 2.d0 * rhoF * aF * (unr - unl)

    ! --- combined mF, mdot ---
    mF = mLP + mRM + Mp
    mdot = aF * ( max(mF, 0.d0) * dltot + min(mF, 0.d0) * drtot )

    ! pressure flux
    phalf = betaLP * pl + betaRM * pr + Pu

    ! --- Momentum and energy fluxes: use mdot upwinding + pressure acoustic term ---
    if (mF >= 0.d0) then
      ! upwind left
      upwind_u = ul
      upwind_v = vl
      upwind_w = wl
      upwind_H = H_L
    else
      upwind_u = ur
      upwind_v = vr
      upwind_w = wr
      upwind_H = H_R
    endif

    ! mass flux output
    F_r = mdot

    ! momentum flux: mdot * upwind velocity + pressure * n
    F_u = mdot * upwind_u + nx * phalf
    F_v = mdot * upwind_v + ny * phalf
    F_w = mdot * upwind_w + nz * phalf

    ! energy flux: mdot * upwind enthalpy (total)
    F_E = mdot * upwind_H

  end subroutine riemann_AUSMp_up


  !> @brief AUSM⁺-up2 Riemann solver (SU2 / Kitamura & Shima pressure flux)
  !>
  !> Computes interface fluxes for the compressible Euler equations using an
  !> AUSM-family flux with AUSM⁺ splitting and the AUSM⁺-up2 pressure flux
  !> modification used in SU2.  This variant implements:
  !>   - Liou-style AUSM⁺ Mach/pressure split functions and Liou's numerical
  !>     interface sound-speed a* (used to form the interface aF),
  !>   - a pressure-diffusion correction in the mass flux (AUSM⁺-up style),
  !>   - the modified pressure-flux formulation proposed by Kitamura & Shima
  !>     (and adopted in SU2's AUSM⁺-up2), which couples the split pressure
  !>     polynomials with the local velocity magnitude to improve shock
  !>     stability and hypersonic heating accuracy.
  !>
  !> @details
  !> The implementation follows the SU2 AUSM⁺-up2 numeric realization:
  !> 1. Compute Liou's critical speed a* (per state) and the numerical
  !>    interface sound speed aF = min(â_L, â_R), with â = a*^2 / max(a*, |u|).
  !> 2. Form interface Mach numbers m_L = u_L / aF and m_R = u_R / aF.
  !> 3. Evaluate dynamic split-polynomials (M^+, M^-, and pressure-coeffs).
  !> 4. Add a pressure-diffusion term Mp to the dimensionless mass coefficient:
  !>       Mp = -(Kp/fa) * φ(Mbar) * (p_R - p_L) / (rho_half * aF^2),
  !>    with φ(x) = max(1 - sigma * x^2, 0) (low-M limiter).
  !> 5. Build the AUSM⁺-up2 pressure flux following Kitamura & Shima:
  !>       p_half = ½(p_L + p_R) + ½(p⁺_L - p⁻_R)(p_L - p_R)
  !>                + sqrt(½(|v|^2_L + |v|^2_R)) (p⁺_L + p⁻_R - 1) ρ_half aF
  !>    This term couples a velocity-scale term (√|v|²) with the split
  !>    pressure coefficients to improve shock stability and heat-flux
  !>    predictions in hypersonic flows.
  !> 6. Mass flux: mdot = aF * [ max(mF,0)*ρ_L + min(mF,0)*ρ_R ] with mF = M^+_L + M^-_R + Mp.
  !> 7. Momentum: mdot * u_upwind + p_half * n ; Energy: mdot * H_upwind.
  !>
  !> @tunable_parameters
  !>  - Kp   : pressure-diffusion coefficient (default 0.25)
  !>  - sigma: limiter parameter for φ (default 1.0)
  !>  - Minf : reference Mach used to define the fa limiter floor (small positive)
  !>
  !> @references
  !>  - Kitamura, K., & Shima, E., "Towards shock-stable and accurate hypersonic
  !>    heating computations: A new pressure flux for AUSM-family schemes."
  !>  - SU2 open-source implementation: AUSM⁺-up2 numerical realization (for
  !>    practical parameter choices and polynomial implementations).
  !>
  !> @param[in]  dl(nsc)   Species densities (left state)
  !> @param[in]  ul,vl,wl  Velocity components (left state)
  !> @param[in]  pl,al     Pressure and sound speed (left state)
  !> @param[in]  dltot     Total density (left state)
  !> @param[in]  dr(nsc)   Species densities (right state)
  !> @param[in]  ur,vr,wr  Velocity components (right state)
  !> @param[in]  pr,ar     Pressure and sound speed (right state)
  !> @param[in]  drtot     Total density (right state)
  !> @param[in]  nx,ny,nz  Unit normal components
  !> @param[out] F_r,F_u,F_v,F_w,F_E Fluxes of mass, momentum, and energy
  !>
  !> @ingroup Lib_RiemannPrivateProcedure
  subroutine riemann_AUSMp_up2(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,switch,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0, f_gamma, f_Rtot
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: switch
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! parameters (SU2 defaults)
    real(R8), parameter :: Kp = 0.25d0
    real(R8), parameter :: sigma = 1.0d0
    real(R8), parameter :: beta = 0.125d0
    real(R8) :: unl, unr, sq_vel
    real(R8) :: H_L, H_R
    real(R8) :: gamma_L, gamma_R, Rgasl, Rgasr
    real(R8) :: astarL2, astarR2, astarL, astarR
    real(R8) :: ahatL, ahatR, aF
    real(R8) :: mL, mR, MFsq, param1, Mrefsq, fa, alpha
    real(R8) :: mLP, mRM, pLP, pRM
    real(R8) :: rhoF, Mp, mF, mdot, phalf
    real(R8) :: upwind_u, upwind_v, upwind_w, upwind_H
    real(R8) :: p1, p2

    ! --- projected velocities and mean-squared magnitude ---
    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz
    sq_vel = 0.5d0*((ul**2 + vl**2 + wl**2) + (ur**2 + vr**2 + wr**2))

    ! --- total enthalpy (per unit mass) ---
    H_L = H0(pl, dl, sqrt(ul**2 + vl**2 + wl**2))
    H_R = H0(pr, dr, sqrt(ur**2 + vr**2 + wr**2))

    ! --- gamma per state ---
    Rgasl = f_Rtot(dl)
    Rgasr = f_Rtot(dr)
    gamma_L = f_gamma(dl, pl, dltot, Rgasl)
    gamma_R = f_gamma(dr, pr, drtot, Rgasr)
    astarL2 = 2d0*(gamma_L - 1d0)/(gamma_L + 1.d0)*H_L
    astarR2 = 2d0*(gamma_R - 1d0)/(gamma_R + 1.d0)*H_R
    astarL = sqrt(astarL2)
    astarR = sqrt(astarR2)

    ! --- numerical interface sound speed (strict AUSM form) ---
    ahatL = astarL2 / max(astarL, unl)
    ahatR = astarR2 / max(astarR, -unr)
    aF = min(ahatL, ahatR)

    ! --- Mach numbers and fa-based dynamic ---
    mL = unl / aF
    mR = unr / aF
    MFsq = 0.5d0*(mL*mL + mR*mR)
    param1 = max(MFsq, obj_riemann%Minf*obj_riemann%Minf)
    Mrefsq = min(1.0d0, param1)
    fa = 2.d0*sqrt(Mrefsq) - Mrefsq
    alpha = 3.d0/16.d0*(-4.d0 + 5.d0*fa*fa)

    ! --- Mach/pressure splitting functions ---
    if (abs(mL) <= 1.d0) then
      p1 = 0.25d0*(mL+1.d0)**2
      p2 = (mL*mL-1.d0)**2
      mLP = p1 + beta*p2
      pLP = p1*(2.d0 - mL) + alpha*mL*p2
    else
      mLP = 0.5d0*(mL + abs(mL))
      pLP = mLP / mL
    endif

    if (abs(mR) <= 1.d0) then
      p1 = 0.25d0*(mR-1.d0)**2
      p2 = (mR*mR-1.d0)**2
      mRM = -p1 - beta*p2
      pRM =  p1*(2.d0 + mR) - alpha*mR*p2
    else
      mRM = 0.5d0*(mR - abs(mR))
      pRM = mRM / mR
    endif

    ! --- mass flux with pressure diffusion ---
    rhoF = 0.5d0*(dltot + drtot)
    Mp = -(Kp/fa)*max(0.d0, 1.d0 - sigma*MFsq)*(pr - pl)/(rhoF*aF*aF)
    mF = mLP + mRM + Mp
    mdot = aF*( max(mF,0.d0)*dltot + min(mF,0.d0)*drtot )

    ! --- modified pressure flux ---
    phalf = 0.5d0*(pr + pl) + 0.5d0*(pLP - pRM)*(pl - pr) + sqrt(sq_vel)*(pLP + pRM - 1.d0)*rhoF*aF

    ! --- Momentum and energy fluxes: use mdot upwinding + pressure acoustic term ---
    if (mF >= 0.d0) then
      ! upwind left
      upwind_u = ul
      upwind_v = vl
      upwind_w = wl
      upwind_H = H_L
    else
      upwind_u = ur
      upwind_v = vr
      upwind_w = wr
      upwind_H = H_R
    endif

    ! mass flux output
    F_r = mdot

    ! momentum flux: mdot * upwind velocity + pressure * n
    F_u = mdot * upwind_u + nx * phalf
    F_v = mdot * upwind_v + ny * phalf
    F_w = mdot * upwind_w + nz * phalf

    ! energy flux: mdot * upwind enthalpy (total)
    F_E = mdot * upwind_H

  end subroutine riemann_AUSMp_up2

end module MOSE_Lib_Riemann_AUSM