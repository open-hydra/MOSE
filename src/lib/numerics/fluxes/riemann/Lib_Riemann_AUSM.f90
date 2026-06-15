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

module MOSE_Lib_Riemann_AUSM
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_AUSMp
  public :: riemann_AUSMp_M

contains

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
  subroutine riemann_AUSMp(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,switch,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0, f_gamma, f_Rtot
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: switch, url, urr
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! specific
    real(R8), parameter :: beta = 0.125d0, alpha = 0.1875d0
    real(R8) :: unl, unr 
    real(R8) :: ML, MR
    real(R8) :: MplusL, MminusR, Mhalf
    real(R8) :: pplusL, pminusR, phalf
    real(R8) :: H_L, H_R
    real(R8) :: a_half, mdotm, mdotp

    ! normal velocities
    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz

    ! total enthalpies
    H_L = H0(pl, dl, sqrt(ul**2 + vl**2 + wl**2))
    H_R = H0(pr, dr, sqrt(ur**2 + vr**2 + wr**2))

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

  !> @brief AUSM+M Riemann solver (Chen, Cai, Xue, Wang & Yan, 2020)
  !>
  !> "An improved AUSM-family scheme with robustness and accuracy for all Mach
  !> number flows", Applied Mathematical Modelling 77 (2020) 1065-1081.
  !>
  !> Three key ingredients over AUSM+/AUSM+-up:
  !>  1. Pressure-diffusion term Mp (Eq.14) whose denominator does NOT contain
  !>     the Mach number, so it has no Kp/fa stagnation singularity and allows a
  !>     larger time step than AUSM+-up at low speed:
  !>        Mp = -1/2 (1-f) (pR-pL)/(rho_1/2 c_1/2^2) (1-g),  f = (1-cos(pi*M))/2
  !>  2. Pressure flux with a scaling Mach function f_o (Eq.19) for low-Mach
  !>     accuracy, plus a multidimensional velocity-diffusion term p_u (Eq.26)
  !>     gated by a pressure-ratio shock sensor g (Eq.25) for carbuncle control.
  !>  3. AUSMPW+ numerical sound speed (Eq.28-29, Kim et al.) for correct oblique
  !>     shocks and no unphysical expansion shocks.
  !>
  !> The shock sensor g uses the local interface pressure ratio (g=1 at a strong
  !> pressure jump, g=0 in smooth flow); the paper's full g additionally scans
  !> all adjacent interfaces of both cells (future multidimensional extension).
  !>
  !> @ingroup Lib_RiemannPrivateProcedure
  subroutine riemann_AUSMp_M(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,switch,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use MOSE_Config_Types_m, only: obj_riemann
    use FLINT_Lib_Thermodynamic, only: H0, H, f_gamma, f_Rtot
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: switch, url, urr
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! parameters
    real(R8), parameter :: beta = 0.125d0, alpha = 0.1875d0   ! 1/8, 3/16
    real(R8), parameter :: pi = 3.141592653589793d0
    ! locals
    real(R8) :: unl, unr, Vl2, Vr2, H_L, H_R
    real(R8) :: gam_L, gam_R, gam, hnorm, cs, c12
    real(R8) :: ML, MR, MLp, MRm, psiLp, psiRm
    real(R8) :: f, g, fo, rho12, Mp, M12, mdot, ps, gfac, pux, puy, puz

    ! --- normal velocities, total enthalpies ---
    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz
    Vl2 = ul*ul + vl*vl + wl*wl
    Vr2 = ur*ur + vr*vr + wr*wr
    H_L = H0(pl, dl, sqrt(Vl2))
    H_R = H0(pr, dr, sqrt(Vr2))

    ! --- AUSMPW+ numerical sound speed (Eq. 28-29) ---
    gam_L = f_gamma(dl, pl, dltot, f_Rtot(dl))
    gam_R = f_gamma(dr, pr, drtot, f_Rtot(dr))
    gam   = 0.5d0*(gam_L + gam_R)
    hnorm = 0.5d0*( H(pl,dl) + H(pr,dr) )                    ! mean static enthalpy
    cs    = sqrt( 2.d0*(gam - 1.d0)/(gam + 1.d0)*hnorm )
    if (unl + unr >= 0.d0) then
      c12 = cs*cs/max(abs(unl), cs)
    else
      c12 = cs*cs/max(abs(unr), cs)
    end if

    ! --- interface Mach numbers and AUSM+ split functions (Eq. 6) ---
    ML = unl/c12 ; MR = unr/c12
    if (abs(ML) < 1.d0) then
      MLp   = 0.25d0*(ML + 1.d0)**2 + beta*(ML*ML - 1.d0)**2
      psiLp = 0.25d0*(ML + 1.d0)**2*(2.d0 - ML) + alpha*ML*(ML*ML - 1.d0)**2
    else
      MLp   = 0.5d0*(ML + abs(ML))
      psiLp = 0.5d0*(1.d0 + sign(1.d0, ML))
    end if
    if (abs(MR) < 1.d0) then
      MRm   = -0.25d0*(MR - 1.d0)**2 - beta*(MR*MR - 1.d0)**2
      psiRm = 0.25d0*(MR - 1.d0)**2*(2.d0 + MR) - alpha*MR*(MR*MR - 1.d0)**2
    else
      MRm   = 0.5d0*(MR - abs(MR))
      psiRm = 0.5d0*(1.d0 - sign(1.d0, MR))
    end if

    ! --- Mach-number function f (Eq. 15, local) and shock sensor g (Eq. 25) ---
    ! 'switch' carries the multidimensional Chen indicator h (h->1 smooth, h->0
    ! at a shock) from the shock-detector framework; map to AUSM+M's gate g.
    f = 0.5d0*(1.d0 - cos(pi*min(1.d0, max(abs(ML), abs(MR)))))
    g = 0.5d0*(1.d0 + cos(pi*switch))

    ! --- pressure-diffusion mass flux (Eq. 14), no Mach in denominator ---
    rho12 = 0.5d0*(dltot + drtot)
    Mp    = -0.5d0*(1.d0 - f)*(pr - pl)/(rho12*c12*c12)*(1.d0 - g)
    M12   = MLp + MRm + Mp
    mdot  = M12*c12*merge(dltot, drtot, M12 >= 0.d0)

    ! --- pressure flux with f_o scaling (Eq. 19-20) ---
    fo = min(1.d0, max(f, obj_riemann%Minf*obj_riemann%Minf))
    ps = 0.5d0*(pl + pr) + 0.5d0*(psiLp - psiRm)*(pl - pr) &
       + fo*(psiLp + psiRm - 1.d0)*0.5d0*(pl + pr)

    ! --- multidimensional velocity diffusion (Eq. 26), global components ---
    gfac = g*gam*(pl + pr)/(2.d0*c12)*psiLp*psiRm
    pux  = -gfac*(ur - ul)
    puy  = -gfac*(vr - vl)
    puz  = -gfac*(wr - wl)

    ! --- assemble flux (Eq. 5, 27) ---
    F_r = mdot
    if (M12 >= 0.d0) then
      F_u = mdot*ul + ps*nx + pux
      F_v = mdot*vl + ps*ny + puy
      F_w = mdot*wl + ps*nz + puz
      F_E = mdot*H_L
    else
      F_u = mdot*ur + ps*nx + pux
      F_v = mdot*vr + ps*ny + puy
      F_w = mdot*wr + ps*nz + puz
      F_E = mdot*H_R
    end if

  end subroutine riemann_AUSMp_M

end module MOSE_Lib_Riemann_AUSM