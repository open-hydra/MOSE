! OVERFLOW2.1 Riemann solvers
! Reference:
!   Tramel, R., Nichols, R., & Buning, P. (2009). Addition of improved shock-capturing schemes to OVERFLOW 2.1. In 19th AIAA Computational Fluid Dynamics (p. 3988).
! Note:
! - HLLC+  | The minimum allowed beta is set to 0.0 instead of 0.4 to allow pure HLLE flows
! - HLLE++ | The Roe lambda_1 is set to un_roe instead of max(un_roe,a_roe). 
!            The adopted forumlation preserves pure Eulerian shear flow, contrary to the reference one.

module MOSE_Lib_Riemann_HLL
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_HLLE, riemann_HLLC
  public :: riemann_HLLCp_Chen
  public :: riemann_HLLCprec
  public :: riemann_HLLEpp
  public :: riemann_HLLCSD
  public :: riemann_HLLCHLLE

contains

  subroutine riemann_HLLC(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    implicit none
    ! Inputs
    real(R8), intent(in) :: dl(nsc), dr(nsc)
    real(R8), intent(in) :: ul, vl, wl, ur, vr, wr
    real(R8), intent(in) :: pl, pr, al, ar
    real(R8), intent(in) :: dltot, drtot
    real(R8), intent(in) :: nx, ny, nz
    real(R8), intent(in) :: dummy, url, urr
    ! Outputs
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_E
    ! Roe + flow variables
    real(R8) :: rho_roe, u_roe, v_roe, w_roe, a_roe, un_roe, h0_roe
    real(R8) :: h0l, h0r, unl, unr
    ! HLLC variables
    real(R8) :: S1, S4, Sstar, pstar
    real(R8) :: inv_denom
    real(R8) :: U1S, U2S, U3S, U4S, U5S
    real(R8) :: e0l, e0r

    unl = ul*nx + vl*ny + wl*nz 
    unr = ur*nx + vr*ny + wr*nz

    !------------------------------------------------
    ! Roe averages
    call roe_averages(nsc, dl, dr, dltot, drtot, &
                      ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, &
                      rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)

    un_roe = u_roe*nx + v_roe*ny + w_roe*nz

    ! Total energy (conservative variable) from total enthalpy
    e0l = h0l - pl / dltot
    e0r = h0r - pr / drtot

    !------------------------------------------------
    ! Batten signal speeds

    S1 = min(0.d0, unl - al, un_roe - a_roe)
    S4 = max(0.d0, unr + ar, un_roe + a_roe)

    inv_denom = 1.d0 / ( dltot*(S1-unl) - drtot*(S4-unr) )

    Sstar = ( pr - pl + dltot*unl*(S1-unl) - drtot*unr*(S4-unr) ) * inv_denom
    pstar = dltot*(unl-S1)*(unl-Sstar) + pl

    !------------------------------------------------
    ! HLLC fluxes
    if (S1 >= 0.d0) then

      F_r = dltot * unl
      F_u = F_r*ul + pl*nx
      F_v = F_r*vl + pl*ny
      F_w = F_r*wl + pl*nz
      F_E = F_r*h0l

    elseif (Sstar >= 0.d0) then

      F_r = dltot * unl
      F_u = F_r*ul + pl*nx
      F_v = F_r*vl + pl*ny
      F_w = F_r*wl + pl*nz
      F_E = F_r*h0l

      inv_denom = 1.d0 / (S1 - Sstar)

      U1S = dltot*(S1-unl)*inv_denom
      U2S = ((S1-unl)*dltot*ul + (pstar-pl)*nx)*inv_denom
      U3S = ((S1-unl)*dltot*vl + (pstar-pl)*ny)*inv_denom
      U4S = ((S1-unl)*dltot*wl + (pstar-pl)*nz)*inv_denom
      U5S = ((S1-unl)*dltot*e0l - pl*unl + pstar*Sstar)*inv_denom

      F_r = F_r + S1*(U1S - dltot)
      F_u = F_u + S1*(U2S - dltot*ul)
      F_v = F_v + S1*(U3S - dltot*vl)
      F_w = F_w + S1*(U4S - dltot*wl)
      F_E = F_E + S1*(U5S - dltot*e0l)

    elseif (S4 > 0.d0) then

      F_r = drtot * unr
      F_u = F_r*ur + pr*nx
      F_v = F_r*vr + pr*ny
      F_w = F_r*wr + pr*nz
      F_E = F_r*h0r

      inv_denom = 1.d0 / (S4 - Sstar)

      U1S = drtot*(S4-unr)*inv_denom
      U2S = ((S4-unr)*drtot*ur + (pstar-pr)*nx)*inv_denom
      U3S = ((S4-unr)*drtot*vr + (pstar-pr)*ny)*inv_denom
      U4S = ((S4-unr)*drtot*wr + (pstar-pr)*nz)*inv_denom
      U5S = ((S4-unr)*drtot*e0r - pr*unr + pstar*Sstar)*inv_denom

      F_r = F_r + S4*(U1S - drtot)
      F_u = F_u + S4*(U2S - drtot*ur)
      F_v = F_v + S4*(U3S - drtot*vr)
      F_w = F_w + S4*(U4S - drtot*wr)
      F_E = F_E + S4*(U5S - drtot*e0r)

    else

      F_r = drtot * unr
      F_u = F_r*ur + pr*nx
      F_v = F_r*vr + pr*ny
      F_w = F_r*wr + pr*nz
      F_E = F_r*h0r

    endif

   end subroutine riemann_HLLC


  !> @brief HLLC+ : low-Mach shock-stable HLLC (Chen, Lin, Li & Yan, 2020)
  !>
  !> "HLLC+: Low-Mach shock-stable HLLC-type Riemann solver for all-speed flows",
  !> SIAM J. Sci. Comput. 42(4) (2020) B921-B950.
  !>
  !> Standard HLLC plus two corrections added in the subsonic star region
  !> (S_L < 0 < S_R), both with prefactor phi = phi_L*phi_R/(phi_R-phi_L):
  !>  - Low-Mach normal-velocity anti-dissipation A_p (Eq. 3.7-3.15):
  !>      momentum += phi*(f*-1)*dU*n ,  energy += phi*(f*-1)*dU*S_HLLC
  !>    with f(M) = M*sqrt(4+(1-M^2)^2)/(1+M^2), M = Thornber clamp (Eq. 3.14).
  !>    This rescales the O(a) pressure-difference dissipation p_d down to the
  !>    convective scale, recovering Ma^2 pressure/density scaling at low speed.
  !>  - Carbuncle shear viscosity (Eq. 3.22-3.25):
  !>      momentum += phi*(S_K/(S_K-S_HLLC))*g*dU_transverse
  !>    gated by the multidimensional pressure-ratio sensor g (Eq. 3.23-24,
  !>    supplied as 'switch' by the shock-detector framework, SD_Chen).
  !>
  !> The near-shock restraint of f* (paper's lambda sonic sensor, Eq. 3.16) is
  !> realised through the same g sensor: f* = f(M) + g*(1-f(M)), so f* -> 1 at a
  !> shock (fix off) and f* -> f(M) in smooth flow (fix on), with a single sensor.
  !>
  !> @ingroup Lib_RiemannPrivateProcedure
  subroutine riemann_HLLCp_Chen(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,switch,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    implicit none
    real(R8), intent(in) :: dl(nsc), dr(nsc)
    real(R8), intent(in) :: ul, vl, wl, ur, vr, wr
    real(R8), intent(in) :: pl, pr, al, ar
    real(R8), intent(in) :: dltot, drtot
    real(R8), intent(in) :: nx, ny, nz
    real(R8), intent(in) :: switch, url, urr
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_E
    ! Roe + flow variables
    real(R8) :: rho_roe, u_roe, v_roe, w_roe, a_roe, un_roe, h0_roe
    real(R8) :: h0l, h0r, unl, unr
    ! HLLC variables
    real(R8) :: S1, S4, Sstar, pstar, inv_denom
    real(R8) :: U1S, U2S, U3S, U4S, U5S, e0l, e0r
    ! Chen low-Mach + shear corrections
    real(R8) :: phiL, phiR, pref, Mloc, fM, g, fstar, SK
    real(R8) :: dU, dut_x, dut_y, dut_z, lm_coef, shear_coef

    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz

    call roe_averages(nsc, dl, dr, dltot, drtot, &
                      ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, &
                      rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)
    un_roe = u_roe*nx + v_roe*ny + w_roe*nz

    e0l = h0l - pl/dltot
    e0r = h0r - pr/drtot

    ! Batten signal speeds
    S1 = min(0.d0, unl - al, un_roe - a_roe)
    S4 = max(0.d0, unr + ar, un_roe + a_roe)

    inv_denom = 1.d0 / ( dltot*(S1-unl) - drtot*(S4-unr) )
    Sstar = ( pr - pl + dltot*unl*(S1-unl) - drtot*unr*(S4-unr) ) * inv_denom
    pstar = dltot*(unl-S1)*(unl-Sstar) + pl

    ! --- Chen all-speed corrections (active only in the subsonic star region) ---
    phiL = dltot*(S1 - unl)
    phiR = drtot*(S4 - unr)
    pref = phiL*phiR / (phiR - phiL)                 ! ~ -1/2 rho a (Eq. 3.20)
    Mloc = min(1.d0, max( sqrt(ul*ul+vl*vl+wl*wl)/al, &
                          sqrt(ur*ur+vr*vr+wr*wr)/ar ))      ! Eq. 3.14
    fM   = Mloc*sqrt(4.d0 + (1.d0 - Mloc*Mloc)**2)/(1.d0 + Mloc*Mloc)   ! Eq. 3.13
    g    = 1.d0 - switch**Mloc                        ! g = 1 - h^M (Eq. 3.23), switch = h
    fstar = fM + g*(1.d0 - fM)                        ! f* -> 1 near shock, f(M) in smooth flow
    dU    = unr - unl                                 ! normal velocity jump
    dut_x = (ur-ul) - dU*nx                           ! transverse velocity jump (Eq. 3.25)
    dut_y = (vr-vl) - dU*ny
    dut_z = (wr-wl) - dU*nz

    !------------------------------------------------
    ! HLLC fluxes (+ Chen corrections in the two star branches)
    if (S1 >= 0.d0) then

      F_r = dltot*unl
      F_u = F_r*ul + pl*nx
      F_v = F_r*vl + pl*ny
      F_w = F_r*wl + pl*nz
      F_E = F_r*h0l

    elseif (Sstar >= 0.d0) then

      F_r = dltot*unl
      F_u = F_r*ul + pl*nx
      F_v = F_r*vl + pl*ny
      F_w = F_r*wl + pl*nz
      F_E = F_r*h0l

      inv_denom = 1.d0 / (S1 - Sstar)
      U1S = dltot*(S1-unl)*inv_denom
      U2S = ((S1-unl)*dltot*ul + (pstar-pl)*nx)*inv_denom
      U3S = ((S1-unl)*dltot*vl + (pstar-pl)*ny)*inv_denom
      U4S = ((S1-unl)*dltot*wl + (pstar-pl)*nz)*inv_denom
      U5S = ((S1-unl)*dltot*e0l - pl*unl + pstar*Sstar)*inv_denom

      F_r = F_r + S1*(U1S - dltot)
      F_u = F_u + S1*(U2S - dltot*ul)
      F_v = F_v + S1*(U3S - dltot*vl)
      F_w = F_w + S1*(U4S - dltot*wl)
      F_E = F_E + S1*(U5S - dltot*e0l)

      SK = S1
      lm_coef    = pref*(fstar - 1.d0)*dU
      shear_coef = pref*( SK/(SK - Sstar) )*g
      F_u = F_u + lm_coef*nx + shear_coef*dut_x
      F_v = F_v + lm_coef*ny + shear_coef*dut_y
      F_w = F_w + lm_coef*nz + shear_coef*dut_z
      F_E = F_E + lm_coef*Sstar

    elseif (S4 > 0.d0) then

      F_r = drtot*unr
      F_u = F_r*ur + pr*nx
      F_v = F_r*vr + pr*ny
      F_w = F_r*wr + pr*nz
      F_E = F_r*h0r

      inv_denom = 1.d0 / (S4 - Sstar)
      U1S = drtot*(S4-unr)*inv_denom
      U2S = ((S4-unr)*drtot*ur + (pstar-pr)*nx)*inv_denom
      U3S = ((S4-unr)*drtot*vr + (pstar-pr)*ny)*inv_denom
      U4S = ((S4-unr)*drtot*wr + (pstar-pr)*nz)*inv_denom
      U5S = ((S4-unr)*drtot*e0r - pr*unr + pstar*Sstar)*inv_denom

      F_r = F_r + S4*(U1S - drtot)
      F_u = F_u + S4*(U2S - drtot*ur)
      F_v = F_v + S4*(U3S - drtot*vr)
      F_w = F_w + S4*(U4S - drtot*wr)
      F_E = F_E + S4*(U5S - drtot*e0r)

      SK = S4
      lm_coef    = pref*(fstar - 1.d0)*dU
      shear_coef = pref*( SK/(SK - Sstar) )*g
      F_u = F_u + lm_coef*nx + shear_coef*dut_x
      F_v = F_v + lm_coef*ny + shear_coef*dut_y
      F_w = F_w + lm_coef*nz + shear_coef*dut_z
      F_E = F_E + lm_coef*Sstar

    else

      F_r = drtot*unr
      F_u = F_r*ur + pr*nx
      F_v = F_r*vr + pr*ny
      F_w = F_r*wr + pr*nz
      F_E = F_r*h0r

    endif

  end subroutine riemann_HLLCp_Chen

  subroutine riemann_HLLE(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    implicit none
    ! Inputs
    real(R8), intent(in) :: dl(nsc), dr(nsc)
    real(R8), intent(in) :: ul, vl, wl, ur, vr, wr
    real(R8), intent(in) :: pl, pr, al, ar
    real(R8), intent(in) :: dltot, drtot
    real(R8), intent(in) :: nx, ny, nz
    real(R8), intent(in) :: dummy, url, urr
    ! Outputs
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_E
    ! Roe data
    real(R8) :: rho_roe, u_roe, v_roe, w_roe, a_roe, un_roe, h0_roe
    real(R8) :: h0l, h0r, unl, unr
    ! HLLE variables
    real(R8) :: S1, S4, inv_denom
    real(R8) :: FrL, FrR

    unl = ul*nx + vl*ny + wl*nz 
    unr = ur*nx + vr*ny + wr*nz

    !------------------------------------------------
    ! Roe averages
    call roe_averages(nsc, dl, dr, dltot, drtot, &
                      ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, &
                      rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)

    un_roe = u_roe*nx + v_roe*ny + w_roe*nz

    !------------------------------------------------
    ! Batten signal speeds

    S1 = min(0.d0, unl - al, un_roe - a_roe)
    S4 = max(0.d0, unr + ar, un_roe + a_roe)

    !------------------------------------------------
    ! HLL fluxes

    inv_denom = 1.d0 / (S4 - S1)

    FrL = dltot*unl
    FrR = drtot*unr

    F_r = ( S4*FrL - S1*FrR + S1*S4*(drtot-dltot) ) * inv_denom

    F_u = ( S4*(FrL*ul + pl*nx) - S1*(FrR*ur + pr*nx) &
          + S1*S4*(drtot*ur - dltot*ul) ) * inv_denom

    F_v = ( S4*(FrL*vl + pl*ny) - S1*(FrR*vr + pr*ny) &
          + S1*S4*(drtot*vr - dltot*vl) ) * inv_denom

    F_w = ( S4*(FrL*wl + pl*nz) - S1*(FrR*wr + pr*nz) &
          + S1*S4*(drtot*wr - dltot*wl) ) * inv_denom

    F_E = ( S4*(FrL*h0l) - S1*(FrR*h0r) &
          + S1*S4*( drtot*(h0r - pr/drtot) &
                    -dltot*(h0l - pl/dltot) ) ) * inv_denom

  end subroutine riemann_HLLE


  !> @brief HLLC Batten preconditioned (Luo–Baum–Löhner style)
  subroutine Riemann_HLLCprec(dl, uL, vL, wL, pL, aL, dltot, dr, uR, vR, wR, pR, aR, drtot, beta, url, urr, nx, ny, nz, F_r, F_u, F_v, F_w, F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H
    implicit none
    ! Inputs
    real(R8), intent(in)  :: dl(nsc), dr(nsc)
    real(R8), intent(in)  :: ul, vl, wl, ur, vr, wr
    real(R8), intent(in)  :: pl, pr, al, ar
    real(R8), intent(in)  :: dltot, drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta, url, urr
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_E
    ! Local
    real(R8) :: unL, unR, eL, eR, h0L, h0R, EtotL, EtotR
    real(R8) :: FrL, FuL, FvL, FwL, FeL, FrR, FuR, FvR, FwR, FeR, hL, hR
    real(R8) :: sqrtRhoL, sqrtRhoR, denom
    real(R8) :: u_Roe, v_Roe, w_Roe, un_Roe, a_Roe
    real(R8) :: S1, S4, Sstar, pstar
    real(R8) :: U1S, U2S, U3S, U4S, U5S
    real(R8) :: VrL, VrR, VrROE
    real(R8) :: alphaL, alphaR, alphaROE
    real(R8) :: vL_star, vR_star, vROE_star
    real(R8) :: cL_star, cR_star, cROE_star

    ! Normal velocities
    unL = uL*nx + vL*ny + wL*nz
    unR = uR*nx + vR*ny + wR*nz

    hL = H(pl,dl)
    hR = H(pr,dr)

    ! Internal energy: e = h - p/rho
    eL = hL - pL/dltot
    eR = hR - pR/drtot

    ! Total enthalpy: h0 = h + 0.5*v²
    h0L = hL + 0.5d0*(uL**2 + vL**2 + wL**2)
    h0R = hR + 0.5d0*(uR**2 + vR**2 + wR**2)

    ! Total energy: E = e + 0.5*v²
    EtotL = eL + 0.5d0*(uL**2 + vL**2 + wL**2)
    EtotR = eR + 0.5d0*(uR**2 + vR**2 + wR**2)

    ! Left/right fluxes
    FrL = dltot * unL
    FuL = FrL * uL + pL * nx
    FvL = FrL * vL + pL * ny
    FwL = FrL * wL + pL * nz
    FeL = FrL * h0L

    FrR = drtot * unR
    FuR = FrR * uR + pR * nx
    FvR = FrR * vR + pR * ny
    FwR = FrR * wR + pR * nz
    FeR = FrR * h0R

    ! Roe averages
    sqrtRhoL = sqrt(dltot)
    sqrtRhoR = sqrt(drtot)
    denom = sqrtRhoL + sqrtRhoR
    u_Roe = (sqrtRhoL*uL + sqrtRhoR*uR) / denom
    v_Roe = (sqrtRhoL*vL + sqrtRhoR*vR) / denom
    w_Roe = (sqrtRhoL*wL + sqrtRhoR*wR) / denom
    a_Roe = (sqrtRhoL*aL + sqrtRhoR*aR) / denom

    ! Reference velocities from reconstructed cell values (includes viscous floor)
    VrL   = min(UrL, aL)
    VrR   = min(UrR, aR)
    VrROE = (sqrtRhoL*VrL + sqrtRhoR*VrR) / denom

    ! Preconditioning parameters
    alphaL   = (1.0d0 - (VrL  **2) / (aL  **2)) / 2.0d0
    alphaR   = (1.0d0 - (VrR  **2) / (aR  **2)) / 2.0d0
    alphaROE = (1.0d0 - (VrROE**2) / (a_Roe**2)) / 2.0d0

    vL_star   = unL   * (1.0d0 - alphaL)
    vR_star   = unR   * (1.0d0 - alphaR)
    vROE_star = un_Roe * (1.0d0 - alphaROE)

    cL_star   = sqrt(max(alphaL  **2 * unL  **2 + VrL  **2, 0.0d0))
    cR_star   = sqrt(max(alphaR  **2 * unR  **2 + VrR  **2, 0.0d0))
    cROE_star = sqrt(max(alphaROE**2 * un_Roe**2 + VrROE**2, 0.0d0))

    ! Preconditioned signal velocities
    S1 = min(0.0d0, vL_star - cL_star, vROE_star - cROE_star)
    S4 = max(0.0d0, vR_star + cR_star, vROE_star + cROE_star)

    ! Star region estimates (Batten)
    Sstar = ( pR - pL + dltot*unL*(S1 - unL) - drtot*unR*(S4 - unR) ) / ( dltot*(S1 - unL) - drtot*(S4 - unR) )
    pstar = dltot*(unL - S1)*(unL - Sstar) + pL

    if (S1 >= 0.0d0) then
      F_r = FrL; F_u = FuL; F_v = FvL; F_w = FwL; F_E = FeL
    else if (S4 <= 0.0d0) then
      F_r = FrR; F_u = FuR; F_v = FvR; F_w = FwR; F_E = FeR
    else if (Sstar >= 0.0d0) then
      U1S = dltot*(S1 - unL)/(S1 - Sstar)
      U2S = ((S1 - unL)*dltot*uL + (pstar - pL)*nx)/(S1 - Sstar)
      U3S = ((S1 - unL)*dltot*vL + (pstar - pL)*ny)/(S1 - Sstar)
      U4S = ((S1 - unL)*dltot*wL + (pstar - pL)*nz)/(S1 - Sstar)
      U5S = ((S1 - unL)*dltot*EtotL - pL*unL + pstar*Sstar)/(S1 - Sstar)
      F_r = FrL + S1*(U1S - dltot)
      F_u = FuL + S1*(U2S - dltot*uL)
      F_v = FvL + S1*(U3S - dltot*vL)
      F_w = FwL + S1*(U4S - dltot*wL)
      F_E = FeL + S1*(U5S - dltot*EtotL)
    else
      U1S = drtot*(S4 - unR)/(S4 - Sstar)
      U2S = ((S4 - unR)*drtot*uR + (pstar - pR)*nx)/(S4 - Sstar)
      U3S = ((S4 - unR)*drtot*vR + (pstar - pR)*ny)/(S4 - Sstar)
      U4S = ((S4 - unR)*drtot*wR + (pstar - pR)*nz)/(S4 - Sstar)
      U5S = ((S4 - unR)*drtot*EtotR - pR*unR + pstar*Sstar)/(S4 - Sstar)
      F_r = FrR + S4*(U1S - drtot)
      F_u = FuR + S4*(U2S - drtot*uR)
      F_v = FvR + S4*(U3S - drtot*vR)
      F_w = FwR + S4*(U4S - drtot*wR)
      F_E = FeR + S4*(U5S - drtot*EtotR)
    end if

  end subroutine Riemann_HLLCprec

  ! HLLE++
  subroutine riemann_HLLEpp(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta, url, urr
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    !----------------------------------------------------
    ! Local variables
    real(R8) :: Sm,Sp,unl,unr
    real(R8) :: Frl,Ful,Fvl,Fwl,Fel
    real(R8) :: Frr,Fur,Fvr,Fwr,Fer
    real(R8) :: h0r,h0l
    real(R8) :: inv_denom
    ! Roe averages
    real(R8) :: rho_ROE,u_ROE,v_ROE,w_ROE
    real(R8) :: a_ROE,un_ROE,h0_roe
    ! Jumps
    real(R8) :: drho,dp,du,dv,dw,dun
    ! Tangential basis
    real(R8) :: t1x,t1y,t1z,t2x,t2y,t2z,normt
    ! Wave strengths
    real(R8) :: alpha_c,alpha_s1,alpha_s2,alpha_m,alpha_p
    ! Eigenvalues
    real(R8) :: lam1_roe, lam2_roe, lam3_roe, eps
    real(R8) :: lam1_hlle, lam2_hlle, lam3_hlle
    real(R8) :: lam1_pp, lam2_pp, lam3_pp

    unl = ul*nx+vl*ny+wl*nz
    unr = ur*nx+vr*ny+wr*nz

    !------------------------------------------------
    ! Roe averages
    call roe_averages(nsc, dl, dr, dltot, drtot, &
                      ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, &
                      rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)

    un_roe = u_roe*nx + v_roe*ny + w_roe*nz

    !----------------------------------------------------
    ! Eigenvalues

    ! HLLE
    Sm = min(0.d0, unl-al, un_ROE-a_ROE)
    Sp = max(0.d0, unr+ar, un_ROE+a_ROE)
    inv_denom = 1d0 / (Sp - Sm)
    lam1_hlle = abs( (Sp+Sm)*un_roe*inv_denom - 2d0*Sp*Sm*inv_denom )
    lam2_hlle = abs( (Sp+Sm)*(un_roe-a_roe)*inv_denom - 2d0*Sp*Sm*inv_denom )
    lam3_hlle = abs( (Sp+Sm)*(un_roe+a_roe)*inv_denom - 2d0*Sp*Sm*inv_denom )

    ! Roe
    eps = 0.3d0*(abs(un_ROE)+a_roe)
    lam1_roe = abs(un_ROE) !max(abs(un_ROE),a_ROE)
    lam2_roe = hartenHyman(un_roe - a_roe, eps)
    lam3_roe = hartenHyman(un_roe + a_roe, eps)

    ! HLLE++
    lam1_pp = beta*lam1_roe + (1.d0 - beta)*lam1_hlle
    lam2_pp = beta*lam2_roe + (1.d0 - beta)*lam2_hlle
    lam3_pp = beta*lam3_roe + (1.d0 - beta)*lam3_hlle

    !----------------------------------------------------
    ! Physical fluxes
    call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,Frl,Ful,Fvl,Fwl,Fel)
    call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,Frr,Fur,Fvr,Fwr,Fer)

    !----------------------------------------------------
    ! Avg fluxes
    F_r = 0.5d0*(Frl + FrR)
    F_u = 0.5d0*(Ful + FuR)
    F_v = 0.5d0*(Fvl + FvR)
    F_w = 0.5d0*(Fwl + FwR)
    F_E = 0.5d0*(Fel + FeR)
    !----------------------------------------------------
    ! Diffusion correction

    ! Jumps
    drho = drtot - dltot
    dp   = pr - pl
    du   = ur - ul
    dv   = vr - vl
    dw   = wr - wl
    dun  = unr - unl

    ! Tangential basis
    if (abs(nx) < 0.9d0) then
      t1x = 0.d0
      t1y = -nz
      t1z =  ny
    else
      t1x = -nz
      t1y = 0.d0
      t1z =  nx
    endif
    normt = dsqrt(t1x*t1x+t1y*t1y+t1z*t1z)
    t1x = t1x/normt ; t1y = t1y/normt ; t1z = t1z/normt

    t2x = ny*t1z - nz*t1y
    t2y = nz*t1x - nx*t1z
    t2z = nx*t1y - ny*t1x

    ! Wave strengths
    alpha_m = (dp - rho_ROE*a_ROE*dun)/(2d0*a_ROE*a_ROE)
    alpha_c  = drho - dp/(a_ROE*a_ROE)
    alpha_s1 = rho_ROE*(du*t1x + dv*t1y + dw*t1z)
    alpha_s2 = rho_ROE*(du*t2x + dv*t2y + dw*t2z)
    alpha_p = (dp + rho_ROE*a_ROE*dun)/(2d0*a_ROE*a_ROE)

    ! Mass
    F_r = F_r - 0.5d0 * ( alpha_c*lam1_pp + alpha_m*lam2_pp + alpha_p*lam3_pp )

    ! Momentum
    ! Contact contribution
    F_u = F_u - 0.5d0*( alpha_c * lam1_pp * u_ROE )
    F_v = F_v - 0.5d0*( alpha_c * lam1_pp * v_ROE )
    F_w = F_w - 0.5d0*( alpha_c * lam1_pp * w_ROE )
    ! Shear contributions
    F_u = F_u - 0.5d0 * lam1_pp * ( alpha_s1*t1x + alpha_s2*t2x )
    F_v = F_v - 0.5d0 * lam1_pp * ( alpha_s1*t1y + alpha_s2*t2y )
    F_w = F_w - 0.5d0 * lam1_pp * ( alpha_s1*t1z + alpha_s2*t2z )
    ! Acoustic contributions
    F_u = F_u - 0.5d0 * ( alpha_m * lam2_pp * (u_roe - a_ROE * nx) + alpha_p * lam3_pp * (u_roe + a_ROE * nx) )
    F_v = F_v - 0.5d0 * ( alpha_m * lam2_pp * (v_roe - a_ROE * ny) + alpha_p * lam3_pp * (v_roe + a_ROE * ny) )
    F_w = F_w - 0.5d0 * ( alpha_m * lam2_pp * (w_roe - a_ROE * nz) + alpha_p * lam3_pp * (w_roe + a_ROE * nz) )

    ! Energy
    ! Contact contribution
    F_E = F_E - 0.5d0 * lam1_pp * ( alpha_c*0.5*(u_ROE*u_ROE + v_ROE*v_ROE + w_ROE*w_ROE) )
    ! Shear contributions
    F_E = F_E - 0.5d0 * lam1_pp * ( alpha_s1*(u_ROE*t1x+v_ROE*t1y+w_ROE*t1z) + alpha_s2*(u_ROE*t2x+v_ROE*t2y+w_ROE*t2z) )
    ! Acoustic contributions
    F_E = F_E - 0.5d0 * ( alpha_m * lam2_pp * (h0_roe - a_ROE * un_ROE) + alpha_p * lam3_pp * (h0_roe + a_ROE * un_ROE) )

  end subroutine riemann_HLLEpp


  ! HLLC - HLLE 
  ! Hybrid solver working with a mulidimensional shock detector
  subroutine riemann_HLLCSD(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta, url, urr
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! 
    real(R8) :: FrE, FuE, FvE, FwE, FeE
    real(R8) :: FrC, FuC, FvC, FwC, FeC
    real(R8) :: beta_

    beta_ = beta
    ! beta_ = max(beta_,min_beta)

    if (beta_ == 0d0) then
      call riemann_HLLE (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    elseif (beta_ >= 1d0) then
      call riemann_HLLC (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    else
      call riemann_HLLE (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,url,urr,nx,ny,nz,FrE,FuE,FvE,FwE,FeE)
      call riemann_HLLC (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,url,urr,nx,ny,nz,FrC,FuC,FvC,FwC,FeC)
      F_r = beta_*FrC + (1d0-beta_)*FrE
      F_u = beta_*FuC + (1d0-beta_)*FuE
      F_v = beta_*FvC + (1d0-beta_)*FvE
      F_w = beta_*FwC + (1d0-beta_)*FwE
      F_e = beta_*FeC + (1d0-beta_)*FeE
    endif

  end subroutine riemann_HLLCSD


  pure function hartenHyman(v,e) result(l)
    implicit none
    real(R8), intent(in) :: v, e
    real(R8) :: l
    real(R8) :: al

    al = abs(v)
    l = al
    if (al<=e) l = (v*v + e*e)/(2d0*e)

  end function hartenHyman


  pure subroutine roe_averages(nsc, dl, dr, dltot, drtot, ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in)   :: nsc
    real(R8), intent(in)  :: dl(nsc), dr(nsc)
    real(R8), intent(in)  :: dltot, drtot
    real(R8), intent(in)  :: ul, vl, wl, ur, vr, wr
    real(R8), intent(in)  :: pl, pr
    real(R8), intent(out) :: h0l, h0r
    real(R8), intent(out) :: rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe
    ! Local 
    integer  :: s, Til, Tir
    real(R8) :: Rl, Rr, Tl, Tr, dTl, dTr
    real(R8) :: inv_dltot, inv_drtot
    real(R8) :: srL, srR, inv_sr
    real(R8) :: cv_roe, sum_ei
    real(R8) :: R_roe, T_roe
    real(R8) :: gam_roe, vel2
    real(R8) :: hl, hr, el, er, invW, d_roe

    !------------------------------------------------
    ! Precompute thermodynamics
    Rl = f_Rtot(dl)
    Rr = f_Rtot(dr)

    inv_dltot = 1.d0 / dltot
    inv_drtot = 1.d0 / drtot

    Tl = pl * inv_dltot / Rl
    Tr = pr * inv_drtot / Rr

    Til = int(Tl)
    Tir = int(Tr)
    dTl = Tl - Til
    dTr = Tr - Tir

    srL = sqrt(dltot)
    srR = sqrt(drtot)
    inv_sr = 1.d0 / (srL + srR)

    !------------------------------------------------
    ! Roe flow averages
    rho_roe = srL * srR

    u_roe = (srR*ur + srL*ul) * inv_sr
    v_roe = (srR*vr + srL*vl) * inv_sr
    w_roe = (srR*wr + srL*wl) * inv_sr

    R_roe = (srR*Rr + srL*Rl) * inv_sr
    T_roe = (srR*Tr + srL*Tl) * inv_sr

    !------------------------------------------------
    ! Species loop
    h0l = 0.d0
    h0r = 0.d0
    cv_roe = 0.d0
    sum_ei = 0.d0

    do s = 1, nsc

      invW = Runiv / Wm_tab(s)

      hl = h_tab(Til, s) + (h_tab(Til+1, s) - h_tab(Til, s)) * dTl
      hr = h_tab(Tir, s) + (h_tab(Tir+1, s) - h_tab(Tir, s)) * dTr

      el = hl - invW * Tl
      er = hr - invW * Tr

      h0l = h0l + hl * dl(s) * inv_dltot
      h0r = h0r + hr * dr(s) * inv_drtot

      d_roe = (srR*dr(s)*inv_drtot + srL*dl(s)*inv_dltot) * inv_sr

      cv_roe = cv_roe + d_roe * ( &
          0.5d0 * ( &
            cp_tab(Til, s) + (cp_tab(Til+1, s)-cp_tab(Til, s))*dTl + &
            cp_tab(Tir, s) + (cp_tab(Tir+1, s)-cp_tab(Tir, s))*dTr ) &
          - invW )

      sum_ei = sum_ei + d_roe * (srR*er + srL*el) * inv_sr

    end do

    !------------------------------------------------
    ! Final scalars
    h0l = h0l + 0.5d0*(ul*ul + vl*vl + wl*wl)
    h0r = h0r + 0.5d0*(ur*ur + vr*vr + wr*wr)

    h0_roe = (srR*h0r + srL*h0l) * inv_sr

    gam_roe = 1.d0 + R_roe / cv_roe

    vel2 = u_roe*u_roe + v_roe*v_roe + w_roe*w_roe

    a_roe = sqrt( (gam_roe - 1.d0) * &
          ( h0_roe - 0.5d0*vel2 + cv_roe*T_roe - sum_ei ) )

  end subroutine roe_averages


  ! Subroutine for computing the conservative fluxes from primitive variables.
  pure subroutine fluxes(p,r,u,v,w,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0
    implicit none
    real(R8), intent(in)  :: p           ! Pressure
    real(R8), intent(in)  :: r(nsc)      ! Density
    real(R8), intent(in)  :: u,v,w       ! Velocity
    real(R8), intent(in)  :: nx,ny,nz    ! Normals
    real(R8), intent(out) :: F_r         ! Flux of mass conservation
    real(R8), intent(out) :: F_u,F_v,F_w ! Flux of momentum conservation
    real(R8), intent(out) :: F_E         ! Flux of energy conservation

    F_r = sum(r)*(u*nx+v*ny+w*nz)
    F_u = F_r*u + p*nx
    F_v = F_r*v + p*ny
    F_w = F_r*w + p*nz
    F_E = F_r*H0(p,r,sqrt(u*u+v*v+w*w))

  end subroutine fluxes


  subroutine riemann_HLLCHLLE  (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc

    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al ! : density(s), velocity, pressure and sound velocity of left state
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar ! : density(s), velocity, pressure and sound velocity of right state
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta, url, urr
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! specific
    real(R8) :: nx1, ny1, nz1, nx2, ny2, nz2
    real(R8) :: alfa1, alfa2
    real(R8) :: F_rHLLE, F_uHLLE, F_vHLLE,F_wHLLE, F_eHLLE,F_rHLLC, F_uHLLC, F_vHLLC, F_wHLLC, F_eHLLC
    real(R8) :: abs_dv, dum

    ! calcolo del modulo del vettore differenza di velocità
    abs_dv=sqrt((ur-ul)**2+(vr-vl)**2+(wr-wl)**2)
    dum=0.5d0*(sqrt(ul**2+vl**2+wl**2)+sqrt(ur**2+vr**2+wr**2))

    if (abs_dv/(dum+1.0d-10)>1.0d-12) then        

      ! calcolo di n1: versore della differenza di velocità
      nx1=(ur-ul)/abs_dv
      ny1=(vr-vl)/abs_dv
      nz1=(wr-wl)/abs_dv

      ! calcolo di alfa1: proiezione di n su n1
      alfa1=nx*nx1+ny*ny1+nz*nz1

      ! rendo alfa1 sempre positivo checkando verso di n1
      nx1=sign(1.d0,alfa1)*nx1
      ny1=sign(1.d0,alfa1)*ny1
      nz1=sign(1.d0,alfa1)*nz1
      alfa1=sign(1.d0,alfa1)*alfa1

      ! dum=sqrt((-nz1-ny1)**2+(-nz1+nx1)**2+(ny1+nx1)**2)
      ! nx2=(-nz1-ny1)/dum
      ! ny2=(-nz1+nx1)/dum
      ! nz2=(nx1+ny1)/dum
      nx2=-ny1
      ny2=nx1
      nz2=nz1
      ! calcolo di alfa2: proiezione di n su n2
      alfa2=nx*nx2+ny*ny2+nz*nz2

      ! rendo alfa2 sempre positivo checkando verso di n2
      nx2=sign(1.d0,alfa2)*nx2
      ny2=sign(1.d0,alfa2)*ny2
      nz2=sign(1.d0,alfa2)*nz2
      alfa2=sign(1.d0,alfa2)*alfa2

      if (nz>0.1d0) then
          alfa1=0.d0
          nx2=nx
          ny2=ny
          nz2=nz
          alfa2=1.d0
      endif
          
      ! chiamata a HLLE in direzione n1
      call riemann_HLLE(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx1,ny1,nz1,F_rHLLE,F_uHLLE,F_vHLLE,F_wHLLE, F_EHLLE)

      ! chiamata a HLLC in direzione n2
      call riemann_HLLC(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx2,ny2,nz2,F_rHLLC,F_uHLLC,F_vHLLC,F_wHLLC, F_EHLLC)

      ! media pesata dei flussi per ottenere flusso rotato
      F_r=F_rHLLE*alfa1+F_rHLLC*alfa2
      F_u=F_uHLLE*alfa1+F_uHLLC*alfa2
      F_v=F_vHLLE*alfa1+F_vHLLC*alfa2
      F_w=F_wHLLE*alfa1+F_wHLLC*alfa2
      if (F_w /= 0.d0 .and. abs(nz) <= 1.d-10) F_w = 0.d0
      F_E=F_EHLLE*alfa1+F_EHLLC*alfa2
    else
      ! chiamata a HLLC in direzione n2, assumendo n1 tangente alla faccia
      call riemann_HLLC(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx,ny,nz,F_rHLLC,F_uHLLC,F_vHLLC,F_wHLLC, F_EHLLC)
      F_r=F_rHLLC
      F_u=F_uHLLC
      F_v=F_vHLLC
      F_w=F_wHLLC
      F_E=F_EHLLC
    endif

  end subroutine riemann_HLLCHLLE

end module MOSE_Lib_Riemann_HLL