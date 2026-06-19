module MOSE_Lib_Riemann_Roe
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_LMRoe

contains

  !> @brief Roe Riemann solver with Rieper low-Mach fix (LMRoe)
  subroutine Riemann_LMRoe(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc), ul, vl, wl, pl, dltot, al
    real(R8), intent(in)  :: dr(nsc), ur, vr, wr, pr, drtot, ar
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: UrL, UrR, beta
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_E

    ! --- locals
    real(R8) :: n1, n2, n3, nmag
    real(R8) :: t1x,t1y,t1z, t2x,t2y,t2z, tmpx,tmpy,tmpz, tmag
    real(R8) :: unL, unR, ut1L, ut1R, ut2L, ut2R
    real(R8) :: EL, ER, H0L, H0R, hL, hR
    real(R8) :: FL1,FL2,FL3,FL4,FL5, FR1,FR2,FR3,FR4,FR5

    real(R8) :: uRoe, vRoe, wRoe, unRoe, utmagRoe
    real(R8) :: HRoe, aRoe, rhoRoe
    real(R8) :: dp, drho, dun, dut1, dut2
    real(R8) :: Ma_loc, phi, dun_mod

    real(R8) :: a2, inva, inva2
    real(R8) :: lam1, lam2, lam3, lam4, lam5
    real(R8) :: abs1, abs2, abs3, abs4, abs5

    real(R8) :: alpha1, alpha2, alpha3, alpha4, alpha5
    real(R8) :: r1_1,r1_2,r1_3,r1_4,r1_5
    real(R8) :: r2_1,r2_2,r2_3,r2_4,r2_5
    real(R8) :: r3_1,r3_2,r3_3,r3_4,r3_5
    real(R8) :: r4_1,r4_2,r4_3,r4_4,r4_5
    real(R8) :: r5_1,r5_2,r5_3,r5_4,r5_5

    real(R8) :: diss1,diss2,diss3,diss4,diss5
    real(R8) :: u2Roe

    ! -------------------------
    ! Normalize face normal
    nmag = sqrt(nx*nx + ny*ny + nz*nz)
    if (nmag <= 0.d0) then
      F_r = 0.d0; F_u = 0.d0; F_v = 0.d0; F_w = 0.d0; F_E = 0.d0
      return
    end if
    n1 = nx/nmag; n2 = ny/nmag; n3 = nz/nmag

    ! Build an orthonormal basis (t1,t2) tangent to n
    ! Choose tmp not parallel to n
    if (abs(n1) < 0.9d0) then
      tmpx = 1.d0; tmpy = 0.d0; tmpz = 0.d0
    else
      tmpx = 0.d0; tmpy = 1.d0; tmpz = 0.d0
    end if
    ! t1 = normalize(tmp - (tmp·n) n)
    t1x = tmpx - (tmpx*n1 + tmpy*n2 + tmpz*n3)*n1
    t1y = tmpy - (tmpx*n1 + tmpy*n2 + tmpz*n3)*n2
    t1z = tmpz - (tmpx*n1 + tmpy*n2 + tmpz*n3)*n3
    tmag = sqrt(t1x*t1x + t1y*t1y + t1z*t1z)
    if (tmag <= 0.d0) then
      ! fallback
      t1x = -n2; t1y =  n1; t1z = 0.d0
      tmag = sqrt(t1x*t1x + t1y*t1y + t1z*t1z)
    end if
    t1x = t1x/tmag; t1y = t1y/tmag; t1z = t1z/tmag
    ! t2 = n x t1
    t2x = n2*t1z - n3*t1y
    t2y = n3*t1x - n1*t1z
    t2z = n1*t1y - n2*t1x

    ! Normal & tangential components
    unL  = uL*n1 + vL*n2 + wL*n3
    unR  = uR*n1 + vR*n2 + wR*n3
    ut1L = uL*t1x + vL*t1y + wL*t1z
    ut1R = uR*t1x + vR*t1y + wR*t1z
    ut2L = uL*t2x + vL*t2y + wL*t2z
    ut2R = uR*t2x + vR*t2y + wR*t2z

    hL = H(pl,dl)
    hR = H(pr,dr)

    ! Total enthalpy H0 = h + 0.5 |u|^2
    H0L = hL + 0.5d0*(uL*uL + vL*vL + wL*wL)
    H0R = hR + 0.5d0*(uR*uR + vR*vR + wR*wR)

    ! Total energy E = h - p/rho + 0.5|u|^2
    EL = H0L - pL/dltot
    ER = H0R - pR/drtot

    ! Physical fluxes projected on n
    FL1 = dltot*unL
    FL2 = dltot*unL*uL + pL*n1
    FL3 = dltot*unL*vL + pL*n2
    FL4 = dltot*unL*wL + pL*n3
    FL5 = dltot*unL*H0L

    FR1 = drtot*unR
    FR2 = drtot*unR*uR + pR*n1
    FR3 = drtot*unR*vR + pR*n2
    FR4 = drtot*unR*wR + pR*n3
    FR5 = drtot*unR*H0R

    call roe_averages(nsc, dl, dr, dltot, drtot, &
                      ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, &
                      rhoroe, uroe, vroe, wroe, aroe, hroe)
    
    unRoe = uRoe*n1 + vRoe*n2 + wRoe*n3

    ! Tangential speed magnitude at Roe state
    utmagRoe = sqrt( max(0.d0, u2Roe - unRoe*unRoe ) )

    ! Low-Mach factor (Rieper): phi = min(1, (|un| + |ut|)/a )
    ! A small floor is retained (as in SU2's LMRoe/L2Roe) to keep a minimum of acoustic upwinding.
    if (aRoe > 0.d0) then
      Ma_loc = (abs(unRoe) + utmagRoe)/aRoe
    else
      Ma_loc = 1.d0
    end if
    phi = max(0.05d0, min(1.d0, Ma_loc))

    ! Jumps
    drho = drtot - dltot
    dp   = pR   - pL
    dun  = unR  - unL
    dut1 = ut1R - ut1L
    dut2 = ut2R - ut2L

    ! Apply LMRoe only to the acoustic part
    dun_mod = phi * dun

    ! Eigenvalues
    lam1 = unRoe - aRoe
    lam2 = unRoe
    lam3 = unRoe
    lam4 = unRoe
    lam5 = unRoe + aRoe

    abs1 = abs(lam1); abs2 = abs(lam2); abs3 = abs(lam3); abs4 = abs(lam4); abs5 = abs(lam5)

    a2    = aRoe*aRoe
    inva  = 1.d0/aRoe
    inva2 = 1.d0/a2

    ! Wave strengths (standard Roe for Euler, but with dun_mod in the acoustic pair)
    ! Entropy wave:
    alpha4 = drho - dp*inva2
    ! Shear waves (two tangential directions):
    alpha2 = rhoRoe * dut1
    alpha3 = rhoRoe * dut2
    ! Acoustic waves:
    alpha1 = 0.5d0*( dp*inva2 - rhoRoe*dun_mod*inva )
    alpha5 = 0.5d0*( dp*inva2 + rhoRoe*dun_mod*inva )

    ! Right eigenvectors in global coordinates (u2Roe computed above)

    ! k=1 (un-a)
    r1_1 = 1.d0
    r1_2 = uRoe - aRoe*n1
    r1_3 = vRoe - aRoe*n2
    r1_4 = wRoe - aRoe*n3
    r1_5 = HRoe - aRoe*unRoe

    ! k=5 (un+a)
    r5_1 = 1.d0
    r5_2 = uRoe + aRoe*n1
    r5_3 = vRoe + aRoe*n2
    r5_4 = wRoe + aRoe*n3
    r5_5 = HRoe + aRoe*unRoe

    ! k=4 entropy (un)
    r4_1 = 1.d0
    r4_2 = uRoe
    r4_3 = vRoe
    r4_4 = wRoe
    r4_5 = 0.5d0*u2Roe

    ! k=2 shear along t1 (un)
    r2_1 = 0.d0
    r2_2 = t1x
    r2_3 = t1y
    r2_4 = t1z
    r2_5 = uRoe*t1x + vRoe*t1y + wRoe*t1z  ! = ut1 at Roe state

    ! k=3 shear along t2 (un)
    r3_1 = 0.d0
    r3_2 = t2x
    r3_3 = t2y
    r3_4 = t2z
    r3_5 = uRoe*t2x + vRoe*t2y + wRoe*t2z  ! = ut2 at Roe state

    ! Dissipation term
    diss1 = abs1*alpha1*r1_1 + abs2*alpha2*r2_1 + abs3*alpha3*r3_1 + abs4*alpha4*r4_1 + abs5*alpha5*r5_1
    diss2 = abs1*alpha1*r1_2 + abs2*alpha2*r2_2 + abs3*alpha3*r3_2 + abs4*alpha4*r4_2 + abs5*alpha5*r5_2
    diss3 = abs1*alpha1*r1_3 + abs2*alpha2*r2_3 + abs3*alpha3*r3_3 + abs4*alpha4*r4_3 + abs5*alpha5*r5_3
    diss4 = abs1*alpha1*r1_4 + abs2*alpha2*r2_4 + abs3*alpha3*r3_4 + abs4*alpha4*r4_4 + abs5*alpha5*r5_4
    diss5 = abs1*alpha1*r1_5 + abs2*alpha2*r2_5 + abs3*alpha3*r3_5 + abs4*alpha4*r4_5 + abs5*alpha5*r5_5

    ! Final flux
    F_r = 0.5d0*(FL1 + FR1) - 0.5d0*diss1
    F_u = 0.5d0*(FL2 + FR2) - 0.5d0*diss2
    F_v = 0.5d0*(FL3 + FR3) - 0.5d0*diss3
    F_w = 0.5d0*(FL4 + FR4) - 0.5d0*diss4
    F_E = 0.5d0*(FL5 + FR5) - 0.5d0*diss5

  end subroutine Riemann_LMRoe



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

end module MOSE_Lib_Riemann_Roe
