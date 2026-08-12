module MOSE_Lib_Riemann_Roe
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_LMRoe
  public :: riemann_MiczekRoe

contains

  !> @brief Roe Riemann solver with Rieper low-Mach fix (LMRoe)
  subroutine Riemann_LMRoe(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use MOSE_Config_Types_m, only: obj_riemann
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
    phi = max(obj_riemann%Mco, min(1.d0, Ma_loc))

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


  !> @brief Miczek preconditioned Roe Riemann solver (low-Mach).
  !!
  !! Implements the low-Mach preconditioned Roe flux of
  !!   Miczek, Roepke & Edelmann, "New numerical solver for flows at various
  !!   Mach numbers", A&A (2015), arXiv:1409.8289v2.
  !!
  !! The numerical flux (their Eq. 35) is
  !!   F = 1/2 (F_L + F_R) - 1/2 P^{-1}|P A|_roe (U_R - U_L),
  !! where only the upwind/dissipation term is preconditioned; the central
  !! (physical) part is unchanged, so the scheme is conservative and reduces
  !! exactly to the standard Roe solver at sonic/supersonic Mach numbers.
  !!
  !! The preconditioning matrix P_V (their Eq. 37) is given in primitive
  !! variables V=(rho,u,v,w,p) with parameter delta = 1/mu - 1, where
  !! mu = min(1, max(M, M_cut)) is the local Mach number limited from above by 1
  !! and from below by M_cut. delta -> 0 (P -> I) for M >= 1.
  !!
  !! Note on non-dimensionalization: the paper is written in non-dimensional
  !! form where a reference Mach number M_r appears in both the flux Jacobian
  !! A_V (the pressure-gradient term scales as 1/M_r^2, eigenvalues u_n +/- c/M_r)
  !! and in P_V. Converting the preconditioned dissipation back to dimensional
  !! variables, every M_r cancels: the dimensional scheme is recovered by setting
  !! M_r = 1 (i.e. dimensional sound speed c, eigenvalues u_n +/- c). Only the
  !! *local* limited Mach number mu (through delta) remains, exactly as it should.
  !!
  !! Implementation strategy (kept parallel to Riemann_LMRoe above):
  !!   * Rotate into a face-aligned orthonormal frame (n,t1,t2). The two
  !!     tangential velocities and any passive scalars are pure advection
  !!     (eigenvalue u_n) and decouple from the preconditioner.
  !!   * The coupled (rho,u_n,p) sub-system reduces to a 3x3 problem. Its
  !!     preconditioned Jacobian PA3 = P3 A3 has eigenvalues u_n, u_n +/- s with
  !!       s = sqrt( (1+delta^2) c^2 - delta^2 u_n^2 )  (>= c for subsonic flow).
  !!   * |PA3| (R|Lambda|R^{-1}) is evaluated through the Lagrange-Sylvester
  !!     matrix-function interpolation, avoiding any explicit (and fragile)
  !!     eigenvector algebra:
  !!       |PA3| = sum_k |lam_k| prod_{j/=k} (PA3 - lam_j I)/(lam_k - lam_j).
  !!   * P3^{-1} has a simple closed form (det P3 = 1 + delta^2).
  !!   * The resulting dissipated primitive jumps are mapped back to conservative
  !!     variables with the standard ideal-gas transform dU = (dU/dV) dV, using
  !!     the Roe-averaged gamma (identical convention to Riemann_LMRoe, so that
  !!     for delta = 0 this routine is bit-for-bit the standard Roe flux).
  subroutine Riemann_MiczekRoe(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use MOSE_Config_Types_m, only: obj_riemann
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc), ul, vl, wl, pl, dltot, al
    real(R8), intent(in)  :: dr(nsc), ur, vr, wr, pr, drtot, ar
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: UrL, UrR, beta
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_E

    !> Lower bound on the local Mach number used to build the preconditioner
    !! (M_cut in the paper). It avoids the singularity of P_V at stagnation
    !! points (where M -> 0 would drive delta = 1/mu - 1 -> infinity) and should
    !! be set on the order of the problem's reference Mach number. It is read
    !! from the `riemann-options-Mco` input key (obj_riemann%Mco); a small
    !! hard floor is kept as a safety net if Mco is left at 0.
    real(R8) :: Mcut

    ! --- locals
    real(R8) :: n1, n2, n3, nmag
    real(R8) :: t1x,t1y,t1z, t2x,t2y,t2z, tmpx,tmpy,tmpz, tmag
    real(R8) :: unL, unR, ut1L, ut1R, ut2L, ut2R
    real(R8) :: EL, ER, H0L, H0R, hL, hR
    real(R8) :: FL1,FL2,FL3,FL4,FL5, FR1,FR2,FR3,FR4,FR5

    real(R8) :: sRL, sRR, denom
    real(R8) :: uRoe, vRoe, wRoe, unRoe, HRoe, aRoe, rhoRoe, u2Roe
    real(R8) :: kappaL, kappaR, kappaRoe, chiL, chiR, chiRoe

    real(R8) :: drho, dp, dun, dut1, dut2
    real(R8) :: cc, Mloc, mu, del, del2, om

    ! 3x3 preconditioned Jacobian PA3 (acting on (rho,u_n,p)), entries
    real(R8) :: p11,p12,p13, p22,p23, p32,p33
    real(R8) :: l_e, l_p, l_m, s, s2
    real(R8) :: va1,va2,va3, vb1,vb2,vb3
    real(R8) :: w1,w2,w3, x1,x2,x3, z1,z2,z3
    real(R8) :: c1,c2,c3, g1,g2,g3
    real(R8) :: dd_rho, dd_un, dd_p, dd_ut1, dd_ut2
    real(R8) :: ddux, dduy, dduz

    ! -------------------------
    ! Normalize face normal
    nmag = sqrt(nx*nx + ny*ny + nz*nz)
    if (nmag <= 0.d0) then
      F_r = 0.d0; F_u = 0.d0; F_v = 0.d0; F_w = 0.d0; F_E = 0.d0
      return
    end if
    n1 = nx/nmag; n2 = ny/nmag; n3 = nz/nmag

    ! Build an orthonormal basis (t1,t2) tangent to n
    if (abs(n1) < 0.9d0) then
      tmpx = 1.d0; tmpy = 0.d0; tmpz = 0.d0
    else
      tmpx = 0.d0; tmpy = 1.d0; tmpz = 0.d0
    end if
    t1x = tmpx - (tmpx*n1 + tmpy*n2 + tmpz*n3)*n1
    t1y = tmpy - (tmpx*n1 + tmpy*n2 + tmpz*n3)*n2
    t1z = tmpz - (tmpx*n1 + tmpy*n2 + tmpz*n3)*n3
    tmag = sqrt(t1x*t1x + t1y*t1y + t1z*t1z)
    if (tmag <= 0.d0) then
      t1x = -n2; t1y =  n1; t1z = 0.d0
      tmag = sqrt(t1x*t1x + t1y*t1y + t1z*t1z)
    end if
    t1x = t1x/tmag; t1y = t1y/tmag; t1z = t1z/tmag
    ! t2 = n x t1
    t2x = n2*t1z - n3*t1y
    t2y = n3*t1x - n1*t1z
    t2z = n1*t1y - n2*t1x

    ! Normal & tangential velocity components
    unL  = uL*n1 + vL*n2 + wL*n3
    unR  = uR*n1 + vR*n2 + wR*n3
    ut1L = uL*t1x + vL*t1y + wL*t1z
    ut1R = uR*t1x + vR*t1y + wR*t1z
    ut2L = uL*t2x + vL*t2y + wL*t2z
    ut2R = uR*t2x + vR*t2y + wR*t2z

    hL = H(pl,dl)
    hR = H(pr,dr)

    ! Total enthalpy and energy
    H0L = hL + 0.5d0*(uL*uL + vL*vL + wL*wL)
    H0R = hR + 0.5d0*(uR*uR + vR*vR + wR*wR)
    EL  = H0L - pL/dltot
    ER  = H0R - pR/drtot

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

    ! Roe averages (density-weighted) -- identical convention to Riemann_LMRoe
    sRL   = sqrt(dltot)
    sRR   = sqrt(drtot)
    denom = sRL + sRR

    uRoe  = (sRL*uL + sRR*uR)/denom
    vRoe  = (sRL*vL + sRR*vR)/denom
    wRoe  = (sRL*wL + sRR*wR)/denom
    HRoe  = (sRL*H0L + sRR*H0R)/denom
    unRoe = uRoe*n1 + vRoe*n2 + wRoe*n3
    rhoRoe = sRL*sRR
    u2Roe  = uRoe*uRoe + vRoe*vRoe + wRoe*wRoe

    ! Roe sound speed (general-gas consistent, see Riemann_LMRoe)
    kappaL   = f_gamma(dl, pL, dltot, f_Rtot(dl)) - 1.d0
    kappaR   = f_gamma(dr, pR, drtot, f_Rtot(dr)) - 1.d0
    chiL     = aL*aL - kappaL*hL
    chiR     = aR*aR - kappaR*hR
    kappaRoe = (sRL*kappaL + sRR*kappaR)/denom
    chiRoe   = (sRL*chiL   + sRR*chiR  )/denom
    aRoe     = sqrt( max(1.d-30, chiRoe + kappaRoe*(HRoe - 0.5d0*u2Roe)) )
    cc       = aRoe

    ! Local limited Mach number and preconditioning parameter delta = 1/mu - 1.
    ! delta -> 0 recovers the standard Roe flux (mu = 1, M >= 1).
    Mcut = max(1.d-6, obj_riemann%Mco)
    Mloc = sqrt(u2Roe)/cc
    mu   = min(1.d0, max(Mloc, Mcut))
    del  = 1.d0/mu - 1.d0
    del2 = del*del
    om   = 1.d0/(1.d0 + del2)            ! 1/det(P3)

    ! Primitive jumps (R - L)
    drho = drtot - dltot
    dp   = pR   - pL
    dun  = unR  - unL
    dut1 = ut1R - ut1L
    dut2 = ut2R - ut2L

    ! Preconditioned Jacobian PA3 = P3*A3 of the coupled (rho,u_n,p) block,
    ! with P3 = [[1,a,0],[0,1,b],[0,e,1]], a=rho*del/c, b=-del/(rho*c), e=rho*c*del,
    ! and A3 = [[u_n,rho,0],[0,u_n,1/rho],[0,rho*c^2,u_n]] (dimensional, M_r=1).
    p11 = unRoe
    p12 = rhoRoe*(1.d0 + del*unRoe/cc)
    p13 = del/cc
    p22 = unRoe - del*cc
    p23 = (1.d0 - del*unRoe/cc)/rhoRoe
    p32 = rhoRoe*cc*(cc + del*unRoe)
    p33 = unRoe + cc*del

    ! Eigenvalues of PA3: u_n, u_n +/- s
    s2  = (1.d0 + del2)*cc*cc - del2*unRoe*unRoe
    s   = sqrt( max(1.d-30, s2) )
    s2  = s*s
    l_e = unRoe
    l_p = unRoe + s
    l_m = unRoe - s

    ! |PA3| * (drho,dun,dp) via Lagrange-Sylvester interpolation of |.|:
    !   |PA3| = |l_e|/((l_e-l_p)(l_e-l_m)) (PA3-l_p)(PA3-l_m)
    !         + |l_p|/((l_p-l_e)(l_p-l_m)) (PA3-l_e)(PA3-l_m)
    !         + |l_m|/((l_m-l_e)(l_m-l_p)) (PA3-l_e)(PA3-l_p)
    ! with denominators -s^2, 2 s^2, 2 s^2 respectively.
    call pa3_shift(l_m, drho,dun,dp,  va1,va2,va3)   ! (PA3 - l_m I) w
    call pa3_shift(l_p, drho,dun,dp,  vb1,vb2,vb3)   ! (PA3 - l_p I) w
    call pa3_shift(l_p, va1,va2,va3,  w1,w2,w3)       ! (PA3-l_p)(PA3-l_m) w
    call pa3_shift(l_e, va1,va2,va3,  x1,x2,x3)       ! (PA3-l_e)(PA3-l_m) w
    call pa3_shift(l_e, vb1,vb2,vb3,  z1,z2,z3)       ! (PA3-l_e)(PA3-l_p) w

    c1 = abs(l_e)/(-s2)
    c2 = abs(l_p)/(2.d0*s2)
    c3 = abs(l_m)/(2.d0*s2)
    g1 = c1*w1 + c2*x1 + c3*z1
    g2 = c1*w2 + c2*x2 + c3*z2
    g3 = c1*w3 + c2*x3 + c3*z3

    ! Apply P3^{-1} (closed form) to obtain the dissipated primitive jumps.
    !   P3^{-1} = [[1, -a*om, a*b*om],[0, om, -b*om],[0, -e*om, om]]
    dd_rho = g1 - (rhoRoe*del/cc)*om*g2 - (del2/(cc*cc))*om*g3
    dd_un  = om*( g2 + (del/(rhoRoe*cc))*g3 )
    dd_p   = om*( -rhoRoe*cc*del*g2 + g3 )

    ! Shear (tangential) waves: pure advection, eigenvalue u_n, unaffected by P.
    dd_ut1 = abs(unRoe)*dut1
    dd_ut2 = abs(unRoe)*dut2

    ! Recompose the dissipated velocity jump in global Cartesian coordinates
    ddux = dd_un*n1 + dd_ut1*t1x + dd_ut2*t2x
    dduy = dd_un*n2 + dd_ut1*t1y + dd_ut2*t2y
    dduz = dd_un*n3 + dd_ut1*t1z + dd_ut2*t2z

    ! Map dissipated primitive jumps -> conservative dissipation, dU = (dU/dV) dV
    ! (ideal-gas energy row with Roe-averaged gamma = kappaRoe + 1; identical to
    !  the conservative eigenvector convention used by Riemann_LMRoe).
    F_r = 0.5d0*(FL1 + FR1) - 0.5d0*( dd_rho )
    F_u = 0.5d0*(FL2 + FR2) - 0.5d0*( uRoe*dd_rho + rhoRoe*ddux )
    F_v = 0.5d0*(FL3 + FR3) - 0.5d0*( vRoe*dd_rho + rhoRoe*dduy )
    F_w = 0.5d0*(FL4 + FR4) - 0.5d0*( wRoe*dd_rho + rhoRoe*dduz )
    F_E = 0.5d0*(FL5 + FR5) - 0.5d0*( 0.5d0*u2Roe*dd_rho &
                                    + rhoRoe*(uRoe*ddux + vRoe*dduy + wRoe*dduz) &
                                    + dd_p/kappaRoe )

  contains

    !> Apply (PA3 - lam I) to a 3-vector x, returning r. PA3 entries are
    !! host-associated; its first column is (u_n,0,0)^T (zero below the diagonal).
    pure subroutine pa3_shift(lam, x1,x2,x3, r1,r2,r3)
      real(R8), intent(in)  :: lam, x1, x2, x3
      real(R8), intent(out) :: r1, r2, r3
      r1 = (p11 - lam)*x1 + p12*x2       + p13*x3
      r2 =                  (p22 - lam)*x2 + p23*x3
      r3 =                  p32*x2       + (p33 - lam)*x3
    end subroutine pa3_shift

  end subroutine Riemann_MiczekRoe

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

    ! Table reads go through the species-contiguous companions: tabT(s,T) holds
    ! exactly tab(T,s), so the values are unchanged, but s is now the fast index
    ! and the loop loads unit-stride instead of gathering with stride Tmax-Tmin.
    do s = 1, nsc

      invW = Ri_tab(s)   ! = Runiv/Wm_tab(s), precomputed once in the FLINT loader

      hl = h_tabT(s, Til) + (h_tabT(s, Til+1) - h_tabT(s, Til)) * dTl
      hr = h_tabT(s, Tir) + (h_tabT(s, Tir+1) - h_tabT(s, Tir)) * dTr

      el = hl - invW * Tl
      er = hr - invW * Tr

      h0l = h0l + hl * dl(s) * inv_dltot
      h0r = h0r + hr * dr(s) * inv_drtot

      d_roe = (srR*dr(s)*inv_drtot + srL*dl(s)*inv_dltot) * inv_sr

      cv_roe = cv_roe + d_roe * ( &
          0.5d0 * ( &
            cp_tabT(s, Til) + (cp_tabT(s, Til+1)-cp_tabT(s, Til))*dTl + &
            cp_tabT(s, Tir) + (cp_tabT(s, Tir+1)-cp_tabT(s, Tir))*dTr ) &
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
