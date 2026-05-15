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
  public :: riemann_HLLE, riemann_HLLEM, riemann_HLLC
  public :: riemann_HLLEpp
  public :: riemann_HLLEMSD, riemann_HLLCSD
  public :: riemann_HLLCHLLE

  real(R8) :: min_beta=0.0d0

contains

  ! HLLEM
  ! Einfeldt, B., Munz, C.C., Roe, P.L., and Sjogreen, B., "On Godunov-type methods near low densities", J. Comput. Phys., 92 (1991), pp. 273–295.
  subroutine riemann_HLLEM(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: dummy
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e

    !----------------------------------------------------
    real(R8) :: h0r,h0l
    ! Local variables
    real(R8) :: S1,S4,C,delta
    real(R8) :: Frl,Ful,Fvl,Fwl,Fel
    real(R8) :: Frr,Fur,Fvr,Fwr,Fer
    ! Roe averages
    real(R8) :: rho_ROE,u_ROE,v_ROE,w_ROE
    real(R8) :: a_ROE,un_ROE,h0_roe
    ! Jumps
    real(R8) :: drho,dp,du,dv,dw
    ! Tangential basis
    real(R8) :: t1x,t1y,t1z,t2x,t2y,t2z,normt
    ! Wave strengths
    real(R8) :: alpha_c,alpha_s1,alpha_s2

    !------------------------------------------------
    ! Roe averages
    call roe_averages(nsc, dl, dr, dltot, drtot, &
                      ul, vl, wl, ur, vr, wr, pl, pr, h0l, h0r, &
                      rho_roe, u_roe, v_roe, w_roe, a_roe, h0_roe)

    un_roe = u_roe*nx + v_roe*ny + w_roe*nz

    !----------------------------------------------------
    ! Wave speed estimates
    S1 = min(0.d0, (ul*nx+vl*ny+wl*nz)-al, un_ROE-a_ROE)
    S4 = max(0.d0, (ur*nx+vr*ny+wr*nz)+ar, un_ROE+a_ROE)

    !----------------------------------------------------
    ! Physical fluxes
    call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,Frl,Ful,Fvl,Fwl,Fel)
    call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,Frr,Fur,Fvr,Fwr,Fer)

    !----------------------------------------------------
    ! HLL flux
    F_r = (S4*Frl - S1*Frr + S1*S4*(drtot-dltot))/(S4-S1)
    F_u = (S4*Ful - S1*Fur + S1*S4*(drtot*ur-dltot*ul))/(S4-S1)
    F_v = (S4*Fvl - S1*Fvr + S1*S4*(drtot*vr-dltot*vl))/(S4-S1)
    F_w = (S4*Fwl - S1*Fwr + S1*S4*(drtot*wr-dltot*wl))/(S4-S1)
    F_E = (S4*Fel - S1*Fer + S1*S4*(drtot*E0(pr,dr,dsqrt(ur*ur+vr*vr+wr*wr)) &
          - dltot*E0(pl,dl,dsqrt(ul*ul+vl*vl+wl*wl))))/(S4-S1)

    !----------------------------------------------------
    ! HLLEM correction (only if S1 < 0 < S4)
    if (S1 < 0.d0 .and. S4 > 0.d0) then

      ! Jumps
      drho = drtot - dltot
      dp   = pr - pl
      du   = ur - ul
      dv   = vr - vl
      dw   = wr - wl

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
      alpha_c  = drho - dp/(a_ROE*a_ROE)
      !alpha_s1 = rho_ROE*(-ny*du + nx*dv)
      alpha_s1 = rho_ROE*(du*t1x + dv*t1y + dw*t1z)
      alpha_s2 = rho_ROE*(du*t2x + dv*t2y + dw*t2z)

      ! Antidiffusion coefficient
      delta = a_ROE/(a_ROE + abs(0.5d0*(S1+S4)))

      ! Prefactor
      C = (S1*S4)/(S4-S1)

      ! Mass
      F_r = F_r - C*delta*alpha_c

      ! Momentum
      ! Contact contribution
      F_u = F_u - C*delta*( alpha_c * u_ROE )
      F_v = F_v - C*delta*( alpha_c * v_ROE )
      F_w = F_w - C*delta*( alpha_c * w_ROE )
      ! Shear contributions
      ! F_u = F_u - C*delta*( -ny*alpha_s1 )
      ! F_v = F_v - C*delta*( nx*alpha_s1 )
      F_u = F_u - C*delta*( alpha_s1*t1x + alpha_s2*t2x )
      F_v = F_v - C*delta*( alpha_s1*t1y + alpha_s2*t2y )
      F_w = F_w - C*delta*( alpha_s1*t1z + alpha_s2*t2z )

      ! Energy
      ! Contact contribution
      F_E = F_E - C*delta*( alpha_c*0.5*(u_ROE*u_ROE + v_ROE*v_ROE + w_ROE*w_ROE) )
      !F_E = F_E - C*delta*( alpha_s1*(-ny*u_ROE + nx*v_ROE) )
      ! Shear contributions
      F_E = F_E - C*delta*( alpha_s1*(u_ROE*t1x+v_ROE*t1y+w_ROE*t1z) &
                          + alpha_s2*(u_ROE*t2x+v_ROE*t2y+w_ROE*t2z) )

    endif

  end subroutine riemann_HLLEM


  subroutine riemann_HLLC(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    implicit none
    ! Inputs
    real(R8), intent(in) :: dl(nsc), dr(nsc)
    real(R8), intent(in) :: ul, vl, wl, ur, vr, wr
    real(R8), intent(in) :: pl, pr, al, ar
    real(R8), intent(in) :: dltot, drtot
    real(R8), intent(in) :: nx, ny, nz
    real(R8), intent(in) :: dummy
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


  subroutine riemann_HLLE(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,dummy,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    implicit none
    ! Inputs
    real(R8), intent(in) :: dl(nsc), dr(nsc)
    real(R8), intent(in) :: ul, vl, wl, ur, vr, wr
    real(R8), intent(in) :: pl, pr, al, ar
    real(R8), intent(in) :: dltot, drtot
    real(R8), intent(in) :: nx, ny, nz
    real(R8), intent(in) :: dummy
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


  ! HLLE++
  subroutine riemann_HLLEpp(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta
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


  ! HLLEM - HLLE 
  ! Hybrid solver working with a mulidimensional shock detector
  subroutine riemann_HLLEMSD(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! 
    real(R8) :: FrE, FuE, FvE, FwE, FeE
    real(R8) :: FrC, FuC, FvC, FwC, FeC
    real(R8) :: beta_

    beta_ = beta
    ! beta_ = max(beta_,min_beta)

    if (beta_ == 0d0) then
      call riemann_HLLE  (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    elseif (beta_ >= 1d0) then
      call riemann_HLLEM (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    else
      call riemann_HLLE  (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,nx,ny,nz,FrE,FuE,FvE,FwE,FeE)
      call riemann_HLLEM (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,nx,ny,nz,FrC,FuC,FvC,FwC,FeC)
      F_r = beta_*FrC + (1d0-beta_)*FrE
      F_u = beta_*FuC + (1d0-beta_)*FuE
      F_v = beta_*FvC + (1d0-beta_)*FvE
      F_w = beta_*FwC + (1d0-beta_)*FwE
      F_e = beta_*FeC + (1d0-beta_)*FeE
    endif

  end subroutine riemann_HLLEMSD


  ! HLLC - HLLE 
  ! Hybrid solver working with a mulidimensional shock detector
  subroutine riemann_HLLCSD(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! 
    real(R8) :: FrE, FuE, FvE, FwE, FeE
    real(R8) :: FrC, FuC, FvC, FwC, FeC
    real(R8) :: beta_

    beta_ = beta
    ! beta_ = max(beta_,min_beta)

    if (beta_ == 0d0) then
      call riemann_HLLE (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    elseif (beta_ >= 1d0) then
      call riemann_HLLC (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    else
      call riemann_HLLE (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,nx,ny,nz,FrE,FuE,FvE,FwE,FeE)
      call riemann_HLLC (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,0d0,nx,ny,nz,FrC,FuC,FvC,FwC,FeC)
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


  subroutine riemann_HLLCHLLE  (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc

    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al ! : density(s), velocity, pressure and sound velocity of left state
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar ! : density(s), velocity, pressure and sound velocity of right state
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta
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
      call riemann_HLLE(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx1,ny1,nz1,F_rHLLE,F_uHLLE,F_vHLLE,F_wHLLE, F_EHLLE)

      ! chiamata a HLLC in direzione n2
      call riemann_HLLC(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx2,ny2,nz2,F_rHLLC,F_uHLLC,F_vHLLC,F_wHLLC, F_EHLLC)

      ! media pesata dei flussi per ottenere flusso rotato
      F_r=F_rHLLE*alfa1+F_rHLLC*alfa2
      F_u=F_uHLLE*alfa1+F_uHLLC*alfa2
      F_v=F_vHLLE*alfa1+F_vHLLC*alfa2
      F_w=F_wHLLE*alfa1+F_wHLLC*alfa2
      if (F_w /= 0.d0 .and. abs(nz) <= 1.d-10) F_w = 0.d0
      F_E=F_EHLLE*alfa1+F_EHLLC*alfa2
    else
      ! chiamata a HLLC in direzione n2, assumendo n1 tangente alla faccia
      call riemann_HLLC(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx,ny,nz,F_rHLLC,F_uHLLC,F_vHLLC,F_wHLLC, F_EHLLC)
      F_r=F_rHLLC
      F_u=F_uHLLC
      F_v=F_vHLLC
      F_w=F_wHLLC
      F_E=F_EHLLC
    endif

  end subroutine riemann_HLLCHLLE

end module MOSE_Lib_Riemann_HLL