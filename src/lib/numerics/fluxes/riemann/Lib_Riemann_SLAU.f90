module MOSE_Lib_Riemann_SLAU
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_SLAU, riemann_SLAU2

contains

  !> @brief Base routine for SLAU-family Riemann solvers (SLAU / SLAU2)
  !>
  !> Computes all common intermediate quantities for the SLAU family:
  !> Mach numbers, Chi, f_rho, Vn weighting, BetaL/BetaR functions, and
  !> the convective (mass) flux. SLAU and SLAU2 differ only in the
  !> pressure flux definition, which is handled in the derived routines.
  !>
  !> @param[in]  dltot,drtot,ul,vl,wl,ur,vr,wr,pl,pr,al,ar,nx,ny,nz
  !> @param[out] F_r,BetaL,BetaR,Chi,sq_veli,sq_velj,aF
  subroutine riemann_SLAU_base(dltot,drtot,ul,vl,wl,ur,vr,wr,pl,pr,al,ar,nx,ny,nz, &
                              F_r,BetaL,BetaR,Chi,sq_veli,sq_velj,aF)
    implicit none
    real(R8), intent(in)  :: dltot,drtot,ul,vl,wl,ur,vr,wr,pl,pr,al,ar,nx,ny,nz
    real(R8), intent(out) :: F_r,BetaL,BetaR,Chi,sq_veli,sq_velj,aF
    ! specific
    real(R8) :: unl, unr, mL, mR, Mach_tilde, f_rho
    real(R8) :: Vn_mag, Vn_magL, Vn_magR

    ! normal velocities and magnitudes
    unl = ul*nx + vl*ny + wl*nz
    unr = ur*nx + vr*ny + wr*nz
    sq_veli = ul**2 + vl**2 + wl**2
    sq_velj = ur**2 + vr**2 + wr**2

    ! mean sound speed
    aF = 0.5d0*(al + ar)

    ! Mach numbers
    mL = unl / aF
    mR = unr / aF

    ! smooth function
    Mach_tilde = min(1.d0, sqrt(0.5d0*(sq_veli + sq_velj)) / aF)
    Chi = (1.d0 - Mach_tilde)**2

    ! density-weighted mean normal velocity magnitude
    Vn_mag = (dltot*abs(unl) + drtot*abs(unr)) / (dltot + drtot)
    f_rho  = -max(min(mL,0.d0), -1.d0) * min(max(mR,0.d0), 1.d0)
    Vn_magL = (1.d0 - f_rho)*Vn_mag + f_rho*abs(unl)
    Vn_magR = (1.d0 - f_rho)*Vn_mag + f_rho*abs(unr)

    ! mass flux
    F_r = 0.5d0 * ( dltot*(unl + Vn_magL) + drtot*(unr - Vn_magR) - (Chi/aF)*(pr - pl) )

    ! BetaL and BetaR functions
    if (abs(mL) < 1.d0) then
      BetaL = 0.25d0*(2.d0 - mL)*(mL + 1.d0)**2
    else if (mL >= 0.d0) then
      BetaL = 1.d0
    else
      BetaL = 0.d0
    endif

    if (abs(mR) < 1.d0) then
      BetaR = 0.25d0*(2.d0 + mR)*(mR - 1.d0)**2
    else if (mR >= 0.d0) then
      BetaR = 0.d0
    else
      BetaR = 1.d0
    endif

  end subroutine riemann_SLAU_base


  subroutine riemann_SLAU(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! Specific
    real(R8) :: BetaL,BetaR,Chi,sq_veli,sq_velj,aF,p_half,H_L,H_R

    call riemann_SLAU_base(dltot,drtot,ul,vl,wl,ur,vr,wr,pl,pr,al,ar,nx,ny,nz,F_r,BetaL,BetaR,Chi,sq_veli,sq_velj,aF)

    ! SLAU pressure flux
    p_half = 0.5d0*(pl + pr) + 0.5d0*(BetaL - BetaR)*(pl - pr) + (1.d0 - Chi)*(BetaL + BetaR - 1.d0)*0.5d0*(pl + pr)

    H_L = H0(pl, dl, sqrt(sq_veli))
    H_R = H0(pr, dr, sqrt(sq_velj))

    if (F_r >= 0.d0) then
      F_u = F_r*ul + nx*p_half
      F_v = F_r*vl + ny*p_half
      F_w = F_r*wl + nz*p_half
      F_E = F_r*H_L
    else
      F_u = F_r*ur + nx*p_half
      F_v = F_r*vr + ny*p_half
      F_w = F_r*wr + nz*p_half
      F_E = F_r*H_R
    endif

  end subroutine riemann_SLAU


  subroutine riemann_SLAU2(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! Specific
    real(R8) :: BetaL,BetaR,Chi,sq_veli,sq_velj,aF,p_half,H_L,H_R

    call riemann_SLAU_base(dltot,drtot,ul,vl,wl,ur,vr,wr,pl,pr,al,ar,nx,ny,nz,F_r,BetaL,BetaR,Chi,sq_veli,sq_velj,aF)

    ! SLAU2 pressure flux (Kitamura & Shima 2013)
    p_half = 0.5d0*(pl + pr) + 0.5d0*(BetaL - BetaR)*(pl - pr) + sqrt(0.5d0*(sq_veli + sq_velj))*(BetaL + BetaR - 1.d0)*aF*0.5d0*(dltot + drtot)

    H_L = H0(pl, dl, sqrt(sq_veli))
    H_R = H0(pr, dr, sqrt(sq_velj))

    if (F_r >= 0.d0) then
      F_u = F_r*ul + nx*p_half
      F_v = F_r*vl + ny*p_half
      F_w = F_r*wl + nz*p_half
      F_E = F_r*H_L
    else
      F_u = F_r*ur + nx*p_half
      F_v = F_r*vr + ny*p_half
      F_w = F_r*wr + nz*p_half
      F_E = F_r*H_R
    endif
    
  end subroutine riemann_SLAU2

end module MOSE_Lib_Riemann_SLAU