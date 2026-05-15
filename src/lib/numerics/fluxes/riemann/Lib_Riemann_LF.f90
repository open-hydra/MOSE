module MOSE_Lib_Riemann_LF
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_PLLF, riemann_LLF

contains

  !> @brief Local Lax-Friedrichs (alias Rusanov)
  subroutine riemann_LLF(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! common
    real(R8) :: Rgasl,Rgasr
    ! specific
    real(R8) :: uln, urn
    real(R8) :: Frr, Fur, Fvr, Fwr, Fer, Frl, Ful, Fvl, Fwl, Fel
    real(R8) :: A

    Rgasl = f_Rtot(dl)
    Rgasr = f_Rtot(dr)

    uln=ul*nx+vl*ny+wl*nz
    urn=ur*nx+vr*ny+wr*nz

    ! Set spectral radius
    A = MAX(ABS(uln-al),ABS(urn-ar),ABS(uln+al),ABS(urn+ar))

    call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,Frr,Fur,Fvr,Fwr,Fer)
    call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,Frl,Ful,Fvl,Fwl,Fel)
    
    F_r = 0.5*(Frl + Frr - A*(drtot - dltot))
    F_u = 0.5*(Fur + Ful - A*(drtot*ur-dltot*ul))
    F_v = 0.5*(Fvr + Fvl - A*(drtot*vr-dltot*vl))
    F_w = 0.5*(Fwr + Fwl - A*(drtot*wr-dltot*wl))
    F_E = 0.5*(Fel + Fer - A*(drtot*E0(pr,dr,sqrt(ur*ur+vr*vr+wr*wr)) - dltot*E0(pl,dl,sqrt(ul*ul+vl*vl+wl*wl))))
        
  end subroutine riemann_LLF

  ! Preconditioned Local Lax-Friedrichs (alias Rusanov)
  subroutine Riemann_PLLF (dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: dl(nsc), ul, vl, wl, pl, al
    real(R8), intent(in)  :: dr(nsc), ur, vr, wr, pr, ar
    real(R8), intent(in)  :: dltot, drtot
    real(R8), intent(in)  :: nx, ny, nz
    real(R8), intent(in)  :: beta
    real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    ! common
    real(R8) :: Rgasl, Rgasr
    ! specific
    real(R8) :: uln, urn, Tl, Tr, rho_Tl, rho_Tr, cpl, cpr, U_rl, U_rr
    real(R8) :: thetal, thetar, alphal, alphar, betal, betar, ul_mod, ur_mod
    real(R8) :: al_mod, ar_mod, A, Frr, Fur, Fvr, Fwr, Fer, Frl, Ful, Fvl, Fwl, Fel

    Rgasl = f_Rtot(dl)
    Rgasr = f_Rtot(dr)

    uln = ul*nx + vl*ny + wl*nz
    urn = ur*nx + vr*ny + wr*nz

    ! Thermodynamic stuff
    Tl = pl / ( dltot * Rgasl )
    Tr = pr / ( drtot * Rgasr )
    rho_Tl = - dltot / Tl
    rho_Tr = - drtot / Tr
    cpl = f_cp ( dl, Tl, dltot )
    cpr = f_cp ( dr, Tr, drtot )

    ! Preconditioning stuff
    U_rl = Max ( 1d-5 * al, Min ( Abs( uln ) , al ) )
    U_rr = Max ( 1d-5 * ar, Min ( Abs( urn ) , ar ) )
    thetal = 1d0 / U_rl**2 - rho_Tl / ( dltot * cpl )
    thetar = 1d0 / U_rr**2 - rho_Tr / ( drtot * cpr )
    betal = 1d0 / al**2
    betar = 1d0 / ar**2
    alphal = 0.5d0 * ( 1d0 - betal * U_rl**2 )
    alphar = 0.5d0 * ( 1d0 - betar * U_rr**2 )
    ul_mod = Abs ( uln * ( 1d0 - alphal ) )
    ur_mod = Abs ( urn * ( 1d0 - alphar ) )
    al_mod = Sqrt ( alphal**2 * uln**2 + U_rl**2 )
    ar_mod = Sqrt ( alphar**2 * urn**2 + U_rr**2 )

    ! Set spectral radius
    A = Max ( ul_mod + al_mod, ur_mod + ar_mod )

    call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,Frr,Fur,Fvr,Fwr,Fer)
    call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,Frl,Ful,Fvl,Fwl,Fel)
    
    F_r = 0.5*(Frl + Frr - A*(drtot - dltot))
    F_u = 0.5*(Fur + Ful - A*(drtot*ur-dltot*ul))
    F_v = 0.5*(Fvr + Fvl - A*(drtot*vr-dltot*vl))
    F_w = 0.5*(Fwr + Fwl - A*(drtot*wr-dltot*wl))
    F_E = 0.5*(Fel + Fer - A*(drtot*E0(pr,dr,sqrt(ur*ur+vr*vr+wr*wr)) - dltot*E0(pl,dl,sqrt(ul*ul+vl*vl+wl*wl))))
        
  end subroutine Riemann_PLLF

  ! Subroutine for computing the conservative fluxes from primitive variables.
  pure subroutine fluxes(p,r,u,v,w,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
    use MOSE_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0
    implicit none
    real(R8), intent(in)::  p           ! Pressure.
    real(R8), intent(in)::  r(nsc)      ! Density.
    real(R8), intent(in)::  u,v,w       ! Velocity.
    real(R8), intent(in)::  nx,ny,nz    ! Normals.
    real(R8), intent(out):: F_r         ! Flux of mass conservation.
    real(R8), intent(out):: F_u,F_v,F_w ! Flux of momentum conservation.
    real(R8), intent(out):: F_E         ! Flux of energy conservation.
    
    F_r = sum(r)*(u*nx+v*ny+w*nz)
    F_u = F_r*u + p*nx
    F_v = F_r*v + p*ny
    F_w = F_r*w + p*nz
    F_E = F_r*H0(p,r,sqrt(u*u+v*v+w*w))
    
  end subroutine fluxes

end module MOSE_Lib_Riemann_LF