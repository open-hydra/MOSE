module MOSE_Mod_Riemann
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Assign_Riemann_Solver

  !> Concrete riemann solver procedure pointing to one of the subroutine realizations
  procedure(riemann_if), pointer, public :: Riemann

  !> Abstract interface relative to the riemann solver procedure
  abstract interface
    subroutine riemann_if(dl,ul,vl,wl,pl,al,dltot,dr,ur,vr,wr,pr,ar,drtot,beta,url,urr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
      use iso_fortran_env, only: I4 => int32, R8 => real64
      use MOSE_Global_m, only: nsc
      use FLINT_Lib_Thermodynamic
      implicit none
      integer s
      real(R8), intent(in)  :: dl(nsc),ul,vl,wl,pl,al ! : density(s), velocity, pressure and sound velocity of left state
      real(R8), intent(in)  :: dr(nsc),ur,vr,wr,pr,ar ! : density(s), velocity, pressure and sound velocity of right state
      real(R8), intent(in)  :: dltot,drtot
      real(R8), intent(in)  :: nx, ny, nz
      real(R8), intent(in)  :: beta
      real(R8), intent(in)  :: url, urr
      real(R8), intent(out) :: F_r, F_u, F_v, F_w, F_e
      ! common
      real(R8) :: Rgasl,Rgasr
    end subroutine riemann_if
  end interface

contains

  subroutine Assign_Riemann_Solver()
    use MOSE_Config_Types_m, only: obj_riemann
    use MOSE_Lib_Riemann_AUSM
    use MOSE_Lib_Riemann_Godunov
    use MOSE_Lib_Riemann_HLL
    use MOSE_Lib_Riemann_LF
    use MOSE_Lib_Riemann_SLAU
    use MOSE_Lib_Riemann_Roe
    implicit none

    nullify(Riemann)

    select case (obj_riemann%description)

      !! AUSM-type solvers
      case ('AUSM+')
        Riemann => riemann_AUSMp
        obj_riemann%description = 'AUSM+'
      case ('AUSM+M')
        Riemann => riemann_AUSMp_M
        obj_riemann%description = 'AUSM+M (Chen et al. 2020)'

      !! Godunov solvers
      case ('exact','Exact')
        Riemann => riemann_exact
        obj_riemann%description = 'Exact Solver'

      !! HLL-type solvers
      case ('HLLE')
        Riemann => riemann_HLLE
        obj_riemann%description = 'HLLE'

      case ('HLLC-PC')
        Riemann => riemann_HLLCprec
        obj_riemann%description = 'HLLC-PC'

      case ('HLLC')
        obj_riemann%description = 'HLLC Batten'
        Riemann => riemann_HLLC

      case ('HLLC+Chen')
        Riemann => riemann_HLLCp_Chen
        obj_riemann%description = 'HLLC+ (Chen et al. 2020)'

      case ('HLLC+Tramel')
        Riemann => riemann_HLLCp_Tramel
        obj_riemann%description = 'HLLC+ (Tramel et al.)'

      case('HLLE++')
        Riemann => riemann_HLLEpp
        obj_riemann%description = 'HLLE++ (Tramel et al.)'

      case ('HLLC Rotated')
        Riemann => riemann_HLLCHLLE
        obj_riemann%description = 'Rotated HLLC Batten / HLLE'

      !! Lax-Friedrichs-type solvers
      case ('LLF','Rusanov')
        Riemann => riemann_LLF
        obj_riemann%description = 'Local Lax-Friedrichs (Rusanov)'

      !! SLAU-type solvers
      case ('SLAU')
        Riemann => riemann_SLAU
        obj_riemann%description = 'SLAU'
      case ('SLAU2')
        Riemann => riemann_SLAU2
        obj_riemann%description = 'SLAU2'

      !! Roe-type solvers
      case ('LMRoe')
        Riemann => riemann_LMRoe
        obj_riemann%description = 'Low-Mach Roe'

      case ('MiczekRoe','Miczek')
        Riemann => riemann_MiczekRoe
        obj_riemann%description = 'Miczek preconditioned Roe (Miczek et al. 2015)'

      !! Default
      case default
        Riemann => riemann_HLLC
        obj_riemann%description = 'HLLC Batten'

    end select

  end subroutine Assign_Riemann_Solver

end module MOSE_Mod_Riemann