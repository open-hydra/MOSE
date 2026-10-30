module MOSE_Lib_Limiters
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

  !> Concrete "rlimiter" procedure pointing to one of the function realizations
  procedure(rlimiter_if), pointer :: rlimiter

  !> Abstract interface relative to the "rlimiter" procedure: the limited slope phi(r)*x of a cell from the slopes x and y
  !> across its two faces along a grid line, r = y/x, with phi bounded for r > 0 by cx and by cy*r, so that the limited slope
  !> is at most cx*|x| and, up to the 1d-30 added to x in ru, cy*|y| (Reconstruction). On a line of equal cells the bounds
  !> are cx = cy = 2, which no limiter exceeds: the slope is the limiter's own.
  abstract interface
    pure function rlimiter_if(x,y,cx,cy) result(result)
      use iso_fortran_env, only: I4 => int32, R8 => real64
      implicit none
      real(R8), intent(in) :: x, y, cx, cy
      real(R8) :: ru
      real(R8) :: phi
      real(R8) :: result
    end function rlimiter_if
  end interface

contains

  !> First order scheme: zero slope, within any bounds
  pure function rlimiter_IORD(x,y,cx,cy) result(result)
    implicit none
    real(R8), intent(in) :: x, y, cx, cy
    real(R8) :: result
    
    result = 0.0d0

  end function rlimiter_IORD


  !> Min-Mod flux limiter: its slope is at most |x| and |y|, within the bounds cx*|x| and cy*|y| of Reconstruction (cx and cy
  !> at least 1), which it therefore does not take
  pure function rlimiter_MINMOD(x,y,cx,cy) result(result)
    implicit none
    real(R8), intent(in) :: x, y, cx, cy
    real(R8) :: result
    
    result = 0.d0
    if (x*y<0.d0) return
    result = x
    if (abs(x)>abs(y)) result = y

  end function rlimiter_MINMOD


  !> Van Leer flux limiter
  pure function rlimiter_VANLEER(x,y,cx,cy) result(result)
    implicit none
    real(R8), intent(in) :: x, y, cx, cy
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    ru = y/(x + 1d-30)
    phi = (ru+abs(ru))/(1d0+abs(ru))
    ! The bounds, for ru > 0: for ru <= 0 both phi and max(cy*ru,0) are zero
    phi = min(phi,cx,max(cy*ru,0d0))
    result = phi*x

  end function rlimiter_VANLEER


  !> Van Albada flux limiter
  pure function rlimiter_VANALBADA(x,y,cx,cy) result(result)
    implicit none
    real(R8), intent(in) :: x, y, cx, cy
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    ru = y/(x+1d-30)
    phi = (ru*ru+ru)/(1d0+ru*ru)
    ! The bounds, for ru > 0 only: this phi is not zero for ru < 0
    if (ru>0d0) phi = min(phi,cx,cy*ru)
    result = phi*x

  end function rlimiter_VANALBADA


  !> MC flux limiter
  pure function rlimiter_MC(x,y,cx,cy) result(result)
    implicit none
    real(R8), intent(in) :: x, y, cx, cy
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,minval([2d0*ru,0.5d0*(1d0+ru),2d0])])
    ! The bounds, for ru > 0: for ru <= 0 both phi and max(cy*ru,0) are zero
    phi = min(phi,cx,max(cy*ru,0d0))
    result = phi*x

  end function rlimiter_MC


  !> Super-bee flux limiter
  pure function rlimiter_SB(x,y,cx,cy) result(result)
    implicit none
    real(R8), intent(in) :: x, y, cx, cy
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,min(2d0*ru,1d0),min(ru,2d0)])
    ! The bounds, for ru > 0: for ru <= 0 both phi and max(cy*ru,0) are zero
    phi = min(phi,cx,max(cy*ru,0d0))
    result = phi*x

  end function rlimiter_SB

end module MOSE_Lib_Limiters