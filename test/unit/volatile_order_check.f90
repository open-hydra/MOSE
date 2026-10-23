!> Unit check of the order of rounding that a VOLATILE variable keeps, on which the geometry of the cells relies
!> (Compute_Metric_Tensor; docs/theory/numerics.md, Cell volume).
!>
!> The release builds of MOSE compile with -ffast-math (gfortran: -fassociative-math; ifx: -fp-model fast), which lets the
!> compiler evaluate a floating-point expression as a mathematically equivalent one, also across statements: with
!> x = q - c1 and a = c2 - c1 in two variables, d = x - a may be evaluated as q - c2. A formula of differences of nodes then
!> rounds at the size of the coordinates instead of the size of the cell. A VOLATILE object is defined and referenced each
!> time as the program says (Fortran 2018, 8.5.19): a difference stored in a VOLATILE variable and read back is rounded once,
!> as written. With q = 1, c1 = 2^54, c2 = 3, read from the command line so that the compiler cannot fold them, the order of
!> the source gives (1 - 2^54) - (3 - 2^54) = -2^54 + (2^54 - 4) = -4 (both differences rounded to even), and q - c2 = -2.
!> The check asserts -4 through VOLATILE variables, a scalar pair and an array as Compute_Metric_Tensor takes its node
!> differences. The value through plain variables is printed, not asserted: -2 with the release flags of gfortran 13 and
!> ifx 2024.1, -4 without -ffast-math, and both are allowed.
program volatile_order_check
  use iso_fortran_env, only: R8 => real64
  implicit none
  real(R8) :: q, c1, c2, dplain, dvol, dvec(3)
  character(len=64) :: arg
  logical :: ok

  call get_command_argument ( 1, arg )
  read(arg,*) q
  call get_command_argument ( 2, arg )
  read(arg,*) c1
  call get_command_argument ( 3, arg )
  read(arg,*) c2
  dplain = Plain ( q, c1, c2 )
  dvol = Washed ( q, c1, c2 )
  dvec = Washed_Array ( [ q, q, q ], [ c1, c1, c1 ], [ c2, c2, c2 ] )
  write(*,'(A,3ES24.16)') ' q, c1, c2:                      ', q, c1, c2
  write(*,'(A,ES24.16)')  ' (q - c1) - (c2 - c1), plain:    ', dplain
  write(*,'(A,ES24.16)')  '                       VOLATILE: ', dvol
  write(*,'(A,3ES24.16)') '         VOLATILE array:         ', dvec
  ok = dvol == -4.0_R8 .and. all ( dvec == -4.0_R8 )
  if ( ok ) then
    write(*,'(A)') 'VolatileOrder: all checks passed'
  else
    write(*,'(A)') 'VolatileOrder: FAILED (a difference read back from a VOLATILE variable was not rounded as written)'
    error stop 1
  end if

contains

  function Plain ( q, c1, c2 ) result ( d )
    real(R8), intent(in) :: q, c1, c2
    real(R8) :: d, x, a
    x = q - c1
    a = c2 - c1
    d = x - a
  end function Plain

  function Washed ( q, c1, c2 ) result ( d )
    real(R8), intent(in) :: q, c1, c2
    real(R8) :: d, x, a
    real(R8), volatile :: vx, va
    vx = q - c1
    va = c2 - c1
    x = vx
    a = va
    d = x - a
  end function Washed

  function Washed_Array ( q, c1, c2 ) result ( d )
    real(R8), intent(in) :: q(3), c1(3), c2(3)
    real(R8) :: d(3), r(3,2)
    real(R8), volatile :: rv(3,2)
    rv(:,1) = q - c1
    rv(:,2) = c2 - c1
    r = rv
    d = r(:,1) - r(:,2)
  end function Washed_Array

end program volatile_order_check
