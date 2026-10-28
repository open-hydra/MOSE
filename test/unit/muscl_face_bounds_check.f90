!> Unit check of the MUSCL reconstruction on lines of cells of any lengths: Reconstruction (MOSE_Lib_Reconstruction) with the
!> limiters of MOSE_Lib_Limiters, called with the distances that Convective_Flux passes (l0, l1, l2 between the centres of the
!> cells, dll and dlr the half lengths of the two cells of the face).
!>
!> A face value outside the values of the two cells of its face is a new extremum; for a variable that is zero in the next
!> cell it is negative. The bounds of the limited slope keep both face values of a cell between the cell and its neighbours on
!> any line (Berger, Aftosmis and Murman, AIAA 2005-0490, eq. 19). Checks:
!>   - three stencils, with van Leer, MC and superbee: A, the cell before the face 1.25 times the cell after it, values 1,
!>     0.01, 0, 0 (left state); B, its mirror image (right state); C, one forward Euler step of the advection to the right at
!>     Courant number 1/2 of the middle cell of h = 0.01, 0.01, 1, 1, 1 and w = 0, 0, 1, 1000, 1000. A and B within the values
!>     of the two cells of the face, C not below zero, within 8 units in the last place (ulp, spacing() of the largest value);
!>   - for every limiter, on 50000 random stencils (lengths with ratios up to 100, values over 14 orders of magnitude and zero,
!>     of one sign or of both, beta = 1 and beta < 1), and on 50000 with values from 1e-33 to 1e-8 times the largest length:
!>     every face value within 8 ulp of the interval of its two cells (van Albada only where the ratio of the slopes is
!>     positive; values below 1e-13 times the largest length are left out: the 1e-30 that the limiters add to the slope moves
!>     their face value by up to about 1e-29 times the length);
!>   - with the value of the cell zero, its length a power of 2 and beta = 1, the face value is the limited slope times the
!>     half length, exactly: the slope has the sign of the limiter's own slope (the call with the bounds 2 and 2 of a line of
!>     equal cells) and is not larger in modulus; minmod and first order take no bounds and give that slope; the two faces of a
!>     cell take the same slope; on a line of equal cells, and for data linear along the line, the slope is the limiter's own;
!>   - the call with the bounds 2 and 2, which the faces of symmetry planes, of walls and of the inlets 405 and 420 make,
!>     gives the limited slope of the formulas of the limiters, copied below. This check compares two compilations of the same
!>     formulas, the library and this program: it is asserted when the build keeps the order of evaluation of the source,
!>     which the three numbers of the command line show as VolatileOrder does ((q - c1) - (c2 - c1) with q = 1, c1 = 2**54,
!>     c2 = 3 gives -4, and -2 when the compiler evaluates it as q - c2, as -ffast-math lets it); with such optimizations the
!>     two compilations need not agree in the last bits, and the differences are printed, not asserted.
!> Values are compared as numbers (a zero of either sign is the same). The random numbers come from the minimal standard
!> generator of Park and Miller with a fixed seed, drawn in the same order in every build.
program muscl_face_bounds_check
  use iso_fortran_env, only: R8 => real64, I8 => int64
  use MOSE_Lib_Limiters
  use MOSE_Lib_Reconstruction, only: Reconstruction
  implicit none
  integer, parameter :: NLIM = 6, IORD = 1, MINMOD = 2, VANLEER = 3, VANALBADA = 4, MC = 5, SB = 6
  character(len=9), parameter :: lname(NLIM) = &
    [character(len=9) :: 'first', 'minmod', 'vanleer', 'vanalbada', 'mc', 'superbee']
  integer, parameter :: NS = 50000
  real(R8), parameter :: ULP_TOL = 8.0_R8, TINY_REL = 1.0e-13_R8
  integer(I8) :: seed
  integer :: nfail, lim
  logical :: source_order

  nfail = 0
  seed = 20261028_I8
  source_order = Source_Order_Kept ()

  call Red_Cases ()
  do lim = 1, NLIM
    call Point ( lim )
    call Face_Values ( lim, .false. )
    call Face_Values ( lim, .true. )
    call Slopes ( lim )
    call Same_Slope ( lim )
    call Equal_Cells ( lim )
    call Linear_Data ( lim )
    call Equal_Bounds ( lim )
  end do

  if ( nfail == 0 ) then
    write(*,'(A)') 'MusclFaceBounds: all checks passed'
  else
    write(*,'(A,I0,A)') 'MusclFaceBounds: FAILED (', nfail, ' checks)'
    error stop 1
  end if

contains

  !> Whether the build keeps the order of evaluation of the source: (q - c1) - (c2 - c1) through plain variables, with q, c1
  !> and c2 read from the command line (1, 2**54 and 3: -4 in the order of the source, -2 as q - c2); without the three numbers,
  !> taken as kept
  function Source_Order_Kept () result ( kept )
    logical :: kept
    character(len=64) :: arg
    real(R8) :: q, c1, c2, x, a, d
    kept = .true.
    if ( command_argument_count() < 3 ) return
    call get_command_argument ( 1, arg )
    read(arg,*) q
    call get_command_argument ( 2, arg )
    read(arg,*) c1
    call get_command_argument ( 3, arg )
    read(arg,*) c2
    x = q - c1
    a = c2 - c1
    d = x - a
    kept = ( d == -4.0_R8 )
    write(*,'(A,ES24.16,A,L1)') ' (q - c1) - (c2 - c1) =', d, ': order of the source kept ', kept
  end function Source_Order_Kept

  !> Uniform random number in (0, 1): minimal standard generator (Park and Miller 1988), 16807*seed mod (2**31 - 1)
  function urand () result ( u )
    real(R8) :: u
    seed = mod ( 16807_I8 * seed, 2147483647_I8 )
    u = real ( seed, R8 ) / 2147483647.0_R8
  end function urand

  subroutine Point ( l )
    integer, intent(in) :: l
    select case ( l )
      case ( IORD )      ; rlimiter => rlimiter_IORD
      case ( MINMOD )    ; rlimiter => rlimiter_MINMOD
      case ( VANLEER )   ; rlimiter => rlimiter_VANLEER
      case ( VANALBADA ) ; rlimiter => rlimiter_VANALBADA
      case ( MC )        ; rlimiter => rlimiter_MC
      case ( SB )        ; rlimiter => rlimiter_SB
    end select
  end subroutine Point

  !> The two states of the face between cells 2 and 3 of a line of four cells of values w and lengths h, with the distances
  !> of Convective_Flux
  subroutine Faces ( w, h, beta, priml, primr )
    real(R8), intent(in)  :: w(4), h(4), beta
    real(R8), intent(out) :: priml, primr
    priml = 0.0_R8
    primr = 0.0_R8
    call Reconstruction ( w(1), w(2), w(3), w(4), 0.5d0 * ( h(1) + h(2) ), 0.5d0 * ( h(2) + h(3) ), &
                          0.5d0 * ( h(3) + h(4) ), 0.5d0 * h(2), 0.5d0 * h(3), beta, priml, primr )
  end subroutine Faces

  !> Excess of s over the interval of a and b, in units of spacing(max(|a|, |b|))
  function Excess ( s, a, b ) result ( e )
    real(R8), intent(in) :: s, a, b
    real(R8) :: e
    e = max ( min ( a, b ) - s, s - max ( a, b ), 0.0_R8 ) / spacing ( max ( abs ( a ), abs ( b ) ) )
  end function Excess

  subroutine Fail ( text )
    character(len=*), intent(in) :: text
    nfail = nfail + 1
    if ( nfail <= 40 ) write(*,'(A)') '  FAIL '//text
  end subroutine Fail

  !> A random line of n cells: lengths 10**U(-1,1); values 0 one time in ten, else 10**U(-7,7) or, with tiny, 10**U(-33,-8)
  !> times the largest length, of one sign or of both; beta 1 or U(0,1)
  subroutine Random_Line ( n, tiny, w, h, beta )
    integer, intent(in) :: n
    logical, intent(in) :: tiny
    real(R8), intent(out) :: w(n), h(n), beta
    integer :: k
    logical :: one_sign
    real(R8) :: u
    one_sign = urand() < 0.5_R8
    do k = 1, n
      h(k) = 10.0_R8**( 2.0_R8 * urand() - 1.0_R8 )
    end do
    do k = 1, n
      if ( urand() < 0.1_R8 ) then
        w(k) = 0.0_R8
      else if ( tiny ) then
        w(k) = 10.0_R8**( 25.0_R8 * urand() - 33.0_R8 ) * maxval ( h )
      else
        w(k) = 10.0_R8**( 14.0_R8 * urand() - 7.0_R8 )
      end if
      u = urand()
      if ( .not. one_sign .and. u < 0.5_R8 ) w(k) = -w(k)
    end do
    beta = 1.0_R8
    u = urand()
    if ( u < 0.5_R8 ) beta = urand()
  end subroutine Random_Line

  !> A length that is a power of 2, from 1/8 to 8
  function Power_Of_Two () result ( h )
    real(R8) :: h
    h = 2.0_R8**( int ( 7.0_R8 * urand() ) - 3 )
  end function Power_Of_Two

  subroutine Red_Cases ()
    integer, parameter :: lims(3) = [ VANLEER, MC, SB ]
    real(R8) :: pl, pr, wm, wp, wn, dummy, e
    integer :: k
    do k = 1, 3
      call Point ( lims(k) )
      ! A: left state of the face between the cells of values 0.01 and 0, the first 1.25 times longer
      call Faces ( [1.0_R8, 0.01_R8, 0.0_R8, 0.0_R8], [1.25_R8, 1.25_R8, 1.0_R8, 1.0_R8], 1.0_R8, pl, pr )
      e = Excess ( pl, 0.01_R8, 0.0_R8 )
      write(*,'(A,A9,A,ES16.8,A,G0.4,A)') ' A ', lname(lims(k)), '  left state  ', pl, ', excess over [0, 0.01] ', e, ' ulp'
      if ( e > ULP_TOL ) call Fail ( 'A, '//trim(lname(lims(k)))//': left state outside the values of its two cells' )
      ! B: the mirror image of A, right state
      call Faces ( [0.0_R8, 0.0_R8, 0.01_R8, 1.0_R8], [1.0_R8, 1.0_R8, 1.25_R8, 1.25_R8], 1.0_R8, pl, pr )
      e = Excess ( pr, 0.0_R8, 0.01_R8 )
      write(*,'(A,A9,A,ES16.8,A,G0.4,A)') ' B ', lname(lims(k)), '  right state ', pr, ', excess over [0, 0.01] ', e, ' ulp'
      if ( e > ULP_TOL ) call Fail ( 'B, '//trim(lname(lims(k)))//': right state outside the values of its two cells' )
      ! C: forward Euler step of the cell of length 1 and value 1 next to cells of 0.01, advection to the right (upwind:
      ! the left states), Courant number 1/2: w = 1 - (1/2) (state at its right face - state of the cell before at its left face)
      call Faces ( [0.0_R8, 0.0_R8, 1.0_R8, 1000.0_R8], [0.01_R8, 0.01_R8, 1.0_R8, 1.0_R8], 1.0_R8, wm, dummy )
      call Faces ( [0.0_R8, 1.0_R8, 1000.0_R8, 1000.0_R8], [0.01_R8, 1.0_R8, 1.0_R8, 1.0_R8], 1.0_R8, wp, dummy )
      wn = 1.0_R8 - 0.5_R8 * ( wp - wm )
      write(*,'(A,A9,A,ES16.8)') ' C ', lname(lims(k)), '  value after the step ', wn
      if ( wn < -ULP_TOL * spacing ( 1000.0_R8 ) ) call Fail ( 'C, '//trim(lname(lims(k)))//': the step makes the cell negative' )
    end do
  end subroutine Red_Cases

  !> Every face value within 8 ulp of the values of its two cells
  subroutine Face_Values ( l, tiny )
    integer, intent(in) :: l
    logical, intent(in) :: tiny
    real(R8) :: w(4), h(4), beta, pl, pr, e, emax, s0, s1, s2
    integer :: n, nin, nout
    character(len=160) :: text
    emax = 0.0_R8
    nin = 0
    nout = 0
    do n = 1, NS
      call Random_Line ( 4, tiny, w, h, beta )
      call Faces ( w, h, beta, pl, pr )
      s0 = ( w(2) - w(1) ) / ( 0.5d0 * ( h(1) + h(2) ) )
      s1 = ( w(3) - w(2) ) / ( 0.5d0 * ( h(2) + h(3) ) )
      s2 = ( w(4) - w(3) ) / ( 0.5d0 * ( h(3) + h(4) ) )
      if ( max ( abs ( w(2) ), abs ( w(3) ) ) < TINY_REL * maxval ( h ) ) then
        nout = nout + 1
        cycle
      end if
      if ( l /= VANALBADA .or. s0 / ( s1 + 1d-30 ) > 0.0_R8 ) then
        nin = nin + 1
        e = Excess ( pl, w(2), w(3) )
        emax = max ( emax, e )
        if ( e > ULP_TOL ) then
          write(text,'(A,4ES11.3,A,4ES11.3,A,ES11.3)') 'w', w, ' h', h, ' left state', pl
          call Fail ( trim(lname(l))//': '//trim(text) )
        end if
      end if
      if ( l /= VANALBADA .or. s1 / ( s2 + 1d-30 ) > 0.0_R8 ) then
        nin = nin + 1
        e = Excess ( pr, w(2), w(3) )
        emax = max ( emax, e )
        if ( e > ULP_TOL ) then
          write(text,'(A,4ES11.3,A,4ES11.3,A,ES11.3)') 'w', w, ' h', h, ' right state', pr
          call Fail ( trim(lname(l))//': '//trim(text) )
        end if
      end if
    end do
    if ( tiny ) then
      text = 'values from 1e-33 to 1e-8 times the largest length'
    else
      text = 'values from 1e-7 to 1e7'
    end if
    write(*,'(1X,A9,A,A,A,I0,A,G0.4,A,I0,A)') lname(l), ' face values, ', trim(text), ': ', nin, &
      ' checked, largest excess over the values of their two cells ', emax, ' ulp; ', nout, &
      ' stencils with values below 1e-13 times the largest length left out'
  end subroutine Face_Values

  !> The limited slope, read from the face value of a cell of value 0, has the sign of the limiter's own slope (the bounds 2
  !> and 2 of a line of equal cells) and is not larger in modulus; minmod and first order take no bounds
  subroutine Slopes ( l )
    integer, intent(in) :: l
    real(R8) :: w(4), h(4), beta, pl, pr, sig, sun, s0, s1, s2
    integer :: n, nsmall
    nsmall = 0
    do n = 1, NS
      call Random_Line ( 4, mod ( n, 2 ) == 0, w, h, beta )
      if ( mod ( n, 4 ) < 2 ) then
        w(2) = 0.0_R8
        h(2) = Power_Of_Two ()
        call Faces ( w, h, 1.0_R8, pl, pr )
        sig = pl / ( 0.5d0 * h(2) )
        s0 = ( w(2) - w(1) ) / ( 0.5d0 * ( h(1) + h(2) ) )
        s1 = ( w(3) - w(2) ) / ( 0.5d0 * ( h(2) + h(3) ) )
        sun = rlimiter ( s1, s0, 2d0, 2d0 )
      else
        w(3) = 0.0_R8
        h(3) = Power_Of_Two ()
        call Faces ( w, h, 1.0_R8, pl, pr )
        sig = - pr / ( 0.5d0 * h(3) )
        s1 = ( w(3) - w(2) ) / ( 0.5d0 * ( h(2) + h(3) ) )
        s2 = ( w(4) - w(3) ) / ( 0.5d0 * ( h(3) + h(4) ) )
        sun = rlimiter ( s2, s1, 2d0, 2d0 )
      end if
      if ( abs ( sig ) < abs ( sun ) ) nsmall = nsmall + 1
      if ( abs ( sig ) > abs ( sun ) .or. ( sig /= 0.0_R8 .and. ( ( sig > 0.0_R8 ) .neqv. ( sun > 0.0_R8 ) ) ) ) &
        call Fail ( trim(lname(l))//': limited slope of another sign or larger than the limiter''s own' )
      if ( ( l == MINMOD .or. l == IORD ) .and. sig /= sun ) call Fail ( trim(lname(l))//': slope changed by the bounds' )
    end do
    write(*,'(1X,A9,A,I0,A,I0,A)') lname(l), ' slopes: ', NS, ' with the sign of the limiter''s own slope and not larger, ', &
      nsmall, ' smaller'
  end subroutine Slopes

  !> The two faces of a cell take the same limited slope: cell 3 of five, of value 0 and length a power of 2
  subroutine Same_Slope ( l )
    integer, intent(in) :: l
    real(R8) :: w(5), h(5), beta, pla, pra, plb, prb
    integer :: n
    do n = 1, NS
      call Random_Line ( 5, mod ( n, 2 ) == 0, w, h, beta )
      w(3) = 0.0_R8
      h(3) = Power_Of_Two ()
      call Faces ( w(1:4), h(1:4), 1.0_R8, pla, pra )
      call Faces ( w(2:5), h(2:5), 1.0_R8, plb, prb )
      if ( plb /= - pra ) call Fail ( trim(lname(l))//': the two faces of a cell take different slopes' )
    end do
    write(*,'(1X,A9,A,I0,A)') lname(l), ' one slope per cell: ', NS, ' cells'
  end subroutine Same_Slope

  !> On a line of equal cells the limited slope is the limiter's own
  subroutine Equal_Cells ( l )
    integer, intent(in) :: l
    real(R8) :: w(4), h(4), beta, pl, pr, sig, sun, s0, s1, s2, hh
    integer :: n
    do n = 1, NS
      call Random_Line ( 4, mod ( n, 2 ) == 0, w, h, beta )
      if ( mod ( n, 3 ) < 2 ) then
        ! lengths a power of 2: the slope is read exactly from the face value
        hh = Power_Of_Two ()
        h = hh
        if ( mod ( n, 3 ) == 0 ) then
          w(2) = 0.0_R8
          call Faces ( w, h, 1.0_R8, pl, pr )
          sig = pl / ( 0.5d0 * hh )
          s0 = ( w(2) - w(1) ) / ( 0.5d0 * ( h(1) + h(2) ) )
          s1 = ( w(3) - w(2) ) / ( 0.5d0 * ( h(2) + h(3) ) )
          sun = rlimiter ( s1, s0, 2d0, 2d0 )
        else
          w(3) = 0.0_R8
          call Faces ( w, h, 1.0_R8, pl, pr )
          sig = - pr / ( 0.5d0 * hh )
          s1 = ( w(3) - w(2) ) / ( 0.5d0 * ( h(2) + h(3) ) )
          s2 = ( w(4) - w(3) ) / ( 0.5d0 * ( h(3) + h(4) ) )
          sun = rlimiter ( s2, s1, 2d0, 2d0 )
        end if
        if ( sig /= sun ) call Fail ( trim(lname(l))//': equal cells, slope changed by the bounds' )
      else
        ! any length: the face value of a cell of value 0 is the limiter's own slope times the half length, rounded once
        hh = h(1)
        h = hh
        w(2) = 0.0_R8
        call Faces ( w, h, 1.0_R8, pl, pr )
        s0 = ( w(2) - w(1) ) / ( 0.5d0 * ( h(1) + h(2) ) )
        s1 = ( w(3) - w(2) ) / ( 0.5d0 * ( h(2) + h(3) ) )
        sun = rlimiter ( s1, s0, 2d0, 2d0 )
        if ( pl /= sun * ( 0.5d0 * hh ) ) call Fail ( trim(lname(l))//': equal cells, face value changed by the bounds' )
      end if
    end do
    write(*,'(1X,A9,A,I0,A)') lname(l), ' equal cells: ', NS, ' stencils with the limiter''s own slope'
  end subroutine Equal_Cells

  !> Data linear along a line of cells of any lengths: the bounds do not act
  subroutine Linear_Data ( l )
    integer, intent(in) :: l
    real(R8) :: w(4), h(4), beta, x(4), c, pl, pr, sig, sun, s0, s1
    integer :: n, k
    do n = 1, NS
      call Random_Line ( 4, .false., w, h, beta )
      h(2) = Power_Of_Two ()
      x(1) = 0.0_R8
      do k = 2, 4
        x(k) = x(k-1) + 0.5d0 * ( h(k-1) + h(k) )
      end do
      c = 10.0_R8**( 14.0_R8 * urand() - 7.0_R8 )
      if ( urand() < 0.5_R8 ) c = -c
      w = c * ( x - x(2) )
      w(2) = 0.0_R8
      call Faces ( w, h, 1.0_R8, pl, pr )
      sig = pl / ( 0.5d0 * h(2) )
      s0 = ( w(2) - w(1) ) / ( 0.5d0 * ( h(1) + h(2) ) )
      s1 = ( w(3) - w(2) ) / ( 0.5d0 * ( h(2) + h(3) ) )
      sun = rlimiter ( s1, s0, 2d0, 2d0 )
      if ( sig /= sun ) call Fail ( trim(lname(l))//': linear data, slope changed by the bounds' )
    end do
    write(*,'(1X,A9,A,I0,A)') lname(l), ' linear data: ', NS, ' stencils with the limiter''s own slope'
  end subroutine Linear_Data

  !> The call with the bounds 2 and 2 of a line of equal cells gives the slope of the formulas of the limiters (asserted when
  !> the build keeps the order of evaluation of the source)
  subroutine Equal_Bounds ( l )
    integer, intent(in) :: l
    real(R8) :: x, y, s, f, dmax
    integer :: n, ndiff
    ndiff = 0
    dmax = 0.0_R8
    do n = 1, NS
      x = 10.0_R8**( 60.0_R8 * urand() - 30.0_R8 )
      y = 10.0_R8**( 60.0_R8 * urand() - 30.0_R8 )
      if ( urand() < 0.5_R8 ) x = -x
      if ( urand() < 0.5_R8 ) y = -y
      if ( urand() < 0.05_R8 ) x = 0.0_R8
      if ( urand() < 0.05_R8 ) y = 0.0_R8
      s = rlimiter ( x, y, 2d0, 2d0 )
      f = Formula ( l, x, y )
      if ( s /= f ) then
        ndiff = ndiff + 1
        dmax = max ( dmax, abs ( s - f ) / spacing ( max ( abs ( s ), abs ( f ) ) ) )
        if ( source_order ) call Fail ( trim(lname(l))//': the call with the bounds 2 and 2 differs from the formula of the limiter' )
      end if
    end do
    if ( source_order .or. ndiff == 0 ) then
      write(*,'(1X,A9,A,I0,A,I0,A)') lname(l), ' bounds 2 and 2: ', NS - ndiff, ' of ', NS, &
        ' pairs with the slope of the formula of the limiter'
    else
      write(*,'(1X,A9,A,I0,A,I0,A,G0.4,A)') lname(l), ' bounds 2 and 2: ', ndiff, ' of ', NS, &
        ' pairs differ from the formula of the limiter compiled in this program, by up to ', dmax, &
        ' ulp (order of evaluation not kept: printed, not asserted)'
    end if
  end subroutine Equal_Bounds

  !> The limited slope of the formulas of the limiters (copied from MOSE_Lib_Limiters, without bounds)
  function Formula ( l, x, y ) result ( s )
    integer, intent(in) :: l
    real(R8), intent(in) :: x, y
    real(R8) :: s, ru, phi
    select case ( l )
      case ( IORD )
        s = 0.0d0
      case ( MINMOD )
        s = 0.d0
        if (x*y<0.d0) return
        s = x
        if (abs(x)>abs(y)) s = y
      case ( VANLEER )
        ru = y/(x + 1d-30)
        phi = (ru+abs(ru))/(1d0+abs(ru))
        s = phi*x
      case ( VANALBADA )
        ru = y/(x+1d-30)
        phi = (ru*ru+ru)/(1d0+ru*ru)
        s = phi*x
      case ( MC )
        ru = y/(x+1d-30)
        phi = maxval([0d0,minval([2d0*ru,0.5d0*(1d0+ru),2d0])])
        s = phi*x
      case ( SB )
        ru = y/(x+1d-30)
        phi = maxval([0d0,min(2d0*ru,1d0),min(ru,2d0)])
        s = phi*x
      case default
        s = 0.0d0
    end select
  end function Formula

end program muscl_face_bounds_check
