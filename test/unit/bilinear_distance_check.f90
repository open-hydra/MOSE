!> Unit check of the wall distance to a bilinear face (docs/theory/turbulence.md, Wall distance): the certified minimum
!> Bilinear_Face_Distance and the fast-path test Face_Fast_Test of MOSE_Lib_Metrics, on the cases of
!> test/unit/bilinear_distance_cases.txt (written by bilinear_distance_cases.py, whose header lists them; every coordinate
!> exact, every reference exact, in rational arithmetic):
!>   D lines: the distance from a point to a face (twisted, plane, non-convex, degenerate, far from the origin, at a focal
!>     point, adverse cases) against the exact distance of the stationarity equations: within tau + 32 u S = 288 u S,
!>     S = |p - c1| + |c2 - c1| + |c4 - c1| + |T|, u = 2^-53 (the proven bound when no piece is cut at the deepest level);
!>     the certificate at most tau = 256 u S, no piece cut at the deepest level, at most 1e5 pieces;
!>   F lines: the fast path taken by the plane faces (convex, thin tapered up to aspect 1e5, triangles of a collapsed edge
!>     also far from the origin, parallelograms and 2-D extruded faces; coplanar to the bit and rounded) and not by the
!>     twisted faces and the plane darts; the upper bound ghat of the gap of the four triangles at least the exact gap.
!> The file is the first argument.
program bilinear_distance_check
  use iso_fortran_env,  only: R8 => real64
  use MOSE_Lib_Metrics, only: Bilinear_Face_Distance, Face_Fast_Test
  implicit none
  real(R8), parameter :: u = epsilon(1.0_R8) / 2.0_R8
  character(len=1024) :: path, line
  character(len=1)    :: tag, want
  character(len=96)   :: label
  real(R8) :: p(3), c(3,4), dref, d, cert, s, t(3), err, gref, ghat
  real(R8) :: worst_err, worst_cert
  integer  :: ios, unit_, ncut, npiece, nd, nf, nbad, maxpiece, nfast, ncuts
  logical  :: fast, fan_ok, ok

  call get_command_argument ( 1, path )
  open ( newunit=unit_, file=trim(path), status='old', action='read', iostat=ios )
  if ( ios /= 0 ) then
    write(*,'(2A)') 'BilinearDistance: FAILED, cannot open ', trim(path)
    error stop 1
  end if
  nd = 0
  nf = 0
  nbad = 0
  nfast = 0
  ncuts = 0
  maxpiece = 0
  worst_err = 0.0_R8
  worst_cert = 0.0_R8
  do
    read(unit_,'(A)',iostat=ios) line
    if ( ios /= 0 ) exit
    if ( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
    read(line,*) tag
    if ( tag == 'D' ) then
      read(line,*) tag, label, p, c, dref
      call Bilinear_Face_Distance ( p, c, d, cert, ncut, npiece )
      t = c(:,1) - c(:,2) + c(:,3) - c(:,4)
      s = norm2 ( p - c(:,1) ) + norm2 ( c(:,2) - c(:,1) ) + norm2 ( c(:,4) - c(:,1) ) + norm2 ( t )
      err = abs ( d - dref ) / ( u * s )
      ok = err <= 288.0_R8 .and. cert <= 256.0_R8 * u * s .and. ncut == 0 .and. npiece <= 100000
      worst_err = max ( worst_err, err )
      worst_cert = max ( worst_cert, cert / ( u * s ) )
      maxpiece = max ( maxpiece, npiece )
      ncuts = ncuts + ncut
      nd = nd + 1
      if ( .not. ok ) then
        nbad = nbad + 1
        write(*,'(3A,ES24.16,A,ES24.16,A,F9.2,A,F9.2,A,I0,A,I0)') ' FAILED ', trim(label), ': distance', d, ', exact', &
          dref, ', error/(u S)', err, ', certificate/(u S)', cert / ( u * s ), ', pieces ', npiece, ', cut ', ncut
      end if
    else if ( tag == 'F' ) then
      read(line,*) tag, label, c, want, gref
      call Face_Fast_Test ( c, fast, ghat, fan_ok )
      ok = ( fast .eqv. want == 'T' ) .and. ghat >= gref
      nf = nf + 1
      if ( fast ) nfast = nfast + 1
      if ( .not. ok ) then
        nbad = nbad + 1
        write(*,'(3A,L1,3A,ES24.16,A,ES24.16)') ' FAILED ', trim(label), ': fast path ', fast, ' (wanted ', want, &
          '), ghat', ghat, ', exact gap', gref
      end if
    end if
  enddo
  close ( unit_ )

  write(*,'(A,I0,A,F7.2,A,F7.2,A,I0,A,I0)') ' distance: ', nd, ' cases, largest error/(u S) ', worst_err, &
    ', largest certificate/(u S) ', worst_cert, ', most pieces ', maxpiece, ', pieces cut at the deepest level ', ncuts
  write(*,'(A,I0,A,I0,A)') ' fast path: ', nf, ' faces, ', nfast, ' taking it'
  if ( nd < 200 .or. nf < 100 ) then
    write(*,'(A)') 'BilinearDistance: FAILED (fewer cases than written)'
    error stop 1
  end if
  if ( nbad == 0 ) then
    write(*,'(A)') 'BilinearDistance: all checks passed'
  else
    write(*,'(A,I0,A)') 'BilinearDistance: FAILED (', nbad, ' cases)'
    error stop 1
  end if
end program bilinear_distance_check
