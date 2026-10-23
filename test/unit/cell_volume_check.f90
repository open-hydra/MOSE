!> Unit check of the volume of a cell, Compute_Metric_Tensor (docs/theory/numerics.md, cell volume).
!>
!> The cell of MOSE is the trilinear map of its 8 nodes; its faces are the bilinear faces of their 4 nodes, the faces whose
!> area vectors Compute_Norm_Area takes. Its volume is the integral of det J over the unit cube: the volume for which the
!> cells sum to the volume enclosed by the boundary faces, and the discrete divergence of a field that grows as x is exact
!> (M. Vinokur, NASA CR-177416, 1986, pp. 14-18). Checks:
!>   1. ten cells against their exact volume: cubes with one twisted face (T_n/h = 1 and 5), a valid non-convex cell with six
!>      twisted faces, a cell with a collapsed edge and two twisted faces, a 2-D prism, two wedge cells of 5 degrees (one on
!>      the axis), the first cube and the prism moved by 1e3, and the first cube mirrored, a left-handed cell whose volume is
!>      stored positive. The exact volumes are the integral of det J in rational arithmetic from the binary64 nodes written
!>      below: test/unit/cell_volume_exact.py prints this data block (--table: det J keeps its sign on every cell but the
!>      last). Tolerance 32 u of the volume, u = 2^-53: a few roundings in each of its six terms, which have the sign of the
!>      volume on these cells (no cancellation);
!>   2. a twisted mesh of 8^3 cells: the sum of the cell volumes is the volume enclosed by the boundary faces, 1/3 of the sum
!>      of xbar_f . S_f over them (the divergence theorem on the bilinear faces), within 1e-13 of it;
!>   3. a cell whose volume has not the sign of det A at its centre, an inverted cell that the centre does not show, raises
!>      bad_metric, the flag of a degenerate or inverted cell after which the set-up stops every rank. Its [ERROR] lines are
!>      expected. bad_metric cannot be reset, so this is the last call.
!> The five tetrahedra of the volume before (each face cut along one diagonal, the same face of the cell next door along the
!> other, and each tetrahedron in absolute value) fail check 1 on the cells with twisted faces, and checks 2 and 3.
program cell_volume_check
  use iso_fortran_env,   only: R8 => real64
  use MOSE_Base_Types_m, only: MOSE_vector_3D_type, MOSE_tensor_3D_type
  use MOSE_Lib_Metrics,  only: Compute_Metric_Tensor, bad_metric
  implicit none
  integer,  parameter :: ncell = 11, nm = 8
  real(R8), parameter :: u = epsilon(1.0_R8) / 2.0_R8, tol_cell = 32.0_R8 * u, tol_mesh = 1.0e-13_R8
  character(len=48)   :: names(ncell)
  real(R8) :: nodes(3,8,ncell), vexact(ncell), x(3,0:nm,0:nm,0:nm), p(3,4), s(3), xb(3)
  real(R8) :: vol, err, vsum, vbound, sgn
  type(MOSE_vector_3D_type) :: nd(8), dl
  type(MOSE_tensor_3D_type) :: m
  integer  :: c, k, i, j, l, d1, d2, d3, d, o1, o2, side, a, b, nfail
  logical  :: ok

  ! ---- data printed by test/unit/cell_volume_exact.py
  names(1) = 'cube, one twisted face, T_n/h = 1'
  nodes(:,1,1) = [ 0.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,2,1) = [ 0.0_R8, 0.0_R8, 1.0_R8 ]
  nodes(:,3,1) = [ 0.0_R8, 1.0_R8, 0.0_R8 ]
  nodes(:,4,1) = [ 0.0_R8, 1.0_R8, 1.0_R8 ]
  nodes(:,5,1) = [ 1.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,6,1) = [ 1.0_R8, 0.0_R8, 1.0_R8 ]
  nodes(:,7,1) = [ 1.0_R8, 1.0_R8, 0.0_R8 ]
  nodes(:,8,1) = [ 1.0_R8, 1.0_R8, 2.0_R8 ]
  vexact(1) = 1.250000000000000000000000e+0_R8
  names(2) = 'cube, one twisted face, T_n/h = 5'
  nodes(:,1,2) = [ 0.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,2,2) = [ 0.0_R8, 0.0_R8, 1.0_R8 ]
  nodes(:,3,2) = [ 0.0_R8, 1.0_R8, 0.0_R8 ]
  nodes(:,4,2) = [ 0.0_R8, 1.0_R8, 1.0_R8 ]
  nodes(:,5,2) = [ 1.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,6,2) = [ 1.0_R8, 0.0_R8, 1.0_R8 ]
  nodes(:,7,2) = [ 1.0_R8, 1.0_R8, 0.0_R8 ]
  nodes(:,8,2) = [ 1.0_R8, 1.0_R8, 6.0_R8 ]
  vexact(2) = 2.250000000000000000000000e+0_R8
  names(3) = 'valid non-convex cell, six twisted faces'
  nodes(:,1,3) = [ 0.05_R8, -0.1_R8, 0.0_R8 ]
  nodes(:,2,3) = [ -0.15_R8, 0.05_R8, 1.1_R8 ]
  nodes(:,3,3) = [ 0.1_R8, 1.05_R8, -0.2_R8 ]
  nodes(:,4,3) = [ 0.0_R8, 0.9_R8, 1.15_R8 ]
  nodes(:,5,3) = [ 1.2_R8, 0.1_R8, 0.05_R8 ]
  nodes(:,6,3) = [ 0.95_R8, -0.15_R8, 0.9_R8 ]
  nodes(:,7,3) = [ 0.9_R8, 1.2_R8, 0.1_R8 ]
  nodes(:,8,3) = [ 1.15_R8, 1.0_R8, 0.95_R8 ]
  vexact(3) = 1.152458333333333279708983e+0_R8
  names(4) = 'collapsed edge, two twisted faces'
  nodes(:,1,4) = [ 0.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,2,4) = [ 0.0_R8, 0.0_R8, 1.0_R8 ]
  nodes(:,3,4) = [ 0.0_R8, 1.0_R8, 0.0_R8 ]
  nodes(:,4,4) = [ 0.0_R8, 1.0_R8, 1.0_R8 ]
  nodes(:,5,4) = [ 0.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,6,4) = [ 1.0_R8, 0.0_R8, 1.0_R8 ]
  nodes(:,7,4) = [ 1.0_R8, 1.0_R8, 0.0_R8 ]
  nodes(:,8,4) = [ 1.0_R8, 1.0_R8, 1.5_R8 ]
  vexact(4) = 8.333333333333333333333333e-1_R8
  names(5) = '2-D prism'
  nodes(:,1,5) = [ 0.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,2,5) = [ 0.0_R8, 0.0_R8, 0.5_R8 ]
  nodes(:,3,5) = [ 0.2_R8, 1.2_R8, 0.0_R8 ]
  nodes(:,4,5) = [ 0.2_R8, 1.2_R8, 0.5_R8 ]
  nodes(:,5,5) = [ 1.1_R8, 0.1_R8, 0.0_R8 ]
  nodes(:,6,5) = [ 1.1_R8, 0.1_R8, 0.5_R8 ]
  nodes(:,7,5) = [ 1.3_R8, 0.9_R8, 0.0_R8 ]
  nodes(:,8,5) = [ 1.3_R8, 0.9_R8, 0.5_R8 ]
  vexact(5) = 5.600000000000000184574578e-1_R8
  names(6) = 'wedge cell, 5 degrees'
  nodes(:,1,6) = [ 0.0_R8, 0.7992385772654863_R8, -0.0348955098922688_R8 ]
  nodes(:,2,6) = [ 0.0_R8, 0.7992385772654863_R8, 0.0348955098922688_R8 ]
  nodes(:,3,6) = [ -0.05_R8, 1.0490006326609507_R8, -0.045800356733602804_R8 ]
  nodes(:,4,6) = [ -0.05_R8, 1.0490006326609507_R8, 0.045800356733602804_R8 ]
  nodes(:,5,6) = [ 0.3_R8, 0.7492861661863933_R8, -0.032714540524002_R8 ]
  nodes(:,6,6) = [ 0.3_R8, 0.7492861661863933_R8, 0.032714540524002_R8 ]
  nodes(:,7,6) = [ 0.35_R8, 1.0989530437400437_R8, -0.0479813261018696_R8 ]
  nodes(:,8,6) = [ 0.35_R8, 1.0989530437400437_R8, 0.0479813261018696_R8 ]
  vexact(6) = 8.530368321427044310086077e-3_R8
  names(7) = 'wedge cell on the axis'
  nodes(:,1,7) = [ 0.0_R8, 0.0_R8, -0.0_R8 ]
  nodes(:,2,7) = [ 0.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,3,7) = [ -0.02_R8, 0.20980012653219013_R8, -0.009160071346720559_R8 ]
  nodes(:,4,7) = [ -0.02_R8, 0.20980012653219013_R8, 0.009160071346720559_R8 ]
  nodes(:,5,7) = [ 0.3_R8, 0.0_R8, -0.0_R8 ]
  nodes(:,6,7) = [ 0.3_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,7,7) = [ 0.32_R8, 0.19980964431637158_R8, -0.0087238774730672_R8 ]
  nodes(:,8,7) = [ 0.32_R8, 0.19980964431637158_R8, 0.0087238774730672_R8 ]
  vexact(7) = 5.983532258769225756659491e-4_R8
  names(8) = 'cube of cell 1 moved by 1e3'
  nodes(:,1,8) = [ 1000.0_R8, -1000.0_R8, 1000.0_R8 ]
  nodes(:,2,8) = [ 1000.0_R8, -1000.0_R8, 1001.0_R8 ]
  nodes(:,3,8) = [ 1000.0_R8, -999.0_R8, 1000.0_R8 ]
  nodes(:,4,8) = [ 1000.0_R8, -999.0_R8, 1001.0_R8 ]
  nodes(:,5,8) = [ 1001.0_R8, -1000.0_R8, 1000.0_R8 ]
  nodes(:,6,8) = [ 1001.0_R8, -1000.0_R8, 1001.0_R8 ]
  nodes(:,7,8) = [ 1001.0_R8, -999.0_R8, 1000.0_R8 ]
  nodes(:,8,8) = [ 1001.0_R8, -999.0_R8, 1002.0_R8 ]
  vexact(8) = 1.250000000000000000000000e+0_R8
  names(9) = '2-D prism moved by 1e3'
  nodes(:,1,9) = [ 1000.0_R8, -1000.0_R8, 1000.0_R8 ]
  nodes(:,2,9) = [ 1000.0_R8, -1000.0_R8, 1000.5_R8 ]
  nodes(:,3,9) = [ 1000.2_R8, -998.8_R8, 1000.0_R8 ]
  nodes(:,4,9) = [ 1000.2_R8, -998.8_R8, 1000.5_R8 ]
  nodes(:,5,9) = [ 1001.1_R8, -999.9_R8, 1000.0_R8 ]
  nodes(:,6,9) = [ 1001.1_R8, -999.9_R8, 1000.5_R8 ]
  nodes(:,7,9) = [ 1001.3_R8, -999.1_R8, 1000.0_R8 ]
  nodes(:,8,9) = [ 1001.3_R8, -999.1_R8, 1000.5_R8 ]
  vexact(9) = 5.599999999999846522769076e-1_R8
  names(10) = 'cube of cell 1 mirrored in x (left-handed)'
  nodes(:,1,10) = [ -0.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,2,10) = [ -0.0_R8, 0.0_R8, 1.0_R8 ]
  nodes(:,3,10) = [ -0.0_R8, 1.0_R8, 0.0_R8 ]
  nodes(:,4,10) = [ -0.0_R8, 1.0_R8, 1.0_R8 ]
  nodes(:,5,10) = [ -1.0_R8, 0.0_R8, 0.0_R8 ]
  nodes(:,6,10) = [ -1.0_R8, 0.0_R8, 1.0_R8 ]
  nodes(:,7,10) = [ -1.0_R8, 1.0_R8, 0.0_R8 ]
  nodes(:,8,10) = [ -1.0_R8, 1.0_R8, 2.0_R8 ]
  vexact(10) = -1.250000000000000000000000e+0_R8
  names(11) = 'inverted cell, det A > 0, V < 0'
  nodes(:,1,11) = [ 0.754373_R8, 0.29754_R8, 0.37438_R8 ]
  nodes(:,2,11) = [ -0.767401_R8, -0.746476_R8, 1.405035_R8 ]
  nodes(:,3,11) = [ 0.295925_R8, 1.736358_R8, -0.344624_R8 ]
  nodes(:,4,11) = [ 0.797913_R8, 0.675172_R8, 0.594738_R8 ]
  nodes(:,5,11) = [ 0.344012_R8, 0.293516_R8, -0.571998_R8 ]
  nodes(:,6,11) = [ 0.262633_R8, 0.047749_R8, 1.744357_R8 ]
  nodes(:,7,11) = [ 0.206212_R8, 1.46736_R8, 0.169779_R8 ]
  nodes(:,8,11) = [ 0.79221_R8, 1.491336_R8, 1.641459_R8 ]
  vexact(11) = -2.781463163047100742504219e-1_R8

  ! ---- end of the data

  nfail = 0

  ! 1. cells against their exact volume
  do c = 1, ncell - 1
    do k = 1, 8
      nd(k) % c = nodes(:,k,c)
    end do
    call Compute_Metric_Tensor ( nd(1), nd(2), nd(3), nd(4), nd(5), nd(6), nd(7), nd(8), m, dl, vol )
    err = abs ( vol - abs ( vexact(c) ) ) / abs ( vexact(c) )
    ok = err <= tol_cell .and. .not. bad_metric
    if ( .not. ok ) nfail = nfail + 1
    write(*,'(A,I3,3A,ES24.16,A,ES24.16,A,ES9.2,2A)') ' cell', c, ' (', trim(names(c)), '): volume', vol, ', exact', &
      abs(vexact(c)), ', error/u', err / u, '  ', merge('ok    ', 'FAILED', ok)
  end do

  ! 2. twisted mesh: a turn about z growing with z and with x y, two sine terms, nodes inside moved by up to 1 % of the side
  do l = 0, nm
    do j = 0, nm
      do i = 0, nm
        x(:,i,j,l) = Mesh_Node ( i, j, l )
      end do
    end do
  end do
  vsum = 0.0_R8
  do l = 0, nm - 1
    do j = 0, nm - 1
      do i = 0, nm - 1
        do d1 = 0, 1
          do d2 = 0, 1
            do d3 = 0, 1
              nd(1 + d3 + 2*d2 + 4*d1) % c = x(:,i+d1,j+d2,l+d3)
            end do
          end do
        end do
        call Compute_Metric_Tensor ( nd(1), nd(2), nd(3), nd(4), nd(5), nd(6), nd(7), nd(8), m, dl, vol )
        vsum = vsum + vol
      end do
    end do
  end do
  ! volume enclosed by the boundary faces: corners in cyclic order in the two other index directions, area vector
  ! (diagonal x diagonal)/2 along +d, outward with the sign of the side
  vbound = 0.0_R8
  do d = 1, 3
    o1 = mod(d, 3) + 1
    o2 = mod(d + 1, 3) + 1
    do side = 0, nm, nm
      sgn = merge(-1.0_R8, 1.0_R8, side == 0)
      do b = 0, nm - 1
        do a = 0, nm - 1
          p(:,1) = Face_Node ( d, side, o1, a,     o2, b     )
          p(:,2) = Face_Node ( d, side, o1, a + 1, o2, b     )
          p(:,3) = Face_Node ( d, side, o1, a + 1, o2, b + 1 )
          p(:,4) = Face_Node ( d, side, o1, a,     o2, b + 1 )
          s = 0.5_R8 * Cross ( p(:,3) - p(:,1), p(:,4) - p(:,2) )
          xb = 0.25_R8 * ( p(:,1) + p(:,2) + p(:,3) + p(:,4) )
          vbound = vbound + sgn * dot_product ( xb, s ) / 3.0_R8
        end do
      end do
    end do
  end do
  err = abs ( vsum - vbound ) / vbound
  ok = err <= tol_mesh .and. .not. bad_metric
  if ( .not. ok ) nfail = nfail + 1
  write(*,'(A,ES24.16,A,ES24.16,A,ES9.2,2A)') ' twisted mesh of 8^3 cells: sum of the volumes', vsum, &
    ', volume enclosed by the boundary', vbound, ', relative difference', err, '  ', merge('ok    ', 'FAILED', ok)

  ! 3. inverted cell (det A > 0, V < 0): bad_metric, no stop
  c = ncell
  do k = 1, 8
    nd(k) % c = nodes(:,k,c)
  end do
  write(*,'(3A)') ' cell', ' 11 (', trim(names(c)) // '): the [ERROR] lines below are expected'
  call Compute_Metric_Tensor ( nd(1), nd(2), nd(3), nd(4), nd(5), nd(6), nd(7), nd(8), m, dl, vol )
  ok = bad_metric
  if ( .not. ok ) nfail = nfail + 1
  write(*,'(A,ES24.16,A,L1,2A)') ' cell 11: exact volume', vexact(c), ', bad_metric ', bad_metric, '  ', &
    merge('ok    ', 'FAILED', ok)

  if ( nfail == 0 ) then
    write(*,'(A)') 'CellVolume: all checks passed'
  else
    write(*,'(A,I0,A)') 'CellVolume: FAILED (', nfail, ' checks): the volume is not the integral of det J of the cell'
    error stop 1
  end if

contains

  !> Node (i, j, l) of the twisted mesh of check 2.
  function Mesh_Node ( i, j, l ) result ( xn )
    integer, intent(in) :: i, j, l
    real(R8) :: xn(3), xx, yy, zz, ang
    integer  :: q
    xx = real(i, R8) / nm
    yy = real(j, R8) / nm
    zz = real(l, R8) / nm
    ang = 0.6_R8 * zz + 0.3_R8 * xx * yy
    xn(1) = xx * cos(ang) - yy * sin(ang) + 0.05_R8 * sin(3.0_R8 * yy + 2.0_R8 * zz)
    xn(2) = xx * sin(ang) + yy * cos(ang) + 0.04_R8 * sin(2.0_R8 * xx * zz)
    xn(3) = zz + 0.1_R8 * xx * yy
    if ( i > 0 .and. i < nm .and. j > 0 .and. j < nm .and. l > 0 .and. l < nm ) then
      do q = 1, 3
        xn(q) = xn(q) + 0.01_R8 / nm * sin(12.9898_R8 * i + 78.233_R8 * j + 37.719_R8 * l + 1.3_R8 * q)
      end do
    end if
  end function Mesh_Node

  !> Node of the boundary face of direction d on the node plane side, at index ia along o1 and ib along o2.
  function Face_Node ( d, side, o1, ia, o2, ib ) result ( xn )
    integer, intent(in) :: d, side, o1, ia, o2, ib
    real(R8) :: xn(3)
    integer  :: idx(3)
    idx(d) = side
    idx(o1) = ia
    idx(o2) = ib
    xn = x(:,idx(1),idx(2),idx(3))
  end function Face_Node

  function Cross ( a, b ) result ( cr )
    real(R8), intent(in) :: a(3), b(3)
    real(R8) :: cr(3)
    cr = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
  end function Cross

end program cell_volume_check
