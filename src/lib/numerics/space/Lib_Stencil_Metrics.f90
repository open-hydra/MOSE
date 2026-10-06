!> Face metrics of the viscous fluxes with diffusive-metric = stencil.
!>
!> The face gradient of MOSE takes, for each quantity, one difference across the face (the two cells) and one per
!> tangential direction (the mean of the central differences of the two cells), D phi = (D_a, D_b, D_c), and maps them to
!> Cartesian components with a face metric: grad phi = M_f^T D phi. Here M_f comes from the same differences of the cell
!> positions: J^T = [D_a x, D_b x, D_c x] (columns), M_f = (J^T)^-1. For a linear field phi = phi0 + g . x every difference
!> is D phi = J g, so the face gradient is g on any mesh: skewed, stretched, with jumps of size, folded at the faces between
!> blocks. "In order to obtain the identity tensor numerically ... the two difference approximations must be the same"
!> (M. Vinokur, An analysis of finite-difference and finite-volume formulations of conservation laws, NASA CR-177416,
!> 1986, p. 26; J. Comput. Phys. 81, 1989). The inverse-mean metric of the two cells, [1/2 (A_L + A_R)]^-T, has
!> the right row across the face but not along it: its tangential rows are the mean edges of the two cells, not the
!> differences of the centres that D_b and D_c take, so on skewed, non-smooth meshes the face gradient of a linear field
!> is wrong.
!>
!> Rows kept from the inverse-mean metric, so that its bits do not change where nothing has to: the row across the face
!> (1/2 (a_L + a_R) is x_R - x_L with the centre as the mean of the eight nodes), and the rows of the index directions with
!> one cell (a 1-D case, the k direction of a 2-D mesh or wedge): there the field does not vary, D phi = 0.
!>
!> Positions: the centre of a cell is the mean of its eight nodes; the first layer of ghost cells takes the position at
!> which its value is right for a linear field: the centre of the cell across a face between blocks (plus the translation of
!> a periodic pair), 3 x1 - 3 x2 + x3 behind the faces whose ghost value is 3 u1 - 3 u2 + u3 (Ghost_Extrapolate: walls,
!> 400, 408, ...; 2 x1 - x2, or x1, in a block with two cells, or one, across the face), the cell itself behind the
!> zero-gradient faces (0, 401-407, 420), its mirror image across a symmetry face (300), and on a 3-D rotational periodic
!> face (200) the cell of the other end rotated by the angle of the sector. Chimera (102) and mapped (410) faces take the
!> extrapolated position too, an approximation: their values come from donor cells or from a file. The centres of every
!> block are computed on every rank, from the nodes every rank holds until the end of the setup of the metrics, so no
!> message is needed; the face metrics are stored for the blocks of the rank only (blk%dir(d)%Mf, on the faces of the
!> block, interior and boundary).
module MOSE_Lib_Stencil_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Stencil_Metrics, Viscous_Face_Metric, METRIC_MEAN, METRIC_INVERSE_MEAN, METRIC_STENCIL

  integer, parameter :: METRIC_MEAN = 0, METRIC_INVERSE_MEAN = 1, METRIC_STENCIL = 2

  !> Cell centres of one block with one layer of ghost cells: x(1:3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)
  type :: centres_t
    real(R8), allocatable :: x(:,:,:,:)
  end type centres_t

contains

  !> Face metric of the interior face of direction d between cell (i,j,k) and the next one along d, for the viscous
  !> fluxes: diffusive-metric = mean (M_L + M_R)/2, inverse-mean [(M_L^-1 + M_R^-1)/2]^-1 (as Diffusive_Flux computed them
  !> before, same operations), stencil the stored metric of the face.
  function Viscous_Face_Metric ( blk, d, i, j, k, metric ) result ( Mf )
    use MOSE_Advanced_Types_m, only: MOSE_block_type
    use MOSE_Lib_Diffusive, only: Inverse3
    implicit none
    type(MOSE_block_type), intent(in) :: blk
    integer, intent(in) :: d, i, j, k, metric
    real(R8) :: Mf(3,3)
    ! Local
    integer :: i2, j2, k2

    i2 = i ; j2 = j ; k2 = k
    select case (d)
      case (1) ; i2 = i + 1
      case (2) ; j2 = j + 1
      case (3) ; k2 = k + 1
    end select
    select case (metric)
      case (METRIC_STENCIL)
        Mf = blk % dir(d) % Mf(i,j,k) % c
      case (METRIC_INVERSE_MEAN)
        Mf = Inverse3 ( 0.5d0 * ( blk % Minv(i,j,k) % c + blk % Minv(i2,j2,k2) % c ) )
      case default
        Mf = 0.5d0 * ( blk % M(i,j,k) % c + blk % M(i2,j2,k2) % c )
    end select

  end function Viscous_Face_Metric


  !> Centre of cell (i,j,k) of a block: the mean of its eight nodes.
  pure function Cell_Centre ( node, i, j, k ) result ( x )
    use MOSE_Base_Types_m, only: MOSE_vector_3D_type
    implicit none
    type(MOSE_vector_3D_type), intent(in) :: node(0:,0:,0:)
    integer, intent(in) :: i, j, k
    real(R8) :: x(3)

    x = 0.125d0 * ( node(i-1,j-1,k-1) % c + node(i,j-1,k-1) % c + node(i-1,j,k-1) % c + node(i,j,k-1) % c &
                  + node(i-1,j-1,k  ) % c + node(i,j-1,k  ) % c + node(i-1,j,k  ) % c + node(i,j,k  ) % c )

  end function Cell_Centre


  !> Centre (mean of the four nodes) and unit normal (direction of growing index) of face Fm of cell (i,j,k).
  pure subroutine Face_Centre_Normal ( node, i, j, k, Fm, xf, nf )
    use MOSE_Base_Types_m, only: MOSE_vector_3D_type
    implicit none
    type(MOSE_vector_3D_type), intent(in) :: node(0:,0:,0:)
    integer, intent(in) :: i, j, k, Fm
    real(R8), intent(out) :: xf(3), nf(3)
    ! Local
    real(R8) :: p1(3), p2(3), p3(3), p4(3), d1(3), d2(3), v(3)
    integer :: a

    select case (Fm)
      case (1:2)
        a = i - mod(Fm,2)
        p1 = node(a,j-1,k-1) % c ; p2 = node(a,j,k-1) % c ; p3 = node(a,j,k) % c ; p4 = node(a,j-1,k) % c
      case (3:4)
        a = j - mod(Fm,2)
        p1 = node(i-1,a,k-1) % c ; p2 = node(i,a,k-1) % c ; p3 = node(i,a,k) % c ; p4 = node(i-1,a,k) % c
      case default
        a = k - mod(Fm,2)
        p1 = node(i-1,j-1,a) % c ; p2 = node(i,j-1,a) % c ; p3 = node(i,j,a) % c ; p4 = node(i-1,j,a) % c
    end select
    xf = 0.25d0 * ( p1 + p2 + p3 + p4 )
    d1 = p3 - p1
    d2 = p4 - p2
    v = [ d1(2)*d2(3) - d1(3)*d2(2), d1(3)*d2(1) - d1(1)*d2(3), d1(1)*d2(2) - d1(2)*d2(1) ]
    nf = v / norm2(v)

  end subroutine Face_Centre_Normal


  !> Mirror image of x across the plane of a face (centre xf, unit normal nf).
  pure function Mirror ( x, xf, nf ) result ( y )
    implicit none
    real(R8), intent(in) :: x(3), xf(3), nf(3)
    real(R8) :: y(3)

    y = x - 2d0 * dot_product ( x - xf, nf ) * nf

  end function Mirror


  !> Positions of the cell centres (every block) and face metrics (blocks of this rank) for diffusive-metric = stencil.
  !> Called by Setup_Metrics for every multigrid level, after the ghost metrics.
  subroutine Setup_Stencil_Metrics ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: gc
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    type(centres_t), allocatable :: xc(:)
    integer :: b, i, j, k, n(3), d, r

    allocate ( xc(domain % nb) )

    ! 1. Centres of the cells of every block
    do b = 1, domain % nb
      n = domain % blk(b) % dim
      allocate ( xc(b) % x(3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc) )
      xc(b) % x = 0d0
      do k = 1, n(3) ; do j = 1, n(2) ; do i = 1, n(1)
        xc(b) % x(:,i,j,k) = Cell_Centre ( domain % blk(b) % node, i, j, k )
      enddo ; enddo ; enddo
    enddo

    ! 2. Ghost cells of the faces that are not between blocks, every block
    do r = 1, domain % nbound
      select case ( domain % bc(r) % type )
        case (101, 201)
          cycle
        case default
          call Ghost_Position ( domain, xc, r )
      end select
    enddo

    ! 3. Ghost cells of the faces between blocks: the centre of the cell across (computed in 1, not a ghost)
    do r = 1, domain % nbound
      select case ( domain % bc(r) % type )
        case (101, 201)
          call Connection_Ghost_Position ( domain, xc, r )
      end select
    enddo

    ! 4. Face metrics of the blocks of this rank
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n = domain % blk(b) % dim
      if ( .not. allocated(domain % blk(b) % dir(1) % Mf) ) then
        allocate ( domain % blk(b) % dir(1) % Mf(0:n(1), n(2), n(3)) )
        allocate ( domain % blk(b) % dir(2) % Mf(n(1), 0:n(2), n(3)) )
        allocate ( domain % blk(b) % dir(3) % Mf(n(1), n(2), 0:n(3)) )
      end if
      do d = 1, 3
        domain % blk(b) % dir(d) % Mf(:,:,:) = MOSE_tensor_3D_type(0d0)
      enddo
      call Interior_Face_Metrics ( b, domain % blk(b), xc(b) )
    enddo
    do r = 1, domain % nbound
      b = domain % bc(r) % b
      if (.not. is_local_block(b)) cycle
      select case ( domain % bc(r) % type )
        case (101, 201)
          call Connection_Face_Metric ( domain, xc, r )
        case (300)
          call Symmetry_Face_Metric ( domain % blk(b), xc(b), domain % bc(r) )
      end select
    enddo

    deallocate ( xc )

  end subroutine Setup_Stencil_Metrics


  !> Position of the first ghost cell of boundary record r (not a face between blocks).
  subroutine Ghost_Position ( domain, xc, r )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: ndir
    use MOSE_Parameters_m, only: guide
    use MOSE_Lib_Metrics, only: delthe
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    type(centres_t), intent(inout) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, Im, Jm, Km, Fm, Ig, Jg, Kg, d, nd, Ks
    integer :: c1(3), c2(3), c3(3)
    real(R8) :: x1(3), x2(3), x3(3), xf(3), nf(3), y(3), ang

    b  = domain % bc(r) % b
    Im = domain % bc(r) % i ; Jm = domain % bc(r) % j ; Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    Ig = Im - guide(Fm,1) ; Jg = Jm - guide(Fm,2) ; Kg = Km - guide(Fm,3)
    d  = (Fm + 1) / 2
    nd = domain % blk(b) % dim(d)
    c1 = [Im, Jm, Km]
    c2 = c1 + guide(Fm,:)
    c3 = c2 + guide(Fm,:)
    x1 = xc(b) % x(:,c1(1),c1(2),c1(3))

    select case ( domain % bc(r) % type )

      case (0, 401:407, 420)                    ! zero gradient: the ghost value is the cell value
        xc(b) % x(:,Ig,Jg,Kg) = x1

      case (300)                                ! symmetry: mirror image across the face
        call Face_Centre_Normal ( domain % blk(b) % node, Im, Jm, Km, Fm, xf, nf )
        xc(b) % x(:,Ig,Jg,Kg) = Mirror ( x1, xf, nf )

      case (200)                                ! wedge (one cell in k: unused) or 3-D rotational periodicity
        if ( ndir == 3 .and. domain % blk(b) % dim(3) > 1 ) then
          Ks = domain % blk(b) % dim(3) - Km + 1
          y = xc(b) % x(:,Im,Jm,Ks)
          ang = delthe
          if ( Fm == 5 ) ang = -delthe
          xc(b) % x(:,Ig,Jg,Kg) = [ y(1), y(2)*cos(ang) - y(3)*sin(ang), y(2)*sin(ang) + y(3)*cos(ang) ]
        else
          xc(b) % x(:,Ig,Jg,Kg) = x1
        end if

      case default                              ! 3 u1 - 3 u2 + u3 (Ghost_Extrapolate), chimera, mapped
        if ( nd >= 3 ) then
          x2 = xc(b) % x(:,c2(1),c2(2),c2(3))
          x3 = xc(b) % x(:,c3(1),c3(2),c3(3))
          xc(b) % x(:,Ig,Jg,Kg) = 3d0 * x1 - 3d0 * x2 + x3
        elseif ( nd == 2 ) then
          x2 = xc(b) % x(:,c2(1),c2(2),c2(3))
          xc(b) % x(:,Ig,Jg,Kg) = 2d0 * x1 - x2
        else
          xc(b) % x(:,Ig,Jg,Kg) = x1
        end if

    end select

  end subroutine Ghost_Position


  !> Position of the first ghost cell of a face between blocks (101) or of a periodic pair (201): the centre of the source
  !> cell, plus the translation from the source face to this face for 201.
  subroutine Connection_Ghost_Position ( domain, xc, r )
    use MOSE_Advanced_Types_m
    use MOSE_Parameters_m, only: guide
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    type(centres_t), intent(inout) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, bs, Ig, Jg, Kg

    b  = domain % bc(r) % b
    bs = domain % bc(r) % bs
    Ig = domain % bc(r) % i - guide(domain % bc(r) % f,1)
    Jg = domain % bc(r) % j - guide(domain % bc(r) % f,2)
    Kg = domain % bc(r) % k - guide(domain % bc(r) % f,3)
    xc(b) % x(:,Ig,Jg,Kg) = xc(bs) % x(:,domain % bc(r) % is, domain % bc(r) % js, domain % bc(r) % ks) &
                          + Periodic_Shift ( domain, r )

  end subroutine Connection_Ghost_Position


  !> Translation of a periodic pair (201): centre of this face minus centre of the source face; zero for 101.
  function Periodic_Shift ( domain, r ) result ( t )
    use MOSE_Advanced_Types_m
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in) :: r
    real(R8) :: t(3)
    ! Local
    real(R8) :: xfm(3), xfs(3), nf(3)
    integer :: b, bs

    t = 0d0
    if ( domain % bc(r) % type /= 201 ) return
    b  = domain % bc(r) % b
    bs = domain % bc(r) % bs
    call Face_Centre_Normal ( domain % blk(b) % node, domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k, &
                              domain % bc(r) % f, xfm, nf )
    call Face_Centre_Normal ( domain % blk(bs) % node, domain % bc(r) % is, domain % bc(r) % js, domain % bc(r) % ks, &
                              domain % bc(r) % fs, xfs, nf )
    t = xfm - xfs

  end function Periodic_Shift


  !> Transposed face Jacobian of the inverse-mean metric of two cells, (A_L^T + A_R^T)/2 with A^T = M^-1: its columns
  !> are the mean edges; the column across the face is x_R - x_L.
  function Mean_Edges ( M1, M2 ) result ( JT )
    use MOSE_Lib_Diffusive, only: Inverse3
    implicit none
    real(R8), intent(in) :: M1(3,3), M2(3,3)
    real(R8) :: JT(3,3)

    JT = 0.5d0 * ( Inverse3 ( M1 ) + Inverse3 ( M2 ) )

  end function Mean_Edges


  !> Face metric M_f = (J^T)^-1 of the face of direction d and indices (i,j,k) of block b. The run stops when the three
  !> differences of the stencil are not independent, |det J^T| <= 1e-12 |D_a x| |D_b x| |D_c x|: a degenerate cell, or a
  !> face whose stencil collapses (a polar line of a 3-D mesh), where this metric is not defined.
  function Face_Inverse ( JT, b, d, i, j, k, kind ) result ( Mf )
    use MOSE_Lib_Diffusive, only: Inverse3
    use MOSE_Mod_MPI, only: mpi_abort_all
    implicit none
    real(R8), intent(in) :: JT(3,3)
    integer, intent(in) :: b, d, i, j, k
    character(len=*), intent(in) :: kind
    real(R8) :: Mf(3,3)
    ! Local
    real(R8) :: det
    integer :: t

    det = JT(1,1) * ( JT(2,2)*JT(3,3) - JT(2,3)*JT(3,2) ) &
        - JT(1,2) * ( JT(2,1)*JT(3,3) - JT(2,3)*JT(3,1) ) &
        + JT(1,3) * ( JT(2,1)*JT(3,2) - JT(2,2)*JT(3,1) )
    if ( .not. ( abs(det) > 1d-12 * norm2(JT(:,1)) * norm2(JT(:,2)) * norm2(JT(:,3)) ) ) then
      write(*,'(A,I0,A,I0,A,3(1X,I0),A)') ' [ERROR] diffusive-metric = stencil: block ', b, ', face of direction ', &
                                          d, ' at', i, j, k, ' ('//kind//')'
      write(*,'(A)') ' [ERROR]   the differences of the cell centres of its stencil are not independent:'
      do t = 1, 3
        write(*,'(A,I0,A,3(1X,ES14.6))') ' [ERROR]   difference ', t, ':', JT(:,t)
      enddo
      write(*,'(A)') ' [ERROR] A degenerate cell, or a stencil that collapses (polar line of a 3-D mesh).'
      write(*,'(A)') ' [ERROR] diffusive-metric = inverse-mean does not need this metric.'
      call mpi_abort_all ( 'face metric of the viscous fluxes not defined (see above)' )
    end if
    Mf = Inverse3 ( JT )

  end function Face_Inverse


  !> Face metrics of the interior faces of block b.
  subroutine Interior_Face_Metrics ( b, blk, xc )
    use MOSE_Advanced_Types_m, only: MOSE_block_type
    implicit none
    integer, intent(in) :: b
    type(MOSE_block_type), intent(inout) :: blk
    type(centres_t), intent(in) :: xc
    ! Local
    integer :: i, j, k, n(3), d, t, i2, j2, k2
    integer :: e(3,3)
    real(R8) :: JT(3,3)

    n = blk % dim
    e = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    do d = 1, 3
      do k = 1, n(3) - e(3,d) ; do j = 1, n(2) - e(2,d) ; do i = 1, n(1) - e(1,d)
        i2 = i + e(1,d) ; j2 = j + e(2,d) ; k2 = k + e(3,d)
        JT = Mean_Edges ( blk % M(i,j,k) % c, blk % M(i2,j2,k2) % c )
        do t = 1, 3
          if ( t == d .or. n(t) == 1 ) cycle
          JT(:,t) = 0.25d0 * ( ( xc % x(:,i +e(1,t),j +e(2,t),k +e(3,t)) - xc % x(:,i -e(1,t),j -e(2,t),k -e(3,t)) ) &
                             + ( xc % x(:,i2+e(1,t),j2+e(2,t),k2+e(3,t)) - xc % x(:,i2-e(1,t),j2-e(2,t),k2-e(3,t)) ) )
        enddo
        blk % dir(d) % Mf(i,j,k) % c = Face_Inverse ( JT, b, d, i, j, k, 'interior' )
      enddo ; enddo ; enddo
    enddo

  end subroutine Interior_Face_Metrics


  !> Face metric of the face of boundary record r between blocks, stored on the face of block b. The two cells in the
  !> order of growing index of block b (L, R) and their tangential neighbours: in block b, and for the ghost cell the
  !> cells of the source block that Fill_BC_Ghost_Connection takes for Pg(:,3:6).
  subroutine Connection_Face_Metric ( domain, xc, r )
    use MOSE_Advanced_Types_m
    use MOSE_Parameters_m, only: guide
    use MOSE_Lib_Diffusive, only: Inverse3
    use MOSE_Lib_BC_Fluxes, only: Face_Index
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(centres_t), intent(in) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, bs, Im, Jm, Km, Fm, Ig, Jg, Kg, dir, fi, fj, fk, t, tb(2), q, n(3)
    integer :: e(3,3), sc(3,4)
    real(R8) :: JT(3,3), xg(3,4), shift(3), dloc(3), dgh(3)

    b  = domain % bc(r) % b
    bs = domain % bc(r) % bs
    Im = domain % bc(r) % i ; Jm = domain % bc(r) % j ; Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    Ig = Im - guide(Fm,1) ; Jg = Jm - guide(Fm,2) ; Kg = Km - guide(Fm,3)
    n  = domain % blk(b) % dim
    e  = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    call Face_Index ( Fm, dir, Im, Jm, Km, fi, fj, fk )
    select case (dir)
      case (1) ; tb = [2, 3]
      case (2) ; tb = [1, 3]
      case default ; tb = [1, 2]
    end select

    ! Cells of the source block that hold Pg(:,3:6) (Fill_BC_Ghost_Connection), and their positions
    call Source_Neighbours ( domain % bc(r), sc )
    shift = Periodic_Shift ( domain, r )
    do q = 1, 4
      xg(:,q) = xc(bs) % x(:,sc(1,q),sc(2,q),sc(3,q)) + shift
    enddo

    ! Face Jacobian: the cells in the order L, R of block b
    if ( mod(Fm,2) == 0 ) then
      JT = Mean_Edges ( domain % blk(b) % M(Im,Jm,Km) % c, domain % blk(b) % M(Ig,Jg,Kg) % c )
    else
      JT = Mean_Edges ( domain % blk(b) % M(Ig,Jg,Kg) % c, domain % blk(b) % M(Im,Jm,Km) % c )
    end if
    do q = 1, 2
      t = tb(q)
      if ( n(t) == 1 ) cycle
      dloc = xc(b) % x(:,Im+e(1,t),Jm+e(2,t),Km+e(3,t)) - xc(b) % x(:,Im-e(1,t),Jm-e(2,t),Km-e(3,t))
      dgh  = xg(:,2*q) - xg(:,2*q-1)
      if ( mod(Fm,2) == 0 ) then
        JT(:,t) = 0.25d0 * ( dloc + dgh )
      else
        JT(:,t) = 0.25d0 * ( dgh + dloc )
      end if
    enddo
    domain % blk(b) % dir(dir) % Mf(fi,fj,fk) % c = Face_Inverse ( JT, b, dir, fi, fj, fk, 'between blocks' )

  end subroutine Connection_Face_Metric


  !> Indices in the source block of the cells of Pg(:,3:6) of a connection record, as in Fill_BC_Ghost_Connection.
  subroutine Source_Neighbours ( rec, sc )
    use MOSE_Advanced_Types_m, only: MOSE_bc_type
    implicit none
    type(MOSE_bc_type), intent(in) :: rec
    integer, intent(out) :: sc(3,4)
    ! Local
    integer :: Is, Js, Ks

    Is = rec % is ; Js = rec % js ; Ks = rec % ks
    select case ( rec % fs )
      case (1:2)
        sc(:,1) = [Is, Js - rec % d11, Ks - rec % d12]
        sc(:,2) = [Is, Js + rec % d11, Ks + rec % d12]
        sc(:,3) = [Is, Js - rec % d21, Ks - rec % d22]
        sc(:,4) = [Is, Js + rec % d21, Ks + rec % d22]
      case (3:4)
        sc(:,1) = [Is - rec % d11, Js, Ks - rec % d12]
        sc(:,2) = [Is + rec % d11, Js, Ks + rec % d12]
        sc(:,3) = [Is - rec % d21, Js, Ks - rec % d22]
        sc(:,4) = [Is + rec % d21, Js, Ks + rec % d22]
      case default
        sc(:,1) = [Is - rec % d11, Js - rec % d12, Ks]
        sc(:,2) = [Is + rec % d11, Js + rec % d12, Ks]
        sc(:,3) = [Is - rec % d21, Js - rec % d22, Ks]
        sc(:,4) = [Is + rec % d21, Js + rec % d22, Ks]
    end select

  end subroutine Source_Neighbours


  !> Face metric of a symmetry face (300) of block b: the cell and its mirror image across the face, D_a = (mirror - cell)
  !> along the direction of growing index, and the tangential differences of the cell and of the mirror images of its
  !> neighbours (BC_Symmetry_Visc). Faces across an index direction with one cell (the k faces of a 2-D mesh) keep the cell
  !> metric: nothing to store.
  subroutine Symmetry_Face_Metric ( blk, xc, rec )
    use MOSE_Advanced_Types_m
    use MOSE_Lib_Diffusive, only: Inverse3
    use MOSE_Lib_BC_Fluxes, only: Face_Index
    implicit none
    type(MOSE_block_type), intent(inout) :: blk
    type(centres_t), intent(in) :: xc
    type(MOSE_bc_type), intent(in) :: rec
    ! Local
    integer :: Im, Jm, Km, Fm, dir, fi, fj, fk, t, n(3), modfm2
    integer :: e(3,3)
    real(R8) :: JT(3,3), xf(3), nf(3), x1(3), xp(3), xm(3)

    Im = rec % i ; Jm = rec % j ; Km = rec % k ; Fm = rec % f
    n = blk % dim
    e = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    call Face_Index ( Fm, dir, Im, Jm, Km, fi, fj, fk )
    if ( n(dir) == 1 ) return
    modfm2 = 1 - 2 * mod(Fm,2)
    call Face_Centre_Normal ( blk % node, Im, Jm, Km, Fm, xf, nf )
    x1 = xc % x(:,Im,Jm,Km)
    JT = Inverse3 ( blk % M(Im,Jm,Km) % c )
    JT(:,dir) = ( Mirror ( x1, xf, nf ) - x1 ) * modfm2
    do t = 1, 3
      if ( t == dir .or. n(t) == 1 ) cycle
      xp = xc % x(:,Im+e(1,t),Jm+e(2,t),Km+e(3,t))
      xm = xc % x(:,Im-e(1,t),Jm-e(2,t),Km-e(3,t))
      JT(:,t) = 0.25d0 * ( ( xp - xm ) + ( Mirror ( xp, xf, nf ) - Mirror ( xm, xf, nf ) ) )
    enddo
    blk % dir(dir) % Mf(fi,fj,fk) % c = Face_Inverse ( JT, rec % b, dir, fi, fj, fk, 'symmetry' )

  end subroutine Symmetry_Face_Metric

end module MOSE_Lib_Stencil_Metrics
