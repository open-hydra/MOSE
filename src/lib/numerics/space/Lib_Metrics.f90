module MOSE_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Base_Types_m
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m

  implicit none
  real(R8), public :: delthe = 0d0  ! grid axisymmetric angle

  ! Wall faces of Compute_Yn, read by Wall_Face_Distance while the k-d tree is queried: for face n, wall_corner(:,1:4,n)
  ! are its corners in cyclic order and wall_corner(:,0,n) their mean; on an axisymmetric wedge (wall_meridional) the
  ! corners are in the meridional plane (x, r, 0), where corners 1 and 2 are the ends of the trace of the face.
  real(R8), allocatable :: wall_corner(:,:,:)
  logical               :: wall_meridional = .false.
  private :: wall_corner, wall_meridional
  private :: Cell_Centre, Meridional, Wall_Face_Distance, Point_Segment_Distance, Point_Triangle_Distance, Cross3

contains

  !> Mesh type, from the nodes of every block.
  !>  - ndir is the number of index directions with more than one cell: 3 if a block has more than one cell in k,
  !>    otherwise 2 if a block has more than one cell in j, otherwise 1.
  !>  - delthe is the angle of the axisymmetric wedge of a 2-D mesh, or of the sectors closed by the rotational
  !>    periodicity (type 200 on faces 5 and 6) of a 3-D mesh, and 0 otherwise. A 2-D mesh is a wedge when, in every
  !>    block, the k = 1 node plane is the k = 0 plane turned about the x axis by the same angle; it is plane when no
  !>    block is so turned. In 3-D the k = dim(3) node plane of every block with a type-200 face is the k = 0 plane turned
  !>    about the x axis by the same angle, or translated (delthe = 0). Every node is checked (Node_Planes); blocks that
  !>    disagree stop the run.
  !>  - The angle is the azimuth turned by the node (0,1) of the block, as it always was, or by the node farthest from
  !>    the axis when that one lies on it. The angles of two blocks agree when they differ by at most ETA of the angle.
  !>  - The checks are relative to the move of the nodes per cell in k, the largest move of a node from one plane to the
  !>    other over the cells between them (the thickness of a 2-D block), not to the size of the coordinates: the
  !>    coordinates of a mesh are written with a finite number of digits, and their round-off, relative to that move, is
  !>    what a turn or a translation leaves over. A check relative to the size of the coordinates, at 1e-6 of it, fails a
  !>    wedge of 5 degrees written with 5 significant digits, which is then taken for a plane mesh, and passes a plane slab
  !>    far from the x axis, which is then taken for a wedge. A 2-D block nearer a turn than a translation, or nearer a
  !>    translation than a turn, but not to within ETA of its thickness, stops the run: its nodes are written with too
  !>    few digits to tell its type, or its k planes are neither parallel nor turned about the x axis.
  subroutine Check_Mesh_Type( domain )
    use MOSE_Mod_MPI, only: mpi_abort_all
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    real(R8), parameter :: TOL = 1d-6, ETA = 1d-2, PI = 3.141592653589793d0
    integer  :: b, l, nb, nturn, nmove
    real(R8) :: sr, ang, ang_b, fit, fit_move
    logical  :: turned, moved, nearer_turn, nearer_move
    logical, allocatable :: sector(:)
    character(len=512) :: msg   ! room for the longest message: 255 characters and the block number

    nb = domain % nb

    ndir = 1
    do b = 1, nb
      if ( domain % Blk(b) % dim(2) > 1 ) ndir = max ( ndir, 2 )
      if ( domain % Blk(b) % dim(3) > 1 ) ndir = 3
    enddo

    ! Size of the distance from the x axis, to tell a probe node on the axis
    sr = 0d0
    do b = 1, nb
      associate ( node => domain % Blk(b) % node )
        sr = max ( sr, maxval ( sqrt ( node(:,:,:) % c(2)**2 + node(:,:,:) % c(3)**2 ) ) )
      end associate
    enddo

    delthe = 0d0
    select case ( ndir )

      case ( 2 )
        ! Wedge if every block is turned by the same angle; plane if none is
        nturn = 0
        do b = 1, nb
          call Node_Planes ( domain % Blk(b), 1, sr, TOL, ETA, turned, ang_b, moved, fit, nearer_turn, fit_move, &
                             nearer_move )
          if ( nearer_turn .and. .not. turned ) then
            write(msg,'(A,I0,A,ES10.3,A,ES9.2,A)') '[ERROR] Mesh type: block ', b, ' is nearer a wedge of angle ', ang_b, &
              ' rad than a plane block, but its nodes match the turn only to ', fit, ' of its thickness (1e-2 needed):'// &
              ' write the node coordinates with more significant digits'
            call mpi_abort_all ( trim(msg) )
          end if
          if ( nearer_move .and. .not. moved ) then
            write(msg,'(A,I0,A,ES9.2,A)') '[ERROR] Mesh type: block ', b, ' is nearer a plane block than a wedge, but its'// &
              ' k = 1 node plane is the k = 0 plane translated only to ', fit_move, ' of its thickness (1e-2 needed):'// &
              ' its k planes must be parallel, or the coordinates written with more significant digits'
            call mpi_abort_all ( trim(msg) )
          end if
          if ( .not. turned ) cycle
          nturn = nturn + 1
          if ( nturn == 1 ) then
            ang = ang_b
          elseif ( abs ( ang_b - ang ) > ETA * abs ( ang ) ) then
            write(msg,'(A,I0,A,ES12.5,A,ES12.5,A)') '[ERROR] Mesh type: block ', b, ' is a wedge of angle ', ang_b, &
              ' rad, another block one of ', ang, ' rad'
            call mpi_abort_all ( trim(msg) )
          end if
        enddo
        if ( nturn == nb ) then
          delthe = ang
        elseif ( nturn > 0 ) then
          write(msg,'(A,I0,A,I0,A)') '[ERROR] Mesh type: ', nturn, ' of the ', nb, &
            ' blocks of the 2-D mesh are axisymmetric wedges (k = 1 plane turned about the x axis), the others are not'
          call mpi_abort_all ( trim(msg) )
        end if

      case ( 3 )
        ! Sectors closed by the rotational periodicity: the blocks with a type-200 face on faces 5-6
        allocate ( sector(nb) )
        sector = .false.
        do l = 1, domain % nbound
          if ( domain % bc(l) % type == 200 .and. domain % bc(l) % f >= 5 ) sector(domain % bc(l) % b) = .true.
        enddo
        nturn = 0
        nmove = 0
        do b = 1, nb
          if ( .not. sector(b) ) cycle
          call Node_Planes ( domain % Blk(b), domain % Blk(b) % dim(3), sr, TOL, ETA, turned, ang_b, moved, fit, &
                             nearer_turn, fit_move, nearer_move )
          if ( turned ) then
            nturn = nturn + 1
            if ( nturn == 1 ) then
              ang = ang_b
            elseif ( abs ( ang_b - ang ) > ETA * abs ( ang ) ) then
              write(msg,'(A,I0,A,ES12.5,A,ES12.5,A)') '[ERROR] Mesh type: block ', b, ' is a sector of angle ', ang_b, &
                ' rad, another block one of ', ang, ' rad'
              call mpi_abort_all ( trim(msg) )
            end if
          elseif ( moved ) then
            nmove = nmove + 1
          else
            write(msg,'(A,I0,A,ES9.2,A,ES9.2,A)') '[ERROR] Mesh type: block ', b, ' has type 200 on its k faces, but its'// &
              ' last k node plane is neither the first one turned about the x axis nor translated, to 1e-2 of the move of'// &
              ' its nodes per cell in k (the turn fits to ', fit, ', the translation to ', fit_move, ')'
            call mpi_abort_all ( trim(msg) )
          end if
        enddo
        if ( nturn > 0 .and. nmove > 0 ) then
          call mpi_abort_all ( '[ERROR] Mesh type: type 200 closes sectors turned about the x axis in some blocks and'// &
                               ' translated ones in others' )
        end if
        if ( nturn > 0 ) delthe = ang
        deallocate ( sector )

    end select

  contains

    !> Compare the node planes k = 0 and k = ktop of a block. The thickness of the block is the largest move of a node
    !> from the plane 0 to the plane ktop (its largest component), and cell = thickness / ktop the move per cell in k (in
    !> 2-D, ktop = 1, the thickness). turned: every node of the plane ktop is the node of the plane 0 turned about the x
    !> axis by ang, with |ang| > tol, to within eta of cell, and closer than as a translation. moved: every node is the
    !> node of the plane 0 translated by the move of the node (0,0), to within eta of cell, and closer than as a turn.
    !> Comparing the two keeps a wedge from passing for a translation and a slab far from the axis, which a small turn
    !> nearly fits, from passing for a wedge; taking cell, not the thickness, keeps a periodic plane of a sector of many
    !> cells that does not match by a fraction of a cell from passing for a turn. fit, fit_move: the misfits of the turn
    !> and of the translation over cell; nearer_turn: |ang| > tol and the turn fits better than the translation;
    !> nearer_move: the translation fits better than the turn.
    subroutine Node_Planes ( blk, ktop, sr, tol, eta, turned, ang, moved, fit, nearer_turn, fit_move, nearer_move )
      type(MOSE_block_type), intent(in) :: blk
      integer, intent(in)               :: ktop
      real(R8), intent(in)              :: sr, tol, eta
      logical, intent(out)              :: turned, moved, nearer_turn, nearer_move
      real(R8), intent(out)             :: ang, fit, fit_move
      ! Local
      integer  :: i, j, i0, j0
      real(R8) :: theta1, theta2, ca, sa, rmax, r, a(3), c(3), shift(3), turn_off, move_off, thick, cell

      associate ( node => blk % node )
        ! Probe node: (0,1), or the node of the plane 0 farthest from the axis when (0,1) lies on it
        i0 = 0
        j0 = 1
        if ( sqrt ( node(0,1,0)%c(2)**2 + node(0,1,0)%c(3)**2 ) <= tol * sr ) then
          rmax = -1d0
          do j = 0, blk % dim(2)
            do i = 0, blk % dim(1)
              r = sqrt ( node(i,j,0)%c(2)**2 + node(i,j,0)%c(3)**2 )
              if ( r > rmax ) then
                rmax = r
                i0 = i
                j0 = j
              end if
            enddo
          enddo
        end if
        theta1 = atan2 ( node(i0,j0,0)%c(3), node(i0,j0,0)%c(2) )
        theta2 = atan2 ( node(i0,j0,ktop)%c(3), node(i0,j0,ktop)%c(2) )
        ang = theta2 - theta1
        if ( ang >  pi ) ang = ang - 2d0 * pi
        if ( ang < -pi ) ang = ang + 2d0 * pi

        ca = cos ( ang )
        sa = sin ( ang )
        shift = node(0,0,ktop)%c - node(0,0,0)%c
        turn_off = 0d0
        move_off = 0d0
        thick    = 0d0
        do j = 0, blk % dim(2)
          do i = 0, blk % dim(1)
            a = node(i,j,0)%c
            c = node(i,j,ktop)%c
            turn_off = max ( turn_off, abs ( c(1) - a(1) ), abs ( c(2) - ( a(2) * ca - a(3) * sa ) ), &
                                       abs ( c(3) - ( a(2) * sa + a(3) * ca ) ) )
            move_off = max ( move_off, maxval ( abs ( c - a - shift ) ) )
            thick    = max ( thick, maxval ( abs ( c - a ) ) )
          enddo
        enddo
        cell = thick / ktop
        nearer_turn = abs ( ang ) > tol .and. turn_off < move_off
        nearer_move = move_off < turn_off
        turned = nearer_turn .and. turn_off <= eta * cell
        moved  = nearer_move .and. move_off <= eta * cell
        fit      = huge ( 1d0 )
        fit_move = huge ( 1d0 )
        if ( cell > 0d0 ) fit = turn_off / cell
        if ( cell > 0d0 ) fit_move = move_off / cell
      end associate

    end subroutine Node_Planes

  end subroutine Check_Mesh_Type


  subroutine Compute_Norm_Area( b )
    implicit none
    class(MOSE_block_type) :: b  
    ! Local
    real(R8) :: Ai, snix, sniy, sniz, Aj, snjx, snjy, snjz, Ak, snkx, snky, snkz
    integer      :: i, j, k, im, jm, km
    real(R8) :: d1(3), d2(3), d3(3)
    real(R8) :: snixx, sniyy, snizz, snjxx, snjyy, snjzz, snkxx, snkyy, snkzz
    real(R8) :: scal, signi, signj, signk
  
    ! Peliminary operations
    im = b%dim(1)
    jm = b%dim(2)
    km = b%dim(3)
  
    ! compute sign of normal vectors to the intefaces
  
    ! i direction
    d1 = b%node(1,1,1)%c - b%node(1,0,0)%c
    d2 = b%node(1,0,1)%c - b%node(1,1,0)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snixx=d3(1)
    sniyy=d3(2)
    snizz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(1,0,1)%c + b%node(1,0,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(0,1,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snixx+d3(2)*sniyy+d3(3)*snizz
  
    signi=sign(1.d0,scal)

    ! j direction
    d1 = b%node(0,1,1)%c - b%node(1,1,0)%c
    d2 = b%node(1,1,1)%c - b%node(0,1,0)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snjxx=d3(1)
    snjyy=d3(2)
    snjzz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,1,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(1,0,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snjxx+d3(2)*snjyy+d3(3)*snjzz
  
    signj=sign(1.d0,scal)

    ! k direction
    d1 = b%node(0,1,1)%c - b%node(1,0,1)%c
    d2 = b%node(0,0,1)%c - b%node(1,1,1)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snkxx=d3(1)
    snkyy=d3(2)
    snkzz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(1,0,1)%c )
    d2 = 0.25d0*( b%node(1,1,0)%c + b%node(0,1,0)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snkxx+d3(2)*snkyy+d3(3)*snkzz
  
    signk=sign(1.d0,scal)
    
    ! compute metrics: n, A
    
    !$omp parallel private (d1,d2,d3,i,j,k,snix,sniy,sniz,Ai,Aj,snjx,snjy,snjz,Ak,snkx,snky,snkz)
    
    ! i direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 1, jm ; do i = 0, im
  
      d1 = b%node(i,j,k)%c - b%node(i,j-1,k-1)%c
      d2 = b%node(i,j-1,k)%c - b%node(i,j,k-1)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Ai = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      if (Ai == 0d0) then
        snix = 0d0
        sniy = 0d0
        sniz = 0d0
      else
        snix = d3(1)/Ai*signi
        sniy = d3(2)/Ai*signi
        sniz = d3(3)/Ai*signi
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(1)%f(i,j,k)%a = Ai
      b%dir(1)%f(i,j,k)%n = [ snix, sniy, sniz ]
      
    enddo; enddo; enddo
  
    ! j direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 0, jm ; do  i = 1, im
  
      d1 = b%node(i-1,j,k)%c - b%node(i,j,k-1)%c
      d2 = b%node(i,j,k)%c - b%node(i-1,j,k-1)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Aj = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      if (Aj == 0d0) then
        snjx = 0d0
        snjy = 0d0
        snjz = 0d0
      else
        snjx = d3(1)/Aj*signj
        snjy = d3(2)/Aj*signj
        snjz = d3(3)/Aj*signj
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(2)%f(i,j,k)%A = Aj
      b%dir(2)%f(i,j,k)%n = [ snjx, snjy, snjz ]
  
    enddo; enddo; enddo
  
    ! k direction
    !$omp do collapse(3)
    do k = 0, km ; do j = 1, jm ; do i = 1, im
  
      d1 = b%node(i-1,j,k)%c - b%node(i,j-1,k)%c
      d2 = b%node(i-1,j-1,k)%c - b%node(i,j,k)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Ak = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      if (Ak == 0d0) then
        snkx = 0d0
        snky = 0d0
        snkz = 0d0
      else
        snkx = d3(1)/Ak*signk
        snky = d3(2)/Ak*signk
        snkz = d3(3)/Ak*signk
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(3)%f(i,j,k)%A = Ak
      b%dir(3)%f(i,j,k)%n = [ snkx, snky, snkz ]
  
    enddo; enddo; enddo
    !$omp end parallel
  
  end subroutine Compute_Norm_Area


  subroutine Compute_Metric_Tensor( N1, N2, N3, N4, N5, N6, N7, N8, M, dl, vol )
    implicit none
    type(MOSE_vector_3D_type), intent(in)  :: N1, N2, N3, N4, N5, N6, N7, N8
    type(MOSE_tensor_3D_type), intent(out) :: M
    type(MOSE_vector_3D_type), intent(out) :: dl
    real(R8), intent(out)             :: vol
    ! Local
    integer  :: h
    real(R8) :: det, A(3,3), cofactor(3,3)
    real(R8) :: vx(8), vy(8), vz(8)

    ! Compute average cell dimension. A is M^-1, inverse of metric tensor.
    ! M1 = [ xcs, ycs, zcs ; xet, yet, zet ; xzi, yzi, zzi ].
    
    A(1,:) = N5 % c - N1 % c + &
             N6 % c - N2 % c + &
             N7 % c - N3 % c + &
             N8 % c - N4 % c
  
    A(2,:) = N3 % c - N1 % c + &
             N4 % c - N2 % c + &
             N7 % c - N5 % c + &
             N8 % c - N6 % c
             
    A(3,:) = N2 % c - N1 % c + &
             N4 % c - N3 % c + &
             N6 % c - N5 % c + &
             N8 % c - N7 % c

    A = 0.25d0 * A

    ! Determinant of M^-1.
    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    if ( abs(det) == 0d0 ) then
      ! Degenerate cell (zero volume): a mesh error, or ghost cells extrapolated from a
      ! wall-thin first cell whose coordinates were read back in single precision.
      ! stop 1, not mpi_abort_all: this can run inside OpenMP threads (MPI is FUNNELED).
      write(*,'(A)') ' [ERROR] Metric tensor det=0: degenerate cell. Check the mesh, or the precision of the restart file'
      write(*,'(A,3(1X,ES14.6))') ' [ERROR]   corner node N1:', N1 % c
      write(*,'(A,3(1X,ES14.6))') ' [ERROR]   corner node N8:', N8 % c
      stop 1

    else

      cofactor(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))
      cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      cofactor(1,3) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))
      cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      cofactor(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))
      cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      cofactor(3,1) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))
      cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      cofactor(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))

      M % c = cofactor / det

    endif

    ! Cell length in i,j,k
    do h = 1, 3
      dl % c(h) = sqrt( A(h,1)**2 + A(h,2)**2 + A(h,3)**2 )
    enddo

    ! Cell volume computation
    vx(1) = N1%c(1)
    vy(1) = N1%c(2)
    vz(1) = N1%c(3)

    vx(2) = N5%c(1)
    vy(2) = N5%c(2)
    vz(2) = N5%c(3)

    vx(3) = N3%c(1)
    vy(3) = N3%c(2)
    vz(3) = N3%c(3)

    vx(4) = N7%c(1)
    vy(4) = N7%c(2)
    vz(4) = N7%c(3)
  
    vx(5) = N2%c(1)
    vy(5) = N2%c(2)
    vz(5) = N2%c(3)

    vx(6) = N6%c(1)
    vy(6) = N6%c(2)
    vz(6) = N6%c(3)

    vx(7) = N4%c(1)
    vy(7) = N4%c(2)
    vz(7) = N4%c(3)
  
    vx(8) = N8%c(1)
    vy(8) = N8%c(2)
    vz(8) = N8%c(3)
  
    vol = tvol(vx,vy,vz,1,2,3,5) + tvol(vx,vy,vz,2,4,3,8) &
        + tvol(vx,vy,vz,5,8,6,2) + tvol(vx,vy,vz,5,7,8,3) &
        + tvol(vx,vy,vz,5,8,2,3)
  
    if( vol <= 0d0 ) then
      write(*,'(A90)') ' [ERROR] Negative volume'
      stop 1
    endif

    contains
  
      pure function tvol(vx, vy, vz, i1, i2, i3, i4) result(volume)
        implicit none
        real(R8), intent(in) :: vx(8), vy(8), vz(8)
        integer, intent(in)  :: i1, i2, i3, i4
        real(R8)             :: volume
  
        volume = abs(((vx(i2)-vx(i1))* &
          ((vy(i3)-vy(i1))*(vz(i4)-vz(i1))-(vy(i4)-vy(i1))*(vz(i3)-vz(i1)))+ &
                                  (vy(i2)-vy(i1))* &
          ((vx(i4)-vx(i1))*(vz(i3)-vz(i1))-(vx(i3)-vx(i1))*(vz(i4)-vz(i1)))+ &
                                  (vz(i2)-vz(i1))* &
          ((vx(i3)-vx(i1))*(vy(i4)-vy(i1))-(vx(i4)-vx(i1))*(vy(i3)-vy(i1)))) &
          /6.d0)
  
      end function tvol
  
  end subroutine Compute_Metric_Tensor

  
  subroutine Compute_Yn ( blk, bc, nb, nbound )
    use MOSE_Lib_KDTree, only: kdtree_t, kdtree_build, kdtree_free
    use MOSE_Mod_MPI,    only: is_local_block, mpi_is_root, mpi_abort_all
    implicit none
    integer, intent(in)                                  :: nb, nbound
    type(MOSE_block_type), dimension(nb), intent(inout)  :: blk
    type(MOSE_bc_type), dimension(nbound), intent(in)    :: bc
    ! Local
    real(R8), allocatable :: facex(:), facey(:), facez(:), facekr(:), facer(:)
    real(R8) :: corner(3,4)
    integer :: Fm, Bm, Im, Jm, Km, l, b, nwall, n
    logical, dimension(nb) :: needs_yn
    type(kdtree_t) :: wall_tree
    ! Periodic images of the wall: image n is the wall turned about the x axis by img_angle(n), then translated by
    ! img_shift(:,n); wall_box(:,1:2) are the low and high corners of the box around the wall faces
    integer               :: nimg
    real(R8), allocatable :: img_shift(:,:), img_angle(:)
    real(R8)              :: wall_box(3,2)

    ! Validation switch: checks every yn against a brute-force scan, roughly
    ! doubling setup cost.  Keep .false. in production.
    logical, parameter  :: DEBUG_KDTREE = .false.
    real(R8), parameter :: KDTREE_TOL   = 1.0e-10_R8

    ! The distance is the minimum distance from the cell centre to the wall surface
    ! (NASA Turbulence Modeling Resource, Spalart-Allmaras page: not the distance to
    ! the nearest wall node or face centre; next to a convex corner, the distance to
    ! the corner). On a 2-D axisymmetric wedge it is measured in the meridional plane
    ! (x, r), from the wall curve: that is the distance to the wall surface of
    ! revolution, whatever the angle of the wedge. Same wedge test as the wedge
    ! stress (Lib_BC_Fluxes_Rotational), delthe of either sign.
    wall_meridional = ndir == 2 .and. .not. ( delthe >= 0d0 .and. delthe <= 0d0 )

    nwall = 0
    do l = 1, nbound
      if ( Is_Wall(bc(l) % type) ) nwall = nwall + 1
    enddo
    allocate ( wall_corner(3, 0:4, max(nwall,1)) )
    allocate ( facex(max(nwall,1)), facey(max(nwall,1)), facez(max(nwall,1)), facekr(max(nwall,1)), facer(max(nwall,1)) )

    nwall=0

    ! bc vector processing.
    do l = 1, nbound

      if ( Is_Wall(bc(l) % type) ) then

        Bm = bc(l) % b
        Im = bc(l) % i
        Jm = bc(l) % j
        Km = bc(l) % k
        Fm = bc(l) % f
        nwall = nwall+1

        ! Boundary face nodes, in cyclic order around the face.
        corner = Face_Corners ( blk(Bm), Im, Jm, Km, Fm )

        if ( wall_meridional ) then
          ! Faces 1-4 of a wedge cell cross the meridional plane along a segment,
          ! whose ends are nodes 1 and 2 (the same x and r as nodes 4 and 3).
          if ( Fm > 4 ) error stop 'Compute_Yn: a wall on a face of the axisymmetric wedge (faces 5-6) is not possible'
          wall_corner(:,1,nwall) = Meridional ( corner(:,1) )
          wall_corner(:,2,nwall) = Meridional ( corner(:,2) )
          wall_corner(:,3,nwall) = wall_corner(:,2,nwall)
          wall_corner(:,4,nwall) = wall_corner(:,1,nwall)
        else
          wall_corner(:,1:4,nwall) = corner
        end if

        ! Face centre (mean of the corners) and radius: every point of the face lies
        ! within the radius of the centre, with a margin for round-off.
        wall_corner(:,0,nwall) = 0.25d0 * ( wall_corner(:,1,nwall) + wall_corner(:,2,nwall) + &
                                            wall_corner(:,3,nwall) + wall_corner(:,4,nwall) )
        facex(nwall) = wall_corner(1,0,nwall)
        facey(nwall) = wall_corner(2,0,nwall)
        facez(nwall) = wall_corner(3,0,nwall)
        facer(nwall) = 0d0
        do n = 1, 4
          facer(nwall) = max ( facer(nwall), norm2 ( wall_corner(:,n,nwall) - wall_corner(:,0,nwall) ) )
        enddo
        facer(nwall) = facer(nwall) * ( 1d0 + 1d-10 )

        ! Face sand-grain roughness, inherited by the cells it is nearest to.
        facekr(nwall) = bc(l) % k_rough

      endif

    enddo

    ! After bc processing nwall is the viscous face counter.
    !
    ! Blocks needing yn here: the locally-owned ones, the remote sources of a
    ! local block-connection BC (Yn_Connection reads blk(Bs)%yn), and on root
    ! all of them, for the eddy viscosity it assembles for the solution file.
    ! The rest keep the zero First_Touch_Block wrote and are read by no one.
    needs_yn = .false.
    do b = 1, nb
      if ( is_local_block(b) ) needs_yn(b) = .true.
    enddo
    do l = 1, nbound
      if ( (bc(l)%type == 101 .or. bc(l)%type == 201) .and. is_local_block(bc(l)%b) ) &
        needs_yn(bc(l)%bs) = .true.
    enddo
    if ( mpi_is_root ) needs_yn = .true.

    ! No wall anywhere in the domain: leave behind the same sentinel the
    ! brute-force scan did, rather than the k-d tree's huge().
    if ( nwall == 0 ) then
      do b = 1, nb
        if (.not. needs_yn(b)) cycle
        blk(b) % yn ( 1:blk(b)%dim(1), 1:blk(b)%dim(2), 1:blk(b)%dim(3) ) = 1d8
      enddo
      deallocate ( wall_corner )
      return
    end if

    ! Periodic images of the wall. Across a periodic connection (201), and across type
    ! 200 on faces 5-6 of a 3-D mesh without sectors, the mesh goes on translated; across
    ! the rotational periodicity of a 3-D sector (200 on faces 5-6) it goes on turned by
    ! delthe about the x axis: the wall nearest to a cell can be an image of the wall
    ! across those boundaries, and the distance is the minimum over the images too.
    call Wall_Images ( )

    ! One k-d tree over the wall-face centres, with the face radii, shared by every
    ! block and thread: it hands the exact point-face distance only the faces that
    ! can still be nearer than the nearest found so far (centre within that
    ! distance plus the face radius).
    call kdtree_build ( wall_tree, facex(1:nwall), facey(1:nwall), facez(1:nwall), nwall, facer(1:nwall) )

    ! Compute distance in block b by querying the shared k-d tree.
    do b = 1, nb
      if (.not. needs_yn(b)) cycle
      call Wall_Distance_Blk ( blk(b), wall_tree, facekr(1:nwall) )
    enddo

    ! Validation, off unless DEBUG_KDTREE: prints max |tree - brute| per block
    ! and aborts past KDTREE_TOL.
    if ( DEBUG_KDTREE ) then
      do b = 1, nb
        if (.not. needs_yn(b)) cycle
        call Validate_Tree_Yn ( blk(b), nwall, b )
      enddo
    end if

    call kdtree_free ( wall_tree )
    deallocate ( wall_corner, img_shift, img_angle )

    contains

      !> Corners of face Fm of cell (Im, Jm, Km) of a block, in cyclic order around the face.
      function Face_Corners ( blk, Im, Jm, Km, Fm ) result ( corner )
        implicit none
        type(MOSE_block_type), intent(in) :: blk
        integer, intent(in)               :: Im, Jm, Km, Fm
        real(R8)                          :: corner(3,4)
        ! Local
        integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4

        select case(Fm)
          case(1:2)
            i1 = Im - mod(Fm,2)
            j1 = Jm - 1
            k1 = Km - 1
            i2 = Im - mod(Fm,2)
            j2 = Jm
            k2 = Km - 1
            i3 = Im - mod(Fm,2)
            j3 = Jm - 1
            k3 = Km
            i4 = Im - mod(Fm,2)
            j4 = Jm
            k4 = Km
          case(3:4)
            i1 = Im - 1
            j1 = Jm - mod(Fm,2)
            k1 = Km - 1
            i2 = Im
            j2 = Jm - mod(Fm,2)
            k2 = Km - 1
            i3 = Im - 1
            j3 = Jm - mod(Fm,2)
            k3 = Km
            i4 = Im
            j4 = Jm - mod(Fm,2)
            k4 = Km
          case(5:6)
            i1 = Im - 1
            j1 = Jm - 1
            k1 = Km - mod(Fm,2)
            i2 = Im
            j2 = Jm - 1
            k2 = Km - mod(Fm,2)
            i3 = Im - 1
            j3 = Jm
            k3 = Km - mod(Fm,2)
            i4 = Im
            j4 = Jm
            k4 = Km - mod(Fm,2)
        end select

        corner(:,1) = blk % node(i1,j1,k1) % c
        corner(:,2) = blk % node(i2,j2,k2) % c
        corner(:,3) = blk % node(i4,j4,k4) % c
        corner(:,4) = blk % node(i3,j3,k3) % c

      end function Face_Corners


      !> Images of the wall across the periodic boundaries (img_shift, img_angle, nimg) and the box around the wall.
      !>  - Translations: each periodic connection (201) gives the vector from its source face to its own face (the
      !>    difference of the face centres), in the meridional plane on a wedge, where it must lie along x. On a 3-D mesh
      !>    without sectors (delthe = 0) so does each block whose k faces are closed by type 200: its last k node plane
      !>    is the first one translated (Check_Mesh_Type), and the cells next to face 5 go on in those next to face 6
      !>    (Ghost_Axysimmetry). The two sides of a periodic boundary match only to within the precision of the mesh
      !>    file, so vectors that differ by less than MATCH times the size of the coordinates are the same (the tolerance
      !>    of Check_Mesh_Type). The distinct vectors, shortest first, are reduced to a basis of the lattice they
      !>    generate; a vector that is not an integer combination of the basis, within that tolerance times the number
      !>    of vectors combined, stops the run. Images: the integer combinations sum k_i B_i with |k_i| <= 2 D |B*_i|,
      !>    D the diagonal of the box around the nodes and B* the dual basis: a farther image is farther than D from
      !>    every cell, beyond the wall of the mesh itself.
      !>  - Rotations: on a 3-D mesh whose sectors are closed by type 200 on faces 5-6 (delthe /= 0), the turns by
      !>    m delthe, |m delthe| <= pi; then every translation must lie along x, the axis of the turns.
      !>  - More than NIMG_MAX images stop the run: the translations would be nearly parallel, or tiny against the mesh.
      subroutine Wall_Images ( )
        implicit none
        ! Local
        integer, parameter    :: MAXV = 3
        real(R8), parameter   :: MATCH = 1d-6, NIMG_MAX = 1d5
        real(R8), allocatable :: raw(:,:), cand(:,:), clen(:)
        integer,  allocatable :: order(:)
        logical,  allocatable :: closed(:)
        real(R8)              :: basis(3,MAXV), ginv(MAXV,MAXV), coef(MAXV), rhs(MAXV)
        real(R8)              :: t(3), lo(3), hi(3), cface(3,4), diag, scale, tol, lim, bound, nimg_r
        integer               :: nraw, nc, nv, kmax(MAXV), mrot, i, j, m, n, l, b, k1, k2, k3, kk(MAXV)
        logical               :: sector, span
        character(len=256)    :: msg

        ! Box around the wall faces, and around every node, for the bound on the images
        wall_box(:,1) = minval ( reshape ( wall_corner(:,1:4,1:nwall), [3, 4*nwall] ), dim=2 )
        wall_box(:,2) = maxval ( reshape ( wall_corner(:,1:4,1:nwall), [3, 4*nwall] ), dim=2 )
        lo = huge(1d0)
        hi = -huge(1d0)
        do b = 1, nb
          do n = 1, 3
            lo(n) = min ( lo(n), minval ( blk(b) % node(:,:,:) % c(n) ) )
            hi(n) = max ( hi(n), maxval ( blk(b) % node(:,:,:) % c(n) ) )
          enddo
        enddo
        scale = maxval ( max ( abs(lo), abs(hi) ) )
        if ( wall_meridional ) then
          ! (x, r): r runs from 0 to the largest |(y, z)|
          hi(2) = sqrt ( max ( lo(2)**2, hi(2)**2 ) + max ( lo(3)**2, hi(3)**2 ) )
          lo(2) = 0d0
          lo(3) = 0d0
          hi(3) = 0d0
        end if
        diag = norm2 ( hi - lo )
        tol  = MATCH * scale

        ! Translations: of the periodic connections, and of the blocks closed in k by type 200 on a 3-D mesh without
        ! sectors
        span = ndir == 3 .and. delthe >= 0d0 .and. delthe <= 0d0
        nraw = 0
        do l = 1, nbound
          if ( bc(l) % type == 201 ) nraw = nraw + 1
        enddo
        if ( span ) nraw = nraw + nb
        allocate ( raw(3, max(nraw,1)), cand(3, max(nraw,1)), clen(max(nraw,1)), order(max(nraw,1)) )
        nraw = 0
        do l = 1, nbound
          if ( bc(l) % type /= 201 ) cycle
          cface = Face_Corners ( blk(bc(l)%b), bc(l)%i, bc(l)%j, bc(l)%k, bc(l)%f )
          t = 0.25d0 * sum ( cface, dim=2 )
          cface = Face_Corners ( blk(bc(l)%bs), bc(l)%is, bc(l)%js, bc(l)%ks, bc(l)%fs )
          nraw = nraw + 1
          raw(:,nraw) = t - 0.25d0 * sum ( cface, dim=2 )
        enddo
        if ( span ) then
          allocate ( closed(nb) )
          closed = .false.
          do l = 1, nbound
            if ( bc(l) % type == 200 .and. bc(l) % f >= 5 ) closed(bc(l) % b) = .true.
          enddo
          do b = 1, nb
            if ( .not. closed(b) ) cycle
            nraw = nraw + 1
            raw(:,nraw) = blk(b) % node(0,0,blk(b)%dim(3)) % c - blk(b) % node(0,0,0) % c
          enddo
          deallocate ( closed )
        end if
        nc = 0
        do l = 1, nraw
          t = raw(:,l)
          if ( norm2 ( t ) <= tol ) cycle
          if ( wall_meridional ) then
            if ( abs ( t(2) ) > tol .or. abs ( t(3) ) > tol ) &
              call mpi_abort_all ( '[ERROR] Wall distance: a periodic connection (201) of an axisymmetric wedge is not '// &
                                   'a translation along the x axis' )
            t(2:3) = 0d0
          end if
          ! One of the vectors t and -t: the largest component positive
          n = maxloc ( abs ( t ), dim=1 )
          if ( t(n) < 0d0 ) t = -t
          do i = 1, nc
            if ( norm2 ( cand(:,i) - t ) <= tol ) exit
          enddo
          if ( i > nc ) then
            nc = nc + 1
            cand(:,nc) = t
            clen(nc) = norm2 ( t )
          end if
        enddo

        ! Basis of the lattice, shortest vectors first
        do i = 1, nc
          order(i) = i
        enddo
        do i = 2, nc
          j = order(i)
          m = i - 1
          do while ( m >= 1 )
            if ( clen(order(m)) <= clen(j) ) exit
            order(m+1) = order(m)
            m = m - 1
          enddo
          order(m+1) = j
        enddo
        nv = 0
        do i = 1, nc
          t = cand(:,order(i))
          if ( nv > 0 ) then
            ! Least-squares coefficients of t on the basis, and the nearest integer ones
            do m = 1, nv
              rhs(m) = dot_product ( basis(:,m), t )
            enddo
            call Invert_Gram ( basis, nv, ginv )
            coef(1:nv) = matmul ( ginv(1:nv,1:nv), rhs(1:nv) )
            lim = tol * ( 1d0 + sum ( abs ( anint ( coef(1:nv) ) ) ) )
            ! An integer combination of the basis: nothing new
            if ( norm2 ( t - matmul ( basis(:,1:nv), anint ( coef(1:nv) ) ) ) <= lim ) cycle
            ! In the span of the basis, between its integer combinations: the periodic boundaries disagree
            if ( norm2 ( t - matmul ( basis(:,1:nv), coef(1:nv) ) ) <= lim ) then
              write(msg,'(A,3ES12.4,A)') '[ERROR] Wall distance: the periodic translation (', t, &
                ') is not an integer combination of the others'
              call mpi_abort_all ( trim(msg) )
            end if
          end if
          if ( nv == MAXV ) call mpi_abort_all ( '[ERROR] Wall distance: more than three independent periodic translations' )
          nv = nv + 1
          basis(:,nv) = t
        enddo

        ! Rotational periodicity of 3-D sectors
        sector = .false.
        if ( ndir == 3 .and. .not. ( delthe >= 0d0 .and. delthe <= 0d0 ) ) then
          do l = 1, nbound
            if ( bc(l) % type == 200 .and. bc(l) % f >= 5 ) sector = .true.
          enddo
        end if
        if ( sector ) then
          do m = 1, nv
            if ( abs ( basis(2,m) ) > tol .or. abs ( basis(3,m) ) > tol ) &
              call mpi_abort_all ( '[ERROR] Wall distance: a periodic connection (201) of a mesh with rotational '// &
                                   'periodicity (200) is not a translation along the x axis' )
          enddo
        end if

        ! Bounds of the integer coefficients, |k_i| <= 2 D |B*_i| with |B*_i|^2 = (G^-1)_ii, and of the turns, and the
        ! number of images (in reals: no integer overflow before the check)
        kmax = 0
        mrot = 0
        nimg_r = 1d0
        if ( nv > 0 ) then
          call Invert_Gram ( basis, nv, ginv )
          do m = 1, nv
            bound = 2d0 * diag * sqrt ( ginv(m,m) )
            if ( .not. ( bound <= NIMG_MAX ) ) bound = NIMG_MAX   ! too large, or not a number (Gram matrix singular)
            kmax(m) = floor ( bound ) + 1
            nimg_r = nimg_r * ( 2 * kmax(m) + 1 )
          enddo
        end if
        if ( sector ) then
          bound = 3.141592653589793d0 / abs ( delthe )
          if ( .not. ( bound <= NIMG_MAX ) ) bound = NIMG_MAX
          mrot = ceiling ( bound )
          nimg_r = nimg_r * ( 2 * mrot + 1 )
        end if
        if ( nimg_r - 1d0 > NIMG_MAX ) then
          write(msg,'(A,ES10.3,A)') '[ERROR] Wall distance: ', nimg_r - 1d0, ' periodic images of the wall: the periodic'// &
            ' translations are nearly parallel, or tiny against the mesh (do the two sides of the periodic boundaries match?)'
          call mpi_abort_all ( trim(msg) )
        end if

        ! The images, all but the identity
        nimg = product ( 2*kmax + 1 ) * ( 2*mrot + 1 ) - 1
        allocate ( img_shift(3, max(nimg,1)), img_angle(max(nimg,1)) )
        img_angle = 0d0
        nimg = 0
        do k3 = -kmax(3), kmax(3)
          do k2 = -kmax(2), kmax(2)
            do k1 = -kmax(1), kmax(1)
              kk = [ k1, k2, k3 ]
              t = 0d0
              do m = 1, nv
                t = t + kk(m) * basis(:,m)
              enddo
              do m = -mrot, mrot
                if ( k1 == 0 .and. k2 == 0 .and. k3 == 0 .and. m == 0 ) cycle
                nimg = nimg + 1
                img_shift(:,nimg) = t
                img_angle(nimg) = m * delthe
              enddo
            enddo
          enddo
        enddo
        deallocate ( raw, cand, clen, order )

      end subroutine Wall_Images


      !> ginv = inverse of the Gram matrix of basis(:,1:nv), nv <= 3 (independent vectors).
      subroutine Invert_Gram ( basis, nv, ginv )
        implicit none
        real(R8), intent(in)  :: basis(:,:)
        integer, intent(in)   :: nv
        real(R8), intent(out) :: ginv(:,:)
        ! Local
        integer  :: i, j
        real(R8) :: gram(3,3), det

        do i = 1, nv
          do j = 1, nv
            gram(i,j) = dot_product ( basis(:,i), basis(:,j) )
          enddo
        enddo
        select case ( nv )
          case ( 1 )
            ginv(1,1) = 1d0 / gram(1,1)
          case ( 2 )
            det = gram(1,1) * gram(2,2) - gram(1,2) * gram(2,1)
            ginv(1,1) =  gram(2,2) / det
            ginv(2,2) =  gram(1,1) / det
            ginv(1,2) = -gram(1,2) / det
            ginv(2,1) = -gram(2,1) / det
          case ( 3 )
            ginv(1,1) = gram(2,2)*gram(3,3) - gram(2,3)*gram(3,2)
            ginv(1,2) = gram(1,3)*gram(3,2) - gram(1,2)*gram(3,3)
            ginv(1,3) = gram(1,2)*gram(2,3) - gram(1,3)*gram(2,2)
            ginv(2,1) = gram(2,3)*gram(3,1) - gram(2,1)*gram(3,3)
            ginv(2,2) = gram(1,1)*gram(3,3) - gram(1,3)*gram(3,1)
            ginv(2,3) = gram(1,3)*gram(2,1) - gram(1,1)*gram(2,3)
            ginv(3,1) = gram(2,1)*gram(3,2) - gram(2,2)*gram(3,1)
            ginv(3,2) = gram(1,2)*gram(3,1) - gram(1,1)*gram(3,2)
            ginv(3,3) = gram(1,1)*gram(2,2) - gram(1,2)*gram(2,1)
            det = gram(1,1)*ginv(1,1) + gram(1,2)*ginv(2,1) + gram(1,3)*ginv(3,1)
            ginv(1:3,1:3) = ginv(1:3,1:3) / det
        end select

      end subroutine Invert_Gram


      !> Distance from p to the box wall_box (0 inside).
      pure function Box_Distance ( p ) result ( d )
        implicit none
        real(R8), intent(in) :: p(3)
        real(R8)             :: d

        d = norm2 ( max ( wall_box(:,1) - p, 0d0, p - wall_box(:,2) ) )

      end function Box_Distance


      !> The point at which the image n of the wall is seen as the wall itself: p translated by -img_shift(:,n) and
      !> turned by -img_angle(n) about the x axis.
      pure function Image_Point ( p, n ) result ( q )
        implicit none
        real(R8), intent(in) :: p(3)
        integer, intent(in)  :: n
        real(R8)             :: q(3)
        ! Local
        real(R8) :: ca, sa

        q = p - img_shift(:,n)
        if ( .not. ( img_angle(n) >= 0d0 .and. img_angle(n) <= 0d0 ) ) then
          ca = cos ( img_angle(n) )
          sa = sin ( img_angle(n) )
          q(2:3) = [ q(2) * ca + q(3) * sa, - q(2) * sa + q(3) * ca ]
        end if

      end function Image_Point


      subroutine Wall_Distance_Blk ( blk, tree, face_kr )
        use MOSE_Lib_KDTree, only: kdtree_t, kdtree_nearest_object
        implicit none
        type(MOSE_block_type), intent(inout) :: blk
        type(kdtree_t), intent(in)           :: tree
        real(R8), intent(in)                 :: face_kr(:)   ! Roughness of each wall face, by original face index
        ! Local
        integer :: i, j, k, iface, n
        real(R8) :: center(3), dist, q(3)

        !$omp parallel
        !$omp do collapse(3) private(i, j, k, center, dist, iface, n, q)
        do k = 1, blk % dim(3)
        do j = 1, blk % dim(2)
        do i = 1, blk % dim(1)

          center = Cell_Centre ( blk, i, j, k )

          ! Nearest wall face via the k-d tree, exact point-face distance.
          call kdtree_nearest_object ( tree, center(1), center(2), center(3), Wall_Face_Distance, dist, iface )

          ! Periodic images: an image can be nearer only if the box of the wall, seen
          ! from the moved point, is not farther than the distance found so far. The
          ! query goes on from that distance and face, so ties keep the lowest face.
          do n = 1, nimg
            q = Image_Point ( center, n )
            if ( Box_Distance ( q ) > dist ) cycle
            call kdtree_nearest_object ( tree, q(1), q(2), q(3), Wall_Face_Distance, dist, iface, resume=.true. )
          enddo

          blk % yn (i,j,k) = dist
          blk % k_rough (i,j,k) = face_kr(iface)

        enddo; enddo; enddo
        !$omp end parallel

      end subroutine Wall_Distance_Blk


      !> Brute-force check of the k-d tree result: recompute each cell's
      !> nearest-wall distance against the full face list and its periodic
      !> images, compare with `blk%yn`, and abort past KDTREE_TOL.
      subroutine Validate_Tree_Yn ( blk, nwall_, b_id )
        use MOSE_Mod_MPI, only: mpi_rank_
        implicit none
        type(MOSE_block_type), intent(in)           :: blk
        integer, intent(in)                         :: nwall_
        integer, intent(in)                         :: b_id
        ! Local
        integer :: i, j, k, ll, n
        real(R8) :: center(3), q(3), dbrute, dummy, err, max_err

        max_err = 0d0

        !$omp parallel
        !$omp do collapse(3) private(i, j, k, ll, n, center, q, dbrute, dummy, err) &
        !$omp             reduction(max:max_err)
        do k = 1, blk % dim(3)
        do j = 1, blk % dim(2)
        do i = 1, blk % dim(1)

          center = Cell_Centre ( blk, i, j, k )

          dbrute = 1d8
          do n = 0, nimg
            q = center
            if ( n > 0 ) q = Image_Point ( center, n )
            do ll = 1, nwall_
              dummy = Wall_Face_Distance ( ll, q(1), q(2), q(3) )
              if (dummy < dbrute) dbrute = dummy
            enddo
          enddo

          err = abs( blk % yn(i,j,k) - dbrute )
          if (err > max_err) max_err = err

        enddo; enddo; enddo
        !$omp end parallel

        write(*,'(A,I0,A,I0,A,ES12.5)') &
          ' [KDTREE] rank ', mpi_rank_, '  block ', b_id, '  max |tree - brute| = ', max_err
        if (max_err > KDTREE_TOL) then
          write(*,'(A,ES12.5,A)') ' [KDTREE] ERROR: exceeds tolerance ', KDTREE_TOL, '. Aborting.'
          stop 1
        end if

      end subroutine Validate_Tree_Yn

  end subroutine Compute_Yn


  !> Point where the wall distance of cell (i,j,k) is measured: the mean of its
  !> eight nodes, in the meridional plane (x, r, 0) on an axisymmetric wedge.
  function Cell_Centre ( blk, i, j, k ) result ( center )
    implicit none
    type(MOSE_block_type), intent(in) :: blk
    integer, intent(in)               :: i, j, k
    real(R8)                          :: center(3)
    ! Local
    integer :: a, b, c

    center = 0d0
    do c = k-1, k
    do b = j-1, j
    do a = i-1, i
      if ( wall_meridional ) then
        center = center + Meridional ( blk % node(a,b,c) % c )
      else
        center = center + blk % node(a,b,c) % c
      end if
    enddo; enddo; enddo
    center = center / 8d0

  end function Cell_Centre


  !> (x, y, z) -> (x, r, 0), with r = sqrt(y^2 + z^2) the distance from the x axis.
  pure function Meridional ( p ) result ( q )
    implicit none
    real(R8), intent(in) :: p(3)
    real(R8)             :: q(3)

    q = [ p(1), sqrt ( p(2)**2 + p(3)**2 ), 0d0 ]

  end function Meridional


  !> Distance from (x,y,z) to wall face i of Compute_Yn (kdtree_object_distance):
  !> the segment of the face in the meridional plane on an axisymmetric wedge;
  !> otherwise the four triangles between the face centre and the sides, which are
  !> the face itself when it is plane, as on every 2-D and wedge mesh.
  function Wall_Face_Distance ( i, x, y, z ) result ( d )
    implicit none
    integer, intent(in)  :: i
    real(R8), intent(in) :: x, y, z
    real(R8)             :: d
    ! Local
    real(R8) :: p(3)

    p = [ x, y, z ]
    if ( wall_meridional ) then
      d = Point_Segment_Distance ( p, wall_corner(:,1,i), wall_corner(:,2,i) )
    else
      d = min ( Point_Triangle_Distance ( p, wall_corner(:,0,i), wall_corner(:,1,i), wall_corner(:,2,i) ), &
                Point_Triangle_Distance ( p, wall_corner(:,0,i), wall_corner(:,2,i), wall_corner(:,3,i) ), &
                Point_Triangle_Distance ( p, wall_corner(:,0,i), wall_corner(:,3,i), wall_corner(:,4,i) ), &
                Point_Triangle_Distance ( p, wall_corner(:,0,i), wall_corner(:,4,i), wall_corner(:,1,i) ) )
    end if

  end function Wall_Face_Distance


  !> Distance from p to the segment ab (a point if a = b).
  pure function Point_Segment_Distance ( p, a, b ) result ( d )
    implicit none
    real(R8), intent(in) :: p(3), a(3), b(3)
    real(R8)             :: d
    ! Local
    real(R8) :: ab(3), ap(3), l2, t

    ab = b - a
    ap = p - a
    l2 = dot_product ( ab, ab )
    t  = 0d0
    if ( l2 > 0d0 ) t = min ( max ( dot_product ( ap, ab ) / l2, 0d0 ), 1d0 )
    d = norm2 ( ap - t * ab )

  end function Point_Segment_Distance


  !> Distance from p to the triangle abc. The triangle is convex: if the projection
  !> of p on its plane falls inside it (on the inner side of the three sides), the
  !> distance is the one from the plane; otherwise the nearest point lies on a side.
  !> A triangle with no area next to its sides (a degenerate face) has no plane and
  !> is treated as its sides.
  pure function Point_Triangle_Distance ( p, a, b, c ) result ( d )
    implicit none
    real(R8), intent(in) :: p(3), a(3), b(3), c(3)
    real(R8)             :: d
    ! Local
    real(R8) :: n(3), n2, e2

    n  = Cross3 ( b - a, c - a )
    n2 = dot_product ( n, n )
    e2 = max ( sum ( (b - a)**2 ), sum ( (c - b)**2 ), sum ( (a - c)**2 ) )
    if ( n2 > ( 1d-12 * e2 )**2 ) then
      if ( dot_product ( Cross3 ( b - a, p - a ), n ) >= 0d0 .and. &
           dot_product ( Cross3 ( c - b, p - b ), n ) >= 0d0 .and. &
           dot_product ( Cross3 ( a - c, p - c ), n ) >= 0d0 ) then
        d = abs ( dot_product ( p - a, n ) ) / sqrt ( n2 )
        return
      end if
    end if
    d = min ( Point_Segment_Distance ( p, a, b ), Point_Segment_Distance ( p, b, c ), &
              Point_Segment_Distance ( p, c, a ) )

  end function Point_Triangle_Distance


  pure function Cross3 ( u, v ) result ( w )
    implicit none
    real(R8), intent(in) :: u(3), v(3)
    real(R8)             :: w(3)

    w = [ u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1) ]

  end function Cross3
  

  subroutine Yn_Connection ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, blks )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm, Is, Js, Ks
    type(MOSE_block_type), intent(inout) :: blks
    type(MOSE_block_type), intent(inout) :: blkm
    ! Local
    integer :: Ig, Jg, Kg

    ! Ghost cell coordinates
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    blkm % yn(Ig,Jg,Kg) = blks % yn(Is,Js,Ks)
    blkm % k_rough(Ig,Jg,Kg) = blks % k_rough(Is,Js,Ks)

  end subroutine Yn_Connection


  subroutine Yn_Chimera ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: blk
    ! Local
    integer :: Ig, Jg, Kg

    ! Ghost cell coordinates
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    ! The chimera ghost cell extends the receiver grid (as its metrics do), so it
    ! takes the wall distance and roughness of the adjacent receiver cell.
    ! Zeroth order on purpose: linear extrapolation can reach yn <= 0 near walls.
    blk % yn(Ig,Jg,Kg) = blk % yn(Im,Jm,Km)
    blk % k_rough(Ig,Jg,Kg) = blk % k_rough(Im,Jm,Km)

  end subroutine Yn_Chimera


  subroutine BC_Extrapolate_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: blk
    type(MOSE_tensor_3D_type), intent(out)    :: Mg(2)
    type(MOSE_vector_3D_type), intent(out)    :: dlg(2)
    real(R8), intent(out)                :: volg(2)
    ! Local
    integer :: g
    type(MOSE_vector_3D_type), dimension(-gc:gc):: N1, N2, N3, N4, N5, N6, N7, N8

    select case (Fm)
      case(1)
        do g = 0, gc
          N1(g)%c = blk % node (Im+g-1,Jm-1,Km-1) % c
          N2(g)%c = blk % node (Im+g-1,Jm-1,Km  ) % c
          N3(g)%c = blk % node (Im+g-1,Jm  ,Km-1) % c
          N4(g)%c = blk % node (Im+g-1,Jm  ,Km  ) % c
        enddo
        do g = -1, -gc, -1
          N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
          N2(g)%c = 3d0*N2(g+1)%c - 3d0*N2(g+2)%c + N2(g+3)%c
          N3(g)%c = 3d0*N3(g+1)%c - 3d0*N3(g+2)%c + N3(g+3)%c
          N4(g)%c = 3d0*N4(g+1)%c - 3d0*N4(g+2)%c + N4(g+3)%c
          call Compute_Metric_Tensor ( N1(g), N2(g), N3(g), N4(g), N1(g+1), N2(g+1), N3(g+1), N4(g+1), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km))
          Mg(-g)   = blk % M(Im+g, Jm, Km)
          dlg(-g)  = blk % dl(Im+g, Jm, Km)
          volg(-g) = blk % vol(Im+g, Jm, Km)
        enddo
        
      case(2)
        do g = 0, -gc, -1
          N5(g)%c = blk % node (Im+g,Jm-1,Km-1) % c
          N6(g)%c = blk % node (Im+g,Jm-1,Km  ) % c
          N7(g)%c = blk % node (Im+g,Jm  ,Km-1) % c
          N8(g)%c = blk % node (Im+g,Jm  ,Km  ) % c
        enddo
        do g = 1, gc
          N5(g)%c = 3d0*N5(g-1)%c - 3d0*N5(g-2)%c + N5(g-3)%c
          N6(g)%c = 3d0*N6(g-1)%c - 3d0*N6(g-2)%c + N6(g-3)%c
          N7(g)%c = 3d0*N7(g-1)%c - 3d0*N7(g-2)%c + N7(g-3)%c
          N8(g)%c = 3d0*N8(g-1)%c - 3d0*N8(g-2)%c + N8(g-3)%c
          call Compute_Metric_Tensor ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km))
          Mg(g)   = blk % M(Im+g, Jm, Km)
          dlg(g)  = blk % dl(Im+g, Jm, Km)
          volg(g) = blk % vol(Im+g, Jm, Km)
        enddo

      case(3)
        if ( blk % dim(2) > 2 ) then
          do g = 0, gc
            N1(g)%c = blk % node (Im-1,Jm+g-1,Km-1) % c
            N2(g)%c = blk % node (Im-1,Jm+g-1,Km  ) % c
            N5(g)%c = blk % node (Im  ,Jm+g-1,Km-1) % c
            N6(g)%c = blk % node (Im  ,Jm+g-1,Km  ) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
            N2(g)%c = 3d0*N2(g+1)%c - 3d0*N2(g+2)%c + N2(g+3)%c
            N5(g)%c = 3d0*N5(g+1)%c - 3d0*N5(g+2)%c + N5(g+3)%c
            N6(g)%c = 3d0*N6(g+1)%c - 3d0*N6(g+2)%c + N6(g+3)%c
            call Compute_Metric_Tensor ( N1(g), N2(g), N1(g+1), N2(g+1), N5(g), N6(g), N5(g+1), N6(g+1), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(-g)   = blk % M(Im, Jm+g, Km)
            dlg(-g)  = blk % dl(Im, Jm+g, Km)
            volg(-g) = blk % vol(Im, Jm+g, Km)
          enddo
        else
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm+g-1,Km-1) % c
            N2(g)%c = blk % node (Im-1,Jm+g-1,Km  ) % c
            N5(g)%c = blk % node (Im  ,Jm+g-1,Km-1) % c
            N6(g)%c = blk % node (Im  ,Jm+g-1,Km  ) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
            N2(g)%c = 2d0*N2(g+1)%c - N2(g+2)%c
            N5(g)%c = 2d0*N5(g+1)%c - N5(g+2)%c
            N6(g)%c = 2d0*N6(g+1)%c - N6(g+2)%c
            call Compute_Metric_Tensor ( N1(g), N2(g), N1(g+1), N2(g+1), N5(g), N6(g), N5(g+1), N6(g+1), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(-g)   = blk % M(Im, Jm+g, Km)
            dlg(-g)  = blk % dl(Im, Jm+g, Km)
            volg(-g) = blk % vol(Im, Jm+g, Km)
          enddo
        endif

      case(4)
        if ( blk % dim(2) > 2 ) then
          do g = 0, -gc, -1
            N3(g)%c = blk % node (Im-1,Jm+g,Km-1) % c
            N4(g)%c = blk % node (Im-1,Jm+g,Km  ) % c
            N7(g)%c = blk % node (Im  ,Jm+g,Km-1) % c
            N8(g)%c = blk % node (Im  ,Jm+g,Km  ) % c
          enddo
          do g = 1, gc
            N3(g)%c = 3d0*N3(g-1)%c - 3d0*N3(g-2)%c + N3(g-3)%c
            N4(g)%c = 3d0*N4(g-1)%c - 3d0*N4(g-2)%c + N4(g-3)%c
            N7(g)%c = 3d0*N7(g-1)%c - 3d0*N7(g-2)%c + N7(g-3)%c
            N8(g)%c = 3d0*N8(g-1)%c - 3d0*N8(g-2)%c + N8(g-3)%c
            call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(g)   = blk % M(Im, Jm+g, Km)
            dlg(g)  = blk % dl(Im, Jm+g, Km)
            volg(g) = blk % vol(Im, Jm+g, Km)
          enddo
        else
          do g = 0, -1, -1
            N3(g)%c = blk % node (Im-1,Jm+g,Km-1) % c
            N4(g)%c = blk % node (Im-1,Jm+g,Km  ) % c
            N7(g)%c = blk % node (Im  ,Jm+g,Km-1) % c
            N8(g)%c = blk % node (Im  ,Jm+g,Km  ) % c
          enddo
          do g = 1, gc
            N3(g)%c = 2d0*N3(g-1)%c - N3(g-2)%c
            N4(g)%c = 2d0*N4(g-1)%c - N4(g-2)%c
            N7(g)%c = 2d0*N7(g-1)%c - N7(g-2)%c
            N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
            call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(g)   = blk % M(Im, Jm+g, Km)
            dlg(g)  = blk % dl(Im, Jm+g, Km)
            volg(g) = blk % vol(Im, Jm+g, Km)
          enddo
        endif

      case(5)
        if ( ndir==1 ) then
          do g = 1, gc
            Mg(g)   = blk % M(Im, Jm, Km)
            dlg(g)  = blk % dl(Im, Jm, Km)
            volg(g) = blk % vol(Im, Jm, Km)
          enddo
        elseif ( ndir==2 .and. delthe==0d0 ) then
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
            N3(g)%c = 2d0*N3(g+1)%c - N3(g+2)%c
            N5(g)%c = 2d0*N5(g+1)%c - N5(g+2)%c
            N7(g)%c = 2d0*N7(g+1)%c - N7(g+2)%c
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo
        elseif (ndir==2 .and. delthe/=0d0) then
          ! 2Dax: extrapolation with a rotation angle delthe (exact).
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c(1) = N1(0) % c(1) 
            N1(g)%c(2) = N1(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N1(g)%c(3) = N1(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N3(g)%c(1) = N3(0) % c(1) 
            N3(g)%c(2) = N3(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N3(g)%c(3) = N3(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N5(g)%c(1) = N5(0) % c(1) 
            N5(g)%c(2) = N5(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N5(g)%c(3) = N5(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N7(g)%c(1) = N7(0) % c(1) 
            N7(g)%c(2) = N7(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N7(g)%c(3) = N7(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo         
        elseif (ndir==3) then
          do g = 0, gc
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
            N3(g)%c = 3d0*N3(g+1)%c - 3d0*N3(g+2)%c + N3(g+3)%c
            N5(g)%c = 3d0*N5(g+1)%c - 3d0*N5(g+2)%c + N5(g+3)%c
            N7(g)%c = 3d0*N7(g+1)%c - 3d0*N7(g+2)%c + N7(g+3)%c
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      case(6)
        if ( ndir==1 ) then
          do g = 1, gc
            Mg(g)   = blk % M(Im, Jm, Km)
            dlg(g)  = blk % dl(Im, Jm, Km)
            volg(g) = blk % vol(Im, Jm, Km)
          enddo
        elseif ( ndir==2 .and. delthe==0d0 ) then
          do g = 0, -1, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c = 2d0*N2(g-1)%c - N2(g-2)%c
            N4(g)%c = 2d0*N4(g-1)%c - N4(g-2)%c
            N6(g)%c = 2d0*N6(g-1)%c - N6(g-2)%c
            N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo
        elseif (ndir==2 .and. delthe/=0d0) then
          ! 2Dax: extrapolation with a rotation angle delthe (exact).
          do g = 0, -1, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c(1) = N2(0) % c(1) 
            N2(g)%c(2) = N2(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N2(g)%c(3) = N2(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N4(g)%c(1) = N4(0) % c(1) 
            N4(g)%c(2) = N4(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N4(g)%c(3) = N4(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N6(g)%c(1) = N6(0) % c(1) 
            N6(g)%c(2) = N6(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N6(g)%c(3) = N6(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N8(g)%c(1) = N8(0) % c(1) 
            N8(g)%c(2) = N8(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N8(g)%c(3) = N8(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo         
        elseif (ndir==3) then
          do g = 0, -gc, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c = 3d0*N2(g-1)%c - 3d0*N2(g-2)%c + N2(g-3)%c
            N4(g)%c = 3d0*N4(g-1)%c - 3d0*N4(g-2)%c + N4(g-3)%c
            N6(g)%c = 3d0*N6(g-1)%c - 3d0*N6(g-2)%c + N6(g-3)%c
            N8(g)%c = 3d0*N8(g-1)%c - 3d0*N8(g-2)%c + N8(g-3)%c
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      end select

  end subroutine BC_Extrapolate_Metrics


  subroutine BC_Symmetry_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout)   :: blk
    type(MOSE_tensor_3D_type), intent(out) :: Mg(2)
    type(MOSE_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    ! Local
    integer :: g

    do g = 1, gc
      select case (Fm)
        case(1)
          blk % M(Im-g,Jm,Km)   = blk % M(Im+g-1,Jm,Km)
          blk % dl(Im-g,Jm,Km)  = blk % dl(Im+g-1,Jm,Km)
          blk % vol(Im-g,Jm,Km) = blk % vol(Im+g-1,Jm,Km)
          Mg(g)   = blk % M(Im-g,Jm,Km)
          dlg(g)  = blk % dl(Im-g,Jm,Km)
          volg(g) = blk % vol(Im-g,Jm,Km)
        case(2)
          blk % M(Im+g,Jm,Km)   = blk % M(Im-g+1,Jm,Km)
          blk % dl(Im+g,Jm,Km)  = blk % dl(Im-g+1,Jm,Km)
          blk % vol(Im+g,Jm,Km) = blk % vol(Im-g+1,Jm,Km)
          Mg(g)   = blk % M(Im+g,Jm,Km)
          dlg(g)  = blk % dl(Im+g,Jm,Km)
          volg(g) = blk % vol(Im+g,Jm,Km)
        case(3)
          blk % M(Im,Jm-g,Km)   = blk % M(Im,Jm+g-1,Km)
          blk % dl(Im,Jm-g,Km)  = blk % dl(Im,Jm+g-1,Km)
          blk % vol(Im,Jm-g,Km) = blk % vol(Im,Jm+g-1,Km)
          Mg(g)   = blk % M(Im,Jm-g,Km)
          dlg(g)  = blk % dl(Im,Jm-g,Km)
          volg(g) = blk % vol(Im,Jm-g,Km)
        case(4)
          blk % M(Im,Jm+g,Km)   = blk % M(Im,Jm-g+1,Km)
          blk % dl(Im,Jm+g,Km)  = blk % dl(Im,Jm-g+1,Km)
          blk % vol(Im,Jm+g,Km) = blk % vol(Im,Jm-g+1,Km)
          Mg(g)   = blk % M(Im,Jm+g,Km)
          dlg(g)  = blk % dl(Im,Jm+g,Km)
          volg(g) = blk % vol(Im,Jm+g,Km)
        case(5)
          blk % M(Im,Jm,Km-g)   = blk % M(Im,Jm,Km+g-1)
          blk % dl(Im,Jm,Km-g)  = blk % dl(Im,Jm,Km+g-1)
          blk % vol(Im,Jm,Km-g) = blk % vol(Im,Jm,Km+g-1)
          Mg(g)   = blk % M(Im,Jm,Km-g)
          dlg(g)  = blk % dl(Im,Jm,Km-g)
          volg(g) = blk % vol(Im,Jm,Km-g)
        case(6)
          blk % M(Im,Jm,Km+g)   = blk % M(Im,Jm,Km-g+1)
          blk % dl(Im,Jm,Km+g)  = blk % dl(Im,Jm,Km-g+1)
          blk % vol(Im,Jm,Km+g) = blk % vol(Im,Jm,Km-g+1)
          Mg(g)   = blk % M(Im,Jm,Km+g)
          dlg(g)  = blk % dl(Im,Jm,Km+g)
          volg(g) = blk % vol(Im,Jm,Km+g)
      end select
    enddo

  end subroutine BC_Symmetry_Metrics


  subroutine BC_Connect_Metrics ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, volg, periodic )
    implicit none
    integer, intent(in)                  :: Im, Jm, Km, Fm, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(MOSE_block_type), intent(in)    :: blks
    type(MOSE_block_type), intent(inout) :: blkm
    type(MOSE_tensor_3D_type), intent(out) :: Mg(2)
    type(MOSE_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    logical, intent(in)                    :: periodic
    ! Local
    integer      :: g, Is1, Js1, Ks1
    integer      :: II, JJ, KK, III, JJJ, KKK
    integer      :: guidem(6,3), guides(6,3), guidem2(6,3), guides2(6,3)
    type(MOSE_vector_3D_type), dimension(0:gc) :: N1, N2, N3, N4, N5, N6, N7, N8
    type(MOSE_vector_3D_type) :: T
    
    ! ---------------------------------------------------------------------------------------------
    ! Preliminary definitions
    guidem(1,1)=-2 ; guidem(1,2)= 0 ; guidem(1,3)= 0 ; guidem(2,1)=1 ; guidem(2,2)=0 ; guidem(2,3)=0
    guidem(3,1)= 0 ; guidem(3,2)=-2 ; guidem(3,3)= 0 ; guidem(4,1)=0 ; guidem(4,2)=1 ; guidem(4,3)=0
    guidem(5,1)= 0 ; guidem(5,2)=0  ; guidem(5,3)=-2 ; guidem(6,1)=0 ; guidem(6,2)=0 ; guidem(6,3)=1

    guidem2(1,1)=-3 ; guidem2(1,2)= 0 ; guidem2(1,3)= 0 ; guidem2(2,1)=2 ; guidem2(2,2)=0 ; guidem2(2,3)=0
    guidem2(3,1)= 0 ; guidem2(3,2)=-3 ; guidem2(3,3)= 0 ; guidem2(4,1)=0 ; guidem2(4,2)=2 ; guidem2(4,3)=0
    guidem2(5,1)= 0 ; guidem2(5,2)= 0 ; guidem2(5,3)=-3 ; guidem2(6,1)=0 ; guidem2(6,2)=0 ; guidem2(6,3)=2

    guides(1,1)= 0 ; guides(1,2)= 0 ; guides(1,3)= 0 ; guides(2,1)=-1 ; guides(2,2)= 0 ; guides(2,3)= 0
    guides(3,1)= 0 ; guides(3,2)= 0 ; guides(3,3)= 0 ; guides(4,1)= 0 ; guides(4,2)=-1 ; guides(4,3)= 0
    guides(5,1)= 0 ; guides(5,2)= 0 ; guides(5,3)= 0 ; guides(6,1)= 0 ; guides(6,2)= 0 ; guides(6,3)=-1

    guides2(1,1)= 1 ; guides2(1,2)= 0 ; guides2(1,3)= 0 ; guides2(2,1)=-2 ; guides2(2,2)= 0 ; guides2(2,3)= 0
    guides2(3,1)= 0 ; guides2(3,2)= 1 ; guides2(3,3)= 0 ; guides2(4,1)= 0 ; guides2(4,2)=-2 ; guides2(4,3)= 0
    guides2(5,1)= 0 ; guides2(5,2)= 0 ; guides2(5,3)= 1 ; guides2(6,1)= 0 ; guides2(6,2)= 0 ; guides2(6,3)=-2
    ! ---------------------------------------------------------------------------------------------

    select case (Fm)
      case(1)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
      case(2)
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(3)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
      case(4)
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(5)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
      case(6)
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
    end select

    T%c = 0._R8
    if (periodic) then
      Is1 = Is + gc*guides(Fs,1) - guides2(Fs,1)
      Js1 = Js + gc*guides(Fs,2) - guides2(Fs,2)
      Ks1 = Ks + gc*guides(Fs,3) - guides2(Fs,3)
      select case(Fs)
        case(1:2)
          II=Is1
          JJ=(2*Js1+d11s+d21s)/2
          KK=(2*Ks1+d12s+d22s)/2
        case(3:4)
          II=(2*Is1+d11s+d21s)/2
          JJ=Js1
          KK=(2*Ks1+d12s+d22s)/2
        case(5:6)
          II=(2*Is1+d11s+d21s)/2
          JJ=(2*Js1+d12s+d22s)/2
          KK=Ks1
      end select
      select case(Fm)
        case(1) ; T%c = N4(0)%c - blks%node(II,JJ,KK)%c
        case(2) ; T%c = N8(0)%c - blks%node(II,JJ,KK)%c
        case(3) ; T%c = N6(0)%c - blks%node(II,JJ,KK)%c
        case(4) ; T%c = N8(0)%c - blks%node(II,JJ,KK)%c
        case(5) ; T%c = N7(0)%c - blks%node(II,JJ,KK)%c
        case(6) ; T%c = N8(0)%c - blks%node(II,JJ,KK)%c
      end select
    end if

    do g = 1, gc

      ! connected cell node index
      Is1 = Is + (gc-g)*guides(Fs,1)+(g-1)*guides2(Fs,1)
      Js1 = Js + (gc-g)*guides(Fs,2)+(g-1)*guides2(Fs,2)
      Ks1 = Ks + (gc-g)*guides(Fs,3)+(g-1)*guides2(Fs,3)

      ! compute the indexes in the connected block corresponding to the ghost nodes to update in the local block 
      select case(Fs)
        case(1:2) ! i-faces
          II=Is1
          JJ=(2*Js1+d11s+d21s)/2
          KK=(2*Ks1+d12s+d22s)/2
        case(3:4) ! j-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=Js1
          KK=(2*Ks1+d12s+d22s)/2
        case(5:6) ! k-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=(2*Js1+d12s+d22s)/2
          KK=Ks1
      end select
      if (Fm == 1) then
          N4(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 2) then
          N8(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 3) then
          N6(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 4) then
          N8(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 5) then
          N7(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 6) then
          N8(g)%c = blks%node(ii,jj,kk)%c + T%c
      endif

      ! update the remaining 3 indexes (only necessary in the corners, but for good measure)
      ! node i1-1,i2 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s
          KKK=KK-d12s
        case(3:4)
          III=II-d11s
          JJJ=JJ
          KKK=KK-d12s
        case(5:6)
          III=II-d11s
          JJJ=JJ-d12s
          KKK=KK
      end select
      if (Fm == 1) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 2) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 3) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 4) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 5) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 6) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      endif

      ! node i1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d21s
          KKK=KK-d22s
        case(3:4)
          III=II-d21s
          JJJ=JJ
          KKK=KK-d22s
        case(5:6)
          III=II-d21s
          JJJ=JJ-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 2) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 3) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 4) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 5) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 6) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      endif

      ! node i1-1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s-d21s
          KKK=KK-d12s-d22s
        case(3:4)
          III=II-d11s-d21s
          JJJ=JJ
          KKK=KK-d12s-d22s
        case(5:6)
          III=II-d11s-d21s
          JJJ=JJ-d12s-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 2) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 3) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 4) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 5) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 6) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      endif

      ! compute metric variables
      select case(Fm)
        case(1)
          call Compute_Metric_Tensor ( N1(g), N2(g), N3(g), N4(g), N1(g-1), N2(g-1), N3(g-1), N4(g-1), &
                                       blkm%M(Im-g, Jm, Km), blkm%dl(Im-g, Jm, Km), blkm%vol(Im-g, Jm, Km))
          Mg(g)   = blkm % M(Im-g, Jm, Km)
          dlg(g)  = blkm % dl(Im-g, Jm, Km)
          volg(g) = blkm % vol(Im-g, Jm, Km)
        case(2)
          call Compute_Metric_Tensor ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                       blkm%M(Im+g, Jm, Km), blkm%dl(Im+g, Jm, Km), blkm%vol(Im+g, Jm, Km))
          Mg(g)   = blkm % M(Im+g, Jm, Km)
          dlg(g)  = blkm % dl(Im+g, Jm, Km)
          volg(g) = blkm % vol(Im+g, Jm, Km)
        case(3)
          call Compute_Metric_Tensor ( N1(g), N2(g), N1(g-1), N2(g-1), N5(g), N6(g), N5(g-1), N6(g-1), &
                                       blkm%M(Im, Jm-g, Km), blkm%dl(Im, Jm-g, Km), blkm%vol(Im, Jm-g, Km))
          Mg(g)   = blkm % M(Im, Jm-g, Km)
          dlg(g)  = blkm % dl(Im, Jm-g, Km)
          volg(g) = blkm % vol(Im, Jm-g, Km)
        case(4)
          call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                       blkm%M(Im, Jm+g, Km), blkm%dl(Im, Jm+g, Km), blkm%vol(Im, Jm+g, Km))
          Mg(g)   = blkm % M(Im, Jm+g, Km)
          dlg(g)  = blkm % dl(Im, Jm+g, Km)
          volg(g) = blkm % vol(Im, Jm+g, Km)
        case(5)
          call Compute_Metric_Tensor ( N1(g), N1(g-1), N3(g), N3(g-1), N5(g), N5(g-1), N7(g), N7(g-1), &
                                       blkm%M(Im, Jm, Km-g), blkm%dl(Im, Jm, Km-g), blkm%vol(Im, Jm, Km-g))
          Mg(g)   = blkm % M(Im, Jm, Km-g)
          dlg(g)  = blkm % dl(Im, Jm, Km-g)
          volg(g) = blkm % vol(Im, Jm, Km-g)
        case(6)
          call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                       blkm%M(Im, Jm, Km+g), blkm%dl(Im, Jm, Km+g), blkm%vol(Im, Jm, Km+g))
          Mg(g)   = blkm % M(Im, Jm, Km+g)
          dlg(g)  = blkm % dl(Im, Jm, Km+g)
          volg(g) = blkm % vol(Im, Jm, Km+g)
        end select

    enddo ! End loop over ghost cell           

  end subroutine BC_Connect_Metrics


  function Is_Wall(type) result (ans)
      implicit none
      integer, intent(in) :: type
      logical :: ans

      select case (type)
        case (103, 104, 301, 302, 503, 504, 505)  ! wall BCs (multi-solver coupling, heat flux, T, GSI)
          ans = .true.
        case default
          ans = .false.
      end select
    end function Is_Wall

end module MOSE_Lib_Metrics