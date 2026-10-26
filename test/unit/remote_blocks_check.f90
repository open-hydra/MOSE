!> Unit check of the blocks whose state an MPI rank holds: blocks_needing_full_arrays and deallocate_remote_computation_data
!> (src/lib/driver/Mod_Allocate_Data.f90).
!>
!> A rank holds the state arrays of the blocks it owns and of the donors of the chimera faces (102) of the receiver blocks it
!> owns: the chimera exchange receives the donor cells of the receivers of the rank only, into the P of the remote donor
!> (build_chimera_schedule and chimera_unpack, src/lib/parallel/Mod_GhostExchange.f90). The donors of a face of type 104 are
!> cells of the other solver, not blocks of MOSE. Four blocks on three ranks, blocks 1 to 4 on ranks 0, 1, 2 and 2:
!>   a face 102 of block 2 (rank 1), its donors in block 3 (rank 2);
!>   a face 102 of block 4 (rank 2), its donors in block 2 (rank 1);
!>   a face 104 of block 1 (rank 0), its donors numbered 4 (a block of the other solver with the number of a block of MOSE).
!> Each rank is set as MOSE_Mod_MPI sets it in a run on three ranks (block_owner, mpi_rank_, mpi_size_, mpi_is_root); no MPI
!> call is made. For each rank:
!>   1. blocks_needing_full_arrays: the blocks whose full arrays the set-up allocates on one grid level;
!>   2. every block allocated in full, as the set-up allocates them with multigrid, then deallocate_remote_computation_data:
!>      the blocks that keep P.
!> Both give block 1 on rank 0, blocks 2 and 3 on rank 1, blocks 2, 3 and 4 on rank 2. Rank 0 holds no chimera donor: it owns
!> no receiver, and the donors of its face 104 are not blocks of MOSE.
program remote_blocks_check
  use MOSE_Advanced_Types_m,  only: MOSE_domain_type, MOSE_bc_type
  use MOSE_Global_m,          only: nprim, model
  use MOSE_Config_Types_m,    only: obj_irs
  use MOSE_Mod_MPI,           only: block_owner, mpi_rank_, mpi_size_, mpi_is_root
  use MOSE_Mod_Allocate_Data, only: Allocate_Block, blocks_needing_full_arrays, deallocate_remote_computation_data
  implicit none
  integer, parameter :: nb = 4, nranks = 3
  ! held(:, r): the blocks whose state rank r holds
  logical, parameter :: held(nb, 0:nranks-1) = reshape( [ .true.,  .false., .false., .false., &
                                                          .false., .true.,  .true.,  .false., &
                                                          .false., .true.,  .true.,  .true. ], [nb, nranks] )
  integer :: r, nfail

  nprim = 5                   ! one species: density, three velocity components, pressure
  model = 0                   ! Euler: no viscous arrays
  obj_irs%enabled = .false.   ! no residual smoothing arrays
  mpi_size_ = nranks
  allocate ( block_owner(nb) )
  block_owner = [0, 1, 2, 2]

  nfail = 0
  do r = 0, nranks - 1
    mpi_rank_ = r
    mpi_is_root = ( r == 0 )
    call check_rank ( r )
  end do

  if ( nfail == 0 ) then
    write(*,'(A)') 'RemoteBlocks: all checks passed'
  else
    write(*,'(A,I0,2A)') 'RemoteBlocks: FAILED (', nfail, ' checks): the blocks whose state a rank holds are not its own', &
      ' and the chimera donors of the receivers it owns'
    error stop 1
  end if

contains

  subroutine check_rank ( r )
    integer, intent(in)    :: r
    type(MOSE_domain_type) :: domain
    logical, allocatable   :: keep(:)
    logical                :: kept(nb)
    integer                :: b

    call build_domain ( domain )
    keep = blocks_needing_full_arrays ( domain )
    call report ( r, 'full arrays at set-up', keep, held(:, r) )
    call deallocate_remote_computation_data ( domain )
    do b = 1, nb
      kept(b) = allocated(domain%blk(b)%P)
    end do
    call report ( r, 'P kept after set-up  ', kept, held(:, r) )
  end subroutine check_rank


  subroutine build_domain ( domain )
    type(MOSE_domain_type), intent(out) :: domain
    integer :: b

    domain%nb = nb
    allocate ( domain%blk(nb) )
    do b = 1, nb
      domain%blk(b)%dim = [2, 2, 1]
      call Allocate_Block ( domain%blk(b), [2, 2, 1] )
    end do
    domain%nbound = 3
    allocate ( domain%bc(domain%nbound) )
    call set_face ( domain%bc(1), 102, 2, 3 )
    call set_face ( domain%bc(2), 102, 4, 2 )
    call set_face ( domain%bc(3), 104, 1, 4 )
  end subroutine build_domain


  !> A boundary record of face 1 of block b, cell (1,1,1), with one donor cell in each of its two lists, both in block donor
  subroutine set_face ( bc, bctype, b, donor )
    type(MOSE_bc_type), intent(inout) :: bc
    integer, intent(in)               :: bctype, b, donor

    bc%type = bctype
    bc%b = b
    bc%f = 1
    bc%i = 1
    bc%j = 1
    bc%k = 1
    bc%ni = [1, 1]
    allocate ( bc%donorID(2, 4) )
    bc%donorID(:, 1) = donor
    bc%donorID(:, 2) = [1, 2]
    bc%donorID(:, 3) = 1
    bc%donorID(:, 4) = 1
  end subroutine set_face


  subroutine report ( r, what, got, want )
    integer, intent(in)          :: r
    character(len=*), intent(in) :: what
    logical, intent(in)          :: got(nb), want(nb)
    character(len=32) :: bad
    integer :: b

    bad = ''
    do b = 1, nb
      if ( got(b) .neqv. want(b) ) then
        nfail = nfail + 1
        write(bad(len_trim(bad)+2:), '(I0)') b
      end if
    end do
    if ( len_trim(bad) == 0 ) then
      write(*,'(A,I0,2A,4(1X,L1),A,4(1X,L1),A)') ' rank ', r, ': ', what, got, ', expected', want, ': ok'
    else
      write(*,'(A,I0,2A,4(1X,L1),A,4(1X,L1),2A)') ' rank ', r, ': ', what, got, ', expected', want, &
        ': FAIL at block', trim(bad)
    end if
  end subroutine report

end program remote_blocks_check
