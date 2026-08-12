module MOSE_Mod_Allocate_Data
  use, intrinsic :: iso_fortran_env, only : iostat_end

  implicit none
  private
  public :: Setup_Data_Structure, Setup_Data_Dimensions, Setup_Data_Fields
  public :: Allocate_Block, deallocate_remote_computation_data
  public :: blocks_needing_full_arrays

contains

  subroutine Setup_Data_Structure ( domain, IOfield )
    !! Dimensions plus a full allocation of every block, for the paths that
    !! cannot partition first.
    use MOSE_Advanced_Types_m
    use Lib_ORION_data

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield

    call Setup_Data_Dimensions ( domain, IOfield )
    call Setup_Data_Fields ( domain, IOfield )

  end subroutine Setup_Data_Structure


  subroutine Setup_Data_Dimensions ( domain, IOfield )
    !! Everything known before any large array is allocated: variable-index
    !! globals, block count and per-block ijk dimensions.  This is what
    !! `partition_blocks` needs, so ownership can be decided before the memory
    !! is committed.
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use Lib_ORION_data

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield
    ! Local
    integer :: b, nblocks

    nu = nsc + 1
    nv = nu + 1
    nw = nv + 1
    np = nw + 1
    nprim = np + nsoot + npass + nrans

    if (nsoot > 0) then
      nc = np + 1
    end if

    if (npass>0) then
      nps = np + nsoot + 1
    endif
    
    if (nrans > 0) then
      nt = np + nsoot + npass + 1
    else
      nt = nprim
    endif
    nres = nprim - nsc + 1

    ! Domain is the MOSE-alias of IOfield
    nblocks = size ( IOfield%block )
    allocate( domain%blk( 1:nblocks ) )
    domain%nb = nblocks

    ! Define the ijk dimensions of each block
    do b = 1, nblocks
      domain%blk(b)%dim(1) = IOfield%block(b)%Ni
      domain%blk(b)%dim(2) = IOfield%block(b)%Nj
      domain%blk(b)%dim(3) = IOfield%block(b)%Nk
    enddo

    ! Check if number of simulation variables in orion-field matches MOSE expectation
    if ( size(IOfield%block(1)%vars, 1) < nprim ) then
      write(*,'(A)')         '[ERROR] Number of variables in IOfield does not match MOSE expectation.'
      write(*,'(A,I0,A,I0)') '        Expected: ', nprim, ', Found: ', size(IOfield%block(1)%vars, 1)
      stop
    end if

  end subroutine Setup_Data_Dimensions


  subroutine Setup_Data_Fields ( domain, IOfield, keep )
    !! Allocate the block arrays and import the initial condition.
    !!
    !! `keep(b)` selects the blocks that get the full set of arrays (~800
    !! B/cell); the others get geometry only, `node` and `yn` (~36 B/cell).
    !! Without `keep`, every block gets everything.
    !!
    !! The excluded blocks still need geometry because set-up reads the mesh
    !! of blocks this rank does not own: `Check_Mesh_Type` reads `blk(1)%node`
    !! everywhere, `BC_Connect_Metrics` builds ghost metrics from a possibly
    !! remote donor's nodes, and `Compute_Yn` minimises over every wall face in
    !! the domain.  `deallocate_remote_computation_data` releases the geometry
    !! too, once metrics are built.
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_io, obj_chemistry
    use MOSE_Global_m
    use Lib_ORION_data

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield
    logical, intent(in), optional         :: keep(:)
    ! Local
    integer :: b, d
    logical :: full

    do b = 1, domain%nb

      full = .true.
      if ( present(keep) ) full = keep(b)

      ! Allocate domain block
      call Allocate_Block ( domain%blk(b), domain%blk(b)%dim, geometry_only = .not. full )

      ! Import domain-block nodes from orion-field
      do d = 1, 3
        domain%blk(b)%node(0:IOfield%block(b)%Ni,0:IOfield%block(b)%Nj,0:IOfield%block(b)%Nk)%c(d) &
        = IOfield%block(b)%mesh(d,0:,0:,0:)
      enddo

      ! Import domain-block primitives from orion-field variables
      if ( full ) then
        do d = 1, nprim
          domain%blk(b)%P(d,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk) &
          = IOfield%block(b)%vars(d,1:,1:,1:)
        enddo
      end if

      ! Flag chemistry-inibited blocks
      if (allocated(obj_chemistry%no_chem_list)) then
        if (any(obj_chemistry%no_chem_list == b)) then
          domain%blk(b)%no_chem = .true.
        end if
      end if

    enddo

    domain%time = obj_io%IOtime

  end subroutine Setup_Data_Fields


  function blocks_needing_full_arrays ( domain ) result ( keep )
    !! Blocks whose state arrays this rank has to hold: the ones it owns, plus
    !! the chimera donors and manifold sources whose `P` it reads directly.
    use MOSE_Advanced_Types_m
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    logical, allocatable               :: keep(:)
    ! Local
    integer :: b, c, i

    allocate(keep(domain%nb))
    do b = 1, domain%nb
      keep(b) = is_local_block(b)
    end do

    do i = 1, domain%nbound
      select case (domain%bc(i)%type)
        case (102) ! chimera
          if (allocated(domain%bc(i)%donorID)) then
            do c = 1, size(domain%bc(i)%donorID, 1)
              b = domain%bc(i)%donorID(c, 1)
              if (b >= 1 .and. b <= domain%nb) keep(b) = .true.
            end do
          end if
        case (501) ! manifold
          b = domain%bc(i)%bs
          if (b >= 1 .and. b <= domain%nb) keep(b) = .true.
      end select
    end do

  end function blocks_needing_full_arrays


  subroutine Allocate_Block( blk, nijk, geometry_only )
    !! `geometry_only` allocates the mesh arrays `node` and `yn` and skips the
    !! state, metric and residual ones -- see Setup_Data_Fields.
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_irs, obj_rans
    use MOSE_Global_m
    implicit none
    integer, intent(in)                  :: nijk(3)
    type(MOSE_block_type), intent(inout) :: blk
    logical, intent(in), optional        :: geometry_only
    ! Local
    integer :: ni, nj, nk
    logical :: geom

    geom = .false.
    if ( present(geometry_only) ) geom = geometry_only

    ni = nijk(1) ; nj = nijk(2) ; nk = nijk(3)

    ! Mesh -- always present, on every rank
    allocate( blk % node ( 0:ni, 0:nj, 0:nk ) )
    allocate ( blk % yn ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )

    if ( geom ) then
      call First_Touch_Block( blk )
      return
    end if

    ! Metrics
    allocate( blk % M    ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % dl   ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % vol  ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % dir(1) % f (0:ni, nj, nk) )
    allocate( blk % dir(2) % f (ni, 0:nj, nk) )
    allocate( blk % dir(3) % f (ni, nj, 0:nk) )

    ! Preconditioning reference velocity
    allocate( blk % Ur ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )

    ! Prim and Residuals
    allocate( blk % P (nprim, 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % PO, blk % R, mold = blk % P )

    ! Dt cell center with no ghost cells
    allocate( blk % dtlocal ( 1:ni, 1:nj, 1:nk) )

    ! Shock-detector flag with no ghost cells
    allocate( blk % beta ( 1:ni, 1:nj, 1:nk) )

    ! Velocity gradient
    if ( model > 1 ) then
      allocate( blk % vel_gradient ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    end if

    ! Rotation/curvature correction terms
    if ( model > 1 ) then
      allocate( blk % rc_term1(ni, nj, nk) )
      allocate( blk % rc_term2(ni, nj, nk) )
    end if

    ! Temp storage for residuals in IRS
    if ( obj_irs%enabled ) then
      allocate( blk % RS1, blk % RS2, mold = blk % R )
    end if

    call First_Touch_Block( blk )

  end subroutine Allocate_Block


  subroutine First_Touch_Block( blk )
    !! Place this block's pages on the NUMA node of the thread that will work
    !! on them: a page binds to a node when it is first written, so without
    !! this the whole block lands wherever the master thread touched it.
    !!
    !! The (k,j) split must match the one the solver uses, or the pages end up
    !! on the wrong node anyway.  Zeroing rather than merely touching also
    !! defines arrays that were previously undefined until first written.
    use MOSE_Advanced_Types_m
    implicit none
    type(MOSE_block_type), intent(inout) :: blk
    ! Local
    integer :: j, k, p, q

    !$omp parallel default(shared) private(j,k,p,q)

    ! Metrics
    !$omp do collapse(2) schedule(static)
    do k = lbound(blk%node,3), ubound(blk%node,3)
      do j = lbound(blk%node,2), ubound(blk%node,2)
        blk % node(:,j,k) % c(1) = 0.0d0
        blk % node(:,j,k) % c(2) = 0.0d0
        blk % node(:,j,k) % c(3) = 0.0d0
      end do
    end do

    !$omp do collapse(2) schedule(static)
    do k = lbound(blk%yn,3), ubound(blk%yn,3)
      do j = lbound(blk%yn,2), ubound(blk%yn,2)
        blk % yn (:,j,k) = 0.0d0
      end do
    end do

    ! Absent on geometry-only blocks.  The conditions are uniform across the
    ! team, so guarding whole `!$omp do` regions is safe.
    if ( allocated(blk%vol) ) then
      !$omp do collapse(2) schedule(static)
      do k = lbound(blk%vol,3), ubound(blk%vol,3)
        do j = lbound(blk%vol,2), ubound(blk%vol,2)
          blk % vol(:,j,k) = 0.0d0
          blk % Ur (:,j,k) = 0.0d0
          do q = 1, 3
            blk % dl(:,j,k) % c(q) = 0.0d0
            do p = 1, 3
              blk % M(:,j,k) % c(p,q) = 0.0d0
            end do
          end do
        end do
      end do
    end if

    ! Primitive state and residuals -- the bulk of the footprint.
    if ( allocated(blk%P) ) then
      !$omp do collapse(2) schedule(static)
      do k = lbound(blk%P,4), ubound(blk%P,4)
        do j = lbound(blk%P,3), ubound(blk%P,3)
          blk % P (:,:,j,k) = 0.0d0
          blk % PO(:,:,j,k) = 0.0d0
          blk % R (:,:,j,k) = 0.0d0
        end do
      end do
    end if

    if ( allocated(blk%RS1) ) then
      !$omp do collapse(2) schedule(static)
      do k = lbound(blk%RS1,4), ubound(blk%RS1,4)
        do j = lbound(blk%RS1,3), ubound(blk%RS1,3)
          blk % RS1(:,:,j,k) = 0.0d0
          blk % RS2(:,:,j,k) = 0.0d0
        end do
      end do
    end if

    if ( allocated(blk%dtlocal) ) then
      !$omp do collapse(2) schedule(static)
      do k = lbound(blk%dtlocal,3), ubound(blk%dtlocal,3)
        do j = lbound(blk%dtlocal,2), ubound(blk%dtlocal,2)
          blk % dtlocal(:,j,k) = 0.0d0
          blk % beta   (:,j,k) = 0.0d0
        end do
      end do
    end if

    if ( allocated(blk%vel_gradient) ) then
      !$omp do collapse(2) schedule(static)
      do k = lbound(blk%vel_gradient,3), ubound(blk%vel_gradient,3)
        do j = lbound(blk%vel_gradient,2), ubound(blk%vel_gradient,2)
          do q = 1, 3
            do p = 1, 3
              blk % vel_gradient(:,j,k) % c(p,q) = 0.0d0
            end do
          end do
        end do
      end do
    end if

    if ( allocated(blk%rc_term1) ) then
      !$omp do collapse(2) schedule(static)
      do k = lbound(blk%rc_term1,3), ubound(blk%rc_term1,3)
        do j = lbound(blk%rc_term1,2), ubound(blk%rc_term1,2)
          blk % rc_term1(:,j,k) = 0.0d0
          blk % rc_term2(:,j,k) = 0.0d0
        end do
      end do
    end if

    ! Face metrics, one loop per direction: different extents, so they cannot
    ! share a (k,j) iteration space.
    if ( allocated(blk%dir(1)%f) ) then
      !$omp do collapse(2) schedule(static)
      do k = lbound(blk%dir(1)%f,3), ubound(blk%dir(1)%f,3)
        do j = lbound(blk%dir(1)%f,2), ubound(blk%dir(1)%f,2)
          blk % dir(1) % f(:,j,k) % A = 0.0d0
          do q = 1, 3
            blk % dir(1) % f(:,j,k) % N(q) = 0.0d0
          end do
        end do
      end do

      !$omp do collapse(2) schedule(static)
      do k = lbound(blk%dir(2)%f,3), ubound(blk%dir(2)%f,3)
        do j = lbound(blk%dir(2)%f,2), ubound(blk%dir(2)%f,2)
          blk % dir(2) % f(:,j,k) % A = 0.0d0
          do q = 1, 3
            blk % dir(2) % f(:,j,k) % N(q) = 0.0d0
          end do
        end do
      end do

      !$omp do collapse(2) schedule(static)
      do k = lbound(blk%dir(3)%f,3), ubound(blk%dir(3)%f,3)
        do j = lbound(blk%dir(3)%f,2), ubound(blk%dir(3)%f,2)
          blk % dir(3) % f(:,j,k) % A = 0.0d0
          do q = 1, 3
            blk % dir(3) % f(:,j,k) % N(q) = 0.0d0
          end do
        end do
      end do
    end if

    !$omp end parallel

  end subroutine First_Touch_Block


  subroutine deallocate_remote_computation_data(domain)
    use MOSE_Advanced_Types_m
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root, mpi_size_, block_owner

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, d, i, c
    logical, allocatable :: needs_remote_P(:)

    call check_donors_are_colocated(domain)

    ! Build mask of remote blocks whose P (and dir) must be kept:
    !  - chimera (102): donorID(:,1) can reference remote blocks
    !  - manifold (501): bc%bs can be remote
    allocate(needs_remote_P(domain%nb))
    needs_remote_P = .false.
    do i = 1, domain%nbound
      select case (domain%bc(i)%type)
        case (102) ! chimera
          if (allocated(domain%bc(i)%donorID)) then
            do c = 1, size(domain%bc(i)%donorID, 1)
              b = domain%bc(i)%donorID(c, 1)
              if (.not. is_local_block(b)) needs_remote_P(b) = .true.
            end do
          end if
        case (501) ! manifold
          b = domain%bc(i)%bs
          if (b > 0 .and. .not. is_local_block(b)) needs_remote_P(b) = .true.
      end select
    end do

    do b = 1, domain%nb
      if (is_local_block(b)) cycle

      ! Computation arrays — free on all ranks
      if (allocated(domain%blk(b)%PO))           deallocate(domain%blk(b)%PO)
      if (allocated(domain%blk(b)%R))            deallocate(domain%blk(b)%R)
      if (allocated(domain%blk(b)%RS1))          deallocate(domain%blk(b)%RS1)
      if (allocated(domain%blk(b)%RS2))          deallocate(domain%blk(b)%RS2)
      if (allocated(domain%blk(b)%dtlocal))      deallocate(domain%blk(b)%dtlocal)
      if (allocated(domain%blk(b)%beta))         deallocate(domain%blk(b)%beta)
      if (allocated(domain%blk(b)%vel_gradient)) deallocate(domain%blk(b)%vel_gradient)
      if (allocated(domain%blk(b)%rc_term1))     deallocate(domain%blk(b)%rc_term1)
      if (allocated(domain%blk(b)%rc_term2))     deallocate(domain%blk(b)%rc_term2)

      ! P — free unless this remote block is a chimera donor or manifold source
      if (.not. needs_remote_P(b)) then
        if (allocated(domain%blk(b)%P)) deallocate(domain%blk(b)%P)
      end if

      ! Metrics — free on non-root ranks only (root needs them for wall I/O)
      ! Keep dir on blocks needed for manifold (BC_Manifold reads blk(Bs)%dir%f%A)
      if (.not. mpi_is_root) then
        if (allocated(domain%blk(b)%node)) deallocate(domain%blk(b)%node)
        if (allocated(domain%blk(b)%M))    deallocate(domain%blk(b)%M)
        if (allocated(domain%blk(b)%dl))   deallocate(domain%blk(b)%dl)
        if (allocated(domain%blk(b)%vol))  deallocate(domain%blk(b)%vol)
        if (allocated(domain%blk(b)%yn))   deallocate(domain%blk(b)%yn)
        if (.not. needs_remote_P(b)) then
          do d = 1, 3
            if (allocated(domain%blk(b)%dir(d)%f)) deallocate(domain%blk(b)%dir(d)%f)
          end do
        end if
      end if
    end do

    deallocate(needs_remote_P)

  end subroutine deallocate_remote_computation_data


  subroutine check_donors_are_colocated(domain)
    use MOSE_Advanced_Types_m
    use MOSE_Mod_MPI, only: mpi_size_, block_owner, mpi_abort_all

    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: i, bm, bs

    if (mpi_size_ <= 1) return
    if (.not. allocated(block_owner)) return

    do i = 1, domain%nbound
      bm = domain%bc(i)%b
      select case (domain%bc(i)%type)
        case (501)
          bs = domain%bc(i)%bs
          call abort_if_split(bm, bs, 'manifold (501)')
      end select
    end do

  end subroutine check_donors_are_colocated


  subroutine abort_if_split(bm, bs, what)
    use MOSE_Mod_MPI, only: block_owner, mpi_abort_all

    implicit none
    integer, intent(in) :: bm, bs
    character(len=*), intent(in) :: what
    character(len=512) :: msg

    if (bs <= 0 .or. bs > size(block_owner)) return
    if (bm <= 0 .or. bm > size(block_owner)) return
    if (block_owner(bs) == block_owner(bm)) return

    write(msg,'(A,A,I0,A,I0,A,I0,A,I0,A)') trim(what), &
      ': donor block ', bs, ' (rank ', block_owner(bs), ') is split from block ', bm, &
      ' (rank ', block_owner(bm), '). Donor P is not refreshed across ranks; co-locate the blocks or use fewer ranks.'
    call mpi_abort_all(trim(msg))

  end subroutine abort_if_split


end module MOSE_Mod_Allocate_Data