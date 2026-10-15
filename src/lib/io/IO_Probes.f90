module MOSE_IO_Probes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Parameters_m
  use MOSE_Base_Types_m, only: MOSE_vector_3D_type
  use MOSE_Config_Types_m, only: obj_io_probes
  use MOSE_Advanced_Types_m, only: MOSE_domain_type

  implicit none
  private
  public :: Setup_Probes, Write_Probes_Data, Close_Probes, Restart_Probes, Drop_Probe_Rows, Drop_Rows_After

  type :: real_ptr
    real(R8), pointer :: p
  end type real_ptr

  type :: obj_probe
    type(MOSE_vector_3D_type)        :: location
    integer, dimension(4)            :: ilocation
    integer                          :: nvar
    character(len=llen), allocatable :: names(:)
    type(real_ptr), allocatable      :: variables(:)
    real(R8), pointer                :: P(:)
    real(R8)                         :: T, M
    real(R8)                         :: dtime
    integer                          :: ntime = 1   ! next multiple of dtime with a time-based row
    integer                          :: diter
    integer                          :: unit
  contains
    private 
    procedure, pass(self) :: Place
    procedure, pass(self) :: Update
  end type obj_probe

  integer, public :: nprobes
  type(obj_probe), allocatable, target :: probe(:)

contains

  subroutine Setup_Probes( domain, newrun )
    use IR_precision
    use strings, only: parse
    use MOSE_Mod_MPI, only: is_local_block, mpi_abort_all
    implicit none
    type(MOSE_domain_type), intent(in)  :: domain
    logical, intent(in)                 :: newrun
    ! Local
    integer             :: i, ii, error
    character(len=clen) :: string(8)

    nprobes = size(obj_io_probes)
    allocate(probe(1:nprobes))

    do i = 1, nprobes

      obj_io_probes(i)%file = 'OUTPUT/'//trim(obj_io_probes(i)%file)//'.txt'

      ! Vars name
      string = ''
      probe(i)%nvar = 0
      call parse(obj_io_probes(i)%varnames,' ',string)
      do ii = 1, size(string)
        if (string(ii) /= '') then
          probe(i)%nvar = probe(i)%nvar + 1
        endif
      enddo
      allocate(probe(i)%names(1:probe(i)%nvar))
      do ii = 1, probe(i)%nvar
        probe(i)%names(ii) = trim(string(ii))
      enddo

      ! Frequency
      probe(i)%dtime = obj_io_probes(i)%dtime
      probe(i)%diter = obj_io_probes(i)%diter

      ! Location
      probe(i)%location%c = obj_io_probes(i)%loc
      probe(i)%ilocation  = obj_io_probes(i)%iloc

      if (sum(obj_io_probes(i)%iloc)==0) call probe(i)%Place(domain)

      ! Only the rank owning this probe's block sets up pointers and opens the file
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle

      call Assign_Variables(probe(i), domain)

      ! A restart appends to the probe file, creating it if missing (e.g. probes added
      ! or files moved away); a new run replaces it. Only the rank that owns the probe
      ! gets here, so a failure aborts all ranks: an unopened unit would write to an
      ! arbitrary file.
      if (.not.newrun) then
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='UNKNOWN', &
             position='APPEND',action='WRITE',iostat=error)
      else
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='REPLACE',iostat=error)
      endif
      if (error/=0) call mpi_abort_all('cannot open probe file '//trim(obj_io_probes(i)%file))

    enddo

  end subroutine Setup_Probes


  !> Find the cell indexes related to the probe starting from space coords: the cell whose centre
  !> (the mean of its eight nodes) is the nearest, among the cells of the blocks (no ghost cells)
  subroutine Place (self, domain)
    use MOSE_Mod_MPI, only: is_local_block, mpi_size_, mpi_rank_, &
                            mpi_allreduce_min_r8, mpi_allreduce_sum_r8_array
    implicit none
    class(obj_probe)                   :: self
    type(MOSE_domain_type), intent(in) :: domain
    ! Local
    integer :: i, j, k, b
    real(8) :: d0, d, dglob, claim, winner, centre(3)
    real(8) :: packed(4)

    d0 = huge(1d0)
    self%ilocation = 0
    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        centre = 0.125d0 * ( domain%blk(b)%node(i-1,j-1,k-1)%c + domain%blk(b)%node(i,j-1,k-1)%c + &
                             domain%blk(b)%node(i-1,j  ,k-1)%c + domain%blk(b)%node(i,j  ,k-1)%c + &
                             domain%blk(b)%node(i-1,j-1,k  )%c + domain%blk(b)%node(i,j-1,k  )%c + &
                             domain%blk(b)%node(i-1,j  ,k  )%c + domain%blk(b)%node(i,j  ,k  )%c )
        d = norm2(self%location%c-centre)
        if (d<d0) then
          d0 = d
          self%ilocation = [b , i, j, k]
        endif
      enddo; enddo; enddo
    enddo

    if (mpi_size_ > 1) then
      call mpi_allreduce_min_r8(d0, dglob)
      claim = huge(1d0)
      if (d0 == dglob) claim = real(mpi_rank_, 8)
      call mpi_allreduce_min_r8(claim, winner)
      packed = 0.d0
      if (real(mpi_rank_, 8) == winner) packed = real(self%ilocation, 8)
      call mpi_allreduce_sum_r8_array(packed, 4)
      self%ilocation = nint(packed)
    endif

  end subroutine Place


  subroutine Assign_Variables(probe, domain)
    use IR_precision
    use MOSE_Global_m, only: nsc, nu, nv, nw, np, nt, nrans
    use MOSE_Mod_MPI,  only: mpi_abort_all
    implicit none
    type(obj_probe), intent(inout), target     :: probe
    type(MOSE_domain_type), intent(in), target :: domain
    ! Local
    integer :: v,s,b,i,j,k

    b = probe%ilocation(1)
    i = probe%ilocation(2)
    j = probe%ilocation(3)
    k = probe%ilocation(4)

    probe%P => domain%blk(b)%P(:,i,j,k)
    allocate(probe%variables(1:probe%nvar))
    do v = 1, probe%nvar
      nullify(probe%variables(v)%p)
    enddo

    do v = 1, probe%nvar
      do s = 1, nsc
        if (probe%names(v)=='rho('//trim(str(.true.,s))//')') probe%variables(v)%p => probe%P(s)
      enddo
      if (probe%names(v)=='u') probe%variables(v)%p => probe%P(nu)
      if (probe%names(v)=='v') probe%variables(v)%p => probe%P(nv)
      if (probe%names(v)=='w') probe%variables(v)%p => probe%P(nw)
      if (probe%names(v)=='p') probe%variables(v)%p => probe%P(np)
      if (probe%names(v)=='T') probe%variables(v)%p => probe%T
      if (probe%names(v)=='M') probe%variables(v)%p => probe%M
      ! First RANS variable, which the solution files name rhoNuTilde (Spalart-Allmaras),
      ! rhoK (two-equation models) or ru'u' (Reynolds-stress model)
      if (probe%names(v)=='mi_tilde' .and. nrans > 0) probe%variables(v)%p => probe%P(nt)
      if (.not. associated(probe%variables(v)%p)) &
        call mpi_abort_all('unknown probe variable "'//trim(probe%names(v))// &
                           '" (accepted: rho(n) u v w p T M, and mi_tilde with a RANS model)')
    enddo

  end subroutine Assign_Variables


  subroutine Update(self)
    use MOSE_Global_m,           only: nsc, nu, nv, nw, np
    use FLINT_Lib_Thermodynamic, only: f_Rtot, f_ss, EOS

    implicit none
    class(obj_probe), intent(inout) :: self
    ! Local
    real(8) :: Rtot, speed_of_sound

    Rtot = f_Rtot( self%P(1:nsc) )
    speed_of_sound = f_ss(self%P(1:nsc),self%P(np),sum(self%P(1:nsc)),Rtot)
    self%T = EOS(p=self%P(np),rho=sum(self%P(1:nsc)),R=Rtot)
    self%M = norm2(self%P(nu:nw))/speed_of_sound

  end subroutine Update


  !> \brief Time counter of the time-based probe rows at a restart.
  !> \details Write_Probes_Data writes a time-based row when time >= dtime*ntime and then counts
  !>          it, so the rows fall on the multiples of dtime, from the first (ntime starts at 1 in a
  !>          new run). A restart that goes on in time (resume) takes for ntime the first multiple
  !>          its time has not reached, computed with the same product as Write_Probes_Data so that an
  !>          exact multiple counts as written; otherwise, and in a steady run (time < 0), the rows
  !>          start again from the first one.
  subroutine Restart_Probes ( time, resume )
    implicit none
    real(R8), intent(in) :: time
    logical,  intent(in) :: resume
    ! Local
    integer :: i

    do i = 1, nprobes
      probe(i)%ntime = 1
      if (.not. resume .or. time < 0.d0) cycle
      if (time / probe(i)%dtime >= 0.5_R8 * real(huge(1), R8)) cycle
      probe(i)%ntime = max(0, int(time / probe(i)%dtime) - 1)
      do while ( time >= probe(i)%dtime * probe(i)%ntime )
        probe(i)%ntime = probe(i)%ntime + 1
      enddo
    enddo

  end subroutine Restart_Probes


  !> \brief Drop from the probe files the rows written after the solution a restart reads.
  !> \details For a restart that goes on from the counters of its solution (see Drop_Rows_After). Each
  !>          rank handles the probes it owns: it closes the file, drops the rows, and opens it again
  !>          at its end, as Setup_Probes does.
  subroutine Drop_Probe_Rows ( iteration, time )
    use MOSE_Mod_MPI, only: is_local_block, mpi_abort_all
    implicit none
    integer,  intent(in) :: iteration
    real(R8), intent(in) :: time
    ! Local
    integer :: i, error

    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      close(probe(i)%unit)
      call Drop_Rows_After ( trim(obj_io_probes(i)%file), 1 + probe(i)%nvar, iteration, time )
      open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='UNKNOWN', &
           position='APPEND',action='WRITE',iostat=error)
      if (error/=0) call mpi_abort_all('cannot open probe file '//trim(obj_io_probes(i)%file))
    enddo

  end subroutine Drop_Probe_Rows


  !> \brief Drop the rows of a history file written after the solution a restart reads.
  !> \details residual-history.dat and the probe files get a row at the iterations, or at the times,
  !>          that call for one. A run that went on past the solution a restart reads (it stopped
  !>          without writing its final solution, and the restart reads an earlier numbered one) left
  !>          in them rows that the restart writes again. The rows kept are those before the first one
  !>          whose first value is an iteration greater than iteration, or a time later than time; that
  !>          row and all the rows after it are dropped, so the file goes on as that of the run in one
  !>          go. A row holds ncol values; written list-directed it may take more than one line, and a
  !>          row starts on a new line: the values are counted to find where each row starts. A time is
  !>          compared as printed, to 1e-12 relative: the row of the step of the solution may print its
  !>          time rounded up in the last digit. A last row with fewer than ncol values (a file cut in
  !>          the middle of a row) is dropped too. A line that starts with '#' (a header) is kept and
  !>          holds no values. The file is truncated and closed, and ends with a newline
  !>          (End_With_Newline); a missing file is left alone.
  subroutine Drop_Rows_After ( file, ncol, iteration, time )
    implicit none
    character(len=*), intent(in) :: file
    integer,          intent(in) :: ncol, iteration
    real(R8),         intent(in) :: time
    ! Local
    character(len=32768) :: line
    integer              :: unit, ios, nkeep, nrow, nval, k, k1, k2, it
    real(R8)             :: t
    logical              :: later

    open(newunit=unit,file=file,status='OLD',action='READWRITE',iostat=ios)
    if (ios /= 0) return
    nkeep = 0
    nrow  = 0
    nval  = 0
    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      ! A header line is kept and holds no values
      if (len_trim(line) > 0) then
        if (line(verify(line, ' '):verify(line, ' ')) == '#') then
          nkeep = nkeep + 1
          cycle
        endif
      endif
      if (mod(nval, ncol) == 0 .and. len_trim(line) > 0) then
        ! The line starts a row: its first value, line(k1:k2), is an iteration or a time
        k1 = verify(line, ' ')
        k2 = k1 + index(line(k1:), ' ') - 2
        later = .false.
        if (scan(line(k1:k2), '.eEdD') > 0) then
          read(line(k1:k2),*,iostat=ios) t
          if (ios == 0) later = t > time + abs(time) * 1.0e-12_R8
        else
          read(line(k1:k2),*,iostat=ios) it
          if (ios == 0) later = it > iteration
        endif
        if (later) exit
        ! Lines kept before the row that starts here
        nrow = nkeep
      endif
      ! Values on the line
      do k = 1, len_trim(line)
        if (line(k:k) /= ' ' .and. (k == 1 .or. line(max(k-1,1):max(k-1,1)) == ' ')) nval = nval + 1
      enddo
      nkeep = nkeep + 1
    enddo
    ! A last row with fewer than ncol values is cut: drop it from its first line
    if (mod(nval, ncol) /= 0) nkeep = nrow
    rewind(unit)
    do k = 1, nkeep
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
    enddo
    endfile(unit)
    close(unit)
    call End_With_Newline ( file )

  end subroutine Drop_Rows_After


  !> \brief End a file with a newline if its last line has none.
  !> \details A history file cut by hand may end without a newline; the rows that a restart then
  !>          writes in append mode would go on that line. Read and written as a stream of bytes.
  subroutine End_With_Newline ( file )
    use iso_fortran_env, only: int64
    implicit none
    character(len=*), intent(in) :: file
    ! Local
    integer        :: unit, ios
    integer(int64) :: nbytes
    character      :: last

    open(newunit=unit,file=file,access='STREAM',form='UNFORMATTED',status='OLD',action='READWRITE',iostat=ios)
    if (ios /= 0) return
    inquire(unit=unit,size=nbytes)
    if (nbytes > 0) then
      read(unit,pos=nbytes,iostat=ios) last
      if (ios == 0 .and. last /= new_line('a')) write(unit,pos=nbytes+1,iostat=ios) new_line('a')
    endif
    close(unit)

  end subroutine End_With_Newline


  !> \brief Write the probe rows due at this iteration.
  !> \details Each rank writes the rows of the probes it owns, and passes each row to the file system
  !>          when it writes it (flush): a run that stops before its end (killed, stopped at the time
  !>          limit of a queue, aborted) leaves in the files every row it wrote, as a restart from its
  !>          last solution expects (Drop_Rows_After); the rows held in the buffer of the unit would be
  !>          lost with the process.
  subroutine Write_Probes_Data( iter, time )
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    integer, intent(in) :: iter
    real(8), intent(in) :: time
    ! Local
    integer :: i, v

    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      if (mod(iter,probe(i)%diter)==0) then
        call probe(i)%update
        write(probe(i)%unit,*) iter, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
        flush(probe(i)%unit)
      elseif (time >= probe(i)%dtime*probe(i)%ntime) then
        call probe(i)%update
        write(probe(i)%unit,*) time, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
        flush(probe(i)%unit)
        probe(i)%ntime = probe(i)%ntime+1
      endif
    enddo

  end subroutine Write_Probes_Data


  !> \brief Close the probe files at the end of the run.
  !> \details Each rank closes the files of the probes it owns, which Setup_Probes or Drop_Probe_Rows
  !>          opened, before the ranks finalize MPI: a unit left open is closed by the Fortran run-time
  !>          library only when the process ends, after MPI_FINALIZE.
  subroutine Close_Probes ()
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    ! Local
    integer :: i

    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      close(probe(i)%unit)
    enddo

  end subroutine Close_Probes


end module MOSE_IO_Probes