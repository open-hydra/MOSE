module MOSE_IO_BC_Q2D
  use iso_fortran_env, only: R8 => real64

  implicit none

  private
  public :: Setup_Q2D_BC_Data

contains

  ! ---------------------------------------------------------------------------
  ! Setup Q2D mapped boundary condition data (BC type 410).
  ! Reads the ATLAS Tecplot file(s) specified in bc.txt (one per face),
  ! parses all solution times, maps variables by name to MOSE primitives,
  ! and populates bc(i)%q2d_map for each type 410 boundary cell.
  ! ---------------------------------------------------------------------------
  subroutine Setup_Q2D_BC_Data ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use MOSE_Mod_MPI, only: mpi_stop_all
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, nq2d_count, nfiles, ifile
    character(len=256), allocatable :: file_list(:)
    character(len=128) :: msg
    logical :: found

    ! Count type 410 BCs
    nq2d_count = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type == 410) nq2d_count = nq2d_count + 1
    enddo
    if (nq2d_count == 0) return

    ! Collect unique filenames from bc(:)%q2d_file
    allocate(file_list(nq2d_count))
    nfiles = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 410) cycle
      if (len_trim(domain%bc(i)%q2d_file) == 0) then
        write(msg,'(A,I0,A)') '[ERROR] BC type 410 entry ', i, ' has no Q2D file specified in bc.txt'
        call mpi_stop_all(trim(msg))
      endif
      found = .false.
      do ifile = 1, nfiles
        if (trim(domain%bc(i)%q2d_file) == trim(file_list(ifile))) then
          found = .true.; exit
        endif
      enddo
      if (.not. found) then
        nfiles = nfiles + 1
        file_list(nfiles) = domain%bc(i)%q2d_file
      endif
    enddo

    ! Process each unique Q2D file
    do ifile = 1, nfiles
      write(*,'(A,A)') '  Q2D file: ', trim(file_list(ifile))
      call Read_Q2D_Tecplot ( domain, trim(file_list(ifile)) )
    enddo

    deallocate(file_list)

  end subroutine Setup_Q2D_BC_Data


  ! ---------------------------------------------------------------------------
  ! Parse the ATLAS Tecplot ASCII file and distribute data to type 410 BCs.
  !
  ! Two-pass reading:
  !   Pass 1 — scan VARIABLES line + all ZONE headers (no data read)
  !   Pass 2 — rewind, read zone data and distribute to matching BCs
  ! The time of a zone, and every value that goes to a type 410 BC, must be
  ! finite: one that is not (NaN or infinity) stops the run at set-up, with the
  ! file, the zone and the line of its header (Not_Finite reads the bits:
  ! every comparison with a NaN is false, and a build that traps the invalid
  ! operations stops on it). The values that the reader skips (coordinates,
  ! variables that map to no primitive, cells of no 410 face) are not used.
  ! A file that cannot be read stops the run at set-up as well, with the file
  ! and, for a zone, the zone and the line of its header: a file that cannot
  ! be opened, with no VARIABLES line or no zone, a value of I, J,
  ! SOLUTIONTIME or of the block number of the title that does not read, zone
  ! data that end before the values that I and J give or hold a value that
  ! is not a number, a line after the data of a zone that is not the header
  ! of the next zone. A read that fails leaves its values undefined, or those
  ! of the zone before (the buffers are reused), and shifts the zones after
  ! it. The distinct times of the zones must increase down the file (the
  ! interpolation in time takes them in that order): a zone whose time no zone
  ! above it has, less than the last such time, stops the run too. Every rank
  ! reads the whole file (Setup_BC runs before the blocks are given to the
  ! ranks), so every rank stops on the same error (mpi_stop_all).
  ! ---------------------------------------------------------------------------
  subroutine Read_Q2D_Tecplot ( domain, filename )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use MOSE_Finite_m, only: Not_Finite
    use MOSE_Mod_MPI,  only: mpi_stop_all
    use strings, only: lowercase
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    character(len=*), intent(in) :: filename
    ! Local
    integer :: u, ios, iz, nvar_file
    integer :: I_nodes, J_nodes, n1, n2, n_nodes, I_max, J_max
    integer :: b_zone, t_zone, iv, ic, i1, i2, ibc, it
    integer :: ncount
    logical :: names_ok
    character(len=32768) :: line
    character(len=:), allocatable :: header
    character(len=32) :: tmp_names(100)
    integer, allocatable :: varmap(:)
    real(R8), allocatable :: node_buf(:), cell_buf(:,:,:)
    integer :: nvar_total
    integer :: nzones_total
    integer, allocatable :: zone_block(:), zone_tidx(:)
    real(R8), allocatable :: zone_time(:), unique_times_tmp(:)
    integer, allocatable :: zone_I(:), zone_J(:)
    integer :: ntimes
    real(R8), allocatable :: unique_times(:)
    integer :: bc_i1, bc_i2
    logical :: found
    integer :: nl                          ! lines of the file read in the pass over the zone headers
    integer, allocatable :: zone_line(:)   ! line of the header of each zone
    character(len=512) :: msg
    character(len=256) :: iomsg
    character(len=64)  :: bad              ! the value of a zone header that does not read (parse_zone_header) ...
    character(len=100) :: badval           ! ... and its text
    character(len=16)  :: tnew, told       ! two times of zones, written for a message

    ! --- Open file ---
    open(newunit=u, file=filename, status='old', iostat=ios, iomsg=iomsg, action='read')
    if (ios /= 0) call mpi_stop_all('[ERROR] Q2D file '//trim(filename)//' cannot be opened: '//trim(iomsg))

    ! --- Pass 1: VARIABLES line + all ZONE headers ---

    ! Read the VARIABLES list: from the VARIABLES line up to the first ZONE record
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) call mpi_stop_all('[ERROR] Q2D file '//trim(filename)//': could not find the VARIABLES line')
      if (index(lowercase(line), 'variables') > 0) exit
    enddo
    header = trim(line(index(line, '=')+1:))
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (Is_Zone_Record(line)) then
        backspace(u)
        exit
      endif
      header = header//' '//trim(line)
    enddo
    call Parse_Variable_Names(header, tmp_names, ncount, names_ok)
    if (.not. names_ok) then
      write(msg,'(3A,I0,A)') '[ERROR] Q2D file ', trim(filename), ': unreadable VARIABLES list (unclosed quote or more than ', &
        size(tmp_names), ' names)'
      call mpi_stop_all(trim(msg))
    endif
    nvar_total = ncount
    nvar_file = nvar_total - 3  ! subtract x, y, z

    if (nvar_file < 1) call mpi_stop_all('[ERROR] Q2D file '//trim(filename)//' has no solution variables (only coordinates)')

    ! Build variable name mapping
    allocate(varmap(nvar_file))
    call build_varmap(tmp_names(4:nvar_total), nvar_file, varmap)

    ! Continue scanning from current position (no rewind): count zones
    nzones_total = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (Is_Zone_Record(line)) nzones_total = nzones_total + 1
    enddo

    if (nzones_total == 0) call mpi_stop_all('[ERROR] Q2D file '//trim(filename)//' has no zones')

    allocate(zone_block(nzones_total), zone_tidx(nzones_total))
    allocate(zone_time(nzones_total), zone_I(nzones_total), zone_J(nzones_total), zone_line(nzones_total))
    allocate(unique_times_tmp(nzones_total))

    ! Rewind once: parse zone headers
    rewind(u)
    iz = 0
    nl = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      nl = nl + 1
      if (Is_Zone_Record(line)) then
        iz = iz + 1
        call parse_zone_header(line, zone_block(iz), &
                               zone_I(iz), zone_J(iz), zone_time(iz), bad, badval)
        zone_line(iz) = nl
        if (len_trim(bad) > 0) then
          write(msg,'(3A,I0,A,I0,4A)') '[ERROR] Q2D file ', trim(filename), ', line ', nl, ' (zone ', iz, &
            '): ', trim(bad), ' cannot be read: ', trim(badval)
          call mpi_stop_all(trim(msg))
        endif
        if (Not_Finite(zone_time(iz))) then
          write(msg,'(3A,I0,A,I0,A)') '[ERROR] Q2D file ', trim(filename), ', line ', nl, ' (zone ', iz, &
            '): SOLUTIONTIME is not finite (NaN or infinity)'
          call mpi_stop_all(trim(msg))
        endif
      endif
    enddo

    ! Extract unique times and build zone → time index mapping (single pass)
    ntimes = 0
    do iz = 1, nzones_total
      found = .false.
      do it = 1, ntimes
        if (abs(zone_time(iz) - unique_times_tmp(it)) < 1.0d-12) then
          found = .true.; zone_tidx(iz) = it; exit
        endif
      enddo
      if (.not. found) then
        ! The distinct times must increase in the order of the file: the interpolation in time (Ghost_Q2D_Mapped and
        ! the update of the time series) takes them so. The times are finite here (tested above), so the comparison holds
        if (ntimes > 0) then
          if (zone_time(iz) < unique_times_tmp(ntimes)) then
            write(tnew,'(ES12.5)') zone_time(iz)
            write(told,'(ES12.5)') unique_times_tmp(ntimes)
            write(msg,'(3A,I0,A,I0,5A)') '[ERROR] Q2D file ', trim(filename), ', line ', zone_line(iz), ' (zone ', iz, &
              '): SOLUTIONTIME ', trim(adjustl(tnew)), ' is less than ', trim(adjustl(told)), &
              ', the time of a zone above it: the zones must come in increasing time'
            call mpi_stop_all(trim(msg))
          endif
        endif
        ntimes = ntimes + 1
        unique_times_tmp(ntimes) = zone_time(iz)
        zone_tidx(iz) = ntimes
      endif
    enddo
    allocate(unique_times(ntimes))
    unique_times = unique_times_tmp(1:ntimes)
    deallocate(unique_times_tmp)

    write(*,'(A,I0,A,I0,A)') '  Q2D data: ', ntimes, ' time steps, ', nvar_file, ' variables'

    ! --- Allocate q2d_map for each type 410 BC ---
    do ibc = 1, domain%nbound
      if (domain%bc(ibc)%type == 410 .and. trim(domain%bc(ibc)%q2d_file) == trim(filename)) then
        allocate(domain%bc(ibc)%q2d_map(nprim))
        do iv = 1, nvar_file
          if (varmap(iv) > 0) then
            ic = varmap(iv)
            domain%bc(ibc)%q2d_map(ic)%exists = .true.
            domain%bc(ibc)%q2d_map(ic)%n = ntimes
            allocate(domain%bc(ibc)%q2d_map(ic)%time(ntimes))
            allocate(domain%bc(ibc)%q2d_map(ic)%var(ntimes))
            domain%bc(ibc)%q2d_map(ic)%time = unique_times
            domain%bc(ibc)%q2d_map(ic)%var = 0.0_R8
            domain%bc(ibc)%q2d_map(ic)%periodic = domain%bc(ibc)%q2d_periodic
          endif
        enddo
        ! Extra species not in file: set to zero for all times
        do ic = 1, nsc
          if (.not. domain%bc(ibc)%q2d_map(ic)%exists) then
            domain%bc(ibc)%q2d_map(ic)%exists = .true.
            domain%bc(ibc)%q2d_map(ic)%n = ntimes
            allocate(domain%bc(ibc)%q2d_map(ic)%time(ntimes))
            allocate(domain%bc(ibc)%q2d_map(ic)%var(ntimes))
            domain%bc(ibc)%q2d_map(ic)%time = unique_times
            domain%bc(ibc)%q2d_map(ic)%var = 0.0_R8
            domain%bc(ibc)%q2d_map(ic)%periodic = domain%bc(ibc)%q2d_periodic
          endif
        enddo
      endif
    enddo

    ! --- Pass 2: read zone data and distribute to BCs ---

    ! Allocate buffers once using max dimensions
    I_max = maxval(zone_I(1:nzones_total))
    J_max = maxval(zone_J(1:nzones_total))
    allocate(node_buf(I_max * J_max))
    allocate(cell_buf(nvar_file, I_max - 1, J_max - 1))

    rewind(u)
    ! Skip to first ZONE line
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) call mpi_stop_all('[ERROR] Q2D file '//trim(filename)//': no ZONE found in the data pass')
      if (Is_Zone_Record(line)) exit
    enddo

    do iz = 1, nzones_total
      ! First zone header already read above; subsequent ones read at end of loop

      I_nodes = zone_I(iz)
      J_nodes = zone_J(iz)
      n1 = I_nodes - 1
      n2 = J_nodes - 1
      n_nodes = I_nodes * J_nodes
      b_zone = zone_block(iz)
      t_zone = zone_tidx(iz)

      ! Read and skip node coordinates (3 * n_nodes values); a read that fails stops the run
      read(u, *, iostat=ios) node_buf(1:n_nodes)  ! x
      if (ios /= 0) call Q2D_Read_Error(filename, iz, zone_line(iz), tmp_names(1), ios)
      read(u, *, iostat=ios) node_buf(1:n_nodes)  ! y
      if (ios /= 0) call Q2D_Read_Error(filename, iz, zone_line(iz), tmp_names(2), ios)
      read(u, *, iostat=ios) node_buf(1:n_nodes)  ! z
      if (ios /= 0) call Q2D_Read_Error(filename, iz, zone_line(iz), tmp_names(3), ios)

      ! Read cell-centered variables: nvar_file blocks of n_cells values
      do iv = 1, nvar_file
        read(u, *, iostat=ios) ((cell_buf(iv, i1, i2), i1=1,n1), i2=1,n2)
        if (ios /= 0) call Q2D_Read_Error(filename, iz, zone_line(iz), tmp_names(3+iv), ios)
      enddo

      ! Distribute to matching type 410 BCs
      do ibc = 1, domain%nbound
        if (domain%bc(ibc)%type /= 410) cycle
        if (domain%bc(ibc)%b /= b_zone) cycle
        if (trim(domain%bc(ibc)%q2d_file) /= trim(filename)) cycle

        ! Face-local indices: (i,j,k,f) -> (i1,i2)
        select case (domain%bc(ibc)%f)
          case (1,2); bc_i1 = domain%bc(ibc)%j; bc_i2 = domain%bc(ibc)%k
          case (3,4); bc_i1 = domain%bc(ibc)%i; bc_i2 = domain%bc(ibc)%k
          case default; bc_i1 = domain%bc(ibc)%i; bc_i2 = domain%bc(ibc)%j
        end select

        ! Check bounds
        if (bc_i1 < 1 .or. bc_i1 > n1 .or. bc_i2 < 1 .or. bc_i2 > n2) then
          write(*,'(A,4I6)') ' [WARNING] Q2D BC index out of range: bc=', ibc, &
            bc_i1, bc_i2, b_zone
          cycle
        endif

        ! Map file variables to MOSE primitives
        do iv = 1, nvar_file
          if (varmap(iv) > 0) then
            if (Not_Finite(cell_buf(iv, bc_i1, bc_i2))) then
              write(msg,'(3A,I0,A,I0,3A,I0,A,I0,A)') '[ERROR] Q2D file ', trim(filename), ', zone ', iz, ' (line ', &
                zone_line(iz), '), variable ', trim(tmp_names(3+iv)), ', cell ', bc_i1, ', ', bc_i2, &
                ': a value that is not finite (NaN or infinity)'
              call mpi_stop_all(trim(msg))
            endif
            domain%bc(ibc)%q2d_map(varmap(iv))%var(t_zone) = cell_buf(iv, bc_i1, bc_i2)
          endif
        enddo
      enddo

      ! Read next zone header (if not the last zone): the line after the data of a zone must be the header of the next
      ! zone, or the zones after it would be read shifted (a zone with more values than I and J give, a blank line)
      if (iz < nzones_total) then
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) line = ''
        if (.not. Is_Zone_Record(line)) then
          write(msg,'(3A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] Q2D file ', trim(filename), ', zone ', iz, ' (line ', &
            zone_line(iz), '): the line after its data is not the header of zone ', iz + 1, ' (line ', &
            zone_line(iz + 1), '): more values than I and J give, or a blank line'
          call mpi_stop_all(trim(msg))
        endif
      endif
    enddo

    close(u)
    deallocate(node_buf, cell_buf)
    deallocate(varmap, zone_block, zone_tidx, zone_time, zone_I, zone_J, zone_line, unique_times)

  end subroutine Read_Q2D_Tecplot


  ! ---------------------------------------------------------------------------
  ! Stop the run on a read of the data of a zone that failed (ios /= 0): the
  ! file ends before all the values that I and J give, or a value is not a
  ! number (also the header of the next zone, met where values were due).
  ! ---------------------------------------------------------------------------
  subroutine Q2D_Read_Error ( filename, iz, nline, name, ios )
    use MOSE_Mod_MPI, only: mpi_stop_all
    implicit none
    character(len=*), intent(in) :: filename, name
    integer, intent(in) :: iz, nline, ios
    ! Local
    character(len=640) :: msg

    if (is_iostat_end(ios)) then
      write(msg,'(3A,I0,A,I0,3A)') '[ERROR] Q2D file ', trim(filename), ', zone ', iz, ' (line ', nline, &
        '): the values of ', trim(name), ': the file ends before all the values that I and J give'
    else
      write(msg,'(3A,I0,A,I0,3A)') '[ERROR] Q2D file ', trim(filename), ', zone ', iz, ' (line ', nline, &
        '): the values of ', trim(name), ': a value that is not a number, or fewer values than I and J give'
    endif
    call mpi_stop_all(trim(msg))

  end subroutine Q2D_Read_Error


  ! ---------------------------------------------------------------------------
  ! Build mapping from file variable index to MOSE primitive index.
  ! Matching:
  !   (a) "rho(N)" pattern — MOSE canonical (IO_Solution, Mod_Probes)
  !   (b) "u","v","w","p"  — case-insensitive
  ! ---------------------------------------------------------------------------
  subroutine build_varmap(names, nvar, varmap)
    use MOSE_Global_m
    use strings, only: lowercase
    use IR_Precision, only: str

    implicit none
    integer, intent(in) :: nvar
    character(len=*), intent(in) :: names(nvar)
    integer, intent(out) :: varmap(nvar)
    ! Local
    integer :: iv, s
    character(len=64) :: lname

    varmap = 0

    do iv = 1, nvar
      lname = lowercase(trim(names(iv)))

      ! (a) "rho(N)" pattern (MOSE canonical: IO_Solution, Mod_Probes)
      do s = 1, nsc
        if (trim(lname) == 'rho('//trim(str(.true.,s))//')') then
          varmap(iv) = s
          exit
        endif
      enddo
      if (varmap(iv) /= 0) cycle

      ! (b) velocity and pressure (case-insensitive)
      select case (trim(lname))
        case ('u');  varmap(iv) = nu
        case ('v');  varmap(iv) = nv
        case ('w');  varmap(iv) = nw
        case ('p');  varmap(iv) = np
      end select
    enddo

    ! Log the mapping
    write(*,'(A)') '  Q2D variable mapping:'
    do iv = 1, nvar
      if (varmap(iv) > 0) then
        write(*,'(A,A,A,I3,A)') '    "', trim(names(iv)), '" -> P(', varmap(iv), ')'
      else
        write(*,'(A,A,A)')    '    "', trim(names(iv)), '" -> (unmapped)'
      endif
    enddo

  end subroutine build_varmap


  ! ---------------------------------------------------------------------------
  ! Split the text after "VARIABLES =" into names. A name is in double quotes,
  ! in single quotes or bare; names are separated by blanks, tabs or commas, and
  ! a quoted name may hold any of them. ok = .false. on an unclosed quote or
  ! more names than size(names).
  ! ---------------------------------------------------------------------------
  subroutine Parse_Variable_Names(text, names, n, ok)
    implicit none
    character(len=*), intent(in)  :: text
    character(len=*), intent(out) :: names(:)
    integer,          intent(out) :: n
    logical,          intent(out) :: ok
    ! Local
    character(len=*), parameter :: seps = ' ,'//achar(9)
    integer :: ip, jp, last

    names = ''; n = 0; ok = .true.
    last = len_trim(text)
    ip = 1
    do while (ip <= last)
      if (index(seps, text(ip:ip)) > 0) then
        ip = ip + 1
        cycle
      endif
      if (n == size(names)) then
        ok = .false.
        return
      endif
      n = n + 1
      if (text(ip:ip) == '"' .or. text(ip:ip) == "'") then
        jp = index(text(ip+1:last), text(ip:ip))
        if (jp == 0) then
          ok = .false.
          return
        endif
        names(n) = text(ip+1:ip+jp-1)
        ip = ip + jp + 1
      else
        jp = scan(text(ip:last), seps//'"'//"'")
        if (jp == 0) jp = last - ip + 2
        names(n) = text(ip:ip+jp-2)
        ip = ip + jp - 1
      endif
    enddo

  end subroutine Parse_Variable_Names


  ! ---------------------------------------------------------------------------
  ! True when the first word of a line is the ZONE keyword (any case), so that
  ! a variable name such as zone_velocity is not taken for a zone record.
  ! ---------------------------------------------------------------------------
  logical function Is_Zone_Record(line)
    use strings, only: lowercase
    implicit none
    character(len=*), intent(in) :: line
    ! Local
    character(len=:), allocatable :: word
    integer :: k

    word = adjustl(line)
    k = scan(word, ' ,'//achar(9))
    if (k > 0) word = word(:k-1)
    Is_Zone_Record = trim(lowercase(word)) == 'zone'

  end function Is_Zone_Record


  ! ---------------------------------------------------------------------------
  ! Parse a Tecplot ZONE header line to extract block number,
  ! I (node count), J (node count), and SOLUTIONTIME.
  ! Zone name format: T="B{block}_T{tidx}"
  ! bad names the first of these values whose text (badval) does not read,
  ! blank when they all read; a value not in the line keeps its default.
  ! ---------------------------------------------------------------------------
  subroutine parse_zone_header(line, b_num, I_nodes, J_nodes, sol_time, bad, badval)
    use strings, only: parse, lowercase
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(out) :: b_num, I_nodes, J_nodes
    real(R8), intent(out) :: sol_time
    character(len=*), intent(out) :: bad, badval
    ! Local
    integer :: i, pos, pos2, ios
    character(len=100) :: args(20), subargs(4)
    character(len=256) :: zone_title

    b_num = 1; I_nodes = 2; J_nodes = 2; sol_time = 0.0_R8
    bad = ''; badval = ''

    ! Split line by commas
    call parse(line, ',', args)

    do i = 1, size(args)
      ! Zone title: T="B{b}_T{t}", T = 'B{b}_T{t}' or T = B{b}_T{t}
      pos = index(args(i), '=')
      if (pos > 1) then
        ! key = last word before '=' (the first argument also holds ZONE)
        zone_title = trim(args(i)(:pos-1))
        zone_title = zone_title(index(trim(zone_title), ' ', back=.true.)+1:)
        if (lowercase(trim(zone_title)) /= 't') pos = 0
      endif
      if (pos > 1) then
        zone_title = adjustl(args(i)(pos+1:))
        pos2 = len_trim(zone_title)
        if (pos2 >= 2 .and. scan(zone_title(1:1), '"'//"'") == 1) then
          if (zone_title(pos2:pos2) == zone_title(1:1)) zone_title = zone_title(2:pos2-1)
        endif
        if (len_trim(zone_title) > 0) then
          ! Parse "B{num}_T{num}"
          pos = index(zone_title, 'B')
          if (pos > 0) then
            pos2 = index(zone_title(pos+1:), '_')
            if (pos2 > 0) then
              read(zone_title(pos+1:pos+pos2-1), *, iostat=ios) b_num
              if (ios /= 0) call Mark_Bad('the block number of the title', zone_title(pos+1:pos+pos2-1))
            endif
          endif
        endif
      endif

      ! I= (avoid matching STRANDID=)
      if (index(args(i), 'I=') > 0 .and. index(args(i), 'STRANDID') == 0) then
        call parse(args(i), '=', subargs)
        read(subargs(2), *, iostat=ios) I_nodes
        if (ios /= 0) call Mark_Bad('the value of I', subargs(2))
      endif

      ! J=
      if (index(args(i), 'J=') > 0) then
        call parse(args(i), '=', subargs)
        read(subargs(2), *, iostat=ios) J_nodes
        if (ios /= 0) call Mark_Bad('the value of J', subargs(2))
      endif

      ! SOLUTIONTIME=
      if (index(lowercase(args(i)), 'solutiontime') > 0) then
        call parse(args(i), '=', subargs)
        read(subargs(2), *, iostat=ios) sol_time
        if (ios /= 0) call Mark_Bad('the value of SOLUTIONTIME', subargs(2))
      endif
    enddo

  contains

    subroutine Mark_Bad ( what, text )
      character(len=*), intent(in) :: what, text
      if (len_trim(bad) == 0) then
        bad = what
        badval = adjustl(text)
      endif
    end subroutine Mark_Bad

  end subroutine parse_zone_header

end module MOSE_IO_BC_Q2D