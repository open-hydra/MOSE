module MOSE_Read_Ini
  
  implicit none

contains

  subroutine Read_Inifile()
    use MOSE_Global_m,       only: nsc
    use Finer,               only: file_ini
    use MOSE_Read_Sim_Param, only: Register_Sim_Param
    use MOSE_Read_IO,        only: Register_IO_Fields, Register_Probes
    use MOSE_Read_Numerics,  only: Register_Numerics
    use MOSE_Read_Physics,   only: Register_Physics
    use MOSE_Backend_INI,    only: Load_Ini, Scan_Ini
    use MOSE_Input_Registry
    use MOSE_Config_Types_m
    implicit none
    ! Local
    type(file_ini) :: fini
    integer :: nprobes, nmgl
    character(len=16), allocatable :: probes_name(:)

    allocate(obj_chemistry%RT(1:nsc+1), obj_chemistry%AT(1:nsc+1))

    ! Load input.ini
    call fini%load(filename='input.ini')

    ! Scan input.ini for unknown number of probes and multigrid levels
    call Scan_Ini(fini, nprobes, probes_name, nmgl)

    ! Build registry entries
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(nprobes, probes_name)
    call Register_Numerics(nmgl)
    call Register_Physics()

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)

    ! The keys of input.ini that the reader does not give
    call Warn_Unread_Keys(fini)

  end subroutine Read_Inifile


  !> The INI reader (FiNeR) stops reading a section at the first line it cannot parse, such as a comment line (one that starts
  !> with !, ; or #) that holds "=": the keys after it are not read, and their parameters keep their defaults. Prints a
  !> warning for every line "key = value" of input.ini whose key the reader does not give in its section.
  subroutine Warn_Unread_Keys(fini)
    use Finer,        only: file_ini
    use MOSE_Mod_MPI, only: mpi_is_root
    implicit none
    type(file_ini), intent(in) :: fini
    ! Local
    character(len=1024) :: line, section, key, val
    integer :: unit, ios, n, eq, error

    if (.not. mpi_is_root) return
    open(newunit=unit, file='input.ini', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    section = ''
    n = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      n = n + 1
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (scan(line(1:1), '!;#') == 1) cycle
      if (line(1:1) == '[') then
        if (index(line, ']') > 2) section = adjustl(line(2:index(line, ']')-1))
        cycle
      end if
      eq = index(line, '=')
      if (eq <= 1) cycle
      key = adjustl(line(1:eq-1))
      call fini%get(section_name=trim(section), option_name=trim(key), val=val, error=error)
      if (error /= 0) write(*,'(A,I0,A)') '[WARNING] input.ini, line ', n, ': '//trim(key)//' of ['//trim(section)// &
        '] is not read (the reader stops a section at a line it cannot parse, such as a comment line that holds "=")'
    end do
    close(unit)

  end subroutine Warn_Unread_Keys


  subroutine Read_Inifile_Runtime()
    use Finer,               only: file_ini
    use MOSE_Backend_INI,    only: Load_Ini
    use MOSE_Input_Registry
    use MOSE_Mod_MPI,        only: mpi_is_root, mpi_abort_all
    implicit none
    ! Local
    type(file_ini) :: fini
    character(len=1024) :: out

    ! Load input.ini
    call fini%load(filename='input.ini')

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)

    ! Validate registry: the values read again keep to the rules of the registry, as at setup
    ! (Check_Input in Wrap_Setup); a value out of its rule stops the run with the message of the rule
    out = Validate_Registry()
    if (index(out,'ERROR') > 0) then
      if (mpi_is_root) then
        write(*,'(A,T35,A)') '   Input file, read again', 'FAIL'
        write(*,'(4X,A)') trim(out)
      end if
      call mpi_abort_all('input file error, input.ini read again during the run (see above)')
    end if

  end subroutine Read_Inifile_Runtime


end module MOSE_Read_Ini