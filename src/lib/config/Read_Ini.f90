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

    ! Load input.ini: a file that cannot be read stops the run
    call Load_Input_File(fini, 'Input file')

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


  !> Prints a warning for every line "key = value" of input.ini whose key the INI reader (FiNeR) does not give in its
  !> section: that key is not read, and its parameter keeps its default. Comment lines (they start with !, ; or #) are
  !> skipped, as the reader skips them since FiNeR 2.1.0, also when they hold "=" or a bracket.
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
        '] is not read (the INI reader does not give it): its parameter keeps its default'
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

    ! Load input.ini: a file that cannot be read stops the run
    call Load_Input_File(fini, 'Input file, read again')

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


  !> Loads input.ini into fini, and stops the run if it cannot be read: a file that cannot be opened (missing, or without
  !> read permission), a file with no line that is not blank, or an error of the INI reader (FiNeR load gives an error
  !> code other than 0: a file with no section, a line that it cannot parse; FiNeR names such a line on stderr, without
  !> its number). Every rank reads input.ini: the ranks take the decision together (the largest of a flag), then every
  !> rank stops with mpi_stop_all (exit code 1). Without this stop, the parameters kept their defaults at set-up and the
  !> values of the previous reading during the run, and an empty file ended the run inside the reader.
  subroutine Load_Input_File(fini, label)
    use Finer,           only: file_ini
    use MOSE_Mod_MPI,    only: mpi_allreduce_max_r8, mpi_stop_all
    use iso_fortran_env, only: R8 => real64
    implicit none
    type(file_ini),   intent(inout) :: fini
    character(len=*), intent(in)    :: label   !< Line of the set-up summary: Input file, or Input file, read again
    ! Local
    character(len=512)  :: reason, msg
    character(len=1024) :: line
    integer  :: unit, ios, error, n
    real(R8) :: flag, flag_all

    reason = ''
    open(newunit=unit, file='input.ini', status='old', action='read', iostat=ios, iomsg=msg)
    if (ios /= 0) then
      reason = 'cannot be opened: '//trim(msg)
    else
      n = 0
      do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) > 0) n = n + 1
      end do
      close(unit)
      if (n == 0) reason = 'is empty'
    end if
    if (len_trim(reason) == 0) then
      call fini%load(filename='input.ini', error=error)
      if (error /= 0) write(reason, '(A,I0,A)') 'cannot be read by the INI reader (FiNeR error ', error, &
        '; a line that FiNeR could not parse is named above, if any)'
    end if

    flag = 0.0_R8
    if (len_trim(reason) > 0) flag = 1.0_R8
    call mpi_allreduce_max_r8(flag, flag_all)
    if (flag_all > 0.0_R8) then
      if (len_trim(reason) > 0) then
        write(*,'(A,T35,A)') '   '//trim(label), 'FAIL'
        write(*,'(4X,A)') '[ERROR] input.ini '//trim(reason)
      end if
      call mpi_stop_all('input.ini cannot be read (see above)')
    end if

  end subroutine Load_Input_File


end module MOSE_Read_Ini