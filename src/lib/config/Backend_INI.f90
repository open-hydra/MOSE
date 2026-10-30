module MOSE_Backend_INI
  implicit none
  private
  public :: Load_Ini, Scan_Ini, Count_Stationary_Faces, Lists_Changed

  ! The numbered lists of input.ini that the set-up read (Scan_Ini), against which a reading during the run checks them
  ! (Lists_Changed): the number of grid levels, whether input.ini gave levels, and the number of probes
  integer :: levels_start = 1
  logical :: levels_given = .false.
  integer :: probes_start = 0

contains

  subroutine Load_Ini(fini)
    use MOSE_Input_Registry, only: reg
    use Finer, only: file_ini, err_option_vals
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: i, error, error_text
    character(len=:), allocatable :: text

    do i = 1, reg%size
      error = 1

      ! Each registry entry has exactly one associated typed pointer.
      if (associated(reg%params(i)%value%i)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%i, error=error)
      else if (associated(reg%params(i)%value%r)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%r, error=error)
      else if (associated(reg%params(i)%value%l)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%l, error=error)
      else if (associated(reg%params(i)%value%s)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%s, error=error)
      else if (associated(reg%params(i)%value%iarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%iarr, error=error)
      else if (associated(reg%params(i)%value%rarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%rarr, error=error)
      end if

      if (error == 0) reg%params(i)%is_set = .true.

      ! A key of the file whose value the reader cannot take as the type of its parameter (a value that does not
      ! convert, a list with more values than the parameter holds, no value: FiNeR gives err_option_vals, and leaves
      ! the parameter as it was) is an error of the file, which Validate_Registry gives with the value as written
      reg%params(i)%unread = ( error == err_option_vals )
      if (allocated(reg%params(i)%text)) deallocate(reg%params(i)%text)
      call fini%get_string(section_name=reg%params(i)%section, option_name=reg%params(i)%name, val=text, error=error_text)
      if (error_text == 0) reg%params(i)%text = text
      ! An integer, or a list of integers, written with a period or an exponent is not one (the input form of an integer is
      ! a signed digit string): gfortran refuses it, ifx takes it and cuts the fraction (2.7 reads 2), so the text decides,
      ! in every build
      if (error == 0 .and. error_text == 0 .and. &
          (associated(reg%params(i)%value%i) .or. associated(reg%params(i)%value%iarr))) then
        if (scan(text, '.eEdDqQ') > 0) reg%params(i)%unread = .true.
      end if
    end do

  end subroutine


  subroutine Scan_Ini(fini, nprobes, probes_name, nmgl)
    use Finer, only: file_ini
    use MOSE_Input_Registry, only: reg
    use MOSE_Lib_RotatingFrame, only: obj_rot
    use MOSE_Global_m
    use IR_precision
    implicit none
    type(file_ini), intent(in) :: fini
    integer, intent(out) :: nprobes, nmgl
    character(len=16), allocatable, intent(out) :: probes_name(:)
    ! Local
    integer :: error, i
    character(len=:), allocatable :: section, msg
  
    ! Read the MGL. A value that is not an integer, or no value: the run would go on with one grid level, without a word.
    ! levels is no key of the registry: Validate_Registry gives this error first (Check_Input), as those of its keys
    nmgl = 1
    call Read_Levels(fini, nmgl, levels_given, section, msg)
    if (len(msg) > 0) reg%load_error = msg
    levels_start = nmgl
    
    ! Count the probes
    nprobes = Count_Probes(fini)
    probes_start = nprobes

    allocate(probes_name(nprobes))
    do i = 1, nprobes
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,i)), val=probes_name(i), error=error)
    enddo
    
    ! Count the stationary faces for the rotating frame: stationary-face-1, -2, ... up to the first missing
    obj_rot%n_stationary = Count_Stationary_Faces(fini)
    if (obj_rot%n_stationary > 0) allocate(obj_rot%stationary_face_str(obj_rot%n_stationary))
    
  end subroutine Scan_Ini


  !> levels of the multigrid section of fini ([MOSE-Multigrid], or [HYDRA-Multigrid] when HYDRA runs MOSE), which section
  !> names: n takes its value, and keeps the one it has when the key is not there (given false). msg is the error of a
  !> value that is not an integer, or of no value, or an empty string: an integer written with a period or an exponent is
  !> not one (gfortran refuses it, ifx takes it and cuts the fraction). The same reading at set-up (Scan_Ini) and when
  !> input.ini is read again during the run (Lists_Changed).
  subroutine Read_Levels(fini, n, given, section, msg)
    use Finer, only: file_ini, err_option_vals
    use MOSE_Config_Types_m, only: obj_sim_param
    use MOSE_Global_m
    implicit none
    type(file_ini), intent(in) :: fini
    integer, intent(inout) :: n
    logical, intent(out) :: given
    character(len=:), allocatable, intent(out) :: section, msg
    ! Local
    integer :: error, error_text
    character(len=:), allocatable :: text

    if (.not. obj_sim_param%HYDRA_MG) then
      section = trim(codename)//'-Multigrid'
    else
      section = 'HYDRA-Multigrid'
    endif
    msg = ''
    call fini%get(section_name=section, option_name='levels', val=n, error=error)
    text = ''
    call fini%get_string(section_name=section, option_name='levels', val=text, error=error_text)
    given = ( error_text == 0 )
    if (error == err_option_vals .or. (error == 0 .and. error_text == 0 .and. scan(text, '.eEdDqQ') > 0)) then
      if (len_trim(text) > 0) then
        msg = '[ERROR] input.ini, ['//section//'] levels = '//trim(text)//': the value is not an integer'
      else
        msg = '[ERROR] input.ini, ['//section//'] levels =: the value is empty'
      endif
    endif

  end subroutine Read_Levels


  !> The number of probes in fini: probe1, probe2, ... of [MOSE-Probes] up to the first missing (at set-up, and when
  !> input.ini is read again during the run).
  function Count_Probes(fini) result(n)
    use Finer, only: file_ini
    use MOSE_Global_m
    use IR_precision
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: n
    ! Local
    character(len=llen) :: wholestring
    integer :: error

    n = 0
    do
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,n+1)), val=wholestring, error=error)
      if (error/=0) exit
      n = n+1
    enddo

  end function Count_Probes


  !> The error of a numbered list of input.ini read again during the run that is not the one of the set-up, or an empty
  !> string. levels, read as at set-up (Read_Levels): a value that is not an integer, or another number than the set-up read
  !> (the grid levels are built at set-up). The probes, counted as at set-up (Count_Probes): more than at the set-up (the
  !> probes are placed and their files opened at set-up). A probe key that names another probe is a key fixed for the run
  !> (Fixed_Key_Changed), and a key that the file leaves out keeps its probe, as every key keeps its value. The stationary
  !> faces of the rotating frame are checked by Reload_RotatingFrame.
  function Lists_Changed(fini) result(out)
    use Finer, only: file_ini
    use MOSE_Global_m
    use IR_precision
    implicit none
    type(file_ini), intent(in) :: fini
    character(len=:), allocatable :: out
    ! Local
    integer :: n
    logical :: given
    character(len=:), allocatable :: section

    n = levels_start
    call Read_Levels(fini, n, given, section, out)
    if (len(out) > 0) return
    if (n /= levels_start) then
      if (levels_given) then
        out = '[ERROR] input.ini, ['//section//'] levels = '//trim(str(.true.,n))//': the run has levels = '// &
              trim(str(.true.,levels_start))//' from its set-up, and levels cannot change during the run'
      else
        out = '[ERROR] input.ini, ['//section//'] levels = '//trim(str(.true.,n))//': the run started without levels in '// &
              'input.ini (its default: 1), and levels cannot change during the run'
      endif
      return
    endif

    n = Count_Probes(fini)
    if (n > probes_start) out = '[ERROR] input.ini, ['//trim(codename)//'-Probes] probe keys (probe1, probe2, ... up to the '// &
      'first missing): '//trim(str(.true.,n))//' read again, '//trim(str(.true.,probes_start))//' at the set-up: the probes '// &
      'cannot change during the run'

  end function Lists_Changed


  !> The number of stationary-face-<k> keys of [MOSE-Rotating-Frame] in fini: stationary-face-1, -2, ... up to the
  !> first missing (at set-up, and when input.ini is read again during the run).
  function Count_Stationary_Faces(fini) result(n)
    use Finer, only: file_ini
    use MOSE_Global_m
    use IR_precision
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: n
    ! Local
    character(len=llen) :: wholestring
    integer :: error

    n = 0
    do
      call fini%get(section_name=trim(codename)//'-Rotating-Frame', &
                    option_name='stationary-face-'//trim(str(.true.,n+1)), val=wholestring, error=error)
      if (error /= 0) exit
      n = n + 1
    enddo

  end function Count_Stationary_Faces

end module MOSE_Backend_INI