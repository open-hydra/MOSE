module MOSE_Backend_INI
  implicit none
  private
  public :: Load_Ini, Scan_Ini

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
    use Finer, only: file_ini, err_option_vals
    use MOSE_Input_Registry, only: reg
    use MOSE_Config_Types_m, only: obj_sim_param
    use MOSE_Lib_RotatingFrame, only: obj_rot
    use MOSE_Global_m
    use IR_precision
    implicit none
    type(file_ini), intent(in) :: fini
    integer, intent(out) :: nprobes, nmgl
    character(len=16), allocatable, intent(out) :: probes_name(:)
    ! Local
    character(len=llen) :: wholestring
    integer :: error, error_text, i
    character(len=:), allocatable :: section, text
  
    ! Read the MGL
    nmgl = 1
    if (.not. obj_sim_param%HYDRA_MG) then
      section = trim(codename)//'-Multigrid'
    else
      section = 'HYDRA-Multigrid'
    endif
    call fini%get(section_name=section, option_name='levels', val=nmgl, error=error)
    ! A value that is not an integer, or no value: the run would go on with one grid level, without a word; an integer
    ! written with a period or an exponent is not one either (gfortran refuses it, ifx takes it and cuts the fraction).
    ! levels is no key of the registry: Validate_Registry gives this error first (Check_Input), as those of its keys
    text = ''
    call fini%get_string(section_name=section, option_name='levels', val=text, error=error_text)
    if (error == err_option_vals .or. (error == 0 .and. error_text == 0 .and. scan(text, '.eEdDqQ') > 0)) then
      if (len_trim(text) > 0) then
        reg%load_error = '[ERROR] input.ini, ['//section//'] levels = '//trim(text)//': the value is not an integer'
      else
        reg%load_error = '[ERROR] input.ini, ['//section//'] levels =: the value is empty'
      endif
    endif
    
    ! Count the probes
    nprobes = 0
    do 
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,nprobes+1)), val=wholestring, error=error)
      if (error/=0) exit
      nprobes = nprobes+1
    enddo

    allocate(probes_name(nprobes))
    do i = 1, nprobes
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,i)), val=probes_name(i), error=error)
    enddo
    
    ! Count the stationary faces for the rotating frame: stationary-face-1, -2, ... up to the first missing
    obj_rot%n_stationary = 0
    do
      call fini%get(section_name=trim(codename)//'-Rotating-Frame', &
                    option_name='stationary-face-'//trim(str(.true.,obj_rot%n_stationary+1)), val=wholestring, error=error)
      if (error /= 0) exit
      obj_rot%n_stationary = obj_rot%n_stationary + 1
    enddo
    if (obj_rot%n_stationary > 0) allocate(obj_rot%stationary_face_str(obj_rot%n_stationary))
    
  end subroutine Scan_Ini

end module MOSE_Backend_INI