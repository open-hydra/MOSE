! FLINT wrapper to load chemistry data
module MOSE_Load_Chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_Chemistry

contains

  subroutine Load_Chemistry()
    use MOSE_Config_Types_m,  only: obj_chemistry
    use FLINT_Load_chemistry, only: read_chemistry
    implicit none
    ! Local
    integer :: error
    
    error = read_chemistry(mech_name = obj_chemistry%mechanism_name)
    if (error/=0) then
      select case (error)
        case (0)
          obj_chemistry%description = 'See next section for details'
        case (1)
          obj_chemistry%description = 'Frozen mixture'
        case (2)
          obj_chemistry%error_message = '[ERROR] Chemistry info file found but could not be read'
        case (3)
          obj_chemistry%error_message = '[ERROR] Chemistry rate table file not found'
        case (4)
          obj_chemistry%error_message = '[ERROR] Chemistry rate table file found but could not be read'
        case default
          obj_chemistry%error_message = '[ERROR] Unknown error loading chemistry'
      end select
    endif

  end subroutine Load_Chemistry

end module MOSE_Load_Chemistry