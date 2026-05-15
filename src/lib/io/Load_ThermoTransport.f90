! FLINT wrapper to load thermodynamic data
module MOSE_Load_ThermoTransport
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_ThermoTransport

contains

  subroutine Load_ThermoTransport()
    use MOSE_Config_Types_m,        only: obj_thermo, obj_transport
    use MOSE_Global_m,              only: nsc, MOSE_phase_prefix
    use FLINT_Load_ThermoTransport, only: read_idealgas_thermo, read_idealgas_transport
    use FLINT_Lib_Thermodynamic,    only: FLINT_phase_prefix, ns
    use MOSE_Parameters_m
    implicit none
    ! Local
    integer :: error

    !! ------------------------------------------------------
    !! Thermo-Transport Tables ------------------------------
    !! ------------------------------------------------------
    obj_thermo%warning_message = 'none'
    obj_thermo%error_message   = 'none'
    obj_thermo%description     = 'none'

    FLINT_phase_prefix = MOSE_phase_prefix
    error = read_idealgas_thermo( 'INPUT' )

    select case (error)
      case (0)
        obj_thermo%description = 'Thermally perfect gas'
      case (5)
        obj_thermo%description = 'Thermally perfect gas'
      case (1)
        obj_thermo%error_message = '[ERROR] Phase file (phase.txt) not found'
      case (2)
        obj_thermo%error_message = '[ERROR] Phase file (phase.txt) found but could not be read'
      case (3)
        obj_thermo%error_message = '[ERROR] Thermo table file not found'
      case (4)
        obj_thermo%error_message = '[ERROR] Thermo table file found but could not be read'
      case (6)
        obj_thermo%error_message = '[ERROR] Composition file found but could not be read'
      case default
        obj_thermo%error_message = '[ERROR] Unknown error loading thermodynamic data'
    end select

    obj_transport%warning_message = 'none'
    obj_transport%error_message   = 'none'
    obj_transport%description     = 'none'

    error = read_idealgas_transport( 'INPUT' )
    select case (error)
      case (0)
        obj_transport%description = 'Mixture-averaged transport properties'
      case (1)
        obj_transport%description = 'Unavailable'
      case (2)
        obj_transport%error_message = '[ERROR] Transport table file found but could not be read'
      case (3)
        obj_transport%error_message = '[ERROR] Transport data temperature range does not match thermo data range'
      case default
        obj_transport%error_message = '[ERROR] Unknown error loading transport data'
    end select

    nsc = ns ! Set the number of species from FLINT
    !! ------------------------------------------------------
    !! ------------------------------------------------------

  end subroutine Load_ThermoTransport

end module MOSE_Load_ThermoTransport