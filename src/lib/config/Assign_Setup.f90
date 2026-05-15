module MOSE_Assign_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ir_precision
  
  implicit none
  private
  public :: Assign_Setup

contains

  subroutine Assign_Setup()
    use MOSE_Config_Types_m
    use MOSE_Global_m,            only: model
    use MOSE_IO_Solution,         only: Setup_Input_Solution
    use MOSE_Mod_Space,           only: Setup_Space_Scheme
    use MOSE_Mod_Riemann,         only: Assign_Riemann_Solver
    use MOSE_Mod_Newstate,        only: Assign_Integration_Variables
    use MOSE_Lib_Strang,          only: Setup_Strang_Splitting
    use MOSE_Lib_Chemistry,       only: Setup_Chemistry
    use MOSE_Mod_Soot,            only: Setup_Soot
    use MOSE_Mod_RANS,            only: Setup_RANS_Model
    use MOSE_Lib_RotatingFrame,   only: Setup_RotatingFrame
    implicit none

    ! Setting simulation type
    if (obj_sim_param%simulation_type=='euler') then
      model = 0
    elseif (obj_sim_param%simulation_type=='navier-stokes') then
      if (obj_rans%model == 'none' .or. obj_rans%model == '') then
        model = 1
        obj_rans%model = 'none'
      else
        model = 2
      end if
    endif

    ! Setting input solution
    call Setup_Input_Solution()

    ! Space
    call Setup_Space_Scheme()
    call Assign_Riemann_Solver()
    if ( obj_riemann%SD .or. obj_space_scheme%SD) then
      obj_shock_detector%SD = .true.
    else
      obj_shock_detector%SD = .false.
    end if

    ! Time
    if (obj_time_scheme%solver_type /= 'euler') then
      read(obj_time_scheme%solver_type(3:3), *) obj_time_scheme%n_rk
    else
      obj_time_scheme%n_rk = 1
    end if
    if (obj_irs%beta>0d0) obj_irs%enabled = .true. 
    call Assign_Integration_Variables()

    ! Assign Chemistry
    call Setup_Chemistry()
    call Setup_Strang_Splitting()
    
    ! Assign soot model
    call Setup_Soot()

    ! Assign RANS model
    call Setup_RANS_Model()

    ! Assign Rotating frame
    call Setup_RotatingFrame()

    !! Descriptions, warnings and errors

    ! Simulation type
    if (obj_sim_param%simulation_type == 'euler') then
      obj_sim_param%description = 'Euler'
    else if (obj_sim_param%simulation_type == 'navier-stokes') then
      obj_sim_param%description = 'Navier-Stokes'
    end if
    ! Time scheme
    if (obj_time_scheme%solver_type == 'euler') then
      obj_time_scheme%description = 'Explicit Euler'
    else if (obj_time_scheme%solver_type == 'RK2') then
      obj_time_scheme%description = 'Second-order Runge-Kutta'
    else if (obj_time_scheme%solver_type == 'RK3') then
      obj_time_scheme%description = 'Third-order Runge-Kutta'
    end if
    if (obj_time_scheme%time_accurate) then
      obj_time_scheme%description = trim(obj_time_scheme%description)//' with time-accurate switch enabled'
    end if
    if (obj_irs%enabled) then
      obj_irs%description = 'Beta set to '//trim(str(.true.,real(obj_irs%beta)))
    end if
    ! Space scheme
    ! ... written in Mod_Space ...
    ! Riemann solver
    if (trim(obj_riemann%description) == 'AUSM+-up' .or. trim(obj_riemann%description) == 'AUSM+-up2') then
      if (obj_riemann%Minf == 0.0d0) then
        obj_riemann%error_message = '[ERROR] AUSM+-up solver selected. Minf must be defined in input.'
      end if
    endif
    ! Transport
    if (model>0 .and. obj_transport%description=='Unavailable') &
    write(*,'(A)') '[ERROR] Transport properties are unavailable for the selected phase: cannot run Navier-Stokes simulation'

  end subroutine Assign_Setup

end module MOSE_Assign_Setup