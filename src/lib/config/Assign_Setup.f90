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

    !! Setting simulation type
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

    !! Setting input solution
    call Setup_Input_Solution()

    !! Space
    call Setup_Space_Scheme()

    !! Riemann
    call Assign_Riemann_Solver()

    !! Shock detector
    obj_shock_detector%id = 0
    if (obj_shock_detector%description=='tramel') then
      obj_shock_detector%id = 1; obj_shock_detector%description='Tramel'
    elseif (obj_shock_detector%description=='chen') then
      obj_shock_detector%id = 2; obj_shock_detector%description='Chen'
    elseif (obj_shock_detector%description=='' .and. index(obj_riemann%description,'Tramel')>0) then
      obj_shock_detector%id = 1; obj_shock_detector%description='Tramel'
    elseif (obj_shock_detector%description=='' .and. index(obj_riemann%description,'Chen')>0) then
      obj_shock_detector%id = 2; obj_shock_detector%description='Chen'
    endif

    !! Time
    ! Time scheme
    if (obj_time_scheme%solver_type /= 'euler') then
      read(obj_time_scheme%solver_type(3:3), *) obj_time_scheme%n_rk
    else
      obj_time_scheme%n_rk = 1
    end if
    ! Implicit residual smoothing
    if (obj_irs%beta>0d0) obj_irs%enabled = .true.
    ! Preconditioning
    if ( trim(obj_time_scheme%integration_variables) == 'prec' ) then
      obj_prec%enabled = .true.
    else
      obj_prec%enabled = .false.
    end if
    ! Integration variables
    call Assign_Integration_Variables()

    !! Assign Chemistry
    call Setup_Chemistry()
    call Setup_Strang_Splitting()
    
    !! Assign soot model
    call Setup_Soot()

    !! Assign RANS model
    call Setup_RANS_Model()

    !! Assign Rotating frame
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
    if (trim(obj_time_scheme%integration_variables) == 'Preconditioned' .and. obj_time_scheme%time_accurate) then
      obj_time_scheme%error_message = '[ERROR] integration-variables=prec is currently supported only for steady/pseudo-time runs.'
    end if
    if (obj_irs%enabled) then
      obj_irs%description = 'Beta set to '//trim(str(.true.,real(obj_irs%beta)))
    end if
    ! Space scheme
    ! ... written in Mod_Space ...
    ! Riemann solver
    if (index(obj_riemann%description, 'AUSM+M')      > 0 .or. &
        index(obj_riemann%description, 'HLLC+ (Chen') > 0 .or. &
        index(obj_riemann%description, 'Low-Mach Roe') > 0) then
      if (obj_riemann%Mco == 0.0d0) then
        obj_riemann%error_message = '[ERROR] Low-Mach Riemann solver selected. Mco (riemann-options-Mco) must be defined in input.'
      end if
    endif
    if (trim(obj_riemann%description) == 'HLLC-PC') then
      if (obj_time_scheme%integration_variables /= 'prec') then
        obj_riemann%error_message = '[ERROR] Preconditioned HLL solver selected. integration-variables must be set to "prec".'
      end if
    endif
    ! Transport
    if (model>0 .and. obj_transport%description=='Unavailable') &
    write(*,'(A)') '[ERROR] Transport properties are unavailable for the selected phase: cannot run Navier-Stokes simulation'

  end subroutine Assign_Setup

end module MOSE_Assign_Setup