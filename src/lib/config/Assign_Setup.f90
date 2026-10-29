module MOSE_Assign_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ir_precision
  
  implicit none
  private
  public :: Assign_Setup, IRS_On, Riemann_Mco_Error, Schmidt_Error

contains

  subroutine Assign_Setup()
    use MOSE_Config_Types_m
    use MOSE_Global_m,            only: model
    use MOSE_Parameters_m,        only: codename
    use MOSE_Input_Registry,      only: Key_Is_Set
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
    ! Local
    character(len=:), allocatable :: check

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
    ! shock-detector takes the values of its rule, Tramel and Chen; left out, the sensor that the name of the Riemann solver
    ! carries, if any
    obj_shock_detector%id = 0
    if (obj_shock_detector%description=='Tramel') then
      obj_shock_detector%id = 1
    elseif (obj_shock_detector%description=='Chen') then
      obj_shock_detector%id = 2
    elseif (obj_shock_detector%description=='' .and. index(obj_riemann%name,'Tramel')>0) then
      obj_shock_detector%id = 1; obj_shock_detector%description='Tramel'
    elseif (obj_shock_detector%description=='' .and. index(obj_riemann%name,'Chen')>0) then
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
    obj_irs%enabled  = IRS_On()
    obj_irs%at_setup = obj_irs%enabled
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
    if (trim(obj_time_scheme%integration_variables) == 'prec' .and. obj_time_scheme%time_accurate) then
      obj_time_scheme%error_message = '[ERROR] integration-variables=prec is currently supported only for steady/pseudo-time runs.'
    end if
    if (obj_irs%enabled) then
      obj_irs%description = 'Beta set to '//trim(str(.true.,real(obj_irs%beta)))// &
                            ', variables: '//trim(obj_irs%variables)
    end if
    ! IRS accelerates steady runs: the smoothed residual changes the transient of a time-accurate run
    if (obj_irs%enabled .and. obj_time_scheme%time_accurate) then
      obj_irs%error_message = '[ERROR] Implicit residual smoothing (irs, irs-beta > 0) is supported only for '// &
                              'steady/pseudo-time runs.'
    end if
    ! Space scheme
    ! ... written in Mod_Space ...
    ! Riemann solver
    check = Riemann_Mco_Error()
    if (len(check) > 0) obj_riemann%error_message = check
    if (trim(obj_riemann%name) == 'HLLC-PC') then
      if (obj_time_scheme%integration_variables /= 'prec') then
        obj_riemann%error_message = '[ERROR] Preconditioned HLL solver selected. integration-variables must be set to "prec".'
      end if
    endif
    ! Transport
    if (model>0 .and. obj_transport%description=='Unavailable') &
    write(*,'(A)') '[ERROR] Transport properties are unavailable for the selected phase: cannot run Navier-Stokes simulation'
    ! Species diffusion model
    check = Schmidt_Error()
    if (len(check) > 0) obj_transport%error_message = check
    ! SAcomp: Compressibility_Correction (Lib_Spalart) takes c from c^2 f(c) = 25 Omega nu~ fv1 / a, a number against a
    ! speed, so the factor it gives the production depends on the units: the set-up refuses the model (Stop_Simulation)
    if (model == 2 .and. obj_rans%SAcomp) obj_rans%error_message = '[ERROR] turbulence = SAcomp: its '// &
      'compressibility correction solves c^2 f(c) = 25 Omega nu_t / a, a number against a speed, so the result '// &
      'depends on the units. Use SA or another turbulence model.'
    ! Turbulence sources in time-accurate runs: explicit, whatever point-implicit says. A warning recalls it when input.ini
    ! sets point-implicit to true; a run that leaves the key out takes its default, which the warning does not concern.
    ! The warning of the rotation-curvature correction (-RC) in time-accurate runs is written by Print_Warnings (Wrap_Setup)
    if (model == 2 .and. obj_time_scheme%time_accurate .and. obj_rans%point_implicit .and. &
        Key_Is_Set(trim(codename)//'-Turbulence', 'point-implicit')) &
      obj_time_scheme%warning_message = '[WARNING] point-implicit is not applied in '// &
        'time-accurate runs: the turbulence sources are explicit, and the time step includes their stiffness.'

  end subroutine Assign_Setup


  !> The check of the set-up between riemann-solver and riemann-options-Mco: a low-Mach Riemann solver (AUSM+M, the HLLC+
  !> of Chen, Low-Mach Roe) needs a cutoff Mach number above zero. Its error, or an empty string. At set-up (Assign_Setup)
  !> and at every reading of input.ini during the run, on the values read (Read_Inifile_Runtime), since a reading takes
  !> riemann-options-Mco: the run never goes on with values that its set-up refuses
  function Riemann_Mco_Error() result(out)
    use MOSE_Config_Types_m, only: obj_riemann
    implicit none
    character(len=:), allocatable :: out

    out = ''
    if (index(obj_riemann%name, 'AUSM+M')      > 0 .or. &
        index(obj_riemann%name, 'HLLC+ (Chen') > 0 .or. &
        index(obj_riemann%name, 'Low-Mach Roe') > 0) then
      if (obj_riemann%Mco == 0.0d0) then
        out = '[ERROR] Low-Mach Riemann solver selected. Mco (riemann-options-Mco) must be defined in input.'
      end if
    endif

  end function Riemann_Mco_Error


  !> The check of the set-up between schmidt and the species diffusion: Sc <= 0 (multicomponent diffusion) needs the table
  !> of the binary diffusion coefficients, INPUT/diffusion.dat, with more than one species in a viscous run (a single
  !> species does not diffuse: co_DS_expr returns D = 0 without the table). Its error, or an empty string. At set-up
  !> (Assign_Setup) and at every reading of input.ini during the run, on the values read (Read_Inifile_Runtime), since a
  !> reading takes schmidt: the run never goes on with values that its set-up refuses
  function Schmidt_Error() result(out)
    use MOSE_Config_Types_m,     only: obj_sim_param
    use MOSE_Global_m,           only: model, nsc
    use FLINT_Lib_Thermodynamic, only: dij_tab
    implicit none
    character(len=:), allocatable :: out

    out = ''
    if (model>0 .and. nsc>1 .and. obj_sim_param%Sc <= 0d0 .and. .not. allocated(dij_tab)) then
      out = '[ERROR] Sc<=0 (multicomponent diffusion) requires a binary diffusion table (INPUT/diffusion.dat)'
    end if

  end function Schmidt_Error


  !> Implicit residual smoothing is on with irs = true or irs-beta > 0: at set-up (Assign_Setup), and from the values of
  !> every reading of input.ini during the run that the run takes (Read_Inifile_Runtime), so that a file read again
  !> unchanged leaves it as it was
  logical function IRS_On()
    use MOSE_Config_Types_m, only: obj_irs
    implicit none

    IRS_On = obj_irs%requested .or. obj_irs%beta > 0d0

  end function IRS_On

end module MOSE_Assign_Setup