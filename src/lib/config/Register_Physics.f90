module MOSE_Read_Physics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Parameters_m
  use MOSE_Config_Types_m, only: obj_chemistry, obj_soot, obj_rans, obj_rot, obj_sim_param
  use MOSE_Input_Registry

  implicit none
  private
  public :: Register_Physics

contains

  subroutine Register_Physics()
    use MOSE_Global_m
    use IR_Precision, only: str
    implicit none
    integer :: k
    character(len=:), allocatable :: section

    !! ------------------------------------------------------
    !! Physical models ------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Physics'
    call reg%add(section, 'equations', obj_sim_param%simulation_type, '', 'Gas dynamics equations', 'euler , navier-stokes', .true., KEY_FIXED)
    call reg%add(section, 'turbulence', obj_rans%model, 'none', 'RANS turbulence model', 'SA, SA-R, SA-RC, SAcomp, SST, SST-RC, Wilcox2006, none', .false., KEY_FIXED)
    call reg%add(section, 'chemistry', obj_chemistry%model, 'frozen', 'Chemistry model', 'frozen, finite-rate, equilibrium', .false., KEY_FIXED)
    call reg%add(section, 'soot-generation', obj_soot%model, 'none', 'Soot generation model', 'LL91, LIN, none', .false., KEY_FIXED)
    call reg%add(section, 'rotational-frame', obj_rot%model, 'none', 'Rotational frame model', 'rigid-body, none', .false., KEY_FIXED_OWN)
    ! Laminar (molecular) transport closure
    call reg%add(section, 'schmidt', obj_sim_param%Sc, '0.0', 'Laminar Schmidt number', '>= 0', .false., KEY_TAKEN)
    call reg%add(section, 'prandtl', obj_sim_param%Prl, '0.0', 'Laminar Prandtl number', '>= 0', .false., KEY_TAKEN)

    !! ------------------------------------------------------
    !! Chemistry --------------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Chemistry'
    ! Exclude blocks from chemistry
    call reg%add(section, 'exclude-blocks', obj_chemistry%exclude_blocks_str, 'none', &
                 'Blocks to exclude from chemistry (finite-rate and equilibrium, on every multigrid level)', &
                 '', .false., KEY_FIXED)
    ! ODE solver selection
    ! Every integrator OSLO's setup_odesolver knows. cvode needs OSLO built with
    ! USE_SUNDIALS; H-dopri5 is explicit and only suited to non-stiff mechanisms.
    call reg%add(section, 'ode-solver', obj_chemistry%ode_name, 'H-radau5', 'ODE solver for chemistry', &
                 'H-radau5, H-rodas, H-sdirk4, H-dopri5, '// &
                 'radau2a, lobatto3c, gauss, radau1a, '// &
                 'ros2, ros3, ros4, rodas3, rodas4, '// &
                 'sdirk2a, sdirk2b, sdirk3a, sdirk4a, sdirk4b, '// &
                 'dodesol, cvode', .false., KEY_FIXED)
    ! ODE solver parameters
    call reg%add(section, 'ode-max-steps', obj_chemistry%max_ode_steps, '100000', 'Maximum ODE integration steps', '> 0', .false., KEY_FIXED)
    ! Only ONERA-7 and Frolov_nopressure ship an analytical Jacobian, and only
    ! the Hairer solvers can use one; anything else falls back to finite
    ! differences.
    call reg%add(section, 'ode-analytical-jacobian', obj_chemistry%analytical_jacobian, '.false.', 'Use the mechanism analytical Jacobian instead of finite differences', 'logical', .false., KEY_FIXED)
    call reg%add(section, 'ode-relative-tol-species', obj_chemistry%RT(1:nsc), '1e-5', 'ODE relative tolerance for species', '> 0', .false., KEY_FIXED)
    call reg%add(section, 'ode-relative-tol-temperature', obj_chemistry%RT(nsc+1), '1e-5', 'ODE relative tolerance for temperature', '> 0', .false., KEY_FIXED)
    call reg%add(section, 'ode-absolute-tol-species', obj_chemistry%AT(1:nsc), '1e-5', 'ODE absolute tolerance for species', '> 0', .false., KEY_FIXED)
    call reg%add(section, 'ode-absolute-tol-temperature', obj_chemistry%AT(nsc+1), '1e-5', 'ODE absolute tolerance for temperature', '> 0', .false., KEY_FIXED)


    !! ------------------------------------------------------
    !! Turbulence -------------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Turbulence'
    call reg%add(section, 'Prt', obj_rans%Prt, '0.85', 'Turbulent Prandtl number', '> 0', .false., KEY_TAKEN)
    call reg%add(section, 'Sct', obj_rans%Sct, '0.90', 'Turbulent Schmidt number', '> 0', .false., KEY_TAKEN)
    call reg%add(section, 'point-implicit', obj_rans%point_implicit, '.true.', 'Point-implicit factor 1/(1 + dt d) on the increment of the turbulence equations, steady runs only (time-accurate runs: explicit sources, stiffness in the time step)', 'logical', .false., KEY_TAKEN)
    call reg%add(section, 'omega-wall-bc', obj_rans%omega_wall_bc, 'practical', 'Omega wall condition for k-omega models: Menter practical (60*nu/(beta1*y^2), beta1 of SST, for both models) or asymptotic (6*nu/(beta*y^2), exact y->0 limit with the beta of the model, beta1 for SST and beta0 for Wilcox2006, for wall-resolved meshes)', 'practical, asymptotic', .false., KEY_FIXED)
    call reg%add(section, 'sst-production', obj_rans%sst_production, 'incompressible', 'SST production: incompressible (P = mu_t*S^2; NASA TMR SST-2003m, or SST-2003s with k-energy-coupling) or compressible (P = tau_ij*du_i/dx_j as OpenFOAM kOmegaSST arranges it: the limiter on mu_t*(S^2 - 2/3*divU^2) only, then -2/3*rho*k*divU added, and -2/3*gamma*rho*omega*divU in the omega equation; TMR SST-2003e limits the whole P and takes gamma/nu_t times it in the omega equation)', 'incompressible, compressible', .false., KEY_FIXED)
    call reg%add(section, 'k-energy-coupling', obj_rans%k_energy_coupling, '.false.', 'Two-equation models (SST, Wilcox2006): k in the total energy, the full compressible RANS form of NASA TMR: -2/3*rho*k*delta_ij in the stress tensor of the momentum and energy equations, and the k-equation source (production - dissipation) taken from the mean-flow energy. .false. = TMR m/e variants', 'logical', .false., KEY_FIXED)

    
    !! ------------------------------------------------------
    !! Rotating Frame ---------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Rotating-Frame'
    call reg%add(section, 'omega',  obj_rot%omega,      '0.0',         'Angular speed [rad/s]', '>= 0', .false., KEY_TAKEN)
    call reg%add(section, 'axis',   obj_rot%axis_str,   '0.0 0.0 1.0', 'Rotation axis direction (3 components)', '', .false., KEY_TAKEN)
    call reg%add(section, 'origin', obj_rot%origin_str, '0.0 0.0 0.0', 'Point on the rotation axis [m]', '', .false., KEY_TAKEN)
    ! Stationary face entries (count scanned in Scan_Ini, strings allocated there)
    do k = 1, obj_rot%n_stationary
      call reg%add(section, 'stationary-face-'//trim(str(.true.,k)), obj_rot%stationary_face_str(k),'', 'Stationary face: block_index face_direction', '', .false., KEY_FIXED_OWN)
    end do

    !! ------------------------------------------------------
    !! GSI --------------------------------------------------
    !! ------------------------------------------------------
    !! TODO: Refactor GSI material properties to registry
    !! For now: HDPE, PP, Paraffin, HTPB properties
    !! These require dedicated sections in config types
    !! ======================================================

  end subroutine Register_Physics

end module MOSE_Read_Physics