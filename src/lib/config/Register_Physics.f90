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
    call reg%add(section, 'equations', obj_sim_param%simulation_type, '', 'Gas dynamics equations', 'euler , navier-stokes', .true.)
    call reg%add(section, 'turbulence', obj_rans%model, 'none', 'RANS turbulence model', 'SA, SA-R, SA-RC, SAcomp, SST, Wilcox2006, none', .false.)
    call reg%add(section, 'chemistry', obj_chemistry%model, 'frozen', 'Chemistry model', 'frozen, finite-rate, equilibrium', .false.)
    call reg%add(section, 'soot-generation', obj_soot%model, 'none', 'Soot generation model', 'LL91, LIN, none', .false.)
    call reg%add(section, 'rotational-frame', obj_rot%model, 'none', 'Rotational frame model', 'rigid-body, none', .false.)

    !! ------------------------------------------------------
    !! Chemistry --------------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Chemistry'
    ! Exclude blocks from chemistry
    call reg%add(section, 'exclude-blocks', obj_chemistry%exclude_blocks_str, 'none', 'Blocks to exclude from chemistry', '', .false.)
    ! ODE solver selection
    call reg%add(section, 'ode-solver', obj_chemistry%ode_name, 'H-radau5', 'ODE solver for chemistry', 'H-radau5, sdirk4b, ros4', .false.)
    ! ODE solver parameters
    call reg%add(section, 'ode-max-steps', obj_chemistry%max_ode_steps, '100000', 'Maximum ODE integration steps', '> 0', .false.)
    call reg%add(section, 'ode-relative-tol-species', obj_chemistry%RT(1:nsc), '1e-5', 'ODE relative tolerance for species', '> 0', .false.)
    call reg%add(section, 'ode-relative-tol-temperature', obj_chemistry%RT(nsc+1), '1e-5', 'ODE relative tolerance for temperature', '> 0', .false.)
    call reg%add(section, 'ode-absolute-tol-species', obj_chemistry%AT(1:nsc), '1e-5', 'ODE absolute tolerance for species', '> 0', .false.)
    call reg%add(section, 'ode-absolute-tol-temperature', obj_chemistry%AT(nsc+1), '1e-5', 'ODE absolute tolerance for temperature', '> 0', .false.)


    !! ------------------------------------------------------
    !! Turbulence -------------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Turbulence'
    call reg%add(section, 'Prt', obj_rans%Prt, '0.85', 'Turbulent Prandtl number', '> 0', .false.)
    call reg%add(section, 'Sct', obj_rans%Sct, '0.90', 'Turbulent Schmidt number', '> 0', .false.)
    call reg%add(section, 'Sc', obj_rans%Sc, '0.7', 'Schmidt number', '> 0', .false.)
    !call reg%add(section, 'k-coupling', obj_rans%k_energy_coupling, '.false.', 'Turbulent kinetic energy coupling', 'logical', .false.)

    
    !! ------------------------------------------------------
    !! Rotating Frame ---------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Rotating-Frame'
    call reg%add(section, 'omega',  obj_rot%omega,      '0.0',         'Angular speed [rad/s]', '>= 0', .false.)
    call reg%add(section, 'axis',   obj_rot%axis_str,   '0.0 0.0 1.0', 'Rotation axis direction (3 components)', '', .false.)
    call reg%add(section, 'origin', obj_rot%origin_str, '0.0 0.0 0.0', 'Point on the rotation axis [m]', '', .false.)
    ! Stationary face entries (count scanned in Scan_Ini, strings allocated there)
    do k = 1, obj_rot%n_stationary
      call reg%add(section, 'stationary-face-'//trim(str(.true.,k)), obj_rot%stationary_face_str(k),'', 'Stationary face: block_index face_direction', '', .false.)
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