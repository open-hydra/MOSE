module MOSE_Read_Numerics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Config_Types_m
  use MOSE_Input_Registry
  
  implicit none
  private
  public :: Register_Numerics

contains

  subroutine Register_Numerics (nmgl)
    use MOSE_Parameters_m
    use IR_precision
    implicit none
    ! Local
    integer, intent(in) :: nmgl
    character(len=llen) :: section

    section = trim(codename)//'-Numerics'

    !! ------------------------------------------------------
    !! Time Scheme ------------------------------------------
    !! ------------------------------------------------------
    obj_time_scheme%warning_message = 'none'
    obj_time_scheme%error_message   = 'none'
    obj_time_scheme%description     = 'none'

    ! Solver-type
    call reg%add( trim(section), 'time-scheme', obj_time_scheme%solver_type, 'euler', 'Time integration solver', 'euler, RK2, RK3', .true., KEY_FIXED )

    ! Stability coefficients and related options
    call reg%add( trim(section), 'cfl', obj_time_scheme%cfl, '0.5', &
                  'CFL number of the convective spectral radii summed over the directions (RK2 <= 1, RK3 <= 1.17; at most 1 '// &
                  'with the multirate integrator on ([MOSE-Multirate] enabled = .true., max-level > 0 and nonlocal-cfl = '// &
                  'redo, levels = recompute and no fixed-macro-step))', '> 0', .true., KEY_TAKEN )
    call reg%add( trim(section), 'vnn', obj_time_scheme%vnn, '0.3', &
                  'Diffusion number of the largest diffusivity summed over the directions (RK2 <= 0.5, RK3 <= 0.62; at most '// &
                  '0.5 in a Navier-Stokes run with the multirate integrator on ([MOSE-Multirate] enabled = .true., '// &
                  'max-level > 0 and nonlocal-cfl = redo, levels = recompute and no fixed-macro-step))', &
                  '> 0', .false., KEY_TAKEN )
    call reg%add( trim(section), 'cfl-rise-threshold', obj_time_scheme%rampa_cfl_iter, '0', 'CFL rise threshold', '>= 0', .false., KEY_TAKEN )

    ! Time-accurate switch
    call reg%add( trim(section), 'time-accurate', obj_time_scheme%time_accurate, '.false.', 'Time accurate switch', 'logical', .true., KEY_FIXED )

    ! New state
    call reg%add( trim(section), 'integration-variables', obj_time_scheme%integration_variables, 'cons', 'Integration variables (cons/prim/prec)', 'cons, prim, prec', .false., KEY_FIXED )


    ! Implicit residual smoothing --------------------------
    obj_irs%description     = 'none'
    obj_irs%warning_message = 'none'
    obj_irs%error_message   = 'none'
    call reg%add( trim(section), 'irs', obj_irs%requested, '.false.', &
          'Implicit Residual Smoothing (steady runs only; a run started without it cannot switch it on)', 'logical', .false., KEY_TAKEN )
    call reg%add( trim(section), 'irs-beta', obj_irs%beta, '0.0', &
          'IRS beta parameter (> 0 also enables IRS; steady runs only)', '>= 0', .false., KEY_TAKEN )
    call reg%add( trim(section), 'irs-variables', obj_irs%variables, 'flow', &
          'Residuals smoothed by IRS: flow (species, momentum, energy) or all (also the RANS, '// &
          'soot and passive-scalar equations; for Spalart-Allmaras: the k-omega models diverge with all)', 'flow, all', .false., KEY_TAKEN )

    ! Preconditioning controls (used when integration-variables = prec)
    obj_prec%description     = 'none'
    obj_prec%warning_message = 'none'
    obj_prec%error_message   = 'none'
    obj_prec%enabled         = .false.
    call reg%add( trim(section), 'preconditioning-Uref', obj_prec%Ur_ref, '0.0', 'Reference velocity for preconditioning', '>= 0', .false., KEY_TAKEN )
    call reg%add( trim(section), 'preconditioning-Mach', obj_prec%Mach_target, '0.0', 'Mach target for preconditioning', '>= 0', .false., KEY_TAKEN )
    call reg%add( trim(section), 'preconditioning-eps-min', obj_prec%eps_min, '0.05', 'Low-Mach cutoff for Ur', 'in (0, 1)', .false., KEY_TAKEN )
    call reg%add( trim(section), 'preconditioning-Ur-min', obj_prec%Ur_min, '0.0', 'Floor on Ur [m/s]', '>= 0', .false., KEY_TAKEN )
    call reg%add( trim(section), 'preconditioning-Ur-factor', obj_prec%Ur_factor, '1.0', 'Weiss-Smith multiplicative factor', '> 0', .false., KEY_TAKEN )
    call reg%add( trim(section), 'preconditioning-Ur-smooth', obj_prec%n_smooth_Ur, '0', 'Ur max-smoothing passes', '>= 0', .false., KEY_TAKEN )

    !! ------------------------------------------------------
    !! ------------------------------------------------------


    !! ------------------------------------------------------
    !! Space Scheme ------------------------------------------
    !! ------------------------------------------------------

    ! Space Discretization ---------------------------------
    obj_space_scheme%description     = 'none'
    obj_space_scheme%warning_message = 'none'
    obj_space_scheme%error_message   = 'none'
    call reg%add( trim(section), 'space-reconstruction', obj_space_scheme%space_reconstruction, '', 'Space reconstruction method', 'MUSCL, first-order', .true., KEY_FIXED )
    call reg%add( trim(section), 'flux-limiter', obj_space_scheme%flux_limiter, '', &
          'Flux limiter for MUSCL reconstruction (van Leer if omitted; ignored with first-order)', &
          'vanalbada, minmod, superbee, vanleer, mc', .false., KEY_FIXED )
    call reg%add( trim(section), 'diffusive-metric', obj_space_scheme%diffusive_metric, 'stencil', &
                  'Face metric of the viscous fluxes: stencil takes the face Jacobian from the same differences of '// &
                  'the cell centres as the face gradient, exact for a linear field on any mesh (at the symmetry faces and the '// &
                  'k faces of type 200 for the linear fields that they allow) except at the chimera faces (102), which take '// &
                  'the inverse of the mean of the inverse metrics of the two cells with every value, and at the symmetry '// &
                  'faces with one cell across them, the walls 502 and the walls 503-505 while they ablate, which take the '// &
                  'metric of the cell; mean takes the mean of the two metrics (the faces between blocks and of type 200 the '// &
                  'inverse of the mean of the inverse metrics), overstates the gradient by (1+r)^2/(4r) at a cell-size ratio '// &
                  'r and does not converge next to a symmetry plane (300) that bounds a sector of a polar mesh', &
                  'stencil, mean', .false., KEY_FIXED )

    ! Shock Detector ---------------------------------------
    obj_shock_detector%warning_message = 'none'
    obj_shock_detector%error_message   = 'none'
    call reg%add( trim(section), 'shock-detector', obj_shock_detector%description, '', 'Shock detector method', 'Tramel, Chen', .false., KEY_FIXED )

    ! Riemann solver ---------------------------------------
    obj_riemann%description     = 'none'
    obj_riemann%warning_message = 'none'
    obj_riemann%error_message   = 'none'
    call reg%add( trim(section), 'riemann-solver', obj_riemann%description, 'HLLC', &
          'Riemann solver', 'HLLC, HLLC+Tramel, HLLC+Chen, HLLC-PC, HLLE, HLLE++, SLAU, SLAU2, LMRoe, MiczekRoe, AUSM+, AUSM+M, LLF, Rusanov, exact', .false., KEY_FIXED )
    call reg%add( trim(section), 'riemann-options-Mco', obj_riemann%Mco, '0.0', 'Low-Mach acoustic-dissipation cutoff Mach (floor)', '>= 0', .false., KEY_TAKEN )


    ! Multigrid levels --------------------------------------
    call Register_Multigrid_Levels(nmgl)

    ! Multirate (local time stepping) ----------------------
    call Register_Multirate()

    !! ------------------------------------------------------
    !! ------------------------------------------------------

  end subroutine Register_Numerics


  subroutine Register_Multirate()
    !> Keys of the multirate integrator (MOSE_Mod_Multirate), section [MOSE-Multirate].
    use MOSE_Config_Types_m
    use MOSE_Parameters_m
    implicit none
    character(len=llen) :: section

    section = trim(codename)//'-Multirate'

    obj_multirate%warning_message = 'none'
    obj_multirate%error_message   = 'none'
    obj_multirate%description     = 'none'

    call reg%add( trim(section), 'enabled', obj_multirate%enabled, '.false.', 'Multirate local time stepping: recursive &
      &Heun (SSP-RK2) with a level per cell, time-accurate runs only', 'logical', .false., KEY_FIXED )
    call reg%add( trim(section), 'max-level', obj_multirate%max_level, '10', &
      'Deepest level allowed (at most 50): the map of the adaptive step takes L = max-level levels below level 0, or fewer '// &
      'when the local steps span fewer octaves or macro-step-target gives a shorter H, so that H stays within the step of '// &
      'the slowest cells divided by safety; the finest cells advance with H/2^L', 'in [0, 50]', .false., KEY_FIXED )
    call reg%add( trim(section), 'safety', obj_multirate%safety, '2.0', &
      'Safety factor sigma: every cell advances with a sub-step h <= dtcap/sigma (dtcap = local CFL/VNN step); '// &
      'the default 2 keeps every cell, at the start of a macro step, within the admissible step of each face neighbour '// &
      'without the interface margin', &
      '>= 1', .false., KEY_FIXED )
    call reg%add( trim(section), 'macro-step-target', obj_multirate%macro_step_target, '0.0', &
      'Target macro step H [s]; 0 = none (H set by the level map, within the step of the slowest cells divided by safety)', &
      '>= 0', .false., KEY_FIXED )
    call reg%add( trim(section), 'levels', obj_multirate%levels, 'recompute', &
      'Level map: computed at the first macro step and kept (fixed) or recomputed from the local steps (recompute)', &
      'fixed, recompute', .false., KEY_FIXED )
    call reg%add( trim(section), 'recompute-every', obj_multirate%recompute_every, '1', &
      'Macro steps between two recomputations of the level map (levels = recompute)', '> 0', .false., KEY_FIXED )
    call reg%add( trim(section), 'fixed-macro-step', obj_multirate%fixed_macro_step, '0.0', 'Fixed macro step as a multiple &
      &of the dt_min of the first macro step; 0 = adaptive (verification runs)', '>= 0', .false., KEY_FIXED )
    call reg%add( trim(section), 'conservative-correction', obj_multirate%conservative_correction, '.true.', 'Conservative &
      &flux correction at the faces between levels (.false. only as a positive control)', 'logical', .false., KEY_FIXED )
    call reg%add( trim(section), 'trace', obj_multirate%trace, '.false.', 'Print the stage sequence, the level map and the &
      &halo lists of the first macro step (standard output)', 'logical', .false., KEY_FIXED )
    call reg%add( trim(section), 'nonlocal-cfl', obj_multirate%nonlocal_cfl, 'redo', &
      'Condition (4) of Gnedin et al. checked at every sub-step start: off; count (counts only); redo (a macro step that '// &
      'fails it is redone from its start with finer levels or a smaller H, with the interface margin; a macro step in '// &
      'which a cell or a face takes a state that is not physical is redone with H halved (with more than one level))', &
      'off, count, redo', .false., KEY_FIXED )
    call reg%add( trim(section), 'interface-margin', obj_multirate%interface_margin, '.false.', &
      'Faces between levels: the coarse step within dtcap/kappa of the fine cell, kappa = max(cfl (1 + theta) of the '// &
      'limiter, 1.1) (always on with nonlocal-cfl = redo)', 'logical', .false., KEY_FIXED )
    call reg%add( trim(section), 'force-redo', obj_multirate%force_redo, '.false.', &
      'Test device of nonlocal-cfl = redo: the first attempt of every macro step is rejected and redone with the same '// &
      'macro step and level map, so the results must not change', 'logical', .false., KEY_FIXED )
    call reg%add( trim(section), 'max-attempts', obj_multirate%max_attempts, '21', &
      'nonlocal-cfl = redo: attempts of a macro step at most (an attempt rejected by force-redo is not counted); when they '// &
      'are all rejected the run stops with the rejections by cause', '>= 1', .false., KEY_FIXED )

  end subroutine Register_Multirate


  subroutine Register_Multigrid_Levels(nmgl)
    use MOSE_Config_Types_m
    use MOSE_Parameters_m
    use IR_precision
    implicit none
    integer, intent(in) :: nmgl
    integer :: m
    character(len=llen) :: option

    ! Allocate iteration threshold array for each level
    obj_multigrid%MGL = max(nmgl, 1)
    allocate( obj_multigrid%iter_threshold(obj_multigrid%MGL) )
    obj_multigrid%iter_threshold = 0

    ! Register per-level iteration parameters (level-2-iter, level-3-iter, ...)
    if (.not. obj_sim_param%HYDRA_MG) then 
      do m = 1, obj_multigrid%MGL
        write(option,'(A5,I0,A5)') 'level', m, '-iter'
        call reg%add( trim(codename)//'-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m)), '>= 0', .false., KEY_TAKEN )
      enddo
      call reg%add( trim(codename)//'-Multigrid', 'prolongated-diter', obj_multigrid%prolongated_diter, '0', &
        'Coarse levels: least level iterations between two writes of field/wall-prolongated, which are always written at the end of a level (0 = with every field write)', '>= 0', .false., KEY_TAKEN )
    else
      do m = 1, obj_multigrid%MGL
        write(option,'(A5,I0,A5)') 'level', m, '-iter'
        call reg%add( 'HYDRA-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m)), '>= 0', .false., KEY_TAKEN )
      enddo
      call reg%add( 'HYDRA-Multigrid', 'prolongated-diter', obj_multigrid%prolongated_diter, '0', &
        'Coarse levels: least level iterations between two writes of field/wall-prolongated, which are always written at the end of a level (0 = with every field write)', '>= 0', .false., KEY_TAKEN )
    endif

  end subroutine Register_Multigrid_Levels  

end module MOSE_Read_Numerics