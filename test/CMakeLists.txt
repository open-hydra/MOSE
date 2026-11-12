# ---------------------------------------------------------------------------
#  CTest definitions of MOSE: one ctest per case, check or part of a check.
#  Tiers, kinds of check and the list of tests: test/TESTING-PLAN.md.
#
#  Every test is registered with mose_test() below: a name, a one-line DESC, the
#  KIND of check and the area LABELS. Its tier is not written here: it comes from
#  test/tiers.txt (one line per test), because the tier follows from the measured
#  time of the test, not from what it checks:
#
#    light    under 5 s     pre-push hook (ctest -L light) and every PR
#    medium   5 s to 100 s  every PR
#    heavy    100 s and up  local, on request (ctest -L heavy)
#
#  Labels of a test: its tier; its kind (exact order identity invariant robustness
#  contract unit reference); its area (physics numerics implementation unit); the
#  dimensions it runs (1D 2D 3D); multirate; and, where it needs a build option,
#  needs-mpi. PR CI: ctest -L "light|medium" -LE needs-mpi, and in a second job
#  built with USE_MPI=ON, ctest -L needs-mpi -LE heavy.
#
#  At configure the registration writes <build>/test/inventory.tsv, from which
#  test/tools/inventory.py renders the tables of the plan.
# ---------------------------------------------------------------------------

# Common environment for all solver tests
set(_MOSE_TEST_ENV
    "OMP_NUM_THREADS=1"
    "KMP_STACKSIZE=100M"
)
# The restart tests also write and read the Tecplot binary format (.szplt), which ORION writes with TecIO
if(USE_TECIO)
    list(APPEND _MOSE_TEST_ENV "MOSE_TECIO=1")
endif()

# ---------------------------------------------------------------------------
#  Tiers: test/tiers.txt, "<Name> <light|medium|heavy> <seconds>" per line; the seconds are the time on the CI runner
#  (tools/ctest_tiers.py keeps them). A change of the file reconfigures.
# ---------------------------------------------------------------------------
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tiers.txt")
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/tiers.txt" _tier_lines)
set(_tier_names "")
foreach(_line IN LISTS _tier_lines)
    if(_line MATCHES "^[ \t]*(#|$)")
        continue()
    endif()
    if(NOT _line MATCHES "^([^ \t]+)[ \t]+(light|medium|heavy)[ \t]+([0-9.]+|-)[ \t]*$")
        message(FATAL_ERROR "test/tiers.txt: cannot read the line '${_line}'")
    endif()
    set(MOSE_TIER_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    set(MOSE_SECONDS_${CMAKE_MATCH_1} ${CMAKE_MATCH_3})
    list(APPEND _tier_names ${CMAKE_MATCH_1})
endforeach()

# ---------------------------------------------------------------------------
#  mose_test(<Name>
#      DESC "<what it checks, one line, no semicolon>"
#      KIND <exact|order|identity|invariant|robustness|contract|unit|reference>
#      LABELS <area and dimension labels>
#      [DIR <folder of test/>] [RUN "<shell command>"]   a script run in DIR (default ./run.sh)
#      [COMMAND <program> <arguments>]                   a program of the build instead (unit checks)
#      [TIMEOUT <s>] [LOCK <name>] [PROCESSORS <n>] [ENV <VAR=value>...]
#      [REQUIRES <MPI|LINUX>...])                         built only with USE_MPI, only on Linux
#  A test that is not built is still listed in the inventory. Two tests that share a folder have the same LOCK.
# ---------------------------------------------------------------------------
function(mose_test NAME)
    cmake_parse_arguments(T "" "DIR;RUN;KIND;DESC;TIMEOUT;LOCK;PROCESSORS" "LABELS;COMMAND;ENV;REQUIRES" ${ARGN})
    if(NOT T_DESC OR NOT T_KIND OR NOT T_LABELS)
        message(FATAL_ERROR "mose_test(${NAME}): DESC, KIND and LABELS are required")
    endif()
    if(NOT T_KIND MATCHES "^(exact|order|identity|invariant|robustness|contract|unit|reference)$")
        message(FATAL_ERROR "mose_test(${NAME}): unknown KIND '${T_KIND}'")
    endif()
    if(NOT DEFINED MOSE_TIER_${NAME})
        message(FATAL_ERROR "mose_test(${NAME}): no line for it in test/tiers.txt")
    endif()
    set(_tier ${MOSE_TIER_${NAME}})
    if(NOT T_TIMEOUT)
        if(_tier STREQUAL "light")
            set(T_TIMEOUT 60)
        elseif(_tier STREQUAL "medium")
            set(T_TIMEOUT 600)
        else()
            set(T_TIMEOUT 3600)
        endif()
    endif()
    set(_labels ${_tier} ${T_KIND} ${T_LABELS})
    set(_env ${_MOSE_TEST_ENV} ${T_ENV})
    set(_build TRUE)
    foreach(_r IN LISTS T_REQUIRES)
        if(_r STREQUAL "MPI")
            list(APPEND _labels needs-mpi)
            list(APPEND _env "MOSE_MPIEXEC=${MPIEXEC_EXECUTABLE}" "MOSE_MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}")
            if(NOT T_PROCESSORS)
                set(T_PROCESSORS 4)
            endif()
            if(NOT USE_MPI)
                set(_build FALSE)
            endif()
        elseif(_r STREQUAL "LINUX")
            if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
                set(_build FALSE)
            endif()
        else()
            message(FATAL_ERROR "mose_test(${NAME}): unknown REQUIRES '${_r}'")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _labels)

    set_property(GLOBAL APPEND PROPERTY MOSE_TEST_NAMES ${NAME})
    set_property(GLOBAL PROPERTY MOSE_TEST_${NAME}_KIND "${T_KIND}")
    set_property(GLOBAL PROPERTY MOSE_TEST_${NAME}_DESC "${T_DESC}")
    set_property(GLOBAL PROPERTY MOSE_TEST_${NAME}_LABELS "${_labels}")
    set_property(GLOBAL PROPERTY MOSE_TEST_${NAME}_BUILT "${_build}")
    if(NOT _build)
        return()
    endif()

    if(T_COMMAND)
        add_test(NAME ${NAME} COMMAND ${T_COMMAND})
    else()
        if(NOT T_RUN)
            set(T_RUN "./run.sh")
        endif()
        add_test(NAME ${NAME} COMMAND bash -c "${T_RUN}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/${T_DIR}")
    endif()
    set_tests_properties(${NAME} PROPERTIES LABELS "${_labels}" TIMEOUT ${T_TIMEOUT} ENVIRONMENT "${_env}")
    if(T_LOCK)
        set_tests_properties(${NAME} PROPERTIES RESOURCE_LOCK ${T_LOCK})
    endif()
    if(T_PROCESSORS)
        set_tests_properties(${NAME} PROPERTIES PROCESSORS ${T_PROCESSORS})
    endif()
endfunction()

# ---------------------------------------------------------------------------
#  Reference cases: 1-D, 2-D and 3-D problems of test/1D, test/2D, test/3D (each folder: MOSE.sh and a verifier)
# ---------------------------------------------------------------------------
mose_test(Sod79
    DIR 1D/Sod79  RUN "./MOSE.sh solve && python3 -B check.py"
    KIND exact  LABELS numerics 1D
    DESC "Sod tube against the exact solution")
mose_test(Einfeldt91
    DIR 1D/Ein91  RUN "./MOSE.sh test && python3 -B verify.py"
    KIND exact  LABELS numerics 1D
    DESC "Einfeldt 1-2-0 tube against the exact solution (ExactPack)")
mose_test(Noh87
    DIR 1D/Noh87  RUN "./MOSE.sh test && python3 -B verify.py"
    KIND exact  LABELS numerics 1D
    DESC "Noh implosion against the exact solution (ExactPack)")
mose_test(Toro99
    DIR 1D/Tor99  RUN "./MOSE.sh test && python3 -B verify.py"
    KIND exact  LABELS numerics 1D
    DESC "Toro tubes against the exact solution (ExactPack)")
mose_test(ShuOsher
    DIR 1D/Shu-Osher  RUN "./MOSE.sh solve && python3 -B verify.py"
    KIND reference  LABELS numerics 1D
    LOCK Shu-Osher
    DESC "Shu-Osher problem, N=200 against the N=1600 reference")
mose_test(ShuOsherConvergence
    DIR 1D/Shu-Osher  RUN "./MOSE.sh test && python3 -B convergence.py"
    KIND order  LABELS numerics 1D
    LOCK Shu-Osher
    DESC "Shu-Osher on N=200, 400, 800: monotone convergence")
mose_test(PrandtlMeyer
    DIR 2D/euler/prandtl-meyer  RUN "./MOSE.sh solve && python3 -B check.py"
    KIND exact  LABELS numerics 2D
    DESC "Prandtl-Meyer expansion against the analytic solution")
mose_test(IsentropicVortex
    DIR 2D/euler/isentropic-vortex  RUN "./MOSE.sh test && python3 -B verify.py"
    KIND order  LABELS numerics 2D
    DESC "isentropic vortex on 32, 64, 128 cells: order of accuracy")
mose_test(MMSViscous
    DIR 2D/viscous/mms  RUN "./MOSE.sh testfast && python3 -B verify.py"
    KIND order  LABELS numerics 2D
    DESC "viscous operator, manufactured solution: order of accuracy")
mose_test(MMSViscousSkewed
    DIR 2D/viscous/mms/skewed  RUN "./MOSE.sh testfast && python3 -B verify.py"
    KIND order  LABELS numerics 2D 3D
    DESC "manufactured solutions on skewed two-block grids: order of accuracy")
mose_test(MMSViscous3D
    DIR 3D/viscous/mms  RUN "./MOSE.sh testfast && python3 -B verify.py"
    KIND order  LABELS numerics 3D
    LOCK mms3d
    DESC "3-D viscous manufactured solutions on two-block grids: order of accuracy")
mose_test(MMSViscous3DFull
    DIR 3D/viscous/mms  RUN "./MOSE.sh test && python3 -B verify.py"
    KIND order  LABELS numerics 3D
    LOCK mms3d  TIMEOUT 259200
    DESC "full form of MMSViscous3D, the runs of the V&V page (18 h)")
mose_test(PlanarNozzleGrids
    DIR 2D/viscous/planar-nozzle  RUN "./MOSE.sh testfast && python3 -B verify.py"
    KIND invariant  LABELS numerics 2D
    DESC "planar nozzle on O and H grids stays mirror-symmetric")

# ---------------------------------------------------------------------------
#  TMR turbulent flat plate: one test per model (restart from the converged 35x25 solution), the criteria of the V&V page, and the full form
# ---------------------------------------------------------------------------
mose_test(FlatPlateTMR.SA
    DIR 2D/viscous/flat-plate-turbulent-tmr  RUN "./MOSE.sh testfast SA && python3 -B verify.py --short SA"
    KIND reference  LABELS physics 2D
    LOCK flat-plate-tmr
    DESC "TMR flat plate, Spalart-Allmaras: restart from the converged 35x25 solution keeps Cf and CD")
mose_test(FlatPlateTMR.SST
    DIR 2D/viscous/flat-plate-turbulent-tmr  RUN "./MOSE.sh testfast SST && python3 -B verify.py --short SST"
    KIND reference  LABELS physics 2D
    LOCK flat-plate-tmr
    DESC "TMR flat plate, SST: restart from the converged 35x25 solution keeps Cf and CD")
mose_test(FlatPlateTMR.Wilcox2006
    DIR 2D/viscous/flat-plate-turbulent-tmr  RUN "./MOSE.sh testfast Wilcox2006 && python3 -B verify.py --short Wilcox2006"
    KIND reference  LABELS physics 2D
    LOCK flat-plate-tmr
    DESC "TMR flat plate, Wilcox 2006: restart from the converged 35x25 solution keeps Cf and CD")
mose_test(FlatPlateTMR.SSGLRR
    DIR 2D/viscous/flat-plate-turbulent-tmr  RUN "./MOSE.sh testfast SSGLRR && python3 -B verify.py --short SSGLRR"
    KIND reference  LABELS physics 2D
    LOCK flat-plate-tmr
    DESC "TMR flat plate, SSG/LRR-omega: restart from the converged 35x25 solution keeps Cf and CD")
mose_test(FlatPlateTMR.Reference
    DIR 2D/viscous/flat-plate-turbulent-tmr  RUN "python3 -B verify.py --reference"
    KIND reference  LABELS physics 2D
    LOCK flat-plate-tmr
    DESC "criteria of the V&V page from the stored results (CFL3D, FUN3D, TAU)")
mose_test(FlatPlateTMRFull
    DIR 2D/viscous/flat-plate-turbulent-tmr  RUN "./MOSE.sh test && python3 -B verify.py"
    KIND reference  LABELS physics 2D
    LOCK flat-plate-tmr  TIMEOUT 370800
    DESC "full form of the TMR flat plate: 4 models, 3 grids, 300000 iterations each (days)")

# ---------------------------------------------------------------------------
#  Cases with their own run.sh in test/2D and test/3D
# ---------------------------------------------------------------------------
mose_test(MMSResiduals
    DIR 2D/viscous/mms/residuals
    KIND exact  LABELS physics 2D 3D
    DESC "RANS right-hand sides on manufactured fields (SA, SST, Wilcox 2006)")
mose_test(PointImplicitSteady
    DIR 2D/viscous/point-implicit-steady
    KIND invariant  LABELS numerics 2D
    DESC "steady state with the point-implicit factor does not depend on the time step")
mose_test(WedgeAxisStep
    DIR 2D/viscous/wedge-axis-step
    KIND robustness  LABELS numerics 2D
    DESC "local time step next to the axis of a wedge is stable at the largest vnn")
mose_test(Preconditioning
    DIR 2D/euler/low-mach-cylinder
    KIND exact  LABELS numerics 2D
    DESC "low-Mach cylinder, Mach 0.01 and 0.001: Cp Mach-independent, near potential flow")
mose_test(AblatingWall
    DIR 2D/viscous/ablating-wall
    KIND exact  LABELS physics 2D
    DESC "ablation 503, 504 and SRM grain 502 against the closed form")
mose_test(RSMDecay
    DIR 2D/viscous/rsm-decay
    KIND exact  LABELS physics 2D
    DESC "SSG/LRR decay of isotropic turbulence against the closed form")

# ---------------------------------------------------------------------------
#  test/checks: small cases that check one mechanism; a check with parts has one test per part (Area.Part)
# ---------------------------------------------------------------------------
mose_test(RiemannSmoke
    DIR checks/riemann-smoke
    KIND robustness  LABELS numerics 1D
    DESC "every Riemann solver stays positive, shock sensor on")
mose_test(NumericsSmoke
    DIR checks/numerics-smoke
    KIND robustness  LABELS numerics 1D
    DESC "limiter, scheme, reconstruction and variable sweep stays positive")
mose_test(Conservation
    DIR checks/conservation
    KIND invariant  LABELS implementation 1D
    DESC "closed slip-wall box conserves mass to round-off")
mose_test(OpenMPEquiv
    DIR checks/openmp-equiv
    KIND identity  LABELS implementation 1D
    DESC "1 against 4 threads, bit for bit: Sod tube and a chimera case")
mose_test(Restart
    DIR checks/restart
    KIND identity  LABELS implementation 1D
    DESC "split run plus restart equals the straight run, every format, probes")
mose_test(Restart.Counters
    DIR checks/restart/restart-counters
    KIND identity  LABELS implementation 1D
    DESC "a restart carries every counter: numbered outputs, histories, probes")
mose_test(TimeIntegration.TimeControl
    DIR checks/time-integration/time-control
    KIND order  LABELS numerics 1D
    DESC "second order in time with a time-dependent inlet, time tables, stop rules")
mose_test(TimeIntegration.LocalTimeStep
    DIR checks/time-integration/local-time-step
    KIND robustness  LABELS numerics 2D 3D
    DESC "local time step stable on square, skewed and cubic cells")
mose_test(TimeIntegration.TurbulenceSources
    DIR checks/time-integration/turbulence-sources
    KIND order  LABELS numerics 1D
    DESC "turbulence sources in time-accurate runs: order and positivity at a wall")
mose_test(TimeIntegration.Quasi2DSymmetry
    DIR checks/time-integration/quasi2d-symmetry
    KIND invariant  LABELS numerics 1D 2D
    DESC "symmetry faces of quasi-2-D and 1-D cases: round-off in w does not grow")
mose_test(BoundaryConditions.Inlets
    DIR checks/boundary-conditions/inlets
    KIND exact  LABELS physics 1D
    DESC "inlets at total conditions, mass flux, normal velocity: fluxes against the face state")
mose_test(BoundaryConditions.Nozzle
    DIR checks/boundary-conditions/nozzle
    KIND exact  LABELS physics 1D
    DESC "choked nozzle 420: exit state below and above the shock pressures")
mose_test(BoundaryConditions.BCInput
    DIR checks/boundary-conditions/bc-input
    KIND contract  LABELS implementation 1D
    DESC "property lines of bc.txt are read as written or refused")
mose_test(BoundaryConditions.RotatingFrame
    DIR checks/boundary-conditions/rotating-frame
    KIND contract  LABELS implementation 1D
    DESC "stationary faces of the rotating frame are read and listed")
mose_test(BoundaryConditions.WallHeatSign
    DIR checks/boundary-conditions/wall-heat-sign
    KIND exact  LABELS physics 1D
    DESC "sign of the wall heat flux and of the wall stress")
mose_test(BoundaryConditions.Extrapolation400
    DIR checks/boundary-conditions/extrapolation-400
    KIND exact  LABELS physics 1D
    DESC "extrapolation faces 400: flux of the cell state, report of the gas entering")
mose_test(BoundaryConditions.SupersonicInlet
    DIR checks/boundary-conditions/supersonic-inlet
    KIND exact  LABELS physics 1D
    DESC "supersonic inlet 405 with the flow direction normal")
mose_test(Connections.ChimeraFace
    DIR checks/connections/chimera-face
    KIND exact  LABELS numerics 2D 3D
    DESC "viscous stencil of chimera faces on a linear field")
mose_test(Connections.ConnectionFlux
    DIR checks/connections/connection-flux
    KIND identity  LABELS numerics 1D
    DESC "faces between blocks: one block equals two, bit for bit, and conserves")
mose_test(Connections.ConnectionOrder
    DIR checks/connections/connection-order
    KIND identity  LABELS implementation 2D
    DESC "sorted and reversed connection records give the same bits")
mose_test(Connections.SectorRadial
    DIR checks/connections/sector-radial
    KIND invariant  LABELS physics 3D
    DESC "3-D sector with rotational periodicity keeps a radial flow its own mirror image")
mose_test(Connections.SectorViscous
    DIR checks/connections/sector-viscous
    KIND exact  LABELS physics 3D
    DESC "viscous flux across the rotational periodicity of a 3-D sector")
mose_test(InputChecks.ExitCodes
    DIR checks/input-checks/exit-codes
    KIND contract  LABELS implementation 1D 2D
    DESC "error stops end with a non-zero exit code and their message")
mose_test(InputChecks.LimiterDefaults
    DIR checks/input-checks/limiter-defaults
    KIND contract  LABELS numerics 1D
    DESC "the declared limiter is the one that runs, set-up warnings say so")
mose_test(InputChecks.IRSTimeAccurate
    DIR checks/input-checks/irs-time-accurate
    KIND contract  LABELS numerics 1D
    DESC "residual smoothing is refused in time-accurate runs")
mose_test(InputChecks.PrecScope
    DIR checks/input-checks/prec-scope
    KIND contract  LABELS numerics 1D
    DESC "prec refused in time-accurate runs, HLLC-PC refused without prec")
mose_test(InputChecks.DecompositionMap
    DIR checks/input-checks/decomposition-map
    KIND contract  LABELS implementation 2D
    DESC "block-to-rank map: valid maps keep the bits, bad maps are refused")
mose_test(InputChecks.RuntimeInput
    DIR checks/input-checks/runtime-input
    KIND contract  LABELS implementation 1D 2D
    DESC "input.ini read again during a run: errors, fixed keys, warn and stop")
mose_test(MeshGeometry.WallDistance
    DIR checks/mesh-geometry/wall-distance
    KIND exact  LABELS physics 2D
    DESC "wall distance of the turbulence models on constructed meshes")
mose_test(MeshGeometry.MeshType
    DIR checks/mesh-geometry/mesh-type
    KIND exact  LABELS implementation 2D
    DESC "type of the mesh told by its nodes")
mose_test(OutputFiles.WallOutput
    DIR checks/output-files/wall-output
    KIND contract  LABELS implementation 2D
    DESC "wall files: extent, points, cell values, one file per wall face")
mose_test(OutputFiles.VtkNames
    DIR checks/output-files/vtk-names
    KIND contract  LABELS implementation 2D
    DESC "VTK outputs are well-formed XML with the right names and times")
mose_test(OutputFiles.VtkRawFormat
    DIR checks/output-files/vtk-raw-format
    KIND contract  LABELS implementation 1D
    DESC "vtk raw solutions hold raw appended data")
mose_test(OutputFiles.WallOutputResidual
    DIR checks/output-files/wall-output-residual
    KIND invariant  LABELS implementation 2D
    DESC "wall output does not change the residual")
mose_test(OutputFiles.ProbePosition
    DIR checks/output-files/probe-position
    KIND exact  LABELS implementation 2D
    DESC "probes given by position take the nearest cell")
mose_test(ViscousPatch.ViscousLinearPatch
    DIR checks/viscous-patch/viscous-linear-patch
    KIND exact  LABELS numerics 2D 3D
    DESC "patch test of the viscous fluxes: linear fields on any mesh")
mose_test(ViscousPatch.ViscousWallPatch
    DIR checks/viscous-patch/viscous-wall-patch
    KIND exact  LABELS numerics 2D 3D
    DESC "patch test next to walls of prescribed temperature")
mose_test(ViscousPatch.ViscousOpenPatch
    DIR checks/viscous-patch/viscous-open-patch
    KIND exact  LABELS numerics 2D 3D
    DESC "patch test next to inlets and outlets")

# ---------------------------------------------------------------------------
#  test/unit: Fortran programs that check one routine
# ---------------------------------------------------------------------------
add_fortran_executable(qcr2000_check unit/qcr2000_check.f90)
add_fortran_executable(cell_volume_check unit/cell_volume_check.f90)
add_fortran_executable(reflective_pressure_check unit/reflective_pressure_check.f90)
add_fortran_executable(volatile_order_check unit/volatile_order_check.f90)
add_fortran_executable(bilinear_distance_check unit/bilinear_distance_check.f90)
add_fortran_executable(remote_blocks_check unit/remote_blocks_check.f90)
add_fortran_executable(muscl_face_bounds_check unit/muscl_face_bounds_check.f90)
add_fortran_executable(prec_transform_check unit/prec_transform_check.f90)

mose_test(QCR2000Stress
    COMMAND qcr2000_check
    KIND unit  LABELS physics unit
    DESC "QCR2000 constitutive relation")
mose_test(CellVolume
    COMMAND cell_volume_check
    KIND unit  LABELS numerics unit
    DESC "cell volume from the metric tensor")
mose_test(ReflectivePressure
    COMMAND reflective_pressure_check
    KIND unit  LABELS numerics unit
    DESC "face pressure of a reflective wall")
mose_test(VolatileOrder
    COMMAND volatile_order_check 1 18014398509481984 3
    KIND unit  LABELS numerics unit
    DESC "rounding order of VOLATILE variables")
mose_test(BilinearDistance
    COMMAND bilinear_distance_check ${CMAKE_CURRENT_SOURCE_DIR}/unit/bilinear_distance_cases.txt
    KIND unit  LABELS numerics unit
    DESC "wall distance to a bilinear face")
mose_test(RemoteBlocks
    COMMAND remote_blocks_check
    KIND unit  LABELS implementation unit
    DESC "blocks held by an MPI rank")
mose_test(MusclFaceBounds
    COMMAND muscl_face_bounds_check 1 18014398509481984 3
    KIND unit  LABELS numerics unit
    DESC "MUSCL face values stay between the neighbours")
mose_test(PrecTransform
    COMMAND prec_transform_check
    KIND unit  LABELS numerics unit
    DESC "Weiss-Smith update of prim and prec, and the reference velocity")

# ---------------------------------------------------------------------------
#  Multirate integrator ([MOSE-Multirate]): one test per check of test/multirate/mr-<check>; the full forms are heavy
# ---------------------------------------------------------------------------
mose_test(Multirate.Smoke
    DIR multirate/mr-smoke  RUN "bash ./run.sh"
    KIND identity  LABELS multirate implementation 1D
    DESC "one level equals the single rate, conservation, restart, refusal at set-up")
mose_test(Multirate.OneLevel
    DIR multirate/mr-onelevel  RUN "bash ./run.sh"
    KIND identity  LABELS multirate implementation numerics 1D 2D 3D
    DESC "one level equals the single-rate step, bit for bit, on 33 configurations")
mose_test(Multirate.Threads
    DIR multirate/mr-threads  RUN "bash ./run.sh"
    KIND identity  LABELS multirate implementation 1D 2D 3D
    DESC "1 against 3 and 4 threads, bit for bit")
mose_test(Multirate.Blocks
    DIR multirate/mr-blocks  RUN "bash ./run.sh"
    KIND identity  LABELS multirate implementation 1D 2D
    DESC "equal levels across block faces, one block equals three")
mose_test(Multirate.Halo
    DIR multirate/mr-halo  RUN "bash ./run.sh"
    KIND exact  LABELS multirate implementation 1D 2D
    DESC "halo of every level against an independent construction")
mose_test(Multirate.Sym
    DIR multirate/mr-sym  RUN "bash ./run.sh"
    KIND invariant  LABELS multirate numerics 2D
    DESC "mirror symmetry of a centred pulse with levels")
mose_test(Multirate.Free
    DIR multirate/mr-free  RUN "bash ./run.sh"
    KIND invariant  LABELS multirate numerics 2D
    DESC "free stream preserved on a curvilinear periodic mesh with levels")
mose_test(Multirate.Conservation
    DIR multirate/mr-conservation  RUN "bash ./run.sh"
    KIND invariant  LABELS multirate numerics 1D 2D 3D
    DESC "mass, momentum and energy conserved with levels, diffusive, turbulence, mixture")
mose_test(Multirate.Restart
    DIR multirate/mr-restart  RUN "bash ./run.sh"
    KIND identity  LABELS multirate implementation 1D 2D
    DESC "exact restart and output instants with levels")
mose_test(Multirate.Order
    DIR multirate/mr-order  RUN "bash ./run.sh"
    KIND order  LABELS multirate numerics 1D 2D
    LOCK mr-order
    DESC "temporal order with levels: waves, vortex, plate, hot square, wall, SA, SST")
mose_test(Multirate.Deep
    DIR multirate/mr-deep  RUN "bash ./run.sh"
    KIND order  LABELS multirate numerics 1D
    LOCK mr-deep
    DESC "9 and 13 levels: order, conservation, 2000-step stability")
mose_test(Multirate.Bc
    DIR multirate/mr-bc  RUN "bash ./run.sh"
    KIND order  LABELS multirate numerics 1D
    DESC "time-dependent inlets 402 and 410 on the finest level: order at the inlet")
mose_test(Multirate.MovingLevels
    DIR multirate/mr-moving-levels  RUN "bash ./run.sh"
    KIND order  LABELS multirate numerics 1D
    DESC "order with level boundaries moving with a simple wave")
mose_test(Multirate.Sod
    DIR multirate/mr-sod  RUN "bash ./run.sh"
    KIND exact  LABELS multirate numerics 1D
    DESC "Sod tube with levels rebuilt every macro step against the exact solution")
mose_test(Multirate.Shock
    DIR multirate/mr-shock  RUN "bash ./run.sh"
    KIND exact  LABELS multirate numerics 2D
    DESC "planar shock through the faces between levels")
mose_test(Multirate.Rpost
    DIR multirate/mr-rpost  RUN "bash ./run.sh"
    KIND invariant  LABELS multirate numerics 1D 2D
    DESC "step taken against the step admitted on the new state")
mose_test(Multirate.Cfl4
    DIR multirate/mr-cfl4  RUN "bash ./run.sh"
    KIND invariant  LABELS multirate numerics 1D 2D
    DESC "non-local CFL condition: certificate, interface margin, redo")
mose_test(Multirate.Positivity
    DIR multirate/mr-positivity  RUN "bash ./run.sh"
    KIND robustness  LABELS multirate physics 1D
    DESC "no negative nu~, k or omega after a stage with levels")
mose_test(Multirate.Reject
    DIR multirate/mr-reject  RUN "bash ./run.sh"
    KIND contract  LABELS multirate implementation 1D 2D 3D
    DESC "34 configurations outside the scope are refused with their message")
mose_test(Multirate.FixedStep
    DIR multirate/mr-fixedstep  RUN "bash ./run.sh"
    KIND contract  LABELS multirate implementation 1D 2D
    DESC "fixed macro step and re-reading of the multirate keys")
mose_test(Multirate.Trace
    DIR multirate/mr-trace  RUN "bash ./run.sh"
    KIND exact  LABELS multirate implementation 1D
    DESC "stage sequence of one macro step against the reference recursion")
mose_test(Multirate.Daru
    DIR multirate/mr-daru  RUN "bash ./run.sh"
    KIND reference  LABELS multirate numerics 2D
    LOCK mr-daru
    DESC "Daru-Tenaud tube on 50x25 cells: multirate against the single rate")
mose_test(Multirate.Nozzle
    DIR multirate/mr-nozzle  RUN "bash ./run.sh"
    KIND reference  LABELS multirate numerics 1D
    DESC "nozzle start-up: moving shock across levels against the single rate")
mose_test(Multirate.OrderFull
    DIR multirate/mr-order  RUN "bash ./run.sh"
    KIND order  LABELS multirate numerics 1D 2D
    LOCK mr-order  TIMEOUT 14400  ENV MR_FORM=full
    DESC "full form of Multirate.Order (about 1 h)")
mose_test(Multirate.DeepFull
    DIR multirate/mr-deep  RUN "bash ./run.sh"
    KIND order  LABELS multirate numerics 1D
    LOCK mr-deep  TIMEOUT 14400  ENV MR_FORM=full
    DESC "full form of Multirate.Deep, with the 13-level order (about 1 h)")
mose_test(Multirate.DaruMedium
    DIR multirate/mr-daru  RUN "bash ./run.sh medium"
    KIND reference  LABELS multirate numerics 2D
    LOCK mr-daru  TIMEOUT 3600
    DESC "Daru-Tenaud tube on 100x50 cells (about 20 min)")
mose_test(Multirate.DaruFull
    DIR multirate/mr-daru  RUN "bash ./run.sh full"
    KIND reference  LABELS multirate numerics 2D
    LOCK mr-daru  TIMEOUT 720000
    DESC "Daru-Tenaud tube on 500x250 cells, as the reference (100 h)")

# ---------------------------------------------------------------------------
#  needs-mpi: built with USE_MPI=ON (the second CI job); up to 4 ranks x threads each, so PROCESSORS says so
# ---------------------------------------------------------------------------
if(USE_MPI AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # the probe library of MPIThreadFunnel and of Multirate.Hybrid: it reports an MPI call from a thread other than the main one
    find_package(Threads REQUIRED)
    add_library(mpi_thread_probe MODULE mpi/thread-funnel/mpi_thread_probe.c)
    target_link_libraries(mpi_thread_probe PRIVATE ${CMAKE_DL_LIBS} Threads::Threads)
    set_target_properties(mpi_thread_probe PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endif()

mose_test(MPIRestart
    DIR mpi/restart
    KIND identity  LABELS implementation 2D
    REQUIRES MPI
    DESC "exact restart on 2 ranks x 2 threads")
mose_test(MPIBoundaryRecords
    DIR mpi/boundary-records
    KIND identity  LABELS implementation 2D
    REQUIRES MPI
    DESC "boundary fluxes identical on any ranks x threads and any record order")
mose_test(MPIVertexStencil
    DIR mpi/vertex-stencil
    KIND identity  LABELS implementation 2D
    REQUIRES MPI
    DESC "viscous stencil at a vertex of four blocks: same bits on any layout")
mose_test(MPIBitIdentity
    DIR mpi/bit-identity
    KIND identity  LABELS implementation 2D 3D
    REQUIRES MPI
    DESC "same output bits on 1x1, 1x4, 2x1, 4x1, 2x2 ranks x threads")
mose_test(MPIThreadFunnel
    DIR mpi/thread-funnel
    KIND contract  LABELS implementation 2D
    REQUIRES MPI LINUX
    DESC "no MPI call from a thread other than the main one")
mose_test(MPIConnectionFlux
    DIR mpi/connection-flux
    KIND identity  LABELS implementation 1D
    REQUIRES MPI
    DESC "faces between ranks give the one-rank bits and conserve")
mose_test(MPIConnectionOrder
    DIR mpi/connection-order
    KIND identity  LABELS implementation 2D
    REQUIRES MPI
    DESC "faces between ranks, sorted and reversed records: same bits")
mose_test(MPIDecompositionMap
    DIR mpi/decomposition-map
    KIND contract  LABELS implementation 2D
    REQUIRES MPI
    DESC "block-to-rank map: valid maps are used, bad maps are refused")
mose_test(MPIExitCodes
    DIR mpi/exit-codes
    KIND contract  LABELS implementation 1D
    REQUIRES MPI
    DESC "an input error found by rank 0 stops every rank with exit code 1")
mose_test(MPIRuntimeInput
    DIR mpi/runtime-input
    KIND contract  LABELS implementation 2D
    TIMEOUT 120  REQUIRES MPI
    DESC "input.ini read again during a run on several ranks")
mose_test(Multirate.Hybrid
    DIR multirate/mr-hybrid  RUN "bash ./run.sh"
    KIND identity  LABELS multirate implementation 1D 2D 3D
    PROCESSORS 8  REQUIRES MPI
    DESC "multirate on MPI ranks and threads: same bits as 1x1, timers")
mose_test(Multirate.RestartHybrid
    DIR multirate/mr-restart-hybrid  RUN "bash ./run.sh"
    KIND identity  LABELS multirate implementation 2D
    PROCESSORS 4  REQUIRES MPI
    DESC "multirate exact restart on 2 ranks x 2 threads")

if(TARGET mpi_thread_probe)
    foreach(_t MPIThreadFunnel Multirate.Hybrid)
        if(TEST ${_t})
            set_property(TEST ${_t} APPEND PROPERTY ENVIRONMENT "MOSE_MPI_THREAD_PROBE=$<TARGET_FILE:mpi_thread_probe>")
        endif()
    endforeach()
endif()
if(TEST MPIThreadFunnel)
    set_tests_properties(MPIThreadFunnel PROPERTIES SKIP_RETURN_CODE 77)   # no probe library: the test skips itself
endif()

# ---------------------------------------------------------------------------
#  The tiers file and the registrations must list the same tests; the inventory goes to the build folder
# ---------------------------------------------------------------------------
get_property(_names GLOBAL PROPERTY MOSE_TEST_NAMES)
foreach(_n IN LISTS _tier_names)
    if(NOT _n IN_LIST _names)
        message(FATAL_ERROR "test/tiers.txt lists '${_n}', which no mose_test() registers")
    endif()
endforeach()
string(ASCII 9 _tab)
set(_inventory "")
foreach(_n IN LISTS _names)
    get_property(_kind GLOBAL PROPERTY MOSE_TEST_${_n}_KIND)
    get_property(_desc GLOBAL PROPERTY MOSE_TEST_${_n}_DESC)
    get_property(_labels GLOBAL PROPERTY MOSE_TEST_${_n}_LABELS)
    get_property(_built GLOBAL PROPERTY MOSE_TEST_${_n}_BUILT)
    string(REPLACE ";" "," _labels "${_labels}")
    string(APPEND _inventory "${_n}${_tab}${MOSE_TIER_${_n}}${_tab}${MOSE_SECONDS_${_n}}${_tab}${_kind}${_tab}${_labels}${_tab}${_built}${_tab}${_desc}\n")
endforeach()
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/inventory.tsv" "${_inventory}")
