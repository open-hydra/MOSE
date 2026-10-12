# ---------------------------------------------------------------------------
#  CTest definitions for MOSE validation cases
#
#  Tiers are expressed through LABELS (see test/TESTING-PLAN.md):
#    fast            -> tiny, dependency-free; run by the pre-push hook
#                       (ctest -L fast)
#    validation      -> analytic / reference comparison
#    needs-cantera   -> requires USE_CANTERA=ON  (finite-rate / equilibrium)
#    needs-sundials  -> requires USE_SUNDIALS=ON (stiff ODE chemistry)
#    needs-mpi       -> requires USE_MPI=ON
#
#  PR CI runs everything except the heavy-dependency labels:
#    ctest -LE "needs-cantera|needs-sundials|needs-mpi"
#  and, in a second job built with USE_MPI=ON (gfortran, Open MPI), the MPI tier:
#    ctest -L needs-mpi
# ---------------------------------------------------------------------------

# Common environment for all solver tests
set(_MOSE_TEST_ENV
    "OMP_NUM_THREADS=1"
    "KMP_STACKSIZE=100M"
)

# Helper: register a case that runs "<run-cmd> && python3 -B <verifier>"
# in a given working directory, with the given labels.
function(mose_add_test)
    cmake_parse_arguments(T "" "NAME;DIR;RUN;VERIFY;TIMEOUT" "LABELS" ${ARGN})
    if(NOT T_TIMEOUT)
        set(T_TIMEOUT 300)
    endif()
    add_test(
        NAME    ${T_NAME}
        COMMAND bash -c "${T_RUN} && python3 -B ${T_VERIFY}"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/${T_DIR}"
    )
    set_tests_properties(${T_NAME} PROPERTIES
        LABELS      "${T_LABELS}"
        TIMEOUT     ${T_TIMEOUT}
        ENVIRONMENT "${_MOSE_TEST_ENV}"
    )
endfunction()

# ===========================================================================
#  FAST tier  (pre-push: ctest -L fast)  — tiny, no heavy dependencies
# ===========================================================================

# 1-D Riemann problems (4-mesh grid-convergence via the "test" verb)
mose_add_test(NAME Sod79  DIR 1D/Sod79
              RUN "./MOSE.sh solve" VERIFY check.py
              LABELS "fast;validation;1D")

mose_add_test(NAME Einfeldt91 DIR 1D/Ein91
              RUN "./MOSE.sh test" VERIFY verify.py
              LABELS "fast;validation;1D")

mose_add_test(NAME Noh87 DIR 1D/Noh87
              RUN "./MOSE.sh test" VERIFY verify.py
              LABELS "fast;validation;1D")

mose_add_test(NAME Toro99 DIR 1D/Tor99
              RUN "./MOSE.sh test" VERIFY verify.py
              LABELS "fast;validation;1D")

# 1-D Shu-Osher (Mach-3 shock / sine interaction; convergence vs N=1600 ref)
#   fast tier : single coarse grid (N=200) vs the reference
mose_add_test(NAME ShuOsher DIR 1D/Shu-Osher
              RUN "./MOSE.sh solve" VERIFY verify.py
              LABELS "fast;validation;numerics;1D")
#   validation tier : 3-grid study (N=200/400/800), monotone convergence + figure
mose_add_test(NAME ShuOsherConvergence DIR 1D/Shu-Osher
              RUN "./MOSE.sh test" VERIFY convergence.py
              LABELS "validation;numerics;1D")

# 2-D steady Euler (analytic Prandtl-Meyer)
mose_add_test(NAME PrandtlMeyer DIR 2D/euler/prandtl-meyer
              RUN "./MOSE.sh solve" VERIFY check.py
              LABELS "fast;validation;2D")

# 2-D isentropic-vortex convection (N1): order-of-accuracy on 32/64/128 grids
mose_add_test(NAME IsentropicVortex DIR 2D/euler/isentropic-vortex
              RUN "./MOSE.sh test" VERIFY verify.py
              LABELS "validation;numerics;2D")

# 2-D viscous Method of Manufactured Solutions (N2): order-of-accuracy of the
# *viscous* operator. Uses the MOSE_mms driver (analytic manufactured source).
# CTest is single-threaded, so it runs the fast 32/64 subset ("testfast"); the
# full 32/64/128 study ("test") is used offline to make the V&V figure.
mose_add_test(NAME MMSViscous DIR 2D/viscous/mms
              RUN "./MOSE.sh testfast" VERIFY verify.py
              LABELS "validation;numerics;2D"
              TIMEOUT 900)

# ===========================================================================
#  Additional fast tests defined in their own directories
#  (added incrementally — see test/TESTING-PLAN.md, section 4)
# ===========================================================================
# Helper: register a self-contained fast test driven by its own run.sh.
function(mose_add_fast_script NAME DIR)
    if(EXISTS "${CMAKE_SOURCE_DIR}/test/${DIR}/run.sh")
        add_test(
            NAME    ${NAME}
            COMMAND bash -c "./run.sh"
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/${DIR}"
        )
        set_tests_properties(${NAME} PROPERTIES
            LABELS      "${ARGN}"
            TIMEOUT     120
            ENVIRONMENT "${_MOSE_TEST_ENV}"
        )
    endif()
endfunction()

# Riemann smoke-matrix (F5): every shock-capturing solver stays positive;
# shock-detector = Tramel turns the shock sensor on.
mose_add_fast_script(RiemannSmoke   fast/riemann-smoke  "fast;numerics;1D")

# Numerics smoke (N4): limiter / time-scheme / reconstruction / int-var sweep.
mose_add_fast_script(NumericsSmoke  fast/numerics-smoke "fast;numerics;1D")

# Time integration (TimeIntegration): the parts below, each in its folder of
# test/fast/time-integration, run one after the other by
# test/fast/time-integration/run.sh.
# Time control: second order in time with a time-dependent inlet (BC 402,
# boundary data at the stage times), time tables before their first instant,
# res-threshold not ending a time-accurate run.
# Local time step: stable on square and skewed 2-D cells, on cubic 3-D cells
# and in diffusion-dominated cells (heat conduction, species), with Heun and RK3;
# heat conduction on 3-D rhombohedra at 115 degrees with Heun at vnn 0.5.
# Turbulence sources in time-accurate runs: nu~, k and omega converge with the
# order of Heun and RK3, and no stage makes them negative near a wall; SST-RC
# prints the start-up warnings of k-omega models and of -RC; the
# vortex-stretching function of Wilcox 2006 on a 3-D flow; the asymptotic omega
# wall value of Wilcox 2006 and SST.
# Quasi-2D case with symmetry k faces: round-off in w must not grow.
mose_add_fast_script(TimeIntegration fast/time-integration "fast;numerics;1D;2D;3D")

# Boundary conditions (BoundaryConditions): the parts below, each in its
# folder of test/fast/boundary-conditions, run one after the other by
# test/fast/boundary-conditions/run.sh.
# Inlet boundary conditions: the far-field inlet (407) with the flow direction
# normal gives the field of the equivalent angles, to the bit, with no
# transverse velocity; inlets at total conditions (401, 407): when the gas
# enters (from rest too), subsonic and choked faces, the slip wall of an inlet the
# reservoir cannot feed, and the outlet of an ambient face it cannot feed with the gas
# moving in; one step, mass fluxes against the face state.
# Boundary input: the property lines of bc.txt are read as written or refused
# (the inlet 402 table under any file name; a missing or decreasing table stops
# the run).
# Rotating frame: stationary-face-1, -2, ... of [MOSE-Rotating-Frame] are read
# and listed, and the run ends (the count of the keys used to loop forever).
mose_add_fast_script(BoundaryConditions fast/boundary-conditions "fast;implementation;numerics;1D;2D")

# Closed-box mass conservation (I4): slip-wall box conserves mass to round-off.
mose_add_fast_script(Conservation   fast/conservation   "fast;implementation;1D")

# Serial vs OpenMP equivalence (I5): 1 vs 4 threads bit-identical, the Sod tube
# and a chimera case whose chimera faces meet at the edges of a block.
mose_add_fast_script(OpenMPEquiv    fast/openmp-equiv   "fast;implementation;1D")

# Faces between blocks (Connections): the parts below, each in its folder of
# test/fast/connections, run one after the other by
# test/fast/connections/run.sh.
# Viscous stencil of chimera faces (F9): one evaluation of the residual of a
# linear field next to chimera faces, also at the edges of the block face, in
# 2-D, in 3-D and on a 2-D axisymmetric wedge (against the same mesh as one
# block).
# Faces between blocks (101, 201): one flux per face. A periodic shock tube as
# one block and as two blocks gives the same bits, and conserves mass,
# momentum and energy, with and without shock sensor and preconditioning.
# Faces between blocks of the same parity (face 4 of one block on face 2 of
# the other, index directions permuted): the record that computes the flux
# of the face does not follow the order of bc.txt; sorted and reversed
# records give the same bits.
# Sector radial: a 3-D sector closed by the rotational periodicity (type 200 on
# the k faces) keeps a radial flow its own mirror image.
mose_add_fast_script(Connections fast/connections "fast;implementation;numerics;1D;2D;3D")

# Restart round-trip (I3/F6): split+restart reproduces the uninterrupted run;
# at a restart a missing probe file is written again, an existing one appended,
# with the probe rows of the uninterrupted run.
mose_add_fast_script(Restart        fast/restart        "fast;implementation;1D")

# Output files (OutputFiles): the parts below, each in its folder of
# test/fast/output-files, run one after the other by
# test/fast/output-files/run.sh.
# Wall output: a vtk raw run with no-slip walls ends, and every wall file has
# the extent and the points of its face.
mose_add_fast_script(OutputFiles fast/output-files "fast;implementation;2D")

# Input checks (InputChecks): the parts below, each in its folder of
# test/fast/input-checks, run one after the other by
# test/fast/input-checks/run.sh.
# Exit codes of the error stops: a run that MOSE stops on an error ends with a
# non-zero exit code and its message: the Sod tube as a viscous mixture of four
# species without a diffusion table, a run that diverges, a NaN in the initial
# condition, an initial condition with too few variables, a rotating frame
# with a zero axis, a probe of an unknown variable; a valid run with chemistry and
# exclude-blocks ends with 0. A steady run with a soot model ends normally and
# its diagnostic output names every variable it writes; a steady single-species
# run without thermo or transport output stops with a message that says so.
# Limiter defaults: the reconstruction that runs is the declared one (MUSCL
# without flux-limiter = van Leer; first-order ignores a given limiter), and
# the set-up warnings say so.
# Implicit residual smoothing is a steady-state accelerator: a time-accurate
# run that asks for it (irs = true or irs-beta > 0) is refused at set-up.
# Block-to-rank assignment from INPUT/decomposition.map on one rank: a map
# with every block on rank 0, also one written by the mesh decomposer (MDB)
# itself, gives the bits of the run without it; maps that do not fit the mesh
# or the run stop it at setup with the fault named.
# input.ini read again during a run (ini-diter): a value out of the rules of
# the registry stops the run with its message, as at setup; read again
# unchanged, with keys left out or written, it gives the run without readings.
mose_add_fast_script(InputChecks fast/input-checks "fast;implementation;numerics;1D;2D")

# Manufactured fields for the RANS right-hand sides (2D/viscous/mms/residuals).
# Spalart-Allmaras with variable density: the right-hand side of the equation
# for rho nu~ converges to the standard conservation form; with SA-RC, the
# right-hand side of a 3-D flow does not depend on how the flow is turned, it
# takes the rotation of the frame in r~, a uniform flow keeps nu~, and air at
# rest next to a wall, or flowing along it on stretched cells, changes nu~ as
# with SA.
# Part 4: the viscous force on the axisymmetric wedge (Homann's flow and a field
# with dilatation and swirl) converges to the Navier-Stokes one with second order.
# Part 5: next to a wall the Spalart-Allmaras sources see the gradient of nu~
# that the wall value nu~ = 0 gives (linear profile at rest).
mose_add_fast_script(MMSResiduals 2D/viscous/mms/residuals "fast;physics;2D;3D")

# Steady state with the point-implicit factor: it must not depend on the time
# step (a short hypersonic plate, Spalart-Allmaras and SST). Validation tier:
# it converges 14 000 iterations per model.
mose_add_fast_script(PointImplicitSteady 2D/viscous/point-implicit-steady "validation;numerics;2D")

# Local time step next to the axis of the wedge (BC 200): the stress on the wedge
# faces stiffens the rows next to the axis. At the largest vnn the step rule
# declares stable (RK3 0.62, Heun 0.5), on square and on elongated cells, and at
# the recommended vnn (RK3 0.5, Heun 0.4) on cells four times as long in x, with
# Prandtl 0.72 and 0.9, no NaN.
# Validation tier: eight runs of 2000 iterations, 155 s of CPU at one thread on
# a loaded machine (163 s of wall time), hence a TIMEOUT of its own.
mose_add_fast_script(WedgeAxisStep 2D/viscous/wedge-axis-step "validation;numerics;2D")
if(TEST WedgeAxisStep)
    set_tests_properties(WedgeAxisStep PROPERTIES TIMEOUT 600)
endif()

# Ablative wall BCs (503 melting / 504 pyrolysis): an isothermal blown column,
# which is an exact steady solution, checked against its closed form. These
# boundary conditions have no other coverage in the suite.
# The melting column also runs without schmidt: the default, 0, needs no
# diffusion table with a single species.
mose_add_fast_script(AblatingWall 2D/viscous/ablating-wall "fast;validation;physics;2D")

# Mesh geometry (MeshGeometry): the parts below, each in its folder of
# test/fast/mesh-geometry, run one after the other by
# test/fast/mesh-geometry/run.sh.
# Wall distance: the distance to the wall of the turbulence models is the
# minimum distance to the wall surface (convex corner, curved wall, nozzle lip
# on an axisymmetric wedge), and to its periodic images, also when the periods
# are much shorter than the mesh or nearly dependent, read back from the
# Spalart-Allmaras source.
# Mesh type: a plane 2-D mesh, an axisymmetric wedge and a 3-D slab periodic in
# k are told apart by every node of every block; blocks that disagree stop.
mose_add_fast_script(MeshGeometry fast/mesh-geometry "fast;implementation;physics;2D")

# Quadratic constitutive relation QCR2000: a unit check of the stress as the
# solver calls it (laminar viscosity, then eddy viscosity), against the
# definition. QCR2000 is not selectable from the input, so no run passes there.
add_fortran_executable(qcr2000_check unit/qcr2000_check.f90)
add_test(NAME QCR2000Stress COMMAND qcr2000_check)
set_tests_properties(QCR2000Stress PROPERTIES LABELS "fast;physics;unit" TIMEOUT 60)

# ===========================================================================
#  MPI tier (needs-mpi): registered only in builds with USE_MPI=ON
#  Each test generates its own small multi-block case and launches the solver
#  itself through the MPI launcher CMake found, with at most 4 ranks; it
#  asserts an invariant of the parallel run.
# ===========================================================================
if(USE_MPI)
    function(mose_add_mpi_script NAME DIR)
        if(EXISTS "${CMAKE_SOURCE_DIR}/test/${DIR}/run.sh")
            add_test(
                NAME    ${NAME}
                COMMAND bash -c "./run.sh"
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/${DIR}"
            )
            set_tests_properties(${NAME} PROPERTIES
                LABELS      "${ARGN}"
                TIMEOUT     600
                PROCESSORS  4
                ENVIRONMENT "${_MOSE_TEST_ENV};MOSE_MPIEXEC=${MPIEXEC_EXECUTABLE};MOSE_MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
            )
        endif()
    endfunction()

    # Boundary records: the boundary-flux loop gives the same bits on any
    # ranks x threads layout and for any order of the records in bc.txt.
    mose_add_mpi_script(MPIBoundaryRecords mpi/boundary-records "needs-mpi;implementation;2D")

    # Viscous connection stencil at a vertex shared by four blocks: the same
    # bits on 1x1, 1x4, 2x1, 4x1 and 2x2 ranks x threads; with equilibrium
    # chemistry, the same bits on 1x1 and 2x1.
    mose_add_mpi_script(MPIVertexStencil mpi/vertex-stencil "needs-mpi;implementation;2D")

    # MPI calls only from the main thread (MPI_THREAD_FUNNELED): a probe
    # library preloaded into the solver counts the MPI calls made by any
    # other thread; it must count none on 2x2, 2x4 and 2x2 with IRS (steady).
    find_package(Threads REQUIRED)
    add_library(mpi_thread_probe MODULE mpi/thread-funnel/mpi_thread_probe.c)
    target_link_libraries(mpi_thread_probe PRIVATE ${CMAKE_DL_LIBS} Threads::Threads)
    set_target_properties(mpi_thread_probe PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    mose_add_mpi_script(MPIThreadFunnel mpi/thread-funnel "needs-mpi;implementation;2D")
    if(TEST MPIThreadFunnel)
        set_property(TEST MPIThreadFunnel APPEND PROPERTY
            ENVIRONMENT "MOSE_MPI_THREAD_PROBE=$<TARGET_FILE:mpi_thread_probe>")
    endif()

    # Faces between ranks: the two-block cases of ConnectionFlux with one block
    # per rank give the one-rank and the one-block bits, and conserve.
    mose_add_mpi_script(MPIConnectionFlux mpi/connection-flux "needs-mpi;implementation;1D")

    # The 2-D case of ConnectionOrder with one block per rank: sorted and
    # reversed records give the one-rank bits.
    mose_add_mpi_script(MPIConnectionOrder mpi/connection-order "needs-mpi;implementation;2D")

    # Block-to-rank assignment from INPUT/decomposition.map: on 2 and 4 ranks
    # every map (MOSE's own assignment as MDB writes it, permutations, an idle
    # rank 0, 110 blocks, fewer blocks than ranks) gives the one-rank field; a
    # map for other ranks stops every rank with the same message.
    mose_add_mpi_script(MPIDecompositionMap mpi/decomposition-map "needs-mpi;implementation;2D")

    # Exit codes under MPI: an input error that rank 0 finds at set-up stops both
    # ranks of a 2-rank run within 60 s, with exit code 1 to 128.
    mose_add_mpi_script(MPIExitCodes mpi/exit-codes "needs-mpi;implementation;1D")
endif()
