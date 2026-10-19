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
# The restart tests also write and read the Tecplot binary format (.szplt), which ORION writes with TecIO
if(USE_TECIO)
    list(APPEND _MOSE_TEST_ENV "MOSE_TECIO=1")
endif()

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
#   the two run in the same folder and write OUTPUT/field.tec there: never at the same time (ctest -j)
set_tests_properties(ShuOsher ShuOsherConvergence PROPERTIES RESOURCE_LOCK Shu-Osher)

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

# Manufactured solutions on general grids (validation tier): the fields of the
# MOSE_mms driver (MOSE_MMS_FIELD) on two-block grid families whose blocks meet
# with permuted indices, zigzag (every cell skewed, its neighbours the other
# way) and iface (a cell-size jump of 1.5 and kinked grid lines): 2-D periodic
# and with walls on N = 10, 20, 3-D zigzag on N = 8, 12 ("testfast"; the full
# studies, "test", make the V&V figures).
mose_add_test(NAME MMSViscousSkewed DIR 2D/viscous/mms/skewed
              RUN "./MOSE.sh testfast" VERIFY verify.py
              LABELS "validation;numerics;2D"
              TIMEOUT 900)
mose_add_test(NAME MMSViscous3D DIR 3D/viscous/mms
              RUN "./MOSE.sh testfast" VERIFY verify.py
              LABELS "validation;numerics;3D"
              TIMEOUT 900)

# Planar nozzle on grid families (validation tier): a turbulent converging-
# diverging nozzle (Spalart-Allmaras, SST) on O and H grids whose four blocks
# have different index orientations stays mirror-symmetric about its centre
# line, and the zig-zag grid gives a finite flow ("testfast": 300 iterations on
# the coarse level; the full study, "test", runs every family to convergence).
mose_add_test(NAME PlanarNozzleGrids DIR 2D/viscous/planar-nozzle
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
# res-threshold not ending a time-accurate run, and ending a steady run when the
# density residual has fallen by that factor from the first iteration, at any
# level of density and across a restart.
# Local time step: stable on square and skewed 2-D cells, on cubic 3-D cells
# and in diffusion-dominated cells (heat conduction, species), with Heun and RK3;
# heat conduction on 3-D rhombohedra at 115 degrees with Heun at vnn 0.5;
# with integration-variables = prec, next to a wall with SST and Wilcox 2006,
# the reference velocity of the step takes the eddy viscosity of the model at
# the strain rate of the cell.
# Turbulence sources in time-accurate runs: nu~, k and omega converge with the
# order of Heun and RK3, and no stage makes them negative near a wall, where
# the step takes the bound of the stiffness of each cell; SST-RC prints the
# start-up warnings of k-omega models and of -RC; the
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
# Sign of the wall heat flux: q_wall of a wall 301 is positive into the fluid on
# every face (the same wall on both sides of a channel heats it alike), and the
# wall output writes qw as k grad T . n along the face normal (301 and 302).
mose_add_fast_script(BoundaryConditions fast/boundary-conditions "fast;implementation;numerics;physics;1D;2D")

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
# Sector viscous: the two k faces of a 3-D sector closed by type 200 are
# interior faces of the annulus, and carry their viscous flux: the sector of
# 90 degrees gives the viscous force of the same cells inside the sector of
# 180 degrees, and the normal stress of symmetry planes on a radial flow.
mose_add_fast_script(Connections fast/connections "fast;implementation;numerics;physics;1D;2D;3D")

# Restart round-trip (I3/F6): split+restart reproduces the uninterrupted run;
# at a restart a missing probe file is written again, an existing one appended,
# with the probe rows of the uninterrupted run.
# Restart counters (I3): with vtk raw, steady, time-based-output and multigrid
# restarts reproduce the uninterrupted runs bit for bit; a multigrid restart from
# Tecplot goes on with the global iteration count; a time-accurate restart from a
# steady solution starts at t = 0; a restart from an older solution keeps no
# residual-history or probe row twice; a first part that ends on an output
# iteration writes the solution and the probe row due there; a first part
# killed at its time limit leaves in the history files every row it wrote.
mose_add_fast_script(Restart        fast/restart        "fast;implementation;1D")

# Output files (OutputFiles): the parts below, each in its folder of
# test/fast/output-files, run one after the other by
# test/fast/output-files/run.sh.
# Wall output: a vtk raw run with no-slip walls ends, and every wall file has
# the extent, the points and the cell values of its face.
# VTK outputs (field, wall, diagnostic) are well-formed XML with the right names.
# vtk raw outputs (field and diagnostic, steady and time-accurate) hold raw
# appended data in every block file, as sol-format asks.
# Wall output and residual: the wall output must not change the residual that
# the diagnostic output prints next (a converged laminar hypersonic plate).
# Probes given by position record the cell whose centre is the nearest, never
# a ghost cell: probes at cell centres hold the values of their cells.
mose_add_fast_script(OutputFiles fast/output-files "fast;implementation;1D;2D")

# Input checks (InputChecks): the parts below, each in its folder of
# test/fast/input-checks, run one after the other by
# test/fast/input-checks/run.sh.
# Exit codes of the error stops: a run that MOSE stops on an error ends with a
# non-zero exit code and its message: the Sod tube as a viscous mixture of four
# species without a diffusion table, a run that diverges, a NaN or an infinity
# in the initial condition (a stop at set-up, exit code 1, in every build), in
# the time table of an inlet 402, in a property line of bc.txt and in the file
# of a face 410, a file of a face 410 that cannot be read (missing, cut, with a
# word among its values, an unreadable zone header or a value too many) or
# whose zones go back in time, an initial condition with too few variables, a
# rotating frame
# with a zero axis, or with an axis or an origin that is not finite,
# chemistry = equilibrium without the elements of the species,
# a value outside the rule of the registry (an interval, a list of reals, a list
# of integers), an input.ini that cannot be read (missing, empty, of comment
# lines only), a key of input.ini whose value does not read as its type, a
# probe of an unknown variable; a valid run with chemistry and
# exclude-blocks ends with 0, and a block excluded from chemistry = equilibrium
# keeps the bits of the run without chemistry, also with multigrid. A steady run
# with a soot model ends normally and its diagnostic output names every variable
# it writes; a steady single-species run without thermo or transport output
# stops with a message that says so. A restart from a solution whose time or a
# counter is not finite stops at set-up.
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
# with dilatation and swirl) converges to the Navier-Stokes one with second order;
# a viscous run with the axis at r = 0 on the wedge or on a 3-D sector stops at setup,
# and on the polar line of a 3-D mesh with diffusive-metric = inverse-mean or mean.
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

# Patch test of the viscous fluxes (ViscousPatch): the parts below, each in
# its folder of test/fast/viscous-patch, run one after the other by
# test/fast/viscous-patch/run.sh.
# Patch test of the viscous fluxes: a linear velocity and a linear temperature
# have zero viscous force and heat conduction in every cell, so the face
# gradients must give zero to round-off on meshes with nodes moved at random,
# folds between blocks, a fan of index lines at a wall (the second layer of
# extrapolated ghost cells turned over), points shared by 3 and 5 blocks, symmetry planes and
# walls with a prescribed heat flux; the old face metrics (inverse-mean, mean)
# against an independent replica of their formulas. ViscousWallPatch: the
# cells next to walls of prescribed temperature varying along the wall.
# Cost: 40 and 61 one-step solver runs, about 9 s and 13 s of CPU time at
# one thread.
mose_add_fast_script(ViscousPatch fast/viscous-patch "fast;numerics;2D;3D")

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
# k are told apart by every node of every block; blocks that disagree stop, and
# so does a cell folded at a corner (2-D, 3-D, and a cell of a real 3-D mesh).
mose_add_fast_script(MeshGeometry fast/mesh-geometry "fast;implementation;physics;2D")

# Quadratic constitutive relation QCR2000: a unit check of the stress as the
# solver calls it (laminar viscosity, then eddy viscosity), against the
# definition. QCR2000 is not selectable from the input, so no run passes there.
add_fortran_executable(qcr2000_check unit/qcr2000_check.f90)
add_test(NAME QCR2000Stress COMMAND qcr2000_check)
set_tests_properties(QCR2000Stress PROPERTIES LABELS "fast;physics;unit" TIMEOUT 60)

# Cell volume (CellVolume): a unit check of the volume of Compute_Metric_Tensor,
# the integral of det J of the trilinear cell, against its exact value on cells
# with twisted faces, a non-convex cell, a collapsed edge, 2-D, wedge and
# left-handed cells, and on a twisted mesh against the volume enclosed by its
# boundary faces; a cell whose volume has not the sign of det A at its centre
# raises the flag of degenerate or inverted cells.
add_fortran_executable(cell_volume_check unit/cell_volume_check.f90)
add_test(NAME CellVolume COMMAND cell_volume_check)
set_tests_properties(CellVolume PROPERTIES LABELS "fast;numerics;unit" TIMEOUT 60)

# Order of rounding (VolatileOrder): a difference stored in a VOLATILE variable
# and read back is rounded as written, also with -ffast-math, which lets the
# compiler evaluate (q - c1) - (c2 - c1) as q - c2; Compute_Metric_Tensor takes
# its node differences so. The values come from the command line (no folding).
add_fortran_executable(volatile_order_check unit/volatile_order_check.f90)
add_test(NAME VolatileOrder COMMAND volatile_order_check 1 18014398509481984 3)
set_tests_properties(VolatileOrder PROPERTIES LABELS "fast;numerics;unit" TIMEOUT 60)

# Wall distance to a bilinear face (BilinearDistance): the certified minimum of
# the wall distance against the exact distance of the stationarity equations on
# twisted, plane, non-convex, degenerate and far faces (within 288 u S, with its
# certificate), and the fast-path test of the four triangles: taken by plane
# faces, not by twisted ones (cases: unit/bilinear_distance_cases.txt).
add_fortran_executable(bilinear_distance_check unit/bilinear_distance_check.f90)
add_test(NAME BilinearDistance
         COMMAND bilinear_distance_check ${CMAKE_CURRENT_SOURCE_DIR}/unit/bilinear_distance_cases.txt)
set_tests_properties(BilinearDistance PROPERTIES LABELS "fast;numerics;unit" TIMEOUT 120)

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

    # Exact restart from vtk raw on 2 ranks x 2 threads: 24 + 16 steps, the
    # first part ending on an output step, give the bits of 40 steps in one
    # go, output numbering and counters included.
    mose_add_mpi_script(MPIRestart mpi/restart "needs-mpi;implementation;2D")

    # Boundary records: the boundary-flux loop gives the same bits on any
    # ranks x threads layout and for any order of the records in bc.txt.
    mose_add_mpi_script(MPIBoundaryRecords mpi/boundary-records "needs-mpi;implementation;2D")

    # Viscous connection stencil at a vertex shared by four blocks: the same
    # bits on 1x1, 1x4, 2x1, 4x1 and 2x2 ranks x threads; with equilibrium
    # chemistry, the same bits on 1x1 and 2x1.
    mose_add_mpi_script(MPIVertexStencil mpi/vertex-stencil "needs-mpi;implementation;2D")

    # The same output bits on 1x1, 1x4, 2x1, 4x1 and 2x2 ranks x threads: a
    # 2x2-block box with walls, inlets and outlets, Spalart-Allmaras and SST
    # steady with implicit residual smoothing, Spalart-Allmaras time-accurate,
    # and the points shared by 3 and by 5 blocks of ViscousLinearPatch.
    mose_add_mpi_script(MPIBitIdentity mpi/bit-identity "needs-mpi;implementation;2D;3D")

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
    # ranks of a 2-rank run within 60 s, with exit code 1 and the STOP line of both
    # ranks; an input.ini that rank 0 cannot read, whose outcome rank 0 passes to
    # every rank, stops them within 60 s, with exit code 1 to 128; a value out of
    # its rule written during the run in the input.ini of rank 0, which rank 0
    # reads again and passes to every rank (each rank in its own folder), stops
    # both ranks with exit code 1.
    # A value of the initial condition that is not finite stops both ranks with exit
    # code 1, also when rank 1 alone holds it; so does an initial condition that
    # rank 1 alone cannot read (two folders).
    mose_add_mpi_script(MPIExitCodes mpi/exit-codes "needs-mpi;implementation;1D")
    # input.ini read again during the run: rank 0 reads it and every rank parses
    # its text. sol-diter = 7 written in the input.ini of rank 0 alone (each
    # rank in its own folder, an MPMD launch) gives every rank the solutions of
    # iterations 21, 28 and 35; sol-diter, res-diter, ini-diter and cfl written
    # in the input.ini of rank 1 alone change nothing: exit code 0 within 60 s
    # and the files of the run without the change, byte for byte. A run whose
    # ranks took different values would hang in its collectives: the test stops
    # at the first case that fails, and its own TIMEOUT stops it in 120 s.
    mose_add_mpi_script(MPIRuntimeInput mpi/runtime-input "needs-mpi;implementation;2D")
    if(TEST MPIRuntimeInput)
        set_tests_properties(MPIRuntimeInput PROPERTIES TIMEOUT 120)
    endif()
endif()
