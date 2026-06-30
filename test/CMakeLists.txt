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

# 2-D steady Euler (analytic Prandtl-Meyer)
mose_add_test(NAME PrandtlMeyer DIR 2D/euler/prandtl-meyer
              RUN "./MOSE.sh solve" VERIFY check.py
              LABELS "fast;validation;2D")

# ===========================================================================
#  Additional fast tests defined in their own directories
#  (added incrementally — see test/TESTING-PLAN.md, section 4)
# ===========================================================================
# Riemann smoke-matrix (F5): exercises every Riemann solver for positivity.
if(EXISTS "${CMAKE_SOURCE_DIR}/test/fast/riemann-smoke/run.sh")
    add_test(
        NAME    RiemannSmoke
        COMMAND bash -c "./run.sh"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/fast/riemann-smoke"
    )
    set_tests_properties(RiemannSmoke PROPERTIES
        LABELS      "fast;numerics;1D"
        TIMEOUT     120
        ENVIRONMENT "${_MOSE_TEST_ENV}"
    )
endif()
