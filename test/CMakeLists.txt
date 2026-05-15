# ---------------------------------------------------------------------------
#  CTest definitions for MOSE validation cases
# ---------------------------------------------------------------------------

# Common environment for all solver tests
set(_MOSE_TEST_ENV
    "OMP_NUM_THREADS=1"
    "KMP_STACKSIZE=100M"
)

# ---- 1-D : Sod shock tube (1 mesh + verification) -------------
add_test(
    NAME    Sod79
    COMMAND bash -c "./MOSE.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/1D/Sod79"
)
set_tests_properties(Sod79 PROPERTIES
    LABELS      "1D;validation"
    TIMEOUT     300
    ENVIRONMENT "${_MOSE_TEST_ENV}"
)

# ---- 2-D : Prandtl-Meyer expansion (steady-state + verification) ----------
add_test(
    NAME    PrandtlMeyer
    COMMAND bash -c "./MOSE.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/prandtl-meyer"
)
set_tests_properties(PrandtlMeyer PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     300
    ENVIRONMENT "${_MOSE_TEST_ENV}"
)
