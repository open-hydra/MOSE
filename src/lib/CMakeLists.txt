file(GLOB_RECURSE MOSE_SRC
    CONFIGURE_DEPENDS
    *.f90
)

add_library(MOSEL ${MOSE_SRC} ${FLINT_SRC})
add_library(MOSEL::MOSEL ALIAS MOSEL)
target_link_libraries(MOSEL PRIVATE FiNeR::FiNeR ORION)
link_FLINT_dependencies(MOSEL)

target_include_directories(MOSEL INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    # Use the imported target so compile *and* link flags propagate
    # transitively to every consumer of MOSEL (PUBLIC).
    # The old SET_TARGET_PROPERTIES(LINK_FLAGS) approach is deprecated and
    # its LINK_FLAGS are silently ignored for static libraries by CMake.
    target_link_libraries(MOSEL PUBLIC OpenMP::OpenMP_Fortran)
ENDIF(USE_OPENMP)

IF(USE_MPI)
    target_compile_definitions(MOSEL PUBLIC USE_MPI)
    target_link_libraries(MOSEL PUBLIC MPI::MPI_Fortran)
ENDIF(USE_MPI)
