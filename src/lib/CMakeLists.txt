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
    # An mpi module paired with a compiler that cannot ignore the types of choice buffers declares no
    # interfaces for the routines with choice buffers (MPI 4.1, Section 19.1.3), as the gfortran modules of
    # Intel MPI do: gfortran >= 10 then stops on calls of one such routine with buffers of different types
    # or ranks. -fallow-argument-mismatch turns those errors into warnings; it is set in MPI builds only,
    # and on Mod_MPI.f90 only, the one file that makes such calls: the builds without MPI, and the other
    # files of MOSE and FLINT in the MPI builds, keep the full check.
    IF(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
        set_source_files_properties(parallel/Mod_MPI.f90 PROPERTIES COMPILE_OPTIONS -fallow-argument-mismatch)
    ENDIF()
ENDIF(USE_MPI)
