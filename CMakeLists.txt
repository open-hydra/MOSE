# CMake project file for 'MOSE'

cmake_minimum_required(VERSION 3.23)

# Name the code
set(solver MOSE)

project (${solver} LANGUAGES Fortran CXX C)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # seach path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # Set some options the user may choose
  # Uncomment the below if you want the user to choose a parallelization library
  OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
  OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)
  # Interprocedural (link-time) optimization. On by default because it is worth
  # ~10% -- see docs/development/performance_notes.md -- but switchable off for
  # toolchains where it misbehaves. It also self-disables if unsupported.
  OPTION(USE_IPO "Use interprocedural optimization in RELEASE builds" ON)

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING. You should  review this file and make sure the flags 
  # are to your liking.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
  # Locate and set parallelization libraries. There are some CMake peculiarities
  # taken care of here, such as the fact that the FindOpenMP routine doesn't know
  # about Fortran.
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # -------- Path fallbacks (standalone builds) --------
  if (NOT DEFINED ORION_PATH)
    set(ORION_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/ORION/")
  endif()
  if (NOT DEFINED FLINT_PATH)
    set(FLINT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/FLINT/")
  endif()
  if (NOT DEFINED FINER_PATH)
    set(FINER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/FiNeR/")
  endif()
  if (NOT DEFINED FINER_DEPS_PATH)
    set(FINER_DEPS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/")
  endif()

  # -------- Helper: initialize a git submodule from the correct git root --------
  macro(ensure_submodule_initialized SUBMOD_PATH)
    if (NOT EXISTS "${SUBMOD_PATH}/CMakeLists.txt")
      message(STATUS "Submodule at ${SUBMOD_PATH} not initialized. Attempting git submodule update --init...")
      get_filename_component(_sm_parent "${SUBMOD_PATH}" DIRECTORY)
      execute_process(
        COMMAND git rev-parse --show-toplevel
        WORKING_DIRECTORY "${_sm_parent}"
        OUTPUT_VARIABLE _sm_git_root
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _sm_git_result)
      if(_sm_git_result EQUAL 0)
        file(RELATIVE_PATH _sm_rel "${_sm_git_root}" "${SUBMOD_PATH}")
        execute_process(
          COMMAND git submodule update --init "${_sm_rel}"
          WORKING_DIRECTORY "${_sm_git_root}"
          RESULT_VARIABLE _sm_init_result)
        if(NOT _sm_init_result EQUAL 0)
          message(FATAL_ERROR "Failed to initialize submodule: ${SUBMOD_PATH}")
        endif()
      else()
        message(FATAL_ERROR "Could not find git root for submodule: ${SUBMOD_PATH}")
      endif()
    endif()
  endmacro()

  # ORION
  if (NOT TARGET ORION)
    ensure_submodule_initialized("${ORION_PATH}")
    set(orionD "${ORION_PATH}" "${CMAKE_BINARY_DIR}/lib/ORION")
    add_subdirectory(${orionD})
  endif()

  # FLINT
  if (NOT TARGET FLINT)
    ensure_submodule_initialized("${FLINT_PATH}")
    set(flintD "${FLINT_PATH}" "${CMAKE_BINARY_DIR}/lib/FLINT")
    add_subdirectory(${flintD})
  endif()

  # Libraries of FiNeR: FiNeR 2.0.7 and later do not carry them as submodules (its CMake build expects them in its
  # src/third_party/, fetched with FoBiS) and add each one only when the target is not defined yet: MOSE provides
  # them, in the order of their dependencies (FLAP needs PENF and FACE, StringiFor needs PENF and BeFoR64)
  foreach (finer_dep PENF FACE FLAP BeFoR64 StringiFor)
    if (NOT TARGET ${finer_dep})
      ensure_submodule_initialized("${FINER_DEPS_PATH}/${finer_dep}")
      add_subdirectory("${FINER_DEPS_PATH}/${finer_dep}" "${CMAKE_BINARY_DIR}/lib/third_party/${finer_dep}")
    endif()
  endforeach()

  # FiNeR
  if (NOT TARGET FiNeR)
    ensure_submodule_initialized("${FINER_PATH}")
    # MOSE needs FiNeR v2.1.0 or later: before it, a comment line that holds "=" or a bracket was taken for an option or
    # a section header, and the reading of its section of input.ini ended there. A clone that pulls without
    # "git submodule update" keeps the FiNeR it had checked out (v2.0.6 before v2.3.0): ensure_submodule_initialized
    # leaves it as it is, and it would build in silence. FiNeR names its version in a VERSION file ("vX.Y.Z") since
    # v2.0.7; a FiNeR without it is older
    set(finer_version "")
    if (EXISTS "${FINER_PATH}/VERSION")
      file(STRINGS "${FINER_PATH}/VERSION" finer_version LIMIT_COUNT 1)
      string(REGEX REPLACE "^[ \t]*v?([0-9][0-9.]*).*$" "\\1" finer_version "${finer_version}")
    endif()
    if (NOT "${finer_version}" VERSION_GREATER_EQUAL "2.1.0")
      if ("${finer_version}" STREQUAL "")
        set(finer_found "a FiNeR without a VERSION file, older than v2.0.7")
      else()
        set(finer_found "FiNeR v${finer_version}")
      endif()
      message(FATAL_ERROR "FiNeR at ${FINER_PATH} is ${finer_found}: MOSE needs FiNeR v2.1.0 or later (it pins v2.3.0), "
                          "as before v2.1.0 a comment line of input.ini that holds \"=\" ends the reading of its section. "
                          "Run: git submodule update --init --recursive")
    endif()
    set(finerD "${FINER_PATH}" "${CMAKE_BINARY_DIR}/lib/third_party/FiNeR")
    add_subdirectory(${finerD})
  endif()

endif()

# Build library
add_subdirectory(src/lib)

# Build application (only if this file has not been called by another CMakeLists)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(src/app)

  # Testing
  enable_testing()
  add_subdirectory(test)
endif()
