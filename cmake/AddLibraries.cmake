# ORION, OSLO and FLINT, as CMakeLists.txt required them with semver_dependency, and FiNeR (the submodule lib/FiNeR,
# with FiNeR's dependencies), each checked against the release MOSE requires, so that a clone of another release stops
# the configure before anything is built. Included only when MOSE is the top-level project: a parent project gives MOSE
# its libraries.

include(${CMAKE_CURRENT_LIST_DIR}/SubmodulePins.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/FiNeRDependencies.cmake)

# ORION, OSLO and FLINT: a release selected by SemVerDeps satisfies the requirement; an external clone (local) only
# warns there, so MOSE checks its version, the newest vX.Y.Z tag of its HEAD as SemVerDeps found it
foreach(_lib ORION OSLO FLINT)
  string(TOUPPER "${_lib}" _LIB)
  get_property(_source GLOBAL PROPERTY SEMVER_${_LIB}_SOURCE)
  if (_source STREQUAL "local")
    get_property(_version GLOBAL PROPERTY SEMVER_${_LIB}_VERSION)
    set(_ok FALSE)
    if (_version)
      _semver_satisfies("${_version}" _ok ${MOSE_${_LIB}_VERSION})
    else()
      set(_version "of no release (no vX.Y.Z tag in its history)")
    endif()
    if (NOT _ok)
      _semver_range("${MOSE_${_LIB}_VERSION}" _min _max)
      message(${MOSE_PIN_CHECK} "${_lib} at ${FETCHCONTENT_SOURCE_DIR_${_LIB}} is ${_version}, but MOSE requires "
                      "${_lib} ${MOSE_${_LIB}_VERSION} (${_min} up to ${_max} excluded). Give a clone of such a release, "
                      "or build the release MOSE selects (without FETCHCONTENT_SOURCE_DIR_${_LIB}), or configure with "
                      "-DMOSE_PIN_CHECK=WARNING to build this ${_lib} on purpose.")
    endif()
  endif()
endforeach()

# FiNeR: its dependencies before its CMakeLists adds them. Those of the submodule are fetched; an external
# clone is built as it is, its dependencies only checked. Its release is in its VERSION file (v2.3.0): its
# project() still declares 2.0.3
FetchContent_Declare(FiNeR SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/FiNeR")
if (FETCHCONTENT_SOURCE_DIR_FINER)
  set(MOSE_FINER_SOURCE_DIR "${FETCHCONTENT_SOURCE_DIR_FINER}")
  ensure_finer_dependencies("${MOSE_FINER_SOURCE_DIR}" FALSE)
else()
  set(MOSE_FINER_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/FiNeR")
  ensure_submodule_initialized("${MOSE_FINER_SOURCE_DIR}")
  mose_check_submodule_pin("${MOSE_FINER_SOURCE_DIR}")
  ensure_finer_dependencies("${MOSE_FINER_SOURCE_DIR}" TRUE)
endif()
set(_finer_version "")
if (EXISTS "${MOSE_FINER_SOURCE_DIR}/VERSION")
  file(STRINGS "${MOSE_FINER_SOURCE_DIR}/VERSION" _finer_version LIMIT_COUNT 1)
  string(REGEX REPLACE "^[ \t]*v" "" _finer_version "${_finer_version}")
endif()
string(REGEX REPLACE "^v" "" _finer_pinned "${MOSE_FINER_TAG}")
string(CONCAT _hint "Give FETCHCONTENT_SOURCE_DIR_FINER a clone of FiNeR ${MOSE_FINER_TAG}, or build the submodule "
                    "lib/FiNeR (without it); if lib/FiNeR was bumped on purpose, set MOSE_FINER_TAG in "
                    "CMakeLists.txt to its release")
mose_check_version(FiNeR "${MOSE_FINER_SOURCE_DIR}" "${_finer_version}" "${_finer_pinned}" "${_hint}")

# FLINT sets FLINT_SRC, the sources src/lib compiles, in this directory
FetchContent_MakeAvailable(OSLO ORION FLINT FiNeR)
