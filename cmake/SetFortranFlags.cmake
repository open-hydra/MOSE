######################################################
# Determine and set the Fortran compiler flags we want 
######################################################

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
INCLUDE(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)

# When using IPO with Intel + MPI/OpenMP, the linker must also receive -ipo.
# CMake does not propagate compile-time IPO flags to the linker automatically.
if(USE_MPI OR USE_OPENMP)
    SET_COMPILE_FLAG(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE}"
                     Fortran "-ipo"
                             "/Qipo"
                    )
endif()

# Make sure the build type is uppercase
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "TESTING")
    SET (CMAKE_BUILD_TYPE TESTING CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG, RELEASE, or TESTING")
ENDIF(BT STREQUAL "RELEASE")

#########################################################
# If the compiler flags have already been set, return now
#########################################################

IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)
    RETURN ()
ENDIF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)

########################################################################
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED 
# flag is given in the call).  This way unknown compiles are supported.
#######################################################################

#####################
### GENERAL FLAGS ###
#####################

# Don't add underscores in symbols for C-compatability
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-underscoring")

# There is some bug where -march=native doesn't work on Mac
IF(APPLE)
    SET(GNUNATIVE "-mtune=native")
ELSE()
    SET(GNUNATIVE "-march=native")
ENDIF()
# Optimize for the host's architecture
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-xHost"        # Intel
                         "/QxHost"       # Intel Windows
                         ${GNUNATIVE}    # GNU
                         "-ta=host"      # Portland Group
                )
# Add preprocessor flag
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-cpp"        # Intel or GNU
                )

SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                  Fortran  "-ffree-line-length-none"     # Intel o GNU
                           "-ffixed-line-length-none"    # Intel
                           "-extend-source"
                )        
###################
### DEBUG FLAGS ###
###################

# Debugging symbols (explicit, do not rely on CMake defaults)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-ggdb3"  # GNU -- max GDB info including macros
                         "-g3"     # generic fallback
                         "-g"      # Intel/PGI
                )

# Disable optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran REQUIRED "-O0" # All compilers not on Windows
                                  "/Od" # Intel Windows
                )

# Turn on all warnings
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-warn all"         # Intel
                         "/warn:all"         # Intel Windows
                         "-Wall"             # GNU
                                             # Portland Group (on by default)
                )

# Extra warnings (implicit interfaces, unused dummies, etc.)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-Wextra"                 # GNU
                         "-Wimplicit-interface"     # GNU
                         "-Wimplicit-procedure"     # GNU
                         "-warn interfaces"         # Intel
                )

# Traceback
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-traceback"   # Intel/Portland Group
                         "/traceback"   # Intel Windows
                         "-fbacktrace"  # GNU (gfortran)
                         "-ftrace=full" # GNU (g95)
                )

# Check everything: bounds, pointers, uninitialized, temporaries
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-check all"    # Intel
                         "/check:all"    # Intel Windows
                         "-fcheck=all"   # GNU (New style) -- replaces -fcheck=bounds
                         "-fbounds-check" # GNU (Old style fallback)
                         "-Mbounds"       # Portland Group
                )

# Initialize reals to signalling NaN to catch use-before-set
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-finit-real=snan"    # GNU
                         "-init=snan,arrays"   # Intel
                )

# Initialize integers to a sentinel value
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-finit-integer=-42"  # GNU
                )

# Trap floating-point exceptions (NaN/Inf/zero-divide crash immediately)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-ffpe-trap=invalid,zero,overflow"  # GNU
                         "-fpe0"                             # Intel (trap all)
                )

# Stack protection
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-fstack-protector-all"  # GNU
                )

#####################
### TESTING FLAGS ###
#####################

# Optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}"
                 Fortran REQUIRED "-O2" # All compilers not on Windows
                                  "/O2" # Intel Windows
                )

# Debug symbols so tests are still debuggable
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}"
                 Fortran "-g"  # GNU/Intel/PGI
                )

# Bounds checking (compatible with -O2)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}"
                 Fortran "-fcheck=bounds"  # GNU
                         "-check bounds"   # Intel
                )

#####################
### RELEASE FLAGS ###
#####################

# Explicit -O3 (do not rely on CMake defaults)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran REQUIRED "-O3"  # GNU/Intel
                                  "/O3"  # Intel Windows
                )

# Unroll loops
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-funroll-loops" # GNU
                         "-unroll"        # Intel
                         "/unroll"        # Intel Windows
                         "-Munroll"       # Portland Group
                )

# Inline functions
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-inline"            # Intel
                         "/Qinline"           # Intel Windows
                         "-finline-functions" # GNU
                         "-Minline"           # Portland Group
                )

# Interprocedural optimization (Intel/PGI only)
# GNU's -flto is excluded: it is fragile with mixed-language static libraries,
# MPI wrappers, and certain linker configurations.
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-ipo"              # Intel
                         "/Qipo"             # Intel Windows
                         "-Mipa=fast,inline" # Portland Group
                )

# Fast math (breaks strict IEEE 754 -- verify numerics are unaffected)
# WARNING: may change floating-point results; disable if results diverge
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-ffast-math"        # GNU
                         "-fp-model fast=2"   # Intel (aggressive)
                )

# Free a register by omitting frame pointer
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-fomit-frame-pointer"  # GNU
                )

#########################################
### DIAGNOSTIC SUMMARY (always prints) ##
#########################################
message(STATUS "----------------------------------------------------")
message(STATUS "Fortran compiler   : ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} (${CMAKE_Fortran_COMPILER})")
message(STATUS "Build type         : ${CMAKE_BUILD_TYPE}")
message(STATUS "Fortran base flags : ${CMAKE_Fortran_FLAGS}")
if(BT STREQUAL "DEBUG")
  message(STATUS "Fortran DEBUG flags: ${CMAKE_Fortran_FLAGS_DEBUG}")
elseif(BT STREQUAL "TESTING")
  message(STATUS "Fortran TEST flags : ${CMAKE_Fortran_FLAGS_TESTING}")
else()
  message(STATUS "Fortran RELEASE flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()
message(STATUS "USE_OPENMP=${USE_OPENMP}  USE_MPI=${USE_MPI}")
message(STATUS "----------------------------------------------------")

