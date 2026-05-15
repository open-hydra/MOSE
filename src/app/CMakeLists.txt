function(add_fortran_executable target source)
  
  add_executable(${target} ${source})
  set_target_properties(${target} PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(${target} PRIVATE MOSEL FiNeR::FiNeR)

  if (USE_TECIO)
    target_link_libraries(${target} PRIVATE -pthread)
  endif()

  target_include_directories(${target} PRIVATE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

  if (USE_OPENMP)
    # OpenMP::OpenMP_Fortran carries both compile and link flags.
    # MOSEL already links it as PUBLIC so this is transitive, but
    # being explicit here documents the dependency clearly.
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  if (USE_MPI)
    target_compile_definitions(${target} PRIVATE USE_MPI)
    target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)
  endif()

endfunction()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_fortran_executable(${solver} "main.f90")
add_fortran_executable(DocGen "docgen.f90")
