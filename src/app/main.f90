program MOSE_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use MOSE_Advanced_Types_m, only: MOSE_simulation_type
  use MOSE_Config_Types_m,   only: obj_sim_param
  use MOSE_Procedures_m,     only: MOSE_type
  use MOSE_Mod_MPI
#ifdef USE_MPI
  use MOSE_Mod_GhostExchange, only: cleanup_ghost_schedule
#endif
  implicit none
  type(MOSE_type)            :: MOSE
  type(MOSE_simulation_type) :: simulation

  ! Initialize MPI environment (no-op if USE_MPI is not defined)
  call mpi_init_env()

#if defined (_OPENMP)
  !$omp parallel
  obj_sim_param%nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  if (mpi_is_root) then
    write(*,'(A)')    ' Parallel execution'
    write(*,'(A)')    ' OpenMP:'
    write(*,'(A,I4)') ' -  Number of threads --> ', obj_sim_param%nthreads
  end if
#else
  if (mpi_is_root) write(*,'(A)')    ' Serial execution'
  obj_sim_param%nthreads = 1
#endif

#ifdef USE_MPI
  if (mpi_is_root) then
    write(*,'(A)')      ' MPI:'
    write(*,'(A,I4)')   ' -  Number of ranks   --> ', mpi_size_
  end if
#endif

  ! Solving with MOSE
  call MOSE%setup( simulation )

  obj_sim_param%TODO = 1
  do while ( obj_sim_param%TODO <= 2 )
    call MOSE%solve( simulation, Dummy_Function )
    if ( obj_sim_param%TODO <= 2 ) call MOSE%postprocess( simulation )
  enddo

  call MOSE%postprocess( simulation )

  ! Free persistent MPI requests before finalizing
#ifdef USE_MPI
  call cleanup_ghost_schedule()
#endif

  ! Finalize MPI environment
  call mpi_finalize_env()

contains

    subroutine Dummy_Function
      ! Empty subroutine to be passed as an argument to MOSE%solve. 
      ! It can be used for user-defined operations during the solution process.
      ! It is used in HYDRA coupling procedures.
    end subroutine Dummy_Function

end program MOSE_program