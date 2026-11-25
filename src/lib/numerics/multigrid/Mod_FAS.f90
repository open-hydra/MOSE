module MOSE_Mod_FAS
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Lib_Limiters, only: rlimiter_if

  implicit none
  private
  public :: FAS_Solve

  !> Limiter of the running level, restored at the end of the cycle
  procedure(rlimiter_if), pointer :: limiter_fine => null()

contains

  !> FAS multigrid cycle on the levels below the running one, after its Explicit_Step (A. Jameson, "Solution of the
  !> Euler equations for two-dimensional transonic flow by a multigrid method", Appl. Math. Comput. 13, 1983).
  !> cycle = v visits each coarser level once, w twice. On the coarse levels the reconstruction is of first order and
  !> only the flow equations are solved, the RANS ones too with cycle-turbulence (Cycle_Turbulence)
  subroutine FAS_Solve ( simulation, External_Function )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_multigrid, obj_time_scheme
    use MOSE_Lib_Limiters,   only: rlimiter
    use MOSE_Mod_Timers,     only: timer_cycle_begin, timer_cycle_end
    use MOSE_Lib_Multigrid,  only: fas_level_iter
    implicit none
    type(MOSE_simulation_type), intent(inout) :: simulation
    external :: External_Function
    ! Local
    integer  :: level, mu, period
    real(R8) :: ramp

    level = obj_multigrid%MG_level
    if ( level == obj_multigrid%MGL ) return

    ! From the iteration cycle-start of the running level, one cycle every pre-iter + post-iter iterations, after its
    ! pre-iter ones
    if ( simulation%domain(level)%iter < obj_multigrid%cycle_start ) return
    period = obj_multigrid%pre_iter + obj_multigrid%post_iter
    if ( mod(simulation%domain(level)%iter, period) /= mod(obj_multigrid%pre_iter, period) ) return

    ! Iteration of the running level, for cycle-turbulence-start (Cycle_Turbulence)
    fas_level_iter = simulation%domain(level)%iter

    ! Visits of each coarser level: v 1, w 2
    mu = 0
    if ( trim(obj_multigrid%cycle) == 'w' ) mu = 1

    ! CFL ramp of the running level, applied to the coarse levels too
    ramp = 1d0
    if ( simulation%domain(level)%iter < obj_time_scheme%rampa_cfl_iter ) &
      ramp = real(simulation%domain(level)%iter, R8) / real(obj_time_scheme%rampa_cfl_iter, R8)

    ! The time of the cycle is counted in the time of the iteration
    call timer_cycle_begin()
    limiter_fine => rlimiter
    call FAS_Level ( simulation%domain, level, mu, ramp, External_Function )
    rlimiter => limiter_fine
    call timer_cycle_end()

  end subroutine FAS_Solve


  !> Level m of the cycle: pre-smoothing, restriction to level m+1, cycle of level m+1 (mu+1 times), correction of
  !> level m, post-smoothing. The running level is smoothed by Explicit_Step, the coarse levels by Smoothing_Step
  recursive subroutine FAS_Level ( domain, m, mu, ramp, External_Function )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_multigrid
    use MOSE_Mod_Explicit,   only: Smoothing_Step
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    integer,  intent(in) :: m, mu
    real(R8), intent(in) :: ramp
    external :: External_Function
    ! Local
    integer :: c, n, imu, next_mu, nsmooth

    c = m + 1

    ! Pre-smoothing: pre-iter iterations, coarse-iter on the coarsest level
    if ( m > obj_multigrid%MG_level ) then
      nsmooth = obj_multigrid%pre_iter
      if ( m == obj_multigrid%MGL ) nsmooth = obj_multigrid%coarse_iter
      call Set_Level_Limiter ( m )
      do n = 1, nsmooth
        call Smoothing_Step ( domain(m), ramp, External_Function )
      enddo
    endif

    if ( m == obj_multigrid%MGL ) return

    call FAS_To_Coarse ( domain, m, ramp, External_Function )

    ! Cycle of level c, mu+1 times (the coarsest level once)
    do imu = 0, mu
      next_mu = mu
      if ( m == obj_multigrid%MGL - 2 ) next_mu = 0
      call FAS_Level ( domain, c, next_mu, ramp, External_Function )
    enddo

    call FAS_Correct ( domain, m )

    ! Post-smoothing: post-iter iterations
    if ( m > obj_multigrid%MG_level ) then
      call Set_Level_Limiter ( m )
      do n = 1, obj_multigrid%post_iter
        call Smoothing_Step ( domain(m), ramp, External_Function )
      enddo
    endif

  end subroutine FAS_Level


  !> Restricted state (saved in fas_q0) and forcing term of level m+1
  subroutine FAS_To_Coarse ( domain, m, ramp, External_Function )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m,      only: obj_multigrid, obj_time_scheme
    use MOSE_Lib_Multigrid,       only: Restrict_State_Blk, fas_store_q0, fas_build_forcing, Add_FAS_Forcing, &
                                        Cycle_Turbulence
    use MOSE_Mod_Explicit,        only: Assemble_Residual
    use MOSE_Mod_dt,              only: Compute_dt
    use MOSE_Lib_Preconditioning, only: update_derived_variables
    use MOSE_Mod_MPI,             only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    integer,  intent(in) :: m
    real(R8), intent(in) :: ramp
    external :: External_Function
    ! Local
    integer :: b, c
    logical :: turb

    c = m + 1
    turb = Cycle_Turbulence()   ! RANS equations in the cycle

    ! Residual of level m for its smoothed state, forcing term included on a coarse level
    call Set_Level_Limiter ( m )
    !$omp parallel
    call Assemble_Residual ( domain(m), External_Function )
    if ( m > obj_multigrid%MG_level ) call Add_FAS_Forcing ( domain(m) )
    !$omp end parallel

    ! Restricted state of level c and its local time step (point-implicit factor of the turbulence sources)
    !$omp parallel
    do b = 1, domain(c) % nb
      if (.not. is_local_block(b)) cycle
      call Restrict_State_Blk ( domain(m)%blk(b)%P,   domain(c)%blk(b)%P,   &
                                domain(m)%blk(b)%vol, domain(c)%blk(b)%vol, &
                                domain(m)%blk(b)%dim, domain(c)%blk(b)%dim )
    enddo
    call Update_Derived_Variables ( domain(c) )
    !$omp end parallel
    domain(c) % dtglobal = 1d5
    call Compute_dt ( domain(c), obj_time_scheme%cfl * ramp, obj_time_scheme%vnn * ramp, 0 )

    ! Forcing term of level c
    call Set_Level_Limiter ( c )
    !$omp parallel
    do b = 1, domain(c) % nb
      if (.not. is_local_block(b)) cycle
      call fas_store_q0 ( domain(c)%blk(b)%P, domain(c)%blk(b)%fas_q0, domain(c)%blk(b)%dim, turb )
    enddo
    call Assemble_Residual ( domain(c), External_Function )
    do b = 1, domain(c) % nb
      if (.not. is_local_block(b)) cycle
      call fas_build_forcing ( domain(m)%blk(b)%R, domain(c)%blk(b)%R, domain(c)%blk(b)%fas_f, &
                               domain(m)%blk(b)%dim, domain(c)%blk(b)%dim, obj_multigrid%damp_restriction, turb )
    enddo
    !$omp end parallel
    domain(c) % fas_forcing = .true.

  end subroutine FAS_To_Coarse


  !> Correction of level m+1, prolongated to level m and added to its state. PO, not used between two steps, holds the
  !> prolongated correction
  subroutine FAS_Correct ( domain, m )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m,      only: obj_multigrid
    use MOSE_Lib_Multigrid,       only: fas_coarse_corr, coarse2fine_corr, coarse2fine_prim, Cycle_Turbulence
    use MOSE_Lib_Preconditioning, only: update_derived_variables
    use MOSE_Mod_MPI,             only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    integer, intent(in) :: m
    ! Local
    integer :: b, c
    logical :: turb

    c = m + 1
    turb = Cycle_Turbulence()   ! RANS equations in the cycle

    !$omp parallel
    do b = 1, domain(c) % nb
      if (.not. is_local_block(b)) cycle
      call fas_coarse_corr ( domain(c)%blk(b)%P, domain(c)%blk(b)%fas_q0, domain(c)%blk(b)%dim, turb )
    enddo
    do b = 1, domain(c) % nb
      if (.not. is_local_block(b)) cycle
      call coarse2fine_prim ( domain(m)%blk(b)%PO, domain(c)%blk(b)%fas_q0, domain(m)%blk(b)%dim, domain(c)%blk(b)%dim )
      call coarse2fine_corr ( domain(m)%blk(b)%P, domain(m)%blk(b)%PO, domain(m)%blk(b)%dim, &
                              obj_multigrid%damp_prolongation, b, turb )
    enddo
    call Update_Derived_Variables ( domain(m) )
    !$omp end parallel
    domain(c) % fas_forcing = .false.

  end subroutine FAS_Correct


  !> Limiter of level m: that of the run on the running level, first order on the coarse levels. Call it outside
  !> parallel regions
  subroutine Set_Level_Limiter ( m )
    use MOSE_Config_Types_m, only: obj_multigrid
    use MOSE_Lib_Limiters,   only: rlimiter, rlimiter_IORD
    implicit none
    integer, intent(in) :: m

    if ( m == obj_multigrid%MG_level ) then
      rlimiter => limiter_fine
    else
      rlimiter => rlimiter_IORD
    endif

  end subroutine Set_Level_Limiter

end module MOSE_Mod_FAS
