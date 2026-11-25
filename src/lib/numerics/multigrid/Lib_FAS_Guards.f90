module MOSE_Lib_FAS_Guards
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use MOSE_Finite_m,   only: Not_Finite

  implicit none
  private
  public :: Admissible_State, Halve_Coarse_Update, Report_FAS_Update_Fallbacks

  !> Coarse-level stage updates halved by this process (Halve_Coarse_Update)
  integer(I8) :: fas_update_fallbacks = 0

  !> Halvings of a coarse update or of a correction whose state is not physical, before it is dropped
  integer, parameter, public :: FAS_HALVINGS = 10

contains

  !> Coarse level of a multigrid cycle: a stage update whose state is not physical is halved, down to
  !> 2**-FAS_HALVINGS, then dropped, instead of stopping the run (the forcing term makes the coarse problem other than a
  !> physical evolution). On input prim holds the state of the full update
  subroutine Halve_Coarse_Update ( prim, cons, cons_stage, consold, update, temperature, irk, n_rk )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_RK
    implicit none
    real(R8), intent(inout) :: prim(nprim)
    real(R8), intent(in)    :: cons(nprim), cons_stage(nprim), consold(nprim), update(nprim), temperature
    integer,  intent(in)    :: irk, n_rk
    ! Local
    integer  :: n
    real(R8) :: theta, cons_h(nprim)

    if ( sum( cons(1:nsc) ) > 0d0 ) then
      if ( Admissible_State( prim ) ) return
    end if

    !$omp atomic
    fas_update_fallbacks = fas_update_fallbacks + 1

    theta = 1d0
    do n = 1, FAS_HALVINGS + 1
      theta = 0.5d0 * theta
      if ( n == FAS_HALVINGS + 1 ) theta = 0d0
      cons_h = RK_stage( irk, n_rk, cons_stage, consold, theta * update )
      if ( sum( cons_h(1:nsc) ) <= 0d0 ) cycle
      prim(1:np) = cons2prim( cons_h(1:np), temperature )
      if ( Admissible_State( prim ) ) exit
    enddo

  end subroutine Halve_Coarse_Update


  !> Physical state for the multigrid cycle: finite primitive variables, positive density and pressure, temperature within
  !> the thermodynamic tables (check_gas_state)
  logical function Admissible_State ( prim )
    use MOSE_Global_m, only: nprim, nsc, np
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in) :: prim(nprim)
    ! Local
    real(R8) :: rhoi(nsc)
    integer  :: error

    Admissible_State = .false.
    if ( any( Not_Finite(prim(1:np)) ) .or. sum( prim(1:nsc) ) <= 0d0 .or. prim(np) <= 0d0 ) return
    rhoi = prim(1:nsc)
    call check_gas_state ( rhoi, prim(np), error )
    Admissible_State = ( error == 0 )

  end function Admissible_State


  !> End of run: number of halved coarse updates over the ranks, if any. Collective
  subroutine Report_FAS_Update_Fallbacks ()
    use MOSE_Mod_MPI, only: mpi_is_root, mpi_allreduce_sum_r8
    implicit none
    real(R8) :: total

    call mpi_allreduce_sum_r8 ( real(fas_update_fallbacks, R8), total )
    if ( mpi_is_root .and. total > 0d0 ) &
      write(*,'(A,I0)') ' Multigrid coarse smoothing updates halved (updated state not physical): ', nint(total, I8)

  end subroutine Report_FAS_Update_Fallbacks

end module MOSE_Lib_FAS_Guards
