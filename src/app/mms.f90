!> @brief MOSE driver for the Method of Manufactured Solutions (MMS) viscous
!>        order-of-accuracy verification (test N2).
!>
!> Identical to the standard MOSE program (src/app/main.f90) except that the
!> external per-stage callback injects the analytic source term that makes a
!> chosen smooth, periodic field an EXACT steady solution of the 2D compressible
!> Navier-Stokes equations (calorically-perfect ideal gas, constant mu/k, Stokes
!> hypothesis).  Grid-refining and measuring the error against that field then
!> verifies the observed order of accuracy of the full (convective + viscous)
!> discretization -- including gradient/metric terms that the inviscid isentropic
!> vortex (N1) cannot exercise.
!>
!> MANUFACTURED PRIMITIVE FIELD (periodic on [0,1]^2, om = 2*pi):
!>     rho = 1.0   + 0.1 *sin(om x) cos(om y)
!>     u   = 40.0  + 8.0 *sin(om x) sin(om y)
!>     v   = 30.0  + 8.0 *cos(om x) cos(om y)
!>     p   = 8.0e3 + 8.0e2*cos(om x) sin(om y)   (low sound speed -> M~0.47)
!> Baked constants (MUST match INPUT/): R = Runiv/W = 8314.51/28.970418,
!> cp = 1004.5, gamma = 1.4, mu = 10.0 (flat transport.dat; Re~5), Pr = 0.72
!> (input.ini Prl), k = mu*cp/Pr.  The S_* expressions below are the analytic flux
!> divergence div(Fc - Fv) generated with sympy and cross-checked against a
!> finite-difference divergence (see test/2D/viscous/mms/mms_gen.py).
!>
!> Residual sign convention (as in Lib_RotatingFrame): V dU/dt = -R, with
!> R = CONV - DIFF - SOURCE, so a physical source S is added by  R -= S * vol.
program MOSE_mms_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use MOSE_Advanced_Types_m, only: MOSE_simulation_type, MOSE_domain_type
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

  if (mpi_is_root) write(*,'(A)') ' MMS driver: injecting manufactured-solution source terms'

  ! Solving with MOSE, injecting the manufactured source each RK stage
  call MOSE%setup( simulation )

  obj_sim_param%TODO = 1
  do while ( obj_sim_param%TODO <= 2 )
    call MOSE%solve( simulation, MMS_Source )
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

  !> @brief External per-stage callback: add the manufactured source to every
  !>        block's residual.  Mirrors RotatingFrame_Source_Terms.
  subroutine MMS_Source ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      call MMS_Source_Blk( domain%blk(b)%R,    &
                           domain%blk(b)%node, &
                           domain%blk(b)%vol,  &
                           domain%blk(b)%dim   )
    end do
  end subroutine MMS_Source

  !> @brief Block kernel: accumulate the analytic source into R (R -= S*vol).
  !>        Uses !$omp do because the caller (Explicit_Step) is already inside an
  !>        active parallel region -- otherwise every thread would add the source.
  subroutine MMS_Source_Blk ( Res, node, vol, n )
    use MOSE_Base_Types_m, only: R8, MOSE_vector_3D_type
    use MOSE_Global_m,     only: nu, nv, np, gc, nprim
    implicit none
    integer,  intent(in) :: n(3)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    type(MOSE_vector_3D_type), dimension(0:n(1), 0:n(2), 0:n(3)), intent(in) :: node
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
    ! Local
    integer  :: i, j, k
    real(R8) :: x, y, S_rho, S_mx, S_my, S_E
    real(R8) :: t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
    real(R8) :: t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38
    real(R8), parameter :: pi = 3.14159265358979323846d0

    !$omp do collapse(3) &
    !$omp private( i,j,k,x,y,S_rho,S_mx,S_my,S_E, &
    !$omp t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19, &
    !$omp t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38 )
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      ! Cell centre (average of the 8 hexahedral corner nodes; z unused in 2D)
      x = 0.125_R8 * ( node(i-1,j-1,k-1)%c(1) + node(i,j-1,k-1)%c(1) &
                     + node(i-1,j,  k-1)%c(1) + node(i,j,  k-1)%c(1) &
                     + node(i-1,j-1,k  )%c(1) + node(i,j-1,k  )%c(1) &
                     + node(i-1,j,  k  )%c(1) + node(i,j,  k  )%c(1) )
      y = 0.125_R8 * ( node(i-1,j-1,k-1)%c(2) + node(i,j-1,k-1)%c(2) &
                     + node(i-1,j,  k-1)%c(2) + node(i,j,  k-1)%c(2) &
                     + node(i-1,j-1,k  )%c(2) + node(i,j-1,k  )%c(2) &
                     + node(i-1,j,  k  )%c(2) + node(i,j,  k  )%c(2) )

      ! ---- analytic source (sympy cse, cross-checked vs FD divergence) ----
      t0 = 2*pi
      t1 = t0*x
      t2 = sin(t1)
      t3 = t0*y
      t4 = sin(t3)
      t5 = t2*t4
      t6 = 8.0d0*t5 + 40.0d0
      t7 = cos(t1)
      t8 = cos(t3)
      t9 = t7*t8
      t10 = t6*t9
      t11 = 8.0d0*t9 + 30.0d0
      t12 = t11*t5
      t13 = 640.0d0*pi
      t14 = t13*t5
      t15 = 0.2d0*t5 + 1
      t16 = t15**2
      t17 = t2*t8
      t18 = 0.1d0*t17
      t19 = t18 + 1.0d0
      t20 = t4*t7
      t21 = 1280.0d0*t15
      t22 = t20*t21
      t23 = 16.0d0*t19
      t24 = t23*t6
      t25 = t11*t23
      t26 = t13*t9
      t27 = 0.266666666666667d0*t9 + 1
      t28 = t27**2
      t29 = 960.0d0*t27
      t30 = t20*t29
      t31 = t4**2
      t32 = pi*t7**2
      t33 = t18 + 1
      t34 = pi/t33**2
      t35 = 800.0d0*t20 + 8000.0d0
      t36 = 3.8888888729863d0*t35/t33**3
      t37 = 160.0d0*t16 + 90.0d0*t28
      t38 = 0.05d0*t17 + 0.5d0
      S_rho = 0.2d0*pi*(t10 - t12)
      S_mx = pi*(-0.2d0*t12*t6 + t14 + 320.0d0*t16*t9 + t17*t25 + t19*t22 - &
      t20*t24 - 1600.0d0*t5)
      S_my = pi*(0.2d0*t10*t11 - t17*t24 - t19*t30 + t20*t25 + t26 - 180.0d0* &
      t28*t5 + 1600.0d0*t9)
      S_E = pi*(-pi*t2**2*t31*t36 + 89288888.8888889d0*pi*t20/( &
      28.7000001173611d0*t17 + 287.000001173611d0) + t11*t26 + t11*( &
      -t37*t5 + t38*(t17*t21 - t30) + 5599.99997710027d0*t9) + t14*t6 - &
      38.888888729863d0*t17*t34*t35 - 10240.0d0*t31*t32 - t32*t36*t8**2 &
      - 62222.2219677808d0*t34*t5*t9 + t6*(t37*t9 + t38*(-t17*t29 + t22 &
      ) - 5599.99997710027d0*t5))

      ! ---- accumulate:  R -= S * vol ----
      Res(1,  i, j, k) = Res(1,  i, j, k) - S_rho * vol(i,j,k)   ! continuity
      Res(nu, i, j, k) = Res(nu, i, j, k) - S_mx  * vol(i,j,k)   ! x-momentum
      Res(nv, i, j, k) = Res(nv, i, j, k) - S_my  * vol(i,j,k)   ! y-momentum
      Res(np, i, j, k) = Res(np, i, j, k) - S_E   * vol(i,j,k)   ! energy

    end do ; end do ; end do
    !$omp end do

  end subroutine MMS_Source_Blk

end program MOSE_mms_program
