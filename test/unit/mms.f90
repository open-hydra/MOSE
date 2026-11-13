!> @brief State of the manufactured-solution driver read by its source callback, kept in a module rather
!>        than in the program: the callback is an internal procedure passed to the solver as an argument,
!>        and reading an allocatable array of the program makes ifx build a trampoline on the stack (an
!>        executable stack in the binary); a module variable needs no access to the frame of the program.
module MOSE_MMS_Driver_State
  use MOSE_Base_Types_m, only: R8
  implicit none
  real(R8), allocatable :: smean(:,:)   ! mean of S removed from the source (equation, grid level)
end module MOSE_MMS_Driver_State


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
!> Baked constants (MUST match INPUT/): R = Runiv/W = 8314.46261815324/28.970418,
!> cp = 1004.5, gamma = 1.4, mu = 10.0 (flat transport.dat; Re~5), Pr = 0.72
!> (input.ini Prl), k = mu*cp/Pr.  The S_* expressions below are the analytic flux
!> divergence div(Fc - Fv) generated with sympy and cross-checked against a
!> finite-difference divergence (see test/2D/viscous/mms/mms_gen.py).
!>
!> Other fields, chosen with the environment variable MOSE_MMS_FIELD (default:
!> periodic2d, the field above):
!>   periodic3d  the field above extended in z (periodic on [0,1]^3), with
!>               w = 20 + 8 cos(om x) sin(om y) sin(om z) and cos(om z) factors;
!>   wall2d      periodic in x, walls at y = 0 and y = 1: u = sin(pi y)(40 + 8 sin(om x)),
!>               v = 8 sin(pi y) cos(om x), T = 28 + 1.4 cos(om x) + 5.6 sin(pi y),
!>               p = 8000 + 800 cos(om x) cos(pi y), rho = p/(R T): zero velocity and
!>               T = T_w(x) on the walls, d2u/dy2 = d2T/dy2 = 0 there;
!>   wall3d      the same, periodic in x and z, with cos(om z) factors and
!>               w = sin(pi y)(20 + 8 cos(om x) sin(om z)).
!> Their sources come from test/2D/viscous/mms/mms_gen_fields.py (sympy, checked
!> against a finite-difference divergence at 20 random points to 1e-6).
!>
!> Mean of the discrete source.  The face fluxes of a conservative scheme cancel in
!> pairs, so for an equation whose total boundary flux is zero for any state the sum
!> of the residuals over the cells is minus the sum of the source, S(c_i)*V_i: the
!> discrete steady problem is solvable only if that sum is zero.  The source at the
!> cell centres sums to zero only to the accuracy of the midpoint rule, so for those
!> equations its mean over the domain is removed, once at setup: all of them with
!> periodic2d and periodic3d (periodic in every direction; the symmetry k faces of a
!> quasi-2D run carry no flux of these equations but the z-momentum, whose residual
!> is zeroed), the mass alone with wall2d and wall3d (wall stress, pressure and heat
!> flux balance the other ones).  The mean and its ratio to the mean of |S| are
!> written to the log.
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
  use MOSE_Base_Types_m,     only: R8
  use MOSE_MMS_Driver_State, only: smean
  use MOSE_Mod_MPI
#ifdef USE_MPI
  use MOSE_Mod_GhostExchange, only: cleanup_ghost_schedule
#endif
  implicit none
  type(MOSE_type)            :: MOSE
  type(MOSE_simulation_type) :: simulation
  character(len=64)          :: mms_field = 'periodic2d'
  character(len=64)          :: env_value
  integer                    :: env_status

  ! Initialize MPI environment (no-op if USE_MPI is not defined)
  call mpi_init_env()

#if defined (_OPENMP)
  ! The size of the team: one thread writes it (every thread of the team writing it would be a data race)
  !$omp parallel
  !$omp single
  obj_sim_param%nthreads = OMP_GET_NUM_THREADS()
  !$omp end single
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

  call get_environment_variable('MOSE_MMS_FIELD', value=env_value, status=env_status)
  if (env_status == 0) then
    mms_field = adjustl(env_value)
  else if (env_status /= 1) then
    ! 1: the variable is not set (default field); -1: value longer than env_value; other: not readable
    if (mpi_is_root) write(*,'(A)') ' MMS driver: MOSE_MMS_FIELD could not be read'
    error stop 1
  end if
  select case (trim(mms_field))
    case ('periodic2d', 'periodic3d', 'wall2d', 'wall3d')
    case default
      if (mpi_is_root) write(*,'(A)') ' MMS driver: unknown MOSE_MMS_FIELD '//trim(mms_field)// &
                                     ' (periodic2d, periodic3d, wall2d, wall3d)'
      error stop 1
  end select
  if (mpi_is_root) write(*,'(A)') ' MMS driver: injecting manufactured-solution source terms, field '//trim(mms_field)

  ! Solving with MOSE, injecting the manufactured source each RK stage
  call MOSE%setup( simulation )
  call MMS_Source_Mean ( simulation )

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
      call MMS_Field_Blk( domain%blk(b)%R,    &
                          domain%blk(b)%node, &
                          domain%blk(b)%vol,  &
                          domain%blk(b)%dim   )
      call MMS_Remove_Mean_Blk( domain%blk(b)%R,   &
                                domain%blk(b)%vol, &
                                domain%blk(b)%dim, &
                                smean(:, domain%mg_level) )
    end do
  end subroutine MMS_Source

  !> @brief Source of the chosen field on one block: R -= S*vol (MMS_Source_Blk for periodic2d,
  !>        MMS_Source_Blk_Field for the others).
  subroutine MMS_Field_Blk ( Res, node, vol, n )
    use MOSE_Base_Types_m, only: R8, MOSE_vector_3D_type
    use MOSE_Global_m,     only: gc, nprim
    implicit none
    integer,  intent(in) :: n(3)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    type(MOSE_vector_3D_type), dimension(0:n(1), 0:n(2), 0:n(3)), intent(in) :: node
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol

    if (trim(mms_field) == 'periodic2d') then
      call MMS_Source_Blk ( Res, node, vol, n )
    else
      call MMS_Source_Blk_Field ( Res, node, vol, n )
    endif
  end subroutine MMS_Field_Blk

  !> @brief Removes the mean of the discrete source (MMS_Source_Mean) on one block: R += sm*vol, so that the
  !>        two calls add -(S - sm)*vol; sm is zero for the equations whose mean is kept.
  subroutine MMS_Remove_Mean_Blk ( Res, vol, n, sm )
    use MOSE_Base_Types_m, only: R8
    use MOSE_Global_m,     only: nu, nv, nw, np, gc, nprim
    implicit none
    integer,  intent(in) :: n(3)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
    real(R8), intent(in) :: sm(5)
    ! Local
    integer :: i, j, k

    !$omp do collapse(3) private( i, j, k )
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      Res(1,  i, j, k) = Res(1,  i, j, k) + sm(1) * vol(i,j,k)   ! continuity
      Res(nu, i, j, k) = Res(nu, i, j, k) + sm(2) * vol(i,j,k)   ! x-momentum
      Res(nv, i, j, k) = Res(nv, i, j, k) + sm(3) * vol(i,j,k)   ! y-momentum
      Res(nw, i, j, k) = Res(nw, i, j, k) + sm(4) * vol(i,j,k)   ! z-momentum
      Res(np, i, j, k) = Res(np, i, j, k) + sm(5) * vol(i,j,k)   ! energy
    end do ; end do ; end do
    !$omp end do

  end subroutine MMS_Remove_Mean_Blk

  !> @brief Mean of the discrete source, sum_i S(c_i)*V_i / sum_i V_i, on every grid level, kept in smean
  !>        for the equations whose total boundary flux is zero (see the head of this file): all of them with
  !>        periodic2d and periodic3d, the mass alone with wall2d and wall3d.  The source of a block is that of
  !>        the kernel of the run (MMS_Field_Blk on a zeroed residual); the cells of a block are summed in a
  !>        fixed order, the blocks reduced over the ranks one by one and summed in their order, so the mean
  !>        does not depend on the decomposition or on the threads.  Called once, after setup, outside any
  !>        parallel region (the kernels then run on one thread).
  subroutine MMS_Source_Mean ( simulation )
    use MOSE_Global_m, only: nu, nv, nw, np
    implicit none
    type(MOSE_simulation_type), intent(in) :: simulation
    ! Local
    integer  :: l, b, i, j, k, q, nb, iq(5)
    logical  :: keep(5)
    real(R8) :: vsum, ssum(5), asum(5)
    real(R8), allocatable :: bs(:,:), ba(:,:), bv(:), work(:,:,:,:)
    character(len=*), parameter :: eq(5) = [ 'mass      ', 'x-momentum', 'y-momentum', 'z-momentum', 'energy    ' ]

    iq = [ 1, nu, nv, nw, np ]
    select case (trim(mms_field))
      case ('wall2d', 'wall3d'); keep = [ .true., .false., .false., .false., .false. ]
      case default;              keep = .true.        ! periodic2d, periodic3d
    end select
    allocate ( smean(5, size(simulation%domain)) )
    smean = 0.0_R8

    do l = 1, size(simulation%domain)
      associate ( domain => simulation%domain(l) )
        nb = domain%nb
        allocate ( bs(5, nb), ba(5, nb), bv(nb) )
        bs = 0.0_R8 ; ba = 0.0_R8 ; bv = 0.0_R8
        do b = 1, nb
          if (.not. is_local_block(b)) cycle
          associate ( blk => domain%blk(b) )
            allocate ( work, mold = blk%R )
            work = 0.0_R8
            call MMS_Field_Blk ( work, blk%node, blk%vol, blk%dim )
            do k = 1, blk%dim(3)
            do j = 1, blk%dim(2)
            do i = 1, blk%dim(1)
              do q = 1, 5
                bs(q, b) = bs(q, b) - work(iq(q), i, j, k)
                ba(q, b) = ba(q, b) + abs(work(iq(q), i, j, k))
              end do
              bv(b) = bv(b) + blk%vol(i,j,k)
            end do ; end do ; end do
            deallocate ( work )
          end associate
        end do
        call mpi_allreduce_sum_r8_array ( bs, 5*nb )
        call mpi_allreduce_sum_r8_array ( ba, 5*nb )
        call mpi_allreduce_sum_r8_array ( bv, nb )
        ssum = 0.0_R8 ; asum = 0.0_R8 ; vsum = 0.0_R8
        do b = 1, nb
          ssum = ssum + bs(:, b)
          asum = asum + ba(:, b)
          vsum = vsum + bv(b)
        end do
        where (keep) smean(:, l) = ssum / vsum
        if (mpi_is_root) then
          do q = 1, 5
            if (asum(q) == 0.0_R8) cycle
            write(*,'(A,I0,A,A,A,ES12.4,A,ES10.2,A)') ' MMS driver: level ', l, ', ', eq(q), &
              ': mean of the source ', ssum(q) / vsum, ', ratio to the mean of |S| ', ssum(q) / asum(q), &
              merge(', removed', ', kept   ', keep(q))
          end do
        end if
        deallocate ( bs, ba, bv )
      end associate
    end do
  end subroutine MMS_Source_Mean

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
    real(R8) :: t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
    real(R8), parameter :: pi = 3.14159265358979323846d0

    !$omp do collapse(3) &
    !$omp private( i,j,k,x,y,S_rho,S_mx,S_my,S_E, &
    !$omp t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19, &
    !$omp t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39 )
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
      t34 = t33**(-2)
      t35 = 800.0d0*t20 + 8000.0d0
      t36 = pi*t35
      t37 = 3.88891103469958d0/t33**3
      t38 = 160.0d0*t16 + 90.0d0*t28
      t39 = 0.05d0*t17 + 0.5d0
      S_rho = 0.2d0*pi*(t10 - t12)
      S_mx = pi*(-0.2d0*t12*t6 + t14 + 320.0d0*t16*t9 + t17*t25 + t19*t22 - &
      t20*t24 - 1600.0d0*t5)
      S_my = pi*(0.2d0*t10*t11 - t17*t24 - t19*t30 + t20*t25 + t26 - 180.0d0* &
      t28*t5 + 1600.0d0*t9)
      S_E = pi*(89288888.8888889d0*pi*t20/(28.6998365648478d0*t17 + &
      286.998365648478d0) - 62222.5765551933d0*pi*t34*t5*t9 + t11*t26 + &
      t11*(-t38*t5 + t39*(t17*t21 - t30) + 5600.0318899674d0*t9) + t14* &
      t6 - 38.8891103469958d0*t17*t34*t36 - t2**2*t31*t36*t37 - &
      10240.0d0*t31*t32 - t32*t35*t37*t8**2 + t6*(t38*t9 + t39*(-t17* &
      t29 + t22) - 5600.0318899674d0*t5))

      ! ---- accumulate:  R -= S * vol ----
      Res(1,  i, j, k) = Res(1,  i, j, k) - S_rho * vol(i,j,k)   ! continuity
      Res(nu, i, j, k) = Res(nu, i, j, k) - S_mx  * vol(i,j,k)   ! x-momentum
      Res(nv, i, j, k) = Res(nv, i, j, k) - S_my  * vol(i,j,k)   ! y-momentum
      Res(np, i, j, k) = Res(np, i, j, k) - S_E   * vol(i,j,k)   ! energy

    end do ; end do ; end do
    !$omp end do

  end subroutine MMS_Source_Blk


  !> @brief Block kernel of the fields other than periodic2d: R -= S*vol for the five
  !>        conservation equations (w-momentum included: in a quasi-2D run with symmetry
  !>        k faces its residual is zeroed afterwards, Mod_Explicit.f90).
  subroutine MMS_Source_Blk_Field ( Res, node, vol, n )
    use MOSE_Base_Types_m, only: R8, MOSE_vector_3D_type
    use MOSE_Global_m,     only: nu, nv, nw, np, gc, nprim
    implicit none
    integer,  intent(in) :: n(3)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    type(MOSE_vector_3D_type), dimension(0:n(1), 0:n(2), 0:n(3)), intent(in) :: node
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
    ! Local
    integer  :: i, j, k, d
    real(R8) :: c(3), S(5)

    !$omp do collapse(3) private( i, j, k, d, c, S )
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      ! Cell centre (average of the 8 hexahedral corner nodes)
      do d = 1, 3
        c(d) = 0.125_R8 * ( node(i-1,j-1,k-1)%c(d) + node(i,j-1,k-1)%c(d) &
                          + node(i-1,j,  k-1)%c(d) + node(i,j,  k-1)%c(d) &
                          + node(i-1,j-1,k  )%c(d) + node(i,j-1,k  )%c(d) &
                          + node(i-1,j,  k  )%c(d) + node(i,j,  k  )%c(d) )
      enddo
      select case (trim(mms_field))
        case ('periodic3d'); call Source_periodic3d ( c(1), c(2), c(3), S )
        case ('wall2d');     call Source_wall2d     ( c(1), c(2), c(3), S )
        case ('wall3d');     call Source_wall3d     ( c(1), c(2), c(3), S )
      end select
      Res(1,  i, j, k) = Res(1,  i, j, k) - S(1) * vol(i,j,k)   ! continuity
      Res(nu, i, j, k) = Res(nu, i, j, k) - S(2) * vol(i,j,k)   ! x-momentum
      Res(nv, i, j, k) = Res(nv, i, j, k) - S(3) * vol(i,j,k)   ! y-momentum
      Res(nw, i, j, k) = Res(nw, i, j, k) - S(4) * vol(i,j,k)   ! z-momentum
      Res(np, i, j, k) = Res(np, i, j, k) - S(5) * vol(i,j,k)   ! energy
    end do ; end do ; end do
    !$omp end do

  end subroutine MMS_Source_Blk_Field

  !> @brief Analytic source of the field periodic3d (mms_gen_fields.py, sympy cse).
  pure subroutine Source_periodic3d ( x, y, z, S )
    use MOSE_Base_Types_m, only: R8
    implicit none
    real(R8), intent(in)  :: x, y, z
    real(R8), intent(out) :: S(5)
    real(R8), parameter   :: pi = 3.14159265358979323846d0
    real(R8) :: S_rho, S_mx, S_my, S_mz, S_E
    real(R8) :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9
    real(R8) :: t10, t11, t12, t13, t14, t15, t16, t17, t18, t19
    real(R8) :: t20, t21, t22, t23, t24, t25, t26, t27, t28, t29
    real(R8) :: t30, t31, t32, t33, t34, t35, t36, t37, t38, t39
    real(R8) :: t40, t41, t42, t43, t44, t45, t46, t47, t48, t49
    real(R8) :: t50, t51, t52, t53, t54, t55, t56, t57, t58, t59

        t0 = 2*pi
        t1 = t0*x
        t2 = sin(t1)
        t3 = t0*z
        t4 = cos(t3)
        t5 = t0*y
        t6 = sin(t5)
        t7 = t4*t6
        t8 = t2*t7
        t9 = t8 + 5
        t10 = cos(t1)
        t11 = cos(t5)
        t12 = t11*t4
        t13 = t10*t12
        t14 = t13*t9
        t15 = t12*t2
        t16 = 0.1d0*t15 + 1
        t17 = t10*t7
        t18 = 16*t17
        t19 = sin(t3)
        t20 = t10*t19
        t21 = t20*t6
        t22 = 2*t21 + 5
        t23 = t19*t2
        t24 = t11*t23
        t25 = t22*t24
        t26 = 4*t13 + 15
        t27 = 0.4d0*t8
        t28 = 1066.66666666667d0*pi
        t29 = t9**2
        t30 = t16*t9
        t31 = 3.2d0*t26
        t32 = t8*t9
        t33 = t23*t6
        t34 = t16*t33
        t35 = 64*t22
        t36 = 32*t16
        t37 = t15*t26
        t38 = pi*t13
        t39 = t26**2
        t40 = t11*t20
        t41 = 1.6d0*t26
        t42 = t22**2
        t43 = t26*t40
        t44 = pi*t6**2
        t45 = t2**2
        t46 = t19**2
        t47 = t45*t46
        t48 = t4**2
        t49 = t10**2*t48
        t50 = t16**(-2)
        t51 = pi*t10*t11*t2*t50*t6
        t52 = t17 + 10
        t53 = pi*t52
        t54 = 3111.12882775967d0/t16**3
        t55 = t11**2*t53*t54
        t56 = 16*t29 + t39 + 4*t42
        t57 = 0.4d0*t56
        t58 = 4*t9
        t59 = 2*t22
        S_rho = pi*(1.6d0*t14 + t16*t18 - 0.8d0*t25 - t26*t27)
        S_mx = pi*(12.8d0*t13*t29 + 256*t17*t30 - 6.4d0*t25*t9 + t28*t8 - t31* &
        t32 - t34*t35 + t36*t37 - 1600*t8)
        S_my = pi*(1600*t13 + t14*t31 - 128*t15*t30 - t16*t35*t40 - t25*t41 + &
        853.333333333333d0*t38 - 0.8d0*t39*t8)
        S_mz = pi*(6.4d0*t14*t22 + 128*t16*t17*t22 + t21*t28 - 1600*t21 - t22* &
        t41*t8 - 3.2d0*t24*t42 - 128*t34*t9 + t36*t43)
        S_E = pi*(133933333.333333d0*pi*t17/(28.6998365648478d0*t15 + &
        286.998365648478d0) + 4266.66666666667d0*pi*t21*t22 + &
        8533.33333333333d0*pi*t32 - 46666.932416395d0*t15*t50*t53 + t18*( &
        2*t16*t56 + 2800.0159449837d0*t17 + 28000.159449837d0) - 4*t22*( &
        5600.0318899674d0*t21 + t24*t57 + t36*(-t17*t59 + t33*t58 + t43 &
        )) + 1706.66666666667d0*t26*t38 + 2*t26*(5600.0318899674d0*t13 - &
        t27*t56 + t36*(t15*t58 - t17*t26 + t40*t59)) - t44*t45*t48*t52* &
        t54 - 10240.0d0*t44*t47 - 13653.3333333333d0*t44*t49 + &
        31111.2882775967d0*t46*t51 - t47*t55 - 62222.5765551933d0*t48*t51 &
        - t49*t55 - 8*t9*(-t13*t57 + t36*(-t17*t58 + t33*t59 + t37) + &
        5600.0318899674d0*t8))
      S = [ S_rho, S_mx, S_my, S_mz, S_E ]
  end subroutine Source_periodic3d

  !> @brief Analytic source of the field wall2d (mms_gen_fields.py, sympy cse).
  pure subroutine Source_wall2d ( x, y, z, S )
    use MOSE_Base_Types_m, only: R8
    implicit none
    real(R8), intent(in)  :: x, y, z
    real(R8), intent(out) :: S(5)
    real(R8), parameter   :: pi = 3.14159265358979323846d0
    real(R8) :: S_rho, S_mx, S_my, S_mz, S_E
    real(R8) :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9
    real(R8) :: t10, t11, t12, t13, t14, t15, t16, t17, t18, t19
    real(R8) :: t20, t21, t22, t23, t24, t25, t26, t27, t28, t29
    real(R8) :: t30, t31, t32, t33, t34, t35, t36, t37, t38, t39
    real(R8) :: t40, t41, t42, t43, t44, t45, t46, t47, t48, t49
    real(R8) :: t50, t51, t52

        t0 = pi*y
        t1 = sin(t0)
        t2 = 2*pi
        t3 = t2*x
        t4 = cos(t3)
        t5 = 1d0/(1607.19084763148d0*t1 + 401.797711907869d0*t4 + &
        8035.95423815738d0)
        t6 = 6400*t5
        t7 = t1**2
        t8 = t4**2
        t9 = t7*t8
        t10 = sin(t3)
        t11 = t10 + 5
        t12 = t1*t11
        t13 = cos(t0)
        t14 = t10*t13
        t15 = t14*t5
        t16 = t1*t4
        t17 = t13*t4
        t18 = t17 + 10
        t19 = t18*t5
        t20 = t17*t18
        t21 = (0.2d0*t1 + 0.05d0*t4 + 1)**(-2)
        t22 = t20*t21
        t23 = t18*t21
        t24 = t1*t10
        t25 = t11*t24
        t26 = 426.666666666667d0*t1
        t27 = pi*t26
        t28 = 53.3333333333333d0*pi
        t29 = 80.0d0*t12
        t30 = 51200*t5
        t31 = t1**3*t30
        t32 = t11**2
        t33 = t32*t7
        t34 = t4*t7
        t35 = 204800*t34
        t36 = t10*t23
        t37 = 0.637136530181885d0*t36
        t38 = 102400*t5
        t39 = t11*t7
        t40 = 102400*t19
        t41 = t1*t13
        t42 = t10*t39
        t43 = 78127.7777777778d0*pi
        t44 = pi*t8
        t45 = t1*t44
        t46 = 8*pi
        t47 = t11*t13
        t48 = 8*t12
        t49 = t32 + t8
        t50 = t19*t49
        t51 = t50*t7
        t52 = t49*t7
        S_rho = pi*(-0.159284132545471d0*t1*t22 - 12800*t12*t15 + 12800*t16*t19 &
        + t20*t6 + 0.0796420662727356d0*t23*t25 - t6*t9)
        S_mx = pi*(pi*t29 + t10*t27 + t11*t19*t35 - t11*t31*t8 + t12*t20*t38 + &
        t14*t28 - 1600*t14 - 102400*t15*t33 - 1.27427306036377d0*t22*t39 &
        + t33*t37)
        S_my = pi*(t11*t34*t37 - 1.27427306036377d0*t13*t23*t9 - 800*t16 - t17* &
        t28 - t17*t38*t42 + t27*t4 - t31*t4**3 + t40*t41*t8 - t40*t42 + &
        t40*t9)
        S_mz = 0
        S_E = pi*(pi*t48*(160.0d0*t14 + t29) + t1*t43 + t13*t44*( &
        853.333333333333d0*t1 - 853.333333333333d0*t13) + 16*t16*( &
        2800.0159449837d0*t17 + 25600*t51 + 28000.159449837d0) + t17*( &
        22400.1275598696d0*t17 + 204800*t51 + 224001.275598696d0) + t2* &
        t24*(-1280.0d0*t24 + 640.0d0*t47) + t25*t46*(-106.666666666667d0* &
        t13 + t26) - t34*(-409600*t13*t50 + 5.09709224145508d0*t23*t41* &
        t49 + t35*t49*t5 + 22400.1275598696d0*t4) + t4*t43 - 16*t45*( &
        213.333333333333d0*t1 - 53.3333333333333d0*t13) + t45*( &
        853.333333333333d0*t1 + 853.333333333333d0*t13) + t45*(2560.0d0* &
        t1 - 1280.0d0*t13) - t46*t47*(-160.0d0*t24 + 80.0d0*t47) + t48*( &
        -t14*t30*t52 - 5600.0318899674d0*t14 + 512000*t19*t34 + &
        0.318568265090942d0*t36*t52))
      S = [ S_rho, S_mx, S_my, S_mz, S_E ]
  end subroutine Source_wall2d

  !> @brief Analytic source of the field wall3d (mms_gen_fields.py, sympy cse).
  pure subroutine Source_wall3d ( x, y, z, S )
    use MOSE_Base_Types_m, only: R8
    implicit none
    real(R8), intent(in)  :: x, y, z
    real(R8), intent(out) :: S(5)
    real(R8), parameter   :: pi = 3.14159265358979323846d0
    real(R8) :: S_rho, S_mx, S_my, S_mz, S_E
    real(R8) :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9
    real(R8) :: t10, t11, t12, t13, t14, t15, t16, t17, t18, t19
    real(R8) :: t20, t21, t22, t23, t24, t25, t26, t27, t28, t29
    real(R8) :: t30, t31, t32, t33, t34, t35, t36, t37, t38, t39
    real(R8) :: t40, t41, t42, t43, t44, t45, t46, t47, t48, t49
    real(R8) :: t50, t51, t52, t53, t54, t55, t56, t57, t58, t59
    real(R8) :: t60, t61, t62, t63, t64, t65, t66, t67, t68, t69
    real(R8) :: t70, t71, t72, t73, t74, t75, t76, t77, t78, t79
    real(R8) :: t80, t81, t82, t83, t84, t85, t86, t87, t88, t89
    real(R8) :: t90

        t0 = 2*pi
        t1 = t0*x
        t2 = cos(t1)
        t3 = t2**2
        t4 = t0*z
        t5 = cos(t4)
        t6 = t5**2
        t7 = t3*t6
        t8 = pi*y
        t9 = sin(t8)
        t10 = t9**2
        t11 = t2*t5
        t12 = 1d0/(401.797711907869d0*t11 + 1607.19084763148d0*t9 + &
        8035.95423815738d0)
        t13 = t10*t12
        t14 = 6400*t13
        t15 = t11*t9
        t16 = cos(t8)
        t17 = t11*t16
        t18 = t17 + 10
        t19 = t12*t18
        t20 = t17*t19
        t21 = sin(t1)
        t22 = t21*t5
        t23 = t22 + 5
        t24 = t16*t23
        t25 = t12*t9
        t26 = 12800*t22
        t27 = sin(t4)
        t28 = t2*t27
        t29 = 2*t28 + 5
        t30 = t16*t29
        t31 = t28*t30
        t32 = (0.05d0*t11 + 0.2d0*t9 + 1)**(-2)
        t33 = t32*t9
        t34 = t18*t33
        t35 = t22*t23
        t36 = t28*t29
        t37 = 1600*t16
        t38 = 853.333333333333d0*t9
        t39 = pi*t22
        t40 = 53.3333333333333d0*t16
        t41 = pi*t9
        t42 = 80.0d0*t23
        t43 = t12*t9**3
        t44 = 51200*t43
        t45 = 102400*t16
        t46 = t23**2
        t47 = t22*t46
        t48 = t10*t11
        t49 = t21*t27
        t50 = t29*t49
        t51 = t13*t18
        t52 = 51200*t51
        t53 = 51200*t13
        t54 = t10*t32
        t55 = t18*t54
        t56 = t20*t9
        t57 = t10*t23
        t58 = t18*t32
        t59 = t57*t58
        t60 = t10*t29
        t61 = t58*t60
        t62 = 0.318568265090942d0*t61
        t63 = t12*t45
        t64 = t2*t21*t6
        t65 = t16*t7
        t66 = t27*t3*t5
        t67 = t30*t53
        t68 = 102400*t51
        t69 = pi*t28
        t70 = 40.0d0*t29
        t71 = t29**2
        t72 = t28*t71
        t73 = t13*t16
        t74 = t11*t29
        t75 = t41*t7
        t76 = t0*t9
        t77 = -106.666666666667d0*t16 + 213.333333333333d0*t9
        t78 = 8*t23
        t79 = t41*t78
        t80 = 4*t29
        t81 = t41*t80
        t82 = 1280.0d0*t9
        t83 = 160.0d0*t16
        t84 = 160.0d0*t9
        t85 = 4*t46 + 4*t7 + t71
        t86 = t18*t85
        t87 = 5600.0318899674d0*t16
        t88 = t73*t85
        t89 = 0.0796420662727356d0*t54*t86
        t90 = 2*t23
        S_rho = pi*(-t14*t7 + 25600*t15*t19 - 0.159284132545471d0*t17*t34 + 6400 &
        *t20 - t24*t25*t26 - 6400*t25*t31 + 0.0796420662727356d0*t34*t35 &
        + 0.0398210331363678d0*t34*t36)
        S_mx = pi*(-t13*t45*t47 - 1.27427306036377d0*t17*t59 + 307200*t19*t23* &
        t48 - t22*t37 + t23*t28*t62 - t23*t31*t53 - t23*t44*t7 + 102400* &
        t23*t56 + t38*t39 + t39*t40 + t41*t42 + 0.637136530181885d0*t47* &
        t55 - t50*t52)
        S_my = pi*(-106.666666666667d0*pi*t17 + 746.666666666667d0*t11*t41 - &
        102400*t13*t24*t64 - 800*t15 + t18*t63*t7*t9 - t2**3*t44*t5**3 - &
        t35*t68 - t36*t52 + 204800*t51*t7 - 1.27427306036377d0*t55*t65 + &
        0.637136530181885d0*t59*t64 + t62*t66 - t66*t67)
        S_mz = pi*(-0.637136530181885d0*t17*t61 - t23*t49*t68 - t28*t37 - 25600* &
        t29*t43*t7 + 51200*t29*t56 + t35*t62 - t35*t67 + t38*t69 + t40* &
        t69 + t41*t70 + 153600*t51*t74 + 0.159284132545471d0*t55*t72 - &
        25600*t72*t73)
        S_E = pi*(-10240.0d0*pi*t10*t21**2*t27**2 + 156255.555555556d0*pi*t11 - &
        8*pi*t24*(t16*t42 - t22*t84) - 4*pi*t30*(t16*t70 - t28*t84) + pi* &
        t65*(-853.333333333333d0*t16 + 1706.66666666667d0*t9) + 32*t15*( &
        t14*t86 + 2800.0159449837d0*t17 + 28000.159449837d0) + t17*( &
        22400.1275598696d0*t17 + t52*t85 + 224001.275598696d0) + t22*t76* &
        (-t22*t82 + 640.0d0*t24) + t22*t77*t79 + t28*t76*(-t28*t82 + &
        320.0d0*t30) + t28*t77*t81 + 5120.0d0*t39*t57 + &
        78127.7777777778d0*t41 - t48*(22400.1275598696d0*t11 + 51200*t12* &
        t48*t85 + 1.27427306036377d0*t16*t33*t86 - t63*t86) + 2560.0d0* &
        t60*t69 + t7*t76*(-1280.0d0*t16 + 2560.0d0*t9) - 32*t75*( &
        -53.3333333333333d0*t16 + 106.666666666667d0*t9) + t75*( &
        1706.66666666667d0*t16 + t38) - t78*t9*(t22*t87 - t22*t89 + t26* &
        t88 + t52*(-t11*t90 + t50 + 2*t64)) + t79*(t22*t83 + t42*t9) - &
        t80*t9*(t28*t87 + 12800*t28*t88 - t28*t89 + t52*(t49*t90 + 2*t66 &
        - t74)) + t81*(t28*t83 + t70*t9))
      S = [ S_rho, S_mx, S_my, S_mz, S_E ]
  end subroutine Source_wall3d

end program MOSE_mms_program
