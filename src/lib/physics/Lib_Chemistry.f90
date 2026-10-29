module MOSE_Lib_Chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Parameters_m, only: llen

  implicit none
  private
  public :: Chemistry_Newstate
  public :: Chemistry_Equilibrium
  public :: Setup_Chemistry
  public :: Chemistry_Failure_Reduce, Chemistry_Failed, Chemistry_Failure_Stop

  ! The integrator of the finite-rate chemistry that Setup_Chemistry prepares in OSLO, once at set-up: a reading of input.ini
  ! during the run does not prepare another one, so the convention of the return code and the name in the message of a failure
  ! come from this name, not from the value of ode-solver read again
  character(len=llen) :: chem_ode_solver = ''

  ! A failed integration of the chemistry on this rank: among the cells of the step whose integration failed, the one of the
  ! smallest block, k, j, i, so that the record does not depend on the number of threads nor on the order in which the
  ! dynamic schedule gives them the cells; with the state of the cell before the chemistry step and the step
  type :: chemistry_failure_t
    logical  :: found = .false.
    integer  :: b = 0, i = 0, j = 0, k = 0, code = 0
    real(R8) :: T = 0.0_R8, p = 0.0_R8, rho = 0.0_R8, dt = 0.0_R8
  end type chemistry_failure_t
  type(chemistry_failure_t) :: chem_failure

  ! The smallest number of a block with a failed integration on any rank, huge(1.0_R8) when there is none
  ! (Chemistry_Failure_Reduce): every thread of every rank reads the same value
  real(R8) :: chem_failed_block = huge(1.0_R8)

contains

  subroutine Setup_Chemistry()
    use MOSE_Config_Types_m,      only: obj_chemistry
    use MOSE_Global_m,            only: nsc, MOSE_phase_prefix
    use FLINT_Lib_Thermodynamic,  only: species_composition
    use FLINT_Lib_Chemistry_wdot, only: Assign_Mechanism
    use FLINT_Lib_Chemistry_rhs,  only: set_analytical_jacobian
    use FLINT_CEA_setup,          only: CEA_initialize_global
    use oslo,                     only: Setup_ODEsolver
    use strings,                  only: parse
    use MOSE_Parameters_m,        only: llen
    implicit none
    integer :: ios, i
    logical :: elements
    character(len=llen) :: args(llen/2 + 1)   ! the words of exclude-blocks: at most one in two characters of the value
    character(len=llen) :: nochem_str
    integer :: nochem_blocks


    ! Initialize chemistry model
    if (obj_chemistry%model=='equilibrium') then

      ! The equilibrium keeps the amount of each element, so it needs the elements of the species, which
      ! INPUT/<prefix>composition.txt gives; when that file is missing the tables load without them
      ! (the number of elements is the first dimension of species_composition, which the file allocates)
      elements = allocated(species_composition)
      if (elements) elements = size(species_composition, 1) > 0
      if (elements) then
        ! Master thread initializes; other threads will lazy-init
        ! inside CEA_solve on first call (threadprivate variables).
        call CEA_initialize_global()
      else
        obj_chemistry%error_message = '[ERROR] chemistry = equilibrium needs the elements of the species, from INPUT/'// &
                                      trim(MOSE_phase_prefix)//'composition.txt: none was read'
      end if
      obj_chemistry%use_strang = .false.
      obj_chemistry%imodel = 2

    elseif (obj_chemistry%model=='finite-rate') then

      ! The finite-rate chemistry needs the reaction mechanism, INPUT/<prefix>chemistry-info.txt: without it read_chemistry
      ! reads no rate table (Load_Chemistry), and FLINT answers every evaluation of the source terms with its bail-out value
      if (.not. obj_chemistry%info_found) then
        obj_chemistry%error_message = '[ERROR] chemistry = finite-rate needs the reaction mechanism, INPUT/'// &
                                      trim(MOSE_phase_prefix)//'chemistry-info.txt: the file cannot be opened'
      else
        call Assign_Mechanism(obj_chemistry%mechanism_name)
        ! Must follow Assign_Mechanism: it resolves against the mechanism just
        ! selected, and falls back to finite differences if that one has no
        ! analytical Jacobian.
        call set_analytical_jacobian(obj_chemistry%analytical_jacobian)
        obj_chemistry%iopt = 0
        obj_chemistry%iopt(1) = obj_chemistry%max_ode_steps
        call Setup_ODEsolver(N=(nsc+1),solver=obj_chemistry%ode_name,RT=obj_chemistry%RT,AT=obj_chemistry%AT,iopt=obj_chemistry%iopt)
        chem_ode_solver = obj_chemistry%ode_name
      end if
      obj_chemistry%use_strang = .true.
      obj_chemistry%imodel = 1

    else

      obj_chemistry%use_strang = .false.
      obj_chemistry%imodel = 0
      return

    endif

    ! exclude-blocks: none, or the numbers of the blocks without chemistry, separated by blanks. Every word must be a
    ! block number, digits only, read whole; any other word stops the run at set-up with the key, the value and the word
    ! (Stop_Simulation, on every rank). The numbers beyond the mesh stop it once the mesh gives the count of its blocks
    ! (Check_Exclude_Blocks in Wrap_Setup). Each word is checked before it is read: digits only, at most nine, read
    ! list-directed and whole, never cut to its first digits.
    if (allocated(obj_chemistry%no_chem_list)) deallocate(obj_chemistry%no_chem_list)
    nochem_str = trim(adjustl(obj_chemistry%exclude_blocks_str))
    if (nochem_str == 'none') then
      allocate(obj_chemistry%no_chem_list(0))
    else
      call parse(nochem_str, ' ', args)
      nochem_blocks = count(args /= '')
      allocate(obj_chemistry%no_chem_list(nochem_blocks))
      do i = 1, nochem_blocks
        ios = 1
        if (verify(trim(args(i)), '0123456789') == 0 .and. len_trim(args(i)) <= 9) &
          read(args(i), *, iostat=ios) obj_chemistry%no_chem_list(i)
        if (ios == 0) then
          if (obj_chemistry%no_chem_list(i) < 1) ios = 1
        end if
        if (ios /= 0) then
          if (index(obj_chemistry%error_message, 'ERROR') == 0) obj_chemistry%error_message = &
            '[ERROR] input.ini, [MOSE-Chemistry] exclude-blocks = '//trim(nochem_str)//': "'//trim(args(i))// &
            '" is not a block number'
          exit
        end if
      end do
    end if

  end subroutine Setup_Chemistry


  subroutine Chemistry_Newstate ( P, dt, n, block_id )
    use MOSE_Global_m
    use MOSE_Parameters_m
    use FLINT_Lib_Chemistry_rhs
    use FLINT_Lib_Thermodynamic, only: f_Rtot, f_tabT, h_tab, &
                                       cp_tab, Ri_tab
    use FLINT_CEA_solver,        only: CEA_solve
    use oslo,                    only: Run_ODESolver

    implicit none
    integer, intent(in) :: n(3)
    real(R8), dimension(n(1),n(2),n(3)), intent(in) :: dt
    real(R8), dimension( nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc ) :: P
    integer, intent(in) :: block_id   ! the number of the block, for the record of a failed integration
    ! Local
    integer :: lz
    real(R8) :: Z(nsc+1), yeq(nsc)
    real(R8) :: dtchem, t0, rho, temp, teq
    real(R8) :: e_old, T_new, e_new, cv_new, departure
    integer :: i, j, k, s, err, iter

    lz = nsc + 1 ! Z vector length: number of species + temperature.

    !$omp do collapse (3) schedule (dynamic, 64) private(i,j,k,s,Z,dtchem,t0,err,rho,temp)
    do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)

      rho = sum( P(1:nsc,i,j,k) )

      Z(1:nsc) = P(1:nsc,i,j,k)
      Z(lz)    = P(np,i,j,k)/(rho*f_Rtot(Z(1:nsc)))
      temp     = Z(lz)
      dtchem   = dt(i,j,k)
      t0 = 0d0

      ! ODE integration according to selected method.  IJAC_chem is 0 unless
      ! an analytical Jacobian was both asked for and available.
      call Run_ODESolver ( lz, t0, dtchem, Z, rhs_native, jac_native, IJAC_chem, err )

      ! An integration that the integrator reports as failed has not reached the end of the step: Z holds the state where it
      ! stopped, or where it started. The cell keeps its state and the failure is recorded; the run stops after the chemistry
      ! of the step, on every thread and rank (Explicit_Step), before any flux reads the cell
      if ( Integration_Failed(err) ) then
        call Record_Failure ( block_id, i, j, k, err, temp, P(np,i,j,k), rho, dt(i,j,k) )
        cycle
      end if

      ! New state from ODE integration in dtchem
      do s = 1, nsc
        if (Z(s)<=1d-20) Z(s) = 1d-20
        ! mass fractions
        P(s,i,j,k) = Z(s)
      end do
      ! pressure
      P(np,i,j,k) = rho*f_Rtot(Z(1:nsc))*Z(lz)

    enddo; enddo; enddo

  end subroutine Chemistry_Newstate


  subroutine Chemistry_Equilibrium ( P, n )
    use MOSE_Global_m
    use FLINT_CEA_solver, only: CEA_solve
    use FLINT_Lib_Thermodynamic, only: f_Rtot
    implicit none
    integer, intent(in) :: n(3)
    real(R8), dimension( nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc ) :: P
    ! Local
    real(R8) :: rho, temp, teq, yeq(nsc), departure
    integer :: i, j, k, s

    !$omp do collapse (3) schedule (dynamic, 64) &
    !$omp private(i,j,k,s,rho,temp,teq,yeq,departure)
    do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)

      rho  = sum( P(1:nsc,i,j,k) )
      temp = P(np,i,j,k)/(rho*f_Rtot(P(1:nsc,i,j,k)))

      call CEA_solve(temp, P(1:nsc,i,j,k), teq, yeq)

      ! Measure departure from equilibrium (max relative change in Y)
      departure = 0.0_R8
      do s = 1, nsc
        if (yeq(s) > 1d-20) then
          departure = max(departure, abs(yeq(s) - P(s,i,j,k)/rho) / yeq(s))
        end if
      end do
      departure = max(departure, abs(teq - temp)/teq)

      ! Skip correction if state is already at equilibrium to avoid noise accumulation
      if (departure < 1d-2) cycle

      ! Replace composition with equilibrium
      do s = 1, nsc
        if (yeq(s) <= 1d-20) yeq(s) = 1d-20
        P(s,i,j,k) = yeq(s)*rho
      end do

      ! Pressure from flow-solver EOS
      P(np,i,j,k) = rho*f_Rtot(P(1:nsc,i,j,k))*teq

    enddo; enddo; enddo

  end subroutine Chemistry_Equilibrium


  !> Whether a return code of the integrator that Setup_Chemistry prepared (chem_ode_solver) reports a failure, by the
  !> convention of that integrator: dodesol of Intel ODE returns 0 on success and a negative code on an error (Return Values
  !> in ODE_manual.pdf of Intel ODE); the integrators of Hairer, of FATODE and CVODE return a positive code on success (CVODE
  !> also 0) and a negative one on a failure (the headers of radau5.f, rodas.f, sdirk4.f, dopri5.f, SDIRK_f90_Integrator.F90,
  !> RK_f90_Integrator.F90, ROS_f90_Integrator.F90 in OSLO; cvodes.h)
  logical function Integration_Failed ( code )
    implicit none
    integer, intent(in) :: code

    if ( chem_ode_solver == 'dodesol' ) then
      Integration_Failed = code /= 0
    else
      Integration_Failed = code < 0
    end if

  end function Integration_Failed


  !> Records a failed integration (inside the OpenMP loop of Chemistry_Newstate: one thread at a time), keeping, among those
  !> of the step on this rank, the one of the smallest block, k, j, i
  subroutine Record_Failure ( block_id, i, j, k, code, T, p, rho, dt )
    implicit none
    integer,  intent(in) :: block_id, i, j, k, code
    real(R8), intent(in) :: T, p, rho, dt   ! the state of the cell before the chemistry step, and the step

    !$omp critical (MOSE_chemistry_failure)
    if ( .not. chem_failure%found .or. Comes_Before([block_id, k, j, i], &
                                                    [chem_failure%b, chem_failure%k, chem_failure%j, chem_failure%i]) ) &
      chem_failure = chemistry_failure_t(.true., block_id, i, j, k, code, T, p, rho, dt)
    !$omp end critical (MOSE_chemistry_failure)

  end subroutine Record_Failure


  !> a before c in the lexicographic order
  pure logical function Comes_Before ( a, c )
    implicit none
    integer, intent(in) :: a(:), c(:)
    integer :: m

    Comes_Before = .false.
    do m = 1, size(a)
      if ( a(m) /= c(m) ) then
        Comes_Before = a(m) < c(m)
        return
      end if
    end do

  end function Comes_Before


  !> The decision of every rank after the chemistry of a step: the smallest number of a block with a failed integration on any
  !> rank. Called inside the parallel region of the step by the primary thread alone (master construct: MPI is initialised
  !> with MPI_THREAD_FUNNELED), after a barrier, so that every thread has integrated its cells; the other threads read the
  !> result (Chemistry_Failed) after the barrier that follows.
  subroutine Chemistry_Failure_Reduce ()
    use MOSE_Mod_MPI, only: mpi_allreduce_min_r8
    implicit none
    real(R8) :: smallest

    smallest = huge(1.0_R8)
    if ( chem_failure%found ) smallest = real(chem_failure%b, R8)
    call mpi_allreduce_min_r8 ( smallest, chem_failed_block )

  end subroutine Chemistry_Failure_Reduce


  !> An integration of the chemistry failed on some rank (Chemistry_Failure_Reduce): the same answer on every thread and rank
  logical function Chemistry_Failed ()
    implicit none

    Chemistry_Failed = chem_failed_block < huge(1.0_R8)

  end function Chemistry_Failed


  !> After the parallel region of a step whose chemistry failed (Chemistry_Failed), on every rank: the rank that holds the
  !> smallest block with a failed integration writes the failure of its smallest cell, once, so that the line is the same
  !> with any number of threads and ranks and with any decomposition, as Check_IC_Values does (Wrap_Setup); then every rank
  !> stops with mpi_stop_all, exit code 1, and no solution of the step is written.
  subroutine Chemistry_Failure_Stop ( level, levels )
    use MOSE_Mod_MPI, only: mpi_size_, mpi_rank_, mpi_stop_all
    implicit none
    integer, intent(in) :: level, levels   ! the grid level of the step and the number of grid levels of the run
    character(len=64) :: place, rank

    if ( chem_failure%found .and. real(chem_failure%b, R8) == chem_failed_block ) then
      write(place,'(A,I0)') 'block ', chem_failure%b
      if ( levels > 1 ) write(place,'(A,I0,A,I0)') 'block ', chem_failure%b, ' of grid level ', level
      rank = ''
      if ( mpi_size_ > 1 ) write(rank,'(A,I0,A)') ' (rank ', mpi_rank_, ')'
      write(*,'(A)') ' [ERROR] Chemistry: the integrator '//trim(chem_ode_solver)//' failed (return code '// &
        Integer_Text(chem_failure%code)//': '//Code_Meaning(chem_failure%code)//') in '//trim(place)// &
        ', cell i, j, k = '//Integer_Text(chem_failure%i)//', '//Integer_Text(chem_failure%j)//', '// &
        Integer_Text(chem_failure%k)//', over a chemistry step of '//Real_Text(chem_failure%dt)//' s from T = '// &
        Real_Text(chem_failure%T)//' K, p = '//Real_Text(chem_failure%p)//' Pa, rho = '//Real_Text(chem_failure%rho)// &
        ' kg/m3'//trim(rank)
    end if
    call mpi_stop_all ( 'chemistry integration error (see above)' )

  end subroutine Chemistry_Failure_Stop


  !> What a return code of the integrator chem_ode_solver means, from its documentation: the headers of the integrators of
  !> Hairer (radau5.f, rodas.f, dopri5.f, sdirk4.f), the error messages of FATODE (SDIRK_ErrorMsg, RK_ErrorMsg, ros_ErrorMsg),
  !> the return values of CVODES (cvodes.h) and those of Intel ODE (Return Values in ODE_manual.pdf). The bound on the steps
  !> of an integration is ode-max-steps (Setup_Chemistry passes it to the integrators of Hairer and of FATODE); dodesol and
  !> cvode do not take it (OSLO gives cvode a bound of 10000 steps).
  function Code_Meaning ( code ) result(text)
    implicit none
    integer, intent(in) :: code
    character(len=:), allocatable :: text
    character(len=*), parameter :: steps = 'more steps are needed than ode-max-steps'
    character(len=20), parameter :: cvode(28) = [character(len=20) :: &
      'CV_TOO_MUCH_WORK', 'CV_TOO_MUCH_ACC', 'CV_ERR_FAILURE', 'CV_CONV_FAILURE', 'CV_LINIT_FAIL', 'CV_LSETUP_FAIL', &
      'CV_LSOLVE_FAIL', 'CV_RHSFUNC_FAIL', 'CV_FIRST_RHSFUNC_ERR', 'CV_REPTD_RHSFUNC_ERR', 'CV_UNREC_RHSFUNC_ERR', &
      'CV_RTFUNC_FAIL', 'CV_NLS_INIT_FAIL', 'CV_NLS_SETUP_FAIL', 'CV_CONSTR_FAIL', 'CV_NLS_FAIL', '', '', '', &
      'CV_MEM_FAIL', 'CV_MEM_NULL', 'CV_ILL_INPUT', 'CV_NO_MALLOC', 'CV_BAD_K', 'CV_BAD_T', 'CV_BAD_DKY', 'CV_TOO_CLOSE', &
      'CV_VECTOROP_ERR']

    text = 'a code that this message does not describe: see the documentation of the integrator'
    select case ( trim(chem_ode_solver) )
    case ('H-radau5', 'H-rodas')
      select case (code)
      case (-1); text = 'input is not consistent'
      case (-2); text = steps
      case (-3); text = 'step size becomes too small'
      case (-4); text = 'matrix is repeatedly singular'
      end select
    case ('H-dopri5')
      select case (code)
      case (-1); text = 'input is not consistent'
      case (-2); text = steps
      case (-3); text = 'step size becomes too small'
      case (-4); text = 'problem is probably stiff'
      end select
    case ('H-sdirk4')
      if ( code == -1 ) text = 'computation unsuccessful'
    case ('sdirk2a', 'sdirk2b', 'sdirk3a', 'sdirk4a', 'sdirk4b')
      select case (code)
      case (-1); text = 'improper value for the maximal number of steps'
      case (-2); text = 'improper value for the maximal number of Newton iterations'
      case (-3); text = 'Hmin, Hmax or Hstart not positive'
      case (-4); text = 'FacMin, FacMax or FacRej not positive'
      case (-5); text = 'improper tolerance values'
      case (-6); text = steps
      case (-7); text = 'step size too small'
      case (-8); text = 'matrix is repeatedly singular'
      end select
    case ('radau2a', 'lobatto3c', 'gauss', 'radau1a')
      select case (code)
      case (-1); text = 'improper value for the maximal number of steps'
      case (-2); text = 'improper value for the maximal number of Newton iterations'
      case (-3); text = 'Hmin, Hmax or Hstart not positive'
      case (-4); text = 'improper values for FacMin, FacMax, FacSafe or FacRej'
      case (-5); text = 'improper value for ThetaMin'
      case (-6); text = 'Newton stopping tolerance too small'
      case (-7); text = 'improper values for Qmin, Qmax'
      case (-8); text = 'tolerances are too small'
      case (-9); text = steps
      case (-10); text = 'step size too small'
      case (-11); text = 'matrix is repeatedly singular'
      case (-12); text = 'non-convergence of the Newton iterations'
      case (-13); text = 'requested Runge-Kutta method not implemented'
      end select
    case ('ros2', 'ros3', 'ros4', 'rodas3', 'rodas4')
      select case (code)
      case (-1); text = 'improper value for the maximal number of steps'
      case (-2); text = 'selected Rosenbrock method not implemented'
      case (-3); text = 'Hmin, Hmax or Hstart not positive'
      case (-4); text = 'FacMin, FacMax or FacRej not positive'
      case (-5); text = 'improper tolerance values'
      case (-6); text = steps
      case (-7); text = 'step size too small'
      case (-8); text = 'matrix is repeatedly singular'
      end select
    case ('cvode')
      if ( code <= -1 .and. code >= -size(cvode) ) then
        if ( len_trim(cvode(-code)) > 0 ) text = trim(cvode(-code))
      end if
    case ('dodesol')
      select case (code)
      case (-100); text = 'the number of equations n is less than zero'
      case (-201); text = 'the initial value of t is greater than the end of the integration interval'
      case (-202); text = 'the initial step size h or the minimal step size hm is not positive'
      case (-203); text = 'the relative error tolerance ep is not positive'
      case (-204); text = 'the threshold for the control of the relative error is not positive'
      end select
    end select

  end function Code_Meaning


  function Integer_Text ( n ) result(text)
    implicit none
    integer, intent(in) :: n
    character(len=:), allocatable :: text
    character(len=16) :: buffer

    write(buffer,'(I0)') n
    text = trim(buffer)

  end function Integer_Text


  function Real_Text ( x ) result(text)
    implicit none
    real(R8), intent(in) :: x
    character(len=:), allocatable :: text
    character(len=32) :: buffer

    write(buffer,'(ES13.6)') x
    text = trim(adjustl(buffer))

  end function Real_Text


end module MOSE_Lib_Chemistry