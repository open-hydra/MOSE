module MOSE_Lib_Chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Chemistry_Newstate
  public :: Chemistry_Equilibrium
  public :: Setup_Chemistry

contains

  subroutine Setup_Chemistry()
    use MOSE_Config_Types_m,      only: obj_chemistry
    use MOSE_Global_m,            only: nsc
    use FLINT_Lib_Chemistry_wdot, only: Assign_Mechanism
    use FLINT_CEA_setup,          only: CEA_initialize_global
    use oslo,                     only: Setup_ODEsolver
    use strings,                  only: parse
    implicit none
    integer :: ios, i
    character(len=4) :: args(100)
    character(len=128) :: nochem_str
    integer :: nochem_blocks


    ! Initialize chemistry model
    if (obj_chemistry%model=='equilibrium') then

      ! Master thread initializes; other threads will lazy-init
      ! inside CEA_solve on first call (threadprivate variables).
      call CEA_initialize_global()
      obj_chemistry%use_strang = .false.
      obj_chemistry%imodel = 2

    elseif (obj_chemistry%model=='finite-rate') then
   
      call Assign_Mechanism(obj_chemistry%mechanism_name)
      obj_chemistry%iopt = 0
      obj_chemistry%iopt(1) = obj_chemistry%max_ode_steps
      call Setup_ODEsolver(N=(nsc+1),solver=obj_chemistry%ode_name,RT=obj_chemistry%RT,AT=obj_chemistry%AT,iopt=obj_chemistry%iopt)
      obj_chemistry%use_strang = .true.
      obj_chemistry%imodel = 1

    else

      obj_chemistry%use_strang = .false.
      obj_chemistry%imodel = 0
      return

    endif

    ! Parse exclude-blocks string
    nochem_str = trim(adjustl(obj_chemistry%exclude_blocks_str))

    if (nochem_str /= 'none') then
      call parse(nochem_str, ' ', args)
      nochem_blocks = count(args /= '')
      if (nochem_blocks > 0) then
        deallocate(obj_chemistry%no_chem_list)
        allocate(obj_chemistry%no_chem_list(nochem_blocks))
        do i = 1, nochem_blocks
          read(args(i), '(I3)', iostat=ios) obj_chemistry%no_chem_list(i)
          if (ios /= 0) then
            obj_chemistry%warning_message = '[WARNING] No-chemistry blocks parsing failed. String: '//trim(nochem_str)
            deallocate(obj_chemistry%no_chem_list)
            allocate(obj_chemistry%no_chem_list(0))
          end if
        end do
      end if
    else
      allocate(obj_chemistry%no_chem_list(0))
    end if

  end subroutine Setup_Chemistry


  subroutine Chemistry_Newstate ( P, dt, n )
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
    ! Local 
    integer :: lz
    real(R8) :: Z(nsc+1), yeq(nsc)
    real(R8) :: dtchem, t0, rho, temp, teq
    real(R8) :: e_old, T_new, e_new, cv_new, departure
    integer :: i, j, k, s, err, iter

    lz = nsc + 1 ! Z vector length: number of species + temperature.

    !$omp do collapse (3) schedule (dynamic) private(i,j,k,s,Z,dtchem,t0,err,rho)
    do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)

      rho = sum( P(1:nsc,i,j,k) )

      Z(1:nsc) = P(1:nsc,i,j,k)
      Z(lz)    = P(np,i,j,k)/(rho*f_Rtot(Z(1:nsc)))
      dtchem   = dt(i,j,k)
      t0 = 0d0

      ! ODE integration according to selected method
      call Run_ODESolver ( lz, t0, dtchem, Z, rhs_native, err )

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

    !$omp do collapse (3) schedule (dynamic) &
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


end module MOSE_Lib_Chemistry