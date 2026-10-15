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
   
      call Assign_Mechanism(obj_chemistry%mechanism_name)
      ! Must follow Assign_Mechanism: it resolves against the mechanism just
      ! selected, and falls back to finite differences if that one has no
      ! analytical Jacobian.
      call set_analytical_jacobian(obj_chemistry%analytical_jacobian)
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

    ! exclude-blocks: none, or the numbers of the blocks without chemistry, separated by blanks. Every word must be a
    ! block number, digits only, read whole; any other word stops the run at set-up with the key, the value and the word
    ! (Stop_Simulation, on every rank). The numbers beyond the mesh stop it once the mesh gives the count of its blocks
    ! (Check_Exclude_Blocks in Wrap_Setup). A word that did not read used to empty the list with a warning, and the read of
    ! the next word wrote beyond the emptied list; a number read with I3 kept its first three digits.
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

    !$omp do collapse (3) schedule (dynamic, 64) private(i,j,k,s,Z,dtchem,t0,err,rho)
    do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)

      rho = sum( P(1:nsc,i,j,k) )

      Z(1:nsc) = P(1:nsc,i,j,k)
      Z(lz)    = P(np,i,j,k)/(rho*f_Rtot(Z(1:nsc)))
      dtchem   = dt(i,j,k)
      t0 = 0d0

      ! ODE integration according to selected method.  IJAC_chem is 0 unless
      ! an analytical Jacobian was both asked for and available.
      call Run_ODESolver ( lz, t0, dtchem, Z, rhs_native, jac_native, IJAC_chem, err )

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


end module MOSE_Lib_Chemistry