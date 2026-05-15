module MOSE_Wrap_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: MOSE_setup

contains

  subroutine MOSE_setup ( simulation )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m
    use MOSE_Global_m
    use MOSE_Read_Ini,             only: Read_Inifile
    use MOSE_Load_ThermoTransport, only: Load_ThermoTransport
    use MOSE_Load_Chemistry,       only: Load_Chemistry
    use MOSE_Assign_Setup,         only: Assign_Setup
    use MOSE_IO_Solution,          only: Read_IC, Setup_Output_Solution
    use MOSE_Mod_Allocate_Data,    only: Setup_Data_Structure, deallocate_remote_computation_data
    use MOSE_Mod_Multigrid,        only: Setup_Multigrid, Restriction
    use MOSE_IO_BC,                only: Setup_BC
    use MOSE_Mod_Metrics,          only: Setup_Metrics
    use MOSE_IO_Probes,            only: Setup_Probes
    use MOSE_IO_Wall,              only: Initialize_Wall_File
    use MOSE_Lib_Ghost,            only: Fill_Ghost_Cell
    use MOSE_Mod_MPI,              only: mpi_is_root, partition_blocks
    use MOSE_Mod_GhostExchange,    only: build_ghost_schedule, build_local_bc_index
    implicit none
    type(MOSE_simulation_type), intent(inout) :: simulation
    ! Local
    integer :: m, ios

    !! ------------------------------------------------------
    !! ------------------------------------------------------
    ! Print Header
    if (mpi_is_root) call Print_Header ()

    ! Load tables (thermo, transport, chemistry)
    call Load_ThermoTransport ()
    call Load_Chemistry ()
    if (mpi_is_root) call Check_Tables ()

    ! Read input.ini
    call Read_Inifile ()
    if (mpi_is_root) call Check_Input ()

    ! Allocate container objects Domain and IOfield.
    allocate ( simulation%domain  ( obj_multigrid%MGL ) )
    allocate ( simulation%IOfield ( obj_multigrid%MGL ) )
    simulation%domain(1)%itermax = obj_sim_param%iter_threshold

    ! Assign setup
    call Assign_Setup ()

    ! Read file for initial solution.
    call Read_IC ( simulation%IOfield(1) )
    if (mpi_is_root) call Check_IC ()

    ! Allocation of data structures and copy solution from IOfield.
    call Setup_Data_Structure ( simulation%domain(1), simulation%IOfield(1) )
    allocate(obj_sim_param%residuotot(nres))

    ! With multigrid, allocate Grid and IOfield 2,...,MGL
    if ( obj_multigrid%MGL > 1 ) then
      call Setup_Multigrid ( simulation )
    end if

    ! Read Boundary Conditions file
    call Setup_BC ( simulation%domain )
    if (mpi_is_root) call Check_BC ()

    ! MPI: partition blocks and build ghost cell communication schedule
    do m = 1, obj_multigrid%MGL
      call partition_blocks(simulation%domain(m)%nb, &
        [(product(simulation%domain(m)%blk(ios)%dim(1:3)), ios=1, simulation%domain(m)%nb)])
      call build_ghost_schedule(simulation%domain(m))
      call build_local_bc_index(simulation%domain(m))
    end do

    ! Read grid file and setup metrics in each grid level.
    do m = 1, obj_multigrid%MGL
      call Setup_Metrics ( simulation%domain(m) )
    end do

    ! Free heavy arrays on non-local blocks to save memory
    do m = 1, obj_multigrid%MGL
      call deallocate_remote_computation_data(simulation%domain(m))
    end do

    ! Setup probes location, output files, and variables to be printed.
    call Setup_Probes ( simulation%domain(1), obj_sim_param%newrun )
 
    ! Solution setup
    call Setup_Output_Solution ( simulation%IOfield )

    ! Initialize Wall file
    if (mpi_is_root .and. obj_io % write_wall) &
      call Initialize_Wall_File ( simulation%domain(1), simulation%IOfield(1)%tec%extension)

    ! Print simulation onto the logfile/shell.
    if (mpi_is_root) call Print_Shell_Info ()

    ! Residual history file
    if (mpi_is_root) then
      open(newunit=obj_io%unitRES,file='OUTPUT/'//trim(MOSE_phase_prefix)//'residual-history.dat',status='unknown',form='formatted')
      if ( .not. obj_sim_param%newrun ) then
        ios = 0
        do while ( ios == 0 )
          read (obj_io%unitRES, *, iostat=ios)
        enddo
        backspace (obj_io%unitRES)
      endif
    end if

    ! Setting the format for writing residual_history file
    write(obj_io%unitRES_format,'(A11,I0,A7)') '(I8,E20.10,', nres, 'E20.10)'

    ! Initialize run variables depending on scheme
    select case (trim(obj_time_scheme%solver_type))

      case ('euler', 'RK2', 'RK3')

        obj_multigrid%MG_level = obj_multigrid%MGL

        do m = 1, obj_multigrid%MGL
          simulation%domain(m)%iter = 0

          if (obj_sim_param%newrun) then
            if (obj_time_scheme%time_accurate) then
              simulation%domain(m)%time = 0.0d0
            else
              simulation%domain(m)%time = -1.0d0
            endif
            obj_sim_param%iter_general = 0
          else
            if (obj_time_scheme%time_accurate) then
              simulation%domain(m)%time = simulation%IOfield(m)%solutiontime
              obj_sim_param%iter_general = 0
            else
              obj_sim_param%iter_general = int(simulation%IOfield(m)%solutiontime)
              simulation%domain(m)%time = -1.d0
            endif
          endif
          obj_sim_param%time_from_call = simulation%domain(m)%time
          obj_sim_param%iter_from_call = 0
        enddo
        call Fill_Ghost_Cell ( simulation%domain(1) )

        ! Interpolate initial solution on coarser domains
        if ( obj_multigrid%MGL > 1 ) then
          do m = 2, obj_multigrid%MGL
            call Restriction ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
            call Fill_Ghost_Cell ( simulation%domain(m) )
          end do
        endif

    end select

    ! Print warnings onto the logfile/shell.
    if (mpi_is_root) call Print_Warnings ()

    ! If errors are found, print error messages and stop the simulation.
    if (mpi_is_root) call Stop_Simulation()

    ! Calculate time at beginning of simulation
    call Cpu_Time ( obj_sim_param%cputime(1) )

  contains

    subroutine Print_Header()

      write(*,*)
      write(*,'(A89)')'_/\\\\____________/\\\\________/\\\\\___________/\\\\\\\\\\\_____/\\\\\\\\\\\\\\\__      '
      write(*,'(A89)')'_\/\\\\\\________/\\\\\\______/\\\///\\\_______/\\\/////////\\\__\/\\\///////////___     '
      write(*,'(A89)')' _\/\\\//\\\____/\\\//\\\____/\\\/__\///\\\____\//\\\______\///___\/\\\______________    '
      write(*,'(A89)')'  _\/\\\\///\\\/\\\/_\/\\\___/\\\______\//\\\____\////\\\__________\/\\\\\\\\\\\______   '
      write(*,'(A89)')'   _\/\\\__\///\\\/___\/\\\__\/\\\_______\/\\\_______\////\\\_______\/\\\///////_______  '
      write(*,'(A89)')'    _\/\\\____\///_____\/\\\__\//\\\______/\\\___________\////\\\____\/\\\______________ '
      write(*,'(A89)')'     _\/\\\_____________\/\\\___\///\\\__/\\\______/\\\______\//\\\___\/\\\______________'
      write(*,'(A89)')'      _\/\\\_____________\/\\\_____\///\\\\\/______\///\\\\\\\\\\\/____\/\\\\\\\\\\\\\\\_'
      write(*,'(A89)')'       _\///______________\///________\/////__________\///////////______\///////////////_'
      write(*,*)
      write(*,*)
      write(*,'(A89)') "  MO•dern  /'mɒd.ən/  n. 1. Relating to the present or recent times as opposed to the    "
      write(*,'(A89)') '  remote past.                                                                           '
      write(*,'(A89)') '  2. Characterized by or using the most up-to-date techniques, ideas, or equipment.      '
      write(*,'(A89)') '  3. Denoting the form of language that is currently used as opposed to any earlier form.'
      write(*,'(A89)') "  SE•quel  /'si:kw(ə)l/  n. 1. A published, broadcast, or recorded work that continues   "
      write(*,'(A89)') '  the story or develops the theme of an earlier one.                                     '
      write(*,'(A89)') '  2. Something that takes place after or as a result of an earlier event.                '
      write(*,*)
      write(*,'(A89)') '                                                               Oxford English Dictionary '
      write(*,*)
                                          
    end subroutine Print_Header


    subroutine Check_Tables()
      implicit none
      logical :: has_error

      has_error = .false.

      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Loading'
      write(*,'(A)') ' ========================================================================================='

      ! Thermodynamics
      if (index(obj_thermo%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Thermodynamics', 'FAIL'
        write(*,'(4X,A)') trim(obj_thermo%error_message)
        has_error = .true.
      else
        write(*,'(A,T35,A)') '   Thermodynamics', 'OK'
      endif

      ! Transport
      if (index(obj_transport%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Transport', 'FAIL'
        write(*,'(4X,A)') trim(obj_transport%error_message)
        has_error = .true.
      else
        write(*,'(A,T35,A)') '   Transport', 'OK'
      endif

      ! Chemistry
      if (index(obj_chemistry%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Chemistry', 'FAIL'
        write(*,'(4X,A)') trim(obj_chemistry%error_message)
        has_error = .true.
      else
        write(*,'(A,T35,A)') '   Chemistry', 'OK'
      endif

      if (has_error) stop

    end subroutine Check_Tables


    subroutine Check_IC()
      implicit none

      if (index(obj_io%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Initial conditions', 'FAIL'
        write(*,'(4X,A)') trim(obj_io%error_message)
        stop
      else
        write(*,'(A,T35,A)') '   Initial conditions', 'OK'
      endif

    end subroutine Check_IC


    subroutine Check_Input()
      use MOSE_Input_Registry
      implicit none
      character(len=hlen) :: out

      out = Validate_Registry()
      if (index(out,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Input file', 'FAIL'
        write(*,'(4X,A)') trim(out)
        stop
      else
        write(*,'(A,T35,A)') '   Input file', 'OK'
      endif

    end subroutine Check_Input


    subroutine Check_BC()
      implicit none

      if (index(obj_io_bc%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Boundary conditions', 'FAIL'
        write(*,'(4X,A)') trim(obj_io_bc%error_message)
        stop
      else
        write(*,'(A,T35,A)') '   Boundary conditions', 'OK'
      endif
      write(*,'(A)') ' ========================================================================================='

    end subroutine Check_BC


    subroutine Print_Warnings
      implicit none

      ! IO warnings
      if (index(obj_io%warning_message,'WARNING')>0) write(*,'(A)') obj_io%warning_message
      ! Physics warnings
      if (index(obj_thermo%warning_message,'WARNING')>0) write(*,'(A)') obj_thermo%warning_message
      if (index(obj_transport%warning_message,'WARNING')>0) write(*,'(A)') obj_transport%warning_message
      if (index(obj_rans%warning_message,'WARNING')>0) write(*,'(A)') obj_rans%warning_message
      if (index(obj_rot%warning_message,'WARNING')>0) write(*,'(A)') obj_rot%warning_message
      if (index(obj_chemistry%warning_message,'WARNING')>0) write(*,'(A)') obj_chemistry%warning_message
      if (index(obj_soot%warning_message,'WARNING')>0) write(*,'(A)') obj_soot%warning_message
      ! Numerical scheme warnings
      if (index(obj_time_scheme%warning_message,'WARNING')>0) write(*,'(A)') obj_time_scheme%warning_message
      if (index(obj_riemann%warning_message,'WARNING')>0) write(*,'(A)') obj_riemann%warning_message
      if (index(obj_space_scheme%warning_message,'WARNING')>0) write(*,'(A)') obj_space_scheme%warning_message
      if (index(obj_irs%warning_message,'WARNING')>0) write(*,'(A)') obj_irs%warning_message

    end subroutine Print_Warnings


    subroutine Stop_Simulation()
      implicit none
      logical :: has_error

      has_error = .false.

      ! Physics errors

      if (index(obj_rans%error_message,'ERROR')>0) then
        write(*,'(A)') obj_rans%error_message;       has_error = .true.
      endif
      if (index(obj_rot%error_message,'ERROR')>0) then
        write(*,'(A)') obj_rot%error_message;       has_error = .true.
      endif
      if (index(obj_soot%error_message,'ERROR')>0) then
        write(*,'(A)') obj_soot%error_message;       has_error = .true.
      endif
      ! Numerical scheme errors
      if (index(obj_time_scheme%error_message,'ERROR')>0) then
        write(*,'(A)') obj_time_scheme%error_message;  has_error = .true.
      endif
      if (index(obj_riemann%error_message,'ERROR')>0) then
        write(*,'(A)') obj_riemann%error_message;      has_error = .true.
      endif
      if (index(obj_space_scheme%error_message,'ERROR')>0) then
        write(*,'(A)') obj_space_scheme%error_message; has_error = .true.
      endif
      if (index(obj_irs%error_message,'ERROR')>0) then
        write(*,'(A)') obj_irs%error_message;          has_error = .true.
      endif

      if (has_error) stop

    end subroutine Stop_Simulation


    subroutine Print_Shell_Info()
      use IR_Precision,   only: str  
      use MOSE_IO_BC,     only: Print_BC_Summary
      use MOSE_IO_Probes, only: nprobes
      implicit none
      ! Local
      character(llen) :: eosword
      integer :: b, total_cells, k

      ! ----- Domain topology -----
      total_cells = 0
      do b = 1, simulation%domain(1)%nb
        total_cells = total_cells + &
          product(simulation%domain(1)%blk(b)%dim(1:3))
      end do

      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Set-up'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Domain'
      write(*,'(A,T35,I0)') '   Blocks', simulation%domain(1)%nb
      write(*,'(A,T35,I0)') '   Cells', total_cells
      write(*,'(A,T35,I0)') '   Boundary faces', simulation%domain(1)%nbound
      if ( obj_multigrid%MGL > 1 ) &
        write(*,'(A,T35,I0)') '   Multigrid levels', obj_multigrid%MGL

      ! ----- Boundary conditions -----
      call Print_BC_Summary ()

      ! ----- IO -----
      write(*,*)
      write(*,'(A)') ' Input/Output'
      write(*,'(A,T35,A)') '   Initial conditions file', trim(obj_io%nameinit)
      write(*,'(A,T35,A)') '   Solution format', trim(obj_io%sol_format)
      write(*,'(A,T35,A)') '   Probes number', str(.true.,nprobes)
      write(*,'(A)') ' ========================================================================================='

      ! ----- Physical model -----
      eosword = 'Ideal'

      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Physical model'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Gas model'
      write(*,'(A,T35,A)') '   Equations', trim(obj_sim_param%description)
      if (model==1) then
        write(*,'(A,T35,A)') '   Viscous model', 'Laminar'
      elseif (model==2) then
        write(*,'(A,T35,A)') '   Viscous model', 'Turbulent, '//trim(obj_rans%model)
      end if
      write(*,'(A,T35,A)') '   Equation of state', trim(eosword)
      write(*,'(A,T35,A)') '   Thermodynamics', trim(obj_thermo%description)
      write(*,'(A,T35,I0)') '   Species', nsc

      if (obj_chemistry%model/='frozen') then
        write(*,*)
        write(*,'(A)') ' Chemistry'
        if (obj_chemistry%model=='finite-rate') then
          write(*,'(A)')       '   Finite-rate '
          write(*,'(A,T35,A)') '   - Mechanism', trim(obj_chemistry%mechanism_name)
          write(*,'(A,T35,A)') '   - Splitting', 'Strang operator'
          write(*,'(A,T35,A)') '   - ODE solver', trim(obj_chemistry%ode_name)
        elseif (obj_chemistry%model=='equilibrium') then
          write(*,'(A)')       '   CEA equilibrium '
        endif
        if ( size(obj_chemistry%no_chem_list) == 0 ) then
          write(*,'(A,T35,A)') '   Active blocks', 'all'
        else
          write(*,'(A,T35,*(I0,:," "))') '   Inactive blocks', &
            obj_chemistry%no_chem_list
        endif
      endif

      if (obj_soot%enabled) then
        write(*,*)
        write(*,'(A)') ' Soot'
        if (trim(obj_soot%model)=='LL91') then
          write(*,'(A,T35,A)') '   Mechanism', 'Leung (1991)'
        elseif (trim(obj_soot%model)=='LIN') then
          write(*,'(A,T35,A)') '   Mechanism', 'Lindsted (1994)'
        end if
        if (nsoot==1) write(*,'(A,T35,A)') '   Model', 'C(gr) => Heavy Gas'
      end if

      if (npass>0) write(*,'(A,T35,I0)') '   Passive scalars', npass

      if (obj_rot%enabled) then
        write(*,*)
        write(*,'(A)')         ' Rotating Frame '
        write(*,'(A,F14.4,A)') '   omega  = ', obj_rot%omega,   ' rad/s'
        write(*,'(A,3F11.5)')  '   axis   = ', obj_rot%axis
        write(*,'(A,3F11.5)')  '   origin = ', obj_rot%origin
        if ( allocated(obj_rot%stationary_face) ) then
          write(*,'(A,I0,A)') '   stationary patches = ', size(obj_rot%stationary_face,2), ':'
          do k = 1, size(obj_rot%stationary_face, 2)
            write(*,'(A,I0,A,I0)') '     block ', obj_rot%stationary_face(1,k), &
                                    '  face-dir ', obj_rot%stationary_face(2,k)
          end do
        else
          write(*,'(A)') '   all type-5/6 walls rotate with the frame'
        end if
      end if

      write(*,'(A)') ' ========================================================================================='

      ! ----- Numerical scheme -----
      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Numerical scheme'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A,T35,A)') ' Space'
      write(*,'(A,T35,A)') '   Reconstruction', trim(obj_space_scheme%description)
      if (len_trim(obj_space_scheme%flux_limiter)>0) &
        write(*,'(A,T35,A)') '   Flux limiter',trim(obj_space_scheme%flux_limiter)
      write(*,*)
      write(*,'(A,T35,A)') ' Time'
      write(*,'(A,T35,A)') '   Scheme', trim(obj_time_scheme%description)
      write(*,'(A,T35,A)') '   Integration variables', trim(obj_time_scheme%integration_variables)
      if (len_trim(obj_irs%description)>0) &
        write(*,'(A,T35,A)') '   Implicit residual smoothing', trim(obj_irs%description)
      write(*,*)
      write(*,'(A,T35,A)') ' Fluxes'
      write(*,'(A,T35,A)') '   Riemann solver', trim(obj_riemann%description)
      write(*,'(A)') ' ========================================================================================='
      write(*,*)

    end subroutine Print_Shell_Info

  end subroutine MOSE_setup

end module MOSE_Wrap_Setup