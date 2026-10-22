module MOSE_Wrap_Postprocess
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: MOSE_postprocess

  character(81), private :: RANS_shell_format  =    "('MOSE | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(69), private :: URANS_shell_format =    "('MOSE | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"
  character(96), private :: RANS_shell_format_MG  = "('MOSE Grid Level', i2, ' | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(84), private :: URANS_shell_format_MG = "('MOSE Grid Level', i2, ' | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"

contains

  subroutine MOSE_postprocess ( simulation )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m
    use MOSE_Global_m
    use MOSE_Parameters_m
    use MOSE_IO_Wall,        only: IOwall, Write_Wall_Solution
    use MOSE_IO_Solution,    only: Write_Solution
    use MOSE_Mod_Diagnostic, only: Write_Diagnostic
    use MOSE_IO_Probes,      only: Write_Probes_Data, Close_Probes
    use MOSE_Lib_RotatingFrame, only: RF_Torque_Walls, RF_Close_Torque_File
    use MOSE_Mod_Multigrid,  only: Prolongation
    use MOSE_Read_Ini,       only: Read_Inifile_Runtime
    use MOSE_Mod_MPI,        only: mpi_is_root
    use MOSE_Mod_Timers,     only: timer_summary
    use MOSE_Lib_Newstate,   only: Report_Negative_Turbulence
    use IR_precision
    implicit none
    type(MOSE_simulation_type), intent(inout) :: simulation
    ! Local
    character(len=llen) :: solfile, dgsfile, wallfile, mgsol, mgwall
    character(len=llen) :: numsol, numwall, numdgs
    integer  :: m, level
    logical  :: iter_based

    level = obj_multigrid%MG_level

    ! Update IOfield solutiontime, on every grid level: a coarse level writes its own solutions
    ! (field-level<L>) with its IOfield
    do m = 1, obj_multigrid%MGL
      if ( obj_time_scheme%time_accurate ) then
        simulation%IOfield(m)%solutiontime = simulation%domain(m)%time
      else
        simulation%IOfield(m)%solutiontime = -real(obj_sim_param%iter_general,8)
      endif
    enddo
    if ( obj_time_scheme%time_accurate ) then
      IOwall%solutiontime = simulation%domain(1)%time
    else
      IOwall%solutiontime = -real(obj_sim_param%iter_general,8)
    endif

    if (level == 1) then
      solfile  = '/'//trim(MOSE_phase_prefix)//'field'
      wallfile = '/'//trim(MOSE_phase_prefix)//'wall'
      dgsfile  = '/'//trim(MOSE_phase_prefix)//'diagnostic'
    else
      solfile  = '/'//trim(MOSE_phase_prefix)//'field-level'//trim(str(.true.,level))
      wallfile = '/'//trim(MOSE_phase_prefix)//'wall-level'//trim(str(.true.,level))
      dgsfile  = '/'//trim(MOSE_phase_prefix)//'diagnostic-level'//trim(str(.true.,level))
      mgsol    = '/'//trim(MOSE_phase_prefix)//'field-prolongated'
      mgwall   = '/'//trim(MOSE_phase_prefix)//'wall-prolongated'
    endif

    ! END OF SIMULATION
    if ( obj_sim_param%TODO == 3 ) then

      if (mpi_is_root) then
        ! Residual history
        if ( mod (simulation%domain(1) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(1)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
        endif
      end if

      ! Probes: the rows due at the last iteration, as at any other
      if (level == 1) call Write_Probes_Data( obj_sim_param%iter_general, simulation%domain(1)%time )

      if (obj_rot%enabled .and. level == 1 .and. mod(simulation%domain(1)%iter, obj_io%res_diter) == 0d0) &
        call RF_Torque_Walls( simulation%domain(1), obj_sim_param%iter_general, simulation%domain(1)%time )

      if (mpi_is_root) then
        ! Shell
        if ( obj_time_scheme%time_accurate ) then
          if (obj_multigrid%MGL > 1 ) then
            write (*,URANS_shell_format_MG) level, obj_sim_param%iter_from_call, simulation%domain(level)%time, simulation%domain(level)%dtglobal
          else
            write (*,URANS_shell_format) obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
          endif
        else
          if (obj_multigrid%MGL > 1 ) then
            write (*,RANS_shell_format_MG) level, obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
          else
            write (*,RANS_shell_format) obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
          endif
        endif
      end if

      ! The solution due at the last iteration (sol-diter, else sol-dtime), as at any other output iteration (TODO == 2):
      ! it takes the next output number and, when due by time only, sets the time of the last output, which the final
      ! solution below then holds for a restart. With sol-overwrite = true the final solution is that output, written
      ! once. Collective: every rank evaluates the same condition. A run ends on the finest level only (Explicit_Step
      ! sets TODO = 3 at level 1), whose domain and file names the block takes: the guard keeps it there
      iter_based = ( mod (simulation%domain(1) % iter, obj_io%sol_diter) == 0 )
      if ( level == 1 .and. ( iter_based .or. &
           simulation%domain(1) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime ) ) then
        obj_io%output_number = obj_io%output_number + 1
        if ( .not. iter_based ) obj_sim_param%time_from_call = simulation%domain(1) % time
        if ( .not. obj_io%sol_overwrite ) then
          if (mpi_is_root) then
            if ( iter_based ) then
              write(*,*) ' ... writing iter-based solution'
            else
              write(*,*) ' ... writing time-based solution'
            endif
          endif
          numsol  = trim(solfile)//trim(str(.true.,obj_io%output_number))
          numwall = trim(wallfile)//trim(str(.true.,obj_io%output_number))
          numdgs  = trim(dgsfile)//trim(str(.true.,obj_io%output_number))
          call Write_Solution ( simulation%domain(1), simulation%IOfield(1), numsol )
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), numwall )
          if (.not. obj_time_scheme%time_accurate) call Write_Diagnostic ( simulation%domain(1), simulation%IOfield(1), numdgs )
        endif
      endif

      ! Timing summary (collective — every rank contributes its accumulators)
      call timer_summary()

      ! Negative turbulence variables after a stage, over the run (collective)
      call Report_Negative_Turbulence()

      ! Write output solution (collective — all ranks participate in gather)
      call Write_Solution ( simulation%domain(1), simulation%IOfield(1), solfile )
      if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), wallfile )
      if (.not. obj_time_scheme%time_accurate) call Write_Diagnostic ( simulation%domain(1), simulation%IOfield(1), dgsfile )

      ! The history files hold their last rows: close them before the ranks finalize MPI (a unit left
      ! open is closed by the Fortran run-time library only when the process ends, after MPI_FINALIZE)
      call Close_Probes()
      if (mpi_is_root) close(obj_io%unitRES)
      call RF_Close_Torque_File()
    

    ! INTERMEDIATE SOLUTION EVALUATION
    elseif ( obj_sim_param%TODO == 2 ) then

      if (level == 1) then
        ! The output number is stored with every vtk solution: a restart goes on with it
        obj_io%output_number = obj_io%output_number + 1
        if (.not.obj_io%sol_overwrite) solfile  = trim(solfile)//trim(str(.true.,obj_io%output_number))
        if (.not.obj_io%sol_overwrite) wallfile = trim(wallfile)//trim(str(.true.,obj_io%output_number))
        if (.not.obj_io%sol_overwrite) dgsfile  = trim(dgsfile)//trim(str(.true.,obj_io%output_number))
      endif

      ! Residuals history
      if (mpi_is_root) then
        if ( mod(simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
        endif
      end if

      ! Probes
      if (level == 1) call Write_Probes_Data( obj_sim_param%iter_general, simulation%domain(1)%time)

      ! Rotating frame torque
      if (level == 1 .and. obj_rot%enabled .and. mod(simulation%domain(1)%iter, obj_io%res_diter) == 0d0) &
        call RF_Torque_Walls(simulation%domain(1), obj_sim_param%iter_general, simulation%domain(1)%time)

      ! Shell (root only)
      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%shell_diter) == 0d0 ) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1 ) then
              write (*,URANS_shell_format_MG) level, obj_sim_param%iter_from_call, simulation%domain(level)%time, simulation%domain(level)%dtglobal
            else
              write (*,URANS_shell_format) obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            endif
          else
            if (obj_multigrid%MGL > 1 ) then
              write (*,RANS_shell_format_MG) level, obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            else
              write (*,RANS_shell_format) obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            endif
          endif
        endif
      end if
      
      ! Solution output (collective — all ranks participate in gather)
      ! Iter-based solution
      if ( mod(simulation%domain(level) % iter, obj_io%sol_diter) == 0d0 ) then
        if (mpi_is_root) write(*,*) ' ... writing iter-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), mgwall )
        endif
        if (level == 1) then
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
          if (.not. obj_time_scheme%time_accurate) call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile )
        endif
      ! Time-based solution
      elseif ( simulation%domain(level) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime ) then
        obj_sim_param%time_from_call = simulation%domain(level) % time
        if (mpi_is_root) write(*,*) ' ... writing time-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), mgwall )
        endif
        if (level == 1) then
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
          if (.not. obj_time_scheme%time_accurate) call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile )
        endif
      endif

      if (obj_multigrid%change_MG) then
        if (mpi_is_root) write(*,*) ' ... writing solution of grid-level ', level
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), mgwall )
        endif
      endif

      if (obj_multigrid%change_MG) then
        obj_multigrid%MG_level = obj_multigrid%MG_level - 1
        obj_sim_param%iter_from_call = 0
        obj_sim_param%time_from_call = 0.0
      endif


    ! AUXILIARY OPERATIONS AT EACH ITERATION
    elseif ( obj_sim_param%TODO == 1 ) then

      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
        end if
      end if

      ! Probes
      if (level == 1) call Write_Probes_Data( obj_sim_param%iter_general, simulation%domain(1)%time )

      ! Rotating frame torque
      if (level == 1 .and. obj_rot%enabled .and. mod(simulation%domain(1)%iter, obj_io%res_diter) == 0d0) &
        call RF_Torque_Walls(simulation%domain(1), obj_sim_param%iter_general, simulation%domain(1)%time)

      ! Shell
      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%shell_diter) == 0d0 ) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1 ) then
              write (*,URANS_shell_format_MG) level, obj_sim_param%iter_from_call, simulation%domain(level)%time, simulation%domain(level)%dtglobal
            else
              write (*,URANS_shell_format) obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            endif
          else
            if (obj_multigrid%MGL > 1 ) then
              write (*,RANS_shell_format_MG) level, obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            else
              write (*,RANS_shell_format) obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            endif
          endif
        endif
      end if

    endif

    ! Update input data
    if ( mod (simulation%domain(level) % iter, obj_io%ini_diter) == 0d0 ) then
      call Read_Inifile_Runtime()
    end if
    

    ! PRINT FOR HYDRA-COUPLING PROBLEMS
    if ( obj_sim_param%HYDRA_postprocess ) then
      if (mpi_is_root) write(*,*) ' ... writing HYDRA-coupling solution'
      call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
      if (obj_multigrid%MGL > 1 .and. level > 1) then
        if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
        do m = obj_multigrid%MG_level, 2, -1
          call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
        enddo
        call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), mgwall )
      endif
      if (obj_io % write_wall .and. level == 1) call Write_Wall_Solution ( simulation%domain(1), wallfile )
    endif

  end subroutine MOSE_postprocess

end module MOSE_Wrap_Postprocess