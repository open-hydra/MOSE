module MOSE_Mod_Diagnostic
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Parameters_m

  implicit none
  private
  public :: Compute_Residual, Write_Diagnostic

contains

  subroutine Compute_Residual ( new, old, dt, n, average, total )
    use MOSE_Global_m
    implicit none
    integer, intent(in)     :: n(3)
    real(R8), intent(in)    :: new ( nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc )
    real(R8), intent(in)    :: old ( nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc )
    real(R8), intent(in)    :: dt ( n(1), n(2), n(3) )
    real(R8), intent(out)   :: average(nres)
    real(R8), intent(inout) :: total(nres)
    ! Local
    integer :: i, j, k
    real(R8) :: resn(nres), residuo(nres)

    residuo = 0d0

    !$omp parallel
    !$omp do private ( resn ) reduction ( + : residuo ) collapse (3)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      resn(1) = abs ( sum(new(1:nsc,i,j,k) - old(1:nsc,i,j,k)) )  ! delta(rho)
      resn(2:nres) = abs ( new(nu:nprim,i,j,k) - old(nu:nprim,i,j,k) ) ! delta(...)
      residuo = residuo + resn*resn  ! sum ( r(i)**2 )

    enddo; enddo; enddo
    !$omp end parallel

    ! do i = 1, nres; resmax_(i) = maxval ( reslocal(i,:,:,:) ); enddo

    ! Average residual of this block
    average = sqrt ( residuo / float ( n(1)*n(2)*n(3) ) )
    ! Overall residual
    total = total + residuo

  end subroutine Compute_Residual


  !> Update the orion-field data and write them accordingly to the chosen format (vtk,tecplot)
  subroutine Write_Diagnostic( domain, IOfield, file )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_io
    use MOSE_Global_m
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use MOSE_Mod_MPI, only: mpi_is_root, mpi_abort_all
    use MOSE_Mod_GhostExchange, only: gather_diagnostic_to_root, mpi_io_barrier
    use strings, only: parse
    use IR_Precision, only: str
    implicit none
    type(MOSE_domain_type), intent(inout)  :: domain
    type(ORION_data), intent(inout)        :: IOfield
    character(llen), intent(in)            :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    character(len=hlen) :: varnames, vtknames
    integer             :: E_IO, b, i, j, k, ndiag
    character(len=clen) :: format(2), extension
    real(R8), allocatable :: full(:,:,:,:)

    ! Residuals of rho, rho*u, rho*v, rho*w, rho*E and of the variables that follow the energy in r: soot, passive
    ! scalars, RANS (named as in the solution files, Setup_Output_Solution), then dt and beta
    ndiag    = nres + 2
    varnames = '"rho" "rhou" "rhov" "rhow" "rhoe"'
    if ( nsoot == 1 ) varnames = trim(varnames)//' "rho*Np"'
    if ( nsoot == 2 ) varnames = trim(varnames)//' "rho*Yp" "rho*Np"'
    do i = 1, npass
      varnames = trim(varnames)//' "Pass('//trim(str(.true.,i))//')"'
    enddo
    varnames = trim(varnames)//trim(obj_io%ORANSname)//' "dt" "beta"'

    ! Gather R, dtlocal, beta from all ranks to root (collective)
    call gather_diagnostic_to_root(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'

      call parse(obj_io%sol_format,' ', format)

      ! The diagnostic reuses the solution-output buffer: nprim+1 slots, two more with thermo, two or three more with
      ! transport. It needs nres+2 = nprim-nsc+3 slots, which are missing only with a single species and neither
      ! thermo nor transport in sol-variables.
      if (size(IOfield%block(1)%vars,1) < ndiag) &
        call mpi_abort_all('the diagnostic output of a steady run with a single species needs sol-variables to '// &
                           'include thermo (the default) or transport')

      ! Update IOfield variables with domain residuals
      do b = 1, size(IOfield%block)
        IOfield%block(b)%vars(1,:,:,:) = sum(domain%blk(b)%r(1:nsc,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk), dim=1)
        IOfield%block(b)%vars(2:5,:,:,:) = domain%blk(b)%r(nu:np,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        if ( nres > 5 ) then ! soot, passive scalars and turbulence variables, in the order of r
          IOfield%block(b)%vars(6:nres,:,:,:) = &
            domain%blk(b)%r(np+1:nprim,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        end if
        ! Auxiliary variables: dt
        do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
              IOfield%block(b)%vars(nres+1,i,j,k) = domain%blk(b)%dtlocal(i,j,k)
        enddo; enddo; enddo
        ! Auxiliary variables: beta
        do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
              IOfield%block(b)%vars(nres+2,i,j,k) = domain%blk(b)%beta(i,j,k)
        enddo; enddo; enddo
      enddo

      ! Write the IOfield accordingly to the solution format
      select case(trim(format(1)))
      case('vtk')
        IOfield%vtk%format = trim(format(2))
        ! The names become XML attributes (Name="..."): the double quotes of the Tecplot list must go
        vtknames = varnames
        do i = 1, len_trim(vtknames)
          if ( vtknames(i:i) == '"' ) vtknames(i:i) = ' '
        enddo
        ! IOfield is the output buffer of the field, with obj_io%Onvar columns, and the vtk
        ! writer writes every column it holds: give it the ndiag diagnostic columns only,
        ! then the field's size back (the field writer fills every column at each write).
        do b = 1, size(IOfield%block)
          call move_alloc(IOfield%block(b)%vars, full)
          allocate(IOfield%block(b)%vars(1:ndiag,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk))
          IOfield%block(b)%vars(:,:,:,:) = full(1:ndiag,:,:,:)
          deallocate(full)
        enddo
        localpath_vtk = trim(path)//'/vtk'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(file), &
                                                             vtmpath=trim(path)//trim(file),varnames=vtknames,time=domain%time)
        do b = 1, size(IOfield%block)
          deallocate(IOfield%block(b)%vars)
          allocate(IOfield%block(b)%vars(1:obj_io%Onvar,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk))
        enddo
      case('tecplot')
        if (format(2)=='binary') then
          extension = '.szplt'
        else
          extension = '.tec'
        end if
        E_IO = tec_write_structured_multiblock(Nvars=ndiag,orion=IOfield,varnames=varnames,filename=trim(path)//trim(file)//trim(extension))
      end select
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_Diagnostic

end module MOSE_Mod_Diagnostic
