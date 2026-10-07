module MOSE_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_BC
  public :: Print_BC_Summary

  ! BC type counters — populated by Read_BCfile, consumed by Print_BC_Summary
  integer :: nconnect, nwall, ngsi, nSRM, nio, nsym, nper, nchimera, ncoupled, next, nKaffs, nstate
  logical :: has_tdep_bc

  ! Empty fields of a property line. A null value of list-directed input (two commas with only blanks between) leaves
  ! its variable as it was: undefined for a field without a default value. Read_BCfile sets those fields to EMPTY
  ! (IEMPTY for integers, a blank token for a direction) before the read, and Need_R, Need_I and Need_Inlet stop the run
  ! on a field still EMPTY after it; a field with a default value keeps it (docs/user/boundary-conditions.md).
  real(R8), parameter :: EMPTY  = -huge(1.0_R8)
  integer,  parameter :: IEMPTY = -huge(1)

contains

  subroutine Setup_BC ( domain )
    use MOSE_Advanced_Types_m
    use FLINT_Lib_Thermodynamic, only: species_names
    use MOSE_Config_Types_m,     only: obj_multigrid, obj_io_bc
    use MOSE_Global_m,           only: model
    use MOSE_IO_BC_Q2D,          only: Setup_Q2D_BC_Data
    use MOSE_Lib_GSI,            only: GSI_Initialize
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    ! Local
    integer :: m, error

    do m = 1, obj_multigrid%MGL

      !! Phase 1: Allocate and check
      call Allocate_BC ( domain(m) )

      if (allocated(obj_io_bc%viscous_flag)) deallocate(obj_io_bc%viscous_flag)
      allocate ( obj_io_bc%viscous_flag( 6 * domain(m) % nb, 6 ) )
      obj_io_bc%viscous_flag = .false.

      if (allocated(obj_io_bc%coupling_flag)) deallocate(obj_io_bc%coupling_flag)
      allocate ( obj_io_bc%coupling_flag( 6 * domain(m) % nb, 6 ) )
      obj_io_bc%coupling_flag = .false.

      error = Check_BC ( domain(m) % nbound, m )
      if (error /= 0) cycle

      !! Phase 2: Read
      call Read_BCfile ( domain(m) % bc, domain(m) % n_bf, m )

      call Pair_Connections ( domain(m), m )

      if (m == 1) call Setup_Q2D_BC_Data ( domain(m) )

      if ((any(obj_io_bc%viscous_flag)) .and. (model==0)) then
        obj_io_bc%error_message = '[ERROR] Viscous wall BCs assigned with Eulerian model'
      end if

    end do

    !! Phase 3: Initialize BCs for specific models
    ! GSI BCs
    if (ngsi > 0) then
      call GSI_Initialize(species_names)
    end if

  end subroutine Setup_BC


  subroutine Allocate_BC ( domain )
    use MOSE_Advanced_Types_m
    implicit none
    class(MOSE_domain_type) :: domain
    ! Local
    integer :: b, n, ni, nj, nk

    n = 0
    do b = 1, domain % nb
      ni = domain % blk(b) % dim(1)
      nj = domain % blk(b) % dim(2)
      nk = domain % blk(b) % dim(3)
      n = n + 2*nj*nk + 2*ni*nk + 2*nj*ni
    enddo

    domain % nbound = n
    allocate ( domain % bc ( domain % nbound ) )
    allocate ( domain % n_bf ( domain % nb, 6 ) )
    
  end subroutine Allocate_BC


  function Check_BC (n, level) result(ios)
    use MOSE_Config_Types_m, only: obj_io_bc
    use MOSE_Global_m,       only: MOSE_phase_prefix
    use IR_Precision,        only: str
    implicit none
    integer, intent(in)  :: n, level
    integer              :: ios
    integer              :: di(5), ti, unitfile, n_proof
    integer              :: c, ci, cii

    ios = 0

    ! Open file
    if (level == 1) then
      open(newunit=unitfile,file='INPUT/'//trim(MOSE_phase_prefix)//'bc.txt',status='old',iostat=ios,action='read')
    else
      open(newunit=unitfile,file='INPUT/'//trim(MOSE_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt',status='old',iostat=ios,action='read')
    endif
    if (ios/=0) then
      obj_io_bc%error_message = '[ERROR] Boundary condition file not found for grid '//trim(str(.true.,level))
      return
    endif

    ! Cheak BC file consistency
    ios = 0; n_proof = -1
    do while (ios==0)
      read( unitfile,*,iostat=ios ) di(1), di(2), di(3), di(4), di(5), ti
      ! Dispatch to skip the right number of extra lines.
      ! This list must stay in step with the property lines actually consumed by
      ! Read_BC below: if a type reads one there but is missing here, the scan
      ! mistakes its property line for a record and the whole file is rejected.
      ! IDs with ONE property line:
      !   101/103 = connection | 201 = periodic | 301,302 = wall | 401-410,420 = inlet/outlet
      !   501 = manifold | 502 = srm | 503-506 = gas-surface interaction
      ! IDs with NO property line:
      !   300 = symmetry | 400 = extrapolation | 409 = forced outlet
      ! IDs with VARIABLE-length property lines:
      !   102 = chimera | 104 = coupled multi-solver chimera
      select case(ti)
      case(101, 103, 201, 301, 302, 401:408, 410, 420, 501:506)
        read( unitfile,*,iostat=ios )
      case(102, 104)
        read( unitfile,*,iostat=ios ) ci, cii
        do c = 1, ci+cii
          read( unitfile,*,iostat=ios )
        enddo
      end select
      n_proof = n_proof + 1
    enddo

    if (n_proof /= n) then
      obj_io_bc%error_message = '[ERROR] Boundary conditions number ('//str(.true.,n_proof)//') is different than the one of the initial conditions ('//str(.true.,n)//')'
      close(unitfile)
      return
    endif

    ! Validation passed: reset ios (non-zero from EOF) to signal success
    ios = 0
    close(unitfile)

  end function Check_BC


  subroutine Read_BCfile ( bc, n_bf, level )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_io, obj_io_bc
    use MOSE_Global_m
    use IR_Precision
    implicit none
    type(MOSE_bc_type), dimension(:), intent(inout) :: bc
    integer, intent(in)                             :: level
    integer, dimension(1:,1:), intent(inout)        :: n_bf
    ! Local
    integer :: cc, i, s, nnozzle, nmanifold
    integer :: unitfile, ios, cios, ip
    integer :: nline                         ! lines of the file read so far (Next_Line)
    integer :: nzero_chim
    character(len=llen) :: p0file
    character(len=32) :: alpha_tok, beta_tok
    character(len=:), allocatable :: pline   ! a property line, read whole
    character(len=256) :: msg402
    character(len=256) :: q2d_line

    cios = 0

    ! Open file
    if (level == 1) then
      open(newunit=unitfile,file='INPUT/'//trim(MOSE_phase_prefix)//'bc.txt',status='old',iostat=ios,action='read')
    else
      open(newunit=unitfile,file='INPUT/'//trim(MOSE_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt',status='old',iostat=ios,action='read')
    endif
    if (ios/=0) then
      obj_io_bc%error_message = '[ERROR] Boundary condition file not found for grid '//trim(str(.true.,level))
      return
    endif

    ! Counters for specific BC types
    if (level == 1) then
      nconnect = 0
      nwall    = 0
      ngsi     = 0
      nSRM     = 0
      nio      = 0
      nsym     = 0
      nper     = 0
      nchimera = 0
      nnozzle  = 0
      ncoupled = 0
      nmanifold = 0
      next     = 0
      nKaffs   = 0
      nstate     = 0
      has_tdep_bc = .false.
    endif

    ! Counter for number of cells per face in each block
    n_bf = 0
    nzero_chim = 0

    ! Read file
    nline = 0
    do i = 1, size(bc)

      ! ── First line: block, ijk, face, ATLAS BC ID ─────────────────────────
      ! Every line is read whole (Next_Line, which counts the lines of the file and skips
      ! the blank ones as the list-directed read did), then as its values.
      call Next_Line ( unitfile, pline, ios, nline )
      if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%k, bc(i)%f, bc(i)%type
      if ( ios /= 0 ) call BC_Line_Error ( i, -1, nline, pline )

      ! n_bf update
      n_bf( bc(i) % b, bc(i) % f ) = n_bf( bc(i) % b, bc(i) % f ) + 1

      ! ── Second line: property data, BC-type dependent ────────────────────
      select case( bc(i) % type )

        ! ─────────────────────────────────────────────────────────────────────
        ! Connection and periodic BCs
        case(101, 201)
          if (level == 1) nconnect = nconnect + 1
          bc(i)%bs = IEMPTY ; bc(i)%is = IEMPTY ; bc(i)%js = IEMPTY ; bc(i)%ks = IEMPTY ; bc(i)%fs = IEMPTY
          bc(i)%d11 = IEMPTY ; bc(i)%d12 = IEMPTY ; bc(i)%d21 = IEMPTY ; bc(i)%d22 = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          if ( ios == 0 ) call Need_I ( i, bc(i)%type, nline, [ bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, &
                                        bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22 ], 'bs is js ks fs d11 d12 d21 d22' )
          allocate ( bc(i) % Pg (nprim, 6) )
          allocate ( bc(i) % Fc (nprim, 2) )
          bc(i) % Fc = 0d0

        ! ─────────────────────────────────────────────────────────────────────
        ! Chimera overlap BC (102) and coupled multi-solver chimera (104, donors are FUSS cells)
        case(102, 104)
          if (level == 1) nchimera = nchimera + 1
          bc(i)%ni = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) (bc(i)%ni(cc),cc=1,2)
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          call Need_I ( i, bc(i)%type, nline, bc(i)%ni, 'n1 n2' )
          allocate(bc(i)%donorID(1:sum(bc(i)%ni),1:4))
          allocate(bc(i)%volume_fraction(1:sum(bc(i)%ni)))
          bc(i)%donorID = IEMPTY
          bc(i)%volume_fraction = EMPTY
          do s = 1, bc(i)%ni(1)
            call Next_Line ( unitfile, pline, ios, nline )
            if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
            if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
            if ( ios == 0 ) call Need_I ( i, bc(i)%type, nline, bc(i)%donorID(s,1:4), 'bD iD jD kD' )
            if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, bc(i)%volume_fraction(s:s), 'vol_frac' )
          enddo
          do s = bc(i)%ni(1)+1, bc(i)%ni(1)+bc(i)%ni(2)
            call Next_Line ( unitfile, pline, ios, nline )
            if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
            if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
            if ( ios == 0 ) call Need_I ( i, bc(i)%type, nline, bc(i)%donorID(s,1:4), 'bD iD jD kD' )
            if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, bc(i)%volume_fraction(s:s), 'vol_frac' )
          enddo
          if (sum(bc(i)%volume_fraction(1:bc(i)%ni(1))) < 0.5d0 .or. &
              sum(bc(i)%volume_fraction(bc(i)%ni(1)+1:sum(bc(i)%ni))) < 0.5d0) then
            nzero_chim = nzero_chim + 1
          endif
          if (bc(i)%type == 102) then
            allocate ( bc(i) % Pg (nprim, 6) )
          else
            if (level == 1) ncoupled = ncoupled + 1
            ! 104: coupled wall as for 103 (solid T stencil, coupling flux, wall output)
            obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
            allocate ( bc(i) % Pg (1, 6) )
            allocate ( bc(i) % ext_flux (nprim) )
            bc(i) % ext_flux = 0.0
          endif

        ! ─────────────────────────────────────────────────────────────────────
        ! Axisymmetry
        case(200)
          if (level == 1) nper = nper + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Euler symmetry
        case(300)
          if (level == 1) nsym = nsym + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Extrapolation
        case(400)
          if (level == 1) next = next + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, prescribed heat flux
        ! Second line: q, roughness_ks, emissivity_eps  (comma-separated reals)
        case(301)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%qw = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%qw, bc(i)%k_rough, bc(i)%eps_wall
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%qw ], 'q_wall' )

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, prescribed temperature
        ! Second line: T, roughness_ks, emissivity_eps
        case(302)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%Tw = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%Tw, bc(i)%k_rough, bc(i)%eps_wall
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%Tw ], 'T_wall' )

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, stag. conditions (T0, p0)
        ! Second line: T0, p0, alpha, beta, rel_fac, massf(1:nsc), turb(1:nrans)
        case(401)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          bc(i)%T0 = EMPTY ; bc(i)%p0 = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%ci = EMPTY
          alpha_tok = ' ' ; beta_tok = ' '
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%T0, bc(i)%p0, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          if ( ios == 0 ) call Need_Inlet ( i, nline, bc(i), [ bc(i)%T0, bc(i)%p0 ], 'T0 p0', &
                                            alpha_tok, beta_tok )
          if ( ios == 0 ) call Need_Direction ( i, bc(i)%type, nline, alpha_tok, beta_tok, bc(i)%alpha, bc(i)%beta )

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, time-varying p0
        ! Second line: T0, p0_time_file, alpha, beta, rel_fac, massf, turb
        case(402)
          if (level == 1) nio = nio + 1
          has_tdep_bc = .true.
          allocate( bc(i) % ci(1 : nsc+nrans) )
          ! Read as one line and split by hand (Read_Inlet_402): read as a list, a file
          ! name with a directory ended at its slash, and the values after it were lost.
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) then
            call Read_Inlet_402 ( pline, bc(i), p0file, alpha_tok, beta_tok, ios, msg402 )
          else
            msg402 = 'the property line is missing'
          end if
          if ( ios /= 0 ) then
            block
              use MOSE_Mod_MPI, only: mpi_abort_all
              character(len=512) :: msg
              write(msg,'(A,I0,A,I0,A)') '[ERROR] bc.txt line ', nline, ' (record ', i, ', type 402): '//trim(msg402)
              call mpi_abort_all(trim(msg))
            end block
          end if
          call bc(i)%p0time%initialize(file=trim(p0file), bar=.true.)
          if ( ios == 0 ) call Need_Direction ( i, bc(i)%type, nline, alpha_tok, beta_tok, bc(i)%alpha, bc(i)%beta )

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, mass-flux g + T0
        ! Second line: T0, g, alpha, beta, rel_fac, massf, turb
        case(403)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          bc(i)%T0 = EMPTY ; bc(i)%mdot = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%ci = EMPTY
          alpha_tok = ' ' ; beta_tok = ' '
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%T0, bc(i)%mdot, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          if ( ios == 0 ) call Need_Inlet ( i, nline, bc(i), [ bc(i)%T0, bc(i)%mdot ], 'T0 g', &
                                            alpha_tok, beta_tok )
          if ( ios == 0 ) call Need_Direction ( i, bc(i)%type, nline, alpha_tok, beta_tok, bc(i)%alpha, bc(i)%beta )

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, mass-flux g + T (static)
        ! Second line: T, g, alpha, beta, rel_fac, massf, turb
        case(404)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          bc(i)%T0 = EMPTY ; bc(i)%mdot = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%ci = EMPTY
          alpha_tok = ' ' ; beta_tok = ' '
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%T0, bc(i)%mdot, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          if ( ios == 0 ) call Need_Inlet ( i, nline, bc(i), [ bc(i)%T0, bc(i)%mdot ], 'T g', &
                                            alpha_tok, beta_tok )
          if ( ios == 0 ) call Need_Direction ( i, bc(i)%type, nline, alpha_tok, beta_tok, bc(i)%alpha, bc(i)%beta )

        ! ─────────────────────────────────────────────────────────────────────
        ! Supersonic inlet, M + T (static) + p (static)
        ! Second line: mach, T, p, alpha, beta, rel_fac, massf, turb
        case(405)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          bc(i)%mach = EMPTY ; bc(i)%T0 = EMPTY ; bc(i)%pamb = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%ci = EMPTY
          alpha_tok = ' ' ; beta_tok = ' '
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%mach, bc(i)%T0, bc(i)%pamb, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          if ( ios == 0 ) call Need_Inlet ( i, nline, bc(i), [ bc(i)%mach, bc(i)%T0, bc(i)%pamb ], 'mach T p', &
                                            alpha_tok, beta_tok )
          if ( ios == 0 ) call Need_Direction ( i, bc(i)%type, nline, alpha_tok, beta_tok, bc(i)%alpha, bc(i)%beta )

        ! ─────────────────────────────────────────────────────────────────────
        ! Outlet, prescribed back pressure (if zero, extrapolation)
        ! Second line: p_back, rel_fac
        case(406)
          if (level == 1) nio = nio + 1
          bc(i)%pamb = EMPTY ; bc(i)%rel_fac = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%pamb, bc(i)%rel_fac
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%pamb, bc(i)%rel_fac ], 'p_back rel_fac' )

        ! ─────────────────────────────────────────────────────────────────────
        ! Ambient, T0 + p0 + back pressure
        ! Second line: T0, p0, p_back, alpha, beta, rel_fac, massf, turb
        case(407)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          bc(i)%T0 = EMPTY ; bc(i)%p0 = EMPTY ; bc(i)%pamb = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%ci = EMPTY
          alpha_tok = ' ' ; beta_tok = ' '
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%T0, bc(i)%p0, bc(i)%pamb, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          if ( ios == 0 ) call Need_Inlet ( i, nline, bc(i), [ bc(i)%T0, bc(i)%p0, bc(i)%pamb ], 'T0 p0 p_back', &
                                            alpha_tok, beta_tok )
          if ( ios == 0 ) call Need_Direction ( i, bc(i)%type, nline, alpha_tok, beta_tok, bc(i)%alpha, bc(i)%beta )

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, mass-flux un + T (static)
        ! Second line: T, un, alpha, beta, rel_fac, massf, turb
        case(408)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          bc(i)%T0 = EMPTY ; bc(i)%un = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%ci = EMPTY
          alpha_tok = ' ' ; beta_tok = ' '
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%T0, bc(i)%un, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          if ( ios == 0 ) call Need_Inlet ( i, nline, bc(i), [ bc(i)%T0, bc(i)%un ], 'T un', &
                                            alpha_tok, beta_tok )
          if ( ios == 0 ) call Need_Direction ( i, bc(i)%type, nline, alpha_tok, beta_tok, bc(i)%alpha, bc(i)%beta )

        ! ─────────────────────────────────────────────────────────────────────
        ! Assigned state via time-varying file
        case(410)
          if (level == 1) nstate = nstate + 1
          allocate ( bc(i) % Pg (nprim, 6) )
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 .and. len_trim(pline) > len(q2d_line) ) ios = 1   ! a name that does not fit q2d_file
          q2d_line = pline
          q2d_line = adjustl(q2d_line)
          ip = index(trim(q2d_line), ',')
          if (ip > 0) then
            bc(i)%q2d_file = q2d_line(1:ip-1)
            if (index(q2d_line(ip:), 'periodic') > 0) bc(i)%q2d_periodic = .true.
          else
            bc(i)%q2d_file = trim(q2d_line)
          endif

        ! ─────────────────────────────────────────────────────────────────────
        ! Choked-nozzle BC
        ! Second line: T0, p0, psub, psup, g, rel_fac, yi, turb
        case(420)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          bc(i)%T0 = EMPTY ; bc(i)%p0 = EMPTY ; bc(i)%mdot = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%ci = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%T0, bc(i)%p0, bc(i)%psub, bc(i)%psup, bc(i)%mdot, bc(i)%rel_fac, &
            (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          if ( ios == 0 ) then
            call Need_R ( i, bc(i)%type, nline, [ bc(i)%T0, bc(i)%p0 ], 'T0 p0' )
            call Need_R ( i, bc(i)%type, nline, [ bc(i)%mdot, bc(i)%rel_fac ], 'g rel_fac' )
            call Need_R ( i, bc(i)%type, nline, bc(i)%ci(nrans+1:nrans+nsc), 'massf' )
            call Need_R ( i, bc(i)%type, nline, bc(i)%ci(1:nrans), 'turb' )
          end if

        ! ─────────────────────────────────────────────────────────────────────
        ! Manifold
        case(501)
          if (level == 1) nmanifold = nmanifold + 1
          bc(i)%bs = IEMPTY ; bc(i)%fs = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%bs, bc(i)%fs
          if ( ios == 0 ) call Need_I ( i, bc(i)%type, nline, [ bc(i)%bs, bc(i)%fs ], 'bs fs' )

        ! ─────────────────────────────────────────────────────────────────────
        ! SRM grain boundary
        ! Second line: Taf, a, n, pRef, rhoGrain, SF, massf(1:nsc) [, turb...]
        case(502)
          if (level == 1) nSRM = nSRM + 1
          allocate( bc(i) % ci(1 : nsc) )
          ! The enthalpy of the injected products is no longer read here: it used to
          ! come from ATLAS and was left at zero when that output dropped it, which
          ! silently injected the propellant gas cold. It is now evaluated from Taf
          ! and the product composition against the thermodynamic tables, in BC_SRM.
          bc(i)%Taf = EMPTY ; bc(i)%aCoeff = EMPTY ; bc(i)%n = EMPTY ; bc(i)%pRef = EMPTY ; bc(i)%ci = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%Taf, bc(i)%aCoeff, bc(i)%n, bc(i)%pRef, &
            bc(i)%rhoGrain, bc(i)%SF_geo, (bc(i)%ci(s), s = 1, nsc)
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%Taf, bc(i)%aCoeff, bc(i)%n, bc(i)%pRef ], &
                                        'T_af a n p_ref' )
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, bc(i)%ci, 'massf' )

        ! ─────────────────────────────────────────────────────────────────────
        ! GSI - Melting
        ! Second line: specific_heat, temperature, initial_temperature, heat_of_fusion, radiative_flux, emissivity, mass fraction(1:nsc)
        case(503)
          if (level == 1) nwall  = nwall  + 1
          if (level == 1) ngsi   = ngsi   + 1
          allocate( bc(i) % ci(1 : nsc) )
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%Tw = EMPTY ; bc(i)%qrad = EMPTY ; bc(i)%ci = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%cp_wall, bc(i)%Tw, bc(i)%Ti_wall, bc(i)%dh_wall, bc(i)%qrad, bc(i)%eps_wall, (bc(i)%ci(s), s = 1, nsc)
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%Tw, bc(i)%qrad ], 'T_wall q_rad' )
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, bc(i)%ci, 'massf' )

        ! ─────────────────────────────────────────────────────────────────────
        ! GSI - Pyrolysis
        ! Second line: pyrolysis model id, radiative_flux, emissivity_eps, mass fraction(1:nsc)
        case(504)
          if (level == 1) nwall  = nwall  + 1
          if (level == 1) ngsi   = ngsi   + 1
          allocate( bc(i) % ci(1 : nsc) )
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%qrad = EMPTY ; bc(i)%ci = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%GSI_pyro_model_id, bc(i)%qrad, bc(i)%eps_wall, (bc(i)%ci(s), s = 1, nsc)
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%qrad ], 'q_rad' )
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, bc(i)%ci, 'massf' )

        ! ─────────────────────────────────────────────────────────────────────
        ! GSI - Surface reactions
        ! Second line: surface reaction id, radiative flux
        case(505)
          if (level == 1) nwall  = nwall  + 1
          if (level == 1) ngsi   = ngsi   + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%qrad = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%GSI_surf_reac_id, bc(i)%qrad
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%qrad ], 'q_rad' )

        ! ─────────────────────────────────────────────────────────────────────
        ! GSI - Surface reactions + pyrolysis
        ! Second line: pyrolysis model id, surface reaction id, radiative flux
        case(506)
          if (level == 1) nwall  = nwall  + 1
          if (level == 1) ngsi   = ngsi   + 1
          allocate( bc(i) % ci(1 : nsc) )
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%qrad = EMPTY ; bc(i)%ci = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%GSI_pyro_model_id, bc(i)%GSI_surf_reac_id, bc(i)%qrad, bc(i)%eps_wall, (bc(i)%ci(s), s = 1, nsc)
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%qrad ], 'q_rad' )
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, bc(i)%ci, 'massf' )

        ! ─────────────────────────────────────────────────────────────────────
        ! Coupled multi-solver wall
        case(103)
          if (level == 1) ncoupled = ncoupled + 1
          obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%bs = IEMPTY ; bc(i)%is = IEMPTY ; bc(i)%js = IEMPTY ; bc(i)%ks = IEMPTY ; bc(i)%fs = IEMPTY
          bc(i)%d11 = IEMPTY ; bc(i)%d12 = IEMPTY ; bc(i)%d21 = IEMPTY ; bc(i)%d22 = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          if ( ios == 0 ) call Need_I ( i, bc(i)%type, nline, [ bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, &
                                        bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22 ], 'bs is js ks fs d11 d12 d21 d22' )
          allocate(bc(i)%ext_flux(nprim))
          bc(i)%ext_flux = 0.0
          allocate ( bc(i) % Pg (1, 6) )

        ! ─────────────────────────────────────────────────────────────────────
        ! ! MOSE 999: Kaffs 1D–multi-D hybrid coupling (unchanged)
        ! case(999)
        !   if (level == 1) nKaffs = nKaffs + 1
        !   read(unitfile,*,iostat=ios) bc(i)%bs, bc(i)%fs

      end select

      ! A property line that cannot be read as the values of its type stops the run
      if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )

    enddo

    if (nwall > 0) obj_io % write_wall = .true.

    close( unitfile )

    if (nzero_chim > 0) then
      block
        use MOSE_Mod_MPI, only: mpi_abort_all
        character(len=256) :: msg
        write(msg,'(I0,A,I0,A)') nzero_chim, ' chimera entries in the level-', level, &
          ' BC file have donor volume fractions summing to ~0 (no donors found by the BC '// &
          'builder). Their ghost states would be garbage; fix the BC file before running.'
        call mpi_abort_all(trim(msg))
      end block
    endif

  end subroutine Read_BCfile


  !> Pairs the two records of each connection face (101, 201): the record of the
  !> same face seen from the other block (partner), and which of the two computes
  !> the flux of the face for both cells (owner, BC_Fluxes). The owner is the
  !> record whose face is the high-index one of its cell (2, 4, 6) when the
  !> other is a low-index one, as on an interior face, where the flux between
  !> cells n and n+1 uses the shock sensor of cell n. When the two faces have the
  !> same parity (blocks that meet with their index directions permuted), the
  !> owner is the record with the smaller (block, face, k, j, i): a choice that
  !> does not depend on the order of the records in bc.txt. A record whose
  !> partner does not point back at it, or that points at itself, is an error.
  !> The search goes through a map, per block face, from the two indices of the
  !> cell along the face to its connection record: one pass over the records,
  !> in any order.
  !> The gc ghost cells of a connection face are the first gc cells of the
  !> block across it (Ghost_Connection), and their metrics come from its first
  !> gc node planes past the face (BC_Connect_Metrics). That block must have at
  !> least gc cells across its face, at every multigrid level (level): with
  !> fewer, a ghost cell would be a ghost cell of that block, which another
  !> record fills at the same time, and a node plane would lie outside the
  !> block. A record whose block across has fewer cells is an error; the
  !> directions along the face, and those of blocks with no connection, need
  !> nothing.
  subroutine Pair_Connections ( domain, level )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_io_bc
    use MOSE_Global_m,       only: gc
    use IR_Precision, only: str
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in)                   :: level
    ! Local
    type :: face_records_type
      integer, allocatable :: r(:,:)
    end type face_records_type
    type(face_records_type), allocatable :: rec(:,:)
    integer :: b, f, i, p, nt(2), t(2), n
    character(len=32) :: at_level

    allocate ( rec(domain % nb, 6) )

    ! Map of the connection records of each block face
    do i = 1, domain % nbound
      if ( domain % bc(i) % type /= 101 .and. domain % bc(i) % type /= 201 ) cycle
      b = domain % bc(i) % b
      f = domain % bc(i) % f
      if ( .not. allocated(rec(b,f) % r) ) then
        nt = Tangential ( f, domain % blk(b) % dim(1), domain % blk(b) % dim(2), domain % blk(b) % dim(3) )
        allocate ( rec(b,f) % r(nt(1), nt(2)) )
        rec(b,f) % r = 0
      end if
      t = Tangential ( f, domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k )
      rec(b,f) % r(t(1), t(2)) = i
    end do

    ! Partner and owner of each connection record
    do i = 1, domain % nbound
      if ( domain % bc(i) % type /= 101 .and. domain % bc(i) % type /= 201 ) cycle
      b = domain % bc(i) % bs
      f = domain % bc(i) % fs
      p = 0
      if ( b >= 1 .and. b <= domain % nb .and. f >= 1 .and. f <= 6 ) then
        if ( allocated(rec(b,f) % r) ) then
          t = Tangential ( f, domain % bc(i) % is, domain % bc(i) % js, domain % bc(i) % ks )
          if ( all(t >= 1) .and. all(t <= shape(rec(b,f) % r)) ) p = rec(b,f) % r(t(1), t(2))
        end if
      end if
      if ( p == 0 ) then
        obj_io_bc % error_message = '[ERROR] Connection record '//trim(str(.true.,i))// &
                                    ': no connection record at the cell and face it points at'
        exit
      end if
      if ( p == i ) then
        obj_io_bc % error_message = '[ERROR] Connection record '//trim(str(.true.,i))// &
                                    ' points at its own cell and face'
        exit
      end if
      if ( domain % bc(p) % type /= domain % bc(i) % type .or. domain % bc(p) % bs /= domain % bc(i) % b .or. &
           domain % bc(p) % fs /= domain % bc(i) % f .or. domain % bc(p) % is /= domain % bc(i) % i .or. &
           domain % bc(p) % js /= domain % bc(i) % j .or. domain % bc(p) % ks /= domain % bc(i) % k ) then
        obj_io_bc % error_message = '[ERROR] Connection records '//trim(str(.true.,i))//' and '//trim(str(.true.,p))// &
                                    ' do not point at each other'
        exit
      end if
      n = domain % blk(b) % dim((f + 1) / 2)
      if ( n < gc ) then
        at_level = ''
        if ( level > 1 ) at_level = ' at multigrid level '//trim(str(.true.,level))
        obj_io_bc % error_message = '[ERROR] Connection record '//trim(str(.true.,i))//' (block '// &
          trim(str(.true.,domain % bc(i) % b))//', face '//trim(str(.true.,domain % bc(i) % f))//'): block '// &
          trim(str(.true.,b))//' across it has '//trim(str(.true.,n))//' cell across its face '//trim(str(.true.,f))// &
          trim(at_level)//', fewer than the '//trim(str(.true.,gc))//' a connection needs'
        exit
      end if
      domain % bc(i) % partner = p
      if ( mod(domain % bc(i) % f, 2) == 0 .and. mod(domain % bc(p) % f, 2) == 1 ) then
        domain % bc(i) % owner = .true.
      else if ( mod(domain % bc(i) % f, 2) == 1 .and. mod(domain % bc(p) % f, 2) == 0 ) then
        domain % bc(i) % owner = .false.
      else
        domain % bc(i) % owner = Key_Less ( domain % bc(i) % b, domain % bc(i) % f, domain % bc(i) % k, &
                                            domain % bc(i) % j, domain % bc(i) % i, &
                                            domain % bc(p) % b, domain % bc(p) % f, domain % bc(p) % k, &
                                            domain % bc(p) % j, domain % bc(p) % i )
      end if
    end do

    contains

      !> Whether the record (b1, f1, k1, j1, i1) comes before (b2, f2, k2, j2, i2):
      !> block, then face, then k, j, i. Every record has its own key, so the two
      !> records of a face always tell which of them is first, whatever their
      !> order in bc.txt.
      pure logical function Key_Less ( b1, f1, k1, j1, i1, b2, f2, k2, j2, i2 )
        implicit none
        integer, intent(in) :: b1, f1, k1, j1, i1, b2, f2, k2, j2, i2
        integer :: key1(5), key2(5), n

        key1 = [ b1, f1, k1, j1, i1 ]
        key2 = [ b2, f2, k2, j2, i2 ]
        Key_Less = .false.
        do n = 1, 5
          if ( key1(n) /= key2(n) ) then
            Key_Less = key1(n) < key2(n)
            return
          end if
        end do

      end function Key_Less

      !> The two indices of a cell (or the two dimensions of a block) along the face f
      pure function Tangential ( f, i, j, k ) result ( t )
        implicit none
        integer, intent(in) :: f, i, j, k
        integer :: t(2)

        select case ( f )
          case (1:2)
            t = [ j, k ]
          case (3:4)
            t = [ i, k ]
          case default
            t = [ i, j ]
        end select

      end function Tangential

  end subroutine Pair_Connections


  subroutine Print_BC_Summary ()
    implicit none

    write(*,*)
    write(*,'(A)') ' Boundary conditions'
    if (nconnect > 0) write(*,'(A,T35,I0)') '   Connection', nconnect
    if (nwall > 0) write(*,'(A,T35,I0)') '   Viscous wall', nwall
    if (ngsi > 0) write(*,'(A,T35,I0)') '   GSI (inc. wall)', ngsi
    if (nio > 0) write(*,'(A,T35,I0)') '   Inflow/outflow', nio
    if (nsym > 0) write(*,'(A,T35,I0)') '   Symmetry', nsym
    if (nper > 0) write(*,'(A,T35,I0)') '   Periodicity', nper
    if (next > 0) write(*,'(A,T35,I0)') '   Extrapolation', next
    if (nchimera > 0) write(*,'(A,T35,I0)') '   Chimera', nchimera
    if (ncoupled > 0) write(*,'(A,T35,I0)') '   Coupled wall', ncoupled
    if (nSRM > 0) write(*,'(A,T35,I0)') '   SRM wall', nSRM
    !if (nKaffs > 0) write(*,'(A,T35,I0)') '   Coupled 1D-MultiD', nKaffs
    if (nstate > 0) write(*,'(A,T35,I0)') '   Q2D Mapped', nstate
    if (has_tdep_bc) write(*,'(A)') '   Time-dependent BC detected'

  end subroutine Print_BC_Summary


  !─────────────────────────────────────────────────────────────────────────────
  ! Next_Line: the next line of the boundary file that is not blank (Read_Line), as the list-directed read skipped
  !   blank lines; nline counts every line read, blank ones included, so it is the line number in the file.
  subroutine Next_Line ( unit, line, ios, nline )
    implicit none
    integer, intent(in)                        :: unit
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out)                       :: ios
    integer, intent(inout)                     :: nline

    do
      call Read_Line ( unit, line, ios )
      if ( ios /= 0 ) exit
      nline = nline + 1
      if ( len_trim(line) > 0 ) exit
    end do

  end subroutine Next_Line


  !─────────────────────────────────────────────────────────────────────────────
  ! BC_Line_Error: stop the run on line nline of the boundary file, of record n, that cannot be read as the values it
  !   must hold (short, or a value that is not a number): its record line (btype < 0) or a property line of type btype;
  !   an empty line is the end of the file.
  subroutine BC_Line_Error ( n, btype, nline, line )
    use MOSE_Mod_MPI, only: mpi_abort_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    character(len=*), intent(in) :: line
    ! Local
    character(len=160) :: head

    if ( len_trim(line) == 0 ) then
      write(head,'(A,I0,A,I0)') '[ERROR] bc.txt ends at line ', nline, ', before the end of record ', n
      call mpi_abort_all ( trim(head) )
    end if
    if ( btype < 0 ) then
      write(head,'(A,I0,A,I0,A)') '[ERROR] bc.txt line ', nline, ' (record ', n, &
        '): the line is not block, i, j, k, face and type'
    else
      write(head,'(A,I0,A,I0,A,I0,A)') '[ERROR] bc.txt line ', nline, ' (record ', n, ', type ', btype, &
        '): the line is short or holds a value that is not a number'
    end if
    call mpi_abort_all ( trim(head)//': "'//trim(line(1:min(len(line),200)))//'"' )

  end subroutine BC_Line_Error


  !─────────────────────────────────────────────────────────────────────────────
  ! Need_R, Need_I: stop the run (Empty_Field) on the first of the values v of record n, type btype, line nline of the
  !   file, still EMPTY after the read of its property line. names: their names, blank-separated, in the order of the
  !   line; the last name stands for the values after it (a list of mass fractions, say).
  subroutine Need_R ( n, btype, nline, v, names )
    implicit none
    integer, intent(in)          :: n, btype, nline
    real(R8), intent(in)         :: v(:)
    character(len=*), intent(in) :: names
    ! Local
    integer :: k

    do k = 1, size(v)
      if ( v(k) == EMPTY ) call Empty_Field ( n, btype, nline, Name_Of ( names, k ) )
    end do

  end subroutine Need_R


  subroutine Need_I ( n, btype, nline, v, names )
    implicit none
    integer, intent(in)          :: n, btype, nline
    integer, intent(in)          :: v(:)
    character(len=*), intent(in) :: names
    ! Local
    integer :: k

    do k = 1, size(v)
      if ( v(k) == IEMPTY ) call Empty_Field ( n, btype, nline, Name_Of ( names, k ) )
    end do

  end subroutine Need_I


  !─────────────────────────────────────────────────────────────────────────────
  ! Need_Inlet: the fields of an inlet line, in its order: the leading values lead (names), alpha, beta, rel_fac, the
  !   mass fractions and the turbulence values.
  subroutine Need_Inlet ( n, nline, bc, lead, names, alpha_tok, beta_tok )
    use MOSE_Advanced_Types_m, only: MOSE_bc_type
    use MOSE_Global_m,         only: nsc, nrans
    implicit none
    integer, intent(in)            :: n, nline
    type(MOSE_bc_type), intent(in) :: bc
    real(R8), intent(in)           :: lead(:)
    character(len=*), intent(in)   :: names, alpha_tok, beta_tok

    call Need_R ( n, bc%type, nline, lead, names )
    if ( len_trim(alpha_tok) == 0 ) call Empty_Field ( n, bc%type, nline, 'alpha' )
    if ( len_trim(beta_tok) == 0 ) call Empty_Field ( n, bc%type, nline, 'beta' )
    call Need_R ( n, bc%type, nline, [ bc%rel_fac ], 'rel_fac' )
    call Need_R ( n, bc%type, nline, bc%ci(nrans+1:nrans+nsc), 'massf' )
    call Need_R ( n, bc%type, nline, bc%ci(1:nrans), 'turb' )

  end subroutine Need_Inlet


  !─────────────────────────────────────────────────────────────────────────────
  ! Empty_Field: stop the run on the field what of the property line of record n, type btype, line nline of the file,
  !   left empty.
  subroutine Empty_Field ( n, btype, nline, what )
    use MOSE_Mod_MPI, only: mpi_abort_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    character(len=*), intent(in) :: what
    ! Local
    character(len=240) :: msg

    write(msg,'(A,I0,A,I0,A,I0,A)') '[ERROR] bc.txt line ', nline, ' (record ', n, ', type ', btype, &
      '): the value of '//trim(what)//' is empty (docs/user/boundary-conditions.md)'
    call mpi_abort_all ( trim(msg) )

  end subroutine Empty_Field


  !─────────────────────────────────────────────────────────────────────────────
  ! Name_Of: the k-th blank-separated word of names, or the last one when there are fewer.
  pure function Name_Of ( names, k ) result ( w )
    implicit none
    character(len=*), intent(in) :: names
    integer, intent(in)          :: k
    character(len=len(names))    :: w
    ! Local
    integer :: p, q, m

    w = ' '
    p = 1
    m = 0
    do
      do while ( p <= len_trim(names) )
        if ( names(p:p) /= ' ' ) exit
        p = p + 1
      end do
      if ( p > len_trim(names) ) exit
      q = index ( names(p:), ' ' )
      if ( q == 0 ) then
        q = len(names) + 1
      else
        q = p + q - 1
      end if
      m = m + 1
      w = names(p:q-1)
      if ( m == k ) exit
      p = q
    end do

  end function Name_Of


  !─────────────────────────────────────────────────────────────────────────────
  ! Read_Line: the next line of a formatted file, whole, whatever its length (read in
  !   pieces without advancing). ios < 0 at the end of the file; a last line without its
  !   end-of-line mark is still a line.
  subroutine Read_Line ( unit, line, ios )
    implicit none
    integer, intent(in)                        :: unit
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out)                       :: ios
    ! Local
    character(len=256) :: piece
    integer :: nread

    line = ''
    do
      read ( unit, '(A)', advance='no', iostat=ios, size=nread ) piece
      line = line // piece(1:nread)
      if ( ios /= 0 ) exit
    end do
    if ( is_iostat_eor(ios) ) ios = 0
    if ( is_iostat_end(ios) .and. len(line) > 0 ) ios = 0

  end subroutine Read_Line


  !─────────────────────────────────────────────────────────────────────────────
  ! Read_Inlet_402: the property line of an inlet 402 (time-varying p0), read as one line and split here into its
  !   values: T0, p0_file, alpha, beta, rel_fac, mass fractions (nsc), turbulence values (nrans). As in list-directed
  !   input, values are separated by a comma or by blanks (blanks around a comma included), and a value in quotes
  !   (' or ") is taken whole, separators included; unlike list-directed input, a slash is part of a value, so that
  !   the file name is taken as written: read as a list, a name with a directory (INPUT/p0.dat) ended at the slash, a
  !   separator that also ends the record, and the values after it were not read. Values after the last one needed
  !   are ignored, as a list-directed read did. ios /= 0, with msg, when the line is short, a value is empty (two
  !   commas with only blanks between: a null value, which would leave the field as it was) or is not a number, or
  !   the file name does not fit p0file.
  subroutine Read_Inlet_402 ( line, bc, p0file, alpha_tok, beta_tok, ios, msg )
    use MOSE_Advanced_Types_m, only: MOSE_bc_type
    use MOSE_Global_m, only: nsc, nrans
    implicit none
    character(len=*), intent(in)       :: line
    type(MOSE_bc_type), intent(inout)  :: bc
    character(len=*), intent(out)      :: p0file, alpha_tok, beta_tok
    integer, intent(out)               :: ios
    character(len=*), intent(out)      :: msg
    ! Local
    character(len=len(line)) :: tok(5 + nsc + nrans)
    integer :: ntok, need, p, n, s
    character(len=1) :: q

    ios = 0
    msg = ' '
    need = 5 + nsc + nrans

    ! Split the line into values, up to the last one needed
    ntok = 0
    p = 1
    n = len_trim(line)
    do while ( ntok < need )
      ! blanks before a value
      do while ( p <= n )
        if ( line(p:p) /= ' ' .and. line(p:p) /= achar(9) ) exit
        p = p + 1
      end do
      if ( p > n ) exit
      ntok = ntok + 1
      tok(ntok) = ' '
      if ( line(p:p) == ',' ) then
        ios = 1
        write(msg,'(A,I0,A)') 'value ', ntok, ' is empty'
        return
      else if ( line(p:p) == '"' .or. line(p:p) == "'" ) then
        q = line(p:p)
        s = index ( line(p+1:n), q )
        if ( s == 0 ) then
          ios = 1
          write(msg,'(A,I0,A)') 'the quotes of value ', ntok, ' are not closed'
          return
        end if
        tok(ntok) = line(p+1:p+s-1)
        p = p + s + 1
      else
        s = p
        do while ( p <= n )
          if ( line(p:p) == ' ' .or. line(p:p) == ',' .or. line(p:p) == achar(9) ) exit
          p = p + 1
        end do
        tok(ntok) = line(s:p-1)
      end if
      ! blanks after the value, then at most one comma
      do while ( p <= n )
        if ( line(p:p) /= ' ' .and. line(p:p) /= achar(9) ) exit
        p = p + 1
      end do
      if ( p <= n ) then
        if ( line(p:p) == ',' ) p = p + 1
      end if
    end do

    if ( ntok < need ) then
      ios = 1
      write(msg,'(A,I0,A,I0,A)') 'it holds ', ntok, ' values, the inlet needs ', need, &
        ' (T0, p0 file, alpha, beta, rel_fac, mass fractions, turbulence values)'
      return
    end if

    read ( tok(1), *, iostat=ios ) bc % T0
    if ( ios /= 0 ) then
      msg = 'T0 is not a number: '//trim(tok(1))
      return
    end if
    if ( len_trim(tok(2)) > len(p0file) ) then
      ios = 1
      write(msg,'(A,I0,A)') 'the name of the p0 file is longer than ', len(p0file), ' characters'
      return
    end if
    p0file    = tok(2)
    alpha_tok = tok(3)
    beta_tok  = tok(4)
    read ( tok(5), *, iostat=ios ) bc % rel_fac
    if ( ios /= 0 ) then
      msg = 'rel_fac is not a number: '//trim(tok(5))
      return
    end if
    do s = 1, nsc
      read ( tok(5+s), *, iostat=ios ) bc % ci(nrans+s)
      if ( ios /= 0 ) then
        msg = 'a mass fraction is not a number: '//trim(tok(5+s))
        return
      end if
    end do
    do s = 1, nrans
      read ( tok(5+nsc+s), *, iostat=ios ) bc % ci(s)
      if ( ios /= 0 ) then
        msg = 'a turbulence value is not a number: '//trim(tok(5+nsc+s))
        return
      end if
    end do

  end subroutine Read_Inlet_402


  !─────────────────────────────────────────────────────────────────────────────
  ! Need_Direction: the flow directions alpha and beta of an inlet line (record n, type btype, line nline of bc.txt)
  !   from their tokens; a token that is neither a number nor normal stops the run.
  subroutine Need_Direction ( n, btype, nline, alpha_tok, beta_tok, alpha, beta )
    use MOSE_Mod_MPI, only: mpi_abort_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    character(len=*), intent(in) :: alpha_tok, beta_tok
    real(R8), intent(out)        :: alpha, beta
    ! Local
    integer            :: ios
    character(len=200) :: msg

    call parse_dir_tok ( alpha_tok, alpha, ios )
    if ( ios /= 0 ) then
      write(msg,'(A,I0,A,I0,A,I0,A)') '[ERROR] bc.txt line ', nline, ' (record ', n, ', type ', btype, &
        "): alpha is neither a number nor normal: '"//trim(adjustl(alpha_tok))//"'"
      call mpi_abort_all ( trim(msg) )
    end if
    call parse_dir_tok ( beta_tok, beta, ios )
    if ( ios /= 0 ) then
      write(msg,'(A,I0,A,I0,A,I0,A)') '[ERROR] bc.txt line ', nline, ' (record ', n, ', type ', btype, &
        "): beta is neither a number nor normal: '"//trim(adjustl(beta_tok))//"'"
      call mpi_abort_all ( trim(msg) )
    end if

  end subroutine Need_Direction


  !─────────────────────────────────────────────────────────────────────────────
  ! parse_dir_tok: convert a direction token (read as character) to real.
  !   If the token contains 'normal', return huge (face-normal flag).
  !   Otherwise parse it as a floating-point number; ios /= 0 when it is not one (it gave 0 in silence).
  pure subroutine parse_dir_tok ( tok, val, ios )
    implicit none
    character(len=*), intent(in) :: tok
    real(R8), intent(out)        :: val
    integer, intent(out)         :: ios
    character(len=len(tok)) :: tok_
    tok_ = adjustl(tok)
    ios = 0
    if (index(trim(tok_), 'normal') > 0) then
      val = huge(1.0_R8)
    else
      read(tok_, *, iostat=ios) val
      if (ios /= 0) val = 0.0_R8
    endif
  end subroutine parse_dir_tok

end module MOSE_IO_BC
