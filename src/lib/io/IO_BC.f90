module MOSE_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_BC
  public :: Print_BC_Summary

  ! BC type counters — populated by Read_BCfile, consumed by Print_BC_Summary
  integer :: nconnect, nwall, ngsi, nSRM, nio, nsym, nper, nchimera, ncoupled, next, nKaffs, nstate
  logical :: has_tdep_bc

contains

  subroutine Setup_BC ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_multigrid, obj_io_bc
    use MOSE_Global_m, only: model
    use MOSE_IO_BC_Q2D, only: Setup_Q2D_BC_Data
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

      if (m == 1) call Setup_Q2D_BC_Data ( domain(m) )

      if ((any(obj_io_bc%viscous_flag)) .and. (model==0)) then
        obj_io_bc%error_message = '[ERROR] Viscous wall BCs assigned with Eulerian model'
      end if

    end do

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
      ! IDs with ONE property line:
      !   101/103 = connection | 201 = periodic | 301-309 = wall | 401-410,420 = inlet/outlet | 501 = manifold | 502 = srm 
      ! IDs with NO property line:
      !   300 = symmetry | 400 = extrapolation | 409 = forced outlet
      ! IDs with VARIABLE-length property lines:
      !   102 = chimera
      select case(ti)
      case(101, 103, 201, 301:309, 401:407, 410, 420, 501:502)
        read( unitfile,*,iostat=ios )
      case(102)
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
    character(len=32) :: p0file
    character(len=32) :: alpha_tok, beta_tok
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

    ! Read file
    do i = 1, size(bc)

      ! ── First line: block, ijk, face, ATLAS BC ID ─────────────────────────
      read( unitfile,*,iostat=ios ) bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%k, bc(i)%f, bc(i)%type
      if (ios/=0) write(*,'(A)') '  Error in BC file'

      ! n_bf update
      n_bf( bc(i) % b, bc(i) % f ) = n_bf( bc(i) % b, bc(i) % f ) + 1

      ! ── Second line: property data, BC-type dependent ────────────────────
      select case( bc(i) % type )

        ! ─────────────────────────────────────────────────────────────────────
        ! Connection and periodic BCs
        case(101, 201)
          if (level == 1) nconnect = nconnect + 1
          read( unitfile,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          allocate ( bc(i) % Pg (nprim, 6) )

        ! ─────────────────────────────────────────────────────────────────────
        ! Chimera overlap BC
        case(102)
          if (level == 1) nchimera = nchimera + 1
          read( unitfile,*,iostat=ios ) (bc(i)%ni(cc),cc=1,2)
          allocate(bc(i)%donorID(1:sum(bc(i)%ni),1:4))
          allocate(bc(i)%volume_fraction(1:sum(bc(i)%ni)))
          do s = 1, bc(i)%ni(1)
            read( unitfile,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
          enddo
          do s = bc(i)%ni(1)+1, bc(i)%ni(1)+bc(i)%ni(2)
            read( unitfile,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
          enddo
          allocate ( bc(i) % Pg (nprim, 6) )

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
          read(unitfile,*,iostat=ios) bc(i)%qw, bc(i)%k_rough, bc(i)%eps_wall

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, prescribed temperature
        ! Second line: T, roughness_ks, emissivity_eps
        case(302)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) bc(i)%Tw, bc(i)%k_rough, bc(i)%eps_wall

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, temperature + radiative flux
        ! Second line: T, qrad, roughness_ks
        case(303)
          if (level == 1) nwall  = nwall  + 1
          if (level == 1) ngsi   = ngsi   + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) bc(i)%Tw, bc(i)%qrad, bc(i)%k_rough

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, radiative flux
        ! Second line: qrad, roughness_ks
        case(304)
          if (level == 1) nwall  = nwall  + 1
          if (level == 1) ngsi   = ngsi   + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) bc(i)%qrad, bc(i)%k_rough

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, stag. conditions (T0, p0)
        ! Second line: T0, p0, alpha, beta, rel_fac, massf(1:nsc), turb(1:nrans)
        case(401)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          read( unitfile,*,iostat=ios ) &
            bc(i)%T0, bc(i)%p0, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, time-varying p0
        ! Second line: T0, p0_time_file, alpha, beta, rel_fac, massf, turb
        case(402)
          if (level == 1) nio = nio + 1
          has_tdep_bc = .true.
          allocate( bc(i) % ci(1 : nsc+nrans) )
          read( unitfile,*,iostat=ios ) &
            bc(i)%T0, p0file, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          call bc(i)%p0time%initialize(file=p0file, bar=.true.)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, mass-flux g + T0
        ! Second line: T0, g, alpha, beta, rel_fac, massf, turb
        case(403)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          read( unitfile,*,iostat=ios ) &
            bc(i)%T0, bc(i)%mdot, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, mass-flux g + T (static)
        ! Second line: T, g, alpha, beta, rel_fac, massf, turb
        case(404)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          read( unitfile,*,iostat=ios ) &
            bc(i)%T0, bc(i)%mdot, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Supersonic inlet, M + T (static) + p (static)
        ! Second line: mach, T, p, alpha, beta, rel_fac, massf, turb
        case(405)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          read( unitfile,*,iostat=ios ) &
            bc(i)%mach, bc(i)%T0, bc(i)%pamb, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Outlet, prescribed back pressure (if zero, extrapolation)
        ! Second line: p_back, rel_fac
        case(406)
          if (level == 1) nio = nio + 1
          read( unitfile,*,iostat=ios ) bc(i)%pamb, bc(i)%rel_fac

        ! ─────────────────────────────────────────────────────────────────────
        ! Ambient, T0 + p0 + back pressure
        ! Second line: T0, p0, p_back, alpha, beta, rel_fac, massf, turb
        case(407)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          read( unitfile,*,iostat=ios ) &
            bc(i)%T0, bc(i)%p0, bc(i)%pamb, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Assigned state via time-varying file
        case(410)
          if (level == 1) nstate = nstate + 1
          allocate ( bc(i) % Pg (nprim, 6) )
          read(unitfile, '(A)', iostat=ios) q2d_line
          q2d_line = adjustl(q2d_line)
          ip = index(trim(q2d_line), ' ')
          if (ip > 0) then
            bc(i)%q2d_file = q2d_line(1:ip-1)
            if (index(q2d_line(ip:), 'periodic') > 0) bc(i)%q2d_periodic = .true.
          else
            bc(i)%q2d_file = trim(q2d_line)
          endif

        ! ─────────────────────────────────────────────────────────────────────
        ! Choked-nozzle BC
        ! Second line: 0(alpha), T0, p0, psub, psup, rt, 0(rel_fac), massf, turb
        case(420)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : nsc+nrans) )
          read( unitfile,*,iostat=ios ) &
            bc(i)%mach, bc(i)%T0, bc(i)%p0, bc(i)%psub, bc(i)%psup, bc(i)%rt_nozzle, bc(i)%rel_fac, &
            (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%pamb = bc(i)%psub   ! pass subsonic pressure as back-pressure to nozzle solver

        ! ─────────────────────────────────────────────────────────────────────
        ! Manifold
        case(501)
          if (level == 1) nmanifold = nmanifold + 1
          read(unitfile,*,iostat=ios) bc(i)%bs, bc(i)%fs

        ! ─────────────────────────────────────────────────────────────────────
        ! SRM grain boundary
        ! Second line: Taf, a, n, pRef, rhoGrain, SF, massf(1:nsc) [, turb...]
        case(502)
          if (level == 1) nSRM = nSRM + 1
          allocate( bc(i) % ci(1 : nsc) )
          bc(i)%haf = 0.0_R8   ! haf no longer in ATLAS output; set to zero
          read(unitfile,*,iostat=ios) bc(i)%Taf, bc(i)%aCoeff, bc(i)%n, bc(i)%pRef, &
            bc(i)%rhoGrain, bc(i)%SF_geo, (bc(i)%ci(s), s = 1, nsc)

        ! ─────────────────────────────────────────────────────────────────────
        ! Coupled multi-solver wall
        case(103)
          if (level == 1) ncoupled = ncoupled + 1
          obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
          read( unitfile,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          allocate(bc(i)%ext_flux(nprim))
          bc(i)%ext_flux = 0.0
          allocate ( bc(i) % Pg (1, 6) )

        ! ─────────────────────────────────────────────────────────────────────
        ! ! MOSE 999: Kaffs 1D–multi-D hybrid coupling (unchanged)
        ! case(999)
        !   if (level == 1) nKaffs = nKaffs + 1
        !   read(unitfile,*,iostat=ios) bc(i)%bs, bc(i)%fs

      end select

    enddo

    if (nwall > 0) obj_io % write_wall = .true.

    close( unitfile )

  end subroutine Read_BCfile


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
  ! parse_dir_tok: convert a direction token (read as character) to real.
  !   If the token contains 'normal', return 1.0e30 (face-normal flag).
  !   Otherwise parse as a floating-point number.
  pure function parse_dir_tok(tok) result(val)
    implicit none
    character(len=*), intent(in) :: tok
    real(R8) :: val
    integer  :: ios_loc
    character(len=len(tok)) :: tok_
    tok_ = adjustl(tok)
    if (index(trim(tok_), 'normal') > 0) then
      val = huge(1.0_R8)
    else
      read(tok_, *, iostat=ios_loc) val
      if (ios_loc /= 0) val = 0.0_R8
    endif
  end function parse_dir_tok

end module MOSE_IO_BC
