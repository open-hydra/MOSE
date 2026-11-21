!> The manifold (BC 501): the faces of a group joined by a lumped plenum. The gas that leaves the domain through the
!> outlet faces of the group (role 1) comes back through its inlet faces (role 2): once per stage, before the boundary
!> fluxes, Manifold_Gather sums over the faces of each group that this rank owns, adds the sums of the other ranks (one
!> reduction per stage), and sets the state the faces of the group take at this stage; BC_Fluxes then applies it with the
!> inlet 403 (mass flux and total temperature, along the face normal, the mean species and turbulence) on the inlet faces
!> and the outlet 406 (a static pressure) on the outlet faces, with the relaxation factor of the record.
!>
!> The means are those of the 2025 manifold of hydra-MF, to reproduce its runs: the area-weighted mean of the conservative
!> state over the outlet faces, from which the mass flux (scaled by the area ratio), the static temperature, the velocity
!> and so the total temperature of the inlet faces follow; the outlet faces take the area-mean static pressure of the
!> inlet faces (K and p0-ratio do not act on this rule). (The mass-flow-weighted means and the equal total pressure of
!> H6_manifold_design.md §7 replace them in the next commit.)
!>
!> The roles are those of the records: no sign of the flow changes them. A cell of an outlet face whose gas moves inward
!> is a slip wall for the stage, as on a 406; when nothing leaves through the outlet faces of a group (the first stages
!> of a run from rest), its inlet faces are slip walls and its outlet faces take the mean static pressure of the inlet
!> side. The record's parameters: K, the loss on the dynamic head of the outlet faces; rel_fac; mdot-ratio, p0-ratio and
!> dh0, the mass flow, total pressure and total enthalpy of the inlet faces against those of the outlet faces (the
!> cycle's topology lumped on one link: one element out of N, a bleed, a turbine).
!>
!> Every rank reads every record, so the groups are built without a reduction (Manifold_Setup, at set-up, per grid
!> level); the faces of a group may lie on any ranks and be cut by the decomposition.
module MOSE_Lib_BC_Fluxes_Manifold
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  private
  public :: Manifold_Setup, Manifold_Gather, MANIFOLD_OUTLET, MANIFOLD_INLET

  integer, parameter :: MANIFOLD_OUTLET = 1   !< role of a face where the gas leaves the domain into the manifold
  integer, parameter :: MANIFOLD_INLET  = 2   !< role of a face where it comes back in

contains

  !─────────────────────────────────────────────────────────────────────────────
  !> The groups of the manifold records of a domain (one grid level), from the records every rank has read: one group
  !> per id, in the order the ids first appear; its parameters from its first record, which every other record of the
  !> group must repeat (ATLAS writes them group-wide); roles 1 or 2 only; at least one record of each role in the whole
  !> domain. Each record gets the index of its group (man_idx). The run stops on every rank on a bad file.
  subroutine Manifold_Setup ( domain, level )
    use MOSE_Mod_MPI, only: mpi_stop_all, mpi_is_root
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer,                intent(in)    :: level
    ! Local
    integer :: i, m, n501, nman
    integer, allocatable :: gids(:)
    character(len=400) :: msg

    domain % nman = 0
    if (allocated(domain % man)) deallocate(domain % man)
    n501 = 0
    do i = 1, domain % nbound
      if (domain % bc(i) % type == 501) n501 = n501 + 1
    end do
    allocate ( domain % man(0) )
    if (n501 == 0) return

    ! The distinct ids, in the order of the file
    allocate ( gids(n501) )
    nman = 0
    do i = 1, domain % nbound
      if (domain % bc(i) % type /= 501) cycle
      if (domain % bc(i) % man_gid < 1) then
        write(msg,'(A,I0,A,I0,A,I0,A)') '[ERROR] manifold (501) record ', i, ' of grid level ', level, &
          ': the group id is ', domain % bc(i) % man_gid, '; the ids start at 1'
        call mpi_stop_all ( trim(msg) )
      end if
      if (domain % bc(i) % man_role /= MANIFOLD_OUTLET .and. domain % bc(i) % man_role /= MANIFOLD_INLET) then
        write(msg,'(A,I0,A,I0,A,I0,A)') '[ERROR] manifold (501) record ', i, ' of grid level ', level, &
          ': the role is ', domain % bc(i) % man_role, '; it is 1 (outlet: the gas leaves the domain) or 2 (inlet)'
        call mpi_stop_all ( trim(msg) )
      end if
      if ( .not. any ( gids(1:nman) == domain % bc(i) % man_gid ) ) then
        nman = nman + 1
        gids(nman) = domain % bc(i) % man_gid
      end if
    end do

    deallocate ( domain % man )
    allocate ( domain % man(nman) )
    domain % nman = nman
    do m = 1, nman
      domain % man(m) % gid = gids(m)
      allocate ( domain % man(m) % ci_in(nsc+nrans) )
      domain % man(m) % ci_in = 0d0
    end do

    ! The records of each group: the parameters of the first, repeated by the others; the roles counted
    do i = 1, domain % nbound
      if (domain % bc(i) % type /= 501) cycle
      do m = 1, nman
        if (gids(m) == domain % bc(i) % man_gid) exit
      end do
      domain % bc(i) % man_idx = m
      associate ( g => domain % man(m), r => domain % bc(i) )
        if (g % n_out + g % n_in == 0) then
          g % K = r % man_K ; g % rf = r % rel_fac ; g % mratio = r % man_mratio ; g % p0ratio = r % man_p0ratio
          g % dh0 = r % man_dh0
          g % first = i
        else
          call Same ( 'K',          r % man_K,       g % K )
          call Same ( 'rel_fac',    r % rel_fac,     g % rf )
          call Same ( 'mdot_ratio', r % man_mratio,  g % mratio )
          call Same ( 'p0_ratio',   r % man_p0ratio, g % p0ratio )
          call Same ( 'dh0',        r % man_dh0,     g % dh0 )
        end if
        if (r % man_role == MANIFOLD_OUTLET) then
          g % n_out = g % n_out + 1
        else
          g % n_in = g % n_in + 1
        end if
      end associate
    end do

    do m = 1, nman
      if (domain % man(m) % n_out == 0 .or. domain % man(m) % n_in == 0) then
        write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] manifold group ', gids(m), ' of grid level ', level, ' has ', &
          domain % man(m) % n_out, ' outlet and ', domain % man(m) % n_in, ' inlet records: a manifold joins at least '// &
          'one face where the gas leaves the domain (role 1) and one where it comes back (role 2)'
        call mpi_stop_all ( trim(msg) )
      end if
    end do

    if (mpi_is_root) then
      write(*,'(A,I0,A)') ' Manifold groups of grid level ', level, ':'
      do m = 1, nman
        associate ( g => domain % man(m) )
          write(*,'(A,I0,A,I0,A,I0,A,ES10.3,A,F6.3,A,ES10.3,A,ES10.3,A,ES12.5,A)') '   group ', g % gid, ': ', g % n_out, &
            ' outlet and ', g % n_in, ' inlet records; K ', g % K, ', rel_fac ', g % rf, ', mdot-ratio ', g % mratio, &
            ', p0-ratio ', g % p0ratio, ', dh0 ', g % dh0, ' J/kg'
        end associate
      end do
    end if

  contains

    subroutine Same ( what, v, v1 )
      character(len=*), intent(in) :: what
      real(R8),         intent(in) :: v, v1
      if ( abs(v - v1) <= 1d-12 * max(1d0, abs(v1)) ) return
      write(msg,'(A,I0,A,I0,A,I0,3A,ES15.8,A,I0,A,ES15.8)') '[ERROR] manifold (501) record ', i, ' of grid level ', &
        level, ' (group ', domain % bc(i) % man_gid, '): its ', trim(what), ' is ', v, ' while record ', &
        domain % man(m) % first, ' of the group has ', v1
      call mpi_stop_all ( trim(msg)//'; the parameters of a manifold are those of its section, the same on every record' )
    end subroutine Same

  end subroutine Manifold_Setup


  !─────────────────────────────────────────────────────────────────────────────
  !> The state the faces of each manifold group take at this stage, from the cells of the group's faces on every rank.
  !> Called by every thread of the parallel region at the top of BC_Fluxes: one thread sums over the records of this rank,
  !> reduces the sums over the ranks (MPI from the primary thread, as the connection exchange) and sets the targets; the
  !> others wait at the barrier.
  !>
  !> Outlet faces, every cell: the area-weighted conservative state (species densities, momentum, total energy; the
  !> turbulence variables likewise) and the area-weighted temperature (the start of the temperature iteration of the
  !> mean state); the cells whose gas leaves also give the mass flow, total enthalpy and dynamic head of the line of the
  !> log. Inlet faces, every cell: the area and the area-weighted static pressure; the dynamic head of the line. Then,
  !> for the group: the mean state of the outlet faces gives the mass flux g_in = |rho u| A_out / A_in times mdot-ratio,
  !> the static temperature T of the mean state and T0_in = T + v^2/(2 cp(T)) of that state, lowered by dh0 through
  !> h(T0); the outlet faces take p_out, the area-mean static pressure of the inlet faces. The species and turbulence
  !> of the inlet faces are those of the mean state.
  subroutine Manifold_Gather ( domain )
    use MOSE_Mod_MPI,        only: mpi_size_, mpi_allreduce_sum_r8_array, mpi_is_root
    use MOSE_Config_Types_m, only: obj_io
    use FLINT_Lib_Thermodynamic
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer  :: nsum, m, ii, i, Bm, Im, Jm, Km, Fm, dir, Face_i, Face_j, Face_k, modfm, modfm1, modfm2, modfm3, s, k
    integer  :: iTA, iAout, iAin, ipA, imdot, igh0, igq, igin, igqin, iturb
    real(R8), allocatable :: sums(:,:), flat(:)
    real(R8) :: A, n(3), Prim(nprim), cons(np), rho, Rgas, Un, g, q, T, hcell, h0cell, cp
    real(R8) :: meancons(np), meanprim(np), v2, A_out, A_in, q_out, q_in, p_in
    real(R8) :: y(nsc)

    if (domain % nman == 0) return

    !$omp master
    ! The slots of the sums of a group: the np conservative variables, then
    iTA   = np + 1    ! sum of T dA over the outlet cells
    iAout = np + 2    ! area of the outlet faces
    iAin  = np + 3    ! area of the inlet faces
    ipA   = np + 4    ! sum of p dA over the inlet cells
    imdot = np + 5    ! mass flow leaving through the outlet cells
    igh0  = np + 6    ! sum of g h0 dA over the leaving outlet cells
    igq   = np + 7    ! sum of g q dA over them, q = v^2/2
    igin  = np + 8    ! sum of g dA over the inlet cells, g = rho |u.n|
    igqin = np + 9    ! sum of g q dA over the inlet cells
    iturb = np + 9    ! then the nrans turbulence variables, rho phi dA
    nsum  = iturb + nrans
    allocate ( sums(nsum, domain % nman) )
    sums = 0d0

    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 501) cycle
      m  = domain % bc(i) % man_idx
      Bm = domain % bc(i) % b ; Im = domain % bc(i) % i ; Jm = domain % bc(i) % j ; Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
      call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
      associate ( blk => domain % blk(Bm) )
        A    = blk % dir(dir) % f(Face_i,Face_j,Face_k) % A
        n    = blk % dir(dir) % f(Face_i,Face_j,Face_k) % n
        Prim = blk % P(:,Im,Jm,Km)
      end associate
      call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
      Un = real(modfm2, R8) * dot_product ( Prim(nu:nw), n )   ! normal velocity out of the domain
      q  = 0.5d0 * sum ( Prim(nu:nw)**2 )
      T  = Prim(np) / ( rho * Rgas )
      select case ( domain % bc(i) % man_role )
      case ( MANIFOLD_OUTLET )
        cons = prim2cons ( Prim(1:np) )
        sums(1:np,m)  = sums(1:np,m)  + cons * A
        sums(iTA,m)   = sums(iTA,m)   + T * A
        sums(iAout,m) = sums(iAout,m) + A
        do k = 1, nrans
          sums(iturb+k,m) = sums(iturb+k,m) + Prim(nt+k-1) * A
        end do
        if (Un > 0d0) then
          g  = rho * Un
          hcell = 0d0
          do s = 1, nsc
            hcell = hcell + Prim(s) / rho * f_tabT ( T, s, h_tab )
          end do
          h0cell = hcell + q
          sums(imdot,m) = sums(imdot,m) + g * A
          sums(igh0,m)  = sums(igh0,m)  + g * h0cell * A
          sums(igq,m)   = sums(igq,m)   + g * q * A
        end if
      case ( MANIFOLD_INLET )
        g = rho * abs(Un)
        sums(iAin,m)  = sums(iAin,m)  + A
        sums(ipA,m)   = sums(ipA,m)   + Prim(np) * A
        sums(igin,m)  = sums(igin,m)  + g * A
        sums(igqin,m) = sums(igqin,m) + g * q * A
      end select
    end do

    if (mpi_size_ > 1) then
      flat = reshape ( sums, [nsum * domain % nman] )
      call mpi_allreduce_sum_r8_array ( flat, nsum * domain % nman )
      sums = reshape ( flat, [nsum, domain % nman] )
    end if

    ! The state of each group
    do m = 1, domain % nman
      associate ( gr => domain % man(m) )
        A_out = sums(iAout,m)
        A_in  = sums(iAin,m)
        p_in  = 0d0
        if (A_in > 0d0) p_in = sums(ipA,m) / A_in
        gr % mdot_out = sums(imdot,m)
        gr % active   = gr % mdot_out > 0d0 .and. A_out > 0d0 .and. A_in > 0d0
        if (gr % active) then
          ! The mean conservative state of the outlet faces and its primitive form (the temperature iteration of
          ! cons2prim starts from the mean temperature)
          meancons = sums(1:np,m) / A_out
          meanprim = cons2prim ( meancons, sums(iTA,m) / A_out )
          call co_rotot_Rtot ( meanprim(1:nsc), rho, Rgas )
          y  = meanprim(1:nsc) / rho
          T  = meanprim(np) / ( rho * Rgas )
          v2 = sum ( meanprim(nu:nw)**2 )
          cp = 0d0
          do s = 1, nsc
            cp = cp + y(s) * f_tabT ( T, s, cp_tab )
          end do
          gr % g_in  = gr % mratio * norm2 ( meancons(nu:nw) ) * A_out / A_in
          gr % T0_in = T + 0.5d0 * v2 / cp
          if (gr % dh0 /= 0d0) gr % T0_in = Total_Temperature ( Mixture_h ( gr % T0_in, y ) - gr % dh0, y, gr % T0_in )
          gr % ci_in(nrans+1:nrans+nsc) = y
          do k = 1, nrans
            gr % ci_in(k) = sums(iturb+k,m) / A_out / rho
          end do
          ! The outlet faces: the area-mean static pressure of the inlet faces, as in 2025 (K and p0-ratio do not act on
          ! this rule); the total pressures of the two sides (mass-weighted dynamic heads q) for the line of the log
          q_out = sums(igq,m) / gr % mdot_out
          q_in  = 0d0
          if (sums(igin,m) > 0d0) q_in = sums(igqin,m) / sums(igin,m)
          gr % p_out  = p_in
          gr % p0_in  = p_in + q_in
          gr % p0_out = p_in + q_out
          ! The line of the log
          gr % mdot_in = gr % g_in * A_in
          gr % h0_out  = sums(igh0,m) / gr % mdot_out
          gr % h0_in   = Mixture_h ( gr % T0_in, y )
        else
          ! Nothing leaves yet: the inlet faces are slip walls, the outlet faces take the pressure of the inlet side
          gr % g_in  = 0d0
          gr % p_out = p_in
          gr % mdot_in = 0d0 ; gr % h0_out = 0d0 ; gr % h0_in = 0d0 ; gr % p0_out = p_in ; gr % p0_in = p_in
        end if
      end associate
    end do

    ! One line per group at the shell interval, once per iteration (the first stage that sees it)
    if ( mpi_is_root .and. obj_io % shell_diter > 0 ) then
      if ( mod ( domain % iter, obj_io % shell_diter ) == 0 ) then
        do m = 1, domain % nman
          associate ( gr => domain % man(m) )
            if ( gr % last_line == domain % iter ) cycle
            gr % last_line = domain % iter
            write(*,'(A,I0,A,I0,A,I0,A,2ES17.9,A,2ES17.9,A,2ES17.9,A,ES11.3,A)') ' manifold ', gr % gid, ', grid level ', &
              domain % mg_level, ', iteration ', domain % iter, ': mdot out/in', gr % mdot_out, gr % mdot_in, &
              ' kg/s, h0 out/in', gr % h0_out, gr % h0_in, ' J/kg, p0 out/in', gr % p0_out, gr % p0_in, ' Pa, taken out', &
              gr % mdot_out * gr % dh0, ' W'
          end associate
        end do
      end if
    end if

    deallocate ( sums )
    !$omp end master
    !$omp barrier

  end subroutine Manifold_Gather


  !─────────────────────────────────────────────────────────────────────────────
  !> The enthalpy of the mixture of mass fractions y at T, from the tables of the gas (f_tabT, h_tab)
  pure real(R8) function Mixture_h ( T, y ) result ( h )
    use FLINT_Lib_Thermodynamic, only: f_tabT, h_tab
    implicit none
    real(R8), intent(in) :: T, y(nsc)
    ! Local
    integer :: s
    h = 0d0
    do s = 1, nsc
      h = h + y(s) * f_tabT ( T, s, h_tab )
    end do
  end function Mixture_h


  !> The temperature at which the mixture of mass fractions y has the enthalpy h0: Newton from T_start with cp of the
  !> tables as the slope (the enthalpy grows with T), to 1e-12 of T_start
  pure real(R8) function Total_Temperature ( h0, y, T_start ) result ( T )
    use FLINT_Lib_Thermodynamic, only: f_tabT, h_tab, cp_tab
    implicit none
    real(R8), intent(in) :: h0, y(nsc), T_start
    ! Local
    integer  :: it, s
    real(R8) :: F, dF, dT
    T = T_start
    do it = 1, 50
      F  = -h0
      dF = 0d0
      do s = 1, nsc
        F  = F  + y(s) * f_tabT ( T, s, h_tab )
        dF = dF + y(s) * f_tabT ( T, s, cp_tab )
      end do
      dT = F / dF
      T  = T - dT
      if ( abs(dT) <= 1d-12 * T_start ) exit
    end do
  end function Total_Temperature

end module MOSE_Lib_BC_Fluxes_Manifold
