module MOSE_Lib_Ghost
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Parameters_m
  use MOSE_Global_m
  use MOSE_Lib_BC_Fluxes, only: Is_Wall

  implicit none

  !> Ghost fills of a wall or 400 face (Ghost_Wall) in which the cut of an extrapolated quantity acted, and in which the
  !> temperature of a 302 wall took the floor T_w/2: reported at the end of a run in which either happened (Report_Wall_Ghosts)
  integer(I8) :: n_ghost_cut = 0, n_ghost_floor = 0

contains

  subroutine Fill_Ghost_Cell ( domain )
    use MOSE_Mod_MPI, only: is_local_block, mpi_size_
    use MOSE_Mod_GhostExchange, only: exchange_ghost_P_post_recv, exchange_ghost_P_pack, &
                                       exchange_ghost_P_post_send, exchange_ghost_P_wait_unpack, &
                                       exchange_ghost_P_wait_send, &
                                       exchange_ghost_chimera_post_recv, exchange_ghost_chimera_pack, &
                                       exchange_ghost_chimera_post_send, exchange_ghost_chimera_wait_recv, &
                                       exchange_ghost_chimera_unpack, exchange_ghost_chimera_wait_send, &
                                       Ghost_Interrank, exchange_ghost_Pg, ghost_sched, &
                                       set_active_mg_level
    use MOSE_Mod_Timers, only: timer_comm_begin, timer_comm_end
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: ii, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    integer :: fg, fdir
    logical :: mpi_on

    ! On a single rank the exchange has nothing to do, and running it anyway
    ! costs a chain of empty master blocks, each one followed by a barrier.  The
    ! test is uniform over the team, so the worksharing constructs and barriers
    ! it guards are still encountered by every thread.
    mpi_on = mpi_size_ > 1

    if (mpi_on) then
      ! MPI: post persistent receives, pack buffer in parallel, then post sends
      ! Chimera donor cells travel in their own non-blocking exchange, started
      ! here so it overlaps with the local BC processing below.
      ! MPI is initialised with MPI_THREAD_FUNNELED: every MPI call is made by
      ! the primary thread, in a master construct (named masked since OpenMP
      ! 5.1; master is kept because older compilers lack masked), followed by
      ! an explicit barrier wherever the other threads need its result. A
      ! single construct would let any thread make the call, which FUNNELED
      ! forbids.
      !$omp master
      call set_active_mg_level(domain%mg_level)
      call exchange_ghost_P_post_recv(domain)
      call exchange_ghost_chimera_post_recv(domain)
      !$omp end master
      !$omp barrier

      ! Pack send buffers in parallel: P over face groups, chimera over donor cells
      !$omp do schedule(static) private(fg)
      do fg = 1, ghost_sched%n_send_faces
        call exchange_ghost_P_pack(domain, fg, fg)
      end do
      !$omp end do nowait
      !$omp do schedule(static) private(ii)
      do ii = 1, ghost_sched%n_chim_send
        call exchange_ghost_chimera_pack(domain, ii, ii)
      end do

      ! Post sends (must wait for all packing to complete — implicit barrier from !$omp do)
      !$omp master
      call exchange_ghost_P_post_send(domain)
      call exchange_ghost_chimera_post_send(domain)
      !$omp end master
    end if

    ! Process LOCAL BC entries while MPI communication is in flight
    ! Uses pre-filtered local_bc_idx to avoid scanning all nbound entries
    !$omp do schedule (dynamic, 64) private(ii, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs)
    do ii = 1, domain % n_local_bc
      i  = domain % local_bc_idx(ii)
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      select case ( domain % bc(i) % type)
        case(101,201) ! block connection
          Bs = domain % bc(i) % bs
          if (.not. is_local_block(Bs)) cycle  ! inter-rank handled after MPI completes
          Is = domain % bc(i) % is
          Js = domain % bc(i) % js
          Ks = domain % bc(i) % ks
          Fs = domain % bc(i) % fs
          call Ghost_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                  Is, Js, Ks, Fs, domain % blk(Bs), &
                                  domain % bc(i) % Pg )
        case(200)
          call Ghost_Axysimmetry ( Im, Jm, Km, Fm, domain % blk(Bm) )
        case(300)
          call Ghost_Symmetry ( Im, Jm, Km, Fm, domain % blk(Bm) )
        case(410)
          call Ghost_Q2D_Mapped ( domain % bc(i), domain % blk(Bm), domain % time )
        case(0,401:407,420) ! inlet/outlet/extrapolation: zero-gradient
          call Ghost_ZG_Extrapolate ( Im, Jm, Km, Fm, domain % blk(Bm) )
        case(103,104,301,302,400,503:505) ! walls (103, 104: coupled to another solver) and extrapolation faces
          call Ghost_Wall ( domain % bc(i), domain % blk(Bm) )
        case(502) ! SRM grain: extrapolation; in a RANS run the turbulence variables from the wall values of the grain
          call Ghost_Extrapolate ( Im, Jm, Km, Fm, domain % blk(Bm) )
          if ( model == 2 ) call Ghost_SRM_Turbulence ( domain % bc(i), domain % blk(Bm) )
        case(102)
          ! Chimera: processed after the donor-cell MPI exchange completes (below)
          continue
        case default
          call Ghost_Extrapolate ( Im, Jm, Km, Fm, domain % blk(Bm) )
      end select
    enddo

    ! Compute Pg for Q2D (410) connections where Bm is local. Two 410 faces of
    ! a block that meet at an edge both write the ghost cell of the edge and
    ! read it back: one loop per direction of the faces, each ending in its
    ! implicit barrier.
    do fdir = 1, 3
      !$omp do schedule (dynamic, 64) private(ii, i, Bm, Im, Jm, Km, Fm)
      do ii = 1, domain % n_local_bc
        i = domain % local_bc_idx(ii)
        Bm = domain % bc(i) % b
        Im = domain % bc(i) % i
        Jm = domain % bc(i) % j
        Km = domain % bc(i) % k
        Fm = domain % bc(i) % f
        if ( (Fm + 1) / 2 /= fdir ) cycle
        select case (domain % bc(i) % type)
          case(410)
            call Fill_BC_Ghost_Q2D ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Pg )
        end select
      enddo
    enddo

    ! MPI: wait for P receives to complete; unpack chimera donor cells into
    ! the kept-alive P arrays of remote donor blocks.  What is timed here is
    ! the communication the local BC work did not hide, plus the wait on
    ! slower neighbours.
    if (mpi_on) then
      !$omp master
      call timer_comm_begin()
      call exchange_ghost_P_wait_unpack(domain)
      call exchange_ghost_chimera_wait_recv(domain)
      call timer_comm_end()
      !$omp end master
      !$omp barrier

      ! Chimera donor cells are distinct, so the unpack has no write conflicts
      !$omp do schedule(static) private(ii)
      do ii = 1, ghost_sched%n_chim_recv
        call exchange_ghost_chimera_unpack(domain, ii, ii)
      end do
    end if

    ! Compute Pg(:,3:6) for type-1 connections where source block is local.
    ! The four tangential neighbours of a source cell on the edge of its face are
    ! ghost cells of Bs, and those can be filled from another rank: this loop must
    ! follow the P unpack above, or it reads them one stage old (MPI-only drift at
    ! block edges and junctions, ~1e-9 after 500 iterations on RL36).
    !$omp do schedule (dynamic, 64) private(ii, i, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do ii = 1, domain % n_local_bs
      i  = domain % local_bs_idx(ii)
      Bs = domain % bc(i) % bs
      Is = domain % bc(i) % is
      Js = domain % bc(i) % js
      Ks = domain % bc(i) % ks
      Fs = domain % bc(i) % fs
      d11s = domain % bc(i) % d11
      d12s = domain % bc(i) % d12
      d21s = domain % bc(i) % d21
      d22s = domain % bc(i) % d22
      call Fill_BC_Ghost_Connection ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, &
                                      domain % blk(Bs), domain % bc(i) % Pg )
    enddo

    ! Chimera ghost fill: all donor data (local and remote) is now current
    !$omp do schedule (dynamic, 64) private(ii, i)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 102) cycle
      call Ghost_Chimera ( domain % nb, domain % blk, domain % bc(i) )
    enddo

    ! Chimera Pg(:,3:6): reads the ghost columns written just above, so it
    ! must stay a separate loop (implicit barrier in between). Each entry
    ! writes its own Pg only: beyond the edges of a block face the ghost row
    ! is extrapolated into Pg, not written in the block (Fill_BC_Ghost_Chimera),
    ! so the entries of two chimera faces of a block that meet at an edge write
    ! no common location, and one loop over all the entries has no data race
    !$omp do schedule (dynamic, 64) private(ii, i, Bm, Im, Jm, Km, Fm)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 102) cycle
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      call Fill_BC_Ghost_Chimera ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Pg )
    enddo

    ! Process INTER-RANK type-1 entries (Bm local, Bs remote)
    if (mpi_on) then
      !$omp do schedule (dynamic, 64) private(ii, i, Bm, Im, Jm, Km, Fm, Bs)
      do ii = 1, domain % n_local_bc
        i  = domain % local_bc_idx(ii)
        if (domain%bc(i)%type /= 101 .and. domain%bc(i)%type /= 201) cycle
        Bs = domain % bc(i) % bs
        if (is_local_block(Bs)) cycle  ! already processed above
        Bm = domain % bc(i) % b
        Im = domain % bc(i) % i
        Jm = domain % bc(i) % j
        Km = domain % bc(i) % k
        Fm = domain % bc(i) % f
        call Ghost_Interrank(Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Pg)
      enddo

      ! Wait for P sends to complete before reusing buffers
      !$omp master
      call exchange_ghost_P_wait_send(domain)
      call exchange_ghost_chimera_wait_send(domain)
      call exchange_ghost_Pg(domain)
      !$omp end master
      !$omp barrier
    end if

  end subroutine Fill_Ghost_Cell


  subroutine Ghost_Connection ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, Pg )

    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(MOSE_block_type), intent(in) :: blks
    type(MOSE_block_type), intent(inout) :: blkm
    real(R8), intent(inout)              :: Pg(nprim,6)
    ! Local
    integer :: g, Ig, Jg, Kg

    do g = 1, gc
      Ig = Im - guide(Fm,1)*(g)
      Jg = Jm - guide(Fm,2)*(g)
      Kg = Km - guide(Fm,3)*(g)
      blkm % P(:,Ig,Jg,Kg) = blks % P (:, Is + guide(Fs,1)*(g-1), &
                                          Js + guide(Fs,2)*(g-1), &
                                          Ks + guide(Fs,3)*(g-1)  )
      ! reference velocity of the preconditioner, read by the connection flux
      blkm % Ur(Ig,Jg,Kg) = blks % Ur (Is + guide(Fs,1)*(g-1), &
                                       Js + guide(Fs,2)*(g-1), &
                                       Ks + guide(Fs,3)*(g-1)  )
      Pg (:,g) = blkm % P (:,Ig,Jg,Kg)
    enddo

  end subroutine Ghost_Connection


  subroutine Ghost_Axysimmetry ( Im, Jm, Km, Fm, blk )
    use MOSE_Lib_Metrics, only: cos_delthe, sin_delthe
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: blk
    ! Local
    integer :: Ig, Jg, Kg, Ks

    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    if (ndir < 3) then 
      ! Rotate v and w of +/- delthe and copy other variables
      blk % P (nv,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Km) * cos_delthe + guide(Fm,3) * blk % P (nw,Im,Jm,Km) * sin_delthe
      blk % P (nw,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Km) * sin_delthe * (-guide(Fm,3)) + blk % P (nw,Im,Jm,Km) * cos_delthe
      blk % P (1:nu,Ig,Jg,Kg) = blk % P (1:nu,Im,Jm,Km)
      blk % P (np:nprim,Ig,Jg,Kg) = blk % P (np:nprim,Im,Jm,Km)
    else
      ! The cell next to the opposite k face, turned about x by -guide(Fm,3)*delthe: back by the
      ! sector angle behind face 5, forward behind face 6 (the same turn as the 2-D branch above)
      Ks = blk % dim(3) - Km + 1
      blk % P (nv,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Ks) * cos_delthe + guide(Fm,3) * blk % P (nw,Im,Jm,Ks) * sin_delthe
      blk % P (nw,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Ks) * sin_delthe * (-guide(Fm,3)) + blk % P (nw,Im,Jm,Ks) * cos_delthe
      blk % P (1:nu,Ig,Jg,Kg) = blk % P (1:nu,Im,Jm,Ks)
      blk % P (np:nprim,Ig,Jg,Kg) = blk % P (np:nprim,Im,Jm,Ks)
    endif
      
  end subroutine Ghost_Axysimmetry


  subroutine Ghost_Symmetry ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: blk
    ! Local
    integer  :: Ig, Jg, Kg, Ic, Jc, Kc
    real(R8) :: normal(3), vel(3), udotn

    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    if (Fm <= 2) then
      Ic = Im - mod(Fm,2)
      normal = blk % dir(1) % f (Ic,Jm,Km) % N
    elseif (Fm <= 4) then
      Jc = Jm - mod(Fm,2)
      normal = blk % dir(2) % f (Im,Jc,Km) % N
    else
      Kc = Km - mod(Fm,2)
      normal = blk % dir(3) % f (Im,Jm,Kc) % N
    endif

    ! Formula to invert normal velocity: u = |u_n|n + |u_t|t -2*|u_n|n = -|u_n|n + |u_t|t
    vel = blk % P (nu:nw,Im,Jm,Km)   ! cell velocity
    udotn = dot_product( vel, normal )  ! normal component
    vel = vel - 2d0*udotn*normal         ! ghost cell velocity
    
    blk % P (nu:nw,Ig,Jg,Kg) = vel
    blk % P (1:nsc,Ig,Jg,Kg) = blk % P (1:nsc,Im,Jm,Km)
    blk % P (np:nprim,Ig,Jg,Kg) = blk % P (np:nprim,Im,Jm,Km)
      
  end subroutine Ghost_Symmetry

  !> Ghost cells of a mapped face (410), both layers, and Pg(:,1:2): the state of the cell across the face that its file
  !> gives, the densities of the species, the velocity and the pressure, interpolated at the time; the quantities after p
  !> (the turbulence variables and the soot), which the file does not give, those of the cell per unit mass, phi_g = phi_1,
  !> that is (rho phi)_g = rho_g (rho phi)_1/rho_1 with rho the sum of the densities of the species. So the specific values
  !> whose differences the viscous fluxes and the sources take, phi = (rho phi)/rho (Visc_State), are those of the cell,
  !> and a uniform field of them stays uniform across the face, whatever the density of the file.
  subroutine Ghost_Q2D_Mapped ( bc, blk, time )
    implicit none
    type(MOSE_bc_type), intent(inout)    :: bc
    type(MOSE_block_type), intent(inout) :: blk
    real(R8), intent(in)                 :: time
    ! Local
    integer  :: Im, Jm, Km, Fm, Ig, Jg, Kg, Ig2, Jg2, Kg2, s
    real(R8) :: Pinterp(nprim)

    Im = bc % i
    Jm = bc % j
    Km = bc % k
    Fm = bc % f

    ! Interpolate Q2D data at current time (per-variable, reusing time_bc_type): the file gives every primitive up to the
    ! pressure (Read_Q2D_Tecplot stops on a file that lacks one; a species that it does not name has zero density)
    do s = 1, np
      if (bc % q2d_map(s) % exists) then
        if (bc % q2d_map(s) % periodic) then
          Pinterp(s) = bc % q2d_map(s) % update(time)
        else
          if (time <= bc % q2d_map(s) % time(1)) then
            Pinterp(s) = bc % q2d_map(s) % var(1)
          elseif (time >= bc % q2d_map(s) % time(bc % q2d_map(s) % n)) then
            Pinterp(s) = bc % q2d_map(s) % var(bc % q2d_map(s) % n)
          else
            Pinterp(s) = bc % q2d_map(s) % update(time)
          endif
        endif
      else
        Pinterp(s) = blk % P(s, Im, Jm, Km)
      endif
    enddo

    ! The quantities after p (the turbulence variables and the soot): the values of the cell per unit mass, at the density
    ! of the file
    if ( nprim > np ) Pinterp(np+1:nprim) = sum ( Pinterp(1:nsc) ) * ( blk % P(np+1:nprim,Im,Jm,Km) &
                                                                       / sum ( blk % P(1:nsc,Im,Jm,Km) ) )

    ! First ghost cell
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)
    blk % P(:, Ig, Jg, Kg) = Pinterp

    ! Second ghost cell: same value (first order)
    Ig2 = Im - guide(Fm,1)*2
    Jg2 = Jm - guide(Fm,2)*2
    Kg2 = Km - guide(Fm,3)*2
    blk % P(:, Ig2, Jg2, Kg2) = Pinterp

    ! Store ghost primitives in Pg for viscous flux
    bc % Pg(:,1) = Pinterp
    bc % Pg(:,2) = Pinterp

  end subroutine Ghost_Q2D_Mapped


  subroutine Ghost_ZG_Extrapolate ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: blk
    ! Local
    integer :: g

    ! Simple <zero gradient> extrapolation
    do g = 1, gc
      select case (Fm)
        case(1)
          blk % P(:,Im-g,Jm,Km) = blk % P (:,Im,Jm,Km)
        case(2)
          blk % P(:,Im+g,Jm,Km) = blk % P (:,Im,Jm,Km)
        case(3)
          blk % P(:,Im,Jm-g,Km) = blk % P (:,Im,Jm,Km)
        case(4)
          blk % P(:,Im,Jm+g,Km) = blk % P (:,Im,Jm,Km)
        case(5)
          blk % P(:,Im,Jm,Km-g) = blk % P (:,Im,Jm,Km)
        case(6)
          blk % P(:,Im,Jm,Km+g) = blk % P (:,Im,Jm,Km)
      end select
    enddo
      
  end subroutine Ghost_ZG_Extrapolate


  subroutine Ghost_Chimera ( nb, Blk, bc )
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in)                  :: nb
    type(MOSE_block_type), intent(inout) :: Blk(nb)
    type(MOSE_bc_type), intent(inout)    :: bc
    ! Local
    integer  :: Bm, Im, Jm, Km, Fm, Ig, Jg, Kg, Ig2, Jg2, Kg2, Bs, Is, Js, Ks, c 
    real(R8) :: rho, Rgas, temperature, consi(nprim), consg(nprim)

    ! Preliminary assignments
    Bm = bc % b
    Im = bc % i
    Jm = bc % j
    Km = bc % k
    Fm = bc % f

    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)
    Ig2 = Im - guide(Fm,1)*2
    Jg2 = Jm - guide(Fm,2)*2
    Kg2 = Km - guide(Fm,3)*2
          
    ! RIEMPIMPENTO GHOST CON VARIABILI CONSERVATE
    call co_rotot_Rtot ( blk(Bm) % P (1:nsc,im,jm,km), rho, Rgas )
    temperature = EOS( blk(Bm) % P (np,im,jm,km), rho, Rgas)
    ! First row of ghost cell coordinates
    consg = 0.d0
    do c = 1, bc % ni(1)
      Bs = bc % donorID(c,1)
      Is = bc % donorID(c,2)
      Js = bc % donorID(c,3)
      Ks = bc % donorID(c,4)
      consi(1:np)       = prim2cons ( blk(Bs) % P (1:np,Is,Js,Ks) )   ! flow: conservative blend
      consi(np+1:nprim) = blk(Bs) % P (np+1:nprim,Is,Js,Ks) &        ! soot/passive/RANS: density-weighted blend
                        * sum ( blk(Bs) % P (1:nsc,Is,Js,Ks) )
      consg = consg + consi * bc % volume_fraction(c)
    enddo
    blk(Bm) % P (1:np,Ig,Jg,Kg)       = cons2prim ( consg(1:np), temperature )
    blk(Bm) % P (np+1:nprim,Ig,Jg,Kg) = consg(np+1:nprim) / sum ( consg(1:nsc) )
    bc % Pg (:,1) = blk(Bm) % P (:, Ig,Jg,Kg)

    ! Second row of ghost cell coordinates
    consg = 0.d0
    do c = bc % ni(1)+1, sum ( bc % ni )
      Bs = bc % donorID(c,1)
      Is = bc % donorID(c,2)
      Js = bc % donorID(c,3)
      Ks = bc % donorID(c,4)
      consi(1:np)       = prim2cons ( blk(Bs) % P (1:np,Is,Js,Ks) )   ! flow: conservative blend
      consi(np+1:nprim) = blk(Bs) % P (np+1:nprim,Is,Js,Ks) &        ! soot/passive/RANS: density-weighted blend
                        * sum ( blk(Bs) % P (1:nsc,Is,Js,Ks) )
      consg = consg + consi * bc % volume_fraction(c)
    enddo
    blk(Bm) % P (1:np,Ig2,Jg2,Kg2)       = cons2prim ( consg(1:np), temperature )
    blk(Bm) % P (np+1:nprim,Ig2,Jg2,Kg2) = consg(np+1:nprim) / sum ( consg(1:nsc) )
    bc % Pg (:,2) = blk(Bm) % P (:,Ig2,Jg2,Kg2)
    
  end subroutine Ghost_Chimera


  subroutine Ghost_Extrapolate ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: blk
    ! Local
    integer :: g

    ! Extrapolation with 2nd order accuracy. Theory: suppose we have a stencil x1,x2,x3,x4 with x4 unknown.
    ! Taylor 2nd order: x4 = x3 + f1(x3) + 1/2 f2(x3) with f1 and f2 1st and 2nd derivative. Approximate 
    ! f1(x3) = (x4 - x2)/2 central and f2(x3) = x3 - 2x2 + x1 left sided => x4 = 3*x3 - 3*x2 + x1
    do g = 1, gc
      select case (Fm)
        case(1)
          blk % P(:,Im-g,Jm,Km) = 3d0 * blk % P (:,Im-g+1,Jm,Km) - &
                                  3d0 * blk % P (:,Im-g+2,Jm,Km) + &
                                        blk % P (:,Im-g+3,Jm,Km)
        case(2)
          blk % P(:,Im+g,Jm,Km) = 3d0 * blk % P (:,Im+g-1,Jm,Km) - &
                                  3d0 * blk % P (:,Im+g-2,Jm,Km) + &
                                        blk % P (:,Im+g-3,Jm,Km)
        case(3)
          blk % P(:,Im,Jm-g,Km) = 3d0 * blk % P (:,Im,Jm-g+1,Km) - &
                                  3d0 * blk % P (:,Im,Jm-g+2,Km) + &
                                        blk % P (:,Im,Jm-g+3,Km)
        case(4)
          blk % P(:,Im,Jm+g,Km) = 3d0 * blk % P (:,Im,Jm+g-1,Km) - &
                                  3d0 * blk % P (:,Im,Jm+g-2,Km) + &
                                        blk % P (:,Im,Jm+g-3,Km)
        case(5)
          blk % P(:,Im,Jm,Km-g) = 3d0 * blk % P (:,Im,Jm,Km-g+1) - &
                                  3d0 * blk % P (:,Im,Jm,Km-g+2) + &
                                        blk % P (:,Im,Jm,Km-g+3)
        case(6)
          blk % P(:,Im,Jm,Km+g) = 3d0 * blk % P (:,Im,Jm,Km+g-1) - &
                                  3d0 * blk % P (:,Im,Jm,Km+g-2) + &
                                        blk % P (:,Im,Jm,Km+g-3)
      end select
    enddo
      
  end subroutine Ghost_Extrapolate


  subroutine Fill_BC_Ghost_Connection ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, blks, Pg )
    implicit none
    integer, intent(in) :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(MOSE_block_type), intent(in) :: blks
    real(R8), intent(inout)           :: Pg(nprim,6)
    ! Local
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4

    select case(Fs)
      case(1:2)
        i1 = Is
        j1 = Js - d11s
        k1 = Ks - d12s
        i2 = Is
        j2 = Js + d11s
        k2 = Ks + d12s
        i3 = Is
        j3 = Js - d21s
        k3 = Ks - d22s
        i4 = Is
        j4 = Js + d21s
        k4 = Ks + d22s
      case(3:4)
        i1 = Is - d11s
        j1 = Js 
        k1 = Ks - d12s
        i2 = Is + d11s
        j2 = Js 
        k2 = Ks + d12s
        i3 = Is - d21s
        j3 = Js
        k3 = Ks - d22s
        i4 = Is + d21s
        j4 = Js 
        k4 = Ks + d22s
      case(5:6)
        i1 = Is - d11s
        j1 = Js - d12s
        k1 = Ks 
        i2 = Is + d11s
        j2 = Js + d12s
        k2 = Ks
        i3 = Is - d21s
        j3 = Js - d22s
        k3 = Ks
        i4 = Is + d21s
        j4 = Js + d22s
        k4 = Ks
    end select

    Pg (:,3) = blks % P (:,i1,j1,k1)
    Pg (:,4) = blks % P (:,i2,j2,k2)
    Pg (:,5) = blks % P (:,i3,j3,k3)
    Pg (:,6) = blks % P (:,i4,j4,k4)

  end subroutine Fill_BC_Ghost_Connection


  !> Pg(:,3:6) of a chimera face: the tangential neighbours of its ghost cell, for the viscous stencil of the face
  !> (BC_Connection_Visc), which takes on the side of the ghost cells the central differences along the face of the row of
  !> the ghost cells, as on the side of the cell those of the row of the cell. Pg(:,3:4) are the lower and the upper
  !> neighbour in the first tangential direction of the face, Pg(:,5:6) in the second (j and k for faces 1-2, i and k for
  !> 3-4, i and j for 5-6). Within the block face a neighbour is the ghost cell of the row, filled from the donors
  !> (Ghost_Chimera). Beyond an edge of the block face the row has no ghost cell: the neighbour is extrapolated linearly
  !> along the row, in the variables whose differences the stencil takes, so that the central difference there is the
  !> one-sided one, exact for a linear field as that of the row of the cell is (the ghost cell copied there halves it, and
  !> the face takes 3/4 of the derivative of a linear field along that direction). A block one cell wide in a direction has
  !> no row to extrapolate along: in the k direction of a 2-D mesh the neighbour is the ghost cell turned by the wedge
  !> angle, as the k neighbours of a cell are (Ghost_Axysimmetry; no turn on a plane mesh), elsewhere the ghost cell
  !> itself. Nothing is written in the block: two chimera faces that meet at an edge of a block share no written location.
  subroutine Fill_BC_Ghost_Chimera ( Im, Jm, Km, Fm, blk, Pg )
    use MOSE_Lib_Metrics, only: delthe
    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(in) :: blk
    real(R8), intent(inout)           :: Pg(nprim,6)
    ! Local
    integer :: c(3), g(3), e(3), t(2), n, s, m

    c = [ Im, Jm, Km ]
    g = c - guide(Fm,:)
    select case (Fm)
      case (1:2)
        t = [ 2, 3 ]
      case (3:4)
        t = [ 1, 3 ]
      case default
        t = [ 1, 2 ]
    end select

    do n = 1, 2
      e = 0
      e(t(n)) = 1
      do s = -1, 1, 2
        m = 2*n + 1 + (s+1)/2   ! Pg(:,3) and Pg(:,4) for n = 1, Pg(:,5) and Pg(:,6) for n = 2
        if ( c(t(n)) + s >= 1 .and. c(t(n)) + s <= blk % dim(t(n)) ) then
          ! Within the block face: the ghost cell of the row
          Pg(:,m) = blk % P (:,g(1)+s*e(1),g(2)+s*e(2),g(3)+s*e(3))
        elseif ( blk % dim(t(n)) > 1 ) then
          ! Beyond an edge of the block face: extrapolated along the row
          Pg(:,m) = Visc_Extrapolated ( blk % P (:,g(1),g(2),g(3)), blk % P (:,g(1)-s*e(1),g(2)-s*e(2),g(3)-s*e(3)) )
        else
          ! One cell in this direction: the ghost cell, in the k direction of a 2-D mesh turned by the wedge angle as by
          ! Ghost_Axysimmetry beyond face 5 (s = -1) or face 6 (s = 1)
          Pg(:,m) = blk % P (:,g(1),g(2),g(3))
          if ( t(n) == 3 .and. ndir < 3 ) then
            Pg(nv,m) = blk % P (nv,g(1),g(2),g(3)) * cos(delthe) - s * blk % P (nw,g(1),g(2),g(3)) * sin(delthe)
            Pg(nw,m) = blk % P (nv,g(1),g(2),g(3)) * sin(delthe) * s + blk % P (nw,g(1),g(2),g(3)) * cos(delthe)
          end if
        end if
      enddo
    enddo

  contains

    !> The state whose viscous variables (mass fractions, velocity, temperature, and the other variables per unit mass:
    !> those of Visc_Variables in BC_Connection_Visc) are twice those of the state a minus those of the state b, with the
    !> density of a, which keeps it positive whatever the jump from b to a
    function Visc_Extrapolated ( a, b ) result ( q )
      use FLINT_Lib_Thermodynamic
      implicit none
      real(R8), intent(in) :: a(nprim), b(nprim)
      real(R8)             :: q(nprim)
      ! Local
      real(R8) :: rho_a, rho_b, rho_q, R_a, R_b, R_q, temperature

      call co_rotot_Rtot ( a(1:nsc), rho_a, R_a )
      call co_rotot_Rtot ( b(1:nsc), rho_b, R_b )
      q(1:nsc)    = ( 2d0 * a(1:nsc) / rho_a - b(1:nsc) / rho_b ) * rho_a
      q(nu:nw)    = 2d0 * a(nu:nw) - b(nu:nw)
      temperature = 2d0 * a(np) / ( rho_a * R_a ) - b(np) / ( rho_b * R_b )
      call co_rotot_Rtot ( q(1:nsc), rho_q, R_q )
      q(np)       = rho_q * R_q * temperature
      if ( nprim > np ) q(np+1:nprim) = ( 2d0 * a(np+1:nprim) / rho_a - b(np+1:nprim) / rho_b ) * rho_a

    end function Visc_Extrapolated

  end subroutine Fill_BC_Ghost_Chimera


  subroutine Fill_BC_Ghost_Q2D ( Im, Jm, Km, Fm, blk, Pg )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: blk
    real(R8), intent(inout)              :: Pg(nprim,6)
    ! Local
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4
    integer :: dim(3), Ig, Jg, Kg

    dim = blk % dim
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    ! Edge handling: extrapolate ghost value to adjacent ghost positions at face boundaries
    select case (Fm)
      case(1:2)
        if ( Jm == 1 )      blk % P (:,Ig, Jm-1, Km) = blk % P (:,Ig, Jm, Km)
        if ( Jm == dim(2) ) blk % P (:,Ig, Jm+1, Km) = blk % P (:,Ig, Jm, Km)
        if ( Km == 1 )      blk % P (:,Ig, Jm, Km-1) = blk % P (:,Ig, Jm, Km)
        if ( Km == dim(3) ) blk % P (:,Ig, Jm, Km+1) = blk % P (:,Ig, Jm, Km)
      case(3:4)
        if ( Im == 1 )      blk % P (:,Im-1, Jg, Km) = blk % P (:,Im, Jg, Km)
        if ( Im == dim(1) ) blk % P (:,Im+1, Jg, Km) = blk % P (:,Im, Jg, Km)
        if ( Km == 1 )      blk % P (:,Im, Jg, Km-1) = blk % P (:,Im, Jg, Km)
        if ( Km == dim(3) ) blk % P (:,Im, Jg, Km+1) = blk % P (:,Im, Jg, Km)
      case(5:6)
        if ( Im == 1 )      blk % P (:,Im-1, Jm, Kg) = blk % P (:,Im, Jm, Kg)
        if ( Im == dim(1) ) blk % P (:,Im+1, Jm, Kg) = blk % P (:,Im, Jm, Kg)
        if ( Jm == 1 )      blk % P (:,Im, Jm-1, Kg) = blk % P (:,Im, Jm, Kg)
        if ( Jm == dim(2) ) blk % P (:,Im, Jm+1, Kg) = blk % P (:,Im, Jm, Kg)
    end select

    ! Fill Pg(:,3:6) from ghost layer (tangential neighbors already populated by Q2D)
    select case(Fm)
      case(1:2)
        i1 = Ig;  j1 = Jm - 1;  k1 = Km
        i2 = Ig;  j2 = Jm + 1;  k2 = Km
        i3 = Ig;  j3 = Jm;      k3 = Km - 1
        i4 = Ig;  j4 = Jm;      k4 = Km + 1
      case(3:4)
        i1 = Im - 1;  j1 = Jg;  k1 = Km
        i2 = Im + 1;  j2 = Jg;  k2 = Km
        i3 = Im;      j3 = Jg;  k3 = Km - 1
        i4 = Im;      j4 = Jg;  k4 = Km + 1
      case(5:6)
        i1 = Im - 1;  j1 = Jm;      k1 = Kg
        i2 = Im + 1;  j2 = Jm;      k2 = Kg
        i3 = Im;      j3 = Jm - 1;  k3 = Kg
        i4 = Im;      j4 = Jm + 1;  k4 = Kg
    end select

    Pg (:,3) = blk % P (:,i1,j1,k1)
    Pg (:,4) = blk % P (:,i2,j2,k2)
    Pg (:,5) = blk % P (:,i3,j3,k3)
    Pg (:,6) = blk % P (:,i4,j4,k4)

  end subroutine Fill_BC_Ghost_Q2D


  !> Ghost cells of a wall (103, 104, 301, 302, 503-505) or of an extrapolation face (400), both layers
  !> (docs/theory/numerics.md, "Ghost cells of walls and 400 faces"). The quantities whose differences the viscous fluxes
  !> and the turbulence sources take, Y_s, T, u and the specific turbulence and scalar variables phi = P/rho, are
  !> extrapolated, and the pressure; the ghost state follows from them, rho = p / (R(Y) T). A quantity with a wall value
  !> q_w takes 2 q_w - q_1: the velocity of a 301 or 302 wall (zero, or the wall velocity of a stationary face of a
  !> rotating frame) and the turbulence variables of the walls (their wall values, as the wall routines take them, so that
  !> the ghost values written here are the only ones: the fluxes and the source terms read the same); the temperature of a
  !> 302 wall, max(2 T_w - T_1, T_w/2). The others q_1 + psi (q_1 - q_2), psi <= s = l0/l1 (rec%sx), cut where needed to keep
  !> each species, the pressure, the temperature and each scalar within [q_1/2, 2 q_1] (Bounded_Extrapolation; the
  !> composition is normalised after the cut), the velocity not cut. Exact for the linear fields whose wall values are
  !> their values at the centre of the face, at the positions of MOSE_Lib_Stencil_Metrics (2 x_f - x_1 for the quantities
  !> with a wall value, x_1 + s (x_1 - x_2) for the extrapolated ones). In an inviscid run the walls take no wall values:
  !> every quantity is extrapolated, as behind a 400 face. With one cell along the direction of the face the extrapolated
  !> quantities are copied.
  subroutine Ghost_Wall ( rec, blk )
    use MOSE_Config_Types_m, only: obj_rot
    use MOSE_Lib_RotatingFrame, only: RF_Is_Stationary_Face, RF_Face_Center, RF_Wall_Velocity
    use MOSE_Lib_RANS, only: RANS_Set_Wall_Values
    use MOSE_Lib_BC_Fluxes, only: Face_Index, Wall_State_Heat, Wall_State_Temperature
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, f_Rtot
    implicit none
    type(MOSE_bc_type), intent(in)       :: rec
    type(MOSE_block_type), intent(inout) :: blk
    ! Local
    integer  :: Im, Jm, Km, Fm, c(3,2), l, s, Dir, Face_i, Face_j, Face_k
    real(R8) :: rho(2), Rgas(2), T(2), Y(nsc,2), u(3,2), phi(nprim-np,2)
    real(R8) :: YG(nsc,2), TG(2), pG(2), uG(3,2), phiG(nprim-np,2), uw(3), phiw(nprim-np), rw(nrans), r_fc(3)
    real(R8) :: rho_w, rhos_w(nsc), mil, kl, temp, dT, Twall, kr, sx, rhoG
    logical  :: wall, cut, cutq

    Im = rec % i ; Jm = rec % j ; Km = rec % k ; Fm = rec % f
    sx = rec % sx
    wall = model > 0 .and. Is_Wall(rec % type)

    ! The boundary cell (1) and the next one along the direction of the face (2; the cell itself in a block of one cell)
    c(:,1) = [Im, Jm, Km]
    c(:,2) = c(:,1)
    if ( blk % dim((Fm+1)/2) > 1 ) c(:,2) = c(:,1) + guide(Fm,:)
    do l = 1, 2
      associate ( P => blk % P(:,c(1,l),c(2,l),c(3,l)) )
      call co_rotot_Rtot ( P(1:nsc), rho(l), Rgas(l) )
      T(l) = P(np) / ( rho(l) * Rgas(l) )
      Y(:,l) = P(1:nsc) / rho(l)
      u(:,l) = P(nu:nw)
      phi(:,l) = P(np+1:nprim) / rho(l)
      end associate
    enddo

    cut = .false.

    ! Species (each cut, then the composition normalised) and pressure
    do s = 1, nsc
      call Bounded_Extrapolation ( Y(s,1), Y(s,2), sx, YG(s,1), YG(s,2), cutq )
      cut = cut .or. cutq
    enddo
    YG(:,1) = YG(:,1) / sum ( YG(:,1) )
    YG(:,2) = YG(:,2) / sum ( YG(:,2) )
    call Bounded_Extrapolation ( blk % P(np,c(1,1),c(2,1),c(3,1)), blk % P(np,c(1,2),c(2,2),c(3,2)), sx, pG(1), pG(2), cutq )
    cut = cut .or. cutq

    ! Temperature
    if ( wall .and. rec % type == 302 ) then
      TG(1) = max ( 2d0 * rec % Tw - T(1), 0.5d0 * rec % Tw )
      TG(2) = max ( 2d0 * rec % Tw - T(2), 0.5d0 * rec % Tw )
      if ( 2d0 * rec % Tw - T(1) < 0.5d0 * rec % Tw ) then
        !$omp atomic
        n_ghost_floor = n_ghost_floor + 1
      end if
    else
      call Bounded_Extrapolation ( T(1), T(2), sx, TG(1), TG(2), cutq )
      cut = cut .or. cutq
    end if

    ! Velocity
    if ( wall .and. rec % type <= 302 ) then
      uw = 0d0
      if ( obj_rot % enabled ) then
        if ( RF_Is_Stationary_Face(rec % b, Fm) ) then
          call RF_Face_Center ( blk % node, Im, Jm, Km, Fm, r_fc )
          uw = RF_Wall_Velocity ( r_fc )
        end if
      end if
      uG(:,1) = 2d0 * uw - u(:,1)
      uG(:,2) = 2d0 * uw - u(:,2)
    else
      uG(:,1) = u(:,1) + sx * ( u(:,1) - u(:,2) )
      uG(:,2) = u(:,1) + 2d0 * sx * ( u(:,1) - u(:,2) )
    end if

    ! Turbulence and scalar variables: the wall values of the turbulence variables, from the wall state of the wall
    ! routines (302: BC_Wall_Temperature; 301: BC_Wall_Heat; 503-505, whose surface model sets the wall state at every
    ! call, and the coupled walls 103 and 104: the state of the cell, as an adiabatic wall)
    if ( wall .and. model == 2 ) then
      kr = 0d0
      select case ( rec % type )
        case (302)
          kr = rec % k_rough
          call Wall_State_Temperature ( blk % P(:,Im,Jm,Km), rho(1), Rgas(1), rec % Tw, rho_w, rhos_w, mil, kl )
        case (301)
          kr = rec % k_rough
          call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )
          call Wall_State_Heat ( blk % P(:,Im,Jm,Km), rho(1), Rgas(1), rec % qw, blk % M(Im,Jm,Km) % c, &
                                 blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n, Dir, temp, mil, kl, dT, Twall, rho_w )
        case default
          call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )
          call Wall_State_Heat ( blk % P(:,Im,Jm,Km), rho(1), Rgas(1), 0d0, blk % M(Im,Jm,Km) % c, &
                                 blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n, Dir, temp, mil, kl, dT, Twall, rho_w )
      end select
      call RANS_Set_Wall_Values ( mil, blk % P(nt:nprim,Im,Jm,Km) * rho_w / rho(1), rw, blk % yn(Im,Jm,Km), kr )
      phiw = 0d0
      phiw(nt-np:) = rw / rho_w
    end if
    do s = 1, nprim - np
      if ( wall .and. model == 2 .and. np + s >= nt ) then
        phiG(s,1) = 2d0 * phiw(s) - phi(s,1)
        phiG(s,2) = 2d0 * phiw(s) - phi(s,2)
      else
        call Bounded_Extrapolation ( phi(s,1), phi(s,2), sx, phiG(s,1), phiG(s,2), cutq )
        cut = cut .or. cutq
      end if
    enddo

    if ( cut ) then
      !$omp atomic
      n_ghost_cut = n_ghost_cut + 1
    end if

    ! Ghost state of the two layers
    do l = 1, 2
      associate ( P => blk % P(:,Im-l*guide(Fm,1),Jm-l*guide(Fm,2),Km-l*guide(Fm,3)) )
      rhoG = pG(l) / ( f_Rtot ( YG(:,l) ) * TG(l) )
      P(1:nsc) = YG(:,l) * rhoG
      P(nu:nw) = uG(:,l)
      P(np) = pG(l)
      P(np+1:nprim) = phiG(:,l) * rhoG
      end associate
    enddo

  end subroutine Ghost_Wall


  !> Turbulence variables of the two layers of ghost cells of an SRM grain face (502) in a RANS run, after
  !> Ghost_Extrapolate: 2 phi_w - phi_l times the density of the ghost cell (phi = Z/rho; l = 1 the boundary cell, l = 2
  !> the next one along the direction of the face, the cell itself in a block of one cell), the rule of the walls
  !> (Ghost_Wall), with the wall values that the boundary flux of the grain takes (SRM_Wall_Turbulence): those of the
  !> blowing wall while the surface burns, of the wall at the flame temperature when it does not. The other variables keep
  !> the extrapolation. MOSE_Lib_Stencil_Metrics puts these values at 2 x_f - x_1.
  subroutine Ghost_SRM_Turbulence ( rec, blk )
    use MOSE_Lib_BC_Fluxes_GSI_SRM, only: SRM_Wall_Turbulence
    implicit none
    type(MOSE_bc_type), intent(in)       :: rec
    type(MOSE_block_type), intent(inout) :: blk
    ! Local
    integer  :: Im, Jm, Km, Fm, c(3), l
    real(R8) :: rho_w, rw(nrans), phiw(nrans)

    Im = rec % i ; Jm = rec % j ; Km = rec % k ; Fm = rec % f
    call SRM_Wall_Turbulence ( Im, Jm, Km, Fm, blk, rec % aCoeff, rec % n, rec % pRef, rec % Taf, rec % ci, rho_w, rw )
    phiw = rw / rho_w
    do l = 1, 2
      c = [Im, Jm, Km]
      if ( l == 2 .and. blk % dim((Fm+1)/2) > 1 ) c = c + guide(Fm,:)
      associate ( G => blk % P(:,Im-l*guide(Fm,1),Jm-l*guide(Fm,2),Km-l*guide(Fm,3)), Pc => blk % P(:,c(1),c(2),c(3)) )
      G(nt:nprim) = sum ( G(1:nsc) ) * ( 2d0 * phiw - Pc(nt:nprim) / sum ( Pc(1:nsc) ) )
      end associate
    enddo

  end subroutine Ghost_SRM_Turbulence


  !> q_1 + psi (q_1 - q_2) on the first ghost layer and q_1 + 2 psi_2 (q_1 - q_2) on the second, with psi <= sx and
  !> psi_2 <= psi the largest values that keep the ghost value within [q_1/2, 2 q_1] (with the sign of q_1; a zero q_1 is
  !> copied). cut: the bound acted on the first layer.
  pure subroutine Bounded_Extrapolation ( q1, q2, sx, g1, g2, cut )
    implicit none
    real(R8), intent(in)  :: q1, q2, sx
    real(R8), intent(out) :: g1, g2
    logical,  intent(out) :: cut
    ! Local
    real(R8) :: d, b

    d = q1 - q2
    b = huge(1d0)
    if ( q1 * d > 0d0 ) then
      b = q1 / d                    ! away from zero: |g| <= 2 |q1|
    elseif ( q1 * d < 0d0 ) then
      b = -0.5d0 * q1 / d           ! towards zero: |g| >= |q1|/2
    elseif ( abs(d) > 0d0 ) then
      b = 0d0                       ! q1 = 0
    end if
    g1 = q1 + min ( sx, b ) * d
    g2 = q1 + 2d0 * min ( sx, 0.5d0 * b ) * d
    cut = b < sx

  end subroutine Bounded_Extrapolation


  !> One line at the end of a run in which the cut of an extrapolated quantity acted in a ghost fill of a wall or 400 face, or
  !> the temperature of a 302 wall took the floor T_w/2: the fills of each kind, over the whole run, its grid levels and ranks.
  !> Either says that the mesh does not resolve the profile next to the face; a run in which neither acted prints no line.
  !> Collective: every rank calls it, and both reductions, in every run, so that no rank takes another MPI path; rank 0
  !> writes the line when one of the two totals is not zero.
  subroutine Report_Wall_Ghosts ()
    use MOSE_Mod_MPI, only: mpi_is_root, mpi_allreduce_sum_r8
    implicit none
    real(R8) :: ncut, nfloor

    call mpi_allreduce_sum_r8 ( real(n_ghost_cut, R8), ncut )
    call mpi_allreduce_sum_r8 ( real(n_ghost_floor, R8), nfloor )
    if ( mpi_is_root .and. ( ncut > 0d0 .or. nfloor > 0d0 ) ) &
      write(*,'(A,I0,A,I0)') ' Ghost cells of walls and 400 faces, fills with the extrapolation cut: ', &
                             nint(ncut, I8), ', with the 302 temperature at T_w/2: ', nint(nfloor, I8)

  end subroutine Report_Wall_Ghosts

end module MOSE_Lib_Ghost