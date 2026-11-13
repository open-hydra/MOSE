!> @brief Multirate (local time stepping) integration of the explicit, time-accurate step.
!>
!> Recursive Heun (SSP-RK2 in Shu-Osher form) with one level per cell: the "Prop. RK2" scheme of Jeanmasson, Mary
!> and Mieussens (2019), with two differences (the coarse corrector reads the final values of the finer cells, where
!> their Table 16 reads the fine predictor; the recursion over several levels, with bands of one cell, is ours), with
!> the coarse-data interpolant of Hoang, Ju, Leng and Wang (2020, eq. 3.2, with the full-stage predictor of their
!> App. B.1, eq. B.4) and the conservative correction of the faces between levels (refluxing).
!> See docs/theory/time-integration.md, section "Multirate time stepping". Notation:
!>
!>   H           macro step; level l advances with sub-steps h_l = H/2**l (l = 0 ... L, L the deepest level of the map);
!>   dtcap_c     admissible local step of cell c (Compute_dt: CFL and VNN limits);
!>   Lhat        depth for which the map of the adaptive step is built: Htry = 2**Lhat dt_min/sigma (at most H_target), Lhat
!>               the largest up to max-level at which the cell with the largest dtcap keeps level 0 (Depth_Cap), so
!>               H <= max_c dtcap_c/sigma; L = Lhat = min(max-level, floor(K)), K = log2(max_c dtcap_c/min_c dtcap_c) the
!>               octaves of the local steps, unless H_target sets Htry (then Lhat = max-level and L the depth H_target gives);
!>   lev_c       = min(max-level, max(0, ceiling(log2(x_c) - 1e-12))), x_c = sigma*Htry/dtcap_c and, with the interface
!>               margin (nonlocal-cfl = redo, the default, or interface-margin), x_c = max(x_c, kappa*Htry/min_v dtcap_v)
!>               over the face neighbours v in the block, kappa = max(cfl (1 + theta), 1.1); then nesting (|lev| jump <= 1 across
!>               internal faces) and equal levels on the two sides of every block connection, iterated;
!>   H           = min(H_target, min_c dtcap_c*2**lev_c/sigma, and with the margin min_c min_v dtcap_v*2**lev_c/kappa), so
!>               every cell advances with h_c <= dtcap_c/sigma.
!>
!>   ADVANCE(l, tau):  PO := P on the level-l cells; STAGE(l, 1, tau);
!>                     if l < L: stash the level-l cells read by finer stages; ADVANCE(l+1, tau);
!>                               ADVANCE(l+1, tau + h_l/2); restore their predictor;
!>                     STAGE(l, 2, tau + h_l).
!>   STAGE(l, k, t):   domain%time = t; the coarser cells read by the stage (halo: union of the stencils of the
!>                     faces of the level, including the sources of the ghost records read) take the linear
!>                     interpolant in conservative variables between the start state and the predictor of their
!>                     own step, theta = (t - tau_m)/h_m, with exact copies at theta = 0 and 1; ghost records read
!>                     by the stage are refreshed (viscous runs: then the tangential ghost values Pg(:,3:6) of the
!>                     connection records of the level, as Fill_Ghost_Cell, and the viscous variables blk%Visc of
!>                     every cell and ghost cell that the viscous fluxes of the stage read, halo and diagonals
!>                     included, from their primitive variables just set, as Viscous_Variables after
!>                     Fill_Ghost_Cell); fluxes on the faces with a level-l cell
!>                     with the routines of the single-rate step, direction by direction the convective then the
!>                     diffusive faces, as Fluxes_blk; conservative correction on the faces between levels, for the
!>                     convective and the diffusive flux each (coarse stage 1: own flux, saved; each of the four fine
!>                     stages: 1/4 of its flux accumulated; coarse stage 2: 2*accumulated - saved); boundary records of
!>                     the level, the connection faces with one flux per face computed by the record that owns it for
!>                     the two cells (Connection_Fluxes); RANS runs: the explicit turbulence sources of the level-l
!>                     cells, which read the ghost cells of the fluxes (Ghost_Wall at the walls), as the
!>                     single-rate stage; planar quasi-2D and 1-D cases with symmetry k faces: no z momentum
!>                     (Zero_W_Residual), 1-D cases with symmetry j faces: no y momentum (Zero_V_Residual);
!>                     Heun update of the level-l cells only, U := U0 + c_k (U - U0 + h_l L), c = (1, 1/2).
!>
!> RANS runs (Spalart-Allmaras, SST): the sources are explicit in time-accurate runs and the local step dtcap_c
!> includes a bound d_hat of the rate of their destruction D_q/q from the state of the cell alone (RANS_Stiffness_Bound,
!> as Compute_dt), so the levels and h_c <= dtcap_c/sigma inherit it: h_c D_q/q <= h_c d_hat <= 1/sigma. With one level
!> (L = 0) the stage calls the block routine of the single-rate stage (RANS_Source_Terms), so that the macro step is
!> the single-rate step to the bit; with more levels it calls its per-cell copies on the cells of the level (the copies
!> are compared with the originals statement by statement by the test mr-onelevel).
!>
!> Several MPI processes (every block on one process): the level map is the same as with one process (the
!> partner levels of the connection records between processes are exchanged at every sweep of Build_Map, and the
!> sweeps stop when no process changes its map); H and L are global. Every process runs the same recursion, so the
!> stages of a level are simultaneous on all processes. A stage of level l exchanges only what it reads from other
!> processes: the ghost cells of the connection records whose partner block is on another process and whose ghost
!> positions the stages of level l read here, and in viscous runs the tangential values Pg(:,3:6) of the level's
!> connection records. The process of the partner block sends them at the stage time: cells of level l or finer as they
!> are, coarser cells interpolated by the same interpolant (they join its halo of level l). The lists ("needs") are
!> agreed between the two processes at every build of the lists, record by record, as bit masks of the levels. Both
!> sides of a connection have the same level, so no face between levels crosses a connection and the conservative
!> correction stays inside the process. MPI calls are made by the master thread only (MPI_THREAD_FUNNELED).
!>
!> Stage times are integers ("ticks") in units of H/2**(L+1): theta is an exact dyadic fraction. The module writes no
!> file of its own: with trace = .true. the stage sequence, the level map and the halo lists of the first macro step go
!> to the standard output (lines MR-TRACE, MR-LEVEL, MR-HALO; with several processes those of the blocks of process 0),
!> and in viscous runs, stage by stage, the viscous variables read by the fluxes checked against their primitive
!> variables (lines MR-VISC).
!> With L = 0 and sigma = 1 one macro step is one single-rate Heun step (bit-identical, test mr-onelevel). Scope: Euler
!> and Navier-Stokes equations, laminar or with the Spalart-Allmaras (SA, SA-R) or SST model (k-energy-coupling and either
!> sst-production), and SA-RC and SST-RC with max-level 0; the SRM grain in its mode 2; one or several MPI processes, each
!> with OpenMP threads; set-up rejections: Multirate_Check_Config (keys) and Check_Domain (first macro step); at the start of
!> every macro step a local step that is not a positive finite number stops the run (Check_Local_Steps).
module MOSE_Mod_Multirate
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Config_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Mod_Timers, only: timer_wtime, timer_comm_begin, timer_comm_end, timer_sync_begin, timer_sync_end, &
                             timer_xstart_enable, timer_xstart_begin, timer_xstart_end

  implicit none
  private
  public :: Multirate_Check_Config, Multirate_Advance, Multirate_Report, Multirate_Reading_Error
  !> The last macro step ended exactly on an output instant t_k = k*sol-dtime (read by Explicit_Step and by the
  !> post-processing, which write the time-based solution at that step when the multirate is enabled)
  logical, public, save :: mr_output_landed = .false.

  integer, parameter :: LMAX_CODE = 50     !< stage times in 64-bit ticks of H/2**(L+1): exact while L+1 <= 52; at L = 50 two
                                           !< successive stage times differ by 2**-51 H, at least 2 units in the last place of H

  !> Cell lists
  type :: cell_list_t
    integer :: n = 0
    integer, allocatable :: b(:), i(:), j(:), k(:)
  end type cell_list_t

  !> Coarser cells (all of one level m) interpolated at the stages of a finer level
  type :: halo_list_t
    integer :: n = 0
    integer, allocatable :: b(:), i(:), j(:), k(:), u(:)   !< u: index in the stash of level m
  end type halo_list_t

  !> Segments of internal faces of one direction: faces a ... e of the row (b, p, q)
  type :: seg_list_t
    integer :: n = 0
    integer, allocatable :: b(:), p(:), q(:), a(:), e(:)
  end type seg_list_t

  type :: level_t
    type(cell_list_t) :: cells                        !< cells of the level
    type(seg_list_t)  :: seg(3)                       !< internal faces with at least one cell of the level
    integer :: nbc = 0, ngrp = 0
    integer, allocatable :: bcr(:)                    !< boundary records of the level's cells, in the order of BC_Fluxes
    integer, allocatable :: grp(:)                    !< start of each (block, face) group in bcr; grp(ngrp+1) = nbc+1
    integer :: ngh = 0
    integer, allocatable :: ghr(:)                    !< ghost records read by the stages of the level
    integer :: npg = 0
    integer, allocatable :: pgr(:)                    !< viscous runs: connection records of the level's cells, whose
                                                      !< tangential ghost values Pg(:,3:6) the stages refresh
    integer :: ncn = 0
    integer, allocatable :: cnr(:)                    !< connection records of the level's cells that own their face:
                                                      !< they compute its flux for both cells (Connection_Record)
    type(cell_list_t) :: visc                         !< viscous runs: cells and ghost cells whose viscous variables
                                                      !< blk%Visc the stages refresh (Visc_Of_Level)
    type(halo_list_t), allocatable :: halo(:)         !< halo(m), m = 0 ... l-1
    type(cell_list_t) :: used                         !< cells of the level read by finer stages
    real(R8), allocatable :: Ws(:,:)                  !< predictor (primitive) of the used cells
    real(R8), allocatable :: Ub(:,:), Us(:,:)         !< start state and predictor, conservative, of the used cells
  end type level_t

  type :: mr_block_t
    integer, allocatable :: lev(:,:,:)                !< level of every cell; -1 on ghost positions
    integer, allocatable :: ifc(:,:,:,:)              !< ifc(d,i,j,k): id of the interface face between (i,j,k) and its +d neighbour
    integer, allocatable :: grec(:,:,:)               !< ghost position -> boundary record that fills it (0: none)
    real(R8), allocatable :: P0(:,:,:,:)              !< state at the start of the macro step (residual iterations)
    logical, allocatable :: cflag(:,:,:)              !< nonlocal-cfl = count or redo: the cell violated the certificate in
                                                      !< this macro step
    real(R8), allocatable :: dmin(:,:,:)              !< nonlocal-cfl = redo: smallest dcap of the cell seen in the attempt
  end type mr_block_t

  !> MPI: the connection records exchanged with one other process. rrec: records whose block is here and whose partner
  !> block is there (we receive their ghost cells, and in viscous runs their Pg(:,3:6)); srec: records whose partner
  !> block is here and whose block is there (we send them). Both lists are in increasing record order, so the rrec list
  !> of one process and the srec list of the other match entry by entry.
  type :: peer_t
    integer :: rank = -1
    integer :: nr = 0, ns = 0
    integer, allocatable :: rrec(:), srec(:)
    integer(I8), allocatable :: rmP(:), rmG(:)        !< on rrec: bit l set when the stages of level l read them here
    integer(I8), allocatable :: smP(:), smG(:)        !< on srec: bit l set when the stages of level l read them there
  end type peer_t

  !> MPI: positions, in the srec or rrec list of one peer, of the records exchanged at the stages of one level
  type :: pos_list_t
    integer :: n = 0
    integer, allocatable :: q(:)
  end type pos_list_t

  !> MPI: what the stages of one level exchange, per peer: ghost cells (P), tangential values Pg(:,3:6) (G) and the
  !> share of the connection flux of the records whose owner record is on the other process (F)
  type :: lev_x_t
    type(pos_list_t), allocatable :: sP(:), rP(:), sG(:), rG(:), sF(:), rF(:)
    integer :: nsP = 0, nrP = 0, nsG = 0, nrG = 0     !< totals over the peers
    integer :: nsF = 0, nrF = 0                       !< shares of the connection fluxes (F): sent, received
  end type lev_x_t

  !> nonlocal-cfl = redo: the first state that is not physical of an attempt on this process, the one with the smallest
  !> key (number of the stage in the attempt, 0 for a face and 1 for a cell, block, direction of the face (0 for a cell),
  !> k, j, i), which does not depend on the number of threads; printed when the macro step falls below 2**-10 (Redo_Stop)
  type :: event_rec_t
    integer(I8) :: key(7) = huge(1_I8)
    integer :: kind = 0                               !< 1 cell, conservative variables; 2 cell, primitive; 3 face
    integer :: lev = 0, kst = 0                       !< level and stage of Heun of the stage
    integer :: flags(2) = 0                           !< face: error flags of check_gas_state, left and right
    real(R8) :: time = 0.0_R8                         !< time of the stage
    real(R8), allocatable :: v(:,:)                   !< cell: v(:,1); face: v(:,1:2) left and right states, v(:,3:6) the
                                                      !< stencil of cells -1 ... 2
  end type event_rec_t

  ! ---- state of the integrator (shared by the threads) ----------------------------------------------------------
  logical, save :: initialised = .false.
  type(multirate_t), save :: cfg0                       !< keys at start-up (fixed for the whole run)
  integer(I8), save :: n_macro = 0                      !< macro steps done in this run
  real(R8), save :: dtmin_first = 0.0_R8                !< dt_min of the first macro step (unit of fixed-macro-step)
  real(R8), save :: max_cap_ratio = 0.0_R8              !< max over the run of h_c/(dtcap_c/sigma)
  real(R8), save :: last_cap_ratio = 0.0_R8
  logical, save :: cap_warned = .false.
  integer(I8), save :: n_ratio_gt = 0                   !< macro steps with h_c/(dtcap_c/sigma) > 1 + 1e-12 at their start
  real(R8), save :: rpost_last = 0.0_R8                 !< largest r_post = h_c(n)/dtcap_c(t_n+1) of the last evaluation
  real(R8), save :: rpost_cnt(2) = 0.0_R8               !<   and its cells with r_post > 1 and > 1/cfl
  real(R8), save :: rpost_max_run = 0.0_R8              !<   largest r_post over the run
  integer(I8), save :: n_rpost_eval = 0, n_rpost_gt1 = 0, n_rpost_gtcfl = 0 !< macro steps evaluated, and with a cell
                                                        !<   above 1 and above 1/cfl
  real(R8), save :: rnb_last = 0.0_R8                   !< largest r_nb = h_c(n)/dtcap_v(t_n+1), v a face neighbour of c in
                                                        !<   its block (Gnedin et al. 2018, condition (4)), last evaluation
  real(R8), save :: rnb_cnt(3) = 0.0_R8                 !<   cells c with r_nb > 1, > 1/cfl, and with r_post or r_nb > 1
  real(R8), save :: rnb_max_run = 0.0_R8                !<   largest r_nb over the run
  integer(I8), save :: n_rnb_gt1 = 0, n_rnb_gtcfl = 0, n_any_gt1 = 0 !< macro steps with a cell: r_nb > 1, r_nb > 1/cfl,
                                                        !<   r_post or r_nb > 1
  ! nonlocal-cfl (count, redo) and interface-margin (docs/theory/time-integration.md, Non-local CFL condition)
  logical, save :: cert_on = .false.                    !< count the violations of the certificate of condition (4)
  logical, save :: redo_on = .false.                    !< nonlocal-cfl = redo with levels = recompute, no fixed-macro-step
  logical, save :: margin_on = .false.                  !< interface margin kappa (interface-margin or redo), max-level > 0
  logical, save :: cert_rejected = .false.              !< the last attempt violated the certificate (all processes)
  logical, save :: cert_bad = .false.                   !< the last attempt gave a dcap that is not a positive number
  integer(I8), save :: cert_nn = 0, cert_nn_step = 0    !< dcap not a positive number: one stage (reduction), the attempt
  integer(I8), save :: n_att_total = 0, n_redone = 0, n_halved = 0 !< attempts, macro steps redone, halvings of H
  integer(I8), save :: n_forced = 0                     !< force-redo (test device): attempts rejected by it
  integer(I8), save :: n_unlanded = 0                   !< macro steps cut to land on an output instant whose H a redo halved
  integer(I8), save :: att_hist(5) = 0                  !< macro steps with 1, 2, 3, 4, 5 or more attempts
  integer, save :: att_last = 1, att_max = 0            !< attempts of the last macro step, largest over the run
  logical, save :: first_attempt = .true.               !< the attempt runs on the dtlocal of Compute_dt (not on dmin)
  ! nonlocal-cfl = redo, a state that is not physical inside an attempt that the redo can redo (redo_on, L > 0): a cell
  ! (Compute_Cons_NoIRS) or a face (Face_Event) that takes one is counted and recorded, the attempt is rejected and redone
  ! with H halved on the same map (docs/theory/time-integration.md, Non-local CFL condition)
  logical, save :: ev_on = .false.                      !< the attempt can be redone
  logical, save :: ev_seen = .false.                    !< this process took a state that is not physical in the attempt
  logical, save :: frozen = .false.                     !< (stage) it took one in an earlier stage of the attempt: the stage
                                                        !<   computes no face of (4) and updates no cell in (6)
  logical, save :: defer_on = .false.                   !< (stage) a face that is not physical is counted (not at tick 0)
  logical, save :: cert_event = .false.                 !< the last attempt took a state that is not physical (all processes)
  logical, save :: cert_ev_face = .false.               !<   one of them at a face
  integer(I8), save :: ev_cells = 0, ev_faces = 0       !< cells and faces that took one, this process, the attempt
  integer(I8), save :: n_ev_rej = 0, n_ev_rej_cell = 0, n_ev_rej_face = 0 !< attempts rejected for it in the run: all, with
                                                        !<   a cell, with a face
  integer(I8), save :: cur_s = 0                       !< (stage) number of the stage in the attempt
  integer, save :: cur_l = 0, cur_k = 0                 !< (stage) its level, its stage of Heun
  real(R8), save :: cur_time = 0.0_R8                   !< (stage) its time
  type(event_rec_t), save :: ev_rec                     !< the first state that is not physical of the attempt, this process
  real(R8), save :: theta_lim = 1.0_R8, kappa_m = 1.0_R8 !< theta of the limiter, kappa = max(cfl (1 + theta), KAPPA_FLOOR)
  ! floor of the interface margin: with it a rejection for a violation of (4) needs a fall of the admissible step of a cell of
  ! more than 1 - 1/kappa within the macro step (docs, Appendix A, Proposition); the value is the fall that the check of
  ! the step of Castro tolerates (change_max), not a value that the scheme gives. Below the default safety 2: with
  ! safety >= kappa the margin adds nothing on the faces between equal levels
  real(R8), parameter :: KAPPA_FLOOR = 1.1_R8
  integer(I8), save :: cert_nv = 0, cert_nm = 0         !< one stage: violations, mismatches with dtlocal at t_n (reductions)
  real(R8), save :: cert_rq = 0.0_R8                    !< one stage: largest hhat/dcap (reduction)
  integer(I8), allocatable, save :: cert_viol(:)        !< violations (cell, sub-step start) of the macro step, per level
  integer(I8), allocatable, save :: cert_cells(:)       !< cells of each level with a violation in the macro step
  integer(I8), allocatable, save :: cert_viol_run(:), cert_cells_run(:) !< the same, summed over the run
  integer(I8), allocatable, save :: ring_cells(:), ring_cells_run(:) !< cells raised by kappa before nesting (last build; summed)
  real(R8), allocatable, save :: cert_wall(:)           !< wall time of the certificate per level (master thread)
  integer(I8), save :: n_cert_rej = 0, n_cert_eval = 0  !< macro steps with a violation (rejected at the first attempt)
  integer(I8), save :: cert_mism_run = 0                !< certificate at t_n against dtlocal: mismatches over 1e-12
  integer(I8), save :: cert_mism_step = 0               !<   in the macro step (this process)
  real(R8), save :: cert_rq_step = 0.0_R8, cert_rq_run = 0.0_R8
  integer(I8), save :: n_rings = 0                      !< builds of the map with kappa
  integer, save :: Lu = -1                             !< deepest level of the current map
  integer, save :: lu_run = 0                           !< deepest level of the maps of the run (lines at the end of the run)
  integer(I8), save :: n_deeper = 0                     !< macro steps whose accepted map a redo made deeper than the first
  logical, save :: have_map = .false.
  integer, save :: nif = 0                              !< interface faces of the current map
  integer, parameter :: FX_CONV = 1, FX_DIFF = 2        !< kind of face flux: convective, diffusive (viscous runs)
  integer, parameter :: RANS_SA = 1, RANS_SST = 2       !< turbulence model of a RANS run (rans_kind)
  integer, save :: rans_kind = 0                        !< 0 (no RANS), RANS_SA or RANS_SST: set by Check_Domain
  real(R8), allocatable, save :: acc(:,:,:), gsave(:,:,:) !< accumulated fine flux and saved coarse flux, per interface
                                                        !< face and kind of flux: acc(:, face, kind)
  type(mr_block_t), allocatable, save :: mb(:)
  type(mr_block_t), allocatable, save :: prev_lev(:)    !< map of the last build (interior cells)
  type(level_t), allocatable, save :: lv(:)
  integer(I8), allocatable, save :: cells_per_level(:)
  real(R8), save :: t_start, H_mac                      !< time at the start of the macro step, macro step
  integer(I8), save :: ticks_H                          !< 2**(Lu+1)
  logical, save :: trace_now = .false.
  real(R8), allocatable, save :: stage_wall(:)          !< wall time per level of the stages (master thread, after the
                                                        !< barriers; timer_wtime of MOSE_Mod_Timers, in every build)
  real(R8), allocatable, save :: xch_wall(:)            !< several processes: wall time per level spent waiting for the
                                                        !< messages of the exchanges of the stages and copying them from
                                                        !< the receive buffers (Exchange_P, Exchange_F, Exchange_G; master
                                                        !< thread), part of stage_wall and of the exchange wait of
                                                        !< MOSE_Mod_Timers
  ! ---- MPI (several processes) ---------------------------------------------------------------------------------------
  integer, save :: npeer = 0                            !< processes this one exchanges connection records with
  type(peer_t), allocatable, save :: peer(:)
  type(lev_x_t), allocatable, save :: xl(:)             !< exchanges of the stages of each level (0:Lu)
  integer, allocatable, save :: rpeer(:), rpos(:)       !< record -> peer and position in its rrec (0: not received)
  real(R8), allocatable, save :: xsbuf(:), xrbuf(:)     !< send and receive buffers of the stage exchanges
  integer, allocatable, save :: xreq(:)                 !< MPI requests of one exchange (2*npeer)
#ifdef USE_MPI
  integer, allocatable, save :: xstat(:,:)              !< their statuses (MPI_STATUS_SIZE, 2*npeer): no MPI_STATUSES_IGNORE,
                                                        !< which gfortran does not recognise with -fno-underscoring
#endif
  ! Tags of the multirate messages, one per exchange, as those of Mod_GhostExchange. Every
  ! receive names the rank of the sender, every exchange sends at most one message per pair of ranks and completes its
  ! messages before it returns, and no two exchanges of the multirate are in flight together (the master thread runs
  ! them one after the other), so a tag has only to tell an exchange from the others (MPI 4.1, Sec. 3.2.4: a receive
  ! takes a message whose source, tag and communicator match its own). They follow those of Mod_GhostExchange and of
  ! the output (0 to 6; the test mr-onelevel checks that no value is shared), and all are below 32767, the least
  ! MPI_TAG_UB the MPI standard allows (Sec. 3.2.3).
  integer, parameter :: TAG_MR_P = 7, TAG_MR_G = 8, TAG_MR_F = 9, TAG_MR_MASK = 10, TAG_MR_LEV = 11

contains

  ! =================================================================================================================
  !  SET-UP REJECTIONS
  ! =================================================================================================================

  !> Checks on the keys, right after input.ini is read (point A): items whose own set-up
  !> could stop MOSE first, or that come from an external program. Every rank runs it.
  subroutine Multirate_Check_Config ()
    use MOSE_Mod_MPI, only: mpi_is_root
    implicit none

    if (.not. obj_multirate%enabled) return

    ! keys of the multirate itself
    if (obj_multirate%max_level < 0 .or. obj_multirate%max_level > LMAX_CODE) &
      call Reject('max-level', 'max-level must be between 0 and 50: the integrator stores its per-level data and its '// &
                  'stage times, integers of H/2**(L+1), for levels 0 to 50')
    if (obj_multirate%safety < 1.0_R8) &
      call Reject('safety', 'safety must be >= 1')
    ! safety, at least 1 (the rule of the registry), divides the admissible step of every cell: no other bound. The default
    ! 2 gives condition (3) of Gnedin et al. on the state at the start of the macro step without the interface margin,
    ! every cell within the admissible step of each face neighbour; with nonlocal-cfl = redo or interface-margin the margin
    ! gives it for every safety >= 1 (docs, Non-local CFL condition)
    if (trim(obj_multirate%levels) /= 'fixed' .and. trim(obj_multirate%levels) /= 'recompute') &
      call Reject('levels', 'levels must be fixed or recompute')
    if (obj_multirate%recompute_every < 1) &
      call Reject('recompute-every', 'recompute-every must be >= 1')
    if (obj_multirate%macro_step_target < 0.0_R8 .or. obj_multirate%fixed_macro_step < 0.0_R8) &
      call Reject('macro step', 'macro-step-target and fixed-macro-step must be >= 0')
    if (trim(obj_multirate%nonlocal_cfl) /= 'off' .and. trim(obj_multirate%nonlocal_cfl) /= 'count' .and. &
        trim(obj_multirate%nonlocal_cfl) /= 'redo') &
      call Reject('nonlocal-cfl', 'nonlocal-cfl must be off, count or redo')
    if (obj_multirate%max_attempts < 1) &
      call Reject('max-attempts', 'max-attempts must be >= 1')
    ! cfl and vnn with levels and nonlocal-cfl = redo: within a macro step the certificate accepts a sub-step up to the whole
    ! step that cfl and vnn admit on its start state (to 1e-12), while the level map keeps every cell within 1/safety of that
    ! step only at the start of the macro step. On the checkerboard datum every limiter gives a zero slope, and Heun
    ! multiplies it by R(z) = 1 + z + z**2/2, z = -2 cfl a - 4 vnn b with a and b the convective and diffusive shares of the
    ! step, above 1 as soon as cfl a + 2 vnn b > 1: cfl <= 1 and vnn <= 1/2 are necessary for every accepted sub-step to be
    ! linearly stable, and sufficient with the first-order reconstruction and the limiters that have a linearisation (van
    ! Leer, MC, van Albada); for minmod and superbee no linear limit is proven (docs, The bounds of cfl and vnn). vnn enters
    ! the step only in Navier-Stokes runs (Compute_dt). With nonlocal-cfl = off or count nothing bounds the sub-steps after
    ! the start of the macro step, where every cell runs within cfl/safety and vnn/safety, and the set-up checks neither;
    ! so it is with redo and levels = fixed or a fixed-macro-step, where the certificate only counts (redo_on, Initialise)
    if (Redo_Keys(obj_multirate) .and. obj_time_scheme%cfl > 1.0_R8) &
      call Reject('cfl', 'cfl must be <= 1 when max-level > 0 and nonlocal-cfl = redo, with levels = recompute and no '// &
                  'fixed-macro-step: within a macro step the '// &
                  'certificate accepts a sub-step up to the whole step that cfl admits on its state (the level map keeps '// &
                  'every cell within 1/safety of that step only at the start of the macro step), and Heun is linearly '// &
                  'stable up to cfl = 1 with the first-order reconstruction and with the limiters that have a linearisation '// &
                  '(vanleer, mc, vanalbada), while above 1 the checkerboard mode grows with every limiter')
    if (Redo_Keys(obj_multirate) .and. Viscous_Keys() .and. obj_time_scheme%vnn > 0.5_R8) &
      call Reject('vnn', 'vnn must be <= 0.5 when max-level > 0 and nonlocal-cfl = redo, with levels = recompute and no '// &
                  'fixed-macro-step, in a Navier-Stokes run: within a '// &
                  'macro step the certificate accepts a sub-step up to the whole step that vnn admits on its state, and '// &
                  'Heun is linearly stable for the central diffusion up to vnn = 0.5, the end of its real stability '// &
                  'interval [-2, 0], while above it the checkerboard mode grows where the diffusion sets the step')
    ! restart (newrun = false): the level map is not written in the outputs, and fixed-macro-step is a multiple of
    ! dt_min of the first macro step of the run, which after a restart is another one
    if (.not. obj_sim_param%newrun .and. trim(obj_multirate%levels) == 'fixed') &
      call Reject('levels fixed', 'levels = fixed cannot restart (newrun = false): the level map is not written in the '// &
                  'outputs; restart with levels = recompute')
    if (.not. obj_sim_param%newrun .and. obj_multirate%fixed_macro_step > 0.0_R8) &
      call Reject('fixed macro step', 'fixed-macro-step cannot restart (newrun = false): H is a multiple of dt_min of '// &
                  'the first macro step of the run, which after a restart is another; set macro-step-target = H [s] instead')

    ! point A
    if (trim(obj_time_scheme%integration_variables) /= 'cons') &
      call Reject('integration-variables = '//trim(obj_time_scheme%integration_variables), &
                  'only integration-variables = cons is supported')
    ! IRS on the keys (irs, irs-beta): Assign_Setup refuses it in time-accurate runs and stops MOSE before point B
    if (obj_irs%enabled .or. obj_irs%beta > 0.0_R8) &
      call Reject('IRS', 'implicit residual smoothing (irs, irs-beta > 0) is not supported in time-accurate runs')
    if (trim(obj_chemistry%model) /= 'frozen' .and. trim(obj_chemistry%model) /= '') &
      call Reject('chemistry', 'chemistry = '//trim(obj_chemistry%model)// &
                  ' (Strang splitting) is not supported: set [MOSE-Physics] chemistry = frozen')
    if (trim(obj_soot%model) /= 'none' .and. trim(obj_soot%model) /= '') &
      call Reject('soot', 'soot-generation = '//trim(obj_soot%model)//' is not supported: set it to none')
    if (trim(obj_rot%model) == 'rigid-body') &
      call Reject('rotating frame', 'rotational-frame = rigid-body is not supported: set it to none')
    if (obj_multigrid%MGL > 1) &
      call Reject('multigrid', 'a grid sequence ([MOSE-Multigrid] levels > 1) is not supported')
    if (obj_sim_param%HYDRA_time_accurate .or. obj_sim_param%HYDRA_MG .or. obj_sim_param%HYDRA_postprocess) &
      call Reject('external time step', 'a time step or grid sequence imposed by an external program (HYDRA) is not supported')
    if (.not. obj_multirate%empty_external_function) &
      call Reject('external function', 'this program calls an External_Function at every stage; the multirate '// &
                  'integrator runs only from the standard MOSE program (src/app/main.f90)')

    if (mpi_is_root) write(*,'(A)') '   Multirate                       enabled (keys checked)'

  end subroutine Multirate_Check_Config


  !> Checks on the final set-up flags and on the mesh, at the first macro step (point B).
  subroutine Check_Domain ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: r, t
    character(len=16) :: tname

    if (obj_time_scheme%n_RK /= 2) then
      if (obj_time_scheme%n_RK == 1) then
        call Reject('time-scheme = euler', 'only time-scheme = RK2 (Heun) is supported')
      else
        call Reject('time-scheme = RK3', 'only time-scheme = RK2 (Heun) is supported; a third-order scheme '// &
                    'needs its own multirate coupling')
      endif
    endif
    ! the set-up keeps integration-variables as input.ini gives it and puts the name of the variables it integrates in
    ! integration_name (Assign_Integration_Variables)
    if (trim(obj_time_scheme%integration_name) /= 'Conservative') &
      call Reject('integration-variables = '//trim(obj_time_scheme%integration_variables), &
                  'only integration-variables = cons is supported')
    ! the key shock-detector turns the sensor on, and so does a Riemann solver whose name carries it (Assign_Setup)
    if (obj_shock_detector%id /= 0) &
      call Reject('shock detector', 'the shock detector ('//trim(obj_shock_detector%description)//') is not '// &
                  'supported: shock-detector turns it on, and so does a Riemann solver that carries it (riemann-solver = '// &
                  trim(obj_riemann%description)//')')
    if (.not. obj_time_scheme%time_accurate) &
      call Reject('time-accurate', 'the multirate integrator needs time-accurate = .true.')
    if (obj_time_scheme%rampa_cfl_iter > 0) &
      call Reject('cfl-rise-threshold', 'the CFL ramp (cfl-rise-threshold > 0) is not supported')
    ! assertions on the items checked at point A
    if (obj_irs%enabled) call Reject('IRS', 'implicit residual smoothing is not supported in time-accurate runs')
    if (obj_chemistry%imodel /= 0 .or. obj_chemistry%N_strang /= 1) call Reject('chemistry', 'chemistry is not supported')
    if (obj_soot%enabled) call Reject('soot', 'soot is not supported')
    if (obj_rot%enabled) call Reject('rotating frame', 'the rotating frame is not supported')
    if (obj_multigrid%MGL /= 1) call Reject('multigrid', 'a grid sequence is not supported')
    ! boundary conditions
    do r = 1, domain%nbound
      if (.not. is_local_block(domain%bc(r)%b)) cycle
      t = domain%bc(r)%type
      select case (t)
      case (102)
        call Reject('chimera', 'chimera boundaries (type 102) are not supported')
      case (103)
        call Reject('coupled wall', 'walls coupled to another solver (type 103) are not supported')
      case (502)
        ! a grain of mode 2 is a record of type 404 from the reading of bc.txt on (IO_BC.f90): only mode 1 comes here
        call Reject('SRM grain', 'an SRM grain in mode 1 (type 502, a wall that blows the products in) is not '// &
                    'supported: BC_Flux_Record has no branch for its flux (BC_SRM); a grain in mode 2, an inflow, runs '// &
                    'as the mass-flux inlet 404')
      case (503:506)
        call Reject('GSI wall', 'gas-surface interaction walls (types 503-506) are not supported')
      case (200)
        if (ndir == 3) call Reject('3D rotational periodicity', &
          'rotational periodic boundaries (type 200) in 3D are not supported (the 2D axisymmetric wedge is)')
      case (0, 101, 201, 300, 301, 302, 400:408, 410, 420)
        ! types with a boundary flux in BC_Flux_Record, or without one as in BC_Fluxes (0 null); 409 and 501 stop
        ! every run at the reading of bc.txt (Unapplied_Type, IO_BC.f90), before this check
      case default
        write(tname,'(I0)') t
        call Reject('boundary type', 'boundary type '//trim(tname)//' has no multirate boundary flux (the multirate '// &
                    'integrator knows the types of BC_Fluxes; a new type needs its branch in BC_Flux_Record)')
      end select
    enddo
    ! turbulence: the stages evaluate the sources of the Spalart-Allmaras (SA, and SA-R with the rotation correction of its
    ! production) and SST models: with one level the block routine of the single-rate stage (RANS_Source_Terms), with
    ! more levels the cell-wise copies of SA_Source_Blk and SST_Blk, k-energy-coupling, both forms of sst-production and
    ! sst-dilatation-dissipation included. Wilcox2006 has no evaluation of its sources in the stages. The rotation and curvature correction (SA-RC,
    ! SST-RC) reads terms that Compute_RC_Terms computes for the whole domain inside the block routine, which the cell-wise
    ! copies do not compute: with max-level 0 the stages call the block routine, with max-level > 0 the run is rejected.
    ! SAcomp (the compressibility correction of SA) is not among them: the cell-wise copy calls the correction, but the
    ! stages are not tested with it, and the run is rejected. The production of SA-R and the dilatation part of the
    ! compressible production of SST can be a sink outside the bound of the stiffness that the step takes (it bounds the
    ! destruction only), as in the single-rate step (docs, Turbulence models)
    rans_kind = 0
    if (model == 2) then
      select case (trim(obj_rans%model))
      case ('SA', 'SA-R', 'SA-RC')
        rans_kind = RANS_SA
      case ('SST', 'SST-RC')
        rans_kind = RANS_SST
      case default
        call Reject('turbulence model', 'turbulence = '//trim(obj_rans%model)//' is not supported: the stages of the '// &
                    'multirate integrator evaluate the sources of SA, SA-R and SST, and of SA-RC and SST-RC with '// &
                    'max-level = 0')
      end select
      if (obj_rans%SpalartShur .and. obj_multirate%max_level > 0) &
        call Reject('turbulence model', 'turbulence = '//trim(obj_rans%model)//' is not supported with max-level > 0: '// &
                    'its rotation and curvature terms are a field of the whole domain (Compute_RC_Terms), which the '// &
                    'cell-wise sources of the stages of a level do not compute')
      ! assertions: options that only model names outside the registry switch on (QCR, blowing correction, Reynolds
      ! stresses); SAcomp is rejected by the select case above
      if (obj_rans%SAcomp .or. obj_rans%QCR2000 .or. obj_rans%RSM .or. obj_rans%blowing_corr) &
        call Reject('turbulence model', 'the options of the turbulence model (compressibility correction, QCR, blowing '// &
                    'correction, Reynolds stresses) are not supported')
    endif

  end subroutine Check_Domain


  !> Every rank that finds the configuration outside the scope prints the message and flushes it before stopping:
  !> with several MPI processes the first error stop makes mpirun kill the others, so a message printed only by the
  !> root process could be lost (seen with 2 processes).
  subroutine Reject ( item, msg )
    use iso_fortran_env, only: output_unit
    use MOSE_Mod_MPI, only: mpi_rank_, mpi_size_
    implicit none
    character(len=*), intent(in) :: item, msg
    character(len=32) :: who
    who = ''
    if (mpi_size_ > 1) write(who,'(A,I0,A)') ' (process ', mpi_rank_, ')'
    write(*,'(A)') ' [ERROR] [MOSE-Multirate] '//trim(item)//': '//trim(msg)//trim(who)
    write(*,'(A)') '         (set [MOSE-Multirate] enabled = .false. to run the single-rate integrator)'
    flush(output_unit)
    error stop 'MOSE-Multirate: configuration not supported'
  end subroutine Reject


  !> input.ini read again during the run (Read_Inifile_Runtime, after the checks of the registry, of the keys fixed for
  !> the run, those of [MOSE-Multirate] among them, and of the set-up between keys): the checks that the multirate
  !> integrator makes at set-up or at its first macro step on the keys of other sections that a reading takes, made
  !> again on the values read: with max-level > 0 and nonlocal-cfl = redo, levels = recompute and no fixed-macro-step,
  !> cfl <= 1 and, in Navier-Stokes runs,
  !> vnn <= 0.5 (Multirate_Check_Config), and cfl-rise-threshold = 0 (Check_Domain). time-accurate, which Check_Domain
  !> checks too, is fixed for the run (KEY_FIXED): a reading cannot change it. The keys of [MOSE-Multirate]
  !> are those of the set-up (cfg0, fixed by Initialise at the first macro step; a reading comes after it,
  !> in the post-processing of an iteration). The error of the first check that fails, with the file, the section, the
  !> key and its value as written (Key_Text), as the other errors of a reading give a key; an empty string when
  !> every check passes. What the error does, ini-error decides (Read_Inifile_Runtime).
  function Multirate_Reading_Error () result ( out )
    use MOSE_Input_Registry, only: Key_Text
    implicit none
    character(len=:), allocatable :: out
    character(len=:), allocatable :: num

    num = trim(codename)//'-Numerics'
    out = ''
    if (Redo_Keys(cfg0) .and. obj_time_scheme%cfl > 1.0_R8) then
      out = '[ERROR] '//Key_Text(num, 'cfl')//': the multirate integrator needs cfl <= 1 (max-level > 0 and '// &
            'nonlocal-cfl = redo, with levels = recompute and no fixed-macro-step, of the set-up): within a macro step '// &
            'the certificate accepts a sub-step up to the whole '// &
            'step that cfl admits on its state, and Heun is linearly stable up to cfl = 1 with the first-order '// &
            'reconstruction and with the limiters that have a linearisation (vanleer, mc, vanalbada)'
    else if (Redo_Keys(cfg0) .and. Viscous_Keys() .and. obj_time_scheme%vnn > 0.5_R8) then
      out = '[ERROR] '//Key_Text(num, 'vnn')//': the multirate integrator needs vnn <= 0.5 in a Navier-Stokes run '// &
            '(max-level > 0 and nonlocal-cfl = redo, with levels = recompute and no fixed-macro-step, of the set-up): '// &
            'within a macro step the certificate accepts a '// &
            'sub-step up to the whole step that vnn admits on its state, and Heun is linearly stable for the central '// &
            'diffusion up to vnn = 0.5'
    else if (obj_time_scheme%rampa_cfl_iter > 0) then
      out = '[ERROR] '//Key_Text(num, 'cfl-rise-threshold')//': the multirate integrator runs without the CFL '// &
            'ramp (cfl-rise-threshold = 0)'
    endif
  end function Multirate_Reading_Error


  ! =================================================================================================================
  !  MACRO STEP
  ! =================================================================================================================

  !> One macro step: local steps, level map, macro step H, recursive sub-steps of all levels, time t_n + H.
  !> Called by Explicit_Step in place of the single-rate step; the residual and the counters stay there.
  subroutine Multirate_Advance ( domain )
    use MOSE_Mod_dt,  only: Compute_dt
    use MOSE_Mod_MPI, only: mpi_is_root, is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    real(R8) :: dtmin, Htry, Hcap, H, ratio, r_loc, hc, H_first, kout, tk, t_land
    logical  :: rebuild, resid_iter, forced, landed
    integer  :: b, i, j, k, l, natt, lu_first
    integer  :: rej_viol, rej_nn, rej_ev, rej_ev_face   ! rejected attempts of the macro step, by cause
    integer(I8) :: taus(0:LMAX_CODE), nst, nh0

    mr_output_landed = .false.
    if (.not. initialised) call Initialise(domain)
    ! interface margin kappa = max(cfl (1 + theta), KAPPA_FLOOR) with the cfl in force: input.ini read again during the run
    ! can change cfl (within 1, Multirate_Reading_Error), and the margin follows it from this macro step on;
    ! theta stays that of the set-up
    kappa_m = max(obj_time_scheme%cfl * (1.0_R8 + theta_lim), KAPPA_FLOOR)

    ! local admissible steps (dtlocal) and their minimum, as in the single-rate step (time-accurate RANS runs: with the
    ! bound of the stiffness of the explicit turbulence sources from the state of each cell, added by Compute_dt)
    domain % dtglobal = 1d5
    call Compute_dt ( domain, obj_time_scheme%cfl, obj_time_scheme%vnn, obj_time_scheme%rampa_cfl_iter )
    ! every local step a positive finite number, or the run stops here, before the map and H are built from them
    call Check_Local_Steps ( domain )
    dtmin = domain % dtglobal
    if (n_macro == 0) dtmin_first = dtmin

    ! The step each cell took in the last macro step against the step admitted on the state it reached (the
    ! dtlocal just computed), with the map of the last macro step, before Build_Map: no extra computation, no effect on H,
    ! on the map or on the results. Not at the first macro step of a run or of a restart: there is no map before it
    if (have_map .and. n_macro > 0) call R_Post ( domain )

    ! level map. levels = recompute: rebuilt at the iterations 1, 1 + k, 1 + 2k, ... of the run counter domain%iter
    ! (k = recompute-every), which a restart from a vtk solution carries on, so a restart on a multiple of k rebuilds
    ! where the uninterrupted run does; a restart rebuilds (the solution holds no map), and the map depends only on the state
    rebuild = .not. have_map
    if (have_map .and. trim(cfg0%levels) == 'recompute') &
      rebuild = mod(int(domain%iter, I8) - 1_I8, int(cfg0%recompute_every, I8)) == 0
    if (rebuild) then
      if (cfg0%fixed_macro_step > 0.0_R8 .and. trim(cfg0%levels) == 'recompute') then
        ! the map of the macro step that is taken: every cell at the lowest level within its limit for H (cells that would
        ! need more than L levels go beyond it, and the run reports the ratio)
        Htry = cfg0%fixed_macro_step * dtmin_first
      else
        ! adaptive step; with levels = fixed also with fixed-macro-step: the map of the adaptive step of the first macro
        ! step, kept, so that every h_c = H/2**lev_c scales with H (studies of the temporal order). Its depth is the
        ! largest up to max-level at which the cell with the largest local step keeps level 0 (Depth_Cap): the macro step
        ! never exceeds the step of the slowest cells divided by sigma, and no coarse level is left empty by the formula
        Htry = Htry_Of ( Depth_Cap ( domain, dtmin ), dtmin )
      endif
      call Build_Map ( domain, Htry )
      ! levels = fixed with a fixed macro step above the step of the map (2**L dt_min/safety, L the depth of the map, or
      ! macro-step-target): the finest cells start the first macro step beyond their limit h <= dtcap/safety, and every later one
      ! until their admissible step grows by the ratio of H to the step of the map, which never rises; the run goes on with the
      ! warning and the run line of the largest h/(dtcap/safety) below, as with levels = recompute where max-level is too low
    endif

    ! macro step: H_cap = min_c dtcap_c*2**lev_c/sigma (every cell within its limit), or the fixed step; with the
    ! interface margin also min_c min_v dtcap_v*2**lev_c/kappa over the face neighbours v in the block
    Hcap = Hcap_Of_Map ( domain )
    if (cfg0%fixed_macro_step > 0.0_R8) then
      H = cfg0%fixed_macro_step * dtmin_first
    else
      H = Hcap
      if (cfg0%macro_step_target > 0.0_R8) H = min(H, cfg0%macro_step_target)
    endif

    ! The macro step lands exactly on the next output instant t_k = k*sol-dtime (k integer, the product in double
    ! precision: the smallest t_k above the time) and on time-threshold. H is cut to reach the nearer of the two; a step
    ! that would end within 1e-12 H of it ends there, so the rounding of H leaves no step of a few units in the last place
    ! of H (H grows by at most 1e-12 H or half a unit in the last place of the time, whichever is larger; where that unit
    ! is above 1e-12 H a step of a few of those units can also remain, one more macro step). The map of the step is kept:
    ! every h_c = H/2**lev_c only shrinks, so every cell stays within
    ! its limit. The time is then assigned, not summed. The instants depend only on the time, so a restart lands where
    ! the uninterrupted run does.
    kout = aint(domain%time / obj_io%sol_dtime) + 1.0_R8
    if ((kout - 1.0_R8) * obj_io%sol_dtime > domain%time) kout = kout - 1.0_R8
    if (kout * obj_io%sol_dtime <= domain%time) kout = kout + 1.0_R8
    tk = kout * obj_io%sol_dtime
    t_land = tk
    if (obj_sim_param%time_threshold > domain%time) t_land = min(t_land, obj_sim_param%time_threshold)
    landed = domain%time + H >= t_land - 1.0e-12_R8 * H
    if (landed) then
      H = t_land - domain%time
      mr_output_landed = (t_land == tk)
    endif

    ! largest h_c/(dtcap_c/sigma) of this macro step
    r_loc = 0.0_R8
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        hc = H / 2.0_R8**mb(b)%lev(i,j,k)
        r_loc = max(r_loc, hc / (domain%blk(b)%dtlocal(i,j,k) / cfg0%safety))
      enddo; enddo; enddo
    enddo
    ratio = Max_R8_Process ( r_loc )
    last_cap_ratio = ratio
    max_cap_ratio = max(max_cap_ratio, ratio)
    if (ratio > 1.0_R8 + 1.0e-12_R8) n_ratio_gt = n_ratio_gt + 1
    if (ratio > 1.0_R8 + 1.0e-12_R8 .and. .not. cap_warned .and. mpi_is_root) then
      write(*,'(A,ES22.15,A)') ' [WARNING] [MOSE-Multirate] the macro step exceeds the limit of some cell: '// &
        'max h/(dtcap/sigma) = ', ratio, ' (fixed-macro-step); reported once'
      cap_warned = .true.
    endif

    ! residual iterations: keep the state at the start of the macro step (Explicit_Step compares P with PO)
    resid_iter = (obj_sim_param%iter_from_call == 1) .or. (mod(domain%iter, obj_io%res_diter) == 0)

    t_start  = domain % time
    trace_now = cfg0%trace .and. n_macro == 0 .and. mpi_is_root
    H_first  = H
    natt     = 0
    ! nonlocal-cfl = redo (condition (4) of Gnedin et al.; docs/theory/time-integration.md): an attempt whose
    ! certificate fails is rejected, the state goes back to t_n, the map is rebuilt for the H of the attempt from the
    ! smallest admissible step seen by every cell (dmin), so that from one redo to the next with the same H the levels
    ! only rise, up to max-level (the map can become deeper than that of the first attempt), and H is halved while a cell
    ! would need a level above max-level. An attempt in which a cell or a face takes a state that is not physical is
    ! rejected too, and redone with H halved once on the same map and dtlocal (the admissible steps seen after that state
    ! come from states that are not the solution). Below 2**-10 of the first attempt the run stops, and so it does when the
    ! attempts of the macro step rejected for one of these causes reach max-attempts (Redo_Stop). Without the redo a macro
    ! step has one attempt
    lu_first = Lu
    rej_viol = 0; rej_nn = 0; rej_ev = 0; rej_ev_face = 0
    nh0 = n_halved
    attempts: do
      natt = natt + 1
      first_attempt = natt == 1
      H_mac   = H
      ticks_H = ishft(1_I8, Lu + 1)
      taus    = 0_I8
      if (cert_on) then
        cert_viol = 0
        cert_rq_step = 0.0_R8
        cert_mism_step = 0
        cert_nn_step = 0
      endif
      if (redo_on) then
        do b = 1, domain % nb
          if (is_local_block(b)) mb(b)%dmin = domain%blk(b)%dtlocal
        enddo
      endif
      ! a state that is not physical: counted and recorded only in an attempt that the redo can redo
      ev_on   = redo_on .and. Lu > 0
      ev_seen = .false.
      frozen  = .false.
      ev_cells = 0
      ev_faces = 0
      ev_rec%key = huge(1_I8)
      cert_event = .false.
      cert_ev_face = .false.

      ! nst: the stages of the attempt, counted by every thread in Advance (the same count on every thread and process)
      !$omp parallel private(nst)
      nst = 0_I8
      if (resid_iter .or. redo_on) call Copy_P0 ( domain, .true. )
      call Advance ( domain, 0, taus, nst )
      !$omp end parallel

      if (cert_on .and. Lu > 0) call Cert_Close ( domain )
      ! force-redo (test device): the first attempt that the certificate accepts is rejected all the same and redone with
      ! the same H and map, so the state must come back to t_n to the bit and the results must not change
      forced = cfg0%force_redo .and. redo_on .and. natt == 1 .and. .not. (Lu > 0 .and. cert_rejected)
      if (.not. (redo_on .and. Lu > 0 .and. cert_rejected) .and. .not. forced) exit attempts

      ! rejected attempt
      !$omp parallel
      call Restore_State ( domain )
      !$omp end parallel
      domain % time = t_start
      trace_now = .false.
      if (forced) then
        n_forced = n_forced + 1
        cycle attempts
      endif
      if (cert_event) then
        ! a state that is not physical (on every process together: cert_event comes from the reduction of Cert_Close): the
        ! map and dtlocal of the attempt kept, H halved once. H_cap does not change (the same map, dtlocal and kappa) and
        ! the H of an attempt is at most H_cap (1 + 1e-12), plus half a unit in the last place of the time after a landing
        ! (it stretches H by no more than the larger of the two), so H/2 is below H_cap: no rebuild, no H_cap to compute
        rej_ev = rej_ev + 1
        if (cert_ev_face) rej_ev_face = rej_ev_face + 1
        H = 0.5_R8 * H
        n_halved = n_halved + 1
      else
        if (cert_bad) then
          rej_nn = rej_nn + 1
        else
          rej_viol = rej_viol + 1
        endif
        do b = 1, domain % nb
          if (is_local_block(b)) domain%blk(b)%dtlocal = min(domain%blk(b)%dtlocal, mb(b)%dmin)
        enddo
        call Build_Map ( domain, H )
        if (cert_bad) then
          H = 0.5_R8 * H
          n_halved = n_halved + 1
        endif
        Hcap = Hcap_Of_Map ( domain )
        ! within the tolerance of the formula of the levels, as the other comparisons of steps: a level kept by that
        ! tolerance after a rebuild asks for no halving
        do while (H > Hcap * (1.0_R8 + 1.0e-12_R8))
          H = 0.5_R8 * H
          n_halved = n_halved + 1
        enddo
      endif
      if (H < scale(H_first, -10)) call Redo_Stop ( H, H_first, natt, rej_viol, rej_nn, rej_ev, .false. )
      ! max-attempts: the attempts of the macro step rejected for a cause (the attempt rejected by force-redo left the loop
      ! iteration above and is not counted), checked after the floor, so that with max-attempts >= 11 a macro step whose
      ! rejections all halve H ends as without the limit; an accepted attempt leaves the loop before both checks
      if (rej_viol + rej_nn + rej_ev >= cfg0%max_attempts) &
        call Redo_Stop ( H, H_first, natt, rej_viol, rej_nn, rej_ev, .true. )
      ! A halved H no longer reaches the output instant or time-threshold that the first attempt landed on; the
      ! instant is reached by a later macro step
      if (landed .and. t_start + H < t_land - 1.0e-12_R8 * H) then
        landed = .false.
        mr_output_landed = .false.
        n_unlanded = n_unlanded + 1
      endif
    enddo attempts

    att_last = natt
    att_max  = max(att_max, natt)
    n_att_total = n_att_total + natt
    if (natt > 1) n_redone = n_redone + 1
    if (Lu > lu_first) n_deeper = n_deeper + 1
    att_hist(min(natt, 5)) = att_hist(min(natt, 5)) + 1
    if (resid_iter) then
      !$omp parallel
      call Copy_P0 ( domain, .false. )
      !$omp end parallel
    endif

    if (landed) then
      domain % time   = t_land               ! the instant itself
    else
      domain % time   = t_start + H
    endif
    domain % dtglobal = H
    n_macro = n_macro + 1

    ! Lines of the macro step, every shell-diter iterations: its measures (H, L, the largest h/(dtcap/sigma), cells and wall
    ! time per level, the largest r_post and r_nb) always; the lines of a situation only when it occurred: the certificate
    ! when the last attempt violated condition (4) or the macro step took more than one attempt, with the attempts rejected
    ! for a state that is not physical when there were some, the margin when it raised cells, the counts of r_post and r_nb
    ! when a cell went above 1 or above 1/cfl
    if (mpi_is_root .and. ((obj_sim_param%iter_from_call == 1) .or. (mod(domain%iter, obj_io%shell_diter) == 0))) then
      if (cert_on .and. Lu > 0 .and. (sum(cert_viol(0:Lu)) > 0 .or. att_last > 1)) then
        write(*,'(A,*(I10))') '   Multirate: certificate, violations per level (last attempt)', (cert_viol(l), l = 0, Lu)
        write(*,'(A,*(I10))') '   Multirate: certificate, cells with a violation per level (last attempt)', &
          (cert_cells(l), l = 0, Lu)
        write(*,'(A,ES22.15)') '   Multirate: certificate, largest hhat/dcap of the last attempt ', cert_rq_step
        write(*,'(A,*(ES10.3))') '   Multirate: certificate wall time per level [s]', (cert_wall(l), l = 0, Lu)
        if (redo_on) write(*,'(A,I0,A,I0,A,I0)') '   Multirate: attempts of the macro step ', att_last, &
          ', macro steps redone so far ', n_redone, ', halvings of H so far ', n_halved
        if (rej_ev > 0) write(*,'(A,I0,A,I0,A,I0,A,I0)') '   Multirate: attempts of the macro step rejected for a state '// &
          'that is not physical ', rej_ev, ' (with a face ', rej_ev_face, '), halvings of H in the macro step ', n_halved - nh0, &
          ', attempts rejected for it in the run ', n_ev_rej
      endif
      if (margin_on .and. sum(ring_cells(0:Lu)) > 0) &
        write(*,'(A,*(I10))') '   Multirate: cells raised by kappa per level (last map)', (ring_cells(l), l = 0, Lu)
      write(*,'(A,ES14.7,A,ES10.3,A,I3,A,ES10.3)') '   Multirate: H = ', H, ' s = ', H/dtmin, ' dt_min, L = ', Lu, &
        ', max h/(dtcap/sigma) = ', ratio
      write(*,'(A,*(I10))') '   Multirate: cells per level', (cells_per_level(l), l = 0, Lu)
      write(*,'(A,*(ES10.3))') '   Multirate: wall time per level [s]', (stage_wall(l), l = 0, Lu)
      if (npeer > 0) write(*,'(A,*(ES10.3))') '   Multirate: exchange time per level [s]', (xch_wall(l), l = 0, Lu)
      if (n_rpost_eval > 0) then
        write(*,'(A,ES22.15,A,ES22.15)') '   Multirate: r_post = h/dtcap(t_n+1) of the last macro step, max ', rpost_last, &
          '; r_nb = h/dtcap(t_n+1) of a face neighbour, max ', rnb_last
        if (any(rpost_cnt > 0.0_R8) .or. any(rnb_cnt > 0.0_R8)) write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') &
          '   Multirate: cells with r_post > 1: ', nint(rpost_cnt(1), I8), ', > 1/cfl: ', nint(rpost_cnt(2), I8), &
          '; with r_nb > 1: ', nint(rnb_cnt(1), I8), ', > 1/cfl: ', nint(rnb_cnt(2), I8), &
          '; with r_post or r_nb > 1: ', nint(rnb_cnt(3), I8)
      endif
    endif

  end subroutine Multirate_Advance


  !> nonlocal-cfl = redo: the macro step fell below 2**-10 of its first attempt (H after natt attempts), or, with capped,
  !> its attempts rejected for a cause reached max-attempts (H_first and H_mac, the H of its first and of its last attempt,
  !> printed). The run stops with the causes of the rejections of the macro step (nviol for a violation of condition (4),
  !> nnn for a local step that is not a positive number, nev for a state that is not physical) and, when the last attempt
  !> took a state that is not physical, the first one of that attempt in the order of the stages: the smallest number of
  !> the stage over the processes, and every process whose record has it prints it (a cell: its values; a face: the lines
  !> of the face of Convective_Flux).
  !> At max-attempts a last attempt without such a state gives its other cause: a local step that is not a positive
  !> number, or its violations of condition (4) and the cells that had one (sums over the processes made by Cert_Close).
  !> Collective: every process gets here together (H, H_first, the counts of the rejections and the decisions to reject
  !> come from the reductions of Cert_Close and Hcap_Of_Map), prints what it holds and waits for the others before the run
  !> ends (mpi_stop_all, exit code 1).
  subroutine Redo_Stop ( H, H_first, natt, nviol, nnn, nev, capped )
    use MOSE_Mod_MPI, only: mpi_is_root, mpi_rank_, mpi_size_, mpi_allreduce_min_r8, mpi_stop_all
    use MOSE_Lib_Convective, only: Interface_State_Lines
    implicit none
    real(R8), intent(in) :: H, H_first
    integer, intent(in) :: natt, nviol, nnn, nev
    logical, intent(in) :: capped
    real(R8) :: smin
    character(len=32) :: who
    character(len=128) :: why
    integer :: loc(5)

    if (mpi_is_root) then
      if (capped) then
        write(*,'(A,I0,A,ES14.7,A,ES14.7,A)') ' [ERROR] [MOSE-Multirate] nonlocal-cfl: the macro step was rejected in all '// &
          'its attempts (max-attempts = ', cfg0%max_attempts, '; H of the first attempt ', H_first, ' s, of the last ', &
          H_mac, ' s)'
      else
        write(*,'(A,ES14.7,A,I0,A)') ' [ERROR] [MOSE-Multirate] nonlocal-cfl: the macro step fell below 2**-10 of its '// &
          'first attempt (H = ', H, ' s after ', natt, ' attempts)'
      endif
      write(*,'(A,I0,A,I0,A,I0,A)') ' [ERROR]   attempts of the macro step rejected: ', nviol, ' for a violation of '// &
        'condition (4), ', nnn, ' for a local step that is not a positive number, ', nev, ' for a state that is not physical'
      if (cert_event) then
        write(*,'(A)') ' [ERROR]   the last attempt took a state that is not physical; the first one, in the order of '// &
          'the stages:'
      else if (.not. capped) then
        write(*,'(A)') ' [ERROR]   the admissible steps collapse within the macro step'
      else if (cert_bad) then
        write(*,'(A)') ' [ERROR]   the last attempt gave a local step that is not a positive number'
      else
        write(*,'(A,I0,A,I0)') ' [ERROR]   the last attempt violated condition (4): violations ', sum(cert_viol), &
          ', cells with a violation ', sum(cert_cells)
      endif
    endif
    if (cert_event) then
      call mpi_allreduce_min_r8 ( real(ev_rec%key(1), R8), smin )
      if (ev_rec%key(1) < huge(1_I8) .and. ev_rec%key(1) == nint(smin, I8)) then
        who = ''
        if (mpi_size_ > 1) write(who,'(A,I0,A)') ' (process ', mpi_rank_, ')'
        write(*,'(A,I0,A,I0,A,I0,A,ES14.7,A,A)') ' [ERROR]   stage ', ev_rec%key(1), ' of the attempt: level ', ev_rec%lev, &
          ', stage ', ev_rec%kst, ' of its sub-step, time ', ev_rec%time, ' s', trim(who)
        select case (ev_rec%kind)
        case (1)
          write(*,'(A,4(1X,I0),A)') ' [ERROR]   cell b, i, j, k:', ev_rec%key(3), ev_rec%key(7), ev_rec%key(6), &
            ev_rec%key(5), ': conservative variables not finite, or total density not positive'
          write(*,'(A,*(1X,ES12.4))') ' [ERROR]   conservative variables:', ev_rec%v(:,1)
        case (2)
          write(*,'(A,4(1X,I0),A)') ' [ERROR]   cell b, i, j, k:', ev_rec%key(3), ev_rec%key(7), ev_rec%key(6), &
            ev_rec%key(5), ': primitive variables not finite, negative pressure or temperature outside the tables'
          write(*,'(A,*(1X,ES12.4))') ' [ERROR]   primitive variables:', ev_rec%v(:,1)
        case (3)
          loc = int([ev_rec%key(3), ev_rec%key(4), ev_rec%key(7), ev_rec%key(6), ev_rec%key(5)])
          call Interface_State_Lines ( ev_rec%flags(1), ev_rec%flags(2), ev_rec%v(:,1), ev_rec%v(:,2), ev_rec%v(:,3:6), &
                                       loc )
        end select
      endif
    endif
    if (capped) then
      write(why,'(A,I0,A)') 'MOSE-Multirate: macro step rejected in all its attempts, max-attempts = ', cfg0%max_attempts, &
        ' (nonlocal-cfl = redo)'
    else
      why = 'MOSE-Multirate: macro step below 2**-10 of its first attempt (nonlocal-cfl = redo)'
    endif
    call mpi_stop_all ( trim(why) )

  end subroutine Redo_Stop


  !> The local steps of the macro step (blk%dtlocal, from Compute_dt) are positive finite numbers in every cell of every
  !> process: the depth of the map, the level map, H_cap and H come from them. A local step that is zero, negative,
  !> infinite or not a number would give a macro step that is not a positive finite number: the attempts would run with it,
  !> the halving of H while H > H_cap would not end with H_cap < 0, and the stop below 2**-10 of the first attempt would
  !> never come with H not a number. Such steps come from states that are not physical: a pressure and a density of opposite signs
  !> give a step that is not a number; a density and a pressure both negative in a viscous run, or a negative eddy viscosity that
  !> makes every diffusivity negative, a negative step where the diffusion rate outweighs the convective one; a zero pressure at
  !> rest in an inviscid run an infinite step; a density so small that the diffusion rate overflows a zero step.
  !> Such a local step stops the run at the start of the macro step: process 0 prints the number of such cells, the process
  !> that holds the first one in the order of the blocks and of k, j, i prints its local step and its primitive variables,
  !> and every process stops with the others (mpi_stop_all, exit code 1). NaN and infinities are found by their bits
  !> (Not_Finite), which -ffast-math and -fp-model fast=2 (the RELEASE build) cannot fold; the finite values are compared
  !> with 0. Collective (one reduction); with every local step positive and finite it changes nothing.
  subroutine Check_Local_Steps ( domain )
    use MOSE_Finite_m, only: Not_Finite
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root, mpi_rank_, mpi_size_, mpi_allreduce_min_r8, mpi_stop_all
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer(I8) :: nbad(1)
    integer  :: b, i, j, k, fb, fi, fj, fk
    real(R8) :: d, bmin
    character(len=32) :: who

    nbad = 0
    fb = 0; fi = 0; fj = 0; fk = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        d = domain%blk(b)%dtlocal(i,j,k)
        if (Not_Finite(d) .or. .not. (d > 0.0_R8)) then
          nbad(1) = nbad(1) + 1
          if (fb == 0) then
            fb = b; fi = i; fj = j; fk = k
          endif
        endif
      enddo; enddo; enddo
    enddo
    call Sum_Process ( nbad )
    if (nbad(1) == 0) return

    ! the first such cell over the processes: the smallest block (each block is on one process), then k, j, i
    call mpi_allreduce_min_r8 ( merge(real(fb, R8), huge(1.0_R8), fb > 0), bmin )
    if (mpi_is_root) write(*,'(A,ES14.7,A,I0,A,A)') ' [ERROR] [MOSE-Multirate] at the start of the macro step (time ', &
      domain%time, ' s) the local step is not a positive finite number in ', nbad(1), &
      trim(merge(' cell ', ' cells', nbad(1) == 1_I8)), '; the first one:'
    if (fb > 0 .and. fb == nint(bmin)) then
      who = ''
      if (mpi_size_ > 1) write(who,'(A,I0,A)') ' (process ', mpi_rank_, ')'
      write(*,'(A,4(1X,I0),A,ES15.7E3,A,A)') ' [ERROR]   cell b, i, j, k:', fb, fi, fj, fk, ': local step ', &
        domain%blk(fb)%dtlocal(fi,fj,fk), ' s', trim(who)
      write(*,'(A,*(1X,ES12.4E3))') ' [ERROR]   primitive variables:', domain%blk(fb)%P(:,fi,fj,fk)
    endif
    call mpi_stop_all ( 'MOSE-Multirate: a local step that is not a positive finite number at the start of the macro step' )

  end subroutine Check_Local_Steps


  !> r_post = h_c(n)/dtcap_c(t_n+1) of every cell, with the macro step H_mac and the level map of the last macro
  !> step and the dtlocal of the state it reached (computed by Compute_dt for the next macro step): largest value and
  !> number of cells above 1 (a step beyond the one admitted, cfl included) and above 1/cfl (a Courant number above 1 where
  !> the convective limit rules), each with the rounding tolerance 1e-12 of the guard on h/(dtcap/sigma): a cell that ran
  !> exactly at its limit gives 1 +/- one unit in the last place. Reads only; every rank takes part in the reductions.
  !> r_nb (Gnedin, Semenov and Kravtsov 2018, section 2, condition (4): every neighbour of a cell ends its step within the
  !> step that the cell admits): h_c(n)/dtcap_v(t_n+1) for the face neighbours v of c in its block (i+-1, j+-1, k+-1),
  !> largest value and cells c above 1 and above 1/cfl with the same tolerance, and cells with r_post or r_nb above 1.
  !> Across a block connection nothing more is needed: the two sides have the same level, so h_c = h_p and
  !> h_c/dtcap_p(t_n+1) is the r_post of the partner p. h_c/min_v dtcap_v is the largest h_c/dtcap_v (division by a
  !> positive number is monotone).
  subroutine R_Post ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    real(R8) :: h, r, rn, dnb, rp_loc, rn_loc, cnt(5), tol1, icfl
    integer(I8) :: c1, c2, c3, c4, c5
    integer  :: b, i, j, k, n1, n2, n3

    tol1 = 1.0_R8 + 1.0e-12_R8
    icfl = tol1 / obj_time_scheme%cfl
    rp_loc = 0.0_R8;  rn_loc = 0.0_R8
    c1 = 0;  c2 = 0;  c3 = 0;  c4 = 0;  c5 = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1);  n2 = domain%blk(b)%dim(2);  n3 = domain%blk(b)%dim(3)
      !$omp parallel do collapse(3) private ( h, r, rn, dnb ) reduction ( max : rp_loc, rn_loc ) &
      !$omp   reduction ( + : c1, c2, c3, c4, c5 )
      do k = 1, n3; do j = 1, n2; do i = 1, n1
        h = H_mac / 2.0_R8**mb(b)%lev(i,j,k)
        r = h / domain%blk(b)%dtlocal(i,j,k)
        ! smallest admissible step of the face neighbours of the cell in the block
        dnb = huge(1.0_R8)
        if (i > 1)  dnb = min(dnb, domain%blk(b)%dtlocal(i-1,j,k))
        if (i < n1) dnb = min(dnb, domain%blk(b)%dtlocal(i+1,j,k))
        if (j > 1)  dnb = min(dnb, domain%blk(b)%dtlocal(i,j-1,k))
        if (j < n2) dnb = min(dnb, domain%blk(b)%dtlocal(i,j+1,k))
        if (k > 1)  dnb = min(dnb, domain%blk(b)%dtlocal(i,j,k-1))
        if (k < n3) dnb = min(dnb, domain%blk(b)%dtlocal(i,j,k+1))
        rn = 0.0_R8
        if (dnb < huge(1.0_R8)) rn = h / dnb
        rp_loc = max(rp_loc, r)
        rn_loc = max(rn_loc, rn)
        if (r  > tol1) c1 = c1 + 1
        if (r  > icfl) c2 = c2 + 1
        if (rn > tol1) c3 = c3 + 1
        if (rn > icfl) c4 = c4 + 1
        if (r > tol1 .or. rn > tol1) c5 = c5 + 1
      enddo; enddo; enddo
      !$omp end parallel do
    enddo
    rpost_last = Max_R8_Process ( rp_loc )
    rnb_last = Max_R8_Process ( rn_loc )
    cnt = real([c1, c2, c3, c4, c5], R8)
    call Sum_R8_Process ( cnt )
    rpost_cnt     = cnt(1:2)
    rnb_cnt       = cnt(3:5)
    rpost_max_run = max(rpost_max_run, rpost_last)
    rnb_max_run   = max(rnb_max_run, rnb_last)
    n_rpost_eval  = n_rpost_eval + 1
    if (cnt(1) > 0.0_R8) n_rpost_gt1   = n_rpost_gt1 + 1
    if (cnt(2) > 0.0_R8) n_rpost_gtcfl = n_rpost_gtcfl + 1
    if (cnt(3) > 0.0_R8) n_rnb_gt1     = n_rnb_gt1 + 1
    if (cnt(4) > 0.0_R8) n_rnb_gtcfl   = n_rnb_gtcfl + 1
    if (cnt(5) > 0.0_R8) n_any_gt1     = n_any_gt1 + 1

  end subroutine R_Post


  !> nonlocal-cfl = count or redo: at the start tau of a
  !> sub-step of level l (stage 1), every cell v of level l against condition (4) of Gnedin et al.: hhat_v <= dcap_v(tau),
  !> with hhat_v = 2 h_l when a face neighbour one level coarser starts its sub-step at tau (tau = taus(l-1)), h_l
  !> otherwise; dcap_v from the per-cell copy of the local step (Dtcap_Cell) on the start state of the sub-step and the
  !> bound of the stiffness on that state (RANS_Stiffness_Bound). A NaN is a violation. With count it only counts; with
  !> redo it also keeps the smallest dcap of the cell seen in the attempt (dmin) and counts a dcap that is not a positive
  !> number, from which the redo rebuilds the map or halves H. At
  !> tau = t_n dcap is also compared with dtlocal of Compute_dt (mismatches over 1e-12). Run by every thread of the
  !> parallel region.
  subroutine Cert_Stage ( domain, l, t, taus )
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l
    integer(I8), intent(in) :: t, taus(0:LMAX_CODE)
    integer  :: q, b, i, j, k
    real(R8) :: hl, hh, dcap, rq, w0
    logical  :: cstart

    ! wall time of the certificate (cert_wall): read by the master thread only, with timer_wtime (MPI_WTIME with MPI,
    ! MPI_THREAD_FUNNELED)
    w0 = 0.0_R8
    !$omp master
    w0 = timer_wtime()
    !$omp end master

    hl = scale(H_mac, -l)
    cstart = .false.
    if (l > 0) cstart = (t == taus(l-1))
    !$omp do schedule(static) private(b, i, j, k, hh, dcap, rq) reduction(+:cert_nv, cert_nm, cert_nn) reduction(max:cert_rq)
    do q = 1, lv(l)%cells%n
      b = lv(l)%cells%b(q); i = lv(l)%cells%i(q); j = lv(l)%cells%j(q); k = lv(l)%cells%k(q)
      dcap = Dtcap_Cell ( domain%blk(b), i, j, k )
      hh = hl
      if (cstart) then
        if (mb(b)%lev(i-1,j,k) == l-1 .or. mb(b)%lev(i+1,j,k) == l-1 .or. mb(b)%lev(i,j-1,k) == l-1 .or. &
            mb(b)%lev(i,j+1,k) == l-1 .or. mb(b)%lev(i,j,k-1) == l-1 .or. mb(b)%lev(i,j,k+1) == l-1) hh = 2.0_R8 * hl
      endif
      ! the limit up to the rounding, the tolerance of the formula of the levels and of the guard on h/(dtcap/sigma): a cell
      ! that starts at its limit gives hh/dcap = 1 within some units in the last place, which no map and no H would change
      if (.not. (hh <= dcap * (1.0_R8 + 1.0e-12_R8))) then
        cert_nv = cert_nv + 1
        mb(b)%cflag(i,j,k) = .true.
      endif
      ! redo: the smallest admissible step of the cell seen in the attempt (a dcap that is not a positive number is
      ! counted: the redo then halves H)
      if (redo_on) then
        if (dcap > 0.0_R8 .and. dcap <= huge(1.0_R8)) then
          mb(b)%dmin(i,j,k) = min(mb(b)%dmin(i,j,k), dcap)
        else
          cert_nn = cert_nn + 1
        endif
      endif
      rq = hh / dcap
      if (.not. (rq <= huge(1.0_R8))) rq = huge(1.0_R8)
      cert_rq = max(cert_rq, rq)
      if (t == 0_I8 .and. first_attempt) then
        if (.not. (abs(dcap - domain%blk(b)%dtlocal(i,j,k)) <= 1.0e-12_R8 * domain%blk(b)%dtlocal(i,j,k))) &
          cert_nm = cert_nm + 1
      endif
    enddo
    !$omp end do
    !$omp single
    cert_viol(l)  = cert_viol(l) + cert_nv
    cert_nv       = 0
    cert_mism_step = cert_mism_step + cert_nm
    cert_nm       = 0
    cert_rq_step  = max(cert_rq_step, cert_rq)
    cert_rq       = 0.0_R8
    cert_nn_step  = cert_nn_step + cert_nn
    cert_nn       = 0
    !$omp end single
    !$omp master
    cert_wall(l) = cert_wall(l) + (timer_wtime() - w0)
    !$omp end master

  end subroutine Cert_Stage


  !> End of an attempt with nonlocal-cfl = count or redo: cells with a violation per level, sums over the processes, macro
  !> steps that the certificate would reject (with a redo, a second attempt at least), and the cells and faces that took a
  !> state that is not physical (redo, L > 0), which reject the attempt too. The counts of the certificate of an attempt
  !> with such a state stay out of the totals of the run (after it the certificate sees states that are not the solution),
  !> except the mismatches with dtlocal, which come from the state at t_n of the first attempt. Collective.
  subroutine Cert_Close ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer  :: b, i, j, k
    integer(I8) :: mm(1), mv(3)
    real(R8) :: r

    cert_cells = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        if (mb(b)%cflag(i,j,k)) then
          cert_cells(mb(b)%lev(i,j,k)) = cert_cells(mb(b)%lev(i,j,k)) + 1
          mb(b)%cflag(i,j,k) = .false.
        endif
      enddo; enddo; enddo
    enddo
    call Sum_Process ( cert_viol )
    call Sum_Process ( cert_cells )
    mm(1) = cert_mism_step
    call Sum_Process ( mm )
    cert_mism_run = cert_mism_run + mm(1)
    ! dcap that is not a positive number, cells and faces with a state that is not physical: one reduction
    mv = [cert_nn_step, ev_cells, ev_faces]
    call Sum_Process ( mv )
    cert_bad = mv(1) > 0
    cert_event = mv(2) + mv(3) > 0
    cert_ev_face = mv(3) > 0
    cert_rejected = sum(cert_viol) > 0 .or. cert_bad .or. cert_event
    r = Max_R8_Process ( cert_rq_step )
    cert_rq_step = r
    n_cert_eval = n_cert_eval + 1
    if (cert_event) then
      n_ev_rej = n_ev_rej + 1
      if (mv(2) > 0) n_ev_rej_cell = n_ev_rej_cell + 1
      if (mv(3) > 0) n_ev_rej_face = n_ev_rej_face + 1
    else
      if (sum(cert_viol) > 0) n_cert_rej = n_cert_rej + 1
      cert_viol_run  = cert_viol_run + cert_viol
      cert_cells_run = cert_cells_run + cert_cells
      cert_rq_run    = max(cert_rq_run, cert_rq_step)
    endif

  end subroutine Cert_Close


  !> Local admissible step of one cell: the internal procedure compute of Compute_dt (Mod_dt.f90) copied statement by
  !> statement (no preconditioning and no CFL ramp, which the multirate rejects), on the current state of the cell and
  !> the bound of the stiffness of the explicit turbulence sources from the state of the cell (time-accurate RANS runs,
  !> RANS_Stiffness_Bound, as Compute_dt gives it to compute), with the distance rc of the
  !> cell centre from the axis and the radial direction er of the 2-D axisymmetric wedge as Compute_dt gives them to
  !> compute. Thread-safe.
  real(R8) function Dtcap_Cell ( blk, i, j, k ) result ( dtmin )
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_RANS, only: Eddy_Viscosity, RANS_Stiffness_Bound
    use MOSE_Lib_Metrics, only: delthe
    implicit none
    type(MOSE_block_type), intent(in) :: blk
    integer, intent(in) :: i, j, k
    real(R8), parameter :: sigma_SA = 2d0/3d0
    integer  :: d, T_i, Tint(2)
    real(R8) :: rho, Rgas, Sound, gn, un, lambda, rate_c, rate_v, temperature, Tdiff, cp, cv, mil, kl, mie, Dmax, stiff
    real(R8) :: cfl, vnn, p
    real(R8) :: rhoi(nsc), vel(3), met(3,3), Dm(nsc), dummy(3,3)
    real(R8) :: rc, er(3), cen(3), nu_e, sum_x, sum_r, s_h
    real(R8) :: cos12, cos13, cos23, ssum, qprod, mroot
    logical  :: wedge

    cfl = obj_time_scheme%cfl
    vnn = obj_time_scheme%vnn
    rhoi = blk%P(1:nsc,i,j,k)
    vel  = blk%P(nu:nw,i,j,k)
    p    = blk%P(np,i,j,k)
    met  = blk%m(i,j,k)%c
    stiff = 0d0
    if ( obj_time_scheme%time_accurate .and. model == 2 ) &
      stiff = RANS_Stiffness_Bound ( blk%P(nt:nt+nRANS-1,i,j,k), sum ( blk%P(1:nsc,i,j,k) ), blk%yn(i,j,k), &
                                     blk%k_rough(i,j,k) )

    ! 2-D axisymmetric wedge: rc and er of Compute_dt (Mod_dt.f90)
    wedge = ndir == 2 .and. .not. ( delthe >= 0d0 .and. delthe <= 0d0 ) .and. model > 0
    rc = 0d0
    er = 0d0
    if ( wedge ) then
      cen = 0.125d0 * ( blk%node(i  ,j  ,k  )%c + blk%node(i-1,j  ,k  )%c + &
                        blk%node(i  ,j-1,k  )%c + blk%node(i-1,j-1,k  )%c + &
                        blk%node(i  ,j  ,k-1)%c + blk%node(i-1,j  ,k-1)%c + &
                        blk%node(i  ,j-1,k-1)%c + blk%node(i-1,j-1,k-1)%c )
      rc = sqrt ( cen(2)**2 + cen(3)**2 )
      if ( rc > 0d0 ) er = [ 0d0, cen(2), cen(3) ] / rc
    end if

    call co_rotot_Rtot ( rhoi, rho, Rgas )
    Sound = f_ss ( rhoi, p, rho, Rgas )

    rate_c = 0d0
    rate_v = 0d0
    do d = 1, ndir
      gn = norm2 ( met(d,:) )
      un = abs ( dot_product ( vel, met(d,:) ) ) / gn
      lambda = un + Sound
      rate_c = rate_c + lambda * gn
      rate_v = rate_v + gn**2
    enddo

    if ( model==0 ) then
      dtmin = cfl / rate_c
      return
    endif

    ! 3-D cells skewed in their three directions: the bound (2 + m + 1/m)/4 of the excess of the spectral radius of the
    ! viscous operator (compute of Compute_dt)
    if ( ndir == 3 ) then
      cos12 = dot_product ( met(1,:), met(2,:) ) / ( norm2 ( met(1,:) ) * norm2 ( met(2,:) ) )
      cos13 = dot_product ( met(1,:), met(3,:) ) / ( norm2 ( met(1,:) ) * norm2 ( met(3,:) ) )
      cos23 = dot_product ( met(2,:), met(3,:) ) / ( norm2 ( met(2,:) ) * norm2 ( met(3,:) ) )
      ssum  = cos12**2 + cos13**2 + cos23**2
      qprod = cos12 * cos13 * cos23
      if ( ssum + 2d0 * qprod > 1d0 ) then
        mroot  = 2d0 * sqrt ( ssum / 3d0 ) * &
                 cos ( acos ( max ( -1d0, min ( 1d0, 3d0 * sqrt ( 3d0 ) * qprod / ssum**1.5d0 ) ) ) / 3d0 )
        rate_v = rate_v * 0.25d0 * ( 2d0 + mroot + 1d0 / mroot )
      end if
    end if

    temperature = p / ( rho * Rgas )
    T_i     = int ( temperature )
    Tdiff   = temperature - T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1
    cp = f_cp_expr ( rhoi, Tint, Tdiff, rho )
    cv = cp - Rgas
    call co_k_mi_lam_Wilke_expr ( rhoi, rho, Tint, Tdiff, mil, kl )
    if ( obj_sim_param%Prl > 0d0 ) kl = mil * cp / obj_sim_param%Prl

    mie = 0.d0
    if (model==2) then
      dummy = 0d0
      call Eddy_Viscosity ( mut=mie, rans_variables=blk%P(nt:nprim,i,j,k), mul=mil, rho=rho, &
                            vel_gradient=dummy, walldist=1d-6, k_rough=0d0 )
    endif

    Dmax = max ( 4d0/3d0 * ( mil + mie ) / rho, ( kl + mie * cp / obj_rans%Prt ) / ( rho * cv ) )

    if ( nsc > 1 ) then
      if ( obj_sim_param%Sc > 0d0 ) then
        Dm = ( mil / obj_sim_param%Sc + mie / obj_rans%Sct ) / rho
      else
        call co_DS_expr ( rhoi, rho, Tint, Tdiff, p, Dm )
        Dm = Dm + mie / ( rho * obj_rans%Sct )
      endif
      Dmax = max ( Dmax, maxval ( Dm ) )
    endif

    if ( model==2 .and. nRANS==1 ) Dmax = max ( Dmax, ( mil + blk%P(nt,i,j,k) ) / ( sigma_SA * rho ) )

    ! Wedge: increase of the row sum from the stress on the wedge faces (compute of Compute_dt)
    s_h = 0d0
    if ( rc > 0d0 ) then
      nu_e  = ( mil + mie ) / rho
      sum_x = 0d0
      sum_r = 0d0
      do d = 1, ndir
        sum_x = sum_x + abs ( met(d,1) )
        sum_r = sum_r + abs ( dot_product ( met(d,:), er ) )
      enddo
      s_h = nu_e * max ( 4d0/3d0 / rc**2 + 2d0/3d0 * ( sum_x + sum_r ) / rc, 1d0 / rc**2 + sum_r / rc )
    end if

    dtmin = 1d0 / ( rate_c / cfl + ( Dmax * rate_v + 0.25d0 * s_h ) / vnn )

    if ( stiff > 0d0 ) dtmin = 1d0 / ( 1d0 / dtmin + stiff )

  end function Dtcap_Cell


  !> theta of the limiter: on cells of equal width the forward-Euler MUSCL step is TVD for a Courant number up to
  !> 1/(1 + theta); on cells of different widths, with phi and phi/r in [0, 2 theta] (all but van Albada), the bracket
  !> of Lib_Reconstruction lies in [1 - 2 theta w+_{j-1}, 1 + 2 theta w-_j], with w+_j = dl_j/(dl_j + dl_{j+1}) and
  !> w-_j = dl_j/(dl_{j-1} + dl_j) (van Albada: [w-_j - 0.207, w-_j + 1.207]). From
  !> u-_{j+1/2} - u-_{j-1/2} = D-u_j [1 + phi(r_j)/(2 r_j) - phi(r_{j-1})/2] (Lib_Reconstruction) and the incremental form of
  !> Harten, the bracket must lie in [0, 1 + theta]: phi and phi/r in [0, 2] (van Leer, MC, superbee) give theta = 1, in
  !> [0, 1] (minmod) 1/2, in [-0.207, 1.207] (van Albada) sqrt(2)/2; first order 0. Docs, Non-local CFL condition.
  !> It reads space-reconstruction and flux-limiter as input.ini gives them, which the set-up keeps (the name of the
  !> limiter in use goes to limiter_name); Initialise takes it once, at the first macro step (theta_lim).
  real(R8) function Theta_Limiter () result ( th )
    implicit none
    if (trim(obj_space_scheme%space_reconstruction) == 'first-order') then
      th = 0.0_R8
    else
      select case (trim(obj_space_scheme%flux_limiter))
      case ('minmod');    th = 0.5_R8
      case ('vanalbada'); th = sqrt(2.0_R8) / 2.0_R8
      case default;       th = 1.0_R8
      end select
    endif
  end function Theta_Limiter


  !> True when the equations of the run take the diffusion, equations = navier-stokes: the test of Assign_Setup, which
  !> gives model > 0 for this value only. It reads the key, so it serves also before Assign_Setup (point A, before
  !> Assign_Setup writes model); the key is fixed for the run (KEY_FIXED), so it is the same at every reading of input.ini.
  logical function Viscous_Keys () result ( v )
    implicit none
    v = obj_sim_param%simulation_type == 'navier-stokes'
  end function Viscous_Keys


  !> True when the keys c of [MOSE-Multirate] make the certificate of condition (4) redo the attempts of a macro step
  !> with levels: max-level > 0 and the condition of redo_on (Initialise), nonlocal-cfl = redo with levels = recompute
  !> and no fixed-macro-step. The keys of the input at set-up (obj_multirate), those of the set-up at a reading (cfg0).
  logical function Redo_Keys ( c ) result ( r )
    implicit none
    type(multirate_t), intent(in) :: c
    r = c%max_level > 0 .and. trim(c%nonlocal_cfl) == 'redo' .and. c%fixed_macro_step <= 0.0_R8 .and. &
        trim(c%levels) == 'recompute'
  end function Redo_Keys


  !> H_cap of the current map from blk%dtlocal: min over the cells of dtcap_c 2**lev_c/sigma and, with the interface margin,
  !> of min_v dtcap_v 2**lev_c/kappa over the face neighbours v in the block (all processes). Collective.
  real(R8) function Hcap_Of_Map ( domain ) result ( Hc )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer  :: b, i, j, k
    real(R8) :: hl
    hl = huge(1.0_R8)
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        hl = min(hl, domain%blk(b)%dtlocal(i,j,k) * 2.0_R8**mb(b)%lev(i,j,k) / cfg0%safety)
        if (margin_on) hl = min(hl, Dtcap_Neighbours(domain%blk(b), i, j, k) * 2.0_R8**mb(b)%lev(i,j,k) / kappa_m)
      enddo; enddo; enddo
    enddo
    Hc = Min_R8_Process ( hl )
  end function Hcap_Of_Map


  !> Rejected attempt (nonlocal-cfl = redo): the state at the start of the macro step back in P (mb%P0, copied at the
  !> start of every attempt), the flux accumulators of the faces between levels zeroed; the time is set by the caller.
  !> Run by every thread of the parallel region.
  subroutine Restore_State ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, i, j, k
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(2) schedule(static)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          domain%blk(b)%P(:,i,j,k)  = mb(b)%P0(:,i,j,k)
          domain%blk(b)%PO(:,i,j,k) = mb(b)%P0(:,i,j,k)
        enddo
      enddo; enddo
      !$omp end do
    enddo
    !$omp single
    if (allocated(acc))   acc   = 0.0_R8
    if (allocated(gsave)) gsave = 0.0_R8
    !$omp end single
  end subroutine Restore_State


  !> Smallest dtlocal of the face neighbours of a cell in its block (huge without any). Interface margin.
  real(R8) function Dtcap_Neighbours ( blk, i, j, k ) result ( d )
    implicit none
    type(MOSE_block_type), intent(in) :: blk
    integer, intent(in) :: i, j, k
    d = huge(1.0_R8)
    if (i > 1)           d = min(d, blk%dtlocal(i-1,j,k))
    if (i < blk%dim(1))  d = min(d, blk%dtlocal(i+1,j,k))
    if (j > 1)           d = min(d, blk%dtlocal(i,j-1,k))
    if (j < blk%dim(2))  d = min(d, blk%dtlocal(i,j+1,k))
    if (k > 1)           d = min(d, blk%dtlocal(i,j,k-1))
    if (k < blk%dim(3))  d = min(d, blk%dtlocal(i,j,k+1))
  end function Dtcap_Neighbours


  !> End of the run (rank 0), the lines of the macro steps of the run, each only when its situation occurred in the run, with
  !> its count and its largest value (15 digits): a macro step that started with a cell above its limit, h/(dtcap/sigma)
  !> > 1 + 1e-12 (with fixed-macro-step, or where max-level levels are not enough); a macro step with a cell above the step
  !> admitted on the state it reached, r_post > 1, or above that of a face neighbour, r_nb > 1; an attempt that violated
  !> condition (4), or a mismatch of the kernel of the certificate with dtlocal (nonlocal-cfl = count or redo); a macro step
  !> redone, a landed macro step whose H a redo halved, an attempt rejected by force-redo, an attempt rejected for a state
  !> that is not physical (redo); a cell raised by the interface margin. The counts of the violations leave out the attempts
  !> rejected for a state that is not physical. A run in which none of them occurred prints none of these lines. With the
  !> certificate the wall times per level of the certificate and of the stages come in every run, as a report of the cost.
  !> A run on several processes ends with the minimum, mean and maximum over the processes, per level, of the wall time
  !> of the stages without the time of their exchanges (stage_wall - xch_wall) and of the time of the exchanges (xch_wall):
  !> every process takes part in their reductions, before process 0 alone writes the lines.
  !> Called by every process from Explicit_Step when the run ends.
  subroutine Multirate_Report ()
    use MOSE_Mod_MPI, only: mpi_is_root, mpi_size_
    implicit none
    integer :: l, lm, n
    real(R8), allocatable :: v(:), vmin(:), vmax(:), vsum(:)

    ! Several processes: the times per level over the processes, reduced by every process before the return below.
    ! lm, the deepest level of the run, is the same on every process (Lu is the largest over the processes, Build_Lists)
    n = 0
    if (initialised .and. mpi_size_ > 1) then
      lm = min(lu_run, LMAX_CODE)
      n = lm + 1
      allocate( v(2*n), vmin(2*n), vmax(2*n), vsum(2*n) )
      v(1:n) = stage_wall(0:lm) - xch_wall(0:lm)
      v(n+1:2*n) = xch_wall(0:lm)
      call Range_R8_Process ( v, vmin, vmax, vsum )
    endif
    if (.not. initialised .or. .not. mpi_is_root) return
    if (n_ratio_gt > 0) &
      write(*,'(A,ES22.15,A,I0,A,I0,A)') '   Multirate run: max h/(dtcap/sigma) at the start of the macro steps ', &
        max_cap_ratio, ' (macro steps above 1 + 1e-12: ', n_ratio_gt, ' of ', n_macro, ')'
    if (n_rpost_gt1 > 0) &
      write(*,'(A,ES22.15,A,I0,A,I0,A,I0,A)') '   Multirate run: max r_post = h/dtcap(t_n+1) ', rpost_max_run, &
        ' (macro steps with a cell above 1: ', n_rpost_gt1, '; above 1/cfl: ', n_rpost_gtcfl, '; evaluated: ', n_rpost_eval, ')'
    if (n_rnb_gt1 > 0) &
      write(*,'(A,ES22.15,A,I0,A,I0,A,I0,A,I0,A)') '   Multirate run: max r_nb = h/dtcap(t_n+1) of a face neighbour ', &
        rnb_max_run, ' (macro steps with a cell above 1: ', n_rnb_gt1, '; above 1/cfl: ', n_rnb_gtcfl, &
        '; with r_post or r_nb above 1: ', n_any_gt1, '; evaluated: ', n_rpost_eval, ')'
    ! certificate (nonlocal-cfl = count or redo) and interface margin, per level up to the deepest level of the run
    lm = min(lu_run, LMAX_CODE)
    if (cert_on) then
      if (n_cert_rej > 0 .or. cert_mism_run > 0) then
        write(*,'(A,A,A,I0,A,I0,A,ES22.15,A,I0)') '   Multirate run: certificate (', trim(cfg0%nonlocal_cfl), &
          '): attempts with a violation ', n_cert_rej, ' of ', n_cert_eval, '; largest hhat/dcap ', cert_rq_run, &
          '; mismatches with dtlocal at t_n: ', cert_mism_run
        write(*,'(A,*(I14))') '   Multirate run: certificate, violations per level, attempts without a state that is not '// &
          'physical', (cert_viol_run(l), l = 0, lm)
        write(*,'(A,*(I14))') '   Multirate run: certificate, cells with a violation per level, summed over the attempts '// &
          'without a state that is not physical', (cert_cells_run(l), l = 0, lm)
      endif
      write(*,'(A,*(ES12.5))') '   Multirate run: certificate wall time per level [s]', (cert_wall(l), l = 0, lm)
      write(*,'(A,*(ES12.5))') '   Multirate run: stage wall time per level [s], all attempts', (stage_wall(l), l = 0, lm)
      if (redo_on .and. n_redone > 0) then
        write(*,'(A,I0,A,I0,A,I0,A,F10.6,A,I0,A,I0,A,I0)') '   Multirate run: redo, macro steps ', n_macro, ', attempts ', &
          n_att_total, ', macro steps redone ', n_redone, ', mean attempts per macro step ', &
          real(n_att_total, R8) / real(max(n_macro, 1_I8), R8), ', largest ', att_max, ', halvings of H ', n_halved, &
          ', maps made deeper by a redo ', n_deeper
        write(*,'(A,5(1X,I0))') '   Multirate run: redo, macro steps with 1, 2, 3, 4, 5 or more attempts', att_hist
      endif
      if (redo_on .and. n_unlanded > 0) &
        write(*,'(A,I0)') '   Multirate run: redo, macro steps cut to land (output instant or time-threshold) whose H '// &
          'a redo halved ', n_unlanded
      if (redo_on .and. cfg0%force_redo .and. n_forced > 0) &
        write(*,'(A,I0)') '   Multirate run: redo, attempts rejected by force-redo (test device) ', n_forced
      if (redo_on .and. n_ev_rej > 0) write(*,'(A,I0,A,I0,A,I0,A)') '   Multirate run: redo, attempts rejected for a '// &
        'state that is not physical ', n_ev_rej, ' (with a cell ', n_ev_rej_cell, ', with a face ', n_ev_rej_face, ')'
    endif
    if (margin_on .and. n_rings > 0) then
      if (sum(ring_cells_run(0:lm)) > 0) then
        write(*,'(A,ES12.5,A,I0,A)') '   Multirate run: interface margin, kappa ', kappa_m, ' (', n_rings, ' maps)'
        write(*,'(A,*(ES12.5))') '   Multirate run: cells raised by kappa per level, mean over the maps', &
          (real(ring_cells_run(l), R8) / real(n_rings, R8), l = 0, lm)
      endif
    endif
    ! several processes: the times per level over the processes (reduced above)
    if (n > 0) then
      write(*,'(A,*(ES12.5))') '   Multirate run: stage wall time per level without the exchanges [s], min over the '// &
        'processes', (vmin(l), l = 1, n)
      write(*,'(A,*(ES12.5))') '   Multirate run: stage wall time per level without the exchanges [s], mean over the '// &
        'processes', (vsum(l) / real(mpi_size_, R8), l = 1, n)
      write(*,'(A,*(ES12.5))') '   Multirate run: stage wall time per level without the exchanges [s], max over the '// &
        'processes', (vmax(l), l = 1, n)
      write(*,'(A,*(ES12.5))') '   Multirate run: exchange time per level [s], min over the processes', (vmin(l), l = n+1, 2*n)
      write(*,'(A,*(ES12.5))') '   Multirate run: exchange time per level [s], mean over the processes', &
        (vsum(l) / real(mpi_size_, R8), l = n+1, 2*n)
      write(*,'(A,*(ES12.5))') '   Multirate run: exchange time per level [s], max over the processes', (vmax(l), l = n+1, 2*n)
    endif

  end subroutine Multirate_Report


  !> First macro step: set-up checks (point B), per-block arrays, ghost-record map, keys frozen.
  subroutine Initialise ( domain )
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root, mpi_size_
    use MOSE_Input_Registry, only: Key_Is_Set
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b, r, g, n1, n2, n3, Ig, Jg, Kg

    call Check_Domain ( domain )
    cfg0 = obj_multirate

    allocate( mb(domain%nb) )
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      allocate( mb(b)%lev (1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc) )
      allocate( mb(b)%grec(1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc) )
      allocate( mb(b)%ifc (3, 1:n1, 1:n2, 1:n3) )
      allocate( mb(b)%P0  (nprim, 1:n1, 1:n2, 1:n3) )
      mb(b)%lev  = -1
      mb(b)%grec = 0
      mb(b)%ifc  = 0
      ! no shock detector: the flux weight is 1 everywhere (Zero_Residuals sets it at every single-rate stage)
      domain%blk(b)%beta = 1.0_R8
    enddo

    ! ghost position -> record that fills it (depths 1 ... gc along the face normal)
    do r = 1, domain % nbound
      b = domain%bc(r)%b
      if (.not. is_local_block(b)) cycle
      do g = 1, gc
        Ig = domain%bc(r)%i - guide(domain%bc(r)%f,1)*g
        Jg = domain%bc(r)%j - guide(domain%bc(r)%f,2)*g
        Kg = domain%bc(r)%k - guide(domain%bc(r)%f,3)*g
        mb(b)%grec(Ig,Jg,Kg) = r
      enddo
    enddo

    allocate( stage_wall(0:LMAX_CODE), xch_wall(0:LMAX_CODE) )
    stage_wall = 0.0_R8
    xch_wall = 0.0_R8
    ! several processes: the start of the exchanges of the module is timed (Exchange_P, _F, _G, _Masks, _Partner_Levels),
    ! and the Timing block gives it; every process comes here at the first macro step, so they all reduce it at the end
    if (mpi_size_ > 1) call timer_xstart_enable()

    ! interface margin and certificate: theta of the limiter from the incremental form of the
    ! forward-Euler MUSCL step (u-_{j+1/2} - u-_{j-1/2} = D-u_j [1 + phi(r_j)/(2 r_j) - phi(r_{j-1})/2] in [0, 1+theta]):
    ! van Leer, MC, superbee 1 (phi, phi/r in [0, 2]); minmod 1/2; van Albada sqrt(2)/2; first order 0
    theta_lim = Theta_Limiter()
    kappa_m = max(obj_time_scheme%cfl * (1.0_R8 + theta_lim), KAPPA_FLOOR)
    cert_on = trim(cfg0%nonlocal_cfl) /= 'off'
    margin_on = (cfg0%interface_margin .or. trim(cfg0%nonlocal_cfl) == 'redo') .and. cfg0%max_level > 0
    ! redo: the certificate imposes condition (4) by rejecting and redoing the macro step; with a fixed macro step or a
    ! fixed map it only counts (verification devices: there (4) is not imposed, docs Scope)
    redo_on = trim(cfg0%nonlocal_cfl) == 'redo' .and. cfg0%fixed_macro_step <= 0.0_R8 .and. trim(cfg0%levels) == 'recompute'
    if (trim(cfg0%nonlocal_cfl) == 'redo' .and. .not. redo_on .and. mpi_is_root) write(*,'(A)') &
      ' [WARNING] [MOSE-Multirate] nonlocal-cfl = redo with fixed-macro-step or levels = fixed: the certificate only '// &
      'counts (condition (4) is not imposed)'
    if (cfg0%force_redo .and. .not. redo_on .and. mpi_is_root) write(*,'(A)') &
      ' [WARNING] [MOSE-Multirate] force-redo acts only with nonlocal-cfl = redo, levels = recompute and no '// &
      'fixed-macro-step: ignored'
    ! max-attempts written in input.ini (a run that leaves it out takes its default, which the warning does not concern)
    if (.not. redo_on .and. Key_Is_Set(trim(codename)//'-Multirate', 'max-attempts') .and. mpi_is_root) write(*,'(A)') &
      ' [WARNING] [MOSE-Multirate] max-attempts acts only with nonlocal-cfl = redo, levels = recompute and no '// &
      'fixed-macro-step: ignored'
    allocate( cert_viol(0:LMAX_CODE), cert_cells(0:LMAX_CODE), cert_viol_run(0:LMAX_CODE), cert_cells_run(0:LMAX_CODE), &
              ring_cells(0:LMAX_CODE), ring_cells_run(0:LMAX_CODE), cert_wall(0:LMAX_CODE) )
    cert_viol = 0; cert_cells = 0; cert_viol_run = 0; cert_cells_run = 0; ring_cells = 0; ring_cells_run = 0
    cert_wall = 0.0_R8
    allocate( ev_rec%v(nprim, 6) )
    ev_rec%v = 0.0_R8
    if (cert_on) then
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        allocate( mb(b)%cflag(domain%blk(b)%dim(1), domain%blk(b)%dim(2), domain%blk(b)%dim(3)) )
        mb(b)%cflag = .false.
        if (redo_on) then
          allocate( mb(b)%dmin(domain%blk(b)%dim(1), domain%blk(b)%dim(2), domain%blk(b)%dim(3)) )
          mb(b)%dmin = huge(1.0_R8)
        endif
      enddo
    endif
    if (mpi_is_root .and. (cert_on .or. margin_on)) write(*,'(A,F8.5,A,F8.5,A,L2,A,A,A,L2)') &
      '   Multirate: theta of the limiter ', theta_lim, ', kappa = max(cfl (1 + theta), 1.1) = ', kappa_m, &
      ', interface margin', margin_on, ', nonlocal-cfl = ', trim(cfg0%nonlocal_cfl), ', redo', redo_on

    call Build_Peers ( domain )
    initialised = .true.

  end subroutine Initialise


  !> P0 := P (start of the macro step) or PO := P0 (end of the macro step) on all local cells.
  subroutine Copy_P0 ( domain, to_p0 )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    logical, intent(in) :: to_p0
    integer :: b, i, j, k
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(2) schedule(static)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          if (to_p0) then
            mb(b)%P0(:,i,j,k) = domain%blk(b)%P(:,i,j,k)
          else
            domain%blk(b)%PO(:,i,j,k) = mb(b)%P0(:,i,j,k)
          endif
        enddo
      enddo; enddo
      !$omp end do
    enddo
  end subroutine Copy_P0


  ! =================================================================================================================
  !  RECURSION (executed by every thread of the parallel region; work sharing inside)
  ! =================================================================================================================

  !> nst: the calls of Stage in the attempt, counted here before every call (Stage returns at once on a process without
  !> work for its level): the same count on every thread and every process, the number of the stage in the record of a
  !> state that is not physical.
  recursive subroutine Advance ( domain, l, taus, nst )
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l
    integer(I8), intent(in) :: taus(0:LMAX_CODE)
    integer(I8), intent(inout) :: nst
    integer(I8) :: tl(0:LMAX_CODE), h

    tl = taus
    h = ishft(1_I8, Lu + 1 - l)                ! ticks of h_l
    call Start_Substep ( domain, l )
    nst = nst + 1_I8
    call Stage ( domain, l, 1, tl(l), tl, nst )
    if (l < Lu) then
      call Stash ( domain, l )
      tl(l+1) = tl(l)
      call Advance ( domain, l+1, tl, nst )
      tl(l+1) = tl(l) + h / 2_I8
      call Advance ( domain, l+1, tl, nst )
      call Restore ( domain, l )
    endif
    nst = nst + 1_I8
    call Stage ( domain, l, 2, tl(l) + h, tl, nst )

  end subroutine Advance


  !> PO := P on the cells of level l (start state of their sub-step).
  subroutine Start_Substep ( domain, l )
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l
    integer :: q
    !$omp do schedule(static)
    do q = 1, lv(l)%cells%n
      associate( b => lv(l)%cells%b(q), i => lv(l)%cells%i(q), j => lv(l)%cells%j(q), k => lv(l)%cells%k(q) )
        domain%blk(b)%PO(:,i,j,k) = domain%blk(b)%P(:,i,j,k)
      end associate
    enddo
    !$omp end do
  end subroutine Start_Substep


  !> Predictor of the level-l cells read by finer stages: primitive copy (restored after the finer sub-steps) and
  !> conservative start state and predictor (the two ends of the interpolant).
  subroutine Stash ( domain, l )
    use FLINT_Lib_Thermodynamic, only: prim2cons
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l
    integer :: q
    !$omp do schedule(static)
    do q = 1, lv(l)%used%n
      associate( b => lv(l)%used%b(q), i => lv(l)%used%i(q), j => lv(l)%used%j(q), k => lv(l)%used%k(q) )
        lv(l)%Ws(:,q) = domain%blk(b)%P(:,i,j,k)
        lv(l)%Us(1:np,q) = prim2cons( domain%blk(b)%P(1:np,i,j,k) )
        lv(l)%Ub(1:np,q) = prim2cons( domain%blk(b)%PO(1:np,i,j,k) )
        if (nprim > np) then
          lv(l)%Us(np+1:nprim,q) = domain%blk(b)%P (np+1:nprim,i,j,k)
          lv(l)%Ub(np+1:nprim,q) = domain%blk(b)%PO(np+1:nprim,i,j,k)
        endif
      end associate
    enddo
    !$omp end do
  end subroutine Stash


  subroutine Restore ( domain, l )
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l
    integer :: q
    !$omp do schedule(static)
    do q = 1, lv(l)%used%n
      associate( b => lv(l)%used%b(q), i => lv(l)%used%i(q), j => lv(l)%used%j(q), k => lv(l)%used%k(q) )
        domain%blk(b)%P(:,i,j,k) = lv(l)%Ws(:,q)
      end associate
    enddo
    !$omp end do
  end subroutine Restore


  !> Stage k (1: predictor, 2: corrector) of level l at tick t, the stage number s of the attempt (Advance).
  subroutine Stage ( domain, l, kst, t, taus, s )
    use MOSE_Lib_RANS,  only: RANS_Source_Terms
    use MOSE_Lib_Diffusive, only: Visc_State
    use MOSE_Lib_BC_Fluxes_Extrapolation, only: Report_400_Entering_Lines
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l, kst
    integer(I8), intent(in) :: t, taus(0:LMAX_CODE), s
    ! Local
    integer :: m, q, d, g
    integer(I8) :: num, den
    real(R8) :: hl, w0
    logical :: fz                                  ! frozen, as this thread read it for the stage

    ! a process without cells of level l still runs the stage when it sends to other processes what their stages of
    ! level l read (interpolated cells, Pg)
    if (lv(l)%cells%n == 0 .and. xl(l)%nsP + xl(l)%nrP + xl(l)%nsG + xl(l)%nrG + xl(l)%nsF + xl(l)%nrF == 0) return

    !$omp single
    domain % time = t_start + H_mac * ( real(t, R8) / real(ticks_H, R8) )
    if (trace_now) call Trace_Line ( l, kst, t, taus )
    ! a state that is not physical (ev_on): the stage, for its record; a process that took one in an earlier stage of the
    ! attempt computes no face of (4) and updates no cell in (6) from here on (the rest of the stage, exchanges included,
    ! goes on, so the processes stay together); the faces of a stage at tick 0 see only the state at t_n, which a shorter H
    ! does not change, so a face that is not physical there stops the run, as in the single-rate step. Fixed here, before
    ! the first barrier of the stage, so every thread reads the same values for the whole stage
    cur_s = s; cur_l = l; cur_k = kst; cur_time = domain % time
    frozen   = ev_seen
    defer_on = ev_on .and. t /= 0_I8
    !$omp end single
    ! each thread reads frozen once, past the barrier of the single: with frozen the stage can end without a barrier after
    ! (5), and the single of the next stage writes frozen again
    fz = frozen
    ! wall time of the stage (stage_wall): read by the master thread only, with timer_wtime (MPI_WTIME with MPI,
    ! MPI_THREAD_FUNNELED)
    w0 = 0.0_R8
    !$omp master
    w0 = timer_wtime()
    !$omp end master

    ! (1) coarser cells read by the stage: interpolant at the stage time
    do m = 0, l - 1
      if (lv(l)%halo(m)%n == 0) cycle
      num = t - taus(m)
      den = ishft(1_I8, Lu + 1 - m)
      !$omp do schedule(static)
      do q = 1, lv(l)%halo(m)%n
        call Interpolate_Cell ( domain, lv(l)%halo(m), q, m, num, den )
      enddo
      !$omp end do
    enddo

    ! (2) ghost records read by the stage, at the stage time
    !$omp do schedule(dynamic, 64)
    do q = 1, lv(l)%ngh
      call Ghost_Record ( domain, lv(l)%ghr(q) )
    enddo
    !$omp end do

    ! (2m) several processes: the ghost cells of the records whose partner block is on another process, from it, at the
    !      stage time (master thread: MPI_THREAD_FUNNELED)
    if (xl(l)%nsP + xl(l)%nrP > 0) then
      !$omp master
      call Exchange_P ( domain, l )
      !$omp end master
      !$omp barrier
    endif

    ! (2b) viscous runs: tangential ghost values Pg(:,3:6) of the connection records of the level-l cells, read by
    !      BC_Connection_Visc, from the partner cells and ghost cells refreshed in (1), (2) and (2m) (Fill_Ghost_Cell: after
    !      the ghost records and the MPI exchange); a partner block on another process: from that process, which also
    !      computes here those of its records whose partner block is here (Exchange_G)
    if (lv(l)%npg > 0) then
      !$omp do schedule(dynamic, 64)
      do q = 1, lv(l)%npg
        if (rpeer(lv(l)%pgr(q)) == 0) call Pg_Record ( domain, lv(l)%pgr(q) )
      enddo
      !$omp end do
    endif
    if (xl(l)%nsG + xl(l)%nrG > 0) then
      !$omp master
      call Exchange_G ( domain, l )
      !$omp end master
      !$omp barrier
    endif

    ! (2v) viscous runs: the viscous variables blk%Visc (Visc_State) of every cell and ghost cell that the viscous fluxes
    !      of the stage read (Visc_Of_Level), from their primitive variables as (1), (2) and (2m) left them: the level-l
    !      cells, the coarser cells interpolated at the stage time, the finer cells as they are, the ghost cells at the
    !      stage time; as Viscous_Variables after Fill_Ghost_Cell in the single-rate stage
    if (lv(l)%visc%n > 0) then
      !$omp do schedule(static)
      do q = 1, lv(l)%visc%n
        associate( b => lv(l)%visc%b(q), i => lv(l)%visc%i(q), j => lv(l)%visc%j(q), k => lv(l)%visc%k(q) )
          call Visc_State ( domain%blk(b)%P(:,i,j,k), domain%blk(b)%Visc(:,i,j,k) )
        end associate
      enddo
      !$omp end do
    endif

    ! (2t) trace (first macro step, viscous runs): the viscous variables that the fluxes of the stage read, against
    !      those of their primitive variables at this point (one thread)
    if (trace_now .and. model > 0) then
      !$omp single
      call Trace_Visc ( domain, l, kst, t )
      !$omp end single
    endif

    ! (3) residuals of the level-l cells
    !$omp do schedule(static)
    do q = 1, lv(l)%cells%n
      domain%blk(lv(l)%cells%b(q))%R(:,lv(l)%cells%i(q),lv(l)%cells%j(q),lv(l)%cells%k(q)) = 0.0_R8
    enddo
    !$omp end do

    ! (4) internal faces with a level-l cell, direction by direction, in each segment the convective then (viscous runs)
    !     the diffusive faces: every cell gets its contributions in the order of Fluxes_blk (none on a process that took a
    !     state that is not physical in an earlier stage of the attempt: frozen)
    if (.not. fz) then
      do d = 1, 3
        !$omp do schedule(static)
        do q = 1, lv(l)%seg(d)%n
          call Face_Segment ( domain, l, kst, d, q )
        enddo
        !$omp end do
      enddo
    endif

    ! (5a) connection faces of the level-l cells: one flux per face, computed at this stage by the record that owns it
    !      for the two cells (Connection_Fluxes at the start of BC_Fluxes); the shares of the records whose owner record
    !      is on another process come from it (Exchange_F, master thread)
    if (lv(l)%ncn > 0) then
      !$omp do schedule(dynamic, 64)
      do q = 1, lv(l)%ncn
        call Connection_Record ( domain, lv(l)%cnr(q) )
      enddo
      !$omp end do
    endif
    if (xl(l)%nsF + xl(l)%nrF > 0) then
      !$omp master
      call Exchange_F ( domain, l )
      !$omp end master
      !$omp barrier
    endif

    ! (5) boundary records of the level-l cells, group by group (block, face) as in BC_Fluxes
    do g = 1, lv(l)%ngrp
      !$omp do schedule(dynamic, 64)
      do q = lv(l)%grp(g), lv(l)%grp(g+1) - 1
        call BC_Flux_Record ( domain, lv(l)%bcr(q) )
      enddo
      !$omp end do
    enddo

    ! (5b) the lines of the gas entering through 400 faces at this stage, as at the end of BC_Fluxes: every thread calls
    !      Report_400_Entering_Lines past the implicit barrier of the last worksharing loop of (5), before the sources and
    !      the update of the cells (a stage without boundary records of level l evaluates no 400 face)
    call Report_400_Entering_Lines ()

    ! (5c) RANS runs: the explicit turbulence sources of the level-l cells (RANS_Source_Terms), in the order of the
    !      single-rate stage: the sources of a cell read its own ghost cells along the axes, those of the fluxes, set in
    !      (2) of this stage (Ghost_Wall at the walls: one set of ghost values for the fluxes and the sources).
    !      One level (every cell at level 0): the block routine of the single-rate stage itself, whose worksharing
    !      loops run in this parallel region, so the macro step is the single-rate step to the bit (the per-cell copies,
    !      compiled in another context, may differ from them in the last bit)
    if (rans_kind /= 0) then
      if (Lu == 0) then
        call RANS_Source_Terms ( domain )
      else
        !$omp do schedule(static)
        do q = 1, lv(l)%cells%n
          call Source_Cell ( domain%blk(lv(l)%cells%b(q)), lv(l)%cells%i(q), lv(l)%cells%j(q), lv(l)%cells%k(q) )
        enddo
        !$omp end do
      endif
    endif

    ! (5e) nonlocal-cfl = count or redo: certificate of condition (4) at the start of the sub-step, on the
    !      state of stage 1 (the start state) and the bound of the stiffness on it
    if (cert_on .and. kst == 1 .and. Lu > 0) call Cert_Stage ( domain, l, t, taus )

    ! (5d) planar quasi-2D and 1-D cases with symmetry k faces (freeze_w) and 1-D cases with symmetry j faces (freeze_v)
    !      (Setup_Quasi2D_W): the z and the y momentum are not integrated, as Zero_W_Residual and Zero_V_Residual do
    !      before the update of the single-rate stage
    if (freeze_w) then
      !$omp do schedule(static)
      do q = 1, lv(l)%cells%n
        domain%blk(lv(l)%cells%b(q))%R(nw,lv(l)%cells%i(q),lv(l)%cells%j(q),lv(l)%cells%k(q)) = 0.0_R8
      enddo
      !$omp end do
    endif
    if (freeze_v) then
      !$omp do schedule(static)
      do q = 1, lv(l)%cells%n
        domain%blk(lv(l)%cells%b(q))%R(nv,lv(l)%cells%i(q),lv(l)%cells%j(q),lv(l)%cells%k(q)) = 0.0_R8
      enddo
      !$omp end do
    endif

    ! (6) Heun update of the level-l cells with h_l (none on a frozen process: its cells keep their states)
    hl = scale(H_mac, -l)
    if (.not. fz) then
      !$omp do schedule(static)
      do q = 1, lv(l)%cells%n
        associate( b => lv(l)%cells%b(q), i => lv(l)%cells%i(q), j => lv(l)%cells%j(q), k => lv(l)%cells%k(q) )
          call Compute_Cons_NoIRS( domain%blk(b)%P(:,i,j,k),    &
                                   domain%blk(b)%PO(:,i,j,k),   &
                                   domain%blk(b)%r(:,i,j,k),    &
                                   hl,                          &
                                   domain%blk(b)%vol(i,j,k),    &
                                   kst, obj_time_scheme%n_RK, b, i, j, k, obj_chemistry%strangcoeff )
        end associate
      enddo
      !$omp end do
    endif

    !$omp master
    stage_wall(l) = stage_wall(l) + (timer_wtime() - w0)
    !$omp end master

  end subroutine Stage


  !> Coarser cell q of level m: P := interpolant at theta = num/den between its start state and its predictor, in
  !> conservative variables; exact copies at theta = 0 (start state) and theta = 1 (predictor).
  subroutine Interpolate_Cell ( domain, hl, q, m, num, den )
    use FLINT_Lib_Thermodynamic, only: cons2prim, co_rotot_Rtot
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    type(halo_list_t), intent(in) :: hl
    integer, intent(in) :: q, m
    integer(I8), intent(in) :: num, den
    integer :: u
    real(R8) :: th, Uc(nprim), rho, Rgas, temperature
    associate( b => hl%b(q), i => hl%i(q), j => hl%j(q), k => hl%k(q) )
      u = hl%u(q)
      if (num == 0_I8) then
        domain%blk(b)%P(:,i,j,k) = domain%blk(b)%PO(:,i,j,k)
      elseif (num == den) then
        domain%blk(b)%P(:,i,j,k) = lv(m)%Ws(:,u)
      else
        th = real(num, R8) / real(den, R8)
        Uc = lv(m)%Ub(:,u) + th * ( lv(m)%Us(:,u) - lv(m)%Ub(:,u) )
        call co_rotot_Rtot( lv(m)%Ws(1:nsc,u), rho, Rgas )
        temperature = lv(m)%Ws(np,u) / ( rho * Rgas )
        domain%blk(b)%P(1:np,i,j,k) = cons2prim( Uc(1:np), temperature )
        if (nprim > np) domain%blk(b)%P(np+1:nprim,i,j,k) = Uc(np+1:nprim)
      endif
    end associate
  end subroutine Interpolate_Cell


  !> Faces a ... e of segment q of direction d (level l, stage kst), convective faces, then (viscous runs) diffusive
  !> faces (Diffusive_Segment). Faces between two cells of the level: the call of Fluxes_blk, unchanged. Faces between
  !> levels: the flux of the face is computed on a zeroed two-cell residual and only the level-l cell receives it; the
  !> coarse side of the correction uses the accumulators.
  subroutine Face_Segment ( domain, l, kst, d, q )
    use MOSE_Lib_Convective, only: Convective_Flux
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l, kst, d, q
    ! Local
    integer :: b, i, j, k, id
    real(R8) :: dl4(-1:2)
    integer  :: fev(2)                    ! a face that is not physical: the error flags of its two states (Face_Event)
    real(R8) :: fst(nprim,2)              !   and the two states

    fev = 0
    b = lv(l)%seg(d)%b(q)
    associate( blk => domain%blk(b) )
    select case (d)
    case (1)
      j = lv(l)%seg(d)%p(q); k = lv(l)%seg(d)%q(q)
      do i = lv(l)%seg(d)%a(q), lv(l)%seg(d)%e(q)
        id = mb(b)%ifc(1,i,j,k)
        dl4(-1) = blk % dl(i-1,j,k) % c(1) ; dl4(0) = blk % dl(i  ,j,k) % c(1)
        dl4( 1) = blk % dl(i+1,j,k) % c(1) ; dl4(2) = blk % dl(i+2,j,k) % c(1)
        if (id == 0) then
          call Convective_Flux ( dl4,                          &
                                 blk % dir(1) % f(i,j,k) % N,  &
                                 blk % dir(1) % f(i,j,k) % A,  &
                                 blk % P(:,i-1:i+2,j,k),       &
                                 blk % R(:,i:i+1,j,k),         &
                                 blk % beta(i,j,k),            &
                                 blk % Ur(i,j,k),              &
                                 blk % Ur(i+1,j,k),            &
                                 [b, 1, i, j, k], fev, fst )
          if (fev(1) /= 0 .or. fev(2) /= 0) call Face_Event ( blk % P(:,i-1:i+2,j,k), [b, 1, i, j, k], fev, fst )
        else
          call Interface_Face ( blk, l, kst, id, dl4, blk % dir(1) % f(i,j,k) % N, blk % dir(1) % f(i,j,k) % A, &
                                blk % P(:,i-1:i+2,j,k), blk % beta(i,j,k), blk % Ur(i,j,k), blk % Ur(i+1,j,k), &
                                mb(b)%lev(i,j,k), mb(b)%lev(i+1,j,k), i, j, k, i+1, j, k, b, 1 )
        endif
      enddo
    case (2)
      i = lv(l)%seg(d)%p(q); k = lv(l)%seg(d)%q(q)
      do j = lv(l)%seg(d)%a(q), lv(l)%seg(d)%e(q)
        id = mb(b)%ifc(2,i,j,k)
        dl4(-1) = blk % dl(i,j-1,k) % c(2) ; dl4(0) = blk % dl(i,j  ,k) % c(2)
        dl4( 1) = blk % dl(i,j+1,k) % c(2) ; dl4(2) = blk % dl(i,j+2,k) % c(2)
        if (id == 0) then
          call Convective_Flux ( dl4,                          &
                                 blk % dir(2) % f(i,j,k) % N,  &
                                 blk % dir(2) % f(i,j,k) % A,  &
                                 blk % P(:,i,j-1:j+2,k),       &
                                 blk % R(:,i,j:j+1,k),         &
                                 blk % beta(i,j,k),            &
                                 blk % Ur(i,j,k),              &
                                 blk % Ur(i,j+1,k),            &
                                 [b, 2, i, j, k], fev, fst )
          if (fev(1) /= 0 .or. fev(2) /= 0) call Face_Event ( blk % P(:,i,j-1:j+2,k), [b, 2, i, j, k], fev, fst )
        else
          call Interface_Face ( blk, l, kst, id, dl4, blk % dir(2) % f(i,j,k) % N, blk % dir(2) % f(i,j,k) % A, &
                                blk % P(:,i,j-1:j+2,k), blk % beta(i,j,k), blk % Ur(i,j,k), blk % Ur(i,j+1,k), &
                                mb(b)%lev(i,j,k), mb(b)%lev(i,j+1,k), i, j, k, i, j+1, k, b, 2 )
        endif
      enddo
    case (3)
      i = lv(l)%seg(d)%p(q); j = lv(l)%seg(d)%q(q)
      do k = lv(l)%seg(d)%a(q), lv(l)%seg(d)%e(q)
        id = mb(b)%ifc(3,i,j,k)
        dl4(-1) = blk % dl(i,j,k-1) % c(3) ; dl4(0) = blk % dl(i,j,k  ) % c(3)
        dl4( 1) = blk % dl(i,j,k+1) % c(3) ; dl4(2) = blk % dl(i,j,k+2) % c(3)
        if (id == 0) then
          call Convective_Flux ( dl4,                          &
                                 blk % dir(3) % f(i,j,k) % N,  &
                                 blk % dir(3) % f(i,j,k) % A,  &
                                 blk % P(:,i,j,k-1:k+2),       &
                                 blk % R(:,i,j,k:k+1),         &
                                 blk % beta(i,j,k),            &
                                 blk % Ur(i,j,k),              &
                                 blk % Ur(i,j,k+1),            &
                                 [b, 3, i, j, k], fev, fst )
          if (fev(1) /= 0 .or. fev(2) /= 0) call Face_Event ( blk % P(:,i,j,k-1:k+2), [b, 3, i, j, k], fev, fst )
        else
          call Interface_Face ( blk, l, kst, id, dl4, blk % dir(3) % f(i,j,k) % N, blk % dir(3) % f(i,j,k) % A, &
                                blk % P(:,i,j,k-1:k+2), blk % beta(i,j,k), blk % Ur(i,j,k), blk % Ur(i,j,k+1), &
                                mb(b)%lev(i,j,k), mb(b)%lev(i,j,k+1), i, j, k, i, j, k+1, b, 3 )
        endif
      enddo
    end select
    end associate

    if (model > 0) call Diffusive_Segment ( domain % blk(b), l, kst, d, q )

  end subroutine Face_Segment


  !> Convective flux of a face between two levels, seen from a stage of level l. Side 0 is the cell (i0,j0,k0), side 1
  !> the cell (i1,j1,k1); Rt(:,s) is what the face adds to the residual of side s. The correction: Interface_Begin and
  !> Interface_End, with the accumulators of the convective flux.
  subroutine Interface_Face ( blk, l, kst, id, dl4, normal, area, Prim, beta, Ur0, Ur1, lev0, lev1, &
                              i0, j0, k0, i1, j1, k1, b, d )
    use MOSE_Lib_Convective, only: Convective_Flux
    implicit none
    type(MOSE_block_type), intent(inout) :: blk
    integer, intent(in) :: l, kst, id, lev0, lev1, i0, j0, k0, i1, j1, k1
    integer, intent(in) :: b, d                    ! block and direction of the face: error messages of Convective_Flux
    real(R8), intent(in) :: dl4(-1:2), normal(3), area, beta, Ur0, Ur1
    real(R8), intent(in) :: Prim(nprim,-1:2)
    ! Local
    real(R8) :: Rt(nprim,0:1), fst(nprim,2)
    integer :: act, fev(2)
    logical :: coarse, compute

    call Interface_Begin ( blk, l, kst, id, FX_CONV, lev0, lev1, i0, j0, k0, i1, j1, k1, act, coarse, compute )
    if (.not. compute) return
    Rt = 0.0_R8
    fev = 0
    call Convective_Flux ( dl4, normal, area, Prim, Rt, beta, Ur0, Ur1, [b, d, i0, j0, k0], fev, fst )
    ! a face that is not physical gives no flux: Rt stays zero for the cell and for the correction
    if (fev(1) /= 0 .or. fev(2) /= 0) call Face_Event ( Prim, [b, d, i0, j0, k0], fev, fst )
    call Interface_End ( blk, id, FX_CONV, act, coarse, Rt, i0, j0, k0, i1, j1, k1 )

  end subroutine Interface_Face


  !> A face of the stage with a reconstructed state that is not physical: Convective_Flux added no flux and returned the
  !> error flags fev of check_gas_state and the two states fst; Prim is the stencil, loc the face (block, direction, i, j,
  !> k). In an attempt that the redo can redo, at a stage after tick 0 (defer_on), the face is counted and recorded and the
  !> attempt will be rejected; otherwise the run stops with the lines of the face, as the single-rate step does. fev back
  !> to zero for the next face.
  subroutine Face_Event ( Prim, loc, fev, fst )
    use MOSE_Lib_Convective, only: Interface_State_Lines
    implicit none
    real(R8), intent(in) :: Prim(nprim,-1:2), fst(nprim,2)
    integer, intent(in) :: loc(5)
    integer, intent(inout) :: fev(2)
    integer(I8) :: key(7)

    if (.not. defer_on) then
      call Interface_State_Lines ( fev(1), fev(2), fst(:,1), fst(:,2), Prim, loc )
      ! stop 1, not mpi_abort_all, as Convective_Flux: this runs inside OpenMP threads and MPI is initialised FUNNELED
      stop 1
    endif
    !$omp atomic
    ev_faces = ev_faces + 1
    key = [cur_s, 0_I8, int(loc(1), I8), int(loc(2), I8), int(loc(5), I8), int(loc(4), I8), int(loc(3), I8)]
    !$omp critical (mr_event)
    ev_seen = .true.
    if (Key_Before(key, ev_rec%key)) then
      ev_rec%key   = key
      ev_rec%kind  = 3
      ev_rec%lev   = cur_l
      ev_rec%kst   = cur_k
      ev_rec%time  = cur_time
      ev_rec%flags = fev
      ev_rec%v(:,1:2) = fst
      ev_rec%v(:,3:6) = Prim
    endif
    !$omp end critical (mr_event)
    fev = 0

  end subroutine Face_Event


  !> A cell of the stage whose update took a state that is not physical (Compute_Cons_NoIRS): vals the state (kind 1: the
  !> conservative variables; kind 2: the primitive variables), counted and recorded; the attempt will be rejected.
  subroutine Cell_Event ( vals, kind, b, i, j, k )
    implicit none
    real(R8), intent(in) :: vals(nprim)
    integer, intent(in) :: kind, b, i, j, k
    integer(I8) :: key(7)

    !$omp atomic
    ev_cells = ev_cells + 1
    key = [cur_s, 1_I8, int(b, I8), 0_I8, int(k, I8), int(j, I8), int(i, I8)]
    !$omp critical (mr_event)
    ev_seen = .true.
    if (Key_Before(key, ev_rec%key)) then
      ev_rec%key   = key
      ev_rec%kind  = kind
      ev_rec%lev   = cur_l
      ev_rec%kst   = cur_k
      ev_rec%time  = cur_time
      ev_rec%flags = 0
      ev_rec%v(:,1) = vals
    endif
    !$omp end critical (mr_event)

  end subroutine Cell_Event


  !> Key a before key b, element by element (the first that differs decides).
  logical function Key_Before ( a, b ) result ( before )
    implicit none
    integer(I8), intent(in) :: a(7), b(7)
    integer :: n
    before = .false.
    do n = 1, 7
      if (a(n) /= b(n)) then
        before = a(n) < b(n)
        return
      endif
    enddo
  end function Key_Before


  !> Conservative correction of a face between two levels, flux of kind fx (convective or diffusive), seen from a stage
  !> of level l: act = side of the level-l cell (0: (i0,j0,k0), 1: (i1,j1,k1)), coarse = the other cell is finer.
  !>   level l fine:   the fine cell receives its flux; acc += 1/4 of the coarse side (weight b*h_l/h_(l-1) = 1/2*1/2)
  !>   level l coarse: stage 1: own flux, saved; stage 2: 2*acc - saved (then acc = saved = 0)
  !> With conservative-correction = .false. both sides use their own fluxes (positive control).
  !> Begin: on stage 2 of the coarse side the cell receives 2*acc - saved and compute = .false.; otherwise compute =
  !> .true. and the caller computes the flux of the face on a zeroed residual Rt(:,0:1), then calls Interface_End.
  subroutine Interface_Begin ( blk, l, kst, id, fx, lev0, lev1, i0, j0, k0, i1, j1, k1, act, coarse, compute )
    implicit none
    type(MOSE_block_type), intent(inout) :: blk
    integer, intent(in)  :: l, kst, id, fx, lev0, lev1, i0, j0, k0, i1, j1, k1
    integer, intent(out) :: act
    logical, intent(out) :: coarse, compute
    integer :: other

    if (lev0 == l) then
      act = 0; other = lev1
    else
      act = 1; other = lev0
    endif
    coarse = other > l
    compute = .true.

    if (coarse .and. kst == 2 .and. cfg0%conservative_correction) then
      if (act == 0) then
        blk%R(:,i0,j0,k0) = blk%R(:,i0,j0,k0) + ( 2.0_R8 * acc(:,id,fx) - gsave(:,id,fx) )
      else
        blk%R(:,i1,j1,k1) = blk%R(:,i1,j1,k1) + ( 2.0_R8 * acc(:,id,fx) - gsave(:,id,fx) )
      endif
      acc(:,id,fx)   = 0.0_R8
      gsave(:,id,fx) = 0.0_R8
      compute = .false.
    endif

  end subroutine Interface_Begin


  !> End of the correction of a face between two levels (see Interface_Begin): the level-l cell receives its side of
  !> the flux Rt; the coarse side saves it (coarse stage 1) or the fine side accumulates 1/4 of the coarse side.
  subroutine Interface_End ( blk, id, fx, act, coarse, Rt, i0, j0, k0, i1, j1, k1 )
    implicit none
    type(MOSE_block_type), intent(inout) :: blk
    integer, intent(in)  :: id, fx, act, i0, j0, k0, i1, j1, k1
    logical, intent(in)  :: coarse
    real(R8), intent(in) :: Rt(nprim,0:1)

    if (act == 0) then
      blk%R(:,i0,j0,k0) = blk%R(:,i0,j0,k0) + Rt(:,0)
    else
      blk%R(:,i1,j1,k1) = blk%R(:,i1,j1,k1) + Rt(:,1)
    endif
    if (.not. cfg0%conservative_correction) return
    if (coarse) then
      gsave(:,id,fx) = Rt(:,act)                              ! stage 1 of the coarse cell
    else
      acc(:,id,fx) = acc(:,id,fx) + 0.25_R8 * Rt(:,1-act)     ! fine stage: impulse seen by the coarse cell
    endif

  end subroutine Interface_End


  !> Diffusive faces a ... e of segment q of direction d (level l, stage kst; viscous runs), after the convective faces
  !> of the segment. Faces between two cells of the level: the call of Fluxes_blk, unchanged (face metric Mfp, as there:
  !> with diffusive-metric = stencil the metric stored for the face, read in place, with mean the one of
  !> Viscous_Face_Metric from the two cell metrics; with a wall ghost cell in the stencil, the metrics per group of
  !> quantities Mg). Faces between levels: the flux is computed on a zeroed two-cell residual, only the level-l cell
  !> receives it, and the correction has its own accumulators (kind FX_DIFF). Stencil of a face (c, c+e_d): the
  !> primitive variables of the two cells, and the viscous variables blk%Visc, which the stage refreshed ((2v) of
  !> Stage), of the two cells and of their neighbours +-1 in the two other directions.
  subroutine Diffusive_Segment ( blk, l, kst, d, q )
    use MOSE_Lib_Diffusive, only: Diffusive_Flux
    use MOSE_Lib_Stencil_Metrics, only: Viscous_Face_Metric, METRIC_STENCIL
    implicit none
    type(MOSE_block_type), intent(inout), target :: blk
    integer, intent(in) :: l, kst, d, q
    ! Local
    real(R8) :: Sc, Sct, Prt, Prl, Rt(nprim,0:1)
    real(R8), target :: Mf(3,3)                  ! face metric computed here (diffusive-metric other than stencil)
    real(R8), pointer, contiguous :: Mfp(:,:)    ! face metric handed to Diffusive_Flux: the stored one (stencil) or Mf
    real(R8), pointer, contiguous :: Mg(:,:,:)   ! metrics per group with a wall ghost cell in the stencil, else absent
    logical  :: soot_enabled, coarse, compute
    integer  :: b, i, j, k, id, act, metric

    ! the arguments of Internal_Fluxes and Fluxes_blk
    Sc  = obj_sim_param%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt
    Prl = obj_sim_param%Prl
    soot_enabled = obj_soot%enabled
    metric = obj_space_scheme%metric_mode

    b = lv(l)%seg(d)%b(q)
    select case (d)
    case (1)
      j = lv(l)%seg(d)%p(q); k = lv(l)%seg(d)%q(q)
      do i = lv(l)%seg(d)%a(q), lv(l)%seg(d)%e(q)
        id = mb(b)%ifc(1,i,j,k)
        if (id == 0) then
          nullify ( Mg )
          if ( metric == METRIC_STENCIL ) then
            Mfp => blk % dir(1) % Mf(i,j,k) % c
            if ( blk % dir(1) % ig(i,j,k) > 0 ) Mg => blk % dir(1) % Mg(:,:,:,blk % dir(1) % ig(i,j,k))
          else
            Mf = Viscous_Face_Metric ( blk, 1, i, j, k, metric )
            Mfp => Mf
          end if
          call Diffusive_Flux ( blk % dir(1) % f(i,j,k) % N,  &
                                blk % dir(1) % f(i,j,k) % A,  &
                                blk % yn(i  ,j,k),            &
                                blk % yn(i+1,j,k),            &
                                blk % k_rough(i  ,j,k),       &
                                blk % k_rough(i+1,j,k),       &
                                blk % P(:,i  ,j,k),           &
                                blk % P(:,i+1,j,k),           &
                                blk % Visc(:,i  ,j,k),        &
                                blk % Visc(:,i+1,j,k),        &
                                blk % Visc(:,i  ,j-1,k),      &
                                blk % Visc(:,i  ,j+1,k),      &
                                blk % Visc(:,i+1,j-1,k),      &
                                blk % Visc(:,i+1,j+1,k),      &
                                blk % Visc(:,i  ,j,k-1),      &
                                blk % Visc(:,i  ,j,k+1),      &
                                blk % Visc(:,i+1,j,k-1),      &
                                blk % Visc(:,i+1,j,k+1),      &
                                Mfp,                          &
                                blk % R(:,i  ,j,k),           &
                                blk % R(:,i+1,j,k),           &
                                1, 2, 3,                &
                                Sc, Sct, Prt, Prl, soot_enabled, Mg = Mg )
        else
          call Interface_Begin ( blk, l, kst, id, FX_DIFF, mb(b)%lev(i,j,k), mb(b)%lev(i+1,j,k), i, j, k, i+1, j, k, &
                                 act, coarse, compute )
          if (.not. compute) cycle
          Rt = 0.0_R8
          nullify ( Mg )
          if ( metric == METRIC_STENCIL ) then
            Mfp => blk % dir(1) % Mf(i,j,k) % c
            if ( blk % dir(1) % ig(i,j,k) > 0 ) Mg => blk % dir(1) % Mg(:,:,:,blk % dir(1) % ig(i,j,k))
          else
            Mf = Viscous_Face_Metric ( blk, 1, i, j, k, metric )
            Mfp => Mf
          end if
          call Diffusive_Flux ( blk % dir(1) % f(i,j,k) % N,  &
                                blk % dir(1) % f(i,j,k) % A,  &
                                blk % yn(i  ,j,k),            &
                                blk % yn(i+1,j,k),            &
                                blk % k_rough(i  ,j,k),       &
                                blk % k_rough(i+1,j,k),       &
                                blk % P(:,i  ,j,k),           &
                                blk % P(:,i+1,j,k),           &
                                blk % Visc(:,i  ,j,k),        &
                                blk % Visc(:,i+1,j,k),        &
                                blk % Visc(:,i  ,j-1,k),      &
                                blk % Visc(:,i  ,j+1,k),      &
                                blk % Visc(:,i+1,j-1,k),      &
                                blk % Visc(:,i+1,j+1,k),      &
                                blk % Visc(:,i  ,j,k-1),      &
                                blk % Visc(:,i  ,j,k+1),      &
                                blk % Visc(:,i+1,j,k-1),      &
                                blk % Visc(:,i+1,j,k+1),      &
                                Mfp,                          &
                                Rt(:,0),                      &
                                Rt(:,1),                      &
                                1, 2, 3,                &
                                Sc, Sct, Prt, Prl, soot_enabled, Mg = Mg )
          call Interface_End ( blk, id, FX_DIFF, act, coarse, Rt, i, j, k, i+1, j, k )
        endif
      enddo
    case (2)
      i = lv(l)%seg(d)%p(q); k = lv(l)%seg(d)%q(q)
      do j = lv(l)%seg(d)%a(q), lv(l)%seg(d)%e(q)
        id = mb(b)%ifc(2,i,j,k)
        if (id == 0) then
          nullify ( Mg )
          if ( metric == METRIC_STENCIL ) then
            Mfp => blk % dir(2) % Mf(i,j,k) % c
            if ( blk % dir(2) % ig(i,j,k) > 0 ) Mg => blk % dir(2) % Mg(:,:,:,blk % dir(2) % ig(i,j,k))
          else
            Mf = Viscous_Face_Metric ( blk, 2, i, j, k, metric )
            Mfp => Mf
          end if
          call Diffusive_Flux ( blk % dir(2) % f(i,j,k) % N,  &
                                blk % dir(2) % f(i,j,k) % A,  &
                                blk % yn(i,j  ,k),            &
                                blk % yn(i,j+1,k),            &
                                blk % k_rough(i,j  ,k),       &
                                blk % k_rough(i,j+1,k),       &
                                blk % P(:,i,j  ,k),           &
                                blk % P(:,i,j+1,k),           &
                                blk % Visc(:,i,j  ,k),        &
                                blk % Visc(:,i,j+1,k),        &
                                blk % Visc(:,i-1,j  ,k),      &
                                blk % Visc(:,i+1,j  ,k),      &
                                blk % Visc(:,i-1,j+1,k),      &
                                blk % Visc(:,i+1,j+1,k),      &
                                blk % Visc(:,i,j  ,k-1),      &
                                blk % Visc(:,i,j  ,k+1),      &
                                blk % Visc(:,i,j+1,k-1),      &
                                blk % Visc(:,i,j+1,k+1),      &
                                Mfp,                          &
                                blk % R(:,i,j  ,k),           &
                                blk % R(:,i,j+1,k),           &
                                2, 1, 3,                &
                                Sc, Sct, Prt, Prl, soot_enabled, Mg = Mg )
        else
          call Interface_Begin ( blk, l, kst, id, FX_DIFF, mb(b)%lev(i,j,k), mb(b)%lev(i,j+1,k), i, j, k, i, j+1, k, &
                                 act, coarse, compute )
          if (.not. compute) cycle
          Rt = 0.0_R8
          nullify ( Mg )
          if ( metric == METRIC_STENCIL ) then
            Mfp => blk % dir(2) % Mf(i,j,k) % c
            if ( blk % dir(2) % ig(i,j,k) > 0 ) Mg => blk % dir(2) % Mg(:,:,:,blk % dir(2) % ig(i,j,k))
          else
            Mf = Viscous_Face_Metric ( blk, 2, i, j, k, metric )
            Mfp => Mf
          end if
          call Diffusive_Flux ( blk % dir(2) % f(i,j,k) % N,  &
                                blk % dir(2) % f(i,j,k) % A,  &
                                blk % yn(i,j  ,k),            &
                                blk % yn(i,j+1,k),            &
                                blk % k_rough(i,j  ,k),       &
                                blk % k_rough(i,j+1,k),       &
                                blk % P(:,i,j  ,k),           &
                                blk % P(:,i,j+1,k),           &
                                blk % Visc(:,i,j  ,k),        &
                                blk % Visc(:,i,j+1,k),        &
                                blk % Visc(:,i-1,j  ,k),      &
                                blk % Visc(:,i+1,j  ,k),      &
                                blk % Visc(:,i-1,j+1,k),      &
                                blk % Visc(:,i+1,j+1,k),      &
                                blk % Visc(:,i,j  ,k-1),      &
                                blk % Visc(:,i,j  ,k+1),      &
                                blk % Visc(:,i,j+1,k-1),      &
                                blk % Visc(:,i,j+1,k+1),      &
                                Mfp,                          &
                                Rt(:,0),                      &
                                Rt(:,1),                      &
                                2, 1, 3,                &
                                Sc, Sct, Prt, Prl, soot_enabled, Mg = Mg )
          call Interface_End ( blk, id, FX_DIFF, act, coarse, Rt, i, j, k, i, j+1, k )
        endif
      enddo
    case (3)
      i = lv(l)%seg(d)%p(q); j = lv(l)%seg(d)%q(q)
      do k = lv(l)%seg(d)%a(q), lv(l)%seg(d)%e(q)
        id = mb(b)%ifc(3,i,j,k)
        if (id == 0) then
          nullify ( Mg )
          if ( metric == METRIC_STENCIL ) then
            Mfp => blk % dir(3) % Mf(i,j,k) % c
            if ( blk % dir(3) % ig(i,j,k) > 0 ) Mg => blk % dir(3) % Mg(:,:,:,blk % dir(3) % ig(i,j,k))
          else
            Mf = Viscous_Face_Metric ( blk, 3, i, j, k, metric )
            Mfp => Mf
          end if
          call Diffusive_Flux ( blk % dir(3) % f(i,j,k) % N,  &
                                blk % dir(3) % f(i,j,k) % A,  &
                                blk % yn(i,j,k  ),            &
                                blk % yn(i,j,k+1),            &
                                blk % k_rough(i,j,k  ),       &
                                blk % k_rough(i,j,k+1),       &
                                blk % P(:,i,j,k  ),           &
                                blk % P(:,i,j,k+1),           &
                                blk % Visc(:,i,j,k  ),        &
                                blk % Visc(:,i,j,k+1),        &
                                blk % Visc(:,i-1,j,k  ),      &
                                blk % Visc(:,i+1,j,k  ),      &
                                blk % Visc(:,i-1,j,k+1),      &
                                blk % Visc(:,i+1,j,k+1),      &
                                blk % Visc(:,i,j-1,k  ),      &
                                blk % Visc(:,i,j+1,k  ),      &
                                blk % Visc(:,i,j-1,k+1),      &
                                blk % Visc(:,i,j+1,k+1),      &
                                Mfp,                          &
                                blk % R(:,i,j,k  ),           &
                                blk % R(:,i,j,k+1),           &
                                3, 1, 2,                &
                                Sc, Sct, Prt, Prl, soot_enabled, Mg = Mg )
        else
          call Interface_Begin ( blk, l, kst, id, FX_DIFF, mb(b)%lev(i,j,k), mb(b)%lev(i,j,k+1), i, j, k, i, j, k+1, &
                                 act, coarse, compute )
          if (.not. compute) cycle
          Rt = 0.0_R8
          nullify ( Mg )
          if ( metric == METRIC_STENCIL ) then
            Mfp => blk % dir(3) % Mf(i,j,k) % c
            if ( blk % dir(3) % ig(i,j,k) > 0 ) Mg => blk % dir(3) % Mg(:,:,:,blk % dir(3) % ig(i,j,k))
          else
            Mf = Viscous_Face_Metric ( blk, 3, i, j, k, metric )
            Mfp => Mf
          end if
          call Diffusive_Flux ( blk % dir(3) % f(i,j,k) % N,  &
                                blk % dir(3) % f(i,j,k) % A,  &
                                blk % yn(i,j,k  ),            &
                                blk % yn(i,j,k+1),            &
                                blk % k_rough(i,j,k  ),       &
                                blk % k_rough(i,j,k+1),       &
                                blk % P(:,i,j,k  ),           &
                                blk % P(:,i,j,k+1),           &
                                blk % Visc(:,i,j,k  ),        &
                                blk % Visc(:,i,j,k+1),        &
                                blk % Visc(:,i-1,j,k  ),      &
                                blk % Visc(:,i+1,j,k  ),      &
                                blk % Visc(:,i-1,j,k+1),      &
                                blk % Visc(:,i+1,j,k+1),      &
                                blk % Visc(:,i,j-1,k  ),      &
                                blk % Visc(:,i,j+1,k  ),      &
                                blk % Visc(:,i,j-1,k+1),      &
                                blk % Visc(:,i,j+1,k+1),      &
                                Mfp,                          &
                                Rt(:,0),                      &
                                Rt(:,1),                      &
                                3, 1, 2,                &
                                Sc, Sct, Prt, Prl, soot_enabled, Mg = Mg )
          call Interface_End ( blk, id, FX_DIFF, act, coarse, Rt, i, j, k, i, j, k+1 )
        endif
      enddo
    end select

  end subroutine Diffusive_Segment


  !> Boundary flux of record r: the body of the select of BC_Fluxes (Mod_BC_Fluxes.f90) for the types the multirate
  !> accepts (the others are rejected at set-up). The rotating frame is rejected, so the wall branches take the
  !> non-rotating form. Connections (101, 201) between blocks: the share of the cell in the one flux of the face, which
  !> the owner record of the face computed at this stage (Connection_Record, as Connection_Fluxes of BC_Fluxes) from the
  !> ghost cells that the stages of the level refresh from the partner block at the stage time (both sides of a
  !> connection have the same level, Build_Map); 402 and 410 read domain%time, the stage time of the level; an SRM grain
  !> in its inflow mode (a 404 with srm_inflow) takes its mass flux and static temperature at the pressure of the boundary
  !> cell in the state of the stage.
  subroutine BC_Flux_Record ( domain, r )
    use MOSE_Lib_BC_Fluxes_Connection
    use MOSE_Lib_BC_Fluxes_Rotational
    use MOSE_Lib_BC_Fluxes_Symmetry
    use MOSE_Lib_BC_Fluxes_Inflow
    use MOSE_Lib_BC_Fluxes_Outflow
    use MOSE_Lib_BC_Fluxes_FarField
    use MOSE_Lib_BC_Fluxes_Inflow_Nozzle
    use MOSE_Lib_BC_Fluxes_Extrapolation
    use MOSE_Lib_BC_Fluxes_Wall_Heat
    use MOSE_Lib_BC_Fluxes_Wall_Temperature
    use MOSE_Lib_BC_Fluxes_GSI_SRM, only: SRM_Static_Temperature
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: r
    ! Local
    integer  :: Bm, Im, Jm, Km, Fm, error
    logical  :: RSM, soot_enabled
    real(R8) :: Sc, Sct, Prt, Prl
    real(R8) :: T0, g                                ! SRM grain of mode 2: static temperature and mass flux of the face
    real(R8) :: dR(nprim)                            ! share of a face in the residual of its cell

    RSM = obj_rans%RSM
    Sc  = obj_sim_param%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt
    Prl = obj_sim_param%Prl
    soot_enabled = obj_soot%enabled

    Bm = domain % bc(r) % b
    Im = domain % bc(r) % i
    Jm = domain % bc(r) % j
    Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    associate( bc => domain % bc(r) )
    select case ( bc % type )

      case (101,201) ! block connection, periodic connection: the share of the cell in the flux of the face
        domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + bc % Fc(:,1)
        if (model>0) &
          domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + bc % Fc(:,2)

      case (200) ! periodic (axisymmetric wedge)
        call BC_Rotational_Periodic_Eul ( Im, Jm, Km, Fm, domain % blk(Bm) )
        if (model>0) &
          call BC_Rotational_Periodic_Visc ( Im, Jm, Km, Fm, domain % blk(Bm), bc % rM, &
                                             Sc, Sct, Prt, Prl, soot_enabled )

      case (300) ! symmetry
        call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
        if (model>0) &
          call BC_Symmetry_Visc ( Im, Jm, Km, Fm, domain % blk(Bm), RSM )

      case (301) ! wall: prescribed heat flux
        call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
        if (model>0) &
          call BC_Wall_Heat ( Im, Jm, Km, Fm, domain % blk(Bm), bc % qw, k_rough = bc % k_rough, rec = bc )

      case (302) ! wall: prescribed temperature
        call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
        if (model>0) &
          call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), bc % Tw, k_rough = bc % k_rough, rec = bc )

      case (400) ! extrapolation
        call BC_Extrapolation ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

      case (401) ! inlet: T0 + p0
        call BC_Inlet_StagnCond ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), bc % T0, bc % p0, &
                                  bc % alpha, bc % beta, bc % ci(nrans+1:nrans+nsc), bc % ci(1:nrans), &
                                  bc % mdot, error )
        if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

      case (402) ! inlet: T0 + time-varying p0 (at the stage time)
        if (bc % p0time % exists) bc % p0 = bc % p0time % update(domain % time)
        call BC_Inlet_StagnCond ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), bc % T0, bc % p0, &
                                  bc % alpha, bc % beta, bc % ci(nrans+1:nrans+nsc), bc % ci(1:nrans), &
                                  bc % mdot, error )
        if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

      case (403) ! inlet: T0 + mass flux
        call BC_Inlet_MassFlux_T0 ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), bc % T0, bc % mdot, bc % rel_fac, &
                                    bc % alpha, bc % beta, bc % ci(nrans+1:nrans+nsc), bc % ci(1:nrans), error )
        if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

      case (404) ! inlet: T + mass flux; an SRM grain of mode 2 (srm_inflow)
        if ( bc % srm_inflow ) then
          g = bc % aCoeff * ( domain % blk(Bm) % P(np,Im,Jm,Km) / bc % pRef ) ** bc % n
          bc % g_srm = g   ! the gas mass flux of this face at this pressure, for the couplers (hydra-MI2)
          T0 = SRM_Static_Temperature ( bc % T0, g, domain % blk(Bm) % P(np,Im,Jm,Km), bc % ci(nrans+1:nrans+nsc) )
          call BC_Inlet_MassFlux_T ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), T0, g, 1d0, bc % alpha, bc % beta, &
                                     bc % ci(nrans+1:nrans+nsc), bc % ci(1:nrans), error, always_in = .true. )
        else
          call BC_Inlet_MassFlux_T ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), bc % T0, bc % mdot, bc % rel_fac, &
                                     bc % alpha, bc % beta, bc % ci(nrans+1:nrans+nsc), bc % ci(1:nrans), error )
        end if
        if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

      case (405) ! inlet: supersonic, Mach + T + p
        call BC_Inlet_Supersonic_Static ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), bc % mach, bc % T0, bc % pamb, &
                                          bc % rel_fac, bc % alpha, bc % beta, bc % ci(nrans+1:nrans+nsc), &
                                          bc % ci(1:nrans), bc % mdot, error )
        if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

      case (406) ! outlet
        call BC_Outflow ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), bc % pAmb, bc % rel_fac, bc % mdot, error )
        if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

      case (407) ! ambient: T0 + p0 + back pressure
        call BC_FarField ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), bc % T0, bc % p0, bc % pamb, &
                           bc % alpha, bc % beta, bc % ci(nrans+1:nrans+nsc), bc % ci(1:nrans), bc % mdot )

      case (408) ! inlet: T + normal velocity
        call BC_Inlet_un_T ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), bc % T0, bc % un, bc % rel_fac, &
                             bc % alpha, bc % beta, bc % ci(nrans+1:nrans+nsc), bc % ci(1:nrans), error )
        if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

      case (410) ! mapped Q2D state (ghost cells at the stage time)
        call BC_Connection_Eul ( Im, Jm, Km, Fm, domain % blk(Bm), dR )
        domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + dR

      case (420) ! choked nozzle
        call BC_Inflow_Nozzle ( Im, Jm, Km, Fm, domain % blk(Bm), bc % T0, bc % p0, bc % rel_fac, bc % psub, &
                                bc % p_shock, bc % mdot, bc % p_design, bc % rho_design, &
                                bc % ci(nrans+1:nrans+nsc), bc % ci(1:nrans), error )
        if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

      case (0)
        ! no boundary flux, as in BC_Fluxes (0 null)

      case default
        ! Check_Domain rejects every other type at the first macro step
        error stop 'MOSE-Multirate: boundary type without a multirate boundary flux'

    end select

    ! Inflow, outflow and extrapolation faces: Reynolds stresses (or -2/3 rho k with k-energy-coupling) of the momentum
    ! and energy equations, as after the select of BC_Fluxes
    if ( model == 2 .and. ( obj_rans%RSM .or. obj_rans%k_energy_coupling ) ) then
      select case ( bc % type )
        case (400:408, 410, 420)
          call BC_Open_Turbulent_Stress ( Im, Jm, Km, Fm, domain % blk(Bm) )
      end select
    end if
    end associate

  end subroutine BC_Flux_Record


  !> One flux for the connection face of the owner record i (101, 201), for the two cells of the face: the body of the
  !> loop of Connection_Fluxes (Mod_BC_Fluxes.f90), unchanged (the record i owns the face, its partner record p gets
  !> the flux with the opposite sign; with the partner block on another process Exchange_F sends that share). Called at
  !> the stages of the level of the two cells (both sides of a connection have the same level).
  subroutine Connection_Record ( domain, i )
    use MOSE_Lib_BC_Fluxes_Connection
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: i
    ! Local
    integer  :: p, Bm
    logical  :: soot_enabled
    real(R8) :: Sc, Sct, Prt, Prl, dR(nprim)

    Sc  = obj_sim_param%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt
    Prl = obj_sim_param%Prl
    soot_enabled = obj_soot%enabled

      Bm = domain % bc(i) % b
      p  = domain % bc(i) % partner
      call BC_Connection_Eul ( domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k, domain % bc(i) % f, &
                               domain % blk(Bm), dR )
      domain % bc(i) % Fc(:,1) = dR
      domain % bc(p) % Fc(:,1) = -dR
      if (model>0) then
        if ( obj_space_scheme % stencil_metric ) then
          call BC_Connection_Visc_Stencil ( domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k, domain % bc(i) % f, &
                                            domain % blk(Bm), domain % bc(i) % Pg, Sc, Sct, Prt, Prl, soot_enabled, dR )
        else
          call BC_Connection_Visc ( domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k, domain % bc(i) % f, &
                                    domain % blk(Bm), domain % bc(i) % Mg(1), domain % bc(i) % Pg, &
                                    Sc, Sct, Prt, Prl, soot_enabled, dR )
        end if
        domain % bc(i) % Fc(:,2) = dR
        domain % bc(p) % Fc(:,2) = -dR
      end if

  end subroutine Connection_Record


  !> Ghost cells of record r: the select of the local-record loop of Fill_Ghost_Cell (Lib_Ghost.f90), every branch,
  !> also those of the boundary types that Check_Domain rejects (the walls coupled to another solver, 103 and 104, the
  !> SRM grain 502 and the gas-surface walls). Fill_BC_Ghost_Q2D (Pg(:,3:6) and the edge ghost cells of the 410 records)
  !> is not called: no flux of a stage reads them (BC_Fluxes gives a 410 face its inviscid flux only).
  subroutine Ghost_Record ( domain, r )
    use MOSE_Lib_Ghost, only: Ghost_Connection, Ghost_Axysimmetry, Ghost_Symmetry, Ghost_Q2D_Mapped, &
                              Ghost_ZG_Extrapolate, Ghost_Wall, Ghost_Extrapolate, Ghost_SRM_Turbulence
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: r
    integer :: Bm, Im, Jm, Km, Fm, Bs

    Bm = domain % bc(r) % b
    Im = domain % bc(r) % i
    Jm = domain % bc(r) % j
    Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    select case ( domain % bc(r) % type )
      case (101,201)
        Bs = domain % bc(r) % bs
        call Ghost_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                domain % bc(r) % is, domain % bc(r) % js, domain % bc(r) % ks, domain % bc(r) % fs, &
                                domain % blk(Bs), domain % bc(r) % Pg )
      case (200)
        call Ghost_Axysimmetry ( Im, Jm, Km, Fm, domain % blk(Bm) )
      case (300)
        call Ghost_Symmetry ( Im, Jm, Km, Fm, domain % blk(Bm) )
      case (410)
        call Ghost_Q2D_Mapped ( domain % bc(r), domain % blk(Bm), domain % time )
      case (0,401:407,420)
        call Ghost_ZG_Extrapolate ( Im, Jm, Km, Fm, domain % blk(Bm) )
      case (103,104,301,302,400,503:506)
        call Ghost_Wall ( domain % bc(r), domain % blk(Bm) )
      case (502)
        call Ghost_Extrapolate ( Im, Jm, Km, Fm, domain % blk(Bm) )
        if ( model == 2 ) call Ghost_SRM_Turbulence ( domain % bc(r), domain % blk(Bm) )
      case (102)
        error stop 'MOSE-Multirate: chimera ghost cells are not supported'
      case default
        call Ghost_Extrapolate ( Im, Jm, Km, Fm, domain % blk(Bm) )
    end select
  end subroutine Ghost_Record


  !> Tangential ghost values Pg(:,3:6) of the connection record r (viscous runs): the call of the loop on the records
  !> with a local partner block in Fill_Ghost_Cell (Lib_Ghost.f90), which copies the four tangential neighbours of the
  !> partner cell (cells or ghost cells of the partner block).
  subroutine Pg_Record ( domain, r )
    use MOSE_Lib_Ghost, only: Fill_BC_Ghost_Connection
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: r
    integer :: Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s

    Bs = domain % bc(r) % bs
    Is = domain % bc(r) % is
    Js = domain % bc(r) % js
    Ks = domain % bc(r) % ks
    Fs = domain % bc(r) % fs
    d11s = domain % bc(r) % d11
    d12s = domain % bc(r) % d12
    d21s = domain % bc(r) % d21
    d22s = domain % bc(r) % d22
    call Fill_BC_Ghost_Connection ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, domain % blk(Bs), domain % bc(r) % Pg )
  end subroutine Pg_Record


  !> Heun stage of one cell: copy of Compute_Cons_NoIRS (contained in Newstate_Conservative, Lib_Newstate.f90),
  !> called with dt = h_l of the cell's level. Same statements, same order, and, in an attempt that the redo can redo
  !> (ev_on), the test of a state that is not physical in two parts: (a) after RK_stage, the conservative variables
  !> finite and the total density positive (cons2prim divides by it, and turns an energy that is not finite into a finite
  !> state); (b) after cons2prim, the primitive variables finite and the test of the faces, check_gas_state, on a copy
  !> (pressure not negative, temperature within the tables after the floor of the partial densities). A cell that fails
  !> it keeps the start state of its sub-step (primold), the copy of Check_And_Fix_State is not called for it (that one
  !> repeats the test of the finite values and of the pressure, and floors), and the attempt will be rejected (Cell_Event).
  !> Outside such an attempt the stop of the copy of Check_And_Fix_State, as in the single-rate step.
  subroutine Compute_Cons_NoIRS( prim, primold, residual, dt, volume, irk, n_rk, b, i, j, k, strangcoeff )
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_RK
    use MOSE_Finite_m, only: Not_Finite
    implicit none
    real(R8), intent(inout) :: prim(nprim), residual(nprim)
    real(R8), intent(in)    :: primold(nprim), dt, volume, strangcoeff
    integer,  intent(in)    :: irk, n_rk, b, i, j, k
    real(R8) :: rho, Rgas, temperature, consold(nprim), cons(nprim)
    real(R8) :: rhoi(nsc)
    integer  :: err
    logical  :: bad

    residual = -residual / volume * dt * strangcoeff

    consold(1:np) = prim2cons( primold(1:np) )
    cons   (1:np) = prim2cons( prim   (1:np) )
    if (nprim > np) then
      consold(np+1:nprim) = primold(np+1:nprim)
      cons   (np+1:nprim) = prim   (np+1:nprim)
    endif

    cons = RK_stage( irk, n_rk, cons, consold, residual )

    ! (a) the conservative variables (the values tested by their bits first: no comparison with a NaN, which a build that
    !     traps the invalid operations would stop on)
    if (ev_on) then
      bad = any( Not_Finite(cons(1:np)) )
      if (.not. bad) bad = .not. ( sum(cons(1:nsc)) > 0.0_R8 )
      if (bad) then
        call Cell_Event ( cons, 1, b, i, j, k )
        prim = primold
        return
      endif
    endif

    call co_rotot_Rtot( prim(1:nsc), rho, Rgas )
    temperature = prim(np) / ( rho * Rgas )

    prim(1:np) = cons2prim( cons(1:np), temperature )
    if (nprim > np) prim(np+1:nprim) = cons(np+1:nprim)

    ! (b) the primitive variables
    if (ev_on) then
      bad = any( Not_Finite(prim(1:nprim)) )
      if (.not. bad) then
        rhoi = prim(1:nsc)
        call check_gas_state ( rhoi, prim(np), err )
        bad = err /= 0
      endif
      if (bad) then
        call Cell_Event ( prim, 2, b, i, j, k )
        prim = primold
        return
      endif
    endif

    call Check_And_Fix_State( prim, b, i, j, k )

  end subroutine Compute_Cons_NoIRS


  !> Copy of Check_And_Fix_State (Lib_Newstate.f90), the message of the stop marked (multirate); compared with the
  !> original statement by statement (verify_copies.py).
  subroutine Check_And_Fix_State ( prim, b, i, j, k )
    use MOSE_Lib_RANS
    use MOSE_Lib_Newstate, only: Count_Negative_Turbulence
    use MOSE_Finite_m,     only: Not_Finite
    implicit none
    real(R8), dimension(:), intent(inout) :: prim
    integer, intent(in) :: b, i, j, k
    integer :: s

    if ( any( Not_Finite(prim(1:nprim)) ) .or. prim(np) < 0d0 ) then
      write(*,'(A,4I4)') "Integration failed at b, i, j, k:", b, i, j, k
      write(*,*) prim(1:nprim)
      write(*,'(A)') "NaN, infinity or p<0 detected (multirate)"
      stop 1
    endif

    do s = 1, nsc
      if ( prim(s) < 1d-20 ) prim(s) = 1d-20
    enddo

    ! Realizability for turbulence, after counting a stage that made nu~, k or omega negative (the counter of
    ! Lib_Newstate: the line at the end of the run counts the stages of the multirate too)
    if (model==2) then
      call Count_Negative_Turbulence( prim )
      call RANS_Enforce_Realizability( prim(nt:) )
    endif

    if (nsoot == 2) then
      if ( prim(nc)   < 1d-20 ) prim(nc)   = 1d-20
      if ( prim(nc+1) < 1d+10 ) prim(nc+1) = 1d+10
    elseif (nsoot == 1) then
      if ( prim(nc)   < 1d+10 ) prim(nc)   = 1d+10
    endif

  end subroutine Check_And_Fix_State


  !> Turbulence source of cell (i,j,k) of block blk (RANS runs with more than one level, stage phase 5c): the routine of
  !> the model for one cell, with the arrays and the options that Spalart_Source_Terms and SST_Source_Terms pass to
  !> SA_Source_Blk and SST_Blk, the cell metrics of the sources included (Source_Metrics). The point-implicit factor
  !> serves steady runs only: in the time-accurate runs of the multirate it is off.
  subroutine Source_Cell ( blk, i, j, k )
    use MOSE_Base_Types_m, only: MOSE_tensor_3D_type
    use MOSE_Lib_Stencil_Metrics, only: Source_Metrics
    implicit none
    type(MOSE_block_type), intent(inout), target :: blk
    integer, intent(in) :: i, j, k
    logical :: point_implicit
    type(MOSE_tensor_3D_type), pointer :: Mi(:,:,:), Mc(:,:,:), Mt(:,:,:)   ! cell metrics of the sources (Source_Metrics)

    point_implicit = obj_rans%point_implicit .and. .not. obj_time_scheme%time_accurate
    select case (rans_kind)
    case (RANS_SA)
      call Source_Metrics ( blk, Mi, Mc, Mt=Mt )
      call SA_Source_Cell ( blk % P, blk % r, Mi, Mc, Mt, blk % yn, blk % k_rough, blk % vol, blk % vel_gradient, &
                            blk % rc_term1, blk % rc_term2, blk % dtlocal, blk % dim, &
                            obj_rans%SpalartShur, obj_rans%SAR, obj_rans%SAcomp, point_implicit, i, j, k )
    case (RANS_SST)
      call Source_Metrics ( blk, Mi, Mt=Mt )
      call SST_Source_Cell ( blk % P, blk % r, Mi, Mt, blk % vol, blk % yn, blk % vel_gradient, &
                             blk % rc_term1, blk % rc_term2, blk % dtlocal, blk % dim, &
                             obj_rans%SpalartShur, obj_rans%k_energy_coupling, point_implicit, &
                             obj_rans%sst_compressible, obj_rans%sst_dd, i, j, k )
    end select

  end subroutine Source_Cell

  !> Spalart-Allmaras source of cell (i,j,k) of a block (RANS runs with more than one level, stage phase 5c): the body
  !> of the cell loop of SA_Source_Blk (Lib_Spalart.f90), copied unchanged, with the constants of that module (private
  !> there) and the indices of the cell as arguments. The test mr-onelevel compares
  !> the body and the constants with the original statement by statement (verify_copies.py); compiled in another
  !> context, the copy may differ from the loop in the last bit, so the one-level stage calls the block routine.
  subroutine SA_Source_Cell ( Prim, Res, M, Mc, Mt, yn, k_rough, vol, gradv, rc1, rc2, dt, n, SpalartShur, SAR, &
                              SAcomp, point_implicit, i, j, k )
    use MOSE_Base_Types_m
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_Fluid
    use MOSE_Lib_RotatingFrame, only: obj_rot
    use MOSE_Lib_Spalart, only: fv1, fv2, Compressibility_Correction
    implicit none
    ! constants of MOSE_Lib_Spalart
    real(R8), parameter :: &
    cb1 = 0.1355d0, &
    cb2 = 0.622d0, &
    cw2 = 0.3d0, &
    cw3 = 2.0d0, &
    cv1 = 7.1d0, &
    cK = 0.41d0, &
    ct1 = 1.d0, &
    ct2 = 2.d0, &
    ct3 = 1.2d0, &
    ct4 = 0.5d0, &
    sigma = 2d0/3d0, &
    c2 = 0.7d0, &
    c3 = 0.9d0, &
    cr1 = 1.d0, &
    cr2 = 12.d0, &
    cr3 = 1.d0, &
    Crot = 2.d0, &
    cw1 = cb1/cK**2 + (1d0+cb2)/sigma, &
    C_cr1 = 0.3d0, &
    C_rough1 = 0.5d0, &   ! SA-rough: chi shift coefficient
    C_rough_d0 = 0.03d0, & ! SA-rough: wall displacement d0 = C_rough_d0 * ks
    tolerance = 1.d-2
    integer, intent(in) :: n(3), i, j, k
    logical, intent(in) :: SpalartShur, SAR, SAcomp, point_implicit
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: yn
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: k_rough
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M    ! velocity
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Mc   ! density
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Mt   ! turbulence variable
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: gradv
    real(R8), dimension(n(1), n(2), n(3)), intent(in) :: rc1, rc2
    real(R8), dimension(n(1), n(2), n(3)), intent(in) :: dt
    ! Local
    real(R8) :: rho, Rgas, nit, mil, Gradvel(3,3), Gradnit(3), omega(3), Om, abs_Eij, Wij(3,3)
    real(R8) :: Eij(3,3), chi, Stilde, Sbar, fw, Gradnit2, r_Fun, g_Fun, fr1, D2, rtilde
    real(R8) :: Production, Diffusion, Destruction, Source, fnu, d, fv2_, Gradrho(3)


      call Co_rotot_Rtot ( Prim(1:nsc,i,j,k), rho, Rgas )
      mil = f_laminarViscosity ( Prim(1:nsc,i,j,k), Prim(np,i,j,k), rho, Rgas )
      nit = Prim (nt,i,j,k) / rho  ! note: its ni-tilde
      chi = Prim (nt,i,j,k) / mil  ! Spalart X variable

      ! Wall distance, displaced below a rough wall (SA-rough)
      d = yn(i,j,k)
      if ( k_rough(i,j,k) > 0d0 ) then
        d    = d + C_rough_d0 * k_rough(i,j,k)
        chi  = chi + C_rough1 * k_rough(i,j,k) / d
        fv2_ = 1d0 - nit / ( mil / rho + nit * fv1(chi) )
      else
        fv2_ = fv2 (chi)
      end if

      if ( SpalartShur ) then
        Gradvel = gradv(i,j,k) % c
      else 
        Gradvel(:,1) = ( Prim(nu:nw,i+1,j,k) - Prim(nu:nw,i-1,j,k) ) * 0.5d0
        Gradvel(:,2) = ( Prim(nu:nw,i,j+1,k) - Prim(nu:nw,i,j-1,k) ) * 0.5d0
        Gradvel(:,3) = ( Prim(nu:nw,i,j,k+1) - Prim(nu:nw,i,j,k-1) ) * 0.5d0
        Gradvel = Matmul ( Gradvel, M(i,j,k) % c )
      end if

      omega = Vorticity_Vector ( Gradvel )    ! vorticity vector (relative frame)
      ! Use absolute vorticity for SA: omega_abs = omega_rel + 2*Omega_frame.
      ! This ensures Om, Stilde, and D2 (SA-RC denominator) are physically correct.
      if (obj_rot%enabled) omega = omega + 2.0_R8 * obj_rot%Omega_vec
      Om    = Sqrt ( sum ( omega**2 ) )     ! vorticity module

      ! Production term computation
      Sbar = nit / ( cK**2 * d**2 ) * fv2_
      if ( Sbar >= (-c2*Om) ) then
        Stilde = Om + Sbar
      else
        Stilde = Om + Om * (c2**2*Om + c3*Sbar) / ((c3 - 2*c2)*Om - Sbar)
      endif
      if ( Stilde == 0d0 ) Stilde = 1d-40 ! Just to avoid division by zero
        
      Eij     = Strain_Tensor ( Gradvel ) ! Strain rate tensor
      abs_Eij = sqrt ( 2*sum (Eij**2) )

      if ( SpalartShur ) then

        Wij = Vorticity_Tensor ( omega )
        ! rc1 is 2 r*/(1 + r*) (Compute_RC_Terms); where strain and vorticity both vanish rc2 is zero, and so is rtilde
        D2 = 0.5d0*(abs_Eij**2 + Om**2)
        rtilde = 0d0
        if ( D2 > 0d0 ) rtilde = rc2(i,j,k) / D2**2
        fr1 = (1.d0 + cr1)*rc1(i,j,k)*(1.d0 - cr3*datan(cr2*rtilde)) - cr1
        Production = fr1 * cb1 * Stilde * nit

      elseif ( SAR ) then
          
        Production = cb1 * ( Stilde + Crot * Min(0d0, abs_Eij - Om) ) * nit

      else

        Production = cb1 * Stilde * nit
        
      end if

      ! finite differences for ni-tilde in computational reference frame
      Gradnit(1) = ( Prim(nt,i+1,j,k) / sum( Prim(1:nsc,i+1,j,k) ) - &
                     Prim(nt,i-1,j,k) / sum( Prim(1:nsc,i-1,j,k) ) ) * 0.5d0
      Gradnit(2) = ( Prim(nt,i,j+1,k) / sum( Prim(1:nsc,i,j+1,k) ) - &
                     Prim(nt,i,j-1,k) / sum( Prim(1:nsc,i,j-1,k) ) ) * 0.5d0
      Gradnit(3) = ( Prim(nt,i,j,k+1) / sum( Prim(1:nsc,i,j,k+1) ) - &
                     Prim(nt,i,j,k-1) / sum( Prim(1:nsc,i,j,k-1) ) ) * 0.5d0

      Gradnit = Matmul ( Gradnit, Mt(i,j,k) % c ) ! ni-tilde gradient

      Gradnit2 = Sum ( Gradnit**2 )              ! ni-tilde gradient squared

      ! Density gradient, same differences
      Gradrho(1) = ( sum( Prim(1:nsc,i+1,j,k) ) - sum( Prim(1:nsc,i-1,j,k) ) ) * 0.5d0
      Gradrho(2) = ( sum( Prim(1:nsc,i,j+1,k) ) - sum( Prim(1:nsc,i,j-1,k) ) ) * 0.5d0
      Gradrho(3) = ( sum( Prim(1:nsc,i,j,k+1) ) - sum( Prim(1:nsc,i,j,k-1) ) ) * 0.5d0
      Gradrho = Matmul ( Gradrho, Mc(i,j,k) % c )   ! behind a wall the density follows the temperature (Ghost_Wall)

      ! Diffusion terms outside the flux div((mu + rho nu~) grad nu~)/sigma: the cb2 term and, with variable density,
      ! -(nu + nu~) grad rho . grad nu~ / sigma, which makes this conservation form equal to the original equation
      ! (S. R. Allmaras, F. T. Johnson, P. R. Spalart, ICCFD7-1902, 2012, eq. 9). The source below is multiplied by
      ! rho, so the second term enters divided by rho.
      Diffusion = cb2 / sigma * Gradnit2 &
                - ( mil / rho + nit ) * Dot_Product ( Gradrho, Gradnit ) / ( sigma * rho )

      ! Destruction term computation
      r_Fun = nit / ( Stilde * cK**2 * d**2 )
      r_Fun = Min ( r_Fun, 10d0 )
      g_Fun = r_Fun + cw2 * ( r_Fun**6 - r_Fun )
      fw = g_Fun * ( ( 1d0 + cw3**6 ) / ( g_Fun**6 + cw3**6 ) )**(1.d0/6.d0)
      Destruction = cw1 * fw * ( nit / d )**2

      if ( SAcomp ) call Compressibility_Correction &
      ( nit, rho, chi, Prim(1:nsc,i,j,k), Prim(np,i,j,k), Rgas, Om, Production )

      ! Point-implicit factor of steady runs. It multiplies the whole increment of the row: the fluxes already in
      ! Res(nt) (Internal_Fluxes, BC_Fluxes) here, the source below. The steady state, zero residual, then does not
      ! depend on the step (point-implicit form of Bussing and Murman: H. C. Yee, J. L. Shinn, NASA TM-89415, 1986,
      ! eq. 3.1; T. H. Pulliam, AIAA Paper 93-3360, eq. 8)
      fnu = 1d0
      if ( point_implicit ) then
        fnu = 1d0 / ( 1d0 + dt(i,j,k) * 2d0 * cw1 * fw * nit / d**2 )
        Res(nt,i,j,k) = Res(nt,i,j,k) * fnu
      end if

      ! Source multiplied by rho since the SA equation is integrated in conservative form
      Source = rho * ( Production + Diffusion - Destruction ) * vol(i,j,k) * fnu
      Res(nt,i,j,k) = Res(nt,i,j,k) - Source

  end subroutine SA_Source_Cell


  !> SST source of cell (i,j,k) of a block (RANS runs with more than one level, stage phase 5c): the body of the cell
  !> loop of SST_Blk (Lib_SST.f90), copied unchanged, with the constants of that module (private there) and the indices
  !> of the cell as arguments. The test mr-onelevel compares the body and the
  !> constants with the original statement by statement (verify_copies.py); compiled in another context, the copy may
  !> differ from the loop in the last bit, so the one-level stage calls the block routine.
  subroutine SST_Source_Cell ( Prim, Res, M, Mt, Volume, WDist, gradv, rc1, rc2, dt, n, SpalartShur, &
                               k_energy_coupling, point_implicit, compressible, dd, i, j, k )
    use MOSE_Base_Types_m
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_Fluid
    use MOSE_Lib_RotatingFrame, only: obj_rot
    use MOSE_Lib_SST, only: compute_F
    implicit none
    ! constants of MOSE_Lib_SST
    real(R8), parameter :: &
  sigma_k1 = 0.85d0, &
  sigma_k2 = 1.0d0, &
  sigma_w1 = 0.5d0, &
  sigma_w2 = 0.856d0, &
  beta_1 = 0.075d0, &
  beta_2 = 0.0828d0, &
  beta_star = 0.09d0, &
  karman = 0.41d0, &
  alpha_1 = 0.31d0, &
  gamma_1 = 5d0/9d0, &
  gamma_2 = 0.44d0, &
  cr1 = 1.d0, &
  cr2 = 2.d0, &
  cr3 = 1.d0
    real(R8), parameter :: dd_xi_w = 2d0, dd_Mt0_w = 0.25d0, dd_xi_z = 0.75d0, dd_Mt0_z = 0.2d0, dd_L_z = 0.66d0
    integer, intent(in) :: n(3), i, j, k
    logical, intent(in) :: SpalartShur, k_energy_coupling, point_implicit, compressible
    integer, intent(in) :: dd   ! dilatation dissipation: 0 none, 1 Wilcox, 2 Zeman (sst-dilatation-dissipation)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Volume
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: WDist
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Mt   ! turbulence variables
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: gradv
    real(R8), dimension(n(1), n(2), n(3)), intent(in) :: rc1, rc2
    real(R8), dimension(n(1), n(2), n(3)), intent(in) :: dt
    ! Local
    real(R8) :: rho, Rgas, mil, kap, ome, dist, Gradvel(3,3), Divel, Sij(3,3), S, Wij(3,3), mi_t, F(2)
    real(R8) :: vort(3), O, D, rtilde, frot, fr1
    real(R8) :: Prod(2), Grad(2,3), dkDotdw, Diff, beta, gamma, Diss(2), Source(2), fk, fw
    real(R8) :: Mach_t, Fdd, beta_star_c, beta_c

      ! Local variables

      ! Local variables
      call Co_rotot_Rtot ( Prim(1:nsc,i,j,k), rho, Rgas )
      mil = f_laminarViscosity ( Prim(1:nsc,i,j,k), Prim(np,i,j,k), rho, Rgas )
      kap = Prim(nt,i,j,k) / rho
      ome = Prim(nt+1,i,j,k) / rho
      dist = WDist(i,j,k)

      ! Velocity gradient computation
      if ( SpalartShur ) then
        Gradvel = gradv(i,j,k) % c
      else
        Gradvel(:,1) = ( Prim(nu:nw,i+1,j,k) - Prim(nu:nw,i-1,j,k) ) * 0.5d0
        Gradvel(:,2) = ( Prim(nu:nw,i,j+1,k) - Prim(nu:nw,i,j-1,k) ) * 0.5d0
        Gradvel(:,3) = ( Prim(nu:nw,i,j,k+1) - Prim(nu:nw,i,j,k-1) ) * 0.5d0
        Gradvel = Matmul ( Gradvel, M(i,j,k) % c )
      end if
      Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)

      ! Gradient of kappa and omega
      Grad(:,1) = ( Prim(nt:nt+1,i+1,j,k) / sum( Prim(1:nsc,i+1,j,k) ) - &
                    Prim(nt:nt+1,i-1,j,k) / sum( Prim(1:nsc,i-1,j,k) ) ) * 0.5d0
      Grad(:,2) = ( Prim(nt:nt+1,i,j+1,k) / sum( Prim(1:nsc,i,j+1,k) ) - &
                    Prim(nt:nt+1,i,j-1,k) / sum( Prim(1:nsc,i,j-1,k) ) ) * 0.5d0
      Grad(:,3) = ( Prim(nt:nt+1,i,j,k+1) / sum( Prim(1:nsc,i,j,k+1) ) - &
                    Prim(nt:nt+1,i,j,k-1) / sum( Prim(1:nsc,i,j,k-1) ) ) * 0.5d0
      Grad = Matmul ( Grad, Mt(i,j,k) % c )

      ! Blending of coefficients
      F = compute_F ( rho=rho, w=ome, k=kap, mil=mil, gradk=Grad(1,:), gradw=Grad(2,:), d=dist )
      gamma = gamma_1*F(1) + gamma_2*(1d0-F(1))
      beta = beta_1*F(1) + beta_2*(1d0-F(1))

      ! Strain and vorticity
      Sij = Strain_Tensor ( Gradvel )
      S = Sqrt ( 2d0 * sum ( Sij**2 ) )

      ! Eddy viscosity
      vort = vorticity_vector ( gradvel )
      ! Use absolute vorticity: omega_abs = omega_rel + 2*Omega_frame.
      ! Needed for O in the SST-RC denominator (rtilde = rc2 / (O*D^3)) and
      ! for Wij passed to the Spalart-Shur fr1 correction.
      if (obj_rot%enabled) vort = vort + 2.0_R8 * obj_rot%Omega_vec
      Wij = Vorticity_Tensor ( vort )
      O = Sqrt ( 2d0 * sum ( Wij**2 ) )
      mi_t = rho*kap*alpha_1 / Max ( alpha_1*ome, S*F(2) )

      ! Production (see the module header of Lib_SST.f90 for the variants and references)
      Prod(1) = mi_t * S**2 ! SST-2003m (TMR), sst-production = incompressible
      ! sst-production = compressible (exact production, OpenFOAM kOmegaSST arrangement), deviatoric part:
      ! mu_t*(S^2 - 2/3 divU^2) = 2 mu_t S_dev:S_dev >= 0
      if ( compressible ) Prod(1) = Prod(1) - 2d0/3d0 * mi_t * Divel**2

      ! Rotation/curvature correction, Smirnov & Menter (2009), TMR SST-RC
      if ( SpalartShur ) then
        ! rc1 is 2 r*/(1 + r*) (Compute_RC_Terms); where the vorticity vanishes rc2 is zero, and so is rtilde
        D = sqrt( max ( S**2, 0.09d0*ome**2) )
        rtilde = 0d0
        if ( O > 0d0 ) rtilde = rc2(i,j,k) / ( O*D**3 )
        frot = (1.d0 + cr1)*rc1(i,j,k)*(1.d0 - cr3*atan(cr2*rtilde)) - cr1
        fr1 = max ( min(frot, 1.25d0), 0d0 )
        Prod(1) = fr1*Prod(1)
      end if

      ! SST-2003 production limiter (10 beta*, Menter et al. 2003); the omega production uses the
      ! limited value, as corrected on the TMR page (the 2003 paper has the unlimited one)
      Prod(1) = Min( Prod(1), 10d0*beta_star*rho*ome*kap )
      Prod(2) = ( gamma*rho / mi_t ) * Prod(1)

      ! Compressible form, dilatation part (-2/3 rho k delta_ij in tau_ij), outside the limiter
      if ( compressible ) then
        Prod(1) = Prod(1) - 2d0/3d0 * rho * kap * Divel
        Prod(2) = Prod(2) - 2d0/3d0 * gamma * rho * ome * Divel
      end if

      ! Cross Diffusion term
      dkDotdw = Dot_Product ( Grad(1,:), Grad(2,:) )
      Diff = 2d0*(1d0-F(1))*rho*sigma_w2/ome*dkDotdw

      ! Dilatation dissipation (sst-dilatation-dissipation; Rumsey, NASA/TM-2009-215705, eqs. 40-43): the destruction
      ! coefficients become beta* (1 + xi F(Mt)) and beta - beta* xi F(Mt), Mt = sqrt(2k)/a the turbulent Mach number.
      !   Wilcox (free shear layers): F = (Mt^2 - Mt0^2) H(Mt - Mt0), xi = 2, Mt0 = 0.25
      !   Zeman (boundary layers):    F = (1 - exp(-((Mt - Mt0)/L)^2)) H(Mt - Mt0), xi = 0.75, Mt0 = 0.2, L = 0.66
      ! beta - beta* xi F is floored at 0 (Wilcox: Mt > 0.65, beyond the Mt of boundary and shear layers): a negative
      ! beta would feed omega. With none, Fdd = 0 and the coefficients are beta* and beta exactly.
      Fdd = 0d0
      if ( dd > 0 ) then
        Mach_t = sqrt ( 2d0 * kap ) / f_ss ( Prim(1:nsc,i,j,k), Prim(np,i,j,k), rho, Rgas )
        if ( dd == 1 ) then
          if ( Mach_t > dd_Mt0_w ) Fdd = dd_xi_w * ( Mach_t**2 - dd_Mt0_w**2 )
        else
          if ( Mach_t > dd_Mt0_z ) Fdd = dd_xi_z * ( 1d0 - exp ( -((Mach_t - dd_Mt0_z) / dd_L_z)**2 ) )
        end if
      end if
      beta_star_c = beta_star * ( 1d0 + Fdd )
      beta_c = max ( beta - beta_star * Fdd, 0d0 )

      ! Dissipation term
      Diss(1) = beta_star_c*rho*kap*ome
      Diss(2) = beta_c*rho*ome**2

      ! Source terms
      Source(1) = Prod(1) - Diss(1)
      Source(2) = Prod(2) - Diss(2) + Diff

      ! Point-implicit factors of steady runs. Each multiplies the whole increment of its row: the fluxes already in
      ! Res (Internal_Fluxes, BC_Fluxes) here, the source below. The steady state, zero residual, then does not depend
      ! on the step (point-implicit form of Bussing and Murman: H. C. Yee, J. L. Shinn, NASA TM-89415, 1986, eq. 3.1;
      ! T. H. Pulliam, AIAA Paper 93-3360, eq. 8)
      ! With k in the total energy (k_energy_coupling) the energy row has the source -Source(1), whose derivative in rho k
      ! is +beta* omega. The point-implicit system of the k and energy rows, D = [[1 + dt beta* omega, 0],
      ! [-dt beta* omega, 1]], gives the energy row its own increment plus dt beta* omega times that of the k row: before
      ! the factor, the energy row takes dt beta* omega fk = 1 - fk times the residual of the k row, and below the source
      ! with the factor. The row of E_m + k then has the fluxes alone at every stage, and the steady state, where the
      ! energy fluxes balance -Source(1) Volume, does not depend on the step. Written 1 - fk, the line shares no product
      ! with the line of fk: written dt beta* omega fk, it made the compilers round fk otherwise (likely no longer fusing
      ! the product with the + 1 in one multiply-add), and the runs without the coupling changed at round-off. 1 - fk
      ! loses relative digits when dt beta* omega is small, where the term is below the round-off of the energy row.
      ! Without the factor (fk = 1) the line does not run
      fk = 1d0 ; fw = 1d0
      if ( point_implicit ) then
        fk = 1d0 / ( 1d0 + dt(i,j,k) * beta_star_c * ome )
        fw = 1d0 / ( 1d0 + dt(i,j,k) * 2d0 * beta_c * ome )
        if ( k_energy_coupling ) Res(np,i,j,k) = Res(np,i,j,k) + ( 1d0 - fk ) * Res(nt,i,j,k)
        Res(nt,  i,j,k) = Res(nt,  i,j,k) * fk
        Res(nt+1,i,j,k) = Res(nt+1,i,j,k) * fw
      end if

      Res(nt,  i,j,k) = Res(nt,  i,j,k) - Source(1) * Volume(i,j,k) * fk
      Res(nt+1,i,j,k) = Res(nt+1,i,j,k) - Source(2) * Volume(i,j,k) * fw
      if ( k_energy_coupling ) Res(np,i,j,k) = Res(np,i,j,k) + Source(1) * Volume(i,j,k) * fk

  end subroutine SST_Source_Cell


  !> Trace line of a stage (first macro step, trace = .true.): level, stage, abscissa t/H and theta of every
  !> coarser level, as exact fractions (format read by the test mr-trace).
  subroutine Trace_Line ( l, kst, t, taus )
    implicit none
    integer, intent(in) :: l, kst
    integer(I8), intent(in) :: t, taus(0:LMAX_CODE)
    character(len=4096) :: line
    character(len=64) :: tok
    integer :: m
    write(line,'(A,I0,A,I0,A,I0,A,I0,A)') 'MR-TRACE l=', l, ' k=', kst, ' t=', t, '/', ticks_H, ' theta='
    do m = 0, l - 1
      write(tok,'(I0,A,I0,A,I0)') m, ':', t - taus(m), '/', ishft(1_I8, Lu + 1 - m)
      line = trim(line)//' '//trim(tok)
    enddo
    write(*,'(A)') trim(line)
  end subroutine Trace_Line


  ! =================================================================================================================
  !  SEVERAL MPI PROCESSES: connection records between processes, needs, exchanges (every block on one process)
  ! =================================================================================================================

  !> The processes this one exchanges connection records with (peer_t) and, for every record received here, its peer
  !> and position (rpeer, rpos). Every process scans all the records in the same order, so the lists of two processes
  !> match without messages. One process: no peers.
  subroutine Build_Peers ( domain )
#ifdef USE_MPI
    use mpi, only: MPI_STATUS_SIZE
#endif
    use MOSE_Mod_MPI, only: mpi_size_, is_local_block, block_owner
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: r, b, bs, p, rk
    integer, allocatable :: nrr(:), nss(:), pid(:)

    allocate( rpeer(domain%nbound), rpos(domain%nbound) )
    rpeer = 0
    rpos = 0
    npeer = 0
    allocate( nrr(0:max(mpi_size_-1,0)), nss(0:max(mpi_size_-1,0)), pid(0:max(mpi_size_-1,0)) )
    nrr = 0; nss = 0; pid = 0
    if (mpi_size_ > 1) then
      do r = 1, domain % nbound
        if (domain%bc(r)%type /= 101 .and. domain%bc(r)%type /= 201) cycle
        b = domain%bc(r)%b
        bs = domain%bc(r)%bs
        if (is_local_block(b) .and. .not. is_local_block(bs)) nrr(block_owner(bs)) = nrr(block_owner(bs)) + 1
        if (is_local_block(bs) .and. .not. is_local_block(b)) nss(block_owner(b)) = nss(block_owner(b)) + 1
      enddo
      do rk = 0, mpi_size_ - 1
        if (nrr(rk) + nss(rk) == 0) cycle
        npeer = npeer + 1
        pid(rk) = npeer
      enddo
    endif
    allocate( peer(npeer) )
    do rk = 0, mpi_size_ - 1
      p = pid(rk)
      if (p == 0) cycle
      peer(p)%rank = rk
      peer(p)%nr = nrr(rk)
      peer(p)%ns = nss(rk)
      allocate( peer(p)%rrec(max(nrr(rk),1)), peer(p)%rmP(max(nrr(rk),1)), peer(p)%rmG(max(nrr(rk),1)) )
      allocate( peer(p)%srec(max(nss(rk),1)), peer(p)%smP(max(nss(rk),1)), peer(p)%smG(max(nss(rk),1)) )
    enddo
    if (npeer > 0) then
      nrr = 0; nss = 0
      do r = 1, domain % nbound
        if (domain%bc(r)%type /= 101 .and. domain%bc(r)%type /= 201) cycle
        b = domain%bc(r)%b
        bs = domain%bc(r)%bs
        if (is_local_block(b) .and. .not. is_local_block(bs)) then
          rk = block_owner(bs)
          nrr(rk) = nrr(rk) + 1
          peer(pid(rk))%rrec(nrr(rk)) = r
          rpeer(r) = pid(rk)
          rpos(r) = nrr(rk)
        endif
        if (is_local_block(bs) .and. .not. is_local_block(b)) then
          rk = block_owner(b)
          nss(rk) = nss(rk) + 1
          peer(pid(rk))%srec(nss(rk)) = r
        endif
      enddo
    endif
    allocate( xreq(max(2*npeer,1)) )
#ifdef USE_MPI
    allocate( xstat(MPI_STATUS_SIZE, max(2*npeer,1)) )
#endif
    call Clear_Masks ()

  end subroutine Build_Peers


  subroutine Clear_Masks ()
    implicit none
    integer :: p
    do p = 1, npeer
      peer(p)%rmP = 0_I8
      peer(p)%rmG = 0_I8
      peer(p)%smP = 0_I8
      peer(p)%smG = 0_I8
    enddo
  end subroutine Clear_Masks


  !> Needs of one kind (1: ghost cells, rmP -> smP; 2: tangential values Pg, rmG -> smG): every process sends to each
  !> peer the masks of the records it receives from it, and receives the masks of the records it sends. Collective among
  !> the peers. The wait for the messages is part of the exchange wait of MOSE_Mod_Timers (the masks are received in
  !> place: no copy from a buffer); posting the messages is the start of the exchange, timed apart as part of the work of
  !> the process (timer_xstart_begin, timer_xstart_end).
  subroutine Exchange_Masks ( kind )
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: check_mpi_error
#endif
    implicit none
    integer, intent(in) :: kind
#ifdef USE_MPI
    integer :: p, nreq, ierr
    if (npeer == 0) return
    call timer_xstart_begin()
    nreq = 0
    do p = 1, npeer
      if (peer(p)%ns > 0) then
        nreq = nreq + 1
        if (kind == 1) then
          call MPI_Irecv( peer(p)%smP, peer(p)%ns, MPI_INTEGER8, peer(p)%rank, TAG_MR_MASK, MPI_COMM_WORLD, &
                          xreq(nreq), ierr )
        else
          call MPI_Irecv( peer(p)%smG, peer(p)%ns, MPI_INTEGER8, peer(p)%rank, TAG_MR_MASK, MPI_COMM_WORLD, &
                          xreq(nreq), ierr )
        endif
        call check_mpi_error(ierr)
      endif
    enddo
    do p = 1, npeer
      if (peer(p)%nr > 0) then
        nreq = nreq + 1
        if (kind == 1) then
          call MPI_Isend( peer(p)%rmP, peer(p)%nr, MPI_INTEGER8, peer(p)%rank, TAG_MR_MASK, MPI_COMM_WORLD, &
                          xreq(nreq), ierr )
        else
          call MPI_Isend( peer(p)%rmG, peer(p)%nr, MPI_INTEGER8, peer(p)%rank, TAG_MR_MASK, MPI_COMM_WORLD, &
                          xreq(nreq), ierr )
        endif
        call check_mpi_error(ierr)
      endif
    enddo
    call timer_xstart_end()
    call timer_comm_begin()
    if (nreq > 0) then
      call MPI_Waitall( nreq, xreq, xstat(:,1:nreq), ierr )
      call check_mpi_error(ierr)
    endif
    call timer_comm_end()
#endif
  end subroutine Exchange_Masks


  !> Build_Map: levels (previous iterate w) of the partner cells of the records received from other processes, wr(r).
  !> Collective among the peers. The wait for the messages and the copy from the receive buffer are part of the exchange
  !> wait of MOSE_Mod_Timers; the copy into the send buffer and the messages posted before the wait are the start of the
  !> exchange, timed apart as part of the work of the process (timer_xstart_begin, timer_xstart_end).
  subroutine Exchange_Partner_Levels ( domain, w, wr )
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: check_mpi_error
#endif
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    type(mr_block_t), intent(in) :: w(:)
    integer, intent(inout) :: wr(:)
#ifdef USE_MPI
    integer, allocatable :: sb(:), rb(:)
    integer :: p, q, r, ns, nr, os, orr, nreq, ierr
    if (npeer == 0) return
    call timer_xstart_begin()
    ns = 0; nr = 0
    do p = 1, npeer
      ns = ns + peer(p)%ns
      nr = nr + peer(p)%nr
    enddo
    allocate( sb(max(ns,1)), rb(max(nr,1)) )
    nreq = 0
    orr = 0
    do p = 1, npeer
      if (peer(p)%nr == 0) cycle
      nreq = nreq + 1
      call MPI_Irecv( rb(orr+1), peer(p)%nr, MPI_INTEGER, peer(p)%rank, TAG_MR_LEV, MPI_COMM_WORLD, xreq(nreq), ierr )
      call check_mpi_error(ierr)
      orr = orr + peer(p)%nr
    enddo
    os = 0
    do p = 1, npeer
      if (peer(p)%ns == 0) cycle
      do q = 1, peer(p)%ns
        r = peer(p)%srec(q)
        sb(os+q) = w(domain%bc(r)%bs)%lev(domain%bc(r)%is, domain%bc(r)%js, domain%bc(r)%ks)
      enddo
      nreq = nreq + 1
      call MPI_Isend( sb(os+1), peer(p)%ns, MPI_INTEGER, peer(p)%rank, TAG_MR_LEV, MPI_COMM_WORLD, xreq(nreq), ierr )
      call check_mpi_error(ierr)
      os = os + peer(p)%ns
    enddo
    call timer_xstart_end()
    call timer_comm_begin()
    if (nreq > 0) then
      call MPI_Waitall( nreq, xreq, xstat(:,1:nreq), ierr )
      call check_mpi_error(ierr)
    endif
    orr = 0
    do p = 1, npeer
      do q = 1, peer(p)%nr
        wr(peer(p)%rrec(q)) = rb(orr+q)
      enddo
      orr = orr + peer(p)%nr
    enddo
    call timer_comm_end()
#endif
  end subroutine Exchange_Partner_Levels


  !> Exchanges of the stages of each level: ghost cells (P) and tangential values (G) from the masks; shares of the
  !> connection fluxes (F) from the owner flags of Pair_Connections and the level of the cells, the same on the two
  !> sides of a connection, so the two processes agree without messages: a record sent here is the partner of an
  !> owner record of the level, a record received here is a record of the level that does not own its face. Positions
  !> in srec/rrec per peer, totals, buffers.
  subroutine Build_Exchange_Lists ( domain )
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: l, p, q, r, o, c, ms, mr
    integer(I8), allocatable :: fm(:)
    if (allocated(xl)) deallocate(xl)
    allocate( xl(0:Lu) )
    ms = 0; mr = 0
    do l = 0, Lu
      allocate( xl(l)%sP(npeer), xl(l)%rP(npeer), xl(l)%sG(npeer), xl(l)%rG(npeer), xl(l)%sF(npeer), xl(l)%rF(npeer) )
      do p = 1, npeer
        call Positions ( peer(p)%smP, peer(p)%ns, l, xl(l)%sP(p) )
        call Positions ( peer(p)%rmP, peer(p)%nr, l, xl(l)%rP(p) )
        call Positions ( peer(p)%smG, peer(p)%ns, l, xl(l)%sG(p) )
        call Positions ( peer(p)%rmG, peer(p)%nr, l, xl(l)%rG(p) )
        ! shares sent: records whose partner (local) owns the face, at the level of the partner's cell
        allocate( fm(max(peer(p)%ns,1)) )
        fm = 0_I8
        do q = 1, peer(p)%ns
          r = peer(p)%srec(q)
          o = domain%bc(r)%partner
          if (.not. domain%bc(o)%owner) cycle
          c = mb(domain%bc(o)%b)%lev(domain%bc(o)%i, domain%bc(o)%j, domain%bc(o)%k)
          fm(q) = ibset( 0_I8, c )
        enddo
        call Positions ( fm, peer(p)%ns, l, xl(l)%sF(p) )
        deallocate( fm )
        ! shares received: records here that do not own their face, at the level of their cell
        allocate( fm(max(peer(p)%nr,1)) )
        fm = 0_I8
        do q = 1, peer(p)%nr
          r = peer(p)%rrec(q)
          if (domain%bc(r)%owner) cycle
          c = mb(domain%bc(r)%b)%lev(domain%bc(r)%i, domain%bc(r)%j, domain%bc(r)%k)
          fm(q) = ibset( 0_I8, c )
        enddo
        call Positions ( fm, peer(p)%nr, l, xl(l)%rF(p) )
        deallocate( fm )
        xl(l)%nsP = xl(l)%nsP + xl(l)%sP(p)%n
        xl(l)%nrP = xl(l)%nrP + xl(l)%rP(p)%n
        xl(l)%nsG = xl(l)%nsG + xl(l)%sG(p)%n
        xl(l)%nrG = xl(l)%nrG + xl(l)%rG(p)%n
        xl(l)%nsF = xl(l)%nsF + xl(l)%sF(p)%n
        xl(l)%nrF = xl(l)%nrF + xl(l)%rF(p)%n
      enddo
      ms = max(ms, xl(l)%nsP * (gc * nprim + 1), xl(l)%nsG * 4 * nprim, xl(l)%nsF * 2 * nprim)
      mr = max(mr, xl(l)%nrP * (gc * nprim + 1), xl(l)%nrG * 4 * nprim, xl(l)%nrF * 2 * nprim)
    enddo
    if (allocated(xsbuf)) deallocate(xsbuf)
    if (allocated(xrbuf)) deallocate(xrbuf)
    allocate( xsbuf(max(ms,1)), xrbuf(max(mr,1)) )
  end subroutine Build_Exchange_Lists


  !> Positions q = 1 ... n whose mask has bit l.
  subroutine Positions ( mask, n, l, pl )
    implicit none
    integer(I8), intent(in) :: mask(:)
    integer, intent(in) :: n, l
    type(pos_list_t), intent(out) :: pl
    integer :: q, c
    c = 0
    do q = 1, n
      if (btest(mask(q), l)) c = c + 1
    enddo
    pl%n = c
    allocate( pl%q(max(c,1)) )
    c = 0
    do q = 1, n
      if (.not. btest(mask(q), l)) cycle
      c = c + 1
      pl%q(c) = q
    enddo
  end subroutine Positions


  !> Stage of level l, several processes (master thread): sends to each peer its partner cells (depth 1 ... gc) of the
  !> records its stages of level l read, at the stage time (refreshed in step 1 of the stage), with the reference
  !> velocity Ur of the cell at depth 1, and receives into the ghost cells of the records whose partner block is on the
  !> peer, as Ghost_Connection and Ghost_Interrank do (also Pg(:,1:gc); Ur of the ghost cell at depth 1, the one the
  !> connection flux reads, as the Pg message of the single-rate exchange). The wait for the messages and the copy from
  !> the receive buffer are timed: the exchange wait of MOSE_Mod_Timers and xch_wall(l); the copy into the send buffer
  !> and the messages posted before the wait are part of the work of the stage, timed apart as the start of the exchange
  !> (timer_xstart_begin, timer_xstart_end: the line Exchange start of the Timing block).
  subroutine Exchange_P ( domain, l )
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: check_mpi_error
#endif
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l
#ifdef USE_MPI
    integer :: p, q, r, g, n, off, nreq, ierr, ent, b, bs, pos
    real(R8) :: dt
    call timer_xstart_begin()
    ent = gc * nprim + 1
    nreq = 0
    off = 0
    do p = 1, npeer
      n = xl(l)%rP(p)%n
      if (n == 0) cycle
      nreq = nreq + 1
      call MPI_Irecv( xrbuf(off+1), n*ent, MPI_DOUBLE_PRECISION, peer(p)%rank, TAG_MR_P, MPI_COMM_WORLD, &
                      xreq(nreq), ierr )
      call check_mpi_error(ierr)
      off = off + n*ent
    enddo
    off = 0
    do p = 1, npeer
      n = xl(l)%sP(p)%n
      if (n == 0) cycle
      do q = 1, n
        r = peer(p)%srec( xl(l)%sP(p)%q(q) )
        bs = domain%bc(r)%bs
        do g = 1, gc
          pos = off + (q-1)*ent + (g-1)*nprim
          xsbuf(pos+1:pos+nprim) = domain%blk(bs)%P(:, domain%bc(r)%is + guide(domain%bc(r)%fs,1)*(g-1), &
                                                       domain%bc(r)%js + guide(domain%bc(r)%fs,2)*(g-1), &
                                                       domain%bc(r)%ks + guide(domain%bc(r)%fs,3)*(g-1))
        enddo
        xsbuf(off + q*ent) = domain%blk(bs)%Ur(domain%bc(r)%is, domain%bc(r)%js, domain%bc(r)%ks)
      enddo
      nreq = nreq + 1
      call MPI_Isend( xsbuf(off+1), n*ent, MPI_DOUBLE_PRECISION, peer(p)%rank, TAG_MR_P, MPI_COMM_WORLD, &
                      xreq(nreq), ierr )
      call check_mpi_error(ierr)
      off = off + n*ent
    enddo
    call timer_xstart_end()
    call timer_comm_begin()
    if (nreq > 0) then
      call MPI_Waitall( nreq, xreq, xstat(:,1:nreq), ierr )
      call check_mpi_error(ierr)
    endif
    off = 0
    do p = 1, npeer
      n = xl(l)%rP(p)%n
      do q = 1, n
        r = peer(p)%rrec( xl(l)%rP(p)%q(q) )
        b = domain%bc(r)%b
        do g = 1, gc
          pos = off + (q-1)*ent + (g-1)*nprim
          domain%blk(b)%P(:, domain%bc(r)%i - guide(domain%bc(r)%f,1)*g, &
                             domain%bc(r)%j - guide(domain%bc(r)%f,2)*g, &
                             domain%bc(r)%k - guide(domain%bc(r)%f,3)*g) = xrbuf(pos+1:pos+nprim)
          domain%bc(r)%Pg(:,g) = xrbuf(pos+1:pos+nprim)
        enddo
        domain%blk(b)%Ur(domain%bc(r)%i - guide(domain%bc(r)%f,1), domain%bc(r)%j - guide(domain%bc(r)%f,2), &
                         domain%bc(r)%k - guide(domain%bc(r)%f,3)) = xrbuf(off + q*ent)
      enddo
      off = off + n*ent
    enddo
    call timer_comm_end(dt)
    xch_wall(l) = xch_wall(l) + dt
#endif
  end subroutine Exchange_P


  !> Stage of level l, several processes (master thread): sends to each peer the shares Fc(:,1:2) of its records whose
  !> owner record is here (computed by Connection_Record at this stage), and receives those of the records here whose
  !> owner record is on the peer, as exchange_connection_flux does. The wait for the messages and the copy from the
  !> receive buffer are timed, as in Exchange_P.
  subroutine Exchange_F ( domain, l )
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: check_mpi_error
#endif
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l
#ifdef USE_MPI
    integer :: p, q, r, n, off, nreq, ierr, ent, pos
    real(R8) :: dt
    call timer_xstart_begin()
    ent = 2 * nprim
    nreq = 0
    off = 0
    do p = 1, npeer
      n = xl(l)%rF(p)%n
      if (n == 0) cycle
      nreq = nreq + 1
      call MPI_Irecv( xrbuf(off+1), n*ent, MPI_DOUBLE_PRECISION, peer(p)%rank, TAG_MR_F, MPI_COMM_WORLD, &
                      xreq(nreq), ierr )
      call check_mpi_error(ierr)
      off = off + n*ent
    enddo
    off = 0
    do p = 1, npeer
      n = xl(l)%sF(p)%n
      if (n == 0) cycle
      do q = 1, n
        r = peer(p)%srec( xl(l)%sF(p)%q(q) )
        pos = off + (q-1)*ent
        xsbuf(pos+1:pos+nprim) = domain%bc(r)%Fc(:,1)
        xsbuf(pos+nprim+1:pos+2*nprim) = domain%bc(r)%Fc(:,2)
      enddo
      nreq = nreq + 1
      call MPI_Isend( xsbuf(off+1), n*ent, MPI_DOUBLE_PRECISION, peer(p)%rank, TAG_MR_F, MPI_COMM_WORLD, &
                      xreq(nreq), ierr )
      call check_mpi_error(ierr)
      off = off + n*ent
    enddo
    call timer_xstart_end()
    call timer_comm_begin()
    if (nreq > 0) then
      call MPI_Waitall( nreq, xreq, xstat(:,1:nreq), ierr )
      call check_mpi_error(ierr)
    endif
    off = 0
    do p = 1, npeer
      n = xl(l)%rF(p)%n
      do q = 1, n
        r = peer(p)%rrec( xl(l)%rF(p)%q(q) )
        pos = off + (q-1)*ent
        domain%bc(r)%Fc(:,1) = xrbuf(pos+1:pos+nprim)
        domain%bc(r)%Fc(:,2) = xrbuf(pos+nprim+1:pos+2*nprim)
      enddo
      off = off + n*ent
    enddo
    call timer_comm_end(dt)
    xch_wall(l) = xch_wall(l) + dt
#endif
  end subroutine Exchange_F


  !> Stage of level l, several processes, viscous runs (master thread): computes the tangential values Pg(:,3:6) of the
  !> records whose partner block is here and whose stages of level l read them on the peer (Pg_Record, from the cells
  !> and ghost cells refreshed in steps 1, 2 and 2m), sends them, and receives those of the records whose partner block
  !> is on the peer, as exchange_ghost_Pg does. The wait for the messages and the copy from the receive buffer are timed,
  !> as in Exchange_P; the values Pg computed for the sends (Pg_Record) are part of the work of the stage and of the
  !> start of the exchange.
  subroutine Exchange_G ( domain, l )
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: check_mpi_error
#endif
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in) :: l
#ifdef USE_MPI
    integer :: p, q, r, c, n, off, nreq, ierr, ent, pos
    real(R8) :: dt
    call timer_xstart_begin()
    ent = 4 * nprim
    nreq = 0
    off = 0
    do p = 1, npeer
      n = xl(l)%rG(p)%n
      if (n == 0) cycle
      nreq = nreq + 1
      call MPI_Irecv( xrbuf(off+1), n*ent, MPI_DOUBLE_PRECISION, peer(p)%rank, TAG_MR_G, MPI_COMM_WORLD, &
                      xreq(nreq), ierr )
      call check_mpi_error(ierr)
      off = off + n*ent
    enddo
    off = 0
    do p = 1, npeer
      n = xl(l)%sG(p)%n
      if (n == 0) cycle
      do q = 1, n
        r = peer(p)%srec( xl(l)%sG(p)%q(q) )
        call Pg_Record ( domain, r )
        do c = 3, 6
          pos = off + (q-1)*ent + (c-3)*nprim
          xsbuf(pos+1:pos+nprim) = domain%bc(r)%Pg(:,c)
        enddo
      enddo
      nreq = nreq + 1
      call MPI_Isend( xsbuf(off+1), n*ent, MPI_DOUBLE_PRECISION, peer(p)%rank, TAG_MR_G, MPI_COMM_WORLD, &
                      xreq(nreq), ierr )
      call check_mpi_error(ierr)
      off = off + n*ent
    enddo
    call timer_xstart_end()
    call timer_comm_begin()
    if (nreq > 0) then
      call MPI_Waitall( nreq, xreq, xstat(:,1:nreq), ierr )
      call check_mpi_error(ierr)
    endif
    off = 0
    do p = 1, npeer
      n = xl(l)%rG(p)%n
      do q = 1, n
        r = peer(p)%rrec( xl(l)%rG(p)%q(q) )
        do c = 3, 6
          pos = off + (q-1)*ent + (c-3)*nprim
          domain%bc(r)%Pg(:,c) = xrbuf(pos+1:pos+nprim)
        enddo
      enddo
      off = off + n*ent
    enddo
    call timer_comm_end(dt)
    xch_wall(l) = xch_wall(l) + dt
#endif
  end subroutine Exchange_G


  ! The reductions of the macro step below are timed, with several processes, as the collective wait of MOSE_Mod_Timers
  ! (timer_sync_begin, timer_sync_end): they are made outside the parallel regions, by the one thread there.

  !> True on every process when x is true on some process. Collective.
  logical function Any_Process ( x )
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: mpi_size_, check_mpi_error
#endif
    implicit none
    logical, intent(in) :: x
#ifdef USE_MPI
    logical :: y
    integer :: ierr
    if (mpi_size_ > 1) then
      call timer_sync_begin()
      call MPI_Allreduce( x, y, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
      call check_mpi_error(ierr)
      call timer_sync_end()
      Any_Process = y
      return
    endif
#endif
    Any_Process = x
  end function Any_Process


  !> Largest n over the processes. Collective.
  integer function Max_Process ( n )
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: mpi_size_, check_mpi_error
#endif
    implicit none
    integer, intent(in) :: n
#ifdef USE_MPI
    integer :: m, ierr
    if (mpi_size_ > 1) then
      call timer_sync_begin()
      call MPI_Allreduce( n, m, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
      call check_mpi_error(ierr)
      call timer_sync_end()
      Max_Process = m
      return
    endif
#endif
    Max_Process = n
  end function Max_Process


  !> a := sum of a over the processes. Collective.
  subroutine Sum_Process ( a )
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: mpi_size_, check_mpi_error
#endif
    implicit none
    integer(I8), intent(inout) :: a(:)
#ifdef USE_MPI
    integer :: ierr
    integer(I8), allocatable :: send(:)
    if (mpi_size_ > 1) then
      ! the send buffer is a copy of a, not MPI_IN_PLACE, which gfortran does not recognise with -fno-underscoring
      send = a
      call timer_sync_begin()
      call MPI_Allreduce( send, a, size(a), MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr )
      call check_mpi_error(ierr)
      call timer_sync_end()
    endif
#endif
  end subroutine Sum_Process


  !> Largest x over the processes (mpi_allreduce_max_r8). Collective.
  real(R8) function Max_R8_Process ( x ) result ( y )
    use MOSE_Mod_MPI, only: mpi_size_, mpi_allreduce_max_r8
    implicit none
    real(R8), intent(in) :: x
    if (mpi_size_ > 1) call timer_sync_begin()
    call mpi_allreduce_max_r8 ( x, y )
    if (mpi_size_ > 1) call timer_sync_end()
  end function Max_R8_Process


  !> Smallest x over the processes (mpi_allreduce_min_r8). Collective.
  real(R8) function Min_R8_Process ( x ) result ( y )
    use MOSE_Mod_MPI, only: mpi_size_, mpi_allreduce_min_r8
    implicit none
    real(R8), intent(in) :: x
    if (mpi_size_ > 1) call timer_sync_begin()
    call mpi_allreduce_min_r8 ( x, y )
    if (mpi_size_ > 1) call timer_sync_end()
  end function Min_R8_Process


  !> a := sum of a over the processes (mpi_allreduce_sum_r8_array). Collective.
  subroutine Sum_R8_Process ( a )
    use MOSE_Mod_MPI, only: mpi_size_, mpi_allreduce_sum_r8_array
    implicit none
    real(R8), intent(inout) :: a(:)
    if (mpi_size_ > 1) call timer_sync_begin()
    call mpi_allreduce_sum_r8_array ( a, size(a) )
    if (mpi_size_ > 1) call timer_sync_end()
  end subroutine Sum_R8_Process


  !> amin, amax, asum := the smallest, the largest and the sum of a over the processes, element by element (a itself on
  !> one process). Collective.
  subroutine Range_R8_Process ( a, amin, amax, asum )
#ifdef USE_MPI
    use mpi
    use MOSE_Mod_MPI, only: mpi_size_, check_mpi_error
#endif
    implicit none
    real(R8), intent(in)  :: a(:)
    real(R8), intent(out) :: amin(:), amax(:), asum(:)
#ifdef USE_MPI
    integer :: ierr
#endif
    amin = a
    amax = a
    asum = a
#ifdef USE_MPI
    if (mpi_size_ > 1) then
      call timer_sync_begin()
      call MPI_Allreduce( a, amin, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
      call check_mpi_error(ierr)
      call MPI_Allreduce( a, amax, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
      call check_mpi_error(ierr)
      call MPI_Allreduce( a, asum, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
      call check_mpi_error(ierr)
      call timer_sync_end()
    endif
#endif
  end subroutine Range_R8_Process


  ! =================================================================================================================
  !  LEVEL MAP AND LISTS (serial, outside the parallel region; lists rebuilt only when the map changes)
  ! =================================================================================================================

  !> Level of a cell from x = sigma*Htry/dtcap_c (with the interface margin, the larger of it and kappa*Htry/min_v
  !> dtcap_v): the least l >= 0 with log2(x) <= l + 1e-12, so that a cell at its limit for a level is not pushed one level
  !> finer by the rounding of the logarithm. The one formula of the levels, used by Build_Map (before the cut at max-level)
  !> and by Depth_Cap.
  pure integer function Level_Of ( x ) result ( l )
    implicit none
    real(R8), intent(in) :: x
    if (x <= 1.0_R8) then
      l = 0
    else
      l = max( ceiling( log(x) / log(2.0_R8) - 1.0e-12_R8 ), 0 )
    endif
  end function Level_Of


  !> Htry of the adaptive step for the depth l: 2**l dt_min/sigma, at most macro-step-target.
  real(R8) function Htry_Of ( l, dtmin ) result ( Ht )
    implicit none
    integer, intent(in)  :: l
    real(R8), intent(in) :: dtmin
    Ht = 2.0_R8**l * dtmin / cfg0%safety
    if (cfg0%macro_step_target > 0.0_R8) Ht = min(Ht, cfg0%macro_step_target)
  end function Htry_Of


  !> Depth for which the map of the adaptive step is built (Htry_Of(l)): the largest l in 0 ... max-level at which the cell
  !> with the largest local step, dtcap_max (every block, every process), has level 0 for Htry_Of(l) in the formula of the
  !> levels alone (Level_Of of sigma*Htry/dtcap_max, without the term of the interface margin). The test is monotone in l
  !> (Htry_Of does not decrease with l) and holds for l = 0 (sigma*Htry_Of(0)/dtcap_max <= dt_min/dtcap_max <= 1), so the
  !> search from max-level down stops at the answer. Without macro-step-target, or with one above dtcap_max/sigma, l =
  !> min(max-level, floor(K + 1e-12)), K = log2(dtcap_max/dt_min) the octaves of the local steps; with a target within
  !> it l = max-level, and the deepest level of the map of Htry_Of(l) is the one the target gives, at most max-level. The
  !> cell with dtcap_max has level 0, so Htry_Of(l) <= dtcap_max/sigma (to 1e-12), and H <= Htry_Of(l): the target bounds
  !> H when it sets Htry_Of(l), otherwise, with kappa <= sigma (the default safety 2 with cfl <= 1), the cell with dt_min
  !> has level l (the margin and the nesting do not raise it), its term in H is Htry_Of(l); so H <= dtcap_max/sigma, the
  !> slowest step over sigma. With kappa > sigma the margin can raise that cell, and a coarse level can stay empty.
  !> Collective (one reduction).
  integer function Depth_Cap ( domain, dtmin ) result ( l )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    real(R8), intent(in) :: dtmin
    real(R8) :: dloc, dmax
    integer  :: b
    dloc = 0.0_R8
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      dloc = max(dloc, maxval(domain%blk(b)%dtlocal(1:domain%blk(b)%dim(1), 1:domain%blk(b)%dim(2), &
                                                     1:domain%blk(b)%dim(3))))
    enddo
    dmax = Max_R8_Process ( dloc )
    l = cfg0%max_level
    do while (l > 0)
      if (Level_Of(cfg0%safety * Htry_Of(l, dtmin) / dmax) == 0) exit
      l = l - 1
    enddo
  end function Depth_Cap


  !> lev_c = min(max-level, Level_Of(sigma*Htry/dtcap_c)) (with the interface margin, Level_Of of the larger of it and
  !> kappa*Htry/min_v dtcap_v), then nesting (|jump| <= 1 across internal
  !> faces) and equal levels on the two sides of every connection, iterated to convergence (Jacobi sweeps: the
  !> least map above the initial one that satisfies both rules). Several processes: at every sweep the levels of the
  !> partner cells on other processes come from them (the previous iterate, as for the local partners), and the sweeps
  !> stop when no process moves, so every process makes the sweeps of a single process. Collective.
  subroutine Build_Map ( domain, Htry )
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: Htry
    ! Local
    type(mr_block_t), allocatable :: w(:)
    integer, allocatable :: wr(:)                        ! partner level of the records received from other processes
    integer :: b, i, j, k, n1, n2, n3, r, lnew, it, lp, ls
    real(R8) :: x, xs, dnb
    logical :: moved, changed, margin

    allocate( w(domain%nb), wr(domain%nbound) )
    wr = -1
    ! interface margin: every cell also within dtcap_v/kappa of each face neighbour v in the block,
    ! kappa = max(cfl (1 + theta), KAPPA_FLOOR): on a face between levels the coarse step is within the TVD bound 1/(1+theta)
    ! of the fine cell (hypothesis (4.15) of Hoang et al.) and within 1/KAPPA_FLOOR of its limit, the margin of condition
    ! (4); with sigma >= kappa it adds nothing on the faces between equal levels or towards coarser cells. Only with
    ! max-level > 0 (one level: no face between levels)
    margin = margin_on
    ring_cells = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      allocate( w(b)%lev(1:n1, 1:n2, 1:n3) )
      do k = 1, n3; do j = 1, n2; do i = 1, n1
        x = cfg0%safety * Htry / domain%blk(b)%dtlocal(i,j,k)
        xs = x
        if (margin) then
          dnb = Dtcap_Neighbours ( domain%blk(b), i, j, k )
          x = max(x, kappa_m * Htry / dnb)
        endif
        w(b)%lev(i,j,k) = min( Level_Of(x), cfg0%max_level )
        if (margin .and. x > xs) then
          ! cell raised by the margin: count it at its level
          ls = min( Level_Of(xs), cfg0%max_level )
          if (w(b)%lev(i,j,k) > ls) ring_cells(w(b)%lev(i,j,k)) = ring_cells(w(b)%lev(i,j,k)) + 1
        endif
      enddo; enddo; enddo
    enddo
    if (margin) then
      call Sum_Process ( ring_cells )
      ring_cells_run = ring_cells_run + ring_cells
      n_rings = n_rings + 1
    endif

    it = 0
    do
      it = it + 1
      moved = .false.
      call Exchange_Partner_Levels ( domain, w, wr )
      ! new map in mb(b)%lev from the previous iterate w(b)%lev
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
        do k = 1, n3; do j = 1, n2; do i = 1, n1
          lnew = w(b)%lev(i,j,k)
          if (i > 1)  lnew = max(lnew, w(b)%lev(i-1,j,k) - 1)
          if (i < n1) lnew = max(lnew, w(b)%lev(i+1,j,k) - 1)
          if (j > 1)  lnew = max(lnew, w(b)%lev(i,j-1,k) - 1)
          if (j < n2) lnew = max(lnew, w(b)%lev(i,j+1,k) - 1)
          if (k > 1)  lnew = max(lnew, w(b)%lev(i,j,k-1) - 1)
          if (k < n3) lnew = max(lnew, w(b)%lev(i,j,k+1) - 1)
          mb(b)%lev(i,j,k) = lnew
        enddo; enddo; enddo
      enddo
      do r = 1, domain % nbound
        if (domain%bc(r)%type /= 101 .and. domain%bc(r)%type /= 201) cycle
        b = domain%bc(r)%b
        if (.not. is_local_block(b)) cycle
        i = domain%bc(r)%i; j = domain%bc(r)%j; k = domain%bc(r)%k
        if (rpeer(r) > 0) then
          lp = wr(r)
        else
          lp = w(domain%bc(r)%bs)%lev(domain%bc(r)%is, domain%bc(r)%js, domain%bc(r)%ks)
        endif
        mb(b)%lev(i,j,k) = max( mb(b)%lev(i,j,k), lp )
      enddo
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
        if (any(mb(b)%lev(1:n1,1:n2,1:n3) /= w(b)%lev)) moved = .true.
        w(b)%lev = mb(b)%lev(1:n1,1:n2,1:n3)
      enddo
      if (.not. Any_Process(moved)) exit
      ! the iterates only rise and, after k sweeps, hold the largest level reachable through at most k faces or partner
      ! records, one less per face: at most (max-level + 1) (m + 1) sweeps, with m the longest chain of partner records of
      ! a cell, that is below this limit for m < 1960 at max-level 50 (a defect of the records otherwise)
      if (it > 100000) error stop 'MOSE-Multirate: level map did not converge within 100000 sweeps (bound (max-level + 1)'// &
                                  '(m + 1) sweeps, m the longest chain of partner records)'
    enddo

    ! compare with the map of the previous build
    changed = .not. have_map
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      if (have_map) then
        if (any(prev_lev(b)%lev /= w(b)%lev)) changed = .true.
      endif
    enddo
    changed = Any_Process(changed)

    if (changed) then
      if (.not. allocated(prev_lev)) allocate( prev_lev(domain%nb) )
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        if (allocated(prev_lev(b)%lev)) deallocate(prev_lev(b)%lev)
        call move_alloc( w(b)%lev, prev_lev(b)%lev )
      enddo
      call Build_Lists ( domain )
      if (mpi_is_root .and. n_macro == 0 .and. cfg0%trace) then
        call Trace_Levels ( domain )
        call Trace_Halo ( domain )
      endif
    endif
    have_map = .true.

  end subroutine Build_Map


  !> Per-level lists of the current map: cells, viscous variables to refresh (viscous runs), face segments, interface
  !> faces, boundary records, ghost records, halos (union of the stencils: offsets +-1, +-2 along each axis for the
  !> convective faces; sources of the ghost records read) and the stash lists. Several processes (collective): the
  !> deepest level is global; the needs of the connection records between processes are agreed level by level (first
  !> the tangential values Pg, whose sources the partner's process then reads, then the ghost cells), and the partner's
  !> process adds the sources it sends to its own halos and ghost records.
  subroutine Build_Lists ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    type(mr_block_t), allocatable :: wk(:)              ! wk(b)%lev: halo mark (last level), wk(b)%grec: stash index
    integer, allocatable :: rflag(:)
    integer :: b, i, j, k, n1, n2, n3, l, m, d, s, r, q, ii, jj, kk, n, nn, a, e, lc, lo, nseg, p1, p2, nfx
    integer :: cnt(0:LMAX_CODE), prevb, prevf, ib
    integer :: ext(3), c(3)

    ! deepest level (of all processes)
    Lu = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      Lu = max(Lu, maxval(mb(b)%lev(1:n1,1:n2,1:n3)))
    enddo
    Lu = Max_Process(Lu)
    lu_run = max(lu_run, Lu)
    if (allocated(lv)) deallocate(lv)
    allocate( lv(0:Lu) )
    if (allocated(cells_per_level)) deallocate(cells_per_level)
    allocate( cells_per_level(0:Lu) )

    ! cells of each level (order b, k, j, i)
    cnt = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        cnt(mb(b)%lev(i,j,k)) = cnt(mb(b)%lev(i,j,k)) + 1
      enddo; enddo; enddo
    enddo
    do l = 0, Lu
      call Alloc_Cells( lv(l)%cells, cnt(l) )
      cells_per_level(l) = cnt(l)
    enddo
    call Sum_Process ( cells_per_level )                  ! all processes (output line of the macro step)
    cnt = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        l = mb(b)%lev(i,j,k)
        cnt(l) = cnt(l) + 1
        lv(l)%cells%b(cnt(l)) = b; lv(l)%cells%i(cnt(l)) = i; lv(l)%cells%j(cnt(l)) = j; lv(l)%cells%k(cnt(l)) = k
      enddo; enddo; enddo
    enddo

    ! viscous runs: the cells and ghost cells whose viscous variables the stages of each level read and refresh
    do l = 0, Lu
      if (model > 0) then
        call Visc_Of_Level ( domain, l )
      else
        call Alloc_Cells ( lv(l)%visc, 0 )
      endif
    enddo

    ! interface faces (internal faces between levels); nesting check
    nif = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n1 = domain%blk(b)%dim(1); n2 = domain%blk(b)%dim(2); n3 = domain%blk(b)%dim(3)
      mb(b)%ifc = 0
      do k = 1, n3; do j = 1, n2; do i = 1, n1
        do d = 1, 3
          ii = i; jj = j; kk = k
          if (d == 1) ii = i + 1
          if (d == 2) jj = j + 1
          if (d == 3) kk = k + 1
          if (ii > n1 .or. jj > n2 .or. kk > n3) cycle
          if (mb(b)%lev(i,j,k) /= mb(b)%lev(ii,jj,kk)) then
            if (abs(mb(b)%lev(i,j,k) - mb(b)%lev(ii,jj,kk)) > 1) error stop 'MOSE-Multirate: nesting violated'
            nif = nif + 1
            mb(b)%ifc(d,i,j,k) = nif
          endif
        enddo
      enddo; enddo; enddo
    enddo
    if (allocated(acc)) deallocate(acc, gsave)
    nfx = FX_CONV
    if (model > 0) nfx = FX_DIFF                           ! viscous runs: accumulators of the diffusive flux too
    allocate( acc(nprim, max(nif,1), nfx), gsave(nprim, max(nif,1), nfx) )
    acc = 0.0_R8
    gsave = 0.0_R8

    ! face segments of each level and direction: maximal runs of faces with at least one cell of the level
    do l = 0, Lu
      do d = 1, 3
        do q = 1, 2                                        ! pass 1 counts, pass 2 fills
          nseg = 0
          do b = 1, domain % nb
            if (.not. is_local_block(b)) cycle
            ext = domain%blk(b)%dim
            if (ext(d) < 2) cycle
            do p2 = 1, ext(Other(d,2)); do p1 = 1, ext(Other(d,1))
              a = 0
              do i = 1, ext(d) - 1
                c = 0
                c(d) = i; c(Other(d,1)) = p1; c(Other(d,2)) = p2
                lc = mb(b)%lev(c(1),c(2),c(3))
                c(d) = i + 1
                lo = mb(b)%lev(c(1),c(2),c(3))
                if (lc == l .or. lo == l) then
                  if (a == 0) a = i
                  e = i
                else
                  if (a /= 0) call Put_Seg(lv(l)%seg(d), q, nseg, b, p1, p2, a, e)
                  a = 0
                endif
              enddo
              if (a /= 0) call Put_Seg(lv(l)%seg(d), q, nseg, b, p1, p2, a, e)
            enddo; enddo
          enddo
          if (q == 1) call Alloc_Seg( lv(l)%seg(d), nseg )
        enddo
      enddo
    enddo

    ! boundary records of the cells of each level, in the order of BC_Fluxes: the list of build_local_bc_index (the
    ! records of the local blocks, block by block and face 1 to 6, in the order of bc.txt within a face), grouped by
    ! (block, face), so that the faces of a corner cell are summed in the order 1..6 whatever the order of bc.txt
    do l = 0, Lu
      n = 0
      do ib = 1, domain % n_local_bc
        r = domain % local_bf_idx(ib)
        b = domain%bc(r)%b
        if (mb(b)%lev(domain%bc(r)%i, domain%bc(r)%j, domain%bc(r)%k) == l) n = n + 1
      enddo
      lv(l)%nbc = n
      allocate( lv(l)%bcr(max(n,1)), lv(l)%grp(n+1) )
      n = 0; nn = 0; prevb = -1; prevf = -1
      do ib = 1, domain % n_local_bc
        r = domain % local_bf_idx(ib)
        b = domain%bc(r)%b
        if (mb(b)%lev(domain%bc(r)%i, domain%bc(r)%j, domain%bc(r)%k) /= l) cycle
        n = n + 1
        lv(l)%bcr(n) = r
        if (b /= prevb .or. domain%bc(r)%f /= prevf) then
          nn = nn + 1
          lv(l)%grp(nn) = n
          prevb = b; prevf = domain%bc(r)%f
        endif
      enddo
      lv(l)%ngrp = nn
      lv(l)%grp(nn+1) = n + 1
    enddo

    ! viscous runs: connection records among the boundary records of each level, read by BC_Connection_Visc with their
    ! tangential ghost values Pg(:,3:6), which the stages of the level refresh (Pg_Record)
    do l = 0, Lu
      n = 0
      if (model > 0) then
        do q = 1, lv(l)%nbc
          r = lv(l)%bcr(q)
          if (domain%bc(r)%type == 101 .or. domain%bc(r)%type == 201) n = n + 1
        enddo
      endif
      lv(l)%npg = n
      allocate( lv(l)%pgr(max(n,1)) )
      n = 0
      if (model > 0) then
        do q = 1, lv(l)%nbc
          r = lv(l)%bcr(q)
          if (domain%bc(r)%type /= 101 .and. domain%bc(r)%type /= 201) cycle
          n = n + 1
          lv(l)%pgr(n) = r
        enddo
      endif
    enddo

    ! connection records among the boundary records of each level that own their face (Pair_Connections): at the stages
    ! of the level they compute the one flux of the face for both cells (Connection_Record)
    do l = 0, Lu
      n = 0
      do q = 1, lv(l)%nbc
        r = lv(l)%bcr(q)
        if (domain%bc(r)%type /= 101 .and. domain%bc(r)%type /= 201) cycle
        if (domain%bc(r)%owner) n = n + 1
      enddo
      lv(l)%ncn = n
      allocate( lv(l)%cnr(max(n,1)) )
      n = 0
      do q = 1, lv(l)%nbc
        r = lv(l)%bcr(q)
        if (domain%bc(r)%type /= 101 .and. domain%bc(r)%type /= 201) cycle
        if (.not. domain%bc(r)%owner) cycle
        n = n + 1
        lv(l)%cnr(n) = r
      enddo
    enddo

    ! halos and ghost records of each level
    allocate( wk(domain%nb), rflag(domain%nbound) )
    rflag = -1
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      allocate( wk(b)%lev (lbound(mb(b)%lev,1):ubound(mb(b)%lev,1), lbound(mb(b)%lev,2):ubound(mb(b)%lev,2), &
                           lbound(mb(b)%lev,3):ubound(mb(b)%lev,3)) )
      allocate( wk(b)%grec(lbound(mb(b)%lev,1):ubound(mb(b)%lev,1), lbound(mb(b)%lev,2):ubound(mb(b)%lev,2), &
                           lbound(mb(b)%lev,3):ubound(mb(b)%lev,3)) )
      wk(b)%lev  = -1          ! last level whose halo took the cell
      wk(b)%grec = 0           ! index of the cell in the stash list of its level
    enddo
    ! several processes: what the stages of each level read from other processes is agreed record by record, as bit
    ! masks of the levels: first the tangential values Pg(:,3:6) of the levels' connection records (the partner's process
    ! computes them, and their sources may be ghost cells of a third process), then the ghost cells
    call Clear_Masks ()
    do l = 0, Lu
      do q = 1, lv(l)%npg
        r = lv(l)%pgr(q)
        if (rpeer(r) > 0) peer(rpeer(r))%rmG(rpos(r)) = ibset( peer(rpeer(r))%rmG(rpos(r)), l )
      enddo
    enddo
    call Exchange_Masks ( 2 )
    do l = 0, Lu
      call Ghost_Records_Of_Level ( domain, l, rflag )
    enddo
    call Exchange_Masks ( 1 )
    do l = 0, Lu
      allocate( lv(l)%halo(0:max(l-1,0)) )
      ! pass 1: halo cells (counts), pass 2: fill
      call Halo_Of_Level ( domain, l, wk )
    enddo
    call Build_Exchange_Lists ( domain )

    ! stash lists: cells of level m read by finer stages; indices into the halo entries
    do m = 0, Lu
      n = 0
      do l = m + 1, Lu
        do q = 1, lv(l)%halo(m)%n
          associate( hb => lv(l)%halo(m)%b(q), hi => lv(l)%halo(m)%i(q), hj => lv(l)%halo(m)%j(q), hk => lv(l)%halo(m)%k(q) )
            if (wk(hb)%grec(hi,hj,hk) == 0) then
              n = n + 1
              wk(hb)%grec(hi,hj,hk) = n
            endif
            lv(l)%halo(m)%u(q) = wk(hb)%grec(hi,hj,hk)
          end associate
        enddo
      enddo
      call Alloc_Cells( lv(m)%used, n )
      allocate( lv(m)%Ws(nprim, max(n,1)), lv(m)%Ub(nprim, max(n,1)), lv(m)%Us(nprim, max(n,1)) )
      do l = m + 1, Lu
        do q = 1, lv(l)%halo(m)%n
          s = lv(l)%halo(m)%u(q)
          lv(m)%used%b(s) = lv(l)%halo(m)%b(q); lv(m)%used%i(s) = lv(l)%halo(m)%i(q)
          lv(m)%used%j(s) = lv(l)%halo(m)%j(q); lv(m)%used%k(s) = lv(l)%halo(m)%k(q)
        enddo
      enddo
    enddo

  end subroutine Build_Lists


  !> Positions read by the stages of level l (the halo of docs/theory/time-integration.md, Recursive Heun scheme),
  !> mk(b)%lev = 1: from every level-l cell the offsets
  !> O_conv = +-1, +-2 along each axis (convective faces) and, in viscous runs, O_diff = the diagonals +-1 +-1 in every
  !> pair of directions (tangential neighbours of the diffusive faces), plus the tangential neighbours of the partner
  !> cells of the level's connection records (Pg(:,3:6), Pg_Record) when the partner block is here. RANS runs: the sources
  !> of a level-l cell read its neighbours +-1 along the axes (central gradients of SA_Source_Blk and SST_Blk), O_src,
  !> a subset of O_conv; the ghost cells they read at a wall (Ghost_Record: Ghost_Wall) belong to the records of the
  !> cell itself, at depth 1, also in O_conv. Several processes: also the tangential neighbours of the partner cells of
  !> the records whose Pg(:,3:6) this process computes for the stages of level l of another process (masks smG).
  subroutine Mark_Level ( domain, l, mk )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in) :: l
    type(mr_block_t), allocatable, intent(out) :: mk(:)
    ! Local
    integer :: b, i, j, k, s, t, q, r, g, p, src(3,4)

    allocate( mk(domain%nb) )
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      allocate( mk(b)%lev(lbound(mb(b)%lev,1):ubound(mb(b)%lev,1), lbound(mb(b)%lev,2):ubound(mb(b)%lev,2), &
                          lbound(mb(b)%lev,3):ubound(mb(b)%lev,3)) )
      mk(b)%lev = 0
    enddo
    do q = 1, lv(l)%cells%n
      b = lv(l)%cells%b(q); i = lv(l)%cells%i(q); j = lv(l)%cells%j(q); k = lv(l)%cells%k(q)
      do s = -gc, gc
        mk(b)%lev(i+s,j,k) = 1
        mk(b)%lev(i,j+s,k) = 1
        mk(b)%lev(i,j,k+s) = 1
      enddo
      if (model > 0) then
        do s = -1, 1, 2
          do t = -1, 1, 2
            mk(b)%lev(i+s,j+t,k) = 1
            mk(b)%lev(i+s,j,k+t) = 1
            mk(b)%lev(i,j+s,k+t) = 1
          enddo
        enddo
      endif
    enddo
    do q = 1, lv(l)%npg
      r = lv(l)%pgr(q)
      if (rpeer(r) > 0) cycle                            ! partner block on another process: that process reads them
      call Pg_Sources ( domain%bc(r), src )
      do g = 1, 4
        mk(domain%bc(r)%bs)%lev(src(1,g), src(2,g), src(3,g)) = 1
      enddo
    enddo
    do p = 1, npeer
      do q = 1, peer(p)%ns
        if (.not. btest(peer(p)%smG(q), l)) cycle
        r = peer(p)%srec(q)
        call Pg_Sources ( domain%bc(r), src )
        do g = 1, 4
          mk(domain%bc(r)%bs)%lev(src(1,g), src(2,g), src(3,g)) = 1
        enddo
      enddo
    enddo

  end subroutine Mark_Level


  !> Ghost records of level l: the records whose ghost positions the stages of level l read (Mark_Level), refreshed by
  !> the stages (Ghost_Record). A connection record whose partner block is on another process is not in the list: its
  !> ghost cells come from that process at the stages of level l (bit l of its mask rmP).
  subroutine Ghost_Records_Of_Level ( domain, l, rflag )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in) :: l
    integer, intent(inout) :: rflag(:)
    ! Local
    type(mr_block_t), allocatable :: mk(:)               ! mk(b)%lev: 1 where read by the stage
    integer :: b, i, j, k, r, ngh

    call Mark_Level ( domain, l, mk )
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do k = lbound(mk(b)%lev,3), ubound(mk(b)%lev,3)
      do j = lbound(mk(b)%lev,2), ubound(mk(b)%lev,2)
      do i = lbound(mk(b)%lev,1), ubound(mk(b)%lev,1)
        if (mk(b)%lev(i,j,k) == 0) cycle
        if (Inside(domain%blk(b)%dim, i, j, k)) cycle
        r = mb(b)%grec(i,j,k)
        if (r == 0) cycle
        if (rpeer(r) > 0) then
          peer(rpeer(r))%rmP(rpos(r)) = ibset( peer(rpeer(r))%rmP(rpos(r)), l )
        else
          rflag(r) = l
        endif
      enddo; enddo; enddo
    enddo
    ngh = count(rflag == l)
    lv(l)%ngh = ngh
    allocate( lv(l)%ghr(max(ngh,1)) )
    ngh = 0
    do r = 1, domain % nbound
      if (rflag(r) /= l) cycle
      ngh = ngh + 1
      lv(l)%ghr(ngh) = r
    enddo

  end subroutine Ghost_Records_Of_Level


  !> Halo of level l: the coarser cells read by the stages of level l, interpolated at their stage times. The interior
  !> positions of Mark_Level, the sources of the level's ghost records (partner cells of a connection up to depth gc;
  !> the boundary cell and the next gc cells inward otherwise) and, with several processes, the partner cells (depth 1
  !> ... gc) of the records whose ghost cells this process sends for the stages of level l of another process (masks
  !> smP); the coarser ones among them.
  subroutine Halo_Of_Level ( domain, l, wk )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in) :: l
    type(mr_block_t), intent(inout) :: wk(:)
    ! Local
    type(mr_block_t), allocatable :: mk(:)               ! mk(b)%lev: 1 where read by the stage
    integer :: b, i, j, k, q, r, g, pass, ii, jj, kk, bs, lsrc, m, p
    integer :: cnt(0:LMAX_CODE)

    if (l == 0) then
      lv(l)%halo(0)%n = 0
      allocate( lv(l)%halo(0)%b(1), lv(l)%halo(0)%i(1), lv(l)%halo(0)%j(1), lv(l)%halo(0)%k(1), lv(l)%halo(0)%u(1) )
      return
    endif

    call Mark_Level ( domain, l, mk )
    do pass = 1, 2
      cnt = 0
      ! interior positions read by the stage
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
          if (mk(b)%lev(i,j,k) == 0) cycle
          call Take ( b, i, j, k )
        enddo; enddo; enddo
      enddo
      ! sources of the ghost records read by the stage
      do q = 1, lv(l)%ngh
        r = lv(l)%ghr(q)
        b = domain%bc(r)%b
        if (domain%bc(r)%type == 101 .or. domain%bc(r)%type == 201) then
          bs = domain%bc(r)%bs
          do g = 1, gc
            ii = domain%bc(r)%is + guide(domain%bc(r)%fs,1)*(g-1)
            jj = domain%bc(r)%js + guide(domain%bc(r)%fs,2)*(g-1)
            kk = domain%bc(r)%ks + guide(domain%bc(r)%fs,3)*(g-1)
            if (Inside(domain%blk(bs)%dim, ii, jj, kk)) call Take ( bs, ii, jj, kk )
          enddo
        else
          do g = 0, gc
            ii = domain%bc(r)%i + guide(domain%bc(r)%f,1)*g
            jj = domain%bc(r)%j + guide(domain%bc(r)%f,2)*g
            kk = domain%bc(r)%k + guide(domain%bc(r)%f,3)*g
            if (Inside(domain%blk(b)%dim, ii, jj, kk)) call Take ( b, ii, jj, kk )
          enddo
        endif
      enddo
      ! several processes: partner cells sent to other processes for their stages of level l
      do p = 1, npeer
        do q = 1, peer(p)%ns
          if (.not. btest(peer(p)%smP(q), l)) cycle
          r = peer(p)%srec(q)
          bs = domain%bc(r)%bs
          do g = 1, gc
            ii = domain%bc(r)%is + guide(domain%bc(r)%fs,1)*(g-1)
            jj = domain%bc(r)%js + guide(domain%bc(r)%fs,2)*(g-1)
            kk = domain%bc(r)%ks + guide(domain%bc(r)%fs,3)*(g-1)
            if (Inside(domain%blk(bs)%dim, ii, jj, kk)) call Take ( bs, ii, jj, kk )
          enddo
        enddo
      enddo
      if (pass == 1) then
        do m = 0, l - 1
          lv(l)%halo(m)%n = cnt(m)
          allocate( lv(l)%halo(m)%b(max(cnt(m),1)), lv(l)%halo(m)%i(max(cnt(m),1)), lv(l)%halo(m)%j(max(cnt(m),1)), &
                    lv(l)%halo(m)%k(max(cnt(m),1)), lv(l)%halo(m)%u(max(cnt(m),1)) )
          lv(l)%halo(m)%u = 0
        enddo
        ! reset the de-duplication marks of this level for pass 2
        do b = 1, domain % nb
          if (.not. is_local_block(b)) cycle
          where (wk(b)%lev == l) wk(b)%lev = -1
        enddo
      endif
    enddo

  contains

    subroutine Take ( bb, ci, cj, ck )
      implicit none
      integer, intent(in) :: bb, ci, cj, ck
      lsrc = mb(bb)%lev(ci,cj,ck)
      if (lsrc < 0 .or. lsrc >= l) return
      if (wk(bb)%lev(ci,cj,ck) == l) return
      wk(bb)%lev(ci,cj,ck) = l
      cnt(lsrc) = cnt(lsrc) + 1
      if (pass == 2) then
        lv(l)%halo(lsrc)%b(cnt(lsrc)) = bb
        lv(l)%halo(lsrc)%i(cnt(lsrc)) = ci
        lv(l)%halo(lsrc)%j(cnt(lsrc)) = cj
        lv(l)%halo(lsrc)%k(cnt(lsrc)) = ck
      endif
    end subroutine Take

  end subroutine Halo_Of_Level


  !> Viscous runs: the cells and ghost cells whose viscous variables blk%Visc the stages of level l read, which they
  !> refresh (Stage, (2v)): from every level-l cell the offsets 0 and +-1 along each axis and the diagonals +-1 +-1 in
  !> every pair of directions, which hold the stencils of the diffusive faces with a level-l cell (the two cells of
  !> the face and their neighbours +-1 in the two other directions, Diffusive_Segment) and the cell and its face
  !> neighbours read by the viscous fluxes of the level's connection and wedge records (BC_Connection_Visc,
  !> BC_Connection_Visc_Stencil, BC_Rotational_Periodic_Visc); the cells and the ghost cells across one face only, as
  !> Viscous_Variables (Mod_Fluxes.f90: no viscous flux reads the ghost cells of the edges and corners). Order b, k, j, i.
  subroutine Visc_Of_Level ( domain, l )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in) :: l
    ! Local
    type(mr_block_t), allocatable :: mk(:)               ! mk(b)%lev: 1 where the stages of level l read blk%Visc
    integer :: b, i, j, k, q, s, t, n(3), pass, cnt

    allocate( mk(domain%nb) )
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n = domain%blk(b)%dim
      allocate( mk(b)%lev(0:n(1)+1, 0:n(2)+1, 0:n(3)+1) )
      mk(b)%lev = 0
    enddo
    do q = 1, lv(l)%cells%n
      b = lv(l)%cells%b(q); i = lv(l)%cells%i(q); j = lv(l)%cells%j(q); k = lv(l)%cells%k(q)
      do s = -1, 1
        mk(b)%lev(i+s,j,k) = 1
        mk(b)%lev(i,j+s,k) = 1
        mk(b)%lev(i,j,k+s) = 1
      enddo
      do s = -1, 1, 2
        do t = -1, 1, 2
          mk(b)%lev(i+s,j+t,k) = 1
          mk(b)%lev(i+s,j,k+t) = 1
          mk(b)%lev(i,j+s,k+t) = 1
        enddo
      enddo
    enddo
    do pass = 1, 2                                         ! pass 1 counts, pass 2 fills
      cnt = 0
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        n = domain%blk(b)%dim
        do k = 0, n(3) + 1; do j = 0, n(2) + 1; do i = 0, n(1) + 1
          if (mk(b)%lev(i,j,k) == 0) cycle
          if ( merge(1,0, i < 1 .or. i > n(1)) + merge(1,0, j < 1 .or. j > n(2)) &
             + merge(1,0, k < 1 .or. k > n(3)) > 1 ) cycle
          cnt = cnt + 1
          if (pass == 2) then
            lv(l)%visc%b(cnt) = b; lv(l)%visc%i(cnt) = i; lv(l)%visc%j(cnt) = j; lv(l)%visc%k(cnt) = k
          endif
        enddo; enddo; enddo
      enddo
      if (pass == 1) call Alloc_Cells ( lv(l)%visc, cnt )
    enddo

  end subroutine Visc_Of_Level


  logical function Inside ( dims, i, j, k )
    implicit none
    integer, intent(in) :: dims(3), i, j, k
    Inside = i >= 1 .and. i <= dims(1) .and. j >= 1 .and. j <= dims(2) .and. k >= 1 .and. k <= dims(3)
  end function Inside


  !> Positions in the partner block read by Fill_BC_Ghost_Connection (Lib_Ghost.f90) for the connection record bc: the
  !> four tangential neighbours of the partner cell (is, js, ks), with the same index rule.
  subroutine Pg_Sources ( bc, src )
    implicit none
    type(MOSE_bc_type), intent(in) :: bc
    integer, intent(out) :: src(3,4)
    associate( Is => bc%is, Js => bc%js, Ks => bc%ks, d11s => bc%d11, d12s => bc%d12, d21s => bc%d21, d22s => bc%d22 )
    select case (bc%fs)
      case (1:2)
        src(:,1) = [Is, Js - d11s, Ks - d12s]
        src(:,2) = [Is, Js + d11s, Ks + d12s]
        src(:,3) = [Is, Js - d21s, Ks - d22s]
        src(:,4) = [Is, Js + d21s, Ks + d22s]
      case (3:4)
        src(:,1) = [Is - d11s, Js, Ks - d12s]
        src(:,2) = [Is + d11s, Js, Ks + d12s]
        src(:,3) = [Is - d21s, Js, Ks - d22s]
        src(:,4) = [Is + d21s, Js, Ks + d22s]
      case default
        src(:,1) = [Is - d11s, Js - d12s, Ks]
        src(:,2) = [Is + d11s, Js + d12s, Ks]
        src(:,3) = [Is - d21s, Js - d22s, Ks]
        src(:,4) = [Is + d21s, Js + d22s, Ks]
    end select
    end associate
  end subroutine Pg_Sources


  integer function Other ( d, n )
    implicit none
    !> the two directions transverse to d: Other(d,1) varies fastest
    integer, intent(in) :: d, n
    integer, parameter :: tab(2,3) = reshape([2,3, 1,3, 1,2], [2,3])
    Other = tab(n, d)
  end function Other


  subroutine Alloc_Cells ( c, n )
    implicit none
    type(cell_list_t), intent(inout) :: c
    integer, intent(in) :: n
    c%n = n
    allocate( c%b(max(n,1)), c%i(max(n,1)), c%j(max(n,1)), c%k(max(n,1)) )
  end subroutine Alloc_Cells


  subroutine Alloc_Seg ( sg, n )
    implicit none
    type(seg_list_t), intent(inout) :: sg
    integer, intent(in) :: n
    sg%n = n
    allocate( sg%b(max(n,1)), sg%p(max(n,1)), sg%q(max(n,1)), sg%a(max(n,1)), sg%e(max(n,1)) )
  end subroutine Alloc_Seg


  subroutine Put_Seg ( sg, pass, nseg, b, p1, p2, a, e )
    implicit none
    type(seg_list_t), intent(inout) :: sg
    integer, intent(in) :: pass, b, p1, p2, a, e
    integer, intent(inout) :: nseg
    nseg = nseg + 1
    if (pass == 2) then
      sg%b(nseg) = b; sg%p(nseg) = p1; sg%q(nseg) = p2; sg%a(nseg) = a; sg%e(nseg) = e
    endif
  end subroutine Put_Seg


  !> Trace (first macro step, trace = .true.): level map of the first macro step on the standard output, one line
  !> "MR-LEVEL b i j k lev" per cell (no file of its own: the outputs of MOSE go through ORION).
  subroutine Trace_Levels ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: b, i, j, k
    write(*,'(A,I0)') 'MR-LEVEL # level map of the first macro step: b i j k lev; deepest level ', Lu
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      do k = 1, domain%blk(b)%dim(3); do j = 1, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        write(*,'(A,5(1X,I0))') 'MR-LEVEL', b, i, j, k, mb(b)%lev(i,j,k)
      enddo; enddo; enddo
    enddo
  end subroutine Trace_Levels


  !> Trace (first macro step, trace = .true.): for every level l its halo cells "MR-HALO H l m b i j k", the ghost
  !> records it refreshes "MR-HALO G l r b i j k f type" (r = index of the record in bc.txt order) and, in viscous runs,
  !> the connection records whose Pg(:,3:6) it refreshes "MR-HALO P l r b i j k f type", on the standard output.
  subroutine Trace_Halo ( domain )
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: l, m, q, r
    write(*,'(A)') 'MR-HALO # halo lists of the first macro step: H l m b i j k | G l r b i j k f type | P l r b i j k f type'
    do l = 0, Lu
      do m = 0, l - 1
        do q = 1, lv(l)%halo(m)%n
          write(*,'(A,7(1X,I0))') 'MR-HALO H', l, m, lv(l)%halo(m)%b(q), lv(l)%halo(m)%i(q), lv(l)%halo(m)%j(q), &
            lv(l)%halo(m)%k(q)
        enddo
      enddo
      do q = 1, lv(l)%ngh
        r = lv(l)%ghr(q)
        write(*,'(A,8(1X,I0))') 'MR-HALO G', l, r, domain%bc(r)%b, domain%bc(r)%i, domain%bc(r)%j, domain%bc(r)%k, &
          domain%bc(r)%f, domain%bc(r)%type
      enddo
      do q = 1, lv(l)%npg
        r = lv(l)%pgr(q)
        write(*,'(A,8(1X,I0))') 'MR-HALO P', l, r, domain%bc(r)%b, domain%bc(r)%i, domain%bc(r)%j, domain%bc(r)%k, &
          domain%bc(r)%f, domain%bc(r)%type
      enddo
    enddo
  end subroutine Trace_Halo


  !> Trace (first macro step, trace = .true., viscous runs), at a stage of level l after its refresh (2v): the viscous
  !> variables blk%Visc that the viscous fluxes of the stage read, against Visc_State of the primitive variables of
  !> the same cell at that point, component by component. Read, as the fluxes of the stage read them: the stencils of
  !> the diffusive faces of the level's segments (the two cells of the face and their neighbours +-1 in the two other
  !> directions, Diffusive_Segment), the cell and its six face neighbours of the level's connection records that own
  !> their face (Connection_Record) and of its wedge records (type 200). A value that differs is stale: computed before
  !> the last change of the primitive variables of its cell. Line "MR-VISC l=.. k=.. t=../.. read=.. stale=.." on the
  !> standard output (read: values read, with repetitions; stale: those that differ). One thread.
  subroutine Trace_Visc ( domain, l, kst, t )
    use MOSE_Lib_Diffusive, only: Visc_State
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in) :: l, kst
    integer(I8), intent(in) :: t
    ! Local
    integer :: d, q, i, s, g, r, b, c(3), x(3), y(3)
    integer(I8) :: nread, nstale

    nread = 0
    nstale = 0
    do d = 1, 3
      do q = 1, lv(l)%seg(d)%n
        b = lv(l)%seg(d)%b(q)
        c(Other(d,1)) = lv(l)%seg(d)%p(q)
        c(Other(d,2)) = lv(l)%seg(d)%q(q)
        do i = lv(l)%seg(d)%a(q), lv(l)%seg(d)%e(q)
          do s = 0, 1                                      ! the two cells of the face (i, i+1 along d)
            x = c
            x(d) = i + s
            call Check ( b, x )
            do g = 1, 2                                    ! their neighbours +-1 in the two other directions
              y = x
              y(Other(d,g)) = x(Other(d,g)) - 1
              call Check ( b, y )
              y(Other(d,g)) = x(Other(d,g)) + 1
              call Check ( b, y )
            enddo
          enddo
        enddo
      enddo
    enddo
    do q = 1, lv(l)%ncn
      call Around ( lv(l)%cnr(q) )
    enddo
    do q = 1, lv(l)%nbc
      r = lv(l)%bcr(q)
      if (domain%bc(r)%type == 200) call Around ( r )
    enddo
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') 'MR-VISC l=', l, ' k=', kst, ' t=', t, '/', ticks_H, ' read=', nread, &
      ' stale=', nstale

  contains

    !> the cell of record rr and its six face neighbours
    subroutine Around ( rr )
      implicit none
      integer, intent(in) :: rr
      integer :: dd, ss, z(3)
      z = [domain%bc(rr)%i, domain%bc(rr)%j, domain%bc(rr)%k]
      call Check ( domain%bc(rr)%b, z )
      do dd = 1, 3
        do ss = -1, 1, 2
          z = [domain%bc(rr)%i, domain%bc(rr)%j, domain%bc(rr)%k]
          z(dd) = z(dd) + ss
          call Check ( domain%bc(rr)%b, z )
        enddo
      enddo
    end subroutine Around

    subroutine Check ( bb, p )
      implicit none
      integer, intent(in) :: bb, p(3)
      real(R8) :: v(nprim)
      call Visc_State ( domain%blk(bb)%P(:,p(1),p(2),p(3)), v )
      nread = nread + 1
      if ( any( v /= domain%blk(bb)%Visc(:,p(1),p(2),p(3)) ) ) nstale = nstale + 1
    end subroutine Check

  end subroutine Trace_Visc

end module MOSE_Mod_Multirate
