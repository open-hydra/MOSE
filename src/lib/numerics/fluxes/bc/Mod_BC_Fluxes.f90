module MOSE_Mod_BC_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: BC_Fluxes

contains

  subroutine BC_Fluxes ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_rans, obj_soot, obj_rot, obj_sim_param
    use MOSE_Global_m, only: model
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Lib_RotatingFrame
    use MOSE_Lib_BC_Fluxes_Connection
    use MOSE_Lib_BC_Fluxes_Rotational
    use MOSE_Lib_BC_Fluxes_Symmetry
    use MOSE_Lib_BC_Fluxes_Manifold
    use MOSE_Lib_BC_Fluxes_Inflow
    use MOSE_Lib_BC_Fluxes_Outflow
    use MOSE_Lib_BC_Fluxes_FarField
    use MOSE_Lib_BC_Fluxes_Inflow_Nozzle
    use MOSE_Lib_BC_Fluxes_Extrapolation
    use MOSE_Lib_BC_Fluxes_Wall_Heat
    use MOSE_Lib_BC_Fluxes_Wall_Temperature
    use MOSE_Lib_BC_Fluxes_GSI_Melting
    use MOSE_Lib_BC_Fluxes_GSI_Pyrolysis
    use MOSE_Lib_BC_Fluxes_GSI_Reactions
    use MOSE_Lib_BC_Fluxes_GSI_SRM
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer  :: f, lower, upper, i, ii, b
    integer  :: error
    integer  :: Bm, Im, Jm, Km, Fm, Bs, Fs
    real(R8) :: T0, g, pstat, rhot(nsc), visct(1:nrans)
    real(R8) :: r_fc(3)                              ! Face centre [m] (stationary-wall BCs)
    logical  :: RSM, soot_enabled
    logical  :: ablating                             ! GSI walls: is the surface consuming material?
    real(R8) :: Tw_surface                           ! GSI walls: wall temperature for the inert fallback
    real(R8) :: Sc, Sct, Prt, Prl
    real(R8) :: dR(nprim)                            ! Share of a face in the residual of its cell

    RSM = obj_rans%RSM
    Sc  = obj_sim_param%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt
    Prl = obj_sim_param%Prl
    soot_enabled = obj_soot%enabled

    ! Connection faces (101, 201): one flux per face, for the two cells
    call Connection_Fluxes ( domain, Sc, Sct, Prt, Prl, soot_enabled )

    ! BC fluxes are computed block by block and, within a block, face by face, over the
    ! records of the blocks this rank owns only (lists built once by build_local_bc_index):
    ! no rank scans the records of the whole domain. The records of one face of one block
    ! write distinct target cells (the periodic BC also writes the partner cell on the
    ! opposite face, one per record), so one worksharing loop over them is race-free; its
    ! implicit barrier keeps the face contributions of a corner cell summed in the order 1..6.
    blocks: do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle blocks
      faces: do f = 1, 6

        lower = domain % local_bf_first(b,f)
        upper = lower + domain % local_bf_count(b,f) - 1

        !$omp do schedule (dynamic, 64) private(ii, i, Bm, Im, Jm, Km, Fm, error, Bs, Fs, rhot, visct, &
        !$omp   T0, g, pstat, r_fc, ablating, Tw_surface, dR)
        do ii = lower, upper
          i  = domain % local_bf_idx(ii)
          Bm = domain % bc(i) % b
          Im = domain % bc(i) % i
          Jm = domain % bc(i) % j
          Km = domain % bc(i) % k
          Fm = domain % bc(i) % f
          select case ( domain % bc(i) % type )

            case (103,104) ! multi-Solver coupling (104 = with chimera)
              call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + domain % bc(i) % ext_flux

            case (101,201) ! block connection, periodic connection: the share of the cell in the flux of the face
              domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + domain % bc(i) % Fc(:,1)
              if (model>0) &
                domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + domain % bc(i) % Fc(:,2)

            case (102) ! chimera
              call BC_Connection_Eul ( Im, Jm, Km, Fm, domain % blk(Bm), dR )
              domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + dR
              if (model>0) then
                call BC_Connection_Visc ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Mg(1), domain % bc(i) % Pg, &
                                          Sc, Sct, Prt, Prl, soot_enabled, dR )
                domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + dR
              end if

            case (200) ! periodic
              call BC_Rotational_Periodic_Eul ( Im, Jm, Km, Fm, domain % blk(Bm) )
              if (model>0) &
                call BC_Rotational_Periodic_Visc ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % rM, &
                                                   Sc, Sct, Prt, Prl, soot_enabled )

            case (300) ! symmetry
              call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              if (model>0) &
                call BC_Symmetry_Visc ( Im, Jm, Km, Fm, domain % blk(Bm), RSM )

            case (301) ! wall: prescribed heat flux
              call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              if (model>0.AND.obj_rot%enabled.AND.RF_Is_Stationary_Face(Bm,Fm)) then
                call RF_Face_Center ( domain % blk(Bm) % node, Im, Jm, Km, Fm, r_fc )
                call BC_Wall_Heat ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % qw, &
                                    w_wall = RF_Wall_Velocity(r_fc), k_rough = domain % bc(i) % k_rough, &
                                    rec = domain % bc(i), w_grad = RF_Wall_Velocity_Gradient() )
              elseif (model>0) then
                call BC_Wall_Heat ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % qw, &
                                    k_rough = domain % bc(i) % k_rough, rec = domain % bc(i) )
              endif

            case (302) ! wall: prescribed temperature
              call BC_Symmetry_Eul (Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              if (model>0.AND.obj_rot%enabled.AND.RF_Is_Stationary_Face(Bm,Fm)) then
                call RF_Face_Center ( domain % blk(Bm) % node, Im, Jm, Km, Fm, r_fc )
                call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tw, &
                                           w_wall = RF_Wall_Velocity(r_fc), k_rough = domain % bc(i) % k_rough, &
                                           rec = domain % bc(i), w_grad = RF_Wall_Velocity_Gradient() )
              elseif (model>0) then
                call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tw, &
                                           k_rough = domain % bc(i) % k_rough, rec = domain % bc(i) )
              endif

            case (400) ! extrapolation
              call BC_Extrapolation ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

            case (401) ! inlet: T0 + p0 (stagnation conditions)
              call BC_Inlet_StagnCond ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), &
                                        domain % bc(i) % T0, domain % bc(i) % p0, &
                                        domain % bc(i) % alpha, domain % bc(i) % beta, &
                                        domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                        domain % bc(i) % ci(1:nrans), &
                                        domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

            case (402) ! inlet: T0 + time-varying p0 (stagnation conditions, subsonic)
              if (domain % bc(i) % p0time % exists) &
                domain % bc(i) % p0 = domain % bc(i) % p0time % update(domain % time)
              call BC_Inlet_StagnCond ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), &
                                        domain % bc(i) % T0, domain % bc(i) % p0, &
                                        domain % bc(i) % alpha, domain % bc(i) % beta, &
                                        domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                        domain % bc(i) % ci(1:nrans), &
                                        domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

            case (403) ! inlet: total temperature T0 + prescribed mass flux
              call BC_Inlet_MassFlux_T0 ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), &
                                           domain % bc(i) % T0, domain % bc(i) % mdot, &
                                           domain % bc(i) % rel_fac, &
                                           domain % bc(i) % alpha, domain % bc(i) % beta, &
                                           domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                           domain % bc(i) % ci(1:nrans), error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

            case (404) ! inlet: static temperature T + prescribed mass flux
              call BC_Inlet_MassFlux_T ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), &
                                          domain % bc(i) % T0, domain % bc(i) % mdot, &
                                          domain % bc(i) % rel_fac, &
                                          domain % bc(i) % alpha, domain % bc(i) % beta, &
                                          domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                          domain % bc(i) % ci(1:nrans), error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

            case (405) ! inlet: supersonic — Mach + static temperature and static pressure
              call BC_Inlet_Supersonic_Static ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), &
                                                domain % bc(i) % mach, &
                                                domain % bc(i) % T0, domain % bc(i) % pamb, &
                                                domain % bc(i) % rel_fac, &
                                                domain % bc(i) % alpha, domain % bc(i) % beta, &
                                                domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                                domain % bc(i) % ci(1:nrans), &
                                                domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

            case (406) ! outlet
              call BC_Outflow ( Bm, Im, Jm, Km, Fm, domain % blk(Bm),domain % bc(i) % pAmb, domain % bc(i) % rel_fac, &
                                       domain % bc(i) % mdot, error )
              if (error == 1) then
                call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              endif

            case (407) ! ambient: T0 + p0 + back pressure (inflow or outflow)
              call BC_FarField ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), &
                                      domain % bc(i) % T0, domain % bc(i) % p0, &
                                      domain % bc(i) % pamb, &
                                      domain % bc(i) % alpha, domain % bc(i) % beta, &
                                      domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                      domain % bc(i) % ci(1:nrans), &
                                      domain % bc(i) % mdot )

            case (408) ! inlet: static temperature T + prescribed velocity (u) normal to the boundary
              call BC_Inlet_un_T ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), &
                                      domain % bc(i) % T0, domain % bc(i) % un, &
                                      domain % bc(i) % rel_fac, &
                                      domain % bc(i) % alpha, domain % bc(i) % beta, &
                                      domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                      domain % bc(i) % ci(1:nrans), error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )

            case (410) ! Mapped BC: Riemann flux with Q2D-interpolated ghost cell
              call BC_Connection_Eul ( Im, Jm, Km, Fm, domain % blk(Bm), dR )
              domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + dR

            case(420) ! choked nozzle
              call BC_Inflow_Nozzle( Im, Jm, Km, Fm, domain % blk(Bm), &
                                     domain % bc(i) % T0, domain % bc(i) % p0, domain % bc(i) % rel_fac, &
                                     domain % bc(i) % psub, domain % bc(i) % psup, domain % bc(i) % mdot, &
                                     domain % bc(i) % ci(nrans+1:nrans+nsc), domain % bc(i) % ci(1:nrans), error )
              if (error == 1) then
                call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              endif

            case (501) ! manifold - mean flow connection
              Bs = domain % bc(i) % bs
              Fs = domain % bc(i) % fs
              call BC_Manifold ( domain % blk(Bm), Fm, domain % blk(Bs), Fs, T0, g, pstat )

            case (502) ! SRM grain combustion
              ! Like the GSI walls below, a burning grain lays down the complete boundary
              ! flux, wall pressure included, so no reflective flux is applied alongside
              ! it. An inhibited surface (a = 0) burns nothing and falls back to an
              ! impermeable wall at the flame temperature. Unlike the GSI walls this one
              ! also runs inviscid: the mass injection is the point, the stress is not.
              ablating = .false.
              Tw_surface = domain % bc(i) % Taf
              call BC_SRM ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % aCoeff, &
                             domain % bc(i) % n, domain % bc(i) % pRef, domain % bc(i) % Taf, &
                             domain % bc(i) % ci, &
                             burning = ablating, T_surface = Tw_surface )
              if (.not. ablating) then
                call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
                if (model>0) &
                  call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), Tw_surface, rec = domain % bc(i) )
              endif

            ! ─────────────────────────────────────────────────────────────────
            ! Gas-surface interaction walls (503-505).
            !
            ! While the surface ablates, the GSI routine lays down the complete
            ! boundary flux -- wall pressure, the momentum and enthalpy convected in
            ! with the blown mass, and the viscous terms -- so no reflective wall
            ! flux is applied alongside it; adding one would impose u.n = 0 at a face
            ! that is deliberately permeable.
            !
            ! When the flow cannot sustain the process the routine touches nothing
            ! and reports back, and the face is closed as an ordinary impermeable
            ! wall: reflective inviscid flux plus a no-slip wall at the temperature
            ! the surface settled on. An inviscid run (model = 0) never ablates
            ! either, and takes the same path.
            case (503) ! GSI - melting
              ablating = .false.
              Tw_surface = domain % bc(i) % Tw
              if (model>0) &
                call BC_Wall_Melting( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % cp_wall, domain % bc(i) % Tw, &
                                                                        domain % bc(i) % Ti_wall, domain % bc(i) % dh_wall, &
                                                                        domain % bc(i) % qrad, domain % bc(i) % ci, &
                                                                        ablating = ablating, T_surface = Tw_surface )
              if (.not. ablating) then
                call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
                if (model>0) &
                  call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), Tw_surface, rec = domain % bc(i) )
              endif

            case (504) ! GSI - pyrolysis
              ablating = .false.
              Tw_surface = 0d0
              if (model>0) &
                call BC_Wall_Pyrolysis( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % GSI_pyro_model_id, &
                                                                          domain % bc(i) % qrad, domain % bc(i) % ci, &
                                                                          ablating = ablating, T_surface = Tw_surface )
              if (.not. ablating) then
                call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
                if (model>0) &
                  call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), Tw_surface, rec = domain % bc(i) )
              endif

            case(505) ! GSI - surface reactions
              ablating = .false.
              Tw_surface = 0d0
              if (model>0) &
                call BC_Wall_Reactions( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % GSI_surf_reac_id, &
                                        domain % bc(i) % qrad, ablating = ablating, T_surface = Tw_surface )
              if (.not. ablating) then
                call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
                if (model>0) &
                  call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), Tw_surface, rec = domain % bc(i) )
              endif

            ! case(506) ! GSI - surface reactions + pyrolysis
            ! TODO: implement this case in the future

          end select

        enddo

      enddo faces
    enddo blocks

    ! The lines of the gas entering through 400 faces at this stage, before the source terms and the update of the cells
    ! (after the implicit barrier of the last loop above: the loops over the records have no nowait)
    call Report_400_Entering_Lines ()

  end subroutine BC_Fluxes


  !> One flux per connection face (101, 201). Computed by each block for its
  !> own cell, the flux of a face between blocks would be the same on the two
  !> sides only when their inputs are: not with the shock sensor, where each side
  !> takes the switch of its own cell, nor with the preconditioned solver, unless
  !> the ghost carries the reference velocity of the other block. Here the record that owns the
  !> face (Pair_Connections: the cell on the low-index side of the face, as on an
  !> interior face) computes it, inviscid and viscous, for its cell, with the
  !> shock sensor of its cell and the reference velocities of both cells (the
  !> ghost carries that of the other block), and gives it with the opposite sign
  !> to the record of the other cell: in the same process directly, between
  !> processes through exchange_connection_flux. BC_Fluxes then adds to each
  !> cell its share, in the order of the records.
  subroutine Connection_Fluxes ( domain, Sc, Sct, Prt, Prl, soot_enabled )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_space_scheme
    use MOSE_Global_m, only: model, nprim
    use MOSE_Mod_MPI, only: mpi_size_
    use MOSE_Mod_GhostExchange, only: exchange_connection_flux
    use MOSE_Lib_BC_Fluxes_Connection
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: Sc, Sct, Prt, Prl
    logical, intent(in)  :: soot_enabled
    ! Local
    integer  :: ii, i, p, Bm
    real(R8) :: dR(nprim)

    !$omp do schedule (dynamic, 64) private(ii, i, p, Bm, dR)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 101 .and. domain % bc(i) % type /= 201) cycle
      if (.not. domain % bc(i) % owner) cycle
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
    enddo

    ! Faces whose other cell is on another process (MPI calls from the primary thread)
    if (mpi_size_ > 1) then
      !$omp master
      call exchange_connection_flux ( domain )
      !$omp end master
      !$omp barrier
    end if

  end subroutine Connection_Fluxes

end module MOSE_Mod_BC_Fluxes