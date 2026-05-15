module MOSE_Mod_BC_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: BC_Fluxes

contains

  subroutine BC_Fluxes ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_rans, obj_space_scheme, obj_riemann, obj_soot, obj_rot
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
    use MOSE_Lib_BC_Fluxes_Wall_Melting
    use MOSE_Lib_BC_Fluxes_Wall_Pyrolysis
    use MOSE_Lib_BC_Fluxes_SRM
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer  :: f, lower, upper, i, b
    integer  :: error
    integer  :: Bm, Im, Jm, Km, Fm, Bs, Fs
    real(R8) :: T0, g, pstat, rhot(nsc), visct(1:nrans)
    real(R8) :: r_fc(3)                              ! Face centre [m] (stationary-wall BCs)
    logical  :: RSM, SD_limiter, SD_riemann, soot_enabled
    real(R8) :: Sc, Sct, Prt

    RSM = obj_rans%RSM
    SD_limiter = obj_space_scheme%SD
    SD_riemann = obj_riemann%SD
    Sc  = obj_rans%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt
    soot_enabled = obj_soot%enabled

    ! BC fluxes are computed in order of block and face type.
    upper = 0

    blocks: do b = 1, domain % nb
      faces: do f = 1, 6

        lower = upper + 1                   ! Update lower bound
        upper = upper + domain % n_bf(b,f)  ! Upper bound: add number of cells on face f of block b

        !$omp do schedule ( dynamic ) private(i, Bm, Im, Jm, Km, Fm, error, Bs, Fs, rhot, visct, T0, g, pstat, r_fc)
        do i = lower, upper
          Bm = domain % bc(i) % b
          if (.not. is_local_block(Bm)) cycle
          Im = domain % bc(i) % i
          Jm = domain % bc(i) % j
          Km = domain % bc(i) % k
          Fm = domain % bc(i) % f
          select case ( domain % bc(i) % type )

            case (103) ! multi-Solver coupling
              call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              domain % blk(Bm) % R(:,Im,Jm,Km) = domain % blk(Bm) % R(:,Im,Jm,Km) + domain % bc(i) % ext_flux

            case (101,102,201) ! connection & chimera (101=block connect, 102=chimera)
              call BC_Connection_Eul ( Im, Jm, Km, Fm, domain % blk(Bm), SD_limiter, SD_riemann )
              if (model>0) &
                call BC_Connection_Visc ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Mg(1), domain % bc(i) % Pg, &
                                          Sc, Sct, Prt, soot_enabled )

            case (200) ! periodic
              call BC_Rotational_Periodic_Eul ( Im, Jm, Km, Fm, domain % blk(Bm) )

            case (300) ! symmetry
              call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              if (model>0) &
                call BC_Symmetry_Visc ( Im, Jm, Km, Fm, domain % blk(Bm), RSM )

            case (301) ! wall: prescribed heat flux
              call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              if (model>0.AND.obj_rot%enabled.AND.RF_Is_Stationary_Face(Bm,Fm)) then
                call RF_Face_Center ( domain % blk(Bm) % node, Im, Jm, Km, Fm, r_fc )
                call BC_Wall_Heat ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % qw, &
                                    w_wall = RF_Wall_Velocity(r_fc) )
              elseif (model>0) then
                call BC_Wall_Heat ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % qw )
              endif

            case (302) ! wall: prescribed temperature
              call BC_Symmetry_Eul (Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              if (model>0.AND.obj_rot%enabled.AND.RF_Is_Stationary_Face(Bm,Fm)) then
                call RF_Face_Center ( domain % blk(Bm) % node, Im, Jm, Km, Fm, r_fc )
                call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tw, &
                                           w_wall = RF_Wall_Velocity(r_fc) )
              elseif (model>0) then
                call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tw )
              endif

            case (303) ! wall: temperature + radiative flux (melting)
              call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              if (model>0) &
                call BC_Wall_Melting( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tw, domain % bc(i) % qrad )
                    
            case (304) ! wall: radiative flux only (pyrolysis/ablation)
              call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              if (model>0) &
                call BC_Wall_Pyrolysis( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % qrad, domain % bc(i) % type )

            case (400) ! extrapolation
              call BC_Extrapolation ( Im, Jm, Km, Fm, domain % blk(Bm) )

            case (401) ! inlet: T0 + p0 (stagnation conditions)
              call BC_Inlet_StagnCond ( Bm, Im, Jm, Km, Fm, domain % blk(Bm), &
                                        domain % bc(i) % T0, domain % bc(i) % p0, &
                                        domain % bc(i) % rel_fac, &
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
                                        domain % bc(i) % rel_fac, &
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
                                      domain % bc(i) % rel_fac, &
                                      domain % bc(i) % alpha, domain % bc(i) % beta, &
                                      domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                      domain % bc(i) % ci(1:nrans), &
                                      domain % bc(i) % mdot )

            case (410) ! Mapped BC: Riemann flux with Q2D-interpolated ghost cell
              call BC_Connection_Eul ( Im, Jm, Km, Fm, domain % blk(Bm), SD_limiter, SD_riemann )

            case(420) ! choked nozzle
              call BC_Inflow_Nozzle( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Mach, &
                                     domain % bc(i) % T0, domain % bc(i) % p0, domain % bc(i) % rel_fac, &
                                     domain % bc(i) % alpha, domain % bc(i) % beta, domain % bc(i) % pAmb, &
                                     domain % bc(i) % ci(nrans+1:nrans+nsc), domain % bc(i) % ci(1:nrans), &
                                     domain % bc(i) % mdot, error )
              if (error == 1) then
                call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              endif

            case (501) ! manifold - mean flow connection
              Bs = domain % bc(i) % bs
              Fs = domain % bc(i) % fs
              call BC_Manifold ( domain % blk(Bm), Fm, domain % blk(Bs), Fs, T0, g, pstat )

            case (502) ! SRM grain combustion
              call BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, domain % blk(Bm) )
              call BC_SRM ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % aCoeff, &
                             domain % bc(i) % n, domain % bc(i) % pRef, domain % bc(i) % Taf, &
                             domain % bc(i) % haf, domain % bc(i) % ci )

          end select

        enddo

      enddo faces
    enddo blocks

  end subroutine BC_Fluxes

end module MOSE_Mod_BC_Fluxes