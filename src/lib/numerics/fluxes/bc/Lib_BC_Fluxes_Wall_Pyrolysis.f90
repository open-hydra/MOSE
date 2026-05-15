module MOSE_Lib_BC_Fluxes_Wall_Pyrolysis
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_Wall_Pyrolysis ( Im, Jm, Km, Fm, Blk, qrad, fuel, Ovar )
    use MOSE_Lib_Fluid
    use MOSE_Lib_GSI
    use MOSE_Lib_RANS
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm, fuel
    real(R8), intent(in) :: qrad
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    ! Local
    integer :: s, modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, Ig, Jg, Kg, iter
    real(R8) :: Normal(3), Area, Dist, dl, M(3,3), Prim(nprim), rho, Rgas, T, Gradient(nprim,3)
    real(R8) :: rho_wall, Prim_Wall(nprim), mil, kl, omega(nsc), mdot
    real(R8) :: Stress(3), Blowing(3), Sum_Omegai_Hwi, Flux(nprim)
    real(R8) :: T_wall, T_wall_old, balance, balanceold, q_conv, q_pyro, q_rad


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    ! Initialization
    Gradient = 0d0

    ! Newton-Raphson loop initialization
    T_wall = T ! First guess: T_wall=T_cell
    rho_wall = Prim(np) / ( Rgas * T_wall )         ! BL: p(wall)=p(cell), R(w)=R(c) and ci(w)=ci(c)
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho ! Partial densities at wall

    select case ( fuel )
      case ( 9, 304 ) ! HDPE / generic ablation (304 is the ATLAS qrad-wall type)
        call GSI_HDPE(T_wall, omega, mdot)
      case ( 10 ) ! PP
        call GSI_PP(T_wall, omega, mdot)
      case ( 12 ) ! Carbon
        call GSI_carbon(Prim_wall(1:nsc), T_wall, omega, mdot)
      case ( 13 ) ! HTPB
        call GSI_HTPB(T_wall, omega, mdot)
    end select
    mdot = mdot * (-modfm2) ! adjust sign
    omega = omega * (-modfm2) ! adjust sign

    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl ) ! transport

    q_conv = 0d0 ! convective heat flux
    select case ( fuel ) ! pyrolysis heat flux
    case ( 9, 304 ) ! HDPE / generic ablation (304 is the ATLAS qrad-wall type)
      q_pyro = mdot * ( Dh_HDPE + cp_HDPE* (T_wall - Ti_HDPE) ) 
    case ( 10 ) ! PP
      q_pyro = mdot * ( Dh_PP + cp_PP* (T_wall - Ti_PP) )
    case ( 12 ) ! Carbon
      q_pyro = 0d0
      do s = 1, nsc
        q_pyro = q_pyro + omega(s)*f_tabT( T_wall,s,h_tab )
      end do    
    case ( 13 ) ! HTPB
      q_pyro = mdot * ( Dh_HTPB + cp_HTPB* (T_wall - Ti_HTPB) )
    end select
    q_rad = qrad * modfm2 ! adjust sign for radiation
    balance = q_conv - q_pyro - q_rad

    do iter = 1, 20 ! start of loop
      ! first part of the iteration
      balanceold = balance
      T_wall_old = T_wall
      T_wall = T_wall * (1d0 + 1d-4) ! increment
      Gradient = 0d0
      Gradient(np,Dir) = ( T - T_wall ) * modfm3
      Gradient(np,:) = matmul ( Gradient(np,:), M )
      rho_wall = Prim(np) / ( Rgas * T_wall )
      Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho
      select case ( fuel )
      case ( 9, 304 ) ! HDPE / generic ablation (304 is the ATLAS qrad-wall type)
        call GSI_HDPE(T_wall, omega, mdot)
      case ( 10 ) ! PP
        call GSI_PP(T_wall, omega, mdot)
      case ( 12 ) ! Carbon
        call GSI_carbon(Prim_wall(1:nsc), T_wall, omega, mdot)
      case ( 13 ) ! HTPB
        call GSI_HTPB(T_wall, omega, mdot)
      end select
      mdot = mdot * (-modfm2)
      omega = omega * (-modfm2)
      call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )
      q_conv = kl * dot_product( Gradient(np,:), normal )
      select case ( fuel ) ! pyrolysis heat flux
      case ( 9, 304 ) ! HDPE / generic ablation (304 is the ATLAS qrad-wall type)
        q_pyro = mdot * ( Dh_HDPE + cp_HDPE* (T_wall - Ti_HDPE) ) 
      case ( 10 ) ! PP
        q_pyro = mdot * ( Dh_PP + cp_PP* (T_wall - Ti_PP) )
      case ( 12 ) ! Carbon
        q_pyro = 0d0
        do s = 1, nsc
          q_pyro = q_pyro + omega(s)*f_tabT( T_wall,s,h_tab )
        end do    
      case ( 13 ) ! HTPB
        q_pyro = mdot * ( Dh_HTPB + cp_HTPB* (T_wall - Ti_HTPB) )
      end select
      balance = q_conv - q_pyro - q_rad

      ! new wall temperature (Newton):
      T_wall = T_wall_old - balanceold * (1d-4*T_wall_old)/(balance-balanceold)

      !if ( Abs ( balance ) < 1d-5 ) exit ! stop criteria

      ! second part of the iteration
      Gradient = 0d0
      Gradient(np,Dir) = ( T - T_wall ) * modfm3
      Gradient(np,:) = matmul ( Gradient(np,:), M )
      rho_wall = Prim(np) / ( Rgas * T_wall )
      Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho
      select case ( fuel )
      case ( 9, 304 ) ! HDPE / generic ablation (304 is the ATLAS qrad-wall type)
        call GSI_HDPE(T_wall, omega, mdot)
      case ( 10 ) ! PP
        call GSI_PP(T_wall, omega, mdot)
      case ( 12 ) ! Carbon
        call GSI_carbon(Prim_wall(1:nsc), T_wall, omega, mdot)
      case ( 13 ) ! HTPB
        call GSI_HTPB(T_wall, omega, mdot)
      end select      
      mdot = mdot * (-modfm2)
      omega = omega * (-modfm2)
      call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )
      q_conv = kl * dot_product( Gradient(np,:), normal )
      select case ( fuel ) ! pyrolysis heat flux
      case ( 9, 304 ) ! HDPE / generic ablation (304 is the ATLAS qrad-wall type)
        q_pyro = mdot * ( Dh_HDPE + cp_HDPE* (T_wall - Ti_HDPE) ) 
      case ( 10 ) ! PP
        q_pyro = mdot * ( Dh_PP + cp_PP* (T_wall - Ti_PP) )
      case ( 12 ) ! Carbon
        q_pyro = 0d0
        do s = 1, nsc
          q_pyro = q_pyro + omega(s)*f_tabT( T_wall,s,h_tab )
        end do    
      case ( 13 ) ! HTPB
        q_pyro = mdot * ( Dh_HTPB + cp_HTPB* (T_wall - Ti_HTPB) )
      end select
      balance = q_conv - q_pyro - q_rad
    end do

    ! From here on T_wall is known
    Gradient = 0d0
    Gradient(np,Dir) = ( T - T_wall ) * modfm3
    Gradient(np,:) = matmul ( Gradient(np,:), M )
    q_conv = kl * dot_product(Gradient(np,:), normal)
    rho_wall = Prim(np) / ( Rgas * T_wall )
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho

    Blowing = mdot / rho_wall * Normal                          ! Blowing velocity
    Gradient(nu:nw,Dir) = ( Prim(nu:nw) - Blowing ) * modfm3    ! Velocity gradient due to blowing in direction normal to the face
    Gradient(nu:nw,:) = matmul ( Gradient(nu:nw,:), M )         ! Velocity gradient transformation
    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, 0d0, Prim(nt:) )

    Sum_Omegai_Hwi = 0d0
    do s = 1, nsc
      Sum_Omegai_Hwi = Sum_Omegai_Hwi + omega(s) * f_tabT( T_wall,s,h_tab )
    end do

    ! Fluxes
    Flux = 0.0d0
    Flux(1:nsc) = - omega * Area
    Flux(nu:nw) = ( Stress - mdot * Blowing ) * Area
    Flux(np) = Area * ( q_conv - Sum_Omegai_Hwi - mdot * 0.5d0 * sum( Blowing**2 ) + &
                        dot_product( Stress, Blowing ) ) ! check sign q_rad!!!

    if (model==2) then
      dl = blk % yn(im,jm,km)
      call RANS_Set_Blowing_Wall ( rho=rho_wall, mil=mil, &
                                   rans_variables=Prim_Wall(nt:nprim), &
                                   tau=stress, mdot=mdot, dist=dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
      Gradient(nt:nprim,:) = matmul ( Gradient(nt:nprim,:), M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=1d-20 )
      ! Ghost-cell extrapolation of RANS variables
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), qw=q_conv, &
                                                    mdot=mdot, y=blk%dl(im,jm,km)%c(dir)*0.5d0, &
                                                    Tw=T_wall, rhow=rho_wall, mu=mil, exit_array=Ovar)

  end subroutine BC_Wall_Pyrolysis

end module MOSE_Lib_BC_Fluxes_Wall_Pyrolysis