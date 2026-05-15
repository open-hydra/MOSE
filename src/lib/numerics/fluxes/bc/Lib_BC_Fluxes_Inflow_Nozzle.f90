module MOSE_Lib_BC_Fluxes_Inflow_Nozzle
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  public
  private :: Trapezoidal

contains

  subroutine BC_Inflow_Nozzle ( Im, Jm, Km, Fm, Blk, BC_Mach, BC_T0, BC_p0, BC_rel_fac, BC_alpha, BC_beta, BC_pAmb, BC_ci, BC_RANS, BC_mdot, error )
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: BC_Mach, BC_T0, BC_p0, BC_rel_fac, BC_alpha, BC_beta, BC_pAmb, BC_ci(nsc), BC_RANS(1:)
    real(R8), intent(out) :: BC_mdot
    type(MOSE_block_type), intent(inout) :: Blk
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k, Dir, Face_i, Face_j, Face_k
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim), Normal(3), Area, t_Vec(3), t_Mod, BC_Sign, Un
    real(R8) :: Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, Inflow_Outflow, Flux(nprim)
    real(R8) :: pSub, pSup, g, alpha_app, beta_app, BC_supersonic
    integer, intent(inout) :: error

    
    error = 0
    BC_supersonic = 0.d0

    pSub = BC_alpha * 1.d5
    pSup = BC_beta * 1.d5
    g = BC_pAmb
    alpha_app = 0.d0
    beta_app = 0.d0

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    ! (im,jm,km): boundary cell; (Int_i,Int_j,Int_k) next cell
    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)
    Int_k = Km + guide(Fm,3)

    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! boundary cell primitive variables
    Bound_Prim = Blk % P(:,Im,Jm,Km)
    Int_Prim = Blk % P(:,Int_i,Int_j,Int_k)

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    ! Other obscure metric stuff 
    Select case ( Fm )
      Case ( 1 : 2 )
        t_Vec = Blk % node(Face_i,Face_j,Face_k) % c - &
                Blk % node(Face_i,Face_j-1,Face_k-1) % c
      Case ( 3 : 4 )
        t_Vec = Blk % node(Face_i,Face_j,Face_k) % c - &
                Blk % node(Face_i-1,Face_j,Face_k-1) % c
      Case ( 5 : 6 )
        t_Vec = Blk % node(Face_i,Face_j,Face_k) % c - &
                Blk % node(Face_i-1,Face_j-1,Face_k) % c
    end select

    t_Mod = norm2 ( t_Vec )
    t_Vec = t_Vec/t_Mod

    BC_Sign = Real ( modfm2 )  ! =-1 for faces 1,3,5 ; =1 for faces 2,4,6
    
    ! Compute stuff in the boundary cell
    Un  = dot_product ( Bound_Prim(nu:nw), Normal )  ! velocity normal to the interface
    call co_rotot_Rtot ( Bound_Prim(1:nsc), Bound_rho, Bound_Rgas )
    Bound_Sound = f_ss( Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas )
    Bound_Gamma = f_gamma ( Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas ) 
    ! riemann invariant a+-delta*u directed towards the interface
    Riem = Bound_Sound + BC_Sign * 0.5d0*( Bound_Gamma - 1d0 ) * Un

    ! check for inflow/outflow. Inflow_Outflow>=0: outflow. Inflow_Outflow<0: inflow
    Inflow_Outflow = BC_Sign * Un

    if ( Inflow_Outflow >= 0d0 .or. Bound_Prim(np) >= BC_P0) then
        ! Force inflow. When the flow is directed outwards, a symmetry BC is employe
        error = 1
        return
      elseif (( Bound_Prim(np) < BC_p0 ) .and. ( Bound_Prim(np) >= pSub )) then
        call BC_Nozzle_Subsonic ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem, alpha_app, beta_app, &
                                  0d0, BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec, Normal, Area, Flux, BC_mdot )
      elseif (( Bound_Prim(np) < pSub) .and. (Bound_Prim(np) >= pSup)) then
        call BC_Nozzle_Subsonic ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem, alpha_app, beta_app, &
                                 -10.d0, BC_ci, BC_RANS, BC_T0, g, BC_rel_fac, BC_Sign, t_Vec, Normal, Area, Flux, BC_mdot )
      elseif ( Bound_Prim(np) < pSup ) then
        if (BC_Mach == 0.d0) then
          BC_supersonic = 1.d0
          call BC_Nozzle_Supersonic ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, g, BC_p0, &
                                      BC_T0, BC_ci, BC_RANS, alpha_app, beta_app, Normal, Area, Flux, BC_mdot, BC_supersonic )
        else
          call BC_Nozzle_Supersonic ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p0, &
                                      BC_T0, BC_ci, BC_RANS, alpha_app, beta_app, Normal, Area, Flux, BC_mdot, BC_supersonic )
        endif
    endif
    ! Residual update
    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + Modfm2 * Flux
    
  end subroutine BC_Inflow_Nozzle


  subroutine BC_Nozzle_Supersonic ( Prim2, Int_Prim, modfm, modfm1, modfm2, g, BC_p0, BC_T0, BC_ci, & 
                                    BC_RANS, alpha_app, beta_app, Normal, Area, Flux, Fmass, BC_supersonic )
    use MOSE_Lib_Limiters, only : rlimiter
    use MOSE_Mod_Riemann
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in) :: modfm, modfm1, modfm2
    real(R8), intent(in) :: Prim2(nprim), Int_Prim(nprim), g, BC_p0, BC_T0, BC_supersonic
    real(R8), intent(in) :: BC_ci(nsc), BC_RANS(1:), alpha_app, beta_app, Normal(3), Area
    real(R8), intent(out) :: Flux(nprim), Fmass
    ! Local
    integer :: s
    real(R8) :: Sup_Prim(nprim), Sup_T, Sup_Rgas, Sup_rho, Sup_Sound, Prim1(nprim), Prim3(nprim)
    real(R8), dimension(nprim) :: Diff32, Diff21, Slope, Face_Prim, Prim_L, Prim_R
    real(R8) :: rho_L, rho_R, Rgas_L, Rgas_R, Sound_L, Sound_R, F_r, F_u, F_v, F_w, F_E
    real(R8) :: su, Sel_L, Sel_R
    real(R8) :: h0_, cpSup, dcpSup, hSup, uSup, Sup_h0_, duSup, Sup_p, dpSup, dh0_Sup, gamSup, dgamSup, daSup
    real(R8) :: fun_T, diff_h0_, alpha, beta

    ! supersonic inflow. BC_Mach enforced at boundary, T0 and p0 total
    if ( ( alpha_app == 0d0 ) .and. ( beta_app == 0d0 ) ) then
      alpha = Atan ( Normal(2)/( Normal(1) + 1d-20 ) )
      beta  = Atan ( Normal(3)/( Normal(1) + 1d-20 ) )
    else
      alpha = alpha_app
      beta = beta_app
    endif

    h0_ = 0d0
    Sup_Rgas = 0d0
    do s = 1, nsc
      h0_ = h0_ + BC_ci(s) * f_tabT( BC_T0,s,h_tab )
      Sup_Rgas = Sup_Rgas + BC_ci(s) * Ri_tab(s)
    enddo
    Sup_T = 0.5d0*BC_T0
    diff_h0_ = 1d0

    do while (diff_h0_ >= 1.d-9)
      hSup = 0d0
      cpSup = 0d0
      dcpSup = 0d0
      do s = 1, nsc
        hSup = hSup + BC_ci(s) * f_tabT( Sup_T,s,h_tab )
        cpSup = cpSup + BC_ci(s) * f_tabT( Sup_T,s,cp_tab )
        dcpSup = dcpSup + BC_ci(s) * f_tabT( Sup_T,s,dcpi_tab )
      enddo
      if (BC_supersonic == 1.d0) then
        call Trapezoidal ( sup_T, BC_T0, BC_ci, Sup_Rgas, fun_T ) 
        Sup_p = BC_P0 / exp (fun_T)
        Sup_rho = Sup_p / ( Sup_Rgas * Sup_T )
        uSup = g / Sup_rho
        Sup_h0_ = hSup + 0.5d0 * uSup ** 2
        diff_h0_ = abs(Sup_h0_-h0_)/abs(h0_)
        dpSup = BC_P0 / exp (fun_T) * cpSup / (Sup_Rgas * Sup_T)
        duSup = g * Sup_Rgas * (Sup_p - dpSup * Sup_T) / Sup_p**2
        dh0_Sup = cpSup + uSup * duSup
      else  
        gamSup = cpSup / (cpSup - Sup_Rgas)
        Sup_Sound = sqrt(gamSup * Sup_Rgas * Sup_T)
        uSup = Sup_Sound * g
        Sup_h0_ = hSup + 0.5d0 * uSup ** 2
        diff_h0_ = abs(Sup_h0_-h0_)/abs(h0_)
        dgamSup = - Sup_Rgas * dcpSup / ( cpSup - Sup_Rgas )**2
        daSup = sqrt(Sup_Rgas) * 0.5d0 / sqrt( gamSup * Sup_T ) * ( gamSup + Sup_T * dgamSup )
        duSup = g * daSup
        dh0_Sup = cpSup + uSup * duSup
      endif
      Sup_T = Sup_T - (Sup_h0_ - h0_)/(dh0_Sup)
    enddo  
    
    ! Numerical Integration (Trapezoidal) to compute p
    call Trapezoidal ( sup_T, BC_T0, BC_ci, Sup_Rgas, fun_T )
    
    ! Compute p from isentropy (II Law of Thermodynamics)
    Sup_Prim(np) = BC_P0 / exp (fun_T)

    Sup_rho = Sup_Prim(np) / ( Sup_Rgas * Sup_T )
    Sup_Prim(1:nsc) = Sup_rho * BC_ci
    Sup_Sound = f_ss ( Sup_Prim(1:nsc), Sup_Prim(np), Sup_rho, Sup_Rgas )
   !        Sup_Prim(nu) = BC_Mach * Sup_Sound / Sqrt( 1d0 + Tan(alpha)**2 + Tan(beta)**2 )
   !        Sup_Prim(nv) = Sup_Prim(nu) * Tan(alpha)
   !        Sup_Prim(nw) = Sup_Prim(nu) * Tan(beta)
    if (BC_supersonic == 1.d0) then
      Sup_Prim(nu:nw) = g / Sup_rho * normal
    else
      Sup_Prim(nu:nw) = g * Sup_Sound * normal
    endif
    if (model==2) Sup_Prim(nt:nprim) = Sup_rho * BC_RANS

    ! building a 3 cell stencil to extrapolate the solution at the interface

    ! stencil cell 2: boundary cell

    ! stencil cell 1 is the ghost cell (odd faces) or the interior cell (even faces)
    Prim1 = modfm*Sup_Prim + modfm1*Int_Prim

    ! stencil cell 3 is the interior cell (odd faces) or the ghost cell (even faces)
    Prim3 = modfm*Int_Prim + modfm1*Sup_Prim
    
    ! The state at the interface from the interior side is reconstructed: (ro,u,p)
    Diff32 = Prim3 - Prim2
    Diff21 = Prim2 - Prim1
    do s = 1, nprim
      Slope(s) = rlimiter ( Diff32(s), Diff21(s) )
    end do
    
    ! Extrapolation at the interface from the boundary cell
    Face_Prim = Prim2 + 0.5d0 * Slope * modfm2

    ! state L of riemann problem is the ghost state for even faces and boundary reconstructed state for odd faces
    Prim_L = modfm*Sup_Prim + modfm1*Face_Prim
    call co_rotot_Rtot ( Prim_L(1:nsc), rho_L, Rgas_L )
    Sound_L = f_ss ( Prim_L(1:nsc), Prim_L(np), rho_L, Rgas_L )

    ! state R of riemann problem is the boundary reconstructed state for odd faces and the ghost state for even faces
    Prim_R = modfm*Face_Prim + modfm1*Sup_Prim
    call co_rotot_Rtot ( Prim_R(1:nsc), rho_R, Rgas_R )
    Sound_R = f_ss ( Prim_R(1:nsc), Prim_R(np), rho_R, Rgas_R )

    call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(np), Sound_L, rho_L, &
                    Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(np), Sound_R, rho_R, &
                    1d0, Normal(1), Normal(2), Normal(3), F_r, F_u, F_v, F_w, F_E)


    su = sign ( 1d0, F_r )
    Sel_L = ( 1d0+su ) / 2d0  ! state left selector
    Sel_R = ( su-1d0 ) / 2d0  ! state right selector

    ! Fluxes
    Flux = 0.0d0
    Fmass = F_r * Area
    Flux(1:nsc) = Fmass * ( Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R )
    Flux(nu) = F_u * Area
    Flux(nv) = F_v * Area
    Flux(nw) = F_w * Area
    Flux(np) = F_E * Area

    if (model==2) &
      Flux(nt:nprim) = Fmass * ( Sel_L*Prim_L(nt:nprim)/rho_L - Sel_R*Prim_R(nt:nprim)/rho_R )

  end subroutine BC_Nozzle_Supersonic
    

  subroutine BC_Nozzle_Subsonic ( Prim, Un, Sound, Gamma, Riem, alpha_app, beta_app, BC_Mach, &
                                   BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec, Normal, Area, Flux, Fmass )
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in) :: Prim(nprim), Un, Sound, Gamma, Riem, alpha_app, beta_app, BC_Mach
     real(R8), intent(in) :: BC_ci(nsc), BC_RANS(1:), BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec(3), Normal(3), Area
    real(R8), intent(out) :: Flux(nprim), Fmass
    ! Local
    integer :: s
    real(R8) :: b_Vec(3), b_Mod, alpha, beta, XA, XB, XC, XD, XE, XF, check, E1, E2
    real(R8) :: Face_Rgas, Un3
    real(R8) :: p3, Ut3, Ub3, mdot, Rgas3, cp3, T3, rho3, h3, h0_3, u3, h0_3_iter, p03, h3_iter
    real(R8) :: Face_T, Face_Prim(nprim), Face_rho, Face_Enthalpy, Face_Un
    real(R8) :: diff_h0_, da, du, dgam, ddel, dT, dp, dh0_, dp0
    real(R8) :: Rgas_acu, T_acu, h_acu, p_acu, u_acu, cp_acu, dcp_acu, gamma_acu, a_acu, T_bound, Bound_ci(nsc),fun_T, diff_T
    real(R8) :: rho_Bound
    ! Vector product n x t (why?)
    b_Vec(1) = Normal(2)*t_Vec(3) - Normal(3)*t_Vec(2)
    b_Vec(2) = Normal(3)*t_Vec(1) - Normal(1)*t_Vec(3)
    b_Vec(3) = Normal(1)*t_Vec(2) - Normal(2)*t_Vec(1)
    b_Mod = norm2 ( b_Vec )
    b_Vec = b_Vec/b_Mod

    if ( ( alpha_app == 0d0 ) .and. ( beta_app == 0d0 ) ) then
      alpha = atan ( Normal(2)/( Normal(1) + 1d-20 ) )
      beta  = atan ( Normal(3)/( Normal(1) + 1d-20 ) )
    else
      alpha = alpha_app
      beta = beta_app
    endif

    ! Obscure esotheric stuff
    XA = Normal(2) - Normal(1) * tan(alpha)
    XB =  t_Vec(2) -  t_Vec(1) * tan(alpha)
    XC =  b_Vec(2) -  b_Vec(1) * tan(alpha)
    XD = Normal(3) - Normal(1) * tan(beta)
    XE =  t_Vec(3) -  t_Vec(1) * tan(beta)
    XF =  b_Vec(3) -  b_vec(1) * tan(beta)
    check = ( XC * XE - XF * XB )
    E1 = ( XF * XA - XC * XD ) / check ! ??
    E2 = ( XD * XB - XE * XA ) / check ! ??

      ! Mach in BC input <0: Subsonic inflow with assigned mass flux and temperature (mdot=po and T=To)
      
      ! Interface state (labeled 3 for some reason)
      mdot = BC_p0
      Rgas3 = sum ( BC_ci * Ri_tab )
      call co_rotot_Rtot ( Prim(1:nsc), rho_Bound, Rgas_acu )
      diff_h0_ = 1.d0
        
      ! Option: total temperature as in BC input is Mach in BC input is ==-10, otherwise To is static temperature.
      ! NB. Heat capacity is extrapolated from the boundary cell. Exact only for gas with constant properties.
      if ( BC_Mach == -10d0 .or. BC_Mach >= 0.d0) then
        h0_3 = 0d0
        do s = 1, nsc
          h0_3 = h0_3 + BC_ci(s) * f_tabT( BC_T0,s,h_tab )
        enddo 
      endif 
        T_Bound = Prim(np)/rho_Bound/Rgas_acu
        T_acu = 1.1d0 * T_Bound
        Bound_ci = Prim(1:nsc)/rho_Bound

        do while (diff_h0_ >= 1.d-6)
          h_acu = 0.d0
          cp_acu = 0.d0
          dcp_acu = 0.d0
          do s = 1, nsc
            h_acu = h_acu + Bound_ci(s) * f_tabT( T_acu,s,h_tab )
            cp_acu = cp_acu + Bound_ci(s) * f_tabT( T_acu,s,cp_tab )
            dcp_acu = dcp_acu + Bound_ci(s) * f_tabT( T_acu,s,dcpi_tab )
          enddo
          gamma_acu = cp_acu / (cp_acu - Rgas_acu)
          a_acu = sqrt(gamma_acu * Rgas_acu * T_acu)
          ! Numerical Integration (Trapezoidal) to compute p
          call Trapezoidal ( T_Bound, T_acu, Bound_ci, Rgas_acu, fun_T )
          ! Compute p from isentropy (II Law of Thermodynamics)
          p_acu = Prim(np) * exp (fun_T)
          u_acu = (Riem - a_acu)/(BC_Sign*(gamma_acu-1)*0.5d0)
          p3 = p_acu
          u3 = u_acu
          dp = Prim(np) * exp(fun_T) * cp_acu / (Rgas_acu * T_acu)
          dgam = - Rgas_acu * dcp_acu / ( cp_acu - Rgas_acu )**2
          ddel = dgam * 0.5d0
          da = sqrt(Rgas_acu) * 0.5d0 / sqrt( gamma_acu * T_acu ) * ( gamma_acu + T_acu * dgam )
          du = (-da*BC_Sign*(gamma_acu-1)*0.5d0 - (Riem-a_acu)*(BC_Sign*ddel))/(BC_Sign*(gamma_acu-1)*0.5d0)**2
          if (BC_Mach >= 0d0) then
            diff_T = 1.d0
            h3 = h0_3 - 0.5d0 * u3**2
            T3 = 0.9*BC_T0
            do while (diff_T >= 1.d-6)
              cp3 = 0.d0
              h3_iter = 0.d0
              do s = 1, nsc
                h3_iter = h3_iter + BC_ci(s) * f_tabT( T3,s,h_tab )
                cp3 = cp3 + BC_ci(s) * f_tabT( T3,s,cp_tab )
              enddo
              diff_T = abs(h3-h3_iter)/abs(h3)
              T3 = T3 - (h3_iter - h3)/(cp3)
            enddo
            cp3 = 0.d0                  
            do s = 1, nsc
              cp3 = cp3 + BC_ci(s) * f_tabT( T3,s,cp_tab )
            enddo                 
            call Trapezoidal ( T3, BC_T0, BC_Ci, Rgas3, fun_T )
            p03 = p3 * exp(fun_T)
            diff_h0_ = abs (p03 - BC_P0)/BC_P0
            dp0 = dp * exp(fun_T) + p3 * exp(fun_T)/(Rgas3*T3)*u3*du
            T_acu = T_acu - (p03 - BC_P0)/dp0
          else
            dT = (1/Rgas3/mdot) * (dp * u3 + du * p3)
            rho3 = mdot/u3
            T3 = p3/(Rgas3*rho3)
            h3 = 0.d0
            cp3 = 0.d0
            do s = 1, nsc
              h3 = h3 + BC_ci(s) * f_tabT( T3,s,h_tab )
              cp3 = cp3 + BC_ci(s) * f_tabT( T3,s,cp_tab )
            enddo
            if ( BC_Mach == -10d0 ) then
              h0_3_iter = h3 + 0.5d0 * u3**2
              diff_h0_ = abs(h0_3 - h0_3_iter) / abs(h0_3)
              dh0_ = cp3 * dT + u3 * du
              T_acu = T_acu - (h0_3_iter - h0_3) / (dh0_)
            else
              diff_h0_ = abs(T3 - BC_T0) / BC_T0
              T_acu = T_acu - (T3 - BC_T0) / dT
            endif
          endif
        enddo
      h_acu = 0.d0
      cp_acu = 0.d0
      dcp_acu = 0.d0
      do s = 1, nsc
        h_acu = h_acu + Bound_ci(s) * f_tabT( T_acu,s,h_tab )
        cp_acu = cp_acu + Bound_ci(s) * f_tabT( T_acu,s,cp_tab )
      enddo
      gamma_acu = cp_acu / (cp_acu - Rgas_acu)
      a_acu = sqrt(gamma_acu * Rgas_acu * T_acu)
      call Trapezoidal ( T_Bound, T_acu, Bound_ci, Rgas_acu, fun_T )
      p_acu = Prim(np) * exp (fun_T)
      u_acu = (Riem - a_acu)/(BC_Sign*(gamma_acu-1)*0.5d0)
      p3 = p_acu
      u3 = u_acu
      if (BC_Mach >= 0.d0) then
        h3 = h0_3 - 0.5d0 * u3**2
        T3 = 0.9*BC_T0
        diff_T = 1.d0
        do while (diff_T >= 1.d-6)
          cp3 = 0.d0
          h3_iter = 0.d0
          do s = 1, nsc
            h3_iter = h3_iter + BC_ci(s) * f_tabT( T3,s,h_tab )
            cp3 = cp3 + BC_ci(s) * f_tabT( T3,s,cp_tab )
          enddo
          diff_T = abs(h3-h3_iter)/abs(h3)
          T3 = T3 - (h3_iter - h3)/(cp3)
        enddo
        rho3 = p3 / (Rgas3 * T3)
      else
        rho3 = mdot/u3
      endif

      Un3 = u3
      Ut3 = Un3 * E1
      Ub3 = Un3 * E2
      Face_Rgas = Rgas3
      Face_rho = rho3
      Face_Prim(1:nsc) = Face_rho * BC_ci
      Face_Prim(nu:nw) = Un3 * Normal + Ut3 * t_Vec + Ub3 * b_Vec
      Face_Prim(np) = p3
      Face_Un = Un3
      if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho
    
    Face_Enthalpy = 0d0
    Face_T = Face_Prim(np) / ( Face_Rgas * Face_rho )
    Flux = 0.0d0
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT( Face_T,s,h_tab )
      Flux(s) = Face_Prim(s) * Face_Un * Area
    end do
    Fmass = Face_rho * Face_Un * Area
    Flux(nu:nw) = Fmass * Face_Prim(nu:nw) + Face_Prim(np) * Area * Normal
    Flux(np) = Fmass * ( 0.5d0 * sum( Face_Prim(nu:nw)**2 ) + Face_Enthalpy )
    if (model==2) Flux(nt:nprim) = Face_Prim(nt:nprim) / Face_rho * Fmass

  end subroutine BC_Nozzle_Subsonic


  subroutine BC_Nozzle_Subsonic_QS ( Prim, Un, Sound, Gamma, Riem, alpha_app, beta_app, BC_Mach, &
                                     BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec, Normal, Area, Flux )
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in) :: Prim(nprim), Un, Sound, Gamma, Riem, alpha_app, beta_app, BC_Mach
    real(R8), intent(in) :: BC_ci(nsc), BC_RANS(1:), BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec(3), Normal(3), Area
    real(R8), intent(out) :: Flux(nprim)
    ! Local
    integer :: s
    real(R8) :: b_Vec(3), alpha, beta
    real(R8) :: Face_Rgas
    real(R8) :: p3, mdot, Rgas3, cp3, T3, rho3, h3, h0_3, u3
    real(R8) :: Face_T, Face_Prim(nprim), Face_rho, Face_Enthalpy, Fmass, Face_Un
    real(R8) :: diff_h0_
    real(R8) :: Rgas_acu, fun_T
    real(R8) :: rho_Bound

    if ( ( alpha_app == 0d0 ) .and. ( beta_app == 0d0 ) ) then
      alpha = Atan ( Normal(2)/( Normal(1) + 1d-20 ) )
      beta  = Atan ( Normal(3)/( Normal(1) + 1d-20 ) )
    else
      alpha = alpha_app
      beta = beta_app
    endif

    ! Mach in BC input <0: Subsonic inflow with assigned mass flux and temperature (mdot=po and T=To)

    ! Interface state (labeled 3 for some reason)
    mdot = BC_p0
    Rgas3 = Sum ( BC_ci * Ri_tab )
    call co_rotot_Rtot ( Prim(1:nsc), rho_Bound, Rgas_acu )
    diff_h0_ = 1.d0

    ! Option: total temperature as in BC input is Mach in BC input is ==-10, otherwise To is static temperature.
    ! NB. Heat capacity is extrapolated from the boundary cell. Exact only for gas with constant properties.
    if ( BC_Mach == -10d0 .or. BC_Mach >= 0.d0) then
      h0_3 = 0d0
      do s = 1, nsc
        h0_3 = h0_3 + BC_ci(s) * f_tabT( BC_T0,s,h_tab )
      enddo 
    endif 

    p3 = Prim(np)
    Face_Prim(np) = Prim(np)
    T3 = 0.9d0 * BC_T0

    do while (diff_h0_ >= 1.d-6)
      if (BC_Mach >= 0d0) then
        cp3 = 0.d0
        do s = 1, nsc
          cp3 = cp3 + BC_ci(s) * f_tabT( T3,s,cp_tab )
        enddo
        call Trapezoidal ( T3, BC_T0, BC_ci, Rgas3, fun_T )
        diff_h0_ = abs(fun_T - log(BC_p0/p3))/abs(log(BC_p0/p3))
        T3 = T3 + (fun_T - log(BC_p0/p3))/(cp3/Rgas3/T3)
        h3 = 0.d0
        do s = 1, nsc
          h3 = h3 + BC_ci(s) * f_tabT( T3,s,h_tab )
        enddo
        u3 = sqrt(2*(h0_3-h3))
      elseif (BC_Mach == -10.d0) then
        h3 = 0.d0
        cp3 = 0.d0
        do s = 1, nsc
          h3 = h3 + BC_ci(s) * f_tabT( T3,s,h_tab )
          cp3 = cp3 + BC_ci(s) * f_tabT( T3,s,cp_tab )
        enddo
          u3 = sqrt(2* (h0_3-h3))
          rho3 = p3 /(Rgas3 * T3)
          diff_h0_ = abs(mdot - rho3 * u3)/abs(mdot)
          T3 = T3 + (rho3*u3-mdot)/(sqrt(2d0) *p3 / Rgas3 * (0.5d0 * cp3 * T3 / sqrt(h0_3-h3)+sqrt(h0_3-h3))/T3**2)
          u3 = mdot / (p3/Rgas3/T3)
        else
          T3 = BC_T0
          rho3 = p3 / (Rgas3 * T3)
          u3 = mdot / rho3
          diff_h0_ = 1.d-7
      endif
    enddo

    rho3 = p3 / (Rgas3 * T3)
    Face_Rgas = Rgas3
    Face_rho = rho3
    Face_Prim(1:nsc) = Face_rho * BC_ci
    Face_Prim(nu:nw) = u3 *Normal
    Face_Un = u3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    Face_Enthalpy = 0d0
    Face_T = Face_Prim(np) / ( Face_Rgas * Face_rho )

    Flux = 0.0d0
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT( Face_T,s,h_tab )
      Flux(s) = Face_Prim(s) * Face_Un * Area
    end do

    Fmass = Face_rho * Face_Un * Area

    Flux(nu:nw) = Fmass * Face_Prim(nu:nw) + Face_Prim(np) * Area * Normal
    Flux(np) = Fmass * ( 0.5d0 * sum( Face_Prim(nu:nw)**2 ) + Face_Enthalpy )

    if (model==2) Flux(nt:nprim) = Face_Prim(nt:nprim) / Face_rho * Fmass


  end subroutine BC_Nozzle_Subsonic_QS


  subroutine Trapezoidal ( sup_T, BC_T0, BC_ci, Rgas, integral )
    use FLINT_Lib_Thermodynamic
    implicit none
    integer :: i, s
    real(R8), intent(in) :: sup_T, BC_T0, Rgas, BC_ci(nsc)
    real(R8), intent(out) :: integral
    real(R8), allocatable :: T_vec (:)
    integer :: int_sx, int_dx, n_points
    real(R8) :: cp_inf, cp_sup

    ! Compute all cps for the integration temperature interval
    if (sup_T > BC_T0) then
      int_sx = ceiling(BC_T0)
      int_dx = floor(Sup_T)
    else
      int_sx = ceiling( Sup_T )
      int_dx = floor( BC_T0 )
    endif

    n_points = int_dx-int_sx + 3
    integral = 0d0
    allocate(T_vec(n_points))

    if (sup_T > BC_T0) then
      T_vec(1) = BC_T0
      T_vec(n_points) = Sup_T
    else
      T_vec(1) = Sup_T
      T_vec(n_points) = BC_T0
    endif

    do i = 2,n_points-1
      T_vec(i) = int_sx + (i-2)
    enddo

    do i = 1,n_points-1
      cp_inf = 0d0
      cp_sup = 0d0
      do s = 1,nsc
        cp_inf = cp_inf + BC_ci(s) * f_tabT(T_vec(i),s,cp_tab )
        cp_sup = cp_sup + BC_ci(s) * f_tabT(T_vec(i+1),s,cp_tab )
      enddo
      integral = integral + 0.5d0 * (T_vec(i+1) - T_vec(i)) * (cp_inf/Rgas/T_vec(i) + cp_sup/Rgas/T_vec(i+1))
    enddo

    if (sup_T > BC_T0) integral = - integral 
  
  end subroutine Trapezoidal

end module MOSE_Lib_BC_Fluxes_Inflow_Nozzle