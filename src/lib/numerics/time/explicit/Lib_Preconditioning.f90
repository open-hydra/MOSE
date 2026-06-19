module MOSE_Lib_Preconditioning
  !-----------------------------------------------------------------------------
  ! Module for preconditioning setup and helper functions.
  !
  ! The actual Weiss-Smith preconditioning for residual transformation is
  ! implemented in Lib_Newstate.f90 (Residual_Cons_to_Prim_Prec).
  ! This module provides:
  !   - INI file reading for PREC flag and preconditioning parameters
  !   - comp_Ur: unified reference velocity computation
  !
  ! Ur modes (mutually exclusive, checked in priority order):
  !   1. Ur_ref > 0:          uniform Ur = Ur_ref (sound if vel or sound > Ur_ref)
  !   2. Mach_target > 0:     Ur from target preconditioned Mach number
  !   3. default:             standard Weiss-Smith (Ur = vel, clamped)
  !
  ! After mode selection, applied in order:
  !   - Viscous floor:  Ur = max(Ur, Ur_visc)       [not applied in Ur_ref mode]
  !   - User floor:     Ur = max(Ur, Ur_min)         [if Ur_min > 0]
  !-----------------------------------------------------------------------------
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Parameters_m
  use MOSE_Config_Types_m, only: obj_prec

  implicit none
  private

  public :: obj_prec
  public :: comp_Ur
  public :: update_derived_variables
 
contains

   function comp_Ur ( vel, sound, Ur_visc ) result ( Ur )
    !---------------------------------------------------------------------------
    ! Unified reference velocity for Weiss-Smith preconditioning.
    !
    ! Modes (priority order, mutually exclusive):
    !   1. Ur_ref > 0:       Ur = Ur_ref, unless vel exceed it (then Ur = sound)
    !   2. Mach_target > 0:  exact Ur from target preconditioned Mach (quadratic)
    !                        fallback to Weiss-Smith if target not achievable
    !   3. default:          Weiss-Smith: Ur = vel, clamped by [eps_min*sound, sound]
    !
    ! Floors (applied after mode selection):
    !   - Ur_visc:  viscous floor (skipped in Ur_ref mode)
    !   - Ur_min:   user-specified global floor
    !
    ! Input:
    !   vel     - local velocity magnitude [m/s]
    !   sound   - local speed of sound [m/s]
    !   Ur_visc - viscous reference velocity [m/s] (0 if not applicable)
    !---------------------------------------------------------------------------
    implicit none
    real(R8), intent(in) :: vel, sound, Ur_visc
    real(R8) :: Ur
    ! Local
    real(R8) :: Mt, m, A_q, B_q, S_q, disc, w

    if (obj_prec%Ur_ref > 0.0_R8) then
      !-- Mode 1: Uniform reference velocity
      if (vel > obj_prec%Ur_ref) then
        Ur = sound
      else
        Ur = obj_prec%Ur_ref
      endif

    elseif (obj_prec%Mach_target > 0.0_R8) then
      !-- Mode 2: Target preconditioned Mach number (exact)
      !   Solve  M_t = vel*(1-alpha)/sqrt(alpha^2*vel^2 + Ur^2)
      !   with   alpha = 0.5*(1 - Ur^2/c^2)
      !   Setting w = (Ur/c)^2, m = vel/c gives quadratic:
      !     A*w^2 + B*w + A = 0   (A = C, roots are reciprocals)
      !   Physical root: w < 1  (Ur < sound)
      !   If no real root: M_t not achievable (M_t < Mach), fallback to Weiss-Smith
      Mt = obj_prec%Mach_target
      m  = vel / sound
      A_q = 0.25_R8 * m*m * (Mt*Mt - 1.0_R8)
      B_q = Mt*Mt - 0.5_R8 * m*m * (Mt*Mt + 1.0_R8)

      if (abs(A_q) > 1.0e-30_R8) then
        S_q  = -B_q / A_q
        disc = S_q*S_q - 4.0_R8
        if (disc > 0.0_R8 .and. S_q > 2.0_R8) then
          ! Numerically stable form: w = 2/(S + sqrt(S^2-4))
          w = 2.0_R8 / (S_q + sqrt(disc))
          Ur = sound * sqrt(w)
        else
          ! Target not achievable (M >= Mt): disable preconditioning
          Ur = sound
        endif
      else
        ! m ~ 0 (stagnation): quadratic is degenerate, use eps_min floor
        Ur = sound * obj_prec%eps_min
      endif

      Ur = max(Ur, obj_prec%eps_min * sound)    ! stagnation floor
      Ur = min(Ur, sound)                       ! supersonic cap
      Ur = max(Ur, Ur_visc)                     ! viscous floor

    else
      !-- Mode 3: Standard Weiss-Smith
      call weiss_smith_standard(vel, sound, Ur)
      Ur = max(Ur, Ur_visc)                     ! viscous floor

    endif

    !-- Global user-specified floor
    if (obj_prec%Ur_min > 0.0_R8) then
      Ur = max(Ur, obj_prec%Ur_min)
    endif

  contains

    pure subroutine weiss_smith_standard(vel, sound, Ur)
      real(R8), intent(in)  :: vel, sound
      real(R8), intent(out) :: Ur
      if (vel < obj_prec%eps_min * sound) then
        Ur = sound * obj_prec%eps_min
      elseif (vel >= sound) then
        Ur = sound
      else
        Ur = vel * obj_prec%Ur_factor
      end if
    end subroutine weiss_smith_standard

  end function comp_Ur


  subroutine Update_Derived_Variables(domain)
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Lib_RANS, only: Eddy_Viscosity
    use MOSE_Config_Types_m, only: obj_rans
    use FLINT_Lib_Thermodynamic

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, j, k, b
    logical :: process_all
    real(R8) :: vel_mag
    real(R8) :: mu_turb_loc, mu_eff_loc, k_eff_loc, dx_min, Ur_visc
    ! Scalar dissipation rate (AMC model)
    real(R8) :: eps_over_k, nu_tilde, nu_lam, chi_sa, fv1_val, S_sq, S_strain, omega_loc
    real(R8) :: Gvel(3,3), rho, rhoi(nsc), R, sound, mu, lambda, T, cp
    real(R8), parameter :: C_d = 2.0_R8, beta_star = 0.09_R8, C_mu = 0.09_R8, cv1 = 7.1_R8
    integer :: ii, jj

    process_all = .false.

    do b = 1, domain % nb
      if (.not. process_all .and. .not. is_local_block(b)) cycle
      !$omp do collapse(3) private(i, j, k, vel_mag, sound, rho, rhoi, R, mu, lambda, T, cp, &
      !$omp   mu_turb_loc, mu_eff_loc, k_eff_loc, dx_min, rho_loc, Ur_visc, &
      !$omp   eps_over_k, nu_tilde, nu_lam, chi_sa, fv1_val, S_sq, S_strain, omega_loc, Gvel, ii, jj)
      do k = 1, domain%blk(b)%dim(3)
      do j = 1, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)

        rhoi = domain%blk(b)%P(1:nsc,i,j,k)

        call co_rotot_Rtot ( rhoi, rho, R )
        T = domain%blk(b)%P(np,i,j,k)/(rho*R)

        !% Compute preconditioning reference velocity Ur
        vel_mag = norm2(domain%blk(b)%P(nu:nw, i, j, k))

        Ur_visc = 0.0_R8
        if (model>0) then
          dx_min  = minval(domain%blk(b)%dl(i,j,k)%c(1:ndir))

          cp = f_cp(domain%blk(b)%P(1:nsc,i,j,k),T,rho)
          call co_k_mi_lam_Wilke(rhoi,rho,T,mu,lambda)

          mu_eff_loc = mu
          k_eff_loc  = lambda
          if (model>1) then
            call Eddy_Viscosity(mut=mu_turb_loc, &
                                rans_variables=domain%blk(b)%P(nt:nt+nrans-1,i,j,k), &
                                mul=mu, rho=rho, &
                                vel_gradient=domain%blk(b)%vel_gradient(i,j,k)%c, &
                                walldist=domain%blk(b)%yn(i,j,k))
            mu_eff_loc = mu_eff_loc + mu_turb_loc
            k_eff_loc  = k_eff_loc  + mu_turb_loc * cp / obj_rans%Prt

            !% Scalar dissipation rate: chi = C_d * (eps/k) * Zv  (AMC model)
            eps_over_k = 0.0_R8
            if (nrans == 1) then
              ! SA: Bradshaw relations for eps/k
              nu_tilde = domain%blk(b)%P(nt, i,j,k) / rho
              nu_lam   = mu / rho
              chi_sa   = nu_tilde / nu_lam
              fv1_val  = chi_sa**3 / (chi_sa**3 + cv1**3)
              Gvel = domain%blk(b)%vel_gradient(i,j,k)%c
              S_sq = 0.0_R8
              do jj = 1, 3; do ii = 1, 3
                S_sq = S_sq + (Gvel(ii,jj) + Gvel(jj,ii))**2
              end do; end do
              S_strain = sqrt(0.5_R8 * S_sq)
              eps_over_k = fv1_val**(5.0_R8/6.0_R8) * sqrt(C_mu) * nu_tilde &
                         * S_strain / (nu_tilde + nu_lam)
            elseif (nrans >= 2) then
              ! k-omega family: eps/k = beta_star * omega
              omega_loc = domain%blk(b)%P(nt+nrans-1, i,j,k) / rho
              eps_over_k = beta_star * omega_loc
            end if
          endif

          Ur_visc = max(mu_eff_loc / (rho * dx_min), &
                        k_eff_loc  / (rho * cp * dx_min))
        endif

        sound = f_ss ( domain%blk(b)%P(1:nsc,i,j,k), domain%blk(b)%P(np,i,j,k), rho, R )
        domain%blk(b)%Ur(i,j,k) = comp_Ur(vel_mag, sound, Ur_visc)

      enddo; enddo; enddo

      !$omp end do

      !% Smooth Ur: max-dilation over face neighbours
      if (obj_prec%n_smooth_Ur > 0) &
        call Smooth_Ur_Block(domain%blk(b)%Ur, &
               domain%blk(b)%dim(1), domain%blk(b)%dim(2), domain%blk(b)%dim(3), &
               obj_prec%n_smooth_Ur)

    enddo

  end subroutine Update_Derived_Variables


  subroutine Smooth_Ur_Block(Ur, ni, nj, nk, n_smooth)
    !> Max-dilation smoothing of Ur: each cell gets max(self, 6 face-neighbours).
    !> Uses a local automatic array as double-buffer. With cell-level OMP each
    !> thread's private copy holds only its own (i,j,k); schedule(static)
    !> guarantees the same thread reads back its own entries in the copy loop.
    use MOSE_Global_m, only: gc
    implicit none
    integer,  intent(in)    :: ni, nj, nk, n_smooth
    real(R8), intent(inout) :: Ur(1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc)
    ! Local
    real(R8) :: Ur_buf(ni, nj, nk)
    integer  :: i, j, k, is

    do is = 1, n_smooth

      !$omp do collapse(3) private(i, j, k)
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        Ur_buf(i,j,k) = max( &
          Ur(i,  j,  k  ), &
          Ur(i-1,j,  k  ), Ur(i+1,j,  k  ), &
          Ur(i,  j-1,k  ), Ur(i,  j+1,k  ), &
          Ur(i,  j,  k-1), Ur(i,  j,  k+1))
      enddo; enddo; enddo
      !$omp end do

      !$omp do collapse(3) private(i, j, k)
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        Ur(i,j,k) = Ur_buf(i,j,k)
      enddo; enddo; enddo
      !$omp end do

    enddo  ! is

  end subroutine Smooth_Ur_Block
























  ! subroutine Residual_Preconditioning(prim, residual, dt, volume, dl, strangcoeff)
  !   use MOSE_Global_m
  !   use FLINT_Lib_Thermodynamic
  !   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  !   implicit none
  !   real(R8), intent(in)    :: prim(nprim), dt, volume, dl(3), strangcoeff
  !   real(R8), intent(inout) :: residual(nprim)
  !   ! Local
  !   integer  :: i, s, T_i, Tint(2)
  !   real(R8) :: rho, Rgas, temp, Tdiff, inv_rho, cp, gamma, sound
  !   real(R8) :: q2, mod_vel, rho_T, rho_p, enthalpy, ur, theta, denom
  !   real(R8) :: t1_inner, t1, t2, t3, t4, lim, log_lim, aa, bb
  !   real(R8) :: drho, drhoN, dyN
  !   real(R8) :: cons_res(5), prec_res(5), invGamma(5,5)
  !   real(R8) :: dPrim(np), Res(np)
  !   real(R8) :: dx
  !   real(R8), parameter :: val_cap = 1.0d100
  !   real(R8), dimension(nsc-1) :: Yi, drdyi, dhdyi, dyi, drhoi
  !   real(R8), dimension(3) :: vel

  !   call co_rotot_Rtot(prim(1:nsc), rho, Rgas)
  !   temp = prim(np) / (rho * Rgas)

  !   T_i = int(temp)
  !   Tdiff = temp - T_i
  !   Tint(1) = T_i
  !   Tint(2) = T_i + 1

  !   inv_rho = 1.0_R8 / rho
  !   vel = prim(nu:nw)

  !   cp = f_cp(prim(1:nsc), temp, rho)
  !   gamma = cp / (cp - Rgas)
  !   sound = sqrt(gamma * Rgas * temp)
  !   mod_vel = norm2(vel)
  !   q2 = mod_vel*mod_vel
  !   rho_T = -rho / max(temp, 1.0d-14)
  !   rho_p = 1.0_R8 / (Rgas * max(temp, 1.0d-14))
  !   enthalpy = H0(prim(np), prim(1:nsc), mod_vel)

  !   dx = min(dl(1), dl(2), dl(3))
  !   ur = comp_Ur(vel=mod_vel, sound=sound)

  !   theta = 1.0_R8 / (ur*ur) - rho_T / (rho * cp)
  !   denom = 1.0_R8 / (rho * cp * theta + rho_T)
  !   theta = max(-val_cap, min(theta, val_cap))
  !   denom = max(-val_cap, min(denom, val_cap))

  !   if (nsc == 1) then

  !     invGamma(1,1) = (rho_T*(enthalpy - q2) + rho*cp) * denom
  !     invGamma(1,2) = rho_T * vel(1) * denom
  !     invGamma(1,3) = rho_T * vel(2) * denom
  !     invGamma(1,4) = rho_T * vel(3) * denom
  !     invGamma(1,5) = -rho_T * denom

  !     invGamma(2,1) = -vel(1) / rho
  !     invGamma(2,2) = 1.0_R8 / rho
  !     invGamma(2,3) = 0.0_R8
  !     invGamma(2,4) = 0.0_R8
  !     invGamma(2,5) = 0.0_R8

  !     invGamma(3,1) = -vel(2) / rho
  !     invGamma(3,2) = 0.0_R8
  !     invGamma(3,3) = 1.0_R8 / rho
  !     invGamma(3,4) = 0.0_R8
  !     invGamma(3,5) = 0.0_R8

  !     invGamma(4,1) = -vel(3) / rho
  !     invGamma(4,2) = 0.0_R8
  !     invGamma(4,3) = 0.0_R8
  !     invGamma(4,4) = 1.0_R8 / rho
  !     invGamma(4,5) = 0.0_R8

  !     invGamma(5,1) = -(theta*(enthalpy - q2) - 1.0_R8) * denom
  !     invGamma(5,2) = -vel(1) * theta * denom
  !     invGamma(5,3) = -vel(2) * theta * denom
  !     invGamma(5,4) = -vel(3) * theta * denom
  !     invGamma(5,5) = theta * denom

  !     cons_res = -[sum(residual(1:nsc)), residual(nu:nw), residual(np)]
  !     prec_res = matmul(invGamma, cons_res)

  !     residual(nu:nw) = prec_res(2:4)
  !     residual(np) = prec_res(1)

  !     drho = rho_p * prec_res(1) + rho_T * prec_res(5)
  !     residual(1:nsc) = drho * prim(1:nsc) / rho

  !   else

  !     Yi = prim(1:nsc-1) * inv_rho
  !     do i = 1, nsc-1
  !       if (Yi(i) <= 1.0d-9) Yi(i) = 0.0_R8
  !       drdyi(i) = -rho * (Ri_tab(i) - Ri_tab(nsc)) / Rgas
  !       dhdyi(i) = f_tabT_expr(i,h_tab,Tint,Tdiff) - f_tabT_expr(nsc,h_tab,Tint,Tdiff)
  !     enddo

  !     Res = -[sum(residual(1:nsc)), residual(nu:np), residual(1:nsc-1)]
  !     Res = max(-val_cap, min(Res, val_cap))
  !     drdyi = max(-val_cap, min(drdyi, val_cap))
  !     dhdyi = max(-val_cap, min(dhdyi, val_cap))

  !     if ( .not.ieee_is_finite(cp) .or. .not.ieee_is_finite(rho) .or. .not.ieee_is_finite(rho_T) .or. &
  !          .not.ieee_is_finite(denom) .or. .not.ieee_is_finite(theta) .or. .not.all(ieee_is_finite(Res)) .or. &
  !          .not.all(ieee_is_finite(drdyi)) .or. .not.all(ieee_is_finite(dhdyi)) .or. .not.all(ieee_is_finite(Yi)) ) then
  !       write(*,*) 'Residual_Preconditioning: non-finite multispecies state'
  !       write(*,*) 'cp,rho,rho_T,theta,denom = ', cp, rho, rho_T, theta, denom
  !       write(*,*) 'temp,p,ur,Rgas = ', temp, prim(np), ur, Rgas
  !       write(*,*) 'prim species sum = ', sum(prim(1:nsc)), '  prim(np)=', prim(np)
  !       write(*,*) 'Res(1:5) = ', Res(1:5)
  !       write(*,*) 'Yi(1:min(5,nsc-1)) = ', Yi(1:min(5,nsc-1))
  !       write(*,*) 'drdyi(1:min(5,nsc-1)) = ', drdyi(1:min(5,nsc-1))
  !       write(*,*) 'dhdyi(1:min(5,nsc-1)) = ', dhdyi(1:min(5,nsc-1))
  !       error stop 'Residual_Preconditioning non-finite diagnostics'
  !     endif

  !     t1_inner = cp*(rho + sum(Yi*drdyi)) - rho_T*(sum(dhdyi*Yi) + q2 - enthalpy)

  !     lim = huge(1.0_R8)
  !     log_lim = log(lim)
  !     aa = abs(Res(1))
  !     bb = abs(t1_inner)
  !     if ( .not.ieee_is_finite(aa) .or. .not.ieee_is_finite(bb) ) then
  !       write(*,*) 'Residual_Preconditioning non-finite in t1 factors'
  !       write(*,*) 'aa, bb = ', aa, bb
  !       write(*,*) 'Res(1), t1_inner = ', Res(1), t1_inner
  !       error stop 'Residual_Preconditioning non-finite t1 factors'
  !     endif
  !     if ( aa > 0.0_R8 .and. bb > 0.0_R8 ) then
  !       if ( log(aa) + log(bb) > log_lim ) then
  !         write(*,*) 'Residual_Preconditioning overflow in t1 = Res(1)*t1_inner'
  !         write(*,*) 'Res(1), t1_inner = ', Res(1), t1_inner
  !         write(*,*) 'cp,rho,rho_T,denom = ', cp, rho, rho_T, denom
  !         error stop 'Residual_Preconditioning overflow t1'
  !       endif
  !     endif
  !     t1 = Res(1) * t1_inner

  !     t2 = 0.0_R8
  !     do i = 6, np
  !       t1_inner = rho_T*dhdyi(i-5) - cp*drdyi(i-5)
  !         aa = abs(Res(i))
  !         bb = abs(t1_inner)
  !         if ( .not.ieee_is_finite(aa) .or. .not.ieee_is_finite(bb) ) then
  !           write(*,*) 'Residual_Preconditioning non-finite in t2 factors'
  !           write(*,*) 'i, aa, bb = ', i, aa, bb
  !           write(*,*) 'Res(i), term = ', Res(i), t1_inner
  !           error stop 'Residual_Preconditioning non-finite t2 factors'
  !         endif
  !         if ( aa > 0.0_R8 .and. bb > 0.0_R8 ) then
  !           if ( log(aa) + log(bb) > log_lim ) then
  !             write(*,*) 'Residual_Preconditioning overflow in t2 species term'
  !             write(*,*) 'i, Res(i), term = ', i, Res(i), t1_inner
  !             write(*,*) 'rho_T,cp,dhdyi,drdyi = ', rho_T, cp, dhdyi(i-5), drdyi(i-5)
  !             error stop 'Residual_Preconditioning overflow t2'
  !           endif
  !       endif
  !       t2 = t2 + Res(i) * t1_inner
  !     enddo

  !     aa = abs(rho_T)
  !     bb = abs(Res(5))
  !     if ( .not.ieee_is_finite(aa) .or. .not.ieee_is_finite(bb) ) then
  !       write(*,*) 'Residual_Preconditioning non-finite in t3 factors'
  !       write(*,*) 'rho_T, Res(5) = ', rho_T, Res(5)
  !       error stop 'Residual_Preconditioning non-finite t3 factors'
  !     endif
  !     if ( aa > 0.0_R8 .and. bb > 0.0_R8 ) then
  !       if ( log(aa) + log(bb) > log_lim ) then
  !         write(*,*) 'Residual_Preconditioning overflow in t3 = -rho_T*Res(5)'
  !         write(*,*) 'rho_T, Res(5) = ', rho_T, Res(5)
  !         error stop 'Residual_Preconditioning overflow t3'
  !       endif
  !     endif
  !     t3 = -rho_T*Res(5)

  !     t4 = rho_T*sum(vel*Res(2:4))
  !     if (.not.ieee_is_finite(t4)) then
  !       write(*,*) 'Residual_Preconditioning non-finite in t4'
  !       write(*,*) 'rho_T = ', rho_T
  !       write(*,*) 'vel = ', vel
  !       write(*,*) 'Res(2:4) = ', Res(2:4)
  !       error stop 'Residual_Preconditioning non-finite t4'
  !     endif

  !     dPrim(1) = (t1 + t2 + t3 + t4) * denom
  !     dPrim(2) = -Res(1)*vel(1)*inv_rho + Res(2)*inv_rho
  !     dPrim(3) = -Res(1)*vel(2)*inv_rho + Res(3)*inv_rho
  !     dPrim(4) = -Res(1)*vel(3)*inv_rho + Res(4)*inv_rho
  !     dPrim(5) = ( Res(1)*(theta*(sum(dhdyi*Yi) - enthalpy + q2) + sum(drdyi*Yi)*inv_rho + 1.0_R8) &
  !                - sum(Res(6:np)*(theta*dhdyi + drdyi*inv_rho)) + theta*(Res(5) - sum(vel*Res(2:4))) ) * denom
  !     do i = 1, nsc-1
  !       dPrim(5+i) = (Res(i+5) - Res(1)*Yi(i)) * inv_rho
  !     enddo
  !     dPrim = max(-val_cap, min(dPrim, val_cap))

  !     dyi = dPrim(6:np)
  !     drho = rho_p*dPrim(1) + rho_T*dPrim(5) + sum(drdyi*dyi)

  !     drhoi = rho*dyi + Yi*drho
  !     dyN = -sum(dyi)
  !     drhoN = rho*dyN + prim(nsc)*inv_rho*drho

  !     residual(1:nsc-1) = drhoi
  !     residual(nsc) = drhoN
  !     residual(nu:nw) = dPrim(2:4)
  !     residual(np) = dPrim(1)

  !   endif

  !   if (nprim > np) residual(np+1:nprim) = -residual(np+1:nprim)
  !   residual = residual * dt / volume * strangcoeff

  !   ! Positivity-oriented limiter for robust turbulent frozen runs.
  !   do s = 1, nsc
  !     residual(s) = max(residual(s), -0.8d0 * prim(s))
  !   enddo
  !   residual(np) = max(residual(np), -0.8d0 * prim(np))

  ! end subroutine Residual_Preconditioning

end module MOSE_Lib_Preconditioning
