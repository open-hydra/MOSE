module MOSE_Lib_Newstate
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Newstate_Conservative, Newstate_Primitive

contains

  ! ===========================================================================
  !  CONSERVATIVE NEWSTATE
  ! ===========================================================================
  subroutine Newstate_Conservative ( domain, irk )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_time_scheme, obj_irs, obj_chemistry
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_RK
    use MOSE_Lib_Strang
    use MOSE_Lib_RANS
    use MOSE_Lib_IRS
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in)                   :: irk
    ! Local
    integer :: i, j, k, b
    integer :: n_rk
    logical :: irs_enabled
    real(R8) :: irs_beta, strangcoeff

    n_rk        = obj_time_scheme%n_RK
    irs_enabled = obj_irs%enabled
    irs_beta    = obj_irs%beta
    strangcoeff = obj_chemistry%strangcoeff

    ! ------------------------------------------------------------------
    ! Branch on IRS once outside all loops – avoids repeated evaluation
    ! ------------------------------------------------------------------
    if ( irs_enabled ) then

      ! PHASE 1: compute residuals only (no state update yet)
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Scale_Residual_Cons( domain%blk(b)%r(:,i,j,k),      &
                                    domain%blk(b)%dtlocal(i,j,k),  &
                                    domain%blk(b)%vol(i,j,k), strangcoeff )
        enddo; enddo; enddo
        !$omp end do
      enddo

      ! PHASE 2: smooth residuals (serial across blocks by design)
      call Residual_Smoothing( domain, irs_beta )

      ! PHASE 3: update state after smoothing
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Update_State_Cons_IRS( domain%blk(b)%P(:,i,j,k),   &
                                      domain%blk(b)%PO(:,i,j,k),  &
                                      domain%blk(b)%r(:,i,j,k),   &
                                      irk, n_rk, b, i, j, k )
        enddo; enddo; enddo
        !$omp end do
      enddo

    else

      ! No IRS: one pass – compute and update in same loop
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_Cons_NoIRS( domain%blk(b)%P(:,i,j,k),    &
                                   domain%blk(b)%PO(:,i,j,k),   &
                                   domain%blk(b)%r(:,i,j,k),    &
                                   domain%blk(b)%dtlocal(i,j,k),&
                                   domain%blk(b)%vol(i,j,k),    &
                                   irk, n_rk, b, i, j, k, strangcoeff )
        enddo; enddo; enddo
        !$omp end do
      enddo

    endif

  contains

    ! Scale residual in-place for conservative integration (IRS path)
    subroutine Scale_Residual_Cons( residual, dt, volume, strangcoeff )
      implicit none
      real(R8), intent(inout) :: residual(nprim)
      real(R8), intent(in)    :: dt, volume, strangcoeff
      residual = -residual / volume * dt * strangcoeff
    end subroutine Scale_Residual_Cons

    ! Apply RK stage + cons2prim after IRS smoothing
    subroutine Update_State_Cons_IRS( prim, primold, residual, irk, n_rk, b, i, j, k )
      implicit none
      real(R8), intent(inout) :: prim(nprim), residual(nprim)
      real(R8), intent(in)    :: primold(nprim)
      integer,  intent(in)    :: irk, n_rk, b, i, j, k
      real(R8) :: rho, Rgas, temperature, consold(nprim), cons(nprim)

      consold(1:np) = prim2cons( primold(1:np) )
      cons   (1:np) = prim2cons( prim   (1:np) )
      if (nprim > np) then
        consold(np+1:nprim) = primold(np+1:nprim)
        cons   (np+1:nprim) = prim   (np+1:nprim)
      endif

      cons = RK_stage( irk, n_rk, cons, consold, residual )

      ! Reuse current prim pressure/density for temperature guess (cheap)
      call co_rotot_Rtot( prim(1:nsc), rho, Rgas )
      temperature = prim(np) / ( rho * Rgas )

      prim(1:np) = cons2prim( cons(1:np), temperature )
      if (nprim > np) prim(np+1:nprim) = cons(np+1:nprim)

      call check_and_fix_state( prim, b, i, j, k )

    end subroutine Update_State_Cons_IRS

    ! Full compute+update without IRS
    subroutine Compute_Cons_NoIRS( prim, primold, residual, dt, volume, irk, n_rk, b, i, j, k, strangcoeff )
      implicit none
      real(R8), intent(inout) :: prim(nprim), residual(nprim)
      real(R8), intent(in)    :: primold(nprim), dt, volume, strangcoeff
      integer,  intent(in)    :: irk, n_rk, b, i, j, k
      real(R8) :: rho, Rgas, temperature, consold(nprim), cons(nprim)

      residual = -residual / volume * dt * strangcoeff

      consold(1:np) = prim2cons( primold(1:np) )
      cons   (1:np) = prim2cons( prim   (1:np) )
      if (nprim > np) then
        consold(np+1:nprim) = primold(np+1:nprim)
        cons   (np+1:nprim) = prim   (np+1:nprim)
      endif

      cons = RK_stage( irk, n_rk, cons, consold, residual )

      call co_rotot_Rtot( prim(1:nsc), rho, Rgas )
      temperature = prim(np) / ( rho * Rgas )

      prim(1:np) = cons2prim( cons(1:np), temperature )
      if (nprim > np) prim(np+1:nprim) = cons(np+1:nprim)

      call Check_And_Fix_State( prim, b, i, j, k )

    end subroutine Compute_Cons_NoIRS

  end subroutine Newstate_Conservative


  ! ===========================================================================
  !  PRIMITIVE NEWSTATE
  ! ===========================================================================
  subroutine Newstate_Primitive ( domain, irk )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_time_scheme, obj_irs, obj_chemistry
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_RK
    use MOSE_Lib_Strang
    use MOSE_Lib_RANS
    use MOSE_Lib_IRS
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer, intent(in)                   :: irk
    ! Local
    integer :: i, j, k, b
    integer :: n_rk
    logical :: irs_enabled
    real(R8) :: irs_beta, strangcoeff

    n_rk        = obj_time_scheme%n_RK
    irs_enabled = obj_irs%enabled
    irs_beta    = obj_irs%beta
    strangcoeff = obj_chemistry%strangcoeff
    
    if ( irs_enabled ) then

    ! ------------------------------------------------------------------
    ! IRS path: compute residuals in conservative space, smooth, then cons2prim + RK update
    ! ------------------------------------------------------------------

      ! PHASE 1: compute residuals only (no state update yet)
      call Residual_Smoothing( domain, irs_beta )

      ! PHASE 2: convert residual to primitive space only
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_Prim_Residual( domain%blk(b)%P(:,i,j,k),    &
                                      domain%blk(b)%r(:,i,j,k),    &
                                      domain%blk(b)%dtlocal(i,j,k),&
                                      domain%blk(b)%vol(i,j,k),    &
                                      domain%blk(b)%Ur(i,j,k),     &
                                      strangcoeff )
        enddo; enddo; enddo
        !$omp end do
      enddo

      ! PHASE 3: apply RK update
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          domain%blk(b)%P(:,i,j,k) = RK_stage( irk, n_rk,                 &
                                               domain%blk(b)%P(:,i,j,k),  &
                                               domain%blk(b)%PO(:,i,j,k), &
                                               domain%blk(b)%r(:,i,j,k)   )
          call Check_And_Fix_State( domain%blk(b)%P(:,i,j,k), b, i, j, k )
        enddo; enddo; enddo
        !$omp end do
      enddo

    else

    ! ------------------------------------------------------------------
    ! No IRS: compute primitive residual + RK update in one shot
    ! ------------------------------------------------------------------

      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(3) schedule(static) private(i,j,k)
        do k = 1, domain%blk(b)%dim(3)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_And_Update_Prim( domain%blk(b)%P(:,i,j,k),    &
                                        domain%blk(b)%PO(:,i,j,k),   &
                                        domain%blk(b)%r(:,i,j,k),    &
                                        domain%blk(b)%dtlocal(i,j,k),&
                                        domain%blk(b)%vol(i,j,k),    &
                                        domain%blk(b)%Ur(i,j,k),     &
                                        irk, n_rk, b, i, j, k, strangcoeff )
        enddo; enddo; enddo
        !$omp end do
      enddo

    endif

  contains

    subroutine Compute_Prim_Residual( prim, residual, dt, volume, Ur, strangcoeff )

      implicit none

      real(R8), intent(in)    :: prim(nprim)
      real(R8), intent(in)    :: dt, volume, Ur, strangcoeff
      real(R8), intent(inout) :: residual(nprim)

      integer  :: T_i, Tint(2), i

      real(R8) :: rho, Rgas, temperature, Tdiff
      real(R8) :: inv_rho, velocity(3), modvel, q2
      real(R8) :: cp, gam, sound, a2
      real(R8) :: rho_T, rho_p, enthalpy

      ! Weiss-Smith preconditioning
      real(R8) :: beta2, Theta, Denom
      real(R8) :: Cons_Res(5), Prec_Res(5), invGamma(5,5), drho

      ! Multi-species coupling
      real(R8), dimension(nsc-1) :: drdyi, dhdyi, Yi, dyi, drhoi
      real(R8) :: drhoN, dyN
      real(R8), dimension(np)    :: Res, dPrim_appo

      call co_rotot_Rtot( prim(1:nsc), rho, Rgas )

      temperature = prim(np)/(rho*Rgas)

      T_i     = int(temperature)
      Tdiff   = temperature - T_i

      Tint(1) = T_i
      Tint(2) = T_i + 1

      inv_rho  = 1.d0/rho
      velocity = prim(nu:nw)
      modvel   = norm2(velocity)
      q2       = modvel*modvel

      !----------------------------------------------------------
      ! Mixture thermodynamics
      !----------------------------------------------------------

      cp       = f_cp_expr( prim(1:nsc), Tint, Tdiff, rho )
      gam      = cp/(cp - Rgas)
      sound    = sqrt( gam*Rgas*temperature )
      a2       = sound*sound

      rho_T    = -rho/temperature           ! (drho/dT)_p  physical
      rho_p    = 1.d0/(Rgas*temperature)     ! (drho/dp)_T  physical
      enthalpy = H0( prim(np), prim(1:nsc), modvel )

      !----------------------------------------------------------
      ! Weiss-Smith preconditioning parameters
      !   Theta = (drho/dp) replaced by pseudo-acoustic value
      !----------------------------------------------------------

      !beta2     = max( Ur*Ur/a2, 1.d-6 )
      beta2     = Ur*Ur/a2

      Theta = 1.d0/(beta2*a2) + 1.d0/(cp*temperature)
      Denom = 1.d0/( rho*cp*Theta + rho_T )

      !----------------------------------------------------------
      ! Apply the full preconditioned cons->prim transform.
      ! Velocity rows are the plain chain rule; pressure (and,
      ! through it, density) carry the preconditioning via Theta.
      !----------------------------------------------------------

      if ( nsc == 1 ) then

        ! invGamma : (drho, drhoV, drhoE) -> (dp, du, dv, dw, dT)
        invGamma(1,1) = ( rho_T*(enthalpy - q2) + rho*cp )*Denom
        invGamma(1,2) = rho_T*velocity(1)*Denom
        invGamma(1,3) = rho_T*velocity(2)*Denom
        invGamma(1,4) = rho_T*velocity(3)*Denom
        invGamma(1,5) = -rho_T*Denom

        invGamma(2,1) = -velocity(1)*inv_rho
        invGamma(2,2) =  inv_rho
        invGamma(2,3) =  0.d0
        invGamma(2,4) =  0.d0
        invGamma(2,5) =  0.d0

        invGamma(3,1) = -velocity(2)*inv_rho
        invGamma(3,2) =  0.d0
        invGamma(3,3) =  inv_rho
        invGamma(3,4) =  0.d0
        invGamma(3,5) =  0.d0

        invGamma(4,1) = -velocity(3)*inv_rho
        invGamma(4,2) =  0.d0
        invGamma(4,3) =  0.d0
        invGamma(4,4) =  inv_rho
        invGamma(4,5) =  0.d0

        invGamma(5,1) = -( Theta*(enthalpy - q2) - 1.d0 )*Denom
        invGamma(5,2) = -velocity(1)*Theta*Denom
        invGamma(5,3) = -velocity(2)*Theta*Denom
        invGamma(5,4) = -velocity(3)*Theta*Denom
        invGamma(5,5) =  Theta*Denom

        Cons_Res = -[ sum(residual(1:nsc)), residual(nu:nw), residual(np) ]
        Prec_Res = matmul( invGamma, Cons_Res )      ! (dp, dV, dT)

        residual(nu:nw) = Prec_Res(2:4)
        residual(np)    = Prec_Res(1)

        ! density from preconditioned dp and dT (physical EOS derivatives)
        residual(1)     = rho_p*Prec_Res(1) + rho_T*Prec_Res(5)

      else

        !----------------------------------------------------------
        ! Multi-species: T and Rgas depend on composition, so the
        ! mass-fraction sensitivities enter the pressure/temperature
        ! rows. Solved algebraically (same structure as nsc==1).
        !----------------------------------------------------------

        Yi = prim(1:nsc-1)*inv_rho
        do i = 1, nsc-1
          if ( Yi(i) <= 1.d-9 ) Yi(i) = 0.d0
          drdyi(i) = -rho*( Ri_tab(i) - Ri_tab(nsc) )/Rgas
          dhdyi(i) =  f_tabT_expr(i,h_tab,Tint,Tdiff) - f_tabT_expr(nsc,h_tab,Tint,Tdiff)
        end do

        ! Res = (total drho, drhoU, drhoV, drhoW, drhoE, drho_1..drho_{nsc-1})
        Res = -[ sum(residual(1:nsc)), residual(nu:np), residual(1:nsc-1) ]

        dPrim_appo(1) = ( Res(1)*( cp*(rho + sum(Yi*drdyi))                         &
                                 - rho_T*( sum(dhdyi*Yi) + q2 - enthalpy ) )        &
                        + sum( Res(6:np)*(rho_T*dhdyi - cp*drdyi) )                 &
                        - rho_T*Res(5) + rho_T*sum(velocity*Res(2:4)) )*Denom
        dPrim_appo(2) = -Res(1)*velocity(1)*inv_rho + Res(2)*inv_rho
        dPrim_appo(3) = -Res(1)*velocity(2)*inv_rho + Res(3)*inv_rho
        dPrim_appo(4) = -Res(1)*velocity(3)*inv_rho + Res(4)*inv_rho
        dPrim_appo(5) = ( Res(1)*( Theta*( sum(dhdyi*Yi) - enthalpy + q2 )          &
                                 + sum(drdyi*Yi)*inv_rho + 1.d0 )                   &
                        - sum( Res(6:np)*( Theta*dhdyi + drdyi*inv_rho ) )          &
                        + Theta*( Res(5) - sum(velocity*Res(2:4)) ) )*Denom
        do i = 1, nsc-1
          dPrim_appo(5+i) = ( Res(i+5) - Res(1)*Yi(i) )*inv_rho
        end do

        dyi   = dPrim_appo(6:np)
        drho  = rho_p*dPrim_appo(1) + rho_T*dPrim_appo(5) + sum(drdyi*dyi)

        drhoi = rho*dyi + Yi*drho
        dyN   = -sum(dyi)
        drhoN = rho*dyN + prim(nsc)*inv_rho*drho

        residual(1:nsc-1) = drhoi
        residual(nsc)     = drhoN
        residual(nu:nw)   = dPrim_appo(2:4)
        residual(np)      = dPrim_appo(1)

      endif

      ! Turbulence variables (not preconditioned)
      if (nprim > np) residual(np+1:nprim) = -residual(np+1:nprim)

      ! Time scaling
      residual = residual * dt / volume * strangcoeff

    end subroutine Compute_Prim_Residual

    ! Full primitive residual + RK update in one shot (no IRS)
    subroutine Compute_And_Update_Prim( prim, primold, residual, dt, volume, Ur, irk, n_rk, b, i, j, k, strangcoeff )
      implicit none
      real(R8), intent(inout) :: prim(nprim), residual(nprim)
      real(R8), intent(in)    :: primold(nprim), dt, volume, Ur, strangcoeff
      integer,  intent(in)    :: irk, n_rk, b, i, j, k

      call Compute_Prim_Residual( prim, residual, dt, volume, Ur, strangcoeff )
      prim = RK_Stage( irk, n_rk, prim, primold, residual )
      call Check_And_Fix_State( prim, b, i, j, k )

    end subroutine Compute_And_Update_Prim

  end subroutine Newstate_Primitive


  ! ===========================================================================
  !  STATE SANITY CHECK AND FLOOR ENFORCEMENT
  ! ===========================================================================
  subroutine Check_And_Fix_State ( prim, b, i, j, k )
    use MOSE_Global_m
    use MOSE_Lib_RANS

    implicit none
    real(R8), dimension(:), intent(inout) :: prim
    integer, intent(in) :: b, i, j, k
    integer :: s

    ! Bail out on NaN or negative pressure (product trick is expensive; use
    ! any(isnan()) + separate pressure check instead)
    if ( any( isnan(prim(1:nprim)) ) .or. prim(np) < 0d0 ) then
      write(*,'(A,4I4)') "Integration failed at b, i, j, k:", b, i, j, k
      write(*,*) prim(1:nprim)
      stop "NaN or p<0 detected"
    endif

    ! Floor species densities
    do s = 1, nsc
      if ( prim(s) < 1d-20 ) prim(s) = 1d-20
    enddo

    ! Realizability for turbulence
    if (model==2) call RANS_Enforce_Realizability( prim(nt:) )

    ! Floor soot variables
    if (nsoot == 2) then
      if ( prim(nc)   < 1d-20 ) prim(nc)   = 1d-20
      if ( prim(nc+1) < 1d+10 ) prim(nc+1) = 1d+10
    elseif (nsoot == 1) then
      if ( prim(nc)   < 1d+10 ) prim(nc)   = 1d+10
    endif

  end subroutine Check_And_Fix_State


end module MOSE_Lib_Newstate
