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
    
    ! ------------------------------------------------------------------
    ! Branch on IRS once outside all loops
    ! ------------------------------------------------------------------
    if ( irs_enabled ) then

      ! PHASE 1: convert residual to primitive space only
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
                                      strangcoeff )
        enddo; enddo; enddo
        !$omp end do
      enddo

      call Residual_Smoothing( domain, irs_beta )

      ! PHASE 3: apply RK stage after smoothing
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
                                        irk, n_rk, b, i, j, k, strangcoeff )
        enddo; enddo; enddo
        !$omp end do
      enddo

    endif

  contains

    ! Converts conservative residual to primitive residual in-place
    ! (no RK update – used when IRS smoothing follows)
    subroutine Compute_Prim_Residual( prim, residual, dt, volume, strangcoeff )
      implicit none
      real(R8), intent(in)    :: prim(nprim), dt, volume, strangcoeff
      real(R8), intent(inout) :: residual(nprim)
      integer  :: T_i, Tint(2), s
      real(R8) :: Tdiff, inv_volume, inv_rhovolume, rho, Rgas, temperature
      real(R8) :: dwrotot, dwE, dwroe, dwTemp, cv, eiroi
      real(R8) :: velocity(3)

      call co_rotot_Rtot( prim(1:nsc), rho, Rgas )
      temperature = prim(np) / ( rho * Rgas )

      T_i     = int( temperature )
      Tdiff   = temperature - T_i
      Tint(1) = T_i
      Tint(2) = T_i + 1

      inv_volume     = 1d0 / volume
      inv_rhovolume  = inv_volume / rho
      velocity       = prim(nu:nw)

      ! Species densities
      residual(1:nsc) = -residual(1:nsc) * inv_volume
      dwrotot = sum( residual(1:nsc) )

      ! Velocity
      residual(nu:nw) = -residual(nu:nw) * inv_rhovolume - velocity / rho * dwrotot

      ! Pressure via energy
      dwE   = -residual(np) * inv_volume
      dwroe = dwE - rho * dot_product( velocity, residual(nu:nw) ) &
                  - 0.5d0 * dwrotot * sum( velocity**2 )
      eiroi = 0d0; cv = 0d0
      do s = 1, nsc
        eiroi = eiroi + ( f_tabT_expr(s,h_tab,Tint,Tdiff) - Ri_tab(s)*temperature ) * residual(s)
        cv    = cv    + prim(s) * ( f_tabT_expr(s,cp_tab,Tint,Tdiff) - Ri_tab(s) )
      enddo
      dwTemp = ( dwroe - eiroi ) / cv

      residual(np) = sum( Ri_tab(1:nsc) * ( temperature*residual(1:nsc) + prim(1:nsc)*dwTemp ) )

      ! Turbulence
      if (nprim > np) residual(np+1:nprim) = -residual(np+1:nprim) * inv_volume

      ! Time scaling
      residual = residual * dt * strangcoeff

    end subroutine Compute_Prim_Residual

    ! Full primitive residual + RK update in one shot (no IRS)
    subroutine Compute_And_Update_Prim( prim, primold, residual, dt, volume, irk, n_rk, b, i, j, k, strangcoeff )
      implicit none
      real(R8), intent(inout) :: prim(nprim), residual(nprim)
      real(R8), intent(in)    :: primold(nprim), dt, volume, strangcoeff
      integer,  intent(in)    :: irk, n_rk, b, i, j, k

      call Compute_Prim_Residual( prim, residual, dt, volume, strangcoeff )
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
