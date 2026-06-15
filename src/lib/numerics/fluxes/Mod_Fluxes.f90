module MOSE_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private :: Fluxes_blk
  public  :: Zero_Residuals, Internal_Fluxes

contains

  subroutine Zero_Residuals ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k, n(3)

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      n = domain % blk(b) % dim

      !$omp do collapse (3)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)
        domain % blk(b) % R(:,i,j,k) = 0.0_R8
        domain % blk(b) % beta(i,j,k) = 1.0_R8
      enddo; enddo; enddo

    enddo

  end subroutine Zero_Residuals


  subroutine Internal_Fluxes ( domain )
    !> Compute internal face fluxes only (no residual zeroing).
    !> Assumes R(:) and beta(:) have been initialized via Zero_Residuals.
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_shock_detector, obj_space_scheme, obj_riemann, obj_rans, obj_soot
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b
    logical  :: SD, SD_limiter, SD_riemann, chen_sensor, soot_enabled
    real(R8) :: Sc, Sct, Prt

    SD           = obj_shock_detector%SD
    SD_limiter   = obj_space_scheme%SD
    SD_riemann   = obj_riemann%SD
    chen_sensor  = obj_riemann%SD_chen
    soot_enabled = obj_soot%enabled
    Sc  = obj_rans%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call Fluxes_blk ( domain % blk(b), &
                        SD, SD_limiter, SD_riemann, chen_sensor, &
                        Sc, Sct, Prt, soot_enabled )
    enddo

  end subroutine Internal_Fluxes


  subroutine Fluxes_blk ( blk, SD, SD_limiter, SD_riemann, chen_sensor, Sc, Sct, Prt, soot_enabled )
    use MOSE_Advanced_Types_m, only: MOSE_block_type
    use MOSE_Global_m, only: model, gc, nprim, np
    use MOSE_Lib_Shock_Detector
    use MOSE_Lib_Convective
    use MOSE_Lib_Diffusive
    implicit none
    ! Inputs
    type(MOSE_block_type), intent(inout) :: blk
    logical, intent(in)  :: SD, SD_limiter, SD_riemann, chen_sensor, soot_enabled
    real(R8), intent(in) :: Sc, Sct, Prt
    ! Local
    integer :: i, j, k, n(3)

    n = blk % dim

    ! -----------------------------------------------------------------
    ! Shock-detector
    if (SD) then
      if (chen_sensor) then
        !$omp do collapse (3)
        do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
              blk % beta(i,j,k) = SD_Chen ( blk % P(np,i-1:i+1,j-1:j+1,k-1:k+1) )
        enddo; enddo; enddo
      else
        !$omp do collapse (3)
        do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
              blk % beta(i,j,k) = SD_Tramel ( blk % P(np,i-1:i+1,j-1:j+1,k-1:k+1) )
        enddo; enddo; enddo
      endif
    endif

    ! -----------------------------------------------------------------
    ! Convective and diffusive fluxes computation
    !$omp do collapse (2)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1) - 1
      call Convective_Flux ( blk % dl(i-1:i+2,j,k) % c(1), &
                             blk % dir(1) % f(i,j,k) % N,  &
                             blk % dir(1) % f(i,j,k) % A,  &
                             blk % P(:,i-1:i+2,j,k),       &
                             blk % R(:,i:i+1,j,k),         &
                             blk % beta(i,j,k),            &
                             SD_limiter,                   &
                             blk % Ur(i,j,k),              &
                             blk % Ur(i+1,j,k) )
    enddo; enddo; enddo

    if (model>0)  then
    !$omp do collapse (2)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1) - 1
      call Diffusive_Flux ( blk % dir(1) % f(i,j,k) % N,  &
                            blk % dir(1) % f(i,j,k) % A,  &
                            blk % yn(i  ,j,k),            &
                            blk % yn(i+1,j,k),            &
                            blk % P(:,i  ,j,k),           &
                            blk % P(:,i+1,j,k),           &
                            blk % P(:,i  ,j-1,k),         &
                            blk % P(:,i  ,j+1,k),         &
                            blk % P(:,i+1,j-1,k),         &
                            blk % P(:,i+1,j+1,k),         &
                            blk % P(:,i  ,j,k-1),         &
                            blk % P(:,i  ,j,k+1),         &
                            blk % P(:,i+1,j,k-1),         &
                            blk % P(:,i+1,j,k+1),         &
                            blk % M(i  ,j,k) % c,         &
                            blk % M(i+1,j,k) % c,         &
                            blk % R(:,i  ,j,k),           &
                            blk % R(:,i+1,j,k),           &
                            1, 2, 3,                &
                            Sc, Sct, Prt, soot_enabled)
      enddo; enddo; enddo
    endif

    !$omp do collapse (2)
    do k = 1, n(3)
    do i = 1, n(1)
    do j = 1, n(2) - 1
      call Convective_Flux ( blk % dl(i,j-1:j+2,k) % c(2), &
                             blk % dir(2) % f(i,j,k) % N,  &
                             blk % dir(2) % f(i,j,k) % A,  &
                             blk % P(:,i,j-1:j+2,k),       &
                             blk % R(:,i,j:j+1,k),         &
                             blk % beta(i,j,k),            &
                             SD_limiter,                   &
                             blk % Ur(i,j,k),              &
                             blk % Ur(i,j+1,k) )
    enddo; enddo; enddo

    if (model>0) then
      !$omp do collapse (2)
      do k = 1, n(3)
      do i = 1, n(1)
      do j = 1, n(2) - 1
      call Diffusive_Flux ( blk % dir(2) % f(i,j,k) % N,  &
                            blk % dir(2) % f(i,j,k) % A,  &
                            blk % yn(i,j  ,k),            &
                            blk % yn(i,j+1,k),            &
                            blk % P(:,i,j  ,k),           &
                            blk % P(:,i,j+1,k),           &
                            blk % P(:,i-1,j  ,k),         &
                            blk % P(:,i+1,j  ,k),         &
                            blk % P(:,i-1,j+1,k),         &
                            blk % P(:,i+1,j+1,k),         &
                            blk % P(:,i,j  ,k-1),         &
                            blk % P(:,i,j  ,k+1),         &
                            blk % P(:,i,j+1,k-1),         &
                            blk % P(:,i,j+1,k+1),         &
                            blk % M(i,j  ,k) % c,         &
                            blk % M(i,j+1,k) % c,         &
                            blk % R(:,i,j  ,k),           &
                            blk % R(:,i,j+1,k),           &
                            2, 1, 3,                &
                            Sc, Sct, Prt, soot_enabled)
      enddo; enddo; enddo
    end if

    !$omp do collapse (2)
    do j = 1, n(2)
    do i = 1, n(1)
    do k = 1, n(3) - 1
      call Convective_Flux ( blk % dl(i,j,k-1:k+2) % c(3), &
                             blk % dir(3) % f(i,j,k) % N,  &
                             blk % dir(3) % f(i,j,k) % A,  &
                             blk % P(:,i,j,k-1:k+2),       &
                             blk % R(:,i,j,k:k+1),         &
                             blk % beta(i,j,k),            &
                             SD_limiter,                   &
                             blk % Ur(i,j,k),              &
                             blk % Ur(i,j,k+1) )
    enddo; enddo; enddo

    if (model>0) then
      !$omp do collapse (2)
      do j = 1, n(2)
      do i = 1, n(1)
      do k = 1, n(3) - 1
      call Diffusive_Flux ( blk % dir(3) % f(i,j,k) % N,  &
                            blk % dir(3) % f(i,j,k) % A,  &
                            blk % yn(i,j,k  ),            &
                            blk % yn(i,j,k+1),            &
                            blk % P(:,i,j,k  ),           &
                            blk % P(:,i,j,k+1),           &
                            blk % P(:,i-1,j,k  ),         &
                            blk % P(:,i+1,j,k  ),         &
                            blk % P(:,i-1,j,k+1),         &
                            blk % P(:,i+1,j,k+1),         &
                            blk % P(:,i,j-1,k  ),         &
                            blk % P(:,i,j+1,k  ),         &
                            blk % P(:,i,j-1,k+1),         &
                            blk % P(:,i,j+1,k+1),         &
                            blk % M(i,j,k  ) % c,         &
                            blk % M(i,j,k+1) % c,         &
                            blk % R(:,i,j,k  ),           &
                            blk % R(:,i,j,k+1),           &
                            3, 1, 2,                &
                            Sc, Sct, Prt, soot_enabled )
      enddo; enddo; enddo
    endif

  end subroutine Fluxes_blk

end module MOSE_Mod_Fluxes