module MOSE_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Fluxes

contains

  subroutine Fluxes ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_shock_detector, obj_space_scheme, obj_riemann, obj_rans, obj_soot
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b
    logical  :: SD, SD_limiter, SD_riemann, soot_enabled
    real(R8) :: Sc, Sct, Prt

    SD           = obj_shock_detector%SD
    SD_limiter   = obj_space_scheme%SD
    SD_riemann   = obj_riemann%SD
    soot_enabled = obj_soot%enabled
    Sc  = obj_rans%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt

    do b = 1, domain % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle
      call Fluxes_blk ( domain % blk(b) % P,   &
                        domain % blk(b) % R,   &
                        domain % blk(b) % dir, &
                        domain % blk(b) % dl,  &
                        domain % blk(b) % beta,  &
                        domain % blk(b) % yn,  &
                        domain % blk(b) % M,   &
                        domain % blk(b) % dim,  &
                        SD, SD_limiter, SD_riemann, &
                        Sc, Sct, Prt, soot_enabled )
    enddo

  end subroutine Fluxes


  subroutine Fluxes_blk ( Prim, Res, Dir, dl, beta, yn, M, n, SD, SD_limiter, SD_riemann, &
                          Sc, Sct, Prt, soot_enabled )
    use MOSE_Global_m
    use MOSE_Base_Types_m
    use MOSE_Lib_Shock_Detector
    use MOSE_Lib_Convective
    use MOSE_Lib_Diffusive
    implicit none
    ! Inputs
    integer, intent(in)  :: n(3)
    logical, intent(in)  :: SD, SD_limiter, SD_riemann, soot_enabled
    real(R8), intent(in) :: Sc, Sct, Prt
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(out) :: Res
    type(MOSE_d_metrics_type), dimension(3), intent(in) :: Dir
    type(MOSE_vector_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: dl
    real(R8), dimension(1:n(1), 1:n(2), 1:n(3)), intent(inout) :: beta
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: yn
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    ! Local
    integer :: i, j, k

    ! -----------------------------------------------------------------
    ! Reset residuals to zero and inizialize shock sensor arrays
    !$omp do collapse (3)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)
        Res(:,i,j,k) = 0d0
        beta(i,j,k)  = 1d0
      enddo; enddo; enddo

    ! -----------------------------------------------------------------
    ! Shock-detector
    if (SD) then
      !$omp do collapse (3)
      do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
            beta(i,j,k) = SD_Tramel ( Prim(np,i-1:i+1,j-1:j+1,k-1:k+1) )
      enddo; enddo; enddo
    endif

    ! -----------------------------------------------------------------
    ! Convective and diffusive fluxes computation
    !$omp do collapse (2)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1) - 1
      call Convective_Flux ( dl(i-1:i+2,j,k) % c(1), &
                             Dir(1) % f(i,j,k) % N,  &
                             Dir(1) % f(i,j,k) % A,  &
                             Prim(:,i-1:i+2,j,k),    &
                             Res (:,i:i+1,j,k),      &
                             beta(i,j,k), SD_limiter )
    enddo; enddo; enddo

    if (model>0)  then
    !$omp do collapse (2)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1) - 1
      call Diffusive_Flux ( Dir(1) % f(i,j,k) % N,  &
                            Dir(1) % f(i,j,k) % A,  &
                            yn(i  ,j,k),            &
                            yn(i+1,j,k),            &
                            Prim(:,i  ,j,k),        &
                            Prim(:,i+1,j,k),        &
                            Prim(:,i  ,j-1,k),      &
                            Prim(:,i  ,j+1,k),      &
                            Prim(:,i+1,j-1,k),      &
                            Prim(:,i+1,j+1,k),      &
                            Prim(:,i  ,j,k-1),      &
                            Prim(:,i  ,j,k+1),      &
                            Prim(:,i+1,j,k-1),      &
                            Prim(:,i+1,j,k+1),      &
                            M(i  ,j,k) % c,         &
                            M(i+1,j,k) % c,         &
                            Res (:,i  ,j,k),        &
                            Res (:,i+1,j,k),        &
                            1, 2, 3,                &
                            Sc, Sct, Prt, soot_enabled)
      enddo; enddo; enddo
    endif

    !$omp do collapse (2)
    do k = 1, n(3)
    do i = 1, n(1)
    do j = 1, n(2) - 1
      call Convective_Flux ( dl(i,j-1:j+2,k) % c(2), &
                             Dir(2) % f(i,j,k) % N,  &
                             Dir(2) % f(i,j,k) % A,  &
                             Prim(:,i,j-1:j+2,k),    &
                             Res (:,i,j:j+1,k),      &
                             beta(i,j,k), SD_limiter )
    enddo; enddo; enddo

    if (model>0) then
      !$omp do collapse (2)
      do k = 1, n(3)
      do i = 1, n(1)
      do j = 1, n(2) - 1
      call Diffusive_Flux ( Dir(2) % f(i,j,k) % N,  &
                            Dir(2) % f(i,j,k) % A,  &
                            yn(i,j  ,k),            &
                            yn(i,j+1,k),            &
                            Prim(:,i,j  ,k),        &
                            Prim(:,i,j+1,k),        &
                            Prim(:,i-1,j  ,k),      &
                            Prim(:,i+1,j  ,k),      &
                            Prim(:,i-1,j+1,k),      &
                            Prim(:,i+1,j+1,k),      &
                            Prim(:,i,j  ,k-1),      &
                            Prim(:,i,j  ,k+1),      &
                            Prim(:,i,j+1,k-1),      &
                            Prim(:,i,j+1,k+1),      &
                            M(i,j  ,k) % c,         &
                            M(i,j+1,k) % c,         &
                            Res (:,i,j  ,k),        &
                            Res (:,i,j+1,k),        &
                            2, 1, 3,                &
                            Sc, Sct, Prt, soot_enabled)
      enddo; enddo; enddo
    end if

    !$omp do collapse (2)
    do j = 1, n(2)
    do i = 1, n(1)
    do k = 1, n(3) - 1
      call Convective_Flux ( dl(i,j,k-1:k+2) % c(3), &
                             Dir(3) % f(i,j,k) % N,  &
                             Dir(3) % f(i,j,k) % A,  &
                             Prim(:,i,j,k-1:k+2),    &
                             Res (:,i,j,k:k+1),      &
                             beta(i,j,k), SD_limiter )
    enddo; enddo; enddo

    if (model>0) then
      !$omp do collapse (2)
      do j = 1, n(2)
      do i = 1, n(1)
      do k = 1, n(3) - 1
      call Diffusive_Flux ( Dir(3) % f(i,j,k) % N,  &
                            Dir(3) % f(i,j,k) % A,  &
                            yn(i,j,k  ),            &
                            yn(i,j,k+1),            &
                            Prim(:,i,j,k  ),        &
                            Prim(:,i,j,k+1),        &
                            Prim(:,i-1,j,k  ),      &
                            Prim(:,i+1,j,k  ),      &
                            Prim(:,i-1,j,k+1),      &
                            Prim(:,i+1,j,k+1),      &
                            Prim(:,i,j-1,k  ),      &
                            Prim(:,i,j+1,k  ),      &
                            Prim(:,i,j-1,k+1),      &
                            Prim(:,i,j+1,k+1),      &
                            M(i,j,k  ) % c,         &
                            M(i,j,k+1) % c,         &
                            Res (:,i,j,k  ),        &
                            Res (:,i,j,k+1),        &
                            3, 1, 2,                &
                            Sc, Sct, Prt, soot_enabled )
      enddo; enddo; enddo
    endif

  end subroutine Fluxes_blk

end module MOSE_Mod_Fluxes