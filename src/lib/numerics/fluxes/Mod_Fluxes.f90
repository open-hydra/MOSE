module MOSE_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private :: Fluxes_blk
  public  :: Zero_Residuals, Internal_Fluxes, Setup_Quasi2D_W, Zero_W_Residual

contains

  subroutine Zero_Residuals ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, j, k, n(3)

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      n = domain % blk(b) % dim

      ! R(:,1:n1,j,k) is contiguous: zeroing it as one array store per (j,k)
      ! is a single large memset instead of n1 short ones.
      !$omp do collapse (2)
      do k = 1, n(3)
      do j = 1, n(2)
        domain % blk(b) % R(:,1:n(1),j,k)  = 0.0_R8
        domain % blk(b) % beta(1:n(1),j,k) = 1.0_R8
      enddo; enddo

    enddo

  end subroutine Zero_Residuals


  !> Planar quasi-2D case (one cell in k, ndir = 2, no wedge angle) whose k faces are all symmetry planes
  !> (300): the two mirror planes one cell apart give w the mode dw/dt = -2 a w / dz (a the speed of
  !> sound), far outside the stability region of the explicit schemes at the time step of the in-plane
  !> cells, so round-off in w grows until the run fails. As in planar codes, which do not solve the third
  !> direction, the z velocity is then set to zero and the z momentum is not integrated: w = 0 at every
  !> step, also when the initial field has w /= 0 (a field restarted from a run hit by the mode, or a 3-D
  !> field), which a frozen rho w would carry as an error on p and u. Null (0) or periodic (200) k faces,
  !> and axisymmetric wedges (delthe /= 0, where the cell-centre w of a cell off the wedge mid-plane
  !> carries part of the radial velocity), are left as they are. The k faces must also be planes normal to z
  !> (the x and y components of the unit normal of every k face at most tol_plane): the planar flow has
  !> w = 0. A layer whose thickness varies in the plane (k faces not normal to z) gives the gas a z velocity
  !> next to an inclined face, and freezing it would change the flow (the mass flow, the pressure): w is then
  !> integrated, with a warning, and its mode across the layer can limit the time step.
  !> Called after the initial field is in place and before the ghost cells and the coarse multigrid levels
  !> are filled from it.
  subroutine Setup_Quasi2D_W ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m,    only: ndir, freeze_w, nw
    use MOSE_Lib_Metrics, only: delthe
    use MOSE_Mod_MPI,     only: mpi_is_root, is_local_block, mpi_allreduce_max_r8
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    real(R8), parameter :: tol_plane = 1.0e-10_R8   ! Largest off-axis component of the unit normal of a face of a plane
    integer  :: i, j, k, b, nk, nsym, n(3)
    real(R8) :: offk, offk_all

    nk = 0
    nsym = 0
    do i = 1, domain % nbound
      if ( domain % bc(i) % f < 5 ) cycle
      nk = nk + 1
      if ( domain % bc(i) % type == 300 ) nsym = nsym + 1
    enddo
    ! No wedge: delthe exactly zero, written as two ordered comparisons (same result, also for -0 and NaN)
    freeze_w = ( ndir == 2 .and. delthe >= 0d0 .and. delthe <= 0d0 .and. nk > 0 .and. nsym == nk )

    ! The k faces must be planes normal to z; the largest value over the blocks of every rank, so that every rank takes
    ! the same decision
    if ( freeze_w ) then
      offk = 0.0_R8
      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        n = domain % blk(b) % dim
        do k = 0, n(3); do j = 1, n(2); do i = 1, n(1)
          offk = max ( offk, abs(domain % blk(b) % dir(3) % f(i,j,k) % N(1)), abs(domain % blk(b) % dir(3) % f(i,j,k) % N(2)) )
        enddo; enddo; enddo
      enddo
      call mpi_allreduce_max_r8 ( offk, offk_all )
      if ( .not. offk_all <= tol_plane ) then
        if ( mpi_is_root ) write(*,'(A,ES9.2,A)') &
          ' [WARNING] Quasi-2D case with symmetry (300) k faces that are not planes normal to z (unit normal off z by up to', &
          offk_all, '): the z velocity is integrated, and its mode across the layer can limit the time step.'
        freeze_w = .false.
      endif
    endif

    if ( .not. freeze_w ) return

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      domain % blk(b) % P(nw,:,:,:) = 0.0_R8
    enddo

    if ( mpi_is_root ) write(*,'(A)') &
      ' [WARNING] Quasi-2D case with symmetry (300) k faces: the z velocity is set to zero and not integrated.'

  end subroutine Setup_Quasi2D_W


  !> With freeze_w, zero the z-momentum residual before the stage update (see Setup_Quasi2D_W).
  subroutine Zero_W_Residual ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: nw
    use MOSE_Mod_MPI,  only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, j, k, n(3)

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      n = domain % blk(b) % dim

      !$omp do collapse (2)
      do k = 1, n(3)
      do j = 1, n(2)
        domain % blk(b) % R(nw,1:n(1),j,k) = 0.0_R8
      enddo; enddo

    enddo

  end subroutine Zero_W_Residual


  subroutine Internal_Fluxes ( domain )
    !> Compute internal face fluxes only (no residual zeroing).
    !> Assumes R(:) and beta(:) have been initialized via Zero_Residuals.
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_shock_detector, obj_rans, obj_soot, obj_sim_param, obj_space_scheme
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer  :: b
    integer  :: SD_id
    logical  :: soot_enabled, inverse_metric
    real(R8) :: Sc, Sct, Prt, Prl

    SD_id = obj_shock_detector%id
    soot_enabled = obj_soot%enabled
    inverse_metric = obj_space_scheme%inverse_metric
    Sc  = obj_sim_param%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt
    Prl = obj_sim_param%Prl

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call Fluxes_blk ( b, domain % blk(b), SD_id, Sc, Sct, Prt, Prl, soot_enabled, inverse_metric )
    enddo

  end subroutine Internal_Fluxes


  subroutine Fluxes_blk ( b, blk, SD_id, Sc, Sct, Prt, Prl, soot_enabled, inverse_metric )
    use MOSE_Advanced_Types_m, only: MOSE_block_type
    use MOSE_Base_Types_m, only: MOSE_tensor_3D_type
    use MOSE_Global_m, only: model, gc, nprim, np
    use MOSE_Lib_Shock_Detector
    use MOSE_Lib_Convective
    use MOSE_Lib_Diffusive
    implicit none
    ! Inputs
    integer, intent(in)  :: b                ! block id, for error messages only
    type(MOSE_block_type), intent(inout), target :: blk
    logical, intent(in)  :: soot_enabled
    logical, intent(in)  :: inverse_metric   ! diffusive-metric = inverse-mean: face metric from blk%Minv
    integer, intent(in)  :: SD_id
    real(R8), intent(in) :: Sc, Sct, Prt, Prl
    ! Local
    integer :: i, j, k, n(3)
    !> Four-cell stencil of cell lengths for the Convective_Flux call.  Filled
    !> by explicit stores: `blk % dl(i-1:i+2,j,k) % c(d)` is a strided section
    !> of a derived type, which the compiler would copy to a temporary on every
    !> face.  Per-thread already — Fluxes_blk runs inside an orphaned parallel
    !> region.
    real(R8) :: dl4(-1:2)
    !> Cell metrics handed to Diffusive_Flux: the inverse metrics (inverse-mean)
    !> or the metrics themselves (mean).
    type(MOSE_tensor_3D_type), pointer :: Mc(:,:,:)

    n = blk % dim
    if ( model > 0 ) then
      if ( inverse_metric ) then
        Mc => blk % Minv
      else
        Mc => blk % M
      end if
    end if

    ! -----------------------------------------------------------------
    ! Shock-detector
    select case(SD_id)
    case(1)
      !$omp do collapse (3)
      do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
            blk % beta(i,j,k) = SD_Tramel ( blk % P(np,i-1:i+1,j-1:j+1,k-1:k+1) )
      enddo; enddo; enddo
    case(2)
      !$omp do collapse (3)
      do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
            blk % beta(i,j,k) = SD_Chen ( blk % P(np,i-1:i+1,j-1:j+1,k-1:k+1) )
      enddo; enddo; enddo
    end select

    ! -----------------------------------------------------------------
    ! Convective and diffusive fluxes computation.
    !
    ! The directional face loops take schedule(runtime): the split is tunable
    ! from OMP_SCHEDULE (`guided` is worth trying on chemistry-imbalanced
    ! meshes) and defaults to static when it is unset.
    !$omp do collapse (2) schedule(runtime)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1) - 1
      dl4(-1) = blk % dl(i-1,j,k) % c(1) ; dl4(0) = blk % dl(i  ,j,k) % c(1)
      dl4( 1) = blk % dl(i+1,j,k) % c(1) ; dl4(2) = blk % dl(i+2,j,k) % c(1)
      call Convective_Flux ( dl4,                          &
                             blk % dir(1) % f(i,j,k) % N,  &
                             blk % dir(1) % f(i,j,k) % A,  &
                             blk % P(:,i-1:i+2,j,k),       &
                             blk % R(:,i:i+1,j,k),         &
                             blk % beta(i,j,k),            &
                             blk % Ur(i,j,k),              &
                             blk % Ur(i+1,j,k),            &
                             [b, 1, i, j, k] )
    enddo; enddo; enddo

    if (model>0)  then
    !$omp do collapse (2) schedule(runtime)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1) - 1
      call Diffusive_Flux ( blk % dir(1) % f(i,j,k) % N,  &
                            blk % dir(1) % f(i,j,k) % A,  &
                            blk % yn(i  ,j,k),            &
                            blk % yn(i+1,j,k),            &
                            blk % k_rough(i  ,j,k),       &
                            blk % k_rough(i+1,j,k),       &
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
                            Mc(i  ,j,k) % c,              &
                            Mc(i+1,j,k) % c,              &
                            blk % R(:,i  ,j,k),           &
                            blk % R(:,i+1,j,k),           &
                            1, 2, 3,                &
                            Sc, Sct, Prt, Prl, soot_enabled, inverse_metric)
      enddo; enddo; enddo
    endif

    !$omp do collapse (2) schedule(runtime)
    do k = 1, n(3)
    do i = 1, n(1)
    do j = 1, n(2) - 1
      dl4(-1) = blk % dl(i,j-1,k) % c(2) ; dl4(0) = blk % dl(i,j  ,k) % c(2)
      dl4( 1) = blk % dl(i,j+1,k) % c(2) ; dl4(2) = blk % dl(i,j+2,k) % c(2)
      call Convective_Flux ( dl4,                          &
                             blk % dir(2) % f(i,j,k) % N,  &
                             blk % dir(2) % f(i,j,k) % A,  &
                             blk % P(:,i,j-1:j+2,k),       &
                             blk % R(:,i,j:j+1,k),         &
                             blk % beta(i,j,k),            &
                             blk % Ur(i,j,k),              &
                             blk % Ur(i,j+1,k),            &
                             [b, 2, i, j, k] )
    enddo; enddo; enddo

    if (model>0) then
      !$omp do collapse (2) schedule(runtime)
      do k = 1, n(3)
      do i = 1, n(1)
      do j = 1, n(2) - 1
      call Diffusive_Flux ( blk % dir(2) % f(i,j,k) % N,  &
                            blk % dir(2) % f(i,j,k) % A,  &
                            blk % yn(i,j  ,k),            &
                            blk % yn(i,j+1,k),            &
                            blk % k_rough(i,j  ,k),       &
                            blk % k_rough(i,j+1,k),       &
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
                            Mc(i,j  ,k) % c,              &
                            Mc(i,j+1,k) % c,              &
                            blk % R(:,i,j  ,k),           &
                            blk % R(:,i,j+1,k),           &
                            2, 1, 3,                &
                            Sc, Sct, Prt, Prl, soot_enabled, inverse_metric)
      enddo; enddo; enddo
    end if

    !$omp do collapse (2) schedule(runtime)
    do j = 1, n(2)
    do i = 1, n(1)
    do k = 1, n(3) - 1
      dl4(-1) = blk % dl(i,j,k-1) % c(3) ; dl4(0) = blk % dl(i,j,k  ) % c(3)
      dl4( 1) = blk % dl(i,j,k+1) % c(3) ; dl4(2) = blk % dl(i,j,k+2) % c(3)
      call Convective_Flux ( dl4,                          &
                             blk % dir(3) % f(i,j,k) % N,  &
                             blk % dir(3) % f(i,j,k) % A,  &
                             blk % P(:,i,j,k-1:k+2),       &
                             blk % R(:,i,j,k:k+1),         &
                             blk % beta(i,j,k),            &
                             blk % Ur(i,j,k),              &
                             blk % Ur(i,j,k+1),            &
                             [b, 3, i, j, k] )
    enddo; enddo; enddo

    if (model>0) then
      !$omp do collapse (2) schedule(runtime)
      do j = 1, n(2)
      do i = 1, n(1)
      do k = 1, n(3) - 1
      call Diffusive_Flux ( blk % dir(3) % f(i,j,k) % N,  &
                            blk % dir(3) % f(i,j,k) % A,  &
                            blk % yn(i,j,k  ),            &
                            blk % yn(i,j,k+1),            &
                            blk % k_rough(i,j,k  ),       &
                            blk % k_rough(i,j,k+1),       &
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
                            Mc(i,j,k  ) % c,              &
                            Mc(i,j,k+1) % c,              &
                            blk % R(:,i,j,k  ),           &
                            blk % R(:,i,j,k+1),           &
                            3, 1, 2,                &
                            Sc, Sct, Prt, Prl, soot_enabled, inverse_metric )
      enddo; enddo; enddo
    endif

  end subroutine Fluxes_blk

end module MOSE_Mod_Fluxes