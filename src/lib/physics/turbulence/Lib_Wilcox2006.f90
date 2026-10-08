!> @brief Wilcox (2006) k-omega model, as specified by the NASA Turbulence Modeling Resource (TMR):
!> https://tmbwg.github.io/turbmodels/wilcox.html (formerly turbmodels.larc.nasa.gov/wilcox.html).
!> Exact production P = tau_ij du_i/dx_j with the k-equation limiter 20 beta* rho omega k
!> (TMR Wilcox2006-klim). Without k-energy-coupling, -2/3 rho k delta_ij is left out of the stress of
!> the momentum and energy equations (TMR Wilcox2006-klim-e); with it, k is in the total energy
!> (TMR Wilcox2006-klim; see Lib_SST.f90 and the TMR compressible-implementation page).
!>
!> References
!>   Wilcox, D. C., "Formulation of the k-omega Turbulence Model Revisited", AIAA J. 46(11), 2008,
!>     pp. 2823-2838, doi:10.2514/1.36541.
!>   Wilcox, D. C., Turbulence Modeling for CFD, 3rd ed., DCW Industries, 2006.
module MOSE_Lib_Wilcox2006
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

  real(R8), parameter, private :: &
  sigma_k = 0.6d0, &
  sigma_w = 0.5d0, &
  beta_star = 0.09d0, &
  beta_0 = 0.0708d0, &
  gamma = 13d0/25d0, &
  Clim = 7d0/8d0

  !> C of the asymptotic omega wall condition, omega_wall = C nu/d**2: the limit 6/beta_0 of the omega equation as y -> 0
  !> (NASA Turbulence Modeling Resource, Wilcox 2006 page)
  real(R8), parameter, public :: Wilcox2006_omega_wall_asymptotic = 6d0/beta_0

contains

  subroutine Wilcox2006_Source_Terms ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_rans, obj_time_scheme
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Base_Types_m, only: MOSE_tensor_3D_type
    use MOSE_Lib_Stencil_Metrics, only: Source_Metrics
    implicit none
    type(MOSE_domain_type), intent(inout), target :: domain
    ! Local
    integer :: b
    type(MOSE_tensor_3D_type), pointer :: Mi(:,:,:), Mt(:,:,:)   ! cell metrics of the sources (Source_Metrics)
    logical :: k_energy_coupling, point_implicit

    k_energy_coupling = obj_rans%k_energy_coupling
    ! The point-implicit factor serves steady runs only: in time-accurate runs the sources are explicit and their
    ! stiffness enters the time step (Compute_dt)
    point_implicit = obj_rans%point_implicit .and. .not. obj_time_scheme%time_accurate

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      call Source_Metrics ( domain % blk(b), Mi, Mt=Mt )
      call Wilcox_Blk ( domain % blk(b) % P,       &
                        domain % blk(b) % r,       &
                        Mi, Mt,                   &
                        domain % blk(b) % vol,     &
                        domain % blk(b) % dtlocal, &
                        domain % blk(b) % stiff,   &
                        domain % blk(b) % dim,     &
                        k_energy_coupling,         &
                        point_implicit )

    enddo

  end subroutine Wilcox2006_Source_Terms


  subroutine Wilcox_Blk ( Prim, Res, M, Mt, Volume, dt, stiff, n, k_energy_coupling, point_implicit )
    use MOSE_Base_Types_m
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_Fluid
    use MOSE_Lib_RotatingFrame, only: obj_rot
    implicit none
    integer, intent(in) :: n(3)
    logical, intent(in) :: k_energy_coupling, point_implicit
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Volume
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Mt   ! turbulence variables
    real(R8), dimension(n(1), n(2), n(3)), intent(in) :: dt
    real(R8), dimension(n(1), n(2), n(3)), intent(out) :: stiff  ! largest destruction rate: max(beta* omega, 2 beta omega)
    ! Local
    integer :: i, j, k, ii, jj, kk
    real(R8) :: rho, kap, ome, Gradvel(3,3), Divel, Sij(3,3), Sij_bar(3,3), diag, ome_hat, mi_t, Tij(3,3), Prod(2)
    real(R8) :: Gradkw(2,3), dkDotdw, sigma_d, Diff, Sij_hat(3,3), Wij(3,3), Num, X_w, f_beta, beta, Diss(2), Source(2)
    real(R8) :: omega_w(3), fk, fw

    !$omp do collapse (3) private (rho, kap, ome, Gradvel, Divel, Sij, Sij_bar, diag, ome_hat, mi_t, Tij, Prod), &
    !$omp private ( Gradkw, dkDotdw, sigma_d, Diff, Sij_hat, Wij, Num, X_w, f_beta, beta, Diss, Source, fk, fw ), &
    !$omp private ( i, j, k, ii, jj, kk, omega_w )

    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      ! Local variables
      rho = Sum ( Prim(1:nsc,i,j,k) )
      kap = Prim(nt,i,j,k) / rho
      ome = Prim(nt+1,i,j,k) / rho

      ! Velocity gradient computation
      Gradvel(:,1) = ( Prim(nu:nw,i+1,j,k) - Prim(nu:nw,i-1,j,k) ) * 0.5d0
      Gradvel(:,2) = ( Prim(nu:nw,i,j+1,k) - Prim(nu:nw,i,j-1,k) ) * 0.5d0
      Gradvel(:,3) = ( Prim(nu:nw,i,j,k+1) - Prim(nu:nw,i,j,k-1) ) * 0.5d0
      Gradvel = Matmul ( Gradvel, M(i,j,k) % c )

      ! Strain tensor and similar stuff
      Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3) 
      Sij = Strain_Tensor ( Gradvel )
      Sij_bar = Sij
      diag = 1d0/3d0 * Divel
      Sij_bar(1,1) = Sij_bar(1,1) - diag
      Sij_bar(2,2) = Sij_bar(2,2) - diag
      Sij_bar(3,3) = Sij_bar(3,3) - diag

      ! Eddy viscosity
      ome_hat = Max ( ome, Clim * Sqrt( 2d0 * Sum ( Sij_bar**2 ) / beta_star ) )
      mi_t = rho * kap / ome_hat

      ! Turbulent stress tensor
      Tij = 2d0 * mi_t * Sij
      diag = 2d0/3d0 * mi_t * Divel + 2d0/3d0 * rho * kap
      Tij(1,1) = Tij(1,1) - diag
      Tij(2,2) = Tij(2,2) - diag
      Tij(3,3) = Tij(3,3) - diag

      ! Production term
      Prod(1) = 0d0
      do ii = 1, 3
      do jj = 1, 3
        Prod(1) = Prod(1) + Tij(ii,jj) * Gradvel(ii,jj)
      end do
      end do
      Prod(2) = ( gamma * ome / kap ) * Prod(1)
      Prod(1) = Min( Prod(1), 20d0*beta_star*rho*ome*kap ) ! k-eqn production limiter

      ! Gradient of kap and ome
      Gradkw(:,1) = ( Prim(nt:nt+1,i+1,j,k) / sum( Prim(1:nsc,i+1,j,k) ) - &
                      Prim(nt:nt+1,i-1,j,k) / sum( Prim(1:nsc,i-1,j,k) ) ) * 0.5d0
      Gradkw(:,2) = ( Prim(nt:nt+1,i,j+1,k) / sum( Prim(1:nsc,i,j+1,k) ) - &
                      Prim(nt:nt+1,i,j-1,k) / sum( Prim(1:nsc,i,j-1,k) ) ) * 0.5d0
      Gradkw(:,3) = ( Prim(nt:nt+1,i,j,k+1) / sum( Prim(1:nsc,i,j,k+1) ) - &
                      Prim(nt:nt+1,i,j,k-1) / sum( Prim(1:nsc,i,j,k-1) ) ) * 0.5d0
      Gradkw = Matmul ( Gradkw, Mt(i,j,k) % c )

      ! Cross Diffusion term
      dkDotdw = Dot_Product ( Gradkw(1,:), Gradkw(2,:) )
      sigma_d = 1d0/8d0
      if ( dkDotdw <= 0d0 ) sigma_d = 0d0
      Diff = sigma_d * rho / ome * dkDotdw

      ! Dissipation term. S^_ki = S_ki - (1/2) du_m/dx_m delta_ki: off the diagonal S^ = S
      Sij_hat = Sij
      Sij_hat(1,1) = Sij(1,1) - 0.5d0 * Divel
      Sij_hat(2,2) = Sij(2,2) - 0.5d0 * Divel
      Sij_hat(3,3) = Sij(3,3) - 0.5d0 * Divel

      omega_w = Vorticity_Vector ( Gradvel )
      ! Use absolute vorticity for f_beta: W_ij^abs = W_ij^rel + Omega_frame tensor.
      if (obj_rot%enabled) omega_w = omega_w + 2.0_R8 * obj_rot%Omega_vec
      Wij = Vorticity_Tensor ( omega_w )

      Num = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        Num = Num + Wij(ii,jj) * Wij(jj,kk) * Sij_hat(kk,ii)
      end do
      end do
      end do

      X_w = Abs ( Num / ( beta_star * ome )**3 )
      f_beta = ( 1d0 + 85d0 * X_w ) / ( 1d0 + 100d0 * X_w )
      beta = beta_0 * f_beta

      Diss(1) = beta_star * rho * kap * ome
      Diss(2) = beta * rho * ome**2

      ! Source terms
      Source(1) = Prod(1) - Diss(1)
      Source(2) = Prod(2) - Diss(2) + Diff

      ! Stiffness of the destruction terms (derivatives with respect to k and omega)
      stiff(i,j,k) = max ( beta_star * ome, 2d0 * beta * ome )

      ! Point-implicit factors of steady runs. Each multiplies the whole increment of its row: the fluxes already in
      ! Res (Internal_Fluxes, BC_Fluxes) here, the source below. The steady state, zero residual, then does not depend
      ! on the step (point-implicit form of Bussing and Murman: H. C. Yee, J. L. Shinn, NASA TM-89415, 1986, eq. 3.1;
      ! T. H. Pulliam, AIAA Paper 93-3360, eq. 8)
      fk = 1d0 ; fw = 1d0
      if ( point_implicit ) then
        fk = 1d0 / ( 1d0 + dt(i,j,k) * beta_star * ome )
        fw = 1d0 / ( 1d0 + dt(i,j,k) * 2d0 * beta * ome )
        Res(nt,  i,j,k) = Res(nt,  i,j,k) * fk
        Res(nt+1,i,j,k) = Res(nt+1,i,j,k) * fw
      end if

      Res(nt,  i,j,k) = Res(nt,  i,j,k) - Source(1) * Volume(i,j,k) * fk
      Res(nt+1,i,j,k) = Res(nt+1,i,j,k) - Source(2) * Volume(i,j,k) * fw
      if (k_energy_coupling) Res(np,i,j,k) = Res(np,i,j,k) + Source(1) * Volume(i,j,k) * fk

    enddo ; enddo ; enddo ! (i, j, k) loop

  end subroutine Wilcox_Blk

  
  subroutine Wilcox2006_Eddy_Viscosity ( mi_t, rkw, mi_l, rho, Gradvel, dist, k_rough )
    use MOSE_Global_m
    use MOSE_Lib_Fluid
    implicit none
    real(R8), intent(in), dimension(nRANS) :: rkw        ! : rkw = [rho*k rho*w]
    real(R8), intent(in)                   :: mi_l       ! : Molecular viscosity
    real(R8), intent(in)                   :: rho        ! : Density
    real(R8), intent(in)                   :: dist       ! : Wall distance (not used)
    real(R8), intent(in)                   :: k_rough    ! : Wall roughness (not modelled)
    real(R8), intent(in), dimension(3,3)   :: Gradvel    ! : Velocity gradient
    real(R8), intent(out)                  :: mi_t       ! : Eddy viscosity
    ! Local
    real(R8) :: w, Divel, Sij(3,3), diag, w_hat

    w = rkw(2) / rho
    Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3) 
    Sij = Strain_Tensor ( Gradvel )
    diag = 1d0/3d0 * Divel
    Sij(1,1) = Sij(1,1) - diag
    Sij(2,2) = Sij(2,2) - diag
    Sij(3,3) = Sij(3,3) - diag
    w_hat = Max ( w, Clim * Sqrt( 2d0 * Sum( Sij**2 ) / beta_star ) )
    mi_t = rkw(1) / w_hat

  end subroutine Wilcox2006_Eddy_Viscosity


  subroutine Wilcox2006_RANS_Diffusive_Flux ( Flux, rkw, mi_l, vGrad, Gradkw, rho, Area, Normal, Dist )
    use MOSE_Global_m
    implicit none
    real(R8), intent(in), dimension(nRANS)   :: rkw                    ! : [ rho*k rho*w ]
    real(R8), intent(in)                     :: mi_l                   ! : laminar viscosity
    real(R8), intent(in), dimension(3,3)     :: vGrad                  ! : Velocity gradient (not used)clear
    real(R8), intent(in), dimension(nRANS,3) :: Gradkw                 ! : Gradient of [ k w ]
    real(R8), intent(in)                     :: rho                    ! : Density
    real(R8), intent(in)                     :: Area, Normal(3), Dist  ! : Metrics
    real(R8), intent(out), dimension(nRANS)  :: Flux                   ! : Diffusive flux
    ! Local
    real(R8) :: rhok, w

    rhok = rkw(1)
    w = rkw(2)/rho
    Flux(1) = ( mi_l + sigma_k * rhok / w ) * Dot_Product ( Gradkw(1,:), Normal )
    Flux(2) = ( mi_l + sigma_w * rhok / w ) * Dot_Product ( Gradkw(2,:), Normal )
    Flux = Flux * Area

  end subroutine Wilcox2006_RANS_Diffusive_Flux


  subroutine Wilcox2006_Set_Wall_Values ( mi_l, rkw_cell, rkw, dist, k_rough )
    use MOSE_Global_m
    use MOSE_Config_Types_m, only: obj_rans
    implicit none
    real(R8), intent(in),  dimension(nRANS)  :: rkw_cell         ! : boundary cell values (not used)
    real(R8), intent(out), dimension(nRANS)  :: rkw              ! : [ rho*k rho*w ] at wall
    real(R8), intent(in)                     :: mi_l             ! : laminar viscosity at wall
    real(R8), intent(in)                     :: dist             ! : cell center distance from wall
    real(R8), intent(in)                     :: k_rough          ! : wall roughness (not modelled)

    rkw(1) = 0d0 ! Solid surface condition on k
    ! Smooth-surface BC: omega_wall = C*mi_l/dist**2, with C set by omega-wall-bc
    ! (800 = Menter practical, 6/beta_0 = asymptotic y->0 limit). See Setup_RANS_Model.
    rkw(2) = obj_rans%omega_wall_coef * mi_l / ( dist**2 )

  end subroutine Wilcox2006_Set_Wall_Values

  subroutine Wilcox2006_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use MOSE_Global_m
    implicit none
    real(R8), intent(in),  dimension(nRANS) :: prim, wall_prim  ! : boundary cell and wall variables
    real(R8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
    real(R8), intent(out), dimension(nRANS) :: ghost_prim ! : ghost cell RANS variables

    ghost_prim(1) = - prim(1) ! enforcing k_wall == 0
    ghost_prim(2) = rho * 2d0 * wall_prim(2)/wall_rho - prim(2) ! extrapolating from omega_wall

  end subroutine Wilcox2006_Extrapolate_Wall


  subroutine Wilcox2006_Enforce_Realizability ( rhokrhow )
    implicit none
    real(R8), intent(inout),  dimension(:) :: rhokrhow

    if ( rhokrhow(1) < 1d-10 ) rhokrhow(1) = 1d-10
    if ( rhokrhow(2) < 1d-10 ) rhokrhow(2) = 1d-10

  end subroutine Wilcox2006_Enforce_Realizability
  
end module MOSE_Lib_Wilcox2006