module MOSE_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_dt, Set_Global_dt, Turbulence_Stiffness

contains

  !> Local time step of every cell (blk%dtlocal) and its minimum over the domain (domain%dtglobal, on input an upper
  !> bound). In time-accurate RANS runs the turbulence sources are explicit and the step includes their stiffness:
  !> Compute_dt first calls Turbulence_Stiffness, so the stiffness is that of the current state. Call it outside OpenMP
  !> parallel regions.
  subroutine Compute_dt ( domain, cfl, vnn, rampa_iter )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use MOSE_Lib_RANS
    use MOSE_Config_Types_m, only: obj_time_scheme
    use MOSE_Mod_MPI, only: is_local_block, mpi_allreduce_min_r8
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: cfl, vnn
    integer, intent(in)  :: rampa_iter
    ! Local
    integer  :: i, j, k, b
    real(R8) :: dtcell, dtglobal, dtglobal_mpi, dstiff
    logical  :: explicit_sources

    dtglobal = domain % dtglobal

    ! Time-accurate RANS runs: the turbulence sources are explicit (the point-implicit factor serves steady runs
    ! only), and their stiffness, evaluated here on the current state, enters the step
    explicit_sources = obj_time_scheme%time_accurate .and. model == 2
    if ( explicit_sources ) call Turbulence_Stiffness ( domain )

    do b = 1, domain % nb ! loop over blocks
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3) private ( dtcell, dstiff ), reduction ( min : dtglobal )
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        
        dstiff = 0d0
        if ( explicit_sources ) dstiff = domain%blk(b)%stiff(i,j,k)

        ! Compute local cell dt according to CFL and VNN numbers
        call compute ( rhoi = domain%blk(b)%p(1:nsc,i,j,k), &
                       Ur = domain%blk(b)%Ur(i,j,k), &
                       vel = domain%blk(b)%p(nu:nw,i,j,k), &
                       p = domain%blk(b)%p(np,i,j,k), &
                       rans_ = domain%blk(b)%P(nt:nprim,i,j,k), &
                       met = domain%blk(b)%m(i,j,k)%c, &
                       stiff = dstiff, &
                       dtmin = dtcell, &
                       cfl = cfl, &
                       vnn = vnn )

        ! Apply CFL reduction if required
        if ( domain%iter < rampa_iter ) dtcell = dtcell * domain%iter / rampa_iter

        ! Update local cell dt and global minimum dt
        domain % blk(b) % dtlocal(i,j,k) = dtcell
        dtglobal = min ( dtcell, dtglobal )

      enddo; enddo; enddo
      !$omp end parallel
    enddo ! end of loop over blocks

    ! MPI: global minimum across all ranks
    call mpi_allreduce_min_r8(dtglobal, dtglobal_mpi)
    domain % dtglobal = dtglobal_mpi

    contains
      
      ! Local time step of a cell, one additive rule for convection and diffusion (Hundsdorfer 2000, notes on the
      ! numerical solution of advection-diffusion-reaction equations, sec. 9.3 and 6; Arnone, Liou, Povinelli 1993,
      ! NASA TM-106373, eq. 26-28):
      !   1/dt = sum_d lambda_d / cfl  +  D_max sum_d |grad xi_d|^2 / vnn,   lambda_d = |u . grad xi_d| + c |grad xi_d|,
      ! with grad xi_d the row d of the cell metric (the normal of the faces in direction d over the cell width) and D_max
      ! the largest diffusivity of the equations. Without dimensional splitting the convective limits of the directions
      ! add up (an explicit scheme on square 2-D cells needs the sum, not the largest of the two); the convective and the
      ! diffusive rates add up as well. For the linearised MUSCL-upwind (Fromm's scheme: the van Leer slope, and the van
      ! Albada and MC ones, in smooth flow; not minmod and superbee) and central-diffusion operators the rule is stable for
      ! cfl <= 1 and vnn <= 1/2 with Heun (RK2), cfl <= 1.17 and vnn <= 0.62 with RK3 (von Neumann limits, proven; at the
      ! worst cell Peclet number the additive law is exact). In 3-D, on a cell skewed in all three pairs of directions, the
      ! viscous rate takes the factor (2 + m + 1/m)/4 below, which bounds the excess of the spectral radius of the viscous
      ! operator over 4 D_max sum_d |grad xi_d|^2 there.
      ! With explicit turbulence sources (time-accurate RANS runs) the stiffness d of their destruction adds to the
      ! rates with the constant 1: 1/dt += d, so that dt d <= 1. A forward-Euler stage then keeps nu~, k and omega
      ! positive, and so do Heun and RK3, convex combinations of forward-Euler stages (Hundsdorfer 2000, sec. 6.2,
      ! theorem 6.4 and the explicit trapezoidal rule; C.-W. Shu, S. Osher, J. Comput. Phys. 77, 1988).
      subroutine compute ( rhoi, Ur, vel, p, rans_, met, stiff, dtmin, cfl, vnn )
        use MOSE_Global_m
        use MOSE_Config_Types_m, only: obj_prec, obj_sim_param, obj_rans
        use FLINT_Lib_Thermodynamic
        use MOSE_Lib_Preconditioning, only: comp_Ur
        implicit none
        real(R8), intent(in)  :: rhoi(nsc), vel(3), p, rans_(:)
        real(R8), intent(in)  :: met(3,3), cfl, vnn, Ur
        real(R8), intent(in)  :: stiff    ! stiffness of the explicit turbulence sources, 0 without them [1/s]
        real(R8), intent(out) :: dtmin
        ! Local
        real(R8), parameter :: sigma_SA = 2d0/3d0   ! diffusion constant of Spalart-Allmaras (Lib_Spalart.f90)
        integer  :: d, T_i, Tint(2)
        real(R8) :: rho, Rgas, Sound, Alpha, Beta, c_star, gn, un, lambda
        real(R8) :: rate_c, rate_v, temperature, Tdiff, cp, cv, mil, kl, mie, Dmax
        real(R8) :: Dm(nsc), dummy(3,3)
        real(R8) :: cos12, cos13, cos23, ssum, qprod, mroot

        call co_rotot_Rtot ( rhoi, rho, Rgas )
        Sound = f_ss ( rhoi, p, rho, Rgas )

        if ( obj_prec%enabled ) then
          Beta = 1.0/Sound**2
          Alpha = 0.5 * ( 1d0 - Beta * Ur**2 )
        endif

        ! Convective spectral radii and squared metrics, summed over the directions
        rate_c = 0d0
        rate_v = 0d0
        do d = 1, ndir
          gn = norm2 ( met(d,:) )                          ! |grad xi_d|: 1 / width of the cell across its d faces
          un = abs ( dot_product ( vel, met(d,:) ) ) / gn  ! velocity normal to the d faces

          if ( obj_prec%enabled ) then
            c_star = Sqrt ( Alpha**2 * un**2 + Ur**2 )
            lambda = un * ( 1d0 - Alpha ) + c_star
          else
            lambda = un + Sound
          endif

          rate_c = rate_c + lambda * gn
          rate_v = rate_v + gn**2
        enddo

        if ( model==0 ) then
          dtmin = cfl / rate_c
          return
        endif

        ! 3-D: on a cell whose three directions are skewed against each other, so that the cosines cos12, cos13, cos23
        ! between the rows grad xi_d of the metric give cos12^2 + cos13^2 + cos23^2 + 2 cos12 cos13 cos23 > 1, the
        ! spectral radius of the viscous operator exceeds 4 D_max sum_d |grad xi_d|^2 by a factor below 9/8. The rate
        ! takes the bound (2 + m + 1/m)/4 of that factor, m + 1 the largest eigenvalue of the matrix of the cosines (m the
        ! largest root of m^3 - s m - 2 cos12 cos13 cos23 = 0, s the sum of the squares of the cosines): exact on
        ! rhombohedra, and not taken where the rate already holds (orthogonal cells, and cells with one direction
        ! orthogonal to the other two)
        if ( ndir == 3 ) then
          cos12 = dot_product ( met(1,:), met(2,:) ) / ( norm2 ( met(1,:) ) * norm2 ( met(2,:) ) )
          cos13 = dot_product ( met(1,:), met(3,:) ) / ( norm2 ( met(1,:) ) * norm2 ( met(3,:) ) )
          cos23 = dot_product ( met(2,:), met(3,:) ) / ( norm2 ( met(2,:) ) * norm2 ( met(3,:) ) )
          ssum  = cos12**2 + cos13**2 + cos23**2
          qprod = cos12 * cos13 * cos23
          if ( ssum + 2d0 * qprod > 1d0 ) then
            mroot  = 2d0 * sqrt ( ssum / 3d0 ) * &
                     cos ( acos ( max ( -1d0, min ( 1d0, 3d0 * sqrt ( 3d0 ) * qprod / ssum**1.5d0 ) ) ) / 3d0 )
            rate_v = rate_v * 0.25d0 * ( 2d0 + mroot + 1d0 / mroot )
          end if
        end if

        ! Largest diffusivity: momentum 4/3 (mu + mu_t)/rho; energy, with the conductivity acting on rho cv T,
        ! (k + mu_t cp/Pr_t)/(rho cv) (gamma nu/Pr for a perfect gas); species; the turbulence variables. The k and omega
        ! diffusivities of SST and Wilcox 2006, (mu + sigma mu_t)/rho with sigma <= 1 and mu_t = rho k/omega at zero
        ! strain, stay below the momentum one.
        temperature = p / ( rho * Rgas )
        T_i     = int ( temperature )
        Tdiff   = temperature - T_i
        Tint(1) = T_i
        Tint(2) = T_i + 1
        cp = f_cp_expr ( rhoi, Tint, Tdiff, rho )
        cv = cp - Rgas
        call co_k_mi_lam_Wilke_expr ( rhoi, rho, Tint, Tdiff, mil, kl )
        if ( obj_sim_param%Prl > 0d0 ) kl = mil * cp / obj_sim_param%Prl

        mie = 0.d0
        if (model==2) then
          ! Note: approximate mit for 2 equation models.
          ! velocity gradient assumed 0 and small distance from wall
          dummy = 0d0
          call Eddy_Viscosity ( mut=mie, rans_variables=rans_, mul=mil, rho=rho, &
                                vel_gradient=dummy, walldist=1d-6, k_rough=0d0 )
        endif

        Dmax = max ( 4d0/3d0 * ( mil + mie ) / rho, ( kl + mie * cp / obj_rans%Prt ) / ( rho * cv ) )

        if ( nsc > 1 ) then
          if ( obj_sim_param%Sc > 0d0 ) then
            Dm = ( mil / obj_sim_param%Sc + mie / obj_rans%Sct ) / rho
          else
            call co_DS_expr ( rhoi, rho, Tint, Tdiff, p, Dm )
            Dm = Dm + mie / ( rho * obj_rans%Sct )
          endif
          Dmax = max ( Dmax, maxval ( Dm ) )
        endif

        if ( model==2 .and. nRANS==1 ) Dmax = max ( Dmax, ( mil + rans_(1) ) / ( sigma_SA * rho ) )

        dtmin = 1d0 / ( rate_c / cfl + Dmax * rate_v / vnn )

        if ( stiff > 0d0 ) dtmin = 1d0 / ( 1d0 / dtmin + stiff )

      end subroutine compute

  end subroutine Compute_dt


  !> Stiffness of the turbulence sources of every cell on the current state of the domain: blk%stiff = d, the
  !> derivative of the destruction (Spalart-Allmaras 2 c_w1 f_w nu~/y^2; SST and Wilcox 2006 max(beta* omega,
  !> 2 beta omega)), written by the source routines. It fills the ghost cells and evaluates the turbulence sources as
  !> the first stage of a step does, without the fluxes, so it is collective (MPI ghost exchange) and must be called
  !> outside OpenMP parallel regions, by all processes. It leaves in blk%r a residual that nobody may use: the next stage
  !> starts by zeroing the residuals, and its own ghost fill replaces the ghost values set here. Compute_dt calls it in
  !> time-accurate RANS runs; call it alone only when the stiffness is needed without the step.
  subroutine Turbulence_Stiffness ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Lib_Ghost, only: Fill_Ghost_Cell, Ghost_Wall_Extrapolation
    use MOSE_Lib_RANS,  only: RANS_Source_Terms
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain

    !$omp parallel
    call Fill_Ghost_Cell ( domain )
    call Ghost_Wall_Extrapolation ( domain )
    call RANS_Source_Terms ( domain )
    !$omp end parallel

  end subroutine Turbulence_Stiffness


  subroutine Set_Global_dt ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        domain % blk(b) % dtlocal(i,j,k) = domain % dtglobal
      enddo; enddo; enddo
      !$omp end parallel
    enddo
    
  end subroutine Set_Global_dt


  subroutine Set_Dt_Control ( domain, control )
    use MOSE_Advanced_Types_m
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    real(R8), intent(in) :: control
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k


    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      !$omp parallel
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        domain % blk(b) % dtlocal(i,j,k) = &
        min( domain % blk(b) % dtlocal(i,j,k), control*domain % dtglobal )
      end do
      end do
      end do
      !$omp end parallel

    end do

  end subroutine Set_Dt_Control

end module MOSE_Mod_dt