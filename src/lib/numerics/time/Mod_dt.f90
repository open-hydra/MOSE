module MOSE_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_dt, Set_Global_dt

contains

  subroutine Compute_dt ( domain, cfl, vnn, rampa_iter )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use MOSE_Lib_RANS
    use MOSE_Mod_MPI, only: is_local_block, mpi_allreduce_min_r8
    use MOSE_Config_Types_m, only: obj_time_scheme
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: cfl, vnn
    integer, intent(in)  :: rampa_iter
    ! Local
    integer  :: i, j, k, b
    real(R8) :: dtcell, dtglobal, dtglobal_mpi
    logical  :: summed, additive

    ! Read every call, so dt-method follows the runtime re-read of input.ini
    summed   = ( trim(obj_time_scheme%dt_method) == 'summed' )
    additive = ( trim(obj_time_scheme%dt_method) == 'additive' )

    dtglobal = domain % dtglobal

    do b = 1, domain % nb ! loop over blocks
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3) private ( dtcell ), reduction ( min : dtglobal )
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        
        ! Compute local cell dt according to CFL and VNN numbers
        if ( summed .or. additive ) then
          call compute_summed ( rhoi = domain%blk(b)%p(1:nsc,i,j,k), &
                                Ur = domain%blk(b)%Ur(i,j,k), &
                                vel = domain%blk(b)%p(nu:nw,i,j,k), &
                                p = domain%blk(b)%p(np,i,j,k), &
                                rans_ = domain%blk(b)%P(nt:nprim,i,j,k), &
                                met = domain%blk(b)%m(i,j,k)%c, &
                                dtmin = dtcell, &
                                cfl = cfl, &
                                vnn = vnn, &
                                additive = additive )
        else
          call compute ( rhoi = domain%blk(b)%p(1:nsc,i,j,k), &
                         Ur = domain%blk(b)%Ur(i,j,k), &
                         vel = domain%blk(b)%p(nu:nw,i,j,k), &
                         p = domain%blk(b)%p(np,i,j,k), &
                         rans_ = domain%blk(b)%P(nt:nprim,i,j,k), &
                         met = domain%blk(b)%m(i,j,k)%c, &
                         dl = domain%blk(b)%dl(i,j,k)%c, &
                         dtmin = dtcell, &
                         cfl = cfl, &
                         vnn = vnn )
        endif

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
      
      subroutine compute ( rhoi, Ur, vel, p, rans_, met, dl, dtmin, cfl, vnn )
        use MOSE_Global_m
        use MOSE_Config_Types_m, only: obj_prec
        use FLINT_Lib_Thermodynamic
        use MOSE_Lib_Preconditioning, only: comp_Ur
        implicit none
        real(R8), intent(in)  :: rhoi(nsc), vel(3), p, rans_(:)
        real(R8), intent(in)  :: met(3,3), dl(3), cfl, vnn, Ur
        real(R8), intent(out) :: dtmin
        ! Local
        integer :: d
        real(R8) :: rho, Rgas, Sound, dt, versor(3), lambda, mie, mil, mi
        real(R8) :: Alpha, Beta, dx
        real(R8) :: c_star
        real(R8) :: vel_, dummy(3,3)=0d0

        call co_rotot_Rtot ( rhoi, rho, Rgas )
        Sound = f_ss ( rhoi, p, rho, Rgas )
       
        if ( obj_prec%enabled ) then
          Beta = 1.0/Sound**2
          Alpha = 0.5 * ( 1d0 - Beta * Ur**2 )
        endif

        dtmin = 1d8

        do d = 1, ndir

          ! CFL condition along d-direction
          versor = met(d,:) / norm2 ( met(d,:) )
          vel_ = abs ( dot_product ( vel, versor ) )

          if ( obj_prec%enabled ) then
            c_star = Sqrt ( Alpha**2 * Vel_**2 + Ur**2 )
            Vel_ = Vel_ * ( 1d0 - Alpha )
            lambda = vel_ + c_star
          else
            lambda = vel_ + Sound
          endif

          dt = dl(d) / lambda * cfl
          dtmin = min ( dt, dtmin )
        
        enddo

        if ( model==0 ) return
        
        mie = 0.d0
        mil = f_laminarViscosity ( rhoi, p, rho, Rgas )
        if (model==2) then 
          ! Note: approximate mit for 2 equation models. 
          ! velocity gradient assumed 0 and small distance from wall
          call Eddy_Viscosity ( mut=mie, rans_variables=rans_, mul=mil, rho=rho, &
                                vel_gradient=dummy, walldist=1d-6, k_rough=0d0 )
        endif
        mi = mie + mil

        do d = 1, ndir
          ! VNN condition along d-direction
          dt = ( rho * dl(d)**2 * vnn ) / mi
          dtmin = min ( dt, dtmin )
        enddo

      end subroutine compute

      !> dt-method = summed / additive: multi-dimensional limits. With g_d = grad(xi_d) (row d of the
      !> cell metric, |g_d| ~ face area / volume), the convective and diffusive spectral radii are
      !>   Lc = sum_d ( |u.g_d| + c |g_d| )   (flow equations; the turbulence equations: Lu = sum_d |u.g_d|)
      !>   Ld = nu * sum_d |g_d|^2
      !> nu_flow is the largest kinematic diffusivity of the flow equations, momentum 4/3 (mul+mut)/rho
      !> or energy gamma (kl/cp + mut/Prt)/rho; nu_turb is the turbulence model's own, which for SA is
      !> (mul + rho*nit)/(sigma*rho) and exceeds the eddy viscosity of the momentum equations.
      !>   summed:    dt = min ( cfl / Lc , vnn / ( max(nu_flow, nu_turb) sum_d |g_d|^2 ) )
      !>   additive:  1/dt = max over equations of ( Lc/cfl + Ld/vnn ), i.e. the convective and diffusive
      !>              fractions of their limits add up to 1 in every cell (the explicit stability limit
      !>              applies to their sum, which min() lets reach cfl/1 + vnn/0.5 where both bind).
      subroutine compute_summed ( rhoi, Ur, vel, p, rans_, met, dtmin, cfl, vnn, additive )
        use MOSE_Global_m
        use MOSE_Config_Types_m, only: obj_prec, obj_rans
        use FLINT_Lib_Thermodynamic
        implicit none
        real(R8), intent(in)  :: rhoi(nsc), vel(3), p, rans_(:)
        real(R8), intent(in)  :: met(3,3), cfl, vnn, Ur
        logical,  intent(in)  :: additive
        real(R8), intent(out) :: dtmin
        ! Local
        real(R8), parameter :: sigma_SA = 2d0/3d0     ! Spalart-Allmaras sigma, as in Lib_Spalart.f90
        integer  :: d
        real(R8) :: rho, Rgas, Sound, T, gam, cp, versor(3), gnorm, vel_, lambda, lam_sum, lamu_sum, g2_sum
        real(R8) :: Alpha, Beta, c_star, mil, kl, mie, nu_eff, nu_flow, nu_turb, rate, dummy(3,3)=0d0

        call co_rotot_Rtot ( rhoi, rho, Rgas )
        Sound = f_ss ( rhoi, p, rho, Rgas )
        if ( obj_prec%enabled ) then
          Beta = 1.0/Sound**2
          Alpha = 0.5 * ( 1d0 - Beta * Ur**2 )
        endif
        lam_sum  = 0d0
        lamu_sum = 0d0
        g2_sum   = 0d0
        do d = 1, ndir
          gnorm = norm2 ( met(d,:) )
          versor = met(d,:) / gnorm
          vel_ = abs ( dot_product ( vel, versor ) )
          if ( obj_prec%enabled ) then
            c_star = Sqrt ( Alpha**2 * vel_**2 + Ur**2 )
            lambda = vel_ * ( 1d0 - Alpha ) + c_star
          else
            lambda = vel_ + Sound
          endif
          lam_sum  = lam_sum + lambda * gnorm
          lamu_sum = lamu_sum + vel_ * gnorm
          g2_sum   = g2_sum + gnorm**2
        enddo
        dtmin = cfl / lam_sum
        if ( model==0 ) return

        T = p / ( rho * Rgas )
        call co_k_mi_lam_Wilke ( rhoi, rho, T, mil, kl )
        gam = f_gamma ( rhoi, p, rho, Rgas )
        cp  = gam * Rgas / ( gam - 1d0 )
        mie = 0.d0
        if ( model==2 ) &
          call Eddy_Viscosity ( mut=mie, rans_variables=rans_, mul=mil, rho=rho, &
                                vel_gradient=dummy, walldist=1d-6, k_rough=0d0 )
        nu_flow = max ( 4d0/3d0 * ( mil + mie ) / rho, gam * ( kl/cp + mie/obj_rans%Prt ) / rho )
        nu_turb = 0d0
        if ( nrans == 1 ) then
          nu_turb = ( mil + max ( rans_(1), 0d0 ) ) / ( sigma_SA * rho )
        elseif ( nrans > 1 ) then
          nu_turb = ( mil + mie ) / rho
        endif
        if ( additive ) then
          rate = lam_sum / cfl + nu_flow * g2_sum / vnn
          if ( nrans > 0 ) rate = max ( rate, lamu_sum / cfl + nu_turb * g2_sum / vnn )
          dtmin = 1d0 / rate
        else
          nu_eff = max ( nu_flow, nu_turb )
          dtmin = min ( dtmin, vnn / ( nu_eff * g2_sum ) )
        endif
      end subroutine compute_summed

  end subroutine Compute_dt


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