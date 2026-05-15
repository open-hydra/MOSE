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
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: cfl, vnn
    integer, intent(in)  :: rampa_iter
    ! Local
    integer  :: i, j, k, b
    real(R8) :: dtcell, dtglobal, dtglobal_mpi

    dtglobal = domain % dtglobal

    do b = 1, domain % nb ! loop over blocks
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3) private ( dtcell ), reduction ( min : dtglobal )
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        
        ! Compute local cell dt according to CFL and VNN numbers
        call compute ( rhoi = domain%blk(b)%p(1:nsc,i,j,k), &
                       vel = domain%blk(b)%p(nu:nw,i,j,k), &
                       p = domain%blk(b)%p(np,i,j,k), &
                       rans_ = domain%blk(b)%P(nt:nprim,i,j,k), &
                       met = domain%blk(b)%m(i,j,k)%c, &
                       dl = domain%blk(b)%dl(i,j,k)%c, &
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
      
      subroutine compute ( rhoi, vel, p, rans_, met, dl, dtmin, cfl, vnn )
        use MOSE_Global_m
        use FLINT_Lib_Thermodynamic
        implicit none
        real(R8), intent(in)  :: rhoi(nsc), vel(3), p, rans_(:)
        real(R8), intent(in)  :: met(3,3), dl(3), cfl, vnn
        real(R8), intent(out) :: dtmin
        ! Local
        integer :: d
        real(R8) :: rho, Rgas, Sound, dt, versor(3), lambda, mie, mil, mi
        real(R8) :: vel_, dummy(3,3)=0d0

        call co_rotot_Rtot ( rhoi, rho, Rgas )
        Sound = f_ss ( rhoi, p, rho, Rgas )
       
        dtmin = 1d8
        do d = 1, ndir
          ! CFL condition along d-direction
          versor = met(d,:) / norm2 ( met(d,:) )
          vel_ = abs ( dot_product ( vel, versor ) )

          lambda  = vel_ + Sound
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
                                vel_gradient=dummy, walldist=1d-6 )
        endif
        mi = mie + mil

        do d = 1, ndir
          ! VNN condition along d-direction
          dt = ( rho * dl(d)**2 * vnn ) / mi
          dtmin = min ( dt, dtmin )
        enddo

      end subroutine compute

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

end module MOSE_Mod_dt