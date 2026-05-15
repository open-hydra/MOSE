!>@brief: Module to evaluate correction terms (depending on DSij/Dt) in -RC turbulence corrections.
module MOSE_Lib_SpalartShur
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

contains
  
  subroutine Compute_Velocity_Gradient ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k, n
    real(R8), dimension(3,3) :: vel_diff

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        ! Cell-centered velocity gradient
        vel_diff(:,1) = ( domain%blk(b)%P(nu:nw,i+1,j,k) - domain%blk(b)%P(nu:nw,i-1,j,k) )/2d0
        vel_diff(:,2) = ( domain%blk(b)%P(nu:nw,i,j+1,k) - domain%blk(b)%P(nu:nw,i,j-1,k) )/2d0
        vel_diff(:,3) = ( domain%blk(b)%P(nu:nw,i,j,k+1) - domain%blk(b)%P(nu:nw,i,j,k-1) )/2d0
        domain%blk(b)%vel_gradient(i,j,k)%c = matmul ( vel_diff + 1d-40, domain%blk(b)%m(i,j,k)%c )
      end do ; end do ; end do

      ! Extrapolation in ghost cells -as second derivatives must be computed
      ! i faces
      !$omp do collapse(2) private(n)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
        n = domain % blk(b) % dim(1) + 1
        domain % blk(b) % vel_gradient(0,j,k) % c = domain % blk(b) % vel_gradient(1,j,k) % c
        domain % blk(b) % vel_gradient(n,j,k) % c = domain % blk(b) % vel_gradient(n-1,j,k) % c
      end do ; end do

      ! j faces
      !$omp do collapse(2) private(n)
      do k = 1, domain % blk(b) % dim(3)
      do i = 1, domain % blk(b) % dim(1)   
        n = domain % blk(b) % dim(2) + 1
        domain % blk(b) % vel_gradient(i,0,k) % c = domain % blk(b) % vel_gradient(i,1,k) % c
        domain % blk(b) % vel_gradient(i,n,k) % c = domain % blk(b) % vel_gradient(i,n-1,k) % c
      end do ; end do

      ! k faces
      !$omp do collapse(2) private(n)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)       
        n = domain % blk(b) % dim(3) + 1
        domain % blk(b) % vel_gradient(i,j,0) % c = domain % blk(b) % vel_gradient(i,j,1) % c
        domain % blk(b) % vel_gradient(i,j,n) % c = domain % blk(b) % vel_gradient(i,j,n-1) % c
      end do ; end do

    end do

  end subroutine Compute_Velocity_Gradient


  ! Computation of rstar and rtilde (missing denominator) for the Spalart-Shur correction.
  ! fr1 (the overall correction term) is not computed here as the model constants and rtilde
  ! denominator definition are different between SA and SST.
  subroutine Compute_RC_Terms ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use MOSE_Lib_Fluid, only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor
    use MOSE_Lib_RotatingFrame, only: obj_rot
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k, ii, jj, kk
    real(R8), dimension(3,3,3) :: strain_gradient
    real(R8), dimension(3,3) :: Sij, Wij, strain_transport
    real(R8) :: omega(3), S, O, rtilde

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3) private (i, j, k, ii, jj, kk, strain_gradient), &
      !$omp private (Sij, Wij, strain_transport, omega, S, O, rtilde)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        
        ! Finite difference strain rate tensor
        strain_gradient(:,:,1) = ( Strain_Tensor(domain%blk(b)%vel_gradient(i+1,j,k)%c) & 
                                 - Strain_Tensor(domain%blk(b)%vel_gradient(i-1,j,k)%c) )*0.5d0
        strain_gradient(:,:,2) = ( Strain_Tensor(domain%blk(b)%vel_gradient(i,j+1,k)%c) & 
                                 - Strain_Tensor(domain%blk(b)%vel_gradient(i,j-1,k)%c) )*0.5d0
        strain_gradient(:,:,3) = ( Strain_Tensor(domain%blk(b)%vel_gradient(i,j,k+1)%c) & 
                                 - Strain_Tensor(domain%blk(b)%vel_gradient(i,j,k-1)%c) )*0.5d0

        ! Transformation component by component
        do ii = 1, 3
        do jj = 1, 3
          strain_gradient(ii,jj,:) = matmul ( strain_gradient(ii,jj,:), domain%blk(b)%m(i,j,k)%c )
        enddo
        enddo

        Sij = Strain_Tensor( domain%blk(b)%vel_gradient(i,j,k)%c )
        omega = Vorticity_Vector(domain%blk(b)%vel_gradient(i,j,k)%c)
        ! In a rotating reference frame the Spalart-Shur correction requires the
        ! absolute vorticity: omega_abs = omega_rel + 2*Omega_frame.
        ! Without this the relative vorticity equals S for uniform shear and r*=1
        ! (no correction) regardless of the frame rotation rate.
        if (obj_rot%enabled) omega = omega + 2.0_R8 * obj_rot%Omega_vec
        Wij = Vorticity_Tensor( omega )
        S = sqrt ( 2d0*sum (Sij**2) )
        O = sqrt ( sum(omega**2) )
        domain%blk(b)%rc_term1(i,j,k) = S/O
        ! Material derivative of strain rate tensor, neglecting the time term
        do ii = 1, 3
        do jj = 1, 3
          strain_transport(ii,jj) = dot_product ( domain%blk(b)%p(nu:nw,i,j,k), &
                                                  strain_gradient(ii,jj,:) )
        enddo
        enddo
        rtilde = 0d0
        do ii = 1, 3
        do jj = 1, 3
        do kk = 1, 3
          rtilde = rtilde + Wij(ii,kk)*Sij(jj,kk)*strain_transport(ii,jj)
        enddo
        enddo
        enddo
        ! Note ! Denominator is missing as its different between models
        domain%blk(b)%rc_term2(i,j,k) = 2d0*rtilde 
        ! Note !
      enddo
      enddo
      enddo
    enddo

  end subroutine Compute_RC_Terms

end module MOSE_Lib_SpalartShur