!> @brief Module for Spalart-Allmaras constants and subroutines, see https://turbmodels.larc.nasa.gov/spalart.html. Note: no ft2 term.
module MOSE_Lib_Spalart
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

  real(R8), parameter, private :: &
    cb1 = 0.1355d0, &
    cb2 = 0.622d0, &
    cw2 = 0.3d0, &
    cw3 = 2.0d0, &
    cv1 = 7.1d0, &
    cK = 0.41d0, &
    ct1 = 1.d0, &
    ct2 = 2.d0, &
    ct3 = 1.2d0, &
    ct4 = 0.5d0, &
    sigma = 2d0/3d0, &
    c2 = 0.7d0, &
    c3 = 0.9d0, &
    cr1 = 1.d0, &
    cr2 = 12.d0, &
    cr3 = 1.d0, &
    Crot = 2.d0, &
    cw1 = cb1/cK**2 + (1d0+cb2)/sigma, &
    C_cr1 = 0.3d0, &
    tolerance = 1.d-2

contains

  subroutine Spalart_Source_Terms ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_rans
    use MOSE_Lib_Fluid, only : Vorticity_Vector, Vorticity_Tensor
    use MOSE_Lib_SpalartShur, only : Compute_Velocity_Gradient, Compute_RC_Terms
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b
    logical :: SpalartShur, SAR, SAcomp

    SpalartShur = obj_rans%SpalartShur
    SAR         = obj_rans%SAR
    SAcomp      = obj_rans%SAcomp

    if ( SpalartShur ) then
      call Compute_Velocity_Gradient ( domain )
      call Compute_RC_Terms ( domain )
    end if
    
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      call SA_Source_Blk ( domain % blk(b) % P,            &
                           domain % blk(b) % r,            &
                           domain % blk(b) % M,            &
                           domain % blk(b) % yn,           &
                           domain % blk(b) % vol,          &
                           domain % blk(b) % vel_gradient, &
                           domain % blk(b) % rc_term1,     &
                           domain % blk(b) % rc_term2,     &
                           domain % blk(b) % dim,          &
                           SpalartShur, SAR, SAcomp        )
    
    end do

  end subroutine Spalart_Source_Terms


  subroutine SA_Source_Blk ( Prim, Res, M, yn, vol, gradv, rc1, rc2, n, SpalartShur, SAR, SAcomp )
    use MOSE_Base_Types_m
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_Fluid
    use MOSE_Lib_RotatingFrame, only: obj_rot
    implicit none
    integer, intent(in) :: n(3)
    logical, intent(in) :: SpalartShur, SAR, SAcomp
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: yn
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: gradv
    real(R8), dimension(n(1), n(2), n(3)), intent(in) :: rc1, rc2
    ! Local
    integer :: i, j, k
    real(R8) :: rho, Rgas, nit, mil, Gradvel(3,3), Gradnit(3), omega(3), Om, abs_Eij, Wij(3,3)
    real(R8) :: Eij(3,3), chi, Stilde, Sbar, fw, Gradnit2, r_Fun, g_Fun, fr1, rstar, D2, rtilde
    real(R8) :: Production, Diffusion, Destruction, Source

    !$omp do collapse (3) private ( rho, Rgas, nit, mil, Gradvel, Gradnit, omega, Om, abs_Eij, Wij ), &
    !$omp private ( Eij, chi, Stilde, Sbar, fw, Gradnit2, r_Fun, g_Fun, fr1, rstar, D2, rtilde ), &
    !$omp private ( Production, Diffusion, Destruction, Source, i, j, k)
    
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      call Co_rotot_Rtot ( Prim(1:nsc,i,j,k), rho, Rgas )
      mil = f_laminarViscosity ( Prim(1:nsc,i,j,k), Prim(np,i,j,k), rho, Rgas )
      nit = Prim (nt,i,j,k) / rho  ! note: its ni-tilde
      chi = Prim (nt,i,j,k) / mil  ! Spalart X variable

      if ( SpalartShur ) then
        Gradvel = gradv(i,j,k) % c
      else 
        Gradvel(:,1) = ( Prim(nu:nw,i+1,j,k) - Prim(nu:nw,i-1,j,k) ) * 0.5d0
        Gradvel(:,2) = ( Prim(nu:nw,i,j+1,k) - Prim(nu:nw,i,j-1,k) ) * 0.5d0
        Gradvel(:,3) = ( Prim(nu:nw,i,j,k+1) - Prim(nu:nw,i,j,k-1) ) * 0.5d0
        Gradvel = Matmul ( Gradvel, M(i,j,k) % c )
      end if

      omega = Vorticity_Vector ( Gradvel )    ! vorticity vector (relative frame)
      ! Use absolute vorticity for SA: omega_abs = omega_rel + 2*Omega_frame.
      ! This ensures Om, Stilde, and D2 (SA-RC denominator) are physically correct.
      if (obj_rot%enabled) omega = omega + 2.0_R8 * obj_rot%Omega_vec
      Om    = Sqrt ( sum ( omega**2 ) )     ! vorticity module

      ! Production term computation
      Sbar = nit / ( cK**2 * yn(i,j,k)**2 ) * fv2 (chi)
      if ( Sbar >= (-c2*Om) ) then
        Stilde = Om + Sbar
      else
        Stilde = Om + Om * (c2**2*Om + c3*Sbar) / ((c3 - 2*c2)*Om - Sbar)
      endif
      if ( Stilde == 0d0 ) Stilde = 1d-40 ! Just to avoid division by zero
        
      Eij     = Strain_Tensor ( Gradvel ) ! Strain rate tensor
      abs_Eij = sqrt ( 2*sum (Eij**2) )

      if ( SpalartShur ) then

        Wij = Vorticity_Tensor ( omega )
        rstar = rc1(i,j,k)
        D2 = 0.5d0*(abs_Eij**2 + Om**2)
        rtilde = rc2(i,j,k) / D2**2
        fr1 = (1.d0 + cr1)*(2.d0*rstar/(1.d0 + rstar))*(1.d0 - cr3*datan(cr2*rtilde)) - cr1
        Production = fr1 * cb1 * Stilde * nit

      elseif ( SAR ) then
          
        Production = cb1 * ( Stilde + Crot * Min(0d0, abs_Eij - Om) ) * nit

      else

        Production = cb1 * Stilde * nit
        
      end if

      ! finite differences for ni-tilde in computational reference frame
      Gradnit(1) = ( Prim(nt,i+1,j,k) / sum( Prim(1:nsc,i+1,j,k) ) - &
                     Prim(nt,i-1,j,k) / sum( Prim(1:nsc,i-1,j,k) ) ) * 0.5d0
      Gradnit(2) = ( Prim(nt,i,j+1,k) / sum( Prim(1:nsc,i,j+1,k) ) - &
                     Prim(nt,i,j-1,k) / sum( Prim(1:nsc,i,j-1,k) ) ) * 0.5d0
      Gradnit(3) = ( Prim(nt,i,j,k+1) / sum( Prim(1:nsc,i,j,k+1) ) - &
                     Prim(nt,i,j,k-1) / sum( Prim(1:nsc,i,j,k-1) ) ) * 0.5d0

      Gradnit = Matmul ( Gradnit, M(i,j,k) % c ) ! ni-tilde gradient

      Gradnit2 = Sum ( Gradnit**2 )              ! ni-tilde gradient squared

      Diffusion = cb2 / sigma * Gradnit2         ! Non-conservative Diffusion term

      ! Destruction term computation
      r_Fun = nit / ( Stilde * cK**2 * yn(i,j,k)**2 )
      r_Fun = Min ( r_Fun, 10d0 ) 
      g_Fun = r_Fun + cw2 * ( r_Fun**6 - r_Fun )
      fw = g_Fun * ( ( 1d0 + cw3**6 ) / ( g_Fun**6 + cw3**6 ) )**(1.d0/6.d0)
      Destruction = cw1 * fw * ( nit / yn(i,j,k) )**2

      if ( SAcomp ) call Compressibility_Correction &
      ( nit, rho, chi, Prim(1:nsc,i,j,k), Prim(np,i,j,k), Rgas, Om, Production )

      ! Source multiplied by rho since the SA equation is integrated in conservative form
      Source = rho * ( Production + Diffusion - Destruction ) * vol(i,j,k)
      Res(nt,i,j,k) = Res(nt,i,j,k) - Source

    enddo ; enddo ; enddo ! (i, j, k) loop

  end subroutine SA_Source_Blk


  function func (yy,rr) result (result)
    implicit none
    real(R8), intent(in) :: yy, rr
    real(R8) :: result, ff
    ff = 0.44d0 / ( 1d0 + 14d0 * yy**5d0 ) + 0.56d0
    result = yy * yy * ff - rr
  end function func

  pure function fv1 (chi) result (result)
    implicit none
    real(R8), intent(in) :: chi
    real(R8) :: result
    result = ( chi**3 ) / ( chi**3 + cv1**3 )
  end function fv1

  function fv2 (chi) result (result)
    implicit none
    real(R8) :: result
    real(R8), intent(in) :: chi
    result = 1d0 - chi / ( 1d0 + chi * fv1(chi) )
  end function fv2


  ! Compressibility correction Paciorri-Sabetta https://doi.org/10.2514/2.3967
  subroutine Compressibility_Correction ( nit, rho, chi, rhoi, p, Rgas, om, prod )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic

    implicit none
    real(R8), intent(in) :: nit, rho, chi, rhoi(nsc), p, Rgas, om
    real(R8), intent(inout) :: prod
    real(R8) :: sound, c_mach, staux, err, deltamach
    real(R8) :: aux_num, aux_den, c_plus, c_minus, res, f_cor1, f_cor2, f_cor
    
    sound = f_ss ( rhoi, p, rho, Rgas )
    c_mach = 1.d0
    staux = 25.d0*(om*nit*fv1(chi))/sound
    err = 1.d0
    do while (err>tolerance)
      deltamach=.1d0*c_mach
      aux_num=func(c_mach,staux)*2.*deltamach
      c_plus=c_mach+deltamach
      c_minus=c_mach-deltamach
      aux_den=func(c_plus,staux)-func(c_minus,staux)
      res=aux_num/aux_den
      err=dabs(res)
      c_mach=c_mach-res
    end do
    ! correction function
    f_cor1=0.6d0*(1.d0/(1.d0+9.d0*(c_mach)**6.))+0.4d0
    f_cor2=0.44d0*(1.d0/(1.d0+14.d0*(c_mach)**5.))+0.56d0
    f_cor=f_cor1*f_cor2
    prod=prod*f_cor

  end subroutine Compressibility_Correction


  subroutine Spalart_Set_Wall_Values ( mi_l, mitilde, dist )
    use MOSE_Global_m

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: mitilde          ! : rho*nit at wall
    real(R8), intent(in)                     :: mi_l             ! : laminar viscosity at wall
    real(R8), intent(in)                     :: dist             ! : cell center distance wall

    mitilde = 0d0 ! Smooth solid surface condition

  end subroutine Spalart_Set_Wall_Values


  subroutine Spalart_Set_Blowing_Wall ( rho, mil, tau, mitilde, mdot, dist )
    use MOSE_Global_m

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: mitilde
    real(R8), intent(in)                     :: rho, mil, mdot, dist
    real(R8), dimension(3), intent(in)       :: tau
    
    mitilde = 0d0 ! Smooth solid surface condition

  end subroutine Spalart_Set_Blowing_Wall


  subroutine Spalart_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use MOSE_Global_m

    implicit none
    real(R8), intent(in),  dimension(nRANS) :: prim, wall_prim  ! : boundary cell and wall variables
    real(R8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
    real(R8), intent(out), dimension(nRANS) :: ghost_prim ! : ghost cell RANS variables

    ghost_prim(1) = - prim(1) ! enforcing mit_wall == 0

  end subroutine Spalart_Extrapolate_Wall


  subroutine Spalart_Enforce_Realizability ( mitilde )

    implicit none
    real(R8), intent(inout),  dimension(:) :: mitilde
    if ( mitilde(1) < 1d-10 ) mitilde(1) = 1d-10

  end subroutine Spalart_Enforce_Realizability


  subroutine Spalart_Eddy_Viscosity ( mut, mitilde, milam, rho, Gradvel, dist )
    use MOSE_Global_m

    implicit none
    real(R8), intent(in), dimension(nRANS) :: mitilde        ! : vector of RANS eqs variables in form (rho*nitilde in this case)
    real(R8), intent(in)                   :: milam          ! : Laminar viscosity
    real(R8), intent(in)                   :: rho            ! : Density
    real(R8), intent(in)                   :: dist           ! : wall distance (not used)
    real(R8), intent(in), dimension(3,3)   :: Gradvel        ! : Velocity gradient
    real(R8), intent(out)                  :: mut            ! : Eddy viscosity
    real(R8) :: chi

    chi = mitilde(1)/milam
    mut = fv1( chi )*mitilde(1)

  end subroutine Spalart_Eddy_Viscosity


  subroutine Spalart_RANS_Diffusive_Flux (flux, mitilde, milam, gradvel, gradnit, rho, area, normal, dist)
    use MOSE_Global_m

    implicit none
    real(R8), intent(in), dimension(nrans)   :: mitilde               ! : vector of RANS eqs, in this case (rho*nitilde)
    real(R8), intent(in)                     :: milam                 ! : laminar viscosity
    real(R8), intent(in), dimension(3,3)     :: gradvel               ! : Velocity gradient (not used)
    real(R8), intent(in), dimension(nrans,3) :: gradnit               ! : Gradient in form dnitilde/dxj
    real(R8), intent(in)                     :: rho                   ! : Density
    real(R8), intent(in)                     :: area, normal(3), dist ! : Metrics
    real(R8), intent(out), dimension(nrans)  :: flux                  ! : Diffusive flux

    flux(1) = area*(milam + mitilde(1))*dot_product(gradnit(1,:), normal)/sigma

  end subroutine Spalart_RANS_Diffusive_Flux
  

end module MOSE_Lib_Spalart
