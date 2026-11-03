!> @brief SSG/LRR-omega differential Reynolds-stress model (Eisfeld, Rumsey & Togiti 2016), as specified
!> by the NASA Turbulence Modeling Resource (TMR): https://tmbwg.github.io/turbmodels/rsm-ssglrr.html
!> (formerly https://turbmodels.larc.nasa.gov/rsm-ssglrr.html).
!>
!> Seven transport equations: the Reynolds stresses rho*R_ij = mean(rho u_i'' u_j''), in the order
!> [R11 R22 R33 R12 R13 R23], and rho*omega. The pressure-strain model and the omega equation blend
!> the Launder-Reece-Rodi (LRR) coefficients near walls with the Speziale-Sarkar-Gatski (SSG) ones
!> away from them, through Menter's function F1.
!>
!> Variants, by the turbulence key (TMR names):
!>   SSGLRR           SSGLRR-RSM-w2012     generalized gradient diffusion (Daly & Harlow)
!>   SSGLRR-SD        SSGLRR-RSM-w2012-SD  simple (eddy-viscosity) diffusion
!>   SSGLRR-w2019     SSGLRR-RSM-w2019     w2012 with the length-scale correction of Eisfeld & Rumsey
!>                                         (2020) on the destruction of omega
!>   SSGLRR-w2019-SD  SSGLRR-RSM-w2019-SD  w2019 with simple diffusion
!>
!> Coupling with the mean flow (TMR, "Implementing Turbulence Models into the Compressible RANS
!> Equations", https://tmbwg.github.io/turbmodels/implementrans.html): the momentum and energy
!> equations carry -rho*R_ij in place of the Boussinesq stress (Stress_Vector_RSM), and the
!> turbulent heat and species fluxes use mu_t = rho*k/omega with Prt and Sct. MOSE advances the
!> energy of the mean flow, E = e + u_i u_i/2, without k: the equation for it is the total-energy
!> equation of TMR minus the k equation (half the trace of the stress equations), so the diffusion
!> of k cancels and the source of k, production - dissipation, is taken from the mean-flow energy.
!>
!> In a rotating frame the stress equations carry the Coriolis term
!> G_ij = -2 Omega_k ( e_ikm R_jm + e_jkm R_im ) and the pressure-strain model the absolute
!> rotation tensor W_ij + e_mji Omega_m (Launder, Tselepidakis & Younis 1987; Speziale 1989).
!>
!> References
!>   Eisfeld, B., Rumsey, C., Togiti, V., "Verification and Validation of a Second-Moment-Closure
!>     Model", AIAA J. 54(5), 2016, pp. 1524-1541, doi:10.2514/1.J054718; erratum AIAA J. 54(9),
!>     2016, p. 2926, doi:10.2514/1.J055336.
!>   Cecora, R.-D., Radespiel, R., Eisfeld, B., Probst, A., "Differential Reynolds-Stress Modeling
!>     for Aeronautics", AIAA J. 53(3), 2015, pp. 739-755, doi:10.2514/1.J053250.
!>   Eisfeld, B., Rumsey, C., "Length-Scale Correction for Reynolds-Stress Modeling", AIAA J. 58(4),
!>     2020, pp. 1518-1528, doi:10.2514/1.J058858 (w2019).
!>   Launder, B. E., Reece, G. J., Rodi, W., "Progress in the Development of a Reynolds-Stress
!>     Turbulence Closure", J. Fluid Mech. 68(3), 1975, pp. 537-566 (LRR).
!>   Speziale, C. G., Sarkar, S., Gatski, T. B., "Modelling the Pressure-Strain Correlation of
!>     Turbulence: an Invariant Dynamical Systems Approach", J. Fluid Mech. 227, 1991, pp. 245-272 (SSG).
!>   Daly, B. J., Harlow, F. H., "Transport Equations in Turbulence", Phys. Fluids 13, 1970,
!>     pp. 2634-2649 (generalized gradient diffusion).
!>   Menter, F. R., "Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications",
!>     AIAA J. 32(8), 1994, pp. 1598-1605 (F1, omega wall value).
!>   Launder, B. E., Tselepidakis, D. P., Younis, B. A., "A Second-Moment Closure Study of Rotating
!>     Channel Flow", J. Fluid Mech. 183, 1987, pp. 63-75 (Coriolis term).
!>   Speziale, C. G., "Turbulence Modeling in Noninertial Frames of Reference", Theor. Comput. Fluid
!>     Dyn. 1, 1989, pp. 3-19 (pressure-strain model in a rotating frame).
module MOSE_Lib_SSGLRR
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  real(R8), dimension(3,3), parameter :: &
  deltaij = reshape((/ 1d0, 0d0, 0d0, &
                       0d0, 1d0, 0d0, &
                       0d0, 0d0, 1d0  /), shape(deltaij), order=(/2,1/) )

  real(R8), parameter, public :: &
  C_mu = 0.09d0, &
  ! omega coefficients
  alphaw_w = 0.5556d0, &
  betaw_w = 0.075d0, &
  sigmaw_w = 0.5d0, &
  sigmad_w = 0d0, &
  C1_w = 1.8d0, &
  C1star_w = 0d0, &
  C2_w = 0d0, &
  C3_w = 0.8d0, &
  C3star_w = 0d0, &
  C2_LRR = 0.52d0, &
  C4_w = 0.5d0*(18d0*C2_LRR + 12d0)/11d0, &
  C5_w = 0.5d0*(-14d0*C2_LRR+20d0)/11d0, &
  D_w = 0.75d0*C_mu, &
  ! epsilon coefficients
  alphaw_e = 0.44d0, &
  betaw_e = 0.0828d0, &
  sigmaw_e = 0.856d0, &
  sigmad_e = 1.712d0, &
  C1_e = 1.7d0, &
  C1star_e = 0.9d0, &
  C2_e = 1.05d0, &
  C3_e = 0.8d0, &
  C3star_e = 0.65d0, &
  C4_e = 0.625d0, &
  C5_e = 0.2d0, &
  D_e = 0.22d0, &
  ! length-scale correction (w2019)
  karman = 0.41d0, &
  A_LSC = 31d0, &
  chiT_LSC = 1d0

contains

  subroutine SSGLRR_Source_Terms ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_rans, obj_time_scheme
    use MOSE_Global_m
    use MOSE_Mod_MPI, only: is_local_block
    use MOSE_Base_Types_m, only: MOSE_tensor_3D_type
    use MOSE_Lib_Stencil_Metrics, only: Source_Metrics
    implicit none
    type(MOSE_domain_type), intent(inout), target :: domain
    ! Local
    integer :: b
    type(MOSE_tensor_3D_type), pointer :: Mi(:,:,:), Mt(:,:,:)   ! cell metrics of the sources (Source_Metrics)
    logical :: point_implicit

    ! The point-implicit factor serves steady runs only: in time-accurate runs the sources are explicit and their
    ! stiffness enters the time step (Compute_dt, SSGLRR_Stiffness_Bound)
    point_implicit = obj_rans%point_implicit .and. .not. obj_time_scheme%time_accurate

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      call Source_Metrics ( domain % blk(b), Mi, Mt=Mt )
      call SSGLRR_Blk ( domain % blk(b) % P,       &
                        domain % blk(b) % r,       &
                        Mi, Mt,                    &
                        domain % blk(b) % vol,     &
                        domain % blk(b) % yn,      &
                        domain % blk(b) % dtlocal, &
                        domain % blk(b) % dim,     &
                        point_implicit, obj_rans%LSC )

    enddo

  end subroutine SSGLRR_Source_Terms



  subroutine SSGLRR_Blk ( Prim, Res, M, Mt, Volume, walldist, dt, n, point_implicit, LSC )
    use MOSE_Base_Types_m
    use MOSE_Global_m
    use MOSE_Lib_Fluid
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_RotatingFrame, only: obj_rot
    implicit none
    integer, intent(in) :: n(3)
    logical, intent(in) :: point_implicit, LSC
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Volume
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: walldist
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    type(MOSE_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Mt   ! turbulence variables
    real(R8), dimension(n(1), n(2), n(3)), intent(in) :: dt
    ! Local
    integer :: i, j, k, ii, jj, kk
    real(R8) :: rho, Rgas, mil, rhoRij(3,3), omega, kappa, rhoe, dist
    real(R8) :: Gradvel(3,3), Grad(8,3), Divel, Sij(3,3), Sij_star(3,3), Wij(3,3), Wframe(3,3)
    real(R8) :: omega_w(3)
    real(R8) :: F1, F2, alphaw, betaw, sigmaw, sigmad, C1, C2, C3, C4, C5, C1star, C3star
    real(R8), dimension(3,3) :: rhoP, rhoG, rhoeij, aij, C1_term, C2_term, C3_term, C4_term, C5_term, aikakj, rhoPI
    real(R8) :: rhoPkk, aklakl, aklSkl, Diff, Source(7), ratio, chi, fR, fw

    !$omp do collapse (3) private (rho, Rgas, mil, rhoRij, omega, kappa, rhoe, dist) &
    !$omp private (Gradvel, Grad, Divel, Sij, Sij_star, Wij, Wframe, omega_w) &
    !$omp private (F1, F2, alphaw, betaw, sigmaw, sigmad, C1, C2, C3, C4, C5, C1star, C3star) &
    !$omp private (rhoP, rhoG, rhoeij, aij, C1_term, C2_term, C3_term, C4_term, C5_term, aikakj, rhoPI) &
    !$omp private (rhoPkk, aklakl, aklSkl, Diff, Source, ratio, chi, fR, fw) &
    !$omp private (i,j,k, ii, jj, kk)

    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      ! Local variables
      call Co_rotot_Rtot ( Prim(1:nsc,i,j,k), rho, Rgas )
      mil = f_laminarViscosity ( Prim(1:nsc,i,j,k), Prim(np,i,j,k), rho, Rgas )
      rhoRij(1,1) = Prim(nt,i,j,k)
      rhoRij(2,2) = Prim(nt+1,i,j,k)
      rhoRij(3,3) = Prim(nt+2,i,j,k)
      rhoRij(1,2) = Prim(nt+3,i,j,k)
      rhoRij(1,3) = Prim(nt+4,i,j,k)
      rhoRij(2,3) = Prim(nt+5,i,j,k)
      rhoRij(2,1) = rhoRij(1,2)
      rhoRij(3,1) = rhoRij(1,3)
      rhoRij(3,2) = rhoRij(2,3)
      omega = Prim(nt+6,i,j,k) / rho
      kappa = 0.5d0 * ( rhoRij(1,1) + rhoRij(2,2) + rhoRij(3,3) ) / rho
      rhoe = rho * C_mu * kappa * omega
      dist = walldist(i,j,k)

      ! Velocity gradient computation
      Gradvel(:,1) = ( Prim(nu:nw,i+1,j,k) - Prim(nu:nw,i-1,j,k) ) * 0.5d0
      Gradvel(:,2) = ( Prim(nu:nw,i,j+1,k) - Prim(nu:nw,i,j-1,k) ) * 0.5d0
      Gradvel(:,3) = ( Prim(nu:nw,i,j,k+1) - Prim(nu:nw,i,j,k-1) ) * 0.5d0
      Gradvel = Matmul ( Gradvel, M(i,j,k) % c )

      ! Gradient of Rij, omega (1-7) and kappa (8)
      Grad(1:7,1) = ( Prim(nt:nt+6,i+1,j,k) / sum( Prim(1:nsc,i+1,j,k) ) - &
                      Prim(nt:nt+6,i-1,j,k) / sum( Prim(1:nsc,i-1,j,k) ) ) * 0.5d0
      Grad(1:7,2) = ( Prim(nt:nt+6,i,j+1,k) / sum( Prim(1:nsc,i,j+1,k) ) - &
                      Prim(nt:nt+6,i,j-1,k) / sum( Prim(1:nsc,i,j-1,k) ) ) * 0.5d0
      Grad(1:7,3) = ( Prim(nt:nt+6,i,j,k+1) / sum( Prim(1:nsc,i,j,k+1) ) - &
                      Prim(nt:nt+6,i,j,k-1) / sum( Prim(1:nsc,i,j,k-1) ) ) * 0.5d0
      Grad(8,1) = 0.5d0 * Sum ( Grad(1:3,1) )
      Grad(8,2) = 0.5d0 * Sum ( Grad(1:3,2) )
      Grad(8,3) = 0.5d0 * Sum ( Grad(1:3,3) )
      Grad = Matmul ( Grad, Mt(i,j,k) % c )

      ! Strain tensor Sij, Sij*, and vorticity tensor Wij
      Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)
      Sij = Strain_Tensor ( Gradvel )
      Sij_star = Sij - 1d0/3d0 * Divel * deltaij
      omega_w = Vorticity_Vector ( Gradvel )
      ! Use absolute vorticity for the C5 pressure-strain term: W_ij^abs captures
      ! the full mean rotation seen in the absolute frame (W_rel + Omega_frame tensor).
      if (obj_rot%enabled) omega_w = omega_w + 2.0_R8 * obj_rot%Omega_vec
      Wij = Vorticity_Tensor ( omega_w )

      ! Blending of coefficients
      F1 = compute_F1 ( rho=rho, w=omega, k=kappa, mil=mil, gradw=Grad(7,:), gradk=Grad(8,:), d=dist )
      F2 = 1d0 - F1
      alphaw = F1*alphaw_w + F2*alphaw_e
      betaw = F1*betaw_w + F2*betaw_e
      sigmaw = F1*sigmaw_w + F2*sigmaw_e
      sigmad = F1*sigmad_w + F2*sigmad_e
      C1 = F1*C1_w + F2*C1_e
      C2 = F1*C2_w + F2*C2_e
      C3 = F1*C3_w + F2*C3_e
      C4 = F1*C4_w + F2*C4_e
      C5 = F1*C5_w + F2*C5_e
      C1star = F1*C1star_w + F2*C1star_e
      C3star = F1*C3star_w + F2*C3star_e

      ! Length-scale correction (w2019): the destruction of omega is switched off where the turbulent
      ! length scale exceeds that of the log layer, i.e. near stagnation and reattachment points
      if ( LSC ) then
        ratio = Sqrt ( kappa ) / ( C_mu**0.25d0 * karman * omega * dist )
        chi = Max ( ( ratio - 1d0 ) * ratio**2, 0d0 )
        betaw = betaw * ( 1d0 - 0.5d0 * ( 1d0 + Tanh ( A_LSC * ( chi - chiT_LSC ) ) ) )
      end if

      ! Production term
      rhoP = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        rhoP(ii,jj) = rhoP(ii,jj) + rhoRij(ii,kk)*Gradvel(jj,kk)
      enddo
      enddo
      enddo

      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        rhoP(ii,jj) = rhoP(ii,jj) + rhoRij(jj,kk)*Gradvel(ii,kk)
      enddo
      enddo
      enddo
      rhoP = -rhoP

      ! Coriolis term in a rotating frame: G_ij = -2 Omega_k ( e_ikm R_jm + e_jkm R_im ), with
      ! e_ikm Omega_k the rotation tensor of the frame. It redistributes the stresses (G_kk = 0).
      rhoG = 0d0
      if (obj_rot%enabled) then
        Wframe = Vorticity_Tensor ( 2.0_R8 * obj_rot%Omega_vec )
        do ii = 1, 3
        do jj = 1, 3
        do kk = 1, 3
          rhoG(ii,jj) = rhoG(ii,jj) - 2d0 * ( Wframe(ii,kk)*rhoRij(jj,kk) + Wframe(jj,kk)*rhoRij(ii,kk) )
        enddo
        enddo
        enddo
      end if

      ! Dissipation term
      rhoeij = 2d0/3d0 * rhoe * deltaij

      ! Pressure-strain correlation
      aij = rhoRij/(rho*kappa) - 2d0/3d0 * deltaij
      aklakl = Sum (aij*aij)

      ! C1 term
      rhoPkk = rhoP(1,1) + rhoP(2,2) + rhoP(3,3)
      C1_term = -( C1*rhoe + 0.5d0*C1star*rhoPkk )*aij

      ! C2 term
      aikakj = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        aikakj(ii,jj) = aikakj(ii,jj) + aij(ii,kk)*aij(kk,jj)
      enddo
      enddo
      enddo
      C2_term = C2*rhoe * ( aikakj - 1d0/3d0*aklakl*deltaij )

      ! C3 term
      C3_term = ( C3 - C3star * Sqrt( aklakl ) ) * rho*kappa*Sij_star

      ! C4 term
      C4_term = 0d0
      aklSkl = Sum ( aij*Sij )
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C4_term(ii,jj) = C4_term(ii,jj) + aij(ii,kk)*Sij(jj,kk)
      enddo
      enddo
      enddo

      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C4_term(ii,jj) = C4_term(ii,jj) + aij(jj,kk)*Sij(ii,kk)
      enddo
      enddo
      enddo

      C4_term = C4_term - 2d0/3d0*aklSkl*deltaij
      C4_term = C4*rho*kappa * C4_term

      ! C5 term
      C5_term = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C5_term(ii,jj) = C5_term(ii,jj) + aij(ii,kk)*Wij(jj,kk)
      enddo
      enddo
      enddo

      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C5_term(ii,jj) = C5_term(ii,jj) + aij(jj,kk)*Wij(ii,kk)
      enddo
      enddo
      enddo
      C5_term = C5*rho*kappa * C5_term

      rhoPI = C1_term + C2_term + C3_term + C4_term + C5_term

      ! Cross Diffusion term
      Diff = sigmad * rho/omega * Max( Dot_Product ( Grad(8,:), Grad(7,:) ), 0d0 )

      ! Source terms
      Source(1) = rhoP(1,1) + rhoG(1,1) + rhoPI(1,1) - rhoeij(1,1)
      Source(2) = rhoP(2,2) + rhoG(2,2) + rhoPI(2,2) - rhoeij(2,2)
      Source(3) = rhoP(3,3) + rhoG(3,3) + rhoPI(3,3) - rhoeij(3,3)
      Source(4) = rhoP(1,2) + rhoG(1,2) + rhoPI(1,2) - rhoeij(1,2)
      Source(5) = rhoP(1,3) + rhoG(1,3) + rhoPI(1,3) - rhoeij(1,3)
      Source(6) = rhoP(2,3) + rhoG(2,3) + rhoPI(2,3) - rhoeij(2,3)
      Source(7) = alphaw*omega/kappa*0.5d0*rhoPkk - betaw*rho*omega**2 + Diff

      ! Point-implicit factors of steady runs, as in SST_Blk: each multiplies the whole increment of its row, the fluxes
      ! already in Res (Internal_Fluxes, BC_Fluxes) here and the source below, so the steady state does not depend on the
      ! step. The rates are those of the destruction terms: the return to isotropy of the stresses, C1*eps/k plus
      ! C1star*P_kk/(2k) where k is produced, and 2*beta*omega for omega.
      ! The energy row (MOSE advances the mean-flow energy, without k) has the source -1/2 the trace of the stress sources.
      ! As the energy row of SST with k-energy-coupling, before the factor it takes 1 - fR times the flux residual of the
      ! k row, half the trace of the first three rows, and below the source with the factor: the row of E_m + k then has
      ! the fluxes alone at every stage. Without the factor (fR = fw = 1) the lines give the bits of the explicit source
      fR = 1d0 ; fw = 1d0
      if ( point_implicit ) then
        fR = 1d0 / ( 1d0 + dt(i,j,k) * ( C1*C_mu*omega + 0.5d0*C1star*Max( rhoPkk, 0d0 )/(rho*kappa) ) )
        fw = 1d0 / ( 1d0 + dt(i,j,k) * 2d0*betaw*omega )
        Res(np,i,j,k) = Res(np,i,j,k) + ( 1d0 - fR ) * 0.5d0 * ( Res(nt,i,j,k) + Res(nt+1,i,j,k) + Res(nt+2,i,j,k) )
        Res(nt:nt+5,i,j,k) = Res(nt:nt+5,i,j,k) * fR
        Res(nt+6,i,j,k)    = Res(nt+6,i,j,k) * fw
      end if

      Res(nt:nt+5,i,j,k) = Res(nt:nt+5,i,j,k) - Source(1:6) * Volume(i,j,k) * fR
      Res(nt+6,i,j,k)    = Res(nt+6,i,j,k)    - Source(7)   * Volume(i,j,k) * fw

      ! Mean-flow energy: the source of k (half the trace: production - dissipation) is taken from it
      Res(np,i,j,k) = Res(np,i,j,k) + 0.5d0 * ( Source(1) + Source(2) + Source(3) ) * Volume(i,j,k) * fR

    enddo ; enddo ; enddo ! (i, j, k) loop

  end subroutine SSGLRR_Blk


  subroutine SSGLRR_Eddy_Viscosity ( mi_t, var, mi_l, rho, Gradvel, dist, k_rough )
    use MOSE_Global_m
    use MOSE_Lib_Fluid
    implicit none
    real(R8), intent(in), dimension(nRANS) :: var        ! : [rho*Rij rho*w]
    real(R8), intent(in)                   :: mi_l       ! : Molecular viscosity
    real(R8), intent(in)                   :: rho        ! : Density
    real(R8), intent(in)                   :: dist       ! : Wall distance (not used)
    real(R8), intent(in)                   :: k_rough    ! : Wall roughness (not modelled)
    real(R8), intent(in), dimension(3,3)   :: Gradvel    ! : Velocity gradient
    real(R8), intent(out)                  :: mi_t       ! : Eddy viscosity
    ! Local
    real(R8) :: rhok, w

    w = var(7) / rho
    rhok = 0.5d0 * Sum( var(1:3) )
    mi_t = rhok / w
  end subroutine SSGLRR_Eddy_Viscosity


  !> Bound of the stiffness of the destruction terms, from the cell alone (time-accurate runs, Compute_dt). With omega
  !> held, the dissipation and the return to isotropy of the stresses, -2/3 rho eps delta_ij - C1 rho eps a_ij with
  !> eps = C_mu k omega, have the rates C_mu omega on the trace (k) and C1 C_mu omega on the deviatoric part; the
  !> destruction of omega, beta rho omega^2, has 2 beta omega. C1 and beta blend the LRR and SSG values through F1, which
  !> takes gradients: their largest values bound them, and the length-scale correction (w2019) only lowers beta. The
  !> part C1star P_kk/(2k) of the return to isotropy needs the velocity gradient and is not in the bound.
  pure function SSGLRR_Stiffness_Bound ( rans_variables, rho, walldist, k_rough ) result ( stiff )
    use MOSE_Global_m
    implicit none
    real(R8), intent(in), dimension(nRANS) :: rans_variables  ! : rho*Rij (6), rho*omega
    real(R8), intent(in)                   :: rho, walldist, k_rough
    real(R8)                               :: stiff

    stiff = max ( max ( C1_w, C1_e ) * C_mu, 2d0 * max ( betaw_w, betaw_e ) ) * max ( rans_variables(7) / rho, 0d0 )

  end function SSGLRR_Stiffness_Bound


  subroutine SSGLRR_RANS_Diffusive_Flux ( Flux, var, mi_l, vGrad, Grad, rho, Area, Normal, Dist )
    use MOSE_Global_m
    implicit none
    real(R8), intent(in), dimension(nRANS)   :: var                    ! : [rho*Rij rho*w]
    real(R8), intent(in)                     :: mi_l                   ! : laminar viscosity
    real(R8), intent(in), dimension(3,3)     :: vGrad                  ! : Velocity gradient (not used)
    real(R8), intent(in), dimension(nRANS,3) :: Grad                   ! : Gradient of [ Rij w ]
    real(R8), intent(in)                     :: rho                    ! : Density
    real(R8), intent(in)                     :: Area, Normal(3), Dist  ! : Metrics
    real(R8), intent(out), dimension(nRANS)  :: Flux                   ! : Diffusive flux
    ! Local
    integer :: i
    real(R8) :: omega, kappa, Gradk(3), F1, F2, D_, sigmaw, Diff(6,3)

    omega = var(7)/rho
    kappa = 0.5d0 * Sum ( var(1:3) ) / rho
    Gradk = 0.5d0 * ( Grad(1,:) + Grad(2,:) + Grad(3,:) )

    ! Blending of coefficients
    F1 = compute_F1 ( rho=rho, w=omega, k=kappa, mil=mi_l, gradw=Grad(7,:), gradk=Gradk, d=Dist )
    F2 = 1d0 - F1
    D_ = ( D_w * F1 + D_e * F2 ) / ( C_mu*omega )
    sigmaw = F1 * sigmaw_w + F2 * sigmaw_e

    Diff(:,1) = mi_l*Grad(1:6,1) + D_*( var(1)*Grad(1:6,1) + var(4)*Grad(1:6,2) + var(5)*Grad(1:6,3) )
    Diff(:,2) = mi_l*Grad(1:6,2) + D_*( var(4)*Grad(1:6,1) + var(2)*Grad(1:6,2) + var(6)*Grad(1:6,3) )
    Diff(:,3) = mi_l*Grad(1:6,3) + D_*( var(5)*Grad(1:6,1) + var(6)*Grad(1:6,2) + var(3)*Grad(1:6,3) )

    do i = 1, 6
      Flux(i) = Dot_Product ( [ Diff(i,1), Diff(i,2), Diff(i,3) ], Normal )
    end do

    Flux(7) = ( mi_l + sigmaw*rho*kappa/omega ) * Dot_Product ( Grad(7,:), Normal )
    Flux = Flux * Area

  end subroutine SSGLRR_RANS_Diffusive_Flux


  subroutine SSGLRR_SD_RANS_Diffusive_Flux ( Flux, var, mi_l, vGrad, Grad, rho, Area, Normal, Dist )
    use MOSE_Global_m
    implicit none
    real(R8), intent(in), dimension(nRANS)   :: var                    ! : [rho*Rij rho*w]
    real(R8), intent(in)                     :: mi_l                   ! : laminar viscosity
    real(R8), intent(in), dimension(3,3)     :: vGrad                  ! : Velocity gradient (not used)
    real(R8), intent(in), dimension(nRANS,3) :: Grad                   ! : Gradient of [ Rij w ]
    real(R8), intent(in)                     :: rho                    ! : Density
    real(R8), intent(in)                     :: Area, Normal(3), Dist  ! : Metrics
    real(R8), intent(out), dimension(nRANS)  :: Flux                   ! : Diffusive flux
    ! Local
    integer :: i
    real(R8) :: omega, kappa, Gradk(3), F1, F2, D_, sigmaw, mi_t, Diff(6,3)

    omega = var(7)/rho
    kappa = 0.5d0 * Sum ( var(1:3) ) / rho
    Gradk = 0.5d0 * ( Grad(1,:) + Grad(2,:) + Grad(3,:) )

    ! Blending of coefficients
    F1 = compute_F1 (rho=rho, w=omega, k=kappa, mil=mi_l, gradw=Grad(7,:), gradk=Gradk, d=Dist)
    F2 = 1d0 - F1
    D_ = 0.5d0*C_mu*F1 + 2d0/3d0*0.22d0*F2
    sigmaw = F1 * sigmaw_w + F2 * sigmaw_e
    mi_t = rho*kappa/omega

    Diff(1:6,:) = (mi_l + D_/C_mu*mi_t) * Grad(1:6,:)
    do i = 1, 6
      Flux(i) = Dot_Product ( [ Diff(i,1), Diff(i,2), Diff(i,3) ], Normal )
    end do
    Flux(7) = ( mi_l + sigmaw*mi_t ) * Dot_Product ( Grad(7,:), Normal )
    Flux = Flux * Area

  end subroutine SSGLRR_SD_RANS_Diffusive_Flux

  subroutine SSGLRR_Set_Wall_Values ( mi_l, var_cell, var, dist, k_rough )
    use MOSE_Global_m
    use MOSE_Config_Types_m, only: obj_rans
    implicit none
    real(R8), intent(in),  dimension(nRANS)  :: var_cell         ! : boundary cell values (not used)
    real(R8), intent(out), dimension(nRANS)  :: var              ! : [ rho*Rij rho*w ] at wall
    real(R8), intent(in)                     :: mi_l             ! : laminar viscosity at wall
    real(R8), intent(in)                     :: dist             ! : cell center distance from wall
    real(R8), intent(in)                     :: k_rough          ! : wall roughness (not modelled)

    var(1:6) = 0d0 ! Solid surface condition on rhoRij
    ! Smooth-surface BC: omega_wall = C*mi_l/dist**2, with C set by omega-wall-bc
    ! (800 = Menter practical, the condition of TMR; 80 = asymptotic y->0 limit). See Setup_RANS_Model.
    var(7) = obj_rans%omega_wall_coef * mi_l / ( dist**2 )

  end subroutine SSGLRR_Set_Wall_Values


  subroutine SSGLRR_Set_Blowing_Wall ( rho, mil, tau, var, mdot, dist )
    use MOSE_Global_m
    use MOSE_Config_Types_m, only: obj_rans
    implicit none
    real(R8), intent(out), dimension(nRANS)  :: var
    real(R8), intent(in)                     :: rho, mil, mdot, dist
    real(R8), dimension(3), intent(in)       :: tau
    ! Local
    real(R8) :: tauw, uT, v_plus, SB

    var(1:6) = 0d0 ! Solid surface condition on rhoRij
    if ( obj_rans%blowing_corr ) then
      tauw = sqrt( sum(tau**2) )
      uT = sqrt ( tauw / rho )
      v_plus = abs(mdot)/rho / uT
      SB = 25d0 / ( v_plus*( 1d0 + 5d0 * v_plus ) )
      var(7) = (rho*uT)**2 * SB / mil ! Wilcox BC for mass injection
    else
      var(7) = obj_rans%omega_wall_coef * mil / ( dist**2 )
    end if

  end subroutine SSGLRR_Set_Blowing_Wall


  subroutine SSGLRR_Extrapolate_Wall ( var, wall_var, rho, wall_rho, ghost_var )
    use MOSE_Global_m
    implicit none
    real(R8), intent(in),  dimension(nRANS) :: var, wall_var  ! : boundary cell and wall variables
    real(R8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
    real(R8), intent(out), dimension(nRANS) :: ghost_var ! : ghost cell RANS variables

    ghost_var(1:6) = - var(1:6) ! enforcing rhoRij_wall == 0
    ghost_var(7) = rho * 2d0 * wall_var(7)/wall_rho - var(7) ! extrapolating from omega_wall

  end subroutine SSGLRR_Extrapolate_Wall


  subroutine SSGLRR_Enforce_Realizability ( rhoRij )
    implicit none
    real(R8), intent(inout),  dimension(:) :: rhoRij

    ! 1st condition: u'u', v'v' and w'w' > 0
    if ( rhoRij(1) < 1d-16 ) rhoRij(1) = 1d-16
    if ( rhoRij(2) < 1d-16 ) rhoRij(2) = 1d-16
    if ( rhoRij(3) < 1d-16 ) rhoRij(3) = 1d-16

    ! 2nd condition: abs( u'v' ) < sqrt ( u'u' * v'v' ), etc
    if ( Abs( rhoRij(4) ) > Sqrt( rhoRij(1)*rhoRij(2) ) ) then
      rhoRij(4) = Sign( 1d0, rhoRij(4) ) * Sqrt( rhoRij(1)*rhoRij(2) )
    end if

    if ( Abs( rhoRij(5) ) > Sqrt( rhoRij(1)*rhoRij(3) ) ) then
      rhoRij(5) = Sign( 1d0, rhoRij(5) ) * Sqrt( rhoRij(1)*rhoRij(3) )
    end if

    if ( Abs( rhoRij(6) ) > Sqrt( rhoRij(2)*rhoRij(3) ) ) then
      rhoRij(6) = Sign( 1d0, rhoRij(6) ) * Sqrt( rhoRij(2)*rhoRij(3) )
    end if

    ! omega > 0
    if ( rhoRij(7) < 1d-10 ) rhoRij(7) = 1d-10

  end subroutine SSGLRR_Enforce_Realizability


  function compute_F1 ( rho, w, k, mil, Gradw, Gradk, d ) result ( F1 )

    implicit none
    real(R8), intent(in) :: rho, w, k, mil, Gradw(3), Gradk(3), d
    real(R8) :: CD, zeta, F1

    CD = sigmad_e * rho / w * Max ( Dot_Product( Gradw, Gradk ), 1d-20 )
    zeta = Min ( Max( Sqrt(k)/(C_mu*w*d), 5d2*mil/(rho*w*d**2) ), &
                 4d0*sigmaw_e*rho*k/(CD*d**2) )
    F1 = Tanh ( zeta**4 )

  end function compute_F1

end module MOSE_Lib_SSGLRR
