module MOSE_Lib_Diffusive
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Diffusive_Flux, Compute_Diffusive_Flux, Face_Metric, Inverse3, Visc_State

contains

  !> Viscous flux of the face between cells 1 and 2, along the normal from 1 to 2: the difference across the face (index
  !> direction a) and the mean central differences of the two cells along the face (b: cells 3..6, c: cells 7..10) of
  !> their viscous variables Visc1..Visc10 (Visc_State, which Viscous_Variables computes once per stage for every cell),
  !> mapped to Cartesian components with the face metric Mf (MOSE_Lib_Stencil_Metrics: Viscous_Face_Metric); the state
  !> of the face is the mean of the primitive variables Prim1 and Prim2 of the two cells. Res1 takes -flux, Res2 +flux.
  !> With a ghost cell of a wall or of a mapped face in the stencil, Mg: one face metric per group of quantities (1 species
  !> and scalars, 2 temperature, 3 velocity, 4 turbulence variables), whose ghost values sit at different positions
  !> (Ghost_Wall, Ghost_Q2D_Mapped). Behind a mapped face the soot, a copy of the cell, sits at x1, while the metric of
  !> group 1 puts it with the species of the file at 2 x_f - x1: no flux reads its gradient (Soot_Diffusive_Flux takes that
  !> of the temperature).
  subroutine Diffusive_Flux ( normal, area, waldis1, waldis2, kr1, kr2, Prim1, Prim2, Visc1, Visc2, Visc3, &
                              Visc4, Visc5, Visc6, Visc7, Visc8, Visc9, Visc10, Mf, Res1, Res2, a, b, c, &
                              Sc, Sct, Prt, Prl, soot_enabled, Mg )
    use MOSE_Global_m
    implicit none
    integer, intent(in)  :: a, b, c
    logical, intent(in)  :: soot_enabled
    real(R8), intent(in) :: Sc, Sct, Prt, Prl
    real(R8), intent(in) :: normal(3), area, waldis1, waldis2
    real(R8), intent(in) :: kr1, kr2      ! nearest-wall roughness of the two cells
    real(R8), intent(in), dimension(nprim) :: Prim1, Prim2
    real(R8), intent(in), dimension(nprim) :: Visc1, Visc2, Visc3, Visc4, Visc5, Visc6
    real(R8), intent(in), dimension(nprim) :: Visc7, Visc8, Visc9, Visc10
    real(R8), intent(in), dimension(3,3) :: Mf       ! face metric
    real(R8), intent(in), optional :: Mg(3,3,4)      ! face metrics per group of quantities
    real(R8), intent(inout), dimension(nprim) :: Res1, Res2
    ! Local
    real(R8) :: Gradient(nprim,3), Prim(nprim), waldis, k_rough, Flux(nprim)
    integer  :: v, g
    real(R8) :: g1, g2, g3

    ! Gradient in the same direction of the face: 1 and 2 (mass fractions, velocity, temperature, turbulence variables)
    Gradient ( :, a ) = Visc2 - Visc1
    if (nprim>np) then
      waldis = 0.5d0 * ( waldis1 + waldis2 ) ! distance to nearest wall
      k_rough = 0.5d0 * ( kr1 + kr2 )        ! roughness of nearest wall
    end if

    ! Gradient in tangential directions: 3-10
    call Tangential_Gradient ( Visc3, Visc4, Visc5, Visc6,  Gradient(:,b) )
    call Tangential_Gradient ( Visc7, Visc8, Visc9, Visc10, Gradient(:,c) )

    !> Transform the gradients from computational to physical space.
    !> Explicit loop, not `Gradient = matmul(Gradient,M)`: Gradient aliases
    !> itself there, so the compiler must build the whole (nprim,3) result in a
    !> temporary on every face.
    if ( present(Mg) ) then
      do v = 1, nprim
        g = 1                                        ! species, soot and passive scalars (behind a mapped face the soot
                                                     ! sits at x1 and this metric at 2 x_f - x1: no flux reads its gradient)
        if ( v == np ) g = 2                         ! temperature
        if ( v >= nu .and. v <= nw ) g = 3           ! velocity
        if ( nrans > 0 .and. v >= nt ) g = 4         ! turbulence variables
        g1 = Gradient(v,1) ; g2 = Gradient(v,2) ; g3 = Gradient(v,3)
        Gradient(v,1) = g1*Mg(1,1,g) + g2*Mg(2,1,g) + g3*Mg(3,1,g)
        Gradient(v,2) = g1*Mg(1,2,g) + g2*Mg(2,2,g) + g3*Mg(3,2,g)
        Gradient(v,3) = g1*Mg(1,3,g) + g2*Mg(2,3,g) + g3*Mg(3,3,g)
      end do
    else
      do v = 1, nprim
        g1 = Gradient(v,1) ; g2 = Gradient(v,2) ; g3 = Gradient(v,3)
        Gradient(v,1) = g1*Mf(1,1) + g2*Mf(2,1) + g3*Mf(3,1)
        Gradient(v,2) = g1*Mf(1,2) + g2*Mf(2,2) + g3*Mf(3,2)
        Gradient(v,3) = g1*Mf(1,3) + g2*Mf(2,3) + g3*Mf(3,3)
      end do
    end if

    Prim = 0.5d0 * ( Prim1 + Prim2 )
    call Compute_Diffusive_Flux ( Prim, Gradient, area, normal, waldis, k_rough, Flux, Sc, Sct, Prt, Prl, soot_enabled )

    Res1 = Res1 - Flux
    Res2 = Res2 + Flux

  end subroutine Diffusive_Flux


  !> Metric tensor at the face between two cells of metrics M1 and M2:
  !> inverse of the mean of the inverse metrics (edge vectors dx/dxi).
  pure function Face_Metric ( M1, M2 ) result ( M )
    implicit none
    real(R8), intent(in) :: M1(3,3), M2(3,3)
    real(R8)             :: M(3,3)

    M = Inverse3 ( 0.5d0 * ( Inverse3 ( M1 ) + Inverse3 ( M2 ) ) )

  end function Face_Metric


  !> Inverse of a 3x3 matrix by its adjugate (transposed cofactors). Works for M
  !> as stored by Compute_Metric_Tensor (cofactors/det = transposed inverse of the
  !> edge vectors): transposition commutes with averaging and inversion.
  pure function Inverse3 ( A ) result ( B )
    implicit none
    real(R8), intent(in) :: A(3,3)
    real(R8)             :: B(3,3)
    real(R8) :: det

    B(1,1) =  ( A(2,2)*A(3,3) - A(2,3)*A(3,2) )
    B(1,2) = -( A(1,2)*A(3,3) - A(1,3)*A(3,2) )
    B(1,3) =  ( A(1,2)*A(2,3) - A(1,3)*A(2,2) )
    B(2,1) = -( A(2,1)*A(3,3) - A(2,3)*A(3,1) )
    B(2,2) =  ( A(1,1)*A(3,3) - A(1,3)*A(3,1) )
    B(2,3) = -( A(1,1)*A(2,3) - A(1,3)*A(2,1) )
    B(3,1) =  ( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
    B(3,2) = -( A(1,1)*A(3,2) - A(1,2)*A(3,1) )
    B(3,3) =  ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )
    det = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
    B = B / det

  end function Inverse3


  !> Mean of the central differences of two cells along a tangential direction, from viscous variables: Visc1 and Visc2
  !> the neighbours of the first cell (at -1 and +1), Visc3 and Visc4 those of the second
  subroutine Tangential_Gradient ( Visc1, Visc2, Visc3, Visc4, Gradient )
    use MOSE_Global_m
    implicit none
    real(R8), intent(in), dimension(nprim) :: Visc1, Visc2, Visc3, Visc4
    real(R8), intent(out), dimension(nprim) :: Gradient

    Gradient = ( Visc2 - Visc1 + Visc4 - Visc3 ) * 0.25d0

  end subroutine Tangential_Gradient


  !> Viscous variables of a cell of primitive variables P: mass fractions, velocity, temperature and turbulence
  !> variables per unit mass, the quantities whose differences make the viscous gradients. Viscous_Variables computes
  !> them once per stage for every cell (blk%Visc); the faces between blocks for the cells of the other block (bc%Pg).
  subroutine Visc_State ( P, Visc )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot
    implicit none
    real(R8), intent(in)  :: P(nprim)
    real(R8), intent(out) :: Visc(nprim)
    ! Local
    real(R8) :: rho, Rgas

    call co_rotot_Rtot ( P(1:nsc), rho, Rgas )
    Visc(1:nsc) = P(1:nsc) / rho
    Visc(nu:nw) = P(nu:nw)
    Visc(np) = P(np) / ( rho * Rgas )
    if (nprim>np) Visc(np+1:nprim) = P(np+1:nprim) / rho

  end subroutine Visc_State


  subroutine Compute_Diffusive_Flux ( Prim, Gradient, area, normal, waldis, k_rough, Flux, Sc, Sct, Prt, Prl, soot_enabled )
    use MOSE_Global_m
    use MOSE_Lib_Fluid
    use MOSE_Lib_RANS
    use MOSE_Mod_Soot, only: Soot_Diffusive_Flux
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: Prim(nprim), Gradient(nprim,3), area, normal(3), waldis, k_rough
    real(R8), intent(in)  :: Sc, Sct, Prt, Prl
    logical, intent(in)   :: soot_enabled
    real(R8), intent(out) :: Flux(nprim)
    ! Local
    integer :: s, T_i, Tint(2)
    real(R8) :: rho, Rgas, T, Tdiff, cp, mil, kl, mie, kappa
    real(R8) :: Dm(nsc), stress(3), DiffHFlux, DmGradYi
    !> Contiguous copy of the velocity-gradient rows, shared by the three
    !> consumers below.  Passing `Gradient(nu:nw,:)` directly is a strided
    !> section, so each consumer would build its own temporary instead.
    real(R8) :: VelGrad(3,3)

    VelGrad(1,1) = Gradient(nu,1) ; VelGrad(1,2) = Gradient(nu,2) ; VelGrad(1,3) = Gradient(nu,3)
    VelGrad(2,1) = Gradient(nv,1) ; VelGrad(2,2) = Gradient(nv,2) ; VelGrad(2,3) = Gradient(nv,3)
    VelGrad(3,1) = Gradient(nw,1) ; VelGrad(3,2) = Gradient(nw,2) ; VelGrad(3,3) = Gradient(nw,3)

    ! Thermodynamic and transport properties at the interface
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    T_i     = int(T)
    Tdiff   = T - T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1
    cp = f_cp_expr ( Prim(1:nsc), Tint, Tdiff, rho )
    call co_k_mi_lam_Wilke_expr ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )
    ! Optional unity/forced laminar Prandtl: override the mixture conductivity with k = mu*cp/Prl
    if ( Prl > 0d0 ) kl = mil * cp / Prl

    ! Eddy viscosity
    mie = 0d0
    if (model==2) then
      call Eddy_Viscosity ( mut=mie, rans_variables=Prim(nt:nprim), &
                            mul=mil, rho=rho, vel_gradient=VelGrad, &
                            walldist=waldis, k_rough=k_rough )
    end if

    ! Species diffusion coefficients. Sc<=0 selects mixture-averaged multicomponent
    ! diffusion (D_k from the binary-diffusion table), mirroring the Prl<=0 idiom above.
    if (Sc <= 0d0) then
      ! Mixture-averaged multicomponent: per-species laminar D_k from binary-diffusion table,
      ! rescaled internally from the table reference pressure to the local pressure Prim(np).
      call co_DS_expr ( Prim(1:nsc), rho, Tint, Tdiff, Prim(np), Dm )
      if (mie > 0d0) Dm(1:nsc) = Dm(1:nsc) + mie/(rho*Sct) ! add turbulent (constant-Sct) part
    else
      ! Constant Schmidt: same coefficient for all species (Dm = mu/(rho*Sc))
      Dm (1:nsc) = ( mil/Sc + mie/Sct ) / rho
    end if

    kappa = kl + mie*cp/Prt ! Laminar + turbulent conductivity

    Stress = stress_vector ( VelGrad, normal, mil, mie, prim(nt:) ) ! Stress tensor in cartesian components

    ! Fluxes computation
    DiffHFlux = 0.0d0
    Flux = 0.0d0

    ! Diffusive mass fluxes
    do s = 1, nsc
      DmGradYi  = Dm(s) * dot_product ( Gradient(s,:), normal )     ! Dm * dyi/dn
      Flux(s) = rho * area * DmGradYi                               ! rho * Dm *dyi/dn * Area
      DiffHFlux = DiffHFlux + DmGradYi * f_tabT_expr ( s, h_tab, Tint, Tdiff )
    enddo

    ! Correcting the mass flux for possible errors (sum must be zero)
    Flux(1:nsc) = Flux(1:nsc) - sum(Flux(1:nsc))*Prim(1:nsc)/rho

    ! Momentum flux
    Flux(nu:nw) = area * Stress

    ! Diffusive enthalpy flux
    Flux(np) = area * ( dot_product ( Stress, Prim(nu:nw) ) + &
                kappa * dot_product ( Gradient(np,:), normal ) + DiffHFlux * rho )

    ! Turbulence variables diffusive flux
    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim(nt:nprim), &
                                 vel_gradient=VelGrad, &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho, &
                                 area=area, normal=normal, dist=waldis )
    end if

    ! Soot diffusive fluxes
    if (soot_enabled) then
      call Soot_Diffusive_Flux ( flux=Flux, &
                                 prim=Prim, &
                                 T_gradient=Gradient(np,:), &
                                 mu=mil+mie, &
                                 rho=rho, &
                                 T=T, &
                                 area=area, & 
                                 normal=normal )
    end if

  end subroutine Compute_Diffusive_Flux

end module MOSE_Lib_Diffusive