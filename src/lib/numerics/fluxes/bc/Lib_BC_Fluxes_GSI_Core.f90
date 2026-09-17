module MOSE_Lib_BC_Fluxes_GSI_Core
  !< Shared machinery for the ablative wall boundary conditions: melting (BC 503),
  !< pyrolysis (BC 504) and surface reactions (BC 505).
  !<
  !< All three model a *blowing* wall, where a gas-surface interaction (GSI) model
  !< injects mass into the domain. They share the same three steps:
  !<
  !<  1. **wall gas state** -- from the boundary layer assumption
  !<     `p(wall) = p(cell)`, `R(wall) = R(cell)` and `ci(wall) = ci(cell)`, so that
  !<     only the density changes with the wall temperature (`Wall_Gas_State`);
  !<  2. **gas-surface interaction** -- a GSI model returns the species source term
  !<     `omega`, the heat absorbed by the surface `q_pyro` and the blown mass flux
  !<     `mdot` (`GSI_Source`);
  !<  3. **flux assembly** -- mass, momentum and energy fluxes are built from the
  !<     blowing velocity and the viscous stress, the RANS variables are set at the
  !<     wall and extrapolated to the ghost cell, and the result is subtracted from
  !<     the residual (`Ablation_Wall_Fluxes`).
  !<
  !< The boundary conditions differ only in how `T_wall` and `mdot` are obtained:
  !<
  !<  * **melting** -- `T_wall` is prescribed and `mdot` follows from the surface
  !<    energy balance, so a single pass through the three steps is enough;
  !<  * **pyrolysis / surface reactions** -- `mdot` follows from an Arrhenius law in
  !<    `T_wall`, and `T_wall` is itself the root of the surface energy balance
  !<    `q_conv - q_pyro - q_rad = 0`, found by the secant iteration in
  !<    `Solve_Wall_Temperature`.
  !<
  !< @note Known divergences between the two families:
  !<
  !<  * melting reports `q_conv + q_rad` as the heat absorbed by the surface, while
  !<    pyrolysis and surface reactions report `q_conv` alone. Either way the gas
  !<    loses only `q_conv` -- see the `q_gas` and `q_wall` arguments of
  !<    `Ablation_Wall_Fluxes` at each call site;
  !<  * melting sets the RANS wall values with `RANS_Set_Wall_Values`, i.e. with the
  !<    no-blowing correlation, whereas pyrolysis and surface reactions use the
  !<    blowing-wall correlation `RANS_Set_Blowing_Wall`.

  use iso_fortran_env, only: I4 => int32, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties
  use MOSE_Lib_Fluid
  use MOSE_Lib_GSI
  use MOSE_Lib_RANS
  use FLINT_Lib_Thermodynamic

  implicit none
  private

  public :: wall_face_t, Get_Wall_Face
  public :: Wall_Gas_State, GSI_Source, Solve_Wall_Temperature, Ablation_Wall_Fluxes
  public :: Surface_Is_Ablating
  public :: BC_Wall_Ablation
  public :: GSI_KIND_PYROLYSIS, GSI_KIND_REACTIONS
  public :: GSI_PYRO_HTPB, GSI_PYRO_HDPE, GSI_PYRO_PP, GSI_REAC_CARBON_BRADLEY
  public :: RANS_WALL_NO_BLOWING, RANS_WALL_BLOWING

  !> Family of gas-surface interaction model, selecting the GSI dispatch.
  integer, parameter :: GSI_KIND_PYROLYSIS = 1  !< Arrhenius pyrolysis of a solid fuel
  integer, parameter :: GSI_KIND_REACTIONS = 2  !< heterogeneous surface reactions

  !> Pyrolysis models, matching `bc % GSI_pyro_model_id`.
  integer, parameter :: GSI_PYRO_HTPB = 1
  integer, parameter :: GSI_PYRO_HDPE = 2
  integer, parameter :: GSI_PYRO_PP   = 3

  !> Surface reaction models, matching `bc % GSI_surf_reac_id`.
  integer, parameter :: GSI_REAC_CARBON_BRADLEY = 1

  !> Correlation used to set the RANS variables at the wall.
  integer, parameter :: RANS_WALL_NO_BLOWING = 1  !< standard wall values
  integer, parameter :: RANS_WALL_BLOWING    = 2  !< corrected for surface blowing

  !> The fluxes are evaluated *on* the wall, so the wall-to-cell distance used by the
  !> RANS diffusive flux degenerates to zero; a tiny positive value avoids the 0/0.
  real(R8), parameter :: WALL_DISTANCE = 1d-20

  !> Relative temperature perturbation used to build the secant slope of the wall
  !> energy balance.
  real(R8), parameter :: T_WALL_PERTURBATION = 1d-4

  !> Iterations of the wall temperature solver. The loop always runs to completion:
  !> the residual-based stopping criterion below is currently disabled, so each
  !> boundary face costs `2 * MAX_T_WALL_ITER` GSI and transport evaluations.
  integer, parameter :: MAX_T_WALL_ITER = 20

  !> Tolerance on the wall energy balance. Unused while the early exit in
  !> `Solve_Wall_Temperature` stays commented out.
  real(R8), parameter :: T_WALL_BALANCE_TOL = 1d-5

  !> Geometry and orientation of a boundary face, gathered once per call.
  type :: wall_face_t
    integer  :: dir        !< index of the face-normal direction (1, 2 or 3)
    integer  :: modfm2     !< +1 on a low face, -1 on a high face: outward sign
    integer  :: modfm3     !< one-sided finite-difference factor for the wall gradient
    integer  :: ghost(3)   !< indices of the ghost cell behind the face
    real(R8) :: normal(3)  !< face normal
    real(R8) :: area       !< face area
    real(R8) :: M(3,3)     !< cell metric, maps computational to physical gradients
  end type wall_face_t

contains

  logical function Surface_Is_Ablating ( face, mdot, T_surface ) result ( ablating )
    !< Whether the gas-surface interaction is actually consuming the wall material.
    !<
    !< `mdot` is the orientation-signed rate the boundary condition works with;
    !< undoing that sign recovers the physical ablation rate, which has to be
    !< strictly positive for there to be anything to apply. It is not, whenever the
    !< flow cannot sustain the process: a melting wall losing more heat to the gas
    !< than it receives gets a non-positive `mdot` straight out of its energy
    !< balance, and a wall temperature solver that failed to converge leaves
    !< non-finite values behind. In either case the surface is inert and the caller
    !< must fall back to an impermeable wall rather than apply a meaningless flux.
    implicit none
    type(wall_face_t), intent(in) :: face
    real(R8), intent(in) :: mdot       !< blown mass flux, signed by face orientation
    real(R8), intent(in) :: T_surface  !< wall temperature the models settled on

    ablating = ieee_is_finite ( mdot ) .and. ieee_is_finite ( T_surface ) .and. &
               T_surface > 0d0 .and. mdot * (-face % modfm2) > 0d0

  end function Surface_Is_Ablating


  subroutine Get_Wall_Face ( Im, Jm, Km, Fm, Blk, face )
    !< Collect the metric and orientation data of boundary face `Fm` of cell `(Im,Jm,Km)`.
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(in) :: Blk
    type(wall_face_t), intent(out) :: face
    ! Local
    integer :: modfm, modfm1, face_i, face_j, face_k

    call Compute_Modfm ( Fm, modfm, modfm1, face % modfm2, face % modfm3 )
    call Face_Index ( Fm, face % dir, Im, Jm, Km, face_i, face_j, face_k )

    face % normal = Blk % dir(face % dir) % f(face_i,face_j,face_k) % n
    face % area   = Blk % dir(face % dir) % f(face_i,face_j,face_k) % a
    face % M      = Blk % M(Im,Jm,Km) % c
    face % ghost  = [ Im, Jm, Km ] - guide(Fm,1:3)

  end subroutine Get_Wall_Face


  subroutine Wall_Gas_State ( face, Prim, rho, Rgas, T_cell, T_wall, &
                              rho_wall, roi_wall, mu_wall, k_wall, grad_T, q_conv, y_wall )
    !< Gas state at the wall for a given wall temperature.
    !<
    !< The boundary layer assumption keeps pressure, gas constant and mass fractions
    !< equal to their boundary cell values, so only the density - and with it the
    !< partial densities and the transport properties - responds to `T_wall`.
    !<
    !< A surface that imposes its own composition on the gas touching it - a burning
    !< propellant grain, whose flame products are what the wall actually sees - passes
    !< `y_wall` to override the inherited mass fractions. The pressure and the gas
    !< constant still come from the boundary cell, so only the partial densities and
    !< the transport properties change.
    implicit none
    type(wall_face_t), intent(in) :: face
    real(R8), intent(in)  :: Prim(nprim)   !< boundary cell primitive variables
    real(R8), intent(in)  :: rho, Rgas     !< boundary cell density and gas constant
    real(R8), intent(in)  :: T_cell        !< boundary cell temperature
    real(R8), intent(in)  :: T_wall        !< wall temperature
    real(R8), intent(out) :: rho_wall      !< wall density
    real(R8), intent(out) :: roi_wall(nsc) !< wall partial densities
    real(R8), intent(out) :: mu_wall       !< wall laminar viscosity
    real(R8), intent(out) :: k_wall        !< wall thermal conductivity
    real(R8), intent(out) :: grad_T(3)     !< wall temperature gradient
    real(R8), intent(out) :: q_conv        !< conductive heat flux at the wall
    real(R8), intent(in), optional :: y_wall(nsc) !< mass fractions imposed by the surface

    rho_wall = Prim(np) / ( Rgas * T_wall )
    if ( present(y_wall) ) then
      roi_wall = y_wall * rho_wall
    else
      roi_wall = Prim(1:nsc) * rho_wall / rho
    end if

    ! One-sided difference across the face; the tangential contributions vanish.
    grad_T = 0d0
    grad_T(face % dir) = ( T_cell - T_wall ) * face % modfm3
    grad_T = matmul ( grad_T, face % M )

    ! The grain combustion wall is the one ablative boundary that also runs inviscid,
    ! where the transport tables are not available and there is no conduction anyway.
    if ( model > 0 ) then
      call co_k_mi_lam_Wilke ( roi_wall, rho_wall, T_wall, mu_wall, k_wall )
      q_conv = k_wall * dot_product ( grad_T, face % normal )
    else
      mu_wall = 0d0
      k_wall  = 0d0
      q_conv  = 0d0
    end if

  end subroutine Wall_Gas_State


  subroutine GSI_Source ( gsi_kind, gsi_id, face, roi_wall, T_wall, omega, q_pyro, mdot, y_wall )
    !< Evaluate the selected gas-surface interaction model.
    !<
    !< The mass flux and the species source term are returned already signed with
    !< respect to the boundary cell, so that a positive `mdot` blows into the domain
    !< whichever side of the block the wall sits on.
    implicit none
    integer, intent(in) :: gsi_kind                !< GSI family, one of `GSI_KIND_*`
    integer, intent(in) :: gsi_id                  !< model within the family
    type(wall_face_t), intent(in) :: face
    real(R8), intent(in)  :: roi_wall(nsc)         !< wall partial densities
    real(R8), intent(in)  :: T_wall                !< wall temperature
    real(R8), intent(out) :: omega(nsc)            !< species source term
    real(R8), intent(out) :: q_pyro                !< heat absorbed by the surface
    real(R8), intent(out) :: mdot                  !< blown mass flux
    real(R8), intent(in), optional :: y_wall(nsc)  !< pyrolysis product mass fractions

    ! Defined even if a model leaves one of them untouched.
    omega  = 0d0
    q_pyro = 0d0
    mdot   = 0d0

    select case ( gsi_kind )

      case ( GSI_KIND_PYROLYSIS )
        if ( .not. present(y_wall) ) &
          error stop 'GSI_Source: pyrolysis requires the product mass fractions y_wall'
        select case ( gsi_id )
          case ( GSI_PYRO_HTPB ) ; call GSI_HTPB ( T_wall, y_wall, omega, q_pyro, mdot )
          case ( GSI_PYRO_HDPE ) ; call GSI_HDPE ( T_wall, y_wall, omega, q_pyro, mdot )
          case ( GSI_PYRO_PP )   ; call GSI_PP   ( T_wall, y_wall, omega, q_pyro, mdot )
          case default
            error stop 'GSI_Source: unknown pyrolysis model id'
        end select

      case ( GSI_KIND_REACTIONS )
        select case ( gsi_id )
          case ( GSI_REAC_CARBON_BRADLEY ) ; call GSI_Carbon_Bradley ( roi_wall, T_wall, omega, q_pyro, mdot )
          case default
            error stop 'GSI_Source: unknown surface reaction model id'
        end select

      case default
        error stop 'GSI_Source: unknown gas-surface interaction kind'

    end select

    ! Orient the gas-surface interaction with respect to the boundary cell.
    !
    ! Everything the model returns is a magnitude, measured positive out of the
    ! solid, so all three quantities flip together with the face orientation.
    ! The energy balance solved in Solve_Wall_Temperature is
    !
    !     q_conv - q_pyro - q_rad = 0
    !
    ! in which q_conv carries the factor (-modfm2), through modfm3 = -2*modfm2,
    ! and -q_rad carries it too, through q_rad = qrad*modfm2. q_pyro must
    ! therefore carry the same factor, otherwise the balance is sign
    ! inconsistent on one of the two orientations and has no physical root
    ! there. Since (-modfm2) is +1 on a low face, this only ever changes the
    ! high faces (Fm = 2, 4, 6).
    mdot   = mdot   * (-face % modfm2)
    omega  = omega  * (-face % modfm2)
    q_pyro = q_pyro * (-face % modfm2)

  end subroutine GSI_Source


  subroutine Solve_Wall_Temperature ( face, gsi_kind, gsi_id, Prim, rho, Rgas, T_cell, q_rad, &
                                      T_wall, rho_wall, roi_wall, mu_wall, k_wall, grad_T, &
                                      q_conv, omega, q_pyro, mdot, y_wall )
    !< Find the wall temperature that closes the surface energy balance
    !< `q_conv - q_pyro - q_rad = 0`, by a secant iteration.
    !<
    !< Each iteration probes the balance at a slightly perturbed temperature to build
    !< a finite-difference slope, extrapolates to the root, and re-evaluates there so
    !< that every returned quantity is consistent with the wall temperature.
    implicit none
    type(wall_face_t), intent(in) :: face
    integer, intent(in) :: gsi_kind, gsi_id
    real(R8), intent(in)  :: Prim(nprim), rho, Rgas, T_cell
    real(R8), intent(in)  :: q_rad          !< radiative flux, signed as seen from the cell
    real(R8), intent(out) :: T_wall
    real(R8), intent(out) :: rho_wall, roi_wall(nsc), mu_wall, k_wall
    real(R8), intent(out) :: grad_T(3), q_conv
    real(R8), intent(out) :: omega(nsc), q_pyro, mdot
    real(R8), intent(in), optional :: y_wall(nsc)
    ! Local
    integer  :: iter
    real(R8) :: balance, balance_old, T_wall_old

    ! First guess: the wall sits at the boundary cell temperature, which makes the
    ! wall temperature gradient - and hence q_conv - vanish identically.
    T_wall = T_cell
    call Evaluate ( T_wall, balance )

    do iter = 1, MAX_T_WALL_ITER

      T_wall_old  = T_wall
      balance_old = balance

      ! Probe the balance slope around the current wall temperature.
      call Evaluate ( T_wall_old * ( 1d0 + T_WALL_PERTURBATION ), balance )

      if ( balance == balance_old ) then
        ! Flat balance: the secant step is undefined. Restore the last consistent
        ! state and give up rather than propagating a division by zero.
        T_wall = T_wall_old
        call Evaluate ( T_wall, balance )
        exit
      end if

      ! Secant update of the wall temperature.
      T_wall = T_wall_old - balance_old * ( T_WALL_PERTURBATION * T_wall_old ) / ( balance - balance_old )

      ! Re-evaluate at the new temperature, both for the next secant step and so that
      ! the state left behind by the loop matches the wall temperature it returns.
      call Evaluate ( T_wall, balance )

      ! TODO: enabling this early exit would cut the cost of the boundary condition
      !       by roughly an order of magnitude, but it changes the converged state
      !       and therefore the results. Left disabled on purpose.
      ! if ( abs(balance) < T_WALL_BALANCE_TOL ) exit

    end do

  contains

    subroutine Evaluate ( T_try, balance )
      !< Wall state, gas-surface interaction and energy balance at a trial temperature.
      implicit none
      real(R8), intent(in)  :: T_try
      real(R8), intent(out) :: balance

      call Wall_Gas_State ( face, Prim, rho, Rgas, T_cell, T_try, &
                            rho_wall, roi_wall, mu_wall, k_wall, grad_T, q_conv )

      call GSI_Source ( gsi_kind, gsi_id, face, roi_wall, T_try, omega, q_pyro, mdot, y_wall )

      balance = q_conv - q_pyro - q_rad

    end subroutine Evaluate

  end subroutine Solve_Wall_Temperature


  subroutine Ablation_Wall_Fluxes ( Im, Jm, Km, face, Blk, Prim, rho, T_wall, rho_wall, Prim_wall, &
                                    mu_wall, grad_T, q_gas, q_wall, omega, mdot, rans_wall_model, Ovar )
    !< Assemble the boundary fluxes of a blowing wall and subtract them from the residual.
    !<
    !< The injected mass leaves the surface at the blowing velocity `mdot/rho_wall`
    !< along the face normal; this drives both a momentum flux and the velocity
    !< gradient that feeds the viscous stress.
    implicit none
    integer, intent(in) :: Im, Jm, Km
    type(wall_face_t), intent(in) :: face
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), intent(in)    :: Prim(nprim)        !< boundary cell primitive variables
    real(R8), intent(in)    :: rho                !< boundary cell density
    real(R8), intent(in)    :: T_wall, rho_wall
    real(R8), intent(inout) :: Prim_wall(nprim)   !< wall variables; RANS slice set here
    real(R8), intent(in)    :: mu_wall
    real(R8), intent(in)    :: grad_T(3)          !< wall temperature gradient
    real(R8), intent(in)    :: q_gas              !< heat the gas loses through the wall
    real(R8), intent(in)    :: q_wall             !< total heat absorbed by the surface (reported)
    real(R8), intent(in)    :: omega(nsc)         !< species source term
    real(R8), intent(in)    :: mdot               !< blown mass flux
    integer,  intent(in)    :: rans_wall_model    !< one of `RANS_WALL_*`
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    ! Local
    integer  :: s
    real(R8) :: Gradient(nprim,3), Blowing(3), Stress(3), Flux(nprim), Sum_Omegai_Hwi, dl

    Gradient = 0d0
    Gradient(np,:) = grad_T

    ! Blowing velocity and the velocity gradient it induces normal to the face.
    Blowing = mdot / rho_wall * face % normal
    Gradient(nu:nw,face % dir) = ( Prim(nu:nw) - Blowing ) * face % modfm3
    Gradient(nu:nw,:) = matmul ( Gradient(nu:nw,:), face % M )

    if ( model > 0 ) then
      Stress = Stress_Vector ( Gradient(nu:nw,:), face % normal, mu_wall, 0d0, Prim(nt:) )
    else
      Stress = 0d0
    end if

    ! Enthalpy carried by the species injected at the wall.
    Sum_Omegai_Hwi = 0d0
    do s = 1, nsc
      Sum_Omegai_Hwi = Sum_Omegai_Hwi + omega(s) * f_tabT ( T_wall, s, h_tab )
    end do

    ! The boundary face carries the *complete* flux: the pressure the surface exerts
    ! on the gas, the momentum and enthalpy convected in with the blown mass, and the
    ! viscous terms. There is deliberately no separate inviscid wall flux alongside
    ! this one: a reflective (symmetry) treatment would impose u.n = 0 at the face and
    ! return the impermeable-wall pressure, which overshoots the true wall pressure by
    ! roughly rho*c*u_blow -- small against p, but comparable to the pressure
    ! differences that drive a low-speed blown boundary layer.
    !
    ! Under the boundary layer assumption used throughout, p(wall) = p(cell).
    Flux = 0d0
    Flux(1:nsc) = - omega * face % area
    Flux(nu:nw) = ( Stress - mdot * Blowing - Prim(np) * face % normal ) * face % area
    ! Energy leaving the gas: what it conducts into the wall, minus the enthalpy the
    ! injected species carry back in. Radiation absorbed by the solid does not appear
    ! here -- it reaches the surface from outside the domain and is spent driving the
    ! ablation, so charging it to the gas would cool the boundary cell for free.
    Flux(np)    = face % area * ( q_gas - Sum_Omegai_Hwi - mdot * 0.5d0 * sum ( Blowing**2 ) + &
                                  dot_product ( Stress, Blowing ) )

    if ( model == 2 ) then
      dl = Blk % yn(Im,Jm,Km)

      select case ( rans_wall_model )
        case ( RANS_WALL_NO_BLOWING )
          call RANS_Set_Wall_Values ( mu_wall, Prim(nt:nprim) * rho_wall / rho, Prim_wall(nt:nprim), dl, &
                                      k_rough=0d0 ) ! GSI walls are smooth
        case ( RANS_WALL_BLOWING )
          call RANS_Set_Blowing_Wall ( rho=rho_wall, mil=mu_wall, &
                                       rans_variables=Prim_wall(nt:nprim), &
                                       tau=Stress, mdot=mdot, dist=dl )
        case default
          error stop 'Ablation_Wall_Fluxes: unknown RANS wall model'
      end select

      Gradient(nt:nprim,face % dir) = ( Prim(nt:nprim)/rho - Prim_wall(nt:nprim)/rho_wall ) * face % modfm3
      Gradient(nt:nprim,:) = matmul ( Gradient(nt:nprim,:), face % M )

      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mu_wall, rho=rho_wall, area=face % area, &
                                 normal=face % normal, dist=WALL_DISTANCE )

      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_wall(nt:nprim), rho, rho_wall, &
                                   Blk % P(nt:nprim, face % ghost(1), face % ghost(2), face % ghost(3)) )
    end if

    ! Residual update
    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) - face % modfm2 * Flux

    if ( present(Ovar) ) &
      call Compute_Wall_Properties ( stress=Stress, pw=Prim(np), qw=q_wall, mdot=mdot, &
                                     y=Blk % dl(Im,Jm,Km) % c(face % dir) * 0.5d0, &
                                     Tw=T_wall, rhow=rho_wall, mu=mu_wall, exit_array=Ovar )

  end subroutine Ablation_Wall_Fluxes


  subroutine BC_Wall_Ablation ( Im, Jm, Km, Fm, Blk, gsi_kind, gsi_id, qrad, Ovar, y_wall, &
                               ablating, T_surface )
    !< Blowing wall whose temperature is set by the surface energy balance.
    !<
    !< Backs both `BC_Wall_Pyrolysis` and `BC_Wall_Reactions`, which differ only in
    !< the gas-surface interaction model they select.
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    integer, intent(in) :: gsi_kind                !< GSI family, one of `GSI_KIND_*`
    integer, intent(in) :: gsi_id                  !< model within the family
    real(R8), intent(in) :: qrad                   !< radiative heat flux at the wall
    type(MOSE_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    real(R8), intent(in), optional :: y_wall(nsc)  !< pyrolysis product mass fractions
    logical,  intent(out), optional :: ablating    !< .false. if the surface is inert
    real(R8), intent(out), optional :: T_surface   !< wall temperature, usable either way
    ! Local
    logical :: ablates
    type(wall_face_t) :: face
    real(R8) :: Prim(nprim), rho, Rgas, T_cell
    real(R8) :: T_wall, rho_wall, Prim_wall(nprim), mu_wall, k_wall
    real(R8) :: grad_T(3), q_conv, q_pyro, q_rad, omega(nsc), mdot

    call Get_Wall_Face ( Im, Jm, Km, Fm, Blk, face )

    ! Boundary cell state
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T_cell = Prim(np) / ( rho * Rgas )

    q_rad = qrad * face % modfm2

    call Solve_Wall_Temperature ( face, gsi_kind, gsi_id, Prim, rho, Rgas, T_cell, q_rad, &
                                  T_wall, rho_wall, Prim_wall(1:nsc), mu_wall, k_wall, grad_T, &
                                  q_conv, omega, q_pyro, mdot, y_wall )

    ! An inert surface gets no ablation flux at all: leave the residual untouched and
    ! report back, so the caller can close the face with an impermeable wall instead.
    ablates = Surface_Is_Ablating ( face, mdot, T_wall )
    if ( present(ablating) ) ablating = ablates
    if ( present(T_surface) ) then
      if ( ieee_is_finite(T_wall) .and. T_wall > 0d0 ) then
        T_surface = T_wall
      else
        ! The energy balance did not converge to anything usable; the boundary cell
        ! temperature at least gives the caller a zero-gradient (adiabatic) wall.
        T_surface = T_cell
      end if
    end if
    if ( .not. ablates ) return

    ! q_rad closes the wall temperature balance above; the gas only ever loses q_conv.
    call Ablation_Wall_Fluxes ( Im, Jm, Km, face, Blk, Prim, rho, T_wall, rho_wall, Prim_wall, &
                                mu_wall, grad_T, q_gas=q_conv, q_wall=q_conv, omega=omega, mdot=mdot, &
                                rans_wall_model=RANS_WALL_BLOWING, Ovar=Ovar )

  end subroutine BC_Wall_Ablation

end module MOSE_Lib_BC_Fluxes_GSI_Core
