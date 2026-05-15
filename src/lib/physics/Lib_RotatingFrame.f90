!> @brief Single Rotating Frame (SRF) source terms for compressible Navier-Stokes.
!>
!> When activated, the full computational domain is assumed to rotate rigidly at
!> constant angular velocity omega [rad/s] about a fixed axis.  The velocity
!> stored in P(nu:nw) is the RELATIVE (rotating-frame) velocity
!>
!>   w = u_abs - Omega x r
!>
!> so that no changes are required to the convective or diffusive flux modules.
!>
!> --------------------------------------------------------------------------
!> GOVERNING EQUATIONS (source terms added to the RHS)
!> --------------------------------------------------------------------------
!>
!>   Continuity / species:  unchanged
!>
!>   Momentum (per unit volume):
!>     S_mom = -2 rho (Omega x w)         [Coriolis]
!>           + rho omega^2 r_perp          [centrifugal = -rho Omega x (Omega x r)]
!>     where  r_perp = r - (r.e) e  is the position vector perpendicular to the
!>     axis,  r = x_cell - x_origin.
!>
!>   Energy (per unit volume, E = cv T + 0.5 |w|^2):
!>     S_E = w . (rho omega^2 r_perp)
!>     [Coriolis does no work:  w . (-2 rho Omega x w) = 0 identically]
!>
!>   RANS equations: unchanged.
!>     For Spalart-Allmaras, activate the Spalart-Shur correction
!>     (turbulence-model = 'SA-RC') to remove spurious production due to
!>     solid-body rotation in the vorticity magnitude.
!>     SST / Wilcox k-omega use the strain-rate for production, which is
!>     already frame-invariant; no correction is needed.
!>
!> --------------------------------------------------------------------------
!> VISCOUS FLUXES
!> --------------------------------------------------------------------------
!>   tau = mu ( grad(w) + grad(w)^T - 2/3 div(w) I )
!>   Viscous energy flux = tau . w  ,  heat conduction q = -k grad(T)
!>   All quantities use the relative velocity w stored in P(nu:nw).
!>   Mod_Diffusive.f90 is ALREADY CORRECT; no changes required.
!>
!> --------------------------------------------------------------------------
!> BOUNDARY CONDITIONS
!> --------------------------------------------------------------------------
!>
!>   Inflow / outflow (type 4):
!>     T0, p0, alpha, beta and mach must be given as RELATIVE total conditions
!>     (rotating frame).  For inflow from a stationary region upstream:
!>
!>       w_inlet  = u_abs_inlet - Omega x r_inlet        (relative velocity)
!>       T0_rel   = T_static + |w_inlet|^2 / (2 cp)
!>       p0_rel   = p_static (T0_rel/T_static)^(gamma/(gamma-1))
!>       mach_rel = |w_inlet| / sqrt(gamma R T_static)
!>
!>     Outflow static-pressure (pamb) is frame-independent; no change needed.
!>
!>   Walls (type 5, 6, 8-14):
!>     No-slip condition w = 0 applies to the relative velocity.  The wall is
!>     stationary in the rotating frame, so the existing BC is correct.
!>
!>   Rotational periodicity (type 1 connection with pitch angle, or type 2):
!>     The relative velocity transforms between the two periodic faces by the
!>     same rotation matrix R(Dtheta) as the position vector.  For a sector with
!>     pitch angle Dtheta about the z-axis:
!>
!>       u_ghost = u_source cos(Dtheta) - v_source sin(Dtheta)
!>       v_ghost = u_source sin(Dtheta) + v_source cos(Dtheta)
!>       w_ghost = w_source   (axial)
!>
!>     This transformation must be embedded in the block connectivity (Mg tensor)
!>     when the mesh is generated.  The connection BC (type 1) then handles it
!>     automatically.  No code changes are needed provided the grid file encodes
!>     the correct rotational mapping between the two periodic faces.
!>
!> --------------------------------------------------------------------------
!> INITIALISATION (called from Read_Ini.f90 after loading input.ini)
!> --------------------------------------------------------------------------
!>   The input.ini section [MOSE-rotating-frame] is read via FiNeR.
!>   The raw fields of obj_rot are filled, then axis is normalised,
!>   Omega_vec computed, and the global flag obj_rot%enabled is set.
!>
!>   [MOSE-rotating-frame]
!>     active = .true.
!>     omega  = 1000.0        ! angular speed [rad/s]; sign follows right-hand rule
!>     axis   = 0.0 0.0 1.0   ! rotation-axis direction (need not be a unit vector)
!>     origin = 0.0 0.0 0.0   ! coordinates of one point on the axis [m]

!>   n-stationary-face  = 4
!>   stationary-face-1 = 1 4
!>   stationary-face-2 = 2 4
!>   stationary-face-3 = 3 4
!>   stationary-face-4 = 4 4
!>
module MOSE_Lib_RotatingFrame
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Config_Types_m, only : obj_rot
  implicit none

contains


  ! ---------------------------------------------------------------------------
  !> @brief Finalise the rotating-frame configuration after registry values
  !>        have been loaded by Load_Ini.
  !>
  !> Actions:
  !>   1. Parse axis_str / origin_str into real(3) arrays.
  !>   2. Parse stationary_face_str(:) into stationary_face(2,N).
  !>   3. Normalise obj_rot%axis to a unit vector.
  !>   4. Compute obj_rot%Omega_vec = omega * axis.
  !>   5. Set the global flag obj_rot%enabled.
  !>   6. Print a summary to stdout (root only).
  !>
  !> Called from Read_Inifile after Load_Ini has populated the registry.
  !> No FiNeR dependency: all values come from the registry.

  subroutine Setup_RotatingFrame()
  
    implicit none
    ! Local
    integer  :: k
    real(R8) :: norm

    if ( obj_rot%model == "" ) then
      obj_rot%enabled = .false.
      return
    elseif ( obj_rot%model == "rigid-body" ) then
      if ( obj_rot%omega < 0.0_R8 ) then
        write(*,'(A)') '[ERROR] Rotational frame model "rigid-body" requires omega >= 0'
        stop
      end if
      obj_rot%enabled = .true.
    end if

    ! Parse axis and origin from registry strings
    read(obj_rot%axis_str, *)   obj_rot%axis
    read(obj_rot%origin_str, *) obj_rot%origin

    ! Normalise axis
    norm = sqrt( sum(obj_rot%axis**2) )
    if ( norm < 1.0e-14_R8 ) then
      write(*,'(A)') '[ERROR] Rotation axis has zero length'
      stop
    end if
    obj_rot%axis      = obj_rot%axis / norm
    obj_rot%Omega_vec = obj_rot%omega * obj_rot%axis

    ! Parse stationary wall patches from registry strings
    if ( obj_rot%n_stationary > 0 ) then
      allocate( obj_rot%stationary_face(2, obj_rot%n_stationary) )
      do k = 1, obj_rot%n_stationary
        read(obj_rot%stationary_face_str(k), *) &
          obj_rot%stationary_face(1,k), obj_rot%stationary_face(2,k)
      end do
    end if

  end subroutine Setup_RotatingFrame


  ! ---------------------------------------------------------------------------
  !> @brief Add Coriolis and centrifugal source terms to the momentum and
  !>        energy residuals for all blocks in the domain.
  subroutine RotatingFrame_Source_Terms ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      call RF_Source_Blk( domain%blk(b)%P,    &
                          domain%blk(b)%R,    &
                          domain%blk(b)%node, &
                          domain%blk(b)%vol,  &
                          domain%blk(b)%dim   )
    end do

  end subroutine RotatingFrame_Source_Terms


  ! ---------------------------------------------------------------------------
  !> @brief Block-level kernel: loops over interior cells and accumulates the
  !>        rotating-frame source terms into the residual array R.
  !>
  !> Residual sign convention (same as throughout MOSE):
  !>   V * dU/dt = -R(U)   with  R = CONV - DIFF - SOURCE
  !> => to add a physical source S, we do  R -= S * vol
  subroutine RF_Source_Blk ( Prim, Res, node, vol, n )
    use MOSE_Base_Types_m
    use MOSE_Global_m

    implicit none
    integer,  intent(in) :: n(3)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    type(MOSE_vector_3D_type), dimension(0:n(1), 0:n(2), 0:n(3)),  intent(in) :: node
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
    ! Local
    integer  :: i, j, k
    real(R8) :: rho, w(3)
    real(R8) :: xc, yc, zc, r_vec(3), r_proj, r_perp(3)
    real(R8) :: Cor(3), Cen(3)
    ! Cache module-level data in local scalars to avoid repeated indirect lookups
    real(R8) :: Omega(3), ax(3), ox(3), omega2

    Omega  = obj_rot%Omega_vec
    ax     = obj_rot%axis
    ox     = obj_rot%origin
    omega2 = obj_rot%omega**2

    !$omp do collapse(3) &
    !$omp private( i, j, k, rho, w, xc, yc, zc, r_vec, r_proj, r_perp, Cor, Cen )
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      rho = sum( Prim(1:nsc, i, j, k) )
      w   = Prim(nu:nw, i, j, k)            ! relative velocity in rotating frame

      ! ---- Cell centre (average of 8 hexahedral corner nodes) ----
      ! Node (i,j,k) in [0..ni,0..nj,0..nk]:
      ! cell (i,j,k) occupies the hexahedron with corners at i-1..i, j-1..j, k-1..k.
      xc = 0.125_R8 * ( node(i-1,j-1,k-1)%c(1) + node(i,j-1,k-1)%c(1) &
                       + node(i-1,j,  k-1)%c(1) + node(i,j,  k-1)%c(1) &
                       + node(i-1,j-1,k  )%c(1) + node(i,j-1,k  )%c(1) &
                       + node(i-1,j,  k  )%c(1) + node(i,j,  k  )%c(1) )
      yc = 0.125_R8 * ( node(i-1,j-1,k-1)%c(2) + node(i,j-1,k-1)%c(2) &
                       + node(i-1,j,  k-1)%c(2) + node(i,j,  k-1)%c(2) &
                       + node(i-1,j-1,k  )%c(2) + node(i,j-1,k  )%c(2) &
                       + node(i-1,j,  k  )%c(2) + node(i,j,  k  )%c(2) )
      zc = 0.125_R8 * ( node(i-1,j-1,k-1)%c(3) + node(i,j-1,k-1)%c(3) &
                       + node(i-1,j,  k-1)%c(3) + node(i,j,  k-1)%c(3) &
                       + node(i-1,j-1,k  )%c(3) + node(i,j-1,k  )%c(3) &
                       + node(i-1,j,  k  )%c(3) + node(i,j,  k  )%c(3) )

      ! ---- Position vector from axis origin to cell centre ----
      r_vec = [ xc - ox(1), yc - ox(2), zc - ox(3) ]

      ! ---- Component of r perpendicular to the rotation axis ----
      !   r_perp = r_vec - (r_vec . axis) * axis
      r_proj = dot_product( r_vec, ax )
      r_perp = r_vec - r_proj * ax

      ! ---- Coriolis force per unit volume:  -2 rho (Omega x w) ----
      Cor(1) = -2.0_R8 * rho * ( Omega(2)*w(3) - Omega(3)*w(2) )
      Cor(2) = -2.0_R8 * rho * ( Omega(3)*w(1) - Omega(1)*w(3) )
      Cor(3) = -2.0_R8 * rho * ( Omega(1)*w(2) - Omega(2)*w(1) )

      ! ---- Centrifugal force per unit volume:  rho omega^2 r_perp ----
      !   (= -rho Omega x (Omega x r)  via the vector triple-product identity)
      Cen = rho * omega2 * r_perp

      ! ---- Accumulate into residuals:  R -= source * vol ----
      !   Momentum equations
      Res(nu, i, j, k) = Res(nu, i, j, k) - ( Cor(1) + Cen(1) ) * vol(i,j,k)
      Res(nv, i, j, k) = Res(nv, i, j, k) - ( Cor(2) + Cen(2) ) * vol(i,j,k)
      Res(nw, i, j, k) = Res(nw, i, j, k) - ( Cor(3) + Cen(3) ) * vol(i,j,k)

      !   Energy equation:  only centrifugal does work on w
      !   w . Coriolis = w . (-2 rho Omega x w) = 0  identically

      Res(np, i, j, k) = Res(np, i, j, k) - dot_product( w, Cen ) * vol(i,j,k)

    end do ; end do ; end do
    !$omp end do

  end subroutine RF_Source_Blk


  ! ---------------------------------------------------------------------------
  !> @brief Return the velocity of a lab-stationary wall in the rotating frame.
  !>
  !> A point fixed in the lab frame moves at  v_frame = Omega x (r - origin)
  !> in the lab frame.  In the rotating frame, that point has velocity
  !> w_wall = -v_frame.
  !>
  !> @param[in]  r_face   Cartesian position of the wall face centre [m]
  !> @return     w_wall   Relative velocity of the stationary wall [m/s]
  function RF_Wall_Velocity(r_face) result(w_wall)
    implicit none
    real(R8), intent(in) :: r_face(3)
    real(R8) :: w_wall(3)
    ! Local
    real(R8) :: r(3)   ! r_face - origin

    r = r_face - obj_rot%origin
    ! w_wall = -(Omega_vec x r)
    w_wall(1) = -( obj_rot%Omega_vec(2)*r(3) - obj_rot%Omega_vec(3)*r(2) )
    w_wall(2) = -( obj_rot%Omega_vec(3)*r(1) - obj_rot%Omega_vec(1)*r(3) )
    w_wall(3) = -( obj_rot%Omega_vec(1)*r(2) - obj_rot%Omega_vec(2)*r(1) )

  end function RF_Wall_Velocity


  ! ---------------------------------------------------------------------------
  !> @brief Compute the centre of a block face as the average of its 4 corner nodes.
  !>

  subroutine RF_Face_Center(node, Im, Jm, Km, Fm, r_face)
    use MOSE_Base_Types_m, only: MOSE_vector_3D_type

    implicit none
    type(MOSE_vector_3D_type), intent(in) :: node(0:,0:,0:)
    integer,  intent(in)  :: Im, Jm, Km, Fm
    real(R8), intent(out) :: r_face(3)
    ! Local
    integer :: Fi, Fj, Fk, d

    ! Face node index (same logic as Face_Index in Lib_BC.f90)
    select case (Fm)
    case(1,2)
      Fi = Im - mod(Fm,2) ; Fj = Jm ; Fk = Km
    case(3,4)
      Fi = Im ; Fj = Jm - mod(Fm,2) ; Fk = Km
    case(5,6)
      Fi = Im ; Fj = Jm ; Fk = Km - mod(Fm,2)
    end select

    ! 4-node average for each Cartesian component
    do d = 1, 3
      select case (Fm)
      case(1,2)   ! i-face: nodes span (Fi, Fj-1:Fj, Fk-1:Fk)
        r_face(d) = 0.25_R8 * ( node(Fi,Fj-1,Fk-1)%c(d) + node(Fi,Fj,Fk-1)%c(d) &
                               + node(Fi,Fj-1,Fk  )%c(d) + node(Fi,Fj,Fk  )%c(d) )
      case(3,4)   ! j-face: nodes span (Fi-1:Fi, Fj, Fk-1:Fk)
        r_face(d) = 0.25_R8 * ( node(Fi-1,Fj,Fk-1)%c(d) + node(Fi,Fj,Fk-1)%c(d) &
                               + node(Fi-1,Fj,Fk  )%c(d) + node(Fi,Fj,Fk  )%c(d) )
      case(5,6)   ! k-face: nodes span (Fi-1:Fi, Fj-1:Fj, Fk)
        r_face(d) = 0.25_R8 * ( node(Fi-1,Fj-1,Fk)%c(d) + node(Fi,Fj-1,Fk)%c(d) &
                               + node(Fi-1,Fj  ,Fk)%c(d) + node(Fi,Fj  ,Fk)%c(d) )
      end select
    end do

  end subroutine RF_Face_Center


  ! ---------------------------------------------------------------------------
  !> @brief Convert inflow BC conditions from the absolute frame to the
  !>        rotating (relative) frame.
  !>
  !> The operator specifies conditions in the lab frame.  This routine
  !> converts them to equivalent conditions in the rotating frame so that
  !> BC_Subsonic_Inflow / BC_Supersonic_Inflow receive relative quantities.
  !>
  !> @param[in]     r_face         Face centre position [m]
  !> @param[in]     BC_Mach_abs    Absolute Mach number (>0)
  !> @param[in]     gamma          Ratio of specific heats at inlet
  !> @param[in]     Rgas           Specific gas constant at inlet [J/(kg K)]
  !> @param[inout]  T0_or_Tstat    In: abs total T or T_static; Out: rel total T or T_static
  !> @param[inout]  p0_or_pstat    In: abs total p or p_static; Out: rel total p or p_static
  !> @param[inout]  alpha          In: abs flow angle (yz-plane); Out: rel flow angle [rad]
  !> @param[inout]  beta           In: abs flow angle (elevation); Out: rel flow angle [rad]
  !> @param[out]    Mach_rel       Relative Mach number
  subroutine RF_Convert_BC_Inflow(Blk, Fm, Im, Jm, Km, BC_Mach_abs, gamma, Rgas, &
                                   T0_or_Tstat, p0_or_pstat, alpha, beta, Mach_rel)
    use MOSE_Advanced_Types_m
    implicit none
    type(MOSE_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Fm, Im, Jm, Km
    real(R8), intent(in)    :: BC_Mach_abs, gamma, Rgas
    real(R8), intent(inout) :: T0_or_Tstat, p0_or_pstat, alpha, beta
    real(R8), intent(out)   :: Mach_rel
    ! Local
    integer  :: Face_i, Face_j, Face_k
    real(R8) :: v_frame(3), r(3), r_fc(3)
    real(R8) :: T_stat, p_stat, a_stat, u_abs(3), w_rel(3), w_norm, cp

    ! Face node index (same convention as Face_Index in Lib_BC.f90)
    select case (Fm)
    case(1,2) ; Face_i = Im - mod(Fm,2) ; Face_j = Jm             ; Face_k = Km
    case(3,4) ; Face_i = Im             ; Face_j = Jm - mod(Fm,2) ; Face_k = Km
    case(5,6) ; Face_i = Im             ; Face_j = Jm             ; Face_k = Km - mod(Fm,2)
    end select

    ! Face centre: 4-node average of the face corner nodes
    select case ( Fm )
        case(1,2)   ! i-face: nodes at (Face_i, Face_j-1:Face_j, Face_k-1:Face_k)
          r_fc(1) = 0.25d0*( Blk%node(Face_i,Face_j-1,Face_k-1)%c(1) + Blk%node(Face_i,Face_j,Face_k-1)%c(1) &
                            + Blk%node(Face_i,Face_j-1,Face_k  )%c(1) + Blk%node(Face_i,Face_j,Face_k  )%c(1) )
          r_fc(2) = 0.25d0*( Blk%node(Face_i,Face_j-1,Face_k-1)%c(2) + Blk%node(Face_i,Face_j,Face_k-1)%c(2) &
                            + Blk%node(Face_i,Face_j-1,Face_k  )%c(2) + Blk%node(Face_i,Face_j,Face_k  )%c(2) )
          r_fc(3) = 0.25d0*( Blk%node(Face_i,Face_j-1,Face_k-1)%c(3) + Blk%node(Face_i,Face_j,Face_k-1)%c(3) &
                            + Blk%node(Face_i,Face_j-1,Face_k  )%c(3) + Blk%node(Face_i,Face_j,Face_k  )%c(3) )
        case(3,4)   ! j-face: nodes at (Face_i-1:Face_i, Face_j, Face_k-1:Face_k)
          r_fc(1) = 0.25d0*( Blk%node(Face_i-1,Face_j,Face_k-1)%c(1) + Blk%node(Face_i,Face_j,Face_k-1)%c(1) &
                            + Blk%node(Face_i-1,Face_j,Face_k  )%c(1) + Blk%node(Face_i,Face_j,Face_k  )%c(1) )
          r_fc(2) = 0.25d0*( Blk%node(Face_i-1,Face_j,Face_k-1)%c(2) + Blk%node(Face_i,Face_j,Face_k-1)%c(2) &
                            + Blk%node(Face_i-1,Face_j,Face_k  )%c(2) + Blk%node(Face_i,Face_j,Face_k  )%c(2) )
          r_fc(3) = 0.25d0*( Blk%node(Face_i-1,Face_j,Face_k-1)%c(3) + Blk%node(Face_i,Face_j,Face_k-1)%c(3) &
                            + Blk%node(Face_i-1,Face_j,Face_k  )%c(3) + Blk%node(Face_i,Face_j,Face_k  )%c(3) )
        case(5,6)   ! k-face: nodes at (Face_i-1:Face_i, Face_j-1:Face_j, Face_k)
          r_fc(1) = 0.25d0*( Blk%node(Face_i-1,Face_j-1,Face_k)%c(1) + Blk%node(Face_i,Face_j-1,Face_k)%c(1) &
                            + Blk%node(Face_i-1,Face_j,  Face_k)%c(1) + Blk%node(Face_i,Face_j,  Face_k)%c(1) )
          r_fc(2) = 0.25d0*( Blk%node(Face_i-1,Face_j-1,Face_k)%c(2) + Blk%node(Face_i,Face_j-1,Face_k)%c(2) &
                            + Blk%node(Face_i-1,Face_j,  Face_k)%c(2) + Blk%node(Face_i,Face_j,  Face_k)%c(2) )
          r_fc(3) = 0.25d0*( Blk%node(Face_i-1,Face_j-1,Face_k)%c(3) + Blk%node(Face_i,Face_j-1,Face_k)%c(3) &
                            + Blk%node(Face_i-1,Face_j,  Face_k)%c(3) + Blk%node(Face_i,Face_j,  Face_k)%c(3) )
    end select

    ! Frame velocity at face centre: v_frame = Omega_vec x (r_fc - origin)
    r = r_fc - obj_rot%origin
    v_frame(1) = obj_rot%Omega_vec(2)*r(3) - obj_rot%Omega_vec(3)*r(2)
    v_frame(2) = obj_rot%Omega_vec(3)*r(1) - obj_rot%Omega_vec(1)*r(3)
    v_frame(3) = obj_rot%Omega_vec(1)*r(2) - obj_rot%Omega_vec(2)*r(1)

    ! Static conditions (frame-invariant)
    if ( BC_Mach_abs > 1.0_R8 ) then
      ! Supersonic: input T0_or_Tstat = T_static, p0_or_pstat = p_static
      T_stat = T0_or_Tstat
      p_stat = p0_or_pstat
    else
      ! Subsonic: input T0_or_Tstat = T0_abs, p0_or_pstat = p0_abs
      T_stat = T0_or_Tstat / ( 1.0_R8 + 0.5_R8*(gamma-1.0_R8)*BC_Mach_abs**2 )
      p_stat = p0_or_pstat * ( T_stat / T0_or_Tstat )**( gamma/(gamma-1.0_R8) )
    end if

    ! Speed of sound (same in both frames)
    a_stat = sqrt( gamma * Rgas * T_stat )

    ! Absolute velocity vector from Mach and flow angles
    !   u = M*a * [cos(alpha)*cos(beta), sin(alpha)*cos(beta), sin(beta)]
    u_abs(1) = BC_Mach_abs * a_stat * cos(alpha) * cos(beta)
    u_abs(2) = BC_Mach_abs * a_stat * sin(alpha) * cos(beta)
    u_abs(3) = BC_Mach_abs * a_stat * sin(beta)

    ! Relative velocity: w = u_abs - v_frame
    w_rel  = u_abs - v_frame
    w_norm = norm2( w_rel )

    ! Relative Mach number and new flow angles
    Mach_rel = w_norm / a_stat
    alpha    = atan2( w_rel(2), w_rel(1) )
    beta     = atan2( w_rel(3), sqrt(w_rel(1)**2 + w_rel(2)**2) )

    if ( BC_Mach_abs > 1.0_R8 ) then
      ! T_static and p_static unchanged (frame-invariant)
      ! Mach_rel and angles already updated above
    else
      ! Relative total conditions from relative velocity
      cp          = gamma * Rgas / ( gamma - 1.0_R8 )
      T0_or_Tstat = T_stat + 0.5_R8 * w_norm**2 / cp
      p0_or_pstat = p_stat * ( T0_or_Tstat / T_stat )**( gamma/(gamma-1.0_R8) )
    end if

  end subroutine RF_Convert_BC_Inflow


  ! ---------------------------------------------------------------------------
  !> @brief Return .true. if the wall patch (block b, face direction f) is
  !>        stationary in the lab frame (casing/hub).
  !>
  !> Called from Mod_BC.f90 for type-5/6 wall BCs when obj_rot%enabled is active.
  !> Patches listed in obj_rot%stationary_face receive w_wall = -(Omega x r);
  !> all other patches (blades) use the standard w_wall = 0.
  logical function RF_Is_Stationary_Face(b, f)
    implicit none
    integer, intent(in) :: b, f
    integer :: k

    RF_Is_Stationary_Face = .false.
    if ( .not. allocated(obj_rot%stationary_face) ) return
    do k = 1, size(obj_rot%stationary_face, 2)
      if ( obj_rot%stationary_face(1,k) == b .and. &
           obj_rot%stationary_face(2,k) == f ) then
        RF_Is_Stationary_Face = .true.
        return
      end if
    end do

  end function RF_Is_Stationary_Face


  ! ---------------------------------------------------------------------------
  !> @brief Estimate the pressure torque on the rotating wall surfaces.
  !>
  !> Loops over every BC entry that is a wall face (types 5,6,8-10,12,13) and
  !> is NOT listed in obj_rot%stationary_face (i.e., the rotating blade/helix
  !> surface).  The pressure force on each face is
  !>
  !>   F_wall = p * (outward-from-fluid unit normal) * area
  !>
  !> and the axial torque contribution is
  !>
  !>   dT = axis . ( (r_face - origin) x F_wall )
  !>
  !> Only pressure contributions are included; viscous shear is neglected
  !> (adequate for a first estimate on a turbine/propeller profile).
  !>
  !> Results are written to OUTPUT/torque.dat at each call, at the same
  !> frequency as the residual history (controlled by the caller).
  subroutine RF_Torque_Walls(domain, iter, time)
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: np
    use MOSE_Lib_Ghost, only: Is_Wall
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root, mpi_reduce_sum_r8

    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer,  intent(in) :: iter
    real(R8), intent(in) :: time
    ! Local
    integer,  save :: unit_torque = -1
    logical,  save :: first_call  = .true.
    integer  :: i, Bm, Im, Jm, Km, Fm, Dir, Fi, Fj, Fk
    real(R8) :: r_face(3), r_rel(3), F_wall(3), torque_vec(3)
    real(R8) :: T_axial, p_face, sign_fm
    real(R8) :: N_face(3), A_face

    T_axial = 0.0_R8

    do i = 1, domain%nbound

      ! Skip non-wall BC types
      if ( .not. Is_Wall(domain%bc(i)%type) ) cycle

      Bm = domain%bc(i)%b

      ! Skip blocks not owned by this rank
      if (.not. is_local_block(Bm)) cycle

      Fm = domain%bc(i)%f

      ! Skip faces that are stationary in the lab frame (casing / hub shroud)
      if ( RF_Is_Stationary_Face(Bm, Fm) ) cycle

      Im = domain%bc(i)%i
      Jm = domain%bc(i)%j
      Km = domain%bc(i)%k

      ! Direction index: 1 for i-faces (Fm=1,2), 2 for j-faces (3,4), 3 for k-faces (5,6)
      Dir = (Fm + 1) / 2

      ! Index of the boundary face in the precomputed metrics arrays
      select case (Fm)
      case (1,2) ; Fi = Im - mod(Fm,2) ; Fj = Jm             ; Fk = Km
      case (3,4) ; Fi = Im             ; Fj = Jm - mod(Fm,2) ; Fk = Km
      case (5,6) ; Fi = Im             ; Fj = Jm             ; Fk = Km - mod(Fm,2)
      end select

      ! Unit normal and area from precomputed face metrics.
      ! dir(d)%f%N always points in the direction of increasing index d.
      N_face = domain%blk(Bm)%dir(Dir)%f(Fi,Fj,Fk)%N
      A_face = domain%blk(Bm)%dir(Dir)%f(Fi,Fj,Fk)%A

      ! Outward-from-fluid sign:  -1 for min faces (Fm odd), +1 for max faces (Fm even)
      sign_fm = real(1 - 2 * mod(Fm, 2), R8)

      ! Pressure at the interior cell adjacent to the wall (zero-order extrapolation)
      p_face = domain%blk(Bm)%P(np, Im, Jm, Km)

      ! Face centre position (4-node average of corner nodes)
      call RF_Face_Center(domain%blk(Bm)%node, Im, Jm, Km, Fm, r_face)

      ! Pressure force on the wall from the fluid:
      !   F_wall = p * (outward-from-fluid normal) * A = p * sign_fm * N * A
      F_wall = p_face * sign_fm * N_face * A_face

      ! Torque contribution about the rotation axis origin:
      !   dT_axial = axis . ( (r_face - origin) x F_wall )
      r_rel        = r_face - obj_rot%origin
      torque_vec(1) = r_rel(2)*F_wall(3) - r_rel(3)*F_wall(2)
      torque_vec(2) = r_rel(3)*F_wall(1) - r_rel(1)*F_wall(3)
      torque_vec(3) = r_rel(1)*F_wall(2) - r_rel(2)*F_wall(1)
      T_axial = T_axial + dot_product(obj_rot%axis, torque_vec)

    end do

    ! Reduce partial torque sums from all ranks to root (blocking collective)
    call mpi_reduce_sum_r8(T_axial)

    ! Only root writes the output
    if (mpi_is_root) then
      if (first_call) then
        open(newunit=unit_torque, file='OUTPUT/torque.dat', &
             status='replace', form='formatted')
        write(unit_torque, '(A)') &
          '#        iter              time   torque_axial (pressure, non-dim)' // &
          '         power (T*omega, non-dim)'
        first_call = .false.
      end if

      write(unit_torque, '(I9,3E20.10)') iter, time, T_axial, T_axial * obj_rot%omega
      flush(unit_torque)
    end if

  end subroutine RF_Torque_Walls


end module MOSE_Lib_RotatingFrame
