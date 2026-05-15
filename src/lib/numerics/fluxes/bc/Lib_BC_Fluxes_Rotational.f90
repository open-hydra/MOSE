module MOSE_Lib_BC_Fluxes_Rotational
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  public

contains

  subroutine BC_Rotational_Symmetry_Eul ( Im, Jm, Km, Fm, Blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: Blk
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, dir, Face_i, Face_j, Face_k
    real(R8) :: Normal(3), Area
    real(R8) :: Prim(nprim), Flux(nprim)

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    normal = Blk % dir(dir) % f(Face_i,Face_j,Face_k) % n
    area = Blk % dir(dir) % f(Face_i,Face_j,Face_k) % a
    Prim = Blk % P(:,Im,Jm,Km)
    Flux(nv) = Prim(np) * Area * Normal(2)
    ! Residual update
    Blk % r(nv,Im,Jm,Km) = Blk % r(nv,Im,Jm,Km) + modFm2 * Flux(nv)

  end subroutine BC_Rotational_Symmetry_Eul


  subroutine BC_Rotational_Periodic_Eul ( Im, Jm, Km, Fm, Blk )
    use MOSE_Lib_Reconstruction, only : reconstruction
    use MOSE_Lib_Metrics, only : delthe
    use MOSE_Mod_Riemann
    use FLINT_Lib_Thermodynamic, only : f_ss, co_rotot_Rtot

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: Blk
    ! Local
    integer :: Dir, Face_i, Face_j, Face_k, FluxSign
    integer :: Fs, Ks, K1, K2, K3, K4, modFm, Modfm1, ModFm2, ModFm3, s
    real(R8) :: Normal(3), Area, Prim(nprim), rho, Face_V, Face_W, Fmass, Flux(nprim)
    real(R8), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4, diff21, diff32, diff43, Slope
    real(R8), dimension(nprim) :: Prim_L, Prim_R
    real(R8) :: v_rot, w_rot, Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34
    real(R8) :: rho_L, rho_R, Rtot_L, Rtot_R, a_L, a_R, F_r, F_u, F_v, F_w, F_e
    real(R8) :: su, sel_L, sel_R


    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    ! boundary cell primitive variables
    Prim = Blk % P(:,Im,Jm,Km)
    rho = sum ( Prim (1:nsc) )
    if (ndir /= 3) then
      ! NB: for face 5 and 6: nx = 0 ; ny = +- sin(delthe/2) ; nz = cos(delthe/2).
      ! Get velocity at the interface by rotating it of an angle = delthe/2 in the y-z plane since the K=1 cell 
      ! spans from -delthe/2 to delthe/2. The component of (uu, vv, ww) on the normal to the face is w because of the rotation.
      ! For reference, the formula are: uf = u ; vf =  v*nz + w*ny ; wf = -v*ny + w*nz
      Face_V =   Prim(nv) * Normal(3) + Prim(nw) * Normal(2)
      Face_W = - Prim(nv) * Normal(2) + Prim(nw) * Normal(3)

      ! Except from f(nsc+2) and f(nsc+3), fluxes from face 5 and 6 are equal and opposite. No need to compute them
      ! (nor the transport of scalar quantities such as roi and mit). Only terms needed: pf*ny (axisymmetry) and fmass*w*ny 
      ! (centrifugal force) in Y direction, -fmass*v*ny (?) in Z direction.
      Fmass  = rho * Area * Prim(nw)     ! written like this because (uu*nx+vv*ny+ww*nz) = w
      Flux = 0d0
      Flux(nv) = Face_V * Fmass + Prim(np) * Area * Normal(2)
      Flux(nw) = Face_W * Fmass

      ! Residual update
      FluxSign = 1 - 2 * mod(Fm,2)
      Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + FluxSign * Flux
    else
      ! Same logic as standard connection
      if (Fm == 5) then
        Fs = 6
        Ks = blk % dim(3)
      elseif (Fm == 6) then
        Fs = 5
        Ks = 1
      else
        write(*,*) '[ERROR] Rotational periodic boundary condition is only for face 5-6'
        stop
      end if

      K1 = modFm*( Ks + guide(Fs,3) ) + modFm1*( Km + guide(Fm,3) )
      K2 = modFm1*Km + modFm*Ks
      K3 = modFm1*Ks + modFm*Km
      K4 = modFm1*(Ks + guide(Fs,3)) + modFm*(Km + guide(Fm,3))

      ! Primitive variables for the stencil
      Prim1 = blk % P(:,Im,Jm,K1)
      Prim2 = blk % P(:,Im,Jm,K2)
      Prim3 = blk % P(:,Im,Jm,K3)
      Prim4 = blk % P(:,Im,Jm,K4)

      ! Rotation
      if (Fm == 5) then
        w_rot =  Prim1(nw) * cos(delthe) - Prim1(nv) * sin(delthe) 
        v_rot =  Prim1(nw) * sin(delthe) + Prim1(nv) * cos(delthe)
        Prim1(nv) = v_rot
        Prim1(nw) = w_rot
        w_rot =  Prim2(nw) * cos(delthe) - Prim2(nv) * sin(delthe) 
        v_rot =  Prim2(nw) * sin(delthe) + Prim2(nv) * cos(delthe)
        Prim2(nv) = v_rot
        Prim2(nw) = w_rot
      elseif (Fm == 6) then
        w_rot =  Prim3(nw) * cos(delthe) + Prim3(nv) * sin(delthe) 
        v_rot = -Prim3(nw) * sin(delthe) + Prim3(nv) * cos(delthe)
        Prim3(nv) = v_rot
        Prim3(nw) = w_rot
        w_rot =  Prim4(nw) * cos(delthe) + Prim4(nv) * sin(delthe) 
        v_rot = -Prim4(nw) * sin(delthe) + Prim4(nv) * cos(delthe)
        Prim4(nv) = v_rot
        Prim4(nw) = w_rot
      end if  

      Dl1 = blk % Dl(Im,Jm,K1) % c(3) * 0.5d0
      Dl2 = blk % Dl(Im,Jm,K2) % c(3) * 0.5d0
      Dl3 = blk % Dl(Im,Jm,K3) % c(3) * 0.5d0
      Dl4 = blk % Dl(Im,Jm,K4) % c(3) * 0.5d0

      Dist23 = Dl2 + Dl3
      Dist12 = Dl1 + Dl2
      Dist34 = Dl3 + Dl4

      ! ! Reconstruction in point (2)
      ! diff32 = ( Prim3 - Prim2 ) / Dist23
      ! diff21 = ( Prim2 - Prim1 ) / Dist12

      ! ! Limiter for point (2) slope => state (1) for the riemann solver
      ! do s = 1, nprim
      !   Slope(s) = rlimiter ( diff32(s), diff21(s) )
      ! end do
    
      ! ! Extrapolation from point (2) to the interface => state (1)
      ! Prim_L = Prim2 + Slope * Dl2

      ! ! Reconstruction in point (3)
      ! diff43 = ( Prim4 - Prim3 ) / Dist34

      ! ! Limiter for point (3) slope => state (4) for the riemann solver
      ! do s = 1, nprim
      !   Slope(s) = rlimiter ( diff43(s), diff32(s) )
      ! end do

      ! ! Extrapolation from point (3) to the interface => state (4)
      ! Prim_R = Prim3 - Slope * Dl3

      call Reconstruction ( Prim1, Prim2, Prim3, Prim4, Dist12, Dist23, Dist34, Dl2, Dl3, 1d0, Prim_L, Prim_R )

      ! Auxiliary variables for state (L)
      call co_rotot_Rtot ( Prim_L(1:nsc), rho_L, Rtot_L )
      a_L = f_ss ( Prim_L(1:nsc), Prim_L(np), rho_L, Rtot_L )

      ! Auxiliary variables for state (R)
      call co_rotot_Rtot ( Prim_R(1:nsc), rho_R, Rtot_R)
      a_R = f_ss ( Prim_R(1:nsc), Prim_R(np), rho_R, Rtot_R )

      call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(np), a_L, rho_L, &
                     Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(np), a_R, rho_R, &
                     1d0, normal(1),normal(2),normal(3), F_r, F_u, F_v, F_w, F_E)

      su = sign ( 0.5d0, F_r )
      Sel_L = 0.5d0 + su
      Sel_R = su - 0.5d0

      ! Fluxes
      Fmass = F_r * area
      Flux(1:nsc) = Fmass*( Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R )
      Flux(nu) = F_u * Area
      Flux(nv) = F_v * Area
      Flux(nw) = F_w * Area
      Flux(np) = F_E * Area

      if (nprim>np) then
        Flux(np+1:nprim) = Fmass*(sel_L*Prim_L(np+1:nprim)/rho_L - sel_R*Prim_R(np+1:nprim)/rho_R)
      end if

      ! Residual update
      if (Fm==6) then
        blk % r(:,Im,Jm,Km) = blk % r(:,Im,Jm,Km) + Flux
        ! Flux is in the face-6 frame (θ=delthe). Rotate momentum back to face-5
        ! frame (θ=0) before applying to Ks, so that (v,w) point in the right
        ! Cartesian directions for a cell sitting at θ=0.
        Flux(nv) = ( F_v * cos(delthe) + F_w * sin(delthe) ) * Area
        Flux(nw) = (-F_v * sin(delthe) + F_w * cos(delthe) ) * Area
        blk % r(:,Im,Jm,Ks) = blk % r(:,Im,Jm,Ks) - Flux
      end if
    end if

  end subroutine BC_Rotational_Periodic_Eul

end module MOSE_Lib_BC_Fluxes_Rotational