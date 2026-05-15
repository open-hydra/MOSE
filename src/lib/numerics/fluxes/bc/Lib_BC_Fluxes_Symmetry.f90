module MOSE_Lib_BC_Fluxes_Symmetry
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm
  use MOSE_Lib_RotatingFrame

  implicit none
  public

contains

  subroutine BC_Symmetry_Eul ( Bm, Im, Jm, Km, Fm, Blk )
    use MOSE_Lib_Limiters, only : rlimiter
    use FLINT_Lib_Thermodynamic
    use MOSE_Mod_Riemann
    use MOSE_Config_Types_m, only: obj_rot
    implicit none
    integer, intent(in) :: Bm, Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: Blk
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k, dir, Face_i, Face_j, Face_k, s
    real(R8) :: Normal(3), Area, Bound_Dl, Int_Dl, Un, Un_wall, Velt(3), Dist12, Dist23, Dist
    real(R8), dimension(nprim) :: Bound_Prim, Int_Prim, Ghost_Prim, Prim1, Prim3
    real(R8), dimension(nprim) :: Slope32, Slope21, Slope, Face_Prim
    real(R8) :: rotot, Rtot, Sound, Un1, F_r, F_u, F_v, F_w, F_E, Flux(nprim), U1(3), U4(3)
    real(R8) :: r_fc(3)


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    ! (im,jm,km): boundary cell; (Int_i,Int_j,Int_k) intern cell for reconstruction
    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)
    Int_k = Km + guide(Fm,3)

    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    normal = Blk % dir(dir) % f(Face_i,Face_j,Face_k) % n
    area = Blk % dir(dir) % f(Face_i,Face_j,Face_k) % a
    Bound_Dl = 0.5d0 * Blk % dl(Im,Jm,Km) % c(dir)
    Int_Dl = 0.5d0 * Blk % dl(Int_i,Int_j,Int_k) % c(dir)

    Bound_Prim = Blk % P(:,Im,Jm,Km)
    Int_Prim = Blk % P(:,Int_i,Int_j,Int_k)

    ! Normal velocity of the wall face in the rotating frame.
    ! For a hub (rotating with the frame) or a symmetry plane: Un_wall = 0.
    ! For a stationary casing/shroud: Un_wall = -(Omega x r_face) . n, which is
    ! non-zero when the surface is not cylindrical around the rotation axis.
    Un_wall = 0.0_R8
    if ( obj_rot%enabled .AND. RF_Is_Stationary_Face(Bm, Fm) ) then
      call RF_Face_Center ( Blk%node, Im, Jm, Km, Fm, r_fc )
      Un_wall = dot_product ( RF_Wall_Velocity(r_fc), Normal )
    end if
    ! boundary cell: normal velocity to the boundary face
    Un = dot_product ( Bound_Prim (nu:nw), Normal )

    ! boundary cell: tangential component
    Velt = Bound_Prim (nu:nw) - Un * Normal

    ! velocity of the boundary cell with opposite normal component => ghost cell
    Ghost_Prim = Bound_Prim
    if ( Un_wall /= 0.0_R8 ) then
      Ghost_Prim (nu:nw) = Ghost_Prim (nu:nw) - 2d0*(Un - Un_wall)*Normal
    else
      Ghost_Prim (nu:nw) = Ghost_Prim (nu:nw) - 2d0 * Un * Normal
    end if
    ! stencil for reconstruction is centered on the boundary cell
    ! for faces 1,3,5 the ghost cell is selected as point (1) of the stencil. Otherwise, the intern cell
    Prim1 = modfm * Ghost_Prim + modfm1 * Int_Prim

    ! distance between point (1) and (2)
    Dist12 = 2.d0 * Bound_Dl * modfm + ( Bound_Dl + Int_Dl ) * modfm1

    ! for faces 1,3,5 the intern point is point (3) of the stencil. Otherwise, the ghost cell
    Prim3 = modfm * Int_Prim + modfm1 * Ghost_Prim

    ! distance between point (2) and (3)
    Dist23 = ( Bound_Dl + Int_Dl ) * modFm + 2d0 * Bound_Dl * modfm1

    ! reconstruction on the boundary face
    Slope32 = ( Prim3 - Bound_Prim ) / Dist23
    Slope21 = ( Bound_Prim - Prim1 ) / Dist12

    Slope = 0d0
    do s = 1, np
      Slope(s) = rlimiter ( Slope32(s), Slope21(s) )
    end do

    ! extrapolation at the boundary face
    Dist = Bound_Dl * modfm2
    Face_Prim = Bound_Prim + Slope * Dist

    ! Riemann problem at the interface
    call co_rotot_Rtot ( Face_Prim(1:nsc), rotot, Rtot )
    Sound = f_ss ( Face_Prim(1:nsc), Face_Prim(np), rotot, Rtot )
    Un  = dot_product ( Face_Prim(nu:nw), Normal )
    Velt = ( Face_Prim(nu:nw) - Un*Normal)
    if ( Un_wall /= 0.0_R8 ) then
      Un1 = (Un - Un_wall) * modFm2
      U1  =  Un_wall * Normal + Un1 * Normal + Velt
      U4  =  Un_wall * Normal - Un1 * Normal + Velt
    else
      Un1 = Un * modFm2
      U1  =  Un1 * Normal + Velt
      U4  = -Un1 * Normal + Velt
    end if
    call Riemann ( Face_Prim(1:nsc), U1(1), U1(2), U1(3), Face_Prim(np), Sound, rotot, &
                   Face_Prim(1:nsc), U4(1), U4(2), U4(3), Face_Prim(np), Sound, rotot, &
                   1d0, Normal(1), Normal(2), Normal(3), F_r, F_u, F_v, F_w, F_e )

    ! Fluxes
    Flux = 0d0
    Flux (nu) = F_u * Area
    Flux (nv) = F_v * Area
    Flux (nw) = F_w * Area

    ! Residual update
    Blk % r(nu:nw,Im,Jm,Km) = Blk % r(nu:nw,Im,Jm,Km) + Modfm2 * Flux(nu:nw)

  end subroutine BC_Symmetry_Eul


  subroutine BC_Symmetry_Visc ( Im, Jm, Km, Fm, Blk, RSM )
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_RANS
    use MOSE_Lib_Fluid
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    logical, intent(in)  :: RSM
    type(MOSE_block_type), intent(inout) :: Blk
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, dir, Face_i, Face_j, Face_k, a, b, c
    integer :: ig, jg, kg
    real(R8) :: Normal(3), Area, M(3,3), Prim(nprim), rho, Rgas, mil, mie
    real(R8) :: Prim_Wall(nprim), Gradrij(nRANS,3), dist
    real(R8) :: Gradient(3,3), Stress(3), Flux(nprim)
    real(R8), dimension(3) :: Vel1, Vel2, Vel3, Vel4, Vel5, Vel6, Vel7, Vel8, Vel9, Vel10, Vel


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    dist = Blk % yn(Im,Jm,Km)
    M = Blk % M (Im,Jm,Km) % c

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    mil = f_laminarViscosity ( Prim(1:nsc), Prim(np), rho, Rgas )

    Gradient = 0d0

    select case (Fm)
      case(1:2)
        Vel1  = Blk % P(nu:nw,im,jm,km)
        Vel3  = Blk % P(nu:nw,im,jm-1,km)
        Vel4  = Blk % P(nu:nw,im,jm+1,km)
        Vel7  = Blk % P(nu:nw,im,jm,km-1)
        Vel8  = Blk % P(nu:nw,im,jm,km+1)
        a = 1
        b = 2
        c = 3
    
      case(3:4)
        Vel1  = Blk % P(nu:nw,im,jm,km)
        Vel3  = Blk % P(nu:nw,im-1,jm,km)
        Vel4  = Blk % P(nu:nw,im+1,jm,km)
        Vel7  = Blk % P(nu:nw,im,jm,km-1)
        Vel8  = Blk % P(nu:nw,im,jm,km+1)
        a = 2
        b = 1
        c = 3
      
      case(5:6)
        Vel1  = Blk % P(nu:nw,im,jm,km)
        Vel3  = Blk % P(nu:nw,im-1,jm,km)
        Vel4  = Blk % P(nu:nw,im+1,jm,km)
        Vel7  = Blk % P(nu:nw,im,jm-1,km)
        Vel8  = Blk % P(nu:nw,im,jm+1,km)
        a = 3
        b = 1
        c = 2

    end select

    ! Symmetry condition
    Vel2  = Vel1 - 2d0 * dot_product ( Vel1, Normal ) * Normal
    Vel5  = Vel3 - 2d0 * dot_product ( Vel3, Normal ) * Normal
    Vel6  = Vel4 - 2d0 * dot_product ( Vel4, Normal ) * Normal
    Vel9  = Vel7 - 2d0 * dot_product ( Vel7, Normal ) * Normal
    Vel10 = Vel8 - 2d0 * dot_product ( Vel8, Normal ) * Normal

    Gradient (:,a) = ( Vel2 - Vel1 ) * Modfm2
    Gradient (:,b) = ( Vel4 - Vel3 + Vel6 - Vel5 ) * 0.25d0
    Gradient (:,c) = ( Vel8 - Vel7 + Vel10 - Vel9 ) * 0.25d0

    Gradient = matmul ( Gradient, M )

    ! Eddy viscosity
    mie = 0d0
    if (model==2) then
      call Eddy_Viscosity ( mut=mie, rans_variables=Prim(nt:nprim), &
                            mul=mil, rho=rho, vel_gradient=Gradient, &
                            walldist=dist )
    end if
    
    Stress = Stress_Vector ( Gradient, Normal, mil, mie, Prim(nt:) )

    Vel = 0.5d0 * ( Vel1 + Vel2 )

    ! Calcolo dei flussi allinterfaccia I+1/2
    Flux = 0d0
    Flux(nu:nw) = Area * Stress
    Flux(np) = Area * dot_product ( Vel, Stress )

    ! Special treatment of Reynolds stress for RSM models
    if ( RSM ) then
      Gradrij = 0d0 ! initialization
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3) ! ghost cell coordinates
      call RSM_Symmetry ( Prim(nt:), Prim_Wall(nt:), Blk%P(nt:,Ig,Jg,Kg), Normal ) ! setting Rij
      Prim_Wall(nprim) = Prim(nprim) 
      Blk%P(nprim,Ig,Jg,Kg) = Prim(nprim) ! setting omega
      Gradrij(:,Dir) = ( Prim(nt:nprim) - Prim_Wall(nt:nprim) )/rho * modfm3 ! gradient normal to BC face
      Gradrij = matmul ( Gradrij, M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient, &
                                 rans_gradient=Gradrij, &
                                 mul=mil, rho=rho, &
                                 area=area, normal=normal, dist=dist )
      Blk % r (nt:nprim,Im,Jm,Km) = Blk % r (nt:nprim,Im,Jm,Km) - Modfm2 * Flux(nt:nprim)
    endif

    ! Residual update
    Blk % r (nu:np,Im,Jm,Km) = Blk % r (nu:np,Im,Jm,Km) - Modfm2 * Flux(nu:np)
  
  end subroutine BC_Symmetry_Visc

end module MOSE_Lib_BC_Fluxes_Symmetry
