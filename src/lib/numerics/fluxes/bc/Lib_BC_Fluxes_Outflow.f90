module MOSE_Lib_BC_Fluxes_Outflow
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none

contains

 subroutine BC_Outflow ( Bm, Im, Jm, Km, Fm, Blk, BC_pexit, BC_rel_fac, BC_mdot, error )
    implicit none
    type(MOSE_block_type), intent(inout) :: Blk
    integer,  intent(in) :: Bm, Im, Jm, Km, Fm
    real(R8), intent(in) :: BC_rel_fac, BC_pexit
    real(R8), intent(out) :: BC_mdot
    integer :: error
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k, Dir, Face_i, Face_j, Face_k
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim), Normal(3), Area, BC_Sign, Un
    real(R8) :: Flux(nprim)

    error = 0

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    ! (im,jm,km): boundary cell; (Int_i,Int_j,Int_k) next cell
    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)
    Int_k = Km + guide(Fm,3)

    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Boundary cell primitive variables
    Bound_Prim = Blk % P(:,Im,Jm,Km)
    Int_Prim = Blk % P(:,Int_i,Int_j,Int_k)

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    BC_Sign = Real ( modfm2 )  ! =-1 for faces 1,3,5 ; =1 for faces 2,4,6
    Un  = dot_product ( Bound_Prim(nu:nw), Normal )  ! velocity normal to the interface

    if (BC_Sign * Un < 0d0) then
      error = 1
      return
    endif

    call Outflow( Bound_Prim, BC_Sign, Un, Normal, Area, BC_pexit, Flux, BC_mdot )

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + Modfm2 * Flux
  
  end subroutine BC_Outflow


  subroutine Outflow( Bound_Prim, BC_Sign, Un, Normal, Area, BC_pexit, Flux, BC_mdot )
    use MOSE_Mod_Riemann
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in)  :: Bound_Prim(nprim), BC_Sign, Un, Normal(3), Area, BC_pexit
    real(R8), intent(out) :: Flux(nprim), BC_mdot
    ! Local
    real(R8) :: Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8) :: Vel_N(3), Vel_T(3), Mach_N, Face_Sound, Isentropic
    real(R8) :: Face_Un, Face_Prim(nprim), Face_rho, Face_Rgas, Face_Enthalpy, Face_T
    real(R8) :: Mach_abs, pstar, den, delta, p_exit
    integer :: s

    call co_rotot_Rtot ( Bound_Prim(1:nsc), Bound_rho, Bound_Rgas )
    Bound_Sound = f_ss( Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas )
    Bound_Gamma = f_gamma ( Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas )

    ! riemann invariant a+-delta*u directed towards the interface
    Riem = Bound_Sound + BC_Sign * 0.5d0*( Bound_Gamma - 1d0 ) * Un

    ! For outflow, interface values (subscript /f) are computed for the fluxes at the end
    Vel_N = Un * Normal                             ! normal velocity vector
    Vel_T = Bound_Prim(nu:nw) - Vel_N               ! tangential velocity vector
    Mach_abs = norm2(Bound_Prim(nu:nw))/Bound_Sound
    Mach_N = abs ( Un / Bound_Sound )               ! Mach number (NB normal to BC face)
  
    if ( ( Mach_N < 1d0 ) .and. ( BC_pexit > 0d0 ) ) then
      ! Subsonic outflow
      p_exit = BC_pexit
      delta = 0.5d0 * (Bound_Gamma - 1.d0)
      den = 1.d0/(1 + delta)**(Bound_Gamma/(Bound_Gamma-1.d0))
      pstar = Bound_Prim(np) * den * (1+delta * Mach_N**2)**(Bound_Gamma/(Bound_Gamma-1.d0))
      if (pstar >= p_exit) p_exit = pstar 
      Isentropic = ( Bound_Gamma - 1d0 ) / ( 2d0 * Bound_Gamma )
      ! isentropic expansion to ambient pressure
      Face_Sound = Bound_Sound * ( ( p_exit/Bound_Prim(np) )** Isentropic )
      ! normal velocity at the interface using riemann invariant
      Face_Un  = ( Riem - Face_Sound ) / ( BC_Sign*0.5d0*(Bound_Gamma - 1d0) )                
      
      Face_Prim(np) = p_exit
      Face_rho = Bound_Gamma * Face_Prim(np) / Face_Sound**2
      Face_Prim(1:nsc) = Bound_Prim(1:nsc) / Bound_rho * Face_rho ! mass fractions extrapolated
      Face_Prim(nu:nw) = Vel_T + Face_Un * Normal ! tangential velocity extrapolated
      if (nprim>np) then 
        Face_Prim(np+1:nprim) = Bound_Prim(np+1:nprim) / Bound_rho * Face_rho ! turbolence quantities (Q/rho) extrapolated
      end if
      Face_Rgas = Bound_Rgas
    
    else
      ! supersonic exit or pambient=0: extrapolation of all variables
      Face_Prim = Bound_Prim
      Face_Un = Un
      Face_Sound = Bound_Sound
      Face_rho = Bound_rho
      Face_Rgas = Bound_Rgas

    end if

    Face_Enthalpy = 0d0
    Face_T = Face_Prim(np) / ( Face_Rgas * Face_rho )
    Flux = 0.0d0
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT( Face_T,s,h_tab )
      Flux(s) = Face_Prim(s) * Area * Face_Un
    end do
    BC_mdot = Face_rho * Area * Face_Un
    Flux(nu:nw) = Face_Prim(nu:nw) * BC_mdot + Face_Prim(np) * Area * Normal
    Flux(np) = BC_mdot * ( 0.5d0 * sum( Face_Prim(nu:nw)**2 ) + Face_Enthalpy )
    if (nprim>np) Flux(np+1:nprim) = Face_Prim(np+1:nprim) / Face_rho * BC_mdot

  end subroutine Outflow

end module MOSE_Lib_BC_Fluxes_Outflow