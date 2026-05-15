module MOSE_Lib_BC_Fluxes_Extrapolation
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index

  implicit none
  public

contains

  subroutine BC_Extrapolation ( Im, Jm, Km, Fm, Blk )
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(MOSE_block_type), intent(inout) :: Blk
    ! Local
    integer :: Dir, Face_i, Face_j, Face_k, FluxSign
    real(R8) :: Normal(3), Area, Prim(nprim), Un, Rgas, Fmass, Flux(nprim)


    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    
    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    ! boundary cell primitive variables
    Prim = Blk % P(:,Im,Jm,Km)
    
    ! boundary cell: normal velocity to the boundary face
    Un = dot_product ( Prim(nu:nw), Normal )

    Rgas = f_Rtot ( Prim(1:nsc) )
    Flux(1:nsc) = Prim(1:nsc) * Un * Area
    Fmass = sum ( Flux(1:nsc) )
    Flux(nu:nw) =  Fmass * Prim(nu:nw) + Prim(np) * Area * Normal
    Flux(np) = Fmass * H0 ( Prim(np), Prim(1:nsc), norm2 ( Prim(nu:nw) ) )
    if (nprim>np) Flux(np+1:nprim) = Prim(np+1:nprim) * Un * Area
    
    ! Residual update
    FluxSign = 1 - 2 * mod(Fm,2)
    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + FluxSign * Flux

  end subroutine BC_Extrapolation  

end module MOSE_Lib_BC_Fluxes_Extrapolation