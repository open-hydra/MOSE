module MOSE_Lib_BC_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  public

contains

  pure subroutine Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    implicit none
    integer, intent(in) :: fm
    integer, intent(out) :: modfm, modfm1, modfm2, modfm3

    modFm  = mod(Fm,2)
    modfm1 = 1-modFm  
    modfm2 = 1-2*modFm
    modfm3 = 2*( 2*modFm-1 )

  end subroutine Compute_Modfm


  pure subroutine Face_Index ( face, dir, cell_i, cell_j, cell_k, face_i, face_j, face_k )
    implicit none
    integer, intent(in) :: face, cell_i, cell_j, cell_k
    integer, intent(out) :: dir, face_i, face_j, face_k

    select case ( face )
      case(1:2)
        dir = 1
        face_i = cell_i - mod ( face, 2 )
        face_j = cell_j
        face_k = cell_k
      case(3:4)
        dir = 2
        face_i = cell_i
        face_j = cell_j - mod ( face, 2 )
        face_k = cell_k
      case(5:6)
        dir = 3
        face_i = cell_i
        face_j = cell_j
        face_k = cell_k - mod ( face, 2 )
    end select

  end subroutine Face_Index


  !> The walls: walls coupled to another solver (103, 104), prescribed heat flux (301), prescribed temperature (302),
  !> gas-surface interaction (503-505).
  elemental logical function Is_Wall ( type )
    implicit none
    integer, intent(in) :: type

    select case ( type )
      case ( 103, 104, 301, 302, 503, 504, 505 )
        Is_Wall = .true.
      case default
        Is_Wall = .false.
    end select

  end function Is_Wall


  !> Gas at the face of a wall of prescribed temperature T_wall, as BC_Wall_Temperature takes it: the pressure of the
  !> boundary cell (boundary layer) and its composition (non-catalytic wall) at the density rho_wall = p / (R T_wall), and
  !> the laminar viscosity and conductivity at T_wall. Prim, rho and Rgas: primitive variables, density and gas constant of
  !> the boundary cell; rhos_wall: the partial densities at the wall.
  subroutine Wall_State_Temperature ( Prim, rho, Rgas, T_wall, rho_wall, rhos_wall, mil, kl )
    use MOSE_Global_m, only: nsc, np
    use FLINT_Lib_Thermodynamic, only: co_k_mi_lam_Wilke
    implicit none
    real(R8), intent(in)  :: Prim(:), rho, Rgas, T_wall
    real(R8), intent(out) :: rho_wall, rhos_wall(nsc), mil, kl

    rho_wall = Prim(np) / ( Rgas * T_wall )
    rhos_wall = Prim(1:nsc) * rho_wall / rho
    call co_k_mi_lam_Wilke ( rhos_wall, rho_wall, T_wall, mil, kl )

  end subroutine Wall_State_Temperature


  !> Gas at the face of a wall of prescribed heat flux, as BC_Wall_Heat takes it: the laminar viscosity and conductivity at
  !> the temperature temp of the boundary cell; the difference dT = T_cell - T_wall across the half cell, from the heat flux
  !> with that conductivity and the metric Mc of the cell along the face direction Dir; the wall temperature Twall = temp - dT
  !> and the density rho_wall = p / (R Twall) (boundary layer, non-catalytic wall).
  subroutine Wall_State_Heat ( Prim, rho, Rgas, Heat_Flux, Mc, normal, Dir, temp, mil, kl, dT, Twall, rho_wall )
    use MOSE_Global_m, only: nsc, np
    use FLINT_Lib_Thermodynamic, only: co_k_mi_lam_Wilke, EOS
    implicit none
    real(R8), intent(in)  :: Prim(:), rho, Rgas, Heat_Flux, Mc(3,3), normal(3)
    integer,  intent(in)  :: Dir
    real(R8), intent(out) :: temp, mil, kl, dT, Twall, rho_wall

    temp = EOS( p=prim(np), rho=rho, R=Rgas )
    call co_k_mi_lam_Wilke ( Prim(1:nsc), rho, temp, mil, kl )
    ! Heat_Flux is positive into the fluid on every face: the wall is hotter than the cell when it heats the fluid
    dT = -Heat_Flux / ( 2d0*kl*dot_product( Mc(Dir,:), normal ) )
    Twall = temp - dT
    rho_wall = Prim(np)/( Rgas * Twall ) ! EOS

  end subroutine Wall_State_Heat


  subroutine Compute_Wall_Properties(stress,pw,Tw,rhow,mu,qw,mdot,y,exit_array)
    implicit none
    real(R8), dimension(:), intent(in)    :: stress
    real(R8), intent(in)                  :: pw
    real(R8), dimension(8), intent(inout) :: exit_array
    real(R8), optional :: Tw, qw, rhow, mu, y, mdot
    ! Local
    real(R8) :: tauw, ut

    ! y+
    tauw = sqrt( sum ( Stress**2 ) )
    ut = sqrt ( tauw/rhow )
    exit_array(1) = rhow*y*ut/mu
    ! Wall stress
    exit_array(2:4) = stress
    ! Wall pressure
    exit_array(5) = pw
    ! Wall temperature
    if (present(Tw)) exit_array(6) = Tw
    ! Heat flux
    if (present(qw)) exit_array(7) = qw
    ! Mass flux
    if (present(mdot)) exit_array(8) = mdot

  end subroutine Compute_Wall_Properties

end module MOSE_Lib_BC_Fluxes