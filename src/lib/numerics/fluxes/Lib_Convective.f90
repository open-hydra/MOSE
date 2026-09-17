module MOSE_Lib_Convective
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Convective_Flux

contains

  subroutine Convective_Flux ( dl, normal, area, Prim, Res, beta, Ur_0, Ur_1 )
    use MOSE_Global_m
    use MOSE_Lib_Reconstruction, only: Reconstruction
    use MOSE_Mod_Riemann
    use FLINT_Lib_Thermodynamic
    implicit none
    real(R8), intent(in), dimension(-1:2) :: dl
    real(R8), intent(in) :: normal(3), area
    real(R8), intent(in), dimension(nprim,-1:2) :: Prim
    real(R8), intent(inout), dimension(nprim,0:1) :: Res
    real(R8), intent(in) :: beta
    real(R8), intent(in) :: Ur_0, Ur_1
    ! Local
    real(R8) :: l0, l1, l2, lm, lp
    real(R8), dimension(nprim) :: Prim1, Prim4
    real(R8) :: rotot1, Rtot1, rotot4, Rtot4, a1, a4
    real(R8) :: F_r, F_u, F_v, F_w, F_E, su, sp, sm, Flux(nprim)
    integer  :: error1, error4

    l0 = 0.5d0 * ( dl(-1) + dl(0) )
    l1 = 0.5d0 * ( dl(0) + dl(1) )
    l2 = 0.5d0 * ( dl(1) + dl(2) )
    lm = 0.5d0 * dl(0)
    lp = 0.5d0 * dl(1)
    
    ! Reconstruction phase. Stencil around interface /i (cells /i and /i+1): (i-1),(i),(i+1),(i+2)
    call Reconstruction ( Prim(:,-1), Prim(:,0), Prim(:,1), Prim(:,2), l0, l1, l2, lm, lp, beta, Prim1, Prim4 )

    ! Check for non-physical states at the interface and correct them if necessary
    call check_gas_state ( Prim1(1:nsc), Prim1(np), error1 )
    call check_gas_state ( Prim4(1:nsc), Prim4(np), error4 )
    if (error1 /= 0 .or. error4 /= 0) then
      write(*,*) '[ERROR] Non-physical state detected at the interface.'
      stop
    end if

    ! Compute rho and Rtot
    call co_rotot_Rtot ( Prim1(1:nsc), rotot1, Rtot1 )
    call co_rotot_Rtot ( Prim4(1:nsc), rotot4, Rtot4 )

    ! Compute speed of sound
    a1   = f_ss ( Prim1(1:nsc), Prim1(np), rotot1, Rtot1 )
    a4   = f_ss ( Prim4(1:nsc), Prim4(np), rotot4, Rtot4 )

    ! Solve the Riemann problem
    call Riemann (prim1(1:nsc), prim1(nu), prim1(nv), prim1(nw), prim1(np), a1, rotot1, &
                  prim4(1:nsc), prim4(nu), prim4(nv), prim4(nw), prim4(np), a4, rotot4, &
                  beta, Ur_0, Ur_1, normal(1), normal(2), normal(3), F_r, F_u, F_v, F_w, F_E)

    ! Sign of velocity at the interface
    su = sign ( 0.5d0, F_r )
    sp = 0.5d0 + su
    sm = su - 0.5d0

    ! Riemann fluxes are multiplied by the interface area
    Flux(1:nsc) = F_r * Area * ( sp*(prim1(1:nsc)/rotot1) - sm*(prim4(1:nsc)/rotot4) )
    Flux(nu)    = F_u * Area
    Flux(nv)    = F_v * Area
    Flux(nw)    = F_w * Area
    Flux(np)    = F_E * Area

    if (nprim>np) then
      Flux(np+1:nprim) = F_r * Area * ( sp * Prim1(np+1:nprim)/rotot1 - sm * Prim4(np+1:nprim)/rotot4 )
    end if

    ! Conservation-form residuals
    Res (:,0) = Res (:,0) + Flux
    Res (:,1) = Res (:,1) - Flux
    
  end subroutine Convective_Flux


  subroutine Sanitize_Primitive_State(state)
    use MOSE_Global_m, only: nprim, nsc, nu, nv, nw, np
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(R8), intent(inout) :: state(nprim)
    integer :: i
    real(R8), parameter :: state_cap = 1.0d100
    real(R8), parameter :: vel_cap = 1.0d100

    do i = 1, nsc
      if ((.not. ieee_is_finite(state(i))) .or. state(i) <= tiny(1.0_R8)) then
        state(i) = tiny(1.0_R8)
      end if
      state(i) = min(state(i), state_cap)
    end do

    if ((.not. ieee_is_finite(state(np))) .or. state(np) <= tiny(1.0_R8)) then
      state(np) = tiny(1.0_R8)
    end if
    state(np) = min(state(np), state_cap)

    if (.not. ieee_is_finite(state(nu))) state(nu) = 0.0_R8
    if (.not. ieee_is_finite(state(nv))) state(nv) = 0.0_R8
    if (.not. ieee_is_finite(state(nw))) state(nw) = 0.0_R8
    state(nu) = max(-vel_cap, min(state(nu), vel_cap))
    state(nv) = max(-vel_cap, min(state(nv), vel_cap))
    state(nw) = max(-vel_cap, min(state(nw), vel_cap))

    if (nprim > np) then
      do i = np + 1, nprim
        if (.not. ieee_is_finite(state(i))) state(i) = 0.0_R8
        state(i) = max(-state_cap, min(state(i), state_cap))
      end do
    end if
  end subroutine Sanitize_Primitive_State


end module MOSE_Lib_Convective