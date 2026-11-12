!> Unit check of the preconditioning of integration-variables = prec (and of the chain rule of = prim).
!>
!> 1. The update of a cell (Compute_Prim_Residual, MOSE_Lib_Newstate) turns the residual of the conservative variables
!>    into the increment of the primitive ones. With prec it solves Gamma dQ = dU, Q = (p, u, v, w, T, Y_1..Y_N-1), dU
!>    the conservative increment, with the Weiss-Smith matrix Gamma: the Jacobian dU/dQ in which dp/drho at constant T,
!>    rho_p = 1/(R T), is replaced by Theta = 1/Ur^2 + 1/(cp T) (Weiss and Smith, AIAA J. 33(11), 1995). With prim
!>    Theta = rho_p, the plain chain rule. This program recovers dQ from the increment of the primitive variables that
!>    the routine returns (linearised equation of state), builds Gamma dQ row by row from the definition, and requires
!>    it to equal dU within 1e-10 of the sum of the magnitudes of the terms of each row: for one species of constant cp
!>    and for three species of different gas constants, heat capacities and heats of formation, at Ur = a, 0.1 a,
!>    1e-2 a, 1e-3 a and 1e-4 a, and with prim. At Ur = a the increment of prec is that of prim (Theta = rho_p),
!>    within 1e-10.
!> 2. The reference velocity (comp_Ur, MOSE_Lib_Preconditioning) in its three modes:
!>    - default (Weiss-Smith): eps_min a below eps_min a, f |v| between, a at and above a; then the viscous floor;
!>    - preconditioning-Uref: Uref where |v| <= Uref, a elsewhere; no viscous floor;
!>    - preconditioning-Mach: the Ur whose preconditioned Mach number u(1 - alpha)/sqrt(alpha^2 u^2 + Ur^2) equals the
!>      target (alpha = (1 - Ur^2/a^2)/2), within 1e-12; a where the flow is already faster than the target;
!>    and preconditioning-Ur-min raises Ur in every mode.
program prec_transform_check
  use iso_fortran_env, only: R8 => real64
  use MOSE_Global_m, only: nsc, nu, nv, nw, np, nt, nrans, nprim
  use MOSE_Config_Types_m, only: obj_prec
  use MOSE_Lib_Newstate, only: Compute_Prim_Residual
  use MOSE_Lib_Preconditioning, only: comp_Ur
  use FLINT_Lib_Thermodynamic, only: ns, Tmin, Tmax, wm_tab, Ri_tab, cp_tab, h_tab, s_tab, Runiv, &
                                     build_transposed_tables, f_Rtot, f_ss
  implicit none
  real(R8), parameter :: tol = 1.0e-10_R8
  integer :: nfail = 0

  call transform_cases(1)
  call transform_cases(3)
  call reference_velocity()

  if (nfail > 0) then
    write(*,'(A,I0,A)') ' PrecTransform: ', nfail, ' check(s) failed'
    stop 1
  end if
  write(*,'(A)') ' PrecTransform: all checks passed'

contains

  !> Species of constant cp: N2, O2, H2O (molar masses, cp in J/(kg K), heat of formation in J/kg)
  subroutine set_gas(n)
    integer, intent(in) :: n
    real(R8), parameter :: wm(3) = [28.0134_R8, 31.9988_R8, 18.0153_R8]
    real(R8), parameter :: cps(3) = [1040.0_R8, 918.0_R8, 1864.0_R8]
    real(R8), parameter :: hf(3) = [0.0_R8, 0.0_R8, -1.3423e7_R8]
    integer :: iT
    ns = n; nsc = n
    nu = nsc + 1; nv = nsc + 2; nw = nsc + 3; np = nsc + 4; nt = np + 1; nrans = 0; nprim = np
    if (allocated(wm_tab)) deallocate(wm_tab, Ri_tab, cp_tab, h_tab, s_tab)
    allocate(wm_tab(ns), Ri_tab(ns))
    wm_tab = wm(1:n)
    Ri_tab = Runiv / wm_tab
    Tmin = 0; Tmax = 5000
    allocate(cp_tab(Tmin:Tmax, ns), h_tab(Tmin:Tmax, ns), s_tab(Tmin:Tmax, ns))
    do iT = Tmin, Tmax
      cp_tab(iT, :) = cps(1:n)
      h_tab(iT, :) = hf(1:n) + cps(1:n) * iT
    end do
    s_tab = 0.0_R8
    call build_transposed_tables()
  end subroutine set_gas

  subroutine transform_cases(n)
    integer, intent(in) :: n
    real(R8), parameter :: ratios(5) = [1.0_R8, 0.1_R8, 1.0e-2_R8, 1.0e-3_R8, 1.0e-4_R8]
    real(R8) :: prim(n + 4), r(n + 4), dprim(n + 4), dprim_a(n + 4), dprim_prim(n + 4)
    real(R8) :: y(n), rho, a, p, T
    integer :: istate, iu
    character(len=40) :: tag

    call set_gas(n)
    do istate = 1, 3
      ! states: p, T, velocity and composition; residuals of mixed sign and scale
      select case (istate)
      case (1); p = 1.0e5_R8; T = 300.0_R8;  prim(nu:nw) = [10.0_R8, -3.0_R8, 2.0_R8]
      case (2); p = 7.0e5_R8; T = 1800.0_R8; prim(nu:nw) = [250.0_R8, 40.0_R8, -90.0_R8]
      case (3); p = 2.0e4_R8; T = 650.0_R8;  prim(nu:nw) = [-0.5_R8, 0.2_R8, 0.0_R8]
      end select
      if (n == 1) then
        y = 1.0_R8
      else
        y = [0.70_R8, 0.20_R8, 0.10_R8]
        if (istate == 2) y = [0.60_R8, 0.05_R8, 0.35_R8]
      end if
      rho = p / (T * sum(y * Ri_tab))
      prim(1:n) = rho * y
      prim(np) = p
      a = f_ss(prim(1:n), p, rho, f_Rtot(prim(1:n)))
      r(1:n) = [(1.0e-3_R8 * rho * (-1.0_R8)**iu * iu, iu = 1, n)]
      r(nu:nw) = rho * [3.0_R8, -1.5_R8, 0.7_R8]
      r(np) = -2.0e-2_R8 * p

      obj_prec%enabled = .false.
      dprim_prim = r
      call Compute_Prim_Residual(prim, dprim_prim, 1.0_R8, 1.0_R8, 0.0_R8, 1.0_R8)
      write(tag,'(A,I0,A,I0,A)') 'ns ', n, ', state ', istate, ', prim'
      call check_gamma(tag, prim, r, dprim_prim, 1.0_R8 / (p / rho))

      obj_prec%enabled = .true.
      do iu = 1, size(ratios)
        dprim = r
        call Compute_Prim_Residual(prim, dprim, 1.0_R8, 1.0_R8, ratios(iu) * a, 1.0_R8)
        write(tag,'(A,I0,A,I0,A,ES7.0,A)') 'ns ', n, ', state ', istate, ', Ur ', ratios(iu), ' a'
        call check_gamma(tag, prim, r, dprim, theta_of(ratios(iu) * a, prim))
        if (iu == 1) dprim_a = dprim
      end do
      if (any(abs(dprim_a - dprim_prim) > tol * (abs(dprim_a) + abs(dprim_prim)) + tiny(1.0_R8))) then
        write(*,'(A,I0,A,I0)') ' [FAIL] prec at Ur = a differs from prim: ns ', n, ', state ', istate
        nfail = nfail + 1
      else
        write(*,'(A,I0,A,I0)') ' [ok  ] prec at Ur = a gives the increment of prim: ns ', n, ', state ', istate
      end if
    end do
  end subroutine transform_cases

  !> Theta of Weiss and Smith: 1/Ur^2 - rho_T/(rho cp), rho_T = -rho/T
  function theta_of(Ur, prim) result(theta)
    real(R8), intent(in) :: Ur, prim(:)
    real(R8) :: theta
    theta = 1.0_R8 / Ur**2 + 1.0_R8 / (cp_mix(prim) * temperature(prim))
  end function theta_of

  function temperature(prim) result(T)
    real(R8), intent(in) :: prim(:)
    real(R8) :: T
    T = prim(np) / (sum(prim(1:nsc)) * sum(prim(1:nsc) * Ri_tab) / sum(prim(1:nsc)))
  end function temperature

  function cp_mix(prim) result(cp)
    real(R8), intent(in) :: prim(:)
    real(R8) :: cp
    cp = sum(prim(1:nsc) * cp_tab(0, :)) / sum(prim(1:nsc))
  end function cp_mix

  !> Recover dQ from the increment of the primitive variables and require Gamma(theta) dQ = -r (the residual is minus
  !> the conservative increment): mass, momentum, energy and the species 1..N-1
  subroutine check_gamma(tag, prim, r, dprim, theta)
    character(len=*), intent(in) :: tag
    real(R8), intent(in) :: prim(:), r(:), dprim(:), theta
    real(R8) :: rho, y(nsc), Rg, T, cp, h(nsc), hmix, u(3), q2, H0
    real(R8) :: drho, dy(nsc), dR, dT, dp, du(3), dmass, rows(nsc + 4), want(nsc + 4), size_(nsc + 4)
    integer  :: i
    logical  :: ok

    rho = sum(prim(1:nsc)); y = prim(1:nsc) / rho; Rg = sum(y * Ri_tab)
    T = prim(np) / (rho * Rg); cp = sum(y * cp_tab(0, :))
    h = h_tab(0, :) + cp_tab(0, :) * T
    hmix = sum(y * h); u = prim(nu:nw); q2 = sum(u * u); H0 = hmix + 0.5_R8 * q2

    ! dQ from the increment of the primitive variables (rho_i, u, v, w, p)
    drho = sum(dprim(1:nsc)); du = dprim(nu:nw); dp = dprim(np)
    dy = (dprim(1:nsc) - y * drho) / rho
    dR = sum(dy * Ri_tab)
    dT = T * (dp / prim(np) - drho / rho - dR / Rg)

    ! Gamma dQ: the mass row with theta in place of rho_p; d rho/dY_i at constant p, T is -rho R_i/R (the sum of dY is 0)
    dmass = theta * dp - rho / T * dT - rho / Rg * sum(Ri_tab * dy)
    rows(1) = dmass
    rows(2:4) = u * dmass + rho * du
    rows(5) = H0 * dmass + rho * (cp * dT + sum(h * dy) + sum(u * du)) - dp
    size_(1) = abs(theta * dp) + abs(rho / T * dT) + rho / Rg * sum(abs(Ri_tab * dy))
    size_(2:4) = abs(u) * size_(1) + rho * abs(du)
    size_(5) = abs(H0) * size_(1) + rho * (abs(cp * dT) + sum(abs(h * dy)) + sum(abs(u * du))) + abs(dp)
    want(1) = -sum(r(1:nsc)); want(2:4) = -r(nu:nw); want(5) = -r(np)
    do i = 1, nsc - 1
      rows(5 + i) = y(i) * dmass + rho * dy(i)
      size_(5 + i) = y(i) * size_(1) + rho * abs(dy(i))
      want(5 + i) = -r(i)
    end do

    ok = all(abs(rows(1:nsc + 4) - want(1:nsc + 4)) <= tol * size_(1:nsc + 4))
    if (ok) then
      write(*,'(A,A)') ' [ok  ] Gamma dQ = dU: ', trim(tag)
    else
      write(*,'(A,A)') ' [FAIL] Gamma dQ /= dU: ', trim(tag)
      write(*,'(A,*(ES11.3))') '        relative error of the rows: ', abs(rows(1:nsc + 4) - want(1:nsc + 4)) / size_(1:nsc + 4)
      nfail = nfail + 1
    end if
  end subroutine check_gamma

  subroutine reference_velocity()
    real(R8), parameter :: a = 340.0_R8
    real(R8) :: Ur, alpha, mach
    integer :: i
    real(R8), parameter :: speeds(4) = [1.0_R8, 50.0_R8, 200.0_R8, 500.0_R8]

    obj_prec%eps_min = 0.05_R8; obj_prec%Ur_factor = 1.0_R8; obj_prec%Ur_min = 0.0_R8
    obj_prec%Ur_ref = 0.0_R8; obj_prec%Mach_target = 0.0_R8
    ! default
    call expect('default, |v| < eps a', comp_Ur(1.0_R8, a, 0.0_R8), 0.05_R8 * a)
    call expect('default, eps a < |v| < a', comp_Ur(50.0_R8, a, 0.0_R8), 50.0_R8)
    call expect('default, |v| > a', comp_Ur(500.0_R8, a, 0.0_R8), a)
    call expect('default, viscous floor', comp_Ur(50.0_R8, a, 80.0_R8), 80.0_R8)
    obj_prec%Ur_factor = 2.0_R8
    call expect('default, Ur-factor 2', comp_Ur(50.0_R8, a, 0.0_R8), 100.0_R8)
    obj_prec%Ur_factor = 1.0_R8
    obj_prec%Ur_min = 30.0_R8
    call expect('default, Ur-min', comp_Ur(1.0_R8, a, 0.0_R8), 30.0_R8)
    obj_prec%Ur_min = 0.0_R8
    ! Uref
    obj_prec%Ur_ref = 20.0_R8
    call expect('Uref, |v| <= Uref', comp_Ur(1.0_R8, a, 0.0_R8), 20.0_R8)
    call expect('Uref, |v| > Uref', comp_Ur(50.0_R8, a, 0.0_R8), a)
    call expect('Uref, no viscous floor', comp_Ur(1.0_R8, a, 80.0_R8), 20.0_R8)
    obj_prec%Ur_min = 30.0_R8
    call expect('Uref, Ur-min', comp_Ur(1.0_R8, a, 0.0_R8), 30.0_R8)
    obj_prec%Ur_min = 0.0_R8; obj_prec%Ur_ref = 0.0_R8
    ! Mach target 0.6: the preconditioned Mach number of the flow is the target
    obj_prec%Mach_target = 0.6_R8
    do i = 1, 2
      Ur = comp_Ur(speeds(i + 1), a, 0.0_R8)
      alpha = 0.5_R8 * (1.0_R8 - Ur**2 / a**2)
      mach = speeds(i + 1) * (1.0_R8 - alpha) / sqrt(alpha**2 * speeds(i + 1)**2 + Ur**2)
      call expect('Mach target, preconditioned Mach number', mach, 0.6_R8)
    end do
    call expect('Mach target, flow faster than the target', comp_Ur(500.0_R8, a, 0.0_R8), a)
    obj_prec%Mach_target = 0.0_R8
  end subroutine reference_velocity

  subroutine expect(what, got, want)
    character(len=*), intent(in) :: what
    real(R8), intent(in) :: got, want
    if (abs(got - want) <= 1.0e-12_R8 * abs(want)) then
      write(*,'(A,A,A,ES14.7)') ' [ok  ] Ur: ', what, ': ', got
    else
      write(*,'(A,A,A,ES14.7,A,ES14.7)') ' [FAIL] Ur: ', what, ': ', got, ', expected ', want
      nfail = nfail + 1
    end if
  end subroutine expect

end program prec_transform_check
