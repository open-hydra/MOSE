!> Unit check of the face pressure of a reflective (symmetry) treatment, the reason Ablation_Wall_Fluxes
!> (MOSE_Lib_BC_Fluxes_GSI_Core) gives for the complete flux of a blowing face, with the pressure of the cell on it:
!> a reflective treatment would return the impermeable-wall pressure, below the pressure of the cell since the gas
!> moves away from the face, by about rho*c*u_blow in the exact reflective Riemann problem.
!>
!> BC_Symmetry_Eul hands the solver of riemann-solver (the pointer Riemann) the state of the face and its mirror image,
!> the normal velocity reversed, with beta = 1. This program does the same on N2 at 1e5 Pa and 300 K moving along the
!> normal at 10 m/s (the moving case of Q2D's SRM-grain-state: the low speed of a blown gas) and at 100 m/s, with the
!> real FLINT functions on one species of
!> constant cp (tables filled from the definition, gamma = 1.4), for the exact solver, the subject of the comment, and
!> for HLLC, the default of riemann-solver:
!>   1. gas moving away from the face (u > 0, the blown gas of a grain, 502, or of an ablating wall, 503-505): the face
!>      pressure below that of the cell; gas moving towards the face (u < 0): above it;
!>   2. the exact solver on the exact solution of the mirror problem, two rarefactions,
!>      p* = p (1 - (gamma-1)/2 u/c)^(2 gamma/(gamma-1)), here p - 0.983 rho c u, within 1e-5 p: the stopping rule of its
!>      Newton iteration, |1 - p2/p3| < 1e-5, leaves the face pressure (one of p2, p3) within 1e-5 of p*, which lies
!>      between them. 1e-5 p is 1 Pa, against 67 Pa between p* and the acoustic estimate p - rho c u at 10 m/s. At
!>      10 m/s the wave is so weak that the relations of a rarefaction and of a shock give the same p* within 6e-7 p
!>      (they agree to the third order in its strength: the left wave of the exact solver taken as a shock moved the
!>      face pressure by 0.06 Pa); at 100 m/s (M = 0.28) they differ by far more than 1e-5 p.
!> The face pressure is read from the momentum flux F_u (normal along x); the momentum of the mass flux, F_r u, must stay
!> below 1e-6 of the pressure difference checked (the mirror pair carries no mass flux).
!> Not checked here: the low-Mach and preconditioned solvers, which take a cutoff Mach number (riemann-options-Mco, which
!> the set-up requires above 0 for AUSM+M, HLLC+Chen and LMRoe) or the preconditioning velocity of the cell (Blk%Ur, the
!> url, urr of the call) that a run sets and this program does not.
program reflective_pressure_check
  use iso_fortran_env, only: R8 => real64
  use MOSE_Global_m, only: nsc
  use MOSE_Config_Types_m, only: obj_riemann
  use MOSE_Mod_Riemann, only: Riemann, Assign_Riemann_Solver
  use FLINT_Lib_Thermodynamic, only: ns, Tmin, Tmax, wm_tab, Ri_tab, cp_tab, h_tab, s_tab, Runiv, &
                                     build_transposed_tables, f_Rtot, f_ss
  implicit none
  integer, parameter :: nsol = 2
  character(len=8), parameter :: solvers(nsol) = [character(len=8) :: 'exact', 'HLLC']
  real(R8), parameter :: gam = 1.4_R8, p = 1.0e5_R8, T = 300.0_R8, speeds(2) = [10.0_R8, 100.0_R8]
  real(R8) :: R, cp, rho, a, d(1), u, speed, Fr, Fu, Fv, Fw, Fe, pf, pstar
  integer :: i, j, k, iT
  logical :: ok, ok_case

  ! One species, N2, of constant cp: gamma = cp / (cp - R) = 1.4 at every temperature
  ns = 1; nsc = 1
  allocate(wm_tab(ns), Ri_tab(ns))
  wm_tab = 28.0134_R8
  Ri_tab = Runiv / wm_tab
  R = Ri_tab(1)
  cp = gam / (gam - 1.0_R8) * R
  Tmin = 0; Tmax = 5000
  allocate(cp_tab(Tmin:Tmax, ns), h_tab(Tmin:Tmax, ns), s_tab(Tmin:Tmax, ns))
  do iT = Tmin, Tmax
    cp_tab(iT, 1) = cp
    h_tab(iT, 1) = cp * iT
  end do
  s_tab = 0.0_R8
  call build_transposed_tables()

  rho = p / (R * T)
  d = rho
  a = f_ss(d, p, rho, f_Rtot(d))
  write(*,'(A,F8.3,A)') ' N2 at 1e5 Pa, 300 K, c = ', a, ' m/s'

  obj_riemann%Mco = 0.0_R8
  ok = .true.
  do j = 1, size(speeds)
  speed = speeds(j)
  pstar = p * (1.0_R8 - 0.5_R8 * (gam - 1.0_R8) * speed / a)**(2.0_R8 * gam / (gam - 1.0_R8))
  write(*,'(A,F6.1,A,ES12.5,A,ES12.5,A)') ' speed ', speed, ' m/s: rho c u = ', rho * a * speed, &
    ' Pa, exact p* - p (away from the face) = ', pstar - p, ' Pa'
  do i = 1, nsol
    obj_riemann%description = solvers(i)
    call Assign_Riemann_Solver()
    do k = 1, 2
      u = merge(speed, -speed, k == 1)            ! the normal velocity of the cell, + away from the face
      ! the face along x, the mirror image on its left, the cell on its right (BC_Symmetry_Eul: U1 = -U4 along n)
      call Riemann(d, -u, 0.0_R8, 0.0_R8, p, a, rho, d, u, 0.0_R8, 0.0_R8, p, a, rho, 1.0_R8, 0.0_R8, 0.0_R8, &
                   1.0_R8, 0.0_R8, 0.0_R8, Fr, Fu, Fv, Fw, Fe)
      pf = Fu
      ok_case = merge(pf < p, pf > p, u > 0.0_R8) .and. abs(Fr * speed) <= 1.0e-6_R8 * abs(pf - p)
      if (trim(solvers(i)) == 'exact' .and. u > 0.0_R8) ok_case = ok_case .and. abs(pf - pstar) <= 1.0e-5_R8 * p
      write(*,'(1X,A8,A,F6.1,A,ES13.5,A,F8.4,A,ES9.2,A)') solvers(i), ' u = ', u, ' m/s: p_face - p = ', pf - p, &
        ' Pa = ', (pf - p) / (rho * a * speed), ' rho c u, F_r = ', Fr, merge('  ok  ', '  FAIL', ok_case)
      ok = ok .and. ok_case
    end do
  end do
  end do

  if (.not. ok) then
    write(*,'(A)') ' FAIL: a reflective face pressure is not below the pressure of the cell for gas moving away from'// &
      ' the face (or not above it for gas moving towards it), or the exact solver is off the exact solution'
    error stop 1
  end if
  write(*,'(A)') ' PASS: exact and HLLC below the cell for gas moving away from the face, above for gas moving'// &
    ' towards it; the exact solver on the exact solution'
end program reflective_pressure_check
