!> Unit check of the surface reactions of a carbon wall (BC 505 and 506), GSI_Carbon_Bradley (MOSE_Lib_GSI), on a gas of
!> nine species (the eight of the model and N2) in an order other than that of the model's list:
!>   1. pure O2 at (Tw, pO2) = (3000 K, 1 atm), (2500 K, 0.1 atm), (2000 K, 0.01 atm): omega_O2, the net mass flux mdot
!>      and the reaction heat q = sum omega_s h_s(Tw) against the two-site form of Nagle and Strickland-Constable,
!>        dotm_O2 = x k5 p/(1 + k6 p) + (1 - x) k7 p,  x = 1/(1 + k8/(k7 p)),  omega_O2 = -(W_O2/W_C) dotm_O2 / 2,
!>      (docs/theory/gsi.md, thesis eq. 2.106), the values printed by an independent Python evaluation (planM
!>      bradley_nsc.py: -6.6752347090332759e-02, -8.4499048439110793e-03, -2.2868110731406573e-03 kg/(m2 s)), within 1e-12
!>      (the two evaluate the same closed form in double precision: a few ulp); the species other than O2 and CO
!>      exactly 0. The form of b3e2ba8, 1 + k6 p x in the denominator, gives 4.69, 3.42 and 1.40 times these values;
!>   2. no O2 at all (p_O2 = 0, the gas of a column of N2 and H2O): every source term finite, omega_O2 = 0.
!> h_s(T) = 1000 s T J/kg (s the index of the species in this list), linear: the 1 K interpolation of f_tabT is exact.
program gsi_bradley_check
  use iso_fortran_env, only: R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use MOSE_Global_m, only: nsc
  use MOSE_Lib_GSI, only: GSI_Initialize, GSI_Carbon_Bradley
  use FLINT_Lib_Thermodynamic, only: ns, Tmin, Tmax, wm_tab, Ri_tab, h_tab, Runiv
  implicit none
  integer, parameter :: N = 9, IO2 = 9, ICO = 5, IH2O = 4
  character(len=8), parameter :: names(N) = [character(len=8) :: 'N2', 'CO2', 'O', 'H2O', 'CO', 'OH', 'H2', 'H', 'O2']
  real(R8), parameter :: W(N) = [28.014_R8, 44.009_R8, 15.999_R8, 18.015_R8, 28.010_R8, 17.007_R8, 2.016_R8, 1.008_R8, &
                                 31.998_R8]
  real(R8), parameter :: pref = 1.01325e5_R8, tol = 1.0e-12_R8
  real(R8), parameter :: Tw(3) = [3000.0_R8, 2500.0_R8, 2000.0_R8], pO2(3) = [1.0_R8, 0.1_R8, 0.01_R8]
  real(R8), parameter :: want_O2(3) = [-6.6752347090332759e-02_R8, -8.4499048439110793e-03_R8, -2.2868110731406573e-03_R8]
  real(R8) :: roi(N), omega(N), q, mdot, d, want_mdot, want_q, e_o2, e_m, e_q
  integer  :: k, s, iT
  logical  :: ok, ok_all

  ns = N; nsc = N
  allocate(wm_tab(N), Ri_tab(N))
  wm_tab = W
  Ri_tab = Runiv / wm_tab
  Tmin = 0; Tmax = 5000
  allocate(h_tab(Tmin:Tmax, N))
  do s = 1, N
    do iT = Tmin, Tmax
      h_tab(iT, s) = 1000.0_R8 * s * iT
    end do
  end do
  call GSI_Initialize ( names )

  ok_all = .true.
  do k = 1, 3
    roi = 0.0_R8
    roi(IO2) = pO2(k) * pref / ( Ri_tab(IO2) * Tw(k) )
    call GSI_Carbon_Bradley ( roi, Tw(k), omega, q, mdot )
    d = -want_O2(k) / ( 0.5_R8 * W(IO2) )                     ! carbon removed, over W_C
    want_mdot = d * ( W(ICO) - 0.5_R8 * W(IO2) )
    want_q = want_O2(k) * 1000.0_R8 * IO2 * Tw(k) + d * W(ICO) * 1000.0_R8 * ICO * Tw(k)
    e_o2 = abs(omega(IO2) / want_O2(k) - 1.0_R8)
    e_m  = abs(mdot / want_mdot - 1.0_R8)
    e_q  = abs(q / want_q - 1.0_R8)
    ok = e_o2 <= tol .and. e_m <= tol .and. e_q <= tol .and. all(omega([1,2,3,4,6,7,8]) == 0.0_R8)
    write(*,'(A,F7.1,A,F5.2,A,ES24.16,A,3ES10.2,A,L1)') ' Tw = ', Tw(k), ' K, pO2 = ', pO2(k), ' atm: omega_O2 = ', &
      omega(IO2), ', |rel| O2 mdot q =', e_o2, e_m, e_q, ', others 0: ', all(omega([1,2,3,4,6,7,8]) == 0.0_R8)
    ok_all = ok_all .and. ok
  end do

  ! no O2: N2 and H2O only
  roi = 0.0_R8
  roi(1) = 0.8_R8 * pref / ( Ri_tab(1) * 2000.0_R8 )
  roi(IH2O) = 0.2_R8 * pref / ( Ri_tab(IH2O) * 2000.0_R8 )
  call GSI_Carbon_Bradley ( roi, 2000.0_R8, omega, q, mdot )
  ok = all(ieee_is_finite(omega)) .and. ieee_is_finite(q) .and. ieee_is_finite(mdot) .and. omega(IO2) == 0.0_R8 .and. &
       mdot > 0.0_R8
  write(*,'(A,ES12.4,A,L1)') ' no O2 (N2 0.8 atm, H2O 0.2 atm, 2000 K): mdot = ', mdot, ', finite and omega_O2 = 0: ', ok
  ok_all = ok_all .and. ok

  if (ok_all) then
    write(*,'(A)') 'GSIBradley: all checks passed'
  else
    write(*,'(A)') 'GSIBradley: FAILED'
    error stop 1
  end if
end program gsi_bradley_check
