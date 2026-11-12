!> Unit check of the pyrolysing ablator (BC 506), GSI_Char_Pyrolysis (MOSE_Lib_GSI), on the nine species of
!> test/common/CHO-N2-viscous (another order) with the pyrolysis-gas table of unit/gsi_ablator_table.dat (the argument;
!> written by unit/gsi_ablator_table.py: the shift CO + H2O = CO2 + H2 of 1 CO + 1 H2O + 0.5 H2, so that every row has
!> the elements of that mixture), at Tw = 1800 K, p = 1e5 Pa, a gas of N2, H2O, CO2, O2, OH and O:
!>   1. phi = 0, h_v0 = 0: the sources, q and mdot of GSI_Carbon_Bradley, bit for bit (the identity 506 = 505);
!>   2. phi = 0, h_v0 = -8.44e5 J/kg: the same sources and mdot, q = q_ch - mdot_ch h_v0 (within 1e-14);
!>   3. phi = 0.25, h_v0 = -8.44e5 J/kg: mdot = (1 + phi) mdot_ch, omega - omega_ch = phi mdot_ch y_pg(p, Tw) with y_pg of
!>      the table (pyrolysis_gas_composition), within 1e-14 of mdot; the elements carried into the gas, sum_s omega_s
!>      n_es W_e / W_s: C = mdot_ch + mdot_pg Y_C, H = mdot_pg Y_H, O = mdot_pg Y_O, Y_e the makeup of the table, within
!>      1e-13 of mdot (the molar masses of the gas are the sums of the element weights to 1.1e-16, planM bradley_nsc.out:
!>      the surface reactions conserve H and O, and their net mass is the carbon removed); q = sum_s omega_s h_s(Tw) -
!>      mdot h_v0 with the linear h_s(T) = 1000 s T J/kg (within 1e-12);
!>   4. the lookup is at the pressure given: at 5e4 and 2e5 Pa (the two rows of the table) the pyrolysis sources differ.
program gsi_ablator_check
  use iso_fortran_env, only: R8 => real64
  use MOSE_Global_m, only: nsc
  use MOSE_Lib_GSI, only: GSI_Initialize, GSI_Carbon_Bradley, GSI_Char_Pyrolysis
  use FLINT_Lib_Thermodynamic, only: ns, Tmin, Tmax, wm_tab, Ri_tab, h_tab, Runiv
  use FLINT_Lib_PyrolysisGas, only: read_pyrolysis_tables, pyrolysis_gas_composition
  implicit none
  integer, parameter :: N = 9
  character(len=8), parameter :: names(N) = [character(len=8) :: 'N2', 'CO2', 'O', 'H2O', 'CO', 'OH', 'H2', 'H', 'O2']
  real(R8), parameter :: W(N) = [28.014_R8, 44.009_R8, 15.999_R8, 18.015_R8, 28.010_R8, 17.007_R8, 2.016_R8, 1.008_R8, &
                                 31.998_R8]
  ! elements C, H, O of each species, and their weights (Cantera)
  integer, parameter :: NE(3,N) = reshape([0,0,0, 1,0,2, 0,0,1, 0,2,1, 1,0,1, 0,1,1, 0,2,0, 0,1,0, 0,0,2], [3, N])
  real(R8), parameter :: WE(3) = [12.011_R8, 1.008_R8, 15.999_R8]
  real(R8), parameter :: Tw = 1800.0_R8, p = 1.0e5_R8, pref = 1.01325e5_R8, hv0 = -8.44e5_R8, phi = 0.25_R8
  real(R8) :: x(N), roi(N), om_ch(N), q_ch, m_ch, om(N), q, m, y_pg(N), y2(N), E(3), Ymk(3), want(3), q_want, mpg
  real(R8) :: om_a(N), om_b(N), q_a, m_a
  character(len=512) :: table_file
  integer  :: s, k, iT, err
  logical  :: ok, ok_all

  call get_command_argument(1, table_file)
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
  err = read_pyrolysis_tables ( trim(table_file), names )
  if ( err /= 0 ) error stop 'GSIAblator: the table cannot be read'
  ! the makeup of the table: 1 C, 3 H, 2 O
  Ymk = [1.0_R8, 3.0_R8, 2.0_R8] * WE / ( WE(1) + 3.0_R8 * WE(2) + 2.0_R8 * WE(3) )

  ! the gas at the wall: mole fractions, partial densities at p and Tw
  x = [0.70_R8, 0.10_R8, 0.002_R8, 0.15_R8, 0.0_R8, 0.003_R8, 0.0_R8, 0.0_R8, 0.045_R8]
  roi = x * p / ( Runiv * Tw ) * W
  call GSI_Carbon_Bradley ( roi, Tw, om_ch, q_ch, m_ch )
  ok_all = m_ch > 0.0_R8
  write(*,'(A,ES24.16)') ' mdot_ch = ', m_ch

  ! 1.
  call GSI_Char_Pyrolysis ( roi, Tw, p, 1, 0.0_R8, 0.0_R8, om, q, m )
  ok = all(om == om_ch) .and. q == q_ch .and. m == m_ch
  write(*,'(A,L1)') ' 1. phi = 0, h_v0 = 0: the sources, q and mdot of the surface reactions, bit for bit: ', ok
  ok_all = ok_all .and. ok

  ! 2.
  call GSI_Char_Pyrolysis ( roi, Tw, p, 1, 0.0_R8, hv0, om, q, m )
  ok = all(om == om_ch) .and. m == m_ch .and. abs(q - (q_ch - m_ch * hv0)) <= 1e-14_R8 * abs(m_ch * hv0)
  write(*,'(A,ES10.2,A,L1)') ' 2. phi = 0, h_v0: q - (q_ch - mdot_ch h_v0) = ', q - (q_ch - m_ch * hv0), ': ', ok
  ok_all = ok_all .and. ok

  ! 3.
  call GSI_Char_Pyrolysis ( roi, Tw, p, 1, phi, hv0, om, q, m )
  mpg = phi * m_ch
  call pyrolysis_gas_composition ( 1, p, Tw, y_pg )
  ok = abs(m / ((1.0_R8 + phi) * m_ch) - 1.0_R8) <= 1e-14_R8 .and. &
       maxval(abs(om - om_ch - mpg * y_pg)) <= 1e-14_R8 * m
  write(*,'(A,ES10.2,A,ES10.2,A,L1)') ' 3. mdot/((1+phi) mdot_ch) - 1 = ', m / ((1.0_R8 + phi) * m_ch) - 1.0_R8, &
    ', max |omega - omega_ch - mdot_pg y_pg| / mdot = ', maxval(abs(om - om_ch - mpg * y_pg)) / m, ': ', ok
  ok_all = ok_all .and. ok
  do k = 1, 3
    E(k) = sum ( om * NE(k,:) * WE(k) / W )
  end do
  want = mpg * Ymk
  want(1) = want(1) + m_ch
  ok = maxval(abs(E - want)) <= 1e-13_R8 * m
  write(*,'(A,3ES24.16)') '    elements C H O into the gas: ', E
  write(*,'(A,3ES24.16)') '    mdot_ch + mdot_pg Y_C, mdot_pg Y_H, mdot_pg Y_O: ', want
  write(*,'(A,ES10.2,A,L1)') '    largest difference / mdot = ', maxval(abs(E - want)) / m, ': ', ok
  ok_all = ok_all .and. ok
  q_want = 0.0_R8
  do s = 1, N
    q_want = q_want + om(s) * 1000.0_R8 * s * Tw
  end do
  q_want = q_want - m * hv0
  ok = abs(q / q_want - 1.0_R8) <= 1e-12_R8
  write(*,'(A,ES10.2,A,L1)') '    q / (sum omega h - mdot h_v0) - 1 = ', q / q_want - 1.0_R8, ': ', ok
  ok_all = ok_all .and. ok

  ! 4.
  call GSI_Char_Pyrolysis ( roi, Tw, 5.0e4_R8, 1, phi, hv0, om_a, q_a, m_a )
  call GSI_Char_Pyrolysis ( roi, Tw, 2.0e5_R8, 1, phi, hv0, om_b, q_a, m_a )
  call pyrolysis_gas_composition ( 1, 5.0e4_R8, Tw, y2 )
  ok = maxval(abs(om_a - om_b)) > 1e-3_R8 * m .and. maxval(abs(om_a - om_ch - mpg * y2)) <= 1e-14_R8 * m
  write(*,'(A,ES10.2,A,L1)') ' 4. the sources at 5e4 and 2e5 Pa differ by ', maxval(abs(om_a - om_b)) / m, &
    ' of mdot, the table read at the pressure given: ', ok
  ok_all = ok_all .and. ok

  if (ok_all) then
    write(*,'(A)') 'GSIAblator: all checks passed'
  else
    write(*,'(A)') 'GSIAblator: FAILED'
    error stop 1
  end if
end program gsi_ablator_check
