module MOSE_Lib_Soot
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Parameters_m

  implicit none
  
  real(R8), parameter :: Mc = 12.011_R8              ! kg/kmol
  real(R8), parameter :: rhos = 1800.0_R8

contains

  subroutine LL91_soot_rates(rho, T, C_C2H2, C_O2, Y_soot, Nvol, &
                             R1, R2, R3, R4)
    implicit none
    real(R8), intent(in)  :: rho, T, C_C2H2, C_O2, Y_soot, Nvol
    real(R8), intent(out) :: R1, R2, R3, R4

    ! Constants
    real(R8), parameter :: pi   = 3.141592653589793_R8
    real(R8), parameter :: NA   = 6.02214076e26_R8
    real(R8), parameter :: kB   = 1.380649e-23_R8
    integer,  parameter :: Cmin = 100
    real(R8), parameter :: Ca   = 9.0_R8

    ! Locals
    real(R8) :: C_mol, d_p, S_area, f_Area
    real(R8) :: k1, k2, k3, k4

    ! Molar conc. of soot carbon
    C_mol = rho*max(Y_soot,0.0_R8)/Mc  ! kmol/m^3

    ! Diameter from volume-based number density
    d_p = (6.0_R8*Mc*C_mol / (pi*rhos*Nvol))**(1.0_R8/3.0_R8)

    ! Surface area per volume
    S_area = pi * d_p**2 * Nvol
    f_Area = sqrt(pi*((6.0_R8*Mc)/(pi*rhos))**(2.0_R8/3.0_R8)) * C_mol**(1.0_R8/3.0_R8) * Nvol**(1.0_R8/6.0_R8)

    ! --- Rate constants (LL91) ---
    k1 = 0.1e5_R8 * exp(-21100.0_R8/T)
    k2 = 0.6e4_R8 * exp(-12100.0_R8/T)
    k3 = 0.1e5_R8 * sqrt(T) * exp(-19680.0_R8/T)
    k4 = 2.0_R8 * Ca * ((6.0_R8*Mc)/(pi*rhos))**(1.0_R8/6.0_R8) * sqrt((6.0_R8*kB*T)/rhos)

    ! --- Reaction rates ---
    R1 = k1 * C_C2H2                   ! nucleation
    R2 = k2 * C_C2H2 * f_Area          ! surface growth
    R3 = k3 * S_area * C_O2            ! oxidation
    R4 = (2.0_R8 * NA / Cmin) * R1 &   ! Nucleation-Coagulation [#/m^3/s]
         - k4 * C_mol**(1.0_R8/6.0_R8) * Nvol**(11.0_R8/6.0_R8) 

  end subroutine LL91_soot_rates

  subroutine LIN_soot_rates(rho, T, C_C2H2, C_O2, Y_soot, Nvol, &
                            R1, R2, R3, R4)
    implicit none
    real(R8), intent(in)  :: rho, T, C_C2H2, C_O2, Y_soot, Nvol
    real(R8), intent(out) :: R1, R2, R3, R4

    ! Constants
    real(R8), parameter :: pi   = 3.141592653589793_R8
    real(R8), parameter :: NA   = 6.02214076e26_R8
    real(R8), parameter :: kB   = 1.380649e-23_R8
    integer,  parameter :: Cmin = 100
    real(R8), parameter :: Ca   = 9.0_R8
    real(R8), parameter :: Mc   = 12.011_R8     ! [kg/kmol]
    real(R8), parameter :: rhos = 1800.0_R8     ! soot density [kg/m^3]

    ! Locals
    real(R8) :: C_mol, d_p, S_area
    real(R8) :: k1, k2, k3, k4

    ! Molar concentration of soot carbon
    C_mol = rho * max(Y_soot, 0.0_R8) / Mc   ! kmol/m^3

    ! Mean particle diameter from number density
    d_p = (6.0_R8 * Mc * C_mol / (pi * rhos * Nvol))**(1.0_R8/3.0_R8)

    ! Surface area per unit volume
    S_area = pi * d_p**2 * Nvol

    ! --- Rate constants ---
    ! Lindstedt (2005) nucleation
    k1 = 0.63e4_R8 * exp(-21100.0_R8 / T)
    ! Lindstedt (1994) surface growth (Bockhorn eq. 27.36)
    k2 = 0.1e-11_R8 * exp(-12100.0_R8 / T)
    ! LL91 oxidation
    k3 = 0.1e5_R8 * sqrt(T) * exp(-19680.0_R8 / T)
    ! LL91 coagulation
    k4 = 2.0_R8 * Ca * ((6.0_R8*Mc)/(pi*rhos))**(1.0_R8/6.0_R8) * sqrt((6.0_R8*kB*T)/rhos)

    ! --- Reaction rates ---
    R1 = k1 * C_C2H2                   ! nucleation [kmol/m^3/s]
    R2 = k2 * C_C2H2  * Nvol           ! surface growth [kmol/m^3/s]
    R3 = k3 * S_area * C_O2            ! Oxidation (LL91)
    R4 = (2.0_R8 * NA / Cmin) * R1 &   ! Nucleation-Coagulation [#/m^3/s]
         - k4 * C_mol**(1.0_R8/6.0_R8) * Nvol**(11.0_R8/6.0_R8) 

  end subroutine LIN_soot_rates

  ! TODO: GIANNI che cazzo è sta merda?
  ! subroutine HACA_soot_rates(rho, T, C_C2H2, C_O2, C_H, C_OH, C_H2, C_H2O, Y_soot, Nvol, &
  !                            R1, R2, R3, R4)
  !   implicit none
  !   real(R8), intent(in)  :: rho, T
  !   real(R8), intent(in)  :: C_C2H2, C_O2, C_H, C_OH, C_H2, C_H2O
  !   real(R8), intent(in)  :: Y_soot, Nvol
  !   real(R8), intent(out) :: R1, R2, R3, R4

  !   ! -------------------------------------------------------------------
  !   ! HACA soot model (Appel-Bockhorn-Frenklach, Combust. Flame 121:122–136, 2000)
  !   ! - Nucleation: Lindstedt (2005)
  !   ! - Growth:     HACA (Frenklach & Wang, 1990; Balthasar & Frenklach, 2005)
  !   ! - Oxidation:  HACA (OH + O2 contributions)
  !   ! - Coagulation: -
  !   ! -------------------------------------------------------------------

  !   ! Constants
  !   real(R8), parameter :: pi   = 3.141592653589793_R8
  !   real(R8), parameter :: NA   = 6.02214076e26_R8
  !   real(R8), parameter :: kB   = 1.380649e-23_R8
  !   integer,  parameter :: Cmin = 100
  !   real(R8), parameter :: Ca   = 9.0_R8
  !   real(R8), parameter :: Mc   = 12.011_R8       ! [kg/kmol]
  !   real(R8), parameter :: rhos = 1800.0_R8       ! [kg/m^3]
  !   real(R8), parameter :: Rg   = 8314.462618_R8  ! [J/kmol/K]

  !   ! Locals
  !   real(R8) :: C_mol
  !   real(R8) :: k1, k3, k4
  !   real(R8) :: RT_kcal, chi_CH, a_param, b_param, alpha, denom, chi_S
  !   real(R8) :: fR1, rR1, fR2, rR2, fR3, fR4, fR5
  !   real(R8) :: rxns_s, rgrowth, RO2, ROH

  !   ! Molar concentration of soot carbon
  !   C_mol = rho * max(Y_soot, 0.0_R8) / Mc   ! kmol/m^3

  !   ! --- Nucleation (Lindstedt 2005)
  !   k1 = 0.63e4_R8 * exp(-21100.0_R8 / T)
  !   R1 = k1 * C_C2H2   ! [kmol/m^3/s]

  !   ! -------------------------------------------------------------------
  !   ! --- HACA growth (Frenklach & Wang 1990, Appel et al. 2000)
  !   ! -------------------------------------------------------------------
  !   RT_kcal = 1.9872036e-3_R8 * T         ! kcal/mol
  !   chi_CH  = 2.3e15_R8                   ! C–H sites/cm^2
  !   a_param = 33.167_R8 - 0.0154_R8 * T
  !   b_param = -2.5786_R8 + 0.00112_R8 * T

  !   ! Raw rates [cm^3/mol/site/s]
  !   fR1 = 4.2e13_R8 * exp(-13.0_R8 / RT_kcal) * C_H    * 1.0e-3_R8
  !   rR1 = 3.9e12_R8 * exp(-11.0_R8 / RT_kcal) * C_H2   * 1.0e-3_R8
  !   fR2 = 1.0e10_R8 * T**0.734_R8 * exp(-1.43_R8 / RT_kcal) * C_OH   * 1.0e-3_R8
  !   rR2 = 3.68e8_R8 * T**1.139_R8 * exp(-17.1_R8 / RT_kcal) * C_H2O  * 1.0e-3_R8
  !   fR3 = 2.0e13_R8 * C_H  * 1.0e-3_R8
  !   fR4 = 8.00e7_R8 * T**1.56_R8 * exp(-3.8_R8 / RT_kcal) * C_C2H2 * 1.0e-3_R8
  !   fR5 = 2.2e12_R8 * exp(-7.5_R8 / RT_kcal) * C_O2 * 1.0e-3_R8

  !   ! Steady-state soot radical site concentration chi_S [sites/cm^2]
  !   denom = rR1 + rR2 + fR3 + fR4 + fR5
  !   chi_S = chi_CH * (fR1 + fR2) / denom

  !   ! Alpha (fraction of active surface sites)
  !   alpha = tanh(a_param / log10((rho*Y_soot)/Nvol) + b_param)
  !   if (alpha < 0.0_R8) alpha = 1.0_R8

  !   ! Growth rate per surface [rxns/s·m²]
  !   rxns_s = fR4 * alpha * chi_S * 1.0e4_R8     ! [1/cm²/s] -> [1/m²/s]
  !   rgrowth = 2.0_R8 * rxns_s * Mc / NA         ! [kg/m²/s]
  !   R2 = 2.0_R8 * rxns_s / Na                   ! [kmol/m³/s]  (proportional to Nvol)

  !   ! -------------------------------------------------------------------
  !   ! --- HACA oxidation (O2 + OH)
  !   ! -------------------------------------------------------------------
  !   RO2 = 2.0_R8 * rxns_s * Mc / NA             ! same structure as growth for O2 sites
  !   ROH = 0.13_R8 * C_OH * sqrt(Rg*T/(2.0_R8*pi*17.0_R8)) * Mc  ! simplified term
  !   R3 = (RO2 + ROH) * Nvol / Mc                ! convert to [kmol/m³/s]

  !   ! --- Nucleation (same as LL91)
  !   R4 = (2.0_R8 * NA / Cmin) * R1

  ! end subroutine HACA_soot_rates

end module MOSE_Lib_Soot
