!> Unit check of the quadratic constitutive relation of Spalart (QCR2000).
!>
!> The solver takes the viscous stress on a face through the pointer Stress_Vector, always as
!> Stress_Vector(velocity gradient, normal, laminar viscosity, eddy viscosity, turbulence variables)
!> (the wall, symmetry and interior fluxes). QCR2000 is not selectable from the input today, so no run
!> passes through it: this program calls it through the same pointer, in the same way, on a 3-D velocity
!> gradient G(i,j) = du_i/dx_j, and compares the stress with the one written from the definition
!> (NASA Turbulence Modeling Resource, Spalart-Allmaras page, SA-QCR2000):
!>   tau_QCR_ij = tau_ij - C_cr1 (O_ik tau_jk + O_jk tau_ik),  O_ik = 2 W_ik / sqrt(G_mn G_mn),  C_cr1 = 0.3,
!> tau_ij = mu_t (2 S_ij - 2/3 G_kk delta_ij) the Boussinesq stress, W_ik = (G_ik - G_ki)/2, plus the laminar
!> stress mu (2 S_ij - 2/3 G_kk delta_ij), which the correction does not touch:
!>   1. laminar viscosity 0, eddy viscosity 1: the stress is tau_QCR . n;
!>   2. laminar viscosity 1, eddy viscosity 0: the stress is the laminar one, without correction.
!> A relation that takes its third argument as the eddy viscosity gives the two cases swapped.
!> The tolerance, 1e-12 of the largest component, covers the order of the floating-point sums.
program qcr2000_check
  use iso_fortran_env, only: R8 => real64
  use MOSE_Lib_Fluid,   only: Stress_Vector
  use MOSE_Lib_QCR2000, only: Stress_Vector_QCR2000
  implicit none
  real(R8), parameter :: c_cr1 = 0.3_R8, tol = 1.0e-12_R8
  real(R8) :: g(3,3), n(3), var(2), tau(3,3), o(3,3), tq(3,3), want(3), got(3), div, den
  integer  :: i, j, k
  logical  :: ok1, ok2

  ! G(i,j) = du_i/dx_j: vorticity and strain along every direction, non-zero divergence
  g = reshape([ 0.30_R8, -1.10_R8,  0.70_R8, &   ! column 1: du_1/dx_1, du_2/dx_1, du_3/dx_1
                0.90_R8,  0.20_R8, -0.40_R8, &   ! column 2
               -0.50_R8,  1.30_R8, -0.10_R8 ], [3, 3])
  n = [0.36_R8, 0.48_R8, 0.80_R8]              ! unit normal
  var = 0.0_R8
  Stress_Vector => Stress_Vector_QCR2000

  div = g(1,1) + g(2,2) + g(3,3)
  do j = 1, 3
    do i = 1, 3
      tau(i,j) = g(i,j) + g(j,i)                 ! 2 S_ij, viscosity 1
    end do
    tau(j,j) = tau(j,j) - 2.0_R8 / 3.0_R8 * div
  end do
  den = sqrt(sum(g**2))
  do k = 1, 3
    do i = 1, 3
      o(i,k) = (g(i,k) - g(k,i)) / den           ! 2 W_ik / sqrt(G_mn G_mn)
    end do
  end do
  do j = 1, 3
    do i = 1, 3
      tq(i,j) = tau(i,j) - c_cr1 * (sum(o(i,:) * tau(j,:)) + sum(o(j,:) * tau(i,:)))
    end do
  end do

  ! 1. eddy viscosity alone
  want = matmul(tq, n)
  got = Stress_Vector(g, n, 0.0_R8, 1.0_R8, var)
  ok1 = maxval(abs(got - want)) <= tol * maxval(abs(want))
  write(*,'(A,3ES24.15)') ' eddy viscosity 1, laminar 0: stress ', got
  write(*,'(A,3ES24.15)') '                     expected (QCR) ', want
  ! 2. laminar viscosity alone
  want = matmul(tau, n)
  got = Stress_Vector(g, n, 1.0_R8, 0.0_R8, var)
  ok2 = maxval(abs(got - want)) <= tol * maxval(abs(want))
  write(*,'(A,3ES24.15)') ' laminar 1, eddy viscosity 0: stress ', got
  write(*,'(A,3ES24.15)') '                 expected (laminar) ', want

  if (ok1 .and. ok2) then
    write(*,'(A)') 'QCR2000Stress: all checks passed'
  else
    write(*,'(A)') 'QCR2000Stress: FAILED (the eddy and the laminar viscosity are not taken in the order of the solver)'
    error stop 1
  end if
end program qcr2000_check
