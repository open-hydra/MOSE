module MOSE_Global_m
  use MOSE_Parameters_m

  implicit none

  character(len=clen) :: MOSE_phase_prefix = ''

  integer :: nsc         ! Number of species
  integer :: nu, nv, nw  ! Velocity/momentum components indexes in prim/residuals vector
  integer :: np          ! Pressure/energy index in prim/residuals vector
  integer :: nt, nrans   ! First RANS variable index in prim/residuals vector and number of RANS variables
  integer :: nc, nsoot   ! Index of first soot variable, Number of soot variables
  integer :: npass=0,nps ! Number of Passive scalar variables
  integer :: nprim       ! Number of primitive variables
  integer :: nres        ! Number of variables whose residuals are printed
  integer :: ndir        ! Number of dimensions of computational frame
  integer :: gc=2        ! ghost cells
  integer :: model       ! Model index (e.g., 0: Euler, 1: laminar NS, 2: RANS, etc.)

end module MOSE_Global_m