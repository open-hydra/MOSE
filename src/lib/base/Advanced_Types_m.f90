module MOSE_Advanced_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Base_Types_m
  use MOSE_Parameters_m
  use MOSE_Series_Data_m
  use Lib_ORION_Data

  implicit none

  !! ------------------------------------------------------
  !! ------------------------------------------------------
  !! FUNDAMENTAL TYPES
  type :: block_type
    integer                                :: dim(3)           ! Number of cells in i-j-k (ghost not included)
    real(R8), allocatable                  :: vol(:,:,:)       ! Cell volume
    type(MOSE_vector_3D_type), allocatable :: node(:,:,:)      ! Mesh grid points (including ghost)
    type(MOSE_tensor_3D_type), allocatable :: M(:,:,:)         ! Metric transformation tensor
    type(MOSE_vector_3D_type), allocatable :: dl(:,:,:)        ! Average cell length (in i/j/k direction). eg: dl%c(1) is sqrt(dx**2+dy**2+dz**2) of the cell in the i direction
    type(MOSE_d_metrics_type)              :: dir(3)           ! Direction object. Contains: i-faces, j-faces, k-faces; eg: dir(1)%face(i,j,k)%n
    real(R8), allocatable                  :: yn(:,:,:)        ! Nearest wall distance
  end type block_type

  type :: bc_type
    integer                    :: i, j, k, b, f                             ! ijk coordinates, block and face in which the boundary element is located
    integer                    :: type                                      ! BC type (1,9)
    integer                    :: bs, is, js, ks, fs, d11, d12, d21, d22    ! BC 1 (connection) specifications
    type(MOSE_tensor_3D_type)  :: Mg(2)                                     ! Ghost cell metric tensor
    type(MOSE_vector_3D_type)  :: dlg(2)                                    ! Ghost cell average cell length
    real(R8)                   :: volg(2)                                   ! Ghost cell volume
    real(R8), allocatable      :: Pg(:,:)                                   ! Ghost cell primitive stencil
    integer                    :: ni(2)                                     ! BC chimera
    integer, allocatable       :: donorID(:,:)                              ! BC chimera
    real(R8), allocatable      :: volume_fraction(:)                        ! BC chimera
    real(R8), allocatable      :: ext_flux(:)                               ! Multi-Solver Coupling
  end type bc_type
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! ------------------------------------------------------
  ! PERFECT GAS EXTENSIONS
  type, extends(block_type) :: MOSE_block_type
    real(R8), dimension(:,:,:,:), allocatable :: P, PO                 ! Primitive variables at time n and n-1: { rho(s) vel p rho*r }, r generic RANS variable
    real(R8), dimension(:,:,:,:), allocatable :: R                     ! Residuals
    real(R8), dimension(:,:,:,:), allocatable :: RS1, RS2              ! Implicit smoothing residuals (temporary storage)
    real(R8), dimension(:,:,:), allocatable   :: dtlocal               ! Local time step
    real(R8), dimension(:,:,:), allocatable   :: beta                  ! Shock detector flag
    type(MOSE_tensor_3D_type), allocatable    :: vel_gradient(:,:,:)   ! Gradient of velocity
    real(R8), allocatable                     :: rc_term1(:,:,:)       ! Spalart-Shur rotation/curvature correction terms
    real(R8), allocatable                     :: rc_term2(:,:,:)       ! Spalart-Shur rotation/curvature correction terms
    logical                                   :: no_chem = .false.     ! Per disattivare la chimica su un unico blocco
  end type MOSE_block_type

  type, extends(bc_type) :: MOSE_bc_type
    real(R8)                            :: qw, Tw, Taw, hg, qrad             ! BC viscous wall specifications   
    real(R8)                            :: T0, p0, alpha, beta, mach, pamb   ! BC 4 (inflow/outflow) specifications
    real(R8)                            :: mdot                              ! BC 4 (inflow/outflow) specifications
    real(R8)                            :: rel_fac                           ! BC 4 (inflow/outflow) specifications
    real(R8), allocatable               :: ci(:)                             ! BC 4 (inflow/outflow) specifications
    real(R8)                            :: SF                                ! BC ? specifications
    real(R8)                            :: aCoeff, n, pRef, Taf, haf         ! BC 14 specification (SRM grain)
    ! Extended fields from ATLAS output format
    real(R8)                            :: k_rough  = 0.0_R8                 ! Surface roughness height [m]  (wall BCs 301-304)
    real(R8)                            :: eps_wall = 0.0_R8                 ! Wall emissivity               (wall BCs 301, 302)
    real(R8)                            :: rhoGrain = 0.0_R8                 ! SRM grain density [kg/m3]     (SRM BC 501)
    real(R8)                            :: SF_geo   = 1.0_R8                 ! SRM geometric scale factor    (SRM BC 501)
    real(R8)                            :: psub     = 0.0_R8                 ! Nozzle subsonic  pressure     (nozzle BC 420)
    real(R8)                            :: psup     = 0.0_R8                 ! Nozzle supersonic pressure    (nozzle BC 420)
    real(R8)                            :: rt_nozzle= 0.0_R8                 ! Nozzle mass flux parameter    (nozzle BC 420)
    type(time_series_type)              :: p0time
    type(time_series_type), allocatable :: q2d_map(:)                        ! BC 667 (Q2D mapped) time-varying primitive data (one per primitive)
    character(len=llen)                 :: q2d_file = ' '                    ! BC 667 (Q2D mapped) Tecplot filename from bc.txt
    logical                             :: q2d_periodic = .false.            ! BC 667 (Q2D mapped) periodic time signal flag
  end type MOSE_bc_type
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! ------------------------------------------------------
  ! COMPUOND TYPES
  type :: MOSE_domain_type
    real(R8)                                         :: time             ! Solution time
    real(R8)                                         :: dtglobal         ! Global dt for time accurate simulation
    integer                                          :: iter, itermax    ! Iteration number
    integer                                          :: nb, nbound       ! Number of blocks, number of boundary faces
    integer, dimension(:,:), allocatable             :: n_bf             ! Number of bc elements per faces per block
    type(MOSE_block_type), dimension(:), allocatable :: blk              ! Allocatable block type
    type(MOSE_bc_type),    dimension(:), allocatable :: bc               ! Allocatable bc object
    ! MPI local BC indices (built by build_local_bc_index)
    integer                                          :: n_local_bc = 0
    integer, dimension(:), allocatable               :: local_bc_idx
    integer                                          :: n_local_bs = 0
    integer, dimension(:), allocatable               :: local_bs_idx
  end type MOSE_domain_type

  type :: MOSE_simulation_type
    type(MOSE_domain_type), dimension(:), allocatable  :: domain
    type(ORION_data), dimension(:), allocatable        :: IOfield
  end type MOSE_simulation_type
  !! ------------------------------------------------------
  !! ------------------------------------------------------

end module MOSE_Advanced_Types_m