module MOSE_Base_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none

  !! ------------------------------------------------------
  !! ------------------------------------------------------
  type :: MOSE_tensor_3D_type
    real(R8)                 :: c(3,3)           !> Metric tensor components.
  end type MOSE_tensor_3D_type


  type :: MOSE_vector_3D_type
    real(R8)                 :: c(3)             !> Average cell length components.
  end type MOSE_vector_3D_type


  type :: MOSE_tensor_3D_R3_type
    real(R8)                 :: c(3,3,3)         !> Metric tensor components.
  end type MOSE_tensor_3D_R3_type


  type :: MOSE_f_metrics_type
    real(R8)                 :: N(3)             !> Unit normal vector.
    real(R8)                 :: A                !> Interface area.
  end type MOSE_f_metrics_type


  type :: MOSE_d_metrics_type
    type(MOSE_f_metrics_type), allocatable  :: f(:,:,:)   
  end type MOSE_d_metrics_type
  !! ------------------------------------------------------
  !! ------------------------------------------------------

end module MOSE_Base_Types_m