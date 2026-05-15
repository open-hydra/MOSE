module MOSE_Procedures_m
  use MOSE_Wrap_Setup
  use MOSE_Wrap_Solve
  use MOSE_Wrap_Postprocess

  implicit none

  type :: MOSE_type

  contains
    procedure, nopass  :: setup => MOSE_setup
    procedure, nopass  :: solve => MOSE_solve
    procedure, nopass  :: postprocess => MOSE_postprocess
  end type MOSE_type
  
end module MOSE_Procedures_m