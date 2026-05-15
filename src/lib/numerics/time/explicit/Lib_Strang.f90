module MOSE_Lib_Strang
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

contains

  subroutine Setup_Strang_Splitting( )
    use MOSE_Config_Types_m, only: obj_chemistry
    implicit none
    
    if (obj_chemistry%use_strang) then
      obj_chemistry%strangcoeff = 0.5_R8
      obj_chemistry%N_strang = 2
    else
      obj_chemistry%strangcoeff = 1._R8
      obj_chemistry%N_strang = 1
    endif

  end subroutine Setup_Strang_Splitting

end module MOSE_Lib_Strang