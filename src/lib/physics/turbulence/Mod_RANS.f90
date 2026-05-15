module MOSE_Mod_RANS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  
contains

  subroutine Setup_RANS_Model ()
    use MOSE_Config_Types_m, only: obj_rans
    use MOSE_Global_m
    use MOSE_Lib_RANS
    use MOSE_Lib_Fluid
    use MOSE_Lib_Spalart
    use MOSE_Lib_SSGLRR
    use MOSE_Lib_SST
    use MOSE_Lib_Wilcox2006
    use MOSE_Lib_QCR2000

    implicit none

    obj_rans%SpalartShur = .false.
    obj_rans%SAR = .false.
    obj_rans%SAcomp = .false.
    obj_rans%QCR2000 = .false.
    obj_rans%blowing_corr = .false.

    if (trim(obj_rans%model) == 'SAcomp') then
      obj_rans%SAcomp = .true.
    endif

    if (index(trim(obj_rans%model), '-R') > 0) then
      if (index(trim(obj_rans%model), '-RC') > 0) then
        obj_rans%SpalartShur = .true.
      else
        obj_rans%SAR = .true.
      end if
    end if

    if (index(trim(obj_rans%model), '-QCR2000') > 0) then
      obj_rans%QCR2000 = .true.
    end if

    if (index(trim(obj_rans%model), '-blowcorr') > 0) then
      obj_rans%blowing_corr = .true.
    end if

    obj_rans%description = 'RANS model: '//trim(obj_rans%model)

    ! Setting RANS or NS model
    if ( index ( trim(obj_rans%model), 'none' ) > 0 ) then
      nRANS = 0
      Eddy_Viscosity => null()
      RANS_Diffusive_Flux => null()
      Stress_Vector => Stress_Vector_Std
      RANS_Enforce_Realizability => null()

    elseif ( index ( trim(obj_rans%model), 'SA' ) > 0 ) then
      nRANS = 1
      Eddy_Viscosity => Spalart_Eddy_Viscosity
      RANS_Diffusive_Flux => Spalart_RANS_Diffusive_Flux
      if ( obj_rans%QCR2000 ) then
        Stress_Vector => Stress_Vector_QCR2000
      else
        Stress_Vector => Stress_Vector_Std
      end if
      RANS_Source_Terms => Spalart_Source_Terms
      RANS_Set_Wall_Values => Spalart_Set_Wall_Values
      RANS_Set_Blowing_Wall => Spalart_Set_Blowing_Wall
      RANS_Extrapolate_Wall => Spalart_Extrapolate_Wall
      RANS_Enforce_Realizability => Spalart_Enforce_Realizability

    elseif ( index ( trim(obj_rans%model), 'Wilcox2006' ) > 0 ) then
      nRANS = 2
      Eddy_Viscosity => Wilcox2006_Eddy_Viscosity
      RANS_Diffusive_Flux => Wilcox2006_RANS_Diffusive_Flux
      if ( obj_rans%k_energy_coupling ) then
        Stress_Vector => Stress_Vector_2eq
      else
        Stress_Vector => Stress_Vector_Std
      end if
      RANS_Source_Terms => Wilcox2006_Source_Terms
      RANS_Set_Wall_Values => Wilcox2006_Set_Wall_Values
      if (obj_rans%blowing_corr) then
        RANS_Set_Blowing_Wall => SST_Blowing_Correction
      else
        RANS_Set_Blowing_Wall => SST_Blowing_noCorrection
      end if
      RANS_Extrapolate_Wall => Wilcox2006_Extrapolate_Wall
      RANS_Enforce_Realizability => Wilcox2006_Enforce_Realizability
    
    elseif ( index ( trim(obj_rans%model), 'SST' ) > 0 ) then
      nRANS = 2
      Eddy_Viscosity => SST_Eddy_Viscosity
      RANS_Diffusive_Flux => SST_RANS_Diffusive_Flux
      if (obj_rans%k_energy_coupling) then
        Stress_Vector => Stress_Vector_2eq
      else
        Stress_Vector => Stress_Vector_Std
      end if
      RANS_Source_Terms => SST_Source_Terms
      RANS_Set_Wall_Values => SST_Set_Wall_Values
      if (obj_rans%blowing_corr) then
        RANS_Set_Blowing_Wall => SST_Blowing_Correction
      else
        RANS_Set_Blowing_Wall => SST_Blowing_noCorrection
      end if
      RANS_Extrapolate_Wall => SST_Extrapolate_Wall
      RANS_Enforce_Realizability => SST_Enforce_Realizability

    elseif ( index ( trim(obj_rans%model), 'SSGLRR' ) > 0 ) then
      nRANS = 7
      obj_rans%RSM = .true.
      Eddy_Viscosity => SSGLRR_Eddy_Viscosity
      if ( index ( trim(obj_rans%model), '-SD' ) > 0 ) then
        RANS_Diffusive_Flux => SSGLRR_SD_RANS_Diffusive_Flux
      else
        RANS_Diffusive_Flux => SSGLRR_RANS_Diffusive_Flux
      end if
      Stress_Vector => Stress_Vector_RSM
      RANS_Source_Terms => SSGLRR_Source_Terms
      RANS_Set_Wall_Values => SSGLRR_Set_Wall_Values
      RANS_Extrapolate_Wall => SSGLRR_Extrapolate_Wall
      RANS_Enforce_Realizability => SSGLRR_Enforce_Realizability

    else
      nRANS = 0
      Eddy_Viscosity => null()
      RANS_Diffusive_Flux => null()
      Stress_Vector => null()
      obj_rans%model = '  Inviscid flow'

    endif

  end subroutine Setup_RANS_Model

end module MOSE_Mod_RANS