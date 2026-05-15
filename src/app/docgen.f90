program MOSE_docgen
    use MOSE_Read_Sim_Param, only: Register_Sim_Param
    use MOSE_Read_IO,        only: Register_IO_Fields, Register_Probes
    use MOSE_Read_Numerics,  only: Register_Numerics
    use MOSE_Read_Physics,   only: Register_Physics
    use MOSE_Backend_INI,    only: Load_Ini, Scan_Ini
    use MOSE_Input_Registry
    use MOSE_Config_Types_m
    implicit none

    allocate(obj_chemistry%RT(1:3+1), obj_chemistry%AT(1:3+1))

    ! Build registry entries
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(1, 'probe-placeholder')  ! Register one probe with placeholder name
    call Register_Numerics(2)
    call Register_Physics()

    call reg%generate_markdown('docs/user/registry.md')

end program MOSE_docgen
    
