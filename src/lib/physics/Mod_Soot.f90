module MOSE_Mod_Soot
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  
  integer, private :: i_C2H2, i_CO2, i_CO, i_H2, i_O2, i_CS, i_C6H6 ! indici specie utili

contains

  subroutine Setup_Soot()
    use MOSE_Config_Types_m, only: obj_soot
    use MOSE_Global_m
    use MOSE_Parameters_m
    use FLINT_Lib_Thermodynamic, only: species_names
    implicit none
    integer :: i

    if (trim(obj_soot%model) == 'none') then
      obj_soot%enabled = .false.
      obj_soot%description = 'Soot model disabled'
    else if (index(obj_soot%model, 'LL91') > 0) then
      obj_soot%enabled = .true.
      obj_soot%LL91 = .true.
      obj_soot%LIN = .false.
      obj_soot%description = 'Leung and Lindstedt (1991) model for soot formation and oxidation'
    else if (index(obj_soot%model, 'LIN') > 0) then
      obj_soot%enabled = .true.
      obj_soot%LL91 = .false.
      obj_soot%LIN = .true.
      obj_soot%description = 'Linear model for soot formation and oxidation'
    else
      obj_soot%error_message = 'Unknown soot model: '//trim(obj_soot%model)
    end if
    
    nsoot = 0
    if (.not.obj_soot%enabled) return

    ! inizializzo
    i_O2   = -1
    i_C2H2 = -1
    i_CO2  = -1
    i_CO   = -1
    i_H2   = -1
    i_CS   = -1
    i_C6H6 = -1

    do i = 1, nsc
      select case (trim(adjustl(species_names(i))))
      case ('O2')
        i_O2 = i
      case ('C2H2')
        i_C2H2 = i
      case ('CO2')
        i_CO2 = i
      case ('CO')
        i_CO = i
      case ('H2')
        i_H2 = i
      case ('C(gr)')
        i_CS = i
      case ('C6H6')
        i_C6H6 = i
      end select
    end do
    
    if (i_CS==-1) then
      nsoot = 2
    else
      nsoot = 1
    end if
    
  end subroutine Setup_Soot


  subroutine Soot_Source_Terms ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_soot
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b
    logical :: LL91, LIN

    LL91 = obj_soot%LL91
    LIN  = obj_soot%LIN

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      call Soot_Source_blk ( domain % blk(b) % P,    &
                             domain % blk(b) % r,    &
                             domain % blk(b) % vol,  &
                             domain % blk(b) % dim,  &
                             LL91, LIN )
    
    enddo

  contains

    subroutine Soot_Source_blk ( prim, res, vol, n, LL91, LIN )
      use MOSE_Global_m
      use MOSE_Lib_Soot
      use FLINT_Lib_Thermodynamic
      implicit none
      integer, intent(in) :: n(3)
      logical, intent(in) :: LL91
      logical, intent(in) :: LIN
      real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: prim
      real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: res
      real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
      ! Local
      integer :: i, j, k, s
      real(R8) :: rho, Rgas, T
      real(R8) :: Yc2h2, Yo2, C_c2h2, C_o2, Y_soot, Nvol
      real(R8) :: R1, R2, R3, R4, S_rhoYs, S_rhoN, omega(nsc), energy

      !$omp do collapse (3) private ( rho, Rgas, T, Yc2h2, Yo2, Y_soot, Nvol, C_c2h2, C_o2 ) &
      !$omp private ( R1, R2, R3, R4, S_rhoYs, S_rhoN, omega, energy, i, j, k, s )
      
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)

        call Co_rotot_Rtot ( Prim(1:nsc,i,j,k), rho, Rgas )
        T = prim(np,i,j,k) / (rho * Rgas)

        ! --- dati specie / soot ---
        Yc2h2  = prim(i_C2H2,i,j,k) / rho
        Yo2    = prim(i_O2  ,i,j,k) / rho
        if (nsoot==1) then
          Y_soot = prim(i_CS,i,j,k) / rho        ! soot mass fraction
          Nvol   = prim(nc  ,i,j,k)              ! number per volume [#/m^3]
        else
          Y_soot = prim(nc  ,i,j,k) / rho        ! soot mass fraction
          Nvol   = prim(nc+1,i,j,k)              ! number per volume [#/m^3]
        end if

        ! concentrazioni molari (kmol/m^3)
        C_c2h2 = rho * Yc2h2 / Wm_tab(i_C2H2)
        C_o2   = rho * Yo2   / Wm_tab(i_O2)

        ! --- reaction rate modello LL91 (R1,R2,R3 in kmol/m^3/s di C, R4 in #/m^3/s) ---
        if (LL91) then
          call LL91_soot_rates(rho, T, C_c2h2, C_o2, Y_soot, Nvol, R1, R2, R3, R4)
        elseif (LIN) then
          call LIN_soot_rates(rho, T, C_c2h2, C_o2, Y_soot, Nvol, R1, R2, R3, R4)
        endif

        ! --- sorgenti per i residui ---
        ! massa soot: d(rho*Ys)/dt = Mc*(2R1 + R2 - R3)  [kg/m^3/s]
        S_rhoYs = Mc * ( 2.0*R1 + R2 - R3 )

        ! numero particelle per volume: d(rhoN)/dt = R4  [#/m^3/s]
        S_rhoN  = R4

        ! Accumulo nei residui
        if (nsoot==2) then
          res(nc  ,i,j,k) = res(nc  ,i,j,k) - S_rhoYs * vol(i,j,k)
          res(nc+1,i,j,k) = res(nc+1,i,j,k) - S_rhoN  * vol(i,j,k)
        else
          omega = 0.0_R8
          ! mass sources, written as production (so negative for consumption)
          omega(i_C2H2) = -(R1+R2) * Wm_tab(i_C2H2)
          omega(i_O2  ) = -0.5*R3  * Wm_tab(i_O2  )
          omega(i_H2  ) = +(R1+R2) * Wm_tab(i_H2  )
          omega(i_CO  ) = +R3      * Wm_tab(i_CO  )
          omega(i_CS  ) = -sum(omega) ! enforcing mass conservation
          ! update species residuals
          res(1:nsc,i,j,k) = res(1:nsc,i,j,k) - omega(1:nsc)*vol(i,j,k)
          ! update number density residual
          res(nc   ,i,j,k) = res(nc   ,i,j,k) - S_rhoN*vol(i,j,k) ! number density
          ! Energy source term
          energy = 0.0_R8
          do s = 1, nsc
            energy = energy + omega(s)*f_tabT(T,s,h_tab)
          enddo
          res(np,i,j,k) = res(np,i,j,k) - energy*vol(i,j,k)
        endif

      enddo; enddo; enddo ! (i, j, k) loop

    end subroutine Soot_Source_blk

  end subroutine Soot_Source_Terms


  subroutine Soot_Diffusive_Flux( flux, prim, T_gradient, &
                                  mu, rho, T, area, normal )
    use MOSE_Global_m
    implicit none

    ! Inputs
    real(R8), intent(in)  :: prim(nprim)
    real(R8), intent(in)  :: T_gradient(3)
    real(R8), intent(in)  :: mu, rho, T, area
    real(R8), intent(in)  :: normal(3)
    ! Outputs
    real(R8), intent(inout) :: flux(nprim)
    ! Constants / model parameters
    real(R8), parameter :: CT = 0.55_R8        ! thermophoretic coeff
    ! Local
    real(R8) :: Ys, Nvol, Thermo
    
    ! Retrieve soot scalars
    Ys   = prim(i_CS)/rho ! adim
    Nvol = prim(nc  )     ! [#/m^3]

    ! Molecular diffusive flux of soot is already computed
    ! Thermophoretic mass flux:
    Thermo = -CT * (mu * Ys) * dot_product(T_gradient, normal)/T * area ! [kg/s]
    flux(i_CS) = flux(i_CS) + Thermo

    ! Equivalent flux on Nv
    if (Ys>1.0e-10_R8) flux(nc  ) = flux(i_CS)/(rho*Ys/Nvol)

  end subroutine Soot_Diffusive_Flux

end module MOSE_Mod_Soot
