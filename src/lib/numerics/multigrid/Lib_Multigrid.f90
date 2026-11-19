module MOSE_Lib_Multigrid
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64

  implicit none
  private
  public :: Check_Multigrid, Coarse_Grid, Coarse_IOfield
  public :: fine2coarse_prim, coarse2fine_prim
  public :: Restrict_State_Blk, fas_store_q0, fas_build_forcing, Add_FAS_Forcing, fas_coarse_corr
  public :: coarse2fine_corr, Cycle_Turbulence
  public :: Report_FAS_Fallbacks

  !> Cell corrections halved or dropped by this process (coarse2fine_corr)
  integer(I8) :: fas_fallbacks = 0

  !> Iteration of the running level at the current cycle (FAS_Solve), for cycle-turbulence-start (Cycle_Turbulence)
  integer, public :: fas_level_iter = 0

contains

  subroutine Check_Multigrid ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Config_Types_m, only: obj_multigrid
    use ir_precision, only: str
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, d, rap, check

    rap = 2**(obj_multigrid%MGL-1)
    do b = 1, domain % nb
      check = 0
      do d = 1, 3
        if ( Mod ( domain % Blk(b) % Dim(d), rap ) == 0 ) check = check + 1
        if (d > 1) then
          if ( domain % Blk(b) % Dim(d) == 1 ) check = check + 1
        endif
      enddo
      if ( check < 3 ) then
        write(*,'(A90)') ' [ERROR] in Check_Multigrid, block: '//trim(str(.true.,b))
        stop 1
      endif
    enddo

  end subroutine Check_Multigrid


  subroutine fine2coarse_prim ( fPrim, cPrim, fVol, cVol, fDim, cDim )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    implicit none
    integer,  intent(in),  dimension(3) :: fDim, cDim
    real(R8), intent(in),  dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(R8), intent(out), dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    real(R8), intent(in),  dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fVol
    real(R8), intent(in),  dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cVol
    ! Local
    integer  :: i, j, k, i2, j2, k2, i2d, j2d, k2d
    real(R8) :: guess(nprim), rho, Rgas, T0
    real(R8), dimension(nprim,fDim(1),fDim(2),fDim(3)) :: fCons
    real(R8), dimension(nprim,cDim(1),cDim(2),cDim(3)) :: cCons

    !$omp do collapse (2)
    do k = 1, fDim(3)
    do j = 1, fDim(2)
    do i = 1, fDim(1)
      fCons(np+1:nprim,i,j,k) = fPrim(np+1:nprim,i,j,k)
      fCons(1:np,i,j,k) = prim2cons ( fPrim(1:np,i,j,k) )
    end do ; end do ; end do

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)

      i2 = 2*i
      j2 = 2*j
      k2 = 2*k
      i2d = i2-1
      j2d = j2-1
      k2d = k2-1

      if ( fDim(2) == 1 .and. fDim(3) == 1 ) then ! 1D
        cCons(:,i,j,k) = fCons(:,i2d,1,1) * fVol(i2d,1,1)  &
                       + fCons(:,i2, 1,1) * fVol(i2, 1,1)
      else if ( fDim(3) == 1 ) then ! 2D
        cCons(:,i,j,k) = fCons(:,i2d,j2d,1) * fVol(i2d,j2d,1)  &
                       + fCons(:,i2, j2d,1) * fVol(i2, j2d,1)  &
                       + fCons(:,i2d,j2, 1) * fVol(i2d,j2, 1)  &
                       + fCons(:,i2, j2, 1) * fVol(i2, j2, 1)
      else ! 3D
        cCons(:,i,j,k) = fCons(:,i2d,j2d,k2d) * fVol(i2d,j2d,k2d)  &
                       + fCons(:,i2, j2d,k2d) * fVol(i2, j2d,k2d)  &
                       + fCons(:,i2d,j2, k2d) * fVol(i2d,j2, k2d)  &
                       + fCons(:,i2, j2, k2d) * fVol(i2, j2, k2d)  &
                       + fCons(:,i2d,j2d,k2 ) * fVol(i2d,j2d,k2 )  &
                       + fCons(:,i2, j2d,k2 ) * fVol(i2, j2d,k2 )  &
                       + fCons(:,i2d,j2, k2 ) * fVol(i2d,j2, k2 )  &
                       + fCons(:,i2, j2, k2 ) * fVol(i2, j2, k2 )
      end if ! 2D

    enddo; enddo; enddo

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)
      cCons(:,i,j,k) = cCons(:,i,j,k) / cVol(i,j,k)
      ! Guess for Newton-Raphson with one fine cell
      i2 = 2*i
      j2 = 2*j
      k2 = 2*k
      if ( fDim(2) == 1 ) j2 = 1
      if ( fDim(3) == 1 ) k2 = 1
      guess = fPrim(:, i2, j2, k2)
      call co_rotot_Rtot( guess(:), rho, Rgas )
      T0 = guess(np) / ( rho * Rgas )
      cPrim(np+1:nprim,i,j,k) = cCons(np+1:nprim,i,j,k) ! turbulence variables
      cPrim(1:np,i,j,k) = cons2prim ( cCons(1:np,i,j,k), T0 ) ! ( rho U p )
    enddo; enddo; enddo

  end subroutine fine2coarse_prim


  subroutine coarse2fine_prim ( fPrim, cPrim, fDim, cDim )
    use MOSE_Global_m
    implicit none
    integer, intent(in), dimension(3) :: fDim, cDim
    real(R8), intent(out), dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(R8), intent(in),  dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    ! Local
    integer :: i, j, k, rap, i2, j2, k2, i2d, j2d, k2d, im, jm, km, ip, jp, kp
    integer :: ii, jj, kk, counter, mask(3), id(6)
    real(R8) :: a1, a2, a3, a4, coeffs(8), interp(nprim)

    a1 = 27d0/64d0
    a2 = 9d0/64d0
    a3 = 3d0/64d0
    a4 = 1d0/64d0
    coeffs(1:8) = [ a1, a2, a2, a2, a3, a3, a3, a4 ]
    rap = 2
    
    if ( fDim(2) == 1 .and. fDim(3) == 1 ) then ! 1D
      a1 = 3.d0/4.d0
      a2 = 1.d0/4.d0
      coeffs(1:8) = [ a1, a2, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 ]
    else if ( fDim(3) == 1 ) then ! 2D
      a1 = 9.d0/16.d0
      a2 = 3.d0/16.d0
      a3 = 1.d0/16.d0
      a4 = 0d0
      coeffs(1:8) = [ a1, a2, a2, a4, a3, a4, a4, a4 ]
    end if

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)            
    do i = 1, cDim(1)
            
      i2 = rap*i
      j2 = rap*j
      k2 = rap*k
      i2d = i2-(rap-1)
      j2d = j2-(rap-1)
      k2d = k2-(rap-1)
          
      im = Max (1,i-1)
      jm = Max (1,j-1)
      km = Max (1,k-1)

      ip = Min (cDim(1), i+1)
      jp = Min (cDim(2), j+1)
      kp = Min (cDim(3), k+1)

      if ( fDim(2) == 1 ) then
        j2 = 1
        j2d = 1
      end if ! 1D

      if ( fDim(3) == 1 ) then
        k2 = 1
        k2d = 1
      end if ! 2D

      id(1:6) = [ im, jm, km, ip, jp, kp ]
          
      counter = 1
          
      do kk = k2d, k2          
      do jj = j2d, j2            
      do ii = i2d, i2

        if ( fDim(3) == 1 ) then ! 2D
          if (counter==1) mask(1:3) = [1,2,3]
          if (counter==2) mask(1:3) = [4,2,3]
          if (counter==3) mask(1:3) = [1,5,3]
          if (counter==4) mask(1:3) = [4,5,3]
        else
          if (counter==1) mask(1:3) = [1,2,3]
          if (counter==2) mask(1:3) = [4,2,3]
          if (counter==3) mask(1:3) = [1,5,3]
          if (counter==4) mask(1:3) = [4,5,3]
          if (counter==5) mask(1:3) = [1,2,6]
          if (counter==6) mask(1:3) = [4,2,6]
          if (counter==7) mask(1:3) = [1,5,6]
          if (counter==8) mask(1:3) = [4,5,6]
        endif  

        interp = coeffs(1) * cPrim(:,i,j,k) &
               + coeffs(2) * cPrim(:,id(mask(1)),j,k) &
               + coeffs(3) * cPrim(:,i,id(mask(2)),k)    &
               + coeffs(4) * cPrim(:,i,j,id(mask(3)))    &
               + coeffs(5) * cPrim(:,id(mask(1)),id(mask(2)),k)   &
               + coeffs(6) * cPrim(:,id(mask(1)),j,id(mask(3)))   &
               + coeffs(7) * cPrim(:,i,id(mask(2)),id(mask(3)))   &
               + coeffs(8) * cPrim(:,id(mask(1)),id(mask(2)),id(mask(3)))

        fPrim(:,ii,jj,kk) = interp

        counter = counter + 1          
        
      enddo; enddo; enddo
    
    enddo; enddo; enddo
      
  end subroutine coarse2fine_prim


  !> Restriction of fine2coarse_prim on one block, inside a parallel region (multigrid cycle)
  subroutine Restrict_State_Blk ( fPrim, cPrim, fVol, cVol, fDim, cDim )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    implicit none
    integer,  intent(in),    dimension(3) :: fDim, cDim
    real(R8), intent(in),    dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(R8), intent(inout), dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    real(R8), intent(in),    dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fVol
    real(R8), intent(in),    dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cVol
    ! Local
    integer  :: i, j, k, i2, j2, k2, i2d, j2d, k2d
    real(R8) :: cCons(nprim), guess(nprim), rho, Rgas, T0

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)

      i2 = 2*i
      j2 = 2*j
      k2 = 2*k
      i2d = i2-1
      j2d = j2-1
      k2d = k2-1

      if ( fDim(2) == 1 .and. fDim(3) == 1 ) then ! 1D
        cCons = Fine_Cons ( i2d,1,1 ) * fVol(i2d,1,1)  &
              + Fine_Cons ( i2, 1,1 ) * fVol(i2, 1,1)
      else if ( fDim(3) == 1 ) then ! 2D
        cCons = Fine_Cons ( i2d,j2d,1 ) * fVol(i2d,j2d,1)  &
              + Fine_Cons ( i2, j2d,1 ) * fVol(i2, j2d,1)  &
              + Fine_Cons ( i2d,j2, 1 ) * fVol(i2d,j2, 1)  &
              + Fine_Cons ( i2, j2, 1 ) * fVol(i2, j2, 1)
      else ! 3D
        cCons = Fine_Cons ( i2d,j2d,k2d ) * fVol(i2d,j2d,k2d)  &
              + Fine_Cons ( i2, j2d,k2d ) * fVol(i2, j2d,k2d)  &
              + Fine_Cons ( i2d,j2, k2d ) * fVol(i2d,j2, k2d)  &
              + Fine_Cons ( i2, j2, k2d ) * fVol(i2, j2, k2d)  &
              + Fine_Cons ( i2d,j2d,k2  ) * fVol(i2d,j2d,k2 )  &
              + Fine_Cons ( i2, j2d,k2  ) * fVol(i2, j2d,k2 )  &
              + Fine_Cons ( i2d,j2, k2  ) * fVol(i2d,j2, k2 )  &
              + Fine_Cons ( i2, j2, k2  ) * fVol(i2, j2, k2 )
      end if

      cCons = cCons / cVol(i,j,k)
      ! Guess for Newton-Raphson with one fine cell
      if ( fDim(2) == 1 ) j2 = 1
      if ( fDim(3) == 1 ) k2 = 1
      guess = fPrim(:, i2, j2, k2)
      call co_rotot_Rtot( guess(:), rho, Rgas )
      T0 = guess(np) / ( rho * Rgas )
      cPrim(np+1:nprim,i,j,k) = cCons(np+1:nprim) ! turbulence variables
      cPrim(1:np,i,j,k) = cons2prim ( cCons(1:np), T0 ) ! ( rho U p )

    enddo; enddo; enddo

  contains

    !> Conservative variables of a fine cell (turbulence variables as they are)
    function Fine_Cons ( ii, jj, kk ) result ( fCons )
      integer, intent(in) :: ii, jj, kk
      real(R8) :: fCons(nprim)

      fCons(np+1:nprim) = fPrim(np+1:nprim,ii,jj,kk)
      fCons(1:np) = prim2cons ( fPrim(1:np,ii,jj,kk) )

    end function Fine_Cons

  end subroutine Restrict_State_Blk


  !> RANS equations on the coarse levels of the cycle: cycle-turbulence, one- and two-equation models, from
  !> cycle-turbulence-start iterations after cycle-start. Otherwise they keep the restricted state, as soot and passive
  !> scalars
  logical function Cycle_Turbulence ()
    use MOSE_Global_m,       only: model, nrans
    use MOSE_Config_Types_m, only: obj_multigrid
    implicit none

    Cycle_Turbulence = obj_multigrid%cycle_turbulence .and. model == 2 .and. nrans <= 2 .and. &
                       fas_level_iter >= obj_multigrid%cycle_start + obj_multigrid%cycle_turbulence_start

  end function Cycle_Turbulence


  !> Restricted conservative state U_2h^(0) of the coarse cells, for the correction U_2h - U_2h^(0) (Jameson 1983):
  !> flow variables, RANS ones with turb
  subroutine fas_store_q0 ( cPrim, Q0, cDim, turb )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    implicit none
    integer,  intent(in),    dimension(3) :: cDim
    real(R8), intent(in),    dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    real(R8), intent(inout), dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: Q0
    logical,  intent(in) :: turb
    ! Local
    integer :: i, j, k

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)
      Q0(1:np,i,j,k) = prim2cons ( cPrim(1:np,i,j,k) )
      Q0(np+1:nprim,i,j,k) = 0d0
      if ( turb ) Q0(nt:nprim,i,j,k) = cPrim(nt:nprim,i,j,k)
    enddo; enddo; enddo

  end subroutine fas_store_q0


  !> Forcing term of the coarse cells, P_2h = a * sum(R_h) - R_2h(U_2h^(0)), a = damp-restriction (a = 1 in Jameson
  !> 1983). The residuals are cell integrals, so those of the fine cells are summed. Flow equations, RANS ones with turb
  subroutine fas_build_forcing ( fRes, cRes, cF, fDim, cDim, damp, turb )
    use MOSE_Global_m
    implicit none
    integer,  intent(in),    dimension(3) :: fDim, cDim
    real(R8), intent(in),    dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fRes
    real(R8), intent(in),    dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cRes
    real(R8), intent(inout), dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cF
    real(R8), intent(in) :: damp
    logical,  intent(in) :: turb
    ! Local
    integer  :: i, j, k, i2, j2, k2, i2d, j2d, k2d
    real(R8) :: rsum(nprim)

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)

      i2 = 2*i
      j2 = 2*j
      k2 = 2*k
      i2d = i2-1
      j2d = j2-1
      k2d = k2-1

      if ( fDim(2) == 1 .and. fDim(3) == 1 ) then ! 1D
        rsum(:) = fRes(:,i2d,1,1) + fRes(:,i2,1,1)
      else if ( fDim(3) == 1 ) then ! 2D
        rsum(:) = fRes(:,i2d,j2d,1) + fRes(:,i2,j2d,1)  &
                + fRes(:,i2d,j2, 1) + fRes(:,i2,j2, 1)
      else ! 3D
        rsum(:) = fRes(:,i2d,j2d,k2d) + fRes(:,i2,j2d,k2d)  &
                + fRes(:,i2d,j2, k2d) + fRes(:,i2,j2, k2d)  &
                + fRes(:,i2d,j2d,k2 ) + fRes(:,i2,j2d,k2 )  &
                + fRes(:,i2d,j2, k2 ) + fRes(:,i2,j2, k2 )
      end if

      cF(1:np,i,j,k) = damp * rsum(1:np) - cRes(1:np,i,j,k)
      cF(np+1:nprim,i,j,k) = 0d0
      if ( turb ) cF(nt:nprim,i,j,k) = damp * rsum(nt:nprim) - cRes(nt:nprim,i,j,k)

    enddo; enddo; enddo

  end subroutine fas_build_forcing


  !> Coarse level of the cycle: forcing term added to the residual of a stage; the residual of the equations not in the
  !> cycle is zeroed
  subroutine Add_FAS_Forcing ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: np, nt, nprim
    use MOSE_Mod_MPI,  only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k, n(3)
    logical :: turb

    turb = Cycle_Turbulence()

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      n = domain % blk(b) % dim

      !$omp do collapse (2)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)
        domain % blk(b) % R(1:np,i,j,k) = domain % blk(b) % R(1:np,i,j,k) + domain % blk(b) % fas_f(1:np,i,j,k)
        if ( turb ) then
          domain % blk(b) % R(np+1:nt-1,i,j,k) = 0.0_R8
          domain % blk(b) % R(nt:nprim,i,j,k) = domain % blk(b) % R(nt:nprim,i,j,k) &
                                              + domain % blk(b) % fas_f(nt:nprim,i,j,k)
        else
          domain % blk(b) % R(np+1:nprim,i,j,k) = 0.0_R8
        end if
      enddo; enddo; enddo

    enddo

  end subroutine Add_FAS_Forcing


  !> Coarse-grid correction U_2h - U_2h^(0), written over U_2h^(0) in Q0 and then prolongated by coarse2fine_prim:
  !> flow variables, RANS ones with turb
  subroutine fas_coarse_corr ( cPrim, Q0, cDim, turb )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    implicit none
    integer,  intent(in),    dimension(3) :: cDim
    real(R8), intent(in),    dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    real(R8), intent(inout), dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: Q0
    logical,  intent(in) :: turb
    ! Local
    integer :: i, j, k

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)
      Q0(1:np,i,j,k) = prim2cons ( cPrim(1:np,i,j,k) ) - Q0(1:np,i,j,k)
      if ( turb ) then
        Q0(np+1:nt-1,i,j,k) = 0d0
        Q0(nt:nprim,i,j,k) = cPrim(nt:nprim,i,j,k) - Q0(nt:nprim,i,j,k)
      else
        Q0(np+1:nprim,i,j,k) = 0d0
      end if
    enddo; enddo; enddo

  end subroutine fas_coarse_corr


  !> Prolongated correction dQ, times damp-prolongation (1 in Jameson 1983), added to the conservative variables of the
  !> fine cells. A cell whose corrected state is not physical (Admissible_State) takes half the correction, down to
  !> 2**-FAS_HALVINGS, then none. With turb the RANS variables are corrected first, in the same way, keeping them positive
  subroutine coarse2fine_corr ( fPrim, dQ, fDim, damp, b, turb )
    use MOSE_Global_m
    use FLINT_Lib_Thermodynamic
    use MOSE_Lib_Newstate,   only: Check_And_Fix_State
    use MOSE_Lib_FAS_Guards, only: Admissible_State, FAS_HALVINGS
    implicit none
    integer,  intent(in),    dimension(3) :: fDim
    real(R8), intent(inout), dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(R8), intent(in),    dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: dQ
    real(R8), intent(in) :: damp
    integer,  intent(in) :: b
    logical,  intent(in) :: turb
    ! Local
    integer  :: i, j, k, n
    logical  :: physical, halved
    real(R8) :: cons0(np), cons(np), prim(nprim), rho, Rgas, temperature, theta
    real(R8) :: rans(nprim)

    !$omp do collapse (2)
    do k = 1, fDim(3)
    do j = 1, fDim(2)
    do i = 1, fDim(1)
      halved = .false.
      if ( turb ) then
        theta = damp
        do n = 0, FAS_HALVINGS
          rans(nt:nprim) = fPrim(nt:nprim,i,j,k) + theta * dQ(nt:nprim,i,j,k)
          physical = all( rans(nt:nprim) > 0d0 )
          if ( physical ) exit
          theta = 0.5d0 * theta
        enddo
        halved = theta < damp
        if ( physical ) fPrim(nt:nprim,i,j,k) = rans(nt:nprim)
      end if
      cons0 = prim2cons ( fPrim(1:np,i,j,k) )
      call co_rotot_Rtot( fPrim(1:nsc,i,j,k), rho, Rgas )
      temperature = fPrim(np,i,j,k) / ( rho * Rgas )
      theta = damp
      do n = 0, FAS_HALVINGS
        cons = cons0 + theta * dQ(1:np,i,j,k)
        physical = sum( cons(1:nsc) ) > 0d0
        if ( physical ) then
          prim = fPrim(:,i,j,k)
          prim(1:np) = cons2prim ( cons, temperature )
          physical = Admissible_State ( prim )
        end if
        if ( physical ) exit
        theta = 0.5d0 * theta
      enddo
      if ( theta < damp .or. halved ) call Count_FAS_Fallback ()   ! once per cell
      if ( .not. physical ) cycle
      fPrim(1:np,i,j,k) = prim(1:np)
      call Check_And_Fix_State ( fPrim(:,i,j,k), b, i, j, k )
    enddo; enddo; enddo

  end subroutine coarse2fine_corr


  !> Counts one halved or dropped cell correction (thread-safe)
  subroutine Count_FAS_Fallback ()
    implicit none

    !$omp atomic
    fas_fallbacks = fas_fallbacks + 1

  end subroutine Count_FAS_Fallback


  !> End of run: number of halved or dropped corrections over the ranks, if any. Collective
  subroutine Report_FAS_Fallbacks ()
    use MOSE_Mod_MPI, only: mpi_is_root, mpi_allreduce_sum_r8
    implicit none
    real(R8) :: total

    call mpi_allreduce_sum_r8 ( real(fas_fallbacks, R8), total )
    if ( mpi_is_root .and. total > 0d0 ) &
      write(*,'(A,I0)') ' Multigrid corrections halved or not applied (corrected state not physical): ', nint(total, I8)

  end subroutine Report_FAS_Fallbacks


  subroutine Coarse_Grid ( Fine, Coarse )
    use MOSE_Advanced_Types_m
    implicit none
    type(MOSE_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b

    do b = 1, Coarse % nb ! Loop over blocks
      call Coarse_Grid_Blk ( Fine % Blk(b) % node, Coarse % Blk(b) % node, &
                             Fine % Blk(b) % dim,  Coarse % Blk(b) % dim )
    enddo

  contains

    subroutine Coarse_Grid_Blk ( fNode, cNode, fDim, cDim )
      use MOSE_Base_Types_m
      implicit none
      integer, intent(in), dimension(3) :: fDim, cDim
      type(MOSE_vector_3D_type), intent(in),  dimension(0:fDim(1),0:fDim(2),0:fDim(3)) :: fNode
      type(MOSE_vector_3D_type), intent(out), dimension(0:cDim(1),0:cDim(2),0:cDim(3)) :: cNode
      ! Local
      integer :: i, j, k, i2, j2, k2

      !$omp do collapse (2)
      do k = 0, fDim(3), 2-mod(fDim(3),2)
      do j = 0, fDim(2), 2-mod(fDim(2),2)
      do i = 0, fDim(1), 2

        i2 = i / 2
        j2 = j / 2
        k2 = k / 2

        if ( fDim(2) == 1 ) then ! 1D
          j2 = j
        end if

        if ( fDim(3) == 1 ) then ! 2D
          k2 = k
        end if

        cNode(i2,j2,k2) % c = fNode(i,j,k) % c
      
      enddo; enddo; enddo
    
    end subroutine Coarse_Grid_Blk

  end subroutine Coarse_Grid


  subroutine Coarse_IOfield ( fIOfield, cIOfield )
    use Lib_ORION_data
    use MOSE_Global_m
    implicit none
    type(ORION_data), intent(in)    :: fIOfield
    type(ORION_data), intent(inout) :: cIOfield
    ! Local
    integer :: b, nb, Ni, Nj, Nk, Onvar

    Onvar = nprim + 3  ! T, gamma, R

    nb = size( fIOfield % block )
    allocate ( cIOfield % block( nb ) )
    
    do b = 1, nb
      Ni = fIOfield % block(b)%Ni /2
      Nj = fIOfield % block(b)%Nj /2
      Nk = fIOfield % block(b)%Nk /2
      Nj = Max ( 1, Nj ) ! 1D case
      Nk = Max ( 1, Nk ) ! 2D case
      cIOfield % block(b) % Ni = Ni
      cIOfield % block(b) % Nj = Nj
      cIOfield % block(b) % Nk = Nk
      allocate( cIOfield%block(b)%mesh(1:3,0:Ni,0:Nj,0:Nk) )
      allocate( cIOfield%block(b)%vars(1:Onvar,Ni,Nj,Nk) )
    enddo

  end subroutine Coarse_IOfield

end module MOSE_Lib_Multigrid