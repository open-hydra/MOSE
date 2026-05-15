module MOSE_Lib_BC_Fluxes_Manifold
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Advanced_Types_m
  use MOSE_Global_m
  use MOSE_Parameters_m
  use MOSE_Lib_BC_Fluxes, only: Face_Index

  implicit none
  public

contains

  subroutine BC_Manifold ( Blkm, Fm, Blks, Fs, T0, g, pstat )
    use FLINT_Lib_Thermodynamic

    type(MOSE_block_type), intent(in) :: Blkm, Blks
    integer, intent(in)   :: Fm, Fs
    real(R8), intent(out) :: T0, g, pstat
    ! Local
    integer  :: im_min, jm_min, km_min, im_max, jm_max, km_max, dirm
    integer  :: is_min, js_min, ks_min, is_max, js_max, ks_max, dirs
    integer  :: im, jm, km, is, js, ks, Face_i, Face_j, Face_k
    real(R8) :: area, area_M, area_S, Trozza, meanTrozza, meanTrozza_S
    real(R8), dimension(1:np) :: cons, meancons, meancons_S, meanprim_S
    real(R8) :: rho, Rgas, v, T, Tdiff, cp
    integer  :: T_i, Tint(2)

    ! Calculate Area of Local Block
    im_min = 1
    jm_min = 1
    km_min = 1
    im_max = Blkm % dim(1)
    jm_max = Blkm % dim(2)
    km_max = Blkm % dim(3)
    select case(Fm)
      case(1)
        im_max = 1
        dirm = 1
      case(2)
        im_min = Blkm % dim(1)
        dirm = 1
      case(3)
        jm_max = 1
        dirm = 2
      case(4)
        jm_min = Blkm % dim(2)
        dirm = 2
      case(5)
        km_max = 1
        dirm = 3
      case(6)
        km_min = Blkm % dim(3)
        dirm = 3
    end select

    area_M = 0.0d0

    do km = km_min, km_max; do jm = jm_min, jm_max; do im = im_min, im_max

      call Face_Index ( Fm, dirm, im, jm, km, Face_i, Face_j, Face_k )

      area = Blkm % dir(dirm) % f(Face_i,Face_j,Face_k) % A
      area_M = area_M + area

    enddo; enddo; enddo

    ! Calculate Mean Properties of Connected Block
    is_min = 1
    js_min = 1
    ks_min = 1
    is_max = Blks % dim(1)
    js_max = Blks % dim(2)
    ks_max = Blks % dim(3)
    select case(Fs)
      case(1)
        is_max = 1
        dirs = 1
      case(2)
        is_min = Blks % dim(1)
        dirs = 1
      case(3)
        js_max = 1
        dirs = 2
      case(4)
        js_min = Blks % dim(2)
        dirs = 2
      case(5)
        ks_max = 1
        dirs = 3
      case(6)
        ks_min = Blks % dim(3)
        dirs = 3
    end select

    area_S = 0.0d0
    meancons = 0.0d0
    meanTrozza = 0.0d0
    do ks = ks_min, ks_max; do js = js_min, js_max; do is = is_min, is_max

      call Face_Index ( Fs, dirs, is, js, ks, Face_i, Face_j, Face_k )

      area = Blks % dir(dirs) % f (Face_i,Face_j,Face_k) % A
      area_S = area_S + area

      cons = prim2cons( Blks % P(1:np,is,js,ks) ) * area

      meancons = meancons + cons

      ! Temperature initial guess for Newton-Raphson
      call co_rotot_Rtot ( Blks % P(1,is,js,ks), rho, Rgas )
      Trozza = Blks % P(np,is,js,ks) / ( rho * Rgas ) * area
      meanTrozza = meanTrozza + Trozza

    enddo; enddo; enddo

    meanTrozza_S = meanTrozza / area_S
    meancons_S = meancons / area_S
    meanprim_S = cons2prim ( meancons_S, meanTrozza_S )
    
    ! Mass FLux
    g = norm2(meancons_S(nu:nw)) * area_S / area_M

    ! Total Temperature
    v = norm2(meanprim_S (nu:nw))
    call co_rotot_Rtot ( meanprim_S(1), rho, Rgas )
    T = meanprim_S(np) / (rho * Rgas)

    T_i = int ( T )
    Tdiff  = T - T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    cp = f_tabT_expr(1,cp_tab,Tint,Tdiff)
    T0 = T + 0.5d0 * v**2 / cp

    ! ! Check T0
    ! meancons_M = meancons_S * area_S / area_M
    ! meanprim_M = cons2prim ( meancons_M, temperature )
    ! v_M = meanprim_M (nu:nw)
    ! call co_rotot_Rtot ( meanprim_M(1:nsc), rho, Rgas )
    ! T_M = meanprim_M(np) / (rho * Rgas)
    ! T_i = int ( T_M )
    ! Tdiff  = T_M - T_i
    ! Tint(1) = T_i
    ! Tint(2) = T_i + 1
    ! cp_M = f_tabT_expr(1,cp_tab,Tint,Tdiff)
    ! T0_M = T_M + 0.5d0 * dot_product (v_M, v_M) / cp_M
  
    ! Static Pressure
    pstat = meanprim_S (np)

  end subroutine BC_Manifold

end module MOSE_Lib_BC_Fluxes_Manifold