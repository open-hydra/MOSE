module MOSE_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

contains

  subroutine Setup_Metrics ( domain )
    use MOSE_Base_Types_m
    use MOSE_Advanced_Types_m
    use MOSE_Global_m
    use MOSE_Lib_Metrics
    use MOSE_Lib_Diffusive, only: Inverse3
    use MOSE_Lib_Stencil_Metrics, only: Setup_Stencil_Metrics
    use MOSE_Config_Types_m, only: obj_space_scheme
    use MOSE_Mod_MPI, only: is_local_block, mpi_is_root
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k, g
    logical :: inverse_metric
    integer :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(MOSE_vector_3D_type) :: N1, N2, N3, N4, N5, N6, N7, N8

    call Check_Mesh_Type ( domain )

    ! Wall roughness is modelled by Spalart-Allmaras only (SA-rough). With any
    ! other model, or laminar, drop it up front so every wall is consistently smooth.
    if ( any( domain % bc(:) % k_rough > 0d0 ) .and. .not. ( model == 2 .and. nRANS == 1 ) ) then
      if ( mpi_is_root .and. domain % mg_level == 1 ) write(*,'(A)') &
        ' [WARNING] Wall roughness (k_rough > 0) is modelled only with the SA turbulence model: all walls are treated as smooth.'
      domain % bc(:) % k_rough = 0d0
    end if

    ! Inverse metric of the interior cells, for the viscous face metric (diffusive-metric =
    ! inverse-mean): the face then needs one 3x3 inversion instead of three.
    inverse_metric = obj_space_scheme % inverse_metric .and. model > 0

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if ( inverse_metric .and. .not. allocated(domain % blk(b) % Minv) ) &
        allocate( domain % blk(b) % Minv ( domain % blk(b) % dim(1), domain % blk(b) % dim(2), domain % blk(b) % dim(3) ) )
      !$omp parallel
      ! Compute metric tensor and cell dimension across i,j,k. => block % M, & block % dl      
      !$omp do collapse(3) private(i, j, k, N1, N2, N3, N4, N5, N6, N7, N8)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        N1 % c = domain % blk(b) % node(i-1,j-1,k-1) % c
        N2 % c = domain % blk(b) % node(i-1,j-1,k  ) % c
        N3 % c = domain % blk(b) % node(i-1,j  ,k-1) % c
        N4 % c = domain % blk(b) % node(i-1,j  ,k  ) % c
        N5 % c = domain % blk(b) % node(i  ,j-1,k-1) % c
        N6 % c = domain % blk(b) % node(i  ,j-1,k  ) % c
        N7 % c = domain % blk(b) % node(i  ,j  ,k-1) % c
        N8 % c = domain % blk(b) % node(i  ,j  ,k  ) % c
        call Compute_Metric_Tensor ( N1, N2, N3, N4, N5, N6, N7, N8, domain % blk(b) % M(i,j,k), domain % blk(b) % dl(i,j,k), domain % blk(b) % vol(i,j,k) )
        if ( inverse_metric ) domain % blk(b) % Minv(i,j,k) % c = Inverse3 ( domain % blk(b) % M(i,j,k) % c )
      enddo; enddo; enddo
      !$omp end parallel

      ! Compute Normal & Area
      call Compute_Norm_Area ( domain % blk(b) )
    enddo

    ! A viscous run on the 2-D axisymmetric wedge stops when the axis is exactly at r = 0 (Check_Wedge_Axis)
    if ( model > 0 ) call Check_Wedge_Axis ( domain )

    if (model>1) then

      ! Compute wall minimum distance => block%yn
      call Compute_Yn ( domain % blk(:), domain % bc, domain % nb, domain % nbound, report = domain % mg_level == 1 )
      
      !$omp parallel
      ! Set wall-distance in connected cells
      !$omp do schedule (dynamic, 64) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks)
      do i = 1, domain % nbound
        select case ( domain % bc(i) % type )
          case(101,201) ! block connection
            Bm = domain % bc(i) % b
            Im = domain % bc(i) % i
            Jm = domain % bc(i) % j
            Km = domain % bc(i) % k
            Fm = domain % bc(i) % f
            Bs = domain % bc(i) % bs
            Is = domain % bc(i) % is
            Js = domain % bc(i) % js
            Ks = domain % bc(i) % ks
            call Yn_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                 Is, Js, Ks, domain % blk(Bs) )
          case(102) ! chimera: no donor copy, extend the receiver
            Bm = domain % bc(i) % b
            if (.not. is_local_block(Bm)) cycle
            Im = domain % bc(i) % i
            Jm = domain % bc(i) % j
            Km = domain % bc(i) % k
            Fm = domain % bc(i) % f
            call Yn_Chimera ( Im, Jm, Km, Fm, domain % blk(Bm) )
        end select
      enddo
      !$omp end parallel
    endif
    
    !$omp parallel private(g)
    ! Create the nodes for gc layers of ghost cell (of the symmetry faces the first layer only: below)
    !$omp do schedule (dynamic, 64) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, domain % nbound
      Bm = domain % bc(i) % b
      ! Ghost metrics land in blk(Bm) and bc(i), read only by the rank owning
      ! Bm, so a non-owning rank needs no blk(Bm)%M/dl/vol at all.
      if (.not. is_local_block(Bm)) cycle
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      select case ( domain % bc(i) % type )
        case(101,201) ! block connection
          Bs = domain % bc(i) % bs
          Is = domain % bc(i) % is
          Js = domain % bc(i) % js
          Ks = domain % bc(i) % ks
          Fs = domain % bc(i) % fs
          d11s = domain % bc(i) % d11
          d12s = domain % bc(i) % d12
          d21s = domain % bc(i) % d21
          d22s = domain % bc(i) % d22
          call BC_Connect_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                    Is, Js, Ks, Fs, domain % blk(Bs), d11s, d12s, d21s, d22s, &
                                    domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, &
                                    periodic = (domain % bc(i) % type == 201) )
        case(300) ! symmetry
          call BC_Symmetry_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                     domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, 1 )
        case default  ! Rientra qui anche chimera
          call BC_Extrapolate_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                        domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg )
      end select
    enddo
    !$omp end do
    ! The further layers of the symmetry faces copy the cells g - 1 behind the face, which in a direction with fewer
    ! cells than g are ghost cells of the opposite face: one layer at a time, after every face has written the previous
    ! one (the barrier at the end of each loop), so that no face reads a ghost cell that another one is writing
    do g = 2, gc
      !$omp do schedule (dynamic, 64) private(i, Bm)
      do i = 1, domain % nbound
        if ( domain % bc(i) % type /= 300 ) cycle
        Bm = domain % bc(i) % b
        if (.not. is_local_block(Bm)) cycle
        call BC_Symmetry_Metrics ( domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k, &
                                   domain % bc(i) % f, domain % blk(Bm), &
                                   domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, g )
      enddo
      !$omp end do
    enddo
    !$omp end parallel

    ! diffusive-metric = stencil: face metrics of the viscous fluxes from the differences of the cell centres (needs the
    ! cell metrics of the ghost cells, above)
    if ( obj_space_scheme % stencil_metric .and. model > 0 ) call Setup_Stencil_Metrics ( domain )

  end subroutine Setup_Metrics


  !> A viscous run on the 2-D axisymmetric wedge (type 200 on the k faces, one cell in k) stops at setup when a boundary
  !> face has zero area, that is when the axis is exactly at r = 0, with every diffusive-metric. The faces of the wedge
  !> (BC_Rotational_Periodic_Visc) take the derivative in r across the axis with the ghost cell one cell away, but behind a
  !> face of zero area the ghost cell is a copy of the cell (Ghost_Symmetry with a zero normal): they see half of
  !> d u_r/dr, and the hoop stress of the first row of cells is wrong (by 1/6 of it in a Homann flow). The meshes of MOSE
  !> put the axis nodes at a small radius, 1e-6 of the radial size of the first cell.
  subroutine Check_Wedge_Axis ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: ndir
    use MOSE_Lib_Metrics, only: delthe
    use MOSE_Lib_BC_Fluxes, only: Face_Index
    use MOSE_Mod_MPI, only: is_local_block, mpi_abort_all
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    ! Local
    integer :: r, b, dir, fi, fj, fk

    ! Not a wedge: the same test as BC_Rotational_Periodic_Visc (delthe exactly zero, as two ordered comparisons)
    if ( ndir == 3 .or. ( delthe >= 0d0 .and. delthe <= 0d0 ) ) return
    do r = 1, domain % nbound
      b = domain % bc(r) % b
      if ( .not. is_local_block(b) ) cycle
      call Face_Index ( domain % bc(r) % f, dir, domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k, fi, fj, fk )
      if ( domain % blk(b) % dir(dir) % f(fi,fj,fk) % a > 0d0 ) cycle
      write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] block ', b, ', cell', domain % bc(r) % i, domain % bc(r) % j, &
                                            domain % bc(r) % k, ', face ', domain % bc(r) % f, &
                                            ': a face of zero area on the 2-D axisymmetric wedge (type 200):'
      write(*,'(A)') ' [ERROR]   the axis is exactly at r = 0, and the hoop stress of the first row of cells would be wrong'
      write(*,'(A)') ' [ERROR]   (by 1/6 of it in a Homann flow). Put the axis nodes at a small radius, 1e-6 of the radial'
      write(*,'(A)') ' [ERROR]   size of the first cell, as the meshes of MOSE do.'
      call mpi_abort_all ( 'axis of the 2-D axisymmetric wedge at r = 0 (see above)' )
    enddo

  end subroutine Check_Wedge_Axis

end module MOSE_Mod_Metrics