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
        call Compute_Metric_Tensor ( N1, N2, N3, N4, N5, N6, N7, N8, domain % blk(b) % M(i,j,k), domain % blk(b) % dl(i,j,k), domain % blk(b) % vol(i,j,k), &
                                     loc = [b, i, j, k, 0] )
        if ( inverse_metric ) domain % blk(b) % Minv(i,j,k) % c = Inverse3 ( domain % blk(b) % M(i,j,k) % c )
      enddo; enddo; enddo
      !$omp end parallel
      ! A degenerate or inverted cell (Compute_Metric_Tensor wrote the first one): no further block here; every rank
      ! stops together below
      if ( bad_metric ) exit

      ! Compute Normal & Area
      call Compute_Norm_Area ( domain % blk(b) )
    enddo
    call Stop_On_Bad_Metric ( )

    ! A viscous run on the 2-D axisymmetric wedge or on a 3-D sector stops when the axis is exactly at r = 0
    ! (Check_Wedge_Axis)
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
                                        domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, Bm )
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
    ! A degenerate or inverted ghost cell (Compute_Metric_Tensor wrote the first one): every rank stops together
    call Stop_On_Bad_Metric ( )

    ! Walls and 400 faces: ratio of the MUSCL distances behind the face, for the ghost values (Ghost_Wall) and their
    ! positions (Setup_Stencil_Metrics); every record, on every rank
    do i = 1, domain % nbound
      select case ( domain % bc(i) % type )
        case (301, 302, 400, 503:505)
          Bm = domain % bc(i) % b
          domain % bc(i) % sx = Extrapolation_Ratio ( domain % blk(Bm) % node, domain % bc(i) % i, domain % bc(i) % j, &
                                                      domain % bc(i) % k, domain % bc(i) % f, domain % blk(Bm) % dim )
      end select
    enddo

    ! diffusive-metric = stencil: face metrics of the viscous fluxes from the differences of the cell centres (needs the
    ! cell metrics of the ghost cells, above)
    if ( obj_space_scheme % stencil_metric .and. model > 0 ) call Setup_Stencil_Metrics ( domain )

  end subroutine Setup_Metrics


  !> A viscous run on the 2-D axisymmetric wedge (type 200 on the k faces, one cell in k), or on a 3-D sector closed by
  !> type 200 on its k faces, stops at setup when a boundary face of a block so closed has zero area, that is when the
  !> axis is exactly at r = 0, with every diffusive-metric. The faces of type 200 (BC_Rotational_Periodic_Visc) take the
  !> derivative in r across the axis with the ghost cell one cell away, but behind a face of zero area the ghost cell is a
  !> copy of the cell (Ghost_Symmetry with a zero normal): they see half of d u_r/dr, and the hoop stress of the first row
  !> of cells is wrong (by 1/6 of it in a Homann flow on the wedge; on a sector of four cells of 1 degree, by 1/12 next to
  !> the faces of type 200 with diffusive-metric = stencil and by 1/6 with inverse-mean). The meshes of MOSE put the axis
  !> nodes at a small radius, 1e-6 of the radial size of the first cell.
  !> On the polar line of a 3-D mesh (a block around the axis, not closed by type 200) the faces on the axis have zero area
  !> too. The stencil metric places their ghost cells at the cells and is exact there for a linear field; with
  !> diffusive-metric = inverse-mean or mean the viscous force of the first row of cells grows as the mesh is refined (for a
  !> linear field, whose exact force is zero, 31.6, 80.6 and 164.5 N/m^3 with inverse-mean on 5 x 4 x 8, 10 x 8 x 16 and
  !> 20 x 16 x 32 cells of the polar case of MMSResiduals): a viscous run stops at setup and names stencil.
  subroutine Check_Wedge_Axis ( domain )
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: ndir
    use MOSE_Lib_Metrics, only: delthe
    use MOSE_Lib_BC_Fluxes, only: Face_Index
    use MOSE_Mod_MPI, only: is_local_block, mpi_abort_all
    use MOSE_Config_Types_m, only: obj_space_scheme
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    ! Local
    integer :: r, b, dir, fi, fj, fk
    logical :: polar
    logical, allocatable :: closed(:)
    character(len=12) :: metric

    ! The blocks closed by type 200: none when delthe is exactly zero (two ordered comparisons, as in
    ! BC_Rotational_Periodic_Visc); in 2-D every block is then a wedge, in 3-D delthe is not zero only with sectors turned
    ! about the x axis (Check_Mesh_Type), the blocks with a type-200 record on faces 5 and 6, as there
    allocate ( closed(domain % nb) )
    closed = .false.
    if ( .not. ( delthe >= 0d0 .and. delthe <= 0d0 ) ) then
      closed = ( ndir < 3 )
      if ( ndir == 3 ) then
        do r = 1, domain % nbound
          if ( domain % bc(r) % type == 200 .and. domain % bc(r) % f >= 5 ) closed(domain % bc(r) % b) = .true.
        enddo
      end if
    end if
    ! The other blocks of a 3-D mesh, with diffusive-metric = inverse-mean or mean: a polar line
    polar = ( ndir == 3 .and. .not. obj_space_scheme % stencil_metric )
    if ( .not. ( any ( closed ) .or. polar ) ) return
    metric = 'mean'
    if ( obj_space_scheme % inverse_metric ) metric = 'inverse-mean'
    do r = 1, domain % nbound
      b = domain % bc(r) % b
      if ( .not. is_local_block(b) .or. .not. ( closed(b) .or. polar ) ) cycle
      call Face_Index ( domain % bc(r) % f, dir, domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k, fi, fj, fk )
      if ( domain % blk(b) % dir(dir) % f(fi,fj,fk) % a > 0d0 ) cycle
      if ( .not. closed(b) ) then
        write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] block ', b, ', cell', domain % bc(r) % i, domain % bc(r) % j, &
                                              domain % bc(r) % k, ', face ', domain % bc(r) % f, &
                                              ': a face of zero area (as on the axis of a 3-D mesh), with diffusive-metric = '// &
                                              trim(metric)//':'
        write(*,'(A)') ' [ERROR]   with this metric the viscous force of the first row of cells next to such a face grows as the'
        write(*,'(A)') ' [ERROR]   mesh is refined. Use diffusive-metric = stencil, exact there, with symmetry (300) on the axis.'
        call mpi_abort_all ( 'polar line of a 3-D mesh with diffusive-metric = '//trim(metric)//' (see above)' )
      else if ( ndir == 3 ) then
        write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] block ', b, ', cell', domain % bc(r) % i, domain % bc(r) % j, &
                                              domain % bc(r) % k, ', face ', domain % bc(r) % f, &
                                              ': a face of zero area on the 3-D sector (type 200 on the k faces):'
        write(*,'(A)') ' [ERROR]   the axis is exactly at r = 0, and the hoop stress of the first row of cells would be wrong'
        write(*,'(A)') ' [ERROR]   (by up to 1/6 of it in a Homann flow). Put the axis nodes at a small radius, 1e-6 of the'
        write(*,'(A)') ' [ERROR]   radial size of the first cell, as the meshes of MOSE do.'
        call mpi_abort_all ( 'axis of the 3-D sector at r = 0 (see above)' )
      else
        write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] block ', b, ', cell', domain % bc(r) % i, domain % bc(r) % j, &
                                              domain % bc(r) % k, ', face ', domain % bc(r) % f, &
                                              ': a face of zero area on the 2-D axisymmetric wedge (type 200):'
        write(*,'(A)') ' [ERROR]   the axis is exactly at r = 0, and the hoop stress of the first row of cells would be wrong'
        write(*,'(A)') ' [ERROR]   (by 1/6 of it in a Homann flow). Put the axis nodes at a small radius, 1e-6 of the radial'
        write(*,'(A)') ' [ERROR]   size of the first cell, as the meshes of MOSE do.'
        call mpi_abort_all ( 'axis of the 2-D axisymmetric wedge at r = 0 (see above)' )
      end if
    enddo

  end subroutine Check_Wedge_Axis

end module MOSE_Mod_Metrics