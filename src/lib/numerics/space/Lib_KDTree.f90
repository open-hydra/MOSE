!> Minimal 3D k-d tree for nearest-neighbor queries.
!> Built once from a point cloud, then queried many times (thread-safe for queries).
!> Used by Compute_Yn to replace the O(N_cells x N_walls) brute-force wall-distance
!> search with an O(N_cells x log N_walls) lookup.
!>
!> The query returns the same distance the brute-force scan did, bit for bit:
!> `sqrt` is monotonic and correctly rounded, so taking the minimum of the squared
!> distances and rooting once is exactly the minimum of the rooted distances.
module MOSE_Lib_KDTree
  use iso_fortran_env, only: R8 => real64, I4 => int32

  implicit none
  private

  public :: kdtree_t
  public :: kdtree_build, kdtree_nearest, kdtree_free

  !> Opaque k-d tree. Tree structure is implicit: for the range [lo..hi], the
  !> root of the subtree is at position (lo+hi)/2, elements with smaller coord
  !> along the current axis live in [lo..mid-1], larger in [mid+1..hi]. Axis
  !> cycles x -> y -> z -> x -> ... with depth.
  type :: kdtree_t
    integer               :: n = 0
    real(R8), allocatable :: px(:), py(:), pz(:)
    integer,  allocatable :: perm(:)        ! perm(i) = original index of the point now at position i
  end type kdtree_t

contains

  !> Build the tree from (fx,fy,fz) of length n. O(N log N) average.
  subroutine kdtree_build(tree, fx, fy, fz, n)
    type(kdtree_t), intent(out) :: tree
    integer,        intent(in)  :: n
    real(R8),       intent(in)  :: fx(n), fy(n), fz(n)
    integer :: i, j
    real(R8) :: r

    tree%n = n
    if (n <= 0) return

    allocate(tree%px(n), tree%py(n), tree%pz(n), tree%perm(n))
    do i = 1, n
      tree%px(i)   = fx(i)
      tree%py(i)   = fy(i)
      tree%pz(i)   = fz(i)
      tree%perm(i) = i
    end do

    ! Randomize input order so quickselect with last-element pivot avoids O(N^2)
    ! adversarial inputs (e.g. already-sorted wall centroids). The shuffle only
    ! changes the tree layout, never the distance a query returns, so ranks that
    ! happen to draw different sequences still agree exactly.
    do i = n, 2, -1
      call random_number(r)
      j = 1 + int(r * real(i, R8))
      if (j > i) j = i
      call swap_point(tree%px, tree%py, tree%pz, tree%perm, i, j)
    end do

    call build_rec(tree%px, tree%py, tree%pz, tree%perm, 1, n, 1)
  end subroutine kdtree_build


  !> Free the tree storage.
  subroutine kdtree_free(tree)
    type(kdtree_t), intent(inout) :: tree
    if (allocated(tree%px))   deallocate(tree%px)
    if (allocated(tree%py))   deallocate(tree%py)
    if (allocated(tree%pz))   deallocate(tree%pz)
    if (allocated(tree%perm)) deallocate(tree%perm)
    tree%n = 0
  end subroutine kdtree_free


  !> Query: return the min Euclidean distance from (x,y,z) to any point in the tree.
  !> Read-only on the tree; safe to call concurrently from OpenMP threads.
  subroutine kdtree_nearest(tree, x, y, z, d)
    type(kdtree_t), intent(in)  :: tree
    real(R8),       intent(in)  :: x, y, z
    real(R8),       intent(out) :: d
    real(R8) :: best2

    if (tree%n <= 0) then
      d = huge(1.0_R8)
      return
    end if
    best2 = huge(1.0_R8)
    call nearest_rec(tree%px, tree%py, tree%pz, 1, tree%n, 1, x, y, z, best2)
    d = sqrt(best2)
  end subroutine kdtree_nearest


  ! -------------------------------------------------------------------------
  ! Internal helpers
  ! -------------------------------------------------------------------------

  recursive subroutine build_rec(px, py, pz, perm, lo, hi, axis)
    real(R8), intent(inout) :: px(:), py(:), pz(:)
    integer,  intent(inout) :: perm(:)
    integer,  intent(in)    :: lo, hi, axis
    integer :: mid, next_axis

    if (lo >= hi) return
    mid = (lo + hi) / 2
    call nth_element(px, py, pz, perm, lo, hi, mid, axis)
    next_axis = mod(axis, 3) + 1
    call build_rec(px, py, pz, perm, lo,    mid - 1, next_axis)
    call build_rec(px, py, pz, perm, mid+1, hi,      next_axis)
  end subroutine build_rec


  recursive subroutine nearest_rec(px, py, pz, lo, hi, axis, qx, qy, qz, best2)
    real(R8), intent(in)    :: px(:), py(:), pz(:)
    real(R8), intent(in)    :: qx, qy, qz
    real(R8), intent(inout) :: best2
    integer,  intent(in)    :: lo, hi, axis
    integer :: mid, next_axis
    real(R8) :: d2, dx, dy, dz, diff

    if (lo > hi) return
    mid = (lo + hi) / 2

    dx = px(mid) - qx
    dy = py(mid) - qy
    dz = pz(mid) - qz
    d2 = dx*dx + dy*dy + dz*dz
    if (d2 < best2) best2 = d2

    select case (axis)
      case (1) ; diff = qx - px(mid)
      case (2) ; diff = qy - py(mid)
      case (3) ; diff = qz - pz(mid)
    end select
    next_axis = mod(axis, 3) + 1

    if (diff <= 0.0_R8) then
      ! query on left side: descend left first
      call nearest_rec(px, py, pz, lo,      mid - 1, next_axis, qx, qy, qz, best2)
      if (diff*diff < best2) &
        call nearest_rec(px, py, pz, mid+1, hi,      next_axis, qx, qy, qz, best2)
    else
      call nearest_rec(px, py, pz, mid+1,   hi,      next_axis, qx, qy, qz, best2)
      if (diff*diff < best2) &
        call nearest_rec(px, py, pz, lo,    mid - 1, next_axis, qx, qy, qz, best2)
    end if
  end subroutine nearest_rec


  !> Rearrange (px, py, pz, perm) in [lo..hi] so that the element at position k
  !> is the k-th smallest along the current axis, with smaller-or-equal to the
  !> left and greater-or-equal to the right. Quickselect, average O(hi-lo).
  recursive subroutine nth_element(px, py, pz, perm, lo, hi, k, axis)
    real(R8), intent(inout) :: px(:), py(:), pz(:)
    integer,  intent(inout) :: perm(:)
    integer,  intent(in)    :: lo, hi, k, axis
    integer :: p

    if (lo >= hi) return
    p = partition(px, py, pz, perm, lo, hi, axis)
    if (k < p) then
      call nth_element(px, py, pz, perm, lo,  p - 1, k, axis)
    else if (k > p) then
      call nth_element(px, py, pz, perm, p+1, hi,    k, axis)
    end if
  end subroutine nth_element


  !> Lomuto partition around the pivot at position hi. Returns final pivot index.
  function partition(px, py, pz, perm, lo, hi, axis) result(p)
    real(R8), intent(inout) :: px(:), py(:), pz(:)
    integer,  intent(inout) :: perm(:)
    integer,  intent(in)    :: lo, hi, axis
    integer :: p, i
    real(R8) :: pivot

    select case (axis)
      case (1) ; pivot = px(hi)
      case (2) ; pivot = py(hi)
      case (3) ; pivot = pz(hi)
    end select

    p = lo
    do i = lo, hi - 1
      select case (axis)
        case (1) ; if (px(i) < pivot) then
                     call swap_point(px, py, pz, perm, p, i); p = p + 1
                   end if
        case (2) ; if (py(i) < pivot) then
                     call swap_point(px, py, pz, perm, p, i); p = p + 1
                   end if
        case (3) ; if (pz(i) < pivot) then
                     call swap_point(px, py, pz, perm, p, i); p = p + 1
                   end if
      end select
    end do
    call swap_point(px, py, pz, perm, p, hi)
  end function partition


  subroutine swap_point(px, py, pz, perm, i, j)
    real(R8), intent(inout) :: px(:), py(:), pz(:)
    integer,  intent(inout) :: perm(:)
    integer,  intent(in)    :: i, j
    real(R8) :: t
    integer  :: it
    if (i == j) return
    t = px(i); px(i) = px(j); px(j) = t
    t = py(i); py(i) = py(j); py(j) = t
    t = pz(i); pz(i) = pz(j); pz(j) = t
    it = perm(i); perm(i) = perm(j); perm(j) = it
  end subroutine swap_point

end module MOSE_Lib_KDTree
