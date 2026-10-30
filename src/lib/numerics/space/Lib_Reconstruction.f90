module MOSE_Lib_Reconstruction
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Reconstruction

contains

  elemental subroutine Reconstruction ( prev, local, next, next2, dl0, dl1, dl2, dll, dlr, beta, priml, primr )
    use MOSE_Lib_Limiters, only : rlimiter
    implicit none
    real(R8), intent(in)    :: prev, local, next, next2
    real(R8), intent(inout) :: priml, primr
    real(R8), intent(in) :: dl0, dl1, dl2, dll, dlr
    real(R8), intent(in) :: beta
    ! Local
    real(R8) :: slope0, slope1, slope2, slopel, sloper

    ! Piecewise Linear Reconstruction
    slope0  = ( local - prev ) / dl0
    slope1  = ( next - local ) / dl1
    slope2  = ( next2 - next ) / dl2
    ! Limited slopes, bounded so that the value of a cell at each of its two faces lies between the value of the cell and that
    ! of the neighbour across the face, on a line of cells of any lengths (the condition of Berger, Aftosmis and Murman, AIAA
    ! 2005-0490, eq. 19): the change over half the cell is at most each jump, |slopel|*dll <= |slope1|*dl1 = |next - local|
    ! and <= |slope0|*dl0 = |local - prev|, |sloper|*dlr <= |slope2|*dl2 and <= |slope1|*dl1. The two faces of a cell take the
    ! same limited slope. On a line of equal cells the bounds are 2 and the slopes are the limiter's own.
    slopel = rlimiter ( slope1, slope0, dl1/dll, dl0/dll )
    sloper = rlimiter ( slope2, slope1, dl2/dlr, dl1/dlr )
    priml = local + beta * slopel * dll
    primr =  next - beta * sloper * dlr

  end subroutine Reconstruction

end module MOSE_Lib_Reconstruction