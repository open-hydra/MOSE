module MOSE_Parameters_m

  implicit none

  integer, parameter  :: hlen=1024
  integer, parameter  :: llen=256
  integer, parameter  :: clen=16
  character(len=clen) :: codename = 'MOSE'

  integer, dimension(6,3), parameter ::  guide  = reshape( [ 1, 0, 0, &
                                                            -1, 0, 0, &
                                                             0, 1, 0, &
                                                             0,-1, 0, &
                                                             0, 0, 1, &
                                                             0, 0,-1  ] , shape(guide), order=[2,1] )

end module MOSE_Parameters_m