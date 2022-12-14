#include "rundeck_opts.h"

!****
!**** Specify horizontal resolution of ocean and
!**** import vertical resolution from olayers module
!****

module oceanres
  use olayers, only : lmo,lmo_min,dZO
  implicit none

  integer, parameter :: & ! number of gridcells in
       imo = 288, & ! longitude
       jmo = 180    ! latitude

end module oceanres
