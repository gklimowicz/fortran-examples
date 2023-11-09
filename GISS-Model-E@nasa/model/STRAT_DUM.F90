#include "rundeck_opts.h"
Subroutine DUMMY_STRAT
!**** Dummy routines in place of STRATDYN
  Entry INIT_GWDRAG
  Entry GWDRAG
  Entry VDIFF
  Entry io_strat
  Entry ALLOC_STRAT_COM
!**** Dummy routines in place of STRAT_DIAG (EP flux calculations)
!**** Note that KEP=0 is set to zero above for the dummy versions.
  Entry EPFLUX
  Entry EPFLXI
  Entry EPFLXP
  Return
End Subroutine DUMMY_STRAT

subroutine get_kep(kep)
  implicit none
  integer :: kep
  kep = 0
end subroutine get_kep
