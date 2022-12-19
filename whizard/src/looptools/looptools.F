! WHIZARD <<Version>> <<Date>>

! F77 LoopTools wrapper
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call LoopTools subroutine CputC to calculate the 3-point integral J0
      double complex function J0_LoopTools (p2, k2, q2, m2c)
      implicit none
#include "looptools.h"
      double precision p2, k2, q2
      double complex m2c
      double complex p2c, k2c, q2c
      double complex res(Ncc)
      p2c = dcmplx(p2)
      k2c = dcmplx(k2)
      q2c = dcmplx(q2)
      call ltini
      call CputC (res, p2c, k2c, q2c, m2c, dcmplx(0d0), m2c)
      call clearcache
      J0_LoopTools = 2d0 * res(1)
      end
