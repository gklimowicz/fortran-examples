#include "rundeck_opts.h"

      subroutine obio_oasimhr(ihr,i0,j0,idm,jdm,lmo)
c
c  Matches up OASIM spectral irradiance data (just above the surface)
c  to the hour.
c

      USE obio_dim
      USE obio_forc, only :Eda2,Esa2,Ed,Es
      USE obio_com,  only :inwst,inwnd,jnwst,jnwnd

!NEED TO CHECK FOR HYCOM (lMO)

      implicit none

      integer, intent(in) :: idm,jdm,lmo
      integer ihr,ihs,nl,i0,j0
      real    tot


c  Obtain surface irradiance above the surface for the hour
      ihs = int((ihr+1)/2)
      inwst = idm
      jnwst = jdm
      inwnd = 0
      jnwnd = 0


       tot = 0.0
       do nl = 1,nlt
        Ed(nl) = Eda2(nl,ihs)
        Es(nl) = Esa2(nl,ihs)
        tot = tot + Ed(nl)+Es(nl)
       enddo

       if (tot .ge. 0.1)then
        inwst = min(i0,inwst)
        jnwst = min(j0,jnwst)

        inwnd = max(i0,inwnd)
        jnwnd = max(j0,jnwnd)
       endif

      return
      end
