#include "rundeck_opts.h"

      subroutine obio_daysetrad(vrbos,i,j,kdm)
c
c  Sets daily parameters for ocean irradiance.
c
      USE obio_dim
      USE obio_incom, only : ac,bc,nl450,excdom
      use ocalbedo_mod, only: aw
      USE obio_com, only : npst,npnd,obio_P,avgq1d
     .                    ,acdom

      implicit none

      integer :: i,j,k
      integer :: nl,nt
      integer, intent(in) :: kdm
      real :: actot450,atot450

      logical vrbos

c  Compute acdom
      if (nl450.eq.0) stop 'obio_daysetrad: nl450=0'
      !!m = indext2
      do k = 1,kdm
        actot450 = 0.0
        atot450  = 0.0
        do nt = 1,nchl
         actot450 = actot450  + obio_P(k,nnut+nt)*ac(nt,nl450)
        enddo
        atot450 = aw(nl450) + actot450
        do nl = npst,npnd
         acdom(k,nl) = 0.2*atot450*excdom(nl)
        enddo
      enddo
 
      return
      end
