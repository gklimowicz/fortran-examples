#include "rundeck_opts.h"

      subroutine obio_sfcirr(noon,rod,ros,vrbos)
 
c  Computes irradiance just below sea surface (i.e., accounts
c  for surface reflectance), determines array portion that is
c  in light, and computes average cosine for direct irradiance.
 
      USE obio_dim
      USE obio_incom, only :rad
      USE obio_forc,  only :rmud,Ed,Es,solz
      USE obio_com,   only :npst,npnd,hour_of_day,day_of_month
      use ocalbedo_mod, only: nlt

      implicit none

      integer iprt,nl
      real rn,sirr,rsza,sinszaw,szaw,rmudl
      real rod(nlt),ros(nlt)
      logical noon,vrbos

      data rn /1.341/  !refractive index of seawater
      data iprt /13/
 
cdiag if (vrbos) write(*,*)'obio_sfcirr, rod= ',rod

c  Compute irradiance for the day, for this hour
      do nl = 1,nlt
       Ed(nl) = Ed(nl)*(1.0-rod(nl))
       Es(nl) = Es(nl)*(1.0-ros(nl))
cdiag    if (vrbos) write(931,'(2i7,2e12.4)')
cdiag.    nstep,nl,Ed(nl),Es(nl)
      enddo

 
c  Noon print-out
      if (noon)then
        sirr = 0.0
        do nl = npst,npnd
         sirr = sirr + Ed(nl) + Es(nl)
        enddo
        write(*,'(a30,2i6,f10.2)')'Surface irradiance W/m2 = ',
     .     day_of_month,hour_of_day,sirr
      endif  !noon
 
c  Compute average cosine for direct irradiance in the water 
c  column given solar zenith angle (in degrees) at surface.
 
      rsza = solz/rad
      sinszaw = sin(rsza)/rn
      szaw = asin(sinszaw)
      rmudl = 1.0/cos(szaw)   !avg cosine direct (1 over)
      rmud = min(rmudl,1.5)
      rmud = max(rmud,0.0)
 
      return
      end
