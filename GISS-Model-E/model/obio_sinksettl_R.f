#include "rundeck_opts.h"
      subroutine obio_sinksettl(vrbos,kmax,errcon,i,j,
     &                    kdm,nstep,dtsrc,ddxypo)

      USE obio_dim
      USE obio_incom,only: wsdeth,mgchltouMC,cchlratio
      USE obio_com, only: P_tend,obio_deltat,D_tend,C_tend
     .                   ,obio_P,det,car
     .                   ,dp1d,wsdet,p1d,obio_ws
     .                   ,rhs,cexp,kzc
#ifdef TRACERS_degC
     .                   ,ndegC1d,Ndeg_tend,cexpdeg !@PL
#endif
      use TimeConstants_mod, only: HOURS_PER_DAY, DAYS_PER_YEAR,
     &                             SECONDS_PER_HOUR


      implicit none

      integer, intent(in) ::kdm,nstep
      real, intent(in) :: dtsrc,ddxypo

      integer :: i,j,k,nt,kmax
      real    :: trnd,sumD,sumD1,sumDdiff
      logical :: vrbos,errcon


!---------------------------------------------------------------
! --- phyto sinking and detrital settling
!---------------------------------------------------------------

!originally: the sinking term is given in units (m/hr)*(mgr,chl/m3)
!in order to be converted into mgr,chl/m3/hr as the tendency
!terms are in the phytoplankton equations, 
!we need to multiply by dz of each layer:
!  dz(k  ) * P_tend(k  ) = dz(k  ) * P_tend(k  ) - trnd
!  dz(k+1) * P_tend(k+1) = dz(k+1) * P_tend(k+1) + trnd
!this way we ensure conservation of tracer after vertical adjustment
!the /hr factor is bcz the obio timestep is in hrs.
! now: all terms in /s, hence tendencies in /s. July 2016

      !phyto sinking
      do nt = nnut+1,ntyp-nzoo
        do k = 1,kmax
        rhs(k,nt,16) = 0.
        enddo
        do k = 1,kmax-1
         trnd = obio_P(k,nt)*obio_ws(k,nt-nnut)
         P_tend(k  ,nt) = P_tend(k,  nt) - trnd/dp1d(k  )
         P_tend(k+1,nt) = P_tend(k+1,nt) + trnd/dp1d(k+1)

         rhs(k  ,nt,16) = rhs(k  ,nt,16) - trnd/dp1d(k  )
         rhs(k+1,nt,16) = rhs(k+1,nt,16) + trnd/dp1d(k+1)
        enddo  ! k
!let phytoplankton that reaches the bottom, disappear in the sediment
!        k = kmax
!        trnd = obio_P(k,nt)*obio_ws(k,nt-nnut)
!        P_tend(k,nt)   = P_tend(k,nt)   - trnd/dp1d(k)
!        rhs(k,nt,16)= - trnd/dp1d(k)
      enddo ! n

     
      !diagnostic for total carbon export at compensation depth
      !total carbon = sinking phyto + settling C detritus
      !term1: sinking phytoplankton

      !detritus settling
      sumD1=sum(D_tend(1:kmax,1)*dp1d(1:kmax))

      do nt = 1,ndet
        do k=1,kmax
        rhs(k,nnut+nchl+nzoo+nt,16) = 0.
        enddo
        do k = 1,kmax-1
         trnd = det(k,nt)*wsdet(k,nt)
         D_tend(k  ,nt) = D_tend(k  ,nt) - trnd/dp1d(k  )
         D_tend(k+1,nt) = D_tend(k+1,nt) + trnd/dp1d(k+1)

         rhs(k  ,nnut+nchl+nzoo+nt,16)= rhs(k  ,nnut+nchl+nzoo+nt,16) 
     .                                - trnd/dp1d(k  )
         rhs(k+1,nnut+nchl+nzoo+nt,16)= rhs(k+1,nnut+nchl+nzoo+nt,16) 
     .                                + trnd/dp1d(k+1)

        enddo  ! k
#ifdef detr_estuarysink
!let detritus that reaches the bottom, disappear in the sediment
       if (dp1d(kmax).le.150.0d0) then
         k = kmax
         trnd = det(k,nt)*wsdet(k,nt)
         D_tend(k,nt)   = D_tend(k,nt)   - trnd/dp1d(k)
         rhs(k,nnut+nchl+nzoo+nt,16)= - trnd/dp1d(k)
       endif
#endif

#ifdef detr_allsink
!let detritus that reaches the bottom, disappear in the sediment
         k = kmax
         trnd = det(k,nt)*wsdet(k,nt)
         D_tend(k,nt)   = D_tend(k,nt)   - trnd/dp1d(k)
         rhs(k,nnut+nchl+nzoo+nt,16)= - trnd/dp1d(k)
#endif
      enddo ! nt

!@PL apply detrital sinking to nondegradable carbon
#ifdef TRACERS_deg
        do k=1,kmax
        rhs(k,ndimndegC,16) = 0.
        enddo
        do k = 1,kmax-1
         trnd = ndegC1d(k)*wsdet(k,1)
         Ndeg_tend(k) = Ndeg_tend(k) - trnd/dp1d(k  )
         Ndeg_tend(k+1) = Ndeg_tend(k+1) + trnd/dp1d(k+1)

         rhs(k  ,ndimndegC,16)= rhs(k  ,ndimndegC,16)
     .                                - trnd/dp1d(k  )
         rhs(k+1,ndimndegC,16)= rhs(k+1,ndimndegC,16)
     .                                + trnd/dp1d(k+1)

        enddo  ! k     

#ifdef detr_estuarysink
!let detritus that reaches the bottom, disappear in the sediment
       if (dp1d(kmax).le.150.0d0) then
         k = kmax
         trnd = ndegC1d(k)*wsdet(k,1)
         Ndeg_tend(k)   = Ndeg_tend(k)   - trnd/dp1d(k)
         rhs(k,ndimndegC,16)= - trnd/dp1d(k)
       endif
#endif

!@PL make entire ocean floor a sink

#ifdef detr_allsink
!let detritus that reaches the bottom, disappear in the sediment
         k = kmax
         trnd = ndegC1d(k)*wsdet(k,1)
         Ndeg_tend(k)   = Ndeg_tend(k)   - trnd/dp1d(k)
         rhs(k,ndimndegC,16)= - trnd/dp1d(k)
#endif
#endif 
 
      sumD=sum(D_tend(1:kmax,1)*dp1d(1:kmax))
      sumDdiff=sumD-sumD1
!     if (vrbos)
!    .write(*,'(a,3i5,3e12.4)')'obio_sinksettl, sumD:'
!    .        ,nstep,i,j,sumD1,sumD,sumDdiff
      
!diagnostic for carbon export at compensation depth
      cexp = 0.

      k=kzc
        do nt=nnut+1,nnut+nchl
           cexp = cexp
     .        + obio_P(k,nt)*obio_ws(k,nt-nnut)   
     .        * cchlratio                       
     .        * SECONDS_PER_HOUR    !July 2016
     .        * HOURS_PER_DAY * DAYS_PER_YEAR         
     .        * 1.d-15 *1.d-3            
     &        * ddxypo              !mgm3 -> PgC/yr              
 
        enddo

      !term2: settling C detritus contribution
      !dont set cexp = 0 here, because adds to before
        nt= 1            !only the for carbon detritus
        cexp = cexp 
     .        + det(k,nt)*wsdet(k,nt)
     .        * SECONDS_PER_HOUR     !July 2016
     .        * HOURS_PER_DAY * DAYS_PER_YEAR
     .        * 1.d-15 *1.d-3                 
     &        * ddxypo               !ugC/l -> PgC/yr
#ifdef TRACERS_degC
!@PL term3 : settling on nondegradable carbon contribution
! only use sinking speed for detrital carbon
        cexp = cexp
     .        + ndegc1d(k)*wsdet(k,nt)
     .        * SECONDS_PER_HOUR     !Jan 2020
     .        * HOURS_PER_DAY * DAYS_PER_YEAR
     .        * 1.d-15 *1.d-3
     &        * ddxypo               !ugC/l -> PgC/yr

        cexpdeg = 0
        cexpdeg = cexpdeg
     .        + ndegc1d(k)*wsdet(k,nt)
     .        * SECONDS_PER_HOUR     !July 2020
     .        * HOURS_PER_DAY * DAYS_PER_YEAR
     .        * 1.d-15 *1.d-3
     &        * ddxypo               !ugC/l -> PgC/yr
#endif
   
      end subroutine obio_sinksettl
