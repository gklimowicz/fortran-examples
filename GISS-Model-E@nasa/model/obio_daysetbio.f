#include "rundeck_opts.h"

      subroutine obio_daysetbio(vrbos,i,j,kdm,nstep,kmax)
c
c  Sets daily parameters for bio.
c

      USE obio_dim
      USE obio_incom, only : rmumax,cchl,cnratio,obio_wsd,obio_wsh
     .                      ,Fescavrate,rik,obio_wss
      USE obio_com,   only : tfac,rmuplsr,rikd,wshc,Fescav,dp1d
     .                      ,avgq1d,gcmax1d,temp1d,obio_P,tzoo,sday


      implicit none

      integer i,j,k

      integer :: nt,ikd
      real :: tfac20,tfact,tfac2,rlightlo,rlighthi,rmidq,fac
     .       ,gcmaxd,wstmp
      real :: cchld       !mean C:chl ratio for the day

      logical vrbos

      integer, intent (in) :: kdm,nstep,kmax

!change: March 15, 2010
      tfac20 = 0.34722*0.851*1.066**20.0
      tfac20 = 0.81*exp(0.0631*20.0)    !Bissinger et al., 2008

c  Compute daily variables
      do k = 1,kdm

c    adjusted to convert to units of /day (instead of doublings/day)
c    (divide by 1.44), and to account for 12-hour photoperiod
c    (divide by 2), making the total factor 0.34722
c    Finally convert to /hr units (instead of /day) by dividing by
c    24.  Normalized to the max growth rate of diatoms
c    Convert to /s  July 016
!change: March 15, 2010
       tfact = 0.34722*0.851*1.066**temp1d(k)     !Eppley, 1972
       tfact = 0.81*exp(0.0631*temp1d(k))    !Bissinger et al., 2008
       tfac(k) = tfact/tfac20
       do nt = 1,nchl
c       rmuplsr(i,k,nt) = rmumax(nt)*rmut*1.12
!       rmuplsr(k,nt) = (rmumax(nt)*tfac(k))/24.0
        rmuplsr(k,nt) = (rmumax(nt)*tfac(k))/sday     !July 2016
       enddo


c  Additional T-dependent factor for cyanobacteria
!#if NCHL_DEFINED > 2
      if (nchl > 2) then
       nt = 3
       tfac2 = temp1d(k)*0.029411 + 0.55823
       tfac2 = min(tfac2,1.0)
       rmuplsr(k,nt) = tfac2*rmuplsr(k,nt)

      endif
!#endif

!@PL option to increase grwoth rate at coasts
#ifdef increase_pp_shelf
       if (dp1d(kmax).le.150.0d0) then
         do nt = 1,nchl
          rmuplsr(k,nt) = rmuplsr(k,nt)*2.0d0
         enddo
       endif
#endif
!@PL

c  Above provides for lower growth of cyanobacteria in cold water.
c  The principle is to keep the delta t difference between cyano to
c  diatoms constant at low T
      enddo


 
c  Zooplankton grazing temperature dependence
c  do not divide by 24: it is a factor
      k = 1
      tzoo = 0.06*exp(0.1*temp1d(k)) + 0.70

cdiag if (vrbos) write(*,'(a,4i5,2d12.4)')'daysetbio: ',
cdiag.   nstep,i,j,k,temp1d(k),tzoo

c
c  Set N:chl ratio and photoadaptation
      rlightlo = 25.0  !threshold for low light adaptation (quanta/m2/s)
      rlighthi = 250.0 !threshold for high light adaptation
      rmidq = (rlighthi-rlightlo)*0.5 + rlightlo
      do k = 1,kdm

       if (avgq1d(k) .gt. 0.0)then

        if (avgq1d(k) .gt. rlighthi)then
         ikd = 3
         !!! cchld = cchl(ikd)    !change: March 15, 2010
         do nt = 1,nchl
          rikd(k,nt) = rik(ikd,nt)
         enddo
        else if (avgq1d(k) .lt. rlightlo)then
         ikd = 1
         !!! cchld = cchl(ikd)    !change: March 15, 2010
         do nt = 1,nchl
          rikd(k,nt) = rik(ikd,nt)
         enddo
        else !avgq<=rlighthi
         ikd = 2
         if (avgq1d(k) .gt. rmidq)then
          fac = (avgq1d(k)-rmidq)/(rlighthi-rmidq)
          !!! cchld = (1.0-fac)*cchl(2) + fac*cchl(3)  !change: March 15, 2010
          do nt = 1,nchl
           rikd(k,nt) = (1.0-fac)*rik(2,nt) + fac*rik(3,nt)
          enddo
         else
          fac = (avgq1d(k)-rlightlo)/(rmidq-rlightlo)
          cchld = (1.0-fac)*cchl(1) + fac*cchl(2)
          do nt = 1,nchl
           rikd(k,nt) = (1.0-fac)*rik(1,nt) + fac*rik(2,nt)
          enddo
         endif
        endif

       else   !avgq<=0
        ikd = 1
        !!! cchld = cchl(ikd)  !change: March 15, 2010
        do nt = 1,nchl
         rikd(k,nt) = rik(ikd,nt)
        enddo
       endif


!change: March15, 2010
!now bn const and in obio_init.f
!      bn(k) = cchld/cnratio  !N:chl ratio (uM/ugl)

      avgq1d(k) = 0.0
      enddo
c


c  Adjustable sinking rate for cocco's: range = 0.3 to 1.4 m/day
!#if NCHL_DEFINED > 3
      if (nchl > 3) then
      nt = 4
      do k = 1,kdm
!      gcmaxd = gcmax1d(k)*24.0
       gcmaxd = gcmax1d(k)*sday     !July 2016

       wstmp = 0.752*gcmaxd + 0.225
       obio_wsd(nt) = max(wstmp,0.3)
       obio_wsd(nt) = min(obio_wsd(nt),1.4)
       obio_wsh(nt) = obio_wsd(nt)/24.0
       obio_wss(nt) = obio_wsd(nt)/sday      !July 2016
!      wshc(k) = obio_wsd(nt)/24.0
       wshc(k) = obio_wsd(nt)/sday           !July 2016
      enddo
      do k = 1,kdm
       gcmax1d(k) = 0.0
      enddo
      endif
!#endif


c  Fe scavenging
      do k = 1,kdm
       if (obio_P(k,4) .lt. 0.6)then
        Fescav(k) = Fescavrate(1)*obio_P(k,4)
       else
        Fescav(k) = Fescavrate(1)* obio_P(k,4)
     .            + Fescavrate(2)*(obio_P(k,4)-0.6)
       endif
      enddo
c
      return
      end
