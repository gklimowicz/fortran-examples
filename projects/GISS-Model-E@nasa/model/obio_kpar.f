#include "rundeck_opts.h"

      subroutine obio_kpar(kmax,vrbos,i,j,kdm,nstep,dtsrc,dxypo,
     &                     mo,g0m,s0m,grav,oAPRESS,fsr,lsrpd)
 
      USE obio_dim
      USE obio_com,   only : npst,npnd,p1d,Kd,Kd_em2d,Kpar,Kpar_em2d
     .                      ,delta_temp1d,temp1d

      USE DOMAIN_DECOMP_1D, only : DIST_GRID
      use ocalbedo_mod, only: nlt


      implicit none

      integer, intent(in) :: kdm,nstep
      integer, intent(in) :: lsrpd
      real, intent(in) :: dtsrc,dxypo,grav,oAPRESS,
     &                mo(kdm),
     &                g0m(kdm),
     &                s0m(kdm) 
      real, intent(inout) :: fsr(lsrpd)


      integer i,j,k
      integer nl,ih,icd,ich,ntr,kmax

      real Ebotq,actot,bctot,bbctot,a,bt,bb
      real acdom450,bbc(10),Etopq,zd,zirrq,chl,chlm,fac
      real*8 temgsp


      real Edz(nlt,kdm),Esz(nlt,kdm)
      real Euz(nlt,kdm)
      real Edtop(nlt),Estop(nlt)
      real fchl(nchl)
      real g,s,pres,delta_g

      logical vrbos

      data bbc / 0.002, 0.00071, 0.0032, 0.00071, 0.0029,
     .           0.0,   0.0,     0.0,    0.0,     0.0/


      pres = oAPRESS    !surface atm. pressure
      do k = 1,kmax
          !integrate kd to get kpar
          Kpar(k) = 0.0d0
          Kpar_em2d(k) = 0.0d0
          delta_temp1d(k) = 0.0d0
          do nl = npst,npnd
             Kpar(k) = Kpar(k) + Kd(nl,k)   !in W/m2
             !Kpar_qm2s(k) = Kpar_qm2s(k) + Kd_qm2s(nl,k)   !in quanta/m2/s
             Kpar_em2d(k) = Kpar_em2d(k) + Kd_em2d(nl,k)   !in Einstein/m2/d
          enddo !nl
          pres=pres+MO(k)*GRAV*.5d0
          g=G0M(k)/(MO(k)*DXYPO)
          s=S0M(k)/(MO(k)*DXYPO)
          !temperature change due to Kpar
          delta_g =      Kpar(k)              ! W/m2
     .                  * 1.0d0               !  -> Joules/s/m2
     .                  / mo(k)               !  -> Joules/kg/s
     .                  * dtsrc               !  -> Joules/kg
          delta_temp1d(k) = temp1d(k) - TEMGSP(g+delta_g,s,pres)
          !add missing pressure to get to the bottom of layer k
          pres=pres+MO(k)*GRAV*.5d0

       if (vrbos)write(*,'(a,5i6,6e12.4)')'kpar: ',
     .     nstep,i,j,k,kmax,p1d(k),pres,g,s,Kpar(k),delta_temp1d(k)

      enddo


      !compute par ratios
      do k = 1,lsrpd
          fsr(k) = kpar(k) / kpar(1)
      enddo

      return
      end
