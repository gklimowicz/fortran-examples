#include "rundeck_opts.h"
      subroutine obio_edeu(kmax,vrbos,i,j,im,jm,kdm,nstep,obio_lambdas)
 
c  Model of irradiance in the water column.  Accounts for three 
c  irradiance streams:
c
c  Edz = direct downwelling irradiance
c  Esz = diffuse downwelling irradiance
c  Euz = diffuse upwelling irradiance
c
c  Propagation is done in energy units, tests are done in quanta,
c  final is quanta for phytoplankton growth.
 

      USE obio_dim
      USE obio_incom, only : ac,bc,facirr,Dmax,bbw
     .                      ,rmus,nl450,zc
      use ocalbedo_mod, only: aw, bw, lam, nlt
      USE obio_forc,  only : Ed,Es,rmud,tirrq,tirrq_critical
      USE obio_com,   only : acdom,npst,npnd,WtoQ,dp1d,avgq1d
     .                      ,obio_P,p1d,Kd,Kd_em2d 
      USE FILEMANAGER, only: openunit,closeunit,file_exists



      implicit none

!      type(DIST_GRID),intent(in) :: ogrid
      integer, intent(in) :: im,jm,kdm,nstep

      real, intent(in) :: obio_lambdas(nlt)

      integer i,j,k
      integer nl,ih,icd,ich,ntr,kmax
      integer :: ichan,iu_bio

      real Ebotq,actot,bctot,bbctot,a,bt,bb
      real acdom450,bbc(10),Etopq,zd,zirrq,chl,chlm,fac
      real*8 temgsp


      real Edz(nlt,kdm),Esz(nlt,kdm) 
      real, ALLOCATABLE, DIMENSION(:)::  delta_temp1d  !change in T due to kpar
      real Euz(nlt,kdm)
      real Edtop(nlt),Estop(nlt)
      real fchl(nchl)
      real g,s,pres,delta_g

      logical vrbos

      data bbc / 0.002, 0.00071, 0.0032, 0.00071, 0.0029,
     .           0.0,   0.0,     0.0,    0.0,     0.0/


      !oAPRESS  => ocnice%APRESS

      !write(913,*) "edeu",i,j,Ed,Es
cddd     &     sum(aw),sum(bw),sum(ac),sum(bc),sum(facirr),sum(lam),
cddd     &     (Dmax),(bbw),
cddd     &     (rmus),(nl450),
cddd     &     sum(Ed),sum(Es),(rmud),sum(tirrq),
cddd     &     sum(acdom),(npst),(npnd),sum(WtoQ),
cddd     &     sum(dp1d),sum(avgq1d),
cddd     &     sum(obio_P),sum(p1d)
 
      if (nl450.eq.0) stop 'obio_edeu: nl450=0'
      do k=1,kdm
      tirrq(k) = 0.0
      enddo
 
       Ebotq = 0.0
c  Convert to quanta: divide by Avos # to get moles quanta; then mult by
c  1E6 to get uM or uEin
       do nl = npst,npnd
        Edtop(nl) = Ed(nl)
        Estop(nl) = Es(nl)
        Ebotq = Ebotq + (Edtop(nl)+Estop(nl))*WtoQ(nl)*1.0E6

!        if (vrbos)
!    .   write(*,'(a,4i7,6e12.4)')'obio_edeu1: ',
!    .        nstep,i,j,nl,
!    .        WtoQ(nl),Ed(nl),Es(nl),Edtop(nl),Estop(nl),Ebotq

       enddo

!       pres = oAPRESS(i,j)    !surface atm. pressure
       do k = 1,kmax
          Etopq = Ebotq
          zd = min(Dmax,p1d(k+1)) 
          zirrq = 0.0
          do nl = npst,npnd
             Edz(nl,k) = 0.0
             Esz(nl,k) = 0.0
             Euz(nl,k) = 0.0
             actot  = 0.0
             bctot  = 0.0
             bbctot = 0.0
             do ntr = 1,nchl
               actot  = actot  + obio_P(k,nnut+ntr)*ac(ntr,nl)
               bctot  = bctot  + obio_P(k,nnut+ntr)*bc(ntr,nl)
               bbctot = bbctot + obio_P(k,nnut+ntr)*bbc(ntr)
     .                          *bc(ntr,nl)
             enddo

             a  = aw(nl) + acdom(k,nl) + actot
             bt = bw(nl) + bctot
             bb = bbw*bw(nl) + bbctot
             bb = max(bb,0.0002)
             if (Edtop(nl) .ge. 1.0E-4 .or. Estop(nl) .ge. 1.0E-4)then
              call radmod(zd,Edtop(nl),Estop(nl),rmud,a,bt,bb,
     .                    Dmax,Edz(nl,k),Esz(nl,k),Euz(nl,k))
             endif
             Edtop(nl) = Edz(nl,k)
             Estop(nl) = Esz(nl,k)
             zirrq = zirrq + (Edz(nl,k)+Esz(nl,k)+Euz(nl,k))
     .                     * WtoQ(nl)*1.0E6

             if (p1d(k+1).le.zd) then
             !Kd_qm2s(nl,k) = (a + bb) / rmus    !in quanta/m2/s
             Kd(nl,k) = Edz(nl,k)+Esz(nl,k)+Euz(nl,k)    !in W/m2

 
             Kd_em2d(nl,k) = Kd(nl,k)*obio_lambdas(nl)*(0.836E-2)*(3600)
     &                       *(24)*(1E-6)  !in Einstein/(m2 day)
             endif

          enddo   !nl

          Ebotq = zirrq
          ih = nint(p1d(k+1))
          acdom450 = acdom(k,nl450)
          ih = min(ih,nh)    !nh=200 set in obio_dim
          ih = max(ih,1)
          icd = nint((alog10(acdom450)*100.0+300.0)/10.0) + 1
          icd = max(icd,1)
          icd = min(icd,ncd)
          chl = 0.0
          do ntr = nnut+1,ntyp-nzoo
             chl = chl + obio_P(k,ntr) 
          enddo
          chlm = max(chl,0.001)
          ich = nint((alog10(chlm)*100.0+300.0)/10.0) + 1
          ich = max(ich,1)
          ich = min(ich,nch)
          chlm = max(chl,0.00001)
          do ntr = 1,nchl
             fchl(ntr) = obio_P(k,nnut+ntr)/chlm
          enddo

          fac = 0.0
          do ntr = 1,nchl
             fac = fac + fchl(ntr)*facirr(ih,ich,ntr,icd)

cdiag        if(vrbos)
cdiag.         write(*,'(a,3i7,e12.4,4i7,2e12.4)')
cdiag.        'edeu diag: ',
cdiag.        nstep,k,ntr,fchl(ntr),ih,ich,ntr,icd,
cdiag.        facirr(ih,ich,ntr,icd),fac
          enddo
 
          tirrq(k) = fac*((Etopq+Ebotq)*0.5)*rmus   !units in quanta/m2/s

cdiag        if (vrbos) then
cdiag        write(*,'(a,4i7,6e12.4)')'obio_edeu2: ',
cdiag.              nstep,i,j,k,fac,Etopq,
cdiag.              Ebotq,rmus,dp1d(k),tirrq(k)
cdiag        endif

       enddo  !k
 
!alternatively, zc is the depth of that light is 1% of top
!     zc=0d0
!     do k=1,kmax
!      !!!definition of zc>=Qtop/100 gives constant zc everywhere
!      !!if (tirrq(1).gt.0. .and. tirrq(k).ge.tirrq(1)/100.) then
!        if (tirrq(1).gt.tirrq_critical) then
!           if (tirrq(k).ge.tirrq_critical) then
!            zc = p1d(k)    !this is really the "shallowest" limit for zc
!           endif
!        else
!            zc = 0d0
!        endif
cdiag write(*,'(a,4i7,3e12.4)')'obio_edeu2: ',
cdiag.    nstep,i,j,k,p1d(k),tirrq(k)
!     enddo
cdiag write(*,'(a,4i7,e12.4)')'obio_edeu3: ',
cdiag.    nstep,i,j,k,zc

c  Irradiance summary loops

      do k = 1,kmax
       avgq1d(k) = avgq1d(k) + tirrq(k)

cdiag write(*,'(a,4i5,3e12.4)')'obio_edeu tirrq: ',
cdiag.       nstep,i,j,k,p1d(k),tirrq(k),zc
      enddo


      return
      end

c-----------------------------------------------------------------------
      subroutine radmod(zd,Edtop,Estop,rmud,a,bt,bb,Dmax,
     .                  Edz,Esz,Euz)
 
c  Model of irradiance in the water column.  Accounts for three 
c  irradiance streams:
c
c  Edz = direct downwelling irradiance
c  Esz = diffuse downwelling irradiance
c  Euz = diffuse upwelling irradiance
c
c  Uses Ackelson's (1994, JGR) mod's to the Aas (1987, AO) 
c  two-stream model.
c
c  Propagation is done in energy units, tests are done in quanta,
c  final is quanta for phytoplankton growth.
c
c  Commented out terms produce a max error of 
c  0.8% in Esz for a > 0.004 and bb > 0.0001 and
c  3.9% in Euz for a > 0.004 and bb > 0.00063
 
      USE obio_incom, only:rmus,rmuu,rd,ru
      USE obio_com, only : EUZ_DEFINED

      implicit none


      real cd,a,bt,rmud,Edz,Edtop
     .    ,zd,au,Bu,bb,Cu,as,Bs,Cs,Bd,Fd,bquad,cquad
     .    ,sqarg,a1,a2,S,SEdz,a2ma1,rM,rN,c2,Ta2z,Esz
     .    ,Eutmp,Euz,Estop,Dmax

 
c  Downwelling irradiance: Edz, Esz
c  Compute irradiance components at depth
      cd = (a+bt)*rmud
      Edz = Edtop*exp(-cd*zd)
      if (Edz.gt.1.e10)stop
      au = a*rmuu
      Bu = ru*bb*rmuu
      Cu = au+Bu
      as = a*rmus
      Bs = rd*bb*rmus
      Cs = as+Bs
      Bd = bb*rmud
      Fd = (bt-bb)*rmud
      bquad = Cs - Cu
      cquad = Bs*Bu - Cs*Cu
!     sqarg = bquad*bquad - 4.0*cquad
      sqarg = max(0.,bquad*bquad - 4.0*cquad)
      a1 = 0.5*(-bquad + sqrt(sqarg))
      a2 = 0.5*(-bquad - sqrt(sqarg))
      S = -(Bu*Bd + Cu*Fd)
      SEdz = S*Edz
      a2ma1 = a2 - a1
      rM = SEdz/(a1*a2ma1)
      rN = SEdz/(a2*a2ma1)

c     ea2Dmax = exp(a2ma1*Dmax)
c     c1 = (rN-rM)*exp(-a1*Dmax) - (Estop-rM+rN)*ea2Dmax
c    *                             /(1.0-ea2Dmax)
c     c2 = Estop - rM + rN - c1
      c2 = Estop - rM + rN


c     a1arg = a1*zd
c     a1arg = min(a1arg,82.8)
c     Ta1z = exp(a1arg)
      Ta2z = exp(a2*zd)


c      Esz = c1*Ta1z + c2*Ta2z + rM - rN
      Esz = c2*Ta2z + rM - rN
      Esz = max(Esz,0.0)


c      Eutmp = ((a1+Cs)*c1)*Ta1z + ((a2+Cs)*c2)*Ta2z + Cs*rM
c     *             - Cs*rN - Fd*Edz


      if (EUZ_DEFINED.eq. 1) then
        Eutmp = ((a2+Cs)*c2)*Ta2z + Cs*rM
     .               - Cs*rN - Fd*Edz
        Euz = Eutmp/Bu
        Euz = max(Euz,0.0)
      else
        Euz = 0.0
      endif

 
      return
      end
c
c
cc ******************************************************************
c      subroutine upeu 
cc
cc  Computes surface upwelling irradiance.  The approximation used for
cc  upwelling irradiance a depth, i.e., that rmud = rmus, is not
cc  valid at the surface, and a full treatment of diffuse and direct
cc  path lengths is required.  
cc
c#include "gloparam.F.h"
c#include "radpth.h"
c#include "comlte.h"
c#include "ostate.F.h"
c#include "combio.h"
c#include "combio.h"
c      real Edtop(nlt),Estop(nlt)
c      real sfceu(nlt)
c      real bbc(10)
c      data ifst /0/
c      data bbc / 0.002, 0.00071, 0.0032, 0.00071, 0.0029,
c     *           0.0,   0.0,     0.0,    0.0,     0.0/
cc
cc  Constants
c      if (ifst .eq. 0)then
c       rmus = 1.0/0.83            !avg cosine diffuse down
c       rmuu = 1.0/0.4             !avg cosine diffuse up
c       bbw = 0.5    !backscattering to forward scattering ratio
c       rbot = 0.0                 !bottom reflectance
c       rd = 1.5   !these are taken from Ackleson, et al. 1994 (JGR)
c       ru = 3.0
c       Dmax = 500.0    !depth at which Ed = 0
c       do nl = npst,npnd
c        sfceu(nl) = 0.0
c       enddo
c       ifst = 1
c      endif
cc
c      m = indext2
c      k = 1
c      do i = inwst,inwnd
c       do nl = npst,npnd
c        Edtop(nl) = Ed(i,nl)
c        Estop(nl) = Es(i,nl)
c       enddo
c       do nl = npst,npnd
c        a = aw(nl) + acdom(nl)
c        bt = bw(nl)
c        bb = bbw*bw(nl)
c        do n = nnut+1,ntyp-nzoo
c         a = a + P(i,k,m,n)*ac(n-nnut,nl)
c         bt = bt + P(i,k,m,n)*bc(n-nnut,nl)
c         bb = bb + P(i,k,m,n)*bbc(n-nnut)*bc(n-nnut,nl)
c        enddo
cc  Compute Eu
cc   Eu from Ed
c        ad = a*rmud(i)
c        bd = rd*bb*rmud(i)
c        au = a*rmuu
c        bu = ru*bb*rmuu
c        cd = ad+bd
c        cu = au+bu
c        bquad = cd - cu
c        cquad = bd*bu - cd*cu
c        sqarg = bquad*bquad - 4.0*cquad
c        a2 = 0.5*(-bquad - sqrt(sqarg))
c        sfceu1 = (a2+cd)/bu
c        EuEd = Etmp*Ta2z
cc
cc   Eu from Es
c        as = a*rmus
c        bs = rd*bb*rmus
c        cs = as+bs
c        bquad = cs - cu
c        cquad = bs*bu - cs*cu
c        sqarg = bquad*bquad - 4.0*cquad
c        a2 = 0.5*(-bquad - sqrt(sqarg))
c        sfceu2 = (a2+cs)/bu
c        sfceu(nl) = Edtop(nl)*sfceu1 + Estop(nl)*sfceu2
c       enddo
cc       if (i .eq. 334)then
cc        write(6,*)ihr,sfceu(10)
cc       endif
c      enddo
cc
c      return
c      end
