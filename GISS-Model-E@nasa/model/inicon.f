#include "rundeck_opts.h"
      subroutine inicon
c
c --- hycom version 0.9

      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS

!      USE FLUXES, only : e0,prec,evapor,flowo,eflowo,dmua,dmva
!     .      ,erunosi,runosi,runpsi,dmui,dmvi,dmsi,dhsi,dssi
!     .      ,gtemp,sss,mlhc,gtempr
!      USE MODEL_COM, only : focean
      USE HYCOM_ARRAYS_GLOB
      USE KPRF_ARRAYS, only : akpar
      use filemanager, only : findunit

      implicit none
c
      integer i,j,k,l,m,n,mm,ia,ja,iu1,iu2,iu3,iu4
      include 'state_eqn.h'
c
      integer totlj(jdm,kdm-1),totl(kdm-1),iz,jz,ni
      character text*24,preambl(5)*79
      real cold,temavg,vol,sst,sigocn,sigstar,sofsig
      real*4 real4(idm,jdm)
      external sigocn,sigstar,sofsig
      character title*80

!!! not sure why I added this line ... IA
!!!      asst(:,:) = 0.d0
c
c --- set minimum salinity for each isopycnic layer
      do k=1,kk
      salmin(k)=sofsig(theta(k),-3.0)
      end do
      write(*,'(a/(10f7.2))') 'minimum salinities:',salmin
      write(*,'('' theta(k)     :'',10f6.2/(15x,10f6.2))')
     .   (theta(k),k=1,kk)
      write(*,'('' dplist(k)     :'',10f6.0/(15x,10f6.0))')
     .   (dplist(k),k=1,kk)
c
      if (nstep0.eq.0) then                ! start from Levitus
        !!call geopar
        delt1=baclin
c
c --- read 3-d temperature
c
        write (*,'(2a)') 'get initial temperature from  ',flnmint
        call findunit(iu1)
        open(unit=iu1,file=flnmint,form='formatted',status='old',
     .     action='read')
        read (iu1,'(a79)') (preambl(n),n=1,5)
        write(*,'(a79)') (preambl(n),n=1,5)
        do k=1,kk
        read (iu1,'(10f8.4)') ((temp(i,j,k),i=1,idm),j=1,jdm)
        enddo
        close (iu1)
        write (*,100) 'temp field read, layers, 1 -',kk
        call zebra(temp(1,1,1),idm,ii1,jj)
c
c --- read 3-d salinity
c
        write (*,'(2a)') 'get initial salinity from  ',flnmins
        call findunit(iu2)
        open(unit=iu2,file=flnmins,form='formatted',status='old',
     .     action='read')
        read (iu2,'(a79)') (preambl(n),n=1,5)
        write(*,'(a79)') (preambl(n),n=1,5)
        do k=1,kk
        read (iu2,'(10f8.4)') ((saln(i,j,k),i=1,idm),j=1,jdm)
        enddo
        close (iu2)
        write (*,100) 'saln field read, layers, 1 -',kk
        call zebra(saln(1,1,1),idm,ii1,jj)
 100    format (a,i4)
c
c --- read interface pressure
c
        write (*,'(2a)') 'get initial pressure from  ',flnminp
        call findunit(iu3)
        open(unit=iu3,file=flnminp,form='formatted',status='old',
     .     action='read')
        read (iu3,'(a79)') (preambl(n),n=1,5)
        write (*,'(a79)') (preambl(n),n=1,5)
        do k=1,kk
        read (iu3,'(10f8.4)') ((p(i,j,k+1),i=1,idm),j=1,jdm)
        enddo
        close (iu3)
        write (*,100) 'pres field read, levels 2 -',kk+1
        call zebra(p(1,1,kk+1),idm,ii1,jj)
c
        do i=1,ii
        do j=1,jj
        if (depths(i,j).le.0.) then	! set "huge" on land points
          temp(i,j,:)=huge
          saln(i,j,:)=huge
             p(i,j,:)=huge
        end if
        end do
        end do

        do 10 j=1,jj
        do 10 l=1,isp(j)
c
        do 15 i=ifp(j,l),ilp(j,l)
 15     p(i,j,1)=0.
c
        do 9 k=1,kk
        do 9 i=ifp(j,l),ilp(j,l)
        p(i,j,k+1)=max(p(i,j,k),min(depths(i,j),p(i,j,k+1))*onem)
        dp(i,j,k)=p(i,j,k+1)-p(i,j,k)
        th3d(i,j,k)=sigocn(temp(i,j,k),saln(i,j,k))
 9      continue
c
        do 17 i=ifp(j,l),ilp(j,l)
 17     pbot(i,j)=p(i,j,kk+1)
c
 10     continue
c
cdiag do k=1,kk,3
cdiag write (text,'(''intf.pressure (m), k='',i3)') k+1
cdiag call prtmsk(ip,p(1,1,k+1),util1,idm,ii1,jj,0.,1./onem,text)
cdiag end do
c
c     call convec(1,1,0,0,1,1)
c
      do 11 j=1,jj
c
      do 12 k=1,kk-1
 12   totlj(j,k)=0
c
      do 11 k=1,kk
      do 11 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
css      tracer(i,j,k)=0.                 ! moved to hycom.f temperarily
      dp  (i,j,k+kk)=dp  (i,j,k)
      th3d(i,j,k+kk)=th3d(i,j,k)
      temp(i,j,k+kk)=temp(i,j,k)
      saln(i,j,k+kk)=saln(i,j,k)
      if (kappa) then
        thstar(i,j,k)=sigstar(temp(i,j,k),saln(i,j,k),p(i,j,k))
      else
        thstar(i,j,k)=th3d(i,j,k)
      end if
c
      if (itest.gt.0.and.jtest.gt.0) then
        if (i.eq.itest.and.j.eq.jtest)
     . write (*,'(2i4,i3,a,3f7.2,2x,2f7.3,f8.1)')
     .  i,j,k,' dens,thstar,kappa,t,s,p=',th3d(i,j,k),thstar(i,j,k)
     .   ,thstar(i,j,k)-th3d(i,j,k)
     .   ,temp(i,j,k),saln(i,j,k),p(i,j,k+1)/onem
      else
        if (i.eq.equatn.and.j.eq.3)
     . write (*,'(2i4,i3,a,3f7.2,2x,2f7.3,f8.1)')
     .  i,j,k,' dens,thstar,kappa,t,s,p=',th3d(i,j,k),thstar(i,j,k)
     .   ,thstar(i,j,k)-th3d(i,j,k)
     .   ,temp(i,j,k),saln(i,j,k),p(i,j,k+1)/onem
      endif
c
      if (k.gt.1) then
        if (thstar(i,j,k).lt.thstar(i,j,k-1))
     .    totlj(j,k-1)=totlj(j,k-1)+1
      endif
 11   continue
c
      write(*,'(a,20i12)') 'ijlist ',((ijlist(i,j),i=30,32),j=4,5)
c
      do k=1,kk
#ifdef HYCOM_UNFINISHED
      write (*,'(i5,a,2i5,a/7x,7(i3,3x),3x,7(i3,3x)/
     .  (/(7(i4,7f6.0,3x,7f6.0/))))')
     .  k,' i,j=',itest,jtest,' input data (t,s,p,depth)'
     . ,        (j,j=jtest-3,jtest+3),(j,j=jtest-3,jtest+3)
     . ,(i,(temp(i,j,k),j=jtest-3,jtest+3)
     . ,   (saln(i,j,k),j=jtest-3,jtest+3),i=itest-3,itest+3)
     . ,(i,(p(i,j,k+1)/onem,j=jtest-3,jtest+3)
     . ,    (depths(i,j),j=jtest-3,jtest+3),i=itest-3,itest+3)
#endif
c
c     if (itest.gt.0.and.jtest.gt.0)
c    . write (*,'(2i4,a,i2/(5(5f7.1,3x,5f7.1/)))')
c    . itest,jtest,' initial t,s,p,depth at k',k,
c    . ((temp(i,j,k),j=jtest-2,jtest+2)
c    . ,(saln(i,j,k),j=jtest-2,jtest+2),i=itest-2,itest+2)
c    . ,((p(i,j,k+1)/onem,j=jtest-2,jtest+2)
c    . , (depths(i,j),j=jtest-2,jtest+2),i=itest-2,itest+2)
c
      enddo
c
      do 18 k=1,kk-1
      totl(k)=0
      do 18 j=1,jj
 18   totl(k)=totl(k)+totlj(j,k)
      write (*,'(a/(10i7))') 'static instability count by layer:',
     .  totl
c
      do 50 j=1,jj
      do 50 l=1,isp(j)
      do 50 i=ifp(j,l),ilp(j,l)
      montg(i,j,1)=0.
      omlhc(i,j)=spcifh*p(i,j,2)/(onem *thref)           ! J/(m2*C)
c
      do 52 k=1,kk-1
 52   montg(i,j,k+1)=montg(i,j,k)-p(i,j,k+1)*(thstar(i,j,k+1)-
     .                                        thstar(i,j,k  ))*thref**2
      thkk(i,j)=thstar(i,j,kk)
 50   psikk(i,j)=montg(i,j,kk)
c
      else                                !  nstep0 > 0
c
c --- start from restart file prescribed
      write (*,'(2a)') 'get initial condition from restart file'
c
      write (*,111) nstep0,time0
 111  format (9x,'chk time step in restart file -',i9,5x,' day ',f9.2)
c
      delt1=baclin+baclin
css   call newbot
c
      do 21 j=1,jj
      do 21 k=1,kk
      do 21 l=1,isp(j)
      do 21 i=ifp(j,l),ilp(j,l)
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
      if (kappa) then
        thstar(i,j,k)=sigstar(temp(i,j,k),saln(i,j,k),p(i,j,k))
      else
        thstar(i,j,k)=th3d(i,j,k)
      end if
 21   continue
c
      end if                                !  nstep0 > 0  or  = 0
c
      call dpthuv
c
      do 16 m=1,2
      mm=(m-1)*kk
c
      do 19 j=1,jj
      do 19 k=1,kk
      do 19 l=1,isp(j)
      do 19 i=ifp(j,l),ilp(j,l)
 19   p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
c
      call dpudpv(mm)
 16   continue
c
c     print *,' focean'
c     call zebra(focean,iia,iia,jja)
c
c     print *,'chk ini. gtemp at nstep=',nstep0
c     call zebra(asst,iia,iia,jja)
c
c     print *,'chk ini. sss at nstep=',nstep0
c     call zebra(sss,iia,iia,jja)
c

      if (itest.gt.0.and.jtest.gt.0) then
        i=itest
        j=jtest
      else
        i=equatn; j=3
        i=154; j=198
      endif

       call pr_9x9(temp(:,:,1),ii,jj,i,j,0.,1.,'sst ini')

      write (*,'(a,2i4,4f8.2)') ' sig=',i,j,temp(i,j,1),saln(i,j,1),
     .   sigocn(temp(i,j,1),saln(i,j,1))
      write (*,103) nstep,i,j,
     .  '  init.profile  temp    saln  thstar   thkns    dpth   montg',
     .  (k,temp(i,j,k),saln(i,j,k),thstar(i,j,k),dp(i,j,k)/onem,
     .  p(i,j,k+1)/onem,montg(i,j,k)/g,k=1,kk)
c
      if (jerlv0.eq.0) then
c ---   read-in monthly kpar file
        write(*,*) 'opening kpar '
        real4=0.
        call findunit(iu4)
        open(iu4,file='kpar',form='unformatted',status='old')
        do k=1,12
        write(*,*) 'reading kpar mo=',k
        read(iu4) title,real4
        write(*,*)'title=',title(1:60)
        akpar(:,:,k)=real4(:,:)
        enddo
        close(iu4)
      endif
c
c     call zebra(akpar,idm,idm,jdm)
 103  format (i7,2i4,a/(24x,i3,2f8.2,f8.2,2f8.1,f8.3))
c
      return
      end
c> Revision history:
c>
c> Mar. 2000 - conversion to SI units
c> Aug. 2000 - added diagnostic count of static instabilities
c> Apr. 2001 - eliminated stmt_funcs.h
c> Sep. 2005 - added EQ refinement
c> JAN. 2008 - no need for EQ refinement - it is done in pre-processing
