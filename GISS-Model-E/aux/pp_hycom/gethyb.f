      subroutine getdat(flnm,day0,day1,succes)
c
c --- read hybrid model fields (binary) and extract portion of global fields.
c --- then  t r a n s f o r m   to   i s o p y c n i c   fields
c
      use hycom_dimen
      use hycom_arrays, only: ubavg,vbavg,srfhgt,dpmixl,thkice,covice,
     .   depths,u,v,dp,p,temp,saln,th3d,tracer,uflx,vflx,diaflx,
     .   uflxav,vflxav,thmix,tmix,smix,umix,vmix,alloc_hycom_arrays
      use const_proc,only: iorign,jorign,g,rho,dz,spval,timav,cnvert,
     .   itest,jtest

!     use gridlayout
!     use constants,only: cnvert,timav,g,rho,iorign,jorign,spval,
!    .                    latij,scvx,scuy,itest,jtest,onem
!     use variables
!     use fields

      implicit none
      character,intent(IN) :: flnm*(*)
      real    ,intent(OUT) :: day0,day1
      logical ,intent(OUT) :: succes
      real*4 real4(idm,jdm),thbase4,time4,theta4(kdm)
      real theta(kdm),pnew(idm,jdm,kdm+1),curl(idm,jdm),ua,ub,va,vb
      integer*4 lenrec,isiz,jsiz,ksiz,nstep,lgth
      integer ni,nt,nrec
!     integer,parameter :: wtlen=12			! old length of 'what'
      integer,parameter :: wtlen=16			! new length of 'what'
      character(len=wtlen) what
      logical,parameter :: vrbos=.true.
      data ni/14/
c
      lenrec=100
      write (*,'(2a/a,i6)') 'open ',trim(flnm),'  with record length',
     .   lenrec
      open (unit=ni,file=flnm,status='old',access='direct',
     .   recl=100,form='unformatted',err=6)
      read (ni,rec=1,err=6) lgth
      close (unit=ni)
c
c --- f90-compiled standalone micom reports file length in bytes
      if (lgth.gt.4*idm*jdm) lgth=lgth/4
ccc   if (lgth.lt.4*idm*jdm) lgth=lgth*4
      write (*,*) '  reopen file with record length',lgth
     .  ,'   timav=',timav
c
      open (unit=ni,file=flnm,status='old',access='direct',
     .   recl=lgth,form='unformatted',err=6)
c
      nrec=1
      read (ni,rec=nrec) lenrec,isiz,jsiz,ksiz,nstep,time4,
     .     thbase4,theta4
      write (*,*) lenrec,isiz,jsiz,ksiz,nstep,time4,thbase4,theta4
      day1=time4
      if (ksiz.ne.kdm)
     .  stop '(error: number of input layers does not match kdm)'
      do 35 k=1,kdm
 35   theta(k)=theta4(k)+thbase4
c
      write (*,'('' nstep,time:'',i11,f11.2,15x,'' sigma values:''/
     .   (1x,11f7.2))') nstep,day1,(theta(k),k=1,kdm)
c
      do j=1,jdm
       do i=1,idm
        thkice(i,j)=0.
        covice(i,j)=0.
       end do
      end do
c
 4    nrec=nrec+1
      read (ni,rec=nrec,err=5) what,k,real4
      what=adjustr(what)		! right-adjust
      print '(a,i4,2a,i3)','record',nrec,'  contains ',what,k 
 100  format (15x,'store ',a,i3,'   as ',a,2(i2,a))
c
      if (timav) then
c
c --- find the 'av_' substring
       l=max(index(what,'av_',.true.),index(what,'av0',.true.))
       if (l.gt.0) then
        l=l+2
c
        if (index(what,'uav_').gt.0 .and. index(what,'avav').eq.0) then
         print 100,what,k,'u(.,.,',k,')'
!        call flipu(real4)
         call extrct(real4,iorign,jorign,    u(1,1,k))
         if (vrbos) call findmx(iu,u(1,1,k),idm,idm,jdm,'uav   ')
        end if
        if (index(what,'vav_').gt.0 .and. index(what,'avav').eq.0) then
         print 100,what,k,'v(.,.,',k,')'
!        call flipv(real4)
         call extrct(real4,iorign,jorign,    v(1,1,k))
         if (vrbos) call findmx(iv,v(1,1,k),idm,idm,jdm,'vav   ')
        end if
        if (index(what,'dpav_').gt.0) then
         print 100,what,k,'dp(.,.,',k,')'
!        call flipp(real4)
         call extrct(real4,iorign,jorign,   dp(1,1,k))
         if (vrbos) call findmx(ip,dp(1,1,k),idm,idm,jdm,'dpav  ')
        end if
        if (index(what,'temav_').gt.0) then
         print 100,what,k,'temp(.,.,',k,')'
!        call flipp(real4)
         call extrct(real4,iorign,jorign, temp(1,1,k))
         if (vrbos) call findmx(ip,temp(1,1,k),idm,idm,jdm,'temav ')
        end if
        if (index(what,'salav_').gt.0) then
         print 100,what,k,'saln(.,.,',k,')'
!        call flipp(real4)
         call extrct(real4,iorign,jorign, saln(1,1,k))
         if (vrbos) call findmx(ip,saln(1,1,k),idm,idm,jdm,'salav ')
        end if
        if (index(what,'th3dav_').gt.0) then
         print 100,what,k,'th3d(.,.,',k,')'
!        call flipp(real4)
         call extrct(real4,iorign,jorign, th3d(1,1,k))
         if (vrbos) call findmx(ip,th3d(1,1,k),idm,idm,jdm,'th3dav')
        end if
        if (index(what,'ubavav_').gt.0) then
         print 100,what,k,'ubavg'
!        call flipu(real4)
         call extrct(real4,iorign,jorign, ubavg)
         if (vrbos) call findmx(iu,ubavg,idm,idm,jdm,'ubavav')
        end if
        if (index(what,'vbavav_').gt.0) then
         print 100,what,k,'vbavg'
!        call flipv(real4)
         call extrct(real4,iorign,jorign, vbavg)
         if (vrbos) call findmx(iv,vbavg,idm,idm,jdm,'vbavav')
        end if
        if (what(l-6:l).eq.'sfhtav_' .or.
     .      what(l-6:l).eq.'srfhav_') then
         print 100,what,k,'srfhgt'
!        call flipp(real4)
         call extrct(real4,iorign,jorign, srfhgt)
         if (vrbos) call findmx(ip,srfhgt,idm,idm,jdm,'sfhtav')
        end if
        if (index(what,'dpmxav_').gt.0) then
         print 100,what,k,'dpmixl'
!        call flipp(real4)
         call extrct(real4,iorign,jorign,dpmixl)
         if (vrbos) call findmx(ip,dpmixl,idm,idm,jdm,'dpmxav')
        end if
        if (index(what,'oiceav_').gt.0) then
         print 100,what,k,'covice'
!        call flipp(real4)
         call extrct(real4,iorign,jorign,covice)
         if (vrbos) call findmx(ip,covice,idm,idm,jdm,'oiceav')
        end if
       if (index(what,'ziceav_').gt.0) then
        print 100,what,k,'thkice'
!       call flipp(real4)
        call extrct(real4,iorign,jorign,thkice)
        if (vrbos) call findmx(ip,thkice,idm,idm,jdm,'ziceav')
       end if
c
!       do nt=1,ntrcr
!        if (index(what,trcid(nt)(1:8)).gt.0) then
!         print 100,what,k,'tracer(.,.,',k,',',nt,')'
!!        call flipp(real4)
!         call extrct(real4,iorign,jorign, tracer(1,1,k,nt))
!         if (vrbos) call findmx(ip,tracer(1,1,k,nt),idm,idm,jdm,
!    .     'tracer')
!        end if
!       end do
c
       end if
      else if (max(index(what,'av_'),index(what,'av0')).eq.0) then

       if (index(what,'  u').gt.0 .and. index(what,'  ub').eq.0) then
        print 100,what,k,'u(.,.,',k,')'
!       call flipu(real4)
        call extrct(real4,iorign,jorign,    u(1,1,k))
        if (vrbos) call findmx(iu,u(1,1,k),idm,idm,jdm,'u     ')
       end if
       if (index(what,'  v').gt.0 .and. index(what,'  vb').eq.0) then
        print 100,what,k,'v(.,.,',k,')'
!       call flipv(real4)
        call extrct(real4,iorign,jorign,    v(1,1,k))
        if (vrbos) call findmx(iv,v(1,1,k),idm,idm,jdm,'v     ')
       end if
       if (index(what,'   dp').gt.0) then
        print 100,what,k,'dp(.,.,',k,')'
!       call flipp(real4)
        call extrct(real4,iorign,jorign,   dp(1,1,k))
        call findmx(ip,dp(1,1,k),idm,idm,jdm,'dp    ')
       end if
       if (index(what,' temp').gt.0) then
        print 100,what,k,'temp(.,.,',k,')'
!       call flipp(real4)
        call extrct(real4,iorign,jorign, temp(1,1,k))
        if (vrbos) call findmx(ip,temp(1,1,k),idm,idm,jdm,'temp  ')
       end if
       if (index(what,' saln').gt.0) then
        print 100,what,k,'saln(.,.,',k,')'
!       call flipp(real4)
        call extrct(real4,iorign,jorign, saln(1,1,k))
        if (vrbos) call findmx(ip,temp(1,1,k),idm,idm,jdm,'saln  ')
       end if
       if (index(what,' th3d').gt.0) then
        print 100,what,k,'th3d(.,.,',k,')'
!       call flipp(real4)
        call extrct(real4,iorign,jorign, th3d(1,1,k))
        if (vrbos) call findmx(ip,th3d(1,1,k),idm,idm,jdm,'th3d  ')
       end if
       if (index(what,'ubavg').gt.0) then
        print 100,what,k,'ubavg'
!       call flipu(real4)
        call extrct(real4,iorign,jorign, ubavg)
        if (vrbos) call findmx(iu,ubavg,idm,idm,jdm,'ubavg ')
       end if
       if (index(what,'vbavg').gt.0) then
        print 100,what,k,'vbavg'
!       call flipv(real4)
        call extrct(real4,iorign,jorign, vbavg)
        if (vrbos) call findmx(iv,vbavg,idm,idm,jdm,'vbavg ')
       end if
       if (index(what,'srfhgt').gt.0) then
        print 100,what,k,'srfhgt'
!       call flipp(real4)
        call extrct(real4,iorign,jorign, srfhgt)
        if (vrbos) call findmx(ip,srfhgt,idm,idm,jdm,'srfhgt ')
       end if
       if (index(what,'mix_dpth').gt.0) then
        print 100,what,k,'dpmixl'
!       call flipp(real4)
        call extrct(real4,iorign,jorign,dpmixl)
        if (vrbos) call findmx(ip,dpmixl,idm,idm,jdm,'dpmixl ')
       end if
       if (index(what,'icethik').gt.0) then
        print 100,what,k,'thkice'
!       call flipp(real4)
        call extrct(real4,iorign,jorign,thkice)
        if (vrbos) call findmx(ip,thkice,idm,idm,jdm,'thkice')
       end if
       if (index(what,'icecover').gt.0) then
        print 100,what,k,'covice'
!       call flipp(real4)
        call extrct(real4,iorign,jorign,covice)
        if (vrbos) call findmx(ip,covice,idm,idm,jdm,'covice')
       end if
      end if			! end-of-month snapshot
c
c --- find the 'av_' substring
      l=max(index(what,'av_',.true.),index(what,'av0',.true.))
      if (l.eq.0) go to 7
      l=l+2
c
! --- time-integrated mass fluxes
      if (index(what,'uflxav_').gt.0) then
       print 100,what,k,'uflx(.,.,',k,')'
!      call flipu(real4)
       call extrct(real4,iorign,jorign,uflx(1,1,k))
       if (vrbos) call findmx(iu,uflx(1,1,k),idm,idm,jdm,
     .    'uflxav')
      end if
      if (index(what,'vflxav_').gt.0) then
       print 100,what,k,'vflx(.,.,',k,')'
!      call flipv(real4)
       call extrct(real4,iorign,jorign,vflx(1,1,k))
       if (vrbos) call findmx(iv,vflx(1,1,k),idm,idm,jdm,
     .    'vflxav')
      end if
      if (index(what,'diaflx_').gt.0) then
       print 100,what,k,'diaflx(.,.,',k,')'
!      call flipp(real4)
       call extrct(real4,iorign,jorign,diaflx(1,1,k))
       if (vrbos) call findmx(ip,diaflx(1,1,k),idm,idm,jdm,
     .    'diaflx')
      end if

! --- surface forcing fields
 7    continue
!     if (max(index(what,'eminp'),index(what,' emp')).gt.0) then
!      print 100,what,k,'surflx1'
!      call flipp(real4)
!      call extrct(real4,iorign,jorign,surflx1)
!      if (vrbos) call findmx(ip,surflx1,idm,idm,jdm,'surflx1')
!      fluxid1=what
!     end if
!     if (max(index(what,'htflx'),index(what,'hflx'),
!    .        index(what,'radfl')).gt.0) then
!      print 100,what,k,'surflx2'
!      call flipp(real4)
!      call extrct(real4,iorign,jorign,surflx2)
!      if (vrbos) call findmx(ip,surflx2,idm,idm,jdm,'surflx2')
!      fluxid2=what
!     end if
!     if (max(index(what,'sflx'),index(what,'sfx'),
!    .        index(what,'snsib')).gt.0) then
!      print 100,what,k,'surflx3'
!      call flipp(real4)
!      call extrct(real4,iorign,jorign,surflx3)
!      if (vrbos) call findmx(ip,surflx3,idm,idm,jdm,'surflx3')
!      fluxid3=what
!     end if
!     if (max(index(what,'brine'),index(what,'brn'),
!    .        index(what,'latnt')).gt.0) then
!      print 100,what,k,'surflx4'
!      call flipp(real4)
!      call extrct(real4,iorign,jorign,surflx4)
!      if (vrbos) call findmx(ip,surflx4,idm,idm,jdm,'surflx4')
!      fluxid4=what
!     end if
!     if (max(index(what,'Tau_x'),index(what,'taux')).gt.0) then
!      print 100,what,k,'surflx5'
!      call flipp(real4)
!      call extrct(real4,iorign,jorign,surflx5)
!      if (vrbos) call findmx(ip,surflx5,idm,idm,jdm,'surflx5')
!      fluxid5=what
!     end if
!     if (max(index(what,'Tau_y'),index(what,'tauy')).gt.0) then
!      print 100,what,k,'surflx6'
!      call flipp(real4)
!      call extrct(real4,iorign,jorign,surflx6)
!      if (vrbos) call findmx(ip,surflx6,idm,idm,jdm,'surflx6')
!      fluxid6=what
!     end if
      go to 4
c
 5    if (timav) then
c --- find substring containing length of averaging interval
       print *,'get averaging interval from string ',what(wtlen-3:wtlen)
       if (wtlen.eq.16) read (what(wtlen-3:wtlen),'(f4.0)') day0
       if (wtlen.eq.12) read (what(wtlen-2:wtlen),'(f3.0)') day0
       day0=day1-day0
       write (*,'(2(a,f9.1))') 'input fields are averaged over days',
     .  day0,' --',day1
      else
        day0=day1
      end if
c
      close (unit=ni)
c
c --- fix physical dimensions
c
!     fluxid1='eminp (m/yr)'
c
      do 20 j=1,jdm
c
      do 21 i=1,idm
      ubavg(i,j)=ubavg(i,j)*100. * iu(i,j)		!  convert to cm/sec
      vbavg(i,j)=vbavg(i,j)*100. * iv(i,j)		!  convert to cm/sec
      p(i,j,1)=0.
      if (ip(i,j).gt.0) then
        dpmixl(i,j,1)=dpmixl(i,j,1)/(g*rho)		!  convert to meters
!       srfhgt(i,j)=srfhgt(i,j)*100./g			!  convert to cm
        srfhgt(i,j)=srfhgt(i,j)*100.			!  convert to cm
!       surflx1(i,j)=surflx1(i,j)*365.*86400.		!  eminp m/s -> m/yr
!       covice(i,j)=covice(i,j)*100.			!  convert to percent
      else
        dpmixl(i,j,1)=spval
        srfhgt(i,j)=spval
        covice(i,j)=spval
!       surflx1(i,j)=spval
!       surflx2(i,j)=spval
!       surflx3(i,j)=spval
!       surflx4(i,j)=spval
      end if
 21   continue
c
c --- combine barotropic and baroclinic velocities
c
      do 20 k=1,kdm
      do 20 i=1,idm
      u(i,j,k)=u(i,j,k)*100./dz*iu(i,j) + ubavg(i,j)	!  convert to cm/sec
      v(i,j,k)=v(i,j,k)*100./dz*iv(i,j) + vbavg(i,j)	!  convert to cm/sec
      uflx(i,j,k)=uflx(i,j,k)*1.e-6/(g*rho*dz**3)	!  convert to Sv*intvl
      vflx(i,j,k)=vflx(i,j,k)*1.e-6/(g*rho*dz**3)	!  convert to Sv*intvl
      if (ip(i,j).gt.0) then
        dp(i,j,k)=dp(i,j,k)/(g*rho*dz)			!  convert to meters
!       th3d(i,j,k)=(th3d(i,j,k)+thbase4)*1000./rho	!  add thbase
      else
        dp(i,j,k)=spval
        temp(i,j,k)=spval
        saln(i,j,k)=spval
        th3d(i,j,k)=spval
      end if
 20   continue
c
      call findmx(ip,dp(1,1,2),idm,idm-1,jdm-1,'layer 2 thknss')
c
      do 39 j=1,jdm
c
c --- save mixed layer fields for future use
      do 38 i=1,idm
      thmix(i,j)=th3d(i,j,1) 
      tmix(i,j)=temp(i,j,1)
      smix(i,j)=saln(i,j,1)
      umix(i,j)=u(i,j,1)
      vmix(i,j)=v(i,j,1)
      p(i,j,2)=dp(i,j,1)
 38   pnew(i,j,1)=0.
c
      do 39 k=2,kdm
      do 39 i=1,idm
      th3d(i,j,k)=max(th3d(i,j,k-1),th3d(i,j,k))
 39   p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
c
c --- transform hybrid fields to isopycnic coord.
c
      if (itest.gt.0 .and. jtest.gt.0) print '(2i5,a/(10f8.1))',
     .  itest,jtest,' (gethyb) pressure:',(p(itest,jtest,k),k=1,kdm+1)

      if (cnvert) then
        write (*,*) 'now transforming input fields to isopyc.coord.'
        call restep(u,v,temp,saln,tracer,th3d,p,u(1,1,kdm+1),
     .    v(1,1,kdm+1),temp(1,1,kdm+1),saln(1,1,kdm+1),
     .    tracer(1,1,kdm+1,1),th3d(1,1,kdm+1),pnew,theta,kdm,kdm)
c
        do 36 j=1,jdm
        do 36 k=1,kdm
        do 36 i=1,idm
        u(i,j,k)=u(i,j,kdm+k)
        v(i,j,k)=v(i,j,kdm+k)
        temp(i,j,k)=temp(i,j,kdm+k)
        saln(i,j,k)=saln(i,j,kdm+k)
        th3d(i,j,k)=th3d(i,j,kdm+k)
        if (ntrcr.gt.0) tracer(i,j,k,:)=tracer(i,j,kdm+k,:)
        pnew(i,j,k+1)=p(i,j,k+1)
        if (ip(i,j).gt.0) dp(i,j,k)=pnew(i,j,k+1)-pnew(i,j,k)
 36     continue
      end if				! cnvert
c
!     write (*,'(a)') 'shown below: mixed layer density'
!     call zebra(thmix,idm,idm-1,jdm-1)
!     write (*,'(a)') 'shown below: sea surface height (cm)'
!     call zebra(srfhgt,idm,idm-1,jdm-1)
      write (*,'(a)') 'shown below: SST'
      call zebra(temp(1,1,2),idm,idm-1,jdm-1)
!     write (*,'(a)') 'shown below: layer 1 thickness (m)'
!     call zebra(dp,idm,idm-1,jdm-1)
!     write (*,'(a)') 'shown below: ice cover (%)'
!     call zebra(covice,idm,idm-1,jdm-1)
c
      succes=.true.
      return
c
 6    succes=.false.
      return
      end
