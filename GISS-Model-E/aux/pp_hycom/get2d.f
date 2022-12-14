      program get2d
c  (1) read in monthly z-output
c  (2) write output in 2d of ij, jk, ik on every z & interval of lat0/lon0
c
c     use hycom_arrays, only: depths,srfhgt,dpmixl,oice,u,v,dp,p
c    .   ,temp,saln,th3d,tracer,alloc_hycom_arrays
      use hycom_dimen
      use const_proc
c     use hycom_o2a
      implicit none
c
      real :: day1
      integer, allocatable :: im(:,:)
c
      integer mo,ny1,ny2,dcd,status,ia,k00,ny,m,nt
      integer, parameter:: nrec=8,lat0=10,lon0=10
      real*8 :: area,avgbot
      character flnmin*80,runid*20,flnmout1*80,flnmout2*80,flnmout3*80
     .         ,ttl*80,ttl1*80,ttl2*80,ttl3*80
      logical timav,cnvert,monave_convert,solo_convert
      character*26 ttlt(ntrcr),ttl0(nrec)
      real :: avg2d(iia,jja),avg3d(iia,jja,k33,nrec-3)
     .  ,sshij(iia,jja),dpmixij(iia,jja),iceij(iia,jja)
     .  ,tz(iia,jja,k33),sz(iia,jja,k33),uz(iia,jja,k33)
     .  ,vz(iia,jja,k33),rz(iia,jja,k33),trz(iia,jja,k33,ntrcr)
     .  ,a2d(iia,jja,12),a3d(iia,jja,k33,12)
     .  ,worka(iia,jja),worko(idm,jdm),depthij(iia,jja),xlat(jja)
     .  ,xlon(iia)
c
      do j=1,jja
      xlat(j)=-90.+(j-0.5)*180./jja
      enddo
c
      do i=1,iia
      xlon(i)=(i-0.5)*360./iia
      enddo
c
      do nt=1,ntrcr
      write(ttlt(nt),'(a,i2.2)') 'tracer No.',nt
      enddo
c
      read(*,*) runid,ny1,ny2,monave_convert,solo_convert
      if (monave_convert=='') monave_convert=.true.
      if (solo_convert=='')   solo_convert=.true.
      if ((ny1==ny2 .and. monave_convert==.true.) .or.
     .    (ny1==ny2 .and. solo_convert==.false. ) ) then
        monave_convert=.false.
        solo_convert=.true.
        write(*,*) 'Note: force monave_convert= F & solo_convert = T for
     . single year'
      endif
      write(*,*) 'get2d ',runid,' from yr ',ny1,' to ',ny2,
     . ' monave_convert=',monave_convert,' solo_convert=',solo_convert

      write(*,'(a,i2)') 'number of tracers =',ntrcr
c
c     call alloc_hycom_arrays
      call alloc_hycom_dimen

      dcd=ny1/10
      if (runid(1:1).eq.' ') stop 'empty runid'
      if (dcd.lt.180 .or. dcd.gt.330) then
        print *,' wrong decade=',dcd
        stop 'wrong decade'
      endif
c
      if (monave_convert) then
        do 152 mo=mo1,mo2+1
        write(flnmin,'(a,i4.4,a,i4.4,2a)') amon(mo),ny1,'_',ny2,'.zout'
     .         ,trim(runid)
        write(flnmout1,'(a,i4.4,a,i4.4,2a)')amon(mo),ny1,'_',ny2,'.oij'
     .         ,trim(runid)
        write(flnmout2,'(a,i4.4,a,i4.4,2a)')amon(mo),ny1,'_',ny2,'.oik'
     .         ,trim(runid)
        write(flnmout3,'(a,i4.4,a,i4.4,2a)')amon(mo),ny1,'_',ny2,'.ojk'
     .         ,trim(runid)
c     write(flnmin,'(5a,i3,a,i3,2a,i4.4,2a)')trim(path1),trim(runid)
c    . ,'/zout',trim(runid),'_',dcd,'0_',dcd,'9/',amon(mo),ny
c    . ,'.zout',trim(runid)

      write(ttl1,'(i3,a,i3,3x,2(1x,a),i4,a,i4)') 
     .         iia,'x',jja,trim(runid),amon(mo),ny1,'_',ny2
      write(ttl2,'(i3,a,i2,2(1x,a),i4,a,i4)') 
     .         iia,'x',k33,trim(runid),amon(mo),ny1,'_',ny2
      write(ttl3,'(i3,a,i2,2(1x,a),i4,a,i4)') 
     .         jja,'x',k33,trim(runid),amon(mo),ny1,'_',ny2

      open(51,file=flnmin,form='unformatted',status='old')
      open(52,file=flnmout1,form='unformatted',status='unknown')
      open(53,file=flnmout2,form='unformatted',status='unknown')
      open(54,file=flnmout3,form='unformatted',status='unknown')

      write(*,'(2a)') 'flnmin=',trim(flnmin)
      write(*,'(2a,1x,a,1x,a)') 'flnmout=',trim(flnmout1),trim(flnmout2)
     .                                    ,trim(flnmout3)
c
c --- read data on z-level
c
      read (51) ttl, sshij
      write(52) ttl, sshij
      read (51) ttl, dpmixij
      write(52) ttl, dpmixij
      read (51) ttl, iceij
      write(52) ttl, iceij
c     
      read(51) ttl, uz
      ttl0(4)=ttl(1:26)
      read(51) ttl, vz
      ttl0(5)=ttl(1:26)
      read(51) ttl, tz
      ttl0(6)=ttl(1:26)
      read(51) ttl, sz
      ttl0(7)=ttl(1:26)
      read(51) ttl, rz
      ttl0(8)=ttl(1:26)
      do nt=1,ntrcr
      read(51) ttl,trz(:,:,:,nt)
      ttlt(nt)=ttl(1:26)
      enddo

      if (diag) then
        i=iatest
        j=jatest
        write(*,'(a,2i4)') 'chk              t          s           u   
     .        v      trc  at ',i,j
        do k=1,k33
        write(*,'(i2,f7.0,5f12.4)')k,z33(k),tz(i,j,k),sz(i,j,k)
     .     ,uz(i,j,k),vz(i,j,k),trz(i,j,k,ntrcr)
        enddo
      endif

      do k=1,k33
      write(ttl,'(2a,i2,1x,a)') ttl0(4),' k=',k,trim(ttl1)
      write(52) ttl, uz(:,:,k)
      write(ttl,'(2a,i2,1x,a)') ttl0(5),' k=',k,trim(ttl1)
      write(52) ttl, vz(:,:,k)
      write(ttl,'(2a,i2,1x,a)') ttl0(6),' k=',k,trim(ttl1)
      write(52) ttl, tz(:,:,k)
      write(ttl,'(2a,i2,1x,a)') ttl0(7),' k=',k,trim(ttl1)
      write(52) ttl, sz(:,:,k)
      write(ttl,'(2a,i2,1x,a)') ttl0(8),' k=',k,trim(ttl1)
      write(52) ttl, rz(:,:,k)
      do nt=1,ntrcr
      write(ttl,'(2a,i2,1x,a)') ttlt(nt),' k=',k,trim(ttl1)
      write(52) ttl,trz(:,:,k,nt)
      enddo
      enddo
c
      do j=1,jja,lat0
      write(ttl,'(2a,f5.1,1x,a)')ttl0(4),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,uz(:,j,:),xlon,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(5),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,vz(:,j,:),xlon,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(6),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,tz(:,j,:),xlon,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(7),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,sz(:,j,:),xlon,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(8),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,rz(:,j,:),xlon,z33
      do nt=1,ntrcr
      write(ttl,'(2a,f5.1,1x,a)')ttlt(nt),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,trz(:,j,:,nt),xlon,z33
      enddo
      enddo
c
      do i=1,iia,lon0
      write(ttl,'(2a,f5.1,1x,a)')ttl0(4),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,uz(i,:,:),xlat,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(5),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,vz(i,:,:),xlat,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(6),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,tz(i,:,:),xlat,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(7),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,sz(i,:,:),xlat,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(8),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,rz(i,:,:),xlat,z33
      do nt=1,ntrcr
      write(ttl,'(2a,f5.1,1x,a)')ttlt(nt),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,trz(i,:,:,nt),xlat,z33
      enddo
      enddo
      close(51)
      close(52)
      close(53)
      close(54)
 152  continue   ! mo=mo1,mo2
      endif      ! monave_convert

      if(solo_convert) then
      do 151 ny=ny1,ny2
      do 151 mo=mo1,mo2+1
c     write(flnmin,'(5a,i3,a,i3,2a,i4.4,2a)')trim(path1),trim(runid)
c    . ,'/zout',trim(runid),'_',dcd,'0_',dcd,'9/',amon(mo),ny
c    . ,'.zout',trim(runid)

        write(flnmin,'(a,i4.4,2a)') amon(mo),ny,'.zout',trim(runid)
        write(flnmout1,'(a,i4.4,2a)')amon(mo),ny,'.oij',trim(runid)
        write(flnmout2,'(a,i4.4,2a)')amon(mo),ny,'.oik',trim(runid)
        write(flnmout3,'(a,i4.4,2a)')amon(mo),ny,'.ojk',trim(runid)

      write(ttl1,'(i3,a,i3,3x,2(1x,a),i4)') 
     .         iia,'x',jja,trim(runid),amon(mo),ny
      write(ttl2,'(i3,a,i2,2(1x,a),i4)') 
     .         iia,'x',k33,trim(runid),amon(mo),ny
      write(ttl3,'(i3,a,i2,2(1x,a),i4)') 
     .         jja,'x',k33,trim(runid),amon(mo),ny

      open(51,file=flnmin,form='unformatted',status='old')
      open(52,file=flnmout1,form='unformatted',status='unknown')
      open(53,file=flnmout2,form='unformatted',status='unknown')
      open(54,file=flnmout3,form='unformatted',status='unknown')

      write(*,'(2a)') 'flnmin=',trim(flnmin)
      write(*,'(2a,1x,a,1x,a)') 'flnmout=',trim(flnmout1),trim(flnmout2)
     .                                    ,trim(flnmout3)
c
c --- read data on z-level
c
      read (51) ttl, sshij
      write(52) ttl, sshij
      read (51) ttl, dpmixij
      write(52) ttl, dpmixij
      read (51) ttl, iceij
      write(52) ttl, iceij
c     
      read(51) ttl, uz
      ttl0(4)=ttl(1:26)
      read(51) ttl, vz
      ttl0(5)=ttl(1:26)
      read(51) ttl, tz
      ttl0(6)=ttl(1:26)
      read(51) ttl, sz
      ttl0(7)=ttl(1:26)
      read(51) ttl, rz
      ttl0(8)=ttl(1:26)
      do nt=1,ntrcr
      read(51) ttl,trz(:,:,:,nt)
      ttlt(nt)=ttl(1:26)
      enddo

      if (diag) then
        i=iatest
        j=jatest
        write(*,'(a,2i4)') 'chk              t          s           u   
     .        v     trc  at ',i,j
        do k=1,k33
        write(*,'(i2,f7.0,5f12.4)')k,z33(k),tz(i,j,k),sz(i,j,k)
     .     ,uz(i,j,k),vz(i,j,k),trz(i,j,k,ntrcr)
        enddo
      endif

      do k=1,k33
      write(ttl,'(2a,i2,1x,a)') ttl0(4),' k=',k,trim(ttl1)
      write(52) ttl, uz(:,:,k)
      write(ttl,'(2a,i2,1x,a)') ttl0(5),' k=',k,trim(ttl1)
      write(52) ttl, vz(:,:,k)
      write(ttl,'(2a,i2,1x,a)') ttl0(6),' k=',k,trim(ttl1)
      write(52) ttl, tz(:,:,k)
      write(ttl,'(2a,i2,1x,a)') ttl0(7),' k=',k,trim(ttl1)
      write(52) ttl, sz(:,:,k)
      write(ttl,'(2a,i2,1x,a)') ttl0(8),' k=',k,trim(ttl1)
      write(52) ttl, rz(:,:,k)
      do nt=1,ntrcr
      write(ttl,'(2a,i2,1x,a)') ttlt(nt),' k=',k,trim(ttl1)
      write(52) ttl,trz(:,:,k,nt)
      enddo
      enddo
c
      do j=1,jja,lat0
      write(ttl,'(2a,f5.1,1x,a)')ttl0(4),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,uz(:,j,:),xlon,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(5),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,vz(:,j,:),xlon,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(6),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,tz(:,j,:),xlon,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(7),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,sz(:,j,:),xlon,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(8),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,rz(:,j,:),xlon,z33
      do nt=1,ntrcr
      write(ttl,'(2a,f5.1,1x,a)')ttlt(nt),' lat=',xlat(j),trim(ttl2)
      write(53) ttl,iia,k33,1,1,trz(:,j,:,nt),xlon,z33
      enddo
      enddo
c
      do i=1,iia,lon0
      write(ttl,'(2a,f5.1,1x,a)')ttl0(4),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,uz(i,:,:),xlat,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(5),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,vz(i,:,:),xlat,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(6),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,tz(i,:,:),xlat,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(7),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,sz(i,:,:),xlat,z33
      write(ttl,'(2a,f5.1,1x,a)')ttl0(8),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,rz(i,:,:),xlat,z33
      do nt=1,ntrcr
      write(ttl,'(2a,f5.1,1x,a)')ttlt(nt),' lon=',xlon(i),trim(ttl3)
      write(54) ttl,jja,k33,1,1,trz(i,:,:,nt),xlat,z33
      enddo
      enddo
      close(51)
      close(52)
      close(53)
      close(54)
 151  continue   ! mo=mo1,mo2, yr=yr1,yr2
      endif      ! solo_convert
c
      stop '(normal finish get2d)'
      end
