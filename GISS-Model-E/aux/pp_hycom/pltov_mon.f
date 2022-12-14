      program pltov_mon
c --- calculate MONTHLY
c  (1) overturning stream function lat vs. rho in "flux(idm,kdm,4)"
c      in mon_ov_[runid]_[decade].tbin: flux(idm,kdm,4)
c  (2) northward heatflux as a function of lat in "heatfl(idm,4)"
c      in avg_hf_[runid]_[decade].txt
c      Last index in flux & heatfl: 1: Atl; 2: Indian; 3: Pac; 4: global
c --- setting rhodot to true will remove model trend during ny1:ny2 period
c
      use hycom_arrays, only: depths,scp2
     .    ,u,v,dp,p,temp,saln,th3d,diaflx
     .    ,uflx,vflx,alloc_hycom_arrays,latij
     .    ,ubavg,vbavg,ubavav,vbavav
      use hycom_dimen
      use const_proc
      use TimeConstants_mod, only: SECONDS_PER_DAY
      implicit none
c
      real :: day0,day1,flxmax,x1,x2,thrufl,tinvrs
      real, allocatable :: pinit(:,:,:),pfinl(:,:,:),lat(:),flux(:,:,:)
     .       ,sunda(:)
      integer, allocatable :: im(:,:)
      real, allocatable :: heatfl(:,:)
      real*4, allocatable :: heatfl_r4(:,:),flux_r4(:,:,:)
c
      integer mo,dcd,mon1,i70,i45,ieq,status,ia,k00,ny
      real*8 :: area,avgbot
      character flnm*128,flnmout*30,title*80
      logical succes
c
      character(len=9) :: ayears  ! string n1-n2 example "1905-1955"
c
      namelist /hdiag_nml/ path0, path1, path2,
     . hycomtopo, latlonij, basinmask, flnmcoso, flnmo2a,
     . runid, ny1, ny2, monave_convert,solo_convert

      open (10,file="hdiag.nml")
      read (10,nml=hdiag_nml)
      write (*,nml=hdiag_nml)
      close(10)
c
      write(*,'(3a,i4,a,i4)')
     .   'processing ',runid,' from yr ',ny1,' to ',ny2
      write(*,'(a,i2)') 'number of tracers =',ntrcr
c
      call alloc_hycom_arrays
      call alloc_hycom_dimen
      allocate (lat(idm),pinit(idm,jdm,kdm+1),pfinl(idm,jdm,kdm+1)
     .    ,flux(idm,kdm,4),im(idm,jdm),heatfl(idm,4),sunda(kdm+1)
     .    ,flux_r4(idm,kdm,4),heatfl_r4(idm,4)
     .    ,stat=status)
      if (status/=0) stop 'wrong allocate1'

c --- determine do-loop limits for u,v,p,q points
      call gtdpth(depths,im)
      call bigrid(depths)
c
c --- determine mesh size
      call meshsz
      avgbot=0.
      area=0.
c
      do 10 j=1,jdm
      do 10 l=1,isp(j)
      do 10 i=ifp(j,l),ilp(j,l)
      avgbot=avgbot+depths(i,j)*scp2(i,j)
 10   area=area+scp2(i,j)
      avgbot=avgbot/area
      write (lp,104) avgbot,area
 104  format(' mean basin depth (m) and area (10^6 km^2):',f9.3,
     .       -12p,f9.3)
c
      dcd=ny1/10
      if (runid(1:1).eq.' ') stop 'empty runid'
      if (dcd.lt.180 .or. dcd.gt.999) then
        print *,' wrong decade=',dcd
        stop 'wrong decade'
      endif
c
      do i=1,idm
      lat(i)=latij(i,1,3)
      enddo
c
      do i=2,idm-1
      if (lat(i+1).lt.70. .and. lat(i).ge.70.) i70=i
      if (lat(i+1).lt.45. .and. lat(i).ge.45.) i45=i
      if (lat(i+1).lt. 0. .and. lat(i).ge. 0.) ieq=i
      enddo
c
      write(ayears,'(i4.4,a1,i4.4)') ny1,'-',ny2
      write(flnmout,'(5a)') 'mon_ov_',trim(runid),'_',ayears,'.tbin'
      open(301,file=trim(path2)//flnmout, 
     +     form='unformatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout) 

C     write(flnmout,'(3a,i3,a)') 'mon_ov_',trim(runid),'_',dcd,'.tbin'
C     open(301,file=flnmout,form='unformatted',status='unknown')
      write(flnmout,'(5a)') 'mon_hf_',trim(runid),'_',ayears,'.tbin'
      open(302,file=trim(path2)//flnmout,
     +     form='unformatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout) 

C     write(flnmout,'(3a,i3,a)') 'mon_hf_',trim(runid),'_',dcd,'.tbin'
C     open(302,file=flnmout,form='unformatted',status='unknown')
c
      do 151 ny=ny1,ny2
      do 152 mo=mo1,mo2
c     write(flnm,'(5a,i3,a,i3,2a,i4.4,2a)')trim(path1),trim(runid)
c    . ,'/out',trim(runid),'_',dcd,'0_',dcd,'9/',amon(mo),ny
c    . ,'.out',trim(runid)

      write(flnm,'(2a,i4.4,2a)')trim(path1),amon(mo),ny
     .  ,'.out',trim(runid)

      write (lp,'(2a)') 'reading: ',trim(flnm)
      mon1=monlg(mod(mo-1,12)+1)
c
c --- read archive data
      timav=.true.
      cnvert=.true.
      call getdat(flnm,day0,day1,succes)

      do 13 i=1,idm
      ia=max(1,i-1)
      do 13 j=1,jdm
      ubavav(i,j)=ubavg(i,j,1)
      vbavav(i,j)=vbavg(i,j,1)
 13   continue
c
      uflx=uflx/(mon1*SECONDS_PER_DAY)  ! Sv
      vflx=vflx/(mon1*SECONDS_PER_DAY)  ! Sv
c
      flux(:,:,:)=0.
      heatfl(:,:)=0.
c --- global domain
      do i=1,idm
      ia=max(1,i-1)
      do j=1,jdm
      do k=1,kdm
        flux(i,k,4)=flux(i,k,4)-uflx(i,j,k)
        heatfl(i,4)=heatfl(i,4)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
      enddo
      enddo
      enddo

c --- each basin
      do 81 i=1,idm
      ia=max(1,i-1)
      do 81 j=1,jdm
      if (im(i,j).eq.1.or.im(i,j).eq.2) then ! Atlantic
        do k=1,kdm
          flux(i,k,1)=flux(i,k,1)-uflx(i,j,k)
          heatfl(i,1)=heatfl(i,1)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
        enddo
      elseif (im(i,j).eq.3.or.im(i,j).eq.4) then ! Indian
        do k=1,kdm
          flux(i,k,2)=flux(i,k,2)-uflx(i,j,k)
          heatfl(i,2)=heatfl(i,2)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
        enddo
      elseif (im(i,j).eq.5.or.im(i,j).eq.6) then ! Pacific
        do k=1,kdm
          flux(i,k,3)=flux(i,k,3)-uflx(i,j,k)
          heatfl(i,3)=heatfl(i,3)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
        enddo
      endif
 81   continue
c
      do 84 l=1,4
      do 84 i=1,idm
c --- convert to petawatt
      heatfl(i,l)=-.5*heatfl(i,l)*spcifh*rho * 1.e-9                ! N-ward > 0
      do 84 k=2,kdm
 84   flux(i,k,l)=flux(i,k,l)+flux(i,k-1,l)
c
      i=i45                ! get max overturning rate at 45N in Atlantic
      flxmax=-999.
      do 85 k=1,kdm
      if (flux(i,k,1).gt.flxmax) then
        flxmax=flux(i,k,1)
        k00=k
      endif
 85   continue
      x1=thrufl(idrk1,jdrk1,idrk2,jdrk2,'(Drake Passage)')
      x2=thrufl(indoi,indoj1,indoi,indoj2,'(Indonesia)')
      x2=-x2               ! take only absolute value
c     write(lp,'(a,i4,f6.2,a,i2,a,2f6.1)') 
c    . 'chk flxmax_i45 =',i45,flxmax,' at k=',k00,'; Drake/Indo=',x1,x2

c --- diagnose indonesian throughflow
c
      sunda=0.
      do 26 k=1,kdm
      i = indoi
      do 35 j=indoj1,indoj2
 35   sunda(k)=sunda(k)+uflx(i,j,k)
 26   sunda(k+1)=sunda(k+1)+sunda(k)
c
c     write(*,'(a,27f4.0)') 'sunda=',(sunda(k),k=1,kdm+1)
c     write(*,'(a,26f4.0)') 'Indo_pa=',(flux(indoi,k,3),k=1,kdm)
c     write(*,'(a,26f4.0)') 'Indo_in=',(flux(indoi,k,2),k=1,kdm)

c --- subtract out portion due to indonesian throughflow
      do 39 k=1,kdm
      do 29 i=indoi+1,idm
 29   flux(i,k,2)=flux(i,k,2)+sunda(k)                !  Indian
      do 39 i=indoi+1,idm
 39   flux(i,k,3)=flux(i,k,3)-sunda(k)                !  Pacific
c     write(*,'(a,26f4.0)') 'pa flux=',(flux(indoi,k,3),k=1,kdm)
c     write(*,'(a,26f4.0)') 'in flux=',(flux(indoi,k,2),k=1,kdm)
c
c     do l=1,4
c     write(301,'(387f6.1)') ((flux(i,k,l),i=1,idm),k=1,kdm)
c     end do
c     do i=1,idm
c     write(302,'(i3,f6.1,8f8.2)') i,lat(i),(heatfl(i,k),k=1,4)
c     end do

      flux_r4=flux
      write(title,'(a,a,i4)') 'overturning rate flux(387,26,4) in Sv '
     .   ,amon(mo),ny
      write(301) title,flux_r4
      heatfl_r4=heatfl
      write(title,'(a,a,i4)') 'northward heatflux(387,4) in PW '
     .   ,amon(mo),ny
      write(302) title,heatfl_r4
 152  continue
 151  continue   ! ny=ny1,ny2
      close (301)
      close (302)

      stop '(normal finish of pltov)'
      end
