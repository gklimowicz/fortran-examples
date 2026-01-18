      program avg
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- output key monthly mean variables in mon_[runid]_[decade].txt
! --- output key annual  mean variables in ann_[runid]_[decade].txt
! --- output MOC in lat/rho space in 4 basins averaged over year ny1:ny2
!      n avg_ov_[runid]_[decade].txt: flux(idm,kdm,4)
! --- output northward heatflux as a function of lat in "heatfl(idm,4)"
!      n 4 basins averaged over year ny1:ny2 in avg_hf_[runid]_[decade].txt
! --- Last index in flux & heatfl: 1: Atl; 2: Indian; 3: Pac; 4: global
! --- Setting rhodot to true will remove model trend during ny1:ny2 period
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
      use hycom_dimen
      use const_proc
!	 ,   only: path0,path1,path2,hycomtopo,latlonij			&
!         ,basinmask,flnmcoso,flnmo2a,runid,ny1,ny2,spcifh,julian	&
!         ,idrk1,idrk2,jdrk1,jdrk2,indoi,indoj1,indoj2,rhodot,amon	&
!         ,monlg,rho,monave_convert,solo_convert,cnvert,timav,mo1,mo2	&
!         ,iberi,jberi,ikuro1,ikuro2,jkuro,igulf1,igulf2,jgulf		&
!         ,g,imed,jmed
      use hycom_arrays, only: srfhgt,dpmixl,covice,depths,scp2 &
          ,u,v,dp,p,temp,saln,th3d,tracer,uflxav,vflxav,diaflx &
          ,uflx,vflx,alloc_hycom_arrays,latij,lonij            &
          ,temav,salav,ubavav,vbavav,scuy,scvx,xlo,ylo,southfl,eastfl
!
      implicit none
!
      integer :: ntime
      integer :: ny,num,m,k00,ia,ja
      real*8 :: tsum,ssum,sstsum,ssssum,arean3 &
       ,annsst,annsss,annssh,anntem,annsal,annhc,annhc_top,thin,ssh0 &
       ,anndh,anndf,anndb,db,trc &
       ,glbsal,glbtem,glbdep,sum1,sum2,sum3,area,avgbot,vol

      real :: iceextn,iceexts,nino3,day0,flxmax,x1,x2,thrufl	&
             ,tinvrs,fl_beri,kuromax,gulfmax,flow_med,htfx0
      real :: flxmax_i26, flxmax_i45
      real, allocatable :: pinit(:,:,:),pfinl(:,:,:),lat(:),lon(:),	&
            flux(:,:,:),sunda(:),heatfl(:,:),ufx_yr(:,:,:),vfx_yr(:,:,:)
      integer, allocatable :: im(:,:)
!
      real :: dpavav(kdm), heatot
      real*8 :: ztop       ! Depth in m interested top layer ( from surface till 700m )
      real*8 :: hc_top      ! Heat content upper  ztop meters
      real*8 :: zdepth      ! Depth in m from surface to current layer (k)
      real*8 :: zdepth_previous ! Depth in m from surface to previous layer (k-1)
      real*8 :: factor1, factor2
      logical :: is_zdepth_shallow_ztop ! is depth of the current layer less than ztop

      real, allocatable :: year(:), dpav(:,:)
      real, allocatable :: fl_kuro(:), fl_gulf(:)
      character flnm*132,flnmout*80
!
      integer mo,dcd,mon1,i70,i45,i26,ieq,status,i0,i36s
      logical :: lexist
      integer :: nino3_i1,nino3_i2,nino3_j1,nino3_j2   ! index bound for nino3
!
      character(len=9) :: ayears  ! string n1-n2 example "1905-1955"

      character(len=*), parameter :: FMT1=              &
       '(a,         t10,a,     t24,a,     t38,a,      t52, a,		&
        t66,a,      t80,a,     t94,a,     t108,a,     t122,a,		&
        t136,a)'
!
      character(len=*), parameter :: FMT2=				&
       '(f7.2,sp,t10,es12.5,t24,es12.5,t38,es12.5, t52, es12.5,		&
        t66,es12.5, t80,es12.5,t94,es12.5,t108,es12.5,t122,es12.5,	&
        t136,es12.5)'
!
      character(len=*), parameter :: FMT3=				&
       '(a,         t10,a,     t24,a,     t38,a,      t52, a,		&
        t66,a,      t80,a,     t94,a,     t108,a,     t122,a,		&
        t136,a,     t150,a,    t164,a,    t178,a,     t192,a,		&
        t206,a)'
!
      character(len=*), parameter :: FMT4=				&
       '(f7.2,sp,t10,es12.5,t24,es12.5,t38,es12.5, t52, es12.5,		&
        t66,es12.5, t80,es12.5,t94,es12.5,t108,es12.5,t122,es12.5,	&
        t136,es12.5,t150,es12.5,t164,es12.5,t178,es12.5,t192,es12.5,	&
        t206,es12.5)'

      character(len=*), parameter :: FMT5=				&
	'(1x,a,t6,a,t16,a,t28,a,t40,a,t52,a)'
      character(len=*), parameter :: FMT6=				&
       '(1x,i3,sp,t6,f6.1,t16,es10.3,t28,es10.3,t40,es10.3,t52,es10.3)'

      namelist /hdiag_nml/ path0, path1, path2,				&
       hycomtopo, latlonij, basinmask, flnmcoso, flnmo2a,		&
       runid, ny1, ny2, monave_convert,solo_convert

      open (10,file="hdiag.nml")
      read (10,nml=hdiag_nml)
      write (*,nml=hdiag_nml)
      close(10)
!
      write(*,'(3a,i4,a,i4)')						&
         'processing RunId=',trim(runid),' from yr ',ny1,' to ',ny2
      write(*,'(a,i2)') 'number of tracers =',ntrcr

      ztop = 10000.
      ztop = 700.

      ntime=(ny2-ny1+1)*12
      allocate(  year(ntime),dpav(ntime,kdm) )
!
      call alloc_hycom_arrays
      call alloc_hycom_dimen
      allocate (lat(idm),pinit(idm,jdm,kdm+1),pfinl(idm,jdm,kdm+1)	&
          ,flux(idm,kdm,4),im(idm,jdm),heatfl(idm,4),sunda(kdm+1)	&
          ,fl_kuro(idm),fl_gulf(idm),lon(jdm),ufx_yr(idm,jdm,kdm) 	&
	  ,vfx_yr(idm,jdm,kdm)						&
          ,stat=status)
      if (status/=0) stop 'wrong allocate1'

! --- determine do-loop limits for u,v,p,q points
      call const
      call gtdpth(depths,im)
      call bigrid(depths)
!
! --- determine mesh size
      call meshsz
      avgbot=0.
      area=0.
!
      do 10 j=1,jdm
      do 10 l=1,isp(j)
      do 10 i=ifp(j,l),ilp(j,l)
      avgbot=avgbot+depths(i,j)*scp2(i,j)
 10   area=area+scp2(i,j)
      avgbot=avgbot/area
      write (*,104) avgbot,area
 104  format(' mean basin depth (m) and area (10^6 km^2):',f9.1,	&
             -12p,f9.1)
!
! --- read archive data
!
      dcd=ny1/10
      if (runid(1:1).eq.' ') stop 'empty runid'
      if (dcd.lt.001 .or. dcd.gt.930) then
        print *,' wrong decade=',dcd
        stop 'wrong decade'
      endif
      write(ayears,'(i4.4,a1,i4.4)') ny1,'-',ny2
!
      do i=1,idm
      lat(i)=latij(i,340,3)
!     write(*,'(a,i3,a,E12.5)') "lat(",i,")=",lat(i)
      enddo
!
      do j=1,jdm
      lon(j)=lonij(idm,j,3)	! take lon at Southern Ocean
!     write(*,'(a,i3,a,f6.1)') "lon(",j,")=",lon(j)
      enddo
!
! --- nino3 index averaged in 5deg S, 5deg N, 150 W - 90 W
      nino3_i1=-99; nino3_i2=-99; nino3_j1=-99; nino3_j2=-99;
      do i=2,idm-1
        if (lat(i+1).lt.70. .and. lat(i).ge.70.) i70=i
        if (lat(i+1).lt.45. .and. lat(i).ge.45.) i45=i
        if (lat(i+1).lt.26. .and. lat(i).ge.26.) i26=i
        if (lat(i+1).lt. 0. .and. lat(i).ge. 0.) ieq=i
        if (lat(i+1).lt.-36..and. lat(i).ge.-36.) i36s=i

        if (lat(i+1).lt. 5. .and. lat(i).ge. 5.) nino3_i1=i	! i index for 5N
        if (lat(i+1).lt.-5. .and. lat(i).ge.-5.) nino3_i2=i	! i index for 5S
      enddo

      do j=2,jdm-1
        if (lon(j).lt.210. .and. lon(j+1).ge.210.) nino3_j1=j	! j index for 150W
        if (lon(j).lt.270. .and. lon(j+1).ge.270.) nino3_j2=j	! j index for  90W
      end do

      write(*,'(a,4i4)')'chk nino3 index:',nino3_i1,nino3_i2,nino3_j1,nino3_j2
!
      write(flnmout,'(5a)') 'mon_',trim(runid),'_',ayears,'.txt'
      open(301,file=trim(path2)//trim(flnmout),				&
           form='formatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout)

      write(flnmout,'(5a)') 'ann_',trim(runid),'_',ayears,'.txt'
      open(302,file=trim(path2)//trim(flnmout),				&
           form='formatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout)

      write(flnmout,'(5a)') 'avg_ov_',trim(runid),'_',ayears,'.txt'
      open(303,file=trim(path2)//trim(flnmout),				&
           form='formatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout)

      write(flnmout,'(5a)') 'avg_hf_',trim(runid),'_',ayears,'.txt'
      open(304,file=trim(path2)//trim(flnmout),				&
           form='formatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout)

      write(301,fmt=FMT1) 						&
       "Time  ","NINO3","Ice_Extent","Ice_Extent","Ocean_Heat",		&
       "Sea_Surface","Sea_Surface", "Global_Ocean","Global_Ocean",	&
       "Top_Ocn_Heat"
      write(301,fmt=FMT1) 						&
       "---- ","Index","Arctic","Antarctic","Content",			&
       "Temperature","Salinity","Temperature","Salinity","Content"
      write(301,fmt=FMT1)						&
       "Year","None","Mln.Sq.km","Mln.Sq.km","Joule(s)",		&
       "degC","PSU","degC","PSU","Joule(s)"
!
      write(302,fmt=FMT3) 						&
       "Time  ","SST","SSS","Tavrg","Savrg","SSH","Ocean_Heat",		&
       "Atl_(45N)","Indonesian","Drake","Bering","Gulf","Kuroshio",	&
       "Top_Ocn_Heat","Med_outflow","Atl_(26N)"
      write(302,fmt=FMT3) &
       "---- ","Surf_Temp","Surf_Saln","Glob_Temp",  &
       "Glob_Saln","Srf_Hght","Content","Max_Overt", &
       "Throughfl","Passage","Strait","Stream","Current","Content", &
       "pos. E-ward","Max_Overt"
      write(302,fmt=FMT3) &
       "Year","degC","PSU","degC","PSU","cm","Joule(s)","Sv","Sv","Sv" &
       ,"Sv","Sv","Sv","Joule(s)",'Sv','Sv'
!
      write(304,'(a)') "North Poleward Ocean Heat Transport"
      write(304,fmt=FMT5) &
       "Num","Latitude","Atlantic","Indian","Pacifiq","Global"
      write(304,fmt=FMT5) "N","degr(N)","pWatts","pWatts", &
        "pWatts","pWatts"

      n=0
      ufx_yr(:,:,:)=0.   ! multi year mean
      vfx_yr(:,:,:)=0.   ! multi year mean
      heatfl(:,:)=0.     ! multi year mean
      do 151 ny=ny1,ny2

      annsst=0.
      annsss=0.
      anntem=0.
      annsal=0.
      annssh=0.
      annhc =0.        ! Annual heat content
      annhc_top =0.    ! Annual heat content of top ocean layer
      ubavav(:,:)=0.   ! annual u-baro
      vbavav(:,:)=0.   ! annual v-baro
      uflxav(:,:,:)=0.
      vflxav(:,:,:)=0.
!
      do 152 mo=mo1,mo2
      n=n+1
      write(flnm,'(2a,i4.4,3a)') trim(path1),amon(mo),ny &
        ,'.out',trim(runid),'.nc'

      inquire(file=trim(flnm),exist=lexist)
      if( .NOT. lexist ) then
        write(*,'(3a)') "!!! ATTENTION !!! file=",trim(flnm),  &
                " does not exist !!!!"
        write(*,'(3a)') " Calculations continue skip this file "
        cycle
      end if

      write (*,'(2a)') 'do 152 reading: ',trim(flnm)
      call get_time(flnm,year(n))
!
!     call getdat(flnm,day0,lexist)
      call getdat_nc(flnm,day0,lexist)
      if (.not.lexist) stop '(file open or read error)'
!
      write(*,'(a,f9.2)') 'cpl model yr =',year(n)
!
! --- compute total icea volume and ice extent (for each hemisphere)
!
      iceextn=0.
      iceexts=0.
      ssh0   =0.
      sstsum =0.
      ssssum =0.
      tsum   =0.
      hc_top =0.
      ssum   =0.
      trc    =0.
      nino3  =0.
      arean3 =0.
      vol=0.
      factor1=1./(grvty*rho)
      factor1=1.          ! in archive already included
!
      do 5 j=1,jdm
      do 5 l=1,isp(j)
      do 5 i=ifp(j,l),ilp(j,l)

      zdepth=0.0
      zdepth_previous=0  ! depth of previous layer( k-1 )
      is_zdepth_shallow_ztop=.true.

      ia=max(1,i-1)
      do 5 k=1,kdm
      if (k.eq.1) then
        if (i.lt.equat) then
          iceextn=iceextn+covice(i,j)*scp2(i,j)
        else
          iceexts=iceexts+covice(i,j)*scp2(i,j)
        end if
!
        ssh0   =   ssh0+srfhgt(i,j)*scp2(i,j)      ! srfhgt in cm
        sstsum = sstsum+temp(i,j,1)*scp2(i,j)
        ssssum = ssssum+saln(i,j,1)*scp2(i,j)
!
! --- nino3 index averaged in 5deg S, 5deg N, 150 W - 90 W
        if (abs(i-.5*(nino3_i2+nino3_i1)) < .5*(nino3_i2-nino3_i1) .and.	&
            abs(j-.5*(nino3_j2+nino3_j1)) < .5*(nino3_j2-nino3_j1)) then
          nino3=nino3+temp(i,j,1)*scp2(i,j)
          arean3=arean3+scp2(i,j)
        endif
      endif  ! k=1

      uflxav(i,j,k)=uflxav(i,j,k)+uflx(i,j,k)*monlg(mo)	! uflx in Sv
      vflxav(i,j,k)=vflxav(i,j,k)+vflx(i,j,k)*monlg(mo)	! vflx in Sv

      heatfl(i,4)=heatfl(i,4)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
      if (im(i,j).eq.1.or.im(i,j).eq.2) then ! Atlantic
        heatfl(i,1)=heatfl(i,1)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
      elseif (im(i,j).eq.3.or.im(i,j).eq.4) then ! Indian
        heatfl(i,2)=heatfl(i,2)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
      elseif (im(i,j).eq.5.or.im(i,j).eq.6) then ! Pacific
        heatfl(i,3)=heatfl(i,3)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
      endif
!
      zdepth=zdepth + dp(i,j,k)*factor1 ! Current depth (in m) from surface to k layer
      if ( ztop >= zdepth  ) then
        hc_top = hc_top + temp(i,j,k)*dp(i,j,k)*scp2(i,j)
      else
         if( is_zdepth_shallow_ztop ) then
           hc_top = hc_top+temp(i,j,k)*(ztop-zdepth_previous)*scp2(i,j)
           is_zdepth_shallow_ztop=.false.
         end if
       end if
!
      tsum=tsum+temp(i,j,k)*dp(i,j,k)*scp2(i,j)
      ssum=ssum+saln(i,j,k)*dp(i,j,k)*scp2(i,j)
      trc=trc+tracer(i,j,k,1)*dp(i,j,k)*scp2(i,j)
      vol=vol+dp(i,j,k)*scp2(i,j)
      zdepth_previous=zdepth    !  Depth of bottom for previous layer
 5    continue
!
      !heatot=spcifh*rho*tsum/area
      heatot=spcifh*rho*tsum
      !hc_top = spcifh*rho*hc_top/area
      hc_top = spcifh*rho*hc_top
      write(301,fmt=FMT2)					&
       year(n),nino3/arean3,iceextn*1.e-12,iceexts*1.e-12,heatot&
       ,sstsum/area,ssssum/area,tsum/vol,ssum/vol,hc_top
!
! --- save annual field
      mon1=monlg(mod(n-1,12)+1)
      annsst=annsst+sstsum*mon1/julian
      annsss=annsss+ssssum*mon1/julian
      annssh=annssh+  ssh0*mon1/julian
      anntem=anntem  +tsum*mon1/julian
      annsal=annsal  +ssum*mon1/julian
      annhc=annhc  +heatot*mon1/julian
      annhc_top=annhc_top  +hc_top*mon1/julian
!
! --- calculate flux
!     write (*,'(2a)') 'reading: ',trim(flnm)
!     call getdat_nc(flnm,day0,lexist)
!     if (.not.lexist) stop 'file open or read error'
!
 152  continue   !mo=mo1,mo2

      uflxav(:,:,:)=uflxav(:,:,:)/365.		! annual mean in Sv
      vflxav(:,:,:)=vflxav(:,:,:)/365.		! annual mean in Sv
      ufx_yr(:,:,:)=ufx_yr(:,:,:)+uflxav(:,:,:)
      vfx_yr(:,:,:)=vfx_yr(:,:,:)+vflxav(:,:,:)
!
      flux(:,:,:)=0.
      do 181 j=1,jdm
      do 181 i=1,idm
      if(im(i,j).eq.1.or.im(i,j).eq.2) then          ! Atlantic
        do k=1,kdm
          flux(i,k,1)=flux(i,k,1)-uflxav(i,j,k)
        enddo
      endif
 181   continue
!
      do 184 k=2,kdm
      do 184 i=1,idm
 184  flux(i,k,1)=flux(i,k,1)+flux(i,k-1,1)
!
      i=i45                ! get max overturning rate at 45N in Atlantic
      flxmax=-999.
      do 185 k=1,kdm
      if (flux(i,k,1).gt.flxmax) then
        flxmax=flux(i,k,1)
        k00=k
      endif
 185  continue
!      rite(*,*) ' yr n=',n,' flxmax_i45 =',i45,flxmax,' at k=',k00
      flxmax_i45=flxmax
!
      i=i26                ! get max overturning rate at 26N in Atlantic
      flxmax=-999.
      do k=1,kdm
      if (flux(i,k,1).gt.flxmax) then
        flxmax=flux(i,k,1)
        k00=k
      endif
      end do
      flxmax_i26=flxmax

      x1=thrufl(idrk1,jdrk1,idrk2,jdrk2,'(Drake Passage)')
      x2=thrufl(indoi,indoj1,indoi,indoj2,'(Indonesia)')
      fl_beri=thrufl(iberi,jberi-1,iberi,jberi+1,'(Bering)')
      flow_med=thrufl(imed-1,jmed,imed+1,jmed,'(Med)')
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- calculate Gulf Stream and Kuroshio transport by storing vertical integral at k=1
!
!$OMP PARALLEL DO PRIVATE(ja)
      do j=1,jdm
       do i=1,idm
        do k=2,kdm
         uflxav(i,j,1)=uflxav(i,j,1)+uflxav(i,j,k)
         vflxav(i,j,1)=vflxav(i,j,1)+vflxav(i,j,k)
        end do
       end do
      end do
!
!$OMP END PARALLEL DO
      do k=1,2
       call usmoo(uflxav)
       call vsmoo(vflxav)
      end do
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call coagflx(uflxav,vflxav,ylo,southfl,xlo,eastfl)
!
      fl_kuro(:)=0.
      do i=ikuro1,ikuro2
      do j=jkuro,jdm
       if (iv(i,j).gt.0) then
!       print *,'southfl at i,j =',i,j,southfl(i,j)
        fl_kuro(i)=-southfl(i,j)
        if (southfl(i,j).lt.0.) exit
       end if
      end do
      write(*,'(a,i3,a,f8.1)') 'i=',i,' kuroshio transport =',fl_kuro(i)

      end do
      kuromax=maxval(fl_kuro)
      print '(a,f8.1)','max. kuroshio transport:',kuromax

      fl_gulf(:)=0.
      do i=igulf1,igulf2
      do j=jgulf,jdm
       if (iv(i,j).gt.0) then
!       print *,'southfl at i,j =',i,j,southfl(i,j)
        fl_gulf(i)=-southfl(i,j)
        if (southfl(i,j).lt.0.) exit
       end if
      end do
      write(*,'(a,i3,a,f8.1)') 'i=',i,' gulf stream transport =',fl_gulf(i)
      end do
      gulfmax=maxval(fl_gulf)
      print '(a,f8.1)','max. gulfstrm transport:',gulfmax
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      write(302,fmt=FMT4) year(n)-0.5       &
         ,annsst/area,annsss/area,anntem/vol,annsal/vol			&
         ,annssh/area,annhc,flxmax_i45,-x2,x1,-fl_beri,gulfmax,kuromax	&
         ,annhc_top,flow_med,flxmax_i26
 151  continue  !do 151 ny=ny1,ny2
!
      uflxav(:,:,:)=ufx_yr(:,:,:)/(ny2-ny1+1)		! => multi-yr mean in Sv
      vflxav(:,:,:)=vfx_yr(:,:,:)/(ny2-ny1+1)		! => multi-yr mean in Sv
      heatfl(:,:) = heatfl(:,:)/(12.*(ny2-ny1+1))	! => multi-yr mean
!
      write(*,'(2(a,i6))') 'overturning stream function is averaged over years',ny1,' -',ny2
!
      flux(:,:,:)=0.
! --- global domain
      do j=1,jdm
      do k=1,kdm
      do i=1,idm
        flux(i,k,4)=flux(i,k,4)-uflxav(i,j,k)
      enddo
      enddo
      enddo

! --- each basin
      do 81 j=1,jdm
      do 81 i=1,idm
      if (im(i,j).eq.1.or.im(i,j).eq.2) then ! Atlantic
        do k=1,kdm
          flux(i,k,1)=flux(i,k,1)-uflxav(i,j,k)
        enddo
      elseif (im(i,j).eq.3.or.im(i,j).eq.4) then ! Indian
        do k=1,kdm
          flux(i,k,2)=flux(i,k,2)-uflxav(i,j,k)
        enddo
      elseif (im(i,j).eq.5.or.im(i,j).eq.6) then ! Pacific
        do k=1,kdm
          flux(i,k,3)=flux(i,k,3)-uflxav(i,j,k)
        enddo
      endif
 81   continue
!
      do 84 l=1,4
      do 84 i=1,idm
! --- convert to petawatt
      heatfl(i,l)=-.5*heatfl(i,l)*spcifh*rho * 1.e-9          !  N-ward > 0
      do 84 k=2,kdm
 84   flux(i,k,l)=flux(i,k,l)+flux(i,k-1,l)		!vertical integral in k
!
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
      x2=abs(x2)               ! take only absolute value
      write(*,'(a,i4,f6.2,a,i2,a,2(i4,a),2f7.1)') 'chk flxmax_i45 =',	&
       i45,flxmax,' at k=',k00,'; Drake/Indo (ave ',ny1,':',ny2,')=',x1,x2

! --- diagnose indonesian throughflow
!
      sunda(:)=0.
      htfx0=heatfl(indoi,3)-heatfl(indoi-1,3)
      do 26 k=1,kdm
      do 35 j=indoj1,indoj2
 35   sunda(k)=sunda(k)+uflxav(i,j,k)
 26   sunda(k+1)=sunda(k+1)+sunda(k)
!
! --- subtract out portion due to indonesian throughflow
      do 38 i=indoi,i36s
      do 39 k=1,kdm
      flux(i,k,2)=flux(i,k,2)+sunda(k)		!  Indian
 39   flux(i,k,3)=flux(i,k,3)-sunda(k)		!  Pacific
!
      heatfl(i,2)=heatfl(i,2)+htfx0		!  Indian
 38   heatfl(i,3)=heatfl(i,3)-htfx0		!  Pacific
      write(*,'(a,2f6.1)') 'chk massflx/htfx @Indo=',sunda(kdm),htfx0

      do i=i36s+1,idm
      do k=1,3
      heatfl(i,k)=0.
      end do
      end do

      if (rhodot) then
        tinvrs=1./((ny2-ny1+1.)*365.*86400.)
!
! --- determine final pressure field (needed for finding rho-dot)
        pfinl(:,:,:)=p(:,:,:)

        write(flnm,'(5a,i3,a,i3,2a,i4.4,2a)')trim(path1),trim(runid)	&
         ,'/out',trim(runid),'_',dcd,'0_',dcd,'9/','DEC',ny1-1		&
         ,'.out',trim(runid)
        write(*,*) ' input 0 file=',trim(flnm)
! --- determine initial pressure field (needed for finding rho-dot)
        write (*,'(2a)') 'reading initial field: ',trim(flnm)
!       call getdat(flnm,day0,lexist)
        call getdat_nc(flnm,day0,lexist)
        if (.not.lexist) stop '(file open or read error)'
!
! --- subtract out portion of flux attributable to interface displacement
        do 34 j=1,jdm
        do 34 k=2,kdm
        do 34 i=1,idm
        flux(i,k,4)=flux(i,k,4)						&
         +(p(i,j,k)-pfinl(i,j,k))*scp2(i,j)*tinvrs*1.e-6	! => Sv
        if (im(i,j).eq.1.or.im(i,j).eq.2) then		! Atlantic
          flux(i,k,1)=flux(i,k,1)					&
           +(p(i,j,k)-pfinl(i,j,k))*scp2(i,j)*tinvrs*1.e-6	! => Sv
        elseif (im(i,j).eq.3.or.im(i,j).eq.4) then	! Indian
          flux(i,k,2)=flux(i,k,2)					&
           +(p(i,j,k)-pfinl(i,j,k))*scp2(i,j)*tinvrs*1.e-6	! => Sv
        elseif (im(i,j).eq.5.or.im(i,j).eq.6) then	! Pacific
          flux(i,k,3)=flux(i,k,3)					&
           +(p(i,j,k)-pfinl(i,j,k))*scp2(i,j)*tinvrs*1.e-6	! => Sv
        endif
 34     continue
      end if
!
      do l=1,4
      if (idm == 387) then
        write(303,'(387f6.1)') ((flux(i,k,l),i=1,idm),k=1,kdm)
      else if (idm == 359) then
        write(303,'(359f6.1)') ((flux(i,k,l),i=1,idm),k=1,kdm)
      end if
      end do
      close (303)
      do i=1,idm
      write(304,fmt=FMT6) i,lat(i),(heatfl(i,k),k=1,4)
      end do
      close (304)

      stop '(normal finish of avg)'
      end


      subroutine get_time(flnm,year)

! --- extract year and julian day from file name

      character*(*),intent(IN)  :: flnm
      real         ,intent(OUT) :: year

      character*3,parameter :: month(12) =				&
        (/'JAN','FEB','MAR','APR','MAY','JUN',				&
          'JUL','AUG','SEP','OCT','NOV','DEC'/)
      real,parameter :: endmon(0:12) =					&
        (/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334.,365./)

      year=0.
      do n=1,12
       i=index(flnm,month(n))
       if (i.gt.0) then
        read(flnm(i+3:i+6),'(f4.0)') year
        year=year+endmon(n)/365.
        exit
       end if
      end do
      print *,'julian day in input file:',year
      return
      end
