!#include "hycom_mpi_hacks.h"
#include "rundeck_opts.h"
      subroutine archiv(n,nn)
c
c --- write archive file for time level n to flnm ( b i n a r y  hycom fmt)
c
      use netcdf
      use cdf_io
      USE JulianCalendar_mod, only : jdendofm
      USE MODEL_COM, only : modelEclock,
     *  itime,iyear1,nday,aMON,xlabel,lrunid,monthi,datei
      USE TimeConstants_mod, only: SECONDS_PER_DAY
      USE HYCOM_SCALARS, only : nstep,time,theta,huge,baclin,onem
     &     ,thref,nhr,g,delt1,pr1d
      USE HYCOM_DIM_GLOB, only : ii1,jj,JDM,kk,isp,ifp,ilp,ntrcr,isu
     &     ,ifu,ilu,isv,ifv,ilv,ii,idm,kdm
      USE HYCOM_ARRAYS_GLOB
      use hycom_arrays_glob_renamer
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: UNPACK_DATA

#ifdef TRACERS_OceanBiology
      USE obio_com, only : pCO2av,ao_co2fluxav,diag_counter
     .                    ,cexpav,pp2tot_dayav
#ifdef TRACERS_Alkalinity
     .                    ,caexpav
#endif
#endif

      use filemanager, only : findunit
c
      implicit none
      integer i,j,k,l,n,nn,kn
c
      integer no,nop,length,nt
      real factor,vol,tts2,temavg,sst,dpsmo(idm,jdm,kdm)
     .    ,icearea,icevol,icearean,icevoln,iceareas,icevols
      character flnm*40,intvl*4,flnm_nc*40
      character what*16
      real*4 real4(idm,jdm),time4,thref4,theta4(kdm),unused,pr1d4(kdm)
      real utotal(idm,jdm,kdm), vtotal(idm,jdm,kdm)
     .       ,dpm(idm,jdm,kdm),dpmixlm(idm,jdm)
      integer*4 length4,idm4,jdm4,kdm4,nstep4
      integer*4 irecl ! specific record lenth, machine dependent
      logical, parameter :: smooth = .false.     ! smooth fields before saving
      data unused/0./
      logical, parameter :: binary_output=.false., ncout=.true.
      integer ncid1
      integer, parameter ::
     . mon_date(13)=(/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer :: year, month, dayOfYear, date, hour
c
      call modelEclock%get(year=year, month=month, date=date,
     .  hour=hour, dayOfYear=dayOfYear)
      call getdte(Itime,Nday,Iyear1,year,month,dayOfYear,date,hour,amon)
      print *,'stamp',Itime,Nday,Iyear1,year,month,dayOfYear,date,hour,
     .  amon

      dpm(:,:,:)=huge
      dpmixlm(:,:)=huge
      do k=1,kk
      do j=1,jj
      do l=1,isp(j)
      do i=ifp(j,l),ilp(j,l)
      dpm(i,j,k)=dp(i,j,k)/onem                 ! convert to m
      if (k==1) dpmixlm(i,j)=dpmixl(i,j,n)/onem ! convert to m
      end do
      end do
      end do
      end do
c --- check if ogcm date matches agcm date
      if (nstep.eq.1) then
        write(flnm,'(a3,i4.4,2a)') amon,0,'.out',xlabel(1:lrunid)
        write(flnm_nc,'(a3,i4.4,3a)')
     .     amon,0,'.out',xlabel(1:lrunid),'.nc'
        temav(:,:,:)=temp(:,:,:)
        salav(:,:,:)=saln(:,:,:)
        th3av(:,:,:)=th3d(:,:,:)
        dpav (:,:,:)= dpm(:,:,:)
        oiceav(:,:) =oice(:,:)
      elseif (abs((itime+1.)/nday-time).gt.1.e-5) then
c --- check if ogcm date matches agcm date
        write(*,*) 'mismatching archive date in agcm/ogcm=',
     .     (itime+1.)/nday,time
c       stop 'mismatching archive date'
      else
        write(flnm,'(a3,i4.4,2a)') amon,year,'.out',xlabel(1:lrunid)
        write(flnm_nc,'(a3,i4.4,3a)')
     .    amon,year,'.out',xlabel(1:lrunid),'.nc'
      endif
c
      if (date.le.999) then ! date = number of days in this month
        write (intvl,'(i4.4)') date
      else
        stop ' wrong date > 999'
      endif
c
      utotal(:,:,:)=huge
      vtotal(:,:,:)=huge
c --- output total velocity
      do k=1,kk
       do j=1,jj
        do l=1,isu(j)
         do i=ifu(j,l),ilu(j,l)
          utotal(i,j,k)=u(i,j,k+nn)+ubavg(i,j,n)
         end do
        end do

        do l=1,isv(j)
         do i=ifv(j,l),ilv(j,l)
          vtotal(i,j,k)=v(i,j,k+nn)+vbavg(i,j,n)
         end do
        end do
       end do      ! do j
      end do       !do k

      if (smooth) then
c
        do j=1,jj
         do k=1,kk
          do l=1,isp(j)
           do i=ifp(j,l),ilp(j,l)
            p(i,j,k+1)=p(i,j,k)+dpm(i,j,k+nn)
           end do
          end do
        end do
       end do
c
       do k=2,kk
        call psmo1(p(1,1,k),pbot)
       end do
c
      dpsmo=huge
       do j=1,jj
        do k=1,kk
         do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            dpsmo(i,j,k)=p(i,j,k+1)-p(i,j,k)
          end do
         end do
        end do
       end do
c
      end if                            !  smooth
c
      if (binary_output) then
      no=4096
      inquire (iolength=irecl)  real4(1,1) ! length of an unformatted real*4
                                           ! irecl=1 on COMPAQ, irecl=4 on SGI
      length=((irecl*idm*JDM+no+15)/no)*no
c
      write (*,'(a/9x,a)') 'storing history data in',flnm
c
      call findunit(nop)
      open (unit=nop,file=flnm,status='unknown',form='unformatted',
     .      access='direct',recl=length)
      no=1
      length4=length
      idm4=idm
      jdm4=jdm
      kdm4=kdm
      nstep4=nstep
      time4=time
      do k=1,kk
        theta4(k)=theta(k)
        pr1d4(k) =pr1d(k)/onem     ! TNL convert in meter
      end do
      write(flnm(1:17),'(1x,2(i2.2,a),i4.4,a,i2)')
     .             month,'/',date,'/',year,' hr ',hour+nhr
      print *,' flnm(1:17)=',flnm(1:17)
      write (nop,rec=no) length4,idm4,jdm4,kdm4,nstep4,time4
!TNL .      ,unused,theta4,      flnm(1:17)
     .      ,unused,theta4,pr1d4,flnm(1:17)
c
      no=no+1
      call r8tor4(srfhgt,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'srfhgt (m)      ',0,real4
      write (*,100)     'srfhgt (m)      ',0,no
      no=no+1
      call r8tor4(dpmixlm,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'mix_dpth(m)     ',0,real4
      write (*,100)     'mix_dpth(m)     ',0,no
      no=no+1
      call r8tor4(oice,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'icecover(%)     ',0,real4
      write (*,100)     'icecover(%)     ',0,no
c
      do 75 k=1,kk
      kn=k+nn
      no=no+1
      call r8tor4(utotal(1,1,kn),real4)
      if (smooth) call usmoo4(real4)
      write (nop,rec=no) 'u               ',k,real4
      write (*,100)     'u               ',k,no
      no=no+1
      call r8tor4(vtotal(1,1,kn),real4)
      if (smooth) call vsmoo4(real4)
      write (nop,rec=no) 'v               ',k,real4
      write (*,100)     'v               ',k,no
      no=no+1
      if (smooth) then
        call r8tor4(dpsmo(1,1,k),real4)
      else
        call r8tor4(dpm(1,1,kn),real4)
      endif
      write (nop,rec=no) 'dp(m)           ',k,real4 ! unit in m
      write (*,100)     'dp(m)           ',k,no
      no=no+1
      call r8tor4(temp(1,1,kn),real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'temp            ',k,real4
      write (*,100)     'temp            ',k,no
      no=no+1
      call r8tor4(saln(1,1,kn),real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'saln            ',k,real4
      write (*,100)     'saln            ',k,no
      no=no+1
      call r8tor4(th3d(1,1,kn),real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'th3d            ',k,real4
      write (*,100)     'th3d            ',k,no
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- ifort miscomputes 'no' in the following loop:
ccc      do nt=1,ntrcr
ccc        no=no+1
ccc        call r8tor4(tracer(1,1,k,nt),real4)
ccc        if (smooth) call psmoo4(real4)
ccc        write (what,'(a6,i2,4x)') 'tracer',nt
ccc        write (nop,rec=no) what,k,real4
ccc      write (*,100)       what,k,no
ccc      end do
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- temporary: store diffusivity as tracer 1
c     tracer(:,:,:,1)=diadff*1.e4
c     write (what,'(a12)') 'diapyc.diffu'
c
c --- code around compiler glitch:
      do nt=1,ntrcr
        call r8tor4(tracer(1,1,k,nt),real4)
        if (smooth) call psmoo4(real4)
        write (what,'(a6,i2,5x)') 'tracer',nt
        write (nop,rec=no+nt) what,k,real4
      write (*,100)       what,k,no+nt
      end do
      no=no+ntrcr
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 75   continue
      end if    ! if binary_output
c
c --- output time-averaged fields
c
      factor=baclin/(date*SECONDS_PER_DAY)
c
      if (nstep > 1) then
      do 55 j=1,jj
      do 55 l=1,isp(j)
      do 55 i=ifp(j,l),ilp(j,l)
      eminpav(i,j)=eminpav(i,j)*factor*1000.*SECONDS_PER_DAY !mm/day
      surflav(i,j)=surflav(i,j)*factor
      salflav(i,j)=salflav(i,j)*factor
      brineav(i,j)=brineav(i,j)*factor
      tauxav(i,j)=tauxav(i,j)*factor
      tauyav(i,j)=tauyav(i,j)*factor
 55   continue
c
#ifdef TRACERS_OceanBiology
      if (diag_counter .ne. 0.d0) then
      do j=1,jj; do l=1,isp(j); do i=ifp(j,l),ilp(j,l)
        ao_co2fluxav(i,j)=ao_co2fluxav(i,j)/diag_counter
        pco2av(i,j)=pco2av(i,j)/diag_counter
        pp2tot_dayav(i,j)=pp2tot_dayav(i,j)/diag_counter
        cexpav(i,j)=cexpav(i,j)/diag_counter
#ifdef TRACERS_Alkalinity
        caexpav(i,j)=caexpav(i,j)/diag_counter
#endif
       enddo; enddo; enddo
        i=243;j=1;
        write(*,'(a,3i5,f10.3,4e12.4)')'2222222222',
     .   nstep,i,j,diag_counter,ao_co2fluxav(i,j),
     .   pco2av(i,j),pp2tot_dayav(i,j),cexpav(i,j)
      endif
#endif

      do 58 k=1,kk
c
      do 59 j=1,jj
      do 591 l=1,isu(j)
      do 591 i=ifu(j,l),ilu(j,l)
      if (dpuav(i,j,k).gt.0) uav(i,j,k)=uav(i,j,k)/dpuav(i,j,k)
      ufxavp(i,j,k)=ufxavp(i,j,k)*baclin*1.e-6/(date*86400.*onem)	! in Sv
 591  uflxav(i,j,k)=uflxav(i,j,k)*baclin*1.e-6/(date*86400.*onem)	! in Sv
      do 592 l=1,isv(j)
      do 592 i=ifv(j,l),ilv(j,l)
      if (dpvav(i,j,k).gt.0) vav(i,j,k)=vav(i,j,k)/dpvav(i,j,k)
      vfxavp(i,j,k)=vfxavp(i,j,k)*baclin*1.e-6/(date*86400.*onem)	! in Sv
 592  vflxav(i,j,k)=vflxav(i,j,k)*baclin*1.e-6/(date*86400.*onem)	! in Sv
      do 59 l=1,isp(j)
      do 59 i=ifp(j,l),ilp(j,l)
      if (dpav(i,j,k).gt.0.) then
        temav(i,j,k)=temav(i,j,k)/dpav(i,j,k)
        salav(i,j,k)=salav(i,j,k)/dpav(i,j,k)
        th3av(i,j,k)=th3av(i,j,k)/dpav(i,j,k)
      end if
      dpav(i,j,k)=dpav(i,j,k)*factor/onem   ! in meter
c
      diaflx(i,j,k)=-diaflx(i,j,k)/(2.*onem)
c --- convert diapycnal thickness changes into actual interface fluxes
      if (k.gt.1) diaflx(i,j,k)=diaflx(i,j,k)+diaflx(i,j,k-1)
 59   continue
ccc      write (*,'(a,i3)') 'shown below: N.Atl. diaflx, bottm of layer',k
ccc      call zebra(diaflx(1,int(.8*jdm),k),idm,idm/3,idm/3)
c
 58   continue
c
      do 56 j=1,jj
      do 56 l=1,isp(j)
      do 56 i=ifp(j,l),ilp(j,l)
      pbavav(i,j)=pbavav(i,j)*factor
      sfhtav(i,j)=sfhtav(i,j)*factor      ! in meter
      dpmxav(i,j)=dpmxav(i,j)*factor/onem ! in meter
      oiceav(i,j)=oiceav(i,j)*factor
 56   continue
c
      end if      ! nstep > 1
c
c     write (*,'(3a,i5)') 'shown below: ',intvl
c    .    ,'- day SSH average step=',nstep
c     call zebra(sfhtav,idm,ii1,jj)
c     write (*,'(3a,i5)') 'shown below: ',intvl
c    .    ,'- day SST average, step=',nstep
c     call zebra(temav,idm,ii1,jj)
c     write (*,'(3a,i5)') 'shown below: ',intvl
c    .    ,'- day ice average, step=',nstep
c     call zebra(oiceav,idm,ii1,jj)
c
      if (binary_output) then
      no=no+1
      call r8tor4(uflxav(1,1,k),real4)
      write (nop,rec=no) '     uflxav_'//intvl,k,real4
      write (*,100)     '     uflxav_'//intvl,k,no
      no=no+1
      call r8tor4(vflxav(1,1,k),real4)
      write (nop,rec=no) '     vflxav_'//intvl,k,real4
      write (*,100)     '     vflxav_'//intvl,k,no
      no=no+1
      call r8tor4(ufxavp(1,1,k),real4)
      write (nop,rec=no) '     ufxavp_'//intvl,k,real4
      write (*,100)     '     ufxavp_'//intvl,k,no
      no=no+1
      call r8tor4(vfxavp(1,1,k),real4)
      write (nop,rec=no) '     vfxavp_'//intvl,k,real4
      write (*,100)     '     vfxavp_'//intvl,k,no
      no=no+1
      call r8tor4(diaflx(1,1,k),real4)
      write (nop,rec=no) '     diaflx_'//intvl,k,real4
      write (*,100)     '     diaflx_'//intvl,k,no
      no=no+1
      call r8tor4(dpmxav,real4)
      write (nop,rec=no) '     dpmxav_'//intvl,0,real4
      write (*,100)     '     dpmxav_'//intvl,0,no
      no=no+1
      call r8tor4(oiceav,real4)
      write (nop,rec=no) '     oiceav_'//intvl,1,real4
      write (*,100)     '     oiceav_'//intvl,0,no
c
#ifdef TRACERS_OceanBiology
      no=no+1
      call r8tor4(ao_co2fluxav,real4)
      write (nop,rec=no) '  ao_co2flux'//intvl,1,real4
      write (*,100)     '  ao_co2flux'//intvl,0,no

      no=no+1
      call r8tor4(pco2av,real4)
      write (nop,rec=no) '      pco2av'//intvl,1,real4
      write (*,100)     '      pco2av'//intvl,0,no

      no=no+1
      call r8tor4(pp2tot_dayav,real4)
      write (nop,rec=no) 'pp2tot_dayav'//intvl,1,real4
      write (*,100)     'pp2tot_dayav'//intvl,0,no

      no=no+1
      write(*,*)'archyb1: ',cexpav(243,1)
      call r8tor4(cexpav,real4)
      write(*,*)'archyb2: ',real4(243,1)
      write (nop,rec=no) '      cexpav'//intvl,1,real4
      write (*,100)     '      cexpav'//intvl,0,no

#ifdef TRACERS_Alkalinity
      no=no+1
      call r8tor4(caexpav,real4)
      write (nop,rec=no) '     caexpav'//intvl,1,real4
      write (*,100)     '     caexpav'//intvl,0,no
#endif
#endif
c
      do 57 k=1,kk
      no=no+1
      call r8tor4(uav(1,1,k),real4)
      write (nop,rec=no) '        uav_'//intvl,k,real4
      write (*,100)     '        uav_'//intvl,k,no
      no=no+1
      call r8tor4(vav(1,1,k),real4)
      write (nop,rec=no) '        vav_'//intvl,k,real4
      write (*,100)     '        vav_'//intvl,k,no
      no=no+1
      call r8tor4(dpav(1,1,k),real4)
      write (nop,rec=no) '       dpav_'//intvl,k,real4
      write (*,100)     '       dpav_'//intvl,k,no
      no=no+1
      call r8tor4(temav(1,1,k),real4)
      write (nop,rec=no) '      temav_'//intvl,k,real4
      write (*,100)     '      temav_'//intvl,k,no
      no=no+1
      call r8tor4(salav(1,1,k),real4)
      write (nop,rec=no) '      salav_'//intvl,k,real4
      write (*,100)     '      salav_'//intvl,k,no
      no=no+1
      call r8tor4(th3av(1,1,k),real4)
      write (nop,rec=no) '     th3dav_'//intvl,k,real4
      write (*,100)     '     th3dav_'//intvl,k,no
 57   continue
c
c --- time-averaged surface fluxes:
      no=no+1
      call r8tor4(eminpav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) '  eminp m/s_'//intvl,0,real4
      write (*,100)     '  eminp m/s_'//intvl,0,no
      no=no+1
      call r8tor4(surflav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) ' htflx W/m2_'//intvl,0,real4
      write (*,100)     ' htflx W/m2_'//intvl,0,no
      no=no+1
      call r8tor4(salflav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) ' sflx g/m2s_'//intvl,0,real4
      write (*,100)     ' sflx g/m2s_'//intvl,0,no
      no=no+1
      call r8tor4(brineav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'brine g/m2s_'//intvl,0,real4
      write (*,100)     'brine g/m2s_'//intvl,0,no
      no=no+1
      call r8tor4(tauxav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) ' Tau_x N/m2_'//intvl,0,real4
      write (*,100)     ' Tau_x N/m2_'//intvl,0,no
      no=no+1
      call r8tor4(tauyav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) ' Tau_y N/m2_'//intvl,0,real4
      write (*,100)     ' Tau_y N/m2_'//intvl,0,no
c
      close (unit=nop)
 100  format (9x,a,' (layer',i3,') archived as record',i5)
      write (*,*) no,' records archived'
c
      end if   ! if binary_output

!<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
      if (ncout) then                   ! archive data in netcdf format
        print *,' flnm=',flnm_nc
        print '(3a)','opening ',trim(flnm_nc),' for netcdf output'
        call errhandl (nf90_create (path=trim(flnm_nc),
     .      cmode=or(NF90_CLOBBER, NF90_64BIT_OFFSET),ncid=ncid1))
        call out2cdf(ncid1,idm,jdm,depths,0.,
     .    'topo','bottom depth','m')
        call out2cdf(ncid1,idm,jdm,scp2,0.,
     .    'scp2','grid cell size','m^2')
        call out2cdf(ncid1,idm,jdm,latij(1,1,3),0.,
     .    'lat','latitude','degrees')
        call out2cdf(ncid1,idm,jdm,lonij(1,1,3),0.,
     .    'lon','longitude','degrees')
        call out1cdf(ncid1,kdm,theta,0.,
     .    'theta','target pot.density, sigma1','kg/m^3')
        call out1cdf(ncid1,kdm,pr1d/onem,0.,
     .    'msf_lvls','z levels for overturning msf','m')
! snapshot
        call out2cdf(ncid1,idm,jdm,srfhgt,time,
     .    'srfht','sea surface height','m')
        call out2cdf(ncid1,idm,jdm,dpmixlm,time,	! unit in m
     .    'zmixl','mixed layer depth','m')
        call out2cdf(ncid1,idm,jdm,oice,time,
     .    'covice','ice coverage','fraction')
        call out2cdf(ncid1,idm,jdm,omlhc,time,
     .    'omlhc','ocean mixed layer heat content','unk')
        call out2cdf(ncid1,idm,jdm,osalt,time,
     .    'osalt','salt flux due to brine rejection','kg/s/m^2')
        call out3cdf(ncid1,idm,jdm,kdm,temp,time,
     .    'temp','potential temperature','deg C')
        call out3cdf(ncid1,idm,jdm,kdm,saln,time,
     .    'saln','salinity','psu')
        call out3cdf(ncid1,idm,jdm,kdm,th3d,time,
     .    'th3d','pot.density, sigma1','kg/m^3')
        call out3cdf(ncid1,idm,jdm,kdm,dpm,time,	! unit in m
     .    'thik','layer thickness','m')
        call out3cdf(ncid1,idm,jdm,kdm,utotal,time,
     .    'utotal','southward velocity','m/sec')
        call out3cdf(ncid1,idm,jdm,kdm,vtotal,time,
     .    'vtotal','eastward velocity','m/sec')
! monthly average
        call out2cdf(ncid1,idm,jdm,sfhtav,time,
     .    'srfhtav','monthly sea surface height','m')
        call out2cdf(ncid1,idm,jdm,dpmxav,time,
     .    'zmixlav','monthly mixed layer depth','m')
!       call out2cdf(ncid1,idm,jdm,ticeav,time,
!    .    'temiceav','monthly ice surface temp','deg C')
        call out2cdf(ncid1,idm,jdm,oiceav,time,
     .    'coviceav','monthly ice coverage','fraction')
        call out2cdf(ncid1,idm,jdm,eminpav,time,
     .    'eminpav','monthly eminp','mm/day')
        call out2cdf(ncid1,idm,jdm,surflav,time,
     .    'surflav','monthly net sfc htflx','W/m2')
        call out2cdf(ncid1,idm,jdm,salflav,time,
     .    'salflav','monthly salflx','g/m2s')
        call out2cdf(ncid1,idm,jdm,brineav,time,
     .    'brineav','monthly brine','g/m2s')
        call out2cdf(ncid1,idm,jdm,tauxav,time,
     .    'tauxav','monthly taux','N/m2')
        call out2cdf(ncid1,idm,jdm,tauyav,time,
     .    'tauyav','monthly tauy','N/m2')
        call out2cdf(ncid1,idm,jdm,diag1,time,
     .    'odmsi','monthly odmsi','kg/m2')
        call out2cdf(ncid1,idm,jdm,diag2,time,
     .    'runsi','monthly runsi','kg/m2')

        call out3cdf(ncid1,idm,jdm,kdm,temav,time,
     .    'tempav','monthly potential temperature','deg C')
        call out3cdf(ncid1,idm,jdm,kdm,salav,time,
     .    'salnav','monthly salinity','psu')
        call out3cdf(ncid1,idm,jdm,kdm,th3av,time,
     .    'th3dav','monthly pot.density, sigma1','kg/m^3')
        call out3cdf(ncid1,idm,jdm,kdm,dpav,time,
     .    'thikav','monthly averaged layer thickness','m')
        call out3cdf(ncid1,idm,jdm,kdm,uav,time,
     .    'uav','monthly southward velocity','m/sec')
        call out3cdf(ncid1,idm,jdm,kdm,vav,time,
     .    'vav','monthly eastward velocity','m/sec')
        if (ntrcr.ge.1)
     .  call out3cdf(ncid1,idm,jdm,kdm,tracer(1,1,1,1),time,
     .    'trc1','passive tracer 1',' ')
        if (ntrcr.ge.2)
     .  call out3cdf(ncid1,idm,jdm,kdm,tracer(1,1,1,2),time,
     .    'trc2','passive tracer 2',' ')
        if (ntrcr.ge.3)
     .  call out3cdf(ncid1,idm,jdm,kdm,tracer(1,1,1,3),time,
     .    'trc3','passive tracer 3',' ')
        if (ntrcr.ge.4)
     .  call out3cdf(ncid1,idm,jdm,kdm,tracer(1,1,1,4),time,
     .    'trc4','passive tracer 4',' ')
        if (ntrcr.ge.5) stop 'stop: need work for ntrcr > 4'
        call out3cdf(ncid1,idm,jdm,kdm,uflxav,time,'uflxav',
     .    'monthly integral of southward isopycnic mass flux','Sv')
        call out3cdf(ncid1,idm,jdm,kdm,vflxav,time,'vflxav',
     .    'monthly integral of eastward isopycnic mass flux','Sv')
        call out3cdf(ncid1,idm,jdm,kdm,ufxavp,time,'ufxavp',
     .    'monthly integral of southward isobaric mass flux','Sv')
        call out3cdf(ncid1,idm,jdm,kdm,vfxavp,time,'vfxavp',
     .    'monthly integral of eastward isobaric mass flux','Sv')
        call out3cdf(ncid1,idm,jdm,kdm,diaflx,time,
     .    'diaflx','monthly integral of diaflx (interlayer mass flux per
     . unit area)','m')
!       call out3cdf(ncid1,idm,jdm,kdm,vctyav,time,
!         'visc','monthly vertical viscosity','m^2/sec')
!       call out3cdf(ncid1,idm,jdm,kdm,diftav,time,
!         'dfft','monthly vert.temperature diffusivity','m^2/sec')
!       call out3cdf(ncid1,idm,jdm,kdm,difsav,time,
!         'dffs','monthly vert.scalar diffusivity','m^2/sec')
!       call out2cdf(ncid1,idm,jdm,radflav,time,
!    . 'radfx','monthly surface net radiative flux, pos.down','W/m^2')
!       call out2cdf(ncid1,idm,jdm,snsibav,time,
!    . 'snsib','monthly surface sensible heat flux, pos.down','W/m^2')
!       call out2cdf(ncid1,idm,jdm,latntav,time,
!    . 'latnt','monthly surface latent heat flux, pos.down','W/m^2')
        call errhandl (nf90_inq_varid(ncid1,'time',i))
        call errhandl (nf90_put_att  (ncid1,i,'calendar','NOLEAP'))
        print '(2a)','closing ',trim(flnm_nc)
        call errhandl (nf90_close (ncid1))
      end if            ! ncout
!<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
#if (defined TRACERS_OceanBiology) || (defined TRACERS_AGE_OCEAN) \
 || (defined TRACERS_OCEAN_WATER_MASSES) || (defined TRACERS_ZEBRA)
!     call obio_archyb(nn,dpav,temav,salav,th3av,dpmxav,oiceav)
#endif

      do 60 j=1,jj
      do 601 l=1,isp(j)
      do 601 i=ifp(j,l),ilp(j,l)
c      eminpav(i,j)=0.
c      surflav(i,j)=0.
c      salflav(i,j)=0.
c      brineav(i,j)=0.
c       tauxav(i,j)=0.
c       tauyav(i,j)=0.
c       pbavav(i,j)=0.
c       dpmxav(i,j)=0.
c       sfhtav(i,j)=0.
c       oiceav(i,j)=0.
c
#ifdef TRACERS_OceanBiology
        diag_counter =0
        ao_co2fluxav(i,j)=0.
        pco2av(i,j)=0.
        pp2tot_dayav(i,j)=0.
        cexpav(i,j)=0.
#ifdef TRACERS_Alkalinity
        caexpav(i,j)=0.
#endif
#endif

 601  continue
c
c     do 60 k=1,kk
c     do 602 l=1,isp(j)
c     do 602 i=ifp(j,l),ilp(j,l)
c     uav(i,j,k)=0.
c     vav(i,j,k)=0.
c     dpuav(i,j,k)=0.
c     dpvav(i,j,k)=0.
c     dpav (i,j,k)=0.
c     temav(i,j,k)=0.
c     salav(i,j,k)=0.
c     th3av(i,j,k)=0.
c     uflxav(i,j,k)=0.
c     vflxav(i,j,k)=0.
c     ufxavp(i,j,k)=0.
c     vfxavp(i,j,k)=0.
c602  diaflx(i,j,k)=0.
 60   continue
c
      return
      end
c
      subroutine r8tor4(real8,real4)
c
      USE HYCOM_DIM_GLOB
      implicit none
      integer i,j
c
      real real8(idm,jdm)
      real*4 real4(idm,jdm)
c
      do 1 j=1,jdm
      do 1 i=1,idm
 1    real4(i,j)=real8(i,j)
c
      return
      end
c
c
c> Revision history:
c>
c> June 2001 - corrected sign error in diaflx
c> Feb. 2005 - added multiple tracer capability
c> June 2005 - reordered output fields, added time-av'ged mxlyr & ice temp/th
