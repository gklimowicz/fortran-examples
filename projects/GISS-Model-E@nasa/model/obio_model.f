#include "rundeck_opts.h"

#define NO_REDIAG_OCNSTATE

      subroutine obio_model(mm)

!@sum  OBIO_MODEL is the main ocean bio-geo-chem routine 
!@auth Natassa Romanou/Watson Gregg

      USE Dictionary_mod
      USE obio_dim
      USE obio_incom
      use bio_inicond_mod, only: bio_inicond
      USE obio_forc, only: solz,tirrq,Ed,Es
     .                    ,rmud,atmFe,stream_atmFe,avgq,sunz
     .                    ,wind
     .                    ,alk
     .                    ,tirrq3d
     .                    ,surfN
#ifdef prescribe_o2sf
     .                    ,o2fc
#endif
#ifdef OBIO_RUNOFF
     .                    ,river_runoff
#endif
      USE obio_com,  only: dobio,gcmax,day_of_month,hour_of_day
     .                    ,temp1d,dp1d,obio_P,det,car,avgq1d
     .                    ,gcmax1d,atmFe_ij,daily_atmFe,covice_ij
     .                    ,P_tend,D_tend,C_tend,saln1d
     .                    ,pCO2_ij,p1d,wsdet,pHsfc
     .                    ,rhs,alk1d
     .                    ,tzoo,tfac,rmuplsr,rikd,wshc,Fescav
     .                    ,tzoo2d,tfac3d,rmuplsr3d,rikd3d
     .                    ,wshc3d,Fescav3d 
     .                    ,acdom,pp2_1d,pp2tot_day,pp2diat_day
     .                    ,pp2chlo_day,pp2cyan_day,pp2cocc_day
     .                    ,tot_chlo,acdom3d
     .                    ,itest,jtest
     .                    ,obio_ws,trmo_unit_factor
     .                    ,cexp,flimit,kzc
     .                    ,rhs_obio,chng_by,Kpar,Kpar_em2d,Edz,Esz,Euz
     .                    ,delta_temp1d,sday
     .                    ,use_qus,num_tracers        
     .                    ,rho_water,rho1d
#ifdef TRACERS_bio_O2 
     .                    ,errchk1,errchk2   !@PLdbg added test for NaNs in o2
#endif
#ifdef TOPAZ_params
     .                    ,ca_det_calc1d
#endif
#ifdef TRACERS_Ocean_O2 
#ifdef TRACERS_bio_O2
     .                    ,o21d,pO2_ij
#ifdef restart_add_o2
     .                    ,o2rst
#endif
#endif
#ifdef TRACERS_abio_O2
     .                    ,abo21d,pabO2_ij
#endif
#endif
#ifdef TRACERS_degC
     .                    ,ndegC1d,cexpdeg
#endif
#ifdef OBIO_RUNOFF
     .                    ,rnitrconc_loc,rdicconc_loc,rdocconc_loc
     .                    ,rsiliconc_loc,rironconc_loc,rpocconc_loc
     .                    ,ralkconc_loc
     .                    ,rnitrconc,rdicconc,rdocconc
     .                    ,rsiliconc,rironconc,rpocconc
     .                    ,ralkconc
#endif
      use obio_com, only: caexp
      use obio_com, only: build_ze

#ifdef STANDALONE_OCEAN 
      USE obio_forc, only: Eda,Esa,Eda2,Esa2
#endif

      use runtimecontrols_mod, only: tracers_alkalinity
      use obio_com, only: tracer,nstep0
      use obio_com, only: obio_deltat
      USE obio_diag, only : ij_pCO2,ij_dic,ij_nitr,ij_diat
     .                 ,ij_amm,ij_sil,ij_chlo,ij_cyan,ij_cocc,ij_herb
     .                 ,ij_doc,ij_iron,ij_alk,ij_Ed,ij_Es,ij_pp,ij_dayl
     .                 ,ij_cexp,ij_lim,ij_sink,ij_setl,ij_ndet,ij_xchl
     .                 ,ij_sunz,ij_solz
     .                 ,ij_pp1,ij_pp2,ij_pp3,ij_pp4
     .                 ,ij_rhs,ij_flux,ij_fca      !,kobio_ij!@PL kobio_ij-> number of obio_ij diagnostics
#ifdef TRACERS_Ocean_O2 
#ifdef TRACERS_bio_O2
     .                 ,ij_o2,ij_oflx,ij_po2
#endif
#ifdef TRACERS_abio_O2
     .                 ,ij_abo2,ij_aboflx,ij_pabo2
#endif
#endif
#ifdef TRACERS_degC
     .                 ,ij_degC,ij_cexpdeg
#endif
#ifdef OBIO_RUNOFF
     .                 ,ij_rnitrconc !,ij_rnitrmflo
     .                 ,ij_rdicconc,ij_rdocconc,ij_rsiliconc
     .                 ,ij_rironconc,ij_rpocconc,ij_ralkconc
#endif

      USE obio_diag, only : oijl=>obio_ijl,ijl_avgq,ijl_kpar,
     .                            ijl_kpar_em2d,ijl_dtemp,ijl_rhs3
     .                           ,ijl_wss,ijl_wsdet
     .                           ,ijl_pp
     .                           ,ijl_lim1,ijl_lim2,ijl_lim3
     .                           ,ijl_lim4,ijl_lim5!hemis_obio_ij !@PL hemis_obio_ij-> hemisphere and global inventories
#ifdef TRACERS_bio_O2 
     .                           ,ijl_cprod
     .                           ,ijl_cdet
     .                           ,ijl_cdoc
     .                           ,ijl_caresp
     .                           ,ijl_chresp
     .                           ,ijl_oprodam
     .                           ,ijl_oprodnit
     .                           ,ijl_odet
     .                           ,ijl_odoc
     .                           ,ijl_oaresp
     .                           ,ijl_ohresp
#endif
      use ocalbedo_mod, only: ocalbedo
      USE MODEL_COM, only: modelEclock
     . ,itime,iyear1,aMON,
     . xlabel,lrunid !,dtsrc
      use JulianCalendar_mod, only: jdendofm
      use TimeConstants_mod, only: HOURS_PER_DAY
      USE FILEMANAGER, only: openunit,closeunit,file_exists
      USE timestream_mod, only: read_stream
      USE obio_com, only : co2flux
#ifdef TRACERS_Ocean_O2 
#ifdef TRACERS_bio_O2
     .               ,o2flux
#endif
#ifdef TRACERS_abio_O2
     .               ,abo2flux
#endif
#endif
      use obio_com, only: ze
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,DIST_GRID,GLOBALSUM
     .                           ,broadcast,sumxpe
      use TimerPackage_mod
      use OCEAN, only: oxyp
      use obio_diag, only: reset_obio_diag
      use runtimecontrols_mod, only: constco2

#ifdef OBIO_ON_GISSocean
      USE MODEL_COM,  only : nstep=>itime,itimei,dtsrc
      USE OCEAN,       only : oLON_DG,oLAT_DG,dlatm,DTS
      USE CONSTANT,   only : grav
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES,   only : kdm=>lmo,dzo
      USE OFLUXES,    only : oice=>oRSI,oAPRESS,ocnatm
#ifdef OBIO_RUNOFF
      USE OFLUXES, only:  oFLOWO
#endif

#ifdef NO_REDIAG_OCNSTATE
      ! Deliberately not taking layer thickness information, for reasons
      ! related to upcoming commits.
      use ocean, only : t3d,s3d,r3d
#else
      USE OCEAN, only : g0m,s0m
#endif
      USE OCEAN,      only : mo=>motr,dxypo,ip=>focean,lmm
     .                      ,trmo,txmo,tymo,tzmo
     .                      ,txxmo,txymo,tzxmo,tyymo,tyzmo,tzzmo
      use ocn_tracer_vector_mod, only:
     &                           vector_ocn_tracer_entry=>vector
      USE OCN_TRACER_COM, only :n_abioDIC,add_ocn_tracer,tracerlist
      USE sw2ocean, only : lsrpd,fsr
      use OCEANRES, only : IM=>IMO,JM=>JMO
      USE EXCHANGE_TYPES, only : iceocn_xchng_vars
      use obio_diag, only: oij=>obio_ij,ij_ph
#ifdef TRACERS_Alkalinity
     .  ,ij_co3
#endif

#else   /* HYCOM */
      USE hycom_dim, only: ogrid,im=>idm,jm=>jdm,kdm,ntrcr
      USE hycom_arrays, only: tracer_h=>tracer,dpinit,temp,saln,oice
     .                            ,p,dpmixl,latij,lonij,scp2
      USE  hycom_arrays_glob, only: latij_glob=>latij,lonij_glob=>lonij
      USE hycom_scalars, only: nstep,onem,time,dtsrc=>baclin,dts=>baclin
      USE obio_com, only: ao_co2fluxav_loc,
     .     pCO2av_loc,pp2tot_dayav_loc,cexpav_loc,caexpav_loc,
     .     pp2diat_dayav_loc,pp2chlo_dayav_loc,pp2cyan_dayav_loc,
     .     pp2cocc_dayav_loc,pHav,pHav_loc
      USE obio_com, only: diag_counter,phav_loc
      use GEOM, only: dlatm
      use hycom_atm, only : ocnatm
#endif

      implicit none

      !molecular weights (gr/mole)
!@PL directives added for Alkalinity/O2 options

      REAL*4, parameter  :: obio_tr_mm(ntrac)= (/ 14.,   !nitrate
     &     14.,      !ammonium
     &     28.055,   !silicate
     &     55.845,   !iron
     &     1., 1., 1., 1., 1.,   !phyto and zooplankton
     &     14.,      !ndet
     &     28.055,   !sdet
     &     55.845,   !idet
     &     12.,      !DOC
     &     12.      !DIC
#ifdef TRACERS_Alkalinity
     &    ,1.       !Alk
#endif
#ifdef TRACERS_Ocean_O2 
#ifdef TRACERS_bio_O2
     &     ,32.    !O2 !@PL I think its better to use MW of O2 = 32 gr O /mole
#endif
#ifdef TRACERS_abio_O2
     &     ,32.
#endif
#endif
#ifdef TRACERS_degC
     &     ,14.     !ndet, nondegradable detritus
#endif
     &      /)

!@PL saturation O2 for restarts from files w/out O2
#ifdef TRACERS_abio_O2 
      real O2sat
#endif

      integer i,j,k,l,km,mm,JMON,khem
      integer ihr,ichan,iyear,nt,ihr0,lgth,kmax
      integer ll,ilim
      real    tot,dummy(6),dummy1
      real    rod(nlt),ros(nlt)
#ifdef OBIO_ON_GISSocean
#ifndef NO_REDIAG_OCNSTATE
      Real*8,External   :: VOLGSP
      real*8 temgs,g,s,temgsp,pres
#endif
      real*8 time,dtr,ftr,totirondep
      REAL*8, DIMENSION(IM,ogrid%J_STRT_HALO:ogrid%J_STOP_HALO)::
     .        iron_atmdep
#ifdef OBIO_RUNOFF
      REAL*8, DIMENSION(IM,ogrid%J_STRT_HALO:ogrid%J_STOP_HALO)::
     .        nitr_river,sili_river,iron_river,poc_river
     .       ,dic_river,doc_river,alk_river
      REAL*8 :: totnitrrivr,totsilirivr,totironrivr,totpocrivr
     .         ,totdocrivr,totdicrivr,totoriver
#endif
#else
      integer :: n_abioDIC,
     &           lmm(im,ogrid%j_strt_halo:ogrid%j_stop_halo)
      real :: dxypo(jm),
     &        MO(im,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm),
     &        trmo(im,ogrid%j_strt_halo:ogrid%j_stop_halo,
     &             kdm,ntrcr),
     &        ip(ogrid%i_strt_halo:ogrid%i_stop_halo,
     &                      ogrid%j_strt_halo:ogrid%j_stop_halo)     
#endif

      integer :: idx_co2
      logical vrbos,noon,errcon
      integer :: year, month, dayOfYear, date, hour
      real, allocatable, dimension(:), save :: obio_lambdas,eda_frac,
     &                                         esa_frac
      integer :: iu_bio
      logical, save :: initialized=.false.
      real, save :: atmco2=-1.
      real :: rlon2D(im,jm),rlat2D(im,jm),ddxypo,mmo
      real :: SDIC!(num_tracers) PLdbg
      integer :: i_0,i_1,j_0,j_1
      integer :: n_co2n
      real :: oij_pH
#ifdef TRACERS_Alkalinity
      real :: oij_co3
#endif

!@PL
!      real*8, dimension(:,:), allocatable :: shnh_loc,shnh
!@PL

      if(.not.dobio) return

#ifdef OBIO_ON_GISSocean
      num_tracers=tracerlist%getsize()
#else
      num_tracers=ntrcr
#endif

      call start(' obio_model')
!--------------------------------------------------------
      diagno_bio=.false.
c
      call modelEclock%get(year=year, month=month, date=date,
     .  hour=hour, dayOfYear=dayOfYear)

      if (JDendOfM(month).eq.dayOfYear.and.hour.eq.12) then
          if (mod(nstep,2).eq.0)
     .    diagno_bio=.true. ! end of month,mid-day
      endif

!Cold initialization

      call start(' obio_init')
      !nstep0=0 : cold initialization
      !nstep0>0 : warm initialization, is the timestep of current restart run
      !nstep    : current timestep   

#ifndef OBIO_QUIET_MODE
      if (AM_I_ROOT()) print*, 'nstep,nstep0 =',
     .         nstep,nstep0
#endif

      call build_ze
      if (nstep0==0) then
        if (AM_I_ROOT()) print*, 'COLD INITIALIZATION....'

        if (AM_I_ROOT()) write(*,'(a)')'BIO:Ocean Biology starts ....'

!PLdbg
!     print*, 'I_0,I_1=',
!     .         ogrid%I_STRT,ogrid%I_STOP

        call obio_init(kdm,dtsrc,ogrid,dlatm)

      !tracer array initialization.
      !note: we do not initialize obio_P,det and car

#ifndef OBIO_ON_GISSocean
        ao_co2fluxav_loc  = 0.
        pCO2av_loc = 0
        pp2tot_dayav_loc = 0
        pp2diat_dayav_loc = 0
        pp2chlo_dayav_loc = 0
        pp2cyan_dayav_loc = 0
        pp2cocc_dayav_loc = 0
        cexpav_loc = 0
        caexpav_loc = 0
        pHav_loc = 0
#endif
        if (AM_I_ROOT())
     .  print*,nstep,'calling bioinit'

#ifdef OBIO_ON_GISSocean
       rlon2D=transpose(spread(oLON_DG(:,1),DIM=1,NCOPIES=jm))
       rlat2D=spread(oLAT_DG(:,1),DIM=1,NCOPIES=im)
#else
       rlon2D=lonij(:,:,3)
       rlat2D=latij(:,:,3)
       !lmm=0d0 this does not work because for GISSocean is
!passed from a module.  
#endif 
       
!PLdbg
!     print*, 'I_0,I_1=',
!     .         ogrid%I_STRT,ogrid%I_STOP
       
!NOTE: ip=>focean is real in GISSocean but ip is integer in
!      hycom. It is passed as a real 2D array

!NOT FOR HYCOM: mo, DXYPO, trmo,n_abioDIC,im,jm, number traces and lmm          passed to subroutine
        call obio_bioinit(kdm,n_abioDIC,DXYPO,ogrid,
     &                    MO,trmo,num_tracers,
     &                    im,jm,rlon2D,rlat2D,ip,lmm)


      endif   !cold restart



!Warm initialization

#ifdef OBIO_ON_GISSocean
      time = float(nstep)
#endif
      if ((nstep0>0).and..not.initialized) then
         write(*,'(a)')'For restart runs.....'
         write(*,'(a,2i9,f10.3)')
     .            'nstep0,nstep,time=',nstep0,nstep,time
         call obio_init(kdm,dtsrc,ogrid,dlatm)
      if (AM_I_ROOT())   print*,'WARM INITIALIZATION'
      endif !for restart only
      call stop(' obio_init')

      i_0=ogrid%I_STRT
      i_1=ogrid%I_STOP
      j_0=ogrid%J_STRT
      j_1=ogrid%J_STOP

#ifdef prescribe_o2sf
      allocate(o2fc(i_0:i_1,j_0:j_1,kdm))

      call bio_inicond('oxygen_inicond',
     &     o2fc(:,:,:),kdm,im,ogrid,ip,lmm)
#endif


!@PL if using a warm ocean tracer restart w/out ocean O2, this 
!@PL will initialize O2 from GLODAP climatology
#ifdef restart_add_o2
      call bio_inicond('oxygen_inicond',
     &     o2rst(:,:,:),kdm,im,ogrid,ip,lmm)
#endif

#ifndef OBIO_QUIET_MODE
      if (AM_I_ROOT()) then
#ifdef OBIO_ON_GISSocean
      write(*,'(a,2i5,2e12.4)')'TEST POINT at: ',itest,jtest,
     .      oLAT_DG(jtest,2),oLON_DG(itest,2)
#else
      write(*,'(a,2i5,2e12.4)')'TEST POINT at: ',itest,jtest,
     .      latij_glob(itest,jtest,3),lonij_glob(itest,jtest,3)
#endif
       endif
#endif /* not OBIO_QUIET_MODE */
       call sync_param( "solFe", solFe)
#ifndef OBIO_QUIET_MODE
       if (AM_I_ROOT()) print*, 'solfe=',solFe

       if (AM_I_ROOT()) print*, 'using QUS?=',USE_QUS
#endif /* not OBIO_QUIET_MODE */
!--------------------------------------------------------

       day_of_month=date
       hour_of_day =hour

       !if (time.kmod(nstep,24*30).eq.0) then
       !   day_of_month=1
       !else
       !   day_of_month=day_of_month+1
       !endif

       !if (mod(nstep,24).eq.0) then
       !  hour_of_day=1
       !else
       !  hour_of_day=hour_of_day+1
       !endif

#ifndef OBIO_QUIET_MODE
      if (AM_I_ROOT()) then
       write(*,'(a,i15,1x,f9.3,2x,3i5)')
     .    'BIO: nstep,time,day_of_month,hour_of_day,dayOfYear=',
     .    nstep,time,day_of_month,hour_of_day,dayOfYear
      endif
#endif /* not OBIO_QUIET_MODE */
        !ihr0 = int(hour_of_day/2)
         ihr0 = nint((hour_of_day+1)/2.)


      if (diagno_bio) then
      endif  !diagno_bio

#ifndef OBIO_QUIET_MODE
#ifdef OBIO_ON_GISSocean
      if (AM_I_ROOT())
     .   write(*,'(/,a,2i5,2e12.4)')'obio_model, test point=',
     .      itest,jtest,oLON_DG(itest,1),oLAT_DG(jtest,1)
#endif
#endif /* not OBIO_QUIET_MODE */
      !print out tracer integrals just before main loop

#ifndef OBIO_ON_GISSocean     /* HYCOM only */
      diag_counter=diag_counter+1
#endif

      if ((dayofyear==1+jdendofm(month-1)).and.
     &       modeleclock%isbeginningofday()) call reset_obio_diag

!#ifdef OBIO_RUNOFF
!!compute riverflow globalsum
      !!don't divide by area here in order to get kg/m2
      !nitr_river = rnitrconc_loc(:,:)*oFLOWO(:,:)
      !sili_river = rsiliconc_loc(:,:)*oFLOWO(:,:)
      !iron_river = rironconc_loc(:,:)*oFLOWO(:,:)
       !poc_river =  rpocconc_loc(:,:)*oFLOWO(:,:)
      !doc_river =  rdocconc_loc(:,:)*oFLOWO(:,:)
       !dic_river =  rdicconc_loc(:,:)*oFLOWO(:,:)
!!#ifdef TRACERS_Alkalinity
!!       alk_river = ralkconc_loc(:,:)*oFLOWO(:,:)
!!#endif
      !call GLOBALSUM(ogrid, nitr_river, totnitrrivr, ALL=.true.)
      !call GLOBALSUM(ogrid, sili_river, totsilirivr, ALL=.true.)
      !call GLOBALSUM(ogrid, iron_river, totironrivr, ALL=.true.)
      !call GLOBALSUM(ogrid, poc_river, totpocrivr,  ALL=.true.)
      !call GLOBALSUM(ogrid, doc_river, totdocrivr,  ALL=.true.)
      !call GLOBALSUM(ogrid, dic_river, totdicrivr,  ALL=.true.)
      !call GLOBALSUM(ogrid, oFLOWO, totoriver,  ALL=.true.)
      !if (AM_I_ROOT()) then
        !write(6,*) "Total Nitr River Mass into Ocean: ",totnitrrivr
        !write(6,*) "Total Sili River Mass into Ocean: ",totsilirivr
        !write(6,*) "Total Iron River Mass into Ocean: ",totironrivr
        !write(6,*) "Total POC  River Mass into Ocean: ",totpocrivr
        !write(6,*) "Total DOC  River Mass into Ocean: ",totdocrivr
        !write(6,*) "Total DIC  River Mass into Ocean: ",totdicrivr
        !write(6,*) "Total River Mass into Ocean: ",totoriver
      !end if
      !call broadcast(ogrid, totnitrrivr)
      !call broadcast(ogrid, totsilirivr)
      !call broadcast(ogrid, totironrivr)
      !call broadcast(ogrid, totpocrivr)
      !call broadcast(ogrid, totdocrivr)
      !call broadcast(ogrid, totdicrivr)
      !call broadcast(ogrid, totoriver)
!#endif

!compute dust iron mass  globalsum
!      iron_atmdep = atmFe(:,:,month)*oXYP(:,:)
!      call GLOBALSUM(ogrid,  iron_atmdep, totirondep,  ALL=.true.)
!      if (AM_I_ROOT())
!     .  write(6,*) "Total Iron Dust deposition into Ocean: ",totirondep
!      call broadcast(ogrid, totirondep)


      call start('  obio main loop')
!@PL allocate north/south hem arrays

!      allocate(shnh_loc(2,kobio_ij),shnh(2,kobio_ij))
!      shnh_loc = 0.
!@PL

!PLdbg
!      print*, 'I_0,I_1,J_0,J_1=',
!     .         i_0,i_1,j_0,j_1

      ocnatm%chl_defined=.true.
       do j=j_0,j_1
       do i=i_0,i_1
!PLdbg
!       print*, 'I,J,I_1,J_1=',I,J,I_1,J_1
!PLdbg
       dp1d(:) = 0.
       IF(ip(I,J)==0) cycle

       vrbos=.false.
       if (i.eq.itest.and.j.eq.jtest) vrbos=.true.

       !fill in reduced rank arrays
       !!covice_ij=covice(i,j)  !for standalone hycom
       covice_ij=oice(i,j)      !for modelE-hycom
       !!!!pCO2_ij=atm%gtracer(atm%n_co2n,i,j)
        if (atmco2<0.) then    ! uninitialized
          if (constco2) then
            call get_param("atmCO2", atmco2)
            if (AM_I_ROOT())print*, 'atmCO2=', atmco2
          else
            atmco2=0.
          endif
        endif
       call compute_pco2_online(nstep,i,j,atmco2,
     .           temp1d(1),saln1d(1),car(1,2),alk1d(1),
     .           obio_P(1,1),obio_P(1,3),(1.-covice_ij),
     .           pCO2_ij,pHsfc,vrbos)
     
#ifdef OBIO_ON_GISSocean
#ifndef NO_REDIAG_OCNSTATE
       pres = oAPRESS(i,j)    !surface atm. pressure
#endif
       do k=1,lmm(i,j)

#ifdef NO_REDIAG_OCNSTATE /* get derived fields from host rather than recomputing */
         temp1d(k) = t3d(k,i,j)         ! in-situ temperature
         saln1d(k) = s3d(k,i,j)*1000.   ! convert to psu (eg. ocean mean salinity=35psu)
         rho_water = r3d(k,i,j)         ! in-situ density
         rho1d(k) = rho_water
#else
         pres=pres+MO(I,J,k)*GRAV*.5
         g=G0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
         s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
!!!!     temp1d(k)=TEMGS(g,s)           !potential temperature
         temp1d(k)=TEMGSP(g,s,pres)     !in situ   temperature
          saln1d(k)=s*1000.             !convert to psu (eg. ocean mean salinity=35psu)
          !dp1d(k)=dzo(k)               !thickenss of each layer in meters
         rho_water = 1d0/VOLGSP(g,s,pres)
         rho1d(k) = rho_water
!!!!!!   rho_water = 1035.
         !add missing part of density to get to the bottom of the layer
         !now pres is at the bottom of the layer
         pres=pres+MO(I,J,k)*GRAV*.5
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! array tracer has units [mole]/m3. convert to/from trmo with units kg as follows:
! trmo = tracer * MB*1e-3/rho_water     * mo * dxypo
! note that occassionally tracer is in milimoles and othertimes in micromoles
! trmo_unit_factor below has units m^3 and depends on which 
! layer you are: trmo = tracer * factor,  tracer=trmo/factor
! txmo,tzmo,tymo are all computed based on trmo
! this should be done at every timestep except when COLD INITIALIZATION
! because then trmo=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         dp1d(k)=MO(I,J,K)/rho_water   !local thickenss of each layer in meters
         if(vrbos.and.k.eq.1)write(*,'(a,4e12.4)')
     .             'obio_model,t,s,p,rho= '
     .             ,temp1d(k),saln1d(k),dp1d(k),rho_water

         do nt=1,ntrac
         
                trmo_unit_factor(k,nt) =  1d-3*1d-3*obio_tr_mm(nt)        ! milimoles/m3=> kg/m3
     .                      *  MO(I,J,k)*DXYPO(J)/rho_water               ! kg/m3 => kg

         !iron and idet
         if (nt.eq.nnut.or.nt.eq.ntyp+ndet) 
     .          trmo_unit_factor(k,nt) =  1d-6*1d-3*obio_tr_mm(nt)        ! nanomoles/lt=> kg/m3
     .                      *  MO(I,J,k)*DXYPO(J)/rho_water               ! kg/m3 => kg

         !ndet
         if (nt.eq.ntyp+1)
     .          trmo_unit_factor(k,nt) =  1d-6 *1d-3/1d-3                 ! micro-grC/lt -> kg/m3
     .                      *  MO(I,J,k)*DXYPO(J)/rho_water               ! kg/m3 => kg
#ifdef TRACERS_degC
         !nondegradable carbon
         if (nt.eq.ndimndegC)
     .          trmo_unit_factor(k,nt) =  1d-6 *1d-3/1d-3                 ! micro-grC/lt -> kg/m3
     .                      *  MO(I,J,k)*DXYPO(J)/rho_water               ! kg/m3 => kg
#endif

         !phyto and zooplankton
         if (nt.ge.nnut+1.and.nt.le.ntyp) 
     .          trmo_unit_factor(k,nt) =  
     .                          cchlratio * 1d-6                          ! miligr,chl/m3=> kg,C/m3
     .                       *  MO(I,J,k)*DXYPO(J)/rho_water              ! kg/m3 => kg

#ifdef TRACERS_Alkalinity
           if (nt.eq.ntyp+ndet+ncar+1)    !factor for alkalinity
     .         trmo_unit_factor(k,nt) = 1d-6*1d-3*obio_tr_mm(nt)      ! umol/kg=micro-mol/kg=> kg,trac/kg,air
     .                                *  MO(I,J,k)*DXYPO(J)           ! kg,trac/kg,air=> kg,trac

#ifdef TOPAZ_params
           if (nt.eq.ntyp+ndet+ncar+2)    !factor for ca_det
     .         trmo_unit_factor(k,nt) = 1d-6*1d-3*obio_tr_mm(nt)      ! umol/kg=micro-mol/kg=> kg,trac/kg,air
     .                                *  MO(I,J,k)*DXYPO(J)           ! kg,trac/kg,air=> kg,trac
#endif
#endif
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
          if (nt.eq.ndimo2)    !factor for oxygen 
     .         trmo_unit_factor(k,nt) = 1d-3*1d-3*obio_tr_mm(nt)      ! mmol/kg=milli-mol/kg=> kg,trac/kg,water
     .                                *  MO(I,J,k)*DXYPO(J)           ! kg,trac/kg,water=> kg,trac
#endif
#ifdef TRACERS_abio_O2
          if (nt.eq.ndimabo2)    !factor for oxygen 
     .         trmo_unit_factor(k,nt) = 1d-3*1d-3*obio_tr_mm(nt)      ! mmol/kg=milli-mol/kg=> kg,trac/kg,water
     .                                *  MO(I,J,k)*DXYPO(J)           ! kg,trac/kg,water=> kg,trac
#endif
#endif

           if (nstep0>0) then
              tracer(i,j,k,nt) = trmo(i,j,k,nt) / trmo_unit_factor(k,nt)
           endif

#ifdef restart_add_o2
          if (.not.initialized) then
#ifdef TRACERS_bio_O2
          if (nt.eq.ndimo2) then
             tracer(i,j,k,nt) = o2rst(i,j,k) !@PL add O2 from GLODAP to tracer array if restart from run w/out O2
          endif
#endif
#ifdef TRACERS_abio_O2
          if (nt.eq.ndimabo2) then
             call init_abo2(temp1d(k),saln1d(k),oAPRESS(i,j)
     &                       ,ocnatm%QSAVG(i,j),O2sat)   !@PL mmol/kg 
             tracer(i,j,k,nt) = O2sat
          endif
#endif
          endif
#endif

#ifdef restart_degC
          if (.not.initialized) then
#ifdef TRACERS_degC
          if (nt.eq.ndimndegC) then
             tracer(i,j,k,nt) = 0 !@PL initalize nondegradable carbon as 0
          endif
#endif
          endif
#endif


         enddo


#ifdef IncreaseN
       tracer(i,j,1,1) = 100.d0 * tracer(i,j,1,1)
#endif
#ifdef Relax2SurfN
#ifdef Relax2maxSurfN
       tracer(i,j,1,1) = max(tracer(i,j,1,1),surfN(i,j))
#else
       tracer(i,j,1,1) = surfN(i,j)
#endif
#endif

#else
       trmo_unit_factor=1
       if (nstep0>0) tracer(i,j,:,:)=tracer_h(i,j,:,:)/trmo_unit_factor
       do k=1,kdm
         km=k+mm
         temp1d(k)=temp(i,j,km)
         saln1d(k)=saln(i,j,km)
         dp1d(k)=dpinit(i,j,k)/onem
#endif
         avgq1d(k)=avgq(i,j,k)
         gcmax1d(k)=gcmax(i,j,k)
         tirrq(k)=tirrq3d(i,j,k)
#ifdef TRACERS_Alkalinity
         nt=1
         alk1d(k)=tracer(i,j,k,ntyp+ndet+ncar+1)
#ifdef TOPAZ_params
         nt=2
         ca_det_calc1d(k)=tracer(i,j,k,ntyp+ndet+ncar+nt)
#endif
#else
         !NOT for INTERACTIVE alk
         alk1d(k)=alk(i,j,k)
#endif

!@PL o21d before bio update
#ifdef TRACERS_Ocean_O2 
#ifdef TRACERS_bio_O2
         o21d(k)=tracer(i,j,k,ntyp+ndet+ncar+nalk+no2)
#endif
#ifdef TRACERS_abio_O2
         abo21d(k) = tracer(i,j,k,ntyp+ndet+ncar+nalk+no2+nabo2)
#endif
#endif

              !----daysetbio/daysetrad arrays----!
         tzoo=tzoo2d(i,j)
         tfac(k)=tfac3d(i,j,k)
         do nt=1,nchl
           rmuplsr(k,nt)=rmuplsr3d(i,j,k,nt)
           rikd(k,nt)=rikd3d(i,j,k,nt)
         enddo
         wshc(k)=wshc3d(i,j,k)
         Fescav(k)=Fescav3d(i,j,k)
         do nt=1,nlt
           acdom(k,nt)=acdom3d(i,j,k,nt)
         enddo
              !----daysetbio arrays----!
         do nt=1,ntyp
           obio_P(k,nt)=tracer(i,j,k,nt)
         enddo
         do nt=1,ndet
           det(k,nt)=tracer(i,j,k,ntyp+nt)
         enddo
         do nt=1,ncar
           car(k,nt)=tracer(i,j,k,ntyp+ndet+nt)
         enddo

#ifdef TRACERS_degC
         ndegC1d(k) = tracer(i,j,k,ndimndegC)
#endif
       enddo  !k=1,kdm or lmm

!@PL prescribe surface o2
#ifdef prescribe_o2sf
         o21d(1) = o2fc(i,j,1)
#endif

       p1d(1)=0.
       do k=2,kdm+1
          p1d(k)=p1d(k-1)+dp1d(k-1)    !in meters
       enddo

#ifdef TRACERS_bio_O2 

      do k=1,kdm

!       if (vrbos) then
!          write(6,'(a,3i7,2e12.4)')'obio_o2(prebio):',
!     .      nstep,i,j,p1d(k),o21d(k)
!
!        endif

      if (ISNAN(o21d(k))) then
        errchk1=1
          write(6,'(a,4i7,3e12.4)')'obio_o2(prebionan):',
     .      nstep,i,j,errchk1,p1d(k),o21d(k),car(k,2)

      endif
      enddo
#endif

      !if(vrbos) write(*,'(a,15e12.4)')'obio_model, strac conc:',
      if(vrbos) write(*,*)'obio_model, strac conc:',
     .    obio_P(1,:),det(1,:),car(1,:)

#ifdef TOPAZ_params
      if(vrbos) write(*,'(a,2e12.4)')'obio_model, sALK  conc:',
     . alk1d(1),ca_det_calc1d(1)
#endif

#ifdef OBIO_ON_GISSocean
      kmax = lmm(i,j)
cdiag if(vrbos)write(*,'(a,4i5)')'nstep,i,j,kmax= ',
cdiag.         nstep,i,j,kmax

#else

! --- find number of non-massless layers in column
      do kmax=1,kdm
c       if (dp1d(kmax).lt.1.e-2) exit
        if (dp1d(kmax).lt.1.) exit
      end do
      kmax=kmax-1
cdiag write(*,'(a,4i5)')'nstep,i,j,kmax= ',nstep,i,j,kmax

      !ensure that all massless layers have same concentration
      !as last mass one
      do k=kmax+1,kdm
       do nt=1,ntyp
        if(obio_P(kmax,nt).le.0.)obio_P(kmax,nt)=1.e-8
        obio_P(k,nt)=obio_P(kmax,nt)
       enddo
       do nt=1,ndet
        det(k,nt)=   det(kmax,nt)
       enddo
       do nt=1,ncar
        car(k,nt)=   car(kmax,nt)
       enddo
#ifdef TRACERS_degC
      ndegC(k) = ndegC(kmax,nt)
#endif
      enddo
#endif

      !find layer index for zc, max depth of sinking phytoplankton
      kzc = 1
      do k=kmax+1,1,-1
           if (p1d(k).gt.zc) kzc = k
      enddo
      if (kzc.lt.1) kzc=1
      if (kzc.gt.kmax) kzc=kmax

      if (vrbos) write(*,*) 'compensation depth, kzc = ',kzc

#ifdef STANDALONE_OCEAN
       !OASIM spectral irradiance data just above the surface
       !Eda and Esa fields have ihr=1:12, ie every 2hrs
       JMON = modelEclock%getMonth()
       do ichan=1,nlt
       Eda2(ichan,ihr0)=Eda(i,j,ichan,ihr0,JMON)
       Esa2(ichan,ihr0)=Esa(i,j,ichan,ihr0,JMON)
       enddo
!      if (vrbos) then
!      write(*,'(a,6i5,2e12.4)') 'obio_model, Eds:',
!    .      nstep,i,j,7,ihr0,JMON,Eda(i,j,7,ihr0,JMON)
!    .     ,Esa(i,j,7,ihr0,JMON)
!      endif
#endif

       !solz is read inside hycom.f and forfun.f
!      do ihr=1,12
!       solz2(ihr)=solz_all(i,j,ihr,l0)*w0
!    .            +solz_all(i,j,ihr,l1)*w1
!    .            +solz_all(i,j,ihr,l2)*w2
!    .            +solz_all(i,j,ihr,l3)*w3

!       sunz2(ihr)=acos(solz2(ihr))*rad   !in degs
!      enddo

!        solz=solz2(ihr0)     !because of bi-hourly values
!        sunz=sunz2(ihr0)     !in degs
         solz=ocnatm%cosz1(i,j) !osolz(i,j)      !read instead from modelE
         sunz=acos(solz)*rad  !in degs


!      wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
!    .     +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3

       wind=ocnatm%wsavg(i,j) !owind(i,j)

!      atmFe_ij=atmFe_all(i,j,l0)*w0 + atmFe_all(i,j,l1)*w1
!    .         +atmFe_all(i,j,l2)*w2 + atmFe_all(i,j,l3)*w3

      if (.not.allocated(daily_atmFe)) then
         allocate(daily_atmFe(ogrid%i_strt_halo:ogrid%i_stop_halo,
     &    ogrid%j_strt_halo:ogrid%j_stop_halo))
      endif

      call read_stream(ogrid,stream_atmFe,year,dayOfYear,
     & daily_atmFe)

       !atmospheric deposition iron
      atmFe_ij=daily_atmFe(i,j)

       if (vrbos) then
         write(*,'(/,a,3i5,4e12.4)')'obio_model, forcing: ',
     .   nstep,i,j,solz,sunz,wind,atmFe_ij
       endif

#ifdef OBIO_RUNOFF
       river_runoff=oFLOWO(i,j)
       rnitrconc = rnitrconc_loc(i,j)
       rsiliconc = rsiliconc_loc(i,j)
       rironconc = rironconc_loc(i,j)
       rpocconc  = rpocconc_loc(i,j)
       rdocconc  = rdocconc_loc(i,j)
       rdicconc  = rdicconc_loc(i,j)
#ifdef TRACERS_Alkalinity
       ralkconc  = ralkconc_loc(i,j)
#endif
#endif

#ifdef Relax2SurfN
       if (vrbos) then
         write(*,'(/,a,3i5,e12.4)')'obio_model, surfN: ',
     .   nstep,i,j,surfN(i,j)
       endif
#endif

       idx_co2=ocnatm%gasex_index%getindex(ocnatm%n_co2n)
       if (idx_co2>0) co2flux=ocnatm%trgasex(idx_co2, i, j)

!@PLtest
!@PL      if (vrbos) then
!@PL          write(*,*) 'idx_co2,co2flux,n_co2n = ',idx_co2,co2flux,
!@PL     .    ocnatm%n_co2n       
!@PL      endif 

!@PL

       !------------------------------------------------------------
       !at the beginning of each day only
          call obio_daysetrad(vrbos,i,j,kdm)
          call obio_daysetbio(vrbos,i,j,kdm,nstep,kmax)

             !fill in the 3d arrays to keep in memory for rest of the day
             !when daysetbio is not called again
             tzoo2d(i,j)=tzoo
             do  k=1,kdm
               tfac3d(i,j,k)=tfac(k)
             do nt=1,nchl
               rmuplsr3d(i,j,k,nt)=rmuplsr(k,nt)
               rikd3d(i,j,k,nt)=rikd(k,nt)
             enddo
               wshc3d(i,j,k)=wshc(k)
               Fescav3d(i,j,k)=Fescav(k)
             do nt=1,nlt
               acdom3d(i,j,k,nt)=acdom(k,nt)
             enddo
             enddo

cdiag    if (vrbos) then
cdiag     write (*,103) nstep,i,j,
cdiag.' aftrsetrad   dpth     dp       nitr    ',
cdiag.'   ammo     sili     iron',
cdiag.    (k,p1d(k+1),dp1d(k),obio_P(k,1),obio_P(k,2),
cdiag.                        obio_P(k,3),obio_P(k,4),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' aftrsetrad   dpth     dp       diat       chlo    ',
cdiag.' cyan       cocc     herb',
cdiag.    (k,p1d(k+1),dp1d(k),obio_P(k,5),obio_P(k,6),obio_P(k,7),
cdiag.                        obio_P(k,8),obio_P(k,9),k=1,kdm)
cdiag    endif
 103     format(i9,2i5,a,a/(25x,i3,6(1x,es9.2)))
 104     format(i9,2i5,a,a/(25x,i3,7(1x,es9.2)))


       do ichan = 1,nlt
         Ed(ichan) = 0.0
         Es(ichan) = 0.0
       enddo
       rmud = 0.0
       iyear=2001


         !compute rod and ros only here. not ocean albedo.
         !ocean albedo is computed in ALBEDO.f
         !have to have hygr =  .true. 
       call ocalbedo(wind,solz,dummy,dummy,dummy1,
     .                      rod,ros,.true.,i,j)


cdiag    if (vrbos)
cdiag.     write(*,105)nstep,
cdiag.  ' channel, surf refl dir, surf refl diff',
cdiag.                        (k,rod(k),ros(k),k=1,nlt)
 105  format(i9,a/(18x,i3,2(1x,es9.2)))

cdiag    do k=1,nlt
cdiag    write(*,'(a,4i5,2e12.4)')'obio_model,dir-diff rad:',
cdiag.    nstep,i,j,k,rod(k),ros(k)
cdiag    enddo

         !only call obio_sfcirr for points in light
       tot = 0.0
       if (.not.allocated(eda_frac)) then
         allocate(obio_lambdas(nlt),eda_frac(nlt), esa_frac(nlt))
         call openunit('eda_esa_ratios',iu_bio,.false.,.false.)
         do ichan=1,nlt
           read(iu_bio,'(3f13.8)')obio_lambdas(ichan),eda_frac(ichan),
     &                            esa_frac(ichan)
          enddo
          close(iu_bio)
        endif


       if (allocated(ocnatm%dirvis)) then
         do ichan = 1,nlt
           if (ichan .le. 18) then     !visible + uv
             Ed(ichan) = ocnatm%dirvis(i,j) * eda_frac(ichan)
             Es(ichan) = ocnatm%difvis(i,j) * esa_frac(ichan)
           else               !nir
             Ed(ichan) = ocnatm%dirnir(i,j) * eda_frac(ichan)
             Es(ichan) = ocnatm%difnir(i,j) * esa_frac(ichan)
           endif
           tot = tot + Ed(ichan)+Es(ichan)
          enddo
       endif

#ifdef STANDALONE_OCEAN
         do ichan = 1,nlt
          Ed(ichan) = Eda2(ichan,ihr0)
          Es(ichan) = Esa2(ichan,ihr0)
          tot = tot + Eda2(ichan,ihr0)+Esa2(ichan,ihr0)
         enddo

#endif

!     if (vrbos) then          
!        do ichan=1,nlt
!         write(*,'(a,5i5,3e12.4)')'obio_model, Eds:',
!    .      nstep,i,j,ichan,ihr0,Ed(ichan),Es(ichan),tot
!        enddo
!     endif

#ifdef OBIO_ON_GISSocean
       !integrate over all ichan
         do ichan = 1,nlt
           OIJ(I,J,IJ_ed) = OIJ(I,J,IJ_ed) + Ed(ichan) ! direct sunlight
           OIJ(I,J,IJ_es) = OIJ(I,J,IJ_es) + Es(ichan) ! diffuse sunlight
         enddo  !ichan
#endif


         noon=.false.
         if (hour_of_day.eq.12)then
          if (i.eq.itest.and.j.eq.jtest)noon=.true.
         endif
         if (tot .ge. 0.1) call obio_sfcirr(noon,rod,ros,vrbos)

      if (vrbos) then
        write(*,'(a,3i9)')
     .       'obio_model: counter days,   i,j,nstep=',i,j,nstep
        write(*,'(3(a,i9))')'hour of day=',hour_of_day,
     .                       ', day of month=',day_of_month,
     .                       ', ihr0=',ihr0

cdiag   write(*,'(a)')
cdiag.'    k     dp          u            v         temp         saln'
cdiag   do k=1,kdm
cdiag   write(*,'(i5,5e12.4)')
cdiag   write(*,*)
cdiag.  k,dp1d(k),u(i,j,k+mm),v(i,j,k+mm),
cdiag.    temp(i,j,k+mm),saln(i,j,k+mm)
cdiag   enddo

!       write(*,'(a)')
!    .'    k     P(1)      P(2)         P(3)       P(4)         P(5) '
!       do k=1,kdm
!       write(*,'(i5,7e12.4)')
!       write(*,*)
!    .   k,obio_P(k,1),obio_P(k,2),obio_P(k,3),obio_P(k,4),
!    .   obio_P(k,5),obio_P(k,6),obio_P(k,7)
!       enddo

!       write(*,'(a)')
!    .'    k     P(8)      P(9)         P(11)      P(12)        P(13)'
!       do k=1,kdm
!       write(*,'(i5,7e12.4)')
!       write(*,*)
!    .   k,obio_P(k,8),obio_P(k,9),det(k,1),det(k,2),
!    .   det(k,3),car(k,1),car(k,2)
!       enddo

        write(*,'(2a)')
cdiag.'    Ed          Es          solz         sunz',
cdiag.'       atmFe       wind'
!       write(*,'(6e12.4)')
cdiag   write(*,*)
cdiag.   Ed(ichan),Es(ichan),solz,sunz,atmFe_ij,wind

cdiag.     write(*,106)nstep,
cdiag.    '      channel, dir dwn irr, diff dwn irr,  tot',
cdiag.                  (ichan,Ed(ichan),Es(ichan),
cdiag.                  tot,ichan=1,nlt)
 106  format(i9,a/(18x,i3,3(1x,es9.2)))

      endif  !vrbos

         !this part decomposes light into gmao bands-- dont need it
         !         m = indext2   !use the indicator of "past"
         !         call gmaoirr(ihr,nwater,Edg,Esg)
         !         call oasimbands(idmo,ihr,nwater,Edg,Esg)

         do k=1,kdm
          tirrq(k) = 0.0
         enddo

         delta_temp1d(:) = 0.d0

         if (tot .ge. 0.1) then
           
        call obio_edeu(kmax,vrbos,i,j,im,jm,kdm,
     &                   nstep,obio_lambdas)
        
#ifdef OBIO_ON_GISSocean
#ifdef KPAR_2_OCEAN
       call obio_kpar(kmax,vrbos,i,j,kdm,nstep,dtsrc,dxypo(j),
     &               mo(i,j,:),g0m(i,j,:),s0m(i,j,:),grav,oAPRESS(i,j),
     &               fsr,lsrpd)
       do k=1,kdm
       OIJL(I,J,k,IJL_kpar) = OIJL(I,J,k,IJL_kpar) + Kpar(k) !  kpar
       OIJL(I,J,k,IJL_kpar_em2d) = OIJL(I,J,k,IJL_kpar_em2d) 
     .                           + Kpar_em2d(k) !  kpar in quanta/m2/s
       OIJL(I,J,k,IJL_dtemp) = OIJL(I,J,k,IJL_dtemp) + delta_temp1d(k) !  change in T due to kpar
       enddo
#endif
#endif
         endif  !tot>=0.1


!        if (vrbos) then
!          write(*,107)nstep,
!    .       '           k   avgq    tirrq',
!    .                 (k,avgq1d(k),tirrq(k),k=1,kdm)
!107  format(i9,a/(18x,i3,2(1x,es9.2)))

!        do k=1,kdm
!        write(*,'(a,4i5,2e12.5)')'obio_model,k,avgq,tirrq:',
!    .       nstep,i,j,k,avgq1d(k),tirrq(k)
!        enddo
!        endif

       !------------------------------------------------------------
       !compute tendency terms on the m level
cdiag     if (vrbos) then
cdiag     write (*,103) nstep,i,j,
cdiag.' bfreptend  dpth       dp       nitr    ',
cdiag.'     ammo     sili     iron',
cdiag.    (k,p(i,j,k+1)/onem,dp1d(k),obio_P(k,1),obio_P(k,2),
cdiag.                       obio_P(k,3),obio_P(k,4),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' bfreptend  dpth     diat       chlo    ',
cdiag.' cyan       cocc     herb',
cdiag.    (k,p(i,j,k+1)/onem,obio_P(k,5),obio_P(k,6),obio_P(k,7),
cdiag.                       obio_P(k,8),obio_P(k,9),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' bfreptend  dpth     nitr_det  sili_det iron_det   doc   ',
cdiag.'    dic ',
cdiag.    (k,p(i,j,k+1)/onem,det(k,1),det(k,2),det(k,3),
cdiag.                       car(k,1),car(k,2),k=1,kdm)
cdiag     endif
       !------------------------------------------------------------

cdiag  if (vrbos)write(*,*)'bfre obio_ptend: ',
cdiag.     nstep,(k,tirrq(k),k=1,kmax)


       n_co2n=ocnatm%n_co2n
#ifdef OBIO_ON_GISSocean
        mmo=mo(i,j,1)
        ddxypo=dxypo(j)
        SDIC=trmo(i,j,1,n_abioDIC)
        oij_pH=oij(i,j,ij_ph)
!       oij_co3=oij(i,j,ij_co3)
#else
        oij_pH=pHav_loc(i,j)
#endif
       
!       call check_sumcarbon(nstep,i,j,kmax,'check1')

        call obio_ptend(vrbos,kmax,i,j,kdm,nstep,n_co2n,
     &                 DTS,mmo,ddxypo,n_abioDIC,num_tracers,SDIC,oij_pH
#ifdef TRACERS_Alkalinity
     &                ,oij_co3
#endif
     &                 )
!       call check_sumcarbon(nstep,i,j,kmax,'check2')
#ifdef OBIO_ON_GISSocean
#ifdef TOPAZ_params
#ifdef TRACERS_Alkalinity
       oij_co3 = oij_co3 + co3_conc
        oij(i,j,ij_co3)=oij_co3
#endif
#endif
        oij_pH=oij_pH+pHsfc
        oij(i,j,ij_ph)=oij_pH
#endif
       !------------------------------------------------------------
cdiag  if (vrbos)then
cdiag   write(*,108)nstep,' aftrptend dpth      dp     P_tend(1:9)',
cdiag.    (k,p(i,j,k+1)/onem,dp1d(k),P_tend(k,1),P_tend(k,2),
cdiag.     P_tend(k,3),P_tend(k,4),P_tend(k,5),P_tend(k,6),
cdiag.     P_tend(k,7),P_tend(k,8),P_tend(k,9),k=1,kdm)
cdiag   write(*,109)nstep,
cdiag.    ' aftrptend dpth      dp     D_tend(1:3) and C_tend(1:2)',
cdiag.    (k,p(i,j,k+1)/onem,dp1d(k),D_tend(k,1),D_tend(k,2),
cdiag.     D_tend(k,3),C_tend(k,1),C_tend(k,2),k=1,kdm)
cdiag  endif
 108  format(i6,a/(12x,i3,11(1x,es9.2)))
 109  format(i6,a/(12x,i3,7(1x,es9.2)))


       !------------------------------------------------------------

#ifdef OBIO_ON_GISSocean
       !update biology to new time level
       !also do phyto sinking and detrital settling here
       !MUST CALL sinksettl BEFORE update
       ddxypo=dxypo(j)
       call obio_sinksettl(vrbos,kmax,errcon,i,j,
     &                kdm,nstep,dtsrc,ddxypo)
!      call check_sumcarbon(nstep,i,j,kmax,'check3')
       call obio_update(vrbos,kmax,i,j,nstep)

!@PLdbg if errchk ~=0, write error status and stop
#ifdef TRACERS_bio_O2
       if (errchk1.eq.1.or.errchk2.eq.1) then
        write(6,'(a,2i2)')'O2 has NaNs, error status=',
     .      errchk1,errchk2
        call stop_model('O2 contains NaNs, aborting',255)

        endif
#endif

#else
       !update biology from m to n level
       !also do phyto sinking and detrital settling here
       !MUST CALL sinksettl AFTER update
       call obio_update(vrbos,kmax,i,j)
       ddxypo=scp2(i,j)
       call obio_sinksettl(vrbos,kmax,errcon,i,j,
     &                kdm,nstep,dtsrc,ddxypo) 
       if (errcon) then
          write (*,'(a,3i5)') 'error update at i,j =',i,j,kmax
          do k=1,kdm
          write (*,'(i5,3e12.4)')
     .           k,dp1d(k),p1d(k),wsdet(k,1)
          enddo
          stop
       endif
#endif

cdiag     if (vrbos) then
cdiag     write (*,*)'     '
cdiag     write (*,103) nstep,i,j,
cdiag.' aftrupdate dpth     dp        nitr      ammo      sili',
cdiag.'      iron',
cdiag.  (k,p(i,j,k+1)/onem,dp1d(k),obio_P(k,1),obio_P(k,2),
cdiag.                     obio_P(k,3),obio_P(k,4),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' aftrupdate dpth     diat      chlo        cyan',
cdiag.'      cocc      herb',
cdiag.  (k,p(i,j,k+1)/onem,obio_P(k,5),obio_P(k,6),obio_P(k,7),
cdiag.                       obio_P(k,8),obio_P(k,9),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' aftrupdate dpth     nitr_det  sili_det    iron_det',
cdiag.'      doc       dic ',
cdiag.    (k,p(i,j,k+1)/onem,det(k,1),det(k,2),det(k,3),
cdiag.                       car(k,1),car(k,2),k=1,kdm)
cdiag     endif


      !------------------------------------------------------------
      !integrate rhs vertically 
      do ll= 1, 17
      do nt= 1, ntrac
      rhs_obio(i,j,nt,ll) = 0.d0
      do k = 1, kdm
          rhs_obio(i,j,nt,ll) = rhs_obio(i,j,nt,ll) +
     .                 rhs(k,nt,ll)*dp1d(k)    
      enddo  !k

c     now all rates are in /s    July 2016
c     if (vrbos) then
c     write(*,'(a,5i5,1x,e20.13)')'rhs_obio (mass,trac/m2/s):',
c    .   nstep,i,j,nt,ll,rhs_obio(i,j,nt,ll)
c     endif
      enddo  !ntrac

      !convert all to mili-mol,C/m2
      rhs_obio(i,j,5,ll) = rhs_obio(i,j,5,ll) * mgchltouMC
      rhs_obio(i,j,6,ll) = rhs_obio(i,j,6,ll) * mgchltouMC
      rhs_obio(i,j,7,ll) = rhs_obio(i,j,7,ll) * mgchltouMC
      rhs_obio(i,j,8,ll) = rhs_obio(i,j,8,ll) * mgchltouMC
      rhs_obio(i,j,9,ll) = rhs_obio(i,j,9,ll) * mgchltouMC
      rhs_obio(i,j,10,ll) = rhs_obio(i,j,10,ll) / uMtomgm3

      chng_by(i,j,1) = (rhs_obio(i,j,5,9)+rhs_obio(i,j,6,9)
     .                  +rhs_obio(i,j,7,9)+rhs_obio(i,j,8,9))
     .               +  rhs_obio(i,j,9,9)
      chng_by(i,j,2) = (rhs_obio(i,j,5,5)+rhs_obio(i,j,6,6)
     .               +  rhs_obio(i,j,7,7)+rhs_obio(i,j,8,8))
     .               +  rhs_obio(i,j,10,5)
      chng_by(i,j,3) = (rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .               +  rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13))
     .               + rhs_obio(i,j,14,6)
      chng_by(i,j,4) = (rhs_obio(i,j,5,14)+rhs_obio(i,j,6,14)
     .               +  rhs_obio(i,j,7,14)+rhs_obio(i,j,8,14))
     .               +  rhs_obio(i,j,13,5)
      chng_by(i,j,5) = (rhs_obio(i,j,5,15)+rhs_obio(i,j,6,15)
     .               +  rhs_obio(i,j,7,15)+rhs_obio(i,j,8,15))
     .               +  rhs_obio(i,j,14,5)
      chng_by(i,j,6) =  rhs_obio(i,j,9,10) 
     .                            + rhs_obio(i,j,13,9)
      chng_by(i,j,7) =rhs_obio(i,j,9,11)
     .                            + rhs_obio(i,j,10,7) 
      chng_by(i,j,8) = rhs_obio(i,j,9,12)
     .                            + rhs_obio(i,j,10,8)
      chng_by(i,j,9) = rhs_obio(i,j,9,14)
     .                            + rhs_obio(i,j,13,12) 
      chng_by(i,j,10)= rhs_obio(i,j,9,15)
     .                            + rhs_obio(i,j,14,15)
      chng_by(i,j,11)= rhs_obio(i,j,10,14)
     .                            + rhs_obio(i,j,13,13)
      chng_by(i,j,12)= rhs_obio(i,j,10,9)
     .                            + rhs_obio(i,j,13,10)
      chng_by(i,j,13)= rhs_obio(i,j,10,10)
     .                            + rhs_obio(i,j,14,10)
      chng_by(i,j,14)= rhs_obio(i,j,13,14)
     .                            + rhs_obio(i,j,14,14)
 
      enddo  !ll

c     call obio_chkbalances(vrbos,nstep,i,j)

#ifdef OBIO_ON_GISSocean
      do nt=1,ntrac
      do ll=1,17
      OIJ(I,J,IJ_rhs(nt,ll)) = OIJ(I,J,IJ_rhs(nt,ll))
     .                                    + rhs_obio(i,j,nt,ll)  ! all terms in rhs
      enddo
      enddo
#endif

#ifdef obio_rhsdiags
      call save_rhs3_diags(nstep,I,J,kdm)
#endif

!     if (vrbos) then
!      print*, 'OBIO TENDENCIES, 1-17, 1,7'
!      do k=1,1
!       write(*,'(17(e9.2,1x))')((rhs(k,nt,ll),ll=1,17),nt=1,7)
!        print*, 'OBIO TENDENCIES, 1-17, 8,14'
!       write(*,'(17(e9.2,1x))')((rhs(k,nt,ll),ll=1,17),nt=8,ntrac-1)
!      enddo
!     endif


       !update 3d tracer array
       do k=1,kmax
        do nt=1,ntyp
         tracer(i,j,k,nt)=obio_P(k,nt)
        enddo
        do nt=1,ndet
         tracer(i,j,k,ntyp+nt)=det(k,nt)
        enddo
        do nt=1,ncar
         tracer(i,j,k,ntyp+ndet+nt)=car(k,nt)
        enddo
#ifdef TRACERS_Alkalinity
         nt=1
         tracer(i,j,k,ntyp+ndet+ncar+nt)=alk1d(k)
#ifdef TOPAZ_params
         nt=2
         tracer(i,j,k,ntyp+ndet+ncar+nt)=ca_det_calc1d(k)
#endif
#endif
#ifdef TRACERS_Ocean_O2 
#ifdef TRACERS_bio_O2
         tracer(i,j,k,ntyp+ndet+ncar+nalk+no2)=o21d(k)
#endif
#ifdef TRACERS_abio_O2
         tracer(i,j,k,ntyp+ndet+ncar+nalk+no2+nabo2)=abo21d(k)
#endif
#endif
#ifdef TRACERS_degC
        tracer(i,j,k,ndimndegC) = ndegC1d(k)
#endif
        !update avgq and gcmax arrays
        avgq(i,j,k)=avgq1d(k)
        OIJL(I,J,k,IJL_avgq)= OIJL(I,J,k,IJL_avgq) + avgq1d(k)
        gcmax(i,j,k)=gcmax1d(k)
        tirrq3d(i,j,k)=tirrq(k)
        alk(i,j,k)=alk1d(k)
       enddo !k

#ifdef OBIO_ON_GISSocean
      !update trmo etc arrays
       do k=1,kmax
       do nt=1,ntrac
          dtr = tracer(i,j,k,nt) * trmo_unit_factor(k,nt) 
     .        - trmo(i,j,k,nt)

        if (dtr.lt.0) then
          ftr = -dtr/trmo(i,j,k,nt)
          txmo(i,j,k,nt)=txmo(i,j,k,nt)*(1.-ftr)
          tymo(i,j,k,nt)=tymo(i,j,k,nt)*(1.-ftr)
          tzmo(i,j,k,nt)=tzmo(i,j,k,nt)*(1.-ftr)       
          if (use_qus==1) then
             txxmo(i,j,k,nt)=txxmo(i,j,k,nt)*(1.-ftr)       
             txymo(i,j,k,nt)=txymo(i,j,k,nt)*(1.-ftr)       
             tzxmo(i,j,k,nt)=tzxmo(i,j,k,nt)*(1.-ftr)       
             tyymo(i,j,k,nt)=tyymo(i,j,k,nt)*(1.-ftr)       
             tyzmo(i,j,k,nt)=tyzmo(i,j,k,nt)*(1.-ftr)       
             tzzmo(i,j,k,nt)=tzzmo(i,j,k,nt)*(1.-ftr)       
          endif
        endif

        trmo(i,j,k,nt)=trmo(i,j,k,nt)+dtr

       enddo
       enddo
#else
       tracer_h(i,j,:,:)=tracer(i,j,:,:)*trmo_unit_factor
#endif

       !compute total chlorophyl at surface layer
       tot_chlo(i,j)=0.
       do nt=1,nchl
          tot_chlo(i,j)=tot_chlo(i,j)+obio_P(1,nnut+nt)
       enddo
       ocnatm%chl(i,j) = tot_chlo(i,j)
       if (vrbos) then
          !!!write(*,'(/,a,3i5,e12.4)')
          write(*,*)
     .       'obio_model, tot_chlo= ',nstep,i,j,tot_chlo(i,j)
       endif

       !compute total primary production per day
       !sum pp over species and depth, NOT over day
       pp2tot_day(i,j)=0.
       pp2diat_day(i,j)=0.
       pp2chlo_day(i,j)=0.
       pp2cyan_day(i,j)=0.
       pp2cocc_day(i,j)=0.
!      if (p1d(kdm+1).gt.200.) then    !if total depth > 200m
          do nt=1,nchl
          do k=1,kdm
          if (nt.eq.1)pp2diat_day(i,j)
     .       =pp2diat_day(i,j)+pp2_1d(k,1)   !mg,C/m2/day     !July 2016: pp2_1d is already mgC/m2/day in ptend.f
          if (nt.eq.2)pp2chlo_day(i,j)
     .       =pp2chlo_day(i,j)+pp2_1d(k,2)  !mg,C/m2/day
          if (nt.eq.3)pp2cyan_day(i,j)
     .       =pp2cyan_day(i,j)+pp2_1d(k,3)   !mg,C/m2/day
          if (nt.eq.4)pp2cocc_day(i,j)
     .       =pp2cocc_day(i,j)+pp2_1d(k,4)   !mg,C/m2/day

              pp2tot_day(i,j)=pp2tot_day(i,j)+pp2_1d(k,nt)     !mg,C/m2/day
!    &                                       * HOURS_PER_DAY   !->mg,C/m2/day
          enddo
          enddo
!       endif
       if (vrbos) then
       write(*,'(a,i8,2i5,7e12.4)')'obio_model, pp:',
     .    nstep,i,j,tot_chlo(i,j),tot,pp2diat_day(i,j),
     .    pp2chlo_day(i,j),pp2cyan_day(i,j),pp2cocc_day(i,j),
     .    pp2tot_day(i,j)
       endif

#ifndef STANDALONE_OCEAN
       ocnatm%gtracer(ocnatm%n_co2n, i,j)=pCO2_ij
#endif

#ifdef OBIO_ON_GISSocean
!diagnostics
       if (tot.gt.0) then
       OIJ(I,J,IJ_dayl) = OIJ(I,J,IJ_dayl) + 1.d0*dtsrc/3600.d0 !number of hrs daylight   
       endif
       OIJ(I,J,IJ_solz) = OIJ(I,J,IJ_solz) + solz  !cos zenith angle
       OIJ(I,J,IJ_sunz) = OIJ(I,J,IJ_sunz) + sunz  !solar zenith angle in degrees  

       OIJ(I,J,IJ_nitr) = OIJ(I,J,IJ_nitr) + tracer(i,j,1,1) ! surf ocean nitrate
       OIJ(I,J,IJ_amm)  = OIJ(I,J,IJ_amm)  + tracer(i,j,1,2) ! surf ocean ammonium
       OIJ(I,J,IJ_sil)  = OIJ(I,J,IJ_sil)  + tracer(i,j,1,3) ! surf ocean silicate
       OIJ(I,J,IJ_iron) = OIJ(I,J,IJ_iron) + tracer(i,j,1,4) ! surf ocean iron    

       OIJ(I,J,IJ_diat) = OIJ(I,J,IJ_diat) + tracer(i,j,1,5) ! surf ocean diatoms
       OIJ(I,J,IJ_chlo) = OIJ(I,J,IJ_chlo) + tracer(i,j,1,6) ! surf ocean chlorophytes
       OIJ(I,J,IJ_cyan) = OIJ(I,J,IJ_cyan) + tracer(i,j,1,7) ! surf ocean cyanobacteria
       OIJ(I,J,IJ_cocc) = OIJ(I,J,IJ_cocc) + tracer(i,j,1,8) ! surf ocean coccolithophores
       OIJ(I,J,IJ_herb) = OIJ(I,J,IJ_herb) + tracer(i,j,1,9) ! surf ocean herbivores
       OIJ(I,J,IJ_pp) = OIJ(I,J,IJ_pp) + pp2tot_day(i,j)     ! depth integrated pp

       OIJ(I,J,IJ_pp1) = OIJ(I,J,IJ_pp1) + pp2diat_day(i,j)  ! ocean pp from diatoms
       OIJ(I,J,IJ_pp2) = OIJ(I,J,IJ_pp2) + pp2chlo_day(i,j)  ! ocean pp from chlorophytes
       OIJ(I,J,IJ_pp3) = OIJ(I,J,IJ_pp3) + pp2cyan_day(i,j)  ! ocean pp from cyanobacteria
       OIJ(I,J,IJ_pp4) = OIJ(I,J,IJ_pp4) + pp2cocc_day(i,j)  ! ocean pp from coccolithophores

       OIJ(I,J,IJ_doc) = OIJ(I,J,IJ_doc) + tracer(i,j,1,13)  ! surf ocean doc
       OIJ(I,J,IJ_dic) = OIJ(I,J,IJ_dic) + tracer(i,j,1,14)  ! surf ocean dic
       OIJ(I,J,IJ_pCO2)= OIJ(I,J,IJ_pCO2)+ pCO2_ij !  *(1.-oice(i,j)) surf ocean pco2 @PL removed sea ice scaling

       OIJ(I,J,IJ_cexp) = OIJ(I,J,IJ_cexp) + cexp             ! export production
       OIJ(I,J,IJ_ndet) = OIJ(I,J,IJ_ndet) + tracer(i,j,kzc,10) ! ndet at zc
#ifdef TRACERS_degC
     .                 + tracer(i,j,kzc,ndimndegC)
       OIJ(I,J,IJ_cexpdeg) = OIJ(I,J,IJ_cexpdeg) + cexpdeg  ! export production of degradable carbon
       OIJ(I,J,IJ_degC) = OIJ(I,J,IJ_degC) + tracer(i,j,kzc,10) ! degradable C at zc
#endif
       OIJ(I,J,IJ_setl)= OIJ(I,J,IJ_setl)+ wsdet(kzc,1)          ! settl vel. n/cdet at kzc
       OIJ(I,J,IJ_sink)= OIJ(I,J,IJ_sink)+ obio_ws(kzc,1)        ! sink. vel. phytoplankton
       if (kzc<=kmax) then
         OIJ(I,J,IJ_xchl) = OIJ(I,J,IJ_xchl)
     .       + tracer(i,j,kzc,5)*trmo_unit_factor(kzc,5)*obio_ws(kzc,1)
     .       + tracer(i,j,kzc,6)*trmo_unit_factor(kzc,6)*obio_ws(kzc,1)
     .       + tracer(i,j,kzc,7)*trmo_unit_factor(kzc,7)*obio_ws(kzc,1)
     .       + tracer(i,j,kzc,8)*trmo_unit_factor(kzc,8)*obio_ws(kzc,1) !total phyto cexp at kzc
       else
         oij(i,j,ij_xchl)=0
       endif

       !limitation diags surface only (for now)
       k = 1
       do nt=1,nchl
       do ilim=1,5
       OIJ(I,J,IJ_lim(nt,ilim)) = OIJ(I,J,IJ_lim(nt,ilim))
     .                          + flimit(k,nt,ilim)
       enddo
       enddo

       !3d pp diags
       do nt=1,nchl
       do k=1,kdm
          OIJL(I,J,k,IJL_pp(nt)) = OIJL(I,J,k,IJL_pp(nt))
     .                           + pp2_1d(k,nt)     !mg,C/m2/day 
       enddo
       enddo

       !3d limitation diags
       do nt=1,nchl
       do k=1,kdm
        OIJL(I,J,k,IJL_lim1(nt))=OIJL(I,J,k,IJL_lim1(nt))+flimit(k,nt,1)
        OIJL(I,J,k,IJL_lim2(nt))=OIJL(I,J,k,IJL_lim2(nt))+flimit(k,nt,2)
        OIJL(I,J,k,IJL_lim3(nt))=OIJL(I,J,k,IJL_lim3(nt))+flimit(k,nt,3)
        OIJL(I,J,k,IJL_lim4(nt))=OIJL(I,J,k,IJL_lim4(nt))+flimit(k,nt,4)
        OIJL(I,J,k,IJL_lim5(nt))=OIJL(I,J,k,IJL_lim5(nt))+flimit(k,nt,5)
       enddo
       enddo


       OIJ(I,J,IJ_flux) = OIJ(I,J,IJ_flux) + co2flux!*sday*365.24      !air-sea CO2 flux(if on ocean grid, this is gr,CO2/m2/yr, if coupled it is in molCO2/m2/yr)
                                                                      !@PL multiply by sday*365.24 so units are molCO2/m2/yr, else it is molCO2/m2/s
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
       OIJ(I,J,IJ_oflx) = OIJ(I,J,IJ_oflx) + o2flux       !air-sea O2 flux (molO2/m2/yr)

!@PL rhs terms withou obiorhs defined, in mmol/kg/day

      do k=1,kdm
      OIJL(I,J,k,IJL_oprodam)=OIJL(I,J,k,IJL_oprodam)
     .                          + rhs(k,ndimo2,6)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar+nalk+no2)
      OIJL(I,J,k,IJL_oprodnit)=OIJL(I,J,k,IJL_oprodnit)+
     .                         rhs(k,ndimo2,7)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar+nalk+no2)
      OIJL(I,J,k,IJL_oaresp)=OIJL(I,J,k,IJL_oaresp) +
     .                          rhs(k,ndimo2,5)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar+nalk+no2)
      OIJL(I,J,k,IJL_ohresp)=OIJL(I,J,k,IJL_ohresp) +
     .                             rhs(k,ndimo2,15)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar+nalk+no2)
      OIJL(I,J,k,IJL_odet)=OIJL(I,J,k,IJL_odet) +
     .                               rhs(k,ndimo2,10)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar+nalk+no2)
      OIJL(I,J,k,IJL_odoc)=OIJL(I,J,k,IJL_odoc) +
     .                          rhs(k,ndimo2,14)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar+nalk+no2)

      OIJL(I,J,k,IJL_cprod)=OIJL(I,J,k,IJL_cprod)
     .                          + rhs(k,ndimc,6)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar)
      OIJL(I,J,k,IJL_caresp)=OIJL(I,J,k,IJL_caresp) +
     .                          rhs(k,ndimc,5)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar)
      OIJL(I,J,k,IJL_chresp)=OIJL(I,J,k,IJL_chresp) +
     .                             rhs(k,ndimc,15)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar)
      OIJL(I,J,k,IJL_cdet)=OIJL(I,J,k,IJL_cdet) +
     .                               rhs(k,ndimc,10)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar)
      OIJL(I,J,k,IJL_cdoc)=OIJL(I,J,k,IJL_cdoc) +
     .                          rhs(k,ndimc,14)*sday
     .           * trmo_unit_factor(k,nnut+nchl+nzoo+ndet+ncar)


!       if (vrbos) then
!          write(6,'(a,3i7,2e12.4)')'obio_o2(bio):',
!     .    nstep,i,j,rhs(k,ndimo2,6),rhs(k,ndimo2,7)
!        endif
!
      enddo

#endif
#ifdef TRACERS_abio_O2
       OIJ(I,J,IJ_aboflx) = OIJ(I,J,IJ_aboflx) + abo2flux       !air-sea abiotic O2 flux (molO2/m2/yr)
#endif
#endif

       !phyt. sinking speed in m/s
       do nt=1,nchl
       do k=1,kmax
       OIJL(I,J,k,IJL_wss(nt)) = OIJL(I,J,k,IJL_wss(nt)) + obio_ws(k,nt)
       enddo
       enddo

       !detr. sinking speed in m/s
       do nt=1,ndet
       do k=1,kmax
       OIJL(I,J,k,IJL_wsdet(nt)) = OIJL(I,J,k,IJL_wsdet(nt))
     .                           + wsdet(k,nt)
       enddo
       enddo

       if (tracers_alkalinity) then
         OIJ(I,J,IJ_alk) = OIJ(I,J,IJ_alk) 
     .          + tracer(i,j,1,nnut+nchl+nzoo+ndet+ncar+nalk)    ! surf ocean alkalinity
         OIJ(I,J,IJ_fca) = OIJ(I,J,IJ_fca) + caexp               ! carbonate export
       else
         OIJ(I,J,IJ_alk) = OIJ(I,J,IJ_alk) + alk(i,j,1)          ! surf ocean alkalinity
       endif
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
         OIJ(I,J,IJ_o2) = OIJ(I,J,IJ_o2)
     .          + tracer(i,j,1,nnut+nchl+nzoo+ndet+ncar+nalk+no2) ! surf ocean oxygen concentration
         OIJ(I,J,IJ_po2) = OIJ(I,J,IJ_po2) + pO2_ij ! *(1.-oice(i,j))  surf ocean partial O2 pressure 
#endif
#ifdef TRACERS_abio_O2
         OIJ(I,J,IJ_abo2) = OIJ(I,J,IJ_abo2)
     .          + tracer(i,j,1,nnut+nchl+nzoo+ndet+ncar+nalk+no2+nabo2) ! surf ocean oxygen concentration
         OIJ(I,J,IJ_pabo2)=OIJ(I,J,IJ_pabo2) + pabO2_ij ! (1.-oice(i,j))  surf ocean partial O2 pressure @PL removed sea ice scaling
#endif
#endif

#ifdef OBIO_RUNOFF
!      OIJ(I,J,IJ_rnitrmflo) = OIJ(I,J,IJ_rnitrmflo)+ rnitrmflo    ! riverine nitr mass flow from dC (kg/s)
       OIJ(I,J,IJ_rnitrconc) = OIJ(I,J,IJ_rnitrconc)+ rnitrconc    ! riverine nitr conc from dC (kg/kg)
       OIJ(I,J,IJ_rdicconc) = OIJ(I,J,IJ_rdicconc)  + rdicconc     ! riverine dic conc from dC (kg/kg)
       OIJ(I,J,IJ_rdocconc) = OIJ(I,J,IJ_rdocconc)  + rdocconc     ! riverine doc conc from dC (kg/kg)
       OIJ(I,J,IJ_rsiliconc) = OIJ(I,J,IJ_rsiliconc)+ rsiliconc    ! riverine silica conc from dC (kg/kg)
       OIJ(I,J,IJ_rironconc) = OIJ(I,J,IJ_rironconc)+ rironconc    ! riverine iron conc from dC (kg/kg)
       OIJ(I,J,IJ_rpocconc) = OIJ(I,J,IJ_rpocconc)  + rpocconc     ! riverine poc conc from dC (kg/kg)
       OIJ(I,J,IJ_ralkconc) = OIJ(I,J,IJ_ralkconc)  + ralkconc     ! riverine alkalinity conc from A-S (mol/kg)
#endif

#endif  /*OBIO_ON_GISSocean*/

#ifndef OBIO_ON_GISSocean    /* HYCOM ACCUMULATED DIAGNOSTICS */
      ao_co2fluxav_loc(i,j)=ao_co2fluxav_loc(i,j) + co2flux
      pCO2av_loc(i,j)=pCO2av_loc(i,j)+pCO2_ij
      pp2tot_dayav_loc(i,j) = pp2tot_dayav_loc(i,j) 
     .                      + pp2tot_day(i,j)
      pp2diat_dayav_loc(i,j) = pp2diat_dayav_loc(i,j)
     .                      + pp2diat_day(i,j)
      pp2chlo_dayav_loc(i,j) = pp2chlo_dayav_loc(i,j)
     .                      + pp2chlo_day(i,j)
      pp2cyan_dayav_loc(i,j) = pp2cyan_dayav_loc(i,j)
     .                      + pp2cyan_day(i,j)
      pp2cocc_dayav_loc(i,j) = pp2cocc_dayav_loc(i,j)
     .                      + pp2cocc_day(i,j)
      cexpav_loc(i,j)=cexpav_loc(i,j)+cexp
#ifdef TRACERS_Alkalinity
            caexpav_loc(i,j)=caexpav_loc(i,j)+caexp
#endif
#endif

      end do
      end do
      call stop('  obio main loop')

      call stop(' obio_model')
      initialized=.true.
      nstep0=nstep0+1
      return

      end subroutine obio_model



      subroutine check_sumcarbon(nstep,i,j,kmax,chkpt)

      USE obio_com,  only: obio_P,det,car,trmo_unit_factor

      implicit none
      integer, intent(in) :: nstep,i,j,kmax
      character*6,intent(in) :: chkpt
      real*8,save :: last = 0
      real*8, dimension(kmax) :: diat,chlo,cyan,cocc,herb
     .                          ,ndet,doc,dic,sumcarbon

      !all units ->kg,C
      !phyto
      diat = obio_P(:,5) * trmo_unit_factor(:,5)
      chlo = obio_P(:,6) * trmo_unit_factor(:,6)
      cyan = obio_P(:,7) * trmo_unit_factor(:,7)
      cocc = obio_P(:,8) * trmo_unit_factor(:,8)
      !herb
      herb = obio_P(:,9) * trmo_unit_factor(:,9)
      !ndet
      ndet = det(:,1) * trmo_unit_factor(:,10)
      !doc
      doc  = car(:,1) * trmo_unit_factor(:,13)
      !dcc
      dic  = car(:,2) * trmo_unit_factor(:,14)

      sumcarbon  = (diat+chlo+cyan+cocc+herb)
     .           + ndet
     .           + (doc+dic)

      write(*,'(2a,3i5,10e12.4)') 'sumcarbon (in kg) at ',chkpt,
     .     nstep,i,j,
     .     sum(diat),sum(chlo),sum(cyan),sum(cocc),sum(herb),
     .     sum(ndet),sum(doc),sum(dic),
     .     sum(sumcarbon),sum(sumcarbon)-last

      last=sum(sumcarbon)

      end subroutine check_sumcarbon

#ifdef obio_rhsdiags
      subroutine save_rhs3_diags(nstep,I,J,kdm)

      USE obio_dim
      USE obio_diag, only : oijl=>rhs_ijl,ijl_rhs3
      USE obio_com,  only: rhs,trmo_unit_factor

      implicit none
      integer, intent(in) :: nstep,i,j,kdm
      integer :: k,nt,ll

!@PL include all or only a few tracer tendancies 
! depending on if oxygen is included
#ifdef TRACERS_Ocean_O2
      do nt=ndimc,ntrac
       if (nt.eq.ndimc
#ifdef TRACERS_bio_O2
     & .or. nt.eq.ndimo2
#endif
#ifdef TRACERS_abio_O2
     & .or. nt.eq.ndimabo2
#endif
     &           ) then
       do ll=1,17
       do k=1,kdm
       OIJL(I,J,k,IJL_rhs3(nt,ll)) = OIJL(I,J,k,IJL_rhs3(nt,ll))
     .                           + rhs(k,nt,ll)*trmo_unit_factor(k,nt) !units in kg/s
       enddo
       enddo
       endif
      enddo
#else


      do nt=1,ntrac
       do ll=1,17
       do k=1,kdm
       OIJL(I,J,k,IJL_rhs3(nt,ll)) = OIJL(I,J,k,IJL_rhs3(nt,ll))
     .                           + rhs(k,nt,ll)*trmo_unit_factor(k,nt) !units in kg/s
       enddo
       enddo
      enddo
#endif
      end subroutine save_rhs3_diags
#endif

