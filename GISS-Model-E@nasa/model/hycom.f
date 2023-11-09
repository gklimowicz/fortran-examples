#include "hycom_mpi_hacks.h"
#include "rundeck_opts.h"

#if defined(CUBED_SPHERE) || defined(NEW_IO)
#else
#define USE_ATM_GLOBAL_ARRAYS
#endif

      subroutine OCEANS(atmocn,iceocn,dynsice)
c
c --- ------------------------------
c --- MICOM-based hybrid ocean model
c
c ---     v e r s i o n    0.9
c --- ------------------------------
c
c                            Program Notes:
c
c                    Online Tracers in MICOM and HYCOM
c                    ---------------------------------
c
c 1. The number of tracers is set by 'ntrcr' in dimensions.h.
c
c 2. Source/sink functions for each tracer must be supplied by the user.
c In one current implementation, ntrcr=2, and there is a call in hycom.f
c to a routine called 'addtrc' which continuously adds tracer at 2
c "pipeline outfall" points on the sea floor. In another implementation,
c ntrcr=21, and hycom calls a biochemical submodel encompassing 21
c biological and chemical tracers.
c
c 3. Tracer transport takes place once every 'trcfrq' time steps; in
c other words, the transport time step is trcfrq*baclin. The transport
c equation is solved in flux form and is based on a time-integrated form
c of the mass continuity equation in which interlayer mass exchange is
c determined as a residual based on cumulative layer thickness change
c and time-integrated horizontal mass flux divergence.
c
c 4. Tracer transport in  n o t  done in leapfrog mode. In the present
c implementation, 'trcfrq' must be even. Only mass fluxes at even time
c steps are used in building up the mass flux time integrals. This is
c to avoid leapfrog-related inconsistencies between layer thickness
c tendencies and horizontal flux fields which would generate spurious
c diapycnal mass fluxes. The single remaining inconsistency is the one
c caused by time smoothing of the thickness field. It spawns small
c diapyncal fluxes (and hence interlayer tracer leakage) even in the
c absence of physically based vertical mixing processes.
c
c 5. Transport is handled by subroutine trcadv.f which has 3 entries:
c
c    -  tradv0 is called at the beginning of each 'long' time interval.
c       It initializes the mass flux time integrals and stores the
c       initial layer thickness.
c
c    -  tradv1 is called every even-numbered 'regular' model time step.
c       It builds up the mass flux time integrals.
c
c    -  tradv2 is called at the end of each 'long' time interval. It
c       diagnoses the diapycnal mass flux components and carries out the
c       actual 3-D tracer transport, using a combination of FCT (flux
c       corrected transport)in the horizontal and PPM (piecewise parabolic
c       method) in the vertical.
c
c 6. Solving transport equations in flux form makes it impossible to
c rigorously conserve total tracer when layer thickness is allowed to go
c to zero. Conservation errors are currently eliminated by applying a
c global correction factor to each tracer field following each transport
c step. With this corrective device in place, tracer is rigorously
c conserved in the model.
c
c 7. Upon exit from tradv2, tracer fields are valid at the current time
c level 'n' and hence can be paired with the current layer thickness
c dp(:,:,k+nn). This temporal consistency presents a window of
c opportunity, which remains open until the next call to tradv0, for
c performing vertical mixing operations on the tracer fields. The logical
c variable 'dotrcr' is set to 'true' during this time window and can be
c used to direct T/S mixing routines to also mix tracer. For numerical
c consistency, the time step for tracer mixing in those routines should
c be set to trcfrq*baclin. Since convective adjustment and Kraus-Turner
c mixed layer routines do not solve a diffusion equation with an explicit
c time step, mathematical consistency of tracer treatment by convec and
c mxlayr is not obvious. These routines operate on the layer thickness
c field every 'regular' time step. Their diffusive effect on the tracer
c fields is therefore captured mainly by diapycnal fluxes spawned by
c layer thickness trends which result from vertical mixing.
c
c 8. Since layer thickness in HYCOM's z-coordinate subdomain does not
c reflect the action of vertical mixing processes the way it does in
c MICOM, vertical mixing by convec and mxlayr is likely to be
c underestimated in HYCOM. This problem is alleviated by using
c vertical mixing schemes like KPP (with time step trcfrq*baclin).
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT, HALO_UPDATE, NORTH,
     &                         haveLatitude, globalsum
      USE DOMAIN_DECOMP_1D, only: pack_data, unpack_data,
     &     band_pack
      USE HYCOM_ATM
!      USE FLUXES, only : e0,prec,eprec,evapor,flowo,eflowo,dmua,dmva
!     . ,erunosi,runosi,srunosi,runpsi,srunpsi,dmui,dmvi,dmsi,dhsi,dssi
!     . ,gtemp,sss,mlhc,ogeoza,uosurf,vosurf,MELTI,EMELTI,SMELTI
!     . ,gmelt,egmelt,solar,gtempr,erunpsi

      !USE SEAICE_COM, only : rsi,msi
      USE SEAICE, only : fsss,tfrez,Ei ! number, functions - ok
      !USE MODEL_COM, only : focean
      USE JulianCalendar_mod, only: jdendofm
      USE MODEL_COM, only: dtsrc
     *  ,itime,iyear1,nday,modelEclock,aMON
      use TimeConstants_mod, only: INT_HOURS_PER_DAY, SECONDS_PER_HOUR,
     &      INT_DAYS_PER_YEAR, INT_SECONDS_PER_HOUR, INT_SECONDS_PER_DAY
     &      ,SECONDS_PER_DAY, DAYS_PER_YEAR
      USE HYCOM_DIM_GLOB
      USE HYCOM_DIM, only : aJ_0, aJ_1, aJ_0H, aJ_1H,
     &                      aI_0, aI_1, aI_0H, aI_1H,
     &                       J_0,  J_1,  J_0H,  J_1H,
     &     ntrcr_age,ntrcr_vent,ntrcr_wm,ntrcr_zebra,ntrcr_obio,
     &      ogrid, isp_loc => isp, ifp_loc => ifp, ilp_loc => ilp
      USE HYCOM_SCALARS
      USE HYCOM_ARRAYS_GLOB
      USE hycom_arrays_glob_renamer
c
      USE KPRF_ARRAYS
      USE HYCOM_CPLER
      USE HYCOM_DYNSI_CPLER
      use mdul_glbdia,only: glbsum,glob2d

      use TimerPackage_mod
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      USE Dictionary_mod, only : get_param
      implicit none
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: iceocn
      type(iceocn_xchng_vars) :: dynsice
c
      integer, parameter :: INT_HOURS_PER_YEAR = INT_HOURS_PER_DAY *
     &                                           INT_DAYS_PER_YEAR
      integer i,j,k,l,m,n,mm,nn,km,kn,k1m,k1n,ia,ja,jb,iam1
!!! afogcm,nsavea should be initialized properly !
      integer :: afogcm=0,nsavea=0,nsaveo
#include "kprf_scalars.h"
c
      real sum_,coord,x,x1,totl,sumice,fusion,saldif,tf,sigocn
     .  ,chk_rho,chk_rhor,chk_rhostar,apehyc,pechg_hyc_bolus
     .  ,hyc_pechg1,hyc_pechg2,q,sum1,sum2,dpini(kdm)
     .  ,thkchg,flxdiv,eflow_gl,sigloc,t4,s4,p4,sigstar
      real*8 sss_restore_dt,sss_restore_dtice
      real totlj(J_0H:J_1H), sumj(J_0H:J_1H), sumicej(J_0H:J_1H)
      integer jj1,no,index,nflip,mo0,mo1,mo2,mo3,rename,iatest,jatest
     .       ,OMP_GET_NUM_THREADS,io,jo,nsub
      integer ipa_loc(aI_0H:aI_1H,aJ_0H:aJ_1H)
#ifdef USE_ATM_GLOBAL_ARRAYS
      integer ipa(iia,jja)
#endif
      integer nt
      external rename
      logical master,slave,diag_ape
c     logical :: chksum=.true.         ! compute global sums after each subr
      logical :: chksum=.false.
      logical :: vrbos=.false.
      character util(idm*jdm+14)*2,charac(20)*1,string*20,
     .          flnm*60
      data charac/'1','2','3','4','5','6','7','8','9','0',
     .            'A','B','C','D','E','F','G','H','I','J'/
      data iatest,jatest/31,41/           ! Iceland
      data nflip/0/
c
      real cnuity_time,tsadvc_time,momtum_time,barotp_time,trcadv_time,
     .     convec_time,thermf_time,enloan_time,mxlayr_time,hybgen_time,
     .     agcm_time,ogcm_time
      integer after,before,rate,bfogcm
      integer :: year, month, dayOfYear, date, hour
c
      real, dimension(idm,jdm) :: osst,osss,odhsi2,oogeoza,usf,vsf
      real, dimension(idm,J_0H:J_1H) :: osst_loc,osss_loc,odhsi2_loc,
     & oogeoza_loc,usf_loc,vsf_loc,tauxi_loc,tauyi_loc,ustari_loc
      real, dimension(aI_0H:aI_1H,aJ_0H:aJ_1H) :: utila
#ifdef TRACERS_GASEXCH_ocean
      real, allocatable :: otrac_loc(:,:,:) !(idm,J_0H:J_1H,ntm)
#endif

#include "state_eqn.h"

      real*8, dimension(:,:), pointer :: rsi_loc,focean_loc,
     &     dmua_loc,dmva_loc,e0_loc,egmelt_loc,eflowo_loc,
     &     prec_loc,eprec_loc,evapor_loc,flowo_loc,gmelt_loc,
     &     melti_loc,emelti_loc,smelti_loc,
     &     runosi_loc,erunosi_loc,   srunosi_loc,
     &     runpsi_loc,erunpsi_loc,srunpsi_loc,mlhc_loc,
     &     sss_loc,ogeoza_loc,uosurf_loc,vosurf_loc,gtemp_loc,gtempr_loc
      real*8, dimension(:,:,:), pointer :: dmsi_loc,dhsi_loc,dssi_loc
      real*8, dimension(:,:), pointer :: cosz1_loc,wsavg_loc,achl_loc
!     real*8, dimension(:,:), allocatable, save ::
!    &     avisdir_loc,avisdif_loc, anirdir_loc, anirdif_loc
!#ifdef TRACERS_OceanBiology /*TNL insert this line */
!      real*8, allocatable, DIMENSION(:,:), save :: asolz_loc, awind_loc
!#endif   /*TNL insert this line */
#ifdef TRACERS_GASEXCH_ocean /*TNL insert this line */
      real*8, dimension(:,:,:), pointer :: GTRACER_loc
#endif   /*TNL insert this line */
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- initiate named-pipe comparison utility
c --- (see pipe.f for instructions)
ccc      inquire(file='master',exist=master)
ccc      inquire(file='slave',exist=slave)
ccc      if ((master .and. slave) .or. (.not.master .and. .not.slave))
ccc     .   stop '(master/slave ambiguity)'
ccc      if (master) call pipe_init(.true.)
ccc      if (slave) call pipe_init(.false.)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      ! move to global atm grid
      call start('hycom')
c
      call modelEclock%get(year=year, month=month, date=date,
     .  hour=hour, dayOfYear=dayOfYear)

      call getdte(Itime,Nday,Iyear1,year,month,dayOfYear,date,hour,amon)
!#ifdef TRACERS_OceanBiology /*TNL insert this line */
!      if (.not.allocated(asolz_loc)) then
!        allocate(asolz_loc(aI_0H:aI_1H,aJ_0H:aJ_1H))
!        allocate(awind_loc(aI_0H:aI_1H,aJ_0H:aJ_1H))
!      endif
!#endif /*TNL insert this line */
      rsi_loc => iceocn%rsi
      focean_loc => atmocn%focean
      dmua_loc => atmocn%dmua
      dmva_loc => atmocn%dmva
      e0_loc => atmocn%e0
      evapor_loc => atmocn%evapor
      flowo_loc => atmocn%flowo
      eflowo_loc => atmocn%eflowo
      gmelt_loc => atmocn%gmelt
      egmelt_loc => atmocn%egmelt
      prec_loc => atmocn%prec
      eprec_loc => atmocn%eprec
      runosi_loc => iceocn%runosi
      erunosi_loc => iceocn%erunosi
      srunosi_loc => iceocn%srunosi
      runpsi_loc => iceocn%runpsi
      erunpsi_loc => iceocn%erunpsi
      srunpsi_loc => iceocn%srunpsi
      melti_loc => iceocn%melti
      emelti_loc => iceocn%emelti
      smelti_loc => iceocn%smelti
      cosz1_loc => atmocn%cosz1
      wsavg_loc => atmocn%wsavg
      achl_loc => atmocn%chl

      gtemp_loc => atmocn%gtemp
      gtempr_loc => atmocn%gtempr
      sss_loc => atmocn%sss
      ogeoza_loc => atmocn%ogeoza
      uosurf_loc => atmocn%uosurf
      vosurf_loc => atmocn%vosurf
        mlhc_loc => atmocn%mlhc
#ifdef TRACERS_ON
      gtracer_loc => atmocn%gtracer
#endif
      dmsi_loc => iceocn%dmsi
      dhsi_loc => iceocn%dhsi
      dssi_loc => iceocn%dssi

      if (atmocn%gasex_index%getsize()>0) trcout=.true.
#ifdef STANDALONE_OCEAN
      ! arrays for salinity restoring
      ! to get the restoring timescales (days) specified in the rundeck:
      call get_param("sss_restore_dt",sss_restore_dt) ! open water
      call get_param("sss_restore_dtice",sss_restore_dtice) ! under ice
#endif

cdiag write(*,'(a,i8,7i5,a)')'chk=',Itime,Nday,Iyear1,Jyear,Jmon
cdiag.        ,Jday,Jdate,Jhour,amon
c
      if (mod(hour,nhr).eq.0 .and.
     &    mod(itime,nday/INT_HOURS_PER_DAY).eq.0) then
        do 28 ja=aJ_0,aJ_1
        do 28 ia=aI_0,aI_1
          ataux_loc(ia,ja)=0.
          atauy_loc(ia,ja)=0.
        aflxa2o_loc(ia,ja)=0.
          aemnp_loc(ia,ja)=0.
          asalt_loc(ia,ja)=0.
           aice_loc(ia,ja)=0.
         austar_loc(ia,ja)=0.
         aswflx_loc(ia,ja)=0.
#ifdef TRACERS_GASEXCH_ocean
        do nt=1,atmocn%gasex_index%getsize()
          atracflx_loc(ia,ja,nt)=0.
        enddo
#endif
!#ifdef TRACERS_OceanBiology
!        awind_loc(ia,ja)=0.
!        asolz_loc(ia,ja)=0.
!#endif
!#ifdef OBIO_RAD_coupling
!        if (allocated(atmocn%dirvis)) then
!          if (.not.allocated(avisdir_loc))
!     &             ALLOCATE(
!     &         avisdir_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
!     &         avisdif_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
!     &         anirdir_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
!     &         anirdif_loc(aI_0H:aI_1H,aJ_0H:aJ_1H) )
!          avisdir_loc(ia,ja)=0.
!          avisdif_loc(ia,ja)=0.
!          anirdir_loc(ia,ja)=0.
!          anirdif_loc(ia,ja)=0.
!        endif
!#endif
 28     continue
#ifdef CUBED_SPHERE
        call reset_dynsi_accum
#endif
      endif
c
c --- accumulate agcm fields over nhr

#ifdef CUBED_SPHERE
! wind and ice stresses will be combined after regridding,
      admui_loc = 0. ! so set atm-grid instance of ice stress to zero
      admvi_loc = 0.
      call do_dynsi_accum(dtsrc,SECONDS_PER_HOUR*real(nhr))
#else
#ifdef STANDALONE_HYCOM
      admui_loc = 0.
      admvi_loc = 0.
#else
c combine wind and ice stresses before regridding
c first redistribute ice-ocean stresses to the atmospheric domain decomp
      call band_pack(dynsice%pack_i2a, dynsice%dmui, admui_loc)
      call band_pack(dynsice%pack_i2a, dynsice%dmvi, admvi_loc)
#endif
#endif
c
CTNL: two next lines from AR5_v2_branch
CTNL  eflow_gl=0.
CTNL  call globalsum_atm(grid, eflowo_loc, eflow_gl, all=.true.)
      eflow_gl = atmocn%eflow_gl
c
      do 29 ia=aI_0,aI_1
#if (defined CUBED_SPHERE) || defined(STANDALONE_HYCOM)
      iam1 = ia ! avoid out-of-bounds indices in admui_loc==0
#else
      iam1=mod(ia-2+iia,iia)+1
#endif
      do 29 ja=aJ_0,aJ_1
      if (focean_loc(ia,ja).eq.0.) then	! set huge to land points
        ipa_loc(ia,ja)=0
       aflxa2o_loc(ia,ja)=huge
         aemnp_loc(ia,ja)=huge
        austar_loc(ia,ja)=huge
        aswflx_loc(ia,ja)=huge
         asalt_loc(ia,ja)=huge
          aice_loc(ia,ja)=huge
         ataux_loc(ia,ja)=huge
         atauy_loc(ia,ja)=huge
      else
        ipa_loc(ia,ja)=1.
c --- accumulate
        aflxa2o_loc(ia,ja)=aflxa2o_loc(ia,ja)                      ! J/m2 => W/m2
     . +((e0_loc(ia,ja)+eprec_loc(ia,ja))*(1.-rsi_loc(ia,ja))      ! ocean water
     . + egmelt_loc(ia,ja)+emelti_loc(ia,ja)/focean_loc(ia,ja)     ! ocn or ice
     . + eflow_gl/area
     . +(erunosi_loc(ia,ja)+erunpsi_loc(ia,ja))*rsi_loc(ia,ja))    ! ice
     . /(SECONDS_PER_HOUR*real(nhr))
        aemnp_loc(ia,ja)=aemnp_loc(ia,ja)				! kg/m2 => m/s
     .  +((prec_loc(ia,ja)-evapor_loc(ia,ja))*(1.-rsi_loc(ia,ja))	! open water
     .  +flowo_loc(ia,ja)+gmelt_loc(ia,ja)				! river/glacial
     .  +melti_loc(ia,ja)/focean_loc(ia,ja)				! meltice
     .  +(runosi_loc(ia,ja)+runpsi_loc(ia,ja))*rsi_loc(ia,ja))*thref	! ice runoff
     .  /(SECONDS_PER_HOUR*real(nhr))
        austar_loc(ia,ja)=austar_loc(ia,ja)+(
     . sqrt(sqrt((dmua_loc(ia,ja)*(1d0-rsi_loc(ia,ja))
     .          +(admui_loc(ia,ja)+admui_loc(iam1,ja))*.5)**2
     .          +(dmva_loc(ia,ja)*(1d0-rsi_loc(ia,ja))
     .          +(admvi_loc(ia,ja)+admvi_loc(ia,max(1,ja-1)))*.5)**2)
     .          /dtsrc*thref))                      ! sqrt(T/r)=>m/s
     .          *dtsrc/(real(nhr)*SECONDS_PER_HOUR)
        aswflx_loc(ia,ja)=aswflx_loc(ia,ja)+(atmocn%solar(ia,ja)*
     .                         (1.-rsi_loc(ia,ja))          !J/m*m=>W/m*m
     .                         +iceocn%solar(ia,ja)*rsi_loc(ia,ja))
     .                         /(SECONDS_PER_HOUR*real(nhr))
        asalt_loc(ia,ja)=asalt_loc(ia,ja)				! salflx kg/m2 => kg/m2/sec
     . +( smelti_loc(ia,ja)/focean_loc(ia,ja)				! meltice
     . +(srunosi_loc(ia,ja)+srunpsi_loc(ia,ja))*rsi_loc(ia,ja))		! ice runoff
     . /(SECONDS_PER_HOUR*real(nhr))
        aice_loc(ia,ja)= aice_loc(ia,ja) + rsi_loc(ia,ja)*
     .                             dtsrc/(real(nhr)*SECONDS_PER_HOUR)
c --- dmua on A-grid, admui on C-grid
        ataux_loc(ia,ja)=ataux_loc(ia,ja)
     .               +(dmua_loc(ia,ja)*(1d0-rsi_loc(ia,ja))
     .               +(admui_loc(ia,ja)+admui_loc(iam1,ja))*.5)          ! scaled by rsi
     .               /(SECONDS_PER_HOUR*real(nhr))                       ! kg/ms => N/m2
        atauy_loc(ia,ja)=atauy_loc(ia,ja)
     .               +(dmva_loc(ia,ja)*(1d0-rsi_loc(ia,ja))
     .               +(admvi_loc(ia,ja)+admvi_loc(ia,max(1,ja-1)))*.5)   ! scaled by rsi
     .               /(SECONDS_PER_HOUR*real(nhr))                       ! kg/ms => N/m2
      end if
#ifdef TRACERS_GASEXCH_ocean
      do nt=1,atmocn%gasex_index%getsize()
        atracflx_loc(ia,ja,nt)=atracflx_loc(ia,ja,nt)
     .       + atmocn%TRGASEX(nt,ia,ja) ! in mol/m2/s
     .       * dtsrc/(real(nhr)*SECONDS_PER_HOUR)

        if (ia == itest .and. ja == jtest) then
          write(*,'(a,4i5,2e12.4)')'hycom, atracflx: ',
     .          itime,hour,ia,ja,atmocn%TRGASEX(nt,ia,ja),
     .                       atracflx_loc(ia,ja,nt)
        end if
      enddo
#endif
!#ifdef TRACERS_OceanBiology
!            asolz_loc(ia,ja)=asolz_loc(ia,ja) !
!     .           +COSZ1_loc(ia,ja)*dtsrc/(SECONDS_PER_HOUR*real(nhr)) !
!            awind_loc(ia,ja)=awind_loc(ia,ja) !
!     .           +wsavg_loc(ia,ja)*dtsrc/(SECONDS_PER_HOUR*real(nhr)) !
!#endif
!#ifdef OBIO_RAD_coupling
!          if (allocated(atmocn%dirvis)) then
!            avisdir_loc(ia,ja)=avisdir_loc(ia,ja) !
!     .         +atmocn%dirvis(ia,ja)*dtsrc/(SECONDS_PER_HOUR*real(nhr)) !
!            avisdif_loc(ia,ja)=avisdif_loc(ia,ja) !
!     .         +atmocn%difvis(ia,ja)*dtsrc/(SECONDS_PER_HOUR*real(nhr)) !
!            anirdir_loc(ia,ja)=anirdir_loc(ia,ja) !
!     .         +atmocn%dirnir(ia,ja)*dtsrc/(SECONDS_PER_HOUR*real(nhr)) !
!            anirdif_loc(ia,ja)=anirdif_loc(ia,ja) !
!     .         +atmocn%difnir(ia,ja)*dtsrc/(SECONDS_PER_HOUR*real(nhr)) !
!          endif
!#endif
 29   continue
c
      nsavea=nsavea+1

      if (mod(hour,nhr).gt.0 .or.
     &    mod(itime,nday/INT_HOURS_PER_DAY).eq.0) goto 9666

      if (nsavea*INT_HOURS_PER_DAY/nday.ne.nhr) then
        write(*,'(a,4i4,i8)')
     .  'nonmatching b.c. accumulation periods: agcm/ogcm:',
     .  nsavea*INT_HOURS_PER_DAY/nday,nhr,nsavea,nday,itime
        stop 'agcm/ogcm accumulating periods differ'
      end if
      nsavea=0
c
c --- pass ice/atmo fields to ocean

      call vec_a2o(  ataux_loc,  atauy_loc,taux_loc,tauy_loc) !wind stress
      call fld_a2o(   aice_loc,   oice_loc,'icecv')	      !ice coverage
      call fld_a2o(aflxa2o_loc,oflxa2o_loc,'surfl')	      !heatflux everywhere
      call fld_a2o(  asalt_loc,  osalt_loc,'osalt')	      !ice->ocn saltflux
      call fld_a2o(  aemnp_loc,  oemnp_loc,'eminp')           !E - P everywhere
      call fld_a2o( austar_loc,  ustar_loc,'ustar')           !friction velocity
      call fld_a2o( aswflx_loc, sswflx_loc,'shtwv')           !shortwave flux
#ifdef CUBED_SPHERE
c combine wind and ice stresses after regridding
      tauxi_loc = 0.; tauyi_loc = 0.; ustari_loc = 0.
      call vec_i2o(taux_dynsi,tauy_dynsi,tauxi_loc,tauyi_loc)
      call fld_i2o(ustar_dynsi,ustari_loc)
      do j=j_0,j_1
        do i=1,idm
          taux_loc(i,j) = taux_loc(i,j) + tauxi_loc(i,j)
          tauy_loc(i,j) = tauy_loc(i,j) + tauyi_loc(i,j)
          ustar_loc(i,j) = ustar_loc(i,j) + ustari_loc(i,j)
        enddo
      enddo
#endif
#ifdef TRACERS_GASEXCH_ocean
      do nt=1,atmocn%gasex_index%getsize()
        !call flxa2o(atracflx_loc(:,:,nt),tracflx_loc(:,:,nt)) !tracer flux
        call flxa2o(atracflx_loc(:,:,nt),ocnatm%work1)
        ocnatm%trgasex(nt,:,:) = ocnatm%work1(:,:)
      enddo
#endif
!#ifdef TRACERS_OceanBiology
!      call fld_a2o(asolz_loc,osolz_loc)
!      call fld_a2o(awind_loc,owind_loc)
!#endif
!#ifdef OBIO_RAD_coupling
!      call fld_a2o(avisdir_loc,ovisdir_loc)
!      call fld_a2o(avisdif_loc,ovisdif_loc)
!      call fld_a2o(anirdir_loc,onirdir_loc)
!      call fld_a2o(anirdif_loc,onirdif_loc)
!#endif
      call scatter1_hycom_arrays ! delete this call, if not the routine
c
      call system_clock(before)
      call system_clock(count_rate=rate)
      bfogcm=before
      agcm_time = real(bfogcm-afogcm)/real(rate)

      if (AM_I_ROOT()) then
      write (*,99009) agcm_time,' sec for AGCM bfor ocn step ',nstep
99009 format (f12.3,a,i8)
      endif ! AM_I_ROOT
c
      if (nstep.eq.0 .or. nstep.eq.nstep0) then
c
      jchunk=50
      if (jchunk.lt.0) then
        if (AM_I_ROOT()) write (*,*) 'you forgot to define jchunk'
        stop   ! proper mpi_abort
      end if
c
c --- compute eqn.of state check values
c
      t4=5.
      s4=36.
      p4=2500.e4			! 2500m
      chk_rho =39.6858			! in_situ rho for t4,s4 & p4
      if (pref.eq.0.) then
        chk_rhor    = 28.469            ! sigma0  for t4 & s4
        chk_rhostar = -9999.            ! sigma0* for t4 & s4
      elseif (pref.eq.1.e7) then
        chk_rhor    = 33.0321           ! sigma1  for t4 & s4
        chk_rhostar = 32.82078          ! sigma1* for t4 & s4
      elseif (pref.eq.2.e7) then
        chk_rhor    = 37.4934           ! sigma2  for t4 & s4
        chk_rhostar = -9999.            ! sigma2* for t4 & s4
      else
        stop 'wrong pref'    ! proper mpi_abort
      endif
c
      if (abs(sigocn(t4,s4)-chk_rhor) .gt. .001) then
      if (AM_I_ROOT())
     . write (*,*) 'error -- sigocn should be',chk_rhor
     . ,', not',sigocn(t4,s4)
        stop
      end if

      if (abs(sigstar(t4,s4,p4)-chk_rhostar).gt. .001) then
      if (AM_I_ROOT())
     .  write (*,*)'error: chk_rhostar should be',
     .  chk_rhostar,', not',sigstar(t4,s4,p4)
        stop
      end if
c
      if (abs(sigloc(t4,s4,p4)-chk_rho) .gt. .001) then
      if (AM_I_ROOT())
     . write (*,*) 'error -- sigloc should be',chk_rho
     . ,', not',sigloc(t4,s4,p4)
        stop
      end if
c
      if (AM_I_ROOT())
     & write (*,109) thkdff,temdff,veldff,viscos,diapyc,vertmx
 109  format (' turb. flux parameters:',1p/
     .  ' thkdff,temdff,veldff =',3e9.2/
     .  ' viscos,diapyc,vertmx =',3e9.2)
c
c --- 'lstep' = number of barotropic time steps per baroclinic time step.
c --- lstep   m u s t   be even.
c
      lstep=baclin/batrop
      lstep=2*((lstep+1)/2)
      dlt=baclin/lstep
      nstepi=real(nhr)*SECONDS_PER_HOUR/baclin + .0001

      if (AM_I_ROOT()) then
c
      write (*,'(i4,'' barotropic steps per baroclinic time step'')')
     .  lstep
      write (*,*) "time0/nstep0=",time0,nstep0
      write (*,'(a,i4,a,i4,a)') "ogcm exchange w. agcm every",nstepi,
     .   " steps, i.e.",nhr," hr"
      if( itest >= lbound(latij,dim=1) .AND.
     .    itest <= ubound(latij,dim=1) .AND.
     .    jtest >= lbound(latij,dim=2) .AND.
     .    jtest <= ubound(latij,dim=2) ) then
      write (*,*) "itest,jtest=",itest,jtest
      write (*,*) "lat,lon=",latij(itest,jtest,3),lonij(itest,jtest,3)
      endif

      endif ! AM_I_ROOT

c
c --- set up parameters defining the geographic environment
c
css   call geopar                ! moved to agcm for zero start or below
css   call inicon                ! moved to agcm
c
      mixfrq=int(43200./baclin+.0001)
      trcfrq=int(43200./baclin+.0001)
      if (AM_I_ROOT())
     &  write (*,'(2(a,i4),a)') 'mixfrq/trcfrq set to'
     &   ,mixfrq,'/',trcfrq,' steps'

      if (trcout) dotrcr=.true.
c
      if (year-iyear1.le.20) then
        if (AM_I_ROOT())
     &   write (*,'(''starting date in ogcm/agcm '',2i5,'' hr '',2i12
     &   /'' iyear1/year='',2i5)')
         ! crashes if nstep0 too big
CTNL .   int((nstep0+nstepi-1)*baclin)/(3600*24),itime/nday
     &   int((nstep0+nstepi-1)/(INT_SECONDS_PER_HOUR*INT_HOURS_PER_DAY)
     &   *baclin),itime/nday,
CTNL .   int((nstep0+nstepi-1)*baclin)/3600,itime*24/nday,iyear1,year ! crashes
     &   int((nstep0+nstepi-1)/INT_SECONDS_PER_HOUR*baclin),
     &   itime*INT_HOURS_PER_DAY/nday,iyear1,year
      else
        if (AM_I_ROOT())
     &   write (*,'(''starting date in ogcm/agcm '',2i12
     &    /'' iyear1/year='',2i5)')
CTNL .    int((nstep0+nstepi-1)*baclin)/
CTNL .    (INT_SECONDS_PER_HOUR*INT_HOURS_PER_DAY),itime/nday,
     &   int((nstep0+nstepi-1)/(INT_SECONDS_PER_HOUR*INT_HOURS_PER_DAY)
     &   *baclin),itime/nday,
     &    iyear1,year,(nstep0+nstepi-1)*baclin/SECONDS_PER_HOUR,
     &    itime*INT_HOURS_PER_DAY/nday
      endif
c
c     if(abs((nstep0+nstepi-1)*baclin/3600.-itime*24./nday).gt.1.e-5)
      if(abs(int((nstep0+nstepi-1)*baclin/SECONDS_PER_HOUR)-
     &   itime*INT_HOURS_PER_DAY/nday).gt.0)  then
        if (AM_I_ROOT()) then
          write (*,'(a,f16.8,i10,f16.8)')'mismatch date found ',
     &     int((nstep0+nstepi-1)*baclin/INT_SECONDS_PER_HOUR),
     &     itime*INT_HOURS_PER_DAY/nday,
     &     (nstep0+nstepi-1)*baclin/INT_SECONDS_PER_HOUR-itime*
     &     INT_HOURS_PER_DAY/nday
          write (*,'(/(a,i9,a,i9,a,f7.1,a))')'chk model start step',
     &     nstep0, ' changed to: ',
     &     int((itime*INT_HOURS_PER_DAY/nday+(iyear1-year)*
     &     INT_HOURS_PER_YEAR)*INT_SECONDS_PER_HOUR/baclin)-nstepi+1,
     &     ' (day ',(int((itime*INT_HOURS_PER_DAY/nday+(iyear1-year)*
     &     INT_HOURS_PER_YEAR)*INT_SECONDS_PER_HOUR/baclin)-nstepi+1)*
     &     baclin/INT_SECONDS_PER_DAY,')'
        end if   ! AM_I_ROOT

        nstep0=int((itime*INT_HOURS_PER_DAY/nday+(iyear1-year)*
     &   INT_HOURS_PER_YEAR)*INT_SECONDS_PER_HOUR/baclin)-nstepi+
     &   INT_SECONDS_PER_HOUR/baclin
        nstep=nstep0
        time0=nstep0*baclin/INT_SECONDS_PER_DAY
        tavini=time0
c
! dpini is used to compute thkchg which is used only for printing
        if (haveLatitude(ogrid, J=jtest)) then
          do 19 k=1,kk
 19         dpini(k)=dp_loc(itest,jtest,k)+dp_loc(itest,jtest,k+kk)
        end if
c
      else
        if (AM_I_ROOT()) then
          write (*,'(/(a,i9))') 'chk model starts at steps',nstep
        end if
      endif
c
      endif       ! if (nstep=0 or nstep=nstep0)
c
      if (nstepi.le.0) stop 'wrong nstepi'
c
c     if (AM_I_ROOT()) then
c     print *,' shown below oice at nstep=',nstep
c     call zebra(oice,iio,iio,jjo)
c     print *,' shown below aflxa2o'
c     call zebra(aflxa2o,iia,iia,jja)
c     end if
c
      osst=0.  ;   osst_loc=0.;
      osss=0.  ;   osss_loc=0.;
      odhsi2=0. ;  odhsi2_loc=0.;

      CALL HALO_UPDATE(ogrid,corio_loc,     FROM=NORTH)
      hekman_loc=0.;
c
      do 202 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 202 l=1,isp_loc(j)
      do 202 i=ifp_loc(j,l),ilp_loc(j,l)
c
      hekman_loc(i,j)=ustar_loc(i,j)*(cekman*4.0)/                ! kpp
     &           (abs(corio_loc(i,j  ))+abs(corio_loc(i+1,j  ))+
     &            abs(corio_loc(i,jb ))+abs(corio_loc(i+1,jb )))
 202  continue
#ifdef TRACERS_GASEXCH_ocean
      allocate(otrac_loc(idm,J_0H:J_1H,atmocn%gasex_index%getsize()))
      otrac_loc(:,:,:) = 0.d0
#endif
c
c --- ---------------------
c --- sub loop starts here
c --- ---------------------
c
c --- letter 'm' refers to mid-time level (example: dp(i,j,km) )
c --- letter 'n' refers to old and new time level
c
c     open(202,file='ip.array',form='unformatted',status='unknown')
c     write(202) ip,iu,iv,ifp,ilp,isp
c     close(202)
c
      nsaveo=0
      do 15 nsub=1,nstepi
      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
      mm=(m-1)*kk
      nn=(n-1)*kk
      k1m=1+mm
      k1n=1+nn
      nstep=nstep+1
      time=time0+(nstep-nstep0)*baclin/SECONDS_PER_DAY
c
      diagno=.false.
      if (JDendOfM(month) .eq. dayOfYear .and.hour .eq. 23 .and
     .  .nsub .eq. nstepi)  diagno = .true. ! end of month
c
      if (nstep.eq.1) then
        call gather_before_archive()
        if (AM_I_ROOT()) then
          call archiv(n,nn)
        endif
        call set_data_after_archiv()
      endif
c
      diag_ape=.false.
c
      trcadv_time = 0.0
c
      if (dotrcr) then
        call system_clock(before)      ! time elapsed since last system_clock
c
c --- initialization of tracer transport arrays (incl. dpinit):
        call tradv0(m,mm)
        dotrcr=.false.
c
c --- inject tracer into ogcm
c     write(*,*) 'nstep=',nstep
c     call findmx(ip,temp(1,1,k1n),idm,ii,jj,'sst')
c     call findmx(ip,saln(1,1,k1n),idm,ii,jj,'sss')
c     if (trcout) call cfcflx(tracer,p,temp(1,1,k1n),saln(1,1,k1n)
c    .           ,latij(1,1,3),scp2,baclin*trcfrq)
CTNL Why is it here?
c     if (AM_I_ROOT())
c    .write(*,'(a,6i3)') 'NTRCR = ',ntrcr,ntrcr_zebra,ntrcr_vent,
c    .                   ntrcr_age,ntrcr_wm,ntrcr_obio
#if defined(TRACERS_HYCOM_Ventilation) || defined(TRACERS_AGE_OCEAN) \
    || defined(TRACERS_OCEAN_WATER_MASSES) || defined(TRACERS_ZEBRA)
        do 12 j=J_0, J_1
        do 12 l=1,isp_loc(j)
        do 12 i=ifp_loc(j,l),ilp_loc(j,l)
#ifdef TRACERS_ZEBRA
           if (nstep.eq.1) then
            do nt=1,ntrcr_zebra
            do k=1,kdm
               if (k.eq.nt) then
                   tracer_loc(i,j,k,nt)=1.   !release in all layers
               else
                   tracer_loc(i,j,k,nt)=0.   !release in all layers
               endif
            enddo
            enddo
           endif
#endif
#ifdef TRACERS_AGE_OCEAN
!for consistency with Gary's ocean need to
!multiply here by trcfrq*baclin/(INT_DAYS_PER_YEAR*SECONDS_PER_DAY)
           nt = ntrcr_zebra + ntrcr_age
           tracer_loc(i,j,1,nt)=0.
           do k=2,kdm
           tracer_loc(i,j,k,nt)=tracer_loc(i,j,k,nt)
     .           + 1.*trcfrq*baclin/(INT_DAYS_PER_YEAR*SECONDS_PER_DAY)
           enddo
#endif
#ifdef TRACERS_OCEAN_WATER_MASSES
            !ACC
            nt=ntrcr_zebra+ntrcr_age+1
            tracer_loc(i,j,1,nt)=0.;
!             if (i.ge.289.and.i.le.386) then
!             if (j.ge.1.and.j.le.360) then
              if (latij(i,j,3).le.-55.d0) then
                tracer_loc(i,j,1,nt)=1.;
              endif
            !North Atlantic;
            nt=ntrcr_zebra+ntrcr_age+2
            tracer_loc(i,j,1,nt)=0.;
!             if (i.ge.100.and.i.le.190) then
!             if (j.ge.250.and.j.le.360) then
              if (latij(i,j,3).ge.40.d0) then
              if (lonij(i,j,3).ge.-90.d0.and.lonij(i,j,3).le.60.d0) then
                      tracer_loc(i,j,1,nt)=1.;
              endif
              endif
#endif
#ifdef TRACERS_HYCOM_Ventilation
        nt=ntrcr_zebra+ntrcr_age+ntrcr_wm+ntrcr_vent
        tracer_loc(i,j,1,nt)=1.              !  surface ventilation tracer
#endif
 12     continue
#endif
        call system_clock(after)        ! time elapsed since last system_clock
        trcadv_time = real(after-before)/real(rate)
      end if ! dotrcr  above if block is done only once
c
#ifdef TRACERS_OceanBiology
      call system_clock(before)      ! time elapsed since last system_clock
      trcout = .true.
c
c      if (dobio) then
c
          !call obio_listDifferences('obio_model', 'before')
        call obio_model(mm,ocnatm)
          !call obio_listDifferences('obio_model', 'after')
c
c      endif
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- available potential energy diagnostics:
      if (diag_ape) then
        call gather_hycom_arrays  !mkb check if all arrays are in
        if (AM_I_ROOT())
     .    write (501,'(f9.1,a,1p,e15.7)') time,'  APE (J):',
     .    hyc_pechg1(dp(1,1,k1n),th3d(1,1,k1n),31)
      end if ! diag_ape
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 103  format (f9.1,a,-12p,f9.3,' TW')
c
      if (chksum) call glbsum(mm,nn,' new timestep',0.,0.)

      call system_clock(before)
      call cnuity(m,n,mm,nn,k1m,k1n)
CTNL  call hycom_arrays_checksum     ! save CPU's time
      call system_clock(after)
      cnuity_time = real(after-before)/real(rate)

c     call sstbud(0,'  initialization',temp(1,1,k1n))
c
ccc      write (string,'(a12,i8)') 'cnuity, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      if (trcout) then
        before = after
c --- long time step tracer advection: build up mass flux time integral
        if (n.eq.oddev) then
          call tradv1(n,nn)
        endif
c
        if (mod(nstep,trcfrq).eq.0) then
          dotrcr=.true.         !  start tracer advection/turb.mixing cycle
          if (AM_I_ROOT())
     .    write (*,'(a)') 'start tracer advection/turb.mixing cycle'
          before = after
c --- tracer transport:
          call tradv2(n,nn)
        end if
        trcadv_time = trcadv_time + real(after-before)/real(rate)
      end if !trcout

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape) then
        call gather_hycom_arrays
        if (AM_I_ROOT())
     .    q=hyc_pechg1(dp(1,1,k1m),th3d(1,1,k1m),32)
      endif  ! diag_ape
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      before = after
      call tsadvc(m,n,mm,nn,k1m,k1n)
      if (chksum) call glbsum(mm,nn,'aft tsadvc',0.,0.)
cdiag if (AM_I_ROOT()) print *,'passed tsadvc'
      call system_clock(after)
      tsadvc_time = real(after-before)/real(rate)
c
c     call sstbud(1,' horiz.advection',temp(1,1,k1n))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape) then
        call gather_hycom_arrays
        if (AM_I_ROOT())
     .    write (501,103) time,'  APE change due to time smoothing: ',
     .    hyc_pechg2(dp(1,1,k1m),th3d(1,1,k1m),32)
      end if ! diag_ape
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      write (string,'(a12,i8)') 'tsadvc, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      call momtum(m,n,mm,nn,k1m,k1n)
cdiag if (AM_I_ROOT()) print *,'passed momtum'
      call system_clock(after)
      momtum_time = real(after-before)/real(rate)
ccc      write (string,'(a12,i8)') 'momtum, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      call barotp(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      barotp_time = real(after-before)/real(rate)
ccc      write (string,'(a12,i8)') 'barotp, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
c     if (iocnmx.eq.0.or.iocnmx.eq.2.or.iocnmx.eq.6)
      call convec(m,n,mm,nn,k1m,k1n)
      if (chksum) call glbsum(mm,nn,'aft convec',0.,0.)
c
      call system_clock(after)
      convec_time = real(after-before)/real(rate)
c     call sstbud(2,' convec.adjustmt',temp(1,1,k1n))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape) then
        call gather_hycom_arrays
        if (AM_I_ROOT())
     .  write (501,103) time,'  APE change due to mass adjustment:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
      end if ! diag_ape
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
ccc      write (string,'(a12,i8)') 'convec, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
c     before = after

      if ((iocnmx.eq.0.or.iocnmx.eq.2.or.iocnmx.eq.6) .and.
     .    nstep*baclin.ge.12.*INT_SECONDS_PER_HOUR)
     .    call diapfl(m,n,mm,nn,k1m,k1n)
      if (chksum) call glbsum(mm,nn,'aft diamix',0.,0.)
c
c     call system_clock(after)
c     diapfl_time = real(after-before)/real(rate)
c     call sstbud(3,'   diapyc.mixing',temp(1,1,k1n))
c
ccc      write (string,'(a12,i8)') 'diapfl, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      call thermf(m,n,mm,nn,k1m,k1n,
     .  sss_restore_dt,sss_restore_dtice)

      if (chksum) call glbsum(mm,nn,'aft thermf',0.,0.)
      call system_clock(after)
      thermf_time = real(after-before)/real(rate)
c
ccc      write (string,'(a12,i8)') 'thermf, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      call enloan(m,n,mm,nn,k1m,k1n)
      if (chksum) call glbsum(mm,nn,'aft enloan',0.,0.)
      call system_clock(after)
      enloan_time = real(after-before)/real(rate)
c
      before = after
c     if (nstep*baclin.le.12.*3600) then
c       call mxlayr(m,n,mm,nn,k1m,k1n)
c     else	! after 12 hrs
         call mxkprf(m,n,mm,nn,k1m,k1n)
         if (chksum) then
           h_glb_cum(m)=h_glb_cum(m)+glob2d(surflx_loc)*g*delt1
           s_glb_cum(m)=s_glb_cum(m)+glob2d(salflx_loc)*g*delt1
           call glbsum(mm,nn,'aft mxkprf',h_glb_cum(m),s_glb_cum(m))
         end if
c     end if
c     if (AM_I_ROOT())  print *,'passed mxkprf'
      call system_clock(after)
      mxlayr_time = real(after-before)/real(rate)

c     call sstbud(4,'  air-sea fluxes',tprime)
c     call sstbud(5,'     entrainment',temp(1,1,k1n))

      if (AM_I_ROOT()) then
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to thermal forcing:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      write (string,'(a12,i8)') 'mxlayr, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      end if ! AM_I_ROOT

      before = after
      call hybgen(m,n,mm,nn,k1m,k1n)
      if (chksum) call glbsum(mm,nn,'aft hybgen',0.,0.)
      call system_clock(after)
      hybgen_time = real(after-before)/real(rate)
c
ccc      write (string,'(a12,i8)') 'hybgrd, step',nstep
ccc      call comparall(m,n,mm,nn,string)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (AM_I_ROOT()) then
      if (diag_ape)
     . write (501,103) time,'  APE change due to vert.regridding:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c     call sstbud(6,' vert.regridding',temp(1,1,k1n))
c
      if (dotrcr)
     .  write (*,'(a)') 'tracer advection/turb.mixing cycle completed'
      end if ! AM_I_ROOT
c
c ---------------------------------------------------------------------------
c
c --- output and diagnostic calculations
c
c ---------------------------------------------------------------------------
c
c --- after 1 day, check terms in time-integrated continuity eqn
c
      if (abs(time-tavini-1.).lt..001) then
        CALL HALO_UPDATE(ogrid, vflxav_loc,  FROM=NORTH)
        i=itest
        j=jtest
        if (haveLatitude(ogrid, J=j)) then
          jb = PERIODIC_INDEX(j+1, jj)
          write (*,'(i9,2i5,a/a)') nstep,i,j,
     . '  time-integrated continuity eqn diagnostics:',
     .  '     thknss_tndcy  horiz.flxdiv   vert.flxdiv      residuum'
          do k=1,kk
            km=k+mm
            kn=k+nn
            thkchg=dp_loc(i,j,km)+dp_loc(i,j,kn)-dpini(k)
            flxdiv=(uflxav_loc(i+1,j,k)-uflxav_loc(i,j,k)
     .             +vflxav_loc(i,jb ,k)-vflxav_loc(i,j,k))*
     .             scp2i_loc(i,j)*delt1
            write (*,102) k,thkchg,flxdiv,-diaflx_loc(i,j,k),
     .      thkchg+flxdiv-diaflx_loc(i,j,k)
 102  format (i3,4f14.1)
          end do
        end if
      end if
c
!------------------------------------------------------------
      if (diagno .or. mod(time+.0001,1.).lt..0002) then    !  once a day
c
c --- diagnose mean sea surface height
        sumj(:)=0.; sumicej(:)=0.;
        do 704 j=J_0,J_1
        do 704 l=1,isp_loc(j)
        do 704 i=ifp_loc(j,l),ilp_loc(j,l)
c --- compute sea surface height (m)
        srfhgt_loc(i,j)=(montg_loc(i,j,1)+thref*pbavg_loc(i,j,m))/g
        sumicej(j)=sumicej(j)+oice_loc(i,j)*scp2_loc(i,j)
 704    sumj(j)=sumj(j)+srfhgt_loc(i,j)*scp2_loc(i,j)
        call GLOBALSUM(ogrid,sumicej,sumice, all=.true.)
        call GLOBALSUM(ogrid,sumj,sum_, all=.true.)

        if (AM_I_ROOT())
     .  write (*,'(i9,'' mean sea srf.hgt. (mm):'',f9.2,f12.0)')nstep,
     .  sum_*1.e3/area,sumice*1.e-6
c
c --- find the largest distance from a tencm layer from bottom
        do 707 j=J_0,J_1
        do 707 l=1,isp_loc(j)
        do 707 i=ifp_loc(j,l),ilp_loc(j,l)
        do 709 k=1,kk
        if (dp_loc(i,j,k+nn).lt.tencm) then
          util2_loc(i,j)=(p_loc(i,j,kk+1)-p_loc(i,j,k))/onem
          goto 708
        endif
 709    continue
 708    continue
 707    continue
c
c       call findmx(ip,util2,idm,ii,jj,'lowest')
c       call findmx(ip,osst,ii,ii,jj,'osst')
c       call findmx(ip,osss,ii,ii,jj,'osss')
c       call findmx(iu,usf,ii,ii,jj,'u_ocn')
c       call findmx(iv,vsf,ii,ii,jj,'v_ocn')
c
c       call findmx(ipa,aflxa2o,iia,iia,jja,'aflxa2o')
c       call findmx(ipa,asst,iia,iia,jja,'asst')
c       call findmx(ipa,sss,iia,iia,jja,'osss')
c       call findmx(ipa,uosurf,iia,iia,jja,'uosurf')
c       call findmx(ipa,vosurf,iia,iia,jja,'vosurf')
c
      end if                      ! once a day
!------------------------------------------------------------
c
c     if (AM_I_ROOT())
c    . write (*,'(2i5,2f8.1,a,f9.3,i4/(5(5f9.2,3x,5f9.2/)))')
c    . itest,jtest,latij(itest,jtest,3),lonij(itest,jtest,3),
c    .  ' output tx,ty,ht,ep,t,s,ice,dp,day=',time,nstep
c    . ,((taux_loc(i,j)*100.,j=jtest-2,jtest+2)
c    . , (tauy_loc(i,j)*100.,j=jtest-2,jtest+2),i=itest-2,itest+2)
c    . ,((oflxa2o_loc(i,j),j=jtest-2,jtest+2)
c    . , (oemnp_loc(i,j)*1.e6,j=jtest-2,jtest+2),i=itest-2,itest+2)
c    . ,((temp_loc(i,j,k1n),j=jtest-2,jtest+2)
c    . , (saln_loc(i,j,k1n),j=jtest-2,jtest+2),i=itest-2,itest+2)
c    . ,((oice_loc(i,j),j=jtest-2,jtest+2)
c    . , (dpmixl_loc(i,j,n)/onem,j=jtest-2,jtest+2),i=itest-2,itest+2)
c
c     write(401,'(i4,f7.3,2f6.1,f6.2,3f8.1)')
c    .  nstep,time,latij(itest,jtest,3),lonij(itest,jtest,3)
c    . ,oice_loc(i,j),dpmixl_loc(i,j,n)/onem,oflxa2o_loc(i,j),
c    .  oemnp_loc(i,j)*1.e6
c     do k=1,kk
c     kn=k+nn
c     write(401,'(i3,4f8.2)')
c    . k,p(i,j,k+1)/onem,temp(i,j,kn),saln(i,j,kn),th3d(i,j,kn)
c     enddo
c
      ogcm_time =
     .    cnuity_time
     .   +tsadvc_time
     .   +trcadv_time
     .   +momtum_time
     .   +barotp_time
     .   +thermf_time
     .   +enloan_time
     .   +mxlayr_time
     .   +hybgen_time

      if (AM_I_ROOT()) then
      write (*,99009) ogcm_time,' sec for OGCM   at ocn step ',nstep
      write (*,'(a/(5(4x,a,i5)))') 'timing (msec) by routine:'
     .   ,'cnuity',int(1000.*cnuity_time)
     .   ,'tsadvc',int(1000.*tsadvc_time)
     .   ,'trcadv',int(1000.*trcadv_time)
     .   ,'momtum',int(1000.*momtum_time)
     .   ,'barotp',int(1000.*barotp_time)
     .   ,'thermf',int(1000.*thermf_time)
     .   ,'enloan',int(1000.*enloan_time)
     .   ,'mxlayr',int(1000.*mxlayr_time)
     .   ,'hybgen',int(1000.*hybgen_time)
c
      if (mod(nstep,5).eq.0) call sys_flush(6)
      end if  ! AM_I_ROOT

c
      if (.not.diagno) go to 23
c
      if (AM_I_ROOT())
     .  write (*,100) nstep,int((time+.001)/DAYS_PER_YEAR),
     &                 mod(time+.001,DAYS_PER_YEAR)
 100  format (' ocn time step',i9,4x,'y e a r',i6,4x,'d a y',f9.2)
c
c --- output to history file
c
      call gather_before_archive()

      if (AM_I_ROOT()) then

        call archiv(n,nn)

c --- diagnose meridional overturning and heat flux
        do 3 j=1,jj
        do 3 k=1,kk
        do 3 l=1,isp(j)
        do 3 i=ifp(j,l),ilp(j,l)
 3      p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
c     call overtn(mm)
c
      write (*,105) nstep
 105  format (' step',i9,' -- archiving completed --')
c
c --- test for cyclic-in-j vs. noncyclic domain
c --- (closed-basin coditions are indicated by a land barrier at j=jj)
c
      jj1=jj-1
      do i=1,ii1
        if (ip(i,jj).gt.0) then
          jj1=jj
        end if
      end do

c
c --- output to line printer
c
ccc      call prtmsk(ip,srfhgt,util3,idm,ii1,jj1,0.,1./(thref*onecm),
ccc     .     'sea surface height (cm)')
ccc      call prtmsk(ip,dpmixl(1,1),util3,idm,ii1,jj1,0.,1./onem,
ccc     .     'mixed layer depth (m)')

      call prtmsk(ip,temp(1,1,k1n),util3,idm,ii1,jj1,0.,10.,
     .     'mix.layer temp. (.1 deg)')
ccc      call prtmsk(ip,saln(1,1,k1n),util3,idm,ii1,jj1,35.,100.,
ccc     .     'mx.lay. salin. (.01 mil)')

ccc      call prtmsk(ip,oice,util3,idm,ii1,jj1,0.,100.,
ccc.     'ice coverage (cm)')
ctmp#ifdef USE_ATM_GLOBAL_ARRAYS
ctmp        call prtmsk(ipa,asst,util3,iia,iia,jja,0.,1.,'asst ')
ctmp#endif

      do 77 j=1,jj1
      do 77 l=1,isu(j)
      do 77 i=ifu(j,l),ilu(j,l)
 77   util1(i,j)=u(i,j,k1n)+ubavg(i,j,n)
ccc      call prtmsk(iu,util1,util3,idm,ii1,jj1,0.,1000.,
ccc     .     'u vel. (mm/s), layer 1')

      do 78 i=1,ii1
      do 78 l=1,jsv(i)
      do 78 j=jfv(i,l),jlv(i,l)
 78   util2(i,j)=v(i,j,k1n)+vbavg(i,j,n)
ccc      call prtmsk(iv,util2,util3,idm,ii1,jj1,0.,1000.,
ccc     .     'v vel. (mm/s), layer 1')
ccc      call prtmsk(iu,ubavg(1,1,n),util3,idm,ii1,jj1,0.,1000.,
ccc     .     'barotrop. u vel. (mm/s)')
ccc      call prtmsk(iv,vbavg(1,1,n),util3,idm,ii1,jj1,0.,1000.,
ccc     .     'barotrop. v vel. (mm/s)')

      endif  !  AM_I_ROOT
      call set_data_after_archiv()

 23   continue

c --- accumulate fields for agcm

      do 201 j=J_0,J_1
      do 201 l=1,isp_loc(j)
      do 201 i=ifp_loc(j,l),ilp_loc(j,l)
      osst_loc(i,j)=osst_loc(i,j)+temp_loc(i,j,k1n)*baclin/
     &              (SECONDS_PER_HOUR*real(nhr))
      osss_loc(i,j)=osss_loc(i,j)+saln_loc(i,j,k1n)*baclin/
     &              (SECONDS_PER_HOUR*real(nhr))
      odhsi2_loc(i,j)=odhsi2_loc(i,j)+odhsi_loc(i,j)*baclin*dtsrc/
     &               (SECONDS_PER_HOUR*real(nhr))          !kg/m2=>kg*.5*hr/m2
      omlhc_loc(i,j)=spcifh*max(dp_loc(i,j,k1n)/onem,thkmin)/thref  ! J/(m2*C)
      oogeoza_loc(i,j)=(montg_loc(i,j,1)+thref*pbavg_loc(i,j,m))*
     &                 g/(thref*onem)                               ! m^2/s^2
 201  continue

#ifdef TRACERS_ON
      do j=J_0,J_1
      do l=1,isp_loc(j)
      do i=ifp_loc(j,l),ilp_loc(j,l)
      !define the tracer that participates in the gas exchange flux.
        do nt=1,atmocn%gasex_index%getsize()
          otrac_loc(i,j,nt)=otrac_loc(i,j,nt)
     .             + ocnatm%gtracer(ocnatm%gasex_index%at(nt),i,j) !pCO2_loc(i,j)  !pCO2 is in ppmv(uatm)
     .             * baclin/(SECONDS_PER_HOUR*real(nhr))
            enddo
      enddo
      enddo
      enddo
#endif
c
      nsaveo=nsaveo+1

      delt1=baclin+baclin

 15   continue

      if (nsaveo*baclin.ne.nhr*INT_SECONDS_PER_HOUR) then
        print *, ' ogcm saved over hr=',nsaveo*baclin/SECONDS_PER_HOUR
        stop ' stop: ogcm saved over hr'
      end if
      nsaveo=0

      call gather6hycom(
     &         osst_loc,osss_loc,odhsi2_loc,oogeoza_loc,
     &         osst,    osss,    odhsi2    ,oogeoza)
c
      do 88 j=J_0,J_1
      do 88 i=1,ii
      usf_loc(i,j)=u_loc(i,j,k1n)+ubavg_loc(i,j,n)
 88   vsf_loc(i,j)=v_loc(i,j,k1n)+vbavg_loc(i,j,n)

      call gather7hycom(usf_loc,vsf_loc,omlhc_loc,usf,vsf,omlhc)

c     if (AM_I_ROOT()) then ! global grid
c
c       call findmx(ip,osst,ii,ii,jj,'osst')
c       call findmx(ip,osss,ii,ii,jj,'osss')
c       call findmx(iu,usf,ii,ii,jj,'u_ocn')
c       call findmx(iv,vsf,ii,ii,jj,'v_ocn')
c
c     endif

c --- pass ocean fields to ice/atmo

      utila(:,:)=0.
      call fld_o2a(osst_loc,atmocn%work1,'s s t')
      call tempr_o2a(osst_loc,atmocn%work2)
      call fld_o2a(osss_loc,sss_loc,'s s s')
      call fld_o2a(omlhc_loc,mlhc_loc,'omlhc')
      call fld_o2a(oogeoza_loc,ogeoza_loc,'geoza')
c --- sign convention: ocn heat gain -> dhsi<0
      call fld_o2a(odhsi2_loc,utila,' dhsi')	! <0 J/m2 per agcm time step
      call vec_o2a(usf_loc,vsf_loc,uosurf_loc,vosurf_loc)
      odhsi_loc=0.;
#ifdef CUBED_SPHERE
      call vec_o2i(usf_loc,vsf_loc,dynsice%uosurf,dynsice%vosurf)
#else
#ifndef STANDALONE_HYCOM
c when atm is lat-lon, dynsi grid == atm grid
      call band_pack(dynsice%pack_a2i, uosurf_loc, dynsice%UOSURF)
      call band_pack(dynsice%pack_a2i, vosurf_loc, dynsice%VOSURF)
#endif
#endif
c
c     call findmx(ipa,asst,iia,iia,jja,'asst')
c     call findmx(ipa,sss,iia,iia,jja,'osss')
c     call findmx(ipa,uosurf,iia,iia,jja,'uosurf')
c     call findmx(ipa,vosurf,iia,iia,jja,'vosurf')
c
cdiag if (time.le.1)
cdiag call vec_a2o(ataux_loc,atauy_loc,taux_loc,tauy_loc) !wind stress
cdiag.write (*,l'(2i5,a,f9.3,i4/(5(5f9.2,3x,5f9.2/)))')
cdiag. iatest,jatest,' output raw uo,vo,ice,focean day=',time,nstep
cdiag. ,((uosurf(i,j)*100.,i=iatest-2,iatest+2)
cdiag. , (vosurf(i,j)*100.,i=iatest-2,iatest+2),j=jatest+2,jatest-2,-1)
cdiag. ,((aice(i,j)*100.,  i=iatest-2,iatest+2)
cdiag. , (focean(i,j)*100.,i=iatest-2,iatest+2),j=jatest+2,jatest-2,-1)
c
      dhsi_loc(:,:,:)=0.
      dmsi_loc(:,:,:)=0.
      dssi_loc(:,:,:)=0.
      do 204 ja=aJ_0,aJ_1
      do 204 ia=aI_0,aI_1
      if (focean_loc(ia,ja).gt.0.) then
        gtemp_loc(ia,ja)=atmocn%work1(ia,ja)
        gtempr_loc(ia,ja)=atmocn%work2(ia,ja)
        tf=tfrez(sss_loc(ia,ja),0.)
c --- 'dmsi' is the new mass (ice & salt) releasing thermal energy 'dhsi'
	dhsi_loc(1,ia,ja)=utila(ia,ja)           !<0 J/m2 per agcm step
c --- sign convention: ice mass gain -> dmsi>0
c --- sign convention: ice salt gain -> dssi>0
 	dmsi_loc(1,ia,ja)=dhsi_loc(1,ia,ja)/
     .    (Ei(tf,fsss*sss_loc(ia,ja))-tf*spcifh) !>0 kg/m2 per agcm step
c	dhsi_loc(1,ia,ja)=utila(ia,ja)*Ei(tf,fsss*sss_loc(ia,ja))
c --- 'dssi' is the salt amount sequestered in new ice (eq.9 in Schmidt et al)
        dssi_loc(1,ia,ja)=1.d-3*dmsi_loc(1,ia,ja)*fsss*sss_loc(ia,ja) !>0 kg/m2 for ice salt gain per agcm step
c --- evenly distribute new ice over open water and sea ice
c --- this is not necessarily a good idea. What about weighting it
c --- with respect to the ice/openwater flux ratio?
        if (aice_loc(ia,ja).gt.1.e-3) then
          dhsi_loc(2,ia,ja)=dhsi_loc(1,ia,ja)
          dmsi_loc(2,ia,ja)=dmsi_loc(1,ia,ja)
          dssi_loc(2,ia,ja)=dssi_loc(1,ia,ja)
        endif	! if aice > 0
      end if	! if focean>0
 204  continue

#ifdef TRACERS_GASEXCH_ocean
      do nt=1,atmocn%gasex_index%getsize()
         call ssto2a(otrac_loc(:,:,nt),atmocn%work1)
         do ja=aJ_0,aJ_1
         do ia=aI_0,aI_1
           if (focean_loc(ia,ja).gt.0.) then
             GTRACER_loc(nt,ia,ja) = atmocn%work1(ia,ja)
            !change units ppmv (uatm) -> kg,CO2/kg,air
     &            *atmocn%vol2mass(nt)* 1.d-6
           endif
         enddo
         enddo
      enddo
      deallocate(otrac_loc)
      if (ocnatm%chl_defined) then
      !call ssto2a(tot_chlo_loc,achl_loc)
        call ssto2a(ocnatm%chl,achl_loc)
        atmocn%chl_defined=.true.
      endif
#endif
c
      call system_clock(afogcm)
      if (abs(time-(itime+1.)/nday).gt..01) then
        write(*,'(a,2f10.3)') 'stop: agcm/ogcm date after ogcm'
     .    ,(itime+1.)/nday,time
        stop ' stop: dates in agcm/ogcm do not match'
      end if
c
 9666 continue

      call stop('hycom')
      return
      end
c
c> Revision history:
c>
c> May  1997 - removed statement "theta(1)=-thbase" after loop 14
c> June 1997 - added loop 60 to fix bug in vertical summation of -diaflx-
c> Mar. 2000 - conversion to SI units
c> June 2000 - added timing logic
c> Aug. 2000 - added code for changing number of OpenMP threads
c> Apr. 2001 - eliminated stmt_funcs.h
c> June 2001 - corrected sign error in diaflx
c> July 2001 - replaced archiving statements by 'call archiv'
c> Oct  2004 - map ice mass to agcm, and then calculate E
c------------------------------------------------------------------
      subroutine gather_before_archive
      USE HYCOM_ARRAYS_GLOB
      use hycom_arrays_glob_renamer
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
      implicit none

#ifdef TRACERS_OceanBiology
      call obio_gather_before_archive()
#endif

      call pack_data( ogrid,  u_loc, u )
      call pack_data( ogrid,  v_loc, v )
      call pack_data( ogrid,  dp_loc, dp )
      call pack_data( ogrid,  p_loc, p )
      call pack_data( ogrid,  temp_loc, temp )
      call pack_data( ogrid,  saln_loc, saln )
      call pack_data( ogrid,  th3d_loc, th3d )
      call pack_data( ogrid,  dpmixl_loc, dpmixl )
      call pack_data( ogrid,  srfhgt_loc, srfhgt )
      call pack_data( ogrid,  ubavg_loc, ubavg )
      call pack_data( ogrid,  vbavg_loc, vbavg )
      call pack_data( ogrid,  uav_loc, uav )
      call pack_data( ogrid,  vav_loc, vav )
      call pack_data( ogrid,  dpuav_loc, dpuav )
      call pack_data( ogrid,  dpvav_loc, dpvav )
      call pack_data( ogrid,  temav_loc, temav )
      call pack_data( ogrid,  salav_loc, salav )
      call pack_data( ogrid,  th3av_loc, th3av )
      call pack_data( ogrid,  dpav_loc, dpav )
      call pack_data( ogrid,  pbavav_loc, pbavav )
      call pack_data( ogrid,  sfhtav_loc, sfhtav )
      call pack_data( ogrid,  uflxav_loc, uflxav )
      call pack_data( ogrid,  vflxav_loc, vflxav )
      call pack_data( ogrid,  ufxavp_loc, ufxavp )
      call pack_data( ogrid,  vfxavp_loc, vfxavp )
      call pack_data( ogrid,  diaflx_loc, diaflx )
      call pack_data( ogrid,  salflav_loc, salflav )
      call pack_data( ogrid,  brineav_loc, brineav )
      call pack_data( ogrid,  eminpav_loc, eminpav )
      call pack_data( ogrid,  surflav_loc, surflav )
      call pack_data( ogrid,  tauxav_loc, tauxav )
      call pack_data( ogrid,  tauyav_loc, tauyav )
      call pack_data( ogrid,  dpmxav_loc, dpmxav )
      call pack_data( ogrid,  oiceav_loc, oiceav )
      call pack_data( ogrid,  pbot_loc, pbot )
      call pack_data( ogrid,  tracer_loc, tracer )
      call pack_data( ogrid,  oice_loc, oice )
      call pack_data( ogrid,  util1_loc, util1 )
      call pack_data( ogrid,  osalt_loc, osalt )

      end subroutine gather_before_archive
c------------------------------------------------------------------
      subroutine set_data_after_archiv

      USE HYCOM_DIM, only : ii, J_0,  J_1, kk, ogrid
      USE HYCOM_ARRAYS_GLOB, only : p_glob => p, util1_glob => util1,
     &                                           util2_glob => util2
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, ONLY: UNPACK_DATA

      implicit none
      integer :: i, j, k

      call unpack_data( ogrid,  p_glob, p )
      call unpack_data( ogrid,  util1_glob, util1 )
      call unpack_data( ogrid,  util2_glob, util2 )

      do 60 j=J_0, J_1
c
      do 601 i=1,ii
      if (dpav(i,j,1).gt.1.) then
        eminpav(i,j)=0.
        surflav(i,j)=0.
        salflav(i,j)=0.
        brineav(i,j)=0.
        tauxav(i,j)=0.
        tauyav(i,j)=0.
        pbavav(i,j)=0.
        dpmxav(i,j)=0.
        sfhtav(i,j)=0.
        oiceav(i,j)=0.
      end if
 601  continue
c
      do 60 k=kk,1,-1
      do 60 i=1,ii
      if (dpav(i,j,1).gt.1.) then
        uav(i,j,k)=0.
        vav(i,j,k)=0.
        dpuav(i,j,k)=0.
        dpvav(i,j,k)=0.
        dpav (i,j,k)=0.
        temav(i,j,k)=0.
        salav(i,j,k)=0.
        th3av(i,j,k)=0.
        uflxav(i,j,k)=0.
        vflxav(i,j,k)=0.
        ufxavp(i,j,k)=0.
        vfxavp(i,j,k)=0.
        diaflx(i,j,k)=0.
      end if
 60   continue

#ifdef TRACERS_OceanBiology
      call obio_set_data_after_archiv()
#endif

      end subroutine set_data_after_archiv
c------------------------------------------------------------------
      subroutine gather6hycom(
     &  osst_loc,osss_loc,odhsi2_loc,oogeoza_loc,  ! in
     &  osst,    osss,    odhsi2,    oogeoza)      ! out

      USE HYCOM_DIM, only : idm, jdm, J_0H,  J_1H, ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
      implicit none
      real,intent(IN) :: osst_loc(idm,J_0H:J_1H),
     &    osss_loc(idm,J_0H:J_1H), odhsi2_loc(idm,J_0H:J_1H),
     &    oogeoza_loc(idm,J_0H:J_1H)
      real,intent(OUT) :: osst(idm,jdm), osss(idm,jdm), odhsi2(idm,jdm),
     &    oogeoza(idm,jdm)

      call pack_data( ogrid,  osst_loc,      osst    )
      call pack_data( ogrid,  osss_loc,      osss    )
      call pack_data( ogrid,  odhsi2_loc,    odhsi2  )
      call pack_data( ogrid,  oogeoza_loc,   oogeoza )

      return
      end subroutine gather6hycom
c------------------------------------------------------------------
      subroutine gather7hycom(usf_loc,vsf_loc,omlhc_loc,usf,vsf,omlhc)

      USE HYCOM_DIM, only : idm, jdm, J_0H,  J_1H, ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
      implicit none
      real,intent(IN)  :: usf_loc(idm,J_0H:J_1H), vsf_loc(idm,J_0H:J_1H)
     .     ,omlhc_loc(idm,J_0H:J_1H)
      real,intent(OUT) :: usf(idm,jdm), vsf(idm,jdm), omlhc(idm,jdm)

      call pack_data( ogrid,  usf_loc,      usf    )
      call pack_data( ogrid,  vsf_loc,      vsf    )
      call pack_data( ogrid,omlhc_loc,    omlhc    )

      return
      end subroutine gather7hycom
c------------------------------------------------------------------
      subroutine hycom_arrays_checksum
      USE HYCOM_DIM, only: ogrid, J_0, J_1
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT, GLOBALSUM

      implicit none
      real :: arraySum
      integer :: iarraySum

      call GLOBALSUM( ogrid, sum(u,dim=3), arraySum )  !sum(u(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(v,dim=3), arraySum )  !sum(v(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(dp,dim=3), arraySum ) !sum(dp(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(dpold,dim=3), arraySum ) !sum(dpold(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(dpu,dim=3), arraySum )   !sum(dpu(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(dpv,dim=3), arraySum )   !sum(dpv(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(p,dim=3), arraySum )   ! sum(p(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(pu,dim=3), arraySum )   ! sum(pu(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(pv,dim=3), arraySum )   ! sum(pv(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(latij,dim=3), arraySum )   ! sum(latij(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(lonij,dim=3), arraySum )   ! sum(lonij(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, corio, arraySum )  ! sum(corio(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, potvor, arraySum )   ! sum(potvor(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(temp,dim=3), arraySum )   ! sum(temp(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(saln,dim=3), arraySum )   ! sum(saln(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(th3d,dim=3), arraySum )   ! sum(th3d(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(thstar,dim=3), arraySum )   ! sum(thstar(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, psikk, arraySum )   ! sum(psikk(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, thkk, arraySum )   ! sum(thkk(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(dpmixl,dim=3), arraySum )   ! sum(dpmixl(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, srfhgt, arraySum )   ! sum(srfhgt(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(montg,dim=3), arraySum )   ! sum(montg(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, defor1, arraySum )   ! sum(defor1(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, defor2, arraySum )   ! sum(defor2(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(ubavg,dim=3), arraySum )   ! sum(ubavg(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(vbavg,dim=3), arraySum )   ! sum(vbavg(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(pbavg,dim=3), arraySum )   ! sum(pbavg(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, ubrhs, arraySum )   ! sum(ubrhs(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, vbrhs, arraySum )   ! sum(vbrhs(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, utotm, arraySum )   ! sum(utotm(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, vtotm, arraySum )   ! sum(vtotm(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, utotn, arraySum )   ! sum(utotn(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, vtotn, arraySum )   ! sum(vtotn(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, uflux, arraySum )   ! sum(uflux(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, vflux, arraysum )   ! sum(vflux(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, uflux1, arraysum )   ! sum(uflux1(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, vflux1, arraysum )   ! sum(vflux1(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, uflux2, arraysum )   ! sum(uflux2(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, vflux2, arraysum )   ! sum(vflux2(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, uflux3, arraysum )   ! sum(uflux3(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, vflux3, arraysum )   ! sum(vflux3(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(uflx,dim=3), arraySum )   ! sum(uflx(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(vflx,dim=3), arraySum )   ! sum(vflx(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(bolusu,dim=3), arraySum )   ! sum(bolusu(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(bolusv,dim=3), arraySum ) !sum(bolusv(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum
c
      call GLOBALSUM( ogrid, sum(uav,dim=3), arraySum ) !sum(uav(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(vav,dim=3), arraySum ) !sum(vav(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(dpuav,dim=3), arraySum ) !sum(dpuav(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(dpvav,dim=3), arraySum ) !sum(dpvav(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(temav,dim=3), arraySum ) !sum(temav(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(salav,dim=3), arraySum ) !sum(salav(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(th3av,dim=3), arraySum ) !sum(th3av(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(dpav,dim=3), arraySum ) !sum(dpav(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, pbavav, arraysum )   ! sum(pbavav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sfhtav, arraysum )   ! sum(sfhtav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(uflxav,dim=3), arraySum ) !sum(uflxav(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(vflxav,dim=3), arraySum ) !sum(vflxav(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(ufxavp,dim=3), arraySum ) !sum(ufxavp(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(vfxavp,dim=3), arraySum ) !sum(vfxavp(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(diaflx,dim=3), arraySum ) !sum(diaflx(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, salflav, arraysum )   ! sum(salflav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, brineav, arraysum )   ! sum(brineav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, eminpav, arraysum )   ! sum(eminpav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, surflav, arraysum )   ! sum(surflav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, tauxav, arraysum )   ! sum(tauxav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, tauyav, arraysum )   ! sum(tauyav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(ufxcum,dim=3), arraySum ) !sum(ufxcum(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(vfxcum,dim=3), arraySum ) !sum(vfxcum(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(dpinit,dim=3), arraySum ) !sum(dpinit(:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, dpmxav, arraysum )   ! sum(dpmxav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, oiceav, arraysum )   ! sum(oiceav(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum
c

      call GLOBALSUM( ogrid, util1, arraysum )   ! sum(util1(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, util2, arraysum )   ! sum(util2(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, util3, arraysum )   ! sum(util3(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, util4, arraysum )   ! sum(util4(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum
c

      call GLOBALSUM( ogrid, scpx, arraysum )   ! sum(scpx(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scpy, arraysum )   ! sum(scpy(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scux, arraysum )   ! sum(scux(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scuy, arraysum )   ! sum(scuy(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scvx, arraysum )   ! sum(scvx(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scvy, arraysum )   ! sum(scvy(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scqx, arraysum )   ! sum(scqx(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scqy, arraysum )   ! sum(scqy(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scu2, arraysum )   ! sum(scu2(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scv2, arraysum )   ! sum(scv2(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scp2, arraysum )   ! sum(scp2(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scq2, arraysum )   ! sum(scq2(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scuxi, arraysum )   ! sum(scuxi(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scvyi, arraysum )   ! sum(scvyi(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scp2i, arraysum )   ! sum(scp2i(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, scq2i, arraysum )   ! sum(scq2i(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum
c

      call GLOBALSUM( ogrid, pgfx, arraysum )   ! sum(pgfx(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, pgfy, arraysum )   ! sum(pgfy(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, gradx, arraysum )   ! sum(gradx(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, grady, arraysum )   ! sum(grady(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, depthu, arraysum )   ! sum(depthu(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, depthv, arraysum )   ! sum(depthv(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, pvtrop, arraysum )   ! sum(pvtrop(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, depths, arraysum )   ! sum(depths(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, drag, arraysum )   ! sum(drag(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, glue, arraysum )   ! sum(glue(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, dampu, arraysum )   ! sum(dampu(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, dampv, arraysum )   ! sum(dampv(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum
c

      call GLOBALSUM( ogrid, uja, arraysum )   ! sum(uja(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, ujb, arraysum )   ! sum(ujb(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, via, arraysum )   !  sum(via(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, vib, arraysum )   ! sum(vib(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, pbot, arraysum )   ! sum(pbot(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, sum(sum(tracer,3),3), arraysum ) !sum(tracer(:,:,:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, tprime, arraysum )   ! sum(tprime(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, sum(sgain(:,:))

      call GLOBALSUM( ogrid, surflx, arraysum )   ! sum(surflx(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, salflx, arraysum )   ! sum(salflx(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, odhsi, arraysum )   ! sum(odhsi(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, dmfz, arraysum )   ! sum(dmfz(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum
c

      ! write(801,*) 'hycom.f ',__LINE__,sum(klist(:,:))
      call GLOBALSUM( ogrid, sum(klist(:,J_0:J_1)),iarraySum )! sum(klist(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, iarraySum
c

      call GLOBALSUM( ogrid, taux, arraysum )   ! sum(taux(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, tauy, arraysum )   ! sum(tauy(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, oemnp, arraysum )   ! sum(oemnp(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, oflxa2o, arraysum )   ! sum(oflxa2o(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, oice, arraysum )   ! sum(oice(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, ustar, arraysum )   ! sum(ustar(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, ustarb, arraysum )   ! sum(ustarb(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, osalt, arraysum )   ! sum(osalt(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum
c
      call GLOBALSUM( ogrid, freshw, arraysum )   ! sum(freshw(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      call GLOBALSUM( ogrid, diafor, arraysum )   ! sum(diafor(:,:))
      if(AM_I_ROOT()) write(801,*) 'hycom.f ',__LINE__, arraySum

      end subroutine hycom_arrays_checksum
c------------------------------------------------------------------
      subroutine scatter1_hycom_arrays
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: UNPACK_DATA
      USE HYCOM_ARRAYS_GLOB
      use hycom_arrays_glob_renamer
      USE KPRF_ARRAYS, ONLY: sswflx, sswflx_loc,
     &                       akpar,  akpar_loc

      implicit none

      call unpack_data( ogrid,  akpar,  akpar_loc )

      end subroutine scatter1_hycom_arrays
c------------------------------------------------------------------
