#include "rundeck_opts.h"
#ifdef SKIP_TRACERS_RAD
#undef TRACERS_ON
#endif
!@sum RAD_DRV contains drivers for the radiation related routines
!@ver  2009/05/11
!@cont init_RAD, RADIA
C**** semi-random cloud overlap (computed opt.d+diagn)
C**** to be used with R99E or later radiation  routines.  carbon/2
C****

      SUBROUTINE CALC_ZENITH_ANGLE
!@sum calculate zenith angle for current time step
!@auth Gavin Schmidt (from RADIA)
      USE CONSTANT, only : twopi
      USE MODEL_COM, only : itime,nday,dtsrc, calendar
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE RAD_COM, only : cosz1
      USE RAD_COSZ0, only : coszt
      USE TimeInterval_mod
      IMPLICIT NONE
      INTEGER JTIME
      REAL*8 ROT1,ROT2
      type (TimeInterval) :: sPerDay

      JTIME=MOD(ITIME,NDAY)
      ROT1=(TWOPI*JTIME)/NDAY
      sPerDay = calendar%getSecondsPerDay()
      ROT2=ROT1+TWOPI*DTsrc/real(sPerDay)
      CALL COSZT (ROT1,ROT2,COSZ1)

      END SUBROUTINE CALC_ZENITH_ANGLE

      SUBROUTINE init_RAD(istart)
!@sum  init_RAD initialises radiation code
!@auth Original Development Team
!@calls RADPAR:RCOMP1, ORBPAR
      USE FILEMANAGER
      use RunTimeControls_mod, only : tracers_minerals
      USE Dictionary_mod
      USE CONSTANT, only : grav,bysha,twopi,planet_name
      USE RESOLUTION, only : jm,lm,psf
      USE ATM_COM, only : t,pk,kradia,lm_req
      USE MODEL_COM, only : dtsrc,iyear1,modelEclock,master_yr
      USE MODEL_COM, only: orbit
      USE ATM_COM, only : pednl00
      USE DOMAIN_DECOMP_ATM, only : grid, write_parallel, am_i_root
     &     ,readt_parallel, getDomainBounds
#ifndef CUBED_SPHERE
      USE GEOM, only : lat_dg
#endif
      USE RADPAR, only : !rcomp1,writer,writet       ! routines
     &      PTLISO ,KTREND ,LMR=>NL, PLB, LS1_loc
     &     ,planck_tmin,planck_tmax
     &     ,transmission_corrections
     *     ,KCLDEM,KSIALB,KSOLAR, SHL, snoage_fac_max, KZSNOW
     *     ,KYEARS,KJDAYS,MADLUV, KYEARG,KJDAYG,MADGHG
     *     ,KYEARO,KJDAYO,MADO3M, KYEARA,KJDAYA,MADAER
     *     ,KYEARD,KJDAYD,MADDST, KYEARV,KJDAYV,MADVOL
     *     ,KYEARE,KJDAYE,MADEPS, KYEARR,KJDAYR
!g95     *     ,FSXAER,FTXAER    ! scaling (on/off) for default aerosols
     *     ,ITR,nraero_aod=>NTRACE ! turning on options for extra aerosols
     *     ,FS8OPX,FT8OPX, TRRDRY,KRHTRA,TRADEN,REFDRY
     *     ,rcomp1, writer, writet
     *     ,FSTASC,FTTASC
#ifdef ALTER_RADF_BY_LAT
     *     ,FS8OPX_orig,FT8OPX_orig
#endif
#ifdef TRACERS_SPECIAL_Shindell
      use photolysis, only: aer2,miedx2,nbfastj
#endif  /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_ON
      use rad_com, only: nraero_rf,nraero_seasalt,
     *                   nraero_koch,nraero_nitrate,nraero_dust,
     *                   nraero_OMA,nraero_AMP,nraero_TOMAS
#endif  /* TRACERS_ON */
      USE RAD_COM, only : rqt, s0x, co2x,n2ox,ch4x,cfc11x,cfc12x,xGHGx
     *     ,o2x,no2x,n2cx,yGHGx,so2x,CH4X_RADoverCHEM,snoage_def
     *     ,s0_yr,s0_day,ghg_yr,ghg_day,volc_yr,volc_day
     *     ,aero_yr,dust_yr,O3_yr
     *     ,H2ObyCH4,dH2O,h2ostratx,O3x,RHfix,CLDx,ref_mult,COSZ1
     *     ,obliq,eccn,omegt,obliq_def,eccn_def,omegt_def
     *     ,CC_cdncx,OD_cdncx,cdncl,pcdnc,vcdnc
     *     ,cloud_rad_forc,TAero_aod_diag,aer_rad_forc
     *     ,PLB0,shl0  ! saved to avoid OMP-copyin of input arrays
     *     ,albsn_yr,dALBsnX,nradfrc
     *     ,rad_interact_aer,clim_interact_chem,rad_forc_lev
     *     ,ntrix_aod,ntrix_rf,wttr
     *     ,variable_orb_par,orb_par_year_bp,orb_par,nrad
     *     ,radiationSetOrbit
     *     ,chl_from_obio,chl_from_seawifs
#ifdef TRACERS_ON
     *     ,njaero,nraero_aod_rsf,nraero_rf_rsf,tau_as,tau_cs,tau_dry
#ifdef CACHED_SUBDD
     *     ,abstau_as,abstau_cs,abstau_dry,swfrc,lwfrc
#endif  /* CACHED_SUBDD */
#endif  /* TRACERS_ON */
#ifdef ALTER_RADF_BY_LAT
     *     ,FULGAS_lat,FS8OPX_lat,FT8OPX_lat
#endif
      use RAD_COSZ0, only : cosz_init
      USE CLOUDS_COM, only : llow
      USE DIAG_COM, only : iwrite,jwrite,itwrite
#ifdef TRACERS_ON
      USE DIAG_COM, only : save3dAOD
      USE TRACER_COM, only: ntm
      USE TRACER_COM, only: n_BCIA, n_BCB, n_NO3p
      USE TRACER_COM, only: n_Clay, n_Silt1, n_Silt2, n_Silt3, n_Silt4,
     &     n_Silt5
      USE TRACER_COM, only: n_SO4, n_Seasalt1, n_Seasalt2
      USE TRACER_COM, only: n_OCB, n_OCIA, n_Isopp1a, n_SO4
      USE TRACER_COM, only: n_vbsAm2
      use RAD_COM, only: diag_fc
      use TRDIAG_COM, only: save_dry_aod
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only: n_ASO4, n_ANACL, n_AECOB, n_AECIL,
     &     n_AOCOB, n_AOCIL, n_ADUST
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      use OldTracer_mod, only: trpdens
      use trdust_mod, only : imDust, nSubClays, dryEffRadMinerals,
     &     subClayWeights
      use trdust_drv, only : calcSubClayWeights
      use tracer_com, only: ntm_clay, ntm_sil1, ntm_sil2, ntm_sil3,
     &     ntm_sil4, ntm_sil5, n_soilDust
      use rad_com, only: nr_soildust
#endif
#ifdef TRACERS_MINERALS
      use tracer_com, only: n_clayilli, n_claykaol, n_claysmec,
     &     n_claycalc, n_clayquar, n_clayfeld, n_clayhema, n_claygyps,
     &     n_clayilhe, n_claykahe, n_claysmhe, n_claycahe, n_clayquhe,
     &     n_clayfehe, n_claygyhe, n_sil1quar, n_sil1feld, n_sil1calc,
     &     n_sil1illi, n_sil1kaol, n_sil1smec, n_sil1hema, n_sil1gyps,
     &     n_sil1quhe, n_sil1fehe, n_sil1cahe, n_sil1gyhe, n_sil1ilhe,
     &     n_sil1kahe, n_sil1smhe, n_sil2quar, n_sil2feld, n_sil2calc,
     &     n_sil2hema, n_sil2gyps, n_sil2illi, n_sil2kaol, n_sil2smec,
     &     n_sil2quhe, n_sil2fehe, n_sil2cahe, n_sil2gyhe, n_sil2ilhe,
     &     n_sil2kahe, n_sil2smhe, n_sil3quar, n_sil3feld, n_sil3calc,
     &     n_sil3hema, n_sil3gyps, n_sil3illi, n_sil3kaol, n_sil3smec,
     &     n_sil3quhe, n_sil3fehe, n_sil3cahe, n_sil3gyhe, n_sil3ilhe,
     &     n_sil3kahe, n_sil3smhe, n_sil4quar, n_sil4feld, n_sil4calc,
     &     n_sil4hema, n_sil4gyps, n_sil4illi, n_sil4kaol, n_sil4smec,
     &     n_sil4quhe, n_sil4fehe, n_sil4cahe, n_sil4gyhe, n_sil4ilhe,
     &     n_sil4kahe, n_sil4smhe, n_sil5quar, n_sil5feld, n_sil5calc,
     &     n_sil5hema, n_sil5gyps, n_sil5illi, n_sil5kaol, n_sil5smec,
     &     n_sil5quhe, n_sil5fehe, n_sil5cahe, n_sil5gyhe, n_sil5ilhe,
     &     n_sil5kahe, n_sil5smhe
#endif
#ifdef TRACERS_AMP
      USE AERO_CONFIG, only: nmodes
      USE TRACER_COM, only:
     *     n_N_AKK_1 ,n_N_ACC_1 ,n_N_DD1_1 ,n_N_DS1_1 ,n_N_DD2_1,
     *     n_N_DS2_1, n_N_SSA_1, n_N_SSC_1, n_N_OCC_1, n_N_BC1_1,
     *     n_N_BC2_1 ,n_N_BC3_1,
     *     n_N_DBC_1, n_N_BOC_1, n_N_BCS_1, n_N_MXX_1
#endif
#ifdef TRACERS_TOMAS
      USE TOMAS_AEROSOL, only: icomp
#endif
      use AerParam_mod, only : aermix
#ifdef OLD_BCdalbsn
      use AerParam_mod, only: depoBC,depoBC_1990
#endif

      use AbstractOrbit_mod, only: AbstractOrbit
      ! begin section for radiation-only SCM
      use constant, only : gasc,tf,mair,mwat,pi,lhe,lhs,mb2kg,kg2mb,kapa
      use atm_com, only : q,p,pmid,pedn,pdsig,pek,ma,byma,ltropo
      use atm_com, only : aml00,byaml00,req_fac,kradia,lm_req
      use resolution, only : im, plbot,ls1=>ls1_nominal
      use resolution, only: mfix,mfrac
#ifndef STDHYB
      use resolution, only: mfixs,mtop
#endif
      use rad_com, only : modrd
      use radpar, only : u0gas,ulgas,set_gases_internally
      use radpar, only : set_aerosols_internally,
     &    sraext,srasct,sragcb,
     &    srdext,srdsct,srdgcb,
     &    srvext,srvsct,srvgcb,
     &    srbext,srbsct,srbgcb,
     &    traalk,trdalk,trvalk,trbalk
      use radpar, only: keepal,srbalb,srxalb,FSTOPX,FTTOPX
      use radpar, only: skip_AOD_in_rad
      use pario, only : par_open,par_close,read_data,read_dist_data
      use fluxes, only : atmsrf,asflx4,focean,fland,flice
      use fluxes, only : atmocn,atmice,atmgla,atmlnd
      use ghy_com, only : fearth
      use lakes_com, only : flake
      use seaice_com, only : si_atm
      use clouds_com, only : svlhx,svlat,rhsav
      !use clouds_com, only : lmid,lhi
      ! end section for radiation-only SCM
      IMPLICIT NONE

      integer, intent(in) :: istart

      INTEGER L,LR,n1,n,nn,iu2 ! LONR,LATR
      REAL*8 PLBx(LM+1),pyear
!@var NRFUN indices of unit numbers for radiation routines
      INTEGER NRFUN(14),IU,DONOTREAD
!@var RUNSTR names of files for radiation routines
      CHARACTER*5 :: RUNSTR(14) = (/"RADN1","RADN2","RADN3",
     *     "RADN4","RADN5","RADN6","RADN7","RADN8",
     *     "RADN9","RADNA","RADNB","RADNC","RADND",
     *     "RADNE"/)
!@var QBIN true if files for radiation input files are binary
      LOGICAL :: QBIN(14)=(/.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.
     *     ,.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)

#ifdef TRACERS_MINERALS
      real(kind=8) :: densclay( 4*ntm_clay ), denssil1( ntm_sil1 ),
     &     denssil2( ntm_sil2 ), denssil3( ntm_sil3 ), denssil4(
     &     ntm_sil4 ), denssil5( ntm_sil5 )
#endif

      character(len=300) :: out_line
      character*6 :: skip

      ! begin section for radiation-only SCM
      real*8 :: cosz_const, mvar
      character(len=6) :: gasnames(13)
      integer :: fid,igas
      real*8 :: szadeg,s0cosz,s0_tmp,cosz_tmp,tloc
      integer :: rad_scm_int
      logical :: rad_scm=.false.
      real*8 qsat ! external
      ! end section for radiation-only SCM

      INTEGER :: I,J
      INTEGER :: I_0,I_1,J_0,J_1
      integer :: I_0H, I_1H
      integer :: J_0H, J_1H

C**** sync radiation parameters from input
      call sync_param( "NRAD", NRAD ) !!
      if (is_set_param("variable_orb_par")) then
        call get_param( "variable_orb_par", variable_orb_par )
      else
        if (master_yr == 0) then
          variable_orb_par=1
        else
          variable_orb_par=0
        endif
      endif
      if (is_set_param("orb_par_year_bp")) then
        call get_param( "orb_par_year_bp", orb_par_year_bp )
      else
        if (master_yr == 0) then
          orb_par_year_bp=0
        else
          orb_par_year_bp=1950-master_yr
        endif
      endif
      call sync_param( "orb_par", orb_par, 3 )
      call sync_param( "S0X", S0X )
      call sync_param( "CO2X", CO2X )     ! fulgas(2)
      call sync_param( "O2X", O2X )       ! fulgas(4)
      call sync_param( "NO2X", NO2X )     ! fulgas(5)
      call sync_param( "N2OX", N2OX )     ! fulgas(6)
      call sync_param( "CH4X", CH4X )     ! fulgas(7)
      call sync_param( "CH4X_RADoverCHEM", CH4X_RADoverCHEM )
      call sync_param( "CFC11X", CFC11X ) ! fulgas(8)
      call sync_param( "CFC12X", CFC12X ) ! fulgas(9)
      call sync_param( "N2CX", N2CX )     ! fulgas(10)
      call sync_param( "XGHGX", XGHGX )   ! fulgas(11)
      call sync_param( "YGHGX", YGHGX )   ! fulgas(12)
      call sync_param( "SO2X", SO2X )     ! fulgas(13)
      call sync_param( "H2OstratX", H2OstratX ) ! fulgas(1)
      call sync_param( "O3X", O3X )       ! fulgas(3)
      call sync_param( "CLDX", CLDX )
      call sync_param( "H2ObyCH4", H2ObyCH4 )
      call get_param( "S0_yr", S0_yr, default=master_yr )
      if (is_set_param("S0_day")) then
        call get_param( "S0_day", S0_day )
      else
        if (s0_yr==0) s0_day=0 ! else use default value
      endif
      call get_param( "ghg_yr", ghg_yr, default=master_yr )
      if (is_set_param("ghg_day")) then
        call get_param( "ghg_day", ghg_day )
      else
        if (ghg_yr==0) ghg_day=0 ! else use default value
      endif
      call get_param( "volc_yr", volc_yr, default=master_yr )
      if (is_set_param("volc_day")) then
        call get_param( "volc_day", volc_day )
      else
        if (volc_yr==0) volc_day=0 ! else use default value
      endif
      call get_param( "aero_yr", aero_yr, default=master_yr )
      call get_param( "dust_yr", dust_yr, default=master_yr )
      call sync_param( "dALBsnX", dALBsnX )
      call get_param( "albsn_yr", albsn_yr, default=master_yr )
      call sync_param( "aermix", aermix , 13 )
      call sync_param( "REFdry", REFdry , 8 )
      call sync_param( "FS8OPX", FS8OPX , 8 )
      call sync_param( "FT8OPX", FT8OPX , 8 )
      call sync_param( "RHfix", RHfix )
      call sync_param( "CC_cdncx", CC_cdncx )
      call sync_param( "OD_cdncx", OD_cdncx )
      call get_param( "O3_yr", O3_yr, default=master_yr )
      if(planet_name.ne.'Earth') then
        PTLISO = .015d0*psf ! reasonable default
      endif
      call sync_param( "PTLISO", PTLISO )
      call sync_param( "KSOLAR", KSOLAR )
      call sync_param( "KSIALB", KSIALB )
      call sync_param( "KZSNOW", KZSNOW )
      call sync_param( "snoage_def", snoage_def )
      call sync_param( "snoage_fac_max", snoage_fac_max )
      call sync_param( "nradfrc", nradfrc )
      if(snoage_fac_max.lt.0. .or. snoage_fac_max.gt.1.) then
        write(out_line,*) 'set 0<snoage_fac_max<1, not',snoage_fac_max
        call write_parallel(trim(out_line),unit=6)
        call stop_model('init_RAD: snoage_fac_max out of range',255)
      end if
      call sync_param( "rad_interact_aer", rad_interact_aer )
      call sync_param( "clim_interact_chem", clim_interact_chem )
      call sync_param( "rad_forc_lev", rad_forc_lev )
      call sync_param( "cloud_rad_forc", cloud_rad_forc )
      call sync_param( "TAero_aod_diag", TAero_aod_diag )
      call sync_param( "aer_rad_forc", aer_rad_forc )
      call sync_param( "ref_mult", ref_mult )
#ifdef TRACERS_ON
      call sync_param( "save3dAOD", save3dAOD)
      CALL sync_param("save_dry_aod",save_dry_aod)
#endif
      REFdry = REFdry*ref_mult

      if(is_set_param('planck_tmin')) then
        call get_param('planck_tmin',planck_tmin)
      endif
      if(is_set_param('planck_tmax')) then
        call get_param('planck_tmax',planck_tmax)
      endif

      call getDomainBounds(grid,
     &     I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      call radiationSetOrbit(orbit)

      if(is_set_param('rad_scm')) then
        call get_param('rad_scm',rad_scm_int)
        rad_scm = im.eq.1 .and. jm.eq.1 .and. rad_scm_int==1
      endif
      if(rad_scm) then
        ! Some of the initializations in this block will be
        ! moved elsewhere once refactorings in the rest of
        ! the GCM are complete.

        modrd = 0 ! just in case
        i = 1
        j = 1

        if(file_exists('TOPO')) then
          fid = par_open(grid,'TOPO','read')
          call read_dist_data(grid,fid,'focean',focean)
          call read_dist_data(grid,fid,'fgice',flice)
          call par_close(grid,fid)
        else
          call get_param('focean',focean(1,1))
          call get_param('flice',flice(1,1))
          call get_param('rsi',si_atm%rsi(1,1))
        endif
        fland = 1d0-focean
        fearth = 1d0 - focean - flice
        flake = 0.

        pednl00(1:lm+1) = plbot ! default
        pednl00(lm+2:lm+lm_req) = req_fac(1:lm_req-1)*pednl00(lm+1)
        pednl00(lm+lm_req+1) = 0.

        if(file_exists('TEMP1D')) then
          fid = par_open(grid,'TEMP1D','read')
          call read_data(grid,fid,'t',t(i,j,:))
          call par_close(grid,fid)
        elseif(file_exists('AIC')) then
          fid = par_open(grid,'AIC','read')
          call read_dist_data(grid,fid,'t',t)
          call read_dist_data(grid,fid,'q',q)
          call read_dist_data(grid,fid,'p',p)  !  surface pressure (mb)
          call read_dist_data(grid,fid,'tsurf',atmsrf%gtempr)
!****     Rescaling: pednl00(1) = ps (= p)
          mvar = p(1,1)*mb2kg
#ifndef STDHYB
     &         - mfixs - mtop
#endif
          pednl00(lm+1) = plbot(lm+1)!mtop*kg2mb
          do l=lm,1,-1
            pednl00(l) = pednl00(l+1) + (mfix(l) + mvar*mfrac(l))*kg2mb
          enddo
          call par_close(grid,fid)
        elseif(is_set_param('temp1d')) then
          call get_param('temp1d',t(i,j,:),lm)
        endif
        do l=1,lm+lm_req
          aml00(l) = mb2kg*(pednl00(l)-pednl00(l+1))
          byaml00(l) = 1d0/aml00(l)
        enddo
        pedn(1:lm+1,i,j) = pednl00(1:lm+1)
        do l=1,lm
          pmid(l,i,j) = .5d0*(pedn(l,i,j)+pedn(l+1,i,j))
          pdsig(l,i,j) = (pedn(l,i,j)-pedn(l+1,i,j))
          ma(l,i,j) = pdsig(l,i,j)*mb2kg
          byma(l,i,j) = 1d0/ma(l,i,j)
        enddo
        pk = pmid**kapa
        pek = pedn**kapa
        p(i,j) = pedn(1,i,j)-pedn(ls1,i,j)

        do l=1,lm
          t(i,j,l) = t(i,j,l)/pk(l,i,j)
        enddo

        if(file_exists('GTEMPR')) then
          fid = par_open(grid,'GTEMPR','read')
          call read_data(grid,fid,'gtempr',atmsrf%gtempr)
          call par_close(grid,fid)
        elseif(is_set_param('gtempr')) then
          call get_param('gtempr',atmsrf%gtempr(1,1))
        endif

        call get_param('wsavg',atmsrf%wsavg(1,1),default=7d0)
        call get_param('bare_soil_wetness',
     &       atmlnd%bare_soil_wetness(1,1),default=1d0)
        call get_param('snow',atmsrf%snow(1,1),default=0d0)
#ifdef RESET_SURFACE_FRACTIONS_ON_SOUTH_POLE
        call get_param('vdata',vdata(1,1,:),size(vdata,3),
     &      default=(/0d0,0d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/))
#endif
        if(file_exists('SUN')) then
          fid = par_open(grid,'SUN','read')
          call read_data(grid,fid,'szadeg',szadeg)
          call read_data(grid,fid,'s0cosz',s0cosz)
          call par_close(grid,fid)
          cosz_tmp = cos(szadeg*pi/180d0)
          s0_tmp = s0cosz/cosz_tmp
          call sync_param('cosz',cosz_tmp)
          call sync_param('s0',s0_tmp)
        endif
        LTROPO = 1

        atmsrf%tsavg = atmsrf%gtempr
        atmocn%gtempr = atmsrf%gtempr
        atmice%gtempr = atmsrf%gtempr
        atmgla%gtempr = atmsrf%gtempr
        atmlnd%gtempr = atmsrf%gtempr
        do n=1,4
          asflx4(n)%gtempr = atmsrf%gtempr
        enddo
        atmgla%snow = atmsrf%snow

      endif

      if(is_set_param('srxalb')) then
        keepal = 1
        call get_param('srxalb',srxalb,size(srxalb))
        call get_param('srbalb',srbalb,size(srbalb))
      endif
      if(is_set_param('cosz')) then
        call get_param('cosz',cosz_const)
        call cosz_init(cosz_const=cosz_const)
      else
        call cosz_init
      endif

      call sync_param("chl_from_obio", chl_from_obio)
      call sync_param("chl_from_seawifs", chl_from_seawifs)
      if ( chl_from_obio>0 .and. chl_from_seawifs>0 ) then
        call stop_model("Make your mind which chl to use",255)
      endif

      if(istart==2) then ! replace with cold vs warm start logic
C**** SET RADIATION EQUILIBRIUM TEMPERATURES FROM LAYER LM TEMPERATURE
        DO J=J_0,J_1
        DO I=I_0,I_1
          RQT(:,I,J)=T(I,J,LM)*PK(LM,I,J)
        ENDDO
        ENDDO
      endif


C****
C**** SET THE CONTROL PARAMETERS FOR THE RADIATION (need mean pressures)
C****
      LMR=LM+LM_REQ
      PLB(1:LMR+1)=PEDNL00(1:LMR+1)
      DO L=1,LM
        PLBx(L)=PLB(L)           ! needed for CH4 prod. H2O
      END DO
      PLBx(LM+1)=0.
      DO LR=LM+1,LMR
        PLB0(LR-LM) = PLB(LR+1)
      END DO
      cdncl=0 ; call reterp(vcdnc,pcdnc,7, cdncl,plb,llow+2)

      KTREND=1   !  GHgas trends are determined by input file
!note KTREND=0 is a possible but virtually obsolete option
C****
C             Model Add-on Data of Extended Climatology Enable Parameter
C     MADO3M  = -1   Reads                      Ozone data the GCM way
C     MADAER  =  1   Reads   Tropospheric Aerosol climatology 1850-2050
C     MADAER  =  3   uses Koch,Bauer 2008 aerosol climatology 1890-2000
C     MADDST  =  1   Reads   Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  1   Reads   Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  1   Reads   Epsilon cloud heterogeniety data     RFILE8
C     MADLUV  =  1   Reads   Lean''s SolarUV 1882-1998 variability RFILE9
C**** Radiative forcings are either constant = obs.value at given yr/day
C****    or time dependent (year=0); if day=0 an annual cycle is used
C****                                         even if the year is fixed
      KYEARS=s0_yr   ; KJDAYS=s0_day ;  MADLUV=1   ! solar 'constant'
      KYEARG=ghg_yr  ; KJDAYG=ghg_day              ! well-mixed GHGases
#ifndef ALTER_RADF_BY_LAT
      if(ghg_yr.gt.0)  MADGHG=0                    ! skip GHG-updating
#endif
      KYEARO=O3_yr   ; KJDAYO=0 ;       MADO3M=-1  ! ozone (ann.cycle)
      if(KYEARO.gt.0) KYEARO=-KYEARO              ! use ONLY KYEARO-data
      KYEARA=Aero_yr ; KJDAYA=0 ! MADAER=1 or 3, trop.aeros (ann.cycle)
      if(KYEARA.gt.0) KYEARA=-KYEARA              ! use ONLY KYEARA-data
      if(file_exists('TAero_SSA')) MADAER=3   ! one of the TAero_XXX set
      KYEARD=Dust_yr
      if(KYEARD.gt.0) KYEARD=-KYEARD              ! use ONLY KYEARD-data
      KYEARV=Volc_yr ; KJDAYV=Volc_day
      if(file_exists('RADN7')) MADVOL=1   ! Volc. Aerosols
      call sync_param( "MADVOL", MADVOL )
!***  KYEARV=0 : use current year
!***  KYEARV<0 : use long term mean stratospheric aerosols (use -1)
!     Hack: KYEARV= -2000 and -2010 were used for 2 specific runs that
!           ended in 2100 and repeated some 20th century volcanos
!***  KYEARV=-2000: use volcanos from 100 yrs ago after 2000
!***  KYEARV=-2010: repeat 2nd half, then first half of 20th century
      if(KYEARV.le.-2000) KYEARV=0   ! use current year (before 2000)
C**** NO time history (yet), except for ann.cycle, for forcings below;
C****  if KJDAY?=day0 (1->365), data from that day are used all year
      KYEARE=0       ; KJDAYE=0 ;
      KYEARR=0       ; KJDAYR=0           ! surf.reflectance (ann.cycle)
      KCLDEM=1                  ! 0:old 1:new LW cloud scattering scheme

      if(file_exists('DUSTaer'))   MADDST=1   ! Desert dust
      if(file_exists('RADN8'))     MADEPS=1   ! cloud Epsln - KCLDEP
      transmission_corrections = file_exists('RADN4')

C**** Aerosols:
C**** Currently there are five different default aerosol controls
C****   1:total 2:background+tracer 3:Climatology 4:dust 5:volcanic
C**** By adjusting FSXAER,FTXAER you can remove the default
C**** aerosols and replace them with your version if required
C**** (through TRACER in RADIA).
C**** FSXAER is for the shortwave,    FTXAER is for the longwave effects
caer  FSXAER = (/ 1.,1.,1.,1.,1. /) ; FTXAER = (/ 1.,1.,1.,1.,1. /)

C**** climatology aerosols are grouped into 6 types from 13 sources:
C****  Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
C****  ---------------------------------  ------------------  --------
C****   1    2    3    4    5    6    7    8    9   10   11   12   13
C****  SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
C**** using the following default scaling/tuning factors  AERMIX(1-13)
C****  1.0, 1.0, .26, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9
C**** The 8 groups are (adding dust and volcanic aerosols as 7. and 8.)
C**** 1. Sulfates (industr and natural), 2. Sea Salt, 3. Nitrates
C**** 4. Organic Carbons, 5. industr Black Carbons(BC), 6. Biomass BC
C**** 7. Dust aerosols, 8. Volcanic aerosols
C**** use FS8OPX and FT8OPX to enhance the optical effect; defaults:
caer  FS8OPX = (/1., 1., 1., 1., 2., 2.,    1.   ,   1./)     solar
caer  FT8OPX = (/1., 1., 1., 1., 1., 1.,    1.3d0,   1./)     thermal
!!!!! Note: FS|T8OPX(7-8) makes FS|TXAER(4-5) redundant.
C**** Particle sizes of the first 4 groups have RelHum dependence

C**** To add up to 8 further aerosols:
C****  1) set nraero_aod to the number of extra aerosol fields
C****  2) ITR defines which set of Mie parameters get used, choose
C****     from the following:
C****     1 SO4,  2 seasalt, 3 nitrate, 4 OCX organic carbons
C****     5 BCI,  6 BCB,     7 dust,    8 H2SO4 volc
C****  2b) set up the indexing array ntrix_aod to map the RADIATION tracers
C****      to the main model tracers
C****  2c) set up the weighting array WTTR to weight main model tracers,
C****      if needed (default value is 1).
C****
C****  3) Use FSTOPX/FTTOPX(1:nraero_aod) to scale them in RADIA
C****  4) Set TRRDRY to dry radius
C****  5) Set KRHTRA=1 if aerosol has RH dependence, 0 if not
C**** Note: whereas FSXAER/FTXAER are global (shared), FSTOPX/FTTOPX
C****       have to be reset for each grid box to allow for the way it
C****       is used in RADIA (TRACERS_AEROSOLS_Koch)
caer   nraero_aod = 0
caer   ITR = (/ 0,0,0,0, 0,0,0,0 /)
caer   TRRDRY=(/ .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0/)
caer   KRHTRA=(/1,1,1,1,1,1,1,1/)

#if (defined TRACERS_ON)

#if defined(TRACERS_AMP)
      nraero_AMP=nmodes
      IF (diag_fc==2) THEN
        nraero_rf=nraero_rf+nraero_AMP
      ELSE IF (diag_fc==1) THEN
        IF (nraero_AMP .gt. 0) nraero_rf=nraero_rf+1
      ENDIF
#elif defined(TRACERS_TOMAS)
!TOMAS does not include NO3 AND VOL, which use its default radiation.
#ifndef TRACERS_NITRATE
      nraero_TOMAS=icomp-2
#else
      nraero_TOMAS=icomp-1
#endif
      IF (diag_fc==2) THEN
        nraero_rf=nraero_rf+nraero_TOMAS
      ELSE IF (diag_fc==1) THEN
        IF (nraero_TOMAS .gt. 0) nraero_rf=nraero_rf+1
      ENDIF
#else
      nraero_OMA=nraero_seasalt+nraero_koch+nraero_nitrate+nraero_dust
      IF (diag_fc==2) THEN
        nraero_rf=nraero_rf+nraero_OMA
      ELSE IF (diag_fc==1) THEN
        IF (nraero_OMA .gt. 0) nraero_rf=nraero_rf+1
      ENDIF
#endif

      nraero_aod=nraero_OMA+nraero_AMP+nraero_TOMAS

      if (nraero_aod_rsf>0) then
        if (nraero_aod_rsf /= nraero_aod) then
          call stop_model('nraero_aod_rsf /= nraero_aod',255)
        endif
      endif

      if (nraero_rf_rsf>0) then
        if (nraero_rf_rsf /= nraero_rf) then
          call stop_model('nraero_rf_rsf /= nraero_rf',255)
        endif
      endif

      if (nraero_aod>0) then
        allocate(ntrix_aod(nraero_aod)) ; ntrix_aod=0
        if (nraero_rf>0) allocate(ntrix_rf(nraero_rf)) ; ntrix_rf=0
        allocate(wttr(nraero_aod))  ; wttr=1.

        if (.not.allocated(tau_as)) then
          allocate(tau_as(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod))
          allocate(tau_cs(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod))
          tau_as = 0.d0
          tau_cs = 0.d0
          if (save_dry_aod>0) then
            allocate(tau_dry(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod))
            tau_dry = 0.d0
          endif
#ifdef CACHED_SUBDD
          allocate(abstau_as(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod))
          allocate(abstau_cs(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod))
          abstau_as = 0.d0
          abstau_cs = 0.d0
          if (save_dry_aod>0) then
            allocate(abstau_dry(I_0H:I_1H,J_0H:J_1H,lm,nraero_aod))
            abstau_dry = 0.d0
          endif
          if (nraero_rf>0) then
            allocate(swfrc(I_0H:I_1H,J_0H:J_1H,nraero_rf))
            allocate(lwfrc(I_0H:I_1H,J_0H:J_1H,nraero_rf))
            swfrc = 0.d0
            lwfrc = 0.d0
          endif
#endif  /* CACHED_SUBDD */
        endif
      endif
#ifdef TRACERS_SPECIAL_Shindell
#if (! defined(TRACERS_AMP)) && (! defined(TRACERS_TOMAS))
      njaero=nraero_aod+2
#else
      njaero=2
#endif
      allocate(miedx2(nbfastj,njaero))
      allocate(aer2(nbfastj,njaero))
#endif  /* TRACERS_SPECIAL_Shindell */

!=======================================================================
! Define indices to map model aerosol tracer arrays to radiation arrays
! and other radiation-related aerosol properties
!=======================================================================
      n=0
!-----------------------------------------------------------------------
#ifdef TRACERS_AEROSOLS_SEASALT
      if (nraero_seasalt > 0) then
        if (rad_interact_aer > 0) then
          FS8OPX(2)=0.d0
          FT8OPX(2)=0.d0
        end if
        ntrix_aod(n+1:n+nraero_seasalt)=(/n_seasalt1,n_seasalt2/)
        trrdry(n+1:n+nraero_seasalt)=(/0.44d0,1.7d0/)
        itr(n+1:n+nraero_seasalt)=(/2,2/)
      endif
      n=n+nraero_seasalt
#endif  /* TRACERS_AEROSOLS_SEASALT */
!-----------------------------------------------------------------------
#ifdef TRACERS_AEROSOLS_Koch
      if (nraero_koch > 0) then
        if (rad_interact_aer > 0) then  ! if BC''s sol.effect are doubled:
          FS8OPX(1)=0.d0
          FT8OPX(1)=0.d0
#ifndef SULF_ONLY_AEROSOLS
          FS8OPX(4:6)=0.d0
          FT8OPX(4:6)=0.d0
#endif
        end if
        ntrix_aod(n+1)=n_SO4
        trrdry(n+1)=0.15d0
        itr(n+1)=1
#ifndef SULF_ONLY_AEROSOLS
        ntrix_aod(n+2:n+nraero_koch)=(/
#ifdef TRACERS_AEROSOLS_VBS
     &                             n_vbsAm2
#else
     &                             n_OCIA,n_OCB
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_SOA
     &                            ,n_isopp1a
#endif  /* TRACERS_AEROSOLS_SOA */
     &                            ,n_BCIA,n_BCB/)
        trrdry(n+2:n+nraero_koch)=(/
#ifdef TRACERS_AEROSOLS_VBS
     &                             0.2d0
#else
     &                             0.2d0,0.2d0
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_SOA
     &                             ,0.2d0
#endif  /* TRACERS_AEROSOLS_SOA */
     &                             ,0.08d0,0.08d0/)
        itr(n+2:n+nraero_koch)=(/
#ifdef TRACERS_AEROSOLS_VBS
     &                             4
#else
     &                             4,4
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_SOA
     &                          ,4
#endif  /* TRACERS_AEROSOLS_SOA */
     &                          ,5,6/)
        krhtra(n+2:n+nraero_koch)=(/
#ifdef TRACERS_AEROSOLS_VBS
     &                             1
#else
     &                             1,1
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_SOA
     &                             ,1
#endif  /* TRACERS_AEROSOLS_SOA */
     &                             ,0,0/)
        fstasc(n+2:n+nraero_koch)=(/
#ifdef TRACERS_AEROSOLS_VBS
     &                             1.d0
#else
     &                             1.d0,1.d0
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_SOA
     &                             ,1.d0
#endif  /* TRACERS_AEROSOLS_SOA */
! augment BC by 50%
     &                             ,1.5d0,1.5d0/)
#endif  /* SULF_ONLY_AEROSOLS */
      endif
      n=n+nraero_koch
#endif  /* TRACERS_AEROSOLS_Koch */
!-----------------------------------------------------------------------
#ifdef TRACERS_NITRATE
      if (nraero_nitrate > 0) then
#ifdef SULF_ONLY_AEROSOLS
        call stop_model('SULF_ONLY_AEROSOLS and TRACERS_NITRATE on',255)
#endif /* OFF: SULF_ONLY_AEROSOLS */
        if (rad_interact_aer > 0) then ! turn off default nitrate
          FS8OPX(3)=0.d0
          FT8OPX(3)=0.d0
        endif
        ntrix_aod(n+1:n+nraero_nitrate)=(/n_NO3p/)
        trrdry(n+1:n+nraero_nitrate)=(/0.15d0/)
        itr(n+1:n+nraero_nitrate) = (/3/)
      endif
      n=n+nraero_nitrate
#endif  /* TRACERS_NITRATE */
!-----------------------------------------------------------------------
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      if (nraero_dust > 0) then
        if (rad_interact_aer > 0) then ! turn off default dust
          FS8OPX(7)=0.d0
          FT8OPX(7)=0.d0
        end if
        nr_soildust=n+1

#ifdef TRACERS_MINERALS

! Adjust if number of dust tracers changes.
        ntrix_aod(n+1:n+nraero_dust)=(/
     &     (n_clayilli, i = 1,nSubClays), (n_claykaol, i = 1,nSubClays),
     &     (n_claysmec, i = 1,nSubClays), (n_claycalc, i = 1,nSubClays),
     &     (n_clayquar, i = 1,nSubClays), (n_clayfeld, i = 1,nSubClays),
     &     (n_clayhema, i = 1,nSubClays), (n_claygyps, i = 1,nSubClays),
     &     (n_clayilhe, i = 1,nSubClays), (n_claykahe, i = 1,nSubClays),
     &     (n_claysmhe, i = 1,nSubClays), (n_claycahe, i = 1,nSubClays),
     &     (n_clayquhe, i = 1,nSubClays), (n_clayfehe, i = 1,nSubClays),
     &     (n_claygyhe, i = 1,nSubClays), n_sil1illi, n_sil1kaol,
     &     n_sil1smec, n_sil1calc, n_sil1quar, n_sil1feld, n_sil1hema,
     &     n_sil1gyps, n_sil1ilhe, n_sil1kahe, n_sil1smhe, n_sil1cahe,
     &     n_sil1quhe, n_sil1fehe, n_sil1gyhe, n_sil2illi, n_sil2kaol,
     &     n_sil2smec, n_sil2calc, n_sil2quar, n_sil2feld, n_sil2hema,
     &     n_sil2gyps, n_sil2ilhe, n_sil2kahe, n_sil2smhe, n_sil2cahe,
     &     n_sil2quhe, n_sil2fehe, n_sil2gyhe, n_sil3illi, n_sil3kaol,
     &     n_sil3smec, n_sil3calc, n_sil3quar, n_sil3feld, n_sil3hema,
     &     n_sil3gyps, n_sil3ilhe, n_sil3kahe, n_sil3smhe, n_sil3cahe,
     &     n_sil3quhe, n_sil3fehe, n_sil3gyhe
#ifdef TRACERS_DUST_Silt4
     &     , n_sil4illi, n_sil4kaol, n_sil4smec, n_sil4calc, n_sil4quar
     &     , n_sil4feld, n_sil4hema, n_sil4gyps, n_sil4ilhe, n_sil4kahe
     &     , n_sil4smhe, n_sil4cahe, n_sil4quhe, n_sil4fehe, n_sil4gyhe
#endif  /* TRACERS_DUST_Silt4 */
#ifdef TRACERS_DUST_Silt5
     &     , n_sil5illi, n_sil5kaol, n_sil5smec, n_sil5calc, n_sil5quar
     &     , n_sil5feld, n_sil5hema, n_sil5gyps, n_sil5ilhe, n_sil5kahe
     &     , n_sil5smhe, n_sil5cahe, n_sil5quhe, n_sil5fehe, n_sil5gyhe
#endif  /* TRACERS_DUST_Silt5 */
     &     /)
        trrdry(n+1:n+nraero_dust)=
     &      (/(dryEffRadMinerals(1:nSubClays), i=1,ntm_clay),
     &        (dryEffRadMinerals(5), i=1,ntm_sil1),
     &        (dryEffRadMinerals(6), i=1,ntm_sil2),
     &        (dryEffRadMinerals(7), i=1,ntm_sil3)
#ifdef TRACERS_DUST_Silt4
     &       ,(dryEffRadMinerals(8), i=1,ntm_sil4)
#endif  /* TRACERS_DUST_Silt4 */
#ifdef TRACERS_DUST_Silt5
     &       ,(dryEffRadMinerals(9), i=1,ntm_sil5)
#endif  /* TRACERS_DUST_Silt5 */
     &       /)

        if ( tracers_minerals ) call calcSubClayWeights

        wttr( n+1:n+nraero_dust ) = (/ ( ( subClayWeights( i, j ), i=1
     &       ,nSubClays ), j=1,ntm_clay ), ( 1.d0, i=1,ntm_sil1+ntm_sil2
     &       +ntm_sil3+ntm_sil4+ntm_sil5 ) /)

        densclay=(/(trpdens(n_clayilli), i=1,nSubClays),
     &             (trpdens(n_claykaol), i=1,nSubClays),
     &             (trpdens(n_claysmec), i=1,nSubClays),
     &             (trpdens(n_claycalc), i=1,nSubClays),
     &             (trpdens(n_clayquar), i=1,nSubClays),
     &             (trpdens(n_clayfeld), i=1,nSubClays),
     &             (trpdens(n_clayhema), i=1,nSubClays),
     &             (trpdens(n_claygyps), i=1,nSubClays),
     &             (trpdens(n_clayilhe), i=1,nSubClays),
     &             (trpdens(n_claykahe), i=1,nSubClays),
     &             (trpdens(n_claysmhe), i=1,nSubClays),
     &             (trpdens(n_claycahe), i=1,nSubClays),
     &             (trpdens(n_clayquhe), i=1,nSubClays),
     &             (trpdens(n_clayfehe), i=1,nSubClays),
     &             (trpdens(n_claygyhe), i=1,nSubClays)/)
        denssil1=(/trpdens(n_sil1illi),trpdens(n_sil1kaol),
     &             trpdens(n_sil1smec),trpdens(n_sil1calc),
     &             trpdens(n_sil1quar),trpdens(n_sil1feld),
     &             trpdens(n_sil1hema),trpdens(n_sil1gyps),
     &             trpdens(n_sil1ilhe),trpdens(n_sil1kahe),
     &             trpdens(n_sil1smhe),trpdens(n_sil1cahe),
     &             trpdens(n_sil1quhe),trpdens(n_sil1fehe),
     &             trpdens(n_sil1gyhe)/)
        denssil2=denssil1
        denssil3=denssil1
#ifdef TRACERS_DUST_Silt4
        denssil4=denssil1
#endif  /* TRACERS_DUST_Silt4 */
#ifdef TRACERS_DUST_Silt5
        denssil5=denssil1
#endif  /* TRACERS_DUST_Silt5 */
        traden(n+1:n+nraero_dust)=(/densclay(:),denssil1(:),denssil2(:)
     &     ,denssil3(:)
#ifdef TRACERS_DUST_Silt4
     &     ,denssil4(:)
#endif  /* TRACERS_DUST_Silt4 */
#ifdef TRACERS_DUST_Silt5
     &     ,denssil5(:)
#endif  /* TRACERS_DUST_Silt5 */
     &       /) * 1d-3 ! Convert from kg/m^3 to g/cm^3

#else  /* not TRACERS_MINERALS */

        ntrix_aod(n+1:n+nraero_dust)=(/ ( n_clay,i = 1,nSubClays ),
     &       n_silt1, n_silt2, n_silt3
#ifdef TRACERS_DUST_Silt4
     &                            ,n_silt4
#endif  /* TRACERS_DUST_Silt4 */
#ifdef TRACERS_DUST_Silt5
     &                            ,n_silt5
#endif  /* TRACERS_DUST_Silt5 */
     &                            /)
        trrdry(n+1:n+nraero_dust)=
     &      (/(dryEffRadMinerals(1:nSubClays), i=1,ntm_clay),
     &        (dryEffRadMinerals(5), i=1,ntm_sil1),
     &        (dryEffRadMinerals(6), i=1,ntm_sil2),
     &        (dryEffRadMinerals(7), i=1,ntm_sil3)
#ifdef TRACERS_DUST_Silt4
     &       ,(dryEffRadMinerals(8), i=1,ntm_sil4)
#endif  /* TRACERS_DUST_Silt4 */
#ifdef TRACERS_DUST_Silt5
     &       ,(dryEffRadMinerals(9), i=1,ntm_sil5)
#endif  /* TRACERS_DUST_Silt5 */
     &                             /)

        if ( imDust >= 4 ) call calcSubClayWeights

        wttr( n+1:n+nraero_dust ) = (/ ( ( subClayWeights( i , j ), i =
     &       1,nSubClays ), j=1,ntm_clay ), ( 1.d0, i=1,ntm_sil1
     &       +ntm_sil2+ntm_sil3+ntm_sil4+ntm_sil5 ) /)

! Particle density of dust
        traden( n+1:n+nraero_dust )=(/ ( trpdens( n_clay ), i = 1
     &       ,nSubClays ), trpdens( n_silt1 ), trpdens( n_silt2 ),
     &       trpdens( n_silt3 )
#ifdef TRACERS_DUST_Silt4
     &                             ,trpdens( n_silt4 )
#endif  /* TRACERS_DUST_Silt4 */
#ifdef TRACERS_DUST_Silt5
     &                             ,trpdens( n_silt5 )
#endif  /* TRACERS_DUST_Silt5 */
     &                             /) * 1d-3 ! Convert from kg/m^3 to g/cm^3

#endif  /* TRACERS_MINERALS */

        itr(n+1:n+nraero_dust) = 7 ! all dust cases, outside ifdefs
        krhtra(n+1:n+nraero_dust) = 0 ! no deliq for dust or minerals
        fttasc(n+1:n+nraero_dust)=1.3d0 ! increase dust AOD by 1.3 in LW
      endif
      n=n+nraero_dust
#endif  /* (defined TRACERS_DUST) || (defined TRACERS_MINERALS) */
!-----------------------------------------------------------------------
!define ntrix_rf, based on the OMA tracers above
      if (n>0) then
        if (diag_fc==2) then
          ntrix_rf(1:nraero_OMA)=ntrix_aod(1:nraero_OMA)
        else if (diag_fc==1) then
          ntrix_rf(1)=ntrix_aod(1)
        endif
      endif
!-----------------------------------------------------------------------
#if (defined TRACERS_AMP) || (defined TRACERS_AMP_M1)
      if (nraero_AMP > 0) then
        if (rad_interact_aer > 0) then
          FS8OPX(1:7)=0.d0
          FT8OPX(1:7)=0.d0
        endif
        ntrix_aod(n+1:n+nraero_AMP)=
     &     (/n_N_AKK_1, n_N_ACC_1, n_N_DD1_1, n_N_DS1_1, n_N_DD2_1,
     &       n_N_DS2_1, n_N_SSA_1, n_N_SSC_1, n_N_OCC_1, n_N_BC1_1,
     &       n_N_BC2_1 ,n_N_BC3_1, n_N_DBC_1, n_N_BOC_1, n_N_BCS_1,
     &       n_N_MXX_1/)
        if (diag_fc==2) then
          ntrix_rf(n+1:n+nraero_AMP)=ntrix_aod(n+1:n+nraero_AMP)
        else if (diag_fc==1) then
          ntrix_rf(n+1)=ntrix_aod(n+1)
        endif
      endif
      n=n+nraero_AMP
#endif  /* (defined TRACERS_AMP) || (defined TRACERS_AMP_M1) */
!-----------------------------------------------------------------------
#ifdef TRACERS_TOMAS
      if (nraero_TOMAS > 0) then
        if (rad_interact_aer > 0) then
          FS8OPX(1:2)=0.d0
          FS8OPX(4:7)=0.d0
          FT8OPX(1:2)=0.d0
          FT8OPX(4:7)=0.d0
#ifdef TRACERS_NITRATE
          FS8OPX(3)=0.d0
          FT8OPX(3)=0.d0
#endif  /* TRACERS_NITRATE */
        endif

        ntrix_aod(n+1:n+nraero_TOMAS)=
     &     (/n_ASO4(1), n_ANACL(1), n_AECOB(1), n_AECIL(1),
     &       n_AOCOB(1), n_AOCIL(1), n_ADUST(1)/)
        itr(n+1:n+nraero_TOMAS) = (/1,2,6,5,4,4,7/)
        krhtra(n+1:n+nraero_TOMAS)=0
! ANUM(1) for internal-mixing case. Others(ncomp-1) for external-mixing case.
        if (diag_fc==2) then
          ntrix_rf(n+1:n+nraero_TOMAS)=ntrix_aod(n+1:n+nraero_TOMAS)
        else if (diag_fc==1) then
          ntrix_rf(n+1)=ntrix_aod(n+1)
        endif
      endif
      n=n+nraero_TOMAS
#endif
!=======================================================================
!=======================================================================
#endif  /* TRACERS_ON */

! set default FSTOPX and FTTOPX values
      if (rad_interact_aer > 0) then
        FSTOPX(:)=1.d0 ; FTTOPX(:)=1.d0
      else
        FSTOPX(:)=0.d0 ; FTTOPX(:)=0.d0
      endif
      skip_AOD_in_rad=rad_interact_aer>0

      if (ktrend.ne.0) then
C****   Read in time history of well-mixed greenhouse gases
        call openunit('GHG',iu,.false.,.true.)
        call ghghst(iu)
        call closeunit(iu)
        if(file_exists('dH2O').and.H2ObyCH4.ne.0..and.Kradia.le.0) then
C****     Read in dH2O: H2O prod.rate in kg/m^2 per day and ppm_CH4
          call openunit('dH2O',iu,.false.,.true.)
#if defined(CUBED_SPHERE)
          call read_qma(iu,plbx)
#else
          call getqma(iu,lat_dg,plbx,dh2o,lm,jm)
#endif
          call closeunit(iu)
        else
          H2ObyCH4 = 0.
        end if
      end if
#ifdef OLD_BCdalbsn
      if(dalbsnX.ne.0.) then
        call updBCd(1990) ; depoBC_1990 = depoBC
      endif
#endif
C**** set up unit numbers for 14 more radiation input files
      donotread = -9999
      nrfun(:) = 0 ! green light
      nrfun(12:13) = donotread ! not used in GCM
      nrfun(10:11) = donotread ! obsolete O3 data
      nrfun(6)     = donotread ! dust read externally now
      if(.not.transmission_corrections) nrfun(4) = donotread
      if(madvol == 0) nrfun(7) = donotread
      if(madeps == 0) nrfun(8) = donotread
!      if(ksolar < 0)  nrfun(9) = donotread
      nrfun(9) = donotread     ! open/read RADN9 inside RCOMP1
      DO IU=1,14
        if(nrfun(iu) == donotread) cycle
        call openunit(RUNSTR(IU),NRFUN(IU),QBIN(IU),.true.)
      END DO

      LS1_loc=1  ! default
C***********************************************************************
C     Main Radiative Initializations
C     ------------------------------------------------------------------
      CALL RCOMP1 (NRFUN)
      if (am_i_root()) CALL WRITER(6,0)  ! print rad. control parameters
C***********************************************************************
      DO IU=1,14
        if(nrfun(iu) == donotread) cycle
        call closeunit(NRFUN(IU))
      END DO
C**** Save initial (currently permanent and global) Q in rad.layers
      do LR=1,LM_REQ
        shl0(LR) = shl(LM+LR)
      end do
      write(out_line,*) 'spec.hum in rad.equ.layers:',shl0
      call write_parallel(trim(out_line),unit=6)

#ifdef ALTER_RADF_BY_LAT
C**** Save initial rad forcing alterations:
      FS8OPX_orig(:)=FS8OPX(:); FT8OPX_orig(:)=FT8OPX(:) ! aerosols

C**** Read in the factors used for alterations:
      call openunit('ALT_GHG_LAT',iu2,.false.,.true.)
      read(iu2,*) ! skip first line
      do n=1,46
        read(iu2,'(a6,13D8.3)') skip,(FULGAS_lat(nn,n),nn=1,13)
      enddo
      call closeunit(iu2)
      call openunit('ALT_AER_LAT',iu2,.false.,.true.)
      read(iu2,*) ! skip first line
      do n=1,46
        read(iu2,'(a6,8D8.3)') skip,(FS8OPX_lat(nn,n),nn=1,8)
      enddo
      read(iu2,*) ! skip first line
      do n=1,46
        read(iu2,'(a6,8D8.3)') skip,(FT8OPX_lat(nn,n),nn=1,8)
      enddo
      call closeunit(iu2)
#endif

c transplanted from main().  needs reviving
c      USE RAD_COM, only : dimrad_sv
c      CHARACTER aDATE*14
c      if (Kradia.ne.0 .and. Kradia<10) then
c        write(aDATE(1:7),'(a3,I4.4)') aMON(1:3),Jyear
c        if (Kradia.gt.0) aDATE(4:7)='    '
c        call openunit(trim('RAD'//aDATE(1:7)),iu_RAD,.true.,.false.)
c        if (Kradia.lt.0) call io_POS(iu_RAD,Itime-1,2*dimrad_sv,Nrad)
c      end if

      if(rad_scm) then
        if(file_exists('GASES')) then
C     GAS NUMBER    1         2    3      4    5         6           7
C                 H2O       CO2   O3     O2  NO2       N2O         CH4
C     GAS NUMBER    8         9   10        11          12          13
C              CCL3F1    CCL2F2   N2     CFC-Y       CFC-Z         SO2

          gasnames =
     &         ['h2o   ','co2   ','o3    ','o2    ','no2   ',
     &          'n2o   ','ch4   ','cfc11 ','cfc12 ','n2    ',
     &          'cfc-y ','cfc-z ','so2   ']
          set_gases_internally = .false.
          u0gas = 0.
          fid = par_open(grid,'GASES','read')
          do igas=1,size(gasnames)
            call read_data(grid,fid,trim(gasnames(igas)),
     &           u0gas(1:lm,igas))
            u0gas(lm+1:,igas) = u0gas(lm,igas) ! fill lm+1:lm+lm_req
            if(trim(gasnames(igas)).eq.'h2o') then
              q(1,1,:) = u0gas(1:lm,igas)*(mwat/mair) ! vol. ratio -> sp. hum.
            endif
            u0gas(:,igas) = aml00*u0gas(:,igas)   ! vol. ratio -> cm-atm
     &           *((1d5/mair)*(gasc*tf/101325d0))
          enddo
          call par_close(grid,fid)

         !fulgas = 1. ! needed?

          ulgas = u0gas

          ! Multiply gas amounts by rundeck scaling factors.
          ! Looping not an option since fulgas array does not yet
          ! contain the factors.

          !ulgas(:, 1) = ulgas(:, 1)*H2OstratX
          ulgas(:, 2) = ulgas(:, 2)*CO2X
          ulgas(:, 3) = ulgas(:, 3)*O3X
          ulgas(:, 4) = ulgas(:, 4)*O2X
          ulgas(:, 5) = ulgas(:, 5)*NO2X
          ulgas(:, 6) = ulgas(:, 6)*N2OX
          ulgas(:, 7) = ulgas(:, 7)*CH4X
          ulgas(:, 8) = ulgas(:, 8)*CFC11X
          ulgas(:, 9) = ulgas(:, 9)*CFC12X
          ulgas(:,10) = ulgas(:,10)*N2CX
          ulgas(:,11) = ulgas(:,11)*XGHGX
          ulgas(:,12) = ulgas(:,12)*YGHGX
          ulgas(:,13) = ulgas(:,13)*SO2X

        endif
        if(file_exists('VISAODangstr')) then
          set_aerosols_internally = .false.
c          fid = par_open(grid,'VISAODangstr','read')
c          not needed for initial CIRC cases which have zero aerosol
c          todo:  read optical depths and scale with Angstrom exponent
c          weighted by solar flux
c          ....
c          call par_close(grid,fid)
          sraext = 0.; srasct = 0.; sragcb = 0.
          srdext = 0.; srdsct = 0.; srdgcb = 0.
          srvext = 0.; srvsct = 0.; srvgcb = 0.
          srbext = 0.; srbsct = 0.; srbgcb = 0.
          traalk = 0.
          trdalk = 0.
          trvalk = 0.
          trbalk = 0.
        endif
        i = 1
        j = 1
        do l=1,lm
          tloc = t(i,j,l)*pk(l,i,j)
          if(tloc.ge.tf) then
            svlhx(l,i,j) = lhe
          else
            svlhx(l,i,j) = lhs
          endif
          svlat(l,i,j) = svlhx(l,i,j)
          rhsav(l,i,j) = q(i,j,l)/qsat(tloc,svlhx(l,i,j),pmid(l,i,j))
        enddo
        !llow=1; lmid=2; lhi=3
      endif

      RETURN
      END SUBROUTINE init_RAD

      subroutine SETATM             ! dummy routine in gcm
      end subroutine SETATM

      subroutine GETVEG(LONR,LATR)  ! dummy routine in gcm
      integer LONR,LATR
      end subroutine GETVEG

      SUBROUTINE daily_RAD(end_of_day)
!@sum  daily_RAD sets radiation parameters that change every day
!@auth G. Schmidt
!@calls RADPAR:RCOMPT
      USE DOMAIN_DECOMP_ATM, only : am_I_root
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      use model_com, only: modelEclock
      USE RADPAR, only : FULGAS,JYEARR=>JYEAR,JDAYR=>JDAY
     *     ,xref,KYEARV
#ifdef ALTER_RADF_BY_LAT
     *     ,FULGAS_orig
#endif
      USE RADPAR, only : rcompt,writet
      USE RAD_COM, only : co2x,n2ox,ch4x,cfc11x,cfc12x,xGHGx,h2ostratx
     *     ,o2x,no2x,n2cx,yghgx,so2x
     *     ,o3x,o3_yr,ghg_yr,co2ppm,Volc_yr,albsn_yr,dalbsnX
     *     ,snoage,snoage_def,chl_from_seawifs
      use DIAG_COM, only : iwrite,jwrite,itwrite,tdiurn
      use geom, only : imaxj
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: end_of_day
      integer :: year, dayOfYear
      integer :: i,j, i_0,i_1,j_0,j_1, itype

      call modelEclock%get(year=year, dayOfYear=dayOfYear)

C**** Update time dependent radiative parameters each day
!     Get black carbon deposition data for the appropriate year
!     (does nothing except at a restart or the beginning of a new year)
      if(dalbsnX.ne.0.) then
        if (albsn_yr.eq.0) then
#ifdef OLD_BCdalbsn
          call updBCd (year)
#else
          call updBCdalbsn (year    ,dayofyear)
#endif
        else
#ifdef OLD_BCdalbsn
          call updBCd (albsn_yr)
#else
          ! as per radiation-code convention, pass -albsn_yr to indicate
          ! perpetual-year mode
          call updBCdalbsn (-albsn_yr,dayofyear)
#endif
        end if
      endif
!     Hack: 2 specific volc. eruption scenarios for 2000-2100 period
      if(volc_yr.eq.-2010) then              ! repeat some old volcanos
         KYEARV=YEAR
         if(YEAR.GT.2010) KYEARV=YEAR-100  ! go back 100 years
      end if
      if(volc_yr.eq.-2000) then
         KYEARV=YEAR
         if(YEAR.GT.2000) KYEARV=YEAR-50   ! go back 50 years til 2050
         if(YEAR.GT.2050) KYEARV=YEAR-150  ! then go back 150 years
      end if

      JDAYR=dayOfYear
      JYEARR=YEAR
      CALL RCOMPT
!     FULGAS(2:) is set only in the first call to RCOMPT unless ghg_yr=0
!     Optional scaling of the observed value only in case it was (re)set
      if(.not. end_of_day .and. H2OstratX.GE.0.)
     *   FULGAS(1)=FULGAS(1)*H2OstratX
      if(.not. end_of_day .or. O3_yr==0.)
     *   FULGAS(3)=FULGAS(3)*O3X
      if(ghg_yr.eq.0 .or. .not. end_of_day) then
         FULGAS(2)=FULGAS(2)*CO2X
         FULGAS(6)=FULGAS(6)*N2OX
         FULGAS(7)=FULGAS(7)*CH4X
         FULGAS(8)=FULGAS(8)*CFC11X
         FULGAS(9)=FULGAS(9)*CFC12X
         FULGAS(11)=FULGAS(11)*XGHGX
         FULGAS(12)=FULGAS(12)*YGHGX
      end if
      if(.not. end_of_day) then
         FULGAS(4)=FULGAS(4)*O2X
         FULGAS(5)=FULGAS(5)*NO2X
         FULGAS(10)=FULGAS(10)*N2CX
         FULGAS(13)=FULGAS(13)*SO2X ! no effect since FULGAS(13)=0.
      end if

C**** write trend table for forcing 'itwrite' for years iwrite->jwrite
C**** itwrite: 1-2=GHG 3=So 4-5=O3 6-9=aerosols: Trop,DesDust,Volc,Total
      if (am_i_root() .and.
     *  jwrite.gt.1500) call writet (6,itwrite,iwrite,jwrite,1,0)

#ifdef ALTER_RADF_BY_LAT
C**** Save initial rad forcing alterations:
      FULGAS_orig(:)=FULGAS(:) ! GHGs
#endif

C**** Define CO2 (ppm) for rest of model
      co2ppm = FULGAS(2)*XREF(1)

      if (chl_from_seawifs>0) call get_chl_from_seawifs

      if (end_of_day) then

        call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
        call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)

        do j=J_0,J_1
        do i=I_0,imaxj(j)
c****
c**** increase snow age depending on snoage_def
c****
          if (snoage_def.eq.0) then ! update indep. of ts
            do itype=1,3
              snoage(itype,i,j)=1.+.98d0*snoage(itype,i,j)
            end do
          elseif (snoage_def.eq.1) then ! update if max T>0
            if (tdiurn(i,j,7).gt.0) snoage(1,i,j)=1.+.98d0
     *           *snoage(1,i,j) ! ocean ice (not currently used)
            if (tdiurn(i,j,8).gt.0) snoage(2,i,j)=1.+.98d0
     *           *snoage(2,i,j) ! land ice
            if (tdiurn(i,j,2).gt.0) snoage(3,i,j)=1.+.98d0
     *           *snoage(3,i,j) ! land
          else
            write(6,*) "This snoage_def is not defined: ",snoage_def
            write(6,*) "Please use: 0 (update indep of T)"
            write(6,*) "            1 (update if T>0)"
            call stop_model('stopped in RAD_DRV.f',255)
          end if
        enddo
        enddo
      endif

      RETURN
      END SUBROUTINE daily_RAD


      subroutine get_chl_from_seawifs

      USE DOMAIN_DECOMP_ATM, only : GRID,REWIND_PARALLEL
     .     ,READT_PARALLEL, getDomainBounds
      USE FLUXES, only : focean,atmocn
      USE CONSTANT, only : by12
      use model_com, only: modelEclock, calendar
      USE RESOLUTION, only : im,jm
      USE FILEMANAGER, only : NAMEUNIT
      USE GEOM, only : imaxj
      USE filemanager, only: openunit
      USE CalendarMonth_mod
      implicit none

      REAL*8 :: TEMP_LOCAL(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     .                     GRID%J_STRT_HALO:GRID%J_STOP_HALO,2)
      integer :: month, date, year
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
      INTEGER :: LSTMON, I, J, J_0, J_1, I_0, I_1
      INTEGER, SAVE :: IMON0 = 0
      integer, save :: iu_chl=-1
!@var ACHL,ECHL1,ECHL0,BCHL,CCHL arrays for the reading in chlorophyll
      REAL*8, ALLOCATABLE, DIMENSION(:,:), save :: ACHL,ECHL1,ECHL0,
     .                   BCHL, CCHL
      REAL*8 :: TIME
      INTEGER :: I_0H, I_1H, J_0H, J_1H
      type (CalendarMonth) :: cMonth

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO
      if (iu_chl.lt.0) then
        call openunit("CHL_DATA",iu_CHL,.true.,.true.)
        allocate(ACHL(I_0H:I_1H,J_0H:J_1H),
     .           ECHL1(I_0H:I_1H,J_0H:J_1H),
     .           ECHL0(I_0H:I_1H,J_0H:J_1H),
     .           BCHL(I_0H:I_1H,J_0H:J_1H),
     .           CCHL(I_0H:I_1H,J_0H:J_1H))
      endif
      call modelEclock%get(month=month, date=date)
      call getDomainBounds(GRID,J_STRT=J_0,J_STOP=J_1,
     .         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     .         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Read in Seawifs files here
      IF (month.NE.IMON0) THEN
      IF (IMON0==0) THEN
C**** READ IN LAST MONTH''S END-OF-MONTH DATA
        LSTMON=month-1
        if (lstmon.eq.0) lstmon = 12
        CALL READT_PARALLEL
     *       (grid,iu_CHL,NAMEUNIT(iu_CHL),TEMP_LOCAL,LSTMON)
        ECHL0 = TEMP_LOCAL(:,:,2)
      ELSE
C**** COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
        ECHL0=ECHL1
      END IF
C**** READ IN CURRENT MONTHS DATA: MEAN AND END-OF-MONTH
      IMON0=month
      if (month.eq.1) CALL REWIND_PARALLEL( iu_CHL )
      CALL READT_PARALLEL
     *     (grid,iu_CHL,NAMEUNIT(iu_CHL),TEMP_LOCAL,1)
      ACHL  = TEMP_LOCAL(:,:,1)
      ECHL1 = TEMP_LOCAL(:,:,2)

C**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          BCHL(I,J)=ECHL1(I,J)-ECHL0(I,J)
          CCHL(I,J)=3.*(ECHL1(I,J)+ECHL0(I,J))-6.*ACHL(I,J)
        END DO
      END DO
      END IF
C**** Calculate CHL for current day
      cMonth = calendar%getCalendarMonth(month, year)
      TIME=(DATE-.5)/cMonth%daysinMonth -.5 ! -.5<TIME<.5
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0) THEN
C**** CHL always uses quadratic fit
            atmocn%CHL(I,J)=ACHL(I,J)+BCHL(I,J)*TIME
     *           +CCHL(I,J)*(TIME**2-BY12)
            if (atmocn%CHL(I,J).lt. 0) atmocn%CHL(I,J)=0. ! just in case
          END IF
        END DO
      END DO
C**** REPLICATE VALUES AT POLE
      IF(HAVE_NORTH_POLE) then
       if (FOCEAN(1,JM).gt.0) atmocn%CHL(2:IM,JM)=atmocn%CHL(1,JM)
      ENDIF
      IF(HAVE_SOUTH_POLE) then
       if (FOCEAN(1, 1).gt.0) atmocn%CHL(2:IM, 1)=atmocn%CHL(1, 1)
      ENDIF
      atmocn%chl_defined=.true.
      return

      end subroutine get_chl_from_seawifs

      SUBROUTINE DAILY_orbit(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@auth Original Development Team
!@calls constant:orbit
      use model_com, only: modelEclock
      USE RAD_COM, only : RSDIST,COSD,SIND,COSZ_day,SUNSET,
     *     variable_orb_par,orb_par_year_bp, useOrbit => orbit
      USE DOMAIN_DECOMP_ATM, only : am_I_root
      use RAD_COSZ0, only : daily_cosz
      use BaseTime_mod
      use TimeInterval_mod
      use Rational_mod
      IMPLICIT NONE
      REAL*8 :: SUNLON,SUNLAT,LAM,EDPY,VEDAY,PYEAR
      LOGICAL, INTENT(IN) :: end_of_day
      integer :: year, dayOfYear
      type (BaseTime) :: t
      real*8 :: declinationAngle
      type (TimeInterval) :: halfDay

      call modelEclock%get(year=year, dayOfYear=dayOfYear)

C**** CALCULATE SOLAR ANGLES AND ORBIT POSITION
C**** This is for noon (GMT) for new day.

C**** The orbital calculation will need to vary depending on the kind
C**** of calendar adopted (i.e. a generic 365 day year, or a transient
C**** calendar including leap years etc.).  For transient calendars the
C**** dayOfYear passed to orbit needs to be adjusted to represent the number
C**** of days from Jan 1 2000AD.
c      EDPY=365.2425d0, VEDAY=79.3125d0  ! YR 2000AD
c      dayOfYear => dayOfYear + 365 * (YEAR-2000) + appropriate number of leaps
C**** Default calculation (no leap, VE=Mar 21 hr 0)
c      EDPY=365d0 ; VEDAY=79d0           ! Generic year
C**** PMIP calculation (no leap, VE=Mar 21 hr 12)
      EDPY=365d0 ; VEDAY=79.5d0           ! Generic year
C**** Update orbital parameters at start of year
      if (dayOfYear == 1) then
         call useOrbit%setYear(real(year,kind=8))
      end if

      ! Use time for the _middle_ of the day to compute
      ! zenith angle:

      halfDay = TimeInterval(useOrbit%getMeanDay() / 2)
      t = newBaseTime(modelEClock%getTimeAtBeginningOfCurrentDay() +
     *                halfDay)

      sinD = useOrbit%getSinDeclinationAngle(t)
      cosD = sqrt(1-sinD**2)
      rsdist = useOrbit%getDistance(t)**2

      call daily_cosz(sind,cosd,cosz_day,sunset)

      RETURN
      END SUBROUTINE DAILY_orbit


      SUBROUTINE DAILY_ch4ox(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@vers 2013/03/27
!@auth Original Development Team
!@calls constant:orbit
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : q
      use model_com, only: modelEclock
      USE MODEL_COM, only : itime
      USE GEOM, only : axyp,imaxj,lat2d
      USE ATM_COM, only : byMA
      USE RADPAR, only : ghgam,ghgyr2,ghgyr1
      USE RAD_COM, only : dh2o,H2ObyCH4,ghg_yr
#ifdef TRACERS_WATER
      use OldTracer_mod, only: tr_wd_type, nWATER,tr_H2ObyCH4, itime_tr0
      USE TRACER_COM, only: trm,NTM
#endif
      USE DIAG_COM, only : ftype,ntype, aij=>aij_loc
      USE DIAG_COM_RAD, only : j_h2och4, ij_h2och4
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, am_I_root
      IMPLICIT NONE
      REAL*8 :: xCH4,xdH2O
      INTEGER i,j,l,iy,it
      LOGICAL, INTENT(IN) :: end_of_day
#ifdef TRACERS_WATER
      INTEGER n
#endif
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      integer :: year, month

      call modelEclock%get(year=year, month=month)

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      IF (.not.end_of_day) RETURN

C**** Tasks to be done at end of day only
      if (H2ObyCH4.gt.0) then
C****   Add obs. H2O generated by CH4(*H2ObyCH4) using a 2 year lag
        iy = year - 2 - ghgyr1 + 1
        if (ghg_yr.gt.0) iy = ghg_yr - 2 - ghgyr1 + 1
        if (iy.lt.1) iy=1
        if (iy.gt.ghgyr2-ghgyr1+1) iy=ghgyr2-ghgyr1+1
        xCH4=ghgam(3,iy)*H2ObyCH4
c        If (AM_I_ROOT())
c     &    write(6,*) 'add in stratosphere: H2O gen. by CH4(ppm)=',xCH4

        do l=1,lm
        do j=J_0,J_1
        do i=I_0,imaxj(j)
#ifdef CUBED_SPHERE
          call lat_interp_qma(lat2d(i,j),l,month,xdH2O)
#else
          xdH2O = dH2O(j,l,month)
#endif
          q(i,j,l)=q(i,j,l)+xCH4*xdH2O*byMA(l,i,j)
#ifdef TRACERS_WATER
C**** Add water to relevant tracers as well
          do n=1,ntm
            if (itime_tr0(n).le.itime) then
              select case (tr_wd_type(n))
              case (nWater)    ! water: add CH4-sourced water to tracers
                trm(i,j,l,n) = trm(i,j,l,n) +
     +                tr_H2ObyCH4(n)*xCH4*xdH2O*axyp(i,j)
              end select
            end if
          end do
#endif
          do it=1,ntype
            call inc_aj(i,j,it,j_h2och4,xCH4*xdH2O*ftype(it,i,j))
          end do
          aij(i,j,ij_h2och4) = aij(i,j,ij_h2och4) + xCH4 * xdH2O
        end do
        end do
        If (HAVE_NORTH_POLE) q(2:im,jm,l)=q(1,jm,l)
        If (HAVE_SOUTH_POLE) q(2:im, 1,l)=q(1, 1,l)
#ifdef TRACERS_WATER
        do n=1,ntm
          If (HAVE_SOUTH_POLE) trm(2:im, 1,l,n)=trm(1, 1,l,n)
          If (HAVE_NORTH_POLE) trm(2:im,jm,l,n)=trm(1,jm,l,n)
        end do
#endif
        end do
      end if

      RETURN
      END SUBROUTINE DAILY_ch4ox

      SUBROUTINE RADIA
!@sum  RADIA adds the radiation heating to the temperatures
!@vers 2013/03/27
!@auth Original Development Team
!@calls tropwmo,coszs,coszt, RADPAR:rcompx ! writer,writet
      USE CONSTANT, only : lhe,lhs,twopi,tf,stbo,rhow,mair,grav
     *     ,bysha,pi,radian,areag
      USE RESOLUTION, only : pmtop
      USE RESOLUTION, only : im,jm,lm
#ifdef TRACERS_SPECIAL_Shindell
      USE RESOLUTION, only : LS1=>ls1_nominal
#endif
      USE ATM_COM, only : kradia,lm_req,p,t,q,iu_rad,req_fac_d
      USE MODEL_COM
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR
      USE ATM_COM, only : byaml00
      USE GEOM, only : imaxj, axyp, byaxyp
     &     ,lat2d,lon2d
      USE RADPAR
     &  , only :  ! routines
     &           lx  ! for threadprivate copyin common block
     &          ,tauwc0,tauic0 ! set in radpar block data
     &          ,writer,rcompx,updghg
C     INPUT DATA         ! not (i,j) dependent
     X          ,S00WM2,RATLS0,S0,JYEARR=>JYEAR,JDAYR=>JDAY,FULGAS
     &          ,use_tracer_chem,FS8OPX,FT8OPX,use_o3_ref,KYEARG,KJDAYG
     &          ,planck_tmin,planck_tmax
#ifdef ALTER_RADF_BY_LAT
     &          ,FS8OPX_orig,FT8OPX_orig,FULGAS_orig
#endif
C     INPUT DATA  (i,j) dependent
     &             ,JLAT46=>JLAT,ILON72=>ILON,JGCM,IGCM
     &             ,L1,LMR=>NL, PLB ,TLB,TLM ,SHL,RHL
     &             ,ltopcl,TAUWC ,TAUIC ,SIZEWC ,SIZEIC, kdeliq
     &             ,POCEAN,PEARTH,POICE,PLICE,PLAKE,COSZ,PVT
     &             ,TGO,TGE,TGOI,TGLI,TSL,WMAG,WEARTH
     &             ,AGESN,SNOWD,SNOWOI,SNOWLI,dALBsn, ZSNWOI,ZOICE
     &             ,zmp,fmp,flags,LS1_loc,snow_frac,zlake
     *             ,TRACER,FSTOPX,FTTOPX,chem_IN
     &             ,nraero_aod=>NTRACE
     *             ,FTAUC,LOC_CHL,FSTASC,FTTASC
#ifdef HEALY_LM_DIAGS
     *             ,VTAULAT
#endif
#ifdef GCC_COUPLE_RAD
     *             ,GCCco2_IN, use_tracer_GCCco2, GCCco2_out
#endif

C     OUTPUT DATA
     &          ,TRDFLB ,TRNFLB ,TRUFLB, TRFCRL ,chem_out
     &          ,SRDFLB ,SRNFLB ,SRUFLB, SRFHRL
     &          ,PLAVIS ,PLANIR ,ALBVIS ,ALBNIR ,FSRNFG
     &          ,SRRVIS ,SRRNIR ,SRAVIS ,SRANIR ,SRXVIS ,SRDVIS
     &          ,BTEMPW ,SRAEXT ,SRASCT ,SRAGCB
     &          ,SRDEXT ,SRDSCT ,SRDGCB ,SRVEXT ,SRVSCT ,SRVGCB
     &          ,aesqex,aesqsc,aesqcb,CO2outCol
     &          ,aesqex_dry,aesqsc_dry,aesqcb_dry
     &          ,SRXNIR,SRDNIR
      USE RAD_COM, only : modrd,nrad
      USE RAD_COM, only : rqt,srhr,trhr,fsf,cosz1,s0x,rsdist,nradfrc
     *     ,CH4X_RADoverCHEM,snoage
     *     ,plb0,shl0,tchg,alb,fsrdir,srvissurf,srdn,cfrac,rcld
     *     ,chem_tracer_save,rad_interact_aer,kliq,RHfix,CLDx
     *     ,ghg_yr,CO2X,N2OX,CH4X,CFC11X,CFC12X,XGHGX,rad_forc_lev
     *     ,ntrix_aod,ntrix_rf,wttr,cloud_rad_forc,CC_cdncx,OD_cdncx
     *     ,cdncl,dALBsnX,rad_to_chem,trsurf,dirvis
     *     ,FSRDIF,DIRNIR,DIFNIR,aer_rad_forc,clim_interact_chem
     *     ,TAUSUMW,TAUSUMI,TAero_aod_diag
     *     ,chl_from_obio,chl_from_seawifs
#ifdef GCC_COUPLE_RAD
     *     ,GCCco2_tracer_save,GCCco2rad_to_chem
#endif
#ifdef mjo_subdd
     *     ,SWHR,LWHR,SWHR_cnt,LWHR_cnt,OLR_acc,OLR_cnt
     *     ,swu_avg,swu_cnt
#endif
#ifdef ALTER_RADF_BY_LAT
     *     ,FULGAS_lat,FS8OPX_lat,FT8OPX_lat
#endif
#ifdef TRACERS_DUST
     &     ,srnflb_save,trnflb_save
#endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &     ,stratO3_tracer_save
#endif
#ifdef TRACERS_ON
      use rad_com, only: tau_as,tau_cs,tau_dry,nraero_rf
#ifdef CACHED_SUBDD
      USE CONSTANT, only : grav,Rgas
      use rad_com, only: abstau_as,abstau_cs,abstau_dry,swfrc,lwfrc
      use RunTimeControls_mod, only: tracers_amp, tracers_tomas
#endif  /* CACHED_SUBDD */
#endif
      USE RANDOM
      USE CLOUDS_COM, only : tauss,taumc,svlhx,rhsav,svlat,cldsav,
     *     cldmc,cldss,csizmc,csizss,llow,lmid,lhi,fss,taussip,csizssip
     *    ,QLss,QIss,QLmc,QImc
     *    ,get_cld_overlap  !  subroutine
      USE DIAG_COM, only : ia_rad,jreg,aij=>aij_loc,aijl=>aijl_loc
     &     ,ntype,ftype,itocean,itlake,itearth,itlandi,itoice,itlkice
     *     ,adiurn=>adiurn_loc,ndiuvar,ia_rad_frc,
#ifdef USE_HDIURN
     *     hdiurn=>hdiurn_loc,
#endif
     *     iwrite,jwrite,itwrite,ndiupt
     *     ,ijdd
     *     ,AFLX_ST, hr_in_day,hr_in_month
      USE DIAG_COM_RAD
#ifdef TRACERS_ON
      USE DIAG_COM, only : adiurn_dust,save3dAOD
      use RAD_COM, only: diag_fc
#endif
      USE ATM_COM, only : pk,pedn,pmid,pdsig,ltropo,MA,byMA
      USE SEAICE_COM, only : si_atm
      USE GHY_COM, only : fearth,snowd_ij=>snowd
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
     &                    ,N_COVERTYPES !YKIM-temp hack
      use ent_drv, only : map_ent2giss  !YKIM-temp hack
      USE LAKES_COM, only : flake,dlake!,mwl
      USE FLUXES, only : asflx4,atmocn,atmice,atmgla,atmlnd,atmsrf
     &     ,flice,fland,focean
      USE DOMAIN_DECOMP_ATM, ONLY: grid, write_parallel
      USE DOMAIN_DECOMP_ATM, ONLY: GLOBALSUM, getDomainBounds
      USE RAD_COSZ0, only : COSZT,COSZS

#ifdef TRACERS_ON
      use OldTracer_mod, only: trname, trpdens
      USE TRACER_COM, only: NTM
     *     ,n_Ox,trm,n_OCB,n_BCII,n_BCIA
     *     ,n_OCIA,N_OCII,n_so4_d2,n_so4_d3,n_SO4,n_stratOx
     *     ,n_N_AKK_1
#ifdef TRACERS_NITRATE
      use OldTracer_mod, only: tr_mm
      use TRACER_COM, only: n_NH4,n_NO3p
#endif
#ifdef TRACERS_AEROSOLS_SOA
      use TRACER_COM, only: n_isopp1a,n_isopp2a
#ifdef TRACERS_TERP
      use TRACER_COM, only: n_apinp1a,n_apinp2a
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_OCEAN
      use TRACER_COM, only: n_ococean
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef GCC_COUPLE_RAD
      USE TRACER_COM, only: n_CO2n
      USE CONSTANT, only: avog
      USE OldTracer_mod, only: tr_mm
#endif
#ifdef TRACERS_AEROSOLS_VBS
      use TRACERS_VBS, only: vbs_tr
#endif
      USE TRDIAG_COM, only: taijs=>taijs_loc,taijls=>taijls_loc,ijts_fc
     *     ,ijts_tau,ijts_tausub,ijts_fcsub
     *     ,ijlt_3dtau,ijlt_3daaod,ijlt_3dtauCS,ijlt_3daaodCS
     *     ,ijlt_3dtauDRY,ijlt_3daaodDRY
     *     ,ijts_sqex
     *     ,ijts_sqexsub,ijts_sqsc,ijts_sqscsub,ijts_sqcb,ijts_sqcbsub
     *     ,diag_rad,diag_aod_3d,save_dry_aod
#ifdef AUXILIARY_OX_RADF
     *     ,ijts_auxfc
#endif /* AUXILIARY_OX_RADF */
#ifdef BC_ALB
     *     ,ijts_alb,ijts_sunlit_snow
#endif  /* BC_ALB */
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: Lmax_rad_O3,Lmax_rad_CH4
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_TOMAS
      USE TOMAS_AEROSOL, only: icomp
#endif
#endif /* TRACERS_ON */
      use AerParam_mod, only: dCDNC_est
#ifdef OLD_BCdalbsn
      use AerParam_mod, only: depoBC,depoBC_1990
#else
      use AerParam_mod, only: BCdalbsn
#endif
      USE TimerPackage_mod, only: startTimer => start, stopTimer => stop
      USE Dictionary_mod, only : get_param, is_set_param
#ifdef CACHED_SUBDD
      use subdd_mod, only : sched_rad, subdd_groups, subdd_type
     &     ,subdd_ngroups,inc_subdd,find_groups, lmaxsubdd
#endif
#ifdef SCM
      use SCM_COM, only : SCMopt,SCMin
      USE CONSTANT, only : SHA
      USE ATM_COM, only : QCL
#endif
      use DIAG_COM, only: ij_nintaerext,ij_nintaersca,ij_nintaerasy
      use RADPAR, only: nintaerext,nintaersca,nintaerasy
      IMPLICIT NONE
      real*8 dz,rho
C
      !@var wtrtau,icetau per-layer opacity for cloud water,ice
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &        wtrtau,icetau
#ifdef SCM
      real*8 q_above(LM+1),q_below(LM+1),Frad(LM+1)
#endif
C     INPUT DATA   partly (i,j) dependent, partly global
      REAL*8 U0GAS,taulim
#ifdef OLD_BCdalbsn
      REAL*8 xdalbs,sumda,tauda,fsnow
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     sumda_psum,tauda_psum
#endif
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     COSZ2,COSZA,TRINCG,BTMPW,WSOIL,fmp_com
      REAL*8, DIMENSION(4,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                    grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     SNFS,TNFS
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     SNFSCRF,TNFSCRF,SNFSCRF2,TNFSCRF2,LWDNCS,
     &     SNFS_AS_noA, TNFS_AS_noA, SNFS_CS_noA, TNFS_CS_noA,
     &     SWUS,CTT,CTP,WTRCLD,ICECLD
      REAL*8, DIMENSION(18,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     SNFSAERRF,TNFSAERRF
#ifdef CFMIP3_SUBDD
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &        swut,swutcs,cfmip_twp,swdcls,swucls,swdt
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &        cfmip_cf,cfmip_qci,cfmip_qcl
#endif
#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups)
      type(subdd_type), pointer :: subdd
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     SDDARR
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     SDDARR3D
#ifdef TRACERS_ON
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,nraero_rf) ::
     &     sddarr3drf
      integer :: f
#endif  /* TRACERS_ON */
#ifdef SCM
C     radiative flux profiles for sub-daily output, generalized
C     for GCM grid but currently limited to SCM use
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     TRDFLB_prof, TRUFLB_prof, SRDFLB_prof, SRUFLB_prof
#endif
#ifdef TRACERS_ON
! types of aods to be saved
! The name will be any combination of {,TRNAME}{as,cs}{,a}aod
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,
     &                  lm,nraero_aod) ::
     &     sddarr4d
      character(len=10), dimension(2) :: sgroups=(/'taijh ','taijlh'/)
      character(len=10), dimension(3) :: ssky=(/'as ','cs ','dry'/)
      character(len=10), dimension(2) :: sabs=(/' ','a'/)
      character(len=10), dimension(2) :: sfrc=(/'swf','lwf'/)
      character(len=10) :: spcname
      character(len=50) :: sname
      integer :: g,s,a
#endif  /* TRACERS_ON */
      !@var CO2out for holding 3D CO2 from rad code for SUBDD
      REAL*8, dimension(LM,grid%i_strt_halo:grid%i_stop_halo,
     & grid%j_strt_halo:grid%j_stop_halo) :: CO2out
#endif /* CACHED_SUBDD */
#if (defined ACCMIP_LIKE_DIAGS) 
#ifndef SKIP_ACCMIP_GHG_RADF_DIAGS
!@var snfs_ghg,tnfs_ghg like SNFS/TNFS but with reference GHG for
!@+   radiative forcing calculations. TOA only.
!@+   index 1=CH4, 2=N2O, 3=CFC11, 4=CFC12
      real*8,dimension(4,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &                   snfs_ghg,tnfs_ghg
      real*8,dimension(4) :: sv_fulgas_ref,sv_fulgas_now
      integer :: nf,GFrefY,GFrefD,GFnowY,GFnowD
!@var nfghg fulgas( ) index of radf diag ghgs:
      integer, dimension(4) :: nfghg=(/7,6,8,9/)
#endif
#endif

! variables for running uncoupled concentration-driven GCC
#ifdef GCC_UNCOUPLE_RAD_CONCEN
      real*8  :: GCCco2_fulgas_ref,GCCco2_fulgas_now
      integer :: GCCco2nowY,GCCco2nowD
#endif

#ifdef HEALY_LM_DIAGS
C  GHG Effective forcing relative to 1850
       real*8 :: ghg_totforc, CO2I=285.2,N2OI=.2754,CH4I=.791  ! 1850  GHG's
       real*8 ::              CO2R=337.9,N2OR=.3012,CH4R=1.547 ! RAD's 1979 Reference values
       real*8 ::              FCO2,FN2O,FCH4                   ! Current Model GHG
       real*8 :: Fe !! Function
#endif
#ifdef TRACERS_ON
!@var SNFST,TNFST like SNFS/TNFS but with/without specific tracers for
!@+   radiative forcing calculations
      REAL*8,DIMENSION(2,nraero_rf,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                             grid%J_STRT_HALO:grid%J_STOP_HALO)::
     *     SNFST,TNFST
!@var SNFST_o3ref,TNFST_o3ref like snfst,tnfst for special case ozone for
!@+   which nraero_rf fields are not defined. Indicies are :
!@+   1=LTROPO,reference, 2=TOA,reference; not saving surface forcing.
!@+   3=LTROPO or LS1-1,auxiliary, 4=TOA,auxiliary; 5=LS1-1,reference
      REAL*8,DIMENSION(5,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     SNFST_o3ref,TNFST_o3ref
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &    ,snfst_stratOx,tnfst_stratOx
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
#ifdef BC_ALB
      REAL*8,DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                 grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     ALBNBC,NFSNBC,
     &     dALBsnBC ! not to be confused with BCdalbsn from an input file
      LOGICAL,DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     bc_snow_present
#endif /* BC_ALB */
#endif /* TRACERS_ON */
      REAL*8, DIMENSION(LM_REQ,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                         grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     TRHRS,SRHRS
      REAL*8, DIMENSION(0:LM+LM_REQ,
     *     grid%I_STRT_HALO:grid%I_STOP_HALO,
     *     grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     TRHRA,SRHRA ! for adj.frc
      REAL*8, DIMENSION(LM) :: TOTCLD,SS_CLD,dcc_cdncl,dod_cdncl
      INTEGER I,J,L,K,KR,LR,JR,IH,IHM,INCH,JK,IT,iy,iend,N,onoff_aer
     *     ,onoff_chem,LFRC,JTIME,n1,moddrf
      REAL*8 ROT1,ROT2,PLAND,CSS,CMC,DEPTH,QSS,TAUSSL,TAUSSLIP
     *     ,TAUMCL,ELHX,CLDCV,X,OPNSKY,CSZ2,tauup,taudn,ptype4(4)
     *     ,taucl,wtlin,MSTRAT,STRATQ,STRJ,MSTJ,optdw,optdi,rsign_aer
     *     ,rsign_chem,tauex5,tauex6,tausct,taugcb,dcdnc
     *     ,QR(LM,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &            grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,CLDinfo(LM,3,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8 tmpS(8),tmpT(8)
      REAL*8 QSAT
#ifdef BC_ALB
      REAL*8 dALBsn1
#endif
      LOGICAL set_clayilli,set_claykaol,set_claysmec,
     &     set_claycalc,set_clayquar
C
      REAL*8  RDSS(LM,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,RDMC(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &           grid%J_STRT_HALO:grid%J_STOP_HALO)

      REAL*8 :: TMP(NDIUVAR)
      INTEGER, PARAMETER :: NLOC_DIU_VAR = 8
      INTEGER :: idx(NLOC_DIU_VAR)
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      INTEGER, PARAMETER :: NLOC_DIU_VARb = 5
#else
      INTEGER, PARAMETER :: NLOC_DIU_VARb = 3
#endif
      INTEGER :: idxb(NLOC_DIU_VARb)

      integer :: aj_alb_inds(8)
      real*8, dimension(lm_req) :: bydpreq

c     INTEGER ICKERR,JCKERR,KCKERR
      INTEGER :: J_0, J_1, I_0, I_1
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      character(len=300) :: out_line

      integer :: nij_before_j0,nij_after_j1,nij_after_i1
      integer :: initial_GHG_setup

      real*8 :: PVT0(N_COVERTYPES), HVT0(N_COVERTYPES)
#ifdef TRACERS_NITRATE
      real*8 :: nh4_on_no3
#endif
#ifdef TRACERS_TOMAS
      real*8 :: qcb_col(6,ICOMP-2),qcb_col_dry(6,ICOMP-2)
#endif

      REAL*8, DIMENSION(:,:), POINTER :: RSI,ZSI,SNOWI,POND_MELT
      LOGICAL, DIMENSION(:,:), POINTER :: FLAG_DSWS
      real*8 :: rhodz ! air density times layer thickness (kg/m2
      integer :: year, dayOfYear, hour, date

#ifdef TRACERS_ON
!@var nsub_ntrix  array of index counters for sub classes of tracers
      integer, dimension( ntm ) :: nsub_ntrix
#endif

#ifdef GCC_COUPLE_RAD
      INTEGER :: Lmax_rad_CO2 = LM
#endif

      call modelEclock%get(year=year, dayOfYear=dayOfYear,
     *     hour=hour, date=date)

      RSI => SI_ATM%RSI
      ZSI => SI_ATM%ZSI
      SNOWI => SI_ATM%SNOWI
      POND_MELT => SI_ATM%POND_MELT
      FLAG_DSWS => SI_ATM%FLAG_DSWS

C
C****
      call startTimer('RADIA()')

      idx = (/ (IDD_CL7+i-1,i=1,7), IDD_CCV /)
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      idxb = (/ IDD_PALB, IDD_GALB, IDD_ABSA, idd_aot, idd_aot2 /)
#else
      idxb = (/ IDD_PALB, IDD_GALB, IDD_ABSA /)
#endif
      call getDomainBounds(grid, HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP


C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** GTEMPR RADIATIVE TEMPERATURE ARRAY OVER ALL SURFACE TYPES (K)
C****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****
C**** VDATA  1-11 RATIOS FOR THE 11 VEGETATION TYPES (1)
C****

C**** limit optical cloud depth from below: taulim
      taulim=min(tauwc0,tauic0) ! currently both .001
      tauwc0 = taulim ; tauic0 = taulim
C**** Calculate mean cosine of zenith angle for the current physics step
      JTIME=MOD(ITIME,NDAY)
      ROT1=(TWOPI*JTIME)/NDAY
c      ROT2=ROT1+TWOPI*DTsrc/SECONDS_PER_DAY
c      CALL COSZT (ROT1,ROT2,COSZ1)
      CALL CALC_ZENITH_ANGLE   ! moved to main loop

      if (kradia.gt.0) then    ! read in all rad. input data (frc.runs)
        iend = 1
        it = itime-1           ! make sure, at least 1 record is read
        do while (mod(itime-it,NDAY*INT_DAYS_PER_YEAR).ne.0)
          read(iu_rad,end=10,err=10) it
C****   input data:          WARNINGS
C****        1 - any changes here also go in later (look for 'iu_rad')
C****        2 - keep "dimrad_sv" up-to-date:         dimrad_sv=IM*JM*{
     *     ,T,RQT,atmsrf%TsAvg                         ! LM+LM_REQ+1+
     *     ,QR,P,CLDinfo,rsi,zsi                       ! LM+1+3*LM+1+1+
!     *     ,(((GTEMPR(k,i,j),k=1,4),i=1,im),j=1,jm)    ! (4+)
     *     ,wsoil,atmsrf%wsavg,snowi,atmgla%snow,atmlnd%snowe ! 1+1+1+1+1+
     *     ,snoage,fmp_com,flag_dsws,ltropo            ! 3+1+.5+.5+
     *     ,atmlnd%fr_snow_rad,dlake,flake                   ! 2+1+1
C****   output data: really needed only if kradia=2
     *     ,srhra,trhra                                ! 2(LM+LM_REQ+1)}
C****   total: dimrad_sv= IM*JM*(7*LM + 3*LM_REQ + 24 (+4)) => RAD_COM.f
     *     ,iy
          if (qcheck) then
            write(out_line,*) 'reading RADfile at Itime',Itime,it,iy
            call write_parallel(trim(out_line),unit=6)
          endif
        end do
        iend = 0
   10   if (it.ne.iy.or.iend.eq.1) then
          write(out_line,*) 'RAD input file bad or too short:',itime,
     *                       it,iy,iend
          call write_parallel(trim(out_line),unit=6)
          call stop_model('RADIA: input file bad or too short',255)
        end if
      end if

      IF (MODRD.NE.0) GO TO 900
      IDACC(ia_rad)=IDACC(ia_rad)+1
      moddrf=1    ! skip rad.forcing diags if nradfrc.le.0
      if (nradfrc>0) moddrf=mod(itime-itimei,nrad*nradfrc)
C****
      if (moddrf==0) IDACC(ia_rad_frc)=IDACC(ia_rad_frc)+1
C**** Interface with radiation routines, done only every NRAD time steps
C****
C**** Calculate mean cosine of zenith angle for the full radiation step
      ROT2=ROT1+TWOPI*NRAD*DTsrc/SECONDS_PER_DAY
      CALL COSZS (ROT1,ROT2,COSZ2,COSZA)
      JDAYR=dayOfYear
      JYEARR=YEAR

      if(is_set_param('s0')) then
        ! typically only used for SCM
        call get_param('s0',s0)
        s00wm2 = s0 ! just in case
      else
        S0=S0X*S00WM2*RATLS0/RSDIST
      endif

#ifdef OLD_BCdalbsn
c**** find scaling factors for surface albedo reduction
      if(dalbsnX.ne.0.) then
      IF (HAVE_SOUTH_POLE) THEN
         sumda_psum(:,1)=axyp(1,1)
         tauda_psum(:,1)=axyp(1,1)*depobc_1990(1,1)
      End If
      do j=J_0S,J_1S
      do i=I_0,I_1
c ilon72, jlat46 are indices w.r.t 72x46 grid
c      JLAT46=INT(1.+(J-1.)*0.25*DLAT_DG+.5)   ! slightly more general
c      ILON72=INT(.5+(I-.5)*72./IM+.5)
        ilon72 = 1 + int( 72d0*lon2d(i,j)/twopi )
        jlat46 = 1 + int( 45d0*(lat2d(i,j)+92d0*radian)/pi )
        fsnow = flice(i,j) + rsi(i,j)*(1-fland(i,j))
        if(atmlnd%SNOWE(I,J).gt.0.) fsnow = fsnow+fearth(i,j)
        sumda_psum(i,j) = axyp(i,j)*fsnow
        tauda_psum(i,j) = axyp(i,j)*fsnow*depobc_1990(i,j)
      end do
      end do
      IF (HAVE_NORTH_POLE) THEN
         sumda_psum(:,JM)=axyp(1,jm)*rsi(1,jm)
         tauda_psum(:,JM)=axyp(1,jm)*rsi(1,jm)*depobc_1990(1,jm)
      END IF
      CALL GLOBALSUM(grid, sumda_psum,sumda,all=.true.)
      CALL GLOBALSUM(grid, tauda_psum,tauda,all=.true.)

      xdalbs=-dalbsnX*sumda/tauda
      IF(QCHECK) write(6,*) 'coeff. for snow alb reduction',xdalbs
      endif ! dalbsnX not zero
#endif

      if(kradia.le.0) then
      IF (QCHECK) THEN
C****   Calculate mean strat water conc
        STRATQ=0.
        MSTRAT=0.
        DO J=J_0,J_1
          STRJ=0.
          MSTJ=0.
          DO I=I_0,IMAXJ(J)
            DO L=LTROPO(I,J)+1,LM
              STRJ = STRJ + Q(I,J,L)*MA(L,I,J)*AXYP(I,J)
              MSTJ = MSTJ + MA(L,I,J)*AXYP(I,J)
            END DO
          END DO
          IF (J.eq.1 .or. J.eq.JM) THEN
            STRJ=STRJ*IM
            MSTJ=MSTJ*IM
          END IF
          STRATQ=STRATQ+STRJ
          MSTRAT=MSTRAT+MSTJ
        END DO
        PRINT*,"Strat water vapour (ppmv), mass (mb)",1d6*STRATQ*mair
     *       /(18.*MSTRAT),PMTOP+1d-2*GRAV*MSTRAT/AREAG
      END IF

C**** Get the random numbers outside openMP parallel regions
C**** but keep MC calculation separate from SS clouds
C**** To get parallel consistency also with mpi, force each process
C**** to generate random numbers for all latitudes (using BURN_RANDOM)

C**** MC clouds are considered as a block for each I,J grid point

      CALL BURN_RANDOM(nij_before_j0(J_0))

      DO J=J_0,J_1                    ! complete overlap
      CALL BURN_RANDOM((I_0-1))
      DO I=I_0,IMAXJ(J)
        RDMC(I,J) = RANDU(X)          ! 1 random number per column
      END DO
      CALL BURN_RANDOM(nij_after_i1(I_1))
      END DO

      CALL BURN_RANDOM((nij_after_j1(J_1)))

C**** SS clouds are considered as a block for each continuous cloud
      CALL BURN_RANDOM(nij_before_j0(j_0)*LM)

      DO J=J_0,J_1                    ! semi-random overlap
      CALL BURN_RANDOM((I_0-1)*LM)
      DO I=I_0,IMAXJ(J)
        ! reverse loop kept only for consistency with previous version
        DO L=LM,1,-1   !   better:  1,LM
          IF(TAUSS(L,I,J).le.taulim) CLDSS(L,I,J)=0.
          IF(TAUMC(L,I,J).le.taulim) CLDMC(L,I,J)=0.
          RDSS(L,I,J) = RANDU(X)
        END DO
      END DO
      CALL BURN_RANDOM(nij_after_i1(I_1)*LM)
      END DO

      CALL BURN_RANDOM(nij_after_j1(j_1)*LM)

      end if                    ! kradia le 0

#if (defined ACCMIP_LIKE_DIAGS) 
#ifndef SKIP_ACCMIP_GHG_RADF_DIAGS
! because of additional updghg calls, these factors won''t apply:
      if(CO2X.ne.1.)  call stop_model('CO2x.ne.1 accmip diags',255)
      if(N2OX.ne.1.)  call stop_model('N2Ox.ne.1 accmip diags',255)
      if(CH4X.ne.1.)  call stop_model('CH4x.ne.1 accmip diags',255)
      if(CFC11X.ne.1.)call stop_model('CFC11x.ne.1 accmip diags',255)
      if(CFC12X.ne.1.)call stop_model('CFC12x.ne.1 accmip diags',255)
      if(XGHGX.ne.1.) call stop_model('XGHGx.ne.1 accmip diags',255)
      GFrefY=1850; GFrefD=182     ! ghg forcing refrnce year, day
      GFnowY=JyearR; GFnowD=JdayR ! ghg current desired year, day
      if(KJDAYG > 0) GFnowD=KJDAYG ! unless presribed in deck
      if(KYEARG > 0) GFnowY=KYEARG !
      call updghg(GFrefY,GFrefD)
      sv_fulgas_ref(1:4)=fulgas(nfghg(1:4))
      call updghg(GFnowY,GFnowD)
      sv_fulgas_now(1:4)=fulgas(nfghg(1:4))
#endif
#endif

! Set variables used in storing reference CO2 for uncoupled runs
#ifdef GCC_UNCOUPLE_RAD_CONCEN
      if(KJDAYG > 0) GCCco2nowD=KJDAYG ! unless presribed in deck
      if(KYEARG > 0) GCCco2nowY=KYEARG !
      call updghg(1850,182)
      GCCco2_fulgas_ref=fulgas(2)
      call updghg(GCCco2nowD,GCCco2nowD)
      GCCco2_fulgas_now=fulgas(2)
#endif

#ifdef HEALY_LM_DIAGS
      FCO2=FULGAS(2)*CO2R
      FN2O=FULGAS(6)*N2OR
      FCH4=FULGAS(7)*CH4R
c
c      write(6,*) 'RJH: GHG: CONC=',
c     * FCO2,FN2O,FCH4
      ghg_totforc=5.35d0*log(FCO2/CO2I)
     *  +0.036d0*(sqrt(FCH4)-sqrt(CH4I))
     *  -(Fe(FCH4,N2OI)-Fe(CH4I,N2OI))
     *  + 0.12d0*(sqrt(FN2O)-sqrt(N2OI))
     *  -(Fe(CH4I,FN2O)-Fe(CH4I,N2OI))
c      write(6,*) 'RJH: GHG: FORC=',ghg_totforc
#endif

      aj_alb_inds = (/ J_PLAVIS, J_PLANIR, J_ALBVIS, J_ALBNIR,
     &                 J_SRRVIS, J_SRRNIR, J_SRAVIS, J_SRANIR /)

      cfrac = 0.
      wtrcld = 0.
      icecld = 0.
      tausumw = 0.
      tausumi = 0.
      ctp = 0.
      ctt = 0.
      swus = 0.
      wtrtau = 0.
      icetau = 0.
#ifdef CFMIP3_SUBDD
      swut = 0.
      swutcs = 0.
      cfmip_twp = 0.
      swdcls = 0.
      swucls = 0.
      swdt = 0.
      cfmip_cf = 0.
      cfmip_qci = 0.
      cfmip_qcl = 0.
#endif
C****
C**** MAIN J LOOP
C****
      DO J=J_0,J_1

c     ICKERR=0
c     JCKERR=0
c     KCKERR=0

C****
C**** MAIN I LOOP
C****
      DO I=I_0,IMAXJ(J)
C**** Radiation input files use a 72x46 grid independent of IM and JM
C**** (ilon72,jlat46) is the 4x5 box containing the center of box (i,j)
c      JLAT46=INT(1.+(J-1.)*45./(JM-1.)+.5)  !  lat_index w.r.to 72x46 grid
c      JLAT46=INT(1.+(J-1.)*0.25*DLAT_DG+.5) ! slightly more general
c      ILON72=INT(.5+(I-.5)*72./IM+.5)  ! lon_index w.r.to 72x46 grid
      igcm = i
      jgcm = j
      ilon72 = 1 + int( 72d0*lon2d(i,j)/twopi )
      jlat46 = 1 + int( 45d0*(lat2d(i,j)+92d0*radian)/pi )
#ifdef ALTER_RADF_BY_LAT
      FULGAS(:)=FULGAS_orig(:)*FULGAS_lat(:,JLAT46)
      FS8OPX(:)=FS8OPX_orig(:)*FS8OPX_lat(:,JLAT46)
      FT8OPX(:)=FT8OPX_orig(:)*FT8OPX_lat(:,JLAT46)
#endif
      L1 = 1                         ! lowest layer above ground
      LMR=LM+LM_REQ                  ! radiation allows var. # of layers
      JR=JREG(I,J)
C**** DETERMINE FRACTIONS FOR SURFACE TYPES AND COLUMN PRESSURE
      PLAND=FLAND(I,J)
      POICE=RSI(I,J)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      PLAKE=FLAKE(I,J)
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
      ptype4(1) = pocean ! open ocean and open lake
      ptype4(2) = poice  ! ocean/lake ice
      ptype4(3) = plice  ! glacial ice
      ptype4(4) = pearth ! non glacial ice covered soil

C**** CHECK SURFACE TEMPERATURES
      DO IT=1,4
        IF(ptype4(IT) > 0.) then
          IF(int(asflx4(it)%GTEMPR(I,J)) .LT. planck_tmin .OR.
     &       int(asflx4(it)%GTEMPR(I,J)) .GE. planck_tmax ) then
            WRITE(6,*) 'In Radia: Time,I,J,IT,TG1',ITime,I,J,IT
     *         ,asflx4(it)%GTEMPR(I,J)
CCC         STOP 'In Radia: Grnd Temp out of range'
c           ICKERR=ICKERR+1
          END IF
        END IF
      END DO

C**** Set Chlorophyll concentration
      if (POCEAN.gt.0) then
        if( (chl_from_seawifs>0 .or. chl_from_obio>0)
     &       .and. atmocn%chl_defined ) then
          LOC_CHL = atmocn%chl(I,J)
          if (ij_chl.gt.0)
     .       AIJ(I,J,IJ_CHL)=AIJ(I,J,IJ_CHL)+atmocn%CHL(I,J)*FOCEAN(I,J)
!         write(*,'(a,3i5,e12.4)')'RAD_DRV:',
!    .    itime,i,j,chl(i,j)
        else
          LOC_CHL = -1.d30
        endif
      endif

      LS1_loc=LTROPO(I,J)+1  ! define stratosphere for radiation
C**** kradia>1: adjusted forcing, i.e. T adjusts in L=LS1_loc->LM+3
      if(kradia>1) LS1_loc=LS1_loc+2-kradia     ! favorite:kradia=3
      if(kradia>3) LS1_loc=1                    ! favorite:kradia=3
      kdeliq=0   ! initialize mainly for L>LM
      if (kradia.gt.0) then     ! rad forcing model
        do l=1,lm
          tlm(l) = T(i,j,l)*pk(l,i,j)
          shl(l) = QR(l,i,j) ; if(shl(l)<0) shl(l)=0
          tauwc(l) = cldx*CLDinfo(l,1,i,j)
          tauic(l) = cldx*CLDinfo(l,2,i,j)
          SIZEWC(L)= CLDinfo(l,3,i,j)
          SIZEIC(L)= SIZEWC(L)
        end do
      else                      ! full model
C****
C**** DETERMINE CLOUDS (AND THEIR OPTICAL DEPTHS) SEEN BY RADIATION
C****
      CSS=0. ; CMC=0. ; CLDCV=0. ; DEPTH=0. ; OPTDW=0. ; OPTDI=0.
      if(cc_cdncx.ne.0. .or. od_cdncx.ne.0.) then
        call dCDNC_EST(i,j,pland, dCDNC)
      else
        dCDNC = 0.
      endif
      dCC_CDNCL = CC_cdncx*dCDNC*CDNCL
      dOD_CDNCL = OD_cdncx*dCDNC*CDNCL

C**** Adjust RDSS for semi-random overlap
      call get_cld_overlap (lm, cldss(:,i,j), randSS=rdss(:,i,j))

      DO L=1,LM
        if(q(i,j,l)<0) then
           WRITE(6,*)'In Radia: Time,I,J,L,Q<0',ITime,I,J,L,Q,'->0'
           Q(I,J,L)=0.
        end if
        QSS=Q(I,J,L)/(RHSAV(L,I,J)+1.D-20)
        shl(L)=QSS
        IF(FSS(L,I,J)*CLDSAV(L,I,J).LT.1.)
     *       shl(L)=(Q(I,J,L)-QSS*FSS(L,I,J)*CLDSAV(L,I,J))/
     /              (1.-FSS(L,I,J)*CLDSAV(L,I,J))
        TLm(L)=T(I,J,L)*PK(L,I,J)
        rhodz=pdsig(l,i,j)*100/grav
        TAUSSL=0.
        TAUSSLIP=0.
        TAUMCL=0.
        TAUWC(L)=0.
        TAUIC(L)=0.
        SIZEWC(L)=0.
        SIZEIC(L)=0.
        TOTCLD(L)=0.
        SS_CLD(L)=0.
C**** Determine large scale and moist convective cloud cover for radia
        IF (CLDSS(L,I,J)*(1.+dcc_cdncl(l)).GT.RDSS(L,I,J)) THEN
          TAUSSL=TAUSS(L,I,J)*(1.+dod_cdncl(l))
          ! tausslip is tau of ice precip in a supercooled water cloud
          TAUSSLIP=TAUSSIP(L,I,J)*(1.+dod_cdncl(l))
          shl(L)=QSS
          CSS=1.
          call inc_ajl(i,j,l,jl_sscld,css)
#ifdef CFMIP3_SUBDD
          ! LS Cloud
          cfmip_cf(i,j,l)=cfmip_cf(i,j,l)+1.
#endif
        END IF
        IF (CLDMC(L,I,J).GT.RDMC(I,J)) THEN
          CMC=1.
          call inc_ajl(i,j,l,jl_mccld,cmc)
#ifdef CFMIP3_SUBDD
          ! MC Cloud
          cfmip_cf(i,j,l)=min(cfmip_cf(i,j,l)+1.,1.)
#endif
          DEPTH=DEPTH+PDSIG(L,I,J)
          IF(TAUMC(L,I,J).GT.TAUSSL+TAUSSLIP) THEN
            TAUMCL=TAUMC(L,I,J)
            ELHX=LHE
            IF(TLm(L).LE.TF) ELHX=LHS
            shl(L)=QSAT(TLm(L),ELHX,PMID(L,I,J))
          END IF
        END IF
        IF(TAUSSL+TAUSSLIP+TAUMCL.GT.0.) THEN
             CLDCV=1.
          TOTCLD(L)=1.
          call inc_ajl(i,j,l,jl_totcld,1d0)
C**** save 3D cloud fraction as seen by radiation
          if(cldx>0) AIJL(I,J,L,IJL_CF)=AIJL(I,J,L,IJL_CF)+1.
          IF(TAUMCL.GT.TAUSSL+TAUSSLIP) THEN
            SIZEWC(L)=CSIZMC(L,I,J)
            SIZEIC(L)=CSIZMC(L,I,J)
            IF(SVLAT(L,I,J).EQ.LHE) THEN
              TAUWC(L)=cldx*TAUMCL
              OPTDW=OPTDW+TAUWC(L)
              call inc_ajl(i,j,l,jl_wcld,1d0)
              call inc_ajl(i,j,l,jl_wcldwt,pdsig(l,i,j))
              aij(i,j,ij_lwprad)=aij(i,j,ij_lwprad)+QLmc(l,i,j)*rhodz
     &                          /cldmc(l,i,j)
              aijl(i,j,l,ijl_QLrad)=aijl(i,j,l,ijl_QLrad)
     &                             +QLmc(l,i,j)*pdsig(l,i,j)
     &                             /cldmc(l,i,j)
#ifdef CFMIP3_SUBDD
              ! MC Cloud Liquid
              cfmip_twp(i,j)=cfmip_twp(i,j)
     &                       +QLmc(l,i,j)*rhodz/cldmc(l,i,j)
              cfmip_qcl(i,j,l)=QLmc(l,i,j)/cldmc(l,i,j)
#endif
            ELSE
              TAUIC(L)=cldx*TAUMCL
              OPTDI=OPTDI+TAUIC(L)
              call inc_ajl(i,j,l,jl_icld,1d0)
              call inc_ajl(i,j,l,jl_icldwt,pdsig(l,i,j))
              aij(i,j,ij_iwprad)=aij(i,j,ij_iwprad)+QImc(l,i,j)*rhodz
     &                          /cldmc(l,i,j)
              aijl(i,j,l,ijl_QIrad)=aijl(i,j,l,ijl_QIrad)
     &                             +QImc(l,i,j)*pdsig(l,i,j)
     &                             /cldmc(l,i,j)
#ifdef CFMIP3_SUBDD
              ! MC Cloud Ice
              cfmip_twp(i,j)=cfmip_twp(i,j)
     &                       +QImc(l,i,j)*rhodz/cldmc(l,i,j)
              cfmip_qci(i,j,l)=QImc(l,i,j)/cldmc(l,i,j)
#endif
            END IF
          ELSE
            SS_CLD(L)=1.
            SIZEWC(L)=CSIZSS(L,I,J)
            SIZEIC(L)=CSIZSS(L,I,J)
            IF(SVLHX(L,I,J).EQ.LHE) THEN
              TAUWC(L)=cldx*TAUSSL
              OPTDW=OPTDW+TAUWC(L)
              call inc_ajl(i,j,l,jl_wcld,1d0)
              call inc_ajl(i,j,l,jl_wcldwt,pdsig(l,i,j))
              aij(i,j,ij_lwprad)=aij(i,j,ij_lwprad)+QLss(l,i,j)*rhodz
     &                          /cldss(l,i,j)
              aijl(i,j,l,ijl_QLrad)=aijl(i,j,l,ijl_QLrad)
     &                             +QLss(l,i,j)*pdsig(l,i,j)
     &                             /cldss(l,i,j)
#ifdef CFMIP3_SUBDD
              ! LS Cloud Liquid
              cfmip_twp(i,j)=cfmip_twp(i,j)
     &                       +QLss(l,i,j)*rhodz/cldss(l,i,j)
              cfmip_qcl(i,j,l)=QLss(l,i,j)/cldss(l,i,j)
#endif
              if(tausslip.gt.0.) then
                SIZEIC(L)=CSIZSSIP(L,I,J)
                TAUIC(L)=cldx*TAUSSLIP
                OPTDI=OPTDI+TAUIC(L)
                call inc_ajl(i,j,l,jl_icld,1d0)
                call inc_ajl(i,j,l,jl_icldwt,pdsig(l,i,j))
                aij(i,j,ij_iwprad)=aij(i,j,ij_iwprad)+QIss(l,i,j)*rhodz
     &                            /cldss(l,i,j)
                aijl(i,j,l,ijl_QIrad)=aijl(i,j,l,ijl_QIrad)
     &                               +QIss(l,i,j)*pdsig(l,i,j)
     &                               /cldss(l,i,j)
#ifdef CFMIP3_SUBDD
               ! LS Snow in supercooled liquid
              cfmip_twp(i,j)=cfmip_twp(i,j)
     &                       +QIss(l,i,j)*rhodz/cldss(l,i,j)
#endif
              endif
            ELSE
              TAUIC(L)=cldx*TAUSSL
              OPTDI=OPTDI+TAUIC(L)
              call inc_ajl(i,j,l,jl_icld,1d0)
              call inc_ajl(i,j,l,jl_icldwt,pdsig(l,i,j))
              aij(i,j,ij_iwprad)=aij(i,j,ij_iwprad)+QIss(l,i,j)*rhodz
     &                          /cldss(l,i,j)
              aijl(i,j,l,ijl_QIrad)=aijl(i,j,l,ijl_QIrad)
     &                             +QIss(l,i,j)*pdsig(l,i,j)
     &                             /cldss(l,i,j)
#ifdef CFMIP3_SUBDD
              ! LS Cloud Ice
              cfmip_twp(i,j)=cfmip_twp(i,j)
     &                       +QIss(l,i,j)*rhodz/cldss(l,i,j)
              cfmip_qci(i,j,l)=QIss(l,i,j)/cldss(l,i,j)
#endif
            END IF
          END IF
          call inc_ajl(i,j,l,jl_wcod,tauwc(l))
          call inc_ajl(i,j,l,jl_icod,tauic(l))
          call inc_ajl(i,j,l,jl_wcsiz,sizewc(l)*tauwc(l))
          call inc_ajl(i,j,l,jl_icsiz,sizeic(l)*tauic(l))
          aijl(i,j,l,ijl_wtrtau)=aijl(i,j,l,ijl_wtrtau)+tauwc(l)
          aijl(i,j,l,ijl_icetau)=aijl(i,j,l,ijl_icetau)+tauic(l)
          wtrtau(i,j,l)=tauwc(l)
          icetau(i,j,l)=tauic(l)
        END IF
C**** save some radiation/cloud fields for wider use
        RCLD(L,I,J)=TAUWC(L)+TAUIC(L)
      END DO
      CFRAC(I,J) = CLDCV    ! cloud fraction consistent with radiation
C**** effective cloud cover diagnostics
         OPNSKY=1.-CLDCV
         DO IT=1,NTYPE
           call inc_aj(i,j,it,J_PCLDSS,CSS  *FTYPE(IT,I,J))
           call inc_aj(i,j,it,J_PCLDMC,CMC  *FTYPE(IT,I,J))
           call inc_aj(i,j,it,J_CLDDEP,DEPTH*FTYPE(IT,I,J))
           call inc_aj(i,j,it,J_PCLD  ,CLDCV*FTYPE(IT,I,J))
         END DO
         call inc_areg(i,j,jr,J_PCLDSS,CSS  )
         call inc_areg(i,j,jr,J_PCLDMC,CMC  )
         call inc_areg(i,j,jr,J_CLDDEP,DEPTH)
         call inc_areg(i,j,jr,J_PCLD  ,CLDCV)
         AIJ(I,J,IJ_PMCCLD)=AIJ(I,J,IJ_PMCCLD)+CMC
         AIJ(I,J,IJ_CLDCV) =AIJ(I,J,IJ_CLDCV) +CLDCV
         DO L=1,LLOW
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDL)=AIJ(I,J,IJ_PCLDL)+1.
           exit
         end do
         DO L=LLOW+1,LMID
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDM)=AIJ(I,J,IJ_PCLDM)+1.
           exit
         end do
         DO L=LMID+1,LHI
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDH)=AIJ(I,J,IJ_PCLDH)+1.
           exit
         end do
         DO L=1,LLOW
           IF (SS_CLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDL_SS)=AIJ(I,J,IJ_PCLDL_SS)+1.
           exit
         end do

         TAUSUMW(I,J) = OPTDW
         TAUSUMI(I,J) = OPTDI
         if(optdw.gt.0.) then
            AIJ(I,J,IJ_optdw)=AIJ(I,J,IJ_optdw)+optdw
            AIJ(I,J,IJ_wtrcld)=AIJ(I,J,IJ_wtrcld)+1.
            WTRCLD(I,J) =   1.
         end if
         if(optdi.gt.0.) then
            AIJ(I,J,IJ_optdi)=AIJ(I,J,IJ_optdi)+optdi
            AIJ(I,J,IJ_icecld)=AIJ(I,J,IJ_icecld)+1.
            ICECLD(I,J) =   1.
         end if

         DO KR=1,NDIUPT
           IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
C**** Warning: this replication may give inaccurate results for hours
C****          1->(NRAD-1)*DTsrc (ADIURN) or skip them (HDIURN)
             TMP(IDD_CL7:IDD_CL7+6)=TOTCLD(1:7)
             TMP(IDD_CCV)=CLDCV
             DO INCH=1,NRAD
               IHM=1+(JTIME+INCH-1)*HR_IN_DAY/NDAY
               IH=IHM
               IF(IH.GT.HR_IN_DAY) IH = IH - HR_IN_DAY
               ADIURN(IDX(:),KR,IH)=ADIURN(IDX(:),KR,IH)+TMP(IDX(:))
#ifdef USE_HDIURN
               IHM = IHM+(DATE-1)*HR_IN_DAY
               IF(IHM.LE.HR_IN_MONTH) THEN
                 HDIURN(IDX(:),KR,IHM)=HDIURN(IDX(:),KR,IHM)+TMP(IDX(:))
               ENDIF
#endif
             ENDDO
           END IF
         END DO
      end if ! kradia le 0 (full model)
C****
C**** SET UP VERTICAL ARRAYS OMITTING THE I AND J INDICES
C****
C**** EVEN PRESSURES
#ifdef TRACERS_TOMAS
      aesqex(:,:,:)=0.0
      aesqsc(:,:,:)=0.0
      aesqcb(:,:,:)=0.0
      aesqex_dry(:,:,:)=0.0
      aesqsc_dry(:,:,:)=0.0
      aesqcb_dry(:,:,:)=0.0
#endif
      PLB(LM+1)=PEDN(LM+1,I,J)
      DO L=1,LM
        PLB(L)=PEDN(L,I,J)
C**** TEMPERATURES
C---- TLm(L)=T(I,J,L)*PK(L,I,J)     ! already defined
        IF(int(TLm(L)) .LT. planck_tmin .OR.
     &     int(TLm(L)) .GE. planck_tmax ) THEN
          WRITE(6,*) 'In Radia: Time,I,J,L,TL',ITime,I,J,L,TLm(L)
          WRITE(6,*) 'GTEMPR:',
     &         asflx4(1)%GTEMPR(I,J),asflx4(2)%GTEMPR(I,J),
     &         asflx4(3)%GTEMPR(I,J),asflx4(4)%GTEMPR(I,J)
CCC       STOP 'In Radia: Temperature out of range'
c         ICKERR=ICKERR+1
        END IF
C**** MOISTURE VARIABLES
C---- shl(L)=Q(I,J,L)        ! already defined and reset to 0 if <0
c       if(shl(l).lt.0.) then
c         WRITE(0,*)'In Radia: Time,I,J,L,QL<0',ITime,I,J,L,shl(L),'->0'
c         KCKERR=KCKERR+1
c         shl(l)=0.
c       end if
        RHL(L) = shl(L)/QSAT(TLm(L),LHE,PMID(L,I,J))
        if(RHfix.ge.0.) RHL(L)=RHfix
C**** Extra aerosol data
C**** For up to nraero_aod aerosols, define the aerosol amount to
C**** be used (kg/m^2)
C**** Only define TRACER if individual tracer is actually defined.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_AEROSOLS_SEASALT)
C**** loop over tracers that are passed to radiation.
C**** Some special cases for black carbon, organic carbon, SOAs where
C**** more than one tracer is lumped together for radiation purposes
      do n=1,nraero_aod
        select case (trname(ntrix_aod(n)))
          case ("OCIA", "vbsAm2")
            TRACER(L,n)=(
#ifdef TRACERS_AEROSOLS_VBS
     *           sum(trm(i,j,l,vbs_tr%iaer))
#else
     *           trm(i,j,l,n_OCII)+trm(i,j,l,n_OCIA)
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_OCEAN
     *          +trm(i,j,l,n_ococean)
#endif  /* TRACERS_AEROSOLS_OCEAN */
     *                  )*BYAXYP(I,J)
          case ("OCB")
#ifdef TRACERS_AEROSOLS_VBS
           TRACER(L,n)=0.d0
#else
           TRACER(L,n)=trm(i,j,l,n_OCB)*BYAXYP(I,J)
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_SOA
          case ("isopp1a")
           TRACER(L,n)=( trm(i,j,l,n_isopp1a)+trm(i,j,l,n_isopp2a)
#ifdef TRACERS_TERP
     &                  +trm(i,j,l,n_apinp1a)+trm(i,j,l,n_apinp2a)
#endif /* TRACERS_TERP */
     &                 )*BYAXYP(I,J)
#endif /* TRACERS_AEROSOLS_SOA */
          case ("BCIA")
           TRACER(L,n)=(trm(i,j,l,n_BCII)+trm(i,j,l,n_BCIA))*BYAXYP(I,J)
          case default
#ifdef TRACERS_NITRATE
! assume full neutralization of NO3p, if NH4 suffice
          select case (trname(ntrix_aod(n)))
          case ("NO3p")
            if (trm(i,j,l,ntrix_aod(n)) > 0.d0) then
              nh4_on_no3=min(trm(i,j,l,n_NO3p)*(tr_mm(n_NO3p)+
     *            tr_mm(n_NH4))/tr_mm(n_NO3p)-trm(i,j,l,n_NO3p),
     *                       trm(i,j,l,n_NH4))
              wttr(n)=(nh4_on_no3+trm(i,j,l,ntrix_aod(n)))/
     *                trm(i,j,l,ntrix_aod(n))
            endif
          case ("SO4")
            if (trm(i,j,l,ntrix_aod(n)) > 0.d0) then
              nh4_on_no3=min(trm(i,j,l,n_NO3p)*(tr_mm(n_NO3p)+
     *            tr_mm(n_NH4))/tr_mm(n_NO3p)-trm(i,j,l,n_NO3p),
     *                       trm(i,j,l,n_NH4))
              wttr(n)=(trm(i,j,l,n_NH4)-nh4_on_no3+
     *                 trm(i,j,l,ntrix_aod(n)))/trm(i,j,l,ntrix_aod(n))
            endif
          end select
#endif
          TRACER(L,n)=wttr(n)*trm(i,j,l,ntrix_aod(n))*BYAXYP(I,J)
        end select
      end do
#endif /* TRACERS_AEROSOLS_Koch/DUST/MINERALS/SEASALT */

#ifdef TRACERS_AMP
      CALL SETAMP_LEV(i,j,l)
#endif
#ifdef TRACERS_TOMAS
      CALL SETTOMAS_LEV(i,j,l)
#endif
      END DO
C**** Radiative Equilibrium Layer data
      DO K=1,LM_REQ
        IF(int(RQT(K,I,J)) .LT. planck_tmin .OR.
     &     int(RQT(K,I,J)) .GE. planck_tmax ) THEN
        WRITE(6,*) 'In RADIA: Time,I,J,L,TL',ITime,I,J,LM+K,RQT(K,I,J)
CCC     STOP 'In Radia: RQT out of range'
c       JCKERR=JCKERR+1
        END IF
        TLm(LM+K)=RQT(K,I,J)
        PLB(LM+k+1) = PLB0(k)
        shl(LM+k)    = shl0(k)
        RHL(LM+k) = shl(LM+k)/QSAT(TLm(LM+k),LHE,
     *                                .5d0*(PLB(LM+k)+PLB(LM+k+1)) )
        tauwc(LM+k) = 0.
        tauic(LM+k) = 0.
        sizewc(LM+k)= 0.
        sizeic(LM+k)= 0.
#ifdef TRACERS_ON
C**** set radiative equilibrium extra tracer amount to zero
        IF (nraero_aod.gt.0) TRACER(LM+k,1:nraero_aod)=0.
#endif
      END DO
      if (kradia.gt.1) then
        do l=1,lm+lm_req
          tlm(l) = tlm(l) + Tchg(l,i,j)
          AFLX_ST(L,I,J,5)=AFLX_ST(L,I,J,5)+Tchg(L,I,J)
        end do
      end if
C**** Zenith angle and GROUND/SURFACE parameters
      COSZ=COSZA(I,J)
      TGO =atmocn%GTEMPR(I,J)
      TGOI=atmice%GTEMPR(I,J)
      TGLI=atmgla%GTEMPR(I,J)
      TGE =atmlnd%GTEMPR(I,J)
      TSL=atmsrf%TSAVG(I,J)
      SNOWOI=SNOWI(I,J)
      SNOWLI=atmgla%SNOW(I,J)
      !SNOWE=atmlnd%SNOWE(I,J)                    ! snow depth (kg/m**2)
      SNOWD(:)=snowd_ij(:,I,J)
      snow_frac(:) = atmlnd%fr_snow_rad(:,i,j)    ! snow cover (1)
      AGESN(1)=SNOAGE(3,I,J)    ! land         ! ? why are these numbers
      AGESN(2)=SNOAGE(1,I,J)    ! ocean ice        so confusing ?
      AGESN(3)=SNOAGE(2,I,J)    ! land ice
c      print*,"snowage",i,j,SNOAGE(1,I,J)
C**** set up parameters for new sea ice and snow albedo
      zsnwoi=atmice%ZSNOWI(I,J)
      if(dalbsnX.ne.0.) then
#ifdef OLD_BCdalbsn
        dALBsn = xdalbs*depobc(i,j)
#else
        dALBsn = dalbsnX*BCdalbsn(i,j)
#endif
      else
        dALBsn = 0.
      endif

c to use on-line tracer albedo impact, set dALBsnX=0. in rundeck
#ifdef BC_ALB
      call GET_BC_DALBEDO(i,j,dALBsn1,bc_snow_present(i,j))
      if (rad_interact_aer > 0) dALBsn=dALBsn1
      dALBsnBC(I,J)=dALBsn1
#endif  /* BC_ALB */
      if (poice.gt.0.) then
        zoice = ZSI(i,j)
        flags=flag_dsws(i,j)
        if (kradia .le. 0) then
          fmp=min(1.6d0*sqrt(pond_melt(i,j)/rhow),1d0)
             AIJ(I,J,IJ_FRMP) = AIJ(I,J,IJ_FRMP) + fmp*POICE
        else
          fmp = fmp_com(i,j)
        end if
        zmp=min(0.8d0*fmp,0.9d0*zoice)
      else
        zoice=0. ; flags=.FALSE. ; fmp=0. ; zmp=0.
      endif
C**** set up new lake depth parameter to incr. albedo for shallow lakes
c      zlake=0.
c      if (plake.gt.0) then
c        zlake = MWL(I,J)/(RHOW*PLAKE*AXYP(I,J))
c      end if
      zlake = dlake(i,j)
C****
      if (kradia .le. 0) then
        !WEARTH=(WEARTH_COM(I,J)+AIEARTH(I,J))/(WFCS(I,J)+1.D-20)
        WEARTH=atmlnd%bare_soil_wetness(i,j)
        if (wearth.gt.1.) wearth=1.
      else                            ! rad.frc. model
        wearth = wsoil(i,j)
      end if
      if ( fearth(i,j) > 0.d0 ) then
        call ent_get_exports( entcells(i,j),
     &       vegetation_fractions=PVT0,
     &       vegetation_heights=HVT0 )
        call map_ent2giss(PVT0,HVT0,PVT) !temp hack: ent pfts->giss veg
      else
        PVT(:) = 0.d0  ! actually PVT is not supposed to be used in this case
      endif
      WMAG=atmsrf%WSAVG(I,J)
C****
C**** Radiative interaction and forcing diagnostics:
C**** If no radiatively active tracers are defined, nothing changes.
C**** Currently this works for aerosols and ozone but should be extended
C**** to cope with all trace gases.
C****
      FTAUC=1. ! deflt (clouds on)
      use_tracer_chem(:) = 0 ! by default use climatological ozone/ch4/co2
C**** Set level for inst. rad. forc. calcs for aerosols/trace gases
C**** This is set from the rundeck.
      LFRC=LM+LM_REQ+1          ! TOA
      if (rad_forc_lev.gt.0) LFRC=LTROPO(I,J) ! TROPOPAUSE
#ifdef ACCMIP_LIKE_DIAGS
      if(rad_forc_lev.gt.0)call stop_model
     &('ACCMIP_LIKE_DIAGS desires TOA RADF diags',255)
#endif
C**** The calculation of the forcing is slightly different.
C**** depending on whether full radiative interaction is turned on
C**** or not.
      onoff_aer=0; onoff_chem=0
      if (rad_interact_aer > 0) onoff_aer=1
      if (clim_interact_chem > 0) onoff_chem=1
      use_o3_ref=0

#ifdef TRACERS_SPECIAL_Shindell
C**** Ozone and Methane:
      CHEM_IN(1,1:LM)=chem_tracer_save(1,1:LM,I,J)
      CHEM_IN(2,1:LM)=chem_tracer_save(2,1:LM,I,J)*CH4X_RADoverCHEM
      if (clim_interact_chem > 0) then
        use_tracer_chem(1)=Lmax_rad_O3  ! O3
        use_tracer_chem(2)=Lmax_rad_CH4 ! CH4
      endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      if(clim_interact_chem<=0)
     &call stop_model("stratOx RADF on, clim_interact_chem<=0",255)
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */

C**** CO2
C Update GCCco2_IN and GCCco2_tracer_save with trm to pass on CO2
C information to RADIATION
#ifdef GCC_COUPLE_RAD
      do L=1,LM
         GCCco2_tracer_save(L,i,j)=(trm(i,j,L,n_CO2n))
     &   * byaxyp(i,j)*avog/(tr_mm(n_CO2n)*2.69e20)
      enddo
      GCCco2_IN(1:LM)=GCCco2_tracer_save(1:LM,I,J)*CO2X
      use_tracer_GCCco2=Lmax_rad_CO2 ! CO2
#endif /* GCC_COUPLE_RAD */

      if (moddrf==0) then
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_AEROSOLS_SEASALT) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
C**** Aerosols (OMA, MATRIX, TOMAS):
        do n=1,nraero_rf
          IF (trname(ntrix_rf(n)).eq."seasalt2") CYCLE ! not for seasalt2
          if (diag_fc==2) then
            FSTOPX(n) = 1-onoff_aer !turns off online tracer
            FTTOPX(n) = 1-onoff_aer !
C**** Warning: small bit of hardcoding assumes that seasalt2 immediately
C****          succeeds seasalt1 in nraero_rf array
            IF (trname(ntrix_rf(n)).eq."seasalt1") THEN          !add seasalt2
              FSTOPX(n+1)=1-onoff_aer;FTTOPX(n+1)=1-onoff_aer !to seasalt1
            END IF
          else if (diag_fc==1) then
            FSTOPX(1:nraero_aod) = 1-onoff_aer !turns off online tracer
            FTTOPX(1:nraero_aod) = 1-onoff_aer !
          endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_AEROSOLS_SEASALT)
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
#endif
          CALL RCOMPX  ! tr.aero.Koch/dust/miner./seasalt
          SNFST(1,n,I,J)=SRNFLB(1) ! surface forcing
          TNFST(1,n,I,J)=TRNFLB(1)
          SNFST(2,n,I,J)=SRNFLB(LFRC) ! Tropopause forcing
          TNFST(2,n,I,J)=TRNFLB(LFRC)
          if (diag_fc==2) then
            FSTOPX(n) = onoff_aer !turns on online tracer
            FTTOPX(n) = onoff_aer !
            IF (trname(ntrix_rf(n)).eq."seasalt1") THEN    ! also for seasalt2
              FSTOPX(n+1)=onoff_aer ; FTTOPX(n+1)=onoff_aer
            END IF
          else if (diag_fc==1) then
            FSTOPX(1:nraero_aod) = onoff_aer !turns on online tracer
            FTTOPX(1:nraero_aod) = onoff_aer !
          endif
        end do
#endif

#ifdef TRACERS_SPECIAL_Shindell
C**** Ozone:
! ozone rad forcing diags now use a constant reference year
! for this first call. And no tracer values...
        use_o3_ref=1 ; use_tracer_chem(1)=0
        kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
        CALL RCOMPX        ! tr_Shindell Ox tracer
        SNFST_o3ref(1,I,J)=SRNFLB(LTROPO(I,J)) ! meteorological tropopause
        TNFST_o3ref(1,I,J)=TRNFLB(LTROPO(I,J))
        SNFST_o3ref(2,I,J)=SRNFLB(LM+LM_REQ+1) ! T.O.A.
        TNFST_o3ref(2,I,J)=TRNFLB(LM+LM_REQ+1)
        SNFST_o3ref(5,I,J)=SRNFLB(LS1-1) ! fixed tropopause
        TNFST_o3ref(5,I,J)=TRNFLB(LS1-1)
#ifdef AUXILIARY_OX_RADF
! if needed, also save the auxiliary ozone field (i.e. climatology
! if tracer is used in final call, tracers if climatology is used.)
#ifdef AUX_OX_RADF_TROP
        ! forces use of tracer from L=1,LS1-1 and reference above that:
        use_o3_ref=1 ; use_tracer_chem(1)=LS1-1
#else
        ! use tracer or climatology, whichever won''t be used in final call:
        use_o3_ref=0 ; use_tracer_chem(1)=(1-onoff_chem)*Lmax_rad_O3
#endif
        kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
        CALL RCOMPX        ! tr_Shindell Ox tracer
#ifdef AUX_OX_RADF_TROP
        SNFST_o3ref(3,I,J)=SRNFLB(LS1-1)       ! fixed tropopause
        TNFST_o3ref(3,I,J)=TRNFLB(LS1-1)
#else
        SNFST_o3ref(3,I,J)=SRNFLB(LTROPO(I,J)) ! meteorological tropopause
        TNFST_o3ref(3,I,J)=TRNFLB(LTROPO(I,J))
#endif
        SNFST_o3ref(4,I,J)=SRNFLB(LM+LM_REQ+1) ! T.O.A.
        TNFST_o3ref(4,I,J)=TRNFLB(LM+LM_REQ+1)
#endif /* AUXILIARY_OX_RADF */
! After AUX call, use either climatological or tracer O3:
        use_o3_ref=0 ; use_tracer_chem(1)=onoff_chem*Lmax_rad_O3
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
! Optional intermediate call with stratOx tracer:
!NEED chem_IN(1,1:LM)=stratO3_tracer_save(1:LM,I,J)
!NEED kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
!NEED CALL RCOMPX        ! stratOx diag tracer
        SNFST_stratOx(1,I,J)=SRNFLB(LTROPO(I,J)) ! tropopause
        TNFST_stratOx(1,I,J)=TRNFLB(LTROPO(I,J))
        SNFST_stratOx(2,I,J)=SRNFLB(LM+LM_REQ+1) ! T.O.A.
        TNFST_stratOx(2,I,J)=TRNFLB(LM+LM_REQ+1)
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
        chem_IN(1,1:LM)=chem_tracer_save(1,1:LM,I,J)  ! Ozone
        chem_IN(2,1:LM)=chem_tracer_save(2,1:LM,I,J)*CH4X_RADoverCHEM  ! Methane
#if (defined ACCMIP_LIKE_DIAGS) 
#ifndef SKIP_ACCMIP_GHG_RADF_DIAGS
! TOA GHG rad forcing: nf=1,4 are CH4, N2O, CFC11, CFC12:
! Initial calls are reference year/day:
        do nf=1,4
          if(nf==1)then ! CH4 reference call must not use tracer
            use_tracer_chem(2)=0
          else ! N2O and CFC call's CH4 should match final call
            use_tracer_chem(2)=onoff_chem*Lmax_rad_CH4
          endif
          fulgas(nfghg(nf))=sv_fulgas_ref(nf)
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
          CALL RCOMPX
          SNFS_ghg(nf,I,J)=SRNFLB(LM+LM_REQ+1)
          TNFS_ghg(nf,I,J)=TRNFLB(LM+LM_REQ+1)
          fulgas(nfghg(nf))=sv_fulgas_now(nf)
        enddo
#endif /* NOT DEFINED SKIP_ACCMIP_GHG_RADF_DIAGS */
#endif /* ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */
      end if ! moddrf=0
#if (defined GCC_COUPLE_RAD)
! final (main) RCOMPX call can use tracer co2 (or not):
      use_tracer_GCCco2=Lmax_rad_CO2
      if (is_set_param('initial_GHG_setup')) then
        call get_param('initial_GHG_setup', initial_GHG_setup)
        if (initial_GHG_setup == 1 .and. itime == itimeI) then
          use_tracer_GCCco2=0  ! special case; model outputs climatology
        end if
      end if
#endif /* GCC_COUPLE_RAD */
#if (defined TRACERS_SPECIAL_Shindell)
! final (main) RCOMPX call can use tracer methane(or not):
      use_tracer_chem(2)=onoff_chem*Lmax_rad_CH4
      if (is_set_param('initial_GHG_setup')) then
        call get_param('initial_GHG_setup', initial_GHG_setup)
        if (initial_GHG_setup == 1 .and. itime == itimeI) then
          use_tracer_chem(2)=0  ! special case; model outputs climatology
        end if
      end if
#endif /* TRACERS_SPECIAL_Shindell */

#ifdef GCC_UNCOUPLE_RAD_CONCEN
! Use reference year CO2 for uncoupling radiation
      fulgas(2) = GCCco2_fulgas_ref  
#endif

      if (moddrf==0) then
#ifdef BC_ALB
        if (rad_interact_aer > 0) dalbsn=0.d0
        CALL RCOMPX
        NFSNBC(I,J)=SRNFLB(LM+LM_REQ+1)
c       NFSNBC(I,J)=SRNFLB(LFRC)
        ALBNBC(I,J)=SRNFLB(1)/(SRDFLB(1)+1.D-20)
c set for BC-albedo effect
        if (rad_interact_aer > 0) dALBsn=dALBsn1
#endif
C**** Optional calculation of CRF using a clear sky calc.
        if (cloud_rad_forc.gt.0) then
          FTAUC=0.   ! turn off cloud tau (tauic +tauwc)
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
          CALL RCOMPX          ! cloud_rad_forc>0 : clr sky
          SNFSCRF(I,J)=SRNFLB(LM+LM_REQ+1)   ! always TOA
          TNFSCRF(I,J)=TRNFLB(LM+LM_REQ+1)   ! always TOA
          LWDNCS(I,J)=STBO*(                 ! clr sky trhr(0)
     &      POCEAN*atmocn%GTEMPR(I,J)**4
     &     + POICE*atmice%GTEMPR(I,J)**4
     &     + PLICE*atmgla%GTEMPR(I,J)**4
     &     +PEARTH*atmlnd%GTEMPR(I,J)**4)
     &     -TRNFLB(1)
C         BEGIN AMIP
          AIJ(I,J,IJ_SWDCLS)=AIJ(I,J,IJ_SWDCLS)+SRDFLB(1)*COSZ2(I,J)
          AIJ(I,J,IJ_SWNCLS)=AIJ(I,J,IJ_SWNCLS)+SRNFLB(1)*COSZ2(I,J)
          AIJ(I,J,IJ_LWDCLS)=AIJ(I,J,IJ_LWDCLS)+TRDFLB(1)
          AIJ(I,J,IJ_SWNCLT)=AIJ(I,J,IJ_SWNCLT)+SRNFLB(LM+LM_REQ+1)
     *     *COSZ2(I,J)
          AIJ(I,J,IJ_LWNCLT)=AIJ(I,J,IJ_LWNCLT)+TRNFLB(LM+LM_REQ+1)
C       END AMIP
#ifdef CFMIP3_SUBDD
          ! SW upward flux at TOA, Csky
          !swutcs(i,j)=sruflb(lm)*csz2
          swutcs(i,j)=sruflb(lm)*cosz2(i,j)
          ! SW downward flux at SFC, Csky
          swdcls(i,j)=srdflb(1)*cosz2(i,j)
          ! SW upward flux at SFC, Csky
          swucls(i,j)=sruflb(1)*cosz2(i,j)
#endif
        end if
        FTAUC=1.     ! default: turn on cloud tau


C**** 2nd Optional calculation of CRF using a clear sky calc. without aerosols and Ox
        if (cloud_rad_forc.eq.2) then
          FTAUC=0.   ! turn off cloud tau (tauic +tauwc)
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
c Including turn off of aerosols and Ox during crf calc.+++++++++++++++++++
#ifdef TRACERS_SPECIAL_Shindell
       use_o3_ref=1 ; use_tracer_chem(1)=0  !turns off ozone
#endif
       FSTOPX(:) = 0 !turns off aerosol tracers
       FTTOPX(:) = 0
        CALL RCOMPX          ! cloud_rad_forc=2 : clr sky
       FSTOPX(:) = onoff_aer !turns on aerosol tracers, if requested
       FTTOPX(:) = onoff_aer !
#ifdef TRACERS_SPECIAL_Shindell
       use_o3_ref=0 ; use_tracer_chem(1)=onoff_chem*Lmax_rad_O3 ! turns on ozone tracers
#endif
          SNFSCRF2(I,J)=SRNFLB(LM+LM_REQ+1)   ! always TOA
          TNFSCRF2(I,J)=TRNFLB(LM+LM_REQ+1)   ! always TOA
        end if
        FTAUC=1.     ! default: turn on cloud tau

        if (cloud_rad_forc.gt.0) then
C**** all sky calc. without aerosol
       kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
       FSTOPX(:) = 0 !turns off aerosol tracers
       FTTOPX(:) = 0
        CALL RCOMPX          !  all sky
       FSTOPX(:) = onoff_aer !turns on aerosol tracers, if requested
       FTTOPX(:) = onoff_aer !

          SNFS_AS_noA(I,J)=SRNFLB(LM+LM_REQ+1)   ! always TOA
          TNFS_AS_noA(I,J)=TRNFLB(LM+LM_REQ+1)   ! always TOA

C**** clear sky calc. without aerosol
          FTAUC=0.   ! turn off cloud tau (tauic +tauwc)
       kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
       FSTOPX(:) = 0 !turns off aerosol tracers
       FTTOPX(:) = 0
        CALL RCOMPX          !  clr sky
       FSTOPX(:) = onoff_aer !turns on aerosol tracers, if requested
       FTTOPX(:) = onoff_aer !

          SNFS_CS_noA(I,J)=SRNFLB(LM+LM_REQ+1)   ! always TOA
          TNFS_CS_noA(I,J)=TRNFLB(LM+LM_REQ+1)   ! always TOA
        FTAUC=1.     ! default: turn on cloud tau
       end if

C**** Optional calculation of the impact of NINT aerosols
        if (aer_rad_forc.gt.0) then
C**** first, separate aerosols
          DO N=1,8
          tmpS(N)=FS8OPX(N)   ; tmpT(N)=FT8OPX(N)
          FS8OPX(N)=0.     ; FT8OPX(N)=0.
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
          CALL RCOMPX           ! aer_rad_forc>0 : no aerosol N
          SNFSAERRF(N,I,J)=SRNFLB(LM+LM_REQ+1) ! TOA
          TNFSAERRF(N,I,J)=TRNFLB(LM+LM_REQ+1) ! TOA
          SNFSAERRF(N+8,I,J)=SRNFLB(1) ! SURF
          TNFSAERRF(N+8,I,J)=TRNFLB(1) ! SURF
          FS8OPX(N)=tmpS(N)   ; FT8OPX(N)=tmpT(N)
          END DO
C**** second, net aerosols
          tmpS(:)=FS8OPX(:)   ; tmpT(:)=FT8OPX(:)
          FS8OPX(:)=0.     ; FT8OPX(:)=0.
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
          CALL RCOMPX             ! aer_rad_forc>0 : no aerosols
          SNFSAERRF(17,I,J)=SRNFLB(LM+LM_REQ+1) ! TOA
          TNFSAERRF(17,I,J)=TRNFLB(LM+LM_REQ+1) ! TOA
          SNFSAERRF(18,I,J)=SRNFLB(1) ! SURF
          TNFSAERRF(18,I,J)=TRNFLB(1) ! SURF
          FS8OPX(:)=tmpS(:)   ; FT8OPX(:)=tmpT(:)
        end if
      end if  ! moddrf=0

C**** End of initial computations for optional forcing diagnostics

C**** Localize fields that are modified by RCOMPX
      kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)

C*****************************************************
C     Main RADIATIVE computations, SOLAR and THERM(A)L
      CALL RCOMPX
C*****************************************************

#ifdef CACHED_SUBDD
      CO2out(1:LM,i,j)=CO2outCol(1:LM)
#endif
#ifdef GCC_UNCOUPLE_RAD_CONCEN
! Put back the actual amount of CO2 to fulgas
      fulgas(2) = GCCco2_fulgas_now   
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)

C**** Save optical depth diags
      nsub_ntrix = 0
      do n=1,nraero_aod
        SELECT CASE (trname(ntrix_aod(n)))
        CASE ('Clay','ClayIlli','ClayKaol','ClaySmec','ClayCalc'
     &         ,'ClayQuar','ClayFeld','ClayHema','ClayGyps'
     &         ,'ClayIlHe','ClayKaHe','ClaySmHe','ClayCaHe'
     &         ,'ClayQuHe','ClayFeHe','ClayGyHe')
          nsub_ntrix(ntrix_aod(n)) = nsub_ntrix(ntrix_aod(n)) + 1

! 3d aod
          if (diag_aod_3d>0 .and. diag_aod_3d<5) then ! valid values are 1-4
            if (ijlt_3Daaod(n).gt.0)
     *           taijls(i,j,1:lm,ijlt_3Daaod(n))
     *           =taijls(i,j,1:lm,ijlt_3Daaod(n))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))
            if (ijlt_3DaaodCS(n).gt.0)
     *           taijls(i,j,1:lm,ijlt_3DaaodCS(n))
     *           =taijls(i,j,1:lm,ijlt_3DaaodCS(n))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))*OPNSKY
            if (ijlt_3DaaodDRY(n).gt.0)
     *           taijls(i,j,1:lm,ijlt_3DaaodDRY(n))
     *           =taijls(i,j,1:lm,ijlt_3DaaodDRY(n))+
     *            (aesqex_dry(1:lm,6,n)-aesqsc_dry(1:lm,6,n))
            if (ijlt_3Dtau(n).gt.0)
     *           taijls(i,j,1:lm,ijlt_3Dtau(n))
     *         =taijls(i,j,1:lm,ijlt_3Dtau(n))+aesqex(1:lm,6,n)
            if (ijlt_3DtauCS(n).gt.0)
     *           taijls(i,j,1:lm,ijlt_3DtauCS(n))
     *         =taijls(i,j,1:lm,ijlt_3DtauCS(n))+aesqex(1:lm,6,n)*OPNSKY
            if (ijlt_3DtauDRY(n).gt.0)
     *           taijls(i,j,1:lm,ijlt_3DtauDRY(n))
     *         =taijls(i,j,1:lm,ijlt_3DtauDRY(n))+
     *          aesqex_dry(1:lm,6,n)
          else if (diag_aod_3d<0 .and. diag_aod_3d>-5) then ! if negative, save total
            if (ijlt_3Daaod(1).gt.0)
     *           taijls(i,j,1:lm,ijlt_3Daaod(1))
     *           =taijls(i,j,1:lm,ijlt_3Daaod(1))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))
            if (ijlt_3DaaodCS(1).gt.0)
     *           taijls(i,j,1:lm,ijlt_3DaaodCS(1))
     *           =taijls(i,j,1:lm,ijlt_3DaaodCS(1))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))*OPNSKY
            if (ijlt_3DaaodDRY(1).gt.0)
     *           taijls(i,j,1:lm,ijlt_3DaaodDRY(1))
     *           =taijls(i,j,1:lm,ijlt_3DaaodDRY(1))+
     *            (aesqex_dry(1:lm,6,n)-aesqsc_dry(1:lm,6,n))
            if (ijlt_3Dtau(1).gt.0)
     *           taijls(i,j,1:lm,ijlt_3Dtau(1))
     *         =taijls(i,j,1:lm,ijlt_3Dtau(1))+aesqex(1:lm,6,n)
            if (ijlt_3DtauCS(1).gt.0)
     *           taijls(i,j,1:lm,ijlt_3DtauCS(1))
     *         =taijls(i,j,1:lm,ijlt_3DtauCS(1))+aesqex(1:lm,6,n)*OPNSKY
            if (ijlt_3DtauDRY(1).gt.0)
     *           taijls(i,j,1:lm,ijlt_3DtauDRY(1))
     *         =taijls(i,j,1:lm,ijlt_3DtauDRY(1))+
     *          aesqex_dry(1:lm,6,n)
          endif ! 0<diag_aod_3d<5 or 0>diag_aod_3d>-5

! 2d aod, per band or just band6, depending on diag_rad
          IF (diag_rad /= 1) THEN
            IF ( ijts_tausub(1,ntrix_aod(n),nsub_ntrix(ntrix_aod(n)))>0
     &           )taijs(i,j,ijts_tausub(1,ntrix_aod(n)
     &           ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j,ijts_tausub(1
     &           ,ntrix_aod(n),nsub_ntrix(ntrix_aod(n)))) +
     &            SUM(aesqex(1:Lm,6,n))
            IF ( ijts_tausub(2,ntrix_aod(n),nsub_ntrix(ntrix_aod(n)))>0
     &           )taijs(i,j,ijts_tausub(2,ntrix_aod(n)
     &           ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j,ijts_tausub(2
     &           ,ntrix_aod(n),nsub_ntrix(ntrix_aod(n)))) +
     &            SUM(aesqex(1:Lm,6,n)) * OPNSKY
            IF ( ijts_tausub(3,ntrix_aod(n),nsub_ntrix(ntrix_aod(n)))>0
     &           )taijs(i,j,ijts_tausub(3,ntrix_aod(n)
     &           ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j,ijts_tausub(3
     &           ,ntrix_aod(n),nsub_ntrix(ntrix_aod(n)))) +
     &            SUM(aesqex_dry(1:Lm,6,n))
          ELSE
            DO kr=1,6
              IF ( ijts_sqexsub(1,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n))) > 0 ) taijs(i,j
     &             ,ijts_sqexsub(1,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j
     &             ,ijts_sqexsub(1,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n))))+SUM(aesqex(1:Lm ,kr,n))
              IF ( ijts_sqexsub(2,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n))) >0 ) taijs(i,j
     &             ,ijts_sqexsub(2,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j
     &             ,ijts_sqexsub(2,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) +SUM(aesqex(1:Lm,kr,n)) *
     &             OPNSKY
              IF ( ijts_sqexsub(3,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n))) >0 ) taijs(i,j
     &             ,ijts_sqexsub(3,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j
     &             ,ijts_sqexsub(3,kr,ntrix_aod(n)
     &           ,nsub_ntrix(ntrix_aod(n)))) +SUM(aesqex_dry(1:Lm,kr,n))
              IF ( ijts_sqscsub(1,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n))) >0 ) taijs(i,j
     &             ,ijts_sqscsub(1,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j
     &             ,ijts_sqscsub(1,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) +SUM(aesqsc(1:Lm,kr,n))
              IF ( ijts_sqscsub(2,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n))) >0 ) taijs(i,j
     &             ,ijts_sqscsub(2,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j
     &             ,ijts_sqscsub(2,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) +SUM(aesqsc(1:Lm,kr,n)) *
     &             OPNSKY
              IF ( ijts_sqscsub(3,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n))) >0 ) taijs(i,j
     &             ,ijts_sqscsub(3,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j
     &             ,ijts_sqscsub(3,kr,ntrix_aod(n)
     &           ,nsub_ntrix(ntrix_aod(n)))) +SUM(aesqsc_dry(1:Lm,kr,n))
              IF ( ijts_sqcbsub(1,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n))) >0 ) taijs(i,j
     &             ,ijts_sqcbsub(1,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j
     &             ,ijts_sqcbsub(1,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) +SUM(aesqcb(1:Lm,kr,n)) /
     &             (SUM(aesqsc(1:Lm,kr,n))+1.D-10)
              IF (ijts_sqcbsub(2,kr,ntrix_aod(n),
     &             nsub_ntrix(ntrix_aod(n)))
     &             >0) taijs(i,j,ijts_sqcbsub(2,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j
     &             ,ijts_sqcbsub(2,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) +SUM(aesqcb(1:Lm,kr,n)) /
     &             (SUM(aesqsc(1:Lm,kr,n))+1.D-10) * OPNSKY
              IF (ijts_sqcbsub(3,kr,ntrix_aod(n),
     &             nsub_ntrix(ntrix_aod(n)))
     &             >0) taijs(i,j,ijts_sqcbsub(3,kr,ntrix_aod(n)
     &             ,nsub_ntrix(ntrix_aod(n)))) = taijs(i,j
     &             ,ijts_sqcbsub(3,kr,ntrix_aod(n)
     &         ,nsub_ntrix(ntrix_aod(n)))) +SUM(aesqcb_dry(1:Lm,kr,n)) /
     &             (SUM(aesqsc_dry(1:Lm,kr,n))+1.D-10)
            END DO
          END IF
        CASE DEFAULT

! 3d aod
          if (diag_aod_3d>0 .and. diag_aod_3d<5) then ! valid values are 1-4
            if (ijlt_3Daaod(n).gt.0)
     &           taijls(i,j,1:lm,ijlt_3Daaod(n))
     &           =taijls(i,j,1:lm,ijlt_3Daaod(n))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))
            if (ijlt_3DaaodCS(n).gt.0)
     &           taijls(i,j,1:lm,ijlt_3DaaodCS(n))
     &           =taijls(i,j,1:lm,ijlt_3DaaodCS(n))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))*OPNSKY
            if (ijlt_3DaaodDRY(n).gt.0)
     &           taijls(i,j,1:lm,ijlt_3DaaodDRY(n))
     &           =taijls(i,j,1:lm,ijlt_3DaaodDRY(n))+
     *            (aesqex_dry(1:lm,6,n)-aesqsc(1:lm,6,n))
            if (ijlt_3Dtau(n).gt.0)
     &           taijls(i,j,1:lm,ijlt_3Dtau(n))
     &         =taijls(i,j,1:lm,ijlt_3Dtau(n))+aesqex(1:lm,6,n)
            if (ijlt_3DtauCS(n).gt.0)
     &           taijls(i,j,1:lm,ijlt_3DtauCS(n))
     &         =taijls(i,j,1:lm,ijlt_3DtauCS(n))+aesqex(1:lm,6,n)*OPNSKY
            if (ijlt_3DtauDRY(n).gt.0)
     &           taijls(i,j,1:lm,ijlt_3DtauDRY(n))
     &         =taijls(i,j,1:lm,ijlt_3DtauDRY(n))+aesqex_dry(1:lm,6,n)
          else if (diag_aod_3d<0 .and. diag_aod_3d>-5) then ! if negative, save total
            if (ijlt_3Daaod(1).gt.0)
     &           taijls(i,j,1:lm,ijlt_3Daaod(1))
     &           =taijls(i,j,1:lm,ijlt_3Daaod(1))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))
            if (ijlt_3DaaodCS(1).gt.0)
     &           taijls(i,j,1:lm,ijlt_3DaaodCS(1))
     &           =taijls(i,j,1:lm,ijlt_3DaaodCS(1))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))*OPNSKY
            if (ijlt_3DaaodDRY(1).gt.0)
     &           taijls(i,j,1:lm,ijlt_3DaaodDRY(1))
     &           =taijls(i,j,1:lm,ijlt_3DaaodDRY(1))+
     *            (aesqex_dry(1:lm,6,n)-aesqsc_dry(1:lm,6,n))
            if (ijlt_3Dtau(1).gt.0)
     &           taijls(i,j,1:lm,ijlt_3Dtau(1))
     &         =taijls(i,j,1:lm,ijlt_3Dtau(1))+aesqex(1:lm,6,n)
            if (ijlt_3DtauCS(1).gt.0)
     &           taijls(i,j,1:lm,ijlt_3DtauCS(1))
     &         =taijls(i,j,1:lm,ijlt_3DtauCS(1))+aesqex(1:lm,6,n)*OPNSKY
            if (ijlt_3DtauDRY(1).gt.0)
     &           taijls(i,j,1:lm,ijlt_3DtauDRY(1))
     &         =taijls(i,j,1:lm,ijlt_3DtauDRY(1))+aesqex_dry(1:lm,6,n)
          endif ! 0<diag_aod_3d<5 or 0>diag_aod_3d>-5

! 2d aod, per band or just band6, depending on diag_rad
          IF (diag_rad /= 1) THEN
            if (ijts_tau(1,ntrix_aod(n)).gt.0)
     &           taijs(i,j,ijts_tau(1,ntrix_aod(n)))
     &         =taijs(i,j,ijts_tau(1,ntrix_aod(n)))+
     &          SUM(aesqex(1:lm,6,n))
            if (ijts_tau(2,ntrix_aod(n)).gt.0)
     &           taijs(i,j,ijts_tau(2,ntrix_aod(n)))
     &           =taijs(i,j,ijts_tau(2,ntrix_aod(n)))
     &           +SUM(aesqex(1:lm,6,n))*OPNSKY
            if (ijts_tau(3,ntrix_aod(n)).gt.0)
     &           taijs(i,j,ijts_tau(3,ntrix_aod(n)))
     &           =taijs(i,j,ijts_tau(3,ntrix_aod(n)))
     &           +SUM(aesqex_dry(1:lm,6,n))
          ELSE
            DO kr=1,6
c               print*,'SUSA  diag',SUM(aesqex(1:Lm,kr,n))
              IF (ijts_sqex(1,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqex(1,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqex(1,kr,ntrix_aod(n)))
     &             +SUM(aesqex(1:Lm,kr,n))
              IF (ijts_sqex(2,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqex(2,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqex(2,kr,ntrix_aod(n)))
     &             +SUM(aesqex(1:Lm,kr,n))*OPNSKY
              IF (ijts_sqex(3,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqex(3,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqex(3,kr,ntrix_aod(n)))
     &             +SUM(aesqex_dry(1:Lm,kr,n))
              IF (ijts_sqsc(1,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqsc(1,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqsc(1,kr,ntrix_aod(n)))
     &             +SUM(aesqsc(1:Lm,kr,n))
              IF (ijts_sqsc(2,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqsc(2,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqsc(2,kr,ntrix_aod(n)))
     &             +SUM(aesqsc(1:Lm,kr,n))*OPNSKY
              IF (ijts_sqsc(3,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqsc(3,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqsc(3,kr,ntrix_aod(n)))
     &             +SUM(aesqsc_dry(1:Lm,kr,n))
#ifndef TRACERS_TOMAS
              IF (ijts_sqcb(1,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqcb(1,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqcb(1,kr,ntrix_aod(n)))
     &             +SUM(aesqcb(1:Lm,kr,n))
     &             /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)
              IF (ijts_sqcb(2,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqcb(2,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqcb(2,kr,ntrix_aod(n)))
     &             +SUM(aesqcb(1:Lm,kr,n))
     &             /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)*OPNSKY
              IF (ijts_sqcb(3,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqcb(3,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqcb(3,kr,ntrix_aod(n)))
     &             +SUM(aesqcb_dry(1:Lm,kr,n))
     &             /(SUM(aesqsc_dry(1:Lm,kr,n))+1.D-10)
#else
              qcb_col(kr,n)=0.d0
              qcb_col_dry(kr,n)=0.d0
              do l=1,lm
                qcb_col(kr,n)=qcb_col(kr,n)+aesqcb(l,kr,n)*
     *               aesqsc(l,kr,n)
                qcb_col_dry(kr,n)=qcb_col(kr,n)+aesqcb_dry(l,kr,n)*
     *               aesqsc_dry(l,kr,n)
              enddo

              IF (ijts_sqcb(1,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqcb(1,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqcb(1,kr,ntrix_aod(n)))
     &             +qcb_col(kr,n)
     &             /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)
              IF (ijts_sqcb(2,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqcb(2,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqcb(2,kr,ntrix_aod(n)))
     &             +qcb_col(kr,n)
     &             /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)*OPNSKY
              IF (ijts_sqcb(3,kr,ntrix_aod(n)) > 0)
     &             taijs(i,j,ijts_sqcb(3,kr,ntrix_aod(n)))
     &             =taijs(i,j,ijts_sqcb(3,kr,ntrix_aod(n)))
     &             +qcb_col_dry(kr,n)
     &             /(SUM(aesqsc_dry(1:Lm,kr,n))+1.D-10)
#endif
            END DO ! kr
          END IF ! diag_rad
        END SELECT ! clay or not
      end do ! nraero_aod

#endif  /* Koch||DUST||MINERALS||AMP||TOMAS||SEASALT */

      if (TAero_aod_diag > 0) then
        do n=1,8 ! 8 radiatively active aerosol tracers
          do kr=1,6 ! 6 bands in the shortwave
            if (TAero_aod_diag == 2 .and. kr /= 6) cycle ! only save band6
            aij(i,j,ij_nintaerext(kr,n))=aij(i,j,ij_nintaerext(kr,n))+
     &        sum(nintaerext(1:lm,kr,n))
            aij(i,j,ij_nintaersca(kr,n))=aij(i,j,ij_nintaersca(kr,n))+
     &        sum(nintaersca(1:lm,kr,n))
            aij(i,j,ij_nintaerasy(kr,n))=aij(i,j,ij_nintaerasy(kr,n))+
     &        sum(nintaerasy(1:lm,kr,n)*nintaersca(1:lm,kr,n))/
     &        (sum(nintaersca(1:lm,kr,n))+1.d-10)
          enddo ! kr
        enddo ! n
      endif


#ifdef TRACERS_ON
      if (nraero_aod>0) then
        tau_as(i,j,1:LM,1:nraero_aod)=aesqex(1:LM,6,1:nraero_aod)
        tau_cs(i,j,1:LM,1:nraero_aod)=aesqex(1:LM,6,1:nraero_aod)*OPNSKY
        if (save_dry_aod>0) then
          tau_dry(i,j,1:LM,1:nraero_aod)=aesqex_dry(1:LM,6,1:nraero_aod)
        endif
#ifdef CACHED_SUBDD
        abstau_as(i,j,1:LM,1:nraero_aod)=
     &    (aesqex(1:LM,6,1:nraero_aod)-aesqsc(1:LM,6,1:nraero_aod))
        abstau_cs(i,j,1:LM,1:nraero_aod)=
     &    (aesqex(1:LM,6,1:nraero_aod)-aesqsc(1:LM,6,1:nraero_aod))*
     &    OPNSKY
        if (save_dry_aod>0) then
          abstau_dry(i,j,1:LM,1:nraero_aod)=
     & (aesqex_dry(1:LM,6,1:nraero_aod)-aesqsc_dry(1:LM,6,1:nraero_aod))
        endif
#endif  /* CACHED_SUBDD */
      endif
#endif /* TRACERS_ON */

      IF (I.EQ.IWRITE .and. J.EQ.JWRITE) CALL WRITER(6,ITWRITE)
      CSZ2=COSZ2(I,J)
      do L=1,LM
#ifdef GCC_COUPLE_RAD
        GCCco2rad_to_chem(L,i,j)=GCCco2_out(L)
#endif
        rad_to_chem(:,L,i,j)=chem_out(L,:)
        rad_to_chem(4,L,i,j)=chem_out(L,4)/CH4X_RADoverCHEM
        do k=1,4
          kliq(L,k,i,j)=kdeliq(L,k) ! save updated flags
        end do
      end do
      if (kradia.gt.0) then  ! rad. forc. model; acc diagn
        do L=1,LM+LM_REQ+1
          AFLX_ST(L,I,J,1)=AFLX_ST(L,I,J,1)+SRUFLB(L)*CSZ2
          AFLX_ST(L,I,J,2)=AFLX_ST(L,I,J,2)+SRDFLB(L)*CSZ2
          AFLX_ST(L,I,J,3)=AFLX_ST(L,I,J,3)+TRUFLB(L)
          AFLX_ST(L,I,J,4)=AFLX_ST(L,I,J,4)+TRDFLB(L)
        end do
        if(kradia.eq.1) then
          tauex6=0. ; tauex5=0. ; tausct=0. ; taugcb=0.
          do L=1,LM
            AFLX_ST(L,I,J,5)=AFLX_ST(L,I,J,5)+1.d2*RHL(L)
            tauex6=tauex6+SRAEXT(L,6)+SRDEXT(L,6)+SRVEXT(L,6)
            tauex5=tauex5+SRAEXT(L,5)+SRDEXT(L,5)+SRVEXT(L,5)
            tausct=tausct+SRASCT(L,6)+SRDSCT(L,6)+SRVSCT(L,6)
            taugcb=taugcb+SRASCT(L,6)*SRAGCB(L,6)+
     +        SRDSCT(L,6)*SRDGCB(L,6)+SRVSCT(L,6)*SRVGCB(L,6)
          end do
          AFLX_ST(LM+1,I,J,5)=AFLX_ST(LM+1,I,J,5)+tauex5
          AFLX_ST(LM+2,I,J,5)=AFLX_ST(LM+2,I,J,5)+tauex6
          AFLX_ST(LM+3,I,J,5)=AFLX_ST(LM+3,I,J,5)+tausct
          AFLX_ST(LM+4,I,J,5)=AFLX_ST(LM+4,I,J,5)+taugcb
          cycle
        end if
        do l=LS1_loc,lm
          tchg(l,i,j) = tchg(l,i,j) + ( srfhrl(l)*csz2-srhra(l,i,j) +
     +      (-trfcrl(l)-trhra(l,i,j)) )*nrad*DTsrc*bysha*byMA(l,i,j)
        end do
        do l=lm+1,lm+lm_req
          tchg(l,i,j) = tchg(l,i,j) + ( srfhrl(l)*csz2-srhra(l,i,j) +
     +      (-trfcrl(l)-trhra(l,i,j)) )*nrad*DTsrc*bysha*byaml00(l)
        end do
        cycle
      else if (kradia.lt.0) then ! save i/o data for frc.runs
        fmp_com(i,j) = fmp                  ! input data
        wsoil(i,j) = wearth
        do L=1,LM
          QR(L,I,J) = shl(L)
          CLDinfo(L,1,I,J) = tauwc(L)
          CLDinfo(L,2,I,J) = tauic(L)
          CLDinfo(L,3,I,J) = sizeic(L)  ! sizeic=sizewc currently
        end do
        SRHRA(0,I,J)=SRNFLB(1)*CSZ2      ! output data (for adj frc)
        TRHRA(0,I,J)=-TRNFLB(1)
        do L=1,LM+LM_REQ
          SRHRA(L,I,J)=SRFHRL(L)*CSZ2
          TRHRA(L,I,J)=-TRFCRL(L)
        end do
      end if
C****
C**** Save relevant output in model arrays
C****
C**** (some generalisation and coherence needed in the rad surf type calc)
      FSF(1,I,J)=FSRNFG(1)   !  ocean
      FSF(2,I,J)=FSRNFG(3)   !  ocean ice
      FSF(3,I,J)=FSRNFG(4)   !  land ice
      FSF(4,I,J)=FSRNFG(2)   !  soil
      SRHR(0,I,J)=SRNFLB(1)
      TRHR(0,I,J)=STBO*(
     &      POCEAN*atmocn%GTEMPR(I,J)**4
     &     + POICE*atmice%GTEMPR(I,J)**4
     &     + PLICE*atmgla%GTEMPR(I,J)**4
     &     +PEARTH*atmlnd%GTEMPR(I,J)**4)
     &     -TRNFLB(1)
      TRSURF(1,I,J) = STBO*atmocn%GTEMPR(I,J)**4  !  ocean
      TRSURF(2,I,J) = STBO*atmice%GTEMPR(I,J)**4  !  ocean ice
      TRSURF(3,I,J) = STBO*atmgla%GTEMPR(I,J)**4  !  land ice
      TRSURF(4,I,J) = STBO*atmlnd%GTEMPR(I,J)**4  !  soil
      DO L=1,LM
        SRHR(L,I,J)=SRFHRL(L)
        TRHR(L,I,J)=-TRFCRL(L)
      END DO
      DO LR=1,LM_REQ
        SRHRS(LR,I,J)= SRFHRL(LM+LR)
        TRHRS(LR,I,J)=-TRFCRL(LM+LR)
      END DO
#ifdef SCM
c**** possibly turn off radiative heating in atmosphere
c**** and use specified profile for thermal heating rate
c**** converting units from K/s to W/m2
      if( SCMopt%Qrad )then
        SRHR(1:LM,I,J)=0.
        TRHR(1:LM,I,J)=0.
        SRHRS(1:LM_REQ,I,J)=0.
        TRHRS(1:LM_REQ,I,J)=0.
        TRHR(1:LM,I,J)=SCMin%Qrad(1:LM)*SHA*MA(1:LM,I,J)
      endif
c**** possibly turn off radiative heating in atmosphere
c**** and use Beers Law for thermal heating rate as
c**** difference of net flux over layer
      if( SCMopt%BeersLaw )then
        SRHR(1:LM,I,J)=0.
        TRHR(1:LM,I,J)=0.
        SRHRS(1:LM_REQ,I,J)=0.
        TRHRS(1:LM_REQ,I,J)=0.
        ! cumulative cloud water paths * extinction coefficient
        q_above(LM+1) = 0.
        do L=LM,1,-1
          q_above(L) = q_above(L+1) +
     &      SCMin%BeersLaw_kappa*MA(L,i,j)*QCL(i,j,L)
        enddo
        q_below(1) = 0.
        do L=1,LM
          q_below(L+1) = q_below(L) +
     &      SCMin%BeersLaw_kappa*MA(L,i,j)*QCL(i,j,L)
        enddo
        ! net upward radiative flux at layer edges
        Frad(:) = SCMin%BeersLaw_f0*exp(-q_above(:)) +
     &            SCMin%BeersLaw_f1*exp(-q_below(:))
        ! radiative flux difference over each layer
        TRHR(1:LM,I,J) = Frad(1:LM) - Frad(2:LM+1)
      endif
c**** save radiative flux profiles for sub-daily output
#ifdef CACHED_SUBDD
      TRDFLB_prof(I,J,1:LM) = TRDFLB(1:LM)
      TRUFLB_prof(I,J,1:LM) = TRUFLB(1:LM)
      SRDFLB_prof(I,J,1:LM) = SRDFLB(1:LM)
      SRUFLB_prof(I,J,1:LM) = SRUFLB(1:LM)
#endif
#endif
C**** Save fluxes at four levels surface, P0, P1, LTROPO
      SNFS(1,I,J)=SRNFLB(1)     ! Surface
      TNFS(1,I,J)=TRNFLB(1)
      SNFS(2,I,J)=SRNFLB(LM+1)  ! P1
      TNFS(2,I,J)=TRNFLB(LM+1)
      SNFS(3,I,J)=SRNFLB(LM+LM_REQ+1) ! P0 = TOA
      TNFS(3,I,J)=TRNFLB(LM+LM_REQ+1)
      SNFS(4,I,J)=SRNFLB(LTROPO(I,J)) ! LTROPO
      TNFS(4,I,J)=TRNFLB(LTROPO(I,J))

C****
      TRINCG(I,J)=TRDFLB(1)
      BTMPW(I,J)=BTEMPW-TF
      ALB(I,J,1)=SRNFLB(1)/(SRDFLB(1)+1.D-20)
      ALB(I,J,2)=PLAVIS
      ALB(I,J,3)=PLANIR
      ALB(I,J,4)=ALBVIS
      ALB(I,J,5)=ALBNIR
      ALB(I,J,6)=SRRVIS
      ALB(I,J,7)=SRRNIR
      ALB(I,J,8)=SRAVIS
      ALB(I,J,9)=SRANIR

#ifdef TRACERS_DUST
      IF (adiurn_dust == 1) THEN
        srnflb_save(i,j,1:lm)=srnflb(1:lm)
        trnflb_save(i,j,1:lm)=trnflb(1:lm)
      END IF
#endif
#ifdef mjo_subdd
      swu_avg(I,J)=swu_avg(I,J)+SRUFLB(1)*CSZ2
#endif

      SWUS(I,J)=SRUFLB(1)*CSZ2
#ifdef CFMIP3_SUBDD
      ! SW upward flux at TOA
      swut(i,j)=sruflb(lm)*csz2
      ! SW downward flux at TOA
      swdt(i,j)=srdflb(lm)*csz2
#endif
      SRDN(I,J) = SRDFLB(1)     ! save total solar flux at surface
C**** SALB(I,J)=ALB(I,J,1)      ! save surface albedo (pointer)
      FSRDIR(I,J)=SRXVIS        ! direct visible solar at surface **coefficient
      SRVISSURF(I,J)=SRDVIS     ! total visible solar at surface
      DIRVIS(I,J)=SRXVIS*SRDVIS  ! direct visible solar at surface
      FSRDIF(I,J)=SRDVIS*(1-SRXVIS) ! diffuse visible solar at surface

      DIRNIR(I,J)=SRXNIR*SRDNIR     ! direct beam nir solar at surface
      DIFNIR(I,J)=SRDNIR*(1-SRXNIR) ! diffuse     nir solar at surface

cdiag write(*,'(a,2i5,6e12.4)')'RAD_DRV: ',
cdiag.    I,J,FSRDIR(I,J),SRVISSURF(I,J),FSRDIF(I,J),
cdiag.        DIRNIR(I,J),SRDNIR,DIFNIR(I,J)
C**** Save clear sky/tropopause diagnostics here
      AIJ(I,J,IJ_CLR_SRINCG)=AIJ(I,J,IJ_CLR_SRINCG)+OPNSKY*
     *     SRDFLB(1)*CSZ2
      AIJ(I,J,IJ_CLR_SRNFG)=AIJ(I,J,IJ_CLR_SRNFG)+OPNSKY*
     *     SRNFLB(1)*CSZ2
      AIJ(I,J,IJ_CLR_TRDNG)=AIJ(I,J,IJ_CLR_TRDNG)+OPNSKY*TRHR(0,I,J)
      AIJ(I,J,IJ_CLR_SRUPTOA)=AIJ(I,J,IJ_CLR_SRUPTOA)+OPNSKY*
     *     SRUFLB(LM+LM_REQ+1)*CSZ2
      AIJ(I,J,IJ_CLR_TRUPTOA)=AIJ(I,J,IJ_CLR_TRUPTOA)+OPNSKY*
     *     TRUFLB(LM+LM_REQ+1)
      AIJ(I,J,IJ_CLR_SRNTP)=AIJ(I,J,IJ_CLR_SRNTP)+OPNSKY*
     *     SRNFLB(LTROPO(I,J))*CSZ2
      AIJ(I,J,IJ_CLR_TRNTP)=AIJ(I,J,IJ_CLR_TRNTP)+OPNSKY*
     *     TRNFLB(LTROPO(I,J))
      AIJ(I,J,IJ_SRNTP)=AIJ(I,J,IJ_SRNTP)+SRNFLB(LTROPO(I,J))*CSZ2
      AIJ(I,J,IJ_TRNTP)=AIJ(I,J,IJ_TRNTP)+TRNFLB(LTROPO(I,J))
      AIJ(I,J,IJ_SISWD)=AIJ(I,J,IJ_SISWD)+POICE*SRDFLB(1)*CSZ2
      AIJ(I,J,IJ_SISWU)=AIJ(I,J,IJ_SISWU)+
     *     POICE*(SRDFLB(1)-FSRNFG(3))*CSZ2

      DO IT=1,NTYPE
         call inc_aj(i,j,it,J_CLRTOA,OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *        *CSZ2-TRNFLB(LM+LM_REQ+1))*FTYPE(IT,I,J))
         call inc_aj(i,j,it,J_CLRTRP,OPNSKY*(SRNFLB(LTROPO(I,J))
     *        *CSZ2-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J))
         call inc_aj(i,j,it,J_TOTTRP,(SRNFLB(LTROPO(I,J))
     *        *CSZ2-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J))
      END DO
      call inc_areg(i,j,jr,J_CLRTOA,OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *     *CSZ2-TRNFLB(LM+LM_REQ+1)))
      call inc_areg(i,j,jr,J_CLRTRP,OPNSKY*(SRNFLB(LTROPO(I,J))
     *     *CSZ2-TRNFLB(LTROPO(I,J))))
      call inc_areg(i,j,jr,J_TOTTRP,(SRNFLB(LTROPO(I,J))
     *     *CSZ2-TRNFLB(LTROPO(I,J))))
C**** Save cloud top diagnostics here
      if (CLDCV.le.0.) go to 590
      AIJ(I,J,IJ_CLDTPPR)=AIJ(I,J,IJ_CLDTPPR)+plb(ltopcl+1)
      AIJ(I,J,IJ_CLDTPT)=AIJ(I,J,IJ_CLDTPT)+(tlb(ltopcl+1) - tf)
      CTT(i,j) = (tlb(ltopcl+1) - tf)
      CTP(i,j) = plb(ltopcl+1)
C**** Save cloud tau=1 related diagnostics here (opt.depth=1 level)
      tauup=0.
      DO L=LM,1,-1
         taucl=tauwc(l)+tauic(l)
         taudn=tauup+taucl
         if (taudn.gt.1.) then
            aij(i,j,ij_cldcv1)=aij(i,j,ij_cldcv1)+1.
            wtlin=(1.-tauup)/taucl
            aij(i,j,ij_cldt1t)=aij(i,j,ij_cldt1t)+( tlb(l+1)-tf +
     +           (tlb(l)-tlb(l+1))*wtlin )
            aij(i,j,ij_cldt1p)=aij(i,j,ij_cldt1p)+( plb(l+1)+
     +           (plb(l)-plb(l+1))*wtlin )
            go to 590
         end if
         tauup=taudn
      end do
 590  continue

      END DO
C****
C**** END OF MAIN LOOP FOR I INDEX
C****

      END DO
C****
C**** END OF MAIN LOOP FOR J INDEX
C****

#ifdef mjo_subdd
       swu_cnt=swu_cnt+1.
#endif

      if(kradia.gt.0) then
         call stopTimer('RADIA()')
         return
      end if
C**** Stop if temperatures were out of range
C**** Now only warning messages are printed for T,Q errors
c     IF(ICKERR.GT.0)
c     &     call stop_model('In Radia: Temperature out of range',11)
c     IF(JCKERR.GT.0)  call stop_model('In Radia: RQT out of range',11)
c     IF(KCKERR.GT.0)  call stop_model('In Radia: Q<0',255)
C**** save all input data to disk if kradia<0
      if (kradia.lt.0) write(iu_rad) itime
     &     ,T,RQT,atmsrf%TsAvg   ! LM+LM_REQ+1+
     &     ,QR,P,CLDinfo,rsi,zsi ! LM+1+3*LM+1+1+
!     &     ,(((GTEMPR(k,i,j),k=1,4),i=1,im),j=1,jm) ! (4+)
     &     ,wsoil,atmsrf%wsavg,snowi,atmgla%snow,atmlnd%snowe ! 1+1+1+1+1+
     &     ,snoage,fmp_com,flag_dsws,ltropo ! 3+1+.5+.5+
     &     ,atmlnd%fr_snow_rad,dlake,flake ! 2+1+1
C**** output data: really needed only if kradia=2
     &     ,srhra,trhra         ! 2(LM+LM_REQ+1)
     &     ,itime
C****
C**** ACCUMULATE THE RADIATION DIAGNOSTICS
C****
      bydpreq(:) = 1d0/(req_fac_d(:)*pmtop)
      DO 780 J=J_0,J_1
         DO 770 I=I_0,IMAXJ(J)
            do l=1,lm
              call inc_ajl(i,j,l,jl_srhr,SRHR(L,I,J)*COSZ2(I,J))
              call inc_ajl(i,j,l,jl_trcr,TRHR(L,I,J))
            enddo
            CSZ2=COSZ2(I,J)
            JR=JREG(I,J)
            DO LR=1,LM_REQ
              call inc_asjl(i,j,lr,3,bydpreq(lr)*SRHRS(LR,I,J)*CSZ2)
              call inc_asjl(i,j,lr,4,bydpreq(lr)*TRHRS(LR,I,J))
            END DO
            DO KR=1,NDIUPT
            IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
              TMP(idd_aot) =SUM(aesqex(1:Lm,6,1:nraero_aod))!*OPNSKY
              TMP(idd_aot2) =SUM(aesqsc(1:Lm,6,1:nraero_aod))!*OPNSKY
#endif
              TMP(IDD_PALB)=(1.-SNFS(3,I,J)/S0)
              TMP(IDD_GALB)=(1.-ALB(I,J,1))
              TMP(IDD_ABSA)=(SNFS(3,I,J)-SRHR(0,I,J))*CSZ2
              DO INCH=1,NRAD
                IHM=1+(JTIME+INCH-1)*HR_IN_DAY/NDAY
                IH=IHM
                IF(IH.GT.HR_IN_DAY) IH = IH - HR_IN_DAY
                ADIURN(IDXB(:),KR,IH)=ADIURN(IDXB(:),KR,IH)+TMP(IDXB(:))
#ifdef USE_HDIURN
                IHM = IHM+(DATE-1)*HR_IN_DAY
                IF(IHM.LE.HR_IN_MONTH) THEN
                  HDIURN(IDXB(:),KR,IHM)=HDIURN(IDXB(:),KR,IHM)+
     &                 TMP(IDXB(:))
                ENDIF
#endif
              ENDDO
            END IF
            END DO

      DO IT=1,NTYPE
         call inc_aj(I,J,IT,J_SRINCP0,(S0*CSZ2)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_SRNFP0 ,(SNFS(3,I,J)*CSZ2)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_SRINCG ,(SRHR(0,I,J)*CSZ2/(ALB(I,J,1)+1.D
     *        -20))*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_BRTEMP ,BTMPW(I,J) *FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_TRINCG ,TRINCG(I,J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_HSURF  ,-(TNFS(3,I,J)-TNFS(1,I,J))
     *        *FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_TRNFP0 ,-TNFS(3,I,J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_TRNFP1 ,-TNFS(2,I,J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_SRNFP1 ,SNFS(2,I,J)*CSZ2*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_HATM   ,-(TNFS(2,I,J)-TNFS(1,I,J))
     *        *FTYPE(IT,I,J))
#ifdef HEALY_LM_DIAGS
         call inc_aj(I,J,IT,j_vtau ,10.*-20*VTAULAT(J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,j_ghg ,10.*ghg_totforc*FTYPE(IT,I,J))
#endif

      END DO
C**** Note: confusing because the types for radiation are a subset
         call inc_aj(I,J,ITOCEAN,J_SRNFG,(FSF(1,I,J)*CSZ2)*FOCEAN(I,J)
     *        *(1.-RSI(I,J)))
         call inc_aj(I,J,ITLAKE ,J_SRNFG,(FSF(1,I,J)*CSZ2)* FLAKE(I,J)
     *        *(1.-RSI(I,J)))
         call inc_aj(I,J,ITEARTH,J_SRNFG,(FSF(4,I,J)*CSZ2)*FEARTH(I,J))
         call inc_aj(I,J,ITLANDI,J_SRNFG,(FSF(3,I,J)*CSZ2)* FLICE(I,J))
         call inc_aj(I,J,ITOICE ,J_SRNFG,(FSF(2,I,J)*CSZ2)*FOCEAN(I,J)
     *        *RSI(I,J))
         call inc_aj(I,J,ITLKICE,J_SRNFG,(FSF(2,I,J)*CSZ2)* FLAKE(I,J)
     *        *RSI(I,J))
C****
         call inc_areg(I,J,JR,J_SRINCP0,(S0*CSZ2))
         call inc_areg(I,J,JR,J_SRNFP0 ,(SNFS(3,I,J)*CSZ2))
         call inc_areg(I,J,JR,J_SRNFP1 ,(SNFS(2,I,J)*CSZ2))
         call inc_areg(I,J,JR,J_SRINCG ,(SRHR(0,I,J)*CSZ2/
     *        (ALB(I,J,1)+1.D-20)))
         call inc_areg(I,J,JR,J_HATM,-(TNFS(2,I,J)-TNFS(1,I,J)))
         call inc_areg(I,J,JR,J_SRNFG,(SRHR(0,I,J)*CSZ2))
         call inc_areg(I,J,JR,J_HSURF,-(TNFS(3,I,J)-TNFS(1,I,J)))
         call inc_areg(I,J,JR,J_BRTEMP, BTMPW(I,J))
         call inc_areg(I,J,JR,J_TRINCG,TRINCG(I,J))
         call inc_areg(I,J,JR,J_TRNFP0,-TNFS(3,I,J))
         call inc_areg(I,J,JR,J_TRNFP1,-TNFS(2,I,J))
         DO K=2,9
           JK=AJ_ALB_INDS(K-1) ! accumulate 8 radiation diags.
           DO IT=1,NTYPE
             call inc_aj(I,J,IT,JK,(S0*CSZ2)*ALB(I,J,K)*FTYPE(IT,I,J))
           END DO
           call inc_areg(I,J,JR,JK,(S0*CSZ2)*ALB(I,J,K))
         END DO
         AIJ(I,J,IJ_SRINCG) =AIJ(I,J,IJ_SRINCG) +(SRHR(0,I,J)*CSZ2/
     /        (ALB(I,J,1)+1.D-20))
         AIJ(I,J,IJ_SRNFG) =AIJ(I,J,IJ_SRNFG) +(SRHR(0,I,J)*CSZ2)
         AIJ(I,J,IJ_BTMPW) =AIJ(I,J,IJ_BTMPW) +BTMPW(I,J)
         AIJ(I,J,IJ_SRREF) =AIJ(I,J,IJ_SRREF) +S0*CSZ2*ALB(I,J,2)
         AIJ(I,J,IJ_SRVIS) =AIJ(I,J,IJ_SRVIS) +S0*CSZ2*ALB(I,J,4)
         AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)-TNFS(3,I,J)
         AIJ(I,J,IJ_SRNFP0)=AIJ(I,J,IJ_SRNFP0)+(SNFS(3,I,J)*CSZ2)
         AIJ(I,J,IJ_RNFP1) =AIJ(I,J,IJ_RNFP1) +(SNFS(2,I,J)*CSZ2
     *                                         -TNFS(2,I,J))
         AIJ(I,J,ij_srvdir)=AIJ(I,J,ij_srvdir)
     &        + FSRDIR(I,J)*SRVISSURF(I,J)
         AIJ(I,J,IJ_SRVISSURF)=AIJ(I,J,IJ_SRVISSURF)+SRVISSURF(I,J)
#ifdef mjo_subdd
         OLR_acc(I,J) = OLR_acc(I,J) -TNFS(3,I,J)
#endif
C**** CRF diags if required
         if (moddrf .ne. 0) go to 770
         if (cloud_rad_forc > 0) then
c    CRF diagnostics
           AIJ(I,J,IJ_SWCRF)=AIJ(I,J,IJ_SWCRF)+
     +          (SNFS(3,I,J)-SNFSCRF(I,J))*CSZ2
           AIJ(I,J,IJ_LWCRF)=AIJ(I,J,IJ_LWCRF)-
     -          (TNFS(3,I,J)-TNFSCRF(I,J))
          endif
         if (cloud_rad_forc.eq.2) then
c    CRF diagnostics without aerosols and Ox
           AIJ(I,J,IJ_SWCRF2)=AIJ(I,J,IJ_SWCRF2)+
     +          (SNFS(3,I,J)-SNFSCRF2(I,J))*CSZ2
           AIJ(I,J,IJ_LWCRF2)=AIJ(I,J,IJ_LWCRF2)-
     -          (TNFS(3,I,J)-TNFSCRF2(I,J))
         end if

C**** AERRF diags if required
         if (aer_rad_forc > 0) then
           do N=1,8
             AIJ(I,J,IJ_SWAERRF+N-1)=AIJ(I,J,IJ_SWAERRF+N-1)+
     *            (SNFS(3,I,J)-SNFSAERRF(N,I,J))*CSZ2
             AIJ(I,J,IJ_LWAERRF+N-1)=AIJ(I,J,IJ_LWAERRF+N-1)-
     *            (TNFS(3,I,J)-TNFSAERRF(N,I,J))
             AIJ(I,J,IJ_SWAERSRF+N-1)=AIJ(I,J,IJ_SWAERSRF+N-1)+
     *            (SNFS(1,I,J)-SNFSAERRF(N+8,I,J))*CSZ2
             AIJ(I,J,IJ_LWAERSRF+N-1)=AIJ(I,J,IJ_LWAERSRF+N-1)-
     *            (TNFS(1,I,J)-TNFSAERRF(N+8,I,J))
           end do
           AIJ(I,J,IJ_SWAERRFNT)=AIJ(I,J,IJ_SWAERRFNT)+
     *          (SNFS(3,I,J)-SNFSAERRF(17,I,J))*CSZ2
           AIJ(I,J,IJ_LWAERRFNT)=AIJ(I,J,IJ_LWAERRFNT)-
     *          (TNFS(3,I,J)-TNFSAERRF(17,I,J))
           AIJ(I,J,IJ_SWAERSRFNT)=AIJ(I,J,IJ_SWAERSRFNT)+
     *          (SNFS(1,I,J)-SNFSAERRF(18,I,J))*CSZ2
           AIJ(I,J,IJ_LWAERSRFNT)=AIJ(I,J,IJ_LWAERSRFNT)-
     *          (TNFS(1,I,J)-TNFSAERRF(18,I,J))
         end if

C***** Clear Sky and All Sky TOA Forcing without aerosol
           AIJ(I,J,IJ_SW_AS_noA)=AIJ(I,J,IJ_SW_AS_noA)+
     +          (SNFS(3,I,J)-SNFS_AS_noA(I,J))*CSZ2
           AIJ(I,J,IJ_LW_AS_noA)=AIJ(I,J,IJ_LW_AS_noA)-
     -          (TNFS(3,I,J)-TNFS_AS_noA(I,J))
           AIJ(I,J,IJ_SW_CS_noA)=AIJ(I,J,IJ_SW_CS_noA)+
     +          (SNFS(3,I,J)-SNFS_CS_noA(I,J))*CSZ2
           AIJ(I,J,IJ_LW_CS_noA)=AIJ(I,J,IJ_LW_CS_noA)-
     -          (TNFS(3,I,J)-TNFS_CS_noA(I,J))


#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS) ||\
    (defined TRACERS_AEROSOLS_SEASALT)
C**** Generic diagnostics for radiative forcing calculations
C**** Depending on whether tracers radiative interaction is turned on,
C**** diagnostic sign changes (for aerosols)
         rsign_aer=1. ; rsign_chem=-1.
         if (rad_interact_aer > 0) rsign_aer=-1.
C**** define SNFS/TNFS level (TOA/TROPO) for calculating forcing
         LFRC=3                 ! TOA
         if (rad_forc_lev.gt.0) LFRC=4 ! TROPOPAUSE
#ifdef BC_ALB
      if(ijts_alb(1).gt.0 .and. bc_snow_present(i,j) .and. csz2>0.) then
        TAIJS(I,J,ijts_sunlit_snow) = TAIJS(I,J,ijts_sunlit_snow) + 1.
        TAIJS(I,J,ijts_alb(1)) = TAIJS(I,J,ijts_alb(1)) + dALBsnBC(I,J)
!     &       + 100.d0*(ALBNBC(I,J)-ALB(I,J,1))
      endif
      if (ijts_alb(2).gt.0)
     & taijs(i,j,ijts_alb(2))
     &     =taijs(i,j,ijts_alb(2))
     &         +(SNFS(3,I,J)-NFSNBC(I,J))*CSZ2
#endif /* BC_ALB */
c     ..........
c     accumulation of forcings for tracers for which nraero_rf fields are
c     defined
c     ..........
           nsub_ntrix = 0
           do n=1,nraero_rf
               SELECT CASE (trname(ntrix_rf(n)))
               CASE ('Clay','ClayIlli','ClayKaol','ClaySmec' ,'ClayCalc'
     &                ,'ClayQuar','ClayFeld','ClayHema' ,'ClayGyps'
     &                ,'ClayIlHe','ClayKaHe','ClaySmHe' ,'ClayCaHe'
     &                ,'ClayQuHe','ClayFeHe','ClayGyHe')
                 nsub_ntrix(ntrix_rf(n)) = nsub_ntrix(ntrix_rf(n)) + 1
c shortwave forcing (TOA or TROPO) of Clay sub size classes
                 if (ijts_fcsub(1,ntrix_rf(n),
     &               nsub_ntrix(ntrix_rf(n))) > 0)
     &                taijs(i,j,ijts_fcsub(1,ntrix_rf(n)
     &                ,nsub_ntrix(ntrix_rf(n)))) =
     &               taijs(i,j,ijts_fcsub(1
     &                ,ntrix_rf(n),nsub_ntrix(ntrix_rf(n)))) +
     &               rsign_aer
     &                *(snfst(2,n,i,j)-snfs(lfrc,i,j))*csz2
c longwave forcing  (TOA or TROPO) of Clay size sub classes
                 if (ijts_fcsub(2,ntrix_rf(n),
     &               nsub_ntrix(ntrix_rf(n))) > 0)
     &                taijs(i,j,ijts_fcsub(2,ntrix_rf(n)
     &                ,nsub_ntrix(ntrix_rf(n)))) =
     &               taijs(i,j,ijts_fcsub(2
     &                ,ntrix_rf(n),nsub_ntrix(ntrix_rf(n)))) -
     &               rsign_aer
     &                *(tnfst(2,n,i,j)-tnfs(lfrc,i,j))
c shortwave forcing (TOA or TROPO) clear sky of Clay sub size classes
                 if (ijts_fcsub(5,ntrix_rf(n),
     &               nsub_ntrix(ntrix_rf(n))) > 0)
     &                taijs(i,j,ijts_fcsub(5,ntrix_rf(n)
     &                ,nsub_ntrix(ntrix_rf(n)))) =
     &               taijs(i,j,ijts_fcsub(5
     &                ,ntrix_rf(n),nsub_ntrix(ntrix_rf(n)))) +
     &               rsign_aer
     &                *(snfst(2,n,i,j)-snfs(lfrc,i,j))*csz2 * (1.D0
     &                -cfrac(i,j))
c longwave forcing  (TOA or TROPO) clear sky of Clay sub size classes
                 if (ijts_fcsub(6,ntrix_rf(n),
     &               nsub_ntrix(ntrix_rf(n))) > 0)
     &                taijs(i,j,ijts_fcsub(6,ntrix_rf(n)
     &                ,nsub_ntrix(ntrix_rf(n)))) =
     &               taijs(i,j,ijts_fcsub(6
     &                ,ntrix_rf(n),nsub_ntrix(ntrix_rf(n)))) -
     &               rsign_aer
     &                *(tnfst(2,n,i,j)-tnfs(lfrc,i,j)) * (1.D0-cfrac(i
     &                ,j))
c shortwave forcing at surface (if required) of Clay sub size classes
                 if (ijts_fcsub(3,ntrix_rf(n),
     &               nsub_ntrix(ntrix_rf(n))) > 0)
     &                taijs(i,j,ijts_fcsub(3,ntrix_rf(n)
     &                ,nsub_ntrix(ntrix_rf(n)))) =
     &               taijs(i,j,ijts_fcsub(3
     &                ,ntrix_rf(n),nsub_ntrix(ntrix_rf(n)))) +
     &               rsign_aer
     &                *(snfst(1,n,i,j)-snfs(1,i,j))*csz2
c longwave forcing at surface (if required) of Clay sub size classes
                 if (ijts_fcsub(4,ntrix_rf(n),
     &               nsub_ntrix(ntrix_rf(n))) > 0)
     &                taijs(i,j,ijts_fcsub(4,ntrix_rf(n)
     &                ,nsub_ntrix(ntrix_rf(n)))) =
     &               taijs(i,j,ijts_fcsub(4
     &                ,ntrix_rf(n),nsub_ntrix(ntrix_rf(n)))) -
     &               rsign_aer
     &                *(tnfst(1,n,i,j)-tnfs(1,i,j))
c shortwave forcing at surface clear sky (if required) of Clay sub size classes
                 if (ijts_fcsub(7,ntrix_rf(n),
     &               nsub_ntrix(ntrix_rf(n))) > 0)
     &                taijs(i,j,ijts_fcsub(7,ntrix_rf(n)
     &                ,nsub_ntrix(ntrix_rf(n)))) =
     &               taijs(i,j,ijts_fcsub(7
     &                ,ntrix_rf(n),nsub_ntrix(ntrix_rf(n)))) +
     &               rsign_aer
     &                *(snfst(1,n,i,j)-snfs(1,i,j))*csz2 * (1.D0-cfrac(i
     &                ,j))
c longwave forcing at surface clear sky (if required) of Clay sub size classes
                 if (ijts_fcsub(8,ntrix_rf(n),
     &               nsub_ntrix(ntrix_rf(n))) > 0)
     &                taijs(i,j,ijts_fcsub(8,ntrix_rf(n)
     &                ,nsub_ntrix(ntrix_rf(n)))) =
     &               taijs(i,j,ijts_fcsub(8
     &                ,ntrix_rf(n),nsub_ntrix(ntrix_rf(n)))) -
     &               rsign_aer
     &                *(tnfst(1,n,i,j)-tnfs(1,i,j)) * (1.D0-cfrac(i,j))
               CASE DEFAULT
                 SELECT CASE (trname(ntrix_rf(n)))
                 CASE ('seasalt2')
                   CYCLE
                 END SELECT
c shortwave forcing (TOA or TROPO)
                 if (ijts_fc(1,ntrix_rf(n)).gt.0)
     &                taijs(i,j,ijts_fc(1,ntrix_rf(n)))
     &                =taijs(i,j,ijts_fc(1,ntrix_rf(n)))
     &                +rsign_aer*(SNFST(2,N,I,J)-SNFS(LFRC,I,J))*CSZ2
c longwave forcing  (TOA or TROPO)
                 if (ijts_fc(2,ntrix_rf(n)).gt.0)
     &                taijs(i,j,ijts_fc(2,ntrix_rf(n)))
     &                =taijs(i,j,ijts_fc(2,ntrix_rf(n)))
     &                -rsign_aer*(TNFST(2,N,I,J)-TNFS(LFRC,I,J))
c shortwave forcing (TOA or TROPO) clear sky
                 if (ijts_fc(5,ntrix_rf(n)).gt.0)
     &                taijs(i,j,ijts_fc(5,ntrix_rf(n)))
     &                =taijs(i,j,ijts_fc(5,ntrix_rf(n)))
     &                +rsign_aer*(SNFST(2,N,I,J)-SNFS(LFRC,I,J))*CSZ2
     &                *(1.d0-CFRAC(I,J))
c longwave forcing  (TOA or TROPO) clear sky
                 if (ijts_fc(6,ntrix_rf(n)).gt.0)
     &                taijs(i,j,ijts_fc(6,ntrix_rf(n)))
     &                =taijs(i,j,ijts_fc(6,ntrix_rf(n)))
     &                -rsign_aer*(TNFST(2,N,I,J)-TNFS(LFRC,I,J))
     &                *(1.d0-CFRAC(I,J))
c shortwave forcing at surface (if required)
                 if (ijts_fc(3,ntrix_rf(n)).gt.0)
     &                taijs(i,j,ijts_fc(3,ntrix_rf(n)))
     &                =taijs(i,j,ijts_fc(3,ntrix_rf(n)))
     &                +rsign_aer*(SNFST(1,N,I,J)-SNFS(1,I,J))*CSZ2
c longwave forcing at surface (if required)
                 if (ijts_fc(4,ntrix_rf(n)).gt.0)
     &                taijs(i,j,ijts_fc(4,ntrix_rf(n)))
     &                =taijs(i,j,ijts_fc(4,ntrix_rf(n)))
     &                -rsign_aer*(TNFST(1,N,I,J)-TNFS(1,I,J))
c shortwave forcing at surface clear sky (if required)
                 if (ijts_fc(7,ntrix_rf(n)).gt.0)
     &                taijs(i,j,ijts_fc(7,ntrix_rf(n)))
     &                =taijs(i,j,ijts_fc(7,ntrix_rf(n)))
     &                +rsign_aer*(SNFST(1,N,I,J)-SNFS(1,I,J))*CSZ2
     &                *(1.d0-CFRAC(I,J))
c longwave forcing at surface clear sky (if required)
                 if (ijts_fc(8,ntrix_rf(n)).gt.0)
     &                taijs(i,j,ijts_fc(8,ntrix_rf(n)))
     &                =taijs(i,j,ijts_fc(8,ntrix_rf(n)))
     &                -rsign_aer*(TNFST(1,N,I,J)-TNFS(1,I,J))
     &                *(1.d0-CFRAC(I,J))
               END SELECT
           end do     ! n=1,nraero_rf

c ..........
c accumulation of forcings for special case ozone (nraero_rf fields
c not defined) Warning: indicies used differently, since we don't
c need CS or Surface, but are doing both TOA and Ltropo:
c ..........
         if(n_Ox > 0) then ! ------ main Ox tracer -------
c shortwave forcing at tropopause
           if (ijts_fc(1,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(1,n_Ox))=taijs(i,j,ijts_fc(1,n_Ox))
     &     +rsign_chem*(SNFST_o3ref(1,I,J)-SNFS(4,I,J))*CSZ2
c longwave forcing at tropopause
           if (ijts_fc(2,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(2,n_Ox))=taijs(i,j,ijts_fc(2,n_Ox))
     &     -rsign_chem*(TNFST_o3ref(1,I,J)-TNFS(4,I,J))
c shortwave forcing at TOA
           if (ijts_fc(3,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(3,n_Ox))=taijs(i,j,ijts_fc(3,n_Ox))
     &     +rsign_chem*(SNFST_o3ref(2,I,J)-SNFS(3,I,J))*CSZ2
c longwave forcing at TOA
           if (ijts_fc(4,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(4,n_Ox))=taijs(i,j,ijts_fc(4,n_Ox))
     &     -rsign_chem*(TNFST_o3ref(2,I,J)-TNFS(3,I,J))
         endif
#ifdef AUXILIARY_OX_RADF
c shortwave forcing at tropopause
         if (ijts_auxfc(1)>0)
     &   taijs(i,j,ijts_auxfc(1))=taijs(i,j,ijts_auxfc(1))
#ifdef AUX_OX_RADF_TROP
     &   +rsign_chem*(SNFST_o3ref(5,I,J)-SNFST_o3ref(3,I,J))*CSZ2
#else
     &   +rsign_chem*(SNFST_o3ref(1,I,J)-SNFST_o3ref(3,I,J))*CSZ2
#endif
c longwave forcing at tropopause
         if (ijts_auxfc(2)>0)
     &   taijs(i,j,ijts_auxfc(2))=taijs(i,j,ijts_auxfc(2))
#ifdef AUX_OX_RADF_TROP
     &   -rsign_chem*(TNFST_o3ref(5,I,J)-TNFST_o3ref(3,I,J))
#else
     &   -rsign_chem*(TNFST_o3ref(1,I,J)-TNFST_o3ref(3,I,J))
#endif
c shortwave forcing at TOA
         if (ijts_auxfc(3)>0)
     &   taijs(i,j,ijts_auxfc(3))=taijs(i,j,ijts_auxfc(3))
     &   +rsign_chem*(SNFST_o3ref(2,I,J)-SNFST_o3ref(4,I,J))*CSZ2
c longwave forcing at TOA
         if (ijts_auxfc(4)>0)
     &   taijs(i,j,ijts_auxfc(4))=taijs(i,j,ijts_auxfc(4))
     &   -rsign_chem*(TNFST_o3ref(2,I,J)-TNFST_o3ref(4,I,J))
#endif /* AUXILIARY_OX_RADF */
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
                      ! ------ diag stratOx tracer -------
! note for now for this diag, there is a failsafe that stops model
! if clim_interact_chem .le. 0 when the below would be wrong:
c shortwave forcing at tropopause
         if (ijts_fc(1,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(1,n_stratOx))=taijs(i,j,ijts_fc(1,n_stratOx))
     &   +rsign_chem*(SNFST_o3ref(1,I,J)-SNFST_stratOx(1,I,J))*CSZ2
c longwave forcing at tropopause
         if (ijts_fc(2,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(2,n_stratOx))=taijs(i,j,ijts_fc(2,n_stratOx))
     &   -rsign_chem*(TNFST_o3ref(1,I,J)-TNFST_stratOx(1,I,J))
c shortwave forcing at TOA
         if (ijts_fc(3,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(3,n_stratOx))=taijs(i,j,ijts_fc(3,n_stratOx))
     &   +rsign_chem*(SNFST_o3ref(2,I,J)-SNFST_stratOx(2,I,J))*CSZ2
c longwave forcing at TOA
         if (ijts_fc(4,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(4,n_stratOx))=taijs(i,j,ijts_fc(4,n_stratOx))
     &   -rsign_chem*(TNFST_o3ref(2,I,J)-TNFST_stratOx(2,I,J))
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS*/
#endif /* any of various tracer groups defined */

#ifdef ACCMIP_LIKE_DIAGS
#ifndef SKIP_ACCMIP_GHG_RADF_DIAGS
         do nf=1,4 ! CH4, N2O, CFC11, and CFC12:
c shortwave GHG forcing at TOA
           if(ij_fcghg(1,nf).gt.0)aij(i,j,ij_fcghg(1,nf))=
     &     aij(i,j,ij_fcghg(1,nf))+(SNFS(3,I,J)-SNFS_ghg(nf,I,J))
     &     *CSZ2
c longwave GHG forcing at TOA
           if(ij_fcghg(2,nf).gt.0)aij(i,j,ij_fcghg(2,nf))=
     &     aij(i,j,ij_fcghg(2,nf))+(TNFS_ghg(nf,I,J)-TNFS(3,I,J))
         enddo
#endif /* NOT DEFINED SKIP_ACCMIP_GHG_RADF_DIAGS */
#endif /* ACCMIP_LIKE_DIAGS */

#ifdef CACHED_SUBDD
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS) ||\
    (defined TRACERS_AEROSOLS_SEASALT)
      if (nraero_rf>0) then
        swfrc(i,j,1:nraero_rf)=
     &    rsign_aer*(SNFST(2,1:nraero_rf,I,J)-SNFS(LFRC,I,J))*CSZ2
        lwfrc(i,j,1:nraero_rf)=
     &    -rsign_aer*(TNFST(2,1:nraero_rf,I,J)-TNFS(LFRC,I,J))
      endif
#endif /* any of various tracer groups defined */
#endif  /* CACHED_SUBDD */

  770    CONTINUE
  780    CONTINUE

#ifdef mjo_subdd
         OLR_cnt = OLR_cnt + 1.
#endif

      DO J=J_0,J_1
      DO I=I_0,I_1
      DO L=1,LM
        AIJL(i,j,l,IJL_RC)=AIJL(i,j,l,IJL_RC)+
     &       (SRHR(L,I,J)*COSZ2(I,J)+TRHR(L,I,J))
#ifdef mjo_subdd
       SWHR(I,J,L)=SWHR(I,J,L)+
     *             SRHR(L,I,J)*COSZ2(I,J)*bysha*byMA(L,I,J)
       LWHR(I,J,L)=LWHR(I,J,L)+TRHR(L,I,J)*bysha*byMA(L,I,J)
#endif
      END DO
      END DO
      END DO
#ifdef mjo_subdd
      SWHR_cnt=SWHR_cnt+1
      LWHR_cnt=LWHR_cnt+1
#endif

#ifdef CACHED_SUBDD
      do k=1,subdd_ngroups
        subdd => subdd_groups(k)
        subdd%nacc(subdd%subdd_period,sched_rad) =
     &  subdd%nacc(subdd%subdd_period,sched_rad) + 1
      enddo
C****
C**** Collect some high-frequency outputs
C****
      call find_groups('rijh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('olrrad')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j) = tnfs(3,i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('olrcs')
        if(cloud_rad_forc.le.0.) call stop_model(
     &       'diagnostic olrcs needs cloud_rad_forc>0',255)
        call inc_subdd(subdd,k,TNFSCRF)
      case ('lwds')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j) = trhr(0,i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('lwdscs')
        if(cloud_rad_forc.le.0.) call stop_model(
     &       'diagnostic lwdscs needs cloud_rad_forc>0',255)
        call inc_subdd(subdd,k,lwdncs)
      case ('lwus')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j) = trhr(0,i,j) + tnfs(1,i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('swus')
        call inc_subdd(subdd,k,SWUS)
      case ('swds')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j)=srdn(i,j)*cosz2(i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('swdf')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j)=fsrdif(i,j)+difnir(i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('swtoa')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j) = snfs(3,i,j)*cosz2(i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      !Net solar flux at surface:
      case ('swns')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j)=snfs(1,i,j)*cosz2(i,j)
        enddo;       enddo
        call inc_subdd(subdd,k,sddarr)
      !Net Longwave flux at surface:
      case ('lwns')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j)=tnfs(1,i,j)
        enddo;       enddo
        call inc_subdd(subdd,k,sddarr)
      case ('totcld')
        call inc_subdd(subdd,k,cfrac)
      case ('totcld_diag')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          call get_cld_overlap(lm,cldss(:,i,j),
     &                         cldmcl=cldmc(:,i,j),
     &                         CldTot=sddarr(i,j))
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('cldss_2d')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          call get_cld_overlap(lm,cldss(:,i,j),
     &                         CldSS=sddarr(i,j))
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('cldmc_2d')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          call get_cld_overlap(lm,cldss(:,i,j),
     &                         cldmcl=cldmc(:,i,j),
     &                         CldMC=sddarr(i,j))
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('wtrcld')
        call inc_subdd(subdd,k,WTRCLD)
      case ('icecld')
        call inc_subdd(subdd,k,ICECLD)
      case ('cod')
        call inc_subdd(subdd,k,TAUSUMW)
      case ('cid')
        call inc_subdd(subdd,k,TAUSUMI)
      case ('ctp')
        call inc_subdd(subdd,k,CTP)
      case ('ctt')
        call inc_subdd(subdd,k,CTT)
#ifdef CFMIP3_SUBDD
      case ('rtmt')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j) = (snfs(3,i,j)*cosz2(i,j))-tnfs(2,i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('swut')
        call inc_subdd(subdd,k,swut)
      case ('swutcs')
        call inc_subdd(subdd,k,swutcs)
      case ('clwvi')
        call inc_subdd(subdd,k,cfmip_twp)
      case ('swdcls')
        call inc_subdd(subdd,k,swdcls)
      case('swucls')
        call inc_subdd(subdd,k,swucls)
      case('swdt')
        call inc_subdd(subdd,k,swdt)
#endif
      end select

      enddo
      enddo

      call find_groups('rijlh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('MRCO2rad')
        do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
          sddarr3d(i,j,l) = CO2out(l,i,j)
        enddo;        enddo;             enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('wtrtau')
        call inc_subdd(subdd,k,wtrtau)
      case ('icetau')
        call inc_subdd(subdd,k,icetau)
      end select
      end do
      end do

#ifdef SCM
      call find_groups('rijlh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('dth_sw')
        do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
          sddarr3d(i,j,l) = SRHR(L,I,J)*bysha*byma(L,I,J)*COSZ2(I,J)/
     &                      PK(L,I,J)
        enddo;                enddo;    enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('dth_lw')
        do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
          sddarr3d(i,j,l) = TRHR(L,I,J)*bysha*byma(L,I,J)/PK(L,I,J)
        enddo;           enddo;         enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('dth_rad')
        do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
          sddarr3d(i,j,l) = (SRHR(L,I,J)*COSZ2(I,J)+TRHR(L,I,J))*
     &                      bysha*byma(L,I,J)/PK(L,I,J)
        enddo;           enddo;         enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('lwdp')
        do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
          sddarr3d(i,j,l) = TRDFLB_prof(i,j,l)
        enddo;           enddo;         enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('lwup')
        do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
          sddarr3d(i,j,l) = TRUFLB_prof(i,j,l)
        enddo;           enddo;         enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('swdp')
        do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
          sddarr3d(i,j,l) = SRDFLB_prof(i,j,l)*COSZ2(I,J)
        enddo;           enddo;         enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('swup')
        do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
          sddarr3d(i,j,l) = SRUFLB_prof(i,j,l)*COSZ2(I,J)
        enddo;           enddo;         enddo
        call inc_subdd(subdd,k,sddarr3d)
      end select
      enddo
      enddo
#endif
#ifdef CFMIP3_SUBDD
      call find_groups('rijlh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('cf')
        call inc_subdd(subdd,k,cfmip_cf)
      case ('qcirad')
        call inc_subdd(subdd,k,cfmip_qci)
      case ('qclrad')
        call inc_subdd(subdd,k,cfmip_qcl)
      end select
      enddo
      enddo
#endif
#ifdef TRACERS_ON

! aod
      do g=1,size(sgroups)
      call find_groups(sgroups(g),grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      do s=1,size(ssky)
      if (trim(ssky(s)).eq.'dry' .and. save_dry_aod==0) cycle
      do a=1,size(sabs)
        select case (trim(ssky(s))//trim(sabs(a)))
          case ('as')
            sddarr4d=tau_as
          case ('cs')
            sddarr4d=tau_cs
          case ('dry')
            sddarr4d=tau_dry
          case ('asa')
            sddarr4d=abstau_as
          case ('csa')
            sddarr4d=abstau_cs
          case ('drya')
            sddarr4d=abstau_dry
          case default
            cycle ! not implemented, silently ignore
        end select
        do n=1,nraero_aod+1 ! +1 for total
          if (n<=nraero_aod) then
            spcname = trim(trname(ntrix_aod(n)))
          else
            spcname = ''
          endif
          !aod
          sname = trim(spcname)//trim(ssky(s))//trim(sabs(a))//'aod'
          if (trim(sgroups(g))=='taijlh') sname=trim(sname)//'3d'
          if (trim(sname)==trim(subdd%name(k))) then ! not select case here
            if (n<=nraero_aod) then
              sddarr3d=sddarr4d(:,:,:,n)
            else
              sddarr3d=sum(sddarr4d,dim=4)
            endif
            select case (trim(sgroups(g)))
              case ('taijh')
                sddarr=sum(sddarr3d,dim=3)
                call inc_subdd(subdd,k,sddarr)
              case ('taijlh')
                call inc_subdd(subdd,k,sddarr3d)
            end select
          endif
          !bext (bcoef) or babs (abcoef)
          if (trim(sgroups(g))=='taijlh') then
            sname = trim(spcname)//trim(ssky(s))//trim(sabs(a))//
     *              'bcoef3d'
            if (trim(sname)==trim(subdd%name(k))) then ! not select case here
              if (n<=nraero_aod) then
                sddarr3d=sddarr4d(:,:,:,n)
              else
                sddarr3d=sum(sddarr4d,dim=4)
              endif
              do j=j_0,j_1
                 do i=i_0,imaxj(j)
                    do l=1,lm
                       tlm(l) = T(i,j,l)*pk(l,i,j)
                       rho = pmid(l,i,j)*100./(Rgas*tlm(l))
                       dz = ma(l,i,j)/rho
                       sddarr3d(i,j,l) = sddarr3d(i,j,l)/dz
                    enddo
                 enddo
              enddo
              call inc_subdd(subdd,k,sddarr3d)
            endif
          endif
        enddo ! n
      enddo ! a
      enddo ! s
      enddo ! k
      enddo ! igrp
      enddo ! g

! rf
      call find_groups('taijh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      do f=1,size(sfrc)
        select case (trim(sfrc(f)))
          case ('swf')
            sddarr3drf=swfrc
          case ('lwf')
            sddarr3drf=lwfrc
          case default
            cycle ! not implemented, silently ignore
        end select
        do n=1,nraero_rf
          if (diag_fc==2) then
            spcname = trim(trname(ntrix_rf(n)))
          else if (diag_fc==1) then
            if (tracers_amp) then
              spcname='AMP'
            elseif (tracers_tomas) then
              spcname='TOMAS'
            else
              spcname='OMA'
            endif
          endif
          sname = trim(sfrc(f))//'_'//trim(spcname)
          if (trim(sname)==trim(subdd%name(k))) then ! not select case here
            call inc_subdd(subdd,k,sddarr3drf(:,:,n))
          endif
        enddo ! n
      enddo ! f
      enddo ! k
      enddo ! igrp

#endif  /* TRACERS_ON */

#endif  /* CACHED_SUBDD */

C****
C**** Update radiative equilibrium temperatures
C****
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          DO LR=1,LM_REQ
            RQT(LR,I,J)=RQT(LR,I,J)+(SRHRS(LR,I,J)*COSZ2(I,J)
     &           +TRHRS(LR,I,J))*NRAD*DTsrc*bysha*byaml00(lr+lm)
          END DO
        END DO
      END DO
C****
C**** Update other temperatures every physics time step
C****
  900 DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          DO L=1,LM
            T(I,J,L)=T(I,J,L)+(SRHR(L,I,J)*COSZ1(I,J)+TRHR(L,I,J))*
     *           DTsrc*bysha*byMA(l,i,j)/PK(L,I,J)
          END DO
          AIJ(I,J,IJ_SRINCP0)=AIJ(I,J,IJ_SRINCP0)+(S0*COSZ1(I,J))
        END DO
      END DO

C**** daily diagnostics
      IH=1+modelEclock%getHour()
      IHM = IH+(modelEclock%getDate()-1)*24
      DO KR=1,NDIUPT
        I = IJDD(1,KR)
        J = IJDD(2,KR)
        IF ((J >= J_0) .AND. (J <= J_1) .AND.
     &      (I >= I_0) .AND. (I <= I_1)) THEN
          ADIURN(IDD_ISW,KR,IH)=ADIURN(IDD_ISW,KR,IH)+S0*COSZ1(I,J)
#ifdef USE_HDIURN
          HDIURN(IDD_ISW,KR,IHM)=HDIURN(IDD_ISW,KR,IHM)+S0*COSZ1(I,J)
#endif
        ENDIF
      ENDDO

      call stopTimer('RADIA()')
      RETURN
      END SUBROUTINE RADIA

      SUBROUTINE RESET_SURF_FLUXES(I,J,ITYPE_OLD,ITYPE_NEW,FTYPE_ORIG,
     *     FTYPE_NOW)
!@sum set incident solar and upward thermal fluxes appropriately
!@+   as fractions change to conserve energy, prevent restart problems
!@auth Gavin Schmidt
      use rad_com, only : fsf,trsurf
      implicit none
!@var itype_old, itype_new indices for the old type turning to new type
      integer, intent(in) :: i,j,itype_old,itype_new
!@var ftype_orig, ftype_now original and current fracs of the 'new' type
      real*8, intent(in) :: ftype_orig, ftype_now
      real*8 :: delf ! change in fraction from old to new

      if (( (ITYPE_OLD==1 .and. ITYPE_NEW==2)
     *  .or.(ITYPE_OLD==2 .and. ITYPE_NEW==1))
     *  .and. (FTYPE_NOW .le. 0. .or. FTYPE_NOW .gt. 1.)) then
        write (6,*) ' RESET_SURF_FLUXES:
     *    I, J, ITYPE_OLD, ITYPE_NEW, FTYPE_ORIG, FTYPE_NOW = ',
     *    I, J, ITYPE_OLD, ITYPE_NEW, FTYPE_ORIG, FTYPE_NOW
        call stop_model('RESET_SURF_FLUXES: INCORRECT RESET',255)
      end if

      delf = FTYPE_NOW-FTYPE_ORIG
C**** Constrain fsf_1*ftype_1+fsf_2*ftype_2 to be constant
      FSF(ITYPE_NEW,I,J)=(FSF(ITYPE_NEW,I,J)*FTYPE_ORIG+
     *     FSF(ITYPE_OLD,I,J)*DELF)/FTYPE_NOW

C**** Same for upward thermal
      TRSURF(ITYPE_NEW,I,J)=(TRSURF(ITYPE_NEW,I,J)*FTYPE_ORIG+
     *     TRSURF(ITYPE_OLD,I,J)*DELF)/FTYPE_NOW

      RETURN
      END SUBROUTINE RESET_SURF_FLUXES

      SUBROUTINE GHGHST(iu)
!@sum  reads history for nghg well-mixed greenhouse gases
!@auth R. Ruedy

      use domain_decomp_atm, only : write_parallel
      USE RADPAR, only : nghg,ghgyr1,ghgyr2,ghgam
      USE RAD_COM, only : ghg_yr
      IMPLICIT NONE
      INTEGER :: iu,n,k,nhead=4,iyr
      CHARACTER*80 title
      character(len=300) :: out_line

      write(out_line,*)  ! print header lines and first data line
      call write_parallel(trim(out_line),unit=6)
      do n=1,nhead+1
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
      end do
      if(title(1:2).eq.'--') then                 ! older format
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
        nhead=5
      end if

!**** find range of table: ghgyr1 - ghgyr2
      read(title,*) ghgyr1
      do ; read(iu,'(a)',end=20) title ; end do
   20 read(title,*) ghgyr2
      rewind iu  !   position to data lines
      do n=1,nhead ; read(iu,'(a)') ; end do

      allocate (ghgam(nghg,ghgyr2-ghgyr1+1))
      do n=1,ghgyr2-ghgyr1+1
        read(iu,*) iyr,(ghgam(k,n),k=1,nghg)
        do k=1,nghg ! replace -999. by reasonable numbers
          if(ghgam(k,n).lt.0.) ghgam(k,n)=ghgam(k,n-1)
        end do
        if(ghg_yr>0 .and. abs(ghg_yr-iyr).le.1) then
          write(out_line,'(i5,6f10.4)') iyr,(ghgam(k,n),k=1,nghg)
          call write_parallel(trim(out_line),unit=6)
        endif
      end do
      write(out_line,*) 'read GHG table for years',ghgyr1,' - ',ghgyr2
      call write_parallel(trim(out_line),unit=6)
      return
      end SUBROUTINE GHGHST

#if defined(CUBED_SPHERE)
      subroutine read_qma (iu,plb)
!@sum  reads H2O production rates induced by CH4 (Tim Hall)
!@auth R. Ruedy
      use domain_decomp_atm, only : write_parallel
      use rad_com, only : dH2O,jma=>jm_dh2o,lat_dh2o
      use resolution, only : lm
      use constant, only : radian
      use TimeConstants_mod, only: DAYS_PER_YEAR
      implicit none
      integer, parameter:: lma=24
      integer m,iu,j,l,ll,ldn(lm),lup(lm)
      real*8 :: plb(lm+1)
      real*4 pb(0:lma+1),h2o(jma,0:lma),z(lma),dz(0:lma)
      character*100 title
      real*4 pdn,pup,dh,fracl
      character(len=300) :: out_line

C**** read headers/latitudes
      read(iu,'(a)') title
      write(out_line,'(''0'',a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
      write(out_line,'(1x,a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
c      write(6,'(1x,a100)') title
      read(title(10:100),*) (lat_dh2o(j),j=1,jma)
      lat_dh2o(:) = lat_dh2o(:)*radian

C**** read heights z(km) and data (kg/km^3/year)
      do m=1,12
        read(iu,'(a)') title
        write(out_line,'(1x,a100)') title
        call write_parallel(trim(out_line),unit=6)
        do l=lma,1,-1
          read(iu,'(a)') title
c          write(6,'(1x,a100)') title
          read(title,*) z(l),(H2O(j,l),j=1,jma)
        end do
        do j=1,jma
          h2o(j,0) = 0.
        end do

C**** Find edge heights and pressures
        dz(0) = 0.
        dz(1) = z(2)-z(1)
        do l=2,lma-1
           dz(l)=.5*(z(l+1)-z(l-1))
        end do
        dz(lma) = z(lma)-z(lma-1)

        pb(0) = plb(1)
        do l=1,lma
           Pb(l)=1000.*10.**(-(z(l)-.5*dz(l))/16.)
        end do
C**** extend both systems vertically to p=0
        pb(lma+1)=0.
        plb(lm+1)=0.

C**** Interpolate vertical resolution to model layers
        ldn(:) = 0
        do l=1,lm
          do while (pb(ldn(l)+1).ge.plb(l) .and. ldn(l).lt.lma)
            ldn(l)=ldn(l)+1
          end do
          lup(l)=ldn(l)
          do while (pb(lup(l)+1).gt.plb(l+1) .and. lup(l).lt.lma)
            lup(l)=lup(l)+1
          end do
        end do

C**** Interpolate (extrapolate) vertically
        do j=1,jma
          do l=1,lm
            dh = 0.
            pdn = plb(l)
            if (lup(l).gt.0) then
              do ll=ldn(l),lup(l)
                pup = max(REAL(pb(ll+1),KIND=8),plb(l+1))
                fracl= (pdn-pup)/(pb(ll)-pb(ll+1))
                dh = dh+h2o(j,ll)*fracl*dz(ll)
                pdn = pup
              end do
            end if
            dh2o(j,l,m) = 1.d-6*dh/1.74d0/DAYS_PER_YEAR !->(kg/m^2/ppm_CH4/day)
          end do
        end do
      end do
      return
      end subroutine read_qma

      subroutine lat_interp_qma (rlat,lev,mon,dh2o_interp)
!@sum  interpolate CH4->H2O production rates in latitude
!@auth R. Ruedy
      use rad_com, only : jma=>jm_dh2o,xlat=>lat_dh2o,dh2o
      implicit none
      real*8 :: rlat ! input latitude (radians)
      integer :: lev,mon ! input level, month
      real*8 :: dh2o_interp  ! output
      real*8 w1,w2
      integer :: j1,j2

C**** Interpolate (extrapolate) horizontally
      j2 = 2+(jma-1)*(rlat-xlat(1))/(xlat(jma)-xlat(1)) ! first guess
      j2 = min(max(2,j2),jma)
      j1 = j2-1
      if(rlat.gt.xlat(j2)) then ! j guess was too low
         do while (j2.lt.jma .and. rlat.gt.xlat(j2))
          j2 = j2+1
         end do
         j1 = j2-1
      elseif(rlat.lt.xlat(j1)) then ! j guess was too high
         do while (j1.gt.1 .and. rlat.lt.xlat(j1))
          j1 = j1-1
         end do
         j2 = j1+1
      endif
C**** coeff. for latitudinal linear inter/extrapolation
      w1 = (xlat(j2)-rlat)/(xlat(j2)-xlat(j1))
C**** for extrapolations, only use half the slope
      if(w1.gt.1.) w1=.5+.5*w1
      if(w1.lt.0.) w1=.5*w1
      w2 = 1.-w1
      dh2o_interp = w1*dh2o(j1,lev,mon)+w2*dh2o(j2,lev,mon)
      return
      end subroutine lat_interp_qma

#endif /* CUBED_SPHERE */

      subroutine getqma (iu,dglat,plb,dh2o,lm,jm)
!@sum  reads H2O production rates induced by CH4 (Tim Hall)
!@auth R. Ruedy
      use domain_decomp_atm, only : grid,getDomainBounds,write_parallel
      use TimeConstants_mod, only: DAYS_PER_YEAR
      implicit none
      integer, parameter:: jma=18,lma=24
      integer m,iu,jm,lm,j,j1,j2,l,ll,ldn(lm),lup(lm)
      real*8 PLB(lm+1),dH2O(grid%j_strt_halo:grid%j_stop_halo,lm,12)
     &     ,dglat(jm)
      real*4 pb(0:lma+1),h2o(jma,0:lma),xlat(jma),z(lma),dz(0:lma)
      character*100 title
      real*4 pdn,pup,w1,w2,dh,fracl
      integer :: j_0,j_1
      character(len=300) :: out_line
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

C**** read headers/latitudes
      read(iu,'(a)') title
      write(out_line,'(''0'',a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
      write(out_line,'(1x,a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
c      write(6,'(1x,a100)') title
      read(title(10:100),*) (xlat(j),j=1,jma)

C**** read heights z(km) and data (kg/km^3/year)
      do m=1,12
        read(iu,'(a)') title
        write(out_line,'(1x,a100)') title
        call write_parallel(trim(out_line),unit=6)
        do l=lma,1,-1
          read(iu,'(a)') title
c          write(6,'(1x,a100)') title
          read(title,*) z(l),(H2O(j,l),j=1,jma)
        end do
        do j=1,jma
          h2o(j,0) = 0.
        end do

C**** Find edge heights and pressures
        dz(0) = 0.
        dz(1) = z(2)-z(1)
        do l=2,lma-1
           dz(l)=.5*(z(l+1)-z(l-1))
        end do
        dz(lma) = z(lma)-z(lma-1)

        pb(0) = plb(1)
        do l=1,lma
           Pb(l)=1000.*10.**(-(z(l)-.5*dz(l))/16.)
        end do
C**** extend both systems vertically to p=0
        pb(lma+1)=0.
        plb(lm+1)=0.

C**** Interpolate vertical resolution to model layers
        ldn(:) = 0
        do l=1,lm
          do while (pb(ldn(l)+1).ge.plb(l) .and. ldn(l).lt.lma)
            ldn(l)=ldn(l)+1
          end do
          lup(l)=ldn(l)
          do while (pb(lup(l)+1).gt.plb(l+1) .and. lup(l).lt.lma)
            lup(l)=lup(l)+1
          end do
        end do

C**** Interpolate (extrapolate) horizontally and vertically
        j2=2
        do j=j_0,j_1
C**** coeff. for latitudinal linear inter/extrapolation
          do while (j2.lt.jma .and. dglat(j).gt.xlat(j2))
            j2 = j2+1
          end do
          j1 = j2-1
          w1 = (xlat(j2)-dglat(j))/(xlat(j2)-xlat(j1))
C**** for extrapolations, only use half the slope
          if(w1.gt.1.) w1=.5+.5*w1
          if(w1.lt.0.) w1=.5*w1
          w2 = 1.-w1
          do l=1,lm
            dh = 0.
            pdn = plb(l)
            if (lup(l).gt.0) then
              do ll=ldn(l),lup(l)
                pup = max(REAL(pb(ll+1),KIND=8),plb(l+1))
                fracl= (pdn-pup)/(pb(ll)-pb(ll+1))
                dh = dh+(w1*h2o(j1,ll)+w2*h2o(j2,ll))*fracl*dz(ll)
                pdn = pup
              end do
            end if
            dh2o(j,l,m) = 1.d-6*dh/1.74d0/DAYS_PER_YEAR !->(kg/m^2/ppm_CH4/day)
          end do
        end do
      end do
      return
      end subroutine getqma


      SUBROUTINE ORBIT (DOBLIQ,ECCEN,DOMEGVP,VEDAY,EDPY, DAY,
     *                  SDIST,SIND,COSD,SUNLON,SUNLAT,EQTIME)
C****
C**** ORBIT receives orbital parameters and time of year, and returns
C**** distance from Sun, declination angle, and Sun's overhead position.
C**** Reference for following caculations is:  V.M.Blanco and
C**** S.W.McCuskey, 1961, "Basic Physics of the Solar System", pages
C**** 135 - 151.  Existence of Moon and heavenly bodies other than
C**** Earth and Sun are ignored.  Earth is assumed to be spherical.
C****
C**** Program author: Gary L. Russell 2004/11/16
C**** Angles, longitude and latitude are measured in radians.
C****
C**** Input: ECCEN  = eccentricity of the orbital ellipse
C****        OBLIQ  = latitude of Tropic of Cancer
C****        OMEGVP = longitude of perihelion (sometimes Pi is added) =
C****               = spatial angle from vernal equinox to perihelion
C****                 with Sun as angle vertex
C****        DAY    = days measured since 2000 January 1, hour 0
C****
C****        EDPY  = Earth days per year
C****                tropical year = 365.2425 (Gregorgian Calendar)
C****                tropical year = 365      (Generic Year)
C****        VEDAY = Vernal equinox
C****                79.0 (Generic year Mar 21 hour 0)
C****                79.5 (Generic year Mar 21 hour 12 - PMIP standard)
C****                79.3125d0 for days from 2000 January 1, hour 0 till vernal
C****                     equinox of year 2000 = 31 + 29 + 19 + 7.5/24
C****
C**** Intermediate quantities:
C****    BSEMI = semi minor axis in units of semi major axis
C****   PERIHE = perihelion in days since 2000 January 1, hour 0
C****            in its annual revolution about Sun
C****       TA = true anomaly = spatial angle from perihelion to
C****            current location with Sun as angle vertex
C****       EA = eccentric anomaly = spatial angle measured along
C****            eccentric circle (that circumscribes Earth's orbit)
C****            from perihelion to point above (or below) Earth's
C****            absisca (where absisca is directed from center of
C****            eccentric circle to perihelion)
C****       MA = mean anomaly = temporal angle from perihelion to
C****            current time in units of 2*Pi per tropical year
C****   TAofVE = TA(VE) = true anomaly of vernal equinox = - OMEGVP
C****   EAofVE = EA(VE) = eccentric anomaly of vernal equinox
C****   MAofVE = MA(VE) = mean anomaly of vernal equinox
C****   SLNORO = longitude of Sun in Earth's nonrotating reference frame
C****   VEQLON = longitude of Greenwich Meridion in Earth's nonrotating
C****            reference frame at vernal equinox
C****   ROTATE = change in longitude in Earth's nonrotating reference
C****            frame from point's location on vernal equinox to its
C****            current location where point is fixed on rotating Earth
C****   SLMEAN = longitude of fictitious mean Sun in Earth's rotating
C****            reference frame (normal longitude and latitude)
C****
C**** Output: SIND = sine of declination angle = sin(SUNLAT)
C****         COSD = cosine of the declination angle = cos(SUNLAT)
C****       SUNDIS = distance to Sun in units of semi major axis
C****       SUNLON = longitude of point on Earth directly beneath Sun
C****       SUNLAT = latitude of point on Earth directly beneath Sun
C****       EQTIME = Equation of Time =
C****              = longitude of fictitious mean Sun minus SUNLON
C****
C**** From the above reference:
C**** (4-54): [1 - ECCEN*cos(EA)]*[1 + ECCEN*cos(TA)] = (1 - ECCEN^2)
C**** (4-55): tan(TA/2) = sqrt[(1+ECCEN)/(1-ECCEN)]*tan(EA/2)
C**** Yield:  tan(EA) = sin(TA)*sqrt(1-ECCEN^2) / [cos(TA) + ECCEN]
C****    or:  tan(TA) = sin(EA)*sqrt(1-ECCEN^2) / [cos(EA) - ECCEN]
C****
      USE CONSTANT, only : twopi,pi,radian
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DOBLIQ,ECCEN,DOMEGVP,DAY,VEDAY,EDPY
      REAL*8, INTENT(OUT) :: SIND,COSD,SDIST,SUNLON,SUNLAT,EQTIME

      REAL*8 MA,OMEGVP,OBLIQ,EA,DEA,BSEMI
     *     ,TAofVE,EAofVE,MAofVE,SUNDIS,TA,SUNX,SUNY,SLNORO
     *     ,VEQLON,ROTATE,SLMEAN
c      REAL*8, PARAMETER :: EDAYzY=365.2425d0, VE2000=79.3125d0
c      REAL*8, PARAMETER :: EDAYzY=365d0, VE2000=79d0  ! original parameters
      REAL*8  EDAYzY,VE2000
C****
      VE2000=VEDAY
      EDAYzY=EDPY
      OMEGVP=DOMEGVP*radian
      OBLIQ=DOBLIQ*radian
C**** Determine EAofVE from geometry: tan(EA) = b*sin(TA) / [e+cos(TA)]
C**** Determine MAofVE from Kepler's equation: MA = EA - e*sin(EA)
C**** Determine MA knowing time from vernal equinox to current day
C****
      BSEMI  = SQRT (1 - ECCEN*ECCEN)
      TAofVE = - OMEGVP
      EAofVE = ATAN2 (BSEMI*SIN(TAofVE), ECCEN+COS(TAofVE))
      MAofVE = EAofVE - ECCEN*SIN(EAofVE)
C     PERIHE = VE2000 - MAofVE*EDAYzY/TWOPI
      MA     = MODULO (TWOPI*(DAY-VE2000)/EDAYzY + MAofVE, TWOPI)
C****
C**** Numerically invert Kepler's equation: MA = EA - e*sin(EA)
C****
      EA  = MA + ECCEN*(SIN(MA) + ECCEN*SIN(2*MA)/2)
   10 dEA = (MA - EA + ECCEN*SIN(EA)) / (1 - ECCEN*COS(EA))
      EA  = EA + dEA
      IF(ABS(dEA).gt.1d-10)  GO TO 10
C****
C**** Calculate distance to Sun and true anomaly
C****
      SUNDIS = 1 - ECCEN*COS(EA)
      TA     = ATAN2 (BSEMI*SIN(EA), COS(EA)-ECCEN)
      SDIST  = SUNDIS*SUNDIS   ! added for compatiblity
C****
C**** Change reference frame to be nonrotating reference frame, angles
C**** fixed according to stars, with Earth at center and positive x
C**** axis be ray from Earth to Sun were Earth at vernal equinox, and
C**** x-y plane be Earth's equatorial plane.  Distance from current Sun
C**** to this x axis is SUNDIS sin(TA-TAofVE).  At vernal equinox, Sun
C**** is located at (SUNDIS,0,0).  At other times, Sun is located at:
C****
C**** SUN = (SUNDIS cos(TA-TAofVE),
C****        SUNDIS sin(TA-TAofVE) cos(OBLIQ),
C****        SUNDIS sin(TA-TAofVE) sin(OBLIQ))
C****
      SIND   = SIN(TA-TAofVE) * SIN(OBLIQ)
      COSD   = SQRT (1 - SIND*SIND)
      SUNX   = COS(TA-TAofVE)
      SUNY   = SIN(TA-TAofVE) * COS(OBLIQ)
      SLNORO = ATAN2 (SUNY,SUNX)
C****
C**** Determine Sun location in Earth's rotating reference frame
C**** (normal longitude and latitude)
C****
      VEQLON = TWOPI*VE2000 - PI + MAofVE - TAofVE  !  modulo 2*Pi
      ROTATE = TWOPI*(DAY-VE2000)*(EDAYzY+1)/EDAYzY
      SUNLON = MODULO (SLNORO-ROTATE-VEQLON, TWOPI)
      IF(SUNLON.gt.PI)  SUNLON = SUNLON - TWOPI
      SUNLAT = ASIN (SIN(TA-TAofVE)*SIN(OBLIQ))
C****
C**** Determine longitude of fictitious mean Sun
C**** Calculate Equation of Time
C****
      SLMEAN = PI - TWOPI*(DAY-FLOOR(DAY))
      EQTIME = MODULO (SLMEAN-SUNLON, TWOPI)
      IF(EQTIME.gt.PI)  EQTIME = EQTIME - TWOPI
C****
      RETURN
      END SUBROUTINE ORBIT

#ifdef HEALY_LM_DIAGS
      real*8 function Fe(M,N)
      real*8 M,N

      Fe=0.47d0*log(1.+2.01d-5*(M*N)**(0.75)
     .   +5.31d-15*M*(M*N)**(1.52))
      return
      end
#endif

#ifdef CACHED_SUBDD
      subroutine rijh_defs(arr,nmax,decl_count)
c
c 2D outputs
c
      use subdd_mod, only : info_type,sched_rad
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)
c
c note: next() is a locally declared function to increment decl_count
c

      decl_count = 0

c
      arr(next()) = info_type_(
     &  sname = 'olrrad',
     &  lname = 'OUTGOING LW RADIATION at TOA (in RADIA)',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'olrcs',
     &  lname = 'OUTGOING LW RADIATION at TOA, CLEAR-SKY',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'lwds',
     &  lname = 'LONGWAVE DOWNWARD FLUX at SURFACE',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'lwdscs',
     &  lname = 'LONGWAVE DOWNWARD FLUX at SURFACE, CLEAR-SKY',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'lwus',
     &  lname = 'LONGWAVE UPWARD FLUX at SURFACE',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'totcld',
     &  lname = 'Total Cloud Cover (as seen by rad)',
     &  units = '%',
     &  scale = 1d2,
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'totcld_diag',
     &  lname = 'Total Cloud Cover (continuous, not seen by rad)',
     &  units = '%',
     &  scale = 1d2,
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cldss_2d',
     &  lname = 'Stratiform Cloud Cover',
     &  units = '%',
     &  scale = 1d2,
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cldmc_2d',
     &  lname = 'Convective Cloud Cover',
     &  units = '%',
     &  scale = 1d2,
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cod',
     &  lname = 'Cloud optical depth warm clouds',
     &  units = '-',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cid',
     &  lname = 'Cloud optical depth ice clouds',
     &  units = '-',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'wtrcld',
     &  lname = 'Water cloud frequency',
     &  units = '-',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'icecld',
     &  lname = 'Ice cloud frequency',
     &  units = '-',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'ctt',
     &  lname = 'Cloud top temperature',
     &  units = 'C',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'ctp',
     &  lname = 'Cloud top pressure',
     &  units = 'hPa',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swds',
     &  lname = 'SOLAR DOWNWARD FLUX at SURFACE',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swus',
     &  lname = 'SOLAR UPWARD FLUX at SURFACE',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swdf',
     &  lname = 'SOLAR DOWNWARD DIFFUSE FLUX at SURFACE',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swtoa',
     &  lname = 'SOLAR NET FLUX, TOA',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swns',
     &  lname = 'Solar net flux at surface',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'lwns',
     &  lname = 'Longwave net flux at surface',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
#ifdef CFMIP3_SUBDD  /* CFMIP3_SUBDD */
      arr(next()) = info_type_(
     &  sname = 'rtmt',
     &  lname = 'Net downward radiative flux, TOA',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swut',
     &  lname = 'TOA outgoing SW',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swutcs',
     &  lname = 'TOA outgoing SW, CSKY',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'clwvi',
     &  lname = 'Total water path',
     &  units = 'kg/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swdcls',
     &  lname = 'SFC downward radiative flux, CSKY',
     &  units = 'kg/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swucls',
     &  lname = 'SFC upward radiative flux, CSKY',
     &  units = 'kg/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swus',
     &  lname = 'SFC upward radiative flux',
     &  units = 'kg/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swdt',
     &  lname = 'TOA incoming SW',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
#endif  /* CFMIP3_SUBDD */
      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine rijh_defs

      subroutine rijlh_defs(arr,nmax,decl_count)
c
c 3D outputs
c
      use subdd_mod, only : info_type,sched_rad
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      use constant, only : kapa
      use TimeConstants_mod, only: SECONDS_PER_DAY
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)
c
c note: next() is a locally declared function to increment decl_count
c
      decl_count = 0
c
      arr(next()) = info_type_(
     &  sname = 'dth_sw',
     &  lname = 'theta tendency from shortwave radiative heating',
     &  units = 'K/day',
     &  scale = 1000.**kapa*SECONDS_PER_DAY,
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dth_lw',
     &  lname = 'theta tendency from longwave radiative heating',
     &  units = 'K/day',
     &  scale = 1000.**kapa*SECONDS_PER_DAY,
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dth_rad',
     &  lname = 'theta tendency from radiative heating',
     &  units = 'K/day',
     &  scale = 1000.**kapa*SECONDS_PER_DAY,
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'lwdp',
     &  lname = 'LONGWAVE DOWNWARD FLUX profile',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'lwup',
     &  lname = 'LONGWAVE UPWARD FLUX profile',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swdp',
     &  lname = 'SHORTWAVE DOWNWARD FLUX profile',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swup',
     &  lname = 'SHORTWAVE UPWARD FLUX profile',
     &  units = 'W/m^2',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'MRCO2rad',
     &  lname = 'radiation code CO2 volume mixing ratio',
     &  units = 'mole species / mole air',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'wtrtau',
     &  lname = 'Cloud Water Opacity Seen by Radiation',
     &  units = '-',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'icetau',
     &  lname = 'Cloud Ice Opacity Seen by Radiation',
     &  units = '-',
     &  sched = sched_rad
     &     )
c
#ifdef CFMIP3_SUBDD
      arr(next()) = info_type_(
     &  sname = 'cf',
     &  lname = 'Cloud Fraction',
     &  units = '%',
     &  scale = 1d2,
     &  sched = sched_rad
     &     )
      arr(next()) = info_type_(
     &  sname = 'qcirad',
     &  lname = 'Ice Water Mass Mixing Ratio Seen by Radiation',
     &  units = 'kg/kg',
     &  sched = sched_rad
     &     )
      arr(next()) = info_type_(
     &  sname = 'qclrad',
     &  lname = 'Liquid Water Mass Mixing Ratio Seen by Radiation',
     &  units = 'kg/kg',
     &  sched = sched_rad
     &     )
#endif
c
      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine rijlh_defs

#endif

      subroutine readIFile(IFile)
! Consolidated duplicate of MODELE.f code snippets that read the
! I-file containing the parameter database and INPUTZ namelist.
! Currently used by radiation-only configuration; to be moved
! to MODELE.f and used by all configurations once full testing
! is completed.
! Note that INPUTZ contains fewer variables in this version.
      USE FILEMANAGER, only : openunit,closeunit
      Use Parser_mod
      use Model_com, only : xlabel, lrunid
      USE MODEL_COM, only : HOURI,DATEI,MONTHI,YEARI,IYEAR1
      use diag_com, only : itwrite
      implicit none
C**** Command line options
      character(len=*), intent(in) :: IFile

      integer :: iu_IFILE
!@nlparam IHRI,TIMEE,IHOURE   end of model run
!@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
      CHARACTER NLREC*80,RLABEL*132
      NAMELIST/INPUTZ/
     *     ITWRITE
C****    List of parameters that are disregarded at restarts
     *     ,HOURI,DATEI,MONTHI,YEARI
      NAMELIST/INPUTZ_cold/
     *     ITWRITE
C****    List of parameters that are disregarded at restarts
     *     ,HOURI,DATEI,MONTHI,YEARI
      character*132 :: bufs
      integer, parameter :: MAXLEN_RUNID = 32
      integer :: lid1,lid2,fid,noff


C****
C**** Reading rundeck (I-file) options
C****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      call parse_params(iu_IFILE)
      call closeunit(iu_IFILE)

C****
C**** Print Header and Label (2 lines) from rundeck
C****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      !if (AM_I_ROOT())
      WRITE (6,'(A,40X,A/)') '0','GISS CLIMATE MODEL'
      READ(iu_IFILE,'(A80)') XLABEL(1:80),NLREC
      NOFF=0
      IF (XLABEL(73:80).EQ.'        ') NOFF=8   ! for 72-column rundecks
      XLABEL(81-NOFF:132)=NLREC(1:52+NOFF)
      !if (AM_I_ROOT())
      WRITE (6,'(A,A/)') '0',XLABEL
      RLABEL = XLABEL !@var RLABEL rundeck-label

      lid1 = INDEX(XLABEL,'(') -1
      if (lid1.lt.1) lid1=MAXLEN_RUNID+1
      lid2 = INDEX(XLABEL,' ') -1
      if (lid2.lt.1) lid2=MAXLEN_RUNID+1
      LRUNID = min(lid1,lid2)
      IF (LRUNID.gt.MAXLEN_RUNID) call stop_model
     *     ('INPUT: Rundeck name too long. Shorten to 32 char or less'
     *     ,255)

C****
C**** Read parameters from the rundeck to database and namelist
C****
      do
        read( iu_IFILE, *, err=910, end=910 ) bufs
        if ( bufs == '&&END_PARAMETERS' ) exit
      enddo

      READ (iu_IFILE,NML=INPUTZ,ERR=900)

      call closeunit(iu_IFILE)

      IF (YearI.lt.0) then
        WRITE(6,*) 'Please choose a proper start year yearI, not',yearI
        call stop_model('INPUT: yearI not provided',255)
      END IF

      return
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP
C****
  900 write (6,*) 'Error in NAMELIST parameters'
      call stop_model('Error in NAMELIST parameters',255)
  910 write (6,*) 'Error readin I-file'
      call stop_model('Error reading I-file',255)

      end subroutine readIFile

      subroutine run_radonly(IFile)
!@sum Call single-column radiation-only model once
      USE Dictionary_mod
      USE DOMAIN_DECOMP_1D, ONLY : init_app
      use Model_com, only: itime, itimeE, master_yr, xlabel, lrunid
      USE MODEL_COM, only : YEARI,IYEAR1
      USE DOMAIN_DECOMP_ATM, ONLY : grid,init_grid
#ifdef CACHED_SUBDD
      use diag_com
      use geom, only : lon_dg,lat_dg
      use cdl_mod
#endif
      use TimerPackage_mod, only: initializeTimerPackage_mod=>initialize
      implicit none
C**** Command line options
      character(len=*), intent(in) :: IFile
c
      integer :: i,j,l,n
      character(len=80) :: filenm

      call init_app()

      call initializeTimerPackage_mod() ! avoid probs when RADIA calls timers

      call readIFile(IFile)

      if (is_set_param("master_yr")) then
        call get_param( "master_yr", master_yr )
      else
        call stop_model('Please define master_yr in the rundeck.',255)
      endif

      Iyear1 = yearI

      call sundial

      itimeE = itime+1  ! for length-1 nominal time axis in diags

      call init_grid(grid, 1, 1, 1, width=0)

      !call alloc_clouds_com(grid)
      call alloc_rad_com(grid)
      !call alloc_veg_com(grid)

      call geom_1pt

      CALL init_RAD(2) ! istart=2
      CALL daily_orbit(.false.)             ! not end_of_day
      CALL daily_RAD(.false.)

      call print_param( 6 )

#ifdef CACHED_SUBDD
      ! Initialize diagnostics framework
      call init_cdl_type('cdl_aij',cdl_ij_template)
      call add_coord(cdl_ij_template,'lon',1,units='degrees_east',
     &     coordvalues=lon_dg(:,1))
      call add_coord(cdl_ij_template,'lat',1,units='degrees_north',
     &     coordvalues=lat_dg(:,1))
      call parse_subdd
      call reset_cached_subdd
      call get_subdd_vinterp_coeffs
      call set_subdd_period()
#endif

      call calc_zenith_angle
      call radia

#ifdef CACHED_SUBDD
      filenm = 'allsteps.subdd'//XLABEL(1:LRUNID)
      call write_subdd_accfile (filenm)
#endif

      call stop_model('Radiation calculations completed.',13)

      contains

      subroutine geom_1pt
      use geom
      use constant, only : pi,twopi,radian
      use dictionary_mod, only : get_param,sync_param
      implicit none
      real*8 :: lon_targ,lat_targ

      ! mandatory rundeck parameters: lon and lat of target point
      call get_param('lon_targ',lon_targ)
      call get_param('lat_targ',lat_targ)

      if(abs(lon_targ).gt.180d0 .or. abs(lat_targ).gt.90d0)
     &     call stop_model(
     &       'geom_atm: invalid lon_targ,lat_targ in rundeck',255)

      lon2d_dg(1,1) = lon_targ
      lat2d_dg(1,1) = lat_targ

      axyp(1,1) = 1.

      byaxyp(1,1) = 1d0/axyp(1,1)

      lon2d(1,1) = lon2d_dg(1,1)*radian
      lat2d(1,1) = lat2d_dg(1,1)*radian

      sinlat2d(1,1) = sin(lat2d(1,1))
      coslat2d(1,1) = cos(lat2d(1,1))
      lon2d(1,1) = lon2d(1,1) + pi ! IDL has a value of zero
      if(lon2d(1,1) .lt. 0.) lon2d(1,1)= lon2d(1,1) + twopi

      imaxj = 1

      lon_dg = lon2d_dg
      lat_dg = lat2d_dg

      return
      end subroutine geom_1pt

      subroutine sundial
! Duplicate of relevant snippets of clock initialization in MODELE.f.
! Currently used by radiation-only configuration; will disappear
! once the clock initialization in MODELE.f has been cleanly isolated
! from other intialization activities.
      USE Dictionary_mod
      USE MODEL_COM, only :
     *      nday,dtsrc,itime,itimei
     *     ,HOURI,DATEI,MONTHI,YEARI
      use MODEL_COM, only: modelEclock, calendar
      use ModelClock_mod, only: ModelClock
      use Time_mod
      use BaseTime_mod
      use Rational_mod
      use TimeInterval_mod
      implicit none

      type (Time) :: modelETime0
      type (Time) :: modelETime
      type (TimeInterval) :: dtSrcUsed
      type (TimeInterval) :: secsPerDay

C**** Get those parameters which are needed in this subroutine
      call get_param( "DTsrc", DTsrc )

!@var NDAY=(1 day)/DTsrc : even integer; adjust DTsrc to be commensurate
      secsPerDay = calendar%getSecondsPerDay()
      NDAY = 2*nint((secsPerDay/(DTsrc*2)))
      dtSrcUsed = TimeInterval(secsPerDay / NDAY)
      DTsrc = real(dtSrcUsed)

      modelETime0 = newTime(calendar)
      modelETime = newTime(calendar)

      call modelEtime%setByDate(yearI, monthI, dateI, hourI)
      call modelEtime0%setByDate(yearI, month=1, date=1, hour=0)

      ITimeI = nint((modelEtime - modelEtime0) / dtSrcUsed)
      Itime = ItimeI

      modelEclock = ModelClock(modelEtime,dtSrcUsed,itime)

      CALL DAILY_cal(.false.)   ! not end_of_day

      end subroutine sundial

      end subroutine run_radonly
