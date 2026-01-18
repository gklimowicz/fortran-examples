#include "rundeck_opts.h"

subroutine CONDSE
!@sum   CONDSE driver for moist convection AND large-scale condensation
!@auth  M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@calls CLOUDS:MSTCNV,CLOUDS:LSCOND
  use CONSTANT, only : bygrav,lhm,rgas,grav,tf,lhe,lhs,sha,deltx &
       ,teeny,undef,bysha
  use TimeConstants_mod, only: SECONDS_PER_DAY, INT_HOURS_PER_DAY, &
                               SECONDS_PER_HOUR
  use RESOLUTION, only : im,jm,lm
  use ATM_COM, only : u,v,t,q,qcl,qci
  use DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds,AM_I_ROOT
  use MODEL_COM, only : DTsrc,itime,modelEclock
  use DOMAIN_DECOMP_ATM, only : GLOBALSUM
  use QUSDEF, only : nmom
  use SOMTQ_COM, only : t3mom=>tmom,q3mom=>qmom
  use GEOM, only : imaxj,axyp,byaxyp
#ifndef SCM
  use GEOM, only : kmaxj
#endif
#ifndef CUBED_SPHERE
#ifndef SCM
  use GEOM, only : ravj
#endif
#endif
  use RANDOM
  use RAD_COM, only : cosz1
  use CLOUDS_COM, only : ttold,qtold,svlhx,svlat,rhsav,cldsav &
       ,isccp_reg2d,ukm,vkm,ncol
#ifdef CLD_AER_CDNC
  use CLOUDS_COM, only : ncl,nci,clwp,cdn3d,cre3d  ! for 3 hrly diag
  use DIAG_COM, only : ij_nclssct
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
  use CLOUDS_COM, only :  ctem,cd3d,cl3d,ci3d  ! for 3 hrly diag
#endif
#ifdef TRACERS_AMP
#ifdef BLK_2MOM
  use CLOUDS_COM, only : NACTC
  use DIAG_COM, only : ij_ccnssct
#endif
#endif
  use CLOUDS_COM, only : tauss,taumc,cldss,cldmc,csizmc,csizss,fss,cldsav1 &
       ,tls,qls,tmc,qmc,ddm1,airx,lmc,taussip,csizssip &
       ,QLss_3d=>QLss,QIss_3d=>QIss,QLmc_3d=>QLmc,QImc_3d=>QImc &
       ,ddms,tdn1,qdn1,ddml
#if (defined mjo_subdd) || (defined etc_subdd)
  use CLOUDS_COM, only : CLWC3D,CIWC3D,TLH3D,LLH3D,SLH3D,DLH3D
#endif
#ifdef etc_subdd
  use CLOUDS_COM, only : LWP2D,IWP2D
#endif
#ifdef mjo_subdd
  use CLOUDS_COM, only : TMCDRY,SMCDRY,DMCDRY,LSCDRY
#endif

  use DIAG_COM, only : ftype,ntype,aij=>aij_loc, &
       aijl=>aijl_loc,adiurn=>adiurn_loc,jreg,ij_pscld, &
       ij_pdcld,ij_scnvfrq,ij_dcnvfrq,ij_cnvfrq,ij_wmsum,ij_snwf,ij_prec, &
       ij_neth,ij_f0oc,j_eprcp,j_prcpmc,j_prcpss,ijl_mc, &
       ijdd,idd_pr,idd_ecnd,idd_mcp,idd_dmc,idd_smc,idd_ssp, &
       jl_mcmflx,jl_sshr,jl_mchr,jl_dammc,jl_rhe,jl_mchphas, &
       jl_mcdtotw,jl_mcldht,jl_mcheat,jl_mcdry,ij_ctpi,ij_taui, &
       ij_lcldi,ij_mcldi,ij_hcldi,ij_tcldi,ij_sstabx,isccp_diags, &
       ndiupt,jl_cldmc,jl_cldss,jl_csizmc,jl_csizss,ij_scldi, &
       jl_mcshlw,jl_mcdeep,ij_mccldtp,ij_mccldbs, &
       ij_mccvtp,ij_mccvbs,ij_precoo,ij_precsi,ij_precli,ij_precgr, &
       saveHCLDI,saveMCLDI,saveLCLDI,saveCTPI,saveTAUI,saveSCLDI, &
       saveTCLDI,saveMCCLDTP
#ifdef USE_HDIURN
  use DIAG_COM, only : hdiurn=>hdiurn_loc
#endif
  use DIAG_COM, only : ntau,npres,aisccp=>aisccp_loc,ij_precmc,ij_cldw,ij_cldi &
       ,ij_tclssct,ij_rclssct &
       ,ij_fwoc,p_acc,pm_acc,ndiuvar,nisccp,adiurn_dust,jl_mcdflx &
       ,lh_diags,ijl_llh,ijl_mctlh,ijl_mcdlh,ijl_mcslh &
       ,ijl_ldry,ijl_tmcdry,ijl_dmcdry,ijl_smcdry &
       ,ijl_cldwtr,ijl_cldice,ijl_MCamFX ! ipcc 3-D model layer diagnostics
!!     ,IJ_CONDLS,IJ_EVAPLS,IJ_CONDSINKMC
#ifdef CLD_AER_CDNC
  use DIAG_COM, only : jl_cnumwm,jl_cnumws,jl_cnumim,jl_cnumis &
       ,ij_dzwm,ij_dzim,ij_dzws,ij_dzis &
       ,ij_3dnwm,ij_3dnws,ij_3dnim,ij_3dnis &
       ,ij_3drwm,ij_3drws,ij_3drim,ij_3dris &
       ,ij_3dlwm,ij_3dlws,ij_3dlim,ij_3dlis &
       ,ijl_rewm,ijl_rews,ijl_cdwm,ijl_cdws,ijl_cwwm,ijl_cwws &
       ,ij_wmclwp,ij_wmctwp &
       ,ijl_reim,ijl_reis,ijl_cdim,ijl_cdis,ijl_cwim,ijl_cwis &
       ,ijl_cfwm,ijl_cfim,ijl_cfws,ijl_cfis,ijl_cdtomas
#endif
#ifdef TRACERS_DUST
  use DIAG_COM, only : idd_wet
#endif
#ifdef TRACERS_AMP
#ifdef BLK_2MOM
  use AERO_CONFIG, only: NMODES
  use AMP_AEROSOL, only: NACTV
#endif
#endif
#if (defined CALCULATE_LIGHTNING) || (defined TRACERS_SPECIAL_Shindell)
      USE LIGHTNING,  only : FLASH_DENS, CG_DENS, FLASH_PERTURB,L440mbM1
#ifdef AUTOTUNE_LIGHTNING
      USE LIGHTNING,  only : TUNE_LT_LAND, TUNE_LT_SEA
      USE LIGHTNING,  only : LAND_FR_LIS, SEA_FR_LIS
      USE LIGHTNING,  only : LAND_FR_UNC, SEA_FR_UNC, FLASH_UNC
      USE LIGHTNING,  only : CNT_FR, NHISTLI
#endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE LIGHTNING, only : ENOx_lgt
#endif
#ifdef TRACERS_ON
  use TRACER_COM, only: TRM,TRMOM,NTM,trdn1
  use OldTracer_mod, only: itime_tr0, trname
#ifdef TRACERS_COSMO
  use TRACER_COM, only: n_Be10,n_Be7
#endif
#ifdef TRACERS_WATER
  use OldTracer_mod, only: trw0, dowetdep
  use TRACER_COM, only: trwm
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
  use TRACER_COM, only: Ntm_dust, n_soilDust
#endif
#endif
#ifdef TRACERS_COSMO
  use COSMO_SOURCES, only : BE7W_acc
#endif
  use trdiag_com, only: taijn=>taijn_loc, tij_prec
#ifndef SKIP_TRACER_DIAGS
  use TRDIAG_COM, only: jlnt_mc,jlnt_lscond,itcon_mc,ijlt_prodSO4aq &
       ,itcon_ss,taijs=>taijs_loc,taijls=>taijls_loc
#ifdef TRACERS_WATER
  use TRDIAG_COM, only: jls_prec,trp_acc
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
  use TRDIAG_COM, only: jls_incloud,ijts_aq
#endif
#ifdef TRDIAG_WETDEPO
  use TRDIAG_COM, only: jls_trdpmc,jls_trdpls,ijts_trdpmc,ijts_trdpls
#endif
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
  use TRDIAG_COM, only: jls_wet,ijts_wet,itcon_wt
#endif
#endif
#endif /*SKIP_TRACER_DIAGS*/

  use CLOUDS, only : tm,tmom,trdnl & ! local  (i,j)
       ,ntx,ntix                     ! global (same for all i,j)
#ifdef TRACERS_WATER
  use CLOUDS, only : trwml,trsvwml,trprmc,trprss
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
  use CLOUDS, only : dt_sulf_mc,dt_sulf_ss
#endif
#ifdef TRDIAG_WETDEPO
  use CLOUDS, only : trcond_mc,trdvap_mc,trflcw_mc,trprcp_mc,trnvap_mc,trwash_mc &
       ,trwash_ls,trevap_ls,trclwc_ls,trprcp_ls,trclwe_ls,trcond_ls &
       ,diag_wetdep
#endif
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
  use CLOUDS, only : tm_dust,tmom_dust,trprc_dust
#endif
#endif
#endif

  use CLOUDS, only : BYDTsrc,mstcnv,lscond & ! glb var & subs
       ,airm,byam,etal,sm,smom,qm,qmom,isc,dxypij,LMCLD,hcndss &
       ,tl,ris,ri1,ri2,mcflx,sshr,dgdsm,dphase,dtotw,dqcond,dctei &
!      ,wml,sdl,u_0,v_0,um,vm,um1,vm1,qs,us,vs,dcl,airxl,prcpss &
       ,qcil,qcll,sdl,u_0,v_0,um,vm,um1,vm1,qs,us,vs,airxl,prcpss &
       ,prcpmc,pearth,ts,taumcl,cldmcl,svwmxl,svlatl,svlhxl,dgdqm &
       ,dcl,zpbl,ppbl &
       ,cldslwij,clddepij,csizel,precnvl,vsubl,lmcmax,lmcmin,wmsum &
!      ,aq,dpdt,th,ql,wmx,ttoldl,rh,taussl,cldssl,cldsavl,rh1,roice &
       ,aq,dpdt,th,ql,qcix,qclx,ttoldl,rh,taussl,cldssl,cldsavl,rh1 &
       ,roice &
       ,kmax,ra,pl,ple,plk,rndssl,lhp,debug,fssl,pland,cldsv1 &
       ,smommc,smomls,qmommc,qmomls,ddmflx,wturb &
       ,tvl,w2l,gzl &
       ,dphashlw,dphadeep,dgshlw,dgdeep,tdnl,qdnl,prebar1 &
       ,DQMTOTAL,DQLSC &
       ,DQMSHLW,DQMDEEP,DQCTOTAL,DQCSHLW,DQCDEEP &
       ,QLss,QIss,QLmc,QImc &
       ,use_vmp,wmpr,tausslip,csizelip
#ifdef CLD_AER_CDNC
       use CLOUDS, only : acdnwm,acdnim,acdnws,acdnis,arews,arewm,areis,areim &
       ,alwim,alwis,alwwm,alwws,nlsw,nlsi,nmcw,nmci &
       ,ncll,ncil,sme &
       ,cdn3dl,cre3dl,smlwp &
       ,wmclwp,wmctwp,CDNC_TOMAS

#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
       use CLOUDS, only : cteml,cd3dl,cl3dl,ci3dl
#endif

#ifdef SCM
  use SCM_COM, only : SCMin
  ! plume diagnostics
  use CLOUDS, only : CUMFLX,DWNFLX,WCUALL,ENTALL,DETALL, &
       MPLUMEALL,PLUME_MAX,PLUME_MIN
#endif
  use PBLCOM, only : dclev,egcm,w2gcm,pblht,pblptop
  use ATM_COM, only : pk,pek,pmid,pedn,gz,PMIDOLD,pdsig,MWs, &
       ua=>ualij,va=>valij,ltropo
  use DYNAMICS, only : wcpsig,bydsig
  use SEAICE_COM, only : si_atm
  use GHY_COM, only : fearth
  use RAD_COM, only : snoage
  use LAKES_COM, only : flake
  use FLUXES, only : prec,eprec,precss,focean,fland,flice, &
       atmocn,atmice,atmgla,atmlnd,atmsrf
#ifdef TRACERS_WATER
  use FLUXES, only : trprec
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
  use FLUXES, only : trprec_dust
#endif
#endif

#ifdef TRACERS_AMP
  use AMP_AEROSOL, only : AQsulfRATE
#endif
#ifdef TRACERS_TOMAS
      USE TOMAS_AEROSOL, only : AQSO4oxid_mc,AQSO4oxid_ls
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
  use tracer_sources, only : n__prec
#endif
  use FILEMANAGER, only: openunit,closeunit
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
  use trdust_mod,only : prelay
#endif
  use TimerPackage_mod, only: startTimer => start, stopTimer => stop
#ifdef CACHED_SUBDD
      use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups, &
           inc_subdd,find_groups
#endif
#ifdef SCM
  use SCM_COM, only : SCMopt
#endif
#ifdef CACHED_SUBDD
#ifdef COSP_SIM
  ! Import modules needed for the COSP simulator
  use atm_com, only : zatmo
  use geom, only : lon_dg,lat_dg
  use clouds, only : ccl_cosp,ccp_cosp,reff_cosp,reffp_cosp
  use clouds_com, only : nsubdd_cosp,npoints_cosp,flag_pfluxes_cosp,flag_re_cosp
  use cosp_drv, only : init_cosp_gbx,run_cosp_sims,free_cosp_sims
  use cosp_drv, only : gbx,i_lscliq,i_lscice,i_lsrain,i_lssnow,i_cvcliq, &
    i_cvcice,i_cvrain,i_cvsnow,i_lsgrpl,cosp_bywc,cosp_byic
#endif
#endif
#ifdef CFMIP3_SUBDD
  use CLOUDS, only : wmctwp,wmclwp
  use ATM_COM, only : ma
#endif

  implicit none

#ifdef TRACERS_ON
!@var tmsave holds tracer value (for diagnostics)
  real*8 tmsave(lm,ntm),tmomsv(nmom,lm,ntm),dtrm(lm)
  integer NX
#endif
#if (defined CALCULATE_LIGHTNING) || (defined TRACERS_SPECIAL_Shindell)
!@var Lfreeze Lowest level where temperature is below freezing (TF)
  integer Lfreeze
#endif

  real*8, dimension(LM,GRID%I_STRT_HALO:GRID%I_STOP_HALO, &
       GRID%J_STRT_HALO:GRID%J_STOP_HALO) &
       :: UASV              ! for U tendency diagnostic

#ifdef SCM
  !  moist convection plume diagnostics:
  !  updraft mass flux, downdraft mass flux, updraft speed, entrainment rate,
  !  detrainment rate, plume mass, plume maximum level, plume minimum level
  !@var mc_mfu_p1,mc_mfu_p2,mc_mfd_p1,mc_mfd_p2,mc_w_p1,mc_w_p2,mc_ent_p1,
  !     mc_ent_p2,mc_det_p1,mc_det_p2,mc_m_p1,mc_m_p2,mc_pl_max_p1,mc_pl_max_p2,
  !     mc_pl_min_p1,mc_pl_min_p2
    real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO, &
            GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM,LM) &
            :: mc_mfu_p1,mc_mfu_p2,mc_mfd_p1,mc_mfd_p2,mc_w_p1,mc_w_p2, &
               mc_ent_p1,mc_ent_p2,mc_det_p1,mc_det_p2,mc_m_p1,mc_m_p2
    real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO, &
            GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) &
            :: mc_pl_max_p1,mc_pl_max_p2, mc_pl_min_p1,mc_pl_min_p2
    integer LMIN
#endif

#ifdef CACHED_SUBDD
#ifdef TRACERS_WATER
   !tracer precipitation variable name
   character(len=20) :: trpname
#endif
#ifdef SCM
   !  isccp diagnostics   save frequency histogram for subdd diagnostics
   !@var save_fq_isccp
    real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO, &
           GRID%J_STRT_HALO:GRID%J_STOP_HALO,NTAU,NPRES) &
           :: save_fq_isccp
#endif
#endif

!@param ENTCON fractional rate of entrainment (km**-1)
  real*8,  parameter :: ENTCON = .2d0
  real*8, parameter :: SLHE=LHE*BYSHA

  integer I,J,K,L,N,LL  !@var I,J,K,L,N loop variables
  integer JR,KR,ITYPE,IT,IH,LP850,LP600,IHM,KMAX_NONPOLAR
!@var JR = JREG(I,J)
!@var KR index for regional diagnostics
!@var ITYPE index for snow age
!@var IT index for surface types
!@var LP850 layer near 850 mb
!@var LP600 layer near 600 mb
!@var LERR,IERR error reporting
  integer :: LERR, IERR

  real*8 :: HCNDMC,PRCP,TPRCP,EPRCP,ENRGP,WMERR,ALPHA1,ALPHA2,ALPHAS
  real*8 :: DTDZ,DTDZS,DUDZ,DVDZ,DUDZS,DVDZS,THSV,THV1,THV2,QG,TGV
  real*8 :: DH1S,BYDH1S,DH12,BYDH12,DTDZG,DUDZG,DVDZG,SSTAB,DIFT,CSC &
       ,TSV,WM1,WMI,PWATER
  !ECON*     ,E,E1,W1,ep,ep1,TSV,q0,q1,q2,WM1,WMI
!@var HCNDMC heating due to moist convection
!@var PRCP precipitation
!@var TPRCP temperature of mc. precip  (deg. C)
!@var EPRCP sensible heat of precip
!@var ENRGP total energy of precip
!@var WMERR DH12,BYDH12,DH1S,BYDH1S,SSTAB dummy variable
!@var THSV,THV1,THV2 vertual potential temperatures
!@var QG,TGV ground humidity,virt.temperature from pbl
!@var ALPHA1,ALPHA2,ALPHAS,DIFT,CSC dummy variables
!@var DTDZ,DTDZS,DTDZG vertical potential temperature gradients
!@var DUDZ,DVDZ,DUDZS,DVDZS,DUDZG,DVDZG vertical wind gradients
!@var TSV virtual surface temperature (K)

  !**** parameters and variables for isccp diags
  real*8, parameter :: bywc = 1./2.56d0 , byic= 1./2.13d0
  real*8 skt(1),conv(lm),qv(lm)
  real*8 pfull(lm),at(lm),cc(lm),dtau_s(lm),dtau_c(lm)
  real*8 dem_s(lm),dem_c(lm),phalf(lm+1)
  real*8 fq_isccp(ntau,npres),ctp(1),tauopt(1)
  real*8 boxtau(ncol),boxptop(ncol)

  integer itau,itrop(1),nbox(1),sunlit,ipres
  !****

  integer ICKERR, JCKERR, JERR, seed, NR
  real*8  RNDSS(3,LM,GRID%I_STRT_HALO:GRID%I_STOP_HALO, &
       GRID%J_STRT_HALO:GRID%J_STOP_HALO),xx
  integer :: nij_before_j0,nij_after_j1,nij_after_i1
  real*8  UKMSP(IM,LM), VKMSP(IM,LM), UKMNP(IM,LM),VKMNP(IM,LM)
  real*4 WCU500(IM,16),SAVWCU(IM,16,LM),SAVEN1(IM,16,LM), &
       SAVEN2(IM,16,LM),W500P1(16),ENTJ(16),SAVWC1(IM,16,LM)
  integer :: J_0,J_1,J_0H,J_1H,J_0S,J_1S,I_0,I_1
  logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

#ifdef CACHED_SUBDD
  integer :: igrp,ngroups,grpids(subdd_ngroups)
  type(subdd_type), pointer :: subdd
!@var sddarr temporary array for passing reordered/derived fields
!@+   to inc_subdd
  real*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
                    grid%j_strt_halo:grid%j_stop_halo) :: sddarr
  real*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
                    grid%j_strt_halo:grid%j_stop_halo,lm) :: sddarr3d, &
                    dth_mc,dq_mc,dth_ss,dq_ss
#endif

  integer, parameter :: n_idx1 = 5
  integer, parameter :: n_idx2 = 3
  integer, parameter :: n_idx3 = 6
#ifdef TRACERS_DUST
  integer,parameter :: n_idxd=1
#endif
  integer :: idx1(n_idx1), idx2(n_idx2), idx3(n_idx3)
#ifdef TRACERS_DUST
  integer :: idxd(n_idxd)
#endif
  real*8 :: tmp(NDIUVAR)
#ifndef TRACERS_WATER
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
  integer :: n1,n_fidx
#endif
#endif
#ifdef CACHED_SUBDD
!@var MCPA moist convective precipitation;
!@var LWPA cloud liquid water path
!@var IWPA cloud ice water path
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
              grid%j_strt_halo:grid%j_stop_halo) :: &
              MCPA,LWPA,IWPA
!@var DCNVF_IJ occurence of deep convecvtio; SCNVF_IJ for shallow.
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
              grid%j_strt_halo:grid%j_stop_halo) :: &
              DCNVF_IJ, SCNVF_IJ  ! , SDDARR
!@var Cloud_daily
!     1:  cdnc Large Scale
!     2:  cdnc Large Scale screened after Ralf Bennartz
!     3:  ctp_mc Cloup top preassure convective clouds
!     4:
!     5:  cdnc conv clouds
!     6:
!     7:  lwp
!     8:  convective lwp
!     9:  reff_w_mc
!     10: reff_w_ls
!     11: reff_i_mc
!     12: reff_i_ls
!     13: dzwm
!     14: dzws
!     15: dzim
!     16: dzis
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
              grid%j_strt_halo:grid%j_stop_halo,16) :: &
              Cloud_daily
!@var Cloud_daily3d

      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
              grid%j_strt_halo:grid%j_stop_halo,LM,15) :: &
              Cloud_daily3d
#endif
#ifdef CFMIP3_SUBDD
    !@var cfmip_ctp_mc cloud top pressure convective clouds for SUBDD output 'ctp_mc'
    !@var cfmip_cbp_mc cloud base pressure convective clouds for SUBDD output 'cbp_mc'
    !@var cfmip_cnvfrq fraction time occurs for SUBDD output 'cnvfrq'
    !@var cfmip_dcnvfrq fraction time deep MC occurs for SUBDD output 'dcnvfrq'
    !@var cfmip_scnvfrq fraction time shallow MC occurs for SUBDD output 'scnvfrq'
    !@var cfmip_mc_twp MC total cloud water path for SUBDD output 'mc_twp'
    !@var cfmip_mc_lwp MC liquid cloud water path for SUBDD output 'mc_lwp'
    !@var cfmip_wvp atmospheric water vapour column for SUBDD output 'qatm'
    REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
                      grid%j_strt_halo:grid%j_stop_halo) :: cfmip_ctp_mc, &
                                                            cfmip_cbp_mc, &
                                                            cfmip_cnvfrq, &
                                                            cfmip_dcnvfrq, &
                                                            cfmip_scnvfrq, &
                                                            cfmip_mc_twp, &
                                                            cfmip_mc_lwp, &
                                                            cfmip_wvp
    REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
                      grid%j_strt_halo:grid%j_stop_halo,lm) :: cfmip_mcamfx
#endif
!@var tau_thresh optical depth threshold for cloud-top diagnostics (-)
  real*8, parameter :: tau_thresh = 1d0
!@var tau_sum cumulative layer-average stratiform cloud optical depth from TOA (-)
!@var rcl_sum cumulative layer-average stratiform cloud optical depth times cloud droplet effective radius    (um)
!@var rci_sum cumulative layer-average stratiform cloud optical depth times cloud ice effective radius (um)
!@var ncl_sum cumulative layer-average stratiform cloud optical depth times cloud droplet num conc (m-3)
!@var nci_sum cumulative layer-average stratiform cloud optical depth times cloud ice num conc (m-3)
!@var ccn_sum cumulative layer-average stratiform cloud optical depth times activated CCN conc (m-3)
  real*8 :: tau_sum,rcl_sum,rci_sum,ncl_sum,nci_sum,ccn_sum
!@var dum1,dum2 dummy variables
  real*8 :: dum1,dum2

#ifdef CLD_AER_CDNC
  real*8 :: cldwt,cldwtdz
#endif
#ifdef SCM
  integer :: kmaxj(1)=1
#endif
  integer :: iThread
  integer :: numThreads
  integer :: I_0thread, I_1thread, imaxj_thread
#ifdef CACHED_SUBDD
#ifdef COSP_SIM
  ! Declare local variables for the COSP simulator
!@var np_cosp grid/column index being passed to COSP
  integer :: np_cosp
  real*8 :: scale_pr_cosp
#endif
#endif
  call startTimer('CONDSE()')
  !**** Initialize
#ifdef TRACERS_SPECIAL_Shindell
  ENOx_lgt(:,:)=0.d0
#endif
  idx1 = (/ IDD_PR, IDD_ECND, IDD_MCP, IDD_DMC, IDD_SMC /)
  idx2 = (/ IDD_PR, IDD_ECND, IDD_SSP /)
  idx3 = (/ IDD_PR, IDD_ECND, IDD_MCP, IDD_DMC, IDD_SMC, IDD_SSP /)
#ifdef TRACERS_DUST
  if (adiurn_dust == 1) idxd=(/idd_wet/)
#endif


  !**** define local grid
  call getDomainBounds(grid, J_STRT=J_0,         J_STOP=J_1, &
       J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S, &
       J_STRT_HALO=J_0H,    J_STOP_HALO=J_1H, &
       HAVE_NORTH_POLE=HAVE_NORTH_POLE, &
       HAVE_SOUTH_POLE=HAVE_SOUTH_POLE        )
  I_0 = grid%I_STRT
  I_1 = grid%I_STOP

  !
  !     OBTAIN RANDOM NUMBERS FOR PARALLEL REGION
  !
  !     Burn some random numbers corresponding to latitudes off
  !     processor
  call BURN_RANDOM(nij_before_j0(J_0)*LMCLD*3)

  do J=J_0,J_1
    call BURN_RANDOM((I_0-1)*LMCLD*3)
    do I=I_0,IMAXJ(J)
      do L=LMCLD,1,-1
        do NR=1,3
          RNDSS(NR,L,I,J) = RANDU(xx)
        end do
      end do
      !     Do not bother to save random numbers for isccp_clouds
    end do
    call BURN_RANDOM(nij_after_i1(I_1)*LMCLD*3)
  end do

  call BURN_RANDOM(nij_after_j1(J_1)*LMCLD*3)

  !     But save the current seed in case isccp_routine is activated
  if (isccp_diags.eq.1) call RFINAL(seed)
  WCU500=0.
  SAVWCU=0.
  SAVWC1=0.
  SAVEN1=0.
  SAVEN2=0.
  W500P1=0.
  ENTJ=0.
#ifdef CACHED_SUBDD
#ifdef SCM
  if( SCMopt%PlumeDiag )then
      ! plume diagnostics
      mc_mfu_p1=0.d0; mc_mfu_p2=0.d0; mc_mfd_p1=0.d0; mc_mfd_p2=0.d0
      mc_w_p1=0.d0; mc_w_p2=0.d0; mc_ent_p1=0.d0; mc_ent_p2=0.d0
      mc_m_p1=0.d0; mc_m_p2=0.d0; mc_det_p1=0.d0; mc_det_p2=0.d0
      mc_pl_max_p1=0.d0; mc_pl_max_p2=0.0; mc_pl_min_p1=0.d0; mc_pl_min_p2=0.d0
  endif
      ! isccp frequency diags
      save_fq_isccp=0.d0
#endif
#endif

  call recalc_agrid_uv ! may not be necessary - check later

  !
  ! collect staggered velocities to be mixed into an A-grid array
  !
  kmax_nonpolar = minval(kmaxj(j_0:j_1))
  call replicate_uv_to_agrid(ukm,vkm,kmax_nonpolar, &
       ukmsp,vkmsp,ukmnp,vkmnp)

  !
  !**** SAVE UC AND VC, AND ZERO OUT CLDSS AND CLDMC
  TLS=T
  QLS=Q
  TMC=T
  QMC=Q
  FSS=1.
  IH=modelEclock%getHour()+1
  IHM = IH+(modelEclock%getDate()-1)*INT_HOURS_PER_DAY
#ifdef TRACERS_ON
  !**** Find the ntx active tracers ntix(1->ntx)
  nx = 0
  do n=1,ntm
    if (itime.lt.itime_tr0(n)) cycle
    nx = nx+1
    ntix(nx) = n
  end do
  ntx = nx

#ifdef TRACERS_AMP
  AQsulfRATE = 0.d0
#endif
#ifdef TRACERS_TOMAS
      AQSO4oxid_mc(:,:,:) = 0.d0
      AQSO4oxid_ls(:,:,:) = 0.d0
#endif
#endif
  saveMCCLDTP(:,:)=undef

#ifdef CACHED_SUBDD
#ifdef COSP_SIM
  if (mod(Itime+1,nsubdd_cosp).eq.0) then
    ! Trigger COSP (CFMIP Observation Simulator Package)
    ! Scaling factor
    scale_pr_cosp = 0.5*(100.0*bygrav)/dtsrc
    ! Reset grid counter for the COSP simulator
    np_cosp = 1
    ! Initialize COSP simulator data structures
    call init_cosp_gbx()
  endif ! nsubdd_cosp
#endif
#endif
#ifdef CFMIP3_SUBDD
  cfmip_ctp_mc = 0.d0
  cfmip_cbp_mc = 0.d0
  cfmip_cnvfrq = 0.d0
  cfmip_dcnvfrq = 0.d0
  cfmip_scnvfrq = 0.d0
  cfmip_mc_twp = 0.d0
  cfmip_mc_lwp = 0.d0
  cfmip_wvp = 0.d0
  cfmip_mcamfx = 0.d0
#endif

  numThreads = 1 ! no openmp

#ifdef AUTOTUNE_LIGHTNING
      IF ( SUM(CNT_FR) .gt. 0 ) THEN
         TUNE_LT_LAND = LAND_FR_LIS*SUM(CNT_FR)/SUM(LAND_FR_UNC)
         TUNE_LT_SEA  =  SEA_FR_LIS*SUM(CNT_FR)/SUM(SEA_FR_UNC)
      ENDIF
      ! Shift in anticipation
      LAND_FR_UNC(2:NHISTLI) = LAND_FR_UNC(1:(NHISTLI-1))
      SEA_FR_UNC(2:NHISTLI)  = SEA_FR_UNC(1:(NHISTLI-1))
      CNT_FR(2:NHISTLI)      = CNT_FR(1:(NHISTLI-1))
      LAND_FR_UNC(1) = 0d0
      SEA_FR_UNC(1) = 0d0
      CNT_FR(1) = 0d0
      IF ( Am_I_Root() )  WRITE(6,'(A50,F15.10,F15.10)') &
               'CLOUDS2_DRV: Calculating land/sea tuning parameters:', &
     	       TUNE_LT_LAND, TUNE_LT_SEA
      !CALL STOP_MODEL( 'LTM Testing',17)
      FLASH_DENS = 0d0
      CG_DENS    = 0d0
      FLASH_UNC  = 0d0
#endif

  !****
  !**** MAIN J LOOP
  !****
  ICKERR=0
  JCKERR=0

  ! Burn random numbers for earlier latitudes here.
  ! Actual generation of random numbers is in CLOUDS2.f::ISCCP_CLOUD_TYPES
  if (isccp_diags.eq.1) then
    call BURN_RANDOM(nij_before_j0(J_0)*NCOL*(LM+1))
  end if

  do J=J_0,J_1

    ! Burn random numbers for earlier longitudes here.
    ! Actual generation of random numbers is in CLOUDS2.f::ISCCP_CLOUD_TYPES
    if (isccp_diags.eq.1) then
      call BURN_RANDOM((I_0-1)*NCOL*(LM+1))
    end if

    do ithread = 0, numThreads - 1
      I_0thread = I_0 + (I_1-I_0+1) * iThread / numThreads
      I_1thread = I_0 + (I_1-I_0+1) * (iThread+1) / numThreads  - 1
      Imaxj_thread = min(IMAXJ(J), I_1thread)
      kmax = kmaxj(j)
      !****
      !**** MAIN I LOOP
      !****
      do I=I_0thread,IMAXJ_thread
        DXYPIJ=AXYP(I,J)
        JR=JREG(I,J)
        !****
        !**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS FOR MSTCNV
        !****
        DEBUG = .false.   ! use for individual box diags in clouds
        PEARTH=FEARTH(I,J)
        PLAND=FLAND(I,J)
        PWATER=1.-PLAND
        ROICE=si_atm%RSI(I,J)*PWATER
        TS=atmsrf%TSAVG(I,J)
        QS=atmsrf%QSAVG(I,J)
        US=atmsrf%USAVG(I,J)
        VS=atmsrf%VSAVG(I,J)
        TGV=atmsrf%TGVAVG(I,J)
        QG=atmsrf%QGAVG(I,J)
        TSV=TS*(1+QS*DELTX)
!!!     DCL=NINT(DCLEV(I,J))   ! prevented by openMP bug
        DCL=int(DCLEV(I,J)+.5)
        ZPBL=PBLHT(I,J)
        PPBL=PBLPTOP(I,J)
#ifdef SCM
        if( SCMopt%PlumeDiag )then
        ! plume diagnostics
          do LL=1,LM
            do L=1,LM
              WCUALL(L,1,LL)=0.
              WCUALL(L,2,LL)=0.
              MPLUMEALL(L,1,LL)=0.
              MPLUMEALL(L,2,LL)=0.
              ENTALL(L,1,LL)=0.
              ENTALL(L,2,LL)=0.
              DETALL(L,1,LL)=0.
              DETALL(L,2,LL)=0.
            enddo
          enddo
        endif
#endif
#ifndef SCM
#ifdef CUBED_SPHERE
        ra = .5d0
#else
#ifdef ALT_CLDMIX_UV
        ra(1:kmax) = 1d0
#else
        do K=1,KMAX
          RA(K)=RAVJ(K,J)
        end do
#endif
#endif
#endif
        !**** PRESSURES, AND PRESSURE TO THE KAPA
        PL(:) =PMID(:,I,J)
        PLE(:)=PEDN(:,I,J)
        PLK(:)=PK(:,I,J)
        AIRM(:)=PDSIG(:,I,J)
        BYAM(:)=1./AIRM(:)
        WTURB(:)=sqrt(.6666667*EGCM(:,I,J))
#ifdef CACHED_SUBDD
        Cloud_daily(I,J,:) = 0.d0
        Cloud_daily3d(I,J,:,:) = 0.d0
        Cloud_daily3d(I,J,:,13) = PL(:)
#endif

        !**** other fields where L is the leading index
        SVLHXL(:)=SVLHX(:,I,J)
        TTOLDL(:)=TTOLD(:,I,J)
        CLDSAVL(:)=CLDSAV(:,I,J)
        CLDSV1(:)=CLDSAV1(:,I,J)
        RH(:)=RHSAV(:,I,J)
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
        CTEML(:) =CTEM(:,I,J)
        CD3DL(:) =CD3D(:,I,J)
        CL3DL(:) =CL3D(:,I,J)
        CI3DL(:) =CI3D(:,I,J)
#endif
#ifdef CLD_AER_CDNC
        NCLL(:)=NCL(:,I,J)
        NCIL(:)=NCI(:,I,J)  ! NCI is for rsf save
        SME(:)  =EGCM(:,I,J)  !saving 3D TKE value
        !       if(l.eq.2)write(6,*)"CTEM_DRV",CTEML(L),SME(L),NCLL(L)
        CDN3DL(:)=CDN3D(:,I,J)
        CRE3DL(:)=CRE3D(:,I,J)
        SMLWP=CLWP(I,J)
#ifdef TRACERS_AMP
        !**not sure if this needed
#ifdef BLK_2MOM
        do n=1,nmodes
          NACTC(:,n)= NACTV(I,J,:,n)
        enddo
#endif
#endif
#endif
        FSSL(:)=FSS(:,I,J)
        DPDT(:) = (PMID(:,I,J)-PMIDOLD(:,I,J))*BYDTsrc
        do L=1,LM
          !**** TEMPERATURES
          SM(L)  =T(I,J,L)*AIRM(L)
          SMOM(:,L) =T3MOM(:,I,J,L)*AIRM(L)
          SMOMMC(:,L) =SMOM(:,L)
          SMOMLS(:,L) =SMOM(:,L)
          TL(L)=T(I,J,L)*PLK(L)
          !**** MOISTURE (SPECIFIC HUMIDITY)
          QM(L)  =Q(I,J,L)*AIRM(L)
          QMOM(:,L) =Q3MOM(:,I,J,L)*AIRM(L)
          QMOMMC(:,L) =QMOM(:,L)
          QMOMLS(:,L) =QMOM(:,L)
          QCIL(L)=QCI(I,J,L)
          QCLL(L)=QCL(I,J,L)
          QL(L) =Q(I,J,L)
          !**** others
#ifdef SCM
          ! apply large-scale divergence if specified
          if( SCMopt%omega .or. SCMopt%w )then
            SDL(L) = SCMin%Omega(L)
          else
            SDL(L) = 0.
          endif
#else
          SDL(L) = MWs(I,J,L)/DTsrc*BYAXYP(I,J) ! averaged SD
#endif
          TVL(L)=TL(L)*(1.+DELTX*QL(L))
          W2L(L)=W2GCM(L,I,J)
          if(L.le.LM-2) &
               ETAL(L+1)=.5*ENTCON*(GZ(I,J,L+2)-GZ(I,J,L))*1.d-3*BYGRAV
          if(L.le.LM-2) GZL(L+1)=ETAL(L+1)/ENTCON
        end do


        ETAL(LM)=ETAL(LM-1)
        ETAL(1)=0.     ! not used
        GZL(LM)=GZL(LM-1)
        GZL(1)=0.
#ifdef TRACERS_ON
        !**** TRACERS: Use only the active ones
        do nx=1,ntx
          do l=1,lm
            tm(l,nx) = trm(i,j,l,ntix(nx))
            tmom(:,l,nx) = trmom(:,i,j,l,ntix(nx))
          end do
        end do
#endif

        !**** SURROUNDING WINDS

        if(j.eq.1 .and. have_south_pole) then
          U_0(1:KMAX,:) = UKMSP(1:KMAX,:)
          V_0(1:KMAX,:) = VKMSP(1:KMAX,:)
        elseif(j.eq.jm .and. have_north_pole) then
          U_0(1:KMAX,:) = UKMNP(1:KMAX,:)
          V_0(1:KMAX,:) = VKMNP(1:KMAX,:)
        else
          U_0(1:KMAX,:) = UKM(1:KMAX,:,I,J)
          V_0(1:KMAX,:) = VKM(1:KMAX,:,I,J)
        endif
        do L=1,LM
          do K=1,KMAX
            UM(K,L) = U_0(K,L)*AIRM(L)
            VM(K,L) = V_0(K,L)*AIRM(L)
            UM1(K,L) = UM(K,L)
            VM1(K,L) = VM(K,L)
          end do
        end do

        !**** INITIALISE PRECIPITATION AND LATENT HEAT
        PRCP=0.
        ENRGP=0.
        !**** temperature of precip is based on pre-mstcnv profile
        TPRCP=T(I,J,1)*PLK(1)-TF
#ifdef TRACERS_WATER
        TRPREC(:,I,J) = 0.
#endif

        !**** SET DEFAULTS FOR AIR MASS FLUX (STRAT MODEL)
        AIRX(I,J)=0.
        DDM1(I,J)=0.
        DDMS(I,J)=0.
        DDML(I,J)=0.
        TDN1(I,J)=0.
        QDN1(I,J)=0.
#ifdef TRACERS_ON
        TRDN1(:,I,J)=0.
#endif
        !****
        !**** Energy conservation note: For future reference the energy function
        !**** for these column calculations (assuming energy reference level
        !**** of 0 K for air, and 0 C liquid for water) is:
        !****  E = SH + LH_vapour + LH_clw + ENRGP
        !****    =  (sum(TL(:)*AIRM(:))*SHA + sum(QM(:))*LHE +
        !****        sum(WML(:)*(LHE-SVLHXL(:))*AIRM(:)))*100.*BYGRAV
        !**** The LH_clw term is slightly different after MSTCNV:
        !****   LH_clw = sum((WML(:)*(LHE-SVLHXL(:))+SVWMXL(:)*(LHE-SVLATL(:)))
        !****                *AIRM(:))*100.*BYGRAV
        !**** After LSCOND, latent heat changes to:
        !****          = sum(WMX(:)*(LHE-SVLHXL(:))*AIRM(:))*100.*BYGRAV
        !****
        !**** Note that the column changes after MSTCNV only apply to the
        !**** moist convective fraction (1-FSSL(:)), and after LSCOND, FSSL(:).
        !**** Condensate is always defined over the whole box.
        !****
        !**** uncomment lines marked ECON to check energy conservation
        !**** uncomment lines marked QCON to check water conservation

        !QCON q0 = sum(QM(:)+WML(:)*AIRM(:))*100.*BYGRAV
        !ECON  E = (sum(TL(:)*AIRM(:))*SHA + sum(QM(:))*LHE +sum(WML(:)*(LHE
        !ECON*     -SVLHXL(:))*AIRM(:)))*100.*BYGRAV

#ifdef CACHED_SUBDD
        !**** save initial values to calculate rates of change below
        do L=1,LM
          dth_mc(I,J,L) = T(I,J,L)
          dq_mc(I,J,L) = Q(I,J,L)
        enddo
#endif

        !**** MOIST CONVECTION
        call MSTCNV(IERR,LERR,i,j)

        !ECON E1 = ( sum( ((T(I,J,:)*PLK(:)-TL(:))*AIRM(:)*SHA + (Q(I,J,:)
        !ECON*     *AIRM(:)-QM(:))*LHE)*(1.-FSSL(:)))-sum(SVWMXL(:)*(LHE
        !ECON*     -SVLATL(:))*AIRM(:)))*100.*BYGRAV
        !QCON q1 = sum( (Q(I,J,:)*AIRM(:)-QM(:))*(1.-FSSL(:))-SVWMXL(:)*AIRM(:))
        !QCON*     *100.*BYGRAV

        !**** Error reports
        if (ierr.gt.0) then
          write(6,*) "Error in moist conv: i,j,l=",i,j,lerr
          if (ierr.eq.2) ickerr = ickerr + 1
        end if

#if (defined CALCULATE_LIGHTNING) || (defined TRACERS_SPECIAL_Shindell)
      FLASH_DENS(i,j) = 0.d0 ! default for subdaily diag
      CG_DENS(i,j)    = 0.d0 ! default for subdaily diag
#ifdef AUTOTUNE_LIGHTNING
      FLASH_UNC(i,j)  = 0.d0 ! default
#endif
      ! Execute lightning parameterization:
      ! first, need the local freezing level:
      if(LMCMAX>0)then
        Lfreeze=1
        do L=1,LMCMAX
          if(T(i,j,L)*plk(L)<TF) then
            Lfreeze=L
            exit
          endif
        enddo
        call calc_lightning(i,j,LMCMAX,Lfreeze,MCFLX(L440mbM1),PRCPMC)
      endif
#endif

        !**** ACCUMULATE MOIST CONVECTION DIAGNOSTICS
        if (LMCMIN.gt.0) then
          AIJ(I,J,IJ_PSCLD)=AIJ(I,J,IJ_PSCLD)+CLDSLWIJ
          AIJ(I,J,IJ_PDCLD)=AIJ(I,J,IJ_PDCLD)+CLDDEPIJ
          if(CLDSLWIJ.gt.1e-6) AIJ(I,J,IJ_SCNVFRQ)=AIJ(I,J,IJ_SCNVFRQ)+1.
          if(CLDDEPIJ.gt.1e-6) AIJ(I,J,IJ_DCNVFRQ)=AIJ(I,J,IJ_DCNVFRQ)+1.
          aij(i,j,ij_cnvfrq)=aij(i,j,ij_cnvfrq)+1.
          AIJ(I,J,IJ_WMSUM)=AIJ(I,J,IJ_WMSUM)+WMSUM
          AIJ(I,J,IJ_MCCLDTP)=AIJ(I,J,IJ_MCCLDTP)+   & ! MC cloud top pressure
               PLE(LMCMAX+1)*CLDMCL(LMCMAX)
          AIJ(I,J,IJ_MCCLDBS)=AIJ(I,J,IJ_MCCLDBS)+   & ! MC cloud base pressure
               PLE(LMCMIN+1)*CLDMCL(LMCMIN+1)
          AIJ(I,J,IJ_MCCVTP)=AIJ(I,J,IJ_MCCVTP)+     & ! MC top cloud cover
               CLDMCL(LMCMAX)
          AIJ(I,J,IJ_MCCVBS)=AIJ(I,J,IJ_MCCVBS)+     & ! MC base cloud cover
               CLDMCL(LMCMIN+1)
#ifdef CLD_AER_CDNC
          AIJ(I,J,IJ_WMCLWP)=AIJ(I,J,IJ_WMCLWP)+WMCLWP
          AIJ(I,J,IJ_WMCTWP)=AIJ(I,J,IJ_WMCTWP)+WMCTWP
#ifdef CACHED_SUBDD
          Cloud_daily(I,J,7) = Cloud_daily(I,J,7)+ WMSUM
          Cloud_daily(I,J,8) = WMCLWP
          Cloud_daily(I,J,3) =  PLE(LMCMAX+1)*CLDMCL(LMCMAX)
#endif
#endif
#ifdef CFMIP3_SUBDD  /* CFMIP3_SUBDD */
          ! MC cloud top pressure
          cfmip_ctp_mc(i,j) = ple(lmcmax+1)
          ! MC cloud base pressure
          cfmip_cbp_mc(i,j) = ple(lmcmin+1)
          ! Fraction of Time MC Occurs
          cfmip_cnvfrq(i,j)=cfmip_cnvfrq(i,j)+1.
          ! Fraction of Time deep MC Occurs
          if (clddepij.gt.1e-6) cfmip_dcnvfrq(i,j)=cfmip_dcnvfrq(i,j)+1.
          ! Fraction of Time Shallow MC Occurs
          if(cldslwij.gt.1e-6) cfmip_scnvfrq(i,j)=cfmip_scnvfrq(i,j)+1.
          ! MC total water path
          cfmip_mc_twp(i,j) = max(0.d0, wmctwp)
          ! MC liquid water path
          cfmip_mc_lwp(i,j) = max(0.d0, wmclwp)
          ! ! Column Water vapor
          ! do l=1,lm
          !   cfmip_wvp(i,j) = cfmip_wvp(i,j)+q(i,j,l)*ma(l,i,j)
          ! enddo
          do l=1,lmcmax
            ! MC air mass flux, limit 100 [hPa/m2]
            if (mcflx(l).gt.0d0 .and. mcflx(l).lt.100d0) then
              cfmip_mcamfx(i,j,l) = max(0.d0, mcflx(l))
            endif
          enddo
#endif  /* CFMIP3_SUBDD */
          ! Also save instantaneous MC cloud top pressure for SUBDDiags:
          saveMCCLDTP(i,j)=PLE(LMCMAX+1)

          HCNDMC=0.
          do L=1,LMCMAX
            HCNDMC=HCNDMC+DGDSM(L)+DPHASE(L)
            call inc_ajl(i,j,l,jl_mchr,DGDSM(L)*BYDSIG(L))
            call inc_ajl(i,j,l,jl_mchphas,DPHASE(L)*BYDSIG(L))
            call inc_ajl(i,j,l,jl_mcdtotw,DTOTW(L)*BYDSIG(L))
            !CC       IF(J.GE.J5S.AND.J.LE.J5N) AIL(I,L,IL_MCEQ)=AIL(I,L,IL_MCEQ)+
            !CC  *         (DGDSM(L)+DPHASE(L))*(AXYP(I,J)*BYDSIG(L))
            AIJL(I,J,L,IJL_MC) = AIJL(I,J,L,IJL_MC) + (DPHASE(L)+DGDSM(L))
            call inc_ajl(i,j,l,jl_mcheat,(DPHASE(L)+DGDSM(L)))
            call inc_ajl(i,j,l,jl_mcdry,(DQCOND(L)-DGDQM(L)))
            call inc_ajl(i,j,l,jl_mcshlw,(DPHASHLW(L)+DGSHLW(L)))
            call inc_ajl(i,j,l,jl_mcdeep,(DPHADEEP(L)+DGDEEP(L)))
            !*** Begin Accumulate 3D convective latent heating
            !*** for monthly diags
            if(lh_diags.eq.1) then
              AIJL(I,J,L,IJL_MCTLH)=AIJL(I,J,L,IJL_MCTLH)+ &
                   (DPHASE(L)+DGDSM(L))
              AIJL(I,J,L,IJL_MCDLH)=AIJL(I,J,L,IJL_MCDLH)+ &
                   (DPHADEEP(L)+DGDEEP(L))
              AIJL(I,J,L,IJL_MCSLH)=AIJL(I,J,L,IJL_MCSLH)+ &
                   (DPHASHLW(L)+DGSHLW(L))
              AIJL(I,J,L,IJL_TMCDRY)=AIJL(I,J,L,IJL_TMCDRY)+ &
                   (DQCTOTAL(L)-DQMTOTAL(L))
              AIJL(I,J,L,IJL_DMCDRY)=AIJL(I,J,L,IJL_DMCDRY)+ &
                   (DQCDEEP(L)-DQMDEEP(L))
              AIJL(I,J,L,IJL_SMCDRY)=AIJL(I,J,L,IJL_SMCDRY)+ &
                   (DQCSHLW(L)-DQMSHLW(L))
            endif
#if (defined mjo_subdd) || (defined etc_subdd)
            !*** For subdaily diags
            TLH3D(L,I,J)=TLH3D(L,I,J)+(DPHASE(L)+DGDSM(L))*BYAM(L)
            SLH3D(L,I,J)=SLH3D(L,I,J)+(DPHASHLW(L)+DGSHLW(L))*BYAM(L)
            DLH3D(L,I,J)=DLH3D(L,I,J)+(DPHADEEP(L)+DGDEEP(L))*BYAM(L)
#endif
#ifdef mjo_subdd
            !*** For subdaily diags
            TMCDRY(L,I,J)=TMCDRY(L,I,J)+(DQCTOTAL(L)-DQMTOTAL(L))
            SMCDRY(L,I,J)=SMCDRY(L,I,J)+(DQCSHLW(L)-DQMSHLW(L))
            DMCDRY(L,I,J)=DMCDRY(L,I,J)+(DQCDEEP(L)-DQMDEEP(L))
#endif
            !*** End Accumulate 3D convective latent heating
            call inc_ajl(i,j,l,jl_mcmflx,MCFLX(L))
            call inc_ajl(i,j,l,jl_cldmc,CLDMCL(L)*AIRM(L))
            call inc_ajl(i,j,l,jl_mcdflx,DDMFLX(L))
            call inc_ajl(i,j,l,jl_csizmc,CSIZEL(L)*CLDMCL(L)*AIRM(L))
            aijl(i,j,l,ijl_MCamFX) = aijl(i,j,l,ijl_MCamFX) + MCFLX(L)
          end do
          do IT=1,NTYPE
            call INC_AJ(I,J,IT,J_PRCPMC,PRCPMC*FTYPE(IT,I,J))
          end do
          call INC_AREG(I,J,JR,J_PRCPMC,PRCPMC)
          do KR=1,NDIUPT
            if(I.eq.IJDD(1,KR).and.J.eq.IJDD(2,KR)) then
              tmp(IDD_PR)  =+PRCPMC
              tmp(IDD_ECND)=+HCNDMC
              tmp(IDD_MCP) =+PRCPMC
              tmp(IDD_DMC) =+CLDDEPIJ
              tmp(IDD_SMC) =+CLDSLWIJ
              ADIURN(IDX1(:),KR,IH)=ADIURN(IDX1(:),KR,IH)+TMP(IDX1(:))
#ifdef USE_HDIURN
              HDIURN(IDX1(:),KR,IHM)=HDIURN(IDX1(:),KR,IHM)+TMP(IDX1(:))
#endif
            end if
          end do
#ifdef CLD_AER_CDNC
          do L =1,LM
            if(SVWMXL(L).le.0.) cycle
            CLDWT = CLDMCL(L)!+teeny
            CLDWTDZ = CLDWT*(RGAS*TL(L)*BYGRAV)*(AIRM(L)/PL(L))
            if(SVLATL(L).eq.LHE) then
              AIJL(I,J,L,IJL_CFWM)= AIJL(I,J,L,IJL_CFWM)+CLDWT
              AIJL(I,J,L,IJL_REWM)= AIJL(I,J,L,IJL_REWM)+AREWM(L)*CLDWT
              AIJL(I,J,L,IJL_CDWM)= AIJL(I,J,L,IJL_CDWM)+ACDNWM(L)*CLDWT
              AIJL(I,J,L,IJL_CWWM)= AIJL(I,J,L,IJL_CWWM)+ALWWM(L)*CLDWT
              AIJ(I,J,IJ_DZWM)=AIJ(I,J,IJ_DZWM)+CLDWTDZ
              AIJ(I,J,IJ_3dNWM)=AIJ(I,J,IJ_3dNWM)+ACDNWM(L)*CLDWTDZ
              AIJ(I,J,IJ_3dRWM)=AIJ(I,J,IJ_3dRWM)+AREWM(L)*CLDWTDZ
              AIJ(I,J,IJ_3dLWM)=AIJ(I,J,IJ_3dLWM)+ALWWM(L)*CLDWTDZ
#ifdef CACHED_SUBDD
      Cloud_daily(I,J,5) =  Cloud_daily(I,J,5) + ACDNWM(L)*CLDWTDZ
      Cloud_daily(I,J,9) = Cloud_daily(I,J,9)+AREWM(L)*CLDWTDZ
      Cloud_daily(I,J,13) = Cloud_daily(I,J,13)+CLDWTDZ
      Cloud_daily3d(I,J,L,2) = ACDNWM(L)*CLDWT
      Cloud_daily3d(I,J,L,11) = CLDWT
      Cloud_daily3d(I,J,L,4) = ALWWM(L)*CLDWT
      Cloud_daily3d(I,J,L,7) = AREWM(L)*CLDWT
#endif
            elseif(SVLATL(L).eq.LHS) then
              AIJL(I,J,L,IJL_CFIM)= AIJL(I,J,L,IJL_CFIM)+CLDWT
              AIJL(I,J,L,IJL_REIM)= AIJL(I,J,L,IJL_REIM)+AREIM(L)*CLDWT
              AIJL(I,J,L,IJL_CDIM)= AIJL(I,J,L,IJL_CDIM)+ACDNIM(L)*CLDWT
              AIJL(I,J,L,IJL_CWIM)= AIJL(I,J,L,IJL_CWIM)+ALWIM(L)*CLDWT
              AIJ(I,J,IJ_DZIM)=AIJ(I,J,IJ_DZIM)+CLDWTDZ
              AIJ(I,J,IJ_3dNIM)=AIJ(I,J,IJ_3dNIM)+ACDNIM(L)*CLDWTDZ
              AIJ(I,J,IJ_3dRIM)=AIJ(I,J,IJ_3dRIM)+AREIM(L)*CLDWTDZ
              AIJ(I,J,IJ_3dLIM)=AIJ(I,J,IJ_3dLIM)+ALWIM(L)*CLDWTDZ
#ifdef CACHED_SUBDD
      Cloud_daily(I,J,11) =  Cloud_daily(I,J,11) + AREIM(L)*CLDWTDZ
      Cloud_daily(I,J,15) =  Cloud_daily(I,J,15)+CLDWTDZ
      Cloud_daily3d(I,J,l,5) = Cloud_daily3d(I,J,l,5)+ALWIM(L)*CLDWT
      Cloud_daily3d(I,J,l,9) = AREIM(L)*CLDWT
#endif
            endif
            if (NMCW.ge.1) then
              call inc_ajl(i,j,l,JL_CNUMWM,ACDNWM(L)*AIRM(L))
              !         write(6,*)"IJL_REWM",AIJL(I,J,L,IJL_REWM),I,J,L,
              !     *   AIJL(I,J,L,IJL_CDWM),AIJL(I,J,L,IJL_CWWM),ALWWM(L)
            endif
            if (NMCI.ge.1) then
              call inc_ajl(i,j,l,JL_CNUMIM,ACDNIM(L)*AIRM(L))
              !         write(6,*)"IJL_REIM",AIJL(I,J,L,IJL_REIM),I,J,L,
              !     *   AIJL(I,J,L,IJL_CDIM),AIJL(I,J,L,IJL_CWIM),ALWIM(L)
            endif

          enddo
#endif
          !**** ACCUMULATE PRECIP
          PRCP=PRCPMC*100.*BYGRAV
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
          precnvl(1)=precnvl(1)+prcpmc*bygrav
#endif
          !**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
          !**** NEED TO TAKE ACCOUNT OF LATENT HEAT THOUGH
          if (TPRCP.gt.0) then
            !         EPRCP=PRCP*TPRCP*SHW
            EPRCP=0.
            ENRGP=ENRGP+EPRCP
            !ECON     ep=0.
          else
            !         EPRCP=PRCP*TPRCP*SHI
            EPRCP=0.
            ENRGP=ENRGP+EPRCP-PRCP*LHM
            !ECON     ep=-PRCP*LHM
            AIJ(I,J,IJ_SNWF)=AIJ(I,J,IJ_SNWF)+PRCP
            AIJ(I,J,atmice%IJ_SISNWF)=AIJ(I,J,atmice%IJ_SISNWF) + &
                 PRCP*FOCEAN(I,J)*si_atm%RSI(I,J)
          end if
          AIJ(I,J,IJ_PRECMC)=AIJ(I,J,IJ_PRECMC)+PRCP

          !**** Uncomment next lines for check on conservation
          !ECON   if (abs(E1-ep).gt.0.01) print*,"energy err0",i,j,(E1-ep)
          !ECON*       *GRAV/100.,E,E1,ep,prcp,tprcp
          !QCON   if (abs(q1-prcp).gt.0.01) print*,"water err0",i,j,(q1-prcp)
          !QCON*       *GRAV/100.,q0,q1,prcp

          do L=1,LMCMAX
            T(I,J,L)=(1.-FSSL(L))*SM(L)*BYAM(L)+FSSL(L)*TLS(I,J,L)
            Q(I,J,L)=(1.-FSSL(L))*QM(L)*BYAM(L)+FSSL(L)*QLS(I,J,L)
            TMC(I,J,L)=SM(L)*BYAM(L)
            QMC(I,J,L)=QM(L)*BYAM(L)
            SMOMMC(:,L)=SMOM(:,L)
            QMOMMC(:,L)=QMOM(:,L)
            do K=1,KMAX
              UM1(K,L)=UM(K,L)
              VM1(K,L)=VM(K,L)
            end do
          end do

          CSIZMC(1:LMCMAX,I,J)=CSIZEL(1:LMCMAX)
          FSS(:,I,J)=FSSL(:)
          AIRX(I,J) = AIRXL*AXYP(I,J)
          do L=1,DCL
            DDML(I,J)=L                    ! the lowest downdraft layer
            if(DDMFLX(L).gt.0.d0) exit
          end do
          if (DDML(I,J).gt.0) then
            TDN1(I,J)=TDNL(DDML(I,J))     ! downdraft temperature
            QDN1(I,J)=QDNL(DDML(I,J))     ! downdraft humidity
            DDMS(I,J)=-100.*DDMFLX(DDML(I,J))/(GRAV*DTsrc) ! downdraft mass flux
#ifdef TRACERS_ON
            do nx=1,ntx
               TRDN1(ntix(nx),I,J)=1d-2*TRDNL(nx,DDML(I,J))*GRAV*BYAXYP(I,J) ! downdraft tracer conc
            end do
#endif
          end if

          !**** level 1 downfdraft mass flux/rho (m/s)
          DDM1(I,J) = DDMFLX(1)*RGAS*TSV/(GRAV*PEDN(1,I,J)*DTSrc)
        end if
#ifdef CFMIP3_SUBDD
          ! Column Water vapor
          do l=1,lm
            cfmip_wvp(i,j) = cfmip_wvp(i,j)+q(i,j,l)*ma(l,i,j)
          enddo
#endif
#ifdef CACHED_SUBDD
        ! calculate profile changes and save initial values again
        do L=1,LM
          dth_mc(I,J,L) = T(I,J,L)-dth_mc(I,J,L)
          dq_mc(I,J,L) = Q(I,J,L)-dq_mc(I,J,L)
          dth_ss(I,J,L) = T(I,J,L)
          dq_ss(I,J,L) = Q(I,J,L)
        enddo
#ifdef SCM
        if( SCMopt%PlumeDiag )then
        ! plume diagnostics
          mc_mfu_p1(I,J,:,:) = CUMFLX(:,1,:)
          mc_mfu_p2(I,J,:,:) = CUMFLX(:,2,:)
          mc_mfd_p1(I,J,:,:) = DWNFLX(:,1,:)
          mc_mfd_p2(I,J,:,:) = DWNFLX(:,2,:)
          mc_w_p1(I,J,:,:) = WCUALL(:,1,:)
          mc_w_p2(I,J,:,:) = WCUALL(:,2,:)
          mc_ent_p1(I,J,:,:) = ENTALL(:,1,:)
          mc_ent_p2(I,J,:,:) = ENTALL(:,2,:)
          mc_det_p1(I,J,:,:) = DETALL(:,1,:)
          mc_det_p2(I,J,:,:) = DETALL(:,2,:)
          mc_m_p1(I,J,:,:) = MPLUMEALL(:,1,:)
          mc_m_p2(I,J,:,:) = MPLUMEALL(:,2,:)
          mc_pl_max_p1(I,J,:) = PLUME_MAX(1,:)
          mc_pl_max_p2(I,J,:) = PLUME_MAX(2,:)
          mc_pl_min_p1(I,J,:) = PLUME_MIN(1,:)
          mc_pl_min_p2(I,J,:) = PLUME_MIN(2,:)
        endif
#endif
#endif

#ifdef TRACERS_ON
        !**** TRACERS: Use only the active ones
        do nx=1,ntx
          n = ntix(nx)

#ifndef SKIP_TRACER_DIAGS
          if(lmcmax > 0) then
            do l=1,lmcmax
              dtrm(l) = (tm(l,nx)-trm(i,j,l,n))*(1.-fssl(l))
#ifdef TRACERS_WATER
              dtrm(l) = dtrm(l) + trsvwml(nx,l)
#else
#endif
            enddo
            if(itcon_mc(n).gt.0) call inc_diagtcb(i,j,sum(dtrm(1:lmcmax)), &
                 itcon_mc(n),n)
            call inc_tajln_column(i,j,1,lmcmax,lm,jlnt_mc,n,dtrm)
          endif
#endif  /*SKIP_TRACER_DIAGS*/

          do l=1,lm

#ifdef TRACERS_WATER
            trwml(nx,l) = trwm(i,j,l,n)+trsvwml(nx,l)
#endif
            tmsave(l,nx) = tm(l,nx) ! save for tajln(large-scale condense)
            tmomsv(:,l,nx) = tmom(:,l,nx)
            tm(l,nx) = trm(i,j,l,n)*fssl(l)   ! kg in lsc fraction only
            tmom(:,l,nx) = trmom(:,i,j,l,n)*fssl(l)
          end do
#ifdef TRACERS_WATER
          trprec(n,i,j) = trprmc(nx)
#endif
        end do
#endif
        LMC(1,I,J) = LMCMIN
        LMC(2,I,J) = LMCMAX+1
        !****
        !**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS FOR LSCOND
        !****
        do L=1,LM
          TL(L)=TLS(I,J,L)*PLK(L)
          TH(L)=TLS(I,J,L)
          QL(L)=QLS(I,J,L)
          SMOM(:,L)=SMOMLS(:,L)
          QMOM(:,L)=QMOMLS(:,L)
        end do
!       WMX(:)=WML(:)+SVWMXL(:)
        QCLX(:)=QCLL(:)
        QCIX(:)=QCIL(:)
        DO L=1,LMCMAX
         IF(SVLATL(L).EQ.LHE) THEN
          QCLX(L)=QCLX(L)+SVWMXL(L)
         ELSEIF(SVLATL(L).EQ.LHS) THEN
          QCIX(L)=QCIX(L)+SVWMXL(L)
         ENDIF
        END DO
        AQ(:)=(QL(:)-QTOLD(:,I,J))*BYDTsrc
        RNDSSL(:,1:LMCLD)=RNDSS(:,1:LMCLD,I,J)
        FSSL(:)=FSS(:,I,J)
        do L=1,LM
          do K=1,KMAX
            UM(K,L) = U_0(K,L)*AIRM(L)
            VM(K,L) = V_0(K,L)*AIRM(L)
          end do
        end do
        !****
        !**** COMPUTE STRATOCUMULUS CLOUDS USING PHILANDER'S FORMULA
        !****
        if (ISC.eq.1.and.FOCEAN(I,J).gt..5) then
          CSC=0.D0
          LP600=LM
          LP850=LM
          do L=2,LM
            if(L.gt.LP600) exit
            if(PL(L).lt.600.) then
              LP600=L
              if(600.-PL(L).gt.PL(L-1)-600.) LP600=L-1
            endif
          enddo
          do L=2,LM
            if(L.gt.LP850) exit
            if(PL(L).lt.850.) then
              LP850=L
              if(850.-PL(L).gt.PL(L-1)-850.) LP850=L-1
            endif
          enddo
          if(SDL(LP600)+SDL(LP600+1).gt.0.) then
            DIFT=TL(LP850)-TGV/(1.+DELTX*QG)
            CSC=.031D0*DIFT+.623D0
            if(CSC.lt.0.) CSC=0.
          endif
          CLDMCL(1)=CLDMCL(1)+CSC
          if(CSC.gt.0.) TAUMCL(1)=AIRM(1)*.08D0
          if(CLDMCL(1).gt.1.) CLDMCL(1)=1.
          !     IF(CSC.GT.0.) WRITE (6,*) I,J,DCL,TL(LP850),TGV/(1.+DELTX*QG),CSC
        endif

        !**** COMPUTE RICHARDSON NUMBER FROM SURFACE CONDITIONS WHEN DEPTH OF
        !**** BOUNDARY LAYER IS AT OR BELOW FIRST LAYER (E.G. AT NIGHT)
        if(DCL.le.1) then
          THSV=TS*(1.+DELTX*QS)/PEK(1,I,J)
          THV1=TH(1)*(1.+DELTX*QL(1))
          THV2=TH(2)*(1.+DELTX*QL(2))
          ALPHAS=2./((TGV/(1.+DELTX*QG)+TS)/PEK(1,I,J))
          ALPHA1=2./(TH(1)+TS/PEK(1,I,J))
          ALPHA2=2./(TH(1)+TH(2))
          DH1S=(PLE(1)-PL(1))*TL(1)*RGAS/(GRAV*PL(1))
          BYDH1S=1./DH1S
          DH12=(GZ(I,J,2)-GZ(I,J,1))*BYGRAV
          BYDH12=1./DH12
          DTDZS=(THV1-THSV)*BYDH1S
          DTDZ=(THV2-THV1)*BYDH12
          DUDZ=(UA(2,I,J)-UA(1,I,J))*BYDH12
          DVDZ=(VA(2,I,J)-VA(1,I,J))*BYDH12
          DUDZS=(UA(1,I,J)-US)*BYDH1S
          DVDZS=(VA(1,I,J)-VS)*BYDH1S
          DUDZG=.1d0*US
          DVDZG=.1d0*VS
          DTDZG=.1d0*(THSV-TGV/PEK(1,I,J))
          RIS=(GRAV*ALPHAS*DTDZG)/(DUDZG*DUDZG+DVDZG*DVDZG+teeny)
          RI1=(GRAV*ALPHA1*DTDZS)/(DUDZS*DUDZS+DVDZS*DVDZS+teeny)
          RI2=(GRAV*ALPHA2*DTDZ)/(DUDZ*DUDZ+DVDZ*DVDZ+teeny)
          !       WRITE (6,*)'I,J,QG,TGV,THSV,RIS,RI1=',I,J,QG,TGV,THSV,RIS,RI1
        endif

        !**** Uncomment next few lines for check on conservation
        !ECON W1 = sum( (WML(1:LMCLD)*(LHE-SVLHXL(1:LMCLD))
        !ECON*     +SVWMXL(1:LMCLD)*(LHE-SVLATL(1:LMCLD)))*AIRM(1:LMCLD))
        !ECON E = ( sum(TL(1:LMCLD)*AIRM(1:LMCLD))*SHA + sum(QL(1:LMCLD)
        !ECON*     *AIRM(1:LMCLD))*LHE )*100.*BYGRAV

        !**** LARGE-SCALE CLOUDS AND PRECIPITATION
        call LSCOND(IERR,WMERR,LERR,i,j)

        !ECON E1 = ( sum( ((TLS(I,J,1:LMCLD)-TH(1:LMCLD))*PLK(1:LMCLD)*AIRM(1:LMCLD)
        !ECON*     *SHA +(QLS(I,J,1:LMCLD)-QL(1:LMCLD))*AIRM(1:LMCLD)*LHE)
        !ECON*     *FSSL(1:LMCLD))+W1-sum(WMX(1:LMCLD)*(LHE-SVLHXL(1:LMCLD))
        !ECON*     *AIRM(1:LMCLD)) )*100.*BYGRAV

        !**** Error reports
        if (IERR.ne.0) write(99,'(I10,3I4,A,D14.5,A)') &
             Itime,I,J,LERR,' CONDSE:H2O<0',WMERR,' ->0'

        !**** Accumulate diagnostics of LSCOND

#ifdef CLD_AER_CDNC
        ! code transplanted from LSCOND
        SMLWP=WMSUM
        CDN3DL=0.
        CRE3DL=0.
        ACDNWS=0.
        ACDNIS=0.
        AREWS=0.
        AREIS=0.
        ALWWS=0.
        ALWIS=0.
        NLSW = 0
        NLSI = 0
!    CDNC_TOMAS=0.
        do l=1,lmcld
          if(CLDSV1(L).gt.1.d-5) then
            if(SVLHXL(L).eq.LHE) then
              ! max(...,20) b/c NCLL had different lower limit than SCDNCW
              ACDNWS(L)= max(NCLL(L),20d0)
              AREWS(L) = CSIZEL(L)
              ALWWS(L) = 1.d5*QCLX(L)*PL(L)/(CLDSV1(L)*TL(L)*RGAS+teeny)
              CDN3DL(L) = NCLL(L)
              CRE3DL(L)=CSIZEL(L)
              NLSW  = NLSW + 1
        !      if(ACDNWS(L).gt.20.d0) write(6,*)"INWCLD",ACDNWS(L),
        !    * SCDNCW,NLSW,AREWS(L),RCLDE,LHX
            elseif(SVLHXL(L).eq.LHS) then
              ACDNIS(L)= NCIL(L)
              AREIS(L) = CSIZEL(L)
              ALWIS(L) = 1.d5*QCIX(L)*PL(L)/(CLDSV1(L)*TL(L)*RGAS+teeny)
              CDN3DL(L) = NCIL(L)
              CRE3DL(L)=CSIZEL(L)
              NLSI  = NLSI + 1
        !      if(ACDNIS(L).gt.0.d0)    write(6,*)"INICLD",ACDNIS(L),
        !    * SCDNCI,NLSI,AREIS(L),RCLDE,LHX
            endif
          endif
        enddo
#endif

        AIJ(I,J,IJ_WMSUM)=AIJ(I,J,IJ_WMSUM)+WMSUM
#ifdef CACHED_SUBDD
      Cloud_daily(I,J,7) = Cloud_daily(I,J,7) + WMSUM
#endif
        do IT=1,NTYPE
          call INC_AJ(I,J,IT,J_PRCPSS,PRCPSS*FTYPE(IT,I,J))
        end do
        call INC_AREG(I,J,JR,J_PRCPSS,PRCPSS)
        do KR=1,NDIUPT
          if(I.eq.IJDD(1,KR).and.J.eq.IJDD(2,KR)) then
            tmp(IDD_PR)  =+PRCPSS
            tmp(IDD_ECND)=+HCNDSS
            tmp(IDD_SSP) =+PRCPSS
            ADIURN(IDX2(:),KR,IH)=ADIURN(IDX2(:),KR,IH)+TMP(IDX2(:))
#ifdef USE_HDIURN
            HDIURN(IDX2(:),KR,IHM)=HDIURN(IDX2(:),KR,IHM)+TMP(IDX2(:))
#endif
          end if
        end do

        !**** TOTAL PRECIPITATION AND AGE OF SNOW
        PRCP=PRCP+PRCPSS*100.*BYGRAV

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
        do l=1,lm
          prelay(i,j,l)=((prebar1(l)*DTsrc*100.+precnvl(l)*100.)+ &
               (prebar1(l+1)*DTsrc*100.+precnvl(l+1)*100.))/2.
        end do
#endif

        !**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
        !**** NEED TO TAKE ACCOUNT OF LATENT HEAT THOUGH
        if (LHP(1).ne.LHS) then
          !       EPRCP=PRCPSS*100.*BYGRAV*TPRCP*SHW
          EPRCP=0.
          ENRGP=ENRGP+EPRCP
          !ECON   ep1=0.
        else
          !       EPRCP=PRCPSS*100.*BYGRAV*TPRCP*SHI
          EPRCP=0.
          ENRGP=ENRGP+EPRCP-PRCPSS*100.*BYGRAV*LHM
          !ECON   ep1=-PRCPSS*100.*BYGRAV*LHM
          AIJ(I,J,IJ_SNWF)=AIJ(I,J,IJ_SNWF)+PRCPSS*100.*BYGRAV
          AIJ(I,J,atmice%IJ_SISNWF)=AIJ(I,J,atmice%IJ_SISNWF) + &
               PRCPSS*FOCEAN(I,J)*si_atm%RSI(I,J)
        end if

        !ECON if (abs(E1-ep1).gt.0.01) print*,"energy err1",i,j,(E1-ep1)
        !ECON*     *GRAV*1d-2,E1-ep1,E1,ep1,prcpss*100.*BYGRAV,lhp(1)

        !**** PRECIPITATION DIAGNOSTICS
        do IT=1,NTYPE
          call INC_AJ(I,J,IT,J_EPRCP,ENRGP*FTYPE(IT,I,J))
        end do
        call INC_AREG(I,J,JR,J_EPRCP,ENRGP)
        AIJ(I,J,IJ_PREC)=AIJ(I,J,IJ_PREC)+PRCP
        AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+ENRGP
        AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+ &
             ENRGP*FOCEAN(I,J)*(1.-si_atm%RSI(I,J))
        AIJ(I,J,IJ_FWOC)=AIJ(I,J,IJ_FWOC)+ &
             PRCP*FOCEAN(I,J)*(1.-si_atm%RSI(I,J))
        AIJ(I,J,IJ_PRECOO)=AIJ(I,J,IJ_PRECOO)+PRCP*PWATER*(1.-si_atm%RSI(I,J))
        AIJ(I,J,IJ_PRECSI)=AIJ(I,J,IJ_PRECSI)+PRCP*PWATER*si_atm%RSI(I,J)
        AIJ(I,J,IJ_PRECLI)=AIJ(I,J,IJ_PRECLI)+PRCP*FLICE(I,J)
        AIJ(I,J,IJ_PRECGR)=AIJ(I,J,IJ_PRECGR)+PRCP*FEARTH(I,J)

        if(ENRGP.lt.0.) then    ! MODIFY SNOW AGES AFTER SNOW FALL
          do ITYPE=1,3
            SNOAGE(ITYPE,I,J)=SNOAGE(ITYPE,I,J)*exp(-PRCP)
          end do
        end if

!**** LS condensation & evaporation diagnostics
!!    AIJ(I,J,IJ_CONDLS)=AIJ(I,J,IJ_CONDLS) &
!!           +SUM(CONDLS_SAVE(:))*100.*BYGRAV
!!    AIJ(I,J,IJ_EVAPLS)=AIJ(I,J,IJ_EVAPLS) &
!!           +SUM(EVAPLS_SAVE(:))*100.*BYGRAV

        !**** cloud water diagnostics
        WM1=0  ; WMI=0
        do L=1,LMCLD
          if(SVLHXL(L).eq.LHE) then
            aijl(i,j,l,ijl_cldwtr) = aijl(i,j,l,ijl_cldwtr) + QCLX(L)*AIRM(L)
#ifdef mjo_subdd
            CLWC3D(L,I,J)=CLWC3D(L,I,J)+QCLX(L)
#endif
#ifdef etc_subdd
            CLWC3D(L,I,J)=QCLX(L)
            LWP2D(I,J)=LWP2D(I,J)+QCLX(L)*AIRM(L)*100.*BYGRAV
#endif
          endif
          if(SVLHXL(L).eq.LHS) then
            aijl(i,j,l,ijl_cldice) = aijl(i,j,l,ijl_cldice) + QCIX(L)*AIRM(L)
#ifdef mjo_subdd
            CIWC3D(L,I,J)=CIWC3D(L,I,J)+QCIX(L)  ! WMX(L)
#endif
#ifdef etc_subdd
            CIWC3D(L,I,J)=QCIX(L)
            IWP2D(I,J)=IWP2D(I,J)+QCIX(L)*AIRM(L)*100.*BYGRAV
#endif
          endif
          WM1=WM1+(QCLX(L)+QCIX(L))*AIRM(L)
          if (SVLHXL(L).eq.LHS) WMI=WMI+QCIX(L)*AIRM(L)
          IF (USE_VMP .AND. LHP(L).eq.LHS) then
           ! Count ice precip generated this timestep as part of IWP.
           ! Todo: diagnose IWP both with and without ice precip and report
           ! the one appropriate for the context.
           ! Todo 2: for timestep independence, introduce a time constant.
            WMI=WMI+WMPR(L)*AIRM(L)
            aijl(i,j,l,ijl_cldice)=aijl(i,j,l,ijl_cldice)+WMPR(L)*AIRM(L)
          ENDIF
        end do
        AIJ(I,J,IJ_CLDW)=AIJ(I,J,IJ_CLDW)+WM1*100.*BYGRAV   ! all condensate
        AIJ(I,J,IJ_CLDI)=AIJ(I,J,IJ_CLDI)+WMI*100.*BYGRAV   ! ice only
        !**** Calculate ISCCP cloud diagnostics if required
        if (isccp_diags.eq.1) then
          fq_isccp = 0.d0
          ctp(1) = 0.d0
          tauopt(1) = 0.d0
          boxtau = 0.d0
          boxptop = 0.d0
          saveCTPI(i,j)=0.d0
          saveTAUI(i,j)=0.d0
          saveLCLDI(i,j)=0.d0
          saveMCLDI(i,j)=0.d0
          saveHCLDI(i,j)=0.d0

          do l=1,lm
            cc(l)=cldmcl(LM+1-L)+cldssl(LM+1-L)
            if(cc(l) .gt. 1.) then
              cc(l)=1.
            endif
            conv(l)=cldmcl(LM+1-L)
            if(conv(l) .gt. 1.) then
              conv(l)=1.
            endif

            dtau_s(l)=taussl(LM+1-L)
            dtau_c(l)=taumcl(LM+1-L)
            pfull(l)=pl(LM+1-L)*100.
            phalf(l)=ple(LM+2-L)*100.
            at(l)=tl(LM+1-L)  ! in situ temperature

            !**** set skt from radiative temperature
            skt=sqrt(sqrt( &
                 (focean(i,j)+flake(i,j))*(1.-si_atm%rsi(i,j))*atmocn%gtempr(i,j)**4+ &
                 (focean(i,j)+flake(i,j))*    si_atm%rsi(i,j) *atmice%gtempr(i,j)**4+ &
                 flice(i,j) *atmgla%gtempr(i,j)**4+ &
                 fearth(i,j)*atmlnd%gtempr(i,j)**4))
            dem_s(l)=0.
            dem_c(l)=0.
            if(svlhxl(LM+1-L) .eq. lhe )   & ! large-scale water cloud
                 dem_s(l)=1.-exp(-taussl(LM+1-L)*bywc)
            if(svlatl(LM+1-L) .eq. lhe )   & ! convective water cloud
                 dem_c(l)=1.-exp(-taumcl(LM+1-L)*bywc)
            if(svlhxl(LM+1-L) .eq. lhs )   & ! large-scale ice cloud
                 dem_s(l)=1.-exp(-taussl(LM+1-L)*byic)
            if(svlatl(LM+1-L) .eq. lhs )   & ! convective ice cloud
                 dem_c(l)=1.-exp(-taumcl(LM+1-L)*byic)

            qv(l)=ql(LM+1-L)
          end do
          phalf(lm+1)=ple(1)*100.
          itrop = LM+1-LTROPO(I,J)
          sunlit=0
          if (cosz1(i,j).gt.0) sunlit=1

          call ISCCP_CLOUD_TYPES([sunlit],pfull,phalf,qv, &
               cc,conv,dtau_s,dtau_c,skt, &
               at,dem_s,dem_c,itrop,fq_isccp,ctp,tauopt, &
               boxtau,boxptop,nbox,jerr)
          if(jerr.ne.0) jckerr = jckerr + 1

          !**** set ISCCP diagnostics
          AIJ(I,J,IJ_SCLDI) = AIJ(I,J,IJ_SCLDI) + sunlit
          saveSCLDI(i,j)=sunlit
          if (nbox(1).gt.0.and.sunlit.gt.0) then
            AIJ(I,J,IJ_CTPI) = AIJ(I,J,IJ_CTPI) + ctp(1)
            AIJ(I,J,IJ_TAUI) = AIJ(I,J,IJ_TAUI) + tauopt(1)
            AIJ(I,J,IJ_TCLDI)= AIJ(I,J,IJ_TCLDI)+ 1.
            saveCTPI(i,j)=ctp(1)    ! saving just the
            saveTAUI(i,j)=tauopt(1) ! current value for
            saveTCLDI(i,j)=1.d0  ! instantaneous SUBDDiags
            !**** note LOW CLOUDS:       ipres=6,7
            !****      MID-LEVEL CLOUDS: ipres=4,5
            !****      HIGH CLOUDS:      ipres=1,2,3
            !**** Sum over itau=2,ntau (itau=1 is no cloud)
            AIJ(I,J,IJ_LCLDI)=AIJ(I,J,IJ_LCLDI)+sum(fq_isccp(2:ntau,6:7))
            AIJ(I,J,IJ_MCLDI)=AIJ(I,J,IJ_MCLDI)+sum(fq_isccp(2:ntau,4:5))
            AIJ(I,J,IJ_HCLDI)=AIJ(I,J,IJ_HCLDI)+sum(fq_isccp(2:ntau,1:3))
            saveLCLDI(i,j)=sum(fq_isccp(2:ntau,6:7)) ! saving just the
            saveMCLDI(i,j)=sum(fq_isccp(2:ntau,4:5)) ! current value for
            saveHCLDI(i,j)=sum(fq_isccp(2:ntau,1:3)) ! instant. SUBDDiags
#ifdef CACHED_SUBDD
#ifdef SCM
            save_fq_isccp(i,j,:,:) = fq_isccp(:,:)
#endif
#endif
            !**** Save area weighted isccp histograms
            n=isccp_reg2d(i,j)
            if (n.gt.0) AISCCP(:,:,n) = AISCCP(:,:,n) + &
                 fq_isccp(:,:)*axyp(i,j)
          end if
        end if
#ifdef CACHED_SUBDD
#ifdef COSP_SIM /* COSP_SIM */
        !@auth Mike Bauer
        !
        ! Populate COSP data structures (Only for nsubdd_cosp)
        ! See init_cosp_gbx() in COSP_drv.F90 for detailed description of
        ! the gbx variables.
        if (mod(Itime+1,nsubdd_cosp).eq.0) then
          ! Geo-location Information for each grid being processed by COSP
          ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          gbx%longitude(np_cosp) = modulo(lon_dg(i,1)+360.0,360.0)
          gbx%latitude(np_cosp) = lat_dg(j,1)
          ! Surface Information for each grid being processed by COSP
          ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          gbx%land(np_cosp)=merge(1,0,pland.ge.0.5)
          gbx%skt(np_cosp)=sqrt(sqrt(                                            &
            (focean(i,j)+flake(i,j))*(1.-si_atm%rsi(i,j))*atmocn%gtempr(i,j)**4+ &
            (focean(i,j)+flake(i,j))*si_atm%rsi(i,j) *atmice%gtempr(i,j)**4+     &
            flice(i,j)*atmgla%gtempr(i,j)**4+fearth(i,j)*atmlnd%gtempr(i,j)**4))
          gbx%psfc(np_cosp)=ple(1)*100.0
          gbx%sunlit(np_cosp)=merge(1,0,cosz1(i,j).gt.0)
          ! Column Based Information for each grid being processed by COSP, SFC -> TOA
          ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          gbx%p(np_cosp,:) = pl(:)*100.0
          gbx%ph(np_cosp,:) = ple(1:lm)*100.
          gbx%t(np_cosp,:) = tl(:)
          gbx%q(np_cosp,:) = rh(:)*100.
          gbx%sh(np_cosp,:) = ql(:)
          gbx%zlev(np_cosp,:) = bygrav*gz(i,j,:)
          ! Some quantities need to be dealt with by layer
          do l=1,lm
            if (l.eq.1) then
              gbx%zlev_half(np_cosp,l) = 0.5*bygrav*(gz(i,j,l)+zatmo(i,j))
            else
              gbx%zlev_half(np_cosp,l) = 0.5*bygrav*(gz(i,j,l)+gz(i,j,l-1))
            endif
            gbx%tca(np_cosp,l) = min(1.d0,cldmcl(l)+cldssl(l))
            if (cldssl(l).le.teeny) then
              ! No LS cloud, ensure consistency
              gbx%dtau_s(np_cosp,l) = 0.0
              gbx%dem_s(np_cosp,l) = 0.0
              gbx%mr_hydro(np_cosp,l,i_lscliq) = 0.0
              gbx%reff(np_cosp,l,i_lscliq) = 0.0
              gbx%mr_hydro(np_cosp,l,i_lscice) = 0.0
              gbx%reff(np_cosp,l,i_lscice) = 0.0
            else
              gbx%dtau_s(np_cosp,l) = taussl(l)
              gbx%dem_s(np_cosp,l) = 0.0

              ! Trigger by latent heat test
              if (svlhxl(l).eq.lhe) then ! LS Cloud Water
                if (qclx(l).gt.teeny) then
                    ! Non-zero LS cloud water mass
                    gbx%mr_hydro(np_cosp,l,i_lscliq) = qclx(l)
                    gbx%dem_s(np_cosp,l) = 1.-exp(-taussl(l)*cosp_bywc)
                    if (flag_re_cosp) then
                      if (reff_cosp(l).le.teeny) then
                        ! Zero EffRad, ensure consistency
                        gbx%mr_hydro(np_cosp,l,i_lscliq) = 0.0
                        gbx%reff(np_cosp,l,i_lscliq) = 0.0
                        gbx%dem_s(np_cosp,l) = 0.0
                        gbx%dtau_s(np_cosp,l) = 0.0
                      else
                        gbx%reff(np_cosp,l,i_lscliq) = 1.e-6*reff_cosp(l)
                      endif ! teeny
                    endif ! flag_re_cosp
                else
                  ! Zero LS cloud water mass, ensure consistency
                  gbx%mr_hydro(np_cosp,l,i_lscliq) = 0.0
                  gbx%dem_s(np_cosp,l) = 0.0
                  gbx%dtau_s(np_cosp,l) = 0.0
                  if (flag_re_cosp) gbx%reff(np_cosp,l,i_lscliq) = 0.0
                endif ! qclx
              endif ! LS Cloud Water

              ! Trigger by latent heat test
              if (svlhxl(l).eq.lhs) then ! LS Cloud Ice
                if (qcix(l).gt.teeny) then
                  ! Non-zero LS cloud ice mass
                  gbx%mr_hydro(np_cosp,l,i_lscice) = qcix(l)
                  gbx%dem_s(np_cosp,l) = 1.-exp(-taussl(l)*cosp_byic)
                  if (flag_re_cosp) then
                    if (reff_cosp(l).le.teeny) then
                      ! Zero EffRad, ensure consistency
                      gbx%mr_hydro(np_cosp,l,i_lscice) = 0.0
                      gbx%reff(np_cosp,l,i_lscice) = 0.0
                      gbx%dem_s(np_cosp,l) = 0.0
                      gbx%dtau_s(np_cosp,l) = 0.0
                    else
                      gbx%reff(np_cosp,l,i_lscice) = 1.e-6*reff_cosp(l)
                    endif ! teeny
                  endif ! flag_re_cosp
                else
                  ! Zero LS cloud ice mass, ensure consistency
                  gbx%mr_hydro(np_cosp,l,i_lscice) = 0.0
                  gbx%dem_s(np_cosp,l) = 0.0
                  gbx%dtau_s(np_cosp,l) = 0.0
                  if (flag_re_cosp) gbx%reff(np_cosp,l,i_lscice) = 0.0
                endif ! qcix
              endif ! LS Cloud Ice
            endif ! cldssl

            ! Question: No need to include flag_pfluxes_cosp as never happens (as in Sandbox)?
            if (lhp(l).eq.lhe) then ! LS Rain
              if (flag_pfluxes_cosp) then
                gbx%mr_hydro(np_cosp,l,i_lsrain) = 50.*(prebar1(l)+prebar1(l+1))
              else
                if (0.5*(prebar1(l)+prebar1(l+1)).ge.teeny) then
                  ! Because if prebar1(l) and prebar1(l+1) zero, the following
                  ! will still be non-zero
                  gbx%mr_hydro(np_cosp,l,i_lsrain) = 50.*dtsrc*grav*           &
                    ((prebar1(l)*byam(l))+(prebar1(l+1)*byam(l+1)))
                endif
              endif ! flag_pfluxes_cosp
              if (gbx%mr_hydro(np_cosp,l,i_lsrain).gt.teeny) then
                ! Non-zero LS rain mass
                if (flag_re_cosp) then
                  if (reffp_cosp(l).le.teeny) then
                    ! Zero EffRad, ensure consistency
                    gbx%reff(np_cosp,l,i_lsrain) = 0.0
                    gbx%mr_hydro(np_cosp,l,i_lsrain) = 0.0
                  else
                   gbx%reff(np_cosp,l,i_lsrain) = 1.e-6*reffp_cosp(l)
                  endif ! teeny
                endif ! flag_re_cosp
              else
                ! Zero LS rain mass, ensure consistency
                gbx%mr_hydro(np_cosp,l,i_lsrain) = 0.0
                if (flag_re_cosp) gbx%reff(np_cosp,l,i_lsrain) = 0.0
              endif ! mr_hydro
            endif  ! LS Rain

            ! Question: Can there be LS Rain and Snow at the same time? if/else if
            if (lhp(l).eq.lhs) then ! LS Snow
              if (flag_pfluxes_cosp) then
                gbx%mr_hydro(np_cosp,l,i_lssnow) = 50.*(prebar1(l)+prebar1(l+1))
              else
                if (0.5*(prebar1(l)+prebar1(l+1)).ge.teeny) then
                  ! Because if prebar1(l) and prebar1(l+1) zero, the following
                  ! will still be non-zero
                  gbx%mr_hydro(np_cosp,l,i_lssnow) = 50.*dtsrc*grav*            &
                    ((prebar1(l)*byam(l))+(prebar1(l+1)*byam(l+1)))
                endif
              endif  ! flag_pfluxes_cosp
              if (gbx%mr_hydro(np_cosp,l,i_lssnow).gt.teeny) then
                ! Non-zero LS snow mass
                if (flag_re_cosp) then
                  if (reffp_cosp(l).le.teeny) then
                      ! Zero EffRad, ensure consistency
                      gbx%reff(np_cosp,l,i_lssnow) = 0.0
                      gbx%mr_hydro(np_cosp,l,i_lssnow) = 0.0
                  else
                    gbx%reff(np_cosp,l,i_lssnow) = 1.e-6*reffp_cosp(l)
                  endif ! teeny
                endif ! flag_re_cosp
              else
                ! Zero LS snow mass, ensure consistency
                gbx%mr_hydro(np_cosp,l,i_lssnow) = 0.0
                if (flag_re_cosp) gbx%reff(np_cosp,l,i_lssnow) = 0.0
              endif ! mr_hydro
            endif ! LS Snow

            ! LS Graupel, N/A modelE
            gbx%mr_hydro(np_cosp,l,i_lsgrpl) = 0.0
            gbx%reff(np_cosp,l,i_lsgrpl) = 0.0

            ! MC Cloud and Precip
            if (cldmcl(l).le.teeny) then
              ! No MC cloud ensure consistency
              gbx%cca(np_cosp,l) = 0.0
              gbx%dtau_c(np_cosp,l) = 0.0
              gbx%dem_c(np_cosp,l) = 0.0
              gbx%mr_hydro(np_cosp,l,i_cvcliq) = 0.0
              gbx%reff(np_cosp,l,i_cvcliq) = 0.0
              gbx%mr_hydro(np_cosp,l,i_cvcice) = 0.0
              gbx%reff(np_cosp,l,i_cvcliq) = 0.0
            else
              gbx%cca(np_cosp,l) = min(1.d0, cldmcl(l))
              gbx%dtau_c(np_cosp,l) = taumcl(l)
              gbx%dem_c(np_cosp,l) = 0.0
              if (svlatl(l).eq.lhe) then ! MC Cloud Water
                if (ccl_cosp(l).gt.teeny) then
                  ! Non-zero MC cloud water mass
                  gbx%mr_hydro(np_cosp,l,i_cvcliq) = ccl_cosp(l)
                  gbx%dem_c(np_cosp,l) = 1.-exp(-taumcl(l)*cosp_bywc)
                  if (flag_re_cosp) then
                    if (csizmc(l,i,j).le.teeny) then
                      ! Zero EffRad, ensure consistency
                      gbx%reff(np_cosp,l,i_cvcliq) = 0.0
                      gbx%mr_hydro(np_cosp,l,i_cvcliq) = 0.0
                      gbx%dem_c(np_cosp,l) = 0.0
                      gbx%dtau_c(np_cosp,l) = 0.0
                    else
                      ! Question: is csizmc MC Cloud Water only?
                      gbx%reff(np_cosp,l,i_cvcliq) = 1.0e-6*csizmc(l,i,j)
                    endif ! teeny
                  endif ! flag_re_cosp
                else
                  ! Zero MC cloud liq mass, ensure consistency
                  gbx%mr_hydro(np_cosp,l,i_cvcliq) = 0.0
                  gbx%dem_c(np_cosp,l) = 0.0
                  gbx%dtau_c(np_cosp,l) = 0.0
                  if (flag_re_cosp) gbx%reff(np_cosp,l,i_cvcliq) = 0.0
                endif ! ccl_cosp
              endif ! MC Cloud Water

              ! Question: Can svlatl test if if/else if?
              if (svlatl(l).eq.lhs) then ! MC Cloud Ice
                if (ccl_cosp(l).gt.teeny) then
                  ! Non-zero MC cloud ice mass
                  gbx%mr_hydro(np_cosp,l,i_cvcice) = ccl_cosp(l)
                  gbx%dem_c(np_cosp,l) = 1.-exp(-taumcl(l)*cosp_byic)
                  if (flag_re_cosp) then
                    if (csizss(l,i,j).le.teeny) then
                      ! Zero MC cloud ice effective radius, ensure consistency
                      gbx%reff(np_cosp,l,i_cvcice) = 0.0
                      gbx%mr_hydro(np_cosp,l,i_cvcice) = 0.0
                      gbx%dem_c(np_cosp,l) = 0.0
                      gbx%dtau_c(np_cosp,l) = 0.0
                    else
                      ! Question: is csizss MC Cloud Ice only?
                      gbx%reff(np_cosp,l,i_cvcice) = 1.0e-6*csizss(l,i,j)
                    endif ! teeny
                  endif ! flag_re_cosp
                else
                  ! Zero MC cloud ice mass, ensure consistency
                  gbx%mr_hydro(np_cosp,l,i_cvcice) = 0.0
                  gbx%dem_c(np_cosp,l) = 0.0
                  gbx%dtau_c(np_cosp,l) = 0.0
                  if (flag_re_cosp) gbx%reff(np_cosp,l,i_cvcice) = 0.0
                endif ! ccl_cosp
              endif ! MC Cloud Ice
            endif ! cldmcl

            ! Question: No need to include flag_pfluxes_cosp as never happens (as in Sandbox)?
            if (lhp(l).eq.lhe) then ! MC Rain
              if (flag_pfluxes_cosp) then
                gbx%mr_hydro(np_cosp,l,i_cvrain) = scale_pr_cosp*(ccp_cosp(l)  &
                  +ccp_cosp(l+1))
              else
                gbx%mr_hydro(np_cosp,l,i_cvrain) = 0.5*(ccp_cosp(l)+ccp_cosp(l+1))
              endif ! flag_pfluxes_cosp
              if (gbx%mr_hydro(np_cosp,l,i_cvrain).gt.teeny) then
                ! Non-zero MC cloud rain mass
                if (flag_re_cosp) then
                  if (csizmc(l,i,j).le.teeny) then
                    ! Zero MC rain effective radius
                    gbx%reff(np_cosp,l,i_cvrain) = 0.0
                    gbx%mr_hydro(np_cosp,l,i_cvrain) = 0.0
                  else
                    ! Rain only?
                    gbx%reff(np_cosp,l,i_cvrain) = 1.0e-6*csizmc(l,i,j)
                  endif ! teeny
                endif ! flag_re_cosp
              else
                ! Zero MC rain mass
                gbx%mr_hydro(np_cosp,l,i_cvrain) = 0.0
                if (flag_re_cosp) gbx%reff(np_cosp,l,i_cvrain) = 0.0
              endif ! mr_hydro
            endif ! MC Rain

            ! Question: Can there be MC Rain and Snow at the same time? if/else if
            if (lhp(l).eq.lhs) then ! MC Snow
              if (flag_pfluxes_cosp) then
                gbx%mr_hydro(np_cosp,l,i_cvsnow) = scale_pr_cosp*              &
                  (ccp_cosp(l)+ccp_cosp(l+1))
              else
                gbx%mr_hydro(np_cosp,l,i_cvsnow) = 0.5*(ccp_cosp(l)+ccp_cosp(l+1))
              endif ! flag_pfluxes_cosp
              if (gbx%mr_hydro(np_cosp,l,i_cvsnow).gt.teeny) then
                ! Non-zero MC cloud snow mass
                if (flag_re_cosp) then
                  if (csizss(l,i,j).le.teeny) then
                    ! Zero MC snow effective radius, ensure consistency
                    gbx%reff(np_cosp,l,i_cvsnow) = 0.0
                    gbx%mr_hydro(np_cosp,l,i_cvsnow) = 0.0
                  else
                    ! Question: is csizss MC Cloud Ice only?
                   gbx%reff(np_cosp,l,i_cvsnow) = 1.0e-6*csizss(l,i,j)
                  endif ! teeny
                endif ! flag_re_cosp
              else
                ! Zero MC snow mass, ensure consistency
                gbx%mr_hydro(np_cosp,l,i_cvsnow) = 0.0
                if (flag_re_cosp) gbx%reff(np_cosp,l,i_cvsnow) = 0.0
              endif ! mr_hydro
            endif ! MC Snow
          end do ! lm
          ! Increment COSP grid counter
          np_cosp = np_cosp + 1
        endif ! nsubdd_cosp
        ! Done Populating COSP data structures
#endif /* COSP_SIM */
#endif
        !**** Peak static stability diagnostic
        SSTAB=-1.d30
        do L=1,DCL
          if(SSTAB.lt.(TH(L+1)-TH(L))/(GZ(I,J,L+1)-GZ(I,J,L))) &
             SSTAB =  (TH(L+1)-TH(L))/(GZ(I,J,L+1)-GZ(I,J,L))
        end do
        AIJ(I,J,ij_sstabx) = AIJ(I,J,ij_sstabx) + SSTAB

        !**** WRITE TO GLOBAL ARRAYS
        TAUMC(:,I,J)=TAUMCL(:)
        CLDMC(:,I,J)=CLDMCL(:)
        SVLAT(:,I,J)=SVLATL(:)

        TAUSS(:,I,J)=TAUSSL(:)
        CLDSS(:,I,J)=CLDSSL(:)
        CLDSAV(:,I,J)=CLDSAVL(:)
        CLDSAV1(:,I,J)=CLDSV1(:)
        SVLHX(:,I,J)=SVLHXL(:)
        CSIZSS(:,I,J)=CSIZEL(:)
        QLss_3d(:,i,j)=QLss(:)
        QIss_3d(:,i,j)=QIss(:)
        QLmc_3d(:,i,j)=QLmc(:)
        QImc_3d(:,i,j)=QImc(:)

        IF(USE_VMP) THEN
          TAUSSIP(:,I,J)=TAUSSLIP(:)
          CSIZSSIP(:,I,J)=CSIZELIP(:)
        ENDIF

        RHSAV(:,I,J)=RH(:)
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
        CTEM(:,I,J) =CTEML(:)
        CD3D(:,I,J) =CD3DL(:)
        CL3D(:,I,J) =CL3DL(:)
        CI3D(:,I,J) =CI3DL(:)
#endif
#ifdef CLD_AER_CDNC
        NCL(:,I,J)=NCLL(:)
        NCI(:,I,J)=NCIL(:)
        EGCM(:,I,J) =SME(:)
        CDN3D(:,I,J)=CDN3DL(:)
        CRE3D(:,I,J)=CRE3DL(:)
        CLWP(I,J) = SMLWP
#ifdef TRACERS_AMP
#ifdef BLK_2MOM
        do n=1,nmodes
          NACTV(I,J,:,n) =  NACTC(:,n)
        enddo
#endif
#endif
#endif
        TTOLD(:,I,J)=TH(:)
        QTOLD(:,I,J)=QL(:)

        PREC(I,J)=PRCP            ! total precip mass (kg/m^2)
        EPREC(I,J)=ENRGP          ! energy of precipitation (J/m^2)
        !**** The PRECSS array is only used if a distinction is being made
        !**** between kinds of rain in the ground hydrology.
        PRECSS(I,J)=PRCPSS*100.*BYGRAV  ! large scale precip (kg/m^2)
        !**** accumulate precip specially for SUBDD
        P_acc(I,J)=P_acc(I,J)+PRCP
        PM_acc(I,J)=PM_acc(I,J)+PRCP-PRECSS(I,J)

#ifdef INTERACTIVE_WETLANDS_CH4
        !**** update running-average of precipitation (in mm/day):
        call running_average(prcp*SECONDS_PER_DAY*byDTsrc,I,J,1.d0,n__prec)
#endif

        do L=1,LM
          call inc_ajl(i,j,l,JL_SSHR,SSHR(L))
          !*** Begin Accumulate 3D heating by large scale condensation --
          if(lh_diags.eq.1) then
            AIJL(I,J,L,IJL_LLH)=AIJL(I,J,L,IJL_LLH)+SSHR(L)
            AIJL(I,J,L,IJL_LDRY)=AIJL(I,J,L,IJL_LDRY)+DQLSC(L)
          endif
#if (defined mjo_subdd) || (defined etc_subdd)
          LLH3D(L,I,J)=LLH3D(L,I,J)+SSHR(L)*BYAM(L)
#endif
#ifdef mjo_subdd
          LSCDRY(L,I,J)=LSCDRY(L,I,J)+DQLSC(L)
#endif
          !*** End Accumulate 3D heating by large scale condensation --
          call inc_ajl(i,j,l,JL_MCLDHT,DCTEI(L))
          call inc_ajl(i,j,l,JL_RHE,RH1(L))
          call inc_ajl(i,j,l,JL_CLDSS,CLDSSL(L)*AIRM(L))
          call inc_ajl(i,j,l,JL_CSIZSS,CSIZEL(L)*CLDSSL(L)*AIRM(L))
          !       write(6,*) "CTEM_DRV",CTEML(L),CTEM(I,J,L),L,I,J

          T(I,J,L)=TH(L)*FSSL(L)+TMC(I,J,L)*(1.-FSSL(L))
          Q(I,J,L)=QL(L)*FSSL(L)+QMC(I,J,L)*(1.-FSSL(L))
          SMOM(:,L)=SMOM(:,L)*FSSL(L)+SMOMMC(:,L)*(1.-FSSL(L))
          QMOM(:,L)=QMOM(:,L)*FSSL(L)+QMOMMC(:,L)*(1.-FSSL(L))
          !**** update moment changes
          T3MOM(:,I,J,L)=SMOM(:,L)*BYAM(L)
          Q3MOM(:,I,J,L)=QMOM(:,L)*BYAM(L)
          QCI(I,J,L)=QCIX(L)
          QCL(I,J,L)=QCLX(L)

          !**** CALCULATE WIND TENDENCIES AND STORE IN UKM,VKM
          if(J.eq.1 .and. HAVE_SOUTH_POLE)  then
            do K=1,IM ! KMAX
              UKMSP(K,L)=(UM(K,L)*FSSL(L)+UM1(K,L)*(1.-FSSL(L)))*BYAM(L) &
                   -UKMSP(K,L)
              VKMSP(K,L)=(VM(K,L)*FSSL(L)+VM1(K,L)*(1.-FSSL(L)))*BYAM(L) &
                   -VKMSP(K,L)
            end do
          else if(J.eq.JM .and. HAVE_NORTH_POLE)  then
            do K=1,IM ! KMAX
              UKMNP(K,L)=(UM(K,L)*FSSL(L)+UM1(K,L)*(1.-FSSL(L)))*BYAM(L) &
                   -UKMNP(K,L)
              VKMNP(K,L)=(VM(K,L)*FSSL(L)+VM1(K,L)*(1.-FSSL(L)))*BYAM(L) &
                   -VKMNP(K,L)
            end do
          else
            do K=1,KMAX
              UKM(K,L,I,J)=(UM(K,L)*FSSL(L)+UM1(K,L)*(1.-FSSL(L)))*BYAM(L) &
                   -UKM(K,L,I,J)
              VKM(K,L,I,J)=(VM(K,L)*FSSL(L)+VM1(K,L)*(1.-FSSL(L)))*BYAM(L) &
                   -VKM(K,L,I,J)
            end do
          end if
        enddo

#ifdef CACHED_SUBDD
        !**** calculate change
        do L=1,LM
          dth_ss(I,J,L) = T(I,J,L)-dth_ss(I,J,L)
          dq_ss(I,J,L) = Q(I,J,L)-dq_ss(I,J,L)
        enddo
#endif

        !**** Uncomment next two lines for check on water conservation
        !QCON q2=sum((Q(I,J,:)+WMX(:))*AIRM(:))*100*BYGRAV+PRCP
        !QCON if (abs(q2-q0).gt.1d-13) print*,"water err1",i,j,q2-q0,q2,q0,q1
        !QCON*     ,prcp

        ! Accumulate stratiform cloud droplet and cloud ice number conc and effective radii
        ! diagnostics near cloud top (down to where cumulative layer-average cloud optical
        ! depth > 1), simplistically applying optional temperature cutoffs only
        ! at layer where threshold is passed and ignoring convective cloud opacity to avoid
        ! overlap complications. This is a crude treatment (all the more so because phases
        ! are treated separately for the sake of answers closer to expectations) for purposes
        ! of comparing model runs but not for rigorous comparisons with, say, satellite-borne
        ! measurements. Note that ncl and nci diagnostics are in-cloud values and thus
        ! cloud fraction appears in denominator but not numerator.

        if( svlhxl(lmcld)==lhe )then ! liquid-phase stratiform cloud
          tau_sum = cldssl(lmcld)*taussl(lmcld)
          rcl_sum = cldssl(lmcld)*taussl(lmcld)*csizel(lmcld)
#ifdef CLD_AER_CDNC
          ncl_sum = taussl(lmcld)*ncll(lmcld)
#ifdef TRACERS_AMP
          ccn_sum = cldssl(lmcld)*taussl(lmcld)*sum(nactc(lmcld,:))
#endif
#endif
        else
          tau_sum = 0.
          rcl_sum = 0.
#ifdef CLD_AER_CDNC
          ncl_sum = 0.
#ifdef TRACERS_AMP
          ccn_sum = 0.
#endif
#endif
        endif
        do L = lmcld-1, 1, -1
          if( svlhxl(L)==lhe )then ! liquid-phase stratiform cloud
            if( tau_sum + cldssl(L)*taussl(L) < tau_thresh )then
              ! accumulate sums vertically
              tau_sum = tau_sum + cldssl(L)*taussl(L)
              rcl_sum = rcl_sum + cldssl(L)*taussl(L)*csizel(L)
#ifdef CLD_AER_CDNC
              ncl_sum = ncl_sum + taussl(L)*ncll(L)
#endif
#ifdef TRACERS_AMP
#ifdef BLK_2MOM
              ccn_sum = ccn_sum + cldssl(L)*taussl(L)*sum(nactc(L,:))
#endif
#endif
            else
              ! cumulative optical depth threshold crossed
              dum1 = tau_sum/tau_thresh
              dum2 = 1d0-dum1
              aij(i,j,ij_tclssct) = aij(i,j,ij_tclssct) + &
                 dum1*tau_sum + dum2*cldssl(L)*taussl(L)
              aij(i,j,ij_rclssct) = aij(i,j,ij_rclssct) + &
                 dum1*rcl_sum + dum2*cldssl(L)*taussl(L)*csizel(L)
#ifdef CLD_AER_CDNC
              aij(i,j,ij_nclssct) = aij(i,j,ij_nclssct) + &
                 dum1*ncl_sum + dum2*taussl(L)*ncll(L)
#endif
#ifdef TRACERS_AMP
#ifdef BLK_2MOM
              aij(i,j,ij_ccnssct) = aij(i,j,ij_ccnssct) + &
                 dum1*ccn_sum + dum2*cldssl(L)*taussl(L)*sum(nactc(L,:))
#endif
#endif
              exit
            endif  ! threshold crossed or not
          endif    ! liquid-phase stratiform cloud
        enddo      ! L

#ifdef CLD_AER_CDNC
        do L=1,LM
          CLDWT = CLDSSL(L)!+teeny
          CLDWTDZ = CLDWT*(RGAS*TL(L)*BYGRAV)*(AIRM(L)/PL(L))
          if(SVLHXL(L).eq.LHE) then
            AIJL(I,J,L,IJL_CFWS)= AIJL(I,J,L,IJL_CFWS)+CLDWT
            AIJL(I,J,L,IJL_REWS)= AIJL(I,J,L,IJL_REWS)+AREWS(L)*CLDWT
            AIJL(I,J,L,IJL_CDWS)= AIJL(I,J,L,IJL_CDWS)+ACDNWS(L)*CLDWT
#ifdef TRACERS_TOMAS
            AIJL(I,J,L,IJL_CDTOMAS)= AIJL(I,J,L,IJL_CDTOMAS)+CDNC_TOMAS(L)*CLDWT
#endif
            AIJL(I,J,L,IJL_CWWS)= AIJL(I,J,L,IJL_CWWS)+ALWWS(L)*CLDWT
! standard cdnc
            AIJ(I,J,IJ_DZWS)=AIJ(I,J,IJ_DZWS)+CLDWTDZ
            AIJ(I,J,IJ_3dNWS)=AIJ(I,J,IJ_3dNWS)+ACDNWS(L)*CLDWTDZ
! screened dcnc
!!          AIJ(I,J,IJ_3dNWSS)=AIJ(I,J,IJ_3dNWSS)+ACDNWSS(L)*CLDWTDZ
!
            AIJ(I,J,IJ_3dRWS)=AIJ(I,J,IJ_3dRWS)+AREWS(L)*CLDWTDZ
            AIJ(I,J,IJ_3dLWS)=AIJ(I,J,IJ_3dLWS)+ALWWS(L)*CLDWTDZ
#ifdef CACHED_SUBDD
      Cloud_daily(I,J,1) =  Cloud_daily(I,J,1)+ACDNWS(L)*CLDWTDZ
!!    Cloud_daily(I,J,2) =  Cloud_daily(I,J,2)+ACDNWSS(L)*CLDWTDZ

      Cloud_daily(I,J,10) = Cloud_daily(I,J,10) +AREWS(L)*CLDWTDZ
      Cloud_daily(I,J,14) = Cloud_daily(I,J,14) +CLDWTDZ

      Cloud_daily3d(i,j,l,1) = ACDNWS(L)*CLDWT
      Cloud_daily3d(i,j,l,12) = CLDWT
      Cloud_daily3d(i,j,l,3) = ALWWS(L)*CLDWT
      Cloud_daily3d(i,j,l,8) = AREWS(L)*CLDWT
#endif
          elseif(SVLHXL(L).eq.LHS) then
            AIJL(I,J,L,IJL_CFIS)= AIJL(I,J,L,IJL_CFIS)+CLDWT
            AIJL(I,J,L,IJL_REIS)= AIJL(I,J,L,IJL_REIS)+AREIS(L)*CLDWT
            AIJL(I,J,L,IJL_CDIS)= AIJL(I,J,L,IJL_CDIS)+ACDNIS(L)*CLDWT
            AIJL(I,J,L,IJL_CWIS)= AIJL(I,J,L,IJL_CWIS)+ALWIS(L)*CLDWT
            AIJ(I,J,IJ_DZIS)=AIJ(I,J,IJ_DZIS)+CLDWTDZ
            AIJ(I,J,IJ_3dNIS)=AIJ(I,J,IJ_3dNIS)+ACDNIS(L)*CLDWTDZ
            AIJ(I,J,IJ_3dRIS)=AIJ(I,J,IJ_3dRIS)+AREIS(L)*CLDWTDZ
            AIJ(I,J,IJ_3dLIS)=AIJ(I,J,IJ_3dLIS)+ALWIS(L)*CLDWTDZ
#ifdef CACHED_SUBDD
      Cloud_daily(I,J,12) = Cloud_daily(I,J,12) +AREIS(L)*CLDWTDZ
      Cloud_daily(I,J,16) = Cloud_daily(I,J,16) +CLDWTDZ

      Cloud_daily3d(I,J,l,5) = Cloud_daily3d(I,J,l,5)+ALWIS(L)*CLDWT
      Cloud_daily3d(I,J,l,10) = AREIS(L)*CLDWT
#endif
          endif
          if (NLSW.ge.1) then
            call inc_ajl(i,j,l,JL_CNUMWS,ACDNWS(L)*AIRM(L))
            !     if(AIJ(I,J,IJ_3dNWS).gt.0.)write(6,*)"OUTDRV",AIJ(I,J,IJ_3dNWS)
            !    * ,ACDNWS(L),NLSW,itime,l
          endif
          if(NLSI.ge.1) then
            call inc_ajl(i,j,l,JL_CNUMIS,ACDNIS(L)*AIRM(L))
          endif
        enddo

#endif
        !QCON q2 = sum((Q(I,J,:)+WMX(:))*AIRM(:))*100.*BYGRAV+PRCP
        !QCON if (abs(q0-q2).gt.1d-13) print*,"pr1",i,j,q0,q1,q2,prcp,prcpss*100
        !QCON*     *bygrav
#ifdef TRACERS_ON
        !**** TRACERS: Use only the active ones
        do nx=1,ntx
          n = ntix(nx)

#ifndef SKIP_TRACER_DIAGS
          do l=1,LMCLD
            dtrm(l) = tm(l,nx)-trm(i,j,l,n)*fssl(l)
#ifdef TRACERS_WATER
            dtrm(l) = dtrm(l) + (trwml(nx,l)-trwm(i,j,l,n)-trsvwml(nx,l))
#endif
          enddo
          if(itcon_ss(n).gt.0) call inc_diagtcb(i,j,sum(dtrm(1:LMCLD)), &
               itcon_ss(n),n)
          call inc_tajln_column(i,j,1,LMCLD,lm,jlnt_lscond,n,dtrm)
#endif  /*SKIP_TRACER_DIAGS*/

          do l=1,LMCLD
#ifdef TRACERS_WATER
            trwm(i,j,l,n) = trwml(nx,l)
#endif
            trm(i,j,l,n) = tm(l,nx)+tmsave(l,nx)*(1.-fssl(l))
            trmom(:,i,j,l,n) = tmom(:,l,nx)+tmomsv(:,l,nx)*(1.-fssl(l))
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
            if (ijts_aq(n).gt.0) then ! use ij mask for jl as well
              call inc_tajls(i,j,l,jls_incloud(1,n), &
                   dt_sulf_mc(n,l)*(1.-fssl(l)))
              call inc_tajls(i,j,l,jls_incloud(2,n),dt_sulf_ss(n,l))
              taijs(i,j,ijts_aq(n))=taijs(i,j,ijts_aq(n))+ &
                   dt_sulf_mc(n,l)*(1.-fssl(l))+dt_sulf_ss(n,l)
            end if
#ifdef ACCMIP_LIKE_DIAGS
            if(trname(n).eq."SO4".and.ijlt_prodSO4aq.gt.0) &
            taijls(i,j,l,ijlt_prodSO4aq)=taijls(i,j,l,ijlt_prodSO4aq)+ &
            (dt_sulf_mc(n,l)*(1.-fssl(l))+dt_sulf_ss(n,l))*byaxyp(i,j)
#endif /* ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_AEROSOLS_Koch or TRACERS_AMP or TRACERS_TOMAS */
#ifdef TRACERS_AMP
            if (trname(n).eq."M_ACC_SU") then
              AQsulfRATE(i,j,l)=dt_sulf_mc(n,l)*(1.-fssl(l))+dt_sulf_ss(n,l)
            endif
#endif
#ifdef TRACERS_TOMAS
           if (trname(n).eq."ASO4__01") then
              AQSO4oxid_mc(i,j,l) = dt_sulf_mc(n,l)*(1.-fssl(l))
              AQSO4oxid_ls(i,j,l) = dt_sulf_ss(n,l)
           endif
#endif
          end do
#ifdef TRACERS_WATER
          trprec(n,i,j) = (trprec(n,i,j)+trprss(nx))*byaxyp(i,j)
#ifndef SKIP_TRACER_DIAGS
          TRP_acc(n,I,J)=TRP_acc(n,I,J)+trprec(n,i,j)
#endif
          !        if (i.eq.64.and.j.eq.7) write(6,'(2i3,a,3f12.2)')
          !     .    n,ntm, ' TRP1::ACC:',trp_acc(n,i,j)*byaxyp(i,j),
          !     .    trprec(n,i,j),trprss(nx)
          !**** diagnostics
          if (dowetdep(n)) then
            taijn(i,j,tij_prec,n) =taijn(i,j,tij_prec,n) + &
                 trprec(n,i,j)
#ifndef SKIP_TRACER_DIAGS
            if (jls_prec(1,n).gt.0) call inc_tajls2(i,j,1,jls_prec(1,n), &
                 trprec(n,i,j))
            if (jls_prec(2,n).gt.0) call inc_tajls2(i,j,1,jls_prec(2,n), &
                 trprec(n,i,j)*focean(i,j))
#ifdef TRACERS_COSMO
            if (n .eq. n_Be7) BE7W_acc(i,j)=BE7W_acc(i,j)+ &
                 trprec(n,i,j)
#endif
#endif /*SKIP_TRACER_DIAGS*/
#ifdef TRDIAG_WETDEPO
            !     ..........
            !     accumulates special wet depo diagnostics
            !     ..........
            if (diag_wetdep == 1) then

              if(jls_trdpmc(1,n)>0) call inc_tajls_column(i,j,1,lmcmax,lm, &
                   jls_trdpmc(1,n),trcond_mc(:,nx))
              if(jls_trdpmc(2,n)>0) call inc_tajls_column(i,j,1,lmcmax,lm, &
                   jls_trdpmc(2,n),trdvap_mc(:,nx))
              if(jls_trdpmc(3,n)>0) call inc_tajls_column(i,j,1,lmcmax,lm, &
                   jls_trdpmc(3,n),trflcw_mc(:,nx))
              if(jls_trdpmc(4,n)>0) call inc_tajls_column(i,j,1,lmcmax,lm, &
                   jls_trdpmc(4,n),trprcp_mc(:,nx))
              if(jls_trdpmc(5,n)>0) call inc_tajls_column(i,j,1,lmcmax,lm, &
                   jls_trdpmc(5,n),trnvap_mc(:,nx))
              if(jls_trdpmc(6,n)>0) call inc_tajls_column(i,j,1,lmcmax,lm, &
                   jls_trdpmc(6,n),trwash_mc(:,nx))

              if (ijts_trdpmc(1,n) > 0) taijs(i,j,ijts_trdpmc(1,n)) &
                   =taijs(i,j,ijts_trdpmc(1,n))+sum(trcond_mc(1:lmcmax,nx))
              if (ijts_trdpmc(2,n) > 0) taijs(i,j,ijts_trdpmc(2,n)) &
                   =taijs(i,j,ijts_trdpmc(2,n))+sum(trdvap_mc(1:lmcmax,nx))
              if (ijts_trdpmc(3,n) > 0) taijs(i,j,ijts_trdpmc(3,n)) &
                   =taijs(i,j,ijts_trdpmc(3,n))+sum(trflcw_mc(1:lmcmax,nx))
              if (ijts_trdpmc(4,n) > 0) taijs(i,j,ijts_trdpmc(4,n)) &
                   =taijs(i,j,ijts_trdpmc(4,n))+sum(trprcp_mc(1:lmcmax,nx))
              if (ijts_trdpmc(5,n) > 0) taijs(i,j,ijts_trdpmc(5,n)) &
                   =taijs(i,j,ijts_trdpmc(5,n))+sum(trnvap_mc(1:lmcmax,nx))
              if (ijts_trdpmc(6,n) > 0) taijs(i,j,ijts_trdpmc(6,n)) &
                   =taijs(i,j,ijts_trdpmc(6,n))+sum(trwash_mc(1:lmcmax,nx))

              if(jls_trdpls(1,n) > 0) call inc_tajls_column(i,j,1,LMCLD,lm, &
                   jls_trdpls(1,n),trwash_ls(:,nx))
              if(jls_trdpls(2,n) > 0) call inc_tajls_column(i,j,1,LMCLD,lm, &
                   jls_trdpls(2,n),trprcp_ls(:,nx))
              if(jls_trdpls(3,n) > 0) call inc_tajls_column(i,j,1,LMCLD,lm, &
                   jls_trdpls(3,n),trclwc_ls(:,nx))
              if(jls_trdpls(4,n) > 0) call inc_tajls_column(i,j,1,LMCLD,lm, &
                   jls_trdpls(4,n),trevap_ls(:,nx))
              if(jls_trdpls(5,n) > 0) call inc_tajls_column(i,j,1,LMCLD,lm, &
                   jls_trdpls(5,n),trclwe_ls(:,nx))
              if(jls_trdpls(6,n) > 0) call inc_tajls_column(i,j,1,LMCLD,lm, &
                   jls_trdpls(6,n),trcond_ls(:,nx))

              if (ijts_trdpls(1,n) > 0) taijs(i,j,ijts_trdpls(1,n)) &
                   =taijs(i,j,ijts_trdpls(1,n))+sum(trwash_ls(1:LMCLD,nx))
              if (ijts_trdpls(2,n) > 0) taijs(i,j,ijts_trdpls(2,n)) &
                   =taijs(i,j,ijts_trdpls(2,n))+sum(trprcp_ls(1:LMCLD,nx))
              if (ijts_trdpls(3,n) > 0) taijs(i,j,ijts_trdpls(3,n)) &
                   =taijs(i,j,ijts_trdpls(3,n))+sum(trclwc_ls(1:LMCLD,nx))
              if (ijts_trdpls(4,n) > 0) taijs(i,j,ijts_trdpls(4,n)) &
                   =taijs(i,j,ijts_trdpls(4,n))+sum(trevap_ls(1:LMCLD,nx))
              if (ijts_trdpls(5,n) > 0) taijs(i,j,ijts_trdpls(5,n)) &
                   =taijs(i,j,ijts_trdpls(5,n))+sum(trclwe_ls(1:LMCLD,nx))
              if (ijts_trdpls(6,n) > 0) taijs(i,j,ijts_trdpls(6,n)) &
                   =taijs(i,j,ijts_trdpls(6,n))+sum(trcond_ls(1:LMCLD,nx))
            end if
#endif
#ifdef TRACERS_DUST
            if (adiurn_dust == 1) then
              do kr=1,Ndiupt
                if(i == ijdd(1,kr) .and. j == ijdd(2,kr)) then
                  select case (trname(n))
                  case ('Clay','Silt1','Silt2','Silt3','Silt4','Silt5')
                    tmp(idd_wet)=+trprec(n,i,j)/Dtsrc
                    ADIURN(IDXD(:),KR,IH)=ADIURN(IDXD(:),KR,IH)+ &
                         TMP(IDXD(:))
#ifdef USE_HDIURN
                    HDIURN(IDXD(:),KR,IHM)=HDIURN(IDXD(:),KR,IHM)+ &
                         TMP(IDXD(:))
#endif
                  end select
                end if
              end do
            end if
#endif
          end if
#endif
        end do
#endif

#ifndef TRACERS_WATER
        !     ..........
        !     call simple wet deposition scheme for dust/mineral tracers
        !     ..........

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
        n_fidx=n_soildust

        do n=1,Ntm_dust
          n1=n_fidx+n-1
          do l=1,Lm
            tm_dust(l,n)=trm(i,j,l,n1)
            tmom_dust(:,l,n)=trmom(:,i,j,l,n1)
          end do
        end do

        call dust_wet(i,j)

        do n=1,Ntm_dust
          n1=n_fidx+n-1
          trprec_dust(n,i,j)=0.D0
          do l=1,Lm
            if (itcon_wt(n).gt.0) call inc_diagtcb(i,j, &
                 tm_dust(l,n)-trm(i,j,l,n1),itcon_wt(n),n)
            trm(i,j,l,n1)=tm_dust(l,n)
            trmom(:,i,j,l,n1)=tmom_dust(:,l,n)
            trprec_dust(n,i,j)=trprec_dust(n,i,j)+trprc_dust(l,n)
            call inc_tajls(i,j,l,jls_wet(n1),trprc_dust(l,n))
            taijs(i,j,ijts_wet(n1))=taijs(i,j,ijts_wet(n1)) &
                 +trprc_dust(l,n)
          end do
        end do

#endif

#ifdef TRACERS_DUST
        if (adiurn_dust == 1) then
          do n=1,Ntm_dust
            do kr=1,Ndiupt
              if(i == ijdd(1,kr) .and. j == ijdd(2,kr)) then
                select case (trname(n))
                case ('Clay','Silt1','Silt2','Silt3','Silt4','Silt5')
                  tmp(idd_wet)=+trprec_dust(n,i,j)*byaxyp(i,j)/Dtsrc
                  ADIURN(IDXD(:),KR,IH)=ADIURN(IDXD(:),KR,IH)+ &
                       TMP(IDXD(:))
#ifdef USE_HDIURN
                  HDIURN(IDXD(:),KR,IHM)=HDIURN(IDXD(:),KR,IHM)+ &
                       TMP(IDXD(:))
#endif
                end select
              end if
            end do
          end do
        end if
#endif
#endif

      end do
      !**** END OF MAIN LOOP FOR INDEX I

    end do ! loop over threads


    ! Burn random numbers for later longitudes here.
    ! Actual generation of random numbers is in CLOUDS2.f::ISCCP_CLOUD_TYPES
    if (isccp_diags.eq.1) then
      call BURN_RANDOM(nij_after_i1(I_1)*NCOL*(LM+1))
    end if

  end do
  !**** END OF MAIN LOOP FOR INDEX J

  !****
  !
  !     WAS THERE AN ERROR IN SUBSID ??
  !
  if(ICKERR.ne.0)  then
    write(6,*)  'SUBSID ERROR: ABS(C) > 1'
    call stop_model('SUBSID ERROR: ABS(C) > 1',255)
  end if
  !
  !     WAS THERE AN ERROR IN ISCCP CLOUD TYPING ??
  !
  if(JCKERR.ne.0)  then
    write(6,*)  'ISCCP CLOUD TYPING ERROR'
    call stop_model('ISCCP CLOUD TYPING ERROR',255)
  end if

#ifdef TRACERS_TOMAS
!C     To fix inconsistent aerosol size distribution and water eqm.
      CALL aeroupdate
#endif
#ifdef CACHED_SUBDD
#ifdef COSP_SIM
  ! Trigger COSP (CFMIP Observation Simulator Package (Only sample nsubdd_cosp)
  if (mod(itime+1,nsubdd_cosp).eq.0) then
    call run_cosp_sims()
  endif
#endif
#endif
  !
  !     NOW UPDATE THE MODEL WINDS
  !

  call avg_replicated_duv_to_vgrid(ukm,vkm,kmax_nonpolar, &
       ukmsp,vkmsp,ukmnp,vkmnp)

  !**** ADD IN CHANGE OF MOMENTUM BY MOIST CONVECTION AND CTEI
  UASV(:,I_0:I_1,J_0S:J_1S) = UA(:,I_0:I_1,J_0S:J_1S)
  call recalc_agrid_uv ! add option for tendency computation?
  do J=J_0S,J_1S
    do I=I_0,I_1
      do L=1,LM
        call inc_ajl(i,j,l,JL_DAMMC, &
             (UA(L,I,J)-UASV(L,I,J))*PDSIG(L,I,J))
      end do
    end do
  end do

  if (isccp_diags.eq.1) call RINIT(seed) ! reset random number sequ.

415 format(1X,'W500 AT I=21 L=5 TIME= ',I10/,1X,10F8.3/,1X,10F8.3)
420 format(1X,'ENT  AT I=21 L=5'/,1X,10F8.2/,1X,10F8.2)

#ifdef CACHED_SUBDD
!****
!**** Collect some high-frequency outputs
!****
  call find_groups('aijh',grpids,ngroups)
  do igrp=1,ngroups
  subdd => subdd_groups(grpids(igrp))
  do k=1,subdd%ndiags
  select case (subdd%name(k))
  case ('prec')
    call inc_subdd(subdd,k,prec)
  case ('ssp')
    call inc_subdd(subdd,k,precss)
  case ('mcp')
    do j=j_0,j_1; do i=i_0,imaxj(j)
      sddarr(i,j) = max(0.,prec(i,j)-precss(i,j))
    enddo;        enddo
    call inc_subdd(subdd,k,sddarr)
  case ('snowfall')
    do j=j_0,j_1; do i=i_0,imaxj(j)
      if(eprec(i,j).ge.0.) then
        sddarr(i,j) = 0.
      else
        sddarr(i,j) = prec(i,j)
      endif
    enddo;        enddo
    call inc_subdd(subdd,k,sddarr)
#ifdef CFMIP3_SUBDD
  case ('mc_lwp')
    call inc_subdd(subdd,k,cfmip_mc_lwp)
#endif
  end select
  enddo
  enddo

  call find_groups('aijlh',grpids,ngroups)
  do igrp=1,ngroups
  subdd => subdd_groups(grpids(igrp))
  do k=1,subdd%ndiags
  select case (subdd%name(k))
  case ('qcl')
    call inc_subdd(subdd,k,qcl)
  case ('qci')
    call inc_subdd(subdd,k,qci)
  case ('cldss')
    call inc_subdd(subdd,k,cldss,jdim=3)
  case ('cldmc')
    call inc_subdd(subdd,k,cldmc,jdim=3)
#ifdef CFMIP3_SUBDD
  case ('mcamfx')
    call inc_subdd(subdd,k,cfmip_mcamfx)
#endif
  end select
  enddo
  enddo

  call find_groups('cijh',grpids,ngroups)
  do igrp=1,ngroups
  subdd => subdd_groups(grpids(igrp))
  do k=1,subdd%ndiags
  select case (subdd%name(k))
  case ('isccp_sunlit')
    call inc_subdd(subdd,k,saveSCLDI)
  case ('isccp_ctp')
    call inc_subdd(subdd,k,saveCTPI)
  case ('isccp_tau')
    call inc_subdd(subdd,k,saveTAUI)
  case ('isccp_lcld')
    call inc_subdd(subdd,k,saveLCLDI)
  case ('isccp_mcld')
    call inc_subdd(subdd,k,saveMCLDI)
  case ('isccp_hcld')
    call inc_subdd(subdd,k,saveHCLDI)
#ifdef CFMIP3_SUBDD  /* CFMIP3_SUBDD */
  case ('ctp_mc')
    call inc_subdd(subdd,k,cfmip_ctp_mc)
  case ('cbp_mc')
    call inc_subdd(subdd,k,cfmip_cbp_mc)
  case ('cnvfrq')
    call inc_subdd(subdd,k,cfmip_cnvfrq)
  case ('dcnvfrq')
    call inc_subdd(subdd,k,cfmip_dcnvfrq)
  case ('scnvfrq')
    call inc_subdd(subdd,k,cfmip_scnvfrq)
  case ('mc_twp')
    call inc_subdd(subdd,k,cfmip_mc_twp)
  ! case ('mc_lwp')
  !   call inc_subdd(subdd,k,cfmip_mc_lwp)
  case ('qatm')
    call inc_subdd(subdd,k,cfmip_wvp)
#endif  /* CFMIP3_SUBDD */
  end select
  enddo
  enddo

  call find_groups('cijlh',grpids,ngroups)
  do igrp=1,ngroups
  subdd => subdd_groups(grpids(igrp))
  do k=1,subdd%ndiags
  select case (subdd%name(k))
  case ('tau_ss')
    do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
      sddarr3d(i,j,l) = tauss(l,i,j)
    enddo;        enddo;        enddo
    call inc_subdd(subdd,k,sddarr3d)
  case ('tau_mc')
    do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
      sddarr3d(i,j,l) = taumc(l,i,j)
    enddo;        enddo;        enddo
    call inc_subdd(subdd,k,sddarr3d)
  case ('dq_mc')
    call inc_subdd(subdd,k,dq_mc)
  case ('dth_mc')
    call inc_subdd(subdd,k,dth_mc)
  case ('dq_ss')
    call inc_subdd(subdd,k,dq_ss)
  case ('dth_ss')
    call inc_subdd(subdd,k,dth_ss)
  end select
  enddo
  enddo
#ifdef CACHED_SUBDD
#ifdef COSP_SIM
  ! Save requested COSP variables to disk (instantaneous)
  ! using the SUBDD framework (Only sample nsubdd_cosp)
  if (mod(Itime+1,nsubdd_cosp).eq.0) then
    ! Save requested COSP variables using the SUBDD framework.
    call save_cosp(grid)
    ! Release COSP data structure memory
    call free_cosp_sims()
  endif
#endif
#endif

#ifdef TRACERS_WATER

  call find_groups('taijh',grpids,ngroups) !2-D tracer variables
  do igrp=1,ngroups
    subdd => subdd_groups(grpids(igrp))
    do k=1,subdd%ndiags
      ntm_loop: do n=1,ntm
        !Set precipitation tracer name (sname):
        trpname = trim(trname(n))//'_in_prec'
        !If name matches subdd name, then add to output and exit loop:
        if(trpname.eq.trim(subdd%name(k))) then
          call inc_subdd(subdd,k,trprec(n,:,:))
          exit ntm_loop
        end if
      end do ntm_loop
    end do !subdd diagnostics/variables
  end do   !subdd groups

#endif

#ifdef SCM
  if (isccp_diags.eq.1) then
      call inc_subdd('isccp_fq',save_fq_isccp,1,.true.,units='fraction', &
           long_name='Cld Fct by ISCCP CldTypes',dim3name='ntau',dim4name='npres')
  endif

  if( SCMopt%PlumeDiag )then
  ! plume diagnostics
    call inc_subdd('mc_mfu_p1',mc_mfu_p1,1,.true.,units='kg/m2/s', &
         long_name='Plume 1 mass flux',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_mfu_p2',mc_mfu_p2,1,.true.,units='kg/m2/s', &
         long_name='Plume 2 mass flux',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_mfd_p1',mc_mfd_p1,1,.true.,units='kg/m2/s', &
         long_name='Plume 1 downdraft mass flux',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_mfd_p2',mc_mfd_p2,1,.true.,units='kg/m2/s', &
         long_name='Plume 2 downdraft mass flux',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_w_p1',mc_w_p1,1,.true.,units='m/s', &
         long_name='Plume 1 updraft speed',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_w_p2',mc_w_p2,1,.true.,units='m/s', &
         long_name='Plume 2 udraft speed',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_ent_p1',mc_ent_p1,1,.true.,units='%/km', &
         long_name='Plume 1 entrainment rate',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_ent_p2',mc_ent_p2,1,.true.,units='%/km', &
         long_name='Plume 2 entrainment rate',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_det_p1',mc_det_p1,1,.true.,units='%/km', &
         long_name='Plume 1 detrainment rate',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_det_p2',mc_det_p2,1,.true.,units='%/km', &
         long_name='Plume 2 detrainment rate',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_m_p1',mc_m_p1,1,.true.,units='mb', &
         long_name='Plume 1 mass',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_m_p2',mc_m_p2,1,.true.,units='mb', &
         long_name='Plume 2 mass',dim3name='level',dim4name='base_level')
    call inc_subdd('mc_pl_max_p1',mc_pl_max_p1,1,.true.,units='mb', &
         long_name='Plume 1 maximum level',dim3name='base_level')
    call inc_subdd('mc_pl_max_p2',mc_pl_max_p2,1,.true.,units='mb', &
         long_name='Plume 2 maximum level',dim3name='base_level')
    call inc_subdd('mc_pl_min_p1',mc_pl_min_p1,1,.true.,units='mb', &
         long_name='Plume 1 minimum level',dim3name='base_level')
    call inc_subdd('mc_pl_min_p2',mc_pl_min_p2,1,.true.,units='mb', &
         long_name='Plume 2 minimum level',dim3name='base_level')
  endif
#endif

#endif

#if (defined CALCULATE_LIGHTNING) || (defined TRACERS_SPECIAL_Shindell)
#ifdef AUTOTUNE_LIGHTNING
      CALL GLOBALSUM( grid, FLASH_UNC*(FOCEAN), SEA_FR_UNC(1), ALL=.true. )
      CALL GLOBALSUM( grid, FLASH_UNC*(1d0-FOCEAN), LAND_FR_UNC(1), ALL=.true. )
      CNT_FR(1) = 1
#endif
#endif

  call stopTimer('CONDSE()')

  return
end subroutine CONDSE

subroutine init_CLD(istart)
!@sum  init_CLD initialises parameters for MSTCNV and LSCOND
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
  use CONSTANT, only : grav,by3,radian,lhe
  use RESOLUTION, only : ls1=>ls1_nominal,plbot
  use RESOLUTION, only : jm,lm
  use MODEL_COM, only : dtsrc
  USE ATM_COM, only : t,q,pmid,pk ! for coldstart istart=2 case
  use DOMAIN_DECOMP_ATM, only : GRID, AM_I_ROOT
  use GEOM, only : lat2d
#ifndef SCM
  use GEOM, only : kmaxj
#endif
#if(defined CALCULATE_LIGHTNING)||(defined TRACERS_SPECIAL_Shindell)
  use LIGHTNING, only : tune_lt_land, tune_lt_sea, FLASH_PERTURB
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE LIGHTNING, only : FLASH_YIELD_MIDLAT, FLASH_YIELD_TROPIC
#endif
  use CLOUDS, only : lmcm,bydtsrc,xmass,brcld,bybr,U00wtrX,U00ice &
       ,U00a,U00b       & ! tuning knobs to replace U00ice and U00wtrX
       ,MAXCTOP         & ! the maximum cloud top (mb)
       ,HRMAX,ISC,LMCLD,RICldX,RWCldOX,xRIcld,do_blU00,tautab,invtau &
       ,funio_denominator,autoconv_multiplier,radius_multiplier &
       ,entrainment_cont1,entrainment_cont2,wmui_multiplier &
       ,RA,UM,VM,UM1,VM1,U_0,V_0
#ifdef TRACERS_ON
  use TRACER_COM, only: NTM
  use QUSDEF, only: nmom

#ifdef TRACERS_WATER
  use CLOUDS, only : trwml,trsvwml,trprmc,trprss
  use CLOUDS, only : TRCOND, TRCONDV
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
  use CLOUDS, only : dt_sulf_mc,dt_sulf_ss
#endif
#endif
#ifdef TRDIAG_WETDEPO
  use CLOUDS, only : &
        trcond_mc,trdvap_mc,trflcw_mc,trprcp_mc,trnvap_mc,trwash_mc &
       ,trcond_ls,trevap_ls,trclwc_ls,trprcp_ls,trclwe_ls,trwash_ls
#endif
  use CLOUDS, only: ntix, TM, TMOM, TRDNL &
       ,DTM, DTMR, TMDNL, DTMOM, DTMOMR, TMOMDNL
#endif

  use CLOUDS_COM, only : llow,lmid,lhi &
       ,isccp_reg2d,UKM,VKM,ttold,qtold,rhsav
  use CLOUDS, only : use_vmp
  use DIAG_COM, only : nisccp,isccp_late &
       ,isccp_diags,ntau,npres
  use ATM_COM, only : pednl00 ! use plbot instead of pednl00
  use Dictionary_Mod
  use FILEMANAGER, only : openunit, closeunit
#ifdef BLK_2MOM
  use resolution, only : im
  USE mo_bulk2m_driver_gcm, ONLY: init_bulk2m_driver
#ifdef TRACERS_AMP
  USE AERO_CONFIG, ONLY: NMODES
#endif
#endif
#ifdef CACHED_SUBDD
#ifdef COSP_SIM
  ! Import modules needed for the COSP simulator
  use MODEL_COM, only : itime
  use domain_decomp_atm, only : hasSouthPole, hasNorthPole
  use resolution, only : im
  use clouds_com, only : npoints_cosp,flag_pfluxes_cosp,flag_re_cosp
  use cosp_drv, only : init_cosp
#endif
#endif
  implicit none
  integer, intent(in) :: istart
  real*8 PLE
  real*8 qsat ! external function
  integer L,I,J,n,iu_ISCCP
  integer :: I_0,I_1,J_0,J_1, I_0H,I_1H,J_0H,J_1H
  character TITLE*80
#ifdef BLK_2MOM
  LOGICAL      :: ldummy=.false.
  INTEGER      :: il0,jl0,kl0,nm0
  CHARACTER*20 :: bname='kuku.txt'
  CHARACTER*15 :: sname='MODELE_mainV3: '
#endif
#ifdef SCM
  integer :: kmaxj(1)=1
#endif

  I_0 =GRID%I_STRT
  I_1 =GRID%I_STOP
  J_0 =GRID%J_STRT
  J_1 =GRID%J_STOP
  I_0H =GRID%I_STRT_HALO
  I_1H =GRID%I_STOP_HALO
  J_0H =GRID%J_STRT_HALO
  J_1H =GRID%J_STOP_HALO

#ifdef BLK_2MOM
  kl0=12;il0=im;jl0=jm ;il0=1;jl0=1;kl0=1
  nm0 = 1 ! or whatever, put a correct value here
#ifdef TRACERS_AMP
  nm0=NMODES
#endif
  ldummy = init_bulk2m_driver(dtsrc,il0,jl0,kl0,nm0,bname)
  if(ldummy) then
    print *,sname,'BLK Initialization is completed...'
  else
    call stop_model("BLK Initialization is not completed: ",255)
  endif
#endif

  if(istart==2) then ! replace with cold vs warm start logic
    do l=1,lm
      ttold(l,:,:)=t(:,:,l)
      qtold(l,:,:)=q(:,:,l)
#ifdef RH_INIT_FROM_TandQ /* define init RH(T,Q) w/r/t liquid instead of using fixed value */
      do j=j_0,j_1
      do i=i_0,i_1
        rhsav(l,i,j)=q(i,j,l)/qsat(t(i,j,l)*pk(l,i,j),lhe,pmid(l,i,j))
      enddo
      enddo
#endif
    end do
  endif

  ! allocate misc variables sized with number of tracers
#ifdef TRACERS_ON
  allocate(ntix(ntm))
  allocate(TM(LM,NTM))
  allocate(TMOM(nmom,lm,ntm))
  allocate(TRDNL(NTM,LM))
#ifdef TRACERS_WATER
  allocate(TRWML(NTM,LM))
  allocate(TRSVWML(NTM,LM))
  allocate(TRPRSS(NTM))
  allocate(TRPRMC(NTM))
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
  ! for diagnostics
    allocate(DT_SULF_MC(NTM,LM))
    allocate(DT_SULF_SS(NTM,LM))
#endif
#ifdef TRDIAG_WETDEPO
    allocate(trcond_mc(LM,NTM))
    allocate(trdvap_mc(LM,NTM))
    allocate(trflcw_mc(LM,NTM))
    allocate(trprcp_mc(LM,NTM))
    allocate(trnvap_mc(LM,NTM))
    allocate(trwash_mc(LM,NTM))

    allocate(trcond_ls(LM,NTM))
    allocate(trevap_ls(LM,NTM))
    allocate(trclwc_ls(LM,NTM))
    allocate(trprcp_ls(LM,NTM))
    allocate(trclwe_ls(LM,NTM))
    allocate(trwash_ls(LM,NTM))
#endif
#else
#endif

    allocate(DTM(LM,NTM)); DTM = 0
    allocate(DTMR(LM,NTM)); DTMR = 0
    allocate(TMDNL(LM,NTM)); TMDNL = 0
    allocate(DTMOM(NMOM,LM,NTM)); DTMOM = 0
    allocate(DTMOMR(NMOM,LM,NTM)); DTMOMR = 0
    allocate(TMOMDNL(NMOM,LM,NTM)); TMOMDNL = 0
#ifdef TRACERS_WATER
    allocate(TRCOND(NTM,LM)); TRCOND = 0
    allocate(TRCONDV(NTM,LM)); TRCONDV = 0
#endif
#endif

  !
  ! allocate space for the varying number of staggered
  ! wind data to be vertically mixed by clouds on the A grid
  !
  n = maxval(kmaxj(j_0:j_1))
  allocate(RA(n))
#ifdef SCM
  RA(:) = 1D0
#endif
  allocate(UM(n,lm),VM(n,lm),UM1(n,lm),VM1(n,lm))
  allocate(U_0(n,lm),V_0(n,lm))
  n = minval(kmaxj(j_0:j_1))
  allocate(UKM(n,lm,i_0h:i_1h,j_0h:j_1h), &
       VKM(n,lm,i_0h:i_1h,j_0h:j_1h))

  call sync_param( 'U00wtrX', U00wtrX )
  call sync_param( 'U00ice', U00ice )
  call sync_param( 'U00a', U00a )
  call sync_param( 'U00b', U00b )
  call sync_param( "LMCM", LMCM )
  call sync_param( "MAXCTOP", MAXCTOP )
  call sync_param( "HRMAX", HRMAX )
  call sync_param( "RICldX", RICldX )
  xRIcld = .001d0*(RICldX-1.d0)
  call sync_param( "RWCldOX", RWCldOX )
  call sync_param( "ISC", ISC)
  call sync_param( "do_blU00", do_blU00)
  call sync_param( "funio_denominator",funio_denominator)
  call sync_param( "autoconv_multiplier",autoconv_multiplier)
  call sync_param( "radius_multiplier",radius_multiplier)
  call sync_param( "wmui_multiplier",wmui_multiplier)
  call sync_param( "entrainment_cont1",entrainment_cont1)
  call sync_param( "entrainment_cont2",entrainment_cont2)
#if(defined CALCULATE_LIGHTNING)||(defined TRACERS_SPECIAL_Shindell)
#ifndef AUTOTUNE_LIGHTNING
  call sync_param( "tune_lt_land", tune_lt_land)
  call sync_param( "tune_lt_sea" , tune_lt_sea )
#endif
#ifdef TRACERS_SPECIAL_Shindell
  call sync_param( "FLASH_YIELD_MIDLAT", flash_yield_midlat )
  call sync_param( "FLASH_YIELD_TROPIC", flash_yield_tropic )
#endif
  call sync_param( "FLASH_PERTURB", flash_perturb )
#endif

  if(LMCM.lt.0) LMCM = LS1-1
  call set_param( "LMCM", LMCM, 'o' )

  i = 0
  call sync_param('use_vmp',i)
  use_vmp = i==1

  BYDTsrc=1./DTsrc
  XMASS=0.1d0*DTsrc*GRAV

  BYBR=((1.-BRCLD)*(1.-2.*BRCLD))**BY3

  !**** SEARCH FOR THE 50 MB LEVEL
  LMCLD=LM
  do L=LM-1,1,-1
    PLE=.25*(PEDNL00(L)+2.*PEDNL00(L+1)+PEDNL00(L+2))
    if (PLE.lt.MAXCTOP) LMCLD=L
  end do
  if (AM_I_ROOT())  write(6,*) &
       "Maximum level for LSCOND calculations (MAXCTOP mb): ",LMCLD

  !**** CLOUD LAYER INDICES USED FOR DIAGNOSTICS (MATCHES ISCCP DEFNs)
  do L=1,LM
    LLOW=L
    if (.5*(PLbot(L+1)+PLbot(L+2)).lt.680.) exit
  end do
  do L=LLOW+1,LM
    LMID=L
    if (.5*(PLbot(L+1)+PLbot(L+2)).lt.440.) exit
  end do
  LHI=LM
  if (LMID+1.gt.LHI) LHI=LMID+1
  if (AM_I_ROOT()) write (6,47) LLOW,LLOW+1,LMID,LMID+1,LHI
47 format (' LOW CLOUDS IN LAYERS 1-',I2,'   MID LEVEL CLOUDS IN', &
       ' LAYERS',I3,'-',I2,'   HIGH CLOUDS IN LAYERS',I3,'-',I2)

  call sync_param( "isccp_diags",isccp_diags)
  if (isccp_diags.eq.1) then

    ! allocate/define distributed 2D ISCCP arrays
    !**** Define regions for ISCCP diagnostics
    allocate(isccp_reg2d(i_0h:i_1h,j_0h:j_1h))
    do j=j_0,j_1
    do i=i_0,i_1
      isccp_reg2d(i,j)=0
      do n=1,nisccp
        if(dble(nint(lat2d(i,j)/radian)).ge.isccp_late(n) .and. &
             dble(nint(lat2d(i,j)/radian)).lt.isccp_late(n+1)) then
          isccp_reg2d(i,j)=n
          exit
        endif
      enddo
    enddo
    enddo

    !**** Read in tau/invtau tables for ISCCP calculations
    call openunit("ISCCP",iu_ISCCP,.true.,.true.)
    read(iu_ISCCP) title,tautab,invtau
    if (AM_I_ROOT())  write(6,*) "Read ISCCP:",trim(title)
    call closeunit(iu_ISCCP)

  endif
#ifdef CACHED_SUBDD
#ifdef COSP_SIM
  ! Number of grids in this domain, pole free, skipping halo
  npoints_cosp = (1+I_1-I_0)*(1+J_1-J_0)
  ! Deal with possible poles
  if (hassouthpole(grid)) then
    ! Domain contains a pole w/ single longitude
    npoints_cosp = npoints_cosp - im + 1
  end if
  if (hasnorthpole(grid)) then
    ! domain contains a pole w/ single longitude
    npoints_cosp = npoints_cosp - im + 1
  end if
  ! Initialize COSP
  call init_cosp(npoints_cosp,im,lm,flag_pfluxes_cosp,flag_re_cosp)
#endif
#endif
end subroutine init_CLD

subroutine qmom_topo_adjustments
  !
  ! Modifies "horizontal" moments of humidity and temperature above steep
  ! topographic slopes to prevent large supersaturations in upslope flow.
  !
  use constant, only : tf,lhe,lhs,bysha
  use resolution, only : ls1=>ls1_nominal
  use resolution, only : im,jm,lm
  use atm_com, only : zatmo,t,q
  use atm_com, only : MUs,MVs,pk,pmid
  use qusdef, only : mx,mxx,my,myy
  use somtq_com, only : tmom,qmom
  use domain_decomp_atm, only : grid,getDomainBounds,halo_update
#ifdef TRACERS_WATER
  use OldTracer_mod, only: tr_wd_type, nWATER
  use tracer_com, only: trm,trmom,NTM
#endif
  use clouds_com, only : svlhx
  implicit none
  integer :: i,j,l
  integer :: iloop_min,iloop_max,jloop_min,jloop_max, &
       ioff_MUs,joff_MVs
  real*8 :: ttmp,qtmp,slh,zthresh,lhx, &
       qe1,qe2,qe1_sv,qe2_sv, &
       te1,te2,te1_sv,te2_sv
  real*8, dimension(lm) :: tl,pl
  real*8 :: qsat ! external function
  real*8, parameter :: qxs=0.1d0
#ifdef TRACERS_WATER
  integer :: n
  real*8 :: ratio
#endif
  integer :: J_0,J_1,J_0S,J_1S,I_0,I_1
  logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

  !**** define local grid
  call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1, &
       J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S, &
       HAVE_NORTH_POLE=HAVE_NORTH_POLE, &
       HAVE_SOUTH_POLE=HAVE_SOUTH_POLE        )
  I_0 = grid%I_STRT
  I_1 = grid%I_STOP

  call halo_update(grid,t)
  call halo_update(grid,pk,jdim=3)   ! already haloed?
  call halo_update(grid,pmid,jdim=3) ! already haloed?
#ifdef CUBED_SPHERE
  ! MVs is already haloed
#else
  call halo_update(grid,MVs)
#endif
  call halo_update(grid,svlhx,jdim=3)

  if(have_south_pole) then
    jloop_min=2 ! horizontal gradients are zero on polar caps
  else
    jloop_min=j_0
  endif
  if(have_north_pole) then
    jloop_max=jm-1 ! horizontal gradients are zero on polar caps
  else
    jloop_max=j_1
  endif
  if(i_0.eq.grid%i_strt_halo) then ! latlon model
    iloop_min = 2                  ! skip IDL
    iloop_max = im-1
    ioff_MUs = 0
    joff_MVs = 0
  else   ! for now, this case is assumed to be the cubed sphere
    iloop_min = i_0
    iloop_max = i_1
    ioff_MUs = 1
    joff_MVs = 1
  endif
  do j=jloop_min,jloop_max
    do i=iloop_min,iloop_max
      zthresh = zatmo(i,j) + 4000. ! ~400 m
      !
      ! north-south
      !
      if(zatmo(i,j-1).gt.zthresh .or. zatmo(i,j+1).gt.zthresh) then
        do l=1,ls1-1
          tl(l) = t(i,j,l)*pk(l,i,j)
          !            if(tl(l).lt.tf) exit ! only apply to relatively moist layers
          if(q(i,j,l).le.0.) cycle
          pl(l) = pmid(l,i,j)
          qe1_sv = q(i,j,l)-qmom(my,i,j,l)+qmom(myy,i,j,l)
          qe2_sv = q(i,j,l)+qmom(my,i,j,l)+qmom(myy,i,j,l)
          qe1 = qe1_sv
          qe2 = qe2_sv
          te1 = t(i,j,l)-tmom(my,i,j,l)+tmom(myy,i,j,l)
          te2 = t(i,j,l)+tmom(my,i,j,l)+tmom(myy,i,j,l)
          if(zatmo(i,j-1).gt.zthresh .and. &
               MVs(i,j-1+joff_MVs,l).lt.0.) then
            ttmp = t(i,j,l)*pk(l,i,j-1); qtmp = q(i,j,l)
            lhx = svlhx(l,i,j-1)
            if(lhx.eq.0.) then
              if(tl(l).gt.tf) then
                lhx = lhe
              else
                lhx = lhs
              endif
            endif
            slh=lhx*bysha
            call moist_adiabat_tq(ttmp,qtmp,lhx,pmid(l,i,j-1))
            ttmp = ttmp-qxs*qtmp*slh
            qtmp = qtmp*(1.+qxs)
            if(qtmp.lt.qe1) then
              qe1 = qtmp
              te1 = ttmp/pk(l,i,j-1)
            endif
          endif
          if(zatmo(i,j+1).gt.zthresh .and. &
               MVs(i,j+joff_MVs,l).gt.0.) then
            ttmp = t(i,j,l)*pk(l,i,j+1); qtmp = q(i,j,l)
            lhx = svlhx(l,i,j+1)
            if(lhx.eq.0.) then
              if(tl(l).gt.tf) then
                lhx = lhe
              else
                lhx = lhs
              endif
            endif
            slh=lhx*bysha
            call moist_adiabat_tq(ttmp,qtmp,lhx,pmid(l,i,j+1))
            ttmp = ttmp-qxs*qtmp*slh
            qtmp = qtmp*(1.+qxs)
            if(qtmp.lt.qe2) then
              qe2 = qtmp
              te2 = ttmp/pk(l,i,j+1)
            endif
          endif
          if(qe1.lt.qe1_sv .or. qe2.lt.qe2_sv) then
            qmom(my ,i,j,l) = .5*(qe2-qe1)
            qmom(myy,i,j,l) = .5*(qe2+qe1)-q(i,j,l)
            tmom(my ,i,j,l) = .5*(te2-te1)
            tmom(myy,i,j,l) = .5*(te2+te1)-t(i,j,l)
#ifdef TRACERS_WATER
            do n=1,ntm
              if(tr_wd_type(n) .eq. nWater) then
                ratio = trm(i,j,l,n)/q(i,j,l)
                trmom(my ,i,j,l,n) = ratio*qmom(my ,i,j,l)
                trmom(myy,i,j,l,n) = ratio*qmom(myy,i,j,l)
              endif
            enddo
#endif
          endif
        enddo
      endif
      !
      ! east-west
      !
      if(zatmo(i-1,j).gt.zthresh .or. zatmo(i+1,j).gt.zthresh) then
        do l=1,ls1-1
          tl(l) = t(i,j,l)*pk(l,i,j)
          !            if(tl(l).lt.tf) exit ! only apply to relatively moist layers
          if(q(i,j,l).le.0.) cycle
          pl(l) = pmid(l,i,j)
          qe1_sv = q(i,j,l)-qmom(mx,i,j,l)+qmom(mxx,i,j,l)
          qe2_sv = q(i,j,l)+qmom(mx,i,j,l)+qmom(mxx,i,j,l)
          qe1 = qe1_sv
          qe2 = qe2_sv
          te1 = t(i,j,l)-tmom(mx,i,j,l)+tmom(mxx,i,j,l)
          te2 = t(i,j,l)+tmom(mx,i,j,l)+tmom(mxx,i,j,l)
          if(zatmo(i-1,j).gt.zthresh .and. &
               MUs(i-1+ioff_MUs,j,l).lt.0.) then
            ttmp = t(i,j,l)*pk(l,i-1,j); qtmp = q(i,j,l)
            lhx = svlhx(l,i-1,j)
            if(lhx.eq.0.) then
              if(tl(l).gt.tf) then
                lhx = lhe
              else
                lhx = lhs
              endif
            endif
            slh=lhx*bysha
            call moist_adiabat_tq(ttmp,qtmp,lhx,pmid(l,i-1,j))
            ttmp = ttmp-qxs*qtmp*slh
            qtmp = qtmp*(1.+qxs)
            if(qtmp.lt.qe1) then
              qe1 = qtmp
              te1 = ttmp/pk(l,i-1,j)
            endif
          endif
          if(zatmo(i+1,j).gt.zthresh .and. &
               MUs(i+ioff_MUs,j,l).gt.0.) then
            ttmp = t(i,j,l)*pk(l,i+1,j); qtmp = q(i,j,l)
            lhx = svlhx(l,i+1,j)
            if(lhx.eq.0.) then
              if(tl(l).gt.tf) then
                lhx = lhe
              else
                lhx = lhs
              endif
            endif
            slh=lhx*bysha
            call moist_adiabat_tq(ttmp,qtmp,lhx,pmid(l,i+1,j))
            ttmp = ttmp-qxs*qtmp*slh
            qtmp = qtmp*(1.+qxs)
            if(qtmp.lt.qe2) then
              qe2 = qtmp
              te2 = ttmp/pk(l,i+1,j)
            endif
          endif
          if(qe1.lt.qe1_sv .or. qe2.lt.qe2_sv) then
            qmom(mx ,i,j,l) = .5*(qe2-qe1)
            qmom(mxx,i,j,l) = .5*(qe2+qe1)-q(i,j,l)
            tmom(mx ,i,j,l) = .5*(te2-te1)
            tmom(mxx,i,j,l) = .5*(te2+te1)-t(i,j,l)
#ifdef TRACERS_WATER
            do n=1,ntm
              if(tr_wd_type(n) .eq. nWater) then
                ratio = trm(i,j,l,n)/q(i,j,l)
                trmom(mx ,i,j,l,n) = ratio*qmom(mx ,i,j,l)
                trmom(mxx,i,j,l,n) = ratio*qmom(mxx,i,j,l)
              endif
            enddo
#endif
          endif
        enddo
      endif
    enddo
  enddo

  return
contains

  subroutine moist_adiabat_tq(t,q,lhx,pl)
    use constant, only : bysha
    implicit none
!@var t,q temperature and specific humidity
    real*8, intent(inout) :: t,q
!@var lhx latent heat of phase change
!@var pl pressure
    real*8, intent(in) :: lhx,pl
    real*8 qst,qsat,dqsatdt,dq,slh
    integer n
    slh=lhx*bysha
    do n=1,3
      qst=qsat(t,lhx,pl)
      dq=(q-qst)/(1.+slh*qst*dqsatdt(t,lhx))
      t=t+slh*dq
      q=q-dq
    end do
    return
  end subroutine moist_adiabat_tq

end subroutine qmom_topo_adjustments

#ifdef CACHED_SUBDD
  subroutine cijh_defs(arr,nmax,decl_count)
  !
  ! 2D outputs
  !
  use subdd_mod, only : info_type
  ! info_type_ is a homemade structure constructor for older compilers
  use subdd_mod, only : info_type_
  implicit none
  integer :: nmax,decl_count
  type(info_type) :: arr(nmax)
  !
  ! note: next() is a locally declared function to increment decl_count
  !

  decl_count = 0

  arr(next()) = info_type_( &
    sname = 'isccp_sunlit', &
    lname = 'Flag for Day (1) or Night (0)', &
    units = '-' &
       )

  arr(next()) = info_type_( &
    sname = 'isccp_ctp', &
    lname = 'Mean Cloud Top Pressure', &
    units = 'mb' &
       )

  arr(next()) = info_type_( &
    sname = 'isccp_tau', &
    lname = 'Mean Optical Thickness', &
    units = '-' &
       )

 arr(next()) = info_type_( &
    sname = 'isccp_lcld', &
    lname = 'Low Cloud Fraction', &
    units = 'fraction' &
       )

 arr(next()) = info_type_( &
    sname = 'isccp_mcld', &
    lname = 'Mid Level Cloud Fraction', &
    units = 'fraction' &
       )

 arr(next()) = info_type_( &
    sname = 'isccp_hcld', &
    lname = 'High Cloud Fraction', &
    units = 'fraction' &
       )
#ifdef CFMIP3_SUBDD  /* CFMIP3_SUBDD */
  arr(next()) = info_type_(                                            &
    sname = 'ctp_mc',                                                  &
    lname = 'Convective cloud top pressure',                           &
    units = 'Pa',                                                      &
    scale = 1d2)

  arr(next()) = info_type_(                                            &
    sname = 'cbp_mc',                                                  &
    lname = 'Convective cloud base pressure',                          &
    units = 'Pa',                                                      &
    scale = 1d2)

  arr(next()) = info_type_(                                            &
     sname = 'cnvfrq',                                                 &
     lname = 'Convective cloud occurrence',                            &
     units = '%',                                                      &
     scale = 1d2 )

  arr(next()) = info_type_(                                            &
     sname = 'dcnvfrq',                                                &
     lname = 'Deep convective cloud occurrence',                       &
     units = '%',                                                      &
     scale = 1d2 )

  arr(next()) = info_type_(                                            &
    sname = 'scnvfrq',                                                 &
    lname = 'Shallow convective cloud occurrence',                     &
    units = '%',                                                       &
    scale = 1d2)

  arr(next()) = info_type_(                                            &
    sname = 'mc_twp',                                                  &
    lname = 'Convective total water path',                             &
    units = 'kg/m^2')

  ! arr(next()) = info_type_(                                            &
  !   sname = 'mc_lwp',                                                  &
  !   lname = 'Convective liquid water path',                            &
  !   units = 'kg/m^2')

  arr(next()) = info_type_(                                            &
    sname = 'qatm',                                                    &
    lname = 'Atmospheric water vapour column',                         &
    units = 'kg/m^2')
#endif  /* CFMIP3_SUBDD */
      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine cijh_defs

  subroutine cijlh_defs(arr,nmax,decl_count)
  !
  ! 3D outputs
  !
  use subdd_mod, only : info_type
  ! info_type_ is a homemade structure constructor for older compilers
  use subdd_mod, only : info_type_
  use constant, only: kapa
  use TimeConstants_mod, only: SECONDS_PER_DAY
  use MODEL_COM, only : dtsrc
  implicit none
  integer :: nmax,decl_count
  type(info_type) :: arr(nmax)
  !
  ! note: next() is a locally declared function to increment decl_count
  !
  decl_count = 0

  arr(next()) = info_type_( &
    sname = 'tau_ss', &
    lname = 'optical depth of stratiform cloud', &
    units = '-' &
       )

  arr(next()) = info_type_( &
    sname = 'tau_mc', &
    lname = 'optical depth of convective cloud', &
    units = '-' &
       )

  arr(next()) = info_type_( &
    sname = 'dq_mc', &
    lname = 'moisture tendency from moist convection', &
    units = 'kg/kg/day', &
    scale = SECONDS_PER_DAY/dtsrc &
       )

  arr(next()) = info_type_( &
    sname = 'dth_mc', &
    lname = 'theta tendency from moist convection', &
    units = 'K/day', &
    scale = 1000.**kapa*SECONDS_PER_DAY/dtsrc &
       )

  arr(next()) = info_type_( &
    sname = 'dq_ss', &
    lname = 'moisture tendency from stratiform cloud', &
    units = 'kg/kg/day', &
    scale = SECONDS_PER_DAY/dtsrc &
       )

  arr(next()) = info_type_( &
    sname = 'dth_ss', &
    lname = 'theta tendency from stratiform cloud', &
    units = 'K/day', &
    scale = 1000.**kapa*SECONDS_PER_DAY/dtsrc &
       )

  return
  contains
  integer function next()
  decl_count = decl_count + 1
  next = decl_count
  end function next
  end subroutine cijlh_defs
#endif

