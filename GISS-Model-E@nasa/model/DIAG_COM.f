#include "rundeck_opts.h"

      MODULE DIAG_COM
!@sum  DIAG_COM Diagnostic model variables
!@auth Original Development Team
!@ver  2010/11/12
      use resolution, only : im,jm,lm,ls1=>ls1_nominal
      USE ATM_COM, only : lm_req
#ifndef SCM
      use diag_zonal, only : jm_budg
#endif
      use socpbl, only : npbl=>n
      use cdl_mod
#ifndef CUBED_SPHERE
#ifndef SCM
      use geom, only : imh,fim,byim
#endif
#endif
      use mdiag_com, only : sname_strlen,units_strlen,lname_strlen
      use mdiag_com, only : ia_cpl
      IMPLICIT NONE
      SAVE
      private

      public LM_REQ,im,jm,lm,jm_budg

#ifdef SCM
      integer, parameter :: jm_budg=1
#endif

#ifndef CUBED_SPHERE
#ifndef SCM
      public :: imh,fim,byim
#endif
#endif

!@var LSTR level of interface between low and mid strat. (approx 10 mb)
      INTEGER, public :: LSTR = LM   ! defaults to model top.

!!  WARNING: if new diagnostics are added, change io_diags/reset_DIAG !!
C**** ACCUMULATING DIAGNOSTIC ARRAYS

!@var LAT_BUDG latitudes of budget grid
!@var DXYP_BUDG area array of budget grid
!@var J_BUDG a mapping array that takes every grid point to the 
!@+   zonal mean budget array
!@var j_0b, j_1b are the min/max zonal budget latitudes for this processor
!@var wtbudg,wtbudg2 area weights for diagnostics on budget grid

      REAL*8, DIMENSION(JM_BUDG), public ::
     &     LAT_BUDG,DXYP_BUDG,DXYP_BUDG_LOC

      integer, public :: j_0b=1, j_1b=1

#ifdef SCM
      integer, public :: J_BUDG(1,1)=1
      real*8, public, dimension(1,1) :: WTBUDG=1d0,WTBUDG2=1d0
#else
      integer, public, allocatable, dimension(:,:) :: J_BUDG
      real*8, public, allocatable, dimension(:,:) :: WTBUDG,WTBUDG2
#endif

!@param KAJ number of accumulated zonal budget diagnostics
      INTEGER, PARAMETER, public :: KAJ=85
#ifdef HEALY_LM_DIAGS
     &                                 + 3
      REAL*8, public :: CROPS_DIAG(IM,JM)
#endif

!@var AJ zonal budget diagnostics for each surface type
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AJ,AJ_loc
     &     ,AJ_out

C**** Define surface types (mostly used for weighting AJ diagnostics)
!@param NTYPE number of different surface types
      INTEGER, PARAMETER, public :: NTYPE=6   ! orig = 3
!@var FTYPE fractions of each surface type
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public   :: FTYPE

!@param ITxx indices of various types (used only when it matters)
      INTEGER, PARAMETER, public :: ITOCEAN=1, ITOICE=2, ITEARTH=3,
     *                              ITLANDI=4, ITLAKE=5, ITLKICE=6

!@var SQRTM moved from DIAG5A where it was a saved local array to this
!@var place so its size could be allocated dynamically and still have
!@var it preserved from call to call of DIAG5A
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: SQRTM

!@param NREG number of regions for budget diagnostics
      INTEGER, PARAMETER, public :: NREG=24
!@var AREG regional budget diagnostics
      REAL*8, DIMENSION(NREG,KAJ), public :: AREG,AREG_loc,
     &     AREG_out

!@var TITREG,NAMREG title and names of regions for AREG diagnostics
      CHARACTER*4, public :: TITREG*80,NAMREG(2,23)
!@var JREG lat/lon array defining regions for AREG diagnostics
      INTEGER, ALLOCATABLE, DIMENSION(:,:), public :: JREG
cmax      INTEGER, DIMENSION(IM,JM), public :: JREG
!@var SAREA_REG areas of the special regions
      REAL*8, DIMENSION(NREG), public :: SAREA_REG
!@var write_regions whether to write regional diags to acc files
      logical, public :: write_regions

!@param KAJL number of AJL diagnostics
      INTEGER, PARAMETER, public :: KAJL=81
!@var AJL latitude/height diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AJL,AJL_loc

!@param KASJL number of ASJL diagnostics
      INTEGER, PARAMETER, public :: KASJL=5
!@var ASJL latitude/height supplementary diagnostics (merge with AJL?)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: ASJL,ASJL_loc

!@param KAIJ number of AIJ diagnostics
      INTEGER, PARAMETER, public :: KAIJ=750
#ifdef ACCMIP_LIKE_DIAGS
     &                                   + 8
#endif
#ifdef HEALY_LM_DIAGS
     &                                   + 1
#endif
#ifdef ENT_DEBUG_DIAGS
     &                                  + 5+256+16*9+2
#endif

!@param KAIJmm maximum number of AIJ min/max diagnostics
      INTEGER, PARAMETER, public :: KAIJmm=10

!@var AIJ latitude/longitude diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AIJ,AIJ_loc,
     &     AIJmm

!@param KAIJL number of AIJL accumulations
      INTEGER, PARAMETER, public :: KAIJL=24
#if (defined mjo_subdd) || (defined etc_subdd)
     &                                  + 8
#endif
#ifdef CLD_AER_CDNC
     &                                  + 16
#endif
#ifdef ACCMIP_LIKE_DIAGS
     &                                  +  1
#endif
#ifdef AIJL_CP_TRANSPORTS
     &                                  +  8
#endif
!@var IJL_xxx,IJK_xxx AIJL diagnostic indices
!@+   IJL/IJK refer to model versus constant-pressure levels
      INTEGER, public ::
     &     IJL_DP,IJK_DP,IJL_U,IJL_V,IJK_TX,IJK_Q,
     &     IJL_W,IJK_RH,IJL_MC
      INTEGER, public :: IJL_MCamFX, IJL_cldwtr,IJL_cldice
     &    ,IJL_LLH,IJL_MCTLH,IJL_MCDLH,IJL_MCSLH
     &    ,IJL_LDRY,IJL_TMCDRY,IJL_DMCDRY,IJL_SMCDRY
     &    ,IJL_REWM,IJL_REWS,IJL_CDWM,IJL_CDWS,IJL_CWWM,IJL_CWWS
     &    ,IJL_REIM,IJL_REIS,IJL_CDIM,IJL_CDIS,IJL_CWIM,IJL_CWIS
     &    ,IJL_CFWM,IJL_CFIM,IJL_CFWS,IJL_CFIS
     &    ,IJL_TEMPL,IJL_GRIDH,IJL_HUSL,IJL_ZL,IJL_CDTOMAS
     &    ,ijl_airmass=0
#ifdef AIJL_CP_TRANSPORTS
      INTEGER, public :: ijk_ucp,ijk_vcp,
     &     ijk_utcp,ijk_vtcp,ijk_uqcp,ijk_vqcp,ijk_uphicp,ijk_vphicp
#endif

!@var AIJL 3D accumulations for longitude/latitude/level diagnostics
      REAL*8, DIMENSION(:,:,:,:), allocatable, public :: AIJL,AIJL_loc

!@var NPTS number of points at which standard conserv. diags are called
      INTEGER, PARAMETER, public :: NPTS = 11
!@param NQUANT Number of conserved quantities in conservation diags
      INTEGER, PARAMETER, public :: NQUANT=24
!@param KCON number of conservation diagnostics
      INTEGER, PARAMETER, public :: KCON=170
!@var CONSRV conservation diagnostics
      REAL*8, DIMENSION(JM_BUDG,KCON), public :: CONSRV_loc,CONSRV
!@var SCALE_CON scales for conservation diagnostics
      REAL*8, DIMENSION(KCON), public :: SCALE_CON
!@var TITLE_CON titles for conservation diagnostics
      CHARACTER*32, DIMENSION(KCON), public :: TITLE_CON
!@var NSUM_CON indices for summation of conservation diagnostics
!@var IA_CON IDACC numbers for conservation diagnostics
      INTEGER, DIMENSION(KCON), public :: NSUM_CON, IA_CON
!@var NOFM indices for CONSRV array
      INTEGER, DIMENSION(NPTS+1,NQUANT), public :: NOFM
!@var icon_xx indexes for conservation quantities
      INTEGER, public ::
     &     icon_AM,icon_KE,icon_MS,icon_TPE,icon_WM,icon_LKM,
     *     icon_LKE,icon_EWM,icon_WTG,icon_HTG,icon_OCE,icon_OMSI,
     *     icon_OHSI,icon_OSSI,icon_LMSI,icon_LHSI,icon_MLI,icon_HLI,
     *     icon_MICB,icon_HICB
!@var KCMX actual number of conservation diagnostics
      INTEGER, public :: KCMX = 25 ! take up first 25 indexes for special cases
!@var CONPT0 default titles for each point where conserv diags. are done
      CHARACTER*10, DIMENSION(NPTS), public :: CONPT0 = (/
     *     "DYNAMICS  ","CONDENSATN","RADIATION ","PRECIPITAT",
     *     "LAND SURFC","SURFACE   ","FILTER    ","OCEAN     ",
     *     "DAILY     ","SRF OCN FL","OCN DYNAM "/)

!@param HR_IN_DAY hours in day
      INTEGER, PARAMETER, public :: HR_IN_DAY=24
!@param lmax_dd2 most upper layer for which multilayer diurnal diagnostics
!+               is written, currently: up to first constant pressure layer
      integer, parameter, public :: lmax_dd2=ls1
!@param NDIUVAR number of diurnal diagnostics
#ifdef TRACERS_AMP
c      INTEGER, PARAMETER, public :: NDIUVAR=73+16+16+100+40+40+40+40
      INTEGER, PARAMETER, public :: NDIUVAR=700
#else
#ifdef TRACERS_DUST
      INTEGER, PARAMETER, public :: NDIUVAR=74+14*lmax_dd2+6*npbl
     &     +4*(npbl-1)
#else
#if (defined TRACERS_MINERALS)
      INTEGER, PARAMETER, public :: NDIUVAR=63
#else
      INTEGER, PARAMETER, public :: NDIUVAR=60
#endif
#endif
#endif
!@param NDIUPT number of points where diurnal diagnostics are kept
#ifdef SCM
      INTEGER, PARAMETER, public :: NDIUPT=1
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      INTEGER, PARAMETER, public :: NDIUPT=34
#else
      INTEGER, PARAMETER, public :: NDIUPT=4
#endif
#endif
!@dbparam adiurn_dust  flag to switch on/off intra daily diagnostics for dust
!@+                    default=0 (off)
      INTEGER, public :: adiurn_dust=0
!@dbparam save3dAOD flag to calculate ttasv_save in rad code even when
!@+ adiurn_dust is off.
      INTEGER, public :: save3dAOD=0
!@dbparam IJDD,NAMDD (i,j)-coord.,names of boxes w/diurnal cycle diag
!@+       defaults set in DIAG_RES (depends on resolution)
      INTEGER, DIMENSION(2,NDIUPT), public :: IJDD
      CHARACTER*4, DIMENSION(NDIUPT), public :: NAMDD
!@dbparam LLDD (lon,lat)-coords (deg) of boxes w/diurnal cycle diag
!@+       defaults set in init_DIAG
      REAL*8, DIMENSION(2,NDIUPT), public :: LLDD

!@var ADIURN diurnal diagnostics (24 hour cycles at selected points)
      REAL*8, DIMENSION(NDIUVAR,NDIUPT,HR_IN_DAY), public :: ADIURN
     &     ,ADIURN_loc
!@param HR_IN_MONTH max hours in month
      INTEGER, public :: HR_IN_MONTH
#ifdef USE_HDIURN
!@var HDIURN hourly diagnostics (hourly value at selected points)
!@+     Same quantities as ADIURN but not averaged over the month
      REAL*8, allocatable, public :: HDIURN(:,:,:), HDIURN_loc(:,:,:)
#endif

!@param KAIJK number of lat/lon constant pressure diagnostics
      INTEGER, PARAMETER, public :: KAIJK=15
!@var AIJK lat/lon constant pressure diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:), public :: AIJK,AIJK_loc

C**** parameters and variables for ISCCP diags
!@param ntau,npress number of ISCCP optical depth,pressure categories
      integer, parameter, public :: ntau=7,npres=7
!@param nisccp number of ISCCP histogram regions
      integer, parameter, public :: nisccp = 5
!@var isccp_press pressure mid points for isccp histogram
      INTEGER, PARAMETER, public ::
     &     isccp_press(npres) = (/ 90, 245, 375, 500,
     *     630, 740, 900 /)
!@var isccp_tau lower bound of optical depth for each isccp tau category
      REAL*8, PARAMETER, public ::
     &     isccp_tau(ntau) = (/ 0d0,.1d0,1.3d0,3.6d0,
     *     9.4d0,23d0,60d0 /)
!@var isccp_taum mid point of optical depth for each isccp tau category
      REAL*8, PARAMETER, public :: isccp_taum(ntau-1) = (/ .05d0,0.7d0,2
     *     .95d0,6.5d0,16.2d0,41.5d0 /)
!@var isccp_late edge latitudes for each isccp lat category (region)
!@var isccp_lat midpoint latitudes for each isccp lat category (region)
! calculation of midpoints is hard-coded until all fortran compilers
! allow array arithmetic in PARAMETERS
      REAL*8, PARAMETER, public ::
     &  isccp_late(nisccp+1)=(/-60,-30,-15,15,30,60/)
     & ,isccp_lat(nisccp)=(/-45.,-22.5,0.,22.5,45./)
!@var AISCCP accumlated array of ISCCP histogram
      real*8, public, dimension(ntau,npres,nisccp) :: AISCCP,AISCCP_loc
!@var WISCCP denominator array for ISCCP histograms
      real*8, public, dimension(nisccp) :: WISCCP

!@param KGZ_max maximum number of pressure levels for some diags
      INTEGER, PARAMETER, public :: KGZ_max = 21
!@param kgz is the actual number of geopotential heights saved
!@+     which will be less than kgz_max for low model tops
      INTEGER, public :: kgz
!@param PMB pressure levels for geopotential heights (extends to strat)
!@param GHT ~mean geopotential heights at PMB level (extends to strat)
!@param PMNAME strings describing PMB pressure levels
      REAL*8, DIMENSION(KGZ_max), PARAMETER, public :: 
     *   GHT = (/     0d0,  900d0, 1500d0, 3000d0, 4500d0, 5600d0,
     *             7800d0, 9500d0,11000d0,12500d0,14500d0,16400d0,
     *            18500d0,20000d0,24000d0,27000d0,30000d0,35000d0,
     *            47000d0,53000d0,64000d0 /),
     *   PMB = (/  1000d0,  925d0,  850d0,  700d0,  600d0,  500d0,
     *              400d0,  300d0,  250d0,  200d0,  150d0,  100d0,
     *               70d0,   50d0,   30d0,   20d0,   10d0,    5d0,
     *                1d0,   .5d0,   .1d0 /)
      CHARACTER*4, DIMENSION(KGZ_max), PARAMETER, public ::
     *  PMNAME= (/ '1000', '925 ', '850 ', '700 ', '600 ', '500 ',
     *             '400 ', '300 ', '250 ', '200 ', '150 ', '100 ',
     *             '70  ', '50  ', '30  ', '20  ', '10  ', '5   ',
     *             '1   ', 'p5  ', 'p1  ' /)
#ifdef TRACERS_SPECIAL_Shindell
!@var O_inst saved instantaneous Ox tracer (at PMB lvls)
!@var X_inst saved instantaneous NOx tracer (at PMB lvls)
!@var N_inst saved instantaneous NO2 non-tracer (at PMB lvls)
!@var M_inst saved instantaneous CO tracer (at PMB lvls)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     & O_inst,N_inst,X_inst,M_inst
#endif
#ifdef TES_LIKE_DIAGS
!@param kgz_max_more is the actual number of TES pressure levels saved
      INTEGER, public :: kgz_max_more
!@param KGZmore more than KGZ number of pressure levels for some diags
      INTEGER, PARAMETER, public :: KGZmore = 32
!@param PMBmore like PMB but more preesure levels for some diags
      REAL*8, DIMENSION(KGZmore), PARAMETER, public ::
     &     PMBmore=(/
     & 1000.0000, 825.40198, 681.29102, 562.34198, 464.16000, 383.11700,
     & 316.22699, 261.01599, 215.44400, 177.82899, 146.77901, 121.15200,
     & 100.00000, 82.540604, 68.129501, 56.233898, 46.415798, 38.311901,
     & 31.622900, 26.101700, 21.544300, 17.782801, 14.678000, 12.115300,
     & 10.000000, 8.2540197, 5.1089802, 3.1622701, 2.1544299, 1.3335201,
     &  0.681292, 0.2154430/)
!@param PMNAMEmore strings describing PMBmore pressure levels
      CHARACTER*4, DIMENSION(KGZmore),PARAMETER,public::PMNAMEmore=(/
     & "1000", "825 ", "681 ", "562 ", "464 ", "383 ",
     & "316 ", "261 ", "215 ", "178 ", "147 ", "121 ",
     & "100 ", "82.5", "68.1", "56.2", "46.4", "38.3",
     & "31.6", "26.1", "21.5", "17.8", "14.7", "12.1",
     & "10.0", "8.25", "5.11", "3.16", "2.15", "1.33",
     & "0.68", "0.22" /)
!@var Q_more saved instantaneous specific hum (at PMBmore lvls)
!@var T_more saved instantaneous temperature(at PMBmore lvls)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: Q_more,T_more
#ifdef TRACERS_SPECIAL_Shindell
!@var O_more saved instantaneous Ox tracer (at PMBmore lvls)
!@var X_more saved instantaneous NOx tracer (at PMBmore lvls)
!@var N_more saved instantaneous NO2 non-tracer (at PMBmore lvls)
!@var M_more saved instantaneous CO tracer (at PMBmore lvls)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     & O_more,N_more,X_more,M_more
#endif
#endif /* TES_LIKE_DIAGS */
!@var save{H,M,L}CLDI,saveCTPI,saveTAUI,save{S,T}CLDI: SUBDDiag
!@+  instantaneous save arrays for ISCCP cloud variables
!@var saveMCCLDTP instnt.SUBDD moist convective cloud top pressure
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: saveHCLDI,
     & saveMCLDI,saveLCLDI,saveCTPI,saveTAUI,saveSCLDI,saveTCLDI,
     & saveMCCLDTP
C**** Instantaneous constant pressure level fields
!@var Z_inst saved instantaneous height field (at PMB levels)
!@var RH_inst saved instantaneous relative hum (at PMB levels)
!@var T_inst saved instantaneous temperature(at PMB levels)
!@var P_acc accumulated precip (special for SUBDD)
!@var PM_acc accumulated moist convective precip (special for SUBDD)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     &     Z_inst,RH_inst,T_inst
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: P_acc,PM_acc
#if (defined ttc_subdd) || (defined etc_subdd)
!@var u_inst saved instantaneous U (at PMB levels)
!@var v_inst saved instantaneous V (at PMB levels)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     &     u_inst,v_inst
#endif
#ifdef ttc_subdd
!@var vt_inst saved instantaneous Vorticity (at PMB levels)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: vt_inst
#endif
#ifdef etc_subdd
!@var omg_inst saved instantaneous omega (at PMB levels)
!@var lwc_inst saved instantaneous cloud liquid water content (at PMB levels)
!@var iwc_inst saved instantaneous cloud ice water content (at PMB levels)
!@var cldmc_inst saved instantaneous convective cloud fraction (at PMB levels)
!@var cldss_inst saved instantaneous stratiform cloud fraction (at PMB levels)
!@var tlh_inst saved instantaneous diabatic heating from moist convective (at PMB levels)
!@var dlh_inst saved instantaneous diabatic heating from deep convective (at PMB levels)
!@var slh_inst saved instantaneous diabatic heating from shallow convective (at PMB levels)
!@var llh_inst saved instantaneous diabatic heating from large-scale condensation (at PMB levels)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     &     omg_inst,lwc_inst,iwc_inst,
     &     cldmc_inst,cldss_inst,tlh_inst,llh_inst,dlh_inst,slh_inst
#endif
#if (defined  mjo_subdd) || (defined etc_subdd)
!@var qsen_avg,qlat_avg,pblht_acc accumulated surface sensible/latent heat flux, PBL height for SUBDD
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public ::
     &        pblht_acc,qlat_avg, qsen_avg
#endif
#ifdef mjo_subdd
!@var E_acc accumulated evaporation (special for SUBDD)
!@var PW_acc accumulated integrated atmospheric water content (precipitable water) (special for SUBDD)
!@var p_avg, lwu_avg,sst_avg, accumulated surface pressure, 
!@     upward LW at sfc, sst, latent heat flux for SUBDD
!@var u_avg,v_avg,w_avg,t_avg,q_avg,r_avg,z_avg accumulated u,v,w,t,q,rh,gh for subdaily period
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: E_acc,PW_acc
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public ::
     &        p_avg,sst_avg,lwu_avg
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     &        u_avg,v_avg,w_avg,t_avg,q_avg,r_avg,z_avg
#endif

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:), public :: AFLX_ST

!@param KTSF number of freezing temperature diagnostics
      integer, parameter, public :: ktsf=4
!@var TSFREZ freezing temperature diagnostics
C****   1  FIRST DAY OF GROWING SEASON (JULIAN DAY)
C****   2  LAST DAY OF GROWING SEASON (JULIAN DAY)
C****   3  LAST DAY OF ICE-FREE LAKE (JULIAN DAY)
C****   4  LAST DAY OF ICED-UP LAKE  (JULIAN DAY)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: TSFREZ,TSFREZ_loc

!@param KTD number of diurnal temperature diagnostics
      INTEGER, PARAMETER, public :: KTD=9
!@var TDIURN diurnal range temperature diagnostics
C****   1  MIN TG1 OVER EARTH FOR CURRENT DAY (C)
C****   2  MAX TG1 OVER EARTH FOR CURRENT DAY (C)
C****   3  MIN TS OVER EARTH FOR CURRENT DAY (K)
C****   4  MAX TS OVER EARTH FOR CURRENT DAY (K)
C****   5  SUM OF COMPOSITE TS OVER TIME FOR CURRENT DAY (C)
C****   6  MAX COMPOSITE TS FOR CURRENT DAY (K)
C****   7  MAX TG1 OVER OCEAN ICE FOR CURRENT DAY (C)
C****   8  MAX TG1 OVER LAND ICE FOR CURRENT DAY (C)
C****   9  MIN COMPOSITE TS FOR CURRENT DAY (K)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: TDIURN
      ! REAL*8 :: TDIURN_glob(IM, JM, KTD)
      REAL*8,allocatable, public :: TDIURN_glob(:,:,:)

!@nlparam KDIAG array of flags to control diagnostics printout
      INTEGER, DIMENSION(13), public :: KDIAG

!@param NKEYNR number of key number diagnostics
      INTEGER, PARAMETER, public :: NKEYNR=43
!@param NKEYMO number of months key diagnostics are saved
      INTEGER, PARAMETER, public :: NKEYMO=50
!@var KEYNR time-series of key numbers
      INTEGER, DIMENSION(NKEYNR,NKEYMO), public :: KEYNR = 0
!@var KEYCT next index in KEYNR to be used (1->nkeymo)
      INTEGER, public :: KEYCT = 1

!@nlparam IWRITE,JWRITE,ITWRITE control rad.debug output (i,j,amount)
      INTEGER, public :: IWRITE = 0, JWRITE = 0, ITWRITE = 23
!@nlparam QDIAG TRUE for outputting binary diagnostics
      LOGICAL, public :: QDIAG = .FALSE.
!@nlparam QDIAG_ratios TRUE for forming ratios if title="q1 x q2"
      LOGICAL, public :: QDIAG_ratios = .TRUE.

!@var OA generic diagnostic array for ocean heat transport calculations
C****
C****       DATA SAVED IN ORDER TO CALCULATE OCEAN TRANSPORTS
C****
C****       1  ACE1I+SNOWOI  (INSTANTANEOUS AT NOON GMT)
C****       2  MSI2   (INSTANTANEOUS AT NOON GMT)
C****       3  HSIT   (INSTANTANEOUS AT NOON GMT)
C****       4  ENRGP  (INTEGRATED OVER THE DAY)
C****       5  SRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       6  TRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       7  SHDT   (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       8  EVHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       9  TRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      10  SHDT   (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      11  EVHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      12  SRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****
C**** Extra array needed for dealing with advected ice
C****      13  HCHSI  (HORIZ CONV SEA ICE ENRG, INTEGRATED OVER THE DAY)
C****
!@param KOA number of diagnostics needed for ocean heat transp. calcs
      INTEGER, public :: iu_VFLXO
      INTEGER, PARAMETER, public :: KOA = 13
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: OA
      ! REAL*8 :: OA_glob    (IM, JM, KOA)
      REAL*8,allocatable, public :: OA_glob(:,:,:)

!**** Diagnostic control parameters
!@dbparam Kvflxo if 1 => vert.fluxes into ocean are saved daily
      integer, public :: Kvflxo=0

C****
C**** Information about acc-arrays:
C****      names, indices, units, idacc-numbers, etc.

!@var iparm/dparm int/double global parameters written to acc-file
      integer, parameter, public :: niparm_max=100
      character(len=20), dimension(niparm_max), public :: iparm_name
      integer, dimension(niparm_max), public :: iparm
      integer, public :: niparm=0
      integer, parameter, public :: ndparm_max=100
      character(len=20), dimension(ndparm_max), public :: dparm_name
      REAL*8, dimension(ndparm_max), public :: dparm
      integer, public :: ndparm=0

!@var J_xxx zonal J diagnostic names
      INTEGER, public ::
     &     J_SRABS,
     *     J_TRHDT, J_RNFP0, J_RNFP1,
     *     J_RHDT, J_SHDT, J_EVHDT, J_HZ1, J_TG2, J_TG1, J_EVAP,
     *     J_PRCP, J_TX, J_TX1, J_TSRF, J_DTSGST, J_DTDGTR, J_RICST,
     *     J_RICTR, J_ROSST, J_ROSTR, J_RSI, J_TYPE, J_RSNOW,
     *     J_OHT, J_DTDJS, J_DTDJT, J_LSTR, J_LTRO, J_EPRCP,
     *     J_RUN, J_ERUN, J_HZ0, J_LWCORR,
     *     J_RVRD,J_ERVR,J_IMPLM, J_IMPLH,
     *     J_WTR1,J_ACE1, J_WTR2,J_ACE2, J_SNOW, J_HZ2,
     *     J_CTOPP, J_PRCPSS, J_PRCPMC, J_QP,
     *     J_GAM,J_GAMM, J_GAMC,J_FTHERM,
     *     J_ALBP0,J_ALBG,J_IRGW,J_IRGWE
#ifdef HEALY_LM_DIAGS
     *     ,J_CROPS
#endif
!@var NAME_J,UNITS_J Names/Units of zonal J diagnostics
      character(len=sname_strlen), dimension(kaj), public :: name_j
      character(len=units_strlen), dimension(kaj), public :: units_j
!@var LNAME_J Long names of zonal J diagnostics
      character(len=lname_strlen), dimension(kaj), public :: lname_j
!@var STITLE_J short titles for print out for zonal J diagnostics
      character(len=16), dimension(kaj), public :: stitle_j
!@var SCALE_J scale for zonal J diagnostics
      real*8, dimension(kaj), public :: scale_j
!@var IA_J IDACC indexes for zonal J diagnostics
      integer, dimension(kaj), public :: ia_j
!@var FMT_J Format strings for zonal J diagnostics
      character(len=30), dimension(kaj), public :: fmt_j,fmt_reg
!@var iden_j denominators for zonal J diagnostics
      integer, dimension(kaj), public :: iden_j,iden_reg
!@var HEMIS_J hemispheric/global averages of AJ
      real*8, dimension(:,:,:), allocatable, public :: hemis_j

      character(len=sname_strlen), dimension(kaj), public :: name_reg

!****
!@var IJ_xxx AIJ diagnostic names
!****
      INTEGER, public ::
!**** Vertical Energy Fluxes
     &  IJ_dSE_Dyn,IJ_dKE_Dyn,IJ_dTE_Dyn,IJ_dHSI_Dyn,
     &     IJ_RSOI, IJ_RSNW, IJ_SNOW, IJ_SHDT, IJ_PREC, IJ_EVAP,
     *  IJ_SSAT, IJ_BETA,  IJ_SLP1,  IJ_P4UV, IJ_PRES,
     *  IJ_PMB1,IJ_TPMB1,IJ_QPMB1,IJ_ZPMB1,IJ_RHPMB1,IJ_UPMB1,IJ_VPMB1,
     *  ij_OMEGAPMB1,
     *     IJ_PBLHT, IJ_DSEV,
     *     IJ_RH1,
     *     IJ_SRTR, IJ_NETH,
     *     IJ_TG1, IJ_RSIT, IJ_TDSL, IJ_TDCOMP, IJ_DTDP,
     *     IJ_RUNE, IJ_TS1, IJ_RUNLI, IJ_WS, IJ_TS, IJ_US, IJ_VS,
     *     IJ_SLP, IJ_UJET, IJ_VJET,
     *     IJ_TOC2, IJ_TAUS, IJ_TAUUS,
     *     IJ_TAUVS, IJ_GWTR, IJ_QS, IJ_STRNGTS, IJ_ARUNU, IJ_DTGDTS,
     *     IJ_PUQ, IJ_PVQ, IJ_TGO, IJ_TGO2, IJ_EVAPO, ij_RHs,
     *     IJ_EVAPI, IJ_EVAPLI,IJ_EVAPE, IJ_F0OC,IJ_F0LI,IJ_F0E,
     *     IJ_F1LI, IJ_SNWF, IJ_TSLI, IJ_SHDTLI, IJ_EVHDT,
     *     IJ_TRHDT, IJ_TMAXE, IJ_TMAXC, IJ_TMINC, IJ_TMNMX, IJ_PEVAP,
     *     IJ_WMSUM, IJ_PSCLD, IJ_PDCLD, IJ_DCNVFRQ, IJ_SCNVFRQ,
     *     IJ_CNVFRQ,
     *     IJ_EMTMOM, IJ_SMTMOM, IJ_FMU, IJ_FMV, IJ_SSTABX,
     *     IJ_FGZU, IJ_FGZV, IJ_ERVR, IJ_MRVR, IJ_SSS, IJ_PRECMC,
     *     IJ_LKON, IJ_LKOFF, IJ_LKICE, IJ_PTROP, IJ_TTROP,
     *     IJ_P850,
     *     IJ_GPP, IJ_IPP, IJ_RAUTO, IJ_CLAB, IJ_DLEAF, IJ_LAI, !VEG DIAGNOSTICS
     *     IJ_SOILRESP, IJ_SOILCPOOLSUM, IJ_ECVF, IJ_LANDCARBON,
     *     IJ_GICE, IJ_GWTR1, IJ_ZSNOW, IJ_AFLMLT, IJ_AERUNS, IJ_AERUNU,
     *     IJ_HTSOIL, IJ_HTSNOW, IJ_AINTRCP, IJ_MCCLDTP, IJ_MCCLDBS,
     *     IJ_TRSDN,
     *     ij_tclssct,ij_rclssct,
#ifdef CLD_AER_CDNC
     *     ij_nclssct,
#endif
#ifdef TRACERS_AMP
#ifdef BLK_2MOM
     *     ij_ccnssct,
#endif
#endif
     *     IJ_TRSUP, IJ_CLDW, IJ_CLDI, IJ_QM, IJ_SSH, IJ_FWOC,
     *     IJ_DSKIN, IJ_DSKINSNOW, IJ_MCCVTP, IJ_MCCVBS,
     *     IJ_LI, IJ_LK,
     &     IJ_FVEG,IJ_GUSTI, IJ_MCCON
     *     ,IJ_WISUM, IJ_SLPQ, IJ_PRESQ
     *     ,ij_dzwm,ij_dzim,ij_dzws,ij_dzis
     *     ,ij_3dnwm,ij_3dnim,ij_3dnws,ij_3dnis
     *     ,ij_3drwm,ij_3drim,ij_3drws,ij_3dris
     *     ,ij_3dlwm,ij_3dlim,ij_3dlws,ij_3dlis
     *     ,ij_ssprec,ij_mcprec,IJ_WMCLWP,IJ_WMCTWP
     &     ,ij_wdry,ij_wtke,ij_wmoist,ij_wsgcm,ij_wspdf
     &     ,ij_flam,ij_CtoG,ij_flash
     *     ,ij_fvden,ij_fireC
     *     ,ij_flam_rh,ij_flam_prec,ij_flam_tsurf
     *     ,ij_nsuppress,ij_cgign,ij_humanign
     *     ,ij_barh1,ij_bawsurf
     *     ,ij_a_tree,ij_a_shrub,ij_a_grass
     *     ,ij_ba_tree,ij_ba_shrub,ij_ba_grass
     *     ,ij_swaerabs,ij_lwaerabs
     *     ,ij_swaerabsnt
     *     ,ij_lwaerabsnt,ij_evapsn,ij_irrW,ij_irrE,ij_irrW_tot
     *     ,ij_mwl,ij_gml,ij_mwlir,ij_gmlir,ij_irrgw,ij_irrgwE
     *     ,ij_rvrflo ,ij_sisnd
     *     ,ij_silwd,ij_silwu,ij_sish
     *     ,ij_impmli,ij_imphli,ij_eicb,ij_micb, IJ_ERVRO, IJ_MRVRO
     *     ,IJ_IMPMGR,IJ_IMPHGR,IJ_IMPMKI,IJ_IMPHKI
     *     ,IJ_MLKtoGR,IJ_HLKtoGR
     *     ,ij_precli,ij_precsi,ij_precoo,ij_precgr
     &     ,ij_ent_debug
#ifdef HEALY_LM_DIAGS
     &     ,IJ_CROPS
#endif
      integer, dimension(6,8), public :: 
     &  ij_nintaerext,ij_nintaersca,ij_nintaerasy
      integer, dimension(:), allocatable, public ::
     &      ij_kw, ij_alpha, ij_gasx
!@var IJ_Gxx names for old AIJG arrays
      INTEGER, public ::
     *     ij_gbsw, ij_gbsbet, ij_gbetpen, ij_gvsw, ij_gbvswt,
     *     ij_gconatm, ij_gconcan, ij_gevppen, ij_gbst, ij_gbsevp,
     *     ij_gdcevp, ij_gwcevp, ij_gvst, ij_gwtbl, ij_gvswet,
     *     ij_gbetat, ij_gbssnd, ij_gvssnd
!@var IJ_GWx names for gravity wave diagnostics
      INTEGER, public ::
     &     IJ_GW1,IJ_GW2,IJ_GW3,IJ_GW4,IJ_GW5,IJ_GW6,IJ_GW7,IJ_GW8
     *     ,IJ_GW9
!@var IJ_xxxI names for ISCCP diagnostics
      INTEGER, public ::
     &     IJ_CTPI,IJ_TAUI,IJ_LCLDI,IJ_MCLDI,IJ_HCLDI,IJ_TCLDI,IJ_SCLDI
c weighting fractions
      INTEGER, public :: IJ_PSOIL,IJ_CLRSKY,IJ_POCEAN,IJ_POPOCN,IJ_VSFR
     *     ,IJ_BSFR,IJ_POPWAT,IJ_PWATER
c derived/composite diagnostics
      INTEGER, public ::
     *  ij_topo, ij_jet, ij_wsmn, ij_jetdir, ij_wsdir, ij_grow,
     *  ij_netrdp, ij_albp, ij_albg, ij_albv, ij_ntdsese, ij_ntdsete,
     *  ij_fland, ij_dzt1, ij_albgv, ij_msutlt,ij_msutmt,ij_msutls,
     *  ij_ssu1, ij_ssu2, ij_ssu3,
     *  ij_Tatm, ij_LOTI, ij_RTSE, ij_HWV, ij_PVS

      integer, public :: ij_tsurfmin,ij_tsurfmax

!@param LEGEND "contour levels" for ij-maps
      CHARACTER(LEN=40), DIMENSION(25), PARAMETER, public :: LEGEND=(/ !
     1  '0=0,1=5...9=45,A=50...K=100             ', ! ir_pct    fac=.2
     2  '0=0...9=90,A=100...I=180...R=270        ', ! ir_angl       .1
     3  '1=.5...9=4.5,A=5...Z=17.5,+=MORE        ', ! ir_0_18        2
     4  '1=.1...9=.9,A=1...Z=3.5,+=MORE          ', ! ir_0_4        10
     5  '1=2...9=18,A=20...Z=70,+=MORE           ', ! ir_0_71       .5
     6  '1=50...9=450,A=500...Z=1750,+=MORE      ', ! ir_0_1775     .02
     7  '1=100...9=900,A=1000...Z=3500,+=MORE    ', ! ir_0_3550     .01
     8  '1=20...9=180,A=200...Z=700,+=MORE       ', ! ir_0_710      .05
     9  'A=1...Z=26,3=30...9=90,+=100-150,*=MORE ', ! ir_0_26_150    1
     O  '0=0,A=.1...Z=2.6,3=3...9=9,+=10-15      ', ! ir_0_3_15     10
     1  '-=LESS,Z=-78...0=0...9=27,+=MORE        ', ! ir_m80_28     .33
     2  '-=LESS,Z=-260...0=0...9=90,+=MORE       ', ! ir_m265_95    .1
     3  '-=LESS,Z=-520...0=0...9=180,+=MORE      ', ! ir_m530_190   .05
     4  '-=LESS,Z=-1300...0=0...9=450,+=MORE     ', ! ir_m1325_475  .02
     5  '-=LESS,Z=-2600...0=0...9=900,+=MORE     ', ! ir_m2650_950  .01
     6  '-=LESS,Z=-3900...0=0...9=1350,+=MORE    ', ! ir_m3975_1425 .007
     7  '-=LESS,Z=-5200...0=0...9=1800,+=MORE    ', ! ir_m5300_1900 .005
     8  '-=LESS,9=-.9...0=0,A=.1...Z=2.6,+=MORE  ', ! ir_m1_3       10
     9  '-=LESS,9=-45...0=0,A=5...I=45...+=MORE  ', ! ir_m45_130    .2
     O  '-=LESS,9=-90...0=0,A=10...Z=260,+=MORE  ', ! ir_m95_265    .1
     1  '-=LESS,9=-180...A=20...Z=520,+=MORE     ', ! ir_m190_530   .05
     2  '-=LESS,9=-9...0=0,A=1...Z=26,+=MORE     ', ! ir_m9_26       1
     3  '-=LESS,9=-36...0=0,A=4...Z=104,+=MORE   ', ! ir_m38_106    .25
     4  '1=5...9=45,A=50...Z=175,+=MORE          ', ! ir_0_180      .2
     5  '9=-512...1=-2,0=0,A=2,B=4,C=8...+=MORE  '/)! ir_log2       1.
!@var ir_xxxx names for indices to LEGEND indicating the (rounded) range
      integer, parameter, public ::
     &     ir_pct=1, ir_angl=2, ir_0_18=3, ir_0_4=4,
     * ir_0_71=5, ir_0_1775=6, ir_0_3550=7, ir_0_710=8, ir_0_26_150=9,
     * ir_0_3_15=10, ir_m80_28=11, ir_m265_95=12, ir_m530_190=13,
     * ir_m1325_475=14, ir_m2650_950=15, ir_m3975_1425=16,
     * ir_m5300_1900=17, ir_m1_3=18, ir_m45_130=19, ir_m95_265=20,
     * ir_m190_530=21, ir_m9_26=22, ir_m38_106=23, ir_0_180=24,
     * ir_log2=25
!@var fac_legnd = 1/(range_of_1_colorbox)
      real*8, dimension(25), public :: fac_legnd=(/
     1      1d0/5,  1d0/10,    2.d0,   10.d0,   1d0/2,
     6     1d0/50, 1d0/100,  1d0/20,    1.d0,   10.d0,
     1      1d0/3,  1d0/10,  1d0/20,  1d0/50, 1d0/100,
     6    1d0/150, 1d0/200,   10.d0,   1d0/5,  1d0/10,
     1     1d0/20,    1.d0,   1d0/4,   1d0/5,     1d0  /)

!@param CBAR "color bars" for ij-maps
      CHARACTER(LEN=38), PARAMETER, DIMENSION(5), public :: CBAR=(/
     &     ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_pos
     &     ' 0123456789ABCDEFGHIJKX               ',  ! ib_pct
     &     '-9876543210ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_npp,ib_ntr
     &     ' 0ABCDEFGHIJKLMNOPQRSTUVWXYZ3456789+* ',  ! ib_hyb
     &     '-ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789+'/) ! ib_nnp
!@var ib_xxx indices for color bars
      integer, parameter, public ::
     &     ib_pos=1,ib_pct=2,ib_npp=3,ib_hyb=4,ib_nnp=5
     &     ,ib_ntr=6

!@dbparam isccp_diags: if 1 accumulate ISCCP cloud data (default 0)
      INTEGER, public :: isccp_diags = 0

!@dbparam lh_diags: if 1 accumulate 3D drying & latent heating profiles (default 0)
      INTEGER, public :: lh_diags = 0

!@var SCALE_IJ scaling for weighted AIJ diagnostics
      REAL*8, DIMENSION(KAIJ), public :: SCALE_IJ
!@var NAME_IJ,UNITS_IJ Names/Units of lat/lon IJ diagnostics
      character(len=sname_strlen), dimension(kaij), public :: name_ij
      character(len=units_strlen), dimension(kaij), public :: units_ij
!@var LNAME_IJ Long names of lat/lon IJ diagnostics
      character(len=lname_strlen), dimension(kaij), public :: lname_ij
!@var HEMIS_IJ hemispheric/global averages of AIJ
      real*8, dimension(:,:,:), allocatable, public :: hemis_ij
!@var nwts_ij = number of weight-ij-arrays used in IJ-diagnostics
      integer, parameter, public :: nwts_ij = 8
!@var wt_ij various weight-arrays use in ij-diagnostics
      real*8, dimension(IM,JM,NWTS_IJ), public :: wt_ij
!@var IW_xxx index for weight-array
      integer, parameter, public :: iw_all=1 , iw_ocn=2 , iw_lake=3,
     *   iw_lice=4 , iw_soil=5 , iw_bare=6 , iw_veg=7, iw_land=8
!@var IR_IJ range indices for IJ diagnostics
      integer, dimension(kaij), public :: ir_ij
!@var IA_IJ IDACC indexes for lat/lon IJ diagnostics
      integer, dimension(kaij), public :: ia_ij
!@var [ij]grid_ij 1=primary grid  2=secondary grid
      integer, dimension(kaij), public :: igrid_ij,jgrid_ij
!@var denom_ij index of AIJ element to use as time/area weight
      integer, dimension(kaij), public :: denom_ij

!@var SCALE_IJmm scale factor for AIJ min/max diagnostics
      REAL*8, DIMENSION(KAIJmm), public :: SCALE_IJmm
!@var [NAME,UNITS,LNAME]_IJmm Names/Units/Longnames of
!@+     min/max IJ diagnostics
      character(len=sname_strlen), dimension(kaijmm), public ::
     &     name_ijmm
      character(len=units_strlen), dimension(kaijmm), public ::
     &     units_ijmm
      character(len=lname_strlen), dimension(kaijmm), public ::
     &     lname_ijmm

!@var JL_xxx, JK_xxx names for AJL indices
!@+   JL/JK refer to model versus constant-pressure levels
      INTEGER, public ::
     &     jl_mcmflx,jl_sshr,jl_trbhr,jl_mchr
     &     ,jl_dtdyn,jl_mcdflx
     &     ,jl_rhe,jl_damdc,jl_dammc,jl_mchphas,jl_mcdtotw
     &     ,jl_mcldht,jl_trbke,jl_trbdlht,jl_mcheat,jl_mcdry
     &     ,jl_mcdeep,jl_mcshlw,jl_cldmc,jl_cldss,jl_csizmc,jl_csizss
     &     ,jl_cnumwm,jl_cnumim,jl_cnumws,jl_cnumis
     &     ,jl_dpa,jl_dpasrc,jl_dwasrc
     &     ,jl_rad_cool
     &     ,jl_epacwt,jl_wpacwt
     &     ,jl_uepac,jl_vepac,jl_wepac,jl_uwpac,jl_vwpac,jl_wwpac
     &     ,jl_dudtsdrg,jl_dtdtsdrg
      INTEGER, public ::
     &      jl_dudfmdrg,jl_dumtndrg,jl_dushrdrg,jl_dumcdrgm10
     &     ,jl_dumcdrgp10,jl_dumcdrgm40,jl_dumcdrgp40,jl_dumcdrgm20
     &     ,jl_dumcdrgp20,jl_sdifcoef,jl_dudtsdif,jl_dudtvdif
     &     ,jl_gwfirst,jl_mcdrgpm10,jl_mcdrgpm40,jl_mcdrgpm20,jl_sumdrg

      INTEGER, public ::
     &     JK_hght, JK_dpwt, JK_tx, JK_q, JK_cldh2o ,JK_rh
     &     ,JK_cldwtr, JK_cldice


!@var JGRID_U, JGRID_KE latitudes at which U-wind and KE diags are defined
!@+   (1 for primary latitudes, 2 for secondary latitudes)
!@+   In future, can be other than the values in GEOM if requested.
      integer, public :: jgrid_u, jgrid_ke
      integer, parameter, public ::
     &     igridc=0,igride=1,jgridc=0,jgride=2,kgridc=0,kgride=4,
     &     ijkgridc=igridc+jgridc+kgridc

!@var SNAME_JL Names of lat-sigma JL diagnostics
      character(len=sname_strlen), dimension(kajl), public :: sname_jl
!@var LNAME_JL,UNITS_JL Descriptions/Units of JL diagnostics
      character(len=lname_strlen), dimension(kajl), public :: lname_jl
      character(len=units_strlen), dimension(kajl), public :: units_jl
!@var SCALE_JL printout scaling factors for JL diagnostics
      REAL*8, dimension(kajl), public :: scale_jl
!@var IA_JL,JGRID_JL,LGRID_JL idacc-numbers,gridtypes for JL diagnostics
      integer, dimension(kajl), public :: ia_jl,jgrid_jl,lgrid_jl
!@var POW_JL printed output scaled by 10**(-pow_jl)
      integer, dimension(kajl), public :: pow_jl,pow_jl_vmean
!@var DENOM_JL index of AJL element to use as weight
      integer, dimension(kajl), public :: denom_jl
!@var HEMIS_JL hemispheric/global averages of AJL
!@var VMEAN_JL vertical sums of AJL
      real*8, dimension(:,:,:), allocatable, public :: hemis_jl,vmean_jl
!@param [CTR,EDG]_[ML,CP] tags for center,edge model-layer,const-pres
!@+   vertical grids
      integer, parameter, public :: ctr_ml=1,edg_ml=2,ctr_cp=3,edg_cp=4

!@var force_jl_vmean a mechanism to force vertical averaging for qtys
!@+   at layer edges or that otherwise lack layer weighting info.
      logical, dimension(kajl), public :: force_jl_vmean

!@var NAME_SJL Names of radiative-layer-only SJL diagnostics
      character(len=sname_strlen), dimension(kasjl), public :: name_sjl
!@var LNAME_SJL,UNITS_SJL Descriptions/Units of SJL diagnostics
      character(len=lname_strlen), dimension(kasjl), public :: lname_sjl
      character(len=units_strlen), dimension(kasjl), public :: units_sjl
!@var SCALE_SJL printout scaling factors for SJL diagnostics
      REAL*8, dimension(kasjl), public :: scale_sjl
!@var IA_SJL idacc-numbers for SJL diagnostics
      integer, dimension(kasjl), public :: ia_sjl

!@var SCALE_IJK scaling for weighted AIJK diagnostics
      REAL*8, DIMENSION(Kaijk), public :: SCALE_IJK
!@var OFF_IJK offset for weighted AIJK diagnostics
      REAL*8, DIMENSION(Kaijk), public :: OFF_IJK
!@var NAME_IJK Names of lon-lat-pressure IJK diagnostics
      character(len=sname_strlen), dimension(kaijk), public :: name_ijk
!@var LNAME_IJK,UNITS_IJK Descriptions/Units of IJK diagnostics
      character(len=lname_strlen), dimension(kaijk), public ::
     &     lname_ijk
      character(len=units_strlen), dimension(kaijk), public ::
     &     units_ijk
!@var jgrid_ijk 1=primary grid  2=secondary grid
      integer, dimension(Kaijk), public :: jgrid_ijk,denom_ijk,ia_ijk

!@var SCALE_IJL scale factor for AIJL diagnostics
      REAL*8, DIMENSION(KAIJL), public :: SCALE_IJL
!@var IA_IJL,DENOM_IJL  idacc-numbers,weights for AIJL diagnostics
      INTEGER, DIMENSION(KAIJL), public :: IA_IJL,DENOM_IJL,LGRID_IJL
     &                                                     ,JGRID_IJL
!@var NAME_IJL Names of lon-lat-level IJL diagnostics
      character(len=sname_strlen), dimension(kaijl), public :: name_ijl
!@var LNAME_IJL,UNITS_IJL Descriptions/Units of IJL diagnostics
      character(len=lname_strlen), dimension(kaijl), public ::
     &     lname_ijl
      character(len=units_strlen), dimension(kaijl), public ::
     &     units_ijl

      character(len=sname_strlen), dimension(kcon), public ::
     &     name_consrv = 'unused'
      character(len=units_strlen), dimension(kcon), public ::
     &     units_consrv
      character(len=lname_strlen), dimension(kcon), public ::
     &     lname_consrv
!@var HEMIS_CONSRV hemispheric/global averages of CONSRV
      real*8, dimension(:,:), allocatable, public :: hemis_consrv

      character(len=sname_strlen), dimension(ndiuvar), public :: name_dd
      character(len=units_strlen), dimension(ndiuvar), public ::
     &     units_dd
      character(len=lname_strlen), dimension(ndiuvar), public ::
     &     lname_dd
      real*8, dimension(ndiuvar), public :: scale_dd
      integer, dimension(ndiuvar), public :: denom_dd

!@var IDD_xxx names for diurnal diagnostics
      INTEGER, public ::
c     standard set of names
     &     IDD_ECND,
     *     IDD_SPR, IDD_PT5, IDD_TS, IDD_TG1, IDD_Q5, IDD_QS,
     *     IDD_QG, IDD_SWG, IDD_LWG, IDD_SH, IDD_LH, IDD_HZ0, IDD_UG,
     *     IDD_VG, IDD_WG, IDD_US, IDD_VS, IDD_WS, IDD_CIA, IDD_RIS,
     *     IDD_RIG, IDD_CM, IDD_CH, IDD_CQ, IDD_EDS, IDD_DBL, IDD_DCF,
     *     IDD_LDC, IDD_PR, IDD_EV, IDD_DMC, IDD_SMC, IDD_W,
     *     IDD_SSP, IDD_MCP ! 56
c     names for one layer dust diagnostics
     &     ,idd_wtke,idd_wd,idd_wm,idd_wsgcm,idd_wspdf,idd_wtrsh
     &     ,idd_emis,idd_emis2,idd_ws2,idd_ustar,idd_us3,idd_stress
     &     ,idd_lmon,idd_rifl,idd_wet,idd_grav,idd_turb ! +17
c     names for llmax_dd2 layers dust diagnostics
     &     ,idd_u1,idd_v1,idd_uv1,idd_t1,idd_qq1,idd_p1,idd_w1,idd_phi1
     &     ,idd_sr1,idd_tr1,idd_load1,idd_conc1,idd_tau1,idd_tau_cs1 ! +14*ls1
c     names for npbl layers dust diagnostics
     &     ,idd_zpbl1,idd_uabl1,idd_vabl1,idd_uvabl1,idd_tabl1
     &     ,idd_qabl1           ! +6*npbl
c     names for npbl-1 layers dust diagnostics
     &     ,idd_zhat1,idd_e1,idd_km1,idd_ri1 ! +4*(npbl-1)
c    hourly AMP diagnostics
     *     ,idd_diam, idd_lwp, idd_ccn, idd_cdnc
     *     ,idd_mass, idd_numb, idd_so2, idd_lwc, idd_ncL, idd_pres

!@var tf_xxx tsfrez diagnostic names
      INTEGER, public :: tf_day1,tf_last,tf_lkon,tf_lkoff
      character(len=sname_strlen), dimension(ktsf), public :: name_tsf
      character(len=units_strlen), dimension(ktsf), public :: units_tsf
      character(len=lname_strlen), dimension(ktsf), public :: lname_tsf

      character(len=8), dimension(ntype), public :: stype_names=
     &     (/ 'OCEAN   ','OCEANICE','EARTH   ',
     &        'LANDICE ','LAKE    ','LAKEICE ' /)

!@param NTYPE_OUT number of output budgets pages
      INTEGER, PARAMETER :: NTYPE_OUT=NTYPE+3  ! to include comp/regio
C**** Expanded version of surfaces (including composites)
!@var TERRAIN name of surface type
      CHARACTER*16, DIMENSION(NTYPE_OUT), PARAMETER :: TERRAIN = (/
     *     '    (GLOBAL)','(OPEN OCEAN)',' (OCEAN ICE)','     (OCEAN)',
     *     '      (LAND)','  (LAND ICE)',' (OPEN LAKE)','  (LAKE ICE)',
     *     '     (LAKES)'/)
C**** weighting functions for surface types
      REAL*8, DIMENSION(NTYPE_OUT,NTYPE), PARAMETER ::
     *     WTJ_COMP=RESHAPE(          ! separate types + composites
     *     (/1.,1.,0.,1.,0.,0.,0.,0.,0., 1.,0.,1.,1.,0.,0.,0.,0.,0.,
     *       1.,0.,0.,0.,1.,0.,0.,0.,0., 1.,0.,0.,0.,0.,1.,0.,0.,0.,
     *       1.,0.,0.,0.,0.,0.,1.,0.,1., 1.,0.,0.,0.,0.,0.,0.,1.,1./),
     *     (/NTYPE_OUT,NTYPE/) )
      public :: ntype_out,terrain,wtj_comp

c idacc-indices of various processes
      integer, parameter, public ::
     &     ia_src=ia_cpl, !=1
     &               ia_rad=2, ia_srf=3, ia_dga=4, ia_d4a=5, ia_d5f=6,
     *     ia_d5d=7, ia_d5s=8, ia_12hr=9, ia_filt=10, ia_rad_frc=11,
     *     ia_inst=12

!@var PLE,PLM, PLE_DN ref pressures at upper, middle and lower edge
      REAL*8, DIMENSION(LM), public :: PLE
      REAL*8, DIMENSION(LM), public :: PLE_DN
      REAL*8, DIMENSION(LM+LM_REQ), public :: PLM
!@var P1000K scaling to change reference pressure from 1mb to 1000mb
      REAL*8, public :: P1000K
CXXXX inci,incj NOT GRID-INDPENDENT
!@var inci,incj print increments for i and j, so maps/tables fit on page
      integer, parameter, public ::
     &     inci=(im+35)/36,incj=(JM+23)/24, jmby2=jm/2
!@var linect = current line on page of print out
      integer, public :: linect

!@var LMOMAX max no. of layers in any ocean
      INTEGER, PARAMETER, public :: LMOMAX=50
!@var ZOC, ZOC1 ocean depths for diagnostics (m) (ONLY FOR DEEP OCEAN)
      REAL*8, public :: ZOC(LMOMAX) = 0. , ZOC1(LMOMAX+1) = 0.

!@param L_ROSSBY_NUMBER length scale for budget-page Rossby number
      real*8, parameter, public :: l_rossby_number=1d6 ! 1000 km

!@dbparam NDAA:   DT_DiagA    =  NDAA*DTsrc + 2*DT(dyn)
!@dbparam NDA5k:  DT_Diag5k   =  NDA5k*DTsrc + 2*DT(dyn) SpAnal KE
!@dbparam NDA5d:  DT_Diag5d   =  NDA5d*DTsrc     Consrv  SpAnal dyn
!@dbparam NDA5s:  DT_Diag5s   =  NDA5s*DTsrc     Consrv  SpAnal src
!@dbparam NDASf:  DT_DiagSrfc =  NDASf*DTsrc + DTsrc/NIsurf
!@dbparam NDA4:   DT_Diag4    =  NDA4 *DTsrc   Energy history
      INTEGER, public ::
     &     NDAa=7, NDA5d=1, NDA5k=7, NDA5s=1, NDASf=1, NDA4=24
!@var MODD5K,MODD5S: if MODxxx=0 do xxx, else skip xxx
      INTEGER, public :: MODD5K, MODD5S

      target :: AIJ_loc,JREG,DXYP_BUDG,CONSRV_loc,NOFM

      type(cdl_type), public :: cdl_latbudg,cdl_heights

!@var CDL_J consolidated metadata for AJ output fields in CDL notation
!@+   CDL_REG                         AREG
      type(cdl_type), public :: cdl_j,cdl_reg
!@var CDL_IJ consolidated metadata for AIJ output fields in CDL notation
      type(cdl_type), public ::
     &     cdl_ij_template,cdl_ij_latlon_template,
     &     cdl_ij,cdl_ij_latlon,cdl_ijmm
!@var CDL_JL consolidated metadata for AJL output fields in CDL notation
      type(cdl_type), public :: cdl_jl,cdl_jl_template
!@var CDL_IJL consolidated metadata for AIJL output fields in CDL notation
      type(cdl_type), public ::
     &     cdl_ijl_template,cdl_ijl_latlon_template,
     &     cdl_ijl,cdl_ijl_latlon
!@var CDL_IJK consolidated metadata for AIJK output fields in CDL notation
!@+   (latlon-only for the moment)
      type(cdl_type), public :: cdl_ijk
!@var CDL_CONSRV consolidated metadata for CONSRV output fields in
!@+   CDL notation
      type(cdl_type), public :: cdl_consrv
!@var CDL_DD consolidated metadata for ADIURN output fields in CDL notation
      type(cdl_type), public :: cdl_dd,cdl_hd


c declarations that facilitate switching between restart and acc
c instances of arrays
      target :: aj,aj_out,areg,areg_out
      REAL*8, dimension(:,:,:), public, pointer ::
     &     AJ_ioptr
      REAL*8, dimension(:,:), public, pointer ::
     &     AREG_ioptr

      logical :: Qbp(NTYPE_OUT+1) ! +1 for regions
      public :: Qbp

      END MODULE DIAG_COM

      SUBROUTINE initDiagj
      use filemanager
      USE DIAG_COM, only : ntype_out, terrain, kdiag, Qbp
      IMPLICIT NONE

      LOGICAL qIbp
      INTEGER :: M,iu_Ibp
      character*80 line

        inquire(file='Ibp',exist=qIbp)
        Qbp=.true.
        if(.not.qIbp) then
          call openunit('Ibp',iu_Ibp,.false.,.false.)
          write (iu_Ibp,'(a)') 'List of budget-pages'
          do m = 1,ntype_out
            write (iu_Ibp,'(i3,1x,a)') m,terrain(m)
          end do
            write (iu_Ibp,'(i3,1x,a)') ntype_out+1,'   (REGIONS)'
        else if(kdiag(1).gt.0) then
          Qbp=.false.
          call openunit('Ibp',iu_Ibp,.false.,.true.)
          read (iu_Ibp,'(a)',end=20) line
   10     read (iu_Ibp,'(a)',end=20) line
          read(line,'(i3)') m
          Qbp(m)=.true.
          go to 10
   20     continue
        end if

      END SUBROUTINE initDiagj

      SUBROUTINE ALLOC_DIAG_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID,getDomainBounds,AM_I_ROOT
      USE RESOLUTION, ONLY : IM,LM
      USE ATM_COM, ONLY : lm_req
      USE DIAG_COM, ONLY : KAJ,KCON,KAJL,KASJL,KAIJ,KAIJK,KAIJmm,
     &                   KGZ_max,KOA,KTSF,nwts_ij,KTD,NREG,KAIJL,JM_BUDG
      USE DIAG_COM, ONLY : SQRTM,AJ_loc,JREG,AJL_loc,ASJL_loc
     *     ,AIJ_loc,AIJK_loc,AIJL_loc,AFLX_ST,ftype,ntype
     *     ,Z_inst,RH_inst,T_inst,TDIURN,TSFREZ_loc,OA,P_acc,PM_acc
#if (defined ttc_subdd) || (defined etc_subdd)
     *     ,u_inst,v_inst
#endif
#ifdef ttc_subdd
     *     ,vt_inst
#endif
#ifdef etc_subdd
     *     ,omg_inst,iwc_inst,lwc_inst
     *     ,cldmc_inst,cldss_inst,tlh_inst,llh_inst,dlh_inst,slh_inst
#endif
#if (defined mjo_subdd) || (defined etc_subdd)
     *     ,qsen_avg,qlat_avg,pblht_acc
#endif
#ifdef mjo_subdd
     *     ,E_acc,PW_acc,p_avg,sst_avg,lwu_avg
     *     ,u_avg,v_avg,w_avg,t_avg,q_avg,r_avg,z_avg
#endif
     *     ,saveHCLDI,saveMCLDI,saveLCLDI,saveCTPI,saveTAUI,saveSCLDI
     *     ,saveTCLDI,saveMCCLDTP
      USE DIAG_COM, ONLY : AJ,AJL,ASJL,AJ_OUT,ntype_out
      USE DIAG_COM, ONLY : hemis_j,hemis_jl,vmean_jl,hemis_consrv
     &     ,hemis_ij
      USE DIAG_COM, only : aijmm
#ifdef TRACERS_SPECIAL_Shindell
      USE DIAG_COM, ONLY : o_inst,n_inst,m_inst,x_inst
#endif
#ifdef TES_LIKE_DIAGS
      USE DIAG_COM, ONLY : T_more,Q_more,KGZmore
#ifdef TRACERS_SPECIAL_Shindell
     &     ,o_more,n_more,m_more,x_more
#endif
#endif
#ifndef SCM
      use diag_zonal, only : get_alloc_bounds
#endif
      use fluxes, only : atmocn
      USE DIAG_COM, only : dxyp_budg,nofm,consrv_loc
#ifdef USE_HDIURN
      USE DIAG_COM, only : hdiurn, hdiurn_loc
#endif
      USE DIAG_COM, only : NDIUVAR, NDIUPT
      use Diag_com, only : HR_IN_MONTH
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR,INT_HOURS_PER_DAY
      use Model_Com, only: calendar
      use CalendarMonth_mod
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER
      LOGICAL, SAVE :: init = .false.
      integer :: j_0budg,j_1budg
      integer :: mnth
      integer :: year
      type (CalendarMonth) :: cMonth

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

      call getDomainBounds( grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H  )
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

#ifdef SCM
      j_0budg = 1
      j_1budg = 1
#else
      call get_alloc_bounds(grid,
     &     j_strt_budg=j_0budg,j_stop_budg=j_1budg)
#endif

      hr_in_month = 0
      year = 2000 ! arbitrary (leap days should not be in longest month!)
      do mnth = 1, INT_MONTHS_PER_YEAR
         cMonth = calendar%getCalendarMonth(mnth, year)
         hr_in_month = 
     &        max(hr_in_month, 
     &        INT_HOURS_PER_DAY*cMonth%daysInMonth)
      end do
#ifdef USE_HDIURN
      allocate(
     &     HDIURN(NDIUVAR, NDIUPT, hr_in_month),
     &     HDIURN_loc(NDIUVAR, NDIUPT, hr_in_month),
     &     STAT=IER)
      HDIURN_loc = 0
      HDIURN = 0
#endif

      ALLOCATE(  JREG(I_0H:I_1H, J_0H:J_1H),
     &         SQRTM(I_0H:I_1H, J_0H:J_1H),
     &         STAT = IER)

      ALLOCATE(
     &         AJ_loc(J_0BUDG:J_1BUDG, KAJ, NTYPE),
     &         AJL_loc(J_0BUDG:J_1BUDG, LM, KAJL),
     &         ASJL_loc(J_0BUDG:J_1BUDG,LM_REQ,KASJL),
     &         AIJ_loc(I_0H:I_1H,J_0H:J_1H,KAIJ),
     &         AIJmm(I_0H:I_1H,J_0H:J_1H,KAIJmm),
     &         Z_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         RH_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         T_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         TSFREZ_loc(I_0H:I_1H,J_0H:J_1H,KTSF),
     &         P_acc(I_0H:I_1H,J_0H:J_1H),
     &         PM_acc(I_0H:I_1H,J_0H:J_1H),
#if (defined ttc_subdd) || (defined etc_subdd)
     &         u_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         v_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
#endif
#ifdef ttc_subdd
     &         vt_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
#endif
#ifdef etc_subdd
     &         lwc_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         iwc_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         omg_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         cldmc_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         cldss_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         tlh_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         dlh_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         slh_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         llh_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
#endif
#if (defined mjo_subdd) || (defined etc_subdd)
     &         qlat_avg(I_0H:I_1H,J_0H:J_1H),
     &         qsen_avg(I_0H:I_1H,J_0H:J_1H),
     &         pblht_acc(I_0H:I_1H,J_0H:J_1H),
#endif
#ifdef mjo_subdd
     &         E_acc(I_0H:I_1H,J_0H:J_1H),
     &         PW_acc(I_0H:I_1H,J_0H:J_1H),
     &         p_avg(I_0H:I_1H,J_0H:J_1H),
     &         sst_avg(I_0H:I_1H,J_0H:J_1H),
     &         lwu_avg(I_0H:I_1H,J_0H:J_1H),
     &         u_avg(I_0H:I_1H,J_0H:J_1H,LM),
     &         v_avg(I_0H:I_1H,J_0H:J_1H,LM),
     &         w_avg(I_0H:I_1H,J_0H:J_1H,LM),
     &         t_avg(I_0H:I_1H,J_0H:J_1H,LM),
     &         q_avg(I_0H:I_1H,J_0H:J_1H,LM),
     &         r_avg(I_0H:I_1H,J_0H:J_1H,LM),
     &         z_avg(I_0H:I_1H,J_0H:J_1H,LM),
#endif
     &         TDIURN(I_0H:I_1H,J_0H:J_1H,KTD),
     &         OA(I_0H:I_1H,J_0H:J_1H,KOA),
     &         saveHCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveMCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveLCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveCTPI(I_0H:I_1H,J_0H:J_1H),
     &         saveTAUI(I_0H:I_1H,J_0H:J_1H),
     &         saveSCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveTCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveMCCLDTP(I_0H:I_1H,J_0H:J_1H),
#ifdef TRACERS_SPECIAL_Shindell
     &         O_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         M_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         N_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
     &         X_inst(KGZ_max,I_0H:I_1H,J_0H:J_1H),
#endif
#ifdef TES_LIKE_DIAGS
     &         Q_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
     &         T_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
#ifdef TRACERS_SPECIAL_Shindell
     &         O_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
     &         N_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
     &         M_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
     &         X_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
#endif
#endif
     &         STAT = IER)

      P_acc=0.d0; PM_acc=0.d0
      OA = 0d0

      ALLOCATE( AIJK_loc(I_0H:I_1H,J_0H:J_1H,LM,KAIJK),
     &         AFLX_ST(LM+LM_REQ+1,I_0H:I_1H,J_0H:J_1H,5),
     &         STAT = IER)


      ALLOCATE( AIJL_loc(I_0H:I_1H,J_0H:J_1H,LM,KAIJL))

c allocate master copies of budget- and JK-arrays on root
      if(am_i_root()) then
        ALLOCATE(AJ(JM_BUDG, KAJ, NTYPE),
     &           AJL(JM_BUDG, LM, KAJL),
     &           ASJL(JM_BUDG,LM_REQ,KASJL),
     &           STAT = IER)
        allocate(aj_out(jm_budg,kaj,ntype_out))
        allocate(hemis_j(3,kaj,ntype_out))
        allocate(hemis_jl(3,lm,kajl))
        allocate(vmean_jl(jm_budg+3,1,kajl))
        allocate(hemis_consrv(3,kcon))
        allocate(hemis_ij(1,3,kaij))
      else
        ALLOCATE(AJ(1,1,1),
     &           AJL(1,1,1),
     &           ASJL(1,1,1),
     &        STAT = IER)
        allocate(aj_out(1,1,1))
        allocate(hemis_j(1,1,1))
        allocate(hemis_jl(1,1,1))
        allocate(vmean_jl(1,1,1))
        allocate(hemis_consrv(1,1))
        allocate(hemis_ij(1,1,1))

      endif

      ALLOCATE(FTYPE(NTYPE,I_0H:I_1H,J_0H:J_1H), STAT = IER)
      FTYPE(:,:,:) = 0.d0

      atmocn%jm_budg = jm_budg
      atmocn%area_of_zone => dxyp_budg
      atmocn%consrv => consrv_loc
      atmocn%nofm => nofm

      call initDiagj

#ifndef SCM
      call alloc_gc_com(grid)
#endif

      RETURN
      END SUBROUTINE ALLOC_DIAG_COM

#ifndef SCM
      module gc_com
      use mdiag_com, only : sname_strlen,units_strlen,lname_strlen
      use resolution, only : jm,lm,ls1=>ls1_nominal,pmtop
      use diag_zonal, only : imlonh,jmlat
      use cdl_mod
      implicit none

!**** Based on model top, determine how much of stratosphere is resolved
!**** ISTRAT = 2:          PMTOP <   1 mb
!**** ISTRAT = 1:  1 mb <= PMTOP <  10 mb
!**** ISTRAT = 0: 10 mb <= PMTOP
      integer, parameter :: ! todo when pmtop no longer a parameter: make
                            ! istrat-related arrays allocatable
     &     istrat = min(2,max(0,ceiling(log10(10d0/pmtop))))

!**** KEP depends on whether stratos. EP flux diagnostics are calculated
!**** If dummy EPFLUX is used set KEP=0, otherwise KEP=21
!@var KEP number of lat/height E-P flux diagnostics
      integer :: KEP

!@param JEQ grid box zone around or immediately north of the equator
      INTEGER, PARAMETER, public :: JEQ=1+JM/2

!@var J50N,J70N,J5NUV,J5SUV,J5S,J5N special latitudes
      INTEGER, PARAMETER, public :: J50N  = (50.+90.)*(JM-1)/180.+1.5
      INTEGER, PARAMETER, public :: J70N  = (70.+90.)*(JM-1)/180.+1.5
      INTEGER, PARAMETER, public :: J5NUV = (90.+5.)*(JM-1.)/180.+2.
      INTEGER, PARAMETER, public :: J5SUV = (90.-5.)*(JM-1.)/180.+2.
      INTEGER, PARAMETER, public :: J5N   = (90.+5.)*(JM-1.)/180.+1.5
      INTEGER, PARAMETER, public :: J5S   = (90.-5.)*(JM-1.)/180.+1.5


!@var KAGC number of latitude-height General Circulation diags
!@param KAGCX number of accumulated+derived GC diagnostics
      INTEGER, public :: KAGC!=82+KEP
      INTEGER, PARAMETER, public :: KAGCX=82+21+100
!@var AGC latitude-height General Circulation diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AGC,AGC_loc
     &     ,AGC_out

C NEHIST=(TROPO/L STRAT/M STRAT/U STRAT)X(ZKE/EKE/SEKE/ZPE/EPE)X(SH/NH)
!@param NED number of different energy history diagnostics
!@param NEHIST,HIST_DAYS number of energy history columns,rows (max)
      INTEGER, PARAMETER, public :: NED=10
      INTEGER, PARAMETER, public :: NEHIST=NED*(2+ISTRAT)
      INTEGER, PARAMETER, public :: HIST_DAYS=100
!@var ENERGY energy diagnostics
      REAL*8, DIMENSION(NEHIST,HIST_DAYS), public :: ENERGY

!@param KSPECA,NSPHER number of spectral diagnostics, and harmonics used
      INTEGER, PARAMETER, public :: KSPECA=20
      INTEGER, PARAMETER, public :: NSPHER=4*(2+ISTRAT)
!@var SPECA spectral diagnostics
      REAL*8, DIMENSION((IMLONH+1),KSPECA,NSPHER), public :: SPECA
     &     ,SPECA_OUT
!@var KLAYER index for dividing up atmosphere into layers for spec.anal.
      INTEGER, DIMENSION(LM), public :: KLAYER
!@param PSPEC pressure levels at which layers are seperated and defined
C**** 1000 - 150: troposphere           150 - 10 : low strat.
C****   10 - 1: mid strat               1 and up : upp strat.
      REAL*8, DIMENSION(4), PARAMETER, public ::
     &     PSPEC = (/ 150., 10., 1., 0. /)
!@var LMAX_SPECA upper layer index of each of the SPECA height zones
      INTEGER, PUBLIC :: LMAX_SPECA(2+ISTRAT)

!@param KTPE number of spectral diagnostics for pot. enthalpy
      INTEGER, PARAMETER, public :: KTPE=8
      integer, parameter, public :: NHEMI=2
!@var ATPE pot. enthalpy spectral diagnostics
      REAL*8, DIMENSION(KTPE,NHEMI), public :: ATPE,ATPE_out

!@param NWAV_DAG number of components in spectral diagnostics
      INTEGER, PARAMETER, public :: NWAV_DAG=min(9,imlonh)
!@param Max12HR_sequ,Min12HR_sequ lengths of time series for wave powers
      INTEGER, PARAMETER, public :: Max12HR_sequ=2*31, Min12HR_sequ=2*28
!@param RE_AND_IM complex components of wave power diagnostics
      INTEGER, PARAMETER, public :: RE_AND_IM=2
!@param KWP number of wave power diagnostics
      INTEGER, PARAMETER, public :: KWP=12
!@var WAVE frequency diagnostics (wave power)
      REAL*8, DIMENSION(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP), public ::
     &     WAVE

      character(len=sname_strlen), dimension(kwp), public :: name_wave
      character(len=units_strlen), dimension(kwp), public :: units_wave
      character(len=lname_strlen), dimension(kwp), public :: lname_wave

!@var SNAME_GC Names of lat-pressure GC diagnostics
      character(len=sname_strlen), dimension(kagcx), public :: sname_gc
!@var LNAME_GC,UNITS_GC Descriptions/Units of GC diagnostics
      character(len=lname_strlen), dimension(kagcx), public :: lname_gc
      character(len=units_strlen), dimension(kagcx), public :: units_gc
!@var SCALE_GC printout scaling factors for GC diagnostics
      REAL*8, dimension(kagcx), public :: scale_gc
!@var IA_GC,JGRID_GC,LGRID_GC idacc-numbers,gridtypes for GC diagnostics
      integer, dimension(kagcx), public :: ia_gc,jgrid_gc,lgrid_gc
!@var DENOM_GC index of AGC element to use as weight
      integer, dimension(kagcx), public :: denom_gc
!@var POW_GC printed output scaled by 10**(-pow_gc)
      integer, dimension(kagcx), public :: pow_gc,pow_gc_vmean
!@var HEMIS_GC hemispheric/global averages of AGC
!@var VMEAN_GC vertical sums of AGC
      real*8, dimension(:,:,:), allocatable, public :: hemis_gc,vmean_gc
!@var lat_gc latitudes of the primary grid for GC diagnostics
      real*8, dimension(jmlat), public :: lat_gc,lat_gc2

!@var CDL_GC consolidated metadata for AGC output fields in CDL notation
      type(cdl_type), public :: cdl_gc

      target :: agc,agc_out
      REAL*8, dimension(:,:,:), public, pointer ::
     &     AGC_ioptr

      target :: speca,speca_out,atpe,atpe_out
      REAL*8, dimension(:,:), public, pointer ::
     &     ATPE_ioptr
      REAL*8, dimension(:,:,:), public, pointer ::
     &     SPECA_ioptr

      end module gc_com

      subroutine alloc_gc_com(grid)
      use domain_decomp_atm, only : dist_grid,am_i_root
      use resolution, only : lm
      use gc_com, only : kep,kagc,jmlat,agc,agc_loc,agc_out
      use gc_com, only : hemis_gc,vmean_gc
      use diag_zonal, only : get_alloc_bounds
      implicit none
      type (dist_grid), intent(in) :: grid
      integer :: ier
      integer :: j_0jk,j_1jk

      call get_alloc_bounds(grid,
     &     j_strt_jk=j_0jk,j_stop_jk=j_1jk)

      call get_kep(kep)
      kagc = 82+kep

      allocate(agc_loc(j_0jk:j_1jk,lm,kagc))

c allocate master copies of budget- and jk-arrays on root
      if(am_i_root()) then
        allocate(agc(jmlat,lm,kagc),
     &           agc_out(jmlat,lm,kagc))
        allocate(hemis_gc(3,lm,kagc))
        allocate(vmean_gc(jmlat+3,1,kagc))
      else
        allocate(agc(1,1,1),
     &           agc_out(1,1,1),
     &        stat = ier)
        allocate(hemis_gc(1,1,1))
        allocate(vmean_gc(1,1,1))
      endif

      return
      end subroutine alloc_gc_com
#endif

      SUBROUTINE ALLOC_ijdiag_glob
!@sum  To allocate large global arrays only when needed
!@auth NCCS (Goddard) Development Team
      USE RESOLUTION, ONLY : IM,JM,LM
      USE DOMAIN_DECOMP_ATM, Only : AM_I_ROOT
      USE DIAG_COM, ONLY : KAIJ,KAIJK,KOA,KTSF,KTD,KAIJL
      USE DIAG_COM, ONLY : AIJ,AIJK,AIJL,TSFREZ,TDIURN_GLOB,OA_GLOB
      IMPLICIT NONE

      if(AM_I_ROOT()) then
         ALLOCATE(AIJ(IM,JM,KAIJ),
     &        TSFREZ(IM,JM,KTSF),
     &        AIJK(IM,JM,LM,KAIJK),
     &        AIJL(IM,JM,LM,KAIJL),
     &        TDIURN_glob(IM, JM, KTD),
     &        OA_glob(IM, JM, KOA))
      else
         ALLOCATE(AIJ(1,1,1),
     &        TSFREZ(1,1,1),
     &        AIJK(1,1,1,1),
     &        AIJL(1,1,1,1),
     &        TDIURN_glob(1,1,1),
     &        OA_glob(1,1,1))
      end if


      RETURN
      END SUBROUTINE ALLOC_ijdiag_glob

      SUBROUTINE DEALLOC_ijdiag_glob
!@sum  To deallocate large global arrays not currently needed
!@auth NCCS (Goddard) Development Team
      USE DOMAIN_DECOMP_ATM, Only : AM_I_ROOT
      USE DIAG_COM, ONLY : AIJ,AIJK,AIJL,TSFREZ,TDIURN_GLOB,OA_GLOB

      IMPLICIT NONE
      DEALLOCATE(AIJ,AIJK,AIJL,TSFREZ,TDIURN_glob,OA_glob)

      RETURN
      END SUBROUTINE DEALLOC_ijdiag_glob

      Subroutine Gather_Diagnostics()
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : pack_data
      use diag_com
      implicit none
      call gather_zonal_diags
      call collect_scalars
      CALL PACK_DATA(grid,  TSFREZ_loc, TSFREZ)
      CALL PACK_DATA(grid,  AIJ_loc,    AIJ)
      CALL PACK_DATA(grid,  AIJK_loc,   AIJK)
      CALL PACK_DATA(grid,  AIJL_loc,   AIJL)
      CALL PACK_DATA(grid,  TDIURN, TDIURN_glob)
      CALL PACK_DATA(grid,  OA, OA_glob)
      return
      End Subroutine Gather_Diagnostics

      Subroutine Collect_Scalars()
      use precision_mod
      use domain_decomp_atm, only : grid,sumxpe,am_i_root
      use diag_com
      implicit none
      integer :: k
      real*8, dimension(:,:), allocatable :: consrv_sv
      CALL SUMXPE(AREG_loc, AREG, increment=.true.)
      AREG_loc(:,:)=0.
      CALL SUMXPE(ADIURN_loc, ADIURN, increment=.true.)
      ADIURN_loc=0
#ifdef USE_HDIURN
      CALL SUMXPE(HDIURN_loc, HDIURN, increment=.true.)
      HDIURN_loc=0
#endif
      CALL SUMXPE(AISCCP_loc, AISCCP, increment=.true.)
      AISCCP_loc=0

      if(am_i_root()) then
c for reproducibility on different numbers of processors
        call reduce_precision(areg,1d-9)
        call reduce_precision(aisccp,1d-9)
      endif

c CONSRV is a mixture of accumulations and instantaneous values, hence the
c more complicated logic
      if(am_i_root()) then
        allocate(consrv_sv(jm_budg,kcon)); consrv_sv = consrv ! pde hack
        do k=1,kcon
          if(ia_con(k).eq.ia_inst) consrv(:,k)=0. ! collect instantaneous values
        enddo
      endif
      CALL SUMXPE(CONSRV_loc, CONSRV, increment=.true.)
      do k=1,kcon
        if(ia_con(k).ne.ia_inst) consrv_loc(:,k)=0. ! keep instantaneous values
      enddo
      if(am_i_root()) then ! pde hack
        do k=1,kcon
          if(ia_con(k).eq.ia_inst .and. all(consrv(:,k)==0.))
     &         consrv(:,k) = consrv_sv(:,k)
        enddo
        deallocate(consrv_sv)
      endif

      return
      End Subroutine Collect_Scalars

      Subroutine Gather_zonal_diags()
      use domain_decomp_atm, only : grid
      use diag_com
#ifdef SCM
      implicit none
      AJ = AJ_loc
      AJL = AJL_loc
      ASJL = ASJL_loc
#else
      use diag_zonal, only : pack_lc
      use gc_com
      implicit none
      call pack_lc(grid, AJ_loc,     AJ)
      call pack_lc(grid, AJL_loc,    AJL)
      call pack_lc(grid, ASJL_loc,   ASJL)
      call pack_lc(grid, AGC_loc,    AGC)
#endif
      return
      End Subroutine Gather_zonal_diags

      Subroutine Scatter_zonal_diags()
      use domain_decomp_atm, only : grid
      use diag_com
#ifdef SCM
      implicit none
      AJ_loc = AJ
      AJL_loc = AJL
      ASJL_loc = ASJL
#else
      use diag_zonal, only : unpack_lc
      use gc_com
      implicit none
      call unpack_lc  (grid, aj,     aj_loc)
      call unpack_lc  (grid, ajl,    ajl_loc)
      call unpack_lc  (grid, asjl,   asjl_loc)
      call unpack_lc  (grid, agc,    agc_loc)
#endif
      return
      End Subroutine Scatter_zonal_diags


C**** Routines associated with the budget grid
#ifndef SCM
      SUBROUTINE SET_J_BUDG
!@sum set_j_budg definition for grid points map to budget-grid zonal means
!@auth Gavin Schmidt
      USE GEOM, only : lat2d_dg
      USE DIAG_COM, only : jm_budg,j_budg,j_0b,j_1b
      USE DOMAIN_DECOMP_ATM, only :GRID,getDomainBounds, halo_update
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
      INTEGER :: I,J,J_0,J_1,I_0,I_1,J_0H,J_1H,I_0H,I_1H
      INTEGER :: IER
      Real*8  :: dLATD  !  latitudinal budget spacing in degrees

C**** define atmospheric grid
      call getDomainBounds(grid,J_STRT=J_0,J_STOP=J_1,J_STRT_HALO=J_0H,
     *     J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO ; I_1H = grid%I_STOP_HALO
      I_0 = grid%I_STRT ; I_1 = grid%I_STOP

      ALLOCATE( J_BUDG(I_0H:I_1H, J_0H:J_1H), STAT = IER)

C**** Define mapping from actual lon/lat point to budget grid
C**** this should be valid for all grids (lat/lon, cubed sphere,...)
      dLATD = 180d0 / JM_BUDG
      If (JM_BUDG == 46)  dLATD = 4
      If (JM_BUDG == 24)  dLATD = 8
      DO J=J_0,J_1
        DO I=I_0,I_1
           J_BUDG(I,J) = Nint (LAT2D_DG(I,J)/dLATD + (JM_BUDG+1)/2d0)
           If (J_BUDG(I,J) < 1)        J_BUDG(I,J) = 1
           If (J_BUDG(I,J) > JM_BUDG)  J_BUDG(I,J) = JM_BUDG
        END DO
      END DO

      call halo_update(grid, J_BUDG)

C**** define limits on budget indices for each processor
      j_0b=MINVAL( J_BUDG(I_0:I_1,J_0:J_1) )
      j_1b=MAXVAL( J_BUDG(I_0:I_1,J_0:J_1) )

      RETURN
      END SUBROUTINE SET_J_BUDG


      subroutine set_wtbudg()
!@sum Precomputes area weights for zonal means on budget grid
!auth Denis Gueyffier
      USE GEOM, only : axyp, imaxj
#ifndef CUBED_SPHERE
      USE DIAG_COM, only : im,fim
#endif
      USE DIAG_COM, only : jm_budg,wtbudg,wtbudg2,lat_budg,dxyp_budg
     &     ,j_budg
      USE DOMAIN_DECOMP_ATM, only :GRID,getDomainBounds
      IMPLICIT NONE
      INTEGER :: I,J,J_0,J_1,I_0,I_1
      INTEGER :: IER
      real*8 :: dlat_budg,fjeq_budg

      call getDomainBounds(grid, J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT ; I_1 = grid%I_STOP

      ALLOCATE( wtbudg(I_0:I_1, J_0:J_1), STAT = IER)  !deallocated near very end, stays in memory all the time
      ALLOCATE( wtbudg2(I_0:I_1, J_0:J_1), STAT = IER)

C**** get the nominal latitudes of the budget grid
      dlat_budg = 180./REAL(JM_budg)  ! for full polar box
      if(jm_budg.eq.46) dlat_budg=4. ! 1/2 box at pole for 4 deg res.
      if(jm_budg.eq.24) dlat_budg=8. ! 1/4 box at pole for 8 deg res.
      lat_budg(1) = -90.; lat_budg(jm_budg) = +90.;
      fjeq_budg = .5*(1+jm_budg)
      do j=2,jm_budg-1
        lat_budg(j) = dlat_budg*(j-fjeq_budg)
      enddo

      call set_budg_area()

C**** Compute area weights of zonal bands, including adjustment for polar
C**** box in lat/lon case.
      do J=J_0,J_1
         do I=I_0,I_1
            wtbudg2(I,J) = axyp(I,J)/dxyp_budg(J_BUDG(I,J))
#ifdef CUBED_SPHERE
            wtbudg (I,J) = wtbudg2(I,J)
#else
            wtbudg (I,J) = wtbudg2(I,J)*fim/real(imaxj(j),kind=8)
#endif
         enddo
      enddo

      RETURN
      END SUBROUTINE set_wtbudg

      subroutine set_budg_area()
!@sum  pre-computes area of budget-grid band (zig-zag band on the cubed sphere)
!@+    accross processors. Several processors can contribute to same band
!@auth Denis Gueyffier
      use GEOM, only: axyp
      use DIAG_COM, only : dxyp_budg,dxyp_budg_loc,JM_BUDG,J_BUDG
      USE DOMAIN_DECOMP_ATM, only :grid,getDomainBounds,sumxpe,broadcast
      IMPLICIT NONE
      INTEGER :: I,J,J_0,J_1,I_0,I_1
      logical :: increment

      call getDomainBounds(grid, J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT ; I_1 = grid%I_STOP

      dxyp_budg_loc(:)=0.0

      do J=J_0,J_1
         do I=I_0,I_1
            dxyp_budg_loc(J_BUDG(I,J))=dxyp_budg_loc(J_BUDG(I,J))
     &       +axyp(I,J)
         enddo
      enddo

      increment = .false.
      call SUMXPE(dxyp_budg_loc,dxyp_budg,increment)   !summing global area
      call broadcast(grid,dxyp_budg)

      return
      end subroutine set_budg_area
#endif /* not SCM */

      SUBROUTINE INC_AJ(I,J,ITYPE,J_DIAG,ACC)
!@sum inc_aj grid dependent incrementer for zonal mean budget diags
!@auth Gavin Schmidt
      USE DIAG_COM, only : aj=>aj_loc, wtbudg, j_budg
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
!@var ITYPE is the surface type
      INTEGER, INTENT(IN) :: I,J,ITYPE
!@var J_DIAG is placing of diagnostic element
      INTEGER, INTENT(IN) :: J_DIAG
!@var ACC value of the diagnostic to be accumulated
      REAL*8, INTENT(IN) :: ACC

C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      !   wtbudg area-weight =1 on lat-lon, <1 on cubed sphere
      AJ(J_BUDG(I,J),J_DIAG,ITYPE) = AJ(J_BUDG(I,J),J_DIAG,ITYPE)
     &     +wtbudg(I,J)*ACC

      RETURN
      END SUBROUTINE INC_AJ

      SUBROUTINE INC_AREG(I,J,JR,J_DIAG,ACC)
!@sum inc_areg incrementer for regional budget diags
!@auth Gavin Schmidt
      USE DIAG_COM, only : areg=>areg_loc
      USE GEOM, only : axyp
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
!@var JR is the region
      INTEGER, INTENT(IN) :: I,J,JR
!@var J_DIAG identifier of diagnostic element
      INTEGER, INTENT(IN) :: J_DIAG
!@var ACC value of the diagnostic to be accumulated
      REAL*8, INTENT(IN) :: ACC

      AREG(JR,J_DIAG) = AREG(JR,J_DIAG)+ACC*AXYP(I,J)

      RETURN
      END SUBROUTINE INC_AREG

      SUBROUTINE INC_AJL(I,J,L,JL_INDEX,ACC)
!@sum inc_ajl adds ACC located at atmospheric gridpoint I,J,L
!@+   to the latitude-height zonal sum AJL(J,L,JL_INDEX).
!@auth M. Kelley
      USE DIAG_COM, only : ajl=>ajl_loc,wtbudg,j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: JL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC

      AJL(J_BUDG(I,J),L,JL_INDEX) = AJL(J_BUDG(I,J),L,JL_INDEX)
     &     +wtbudg(I,J)*ACC

      RETURN
      END SUBROUTINE INC_AJL

      SUBROUTINE INC_AJL2(I,J,L,JL_INDEX,ACC)
c temporary variant of inc_ajl without any weighting
      USE DIAG_COM, only : ajl=>ajl_loc,j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: JL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC

      AJL(J_BUDG(I,J),L,JL_INDEX) = AJL(J_BUDG(I,J),L,JL_INDEX)
     &     +ACC

      RETURN
      END SUBROUTINE INC_AJL2

      SUBROUTINE INC_ASJL(I,J,L,SJL_INDEX,ACC)
!@sum inc_asjl adds ACC located at atmospheric gridpoint I,J,L
!@+   to the latitude-height zonal sum ASJL(J,L,JL_INDEX).
!@+   This is a trivial version for the latlon grid.
!@auth M. Kelley
      USE DIAG_COM, only : asjl=>asjl_loc,wtbudg,j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var SJL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: SJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC

      ASJL(J_BUDG(I,J),L,SJL_INDEX) = ASJL(J_BUDG(I,J),L,SJL_INDEX)
     *     +wtbudg(I,J)*ACC

      RETURN
      END SUBROUTINE INC_ASJL

      subroutine def_rsf_acc(fid,r4_on_disk)
!@sum  def_rsf_acc defines accumulation array structure in restart/acc files
!@auth M. Kelley
!@ver  beta
      use model_com, only : idacc
      use mdiag_com, only : monacc
      use diag_com, only :
     &     aj=>aj_ioptr,areg=>areg_ioptr,
     &     aij=>aij_loc,aijl=>aijl_loc,aijk=>aijk_loc, ! dist
     &     oa,tdiurn,aijmm,                            ! dist
     &     ajl,asjl,consrv,
     &     adiurn,aisccp
#ifdef USE_HDIURN
      use diag_com, only :  hdiurn
#endif
#ifndef SCM
      use gc_com, only : agc=>agc_ioptr,energy,wave,
     &     speca=>speca_ioptr,atpe=>atpe_ioptr
#endif
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      implicit none
      integer fid   !@var fid file id
      logical :: r4_on_disk  !@var r4_on_disk if true, real*8 stored as real*4
      real*8 :: r8dum
      character(len=20) :: tmpstr

      call defvar(grid,fid,idacc,'monacc(twelve)')
      call defvar(grid,fid,idacc,'idacc(nsampl)')
      call defvar(grid,fid,tdiurn,'tdiurn(dist_im,dist_jm,ktd)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,oa,'oa(dist_im,dist_jm,koa)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,aij,'aij(dist_im,dist_jm,kaij)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,aijmm,'aijmm(dist_im,dist_jm,kaijmm)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,aijl,'aijl(dist_im,dist_jm,lm,kaijl)',
     &     r4_on_disk=r4_on_disk)
#ifndef CUBED_SPHERE
#ifndef SCM
      call defvar(grid,fid,aijk,'aijk(dist_im,dist_jm,lm,kaijk)',
     &     r4_on_disk=r4_on_disk)
#endif
#endif

#ifdef CACHED_SUBDD
      if(.not.r4_on_disk) call def_rsf_subdd_acc(fid,r4_on_disk)
#endif

      call defvar(grid,fid,aj,'aj(jm_budg,kaj,ntype)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,ajl,'ajl(jm_budg,lm,kajl)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,asjl,'asjl(jm_budg,lm_req,kasjl)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,consrv,'consrv(jm_budg,kcon)',
     &     r4_on_disk=r4_on_disk)

      call defvar(grid,fid,areg,'areg(nreg,kaj)',r4_on_disk=r4_on_disk)

      call defvar(grid,fid,aisccp,'aisccp(ntau,npres,nisccp)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,adiurn,
     &     'adiurn(ndiuvar,ndiupt,hr_in_day)',r4_on_disk=r4_on_disk)
#ifdef USE_HDIURN
      call defvar(grid,fid,hdiurn,
     &     'hdiurn(ndiuvar,ndiupt,hr_in_month)',r4_on_disk=r4_on_disk)
#endif

#ifndef SCM
      call defvar(grid,fid,agc,'agc(jmlat,lm,kagc)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,energy,'energy(nehist,hist_days)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,speca,
     &     'speca(imlonh_plus_1,kspeca,nspher)',r4_on_disk=r4_on_disk)
      call defvar(grid,fid,atpe,'atpe(ktpe,nhemi)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,wave,
     &     'wave(re_and_im,max12hr_sequ,nwav_dag,kwp)',
     &     r4_on_disk=r4_on_disk)
#endif

      if(.not.r4_on_disk) then          ! reproducibility info
        tmpstr = 'is_npes_reproducible' ! for checkpoint files
        call defvar(grid,fid,r8dum,tmpstr)
        call write_attr(grid,fid,tmpstr,'aisccp' ,'no')
        call write_attr(grid,fid,tmpstr,'areg'   ,'no')
        call write_attr(grid,fid,tmpstr,'speca'  ,'no')
        call write_attr(grid,fid,tmpstr,'energy' ,'no')
c       if(using Ocean R not on atm grid)
c       call write_attr(grid,fid,tmpstr,'consrv' ,'no')
#ifdef CUBED_SPHERE
        call write_attr(grid,fid,tmpstr,'aj'     ,'no')
        call write_attr(grid,fid,tmpstr,'ajl'    ,'no')
        call write_attr(grid,fid,tmpstr,'asjl'   ,'no')
        call write_attr(grid,fid,tmpstr,'agc'    ,'no')
        call write_attr(grid,fid,tmpstr,'consrv' ,'no')
#ifdef TRACERS_ON /* move these declarations to the tracer code? */
        call write_attr(grid,fid,tmpstr,'tajln'  ,'no')
        call write_attr(grid,fid,tmpstr,'tajls'  ,'no')
        call write_attr(grid,fid,tmpstr,'tconsrv','no')
#endif
#endif
      endif

#ifdef TRACERS_ON
      call def_rsf_trdiag(fid,r4_on_disk)
#endif

      return
      end subroutine def_rsf_acc

      subroutine new_io_acc(fid,iaction)
!@sum  new_io_acc read/write accumulation arrays from/to restart/acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use resolution, only : im,jm,lm
      use model_com, only : ioread,iowrite,iowrite_single,iowrite_mon
     &     ,idacc
c i/o pointers point to:
c    primary instances of arrays when writing restart files
c    extended/rescaled instances of arrays when writing acc files
      use mdiag_com, only : monacc
      use diag_com, only : kaijl,
     &     aj=>aj_ioptr,areg=>areg_ioptr,
     &     aij=>aij_loc,aijl=>aijl_loc,aijk=>aijk_loc, ! dist
     &     oa,tdiurn,aijmm,                            ! dist
     &     ajl,asjl,consrv,
     &     adiurn,aisccp
#ifdef USE_HDIURN
      use diag_com, only :  hdiurn
#endif
#ifndef SCM
      use gc_com, only : agc=>agc_ioptr,energy,wave,
     &     speca=>speca_ioptr,atpe=>atpe_ioptr
#endif
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : hasNorthPole, hasSouthPole
      use pario, only : write_dist_data,read_dist_data
     &     ,write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: k,l
      select case (iaction)
      case (iowrite,iowrite_single) ! output to restart or acc file
        if(iaction.eq.iowrite) then ! already done for iowrite_single
          call gather_zonal_diags
          call collect_scalars
        else  ! pole fills needed for acc-files
          do k=1,kaijl
            do l=1,lm
              if(hasSouthPole(grid)) then
                aijl(2:im, 1,l,k) = aijl(1, 1,l,k)
              endif
              if(hasNorthPole(grid)) then
                aijl(2:im,jm,l,k) = aijl(1,jm,l,k)
              endif
            enddo
          enddo
        endif
        call write_data(grid,fid,'monacc',monacc)
        call write_data(grid,fid,'idacc',idacc)
        call write_data(grid,fid,'aisccp',aisccp)
        call write_data(grid,fid,'adiurn',adiurn)
#ifdef USE_HDIURN
        call write_data(grid,fid,'hdiurn',hdiurn)
#endif
        call write_dist_data(grid,fid,'tdiurn',tdiurn)
        call write_dist_data(grid,fid,'oa',oa)
        call write_dist_data(grid,fid,'aij',aij)
        call write_dist_data(grid,fid,'aijmm',aijmm)
        call write_dist_data(grid,fid,'aijl',aijl)
#ifndef CUBED_SPHERE
#ifndef SCM
        call write_dist_data(grid,fid,'aijk',aijk)
#endif
#endif

        call write_data(grid,fid,'aj',aj)
        call write_data(grid,fid,'ajl',ajl)
        call write_data(grid,fid,'asjl',asjl)
        call write_data(grid,fid,'consrv',consrv)

        call write_data(grid,fid,'areg',areg)

#ifndef SCM
        call write_data(grid,fid,'agc',agc)
        call write_data(grid,fid,'energy',energy)
        call write_data(grid,fid,'speca',speca)
        call write_data(grid,fid,'atpe',atpe)
        call write_data(grid,fid,'wave',wave)
#endif

      case (ioread)            ! input from restart or acc file
c for which scalars is bcast_all=.true. necessary?
        call read_data(grid,fid,'monacc',monacc,bcast_all=.true.)
        call read_data(grid,fid,'idacc',idacc,bcast_all=.true.)
        call read_data(grid,fid,'aisccp',aisccp,bcast_all=.true.)
        call read_data(grid,fid,'adiurn',adiurn,bcast_all=.true.)
#ifdef USE_HDIURN
        call read_data(grid,fid,'hdiurn',hdiurn,bcast_all=.true.)
#endif
        call read_dist_data(grid,fid,'tdiurn',tdiurn)
        call read_dist_data(grid,fid,'oa',oa)
        call read_dist_data(grid,fid,'aij',aij)
        call read_dist_data(grid,fid,'aijmm',aijmm)
        call read_dist_data(grid,fid,'aijl',aijl)
#ifndef CUBED_SPHERE
#ifndef SCM
        call read_dist_data(grid,fid,'aijk',aijk)
#endif
#endif

        call read_data(grid,fid,'aj',aj)
        call read_data(grid,fid,'ajl',ajl)
        call read_data(grid,fid,'asjl',asjl)
        call read_data(grid,fid,'consrv',consrv)

        call read_data(grid,fid,'areg',areg)

#ifndef SCM
        call read_data(grid,fid,'agc',agc)
        call read_data(grid,fid,'energy',energy,bcast_all=.true.)
        call read_data(grid,fid,'speca',speca,bcast_all=.true.)
        call read_data(grid,fid,'atpe',atpe,bcast_all=.true.)
        call read_data(grid,fid,'wave',wave,bcast_all=.true.)
#endif

        call scatter_zonal_diags

      case (iowrite_mon)            ! specials
        call stop_model('new_io_acc: fix io_oda call',255)
c            If (AM_I_ROOT())
c     *           call openunit(aDATE(1:7)//'.oda'//XLABEL(1:LRUNID)
c     *           ,iu_ODA,.true.,.false.)
c            call io_oda(iu_ODA,Itime,iowrite,ioerr)
c            IF (AM_I_ROOT()) call closeunit(iu_ODA)
      end select

#ifdef TRACERS_ON
        call new_io_trdiag (fid,iaction)
#endif

#ifdef CACHED_SUBDD
      if(iaction.eq.iowrite) then
        call write_subdd_accdata(fid,iaction)
      endif
#endif

      return
      end subroutine new_io_acc

      subroutine def_rsf_longacc(fid,r4_on_disk)
!@sum  def_rsf_longacc defines accumulation array structure in restart/acc files
!@auth M. Kelley
!@ver  beta
      use diag_com, only : tsfrez=>tsfrez_loc,keyct,keynr
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid           !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
      call defvar(grid,fid,keyct,'keyct')
      call defvar(grid,fid,keynr,'keynr(nkeynr,nkeymo)')
      call defvar(grid,fid,tsfrez,'tsfrez(dist_im,dist_jm,ktsf)',
     &     r4_on_disk=r4_on_disk)
      return
      end subroutine def_rsf_longacc

      subroutine new_io_longacc(fid,iaction)
!@sum  new_io_longacc read/write accumulation arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use diag_com, only : tsfrez=>tsfrez_loc,keyct,keynr
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
     &     ,write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart or acc file
        call write_data(grid,fid,'keyct',keyct)
        call write_data(grid,fid,'keynr',keynr)
        call write_dist_data(grid,fid,'tsfrez',tsfrez)
      case (ioread)            ! input from restart or acc file
c for which scalars is bcast_all=.true. necessary?
        call read_data(grid,fid,'keyct',keyct,bcast_all=.true.)
        call read_data(grid,fid,'keynr',keynr,bcast_all=.true.)
        call read_dist_data(grid,fid,'tsfrez',tsfrez)
      end select
      return
      end subroutine new_io_longacc

c def_rsf_subdd
      subroutine def_rsf_subdd(fid)
!@sum def_rsf_subdd defines write/read structure of accumulation arrays
!@+                 for subdaily diagnostics to/from restart files
!@auth Jan Perlwitz

      use domain_decomp_atm, only: grid
      use pario, only: defvar
      use diag_com, only: P_acc,PM_acc

      implicit none

      integer,intent(in) :: fid

      call defvar(grid,fid,P_acc,'P_acc(dist_im,dist_jm)')
      call defvar(grid,fid,PM_acc,'PM_acc(dist_im,dist_jm)')

      return
      end subroutine def_rsf_subdd

c new_io_subdd
      subroutine new_io_subdd(fid,iaction)
!@sum new_io_subdd write/read of accumulation arrays for subdaily diagnostics
!+                 to/from restart files
!@auth Jan Perlwitz

      use domain_decomp_atm, only: grid
      use pario, only: read_dist_data,write_dist_data
      use model_com, only: ioread,iowrite
      use diag_com, only: P_acc,PM_acc

      implicit none

      integer,intent(in) :: iaction,fid

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'P_acc',P_acc)
        call write_dist_data(grid,fid,'PM_acc',PM_acc)

      case (ioread)             ! input from restart file
        call read_dist_data(grid,fid,'P_acc',P_acc)
        call read_dist_data(grid,fid,'PM_acc',PM_acc)

      end select

      return
      end subroutine new_io_subdd

      subroutine def_meta_atmacc(fid)
!@sum  def_meta_atmacc defines metadata in atm acc files
!@auth M. Kelley
!@ver  beta
      use diag_com, only :
     &     ia_j,ia_jl,ia_ij,ia_ijl,ia_con,ia_ijk,
     &     name_j,name_reg,sname_jl,name_ij,name_ijl,name_dd,
     &     name_consrv,name_ijk,
     &     cdl_j,cdl_reg,cdl_jl,
     &     cdl_ij,cdl_ijl,cdl_ij_latlon,cdl_ijl_latlon,
     &     cdl_dd,cdl_hd,cdl_consrv,cdl_ijk,
     &     hemis_j,hemis_jl,vmean_jl,hemis_consrv,
     &     hemis_ij,
     &     scale_j,scale_jl,scale_ij,scale_ijl,scale_dd,scale_con,
     &     scale_ijk,
     &     iden_j,iden_reg,denom_jl,denom_ijl,denom_ij,denom_dd,
     &     denom_ijk,
     &     lm,
     &     isccp_diags,isccp_press,isccp_tau,isccp_late,wisccp,
     &     scale_ijmm,name_ijmm,cdl_ijmm,
     &     write_regions
#ifndef SCM
      use gc_com, only : kagc,ia_gc,sname_gc,cdl_gc,hemis_gc,
     &     vmean_gc,scale_gc,denom_gc,lmax_speca
#endif
      use geom, only : axyp
#ifdef CUBED_SPHERE
      use geom, only : lon2d_dg,lat2d_dg,lonbds,latbds
#endif
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id
      integer :: int_dummy
      real*8 :: r8dum

#ifdef CUBED_SPHERE
      call defvar(grid,fid,lon2d_dg,'lon(dist_im,dist_jm)')
      call defvar(grid,fid,lat2d_dg,'lat(dist_im,dist_jm)')
      call defvar(grid,fid,lonbds,'lonbds(four,dist_im,dist_jm)')
      call defvar(grid,fid,latbds,'latbds(four,dist_im,dist_jm)')
#endif

      call defvar(grid,fid,axyp,'axyp(dist_im,dist_jm)')
      call defvar(grid,fid,r8dum,'time')
      call write_attr(grid,fid,'time','reduction','avg')

      call write_attr(grid,fid,'aj','reduction','sum')
      call write_attr(grid,fid,'aj','split_dim',2)
#ifndef SCM
      call defvar(grid,fid,hemis_j,'hemis_aj(shnhgm,kaj,ntype)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_aj','reduction','sum')
#endif
      call defvar(grid,fid,ia_j,'ia_aj(kaj)')
      call defvar(grid,fid,scale_j,'scale_aj(kaj)')
      call defvar(grid,fid,iden_j,'denom_aj(kaj)')
      call defvar(grid,fid,name_j,'sname_aj(sname_strlen,kaj)')
      call defvar_cdl(grid,fid,cdl_j,'cdl_aj(cdl_strlen,kcdl_aj)')

      if(write_regions) then
      call write_attr(grid,fid,'areg','reduction','sum')
      call write_attr(grid,fid,'areg','split_dim',2)
      call defvar(grid,fid,ia_j,'ia_areg(kaj)')
      call defvar(grid,fid,scale_j,'scale_areg(kaj)')
      call defvar(grid,fid,iden_reg,'denom_areg(kaj)')
      call defvar(grid,fid,name_reg,'sname_areg(sname_strlen,kaj)')
      call defvar_cdl(grid,fid,cdl_reg,
     &     'cdl_areg(cdl_strlen,kcdl_areg)')
      endif

      call write_attr(grid,fid,'consrv','reduction','sum')
      call write_attr(grid,fid,'consrv','split_dim',2)
#ifndef SCM
      call defvar(grid,fid,hemis_consrv,'hemis_consrv(shnhgm,kcon)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_consrv','reduction','sum')
#endif
      call defvar(grid,fid,ia_con,'ia_consrv(kcon)')
      call defvar(grid,fid,scale_con,'scale_consrv(kcon)')
      call defvar(grid,fid,name_consrv,
     &     'sname_consrv(sname_strlen,kcon)')
      call defvar_cdl(grid,fid,cdl_consrv,
     &     'cdl_consrv(cdl_strlen,kcdl_consrv)')

      call write_attr(grid,fid,'ajl','reduction','sum')
      call write_attr(grid,fid,'ajl','split_dim',3)
#ifndef SCM
      call defvar(grid,fid,hemis_jl,'hemis_ajl(shnhgm,lm,kajl)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_ajl','reduction','sum')
      call defvar(grid,fid,vmean_jl,'vmean_ajl(jm_budg_plus3,one,kajl)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'vmean_ajl','reduction','sum')
#endif
      call defvar(grid,fid,ia_jl,'ia_ajl(kajl)')
      call defvar(grid,fid,scale_jl,'scale_ajl(kajl)')
      call defvar(grid,fid,denom_jl,'denom_ajl(kajl)')
      call defvar(grid,fid,sname_jl,'sname_ajl(sname_strlen,kajl)')
      call defvar_cdl(grid,fid,cdl_jl,'cdl_ajl(cdl_strlen,kcdl_ajl)')

#ifndef SCM
      call write_attr(grid,fid,'agc','reduction','sum')
      call write_attr(grid,fid,'agc','split_dim',3)
      call defvar(grid,fid,hemis_gc,'hemis_agc(shnhgm,lm,kagc)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_agc','reduction','sum')
      call defvar(grid,fid,vmean_gc,'vmean_agc(jmlat_plus3,one,kagc)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'vmean_agc','reduction','sum')
      call defvar(grid,fid,ia_gc(1:kagc),'ia_agc(kagc)')
      call defvar(grid,fid,scale_gc(1:kagc),'scale_agc(kagc)')
      call defvar(grid,fid,denom_gc(1:kagc),'denom_agc(kagc)')
      call defvar(grid,fid,sname_gc(1:kagc),
     &     'sname_agc(sname_strlen,kagc)')
      call defvar_cdl(grid,fid,cdl_gc,
     &     'cdl_agc(cdl_strlen,kcdl_agc)')
#endif

      call write_attr(grid,fid,'aij','reduction','sum')
      call write_attr(grid,fid,'aij','split_dim',3)
      call defvar(grid,fid,ia_ij,'ia_aij(kaij)')
      call defvar(grid,fid,scale_ij,'scale_aij(kaij)')
      call defvar(grid,fid,denom_ij,'denom_aij(kaij)')
      call defvar(grid,fid,name_ij,'sname_aij(sname_strlen,kaij)')
      call defvar_cdl(grid,fid,cdl_ij,'cdl_aij(cdl_strlen,kcdl_aij)')
#ifndef SCM
      call defvar(grid,fid,hemis_ij,'hemis_aij(one,shnhgm,kaij)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_aij','reduction','sum')
#endif
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_ij_latlon,
     &     'cdl_aij_latlon(cdl_strlen,kcdl_aij_latlon)')
#endif

      call write_attr(grid,fid,'aijmm','reduction','max')
      call write_attr(grid,fid,'aijmm','split_dim',3)
      call defvar(grid,fid,scale_ijmm,'scale_aijmm(kaijmm)')
      call defvar(grid,fid,name_ijmm,'sname_aijmm(sname_strlen,kaijmm)')
      call defvar_cdl(grid,fid,cdl_ijmm,
     &     'cdl_aijmm(cdl_strlen,kcdl_aijmm)')

      call write_attr(grid,fid,'aijl','reduction','sum')
      call write_attr(grid,fid,'aijl','split_dim',4)
      call defvar(grid,fid,ia_ijl,'ia_aijl(kaijl)')
      call defvar(grid,fid,scale_ijl,'scale_aijl(kaijl)')
      call defvar(grid,fid,denom_ijl,'denom_aijl(kaijl)')
      call defvar(grid,fid,name_ijl,'sname_aijl(sname_strlen,kaijl)')
      call defvar_cdl(grid,fid,cdl_ijl,
     &     'cdl_aijl(cdl_strlen,kcdl_aijl)')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_ijl_latlon,
     &     'cdl_aijl_latlon(cdl_strlen,kcdl_aijl_latlon)')
#endif

#ifndef CUBED_SPHERE
#ifndef SCM
      call write_attr(grid,fid,'aijk','reduction','sum')
      call write_attr(grid,fid,'aijk','split_dim',4)
      call defvar(grid,fid,ia_ijk,'ia_aijk(kaijk)')
      call defvar(grid,fid,scale_ijk,'scale_aijk(kaijk)')
      call defvar(grid,fid,denom_ijk,'denom_aijk(kaijk)')
      call defvar(grid,fid,name_ijk,'sname_aijk(sname_strlen,kaijk)')
      call defvar_cdl(grid,fid,cdl_ijk,
     &     'cdl_aijk(cdl_strlen,kcdl_aijk)')
#endif
#endif

      call write_attr(grid,fid,'adiurn','reduction','sum')
      call write_attr(grid,fid,'adiurn','split_dim',1)
      call defvar(grid,fid,int_dummy,'ntime_adiurn')
      call write_attr(grid,fid,'ntime_adiurn','reduction','sum')
      call defvar(grid,fid,denom_dd,'denom_adiurn(ndiuvar)')
      call defvar(grid,fid,scale_dd,'scale_adiurn(ndiuvar)')
      call defvar(grid,fid,name_dd,'sname_adiurn(sname_strlen,ndiuvar)')
      call defvar_cdl(grid,fid,cdl_dd,
     &     'cdl_adiurn(cdl_strlen,kcdl_adiurn)')
#ifdef USE_HDIURN
      call write_attr(grid,fid,'hdiurn','split_dim',1)
      call defvar(grid,fid,int_dummy,'ntime_hdiurn')
      call defvar(grid,fid,denom_dd,'denom_hdiurn(ndiuvar)')
      call defvar(grid,fid,scale_dd,'scale_hdiurn(ndiuvar)')
      call defvar(grid,fid,name_dd,'sname_hdiurn(sname_strlen,ndiuvar)')
      call defvar_cdl(grid,fid,cdl_hd,
     &     'cdl_hdiurn(cdl_strlen,kcdl_hdiurn)')
#endif

      if(isccp_diags.eq.1) then
      call write_attr(grid,fid,'aisccp','reduction','sum')
      call defvar(grid,fid,wisccp,'wisccp(nisccp)')
      call write_attr(grid,fid,'wisccp','reduction','sum')
      call defvar(grid,fid,isccp_press,'isccp_press(npres)')
      call defvar(grid,fid,isccp_tau,'isccp_tau(ntau)')
      call defvar(grid,fid,isccp_late,'isccp_late(nisccp_plus_1)')
      endif

      call def_meta_rvracc(fid)

#ifdef TRACERS_ON
      call def_meta_trdiag(fid)
#endif

#ifndef SCM
      call write_attr(grid,fid,'speca','reduction','sum')
      call defvar(grid,fid,lmax_speca,'lmax_speca(nlspeca)')
      call write_attr(grid,fid,'atpe','reduction','sum')
#endif

      return
      end subroutine def_meta_atmacc

      subroutine write_meta_atmacc(fid)
!@sum  write_meta_atmacc write atm accumulation metadata to file
!@auth M. Kelley
      use model_com, only : nday,idacc,jyear0,jmon0
      use diag_com, only :
     &     ia_j,ia_jl,ia_ij,ia_ijl,ia_con,ia_ijk,
     &     name_j,name_reg,sname_jl,name_ij,name_ijl,name_dd,
     &     name_consrv,name_ijk,nisccp,ntau,
     &     cdl_j,cdl_reg,cdl_jl,
     &     cdl_ij,cdl_ijl,cdl_ij_latlon,cdl_ijl_latlon,
     &     cdl_dd,cdl_hd,cdl_consrv,cdl_ijk,
     &     hemis_j,hemis_jl,vmean_jl,hemis_consrv,
     &     hemis_ij,
     &     scale_j,scale_jl,scale_ij,scale_ijl,scale_dd,scale_con,
     &     scale_ijk,
     &     iden_j,iden_reg,denom_jl,denom_ij,denom_ijl,denom_dd,
     &     denom_ijk,
     &     lm,ia_12hr,
     &     isccp_diags,isccp_press,isccp_tau,isccp_late,wisccp,
     &     scale_ijmm,name_ijmm,cdl_ijmm,
     &     write_regions
#ifndef SCM
      use gc_com, only : kagc,ia_gc,sname_gc,cdl_gc,hemis_gc,
     &     vmean_gc,scale_gc,denom_gc,lmax_speca
#endif
      use geom, only : axyp
#ifdef CUBED_SPHERE
      use geom, only : lon2d_dg,lat2d_dg,lonbds,latbds
#endif
      use domain_decomp_atm, only : grid
      use pario, only : write_data,write_dist_data
      use cdl_mod, only : write_cdl
      implicit none
      integer fid   !@var fid unit number of read/write
      integer :: ntime_dd,ntime_hd
      real*8, allocatable :: tmpArr(:)
      real*8 :: r8time

#ifdef CUBED_SPHERE
      call write_dist_data(grid,fid,'lon',lon2d_dg)
      call write_dist_data(grid,fid,'lat',lat2d_dg)
      call write_dist_data(grid,fid,'lonbds',lonbds,jdim=3)
      call write_dist_data(grid,fid,'latbds',latbds,jdim=3)
#endif

      call write_dist_data(grid,fid,'axyp',axyp)

      r8time = real(jyear0,kind=8)+(real(jmon0,kind=8)-.5)/12d0
      call write_data(grid,fid,'time',r8time)

#ifndef SCM
      call write_data(grid,fid,'hemis_aj',hemis_j)
#endif
      call write_data(grid,fid,'ia_aj',ia_j)
      call write_data(grid,fid,'scale_aj',scale_j)
      call write_data(grid,fid,'denom_aj',iden_j)
      call write_data(grid,fid,'sname_aj',name_j)
      call write_cdl(grid,fid,'cdl_aj',cdl_j)

      if(write_regions) then
      call write_data(grid,fid,'ia_areg',ia_j)
      call write_data(grid,fid,'scale_areg',scale_j)
      call write_data(grid,fid,'denom_areg',iden_reg)
      call write_data(grid,fid,'sname_areg',name_reg)
      call write_cdl(grid,fid,'cdl_areg',cdl_reg)
      endif

#ifndef SCM
      call write_data(grid,fid,'hemis_consrv',hemis_consrv)
#endif
      call write_data(grid,fid,'ia_consrv',ia_con)
      call write_data(grid,fid,'scale_consrv',scale_con)
      call write_data(grid,fid,'sname_consrv',name_consrv)
      call write_cdl(grid,fid,'cdl_consrv',cdl_consrv)

#ifndef SCM
      call write_data(grid,fid,'hemis_ajl',hemis_jl)
      call write_data(grid,fid,'vmean_ajl',vmean_jl)
#endif
      call write_data(grid,fid,'ia_ajl',ia_jl)
      call write_data(grid,fid,'scale_ajl',scale_jl)
      call write_data(grid,fid,'denom_ajl',denom_jl)
      call write_data(grid,fid,'sname_ajl',sname_jl)
      call write_cdl(grid,fid,'cdl_ajl',cdl_jl)

#ifndef SCM
      call write_data(grid,fid,'hemis_agc',hemis_gc)
      call write_data(grid,fid,'vmean_agc',vmean_gc)
      call write_data(grid,fid,'ia_agc',ia_gc(1:kagc))
      call write_data(grid,fid,'scale_agc',scale_gc(1:kagc))
      call write_data(grid,fid,'denom_agc',denom_gc(1:kagc))
      call write_data(grid,fid,'sname_agc',sname_gc(1:kagc))
      call write_cdl(grid,fid,'cdl_agc',cdl_gc)
#endif

#ifndef SCM
      call write_data(grid,fid,'hemis_aij',hemis_ij)
#endif
      call write_data(grid,fid,'ia_aij',ia_ij)
      call write_data(grid,fid,'scale_aij',scale_ij)
      call write_data(grid,fid,'denom_aij',denom_ij)
      call write_data(grid,fid,'sname_aij',name_ij)
      call write_cdl(grid,fid,'cdl_aij',cdl_ij)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_aij_latlon',cdl_ij_latlon)
#endif

      call write_data(grid,fid,'scale_aijmm',scale_ijmm)
      call write_data(grid,fid,'sname_aijmm',name_ijmm)
      call write_cdl(grid,fid,'cdl_aijmm',cdl_ijmm)

      call write_data(grid,fid,'ia_aijl',ia_ijl)
      call write_data(grid,fid,'scale_aijl',scale_ijl)
      call write_data(grid,fid,'denom_aijl',denom_ijl)
      call write_data(grid,fid,'sname_aijl',name_ijl)
      call write_cdl(grid,fid,'cdl_aijl',cdl_ijl)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_aijl_latlon',cdl_ijl_latlon)
#endif

#ifndef CUBED_SPHERE
#ifndef SCM
      call write_data(grid,fid,'ia_aijk',ia_ijk)
      call write_data(grid,fid,'scale_aijk',scale_ijk)
      call write_data(grid,fid,'denom_aijk',denom_ijk)
      call write_data(grid,fid,'sname_aijk',name_ijk)
      call write_cdl(grid,fid,'cdl_aijk',cdl_ijk)
#endif
#endif

      ntime_dd = (idacc(ia_12hr)/2)*(nday/24)
      call write_data(grid,fid,'ntime_adiurn',ntime_dd)
      call write_data(grid,fid,'scale_adiurn',scale_dd)
      call write_data(grid,fid,'denom_adiurn',denom_dd)
      call write_data(grid,fid,'sname_adiurn',name_dd)
      call write_cdl(grid,fid,'cdl_adiurn',cdl_dd)

#ifdef USE_HDIURN
      ntime_hd = nday/24
      call write_data(grid,fid,'ntime_hdiurn',ntime_hd)
      call write_data(grid,fid,'scale_hdiurn',scale_dd)
      call write_data(grid,fid,'denom_hdiurn',denom_dd)
      call write_data(grid,fid,'sname_hdiurn',name_dd)
      call write_cdl(grid,fid,'cdl_hdiurn',cdl_hd)
#endif

      if(isccp_diags.eq.1) then
      call write_data(grid,fid,'isccp_press',(isccp_press))
! tmpArr is a workaround for pnetcdf under gfortran.  Apparently
! some issue with arrays with the PARAMETER attribute
      allocate(tmpArr(ntau))
      tmpArr = isccp_tau
      call write_data(grid,fid,'isccp_tau',tmpArr)
      deallocate(tmpArr)
      allocate(tmpArr(nisccp+1))
      tmpArr = isccp_late
      call write_data(grid,fid,'isccp_late',tmpArr)
      call write_data(grid,fid,'wisccp',wisccp)
      endif

      call write_meta_rvracc(fid)

#ifdef TRACERS_ON
      call write_meta_trdiag(fid)
#endif

#ifndef SCM
      call write_data(grid,fid,'lmax_speca',lmax_speca)
#endif

      return
      end subroutine write_meta_atmacc

      subroutine set_ioptrs_atmacc_default
c point i/o pointers for diagnostic accumlations to the
c instances of the arrays used during normal operation.
      use diag_com
#ifndef SCM
      use gc_com
#endif
      implicit none
      aj_ioptr     => aj
      areg_ioptr   => areg
#ifndef SCM
      agc_ioptr    => agc
      speca_ioptr  => speca
      atpe_ioptr   => atpe
#endif
      return
      end subroutine set_ioptrs_atmacc_default

      subroutine set_ioptrs_atmacc_extended
c point i/o pointers for diagnostic accumlations to the
c instances of the arrays containing derived outputs
      use diag_com
#ifndef SCM
      use gc_com
#endif
      implicit none
      aj_ioptr     => aj_out
      areg_ioptr   => areg_out
#ifndef SCM
      agc_ioptr    => agc_out
      speca_ioptr  => speca_out
      atpe_ioptr   => atpe_out
#endif
      return
      end subroutine set_ioptrs_atmacc_extended
