#include "rundeck_opts.h"

#ifndef SWFIX_20151201
#define SWFIX_20151201
#endif

      MODULE RADPAR
!@sum radiation module based originally on rad00b.radcode1.F
!@auth A. Lacis/V. Oinas/R. Ruedy
#ifndef USE_RAD_OFFLINE
      use constant, only: pO2, avog, mair, grav, loschmidt_constant
      use atm_com, only : lm_req
      use resolution, only : lm_gcm=>lm
#endif
#ifdef HEALY_LM_DIAGS
      USE RESOLUTION, only : JM_DIAG=>jm
#endif
      IMPLICIT NONE

C--------------------------------------------------
C     Grid parameters: Vertical resolution/profiles
C--------------------------------------------------

!@var LX max.number of vertical layers of the radiation (1D)-model
!@+
!@+   The Radiation Model can accomodate  arbitrary vertical resolution,
!@+              the number of layers may be time or location dependent,
!@+              but it cannot exceed LX.
#ifndef USE_RAD_OFFLINE
!@+   The GCM uses LM_REQ radiative equilibrium layers on top of the LM
!@+   atmospheric layers
      INTEGER, PARAMETER :: LX = lm_gcm+lm_req
#else
      INTEGER, PARAMETER :: LX = 57
#endif
!     optional repartitioning of gases - OFFLINE use only
!@var MRELAY if not 0, gases/aerosols are repartitioned to new layering
!@var KEEP10 if =10 N2 is kept, not repartitioned    (only if MRELAY>0)
!@+           n=1-9 N2 not repartitioned and replaces gas n
!@+         n=11-19 N2 not repartitioned and added to gas n-10
!@var NO3COL if >0 ozone is rescaled before repartitioning if MRELAY>0
!@var RO3COL = rescaled column amount of O3 if NO3COL>0   (if MRELAY>0)
      INTEGER :: MRELAY=0, KEEP10=0, NO3COL=0 ; REAL*8 :: RO3COL=1.

!     temperature profile within a layer: TLB,TLM,TLT bottom,mid,top T
!@var TLGRAD if >=0 tlt=tlm+dT*TLGRAD, tlb=tlm-dT*TLGRAD where
!@+        dT is chosen to try to minimize discontinuities if TLGRAD=1
!@+        if TLGRAD<0 tlt,tlm,tlb are all inputs        (OFFLINE use)
!@var PTLISO tlt=tlb=tlm above PTLISO mb independent of TLGRAD
      REAL*8 :: TLGRAD=1.             !     control param
      REAL*8 :: PTLISO=0d0            ! GCM control param

C-------------------------------------------
C     Grid parameters: Horizontal resolution
C-------------------------------------------

!@var MLAT46,MLON72 horizontal grid dimensions referred to in this model
!@+   The Radiation Model utilizes Data with 72x46 (lon,lat) resolution.
!@+               For GCM resolution other than 72x46, set JLAT and ILON
!@+               to appropriately Sample  (rather than interpolate) the
!@+               72x46 aerosol, ozone, cloud heterogeneity data sets
      INTEGER, PARAMETER ::  MLAT46=46,MLON72=72

!@var JNORTH latitude index defining northern hemisphere : jlat>jnorth
      INTEGER, PARAMETER ::  JNORTH=MLAT46/2

!     longitudes of box centers (degrees): -177.5,-172.5., ... ,177.5
!@var DLAT46 latitudes of box centers (degrees)
      REAL*8, PARAMETER :: DLAT46(46)=(/
     A    -90.,-86.,-82.,-78.,-74.,-70.,-66.,-62.,-58.,-54.,-50.,-46.,
     B    -42.,-38.,-34.,-30.,-26.,-22.,-18.,-14.,-10., -6., -2.,  2.,
     C      6., 10., 14., 18., 22., 26., 30., 34., 38., 42., 46., 50.,
     D     54., 58., 62., 66., 70., 74., 78., 82., 86., 90./)

C----------------
C     Input data               for the 1-d radiation
C----------------

!@var LASTVC if >= 0 picks sample atmosph. and ground data, OFFLINE only
      INTEGER :: LASTVC=-123456

!@var COSZ          cosine of zenith angle  (1)
      REAL*8 cosz
!@var JLAT,ILON     lat,lon index  w.r.to 72x46 lon-lat grid
!@var JGCM,IGCM     host GCM grid indices
!@var NL,L1         highest and lowest above ground layer
!@var LS1_loc       local tropopause level, used to limit H2O-scaling
      INTEGER   :: JLAT,ILON, NL,L1=1, LS1_loc ! Offline deflts L1=LS1_loc=1
      INTEGER   :: JGCM, IGCM
!@var JYEAR,JDAY    current year, Julian date
      INTEGER :: JYEAR=1980, JDAY=1

!@var PLB           layer pressure (mb) at bottom of layer
!@var HLB           height (km) at bottom of layer - currently NOT Used
!@var TLm           mean layer temperature (K)
!@var TLb,TLt       bottom,top layer temperature (K) - derived from TLm
!@+                                                    (unless TLGRAD<0)
!@var SHL,RHL       layer specific,relative humidity (1)
      REAL*8, dimension(LX+1) :: PLB,HLB,TLB
      REAL*8, dimension(LX)   :: TLT,TLM,SHL,RHL
!@var KEEPRH  if 0: find RH from SH, 1: find SH from RH, 2: keep both
      INTEGER :: KEEPRH=2

!@var ULGAS         current gas amounts, 13 types  (cm atm) (in getgas)
!@var TAUWC,TAUIC   opt.depth of water,ice cloud layer (1)
!@var SIZEWC,SIZEIC particle size of water,ice clouds (micron)
!@var CLDEPS        cloud heterogeneity; is computed using KCLDEP,EPSCON
      REAL*8 :: ULGAS(LX,13),TAUWC(LX),TAUIC(LX),SIZEWC(LX),SIZEIC(LX)
     *     ,CLDEPS(LX)
!@var       EPSCON  cldeps=EPSCON if KCLDEP=1
!@var       KCLDEP  KCLDEP=0->CLDEPS=0, 1->=EPSCON, 2->as is, 3,4->isccp
      REAL*8 :: EPSCON=0. ; INTEGER :: KCLDEP=4 ! control param

!@var KDELIQ Flag for dry(0) or wet(1) air deliquescence
      INTEGER :: KDELIQ(LX,4)
!@var KRHDTK if 1, RHlevel for deliquescence is temperature dependent
      INTEGER :: KRHDTK=1    !  control parameter

!@var SRBALB,SRXALB diffuse,direct surface albedo (1); see KEEPAL
      REAL*8 :: SRBALB(6),SRXALB(6),dalbsn ! prescr change in snowalbedo
!@var       KEEPAL  if 0, SRBALB,SRXALB are computed in SET/GETSUR
      INTEGER :: KEEPAL=0       ! control param
!@dbparm    KSIALB  sea ice albedo computation flag: 0=Hansen 1=Lacis
      INTEGER :: KSIALB=0
!@var PVT           frac. of surf.type (bareWhite+veg*8+bareDark+ocn)(1)
!@var AGESN 1-3     age of snow    (over soil,oice,land ice) (days)
!@var SNOWLI  amount of snow (over land ice)   (kg/m^2)
!@var SNOWD  amount of snow (over soil)   (m)
!@var SNOWOI        amount of snow (over ocean/lake ice)  (kg/m^2)
!@var WEARTH        soil wetness (1)
!@var WMAG          wind speed (m/s)
!@var POCEAN        fraction of box covered by ocean or lake  (1)
!@var PLAKE         fraction of box covered by lake           (1)
!@var PEARTH        fraction of box covered by soil           (1)
!@var POICE         fraction of box covered by ocean/lakeice  (1)
!@var PLICE         fraction of box covered by glacial ice    (1)
!@var TGO           top layer water temperature (K) of ocean/lake
!@var TGE,TGOI,TGLI top layer ground temperature (K) soil,seaice,landice
!@var TSL           surface air temperature (K)
      REAL*8 PVT(12),AGESN(3),SNOWD(2),SNOWOI,SNOWLI,WEARTH,WMAG,POCEAN
     *     ,PEARTH,POICE,PLICE,PLAKE,TGO,TGE,TGOI,TGLI,TSL
!@var KZSNOW        =1 for snow/ice albedo zenith angle dependence
      INTEGER :: KZSNOW=1
!     Additional info for Schramm/Schmidt/Hansen sea ice albedo KSIALB=0
!@var ZSNWOI        depth of snow over ocean ice (m)
!@var zoice         depth of ocean ice (m)
!@var zmp           depth of melt pond (m)
!@var fmp           fraction of melt pond area (1)
!@var zlake         lake depth (m)
!@var flags         true if snow is wet
!@var snow_frac(2)  fraction of snow over bare(1),vegetated(2) soil (1)
!@var snoage_fac_max  max snow age reducing-factor for sea ice albedo
      REAL*8 :: zsnwoi,zoice,zmp,fmp,zlake,snow_frac(2)
      REAL*8 :: snoage_fac_max=.5d0

!@var ITRMAX maximum number of optional tracers
      INTEGER, PARAMETER :: ITRMAX=150
!@var TRACER array to add up to ITRMAX additional aerosol species
      REAL*8    :: TRACER(LX,ITRMAX)
!@var FSTOPX,FTTOPX switches on/off aerosol for diagnostics (solar,thermal component)
!@var FSTASC,FTTASC scales optional aerosols (solar,thermal component)
      REAL*8    :: FSTOPX(ITRMAX),FTTOPX(ITRMAX)
!@var skip_AOD_in_rad If true, no optical depth calculations in RADIATION.f
      logical :: skip_AOD_in_rad
!@var chem_IN column variable for importing ozone(1) and methane(2)
!@+   fields from rest of model
!@var use_tracer_chem:set U0GAS(L, )=chem_IN( ,L), L=L1,use_tracer_chem( )
!@var GCCco2_IN column variable for importing CO2 and use_tracer_GCCco2 variable
#ifdef GCC_COUPLE_RAD
      REAL*8 :: GCCco2_IN(LX)
      INTEGER :: use_tracer_GCCco2
#endif
      REAL*8 :: chem_IN(2,LX)
      INTEGER :: use_tracer_chem(2),use_o3_ref=0
      LOGICAL*4 :: flags
!@var LOC_CHL local chlorophyll value (unit?) for albedo calculation (optional)
      REAL*8    :: LOC_CHL
#ifdef HEALY_LM_DIAGS
      REAL*8 :: VTAULAT(JM_DIAG)
#endif

      logical :: set_gases_internally = .true.,
     &           set_aerosols_internally = .true.

!@var U0GAS   reference gas amounts, 13 types  (cm atm)      (in setgas)
C     array with local and global entries: repeat this section in driver
      REAL*8 U0GAS(LX,13)
C     end of section to be repeated in driver (needed for 'copyin')

C--------------------------------------------------------
C     Output data     (from RCOMPX)  grid point dependent
C--------------------------------------------------------

!@var TRDFLB,TRUFLB,TRNFLB  Thrml down,up,net Flux at Layr Bottom (W/m2)
!@var SRDFLB,SRUFLB,SRNFLB  Solar down,up,net Flux at Layr Bottom (W/m2)
!@var TRFCRL,SRFHRL         layer LW Cooling Rate,SW Heating Rate (W/m2)
!@var SR.VIS,SR.NIR         SW fluxes in vis,near-IR domain       (W/m2)
!@var PLA...,ALB...         planetary and surface albedos            (1)
!@var TR...W,WINDZF            fluxes in the window region        (W/m2)
!@var BTEMPW,WINDZT      Brightness temperature in the window region (K)
!@var SK...,SRK...       Spectral breakdown of fluxes/heat.rates  (W/m2)
!@var FSRNFG,FTRUFG      surface type fractions of SW,LW fluxes   (W/m2)
!@var DTRUFG               not used                               (W/m2)
!sl!@var FTAUSL,TAUSL,...  surface layer computations commented out: !sl
!@var LBOTCL,LTOPCL  bottom and top cloud level (lbot < ltop)
!@var chem_out column variable for exporting radiation code quantities
!@    1=Ozone, 2=aerosol ext, 3=N2O, 4=CH4,5=CFC11+CFC12
!@var CO2outCol column CO2 export [mole mole-1] for SUBDD
!@var aesqex saves extinction aerosol optical thickness
!@var aesqsc saves scattering aerosol optical thickness
!@var aesqcb saves aerosol scattering asymmetry factor
!@var aesqex_dry saves dry extinction aerosol optical thickness
!@var aesqsc_dry saves dry scattering aerosol optical thickness
!@var aesqcb_dry saves dry aerosol scattering asymmetry factor

      REAL*8 TRDFLB(LX+1),TRUFLB(LX+1),TRNFLB(LX+1), TRFCRL(LX)
      REAL*8 SRDFLB(LX+1),SRUFLB(LX+1),SRNFLB(LX+1), SRFHRL(LX)
!@var GCCco2_out column CO2 for exporting
#ifdef GCC_COUPLE_RAD
      REAL*8 :: GCCco2_out(LX)=0d0
#endif
      REAL*8 :: chem_out(LX,5)=0d0
      REAL*8 :: CO2outCol(LX)=0.d0
      REAL*8 SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR,
     *       SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,
     *       SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR
      REAL*8 TRDFGW,TRUFGW,TRUFTW,BTEMPW,SRXVIS,SRXNIR
      REAL*8 WINDZF(3),WINDZT(3),TOTLZF(3),TOTLZT(3)
      REAL*8 SRKINC(16),SRKALB(16),SRKGAX(16,4),SRKGAD(16,4)
      REAL*8, dimension(LX,17) ::   SKFHRL
      REAL*8, dimension(LX+1,17) :: SKDFLB,SKUFLB,SKNFLB
      REAL*8 FSRNFG(4),FTRUFG(4),DTRUFG(4) ! ,SRXATM(4)
!sl   REAL*8 FTAUSL(33),TAUSL(33)             ! surf.layer input data
!nu  K      ,TRDFSL,TRUFSL,TRSLCR,SRSLHR,TRSLWV  !nu = not (yet) used
!sl  K      ,TRSLTS,TRSLTG,TRSLBS
      REAL*8 aesqex(lx,6,itrmax),aesqsc(lx,6,itrmax),aesqcb(lx,6,itrmax)
      REAL*8 aesqex_dry(lx,6,itrmax),aesqsc_dry(lx,6,itrmax),
     &       aesqcb_dry(lx,6,itrmax)
      INTEGER :: LBOTCL,LTOPCL

C----------------   scratch pad for temporary arrays that are passed to
C     Work arrays   other routines while working on a lat/lon point;
C----------------   but with openMP, each cpu needs its own copy !!

      real*8, dimension(LX,6,8) ::
     &  nintaerext,nintaersca,nintaerasy
      REAL*8, dimension(LX,6) ::
     *     SRAEXT,SRASCT,SRAGCB,SRBEXT,SRBSCT,SRBGCB,
     *     SRDEXT,SRDSCT,SRDGCB,SRVEXT,SRVSCT,SRVGCB,
     *     SRCEXT,SRCSCT,SRCGCB,SRCPI0
      REAL*8, dimension(LX+1,6) :: DBLEXT,DBLSCT,DBLGCB,DBLPI0
      REAL*8, dimension(LX ,33) :: TRTAUK,TRGXLK
     *        ,TRCALK,TRAALK,TRBALK,TRDALK,TRVALK
      REAL*8 DFLB(LX+1,33),UFLB(LX+1,33)
      REAL*8, dimension(33) :: TRCTCA,DFSL,UFSL,TXCTPG,TSCTPG
     *     ,TGCTPG,AVH2S,TRGALB,BGFEMT,BGFEMD
      REAL*8, dimension(LX) :: PL,DPL,O2FHRL,SRAXNL,SRASNL,SRAGNL,O2FHRB
      REAL*8 BXA(7),PRNB(6,4),PRNX(6,4),Q55H2S
     *     ,QVH2S(6),SVH2S(6),GVH2S(6),XTRU(LX,4),XTRD(LX,4)
     *     ,DXAERU(LX,4,4,LX+4),DXAERD(LX,4,4,LX+4)
      INTEGER IP24C9(LX)
C**** local except for special radiative aerosol diagnostics aadiag

      REAL*8 ::  SRCQPI(6,15),TRCQPI(33,15)       !??? to setcld/getcld
                 !  Temp data used by WRITER, WRITET
      REAL*8  :: TRAQAB(33,11),TRBQAB(33,10),TRCQAB(33,15),TRDQAB(33,25)
      REAL*8  :: AMP_TAB_SPEC(33,ITRMAX)
      INTEGER :: NORDER(16),NMWAVA(16),NMWAVB(16)

C------------------------------------------
C     Reference data, Tables, Climatologies
C------------------------------------------

      REAL*8, PARAMETER :: DKS0(16)=(/
     *           .010, .030, .040, .040, .040, .002, .004, .013,
     +           .002, .003, .003, .072, .200, .480, .050, .011/)

      INTEGER ::  NKSLAM=14
      INTEGER,parameter :: KSLAM(16)=(/1,1,2,2,5,5,5,5,1,1,1,3,4,6,6,1/)

      REAL*8 ::                 !   Model parameters generated by RCOMP1
     H              HLB0(LX+1),PLB0(LX+1),TLM0(LX),U0GAS3(LX)
     A             ,TKPFW(630),TKPFT(900),AO3(460)
     D             ,FPXCO2(LX),FPXOZO(LX) !nu ,PIAERO(10)
!     E ,QXDUST(6,8),QSDUST(6,8),QCDUST(6,8),ATDUST(33,8),QDST55(8) !?DST   !ron
     D             ,TRAX(LX,33,5),DBLN(30),TCLMIN

      logical :: dust_optics_initialized=.false.
      real*8, dimension(:,:), allocatable :: QXDUST, QSDUST, QCDUST, !ron
     *     ATDUST                                                    !ron
      real*8, dimension(  :), allocatable :: QDST55                  !ron
      real*8, dimension(:), allocatable :: taucon_dust

!@dbparam planck_tmin, planck_tmax temperature range for Planck function
!@+       lookup table.  If the requested tmin is less than the default
!@+       value of 124 K, the lookup table is extrapolated at startup to
!@+       cover the requested range (same for tmax exceeding 373 K).
      integer :: planck_tmin=1, planck_tmax=800

!@var transmission_corrections whether to apply correction factors
!@+   to longwave transmission
      logical :: transmission_corrections
C            RADDAT_TR_SGP_TABLES          read from  radfile1, radfile2
      INTEGER, PARAMETER :: NGUX=1024, NTX=8, NPX=19
      REAL*8, dimension(NGUX,NTX,NPX) :: TAUTBL,TAUWV0,TAUCD0,TAUO30
      REAL*8  H2O(100),FCO2(100)
      REAL*8, DIMENSION(1:800,33) :: PLANCK
      REAL*8  XKCFC(12,8,17:20),ULOX(19,16),DUX(19,16), XTFAC(11,9)
! Correction-factor lookup-table sizes
! NLCF  : number of layers in ref. atm. used to compute the table
! NWVCF : number of H2O vapor column amounts
! NUCF  : number of column amounts for absorbers other than H2O
! NRCF  : number of principal absorber regions (H2O, CO2, O3)
      INTEGER, PARAMETER, PRIVATE :: NLCF=43,NWVCF=9,NUCF=7,NRCF=3
      REAL*8, dimension(NLCF,NRCF) :: XTU0,XTD0
      REAL*8, dimension(NLCF,NWVCF,NRCF) :: XTRUP,XTRDN,DXUP13,DXDN13
      REAL*8, dimension(NLCF,NWVCF,NUCF,NRCF) ::
     *        DXUP2,DXUP3,DXUP6,DXUP7,DXUP8,DXUP9
     *       ,DXDN2,DXDN3,DXDN6,DXDN7,DXDN8,DXDN9
c---------------------------------------------------------------------
C         Default h2o continuum is Ma 2000.  Other options: Ma 2004
C         Roberts, MT_CKD model (Mlawer/Tobin_Clough/Kneizys/Davies)
C---------------------------------------------------------------------
      REAL*8 H2OCN8(33,8,14),H2OCF8(33,8,5)

C            RADDAT_AERCLD_MIEPAR          read from            radfile3
      REAL*8 ::
     A              SRAQEX( 6,11),SRAQSC( 6,11),SRAQCB( 6,11),Q55A11(11)
     B             ,TRAQEX(33,11),TRAQSC(33,11),TRAQCB(33,11),REFA11(11)
     C             ,SRBQEX( 6,10),SRBQSC( 6,10),SRBQCB( 6,10),Q55B10(10)
     D             ,TRBQEX(33,10),TRBQSC(33,10),TRBQCB(33,10),REFB10(10)
     E             ,SRCQEX( 6,15),SRCQSC( 6,15),SRCQCB( 6,15),Q55C15(15)
     F             ,TRCQEX(33,15),TRCQSC(33,15),TRCQCB(33,15),REFC15(15)
     G             ,TRCQAL(33,15),VEFC15(15)   ,VEFA11(   11),VEFB10(10)
     H             ,SRDQEX( 6,25),SRDQSC( 6,25),SRDQCB( 6,25),Q55D25(25)
     .             ,YRDQEX( 6,25),YRDQSC( 6,25),YRDQCB( 6,25),Y55D25(25)
     I             ,TRDQEX(33,25),TRDQSC(33,25),TRDQCB(33,25),REFD25(25)
     J             ,TRDQAL(33,25),VEFD25(25)
     K         ,SRVQEX( 6,20,6),SRVQSC( 6,20,6),SRVQCB( 6,20,6)
     L         ,TRVQEX(33,20,6),TRVQSC(33,20,6),TRVQCB(33,20,6)
     M         ,TRVQAL(33,20,6),Q55V20(20,6),REFV20(20,6),VEFV20(20,6)
     N         ,SRUQEX( 6,120),SRUQSC( 6,120),SRUQCB( 6,120),Q55U22(120)
     O         ,TRUQEX(33,120),TRUQSC(33,120),TRUQCB(33,120),REFU22(120)
     P         ,TRUQAL(33,120),VEFU22(120),TRSQAL(33,25),VEFS25(25)
     Q             ,SRSQEX( 6,25),SRSQSC( 6,25),SRSQCB( 6,25),Q55S25(25)
     R             ,TRSQEX(33,25),TRSQSC(33,25),TRSQCB(33,25),REFS25(25)

      REAL*8    SRQV( 6,20),SRSV( 6,20),SRGV( 6,20),Q55V(   20),REFV(20)
      REAL*8    TRQV(33,20),TRSV(33,20),TRGV(33,20),TRAV(33,20),VEFV(20)
      EQUIVALENCE (SRVQEX(1,1,6),SRQV(1,1)), (SRVQSC(1,1,6),SRSV(1,1))
      EQUIVALENCE (SRVQCB(1,1,6),SRGV(1,1)),   (Q55V20(1,6),Q55V(1))
      EQUIVALENCE (TRVQEX(1,1,6),TRQV(1,1)), (TRVQSC(1,1,6),TRSV(1,1))
      EQUIVALENCE (TRVQCB(1,1,6),TRGV(1,1)), (TRVQAL(1,1,6),TRAV(1,1))
      EQUIVALENCE   (REFV20(1,6),REFV(1)),     (VEFV20(1,6),VEFV(1))

C            RADDAT_CLDCOR_TRSCAT           read from           radfileE
      REAL*8 :: RIJTPG(6,49,17,21),FDXTPG(3,49,17,21),FEMTPG(3,49,17,21)

!@var ppmv_to_cm_at_stp Conversion factor for conversion from PPMV to cm at
!                       STP. Also needs an additional factor dP for the
!                       conversion.
      REAL*8, PARAMETER :: ppmv_to_cm_at_stp = 1.0D-05*avog/
     *      (grav*mair*loschmidt_constant)
!@var h2o_mmr_to_cm_at_stp Conversion factor for conversion from mass
!                          mixing ratio to cm at STP for water vapor.
!                          Also needs an additional factor dP for the
!                          conversion.
      REAL*8, PARAMETER :: h2o_mmr_to_cm_at_stp = ppmv_to_cm_at_stp*
     *      1.0D+06*mair/18.0153D0


C--------------------------------------   This also should be moved out
C     History files (+ control options)   of RADPAR, which should just
C--------------------------------------   have to handle 1 point in time

!     -------------------------------------------------------i/o control
!@var MADxxx  Model Add-on Data of Extended Climatology Enable Parameter
!@+   ------   if 0   input process is skipped
!@+ 2 MADAER   =  1   Reads  Aerosol tropospheric climatology
!@+ 3 MADDST   =  1   Reads  Dust-windblown mineral climatology   RFILE6
!@+ 4 MADVOL   =  1   Reads  Volcanic 1950-00 aerosol climatology RFILE7
!@+ 5 MADEPS   =  1   Reads  Epsilon cloud heterogeneity data     RFILE8
!@+ 6 MADLUV   =  1   Reads  Lean format Spectral Solar Irrad.    RFILE9
!@+   MADGHG   =  1          Enables UPDGHG update. MADGHG=0: no update
!@+   MADSUR   =  1   Reads  Vegetation,Topography data    RFILEC,RFILED
!@+   MADBAK   if 1          Adds background aerosols
!@+   MADO2A   if > 0      call set/geto2a,  activating O2 solar heating
!     ------------------------------------------------------------------
      INTEGER :: MADO3M=1,MADAER=0,MADDST=0,MADVOL=0,MADEPS=0,MADLUV=1
      INTEGER :: MADGHG=1,MADSUR=0,MADBAK=0 ! MADSUR=1 for OFF-line use
      INTEGER :: MADO2A=1

!     ------------------------------------------------------time control
!@var KYEARx,KJDAYx if both are 0   : data are updated to current yr/day
!@+   -------------    only KJDAYx=0: data cycle through year KYEARx
!@+                    neither is 0 : yr/day=KYEARx/KJDAYx data are used
!@+   KYEARS,KJDAYS: Solar Trend
!@+   KYEARO,KJDAYO: Ozone Trend
!@+   KYEARD,KJDAYD: Dust Trend
!@+   KYEARE,KJDAYE: CldEps Trend
!@+   KYEARG,KJDAYG: GHG  Trend
!@+   KYEARR,KJDAYR: RVegeTrend (Ground Albedo)
!@+   KYEARV,KJDAYV: Volc.Aerosol Trend
!@+   KYEARA,KJDAYA: trop.Aerosol Trend
!     ------------------------------------------------------------------
      INTEGER ::                    KYEARS=0,KJDAYS=0, KYEARG=0,KJDAYG=0
     *          ,KYEARO=0,KJDAYO=0, KYEARA=0,KJDAYA=0, KYEARD=0,KJDAYD=0
     *          ,KYEARV=0,KJDAYV=0, KYEARE=0,KJDAYE=0, KYEARR=0,KJDAYR=0

      REAL*8, dimension(:,:,:), pointer :: o3jday,o3jref
#ifdef HIGH_FREQUENCY_O3_INPUT
      REAL*8, dimension(:,:,:), pointer :: o3jday_HF_modelLevels
#endif

!@var PLBA21 Vert. Layering for tropospheric aerosols (reference)
      REAL*8, PARAMETER :: PLBA20(21)=(/
     *  984.,964.,934.,884.,810.,710.,550.,390.,285.,210.,
     *  150.,110., 80., 55., 35., 20., 10., 3.,  1.,0.3,0.1/)
!@var PLBA09 Vert. Layering for tropospheric aerosols/dust (reference)
      REAL*8, PARAMETER :: PLBA09(10)=(/
     *  1010.,934.,854.,720.,550.,390.,255.,150., 70., 10./)
      real*8, dimension(:), pointer ::  plbaer => null()
      real*8, dimension(:,:,:,:), pointer :: A6JDAY => null()


C            RADMAD3_DUST_SEASONAL            (user SETDST)     radfile6
!      REAL*4 TDUST(72,46,9,8,12)                                   !ron
!      REAL*8 DDJDAY(9,8,72,46)                                     !ron

C            RADMAD4_VOLCAER_DECADAL          (user SETVOL)     radfile7
      INTEGER JVOLYI,JVOLYE,NVOLMON,NVOLLAT,NVOLK
      real*8, dimension(:), allocatable :: ELATVOL,HVOLKM
      real*8, dimension(:,:,:), allocatable :: VTauTJK ! (NVOLMON,NVOLLAT,NVOLK)
      real*8, dimension(:,:), allocatable :: VReffTJ   ! (NVOLMON,NVOLLAT)


C            RADMAD5_CLDEPS_3D_SEASONAL       (user SETCLD)     radfile8
      REAL*4 EPLMHC(72,46,12,4)
      REAL*8 EPLOW(72,46),EPMID(72,46),EPHIG(72,46),EPCOL(72,46)

C            RADMAD6_SOLARUV_DECADAL          (user SETSOL)     radfile9
!@var iy1S0,MS0X first year, max.number of months for S0 history
!@var icycs0  solar cycle in yrs used to extend S0 history before 2000
!@var icycs0f solar cycle in yrs used to extend S0 history after 2000
!@var KSOLAR controls which data are used: <0 Thekaekara, else Lean:
!@+          1: use monthly data, 2: use annual data, 0: constant data
!@+          9: use annual data from file but with Thekaekara bins
      INTEGER :: KSOLAR=2       ! MADLUV=KSOLAR=0 only possible OFF-line

      INTEGER, PARAMETER :: iy1S0=1882, MS0X=12*(1998-iy1S0+1)
      INTEGER, PARAMETER :: icycs0=11,  icycs0f=12
      INTEGER  iMS0X
      REAL*4 yr1S0,yr2S0
      real,    ALLOCATABLE, DIMENSION(:,:):: UV_SSI
      real,    ALLOCATABLE, DIMENSION(:)  :: TSI1,TSI2
      REAL*8 FS_SSI(190),W1_SSI(190)

      REAL*8 :: S00WM2=1366.2911d0, S0=1366.d0, RATLS0=1.

      REAL*8 :: WSOLAR(190),FSOLAR(190)

C***  alternate sources to get WSOLAR,FSOLAR:
      REAL*8, dimension(190) :: WS_SSI,DS_SSI,FR_SSI
#ifdef USE_RAD_OFFLINE
      common/lean1950/ WS_SSI,DS_SSI,FR_SSI ! for MADLUV=0 uses block data
#endif
      REAL*8, PARAMETER :: WTHEK(190)=(/        ! if KSOLAR<0
     *           .115,.120,.125,.130,.140,.150,.160,.170,.180,.190,.200,
     1 .210,.220,.225,.230,.235,.240,.245,.250,.255,.260,.265,.270,.275,
     2      .280,.285,.290,.295,.300,.305,.310,.315,.320,.325,.330,.335,
     3           .340,.345,.350,.355,.360,.365,.370,.375,.380,.385,.390,
     4           .395,.400,.405,.410,.415,.420,.425,.430,.435,.440,.445,
     5           .450,.455,.460,.465,.470,.475,.480,.485,.490,.495,.500,
     6           .505,.510,.515,.520,.525,.530,.535,.540,.545,.550,.555,
     7           .560,.565,.570,.575,.580,.585,.590,.595,.600,.605,.610,
     8           .620,.630,.640,.650,.660,.670,.680,.690,.700,.710,.720,
     9           .730,.740,.750,.760,.770,.780,.790,.800,.810,.820,.830,
     A .840,.850,.860,.870,.880,.890,.900,.910,.920,.930,.940,.950,.960,
     B 0.97,0.98,0.99,1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,
     C 1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00,2.10,2.20,
     D 2.30,2.40,2.50,2.60,2.70,2.80,2.90,3.00,3.10,3.20,3.30,3.40,3.50,
     E 3.60,3.70,3.80,3.90,4.00,4.10,4.20,4.30,4.40,4.50,4.60,4.70,4.80,
     F  4.9, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,11.0,12.0,13.0,14.0,15.00/)

      REAL*8, PARAMETER :: FTHEK(190)=(/
     *         .007,.900,.007,.007,.030,.070,.230,.630,1.25,2.71,10.7,
     1 22.9,57.5,64.9,66.7,59.3,63.0,72.3,70.4,104.,130.,185.,232.,204.,
     2    222.,315.,482.,584.,514.,603.,689.,764.,830.,975.,1059.,1081.,
     31074.,1069.,1093.,1083.,1068.,1132.,1181.,1157.,1120.,1098.,1098.,
     41189.,1429.,1644.,1751.,1774.,1747.,1693.,1639.,1663.,1810.,1922.,
     52006.,2057.,2066.,2048.,2033.,2044.,2074.,1976.,1950.,1960.,1942.,
     61920.,1882.,1833.,1833.,1852.,1842.,1818.,1783.,1754.,1725.,1720.,
     71695.,1705.,1712.,1719.,1715.,1712.,1700.,1682.,1666.,1647.,1635.,
     81602.,1570.,1544.,1511.,1486.,1456.,1427.,1402.,1389.,1344.,1314.,
     91290.,1260.,1235.,1211.,1185.,1159.,1134.,1109.,1085.,1060.,1036.,
     A1013.,990.,968.,947.,926.,908.,891.,880.,869.,858.,847.,837.,820.,
     B 803.,785.,767.,748.,668.,593.,535.,485.,438.,397.,358.,337.,312.,
     C 288.,267.,245.,223.,202.,180.,159.,142.,126.,114.,103., 90., 79.,
     D 69.0,62.0,55.0,48.0,43.0,39.0,35.0,31.0,26.0,22.6,19.2,16.6,14.6,
     E 13.5,12.3,11.1,10.3, 9.5,8.70,7.80,7.10,6.50,5.92,5.35,4.86,4.47,
     F  4.11,3.79,1.82,0.99,.585,.367,.241,.165,.117,.0851,.0634,.0481/)

!icb         RADMAD7_VEG_TOPOG          (user SETSUR)  radfileC,radfileD
!icb                 FVEG11(72,46,11),FOLGIZ(72,46,9)

C            RADMAD8_RELHUM_AERDATA     (user SETAER,SETREL)    radfileH
!nu   KRHAER(4) -1/0/1 flag to base aeros.sizes on 70%/0%/model rel.humi
!nu   INTEGER :: KRHAER(4)=(/1,1,1,1/) ! SO4,SSalt,NO3,OC
!@var KRHTRA(ITRMAX) 0/1 to make tracer aerosols rel.humid dependent
      INTEGER :: KRHTRA(ITRMAX)= 1
      REAL*8 ::
     A               SRHQEX(6,190,4),SRHQSC(6,190,4),SRHQCB( 6,190,4)
     B              ,TRHQAB(33,190,4),RHINFO(190,15,4)
     C   ,SRTQEX(6,190,ITRMAX),SRTQSC(6,190,ITRMAX),SRTQCB(6,190,ITRMAX)
     D   ,TRTQAB(33,190,ITRMAX),RTINFO(190,15,ITRMAX)

!new
!new  save TSOIL,TVEGE                  (not implemented)
!nu   DIMENSION PI0TRA(11)
!new  save FTRUFS,FTRUFV,DTRUFS,DTRUFV  (not implemented)

C     -----------------------
C     Ozone absorption tables
C     -----------------------
      REAL*8, PARAMETER ::        XWAVO3(226)=(/
     *            .2002,.2012,.2022,.2032,.2042,.2052,.2062,.2072,.2082,
     A.2092,.2102,.2112,.2122,.2132,.2142,.2152,.2162,.2172,.2182,.2192,
     B.2202,.2212,.2222,.2232,.2242,.2252,.2262,.2272,.2282,.2292,.2302,
     C.2312,.2322,.2332,.2342,.2352,.2362,.2372,.2382,.2392,.2400,.2402,
     D.2412,.2422,.2432,.2438,.2444,.2452,.2458,.2463,.2472,.2478,.2482,
     E.2490,.2492,.2500,.2508,.2519,.2527,.2539,.2543,.2553,.2562,.2566,
     F.2571,.2575,.2579,.2587,.2597,.2604,.2617,.2624,.2635,.2643,.2650,
     G.2654,.2662,.2669,.2675,.2682,.2692,.2695,.2702,.2712,.2718,.2722,
     H.2732,.2742,.2746,.2752,.2762,.2772,.2782,.2792,.2802,.2812,.2822,
     I.2830,.2842,.2852,.2862,.2872,.2882,.2892,.2902,.2912,.2922,.2932,
     J.2942,.2952,.2962,.2972,.2982,.2992,.2998,
     &            .3004,.3016,.3021,.3029,.3036,.3037,.3051,.3053,.3059,
     A.3061,.3066,.3075,.3077,.3083,.3085,.3092,.3098,.3100,.3104,.3106,
     B.3109,.3112,.3130,.3135,.3146,.3148,.3151,.3154,.3167,.3170,.3173,
     C.3176,.3190,.3194,.3199,.3200,.3209,.3210,.3216,.3220,.3223,.3226,
     D.3239,.3242,.3245,.3248,.3253,.3255,.3269,.3272,.3275,.3279,.3292,
     E.3295,.3299,.3303,.3309,.3312,.3328,.3332,.3334,.3338,.3357,.3365,
     F.3369,.3372,.3391,.3395,.3398,.3401,.3417,.3421,.3426,.3430,.3437,
     G.3439,.3451,.3455,.3460,.3463,.3466,.3472,.3481,.3485,.3489,.3493,
     H.3499,.3501,.3506,.3514,.3521,.3523,.3546,.3550,.3554,.3556,.3561,
     I.3567,.3572,.3573,.3588,.3594,.3599,.3600,.3604,.3606,.3639,.3647,
     J.3650,.3654,.3660/)
      REAL*8 ::  UVA(226)
      REAL*8, PARAMETER ::  FUVKO3(226)=(/
     *             8.3,  8.3,  8.1,  8.3,  8.6,  9.0,  9.7, 10.8, 11.7,
     A 13.0, 14.3, 16.0, 18.0, 20.6, 23.0, 26.1, 29.3, 32.6, 36.9, 40.8,
     B 46.9, 51.4, 56.7, 63.4, 69.1, 76.6, 84.0, 91.4, 99.9,110.0,118.0,
     C126.0,136.0,145.0,154.0,164.0,175.0,186.0,192.0,201.0,210.0,212.0,
     D221.0,230.0,239.0,248.0,250.0,259.0,264.0,264.0,273.0,277.0,275.0,
     E283.0,283.0,290.0,283.0,297.0,290.0,300.0,290.0,302.0,295.0,283.0,
     F293.0,290.0,286.0,297.0,281.0,280.0,271.0,275.0,254.0,264.0,250.0,
     G248.0,242.0,228.0,230.0,216.0,213.0,211.0,199.0,188.0,188.0,178.0,
     H169.0,153.0,155.0,148.0,136.0,127.0,117.0,108.0, 97.0, 88.7, 81.3,
     I 78.7, 67.9, 61.4, 54.3, 49.6, 43.1, 38.9, 34.6, 30.2, 27.5, 23.9,
     J 21.0, 18.6, 16.2, 14.2, 12.3, 10.7,  9.5,
     &            8.880,7.520,6.960,6.160,5.810,5.910,4.310,4.430,4.130,
     A4.310,4.020,3.330,3.390,3.060,3.100,2.830,2.400,2.490,2.330,2.320,
     B2.120,2.200,1.436,1.595,1.074,1.138,1.068,1.262,0.818,0.948,0.860,
     C1.001,0.543,0.763,0.665,0.781,0.382,0.406,0.373,0.608,0.484,0.601,
     D0.209,0.276,0.259,0.470,0.319,0.354,0.131,0.223,0.185,0.339,0.080,
     E0.093,0.079,0.184,0.139,0.214,0.053,0.074,0.068,0.152,0.038,0.070,
     F.0540000,.1030000,.0240000,.0382500,.0292500,.0550000,.0135000,
     G.0155250,.0127500,.0188250,.0167250,.0262500,.0115500,.0140250,
     H.0099750,.0115500,.0081000,.0104250,.0050100,.0057000,.0046650,
     I.0073425,.0051825,.0055275,.0040575,.0077700,.0048900,.0054600,
     J.0015375,.0017775,.0013275,.0014100,.0011550,.0023325,.0018825,
     K.0019650,.0009600,.0013650,.0011925,.0013200,.0008925,.0009825,
     L.0001350,.0006300,.0004500,.0006225,0.0/)

C     ------------------------------------------------------------------
C          NO2 Trace Gas Vertical Distribution and Concentration Profile
C     ------------------------------------------------------------------

      REAL*8, PARAMETER ::
     *     CMANO2(42)=(/            ! every 2 km starting at 0km
     1  8.66E-06,5.15E-06,2.85E-06,1.50E-06,9.89E-07,6.91E-07,7.17E-07,
     2  8.96E-07,3.67E-06,4.85E-06,5.82E-06,6.72E-06,7.77E-06,8.63E-06,
     3  8.77E-06,8.14E-06,6.91E-06,5.45E-06,4.00E-06,2.67E-06,1.60E-06,
     4  8.36E-07,3.81E-07,1.58E-07,6.35E-08,2.57E-08,1.03E-08,4.18E-09,
     5  1.66E-09,6.57E-10,2.58E-10,1.02E-10,4.11E-11,1.71E-11,7.73E-12,
     6  9.07E-12,4.63E-12,2.66E-12,1.73E-12,1.28E-12,1.02E-12,1.00E-30/)

C     ------------------------------------------------------------------
C     TRACE GAS REFERENCE AMOUNTS & DISTRIBUTIONS ARE DEFINED IN  SETGAS
C     ------------------------------------------------------------------

C-------------------------
C     Scaling/kill factors
C-------------------------

!@var FULGAS scales the various atmospheric constituents:
!@+         H2O CO2 O3 O2 NO2 N2O CH4 F11 F12 N2C CFC11 CFC12 SO2
!@+   Note: FULGAS(1) only acts in the stratosphere (unless LS1_loc=1)
      REAL*8 :: FULGAS(13) = (/    ! scales ULGAS

C      H2O CO2  O3  O2 NO2 N2O CH4 F11 F12 N2C CFC11+ CFC12+ SO2
C        1   2   3   4   5   6   7   8   9  10    11     12   13
     +   1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,   1.,    1.,  0./)
#ifdef ALTER_RADF_BY_LAT
!@var FULGAS_orig saves initial FULGAS values
      REAL*8, dimension(13) :: FULGAS_orig
#endif

!@var FGOLDH scales background aerosols for Glb Ocn Land Desert Haze
C                         GLOBAL  OCEAN   LAND  DESERT    HAZE
C     for setbak/getbak only   1      2      3       4       5
      REAL*8 :: FGOLDH(5)=(/ 1d0, .68d0, .32d0, 1.d-20, 1.d-20 /)

!@var FSxAER,FTxAER scales solar,thermal opt.depth for var. aerosols:
!@+     x =    T:total B:background A:atmClim  D:dust  V:volcanic
      REAL*8 :: FSTAER=1.,FSBAER=1.,FSAAER=1.,FSDAER=1.,FSVAER=1.
     *         ,FTTAER=1.,FTBAER=1.,FTAAER=1.,FTDAER=1.,FTVAER=1.

!@var FTAUC factor to control cloud optical depth in radiation calc.
!@+   =1 for full expression, =0 for clear sky calculation.
      REAL*8 :: FTAUC ! to be set in calling routine, thread-private ! deflt=1

!@var PIVMAX limits PI0 of volcanic aerosols
      REAL*8 :: PIVMAX=1.0
!@var ECLTRA,KCLDEM scales,enables full cloud scattering correction
      REAL*8 :: ECLTRA=1. ; INTEGER :: KCLDEM=1
!@var FCLDTR,FCLDSR scales opt.depth of clouds - not used (yet)
!@var FRAYLE        scales Rayleigh parameter
      REAL*8 ::   FCLDTR=1.,  FCLDSR=1.,  FRAYLE=1.

!@var KUVFAC,UVFACT,UVWAVL,KSNORM rescale UV spectral flux distribution
      INTEGER :: KUVFAC=0,  KSNORM=0  ! no rescaling
      REAL*8  :: UVWAVL(3)=(/0.295d0, 0.310d0, 0.366d0/)
      REAL*8  :: UVFACT(3)=(/0.98011d0, 0.99467d0, 0.99795d0/)

!@var SRCGSF Scaling Factors for Cloud Asymmetry Parameter for
!@+                            Water    Ice    MieIce
      REAL*8  ::  SRCGSF(3)=(/ 1.000,  1.000,  1.000/)

!@var TAUWC0,TAUIC0 lower limits for water/ice cloud opt.depths
      REAL*8 ::  TAUWC0=1d-3, TAUIC0=1d-3

!@var KFPCO2,KPFOZO if > 0 scale CO2,O3 vertical profile
      INTEGER :: KFPCO2=-1, KPFOZO=0

!@var KANORM,KCNORM if > 0 renormalize aerosols,cloud albedos
      INTEGER :: KANORM=0, KCNORM=0

!@var KWVCON        ON/OFF flag for water vapor continuum absorption
!@var KUFH2O,KUFCO2 H2O,CO2 column absorb.scaling
!@var KCSELF,KCFORN H2O_ContSelf-Broadening,CO2_ContForeign-Broadening
      INTEGER :: KWVCON=1, KUFH2O=1,  KUFCO2=1,  KCSELF=1,  KCFORN=1
!@var XCSELF,XCFORN scaling factors for Cont.Broadening (Deflt: Ma 2000)
      REAL*8 :: XCSELF=1. , XCFORN=1.

!@var ICE012 pick ice droplet type: 0 liquid, 1 ice non-spher, 2 ice Mie
      INTEGER :: ICE012=1

!@var VEFF0 effective volc. aerosol size distribution variance
      REAL*8  :: VEFF0=0.35d0,  REFF0=0.30d0      ! REFF0 not used

!@var NORMS0 if =1, Incident (TOA) Solar flux is normalized to equal S0
      INTEGER :: NORMS0=1

!@var fOnOff if =1 fully turns on SW long-path H2O absorption correction
      REAL*8 :: fOnOff=1. ! if =0. disables SW-H2O correction (tunable)

!@var KORDER,KWTRAB controls WRITER-output (Mie-scattering info)
      INTEGER :: KWTRAB=0, KORDER=0

C-----------------------------------------------------------------------
C      COMPOSITION & VERTICAL DISTRIBUTION FOR 5 SPECIFIED AEROSOL TYPES
C-----------------------------------------------------------------------
C TYPE
C    1   STRATOSPHERIC GLOBAL AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES
C    2    TROPOSPHERIC  OCEAN AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES
C    3    TROPOSPHERIC   LAND AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES
C    4    TROPOSPHERIC DESERT AEROSOL  A,B,C ARE  LOCAL AVERAGE VALUES
C    5    TROPOSPHERIC   HAZE AEROSOL  A,B,C ARE  LOCAL AVERAGE VALUES

C        1     2     3     4     5     6     7     8     9    10    11
C      ACID1 SSALT SLFT1 SLFT2 BSLT1 BSLT2 DUST1 DUST2 DUST3 CARB1 CARB2
      REAL*8, dimension(11,5) :: AGOLDH=reshape( (/
     1 .005,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,
     2   .0, .020, .010, .010, .005,   .0, .010,   .0,   .0, .005,   .0,
     3   .0,   .0,   .0, .020, .005,   .0, .010, .010,   .0,   .0, .015,
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, .020, .010,   .0,   .0,
     5   .0,   .0,   .0, .010,   .0,   .0,   .0,   .0,   .0,   .0, .005/
     *  ),(/11,5/) )
      REAL*8, dimension(11,5) :: BGOLDH=reshape( (/
     1 20.0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,
     2   .0, 1.00, 4.00, 1.00, 4.00, 1.00, 4.00,   .0,   .0, 1.00,   .0,
     3   .0,   .0,   .0, 0.00, 2.00,   .0, 4.00, 2.00,   .0,   .0, 0.00,
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, 2.00, 0.00,   .0,   .0,
     5   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0, 0.00/
     *  ),(/11,5/) )
      REAL*8, dimension(11,5) :: CGOLDH=reshape( (/
     1 3.00,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,
     2   .0, 1.00, 3.00, 2.00, 3.00, 1.00, 2.00,   .0,   .0, 1.00,   .0,
     3   .0,   .0,   .0, 1.00, 3.00,   .0, 1.00, 1.00,   .0,   .0, 1.00,
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, 1.00, 1.00,   .0,   .0,
     5   .0,   .0,   .0, 1.00,   .0,   .0,   .0,   .0,   .0,   .0, 1.00/
     *  ),(/11,5/) )

!nu   REAL*8, dimension(11) :: PI0VIS=(/
!nu         1          2          3          4          5          6
!nu       ACID1      SSALT      SLFT1      SLFT2      BSLT1      BSLT2
!nu  1   1.00000,   1.00000,   1.00000,   1.00000,   0.98929,   0.95609,
!nu
!nu         7          8          9         10         11
!nu       DUST1      DUST2      DUST3      CARB1      CARB2
!nu  2   0.91995,   0.78495,   0.63594,   0.31482,   0.47513/)

        REAL*8, dimension(8) ::
C                TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
C                  SO4    SEA    ANT    OCX    BCI    BCB    DST   VOL
     *  REFDRY=(/0.150, 1.000, 0.300, 0.200, 0.080, 0.080, 1.000,1.000/)

!nu  * ,REFWET=(/0.272, 1.808, 0.398, 0.318, 0.100, 0.100, 1.000,1.000/)
CKoch   DRYM2G=(/5.000, 2.866, 8.000, 8.000, 9.000, 9.000, 1.000,1.000/)
!nu     RHTMAG=(/1.788, 3.310, 1.756, 1.163, 1.000, 1.000, 1.000,1.000/)
!nu alt RHTMAG=(/1.982, 3.042, 1.708, 1.033, 1.000, 1.000, 1.000,1.000/)
!old *  WETM2G=(/8.345, 2.866, 7.811, 5.836, 9.000, 9.000, 1.000,1.000/)
!nu  * ,WETM2G=(/9.250, 2.634, 7.598, 5.180, 9.000, 9.000, 1.000,1.000/)
     * ,Q55DRY=(/2.191, 2.499, 3.069, 3.010, 1.560, 1.560, 1.000,1.000/)

     * ,DENAER=(/1.760, 2.165, 1.725, 1.500, 1.300, 1.300, 2.000,2.000/)

C     TROP AEROSOL 1850 BACKGROUND, INDUSTRIAL & BIO-BURNING PARAMETERS
      REAL*8, dimension(8) ::
C                TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
C                  SO4    SEA    ANT    OCX    BCI    BCB    DST   VOL
     *  FS8OPX=(/1.000, 1.000, 1.000, 1.000, 1.500, 1.500, 1.000, 1.00/)

     * ,FT8OPX=(/1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.300, 1.00/)

     * ,FRSULF=(/0.000, 0.000, 0.000, 0.330, 0.000, 0.000, 0.000, 1.00/)

     * ,PI0MAX=(/1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.00/)

!nu  * ,A8VEFF=(/ .200,  .200,  .200,  .200,  .200,  .200,  .200, .200/)

#ifdef ALTER_RADF_BY_LAT
!@var FS8OPX_orig saves initial FS8OPX values
!@var FT8OPX_orig saves initial FT8OPX values
      REAL*8, dimension(8) :: FS8OPX_orig, FT8OPX_orig
#endif

!      REAL*8, dimension(8) ::                                                    !ron
C                          MINERAL DUST PARAMETERS
C                         CLAY                  SILT
!     *REDUST=(/0.132D0,0.23D0,0.416D0,0.766D0,1.386D0,2.773D0,5.545D0,           !ron
!     &                                        8D0/) ! <- not used; 3 silt only   !ron
!nu  *  ,VEDUST=(/ 0.2, 0.2, 0.2, 0.2,   0.2, 0.2, 0.2, 0.2/)
!     *  ,RODUST=(/2.5D0,2.5D0,2.5D0,2.5D0,2.65D0,2.65D0,2.65D0,                  !ron
!     &                                        2.65D0/)! <- not used; 3 silt only !ron
!nu  *  ,FSDUST=(/ 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0/)
!nu  *  ,FTDUST=(/ 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0/)

!@var DUSTAB: specifies relative mixture of particles with Sinyuk 2003
!@+   and Patterson 1977 SW properties.
!@+   DUSTAB=1.0: all particles have Sinyuk 2003 properties
!@+   DUSTAB=0.0: all particles have Patterson 1977 properties
      real*8, parameter :: DUSTAB = 0.5

C-----------------------------------------------------------------------
C     GHG 1980 Reference Concentrations and Vertical Profile Definitions
C-----------------------------------------------------------------------

!@var KTREND if > 0 table GHG concentrations (Trend G) are used for
!@+             yr/day KYEARG/KJDAYG; if KTREND=0, GHG are set to PPMVK0
      INTEGER :: KTREND=1

!@var PPMV80  reference GHG concentrations (ppm)
      REAL*8, dimension(13) ::
C     GAS NUMBER    1         2    3      4    5         6           7
C                 H2O       CO2   O3     O2  NO2       N2O         CH4
#ifdef V2_O2_MODE /* temporary option to exactly match v2_branch */
     *   PPMV80=(/0d0, 337.90d0, 0d0,    21d4,0d0,  .3012d0,   1.5470d0
#else
     *   PPMV80=(/0d0, 337.90d0, 0d0,pO2*1.d6,0d0,  .3012d0,   1.5470d0
#endif
     *     ,.1666d-03,.3003d-03, 0d0,   .978D-04,  .0010D-10,  .0420d0/)
C              CCL3F1    CCL2F2   N2     CFC-Y       CFC-Z         SO2
C     GAS NUMBER    8         9   10        11          12          13

!@var PPMVK0  user set  GHG concentrations (ppm), used if KTREND=0
      REAL*8, dimension(12) ::
C     GAS  NUMBER   1         2    3      4    5         6           7
C                 H2O       CO2   O3     O2  NO2       N2O         CH4
     *   PPMVK0=(/0d0, 337.90d0, 0d0, 21.d4, 0d0,  .3012d0,   1.5470d0
     *               ,.1666d-03,  .3003d-03, 0d0, .978D-04, 0.0010D-10/)
C                        CCL3F1      CCL2F2   N2     CFC-Y       CFC-Z
C     GAS  NUMBER             8           9   10        11          12

C     Makiko GHG Trend Compilation  GHG.1850-2050.Dec1999 in GTREND
C     ---------------------------------------------------------------
!@var nghg nr. of well-mixed GHgases: CO2 N2O CH4 CFC-11 CFC-12 others
!@var nyrsghg max.number of years of prescr. greenhouse gas history
      INTEGER, PARAMETER :: nghg=6

!@var ghgyr1,ghgyr2 first and last year of GHG history
      INTEGER ghgyr1,ghgyr2
!@var ghgam,xref,xnow     GHG-mixing ratios in ppm,ppm,ppm,ppb,ppb,ppb
      REAL*8 XREF(nghg+1),XNOW(nghg+1)
      real*8, allocatable :: ghgam(:,:)

C     GTREND:  1980.,  337.9,  .3012,  1.547,  .1666,  .3003,  .0978,
C     ---------------------------------------------------------------

!@var KGGVDF,KPGRAD,KLATZ0 control parameters for vertical GHG profiles
!@+   -----------------------------------------------------------------
!@+   Minschwaner et al JGR (1998) CH4, N2O, CFC-12 Vertical profiles
!@+   IF(KGGVDF > 0) Then:
!@+      Gas decreases are linear with pressure, from unity at ground to
!@+      the fractional value PPMVDF(NGAS) at the top of the atmosphere.
!@+   Exponential decrease by EXP(-(Z-Z0)/H) is superimposed on this.
!@+   IF(KLATZ0 > 0) Then: Z0 depends on latitude, KGGVDF not used
!@+   KPGRAD>0: Pole-to-Pole lat. gradient (PPGRAD) is also superimposed
!@+   ------------------------------------------------------------------
!@var Z0,ZH   scale heights used for vertical profile (km)
!@var PPMVDF  frac. value at top of atmosphere (used if KGGVDF > 0)
!@var PPGRAD  Pole-to-Pole latitud.gradient for GHG (used if KPGRAD > 0)
      INTEGER :: KGGVDF=0, KPGRAD=1, KLATZ0=1

      REAL*8, dimension(12) ::
C     NUMBER   1    2    3    4  5    6    7    8     9   10   11  12
C             H2O  CO2  O3   O2 NO2  N2O  CH4 CFC11 CFC12 N2 CF-Y  CF-Z
     *   Z0=(/0.0, 0.0,0.0, 0.0,0.0, 16., 16., 16., 16., 0.0, 16., 16./)
     *  ,ZH=(/8.0, 8.0,8.0, 8.0,8.0, 30., 50., 30., 30., 0.0, 30., 30./)

C     GAS NUMBER    1     2    3    4    5         6         7
C                 H2O   CO2   O3   O2  NO2       N2O       CH4
     *  ,PPMVDF=(/1.0,  1.0, 1.0, 1.0, 1.0,  0.88888,  0.88888,
     *              0.88888,  0.88888, 1.0,  0.88888,  0.88888/)
C                    CCL3F1    CCL2F2   N2     CFC-Y     CFC-Z
C     GAS NUMBER          8         9   10        11        12

C     GAS  NUMBER   1     2    3    4    5         6         7
C                 H2O   CO2   O3   O2  NO2       N2O       CH4
     *  ,PPGRAD=(/0.0,  0.0, 0.0, 0.0, 0.0,   0.0100,   0.0900,
     *               0.0600,   0.0600, 0.0,   0.0600,   0.0600/)
C                    CCL3F1    CCL2F2   N2     CFC-Y     CFC-Z
C     GAS  NUMBER         8         9   10        11        12

C---------------------
C     Optional Tracers    used via setbak/getbak
C---------------------
      INTEGER, dimension(ITRMAX) :: ITR=1
      INTEGER :: NTRACE=0

      REAL*8, dimension(ITRMAX) ::
C                TRACER AEROSOL COMPOSITIONAL/TYPE PARAMETERS
     *  TRRDRY= .1d0
!nu  * ,TRVEFF= .2d0
     * ,TRADEN= 1.d0
!loc * ,FSTOPX= 1.d0
!loc * ,FTTOPX= 1.d0
     * ,FSTASC= 1.d0
     * ,FTTASC= 1.d0

      SAVE

      CONTAINS

      SUBROUTINE RCOMP1(NRFUN)
      use DOMAIN_DECOMP_ATM, only: AM_I_ROOT, grid
      use pario, only : par_open,par_close, variable_exists
     &                 ,get_dimlen,read_data
      use filemanager, only : file_exists

      IMPLICIT NONE
C     ------------------------------------------------------------------
C     Solar,GHG Trend, VolcAer Size Selection Parameters:    Defaults
C                                           Process       KYEARX  KJDAYX
c                                         SolarCon, UV       0       0
c                                         GH Gas Trend       0       0
c                                                         REFF0= 0.3
c                                                         VEFF0= 0.35
C     ------------------------------------------------------------------

c     NRFUN is now set as an argument from calling routine so that unit
c     numbers can be set automatically
      INTEGER :: NRFUN(14)
C          radfile1   2   3   4   5   6   7   8   9   A   B   C   D   E
!?    DATA NRFN0/71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84/

      INTEGER, SAVE :: IFIRST=1 ! ,NRFN0
      CHARACTER*80 EPSTAG,TITLE

      REAL*4 OZONLJ(44,46),R72X46(72,46)
      real*4, dimension(:,:), allocatable :: VTAUR4              !rjh
      real*4, allocatable :: vtau4(:,:,:),vreff4(:,:),hv4(:),lat4(:)

      INTEGER :: I,J,K,L,M,N,N1,N2,NRFU,KK,NN,IYEAR,IMONTH,JJDAYS,JYEARS
     *     ,JJDAYG,JYEARG,yr2S0i
      REAL*8 :: WAVNA,WAVNB,PFWI,TKOFPF,SUMV,EPK,EPL,DEP,SFNORM,D,O,Q,S
     *     ,OCM,WCM,YQSCCB
!@var GTAU,TGDATA temporary array to read data and pass it to RAD_UTILS
      REAL*8 :: GTAU(51,11,143),TGDATA(122,13)

      INTEGER :: N_BIN,fid
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SSI_IN
      REAL*8, ALLOCATABLE, DIMENSION(:) :: calyear,WS_IN,DS_IN,TSI_IN
      logical :: have_RADN9_file

!?    IF(LASTVC > 0) NRFUN=NRFN0
      IF(IFIRST < 1) GO TO 9999

C     ------------------------------------------------------------------
C     Input data are read as specified in the first CALL RCOMP1 (NRFUN).
C     Subsequent calls to RCOMP1 can be used to re-initialize parameters
C     in SETXXX subroutines to different values, but no new data is read
C     ------------------------------------------------------------------

C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C             e.g. MADVEL=123456   (zero digit skips input process)
C
C     MADO3M   =  1   Reads  Decadal Ozone files and Ozone trend file
C     MADAER   =  2   Reads  Aerosol 50y tropospheric climatology RFILE5
C     MADDST   =  3   Reads  Dust-windblown mineral climatology   RFILE6
C     MADVOL   =  4   Reads  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS   =  5   Reads  Epsilon cloud heterogeneity data     RFILE8
C     MADLUV   =  6   Reads  Lean formar Solar Spectral Irrad.    RFILE9
C
C                 Related Model Add-on Data Parameters set in RADPAR
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1   Reads  V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------


C              Initialize variables that might not otherwise get defined
C              ---------------------------------------------------------

       TAUWC(:) = 0  ;   TAUIC(:) = 0
      SIZEWC(:) = 0  ;  SIZEIC(:) = 0
      CLDEPS(:) = 0
      FPXCO2(:) = 1  ;  FPXOZO(:) = 1
       TLB(:) = 250  ;   TLT(:) = 250  ;   TLM(:) = 250
       SHL(:) =   0  ;   RHL(:) =   0
      SRAEXT(:,:) = 0  ;  SRASCT(:,:) = 0  ;  SRAGCB(:,:) = 0
      SRBEXT(:,:) = 0  ;  SRBSCT(:,:) = 0  ;  SRBGCB(:,:) = 0
      SRDEXT(:,:) = 0  ;  SRDSCT(:,:) = 0  ;  SRDGCB(:,:) = 0
      SRVEXT(:,:) = 0  ;  SRVSCT(:,:) = 0  ;  SRVGCB(:,:) = 0
      SRCEXT(:,:) = 0  ;  SRCSCT(:,:) = 0  ;  SRCGCB(:,:) = 0
      SRCPI0(:,:) = 0  ;  DBLPI0(:,:) = 0
      DBLEXT(:,:) = 0  ;  DBLSCT(:,:) = 0  ;  DBLGCB(:,:) = 0
      TRAALK(:,:) = 0  ;  TRBALK(:,:) = 0  ;  TRDALK(:,:) = 0
      TRVALK(:,:) = 0  ;  TRCALK(:,:) = 0  ;  TRGXLK(:,:) = 0
       U0GAS(:,:) = 0  ;   ULGAS(:,:) = 0
      TRACER(:,:) = 0
       EPLOW(:,:)  = 0 ;   EPMID(:,:) = 0  ;   EPHIG(:,:) = 0

      IF(LASTVC > 0) CALL SETATM
      IF(NL > LX)   call stop_model('rcomp1: increase LX',255)

C**** Use (global mean) pressures to get standard mid-latitude summer
C**** values for height, density, temperature, ozone, water vapor
      DO 120 L=1,NL+1
      PLB0(L)=PLB(L)
      CALL PHATMO(PLB0(L),HLB0(L),D,TLB(L),O,Q,S,OCM,WCM,1,2)
  120 CONTINUE
      DO 121 L=1,NL
      TLT(L)=TLB(L+1)
      TLM(L)=0.5D0*(TLB(L)+TLT(L))
  121 CONTINUE

!sl   De-activate surface layer computations
!sl   TAUSL(:)=0.0
!sl   FTAUSL(:)=0.0

C-----------------------------------------------------------------------
CR(1) Reads GTAU Asymmetry Parameter Conversion Table used within SGPGXG
C
C       (SGPGXG does Multiple Scattering Parameterization used in SOLAR)
C       ----------------------------------------------------------------

      NRFU=NRFUN(1)
      READ (NRFU) GTAU,TGDATA
      CALL SETGTS(TGDATA)
      CALL SET_SGPGXG(GTAU)


C-----------------------------------------------------------------------
CR(2)    Reads in Merged k-Distribution Tau Tables for Thermal Radiation
C        CFCs, H2O Continuum Tau Table, Merged k-Distr Planck Flux Table
C
C        (Reads: TAUCD0,TAUTBL,TAUWV0,TAUO30,PLANCK,XKCFC,H2OCN8,H2OCF8
c                DUCH4,SDUCH4,DUN2O,SDUN2O,ULOX,DUX      used in TAUGAS)
C       ----------------------------------------------------------------

      NRFU=NRFUN(2)
      READ(NRFU) title,TAUTBL
      READ(NRFU) title,TAUWV0
      READ(NRFU) title,TAUCD0
      READ(NRFU) title,TAUO30
      READ(NRFU) title,PLANCK
      READ(NRFU) title,XKCFC
      READ(NRFU) title,ULOX,DUX

      if(transmission_corrections) then
      NRFU=NRFUN(4)
      READ(NRFU) title,XTRUP,XTRDN,XTU0,XTD0
      READ(NRFU) title,XTFAC
      READ(NRFU) title,DXUP2,DXDN2    ! CO2
      READ(NRFU) title,DXUP3,DXDN3    ! O3
      READ(NRFU) title,DXUP6,DXDN6    ! N2O
      READ(NRFU) title,DXUP7,DXDN7    ! CH4
      READ(NRFU) title,DXUP8,DXDN8    ! CFC11
      READ(NRFU) title,DXUP9,DXDN9    ! CFC12
      READ(NRFU) title,DXUP13,DXDN13  ! SO2
      endif

C**** H2O Continuum Tau Tables (Ma_2000 or Ma_2004,Roberts,MT_CKD)
      NRFU=NRFUN(5)
      READ(NRFU) title,H2OCN8,XCSELF
      if(Am_I_Root()) write(6,*) title,' scaling factor:',XCSELF
      READ(NRFU) title,H2OCF8,XCFORN
      if(Am_I_Root()) write(6,*) title,' scaling factor:',XCFORN

C        Define Window Flux to Brightness Temperature Conversion Factors
C        ---------------------------------------------------------------

      do i=1,100 ; TKPFW(i    ) = TKofPF(85d1,9d2,    .001d0*I) ; end do
      do i=1, 90 ; TKPFW(i+100) = TKofPF(85d1,9d2,.1d0+.01d0*I) ; end do
      do i=1,440 ; TKPFW(i+190) = TKofPF(85d1,9d2,1.d0+ .1d0*I) ; end do
      do i=1,900 ; TKPFT(i    ) = TKofPF( 0d0,1d4,     dble(I)) ; end do

C                            PLANCK Table interpolation limit parameters
C                            -------------------------------------------
C-----------------------------------------------------------------------
CR(3)        Read Mie Scattering Parameters [Qext, Qscat, AsymParameter]
C            (1) Tropospheric Aerosols [11 Background, 8 Trop8 Aerosols]
C            (2) Clouds [5 Water, 5 non-spherical Ice, 5 Mie Ice Clouds]
C            (3) Desert Dust Aerosols  [25 particle sizes - to select 8]
C            (4) Volcanic Aerosols [20 particle sizes, 5 size variances]
C            (5) Sulfate  Aerosols [22 particle sizes, 0.1 - 10. micron]
C            (6) Soot   Aerosols [25 particle sizes, 0.001 - 5.0 micron]
C            -----------------------------------------------------------

      NRFU=NRFUN(3)

C                               GCM 11 background aerosol Mie parameters
C                               ----------------------------------------
      DO 301 N=1,11
      READ (NRFU,3000) TITLE
 3000 FORMAT(A80)
      READ (NRFU,3001) (SRAQEX(K,N),K=1,6)
 3001 FORMAT( 18X,6(F7.5,1X))
      READ (NRFU,3001) (SRAQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRAQCB(K,N),K=1,6)
  301 CONTINUE
      READ (NRFU,3002) (Q55A11(N),N=1,11)
 3002 FORMAT( 18X,6(F7.5,1X)/18X,6(F7.5,1X))
      READ (NRFU,3003) (REFA11(N),N=1,11)
 3003 FORMAT( 18X,6(F7.3,1X)/18X,5(F7.3,1X))
      READ (NRFU,3003) (VEFA11(N),N=1,11)
      DO 302 N=1,11
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRAQEX(K,N),K=1,33)
 3004 FORMAT( 14X,7(F7.5,1X),4(/14X,7(F7.5,1X)))
 3005 FORMAT(/14X,7(F7.5,1X),4(/14X,7(F7.5,1X)))
      READ (NRFU,3005) (TRAQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRAQCB(K,N),K=1,33)
  302 CONTINUE

C                       GCM 9 (of 10) climatology aerosol Mie parameters
C                       ------------------------------------------------
      DO 303 N=1,10
      IF(N==6) GO TO 303
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRBQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRBQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRBQCB(K,N),K=1,6)
  303 CONTINUE
      READ (NRFU,3002) (Q55B10(N),N=1,5),(Q55B10(N),N=7,10)
      READ (NRFU,3003) (REFB10(N),N=1,5),(REFB10(N),N=7,10)
      READ (NRFU,3003) (VEFB10(N),N=1,5),(VEFB10(N),N=7,10)
      DO 304 N=1,10
      IF(N==6) GO TO 304
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRBQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRBQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRBQCB(K,N),K=1,33)
  304 CONTINUE


C                               Cloud Water, Ice-non, Ice-Mie parameters
C                               ----------------------------------------
      DO 305 N=1,15
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRCQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRCQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRCQCB(K,N),K=1,6)
  305 CONTINUE
      READ (NRFU,3006) (Q55C15(N),N=1,15)
 3006 FORMAT( 18X,6(F7.5,1X)/18X,6(F7.5,1X)/18X,6(F7.5,1X))
      READ (NRFU,3007) (REFC15(N),N=1,15)
 3007 FORMAT( 18X,6(F7.3,1X)/18X,6(F7.3,1X)/18X,6(F7.3,1X))
      READ (NRFU,3007) (VEFC15(N),N=1,15)
      DO 306 N=1,15
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRCQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRCQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRCQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRCQAL(K,N),K=1,33)
  306 CONTINUE

C                               Desert Dust 25 sizes, Mie parameter data
C                               ----------------------------------------
      DO 307 N=1,25
      READ (NRFU,3001) (SRDQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRDQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRDQCB(K,N),K=1,6)
  307 CONTINUE
      READ (NRFU,3008) (Q55D25(N),N=1,25)
 3008 FORMAT( 18X,5(F7.5,1X),4(/18X,5(F7.5,1X)))
      READ (NRFU,3009) (REFD25(N),N=1,25)
 3009 FORMAT( 18X,12(F3.1,1X)/18X,12(F3.1,1X)/18X,F3.0)
      READ (NRFU,3010) (VEFD25(N),N=1,25)
 3010 FORMAT( 18X,12(F3.1,1X)/18X,12(F3.1,1X)/18X,F3.1)
      DO 308 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRDQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRDQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRDQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRDQAL(K,N),K=1,33)
  308 CONTINUE

      TRDQAB(:,:)=TRDQEX(:,:)-TRDQSC(:,:)  !  used in writer only

C                               Volcanic aerosol Mie size, variance data
C                               ----------------------------------------
      DO 313 M=1,5
      IF(M==4) GO TO 313
      DO 311 N=1,20
      READ (NRFU,3001) (SRVQEX(K,N,M),K=1,6)
      READ (NRFU,3001) (SRVQSC(K,N,M),K=1,6)
      READ (NRFU,3001) (SRVQCB(K,N,M),K=1,6)
  311 CONTINUE
      READ (NRFU,3011) (Q55V20(N,M),N=1,20)
 3011 FORMAT( 18X,5(F7.5,1X),3(/18X,5(F7.5,1X)))
      READ (NRFU,3012) (REFV20(N,M),N=1,20)
 3012 FORMAT( 18X,12(F3.1,1X)/18X,8(F3.1,1X))
      READ (NRFU,3012) (VEFV20(N,M),N=1,20)
      DO 312 N=1,20
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRVQEX(K,N,M),K=1,33)
      READ (NRFU,3005) (TRVQSC(K,N,M),K=1,33)
      READ (NRFU,3005) (TRVQCB(K,N,M),K=1,33)
      READ (NRFU,3005) (TRVQAL(K,N,M),K=1,33)
  312 CONTINUE
  313 CONTINUE
      DO 316 N=1,20
      DO 314 K=1,6
      SRVQEX(K,N,4)=(SRVQEX(K,N,3)+SRVQEX(K,N,5))/2.D0
      SRVQSC(K,N,4)=(SRVQSC(K,N,3)+SRVQSC(K,N,5))/2.D0
      SRVQCB(K,N,4)=(SRVQCB(K,N,3)+SRVQCB(K,N,5))/2.D0
  314 CONTINUE
      Q55V20(N,4)=(Q55V20(N,3)+Q55V20(N,5))/2.D0
      REFV20(N,4)=(REFV20(N,3)+REFV20(N,5))/2.D0
      VEFV20(N,4)=(VEFV20(N,3)+VEFV20(N,5))/2.D0
      DO 315 K=1,33
      TRVQEX(K,N,4)=(TRVQEX(K,N,3)+TRVQEX(K,N,5))/2.D0
      TRVQSC(K,N,4)=(TRVQSC(K,N,3)+TRVQSC(K,N,5))/2.D0
      TRVQCB(K,N,4)=(TRVQCB(K,N,3)+TRVQCB(K,N,5))/2.D0
      TRVQAL(K,N,4)=(TRVQAL(K,N,3)+TRVQAL(K,N,5))/2.D0
  315 CONTINUE
  316 CONTINUE

C                            Sulfate aerosol, Mie parameter 22-size data
C                            -------------------------------------------
      DO 321 N=1,22
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRUQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRUQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRUQCB(K,N),K=1,6)
  321 CONTINUE
      READ (NRFU,3008) (Q55U22(N),N=1,22)
      READ (NRFU,3013) (REFU22(N),N=1,22)
 3013 FORMAT( 18X,5(F7.3,1X),4(/18X,5(F7.3,1X)))
      READ (NRFU,3013) (VEFU22(N),N=1,22)
      DO 322 N=1,22
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRUQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRUQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRUQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRUQAL(K,N),K=1,33)
  322 CONTINUE

C                               Soot aerosol, Mie parameter 25-size data
C                               ----------------------------------------
      DO 323 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRSQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRSQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRSQCB(K,N),K=1,6)
  323 CONTINUE
      READ (NRFU,3008) (Q55S25(N),N=1,25)
      READ (NRFU,3013) (REFS25(N),N=1,25)
      READ (NRFU,3013) (VEFS25(N),N=1,25)
      DO 324 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRSQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRSQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRSQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRSQAL(K,N),K=1,33)
  324 CONTINUE

C                            Seasalt aerosol, Mie parameter 22-size data
C                            Nitrate aerosol, Mie parameter 22-size data
C                            (Water) aerosol, Mie parameter 22-size data
C                            Organic aerosol, Mie parameter 22-size data
C                            -------------------------------------------
      N1=23
      DO 326 KK=1,4
      N2=N1+21
      DO 325 N=N1,N2
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRUQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRUQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRUQCB(K,N),K=1,6)
  325 CONTINUE
      READ (NRFU,3008) (Q55U22(N),N=N1,N2)
      READ (NRFU,3013) (REFU22(N),N=N1,N2)
      READ (NRFU,3013) (VEFU22(N),N=N1,N2)
      N1=N2+1
  326 CONTINUE
      N1=23
      DO 328 KK=1,4
      N2=N1+21
      DO 327 N=N1,N2
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRUQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRUQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRUQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRUQAL(K,N),K=1,33)
  327 CONTINUE
      N1=N2+1
  328 CONTINUE

C                        Sinyuk Desert Dust 25 sizes, Mie parameter data
C                        -----------------------------------------------

      DO 347 N=1,25
      READ (NRFU,3001) (YRDQEX(K,N),K=1,6)
      READ (NRFU,3001) (YRDQSC(K,N),K=1,6)
      READ (NRFU,3001) (YRDQCB(K,N),K=1,6)
  347 CONTINUE
      READ (NRFU,3008) (Y55D25(N),N=1,25)
      READ (NRFU,3009) (REFD25(N),N=1,25)
      READ (NRFU,3010) (VEFD25(N),N=1,25)
      DO 348 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRDQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRDQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRDQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRDQAL(K,N),K=1,33)
  348 CONTINUE

      TRDQAB(:,:)=TRDQEX(:,:)-TRDQSC(:,:)  !  used in writer only

!-----------------------------------------------------------------------
!     Create external mixture of Patterson and Sinyuk dust particles
      DO 352 N=1,25
      DO 351 K=1,6
      YQSCCB= YRDQSC(K,N)*YRDQCB(K,N)*DUSTAB
     *      + SRDQSC(K,N)*SRDQCB(K,N)*(1.d0-DUSTAB)
      SRDQEX(K,N)=DUSTAB*YRDQEX(K,N)+(1.d0-DUSTAB)*SRDQEX(K,N)
      SRDQSC(K,N)=DUSTAB*YRDQSC(K,N)+(1.d0-DUSTAB)*SRDQSC(K,N)
      SRDQCB(K,N)=YQSCCB/(1.d-10+SRDQSC(K,N))
 351  CONTINUE
      Q55D25(N)=DUSTAB*Y55D25(N)+(1.d0-DUSTAB)*Q55D25(N)
 352  CONTINUE

C-----------------------------------------------------------------------
CR(7)        Read Stratospheric Volcanic binary data
C            (NVOLMON months (years JVOLYI to JVOLYE) x NVOLLAT latitudes)
C            If KyearV<0 use the NVOLMON-month mean as background aerosol
C            ---------------------------------------------------------
      NRFU=NRFUN(7)
      if(madvol==1) then
        READ (NRFU) TITLE,NVOLMON,JVOLYI,JVOLYE
        NVolLat = 24 ; NVolK = 4
        IF(TITLE(1:9).ne.'OD Header')
     &      call stop_model('rcomp1: use new RADN7 header file',255)
        ALLOCATE (VTauTJK(NVOLMON,NVolLat,NVolK), HVolKM(NVolK+1))
        ALLOCATE (VReffTJ(NVOLMON,NVolLat), VTAUR4(NVOLMON,NVolLat))
        DO K=1,NVolK
          READ (NRFU) TITLE,VTAUR4
          DO J=1,NVolLat
            SUMV=0.
            DO I=1,NVOLMON
              VTauTJK(I,J,K)=VTAUR4(I,J)
              SUMV=SUMV+VTAUR4(I,J)
            END DO
            if(kyearv < 0) VTauTJK(1,J,K)=SUMV/NVOLMON
          END DO
        END DO
        READ (NRFU) TITLE,VTAUR4
        DO J=1,NVolLat
          SUMV=0.
          DO I=1,NVOLMON
            VReffTJ(I,J)=VTAUR4(I,J)
            SUMV=SUMV+VTAUR4(I,J)
          END DO
          if(kyearv < 0) VReffTJ(1,J)=SUMV/NVOLMON
        END DO
        deallocate (VTAUR4)
      else if(madvol==2) then
        READ (NRFU) TITLE ; rewind NRFU
        IF(TITLE(1:12).ne.'CMIP6 Header')
     &      call stop_model('rcomp1: use CMIP6 RADN7 header file',255)
        READ (NRFU) TITLE,NVOLMON,JVOLYI,JVOLYE,NVolLat,NVolK
        ALLOCATE (VTauTJK(NVOLMON,NVolLat,NVolK))
        ALLOCATE (VTau4(NVOLMON,NVolLat,NVolK))
        ALLOCATE (VReff4(NVOLMON,NVolLat))
        ALLOCATE (VReffTJ(NVOLMON,NVolLat))
        ALLOCATE (HVOLKM(NVolK+1),ELATVOL(NVolLat+1))
        ALLOCATE (hv4(NVolK+1),LAT4(NVolLat+1))
        READ (NRFU) TITLE,VTau4  ; VTauTJK = VTau4
        READ (NRFU) TITLE,VReff4 ; VReffTJ = VReff4
        READ (NRFU) TITLE,hv4    ; HVOLKM  = hv4
        READ (NRFU) TITLE,LAT4   ; ELATVOL = lat4
        deallocate (VTau4,VReff4,hv4,LAT4)
        if(kyearv < 0) then
          do j=1,NVolLat
          VReffTJ(1,J)=sum(VReffTJ(:,j))/nvolmon
            do k=1,NVolK
            VTauTJK(1,J,K)=sum(VTauTJK(:,j,k))/nvolmon
            end do
          end do
        end if
      end if
C-----------------------------------------------------------------------
CR(8)  ISCCP Derived Cloud Variance (EPSILON) Cloud Optical Depth Factor
C      Low, Mid, High  Cloud Optical Depths are Reduced by (1 - EPSILON)
C
C               INPUT DATA FILE:  UNIT = INFILE
C                                 TAG  = EPSTAG  (CHARACTER*80)
C                                 DATA = EPLMHC  (72,46,12,4) REAL*4
C
C      Data are 72X46 Monthly Mean Low, Mid, High, Column EPSILON Values
C      Cloud Heterogeneity selections used in UPDEPS, GETEPS (in SETCLD)
C
C             EPSCON  Column Cloud Inhomogeneity EPSILON (when KCLDEP=1)
C             KCLDEP  Selects Cloud Inhomogeneity Option (0-4):
C                     KCLDEP =  0  Sets Column CLDEPS to Zero
C                     KCLDEP =  1  Sets Column CLDEPS to EPSCON
C                     KCLDEP =  2  Keeps whatever is specified in CLDEPS
C                     KCLDEP =  3  Uses: Column EPCOL(72,46) Climatology
C                     KCLDEP =  4  Uses: Ht Dep EPLOW, EPMID, EPHIG Data
C               --------------------------------------------------------

      IF(MADEPS < 1) GO TO 899
      NRFU=NRFUN(8)
      READ (NRFU) EPSTAG,EPLMHC


      DO 810 N=1,4
      DO 810 M=1,12
      DO 810 I=1,72
!**** extend northern-most non-neg.value to N.Pole
      J=46 ! MLAT46
      do while (EPLMHC(I,J,M,N) < 0) ; J=J-1 ; end do
      IF(J < 46) EPLMHC(I,J+1:46,M,N) = EPLMHC(I,J,M,N)
!**** extend southern-most non-neg.value to S.Pole
      J=1
      do while (EPLMHC(I,J,M,N) < 0) ; J=J+1 ; end do
      IF (J > 1) EPLMHC(I,1:J-1,M,N) = EPLMHC(I,J,M,N)
      IF (J==46) GO TO 810
!**** linearly interpolate across remaining intervals with EP<0
  805 J=J+1                                    ! find start of interval:
      do while (EPLMHC(I,J,M,N) >= 0)
        J=J+1 ; IF (J > 46) GO TO 810
      end do
      K=J-1
      EPK=EPLMHC(I,K,M,N)
      J=J+1                                    ! find  end  of interval:
      do while (EPLMHC(I,J,M,N) < 0) ; J=J+1 ; end do
      L=J
      EPL=EPLMHC(I,L,M,N)            ! EPk>=0, EPk+1,...,EPl-1<0, EPl>=0
      DEP=(EPL-EPK)/(L-K)
      do NN=1,L-1-K                  ! replace EP(k+1)...EP(l-1)
        EPLMHC(I,K+NN,M,N)=EPK+NN*DEP
      end do
      IF (J < 46) GO TO 805
  810 CONTINUE

  899 CONTINUE


C-----------------------------------------------------------------------
CR(E)
C             KCLDEM  Selects: Top-Cloud (Thermal) Scattering Correction
C                     KCLDEM =  0  Utilizes Non-scattering approximation
C                     KCLDEM =  1  Modifies emission and transmission by
C                                  top cloud (over-rides old correction)
C             ----------------------------------------------------------

      NRFU=NRFUN(14)
      READ (NRFU) RIJTPG,FDXTPG,FEMTPG


C-----------------------------------------------------------------------
CR(9)         Read Judith Lean Solar UV and Solar Constant Variability
C                                                Monthly-Mean Solar UV
C                                      ---------------------------------
      iMS0X = MS0X

      IF(KSOLAR < 0) GO TO 949
      IF(MADLUV < 1) THEN
        WS_SSI(:)=WS_SSI(:)/1000.D0
        DS_SSI(:)=DS_SSI(:)/1000.D0
        W1_SSI(:)=WS_SSI(:)-0.5D0*DS_SSI(:)
        GO TO 949
      END IF
!      NRFU=NRFUN(9)

!      IF(KSOLAR.ne.9) THEN
!        READ(NRFU,'(a80)') TITLE
!        if(ksolar >= 2 .and. TITLE(1:3).ne.'ANN')
!     &    call stop_model('rcomp1: change RADN9 to ann.file',255)
!        if(ksolar < 2 .and. TITLE(1:3)=='ANN')
!     &    call stop_model('rcomp1: change RADN9 to monthly file',255)
!        READ(NRFU,'(5F14.2)') WSLEAN   !  1:190
!        READ(NRFU,'(a80)') TITLE
!        READ(NRFU,'(5E14.3)') DSLEAN   !  1:190

      have_RADN9_file = file_exists('RADN9')

      if(have_RADN9_file) then
        fid=par_open(grid,'RADN9','read')
        iMs0X=get_dimlen(grid,fid,'time')
        N_BIN=get_dimlen(grid,fid,'wlen')
        ALLOCATE (TSI_IN(iMS0X),calyear(iMS0X))
        ALLOCATE (WS_IN(N_BIN),DS_IN(N_BIN),SSI_IN(N_BIN,iMS0X))
        if(variable_exists(grid,fid,'calyear'))then
          call read_data(grid,fid,'calyear',calyear,bcast_all=.true.)
        else
          call stop_model('missing calyear in RADN9 file',255)
        endif
        if(variable_exists(grid,fid,'wlen'))then
          call read_data(grid,fid,'wlen',WS_IN,bcast_all=.true.)
        else
          call stop_model('missing the wlen variable in RADN9 file',255)
        endif
        if(variable_exists(grid,fid,'wlenbinsize'))then
          call read_data(grid,fid,'wlenbinsize',DS_IN,bcast_all=.true.)
        else
          call stop_model('missing wlenbinsize in RADN9 file',255)
        endif
          if(variable_exists(grid,fid,'ssi'))then
          call read_data(grid,fid,'ssi',SSI_IN,bcast_all=.true.)
        else
          call stop_model('missing the ssi variable in RADN9 file',255)
        endif
        if(variable_exists(grid,fid,'tsi'))then
          call read_data(grid,fid,'tsi',TSI_IN,bcast_all=.true.)
        else
          call stop_model('missing the tsi variable in RADN9 file',255)
        endif
        call par_close(grid,fid)
      else
        call stop_model('missing the RADN9 file',255)
      endif

      WS_SSI(:)=WS_IN(N_BIN-189:N_BIN)/1000.D0
      DS_SSI(:)=DS_IN(N_BIN-189:N_BIN)/1000.D0
      W1_SSI(:)=WS_SSI(:)-0.5D0*DS_SSI(:)

!        WSLEAN(:)=WSLEAN(:)/1000.D0
!        DSLEAN(:)=DSLEAN(:)/1000.D0
!        W1LEAN(:)=WSLEAN(:)-0.5D0*DSLEAN(:)

!        READ(NRFU,'(a80)') TITLE
!        READ(NRFU,'(a80)') TITLE
!        READ(NRFU,'(a80)') TITLE
!        if(TITLE(1:5).ne.'MS0X=') then  ! old no_header file
!          backspace (NRFU)
!        else
!          read (title(6:80),*) iMs0X
!        endif
!      END IF
      ALLOCATE (UV_SSI(iMS0X,190),TSI1(iMS0X),TSI2(iMS0X))
      UV_SSI(:,:)=TRANSPOSE(SSI_IN(N_BIN-189:N_BIN,:))
      TSI1(:)=TSI_IN(:)
      TSI2(:)=TSI_IN(:)
      yr1S0=calyear(1)
      yr2S0=calyear(iMS0X)
      DEALLOCATE(WS_IN,DS_IN,SSI_IN,TSI_IN,calyear)
!      IF(KSOLAR < 2) THEN
C****   Read in monthly-mean data
!        DO I=1,iMs0X
!          READ(NRFU,'(2I6,3F17.6)') IYEAR,IMONTH,TSI1(I),TSI2(I)
!          READ(NRFU,'(5E14.6)')     FSLEAN    ! 1:190
!          SFNORM = TSI1(I) / SUM(FSLEAN(:)*DSLEAN(:))
!          UVLEAN(I,:)=FSLEAN(:)*SFNORM
!        END DO
!      ELSE
C****   Read in annual-mean data
!        DO I=1,iMs0X
!          IF(KSOLAR.ne.9) THEN
!            READ(NRFU,'(F12.1,2F15.4)',end=908) yr2S0,TSI1(I),TSI2(I)
!          ELSE
!            READ(NRFU,'(I6,2F17.6)',end=908) yr2S0i,TSI1(I),TSI2(I)
!            yr2S0=real(yr2S0i)+0.5
!          END IF
!          if(I==1) yr1S0 = yr2S0
!          IF(KSOLAR.ne.9) THEN
!            READ(NRFU,'(5E14.6)')   FSLEAN    ! 1:190
!            SFNORM=TSI1(I) / SUM(FSLEAN(:)*DSLEAN(:))
!            UVLEAN(I,:)=FSLEAN(:)*SFNORM
!          ELSE ! ksolar=9
!            READ(NRFU,'(5E14.6)')               (UVLEAN(I,K),K=1,190)
!          ENDIF
!        END DO
  908   if(Am_I_Root()) write(6,*) 'read S0-history: ',yr1S0,' - ',yr2S0
!      END IF

  949 CONTINUE

C-----------------------------------------------------------------------
CR(C)     Read:    Elaine Mathews 10 Fractional Vegetation Distributions
C         10 global maps (72x46) depict fractional vegetation/soil types
C         Map-1 (bright sand) + Map-10 (black dirt) define desert albedo
C         (sum of Maps 1-10 over land-area (ILON,JLAT) grid boxes = 1.0)
C
C         Map-11 refers to plankton concentrations over ocean areas that
C         are yet to be implemented.
C         --------------------------------------------------------------





C-----------------------------------------------------------------------
CR(D)      Read:   1   FOCEAN   72x46 ocean fraction   (FOCEAN = 0 or 1)
C                  2   FLAKE    72x46 lake  fraction
C                  3   FGRND    72x46 lake  fraction
C                  4   FGICE    72x46 glacial ice fraction
C                               (FLAKE + FGRND + FGICE + FOCEAN = 1.000)
C
C                  5   ZATMO    72x46 topography (ocean = 0.0)
C                  6   HOCEAN   72x46 ocean depth
C                  7   HLAKE    72x46 lake  depth
C                  8   HGICE    72x46 glice depth
C                  9   ZSOLID   72x46 topography of solid ground surface
C                  -----------------------------------------------------
C
C     FOLGIZ is for off-line use only, and is not used in GCM radiation.
C     GCM supplies dynamically changing POCEAN,POICE,PEARTH,PLICE values
C     ------------------------------------------------------------------



  999 CONTINUE

      IFIRST=0
 9999 CONTINUE


C     ---------------------------------------------------------------
C     LASTVC     Initialize:  Default Atmospheric Layering, Structure
C                (for Off-Line use)  as Specified by LASTVC Parameter
C                If LASTVC < 0, GCM defines all Radiation Model Input
C     otherwise:
C                Each LASTVC digit(6) specifies a model configuration
C          e.g.:    LASTVC= 123456
C                L=0,1,..9  Layers NL=  Any,GCM12,GCM23,Pset,Hset,etc
C                A=0,1,..6  Atmosphere  Any,Trop,MLS,MLW,SAS,SAW,Std
C                S=0,1,..9  Surf Types  POCEAN=1,PEARTH=1,POICE=1,etc
C                T=0,1,..9  Tracer Aer  Tau=0,  Tau=0.1 Aer Comp(1-9)
C                V=0,1,..9  Vegetation  Sand,Tundra,Grass,Shrubs, etc
C                C=0,1,..9  Cloud,R=10  Clim Cloud Tau in Layer(1,-9)
C                ----------------------------------------------------

      IF(LASTVC >= 0) CALL SETATM

C             -------------------------------------------------------
C             Set Solar Constant for Default Reference Time: Jan 1950
C             Default used for KSOLAR(=1) is that specified in RADPAR
C             -------------------------------------------------------

      JJDAYS=1
      JYEARS=1950
      IF(KJDAYS > 0)             JJDAYS=KJDAYS
      IF(KYEARS > 0)             JYEARS=KYEARS
C----------------------------------------------
                      CALL SETSOL(JYEARS,JJDAYS)
C----------------------------------------------


C             -------------------------------------------------------
C             Set Default Greenhouse Gas Reference Year to:  Mid 1980
C             Default used for KTREND(=1) is that specified in RADPAR
C             -------------------------------------------------------

      JJDAYG=184
      JYEARG=1980
C----------------------------------------------
                      CALL SETGHG(JYEARG,JJDAYG)
C----------------------------------------------
      IF(KJDAYG > 0)             JJDAYG=KJDAYG
      IF(KYEARG > 0)             JYEARG=KYEARG
C----------------------------------------------
                      CALL UPDGHG(JYEARG,JJDAYG)
C----------------------------------------------

C--------------------------------
                      CALL SETGAS
C
                                   CALL SETBAK
      IF(MADAER > 0.or.NTRACE > 0) CALL SETAER
      ! SETDST ops deferred to first call to GETDST once dust info known
      !IF(MADDST > 0) CALL SETDST
C--------------------------------


C               -----------------------------------------------------
C               Set Volcanic Aerosol Effective Variance Default Value
C                   Particle Size(REFF0=0.3) when not known from data
C                   (VEFF0=0.35 is value based on thermal ISAMS data)
C                   -------------------------------------------------

C----------------------------------------------
      IF(MADVOL > 0) CALL SETVOL
C----------------------------------------------

C--------------------------------
                      CALL SETCLD

C--------------------------------

      CALL SOLAR0

      RETURN
      END SUBROUTINE RCOMP1

      SUBROUTINE RCOMPT
      use SURF_ALBEDO, only : UPDSUR
      use AerParam_mod, only : updateAerosol,updateAerosol2
      use DustParam_mod, only : upddst2
      use O3mod, only : updO3d,updO3d_solar,plbo3,nlo3
#ifdef HIGH_FREQUENCY_O3_INPUT
      use O3mod, only : UPDO3D_highFrequency
#endif
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C     Time Trend Selection Parameters and Options:
C     -------------------------------------------
C
C             The Nominal Default Values are KYEARX = 0, and KJDAYX = 0,
C             in which case RADPAR supplied Time JYEAR and JDAY are used
C
C             When Non-Zero Values are specified for  KYEARX and KJDAYX,
C             the JYEAR,JDAY Time Dependence of the Specified Process is
C             over-ridden by the Non-Zero KYEARX and KJDAYX Value.
C             ----------------------------------------------------------
C                                  Process       KYEARX  KJDAYX
c             KYEARS,KJDAYS        SolarCon, UV       0       0
c             KYEARG,KJDAYG        GH Gas Trend       0       0
c             KYEARO,KJDAYO        Ozone Distr        0       0
c             KYEARA,KJDAYA        AerClimtolgy       0       0
c             KYEARD,KJDAYD        Desert Dust        0       0
c             KYEARV,KJDAYV        Volcanic Aer       0       0
c             KYEARE,KJDAYE        Epsilon Clds       0       0
c             KYEARR,KJDAYR        Refl Surface       0       0

C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C             e.g. MADVEL=123456   (zero digit skips input process)
C
C     MADAER  =  2  Updates  Aerosol 50y tropospheric climatology RFILE5
C     MADDST  =  3  Updates  Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  4  Updates  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  5  Updates  Epsilon cloud heterogeneity data     RFILE8
C     MADLUV  =  6  Updates  Lean format Spectral Solar Irrad.    RFILE9
C
C                 Related Model Add-on Data Parameters set in RADPAR
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1          V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------
      INTEGER JJDAYS,JYEARS,JJDAYG,JYEARG,JJDAYO,JYEARO,JJDAYA,JYEARA
     *     ,JJDAYD,JYEARD,JJDAYV,JYEARV,JJDAYE,JYEARE,JJDAYR,JYEARR

C                      -------------------------------------------------
C                      Set Seasonal and Time (JDAY) Dependent Quantities
C                      -------------------------------------------------

      JJDAYS=JDAY
      JYEARS=JYEAR
      IF(KJDAYS > 0)             JJDAYS=KJDAYS
      IF(KYEARS > 0)             JYEARS=KYEARS
C----------------------------------------------
      IF(MADLUV > 0) CALL UPDSOL(JYEARS,JJDAYS)
C----------------------------------------------

      JJDAYG=JDAY
      JYEARG=JYEAR
      IF(KJDAYG > 0)             JJDAYG=KJDAYG
      IF(KYEARG > 0)             JYEARG=KYEARG
C----------------------------------------------
      IF(MADGHG > 0) CALL UPDGHG(JYEARG,JJDAYG)
C----------------------------------------------

      JJDAYO=JDAY
      JYEARO=JYEAR
      IF(KJDAYO.ne.0)             JJDAYO=KJDAYO
      IF(KYEARO.ne.0)             JYEARO=KYEARO
C----------------------------------------------
      CALL UPDO3D(JYEARO,JJDAYO,O3JDAY,O3JREF)
#ifdef HIGH_FREQUENCY_O3_INPUT
      CALL UPDO3D_highFrequency(JYEARO,JJDAYO,O3JDAY_HF_modelLevels)
#endif
      CALL UPDO3D_solar(JJDAYO,S00WM2*RATLS0,O3JDAY)
C----------------------------------------------

      JJDAYA=JDAY
      JYEARA=JYEAR
      IF(KJDAYA > 0)             JJDAYA=KJDAYA
      IF(KYEARA.ne.0)            JYEARA=KYEARA
C----------------------------------------------
      IF(MADAER.eq.3) THEN
        CALL updateAerosol2(JYEARA,JJDAYA, a6jday, plbaer)
      ELSE IF(MADAER.ne.0) THEN
        CALL updateAerosol(JYEARA,JJDAYA, a6jday, plbaer)
      ENDIF
C----------------------------------------------

      JJDAYD=JDAY
      JYEARD=JYEAR
      IF(KJDAYD > 0)             JJDAYD=KJDAYD
      IF(KYEARD.ne.0)            JYEARD=KYEARD
C----------------------------------------------
      IF(MADDST > 0) CALL UPDDST2(JYEARD,JJDAYD)
C----------------------------------------------

      JJDAYV=JDAY
      JYEARV=JYEAR
      IF(KJDAYV > 0)             JJDAYV=KJDAYV
      IF(KYEARV.ne.0)             JYEARV=KYEARV
C----------------------------------------------
      IF(MADVOL > 0) CALL UPDVOL(JYEARV,JJDAYV)
C----------------------------------------------

      JJDAYE=JDAY
      JYEARE=JYEAR
      IF(KJDAYE > 0)             JJDAYE=KJDAYE
      IF(KYEARE > 0)             JYEARE=KYEARE
C----------------------------------------------
      IF(MADEPS > 0) CALL UPDEPS(JYEARE,JJDAYE)
C----------------------------------------------

      JJDAYR=JDAY
      JYEARR=JYEAR
      IF(KJDAYR > 0)             JJDAYR=KJDAYR
      IF(KYEARR > 0)             JYEARR=KYEARR
C----------------------------------------------
                      CALL UPDSUR(JYEARR,JJDAYR)
C----------------------------------------------

      RETURN
      END SUBROUTINE RCOMPT

      SUBROUTINE RCOMPX
      use SURF_ALBEDO, only : getsur
      use O3mod, only : plbo3,nlo3,plbo3_traditional,nlo3_traditional
#ifdef SCM
      use SCM_COM, only : SCMopt,SCMin
#endif
      IMPLICIT NONE
      integer k
C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C             e.g. MADVEL=123456   (zero digit skips process)
C
C     MADO3M  =  1           Makiko   1951-1997 Ozone climatology RFILEA
C     MADAER  =  2  Updates  Aerosol 50y tropospheric climatology RFILE5
C     MADDST  =  3  Updates  Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  4  Updates  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  5           Epsilon cloud heterogeneity data     RFILE8
C     MADLUV  =  6           Lean format Spectral Solar Irrad.    RFILE9
C
C                 Related Model Add-on Data Parameters set in RADPAR
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1          V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------
C
C      -----------------------------------------------------------------
C      Get Surface, Atmosphere, Sun Angle, Radiative Forcing, etc. Input
C      to compute Solar/Thermal Radiation for given (JLAT,ILON) Grid-box
C
C      The Radiation Model utilizes Data with 72x46 (lon,lat) resolution
C                 for GCM resolution other than 72x46, set JLAT and ILON
C                 to appropriately Sample  (rather than interpolate) the
C                 72x46 aerosol, ozone, cloud heterogeneity data sets
C
C      The Radiation Model can accommodate arbitrary vertical resolution
C      -----------------------------------------------------------------


C--------------------------------
      if(set_gases_internally) then
!!!                   CALL GETO3D(ILON,JLAT) ! may have to be changed ??
      if(use_o3_ref > 0 )then
        CALL REPART (O3JREF(1,IGCM,JGCM),
     *          PLBO3_traditional,NLO3_traditional+1, ! in
     *                        U0GAS(1,3),PLB0, NL+1)  ! out, ok if L1>1 ?
        ! next block may seem weird but it is here to allow RCOMPX calls with
        ! reference ozone in part of the atmosphere and tracer below:
        if(use_tracer_chem(1) > 0) then
          U0GAS(1:use_tracer_chem(1),3)=chem_IN(1,1:use_tracer_chem(1))
        endif
        FULGAS(3)=1.d0
      else
        CALL REPART (O3JDAY(1,IGCM,JGCM),PLBO3,NLO3+1, ! in
     *                        U0GAS(1,3),PLB0, NL+1)   ! out, ok if L1>1 ?
#ifdef HIGH_FREQUENCY_O3_INPUT
        ! Overwrite the lm_gcm levels with higher frequency ozone, leaving
        ! climatology above those levels:
        U0GAS(1:lm_gcm,3)=O3JDAY_HF_modelLevels(1:lm_gcm,IGCM,JGCM)
        FULGAS(3)=1.d0
#endif
#ifdef SCM
        if(SCMopt%ozone)then
        ! Overwrite specified SCM levels (indicated by non-zero values),
        ! leaving climatology above those levels:
          do k = 1,lm_gcm
            if(SCMin%O3(k) > 0.) U0GAS(k,3)=SCMin%O3(k)
          enddo
          FULGAS(3)=1.d0
        endif
#endif
        ! considering this move to here from setgas:
        ! chem_out(:,1)=U0GAS(:,3)*FULGAS(3) ! save climatology O3 for chem
        ! and might then need something like:
        ! IF(KPFOZO==1)chem_out(1:NL0,1)=chem_out(1:NL0,1)*FPXOZO(1:NL0)
        if(use_tracer_chem(1) > 0) then
          U0GAS(1:use_tracer_chem(1),3)=chem_IN(1,1:use_tracer_chem(1))
          FULGAS(3)=1.d0
        endif
      endif
                      CALL GETGAS
      else
        CALL TAUGAS
      endif
C--------------------------------


C--------------------------------
      if(set_aerosols_internally) then
      SRBEXT=1.d-20 ; SRBSCT=0. ; SRBGCB=0. ; TRBALK=0.
      IF(MADBAK > 0) CALL GETBAK

      IF(MADAER.ne.0.OR.NTRACE > 0) THEN ; CALL GETAER
       ELSE ; SRAEXT=0.     ; SRASCT=0. ; SRAGCB=0. ; TRAALK=0. ; END IF
      IF(MADDST > 0) THEN ; CALL GETDST
       ELSE ; SRDEXT=0.     ; SRDSCT=0. ; SRDGCB=0. ; TRDALK=0. ; END IF
      IF(MADVOL > 0) THEN ; CALL GETVOL
       ELSE ; SRVEXT=0.     ; SRVSCT=0. ; SRVGCB=0. ; TRVALK=0. ; END IF
      chem_out(:,2)=SRVEXT(:,6) ! save 3D aerosol extinction in SUB RADIA
      endif
C--------------------------------


C--------------------------------  (GETSUR sets albedo needed by GETCLD)
                      CALL GETSUR(
     i     snoage_fac_max,
     i     MLAT46,jnorth,KEEPAL,KSIALB,KZSNOW,MADSUR,
     i     COSZ,PLANCK,PLANCK_TMIN,PLANCK_TMAX,
     i     ILON,JLAT,
     i     AGESN,POCEAN,POICE,PEARTH,PLICE,PLAKE,zlake,
     i     TGO,TGOI,TGE,TGLI,ZOICE,FMP,ZSNWOI,zmp,
     i     SNOWOI,SNOWD,SNOWLI,SNOW_FRAC,WEARTH,WMAG,PVT,dalbsn,
     i     flags,LOC_CHL,
     o     BXA,PRNB,PRNX,SRBALB,SRXALB,TRGALB,
     o     BGFEMD,BGFEMT,
     o     DTRUFG,FTRUFG
     &     )
                      CALL GETEPS
                      CALL GETCLD
C--------------------------------

C--------------------------------
                      CALL THERML

                      CALL SOLARM
C--------------------------------
      RETURN
      END SUBROUTINE RCOMPX


      subroutine UPDSOL(JYEARS,JJDAYS)
      INTEGER, INTENT(IN) :: JYEARS,JJDAYS
      call SETSOL(JYEARS,JJDAYS,1)
      end subroutine UPDSOL

      SUBROUTINE SETSOL(JYEARS,JJDAYS,UPDSOL_flag)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C     SETSOL Parameters:
C----------------------
C            KSOLAR    Selects Solar Spectrum, (Lean vs Thekaekara Flux)
C            JYEARS    JYEAR Proxy:  Sets: Solar Constant Reference Year
C            JJDAYS    JDAY  Proxy:  Sets Reference Year Month JDAY/30.5
C                      (Nominal Reference: JYEARS= 1950 JJDAYS= January)
C
C-----------------------------------------------------------------------
C            KSOLAR   SOLSPEC     UVWAVLs       UVFACTs         KUVFAC
C-----------------------------------------------------------------------
C              -1     THEK      Can be set    Can be set   (if KUVFAC=1)
C-----------------------------------------------------------------------
C               0     SSI       Can be set    Can be set   (if KUVFAC=1)
C-----------------------------------------------------------------------
C               1     SSI       Can be set    Can be set   (if KUVFAC=1)
C-----------------------------------------------------------------------
C
C                               (Option to Modify Solar UV Fluxes)
C            UVWAVL    Specified Edges of UV Flux Variation SubIntervals
C            UVFACT    Factors to Change the Amplitude of UV Variability
C
C            KUVFAC    ON/OFF switch for activating UV Flux Modification
C            KSNORM    Re-Normalize S0 (VIS) (after UV Amplitude Change)
C                         (Nominal UVWAVLs are: 0.295,0.310,0.366)
C
C-----------------------------------------------------------------------
C     SETSOL Output:
C------------------
C
C        AO3 =   Ozone Absorption Table AO3(460)
C                (Solar UV Flux Weighted Absorption Table is used by the
C                FUNCTION AO3ABS(OCM) in SOLAR to compute Ozone Heating)
C                AO3 is the fraction of total Solar Flux absorbed by O3.
C
C     S00WM2 =   Solar Constant Reference Value for Time = JYEARS,JJDAYS
C                (Thekaekara, if KSOLAR=-1, Reference = 1367 WATTS/M**2)
C
C
C     SETSOL  is Generally Called once at Model Initialization to Select
C                Solar Flux (SSI,THEK), and to Define S00WM2 (RATLS0=1)
C
C-----------------------------------------------------------------------
C NOTE:
C-----
C     S00WM2 = Nominal Reference Solar Constant 1366.448785D0 WATTS/M**2
C                (Spectral Integral: Lean99 Solar Flux for January 1950)
C
C     KSOLAR=-1  Reproduces Thekaekhara Ozone Absorption, e.g., XRAD83XX
C     KSOLAR= 0  Uses Lean99 Solar Flux as set for Time= (JYEARS,JJDAYS)
C     KSOLAR= 1  Sets Lean99 Solar Flux to Current Time= (JYEARS,JJDAYS)
C     KSOLAR= 2  same as 1 but based on annual (not monthly) data
C                (JJDAYS used to select the specified Monthly-Mean Flux)
C     KSOLAR= 9  annual data for current time from file, but Thekaekhara
C                wavelength bins
C
C-----------------------------------------------------------------------
C
C     UPDSOL Parameters:
C----------------------
C            JYEARS    JYEAR Proxy:  Selects Solar Constant Current Year
C            JJDAYS    JDAY  Proxy:  Selects Lean Data Month JJDAYS/30.5
C
C     UPDSOL Output:
C------------------
C
C        AO3 =   Ozone Absorption Table AO3(460)
C                (Solar UV Flux Weighted Absorption Table is used by the
C                FUNCTION AO3ABS(OCM) in SOLAR to compute Ozone Heating)
C                AO3 is the fraction of total Solar Flux absorbed by O3.
C
C     RATLS0 =   Ratio:  Current-Time Solar Constant to Reference S00WM2
C
C-----------------------------------------------------------------------
C     Remark:
C
C     UPDSOL  is Called in RCOMPT to Update Solar Constant and Ozone AO3
C                Solar UV Absorption Dependence.  (Monthly-Mean Data are
C                NOT Interpolated in Time, but get Updated with Changing
C                Month, i.e., whenever JDAY/30.5 Reaches Integer Value.)
C
C-----------------------------------------------------------------------
      REAL*8, PARAMETER :: CORFAC=1366.2911D0/1366.4487855D0
      INTEGER, INTENT(IN) :: JYEARS,JJDAYS
      INTEGER, INTENT(IN), optional :: UPDSOL_flag
      INTEGER, SAVE :: LMOREF=0
      INTEGER JMO,LMO,Is0x,K,I,NWSUV,II,J,NUV,icyc
      REAL*8 FLXSUM,FFLUX(3),UVNORM,XX,OCM,TAUK,UVWAVA,UVWAVB,AO33

      if ( present(UPDSOL_flag) ) goto 777

C                                          Thekaekhara Solar Flux Option
C                                          -----------------------------
      IF(KSOLAR < 0) THEN
        WSOLAR(1:190)=WTHEK(1:190)
        FSOLAR(1:190)=FTHEK(1:190)
        S00WM2=1367.D0
        LMOREF=-111
        NWSUV=190
        GO TO 130
      END IF
C                                           Lean99 Solar Flux, UV Option
C                                           ----------------------------
      if(jyears > 2000) then
        icyc = icycs0f
      else
        icyc = icycs0
      end if
      if(Ksolar < 2) then    ! monthly data
        icyc = icyc*12
        JMO=1+JJDAYS/30.5D0
        IF(JMO > 12) JMO=12
        LMO=(JYEARS-iy1S0)*12+JMO
        IF(LMO > iMs0X) LMO=LMO-icyc*((LMO-iMs0X+icyc-1)/icyc)
        IF(LMO < 1) LMO=LMO+icyc*((icyc-lmo)/icyc)
      else                    ! annual data
        Is0x = nint( yr2s0-yr1s0+1 )
        lmo = nint( jyears - yr1s0 + 1.5 )
        IF(LMO > Is0X) LMO=LMO-icyc*((LMO-Is0X+icyc-1)/icyc)
        IF(LMO < 1) LMO=LMO+icyc*((icyc-lmo)/icyc)
      end if
      LMOREF=LMO

C                        IF(MADLUV==0) Default Option is then in force
C                        Default (FR_SSI) = Lean 1950 Jan Solar, UV flux
C                        CORFAC accounts for DS_SSI units in BLOCK DATA,
C                        and TSI1/TSI2 normalization of Lean input data.
C                        -----------------------------------------------

c      CORFAC=1366.2911D0/1366.4487855D0
      IF(KSOLAR.ne.9) THEN
        IF(MADLUV == 0) S00WM2 = SUM(FR_SSI(:)*DS_SSI(:)*CORFAC)
        IF(MADLUV >  0) S00WM2 = SUM(UV_SSI(LMO,:)*DS_SSI(:))
      ELSE
        S00WM2=TSI2(LMO)
      END IF

      IF(KSOLAR.ne.9) THEN
        I=0
        DO K=1,50
          I=I+1
          WSOLAR(I)=W1_SSI(K)
          IF(MADLUV == 0) FSOLAR(I)=FR_SSI(K)
          IF(MADLUV >  0) FSOLAR(I)=UV_SSI(LMO,K)
          I=I+1
          WSOLAR(I)=W1_SSI(K+1)
          FSOLAR(I)=FSOLAR(I-1)
        END DO
        NWSUV=100
      ELSE
        IF(MADLUV==0)call stop_model("invalid MADLUV for KSOLAR=9",255)
        WSOLAR(1:190)=WTHEK(1:190)
        FSOLAR(1:190)=UV_SSI(LMO,1:190)
        NWSUV=190
      END IF

C                                         Option to Modify Solar UV Flux
C                                         ------------------------------
  130 IF(KUVFAC==1) THEN
        FFLUX(:)=0.D0
        NUV=1
        DO I=1,NWSUV,2          ! by twos to account for histogram
 140      IF(WSOLAR(I+1) <= UVWAVL(NUV)) THEN
            FFLUX(NUV)=FFLUX(NUV)+FSOLAR(I)*(WSOLAR(I+1)-WSOLAR(I))
            FSOLAR(I:I+1)=FSOLAR(I:I+1)*UVFACT(NUV)
          ELSE
            NUV=NUV+1
            IF (NUV.LE.3) GOTO 140
            EXIT
          END IF
        END DO
        UVNORM = SUM(FFLUX(:)*(1d0-UVFACT(:)))
        IF (MADLUV==0) UVNORM=UVNORM*CORFAC
        IF (KSNORM==0) S00WM2=S00WM2-UVNORM
      ENDIF
C                  -----------------------------------------------------
C                  When KUVFAC=1 option multiplicative factors UVFACT(I)
C                  are used to change the UV spectral flux distribution,
C                  KSNORM=1 provides the option to keep S00WM2 constant.
C                  -----------------------------------------------------

      RATLS0=1.D0

      DO 190 I=1,460
      II=(I-10)/90-4
      XX=I-((I-10)/90)*90
      OCM=XX*10.D0**II
      DO 180 J=1,226
      TAUK=FUVKO3(J)*OCM
      IF(TAUK > 35.D0) TAUK=35.D0
      UVA(J)=1.D0-EXP(-TAUK)
  180 CONTINUE
      UVWAVA=0.100D0
      UVWAVB=0.400D0
      CALL FXGINT(UVA,XWAVO3,226,FSOLAR,WSOLAR,NWSUV,UVWAVA,UVWAVB,AO33)
      AO3(I)=AO33/S00WM2
  190 CONTINUE

C                       ------------------------------------------------
C                       NOTE:  AO3 is the Ozone-path Absorption Function
C                              AO3 convolves O3 asborption with solar UV
C                              spectral variations by FXGINT integration
C                              AO3 is expressed as the absorbed fraction
C                              of the total solar flux (S00WM2=1366W/m2)
C                              -----------------------------------------

      RETURN

C--------------------------------
!      ENTRY UPDSOL(JYEARS,JJDAYS)
C--------------------------------
 777  continue

      IF(KSOLAR < 1) RETURN   ! solar constant not time dependent
      IF(JYEARS < 1) RETURN   ! solar constant not time dependent

      if(jyears > 2000) then
        icyc = icycs0f
      else
        icyc = icycs0
      end if
      if(Ksolar == 1) then    ! monthly data
        icyc = icyc*12
        JMO=1+JJDAYS/30.5D0
        IF(JMO > 12) JMO=12
        LMO=(JYEARS-iy1S0)*12+JMO
        IF(LMO > iMs0X) LMO=LMO-icyc*((LMO-iMs0X+icyc-1)/icyc)
        IF(LMO < 1) LMO=LMO+icyc*((icyc-lmo)/icyc)
      else                    ! annual data        ksolar=2,9
        Is0x = nint( yr2s0-yr1s0+1 )
        lmo = nint( jyears - yr1s0 + 1.5 )
        IF(LMO > Is0X) LMO=LMO-icyc*((LMO-Is0X+icyc-1)/icyc)
        IF(LMO < 1) LMO=LMO+icyc*((icyc-lmo)/icyc)
      end if

      IF(LMO==LMOREF) RETURN  ! solar constant up-to-date
      LMOREF=LMO

C                                               Select Lean99 Solar Flux
C                                               ------------------------
      IF(KSOLAR.ne.9)THEN
        FLXSUM = SUM(UV_SSI(LMO,1:190)*DS_SSI(1:190))
      ELSE
        FLXSUM=TSI2(LMO)
      END IF
c        write(6,*) 'UPDSOLAR::FLXSUM::',FLXSUM

      IF(KSOLAR.ne.9) THEN
        I=0
        DO K=1,50
          I=I+1
          WSOLAR(I)=W1_SSI(K)
          FSOLAR(I)=UV_SSI(LMO,K)
          I=I+1
          WSOLAR(I)=W1_SSI(K+1)
          FSOLAR(I)=FSOLAR(I-1)
        END DO
        NWSUV=100
      ELSE
C                                          Select Thekaekhara Solar Flux
C                                          -----------------------------
        WSOLAR(1:190)=WTHEK(1:190)
        FSOLAR(1:190)=UV_SSI(LMO,1:190)
        NWSUV=190
      END IF
C                                         Option to Modify Solar UV Flux
C                                         ------------------------------
      IF(KUVFAC==1) THEN
        FFLUX(:)=0.D0
        NUV=1
        DO I=1,NWSUV,2          ! by twos to account for histogram
 240      IF(WSOLAR(I+1) <= UVWAVL(NUV)) THEN
            FFLUX(NUV)=FFLUX(NUV)+FSOLAR(I)*(WSOLAR(I+1)-WSOLAR(I))
            FSOLAR(I:I+1)=FSOLAR(I:I+1)*UVFACT(NUV)
          ELSE
            NUV=NUV+1
            IF (NUV.LE.3) GOTO 240
            EXIT
          END IF
        END DO
        UVNORM = SUM(FFLUX(:)*(1d0-UVFACT(:)))
        IF (MADLUV==0) UVNORM=UVNORM*CORFAC
        IF (KSNORM==0) FLXSUM=FLXSUM-UVNORM
      ENDIF

      RATLS0=FLXSUM/S00WM2

      DO 290 I=1,460
      II=(I-10)/90-4
      XX=I-((I-10)/90)*90
      OCM=XX*10.D0**II
      DO 280 J=1,226
      TAUK=FUVKO3(J)*OCM
      IF(TAUK > 35.D0) TAUK=35.D0
      UVA(J)=1.D0-EXP(-TAUK)
  280 CONTINUE
      UVWAVA=0.100D0
      UVWAVB=0.400D0
      CALL FXGINT(UVA,XWAVO3,226,FSOLAR,WSOLAR,NWSUV,UVWAVA,UVWAVB,AO33)
      AO3(I)=AO33/FLXSUM
  290 CONTINUE

      RETURN
      END SUBROUTINE SETSOL


      SUBROUTINE SETGHG(JYEARG,JJDAYG)
      IMPLICIT NONE
C
C
C     ---------------------------------------------------------------
C     SETGHG  Sets Default Greenhouse Gas Reference Year (for FULGAS)
C
C     Control Parameter:
C                     KTREND (specified in RADPAR) activates GH Trend
C               Default
C     KTREND  =   1
C     Selects   GTREND
C     ---------------------------------------------------------------
      INTEGER, INTENT(IN) :: JYEARG,JJDAYG
      REAL*8 TREF,TNOW
      INTEGER I
C
      TREF=JYEARG+(JJDAYG-0.999D0)/366.D0
C
      IF(KTREND==0) THEN
        XREF(1)=PPMV80(2)
        XREF(2)=PPMV80(6)
        XREF(3)=PPMV80(7)
        XREF(4)=PPMV80(8)*1000.D0
        XREF(5)=PPMV80(9)*1000.D0
        XREF(6)=PPMV80(11)*1000.D0  ! YREF11=PPMV80(11)*1000.D0
        XREF(7)=PPMV80(12)*1000.D0  ! ZREF12=PPMV80(12)*1000.D0
        RETURN
      END IF

      CALL GTREND(XREF,TREF)     ! finds xref 1-6 (yref11=xx6=xref(6))
      XREF(7)=1.D-13             ! ZREF12=1.D-13
      DO 120 I=1,NGHG
      IF(XREF(I) < 1.D-06) XREF(I)=1.D-06
  120 CONTINUE
      PPMV80(2)=XREF(1)
      PPMV80(6)=XREF(2)
      PPMV80(7)=XREF(3)
      PPMV80(8)=XREF(4)/1000.D0
      PPMV80(9)=XREF(5)/1000.D0
      PPMV80(11)=XREF(6)/1000.d0   ! YREF11/1000.D0
      PPMV80(12)=XREF(7)/1000.D0   ! ZREF12/1000.D0
      RETURN
      end SUBROUTINE SETGHG
C
C--------------------------------
!      ENTRY UPDGHG(JYEARG,JJDAYG)
C--------------------------------
      subroutine UPDGHG(JYEARG,JJDAYG)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JYEARG,JJDAYG
      REAL*8 TREF,TNOW
      INTEGER I
C
      TNOW=JYEARG+(JJDAYG-0.999D0)/366.D0
C
      IF(KTREND==0) THEN
        FULGAS(2)=PPMVK0(2)/XREF(1)
        FULGAS(6)=PPMVK0(6)/XREF(2)
        FULGAS(7)=PPMVK0(7)/XREF(3)
        FULGAS(8)=PPMVK0(8)/XREF(4)
        FULGAS(9)=PPMVK0(9)/XREF(5)
        FULGAS(11)=PPMVK0(11)/XREF(6) ! YREF11
        FULGAS(12)=PPMVK0(12)/XREF(7) ! .../ZREF12
        RETURN
      END IF

      CALL GTREND(XNOW,TNOW) ! finds xnow 1-6 (ynow11=xx6=xnow(6))
      XNOW(7)=1.D-20         ! ZNOW12=1.D-20
      FULGAS(2)=XNOW(1)/XREF(1)
      FULGAS(6)=XNOW(2)/XREF(2)
      FULGAS(7)=XNOW(3)/XREF(3)
      FULGAS(8)=XNOW(4)/XREF(4)
      FULGAS(9)=XNOW(5)/XREF(5)
      FULGAS(11)=XNOW(6)/XREF(6) ! YNOW11/YREF11
      FULGAS(12)=XNOW(7)/XREF(7) ! ZNOW12/ZREF12
C
      RETURN
      end subroutine UPDGHG

      subroutine GETGAS
      call SETGAS(1)
      end subroutine GETGAS

      SUBROUTINE SETGAS( GETGAS_flag)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Global   U.S. (1976) Standard Atmosphere  P, T, Geo Ht  Parameters
C-----------------------------------------------------------------------
      INTEGER, optional :: GETGAS_flag
      REAL*8, PARAMETER :: HPCON=34.16319d0,P0=1013.25d0,
     *     PI=3.141592653589793D0
      REAL*8, SAVE :: SINLAT(46)
      INTEGER, SAVE :: IFIRST=1, NL0
      INTEGER I,NLAY,NATM,L,J,K,N
      REAL*8 RHP,EST,FWB,FWT,PLT,DP,EQ,ES,ACM,HI,FI,HL,HJ,FJ,DH
     *     ,FF,GGVDF,ZT,ZB,EXPZT,EXPZB,PARTTR,PARTTG,PTRO,DL,DLS,DLN
     *     ,Z0LAT,SUMCOL,ULGASL,UGAS0(LX),UGASR(LX)

      if ( present(GETGAS_flag) ) goto 777

      IF(IFIRST==1) THEN
        SINLAT(:) = SIN(DLAT46(:)*PI/180.D0)
        NL0=NL
        IFIRST=0
      ENDIF
C                  -----------------------------------------------------
C                  Use PLB to fix Standard Heights for Gas Distributions
C                  -----------------------------------------------------

!nu   PS0=PLB0(1)

      DO 100 L=1,NL0
      DPL(L)=PLB0(L)-PLB0(L+1)
      PL(L)=(PLB0(L)+PLB0(L+1))*0.5D0
!nu   HLB(L)=HLB0(L)
  100 CONTINUE
!nu   HLB(NL0+1)=HLB0(NL0+1)
ccc      CALL RETERP(UFAC36,P36,36,FPXCO2,PL,NL0)
      CALL SET_FPXCO2(PL,FPXCO2,NL0,KFPCO2)
cc    IUFAC=1
cc    IF(IUFAC==0) FPXCO2(:)=1

      NLAY=LASTVC/100000
      NATM=(LASTVC-NLAY*100000)/10000
      IF(NATM > 0) GO TO 112

C     ----------------------------------------------------------------
C     Define Default Global Mean Gas Amounts for Off-Line Use Purposes
C
C     IGAS=1                              Global Mean H2O Distribution
C                                         ----------------------------
      RHP=0.77D0
      EST=10.D0**(9.4051D0-2353.D0/TLB(1))
      FWB=0.662D0*RHP*EST/(PLB0(1)-RHP*EST)
      DO 111 L=1,NL0
      PLT=PLB0(L+1)
      DP=PLB0(L)-PLT
      RHP=0.77D0*(PLT/P0-0.02D0)/.98D0
      EST=10.D0**(9.4051D0-2353.D0/TLT(L))
      FWT=0.662D0*RHP*EST/(PLT-RHP*EST)
      IF(FWT > 3.D-06) GO TO 110
      FWT=3.D-06
      RHP=FWT*PLT/(EST*(FWT+0.662D0))
  110 CONTINUE
      ULGASL=0.5D0*(FWB+FWT)*DP*h2o_mmr_to_cm_at_stp
      U0GAS(L,1)=ULGASL
      SHL(L)=ULGASL/(ULGASL+h2o_mmr_to_cm_at_stp*DP)
      EQ=0.5D0*(PLB0(L)+PLT)*SHL(L)/(0.662D0+0.378D0*SHL(L))
      ES=10.D0**(9.4051D0-2353.D0/TLM(L))
      RHL(L)=EQ/ES
      FWB=FWT
  111 CONTINUE
  112 CONTINUE

C                                         ----------------------------
C     IGAS=5                              Global Mean NO2 Distribution
C                                         ----------------------------
      ACM=0.D0
      HI=0.D0
      FI=CMANO2(1)
      HL=HLB0(2)
      L=1
      J=1
  130 CONTINUE
      J=J+1
      IF(J > 42) GO TO 133
      HJ=HI+2.D0
      FJ=CMANO2(J)
  131 CONTINUE
      DH=HJ-HI
      IF(HJ > HL) GO TO 132
      ACM=ACM+(FI+FJ)*DH*0.5D0
      HI=HJ
      FI=FJ
      GO TO 130
  132 CONTINUE
      FF=FI+(FJ-FI)*(HL-HI)/DH
      DH=HL-HI
      ACM=ACM+(FI+FJ)*DH*0.5D0
      U0GAS(L,5)=ACM
      ACM=0.D0
      HI=HL
      FI=FF
      IF(L==NL0) GO TO 133
      L=L+1
      HL=HLB0(L+1)
      GO TO 131
  133 CONTINUE
      U0GAS(L,5)=ACM
      ACM=0.D0
      L=L+1
      IF(L < NL0+1) GO TO 133
C                            -----------------------------------------
C     IGAS=2 and 4           (CO2,O2) Uniformly Mixed Gas Distribution
C                            -----------------------------------------
      DO 140 K=2,4,2
      U0GAS(1:NL0,K)=PPMV80(K)*ppmv_to_cm_at_stp*DPL(1:NL0)
  140 CONTINUE
C                -----------------------------------------------------
C     IGAS=6-12  (N20,CH4,F11,F12) Specified Vertical Gas Distribution
C                -----------------------------------------------------
      DO 151 K=6,12
      IF(K==10) GO TO 151
      DO 150 N=1,NL0
      GGVDF=1.D0-(1.D0-PPMVDF(K))*(1.D0-PLB0(N)/PLB0(1))
      IF(KGGVDF < 1) GGVDF=1.D0
      U0GAS(N,K)=PPMV80(K)*ppmv_to_cm_at_stp*DPL(N)*GGVDF
      ZT=(HLB0(N+1)-Z0(K))/ZH(K)
      IF(ZT <= 0.D0) GO TO 150
      ZB=(HLB0(N)-Z0(K))/ZH(K)
      EXPZT=EXP(-ZT)
      EXPZB=EXP(-ZB)
      IF(ZB < 0.D0) EXPZB=1.D0-ZB
      U0GAS(N,K)=U0GAS(N,K)*(EXPZB-EXPZT)/max(ZT-ZB,1d-6)
  150 CONTINUE
  151 CONTINUE
C                         --------------------------------------------
C                         Specification of  FULGAS  Scaled Gas Amounts
C                         --------------------------------------------

Cc*** Adjust water vapor in ALL layers                        ! IGAS=1
cc    ULGAS(1:NL0,1)=U0GAS(1:NL0,1)*FULGAS(1)
c**** Only adjust stratospheric levels (above LS1_loc)
      ULGAS(1:LS1_loc-1,1)=U0GAS(1:LS1_loc-1,1)
      ULGAS(LS1_loc:NL0,1)=U0GAS(LS1_loc:NL0,1)*FULGAS(1)
C****
      ULGAS(1:NL0,3)=U0GAS(1:NL0,3)*FULGAS(3)                 ! IGAS=3
      IF(KPFOZO==1) ULGAS(1:NL0,3)=ULGAS(1:NL0,3)*FPXOZO(1:NL0)

      DO 240 L=1,NL0                                     ! IGAS=2,4-13
!!!   PARTTR = (PLB(L)-PLB(L+1)) / (PLB0(L)-PLB0(L+1))   ! PLB=PLB0 ??
      DO 240 K=2,12
      IF(K==3) GO TO 240
!!!   PARTTG=PARTTR  ! next line not possible at this point (jlat=???)
!!!   IF(KPGRAD > 0) PARTTG=PARTTG*(1.D0+0.5D0*PPGRAD(K)*SINLAT(JLAT))
      ULGAS(L,K)=U0GAS(L,K)*FULGAS(K) !!! *PARTTG
  240 CONTINUE
      ULGAS(1:NL0,13)=U0GAS(1:NL0,13)*FULGAS(13)

      ULGAS(1:NL0,2)=ULGAS(1:NL0,2)*FPXCO2(1:NL0)

      RETURN


C-----------------
!      ENTRY GETGAS
C-----------------
 777  continue
C                        ---------------------------------------------
C                        Specify ULGAS: Get Gas Absorption from TAUGAS
C                        ---------------------------------------------

C                -----------------------------------------------------
C                N20,CH4,F11,F12 Specified Latitudinal Z0 Distribution
C                -----------------------------------------------------

      IF(KLATZ0 > 0) THEN
        PTRO=100.D0
        DL=DLAT46(JLAT)
        DLS=-40.D0
        DLN= 40.D0
        IF(DL < DLS) PTRO=189.D0-(DL+40.D0)*2.22D0
        IF(DL > DLN) PTRO=189.D0+(DL-40.D0)*2.22D0
        DO L=1,NL0
          IF(PLB0(L) >= PTRO) Z0LAT=HLB0(L)  ! orig. hlb not hlb0
        END DO
        DO 251 K=6,12
        IF(K==10) GO TO 251
        DO 250 L=1,NL0
        U0GAS(L,K)=PPMV80(K)*ppmv_to_cm_at_stp*(PLB0(L)-PLB0(L+1))
        IF(PLB0(1) >= PTRO) THEN ! safety check until P,H hard-coding removed
        ZT=(HLB0(L+1)-Z0LAT)/ZH(K)           ! orig. hlb not hlb0
        IF(ZT <= 0.D0) GO TO 250
        ZB=(HLB0(L)-Z0LAT)/ZH(K)             ! orig. hlb not hlb0
        EXPZT=EXP(-ZT)
        EXPZB=EXP(-ZB)
        IF(ZB < 0.D0) EXPZB=1.D0-ZB
        U0GAS(L,K)=U0GAS(L,K)*(EXPZB-EXPZT)/max(ZT-ZB,1d-6)
        ENDIF                    ! safety check
  250   CONTINUE
  251   CONTINUE
      ENDIF

      DO 300 L=L1,NL
      DPL(L)=PLB(L)-PLB(L+1)
      PL(L)=(PLB(L)+PLB(L+1))*0.5D0
  300 CONTINUE

      IF(KEEPRH==2) GO TO 313                  ! keep RH,SH
      IF(KEEPRH==1) GO TO 311                  ! find SH from RH
      DO 310 L=L1,NL                           ! find RH from SH
      EQ=PL(L)*SHL(L)/(0.662D0+0.378D0*SHL(L))
      ES=10.D0**(9.4051D0-2353.D0/TLM(L))
      RHL(L)=EQ/ES
  310 CONTINUE
      GO TO 313
  311 CONTINUE
      DO 312 L=L1,NL
      ES=10.D0**(9.4051D0-2353.D0/TLM(L))
      SHL(L)=0.622D0*(RHL(L)*ES)/(PL(L)-0.378D0*(RHL(L)*ES))
  312 CONTINUE
  313 CONTINUE

      U0GAS(L1:NL,1)=h2o_mmr_to_cm_at_stp*DPL(L1:NL)*
     *    SHL(L1:NL)/(1-SHL(L1:NL))
Cc*** Adjust water vapor in ALL layers
cc    ULGAS(L1:NL,1)=U0GAS(L1:NL,1)*FULGAS(1)
c**** Only adjust stratospheric levels (above LS1_loc)
      ULGAS(L1:LS1_loc-1,1)=U0GAS(L1:LS1_loc-1,1)
      ULGAS(LS1_loc:NL,1)=U0GAS(LS1_loc:NL,1)*FULGAS(1)
C****
      ULGAS(1:NL0,3)=U0GAS(1:NL0,3)*FULGAS(3)
      IF(KPFOZO==1) ULGAS(1:NL0,3)=ULGAS(1:NL0,3)*FPXOZO(1:NL0)

      DO 340 L=L1,NL0   ! =L1,NL for GCM use, =1,NL0 for offline use
      PARTTR = (PLB(L)-PLB(L+1)) / (PLB0(L)-PLB0(L+1))
      DO 339 K=2,12
      IF(K==3) GO TO 339
      PARTTG=PARTTR
      IF(KPGRAD > 0) PARTTG=PARTTG*(1.D0+0.5D0*PPGRAD(K)*SINLAT(JLAT))
      ULGAS(L,K)=U0GAS(L,K)*FULGAS(K)*PARTTG
  339 CONTINUE
      ULGAS(L,13)=U0GAS(L,13)*FULGAS(13)
  340 CONTINUE

      chem_out(:,4)=ULGAS(:,7) ! climatological CH4 saved for chemistry
      if(use_tracer_chem(2) > 0) ! allow use of tracer CH4.
     * ULGAS(1:use_tracer_chem(2),7)=chem_IN(2,1:use_tracer_chem(2))
#ifdef GCC_COUPLE_RAD
      GCCco2_out(:)=ULGAS(:,2) 
      if(use_tracer_GCCco2 > 0) then 
       ULGAS(1:use_tracer_GCCco2,2)=GCCco2_IN(1:use_tracer_GCCco2)
       ULGAS(use_tracer_GCCco2+1:NL,2)=GCCco2_IN(use_tracer_GCCco2)
      endif
#endif

      IF(MRELAY > 0) THEN          ! for offline use only
        IF(NO3COL > 0)             ! rescale ozone to col.amount RO3COL
     *    ULGAS(1:NL0,3) = U0GAS(1:NL0,3)*RO3COL/SUM( U0GAS(1:NL0,3) )
        DO 450 K=2,12                      ! repartition to new layering
          IF(K==10.and.KEEP10 > 0) GO TO 450
          UGAS0(1:NL0) = ULGAS(1:NL0,K)
          CALL REPART(UGAS0,PLB0,NL0+1, UGASR,PLB,NL+1)
          ULGAS(1:NL,K)=UGASR(1:NL)
  450   CONTINUE
        IF (KEEP10 > 0 .and. KEEP10 < 10)
     *    ULGAS(1:NL,KEEP10) = ULGAS(1:NL,10)
        IF (KEEP10 > 10)
     *    ULGAS(1:NL,KEEP10-10)=ULGAS(1:NL,KEEP10-10)+ULGAS(L,10)
      ENDIF

      if (NL > 40) then
        if(kfpco2.ge.3)
     *     call get_FPXCO2_105(FPXCO2(NL-38:NL),jlat,mlat46,jday)
      end if
      ULGAS(1:NL0,2)=ULGAS(1:NL0,2)*FPXCO2(1:NL0)

      chem_out(:,1)=ULGAS(:,3)!O3 considering move to RCOMPX; see above
C     chem_out(:,2)= _________              ! set in RCOMPX
      chem_out(:,3)=ULGAS(:,6)              ! N2O
C     chem_out(:,4)=ULGAS(:,7) ! CH4 (moved above before tracer option)
      chem_out(:,5)=ULGAS(:,8)+ULGAS(:,9)   ! CFC11(+)   +  CFC12(+)
      ! output CO2 in mole CO2 per mole air:
      CO2outCol(1:NL)=1.d-6*ULGAS(1:NL,2)/(ppmv_to_cm_at_stp*DPL(1:NL))

C-----------------
      CALL  TAUGAS
C-----------------

      RETURN
      END SUBROUTINE SETGAS


      subroutine GETO2A
      call SETO2A(1)
      end subroutine GETO2A

      SUBROUTINE SETO2A( GETO2A_flag )
      IMPLICIT NONE
      INTEGER,optional :: GETO2A_flag

      INTEGER, PARAMETER :: NW=18, NZ=11, NKO2=6
      REAL*8, PARAMETER ::
     *SFWM2(NW) = (/
     A     2.196E-3, 0.817E-3, 1.163E-3, 1.331E-3, 1.735E-3, 1.310E-3,
     B     1.311E-3, 2.584E-3, 2.864E-3, 4.162E-3, 5.044E-3, 6.922E-3,
     C     6.906E-3,10.454E-3, 5.710E-3, 6.910E-3,14.130E-3,18.080E-3/),

     *SIGMA(NW,NKO2) = RESHAPE( (/
     A     2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19,
     B     4.33E-21, 4.89E-21, 6.63E-21, 1.60E-20, 7.20E-20, 1.59E-18,
     C     2.10E-21, 2.32E-21, 3.02E-21, 6.30E-21, 3.46E-20, 7.52E-19,
     D     5.95E-22, 9.72E-22, 2.53E-21, 7.57E-21, 7.38E-20, 7.44E-19,
     E     3.33E-22, 1.02E-22, 4.09E-21, 1.63E-20, 8.79E-20, 3.81E-19,
     F     1.09E-21, 1.16E-21, 1.45E-21, 3.32E-21, 2.00E-20, 4.04E-19,
     G     1.15E-21, 1.30E-21, 1.90E-21, 4.89E-21, 2.62E-20, 4.08E-19,
     H     3.90E-22, 4.90E-22, 9.49E-22, 3.33E-21, 2.14E-20, 2.39E-19,
     I     1.29E-22, 2.18E-22, 8.28E-22, 3.46E-21, 1.94E-20, 1.06E-19,
     J     6.26E-23, 7.80E-23, 2.62E-22, 1.83E-21, 1.25E-20, 3.95E-20,
     K     2.74E-23, 3.58E-23, 8.64E-23, 4.03E-22, 2.13E-21, 1.95E-20,
     L     1.95E-23, 2.44E-23, 4.89E-23, 2.87E-22, 1.95E-21, 1.36E-20,
     M     1.84E-23, 1.96E-23, 2.71E-23, 8.52E-23, 6.48E-22, 3.89E-21,
     N     1.80E-23, 1.81E-23, 1.87E-23, 2.69E-23, 1.34E-22, 1.52E-21,
     O     1.80E-23, 1.80E-23, 1.82E-23, 2.40E-23, 5.71E-23, 5.70E-22,
     P     1.76E-23, 1.76E-23, 1.76E-23, 1.76E-23, 1.76E-23, 3.50E-23,
     Q     1.71E-23, 1.71E-23, 1.71E-23, 1.71E-23, 1.71E-23, 2.68E-23,
     R     1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23/)
     *     , (/ NW, NKO2 /) ),
     *WTKO2(NKO2) = (/0.05,0.20,0.25,0.25,0.20,0.05/),
     *STPMOL=2.68714D+19

      REAL*8 , SAVE :: ZTABLE(LX+1,11)
      INTEGER, SAVE :: NL0
      INTEGER, SAVE :: IFIRST=1
      REAL*8 FSUM,SUMMOL,ZCOS,WSUM,TAU,DLFLUX,WTI,WTJ
      INTEGER I,J,K,L,JI,JJ,N,LL

      if ( present(GETO2A_flag) ) goto 777

      if (mado2a==0) then
        ZTABLE(:,:)=0.
        return
      end if

      IF(IFIRST==1) THEN
        NL0=NL
        DO N=1,NL0
          ULGAS(N,4)=PPMV80(4)*ppmv_to_cm_at_stp*(PLB0(N)-PLB0(N+1))
        END DO
        IFIRST=0
      ENDIF

      FSUM = SUM(SFWM2(:))
      ZTABLE(NL0+1,:) = FSUM

      SUMMOL=0.D0
      DO 150 L=NL0,1,-1
      SUMMOL=SUMMOL+ULGAS(L,4)*STPMOL
      DO 140 J=1,NZ
      ZCOS=0.01D0*(1/J)+0.1D0*(J-1)
      FSUM=0.D0
      DO 130 I=1,NW
      WSUM=0.D0
      DO 120 K=1,NKO2
      TAU = SIGMA(I,K)*SUMMOL/ZCOS ; IF (TAU > 30) TAU=30
      WSUM=WSUM+WTKO2(K)*EXP(-TAU)
  120 CONTINUE
      FSUM=FSUM+WSUM*SFWM2(I)
  130 CONTINUE
      ZTABLE(L,J)=FSUM
  140 CONTINUE
  150 CONTINUE
      DO 160 J=1,NZ
      DO 160 L=1,NL0
      DLFLUX=ZTABLE(L+1,J)-ZTABLE(L,J)
      ZTABLE(L,J)=DLFLUX/1366.D0
  160 CONTINUE

      RETURN

C-----------------
!      ENTRY GETO2A
C-----------------
 777  continue

C              ---------------------------------------------------------
C              UV absorption by Oxygen is expressed as a fraction of the
C              total solar flux S0. Hence, O2FHRL(L)=ZTABLE(L,J) must be
C              normalized within SOLARM, dividing the GETO2A absorptions
C              O2FHRL(L) and O2FHRB(L) by the fraction of the solar flux
C              within the spectral interval DKS0(15), nominally by 0.05.
C              ---------------------------------------------------------
                          ! offline: may not yet work properly if NL>NL0
      ZCOS=1.D0+10.D0*COSZ
      JI=ZCOS ; IF (JI > 10) JI=10
      JJ=JI+1
      WTJ=ZCOS-JI
      WTI=1.0-WTJ
      O2FHRL(L1:NL) = WTI*ZTABLE(L1:NL,JI) + WTJ*ZTABLE(L1:NL,JJ)
      O2FHRB(L1:NL) =     ZTABLE(L1:NL,6)

      RETURN
      END SUBROUTINE SETO2A


      subroutine GETBAK
      call SETBAK(1)
      end subroutine GETBAK

      SUBROUTINE SETBAK( GETBAK_flag )
      IMPLICIT NONE
      INTEGER, optional :: GETBAK_flag
C     ------------------------------------------------------------------
C     SETBAK,GETBAK  Initializes Background Aerosol Specification, i.e.,
C                    Aerosol Composition and Distribution that is set in
C                    RADPAR by AGOLDH, BGOLDH, CGOLDH Factors
C                    and controlled by FGOLDH ON/OFF Scaling Parameters.
C                    Optional tracers may be added in SETAER/GETAER
C     ------------------------------------------------------------------
C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSBAER   FTBAER  ! Bgrnd Aerosol
C
C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW   (All-type) Aerosol Optical Depth
C                        FTTAER    LW   (All-type) Aerosol Optical Depth
C                        FSBAER    SW   SETBAKonly Aerosol Optical Depth
C                        FTBAER    LW   SETBAKonly Aerosol Optical Depth
C                        -----------------------------------------------

      REAL*8,  SAVE :: SRAX(LX,6,5),SRAS(LX,6,5),SRAC(LX,6,5)
      INTEGER, SAVE :: IFIRST=1
      INTEGER, SAVE :: NL0=0

      REAL*8 SGOLDH(5),TGOLDH(5),C,BC,ABC,HXPB,HXPT,ABCD
      INTEGER I,J,K,L,JJ

      if ( present(GETBAK_flag) ) goto 777

C**** Background aerosols
C     ------------------------------------------------------------------
C     Thermal: Set (5) Aerosol Type Compositions & Vertical Distribution
C     ------------------------------------------------------------------
      IF(IFIRST==1) THEN
        NL0=NL
        IFIRST=0
      ENDIF

      TRAX(:,:,:)=0                          ! 1:NL0,1:NKBAND,1:5

      DO 105 I=1,11
      DO 103 J=1,5
      IF(AGOLDH(I,J) < 1.D-06) GO TO 103
      C=CGOLDH(I,J)
      BC=EXP(-BGOLDH(I,J)/C)
      ABC=AGOLDH(I,J)*(1.D0+BC)

      HXPB=1.D0
      DO 102 L=1,NL0
      HXPT=HLB0(L+1)/C                       ! orig. hlb not hlb0
      IF(HXPT > 80.D0) GO TO 102
      HXPT=EXP(HXPT)
      ABCD=ABC/(1.D0+BC*HXPB)
     +    -ABC/(1.D0+BC*HXPT)
      HXPB=HXPT
      TRAX(L,:,J)=TRAX(L,:,J)+ABCD*(TRAQEX(:,I)-TRAQSC(:,I)) ! 1:NKBAND
  102 CONTINUE
  103 CONTINUE
      TRAQAB(:,I)=TRAQEX(:,I)-TRAQSC(:,I)
  105 CONTINUE

      TRBALK(:,:)=0                          ! 1:NL0,1:NKBAND

C-----------------------------------------------------------------------
C     SOLAR:   Set (5) Aerosol Type Compositions & Vertical Distribution
C-----------------------------------------------------------------------

      SRAX(:,:,:) = 1.D-20                   ! 1:NL0,1:6,1:5
      SRAS(:,:,:) = 1.D-30
      SRAC(:,:,:) = 0

      DO 114 I=1,11
      DO 113 J=1,5
      IF(AGOLDH(I,J) < 1.D-06) GO TO 113
      C=CGOLDH(I,J)
      BC=EXP(-BGOLDH(I,J)/C)
      ABC=AGOLDH(I,J)*(1.D0+BC)

      HXPB=1.D0
      DO 112 L=1,NL0
      HXPT=HLB0(L+1)/C                       ! orig. hlb not hlb0
      IF(HXPT > 80.D0) GO TO 112
      HXPT=EXP(HXPT)
      ABCD=ABC/(1.D0+BC*HXPB)
     +    -ABC/(1.D0+BC*HXPT)
      HXPB=HXPT
      SRAX(L,:,J) = SRAX(L,:,J) + ABCD*SRAQEX(:,I)
      SRAS(L,:,J) = SRAS(L,:,J) + ABCD*SRAQSC(:,I)
      SRAC(L,:,J) = SRAC(L,:,J) + ABCD*SRAQCB(:,I)*SRAQSC(:,I)
  112 CONTINUE
  113 CONTINUE
  114 CONTINUE

      SRAC(:,:,:) = SRAC(:,:,:)/SRAS(:,:,:)  ! 1:NL0,1:6,1:5

      SRBEXT(:,:) = 1.D-20                   ! 1:NL0,1:6
      SRBSCT(:,:) = 0
      SRBGCB(:,:) = 0
!nu   SRBPI0(:,:) = 0

      RETURN

C-----------------
!      ENTRY GETBAK
C-----------------
 777  continue
C     ------------------------------------------------------------------
C     GETBAK   Specifies Background Aerosol Contribution and Initializes
C                    (1) Thermal Radiation Aerosol Coefficient Table:
C                        TRAALK(L,K), for (L=1,NL), (K=1,33)
C
C                    (2) Solar Radiation Coefficient Tables:
C                        SRAEXT(L,K),SRASCT(L,K),SRAGCB(L,K) for (K=1,6)
C                    ---------------------------------------------------
C     Warning: MRELAY-section missing: not ready if NL.ne.NL0
C                                                              (Thermal)
C                                                              ---------
      TGOLDH(:)=FTTAER*FTBAER*FGOLDH(:)      ! 1:5
      DO 202 K=1,33
      DO 202 L=L1,NL0
      TRBALK(L,K) = SUM( TGOLDH(:)*TRAX(L,K,:) ) + 1.D-20
  202 CONTINUE

C                                                                (Solar)
C                                                                -------

      SGOLDH(:)=FSTAER*FSBAER*FGOLDH(:)      ! 1:5
      DO 212 K=1,6
      DO 212 L=L1,NL0
      SRBEXT(L,K) = SUM(SGOLDH(:)*SRAX(L,K,:)) + 1.D-20
      SRBSCT(L,K) = SUM(SGOLDH(:)*SRAS(L,K,:)) + 1.D-30
      SRBGCB(L,K) = SUM(SGOLDH(:)*SRAS(L,K,:)*SRAC(L,K,:)) / SRBSCT(L,K)
  212 CONTINUE


      RETURN
      END SUBROUTINE SETBAK


      subroutine GETAER
      call SETAER(1)
      end subroutine GETAER

      SUBROUTINE SETAER( GETAER_flag )
cc    INCLUDE  'rad00def.radCOMMON.f'
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      USE RESOLUTION, only :LM
#endif
      use AerParam_mod, only : DRYM2G
      use AerParam_mod, only : LMA
      IMPLICIT NONE
      INTEGER, optional :: GETAER_flag
C     ---------------------------------------------------------------
C     GISS MONTHLY-MEAN (1850-2050)  TROPOSPHERIC AEROSOL CLIMATOLOGY
C     ---------------------------------------------------------------

C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSAAER   FTAAER  ! AClim Aerosol

C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW   (All-type) Aerosol Optical Depth
C                        FTTAER    LW   (All-type) Aerosol Optical Depth
C                        FSAAER    SW   AClim Aer  Aerosol Optical Depth
C                        FTAAER    LW   AClim Aer  Aerosol Optical Depth
C                        -----------------------------------------------

!nu   DIMENSION ATAU09(9)
cc    DIMENSION PLBA09(10)          !       Aerosol data pressure levels
cc    DATA PLBA09/1010.,934.,854.,720.,550.,390.,255.,150.,70.,10./
      REAL*8, PARAMETER, dimension(4) ::
C              Crystallization RH               Deliquescence RH
     *  RHC=(/.38d0,.47d0,.28d0,.38d0/), RHD=(/.80d0,.75d0,.62d0,.80d0/)

C     ------------------------------------------------------------------
C                  Define aerosol size according to REFDRY specification
C                                      (if KRHAER(NA)=0, REFWET is used)
C                  FRSULF= Sulfate fraction of basic aerosol composition
C
C          Set size SO4 (NA=1) = Sulfate aerosol  (Nominal dry Reff=0.2)
C          Set size SEA (NA=2) = SeaSalt aerosol  (Nominal dry Reff=1.0)
C          Set size ANT (NA=3) = Nitrate aerosol  (Nominal dry Reff=0.3)
C          Set size OCX (NA=4) = Organic aerosol  (Nominal dry Reff=0.3)
C     ------------------------------------------------------------------
      REAL*8 AREFF, XRH,FSXTAU,FTXTAU,SRAGQL,RHFTAU,q55,RHDNA,RHDTNA
      REAL*8 ATAULX(LX,6),TTAULX(LX,ITRMAX),SRBGQL,FAC,RHFTAU_dry
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      REAL*8, DIMENSION(LM,6)  :: EXT,SCT,GCB
      REAL*8, DIMENSION(LM,33) :: TAB
#endif
      INTEGER NRHNAN(LX,8),K,L,NA,N,NRH,M,KDREAD,NT


#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      if (skip_AOD_in_rad) then ! rad_interact_aer>0
#ifdef TRACERS_AMP
        CALL SETAMP(EXT,SCT,GCB,TAB)
#endif
#ifdef TRACERS_TOMAS
        CALL SETTOMAS(EXT,SCT,GCB,TAB)
#endif
!radiation has 3 extra levels on the top - aerosols are zero
c SW
        SRBEXT(L1:LM,:) = EXT(L1:LM,:)
        SRBSCT(L1:LM,:) = SCT(L1:LM,:)
        SRBGCB(L1:LM,:) = GCB(L1:LM,:)
c LW
        TRBALK(L1:LM,:) = TAB(L1:LM,:)

        return ! nothing else to do here, everything handled in SETAMP/SETTOMAS
      endif
#endif

      if ( present(GETAER_flag) ) goto 200

      IF(MADAER <= 0) GO TO 150
      DO 110 NA=1,4
      AREFF=REFDRY(NA)
!nu   IF(KRHAER(NA) < 0) AREFF=REFWET(NA)
      CALL GETMIE(NA,AREFF,SRHQEX(1,1,NA),SRHQSC(1,1,NA),SRHQCB(1,1,NA)
     +                    ,TRHQAB(1,1,NA),Q55DRY(NA))
      DRYM2G(NA)=0.75D0/DENAER(NA)*Q55DRY(NA)/AREFF
!nu   IF(KRHAER(NA) < 0) DRYM2G(NA)=WETM2G(NA)
      RHINFO(1,1,NA)=0.D0                                     !  Rel Hum
      RHINFO(1,2,NA)=1.D0                                     !  TAUFAC
      RHINFO(1,3,NA)=AREFF                                    !  AerSize
      RHINFO(1,4,NA)=0.D0                                     !  LW g/m2
      RHINFO(1,5,NA)=1.33333333D0*AREFF*DENAER(NA)/Q55DRY(NA) !  Dryg/m2
      RHINFO(1,6,NA)=1.33333333D0*AREFF*DENAER(NA)/Q55DRY(NA) !  Totg/m2
      RHINFO(1,7,NA)=1.D0                                     !  Xmas fr
      RHINFO(1,8,NA)=DENAER(NA)                               !  Density
      RHINFO(1,9,NA)=Q55DRY(NA)                               !  Q55 Ext
  110 CONTINUE

C     Set size BCI (NA=5) = Black Carbon (Industrial) (Nominal Reff=0.1)
C     Set size BCB (NA=6) = Black Carbon (BioBurning) (Nominal Reff=0.1)
C     ------------------------------------------------------------------
      DO 120 NA=5,6
      AREFF=REFDRY(NA)
      CALL GETMIE(NA,AREFF,SRBQEX(1,NA),SRBQSC(1,NA),SRBQCB(1,NA)
     +                    ,TRBQAB(1,NA),Q55DRY(NA))
      DRYM2G(NA)=0.75D0/DENAER(NA)*Q55DRY(NA)/AREFF
  120 CONTINUE

              !      Extend default dry aerosol coefficients for N=2,190
      DO 135 N=2,190
      DO 135 NA=1,4
      SRHQEX(:,N,NA) = SRHQEX(:,1,NA)   !  1:6
      SRHQSC(:,N,NA) = SRHQSC(:,1,NA)   !  1:6
      SRHQCB(:,N,NA) = SRHQCB(:,1,NA)   !  1:6
      TRHQAB(:,N,NA) = TRHQAB(:,1,NA)   !  1:33
      RHINFO(N,1:9,NA) = RHINFO(1,1:9,NA)
  135 CONTINUE
                          !  Over-write dry coefficients if KRHAER(NA)=1
      KDREAD=71           !  default unit number for offline use only
      DO 140 NA=1,4
!nu   IF(KRHAER(NA) > 0) THEN
      CALL SETREL(REFDRY(NA),NA ,kdread
     A           ,SRUQEX,SRUQSC,SRUQCB
     B           ,TRUQEX,TRUQSC,TRUQCB
     C           ,REFU22,Q55U22,FRSULF
     D      ,SRHQEX(1,1,NA),SRHQSC(1,1,NA),SRHQCB(1,1,NA)
     E      ,TRHQAB(1,1,NA)
     F      ,RHINFO(1,1,NA))
!nu   ENDIF
  140 CONTINUE

  150 CONTINUE
      IF(NTRACE <= 0) RETURN

C**** Optional Tracer aerosols initializations
      DO NT=1,NTRACE
      NA=ITR(NT)
      AREFF=TRRDRY(NT)
      CALL GETMIE(NA,AREFF,SRTQEX(1,1,NT),SRTQSC(1,1,NT),SRTQCB(1,1,NT)
     +                    ,TRTQAB(1,1,NT),Q55)
      RTINFO(1,1,NT)=0.0
      RTINFO(1,2,NT)=1.0
      RTINFO(1,3,NT)=AREFF
      RTINFO(1,4,NT)=0.0
      RTINFO(1,5,NT)=1.33333333D0*AREFF*DENAER(NA)/Q55
      RTINFO(1,6,NT)=1.33333333D0*AREFF*DENAER(NA)/Q55
      RTINFO(1,7,NT)=1.0
      RTINFO(1,8,NT)=DENAER(NA)
      RTINFO(1,9,NT)=Q55
      END DO
            !      Define default dry aerosol coefficients for N=2,190
      DO N=2,190
      DO NT=1,NTRACE
        SRTQEX(:,N,NT) = SRTQEX(:,1,NT)    ! 1:6
        SRTQSC(:,N,NT) = SRTQSC(:,1,NT)    ! 1:6
        SRTQCB(:,N,NT) = SRTQCB(:,1,NT)    ! 1:6
        TRTQAB(:,N,NT) = TRTQAB(:,1,NT)    ! 1:33
        RTINFO(N,1:9,NT) = RTINFO(1,1:9,NT)
      END DO
      END DO
                          !  Over-write dry coefficients if KRHTRA(NT)=1
      KDREAD=71           !  default unit number for offline use only
      DO NT=1,NTRACE
      NA=ITR(NT)
      IF (KRHTRA(NT) > 0 .and. NA <= 4) THEN
      CALL SETREL(TRRDRY(NT),NA,KDREAD
     A           ,SRUQEX,SRUQSC,SRUQCB
     B           ,TRUQEX,TRUQSC,TRUQCB
     C           ,REFU22,Q55U22,FRSULF
     D      ,SRTQEX(1,1,NT),SRTQSC(1,1,NT),SRTQCB(1,1,NT)
     E      ,TRTQAB(1,1,NT)
     F      ,RTINFO(1,1,NT))
      ENDIF
      END DO

      RETURN


C-----------------
!      ENTRY GETAER
C-----------------
  200 continue

      NRHNAN(:,:) = 1
      DO 230 L=L1,NL
      if (RHL(L) > 0.9005D0) then
        XRH = (RHL(L)-0.899499D0)*1000.D0
        NRH = XRH+90 ; if (NRH > 189) NRH=189
      else
        XRH=RHL(L)*100.D0+0.5D0
        NRH=XRH  ;   if (NRH < 0) NRH=0
      endif
      DO 220 NA=1,4
      if (KDELIQ(L,NA)==0) then
        RHDNA = RHD(NA) ; if (KRHDTK==1) RHDNA = RHDTNA(TLM(L),NA)
        if (RHL(L) > RHDNA)   KDELIQ(L,NA)=1
      else
        if (RHL(L) < RHC(NA)) KDELIQ(L,NA)=0
      endif
      NRHNAN(L,NA)=NRH*KDELIQ(L,NA)+1
  220 CONTINUE
  230 CONTINUE

      IF(MADAER <= 0) GO TO 500

      DO NA=1,6
      IF(MADAER.eq.3) THEN
      CALL REPART (A6JDAY(1,NA,IGCM,JGCM),PLBAER,lma+1,    ! in
#ifdef REPART_AER_FIX
      ! passing plb0 instead of plb for approximate consistency with input
     *             ATAULX(1,NA),PLB0,NL+1)              ! out
#else
     *             ATAULX(1,NA),PLB,NL+1)               ! out
#endif
      ELSE
      CALL REPART (A6JDAY(1,NA,ILON,JLAT),PLBA09,10,    ! in
     *             ATAULX(1,NA),PLB,NL+1)               ! out
      ENDIF
      END DO

      FSXTAU=FSTAER*FSAAER+1.D-10
      FTXTAU=FTTAER*FTAAER
                           !            (Solar BCI,BCB components)
      DO 250 L=L1,NL
      nintaerext(L,:,5)=SRBQEX(:,5)*ATAULX(L,5)*FSXTAU*FS8OPX(5)
      nintaerext(L,:,6)=SRBQEX(:,6)*ATAULX(L,6)*FSXTAU*FS8OPX(6)
      nintaersca(L,:,5)=SRBQSC(:,5)*ATAULX(L,5)*FSXTAU*FS8OPX(5)
      nintaersca(L,:,6)=SRBQSC(:,6)*ATAULX(L,6)*FSXTAU*FS8OPX(6)
      nintaerasy(L,:,5)=SRBQCB(:,5)
      nintaerasy(L,:,6)=SRBQCB(:,6)
      SRAEXT(L,:)=nintaerext(L,:,5)+nintaerext(L,:,6)
      SRASCT(L,:)=nintaersca(L,:,5)+nintaersca(L,:,6)
      SRAGCB(L,:)=(nintaersca(L,:,5)*nintaerasy(L,:,5)
     &            +nintaersca(L,:,6)*nintaerasy(L,:,6))
     /            /(SRASCT(L,:)+1.D-10)
  250 CONTINUE
                          !           (Thermal BCI,BCB components)
      DO 260 L=L1,NL
      TRAALK(L,:) = TRBQAB(:,5)*ATAULX(L,5)*FTXTAU*FT8OPX(5) +   ! 1:33
     +              TRBQAB(:,6)*ATAULX(L,6)*FTXTAU*FT8OPX(6)
      IF(PLB(L) > 10) GO TO 260
      TRAALK(L,:)=0
  260 CONTINUE

      DO 330 NA=1,4
      DO 330 L=L1,NL
      RHFTAU=RHINFO(NRHNAN(L,NA),2,NA)*ATAULX(L,NA)*FSXTAU*FS8OPX(NA)
      DO 330 K=1,6
      nintaerext(L,K,NA)=SRHQEX(K,NRHNAN(L,NA),NA)*RHFTAU
      nintaersca(L,K,NA)=SRHQSC(K,NRHNAN(L,NA),NA)*RHFTAU
      nintaerasy(L,K,NA)=SRHQCB(K,NRHNAN(L,NA),NA)
      SRAEXT(L,K)=SRAEXT(L,K)+nintaerext(L,K,NA)
      SRAGQL     =SRAGCB(L,K)*SRASCT(L,K)
     &           +nintaerasy(L,K,NA)*nintaersca(L,K,NA)
      SRASCT(L,K)=SRASCT(L,K)+nintaersca(L,K,NA)
      SRAGCB(L,K)=SRAGQL/(SRASCT(L,K)+1.D-10)
  330 CONTINUE

      DO 360 NA=1,4
      DO 360 L=L1,NL
      RHFTAU=RHINFO(NRHNAN(L,NA),2,NA)*ATAULX(L,NA)*FTXTAU*FT8OPX(NA)
      TRAALK(L,:)=TRAALK(L,:)+TRHQAB(:,NRHNAN(L,NA),NA)*RHFTAU    ! 1:33
  360 CONTINUE

  500 CONTINUE

      IF(NTRACE <= 0) RETURN

C     ------------------------------------------------------------------
C     Option to add on Tracer Type aerosol thermal & solar contributions
C
C     NOTE:  Aerosol carried as a tracer is assumed to be in kg/m2 units
C     ------------------------------------------------------------------

      DO NT=1,NTRACE
        IF (ITR(NT) == 7) THEN
          FAC = 1d3*.75d0/TRADEN(NT)*RTINFO(1,9,NT)/TRRDRY(NT)
        ELSE
          FAC = 1d3*.75d0/DENAER(ITR(NT))*Q55DRY(ITR(NT))/TRRDRY(NT)
        END IF
        TTAULX(L1:NL,NT) = TRACER(L1:NL,NT) * FAC
      END DO

      FSXTAU=FSTAER*FSBAER+1.D-10
      FTXTAU=FTTAER*FTBAER

      DO 700 NT=1,NTRACE
      NA=ITR(NT)
      DO 700 L=L1,NL
      RHFTAU=RTINFO(NRHNAN(L,NA),2,NT)*TTAULX(L,NT)*FSXTAU
      RHFTAU_dry=RTINFO(1,2,NT)*TTAULX(L,NT)*FSXTAU
      IF (FSTOPX(NT) > 0) THEN
        RHFTAU=RHFTAU*FSTOPX(NT)*FSTASC(NT)
        RHFTAU_dry=RHFTAU_dry*FSTOPX(NT)*FSTASC(NT)
        DO K=1,6
          SRBEXT(L,K)=SRBEXT(L,K)+SRTQEX(K,NRHNAN(L,NA),NT)*RHFTAU
          SRBGQL =SRBGCB(L,K)*SRBSCT(L,K)+SRTQCB(K,NRHNAN(L,NA),NT)
     +           *SRTQSC(K,NRHNAN(L,NA),NT)*RHFTAU
          SRBSCT(L,K)=SRBSCT(L,K)+SRTQSC(K,NRHNAN(L,NA),NT)*RHFTAU
          SRBGCB(L,K)=SRBGQL/(SRBSCT(L,K)+1.D-10)
        END DO
      END IF
      aesqex(L,:,nt)=srtqex(:,nrhnan(L,na),nt)*rhftau           ! 1:6
      aesqsc(L,:,nt)=srtqsc(:,nrhnan(L,na),nt)*rhftau
      aesqcb(L,:,nt)=srtqcb(:,nrhnan(L,na),nt)*aesqsc(L,:,nt)
      aesqex_dry(L,:,nt)=srtqex(:,1,nt)*rhftau_dry           ! 1:6
      aesqsc_dry(L,:,nt)=srtqsc(:,1,nt)*rhftau_dry
      aesqcb_dry(L,:,nt)=srtqcb(:,1,nt)*aesqsc_dry(L,:,nt)
  700 CONTINUE

      DO 750 NT=1,NTRACE
      NA=ITR(NT)
      DO 750 L=L1,NL
      RHFTAU=RTINFO(NRHNAN(L,NA),2,NT)*TTAULX(L,NT)*FTXTAU*FTTOPX(NT)
     *  *FTTASC(NT)
      TRBALK(L,:)=TRBALK(L,:)+TRTQAB(:,NRHNAN(L,NA),NT)*RHFTAU ! 1:33
  750 CONTINUE

      RETURN
      END SUBROUTINE SETAER

C-----------------
!      ENTRY GETDST
C-----------------
      subroutine GETDST
C     ---------------------------------------------------------------
C     MONTHLY-MEAN DESERT DUST CLIMATOLOGY
C     ---------------------------------------------------------------

C           OUTPUT: via SRDEXT(L,K)   D Dust Extinction Optical Depth
C                       SRDSCT(L,K)   D Dust Scattering Optical Depth
C                       SRDGCB(L,K)   D Dust Asymmetry Parameter  g
C                       TRDALK(L,K)  Thermal Absorption Optical Depth

C     Tau Scaling Factors:    Solar    Thermal    apply to:
C                             FSTAER   FTTAER  ! Total Aerosol
C                             FSDAER   FTDAER  ! Dust  Aerosol
C
C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW   (All-type) Aerosol Optical Depth
C                        FTTAER    LW   (All-type) Aerosol Optical Depth
C                        FSDAER    SW   Dust Aer   Aerosol Optical Depth
C                        FTDAER    LW   Dust Aer   Aerosol Optical Depth
C                        -----------------------------------------------
      use DustParam_mod
      IMPLICIT NONE
      REAL*8 FSXTAU,FTXTAU,DTAULX(LX+1,nsized) !ron
      INTEGER K,L,N
      real*8 :: TDUST_col(lmd)

      if(.not.dust_optics_initialized) then
        dust_optics_initialized = .true.
        allocate( QXDUST(6,nsized), QSDUST(6,nsized), QCDUST(6,nsized),
     *       ATDUST(33,nsized), QDST55(nsized) )

        allocate(taucon_dust(nsized))
        DO N=1,nsized
          CALL GETMIE(7,REDUST(N),QXDUST(1,N),QSDUST(1,N),QCDUST(1,N)
     +         ,ATDUST(1,N),QDST55(N))
          ! save the factor for converting from concentration to AOT
          TAUCON_dust(N)=0.75E+03*QDST55(N)/(RODUST(N)*REDUST(N))
        ENDDO
      endif

      DO N=1,nsized
        TDUST_col(:) = DDJDAY(:,N,IGCM,JGCM)*taucon_dust(n) ! kg/m2 -> tau
#ifdef REPART_AER_FIX
        ! passing plb0 instead of plb for approximate consistency with input
        CALL REPART(TDUST_col,PLBdust,lmd+1,DTAULX(1,N),PLB0,NL+1)
#else
        CALL REPART(TDUST_col,PLBdust,lmd+1,DTAULX(1,N),PLB,NL+1)
#endif
      ENDDO

C                     Apply Solar/Thermal Optical Depth Scaling Factors
C                              Dust Aerosol  Solar   FSXD=FSTAER*FSDAER
C                              Dust Aerosol Thermal  FTXD=FSTAER*FTDAER
C                              ----------------------------------------

      FSXTAU=FSTAER*FSDAER+1.D-10
      FTXTAU=FTTAER*FTDAER

      DO 220 K=1,6
      DO 210 L=L1,NL
      SRDEXT(L,K)=2.D-10
      SRDSCT(L,K)=1.D-10
      SRDGCB(L,K)=0.D0
  210 CONTINUE
  220 CONTINUE

      DO 270 L=L1,NL
      DO 260 K=1,6
      nintaerext(L,K,7)=sum( QXDUST(K,:)*DTAULX(L,:) )*FSXTAU*FS8OPX(7)
      nintaersca(L,K,7)=sum( QSDUST(K,:)*DTAULX(L,:) )*FSXTAU*FS8OPX(7)
      SRDEXT(L,K)= SRDEXT(L,K) + nintaerext(L,K,7)
      SRDSCT(L,K)= SRDSCT(L,K) + nintaersca(L,K,7)
      nintaerasy(L,K,7)=sum( QCDUST(K,:)*QSDUST(K,:)*DTAULX(L,:) )*
     &                  FSXTAU*FS8OPX(7)/(SRDSCT(L,K)+1.D-10)
      SRDGCB(L,K)= nintaerasy(L,K,7)
  260 CONTINUE
  270 CONTINUE

      DO 280 L=L1,NL
      DO 280 K=1,33
      TRDALK(L,K)= sum (ATDUST(K,:)*DTAULX(L,:)*FTXTAU*FT8OPX(7)) ! 1:nsized !ron
  280 CONTINUE

      RETURN
      end subroutine GETDST


      subroutine UPDVOL(JYEARV,JDAYVA)
      INTEGER, INTENT(IN) :: JYEARV,JDAYVA
      call SETVOL(JYEARV,JDAYVA)
      end subroutine UPDVOL

      subroutine GETVOL
      call SETVOL(GETVOL_flag=1)
      end subroutine GETVOL

      SUBROUTINE SETVOL(JYEARV,JDAYVA,GETVOL_flag)
      IMPLICIT NONE


      REAL*8, SAVE :: E46LAT(47),SIZLAT(46),TAULAT(46)
      INTEGER, SAVE :: NJ46
      REAL*8, PARAMETER :: HVOL00(5) = (/15.0, 20.0, 25.0, 30.0, 35.0/)
cx    INTEGER, SAVE :: LATVOL = 0   ! not ok for grids finer than 72x46

!nu   REAL*8, PARAMETER :: htplim=1.d-3
      REAL*8, SAVE :: FSXTAU,FTXTAU
      INTEGER, INTENT(IN), optional :: JYEARV,JDAYVA,GETVOL_flag
      INTEGER J,L,MI,MJ,K
      REAL*8 XYYEAR,XYI,WMI,WMJ,SIZVOL !nu ,SUMHTF
      REAL*8, save, allocatable :: gdata(:), hlattf(:), HTFLAT(:,:)
      REAL*8, save, allocatable :: HTPROF(:)
#ifdef HEALY_LM_DIAGS
      INTEGER, SAVE :: NJDG
      REAL*8, SAVE :: EDGLAT(JM_DIAG)
#endif

C     ------------------------------------------------------------------
C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSVAER   FTVAER  ! SETVOL Aer

C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW  (All-type) Aerosol Optical Depth
C                        FTTAER    LW  (All-type) Aerosol Optical Depth
C                        FSVAER    SW  SETVOLonly Aerosol Optical Depth
C                        FTVAER    LW  SETVOLonly Aerosol Optical Depth
C                        -----------------------------------------------

C     -----------------------------------------------------------------
C     VEFF0   Selects Size Distribution Variance (this affects Thermal)
C     REFF0   Selects Effective Particle Size for Archive Volcanic Data
C     -----------------------------------------------------------------

      if ( present(JYEARV) ) goto 777 ! UPDVOL
      if ( present (GETVOL_flag) ) goto 778 ! GETVOL

      FSXTAU=FSTAER*FSVAER
      FTXTAU=FTTAER*FTVAER

C                   Set Grid-Box Edge Latitudes for Data Repartitioning
C                   ---------------------------------------------------
      if(madvol==1) then
        allocate(ELATVOL(NVolLat+1))
        DO J=2,24 ! NVolLat
          ELATVOL(J)=-90.D0+(J-1.5D0)*180.D0/23.D0
        END DO
        ELATVol( 1)=-90.D0
        ELATVol(25)= 90.D0
        HVolKM = HVol00
      end if
      NJ46=46+1
      DO J=2,46
        E46LAT(J)=-90.D0+(J-1.5D0)*180.D0/(MLAT46-1)
      END DO
      E46LAT(   1)=-90.D0
      E46LAT(NJ46)= 90.D0
#ifdef  HEALY_LM_DIAGS
      NJDG=JM_DIAG+1
      DO J=2,JM_DIAG
      EDGLAT(J)=-90.D0+(J-1.5D0)*180.D0/(JM_DIAG-1)
      END DO
      EDGLAT(   1)=-90.D0
      EDGLAT(NJDG)= 90.D0
#endif
      allocate (gdata(NVolLat),hlattf(NVolK),HTFLAT(NJ46,NVOLK))
      allocate (HTPROF(NL))

      HTPROF(:)=0

C                       -----------------------------------------------
C                       Initialize H2SO4 Q,S,C,A Tables for Input VEFF0
C                       -----------------------------------------------
C     ------------------
      CALL SETQVA(VEFF0)
C     ------------------

      RETURN


C--------------------------------
!      ENTRY UPDVOL(JYEARV,JDAYVA)
C--------------------------------
 777  continue

C                                          (Volcanic data)
C                                          -------------------------
      XYYEAR=JYEARV+JDAYVA/366.D0
      IF(XYYEAR < JVOLYI) XYYEAR=JVOLYI
      XYI=(XYYEAR-JVOLYI)*12.D0+1.D0
      IF(XYI > NVOLMON - .001D0) XYI=NVOLMON-.001D0
!!    write(6,'(a,2f9.1,3i7)') 'VOLCYEAR=',
!!   .   XYI,XYYEAR,JVOLYI,JYEARV,JDAYVA
      MI=XYI
      WMJ=XYI-MI
      WMI=1.D0-WMJ
      MJ=MI+1
      DO 250 J=1,NVolLat
      GDATA(J)=WMI*VReffTJ(MI,J)+WMJ*VReffTJ(MJ,J)
!!    write(6,'(a,2I7,2f8.1,2f10.4)')'VOLCREFF:: ',MI,MJ,
!!   . XYYEAR,XYI,VReffTJ(MI,J),VReffTJ(MJ,J)
  250 CONTINUE
      CALL RETERP(GDATA,ELATVol,NVolLat+1,SIZLAT,E46LAT,NJ46)
      DO 270 K=1,NVOLK
      DO 260 J=1,NVolLat
      GDATA(J)=WMI*VTauTJK(MI,J,K)+WMJ*VTauTJK(MJ,J,K)
!!    write(6,'(a,3I7,2f8.1,2F10.4)')'VOLCAER:: ',K,MI,MJ,
!!   . XYYEAR,XYI,VTauTJK(MI,J,K),VTauTJK(MJ,J,K)
  260 CONTINUE
      CALL RETERP(GDATA,ELATVOL,NVolLat+1,HTFLAT(1,K),E46LAT,NJ46)
  270 CONTINUE
#ifdef HEALY_LM_DIAGS
      DO J=1,46
      TAULAT(J) = SUM (HTFLAT(J,:))
      END DO
      CALL RETERP(TAULAT,E46LAT,NJ46,VTAULAT,EDGLAT,NJDG)
#endif


      RETURN


C-----------------
!      ENTRY GETVOL
C-----------------
 778   continue
cx    IF(MRELAY > 0)    GO TO 300
cx    IF(JLAT==LATVOL) GO TO 350  ! not ok for grids finer than 72x46

C                      Set JLAT Dependent Aerosol Distribution and Size
C                      ------------------------------------------------
cx300 CONTINUE

      HLATTF(1:NVolK)=HTFLAT(JLAT,1:NVolK)
      CALL REPART(HLATTF,HVOLKM,NVolK+1,HTPROF,HLB0,NL+1)
!nu   LHPMAX=0       ! not used
!nu   LHPMIN=NL      ! not used
!nu   DO L=L1,NL
!nu     N=NL+1-L
!nu     IF(HTPROF(L) >= HTPLIM) LHPMAX=L
!nu     IF(HTPROF(N) >= HTPLIM) LHPMIN=N
!nu   END DO
!nu   SUMHTF=1.D-10
      DO 330 L=L1,NL
      IF(HTPROF(L) < 0.) HTPROF(L)=0.D0
!nu   SUMHTF=SUMHTF+HTPROF(L)
  330 CONTINUE

      SIZVOL=SIZLAT(JLAT)

C                        Select H2SO4 Q,S,C,A Tables for  Size = SIZVOL
C                        ----------------------------------------------

C------------------------
      CALL GETQVA(SIZVOL)
C------------------------

cx    LATVOL=JLAT
cx350 CONTINUE
C                                  ------------------------------------
C                                  H2SO4 Thermal Contribution in TRVALK
C                                  ------------------------------------
      DO 420 K=1,33
      TRVALK(L1:NL,K)=HTPROF(L1:NL)*AVH2S(K)*FTXTAU*FT8OPX(8)
  420 CONTINUE

C                      H2SO4 Solar Contribution in SRVEXT,SRVSCT,SRVGCB
C                      ------------------------------------------------

      DO 440 K=1,6
      nintaerext(L1:NL,K,8)=QVH2S(K)*HTPROF(L1:NL)*FSXTAU*FS8OPX(8)
      nintaersca(L1:NL,K,8)=SVH2S(K)*HTPROF(L1:NL)*FSXTAU*PIVMAX
     &                     *FS8OPX(8)
      nintaerasy(L1:NL,K,8)=GVH2S(K)
      SRVEXT(L1:NL,K)=nintaerext(L1:NL,K,8)
      SRVSCT(L1:NL,K)=nintaersca(L1:NL,K,8)
      SRVGCB(L1:NL,K)=nintaerasy(L1:NL,K,8)
  440 CONTINUE

      RETURN
      END SUBROUTINE SETVOL


      subroutine GETQVA(SIZVOL)
      REAL*8, INTENT(IN) :: SIZVOL
      call SETQVA(SIZVOL=SIZVOL)
      end subroutine GETQVA

      SUBROUTINE SETQVA(VEFF,SIZVOL)
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SETQVA   Selects (interpolates) H2SO4 Mie Parameters for specified
C              Variance VEFF for subsequent Size interpolation by GETQVA
C     ------------------------------------------------------------------

ceq   REAL*8 SRQV( 6,20),SRSV( 6,20),SRGV( 6,20),Q55V(   20),REFV(20)
ceq   REAL*8 TRQV(33,20),TRSV(33,20),TRGV(33,20),TRAV(33,20),VEFV(20)
      REAL*8 TRAB(33,20),Q5(5),RV20(20),QV20(20)
      REAL*8, PARAMETER ::  V5(5)=(/ .1d0, .2d0, .3d0, .4d0, .5d0/)
      SAVE TRAB

C     ------------------------------------------------------------------
C     SRVQEX Volcanic Aerosol sizes (Reff) range from 0.1 to 5.0 microns
C     To utilize equal interval interpolation, Reff N=9,20 are redefined
C     so Volcanic Aerosol sizes have effective range of 0.1-2.0 microns.
C     ------------------------------------------------------------------
      REAL*8, INTENT(IN), optional :: SIZVOL,VEFF
      REAL*8 REFN,RADX,WTJHI,WTJLO
      INTEGER I,K,N,JRXLO,JRXHI

      if ( present(SIZVOL) ) goto 777

      DO 130 N=1,20
      RV20(N)=REFV20(N,1)
      VEFV(N)=VEFF
      REFV(N)=N/10.D0
      Q5(:) = Q55V20(N,:5)
      CALL SPLINE(V5,Q5,5,VEFF,Q55V(N),1.D0,1.D0,1)
      DO 115 K=1,6
      Q5(:) = SRVQEX(K,N,:5)
      CALL SPLINE(V5,Q5,5,VEFF,SRQV(K,N),1.D0,1.D0,1)
      Q5(:) = SRVQSC(K,N,:5)
      CALL SPLINE(V5,Q5,5,VEFF,SRSV(K,N),1.D0,1.D0,1)
      Q5(:) = SRVQCB(K,N,:5)
      CALL SPLINE(V5,Q5,5,VEFF,SRGV(K,N),1.D0,1.D0,1)
  115 CONTINUE
      DO 120 K=1,33
      Q5(:) = TRVQEX(K,N,:5)
      CALL SPLINE(V5,Q5,5,VEFF,TRQV(K,N),1.D0,1.D0,1)
      Q5(:) = TRVQSC(K,N,:5)
      CALL SPLINE(V5,Q5,5,VEFF,TRSV(K,N),1.D0,1.D0,1)
      Q5(:) = TRVQCB(K,N,:5)
      CALL SPLINE(V5,Q5,5,VEFF,TRGV(K,N),1.D0,1.D0,1)
      Q5(:) = TRVQAL(K,N,:5)
      CALL SPLINE(V5,Q5,5,VEFF,TRAV(K,N),1.D0,1.D0,1)
  120 CONTINUE
      TRAB(:,N) = TRQV(:,N)-TRSV(:,N)    ! 1:33
  130 CONTINUE

      QV20(:) = Q55V(:)                  ! 1:20
      DO 132 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,Q55V(N),1.D0,1.D0,1)
  132 CONTINUE
      DO 140 K=1,6
      QV20(:) = SRQV(K,:)
      DO 134 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,SRQV(K,N),1.D0,1.D0,1)
  134 CONTINUE
      QV20(:) = SRSV(K,:)
      DO 136 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,SRSV(K,N),1.D0,1.D0,1)
  136 CONTINUE
      QV20(:) = SRGV(K,:)
      DO 138 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,SRGV(K,N),1.D0,1.D0,1)
  138 CONTINUE
  140 CONTINUE
      DO 150 K=1,33
      QV20(:) = TRQV(K,:)
      DO 142 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRQV(K,N),1.D0,1.D0,1)
  142 CONTINUE
      QV20(:) = TRSV(K,:)
      DO 144 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRSV(K,N),1.D0,1.D0,1)
  144 CONTINUE
      QV20(:) = TRGV(K,:)
      DO 146 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRGV(K,N),1.D0,1.D0,1)
  146 CONTINUE
      QV20(:) = TRAV(K,:)
      DO 148 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRAV(K,N),1.D0,1.D0,1)
  148 CONTINUE
      DO 149 N=9,20
      TRAB(K,N)=TRQV(K,N)-TRSV(K,N)
  149 CONTINUE
  150 CONTINUE

      RETURN

C-------------------------
!      ENTRY GETQVA(SIZVOL)
C-------------------------
 777  continue
C     ------------------------------------------------------------------
C     Volcanic Aerosol sizes have effective range of  0.1 - 2.0 microns.
C     ------------------------------------------------------------------

      RADX=SIZVOL*10.D0
      IF(RADX < 1.000001D0) RADX=1.000001D0
      IF(RADX > 19.99999D0) RADX=19.99999D0
      JRXLO=RADX
      WTJHI=RADX-JRXLO
      WTJLO=1.D0-WTJHI
      JRXHI=JRXLO+1

      QVH2S(:) = WTJLO*SRQV(:,JRXLO) + WTJHI*SRQV(:,JRXHI)   ! 1:6
      SVH2S(:) = WTJLO*SRSV(:,JRXLO) + WTJHI*SRSV(:,JRXHI)
      GVH2S(:) = WTJLO*SRGV(:,JRXLO) + WTJHI*SRGV(:,JRXHI)

      Q55H2S=WTJLO*Q55V(JRXLO)+WTJHI*Q55V(JRXHI)

      AVH2S(:) = WTJLO*TRAB(:,JRXLO) + WTJHI*TRAB(:,JRXHI)   ! 1:33

      RETURN
      END SUBROUTINE SETQVA

      SUBROUTINE SETCLD
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Control Parameters used in SETCLD,GETCLD,GETEPS: defined in RADPAR
C
C             ICE012  Selects Water, Non-Mie, Mie  Ice Cloud Qex,Qsc,Pi0
C             TAUWC0  Minimum Optical Depth for Water Clouds
C             TAUIC0  Minimum Optical Depth for   Ice Clouds
C             FCLDTR  Scaling Factor for Thermal Cloud Optical Depth
C             FCLDSR  Scaling Factor for  Solar  Cloud Optical Depth
C             EPSCON  Column Cloud Inhomogeneity EPSILON (when KCLDEP=1)
C             KCLDEP  Selects Cloud Inhomogeneity Option (0-4):
C                     KCLDEP =  0  Sets Column CLDEPS to Zero
C                     KCLDEP =  1  Sets Column CLDEPS to EPSCON
C                     KCLDEP =  2  Keeps whatever is specified in CLDEPS
C                     KCLDEP =  3  Uses: Column EPCOL(72,46) Climatology
C                     KCLDEP =  4  Uses: Ht Dep EPLOW, EPMID, EPHIG Data
C
C-----------------------------------------------------------------------
C                    Define Cloud Absorption Cross-Sections
C
C     Selected by:   ICE012 = 0    Liquid Water Droplets   (N =  1 -  5)
C                    ICE012 = 1    Ice - Non-Spherical     (N =  6 - 10)
C                    ICE012 = 2    Ice - Mie (Spherical)   (N = 11 - 15)
C
C     Define Solar,Thermal Cloud Single Scattering Albedo: SRCQPI( 6,15)
C                                                          TRCQPI(33,15)
C-----------------------------------------------------------------------

      TRCQAB(:,:)=TRCQEX(:,:)-TRCQSC(:,:)   ! 1:33,1:15
      TRCQPI(:,:)=TRCQSC(:,:)/TRCQEX(:,:)

      SRCQPI(:,:)=SRCQSC(:,:)/SRCQEX(:,:)   ! 1:6,1:15

C                          Initialize  GETCLD Output Parameters to Zero
C                          --------------------------------------------
      TRCTCA(:)=0          ! 1:33
      TRCALK(:,:)=0        ! 1:NL,1:33
      SRCEXT(:,:)=1.D-20   ! 1:NL,1:6
      SRCSCT(:,:)=0
      SRCGCB(:,:)=0

      RETURN
      end SUBROUTINE SETCLD
C-----------------
!      ENTRY GETCLD
C-----------------
      subroutine GETCLD
      IMPLICIT NONE
      REAL*8 SIZWCL,SIZICL,XRW,XMW,XPW,EPS,VEP,VEP1,VEP2,VEPP,TAUWCL
     *     ,TAUICL,QAWATK,QPWATK,SRCGFW,QXWATK,QSWATK,QGWATK,XRI,XMI,XPI
     *     ,QAICEK,QPICEK,SRCGFC,QXICEK,QSICEK,QGICEK,SCTTAU,GCBICE
     *     ,SCTGCB,TCTAUW,TCTAUC,ALWATK,WTI,WTW,ALICEK,TRCTCI,XJDAY,XMO
     *     ,WTMJ,WTMI
      INTEGER I,J,N,K,L,LBOTCW,LTOPCW,LBOTCI,LTOPCI,IRWAT,IRICE,MI,MJ

C-----------------------------------------------------------------------
C           Define:    TRCALK(LX,33)  Thermal Radiation Cloud Absorption
C                      TRCTCA(33)     Thermal Radiation Top Cloud Albedo
C
C                      SRCEXT(LX,6)   Solar Radiation Cloud Ext Op Depth
C                      SRCSCT(LX,6)   Solar Radiation Cloud Sct Op Depth
C                      SRCGCB(LX,6)   Solar Radiation Cloud Asym Param g
C
C                         LTOPCL      Top Cloud Layer Location
C                         LBOTCL      Bot Cloud Layer Location
C
C                         LTOPCW      Top Water Cloud Layer Location
C                         LBOTCW      Bot Water Cloud Layer Location
C
C                         LTOPCI      Top Ice Cloud Layer Location
C                         LBOTCI      Bot Ice Cloud Layer Location
C
C-----------------------------------------------------------------------

      LBOTCW=0
      LTOPCW=0
      LBOTCI=0
      LTOPCI=0
      TRCTCA(:)=0              ! 1:33
      DO 280 L=L1,NL
      TRCALK(L,:)=0
      SRCEXT(L,:)=1.D-20       ! 1:6
      SRCSCT(L,:)=1.D-30
      SRCGCB(L,:)=0
      SRCPI0(L,:)=0
C                                         Water Cloud Size Interpolation
C                                         ------------------------------

      IF(FTAUC*TAUWC(L) > TAUWC0) THEN
      SIZWCL=SIZEWC(L)
      LTOPCW=L
      IF(LBOTCW==0) LBOTCW=L
      IF(SIZWCL < 15.D0) THEN
      IF(SIZWCL < 3.0D0) SIZWCL=3.0D0
      IRWAT=2
      XRW=SIZWCL/10.0D0-1.00D0
      ELSE
      IF(SIZWCL > 25.D0) SIZWCL=25.D0
      IRWAT=4
      XRW=SIZWCL/10.0D0-2.00D0
      ENDIF
      XMW=1.D0-XRW-XRW
      XPW=1.D0+XRW+XRW
      EPS=CLDEPS(L)
      VEP=EPS/(1.D0-EPS)
      VEP1=1.D0+VEP
      TAUWCL=FTAUC*TAUWC(L)
      DO 240 K=1,33
      QAWATK=XMW*XPW*TRCQAB(K,IRWAT)
     +      -XMW*XRW*TRCQAB(K,IRWAT-1)+XPW*XRW*TRCQAB(K,IRWAT+1)
      QPWATK=XMW*XPW*TRCQPI(K,IRWAT)
     +      -XMW*XRW*TRCQPI(K,IRWAT-1)+XPW*XRW*TRCQPI(K,IRWAT+1)
      VEPP=VEP*QPWATK
      TRCALK(L,K)=TRCALK(L,K)+TAUWCL*QAWATK/(VEP1-VEPP)
  240 CONTINUE
      SRCGFW=SRCGSF(1)
      DO 250 K=1,6
      QXWATK=XMW*XPW*SRCQEX(K,IRWAT)
     +      -XMW*XRW*SRCQEX(K,IRWAT-1)+XPW*XRW*SRCQEX(K,IRWAT+1)
      QSWATK=XMW*XPW*SRCQSC(K,IRWAT)
     +      -XMW*XRW*SRCQSC(K,IRWAT-1)+XPW*XRW*SRCQSC(K,IRWAT+1)
      QGWATK=XMW*XPW*SRCQCB(K,IRWAT)
     +      -XMW*XRW*SRCQCB(K,IRWAT-1)+XPW*XRW*SRCQCB(K,IRWAT+1)
      QPWATK=XMW*XPW*SRCQPI(K,IRWAT)
     +      -XMW*XRW*SRCQPI(K,IRWAT-1)+XPW*XRW*SRCQPI(K,IRWAT+1)
      QGWATK=QGWATK*SRCGFW
      VEPP=VEP*QPWATK
      VEP2=VEP1-VEPP
      SRCEXT(L,K)=SRCEXT(L,K)+TAUWCL*QXWATK/VEP1
      SRCSCT(L,K)=TAUWCL*QSWATK/(VEP1*VEP2)
      SRCGCB(L,K)=QGWATK*VEP2/(VEP1-VEPP*QGWATK)
  250 CONTINUE
      ENDIF
C                                           Ice Cloud Size Interpolation
C                                           ----------------------------
      IF(FTAUC*TAUIC(L) > TAUIC0) THEN
      SIZICL=SIZEIC(L)
      LTOPCI=L
      IF(LBOTCI==0) LBOTCI=L
      IF(SIZICL < 25.D0) THEN
      IF(SIZICL < 3.0D0) SIZICL=3.0D0
      IRICE=2+ICE012*5
      XRI=SIZICL/20.D0-0.75D0
      ELSE
      IF(SIZICL > 75.D0) SIZICL=75.D0
      IRICE=4+ICE012*5
      XRI=SIZICL/50.D0-1.00D0
      ENDIF
      XMI=1.D0-XRI-XRI
      XPI=1.D0+XRI+XRI
      EPS=CLDEPS(L)
      VEP=EPS/(1.D0-EPS)
      VEP1=1.D0+VEP
      TAUICL=FTAUC*TAUIC(L)
      DO 260 K=1,33
      QAICEK=XMI*XPI*TRCQAB(K,IRICE)
     +      -XMI*XRI*TRCQAB(K,IRICE-1)+XPI*XRI*TRCQAB(K,IRICE+1)
      QPICEK=XMI*XPI*TRCQPI(K,IRICE)
     +      -XMI*XRI*TRCQPI(K,IRICE-1)+XPI*XRI*TRCQPI(K,IRICE+1)
      VEPP=VEP*QPICEK
      TRCALK(L,K)=TRCALK(L,K)+TAUICL*QAICEK/(VEP1-VEPP)
  260 CONTINUE

      SRCGFC=SRCGSF(2)
      IF(ICE012==2) SRCGFC=SRCGSF(3)
      DO 270 K=1,6
      QXICEK=XMI*XPI*SRCQEX(K,IRICE)
     +      -XMI*XRI*SRCQEX(K,IRICE-1)+XPI*XRI*SRCQEX(K,IRICE+1)
      QSICEK=XMI*XPI*SRCQSC(K,IRICE)
     +      -XMI*XRI*SRCQSC(K,IRICE-1)+XPI*XRI*SRCQSC(K,IRICE+1)
      QGICEK=XMI*XPI*SRCQCB(K,IRICE)
     +      -XMI*XRI*SRCQCB(K,IRICE-1)+XPI*XRI*SRCQCB(K,IRICE+1)
      QPICEK=XMI*XPI*SRCQPI(K,IRICE)
     +      -XMI*XRI*SRCQPI(K,IRICE-1)+XPI*XRI*SRCQPI(K,IRICE+1)
      QGICEK=QGICEK*SRCGFC
      VEPP=VEP*QPICEK
      VEP2=VEP1-VEPP
      SRCEXT(L,K)=SRCEXT(L,K)+TAUICL*QXICEK/VEP1
      SCTTAU=TAUICL*QSICEK/(VEP1*VEP2)
      GCBICE=QGICEK*VEP2/(VEP1-VEPP*QGICEK)
      SCTGCB=SRCSCT(L,K)*SRCGCB(L,K)+SCTTAU*GCBICE
      SRCSCT(L,K)=SRCSCT(L,K)+SCTTAU
      SRCGCB(L,K)=SCTGCB/SRCSCT(L,K)
  270 CONTINUE
      ENDIF
  280 CONTINUE

C     ------------------------------------------------------------------
C     Identify Top Cloud (LTOPCL) and define top cloud albedo correction
C
C     Full Scattering Correction:      KCLDEM=1   ECLTRA=1.0   (default)
C     Partial(rad99a) Correction:      KCLDEM=0   ECLTRA=1.0
C       No Scattering Correction:      KCLDEM=0   ECLTRA=0.0
C
C     KCLDEM=1 Top-cloud scattering correction uses TXCTPG,TSCTPG,TGCTPG
C              to generate correction (over-rides old ECLTRA correction)
C              (KCLDEM correction is computed in THRMAL at LTOPCL level)
C     ------------------------------------------------------------------

      LTOPCL=LTOPCI
      IF(LTOPCI > LTOPCW) GO TO 330
      IF(LTOPCW < 1) GO TO 350
      LTOPCL=LTOPCW
      TCTAUW=FTAUC*TAUWC(LTOPCL)
      DO 310 K=1,33
      ALWATK=XMW*XPW*TRCQAL(K,IRWAT)
     +      -XMW*XRW*TRCQAL(K,IRWAT-1)+XPW*XRW*TRCQAL(K,IRWAT+1)
      QXWATK=XMW*XPW*TRCQEX(K,IRWAT)
     +      -XMW*XRW*TRCQEX(K,IRWAT-1)+XPW*XRW*TRCQEX(K,IRWAT+1)
      TRCTCA(K)=(1.D0-EXP(-FTAUC*TAUWC(LTOPCL)*QXWATK))*ALWATK*ECLTRA
      QSWATK=XMW*XPW*TRCQSC(K,IRWAT)
     +      -XMW*XRW*TRCQSC(K,IRWAT-1)+XPW*XRW*TRCQSC(K,IRWAT+1)
      QGWATK=XMW*XPW*TRCQCB(K,IRWAT)
     +      -XMW*XRW*TRCQCB(K,IRWAT-1)+XPW*XRW*TRCQCB(K,IRWAT+1)
      TXCTPG(K)=QXWATK*TCTAUW
      TSCTPG(K)=QSWATK*TCTAUW
      TGCTPG(K)=QGWATK
  310 CONTINUE
      LBOTCL=LBOTCW
      IF(LBOTCI < 1) GO TO 360
      IF(LBOTCI <= LBOTCW) LBOTCL=LBOTCI
      IF(LTOPCI==LTOPCW) THEN
      TCTAUW=FTAUC*TAUWC(LTOPCL)
      TCTAUC=FTAUC*TAUIC(LTOPCL)
      WTI=TAUIC(LTOPCL)/(TAUIC(LTOPCL)+TAUWC(LTOPCL))
      WTW=TAUWC(LTOPCL)/(TAUIC(LTOPCL)+TAUWC(LTOPCL))
      DO 320 K=1,33
      ALICEK=XMI*XPI*TRCQAL(K,IRICE)
     +      -XMI*XRI*TRCQAL(K,IRICE-1)+XPI*XRI*TRCQAL(K,IRICE+1)
      QXICEK=XMI*XPI*TRCQEX(K,IRICE)
     +      -XMI*XRI*TRCQEX(K,IRICE-1)+XPI*XRI*TRCQEX(K,IRICE+1)
      TRCTCI=(1.D0-EXP(-FTAUC*TAUIC(LTOPCL)*QXICEK))*ALICEK*ECLTRA
      TRCTCA(K)=WTW*TRCTCA(K)+WTI*TRCTCI
      QSICEK=XMI*XPI*TRCQSC(K,IRICE)
     +      -XMI*XRI*TRCQSC(K,IRICE-1)+XPI*XRI*TRCQSC(K,IRICE+1)
      QGICEK=XMI*XPI*TRCQCB(K,IRICE)
     +      -XMI*XRI*TRCQCB(K,IRICE-1)+XPI*XRI*TRCQCB(K,IRICE+1)
      TXCTPG(K)=TXCTPG(K)+QXICEK*TCTAUC
      SCTGCB=TSCTPG(K)*TGCTPG(K)+QSICEK*TCTAUC*QGICEK
      TSCTPG(K)=TSCTPG(K)+QSICEK*TCTAUC
      TGCTPG(K)=SCTGCB/(1.D-10+TSCTPG(K))
  320 CONTINUE
      ENDIF
      GO TO 360
  330 CONTINUE
      LTOPCL=LTOPCI
      TCTAUC=FTAUC*TAUIC(LTOPCL)
      DO 340 K=1,33
      ALICEK=XMI*XPI*TRCQAL(K,IRICE)
     +      -XMI*XRI*TRCQAL(K,IRICE-1)+XPI*XRI*TRCQAL(K,IRICE+1)
      QXICEK=XMI*XPI*TRCQEX(K,IRICE)
     +      -XMI*XRI*TRCQEX(K,IRICE-1)+XPI*XRI*TRCQEX(K,IRICE+1)
      TRCTCA(K)=(1.D0-EXP(-FTAUC*TAUIC(LTOPCL)*QXICEK))*ALICEK*ECLTRA
      QSICEK=XMI*XPI*TRCQSC(K,IRICE)
     +      -XMI*XRI*TRCQSC(K,IRICE-1)+XPI*XRI*TRCQSC(K,IRICE+1)
      QGICEK=XMI*XPI*TRCQCB(K,IRICE)
     +      -XMI*XRI*TRCQCB(K,IRICE-1)+XPI*XRI*TRCQCB(K,IRICE+1)
      TXCTPG(K)=QXICEK*TCTAUC
      TSCTPG(K)=QSICEK*TCTAUC
      TGCTPG(K)=QGICEK
  340 CONTINUE
      LBOTCL=LBOTCI
      IF(LBOTCW==0) GO TO 360
      LBOTCL=LBOTCW
      GO TO 360
  350 CONTINUE
      LBOTCL=0
      LTOPCL=0
  360 CONTINUE

      RETURN
      end subroutine GETCLD

C--------------------------------
!      ENTRY UPDEPS(JYEARE,JJDAYE)
C--------------------------------
      subroutine UPDEPS(JYEARE,JJDAYE)
C                 Select ISCCP-Based Cloud Heterogeneity Time Dependence
C                 ------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JYEARE,JJDAYE
      REAL*8 XJDAY,XMO,WTMJ,WTMI
      INTEGER MI,MJ

      XJDAY=JJDAYE-0.999D0
      XMO=XJDAY/30.5D0+.5D0
      MI=XMO
      WTMJ=XMO-MI
      WTMI=1.D0-WTMJ
      IF(MI < 1) MI=12
      MJ=MI+1
      IF(MJ > 12) MJ=1

      EPLOW(:,:) = WTMI*EPLMHC(:,:,MI,1) + WTMJ*EPLMHC(:,:,MJ,1) ! 72,46
      EPMID(:,:) = WTMI*EPLMHC(:,:,MI,2) + WTMJ*EPLMHC(:,:,MJ,2)
      EPHIG(:,:) = WTMI*EPLMHC(:,:,MI,3) + WTMJ*EPLMHC(:,:,MJ,3)
      EPCOL(:,:) = WTMI*EPLMHC(:,:,MI,4) + WTMJ*EPLMHC(:,:,MJ,4)

      RETURN
      end subroutine UPDEPS

C-----------------
!      ENTRY GETEPS
C-----------------
      subroutine GETEPS
C             ----------------------------------------------------------
C                     Select Cloud Heterogeneity CLDEPS Options
C             EPSCON  Column Cloud Inhomogeneity EPSILON (when KCLDEP=1)
C             KCLDEP  Selects Cloud Inhomogeneity Option (0-4):
C                     KCLDEP =  0  Sets Column CLDEPS to Zero
C                     KCLDEP =  1  Sets Column CLDEPS to EPSCON
C                     KCLDEP =  2  Keeps whatever is specified in CLDEPS
C                     KCLDEP =  3  Uses: Column EPCOL(72,46) Climatology
C                     KCLDEP =  4  Uses: Ht Dep EPLOW, EPMID, EPHIG Data
C                     --------------------------------------------------
      IMPLICIT NONE
      INTEGER L

      IF(KCLDEP == 0)  CLDEPS(L1:NL) = 0
      IF(KCLDEP == 1)  CLDEPS(L1:NL) = EPSCON
      IF(KCLDEP == 3)  CLDEPS(L1:NL) = EPCOL(ILON,JLAT)
      IF(KCLDEP == 4)  then
        DO L=L1,NL
          CLDEPS(L) = EPMID(ILON,JLAT)
          IF(PLB(L) > 750)  CLDEPS(L) = EPLOW(ILON,JLAT)
          IF(PLB(L) < 430)  CLDEPS(L) = EPHIG(ILON,JLAT)
        END DO
      ENDIF

      RETURN
      end subroutine GETEPS

      SUBROUTINE TAUGAS
      IMPLICIT NONE
C     -------------------------------------------------------------
C     TAUGAS INPUT REQUIRES:  L1,NL,PL,DPL,TLM,ULGAS, TAUTBL,TAUWV0
C                             TAUCD0,TAUO30, XKCFC,H2OCN8,H2OCF8
C                             ULOX,DUX,XTRUP,XTU0,XTRDN,XTD0
C                             DXUP2,DXDN2,DXUP3,DXDN3,DXUP6,DXDN6
C                             DXUP7,DXDN7,DXUP8,DXDN8,DXUP9,DXDN9
C                             DXUP13,DXDN13
C     TAUGAS OUTPUT DATA IS:  TRGXLK,XTRU,XTRD
C     ----------------------------------------------------------

      INTEGER, PARAMETER :: NPU2=14, NPU=5
      REAL*8, PARAMETER :: TLOX=181.d0, DTX=23.d0, P0=1013.25d0

      REAL*8, PARAMETER :: PX(NPX)= (/1000d0, 750d0, 500d0, 300d0,
     *        200d0, 100d0, 50d0, 20d0, 10d0,   5d0,   2d0,   1d0,
     *         .5d0,  .2d0, .1d0,.03d0,.01d0,.003d0,.001d0/)

      INTEGER, PARAMETER :: NGX(4) = (/12,12, 8,33/),
     *                     IG1X(4) = (/ 2,14,26, 1/)
      REAL*8, PARAMETER :: PDPU2(NPU2) = (/1.d4, 1.d5,2.d5,5.d5,
     *           1.d6,2.d6,5.d6, 1.d7,2.d7,5.d7, 1.d8,2.d8,5.d8, 1.d9/)
      REAL*8, PARAMETER ::  PU(NPU) = (/  50.,200.,800.,3200.,12800./)
      INTEGER, PARAMETER :: IGASX(21) = (/ 1, 2, 3, 1, 1, 2, 2, 3, 3, 6,
     *     6, 6, 7, 7,13,13, 8, 8, 9, 9, 1/)
      INTEGER, PARAMETER :: KGX(21) =   (/ 1, 2, 3, 2, 3, 1, 3, 1, 2, 1,
     *     2, 3, 1, 3, 1, 3, 2, 3, 2, 3, 4/)
      INTEGER, PARAMETER :: NUX(16) = (/25, 9, 9, 9, 9, 5, 5, 5, 5, 2, 2
     *     ,2 , 2, 2, 2,2/)
      INTEGER, PARAMETER :: IGUX(16) = (/ 0,300,408,480,588,660,720,760
     *     ,820,880,904,928,944,968,984,1008/)


      REAL*8, PARAMETER :: XKH2OW(8) = (/
     & .432D-5,.943D-5,.188D-4,.352D-4,.623D-4,.105D-3,.170D-3,.262D-3
     &     /)

      REAL*8, PARAMETER ::  XKCFCW(8,2) = RESHAPE( (/
     + 11.0, 11.7, 11.5, 10.9, 10.3, 9.90, 9.90, 9.90,
     + 5.75, 5.72, 5.95, 5.95, 5.90, 6.51, 6.51, 6.51 /)
     *     , (/8,2/) )

      REAL*8, PARAMETER ::  PCF(NLCF) = (/
     & 0.98981D+03,0.96840D+03,0.94446D+03,0.91796D+03,0.88891D+03,
     & 0.85579D+03,0.81757D+03,0.77425D+03,0.72686D+03,0.67692D+03,
     & 0.62545D+03,0.57296D+03,0.52098D+03,0.47104D+03,0.42365D+03,
     & 0.37932D+03,0.33855D+03,0.30186D+03,0.26874D+03,0.23867D+03,
     & 0.21115D+03,0.18567D+03,0.16172D+03,0.13900D+03,0.11800D+03,
     & 0.99000D+02,0.81500D+02,0.65000D+02,0.50000D+02,0.37000D+02,
     & 0.25500D+02,0.15000D+02,0.78100D+01,0.43900D+01,0.24700D+01,
     & 0.13900D+01,0.78100D+00,0.43900D+00,0.24700D+00,0.13900D+00,
     & 0.75000D-01,0.35000D-01,0.10000D-01/)

      REAL*8, PARAMETER ::  DPCF(NLCF) = (/
     & 0.20380D+02,0.22430D+02,0.25470D+02,0.27520D+02,0.30580D+02,
     & 0.35670D+02,0.40770D+02,0.45860D+02,0.48920D+02,0.50960D+02,
     & 0.51980D+02,0.53000D+02,0.50960D+02,0.48920D+02,0.45860D+02,
     & 0.42810D+02,0.38720D+02,0.34660D+02,0.31590D+02,0.28540D+02,
     & 0.26500D+02,0.24460D+02,0.23440D+02,0.22000D+02,0.20000D+02,
     & 0.18000D+02,0.17000D+02,0.16000D+02,0.14000D+02,0.12000D+02,
     & 0.11000D+02,0.10000D+02,0.43800D+01,0.24600D+01,0.13800D+01,
     & 0.78000D+00,0.43800D+00,0.24600D+00,0.13800D+00,0.78000D-01,
     & 0.50000D-01,0.30000D-01,0.20000D-01/)

      REAL*8, PARAMETER ::  PLBCF(NLCF+1) = (/
     & 0.10000D+04,0.97962D+03,0.95719D+03,0.93172D+03,0.90420D+03,
     & 0.87362D+03,0.83795D+03,0.79718D+03,0.75132D+03,0.70240D+03,
     & 0.65144D+03,0.59946D+03,0.54646D+03,0.49550D+03,0.44658D+03,
     & 0.40072D+03,0.35791D+03,0.31919D+03,0.28453D+03,0.25294D+03,
     & 0.22440D+03,0.19790D+03,0.17344D+03,0.15000D+03,0.12800D+03,
     & 0.10800D+03,0.90000D+02,0.73000D+02,0.57000D+02,0.43000D+02,
     & 0.31000D+02,0.20000D+02,0.10000D+02,0.56200D+01,0.31600D+01,
     & 0.17800D+01,0.10000D+01,0.56200D+00,0.31600D+00,0.17800D+00,
     & 0.10000D+00,0.50000D-01,0.20000D-01,0.00000D+00/)


      REAL*8, PARAMETER ::  DLOG2 =.30103d0,  ULMNH2=1.85124d0,
     *     ULMNCH=-.8160d0, ULMNN2=-1.527d0, ULMNF1=-4.780d0,
     *     ULMNO3=-1.368d0, ULMNCO= 1.523d0, ULMNF2=-4.524d0,
     *     USO2S=.042d0
      REAL*8, DIMENSION(NLCF,NRCF) :: XTU,XTD
      REAL*8, DIMENSION(NLCF,NWVCF,NRCF) :: DXUP,DXDN
      REAL*8 PRATCF(NLCF)
      INTEGER MLGAS(21)
      INTEGER I,IM,L,LL,LCF,LCFdn,LCFup,NLPrat,IULOW,IU,IPX,IAA,ITX
     *     ,IGAS,NG,KK,IK1,IK2,IPU,IK,NU,IUA,nsum
     *     ,IUB,IH2O0,IG  ,ICDlow,ICO20, IO3low,IO30, IUW,IU1,IU2
     *     ,i2u1,i2u2,i3u1,i3u2,i6u1,i6u2,i7u1,i7u2,i8u1,i8u2,i9u1,i9u2

      REAL*8 UH2O,UCO2L,UO3LL,UCH4L,UN2OL,UCF1L,UCF2L,USO2
     *                       ,UCH4L1,CH4RAT
     *     ,DUH2,DU1,DU2,DUCO,D2U1,D2U2,DUO3,D3U1,D3U2,DUCH,D7U1,D7U2
     *     ,DUN2,D6U1,D6U2,DUF1,D8U1,D8U2,DUF2,D9U1,D9U2,SUM1,SUM2,sumPR
     *     ,TAUT1,TAUT2,TAUHFB,TAUCF,TAUIPG,TAUSUM,TAU11,TAU12
     *     ,QAA,QAB,QBA,QBB, PLL,FPL,PU2, U,UP,UGAS, FNU1
     *     ,UAA,UAB,UBA,UBB, WPB, WTB,WTPU, XA,XB,XK,XUA,XUB
     *     ,WAA,WAB,WBA,WBB,WAAA,WAAB,WABA,WABB,WBAA,WBAB,WBBA,WBBB
      REAL*8 PRAT(LX),WT(LX)
      INTEGER LCFofL(LX)

      ! The variation of correction factors with the water vapor
      ! profile is determined via a two-step procedure.
      ! As for the other absorbers, a lookup-table dependence upon
      ! total column absorber amount is constructed by multiplying
      ! a reference vapor profile by a set of powers of 2 that index
      ! the tables.
      ! The actual shape of the water vapor profile for which fluxes
      ! are being computed is then folded into this column-oriented
      ! framework via per-layer interpolations for downward/upward flux
      ! correction factors that select the reference profile having
      ! the same column amount above/below each layer.
      ! In rare cases for which water vapor mixing ratios increase upward
      ! in the lower troposphere, an additional correction is performed.
      ! Since the interpolations in absorber amount are performed
      ! on the layers of the reference atmosphere, a call to REPART
      ! is needed to regrid the GCM water vapor to the reference layers.

      ! Interp. weights/indices and their prerequisites for downward-flux
      ! correction factors
      real*8, dimension(nlcf) :: qabove,uh2otl,duh2o1dn,duh2o2dn
      integer, dimension(nlcf) :: iuh2o1dn,iuh2o2dn
      REAL*8, PARAMETER ::  UCMRCF(NLCF) = (/
     & 0.10000D+01,0.11031D+01,0.12329D+01,0.14042D+01,0.16222D+01,
     & 0.19138D+01,0.23456D+01,0.29937D+01,0.40048D+01,0.55959D+01,
     & 0.81519D+01,0.12375D+02,0.19808D+02,0.32524D+02,0.54780D+02,
     & 0.93174D+02,0.15748D+03,0.25287D+03,0.37309D+03,0.50049D+03,
     & 0.60868D+03,0.71320D+03,0.84524D+03,0.10234D+04,0.12685D+04,
     & 0.16078D+04,0.20909D+04,0.28557D+04,0.41673D+04,0.66299D+04,
     & 0.12554D+05,0.26315D+05,0.59589D+05,0.10604D+06,0.18859D+06,
     & 0.33480D+06,0.59602D+06,0.10608D+07,0.18872D+07,0.33520D+07,
     & 0.59729D+07,0.11972D+08,0.30143D+08/)


      ! Interp. weights/indices and their prerequisites for upward-flux
      ! correction factors
      real*8, dimension(nlcf) :: qbelow,uh2oul,duh2o1up,duh2o2up
      integer, dimension(nlcf) :: iuh2o1up,iuh2o2up
      REAL*8, PARAMETER ::  UCMUCF(NLCF) = (/
     & 0.10701D+02,0.52946D+01,0.34738D+01,0.26072D+01,0.20943D+01,
     & 0.17432D+01,0.15016D+01,0.13328D+01,0.12176D+01,0.11398D+01,
     & 0.10879D+01,0.10532D+01,0.10317D+01,0.10186D+01,0.10108D+01,
     & 0.10064D+01,0.10040D+01,0.10027D+01,0.10020D+01,0.10016D+01,
     & 0.10014D+01,0.10012D+01,0.10010D+01,0.10008D+01,0.10006D+01,
     & 0.10005D+01,0.10004D+01,0.10002D+01,0.10002D+01,0.10001D+01,
     & 0.10000D+01,0.10000D+01,0.10000D+01,0.10000D+01,0.10000D+01,
     & 0.10000D+01,0.10000D+01,0.10000D+01,0.10000D+01,0.10000D+01,
     & 0.10000D+01,0.10000D+01,0.10000D+01/)

      real*8 :: dudp(lx),ddudp ! ddudp is vertical gradient of water vapor
      real*8 :: dxtru3_10(10) ! optional correction of top 10 layers <.2mb

#ifdef TAPER_UTCF
      real*8 :: pcen_tap,wt_one,pwid_tap
#endif

      ! compute QABOVE/QBELOW, the WV amount above/below each reference level
      CALL REPART (ULGAS(1,1),PLB,  NL+1,
     &             QABOVE,    PLBCF, NLCF+1)
      QBELOW = QABOVE
      DO L=2,NLCF
        QBELOW(L) = QBELOW(L) + QBELOW(L-1)
      ENDDO
      DO L=NLCF-1,1,-1
        QABOVE(L) = QABOVE(L) + QABOVE(L+1)
      ENDDO
      DO L=1,NLCF
        IF(QABOVE(L).GT.0.) THEN
          UH2OTL(L) = LOG10(UCMRCF(L)*QABOVE(L))
        ELSE
          UH2OTL(L) = 0. ! should not happen
        ENDIF
        IF(QBELOW(L).GT.0.) THEN
          UH2OUL(L) = LOG10(UCMUCF(L)*QBELOW(L))
        ELSE
          UH2OUL(L) = 0. ! below ground
        ENDIF
      ENDDO

C                          MLGAS DEF.
C                          ----------
C     H2O: 1,4,5   CO2: 2,6,7   O3: 3,8,9   N2O: 10,11,12   CH4: 13,14
C     SO2: 15,16   CFC: 17-20   WVCON: 21

      MLGAS(:)=1  !  1:21

C              KWVCON = ON/OFF flag for water vapor continuum absorption
C              ---------------------------------------------------------
      IF(KWVCON < 1) MLGAS(21)=0

C**** Find correction factors XTU and XTD
!     Prepare interpolation from PL to PCF pressure levels
      LCF=2
      NLPrat = NL
      DO L=L1,NL
        PLL = PL(L)
!       Find LCF s.t. PLLmid is between PCF(LCF) and PCF(LCF-1)
        DO WHILE (PLL < PCF(LCF))
          LCF = LCF + 1
          IF (LCF > NLCF) THEN      ! PL-levels higher than PCF_top
            NLPrat = L-1
            GO TO 100
          END IF
        END DO
        LCFofL(L) = LCF
        WT(L) = (PLL - PCF(LCF))/(PCF(LCF-1) - PCF(LCF))
        WT(L) = MIN(WT(L), 1D0)
        Prat(L) = DPL(L) / (DPCF(LCF-1)*WT(L) + DPCF(LCF)*(1-WT(L)))
      END DO
  100 CONTINUE

      ICDlow = 0 ! default: CO2 not low
      IO3low = 0 ! default: O3  not low
      IUlow  = 0 ! water vapor  not low

      UH2O = 1d-10 + SUM(ULGAS(L1:NL,1))
      IF (UH2O < 1.1d-10) THEN  ! low water vapor
        IUlow    = 1
        XTU(:,:) = XTU0(:,:)                                  ! 1:NLCF,1:NRCF
        XTD(:,:) = XTD0(:,:)                                  ! 1:NLCF,1:NRCF
        GO TO 180               ! if no water vapor
      END IF

      UCO2L = LOG10 (1d-10 + SUM(ULGAS(L1:NL,2)))
      UO3LL = LOG10 (1d-10 + SUM(ULGAS(L1:NL,3)))
      UCH4L = LOG10 (1d-10 + SUM(ULGAS(L1:NL,7)))
      UN2OL = LOG10 (1d-10 + SUM(ULGAS(L1:NL,6)))
      UCF1L = LOG10 (1d-10 + SUM(ULGAS(L1:NL,8)))
      UCF2L = LOG10 (1d-10 + SUM(ULGAS(L1:NL,9)))
      USO2  = SUM(ULGAS(L1:NL,13))

      CH4RAT = 1.
      if(UCH4L>.7) then                    ! high CH4 concentration case
         if(UCH4L<1.1) then
            UCH4L1 = 1.15*UCH4L - .1
         else if(UCH4L<1.7) then
            UCH4L1 = 0.70*UCH4L + .4
         else
            UCH4L1 = 0.375*UCH4L + .95
         end if
         CH4RAT = 10**UCH4L1/10**UCH4L
         UCH4L = UCH4L1
      end if

      if(UCO2L < -9.958607315d0) ICDlow = 1  ! if UCO2<1.1d-10 (low CO2)
      IF(UO3LL < -9.6)           IO3LOW = 1  ! low ozone

      DUCO=UCO2L-ULMNCO
      IF(DUCO.LT.0.) DUCO=0.
      I2U1=DUCO/DLOG2+1
      IF(I2U1.LT.1) I2U1=1
      IF(I2U1.GT.NUCF-1) I2U1=NUCF-1
      I2U2=I2U1+1
      D2U1=DUCO-(I2U1-1)*DLOG2
      D2U2=DLOG2-D2U1

      DUO3=UO3LL-ULMNO3
      IF(DUO3.LT.0.) DUO3=0.
      I3U1=DUO3/DLOG2+1
      IF(I3U1.LT.1) I3U1=1
      IF(I3U1.GT.NUCF-1) I3U1=NUCF-1
      I3U2=I3U1+1
      D3U1=DUO3-(I3U1-1)*DLOG2
      D3U2=DLOG2-D3U1

      DUCH=UCH4L-ULMNCH
      I7U1=DUCH/DLOG2+1
      IF(I7U1.LT.1) I7U1=1
      IF(I7U1.GT.NUCF-1) I7U1=NUCF-1
      I7U2=I7U1+1
      D7U1=DUCH-(I7U1-1)*DLOG2
      D7U2=DLOG2-D7U1

      DUN2=UN2OL-ULMNN2
      IF(DUN2.LT.0.) DUN2=DUN2*.5
      IF(DUN2.LT.-.56) DUN2=-.56
      I6U1=DUN2/DLOG2+1
      IF(I6U1.LT.1) I6U1=1
      IF(I6U1.GT.NUCF-1) I6U1=NUCF-1
      I6U2=I6U1+1
      D6U1=DUN2-(I6U1-1)*DLOG2
      D6U2=DLOG2-D6U1

      DUF1=UCF1L-ULMNF1
      IF(DUF1.LT.-.25) DUF1=-.25
      I8U1=DUF1/DLOG2+1
      IF(I8U1.LT.1) I8U1=1
      IF(I8U1.GT.NUCF-1) I8U1=NUCF-1
      I8U2=I8U1+1
      D8U1=DUF1-(I8U1-1)*DLOG2
      D8U2=DLOG2-D8U1

      DUF2=UCF2L-ULMNF2
      IF(DUF2.LT.-.2) DUF2=-.2
      I9U1=DUF2/DLOG2+1
      IF(I9U1.LT.1) I9U1=1
      IF(I9U1.GT.NUCF-1) I9U1=NUCF-1
      I9U2=I9U1+1
      D9U1=DUF2-(I9U1-1)*DLOG2
      D9U2=DLOG2-D9U1
c      IF(I9U1.GT.9xxxfixthis) THEN
c        I9U1=1
c        I9U2=2
c        D9U1=0.
c        D9U2=0.
c      end if

!     Find pressure ratios on PCF levels by averaging Prat
!     Fill missed layers copying from the nearest layer above
      LCFdn = 1                ! bottom of current segment
      LCFup = LCFofL(L1)       ! top of current segment
      sumPR = Prat(L1)
      PratCF(LCFdn:LCFup) = sumPR
      nsum = 1
      DO L=L1+1,NLPrat
        LCF=LCFofL(L)
        IF(LCF == LCFup) THEN  ! update the current PratLCF segment
          sumPR=sumPR+Prat(L)
          NSUM=NSUM+1
          PratCF(LCFdn:LCFup)=sumPR/DFLOAT(NSUM)
        ELSE                   ! start next PratLCF segment
          sumPR=Prat(L)
          PratCF(LCFup+1:LCF)=sumPR
          NSUM=1
          LCFdn=LCFup+1
          LCFup=LCF
        END IF
      END DO
      PratCF(LCFup+1:NLCF)=PratCF(LCFup) ! at top fill from below

      DO I=1,NLCF
        DUH2=UH2OUL(I)-ULMNH2
        IF(DUH2.LT.0.) DUH2=0.
        IU1=DUH2/DLOG2+1.
        IF(IU1.LT.1) IU1=1
        IF(IU1.GT.NWVCF-1) IU1=NWVCF-1
        IU2=IU1+1
        DUH2O1up(I)=DUH2-(IU1-1)*DLOG2
        DUH2O2up(I)=DLOG2-DUH2O1up(I)
        IUH2O1up(I)=IU1
        IUH2O2up(I)=IU2
      ENDDO

      DO I=1,NLCF
        DUH2=UH2OTL(I)-ULMNH2
        IF(DUH2.LT.0.) DUH2=0.
        IU1=DUH2/DLOG2+1.
        IF(IU1.LT.1) IU1=1
        IF(IU1.GT.NWVCF-1) IU1=NWVCF-1
        IU2=IU1+1
        DUH2O1dn(I)=DUH2-(IU1-1)*DLOG2
        DUH2O2dn(I)=DLOG2-DUH2O1dn(I)
        IUH2O1dn(I)=IU1
        IUH2O2dn(I)=IU2
      ENDDO

      DO IM=1,NRCF
      DO I=1,NLCF
      DO IUW=IUH2O1up(I),IUH2O2up(I)
      SUM1=(DXUP2(I,IUW,I2U2,IM)*D2U1+DXUP2(I,IUW,I2U1,IM)*D2U2)+
     $     (DXUP3(I,IUW,I3U2,IM)*D3U1+DXUP3(I,IUW,I3U1,IM)*D3U2)+

     $     PratCF(I)*(
     $     (DXUP7(I,IUW,I7U2,IM)*D7U1+DXUP7(I,IUW,I7U1,IM)*D7U2)+
     $     (DXUP6(I,IUW,I6U2,IM)*D6U1+DXUP6(I,IUW,I6U1,IM)*D6U2) )
     $    +(DXUP8(I,IUW,I8U2,IM)*D8U1+DXUP8(I,IUW,I8U1,IM)*D8U2)
     $    +(DXUP9(I,IUW,I9U2,IM)*D9U1+DXUP9(I,IUW,I9U1,IM)*D9U2)
      DXUP(I,IUW,IM)=SUM1/DLOG2+DXUP13(I,IUW,IM)*USO2/USO2S
      ENDDO
      DO IUW=IUH2O1dn(I),IUH2O2dn(I)
      SUM2=(DXDN2(I,IUW,I2U2,IM)*D2U1+DXDN2(I,IUW,I2U1,IM)*D2U2)+
     $     (DXDN3(I,IUW,I3U2,IM)*D3U1+DXDN3(I,IUW,I3U1,IM)*D3U2)+

     $     PratCF(I)*(
     $     (DXDN7(I,IUW,I7U2,IM)*D7U1+DXDN7(I,IUW,I7U1,IM)*D7U2)+
     $     (DXDN6(I,IUW,I6U2,IM)*D6U1+DXDN6(I,IUW,I6U1,IM)*D6U2) )
     $    +(DXDN8(I,IUW,I8U2,IM)*D8U1+DXDN8(I,IUW,I8U1,IM)*D8U2)
     $    +(DXDN9(I,IUW,I9U2,IM)*D9U1+DXDN9(I,IUW,I9U1,IM)*D9U2)
      DXDN(I,IUW,IM)=SUM2/DLOG2+DXDN13(I,IUW,IM)*USO2/USO2S
      ENDDO
      ENDDO ! LAYER
      ENDDO ! IM

      DO IM=1,NRCF
      DO I=1,NLCF
      DU1=DUH2O1up(I)
      DU2=DUH2O2up(I)
      IU1=IUH2O1up(I)
      IU2=IUH2O2up(I)
      XTU(I,IM)=((XTRUP(I,IU2,IM)+DXUP(I,IU2,IM))*DU1+
     $           (XTRUP(I,IU1,IM)+DXUP(I,IU1,IM))*DU2)/DLOG2
      DU1=DUH2O1dn(I)
      DU2=DUH2O2dn(I)
      IU1=IUH2O1dn(I)
      IU2=IUH2O2dn(I)
      XTD(I,IM)=((XTRDN(I,IU2,IM)+DXDN(I,IU2,IM))*DU1+
     $           (XTRDN(I,IU1,IM)+DXDN(I,IU1,IM))*DU2)/DLOG2
      ENDDO
      ENDDO

!**** Interpolate correction factors to model grid: XTU/D=>XTRU/D
  180 CONTINUE

      if(transmission_corrections) then
      ! note window region is position 1 in XTRU, XTRD
        XTRU(:,1)=1.
        XTRD(:,1)=1.

        DO L=L1,min(NLPrat,NL-1)
          LCF=LCFofL(L)
          XTRU(L,2:NRCF+1)=
     *    1.-PRAT(L)*(1.-XTU(LCF-1,:)*WT(L)-XTU(LCF,:)*(1.-WT(L)))
          XTRD(L,2:NRCF+1)=
     *    1.-PRAT(L)*(1.-XTD(LCF-1,:)*WT(L)-XTD(LCF,:)*(1.-WT(L)))
        END DO

        DO L=NLPrat+1,NL-1
          XTRU(L,2:NRCF+1)=XTU(NLCF,:)
          XTRD(L,2:NRCF+1)=XTD(NLCF,:)
        END DO

        XTRU(NL,2:NRCF+1) = 1.
        XTRD(NL,2:NRCF+1) = 1.

#ifdef TAPER_UTCF
        ! force upward transmission correction factors to 1 near the model top
        !pcen_tap = 1d0 ! center pressure (mb) of blending region
        pcen_tap = .1d0 ! center pressure (mb) of blending region
        pwid_tap = .5d0*pcen_tap ! width (mb) of blending region
        do l=nl,1,-1
          wt_one = .5d0*(1d0+tanh((pcen_tap-plb(l))/pwid_tap)) ! blending weight
          xtru(l,2:nrcf+1) = wt_one*1d0 + (1d0-wt_one)*xtru(l,2:nrcf+1)
          xtrd(l,2:nrcf+1) = wt_one*1d0 + (1d0-wt_one)*xtrd(l,2:nrcf+1)
          if(wt_one.lt.1d-3) exit ! far from model top
        enddo
#endif

      ! correction for cases when water vapor mixing ratio increases upward
        DO L=1,NL
          DUDP(L)=ULGAS(L,1)/(PLB(L)-PLB(L+1))
        ENDDO
        DO L=2,NL-1
          IF(PLB(L).LT.600.) EXIT
        !DDUDP=(DUDP(L)-DUDP(L+1))/(PLB(L)-PLB(L+1))
          DDUDP=(DUDP(L-1)-DUDP(L))/(PL(L-1)-PL(L))
          IF(DDUDP.GE.0.) CYCLE
          IF(DDUDP .GT. -.00037D0) THEN ! avoid nonzero effect for DDUDP==0
            XTRD(L,2) = XTRD(L,2) - 100d0*DDUDP
          ELSE
            XTRD(L,2) = XTRD(L,2) + (.035d0-5.25d0*DDUDP)
          ENDIF
        ENDDO

      else
        XTRU(:,:)=1.
        XTRD(:,:)=1.
      endif

C**** Find TRGXLK
  200 TRGXLK(L1:NL,1:33)=0.D0
      IPX=2
      DO 600 L=L1,NL
C         Locate model layer pressure between IPX and IPX-1
  280 CONTINUE
      WPB = (PL(L)-PX(IPX))/(PX(IPX-1)-PX(IPX))
      IF(WPB >= 0 .or. IPX >= NPX) GO TO 290
      IPX = IPX+1
      GO TO 280
C         Locate model layer temperature between ITX and ITX+1
  290 CONTINUE
      WTB = (TLM(L)-TLOX)/DTX + 1
      ITX = WTB  ;  IF(ITX < 1) ITX=1  ;  IF(ITX >= NTX) ITX=NTX-1
      WTB = WTB-ITX

      WBB = WPB*WTB
      WBA = WPB-WBB
      WAB = WTB-WBB
      WAA = 1 - (WBB+WBA+WAB)

      DO 500 IGAS=1,21
      IF(MLGAS(IGAS) < 1) GO TO 500
      KK   = IG1X(KGX(IGAS))
      NG   = NGX (KGX(IGAS))
      UGAS = ULGAS(L,IGASX(IGAS))
      IF(IGAS == 13.OR.IGAS == 14) UGAS = UGAS*CH4RAT
      IF(IGAS == 17.OR.IGAS == 18) UGAS = UGAS + ULGAS(L,11)

      IF(IGAS < 21) GO TO 375

C     IGAS = 21                   Apply water vapor continuum absorption
C     ---------                   --------------------------------------
C                  KCSELF = ON/FF flag for H2O self broadening continuum
C                  -----------------------------------------------------
      IF(KCSELF <= 0) GO TO 335

      DO 330 IK1=1,2
      if (IK1==1) then
          IK2=1 ;       U = UGAS*1.15d0  ! thermal K-domain 1
      else ! IK1=2
             IK2=33 ;   U = UGAS*XCSELF  ! thermal K-domain 2-33
      end if
      PU2 = PL(L)/DPL(L) * U**2
      IF (PU2 > PDPU2(1)) THEN
        IPU=2
        do while (PU2>PDPU2(IPU) .and. IPU<NPU2) ; IPU=IPU+1 ; end do
        WTPU   = (PU2-PDPU2(IPU-1))/(PDPU2(IPU)-PDPU2(IPU-1))
        DO IK=IK1,IK2
          TAUT1 = WTPU*(H2OCN8(IK,ITX,IPU)-H2OCN8(IK,ITX,IPU-1))+
     +            H2OCN8(IK,ITX,IPU-1)
          TAUT2 = WTPU*(H2OCN8(IK,ITX+1,IPU)-H2OCN8(IK,ITX+1,IPU-1))+
     +            H2OCN8(IK,ITX+1,IPU-1)
          TRGXLK(L,KK) = TRGXLK(L,KK) + (WTB*(TAUT2-TAUT1) + TAUT1)
          KK=KK+1
        END DO
      ELSE
        WTPU   = PU2/PDPU2(1)
        DO IK=IK1,IK2
          TAUT1 = WTPU*H2OCN8(IK,ITX,1)
          TAUT2 = WTPU*H2OCN8(IK,ITX+1,1)
          TRGXLK(L,KK) = TRGXLK(L,KK) + (WTB*(TAUT2-TAUT1) + TAUT1)
          KK=KK+1
        END DO
      END IF
  330 CONTINUE

C               KCFORN = ON/FF flag for H2O foreign broadening continuum
C               --------------------------------------------------------
  335 IF(KCFORN < 1) GO TO 500
      KK=IG1X(KGX(IGAS))
      DO 370 IK1=1,2
      if(IK1==1) then
         IK2=1      ;   U = UGAS*1.15d0
      else ! IK1=2
             IK2=33 ;   U = UGAS*XCFORN
      end if
      UP=PL(L)/P0*U
      IF(UP > PU(1)) THEN
        IPU=2
        DO WHILE (UP>PU(IPU) .and. IPU<NPU) ; IPU=IPU+1 ; END DO
        WTPU=(UP-PU(IPU-1))/(PU(IPU)-PU(IPU-1))
        DO IK=IK1,IK2
          TAUT1=WTPU*(H2OCF8(IK,ITX,IPU)-H2OCF8(IK,ITX,IPU-1))+
     +          H2OCF8(IK,ITX,IPU-1)
          TAUT2=WTPU*(H2OCF8(IK,ITX+1,IPU)-H2OCF8(IK,ITX+1,IPU-1))+
     +          H2OCF8(IK,ITX+1,IPU-1)
          TRGXLK(L,KK) = TRGXLK(L,KK) + (WTB*(TAUT2-TAUT1) + TAUT1)
          KK=KK+1
        END DO
      ELSE
        DO IK=IK1,IK2
          WTPU=UP/PU(1)
          TAUT1=WTPU*H2OCF8(IK,ITX,1)
          TAUT2=WTPU*H2OCF8(IK,ITX+1,1)
          TRGXLK(L,KK) = TRGXLK(L,KK) + (WTB*(TAUT2-TAUT1) + TAUT1)
          KK=KK+1
        END DO
      END IF
  370 CONTINUE
      GO TO 500

  375 IF(IGAS < 17) GO TO 385
C                               IGAS=17-20       Chloro Fluoro Carbons
C                               ----------       ---------------------
      DO IK=1,NG
        XA=WTB*(XKCFC(IK,ITX+1,IGAS)-XKCFC(IK,ITX,IGAS))+
     +     XKCFC(IK,ITX,IGAS)
        XB=WTB*(XKCFC(IK,ITX+1,IGAS)-XKCFC(IK,ITX,IGAS))+
     +     XKCFC(IK,ITX,IGAS)
        XK=WPB*(XA-XB)+XB
        TAUCF=XK*UGAS
        TRGXLK(L,KK)=TRGXLK(L,KK)+TAUCF
        KK=KK+1
      END DO
      GO TO 500

  385 CONTINUE               !  IGAS=1-16        H2O,CO2,O3,N2O,CH4,SO2
C                               ---------        ----------------------
      NU = NUX(IGAS)
      XUA = (UGAS-ULOX(IPX  ,IGAS)) / DUX(IPX  ,IGAS)
      XUB = (UGAS-ULOX(IPX-1,IGAS)) / DUX(IPX-1,IGAS)
C     IF(NU <= 1) then  ;  XUA = 0  ;  XUB = 0  ;  endif
      IUA = XUA
      IUB = XUB

      QAA = 1
      QAB = 1
      IF(XUA <= 0)  then
         XUA = 0
         IUA = 0
         QAA = UGAS /  ULOX(IPX,IGAS)
         QAB = UGAS / (ULOX(IPX,IGAS)+DUX(IPX,IGAS))
      endif
      IF(XUA >= NU-1)  then
         XUA = NU-1
         IUA = NU-2
         QAA = UGAS / (ULOX(IPX,IGAS)+DUX(IPX,IGAS)*(NU-2))
         QAB = UGAS / (ULOX(IPX,IGAS)+DUX(IPX,IGAS)*(NU-1))
      endif
      QBA = 1
      QBB = 1
      IF(XUB <= 0)  then
         XUB = 0
         IUB = 0
         QBA = UGAS /  ULOX(IPX-1,IGAS)
         QBB = UGAS / (ULOX(IPX-1,IGAS)+DUX(IPX-1,IGAS))
      endif
      IF(XUB >= NU-1)  then
         XUB = NU-1
         IUB = NU-2
         QBA = UGAS / (ULOX(IPX-1,IGAS)+DUX(IPX-1,IGAS)*(NU-2))
         QBB = UGAS / (ULOX(IPX-1,IGAS)+DUX(IPX-1,IGAS)*(NU-1))
      endif
      UAB = XUA-IUA
      UBB = XUB-IUB
      UAA = 1-UAB
      UBA = 1-UBB

      WAAA = WAA*UAA*QAA
      WAAB = WAA*UAB*QAB
      WABA = WAB*UAA*QAA
      WABB = WAB*UAB*QAB
      WBAA = WBA*UBA*QBA
      WBAB = WBA*UBB*QBB
      WBBA = WBB*UBA*QBA
      WBBB = WBB*UBB*QBB

      IH2O0=0
      IF( (IGAS==6.OR.IGAS==8.OR.IGAS==10.OR.IGAS==13.OR.IGAS==15)
     +   .and. IULOW==1 ) IH2O0=1

      ICO20=0
      IF( (IGAS==4.OR.IGAS==9.OR.IGAS==11) .and. ICDLOW==1 ) ICO20=1

      IO30=0
      IF( (IGAS==5.or.IGAS==7.or.IGAS==12.or.IGAS==14.or.IGAS==16)
     +   .and. IO3LOW==1 ) IO30=1

!!!   WARNING: If IH2O0+ICO20+IO30=2 accuracy is reduced
!!!   WARNING: If IH2O0+ICO20+IO30=3 result is unusable

      DO 430 IG=1,NG
      IF      (IH2O0 == 1) THEN
      TAUIPG = WAAA*TAUWV0(IG+IGUX(IGAS)+NG* IUA   ,ITX  ,IPX) ! low H2O
     +       + WAAB*TAUWV0(IG+IGUX(IGAS)+NG*(IUA+1),ITX  ,IPX)
     +       + WABA*TAUWV0(IG+IGUX(IGAS)+NG* IUA   ,ITX+1,IPX)
     +       + WABB*TAUWV0(IG+IGUX(IGAS)+NG*(IUA+1),ITX+1,IPX)
     +       + WBAA*TAUWV0(IG+IGUX(IGAS)+NG* IUB   ,ITX  ,IPX-1)
     +       + WBAB*TAUWV0(IG+IGUX(IGAS)+NG*(IUB+1),ITX  ,IPX-1)
     +       + WBBA*TAUWV0(IG+IGUX(IGAS)+NG* IUB   ,ITX+1,IPX-1)
     +       + WBBB*TAUWV0(IG+IGUX(IGAS)+NG*(IUB+1),ITX+1,IPX-1)
      ELSE IF (ICO20 == 1) THEN
      TAUIPG = WAAA*TAUCD0(IG+IGUX(IGAS)+NG* IUA   ,ITX  ,IPX) ! low CO2
     +       + WAAB*TAUCD0(IG+IGUX(IGAS)+NG*(IUA+1),ITX  ,IPX)
     +       + WABA*TAUCD0(IG+IGUX(IGAS)+NG* IUA   ,ITX+1,IPX)
     +       + WABB*TAUCD0(IG+IGUX(IGAS)+NG*(IUA+1),ITX+1,IPX)
     +       + WBAA*TAUCD0(IG+IGUX(IGAS)+NG* IUB   ,ITX  ,IPX-1)
     +       + WBAB*TAUCD0(IG+IGUX(IGAS)+NG*(IUB+1),ITX  ,IPX-1)
     +       + WBBA*TAUCD0(IG+IGUX(IGAS)+NG* IUB   ,ITX+1,IPX-1)
     +       + WBBB*TAUCD0(IG+IGUX(IGAS)+NG*(IUB+1),ITX+1,IPX-1)
      ELSE IF (IO30 == 1) THEN
      TAUIPG = WAAA*TAUO30(IG+IGUX(IGAS)+NG* IUA   ,ITX  ,IPX) ! low O3
     +       + WAAB*TAUO30(IG+IGUX(IGAS)+NG*(IUA+1),ITX  ,IPX)
     +       + WABA*TAUO30(IG+IGUX(IGAS)+NG* IUA   ,ITX+1,IPX)
     +       + WABB*TAUO30(IG+IGUX(IGAS)+NG*(IUA+1),ITX+1,IPX)
     +       + WBAA*TAUO30(IG+IGUX(IGAS)+NG* IUB   ,ITX  ,IPX-1)
     +       + WBAB*TAUO30(IG+IGUX(IGAS)+NG*(IUB+1),ITX  ,IPX-1)
     +       + WBBA*TAUO30(IG+IGUX(IGAS)+NG* IUB   ,ITX+1,IPX-1)
     +       + WBBB*TAUO30(IG+IGUX(IGAS)+NG*(IUB+1),ITX+1,IPX-1)
      ELSE   !! if H2O, CO2, O3 are present (I..0=0)
      TAUIPG = WAAA*TAUTBL(IG+IGUX(IGAS)+NG* IUA   ,ITX  ,IPX)
     +       + WAAB*TAUTBL(IG+IGUX(IGAS)+NG*(IUA+1),ITX  ,IPX)
     +       + WABA*TAUTBL(IG+IGUX(IGAS)+NG* IUA   ,ITX+1,IPX)
     +       + WABB*TAUTBL(IG+IGUX(IGAS)+NG*(IUA+1),ITX+1,IPX)
     +       + WBAA*TAUTBL(IG+IGUX(IGAS)+NG* IUB   ,ITX  ,IPX-1)
     +       + WBAB*TAUTBL(IG+IGUX(IGAS)+NG*(IUB+1),ITX  ,IPX-1)
     +       + WBBA*TAUTBL(IG+IGUX(IGAS)+NG* IUB   ,ITX+1,IPX-1)
     +       + WBBB*TAUTBL(IG+IGUX(IGAS)+NG*(IUB+1),ITX+1,IPX-1)
      ENDIF

      TAUSUM=TRGXLK(L,KK)+TAUIPG
      IF(TAUSUM > 0) TRGXLK(L,KK)=TAUSUM
      KK=KK+1
  430 CONTINUE
  500 CONTINUE


C-------------------------------------------------------------------
C            H2O WINDOW ABSORPTION (2013)
C-------------------------------------------------------------------
      IF(MLGAS(1) == 1) THEN
        XK=WTB*(XKH2OW(ITX+1)-XKH2OW(ITX))+XKH2OW(ITX)
        TRGXLK(L,1)=TRGXLK(L,1)+XK*ULGAS(L,1)
      END IF

C                               CFC11 and CFC12 Window Absorption (1997)
C                               ----------------------------------------

      IF(MLGAS(17) == 1.OR.MLGAS(18) == 1) THEN
        XK=WTB*(XKCFCW(ITX+1,1)-XKCFCW(ITX,1))+XKCFCW(ITX,1)
        TAU11=XK*(ULGAS(L,8)+ULGAS(L,11))
        TRGXLK(L,1)=TRGXLK(L,1)+TAU11
      ENDIF
      IF(MLGAS(19) == 1.OR.MLGAS(20) == 1) THEN
        XK=WTB*(XKCFCW(ITX+1,2)-XKCFCW(ITX,2))+XKCFCW(ITX,2)
        TAU12=XK*ULGAS(L,9)
        TRGXLK(L,1)=TRGXLK(L,1)+TAU12
      ENDIF
  600 CONTINUE

!     Optional LW up-flux correction for top 10 layers above 0.2 mb
      if (kfpco2==4) then
         call get_dxtru3_corr(dxtru3_10,jlat,mlat46,jday)
         xtru(nl-9:nl,3) = 1.d0+dxtru3_10(1:10)
      end if

      RETURN
      END SUBROUTINE TAUGAS

      SUBROUTINE THERML
#ifdef PLANET_PARAMS
      use constant, only : kapa ! exceptional use of external module
#endif
      IMPLICIT NONE
C     ------------------------------------------------------------------
C             Top-cloud Thermal Scattering Correction Control Parameters
C             ----------------------------------------------------------
C
C             ECLTRA = 1.0  Scattering correction is enabled
C             with KCLDEM = 1, Rigorous scattering correction is applied
C             with KCLDEM = 0, Approximate scattering correction is used
C
C             ECLTRA = 0.0  No scattering correction is used
C                                          (Independent of KCLDEM value)
C
C     ------------------------------------------------------------------
C                                   Lower Edge Temperature Interpolation
C                                   ------------------------------------
C     TLGRAD=1.0  (Default)
C                 Layer-mean temperatures (TLM) supplied by GCM are used
C                 to define the layer edge temperature TLT (top) and TLB
C                 (bottom) using overall atmospheric temperature profile
C                 to establish temperature gradient within each layer so
C                 as to minimize the temperature discontinuities between
C                 layer edges and to conserve layer thermal energy.
C
C     TLGRAD=0.0  This results in isothermal layers with TLT = TLB = TLM
C
C     TLGRAD<0.0  TLT and TLB are used as specified, without any further
C                 adjustments.  This is mainly for off-line use when the
C                 temperature profile (TLM,TLT,TLB) can be fully defined
C                 from a continuous temperature profile.
C
C     NOTE:       TLGRAD can also accommodate values between 0.0 and 1.0
C
C     PTLISO      (Default PTLISO=2.5mb)
C                 Pressure level above which model layers are defined to
C                 be isothermal.  This is appropriate for optically thin
C                 layers where emitted flux depends on mean temperature.
C     ------------------------------------------------------------------
      REAL*8 :: PX(9)=(/1001.,973.,934.,865.,752.,603.,439.,283.,156./)
      REAL*8 :: ALG2=.30103d0, TAUMNL=-2.20412d0

      REAL*8, PARAMETER :: R6=.16666667D0, R24=4.1666667D-02
      REAL*8, PARAMETER :: A=0.3825D0,B=0.5742D0,C=0.0433D0

#ifndef PLANET_PARAMS
      real*8, parameter :: kapa = .286d0
#endif

      REAL*8 TA,TB,TC,P1,P2,P3,P4,DT1CPT,DTHALF,CLTAUX,CLTAUS,CLCOSB
     *     ,CTX,DT2,DT1,CTG,DG2,DG1,WT1,WT2,WT3,WT4,WT5,WT6,WT7
     *     ,WT8,BG,DNACUM,DNBCUM,DNCCUM,TAUAG,TAUAP,TAUBP,TAUCP,TAUAX
     *     ,TAUBX,TAUCX,XTRDL,BTOP,BBOT,BBAR,TX,PLBN,F,TAUA,TAUB,TAUC
     *     ,BDIF,BBTA,BBTB,BBTC,TRANA,TRANB,TRANC,DEC
     *     ,DEB,DEA,COALB1,COALB2,COALB3,FDNABC,UNA,UNB,UNC,FUNABC
     *     ,PFW,DPF,CTP,DP1,DP2,TAUBG,TAUCG,DDFLUX,XTRUL
     *     ,FSUM,XFSUM,PLL,DTAU0,TAUPLG,AP1,AP2,XTF,XTFACN
      REAL*8 ENA(LX),ENB(LX),ENC(LX), TRA(LX),TRB(LX),TRC(LX)
      REAL*8 DNA(LX),DNB(LX),DNC(LX), WTLB(LX),WTLT(LX)
      REAL*8 RIJTCK(6,33), FDXTCK(3,33),FEMTCK(3,33),ALBTCK(3,33)
      REAL*8 CLPI0(33),CLPI0K
      INTEGER K,L,LL,II,ITL,ICT,IT1,IT2,IP1,IP2,ICG,IG1,IG2,IMOL
     *     ,IPF,ICP,ITLT(LX),ITLB(LX),IP,IPX0,ITAU1,ITAU2,LTOPA
     *     ,LCL(LX),ia,iaa,ic,iu,lvlo,lvhi,lskip,lcbot,nclds,icomb

C-----------------------------------------------------------------------
C                                   Layer edge temperature interpolation
C-----------------------------------------------------------------------
      if (TLGRAD < 0.D0) GO TO 130
      TA = TLM(L1)
      TB = TLM(L1+1)
      P1 = PLB(L1)
      P2 = PLB(L1+1)
      P3 = PLB(L1+2)
      DT1CPT = .5*TA*(P1**kapa-P2**kapa) / PL(L1)**kapa
      DTHALF = (TA-TB)*(P1-P2)/(P1-P3)
      if (DTHALF > DT1CPT) DTHALF = DT1CPT
      TLB(L1) = TA+DTHALF*TLGRAD
      TLT(L1) = TA-DTHALF*TLGRAD
      DO L = L1+1,NL-1
        TC = TLM(L+1)
        P4 = PLB(L+2)
        DTHALF = .5*((TA-TB)/(P1-P3)+(TB-TC)/(P2-P4))*(P2-P3)*TLGRAD
        TLB(L) = TB+DTHALF
        TLT(L) = TB-DTHALF
        TA = TB
        TB = TC
        P1 = P2
        P2 = P3
        P3 = P4
      END DO
      DTHALF = (TA-TB)*(P2-P3)/(P1-P3)*TLGRAD
      TLB(NL) = TC+DTHALF
      TLT(NL) = TC-DTHALF
      DO L = NL,L1,-1
        if (PLB(L) > PTLISO) GO TO 130
        TLT(L) = TLM(L)
        TLB(L) = TLM(L)
      END DO
  130 CONTINUE
      TLB(NL+1) = TLT(NL)

C     ------------------------------------------------------------------
C     weight assignments for Planck function interpolation
C     (Effective range (K) is from TK = planck_tmin to TK = planck_tmax)
C     ------------------------------------------------------------------

      DO 140 L=L1,NL
      ITLB(L) = TLB(L)
      WTLB(L) = TLB(L)-ITLB(L)
      if (ITLB(L) < planck_tmin  ) ITLB(L) = planck_tmin
      if (ITLB(L) > planck_tmax-1) ITLB(L) = planck_tmax-1
      ITLT(L) = TLT(L)
      WTLT(L) = TLT(L)-ITLT(L)
      if (ITLT(L) < planck_tmin  ) ITLT(L) = planck_tmin
      if (ITLT(L) > planck_tmax-1) ITLT(L) = planck_tmax-1
  140 CONTINUE

      if (LTOPCL==0) GO TO 180

      DO 170 K=1,33
      CLTAUX=TXCTPG(K)+TRGXLK(LTOPCL,K)+1d-10
      CLTAUS=TSCTPG(K)
      CLCOSB=TGCTPG(K)
      CLPI0K=CLTAUS*ECLTRA/CLTAUX
      CLPI0(K)=CLPI0K
      CTX=CLTAUX*10.D0
      if (CLTAUX >= 3.D0) then
        CTX=CLTAUX*2 + 24
        if (CTX > 47.999999D0) CTX=47.999999D0
      end if
      ICT=CTX
      DT2=CTX-ICT
      DT1=1.D0-DT2
      IT1=ICT+1
      IT2=ICT+2
      CTP=CLPI0K*20.D0
      ICP=CTP
      DP2=CTP-ICP
      DP1=1.D0-DP2
      IP1=ICP+1
      IP2=ICP+2
      CTG=CLCOSB*20.D0
      ICG=CTG
      DG2=CTG-ICG
      DG1=1.D0-DG2
      IG1=ICG+1
      IG2=ICG+2
      WT1=DT1*DP1*DG1
      WT2=DT2*DP1*DG1
      WT3=DT2*DP2*DG1
      WT4=DT1*DP2*DG1
      WT5=DT1*DP1*DG2
      WT6=DT2*DP1*DG2
      WT7=DT2*DP2*DG2
      WT8=DT1*DP2*DG2
      RIJTCK(:,K)=WT1*RIJTPG(:,IT1,IP1,IG1)+WT2*RIJTPG(:,IT2,IP1,IG1) ! 1:6
     +           +WT3*RIJTPG(:,IT2,IP2,IG1)+WT4*RIJTPG(:,IT1,IP2,IG1)
     +           +WT5*RIJTPG(:,IT1,IP1,IG2)+WT6*RIJTPG(:,IT2,IP1,IG2)
     +           +WT7*RIJTPG(:,IT2,IP2,IG2)+WT8*RIJTPG(:,IT1,IP2,IG2)
      FEMTCK(:,K)=WT1*FEMTPG(:,IT1,IP1,IG1)+WT2*FEMTPG(:,IT2,IP1,IG1) ! 1:3
     +           +WT3*FEMTPG(:,IT2,IP2,IG1)+WT4*FEMTPG(:,IT1,IP2,IG1)
     +           +WT5*FEMTPG(:,IT1,IP1,IG2)+WT6*FEMTPG(:,IT2,IP1,IG2)
     +           +WT7*FEMTPG(:,IT2,IP2,IG2)+WT8*FEMTPG(:,IT1,IP2,IG2)
      FDXTCK(:,K)=WT1*FDXTPG(:,IT1,IP1,IG1)+WT2*FDXTPG(:,IT2,IP1,IG1)
     +           +WT3*FDXTPG(:,IT2,IP2,IG1)+WT4*FDXTPG(:,IT1,IP2,IG1)
     +           +WT5*FDXTPG(:,IT1,IP1,IG2)+WT6*FDXTPG(:,IT2,IP1,IG2)
     +           +WT7*FDXTPG(:,IT2,IP2,IG2)+WT8*FDXTPG(:,IT1,IP2,IG2)
  170 CONTINUE

  180 CONTINUE
      TRDFLB(:)=0.D0
      TRUFLB(:)=0.D0

      BG=BGFEMT(1)
      TOTLZF(1:3)=0.D0
!sl   TRSLTS=0.D0
!sl   TRSLTG=0.D0
!sl   TRSLBS=0.D0

C     ------------------------------------------------------------------
C                                                      LOOP OVER K-BANDS
C     ------------------------------------------------------------------
      K=0
      IMOL=0
  200 CONTINUE
      K=K+1
      if (K > 33) GO TO 300
      BG=BGFEMT(K)
      if (K > 1 .and. K < 14) IMOL=1
      if (K > 13 .and. K < 26) IMOL=2
      if (K > 25) IMOL=3
      DFLB(NL+1,K)=0.D0
      DNACUM=0.D0
      DNBCUM=0.D0
      DNCCUM=0.D0
C**** Find top layer with absorbers: LtopA
      DO 210 L=NL,L1,-1
      LTOPA=L
      TAUAG=TRGXLK(L,K)
      TAUAP=TRCALK(L,K)+TRAALK(L,K)+TRBALK(L,K)+TRDALK(L,K)+TRVALK(L,K)
      TAUAX=TAUAG+TAUAP
      if (TAUAX > 1.D-06) GO TO 211
      DFLB(L,K)=0.D0
      ENA(L)=0.D0
      DNA(L)=0.D0
      TRA(L)=1.D0
      ENB(L)=0.D0
      DNB(L)=0.D0
      TRB(L)=1.D0
      ENC(L)=0.D0
      DNC(L)=0.D0
      TRC(L)=1.D0
  210 CONTINUE
      UFLB(L1:NL+1,K)=BG                 ! no absorbers in whole column
      TRUFLB(L1:NL+1)=TRUFLB(L1:NL+1)+BG
      TOTLZF(1)=TOTLZF(1)+BG
      TOTLZF(2)=TOTLZF(2)+BG
      TOTLZF(3)=TOTLZF(3)+BG
      GO TO 200                        ! next K

  211 CONTINUE
      FSUM=0.
      XFSUM=0.
      XTFACN=0.
      IPX0=9
C     ------------------------------------------------------------------
C                                              DOWNWARD FLUX COMPUTATION
C     ------------------------------------------------------------------
      DO 250 L=LTOPA,L1,-1
      BTOP = PLANCK(ITLT(L),K)-
     -      (PLANCK(ITLT(L),K)-PLANCK(ITLT(L)+1,K))*WTLT(L)
      BBOT = PLANCK(ITLB(L),K)-
     -      (PLANCK(ITLB(L),K)-PLANCK(ITLB(L)+1,K))*WTLB(L)
      TAUAG=TRGXLK(L,K)
      TAUAP=TRCALK(L,K)+TRAALK(L,K)+TRBALK(L,K)+TRDALK(L,K)+TRVALK(L,K)
      TAUAX=TAUAG+TAUAP
      IF(TAUAP.LT..003) GO TO 219
      PLL=PL(L)
      DO 217 IP=IPX0,1,-1
      IP1=IP
      IF(PLL.LT.PX(IP)) GO TO 218
  217 CONTINUE
  218 CONTINUE
      IF(IP1.EQ.9) IP1=8
      IP2=IP1+1
      IPX0=IP2
      TAUPLG=DLOG10(TAUAP)
      DTAU0=TAUPLG-TAUMNL
C     IF(DTAU0.LT.0.) DTAU0=0.
      ITAU1=DTAU0/ALG2+1
      IF(ITAU1.LT.1) ITAU1=1
      IF(ITAU1.GT.10) ITAU1=10
      ITAU2=ITAU1+1
      DT1=DTAU0-(ITAU1-1)*ALG2
      DT2=ALG2-DT1
      AP1=(XTFAC(ITAU2,IP1)*DT1+XTFAC(ITAU1,IP1)*DT2)/ALG2
      AP2=(XTFAC(ITAU2,IP2)*DT1+XTFAC(ITAU1,IP2)*DT2)/ALG2
      XTF=(AP2*(PLL-PX(IP1))+AP1*(PX(IP2)-PLL))/(PX(IP2)-PX(IP1))
      FSUM=FSUM+XTF/(1.+1.75*XFSUM**2)**2
      XTFACN=FSUM
      IF(XTFACN.GT.1.) XTFACN=1.
      IF(XTFACN.LT.0.) XTFACN=0.
      XFSUM=XFSUM+XTF
  219 CONTINUE

      XTRDL=XTRD(L,IMOL+1)
      XTRDL=XTRDL+XTFACN*(1.-XTRDL)

C               Optically thin limit emission/transmission approximation
C               --------------------------------------------------------

      IF (TAUAX >= 1.D-04) GO TO 220
      TAUBX=TAUAX+TAUAX
      TAUCX=10.D0*TAUAX
      BBAR=0.5D0*(BTOP+BBOT)
      TRA(L)=1.D0-TAUAX
      ENA(L)=BBAR*TAUAX
      DNA(L)=ENA(L)
      TX=TRA(L)*XTRDL ! ; if(TX > 1) TX=1
      DNACUM=DNACUM*TX+DNA(L)
      TRB(L)=1.D0-TAUBX
      ENB(L)=BBAR*TAUBX
      DNB(L)=ENB(L)
      TX=TRB(L)*XTRDL ! ; if(TX > 1) TX=1
      DNBCUM=DNBCUM*TX+DNB(L)
      TRC(L)=1.D0-TAUCX
      ENC(L)=BBAR*TAUCX
      DNC(L)=ENC(L)
      TX=TRC(L)*XTRDL ! ; if(TX > 1) TX=1
      DNCCUM=DNCCUM*TX+DNC(L)
      GO TO 230

C                     TAUB absorber-dependent extinction path adjustment
C                     --------------------------------------------------

  220 PLBN=PLB(L)
      ICOMB=0
      IF (TAUAG > TAUAP) THEN
        ICOMB=1
        TAUAG=TAUAX
      END IF
      TAUBG=TAUAG+TAUAG
      TAUCG=10.D0*TAUAG

      F=1
      if (IMOL==3 .and. PLBN>500 .and. TAUAG>.05d0 .and. TAUAG<.25) then
        F=23.71D0*TAUAG**2-7.113D0*TAUAG+1.296D0
        GO TO 221
      end if

      if (TAUAG > .1D0) then
        if      (IMOL==1) then
          if (PLBN > 250.D0) then
            F=.761D0
            if (TAUAG < 3.D0) F=.92D0-.053D0*TAUAG
            if (TAUAG < .2D0) F=1.091D0-.906D0*TAUAG
          else
            F=.718D0
            if (TAUAG < 2.5D0) F=.90D0-.073D0*TAUAG
            if (TAUAG < .2D0) F=1.115D0-1.146D0*TAUAG
          end if
        else if (IMOL==2) then
          if (PLBN > 250.D0) then
            F=.590D0
            if (TAUAG < 3.5D0) F=.93D0-.097D0*TAUAG
            if (TAUAG < .2D0) F=1.089D0-.894D0*TAUAG
          else
            F=.703D0
            if (TAUAG < 3.5D0) F=.92D0-.062D0*TAUAG
            if (TAUAG < .2D0) F=1.092D0-.924D0*TAUAG
          end if
        else if (IMOL==3) then
          if (PLBN > 250.D0) then
            F=.982D0
            if (TAUAG < .5D0) F=.99D0-.016D0*TAUAG
            if (TAUAG < .2D0) F=1.013D0-.132D0*TAUAG
          else
            F=.748D0
            if (TAUAG < 3.7D0) F=.97D0-.060D0*TAUAG
            if (TAUAG < .2D0) F=1.042D0-.420D0*TAUAG
          end if
        end if
      end if
  221 TAUBG=TAUBG*F

C                     TAUC absorber-dependent extinction path adjustment
C                     --------------------------------------------------
      F=1
      if (IMOL==3 .and. PLBN>500 .and. TAUAG>.01d0 .and. TAUAG<.25) then
        F=26.14D0*TAUAG**2-6.796D0*TAUAG+1.065D0
        GO TO 222
      end if

      if (TAUAG > .01D0) then
        if      (IMOL==1) then
          if (PLBN > 250.D0) then
            F=.712D0
            if (TAUAG < .37D0) F=.96D0-.67D0*TAUAG
            if (TAUAG < .02D0) F=1.053D0-5.34D0*TAUAG
          else
            F=.536D0
            if (TAUAG < .47D0) F=.87D0-.71D0*TAUAG
            if (TAUAG < .02D0) F=1.144D0-14.42D0*TAUAG
          end if
        else if (IMOL==2) then
          if (PLBN > 250.D0) then
            F=.710D0
            if (TAUAG < .75D0) F=.95D0-.32D0*TAUAG
            if (TAUAG < .02D0) F=1.056D0-5.64D0*TAUAG
          else
            F=.487D0
            if (TAUAG < .70D0) F=.90D0-.59D0*TAUAG
            if (TAUAG < .02D0) F=1.112D0-11.18D0*TAUAG
          end if
        else if (IMOL==3) then
          if (PLBN > 250.D0) then
            F=.961D0
            if (TAUAG < .5D0) F=.98D0-.039D0*TAUAG
            if (TAUAG < .02D0) F=1.021D0-2.08D0*TAUAG
          else
            F=.777D0
            if (TAUAG < .70D0) F=.98D0-.29D0*TAUAG
            if (TAUAG < .02D0) F=1.026D0-2.58D0*TAUAG
          end if
        end if
      end if
  222 TAUCG=TAUCG*F

      IF (ICOMB==0) THEN
        TAUBP=TAUAP+TAUAP
        TAUCP=10.D0*TAUAP
        TAUA=TAUAG+TAUAP
        TAUB=TAUBG+TAUBP
        TAUC=TAUCG+TAUCP
      ELSE
        TAUA=TAUAG
        TAUB=TAUBG
        TAUC=TAUCG
      END IF

      if (L==LTOPCL .and. KCLDEM==1) GO TO 225

      BDIF=BBOT-BTOP
      BBTA=BDIF/TAUA
      BBTB=BDIF/TAUB
      BBTC=BDIF/TAUC

C            Optically thick limit non-scattering emission approximation
C            -----------------------------------------------------------

      if (TAUA > 9.D0) then
        TRA(L)=0.D0
        TRB(L)=0.D0
        TRC(L)=0.D0
        ENA(L)=BTOP+BBTA
        ENB(L)=BTOP+BBTB
        ENC(L)=BTOP+BBTC
        DNA(L)=BBOT-BBTA
        DNB(L)=BBOT-BBTB
        DNC(L)=BBOT-BBTC
        DNACUM=BBOT-BBTA
        DNBCUM=BBOT-BBTB
        DNCCUM=BBOT-BBTC
        GO TO 230
      end if

      if (TAUA < 0.5D0) then
        TRANA = 1 - TAUA + (.5 - R6*TAUA + R24*(TAUA*TAUA))*(TAUA*TAUA)
      else
        TRANA = EXP(-TAUA)
      end if
      if (TAUB < 0.5D0) then
        TRANB = 1 - TAUB + (.5 - R6*TAUB + R24*(TAUB*TAUB))*(TAUB*TAUB)
      else
        TRANB = EXP(-TAUB)
      end if
      if (TAUC < 0.5D0) then
        TRANC = 1 - TAUC + (.5 - R6*TAUC + R24*(TAUC*TAUC))*(TAUC*TAUC)
      else
        TRANC = EXP(-TAUC)
      end if

      TRA(L)=TRANA
      ENA(L)=BTOP+BBTA-(BBOT+BBTA)*TRANA
      DNA(L)=BBOT-BBTA-(BTOP-BBTA)*TRANA
      TX=TRANA*XTRDL ! ; if(TX > 1) TX=1
      DNACUM=DNACUM*TX+DNA(L)
      TRB(L)=TRANB
      ENB(L)=BTOP+BBTB-(BBOT+BBTB)*TRANB
      DNB(L)=BBOT-BBTB-(BTOP-BBTB)*TRANB
      TX=TRANB*XTRDL ! ; if(TX > 1) TX=1
      DNBCUM=DNBCUM*TX+DNB(L)
      TRC(L)=TRANC
      ENC(L)=BTOP+BBTC-(BBOT+BBTC)*TRANC
      DNC(L)=BBOT-BBTC-(BTOP-BBTC)*TRANC
      TX=TRANC*XTRDL ! ; if(TX > 1) TX=1
      DNCCUM=DNCCUM*TX+DNC(L)
      GO TO 230

C                          ---------------------------------------------
C                          Top-cloud multiple scattering corrections for
C                          emitted, transmitted, and reflected radiances
C                          and fluxes at the top-cloud (L=LTOPCL) level.
C                          ---------------------------------------------

  225 CONTINUE
      IF (ICOMB==1) THEN
        TAUBP=TAUAP*(TAUBG/TAUAG)
        TAUCP=TAUAP*(TAUCG/TAUAG)
        TAUBG=TRGXLK(L,K)*(TAUBG/TAUAG)
        TAUCG=TRGXLK(L,K)*(TAUCG/TAUAG)
        TAUAG=TAUAG-TAUAP
      END IF
      TRA(L)=EXP(-TAUAG-TAUAP*FDXTCK(3,K))
      TRB(L)=EXP(-TAUBG-TAUBP*FDXTCK(2,K))
      TRC(L)=EXP(-TAUCG-TAUCP*FDXTCK(1,K))
      DEC=C*DNCCUM*RIJTCK(1,K)+B*DNBCUM*RIJTCK(2,K)+A*DNACUM*RIJTCK(3,K)
      DEB=C*DNCCUM*RIJTCK(2,K)+B*DNBCUM*RIJTCK(4,K)+A*DNACUM*RIJTCK(5,K)
      DEA=C*DNCCUM*RIJTCK(3,K)+B*DNBCUM*RIJTCK(5,K)+A*DNACUM*RIJTCK(6,K)
      ALBTCK(1,K)=C*RIJTCK(1,K)+B*RIJTCK(2,K)+A*RIJTCK(3,K)
      ALBTCK(2,K)=C*RIJTCK(2,K)+B*RIJTCK(4,K)+A*RIJTCK(5,K)
      ALBTCK(3,K)=C*RIJTCK(3,K)+B*RIJTCK(5,K)+A*RIJTCK(6,K)
      COALB1=1.D0-ALBTCK(1,K)
      COALB2=1.D0-ALBTCK(2,K)
      COALB3=1.D0-ALBTCK(3,K)
      TAUA=TAUAG+TAUAP*FEMTCK(3,K)
      TAUB=TAUBG+TAUBP*FEMTCK(2,K)
      TAUC=TAUCG+TAUCP*FEMTCK(1,K)
      TRANA=EXP(-TAUA)
      TRANB=EXP(-TAUB)
      TRANC=EXP(-TAUC)
      BDIF=BBOT-BTOP
      BBTA=BDIF/TAUA
      BBTB=BDIF/TAUB
      BBTC=BDIF/TAUC
      ENA(L)=(BTOP+BBTA-(BBOT+BBTA)*TRANA)*COALB3
      DNA(L)=(BBOT-BBTA-(BTOP-BBTA)*TRANA)*COALB3
      TX=TRA(L)*XTRDL ! ; if(TX > 1) TX=1
      DNACUM=DNACUM*TX+DNA(L)
      ENB(L)=(BTOP+BBTB-(BBOT+BBTB)*TRANB)*COALB2
      DNB(L)=(BBOT-BBTB-(BTOP-BBTB)*TRANB)*COALB2
      TX=TRB(L)*XTRDL ! ; if(TX > 1) TX=1
      DNBCUM=DNBCUM*TX+DNB(L)
      ENC(L)=(BTOP+BBTC-(BBOT+BBTC)*TRANC)*COALB1
      DNC(L)=(BBOT-BBTC-(BTOP-BBTC)*TRANC)*COALB1
      TX=TRC(L)*XTRDL ! ; if(TX > 1) TX=1
      DNCCUM=DNCCUM*TX+DNC(L)
      ENC(L)=ENC(L)+DEC
      ENB(L)=ENB(L)+DEB
      ENA(L)=ENA(L)+DEA
  230 CONTINUE
      FDNABC=A*DNACUM+B*DNBCUM+C*DNCCUM
      TRDFLB(L)=TRDFLB(L)+FDNABC
      DFLB(L,K)=FDNABC
  250 CONTINUE

C             Old form of scattering correction is skipped when KCLDEM=1
C             ----------------------------------------------------------

      if (KCLDEM==0 .and. LTOPCL > 0) then
        ENA(LTOPCL)=ENA(LTOPCL)*(1-TRCTCA(K))+TRCTCA(K)*DFLB(LTOPCL+1,K)
        ENB(LTOPCL)=ENB(LTOPCL)*(1-TRCTCA(K))+TRCTCA(K)*DFLB(LTOPCL+1,K)
        ENC(LTOPCL)=ENC(LTOPCL)*(1-TRCTCA(K))+TRCTCA(K)*DFLB(LTOPCL+1,K)
      end if

!sl   ------------------------------------------------------------------
!sl                                       SURFACE LAYER FLUX COMPUTATION
!sl   with TAUSL,FTAUSL=0 defaults, surface layer calculation is skipped
!sl   ------------------------------------------------------------------

      DFSL(K)=FDNABC
!sl   TAUA=TAUSL(K)+FTAUSL(K)
!sl   if (TAUA > 1.D-06) GO TO 24
      BG=BG+FDNABC*TRGALB(K)
      UNA=BG
      UNB=BG
      UNC=BG
      FUNABC=BG
!sl   GO TO 245
!sl24 CONTINUE
!sl   ITS=TSL
!sl   WTS=TSL-ITS
!sl   WTS1=1-WTS
!sl   BS = PLANCK(ITS,K)*WTS1 + PLANCK(ITS+1,K)*WTS
!sl   TA=EXP(-TAUA)
!sl   TB=TA*TA
!sl   TC=(TB*TB*TA)**2
!sl   DNA(1)=(DNA(1)-BS)*TA+BS
!sl   DNB(1)=(DNB(1)-BS)*TB+BS
!sl   DNC(1)=(DNC(1)-BS)*TC+BS
!sl   FDNABC=A*DNA(1)+B*DNB(1)+C*DNC(1)
!sl   BG=BGFEMT(K)+FDNABC*TRGALB(K)
!sl   UNA=(BG-BS)*TA+BS
!sl   UNB=(BG-BS)*TB+BS
!sl   UNC=(BG-BS)*TC+BS
!sl   FUNABC=A*UNA+B*UNB+C*UNC
!sl   BSP = PLANCK(ITS+1,K)*WTS1 + PLANCK(ITS+2,K)*WTS
!sl   BSM = PLANCK(ITS-1,K)*WTS1 + PLANCK(ITS  ,K)*WTS
!sl   SLABS=1.D0-A*TA-B*TB-C*TC
!sl   TRSLTS=TRSLTS+(BSP-BSM)*SLABS
!sl   TRSLTG=TRSLTG+BGFEMD(K)*SLABS
!sl   TRSLBS=TRSLBS+BS*SLABS

C     ------------------------------------------------------------------
C                                                UPWARD FLUX COMPUTATION
C     ------------------------------------------------------------------

  245 DO 260 L=L1,NL
      TRUFLB(L)=TRUFLB(L)+FUNABC
      UFLB(L,K)=FUNABC

C       ----------------------------------------------------------------
C       At top-cloud level, find  component of  upwelling flux reflected
C       downward by cloud bottom and add to downwelling flux below cloud
C       ----------------------------------------------------------------

      if (L==LTOPCL .and. KCLDEM==1) then
        DEC=C*UNC*RIJTCK(1,K)+B*UNB*RIJTCK(2,K)+A*UNA*RIJTCK(3,K)
        DEB=C*UNC*RIJTCK(2,K)+B*UNB*RIJTCK(4,K)+A*UNA*RIJTCK(5,K)
        DEA=C*UNC*RIJTCK(3,K)+B*UNB*RIJTCK(5,K)+A*UNA*RIJTCK(6,K)
        DO LL=L,L1,-1
          DNA(LL)=DNA(LL)+DEA
          DNB(LL)=DNB(LL)+DEB
          DNC(LL)=DNC(LL)+DEC
          DDFLUX=A*DEA+B*DEB+C*DEC
          TRDFLB(LL)=TRDFLB(LL)+DDFLUX
          DFLB(LL,K)=DFLB(LL,K)+DDFLUX
          if (LL == L1) exit ! LL-loop
          DEA=DEA*TRA(LL-1)
          DEB=DEB*TRB(LL-1)
          DEC=DEC*TRC(LL-1)
        END DO
      end if
      XTRUL=XTRU(L,IMOL+1)
      TX=TRA(L)*XTRUL ! ; if(TX > 1) TX=1
      UNA=UNA*TX+ENA(L)
      TX=TRB(L)*XTRUL ! ; if(TX > 1) TX=1
      UNB=UNB*TX+ENB(L)
      TX=TRC(L)*XTRUL ! ; if(TX > 1) TX=1
      UNC=UNC*TX+ENC(L)
      FUNABC=A*UNA+B*UNB+C*UNC
  260 CONTINUE

      if (K==1) then
        TRUFTW=FUNABC
        TRDFGW=TRDFLB(1)
        TRUFGW=BG
        WINDZF(1)=UNA
        WINDZF(2)=UNB
        WINDZF(3)=UNC
      end if

      TRUFLB(NL+1)=TRUFLB(NL+1)+FUNABC
      UFLB(NL+1,K)=FUNABC
      UFSL(K)=UFLB(1,K)
      TOTLZF(1)=TOTLZF(1)+UNA
      TOTLZF(2)=TOTLZF(2)+UNB
      TOTLZF(3)=TOTLZF(3)+UNC

      GO TO 200  !  next K
  300 CONTINUE

      TRNFLB(L1:NL+1) = TRUFLB(L1:NL+1) - TRDFLB(L1:NL+1)
      TRFCRL(L1:NL)   = TRNFLB(L1+1:NL+1) - TRNFLB(L1:NL)

C**** Window region and spectr. integrated total flux diagnostics
      DO 390 II=0,3
      if (II > 0) then
        PFW = TOTLZF(II) ; IF (PFW < 1) PFW=1
        if (PFW > 899.999d0) PFW=899.999d0
        IPF=PFW
        TOTLZT(II) = TKPFT(IPF) + (PFW-IPF)*(TKPFT(IPF+1)-TKPFT(IPF))

        PFW = 10*WINDZF(II)
      else
        PFW = 10*TRUFTW
      end if
        IF (PFW < 1.0001d-2) PFW=1.0001d-2
        IF (PFW > 719.999d0) PFW=719.999d0
        IPF=PFW
        IF (PFW < 1) THEN
          PFW = 100.*PFW
          IPF = PFW ; DPF = PFW-IPF         ! IPF=  1- 99
        ELSE IF (PFW < 10) THEN
          PFW = 10.*PFW
          IPF = PFW ; DPF = PFW-IPF
          IPF = IPF + 90                    ! IPF=100-189
        ELSE
          IPF = PFW ; DPF = PFW-IPF
          IPF = IPF + 180                   ! IPF=190-899
        END IF
      if (II > 0) then
        WINDZT(II) = TKPFW(IPF) + DPF*(TKPFW(IPF+1)-TKPFW(IPF))
      else
        BTEMPW     = TKPFW(IPF) + DPF*(TKPFW(IPF+1)-TKPFW(IPF))
      end if
  390 CONTINUE

      RETURN
      END SUBROUTINE THERML

      SUBROUTINE SOLAR0
      IMPLICIT NONE

      INTEGER, PARAMETER, DIMENSION(17) :: NMKWAV = (/ 200, 360, 770,
     *     795, 805, 810, 860,1250,1500,1740,2200,3000,3400,3600,3800
     *     ,4000,9999/)
      INTEGER, PARAMETER, DIMENSION(16) :: LORDER = (/15,14, 8, 7, 6, 5,
     *     13, 12, 4, 3, 2, 1, 11, 10, 9,16/)
      INTEGER N2,I

      DO I=1,30
        DBLN(I) = 2**I
      END DO

      NORDER(1:16)=LORDER(1:16)
      NMWAVA(1:16)=NMKWAV(1:16)
      NMWAVB(1:16)=NMKWAV(2:17)

      TCLMIN=MIN(TAUIC0,TAUWC0)

      call SETO2A

      RETURN
      END SUBROUTINE SOLAR0

      SUBROUTINE SOLARM
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SOLARM Returns:
C                      SRDFLB   Solar downward flux at layer bottom edge
C                      SRUFLB   Solar  upward  flux at layer bottom edge
C                      SRNFLB   Solar net downward flux  (in Watts/m**2)
C                      SRFHRL   Solar heating rate/layer (in Watts/m**2)
C                      FSRNFG   Solar flux abs at ground by surface-type
C                               (see explanatory note at end of SOLARM)
C               Also:
C                      TOA:   SRIVIS SROVIS PLAVIS  SRINIR SRONIR PLANIR
C                      BOA:   SRDVIS SRUVIS ALBVIS  SRDNIR SRUNIR ALBNIR
C                      ATM:   SRTVIS SRRVIS SRAVIS  SRTNIR SRRNIR SRANIR
C                             SRXVIS SRXNIR (Direct beam only at ground)
C
C           Spectral:  (by k-distribution/pseudo-spectral) breakdown:
C                      SKDFLB   Solar downward flux at layer bottom edge
C                      SKUFLB   Solar  upward  flux at layer bottom edge
C                      SKNFLB   Solar net downward flux  (in Watts/m**2)
C                      SKFHRL   Solar heating rate/layer (in Watts/m**2)
C
C                      SRKALB   Planetary albedo (by spectral breakdown)
C                      SRKINC   Incident fluxedo (by spectral breakdown)
C                      SRKGAX   Direct  k-d flux absorbed by ground-type
C                      SRKGAD   Diffuse k-d flux absorbed by ground-type
C     ------------------------------------------------------------------
C     Remarks:
C              NORMS0=1 Incident (TOA) Solar flux normalized to equal S0
C                       (COSZ dependence included in calculated results)
C                        The returned solar fluxes have to be multiplied
C                       by COSZ to yield actual atmospheric heating rate
C
C              NMKWAV   Spectral/k-distribution subdivisions are nominal
C                       (due to spectral trading of absorption features)
C
C                 VIS   Designates solar visible wavelengths (   <770nm)
C                 NIR   Designates solar near-IR wavelengths (770>   nm)
C                       VIS comprises .53 of S0, NIR comprises .47 of S0
C     ------------------------------------------------------------------
C
C     ------------------------------------------------------------------
C         Fractional solar flux k-distribution/pseudo-spectral intervals
C
C         KSLAM=    1     1     2     2     5     5     5     5
C             K=    1     2     3     4     5     6     7     8
C     DATA DKS0/ .010, .030, .040, .040, .040, .002, .004, .013,
C         KSLAM=    1     1     1     3     4     6     6     1
C             K=    9    10    11    12    13    14    15    16
C    +           .002, .003, .003, .072, .200, .480, .050, .011/
C
C     ------------------------------------------------------------------
C     The nominal spectral order for k-dist/pseudo-spectral intervals is
C     (WavA and WavB designate approximate spectral interval boundaries)
C
C             L=   12    11    10     9     6     5     4     3
C     WavA (nm)= 3000  2200  1740  1500   810   805   795   770
C     WavB (nm)= 3400  3000  2200  1740   860   810   805   795
C             K=    1     2     3     4     5     6     7     8
C     DATA DKS0/ .010, .030, .040, .040, .040, .002, .004, .013,
C
C             L=   15    14    13     8     7     2     1    16
C     WavA (nm)= 3800  3500  3400  1250   860   360   200  4000
C     WavB (nm)= 4000  3800  3600  1500  1250   770   360  9999
C             K=    9    10    11    12    13    14    15    16
C    +           .002, .003, .003, .072, .200, .480, .050, .011/
C
C     ------------------------------------------------------------------
C     6 spectral intervals overlap the 16 solar k-distribution intervals
C
C     Cloud and aerosol Mie scattering parameters (also surface albedos)
C     are averaged over these spectral intervals. These intervals are in
C     reverse spectral order. Thus spectral interval 6 refers to visible
C     (VIS) wavelengths, intervals 1-5 refer to nearIR (NIR) wavelengths
C     KSLAM designates the spectral interval of first 14 k-distributions
C     (K=15 for UV ozone absorption refers to (VIS) spectral interval 6)
C     (K=16 represents strong absorbing spectral regions via interval 1)
C
C     The nominal Mie scattering spectral band subdivisions are:
C
C                         -------------NIR------------      VIS
C                     L=    1     2     3     4     5        6
C             WavA (nm)=  2200  1500  1250   860   770      300
C             WavB (nm)=  4000  2200  1500  1250   860      770
C
C     ------------------------------------------------------------------

      REAL*8 COLEXT(6),COLSCT(6),COLGCB(6)   ! ,ALLGCB(6)

C                            -------------------------------------------
C                            NO2, O3 Chappuis Band, Rayleigh, parameters
C                            -------------------------------------------
      REAL*8, PARAMETER :: XCMNO2=5.465d0, XCMO3=.0399623d0
      REAL*8, PARAMETER ::
     &    SIGMA_RAY= 4.4028450689125004d-07
C           Rayleigh scattering cross-section [m2/mol]
      REAL*8 RNB(LX),RNX(LX), TNB(LX),TNX(LX), XNB(LX),XNX(LX)
      REAL*8 SRB(LX),SRX(LX), VRU(LX+1),VRD(LX+1),FAC(LX+1)
      REAL*8 AO3D(LX),AO3U(LX),AO3X(LX)
      REAL*8 S0COSZ,COSMAG,SECZ,TAURAY,RTAU,SUMEXT,SUMCST,SUMCGB,COLPFG
     *     ,SURFBB,TAUSBB,ALLTAU,TAULAY,GCBLAY,RTAUL,DKS0X,RBNB,RBNX
     *     ,RCNB,RCNX,TLN,PLN,ULN,TERMA,TERMB,TAU1,TAU,PIZERO,PR,PT,DBLS
     *     ,XANB,XANX,TANB,TANX,XXT,RASB,RASX,BNORM,XNORM,RARB,RARX,XATB
     *     ,DENOM,DB,DX,UB,UX,RBXTOA,ATOPX,ATOPD,O3CMX,O3CMD
     *     ,SUMSCT,SUMGCB,XXG,SURX,PFF,XATC,XBNB,XBNX,TBNB,TBNX,XBTB
     *     ,ABOTX,ABOTD,AO3UXN,AO3UDN,SRKA16,DKS0XX,TRNC,CLX,TRNU,TRN1
     *     ,TRN2,TRN3,TAUG,TAU2,TAU3,S0VIS,S0NIR,SUMX,SUMD,SUMU,SUMN
     *     ,SUMH,SGPG
      INTEGER I,K,KK,L,M,N,NN,KLAM,NDBLS
      REAL*8 :: WVCOL,ZWPATH,ALPH,BETA,FACK12,ROOT,PTROOT,TAUK,FACK13

      S0COSZ=S0 ; IF (NORMS0==0) S0COSZ=S0*COSZ

      SRDFLB(L1:NL+1)=0    ; SRUFLB(L1:NL+1)=0  ;  SRNFLB(L1:NL+1)=0
      SRFHRL(L1:NL )=0

      SKDFLB(L1:NL+1,16)=0 ; SKUFLB(L1:NL+1,16)=0

      SRKALB(1:16)=0.D0                             ! for WRITER only
      dblext=0. ; dblsct=0. ; dblgcb=0. ; dblpi0=0. ! for writer only
      skdflb=0. ; sknflb=0. ; skuflb=0.             ! for writer only
      skfhrl=0. ; srkgax=0. ; srkgad=0.             ! for writer only

C                     TOA solar flux VIS/NIR subdivision
C                     (incident, outgoing, plane albedo)
C                     ----------------------------------
      SRIVIS=0.D0
      SROVIS=0.D0
      PLAVIS=1.D0
      SRINIR=0.D0
      SRONIR=0.D0
      PLANIR=1.D0
C                     BOA solar flux VIS/NIR subdivision
C                     (incident, upward, surface albedo)
C                     ----------------------------------
      SRDVIS=0.D0
      SRUVIS=0.D0
      ALBVIS=1.D0
      SRDNIR=0.D0
      SRUNIR=0.D0
      ALBNIR=1.D0
C                     Fractional atmos only flux VIS/NIR subdivision
C                     (fractions reflected, transmitted, & absorbed)
C                     ----------------------------------------------
      SRRVIS=1.D0
      SRTVIS=0.D0
      SRAVIS=0.D0
      SRRNIR=0.D0
      SRTNIR=0.D0
      SRANIR=0.D0
C                     Direct beam, fractional S0 VIS/NIR subdivision
C                     ----------------------------------------------
      SRXVIS=0.D0
      SRXNIR=0.D0
C                     Ground surface absorbed solar flux subdivision
C                     according to 4 fractional surface-type albedos
C                     ----------------------------------------------
      FSRNFG(1:4)=0

      IF(COSZ < 0.001D0) RETURN
      COSMAG=35.D0/SQRT(1224.D0*COSZ*COSZ+1.D0)
      SECZ=1.D0/COSZ

C     Compute Rayleigh optical depth, still missing dP in units of mbar
      TAURAY=SIGMA_RAY/(grav*mair*1d-3)*1d+2*FRAYLE

      DO 90 K=1,6
      RTAU=1.D-10
      IF(K==6) RTAU=TAURAY
      COLEXT(K)=0.D0
      COLSCT(K)=0.D0
      COLGCB(K)=0.D0
      DO 30 L=L1,NL
      RTAUL=RTAU*(PLB(L)-PLB(L+1))
      SUMEXT=RTAUL+SRCEXT(L,K)+SRAEXT(L,K)+SRBEXT(L,K)
     +                        +SRDEXT(L,K)+SRVEXT(L,K)
      SUMSCT=RTAUL+SRCSCT(L,K)+SRASCT(L,K)+SRBSCT(L,K)
     +                        +SRDSCT(L,K)+SRVSCT(L,K)
      SUMGCB=      SRCSCT(L,K)*SRCGCB(L,K)+SRASCT(L,K)*SRAGCB(L,K)
     +            +SRBSCT(L,K)*SRBGCB(L,K)+SRDSCT(L,K)*SRDGCB(L,K)
     +            +SRVSCT(L,K)*SRVGCB(L,K)
      DBLEXT(L,K)=SUMEXT
      DBLSCT(L,K)=SUMSCT
      DBLGCB(L,K)=SUMGCB/(SUMSCT+1.D-10)
      DBLPI0(L,K)=SUMSCT/(SUMEXT+1.D-10)
      COLEXT(K)=COLEXT(K)+DBLEXT(L,K)
      COLSCT(K)=COLSCT(K)+DBLSCT(L,K)
      COLGCB(K)=COLGCB(K)+DBLSCT(L,K)*DBLGCB(L,K)
   30 CONTINUE
      COLGCB(K)=COLGCB(K)/(COLSCT(K)+1.D-10)

      IF(KANORM > 0) THEN
C     -----------------------------------------------------------------
C     KANORM (default = 0)  Option to renormalize aerosol column albedo
C                           to make column albedo less dependent on the
C                           number of model layers due to SGP treatment
C
C                           KANORM=1  aerosol column only is normalized
C
C                           KANORM=2  aerosol plus ground is normalized
C                                     with Tau equivalent ground albedo
C                                     ---------------------------------
      COLPFG=COLGCB(K)
      SURFBB=SRBALB(K)
      TAUSBB=0.D0
      IF(KANORM > 1) CALL GTSALB(XXG,XXT,SURX,SURFBB,COLPFG,TAUSBB,2)
      DBLEXT(NL+1,K)=TAUSBB
      ALLTAU=TAUSBB+COLEXT(K)
      CALL SGPGXG(COSZ,ALLTAU,COLPFG,SGPG)
cc    ALLGCB(K)=SGPG
      DBLGCB(L1:NL,K)=SGPG
      ELSE

      DO 50 L=L1,NL
      TAULAY=DBLEXT(L,K)
      GCBLAY=DBLGCB(L,K)
      CALL SGPGXG(COSZ,TAULAY,GCBLAY,SGPG)
      DBLGCB(L,K)=SGPG
   50 CONTINUE
      ENDIF

      IF(LTOPCL == 0) GO TO 90
      RTAU=1.D-10
      IF(K==6) RTAU=TAURAY
      COLEXT(K)=0.D0
      COLSCT(K)=0.D0
      COLGCB(K)=0.D0
      DO 60 L=L1,NL
      IF(SRCEXT(L,K) < TCLMIN) GO TO 60
      RTAUL=RTAU*(PLB(L)-PLB(L+1))
      SUMEXT=RTAUL+SRCEXT(L,K)+SRAEXT(L,K)+SRBEXT(L,K)
     +                        +SRDEXT(L,K)+SRVEXT(L,K)
      SUMSCT=RTAUL+SRCSCT(L,K)+SRASCT(L,K)+SRBSCT(L,K)
     +                        +SRDSCT(L,K)+SRVSCT(L,K)
      SUMGCB=      SRCSCT(L,K)*SRCGCB(L,K)+SRASCT(L,K)*SRAGCB(L,K)
     +            +SRBSCT(L,K)*SRBGCB(L,K)+SRDSCT(L,K)*SRDGCB(L,K)
     +            +SRVSCT(L,K)*SRVGCB(L,K)
      DBLEXT(L,K)=SUMEXT
      DBLSCT(L,K)=SUMSCT
      DBLGCB(L,K)=SUMGCB/(SUMSCT+1.D-10)
      DBLPI0(L,K)=SUMSCT/(SUMEXT+1.D-10)
      COLEXT(K)=COLEXT(K)+DBLEXT(L,K)
      COLSCT(K)=COLSCT(K)+DBLSCT(L,K)
      COLGCB(K)=COLGCB(K)+DBLSCT(L,K)*DBLGCB(L,K)
   60 CONTINUE
      COLGCB(K)=COLGCB(K)/(COLSCT(K)+1.D-10)

C     -----------------------------------------------------------------
C     KCNORM (default = 0)  Option to renormalize  cloud  column albedo
C                           to make column albedo less dependent on the
C                           number of model layers due to SGP treatment
C
C                           KCNORM=1    cloud column only is normalized
C
C                           KCNORM=2    cloud plus ground is normalized
C                                     with Tau equivalent ground albedo
C                                     ---------------------------------
      IF(KCNORM > 0) THEN
        COLPFG=COLGCB(K)
        SURFBB=SRBALB(K)
        TAUSBB=0.D0
        IF(KCNORM > 1) CALL GTSALB(XXG,XXT,SURX,SURFBB,COLPFG,TAUSBB,2)
        DBLEXT(NL+1,K)=TAUSBB
        ALLTAU=TAUSBB+COLEXT(K)
        CALL SGPGXG(COSZ,ALLTAU,COLPFG,SGPG)
cc      ALLGCB(K)=SGPG
        DO 70 L=L1,NL
        IF(SRCEXT(L,K) < TCLMIN) GO TO 70
        DBLGCB(L,K)=SGPG
   70   CONTINUE
      ELSE
        DO 80 L=L1,NL
        IF(SRCEXT(L,K) < TCLMIN) GO TO 80
        TAULAY=DBLEXT(L,K)
        GCBLAY=DBLGCB(L,K)
        CALL SGPGXG(COSZ,TAULAY,GCBLAY,SGPG)
        DBLGCB(L,K)=SGPG
   80   CONTINUE
      ENDIF
   90 CONTINUE

      WVCOL = SUM(ULGAS(:,1))
      ZWPATH = WVCOL*(1d0/COSZ + 2d0*srbalb(6))

#ifdef SWFIX_20151201
      FACK12 = 0.09325D0*
     &     ((ZWPATH**0.97D0)/(1.D0+5.D-4*(ZWPATH**1.31D0)))*0.462D-05
      FACK13 = 0.0001982D0*
     &     ((WVCOL**1.08D0)*(1.D0+6.D-5*(WVCOL**0.93D0)))*0.277D-05
#endif

      K = 0
  300 CONTINUE    !   DO K=1,NKSLAM
      K = K+1

      KLAM=KSLAM(K)
      DKS0X=DKS0(K)*S0COSZ

!     write(*,'(a,3i5,3(e12.4,1x))')'RADIATION1: ',
!    .     ILON,JLAT,K,DKS0(K),S0COSZ,DKS0X

      RBNB=SRBALB(KLAM)
      RBNX=SRXALB(KLAM)
      RCNB=0.D0
      RCNX=0.D0
      SRKINC(K)=DKS0X

      DO 200 N=L1,NL

      SRB(N)=RBNB
      SRX(N)=RBNX
      TLN=TLM(N)
      PLN=PL(N)
      ULN=ULGAS(N,1)

C     Select parameterized k-distribution gas absorption by H2O, O2, CO2
C     ------------------------------------------------------------------

      SELECT CASE (K)
      CASE (1)
C--------K=6-------H2O       DS0=.01
      TERMA=(35.66+TLN*(.0416-.0004622*TLN+.001057*PLN))*(1.+.04286*PLN)
      TERMB=(1.+.00171*ULN)*(1.+PLN*(189.088+.1316*PLN))
      IF(TERMB < 1000.) TERMB = 1000.
      TAU1 =TERMA/TERMB
      !IF(TAU1 > 0.02343) TAU1=0.02343
      IF(TAU1 > .05) TAU1 = .05
      TAU=TAU1*ULN

      CASE (2)
C--------K=5-------H2O       DS0=.03
      TERMA=(2.792+TLN*(.0914-.0002848*TLN+.0003395*PLN))
     +     *(1.+.02964*PLN)
      TERMB=(1.0+.000657*ULN)*(1.+PLN*(240.70+.13847*PLN))
      IF(TERMB < 1000.) TERMB = 1000.
      TAU1 =TERMA/TERMB
      !IF(TAU1 > 0.00520) TAU1=0.00520
      IF(TAU1 > .01) TAU1 = .01
      TAU=TAU1*ULN

      CASE (3)
C--------K=4-------H2O       DS0=.04
      TERMA=(.4768+.467E-04*PLN*TLN)*(1.+TLN*(.00191-.719E-05*TLN))
      TERMB=(1.+.717E-04*ULN)*(1.+PLN*(130.56+.0876*PLN))/(1.+.0266*PLN)
      IF(TERMB < 1000.) TERMB = 1000.
      TAU1 =TERMA/TERMB
      !IF(TAU1 > 0.00150) TAU1=0.0015
      IF(TAU1 > .01) TAU1 = .01
      TAU=TAU1*ULN

      CASE (4)
C--------K=3-------H2O       DS0=.04
      TERMA=(.000247*TLN-.091+PLN*(.00035+.78E-06*TLN))*(1.+.2847*PLN)
      TERMB=(1.+.2066E-04*ULN)*(1.+PLN*(137.17+.16132*PLN))
      IF(TERMA < 20.) TERMA = 20.
      IF(TERMB < 1000.) TERMB = 1000.
      TAU  =(TERMA/TERMB)*ULN

      CASE (5)
C--------K=2-------H2O       DS0=.04
      TERMA=(PLN*(1.974/TLN+.0001117*TLN)-10.713)*(1.+.005788*TLN)
     +     *(1.+.001517*PLN)
      TERMB=(1.+.3218E-04*ULN)*(1.+PLN*(863.44+.2048*PLN))
      IF(TERMA < 20.) TERMA = 20.
      IF(TERMB < 1000.) TERMB = 1000.
      TAU  =(TERMA/TERMB)*ULN

      CASE (6)
C--------K=4-------O2        DS0=.002
      ULN=ULGAS(N,4)
      TERMA=(.2236E-05-.1181E-09*TLN)*(1.+PLN*(.6364E-05*PLN+.001168))
      TERMB=1.+.1521E-05*ULN
      TAU  =(TERMA/TERMB)*ULN

      CASE (7)
C--------K=3-------O2        DS0=.004
      ULN=ULGAS(N,4)
      TERMA=(.3179E-06-.9263E-11*TLN)*(1.+PLN*(.8832E-05*PLN+.0005292))
      TERMB=1.+.1968E-06*ULN
      TAU  =(TERMA/TERMB)*ULN

      CASE (8)
C--------K=2-------O2        DS0=.013
      ULN=ULGAS(N,4)
      TERMA=(.2801E-07-.1638E-12*TLN)*(1.+PLN*(.1683E-04*PLN-.001721))
      TERMB=1.+.8097E-07*ULN
      TAU  =(TERMA/TERMB)*ULN

      CASE (9)
C--------K=4-------CO2       DS0=.002
      ULN=ULGAS(N,2)
      TERMA=(50.73-.03155*TLN-PLN*(.5543+.00091*TLN))*(1.-.1004*PLN)
      TERMB=(1.+.006468*ULN)*(1.+PLN*(49.51+.8285*PLN))
      TAU  =(TERMA/TERMB)*ULN
      IF(PLN < 175.0) TAU=(.00018*PLN+0.00001)*ULN

      CASE (10)
C--------K=3-------CO2       DS0=.003
      ULN=ULGAS(N,2)
      TERMA=(1.+.01319*TLN)*(PLN*(.008001*ULN+.4589E-03)-.8396*ULN)
      TERMB=ULN*(PLN+295.7+1.967*ULN)+.15126*PLN
      TAU  =(TERMA/TERMB)*ULN

      CASE (11)
C--------K=2-------CO2       DS0=.003
      ULN=ULGAS(N,2)
      TERMA=(1.+.02257*TLN)*(PLN*(.002295*ULN-.5489E-04)-.7571*ULN)
      TERMB=ULN*(PLN+803.9+2.477*ULN)-.09899*PLN
      TAU  =(TERMA/TERMB)*ULN

C     -------------------------------------------------------------------
C     fOnOff (default = 1.) scales SW long-path H2O absorption correction
C            =1. fully turned on, =0. disables correction (older version)
C            fOnOff is 'tunable' from 0. to 1. (introduced 7/3/2014)

      CASE (12)
      !ULN=ULGAS(N,1) ! not needed because uln set to this before select case
#ifdef SWFIX_20151201
      PTROOT=(((PLN+10.0)/1000.0)**0.5D0)/SQRT(TLN/296.D0)
      TAUK=PTROOT*ULN
      TAU=TAUK*FACK12*FONOFF
#else
      ALPH=0.002d0
      BETA=0.200d0
      !FACK12=1.05D-04*ZWPATH/(1.D0-1.D-05*ZWPATH)
      FACK12=.525D-04*ZWPATH/(1.D0+2.73D-04*ZWPATH)
      ROOT=SQRT(((PLN+50.0)/1000.0)**2+1000.0*BETA*ULN/(PLN+50.0))
      TAUK=ALPH*(ROOT-(PLN+50.0)/1000.0)
      TAU=TAUK*FACK12*fOnOff
#endif

      CASE (13)
      !ULN=ULGAS(N,1) ! not needed because uln set to this before select case
#ifdef SWFIX_20151201
      PTROOT=(((PLN+10.0)/1000.0)**0.5D0)/SQRT(TLN/296.D0)
      TAUK=PTROOT*ULN
      TAU=TAUK*FACK13*FONOFF
#else
      ALPH=0.004d0
      BETA=0.200d0
      !FACK13=1.05D-04*ZWPATH/(1.D0-1.D-05*ZWPATH)
      FACK13=.525D-04*ZWPATH/(1.D0+2.73D-04*ZWPATH)
      ROOT=SQRT(((PLN+50.0)/1000.0)**2+1000.0*BETA*ULN/(PLN+50.0))
      TAUK=ALPH*(ROOT-(PLN+50.0)/1000.0)
      TAU=TAUK*FACK13*fOnOff
#endif

      CASE (14)
        TAU=XCMNO2*ULGAS(N,5)+XCMO3*ULGAS(N,3)
      END SELECT

C     With 10 doublings to get to Tau=1.0, maximum seed tau is < 1/1024.
C     ------------------------------------------------------------------

      IF(TAU < 0.D0) TAU=0.D0

      TAU=TAU+DBLEXT(N,KLAM)
      IF(TAU < 1.D-06) GO TO 180
      PIZERO=DBLSCT(N,KLAM)/TAU
      IF(PIZERO < 0.001D0) GO TO 180

      PFF=DBLGCB(N,KLAM)

      NDBLS=0
      PR=1.D0-PFF
      PT=1.D0+PFF
      IF(TAU > 0.0019531D0) THEN
        DBLS=10.D0+1.44269D0*LOG(TAU)
        NDBLS=DBLS
        TAU=TAU/DBLN(NDBLS)
      ENDIF

C     Set optically thin limit values of R,T,X using PI0 renormalization
C     ------------------------------------------------------------------

      XANB=EXP(-TAU-TAU)
      XANX=EXP(-TAU*SECZ)
      TANB=PT*XANB
      XXT=(SECZ-2.D0)*TAU
      TANX=PT*SECZ
     +    *(.5D0+XXT*(.25D0+XXT*(.0833333D0+XXT*(.0208333D0+XXT))))*XANX
      RASB=PR*(1.D0-TAU*(2.D0-2.66667D0*TAU*(1.D0-TAU)))
      XXT=(SECZ+2.D0)*TAU
      RASX=PR*SECZ
     +    *(.5D0-XXT*(.25D0-XXT*(.0833333D0-XXT*(.0208333D0-XXT))))
      BNORM=(1.D0-XANB)/(RASB+TANB)*PIZERO
      XNORM=(1.D0-XANX)/(RASX+TANX)*PIZERO
      RASB=RASB*BNORM
      RASX=RASX*XNORM
      TANB=TANB*BNORM
      TANX=TANX*XNORM

C     Compute and record R,T,X atmospheric layer doubling/adding results
C     ------------------------------------------------------------------

      IF(NDBLS < 1) GO TO 170
      DO 160 NN=1,NDBLS
      RARB=RASB*RASB
      RARX=XANX*RASX
      XATB=XANB+TANB
      DENOM=1.D0-RARB
      DB=(TANB+XANB*RARB)/DENOM
      DX=(TANX+RARX*RASB)/DENOM
      UB=RASB*(XANB+DB)
      UX=RARX+RASB*DX
      RASB=RASB+XATB*UB
      RASX=RASX+XATB*UX
      TANB=XANB*TANB+XATB*DB
      TANX=XANX*TANX+XATB*DX
      XANB=XANB*XANB
      XANX=XANX*XANX
  160 CONTINUE
  170 CONTINUE
      RARB=RASB*RBNB
      RARX=RASB*RBNX
      XATB=XANB+TANB
      DENOM=1.D0-RARB
      DB=(TANB+XANB*RARB)/DENOM
      DX=(TANX+XANX*RARX)/DENOM
      UB=RBNB*(XANB+DB)
      UX=RBNX*XANX+RBNB*DX
      RBNB=RASB+XATB*UB
      RBNX=RASX+XATB*UX
      XATC=XATB/(1.D0-RASB*RCNB)
      RCNX=RASX+(XANX*RCNX+TANX*RCNB)*XATC
      RCNB=RASB+RCNB*XATB*XATC
      GO TO 190
  180 CONTINUE
      RASB=0.D0
      RASX=0.D0
      TANB=0.D0
      TANX=0.D0
      XANB=EXP(-TAU-TAU)
      XANX=EXP(-TAU*SECZ)
      DX=0.D0
      UX=RBNX*XANX
      RBNB=RBNB*XANB*XANB
      RBNX=UX*XANB
      RCNB=RCNB*XANB*XANB
      RCNX=RCNX*XANX*XANB
  190 CONTINUE
      RNB(N)=RASB
      RNX(N)=RASX
      TNB(N)=TANB
      TNX(N)=TANX
      XNB(N)=XANB
      XNX(N)=XANX
  200 CONTINUE

C     Record fluxes, spectral components at TOA, & top-layer bottom edge
C     ------------------------------------------------------------------

      SRDFLB(NL+1)=SRDFLB(NL+1)+DKS0X
      SRUFLB(NL+1)=SRUFLB(NL+1)+DKS0X*RBNX
      SRDFLB(NL)=SRDFLB(NL)+DKS0X*(XANX+DX)
      SRUFLB(NL)=SRUFLB(NL)+DKS0X*UX
      SKDFLB(NL+1,K)=DKS0X
      SKUFLB(NL+1,K)=DKS0X*RBNX
      SKDFLB(NL,K)=DKS0X*(XANX+DX)
      SKUFLB(NL,K)=DKS0X*UX
      RBXTOA=RBNX
      SRKALB(K)=RBNX

C     Add successively layer N (at bottom) to form upper composite layer
C     ------------------------------------------------------------------

      DO 230 N=NL-1,L1,-1
      XBNB=XNB(N)
      XBNX=XNX(N)
      RBNX=RNX(N)
      IF(RBNX > 1.D-05) GO TO 210
      RASB=RASB*XBNB*XBNB
      TANX=TANX*XBNB
      GO TO 220
  210 CONTINUE
      RBNB=RNB(N)
      TBNB=TNB(N)
      TBNX=TNX(N)
      RARB=RASB*RBNB
      XBTB=XBNB+TBNB
      DENOM=1.D0-RARB
      TANX=TBNX*XANX+XBTB*(TANX+XANX*RBNX*RASB)/DENOM
      RASB=RBNB+XBTB*XBTB*RASB/DENOM
  220 CONTINUE
      XANX=XANX*XBNX
      RBNB=SRB(N)
      RBNX=SRX(N)
      DX=(TANX+XANX*RBNX*RASB)/(1.D0-RASB*RBNB)
      UX=RBNX*XANX+RBNB*DX
      SRUFLB(N)=SRUFLB(N)+DKS0X*UX
      SRDFLB(N)=SRDFLB(N)+DKS0X*(XANX+DX)
      SKUFLB(N,K)=DKS0X*UX
      SKDFLB(N,K)=DKS0X*(XANX+DX)
  230 CONTINUE

C     Record absorbed spectral flux at ground for surface type fractions
C     ------------------------------------------------------------------

      SRKGAX(K,1:4)=DKS0X*XANX*(1.D0-PRNX(KLAM,1:4))
      SRKGAD(K,1:4)=DKS0X*  DX*(1.D0-PRNB(KLAM,1:4))

      IF(K==NKSLAM) GO TO 301
      SRINIR=SRINIR+DKS0X
      SRONIR=SRONIR+DKS0X*RBXTOA
      SRDNIR=SRDNIR+SKDFLB(1,K)
      SRUNIR=SRUNIR+SKUFLB(1,K)
      SRRNIR=SRRNIR+DKS0X*RCNX
      SRTNIR=SRTNIR+DKS0X*(TANX+XANX)
      SRXNIR=SRXNIR+DKS0X*XANX
      GO TO 300

  301 CONTINUE
      SRIVIS=DKS0X
      SROVIS=DKS0X*RBXTOA
      SRDVIS=SKDFLB(1,K)
      SRUVIS=SKUFLB(1,K)
      SRRVIS=DKS0X*RCNX
      SRTVIS=DKS0X*(TANX+XANX)
      SRXVIS=SRXVIS+DKS0X*XANX

!     write(*,'(a,3i5,3(e12.4,1x))')'RADIATION2: ',
!    . ILON,JLAT,K,XANX,DKS0X,SRXVIS

C     ------------------------------------------------------------------
C     UV absorption by O3 and O2 within solar spectral band DKS0(15)=.05
C     ------------------------------------------------------------------

      K=15
      DKS0X=DKS0(K)*S0COSZ
      SRKINC(K)=DKS0X
!     write(*,'(a,3i5,3(e12.4,1x))')'RADIATION3: ',
!    . ILON,JLAT,K,DKS0(K),S0COSZ,DKS0X

      N=NL+1
      ATOPX=0.D0
      ATOPD=0.D0
      O3CMX=0.D0
      O3CMD=0.D0
  302 CONTINUE
      N=N-1
      O3CMX=O3CMX+COSMAG*ULGAS(N,3)
      O3CMD=O3CMD+1.90D0*ULGAS(N,3)
      CALL AO3ABS(O3CMX,ABOTX)
      CALL AO3ABS(O3CMD,ABOTD)
      AO3X(N)=(ABOTX-ATOPX)/DKS0(15)
      AO3D(N)=(ABOTD-ATOPD)/DKS0(15)
      ATOPX=ABOTX
      ATOPD=ABOTD
      IF(N > L1) GO TO 302
  303 CONTINUE
      O3CMX=O3CMX+1.90D0*ULGAS(N,3)
      O3CMD=O3CMD+1.90D0*ULGAS(N,3)
      CALL AO3ABS(O3CMX,ATOPX)
      CALL AO3ABS(O3CMD,ATOPD)
      AO3UXN=(ATOPX-ABOTX)/DKS0(15)
      AO3UDN=(ATOPD-ABOTD)/DKS0(15)
      AO3U(N)=XNX(N)*AO3UXN+(1.D0-XNX(N))*AO3UDN
      ABOTX=ATOPX
      ABOTD=ATOPD
      N=N+1
      IF(N < NL+1) GO TO 303
      RBNB=SRBALB(KLAM)
      RBNX=SRXALB(KLAM)
      RCNB=0.D0
      RCNX=0.D0
C                                  -------------------------------------
C                                  Get Oxygen UV absorption contribution
C                                  -------------------------------------
C----------------
      CALL GETO2A
C----------------
C                 ------------------------------------------------------
C                 Add Layers from Ground up. Retain Composite RBNB, RBNX
C             R,T,X of "A" (above) layer are corrected for O3 absorption
C             ----------------------------------------------------------

      DO 304 N=L1,NL
      O2FHRL(N)=O2FHRL(N)/DKS0(15)*FULGAS(4)
      O2FHRB(N)=O2FHRB(N)/DKS0(15)*FULGAS(4)
      SRB(N)=RBNB
      SRX(N)=RBNX
      XANX=XNX(N)*(1.D0-AO3X(N)-O2FHRL(N))
      XANB=XNB(N)*(1.D0-AO3D(N)-O2FHRB(N))
      RASX=RNX(N)*(1.D0-AO3U(N))
      RASB=RNB(N)*(1.D0-AO3U(N))
      TANX=TNX(N)*(1.D0-AO3D(N))
      TANB=TNB(N)*(1.D0-AO3D(N))
!nu      ABSRTX=1.D0-XANX-TANX-RASX
!nu      ABSRTB=1.D0-XANB-TANB-RASB
      RARB=RASB*RBNB
      RARX=RASB*RBNX
      XATB=XANB+TANB
      DENOM=1.D0-RARB
      DB=(TANB+XANB*RARB)/DENOM
      DX=(TANX+XANX*RARX)/DENOM
      UB=RBNB*(XANB+DB)
      UX=RBNX*XANX+RBNB*DX
      RBNB=RASB+XATB*UB
      RBNX=RASX+XATB*UX
      XATC=XATB/(1.D0-RASB*RCNB)
      RCNX=RASX+(XANX*RCNX+TANX*RCNB)*XATC
      RCNB=RASB+RCNB*XATB*XATC
  304 CONTINUE
      VRD(NL+1)=1.D0
      VRU(NL+1)=RBNX
      SRKALB(15)=RBNX
      N=NL
      VRD(N)=XANX+DX
      VRU(N)=UX
  310 CONTINUE
      N=N-1
      XBNX=XNX(N)*(1.D0-AO3X(N)-O2FHRL(N))
      XBNB=XNB(N)*(1.D0-AO3D(N)-O2FHRB(N))
      RBNX=RNX(N)*(1.D0-AO3U(N))
      RBNB=RNB(N)*(1.D0-AO3U(N))
      TBNX=TNX(N)*(1.D0-AO3D(N))
      TBNB=TNB(N)*(1.D0-AO3D(N))

C     Add successively layer N (at bottom) to form upper composite layer
C     ------------------------------------------------------------------

      RARB=RASB*RBNB
      XBTB=XBNB+TBNB
      DENOM=1.D0/(1.D0-RARB)
      TANX=TBNX*XANX+XBTB*(TANX+XANX*RBNX*RASB)*DENOM
      RASB=RBNB+XBTB*XBTB*RASB*DENOM
      XANX=XANX*XBNX

C     Add upper & bottom composite layers to get flux at layer interface
C     ------------------------------------------------------------------

      RBNB=SRB(N)
      RBNX=SRX(N)
      DX=(TANX+XANX*RBNX*RASB)/(1.D0-RASB*RBNB)
      UX=RBNX*XANX+RBNB*DX
      VRD(N)=XANX+DX
      VRU(N)=UX
      IF(N > 1) GO TO 310
      SRKGAX(15,1:4)=DKS0X*XANX*(1-PRNX(6,1:4))
      SRKGAD(15,1:4)=DKS0X*  DX*(1-PRNB(6,1:4))

      DO 325 N=L1,NL+1
      VRD(N)=VRD(N)*DKS0X
      VRU(N)=VRU(N)*DKS0X
      SKDFLB(N,K)=VRD(N)
      SKUFLB(N,K)=VRU(N)
  325 CONTINUE
      SRIVIS=SRIVIS+VRD(NL+1)
      SROVIS=SROVIS+VRU(NL+1)
      PLAVIS=SROVIS/SRIVIS
      SRDVIS=SRDVIS+VRD(L1)
      SRUVIS=SRUVIS+VRU(L1)
      ALBVIS=SRUVIS/(SRDVIS+1.D-10)
      SRRVIS=SRRVIS+DKS0X*RCNX
      SRTVIS=SRTVIS+DKS0X*(TANX+XANX)
      SRXVIS=SRXVIS+DKS0X*XANX
      SRAVIS=1.D0-SRRVIS-SRTVIS

C     K16 strong absorbing contributions are computed without scattering
C     ------------------------------------------------------------------

      K=16
      DKS0X=DKS0(16)*S0COSZ
      SRKINC(16)=DKS0X
      SRKA16=0.D0
      SRKGAX(16,1:4)=0.D0
      SRKGAD(16,1:4)=0.D0
      DO 345 KK=1,3
      IF(KK==1) DKS0XX=DKS0X*0.002D0/0.011D0
      IF(KK==2) DKS0XX=DKS0X*0.008D0/0.011D0
      IF(KK==3) DKS0XX=DKS0X*0.001D0/0.011D0
      TRNC=1.D0
      DO 330 N=NL,L1,-1
      PLN=PL(N)
      CLX=DBLEXT(N,1)-DBLSCT(N,1)

C--------K=5-------CO2       DS0=.002
      IF(KK==1) THEN
      TRN1=0.D0
      ULN=ULGAS(N,2)*SECZ
      IF(ULN > 7.D0) ULN=7.D0
      TERMA=.003488*PLN*(1.+39.59*EXP(-8.769*ULN/(1.+4.419*ULN)))
     +     *(1.+ULN*(.001938*PLN-.00503*ULN))
      TERMB=(1.+.04712*PLN*(1.+.4877*ULN))
      TAUG=TERMA/TERMB*ULN
      TAU1=TAUG+CLX*SECZ
      IF(TAU1 < 10.0) TRN1=EXP(-TAU1)
      FAC(N)=TRN1
      ENDIF

C--------K=7-------H2O       DS0=.008
      IF(KK==2) THEN
      TRN2=0.D0
      ULN=ULGAS(N,1)*SECZ
      TERMA=.001582*PLN*(1.+6.769*EXP(-9.59*ULN/(1.+5.026*ULN)))
     +     *(1.+ULN*(.2757E-03*PLN+.001429*ULN))
      TERMB=(1.+.003683*PLN*(1.+1.187*ULN))
      TAUG=TERMA/TERMB*ULN
      TAU2=TAUG+CLX*SECZ
      IF(TAU2 < 10.0) TRN2=EXP(-TAU2)
      FAC(N)=TRN2
      ENDIF

C--------K=5-------O2        DS0=.001
      IF(KK==3) THEN
      TRN3=0.D0
      ULN=ULGAS(N,4)*SECZ
      TERMA=(.1366E-03-.2203E-07*TLN)*(1.+PLN*(.1497E-06*ULN+.001261))
      TERMB=(1.+.3867E-03*ULN)/(1.+.2075E-04*ULN)
      TAUG=TERMA/TERMB*ULN
      TAU3=TAUG+CLX*SECZ
      IF(TAU3 < 10.0) TRN3=EXP(-TAU3)
      FAC(N)=TRN3
      ENDIF

      TRNC=TRNC*FAC(N)
      SRDFLB(N)=SRDFLB(N)+DKS0XX*TRNC
      SKDFLB(N,K)=SKDFLB(N,K)+DKS0XX*TRNC
  330 CONTINUE
      SRDFLB(NL+1)=SRDFLB(NL+1)+DKS0XX
      SRUFLB(L1)=SRUFLB(L1)+DKS0XX*TRNC*SRXALB(1)
      SKDFLB(NL+1,K)=SKDFLB(NL+1,K)+DKS0XX
      SKUFLB(L1,K)=SKUFLB(L1,K)+DKS0XX*TRNC*SRXALB(1)

C     For completeness, any incident flux at ground is relflected upward
C     ------------------------------------------------------------------

      TRNU=TRNC
      DO 340 N=L1+1,NL+1
      TRNU=TRNU*FAC(N-1)
      SRUFLB(N)=SRUFLB(N)+DKS0XX*TRNC*SRXALB(1)*TRNU
      SKUFLB(N,K)=SKUFLB(N,K)+DKS0XX*TRNC*SRXALB(1)*TRNU
  340 CONTINUE
      SRKGAX(16,1:4)=SRKGAX(16,1:4)+DKS0XX*TRNC*(1-PRNX(1,1:4))
      SRKA16=SRKA16+TRNU*SRXALB(1)

      SRINIR=SRINIR+DKS0XX
      SRONIR=SRONIR+DKS0XX*TRNU*SRXALB(1)
      SRDNIR=SRDNIR+SKDFLB(L1,K)
      SRUNIR=SRUNIR+SKUFLB(L1,K)
  345 CONTINUE
      PLANIR=SRONIR/SRINIR
      ALBNIR=SRUNIR/(SRDNIR+1.D-10)
      SRKALB(16)=SRKA16/DKS0X

      SRDFLB(L1:NL+1) = SRDFLB(L1:NL+1) + VRD(L1:NL+1)
      SRUFLB(L1:NL+1) = SRUFLB(L1:NL+1) + VRU(L1:NL+1)
      SRNFLB(L1:NL+1) = SRDFLB(L1:NL+1) - SRUFLB(L1:NL+1)
      SRFHRL(L1:NL)   = SRNFLB(L1+1:NL+1) - SRNFLB(L1:NL)
      SRRNIR=SRRNIR+DKS0X*RCNX
      SRTNIR=SRTNIR+DKS0X*(TANX+XANX)
      SRXNIR=SRXNIR+DKS0X*XANX

      S0VIS=0.53D0*S0
      SRTVIS=SRTVIS/S0VIS
      SRRVIS=SRRVIS/S0VIS
      SRXVIS=SRXVIS/S0VIS
      SRAVIS=1.D0-SRTVIS-SRRVIS

      S0NIR=0.47D0*S0
      SRTNIR=SRTNIR/S0NIR
      SRRNIR=SRRNIR/S0NIR
      SRXNIR=SRXNIR/S0NIR
      SRANIR=1.D0-SRTNIR-SRRNIR


C     ------------------------------------------------------------------
C     FSRNFG defines the total solar flux absorbed at the ground surface
C              taking into account the albedo of different surface types
C     Thus:
C              SRNFLB(1)=POCEAN*FSRNFG(1)+PEARTH*FSRNFG(2)
C                       + POICE*FSRNFG(3)+ PLICE*FSRNFG(4)
C
C     NOTE:    If any surface type POCEAN, PEARTH, POICE, PLICE are Zero
C              the corresponding FSRNFG(I) absorbed solar flux at ground
C              is computed with that surface-type albedo set equal to 0.
C              ---------------------------------------------------------

      DO 420 I=1,4
      FSRNFG(I) = sum(SRKGAX(1:16,I)) + sum(SRKGAD(1:16,I))
  420 CONTINUE


      DO 510 K=1,16
      SKNFLB(L1:NL+1,K)=SKDFLB(L1:NL+1,K)-SKUFLB(L1:NL+1,K)
  510 CONTINUE

      DO 530 K=1,16
      SKFHRL(L1:NL,K)=SKNFLB(L1+1:NL+1,K)-SKNFLB(L1:NL,K)
  530 CONTINUE

      DO 560 L=L1,NL+1
      SKDFLB(L,17) = sum (SKDFLB(L,1:16))
      SKUFLB(L,17) = sum (SKUFLB(L,1:16))
      SKNFLB(L,17) = sum (SKNFLB(L,1:16))
  560 CONTINUE
      DO L=L1,NL
        SKFHRL(L,17) = sum (SKFHRL(L,1:16))
      END DO

      RETURN
      END SUBROUTINE SOLARM




      SUBROUTINE GETMIE(NA,AREFF,SQEX,SQSC,SQCB,TQAB,Q55)

c     INCLUDE 'rad00def.radCOMMON.f'

      INTEGER, INTENT(IN) :: NA
      REAL*8,  intent(in) :: areff
      REAL*8   SQEX(6),SQSC(6),SQCB(6),TQEX(33),TQSC(33),TQAB(33),Q55
      REAL*8   QXAERN(25),QSAERN(25),QGAERN(25),Q55AER(25)

      REAL*8 wts,wta,QGAERX,pi,vreff
      INTEGER n0,k,n,nn
                         !                               1   2   3   4
      IF(NA < 5) THEN    !    NA : Aerosol compositions SO4,SEA,ANT,OCX
        N0=0
        IF(NA==2) N0=22
        IF(NA==3) N0=44
        IF(NA==4) N0=88
        DO 112 K=1,6
        DO 111 N=1,22
        NN=N0+N
        WTS=FRSULF(NA)
        WTA=1.D0-WTS
        QXAERN(N)=SRUQEX(K,NN)*WTA+SRUQEX(K,N)*WTS
        QSAERN(N)=SRUQSC(K,NN)*WTA+SRUQSC(K,N)*WTS
        QGAERX=SRUQCB(K,NN)*SRUQSC(K,NN)*WTA+SRUQCB(K,N)*SRUQSC(K,N)*WTS
        QGAERN(N)=QGAERX/QSAERN(N)
  111   CONTINUE
        CALL SPLINE(REFU22,QXAERN,22,AREFF,SQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFU22,QSAERN,22,AREFF,SQSC(K),1.D0,1.D0,1)
        CALL SPLINE(REFU22,QGAERN,22,AREFF,SQCB(K),1.D0,1.D0,1)

        PI=SQSC(K)/SQEX(K)
        IF(PI > PI0MAX(NA)) SQSC(K)=SQSC(K)*PI0MAX(NA)/PI
  112   CONTINUE
        DO 114 K=1,33
        DO 113 N=1,22
        NN=N0+N
        WTS=FRSULF(NA)
        WTA=1.D0-WTS
        QXAERN(N)=TRUQEX(K,NN)*WTA+TRUQEX(K,N)*WTS
        QSAERN(N)=TRUQSC(K,NN)*WTA+TRUQSC(K,N)*WTS
        QGAERX=TRUQCB(K,NN)*TRUQSC(K,NN)*WTA+TRUQCB(K,N)*TRUQSC(K,N)*WTS
        QGAERN(N)=QGAERX/(QSAERN(N)+1.d-20)
  113   CONTINUE
        CALL SPLINE(REFU22,QXAERN,22,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFU22,QSAERN,22,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
  114   CONTINUE
        DO 115 N=1,22
        NN=N0+N
        WTS=FRSULF(NA)
        WTA=1.D0-WTS
        Q55AER(N)=Q55U22(NN)*WTA+Q55U22(N)*WTS
  115   CONTINUE
        CALL SPLINE(REFU22,Q55U22,22,AREFF,Q55,1.D0,1.D0,1)
      ENDIF

                               !                              5   6
      IF(NA==5.OR.NA==6) THEN  !   NA : Aerosol compositions BIC,BCB
cc      AREFF=REFDRY(NA)
        DO 122 K=1,6
        QXAERN(:)=SRSQEX(K,:)    ! 1:25
        QSAERN(:)=SRSQSC(K,:)    ! 1:25
        QGAERN(:)=SRSQCB(K,:)    ! 1:25
        CALL SPLINE(REFS25,QXAERN,25,AREFF,SQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFS25,QSAERN,25,AREFF,SQSC(K),1.D0,1.D0,1)
        CALL SPLINE(REFS25,QGAERN,25,AREFF,SQCB(K),1.D0,1.D0,1)
  122   CONTINUE
        DO 124 K=1,33
        QXAERN(:)=TRSQEX(K,:)    ! 1:25
        QSAERN(:)=TRSQSC(K,:)    ! 1:25
        QGAERN(:)=TRSQCB(K,:)    ! 1:25
        CALL SPLINE(REFS25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFS25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
  124   CONTINUE
        CALL SPLINE(REFS25,Q55S25,25,AREFF,Q55,1.D0,1.D0,1)
      ENDIF

                                      !                             7
      IF(NA==7) THEN                  !   NA : Aerosol composition DST
cc      AREFF=REFDRY(NA)
        DO 132 K=1,6
        QXAERN(:)=SRDQEX(K,:)    ! 1:25
        QSAERN(:)=SRDQSC(K,:)    ! 1:25
        QGAERN(:)=SRDQCB(K,:)    ! 1:25
        CALL SPLINE(REFD25,QXAERN,25,AREFF,SQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFD25,QSAERN,25,AREFF,SQSC(K),1.D0,1.D0,1)
        CALL SPLINE(REFD25,QGAERN,25,AREFF,SQCB(K),1.D0,1.D0,1)
  132   CONTINUE
        DO 134 K=1,33
        QXAERN(:)=TRDQEX(K,:)    ! 1:25
        QSAERN(:)=TRDQSC(K,:)    ! 1:25
        QGAERN(:)=TRDQCB(K,:)    ! 1:25
        CALL SPLINE(REFD25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFD25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
  134   CONTINUE
        CALL SPLINE(REFD25,Q55D25,25,AREFF,Q55,1.D0,1.D0,1)
      ENDIF

                             !                                      8
      IF(NA==8) THEN         !     NA : Aerosol composition(H2SO4) VOL
        VREFF=AREFF
        IF(VREFF < 0.1D0) VREFF=0.1D0
        IF(VREFF > 2.0D0) VREFF=2.0D0
        CALL GETQVA(VREFF)
        SQEX(:)=QVH2S(:)       ! 1:6
        SQSC(:)=SVH2S(:)       ! 1:6
        SQCB(:)=GVH2S(:)       ! 1:6
        TQAB(:)=AVH2S(:)       ! 1:33
        Q55=Q55H2S
      ENDIF
      RETURN
      END SUBROUTINE GETMIE

      SUBROUTINE AO3ABS(OCM,O3ABS)
      IMPLICIT NONE
C              ---------------------------------------------------------
C              UV absorption by Ozone  is expressed as a fraction of the
C              total solar flux S0. Hence O3ABS (fraction of total solar
C              flux absored by OCM cm ofozone) must be normalized within
C              SOLARM by dividing O3ABS by the corresponding fraction of
C              the solar flux within the spectral interval DKS0(15)=0.05
C              ---------------------------------------------------------
      REAL*8, INTENT(IN) :: OCM
      REAL*8, INTENT(OUT) :: O3ABS
      REAL*8 XX,DX
      INTEGER IP,IX

      O3ABS=AO3(460)
      IP=0
      XX=OCM*1.D+04
      IX=XX
      IF(IX > 99) GO TO 110
      IF(IX < 1 ) GO TO 130
      GO TO 120
  110 CONTINUE
      IP=IP+90
      XX=XX*0.1D0
      IX=XX
      IF(IX > 99) GO TO 110
  120 CONTINUE
      DX=XX-IX
      IX=IX+IP
      IF(IX > 459) GO TO 140
      O3ABS=AO3(IX)+DX*(AO3(IX+1)-AO3(IX))
      GO TO 140
  130 CONTINUE
      O3ABS=XX*AO3(1)
  140 CONTINUE

      RETURN
      END SUBROUTINE AO3ABS

      SUBROUTINE WRITER(KWRU,INDEX)
C
c      USE SURF_ALBEDO, only : AVSCAT, ANSCAT, AVFOAM, ANFOAM,
c     *     WETTRA, WETSRA, ZOCSRA, ZSNSRA, ZICSRA, ZDSSRA, ZVGSRA,
c     *     EOCTRA, ESNTRA, EICTRA, EDSTRA, EVGTRA, AGEXPF, ALBDIF
      USE SURF_ALBEDO, only : get_albedo_data
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE DustParam_mod, only : redust
      IMPLICIT NONE
C
C     ------------------------------------------------------------------
C     WRITER Radiative Input/Output Cloud/Aerosol Data/Conrol Parameters
C
C         INDEX
C           0     control parameter defaults in RADPAR
C           1     RADPAR Radiative control/scaling params; GHG defaults
C           2     RADPAR Atmospheric composition P,H,T,Cld,Aer  profiles
C           3     RADPAR Computed LW SW fluxes cooling and heating rates
C           4     Aerosol and Cloud: Mie scattering radiative parameters
C                 A  SW aerosol Mie scattering Qx,Qs,g in use parameters
C                 B  SW  cloud  Mie scattering Qx,Qs,g in use parameters
C                 C  SW cld+aer Mie scattering Qx,Qs,g in use parameters
C                 D  SW LW aerosol 11-compositon  Mie Qx,Qs,g parameters
C                 E  SW LW aerosol  6-compositon  Mie Qx,Qs,g parameters
C                 F  SW LW aerosol 8-size D dust  Mie Qx,Qs,g parameters
C                 G  SW LW  cloud  15-size/phase  Mie Qx,Qs,g parameters
C           5     LW cld,aer,gas total optical k-distribution extinction
C           6     LW gas absorb: total optical k-distribution extinction
C           7     A  LW cloud   TRCALK optical k-distribution extinction
C                 B  LW aerosol TRAALK optical k-distribution extinction
C           8     SW Spectral/k-dist flux, albedo, absorption components
C                 A  Spectral components of downward & upward solar flux
C                 B  Spectral components of net solar flux, heating rate
C           9     LW flux contribution from each k-distribution interval
C                 1  Downward LW flux  from each k-distribution interval
C                 2  Upward   LW flux  from each k-distribution interval
C                 3  Net (Up) LW flux  from each k-distribution interval
C                 4  Flux cooling rate from each k-distribution interval
C                 5  Fraction coolrate from each k-distribution interval
C        NOTE:
C                 KWTRAB sets LW Mie parameters in 4-D,E,F,G
C                        KWTRAB=0 (default) sets LW output to be Mie Qab
C                        KWTRAB=1           sets LW output to be Mie Qex
C                        KWTRAB=2           sets LW output to be Mie Qsc
C                        KWTRAB=3           sets LW output to be Mie Qcb
C                        KWTRAB=4           sets LW output to be Mie Pi0
C
C                 INDEX  0-9 : show item 'INDEX' only
C                 INDEX 11-19: show items 1->last digit of 'INDEX'
C                 INDEX 21-29: show items 0->last digit of 'INDEX'
C                 KWRU directs the output to selected (KWRU) file number
C     ------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: INDEX
      REAL*8 AVSCAT, ANSCAT, AVFOAM, ANFOAM,
     *     WETTRA, WETSRA, ZOCSRA, ZSNSRA, ZICSRA, ZDSSRA, ZVGSRA,
     *     EOCTRA, ESNTRA, EICTRA, EDSTRA, EVGTRA,
     &     AGEXPF(3,2), ALBDIF(3,2)
      REAL*8, dimension(10) ::       ! no longer needed except in writer
C!nu               TROPOSPHERIC AEROSOL effective radius
C!nu               BCI  OCI  SUI  SEA  SUN    ANT  OCN  OCB  BCB  SSB
     *   REAERO=(/ 0.1, 0.3, 0.3, 2.0, 0.3,   1.0, 0.3, 0.3, 0.2, 0.5/)
      CHARACTER*8, PARAMETER :: FTYPE(5) =
     *     (/'DOWNWARD','  UPWARD','UPWD NET','COOLRATE','FRACTION'/)
      CHARACTER*6, PARAMETER :: GHG(12) =
     *     (/'   H2O','   CO2','    O3','    O2','   NO2','   N2O'
     +      ,'   CH4','CCL3P1','CCL2P2','    N2',' CFC-Y',' CFC-Z'/)
      CHARACTER*3 TRABCD(5),TRAXSG(5),snotyp
      DATA TRABCD/'TRA','TRB','TRC','TRD','TRE'/
      DATA TRAXSG/'QAB','QEX','QSC','QCB','PI0'/

      REAL*8 TKEFF(3),TRPI0K(25)
      REAL*8 WFLB(LX,33),WFSL(33),UXGAS(LX,9)
      REAL*8 BGFLUX(33),BGFRAC(33),TAUSUM(33)
      REAL*8 SUM0(20),SUM1(LX+1),SUM2(LX+1),SUM3(LX+1)
      REAL*8, DIMENSION(LX,6) :: WSREXT,WSRSCT,WSRGCB,WSRPI0
      REAL*8 FSR1(17),FSR2(17)
      INTEGER :: ISR1(16), KWRU
      INTEGER, PARAMETER ::
     *  KSLAMW(16) = (/ 1, 1, 2, 2, 5, 5, 5, 5, 1, 1, 1, 3, 4, 6, 6, 1/)
     * ,IORDER(16) = (/12,11,10, 9, 6, 5, 4, 3,15,14,13, 8, 7, 2, 1,16/)

      character*1,parameter :: AUXGAS(4) = (/'0','L','X','X'/)
      REAL*8, PARAMETER :: P0=1013.25,SIGMA=5.6697D-08
      REAL*8 ACOLX,BCOLX,DCOLX,VCOLX,TCOLX,FACTOR,PPMCO2,PPMO2
     *     ,PPMN2O,PPMCH4,PPMF11,PPMF12,PPMY11,PPMZ12,EPS,TAER,HLM,TLAPS
     *     ,TAU55,TGMEAN,PSUM,SRALB,STNFLB,CRHRF,STFHR,TRDCR,SRDHR,STDHR
     *     ,PFW,DPF,FRACSL,SIGT4,WTG,SUMK,SUMT,SUMK1,SUMK2
     *     ,ASUM1,BSUM1,CSUM1,DSUM1,ESUM1,FSUM1,ASUM2,BSUM2,CSUM2,DSUM2
     *     ,ESUM2,FSUM2,ASUM3,BSUM3,CSUM3,DSUM3,ESUM3,FSUM3,SUML,SUMA
     *     ,SUMB,SUMC,SUMD,SUME,SUMF
      INTEGER I,J,K,L,KW,INDJ,INDI,INDX,KPAGE,NPAGE,LUXGAS,LGS,IPI0,IRHL
     *     ,N,II,IPF,ITG,LK,KK,NW,LINFIL

      call get_albedo_data( AVSCAT, ANSCAT, AVFOAM, ANFOAM,
     *     WETTRA, WETSRA, ZOCSRA, ZSNSRA, ZICSRA, ZDSSRA, ZVGSRA,
     *     EOCTRA, ESNTRA, EICTRA, EDSTRA, EVGTRA, AGEXPF, ALBDIF )

      KW=KWRU
      INDJ=MOD(INDEX,10)
      IF(INDJ < 1.and.INDEX > 0) INDJ=10
      INDI=1
      IF(INDEX > 20.or.INDEX==0) INDI=0
      IF(INDEX < 11) INDI=INDJ

      if (INDJ > 0) then
        DO K=1,6
        DO L=L1,NL
          WSREXT(L,K)=SRAEXT(L,K)+SRBEXT(L,K)
     +               +SRDEXT(L,K)+SRVEXT(L,K)
          WSRSCT(L,K)=SRASCT(L,K)+SRBSCT(L,K)
     +               +SRDSCT(L,K)+SRVSCT(L,K)
          WSRGCB(L,K)=SRASCT(L,K)*SRAGCB(L,K)+SRBSCT(L,K)*SRBGCB(L,K)
     +               +SRDSCT(L,K)*SRDGCB(L,K)+SRVSCT(L,K)*SRVGCB(L,K)
          WSRPI0(L,K)=WSRSCT(L,K)/(WSREXT(L,K)+1.E-10)
          WSRGCB(L,K)=WSRGCB(L,K)/(WSRSCT(L,K)+1.D-10)
        END DO
        END DO
C
        ACOLX = sum (SRAEXT(L1:NL,6))
        BCOLX = sum (SRBEXT(L1:NL,6))
        DCOLX = sum (SRDEXT(L1:NL,6))
        VCOLX = sum (SRVEXT(L1:NL,6))
        TCOLX = ACOLX+BCOLX+DCOLX+VCOLX
      end if

      DO 9999 INDX=INDI,INDJ

      KPAGE=1
      IF(INDX==0) GO TO 90

      GO TO (100,200,300,400,500,600,700,800,900,1000),INDX
C
C-------------
   90 CONTINUE
C-------------
      IF (AM_I_ROOT()) THEN
        WRITE(KW,6000)
 6000 FORMAT(' CALL WRITER(KW,0) :',2X,'PAGE 1/2  '
     +          ,'CONTROL PARAMS   DEFINITIONS'/
     +      /' CONTROL PARAMTER      DEFAULT  PARAMETER DESCRIPTION')

       WRITE(KW,6001)                              KUVFAC,KSNORM
     + ,KWTRAB,KGGVDF,KPGRAD,KLATZ0,KCLDEM,KANORM,KFPCO2,KPFOZO,KSIALB
     + ,KORDER,KUFH2O,KUFCO2,KCSELF,KCFORN
 6001 FORMAT( ! 7X,'   KVRAER = ',I1,'     1      Repartition Aer VDist'
!nu  2    ! /7X,'   MEANAC = ',I1,'     0      Use Ann-Mean Aer Clim'
!nu  3    ! /7X,'   MEANDD = ',I1,'     0      Use Ann-Mean Des Dust'
!nu  4    ! /7X,'   MEANVA = ',I1,'     0      Use Ann-Mean Volc Aer'
!nu  5      /7X,'   NCARO3 = ',I1,'     0      NCAR London 1976 Ozon'/
     6       7X,'   KUVFAC = ',I1,'     0      ON/OFF UV Mult Factor'
     7      /7X,'   KSNORM = ',I1,'     0      Norm S0 when KUVFAC=1'
     8      /7X,'   KWTRAB = ',I1,'     0      WRITER: Qab,Qex,Qsc,g'
     9      /7X,'   KGGVDF = ',I1,'     0      Use GHG VertProf Grad'
     A      /7X,'   KPGRAD = ',I1,'     1      Pole-to-Pole GHG Grad'
     1      /7X,'   KLATZ0 = ',I1,'     1      Use GHG VDist Lat Dep'
     2      /7X,'   KCLDEM = ',I1,'     1      Use TopCloud Scat Cor'
     3      /7X,'   KANORM = ',I1,'     0      Use SGP Atmo Col Norm'
     4      /7X,'   KFPCO2 = ',I1,'     0      1=MOD CO2PROF: FPXCO2'
     5      /7X,'   KPFOZO = ',I1,'     0      1=MOD O3 PROF: FPXOZO'
     6      /7X,'   KSIALB = ',I1,'     0      Schramm"s ocn ice alb'
     7      /7X,'   KORDER = ',I1,'     0      WRITER k-d spec order'
     8      /7X,'   KUFH2O = ',I1,'     1      Col Absorber Scal H2O'
     +      /7X,'   KUFCO2 = ',I1,'     1      Col Absorber Scal CO2'
     1      /7X,'   KCSELF = ',I1,'     1      H2O Cont Self-Broaden'
     2      /7X,'   KCFORN = ',I1,'     1      H2O Con Foreign-Broad'
     R      )
C
      WRITE(KW,6004)
 6004 FORMAT(/' CONTROL PARAMTER    DEFAULT       SNOW/ICE FACTORS')

      WRITE(KW,6005) agexpf,albdif
 6005 FORMAT(7X   ,'   AGEXPF = ',F7.3,'      SNOWAGE XPFACTOR SH EARTH'
     A      /7X   ,'    SH O  = ',F7.3,'         "       "     "  OCICE'
     B      /7X   ,'    SH L  = ',F7.3,'         "       "     "  LDICE'
     C      /7X   ,'    NH E  = ',F7.3,'         "       "     NH EARTH'
     D      /7X   ,'    NH O  = ',F7.3,'         "       "     "  OCICE'
     E      /7X   ,'    NH L  = ',F7.3,'         "       "     "  LDICE'
     F      /7X   ,'   ALBDIF = ',F7.3,'      SNOW/ICE ALBDIF  SH EARTH'
     G      /7X   ,'    SH O  = ',F7.3,'         "       "     "  OCICE'
     H      /7X   ,'    SH L  = ',F7.3,'         "       "     "  LDICE'
     I      /7X   ,'    NH E  = ',F7.3,'         "       "     NH EARTH'
     J      /7X   ,'    NH O  = ',F7.3,'         "       "     "  OCICE'
     K      /7X   ,'    NH L  = ',F7.3,'         "       "     "  LDICE'
     L      )
C
      WRITE(KW,6006)
 6006 FORMAT('0CONTROL PARAMTER    VALUE',16X,' DEFAULT')
      WRITE(KW,6007) REFF0,VEFF0, AVSCAT,ANSCAT,AVFOAM,ANFOAM
 6007 FORMAT(7X,'   REFF0  = ',F7.3,'                 0.300         '
     B      /7X,'   VEFF0  = ',F7.3,'                 0.350         '
     H      /7X,'   AVSCAT = ',F7.5,'                 0.01560       '
     I      /7X,'   ANSCAT = ',F7.5,'                 0.00020       '
     J      /7X,'   AVFOAM = ',F7.5,'                 0.21970       '
     K      /7X,'   ANFOAM = ',F7.5,'                 0.15140       '
     X      )
      WRITE(KW,6008)
 6008 FORMAT(/10X,'UV Solar Flux Spectral Partitions and Factors')
      WRITE(KW,6009) UVWAVL,UVFACT
 6009 FORMAT(10X,'UVWAVL = ',F7.5,2F8.5/10X,'UVFACT = ',F7.5,2F8.5)
C
!nu   WRITE(KW,6013)
 6013 FORMAT(/' CONTROL PARAMETER  PI0VIS     PI0TRA      DEFAULT')
!nu   WRITE(KW,6014) PI0VIS,PI0TRA
 6014 FORMAT(7X,'   ACID1 = ',F8.6,F11.6,'     1.0                 '
     A      /7X,'   SSALT = ',F8.6,F11.6,'     1.0                 '
     B      /7X,'   SLFT1 = ',F8.6,F11.6,'     1.0                 '
     C      /7X,'   SLFT2 = ',F8.6,F11.6,'     1.0                 '
     D      /7X,'   BSLT1 = ',F8.6,F11.6,'      .98929             '
     E      /7X,'   BSLT2 = ',F8.6,F11.6,'      .95609             '
     F      /7X,'   DUST1 = ',F8.6,F11.6,'      .91995             '
     G      /7X,'   DUST2 = ',F8.6,F11.6,'      .78495             '
     H      /7X,'   DUST3 = ',F8.6,F11.6,'      .63576             '
     I      /7X,'   CARB1 = ',F8.6,F11.6,'      .31482             '
     J      /7X,'   CARB2 = ',F8.6,F11.6,'      .47513             '
     K      )

      WRITE(KW,6019)
 6019 FORMAT(/'  GHGAS',9X,'PPMVK0    PPMVDF    PPGRAD')
      WRITE(KW,6020) (ghg(I),PPMVK0(I),PPMVDF(I),PPGRAD(I),I=1,12)
 6020 FORMAT(1X,a6,' ',F15.7,F10.5,F10.5)
      END IF
      GO TO 9999
C
C-------------
  100 CONTINUE
C-------------
C
      NPAGE=1
      IF(INDEX < 11) NPAGE=KPAGE
      WRITE(KW,6101)
      WRITE(KW,6102)
      FACTOR=1D0/((PLB(L1)-PLB(L1+1))*ppmv_to_cm_at_stp)
      PPMCO2=ULGAS(L1,2)*FACTOR
      PPMO2 =ULGAS(L1,4)*FACTOR
      PPMN2O=ULGAS(L1,6)*FACTOR
      PPMCH4=ULGAS(L1,7)*FACTOR
      PPMF11=ULGAS(L1,8)*FACTOR
      PPMF12=ULGAS(L1,9)*FACTOR
      PPMY11=ULGAS(L1,11)*FACTOR
      PPMZ12=ULGAS(L1,12)*FACTOR
      WRITE(KW,6103) (FULGAS(I),I=1,3),(FULGAS(I),I=6,9)
     +              ,FULGAS(11),FULGAS(12),    (FGOLDH(I),I=1,5)
C     IF(KGASSR > 0)
C    +WRITE(KW,6104) (FULGAS(I+9),I=1,2),(FULGAS(I+9),I=4,9)
C    +              ,FULGAS(11),FULGAS(12),   (FGOLDH(I+9),I=1,5)
      WRITE(KW,6105) PPMCO2,PPMN2O,PPMCH4,PPMF11,PPMF12,PPMY11,PPMZ12
     +             ,(FSTOPX(I),I=1,4),PPMV80(2),(PPMV80(I),I=6,9)
     +             ,(PPMV80(I),I=11,12),KTREND,JYEAR,JDAY,LASTVC
      WRITE(KW,6106) TAUWC0,FCLDTR,EOCTRA,ZOCSRA,KZSNOW,KCLDEM,NTRACE
     +             ,FSAAER,FTTAER,KCLDEP,MADO3M,L1
      WRITE(KW,6107) TAUIC0,FCLDSR,ESNTRA,ZSNSRA,WETTRA,KSIALB,ITR(1)
     +             ,ITR(5),FSBAER,FTBAER,KEEPAL,NL
      WRITE(KW,6108)       FRAYLE,EICTRA,ZICSRA,WETSRA,KCNORM,ITR(2)
     +             ,ITR(6),FSAAER,FTAAER,       KEEP10,MLAT46
      WRITE(KW,6109) TLGRAD,ECLTRA,EDSTRA,ZDSSRA,KANORM,KPGRAD,ITR(3)
     +             ,ITR(7),FSDAER,FTDAER,KWVCON,ICE012,MLON72
      WRITE(KW,6110) PTLISO       ,EVGTRA,ZVGSRA,KEEPRH,KLATZ0,ITR(4)
     +             ,ITR(8),FSVAER,FTVAER,KSOLAR,NORMS0
C
 6101 FORMAT(' (1)FUL:  1',7X,'2',8X,'3',7X,'6',7X,'7',8X,'8',8X,'9'
     +      ,8X,'11',7X,'12',4X,'RADPAR 1/F: (Control/Default'
     +      ,'/Scaling Parameters)')
 6102 FORMAT(4X,'GAS: ','H2O',5X,'CO2',7X,'O3'
     +      ,5X,'N2O',5X,'CH4',5X,'CFC-11',3X,'CFC-12'
     +      ,3X,'CFY-11',3X,'CFZ-12'
     +      ,2X,'Aerosol Global    Ocean     Land  Desert    Haze')
 6103 FORMAT(1X,'FULGAS=',F5.3,F10.5,F7.3,F9.5,F8.5,4F9.5
     +      ,2X,'FGOLDH=',F7.5,2F9.6,2F8.5)
C6104 FORMAT('+',T84,'T'
C    +      /1X,'FULGAS=',1P,1E7.1,1P,2E8.1,1P,2E8.1,1P,4E9.1
C    +      ,' S','FGOLDH=',1P,1E7.1,1P,2E9.2,1P,2E8.1)
 6105 FORMAT(1X,'PPM(1)=(now)',2X,F8.3,8X,F8.5,F8.5,4(1X,F8.7)
     +      ,2X,'TRACER=',F7.5,2F9.6,F8.5
     +      /' PPMV80=(ref)=',0P,F9.3,8X,2F8.5,4(1X,F8.7),2X
     +      ,'KTREND=',I1,2X,'JYEAR=',I4,' JDAY=',I3,5X,'LASTVC=',I7)
 6106 FORMAT(1X,'TAUWC0=',1P,E6.0,' FCLDTR=',0P,F4.2,' EOCTRA=',F3.1
     +      ,1X,'ZOCSRA=',   F3.1,' KZSNOW=',     I4,' KCLDEM=',  I3
     +      ,1X,'NTRACE=',    I3,2X,'FSTAER=',  F3.1,' FTTAER=',F3.1
     +      ,1X,'KCLDEP=',    I1,1X,'MADO3M=',    I2, '   L1=',  I3)
 6107 FORMAT(1X,'TAUIC0=',1P,E6.0,' FCLDSR=',0P,F4.2,' ESNTRA=',F3.1
     +      ,1X,'ZSNSRA=',  F3.1,1X,'WETTRA=',  F4.2,' KSIALB=',  I3
     +      ,1X,'ITR(1)=',    2I2,1X,'FSBAER=', F3.1,' FTBAER=',F3.1
     +      ,1X,'KEEPAL=',    I1,1X,'       ',   ' ','    NL=',  I3)
 6108 FORMAT(1X,'       ',  6X ,  ' FRAYLE=',0P,F4.1,' EICTRA=',F3.1
     +      ,1X,'ZICSRA=',  F3.1,1X,'WETSRA=',  F4.2,' KCNORM=',  I3
     +      ,1X,'ITR(2)=',    2I2,1X,'FSAAER=', F3.1,' FTAAER=',F3.1
     +      ,1X,'KEEP10=',    I1,1X,'       ',   ' ',' MLAT46=',  I2)
 6109 FORMAT(1X,'TLGRAD=',   F6.2,' ECLTRA=',0P,F4.2,' EDSTRA=',F3.1
     +      ,1X,'ZDSSRA=',  F3.1,1X,'KANORM=',  I4  ,' KPGRAD=',  I3
     +      ,1X,'ITR(3)=',    2I2,1X,'FSDAER=', F3.1,' FTDAER=',F3.1
     +      ,1X,'KWVCON=',    I1,1X,'ICE012=',    I1,' MLON72=',  I2)
 6110 FORMAT(1X,'PTLISO=',  F6.1,1X,'           '   ,' EVGTRA=',F3.1
     +      ,1X,'ZVGSRA=',  F3.1,1X,'KEEPRH=',  I4  ,' KLATZ0=',  I3
     +      ,1X,'ITR(4)=',    2I2,1X,'FSVAER=', F3.1,' FTVAER=',F3.1
     +      ,1X,'KSOLAR=',    I1,1X,'NORMS0=',    I1,'        ')
      GO TO 9999
C
C-------------
  200 CONTINUE
C-------------
C
      NPAGE=0
      LUXGAS=0
      IF(INDEX < 11) NPAGE=KPAGE
      WRITE(KW,6201) AUXGAS(LUXGAS+1),S00WM2,S0,COSZ
      DO 202 K=1,9
      DO 201 L=L1,NL
      UXGAS(L,K)=ULGAS(L,K)
  201 CONTINUE
  202 CONTINUE
      IF(LUXGAS < 2) GO TO 205
      LGS=(LUXGAS-2)*9
      DO 203 L=L1,NL
      UXGAS(L,1)=U0GAS(L,1)*FULGAS(1+LGS)
      UXGAS(L,3)=U0GAS(L,3)*FULGAS(3+LGS)
  203 UXGAS(L,5)=U0GAS(L,5)*FULGAS(5+LGS)
C
      DO 204 L=L1,NL
      UXGAS(L,2)=U0GAS(L,2)*FULGAS(2+LGS)
      UXGAS(L,4)=U0GAS(L,4)*FULGAS(4+LGS)
      UXGAS(L,6)=U0GAS(L,6)*FULGAS(6+LGS)
      UXGAS(L,7)=U0GAS(L,7)*FULGAS(7+LGS)
      UXGAS(L,8)=U0GAS(L,8)*FULGAS(8+LGS)
  204 UXGAS(L,9)=U0GAS(L,9)*FULGAS(9+LGS)
  205 CONTINUE
      DO 206 L=NL,L1,-1
      EPS=CLDEPS(L)
      TAER=WSREXT(L,6)
      IPI0=WSRPI0(L,6)*1000.D0+1.D-05
      HLM=0.5D0*(HLB0(L+1)+HLB0(L))
      TLAPS=(TLT(L)-TLB(L))/max(1d-3,HLB0(L+1)-HLB0(L))
      IRHL=RHL(L)*100.0
      IF(PL(L) < 1.D0) THEN
      WRITE(KW,6212) L,PL(L),HLM,TLM(L),TLAPS,SHL(L),IRHL
     +       ,(UXGAS(L,K),K=1,3),(UXGAS(L,K),K=6,9),UXGAS(L,5)
     +       ,SIZEWC(L),SIZEIC(L),FTAUC*TAUWC(L),FTAUC*TAUIC(L),EPS,TAER
     *       ,IPI0
      ELSE
      IF(UXGAS(L,1) >= 1.D0) THEN
      WRITE(KW,6202) L,PL(L),HLM,TLM(L),TLAPS,SHL(L),IRHL
     +       ,(UXGAS(L,K),K=1,3),(UXGAS(L,K),K=6,9),UXGAS(L,5)
     +       ,SIZEWC(L),SIZEIC(L),FTAUC*TAUWC(L),FTAUC*TAUIC(L),EPS,TAER
     *       ,IPI0
      ELSE
      WRITE(KW,6211) L,PL(L),HLM,TLM(L),TLAPS,SHL(L),IRHL
     +       ,(UXGAS(L,K),K=1,3),(UXGAS(L,K),K=6,9),UXGAS(L,5)
     +       ,SIZEWC(L),SIZEIC(L),FTAUC*TAUWC(L),FTAUC*TAUIC(L),EPS,TAER
     *       ,IPI0
      ENDIF
      ENDIF
  206 CONTINUE
      DO 207 I=1,16
  207 SUM0(I)=0.
      DO 210 L=L1,NL
      DO 208 I=1,9
      SUM0(I)=SUM0(I)+ULGAS(L,I)
  208 CONTINUE
      DO 209 I=1,4
      SUM0(12+I)=SUM0(12+I)+TRACER(L,I)*
     *  1d3*.75d0/DENAER(ITR(I))*Q55DRY(ITR(I))/TRRDRY(I)
  209 CONTINUE
      SUM0(10)=SUM0(10)+FTAUC*TAUWC(L)
      SUM0(11)=SUM0(11)+FTAUC*TAUIC(L)
  210 CONTINUE
      TAU55=0.0
      DO 211 L=L1,NL
      TAU55=TAU55+WSREXT(L,6)
  211 CONTINUE
      SUM0(12)=TAU55
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      WRITE(KW,6203) (SUM0(I),I=1,3),(SUM0(I),I=6,9),SUM0(5)
     +               ,SUM0(10),SUM0(11),SUM0(12)
      WRITE(KW,6204) POCEAN,TGO,PLAKE,zlake,SUM0(13),JYEAR
     +             ,BXA(4:5),LASTVC
      WRITE(KW,6205) PEARTH,TGE,SNOWD,ZSNWOI,SUM0(14),JDAY,BXA(6:7)
      WRITE(KW,6206) POICE,TGOI,SNOWOI,ZOICE,SUM0(15),JLAT
     +             ,(SRBALB(I),I=1,6)
      WRITE(KW,6207) PLICE,TGLI,SNOWLI,zmp,SUM0(16),ILON
     +             ,(SRXALB(I),I=1,6)
      PSUM=POCEAN+PEARTH+POICE+PLICE
      snotyp='DRY' ; if(flags) snotyp='WET'
      WRITE(KW,6208) TGMEAN,snotyp,fmp
     +               ,PSUM,TSL,WMAG,LS1_loc,(PVT(I),I=1,11)
      write(kw,6213) snow_frac(1),snow_frac(2),agesn(1),
     +  agesn(2),agesn(3),wearth,fulgas(4),fulgas(5),fulgas(10)
 6213 FORMAT(1X,'FSNWds=',F6.4,' FSNWvg=',F6.4,'  AGESN=[EA:',F6.3,
     +      ' OI:',F6.3,' LI:',F6.3,'] WEARTH=',F6.4,1X,
     +      ' FULGAS[ 4=O2:',F3.1,' 5=NO2:',F3.1,' 10=N2C:',F3.1,']')
      WRITE(KW,6209) (PRNB(1:2,I),PRNX(1:2,I),I=1,4),BXA(1:3)
      WRITE(KW,6210)
 6201 FORMAT(' (2) RADPAR G/L: (Input Data)'
     +      ,2X,'Absorber Amount per Layer:'
     +      ,'  U',1A1,'GAS(L,K) in cm**3(STP)/cm**2',2X,'S00WM2=',F9.4
     +      ,1X,'S0=',F9.4,2X,'COSZ=',F6.4
     +      /' LN     PL  HLM   TLM  TLAP   SHL  .RH     '
     +      ,'H2O   CO2    O3    N2O    CH4  CFC-11'
     +      ,'  CFC-12   NO2   WC.SIZ.IC  WC.TAU.IC CLEP  A TAU PI0')
 6202 FORMAT(1X,I2,F7.2,F5.1,F7.2,F5.1,1X,F7.6,I3
     +      ,F8.2,F6.2,1X,F6.5,1X,F5.4
     +      ,F7.4,1P,3E8.1,0P,2F5.1,F6.2,F5.2,1X,F4.3,F6.3,I5)
 6203 FORMAT(24X,' Column Amount',F7.1,F7.2,1X,F6.5
     +       ,1X,F5.4,F7.4,1P,3E8.1,0P,10X,F6.2,F5.2,5X,F6.3)
 6204 FORMAT( 1X,'PWATER=',F6.4,'    TGO=' ,F6.2,1X,' PLAKE=',F6.3
     +      , 1X,' ZLAKE=',F6.3,' TRACER 1=',F5.3,' JYEAR=',I4
     +      , 3X,'BSNVIS=',F6.4,' BSNNIR=' ,F6.4,7X,'LASTVC=',I7)
 6205 FORMAT(    ' PEARTH=',F6.4,'    TGE=',F6.2,'  SNOWD=',2F6.3
     +      ,    '  ZSNOW=',F6.3,'  Sums: 2=',F5.3
     +      ,     '  JDAY=',I4  ,2X,' XSNVIS=',F6.4,' XSNNIR=',F6.4
     +      , 8X,'NIRALB VISALB')
 6206 FORMAT(    '  POICE=',F6.4,'   TGOI=',F6.2,' SNOWOI=',F6.3
     +      ,    '  ZOICE=',F6.3,'        3=',F5.3
     +      ,     '  JLAT=',I4,  2X,' SRBALB=',F6.4
     +      ,4F7.4,F7.4)
 6207 FORMAT(    '  PLICE=',F6.4,'   TGLI=',F6.2,' SNOWLI=',F6.3
     +      ,    '  ZMLTP=',F6.3,'        4=',F5.3
     +      ,     '  ILON=',I4,  2X,' SRXALB=',F6.4
     +      ,4F7.4,F7.4)
 6208 FORMAT(8X,6('-'),' TGMEAN=',F6.2,'    SNOW : ',a3,'  FMLTP='
     +      ,F6.3,'  BSAND TUNDRA GRASSL SHRUBS  TREES DECIDF'
     +       ,' EVERGF','  RAINF','  CROPS','  BDIRT','  ALGAE'
     +      /    '   PSUM=',F6.4,'    TSL=',F6.2,' WINDSP=',F6.3
     +       ,'  LS1L=',I2,T54,'PVT=',F6.4,10F7.4)
 6209 FORMAT(' BOCVIS BOCNIR XOCVIS XOCNIR BEAVIS BEANIR XEAVIS XEANIR'
     +      ,' BOIVIS BOINIR XOIVIS XOINIR BLIVIS BLINIR XLIVIS XLINIR'
     +      ,' EXPSNE EXPSNO EXPSNL'/1X,F6.4,18F7.4)
 6210 FORMAT(' ')
 6211 FORMAT(1X,I2,F7.2,F5.1,F7.2,F5.1,1X,F7.6,I3
     +      ,F8.5,F6.2,1X,F6.5,1X,F5.4
     +      ,F7.4,1P,3E8.1,0P,2F5.1,F6.2,F5.2,1X,F4.3,F6.3,I5)
 6212 FORMAT(1X,I2,F7.4,F5.1,F7.2,F5.1,1X,F7.6,I3
     +      ,F8.5,F6.2,1X,F6.5,1X,F5.4
     +      ,F7.4,1P,3E8.1,0P,2F5.1,F6.2,F5.2,1X,F4.3,F6.3,I5)
C
      GO TO 9999
C
C-------------
  300 CONTINUE
C-------------
C
      NPAGE=0
      IF(INDEX < 11) NPAGE=KPAGE
      IF(NL > 13) NPAGE=1
      L=NL+1
      SRALB =SRUFLB(L)/(SRDFLB(L)+1.E-10)
      STNFLB=SRNFLB(L)-TRNFLB(L)
      WRITE(KW,6301) NORMS0
      WRITE(KW,6302) L,PLB(L),HLB0(L),TLT(L-1) ! TLB(LN+1) unused/set
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L)
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),STNFLB,SRALB
      DO 301 L=NL,L1,-1
      CRHRF=8.4167/(PLB(L)-PLB(L+1))
      STNFLB=SRNFLB(L)-TRNFLB(L)
      STFHR =SRFHRL(L)-TRFCRL(L)
      TRDCR =TRFCRL(L)*CRHRF
      SRDHR =SRFHRL(L)*CRHRF
      STDHR=STFHR*CRHRF
      SRALB =SRUFLB(L)/(SRDFLB(L)+1.E-10)
!eq   SRXVIS=SRXATM(1)
!eq   SRXNIR=SRXATM(2)
      IF(PLB(L) < 1.D0) THEN
      WRITE(KW,6313) L,PLB(L),HLB0(L),TLB(L),TLT(L)
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L),TRFCRL(L)
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),SRFHRL(L)
     +             ,STNFLB,STFHR,STDHR,TRDCR,SRDHR,SRALB
      ELSE
      WRITE(KW,6303) L,PLB(L),HLB0(L),TLB(L),TLT(L)
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L),TRFCRL(L)
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),SRFHRL(L)
     +             ,STNFLB,STFHR,STDHR,TRDCR,SRDHR,SRALB
      ENDIF
  301 CONTINUE
C
      DO 302 II=1,3
      PFW=TRDFLB(L1)
      IF(II==2) PFW=TRUFLB(L1)
      IF(II==3) PFW=TRUFLB(NL+1)
      IPF=PFW
      DPF=PFW-IPF
      IF(IPF < 1) IPF=1
      IF(IPF > 899) IPF=899
  302 TKEFF(II)=TKPFT(IPF)+DPF*(TKPFT(IPF+1)-TKPFT(IPF))
C
      WRITE(KW,6304) WINDZF(1),WINDZT(1),TOTLZF(1),TOTLZT(1),
     +              (FSRNFG(I),I=1,4),LTOPCL,JLAT,JYEAR
      WRITE(KW,6305) WINDZF(2),WINDZT(2),TOTLZF(2),TOTLZT(2),
     +              (FTRUFG(I),I=1,4),LBOTCL,ILON,JDAY
      IF(KORDER==0) WRITE(KW,6306) WINDZF(3),WINDZT(3),TOTLZF(3)
     +              ,TOTLZT(3),(I,I=1,16)
      IF(KORDER==1) WRITE(KW,6307) WINDZF(3),WINDZT(3),TOTLZF(3)
     +              ,TOTLZT(3),(I,I=1,16)
      FRACSL=0.D0
      IF(KORDER==0) WRITE(KW,6308) TKEFF(1),TKEFF(2),TKEFF(3)
     +               ,(SRKALB(NORDER(I)),I=1,16),BTEMPW,TRUFTW,SRIVIS
     +               ,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR
      IF(KORDER==1) WRITE(KW,6308) TKEFF(1),TKEFF(2),TKEFF(3)
     +               ,(SRKALB(I),I=1,16),BTEMPW,TRUFTW,SRIVIS
     +               ,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR
      WRITE(KW,6309) TRDFGW,TRUFGW,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR
     +             ,ALBNIR
      WRITE(KW,6310) SRXVIS,SRXNIR,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR
     +             ,SRANIR
C
C
 6301 FORMAT(/' (3) RADPAR M/S: (Output Data)'
     +      ,T37,'Thermal Fluxes (W/M**2)',4X,'Solar Fluxes (W/M**2)'
     +      ,1X,'NORMS0=',I1,'  Energy  Input  Heat/Cool Deg/Day Alb'
     +      ,'do'/' LN     PLB   HLB    TLB    TLT '
     +      ,'  TRDFLB TRUFLB TRNFLB TRFCRL   SRDFLB  SRUFLB  SRNFLB'
     +      ,' SRFHRL  STNFLB  STFHR  SR-TR TR=CR SR=HR SRALB')
 6302 FORMAT(1X,I2,F9.3,F6.2,1X,F6.2,8X,3F7.2,8X,3F8.2,7X,F8.2,26X,F6.4)
 6303 FORMAT(1X,I2,F9.3,F6.2,2F7.2,1X,3F7.2,F7.2,1X,3F8.2,F7.2,1X,F7.2
     +      ,1X,F6.2,1X,3F6.2,1X,F5.4)
 6304 FORMAT( 1X,'XMU  WINDZF WINDZT   TOTLZF TOTLZT'/
     +        1X,'1.0',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'FR.SRNLB1'
     +      ,' OCEAN=',F7.2,' EARTH=',F7.2,'  OICE=',F7.2,'   LICE='
     +      ,F7.2,1X,' LTOPCL=',I2,' JLAT=',I2,' JYEAR=',I4)
 6305 FORMAT( 1X,'0.5',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'FR.TRULB1'
     +      ,' OCEAN=',F7.4,' EARTH=',F7.4,'  OICE=',F7.4,'   LICE='
     +      ,F7.4,1X,' LBOTCL=',I2,' ILON=',I2,'  JDAY=',I4)
 6306 FORMAT( 1X,'0.1',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'L=',I3,15I6)
 6307 FORMAT( 1X,'0.1',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'K=',I3,15I6)
 6308 FORMAT(' TKeff= ',F6.2,2F7.2,'  SRKALB=',16F6.4/
     +        1X,'At Top of Atm: ',' BTEMPW=',F6.2,1X,' TRUFTW=',F6.3
     +      , 2X,' SRIVIS=',F6.2,' SROVIS=',F6.2,   ' PLAVIS=',F6.4
     +      , 2X,' SRINIR=',F6.2,' SRONIR=',F6.2,   ' PLANIR=',F6.4)
 6309 FORMAT( 1X,'At Bot of Atm: ',' TRDFGW=',F6.3,1X,' TRUFGW=',F6.3
     +      , 2X,' SRDVIS=',F6.2,' SRUVIS=',F6.2,   ' ALBVIS=',F6.4
     +      , 2X,' SRDNIR=',F6.2,' SRUNIR=',F6.2,   ' ALBNIR=',F6.4)
 6310 FORMAT( 1X,'In Atmosphere: ',' SRXVIS=',F6.4,1X,' SRXNIR=',F6.4
     +      , 2X,' SRTVIS=',F6.4,' SRRVIS=',F6.4,   ' SRAVIS=',F6.4
     +      , 2X,' SRTNIR=',F6.4,' SRRNIR=',F6.4,   ' SRANIR=',F6.4)
 6311 FORMAT(' ')
 6313 FORMAT(1X,I2,F9.5,F6.2,2F7.2,1X,F7.4,2F7.2,F7.4,1X,3F8.2,F7.4
     +      ,1X,F7.2,F7.4,1X,3F6.2,1X,F5.4)
      GO TO 9999
C
C-------------
  400 CONTINUE
C-------------
C
C                                (4A)  Total Aerosol Qx, Qs, g, Pi0
C                                ----------------------------------
      NPAGE=1
      IF(INDEX < 11) NPAGE=KPAGE
      WRITE(KW,6401)
      DO 402 K=1,6
      SUM1(K)=0.
      SUM2(K)=0.
      SUM3(K)=0.
      DO 401 L=L1,NL
      SUM1(K)=SUM1(K)+WSREXT(L,K)
      SUM2(K)=SUM2(K)+WSRSCT(L,K)
      SUM3(K)=SUM3(K)+WSRSCT(L,K)*WSRGCB(L,K)
  401 CONTINUE
      SUM3(K)=SUM3(K)/(SUM2(K)+1.D-10)
      SUM0(K)=SUM2(K)/(SUM1(K)+1.D-10)
  402 CONTINUE
      WRITE(KW,6402) (K,K=1,6),(K,K=1,6)
      DO 403 L=NL,L1,-1
      WRITE(KW,6403) L,PLB(L),HLB0(L)
     +              ,(WSREXT(L,J),J=1,6),(WSRSCT(L,J),J=1,6)
  403 CONTINUE
      WRITE(KW,6404) (SUM1(K),K=1,6),(SUM2(K),K=1,6)
      NPAGE=0
      IF(NL > 13) NPAGE=1
      WRITE(KW,6405) KANORM
      WRITE(KW,6406) (K,K=1,6),(K,K=1,6)
      DO 404 L=NL,L1,-1
      WRITE(KW,6407) L,PL(L),DPL(L)
     +              ,(WSRGCB(L,J),J=1,6),(WSRPI0(L,J),J=1,6)
  404 CONTINUE
      WRITE(KW,6408) (SUM3(K),K=1,6),(SUM0(K),K=1,6)
C     WRITE(KW,6420) (SRBALB(K),K=1,6)
C     WRITE(KW,6421) (SRXALB(K),K=1,6)
C     WRITE(KW,6422)
      SUMT=0.
      DO 406 J=1,5
      TAU55=0.
      DO 405 I=1,11  !  NAERO
  405 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)
      WRITE(KW,6423) J,FGOLDH(J),TAU55
  406 SUMT=SUMT+TAU55
      WRITE(KW,6438) SUMT
      WRITE(KW,6424) BCOLX,ACOLX,DCOLX,VCOLX,TCOLX
      DO 407 I=1,8
      WRITE(KW,6425)
  407 CONTINUE
C
C                                (4B)  Water/Ice Cloud Qx, Qs, g, Pi0
C                                ------------------------------------
      NPAGE=1
      IF(INDEX < 11) NPAGE=KPAGE
      WRITE(KW,6411)
      DO 412 K=1,6
      SUM1(K)=0.
      SUM2(K)=0.
      SUM3(K)=0.
      DO 411 L=L1,NL
      SUM1(K)=SUM1(K)+SRCEXT(L,K)
      SUM2(K)=SUM2(K)+SRCSCT(L,K)
      SUM3(K)=SUM3(K)+SRCSCT(L,K)*SRCGCB(L,K)
  411 SRCPI0(L,K)=SRCSCT(L,K)/(SRCEXT(L,K)+1.D-10)
      SUM3(K)=SUM3(K)/(SUM2(K)+1.D-10)
      SUM0(K)=SUM2(K)/(SUM1(K)+1.D-10)
  412 CONTINUE
      WRITE(KW,6412) (K,K=1,6),(K,K=1,6)
      DO 413 L=NL,L1,-1
      WRITE(KW,6413) L,PLB(L),HLB0(L)
     +              ,(SRCEXT(L,J),J=1,6),(SRCSCT(L,J),J=1,6)
  413 CONTINUE
      WRITE(KW,6414) (SUM1(K),K=1,6),(SUM2(K),K=1,6)
      NPAGE=0
      IF(NL > 13) NPAGE=1
      WRITE(KW,6415) KANORM
      WRITE(KW,6416) (K,K=1,6),(K,K=1,6)
      DO 414 L=NL,L1,-1
      WRITE(KW,6417) L,PL(L),DPL(L)
     +              ,(SRCGCB(L,J),J=1,6),(SRCPI0(L,J),J=1,6)
  414 CONTINUE
      WRITE(KW,6418) (SUM3(K),K=1,6),(SUM0(K),K=1,6)
      WRITE(KW,6420) (SRBALB(K),K=1,6)
      WRITE(KW,6421) (SRXALB(K),K=1,6)
      WRITE(KW,6422)
      SUMT=0.
      DO 416 J=1,5
      TAU55=0.
      DO 415 I=1,11  !  NAERO
  415 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)
      WRITE(KW,6423) J,FGOLDH(J),TAU55
  416 SUMT=SUMT+TAU55
      WRITE(KW,6438) SUMT
      DO 417 I=1,2
      WRITE(KW,6425)
  417 CONTINUE
C
C                                (4C)  Aerosol + Cloud Qx, Qs, g, Pi0
C                                ------------------------------------
      NPAGE=1
      IF(INDEX < 11) NPAGE=KPAGE
      WRITE(KW,6426)
      DO 419 K=1,6
      SUM1(K)=0.
      SUM2(K)=0.
      SUM3(K)=0.
      DO 418 L=L1,NL
      SUM1(K)=SUM1(K)+DBLEXT(L,K)
      SUM2(K)=SUM2(K)+DBLSCT(L,K)
      SUM3(K)=SUM3(K)+DBLSCT(L,K)*DBLGCB(L,K)
  418 DBLPI0(L,K)=DBLSCT(L,K)/(DBLEXT(L,K)+1.E-10)
      SUM3(K)=SUM3(K)/(SUM2(K)+1.E-10)
      SUM0(K)=SUM2(K)/(SUM1(K)+1.E-10)
  419 CONTINUE
      WRITE(KW,6427) (K,K=1,6),(K,K=1,6)
      DO 420 L=NL,L1,-1
      WRITE(KW,6428) L,PLB(L),HLB0(L)
     +              ,(DBLEXT(L,J),J=1,6),(DBLSCT(L,J),J=1,6)
  420 CONTINUE
      WRITE(KW,6429) (SUM1(K),K=1,6),(SUM2(K),K=1,6)
      NPAGE=0
      IF(NL > 13) NPAGE=1
      WRITE(KW,6430) KANORM
      WRITE(KW,6431) (K,K=1,6),(K,K=1,6)
      DO 421 L=NL,L1,-1
      WRITE(KW,6432) L,PL(L),DPL(L)
     +              ,(DBLGCB(L,J),J=1,6),(DBLPI0(L,J),J=1,6)
  421 CONTINUE
      WRITE(KW,6433) (SUM3(K),K=1,6),(SUM0(K),K=1,6)
      WRITE(KW,6434) (SRBALB(K),K=1,6)
      WRITE(KW,6435) (SRXALB(K),K=1,6)
      WRITE(KW,6436)
      SUMT=0.
      DO 423 J=1,5
      TAU55=0.
      DO 422 I=1,11  !  NAERO
  422 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)
      WRITE(KW,6437) J,FGOLDH(J),TAU55
  423 SUMT=SUMT+TAU55
      WRITE(KW,6438) SUMT
      DO 424 I=1,2
      WRITE(KW,6439)
  424 CONTINUE
C
C                                (4D)  11-Comp Aerosol Qx, Qs, g, Pi0
C                                ------------------------------------
      NPAGE=1
      IF(INDEX < 11) NPAGE=KPAGE
      WRITE(KW,6440)  KWTRAB,(N,N=1,11)
      WRITE(KW,6441)
      DO 425 K=1,6
      WRITE(KW,6442) K,(SRAQEX(K,N),N=1,11)
  425 CONTINUE
      WRITE(KW,6443)
      DO 426 K=1,6
      WRITE(KW,6442) K,(SRAQSC(K,N),N=1,11)
  426 CONTINUE
      WRITE(KW,6444)
      DO 427 K=1,6
      WRITE(KW,6442) K,(SRAQCB(K,N),N=1,11)
  427 CONTINUE
      WRITE(KW,6445) TRABCD(1),TRAXSG(KWTRAB+1)
      DO 428 K=1,33
      IF(KWTRAB==0) WRITE(KW,6442) K,(TRAQAB(K,N),N=1,11)
      IF(KWTRAB==1) WRITE(KW,6442) K,(TRAQEX(K,N),N=1,11)
      IF(KWTRAB==2) WRITE(KW,6442) K,(TRAQSC(K,N),N=1,11)
      IF(KWTRAB==3) WRITE(KW,6442) K,(TRAQCB(K,N),N=1,11)
      IF(KWTRAB==4) THEN
      DO N=1,11
      TRPI0K(N)=TRAQSC(K,N)/(1.D-10+TRAQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,11)
      ENDIF
  428 CONTINUE
      DO 429 I=1,1
      WRITE(KW,6446)
  429 CONTINUE
C
C
C                                (4E)  10-Comp Aerosol Qx, Qs, g, Pi0
C                                ------------------------------------
      WRITE(KW,6450) KWTRAB,(N,N=1, 6),(REFDRY(N),N=1, 6)
      WRITE(KW,6451)
      DO 435 K=1,6
      WRITE(KW,6452) K,(SRHQEX(K,1,N),N=1, 4),(SRBQEX(K,N),N=5, 6)
  435 CONTINUE
      WRITE(KW,6453)
      DO 436 K=1,6
      WRITE(KW,6452) K,(SRHQSC(K,1,N),N=1, 4),(SRBQSC(K,N),N=5, 6)
  436 CONTINUE
      WRITE(KW,6454)
      DO 437 K=1,6
      WRITE(KW,6452) K,(SRHQCB(K,1,N),N=1, 4),(SRBQCB(K,N),N=5, 6)
  437 CONTINUE
      WRITE(KW,6455) TRABCD(2),TRAXSG(1) !obs TRAXSG(KWTRAB+1)
      DO 438 K=1,33
      IF(KWTRAB==0) WRITE(KW,6442) K,(TRHQAB(K,1,N),N=1, 4),
     *                                 (TRBQAB(K,N),N=5, 6)
!obs  IF(KWTRAB==1) WRITE(KW,6442) K,(TRHQEX(K,1,N),N=1, 4),
!obs *                                 (TRBQEX(K,N),N=5, 6)
!obs  IF(KWTRAB==2) WRITE(KW,6442) K,(TRHQSC(K,1,N),N=1, 4),
!obs *                                 (TRBQSC(K,N),N=5, 6)
!obs  IF(KWTRAB==3) WRITE(KW,6442) K,(TRHQCB(K,1,N),N=1, 4),
!obs *                                 (TRBQCB(K,N),N=5, 6)
!obs  IF(KWTRAB==4) THEN
!obs  DO N=1,4
!obs  TRPI0K(N)=TRHQCB(K,1,N)/(1.D-10+TRHQEX(K,1,N))
!obs  END DO
!obs  DO N=5,6 ! 10
!obs  TRPI0K(N)=TRBQSC(K,N)/(1.D-10+TRBQEX(K,N))
!obs  END DO
!obs  WRITE(KW,6442) K,(TRPI0K(N),N=1, 6)
!obs  ENDIF
  438 CONTINUE
      DO 439 I=1,1
      WRITE(KW,6456)
  439 CONTINUE
C
C                             (4F  8-size Dust Aerosol Qx, Qs, g, Pi0
C                             ---------------------------------------
      NPAGE=1
      IF(INDEX < 11) NPAGE=KPAGE
      WRITE(KW,6460) KWTRAB,(N,N=1,8),(REDUST(N),N=1,8)
      WRITE(KW,6461)
      DO 445 K=1,6
      WRITE(KW,6462) K,(SRAQEX(K,N),N=1,8)
  445 CONTINUE
      WRITE(KW,6463)
      DO 446 K=1,6
      WRITE(KW,6462) K,(SRAQSC(K,N),N=1,8)
  446 CONTINUE
      WRITE(KW,6464)
      DO 447 K=1,6
      WRITE(KW,6462) K,(SRAQCB(K,N),N=1,8)
  447 CONTINUE
      WRITE(KW,6465) TRABCD(4),TRAXSG(KWTRAB+1)
      DO 448 K=1,33
      IF(KWTRAB==0) WRITE(KW,6442) K,(TRDQAB(K,N),N=1,8)
      IF(KWTRAB==1) WRITE(KW,6442) K,(TRDQEX(K,N),N=1,8)
      IF(KWTRAB==2) WRITE(KW,6442) K,(TRDQSC(K,N),N=1,8)
      IF(KWTRAB==3) WRITE(KW,6442) K,(TRDQCB(K,N),N=1,8)
      IF(KWTRAB==4) THEN
      DO N=1,8
      TRPI0K(N)=TRDQSC(K,N)/(1.D-10+TRDQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,8)
      ENDIF
  448 CONTINUE
      DO 449 I=1,1
      WRITE(KW,6466)
  449 CONTINUE
C
C                           (4G  15-Size/phase Cloud Qx, Qs, g, Pi0
C                           ---------------------------------------
      NPAGE=1
      IF(INDEX < 11) NPAGE=KPAGE
      WRITE(KW,6470) KWTRAB,(N,N=1,15)
      WRITE(KW,6471)
      DO 430 K=1,6
      WRITE(KW,6472) K,(SRCQEX(K,N),N=1,15)
  430 CONTINUE
      WRITE(KW,6473)
      DO 431 K=1,6
      WRITE(KW,6472) K,(SRCQSC(K,N),N=1,15)
  431 CONTINUE
      WRITE(KW,6474)
      DO 432 K=1,6
      WRITE(KW,6472) K,(SRCQCB(K,N),N=1,15)
  432 CONTINUE
      WRITE(KW,6475) TRABCD(3),TRAXSG(KWTRAB+1)
      DO 433 K=1,33
      IF(KWTRAB==0) WRITE(KW,6472) K,(TRCQAB(K,N),N=1,15)
      IF(KWTRAB==1) WRITE(KW,6472) K,(TRCQEX(K,N),N=1,15)
      IF(KWTRAB==2) WRITE(KW,6472) K,(TRCQSC(K,N),N=1,15)
      IF(KWTRAB==3) WRITE(KW,6472) K,(TRCQCB(K,N),N=1,15)
      IF(KWTRAB==4) THEN
      DO N=1,15
      TRPI0K(N)=TRCQSC(K,N)/(1.D-10+TRCQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,15)
      ENDIF
  433 CONTINUE
      DO 434 I=1,2
      WRITE(KW,6476)
  434 CONTINUE
C
 6401 FORMAT(' (4A) Aerosol Input for Solar Radiation:'
     +      ,' Aerosol Radiative Parameters'
     +      ,T81,'LIST: SRAEXT(L,K),SRASCT(L,K),SRAGCB(L,K),SRAPI0(L,K)'
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'
     +      ,/T24,53('-'),4X,53('-'))
 6402 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)
 6403 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6404 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)
 6405 FORMAT(6X,'KANORM=',1I1/T48,'COSBAR',T105,'PIZERO'
     +      ,/T24,53('-'),4X,53('-'))
 6406 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)
 6407 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6408 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)
C
 6411 FORMAT(' (4B) Cloud Input for Solar Radiation:'
     +      ,'   Cloud Radiative Parameters'
     +      ,T81,'LIST: SRCEXT(L,K),SRCSCT(L,K),SRCGCB(L,K),SRCPI0(L,K)'
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'
     +      ,/T24,53('-'),4X,53('-'))
 6412 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)
 6413 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6414 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)
 6415 FORMAT(6X,'KANORM=',1I1/T48,'COSBAR',T105,'PIZERO'
     +      ,/T24,53('-'),4X,53('-'))
 6416 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)
 6417 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6418 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)
C
 6420 FORMAT(/1X,T7,'ALBEDO RSURFB=',2X,6F9.6,3X,6F9.6)
 6421 FORMAT( 1X,T7,'ALBEDO RSURFX=',2X,6F9.6,3X,6F9.6)
 6422 FORMAT(///T44,'AEROSOL COMPOSITION AND TYPE MIX:'
     +      ,T81,'FACTOR',6X,'VALUE',T107,'TAU(0.55)'/)
 6423 FORMAT(T81,'FGOLDH(',I1,') =',1P,E9.2,5X,0P,F7.4)
 6424 FORMAT(/T11,'SUM COLUMN TAU(0.55) =    BkGrnd    ClimAer   D Dust'
     +       ,'    VolAer    TotAer'/T33,5F10.5)
 6425 FORMAT(' ')
C
 6426 FORMAT(' (4C) Cloud+Aerosol  Output from SOLARM/SGPGXG:'
     +      ,'   Cloud+Aerosol Rad Parameters'
     +      ,T81,'LIST: DBLEXT(L,K),DBLSCT(L,K),DBLGCB(L,K),DBLPI0(L,K)'
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'
     +      ,/T24,53('-'),4X,53('-'))
 6427 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)
 6428 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6429 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)
 6430 FORMAT(6X,'KANORM=',1I1/T48,'COSBAR',T105,'PIZERO'
     +      ,/T24,53('-'),4X,53('-'))
 6431 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)
 6432 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6433 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)
 6434 FORMAT(/1X,T7,'ALBEDO RSURFB=',2X,6F9.6,3X,6F9.6)
 6435 FORMAT( 1X,T7,'ALBEDO RSURFX=',2X,6F9.6,3X,6F9.6)
 6436 FORMAT(///T44,'AEROSOL COMPOSITION AND TYPE MIX:'
     +      ,T81,'FACTOR',6X,'VALUE',T107,'TAU(0.55)'/)
 6437 FORMAT(T81,'FGOLDH(',I1,') =',1P,E9.2,5X,0P,F7.4)
 6438 FORMAT(/T81,'SUM COLUMN TAU(0.55) =',F10.4)
 6439 FORMAT(' ')
C
 6440 FORMAT(' (4D) Background Aerosol Solar and Thermal Mie '
     +      ,'Scattering Parameters:'
     +      ,T81,'List: SRAQEX(L,K),SRAQST(L,K),SRAQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,11I8/
     +        '   AEROSOL  ACID1   SSALT   SLFT1   SLFT2   BSLT1'
     +        ,'   BSLT2   DUST1   DUST2   DUST3   CARB1   CARB2'/
     +        '   SIZE      0.5     2.0     0.3     1.0     0.5 '
     +        ,'    2.0     0.5     2.0     8.0     0.1     0.5 ')
 6441 FORMAT('  K  SRAQEX')
 6442 FORMAT(I3,6X,15F8.5)
 6443 FORMAT('  K  SRAQSC')
 6444 FORMAT('  K  SRAQCB')
 6445 FORMAT('  K  ',2A3)
 6446 FORMAT(' ')
C
 6450 FORMAT(' (4E) Climatology Aerosol Solar and Thermal Mie '
     +      ,'Scattering Parameters:'
     +      ,T81,'List: SRBQEX(L,K),SRBQST(L,K),SRBQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X, 6I8/
     +        '   AEROSOL   SO4     SEA     ANT     OCX     BCI '
     +        ,'    BCB'/ ! OCN     OCB     BCB     SSB
     +        '   SIZE ', 6F8.1)
 6451 FORMAT('  K  SRBQEX - DRY')
 6452 FORMAT(I3,6X,15F8.5)
 6453 FORMAT('  K  SRBQSC - DRY')
 6454 FORMAT('  K  SRBQCB - DRY')
 6455 FORMAT('  K  ',2A3,' - DRY')
 6456 FORMAT(' ')
C
 6460 FORMAT(' (4F) Desert Dust Aerosol Solar and Thermal Mie '
     +      ,'Scattering Parameters:'
     +      ,T81,'List: SRDQEX(L,K),SRDQST(L,K),SRDQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,8I8/
     +        '   AEROSOL  CLAY1   CLAY2   CLAY3   CLAY4   SILT1'
     +        ,'   SILT2   SILT3   SILT4                        '/
     +        '   SIZE ',8F8.1)
 6461 FORMAT('  K  SRDQEX')
 6462 FORMAT(I3,6X,15F8.5)
 6463 FORMAT('  K  SRDQSC')
 6464 FORMAT('  K  SRDQCB')
 6465 FORMAT('  K  ',2A3)
 6466 FORMAT(' ')
C
 6470 FORMAT(' (4G) Cloud Input for Solar, Thermal Radiation:'
     +      ,'  Mie Cloud Radiative Properties'
     +      ,T81,'List: SRCQEX(L,K),SRCQST(L,K),SRCQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,15I8/
     +        ' WIM CLOUD  WAT05   WAT10   WAT15   WAT20   WAT25'
     +                ,'   ICE05   ICE15   ICE25   ICE50   ICE75'
     +                ,'   MIC05   MIC15   MIC25   MIC50   MIC75')
 6471 FORMAT('  K  SRCQEX')
 6472 FORMAT(I3,6X,15F8.5)
 6473 FORMAT('  K SRCQSC')
 6474 FORMAT('  K SRCQCB')
 6475 FORMAT('  K ',2A3)
 6476 FORMAT(' ')
      GO TO 9999
C
C-------------
  500 CONTINUE
C-------------
C
      NPAGE=1
      IF(INDEX < 11) NPAGE=KPAGE
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      SUMK=0.0
      DO 501 K=1,33
      BGFLUX(K) = PLANCK(ITG,K) - (PLANCK(ITG,K)-PLANCK(ITG+1,K))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
  501 CONTINUE
      LK=0
      DO 503 K=1,33
      TAUSUM(K)=0. !!sl TAUSL(K)
      DO 502 L=L1,NL
      TRTAUK(L,K)=TRGXLK(L,K)+TRCALK(L,K)+TRAALK(L,K)
      TAUSUM(K)=TAUSUM(K)+TRGXLK(L,K)+TRCALK(L,K)+TRAALK(L,K)
  502 CONTINUE
  503 CONTINUE
      WRITE(KW,6501)
      WRITE(KW,6502) (K,K=1,13)
      DO 504 L=NL,L1,-1
      WRITE(KW,6503) L,PL(L),TLM(L),(TRTAUK(L,K),K=1,13)
  504 CONTINUE
!sl   WRITE(KW,6504) (TAUSL(K),K=1,13)
      WRITE(KW,6505) (TAUSUM(K),K=1,13)
      WRITE(KW,6506) SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6507) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
      NPAGE=0
      IF(NL > 13)  NPAGE=1
      WRITE(KW,6508) NPAGE
      WRITE(KW,6509) (K,K=14,33)
      DO 505 L=NL,L1,-1
      WRITE(KW,6510) L,(TRTAUK(L,K),K=14,33)
  505 CONTINUE
!sl   WRITE(KW,6511) ( TAUSL(K),K=14,33)
      WRITE(KW,6512) (TAUSUM(K),K=14,33)
      WRITE(KW,6513) (BGFLUX(K),K=14,33)
      WRITE(KW,6514) (BGFRAC(K),K=14,33)
      DO 506 I=1,10
      WRITE(KW,6515)
  506 CONTINUE
C
 6501 FORMAT(' (5) TAU TABLE FOR THERMAL RADIATION: CONTAINS'
     +      ,' TOTAL SPECIFIED GAS, CLOUD & AEROSOL ABSORPTION'
     +      ,T99,'TRGXLK(L,K),TRCALK(L,K),TRCAAK(L,K)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6502 FORMAT(' LN     PL     TLM      K=',I1,4X,'K=',I2,9I9,3I8)
 6503 FORMAT(1X,I2,F8.3,F7.2,1X,10F9.4,3F8.3)
!sl6504 FORMAT(/4X,'SURFACE LAYER= ',10F9.4,3F8.3)
 6505 FORMAT( 4X,'COLUMN AMOUNT= ',10F9.4,3F8.3)
 6506 FORMAT(/1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6507 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6508 FORMAT(1I1/4X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +      ,T83,'OZONE:',T100,'PRINCIPAL ABSORBER REGION'
     +      /4X,76('-'),2X,50('-'))
 6509 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6510 FORMAT( 1X,  I2,6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
!sl6511 FORMAT(/1X,'SL',6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
 6512 FORMAT( 1X,'CA',5F7.4,1F7.3,3F6.2,2F6.1,1F6.0,4F7.4,2F6.3,2F6.2)
 6513 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6514 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6515 FORMAT(' ')
      GO TO 9999
C
C-------------
  600 CONTINUE
C-------------
C
      NPAGE=1
      IF(INDEX < 11) NPAGE=KPAGE
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      SUMK=0.0
      DO 601 K=1,33
      BGFLUX(K) = PLANCK(ITG,K) - (PLANCK(ITG,K)-PLANCK(ITG+1,K))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
  601 CONTINUE
      WRITE(KW,6601)
      WRITE(KW,6602) (K,K=1,13)
      DO 602 L=NL,L1,-1
      WRITE(KW,6603) L,PL(L),TLM(L),(TRGXLK(L,K),K=1,13)
  602 CONTINUE
      LK=0
      DO 604 K=1,33
      TAUSUM(K)=0. !!sl  TAUSL(K)
      DO 603 L=L1,NL
  603 TAUSUM(K)=TAUSUM(K)+TRGXLK(L,K)
  604 CONTINUE
!sl   WRITE(KW,6604) (TAUSL(K),K=1,13)
      WRITE(KW,6605) (TAUSUM(K),K=1,13)
      WRITE(KW,6606) SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6607) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
      NPAGE=0
      IF(NL > 13)  NPAGE=1
      WRITE(KW,6608)
      WRITE(KW,6609) (K,K=14,33)
      DO 605 L=NL,L1,-1
      WRITE(KW,6610) L,(TRGXLK(L,K),K=14,33)
  605 CONTINUE
!sl   WRITE(KW,6611) ( TAUSL(K),K=14,33)
      WRITE(KW,6612) (TAUSUM(K),K=14,33)
      WRITE(KW,6613) (BGFLUX(K),K=14,33)
      WRITE(KW,6614) (BGFRAC(K),K=14,33)
      DO 606 I=1,10
      WRITE(KW,6615)
  606 CONTINUE
C
 6601 FORMAT(' (6) TAU TABLE FOR THERMAL RADIATION: INCLUDES ANY'
     +      ,' SPECIFIED OVERLAP, CLOUD & AEROSOL ABSORPTION'
     +      ,T114,'TRGXLK(L,K),TAUSL(L)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6602 FORMAT(' LN     PL     TLM      K=',I1,4X,'K=',I2,9I9,3I8)
 6603 FORMAT(1X,I2,F8.3,F7.2,1X,10F9.4,3F8.3)
 6604 FORMAT(/4X,'SURFACE LAYER= ',10F9.4,3F8.3)
 6605 FORMAT( 4X,'COLUMN AMOUNT= ',10F9.4,3F8.3)
 6606 FORMAT(/1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6607 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6608 FORMAT(/4X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +      ,T83,'OZONE:',T100,'PRINCIPAL ABSORBER REGION'
     +      /4X,76('-'),2X,50('-'))
 6609 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6610 FORMAT( 1X,  I2,6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
!sl6611 FORMAT(/1X,'SL',6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
 6612 FORMAT( 1X,'CA',5F7.4,1F7.3,3F6.2,2F6.1,1F6.0,4F7.4,2F6.3,2F6.2)
 6613 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6614 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6615 FORMAT(' ')
      GO TO 9999
C
C-------------
  700 CONTINUE
C-------------
C
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      SUMK=0.0
      DO 701 K=1,33
      BGFLUX(K) = PLANCK(ITG,K) - (PLANCK(ITG,K)-PLANCK(ITG+1,K))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
  701 CONTINUE
      WRITE(KW,6701)
      WRITE(KW,6702) (K,K=1,13)
      DO 702 L=NL,L1,-1
      WRITE(KW,6703) L,PL(L),TLM(L),(TRCALK(L,K),K=1,13)
  702 CONTINUE
      LK=0
      DO 704 K=1,33
      TAUSUM(K)=0.0
      DO 703 L=L1,NL
      LK=LK+1
  703 TAUSUM(K)=TAUSUM(K)+TRCALK(L,K)
  704 CONTINUE
      WRITE(KW,6704) (TAUSUM(K),K=1,13),(TRCTCA(K),K=1,13)
      WRITE(KW,6705)
      WRITE(KW,6706)  SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6707) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
C
      WRITE(KW,6708)
      WRITE(KW,6709) (K,K=14,33)
      DO 705 L=NL,L1,-1
      WRITE(KW,6710) L,(TRCALK(L,K),K=14,33)
  705 CONTINUE
      WRITE(KW,6711) (TAUSUM(K),K=14,33),(TRCTCA(K),K=14,33)
      WRITE(KW,6712) (BGFLUX(K),K=14,33)
      WRITE(KW,6713) (BGFRAC(K),K=14,33)
      DO 706 I=1,8
      WRITE(KW,6714)
  706 CONTINUE
C
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      SUMK=0.0
      DO 711 K=1,33
      BGFLUX(K) = PLANCK(ITG,K) - (PLANCK(ITG,K)-PLANCK(ITG+1,K))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
  711 CONTINUE
      WRITE(KW,6721)
      WRITE(KW,6722) (K,K=1,13)
      DO 712 L=NL,L1,-1
      WRITE(KW,6723) L,PL(L),TLM(L),(TRAALK(L,K),K=1,13)
  712 CONTINUE
      DO 714 K=1,33
      TAUSUM(K)=0.0
      DO 713 L=L1,NL
  713 TAUSUM(K)=TAUSUM(K)+TRAALK(L,K)
  714 CONTINUE
      WRITE(KW,6724) (TAUSUM(K),K=1,13)
      WRITE(KW,6725)
      WRITE(KW,6726)         SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6727) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
      NPAGE=0
      IF(NL > 13)  NPAGE=1
      WRITE(KW,6728) NPAGE
      WRITE(KW,6729) (K,K=14,33)
      DO 715 L=NL,L1,-1
      WRITE(KW,6730) L,(TRAALK(L,K),K=14,33)
  715 CONTINUE
      WRITE(KW,6731) (TAUSUM(K),K=14,33)
      WRITE(KW,6732) (BGFLUX(K),K=14,33)
      WRITE(KW,6733) (BGFRAC(K),K=14,33)
      DO 716 I=1,12
      WRITE(KW,6734)
  716 CONTINUE
C
 6701 FORMAT(' (7A) TRCALK TABLE FOR THERMAL RADIATION: CONTAINS'
     +      ,' 33 KD CLOUD ABSORPTION OPTICAL DEPTHS AT'
     +      ,' THERMAL WAVELENGTHS ',T117,'LIST: TRCALK(L,K)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6702 FORMAT(' LN     PL     TLM      K=',I1,6X,I2,9I9,3I8)
 6703 FORMAT(1X,I2,F8.3,F7.2,1X,9F9.5,4F8.5)
 6704 FORMAT(/4X,'COLUMN AMOUNT= ',9F9.4,4F8.5
     +       /4X,'TOPCLD ALBEDO= ',9F9.4,4F8.5)
 6705 FORMAT(/' K-INTERVAL CONTRIBUTIONS:'/' COMPARE WITH GROUND FLUX:')
 6706 FORMAT( 1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6707 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6708 FORMAT(/T25,'CARBON DIOXIDE:   PRINCIPAL ABSORBER REGION'
     +      ,T93,'OZONE:   PRINCIPAL ABSORBER REGION'
     +      /4X,77('-'),1X,51('-'))
 6709 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6710 FORMAT( 1X,  I2,6F7.4,5F6.4,F6.4,4F7.4,3F6.4,F6.4)
 6711 FORMAT(/1X,'CA',5F7.4,1F7.3,3F6.2,2F6.1,1F6.0,4F7.4,2F6.3,2F6.2
     +       /1X,'TA',5F7.4,1F7.4,3F6.3,2F6.3,1F6.3,4F7.4,2F6.3,2F6.3)
 6712 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6713 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6714 FORMAT(' ')
C
 6721 FORMAT(' (7B) AEROSOL TAU TABLE FOR THERMAL RADIATION:'
     +      ,'  AEROSOL ABSORPTION OPTICAL DEPTH AT THERMAL WAVELENGTHS'
     +      ,T116,'LIST:  TRAALK(L,K)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6722 FORMAT(' LN     PL     TLM      K=',I1,6X,I2,9I9,3I8)
 6723 FORMAT(1X,I2,F8.3,F7.2,1X,10F9.5,3F8.5)
 6724 FORMAT(/4X,'COLUMN AMOUNT= ',10F9.5,3F8.5)
 6725 FORMAT(' K-INTERVAL CONTRIBUTIONS:'/' COMPARE WITH GROUND FLUX:')
 6726 FORMAT( 1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6727 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6728 FORMAT(1I1/4X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +      ,T83,'OZONE:',T100,'PRINCIPAL ABSORBER REGION'
     +      /4X,76('-'),2X,50('-'))
 6729 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6730 FORMAT( 1X,  I2,6F7.5,2F6.4,3F6.4,F6.4,4F7.4,3F6.4,F6.4)
 6731 FORMAT( 1X,'CA',5F7.5,1F7.5,3F6.4,2F6.4,1F6.4,4F7.4,2F6.4,2F6.4)
 6732 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6733 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6734 FORMAT(' ')
      GO TO 9999
C
C-------------
  800 CONTINUE
C-------------
C
      WRITE(KW,6800)
      DO 801 K=1,16
      ISR1(K)=NORDER(K)
      IF(KORDER==1) ISR1(K)=K
  801 CONTINUE
      WRITE(KW,6801) (ISR1(K),K=1,16)
      SUMK=0.0
      DO 802 K=1,16
      FSR1(K)=DKS0(NORDER(K))
      IF(KORDER==1) FSR1(K)=DKS0(K)
      SUMK=SUMK+FSR1(K)
  802 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6802) (FSR1(K),K=1,17)
      DO 803 K=1,16
      ISR1(K)=NMWAVA(K)
      IF(KORDER==1) ISR1(K)=NMWAVA(IORDER(K))
  803 CONTINUE
      WRITE(KW,6803) (ISR1(K),K=1,16)
      DO 804 K=1,16
      ISR1(K)=NMWAVB(K)
      IF(KORDER==1) ISR1(K)=NMWAVB(IORDER(K))
  804 CONTINUE
      WRITE(KW,6804) (ISR1(K),K=1,16)
      IF(KORDER==0) WRITE(KW,6805)
      IF(KORDER==1) WRITE(KW,6806)
      DO 805 K=1,16
      ISR1(K)=K
      IF(KORDER==1)ISR1(K)=IORDER(K)
  805 CONTINUE
      WRITE(KW,6807) (ISR1(K),K=1,16)
      DO 807 L=NL+1,L1,-1
      SUMK=0.0
      DO 806 K=1,16
      FSR1(K)=SKDFLB(L,NORDER(K))
      IF(KORDER==1) FSR1(K)=SKDFLB(L,K)
      SUMK=SUMK+FSR1(K)
  806 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6808) L,(FSR1(K),K=1,17)
  807 CONTINUE
      DO 808 K=1,16
      ISR1(K)=K
      IF(KORDER==1) ISR1(K)=IORDER(K)
  808 CONTINUE
      WRITE(KW,6809) (ISR1(K),K=1,16)
      DO 810 L=NL+1,L1,-1
      SUMK=0.0
      DO 809 K=1,16
      FSR1(K)=SKUFLB(L,NORDER(K))
      IF(KORDER==1) FSR1(K)=SKUFLB(L,K)
      SUMK=SUMK+FSR1(K)
  809 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6810) L,(FSR1(K),K=1,17)
  810 CONTINUE
      DO 811 K=1,16
      ISR1(K)=K
      IF(KORDER==1) ISR1(K)=IORDER(K)
  811 CONTINUE
      WRITE(KW,6811) (ISR1(K),K=1,16)
      SUMT=0.D0
      SUMK=0.D0
      DO 812 K=1,16
      FSR1(K)=SRKALB(NORDER(K))
      FSR2(K)=DKS0(NORDER(K))
      IF(KORDER==1) FSR1(K)=SRKALB(K)
      IF(KORDER==1) FSR2(K)=DKS0(K)
      SUMK=SUMK+FSR1(K)*FSR2(K)
  812 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6812) (FSR1(K),K=1,17)
      SUMT=SUMT+FSR1(17)
      SUMK1=0.D0
      SUMK2=0.D0
      DO 813 K=1,16
      FSR1(K)=SKNFLB(NL+1,NORDER(K))-SKNFLB(L1,NORDER(K))
      FSR2(K)=SKDFLB(NL+1,NORDER(K))
      IF(KORDER==1) FSR1(K)=SKNFLB(NL+1,K)-SKNFLB(L1,K)
      IF(KORDER==1) FSR2(K)=SKDFLB(NL+1,K)
      SUMK1=SUMK1+FSR1(K)
      SUMK2=SUMK2+FSR2(K)
      FSR1(K)=FSR1(K)/(FSR2(K)+1.d-20)
  813 CONTINUE
      FSR1(17)=SUMK1/(SUMK2+1.d-20)
      WRITE(KW,6813) (FSR1(K),K=1,17)
      SUMT=SUMT+FSR1(17)
      SUMK1=0.D0
      SUMK2=0.D0
      DO 814 K=1,16
      FSR1(K)=SKNFLB(L1,NORDER(K))
      FSR2(K)=SKDFLB(NL+1,NORDER(K))
      IF(KORDER==1) FSR1(K)=SKNFLB(L1,K)
      IF(KORDER==1) FSR2(K)=SKDFLB(NL+1,K)
      SUMK1=SUMK1+FSR1(K)
      SUMK2=SUMK2+FSR2(K)
      FSR1(K)=FSR1(K)/(FSR2(K)+1.d-20)
  814 CONTINUE
      FSR1(17)=SUMK1/(SUMK2+1.d-20)
      WRITE(KW,6814) (FSR1(K),K=1,17)
      SUMT=SUMT+FSR1(17)
      DO 815 K=1,16
      ISR1(K)=KSLAMW(NORDER(K))
      IF(KORDER==1) ISR1(K)=KSLAMW(K)
  815 CONTINUE
      WRITE(KW,6815) SUMT,(ISR1(K),K=1,16)
      SUMK=0.D0
      DO 816 K=1,16
      KK=KSLAMW(NORDER(K))
      IF(KORDER==1) KK=KSLAMW(K)
      FSR1(K)=SRBALB(KK)
      FSR2(K)=SRXALB(KK)
  816 CONTINUE
      WRITE(KW,6816) (FSR1(K),K=1,16)
      WRITE(KW,6817) (FSR2(K),K=1,16)
      WRITE(KW,6818)   COSZ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR
      WRITE(KW,6819) SRXVIS,SRXNIR,SRDVIS,SRUVIS,ALBVIS,SRDNIR
     +                            ,SRUNIR,ALBNIR
      WRITE(KW,6820) SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR
      DO 817 I=1,1
      IF(KORDER==1) WRITE(KW,6821)
  817 CONTINUE
C
      WRITE(KW,6840)
      DO 821 K=1,16
      ISR1(K)=NORDER(K)
      IF(KORDER==1) ISR1(K)=K
  821 CONTINUE
      WRITE(KW,6841) (ISR1(K),K=1,16)
      SUMK=0.0
      DO 822 K=1,16
      FSR1(K)=DKS0(NORDER(K))
      IF(KORDER==1) FSR1(K)=DKS0(K)
      SUMK=SUMK+FSR1(K)
  822 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6842) (FSR1(K),K=1,17)
      IF(KORDER==0) WRITE(KW,6843)
      IF(KORDER==1) WRITE(KW,6844)
      DO 825 K=1,16
      ISR1(K)=K
      IF(KORDER==1)ISR1(K)=IORDER(K)
  825 CONTINUE
      WRITE(KW,6845) (ISR1(K),K=1,16)
      DO 827 L=NL+1,L1,-1
      SUMK=0.0
      DO 826 K=1,16
      FSR1(K)=SKNFLB(L,NORDER(K))
      IF(KORDER==1) FSR1(K)=SKNFLB(L,K)
      SUMK=SUMK+FSR1(K)
  826 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6846) L,(FSR1(K),K=1,17)
  827 CONTINUE
      DO 828 K=1,16
      ISR1(K)=K
      IF(KORDER==1) ISR1(K)=IORDER(K)
  828 CONTINUE
      WRITE(KW,6847) (ISR1(K),K=1,16)
      DO 830 L=NL,L1,-1
      SUMK=0.0
      DO 829 K=1,16
      FSR1(K)=SKFHRL(L,NORDER(K))
      IF(KORDER==1) FSR1(K)=SKFHRL(L,K)
      SUMK=SUMK+FSR1(K)
  829 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6848) L,(FSR1(K),K=1,17)
  830 CONTINUE
      DO 831 K=1,16
      ISR1(K)=K
      IF(KORDER==1)ISR1(K)=IORDER(K)
  831 CONTINUE
      WRITE(KW,6849) (ISR1(K),K=1,16)
      DO 833 N=1,4
      SUMK=0.0
      DO 832 K=1,16
      FSR1(K)=SRKGAX(NORDER(K),N)
      IF(KORDER==1) FSR1(K)=SRKGAX(K,N)
      SUMK=SUMK+FSR1(K)
  832 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6850) 0,(FSR1(K),K=1,17)
  833 CONTINUE
      DO 834 K=1,16
      ISR1(K)=K
      IF(KORDER==1)ISR1(K)=IORDER(K)
  834 CONTINUE
      WRITE(KW,6851)
      DO 836 N=1,4
      SUMK=0.0
      DO 835 K=1,16
      FSR1(K)=SRKGAD(NORDER(K),N)
      IF(KORDER==1) FSR1(K)=SRKGAD(K,N)
      SUMK=SUMK+FSR1(K)
  835 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6852) N,(FSR1(K),K=1,17)
  836 CONTINUE
      WRITE(KW,6853)
      DO 838 N=1,4
      SUMK=0.0
      DO 837 K=1,16
      FSR1(K)=SRKGAX(NORDER(K),N)
      FSR2(K)=SRKGAD(NORDER(K),N)
      IF(KORDER==1) FSR1(K)=SRKGAX(K,N)
      IF(KORDER==1) FSR2(K)=SRKGAD(K,N)
      FSR1(K)=FSR1(K)+FSR2(K)
      SUMK=SUMK+FSR1(K)
  837 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6854) N,(FSR1(K),K=1,17)
  838 CONTINUE
      WRITE(KW,6855)SRNFLB(L1),POCEAN,FSRNFG(1),PEARTH,FSRNFG(2)
     +                        ,POICE ,FSRNFG(3),PLICE ,FSRNFG(4)
C
 6800 FORMAT(' (8A)     SPECTRAL/k-DISTRIBUTION COMPONENT BREAKDOWN'
     +      ,' FOR DOWNWARD AND UPWARD SOLAR RADIATIVE FLUXES'
     +      ,T108,'SKDFLB(L,K)  SKUFLB(L,K)  SRKALB(K)'/)
 6801 FORMAT('     K=',I5,2I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6802 FORMAT('  DKS0=',F6.3,2F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6803 FORMAT(' NMWAVA=',I6,2I8,2I7,I8,I9,7I8,I7,I8)
 6804 FORMAT(' NMWAVB=',I6,2I8,2I7,I8,I9,7I8,I7,I8)
 6805 FORMAT(' ABSORB'/'   GAS= O3,O2  O3,NO2      O2     O2     O2'
     +        ,'     H2O',22X,'H2O     H2O     H2O     H2O     CO2'
     +        ,'     CO2    CO2  CO2,H2O,O2'
     +        /' SKDFLB (Downward Spectral Flux)'
     +        /6X,6('-'),'VIS',6('-'),2X,46('-'),'NIR',59('-'))
 6806 FORMAT(' ABSORB'/'   GAS=   H2O     H2O     H2O    H2O    H2O'
     +        ,'      O2       O2      O2     CO2     CO2     CO2',18X
     +        ,'O3,NO2  O3,O2 CO2,H2O,O2'/'SKDFLB  (Downard Spectral'
     +        ,' Flux)',T110,6('-'),'VIS',5('-'))
 6807 FORMAT('  N  L=',I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6808 FORMAT(I3,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6809 FORMAT(/' SKUFLB (Upward Spectral Flux)'/'  N  L='
     +         ,I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6810 FORMAT(I3,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6811 FORMAT(/' SRKALB '
     +         ,I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6812 FORMAT('   TOA='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6813 FORMAT(' ABSORB'/'  ATMO='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6814 FORMAT(' ABSORB'/'  SURF='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6815 FORMAT(' ALSURF',T133,'Sum=',F6.4
     +       /' KSLAM= ',I3,I9,I8,2I7,I8,I9,7I8,I7,I8)
 6816 FORMAT('   SRX='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6817 FORMAT('   SRB='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6818 FORMAT(/' At Top of Atm:  ',' COSZ  =',F6.4,14X
     +      , 2X,' SRIVIS=',F7.3,'  SROVIS=',F7.3,   '   PLAVIS=',F6.4
     +      , 2X,' SRINIR=',F7.3,'  SRONIR=',F7.3,   '   PLANIR=',F6.4)
 6819 FORMAT( ' At Bot of Atm:  ',' SRXVIS=',F6.4,1X,' SRXNIR=',F6.4
     +      , 1X,' SRDVIS=',F7.3,'  SRUVIS=',F7.3,   '   ALBVIS=',F6.4
     +      , 2X,' SRDNIR=',F7.3,'  SRUNIR=',F7.3,   '   ALBNIR=',F6.4)
 6820 FORMAT( ' In Atmosphere:  ',' (VIS=0.53*S0)',2X,'(NIR=0.47*S0)'
     +      , 1X,' SRTVIS=',F7.5,'  SRRVIS=',F7.5,   '   SRAVIS=',F6.4
     +      , 2X,' SRTNIR=',F7.5,'  SRRNIR=',F7.5,   '   SRANIR=',F6.4)
 6821 FORMAT(' ')

 6840 FORMAT(' (8B)  SPECTRAL/k-DISTRIBUTION COMPONENT BREAKDOWN'
     +      ,' FOR NET DOWNWARD SOLAR FLUX & HEATING RATE'
     +      ,T106,'SKNFLB(L,K)  SKFHRL(L,K)  SRKGAX(L,I)'/)
 6841 FORMAT('     K=',I5,2I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6842 FORMAT('  DKS0=',F6.3,2F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6843 FORMAT('   GAS= O3,O2  O3,NO2      O2     O2     O2'
     +        ,'     H2O',22X,'H2O     H2O     H2O     H2O     CO2',5X
     +        ,'CO2    CO2  CO3,H2O,O2'/' SKNFLB (Spectral Net Flux)')
 6844 FORMAT('   GAS=   H2O     H2O     H2O    H2O    H2O'
     +        ,'      O2       O2      O2     CO2     CO2     CO2',18X
     +        ,'O3,NO2  O3,O2 CO2,H2O,O2'
     +        /' SKDFLB (Spectral Net Flux)',T110,6('-'),'VIS',5('-'))
 6845 FORMAT('  N  L=',I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6846 FORMAT(I3,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6847 FORMAT(/' SKFHRL (Spectral Heating Rate)'/'  N  L='
     +       ,I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6848 FORMAT(I3,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6849 FORMAT(/' SRKGAX (Direct Beam Spectral Absorption at Ground)'
     +       /' N   L=',I4,8I8,I9,4I8,I7,I8,'      Total')
 6850 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6851 FORMAT( ' SRKGAD (Diffuse Spectral Absorption at Ground)')
 6852 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6853 FORMAT( ' SRKGAD (Total Spectral Absorption at Ground)')
 6854 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6855 FORMAT(/' Absorption at Ground by Surface-type'
     +       ,T39,'SRNFLB(1) = POCEAN * FSRNFG(1) + PEARTH * FSRNFG(2) '
     +                     ,'+  POICE * FSRNFG(3) +  PLICE * FSRNFG(4) '
     +       /T39,F7.3,'   = ',F6.4,' *',F8.3,'   + ',F6.4,' *',F8.3
     +                ,'   + ',F6.4,' *',F8.3,'   + ',F6.4,' *',F8.3)
      GO TO 9999
C-------------
  900 CONTINUE
C-------------
C
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      DO 901 K=1,33
      BGFLUX(K) = PLANCK(ITG,K) - (PLANCK(ITG,K)-PLANCK(ITG+1,K))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
  901 CONTINUE
      DO 910 NW=1,5
      DO 903 K=1,33
      DO 902 L=L1,NL+1
      IF(NW==1) WFLB(L,K)=DFLB(L,K)
      IF(NW==2) WFLB(L,K)=UFLB(L,K)
      IF(NW==3) WFLB(L,K)=UFLB(L,K)-DFLB(L,K)
      IF(NW > 3.and.L > NL) GO TO 902
      IF(NW==4) WFLB(L,K)=WFLB(L+1,K)-WFLB(L,K)
      IF(NW==5.and.ABS(TRFCRL(L)) < 1.E-10) WFLB(L,K)=1.E-30
      IF(NW==5) WFLB(L,K)=WFLB(L,K)/(ABS(TRFCRL(L))+1.E-10)
  902 CONTINUE
      IF(NW==1) WFSL(K)=DFSL(K)
      IF(NW==2) WFSL(K)=UFSL(K)
      IF(NW==3) WFSL(K)=UFSL(K)-DFSL(K)
      IF(NW==4) WFSL(K)=WFSL(K)-UFLB(L1,K)+DFLB(L1,K)
!sl   IF(NW==5.and.ABS(TRSLCR) < 1.E-10) WFSL(K)=1.E-30
      IF(NW==5) WFSL(K)=0.   !nu =WFSL(K)/(ABS(TRSLCR)+1.E-10)
  903 CONTINUE
      DO 907 L=L1,NL+1
      IF(L > NL .and. NW > 3) GO TO 907
      ASUM1=0.
      BSUM1=0.
      CSUM1=0.
      DSUM1=0.
      ESUM1=0.
      FSUM1=0.
      SUMF=0.
      DO 904 K=2,13
      ASUM1=ASUM1+  WFSL(K)
      BSUM1=BSUM1+BGFEMT(K)
      CSUM1=CSUM1+BGFLUX(K)
      DSUM1=DSUM1+BGFRAC(K)
      ESUM1=ESUM1+TRCTCA(K)
      FSUM1=FSUM1+TRGALB(K)
  904 SUMF=SUMF+WFLB(L,K)
      SUM1(L)=SUMF
      ASUM2=0.
      BSUM2=0.
      CSUM2=0.
      DSUM2=0.
      ESUM2=0.
      FSUM2=0.
      SUMF=0.
      DO 905 K=14,25
      ASUM2=ASUM2+  WFSL(K)
      BSUM2=BSUM2+BGFEMT(K)
      CSUM2=CSUM2+BGFLUX(K)
      DSUM2=DSUM2+BGFRAC(K)
      ESUM2=ESUM2+TRCTCA(K)
      FSUM2=FSUM2+TRGALB(K)
  905 SUMF=SUMF+WFLB(L,K)
      SUM2(L)=SUMF
      ASUM3=0.
      BSUM3=0.
      CSUM3=0.
      DSUM3=0.
      ESUM3=0.
      FSUM3=0.
      SUMF=0.
      DO 906 K=26,33
      ASUM3=ASUM3+  WFSL(K)
      BSUM3=BSUM3+BGFEMT(K)
      CSUM3=CSUM3+BGFLUX(K)
      DSUM3=DSUM3+BGFRAC(K)
      ESUM3=ESUM3+TRCTCA(K)
      FSUM3=FSUM3+TRGALB(K)
  906 SUMF=SUMF+WFLB(L,K)
      SUM3(L)=SUMF
  907 CONTINUE
C
      NPAGE=1
      WRITE(KW,6901) NW,FTYPE(NW)
      WRITE(KW,6902) (K,K=1,13)
      DO 908 L=NL+1,L1,-1
      IF(L > NL .and. NW > 3) GO TO 908
      SUML=SUM1(L)+SUM2(L)+SUM3(L)+WFLB(L,1)
      WRITE(KW,6903) L,PL(L),SUML,SUM1(L),SUM2(L),SUM3(L)
     +               ,(WFLB(L,K),K=1,13)
  908 CONTINUE
      SUMA=ASUM1+ASUM2+ASUM3+  WFSL(1)
      SUMB=BSUM1+BSUM2+BSUM3+BGFEMT(1)
      SUMC=CSUM1+CSUM2+CSUM3+BGFLUX(1)
      SUMD=DSUM1+DSUM2+DSUM3+BGFRAC(1)
      SUME=ESUM1+ESUM2+ESUM3+TRCTCA(1)
      SUMF=FSUM1+FSUM2+FSUM3+TRGALB(1)
      WRITE(KW,6904) SUMA,ASUM1,ASUM2,ASUM3,(  WFSL(K),K=1,13)
      WRITE(KW,6905) SUMB,BSUM1,BSUM2,BSUM3,(BGFEMT(K),K=1,13)
      WRITE(KW,6906) SUMC,CSUM1,CSUM2,CSUM3,(BGFLUX(K),K=1,13)
      WRITE(KW,6907) SUMD,DSUM1,DSUM2,DSUM3,(BGFRAC(K),K=1,13)
      WRITE(KW,6908) SUME,ESUM1,ESUM2,ESUM3,(TRCTCA(K),K=1,13)
      WRITE(KW,6909) SUMF,FSUM1,FSUM2,FSUM3,(TRGALB(K),K=1,13)
      NPAGE=0
      WRITE(KW,6910) NPAGE
      WRITE(KW,6911) (K,K=14,33)
      DO 909 L=NL+1,L1,-1
      IF(L > NL.and.NW > 3) GO TO 909
      WRITE(KW,6912) L,(WFLB(L,K),K=14,33)
  909 CONTINUE
      WRITE(KW,6913) (  WFSL(K),K=14,33)
      WRITE(KW,6914) (BGFEMT(K),K=14,33)
      WRITE(KW,6915) (BGFLUX(K),K=14,33)
      WRITE(KW,6916) (BGFRAC(K),K=14,33)
      WRITE(KW,6917) (TRCTCA(K),K=14,33)
      WRITE(KW,6918) (TRGALB(K),K=14,33)
      LINFIL=2
      IF(NW > 3) LINFIL=4
      DO 911 I=1,LINFIL
      WRITE(KW,6919)
  911 CONTINUE
  910 CONTINUE
C
 6901 FORMAT(' (9.',I1,') THERMAL RADIATION: K-DISTRIBUTION'
     +       ,' BREAKDOWN FOR  ',1A8,' FLUX'/
     +       /T21,'PRINCIPAL REGION SUM',2X,'WINDOW'
     +       ,T52,'WATER VAPOR:',T76,'PRINCIPAL ABSORBER REGION'
     +       /20X,20('-'),2X,6('-'),3X,81('-'))
 6902 FORMAT(1X,'LN     PL    TOTAL    H2O   CO2    O3     K='
     +       ,I2,4X,'K=',I2,12I7)
 6903 FORMAT(1X,I2,2F8.2,3F7.2,F8.3,1X,12F7.3)
 6904 FORMAT(/' SL',  9X,4F7.2,F8.3,1X,12F7.3)
 6905 FORMAT(/' BG',  9X,4F7.2,F8.3,1X,12F7.3)
 6906 FORMAT( ' PF',  9X,4F7.2,F8.3,1X,12F7.3)
 6907 FORMAT( ' FR',  9X,4F7.2,F8.3,1X,12F7.3)
 6908 FORMAT(/' AC',  9X,4F7.2,F8.3,1X,12F7.3)
 6909 FORMAT( ' AG',  9X,4F7.2,F8.3,1X,12F7.3)
 6910 FORMAT(1I1/5X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +       ,T85,'OZONE:',T101,'PRINCIPAL ABSORBER REGION'
     +       /5X,76('-'),3X,48('-'))
 6911 FORMAT(1X,'LN  K=',I2,6I7,5I6,4X,'K=',I2,1I7,6I6)
 6912 FORMAT( 1X,I2,7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6913 FORMAT(/' SL',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6914 FORMAT(/' BG',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6915 FORMAT( ' PF',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6916 FORMAT( ' FR',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6917 FORMAT(/' AC',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6918 FORMAT( ' AG',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6919 FORMAT(' ')
      RETURN
C-------------
 1000 CONTINUE
C-------------
C
 9999 CONTINUE
      RETURN
      END SUBROUTINE WRITER

      SUBROUTINE WRITET(KWRU,INDEX,JYRREF,JYRNOW,JMONTH,KLIMIT)
      use AerParam_mod, only : updateAerosol,updateAerosol2
      use DustParam_mod, only : upddst2
      use O3mod, only : updO3d,updO3d_solar,plbo3,nlo3
#ifdef HIGH_FREQUENCY_O3_INPUT
      use O3mod, only : UPDO3D_highFrequency
#endif
      IMPLICIT NONE
C
C
C     ------------------------------------------------------------------
C     WRITET  GHG, Solar UV, Ozone, Aerosol Trend Diagnostic Information
C
C         INDEX
C           1     GHG DT0 Trends / FULGAS Ratios for CO2,NO2,CH4,F11,F12
C           2     GHG DF  Change / Ann Increase Rate CO2,NO2,CH4,F11,F12
C           3     Lean Solar Constant, UV Spectral Variation Time Trends
C           4     Ozone Zonal-mean (Latitude and Vertical) Distributions
C           5     Ozone Surface-150mb, 150mb-TOA, Column Longitude Distr
C                 A  O3 (Wang-Jacobs) Relative Longitudinal Distribution
C                 B  O3 (London-NCAR) Relative Longitudinal Distribution
C                 C  O3 (W-J, London) Relative Longitudinal Distribution
C           6     Tropospheric Climatology Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C           7     Tropospheric Desert Dust Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C           8     Stratospheric (Volcanic) Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C           9     Total Column Atmospheric Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C        NOTE:
C                 Time Trend (year) Specification is by JYRREF to JYRNOW
C                 Time Specification (O3,Aerosol) is by JYRREF to JMONTH
C                      (If JMONTH = 0, JDAY is used)
C
C                 INDEX < 10 is selective, INDEX > 10 is digit inclusive
C                 KLIMIT = 0 full output,  KLIMIT > 0 abbreviated output
C                 KWRU directs the output to selected (KWRU) file number
C     ------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: KWRU,INDEX,JYRREF,JYRNOW,JMONTH,KLIMIT

      REAL*8 WREF(7),WDAT(7),WPPM(7),XRAT(5)
      REAL*8, DIMENSION(49,LX) :: QX,QS,QG,QP,O3
      REAL*8, DIMENSION(49) :: QXCOL,QSCOL,QGCOL,QPCOL,O3COL
      REAL*8 SFL0(5),SFLX(5),DFLX(5),RFLX(5),O3L(46,72)
      INTEGER :: LO3(36)
C
      INTEGER, PARAMETER :: NSW1=24, NSW2=32, NSW3=40, NSW4=48
C
      CHARACTER*32, PARAMETER :: CHAER(4) = (/
     *     'Tropospheric Climatology Aerosol',
     +     'Tropospheric Desert Dust Aerosol',
     +     'Stratospheric (Volcanic) Aerosol',
     +     'Total Column Atmospheric Aerosol'/)

      REAL*8 YREF11,ZREF12,SUMO3,QOSH,QONH,QOGL
     *     ,SUMXL,SUMGL,SUMSL,QXSH,QXNH,QXGL,QSSH,QSNH,QSGL,QPSH,QPNH
     *     ,QPGL,QGSH,QGNH,QGGL
      INTEGER KW,INDJ,INDI,INDX,KINDEX,I,JJDAYG,JYEARG,KWSKIP,J,IYEAR
     *     ,NSPACE,LMO,K,mavg,iyr1,lmax,icyc,JYEARS,M,JJDAYO,L,JJ,N,N1
     *     ,N2,II,KAEROS,LL1,KA,JJDAY,icycf
C
      KW=KWRU
      INDJ=MOD(INDEX,10)
      IF(INDJ < 1) INDJ=10
      INDI=1
      IF(INDEX==0)  INDJ=1
      IF(INDEX < 11) INDI=INDJ
      DO 9999 INDX=INDI,INDJ
C
      GO TO (100,100,300,400,500,600,600,600,600,1000),INDX
C
C-------------
  100 CONTINUE
C-------------
C
      KINDEX=INDX
      DO 110 I=1,5
      WREF(I)=XREF(I)
      WPPM(I)=PPMV80(I+4)
  110 CONTINUE
      WREF(6)=PPMV80(11)*1000.D0
      WREF(7)=PPMV80(12)*1000.D0
      WPPM(1)=PPMV80(2)
      WPPM(6)=PPMV80(11)
      WPPM(7)=PPMV80(12)
C
      YREF11=PPMV80(11)
      ZREF12=PPMV80(12)
C
      JJDAYG=184
C
      IF(KINDEX==1) THEN
      WRITE(KW,6101) JJDAYG
 6101 FORMAT(/1X,'(1)=INDEX'
     +      ,T12,'JDAY=',I3,'   RCM RAD EQUIL NO-FEEDBACK DT0'
     +      ,T55,'PRESENT TREND UPDGHG INPUT DATA TO GCM'
     +      ,T96,'FULGAS FACTOR RELATIVE TO 1980 AMOUNTS')
       WRITE(KW,6102) KTREND
 6102 FORMAT(1X,'KTREND=',I2,1X,40('-'),3X,38('-'),3X,38('-')
     +      /1X,'YEAR DTSUM  *DTCO2   DTN2O   DTCH4   DTF11   DTF12'
     +         ,         '   PPMCO2  PPMN20  PPMCH4  PPBF11  PPBF12'
     +         ,         '   FULCO2  FULN2O  FULCH4  FULF11  FULF12')
      ENDIF
C
      IF(KINDEX==2) THEN
      WRITE(KW,6201) JJDAYG
 6201 FORMAT(/1X,'(2)=INDEX'
     +      ,T12,'JDAY=',I3,'   RCM EQ NO-FEEDBACK DFLUX W/M2'
     +      ,T55,'PRESENT TREND UPDGHG INPUT DATA TO GCM'
     +      ,T96,'ANNUAL CHANGE RATE OF TRACE GAS AMOUNT')
      WRITE(KW,6202) KTREND
 6202 FORMAT(1X,'KTREND=',I2,1X,40('-'),3X,38('-'),3X,38('-')
     +      /1X,'YEAR DTSUM  *DTCO2   DTN2O   DTCH4   DTF11   DTF12'
     +         ,         '   PPMCO2  PPMN20  PPMCH4  PPBF11  PPBF12'
     +         ,         '   RATCO2  RATN2O  RATCH4  RATF11  RATF12')
       ENDIF
C
      JYEARG=JYRREF-1
      CALL UPDGHG(JYEARG,JJDAYG)
C
      DO 120 I=1,5
      WDAT(I)=XNOW(I)
  120 CONTINUE
C
      DO 230 J=JYRREF,JYRNOW
      KWSKIP=0
      IF(J > JYRREF) KWSKIP=KLIMIT
      IF(J==1980)   KWSKIP=0
      IF(J==JYRNOW) KWSKIP=0
      JYEARG=J
      CALL UPDGHG(JYEARG,JJDAYG)
      DO 220 I=1,5
      XRAT(I)=(XNOW(I)-WDAT(I))/(1.D-10+WDAT(I))
      IF(XRAT(I) > 9.9999) XRAT(I)=9.9999
      WDAT(I)=XNOW(I)
  220 CONTINUE
      IYEAR=JYEARG
      IF(KINDEX==1) THEN
      IF(KWSKIP==0)
     +WRITE(KW,6103) IYEAR,(XNOW(I),I=1,5),FULGAS(2),(FULGAS(I),I=6,9)
 6103 FORMAT(1X,I4,1X,F8.2,4F8.4,1X,5F8.4)
      ENDIF
      IF(KINDEX==2) THEN
      IF(KWSKIP==0)
     +WRITE(KW,6203) IYEAR,(XNOW(I),I=1,5),(XRAT(I),I=1,5)
 6203 FORMAT(1X,I4,1X,F8.2,4F8.4,1X,5F8.4)
      ENDIF
      NSPACE=IYEAR-(IYEAR/10)*10
      IF(KLIMIT > 0) GO TO 230
      IF(NSPACE==0) WRITE(KW,6104)
 6104 FORMAT(' ')
  230 CONTINUE
      GO TO 9999
C
C-------------
  300 CONTINUE
C-------------
C
      if(ksolar < 0) go to 9999
      LMO=(1950-iy1S0)*12+1
      if(ksolar > 1) LMO=nint(1950 - yr1s0 + 1.5)
      DO 310 I=1,5
      SFL0(I)=0.D0
  310 CONTINUE
      DO 320 K=1,190
      IF(K <= NSW1)              SFL0(1)=SFL0(1)+UV_SSI(LMO,K)*DS_SSI(K)
      IF(K > NSW1.and.K <= NSW2)SFL0(2)=SFL0(2)+UV_SSI(LMO,K)*DS_SSI(K)
      IF(K > NSW2.and.K <= NSW3)SFL0(3)=SFL0(3)+UV_SSI(LMO,K)*DS_SSI(K)
      IF(K > NSW3.and.K <= NSW4)SFL0(4)=SFL0(4)+UV_SSI(LMO,K)*DS_SSI(K)
                                 SFL0(5)=SFL0(5)+UV_SSI(LMO,K)*DS_SSI(K)
  320 CONTINUE
C
      if(ksolar==2.or.ksolar==9)
     *   WRITE(KW,6299) int(yr1s0),int(yr2s0),JYRREF,JYRNOW,SFL0(5)
      if(ksolar < 2) WRITE(KW,6300) JYRREF,JYRNOW,SFL0(5)
 6299 FORMAT(/' (3)=INDEX  Annual-mean Solar flux (from ann. SSI input'
     +      ,I6,'-',I4,' data) for JYRREF=',I4,' to JYRNOW=',I4,'  mid'
     +      ,' 1950 Ref S00WM2=',F9.4/12X,'Solar UV Spectral Flux W/m2'
     +      ,T57,'Delta Solar UV Spectral Flux W/m2'
     +      ,T97,'Solar UV Spectral Flux Ratios'
     +      /'  YEAR    0-280 280-320 320-360 360-400   Total '
     +      ,6X,'0-280 280-320 320-360 360-400   Total '
     +      ,4X,'0-280 280-320 320-360 360-400   Total ')
 6300 FORMAT(/' (3)=INDEX  Annual-mean Solar flux (from J.Lean monthly'
     +      ,' 1882-1998 data) for JYRREF=',I4,' to JYRNOW=',I4,'  Jan'
     +      ,' 1950 Ref S00WM2=',F9.4/12X,'Solar UV Spectral Flux W/m2'
     +      ,T57,'Delta Solar UV Spectral Flux W/m2'
     +      ,T97,'Solar UV Spectral Flux Ratios'
     +      /'  YEAR    0-280 280-320 320-360 360-400   Total '
     +      ,6X,'0-280 280-320 320-360 360-400   Total '
     +      ,4X,'0-280 280-320 320-360 360-400   Total ')
C
      if(ksolar < 2) then
        mavg = 12
        iyr1 = iy1s0
        lmax = ms0x
      else
        mavg = 1
        iyr1 = yr1s0
        lmax = nint( yr2s0-yr1s0+1 )
      end if
      icyc  = mavg*icycs0
      icycf = mavg*icycs0f
      DO 370 J=JYRREF,JYRNOW
      if(j > 2000) icyc = icycf
      KWSKIP=0
      IF(J > JYRREF) KWSKIP=KLIMIT
      IF(J==JYRNOW) KWSKIP=0
      JYEARS=J
      DO 330 I=1,5
      SFLX(I)=0.D0
  330 CONTINUE
      LMO=(JYEARS-iyr1)*mavg
      DO 350 M=1,mavg
      LMO=LMO+1
      IF(LMO > lmax) LMO=LMO-icyc*((LMO-lmax+icyc-1)/icyc)
      IF(LMO < 1) LMO=LMO+icyc*((icyc-LMO)/icyc)
      DO 340 K=1,190
      IF(K <= NSW1)              SFLX(1)=SFLX(1)+UV_SSI(LMO,K)*DS_SSI(K)
      IF(K > NSW1.and.K <= NSW2)SFLX(2)=SFLX(2)+UV_SSI(LMO,K)*DS_SSI(K)
      IF(K > NSW2.and.K <= NSW3)SFLX(3)=SFLX(3)+UV_SSI(LMO,K)*DS_SSI(K)
      IF(K > NSW3.and.K <= NSW4)SFLX(4)=SFLX(4)+UV_SSI(LMO,K)*DS_SSI(K)
                                 SFLX(5)=SFLX(5)+UV_SSI(LMO,K)*DS_SSI(K)
  340 CONTINUE
  350 CONTINUE
      DO 360 I=1,5
      SFLX(I)=SFLX(I)/mavg
      DFLX(I)=SFLX(I)-SFL0(I)
      RFLX(I)=SFLX(I)/SFL0(I)
  360 CONTINUE
      IF(KWSKIP==0)
     +WRITE(KW,6301) JYEARS,(SFLX(I),I=1,5),(DFLX(I),I=1,5)
     +                     ,(RFLX(I),I=1,5)
 6301 FORMAT(2X,I4,1X,4F8.4,F10.4,2X,5F8.4,2X,5F8.5)
      NSPACE=JYEARS-(JYEARS/10)*10
      IF(KLIMIT > 0) GO TO 370
      IF(NSPACE==0) WRITE(KW,6302)
 6302 FORMAT(' ')
  370 CONTINUE
      GO TO 9999
C
C-------------
  400 CONTINUE
C-------------
C
      JJDAYO=JMONTH*30-15
      IF(JMONTH < 1) JJDAYO=JDAY
      CALL UPDO3D(JYRREF,JJDAYO,O3JDAY,O3JREF)
#ifdef HIGH_FREQUENCY_O3_INPUT
      CALL UPDO3D_highFrequency(JYRREF,JJDAYO,O3JDAY_HF_modelLevels)
#endif
      CALL UPDO3D_solar(JJDAYO,S00WM2*RATLS0,O3JDAY)
      DO 450 J=1,46
      DO 410 L=1,NL
      O3(J,L)=0.D0
  410 CONTINUE
      JLAT=J
      DO 430 I=1,72
      ILON=I
!!!   CALL GETO3D(ILON,JLAT)
      CALL REPART(O3JDAY(1,IGCM,JGCM),PLBO3,NLO3+1,U0GAS(1,3),PLB0,NL+1)
      DO 420 L=1,NL
      O3(J,L)=O3(J,L)+U0GAS(L,3)/72.D0
  420 CONTINUE
  430 CONTINUE
      SUMO3=0.D0
      DO 440 L=1,NL
      SUMO3=SUMO3+O3(J,L)
  440 CONTINUE
      O3COL(J)=SUMO3
  450 CONTINUE
      CALL BOXAV1(DLAT46,O3COL,46, 1,23,QOSH)
      CALL BOXAV1(DLAT46,O3COL,46,24,46,QONH)
      CALL BOXAV1(DLAT46,O3COL,46, 1,46,QOGL)
      O3COL(47)=QOSH
      O3COL(48)=QONH
      O3COL(49)=QOGL
      DO 460 L=1,NL
      CALL BOXAV1(DLAT46,O3(1,L),46, 1,23,QOSH)
      CALL BOXAV1(DLAT46,O3(1,L),46,24,46,QONH)
      CALL BOXAV1(DLAT46,O3(1,L),46, 1,46,QOGL)
      O3(47,L)=QOSH
      O3(48,L)=QONH
      O3(49,L)=QOGL
  460 CONTINUE
C
      IF(KLIMIT > 0)
     +WRITE(KW,6400) JYRREF,JJDAYO,JMONTH,MADO3M,(L,L=2,NL)
 6400 FORMAT(/' (4)=INDEX  JYRREF=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone: Zonal-mean Vertical Distribution (cmSTP)'
     +      ,T126,'MADO3=',I2/'  JLAT DLAT46   COLUMN  L =   1',14I7/
     +    I31,14I7)
      IF(KLIMIT < 1)
     +WRITE(KW,7400) JYRREF,JJDAYO,JMONTH,MADO3M
     +      ,(PLB0(I),I=1,15),(L,L=2,15)
      IF(KLIMIT < 1.and.nl > 15) then
        write(KW,'(F33.2,14F7.2)') (PLB0(I),I=16,NL)
        write(KW,'(I31,14I7)') (L,L=16,NL)
      END IF
 7400 FORMAT(/' (4)=INDEX  JYRREF=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone: Zonal-mean Vertical Distribution (cmSTP)'
     +      ,T126,'MADO3=',I2//21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =   1',14I7)
C
      DO 470 JJ=1,46
      J=47-JJ
      IF(KLIMIT > 0) GO TO 470
      WRITE(KW,6401) J,DLAT46(J),O3COL(J),(O3(J,L),L=1,NL)
 6401 FORMAT(I5,F8.2,F9.5,4X,15(1x,F6.5)/26X,15(1x,F6.5))
  470 CONTINUE
      IF(KLIMIT < 1) WRITE(KW,6402)
 6402 FORMAT(' ')
      WRITE(KW,6403) O3COL(48),(O3(48,L),L=1,NL)
 6403 FORMAT(11X,'NH',F9.5,4X,15(1x,F6.5)/26X,15(1x,F6.5))
      IF(KLIMIT < 1) WRITE(KW,6402)
      WRITE(KW,6404) O3COL(47),(O3(47,L),L=1,NL)
 6404 FORMAT(11X,'SH',F9.5,4X,15(1x,F6.5)/26X,15(1x,F6.5))
      IF(KLIMIT < 1) WRITE(KW,6402)
      WRITE(KW,6405) O3COL(49),(O3(49,L),L=1,NL)
 6405 FORMAT( 7X,'GLOBAL',F9.5,4X,15(1x,F6.5)/26X,15(1x,F6.5))
      GO TO 9999
C
C
C-------------
  500 CONTINUE
C-------------
C
      JJDAYO=JMONTH*30-15
      IF(JMONTH < 1) JJDAYO=JDAY
      CALL UPDO3D(JYRREF,JJDAYO,O3JDAY,O3JREF)
#ifdef HIGH_FREQUENCY_O3_INPUT
      CALL UPDO3D_highFrequency(JYRREF,JJDAYO,O3JDAY_HF_modelLevels)
#endif
      CALL UPDO3D_solar(JJDAYO,S00WM2*RATLS0,O3JDAY)
      DO 590 N=1,3
      N1=1
      N2=8
      IF(N==2) N1=9
      IF(N > 1) N2=NL
      DO 530 J=1,46
      JLAT=J
      DO 520 I=1,72
      ILON=I
!!!   CALL GETO3D(ILON,JLAT)
      CALL REPART(O3JDAY(1,IGCM,JGCM),PLBO3,NLO3+1,U0GAS(1,3),PLB0,NL+1)
      SUMO3=0.D0
      DO 510 L=N1,N2
      SUMO3=SUMO3+U0GAS(L,3)
  510 CONTINUE
      O3L(J,I)=SUMO3
  520 CONTINUE
  530 CONTINUE
      DO 560 J=1,46
      SUMO3=0.D0
      DO 540 I=1,72
      SUMO3=SUMO3+O3L(J,I)/72.D0
  540 CONTINUE
      DO 550 I=1,72
      O3L(J,I)=O3L(J,I)/SUMO3
  550 CONTINUE
  560 CONTINUE
C
      IF(N==1)
     +WRITE(KW,6510) JYRREF,JJDAYO,JMONTH,MADO3M,(I,I=10,310,10)
 6510 FORMAT(/'  5A=INDEX  JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone Longitudinal Variation:  Troposphere'
     +      ,' (Wang-Jacobs) Surf to 150 mb',T126,'MADO3=',I2
     +      /'  J LON=0',31I4)
      IF(N==2)
     +WRITE(KW,6520) JYRREF,JJDAYO,JMONTH,MADO3M,(I,I=10,310,10)
 6520 FORMAT(/'  5B=INDEX  JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone Longitudinal Variation:  Stratosphere'
     +      ,' (London-NCAR) 150 mb to TOA',T126,'MADO3=',I2
     +      /'  J LON=0',31I4)
      IF(N==3.and.KLIMIT < 1)
     +WRITE(KW,6530) JYRREF,JJDAYO,JMONTH,MADO3M,(I,I=10,310,10)
 6530 FORMAT(/'  5C=INDEX  JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone Longitudinal Variation:  Total Column'
     +      ,' (W-J/London) Surface to TOA',T126,'MADO3=',I2
     +      /'  J LON=0',31I4)
      IF(KLIMIT < 1) WRITE(KW,6540)
 6540 FORMAT(' ')
C
      IF(N==3.and.KLIMIT > 0) GO TO 590
      DO 580 JJ=1,46
      J=47-JJ
      KWSKIP=KLIMIT
      IF(J==36) KWSKIP=0
      IF(J==24) KWSKIP=0
      IF(J==12) KWSKIP=0
      DO 570 I=1,36
      II=I*2-1
      LO3(I)=O3L(J,II)*100.D0+0.5D0
  570 CONTINUE
      IF(KWSKIP==0)
     +WRITE(KW,6501) J,(LO3(I),I=1,32)
 6501 FORMAT(I4,1X,36I4)
  580 CONTINUE
  590 CONTINUE
C
      GO TO 9999
C
C-------------
  600 CONTINUE
C-------------
C
      KAEROS=4
      IF(INDX==6) KAEROS=1
      IF(INDX==7) KAEROS=2
      IF(INDX==8) KAEROS=3
      LL1=1
      IF(INDX==8 .and. NL > 15) LL1=NL-14
      JJDAY=JMONTH*30-15
      IF(JMONTH < 1) JJDAY=JDAY
      K=6
      IF (MADAER.eq.3) THEN ! newer aerosol fields
      IF(KAEROS==1.OR.KAEROS > 3)
     &       CALL updateAerosol2(JYRREF,JJDAY,a6jday, plbaer)
      ELSE
      IF(KAEROS==1.OR.KAEROS > 3)
     &       CALL updateAerosol(JYRREF,JJDAY, a6jday, plbaer)
      ENDIF
      IF(KAEROS==2.OR.KAEROS > 3) CALL UPDDST2(JYRREF,JJDAY)
      IF(KAEROS==3.OR.KAEROS > 3) CALL UPDVOL(JYRREF,JJDAY)
C
      DO 650 J=1,46
      DO 610 L=1,NL
      QX(J,L)=0.D0
      QS(J,L)=0.D0
      QG(J,L)=0.D0
  610 CONTINUE
      JLAT=J
      DO 630 I=1,72
      ILON=I
      IF(KAEROS==1.OR.KAEROS > 3) CALL GETAER
      IF(KAEROS==2.OR.KAEROS > 3) CALL GETDST
      IF(KAEROS==3.OR.KAEROS > 3) CALL GETVOL
      DO 620 L=1,NL
      IF(KAEROS==1.OR.KAEROS > 3) QX(J,L)=QX(J,L)+SRAEXT(L,K)/72.D0
      IF(KAEROS==2.OR.KAEROS > 3) QX(J,L)=QX(J,L)+SRDEXT(L,K)/72.D0
      IF(KAEROS==3.OR.KAEROS > 3) QX(J,L)=QX(J,L)+SRVEXT(L,K)/72.D0
      IF(KAEROS==1.OR.KAEROS > 3) QS(J,L)=QS(J,L)+SRASCT(L,K)/72.D0
      IF(KAEROS==2.OR.KAEROS > 3) QS(J,L)=QS(J,L)+SRDSCT(L,K)/72.D0
      IF(KAEROS==3.OR.KAEROS > 3) QS(J,L)=QS(J,L)+SRVSCT(L,K)/72.D0
      IF(KAEROS==1.OR.KAEROS > 3) QG(J,L)=QG(J,L)+SRAGCB(L,K)
     +                                              *SRASCT(L,K)/72.D0
      IF(KAEROS==2.OR.KAEROS > 3) QG(J,L)=QG(J,L)+SRDGCB(L,K)
     +                                              *SRDSCT(L,K)/72.D0
      IF(KAEROS==3.OR.KAEROS > 3) QG(J,L)=QG(J,L)+SRVGCB(L,K)
     +                                              *SRVSCT(L,K)/72.D0
  620 CONTINUE
  630 CONTINUE
      SUMXL=1.D-10
      SUMSL=1.D-20
      SUMGL=1.D-20
      DO 640 L=1,NL
      SUMXL=SUMXL+QX(J,L)
      SUMSL=SUMSL+QS(J,L)
      SUMGL=SUMGL+QG(J,L)
      QG(J,L)=(1.D-20+QG(J,L))/(1.D-10+QS(J,L))
      QP(J,L)=(1.D-20+QS(J,L))/(1.D-10+QX(J,L))
      IF(QP(J,L) > 0.99999D0) QP(J,L)=0.99999D0
  640 CONTINUE
      QXCOL(J)=SUMXL
      QSCOL(J)=SUMSL
      QGCOL(J)=(1.D-15+SUMGL)/(1.D-05+SUMSL)
      QPCOL(J)=(1.D-20+SUMSL)/(1.D-10+SUMXL)
  650 CONTINUE
      CALL BOXAV1(DLAT46,QXCOL,46, 1,23,QXSH)
      CALL BOXAV1(DLAT46,QXCOL,46,24,46,QXNH)
      CALL BOXAV1(DLAT46,QXCOL,46, 1,46,QXGL)
      QXCOL(47)=QXSH
      QXCOL(48)=QXNH
      QXCOL(49)=QXGL
      CALL BOXAV1(DLAT46,QSCOL,46, 1,23,QSSH)
      CALL BOXAV1(DLAT46,QSCOL,46,24,46,QSNH)
      CALL BOXAV1(DLAT46,QSCOL,46, 1,46,QSGL)
      QSCOL(47)=QSSH
      QSCOL(48)=QSNH
      QSCOL(49)=QSGL
      CALL BOXAV2(DLAT46,QXCOL,QPCOL,46, 1,23,QPSH)
      CALL BOXAV2(DLAT46,QXCOL,QPCOL,46,24,46,QPNH)
      CALL BOXAV2(DLAT46,QXCOL,QPCOL,46, 1,46,QPGL)
      QPCOL(47)=QPSH
      QPCOL(48)=QPNH
      QPCOL(49)=QPGL
      CALL BOXAV2(DLAT46,QSCOL,QGCOL,46, 1,23,QGSH)
      CALL BOXAV2(DLAT46,QSCOL,QGCOL,46,24,46,QGNH)
      CALL BOXAV2(DLAT46,QSCOL,QGCOL,46, 1,46,QGGL)
      QGCOL(47)=QGSH
      QGCOL(48)=QGNH
      QGCOL(49)=QGGL
      DO 660 L=1,NL
      CALL BOXAV1(DLAT46,QX(1,L),46, 1,23,QXSH)
      CALL BOXAV1(DLAT46,QX(1,L),46,24,46,QXNH)
      CALL BOXAV1(DLAT46,QX(1,L),46, 1,46,QXGL)
      QX(47,L)=QXSH
      QX(48,L)=QXNH
      QX(49,L)=QXGL
      CALL BOXAV1(DLAT46,QS(1,L),46, 1,23,QSSH)
      CALL BOXAV1(DLAT46,QS(1,L),46,24,46,QSNH)
      CALL BOXAV1(DLAT46,QS(1,L),46, 1,46,QSGL)
      QS(47,L)=QSSH
      QS(48,L)=QSNH
      QS(49,L)=QSGL
      CALL BOXAV2(DLAT46,QX(1,L),QP(1,L),46, 1,23,QPSH)
      CALL BOXAV2(DLAT46,QX(1,L),QP(1,L),46,24,46,QPNH)
      CALL BOXAV2(DLAT46,QX(1,L),QP(1,L),46, 1,46,QPGL)
      QP(47,L)=QPSH
      QP(48,L)=QPNH
      QP(49,L)=QPGL
      CALL BOXAV2(DLAT46,QS(1,L),QG(1,L),46, 1,23,QGSH)
      CALL BOXAV2(DLAT46,QS(1,L),QG(1,L),46,24,46,QGNH)
      CALL BOXAV2(DLAT46,QS(1,L),QG(1,L),46, 1,46,QGGL)
      QG(47,L)=QGSH
      QG(48,L)=QGNH
      QG(49,L)=QGGL
  660 CONTINUE
C
      KA=KAEROS
      IF(KLIMIT > 0)
     +WRITE(KW,6600) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA),(L,L=LL1,LL1+14)
 6600 FORMAT(/I3,'A=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL OPTICAL DEPTH',T100,A32
     +      /'  JLAT DLAT46   COLUMN  L =',I4,14I7)
      IF(KLIMIT < 1)
     +WRITE(KW,7600) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA)
     +              ,(PLB0(I),I=LL1,LL1+14),(L,L=LL1,LL1+14)
 7600 FORMAT(/I3,'A=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL OPTICAL DEPTH',T100,A32
     +      //21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =',I4,14I7)
C
      IF(KLIMIT < 1) THEN
      DO 670 JJ=1,46
      J=47-JJ
      WRITE(KW,6601) J,DLAT46(J),QXCOL(J),(QX(J,L),L=LL1,LL1+14)
 6601 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  670 CONTINUE
      WRITE(KW,6602) QXCOL(48),(QX(48,L),L=LL1,LL1+14)
 6602 FORMAT(/11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6603) QXCOL(47),(QX(47,L),L=LL1,LL1+14)
 6603 FORMAT(/11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6604) QXCOL(49),(QX(49,L),L=LL1,LL1+14)
 6604 FORMAT(/7X,'GLOBAL',F9.5,4X,15F7.5)
      ENDIF
      IF(KLIMIT > 0) THEN
      WRITE(KW,6605) QXCOL(48),(QX(48,L),L=LL1,LL1+14)
 6605 FORMAT( 11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6606) QXCOL(47),(QX(47,L),L=LL1,LL1+14)
 6606 FORMAT( 11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6607) QXCOL(49),(QX(49,L),L=LL1,LL1+14)
 6607 FORMAT( 7X,'GLOBAL',F9.5,4X,15F7.5)
      ENDIF
C
      IF(KLIMIT > 0) GO TO 699
      WRITE(KW,6610) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA)
     +              ,(PLB0(I),I=LL1,LL1+14),(L,L=LL1,LL1+14)
 6610 FORMAT(/I3,'B=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL SINGLE SCATTERING ALBEDO'
     +      ,T100,A32//21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =',I4,14I7)
C
      DO 680 JJ=1,46
      J=47-JJ
      WRITE(KW,6611) J,DLAT46(J),QPCOL(J),(QP(J,L),L=LL1,LL1+14)
 6611 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  680 CONTINUE
      WRITE(KW,6612) QPCOL(48),(QP(48,L),L=LL1,LL1+14)
 6612 FORMAT(/11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6613) QPCOL(47),(QP(47,L),L=LL1,LL1+14)
 6613 FORMAT(/11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6614) QPCOL(49),(QP(49,L),L=LL1,LL1+14)
 6614 FORMAT(/7X,'GLOBAL',F9.5,4X,15F7.5)
C
      WRITE(KW,6620) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA)
     +              ,(PLB0(I),I=LL1,LL1+14),(L,L=LL1,LL1+14)
 6620 FORMAT(/I3,'C=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL ASYMMETRY PARAMETER'
     +      ,T100,A32//21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =',I4,14I7)
C
      DO 690 JJ=1,46
      J=47-JJ
      WRITE(KW,6621) J,DLAT46(J),QGCOL(J),(QG(J,L),L=LL1,LL1+14)
 6621 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  690 CONTINUE
      WRITE(KW,6622) QGCOL(48),(QG(48,L),L=LL1,LL1+14)
 6622 FORMAT(/11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6623) QGCOL(47),(QG(47,L),L=LL1,LL1+14)
 6623 FORMAT(/11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6624) QGCOL(49),(QG(49,L),L=LL1,LL1+14)
 6624 FORMAT(/7X,'GLOBAL',F9.5,4X,15F7.5)
  699 CONTINUE
      GO TO 9999
C
 1000 CONTINUE
C
 9999 CONTINUE
C
      RETURN
      END SUBROUTINE WRITET

      END MODULE RADPAR


      SUBROUTINE GTREND(XNOW,TNOW)
C
      USE RADPAR, only: nghg,ghgyr1,ghgyr2,ghgam
      IMPLICIT NONE
      REAL*8 xnow(nghg),tnow,year,dy,frac
      INTEGER iy,n
C
C-------------------------------------------------------------
C        Makiko GHG Trend Compilation  GHG.1850-2050.Dec1999
C
C        Annual-Mean      Greenhouse Gas Mixing Ratios
C-------------------------------------------------------------
C                 CO2     N2O     CH4   CFC-11  CFC-12  others
C        Year     ppm     ppm     ppm     ppb     ppb     ppb
C-------------------------------------------------------------
C     Read from external file - outside table: use value from
C                                      years ghgyr1 or ghgyr2
      YEAR=TNOW
      IF(TNOW <= ghgyr1+.5D0) YEAR=ghgyr1+.5D0
      IF(TNOW >= ghgyr2+.49999D0) YEAR=ghgyr2+.49999D0
      DY=YEAR-(ghgyr1+.5D0)
      IY=DY
      frac=DY-IY
      IY=IY+1
C
C     CO2 N2O CH4 CFC-11 CFC-12 other_GHG  SCENARIO
C--------------------------------------------------
C
      do n=1,nghg
        XNOW(N)=GHGAM(N,IY)+frac*(GHGAM(N,IY+1)-GHGAM(N,IY))
      end do
C
      RETURN
      END SUBROUTINE GTREND

