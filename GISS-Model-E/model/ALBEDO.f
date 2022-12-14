#include "rundeck_opts.h"       
      MODULE SURF_ALBEDO
!@sum SURF_ALBEDO contains parameters/variables needed for albedo calc
!@auth A. Lacis/V. Oinas (modifications by I. Aleinov/G. Schmidt)
      implicit none
      save
      private

      public get_albedo_data, getsur, updsur, albvnh
!@var NKBAND number of K-bands
      INTEGER, PARAMETER :: NKBAND=33

!@var NV total number of vegetation types
      INTEGER, PARAMETER :: NV=12

!@var SRFOAM look up table for ocean foam as a function of wind speed
      REAL*8, PARAMETER :: SRFOAM(25) = (/
     *     0.000,0.000,0.000,0.000,0.001,0.002,0.003,0.005,0.007,0.010,
     *     0.014,0.019,0.025,0.032,0.041,0.051,0.063,0.077,0.094,0.112,
     *     0.138,0.164,0.191,0.218,0.246/)

!@var SEASON julian day for start of season (used for veg albedo calc)
C                                        WINTER SPRING SUMMER AUTUMN
      REAL*8, PARAMETER :: SEASON(4) = (/ 15.00, 105.0, 196.0, 288.0/)
C**** parameters used for vegetation albedo
!@var albvnd veg alb by veg type, season and band
      REAL*8, PARAMETER :: ALBVND(NV,4,6) = RESHAPE( (/
C     (1)  >SRBALB(6) = VIS  (300 - 770 nm)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.067,.089,.089,.078,.100,.067,.061,.089,.000,.200,.089,
     2 .500,.062,.100,.100,.073,.055,.067,.061,.100,.000,.200,.100,
     3 .500,.085,.091,.139,.085,.058,.083,.061,.091,.000,.200,.091,
     4 .500,.080,.090,.111,.064,.055,.061,.061,.090,.000,.200,.090,
C
C     (2)  >SRBALB(5) = NIR  (770 - 860 nm)    (ANIR=Ref)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.200,.267,.267,.233,.300,.200,.183,.267,.000,.200,.267,
     2 .500,.206,.350,.300,.241,.218,.200,.183,.350,.000,.200,.350,
     3 .500,.297,.364,.417,.297,.288,.250,.183,.364,.000,.200,.364,
     4 .500,.255,.315,.333,.204,.218,.183,.183,.315,.000,.200,.315,
C
C     (3)  >SRBALB(4) = NIR  (860 -1250 nm)    (ANIR*1.0)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.200,.267,.267,.233,.300,.200,.183,.267,.000,.200,.267,
     2 .500,.206,.350,.300,.241,.218,.200,.183,.350,.000,.200,.350,
     3 .500,.297,.364,.417,.297,.288,.250,.183,.364,.000,.200,.364,
     4 .500,.255,.315,.333,.204,.218,.183,.183,.315,.000,.200,.315,
C
C     (4)  >SRBALB(3) = NIR  (1250-1500 nm)    (ANIR*0.4)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.080,.107,.107,.093,.120,.080,.073,.107,.000,.200,.107,
     2 .500,.082,.140,.120,.096,.083,.080,.073,.140,.000,.200,.140,
     3 .500,.119,.145,.167,.119,.115,.100,.073,.145,.000,.200,.145,
     4 .500,.102,.126,.132,.081,.087,.073,.073,.126,.000,.200,.126,
C
C     (5)  >SRBALB(2) = NIR  (1500-2200 nm)    (ANIR*0.5)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.100,.133,.133,.116,.150,.100,.091,.133,.000,.200,.133,
     2 .500,.103,.175,.150,.120,.109,.100,.091,.175,.000,.200,.175,
     3 .500,.148,.182,.208,.148,.144,.125,.091,.182,.000,.200,.182,
     4 .500,.127,.157,.166,.102,.109,.091,.091,.157,.000,.200,.157,
C
C     (6)  >SRBALB(1) = NIR  (2200-4000 nm)    (ANIR*0.1)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.020,.027,.027,.023,.030,.020,.018,.027,.000,.200,.027,
     2 .500,.021,.035,.030,.024,.022,.020,.018,.035,.000,.200,.035,
     3 .500,.030,.036,.042,.030,.029,.025,.018,.036,.000,.200,.036,
     4 .500,.026,.032,.033,.020,.022,.018,.018,.032,.000,.200,.032
     *     /),(/NV,4,6/) )
C
!@var VTMASK vegetation depth mask by type (m)
      REAL*8, PARAMETER :: VTMASK(NV) = (/
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
!!!  * 1d1, 2d1, 2d1, 5d1, 2d2, 5d2, 1d3, 25d2,2d1, 1d1,.001d0,2d1 ! (kg/m^2)
     * .1d0,.2d0,.2d0,.5d0,2d0, 5d0, 10d0,25d0,.2d0,.1d0,1d-5,.2d0
     *     /)

!@var ASHZOI,ANHZOI hemisph.Ice Albedo half-max depth (m) (orig.version)
      REAL*8 :: ASHZOI=.1d0, ANHZOI=.1d0             ! tuning parameters
!@var DMOICE  masking depth for snow on sea ice           (orig.version)
      REAL*8 :: DMOICE = 10.                         ! tuning parameter
!@var DMLICE  masking depth for snow on land ice
      REAL*8 :: DMLICE = 10.
!@var AVSCAT,ANSCAT,AVFOAM,ANFOAM for ocean albedo calc
      REAL*8 ::                                      ! tuning parameters
     *     AVSCAT=.0156d0, ANSCAT=0d0, AVFOAM=.2197d0, ANFOAM=.1514d0

!@var ASNALB snow albedo for old snow
!@var AOIALB seaice albedo                            (original version)
!@var ALIALB land ice albedo
      REAL*8, parameter ::
C                        VIS  NIR1  NIR2  NIR3  NIR4  NIR5    NIR
     *     ASNALB(7)=(/.60d0,.55d0,.55d0,.30d0,.10d0,.05d0, .35d0/),
     *     AOIALB(7)=(/.55d0,.50d0,.45d0,.25d0,.10d0,.05d0, .30d0/),
     *     ALIALB(7)=(/.60d0,.55d0,.50d0,.30d0,.10d0,.05d0, .35d0/)
C**** shorthand for the 2 band version
!      do we need these? they are not used
!      REAL*8, parameter :: ASNVIS=ASNALB(1), ASNNIR=ASNALB(7),
!     &     AOIVIS=AOIALB(1), AOINIR=AOIALB(7),
!     &     ALIVIS=ALIALB(1), ALINIR=ALIALB(7)

C**** variables that control snow aging calculation (over land)
!@var AGEXPF exponent in snowage calculation depends on hemi/surf type
!@var ALBDIF difference in albedo as function of snowage
      REAL*8 ::
     *     AGEXPF(3,2) = RESHAPE( (/
C          SH EA   SH OC   SH LI   NH EA   NH OC   NH LI
     *     0.2d0,  0.2d0,  0.2d0,  0.2d0,  0.2d0,  0.2d0 /), (/3,2/) ),
     *     ALBDIF(3,2) = RESHAPE( (/
C          SH EA   SH OC   SH LI   NH EA   NH OC   NH LI
     *     0.35d0, 0.35d0, 0.35d0, 0.35d0, 0.35d0, 0.35d0/), (/3,2/) )

C**** parameters used for Schramm sea ice albedo scheme (Hansen)
!@var AOImin,AOImax           range for seaice albedo
!@var ASNwet,ASNdry           wet,dry snow albedo over sea ice
!@var AMPmin                  mininimal melt pond albedo
      REAL*8 ::
C                         VIS   NIR1   NIR2   NIR3   NIR4   NIR5
     *     AOImin(6)=(/ .05d0, .05d0, .05d0, .050d0, .05d0, .03d0/),
     *     AOImax(6)=(/ .62d0, .42d0, .30d0, .120d0, .05d0, .03d0/),
     *     ASNwet(6)=(/ .85d0, .75d0, .50d0, .175d0, .03d0, .01d0/),
     *     ASNdry(6)=(/ .90d0, .85d0, .65d0, .450d0, .10d0, .10d0/),
     *     AMPmin(6)=(/ .10d0, .05d0, .05d0, .050d0, .05d0, .03d0/)

!@var AOCEAN K-band dependent Thermal radiation characteristics for ocn
      REAL*8, PARAMETER :: AOCEAN(NKBAND) = (/
     +        0.04000,0.09566,0.10273,0.10389,0.10464,0.10555,0.10637,
     +        0.10666,0.10697,0.10665,0.10719,0.10728,0.11007,0.04009,
     +        0.04553,0.05554,0.08178,0.09012,0.09464,0.09548,0.09532,
     +        0.09558,0.09558,0.09568,0.09565,0.05771,0.04985,0.04670,
     +        0.04630,0.04575,0.04474,0.04468,0.04500/)

!@var AGSIDV K-band dependent Thermal radiation for other types
C                     AGSNOW  AGLICE  AGROCK  AGVEG
      REAL*8, PARAMETER :: AGSIDV(NKBAND,4) = RESHAPE( (/
     +        0.01400,0.09262,0.09170,0.07767,0.07130,0.06603,0.06540,
     +        0.06397,0.06358,0.06361,0.06365,0.06386,0.06564,0.01354,
     +        0.01537,0.02320,0.04156,0.03702,0.03633,0.03417,0.03346,
     +        0.03342,0.03322,0.03350,0.03170,0.01967,0.01845,0.01977,
     +        0.01986,0.01994,0.02013,0.02041,0.02100,
     +        0.01400,0.09262,0.09170,0.07767,0.07130,0.06603,0.06540,
     +        0.06397,0.06358,0.06361,0.06365,0.06386,0.06564,0.01354,
     +        0.01537,0.02320,0.04156,0.03702,0.03633,0.03417,0.03346,
     +        0.03342,0.03322,0.03350,0.03170,0.01967,0.01845,0.01977,
     +        0.01986,0.01994,0.02013,0.02041,0.02100,
     +        0.04500,0.10209,0.08806,0.05856,0.04835,0.04052,0.04001,
     +        0.03775,0.03687,0.03740,0.03637,0.03692,0.03570,0.07001,
     +        0.05665,0.05326,0.05349,0.04356,0.03845,0.03589,0.03615,
     +        0.03610,0.03602,0.03613,0.03471,0.13687,0.14927,0.16484,
     +        0.16649,0.16820,0.17199,0.17484,0.18000,
     +        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     *        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /),
     *        (/ NKBAND,4 /) )

C**** miscellaneous tuning parameters
      REAL*8 ::
!@var WETTRA,WETSRA adjustment factors for wet earth albedo calc
     *     WETTRA=1.0, WETSRA=1.0,
!@var ZOCSRA,ZSNSRA,ZICSRA,ZDSSRA,ZVGSRA adjustment factors for
!@+   solar zenith angle effects (not fully enabled)
     *     ZOCSRA=1.0, ZSNSRA=1.0, ZICSRA=1.0, ZDSSRA=1.0, ZVGSRA=1.0,
!@var EOCTRA,ESNTRA,EICTRA,EDSTRA,EVGTRA adjustment factors for
!@+   thermal radiation effects  (not fully enabled)
     *     EOCTRA=1.0, EDSTRA=1.0, ESNTRA=1.0, EICTRA=1.0, EVGTRA=1.0

C**** ALBVNH is set only once a day then saved
!@var ALBVNH hemispherically varying vegetation albedo
      REAL*8 ALBVNH(NV,6,2)

!@var GZSNOW asymmetry parameter for snow over three types
!@+   from Wiscombe and Warren (1980) JAS
!@+   Note this is used for ice + melt-ponds as well.
      REAL*8, PARAMETER :: GZSNOW(7,3,2) = RESHAPE( (/
C       VIS     NIR1    NIR2     NIR3     NIR4     NIR5    NIRT
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,! SH EA
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,!    OI
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,!    LI
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,! NH EA
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,!    OI
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0 !    LI
     *     /), (/7,3,2/) )

      contains

      subroutine get_albedo_data(  AVSCAT_out, ANSCAT_out, AVFOAM_out,
     &     ANFOAM_out, WETTRA_out, WETSRA_out, ZOCSRA_out, ZSNSRA_out,
     &     ZICSRA_out, ZDSSRA_out, ZVGSRA_out, EOCTRA_out, ESNTRA_out,
     &     EICTRA_out, EDSTRA_out, EVGTRA_out, AGEXPF_out, ALBDIF_out )
!@sum returns some internal SURF_ALBEDO data. Needed to provide
!@+   corresponding data to WRITER in RADPAR.
      implicit none
      REAL*8, intent(out) ::       AVSCAT_out, ANSCAT_out, AVFOAM_out,
     &     ANFOAM_out, WETTRA_out, WETSRA_out, ZOCSRA_out, ZSNSRA_out,
     &     ZICSRA_out, ZDSSRA_out, ZVGSRA_out, EOCTRA_out, ESNTRA_out,
     &     EICTRA_out, EDSTRA_out, EVGTRA_out,
     &     AGEXPF_out(3,2), ALBDIF_out(3,2)

      AVSCAT_out = AVSCAT
      ANSCAT_out = ANSCAT
      AVFOAM_out = AVFOAM
      ANFOAM_out = ANFOAM
      WETTRA_out = WETTRA
      WETSRA_out = WETSRA
      ZOCSRA_out = ZOCSRA
      ZSNSRA_out = ZSNSRA
      ZICSRA_out = ZICSRA
      ZDSSRA_out = ZDSSRA
      ZVGSRA_out = ZVGSRA
      EOCTRA_out = EOCTRA
      ESNTRA_out = ESNTRA
      EICTRA_out = EICTRA
      EDSTRA_out = EDSTRA
      EVGTRA_out = EVGTRA
      AGEXPF_out = AGEXPF
      ALBDIF_out = ALBDIF

      return
      end subroutine get_albedo_data

      SUBROUTINE UPDSUR(JYEARR,JJDAYR)
!@sum UPDSUR updates variables for surface albedo once a day
!@auth A. Lacis/V. Oinas (modifications by I. Aleinov/G. Schmidt)
  !    use SURF_ALBEDO, only : albvnh,albvnd,season,nv
      implicit none
!@var jyearr radiation year (not used for anything yet)
!@var jjdayr julian day (used for seasonality of veg albedo)
      INTEGER, intent(in) :: jyearr, jjdayr
      INTEGER K,KS1,KS2,KN1,KN2,L
      REAL*8 XJDAY,SEASN1,SEASN2,WT2,WT1
C
C                      Define Seasonal Albedo Dependence
C                      ---------------------------------
C
      XJDAY=JJDAYR
      SEASN1=-77.0D0
      DO K=1,4
        SEASN2=SEASON(K)
        IF(XJDAY <= SEASN2) GO TO 120
        SEASN1=SEASN2
      END DO
      K=1
      SEASN2=380.0D0
  120 CONTINUE
      WT2=(XJDAY-SEASN1)/(SEASN2-SEASN1)
      WT1=1.D0-WT2
      KS1=1+MOD(K,4)
      KS2=1+MOD(K+1,4)
      KN1=1+MOD(K+2,4)
      KN2=K

      DO K=1,NV
      DO L=1,6
C     -------------------
C     Southern Hemisphere
C     -------------------
        ALBVNH(K,L,1)=WT1*ALBVND(K,KS1,L)+WT2*ALBVND(K,KS2,L)
C     -------------------
C     Northern Hemisphere
C     -------------------
        ALBVNH(K,L,2)=WT1*ALBVND(K,KN1,L)+WT2*ALBVND(K,KN2,L)
      END DO
      END DO
      RETURN
      END SUBROUTINE UPDSUR


      SUBROUTINE GETSUR(
     i     SNOAGE_FAC_MAX,
     i     MLAT46,jnorth,KEEPAL,KSIALB,KZSNOW,MADSUR,
     i     COSZ,PLANCK,PLANCK_TMIN,PLANCK_TMAX,
     i     ILON,JLAT,
     i     AGESN,POCEAN,POICE,PEARTH,PLICE,PLAKE,ZLAKE,
     i     TGO,TGOI,TGE,TGLI,ZOICE,FMP,ZSNWOI,ZMP,
     i     SNOWOI,SNOWD,SNOWLI,SNOW_FRAC,WEARTH,WMAG,PVT,dalbsn,
     i     FLAGS,LOC_CHL,
     o     BXA,PRNB,PRNX,SRBALB,SRXALB,TRGALB,
     o     BGFEMD,BGFEMT,
     o     DTRUFG,FTRUFG
     &     )
!@sum GETSUR computes surface albedo for each grid box
!@auth A. Lacis/V. Oinas (modifications by I. Aleinov/G. Schmidt)

      use ocalbedo_mod, only: ocalbedo

#ifdef SCM
      USE SCM_COM, only : SCMopt,SCMin
#endif

      implicit none

!********* start  in/out *****************************
C**** config data
      REAL*8 SNOAGE_FAC_MAX
      INTEGER  MLAT46,jnorth,KEEPAL,KSIALB,KZSNOW,MADSUR
C**** input from radiation
      INTEGER :: PLANCK_TMIN,PLANCK_TMAX
      REAL*8 COSZ,PLANCK(PLANCK_TMIN:PLANCK_TMAX,NKBAND)  ! T-range(K)
C**** input from driver
      INTEGER ILON,JLAT
      REAL*8 AGESN(3),POCEAN,POICE,PEARTH,PLICE,PLAKE,ZLAKE,
     *     TGO,TGOI,TGE,TGLI,ZOICE,FMP,ZSNWOI,ZMP,
     *     SNOWOI,SNOWD(2),SNOWLI,SNOW_FRAC(2),WEARTH,WMAG,PVT(12),
     &     dalbsn,LOC_CHL
      LOGICAL*4 :: FLAGS
C**** output
      REAL*8 BXA(7),PRNB(6,4),PRNX(6,4),SRBALB(6),SRXALB(6),TRGALB(33),
     &     BGFEMD(33),BGFEMT(33),
     *     DTRUFG(4),FTRUFG(4)
!********* end in/out *****************************
      INTEGER K,J,L,JH,IWM,JWM,ITOC,ITEA,ITOI,ITLI
      REAL*8 ASNAGE,FSNAGE,BOCSUM,BEASUM,BOISUM,BLISUM,AVSCUM,ANSCUM,WMJ
     *   ,WMI,FRFOAM,AV,BV,WTOC,BOCM,BOCP,TRAPOC,BOCM1,BOCP1,BOC
     *   ,DSFRAC,VGFRAC,VTFRAC,WTEA,BEAM,BEAP,TRAPEA
     *   ,BEAM1,BEAP1,BEA,WTOI,BOIM,BOIP,TRAPOI,BOIM1,BOIP1,BOI,AMEAN
     *   ,WTLI,BLIM,BLIP,BGF,TRAPLI,BLIM1,BLIP1,BLI,KKZSNO,FDZICE,AHMZOI

C**** local variables for albedos for each surface type
      REAL*8 BOCVIS,BOCNIR,BLIVIS,BLINIR,XOCVIS,XOCNIR,
     &     EXPSNE,EXPSNO,EXPSNL

C**** arrays needed if 6 band albedo is used
      REAL*8, dimension(6) :: BOCVN,BEAVN,BOIVN,BLIVN,BSNVN,BVNSUR,
     *                        XOCVN,XEAVN,XOIVN,XLIVN,XSNVN,XVNSUR

C**** variables used for sea ice albedo calculation (6 bands, Hansen)
      REAL*8, dimension(6) :: almp6,alsf6
      REAL*8 :: patchy,snagfac,dummy1(33),dummy2(33)

C     -----------------------------------------------------------------
C     Ocean Albedo Dependence on Zenith Angle and Wind Speed
C
      REAL*8 BVH2O, XVH2O, WMAG1, X
      BVH2O(WMAG1)=.0488D0+.0974D0/(5.679D0+WMAG1)+
     &     .0004D0/(.3333D0+WMAG1)
      XVH2O(WMAG1,X)=.021D0+X*X*(.0421D0+X*(.1283D0+X*(-.04D0+X*(3.117D0
     +              /(5.679D0+WMAG1)+X*.025D0/(.3333D0+WMAG1)))))
C     -----------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     Select albedo computation using KSIALB
C     KSIALB= 0  Schramm oi.alb, Antarc/Greenl alb=.8 (J.Hansen)
C     KSIALB= 1  6-band original albedo - no 'fixups' (Andy Lacis)
C       else     Schramm oi.alb but no land ice 'fixups'

C     For offline use: if MADSUR=1 get vegetation fractions from ij-map
      if(MADSUR==1) call getveg(ilon,jlat)

C           Get Albedo, Thermal Flux, Flux Derivative for each Surf Type
C           ------------------------------------------------------------

      JH=1
      IF(JLAT > JNORTH) JH=2
      KKZSNO=KZSNOW
      IF(COSZ < 0.001) KKZSNO=0
C
      EXPSNE=1.D0 ; EXPSNO=1.D0 ;  EXPSNL=1.D0
C
      DO K=1,NKBAND
        TRGALB(K)=0.D0
        BGFEMD(K)=0.D0
        BGFEMT(K)=0.D0
      END DO
C
      BOCSUM=0.D0
      BEASUM=0.D0
      BOISUM=0.D0
      BLISUM=0.D0
      DO K=1,4
        DTRUFG(K)=0.D0
      END DO
C
      AVSCUM=0.D0
      ANSCUM=0.D0
C
      BOCVN=0. ; BEAVN=0. ; BOIVN=0. ; BLIVN=0. ; BSNVN=0.
      XOCVN=0. ; XEAVN=0. ; XOIVN=0. ; XLIVN=0. ; XSNVN=0. ; BXA=0.
C
C                                             --------------------------
C                                             Ocean Albedo Specification
C                                             --------------------------
C
      IF(POCEAN < 1.D-04) GO TO 400
      X=0.5D0+(0.5D0-COSZ)*ZOCSRA
      BOCVIS=BVH2O(WMAG)  +AVSCAT+AVSCUM
      XOCVIS=XVH2O(WMAG,X)+AVSCAT+AVSCUM
      BOCNIR=BVH2O(WMAG)  +ANSCAT+ANSCUM
      XOCNIR=XVH2O(WMAG,X)+ANSCAT+ANSCUM
C
      IWM=WMAG
      IF(IWM < 1) IWM=1
      IF(IWM > 24) IWM=24
      JWM=IWM+1
      WMJ=WMAG-IWM
      WMI=1.D0-WMJ
      FRFOAM=WMI*SRFOAM(IWM)+WMJ*SRFOAM(JWM)
C
      BOCVIS=BOCVIS*(1.D0-FRFOAM)+FRFOAM*AVFOAM
      XOCVIS=XOCVIS*(1.D0-FRFOAM)+FRFOAM*AVFOAM
      BOCNIR=BOCNIR*(1.D0-FRFOAM)+FRFOAM*ANFOAM
      XOCNIR=XOCNIR*(1.D0-FRFOAM)+FRFOAM*ANFOAM
      BOCVN(1)=BOCVIS
      XOCVN(1)=XOCVIS
      DO L=2,6                  ! fill in higher bands
        BOCVN(L)=BOCNIR         ! 1/2 already equivalenced
        XOCVN(L)=XOCNIR
      END DO

      if ( LOC_CHL >= 0.d0 ) then
C**** chlorophyl modification of albedo
! bocvn is the diffuse albedo (function of wind speed)
! xocnv is the direct albedo (function of the solar zenith angle)
! chlorophyll changes only the diffuse albedo
! however, direct albedo calculation in ocalbedo is 
! slightly different than in the default model.


      !call ocalbedo with hycgr=.false. because the
      !calculation is done on the amtos grid here and we
      !need to return bocvn,xocvn but dont return rod and ros
        call  ocalbedo(WMAG,COSZ,BOCVN,XOCVN,
     .     LOC_CHL,dummy1,dummy2,.false.,ILON,JLAT)
      endif

C**** For lakes increase albedo if lakes are very shallow
C**** This is a fix to prevent lakes from overheating when they
C**** are not allowed to evaporate away. We assume that departing
C**** lakes would leave bare soil behind.
C**** Reduce effect by a half...
      if (PLAKE > 0 .and. ZLAKE < 1) then  ! < 1m
        DO L=1,6
          BOCVN(L) = (BOCVN(L)*(ZLAKE-.4d0) +
     *         (1-ZLAKE)*.5*ALBVNH(1,L,1))/.6d0
          XOCVN(L) = (XOCVN(L)*(ZLAKE-.4d0) +
     *         (1-ZLAKE)*.5*ALBVNH(1,L,1))/.6d0
        END DO
      end if
C
      X=1.D0/(1.D0+WMAG)
      AV=(-.0147087D0*X*X+.0292266D0*X-.0081079D0)*EOCTRA
      BV=(1.01673D0-0.0083652D0*WMAG)*EOCTRA
C
      ITOC=TGO
      WTOC=TGO-ITOC
c     IF(ITOC < 124) ITOC=124   ! if this is necessary than something bad has happened!
      BOCSUM=0
      BOCM=0
      BOCP=0
C
      DO K=1,NKBAND
        TRAPOC=AV+BV*AOCEAN(K)
        BOCM1 = (PLANCK(ITOC-1,K)-
     -          (PLANCK(ITOC-1,K)-PLANCK(ITOC,K))*WTOC) * (1-TRAPOC)
        BOCM  =BOCM+BOCM1
        BOCP1 = (PLANCK(ITOC+1,K)-
     -          (PLANCK(ITOC+1,K)-PLANCK(ITOC+2,K))*WTOC) * (1-TRAPOC)
        BOCP  =BOCP+BOCP1
        BOC   = (PLANCK(ITOC,K)-
     -          (PLANCK(ITOC,K)-PLANCK(ITOC+1,K))*WTOC) * (1-TRAPOC)
        BOCSUM=BOCSUM+BOC
C
        TRGALB(K)=TRGALB(K)+POCEAN*TRAPOC
        BGFEMD(K)=BGFEMD(K)+POCEAN*(BOCP1-BOCM1)
        BGFEMT(K)=BGFEMT(K)+POCEAN*BOC
      END DO
      DTRUFG(1)=.5D0*(BOCP-BOCM)
C
  400 CONTINUE
      IF(PEARTH < 1.D-04) GO TO 500
C
C                                         ------------------------------
C                                         Land Snow Albedo Specification
C                                         ------------------------------
      ASNAGE=ALBDIF(1,JH)*EXP(-AGEXPF(1,JH)*AGESN(1))
      DO L=1,6
        FSNAGE=1.D0
        IF(L > 2) FSNAGE=2.0D0/L
        BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
C**** Set zenith angle dependence if required
        IF (KKZSNO > 0) THEN
          CALL RXSNOW(BSNVN(L),COSZ,GZSNOW(L,1,JH),XSNVN(L))
        ELSE
          XSNVN(L)=BSNVN(L)
        END IF
      END DO

C                                          -----------------------------
C                                          Soil/Veg Albedo Specification
C                                          -----------------------------
c**** In the following code when computing albedo for snow covered soil
c**** we use snow_frac(1:2) which is the snow fraction cover for
c**** bare/vegetated soil. It is computed in GHY_DRV.f in accordance
c**** with the surface topography.
c**** The final snow cover is minimum of snow_frac and the snow fraction
c**** obtained using the vegetation masking.
      DSFRAC=PVT(1)+PVT(10)
      VGFRAC=1.D0-DSFRAC
      IF(SNOWD(1)+SNOWD(2)  <= 1.D-04 * 10.d0) THEN
        DO L=1,6
          BEAVN(L)=PVT(1)*ALBVNH(1,L,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        END DO
!nu     BVSOIL=BEAVN(1)
!nu     BNSOIL=BEAVN(2)
        DO K=2,NV
          if ( k==10 .or. k==11 ) cycle
          DO L=1,6
            BEAVN(L)=BEAVN(L)+PVT(K)*ALBVNH(K,L,JH)
          END DO
        END DO
!nu     SEAVIS=BEAVN(1)
!nu     SEANIR=BEAVN(2)
!nu     BVVEGE=BVSOIL
!nu     BNVEGE=BNSOIL
!nu     IF(VGFRAC > 0.001D0) THEN
!nu       BVVEGE=(BEAVN(1)-BVSOIL*DSFRAC)/VGFRAC
!nu       BNVEGE=(BEAVN(2)-BNSOIL*DSFRAC)/VGFRAC
!nu     ENDIF
      ELSE
        VTFRAC=PVT(1)*( 1.d0
     &       - snow_frac(1)*( 1.d0 - EXP(-SNOWD(1)/VTMASK(1)) ) )
        EXPSNE=VTFRAC +
     &       PVT(10)*( 1.d0
     &       - snow_frac(1)*( 1.d0 - EXP(-SNOWD(1)/VTMASK(10)) ) )
        DSFRAC=EXPSNE
        DO L=1,6
          BEAVN(L)=VTFRAC*ALBVNH(1,L,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        END DO
        DO K=2,NV
          if ( k==10 .or. k==11 ) cycle
          VTFRAC=PVT(K)*( 1.d0
     &         - snow_frac(2)*( 1.d0 - EXP(-SNOWD(2)/VTMASK(K)) ) )
          DO L=1,6
            BEAVN(L)=BEAVN(L)+VTFRAC*ALBVNH(K,L,JH)
          END DO
          EXPSNE=EXPSNE+VTFRAC
        END DO
      END IF
      DO L=1,6
        XEAVN(L)=BEAVN(L)
      END DO
      DO L=1,2
        BEAVN(L)=BEAVN(L)+max(0.d0,BSNVN(L)+dalbsn/L)*(1.D0-EXPSNE)
        XEAVN(L)=XEAVN(L)+max(0.d0,XSNVN(L)+dalbsn/L)*(1.D0-EXPSNE)
      END DO
      DO L=3,6
        BEAVN(L)=BEAVN(L)+BSNVN(L)*(1.D0-EXPSNE)
        XEAVN(L)=XEAVN(L)+XSNVN(L)*(1.D0-EXPSNE)
      END DO
      VGFRAC=EXPSNE-DSFRAC
!nu   XVSOIL=BVSOIL
!nu   XNSOIL=BNSOIL
!nu   XVVEGE=BVVEGE
!nu   XNVEGE=BNVEGE

      ITEA=TGE
      WTEA=TGE-ITEA
      BEASUM=0.D0
      BEAM=0.D0
      BEAP=0.D0

      DO K=1,NKBAND
        TRAPEA=AGSIDV(K,1)*(1.D0-EXPSNE)
     +        +AGSIDV(K,3)*DSFRAC*EDSTRA*(1.D0-WETTRA*WEARTH)
     +        +AGSIDV(K,4)*VGFRAC*EVGTRA
        BEAM1 = (PLANCK(ITEA-1,K)-
     -          (PLANCK(ITEA-1,K)-PLANCK(ITEA,K))*WTEA) * (1-TRAPEA)
        BEAM  = BEAM+BEAM1
        BEAP1 = (PLANCK(ITEA+1,K)-
     -          (PLANCK(ITEA+1,K)-PLANCK(ITEA+2,K))*WTEA) * (1-TRAPEA)
        BEAP  = BEAP+BEAP1
        BEA   = (PLANCK(ITEA,K)-
     -          (PLANCK(ITEA,K)-PLANCK(ITEA+1,K))*WTEA) * (1-TRAPEA)
        BEASUM=BEASUM+BEA
C
        TRGALB(K)=TRGALB(K)+PEARTH*TRAPEA
        BGFEMD(K)=BGFEMD(K)+PEARTH*(BEAP1-BEAM1)
        BGFEMT(K)=BGFEMT(K)+PEARTH*BEA
      END DO
      DTRUFG(2)=0.5D0*(BEAP-BEAM)
C
  500 CONTINUE
      IF(POICE < 1.D-04) GO TO 600
C
C                                         ------------------------------
C                                         Ocean Ice Albedo Specification
C                                         ------------------------------
      IF(KSIALB==1) THEN
C****                      original version (6-band)
        EXPSNO=EXP(-SNOWOI/DMOICE)
C**** Set snow albedo over sea ice
        ASNAGE=ALBDIF(2,JH)*EXP(-AGEXPF(2,JH)*AGESN(2))
        DO L=1,6
          FSNAGE=1.D0
          IF(L > 2) FSNAGE=2.0D0/L
          BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
        END DO

C**** set ice albedo
        AHMZOI=ASHZOI
        IF(JLAT > JNORTH) AHMZOI=ANHZOI
        FDZICE=ZOICE/(ZOICE+AHMZOI)   ! ZOICE = ice depth (m)
        DO L=1,6
          BOIVN(L)=FDZICE*AOIALB(L)*EXPSNO+BSNVN(L)*(1.D0-EXPSNO)
C**** Set zenith angle dependence if required
          IF (KKZSNO > 0) THEN
            CALL RXSNOW(BOIVN(L),COSZ,GZSNOW(L,2,JH),XOIVN(L))
          ELSE
            XOIVN(L)=BOIVN(L)
          END IF
        END DO
C**** end of original version (6-band)

      else        ! KSIALB not 1
C**** Schramm/J. Hansen's sea ice albedo formulas (6 spectral bands)
C**** Bare ice:
        BOIVN(1:6) = aoimax(1:6)
        if(ZOICE < 1.)then       ! ZOICE: ice depth (at least Z1I=.1m)
          BOIVN(1:4)=aoimin(1:4)+(aoimax(1:4)-aoimin(1:4))*sqrt(ZOICE)
        endif
C**** Snow:    patchy: snow_cover_fraction (<1 if snow depth < .1m)
        patchy = 0.       !  max(0.d0 , min(1.d0, 10.d0*ZSNWOI) )
        if(ZSNWOI > 0.)then
          if(ZSNWOI >= 0.1d0)then      ! snow deeper than .1m
            patchy=1d0
          else
            patchy=ZSNWOI/0.1d0
          endif
          if(FLAGS)then         ! wet snow
            alsf6(1:6)=asnwet(1:6)
          else                  ! dry snow
            alsf6(1:6)=asndry(1:6)
          endif
C       snow aging based on Loth and Graf (1998)
C       Dry, Wet(thick), Wet(thin) snow decreases by
C       0.006,  0.015 and 0.071 per day, respectively (for mean)
C       assume decrease for each band is proportional
          if (FLAGS) then
            if (ZSNWOI > 0.25) then
              snagfac = 0.015d0/0.7d0 * AGESN(2)
            else
              snagfac = 0.071d0/0.7d0 * AGESN(2)
            end if
          else
            snagfac = 0.006d0/0.82d0  * AGESN(2)
          end if
C       make sure snow albedo doesn't get too low!
          snagfac=min(SNOAGE_FAC_MAX,snagfac)
          alsf6(1:6)=alsf6(1:6)*(1.-snagfac)
C       combine bare ice and snow albedos
          BOIVN(1:6)=BOIVN(1:6)*(1.-patchy)+alsf6(1:6)*patchy
        endif
C**** Melt ponds:
        almp6(1:6)=ampmin(1:6)    ! used if melt pond deeper than .5m
        if(ZMP > 0. .and. ZMP < 0.5)then
          almp6(1:6)=almp6(1:6)+(aoimax(1:6)-almp6(1:6))*(1.-2.*ZMP)**2
        end if
c**** combined sea ice albedo
        BOIVN(1:6)=BOIVN(1:6)*(1.-FMP)+almp6(1:6)*FMP
C**** set zenith angle dependence
        IF (KKZSNO > 0) THEN     ! for all surface types
          DO L=1,6
            CALL RXSNOW(BOIVN(L),COSZ,GZSNOW(L,2,JH),XOIVN(L))
          END DO
        ELSE
          XOIVN(1:6)=BOIVN(1:6)
        END IF
        EXPSNO=1.-patchy

c**** Reduce the Ocean Ice albedo by dalbsn
        DO L=1,2
          BOIVN(L) = max(0.d0,BOIVN(L)+dalbsn/L)
          XOIVN(L) = max(0.d0,XOIVN(L)+dalbsn/L)
        END DO
      end if  !  KSIALB not 1:  Schramm/Hansen
C*
      ITOI=TGOI
      WTOI=TGOI-ITOI
      BOISUM=0.D0
      BOIM=0.D0
      BOIP=0.D0
C
      DO K=1,NKBAND
        TRAPOI=AGSIDV(K,1)*ESNTRA*(1.-EXPSNO)
     +        +AGSIDV(K,2)*EICTRA*EXPSNO
        BOIM1 = (PLANCK(ITOI-1,K)-
     -          (PLANCK(ITOI-1,K)-PLANCK(ITOI,K))*WTOI) * (1-TRAPOI)
        BOIM  =BOIM+BOIM1
        BOIP1 = (PLANCK(ITOI+1,K)-
     -          (PLANCK(ITOI+1,K)-PLANCK(ITOI+2,K))*WTOI) * (1-TRAPOI)
        BOIP  =BOIP+BOIP1
        BOI   = (PLANCK(ITOI,K)-
     -          (PLANCK(ITOI,K)-PLANCK(ITOI+1,K))*WTOI) * (1-TRAPOI)
        BOISUM=BOISUM+BOI
C
        TRGALB(K)=TRGALB(K)+POICE*TRAPOI
        BGFEMD(K)=BGFEMD(K)+POICE*(BOIP1-BOIM1)
        BGFEMT(K)=BGFEMT(K)+POICE*BOI
      END DO
      DTRUFG(3)=0.5D0*(BOIP-BOIM)
C
  600 CONTINUE
      IF(PLICE < 1.E-04) GO TO 700
C                                          -----------------------------
C                                          Land Ice Albedo Specification
C                                          -----------------------------
C**** Set snow albedo over land ice
      ASNAGE=ALBDIF(3,JH)*EXP(-AGEXPF(3,JH)*AGESN(3))
      DO L=1,6
        FSNAGE=1.D0
        IF(L > 2) FSNAGE=2.0D0/L
        BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
      END DO
      XSNVN(1:6)=BSNVN(1:6)

      EXPSNL=EXP(-SNOWLI/DMLICE)
      DO L=1,6
        BLIVN(L)=ALIALB(L)*EXPSNL+BSNVN(L)*(1.D0-EXPSNL)
      END DO

      if (ksialb.eq.0) then                  ! use also land ice fixup:
C****   Specify the Albedo for Antarctica and Greenland: vis.alb = 95%
C****   and mean albedo=80%, i.e. AMEAN = .585*BLIVIS+.415*BLINIR = .80
        IF( JLAT < NINT(MLAT46/6.) .or.
     *     (JLAT < 45.and.JLAT > 38.and.ILON < 33.and.ILON > 23)) THEN
          IF (KKZSNO > 0) THEN
           AMEAN=.78d0 ! compensate for zenith angle effects in the mean
          ELSE
           AMEAN=.8d0
          END IF
          BLIVIS=.95d0
          BLINIR=(AMEAN-.585d0*BLIVIS)/.415d0
          BLIVN(1)=BLIVIS
          DO L=2,6  ! fill in higher bands
            BLIVN(L)=BLINIR
          END DO
        END IF
      end if

C**** zenith angle dependence if required
      DO L=1,6
        IF (KKZSNO > 0) THEN
          CALL RXSNOW(BLIVN(L),COSZ,GZSNOW(L,3,JH),XLIVN(L))
        ELSE
          XLIVN(L)=BLIVN(L)
        END IF
      END DO

c**** Reduce the Land Ice albedo by dalbsn
      DO L=1,2
        BLIVN(L) = max(0.d0,BLIVN(L)+dalbsn/L)
        XLIVN(L) = max(0.d0,XLIVN(L)+dalbsn/L)
      END DO

      ITLI=TGLI
      WTLI=TGLI-ITLI
      BLISUM=0.D0
      BLIM=0.D0
      BLIP=0.D0
      BGF=0.D0
      DO K=1,NKBAND
        TRAPLI=AGSIDV(K,1)*ESNTRA*(1.-EXPSNL) +
     +         AGSIDV(K,2)*EICTRA*EXPSNL
        BLIM1 = (PLANCK(ITLI-1,K)-
     -          (PLANCK(ITLI-1,K)-PLANCK(ITLI,K))*WTLI) * (1-TRAPLI)
        BLIM  =BLIM+BLIM1
        BLIP1 = (PLANCK(ITLI+1,K)-
     -          (PLANCK(ITLI+1,K)-PLANCK(ITLI+2,K))*WTLI) * (1-TRAPLI)
        BLIP  =BLIP+BLIP1
        BLI   = (PLANCK(ITLI,K)-
     -          (PLANCK(ITLI,K)-PLANCK(ITLI+1,K))*WTLI) * (1-TRAPLI)
        BLISUM=BLISUM+BLI
        TRGALB(K)=TRGALB(K)+PLICE*TRAPLI
        BGFEMD(K)=BGFEMD(K)+PLICE*(BLIP1-BLIM1)
        BGFEMT(K)=BGFEMT(K)+PLICE*BLI
      END DO
      DTRUFG(4)=0.5D0*(BLIP-BLIM)
  700 CONTINUE

C**** write some BXA for diagnostic output (in WRITER) (replaces equiv.)
      BXA(1)=EXPSNE ; BXA(2)=EXPSNO ; BXA(3)=EXPSNL
      BXA(4)=BSNVN(1); BXA(5)=BSNVN(2); BXA(6)=XSNVN(1); BXA(7)=XSNVN(2)

C**** calculate final variables always over 6-bands
      DO L=1,6
        BVNSUR(L)=POCEAN*BOCVN(L)+PEARTH*BEAVN(L)
     +           + POICE*BOIVN(L)+ PLICE*BLIVN(L)
        XVNSUR(L)=POCEAN*XOCVN(L)+PEARTH*XEAVN(L)
     +           + POICE*XOIVN(L)+ PLICE*XLIVN(L)
      END DO
      DO L=1,6
        J=7-L
        PRNB(J,1)=BOCVN(L)
        PRNB(J,2)=BEAVN(L)
        PRNB(J,3)=BOIVN(L)
        PRNB(J,4)=BLIVN(L)
        PRNX(J,1)=XOCVN(L)
        PRNX(J,2)=XEAVN(L)
        PRNX(J,3)=XOIVN(L)
        PRNX(J,4)=XLIVN(L)
      END DO
      IF(KEEPAL.NE.1) THEN
        DO J=1,6
          L=7-J
          SRBALB(J)=BVNSUR(L)
          SRXALB(J)=XVNSUR(L)
        END DO
#ifdef SCM
      ELSE
        DO J=1,6
          L=7-J
          if( SCMopt%alb )then
            SRBALB(j)=SCMin%alb
            SRXALB(j)=SCMin%alb
          endif
        END DO
#endif
      ENDIF
C
C                     --------------------------------------------------
C                     Define each Surface Flux Factors, Flux Derivatives
C                     --------------------------------------------------
      BGF=0.D0
      DO K=1,NKBAND
        BGFEMD(K)=BGFEMD(K)*0.5D0
        BGF=BGF+BGFEMT(K)
      END DO
C
!nu   BGM=BOCM*POCEAN+BEAM*PEARTH+BOIM*POICE+BLIM*PLICE
!nu   BGP=BOCP*POCEAN+BEAP*PEARTH+BOIP*POICE+BLIP*PLICE
!nu   TTRUFG=0.5D0*(BGP-BGM)
      FTRUFG(1)=BOCSUM/BGF
      FTRUFG(2)=BEASUM/BGF
      FTRUFG(3)=BOISUM/BGF
      FTRUFG(4)=BLISUM/BGF

      RETURN
      END SUBROUTINE GETSUR

      END MODULE SURF_ALBEDO

