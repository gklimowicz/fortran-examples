#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
      MODULE SEAICE
!@sum  SEAICE contains all the sea ice related subroutines
!@auth Original Development Team
!@cont PREC_SI,SEA_ICE,ADDICE,SIMELT
      USE CONSTANT, only: lhm,rhoi,byrhoi,rhow,shi,shw,byshi,bylhm,rhows
      use TimeConstants_mod, only: SECONDS_PER_DAY
      IMPLICIT NONE
      SAVE

!@param dtdp Clausius-Clapeyron Constant (dT/dp of ice) [K Pa-1]
      real*8, parameter :: dtdp = -7.5d-8

#ifdef TRACERS_WATER
      INTEGER :: ntm
#endif
!@param LMI number of temperature layers in ice
      INTEGER, PARAMETER :: LMI = 4
!@param XSI fractions of mass layer in each temp. layer
!@param BYXSI recipricol of XSI
      REAL*8, PARAMETER, DIMENSION(LMI) ::
     *     XSI= (/0.5d0, 0.5d0, 0.5d0, 0.5d0/),
     *     BYXSI= (/ 1./XSI(1), 1./XSI(2), 1./XSI(3), 1./XSI(4) /)
!@param Z1I thickness of first layer ice [m]
      REAL*8, PARAMETER :: Z1I = .1d0
!@param ACE1I ice mass first layer [kg/m^2]
      REAL*8, PARAMETER :: ACE1I = Z1I*RHOI
!@param Z2OIM min. thickness of 2nd layer ice [m] (if > 0)
      REAL*8, PARAMETER :: Z2OIM = .1d0    ! .4d0
!@param AC2OIM min. ice mass 2nd layer [kg/m^2]   (if > 0)
      REAL*8, PARAMETER :: AC2OIM = Z2OIM*RHOI
C**** snow/ice thermal diffusivity (Pringle et al, 2007)
!@param alami0 lambda coefficient for ice [J m-1 K-1 s-1]
      REAL*8, PARAMETER :: alami0=2.11d0
!@param alams lambda coefficient for snow [J m-1 K-1 s-1]
      REAL*8, PARAMETER :: alams=0.35d0
!@param alamdS salinity/temp coefficient for conductivity (J/(m*s)/psu)
!@param alamdT temp coefficient for conductivity (J/(m*s)/degC^2)
      real*8, parameter :: alamdS=0.09d0, alamdT=-0.011d0
!@param RHOS density of snow [kg/m^3]
      REAL*8, PARAMETER :: RHOS = 300.0
!@var FLEADOC lead fraction for ocean ice (1) for mean ice thickness of 1 m
      REAL*8, PARAMETER :: FLEADOC = 0.06d0
!@var FLEADLK lead fraction for lakes [%]
      REAL*8, PARAMETER :: FLEADLK = 0.
!@var FLEADMX maximum thickness for lead fraction [m]
      REAL*8, PARAMETER :: FLEADMX = 5.
!@param BYHREF (1/m) reciprocal of scale depth for calculating open water
!@+     fraction as a function of mean ice thickness
      REAL*8, PARAMETER :: BYHREF = 1.1d0
!@param BYRLS reciprocal of snow density*lambda
      REAL*8, PARAMETER :: BYRLS = 1./(RHOS*ALAMS)
!@param MU coefficient of seawater freezing point w.r.t. salinity
      REAL*8, PARAMETER :: MU = 0.054d0   ! C/ppt
!@param SSI0 default value for sea ice salinity (=3.2ppt) [kg/kg]
      REAL*8, PARAMETER :: SSI0 = 0.0032d0
!@param FSSS fraction of ocean salinity found in new-formed ice [1]
      REAL*8, PARAMETER :: FSSS = 8d0/35d0
!@var qsfix flag is true if salinity of sea ice is constant
      LOGICAL :: QSFIX = .false.
!@var alpha implicity for heat diffusion in sea ice (1=fully implicit)
      REAL*8, PARAMETER :: ALPHA = 1.0
!@dbparam oi_ustar0 default ice-ocean friction velocity [m/s]
      REAL*8 :: oi_ustar0 = 1d-3  ! 5d-3 ! not used if ice dynamics is
!@dbparam silmfac factor controlling lateral melt of ocean ice
      REAL*8 :: silmfac = 1.d-7 ! = pi*(3d-6)/0.66/1000
!@var silmpow exponent for temperature dependence of lateral melt
      REAL*8 :: silmpow = 1.36d0
!@dbparam snow_ice =1 to allow for snow ice formation (default=1)
      INTEGER :: snow_ice = 1
!@var osurf_tilt controls calc. of ocean surface tilt for ice dyn:
!@+       from geostrophy (=0) or from free surface (=1, default)
      INTEGER :: osurf_tilt = 1
!@param ssimin critical cutoff for salt amount [kg/kg]
      REAL*8 :: ssimin = 1d-4
!@var DEBUG flag
      LOGICAL DEBUG
!@param seaice_thermo formulation of sea ice thermodynamics (BP or SI)
      CHARACTER*2 :: seaice_thermo = "BP"  ! default is BP

      CONTAINS

      SUBROUTINE PREC_SI(SNOW,MSI2,HSIL,TSIL,SSIL,PRCP,ENRGP,RUN0,SRUN0
     *     ,ERUN0,
#ifdef TRACERS_WATER
     &     TRSIL,TRPRCP,TRUN0,
#endif
     &     WETSNOW,CMPRS)
!@sum  PREC_SI Adds the precipitation to sea/lake ice
!@auth Gary Russell
      IMPLICIT NONE
!@param SNOMAX maximum allowed snowdepth (1m equivalent) (kg/m^2)
!@param dSNdRN rate of conversion of snow to ice as a function of rain
      REAL*8, PARAMETER :: SNOMAX=1d0*RHOS, dSNdRN=0.
!@var ENRGP total energy of precip (C),(J/m^2)
!@var PRCP amount of precip (kg/m^2)
      REAL*8, INTENT(IN) :: ENRGP, PRCP
!@var SNOW snow mass (kg/m^2)
!@var MSI1 first layer ice mass (= SNOW + ACE1I) (kg/m^2)
!@var MSI2 second layer ice mass (kg/m^2)
      REAL*8, INTENT(INOUT) :: SNOW, MSI2
!@var HSIL enthalpy of ice layers (J/m^2)
!@var SSIL salt in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL, SSIL
!@var TSIL temperature of ice layers (C) (DIAGNOSTIC ONLY)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var WETSNOW true if snow is wet (i.e. has been rained on)
      LOGICAL, INTENT(OUT) :: WETSNOW
!@var RUN0 runoff from ice (kg/m^2)
!@var SRUN0 salt in runoff from ice (kg/m^2)
!@var ERUN0 energy in runoff from ice (J/m^2)
      REAL*8, INTENT(OUT) :: RUN0, SRUN0, ERUN0
      REAL*8 :: FMSI2, CMPRS, SNWF, RAIN, FREZI, MSI1, FSSI2,
     *     MELTI, SMELTI, HCMPRS, HMELTI, MELTS
      REAL*8 :: MICE(LMI),SICE(LMI),HICE(LMI)
      REAL*8 :: HSNOW(2),SNOWL(2),TSNW(2)
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, DIMENSION(NTM,LMI), INTENT(INOUT) :: TRSIL
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN0
      REAL*8, DIMENSION(NTM) :: FTRSI2,FTRSI3,TRRMF
     *     ,TRMELTI,TRMELTS,TRCMPRS
      REAL*8, DIMENSION(NTM,2) :: TRSNOW
      REAL*8, DIMENSION(NTM,LMI) :: TRICE
#endif

      WETSNOW=.FALSE.
      RUN0=0. ; SRUN0=0. ; ERUN0=0.
#ifdef TRACERS_WATER
      TRUN0=0.
#endif
      MSI1=SNOW+ACE1I

      IF (PRCP.gt.0) THEN
C**** Snowfall is calculated from precip energy (0 deg or colder)
        SNWF = MAX(0d0,MIN(PRCP,-ENRGP*BYLHM))
C**** Rain is remaining precip (0 deg or warmer)
        RAIN = PRCP-SNWF
        WETSNOW = RAIN.GT.1d-5*PRCP ! i.e. a noticeable fraction of prec
C**** set defaults
        MELTS = 0.  ! snow melt
        MELTI = 0. ; SMELTI = 0. ; HMELTI = 0. ! ice melt
        FREZI = 0.               ! rain freezing
        CMPRS = 0.               ! snow compression to ice
#ifdef TRACERS_WATER
        TRMELTS(:) = 0.
        TRMELTI(:) = 0.
        TRCMPRS(:)= 0.
        TRRMF(:) =0.
#endif
        
C**** separate out snow and ice components
        call get_snow_ice_layer(SNOW,MSI2,HSIL,SSIL,
#ifdef TRACERS_WATER
     *       TRSIL,TRSNOW,TRICE, 
#endif 
     *       SNOWL,HSNOW,HICE,SICE,TSNW,TSIL,MICE,.false.)

C**** Steps: 
C****    1) If SNWF > 0, add SNWF to SNOW1
C****    ) If RAIN > 0, calculate MELTS/MELTI or FREZI
C****    ) If SNOW1+SNOW2 > SNOMAX, calculate CMPRS 
C****    ) If RAIN > 0, calculate CMPRS 
C****    ) Calculate resulting fluxes to/from layer 3, then relayer
C****
        IF (SNOW.gt.0) THEN     ! apply fluxes to snow portion first
          MELTS=MAX(0d0,(HSNOW(1)+ENRGP)*BYLHM+SNOWL(1)+SNWF)
          IF (MELTS.gt.SNOWL(1)+SNWF) THEN ! all level 1 snow melts 
            MELTS = SNOWL(1)+SNWF
#ifdef TRACERS_WATER
            TRMELTS(:)= TRSNOW(:,1) + TRPRCP(:)*SNWF/PRCP
#endif
            IF (MICE(1).gt.0) THEN ! check for layer 1 ice melt
              MELTI = Mi(HICE(1)+HSNOW(1)+ENRGP,SICE(1),MICE(1))
              SMELTI = MELTI*SICE(1)/MICE(1)
              HMELTI = MELTI*Em(1d3*SICE(1)/MICE(1))
#ifdef TRACERS_WATER
              TRMELTI(:)= (MELTI-SMELTI)*TRICE(:,1)/(MICE(1)-SICE(1))
#endif
              HICE(1) = HICE(1)+HSNOW(1)+ENRGP-HMELTI
              SICE(1) = SICE(1) - SMELTI
              MICE(1) = MICE(1) - MELTI
#ifdef TRACERS_WATER
              TRICE(:,1)  = TRICE(:,1) - TRMELTI(:)
#endif
            END IF
            HSNOW(1) = 0.
            SNOWL(1) = 0.
#ifdef TRACERS_WATER
            TRSNOW(:,1) = 0.
#endif
          ELSE                   ! some level 1 snow remains
C**** calculate amount of rain that can freeze in the snow and add to ice
            FREZI=MIN(RAIN,MAX(-(HSNOW(1)+ENRGP)*BYLHM-SNOWL(1)-SNWF
     $           +MELTS,0d0))
            HSNOW(1)= HSNOW(1) + ENRGP + LHM*FREZI
            IF (SNOWL(2).gt.0) THEN
               HICE(2) = HICE(2) - LHM*FREZI
               MICE(2) = MICE(2) + FREZI
#ifdef TRACERS_WATER
               TRICE(:,2)  = TRICE(:,2) + TRPRCP(:)*FREZI/PRCP
#endif
            ELSE
               HICE(1) = HICE(1) - LHM*FREZI
               MICE(1) = MICE(1) + FREZI
#ifdef TRACERS_WATER
               TRICE(:,1)  = TRICE(:,1) + TRPRCP(:)*FREZI/PRCP
#endif
            END IF
#ifdef TRACERS_WATER
            TRRMF(:)  = TRPRCP(:)*(1.-(SNWF+FREZI)/PRCP)
            TRSNOW(:,1) = TRSNOW(:,1)+ TRPRCP(:)*SNWF/PRCP
            TRMELTS(:)= MELTS*TRSNOW(:,1)/(SNOWL(1)+SNWF)
            TRSNOW(:,1) = TRSNOW(:,1) - TRMELTS(:)
#endif
            SNOWL(1)=SNOWL(1)+SNWF-MELTS

C**** Calculate potential compression of snow to ice
            IF (SNWF.ge.PRCP .and. SNOW+SNWF .GT. SNOMAX) THEN
C**** TOO MUCH SNOW HAS ACCUMULATED, SOME SNOW IS COMPACTED INTO ICE
              CMPRS = SNOW+SNWF-0.9d0*SNOMAX
            ELSE
C**** RAIN and MELT COMPRESSES SNOW INTO ICE
              CMPRS = MIN(dSNdRN*(RAIN+MELTS), SNOW+SNWF-MELTS)
            END IF
            IF (CMPRS.LT.SNOWL(2)) THEN
             IF(SNOWL(2).NE.0.) THEN
               HICE(2) = HICE(2)+HSNOW(2)*CMPRS/SNOWL(2) 
               HSNOW(2)= HSNOW(2)*(1.-CMPRS/SNOWL(2))
               MICE(2) = MICE(2)+CMPRS
               SNOWL(2)=SNOWL(2)-CMPRS
#ifdef TRACERS_WATER
               TRCMPRS(:)= CMPRS*(TRSNOW(:,1)-TRMELTS(:))/(SNOWL(1)+SNWF
     $              -MELTS) 
               TRICE(:,2)  = TRICE(:,2)  + TRCMPRS(:)
               TRSNOW(:,2) = TRSNOW(:,2) - TRCMPRS(:)
#endif
             END IF
            ELSE
             IF(SNOWL(1).NE.0.) THEN
               HICE(1) = HICE(1)+HSNOW(1)*(CMPRS-SNOWL(2))/SNOWL(1)
               HICE(2) = HICE(2)+HSNOW(2)
               HSNOW(1)= HSNOW(1)*(1.-(CMPRS-SNOWL(2))/SNOWL(1))
               HSNOW(2)= 0. 
               MICE(1) = MICE(1)+CMPRS-SNOWL(2)
               MICE(2) = MICE(2)+SNOWL(2)
               SNOWL(1)= SNOWL(1)-(CMPRS-SNOWL(2))
               SNOWL(2)= 0.
#ifdef TRACERS_WATER
               TRCMPRS(:)= CMPRS*(TRSNOW(:,1)-TRMELTS(:))/(SNOWL(1)+SNWF
     $              -MELTS)
               TRICE(:,1)  = TRICE(:,1) + TRCMPRS(:) - TRSNOW(:,2)
               TRICE(:,2)  = TRICE(:,2) + TRSNOW(:,2) 
               TRSNOW(:,1) = TRSNOW(:,1)- (TRCMPRS(:) - TRSNOW(:,2))
               TRSNOW(:,2) = 0.
#endif
             END IF
            END IF
          END IF
        ELSE                     ! no existing snow
          SNOWL(1) = SNWF
          HSNOW(1) = MIN(ENRGP,0d0) ! new snow
          HICE(1)  = HICE(1)+ENRGP-HSNOW(1)
C**** Calculate possible melt from ice layer 1
          MELTI = Mi(HICE(1),SICE(1),MICE(1))
          SMELTI= MELTI*SICE(1)/MICE(1)
          HMELTI= MELTI*Em(1d3*SICE(1)/MICE(1))
          HICE(1) = HICE(1) - HMELTI
          SICE(1) = SICE(1) - SMELTI
C**** Calculate possible freezing 
          FREZI = Fi(RAIN,HICE(1),SICE(1),MICE(1))
          MICE(1) = MICE(1) - MELTI + FREZI
#ifdef TRACERS_WATER
          TRSNOW(:,1) = TRPRCP(:)*SNWF/PRCP
          TRICE(:,1)  = TRICE(:,1) + TRPRCP(:)*FREZI/PRCP
          TRMELTI(:)= (MELTI-SMELTI)*TRICE(:,1)/(MICE(1)-SICE(1))
          TRICE(:,1)  = TRICE(:,1) - TRMELTI(:)
#endif
        END IF

#ifdef TRACERS_WATER
          TRRMF(:)  = TRPRCP(:)*(1.-(SNWF+FREZI)/PRCP)
#endif
 10     continue

C**** Mass fluxes required to keep first layer ice = ACE1I
        FMSI2 = CMPRS+FREZI-MELTI ! either up or down

C**** relayer lower ice layers
        call relayer(FMSI2,MICE,HICE,SICE
#ifdef TRACERS_WATER
     *       ,TRICE
#endif
     *       )

C**** relayer upper two layers
        call relayer_12(HSNOW,HICE,SICE,MICE,SNOWL
#ifdef TRACERS_WATER
     *       ,TRSNOW,TRICE 
#endif 
     *       )

C**** reconstitute upper snow and ice layers
        call set_snow_ice_layer(HSNOW,HICE,SICE,MICE,SNOWL,
#ifdef TRACERS_WATER
     *       TRSNOW,TRICE,TRSIL, 
#endif 
     *       SNOW,MSI1,MSI2,HSIL,SSIL)

C**** Diagnostics for output
        RUN0 = MELTS+MELTI+(RAIN-FREZI) ! runoff mass to ocean
        IF (RUN0.lt.1d-13) RUN0=0. ! clean up roundoff errors
        SRUN0 = SMELTI          ! salt in runoff
        ERUN0 = HMELTI          ! energy of runoff
        
#ifdef TRACERS_WATER
        TRUN0(:)=TRMELTS(:)+TRMELTI(:)+TRRMF(:) ! tracer flux to ocean
#endif
      END IF

 ! ice temp L=1,..LMI
      call tice(hsil,ssil,msi1,msi2,tsil)

      RETURN
      END SUBROUTINE PREC_SI

      SUBROUTINE SEA_ICE(DTSRCE,SNOW,ROICE,HSIL,SSIL,MSI2,F0DT,F1DT,EVAP
     &     ,SROX,
#ifdef TRACERS_WATER
     &     TRSIL,TREVAP,FTROC,TRUN,
#endif
     *     FMOC,FHOC,FSOC,RUN,ERUN,SRUN,WETSNOW,MELT12,CMPRS)
!@sum  SEA_ICE applies surface fluxes to ice covered areas
!@auth Gary Russell
      IMPLICIT NONE

      REAL*8, PARAMETER :: dSNdML =0.
!@var DTSRCE source time step (s)
      REAL*8, INTENT(IN) :: DTSRCE
!@var F0DT heat flux on the ice top surface (J/m^2)
      REAL*8, INTENT(IN) :: F0DT
!@var F1DT heat flux between the 1st and 2nd ice layers (J/m^2)
      REAL*8, INTENT(IN) :: F1DT
!@var SROX solar radiation at top and bottom ice surface (J/m^2)
      REAL*8, INTENT(INOUT) :: SROX(2)
!@var EVAP evaporation/dew on the top ice surface (kg/m^2)
      REAL*8, INTENT(IN) :: EVAP
!@var FSRI fraction of solar radiation that goes through the bottom
      REAL*8 :: FSRI(LMI)
!@var FMOC,FHOC,FSOC basal fluxes down of mass,heat,salt (J or kg/m^2)
      REAL*8, INTENT(IN) :: FMOC, FHOC, FSOC
!@var RUN,ERUN,SRUN runoff fluxes down of mass,heat,salt (J or kg/m^2)
      REAL*8, INTENT(OUT) :: RUN,ERUN,SRUN
!@var WETSNOW whether snow is wet or not (Used for albedo calc)
      LOGICAL, INTENT(INOUT) :: WETSNOW
!@var MELT12 save surface melt (Used for albedo calc)
      REAL*8, INTENT(OUT) :: MELT12
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, DIMENSION(NTM,LMI), INTENT(INOUT) :: TRSIL
!@var TREVAP tracer amount in evap (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TREVAP
!@var FTROC basal tracer flux to ocean (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: FTROC
!@var TRUN tracer runoff flux to ocean (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN
      REAL*8, DIMENSION(NTM) :: FTRSI2,FTRSI3,TRCMPRS,TRMELTS,TRDEW
      REAL*8, DIMENSION(NTM,LMI) :: TRMELTI,TRICE
      REAL*8, DIMENSION(NTM,2) :: TRSNOW
      INTEGER N
#endif
!@var ROICE sea ice fraction of open water
      REAL*8, INTENT(IN) :: ROICE
!@var HSIL enthalpy of ice layers (J/m^2)
!@var SSIL salt in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL, SSIL
!@var SNOW snow mass (kg/m^2)
!@var MSI1 first layer ice mass (= SNOW + ACE1I) (kg/m^2)
!@var MSI2 second layer ice mass (kg/m^2)
      REAL*8, INTENT(INOUT) :: SNOW, MSI2
!@var TSIL temperature of ice layers (C)
      REAL*8, DIMENSION(LMI) :: TSIL,MICE,HICE,SICE,DMSI,ALAM, MELTI
     $     ,SMELTI,HMELTI,HC,dFdTI,F
      REAL*8 SNOWL(2),HSNOW(2),TSNW(2), DEWI(2)
      REAL*8 MSI1, MELTS, MELTS2, DEW, CMPRS, DEWS, HCMPRS
      REAL*8 FMSI2, FMSI3, FSSI2, FSSI3, FHSI2, FHSI3
      INTEGER L

      FMSI2=0. ; FHSI2=0. ; FHSI3=0. ; FSSI2=0. ; FSSI3=0.
      MELTS=0. ; MELTS2=0. ; MELTI(:)=0. ; SMELTI(:)=0 ; HMELTI(:)=0  
#ifdef TRACERS_WATER
      FTRSI2=0. ; FTRSI3=0.
      TRMELTS=0. ; TRMELTI(:,:)=0. 
#endif
C****
      MSI1 = SNOW+ACE1I ! snow and first (physical) layer ice mass

C**** Calculate solar fractions
      IF (SROX(1).gt.0) THEN
        call solar_ice_frac(snow,msi2,wetsnow,fsri,lmi)
      ELSE
        FSRI = 0.
      END IF
C**** solar flux through bottom
      SROX(2) = SROX(1)*FSRI(LMI)
      F(1)=F1DT  ! accounting

C**** separate out snow and ice components
      call get_snow_ice_layer(SNOW,MSI2,HSIL,SSIL,
#ifdef TRACERS_WATER
     *     TRSIL,TRSNOW,TRICE, 
#endif 
     *     SNOWL,HSNOW,HICE,SICE,TSNW,TSIL,MICE,.true.)

C**** diffusive heat flux through lower ice levels
      DO L=2,LMI-1
 ! effective heat capacities of ice (J/(degC*m^2))
        HC(L) = dEidTi(tsil(L),1d3*(SICE(L)/MICE(L)))*MICE(L)
 ! thermal diffusion at grid edge
        ALAM(L) = Alami((TSIL(L)*MICE(L+1)+TSIL(L+1)*MICE(L))/(MICE(L+1
     $        )+MICE(L)),1d3*(SSIL(L)*MICE(L+1)/MICE(L)+SSIL(L+1)*MICE(L
     $        )/MICE(L+1))/(MICE(L+1)+MICE(L)))
! diffusive heat flux at grid edge
        dFdTI(L) = 2.*ALAM(L)*RHOI*DTSRCE/(MICE(L)+MICE(L+1))

C**** EXPLCIIT DIFFUSIVE FLUXES
C       F(L) = dFdTI(L)*(TSIL(L)-TSIL(L+1))+SROX(1)*FSRI(L)

C**** PARTIALLY IMPLICIT METHOD
        F(L)=(dFdTI(L)*(HC(L)*(TSIL(L)-TSIL(L+1))+ALPHA*F(L-1))+HC(L)
     $       *SROX(1)*FSRI(L))/(HC(L)+ALPHA*dFdTI(L))
      END DO

C**** Update top two thermal layers with surface fluxes
      HSIL(1) = HSIL(1)+(F0DT-F1DT)
      HSIL(2) = HSIL(2)+ F1DT

C**** redo separation into snow and ice layers
      call get_snow_ice_layer(SNOW,MSI2,HSIL,SSIL,
#ifdef TRACERS_WATER
     *     TRSIL,TRSNOW,TRICE,
#endif 
     *     SNOWL,HSNOW,HICE,SICE,TSNW,TSIL,MICE,.false.)

C**** Update ice heat      
      DO L=2,LMI-1
        HICE(L)=HICE(L)-F(L)
        HICE(L+1)=HICE(L+1)+F(L)
      END DO

C**** add bottom fluxes
      HICE(LMI) = HICE(LMI)-SROX(2)-FHOC
      SICE(LMI) = SICE(LMI)-FSOC
#ifdef TRACERS_WATER
      TRICE(:,LMI) = TRICE(:,LMI)-FTROC(:)
#endif
      MICE(LMI)=MICE(LMI)-FMOC

C**** Apply dew/evap to the surface
      DEW = -EVAP
#ifdef TRACERS_WATER
      TRDEW(:) = -TREVAP(:)
#endif

C**** Add DEW directly to ice (which can be in either layer)
C**** DEWI adds to ice, DEWS adds to snow (if applicable)
      DEWI(1) = MAX(0d0,DEW)       ! >0 i.e. dew to ice
      DEWS = DEW - DEWI(1)         ! <0 remaining evap applied to snow
      DEWI(2) = 0.
      IF (SNOWL(1)+DEWS.le.0) THEN ! not enough snow in layer 1
        DEWS=-SNOWL(1)       ! <0 evap up to snow amount
        DEWI(1) = DEW - DEWS       ! either sign, evap/dew applied to ice
c        SNOWL(1)=0.    ! set below
        MICE(1) =MICE(1)+DEWI(1)
      ELSE
c        SNOWL(1)=SNOWL(1)+DEWS   ! set below
        IF (MICE(1).eq.0) THEN  ! first layer is all snow
          DEWI(2)=DEWI(1)
          DEWI(1)=0
          MICE(2) =MICE(2) +DEWI(2)
        ELSE
          MICE(1) =MICE(1) +DEWI(1)
        END IF
      END IF

#ifdef TRACERS_WATER
      DO N=1,NTM
C**** Take into account occasional sign differences when summing
C**** over multiple surface time steps
        IF (DEWI(1)+DEWI(2).gt.0. .and. TRDEW(N).gt.0) THEN
          IF (MICE(1).gt.0) THEN 
            TRICE(N,1) = TRICE(N,1)+TRDEW(N)
          ELSE
            TRICE(N,2) = TRICE(N,2)+TRDEW(N)
          END IF
        ELSEIF (TRSNOW(N,1)+TRDEW(N).gt.0) THEN
          TRSNOW(N,1) = TRSNOW(N,1)+TRDEW(N)
        ELSE
          TRICE(N,1) = TRICE(N,1)+TRSNOW(N,1)+TRDEW(N)
          TRSNOW(N,1) = 0.
        END IF
      END DO
#endif

C**** calculate melt in snow and ice layers
      MELTS = MAX(0d0,HSNOW(1)*BYLHM+SNOWL(1)+DEWS) ! snow melt amount
      MELTS2= MAX(0d0,HSNOW(2)*BYLHM+SNOWL(2))
      DO L=1,LMI
        IF (MICE(L).gt.0) THEN
           MELTI(L)= Mi(HICE(L),SICE(L),MICE(L))
          SMELTI(L)= MELTI(L)*SICE(L)/MICE(L)
          HMELTI(L)= MELTI(L)*Em(1d3*SICE(L)/MICE(L))
#ifdef TRACERS_WATER
          TRMELTI(:,L)=MELTI(L)*TRICE(:,L)/MICE(L)
          TRICE(:,L) = TRICE(:,L) - TRMELTI(:,L)
#endif
          MICE(L) = MICE(L)- MELTI(L)
          IF(MICE(L).EQ.0..AND.dabs(SICE(L)-SMELTI(L)).LT.1.d-12) THEN
            SICE(L) = 0.
          ELSE
            SICE(L) = SICE(L)-SMELTI(L)
          ENDIF
          IF(MICE(L).EQ.0..AND.dabs(HICE(L)-HMELTI(L)).LT.1.d-12) THEN
            HICE(L) = 0.
          ELSE
            HICE(L) = HICE(L)-HMELTI(L)
          ENDIF
        END IF
      END DO

C**** CMPRS is amount of snow turned to ice during melting
C**** Calculate remaining snow and necessary mass flux using
C****        SNOW1 =SNOW1 - MELTS + DEWS - CMPRS
C**** either SNOW2 =SNOW2 - MELTS2               
C****  or    MICE1=MICE1- MELTI1+ DEWI1+ CMPRS 
C****        MICE2=MICE2- MELTI2+ DEWI2+ CMPRS - FMSI2
C****
      IF (SNOWL(1).GT.MELTS-DEWS) THEN ! some L1 snow remains
        CMPRS = MIN(dSNdML*MELTS, SNOWL(1)+DEWS-MELTS) ! > 0.
        HCMPRS = HSNOW(1)*CMPRS/(SNOWL(1)+DEWS-MELTS)
#ifdef TRACERS_WATER
        TRMELTS(:)= MELTS*TRSNOW(:,1)/(SNOWL(1)+DEWS)
        TRCMPRS(:)= CMPRS*(TRSNOW(:,1)-TRMELTS(:))/(SNOWL(1)+DEWS-MELTS)
        TRSNOW(:,1)=TRSNOW(:,1)-TRMELTS(:)-TRCMPRS(:)
#endif
        SNOWL(1)=SNOWL(1)+DEWS-MELTS-CMPRS
        HSNOW(1)=HSNOW(1)-HCMPRS
        IF (MICE(1).gt.0) THEN
          MICE(1)=MICE(1) +  CMPRS
          HICE(1)=HICE(1) + HCMPRS
#ifdef TRACERS_WATER
          TRICE(:,1)=TRICE(:,1) + TRCMPRS(:)
#endif
        ELSE
          MICE(2)=MICE(2) +  CMPRS
          HICE(2)=HICE(2) + HCMPRS
#ifdef TRACERS_WATER
          TRICE(:,2)=TRICE(:,2) + TRCMPRS(:)
#endif
        END IF
      ELSE ! all snow L1 and some ice melts or evaporates
        CMPRS = 0.
#ifdef TRACERS_WATER
        TRMELTS(:)= TRSNOW(:,1)
        TRSNOW(:,1)=0.
#endif
c        MICE(1)=MICE(1) + SNOWL(1)+DEWS-MELTS
        HICE(1)=HICE(1) + HSNOW(1) 
        SNOWL(1)=0.
        HSNOW(1)=0.
      END IF
      IF (MELTS2.gt.0) THEN
        IF (SNOWL(2).GT.MELTS2) THEN
#ifdef TRACERS_WATER
          TRMELTS(:)= TRMELTS(:) + MELTS2*TRSNOW(:,2)/(SNOWL(2))
          TRSNOW(:,2)=TRSNOW(:,2)- MELTS2*TRSNOW(:,2)/(SNOWL(2))
#endif
          SNOWL(2)=SNOWL(2)-MELTS2
        ELSE
#ifdef TRACERS_WATER
          TRMELTS(:)= TRMELTS(:) + TRSNOW(:,2)
          TRSNOW(:,2)=0.
#endif
c         MICE(2)=MICE(2)+SNOWL(2)-MELTS2
          HICE(2)=HICE(2)+HSNOW(2)
          SNOWL(2)=0.
          HSNOW(2)=0.
        END IF
      END IF

C**** Mass fluxes required to keep first layer ice = ACE1I
      FMSI2 = -MELTI(1)-MELTI(2)+DEWI(1)+DEWI(2)+CMPRS  ! either up or down
c      MICE(2)=MICE(2)-FMSI2
c      MSI2=MSI2-SUM(MELTI(3:LMI))+FMSI2-FMOC=SUM(MICE(3:LMI))+FMSI2

C**** relayer lower ice layers
      call relayer(FMSI2,MICE,HICE,SICE
#ifdef TRACERS_WATER
     *     ,TRICE
#endif
     *     )

C**** relayer upper two layers
      call relayer_12(HSNOW,HICE,SICE,MICE,SNOWL
#ifdef TRACERS_WATER
     *     ,TRSNOW,TRICE 
#endif 
     *     )

C**** reconstitute snow and ice layers
      call set_snow_ice_layer(HSNOW,HICE,SICE,MICE,SNOWL,
#ifdef TRACERS_WATER
     *     TRSNOW,TRICE,TRSIL, 
#endif 
     *     SNOW,MSI1,MSI2,HSIL,SSIL)
      
C**** Calculate additional runoff output fluxes
      RUN = MELTS + MELTS2 + SUM( MELTI(:)) ! mass flux to ocean
      SRUN=                  SUM(SMELTI(:)) ! salt flux to ocean
      ERUN= SROX(2)        + SUM(HMELTI(:)) ! energy flux to ocean
#ifdef TRACERS_WATER
      TRUN(:)=TRMELTS(:)+TRMELTI(:,1)+TRMELTI(:,2)+TRMELTI(:,3)
     $     +TRMELTI(:,4)        ! tracer flux to ocean
#endif

C**** update temperatures
      call tice(hsil,ssil,msi1,msi2,tsil)
      
c**** Decide WETSNOW for albedo calculations and save MELT12 for
c**** possible melt ponds
      MELT12=MELTS+MELTI(1)
      WETSNOW = WETSNOW .or. (MELT12.gt.0)
C**** 

      RETURN
      END SUBROUTINE SEA_ICE

      SUBROUTINE ADDICE (SNOW,ROICE,HSIL,SSIL,MSI2,TSIL,ENRGFO
     &     ,ACEFO,ACEFI,ENRGFI,SALTO,SALTI,
#ifdef TRACERS_WATER
     &     TRSIL,TRO,TRI,DTRIMP,
#endif
     *     DMIMP,DHIMP,DSIMP,FLEAD,QFIXR,debug)
!@sum  ADDICE adds ice formed in the ocean to ice variables
!@auth Gary Russell/Gavin Schmidt
      IMPLICIT NONE

c     REAL*8, PARAMETER :: Z2OIX = 5.d0-Z1I,     ! max ice thickness l=2
c    *                     BYZICX=1./(Z1I+Z2OIX)
!@var QFIXR true if RSI and MSI2 are fixed (ie. for fixed SST run)
      LOGICAL, INTENT(IN) :: QFIXR,debug
!@var FLEAD minimum lead fraction for ice (%)
      REAL*8, INTENT(IN) :: FLEAD
!@var ACEFO ice mass formed in ocean for open water fraction (kg/m^2)
!@var ACEFI ice mass formed in ocean for ice covered fraction (kg/m^2)
!@var ENRGFO energy of ice formed in ocean for open water frac (J/m^2)
!@var ENRGFI energy of ice formed in ocean for ice covered frac (J/m^2)
!@var SALTO salt in ice formed in ocean for open water frac (kg/m^2)
!@var SALTI salt in ice formed in ocean for ice covered frac (kg/m^2)
      REAL*8, INTENT(IN) ::  ENRGFI, ENRGFO, ACEFO, ACEFI, SALTO, SALTI
!@var ROICE  ice fraction over open water
!@var SNOW snow amount on ice (kg/m^2)
!@var MSI2 second mass layer ice thickness (kg/m^2)
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
!@var HSIL enthalpy of ice layers (J/m^2)
!@var SSIL salt in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL,SSIL
!@var TSIL temperature of ice layers (C)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var DMIMP,DSIMP,DHIMP implicit mass,salt and heat fluxes required
!@+   to maintain minimum ice thickness if ice fraction is fixed
      REAL*8, INTENT(OUT) :: DMIMP,DSIMP,DHIMP
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: trsil
!@var TRO tracer in ice formed in ocean for open water frac (kg/m^2)
!@var TRI tracer in ice formed in ocean for ice covered frac (kg/m^2)
      REAL*8, INTENT(IN), DIMENSION(NTM) :: tro,tri
!@var DTRIMP implicit tracer flux required to maintain minimum ice
!@+   thickness if ice fraction is fixed
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: DTRIMP
      REAL*8, DIMENSION(NTM) :: FTRSI1,FTRSI2,FTRSI3,FTRSI4
      REAL*8 :: TRSNOW(NTM,2),TRICE(NTM,LMI)
      INTEGER N
#endif
      REAL*8, DIMENSION(LMI) :: FRI
      REAL*8 FMSI1,FMSI2,FMSI3,FMSI4, FHSI1,FHSI2,FHSI3,FHSI4,
     &       FSSI1,FSSI2,FSSI3,FSSI4
      REAL*8 HSNOW(2),SNOWL(2),TSNW(2),HICE(LMI),SICE(LMI),MICE(LMI)
      REAL*8 ROICEN, OPNOCN, DRSI, MSI1, MSI2xx, HAVG, MSI2NO
      integer l

      DMIMP=0. ; DHIMP=0. ; DSIMP=0.
#ifdef TRACERS_WATER
      DTRIMP=0.
#endif
c      if (debug) print*,"add0",SNOW+ACE1I-SSIL(1)-SSIL(2),TRSIL(1,1)+
c     *     TRSIL(1,2),
c     *     ACEFO-SALTO,TRO(1),ACEFI-SALTI,TRI(1),MSI2-SSIL(3)-SSIL(4),
c     *     TRSIL(1,3)+TRSIL(1,4),DMIMP-DSIMP,DTRIMP(1)
      MSI1=SNOW+ACE1I
      IF (.not. QFIXR) THEN
      IF (ROICE.LE.0. .and. ACEFO.gt.0) THEN
C**** Create new ice in ice-free ocean
       ROICE=MIN(1d0,ACEFO/(ACE1I+AC2OIM))
       IF(ACEFO/(ACE1I+AC2OIM).GT.1.) THEN
         write(6,*) "Warning: ADDICE: Excessive sea ice created:",
     *     ACEFO/(ACE1I+AC2OIM)
       ENDIF
       MSI1=ACE1I
       MSI2=MAX(AC2OIM,ACEFO-ACE1I)
       SNOW=0.

       HSIL(1:2) =(ENRGFO/ACEFO)*XSI(1:2)*MSI1
       HSIL(3:4) =(ENRGFO/ACEFO)*XSI(3:4)*MSI2
       SSIL(1:2) = (SALTO/ACEFO)*XSI(1:2)*MSI1
       SSIL(3:4) = (SALTO/ACEFO)*XSI(3:4)*MSI2
#ifdef TRACERS_WATER
       DO N=1,NTM
         TRSIL(N,1:2) = (TRO(N)/ACEFO)*XSI(1:2)*MSI1
         TRSIL(N,3:4) = (TRO(N)/ACEFO)*XSI(3:4)*MSI2
       END DO
#endif
       
      ELSEIF (ROICE.gt.0) THEN
C**** Create new ice in partially ice-covered ocean
      IF (ACEFI.gt.0) THEN
C**** CALCULATE ADVECTIVE HEAT FLUX FROM LAYER 3 TO LAYER 4 OF ICE
CC      FMSI3 = -XSI(3)*ACEFI ! < 0.
CC      FMSI4 = -ACEFI
        IF (XSI(3)*ACEFI.gt.XSI(4)*MSI2) THEN ! exceptional case
          FHSI3 = -HSIL(4)-(XSI(3)*ACEFI-XSI(4)*MSI2)*ENRGFI/ACEFI
          FSSI3 = -SSIL(4)-(XSI(3)*ACEFI-XSI(4)*MSI2)*SALTI/ACEFI
#ifdef TRACERS_WATER
          FTRSI3(:) = -TRSIL(:,4)-(XSI(3)*ACEFI-XSI(4)*MSI2)*TRI(:)
     $         /ACEFI
#endif
        ELSE
          FHSI3 = -HSIL(4)*ACEFI*(XSI(3)/XSI(4))/MSI2
          FSSI3 = -SSIL(4)*ACEFI*(XSI(3)/XSI(4))/MSI2
#ifdef TRACERS_WATER
          FTRSI3(:) = -TRSIL(:,4)*ACEFI*(XSI(3)/XSI(4))/MSI2
#endif
        END IF
      ELSE
        FHSI3 = 0. ; FSSI3 = 0.
#ifdef TRACERS_WATER
        FTRSI3(:) = 0.
#endif
      END IF

      IF (ACEFO .EQ. 0.) THEN
C**** NEW ICE IS ONLY FORMED BELOW OLD SEA ICE
        HSIL(3) = HSIL(3)-FHSI3
        HSIL(4) = HSIL(4)+(FHSI3+ENRGFI)
        SSIL(3) = SSIL(3)-FSSI3
        SSIL(4) = SSIL(4)+(FSSI3+SALTI)
#ifdef TRACERS_WATER
        TRSIL(:,3) = TRSIL(:,3) - FTRSI3(:)
        TRSIL(:,4) = TRSIL(:,4) +(FTRSI3(:)+TRI(:))
#endif
        MSI2 = MSI2+ACEFI       ! new ice mass of physical layer 2
      ELSE
C**** NEW ICE IS FORMED ON OPEN OCEAN AND POSSIBLY BELOW OLD SEA ICE
        DRSI = MIN((1.-ROICE)*ACEFO/(ACE1I+AC2OIM),1.-ROICE) ! new ice on the open oc.
        ROICEN=ROICE+DRSI
        IF(ACEFO/(ACE1I+AC2OIM).GT.1.) THEN
          write(6,*) "Warning: ADDICE: Excessive sea ice created:",
     *      ACEFO/(ACE1I+AC2OIM)
        ENDIF       

c**** new 2nd layer ice thickness on open water fraction
        MSI2NO = MAX(AC2OIM,ACEFO-ACE1I)

C**** separate out snow and ice components
      call get_snow_ice_layer(SNOW,MSI2,HSIL,SSIL,
#ifdef TRACERS_WATER
     *     TRSIL,TRSNOW,TRICE, 
#endif 
     *     SNOWL,HSNOW,HICE,SICE,TSNW,TSIL,MICE,.false.)

C**** distribute snow variables over new ice extent
        SNOWL(:) = SNOWL(:)*(ROICE/ROICEN)
        HSNOW(:) = HSNOW(:)*(ROICE/ROICEN)
CC      SSNOW = SSNOW*(ROICE/ROICEN)  ! always zero
#ifdef TRACERS_WATER
        TRSNOW(:,1:2) = TRSNOW(:,1:2)*(ROICE/ROICEN)
#endif

C**** Add new ice to ice variables
        HICE(1:2)=((1.-ROICE)*ENRGFO*XSI(1:2)*ACE1I/(ACE1I+MSI2NO)+ROICE
     $       *HICE(1:2))/ROICEN
        SICE(1:2)=((1.-ROICE)*SALTO *XSI(1:2)*ACE1I/(ACE1I+MSI2NO)+ROICE
     $       *SICE(1:2))/ROICEN
#ifdef TRACERS_WATER
        TRICE(:,1)=((1.-ROICE)*TRO(:)*XSI(1)*ACE1I/(ACE1I+MSI2NO)+ROICE
     $       *TRICE(:,1))/ROICEN
        TRICE(:,2)=((1.-ROICE)*TRO(:)*XSI(2)*ACE1I/(ACE1I+MSI2NO)+ROICE
     $       *TRICE(:,2))/ROICEN
#endif
        MICE(1:2)=((1.-ROICE)*ACEFO*XSI(1:2)*ACE1I/(ACE1I+MSI2NO)+ROICE
     $       *MICE(1:2))/ROICEN

C**** relayer upper two layers
        call relayer_12(HSNOW,HICE,SICE,MICE,SNOWL
#ifdef TRACERS_WATER
     *       ,TRSNOW,TRICE 
#endif 
     *       )

C**** reconstitute snow and ice layers
        call set_snow_ice_layer(HSNOW,HICE,SICE,MICE,SNOWL,
#ifdef TRACERS_WATER
     *       TRSNOW,TRICE,TRSIL, 
#endif 
     *       SNOW,MSI1,MSI2,HSIL,SSIL)

C**** add new ice to old ice (not including snow)
        MSI2 = (DRSI*MAX(AC2OIM,ACEFO-ACE1I)+ROICE*(MSI2+ACEFI))/ROICEN ! layer 2
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
        HSIL(3)=((1.-ROICE)*ENRGFO*XSI(3)*MSI2NO/(ACE1I+MSI2NO)
     *           +ROICE*(HSIL(3)-FHSI3))/ROICEN
        HSIL(4)=((1.-ROICE)*ENRGFO*XSI(4)*MSI2NO/(ACE1I+MSI2NO)
     *           +ROICE*(HSIL(4)+FHSI3+ENRGFI)
     *       )/ROICEN
        SSIL(3)=((1.-ROICE)*SALTO*XSI(3)*MSI2NO/(ACE1I+MSI2NO)
     *           +ROICE*(SSIL(3)-FSSI3))/ROICEN
        SSIL(4)=((1.-ROICE)*SALTO*XSI(4)*MSI2NO/(ACE1I+MSI2NO)
     *           +ROICE*(SSIL(4)+FSSI3+SALTI))
     *       /ROICEN
#ifdef TRACERS_WATER
        DO N=1,NTM
          TRSIL(N,3)=((1.-ROICE)*TRO(N)*XSI(3)*MSI2NO/(ACE1I+MSI2NO)
     *                +ROICE*(TRSIL(N,3)
     *         -FTRSI3(N)))/ROICEN
          TRSIL(N,4)=((1.-ROICE)*TRO(N)*XSI(4)*MSI2NO/(ACE1I+MSI2NO)
     *                +ROICE*(TRSIL(N,4)
     *         +FTRSI3(N)+TRI(N)))/ROICEN
        END DO
#endif
        ROICE = ROICEN          ! new ice concentration
      END IF

C**** COMPRESS THE ICE HORIZONTALLY IF TOO THIN OR LEAD FRAC. TOO SMALL
      HAVG = ROICE*(ACE1I+MSI2)/RHOI ! average ice thickness in meters
      ! using exp() form to reduce wintertime open ocean heat loss in the
      ! presence of thick ice, keeping same lead fraction at 1 m thickness
      !OPNOCN=MIN(0.1d0,FLEAD/HAVG)
      OPNOCN=MIN(0.0D0,FLEAD*EXP(-BYHREF*(HAVG-1D0)))
      IF ((ROICE*(ACE1I+MSI2)).gt.FLEADMX*RHOI) OPNOCN=0. ! no leads for h>mx
      IF (MSI2.LT.AC2OIM .or. ROICE.GT.1.-OPNOCN) THEN

C**** separate out snow and ice components
      call get_snow_ice_layer(SNOW,MSI2,HSIL,SSIL,
#ifdef TRACERS_WATER
     *       TRSIL,TRSNOW,TRICE, 
#endif 
     *       SNOWL,HSNOW,HICE,SICE,TSNW,TSIL,MICE,.false.)

      ROICEN = MIN(ROICE*(ACE1I+MSI2)/(ACE1I+AC2OIM), 1.-OPNOCN)
      DRSI  = ROICEN - ROICE ! < 0 compressed ice concentration
      FMSI1 = - MICE(1) * DRSI / ROICEN
      FMSI2 = - (MICE(1)+MICE(2)) * DRSI / ROICEN
      FMSI3 = FMSI2 * XSI(4)
C     FMSI4 = 0
      FHSI1 = HICE(1) * FMSI1 / (MICE(1) + 1d-30)
      FHSI2 = HICE(2) * FMSI2 /  MICE(2)
      FHSI3 = HICE(3) * FMSI3 / (MSI2*XSI(3))
      FSSI1 = SICE(1) * FMSI1 / (MICE(1) + 1d-30)
      FSSI2 = SICE(2) * FMSI2 /  MICE(2)
      FSSI3 = SICE(3) * FMSI3 / (MSI2*XSI(3))
      HICE(2) = HICE(2)*ROICE/ROICEN + FHSI1-FHSI2
      HICE(3) = HICE(3)*ROICE/ROICEN + FHSI2-FHSI3
      HICE(4) = HICE(4)*ROICE/ROICEN + FHSI3
      SICE(2) = SICE(2)*ROICE/ROICEN + FSSI1-FSSI2
      SICE(3) = SICE(3)*ROICE/ROICEN + FSSI2-FSSI3
      SICE(4) = SICE(4)*ROICE/ROICEN + FSSI3
#ifdef TRACERS_WATER
      FTRSI1(:) = TRICE(:,1) * FMSI1 / (MICE(1) + 1d-30)
      FTRSI2(:) = TRICE(:,2) * FMSI2 /  MICE(2)
      FTRSI3(:) = TRICE(:,3) * FMSI3 / (MSI2*XSI(3))
      TRICE(:,2) = TRICE(:,2)*ROICE/ROICEN + FTRSI1(:)-FTRSI2(:)
      TRICE(:,3) = TRICE(:,3)*ROICE/ROICEN + FTRSI2(:)-FTRSI3(:)
      TRICE(:,4) = TRICE(:,4)*ROICE/ROICEN + FTRSI3(:)
      TRSNOW(:,:) = TRSNOW(:,:)*ROICE/ROICEN
#endif
      MSI2 = MSI2*ROICE/ROICEN + FMSI2
C     SNOW*ROICEN = SNOW*roice  ! snow mass (kg) is conserved
      SNOWL(:) = SNOWL(:)*ROICE/ROICEN
      HSNOW(:) = HSNOW(:)*ROICE/ROICEN
      ROICE = ROICEN
      
C**** relayer upper two layers
        call relayer_12(HSNOW,HICE,SICE,MICE,SNOWL
#ifdef TRACERS_WATER
     *       ,TRSNOW,TRICE 
#endif 
     *       )

C**** reconstitute snow and ice layers
      call set_snow_ice_layer(HSNOW,HICE,SICE,MICE,SNOWL,
#ifdef TRACERS_WATER
     *       TRSNOW,TRICE,TRSIL, 
#endif 
     *       SNOW,MSI1,MSI2xx,HSIL,SSIL)
      END IF
C****
      END IF

C**** Clean up ice fraction (if rsi>(1-OPNOCN)-1d-3) => rsi=(1-OPNOCN))
      IF (ROICE.gt.0) THEN
      HAVG = ROICE*(ACE1I+MSI2)/RHOI ! average ice thickness in meters
      !OPNOCN=MIN(0.1d0,FLEAD/HAVG)
      OPNOCN=MIN(0.0D0,FLEAD*EXP(-BYHREF*(HAVG-1D0)))
      IF ((ROICE*(ACE1I+MSI2)).gt.FLEADMX*RHOI) OPNOCN=0.  ! no leads for h>mx
      IF (ROICE.gt.(1.-OPNOCN)-1d-3) THEN
        ROICEN = 1.-OPNOCN
        DRSI = MAX(0d0,ROICEN-ROICE) ! +ve
        IF (DRSI.gt.0) THEN
        FMSI4 = (ACE1I+MSI2)*(DRSI/ROICEN) ! >0 ice mass out of layer 4
        FHSI4 = HSIL(4)*FMSI4/(XSI(4)*MSI2)
        FSSI4 = SSIL(4)*FMSI4/(XSI(4)*MSI2)
        FHSI3 = HSIL(3)*FMSI4/MSI2 ! downward heat flux into layer 4
        FSSI3 = SSIL(3)*FMSI4/MSI2 ! downward salt flux into layer 4
#ifdef TRACERS_WATER
        FTRSI4(:) = TRSIL(:,4)*FMSI4/(XSI(4)*MSI2)
        FTRSI3(:) = TRSIL(:,3)*FMSI4/MSI2
#endif
        MSI2=MSI2-FMSI4         ! new ice mass of second physical layer
        
C**** separate out snow and ice components
        call get_snow_ice_layer(SNOW,MSI2,HSIL,SSIL,
#ifdef TRACERS_WATER
     *       TRSIL,TRSNOW,TRICE, 
#endif 
     *       SNOWL,HSNOW,HICE,SICE,TSNW,TSIL,MICE,.false.)

C**** fractions of new ice in each layer (deals with thick snow case)
        FRI(1:2)=MICE(1:2)/(ACE1I+MSI2)
        FRI(3:4)=XSI(3:4)*MSI2/(ACE1I+MSI2)

C**** distribute snow variables over new ice extent
        SNOWL(:) = SNOWL(:)*(ROICE/ROICEN)
        HSNOW(:) = HSNOW(:)*(ROICE/ROICEN)
CC      SSNOW = SSNOW*(ROICE/ROICEN)  ! always zero
#ifdef TRACERS_WATER
        TRSNOW(:,1:2) = TRSNOW(:,1:2)*(ROICE/ROICEN)
#endif

C**** Add new ice to ice variables
        HICE(1:2) = (ROICE/ROICEN)*(FHSI4*FRI(1:2)+HICE(1:2))
        SICE(1:2) = (ROICE/ROICEN)*(FSSI4*FRI(1:2)+SICE(1:2))
#ifdef TRACERS_WATER
        TRICE(:,1)=(ROICE/ROICEN)*(FTRSI4(:)*FRI(1)+TRICE(:,1))
        TRICE(:,2)=(ROICE/ROICEN)*(FTRSI4(:)*FRI(2)+TRICE(:,2))
#endif
        MICE(1:2)=(ROICE/ROICEN)*(FMSI4*FRI(1:2)+MICE(1:2))

C**** relayer upper two layers
        call relayer_12(HSNOW,HICE,SICE,MICE,SNOWL
#ifdef TRACERS_WATER
     *       ,TRSNOW,TRICE 
#endif 
     *       )

C**** reconstitute snow and ice layers
        call set_snow_ice_layer(HSNOW,HICE,SICE,MICE,SNOWL,
#ifdef TRACERS_WATER
     *       TRSNOW,TRICE,TRSIL, 
#endif 
     *       SNOW,MSI1,MSI2,HSIL,SSIL)

C**** lower layer adjustmensts
        HSIL(3) = (ROICE/ROICEN)*(HSIL(3)+FHSI4*FRI(3)-FHSI3)
        HSIL(4) = (ROICE/ROICEN)*(HSIL(4)+FHSI4*FRI(4)+FHSI3-FHSI4)
        SSIL(3) = (ROICE/ROICEN)*(SSIL(3)+FSSI4*FRI(3)-FSSI3)
        SSIL(4) = (ROICE/ROICEN)*(SSIL(4)+FSSI4*FRI(4)+FSSI3-FSSI4)
#ifdef TRACERS_WATER
        DO N=1,NTM
          TRSIL(N,3) = (ROICE/ROICEN)*(TRSIL(N,3)+FTRSI4(N)*FRI(3)
     *         -FTRSI3(N))
          TRSIL(N,4) = (ROICE/ROICEN)*(TRSIL(N,4)+FTRSI4(N)*FRI(4)
     *         +FTRSI3(N)-FTRSI4(N))
        END DO
#endif
        ROICE = ROICEN          ! new ice concentration
        END IF
      END IF
      END IF

      ELSE
C**** Ensure that MSI2 does not get too small for fixed-SST case.
        IF (ROICE.gt.0. and. MSI2.lt.AC2OIM) then
          DMIMP=AC2OIM-MSI2
          DHIMP=SUM(HSIL(3:4)*XSI(3:4)*DMIMP)
          DSIMP=SUM(SSIL(3:4)*XSI(3:4)*DMIMP)
          HSIL(3:4)=HSIL(3:4)*AC2OIM/MSI2
          SSIL(3:4)=SSIL(3:4)*AC2OIM/MSI2
#ifdef TRACERS_WATER
          DTRIMP(:)=(TRSIL(:,3)*XSI(3)+TRSIL(:,4)*XSI(4))*DMIMP
          TRSIL(:,3:4)=TRSIL(:,3:4)*AC2OIM/MSI2
#endif
          MSI2=AC2OIM
       END IF
      END IF
C**** Calculate temperatures for diagnostics and radiation
      call tice(hsil,ssil,msi1,msi2,tsil)

      RETURN
      END SUBROUTINE ADDICE

      SUBROUTINE SIMELT(DT,ROICE,SNOW,MSI2,HSIL,SSIL,POCEAN,Tm,TFO,TSIL,
#ifdef TRACERS_WATER
     &     TRSIL,TRUN0,
#endif
     &     ENRGMAX,ENRGUSED,RUN0,SALT)
!@sum  SIMELT melts sea ice laterally and if it is too small
!@+    Note: all amounts are with respect to the ocean/lake fraction
!@auth Original Development Team
!@ver  1.1
      IMPLICIT NONE
!@var POCEAN ocean fraction (zero if lake)
      REAL*8, INTENT(IN) :: POCEAN
!@var TFO freezing temperature of water (C)
      REAL*8, INTENT(IN) :: TFO
!@var ROICE,SNOW,MSI2 ice variables (%,kg/m^2,kg/m^2)
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
!@var Tm mixed layer ocean temperature (C)
      REAL*8, INTENT(IN) :: Tm
!@var DT time step (s)
      REAL*8, INTENT(IN) :: DT
!@var HSIL ice enthalpy  (J/m^2)
!@var SSIL ice salt (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL, SSIL
!@var TSIL ice temperature (C)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var ENRGMAX max energy available for melting ice (J/m^2)
      REAL*8, INTENT(IN) :: ENRGMAX
!@var ENRGUSED energy used to melt ice (J/m^2)
      REAL*8, INTENT(OUT) :: ENRGUSED
!@var RUN0 amount of sea ice melt (kg/m^2)
!@var SALT amount of salt in sea ice melt (kg/m^2)
      REAL*8, INTENT(OUT) :: RUN0, SALT
!@var DRSI change in RSI fraction
!@var DTEMP temperature diff. used in lateral melt calculation.
      REAL*8 :: DRSI,DTEMP
#ifdef TRACERS_WATER
!@var TRSIL tracer ice amounts (kg/m^2)
!@var TRUN0 tracer amount in sea ice melt (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: TRUN0
#endif
      INTEGER L

C**** Estimate DRSI
      DRSI=0.
      IF (ROICE.lt.1d-3) THEN
        DRSI=ROICE
      ELSE      ! IF (POCEAN.gt.0) THEN ! now for lakes too
C**** Estimate lateral melt using parameterisation from Makyut/Steele
C**** (via C. Bitz): Rside=dt*pi/(floesize*eta)*(3d-6)*(delT)^(1.36)
        dtemp=MAX(Tm-TFO,0d0)
        DRSI=DT*SILMFAC*dtemp**SILMPOW
        IF (ROICE-DRSI < 1d-3) DRSI=ROICE
        IF (ENRGMAX+DRSI*SUM(HSIL).lt.0) DRSI=-ENRGMAX/SUM(HSIL)
        IF (ROICE-DRSI > 1)    DRSI = 1 - ROICE
      END IF
C**** Remove DRSI amount of ice
      ENRGUSED=-DRSI*SUM(HSIL) !  [J/m^2]
      RUN0=DRSI*(SNOW + ACE1I + MSI2)
      SALT=DRSI*SUM(SSIL)
#ifdef TRACERS_WATER
      TRUN0(:)=DRSI*(TRSIL(:,1)+TRSIL(:,2)+TRSIL(:,3)+TRSIL(:,4))
#endif
      ROICE=min(1d0,ROICE-DRSI) ! possible round off error
      IF (ROICE.lt.1d-10) THEN
        ROICE=0.                ! deal with possible round off err
C**** set defaults if no ice is left
        SNOW=0.
        MSI2=AC2OIM
        IF (POCEAN.gt.0) THEN
          SSIL(1:2)=SSI0*XSI(1:2)*ACE1I
          SSIL(3:LMI)=SSI0*XSI(3:LMI)*AC2OIM
        ELSE
          SSIL(:) = 0.
        END IF
        DO L=1,2
          HSIL(L)=(XSI(L)*ACE1I)*Ei(TFO,1d3*SSIL(L)/(XSI(L)*ACE1I))
        END DO
        DO L=3,LMI
          HSIL(L)=(XSI(L)*AC2OIM)*Ei(TFO,1d3*SSIL(L)/(XSI(L)*AC2OIM))
        END DO

        TSIL=TFO

        IF(TFO.GT.Ti(0d0,1d3*SSI0).AND.TFO.NE.0.) HSIL=0.

#ifdef TRACERS_WATER
        TRSIL(:,:)=0.
#endif
      END IF
C****
      RETURN
      END SUBROUTINE SIMELT

      SUBROUTINE SSIDEC(SNOW,MSI2,HSIL,SSIL,DT,MELT12,
#ifdef TRACERS_WATER
     &     TRSIL,TRFLUX,
#endif
     &     MFLUX,HFLUX,SFLUX,debug)
!@sum  SSIDEC decays salinity in sea ice
!@auth Jiping Liu
      IMPLICIT NONE

!@var HSIL enthalpy of ice layers (J/m^2)
!@var SSIL salt in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: SSIL,HSIL
!@var MSI2 second mass layer ice thickness (kg/m^2)
!@var MSI1 first mass layer ice thickness (kg/m^2)
!@var DT source time step (s)
!@var MELT12 surface melt this time step (kg/m2)
      REAL*8, INTENT(INOUT) :: MSI2, MELT12
      REAL*8 :: DT, SNOW
!@var MFLUX,SFLUX,HFLUX mass, salt and heat flux arising from
!@+   sea salinity decay
      REAL*8, INTENT(OUT) :: MFLUX,HFLUX,SFLUX
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: TRSIL
!@var TRFLUX tracer flux arising from sea salinity decay
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: TRFLUX
      REAL*8, DIMENSION(NTM,2) :: TRSNOW
      REAL*8, DIMENSION(NTM,LMI-1) :: FTRSI
      REAL*8, DIMENSION(NTM,LMI) :: DTRSI,TRICE
      INTEGER N
#endif
      LOGICAL, INTENT(IN) :: debug
      REAL*8 DSSI(LMI),DMSI(LMI),DHSI(LMI)
      REAL*8 FMSI(LMI-1),FSSI(LMI-1),FHSI(LMI-1)
      REAL*8 HSNOW(2),HICE(LMI),SICE(LMI),TSIL(LMI),MICE(LMI),MSI1
     $     ,SNOWL(2),TSNW(LMI) 
      INTEGER L
!@var dtssi decay time scale for sea ice salinity (days)
!@var bydtssi decay constant for sea ice salinity (1/s)
      REAL*8, parameter :: dtssi = 30.d0
      REAL*8, parameter :: bydtssi = 1./(dtssi*SECONDS_PER_DAY)
      REAL*8 :: rate,brine_frac

      DSSI = 0. ; DMSI = 0. ; DHSI = 0.
#ifdef TRACERS_WATER
      DTRSI(:,:)= 0.
#endif

C**** separate out snow and ice components
      call get_snow_ice_layer(SNOW,MSI2,HSIL,SSIL,
#ifdef TRACERS_WATER
     *     TRSIL,TRSNOW,TRICE, 
#endif 
     *     SNOWL,HSNOW,HICE,SICE,TSNW,TSIL,MICE,.true.)

C**** brine removal is very sensitive to seaice thermo formulation
      SELECT CASE (seaice_thermo)
      CASE ("SI")               ! salinity only affects mass
C**** calculate removal of excess salinity using simple decay
        DO L=1,LMI
          IF (SICE(L).GT.ssi0*MICE(L)) THEN
            DSSI(L) = (SICE(L)-MICE(L)*ssi0)*DT*BYDTSSI
            DMSI(L) = DSSI(L)
            DHSI(L) = -DMSI(L)*TSIL(L)*SHI
#ifdef TRACERS_WATER
C**** no tracer removed if pure salt is lost
c           DTRSI(:,L) = (DMSI(L)-DSSI(L))*TRICE(:,L)/(MICE(L)-SICE(L))
c           TRICE(:,L) = TRICE(:,L)-DTRSI(:,L)
#endif
            SICE(L)=SICE(L)-DSSI(L)
            HICE(L)=HICE(L)-DHSI(L)
            MICE(L)=MICE(L)-DMSI(L)
          END IF
        END DO

      CASE ("BP")               ! Brine pocket formulation
C**** calculate removal of excess salinity using flushing and brine pocket limit
        DO L=1,LMI
          IF (SICE(L).gt.0) THEN
            IF(1d3*SICE(L)/MICE(L).GT.1d-10) THEN
              brine_frac=-mu*1d3*(SICE(L)/TSIL(L))/MICE(L)
            ELSE
              brine_frac=0.
            ENDIF
C**** flushing (30% of MELT12 pushes out an equivalent mass of brine)
            IF(brine_frac.NE.0.) THEN
              rate = min(1d0,0.3d0*MELT12/(MICE(L)*brine_frac)) ! fractional loss
            ELSE
              rate = 1.
            ENDIF
C**** basic gravity drainage (3 day timescale)
            if (brine_frac.gt.0.01d0) rate =
     *           min(rate + DT*BYDTSSI*100.*(brine_frac-0.01d0),1d0)
C**** remove very small amounts of salt
            if (SICE(L).lt.ssimin*MICE(L)) rate=1d0
            DMSI(L) = rate*brine_frac*MICE(L)
            DSSI(L) = rate*SICE(L)
c            DHSI(L) = rate*brine_frac*MICE(L)*shw*TSIL(L)
            DHSI(L) = -rate*mu*1d3*SICE(L)*shw
#ifdef TRACERS_WATER
C**** tracers may be fractionated in the brine....(assume not for now)
            DTRSI(:,L) = (DMSI(L)-DSSI(L))*TRICE(:,L)/(MICE(L)-SICE(L))
            TRICE(:,L) = TRICE(:,L)-DTRSI(:,L)
#endif
            SICE(L)=max(0d0,SICE(L)-DSSI(L))
            HICE(L)=HICE(L)-DHSI(L)
            MICE(L)=MICE(L)-DMSI(L)

C**** add brine expulsed to melt ponds if snow is gone
            if (SNOW.eq.0. .and. L.le.2) MELT12=MELT12+DMSI(L)
          END IF
        END DO

      END SELECT

C**** Calculate fluxes required to rebalance layers
C**** Force no upward salt flux to permit zero salt in upper layers
C**** under flushing conditions.

C**** Mass/heat moves from layer 2 to 1
      FMSI(1) = -DMSI(1)  ! <= 0
      FSSI(1) = 0.   ! no salt up
      FHSI(1) = FMSI(1)*(TSIL(2)*shi-lhm) ! energy of pure ice
#ifdef TRACERS_WATER
      FTRSI(:,1) = FMSI(1)*TRICE(:,2)/(MICE(2)-SICE(2))
#endif
C**** Mass/heat/salt moves from layer 3 to 2
      FMSI(2) = -DMSI(1)-DMSI(2)
      FSSI(2) = 0.   ! no salt up
      FHSI(2) = FMSI(2)*(TSIL(3)*shi-lhm) ! energy of pure ice
#ifdef TRACERS_WATER
      FTRSI(:,2) = FMSI(2)*TRICE(:,3)/(MICE(3)-SICE(3))
#endif

C**** add fluxes in upper layers
      IF (DMSI(1)+DMSI(2).gt.0) THEN
        HICE(1) = HICE(1)-FHSI(1)
        SICE(1) = SICE(1)-FSSI(1)
        HICE(2) = HICE(2) + FHSI(1) - FHSI(2)
        SICE(2) = SICE(2) + FSSI(1) - FSSI(2)
        MICE(1) = MICE(1) -FMSI(1) 
        MICE(2) = MICE(2) +FMSI(1) - FMSI(2)
#ifdef TRACERS_WATER
        TRICE(:,1) = TRICE(:,1)- FTRSI(:,1)
        TRICE(:,2) = TRICE(:,2)+ FTRSI(:,1) - FTRSI(:,2)
#endif
      END IF

C**** reconstitute snow and ice layers
      call set_snow_ice_layer(HSNOW,HICE,SICE,MICE,SNOWL,
#ifdef TRACERS_WATER
     *     TRSNOW,TRICE,TRSIL, 
#endif 
     *     SNOW,MSI1,MSI2,HSIL,SSIL)

C**** Mass/heat moves between layers 3 to 4
      FMSI(3) = XSI(3)*DMSI(4)+XSI(4)*(FMSI(2)-DMSI(3))
      IF (FMSI(3).gt.0) THEN      ! downward flux to layer 4
       FHSI(3) = FMSI(3)*HSIL(3)/MICE(3)
       FSSI(3) = FMSI(3)*SSIL(3)/MICE(3)
#ifdef TRACERS_WATER
       FTRSI(:,3) = (FMSI(3)-FSSI(3))*TRSIL(:,3)/(MICE(3)-SICE(3)) 
#endif
      ELSE                      ! upward flux
       FHSI(3) = FMSI(3)*(TSIL(4)*shi-lhm)    ! energy of pure ice only
       FSSI(3) = 0.                           ! no salt up
#ifdef TRACERS_WATER
       FTRSI(:,3) = FMSI(3)*TRSIL(:,4)/(MICE(4)-SICE(4))
#endif
      END IF

C**** Apply the second mass layer fluxes
      HSIL(3)=HSIL(3)+ FHSI(2)-FHSI(3)
      HSIL(4)=HSIL(4)+ FHSI(3)
      SSIL(3)=SSIL(3)+ FSSI(2)-FSSI(3)
      SSIL(4)=SSIL(4)+ FSSI(3)
c      if (SSIL(3).lt.ssimin*MICE(3)) SSIL(3)=0.
c      if (SSIL(4).lt.ssimin*MICE(4)) SSIL(4)=0.
#ifdef TRACERS_WATER
C**** Apply the tracer fluxes
      TRSIL(:,3)=TRSIL(:,3)+ FTRSI(:,2)-FTRSI(:,3)
      TRSIL(:,4)=TRSIL(:,4)+ FTRSI(:,3)
#endif
c     MSI1 = MSI1 - (DMSI(1)+DMSI(2)) - FMSI(2) ! stays fixed
      MSI2 = MSI2 + FMSI(2)

C**** output fluxes and diagnostics
      MFLUX = SUM(DMSI(1:LMI))         ! mass flux to ocean
      HFLUX = SUM(DHSI(1:LMI))         ! energy flux to ocean
      SFLUX = SUM(DSSI(1:LMI))         ! salt flux to ocean
#ifdef TRACERS_WATER
      DO N=1,NTM
        TRFLUX(N)=SUM(DTRSI(N,1:LMI)) ! tracer flux to ocean
      END DO
#endif
C****
      call tice(hsil,ssil,msi1,msi2,tsil)

      RETURN
      END SUBROUTINE SSIDEC

      subroutine iceocean_fluxes(Ti,Si,Tm,Sm,dh,ustar,Coriol,dtsrc,mlsh,
#ifdef TRACERS_WATER
     &     Tri,Trm,trflux,tralpha,
#endif
     &     mflux,sflux,hflux)
!@sum  iceocean calculates fluxes at base of sea ice
!@auth Gavin Schmidt
!@usage At the ice-ocean interface the following equations hold:
!@+                                              Tb = -mu Sb         (1)
!@+ -lam_i(Ti,Si)(Ti-Tb)/dh + rho_m shw g_T (Tb-Tm) = -m Lh(Tib,Sib)
!@+                                                  + m shw (Tib-Tb)(2)
!@+                           rho_m     g_S (Sb-Sm) =  m (Sib-Sb )   (3)
!@+  with  Sib=Si, Tib=Ti m>0,  or Sib = fsss Sb, Tib=Tb m<0        (4)
!@+ This routine iteratively solves for m, Tb, Sb and Sib, Tib using
!@+ Newton's method. The form of the latent heat and conductivity
!@+ can optionally include salinity effects.
!@+ The actual fluxes (positive down) are then:
!@+  mflux = m   ; sflux = 1d-3 m Sib (kg/m^2/s)
!@+  hflux = lam_i(Ti,Si) (Ti-Tb)/dh -m Lh(Tib,Sib)+m shw Tib  (J/m^2/s)
!@+
!@+  Formulations of salinty effects PI/SI/BP from Schmidt et al (2004)
!@+
      implicit none
!@var G_mole_T,G_mole_S molecular diffusion terms for heat/salt
C****  G_mole_T = 12.5 * 13.8d0**(2d0/3d0) - 6. = 65.9d0
C****  G_mole_S = 12.5 * 2432d0**(2d0/3d0) - 6. = 2255d0
      real*8, parameter :: G_mole_T = 65.9d0 , G_mole_S = 2255d0
!@var Si,Sm salinity in lowest ice layer and mixed layer (psu)
!@var Ti,Tm temperatures in ice and mixed layer (C)
!@var dh distance from center of bottom ice layer to base of ice (m)
      real*8, intent(in) :: Ti,Si,Tm,Sm,dh
!@var coriol Coriolis parameter (used in turbulent flux calc) (1/s)
!@var ustar friction velocity at ice-ocean interface (m/s)
      real*8, intent(in) :: ustar,Coriol
!@var dtsrc source time step (s)
!@var mlsh mixed layer specific heat capactity (J/m^2 C) (UNUSED)
      real*8, intent(in) :: dtsrc,mlsh
!@var mflux,sflux,hflux mass, salt and heat fluxes at base of ice
      real*8, intent(out) :: mflux,sflux,hflux
#ifdef TRACERS_WATER
!@var Tri,Trm tracer concentration in ice and mixed layer (kg/kg)
!@var tralpha tracer fraction going into ice (1)
      real*8, dimension(ntm), intent(in) :: Trm,Tri,tralpha
!@var trflux tracer mass flux at base
      real*8, dimension(ntm), intent(out) :: trflux
!@var Trib tracer concentration at interface ice (kg/kg)
      real*8, dimension(ntm) :: Trib
#endif

!@var g_T,g_S turbulent exchange velocities (m/s)
!@var G_turb turbulent diffusion term
!@var Sb,Sb0 final and initial estimate for salinity at interface (psu)
      real*8 :: G_turb,g_T,g_S,Sb,Sb0
!@var dSbdSb differential of Sb with respect to initial Sb0
      real*8 :: dSbdSb
!@var m melt rate (negative implies freezing) (kg/m^2/s)
!@var Tb temperature at interface (C)
!@var Tib,Sib actual temperature/salinity in melting/freezing ice (psu)
      real*8 :: m,Tb,Sib  ! ,Tib
      real*8 :: lh,left2,df3dm,dmdTb,dmdSi,alamdh,f0,df,rsg
      integer, parameter :: niter=5 !@param niter number of iterations
      integer :: i

C**** Calculate turbulent exchange velocities g_T, g_S

C**** G_turb = (1/k) ln (u* E n^2 / f h) + 1 / (2 E n) - (1/k)
C****        = 2.5 ln ( 5300*(u*)^2/coriol ) + 9.62 - 2.5 (assuming n=1)
C**** Note: n should theoretically depend on the buoyancy flux and
C****       therfore this should be included in the iteration below.
      G_turb = 2.5d0 * log ( 5300.*ustar*ustar/coriol ) + 7.12d0
C****  g  = u* / ( G_turb + G_mole )
      g_T = ustar / ( G_turb + G_mole_T )
      g_S = ustar / ( G_turb + G_mole_S )

C**** set conductivity term
C**** Diffusive flux is implicit in ice !  + ml temperature (not used)
      rsg=rhows*shw*g_T    ! /(1.+alpha*dtsrc*rhows*shw*g_T/mlsh)

      if (seaice_thermo.eq."SI" .or. Si.eq.0.) then
                                ! no salinity effects on diffusion
        alamdh = alami0/(dh+alpha*dtsrc*alami0*byshi/(2.*dh*rhoi))
      else                      ! brine pocket impact on diffusion
        alamdh = alami(Ti,Si)/(dh+alpha*dtsrc*
     *           alami(Ti,Si)/(dEidTi(Ti,Si)*2.*rhoi*dh))
      end if

C**** solve for boundary values (uses Newton-Rapheson with bounds)
C**** Estimate initial boundary salinity
      Sb0 = 0.75d0*Sm+0.25d0*Si
      do i=1,niter
C**** we use the full tfrez calc here, but assume for simplicity that
C**** the gradient w.r.t. T is still = -mu below.
        Tb = tfrez(Sb0)   ! -mu*Sb0  ! freezing point at salinity Sb0
C**** calculate left hand side of equation 2
        left2 = -alamdh*(Ti-Tb) + rsg*(Tb-Tm)
C**** depending on whether there is freezing or melting, the sea ice
C**** salinity in the melt/freeze fraction is set to be the original
C**** ice salinity or is a function of the boundary salinity.
C**** Similarly for temperature Tib
        if (left2.gt.0) then    ! freezing
          if (qsfix) then   ! keep salinity in ice constant
            Sib = ssi0*1d3
          else              ! it is a function of boundary value
            Sib = Sb0*fsss
          end if

          select case (seaice_thermo)
c          case ("PI")           ! pure ice
c            lh = lhm + Tb*(shw-shi)
          case ("SI")           ! salinity effects only mass
            lh = lhm*(1.-Sib*1d-3) + Tb*(shw-shi)
          case ("BP")           ! brine pockets
            if (Sib.gt.0) then
              lh = lhm*(1.+mu*Sib/Tb) + (Tb+mu*Sib)*(shw-shi)
            else
              lh = lhm + Tb*(shw-shi)
            end if
          end select
c
          m = -left2/lh

          select case (seaice_thermo)
c          case ("PI")           ! pure ice
c            dmdTb = (left2*(shw-shi)-lh*(alamdh + rsg))/(lh*lh)
c            dmdSi = 0.
          case ("SI")           ! salinity effects only mass
            dmdTb = (left2*(shw-shi)-lh*(alamdh + rsg))/(lh*lh)
            dmdSi = -left2*lhm*1d-3/(lh*lh)
          case ("BP")           ! brine pockets
            if (Sib.gt.0) then
              dmdTb = (left2*(-lhm*mu*Sib/Tb**2+shw-shi)-lh*(alamdh+rsg)
     *             )/(lh*lh)
              dmdSi = (left2*mu*(lhm/Tb+shw-shi))/(lh*lh)
            else
              dmdTb = (left2*(shw-shi)-lh*(alamdh + rsg))/(lh*lh)
              dmdSi = 0.
            end if
          end select

          if (qsfix) then   ! keep salinity in ice constant
            Sb = (rhows*g_S*Sm+m*Sib)/(rhows*g_S+m)
            df3dm= rhows*g_S*(Sib-Sm)/(rhows*g_S+m)**2
            dSbdSb= -mu*dmdTb*df3dm
          else              ! it is a function of boundary value
            Sb = rhows*g_S*Sm/(rhows*g_S+m*(1.-fsss))
            df3dm= -rhows*g_S*Sm*(1.-fsss)/(rhows*g_S+m*(1.-fsss))**2
            dSbdSb= (fsss*dmdSi-mu*dmdTb)*df3dm
          end if
        else                    ! melting
          Sib = Si

          select case (seaice_thermo)
c          case ("PI")           ! pure ice
c            lh = lhm + Tb*shw - Ti*shi
          case ("SI")           ! salinity effects only mass
            lh = lhm*(1.-Sib*1d-3) + Tb*shw - Ti*shi
          case ("BP")           ! brine pockets
            if (Sib.gt.0.AND.Ti.NE.0) then
              lh = lhm*(1.+mu*Sib/Ti)
     *            +(Ti+mu*Sib)*(shw-shi)-shw*(Ti-Tb)
            else
              lh = lhm + Tb*shw - Ti*shi
            end if
          end select
c
          m = -left2/lh

          Sb = (m*Sib+rhows*g_S*Sm)/(rhows*g_S+m)
          df3dm= (Sib*(rhows*g_S+m)-(m*Sib+rhows*g_S*Sm))/(rhows*g_S+m)
     *         **2
          dmdTb = -(alamdh + rsg)/lh + left2*shw/lh**2
          dSbdSb= -mu*dmdTb*df3dm
        end if

        f0=Sb-Sb0
        df=dSbdSb-1.+1d-20
        Sb = min(max(Sib,Sb0 - f0/df),40d0)
        Sb0=Sb
      end do
C**** define fluxes (positive down)
C**** Cap mass flux at at 90% of bottom layer
      m=max(min(0.9d0*2.*dh*rhoi/dtsrc,m),-0.9d0*2.*dh*rhoi/dtsrc)

      mflux = m                       ! (kg/m^2 s)
      sflux = 1d-3*m*Sib              ! (kg/m^2 s)
      hflux = alamdh*(Ti-Tb) - m*lh +m*shw*Tb ! (J/m^2 s)
#ifdef TRACERS_WATER
C**** Tracers use salinity turbulent diffusion term
      if (m.gt.0) then
        Trib(:)=Tri(:)
      else
        Trib(:)=tralpha(:)*rhows*g_S*Trm(:)/(rhows*g_S+(m-sflux)*(1.
     *       -tralpha(:)))
      end if
      trflux(:) = (m - sflux) * Trib(:)  ! (kg/m^2 s)
#endif
C****
      return
      end subroutine iceocean_fluxes

      subroutine icelake_fluxes(Ti,Tm,dh,dtsrc,mlsh,
#ifdef TRACERS_WATER
     &     Tri,Trm,trflux,tralpha,
#endif
     &     mflux,hflux)
!@sum  icelake calculates fluxes at base of lake ice (no salinity)
!@auth Gavin Schmidt
!@usage At the ice-lake interface the following equations hold:
!@+   interface is at freezing point (Tb=0.)   (1)
!@+   -lam_i Ti/dh - rho_m shw g_T Tm = -m Lh(Tib) + m shw Tib   (2)
!@+     with  Tib=Ti m>0,  or Tib=0. m<0        (4)
      implicit none
!@var mflux,hflux mass and heat fluxes at base of ice
      real*8, intent(out) :: mflux,hflux
!@var dtsrc source time step (s)
!@var Ti,Tm temperatures in ice and upper lake layer (C)
!@var dh distance from center of bottom ice layer to base of ice (m)
!@var mlsh mixed layer specific heat capactity (J/m^2 C) (UNUSED)
      real*8, intent(in) :: Ti,Tm,dh,dtsrc,mlsh
C**** Assume constant g_T = 1.3d-5, g_S = 0.025 * g_T m/s
!@var rsg = rhow * shw * g_T turbulent energy flux (J/m^2 K s)
!@var rgS = rhow * g_S turbulent tracer flux (kg/m^2 s)
      real*8, parameter ::  rsg = rhow*shw*1.3d-5, rgS=rhow*3.2d-7
      real*8 left2, lh, m, alamdh
#ifdef TRACERS_WATER
!@var Tri,Trm tracer concentration in ice and upper lake layer (kg/kg)
!@var tralpha tracer fraction going into ice (1)
      real*8, dimension(ntm), intent(in) :: Trm,Tri,tralpha
!@var trflux tracer mass flux at base
      real*8, dimension(ntm), intent(out) :: trflux
!@var Trib tracer concentration at interface ice (kg/kg)
      real*8, dimension(ntm) :: Trib
#endif

C**** Diffusive flux is implicit in ice temperature
      alamdh = alami0/(dh+alpha*dtsrc*byshi*alami0/(2.*dh*rhoi))
C**** calculate left hand side of equation 2
      left2 = -alamdh*Ti - rsg*Tm    !/(1.+alpha*dtsrc*rsg/mlsh)
      if (left2.gt.0) then      ! freezing
        lh = lhm
      else                      ! melting
        lh = lhm -  Ti*shi
      end if
C**** define fluxes (positive down)
      m = -left2/lh
C**** Cap mass flux at 90% of bottom layer
      m=max(min(0.9d0*2.*dh*rhoi/dtsrc,m),-0.9d0*2.*dh*rhoi/dtsrc)

      mflux = m                 ! (kg/m^2 s)
      hflux = alamdh*Ti - m*lh  ! (J/m^2 s)
#ifdef TRACERS_WATER
C**** Tracers use tracer turbulent diffusion term
      if (m.gt.0) then
        Trib(:) = Tri(:)
      else
        Trib(:) = tralpha(:)*rgS*Trm(:)/(rgS+m*(1.-tralpha(:)))
      end if
      trflux(:) = m*Trib(:)     ! (kg/m^2 s)
#endif
C****
      return
      end subroutine icelake_fluxes

      subroutine solar_ice_frac(snow,msi2,wetsnow,fsri,lmax)
!@sum solar_ice_frac calculates the fraction of solar radiation
!@+   transmitted through the ice, taking albedo into account, but
!@+   assuming constant proportions of visible and near-ir.
!@+   Based on Ebert et al, 1995.
!@auth Gavin Schmidt
      implicit none
!@var kiextvis,kiextnir1 ice extinction coeff. for vis and nir (1/m)
      real*8, parameter :: kiextvis=1.5d0, kiextnir1=18.d0
!@var snow,msi2 snow and second layer ice amounts (kg/m^2)
      real*8, intent(in) :: snow,msi2
!@var wetsnow =true if snow is wet
      logical, intent(in) :: wetsnow
!@var lmax number of ice layers to do calculation for
      integer, intent(in) :: lmax
!@var fsri fraction of solar energy reaching bottom of layer
      real*8, dimension(lmax), intent(out) :: fsri
      real*8, dimension(lmax) :: fsrvis, fsrnir1, hice
!@var ksextvis,ksextnir1 snow extinction coeff. for vis and nir (1/m)
!@var fracvis,fracnir fraction of absorbed solar in vis and nir
      real*8 fracvis,fracnir1,ksextvis,ksextnir1
      real*8 dsnow,dsnow1,hice12
      integer l

      dsnow = snow/rhos
      hice12= ace1i*byrhoi

c**** Set fractions of visible and near-ir bands based on weighting the
c**** solar input by the approximate co-albedo (ignoring melt ponds etc)
c**** VIS:  250 - 690 nm,       49.3% of incoming solar at ground
c**** NIR1: 690 - 1190 nm       34.9%       "       "   "    "
c**** NIR2/3 (> 1190 nm assumed not to be transmitted)
c****
      if (dsnow.gt.0.02) then    ! same cutoff as for albedo
c**** Extinction coefficients for snow depend on wetness
c**** Since albedo does as well, absorbed fractions also vary
        if (wetsnow) then
          fracvis =0.20d0 ;  fracnir1=0.33d0
          ksextvis=10.7d0 ; ksextnir1=118.d0
        else
          fracvis =0.06d0 ;  fracnir1=0.31d0
          ksextvis=19.6d0 ; ksextnir1=196.d0
        end if
      else
        fracvis=0.24d0 ; fracnir1=0.43d0
        ksextvis=10.7d0 ; ksextnir1=118.d0
      end if

C**** calculate sucessive fractions assuming each band is attenuated
C**** independently.
      if (ace1i*xsi(1).gt.snow*xsi(2)) then
c**** First thermal layer is snow and ice
        hice(1)    = (ace1i-xsi(2)*(snow+ace1i))*byrhoi
        fsrvis (1) = exp(-ksextvis *dsnow - kiextvis *hice(1))
        fsrnir1(1) = exp(-ksextnir1*dsnow - kiextnir1*hice(1))
        fsri(1)    = fracvis*fsrvis(1) + fracnir1*fsrnir1(1)
      else                      ! all snow
        dsnow1     = (ace1i+snow)*xsi(1)/rhos
        fsrvis (1) = exp(-ksextvis *dsnow1)
        fsrnir1(1) = exp(-ksextnir1*dsnow1)
        fsri(1)    = fracvis*fsrvis(1) + fracnir1*fsrnir1(1)
      end if
      fsrvis (2) = exp(-ksextvis *dsnow - kiextvis *hice12)
      fsrnir1(2) = exp(-ksextnir1*dsnow - kiextnir1*hice12)
      fsri(2)    = fracvis*fsrvis(2) + fracnir1*fsrnir1(2)

      do l=3,lmax   ! only done if lmax>2
        hice(l)    = xsi(l)*msi2*byrhoi
        fsrvis (l) = exp(-kiextvis *hice(l)) * fsrvis (l-1)
        fsrnir1(l) = exp(-kiextnir1*hice(l)) * fsrnir1(l-1)
        fsri(l)    = fracvis*fsrvis(l) + fracnir1*fsrnir1(l)
      end do
c****
      return
      end subroutine solar_ice_frac

      function tfrez(sss,press)
!@sum tfrez calculates freezing temperature of sea water
!@auth Gavin Schmidt
      implicit none
!@var sss sea surface salinity (psu)
      real*8, intent(in) :: sss
!@var press gauge pressure (default=0) (Pa)
      real*8, intent(in), optional :: press
!@var tfrez approx. freezing point of sea water (C)
      real*8 tfrez,pr
      real*8 :: a01 = -.0575d0, a02 = -2.154996d-4, a03 =1.710523d-3

      pr=0.
      if (present(press)) then
        pr=press
      end if
C**** linear approximation
c      tfrez = -mu*sss + dtdp*pr
C**** UNESCO formula (1983)
      tfrez = (a01 + a02*sss)*sss + a03*sss*sqrt(sss)+ dtdp*pr

      return
      end function tfrez

      subroutine snowice(Tm,Sm,SNOW,MSI2,HSIL,SSIL,qsfix,
#ifdef TRACERS_WATER
     *         Trm,Tralpha,TRSIL,TRSNWIC,
#endif
     *         MSNWIC,HSNWIC,SSNWIC,DSNOW)
      IMPLICIT NONE
!@var Tm ocean mixed layer temperature (C)
!@var Sm ocean mixed layer salinity (psu)
      REAL*8, INTENT(IN) :: Tm, Sm
!@var MSI1, MSI2 snow and ice mass in two mass layers (kg/m^2)
      REAL*8, INTENT(INOUT) :: SNOW,MSI2
!@var HSIL, SSIL energy and salt profiles in ice (J, kg)/m^2
      REAL*8, DIMENSION(LMI), INTENT(INOUT) :: HSIL,SSIL
!@var MSNWIC, HSNWIC, SSNWIC mass, energy and salt flux to ocean (<0)
      REAL*8, INTENT(OUT) :: MSNWIC, HSNWIC, SSNWIC
!@var DSNOW mass of snow changed to seaice (kg/m^2)
      REAL*8, INTENT(OUT) :: DSNOW
!@var qsfix true if sea ice salinity is fixed
      LOGICAL, INTENT(IN) :: qsfix
#ifdef TRACERS_WATER
!@var Trm tracer concentration in mixed layer
!@var Tralpha fractionation for snow ice
C**** Be careful to avoid taking too much tracer from ocean box
      REAL*8, DIMENSION(NTM), INTENT(IN) :: Trm, Tralpha
!@var TRSIL tracer profile in ice (kg/m^2)
      REAL*8, DIMENSION(NTM,LMI), INTENT(INOUT) :: TRSIL
!@var TRSNWIC tracer flux to ocean (<0)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRSNWIC
!@var Tri tracer concentration in snow ice
      REAL*8 :: Tri(NTM)   
      REAL*8 :: TRICE(NTM,LMI), TRSNOW(NTM,2)
#endif
      REAL*8 Z0,MAXM,MAXME,Eoc,Esnow1,Esnow2,Eic,Erat1
     $     ,Erat2,Si,MSI1,Tf,Eratd
      REAL*8 SICE(LMI),HICE(LMI),HSNOW(2),MICE(LMI),SNOWL(2),FMSI2,
     $     TSNW(2),TSIL(LMI)

C**** test for snow ice possibility
      IF (RHOI*SNOW.gt.(ACE1I+MSI2)*(RHOWS-RHOI)) THEN
        MSI1=SNOW+ACE1I
        Z0=(MSI1+MSI2)/RHOWS-(ACE1I+MSI2)/RHOI
        MAXM=Z0*RHOWS*(RHOI-RHOS)/(RHOWS+RHOS-RHOI)

C**** separate out snow and ice components
        call get_snow_ice_layer(SNOW,MSI2,HSIL,SSIL,
#ifdef TRACERS_WATER
     *       TRSIL,TRSNOW,TRICE, 
#endif 
     *       SNOWL,HSNOW,HICE,SICE,TSNW,TSIL,MICE,.true.)

C**** Check energy constraint on mass flux
C**** Si is the salinity only in the frozen sea water part
        if (qsfix) then
          Si=1d3*ssi0
        else
          Si=fsss*Sm
        end if
#ifdef TRACERS_WATER
C**** need to be careful with fractionation for tracers
        Tri(:)=Tralpha(:)*Trm(:)
#endif
        Tf=Tfrez(Sm)
        Eic=Ei(Tf,Si)
        Eoc=Tm*SHW
C**** There is a small error in Schmidt et al (2004) equation (12):
C**** (Esnow-Ei(Tm,Si))/(Ei(Tm,Si)-Eoc) should be
C**** (Esnow-Ei(Tm, 0))/(Ei(Tm,Si)-Eoc)
        Esnow1=HSNOW(1)/SNOWL(1)
        ERAT1=(Esnow1-Ei(Tf,0d0))/(Eic-Eoc)
        if (SNOWL(2).gt.0) then
          Esnow2=HSNOW(2)/SNOWL(2)
          ERAT2=(Esnow2-Ei(Tf,0d0))/(Eic-Eoc)
        else
          Esnow2=0.
          ERAT2=0.
        end if
        ERATD=(Esnow1-Esnow2)/(Eic-Eoc)
         
        MAXME=(Z0*RHOI*ERAT1-SNOWL(2)*ERATD)/(1.+ERAT1*(RHOWS-RHOI)
     $       /RHOWS)
        MAXM=MAX(0d0,MIN(MAXME,MAXM))   ! mass of sea water to be added
        MAXM=MAX(0d0,MIN(0.9*RHOWS*(MICE(2)/RHOI-Z0),MAXM)) ! limit for advection
        DSNOW=Z0*RHOI-MAXM*(RHOWS-RHOI)/RHOWS  ! total loss of snow
        IF (MAXM.eq.0.) DSNOW=MIN(0.9*MICE(2),DSNOW) ! limit for advection 

C**** distribute changes over ice and snow
        if (DSNOW.gt.SNOWL(2)) THEN ! all of snow layer 2 and some snow layer 1 go to ice
           HSNOW(1)=HSNOW(1)*(SNOWL(1)+SNOWL(2)-DSNOW)/SNOWL(1) ! temperature constant in remaining snow
           HSNOW(2)=0.
           MICE(1)=MICE(1)+MAXM+DSNOW
           SICE(1)=SICE(1)+MAXM*0.001d0*Si
           HICE(1)=HICE(1)+MAXM*Eoc+SNOWL(2)*Esnow2+(DSNOW-SNOWL(2))
     $          *Esnow1
#ifdef TRACERS_WATER
           TRICE(:,1)=TRICE(:,1)+MAXM*(1.-0.001d0*Si)*Tri(:)+TRSNOW(:,2)
     $          +(DSNOW-SNOWL(2))*TRSNOW(:,1)/SNOWL(1)
           TRSNOW(:,1)=TRSNOW(:,1)*(SNOWL(1)+SNOWL(2)-DSNOW)/SNOWL(1)
           TRSNOW(:,2)=0.
#endif
           SNOWL(1)=SNOWL(1)+SNOWL(2)-DSNOW
           SNOWL(2)=0.
        else
           HSNOW(2)=HSNOW(2)*(SNOWL(2)-DSNOW)/SNOWL(2)
           MICE(2)=MICE(2)+DSNOW+MAXM
           SICE(2)=SICE(2)+MAXM*0.001d0*Si
           HICE(2)=HICE(2)+MAXM*Eoc+DSNOW*Esnow2
#ifdef TRACERS_WATER
           TRICE(:,2)=TRICE(:,2)+MAXM*(1.-0.001d0*Si)*Tri(:)+DSNOW
     $          *TRSNOW(:,2)/SNOWL(2)
           TRSNOW(:,2)=TRSNOW(:,2)*(SNOWL(2)-DSNOW)/SNOWL(2)
#endif
           SNOWL(2)=SNOWL(2)-DSNOW
        endif

        FMSI2 = MAXM+DSNOW      ! mass of ice pushed down to layer 3

C**** relayer lower ice layers
        call relayer(FMSI2,MICE,HICE,SICE
#ifdef TRACERS_WATER
     *       ,TRICE
#endif
     *       )

C**** relayer upper two layers
        call relayer_12(HSNOW,HICE,SICE,MICE,SNOWL
#ifdef TRACERS_WATER
     *       ,TRSNOW,TRICE 
#endif 
     *       )

C**** reconstitute upper snow and ice layers
        call set_snow_ice_layer(HSNOW,HICE,SICE,MICE,SNOWL,
#ifdef TRACERS_WATER
     *       TRSNOW,TRICE,TRSIL, 
#endif 
     *       SNOW,MSI1,MSI2,HSIL,SSIL)

C**** output flux (positive down)
        MSNWIC = -MAXM
        HSNWIC = MSNWIC*Eoc
        SSNWIC = 0.001d0*Si*MSNWIC
#ifdef TRACERS_WATER
        TRSNWIC(:) = Tri(:)*MSNWIC*(1.-0.001d0*Si)
#endif
      ELSE
        MSNWIC = 0. ; HSNWIC = 0. ; SSNWIC = 0. ; DSNOW=0.
#ifdef TRACERS_WATER
        TRSNWIC = 0.
#endif
      END IF
      RETURN
      end subroutine snowice

      subroutine get_snow_ice_layer(SNOW,MSI2,HSIL,SSIL,
#ifdef TRACERS_WATER
     *    TRSIL,TRSNOW,TRICE, 
#endif 
     *    SNOWL,HSNOW,HICE,SICE,TSNW,TSIL,MICE,needtemp)
!@sum Split thermal layer fields into ice and snow components
!@+   in each thermal layer: i.e. SNOW in L1/L2, ICE in L1/L2
      REAL*8, INTENT(IN) :: SNOW,MSI2,HSIL(LMI),SSIL(LMI)
      REAL*8, INTENT(OUT) :: HSNOW(2),HICE(LMI),SICE(LMI),TSIL(LMI)
     *     ,MICE(LMI),SNOWL(2),TSNW(2)
      LOGICAL, INTENT(IN) :: needtemp
#ifdef TRACERS_WATER
      REAL*8, INTENT(IN) :: TRSIL(NTM,LMI)
      REAL*8, INTENT(OUT) :: TRSNOW(NTM,2),TRICE(NTM,LMI)
#endif 
      REAL*8 MSI1,Ti1,Si1
      INTEGER L

C**** Assume equal temperatures over snow and ice in single thermal layer
      MSI1=SNOW+ACE1I

      IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
        MICE(1) = ACE1I-XSI(2)*MSI1
        MICE(2) = XSI(2)*MSI1
        SNOWL(1)= SNOW
        SNOWL(2)= 0.

C**** Calculcate actual salinity (not bulk) and temperature of the ice in layer 1
        Si1 = 1d3*SSIL(1)/MICE(1)
        Ti1 = Ti2b(HSIL(1)/(XSI(1)*MSI1),Si1,SNOWL(1),MICE(1))

C**** We can’t have more energy in HICE than MICE*Em(Si) = -MICE*mu*Si*shw.
C**** All excess energy then has to be in HSNOW.
        HICE(1)  = min(max(MICE(1)*Ei(Ti1,Si1),HSIL(1)),MICE(1)*Em(Si1))

        HSNOW(1) = HSIL(1) - HICE(1)
        IF(SNOWL(1).EQ.0..OR.dabs(HSIL(1)-HICE(1)).LT.1.d-8) HSNOW(1)=0. ! clean up roundoff errors
        HSNOW(2) = 0.

        HICE(2) = HSIL(2)
        SICE(1) = SSIL(1)
        SICE(2) = SSIL(2)
#ifdef TRACERS_WATER
        TRSNOW(:,1) = TRSIL(:,1)*SNOWL(1)/(XSI(1)*MSI1-SSIL(1))
        TRSNOW(:,2) = 0.
        TRICE(:,1) = TRSIL(:,1)-TRSNOW(:,1)
        TRICE(:,2) = TRSIL(:,2)
#endif
      ELSE  ! some snow in second layer
        MICE(1) = 0.
        MICE(2) = ACE1I
        SNOWL(1)= XSI(1)*MSI1
        SNOWL(2)= XSI(2)*MSI1-ACE1I
        HSNOW(1) = HSIL(1)

C**** Calculcate actual salinity (not bulk) and temperature of the ice in layer 2
        Si1 = 1d3*SSIL(2)/MICE(2)
        Ti1 = Ti2b(HSIL(2)/(XSI(2)*MSI1),Si1,SNOWL(2),MICE(2))

        HICE(1) = 0.
C**** We can’t have more energy in HICE than MICE*Em(Si) = -MICE*mu*Si*shw.
C**** All excess energy then has to be in HSNOW.
        HICE(2)  = min(max(MICE(2)*Ei(Ti1,Si1),HSIL(2)),MICE(2)*Em(Si1))

        HSNOW(2) = HSIL(2) - HICE(2)

        IF(SNOWL(2).EQ.0..OR.dabs(HSIL(2)-HICE(2)).LT.1.d-8) HSNOW(2)=0. ! clean up roundoff errors

        SICE(1) = 0.
        SICE(2) = SSIL(2)
#ifdef TRACERS_WATER
        TRSNOW(:,1) = TRSIL(:,1)
        TRSNOW(:,2) = TRSIL(:,2)*SNOWL(2)/(XSI(2)*MSI1-SSIL(2))
        TRICE(:,1) = 0.
        TRICE(:,2) = TRSIL(:,2)-TRSNOW(:,2)
#endif
      END IF

C**** lower levels
      do L=3,LMI
        MICE(L)=XSI(L)*MSI2
        HICE(L)=HSIL(L)
        SICE(L)=SSIL(L)
#ifdef TRACERS_WATER
        TRICE(:,L)=TRSIL(:,L)
#endif
      end do

C**** snow/ice temperatures
      if (needtemp) then
        do L=1,2
          IF (SNOWL(L).gt.0) THEN
            TSNW(L) = Ti(HSNOW(L)/SNOWL(L),0d0)
          ELSE
            TSNW(L) = 0.
          END IF
        end do
        do L=1,LMI
          IF (MICE(L).gt.0) THEN
            TSIL(L) = Ti(HICE(L)/MICE(L),1d3*SICE(L)/MICE(L))
          ELSE
            TSIL(L)=0.
          END IF
        end do
      end if

      return
      end subroutine get_snow_ice_layer

      subroutine set_snow_ice_layer(HSNOW,HICE,SICE,MICE,SNOWL,
#ifdef TRACERS_WATER
     *    TRSNOW,TRICE,TRSIL, 
#endif 
     *    SNOW,MSI1,MSI2,HSIL,SSIL)
!@sum Collect snow and ice mass layer fields into thermal layers
!@+   SNOW in L1 and L2, ICE in L1 , ICE in L2
      REAL*8, INTENT(IN) :: HSNOW(2),HICE(LMI),SICE(LMI),MICE(LMI)
     $     ,SNOWL(2)
      REAL*8, INTENT(OUT) :: SNOW,MSI1,MSI2,HSIL(LMI),SSIL(LMI)
#ifdef TRACERS_WATER
      REAL*8, INTENT(IN) :: TRSNOW(NTM,2),TRICE(NTM,LMI)
      REAL*8, INTENT(OUT) :: TRSIL(NTM,LMI)
#endif 
      INTEGER L

      SNOW=SNOWL(1)+SNOWL(2)
      MSI1=SNOW+ACE1I
      MSI2=SUM(MICE(3:LMI))
      
      HSIL(1) = HSNOW(1) + HICE(1)
      SSIL(1) =            SICE(1)
      HSIL(2) = HSNOW(2) + HICE(2)
      SSIL(2) =            SICE(2)
#ifdef TRACERS_WATER
      TRSIL(:,1) = TRSNOW(:,1) + TRICE(:,1)
      TRSIL(:,2) = TRSNOW(:,2) + TRICE(:,2)
#endif

C**** lower levels
      do L=3,LMI
        HSIL(L)=HICE(L)
        SSIL(L)=SICE(L)
#ifdef TRACERS_WATER
        TRSIL(:,L)=TRICE(:,L)
#endif
      end do
      
      return
      end subroutine set_snow_ice_layer

      subroutine relayer_12(HSNOW,HICE,SICE,MICE,SNOWL
#ifdef TRACERS_WATER
     *    ,TRSNOW,TRICE 
#endif 
     *    )
!@sum Relayer upper 2 snow/ice layers if required
      REAL*8, INTENT(INOUT) :: HSNOW(2),HICE(LMI),SICE(LMI),MICE(LMI)
     $     ,SNOWL(2)
#ifdef TRACERS_WATER
      REAL*8, INTENT(INOUT) :: TRSNOW(NTM,2),TRICE(NTM,LMI)
      REAL*8 FTRSI1(NTM)
#endif
      REAL*8 FMSI1,FHSI1,FSSI1,FMSI0

      FMSI1 = SNOWL(1)+MICE(1)-XSI(1)*(SNOWL(1)+SNOWL(2)+ACE1I)

      IF(dabs(FMSI1).LT.1.d-10) FMSI1=0.

      FMSI0=XSI(1)*(SNOWL(1)+SNOWL(2)+ACE1I)
      IF(SNOWL(2).GT.0..AND.MICE(1)+MICE(2).GT.FMSI0) THEN
        IF(FMSI0.GT.MICE(2)) THEN
          FSSI1 = (FMSI0 - MICE(2))*SICE(1)/MICE(1)
          FHSI1 = (FMSI0 - MICE(2))*HICE(1)/MICE(1)
          SICE(1)=SICE(1) - FSSI1
          SICE(2)=SICE(2) + FSSI1
          HICE(1)=HICE(1) - FHSI1
          HICE(2)=HICE(2) + FHSI1
#ifdef TRACERS_WATER
          FTRSI1(:) = (FMSI0 - MICE(2))*TRICE(:,1)/MICE(1)
          TRICE(:,1)=TRICE(:,1)-FTRSI1(:)
          TRICE(:,2)=TRICE(:,2)+FTRSI1(:)
#endif
        ELSE
          FSSI1 = (MICE(2) - FMSI0)*SICE(2)/MICE(2)
          FHSI1 = (MICE(2) - FMSI0)*HICE(2)/MICE(2)
          SICE(1)=SICE(1) + FSSI1
          SICE(2)=SICE(2) - FSSI1
          HICE(1)=HICE(1) + FHSI1
          HICE(2)=HICE(2) - FHSI1
#ifdef TRACERS_WATER
          FTRSI1(:) = (MICE(2)-FMSI0)*TRICE(:,2)/MICE(2)
          TRICE(:,1)=TRICE(:,1)+FTRSI1(:)
          TRICE(:,2)=TRICE(:,2)-FTRSI1(:)
#endif
        ENDIF
        MICE(1)=MICE(1)+MICE(2)-FMSI0
        MICE(2)=FMSI0
        SNOWL(1)=SNOWL(1)+SNOWL(2)
        SNOWL(2)=0.
        HSNOW(1)=HSNOW(1)+HSNOW(2)
        HSNOW(2)=0.
#ifdef TRACERS_WATER
        TRSNOW(:,1)=TRSNOW(:,1)+TRSNOW(:,2)
        TRSNOW(:,2)=0.
#endif
      ELSEIF(MICE(1).GT.0..AND.MICE(1)+MICE(2).LT.FMSI0) THEN
        IF(FMSI0.LT.SNOWL(1)) THEN
          FHSI1 = (SNOWL(1)-FMSI0)*HSNOW(1)/SNOWL(1)
          HSNOW(2)=HSNOW(2)+FHSI1
          HSNOW(1)=HSNOW(1)-FHSI1
#ifdef TRACERS_WATER
          FTRSI1(:) = (SNOWL(1)-FMSI0)*TRSNOW(:,1)/SNOWL(1)
          TRSNOW(:,2)=TRSNOW(:,2)+FTRSI1(:)
          TRSNOW(:,1)=TRSNOW(:,1)-FTRSI1(:)
#endif
        ELSE
          FHSI1 = (FMSI0-SNOWL(1))*HSNOW(2)/SNOWL(2)
          HSNOW(2)=HSNOW(2)-FHSI1
          HSNOW(1)=HSNOW(1)+FHSI1
#ifdef TRACERS_WATER
          FTRSI1(:) = (FMSI0-SNOWL(1))*TRSNOW(:,2)/SNOWL(2)
          TRSNOW(:,2)=TRSNOW(:,2)-FTRSI1(:)
          TRSNOW(:,1)=TRSNOW(:,1)+FTRSI1(:)
#endif
        ENDIF
        SNOWL(2)=SNOWL(1)+SNOWL(2)-FMSI0
        SNOWL(1)=FMSI0
        HICE(2)=HICE(2)+HICE(1)
        HICE(1)=0.
        SICE(2)=SICE(2)+SICE(1)
        SICE(1)=0.
        MICE(2)=MICE(2)+MICE(1)
        MICE(1)=0.
#ifdef TRACERS_WATER
        TRICE(:,2) = TRICE(:,1)+TRICE(:,2)
        TRICE(:,1) = 0.
#endif
      ELSE

      IF (FMSI1.gt.0) THEN  ! flux from layer 1 to layer 2
        IF (MICE(1).gt.0) THEN ! flux ice and check for enough
          IF (FMSI1.gt.MICE(1)) THEN ! flux all ice and some snow  
#ifdef TRACERS_WATER
            FTRSI1(:) = (FMSI1-MICE(1))*TRSNOW(:,1)/SNOWL(1)
            TRSNOW(:,2)=TRSNOW(:,2)+FTRSI1(:)
            TRSNOW(:,1)=TRSNOW(:,1)-FTRSI1(:)
            TRICE(:,2) = TRICE(:,1)+TRICE(:,2)
            TRICE(:,1) = 0.
#endif 
            HSNOW(2)=HSNOW(2)+(FMSI1-MICE(1))*HSNOW(1)/SNOWL(1)
            HSNOW(1)=HSNOW(1)-(FMSI1-MICE(1))*HSNOW(1)/SNOWL(1)
            SNOWL(2)=SNOWL(2)+(FMSI1-MICE(1))
            SNOWL(1)=SNOWL(1)-(FMSI1-MICE(1))
            HICE(2)=HICE(2)+HICE(1)
            HICE(1)=0.
            SICE(2)=SICE(2)+SICE(1)
            SICE(1)=0.
            MICE(2)=MICE(2)+MICE(1)
            MICE(1)=0.
          ELSE  ! only flux ice
            FHSI1 = FMSI1*HICE(1)/MICE(1)
            FSSI1 = FMSI1*SICE(1)/MICE(1)
#ifdef TRACERS_WATER
            FTRSI1(:) = FMSI1*TRICE(:,1)/MICE(1)
            TRICE(:,1)=TRICE(:,1)-FTRSI1(:)
            TRICE(:,2)=TRICE(:,2)+FTRSI1(:)
#endif
            MICE(1)=MICE(1)-FMSI1
            MICE(2)=MICE(2)+FMSI1
            HICE(1)=HICE(1)-FHSI1
            HICE(2)=HICE(2)+FHSI1
            SICE(1)=SICE(1)-FSSI1
            SICE(2)=SICE(2)+FSSI1
          END IF
        ELSE ! flux snow 
          FHSI1 = FMSI1*HSNOW(1)/SNOWL(1)
c         FSSI1 = 0.
#ifdef TRACERS_WATER
          FTRSI1(:) = FMSI1*TRSNOW(:,1)/SNOWL(1)
          TRSNOW(:,1)=TRSNOW(:,1)-FTRSI1(:)
          TRSNOW(:,2)=TRSNOW(:,2)+FTRSI1(:)
#endif 
          SNOWL(1)=SNOWL(1)-FMSI1
          SNOWL(2)=SNOWL(2)+FMSI1
          HSNOW(1)=HSNOW(1)-FHSI1
          HSNOW(2)=HSNOW(2)+FHSI1
        END IF
      ELSE  ! < 0 up from layer 2 to layer 1
        IF (MICE(1).gt.0) THEN ! flux ice 
          FHSI1 = FMSI1*HICE(2)/MICE(2)
          FSSI1 = FMSI1*SICE(2)/MICE(2)
#ifdef TRACERS_WATER
          FTRSI1(:) = FMSI1*TRICE(:,2)/MICE(2)
          TRICE(:,1)=TRICE(:,1)-FTRSI1(:)
          TRICE(:,2)=TRICE(:,2)+FTRSI1(:)
#endif 
          MICE(1)=MICE(1)-FMSI1
          MICE(2)=MICE(2)+FMSI1
          HICE(1)=HICE(1)-FHSI1
          HICE(2)=HICE(2)+FHSI1
          SICE(1)=SICE(1)-FSSI1
          SICE(2)=SICE(2)+FSSI1
        ELSE ! flux snow and check if enough
          IF (SNOWL(2)+FMSI1.lt.0) THEN ! not enough snow
#ifdef TRACERS_WATER
            TRICE(:,1)= -(SNOWL(2)+FMSI1)*TRICE(:,2)/MICE(2)
            TRICE(:,2)=TRICE(:,2)-TRICE(:,1)
            TRSNOW(:,1)=TRSNOW(:,1)+TRSNOW(:,2)
            TRSNOW(:,2)=0.
#endif
            HICE(1)= -(SNOWL(2)+FMSI1)*HICE(2)/MICE(2)
            HICE(2)= HICE(2)-HICE(1)
            SICE(1)= -(SNOWL(2)+FMSI1)*SICE(2)/MICE(2)
            SICE(2)= SICE(2)-SICE(1)
            MICE(1)= -(SNOWL(2)+FMSI1)
            MICE(2)= MICE(2)-MICE(1)
            HSNOW(1)=HSNOW(1)+HSNOW(2)
            HSNOW(2)=0.
            SNOWL(1)=SNOWL(1)+SNOWL(2)
            SNOWL(2)=0.
          ELSE ! only flux snow
            FHSI1 = FMSI1*HSNOW(2)/SNOWL(2)
c           FSSI1 = 0.
#ifdef TRACERS_WATER
            FTRSI1(:) = FMSI1*TRSNOW(:,2)/SNOWL(2)
            TRSNOW(:,1)=TRSNOW(:,1)-FTRSI1(:)
            TRSNOW(:,2)=TRSNOW(:,2)+FTRSI1(:)
#endif 
            SNOWL(1)=SNOWL(1)-FMSI1
            SNOWL(2)=SNOWL(2)+FMSI1
            HSNOW(1)=HSNOW(1)-FHSI1
            HSNOW(2)=HSNOW(2)+FHSI1
          END IF
        END IF
      END IF
      
      END IF
      
      return
      end subroutine relayer_12

      subroutine relayer(FMSI2,MICE,HICE,SICE
#ifdef TRACERS_WATER
     *     ,TRICE
#endif
     *     )
!@sum relayer lower level ice because of mass flux btw layer 2 and 3
!@+   and possible changes in mass in each layer
      REAL*8, INTENT(IN) :: FMSI2
      REAL*8, INTENT(INOUT) :: MICE(LMI),HICE(LMI),SICE(LMI)
#ifdef TRACERS_WATER
      REAL*8, INTENT(INOUT) :: TRICE(NTM,LMI)
      REAL*8 :: FTRSI(NTM,LMI)
#endif
      REAL*8 :: FMSI(LMI),FHSI(LMI),FSSI(LMI),MSI2
      INTEGER L

      MSI2=SUM(MICE(3:LMI))

C**** Calculate consequential heat/salt flux between layers
      DO L=2,LMI-1
C**** Calculate mass/heat flux between layers L and L+1
        FMSI(L) = SUM(XSI(L+1:LMI))*(MSI2+FMSI2)-SUM(MICE(L+1:LMI))
        IF (FMSI(L).gt.0) THEN  ! downward flux to layer L+1
          IF (FMSI(L).gt.MICE(L)) THEN ! exceptional case
            FHSI(L) = HICE(L)+(FMSI(L)-MICE(L))*HICE(L-1)/MICE(L-1)
            FSSI(L) = SICE(L)+(FMSI(L)-MICE(L))*SICE(L-1)/MICE(L-1)
#ifdef TRACERS_WATER
            FTRSI(:,L) = TRICE(:,L)+(FMSI(L)-MICE(L))*TRICE(:,L-2)
     $           /MICE(L-1)
#endif
          ELSE
            FHSI(L) = HICE(L)*FMSI(L)/MICE(L)
            FSSI(L) = SICE(L)*FMSI(L)/MICE(L)
#ifdef TRACERS_WATER
            FTRSI(:,L) = TRICE(:,L)*FMSI(L)/MICE(L)
#endif
          END IF
        ELSE                    ! upward flux
          FHSI(L) = HICE(L+1)*FMSI(L)/MICE(L+1)
          FSSI(L) = SICE(L+1)*FMSI(L)/MICE(L+1)
#ifdef TRACERS_WATER
          FTRSI(:,L) = TRICE(:,L+1)*FMSI(L)/MICE(L+1)
#endif
        END IF
      END DO

C**** Adjust amounts resulting from fluxes.

      DO L=2,LMI
        IF (L.gt.2) THEN
          MICE(L)=MICE(L)+FMSI(L-1)
          HICE(L)=HICE(L)+FHSI(L-1)
          SICE(L)=SICE(L)+FSSI(L-1)
#ifdef TRACERS_WATER
          TRICE(:,L)=TRICE(:,L)+FTRSI(:,L-1)
#endif
        END IF
        IF (L.lt.LMI) THEN
c          IF (FMSI(L).gt.MICE(L)) THEN ! exceptional case 
c            MICE(L-1)=MICE(L-1)-FMSI(L)
c            HICE(L)=HICE(L-1)
c            HICE(L-1)=HICE(L-1)-FHSI(L)
c            SICE(L)=SICE(L)-FSSI(L)
c#ifdef TRACERS_WATER
c            TRICE(:,L)=TRICE(:,L)-FTRSI(:,L)
c#endif
c          ELSE
            MICE(L)=MICE(L)-FMSI(L)
            HICE(L)=HICE(L)-FHSI(L)
            SICE(L)=SICE(L)-FSSI(L)
#ifdef TRACERS_WATER
            TRICE(:,L)=TRICE(:,L)-FTRSI(:,L)
#endif
c          END IF
        END IF
      END DO

      return
      end subroutine relayer

      SUBROUTINE TICE(HSIL,SSIL,MSI1,MSI2,TSIL)
!@sum TICE returns array of ice temperatures from model variables
!@auth Gavin Schmidt
      IMPLICIT NONE
      REAL*8, DIMENSION(LMI), INTENT(IN) :: HSIL, SSIL  ! J/m2, kg/m2
      REAL*8, INTENT(IN) :: MSI1,MSI2  ! kg/m2
      REAL*8, DIMENSION(LMI), INTENT(OUT) :: TSIL  ! deg C
      REAL*8, DIMENSION(2) :: SNOWL,MICE
      INTEGER l

      IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
        MICE(1) = ACE1I-XSI(2)*MSI1
        MICE(2) = XSI(2)*MSI1
c       SNOWL(1)= SNOW
        SNOWL(1)= MSI1-ACE1I
        SNOWL(2)= 0.
      ELSE  ! some snow in second layer
        MICE(1) = 0.
        MICE(2) = ACE1I
        SNOWL(1)= XSI(1)*MSI1
        SNOWL(2)= XSI(2)*MSI1-ACE1I
      ENDIF

      IF(MICE(1).NE.0.) THEN
        TSIL(1)=Ti2b(HSIL(1)/(XSI(1)*MSI1),1d3*SSIL(1)/MICE(1),
     *             SNOWL(1),MICE(1))
      ELSE
        TSIL(1)=Ti(HSIL(1)/(XSI(1)*MSI1),0d0)
      ENDIF
        TSIL(2)=Ti2b(HSIL(2)/(XSI(2)*MSI1),1d3*SSIL(2)/MICE(2),
     *             SNOWL(2),MICE(2))
      do l=3,lmi
        TSIL(l)=Ti(HSIL(l)/(XSI(l)*MSI2),1d3*SSIL(l)/(XSI(l)*MSI2))
      end do

      RETURN
      END SUBROUTINE TICE

      REAL*8 FUNCTION Ti(Ei,Si)
!@sum Ti calculates sea ice temperature as a function of internal
!@+   energy and salinity content depending on ice thermo formulation
!@auth Jiping Liu/Gavin Schmidt
      USE CONSTANT, only : lhm,shw,shi,byshi
      IMPLICIT NONE
!@var Ei internal energy of sea ice (J/kg)
!@var Si salinity  of sea ice (10-3 kg/kg = ppt)
!@var TI temperature (deg C)
      REAL*8, INTENT(IN) :: Ei, Si
      real*8 b,c,det,tm

      select case (seaice_thermo)
c      case ("PI")               ! pure ice
cc**** solve Ei = shi*Ti-lhm
c        Ti=(Ei+lhm)*byshi
      case ("SI")               ! salinity mass effect
c**** solve Ei = shi*Ti-lhm*(1-1d-3*Si)
        Ti=(Ei+lhm*(1.-1d-3*Si))*byshi
      case ("BP")               ! Brine pocket formulation
        if (Si.gt.1d-10) then
c**** solve Ei = shi*(Ti+mu*Si)-lhm*(1+mu*Si/Ti)-mu*Si*shw
          Tm=-mu*Si
          if(Ei .ge. shw*Tm) then ! at or above melting point
            Ti=Tm
          else                  ! solve quadratic, take most negative root
            b=Tm*(shw-shi)-(Ei+lhm)
            c=lhm*Tm
            det=b*b-4d0*shi*c   ! > 0
            Ti=-5d-1*(b + sqrt(det))*byshi
          end if
        else
          Ti=(Ei+lhm)*byshi  ! pure ice case
          IF (dabs(Ei+lhm).LT.1.d-10) Ti=0. ! clean up roundoff errors
        end if
      end select

      return
      END FUNCTION Ti

      REAL*8 FUNCTION Ti2b(Eit,Si,SNOWL,MICE)
!@sum Ti2b calculates sea ice temperature as a function of internal
!@+   energy and salinity content depending on ice thermo formulation
!@+   when ice and snow are present in the same layer
!@auth Anthony Leboissetier/Gavin Schmidt
      USE CONSTANT, only : lhm,shw,shi,byshi
      IMPLICIT NONE
!@var Ei internal energy of sea ice (J/kg)
!@var Si salinity  of sea ice (10-3 kg/kg = ppt)
!@var SNOWL snow mass (kg/m^2) (use the actual salinity of the ice not the bulk salinity)
!@var MICE ice mass (kg/m^2)
!@var TI2b temperature (deg C)
      REAL*8, INTENT(IN) :: Eit, Si,SNOWL, MICE
      real*8 b,c,det,tm, Eil

      select case (seaice_thermo)
c      case ("PI")               ! pure ice
cc**** solve Ei = shi*Ti-lhm
c        Ti=(Ei+lhm)*byshi
      case ("SI")               ! salinity mass effect
c**** solve Ei = shi*Ti-lhm*(1-1d-3*Si)
c**** in this particular case, Si is the bulk salinity (like in the SI case for Ti)
c**** so it has to be multiplied by MICE/(MICE+SNOWL)
        Ti2b=(Eit+lhm*(1.-1d-3*Si*(MICE/(MICE+SNOWL))))*byshi
      case ("BP")               ! Brine pocket formulation
        if (Si.gt.1.d-10) then
c**** solve Ei = shi*(Ti+mu*Si)-lhm*(1+mu*Si/Ti)-mu*Si*shw
          Tm=-mu*Si
          b=(MICE/(MICE+SNOWL))*Tm*(shw-shi)-(Eit+lhm)
          c=(MICE/(MICE+SNOWL))*lhm*Tm
          det=b*b-4d0*shi*c   ! > 0
          Ti2b=-5d-1*(b + sqrt(det))*byshi
        else
          Ti2b=(Eit+lhm)*byshi  ! pure ice case
          IF (dabs(Eit+lhm).LT.1.d-10) Ti2b=0. ! clean up roundoff errors
        end if
      end select

      return
      END FUNCTION Ti2b

      REAL*8 FUNCTION Ei(Ti,Si)
!@sum Calculation of energy of ice (J/kg) as a function of temp and salinity
!@auth Gavin Schmidt
      USE CONSTANT, only : shi, lhm, shw
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: Ti,Si   ! deg C and psu

      select case (seaice_thermo)
c      case ("PI")               ! no salinity effect
c        Ei=Ti*shi-lhm
      case ("SI")               ! salinity affects only mass
        Ei=Ti*shi-lhm*(1.-1d-3*Si)
      case ("BP")               ! brine pocket formulation
        if (Si.gt.1d-10.AND.Ti.NE.0.) then  ! is this safe from T=0? (or T>-muS?
          Ei= (Ti+mu*Si)*shi-lhm*(1.+mu*Si/Ti)-shw*mu*Si
        else
          Ei= Ti*shi-lhm
        end if
      end select

      RETURN
      END FUNCTION Ei

      REAL*8 FUNCTION Mi(HSI,SSI,MSI)
!@sum Calculation of ice melt (kg/m2) as a function of hsi and salinity
!@auth Gavin Schmidt
      USE CONSTANT, only : lhm, shw, bylhm
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: HSI,SSI,MSI   ! J/m2, kg/m2, kg/m2

      select case (seaice_thermo)
c      case ("PI")               ! no salinity effect
c        Mi=max(0d0,msi+hsi*bylhm)
      case ("SI")               ! salinity affects only mass
        Mi=max(0d0,msi+hsi*bylhm/(1.-ssi/msi))
      case ("BP")               ! brine pocket formulation
        if (1d3*ssi/msi.GT.1d-10) then
          Mi=0.
          if (hsi+shw*mu*1d3*ssi.gt.0.OR.
     *        dabs(hsi+shw*mu*1d3*ssi).LT.1.d-12) Mi=msi
        else
          Mi=max(0d0,msi+hsi*bylhm)
        end if
      end select

      RETURN
      END FUNCTION Mi

      REAL*8 FUNCTION Fi(WAT,HSI,SSI,MSI)
!@sum Calculation of liquid water that can be frozen in ice (kg/m2) 
!@auth Gavin Schmidt
      USE CONSTANT, only : lhm, shw, bylhm
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: WAT,HSI,SSI,MSI   ! kg/m2,J/m2,kg/m2,kg/m2

      select case (seaice_thermo)
c      case ("PI")               ! no salinity effect
c        Fi=min(wat,max(-hsi*bylhm-msi,0d0))
      case ("SI")               ! salinity affects only mass
        Fi=min(wat,max(-hsi*bylhm-msi+ssi,0d0))
      case ("BP")               ! brine pocket formulation
        if (1d3*ssi/msi.GT.1d-10) then
           Fi=wat
        else
           Fi=min(wat,max(-hsi*bylhm-msi,0d0))
        end if
      end select

      RETURN
      END FUNCTION Fi


      REAL*8 FUNCTION Em(Si)
!@sum Calculation of the energy of ice melt (J/kg) as a function of salinity
!@auth Gavin Schmidt
      USE CONSTANT, only : shw
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: Si ! psu

      select case (seaice_thermo)
c      case ("PI")               ! no salinity effect
c        Em= 0.
      case ("SI")               ! salinity affects only mass
        Em= 0.
      case ("BP")               ! brine pocket formulation
        Em= -mu*Si*shw
      end select

      RETURN
      END FUNCTION Em

      REAL*8 FUNCTION Alami(Ti,Si)
!@sum Calculation of thermal diffusion of ice J/(m*degC*sec) =f(T,S)
!@+   (based on Pringle et al, 2007)
!@auth Gavin Schmidt
      USE CONSTANT, only : shi, lhm, shw
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: Ti,Si   ! deg C and psu

      if (seaice_thermo.eq."SI" .or. Si.lt.1d-10) then ! pure ice value
        alami=alami0
      else                      ! use brine fraction
        IF (Ti.NE.0.) THEN
        alami=alami0 + alamdT*Ti + alamdS*Si/Ti
        IF(alami.LE.0.) alami=alami0
        ELSE
        alami=alami0
        ENDIF
      end if

      RETURN
      END FUNCTION Alami

      REAL*8 FUNCTION dEidTi(Ti,Si)
!@sum Calculation of effective specific heat of ice J/(kg degC) =f(T,S)
!@auth Gavin Schmidt
      USE CONSTANT, only : shi, lhm
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: Ti,Si   ! deg C and psu

      if (seaice_thermo.eq."SI" .or. Si.lt.1d-10) then ! pure ice value
        dEidTi=shi
      else                      ! use brine fraction
        IF(Ti.NE.0.) THEN
        dEidTi=shi+lhm*mu*Si/(Ti*Ti)
        ELSE
        dEidTi=shi
        ENDIF
      end if

      RETURN
      END FUNCTION dEidTi

      REAL*8 FUNCTION dEidTiws(Ti,Si,SNOWL,MICE)
!@sum Calculation of effective specific heat of ice J/(kg degC) =f(T,S)
!@auth Gavin Schmidt
      USE CONSTANT, only : shi, lhm
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: Ti,Si,SNOWL,MICE   ! deg C and psu

      if (seaice_thermo.eq."SI" .or. Si.lt.1d-10) then ! pure ice value
        dEidTiws=shi
      else                      ! use brine fraction
        IF(Ti.NE.0..OR.MICE.EQ.0.) THEN
        dEidTiws=shi+(MICE/(MICE+SNOWL))*lhm*mu*Si/(Ti*Ti)
        ELSE
        dEidTiws=shi
        ENDIF
      end if

      RETURN 
      END FUNCTION dEidTiws

      END MODULE SEAICE

      MODULE SEAICE_COM
!@sum  SEAICE_COM contains the model arrays for seaice
!@auth Gavin Schmidt
      USE SEAICE, only : lmi
      USE EXCHANGE_TYPES, only : iceocn_xchng_vars
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      USE TIMESTREAM_MOD, only : timestream
      IMPLICIT NONE

      type icestate
!@var grid a pointer to the grid object whose domain bounds were
!@+   used to allocate an instance of this type
         type(dist_grid), pointer :: grid
         CHARACTER(LEN=8) :: DOMAIN
         INTEGER :: I_0,I_1, J_0,J_1
         INTEGER :: I_0H,I_1H, J_0H,J_1H
         LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
         INTEGER, DIMENSION(:), POINTER :: IMAXJ
         REAL*8, DIMENSION(:,:), POINTER ::
!@var FWATER water fraction of gridbox
     &        FWATER
!@var RSI fraction of water area covered in ice
     &        ,RSI
!@var SNOWI snow amount on sea ice (kg/m^2)
     &        ,SNOWI
!@var MSI mass of ice second layer (layer 1=const) (kg/m^2)
C**** Note that MSI includes the mass of salt in sea ice
     &        ,MSI
!@var ZSI total ice thickness (m)
     &        ,ZSI
!@var pond_melt amount of melt pond mass (kg/m^2)
C**** Note this is a virtual meltpond and is only used for
C**** albedo calculations
     &        ,pond_melt
!@var RSIX,RSIY first order moments of ice concentration (for advection)
     &        ,RSIX,RSIY
!@var RSISAVE saved value of sea ice concentration before DYNSI
     &        ,RSISAVE
!@var HSI enthalpy of each ice layer (J/m^2)
!@var SSI sea ice salt content (kg/m^2)
         REAL*8, DIMENSION(:,:,:), POINTER ::
     &        HSI,SSI
!@var flag_dsws true if snow on ice is wet (ie. rain or surface melt)
         LOGICAL, DIMENSION(:,:), POINTER ::
     &        flag_dsws
#ifdef TRACERS_WATER
         INTEGER :: ntm=0
!@var TRSI tracer amount in sea ice (kg/m^2)
         REAL*8, DIMENSION(:,:,:,:), POINTER :: TRSI
#endif
      end type icestate

!@var si_ocn ocean instance of ice variables
      type(icestate) :: si_ocn
!@var si_atm is the atmos-grid instance of seaice variables.
!@+   It contains lake ice variables in ocean-free atmos gridcells.
      type(icestate) :: si_atm

!@var iceocn derived-type strucure containing variables
!@+   (or pointers thereto) needed for seaice interactions with ocean.
      type(iceocn_xchng_vars) :: iceocn

      type(dist_grid), pointer :: sigrid

!@var RSIstream,ZSIstream interfaces for reading and time-interpolating
!@+   ice cover and ice thickness files when kocean=0
!@+   See general usage notes in timestream_mod.
!@+   Selection of time period occurs as per the comments concerning
!@+   seaice_yr in init_oceanice.
!@+   Note regarding comparison of results to runs that use traditional I/O:
!@+   if SICE datafiles contain monthly means and the piecewise parabolic
!@+   method is used for monthly->daily interpolation, the presence of OSST_eom
!@+   in the rundeck will prompt read_stream to read end-of-month values from
!@+   OSST_eom rather than computing them on the fly.  SICE and SICE_eom may
!@+   refer to the same file or directory.  The on-the-fly result will differ
!@+   due to roundoff effects.
      type(timestream) :: RSIstream,ZSIstream
!@var ZSI_exists flag indicating availability of an input file ZSI containing
!@+   ice thickness (only applicable to kocean=0 netcdf-input case).
!@var dm empirical factor relating ice thickness and ice fraction for kocean=0
!@+   (only applicable when ZSI_exists=.false.)
!@+   The presence of ZSIFAC in a rundeck indicates that the factor dm should
!@+   be read from that file and used; otherwise ZSI must be present in the
!@+   rundeck.
      logical :: ZSI_exists=.true.
!@var RSI_exists flag indicating availability of an input file SICE containing
!@+   ice cover
      logical :: RSI_exists=.true.
      real*8, allocatable, dimension(:,:) :: dm

      CONTAINS

      subroutine alloc_icestate_type(grd_dum,state,domain)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID,getDomainBounds
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(icestate) :: state
      character(len=*), intent(in) :: domain
      INTEGER :: I_0H,I_1H, J_0H,J_1H, I_0,I_1, J_0,J_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
      INTEGER :: IER
      INTEGER :: NTM

      call getDomainBounds(grd_dum,
     &     I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1,
     &     I_STRT_HALO=I_0H, I_STOP_HALO=I_1H,
     &     J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      if(trim(domain).ne.'OCEAN' .and.
     &   trim(domain).ne.'LAKES') then
        call stop_model('alloc_icestate_type: unrecognized domain',255)
      endif

      state % domain = trim(domain)
      state % I_0H = I_0H
      state % I_1H = I_1H
      state % J_0H = J_0H
      state % J_1H = J_1H

      state % I_0 = I_0
      state % I_1 = I_1
      state % J_0 = J_0
      state % J_1 = J_1

      state % HAVE_SOUTH_POLE = HAVE_SOUTH_POLE
      state % HAVE_NORTH_POLE = HAVE_NORTH_POLE

      ALLOCATE(
     &     state % IMAXJ(J_0H:J_1H),
     &     state % FWATER(I_0H:I_1H, J_0H:J_1H),
     &     state % RSI(I_0H:I_1H, J_0H:J_1H),
     &     state % SNOWI(I_0H:I_1H, J_0H:J_1H),
     &     state % MSI(I_0H:I_1H, J_0H:J_1H),
     &     state % ZSI(I_0H:I_1H, J_0H:J_1H),
     &     state % pond_melt(I_0H:I_1H, J_0H:J_1H),
     &     state % flag_dsws(I_0H:I_1H, J_0H:J_1H),
     &     state % HSI(LMI, I_0H:I_1H, J_0H:J_1H),
     &     state % SSI(LMI, I_0H:I_1H, J_0H:J_1H),
     &     state % RSIX(I_0H:I_1H, J_0H:J_1H),
     &     state % RSIY(I_0H:I_1H, J_0H:J_1H),
     &     state % RSISAVE(I_0H:I_1H, J_0H:J_1H),
     &     STAT=IER)

      state % IMAXJ(:) = I_1
      IF(HAVE_SOUTH_POLE) state%IMAXJ(J_0) = 1
      IF(HAVE_NORTH_POLE) state%IMAXJ(J_1) = 1

      state % rsi = 0.d0
      state % SNOWI = 0.d0
      state % MSI = 0.d0
      state % ZSI = 0.d0
      state % pond_melt = 0.d0
      state % flag_dsws  = .false.
      state % hsi = 0.d0
      state % ssi = 0.d0

#ifdef TRACERS_WATER
      NTM = state%NTM
      ALLOCATE( state % TRSI(NTM, LMI, I_0H:I_1H, J_0H:J_1H),
     *     STAT=IER)
      state % TRSI = 0.  ! default to prevent unecessary crash
#endif

      return
      end subroutine alloc_icestate_type

      END MODULE SEAICE_COM

      SUBROUTINE ALLOC_SEAICE_COM(grid) !(IM_in,JM_in)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run-time
!@auth Rodger Abel
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      USE SEAICE_COM, ONLY : iceocn, si_ocn, alloc_icestate_type
      USE EXCHANGE_TYPES, only : alloc_xchng_vars
      use seaice_com, only : dm
      IMPLICIT NONE
      TYPE(DIST_GRID), INTENT(IN) :: GRID
      !INTEGER, INTENT(IN) :: IM_in,JM_in

      !CALL INIT_GRID(sigrid,im_in,jm_in,1)

      call alloc_icestate_type(grid,si_ocn,'OCEAN')

#ifdef TRACERS_WATER
      iceocn%ntm = si_ocn%ntm
#endif
      call alloc_xchng_vars(grid,iceocn)
      deallocate(iceocn%rsi); iceocn%rsi => si_ocn%rsi

      allocate(dm(grid%i_strt_halo:grid%i_stop_halo,
     &            grid%j_strt_halo:grid%j_stop_halo))

      RETURN
      END SUBROUTINE ALLOC_SEAICE_COM

      subroutine read_seaice_ic
!@sum   read_seaice_ic read sea ice initial conditions file.
      use model_com, only : ioread
      use seaice_com, only : grid=>sigrid
      use pario, only : par_open,par_close,variable_exists
      use filemanager, only : file_exists
      implicit none
      integer :: fid
      logical :: altstyle_coldstart=.false.
      if(file_exists('SICEIC')) then
        fid = par_open(grid,'SICEIC','read')
        altstyle_coldstart = variable_exists(grid,fid,'tsi_top')
      elseif(file_exists('GIC')) then
        fid = par_open(grid,'GIC','read')
      else
        return
      endif
      if(altstyle_coldstart) then
        ! IC file has T,S
        call altstyle_seaice_ic(fid)
      else
        ! IC file has per-layer heat,salt
        call new_io_seaice (fid,ioread)
      endif
      call par_close(grid,fid)
      return
      end subroutine read_seaice_ic

      subroutine def_rsf_seaice(fid)
!@sum  def_rsf_seaice defines seaice array structure in restart files
!@auth M. Kelley
!@ver  beta
      use seaice_com, only : grid=>sigrid,si_ocn,iceocn
      use pario, only : defvar
      use domain_decomp_1d, only : getDomainBounds
      implicit none
      integer fid   !@var fid file id
      integer :: i_0h,i_1h, j_0h,j_1h
      real*8, dimension(:,:), allocatable :: arrdum
      call defvar(grid,fid,si_ocn%rsi,'rsi(dist_im,dist_jm)')
      call defvar(grid,fid,si_ocn%snowi,'snowi(dist_im,dist_jm)')
      call defvar(grid,fid,si_ocn%msi,'msi(dist_im,dist_jm)')
      call defvar(grid,fid,si_ocn%pond_melt,
     &     'pond_melt(dist_im,dist_jm)')
      call defvar(grid,fid,si_ocn%flag_dsws,
     &     'flag_dsws(dist_im,dist_jm)')
      call defvar(grid,fid,si_ocn%hsi,'hsi(lmi,dist_im,dist_jm)')
      call defvar(grid,fid,si_ocn%ssi,'ssi(lmi,dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid,fid,si_ocn%trsi,'trsi(ntm,lmi,dist_im,dist_jm)')
#endif
      call defvar(grid,fid,si_ocn%rsix,'rsix(dist_im,dist_jm)')
      call defvar(grid,fid,si_ocn%rsiy,'rsiy(dist_im,dist_jm)')
      call defvar(grid,fid,iceocn%dhsi,'dhsi(two,dist_im,dist_jm)')
      call defvar(grid,fid,iceocn%dmsi,'dmsi(two,dist_im,dist_jm)')
      call defvar(grid,fid,iceocn%dssi,'dssi(two,dist_im,dist_jm)')

      call getDomainBounds(grid, i_strt_halo=i_0h,i_stop_halo=i_1h,
     &               j_strt_halo=j_0h,j_stop_halo=j_1h)
      allocate(arrdum(i_0h:i_1h,j_0h:j_1h))
      call defvar(grid, fid, arrdum, 'wseai(dist_im,dist_jm)' )
      call defvar(grid, fid, arrdum, 'eseai(dist_im,dist_jm)' )
      call defvar(grid, fid, arrdum, 'sseai(dist_im,dist_jm)' )
      deallocate(arrdum)

      return
      end subroutine def_rsf_seaice

      subroutine new_io_seaice(fid,iaction)
!@sum  new_io_seaice read/write seaice arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use seaice_com, only : grid=>sigrid,si_ocn,iceocn
      use pario, only : write_dist_data,read_dist_data
      use domain_decomp_1d, only : getDomainBounds
      implicit none
      integer fid      !@var fid unit number of read/write
      integer iaction  !@var iaction flag for reading or writing to file
      integer :: i_0h,i_1h, j_0h,j_1h
      real*8, dimension(:,:), allocatable :: arrdum
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid, fid, 'rsi', si_ocn%rsi)
        call write_dist_data(grid, fid, 'snowi', si_ocn%snowi)
        call write_dist_data(grid, fid, 'msi', si_ocn%msi)
        call write_dist_data(grid, fid, 'pond_melt', si_ocn%pond_melt)
        call write_dist_data(grid, fid, 'flag_dsws', si_ocn%flag_dsws)
        call write_dist_data(grid, fid, 'hsi', si_ocn%hsi, jdim=3)
        call write_dist_data(grid, fid, 'ssi', si_ocn%ssi, jdim=3)
#ifdef TRACERS_WATER
        call write_dist_data(grid, fid, 'trsi', si_ocn%trsi, jdim=4)
#endif
        call write_dist_data(grid, fid, 'rsix', si_ocn%rsix)
        call write_dist_data(grid, fid, 'rsiy', si_ocn%rsiy)
        call write_dist_data(grid, fid, 'dhsi', iceocn%dhsi, jdim=3)
        call write_dist_data(grid, fid, 'dmsi', iceocn%dmsi, jdim=3)
        call write_dist_data(grid, fid, 'dssi', iceocn%dssi, jdim=3)

        call getDomainBounds(grid, i_strt_halo=i_0h,i_stop_halo=i_1h,
     &                 j_strt_halo=j_0h,j_stop_halo=j_1h)
        allocate(arrdum(i_0h:i_1h,j_0h:j_1h))
        call conserv_OMSI(arrdum)
        call write_dist_data( grid, fid, 'wseai', arrdum)
        call conserv_OHSI(arrdum)
        call write_dist_data( grid, fid, 'eseai', arrdum)
        call conserv_OSSI(arrdum)
        call write_dist_data( grid, fid, 'sseai', arrdum)
        deallocate(arrdum)
      case (ioread)             ! input from restart file
        call read_dist_data(grid, fid, 'rsi', si_ocn%rsi)
        call read_dist_data(grid, fid, 'snowi', si_ocn%snowi)
        call read_dist_data(grid, fid, 'msi', si_ocn%msi)
        call read_dist_data(grid, fid, 'pond_melt', si_ocn%pond_melt)
        call read_dist_data(grid, fid, 'flag_dsws', si_ocn%flag_dsws)
        call read_dist_data(grid, fid, 'hsi', si_ocn%hsi, jdim=3)
        call read_dist_data(grid, fid, 'ssi', si_ocn%ssi, jdim=3)
#ifdef TRACERS_WATER
        call read_dist_data(grid, fid, 'trsi', si_ocn%trsi, jdim=4)
#endif
        call read_dist_data(grid, fid, 'rsix', si_ocn%rsix)
        call read_dist_data(grid, fid, 'rsiy', si_ocn%rsiy)
        call read_dist_data(grid, fid, 'dhsi', iceocn%dhsi, jdim=3)
        call read_dist_data(grid, fid, 'dmsi', iceocn%dmsi, jdim=3)
        call read_dist_data(grid, fid, 'dssi', iceocn%dssi, jdim=3)
      end select
      return
      end subroutine new_io_seaice

      subroutine altstyle_seaice_ic(fid)
!@sum  altstyle_seaice_ic reads sea ice state from an alternate
!@+    style IC file with temperature,salinity (intensive units)
!@+    and converts to the per-layer heat,salt (extensive units
!@+    J/m2, kg/m2) provided by the traditional IC file.
!@+    Other fields (ice mass, meltponds etc.) in the alternate file
!@+    are the same as in the traditional file.
!@+    For the moment, temperature and salinity in the IC file are
!@+    specified via
!@+      1. tsnow_top (C), the temperature at the top of any
!@+         snow existing on the ice
!@+      2. tsi_top (C), the temperature at the top of the ice
!@+         (i.e. snow base)
!@+      3. salt (psu), the mean salinity of the ice
!@+    Temperature is piecewise linearly interpolated between
!@+      1. tsnow_top and tsi_top
!@+      2. tsi_top and the ocean temperature tfo
!@+    Before trying to understand this subroutine, see module
!@+    SEAICE for explanation of
!@+      1. the distinction between mass layers and thermal layers
!@+      2. the blending of ice and snow in thermal layers
!@+      3. ACE1I
!@auth M. Kelley
      use constant, only : tf
      use seaice_com, only : grid=>sigrid,si_ocn
      use pario, only : read_dist_data
      use domain_decomp_1d, only : getDomainBounds
      use seaice, only : lmi,Ei,ace1i,xsi
      implicit none
      integer fid      !@var fid unit number of read
      integer :: i_0h,i_1h, j_0h,j_1h
      real*8, dimension(:,:), allocatable ::
     &     tsi_top,tsnow_top,salt
      real*8, dimension(lmi) :: tl,ml
      real*8 :: mice,msnow,msi1,msi2,icefr(2),zbot,ztop,dtdz,tfo
      integer :: i,j,l
      
      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo
      allocate(tsi_top(i_0h:i_1h,j_0h:j_1h))
      allocate(tsnow_top(i_0h:i_1h,j_0h:j_1h))
      allocate(salt(i_0h:i_1h,j_0h:j_1h))

      call read_dist_data(grid, fid, 'rsi', si_ocn%rsi)
      call read_dist_data(grid, fid, 'snow', si_ocn%snowi)
      call read_dist_data(grid, fid, 'msi', si_ocn%msi)
      call read_dist_data(grid, fid, 'tsi_top', tsi_top)
      call read_dist_data(grid, fid, 'tsnow_top', tsnow_top)
      call read_dist_data(grid, fid, 'salt', salt)
      call read_dist_data(grid, fid, 'pond_melt', si_ocn%pond_melt)
      call read_dist_data(grid, fid, 'flag_dsws', si_ocn%flag_dsws)


      !tsi_top = tsi_top - tf
      !tsnow_top = tsnow_top - tf

      tfo = -1.9d0
      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,grid%i_stop
        if(si_ocn%rsi(i,j).le.0d0) then
          si_ocn%hsi(:,i,j) = 0d0
          si_ocn%ssi(:,i,j) = 0d0
          cycle
        endif
        msnow = si_ocn%snowi(i,j)
        mice = si_ocn%msi(i,j)
        msi2 = max(0d0,mice - ace1i)
        msi1 = si_ocn%snowi(i,j) + ace1i
        ml(1:2) = xsi(1:2)*msi1
        ml(3:lmi) = xsi(3:lmi)*msi2
        si_ocn%msi(i,j) = msi2

        icefr(1) = max(0d0,ace1i-ml(2))/ml(1)
        icefr(2) = min(1d0,ace1i/ml(2))

        dtdz = (tsi_top(i,j)-tfo)/(mice+1d-30)
        zbot = 0d0
        do l=1,lmi
          ztop = zbot
          if(l.le.2) then
            zbot = zbot + icefr(l)*ml(l)
          else
            zbot = zbot + ml(l)
          endif
          tl(l) = tsi_top(i,j) - dtdz*.5d0*(zbot+ztop)
        enddo

        dtdz = (tsnow_top(i,j)-tsi_top(i,j))/(msnow+1d-30)
        zbot = 0d0
        do l=1,2
          ztop = zbot
          zbot = zbot + (1d0-icefr(l))*ml(l)
          tl(l) = icefr(l)*tl(l) + (1d0-icefr(l))*
     &       (tsnow_top(i,j) - dtdz*.5d0*(zbot+ztop))
        enddo
        do l=1,2
          si_ocn%hsi(l,i,j) = Ei(tl(l),icefr(l)*salt(i,j))*ml(l)
          si_ocn%ssi(l,i,j) = (salt(i,j)*1d-3)*ml(l)*icefr(l)
        enddo
        do l=3,lmi
          si_ocn%hsi(l,i,j) = Ei(tl(l),salt(i,j))*ml(l)
          si_ocn%ssi(l,i,j) = (salt(i,j)*1d-3)*ml(l)
        enddo

      enddo
      enddo

      return
      end subroutine altstyle_seaice_ic
