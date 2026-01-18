#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
      MODULE LANDICE
!@sum  LANDICE contains variables and routines dealing with land ice
!@auth Original Development Team
!@ver  2010/10/13
!@cont PRECLI,LNDICE
      USE CONSTANT, only : lhm,rhoi,shi
#ifdef TRACERS_WATER
      USE TRACER_COM, only : NTM
#endif
      IMPLICIT NONE
      SAVE
!@var Z1E,Z2LI first and second layer thicknesses for land ice
      REAL*8, PARAMETER :: Z1E=.1d0, Z2LI=2.9d0
!@var ACE1LI first layer land ice mass (kg/m^2)
      REAL*8, PARAMETER :: ACE1LI = Z1E*RHOI
!@var HC1LI heat capacity of first layer land ice (J/m^2)
      REAL*8, PARAMETER :: HC1LI = ACE1LI*SHI
!@var ACE2LI second layer land ice mass (kg/m^2)
      REAL*8, PARAMETER :: ACE2LI = Z2LI*RHOI
!@var HC2LI heat capacity of second layer land ice (J/m^2)
      REAL*8, PARAMETER :: HC2LI = ACE2LI*SHI
!@param SNMIN snow minimum used in some tracer calculations (kg/m2)
      REAL*8, PARAMETER :: SNMIN = 1d-5
!@dbparam glmelt_on determines whether glacial melt is used for oceans
      INTEGER :: glmelt_on = 1 ! yes w/ann adjustment (2: yes but fixed)
!@dbparam glmelt_fac_nh is a factor to multiply glacial melt by in NH
      REAL*8  :: glmelt_fac_nh = 1.d0
!@dbparam glmelt_fac_sh is a factor to multiply glacial melt by in SH
      REAL*8  :: glmelt_fac_sh = 1.d0

!@var FWAREA_SH, FWAREA_NH area in each hemisphere for glmelt
      REAL*8  :: FWAREA_SH, FWAREA_NH

C**** The net fluxes of glacier mass into the ocean (per hemisphere)
!@var ACCPDA total accumulation per year for Antarctica (kg/yr)
!@var ACCPDG total accumulation per year for Greenland (kg/yr)
!@var EACCPDA total energy flux per year for Antarctica (kg/yr)
!@var EACCPDG total energy flux per year for Greenland (kg/yr)
      REAL*8 :: ACCPDA, ACCPDG,  EACCPDA, EACCPDG

C**** NOTE: ICB arrays are purely diagnostic and kept for
C**** accounting purposes only. They are not real icebergs!
!@var MICBIMP( ) + Sum[MDWNIMP(over hemisphere)] = instantaneous
!@+     mass of icebergs (implicit) in each hemisphere (kg)
!@var EICBIMP( ) + Sum[EDWNIMP(over hemisphere)] = instantaneous
!@+     energy of icebergs (implicit) in each hemisphere (J)
      Real*8,Dimension(2) :: MICBIMP,EICBIMP

#ifdef TRACERS_WATER /* TNL: inserted */
#ifdef TRACERS_OCEAN
!@var TRACCPDA total tracer flux per year for Antarctica (kg/yr)
!@var TRACCPDG total tracer flux per year for Greenland (kg/yr)
      !REAL*8 :: TRACCPDA(NTM), TRACCPDG(NTM)
      ! landice_com should really be the owner of *ACC[PDA,PDG]
      ! and other "iceberg" variables, since these variables are
      ! not intrinsic to the 1D ice physics.
      ! Allocation using runtime-determined NTM is in alloc_landice_com
      REAL*8, ALLOCATABLE :: TRACCPDA(:), TRACCPDG(:)
!@var TRICBIMP( ) + Sum[TRDWNIMP(over hemisphere)] = instantaneous
!@var   tracer content of icebergs (implicit) in each hemisphere (kg)
c      REAL*8, DIMENSION(NTM,2) :: TRICBIMP
#endif
#endif   /* TNL: inserted */

      CONTAINS

      SUBROUTINE PRECLI(SNOW,TG1,TG2,PRCP,ENRGP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TRPRCP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,ERUN2,RUN0)
!@sum  PRECLI apply precipitation to land ice fraction
!@auth Original Development Team
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: SNOW,TG1,TG2
      REAL*8, INTENT(IN) :: PRCP,ENRGP
      REAL*8, INTENT(OUT) :: EDIFS,DIFS,ERUN2,RUN0
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRLI
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRDIFS
#endif

      REAL*8 HC1,DWATER
C**** initiallize output
      EDIFS=0. ; DIFS=0. ; ERUN2=0. ; RUN0=0.
#ifdef TRACERS_WATER
      TRDIFS(:)=0. ; TRUN0(:)=0.
#endif

      HC1=HC1LI+SNOW*SHI
      IF (ENRGP.GE.0.) THEN
        IF (ENRGP.GE.-TG1*HC1) THEN
C**** RAIN HEATS UP TG1 TO FREEZING POINT AND MELTS SOME SNOW OR ICE
          DWATER=(TG1*HC1+ENRGP)/LHM
          TG1=0.
          RUN0=DWATER+PRCP
          IF (DWATER.LT.SNOW) THEN
C**** RAIN MELTS SOME SNOW
#ifdef TRACERS_WATER
            TRUN0(:)=TRPRCP(:)+DWATER*TRSNOW(:)/SNOW
            TRSNOW(:)=TRSNOW(:)*(1.-DWATER/SNOW)
#endif
            SNOW=SNOW-DWATER
          ELSE
C**** RAIN MELTS ALL SNOW AND SOME ICE, ICE MOVES UP THROUGH THE LAYERS
            DIFS=SNOW-DWATER  ! < 0
            SNOW=0.
            TG1=-TG2*DIFS/ACE1LI
            EDIFS=DIFS*(TG2*SHI-LHM)
            ERUN2=EDIFS
#ifdef TRACERS_WATER
            TRUN0(:)=TRPRCP(:)+TRSNOW(:)-DIFS*TRLI(:)/(ACE2LI+ACE1LI)
            TRDIFS(:)=DIFS*TRLI(:)/(ACE2LI+ACE1LI)
            TRSNOW(:)=0.
#endif
          END IF
        ELSE
C**** RAIN COOLS TO FREEZING POINT AND HEATS UP TG1
          TG1=TG1+ENRGP/HC1
          RUN0=PRCP
#ifdef TRACERS_WATER
          TRUN0(:)=TRPRCP(:)
#endif
        END IF
      ELSE
C**** SNOW INCREASES SNOW AMOUNT AND SNOW TEMPERATURE RECOMPUTES TG1
        TG1=(TG1*HC1+ENRGP+LHM*PRCP)/(HC1+PRCP*SHI)
        SNOW=SNOW+PRCP
#ifdef TRACERS_WATER
        TRSNOW(:)=TRSNOW(:)+TRPRCP(:)
#endif
        IF (SNOW.gt.ACE1LI) THEN
C**** SNOW IS COMPACTED INTO ICE, ICE MOVES DOWN THROUGH THE LAYERS
          DIFS=SNOW-.9d0*ACE1LI
#ifdef TRACERS_WATER
          TRDIFS(:)=TRLI(:)*DIFS/(ACE2LI+ACE1LI)
          TRLI(:)=TRLI(:)+DIFS*(TRSNOW(:)/SNOW-TRLI(:)/(ACE2LI+ACE1LI))
          TRSNOW(:)=TRSNOW(:)*(1.-DIFS/SNOW)
#endif
          SNOW=.9d0*ACE1LI
          EDIFS=DIFS*(TG1*SHI-LHM)
          ERUN2=DIFS*(TG2*SHI-LHM)
          TG2=TG2+(TG1-TG2)*DIFS/ACE2LI
        END IF
      END IF
      RETURN
      END SUBROUTINE PRECLI

      SUBROUTINE LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TREVAP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,RUN0)
!@sum  LNDICE apply surface fluxes to land ice fraction
!@auth Original Development team
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: SNOW,TG1,TG2
      REAL*8, INTENT(IN) :: F0DT,F1DT,EVAP
      REAL*8, INTENT(OUT) :: EDIFS,DIFS,RUN0
#ifdef TRACERS_WATER
C**** Tracer concentration in ice assumed constant
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRLI
!@var TREVAP tracer amount in evaporation (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TREVAP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRDIFS
      INTEGER N
#endif
      REAL*8 :: SNANDI,HC1,ENRG1

C**** initiallize output
      EDIFS=0. ; DIFS=0. ; RUN0=0.
#ifdef TRACERS_WATER
      TRDIFS(:)=0. ; TRUN0(:)=0.
#endif

C**** CALCULATE TG1
      SNANDI=SNOW+ACE1LI-EVAP
      HC1=SNANDI*SHI
      ENRG1=F0DT+EVAP*(TG1*SHI-LHM)-F1DT
      IF (ENRG1.GT.-TG1*HC1) THEN
C**** FLUXES HEAT UP TG1 TO FREEZING POINT AND MELT SOME SNOW AND ICE
        RUN0=(ENRG1+TG1*HC1)/LHM
        TG1=0.
        SNANDI=SNANDI-RUN0
#ifdef TRACERS_WATER
        IF (SNOW-EVAP.gt.RUN0) THEN  ! only snow melts
          TRUN0(:)=RUN0*(TRSNOW(:)-TREVAP(:))/(SNOW-EVAP)
          TRSNOW(:)=TRSNOW(:)-TREVAP(:)-TRUN0(:)
        ELSE ! all snow + some ice melts
          TRUN0(:)=TRSNOW(:)-TREVAP(:)+(ACE1LI-SNANDI)*TRLI(:)/(ACE2LI
     *         +ACE1LI)
c              ! done below
c          TRSNOW(:)=0.
c          TRLI(:)=TRLI(:)*(1.- ((ACE1LI-SNANDI)/(ACE2LI+ACE1LI)))
          do n=1,ntm
            TRUN0(N)=MAX(TRUN0(N),0d0)
          end do
        END IF
#endif
      ELSE
C**** FLUXES RECOMPUTE TG1 WHICH IS BELOW FREEZING POINT
        TG1=TG1+ENRG1/HC1
#ifdef TRACERS_WATER
        if (SNOW.gt.EVAP) THEN
          TRSNOW(:)=TRSNOW(:)-TREVAP(:)
c        else   ! done below
c          TRLI(:)=TRLI(:)-(TREVAP(:)-TRSNOW(:))
c          TRSNOW(:)=0.
        end if
#endif
      END IF
      IF (SNANDI.LT.ACE1LI) THEN
C**** SOME ICE HAS MELTED OR EVAPORATED, TAKE IT FROM G2
        SNOW=0.
        DIFS=SNANDI-ACE1LI
        TG1=(TG1*SNANDI-TG2*DIFS)/ACE1LI
        EDIFS=DIFS*(TG2*SHI-LHM)
#ifdef TRACERS_WATER
        TRDIFS(:)=TRLI(:)*(DIFS/(ACE2LI+ACE1LI))
        TRLI(:)=TRLI(:)-DIFS*TRLI(:)/(ACE2LI+ACE1LI)+
     *       TRSNOW(:)-TREVAP(:)-TRUN0(:)
        TRSNOW(:)=0.
#endif
      ELSE
        SNOW=SNANDI-ACE1LI
      END IF
C**** CALCULATE TG2
      TG2=TG2+F1DT/HC2LI

      RETURN
      END SUBROUTINE LNDICE

      END MODULE LANDICE

