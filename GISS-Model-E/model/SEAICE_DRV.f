#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
!@sum  SEAICE_DRV contains drivers for SEAICE related routines
!@auth Gavin Schmidt
!@cont PRECIP_SI,GROUND_SI

      SUBROUTINE CALC_APRESS(atmice)
C**** Calculate pressure anomaly at ocean surface
      USE CONSTANT, only : grav
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : iceocn,si_ocn
      USE EXCHANGE_TYPES, only : atmice_xchng_vars
      IMPLICIT NONE
      type(atmice_xchng_vars) :: atmice
c
      INTEGER I,J
C****
C**** Extract useful local domain parameters from "grid"
C****
      integer :: J_0, J_1, I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      I_0 = atmice%I_0
      I_1 = atmice%I_1
      J_0 = atmice%J_0
      J_1 = atmice%J_1

      DO J=J_0, J_1
      DO I=I_0,atmice%IMAXJ(J)
        iceocn%APRESS(I,J) = 100.*(atmice%SRFP(I,J)-1013.25d0)+
     *       si_ocn%RSI(I,J)*
     *       (si_ocn%SNOWI(I,J)+ACE1I+si_ocn%MSI(I,J))*GRAV
      END DO
      END DO
      IF (atmice%HAVE_SOUTH_POLE)
     &     iceocn%APRESS(2:I_1,1)   = iceocn%APRESS(1,1)
      IF (atmice%HAVE_NORTH_POLE)
     &     iceocn%APRESS(2:I_1,J_1) = iceocn%APRESS(1,J_1)
      RETURN
      END SUBROUTINE CALC_APRESS

      SUBROUTINE PRECIP_SI(si_state,iceocn,atmice)
!@sum  PRECIP_SI driver for applying precipitation to sea ice fraction
!@auth Original Development team
!@calls seaice:prec_si
      USE CONSTANT, only : teeny,grav,tf,bylhm
      USE EXCHANGE_TYPES, only : atmice_xchng_vars,iceocn_xchng_vars
      USE SEAICE_COM, only : icestate
#ifdef TRACERS_WATER
      USE SEAICE, only : ntm
#endif
      USE SEAICE, only : prec_si, ace1i, lmi,xsi,debug
      IMPLICIT NONE
      type(icestate) :: si_state
      type(iceocn_xchng_vars) :: iceocn
      type(atmice_xchng_vars) :: atmice
      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 SNOW,MSI2,PRCP,ENRGP,RUN0,POICE,SRUN0,ERUN0
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, DIMENSION(NTM) :: TRUN0,TRPRCP
#endif
      CHARACTER(LEN=8) :: DOMAIN
      INTEGER I,J
      LOGICAL WETSNOW,DOPOINT
      real*8, dimension(:,:), pointer :: runpsi,srunpsi,erunpsi
     &     ,prec,eprec,fwater
      real*8, dimension(:,:), pointer :: rsi,msi,snowi,pond_melt
      logical, dimension(:,:), pointer :: flag_dsws
      real*8, dimension(:,:,:), pointer :: hsi,ssi
#ifdef TRACERS_WATER
      real*8, dimension(:,:,:,:), pointer :: trsi
      real*8, dimension(:,:,:), pointer :: trprec,trunpsi
#endif
      integer :: J_0, J_1, I_0,I_1
      REAL*8 CMPRS
      real*8, dimension(:,:,:), pointer :: aij

c#ifdef STANDALONE_OCEAN
c      call stop_model('add snow aging in here',255)
c#endif

      I_0 = si_state%I_0
      I_1 = si_state%I_1
      J_0 = si_state%J_0
      J_1 = si_state%J_1

      domain = si_state%domain

c
c set pointers
c
      rsi => si_state%rsi
      msi => si_state%msi
      hsi => si_state%hsi
      ssi => si_state%ssi
      snowi => si_state%snowi
      pond_melt => si_state%pond_melt
      flag_dsws => si_state%flag_dsws
#ifdef TRACERS_WATER
      trsi => si_state%trsi
#endif

       fwater => iceocn%fwater
       runpsi => iceocn%runpsi
      srunpsi => iceocn%srunpsi
      erunpsi => iceocn%erunpsi
       prec => atmice%prec
      eprec => atmice%eprec
#ifdef TRACERS_WATER
      trprec => atmice%trprec
      trunpsi => iceocn%trunpsi
#endif

      aij => atmice%aij

      DO J=J_0, J_1
      DO I=I_0,si_state%IMAXJ(J)
      POICE=    RSI(I,J) *FWATER(I,J)
      IF(POICE.EQ.0.) THEN
        RUNPSI(I,J)=0
        SRUNPSI(I,J)=0
        ERUNPSI(I,J)=0
#ifdef TRACERS_WATER
        TRUNPSI(:,I,J)=0
#endif
      ENDIF
      DOPOINT = POICE.GT.0.
      IF (DOPOINT) THEN ! todo: no more dopoint
        PRCP=PREC(I,J)
        ENRGP=EPREC(I,J)      ! energy of precip
        SNOW=SNOWI(I,J)
        MSI2=MSI(I,J)
        HSIL(:) = HSI(:,I,J)      ! sea ice temperatures
        SSIL(:) = SSI(:,I,J)      ! sea ice salt
#ifdef TRACERS_WATER
        TRSIL(:,:)=TRSI(:,:,I,J)  ! sea ice tracers
        TRPRCP(:)=TRPREC(:,I,J)   ! tracer in precip
#endif
        CMPRS=0d0

C**** CALL SUBROUTINE FOR CALCULATION OF PRECIPITATION OVER SEA ICE

        CALL PREC_SI(SNOW,MSI2,HSIL,TSIL,SSIL,PRCP,ENRGP,RUN0,SRUN0
     *       ,ERUN0,
#ifdef TRACERS_WATER
     *       TRSIL,TRPRCP,TRUN0,
#endif
     *       WETSNOW,CMPRS)

        SNOWI(I,J)  =SNOW
        RUNPSI(I,J) =RUN0
        SRUNPSI(I,J)=SRUN0
        ERUNPSI(I,J)=ERUN0
        HSI(:,I,J)=HSIL(:)
        SSI(:,I,J)=SSIL(:)
#ifdef TRACERS_WATER
        TRSI(:,:,I,J)=TRSIL(:,:)
        TRUNPSI(:,I,J)=TRUN0(:) ! tracer in runoff
#endif
        FLAG_DSWS(I,J)=FLAG_DSWS(I,J).or.WETSNOW
C**** reset flag if there was fresh snow (i.e. prcp but no rain!)
        IF (.not. WETSNOW .and. PRCP.gt.0.) FLAG_DSWS(I,J)=.FALSE.
C**** pond_melt accumulation
        pond_melt(i,j)=pond_melt(i,j)+0.3d0*RUN0

        MSI(I,J)=MSI2

#ifndef STANDALONE_OCEAN
C**** snow-to-ice accumulation
        AIJ(I,J,ATMICE%IJ_SNTOSI)=AIJ(I,J,ATMICE%IJ_SNTOSI)
     &         +POICE*CMPRS
C****
#endif

      END IF

      END DO
      END DO

C****
      END SUBROUTINE PRECIP_SI

      SUBROUTINE UNDERICE(si_state,iceocn,atmocn)
!@sum  underice calculates basal fluxes under sea and lake ice
!@+    saves the resulting fluxes
!@auth Gavin Schmidt
!@calls iceocean_fluxes,icelake_fluxes
      USE CONSTANT, only : rhow,rhows,omega,rhoi,shw
      USE MODEL_COM, only : dtsrc,qcheck,kocean
#ifdef TRACERS_WATER
      USE SEAICE, only : ntm
#endif
      USE SEAICE, only : lmi,xsi,icelake_fluxes,iceocean_fluxes,
     *     ac2oim,alpha,tfrez,debug,Ti,dEidTi,alami
      USE EXCHANGE_TYPES, only : iceocn_xchng_vars,atmocn_xchng_vars
      USE SEAICE_COM, only : icestate
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      IMPLICIT NONE
      type(icestate) :: si_state
      type(iceocn_xchng_vars) :: iceocn
      type(atmocn_xchng_vars) :: atmocn
      CHARACTER(LEN=8) :: DOMAIN
      INTEGER I,J
      LOGICAL :: DOPOINT
      REAL*8 coriol,ustar,Tm,Sm,Si,Tic,dh,mflux,hflux,sflux,fluxlim
     *     ,mlsh,icefrac   !,mfluxmax
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM) :: Trm,Tri,trflux,tralpha
#ifdef TRACERS_SPECIAL_O18
      REAL*8 fracls
      INTEGER N
#endif
#endif
      real*8, dimension(:,:), pointer :: fmsi_io,fhsi_io,fssi_io,ui2rho
     &     ,sss,mlhc,fwater
     &     ,mldlk,dlake,glake
      real*8, dimension(:,:), pointer :: rsi,msi
      real*8, dimension(:,:,:), pointer :: hsi,ssi
#ifdef TRACERS_WATER
      real*8, dimension(:,:,:,:), pointer :: trsi
      real*8, dimension(:,:,:), pointer :: mltracer,ftrsi_io
#endif
      integer :: J_0, J_1 ,I_0,I_1

      call startTimer('UNDERICE()')

      domain = si_state%domain

      I_0 = si_state%I_0
      I_1 = si_state%I_1
      J_0 = si_state%J_0
      J_1 = si_state%J_1

c
c set pointers
c
      rsi => si_state%rsi
      msi => si_state%msi
      hsi => si_state%hsi
      ssi => si_state%ssi
#ifdef TRACERS_WATER
      trsi => si_state%trsi
#endif
       fwater => iceocn%fwater
      fmsi_io => iceocn%fmsi_io
      fhsi_io => iceocn%fhsi_io
      fssi_io => iceocn%fssi_io
      ui2rho  => iceocn%ui2rho
      sss     => atmocn%sss
      mlhc    => atmocn%mlhc
      if (domain.eq.'LAKES') then
        mldlk => iceocn%mldlk
        dlake => iceocn%dlake
        glake => iceocn%glake
      endif
#ifdef TRACERS_WATER
      mltracer => atmocn%gtracer
      ftrsi_io => iceocn%ftrsi_io
#endif

      DO J=J_0, J_1
C**** Coriolis parameter (on tracer grid)
        DO I=I_0,si_state%IMAXJ(J)
          coriol = iceocn%coriol(i,j) !ABS(2.*OMEGA*SINLAT2D(I,J))
          icefrac = RSI(I,J)*FWATER(I,J)
          if(icefrac.eq.0.) then
            FMSI_IO(I,J) = 0.
            FHSI_IO(I,J) = 0.
            FSSI_IO(I,J) = 0.
#ifdef TRACERS_WATER
            FTRSI_IO(:,I,J)=0.
#endif
          endif
          DOPOINT = ICEFRAC.GT.0.
          IF (DOPOINT) THEN  ! todo: no more dopoint
C**** Set mixed layer conditions
            Tm = atmocn%gtemp(i,j)
#ifdef TRACERS_WATER
            Trm(:)=MLTRACER(:,I,J)
            Tri(:)=TRSI(:,LMI,I,J)/(XSI(LMI)*MSI(I,J)-SSI(LMI,I,J))
#ifdef TRACERS_SPECIAL_O18
            do n=1,ntm
              tralpha(n)=fracls(n)
            end do
#else
            tralpha(:)=1.
#endif
#endif
            dh = 0.5*(XSI(LMI)*MSI(I,J))/RHOI
c            mfluxmax = (MSI(I,J)-AC2OIM)/dtsrc
            IF (DOMAIN.EQ.'OCEAN') THEN
C**** Ice lowest layer conditions
              Si = 1d3*SSI(LMI,I,J)/(XSI(LMI)*MSI(I,J))
              Tic = Ti(HSI(LMI,I,J)/(XSI(LMI)*MSI(I,J)),Si)
              IF (KOCEAN.ge.1) THEN
C**** should we calculate ocean rho(Tm,Sm) here?
                Ustar = MAX(5d-4,SQRT(UI2rho(I,J)/RHOWS))
                Sm = SSS(I,J)
                mlsh = MLHC(I,J)
                call iceocean_fluxes(
     *               Tic,Si,Tm,Sm,dh,Ustar,Coriol,dtsrc,mlsh,
#ifdef TRACERS_WATER
     *               Tri,Trm,trflux,tralpha,
#endif
     *               mflux,sflux,hflux)
              ELSE ! for fixed SST assume freezing temp at base,implicit
                hflux=alami(Tic,Si)*(Tic-tfrez(sss(i,j)))/(dh+alpha
     *               *dtsrc*alami(Tic,Si)/(dEidTi(Tic,Si)*XSI(LMI)*MSI(I
     *               ,J)))
                mflux=0.
                sflux=0.
#ifdef TRACERS_WATER
                trflux= 0.
#endif
              END IF
            ELSE   ! for lakes (no salinity so solve directly)
              Tic =Ti(HSI(LMI,I,J)/(XSI(LMI)*MSI(I,J)),0d0)
              mlsh=SHW*MLDLK(I,J)*RHOW
              sflux = 0.
              call icelake_fluxes(Tic,Tm,dh,dtsrc,mlsh,
#ifdef TRACERS_WATER
     *             Tri,Trm,trflux,tralpha,
#endif
     *             mflux,hflux)

C**** Limit lake-to-ice flux if lake is too shallow (< 40cm)
              IF (DLAKE(I,J).lt.0.4d0) THEN
                FLUXLIM=-GLAKE(I,J)/DTSRC
                IF (hflux.lt.FLUXLIM) hflux = FLUXLIM
                if (mflux.lt.0) then
                  mflux = 0.
#ifdef TRACERS_WATER
                  trflux= 0.
#endif
                end if
                if (qcheck) print*,"Flux limiting",I,J,DLAKE(I,J),
     *               FLUXLIM*DTSRC
              END IF
            END IF
            FMSI_IO(I,J) = mflux*dtsrc   ! positive down
            FHSI_IO(I,J) = hflux*dtsrc
            FSSI_IO(I,J) = sflux*dtsrc
#ifdef TRACERS_WATER
            FTRSI_IO(:,I,J)=trflux(:)*dtsrc
#endif
          END IF
        END DO
      END DO

C****
      call stopTimer('UNDERICE()')
      RETURN
      END SUBROUTINE UNDERICE

      SUBROUTINE MELT_SI(si_state,iceocn,atmocn,atmice)
!@sum  MELT_SI driver for lateral melt of sea ice
!@auth Gary Russell/Gavin Schmidt
!@calls SEAICE:SIMELT
      USE CONSTANT, only : TF
      USE MODEL_COM, only : kocean,dtsrc
      USE SEAICE, only : lmi,simelt,tfrez,xsi,Ti,ace1i,debug
      USE EXCHANGE_TYPES, only : iceocn_xchng_vars,
     &     atmocn_xchng_vars,atmice_xchng_vars
      USE SEAICE_COM, only : icestate
#ifdef TRACERS_WATER
      USE SEAICE, only : ntm
#endif
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      IMPLICIT NONE
      type(icestate) :: si_state
      type(iceocn_xchng_vars) :: iceocn
      type(atmocn_xchng_vars) :: atmocn
      type(atmice_xchng_vars) :: atmice
      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 MSI2,ROICE,SNOW,ENRGUSED,RUN0,SALT,POCEAN,TFO
     *     ,PWATER,Tm,DT,ENRGMAX
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, DIMENSION(NTM) :: TRUN0
#endif
      CHARACTER(LEN=8) :: DOMAIN
      LOGICAL :: DOPOINT
      INTEGER I,J
      real*8, dimension(:,:), pointer :: melti,emelti,smelti,sss,mlhc
     &     ,fwater,rsistart
      real*8, dimension(:,:), pointer :: rsi,msi,snowi
      real*8, dimension(:,:,:), pointer :: hsi,ssi
#ifdef TRACERS_WATER
      real*8, dimension(:,:,:,:), pointer :: trsi
      real*8, dimension(:,:,:), pointer :: trmelti
#endif
      integer :: J_0, J_1 ,I_0,I_1

      call startTimer('MELT_SI()')

      I_0 = si_state%I_0
      I_1 = si_state%I_1
      J_0 = si_state%J_0
      J_1 = si_state%J_1

      domain = si_state%domain
c
c set pointers
c
      rsi => si_state%rsi
      msi => si_state%msi
      hsi => si_state%hsi
      ssi => si_state%ssi
      snowi => si_state%snowi
#ifdef TRACERS_WATER
      trsi => si_state%trsi
#endif
      fwater => iceocn%fwater
       melti => iceocn%melti
      emelti => iceocn%emelti
      smelti => iceocn%smelti
      sss    => atmocn%sss
      mlhc   => atmocn%mlhc
      rsistart => atmice%rsistart
#ifdef TRACERS_WATER
      trmelti => iceocn%trmelti
#endif

      DO J=J_0, J_1
      DO I=I_0,si_state%IMAXJ(J)
        IF(FWATER(I,J).GT.0.) THEN
          RSIstart(I,J)=RSI(I,J)
          MELTI(I,J) = 0.
          EMELTI(I,J)= 0.
          SMELTI(I,J)= 0.
#ifdef TRACERS_WATER
          TRMELTI(:,I,J)=0.
#endif
        ENDIF
      ENDDO
      ENDDO

C**** CALCULATE LATERAL MELT (ALSO ELIMINATE SMALL AMOUNTS)
C**** EVERY PHYSICS TIME STEP
      DT=DTsrc
      DO J=J_0, J_1
        DO I=I_0,si_state%IMAXJ(J)
          PWATER=FWATER(I,J)
          ROICE=RSI(I,J)
C**** Call simelt if (lake and v. small ice) or (q-flux ocean, some ice)
C**** now include lat melt for lakes and any RSI < 1
          DOPOINT = PWATER*ROICE.GT.0.
          IF(DOMAIN.EQ.'OCEAN') THEN
            DOPOINT = DOPOINT .AND. KOCEAN.GE.1
            POCEAN = PWATER
          ELSE
            DOPOINT = DOPOINT .AND. ROICE.LT.1.
            POCEAN = 0.
          ENDIF
          IF(DOPOINT) THEN
            IF (POCEAN.gt.0) THEN
              TFO = tfrez(sss(i,j))
            ELSE
              TFO = 0.
            END IF
            Tm=atmocn%gtemp(i,j)
            MSI2=MSI(I,J)
            SNOW=SNOWI(I,J)     ! snow mass
            HSIL(:)= HSI(:,I,J) ! sea ice enthalpy
            SSIL(:)= SSI(:,I,J) ! sea ice salt
            ENRGMAX= MAX(Tm-TFO,0d0)*MLHC(I,J) ! max energy for melt
#ifdef TRACERS_WATER
            TRSIL(:,:)=TRSI(:,:,I,J) ! tracer content of sea ice
#endif
            CALL SIMELT(DT,ROICE,SNOW,MSI2,HSIL,SSIL,POCEAN,Tm,TFO,TSIL
#ifdef TRACERS_WATER
     *           ,TRSIL,TRUN0
#endif
     *           ,ENRGMAX,ENRGUSED,RUN0,SALT)

#ifndef STANDALONE_OCEAN
C**** Update prognostic sea ice variables + correction for rad. fluxes
          IF(DOMAIN.EQ.'LAKES') THEN
            if (roice.gt.rsi(i,j)) ! ice from ocean
     *          call RESET_SURF_FLUXES(I,J,1,2,RSI(I,J),ROICE)
            if (roice.lt.rsi(i,j)) ! ocean from ice
     *          call RESET_SURF_FLUXES(I,J,2,1,1.-RSI(I,J),1.-ROICE)
          ENDIF
C****
#endif
            RSI(I,J)=ROICE
            RSIstart(I,J)=ROICE
            MSI(I,J)=MSI2
            SNOWI(I,J)=SNOW
            HSI(:,I,J)=HSIL(:)
            SSI(:,I,J)=SSIL(:)
#ifdef TRACERS_WATER
            TRSI(:,:,I,J)=TRSIL(:,:)
#endif
C**** Save fluxes (in kg, J etc.), positive into ocean
          MELTI(I,J) = RUN0*PWATER
          EMELTI(I,J)=-ENRGUSED*PWATER
          SMELTI(I,J)= SALT*PWATER
#ifdef TRACERS_WATER
          TRMELTI(:,I,J)=TRUN0(:)*PWATER
#endif
          END IF

C****
        END DO
      END DO
C****

C**** replicate ice values at the poles
      IF(DOMAIN.EQ.'OCEAN') THEN
      IF (si_state%HAVE_NORTH_POLE) THEN
        DO I=2,I_1
          RSI(I,J_1)=RSI(1,J_1)
          MSI(I,J_1)=MSI(1,J_1)
          HSI(:,I,J_1)=HSI(:,1,J_1)
          SSI(:,I,J_1)=SSI(:,1,J_1)
          SNOWI(I,J_1)=SNOWI(1,J_1)
#ifdef TRACERS_WATER
          TRSI(:,:,I,J_1) = TRSI(:,:,1,J_1)
#endif
        END DO
      END IF
      IF (si_state%HAVE_SOUTH_POLE) THEN
        DO I=2,I_1
          RSI(I,1)=RSI(1,1)
          MSI(I,1)=MSI(1,1)
          HSI(:,I,1)=HSI(:,1,1)
          SSI(:,I,1)=SSI(:,1,1)
          SNOWI(I,1)=SNOWI(1,1)
#ifdef TRACERS_WATER
          TRSI(:,:,I,1) = TRSI(:,:,1,1)
#endif
        END DO
      END IF
      END IF

      call stopTimer('MELT_SI()')
      RETURN
      END SUBROUTINE MELT_SI

      SUBROUTINE GROUND_SI(si_state,iceocn,atmice,atmocn)
!@sum  GROUND_SI driver for applying surface + base fluxes to sea ice
!@auth Gary Russell/Gavin Schmidt
!@ver  2010/11/12
!@calls SEAICE:SEA_ICE
      USE CONSTANT, only : grav,rhows,rhow
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE MODEL_COM, only : dtsrc
      USE EXCHANGE_TYPES, only :
     &     atmice_xchng_vars,iceocn_xchng_vars,atmocn_xchng_vars
      USE SEAICE_COM, only : icestate
#ifdef TRACERS_WATER
      USE SEAICE, only: ntm
#endif
      USE SEAICE, only : sea_ice,ssidec,lmi,xsi,ace1i,qsfix,debug
     *     ,snowice, snow_ice, rhos, Ti, Ti2b
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      IMPLICIT NONE
      type(icestate) :: si_state
      type(iceocn_xchng_vars) :: iceocn
      type(atmice_xchng_vars) :: atmice
      type(atmocn_xchng_vars) :: atmocn
      REAL*8, DIMENSION(LMI) :: HSIL,SSIL
      REAL*8 SNOW,ROICE,MSI2,F0DT,F1DT,EVAP,SROX(2)
     *     ,FMOC,FHOC,FSOC,POICE,PWATER,SCOVI
      REAL*8 MFLUX,HFLUX,SFLUX,RUN,ERUN,SRUN,MELT12
      REAL*8 MSNWIC,HSNWIC,SSNWIC,SM,TM,Ti1,DSNOW
      CHARACTER(LEN=8) :: DOMAIN
      INTEGER I,J
      LOGICAL WETSNOW,DOPOINT
      real*8, dimension(:,:), pointer :: runosi,erunosi,srunosi,
     *     fmsi_io,fhsi_io,fssi_io,solar_io,solar,e0,e1,evapor,sss
     *     ,fwater
     &     ,MSIsave,SNTOSI,SITOPMLT,MSNFLOOD,HSNFLOOD,TI1save,SIHC
     &     ,SNOWsave
      real*8, dimension(:,:), pointer :: rsi,msi,snowi,pond_melt
      logical, dimension(:,:), pointer :: flag_dsws
      real*8, dimension(:,:,:), pointer :: hsi,ssi
#ifdef TRACERS_WATER
      real*8, dimension(:,:,:,:), pointer :: trsi
      real*8, dimension(:,:,:), pointer :: trevapor,mltracer
     &     ,ftrsi_io,trunosi
#ifdef TRACERS_DRYDEP
      real*8, dimension(:,:,:), pointer :: trdrydep
#endif
#endif
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: trsil
      REAL*8, DIMENSION(NTM) :: trflux,ftroc,trevap,trrun,trsnwic,trm
     *     ,tralpha
#endif
      integer :: J_0, J_1, I_0,I_1
      real*8 MSI1,SNOWL1,MICE1
      REAL*8 CMPRS

      call startTimer('GROUND_SI()')

      I_0 = si_state%I_0
      I_1 = si_state%I_1
      J_0 = si_state%J_0
      J_1 = si_state%J_1

      domain = si_state%domain
c
c set pointers
c
      rsi => si_state%rsi
      msi => si_state%msi
      hsi => si_state%hsi
      ssi => si_state%ssi
      snowi => si_state%snowi
      pond_melt => si_state%pond_melt
      flag_dsws => si_state%flag_dsws
#ifdef TRACERS_WATER
      trsi => si_state%trsi
#endif
       fwater => iceocn%fwater
       runosi => iceocn%runosi
      srunosi => iceocn%srunosi
      erunosi => iceocn%erunosi
      fmsi_io => iceocn%fmsi_io
      fhsi_io => iceocn%fhsi_io
      fssi_io => iceocn%fssi_io
      solar_io=> iceocn%solar
      sss     => atmocn%sss
#ifdef TRACERS_WATER
      mltracer => atmocn%gtracer
      trunosi => iceocn%trunosi
      ftrsi_io => iceocn%ftrsi_io
#endif
      e0     => atmice%e0
      e1     => atmice%e1
      evapor => atmice%evapor
      solar  => atmice%solar
#ifdef TRACERS_WATER
      trevapor => atmice%trevapor
#ifdef TRACERS_DRYDEP
      trdrydep => atmice%trdrydep
#endif
#endif

      MSIsave  => atmice%MSIsave
      SNTOSI   => atmice%SNTOSI
      SITOPMLT => atmice%SITOPMLT
      MSNFLOOD => atmice%MSNFLOOD
      HSNFLOOD => atmice%HSNFLOOD
      TI1save  => atmice%TI1save
      SIHC     => atmice%SIHC
      SNOWsave => atmice%SNOWsave

      debug=.false.

      DO J=J_0, J_1
      DO I=I_0,si_state%IMAXJ(J)
c        debug = i.eq.52.and.j.eq.87
c      if (debug) print*,"gs0",trsi(1,1,i,j)-(snowi(i,j)+ace1i)*xsi(1)+
c     &             ssi(1,i,j), trsi(1,2,i,j)-(snowi(i,j)+ace1i)
c     *     *xsi(2)+ssi(2,i,j),trsi(1,3,i,j)-msi(i,j)*xsi(3)+
c     *     ssi(3,i,j),trsi(1,4,i,j)-msi(i,j)*xsi(4)+ssi(4,i,j)

      PWATER=FWATER(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      IF (POICE.eq.0) THEN
        SOLAR_IO(I,J)=0
        RUNOSI(I,J)=0
        ERUNOSI(I,J)=0
        SRUNOSI(I,J)=0
#ifdef TRACERS_WATER
        TRUNOSI(:,I,J) = 0.
#endif
      ENDIF
      DOPOINT = POICE.GT.0.
      IF(DOPOINT) THEN ! todo: no more dopoint
c      IF (POICE.gt.0) THEN
        F0DT=E0(I,J)   ! heat flux to the top ice surface (J/m^2)
        F1DT=E1(I,J)   ! heat flux between 1st and 2nd ice layer (J/m^2)
        EVAP=EVAPOR(I,J)   ! evaporation/dew at the ice surface (kg/m^2)
        SROX(1)=SOLAR(I,J) ! solar radiation absrbd by sea ice (J/m^2)
        FMOC=fmsi_io(i,j)  ! mass flux at base (kg/m^2)
        FSOC=fssi_io(i,j)  ! salt flux at base (kg/m^2)
        FHOC=fhsi_io(i,j)  ! heat flux at base (J/M^2)
        SNOW= SNOWI(I,J)  ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
        HSIL(:) = HSI(:,I,J)  ! sea ice enthalpy
        SSIL(:) = SSI(:,I,J)  ! sea ice salt
        WETSNOW=FLAG_DSWS(I,J)  ! wetness of snow
        Tm=atmocn%gtemp(i,j) ! ocean mixed layer temperature (C)
#ifdef TRACERS_WATER
        TREVAP(:) = TREVAPOR(:,I,J)
#ifdef TRACERS_DRYDEP
     *       -trdrydep(:,i,j)
#endif
        FTROC(:)  = ftrsi_io(:,i,j)
        TRSIL(:,:)= TRSI(:,:,I,J)
        Trm(:)=MLTRACER(:,I,J)
        Tralpha(:)=1d0  ! no fractionation for snow ice formation
#endif
        IF(DOMAIN.EQ.'OCEAN') THEN
          Sm=SSS(I,J)           ! ocean mixed layer salinity (psu)
        ELSE
          Sm=0.                 ! lakes always fresh
        END IF

        MSIsave(I,J) = (ACE1I+MSI2)*POICE

        CMPRS=0d0

c        if (debug) write(6,'(A,2I4,4F11.6)') "si0",i,j,SNOW,1d3*SSIL(1)
c     $       /(XSI(1)*(ACE1I+SNOW)),1d3*SSIL(2)/(XSI(2)*(ACE1I+SNOW))
c     $       ,1d3*(SSIL(1)+SSIL(2))/ACE1I
c        if (debug) print*,"si0",i,j,SNOW,ROICE,HSIL,SSIL,MSI2,F0DT,F1DT,
c     *       EVAP,SROX,SNOW+ACE1I-SSIL(1)-SSIL(2) !, TRSIL(1,1)+TRSIL(1,2)

        CALL SEA_ICE(DTSRC,SNOW,ROICE,HSIL,SSIL,MSI2,F0DT,F1DT,EVAP,SROX
#ifdef TRACERS_WATER
     *       ,TRSIL,TREVAP,FTROC,TRRUN
#endif
     *       ,FMOC,FHOC,FSOC,RUN,ERUN,SRUN,WETSNOW,MELT12,CMPRS)
c        if (debug)  write(6,'(A,2I4,4F11.6)') "si1",i,j,SNOW,1d3*SSIL(1)
c     $       /(XSI(1)*(ACE1I+SNOW)),1d3*SSIL(2)/(XSI(2)*(ACE1I+SNOW))
c     $       ,1d3*(SSIL(1)+SSIL(2))/ACE1I
c        if (debug) print*,"si1",i,j,SNOW,HSIL,SSIL,MSI2,FMOC,FHOC,FSOC,
c     *       RUN,ERUN,SRUN,WETSNOW,MELT12,SNOW
c     $       +ACE1I-SSIL(1)-SSIL(2) !, TRSIL(1,1)+TRSIL(1,2)
C**** Decay sea ice salinity
        IF(DOMAIN.EQ.'OCEAN') THEN
          CALL SSIDEC(SNOW,MSI2,HSIL,SSIL,DTsrc,MELT12,
#ifdef TRACERS_WATER
     *         TRSIL,TRFLUX,
#endif
     *         MFLUX,HFLUX,SFLUX,debug)
        else
          MFLUX=0. ; SFLUX=0. ; HFLUX=0.
#ifdef TRACERS_WATER
          TRFLUX = 0.
#endif
        end if

c      if (debug) print*,"gs4",trsil(1,1)-(snow+ace1i)*xsi(1)+
c     &             ssil(1), trsil(1,2)-(snow+ace1i)
c     *     *xsi(2)+ssil(2),trsil(1,3)-msi2*xsi(3)+
c     *     ssil(3),trsil(1,4)-msi2*xsi(4)+ssil(4),TRFLUX(1)-MFLUX+SFLUX
c        if (debug)  write(6,'(A,2I4,4F11.6)') "si2",i,j,SNOW,1d3*SSIL(1)
c     $       /(XSI(1)*(ACE1I+SNOW)),1d3*SSIL(2)/(XSI(2)*(ACE1I+SNOW))
c     $       ,1d3*(SSIL(1)+SSIL(2))/ACE1I
c        if (debug) print*,"si2",i,j,SNOW,HSIL,SSIL,MSI2,MFLUX,HFLUX,
c     *       SFLUX,SNOW+ACE1I-SSIL(1)-SSIL(2) !, TRSIL(1,1)+TRSIL(1,2)
C**** Calculate snow-ice possibility
        IF(snow_ice .eq. 1 .and. DOMAIN.EQ.'OCEAN') THEN
          call snowice(Tm,Sm,SNOW,MSI2,HSIL,SSIL,qsfix,
#ifdef TRACERS_WATER
     *         Trm,Tralpha,TRSIL,TRSNWIC,
#endif
     *         MSNWIC,HSNWIC,SSNWIC,DSNOW)
        else
          MSNWIC=0. ; SSNWIC=0. ; HSNWIC=0. ; DSNOW=0.
#ifdef TRACERS_WATER
          TRSNWIC = 0.
#endif
        end if
c     if (debug)  write(6,'(A,2I4,4F11.6)') "si3",i,j,SNOW,1d3*SSIL(1)
c     $       /(XSI(1)*(ACE1I+SNOW)),1d3*SSIL(2)/(XSI(2)*(ACE1I+SNOW))
c     $       ,1d3*(SSIL(1)+SSIL(2))/ACE1I
c        if (debug) print*,"si3",i,j,SNOW,HSIL,SSIL,MSI2, MSNWIC,HSNWIC,
c     *       SSNWIC,SNOW+ACE1I-SSIL(1)-SSIL(2) ! ,TRSIL(1,1)+TRSIL(1,2),

C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J)=SNOW
        HSI(:,I,J)=HSIL(:)
        SSI(:,I,J)=SSIL(:)
        MSI(I,J) = MSI2
#ifdef TRACERS_WATER
        TRSI(:,:,I,J) = TRSIL(:,:)
#endif
        FLAG_DSWS(I,J)=WETSNOW

        MSI1 = ACE1I + SNOW
        IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
          MICE1 = ACE1I-XSI(2)*MSI1
          SNOWL1= SNOW
        ELSE  ! some snow in second layer
          MICE1 = 0.
          SNOWL1= XSI(1)*MSI1
        ENDIF
        IF (MICE1.NE.0.) THEN
        Ti1 = Ti2b(HSIL(1)/(XSI(1)*MSI1),1d3*SSIL(1)/MICE1,SNOWL1,MICE1)
        ELSE
        Ti1 = Ti(HSIL(1)/(XSI(1)*(SNOW+ACE1I)),0d0)
        ENDIF

        TI1save(I,J) = Ti1

        SIHC(I,J) = SUM(HSIL(:))
        SNOWsave(I,J) = SNOW

C**** pond_melt accumulation
        pond_melt(i,j)=pond_melt(i,j)+0.3d0*MELT12
        pond_melt(i,j)=MIN(pond_melt(i,j),0.5*(MSI2+SNOW+ACE1I))

C**** decay is slow if there is some melting, faster otherwise
        if (MELT12.gt.0) then   ! 30 day decay
          pond_melt(i,j)=pond_melt(i,j)*(1.-dtsrc/(30.*SECONDS_PER_DAY))
        else                    ! 10 day decay
          pond_melt(i,j)=pond_melt(i,j)*(1.-dtsrc/(10.*SECONDS_PER_DAY))
        end if

C**** saftey valve to ensure that melt ponds eventually disappear (Ti<-10)
        if (Ti1 .lt.-10.) pond_melt(i,j)=0.  ! refreeze

c        if (debug) write(6,'(A,4I4,4F11.6)') "ponds",i,j,jday,jhour,
c     *       melt12,pond_melt(i,j),ti(HSIL(1)/(XSI(1)*(SNOW+ACE1I)),
c     *         1d3*SSIL(1)/(XSI(1)*(SNOW+ACE1I))),
c     *         1d3*(SSIL(1)+SSIL(2))/ACE1I
        
C**** Net fluxes to ocean
        RUNOSI(I,J) = FMOC + RUN  + MFLUX + MSNWIC
        ERUNOSI(I,J)= FHOC + ERUN + HFLUX + HSNWIC
        SRUNOSI(I,J)= FSOC + SRUN + SFLUX + SSNWIC
        SOLAR_IO(I,J)= SROX(2)
#ifdef TRACERS_WATER
        TRUNOSI(:,I,J) = FTROC(:) + TRRUN(:) + TRFLUX(:) + TRSNWIC(:)
#endif
        SNTOSI(I,J) = POICE*(DSNOW-MSNWIC + CMPRS)
        SITOPMLT(I,J) = POICE*(RUN+MFLUX)
        MSNFLOOD(I,J) = -POICE*MSNWIC
        HSNFLOOD(I,J) = -POICE*HSNWIC

      END IF

      END DO
      END DO

C****
      call stopTimer('GROUND_SI()')
      END SUBROUTINE GROUND_SI

      SUBROUTINE FORM_SI(si_state,iceocn,atmice)
!@sum  FORM_SI driver for adding new sea ice
!@auth Original Development team
!@calls seaice:addice
      USE CONSTANT, only : tf
      USE MODEL_COM, only : kocean
      USE SEAICE, only : ace1i,addice,lmi,fleadoc,fleadlk,xsi,debug
      USE SEAICE_COM, only : icestate
      USE EXCHANGE_TYPES, only : iceocn_xchng_vars,atmice_xchng_vars
#ifdef TRACERS_WATER
      USE SEAICE, only : ntm
#endif
      IMPLICIT NONE
      type(icestate) :: si_state
      type(iceocn_xchng_vars) :: iceocn
      type(atmice_xchng_vars) :: atmice
      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 SNOW,ROICE,MSI2,ENRGFO,ACEFO,ACEFI,ENRGFI,SALTO,SALTI
     *     ,POICE,PWATER,FLEAD,POCEAN,DMIMP,DHIMP,DSIMP
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: trsil
      REAL*8, DIMENSION(NTM) :: tro,tri,dtrimp
#endif
      CHARACTER(LEN=8) :: DOMAIN
      LOGICAL QFIXR,DOPOINT
      INTEGER I,J,N
      real*8, dimension(:,:), pointer :: fwater
      real*8, dimension(:,:), pointer :: ticesave,ssi1save,ssi2save
      real*8, dimension(:,:,:), pointer :: dmsi,dhsi,dssi,aij
      real*8, dimension(:,:), pointer :: rsi,msi,snowi
      real*8, dimension(:,:,:), pointer :: hsi,ssi
#ifdef TRACERS_WATER
      real*8, dimension(:,:,:,:), pointer :: trsi,dtrsi
#endif

      integer :: J_0, J_1, I_0,I_1

      I_0 = si_state%I_0
      I_1 = si_state%I_1
      J_0 = si_state%J_0
      J_1 = si_state%J_1

      domain = si_state%domain

      debug=.false.

c
c set pointers
c
      rsi => si_state%rsi
      msi => si_state%msi
      hsi => si_state%hsi
      ssi => si_state%ssi
      snowi => si_state%snowi
#ifdef TRACERS_WATER
      trsi => si_state%trsi
#endif
      fwater => iceocn%fwater
      dmsi => iceocn%dmsi
      dhsi => iceocn%dhsi
      dssi => iceocn%dssi
#ifdef TRACERS_WATER
      dtrsi => iceocn%dtrsi
#endif

      ticesave => atmice%ticesave
      SSI1save => atmice%SSI1save
      SSI2save => atmice%SSI2save
      aij => atmice%aij

      DO J=J_0, J_1
      DO I=I_0,si_state%IMAXJ(J)

c        debug = i.eq.52.and.j.eq.87
        
      PWATER=FWATER(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      POCEAN=(1.-ROICE)*PWATER
      DOPOINT = PWATER.GT.0.
      IF(DOPOINT) THEN ! todo: no more dopoint
c      IF (PWATER.gt.0) THEN
        SNOW= SNOWI(I,J)      ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
        HSIL(:) = HSI(:,I,J)      ! sea ice enthalpy
        SSIL(:) = SSI(:,I,J)      ! sea ice salt
#ifdef TRACERS_WATER
        TRSIL(:,:)= TRSI(:,:,I,J)
#endif

        IF (DOMAIN.EQ.'OCEAN') THEN ! todo: move outside loop
          FLEAD=FLEADOC
          IF (KOCEAN.ge.1) THEN
            QFIXR=.FALSE.
          ELSE
            QFIXR=.TRUE.
          END IF
        ELSE
          FLEAD=FLEADLK
          QFIXR=.FALSE.
        END IF

        ACEFO=DMSI(1,I,J)
        ACEFI=DMSI(2,I,J)
        ENRGFO=DHSI(1,I,J)
        ENRGFI=DHSI(2,I,J)
        SALTO=DSSI(1,I,J)
        SALTI=DSSI(2,I,J)
#ifdef TRACERS_WATER
        TRO(:) = DTRSI(:,1,I,J)
        TRI(:) = DTRSI(:,2,I,J)
#endif

        CALL ADDICE (SNOW,ROICE,HSIL,SSIL,MSI2,TSIL,ENRGFO,ACEFO,ACEFI,
     *       ENRGFI,SALTO,SALTI,
#ifdef TRACERS_WATER
     *       TRSIL,TRO,TRI,DTRIMP,
#endif
     *       DMIMP,DHIMP,DSIMP,FLEAD,QFIXR,debug)

C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J) = SNOW
        MSI(I,J)=MSI2
        HSI(:,I,J) = HSIL(:)
        SSI(:,I,J) = SSIL(:)
#ifdef TRACERS_WATER
        TRSI(:,:,I,J) = TRSIL(:,:)
#endif
#ifndef STANDALONE_OCEAN
        IF (.not. QFIXR) THEN
          IF(DOMAIN.EQ.'LAKES') THEN
            if (roice.gt.rsi(i,j)) ! ice from ocean
     *         call RESET_SURF_FLUXES(I,J,1,2,RSI(I,J),ROICE)
            if (roice.lt.rsi(i,j)) ! ocean from ice
     *         call RESET_SURF_FLUXES(I,J,2,1,1.-RSI(I,J),1.-ROICE)
          ENDIF
        ELSE
C**** save implicit mass-flux diagnostics
          AIJ(I,J,atmice%IJ_SMFX)=AIJ(I,J,atmice%IJ_SMFX)+ROICE*DMIMP
          CALL INC_AJ(I,J,atmice%ITOICE,atmice%J_IMPLM,
     &         -(DMIMP-DSIMP)*POICE)
          CALL INC_AJ(I,J,atmice%ITOICE,atmice%J_IMPLH,
     &                -DHIMP *POICE)
        END IF
#endif
        IF (.not. QFIXR) RSI(I,J)=ROICE

        TICEsave(I,J)=(XSI(3)*TSIL(3)+XSI(4)*TSIL(4))
        SSI1save(I,J)=(SSI(1,I,J)+SSI(2,I,J))/ACE1I
        SSI2save(I,J)=(SSI(3,I,J)+SSI(4,I,J))/MSI(I,J)
        atmice%SNOWsave2(I,J)=SNOWI(I,J)
        atmice%MSIsave2(I,J)=MSI(I,J)
#ifdef TRACERS_WATER
        atmice%TRSIsum(:,I,J) = sum(TRSI(:,:,I,J),dim=2)
#endif

      if (TSIL(1).lt.-100.) then
         write(6,*) "Seaice: T < -100. i,j,TSI = ",i,j,TSIL(1:LMI)
         call stop_model("Seaice too cold after ADDICE",255)
      end if

      END IF

      END DO
      END DO

C**** replicate ice values at the poles
      IF(DOMAIN.EQ.'OCEAN') THEN
      IF (si_state%HAVE_NORTH_POLE) THEN
        DO I=2,I_1
          RSI(I,J_1)=RSI(1,J_1)
          MSI(I,J_1)=MSI(1,J_1)
          HSI(:,I,J_1)=HSI(:,1,J_1)
          SSI(:,I,J_1)=SSI(:,1,J_1)
          SNOWI(I,J_1)=SNOWI(1,J_1)
#ifdef TRACERS_WATER
          TRSI(:,:,I,J_1) = TRSI(:,:,1,J_1)
#endif
        END DO
      END IF
      IF (si_state%HAVE_SOUTH_POLE) THEN
        DO I=2,I_1
          RSI(I,1)=RSI(1,1)
          MSI(I,1)=MSI(1,1)
          HSI(:,I,1)=HSI(:,1,1)
          SSI(:,I,1)=SSI(:,1,1)
          SNOWI(I,1)=SNOWI(1,1)
#ifdef TRACERS_WATER
          TRSI(:,:,I,1) = TRSI(:,:,1,1)
#endif
        END DO
      END IF
      END IF
C****
      END SUBROUTINE FORM_SI


#ifndef STANDALONE_OCEAN
      SUBROUTINE SI_diags(si_state,iceocn,atmice)
      USE MODEL_COM, only : kocean,itime
      USE CONSTANT, only : rhows,rhow,bylhm
      USE EXCHANGE_TYPES, only : atmice_xchng_vars,iceocn_xchng_vars
      USE SEAICE_COM, only : icestate
      USE SEAICE, only : rhos,ace1i
#ifdef TRACERS_WATER
      USE SEAICE, only : ntm
      !USE TRACER_COM, only : itime_tr0,tr_wd_type,nWater,nPART
#endif
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      IMPLICIT NONE
      type(icestate) :: si_state
      type(iceocn_xchng_vars) :: iceocn
      type(atmice_xchng_vars) :: atmice
      REAL*8 POCEAN,ROICE,POICE,PWATER,SCOVI
      REAL*8 :: IMLT,HMLT,SMLT
      INTEGER I,J,JR,ITYPE,ITYPEO,N
      LOGICAL DOMELT
      CHARACTER(LEN=8) :: DOMAIN
      real*8, dimension(:,:), pointer ::
     &     runpsi,srunpsi,erunpsi,
     &     runosi,srunosi,erunosi,
     &     melti,emelti,smelti,
     &     fmsi_io,fhsi_io,fssi_io,
     &     prec,eprec,fwater
     &     ,RSIstart,MSIsave,SNTOSI,SITOPMLT,MSNFLOOD,HSNFLOOD
     &     ,ticesave,Ti1save,SNOWsave,ssi1save,ssi2save,SIHC
      real*8, dimension(:,:,:), pointer :: dmsi,dhsi,dssi,aij
      real*8, dimension(:,:), pointer :: msi,snowi,pond_melt
#ifdef TRACERS_WATER
      real*8, dimension(:,:,:,:), pointer :: trsi,dtrsi,taijn
      real*8, dimension(:,:,:), pointer ::
     &     trunpsi,ftrsi_io,trunosi,trmelti
#endif
      integer :: J_0, J_1, I_0,I_1

      call startTimer('Diagnostics')

      domain = si_state%domain

      I_0 = si_state%I_0
      I_1 = si_state%I_1
      J_0 = si_state%J_0
      J_1 = si_state%J_1

c
c set pointers
c
      msi => atmice%msisave2 !si_state%msi
      snowi => atmice%snowsave2 !si_state%snowi
      pond_melt => si_state%pond_melt
#ifdef TRACERS_WATER
      trsi => si_state%trsi
#endif
       fwater => iceocn%fwater
       runpsi => iceocn%runpsi
      srunpsi => iceocn%srunpsi
      erunpsi => iceocn%erunpsi
      fmsi_io => iceocn%fmsi_io
      fhsi_io => iceocn%fhsi_io
      fssi_io => iceocn%fssi_io
       melti => iceocn%melti
      emelti => iceocn%emelti
      smelti => iceocn%smelti
       runosi => iceocn%runosi
      srunosi => iceocn%srunosi
      erunosi => iceocn%erunosi
      dmsi => iceocn%dmsi
      dhsi => iceocn%dhsi
      dssi => iceocn%dssi
       prec => atmice%prec
      eprec => atmice%eprec
#ifdef TRACERS_WATER
      trunosi => iceocn%trunosi
      trunpsi => iceocn%trunpsi
      trmelti => iceocn%trmelti
      dtrsi => iceocn%dtrsi
      ftrsi_io => iceocn%ftrsi_io
      taijn => atmice%taijn
#endif

      RSIstart => atmice%RSIstart
      ticesave => atmice%ticesave
      MSIsave  => atmice%MSIsave
      SSI1save => atmice%SSI1save
      SSI2save => atmice%SSI2save
      SNTOSI   => atmice%SNTOSI
      SITOPMLT => atmice%SITOPMLT
      MSNFLOOD => atmice%MSNFLOOD
      HSNFLOOD => atmice%HSNFLOOD
      TI1save  => atmice%TI1save
      SIHC     => atmice%SIHC
      SNOWsave => atmice%SNOWsave

      aij => atmice%aij

      IF(DOMAIN.EQ.'OCEAN') THEN
        DO J=J_0, J_1
        DO I=I_0, I_1
          IF(atmice%FOCEAN(I,J).LE.0.) cycle
          AIJ(I,J,ATMICE%IJ_MUSI)=AIJ(I,J,ATMICE%IJ_MUSI)
     &         +ATMICE%MUSI(I,J)
          AIJ(I,J,ATMICE%IJ_HUSI)=AIJ(I,J,ATMICE%IJ_HUSI)
     &         +ATMICE%HUSI(I,J)
          AIJ(I,J,ATMICE%IJ_SUSI)=AIJ(I,J,ATMICE%IJ_SUSI)
     &         +ATMICE%SUSI(I,J)
          AIJ(I,J,ATMICE%IJ_MVSI)=AIJ(I,J,ATMICE%IJ_MVSI)
     &         +ATMICE%MVSI(I,J)
          AIJ(I,J,ATMICE%IJ_HVSI)=AIJ(I,J,ATMICE%IJ_HVSI)
     &         +ATMICE%HVSI(I,J)
          AIJ(I,J,ATMICE%IJ_SVSI)=AIJ(I,J,ATMICE%IJ_SVSI)
     &         +ATMICE%SVSI(I,J)
          AIJ(I,J,ATMICE%IJ_dHSI_Dyn) = AIJ(I,J,ATMICE%IJ_dHSI_Dyn) +
     +          ATMICE%HSICNV(I,J)
        ENDDO
        ENDDO
#ifdef TRACERS_WATER
        DO N=1,NTM
        DO J=J_0, J_1
        DO I=I_0, I_1
          IF(atmice%FOCEAN(I,J).LE.0.) cycle
          TAIJN(I,J,atmice%TIJ_TUSI,N)=TAIJN(I,J,atmice%TIJ_TUSI,N)
     &         +ATMICE%TUSI(I,J,N)
          TAIJN(I,J,atmice%TIJ_TVSI,N)=TAIJN(I,J,atmice%TIJ_TVSI,N)
     &         +ATMICE%TVSI(I,J,N)
        ENDDO
        ENDDO
        ENDDO
#endif
      ENDIF

      IF(DOMAIN.EQ.'OCEAN') THEN
        ITYPE=atmice%ITOICE
        ITYPEO=atmice%ITOCEAN
        DOMELT=KOCEAN.GE.1
      ELSEIF(DOMAIN.EQ.'LAKES') THEN
        ITYPE=atmice%ITLKICE
        ITYPEO=atmice%ITLAKE
        DOMELT=.TRUE.
      ENDIF

      DO J=J_0, J_1
      DO I=I_0,si_state%IMAXJ(J)
        PWATER=FWATER(I,J)
        IF(PWATER.LE.0.) CYCLE
        ROICE=RSIstart(I,J)
        POCEAN=(1.-ROICE)*PWATER
        POICE=ROICE*PWATER
        JR=atmice%JREG(I,J)

C**** MELT_SI diags
        IF(DOMELT .and. MELTI(I,J).NE.0.) THEN
          AIJ(I,J,atmice%IJ_SIGRLT)=AIJ(I,J,atmice%IJ_SIGRLT)-MELTI(I,J)
          IF (DOMAIN.EQ.'OCEAN') THEN
            AIJ(I,J,ATMICE%IJ_FWIO)=AIJ(I,J,ATMICE%IJ_FWIO)
     &           +(MELTI(I,J)-SMELTI(I,J))
            AIJ(I,J,ATMICE%IJ_HTIO)=AIJ(I,J,ATMICE%IJ_HTIO)+EMELTI(I,J)
            AIJ(I,J,ATMICE%IJ_STIO)=AIJ(I,J,ATMICE%IJ_STIO)+SMELTI(I,J)
#ifdef TRACERS_WATER
            TAIJN(I,J,atmice%TIJ_ICOCFLX,:)=
     &           TAIJN(I,J,atmice%TIJ_ICOCFLX,:) +TRMELTI(:,I,J)
#endif
          END IF

          CALL INC_AJ(I,J,ITYPE,ATMICE%J_HMELT,EMELTI(I,J)*ROICE)
          CALL INC_AJ(I,J,ITYPE,ATMICE%J_SMELT,SMELTI(I,J)*ROICE)
          CALL INC_AJ(I,J,ITYPE,ATMICE%J_IMELT, MELTI(I,J)*ROICE)
          CALL INC_AJ(I,J,ITYPEO,ATMICE%J_HMELT,EMELTI(I,J)*(1.-ROICE))
          CALL INC_AJ(I,J,ITYPEO,ATMICE%J_SMELT,SMELTI(I,J)*(1.-ROICE))
          CALL INC_AJ(I,J,ITYPEO,ATMICE%J_IMELT,MELTI(I,J)*(1.-ROICE))

          CALL INC_AREG(I,J,JR,ATMICE%J_HMELT,EMELTI(I,J))
          CALL INC_AREG(I,J,JR,ATMICE%J_SMELT,SMELTI(I,J))
          CALL INC_AREG(I,J,JR,ATMICE%J_IMELT, MELTI(I,J))

        END IF ! DOMELT .and. MELTI(I,J).NE.0.

C****
C**** Open water ice formation diagnostics
C****
        IF(POCEAN*DMSI(1,I,J).NE.0.) THEN

C**** ice formation diagnostics on the atmospheric grid
! define open ocean ice formation as frazil ice growth
          AIJ(I,J,ATMICE%IJ_SIGRFR)=AIJ(I,J,ATMICE%IJ_SIGRFR)
     &         +POCEAN*DMSI(1,I,J)

          IF (DOMAIN.EQ.'OCEAN') THEN
            AIJ(I,J,ATMICE%IJ_FWIO)=AIJ(I,J,ATMICE%IJ_FWIO)
     &           -POCEAN*(DMSI(1,I,J)-DSSI(1,I,J))
            AIJ(I,J,ATMICE%IJ_HTIO)=AIJ(I,J,ATMICE%IJ_HTIO)
     &           -POCEAN* DHSI(1,I,J)
            AIJ(I,J,ATMICE%IJ_STIO)=AIJ(I,J,ATMICE%IJ_STIO)
     &           -POCEAN* DSSI(1,I,J)
#ifdef TRACERS_WATER
            TAIJN(I,J,atmice%TIJ_ICOCFLX,:)=
     &           TAIJN(I,J,atmice%TIJ_ICOCFLX,:)-POCEAN*DTRSI(:,1,I,J)
#endif
          END IF

          CALL INC_AJ(I,J,ITYPEO,ATMICE%J_IMELT,-DMSI(1,I,J)*POCEAN)
          CALL INC_AJ(I,J,ITYPEO,ATMICE%J_HMELT,-DHSI(1,I,J)*POCEAN)
          CALL INC_AJ(I,J,ITYPEO,ATMICE%J_SMELT,-DSSI(1,I,J)*POCEAN)

C**** regional diagnostics
          CALL INC_AREG(I,J,JR,ATMICE%J_IMELT,-DMSI(1,I,J)*POCEAN)
          CALL INC_AREG(I,J,JR,ATMICE%J_HMELT,-DHSI(1,I,J)*POCEAN)
          CALL INC_AREG(I,J,JR,ATMICE%J_SMELT,-DSSI(1,I,J)*POCEAN)

        ENDIF ! POCEAN*DMSI(1,I,J).NE.0.

C**** Ice-covered ocean diagnostics

        IF(POICE.GT.0.) THEN

! define under ice formation as congelation ice growth
          AIJ(I,J,ATMICE%IJ_SIGRCG)=AIJ(I,J,ATMICE%IJ_SIGRCG)
     &         +POICE*DMSI(2,I,J)

          AIJ(I,J,ATMICE%IJ_SSI1)=AIJ(I,J,ATMICE%IJ_SSI1)
     &         +POICE*SSI1save(I,J)
          AIJ(I,J,ATMICE%IJ_SSI2)=AIJ(I,J,ATMICE%IJ_SSI2)
     &         +POICE*SSI2save(I,J)
          AIJ(I,J,ATMICE%IJ_TSI)=AIJ(I,J,ATMICE%IJ_TSI)
     &         +POICE*TICEsave(I,J)

#ifdef TRACERS_WATER
C**** Save sea ice tracer amount
          do n=1,ntm
c            if (itime_tr0(n).le.itime .and.
c     &       (tr_wd_TYPE(n).eq.nWater .or. tr_wd_TYPE(n).eq.nPART)) then
            if(atmice%do_accum(n)) then
              taijn(i,j,atmice%tij_seaice,n)=
     &             taijn(i,j,atmice%tij_seaice,n)+
     &             POICE*atmice%TRSIsum(n,I,J)
            end if
          end do
#endif

          AIJ(I,J,ATMICE%IJ_F0OI) = AIJ(I,J,ATMICE%IJ_F0OI)
     &         +EPREC(I,J)*POICE
          AIJ(I,J,ATMICE%IJ_SISNWF) = AIJ(I,J,ATMICE%IJ_SISNWF)
     &         -MIN(EPREC(I,J)*BYLHM,0d0)*POICE
          AIJ(I,J,ATMICE%IJ_RSOI) =AIJ(I,J,ATMICE%IJ_RSOI) +POICE
          IF(POICE.GT.0d0)
C**** SITF: Fraction of time steps of the averaging period during which
c**** sea ice is present (siconc >0 ) in a grid cell
     &    AIJ(I,J,ATMICE%IJ_SITF) =AIJ(I,J,ATMICE%IJ_SITF) + 1
          AIJ(I,J,ATMICE%IJ_MSI) =AIJ(I,J,ATMICE%IJ_MSI) + MSIsave(I,J)
          AIJ(I,J,ATMICE%IJ_SIMASS) =AIJ(I,J,ATMICE%IJ_SIMASS)
     *         + MSIsave(I,J)
          AIJ(I,J,ATMICE%IJ_SIVOL) =AIJ(I,J,ATMICE%IJ_SIVOL)
     *         + MSIsave(I,J)

          AIJ(I,J,ATMICE%IJ_SITOPMLT)=AIJ(I,J,ATMICE%IJ_SITOPMLT)
     &         +RUNPSI(I,J)*POICE

          IMLT =  RUNPSI(I,J)+ RUNOSI(I,J)-DMSI(2,I,J)
          HMLT = ERUNPSI(I,J)+ERUNOSI(I,J)-DHSI(2,I,J)
          SMLT = SRUNPSI(I,J)+SRUNOSI(I,J)-DSSI(2,I,J)

          IF (DOMAIN.EQ.'OCEAN') THEN
            AIJ(I,J,ATMICE%IJ_FWIO)=AIJ(I,J,ATMICE%IJ_FWIO)
     &           +(IMLT-SMLT)*POICE
            AIJ(I,J,ATMICE%IJ_HTIO)=AIJ(I,J,ATMICE%IJ_HTIO) +HMLT*POICE
            AIJ(I,J,ATMICE%IJ_STIO)=AIJ(I,J,ATMICE%IJ_STIO) +SMLT*POICE
#ifdef TRACERS_WATER
            TAIJN(I,J,atmice%TIJ_ICOCFLX,:)=
     &           TAIJN(I,J,atmice%TIJ_ICOCFLX,:)
     &           +(TRUNPSI(:,I,J)+TRUNOSI(:,I,J)-DTRSI(:,2,I,J))*POICE
#endif
          END IF

C**** snow cover diagnostic now matches that seen by the radiation
          IF (SNOWsave(I,J).GT.0) THEN
            SCOVI=MIN(1d0,SNOWsave(I,J)/(RHOS*0.1d0))*POICE
          ELSE
            SCOVI=0.
          ENDIF

          AIJ(I,J,ATMICE%IJ_RSNW)=AIJ(I,J,ATMICE%IJ_RSNW)+SCOVI
          AIJ(I,J,ATMICE%IJ_SNOW)=AIJ(I,J,ATMICE%IJ_SNOW)
     &         +SNOWsave(I,J)*POICE
          AIJ(I,J,ATMICE%IJ_RSIT)=AIJ(I,J,ATMICE%IJ_RSIT)+POICE
          AIJ(I,J,ATMICE%IJ_MLTP)=AIJ(I,J,ATMICE%IJ_MLTP)
     &         +pond_melt(i,j)*POICE
          AIJ(I,J,ATMICE%IJ_ZSNOW)=AIJ(I,J,ATMICE%IJ_ZSNOW)
     &         +POICE*SNOWsave(I,J)/RHOS
          AIJ(I,J,ATMICE%IJ_TSICE)=AIJ(I,J,ATMICE%IJ_TSICE)
     &         +Ti1save(I,J)*POICE
          AIJ(I,J,ATMICE%IJ_SIHC)=AIJ(I,J,ATMICE%IJ_SIHC)
     &         +SIHC(I,J)*POICE
          AIJ(I,J,ATMICE%IJ_SNTOSI)=AIJ(I,J,ATMICE%IJ_SNTOSI)
     &         +SNTOSI(I,J)
          AIJ(I,J,ATMICE%IJ_SITOPMLT)=AIJ(I,J,ATMICE%IJ_SITOPMLT)
     &         +SITOPMLT(I,J)
          AIJ(I,J,ATMICE%IJ_MSNFLOOD)=AIJ(I,J,ATMICE%IJ_MSNFLOOD)
     &         +MSNFLOOD(I,J)
          AIJ(I,J,ATMICE%IJ_HSNFLOOD)=AIJ(I,J,ATMICE%IJ_HSNFLOOD)
     &         +HSNFLOOD(I,J)

          IF (fmsi_io(i,j).lt.0) THEN   ! define as congelation growth
            AIJ(I,J,ATMICE%IJ_SIGRCG)=AIJ(I,J,ATMICE%IJ_SIGRCG)
     &           -POICE*fmsi_io(i,j)
          ELSE                  ! basal melt
            AIJ(I,J,ATMICE%IJ_SIBOTMLT)=AIJ(I,J,ATMICE%IJ_SIBOTMLT)
     &           +POICE*fmsi_io(i,j)
          END IF

          CALL INC_AJ(I,J,ITYPE,atmice%J_RSNOW,SCOVI)
          CALL INC_AJ(I,J,ITYPE,atmice%J_RSI ,      POICE)
          CALL INC_AJ(I,J,ITYPE,atmice%J_ACE1, ACE1I*POICE)
          CALL INC_AJ(I,J,ITYPE,atmice%J_ACE2, MSI(I,J)*POICE)
          CALL INC_AJ(I,J,ITYPE,atmice%J_SNOW, SNOWI(I,J)*POICE)

          CALL INC_AJ(I,J,ITYPE,ATMICE%J_IMELT, IMLT*POICE)
          CALL INC_AJ(I,J,ITYPE,ATMICE%J_HMELT, HMLT*POICE)
          CALL INC_AJ(I,J,ITYPE,ATMICE%J_SMELT, SMLT*POICE)

C**** Accumulate regional diagnostics
          IF (JR.ne.24) THEN
            CALL INC_AREG(I,J,JR,atmice%J_RSI ,      POICE)
            CALL INC_AREG(I,J,JR,atmice%J_SNOW, SNOWI(I,J)*POICE)
            CALL INC_AREG(I,J,JR,atmice%J_ACE1, ACE1I*POICE)
            CALL INC_AREG(I,J,JR,atmice%J_ACE2, MSI(I,J)*POICE)
          END IF
          CALL INC_AREG(I,J,JR,atmice%J_RSNOW,SCOVI)

          CALL INC_AREG(I,J,JR,ATMICE%J_IMELT, IMLT*POICE)
          CALL INC_AREG(I,J,JR,ATMICE%J_HMELT, HMLT*POICE)
          CALL INC_AREG(I,J,JR,ATMICE%J_SMELT, SMLT*POICE)

        END IF ! POICE.GT.0.

      END DO
      END DO

      call stopTimer('Diagnostics')

      RETURN
      END SUBROUTINE SI_diags
#endif

      subroutine set_noice_defaults(si_state,iceocn)
!@sum  set_noice_defaults sets defaults for ice-free conditions
!@auth Original Development Team
      USE SEAICE, only : xsi,ace1i,ac2oim,ssi0,Ti,Ei
      USE SEAICE_COM, only : icestate
      USE EXCHANGE_TYPES, only : iceocn_xchng_vars
      IMPLICIT NONE
      type(icestate) :: si_state
      type(iceocn_xchng_vars) :: iceocn
      INTEGER I,J
      REAL*8 TFO,SAL
      integer :: I_0, I_1, J_0, J_1

      I_0 = si_state%I_0
      I_1 = si_state%I_1
      J_0 = si_state%J_0
      J_1 = si_state%J_1

      if(si_state%domain.eq.'OCEAN') then
        SAL = SSI0
        TFO = -1.87d0       ! reasonable value, doesn't really matter
      else
        SAL = 0.
        TFO = 0.
      endif

      DO J=J_0, J_1
      DO I=I_0, I_1
        IF(si_state%RSI(I,J).le.0) THEN
          si_state%MSI(I,J) = AC2OIM
          si_state%SNOWI(I,J) = 0.
          if(iceocn%fwater(i,j).gt.0.) THEN
            si_state%SSI(1:2,I,J)=SAL*XSI(1:2)*ACE1I
            si_state%SSI(3:4,I,J)=SAL*XSI(3:4)*AC2OIM
            si_state%HSI(1:2,I,J)=Ei(TFO,1d3*SAL)*XSI(1:2)*ACE1I
            si_state%HSI(3:4,I,J)=Ei(TFO,1d3*SAL)*XSI(3:4)*AC2OIM
          else
            si_state%SSI(1:2,I,J)=0.
            si_state%SSI(3:4,I,J)=0.
            si_state%HSI(1:2,I,J)=Ei(0d0,0d0)*XSI(1:2)*ACE1I
            si_state%HSI(3:4,I,J)=Ei(0d0,0d0)*XSI(3:4)*AC2OIM
          end if
#ifdef TRACERS_WATER
          si_state%TRSI(:,:,I,J)=0.
#endif
          si_state%pond_melt(i,j) = 0.
          si_state%flag_dsws(i,j) = .FALSE.
        END IF
      END DO
      END DO

      return
      end subroutine set_noice_defaults

#ifdef TRACERS_WATER
      subroutine init_single_seaice_tracer(si_state,n,conc)
!@sum  
      USE SEAICE, only : xsi,ace1i
      USE SEAICE_COM, only : icestate
      IMPLICIT NONE
      type(icestate) :: si_state
      integer :: n
      real*8, intent(in) :: conc
      integer i,j
      do j=si_state%j_0,si_state%j_1
      do i=si_state%i_0,si_state%i_1
        if (si_state%msi(i,j).gt.0) then
          si_state%trsi(n,1,i,j)=conc*
     &         (xsi(1)*(si_state%snowi(i,j)+ace1i)-si_state%ssi(1,i,j))
          si_state%trsi(n,2,i,j)=conc*
     &         (xsi(2)*(si_state%snowi(i,j)+ace1i)-si_state%ssi(2,i,j))
          si_state%trsi(n,3,i,j)=conc*
     &         (xsi(3)*si_state%msi(i,j)-si_state%ssi(3,i,j))
          si_state%trsi(n,4,i,j)=conc*
     &         (xsi(4)*si_state%msi(i,j)-si_state%ssi(4,i,j))
        end if
      enddo
      enddo
      return
      end subroutine init_single_seaice_tracer
#endif

      SUBROUTINE init_oceanice(iniOCEAN,do_IC_fixups,atmocn)
!@sum  init_ice initialises ice arrays
!@auth Original Development Team
      USE CONSTANT, only : rhows,omega
      USE MODEL_COM, only : kocean, master_yr
      USE SEAICE, only : oi_ustar0,silmfac,snow_ice,seaice_thermo
      USE SEAICE_COM, only : si_ocn,iceocn
      USE Dictionary_mod
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      USE MODEL_COM, only : modelEclock
      USE TIMESTREAM_MOD, only : init_stream
      use pario, only : par_open, par_close, read_dist_data
      USE SEAICE_COM, only : RSIstream,ZSIstream
      USE SEAICE_COM, only : dm,rsi_exists,zsi_exists
      use seaice_com, only : grid=>sigrid
      use filemanager, only : file_exists
      IMPLICIT NONE
      LOGICAL :: iniOCEAN
      INTEGER do_IC_fixups
      type(atmocn_xchng_vars) :: atmocn
c
      INTEGER I,J
      integer :: I_0, I_1, J_0, J_1
      integer :: fid,jyear,jday,seaice_yr
      logical :: cyclic
      logical :: exists

      I_0 = iceocn%I_0
      I_1 = iceocn%I_1
      J_0 = iceocn%J_0
      J_1 = iceocn%J_1

      si_ocn%fwater = atmocn%focean
      iceocn%fwater = atmocn%focean
      deallocate(iceocn%ogeoza);   iceocn%ogeoza => atmocn%ogeoza

      DO J=J_0,J_1
      DO I=I_0,I_1
        iceocn%coriol(i,j) = ABS(2.*OMEGA*SIN(atmocn%LAT(I,J)))
      ENDDO
      ENDDO

C**** set up a default ice-ocean stress field. This can be changed by
C**** adjusting oi_ustar0 in the parameter list. If ice dynamics
C**** is used, this is overwritten.
      call sync_param("oi_ustar0",oi_ustar0)
      iceocn%UI2rho = rhows*(oi_ustar0)**2

C**** Adjust degree of lateral melt by changing silmfac
C**** Default is 2.5d-8, but could be changed by a factor of 2.
      call sync_param("silmfac",silmfac)

C**** Decide whether snow_ice formation is allowed
      call sync_param("snow_ice",snow_ice)

C**** Define the ice thermodynamics (SI or BP)
      call sync_param("seaice_thermo",seaice_thermo)

C**** clean up ice fraction/sea ice salinity possibly incorrect in I.C.
      if (do_IC_fixups == 1) then
        DO J=J_0, J_1
        DO I=I_0, I_1
          IF (si_ocn%FWATER(I,J).eq.0 .and. si_ocn%RSI(i,j).gt.0)
     &         si_ocn%RSI(I,J)=0
        END DO
        END DO
      end if

      IF (iniOCEAN) THEN
        si_ocn%rsix(:,:) = 0.
        si_ocn%rsiy(:,:) = 0.
        if(KOCEAN.EQ.0) call set_noice_defaults(si_ocn,iceocn)
      END IF

#ifndef STANDALONE_OCEAN
C**** Set conservation diagnostics for ice mass, energy, salt
      CALL DECLARE_SEAICE_CONSRV
#endif

      rsi_exists = (kocean.eq.0 .and. file_exists('SICE'))
      if(rsi_exists) then
        if(is_set_param('seaice_yr')) then
          ! If parameter seaice_yr exists, ice data from that year is
          ! selected (only relevant if SICE is a multi-year dataset).
          call get_param( 'seaice_yr', seaice_yr )
        else
          ! Otherwise, seaice_yr is set to ocean_yr or master_yr.
          call get_param( 'ocean_yr', seaice_yr, default=master_yr )
        endif
        cyclic = seaice_yr /= 0 ! seaice_yr==0 implies transient mode.
        seaice_yr = abs(seaice_yr)
        call modelEclock%get(year=jyear, dayOfYear=jday)
        if(cyclic) jyear = seaice_yr
        call init_stream(grid,RSIstream,'SICE','rsi',0d0,1d0,'ppm',
     &       jyear,jday,msk=atmocn%focean,cyclic=cyclic)
        inquire(file='ZSIFAC',exist=exists)
        if(exists) then
          zsi_exists = .false.
          fid = par_open(grid,'ZSIFAC','read')
          call read_dist_data(grid,fid,'dm',dm)
          call par_close(grid,fid)
        else
          zsi_exists = .true.
        endif
        if(zsi_exists) then
          call init_stream(grid,ZSIstream,'ZSI','ZSI',0d0,1d30,'ppm',
     &         jyear,jday,msk=atmocn%focean,cyclic=cyclic)
        endif
      endif

      END SUBROUTINE init_oceanice

      SUBROUTINE conserv_OMSI(ICE)
!@sum  conserv_MSI calculates total amount of snow and ice over ocean
!@auth Gavin Schmidt
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : si_ocn
      IMPLICIT NONE
!@var ICE total ocean snow and ice mass (kg/m^2)
      REAL*8, DIMENSION(si_ocn%I_0H:si_ocn%I_1H,
     &                  si_ocn%J_0H:si_ocn%J_1H) :: ICE
      INTEGER I,J
      INTEGER J_0, J_1, I_0,I_1

      I_0 = si_ocn%I_0
      I_1 = si_ocn%I_1
      J_0 = si_ocn%J_0
      J_1 = si_ocn%J_1

      DO J=J_0,J_1
      DO I=I_0,si_ocn%IMAXJ(J)
        ICE(I,J)=si_ocn%RSI(I,J)*
     &      (si_ocn%MSI(I,J)+ACE1I+si_ocn%SNOWI(I,J))*si_ocn%FWATER(I,J)
      END DO
      END DO
      IF (si_ocn%HAVE_SOUTH_POLE) ICE(2:I_1,  1)=ICE(1,  1)
      IF (si_ocn%HAVE_NORTH_POLE) ICE(2:I_1,J_1)=ICE(1,J_1)
      RETURN
C****
      END SUBROUTINE conserv_OMSI

      SUBROUTINE conserv_OHSI(EICE)
!@sum  conserv_HSI calculates total ice energy over ocean
!@auth Gavin Schmidt
      USE SEAICE_COM, only : si_ocn
      IMPLICIT NONE
!@var EICE total ocean snow and ice energy (J/m^2)
      REAL*8, DIMENSION(si_ocn%I_0H:si_ocn%I_1H,
     &                  si_ocn%J_0H:si_ocn%J_1H) :: EICE
      INTEGER I,J

      INTEGER J_0, J_1, I_0,I_1

      I_0 = si_ocn%I_0
      I_1 = si_ocn%I_1
      J_0 = si_ocn%J_0
      J_1 = si_ocn%J_1

      DO J=J_0,J_1
      DO I=I_0,si_ocn%IMAXJ(J)
        EICE(I,J)=si_ocn%RSI(I,J)*si_ocn%FWATER(I,J)*
     &       SUM(si_ocn%HSI(:,I,J))
      END DO
      END DO
      IF (si_ocn%HAVE_SOUTH_POLE) EICE(2:I_1,  1)=EICE(1,  1)
      IF (si_ocn%HAVE_NORTH_POLE) EICE(2:I_1,J_1)=EICE(1,J_1)
      RETURN
C****
      END SUBROUTINE conserv_OHSI

      SUBROUTINE conserv_OSSI(SALT)
!@sum  conserv_SSI calculates total amount of salt in ocean ice
!@auth Gavin Schmidt
      USE SEAICE_COM, only : si_ocn
      IMPLICIT NONE
!@var SALT total salt in ocean ice (kg/m^2)
      REAL*8, DIMENSION(si_ocn%I_0H:si_ocn%I_1H,
     &                  si_ocn%J_0H:si_ocn%J_1H) :: SALT
      INTEGER I,J
      INTEGER J_0, J_1, I_0,I_1

      I_0 = si_ocn%I_0
      I_1 = si_ocn%I_1
      J_0 = si_ocn%J_0
      J_1 = si_ocn%J_1

      DO J=J_0,J_1
      DO I=I_0,si_ocn%IMAXJ(J)
        IF (SI_OCN%FWATER(I,J).gt.0) THEN
          SALT(I,J)=SI_OCN%FWATER(I,J)*si_ocn%RSI(I,J)*
     &         SUM(si_ocn%SSI(:,I,J))
        ELSE
          SALT(I,J)=0
        END IF
      END DO
      END DO
      IF (si_ocn%HAVE_SOUTH_POLE) SALT(2:I_1,  1)=SALT(1,  1)
      IF (si_ocn%HAVE_NORTH_POLE) SALT(2:I_1,J_1)=SALT(1,J_1)
      RETURN
C****
      END SUBROUTINE conserv_OSSI

      SUBROUTINE seaice_to_atmgrid(atmice)
!@sum seaice_to_atmgrid set sea ice properties on the atm grid
!@auth Gavin Schmidt
      USE CONSTANT, only : tf
#ifdef SCM
      USE SCM_COM, only : SCMopt,SCMin
#endif
      USE SEAICE_COM, only : si_atm,si_ocn
      USE SEAICE, only : ace1i,xsi,lmi,Ti,rhoi,rhos,Ti2b
      USE EXCHANGE_TYPES, only : atmice_xchng_vars
      IMPLICIT NONE
      type(atmice_xchng_vars) :: atmice
c
      INTEGER I,J, J_0, J_1 ,I_0,I_1
      REAL*8 MSI1,SNOWL(2),MICE(2)

      I_0 = atmice%I_0
      I_1 = atmice%I_1
      J_0 = atmice%J_0
      J_1 = atmice%J_1

      DO J=J_0, J_1
      DO I=I_0, I_1
        IF(atmice%FOCEAN(I,J).GT.0.) THEN
          ! while ocean ice still on atm grid:
          si_atm%rsi(i,j) = si_ocn%rsi(i,j)
          si_atm%snowi(i,j) = si_ocn%snowi(i,j)
          si_atm%msi(i,j) = si_ocn%msi(i,j)
          si_atm%pond_melt(i,j) = si_ocn%pond_melt(i,j)
          si_atm%flag_dsws(i,j) = si_ocn%flag_dsws(i,j)
          si_atm%hsi(:,i,j) = si_ocn%hsi(:,i,j)
          si_atm%ssi(:,i,j) = si_ocn%ssi(:,i,j)
#ifdef TRACERS_WATER
          si_atm%trsi(:,:,i,j) = si_ocn%trsi(:,:,i,j)
#endif
        ENDIF
      END DO
      END DO

      DO J=J_0, J_1
      DO I=I_0, atmice%IMAXJ(J)
C**** set GTEMP etc. array for ice
        MSI1=si_atm%SNOWI(I,J)+ACE1I

        IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
          MICE(1) = ACE1I-XSI(2)*MSI1
          MICE(2) = XSI(2)*MSI1
c         SNOWL(1)= SNOW
          SNOWL(1)= MSI1-ACE1I
          SNOWL(2)= 0.
        ELSE  ! some snow in second layer
          MICE(1) = 0.
          MICE(2) = ACE1I
          SNOWL(1)= XSI(1)*MSI1
          SNOWL(2)= XSI(2)*MSI1-ACE1I
        ENDIF
        IF(MICE(1).NE.0.) THEN
        atmice%GTEMP(I,J)=Ti2b(si_atm%HSI(1,I,J)/(XSI(1)*MSI1),
     *           1d3*si_atm%SSI(1,I,J)/MICE(1),SNOWL(1),MICE(1))
        ELSE
        atmice%GTEMP(I,J)=Ti(si_atm%HSI(1,I,J)/(XSI(1)*MSI1),0d0)
        ENDIF
        atmice%GTEMP2(I,J)=Ti2b(si_atm%HSI(2,I,J)/(XSI(2)*MSI1),
     *           1d3*si_atm%SSI(2,I,J)/MICE(2),SNOWL(2),MICE(2))

        atmice%GTEMPR(I,J) = atmice%GTEMP(I,J)+TF
        atmice%ZSNOWI(I,J)=si_atm%SNOWI(I,J)/rhos
        si_atm%ZSI(I,J)=(ace1i+si_atm%msi(i,j))/rhoi
#ifdef SCM
        if( SCMopt%Tskin )then
          atmice%GTEMP(I,J) = SCMin%Tskin - TF
          atmice%GTEMP2(I,J) = SCMin%Tskin - TF
          atmice%GTEMPR(I,J) = SCMin%Tskin
        endif
#endif
#ifdef TRACERS_WATER
        atmice%GTRACER(:,I,J) = si_atm%TRSI(:,1,I,J)/
     &       (XSI(1)*MSI1-si_atm%SSI(1,I,J))
#endif
        atmice%FWSIM(I,J) = si_atm%RSI(I,J)*
     &       (MSI1+si_atm%MSI(I,J)-SUM(si_atm%SSI(1:LMI,I,J)))
      END DO
      END DO

      DO J=J_0, J_1
      DO I=I_0, atmice%IMAXJ(J)
#ifndef STANDALONE_OCEAN
        IF(atmice%FOCEAN(I,J).GT.0.) THEN
C**** adjust rad fluxes for change in ice fraction
        if (si_atm%rsi(i,j).gt.si_atm%rsisave(i,j)) then ! ice from ocean
          call RESET_SURF_FLUXES(I,J,1,2,
     &         si_atm%RSISAVE(I,J),si_atm%RSI(I,J))
        elseif (si_atm%rsi(i,j).lt.si_atm%rsisave(i,j)) then ! ocean from ice
          call RESET_SURF_FLUXES(I,J,2,1,
     &         1.-si_atm%RSISAVE(I,J),1.-si_atm%RSI(I,J))
        endif
        ENDIF
#endif
        si_atm%RSISAVE(i,j) = si_atm%RSI(i,j)
      END DO
      END DO

      RETURN
      END SUBROUTINE seaice_to_atmgrid

#ifndef STANDALONE_OCEAN /* remainder of file skipped */
      subroutine daily_seaice(end_of_day,atmocn,atmice)
      use model_com, only : kocean
      use exchange_types, only : atmocn_xchng_vars,atmice_xchng_vars
      implicit none
      logical :: end_of_day
      type(atmocn_xchng_vars) :: atmocn
      type(atmice_xchng_vars) :: atmice

      if(kocean.eq.0) then
        CALL read_seaice(end_of_day,atmocn,atmice)
        CALL reset_gtemp_noice(atmocn,atmice)
      endif

      !if(.not.end_of_day) return

      return
      end subroutine daily_seaice

      SUBROUTINE READ_SEAICE(END_OF_DAY,atmocn,atmice)
!@sum read_seaice reads sea ice concentration (+thickness if provided),
!@+   adjusts ice heat/salt/tracers for conservation purposes, and
!@+   updates diagnostics
!@+   This version is based on the timestream module for netcdf files.
!@ver  beta
!@auth Original Development Team
!@auth M. Kelley restructuring and netcdf-based input options
      USE CONSTANT, only : tf,rhoi
      USE MODEL_COM, only : itime,itimei
      USE MODEL_COM, only : modelEclock
      USE SEAICE_COM, only : RSIstream,ZSIstream
      USE SEAICE, only : xsi,ace1i,ac2oim,ssi0,tfrez,lmi, Ei
#ifdef TRACERS_WATER
      use OldTracer_mod, only: trsi0
      USE SEAICE, only : ntm
#endif
      use seaice_com, only : rsi_exists,zsi_exists,dm
      USE SEAICE, only : z1i,z2oim,fleadoc
      use timestream_mod, only : read_stream
      USE SEAICE_COM, only : si_ocn
      use seaice_com, only : grid=>sigrid
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,atmice_xchng_vars
      IMPLICIT NONE
      LOGICAL :: END_OF_DAY
      type(atmocn_xchng_vars) :: atmocn
      type(atmice_xchng_vars) :: atmice
c
      INTEGER n,I,J,JR
      REAL*8 TFO,MSIrat,RSIMSI,RSIMSIold,dRSI

      REAL*8 ZIMIN,ZIMAX,RSINEW,MSINEW,OPNOCN

      real*8, dimension(2,grid%i_strt_halo:grid%i_stop_halo,
     &                    grid%j_strt_halo:grid%j_stop_halo)
     &     :: TLIM
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo)
     &     :: ZSI


      INTEGER :: IJ_SMFX,IJ_FWIO
      REAL*8, DIMENSION(:,:), POINTER :: RSI,MSI,SNOWI,FWSIM,SSS
      REAL*8, DIMENSION(:,:,:), POINTER :: HSI,SSI,AIJ
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:,:), POINTER :: TRSI,taijn
      INTEGER :: TIJ_ICOCFLX
#endif

      INTEGER :: J_0,J_1, I_0,I_1
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     MSIold,RSIold
      integer :: jyear,jday,itocean,itoice,j_implh,j_implm

      if(.not.rsi_exists) then
        SI_OCN%RSI(:,:) = 0.
        return
      endif

      call modelEclock%get(year=jyear, dayOfYear=jday)

      itocean = atmice%itocean
      itoice = atmice%itoice
      j_implh = atmice%j_implh
      j_implm = atmice%j_implm
      aij => atmice%aij

      FWSIM => ATMICE%FWSIM
      SSS => ATMOCN%SSS
      RSI => SI_OCN%RSI
      MSI => SI_OCN%MSI
      HSI => SI_OCN%HSI
      SSI => SI_OCN%SSI
      SNOWI => SI_OCN%SNOWI
#ifdef TRACERS_WATER
      TRSI => SI_OCN%TRSI
      taijn => atmice%taijn
      TIJ_ICOCFLX = atmice%TIJ_ICOCFLX
#endif
      IJ_SMFX = atmice%IJ_SMFX
      IJ_FWIO = atmice%IJ_FWIO

      I_0 = atmice%I_0
      I_1 = atmice%I_1
      J_0 = atmice%J_0
      J_1 = atmice%J_1

      if(.not.(end_of_day.or.itime.eq.itimei)) return

C**** Get OST, RSI and MSI for current day (save prev. values first)
      MSIold = MSI
      RSIold = RSI

C**** Calculate RSI and MSI for current day

      ZIMIN=Z1I+Z2OIM

      if(zsi_exists) then
        call read_stream(grid,RSIstream,jyear,jday,RSI,TLIM=TLIM)
        call read_stream(grid,ZSIstream,jyear,jday,ZSI,TLIM=TLIM)
        DO J=J_0,J_1
        DO I=I_0,si_ocn%IMAXJ(J)
          if(atmice%focean(i,j).le.0.) cycle
          IF(RSI(I,J).EQ.0.) ZSI(I,J)=ZIMIN
          ZSI(I,J) = max(ZIMIN,ZSI(I,J))
          MSI(I,J) = RHOI*(ZSI(I,J)-Z1I)
        END DO
        END DO
      else
        call read_stream(grid,RSIstream,jyear,jday,RSI)
        DO J=J_0,J_1
        DO I=I_0,si_ocn%IMAXJ(J)
          if(atmice%focean(i,j).le.0.) cycle
          ZIMAX=2d0
          IF(atmice%lat(i,j).GT.0.) ZIMAX=3.5d0 ! northern hemisphere
          rsinew = rsi(i,j)
          MSINEW=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
C**** Ensure that lead fraction is consistent with kocean=1 case
          IF (RSINEW.gt.0) THEN
            OPNOCN=MIN(0.0d0,FLEADOC*RHOI/(RSINEW*(ACE1I+MSINEW)))
            IF (RSINEW.GT.1.-OPNOCN) THEN
              RSINEW = 1.-OPNOCN
              MSINEW=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
            END IF
          END IF
          RSI(I,J)=RSINEW
          MSI(I,J)=MSINEW
        END DO
        END DO
      endif

      where(MSIold.eq.0.) MSIold=MSI  ! does this happen?

C**** Accumulate diagnostics
      IF (end_of_day) THEN
        DO J=J_0,J_1
        DO I=I_0,si_ocn%IMAXJ(J)
          IF (atmice%FOCEAN(I,J).le.0.) cycle
          MSIrat = MSI(I,J)/MSIold(I,J)
          RSIMSI = RSI(I,J)*MSI(I,J)
          RSIMSIold = RSIold(I,J)*MSIold(I,J)
          dRSI = RSI(I,J)-RSIold(I,J)
          AIJ(I,J,IJ_SMFX)=AIJ(I,J,IJ_SMFX)+
     &         (SNOWI(I,J)+ACE1I)*dRSI+RSIMSI-RSIMSIold
          AIJ(I,J,IJ_FWIO)=AIJ(I,J,IJ_FWIO)
     &         -(SNOWI(I,J)+ACE1I-SUM(SSI(1:2,I,J)))*dRSI
     &         -(RSIMSI-RSIMSIold)*(1-SUM(SSI(3:4,I,J))/MSIold(I,J))
          CALL INC_AJ(I,J,ITOICE,J_IMPLM,-atmice%FOCEAN(I,J)*RSIold(I,J)
     &         *(MSI(I,J)-MSIold(I,J))
     &         *(1.-SUM(SSI(3:4,I,J))/MSIold(I,J)))
          CALL INC_AJ(I,J,ITOICE,J_IMPLH,-atmice%FOCEAN(I,J)*RSIold(I,J)
     &         *SUM(HSI(3:4,I,J))*(MSIrat-1.)) 
          CALL INC_AJ(I,J,ITOCEAN,J_IMPLM,-atmice%FOCEAN(I,J)*dRSI
     &         *(MSI(I,J)+ACE1I+SNOWI(I,J)-SUM(SSI(1:2,I,J))
     &         -SUM(SSI(3:4,I,J))*MSIrat))
          CALL INC_AJ(I,J,ITOCEAN,J_IMPLH,-atmice%FOCEAN(I,J)*dRSI
     &         *(SUM(HSI(1:2,I,J))+SUM(HSI(3:4,I,J))*MSIrat))
#ifdef TRACERS_WATER
          DO N=1,NTM
            TAIJN(I,J,TIJ_ICOCFLX,N)=TAIJN(I,J,TIJ_ICOCFLX,N)-
     &           SUM(TRSI(N,1:2,I,J))*(RSI(I,J)       -RSIold(I,J))-
     &           SUM(TRSI(N,3:4,I,J))*(RSI(I,J)*MSIrat-RSIold(I,J))
          END DO
#endif
        END DO
        END DO
      END IF

C****
C**** Update sea ice internal properties (heat,salt,tracers)
C****
      DO J=J_0,J_1
      DO I=I_0,si_ocn%IMAXJ(J)
        IF (atmice%FOCEAN(I,J).le.0.) cycle
        MSIrat = MSI(I,J)/MSIold(I,J)
C**** adjust enthalpy and salt so temperature/salinity remain constant
        HSI(3:4,I,J)=HSI(3:4,I,J)*MSIrat
        SSI(3:4,I,J)=SSI(3:4,I,J)*MSIrat
#ifdef TRACERS_WATER
        TRSI(:,3:4,I,J)=TRSI(:,3:4,I,J)*MSIrat
#endif
C**** adjust some radiative fluxes for changes in ice fraction
c        if (rsi(i,j).gt.rsiold(i,j)) ! ice from ocean
c     *       call RESET_SURF_FLUXES(I,J,1,2,RSIold(I,J),RSI(I,J))
c        if (rsi(i,j).lt.rsiold(i,j)) ! ocean from ice
c     *       call RESET_SURF_FLUXES(I,J,2,1,1.-RSIold(I,J),1.-RSI(I,J))
C****
        TFO=tfrez(sss(i,j))
C**** SET DEFAULTS IF NO OCEAN ICE
        IF (RSI(I,J).LE.0.) THEN
          HSI(1:2  ,I,J)=Ei(TFO,1d3*SSI0)*XSI(1:2)*ACE1I
          HSI(3:LMI,I,J)=Ei(TFO,1d3*SSI0)*XSI(3:LMI)*AC2OIM
          SSI(1:2,I,J)=SSI0*XSI(1:2)*ACE1I
          SSI(3:LMI,I,J)=SSI0*XSI(3:LMI)*AC2OIM
#ifdef TRACERS_WATER
          DO N=1,NTM
            TRSI(N,1:2,I,J)=TRSI0(N)*(1.-SSI0)*XSI(1:2)*ACE1I
            TRSI(N,3:LMI,I,J)=TRSI0(N)*(1.-SSI0)*XSI(3:LMI)*AC2OIM
          END DO
#endif
          SNOWI(I,J)=0.
        END IF
        atmice%FWSIM(I,J)=RSI(I,J)*
     &       (ACE1I+SNOWI(I,J)+MSI(I,J)-SUM(SSI(1:LMI,I,J)))
      END DO
      END DO

C**** REPLICATE VALUES AT POLE
      IF(atmice%HAVE_NORTH_POLE) THEN
        IF (atmice%FOCEAN(1,J_1).gt.0) THEN
          DO I=2,I_1
            RSI(I,J_1)=RSI(1,J_1)
            MSI(I,J_1)=MSI(1,J_1)
            SNOWI(I,J_1)=SNOWI(1,J_1)
            HSI(:,I,J_1)=HSI(:,1,J_1)
            SSI(:,I,J_1)=SSI(:,1,J_1)
#ifdef TRACERS_WATER
            TRSI(:,:,I,J_1)=TRSI(:,:,1,J_1)
#endif
            atmice%FWSIM(I,J_1)=atmice%FWSIM(1,J_1)
          END DO
        END IF
      END IF
      IF(atmice%HAVE_SOUTH_POLE) THEN
        IF (atmice%FOCEAN(1,1).gt.0) THEN
          DO I=2,I_1
            RSI(I,1)=RSI(1,1)
            MSI(I,1)=MSI(1,1)
            SNOWI(I,1)=SNOWI(1,1)
            HSI(:,I,1)=HSI(:,1,1)
            SSI(:,I,1)=SSI(:,1,1)
#ifdef TRACERS_WATER
            TRSI(:,:,I,1)=TRSI(:,:,1,1)
#endif
            atmice%FWSIM(I,1)=atmice%FWSIM(1,1)
          END DO
        END IF
      END IF
      RETURN

      RETURN
      END SUBROUTINE READ_SEAICE

      SUBROUTINE reset_gtemp_noice(atmocn,atmice)
C**** SET DEFAULTS IF NO OCEAN ICE
      USE CONSTANT, only : tf
      USE SEAICE_COM, only : si_ocn
      USE SEAICE, only : tfrez
#ifdef SCM
      USE SCM_COM, only : SCMopt,SCMin
#endif
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,atmice_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atmocn
      type(atmice_xchng_vars) :: atmice
c
      REAL*8 :: TFO
      INTEGER :: I,J,J_0,J_1, I_0,I_1
      I_0 = atmice%I_0
      I_1 = atmice%I_1
      J_0 = atmice%J_0
      J_1 = atmice%J_1
      DO J=J_0,J_1
      DO I=I_0,si_ocn%IMAXJ(J)
        IF (atmice%FOCEAN(I,J).gt.0. .and. si_ocn%RSI(I,J).LE.0.) THEN
          TFO=tfrez(atmocn%sss(i,j))
          atmice%GTEMP(I,J)=TFO
          atmice%GTEMP2(I,J)=TFO
          atmice%GTEMPR(I,J) = TFO+TF
#ifdef SCM
          if( SCMopt%Tskin )then
            atmice%GTEMP(I,J) = SCMin%Tskin - TF
            atmice%GTEMPR(I,J) = SCMin%Tskin
          endif
#endif
        ENDIF
      ENDDO
      ENDDO
      END SUBROUTINE reset_gtemp_noice

      SUBROUTINE ADVSI_DIAG_OCNML(Z1O,Z12O,atmocn,atmice)
!@sum  ADVSI_DIAG_OCNML adjust diagnostics + mlhc for qflux
!@auth Gavin Schmidt
      USE CONSTANT, only : shw,rhows
      USE MODEL_COM, only : kocean,itime
      USE SEAICE, only : ace1i,lmi
      USE SEAICE_COM, only : si_ocn
      use seaice_com, only : grid=>sigrid
#ifdef TRACERS_WATER
      USE SEAICE, only : ntm
#endif
      USE MODEL_COM, only : modelEclock
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,atmice_xchng_vars
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     Z1O,Z12O
      type(atmocn_xchng_vars) :: atmocn
      type(atmice_xchng_vars) :: atmice

      INTEGER I,J,JR,N
      REAL*8 RUN4,ERUN4,TGW,POICE,POCEAN,Z1OMIN,MSINEW
      INTEGER :: J_0,J_1, I_0,I_1

      REAL*8, DIMENSION(:,:), POINTER :: RSI,MSI,SNOWI,FWSIM,FOCEAN
      REAL*8, DIMENSION(:,:,:), POINTER :: HSI,SSI,AIJ
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:,:), POINTER :: TRSI,taijn
      INTEGER :: TIJ_ICOCFLX
#endif
      INTEGER :: IJ_FWIO,J_IMELT,J_HMELT,J_SMELT
      INTEGER :: JMON,itocean,itoice,J_IMPLM,J_IMPLH

      jmon = modelEclock%getMonth()

      itocean = atmice%itocean
      itoice = atmice%itoice
      j_implh = atmice%j_implh
      j_implm = atmice%j_implm
      aij => atmice%aij

      I_0 = atmice%I_0
      I_1 = atmice%I_1
      J_0 = atmice%J_0
      J_1 = atmice%J_1

      FOCEAN => ATMICE%FOCEAN
      FWSIM => ATMICE%FWSIM
      RSI => SI_OCN%RSI
      MSI => SI_OCN%MSI
      HSI => SI_OCN%HSI
      SSI => SI_OCN%SSI
      SNOWI => SI_OCN%SNOWI
#ifdef TRACERS_WATER
      TRSI => SI_OCN%TRSI
      taijn => atmice%taijn
      TIJ_ICOCFLX = atmice%TIJ_ICOCFLX
#endif
      IJ_FWIO = atmice%IJ_FWIO
      J_IMELT = atmice%J_IMELT
      J_HMELT = atmice%J_HMELT
      J_SMELT = atmice%J_SMELT

      DO J=J_0,J_1
      DO I=I_0,si_ocn%IMAXJ(J)
        JR=atmice%JREG(I,J)
        POICE=FOCEAN(I,J)*RSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-RSI(I,J))
        IF (FOCEAN(I,J).gt.0) THEN
          !TGW  = TOCEAN(1,I,J)
          RUN4  = atmice%MSICNV(I,J)
          ERUN4 = 0.  ! TGW*SHW*RUN4 ! force energy conservation
C**** Ensure that we don't run out of ocean if ice gets too thick
          IF (POICE.GT.0) THEN
            Z1OMIN=1.+FWSIM(I,J)/(RHOWS*RSI(I,J))
            IF (Z1OMIN.GT.Z1O(I,J)) THEN
C**** MIXED LAYER DEPTH IS INCREASED TO OCEAN ICE DEPTH + 1 METER
              WRITE(6,602) ITime,I,J,JMON,Z1O(I,J),Z1OMIN,z12o(i,j)
 602          FORMAT (' INCREASE OF MIXED LAYER DEPTH ',I10,3I4,3F10.3)
              Z1O(I,J)=MIN(Z1OMIN, z12o(i,j))
              IF (Z1OMIN.GT.Z12O(I,J)) THEN
C****       ICE DEPTH+1>MAX MIXED LAYER DEPTH :
C****       lose the excess mass to the deep ocean
C**** Calculate freshwater mass to be removed, and then any energy/salt
                MSINEW=MSI(I,J)*(1.-RHOWS*(Z1OMIN-Z12O(I,J))*RSI(I,J)/
     *       (FWSIM(I,J)-RSI(I,J)*(ACE1I+SNOWI(I,J)-SUM(SSI(1:2,I,J))))) 
C**** save diagnostics
                AIJ(I,J,IJ_FWIO)=AIJ(I,J,IJ_FWIO)+RSI(I,J)*(MSI(I,J)
     *               -MSINEW)*(1-SUM(SSI(3:4,I,J))/MSI(I,J))
                CALL INC_AJ(I,J,ITOICE,J_IMELT,-FOCEAN(I,J)*RSI(I,J)
     *               *(MSINEW-MSI(I,J)))
                CALL INC_AJ(I,J,ITOICE,J_HMELT,-FOCEAN(I,J)*RSI(I,J)
     *               *SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)) 
                CALL INC_AJ(I,J,ITOICE,J_SMELT,-FOCEAN(I,J)*RSI(I,J)
     *               *SUM(SSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)) 
                CALL INC_AJ(I,J,ITOICE,J_IMPLM,-FOCEAN(I,J)*RSI(I,J)
     *               *(MSINEW-MSI(I,J))*(1.-SUM(SSI(3:4,I,J))/MSI(I,J)))
                CALL INC_AJ(I,J,ITOICE,J_IMPLH,-FOCEAN(I,J)*RSI(I,J)
     *               *SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)) 
                CALL INC_AREG(I,J,JR,J_IMPLM,-FOCEAN(I,J)*RSI(I,J)
     *               *(MSINEW-MSI(I,J))*(1.-SUM(SSI(3:4,I,J))
     *               /MSI(I,J)))
                CALL INC_AREG(I,J,JR,J_IMPLH,-FOCEAN(I,J)*RSI(I,J)
     *               *SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.))
#ifdef TRACERS_WATER
                DO N=1,NTM
                  TAIJN(I,J,TIJ_ICOCFLX,N)=TAIJN(I,J,TIJ_ICOCFLX,N)-
     *              SUM(TRSI(N,3:4,I,J))*RSI(I,J)*(MSINEW/MSI(I,J)-1.)
                END DO
#endif
C**** update heat and salt
                HSI(3:4,I,J) = HSI(3:4,I,J)*(MSINEW/MSI(I,J))
                SSI(3:4,I,J) = SSI(3:4,I,J)*(MSINEW/MSI(I,J))
#ifdef TRACERS_WATER
                TRSI(:,3:4,I,J) = TRSI(:,3:4,I,J)*(MSINEW/MSI(I,J))
#endif
                MSI(I,J)=MSINEW
                FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)
     *               -SUM(SSI(1:LMI,I,J)))
              END IF
            END IF
          END IF
          atmocn%MLHC(I,J) = SHW*(Z1O(I,J)*RHOWS-FWSIM(I,J))
C**** Open Ocean diagnostics
          CALL INC_AJ(I,J,ITOCEAN,J_IMPLM, RUN4*POCEAN)
          CALL INC_AJ(I,J,ITOCEAN,J_IMPLH,ERUN4*POCEAN)
C**** Ice-covered ocean diagnostics
          CALL INC_AJ(I,J,ITOICE,J_IMPLM, RUN4*POICE)
          CALL INC_AJ(I,J,ITOICE,J_IMPLH,ERUN4*POICE)
C**** regional diagnostics
          CALL INC_AREG(I,J,JR,J_IMPLM, RUN4*FOCEAN(I,J))
          CALL INC_AREG(I,J,JR,J_IMPLH,ERUN4*FOCEAN(I,J))
        END IF
      END DO
      END DO

      RETURN
      END SUBROUTINE ADVSI_DIAG_OCNML

#endif /* ifndef STANDALONE_OCEAN */
