      MODULE AERO_DISCRETE
      USE AERO_PARAM, ONLY: DENSP, PI6
      USE AERO_COAG,  ONLY: TOTAL_COAG_COEF, BROWNIAN_COAG_COEF
      IMPLICIT NONE
      INTEGER, PARAMETER :: NBINS = 200
      INTEGER, PARAMETER :: NSPCS =   2
      REAL(8) :: APDF(NBINS,1+NSPCS)       ! mode A: number conc. [#/m^3] and mass concs. [ug/m^3] for each bin. 
      REAL(8) :: BPDF(NBINS,1+NSPCS)       ! mode B: number conc. [#/m^3] and mass concs. [ug/m^3] for each bin. 
      REAL(8) :: CPDF(NBINS,1+NSPCS)       ! mode C: number conc. [#/m^3] and mass concs. [ug/m^3] for each bin. 
      REAL(8) :: DGRID(NBINS)              ! fixed diameter        grid [um]
      REAL(8) :: VGRID(NBINS)              ! fixed volume/particle grid [um^3/particle]
      REAL(8) :: MGRID(NBINS)              ! fixed mass/particle   grid [ug/particle]
      REAL(8) :: D3GRID(NBINS)             ! fixed diameter-cubed  grid [um^3/particle]
      REAL(8) :: DLOWER(NBINS)             ! lower boundary fixed diameter grid [um]
      REAL(8) :: DUPPER(NBINS)             ! upper boundary fixed diameter grid [um]
      REAL(8), PARAMETER :: DMIN     =  0.001D+00       ! smallest particle diameter of the discrete grid [um]
      REAL(8), PARAMETER :: DMAX     = 20.000D+00       ! largest  particle diameter of the discrete grid [um]
      REAL(8), SAVE      :: RDMIN    =  0.000D+00       ! reciprocal of DMIN to optimize coagulation [1/um]
      REAL(8), SAVE      :: RDLOGDSC =  0.000D+00       ! reciprocal of log10 of the grid spacing [1]
      REAL(8), PARAMETER :: XKB      =  1.3806505D-23   ! [J/K] http://en.wikipedia.org/wiki/Boltzmann_constant
      !-----------------------------------------------------------------------------------------------------------------
      ! ISPCA, ISPCB, ISPCC: 1=SULF, 2=BCAR, 3=OCAR, 4=DUST, 5=SEAS
      !-----------------------------------------------------------------------------------------------------------------
      INTEGER, SAVE :: ISPCA               ! index of the chemical species for mode A: set in AERO_INIT
      INTEGER, SAVE :: ISPCB               ! index of the chemical species for mode B: set in AERO_INIT
      INTEGER, SAVE :: ISPCC               ! index of the chemical species for mode C: set in AERO_INIT
      INTEGER, PARAMETER :: ITCOAG = 10    ! number of subdivisions of the time step for integration of coagulation 
      REAL(8) :: DTCOAG = 0.0D+00
      REAL(8) :: KIJ_DISCRETE(NBINS,NBINS) ! fixed coagulation coefficients for the discrete grid [m^3/s/particle]
      REAL(8) :: FRAC_LO (NBINS,NBINS)     ! used in DISCRETE_INTRACOAG [1]
      REAL(8) :: FRAC_HI (NBINS,NBINS)     ! used in DISCRETE_INTRACOAG [1]
      INTEGER :: NLO_GRID(NBINS,NBINS)     ! used in DISCRETE_INTRACOAG [1]
      INTEGER :: NHI_GRID(NBINS,NBINS)     ! used in DISCRETE_INTRACOAG [1]
      

      CONTAINS


      SUBROUTINE DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
!-----------------------------------------------------------------------------------------------------------------------
!@auth    Susanne Bauer/Doug Wright
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.
 
      INTEGER, INTENT(   IN) :: ICSET                      ! identifies test case [1]
      REAL(8), INTENT(   IN) :: NA,      NB,      NC       ! number concentrations for modes A, B, and C [#/m^3]
      REAL(8), INTENT(   IN) :: DGA,     DGB,     DGC      ! geo. mean diameters   for modes A, B, and C [um]
      REAL(8), INTENT(   IN) :: SIGMAGA, SIGMAGB, SIGMAGC  ! geo. std. deviations  for modes A, B, and C [1]
      REAL(8), INTENT(   IN) :: TEMP                       ! ambient temperature [K]
      REAL(8), INTENT(   IN) :: PRES                       ! ambient pressure    [Pa]
      REAL(8), INTENT(  OUT) :: MASSA,   MASSB,   MASSC    ! mass concentrations   for modes A, B, and C [ug/m^3]

      ! Local variables.

      INTEGER :: I, J, NHI, NLO
      REAL(8) :: SCALE, FA, FB, FC, WA(NBINS), WB(NBINS), WC(NBINS), WTOTA, WTOTB, WTOTC, DNEWCUB
      REAL(8) :: DMINL, DMAXL, ETAA, VP_PDF, MFRAC(NBINS,NSPCS)

      ! For the call to the coagulation coefficient routines

      REAL(4) :: DI, DJ    ! ambient particle diameters [um]
      REAL(4) :: BETAIJ    ! coagulation coefficients [m^2/s/particle]
      REAL(4) :: TEMPL     ! ambient temperature [K]
      REAL(4) :: PRESL     ! ambient pressure [Pa]

      REAL(8), PARAMETER :: ONETHIRD = 1.0D+00 / 3.0D+00 


      APDF(:,:) = 0.0D+00
      BPDF(:,:) = 0.0D+00
      CPDF(:,:) = 0.0D+00
      
      ! Set up discrete grid.

      SELECT CASE( ICSET )
      CASE ( 10 )
        DMAXL = 0.600D+00
        DMINL = 0.006D+00
      CASE ( 11 )
        DMAXL = 1.000D+00
        DMINL = 0.001D+00
      CASE ( 12 )
        DMAXL = DMAX
        DMINL = DMIN
      CASE ( 13 )
        DMAXL = DMAX
        DMINL = DMIN
      CASE ( 14 )
        DMAXL = 100.000D+00
        DMINL =   0.001D+00
      CASE ( 15 )
        DMAXL = DMAX
        DMINL = DMIN
      CASE ( 16 )
        DMAXL = DMAX
        DMINL = DMIN
      CASE ( 17 )
        DMAXL = 200.000D+00
        DMINL =   0.001D+00
      CASE ( 18 )
        DMAXL = 8.000D+00
        DMINL = 0.001D+00
      CASE ( 19 )
        DMAXL = DMAX
        DMINL = DMIN
      CASE ( 20 )
        DMAXL = DMAX
        DMINL = DMIN
      CASE DEFAULT
        WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete_init (A).'
        STOP
      END SELECT
      SCALE    = ( DMAXL / DMINL )**(1.0D+00/REAL(NBINS-1))
      RDLOGDSC = 1.0D+00 / LOG10( SCALE )
      RDMIN    = 1.0D+00 / DMINL
      ! WRITE(*,'(I4,5D14.5)') ICSET,DMINL,DMAXL,SCALE,RDLOGDSC,RDMIN
      WRITE(36,'(/A,I6,/)') 'ICSET = ', ICSET
      WRITE(36,*)'I, DLOWER(I) [um], DGRID(I) [um], DUPPER(I) [um], VGRID(I) [um^3/particle], MGRID(I) [ug/particle]'
      DO I=1, NBINS
        DGRID(I)  = DMINL * SCALE**(I-1)                  ! [um]
        DLOWER(I) = DGRID(I) / SCALE**0.5D+00             ! [um]
        DUPPER(I) = DGRID(I) * SCALE**0.5D+00             ! [um]
        D3GRID(I) = DGRID(I)**3                           ! [um^3/particle]
        VGRID(I)  = PI6 * DGRID(I)**3                     ! [um^3/particle]
        MGRID(I)  = 1.0D-06 * DENSP * VGRID(I)            ! [ug/particle]
        WRITE(36,'(I4,3F12.7,2D16.7)') I, DLOWER(I), DGRID(I), DUPPER(I), VGRID(I), MGRID(I)
      ENDDO

      WTOTA = 0.0D+00
      WTOTB = 0.0D+00
      WTOTC = 0.0D+00
      FA = 0.0D+00
      FB = 0.0D+00
      FC = 0.0D+00
      WRITE(36,'(/A,I6,/)') 'ICSET = ', ICSET
      WRITE(36,*)'I, DGRID(I) [um], FLN [/m^3/um], W(I) [#/m^3], WTOT [#/m^3]'
      WRITE(36,*)'I, DGRID(I), FA, FB, FC, WA(I), WB(I), WC(I), WTOTA, WTOTB, WTOTC'
      SELECT CASE( ICSET )
      CASE ( 10 )
        WA(:) = 0.0D+00
        WB(:) = 0.0D+00
        WC(:) = 0.0D+00
        WA(1) = NA
        WTOTA = SUM( WA(:) ) 
        WTOTB = SUM( WB(:) ) 
        WTOTC = SUM( WC(:) )
        DO I=1, NBINS 
          WRITE(36,'(I4,F10.5,9D14.6)') I, DGRID(I), 0D0, 0D0, 0D0, WA(I), WB(I), WC(I), WTOTA, WTOTB, WTOTC
        ENDDO
      CASE ( 11 )
        VP_PDF = 0.0005236D+00
        WB(:) = 0.0D+00
        WC(:) = 0.0D+00
        WTOTB = 0.0D+00
        WTOTC = 0.0D+00
        DO I=1, NBINS
          WA(I) = NA * ( (PI6*(DUPPER(I)**3-DLOWER(I)**3))/VP_PDF ) * EXP( -VGRID(I)/VP_PDF ) 
          WTOTA = WTOTA + WA(I) 
          WRITE(36,'(I4,F10.5,9D14.6)') I, DGRID(I), 0D0, 0D0, 0D0, WA(I), WB(I), WC(I), WTOTA, WTOTB, WTOTC
        ENDDO
      CASE ( 12, 13, 14, 15, 16, 17, 18, 19, 20 )
        DO I=1, NBINS
          FA = NA * FLN( DGRID(I), DGA, SIGMAGA )
          FB = NB * FLN( DGRID(I), DGB, SIGMAGB )
          FC = NC * FLN( DGRID(I), DGC, SIGMAGC )
          WA(I) = FA * ( DUPPER(I) - DLOWER(I) )
          WB(I) = FB * ( DUPPER(I) - DLOWER(I) )
          WC(I) = FC * ( DUPPER(I) - DLOWER(I) )
          WTOTA = WTOTA + WA(I) 
          WTOTB = WTOTB + WB(I) 
          WTOTC = WTOTC + WC(I) 
          WRITE(36,'(I4,F10.5,9D14.6)') I, DGRID(I), FA, FB, FC, WA(I), WB(I), WC(I), WTOTA, WTOTB, WTOTC
        ENDDO
      CASE DEFAULT
        WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete_init (B).'
        STOP
      END SELECT

      WA(:) = WA(:) * ( NA / MAX( WTOTA, 1.0D-30) )  ! Renormalize to the precise input number concentration.
      WB(:) = WB(:) * ( NB / MAX( WTOTB, 1.0D-30) )  ! Renormalize to the precise input number concentration.
      WC(:) = WC(:) * ( NC / MAX( WTOTC, 1.0D-30) )  ! Renormalize to the precise input number concentration.

      WRITE(36,*)'Check sums on WA(:), WB(:), WC(:)'
      WRITE(36,'(3D22.14)') NA,        NB,        NC
      WRITE(36,'(3D22.14)') SUM(WA(:)),SUM(WB(:)),SUM(WC(:))

      ! Get total mass concentration summed over all bins. 

      MASSA = SUM( WA(:) * MGRID(:) ) 
      MASSB = SUM( WB(:) * MGRID(:) ) 
      MASSC = SUM( WC(:) * MGRID(:) ) 

      ! Initialize the main discrete grids for number and mass concentrations.
      ! WRITE(*,*)ISPCA,ISPCB,ISPCC

      APDF(:,1) = WA(:)                                   ! [#/m^3]
      BPDF(:,1) = WB(:)                                   ! [#/m^3]
      CPDF(:,1) = WC(:)                                   ! [#/m^3]
      MFRAC(:,:) = 0.0D+00                                ! [1]
      SELECT CASE( ICSET )
      CASE( 10, 11 )                                  ! one-component for all modes 
        DO I=1, NBINS
          MFRAC(I,2) = REAL(I) / REAL(NBINS)              ! [1]
          MFRAC(I,1) = 1.0D+00 - MFRAC(I,2)               ! [1]
        ENDDO      
        DO I=1, NBINS
          DO J=1, ISPCA
            APDF(I,1+J) = WA(I) * MGRID(I) * MFRAC(I,J)   ! [ug/m^3]
            BPDF(I,1+J) = WB(I) * MGRID(I) * MFRAC(I,J)   ! [ug/m^3]
            CPDF(I,1+J) = WC(I) * MGRID(I) * MFRAC(I,J)   ! [ug/m^3]
          ENDDO
        ENDDO
      CASE( 12 )
        APDF(:,2:NSPCS+1) = 0.0D+00                       ! [ug/m^3]
        BPDF(:,2:NSPCS+1) = 0.0D+00                       ! [ug/m^3]
        CPDF(:,2:NSPCS+1) = 0.0D+00                       ! [ug/m^3]
        APDF(:,2) = WA(:) * MGRID(:)                      ! [ug/m^3] pure sulfate
        BPDF(:,2) = WB(:) * MGRID(:)                      ! [ug/m^3] pure BC
      CASE( 13, 14, 15, 16, 17, 18, 19, 20 )
        APDF(:,2:NSPCS+1) = 0.0D+00                       ! [ug/m^3]
        BPDF(:,2:NSPCS+1) = 0.0D+00                       ! [ug/m^3]
        CPDF(:,2:NSPCS+1) = 0.0D+00                       ! [ug/m^3]
        APDF(:,2) = WA(:) * MGRID(:)                      ! [ug/m^3] pure sulfate
        BPDF(:,3) = WB(:) * MGRID(:)                      ! [ug/m^3] pure BC
      CASE DEFAULT
        WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete_init (C).'
        STOP
      END SELECT

      ! Calculate particle partitioning between adjacient grid points to optimize the coagulation integration.

      FRAC_LO(:,:) = 0.0D+00
      FRAC_HI(:,:) = 0.0D+00
      NLO_GRID(:,:) = 1
      NHI_GRID(:,:) = 1
      DO I=1, NBINS
      DO J=1, NBINS
        DNEWCUB = D3GRID(I) + D3GRID(J)
        NLO  = INT( LOG10( RDMIN * DNEWCUB**ONETHIRD ) * RDLOGDSC ) + 1
        NLO  = MIN( NBINS-1, NLO )
        NHI  = NLO + 1
        NHI  = MIN( NBINS, NHI )
        NLO_GRID(I,J) = NLO
        NHI_GRID(I,J) = NHI
        IF(     1.000000001D+00*DNEWCUB .LT. D3GRID(NLO) ) THEN
          FRAC_LO(I,J) = 1.0D+00
          FRAC_HI(I,J) = 1.0D+00 - FRAC_LO(I,J)
          IF( NLO .EQ. 1 .OR. NLO .EQ. NBINS-1 ) GOTO 101
          IF( NHI .EQ. 2 .OR. NHI .EQ. NBINS   ) GOTO 101
          WRITE(*,90) NLO, NHI, D3GRID(NLO), DNEWCUB, D3GRID(NHI)
          STOP
        ELSEIF( 0.999999999D+00*DNEWCUB .GT. D3GRID(NHI) ) THEN
          FRAC_LO(I,J) = 0.0D+00
          FRAC_HI(I,J) = 1.0D+00 - FRAC_LO(I,J)
          IF( NLO .EQ. 1 .OR. NLO .EQ. NBINS-1 ) GOTO 101
          IF( NHI .EQ. 2 .OR. NHI .EQ. NBINS   ) GOTO 101
          WRITE(*,90) NLO, NHI, D3GRID(NLO), DNEWCUB, D3GRID(NHI)
          STOP
        ELSE 
          FRAC_LO(I,J) = ( D3GRID(NHI) - DNEWCUB ) / ( D3GRID(NHI) - D3GRID(NLO) )
          FRAC_HI(I,J) = 1.0D+00 - FRAC_LO(I,J)
        ENDIF
101     CONTINUE
        ! WRITE(36,'(4I5,2D15.5)') I,J,NLO,NHI,FRAC_LO(I,J), FRAC_HI(I,J)
      ENDDO
      ENDDO

      ! Calculate the fixed grid coagulation coefficients.

      SELECT CASE( ICSET )
      CASE ( 10, 11 )    
        !---------------------------------------------------------------------------------------------------------------
        ! ETAA: Dynamic viscosity of air [kg/m/s], Jacobson, 1999, eq.(4.55)
        !---------------------------------------------------------------------------------------------------------------
        ETAA = 1.832D-05*(416.16D+00/(TEMP+120.0D+00))*(TEMP/296.16D+00)**1.5D+00
        KIJ_DISCRETE(:,:) = 8.0D+00 * XKB * TEMP / ( 3.0D+00 * ETAA )
      CASE( 12, 13, 14, 15, 16, 17, 18, 19, 20 )
        KIJ_DISCRETE(:,:) = 0.0D+00
        TEMPL = REAL( TEMP )                                         ! convert to single precision               
        PRESL = REAL( PRES )                                         ! convert to single precision
        DO I=1, NBINS
        DO J=I, NBINS
          DI = REAL( DGRID(I) )                                      ! convert to single precision
          DJ = REAL( DGRID(J) )                                      ! convert to single precision
!         CALL BROWNIAN_COAG_COEF( DI, DJ, TEMPL, PRESL, BETAIJ )    ! all variables single precision
          CALL TOTAL_COAG_COEF   ( DI, DJ, TEMPL, PRESL, BETAIJ )    ! all variables single precision
          KIJ_DISCRETE(I,J) = REAL( BETAIJ )                         ! convert to double precision
          KIJ_DISCRETE(J,I) = REAL( BETAIJ )                         ! convert to double precision
        ENDDO
        ENDDO
      CASE DEFAULT
        WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete_init (D).'
        STOP
        !---------------------------------------------------------------------------------------------------------------
        ! ETAA = 1.832D-05*(416.16D+00/(TEMP+120.0D+00))*(TEMP/296.16D+00)**1.5D+00
        ! KIJ_DISCRETE(:,:) = 8.0D+00 * XKB * TEMP / ( 3.0D+00 * ETAA )
        ! KIJ_DISCRETE(:,:) = 1.0D-14
        !---------------------------------------------------------------------------------------------------------------
      END SELECT

      WRITE(37,'(F8.2,11D13.5)') 0D0, 
     &                          1.0D-06*SUM(APDF(:,1)), 1.0D-06*SUM(BPDF(:,1)), 1.0D-06*SUM(CPDF(:,1)),
     &                          SUM( APDF(:,2) ),  SUM( BPDF(:,2) ),  SUM( CPDF(:,2) ),
     &                          SUM( APDF(:,3) ),  SUM( BPDF(:,3) ),  SUM( CPDF(:,3) ),
     &                          SUM( APDF(:,2) ) + SUM( BPDF(:,2) ) + SUM( CPDF(:,2) ),
     &                          SUM( APDF(:,3) ) + SUM( BPDF(:,3) ) + SUM( CPDF(:,3) )

      CALL DISCRETE_OUT(38,ICSET,0D0)      
      
90    FORMAT('DISCRETE_INIT: Bad NLO or NHI: NLO, NHI, D3GRID(NLO), DNEWCUB, D3GRID(NHI) = ',/,2I10,3D22.14)
      RETURN
      END SUBROUTINE DISCRETE_INIT


      SUBROUTINE DISCRETE( ICSET, TIMEH, TSTEP )
!-----------------------------------------------------------------------------------------------------------------------
!     10-30-06, DLW: Discrete model of the pdf used to obtain results for evaluation of MATRIX. 
!-----------------------------------------------------------------------------------------------------------------------

      ! Arguments.
 
      INTEGER             :: ICSET                 ! label for set of initial conditions [1] 
      REAL(8), INTENT(IN) :: TIMEH                 ! model time [h]
      REAL(8), INTENT(IN) :: TSTEP                 ! model physics time step [s]

      ! Local variables.

      INTEGER :: IT
      LOGICAL, SAVE      :: FIRSTIME           = .TRUE.
      LOGICAL, PARAMETER :: DISCRETE_COAG_FLAG = .TRUE.

!----------------------------------------------------------------------------------------------------------------------
!     Begin execution.
!----------------------------------------------------------------------------------------------------------------------
      IF ( FIRSTIME ) THEN
        FIRSTIME = .FALSE.
        DTCOAG = TSTEP / REAL( ITCOAG )
      ENDIF                           


      IF( DISCRETE_COAG_FLAG ) THEN
        SELECT CASE( ICSET )
        CASE ( 10, 11, 18 )
          DO IT=1, ITCOAG
            CALL DISCRETE_INTRACOAG    ( APDF )
          ENDDO
        CASE ( 12, 20 )
          DO IT=1, ITCOAG
            CALL DISCRETE_INTRACOAG    ( APDF )
            CALL DISCRETE_INTRACOAG    ( BPDF )
            CALL DISCRETE_INTERCOAG_AB ( APDF, BPDF )
          ENDDO
        CASE ( 13, 14, 15, 16, 17, 19 )
          DO IT=1, ITCOAG
            CALL DISCRETE_INTRACOAG    ( APDF )
            CALL DISCRETE_INTRACOAG    ( BPDF )
            CALL DISCRETE_INTRACOAG    ( CPDF )
            CALL DISCRETE_INTERCOAG_AB ( APDF, CPDF )
            CALL DISCRETE_INTERCOAG_AB ( BPDF, CPDF )
            CALL DISCRETE_INTERCOAG_ABC( APDF, BPDF, CPDF )
          ENDDO
        CASE DEFAULT
          WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete.'
          STOP
        END SELECT
      ENDIF

      WRITE(37,'(F8.2,11D13.5)') TIMEH+TSTEP/3.6D+03,
     &                          1.0D-06*SUM(APDF(:,1)), 1.0D-06*SUM(BPDF(:,1)), 1.0D-06*SUM(CPDF(:,1)),
     &                          SUM( APDF(:,2) ),  SUM( BPDF(:,2) ),  SUM( CPDF(:,2) ),
     &                          SUM( APDF(:,3) ),  SUM( BPDF(:,3) ),  SUM( CPDF(:,3) ),
     &                          SUM( APDF(:,2) ) + SUM( BPDF(:,2) ) + SUM( CPDF(:,2) ),
     &                          SUM( APDF(:,3) ) + SUM( BPDF(:,3) ) + SUM( CPDF(:,3) )

      RETURN
      END SUBROUTINE DISCRETE


      SUBROUTINE DISCRETE_INTRACOAG( PDF )
!-----------------------------------------------------------------------------------------------------------------------
!     10-30-06, DLW: Discrete model of the pdf used to obtain results for evaluation of MATRIX. 
!-----------------------------------------------------------------------------------------------------------------------

      ! Arguments.
 
      REAL(8), INTENT(INOUT) :: PDF(NBINS,1+NSPCS)    ! working array for discrete variables [#/m^3], [ug/m^3]

      ! Local variables.

      INTEGER :: I, J, NLO, NHI 
      REAL(8) :: W(NBINS), DELW(NBINS), DW, DELM(NBINS,NSPCS), MP(NBINS,NSPCS), DM(NSPCS)

      W(:)      = PDF(:,1)   ! load number concentrations in work array [#/m^3]
      DELW(:)   = 0.0D+00
      DELM(:,:) = 0.0D+00
      DO I=1, NBINS
        MP(I,:) = PDF(I,2:NSPCS+1) / MAX( W(I), 1.0D-30 )   ! load mass per particle for each bin and species [ug]
      ENDDO
      DO I=1, NBINS
      DO J=I, NBINS
        IF(I.EQ.J) THEN
          DW        = 0.5D+00*KIJ_DISCRETE(I,J)*W(I)*W(J)*DTCOAG
          DM(:)     = 2.0D+00 * DW * MP(I,:)
          DELW(I)   = DELW(I)   - 2.0D+00*DW
          DELM(I,:) = DELM(I,:) - DM(:)
        ELSE
          DW = KIJ_DISCRETE(I,J)*W(I)*W(J)*DTCOAG
          DM(:)     = DW * ( MP(I,:) + MP(J,:) )
          DELW(I)   = DELW(I)   - DW
          DELW(J)   = DELW(J)   - DW
          DELM(I,:) = DELM(I,:) - DW*MP(I,:)
          DELM(J,:) = DELM(J,:) - DW*MP(J,:)
        ENDIF
        ! WRITE(*,'(2I4,7D12.3)')I,J,DW,WLO,WHI,W(I),W(J),FRAC_LO(I,J),FRAC_HI(I,J)
        NLO         = NLO_GRID(I,J)
        NHI         = NHI_GRID(I,J)
        DELW(NLO)   = DELW(NLO)   + DW    * FRAC_LO(I,J)
        DELW(NHI)   = DELW(NHI)   + DW    * FRAC_HI(I,J)
        DELM(NLO,:) = DELM(NLO,:) + DM(:) * FRAC_LO(I,J)
        DELM(NHI,:) = DELM(NHI,:) + DM(:) * FRAC_HI(I,J)
      ENDDO
      ENDDO
      PDF(:,1)         = W(:)             + DELW(:)          ! update output array [#/m^3]
      PDF(:,2:NSPCS+1) = PDF(:,2:NSPCS+1) + DELM(:,:)        ! update output array [ug/m^3]
      DO I=1, NSPCS+1
        PDF(:,I) = MAX( PDF(:,I), 0.0D-30 )
      ENDDO


      RETURN
      END SUBROUTINE DISCRETE_INTRACOAG


      SUBROUTINE DISCRETE_INTERCOAG_AB( PDFA, PDFB )
!-----------------------------------------------------------------------------------------------------------------------
!     10-30-06, DLW: Discrete model of the pdf used to obtain results for evaluation of MATRIX. 
!
!     PDFA coagulates with PDFB to produce additional particles in PDFB. 
!-----------------------------------------------------------------------------------------------------------------------

      ! Arguments.

      REAL(8), INTENT(INOUT) :: PDFA(NBINS,1+NSPCS)    ! working array for discrete variables [#/m^3], [ug/m^3]
      REAL(8), INTENT(INOUT) :: PDFB(NBINS,1+NSPCS)    ! working array for discrete variables [#/m^3], [ug/m^3]

      ! Local variables.

      INTEGER :: I, J, NLO, NHI 
      REAL(8) :: DW, DM(NSPCS)
      REAL(8) :: WA(NBINS), DELWA(NBINS), DELMA(NBINS,NSPCS), MPA(NBINS,NSPCS), DMA(NSPCS)
      REAL(8) :: WB(NBINS), DELWB(NBINS), DELMB(NBINS,NSPCS), MPB(NBINS,NSPCS), DMB(NSPCS)

      WA(:)      = PDFA(:,1)             ! load number concentrations in work array [#/m^3]
      WB(:)      = PDFB(:,1)             ! load number concentrations in work array [#/m^3]
      DELWA(:)   = 0.0D+00
      DELWB(:)   = 0.0D+00
      DELMA(:,:) = 0.0D+00
      DELMB(:,:) = 0.0D+00
      DO I=1, NBINS
        MPA(I,:) = PDFA(I,2:NSPCS+1) / MAX( WA(I), 1.0D-30 )   ! load mass per particle for each bin and species [ug]
        MPB(I,:) = PDFB(I,2:NSPCS+1) / MAX( WB(I), 1.0D-30 )   ! load mass per particle for each bin and species [ug]
      ENDDO
      DO I=1, NBINS
      DO J=1, NBINS
        DW           = KIJ_DISCRETE(I,J)*WA(I)*WB(J)*DTCOAG
        DMA(:)       = DW * MPA(I,:)
        DMB(:)       = DW * MPB(J,:)
        DM(:)        = DMA(:) + DMB(:)
        DELWA(I)     = DELWA(I)   - DW
        DELWB(J)     = DELWB(J)   - DW
        DELMA(I,:)   = DELMA(I,:) - DMA(:)
        DELMB(J,:)   = DELMB(J,:) - DMB(:)
        NLO          = NLO_GRID(I,J)
        NHI          = NHI_GRID(I,J)
        DELWB(NLO)   = DELWB(NLO)   + DW    * FRAC_LO(I,J)
        DELWB(NHI)   = DELWB(NHI)   + DW    * FRAC_HI(I,J)
        DELMB(NLO,:) = DELMB(NLO,:) + DM(:) * FRAC_LO(I,J)
        DELMB(NHI,:) = DELMB(NHI,:) + DM(:) * FRAC_HI(I,J)
      ENDDO
      ENDDO
      PDFA(:,1)         = WA(:)             + DELWA(:)          ! update output array [#/m^3]
      PDFB(:,1)         = WB(:)             + DELWB(:)          ! update output array [#/m^3]
      PDFA(:,2:NSPCS+1) = PDFA(:,2:NSPCS+1) + DELMA(:,:)        ! update output array [ug/m^3]
      PDFB(:,2:NSPCS+1) = PDFB(:,2:NSPCS+1) + DELMB(:,:)        ! update output array [ug/m^3]

      DO I=1, NSPCS+1
        PDFA(:,I) = MAX( PDFA(:,I), 0.0D-30 )
        PDFB(:,I) = MAX( PDFB(:,I), 0.0D-30 )
      ENDDO

      RETURN
      END SUBROUTINE DISCRETE_INTERCOAG_AB


      SUBROUTINE DISCRETE_INTERCOAG_ABC( PDFA, PDFB, PDFC )
!-----------------------------------------------------------------------------------------------------------------------
!     10-30-06, DLW: Discrete model of the pdf used to obtain results for evaluation of MATRIX. 
!
!     PDFA coagulates with PDFB to produce additional particles in PDFC. 
!     Either PDFA or PDFB may be identical with PDFC, but PDFA and PDFB cannot be identical.
!
!     IF PDFA is not PDFC, and PDFB is not PDFC --> ICASE = 0
!     IF PDFA is. PDFC                          --> ICASE = 1
!     IF PDFB is. PDFC                          --> ICASE = 2
!-----------------------------------------------------------------------------------------------------------------------

      ! Arguments.

      REAL(8), INTENT(INOUT) :: PDFA(NBINS,1+NSPCS)    ! working array for discrete variables [#/m^3], [ug/m^3]
      REAL(8), INTENT(INOUT) :: PDFB(NBINS,1+NSPCS)    ! working array for discrete variables [#/m^3], [ug/m^3]
      REAL(8), INTENT(INOUT) :: PDFC(NBINS,1+NSPCS)    ! working array for discrete variables [#/m^3], [ug/m^3]

      ! Local variables.

      INTEGER :: I, J, NLO, NHI 
      REAL(8) :: DW, DM(NSPCS)
      REAL(8) :: WA(NBINS), DELWA(NBINS), DELMA(NBINS,NSPCS), MPA(NBINS,NSPCS), DMA(NSPCS)
      REAL(8) :: WB(NBINS), DELWB(NBINS), DELMB(NBINS,NSPCS), MPB(NBINS,NSPCS), DMB(NSPCS)
      REAL(8) ::            DELWC(NBINS), DELMC(NBINS,NSPCS)

      WA(:)      = PDFA(:,1)             ! load number concentrations in work array [#/m^3]
      WB(:)      = PDFB(:,1)             ! load number concentrations in work array [#/m^3]
      DELWA(:)   = 0.0D+00
      DELWB(:)   = 0.0D+00
      DELWC(:)   = 0.0D+00
      DELMA(:,:) = 0.0D+00
      DELMB(:,:) = 0.0D+00
      DELMC(:,:) = 0.0D+00
      DO I=1, NBINS
        MPA(I,:) = PDFA(I,2:NSPCS+1) / MAX( WA(I), 1.0D-30 )   ! load mass per particle for each bin and species [ug]
        MPB(I,:) = PDFB(I,2:NSPCS+1) / MAX( WB(I), 1.0D-30 )   ! load mass per particle for each bin and species [ug]
        ! WRITE(39,'(I6,8D15.5)') I, MPB(I,:), PDFB(I,2:NSPCS+1), WB(I)
      ENDDO
      ! WRITE(39,*)'  '
      DO I=1, NBINS
      DO J=1, NBINS
        DW           = KIJ_DISCRETE(I,J)*WA(I)*WB(J)*DTCOAG
!       IF( DW .GT. MIN( WA(I), WB(J) ) ) THEN
!         WRITE(*,*)'DW .GT. MIN( WA(I), WB(J) ): I,J, DW, WA, WB = ', I, J, DW, WA(I), WB(J)
!         STOP
!       ENDIF
        DMA(:)       = DW * MPA(I,:)
        DMB(:)       = DW * MPB(J,:)
        DM(:)        = DMA(:) + DMB(:)
        DELWA(I)     = DELWA(I)   - DW
        DELWB(J)     = DELWB(J)   - DW
        DELMA(I,:)   = DELMA(I,:) - DMA(:)
        DELMB(J,:)   = DELMB(J,:) - DMB(:)
        NLO          = NLO_GRID(I,J)
        NHI          = NHI_GRID(I,J)
        DELWC(NLO)   = DELWC(NLO)   + DW    * FRAC_LO(I,J)
        DELWC(NHI)   = DELWC(NHI)   + DW    * FRAC_HI(I,J)
        DELMC(NLO,:) = DELMC(NLO,:) + DM(:) * FRAC_LO(I,J)
        DELMC(NHI,:) = DELMC(NHI,:) + DM(:) * FRAC_HI(I,J)
      ENDDO
      ENDDO
      PDFA(:,1)         = WA(:)             + DELWA(:)          ! update output array [#/m^3]
      PDFB(:,1)         = WB(:)             + DELWB(:)          ! update output array [#/m^3]
      PDFC(:,1)         = PDFC(:,1)         + DELWC(:)          ! update output array [#/m^3]
      PDFA(:,2:NSPCS+1) = PDFA(:,2:NSPCS+1) + DELMA(:,:)        ! update output array [ug/m^3]
      PDFB(:,2:NSPCS+1) = PDFB(:,2:NSPCS+1) + DELMB(:,:)        ! update output array [ug/m^3]
      PDFC(:,2:NSPCS+1) = PDFC(:,2:NSPCS+1) + DELMC(:,:)        ! update output array [ug/m^3]

      DO I=1, NSPCS+1
        PDFA(:,I) = MAX( PDFA(:,I), 0.0D-30 )
        PDFB(:,I) = MAX( PDFB(:,I), 0.0D-30 )
      ENDDO

      RETURN
      END SUBROUTINE DISCRETE_INTERCOAG_ABC


      SUBROUTINE DISCRETE_OUT(IUNIT,ICSET,TIMEH)
!-----------------------------------------------------------------------------------------------------------------------
!     11-01-06, DLW
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.
 
      INTEGER, INTENT(   IN) :: IUNIT                 ! output logical unit number [1]
      INTEGER, INTENT(   IN) :: ICSET                 ! identifies test case [1]
      REAL(8), INTENT(   IN) :: TIMEH                 ! model time [h]

      ! Local variables.

      INTEGER :: I
      REAL(8) :: ADNDLOGD,  BDNDLOGD,  CDNDLOGD
      REAL(8) :: ADM1DLOGD, BDM1DLOGD, CDM1DLOGD
      REAL(8) :: ADM2DLOGD, BDM2DLOGD, CDM2DLOGD

      DO I=1, NBINS
        ADNDLOGD  = APDF(I,1) * RDLOGDSC * 1.0D-06           ! convert from [#/m^3] to [#/cm^3]
        BDNDLOGD  = BPDF(I,1) * RDLOGDSC * 1.0D-06           ! convert from [#/m^3] to [#/cm^3]
        CDNDLOGD  = CPDF(I,1) * RDLOGDSC * 1.0D-06           ! convert from [#/m^3] to [#/cm^3]
        ADNDLOGD  = MAX( ADNDLOGD, 1.0D-30 )
        BDNDLOGD  = MAX( BDNDLOGD, 1.0D-30 )
        CDNDLOGD  = MAX( CDNDLOGD, 1.0D-30 )
        ADM1DLOGD = APDF(I,2) * RDLOGDSC                     ! [ug/m^3]
        BDM1DLOGD = BPDF(I,2) * RDLOGDSC                     ! [ug/m^3]
        CDM1DLOGD = CPDF(I,2) * RDLOGDSC                     ! [ug/m^3]
        ADM1DLOGD = MAX( ADM1DLOGD, 1.0D-30 )
        BDM1DLOGD = MAX( BDM1DLOGD, 1.0D-30 )
        CDM1DLOGD = MAX( CDM1DLOGD, 1.0D-30 )
        ADM2DLOGD = APDF(I,3) * RDLOGDSC                     ! [ug/m^3]
        BDM2DLOGD = BPDF(I,3) * RDLOGDSC                     ! [ug/m^3]
        CDM2DLOGD = CPDF(I,3) * RDLOGDSC                     ! [ug/m^3]
        ADM2DLOGD = MAX( ADM2DLOGD, 1.0D-30 )
        BDM2DLOGD = MAX( BDM2DLOGD, 1.0D-30 )
        CDM2DLOGD = MAX( CDM2DLOGD, 1.0D-30 )
        WRITE(IUNIT,91) I, DGRID(I), ADNDLOGD,  BDNDLOGD,  CDNDLOGD,
     &                               ADM1DLOGD, BDM1DLOGD, CDM1DLOGD,
     &                               ADM2DLOGD, BDM2DLOGD, CDM2DLOGD
         
      ENDDO
      
91    FORMAT(I5,10D14.6)
      RETURN
      END SUBROUTINE DISCRETE_OUT


      REAL(8) FUNCTION FLN(X,XG,SIGMAG)
      REAL(8) :: X      ! particle radius or diameter [any units]
      REAL(8) :: XG     ! geometric mean radius or diameter [any units]
      REAL(8) :: SIGMAG ! geometric standard deviation [monodisperse = 1.0]
      REAL(8), PARAMETER :: SQRTTWOPI = 2.506628275D+00
      FLN = EXP(-0.5D+00*(LOG(X/XG)/LOG(SIGMAG))**2) / (X*LOG(SIGMAG)*SQRTTWOPI)
      RETURN
      END FUNCTION FLN

 
      END MODULE AERO_DISCRETE

