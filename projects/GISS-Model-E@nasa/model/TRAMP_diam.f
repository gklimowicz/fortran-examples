      MODULE AERO_DIAM
      USE AERO_PARAM,  ONLY: NLAYS, AUNIT1, WRITE_LOG
      USE AERO_CONFIG, ONLY: NMODES, MECH
      USE AERO_SETUP,  ONLY: N_DP_CONDTABLE, DP_CONDTABLE, MODE_NAME
      USE AMP_AEROSOL, ONLY: DIAM
!-------------------------------------------------------------------------------------------------------------------------
!     The array DIAM(x,y,z) contains current values of some measure of average ambient mode diameter for each mode  
!       for use outside of the MATRIX microphysical module where it is calculated. 
!       Values in DIAM are saved at the top of the subr. MATRIX before microphysical evolution 
!       for the current time step is done. 
!
!     The current measure of particle diameter is the diameter of average mass:
! 
!        DIAM(x,y,z) = [ (6/pi) * (Mi/Ni) * (1/D) ]^(1/3)
!
!     with Mi the total mass concentration (including water) in mode i, Ni the number concentration in mode i, and
!     D is the particle density. 
!-------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
!      REAL(8) :: DIAM_HISTOGRAM(NMODES,N_DP_CONDTABLE,2)    ! [1]

      CONTAINS

      SUBROUTINE SETUP_DIAM
!-------------------------------------------------------------------------------------------------------------------------
!@sum     Routine to initialize DIAM before model time stepping is begun.
!@+     This routine should be called only once. 
!@+     The diameter of average mass is the cube root of the normalized third diameter moment: (M3/M0)**(1/3) = Dg*Sg**3
!@auth Susanne Bauer/Doug Wright
!-------------------------------------------------------------------------------------------------------------------------
      USE AERO_SETUP, ONLY: DGN0, SIG0
      INTEGER :: I
      REAL(8) :: SG, D_NUMBER_MEAN
      IF( WRITE_LOG ) WRITE(AUNIT1,'(/A/)') 'I,DGN0(I),SIG0(I),SG,DIAM(1,1,1,I)*1.0D+06,D_NUMBER_MEAN'
      DO I=1, NMODES
        SG = EXP( 0.5d+00*( LOG(SIG0(I)) )**2 )
        DIAM(:,:,:,I) = 1.0D-06 * DGN0(I)*SG**3           ! convert from [um] to [m]
        D_NUMBER_MEAN =           DGN0(I)*SG             
        IF( WRITE_LOG ) WRITE(AUNIT1,90) I,DGN0(I),SIG0(I),SG,DIAM(1,1,1,I)*1.0D+06,D_NUMBER_MEAN
      ENDDO
      WRITE(AUNIT1,'(A)') '  '

      ! Zero histogram for DIAM values.
      
!      DIAM_HISTOGRAM(:,:,:) = 0.0D+00 

90    FORMAT(I5,6F12.5)
      RETURN
      END SUBROUTINE SETUP_DIAM


!      SUBROUTINE WRITE_DIAM_HISTOGRAM
!!-------------------------------------------------------------------------------------------------------------------------
!!     Call this routine at the end of all time stepping.
!!-------------------------------------------------------------------------------------------------------------------------
!      INTEGER :: I, N
!      INTEGER, PARAMETER :: OUTUNIT  = 93
!      INTEGER, PARAMETER :: NBINSSUM =  5
!      REAL(8) :: SUM_HIST(NMODES),SUM_HIST_DRY(NMODES),D,F16,F18,F20
!
!      SELECT CASE( MECH ) 
!      CASE( 1 )
!        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech1_diam_histogram.plt')
!      CASE( 2 )
!        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech2_diam_histogram.plt')
!      CASE( 3 )
!        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech3_diam_histogram.plt')
!      CASE( 4 )
!        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech4_diam_histogram.plt')
!      CASE( 5 )
!        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech5_diam_histogram.plt')
!      CASE( 6 )
!        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech6_diam_histogram.plt')
!      CASE( 7 )
!        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech7_diam_histogram.plt')
!      CASE( 8 )
!        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech8_diam_histogram.plt')
!      CASE DEFAULT
!        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mechx_diam_histogram.plt')
!      END SELECT 
!
!      !-------------------------------------------------------------------------------------------------------------------
!      ! F16, F18, and F20 convert the diameter of average mass (D) to the number mean diameter for 
!      ! geometric standard deviations of 1.6, 1.8, and 2.0, respectively. 
!      !-------------------------------------------------------------------------------------------------------------------
!      F16 = 1.0D+00 / ( EXP( ( LOG(1.6D+00) )**2 ) ) 
!      F18 = 1.0D+00 / ( EXP( ( LOG(1.8D+00) )**2 ) ) 
!      F20 = 1.0D+00 / ( EXP( ( LOG(2.0D+00) )**2 ) ) 
!
!      !-------------------------------------------------------------------------------------------------------------------
!      ! For each mode, normalize the histogram to unity. 
!      !-------------------------------------------------------------------------------------------------------------------
!      DO I=1, NMODES
!        SUM_HIST(I) = 0.0D+00
!        DO N=1, N_DP_CONDTABLE
!          SUM_HIST(I) = SUM_HIST(I) + DIAM_HISTOGRAM(I,N,1)
!        ENDDO
!        DIAM_HISTOGRAM(I,:,1) = DIAM_HISTOGRAM(I,:,1) / ( SUM_HIST(I) + 1.0D-30 ) 
!        DIAM_HISTOGRAM(I,:,2) = DIAM_HISTOGRAM(I,:,2) / ( SUM_HIST(I) + 1.0D-30 ) 
!        WRITE(OUTUNIT,'(A,I3,1X,A,F15.1)') 'Total count for mode I = ',I,' is ', SUM_HIST(I)
!      ENDDO
!
!      WRITE(OUTUNIT,91) '         DAM','       DNM16','       DNM18','       DNM20',
!     &                                 '       DGN16','       DGN18','       DGN20',MODE_NAME(:),MODE_NAME(:)
!      SUM_HIST    (:) = 0.0D+00
!      SUM_HIST_DRY(:) = 0.0D+00
!      D = 1.0D+00
!      N = 0
!      DO I=1, N_DP_CONDTABLE
!        SUM_HIST    (1:NMODES) = SUM_HIST    (1:NMODES) + DIAM_HISTOGRAM(1:NMODES,I,1) + 1.0D-20
!        SUM_HIST_DRY(1:NMODES) = SUM_HIST_DRY(1:NMODES) + DIAM_HISTOGRAM(1:NMODES,I,2) + 1.0D-20
!        D = D*DP_CONDTABLE(I)
!        IF( MOD( I, NBINSSUM ) .EQ. 0 ) THEN
!          N = N + 1
!          D = 1.0D+06 * D**(1.0D+00/REAL(NBINSSUM))
!          WRITE(OUTUNIT,90)N,I,D,D*F16,D*F18,D*F20,D*F16**1.5,D*F18**1.5,D*F20**1.5,SUM_HIST(1:NMODES),SUM_HIST_DRY(1:NMODES)
!          SUM_HIST    (:) = 0.0D+00
!          SUM_HIST_DRY(:) = 0.0D+00
!          D = 1.0D+00
!        ENDIF
!      ENDDO
!      CLOSE(OUTUNIT)
!90    FORMAT(2I5,7F12.7,32F14.10)
!91    FORMAT(10X,7A12,32A14)
!      RETURN
!      END SUBROUTINE WRITE_DIAM_HISTOGRAM

      END MODULE AERO_DIAM
 
