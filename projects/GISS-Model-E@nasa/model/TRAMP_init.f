      MODULE AERO_INIT  
!----------------------------------------------------------------------------------------------------------------------
!
!@sum     Defines initial values for the aerosol and gas-phase species for the 
!@+     stand-alone version of the MATRIX microphysical module. 
!@auth  Susanne Bauer/Doug Wright
!----------------------------------------------------------------------------------------------------------------------
      USE AERO_PARAM,  ONLY: CONV_DP_TO_MASS, MINCONC, DENSP, AUNIT1, WRITE_LOG
      USE AERO_CONFIG, ONLY: NWEIGHTS
      USE AERO_SETUP 
      USE AERO_DISCRETE, ONLY: DISCRETE_INIT, ISPCA, ISPCB, ISPCC   
      IMPLICIT NONE
      REAL(8), SAVE :: AERO_IN(NAEROBOX)
      REAL(8), SAVE :: GAS_IN(NGASES) = 0.0D-30   

      CONTAINS


      SUBROUTINE INIT_AERO( ICSET, TEMP, PRES )
!----------------------------------------------------------------------------------------------------------------------
!     Routine to set initial concentrations for various test cases.
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I, INDEX, INDEX1, INDEX2
      INTEGER :: ICSET           ! identifier of the desired set of initial concentrations
      REAL(8) :: TEMP            ! ambient temperature [K]
      REAL(8) :: PRES            ! ambient pressure [Pa]
      REAL(8) :: MPP(NWEIGHTS)   ! mass per particle in mode I [ug]
      REAL(8), PARAMETER :: DEFAULT_N_AKK = 1.0D+10  ! 1.0D+10
      REAL(8), PARAMETER :: DEFAULT_N_ACC = 1.0D+09  ! 1.0D+09
      REAL(8), PARAMETER :: DEFAULT_N_DD1 = 1.0D+07  ! 1.0D+07
      REAL(8), PARAMETER :: DEFAULT_N_DD2 = 1.0D+06  ! 1.0D+06
      REAL(8), PARAMETER :: DEFAULT_N_SSA = 1.0D+08  ! 1.0D+08
      REAL(8), PARAMETER :: DEFAULT_N_SSC = 1.0D+05  ! 1.0D+05
      REAL(8), PARAMETER :: DEFAULT_N_SSS = 1.0D+05  ! 1.0D+05
      REAL(8), PARAMETER :: DEFAULT_N_OCC = 1.0D+08  ! 1.0D+08
      REAL(8), PARAMETER :: DEFAULT_N_BC1 = 1.0D+08  ! 1.0D+08
      
      ! Varibles for the discrete pdf model. 
      
      REAL(8) :: NA,      NB,      NC       ! number concentrations for modes A, B, and C [#/m^3]
      REAL(8) :: DGA,     DGB,     DGC      ! geo. mean diameters   for modes A, B, and C [um]
      REAL(8) :: SIGMAGA, SIGMAGB, SIGMAGC  ! geo. std. deviations  for modes A, B, and C [1]
      REAL(8) :: MASSA,   MASSB,   MASSC    ! mass concentrations   for modes A, B, and C [ug/m^3]
      REAL(8) :: MAFORM,  MBFORM,  MCFORM   ! formula mass concentrations for modes A, B, and C [ug/m^3]
      REAL(8), PARAMETER :: DG_DEFAULT     = 0.08D+00       ! [um]
      REAL(8), PARAMETER :: SIGMAG_DEFAULT = 1.80D+00       ! [1]


      AERO_IN(:) = MINCONC

      !----------------------------------------------------------------------------------------------------------------
      ! Use the default mode lognormal parameters and default particle densities for each mode to get a mean mass 
      !   per particle to obtain initial mass concentrations from initial number concentrations.
      !
      ! DENSPI(:) contains the default density for mode I based on its principal chemical component. 
      !   Parameter CONV_DP_TO_MASS contains the default ambient aerosol density DENSP, which is
      !   divided out in favor of DENSPI(:).  
      !----------------------------------------------------------------------------------------------------------------
      MPP(:) = DENSPI(:) * ( CONV_DP_TO_MASS / DENSP ) 
     &       * ( 1.0D-06 * DGN0(:) )**3 * EXP( 4.5D+00 * ( LOG( SIG0(:) ) )**2 ) ! [ug/particle]

      IF ( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A/)') '   I    MPP(I) [ug]  <--  Initial mean mass per particle in subr. INIT_AERO'
        DO I=1, NWEIGHTS
          WRITE(AUNIT1,'(I4,D15.5,F12.5)') I, MPP(I), DENSPI(I)
        ENDDO
      ENDIF

      !----------------------------------------------------------------------------------------------------------------
      ! Set number concentrations for each mode.
      !----------------------------------------------------------------------------------------------------------------
      SELECT CASE ( ICSET ) 

      CASE( 0 ) 

      CASE( 1 ) 

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK

      CASE( 2 )

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK
        AERO_IN( NUMB_ACC_1 ) = DEFAULT_N_ACC

      CASE( 3 ) 

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK 
        AERO_IN( NUMB_ACC_1 ) = DEFAULT_N_ACC
        AERO_IN( NUMB_DD1_1 ) = DEFAULT_N_DD1
        IF ( NUMB_DD2_1 .GT. 0 ) AERO_IN( NUMB_DD2_1 ) = DEFAULT_N_DD2

      CASE( 4 )  

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK
        AERO_IN( NUMB_ACC_1 ) = DEFAULT_N_ACC
        AERO_IN( NUMB_DD1_1 ) = DEFAULT_N_DD1
        IF ( NUMB_DD2_1 .GT. 0 ) AERO_IN( NUMB_DD2_1 ) = DEFAULT_N_DD2
        IF ( NUMB_SSA_1 .GT. 0 ) AERO_IN( NUMB_SSA_1 ) = DEFAULT_N_SSA
        IF ( NUMB_SSS_1 .GT. 0 ) AERO_IN( NUMB_SSS_1 ) = DEFAULT_N_SSS

      CASE( 5 )  

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK
        AERO_IN( NUMB_ACC_1 ) = DEFAULT_N_ACC
        AERO_IN( NUMB_DD1_1 ) = DEFAULT_N_DD1
        IF ( NUMB_DD2_1 .GT. 0 ) AERO_IN( NUMB_DD2_1 ) = DEFAULT_N_DD2
        IF ( NUMB_SSA_1 .GT. 0 ) AERO_IN( NUMB_SSA_1 ) = DEFAULT_N_SSA
        IF ( NUMB_SSC_1 .GT. 0 ) AERO_IN( NUMB_SSC_1 ) = DEFAULT_N_SSC
        IF ( NUMB_SSS_1 .GT. 0 ) AERO_IN( NUMB_SSS_1 ) = DEFAULT_N_SSS

      CASE( 6 )  

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK
        AERO_IN( NUMB_ACC_1 ) = DEFAULT_N_ACC
        AERO_IN( NUMB_DD1_1 ) = DEFAULT_N_DD1
        IF ( NUMB_DD2_1 .GT. 0 ) AERO_IN( NUMB_DD2_1 ) = DEFAULT_N_DD2
        IF ( NUMB_SSA_1 .GT. 0 ) AERO_IN( NUMB_SSA_1 ) = DEFAULT_N_SSA
        IF ( NUMB_SSC_1 .GT. 0 ) AERO_IN( NUMB_SSC_1 ) = DEFAULT_N_SSC
        IF ( NUMB_SSS_1 .GT. 0 ) AERO_IN( NUMB_SSS_1 ) = DEFAULT_N_SSS
        AERO_IN( NUMB_OCC_1 ) = DEFAULT_N_OCC

      CASE( 7 )  

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK
        AERO_IN( NUMB_ACC_1 ) = DEFAULT_N_ACC
        AERO_IN( NUMB_DD1_1 ) = DEFAULT_N_DD1
        IF ( NUMB_DD2_1 .GT. 0 ) AERO_IN( NUMB_DD2_1 ) = DEFAULT_N_DD2
        IF ( NUMB_SSA_1 .GT. 0 ) AERO_IN( NUMB_SSA_1 ) = DEFAULT_N_SSA
        IF ( NUMB_SSC_1 .GT. 0 ) AERO_IN( NUMB_SSC_1 ) = DEFAULT_N_SSC
        IF ( NUMB_SSS_1 .GT. 0 ) AERO_IN( NUMB_SSS_1 ) = DEFAULT_N_SSS
        AERO_IN( NUMB_OCC_1 ) = DEFAULT_N_OCC
        AERO_IN( NUMB_BC1_1 ) = DEFAULT_N_BC1

      CASE( 8 )  

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK
        AERO_IN( NUMB_ACC_1 ) = DEFAULT_N_ACC * 1.0D+01
        AERO_IN( NUMB_DD1_1 ) = DEFAULT_N_DD1 * 1.0D+01
        IF ( NUMB_DD2_1 .GT. 0 ) AERO_IN( NUMB_DD2_1 ) = DEFAULT_N_DD2 * 1.0D+01
        IF ( NUMB_SSA_1 .GT. 0 ) AERO_IN( NUMB_SSA_1 ) = DEFAULT_N_SSA
        IF ( NUMB_SSC_1 .GT. 0 ) AERO_IN( NUMB_SSC_1 ) = DEFAULT_N_SSC
        IF ( NUMB_SSS_1 .GT. 0 ) AERO_IN( NUMB_SSS_1 ) = DEFAULT_N_SSS
        AERO_IN( NUMB_OCC_1 ) = DEFAULT_N_OCC * 1.0D+01
        AERO_IN( NUMB_BC1_1 ) = DEFAULT_N_BC1 * 1.0D+01

      CASE( 9 )  

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK
        AERO_IN( NUMB_ACC_1 ) = DEFAULT_N_ACC
        AERO_IN( NUMB_DD1_1 ) = DEFAULT_N_DD1
        IF ( NUMB_DD2_1 .GT. 0 ) AERO_IN( NUMB_DD2_1 ) = DEFAULT_N_DD2
        IF ( NUMB_SSA_1 .GT. 0 ) AERO_IN( NUMB_SSA_1 ) = DEFAULT_N_SSA
        IF ( NUMB_SSC_1 .GT. 0 ) AERO_IN( NUMB_SSC_1 ) = DEFAULT_N_SSC
        IF ( NUMB_SSS_1 .GT. 0 ) AERO_IN( NUMB_SSS_1 ) = DEFAULT_N_SSS
        AERO_IN( NUMB_OCC_1 ) = DEFAULT_N_OCC
        AERO_IN( NUMB_BC1_1 ) = DEFAULT_N_BC1

      !----------------------------------------------------------------------------------------------------------------
      ! Case 10 is for testing the discrete model of the PDF. MZJ 2005, Figure 15.2.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 10 ) 

        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        NA      = 1.0D+12
        DGA     = DG_DEFAULT
        SIGMAGA = SIGMAG_DEFAULT
        NB      = 0.0D+00
        DGB     = DG_DEFAULT
        SIGMAGB = SIGMAG_DEFAULT
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        AERO_IN( MASS_AKK_SULF ) = AERO_IN( NUMB_MAP(1) ) * MPP(1)
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 11 is for testing the discrete model of the PDF. MZJ 2005, Figure 15.3.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 11 ) 

        ISPCA   = 1
        ISPCB   = 1
        ISPCC   = 1
        NA      = 1.0D+11
        DGA     = DG_DEFAULT
        SIGMAGA = SIGMAG_DEFAULT
        NB      = 0.0D+00
        DGB     = DG_DEFAULT
        SIGMAGB = SIGMAG_DEFAULT
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        AERO_IN( MASS_AKK_SULF ) = AERO_IN( NUMB_MAP(1) ) * MPP(1)
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 12 is for comparison with results for the discrete model of the PDF. Modes AKK and ACC only. 
      !----------------------------------------------------------------------------------------------------------------
      CASE( 12 ) 

        AERO_IN( NUMB_AKK_1 ) = 1.0D+10     ! [#/m^3]
        AERO_IN( NUMB_ACC_1 ) = 1.0D+09     ! [#/m^3]
        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        INDEX1  = 1
        INDEX2  = 2
        NA      = AERO_IN( NUMB_AKK_1 )
        DGA     = DGN0(INDEX1)
        SIGMAGA = SIG0(INDEX1)
        NB      = AERO_IN( NUMB_ACC_1 )
        DGB     = DGN0(INDEX2)
        SIGMAGB = SIG0(INDEX2)
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        MAFORM = AERO_IN( NUMB_MAP(INDEX1) ) * MPP(INDEX1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', MASSA
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', MAFORM
        AERO_IN( MASS_AKK_SULF ) = MASSA
        MBFORM = AERO_IN( NUMB_MAP(INDEX2) ) * MPP(INDEX2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', MASSB
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', MBFORM
        AERO_IN( MASS_ACC_SULF ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 13 is for comparison with results for the discrete model of the PDF. ACC + BC1 --> BCS.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 13 ) 

        AERO_IN( NUMB_ACC_1 ) = 1.0D+09     ! [#/m^3]
        AERO_IN( NUMB_BC1_1 ) = 1.0D+09     ! [#/m^3]
        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        INDEX1  = 2
        INDEX2  = 10
        NA      = AERO_IN( NUMB_ACC_1 )
        DGA     = DGN0(INDEX1)
        SIGMAGA = SIG0(INDEX1)
        NB      = AERO_IN( NUMB_BC1_1 )
        DGB     = DGN0(INDEX2)
        SIGMAGB = SIG0(INDEX2)
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        MAFORM = AERO_IN( NUMB_MAP(INDEX1) ) * MPP(INDEX1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', MASSA
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', MAFORM
        AERO_IN( MASS_ACC_SULF ) = MASSA
        MBFORM = AERO_IN( NUMB_MAP(INDEX2) ) * MPP(INDEX2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', MASSB
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', MBFORM
        AERO_IN( MASS_BC1_BCAR ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 14 is for comparison with results for the discrete model of the PDF. DD2 + OCC --> MXX.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 14 ) 

        IF ( NUMB_DD2_1 .GT. 0 ) THEN
          AERO_IN( NUMB_DD2_1 ) = 1.0D+07     ! [#/m^3]
        ELSE
          WRITE(*,*)'Cannot use ICSET = 14 with this mechanism - must have mode DD2 to use this ICSET.'
          STOP
        ENDIF
        AERO_IN( NUMB_OCC_1 ) = 1.0D+09     ! [#/m^3]
        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        INDEX1  = 5
        INDEX2  = 9
        NA      = AERO_IN( NUMB_DD2_1 )
        DGA     = DGN0(INDEX1)
        SIGMAGA = SIG0(INDEX1)
        NB      = AERO_IN( NUMB_OCC_1 )
        DGB     = DGN0(INDEX2)
        SIGMAGB = SIG0(INDEX2)
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        MAFORM = AERO_IN( NUMB_MAP(INDEX1) ) * MPP(INDEX1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', MASSA
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', MAFORM
        AERO_IN( MASS_DD2_DUST ) = MASSA
        MBFORM = AERO_IN( NUMB_MAP(INDEX2) ) * MPP(INDEX2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', MASSB
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', MBFORM
        AERO_IN( MASS_OCC_OCAR ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 15 is for comparison with results for the discrete model of the PDF. AKK + BC1 --> BCS.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 15 ) 

        AERO_IN( NUMB_AKK_1 ) = 1.0D+10     ! [#/m^3]
        AERO_IN( NUMB_BC1_1 ) = 1.0D+09     ! [#/m^3]
        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        INDEX1  = 1
        INDEX2  = 10
        NA      = AERO_IN( NUMB_AKK_1 )
        DGA     = DGN0(INDEX1)
        SIGMAGA = SIG0(INDEX1)
        NB      = AERO_IN( NUMB_BC1_1 )
        DGB     = DGN0(INDEX2)
        SIGMAGB = SIG0(INDEX2)
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        MAFORM = AERO_IN( NUMB_MAP(INDEX1) ) * MPP(INDEX1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', MASSA
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', MAFORM
        AERO_IN( MASS_AKK_SULF ) = MASSA
        MBFORM = AERO_IN( NUMB_MAP(INDEX2) ) * MPP(INDEX2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', MASSB
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', MBFORM
        AERO_IN( MASS_BC1_BCAR ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 16 is for comparison with results for the discrete model of the PDF. DD1 + OCC --> MXX.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 16 ) 

        AERO_IN( NUMB_DD1_1 ) = 1.0D+09     ! [#/m^3]
        AERO_IN( NUMB_OCC_1 ) = 1.0D+09     ! [#/m^3]
        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        INDEX1  = 3
        INDEX2  = 9
        NA      = AERO_IN( NUMB_DD1_1 )
        DGA     = DGN0(INDEX1)
        SIGMAGA = SIG0(INDEX1)
        NB      = AERO_IN( NUMB_OCC_1 )
        DGB     = DGN0(INDEX2)
        SIGMAGB = SIG0(INDEX2)
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        MAFORM = AERO_IN( NUMB_MAP(INDEX1) ) * MPP(INDEX1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', MASSA
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', MAFORM
        AERO_IN( MASS_DD1_DUST ) = MASSA
        MBFORM = AERO_IN( NUMB_MAP(INDEX2) ) * MPP(INDEX2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', MASSB
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', MBFORM
        AERO_IN( MASS_OCC_OCAR ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 17 is for comparison with results for the discrete model of the PDF. DD1 + SSC --> MXX.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 17 ) 

        AERO_IN( NUMB_DD1_1 ) = 1.00D+09     ! [#/m^3]
        AERO_IN( NUMB_SSC_1 ) = 0.40D+06     ! [#/m^3]
        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        INDEX1  = 3
        INDEX2  = 8
        NA      = AERO_IN( NUMB_DD1_1 )
        DGA     = DGN0(INDEX1)
        SIGMAGA = SIG0(INDEX1)
        NB      = AERO_IN( NUMB_SSC_1 )
        DGB     = DGN0(INDEX2)
        SIGMAGB = SIG0(INDEX2)
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        MAFORM = AERO_IN( NUMB_MAP(INDEX1) ) * MPP(INDEX1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', MASSA
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', MAFORM
        AERO_IN( MASS_DD1_DUST ) = MASSA
        MBFORM = AERO_IN( NUMB_MAP(INDEX2) ) * MPP(INDEX2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', MASSB
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', MBFORM
        AERO_IN( MASS_SSC_SEAS ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 18 is for comparison with the discrete model of the PDF. Mode AKK. Intramodal coagulation only.  
      !----------------------------------------------------------------------------------------------------------------
      CASE( 18 ) 

        IF ( NUMB_AKK_1 .GT. 0 ) AERO_IN( NUMB_AKK_1 ) = DEFAULT_N_AKK
        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        INDEX1  = 1
        INDEX2  = 2
        NA      = AERO_IN( NUMB_AKK_1 )
        DGA     = DGN0(INDEX1)
        SIGMAGA = SIG0(INDEX1)
        NB      = AERO_IN( NUMB_ACC_1 )
        DGB     = DGN0(INDEX2)
        SIGMAGB = SIG0(INDEX2)
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        MAFORM = AERO_IN( NUMB_MAP(INDEX1) ) * MPP(INDEX1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', MASSA
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', MAFORM
        AERO_IN( MASS_AKK_SULF ) = MASSA
        MBFORM = AERO_IN( NUMB_MAP(INDEX2) ) * MPP(INDEX2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', MASSB
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', MBFORM
        AERO_IN( MASS_ACC_SULF ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 19 is for comparison with results for the discrete model of the PDF. BC1 + OCC --> BOC.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 19 ) 

        AERO_IN( NUMB_OCC_1 ) = 1.0D+09     ! [#/m^3]
        AERO_IN( NUMB_BC1_1 ) = 1.0D+09     ! [#/m^3]
        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        INDEX1  = 9
        INDEX2  = 10
        NA      = AERO_IN( NUMB_OCC_1 )
        DGA     = DGN0(INDEX1)
        SIGMAGA = SIG0(INDEX1)
        NB      = AERO_IN( NUMB_BC1_1 )
        DGB     = DGN0(INDEX2)
        SIGMAGB = SIG0(INDEX2)
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        MAFORM = AERO_IN( NUMB_MAP(INDEX1) ) * MPP(INDEX1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', MASSA
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', MAFORM
        AERO_IN( MASS_OCC_OCAR ) = MASSA
        MBFORM = AERO_IN( NUMB_MAP(INDEX2) ) * MPP(INDEX2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', MASSB
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', MBFORM
        AERO_IN( MASS_BC1_BCAR ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 20 is for comparison with results for the discrete model of the PDF. ACC + DD1 --> DD1. (no DS1 transfer)
      !----------------------------------------------------------------------------------------------------------------
      CASE( 20 ) 

        AERO_IN( NUMB_ACC_1 ) = 1.0D+09     ! [#/m^3]
        AERO_IN( NUMB_DD1_1 ) = 1.0D+09     ! [#/m^3]
        ISPCA   = 2
        ISPCB   = 2
        ISPCC   = 2
        INDEX1  = 2
        INDEX2  = 3
        NA      = AERO_IN( NUMB_ACC_1 )
        DGA     = DGN0(INDEX1)
        SIGMAGA = SIG0(INDEX1)
        NB      = AERO_IN( NUMB_DD1_1 )
        DGB     = DGN0(INDEX2)
        SIGMAGB = SIG0(INDEX2)
        NC      = 0.0D+00
        DGC     = DG_DEFAULT
        SIGMAGC = SIGMAG_DEFAULT
        CALL DISCRETE_INIT(ICSET,NA,DGA,SIGMAGA,NB,DGB,SIGMAGB,NC,DGC,SIGMAGC,TEMP,PRES,MASSA,MASSB,MASSC)
        MAFORM = AERO_IN( NUMB_MAP(INDEX1) ) * MPP(INDEX1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', MASSA
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', MAFORM
        AERO_IN( MASS_ACC_SULF ) = MASSA
        MBFORM = AERO_IN( NUMB_MAP(INDEX2) ) * MPP(INDEX2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', MASSB
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', MBFORM
        AERO_IN( MASS_DD1_DUST ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 21 is for illustration of aerosol activation. 
      !----------------------------------------------------------------------------------------------------------------
      CASE( 21 )  

        IF ( NUMB_AKK_1 .GT. 0 ) THEN 
          AERO_IN( NUMB_AKK_1 ) = 1.0D+09 * 2.0D+00
          AERO_IN( NUMB_ACC_1 ) = 1.0D+08 * 2.0D+00
        ELSE
          AERO_IN( NUMB_ACC_1 ) = 1.1D+09 * 2.0D+00
        ENDIF
        IF ( NUMB_DD2_1 .GT. 0 ) THEN
          AERO_IN( NUMB_DD1_1 ) = 1.0D+08 * 0.1D+00
          AERO_IN( NUMB_DD2_1 ) = 1.0D+07 * 0.1D+00
        ELSE
          AERO_IN( NUMB_DD1_1 ) = 1.1D+08 * 0.1D+00
        ENDIF
        IF ( NUMB_SSA_1 .GT. 0 ) AERO_IN( NUMB_SSA_1 ) = 1.000D+08 * 0.2D+00
        IF ( NUMB_SSC_1 .GT. 0 ) AERO_IN( NUMB_SSC_1 ) = 0.004D+08 * 0.2D+00
        IF ( NUMB_SSS_1 .GT. 0 ) AERO_IN( NUMB_SSS_1 ) = 1.004D+08 * 0.2D+00
        AERO_IN( NUMB_OCC_1 ) = 1.0D+08 * 1.0D+01
        AERO_IN( NUMB_BC1_1 ) = 1.0D+08 * 1.0D+01

      CASE DEFAULT

        WRITE(*,*)'BAD VALUE OF ICSET IN SUBR. INIT_AERO'
        STOP

      END SELECT

      !----------------------------------------------------------------------------------------------------------------
      ! Calculate masses for all primary modes.
      !----------------------------------------------------------------------------------------------------------------
      DO I=1, NMODES
 
        ! Set the defining chemical species for each mode to receive the initial concentrations. 

        INDEX = 0
        IF( MODE_NAME(I) .EQ. 'AKK' ) INDEX = MASS_AKK_SULF
        IF( MODE_NAME(I) .EQ. 'ACC' ) INDEX = MASS_ACC_SULF
        IF( MODE_NAME(I) .EQ. 'DD1' ) INDEX = MASS_DD1_DUST
        IF( MODE_NAME(I) .EQ. 'DD2' ) INDEX = MASS_DD2_DUST
        IF( MODE_NAME(I) .EQ. 'SSA' ) INDEX = MASS_SSA_SEAS
        IF( MODE_NAME(I) .EQ. 'SSC' ) INDEX = MASS_SSC_SEAS
        IF( MODE_NAME(I) .EQ. 'SSS' ) INDEX = MASS_SSS_SEAS
        IF( MODE_NAME(I) .EQ. 'OCC' ) INDEX = MASS_OCC_OCAR
        IF( MODE_NAME(I) .EQ. 'BC1' ) INDEX = MASS_BC1_BCAR

        IF( INDEX .GT. 0 ) THEN

          ! WRITE(*,*)'INDEX = ', INDEX
          AERO_IN( INDEX ) = AERO_IN( NUMB_MAP(I) ) * MPP(I)

          ! For some IC sets, initial sulfate concentrations are also given to these modes.

          IF( ICSET.EQ.9 .AND. MODE_NAME(I).NE.'SSC' ) AERO_IN( SULF_MAP(I) ) = AERO_IN( NUMB_MAP(I) ) * MPP(I)

        ENDIF
      ENDDO

      ! Set ammonium assuming ammonium sulfate. 

      IF( ICSET .EQ. 21 ) AERO_IN( 2 ) = 0.3755532793D+00 * SUM( AERO_IN(SULF_MAP(:)) ) 

      RETURN
      END SUBROUTINE INIT_AERO


      END MODULE AERO_INIT

