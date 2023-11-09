      MODULE AERO_NPF
!-------------------------------------------------------------------------------------------------------------------
!
!@sum     This module contains all sub-programs to calculate nucleation and 
!@+     new particle formation rates.
!@auth  Susanne Bauer/Doug Wright
!-------------------------------------------------------------------------------------------------------------------
      USE AERO_PARAM
      USE AERO_SETUP, ONLY: DIFFCOEF_M2S, AVG_DP_OF_AVG_MASS_METERS
      IMPLICIT NONE
      !-------------------------------------------------------------------------------------------------------------
      ! Select the nucleation scheme: INUC = 1, JVM 
      !                               INUC = 2, VEHKAMAKI
      !                               INUC = 3, NAPARI
      !                               INUC = 4, EISELE AND MCMURRY, 1997
      !-------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: INUC = 2
      LOGICAL, PARAMETER :: INCLUDE_ION_ION = .TRUE.      ! include Turco scheme

      !-------------------------------------------------------------------------------------------------------------
      ! Flag for use of the Kerminem and Kulmala (2002) parameterization for
      ! conversion of a nucleation rate to a particle formation rate at a
      ! larger user-selected size.
      !-------------------------------------------------------------------------------------------------------------
      LOGICAL, PARAMETER :: KK02 = .TRUE.
      LOGICAL, PARAMETER :: WRITE_F_KK02 = .FALSE.

      !-------------------------------------------------------------------------------------------------------------
      ! New particle parameters.
      !-------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DSTAR_NM =  1.0D+00   ! diameter of critical nucleus [nm]
      REAL(8), PARAMETER :: DNPF_NM  =  3.0D+00   ! diameter of E&M(1997) new particles [nm]
      REAL(8), PARAMETER :: DNU_NM   =  3.0D+00   ! diameter of a new particle [nm]
      REAL(8), PARAMETER :: RNU_NM = DNU_NM*0.5D+00 ! radius of a new particle [m]
      REAL(8), PARAMETER :: DNU = DNU_NM*1.0D-09  ! diameter of a new particle [m]
      REAL(8), PARAMETER :: VNU = PI6*DNU*DNU*DNU ! volume of a new particle   [m^3]
      !-------------------------------------------------------------------------------------------------------------
      ! For conversion from 1-nm particles to DNU_NM-nm particles.
      !-------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: CONV_NUCL_TO_NPF = 1.0D+06    ! [#/cm^3/s] to [#/m^3/s]
     &         * DSTAR_NM * DSTAR_NM * DSTAR_NM / ( DNU_NM * DNU_NM * DNU_NM )
      !-------------------------------------------------------------------------------------------------------------
      ! For conversion from DNPF_NM-nm particles to DNU_NM-nm particles.
      !-------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: CONV_EMNPF_TO_NPF = 1.0D+06   ! [#/cm^3/s] to [#/m^3/s]
     &         * DNPF_NM  * DNPF_NM  * DNPF_NM  / ( DNU_NM * DNU_NM * DNU_NM )

      !-------------------------------------------------------------------------------------------------------------
      ! Limits on the nucleation rate applied to all parameterizations,
      ! except that the upper limit is not applied to the Vehkamaki et al. 2002
      ! parameterization for binary H2SO4-H2O nucleation.
      !-------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: J_LOWER = 1.0D-07    ! [#/cm^3/s]
      REAL(8), PARAMETER :: J_UPPER = 1.0D+07    ! [#/cm^3/s]

      !-------------------------------------------------------------------------------------------------------------
      ! NPFMASS(I,N) is the sulfate mass [ugSO4] per new particle of volume VNU for I=NINT[RH(%)].
      !
      ! Particle composition is presently divided into three regimes:
      !   (1) ammonium sulfate   - second index set to 2
      !   (2) ammonium bisulfate - second index set to 1
      !   (3) sulfuric acid      - second index set to 0
      !-------------------------------------------------------------------------------------------------------------
      REAL(8), SAVE :: NPFMASS(0:100,0:2)   ! [ugSO4/particle]
      INTEGER, SAVE :: NPFMASS_REGIME = 2   ! second index to NPFMASS
      CONTAINS
 

      SUBROUTINE NPFRATE(PRS,RH,TEMP,XH2SO4,SO4RATE,XNH3,KC,DNDT,DMDT_SO4,ICALL)
!-------------------------------------------------------------------------------------------------------------------
!     DLW 2006.           
!     Routine to calculate the rate of production of new particles and the 
!     corresponding rate of production of particulate sulfate mass.
!-------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8), INTENT(IN)   :: PRS       ! pressure [Pa]
      REAL(8), INTENT(IN)   :: RH        ! fractional relative humidity [1]
      REAL(8), INTENT(IN)   :: TEMP      ! ambient temperature [K]
      REAL(8), INTENT(IN)   :: XH2SO4    ! sulfuric acid (as SO4) concentration [ugSO4/m^3]
      REAL(8), INTENT(IN)   :: SO4RATE   ! gas-phase H2SO4 (as SO4) production rate [ugSO4/m^3 s]
      REAL(8), INTENT(IN)   :: XNH3      ! ammonia mixing ratio [ppmV]
      REAL(8), INTENT(IN)   :: KC        ! condensational sink [1/s]
      INTEGER, INTENT(IN)   :: ICALL     ! flag signalling type of call

      ! Output arguments.

      REAL(8), INTENT(OUT)  :: DNDT      ! particle number production rate [m^-3 s^-1]
      REAL(8), INTENT(OUT)  :: DMDT_SO4  ! SO4 mass production rate        [ugSO4/m^3 s]

      ! Scratch local variables.

      REAL(8) :: JTOT        ! total nucleation rate [cm^-3 s^-1]
      REAL(8) :: JGAS        ! homogeneous nucleation rate [cm^-3 s^-1]
      REAL(8) :: JION        ! ion-ion recombination nucleation rate [cm^-3 s^-1]
      REAL(8) :: SO4MASS     ! mass of SO4 per new particle [ugSO4]
      REAL(8) :: H2SO4_TMP   ! [H2SO4] scratch variable [molecules/cm^3]
      REAL(8) :: NH3_TMP     ! [NH3]   scratch variable [ppt]                 
      REAL(8) :: JP_TO_J     ! ratio of the particle formation rate to the nucleation rate [1]

      !-------------------------------------------------------------------------------------------------------------
      ! In the call to NUCLEATION_RATE, RH is converted to [%], the sulfuric
      ! acid (SO4) concentration is converted from ugSO4/m^3 to molecules/cm^3,
      ! and ammonia is converted from ppm to ppt.
      !-------------------------------------------------------------------------------------------------------------
      H2SO4_TMP = UGM3_NCM3 * MAX ( XH2SO4, 1.0D-30 )     ! [molecule/cm^3]
      NH3_TMP   =   1.0D+06 * MAX ( XNH3  , 1.0D-30 )     ! [ppt]

      CALL NUCLEATION_RATE( TEMP, 1.0D+02*RH, H2SO4_TMP, NH3_TMP,
     &                      ZHEIGHT(ILAY), JTOT, JGAS, JION )

      ! WRITE(78,*)'TEMP, 100.0D+00*RH, H2SO4_TMP, NH3_TMP, ZHEIGHT(ILAY), JTOT, JGAS, JION'
      ! WRITE(78,*) TEMP, 100.0D+00*RH, H2SO4_TMP, NH3_TMP, ZHEIGHT(ILAY), JTOT, JGAS, JION  

      !-------------------------------------------------------------------------------------------------------------
      ! Convert the nucleation rate in [#/cm^3/s] into a new particle formation rate for 
      !   particles of diameter DNU_NM in [#/m^3/s].
      !
      ! There are two options, either the Kerminen and Kulmala (2002) parameterization, 
      !   simply reduce the nucleation rate based upon mass conservation as 
      !   1-nm diameter particles to converted to particles at the selected size. 
      !
      ! If the Eisele and McMurry (1997) curves are used, the conversion is
      !   from 3.0 nm (rather than 1 nm) to the selected size. 
      !-------------------------------------------------------------------------------------------------------------
      IF( INUC.EQ.1 .OR. INUC.EQ.2 .OR. INUC.EQ.3 ) THEN
        IF( KK02 ) THEN
          CALL F_KK02( PRS, TEMP, 1.0D+02*RH, H2SO4_TMP, NH3_TMP, KC, DNU_NM, DSTAR_NM, JP_TO_J )
          DNDT = ( 1.0D+06 * JTOT ) * JP_TO_J ! Convert [#/cm^3/s] to [#/m^3/s].
        ELSE
          DNDT = CONV_NUCL_TO_NPF * JTOT      ! CONV_NUCL_TO_NPF includes the [#/cm^3/s] to [#/m^3/s] conversion. 
        ENDIF
      ELSEIF( INUC.EQ.4 ) THEN    ! Eisele and McMurry (1997) curves are used.
        IF( KK02 ) THEN
          CALL F_KK02( PRS, TEMP, 1.0D+02*RH, H2SO4_TMP, NH3_TMP, KC, DNU_NM, DNPF_NM, JP_TO_J )
          DNDT = ( 1.0D+06 * JTOT ) * JP_TO_J ! Convert [#/cm^3/s] to [#/m^3/s].
        ELSE
          DNDT = CONV_EMNPF_TO_NPF * JTOT     ! CONV_EMNPF_TO_NPF includes the [#/cm^3/s] to [#/m^3/s] conversion.
        ENDIF
      ENDIF
      ! WRITE(34,'(A,4D13.5)')'IN NPFRATE: JTOT, JP_TO_J, DNDT, H2SO4 = ', JTOT, JP_TO_J, DNDT, H2SO4_TMP

      !-------------------------------------------------------------------------------------------------------------
      ! Calculate mass production rate [ugSO4/m^3/s)], limited by the
      ! production rate by the production rate of H2SO4.
      ! Adjust the number production rate if necessary.
      !
      ! NPFMASS(I,N) is the sulfate mass [ugSO4] per new particle of volume VNU for I=NINT[RH(%)].
      !
      ! Particle composition is presently divided into three regimes:--> 
      !   (1) ammonium sulfate   --> second index set to 2 -->  NPFMASS_REGIME = 2
      !   (2) ammonium bisulfate --> second index set to 1 -->  NPFMASS_REGIME = 1
      !   (3) sulfuric acid      --> second index set to 0 -->  NPFMASS_REGIME = 0
      !-------------------------------------------------------------------------------------------------------------
      SO4MASS = NPFMASS( NINT( 100.0D+00*RH ), NPFMASS_REGIME )   ! [ugSO4]
      DMDT_SO4 = SO4MASS * DNDT                                   ! [ugSO4/m^3/s]
      IF( ICALL .GT. 0 ) RETURN                                   ! do not impose mass limitation for this call
      IF ( DMDT_SO4 .GT. SO4RATE ) THEN                           ! [ugSO4/m^3/s]
        ! IF ( DMDT_SO4 .GT. 1.0D-10 ) WRITE(34,*)'NPFRATE MASS-LIMIT IMPOSED: DMDT_SO4, SO4RATE=',DMDT_SO4,SO4RATE
        DMDT_SO4 = SO4RATE                                        ! [ugSO4/m^3/s]
        DNDT = DMDT_SO4 / SO4MASS                                 ! [   # /m^3/s]
      ENDIF

      RETURN
      END SUBROUTINE NPFRATE


      SUBROUTINE SETUP_NPFMASS
!-------------------------------------------------------------------------------------------------------------------
!     DLW 2006.
!     Routine to pre-calculate the SO4 mass in a single new particle as
!     a function of RH and new particle volume, the array NPFMASS(0:100,0:2).
!-------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I 
      REAL(8) :: H, XI_KE_SULFATE, XI_KE_H2SO4
      !-------------------------------------------------------------------------------------------------------------
      ! CONV_V_TO_SO4_XXX is the mass of sulfate (MW=96g/mol) in [ugSO4] for either XXX=SULFATE or XXX=H2SO4.
      !   in a dry particle of volume VNU [m^3]. Checked on 7-27-06.
      !   In CONV_V_TO_SO4, the 1.0D+12 converts [gSO4/cm^3] to [ugSO4/m^3].
      !
      ! RHC_NH42SO4_RNU_NM is the crystallization RH for a 5-nm (radius) dry particle (E. Lewis).
      !-------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: CONV_V_TO_SO4_SULFATE = RHO_NH42SO4 * 1.0D+12 * MW_SO4 / MW_NH42SO4
      REAL(8), PARAMETER :: CONV_V_TO_SO4_H2SO4   = RHO_H2SO4   * 1.0D+12 * MW_SO4 / MW_H2SO4
      REAL(8), PARAMETER :: RHC_NH42SO4_DNU_03NM  = 0.00D+00
      REAL(8), PARAMETER :: RHC_NH42SO4_DNU_10NM  = 0.50D+00
      REAL(8), PARAMETER :: RHC_NH42SO4_DNU_20NM  = 0.50D+00
      REAL(8), PARAMETER :: RHC_H2SO4_DNU_03NM    = 0.10D+00
      REAL(8), PARAMETER :: RHC_H2SO4_DNU_10NM    = 0.05D+00
      REAL(8), PARAMETER :: RHC_H2SO4_DNU_20NM    = 0.05D+00

      !-------------------------------------------------------------------------------------------------------------
      ! NPFMASS(I,J) is the sulfate mass [ugSO4] per new particle of volume VNU for I=NINT[RH(%)].
      !
      ! Particle composition is presently divided into three regimes:
      !   (1) ammonium sulfate   - second index J = 2
      !   (2) ammonium bisulfate - second index J = 1
      !   (3) sulfuric acid      - second index J = 0
      ! 
      ! The Kelvin Effect is taken into account.
      !
      ! XI_KE is the radius ratio r_ambient/r_dry, including the Kelvin effect, and is used to convert
      !   ambient particle volume to dry particle volume.
      !
      ! CONV_V_TO_SO4 converts dry volume ammonium sulfate [m^3] to dry mass sulfate (MW=96g/mol) [ugSO4].
      !-------------------------------------------------------------------------------------------------------------
      IF( WRITE_LOG ) WRITE(AUNIT1,90)
      DO I=0, 100
        H = MIN( 1.0D-02 * DBLE(I), 0.999D+00 )
        !-----------------------------------------------------------------------------------------------------------
        ! Ammonium sulfate. 
        ! Ammonium bisulfate. The sulfate value is also used because of lack of data for bisulfate solutions.
        !-----------------------------------------------------------------------------------------------------------
        IF    ( DNU_NM .EQ.  3.0D+00 ) THEN
          IF ( H .GE. RHC_NH42SO4_DNU_03NM ) THEN     ! The particle is wet.
            XI_KE_SULFATE = 1.0D+00 + 0.2D+00*H       ! Linear fit over the entire RH range.
          ELSE
            XI_KE_SULFATE = 1.0D+00 + 0.2D+00*H       ! Linear fit over the entire RH range.
          ENDIF
        ELSEIF( DNU_NM .EQ. 10.0D+00 ) THEN
          IF ( H .GE. RHC_NH42SO4_DNU_10NM ) THEN     ! The particle is wet.
            XI_KE_SULFATE = 0.677D+00 + H*(1.816D+00 + H*( -2.345D+00 + H*1.296D+00 ) )
          ELSE                                        ! The particle is dry. 
            XI_KE_SULFATE = 1.0D+00 + 0.16075D+00*(H/RHC_NH42SO4_DNU_10NM)
          ENDIF
        ELSEIF( DNU_NM .EQ. 20.0D+00 ) THEN
          IF ( H .GE. RHC_NH42SO4_DNU_20NM ) THEN     ! The particle is wet.
            XI_KE_SULFATE = 0.175D+00 + H*(4.532D+00 + H*( -6.894D+00 + H*3.856D+00 ) )
          ELSE                                        ! The particle is dry. 
            XI_KE_SULFATE = 1.0D+00 + 0.1995D+00*(H/RHC_NH42SO4_DNU_20NM)
          ENDIF
        ELSE
          WRITE(*,*)'Bad value of DNU_NM in subr. SETUP_NPFMASS: DNU_NM = ', DNU_NM
          STOP
        ENDIF
        NPFMASS(I,2) = ( VNU / XI_KE_SULFATE**3 ) * CONV_V_TO_SO4_SULFATE
        NPFMASS(I,1) = NPFMASS(I,2)
        !-----------------------------------------------------------------------------------------------------------
        ! Sulfuric acid.
        !-----------------------------------------------------------------------------------------------------------
        IF    ( DNU_NM .EQ.  3.0D+00 ) THEN
          IF ( H .GE. RHC_H2SO4_DNU_03NM ) THEN                 
            XI_KE_H2SO4 = 1.14D+00 + H*(0.464D+00 + H*( -0.336D+00 + H*0.189D+00 ) )
          ELSE
            XI_KE_H2SO4 = 1.0D+00 + 0.179D+00*(H/RHC_H2SO4_DNU_03NM) ! Linear fit over this range.
          ENDIF
        ELSEIF( DNU_NM .EQ. 10.0D+00 ) THEN
          IF ( H .GE. RHC_H2SO4_DNU_10NM ) THEN                
            XI_KE_H2SO4 = 1.14D+00 + H*(0.765D+00 + H*( -0.850D+00 + H*0.745D+00 ) )
          ELSE                                    
            XI_KE_H2SO4 = 1.0D+00 + 0.211D+00*(H/RHC_H2SO4_DNU_10NM) ! Linear fit over this range.
          ENDIF
        ELSEIF( DNU_NM .EQ. 20.0D+00 ) THEN
          IF ( H .GE. RHC_H2SO4_DNU_20NM ) THEN                
            XI_KE_H2SO4 = 1.113D+00 + H*(1.190D+00 + H*( -2.001D+00 + H*1.750D+00 ) )
          ELSE                                    
            XI_KE_H2SO4 = 1.0D+00 + 0.168D+00*(H/RHC_H2SO4_DNU_20NM) ! Linear fit over this range.
          ENDIF
        ELSE
          WRITE(*,*)'Bad value of DNU_NM in subr. SETUP_NPFMASS: DNU_NM = ', DNU_NM
          STOP
        ENDIF
        NPFMASS(I,0) = ( VNU / XI_KE_H2SO4**3 ) * CONV_V_TO_SO4_H2SO4
        ! IF( WRITE_LOG ) WRITE(AUNIT1,'(I6,2F25.6,3D15.6)') I, XI_KE_H2SO4, XI_KE_SULFATE, NPFMASS(I,0:2)
      ENDDO

90    FORMAT(/,' RH[%]','   R_ambient/R_dry[H2SO4]',' R_ambient/R_dry[SULFATE]',
     &         '   NPFMASS(I,0:2)[ugSO4]')
      RETURN
      END SUBROUTINE SETUP_NPFMASS


      SUBROUTINE NUCLEATION_RATE (T,RH,NA,MB,Z,JTOT,JGAS,JION)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8), INTENT( IN ) :: T      ! temperature [K]
      REAL(8), INTENT( IN ) :: RH     ! relative humidity [%]
      REAL(8), INTENT( IN ) :: NA     ! H2SO4 concentration [molecules/cm^3]
      REAL(8), INTENT( IN ) :: MB     ! NH3 concentration [ppt]
      REAL(8), INTENT( IN ) :: Z      ! height above the Earth's surface [km]

      ! Output arguments.

      REAL(8) :: JGAS ! homogeneous nucleation rate [#/cm^3/s]
      REAL(8) :: JION ! ion-ion     nucleation rate [#/cm^3/s]
      REAL(8) :: JTOT ! total       nucleation rate [#/cm^3/s]


      SELECT CASE (INUC)
     
      CASE (1)     ! INUC=1: JVM BINARY NUCLEATION SCHEME

        CALL NUCL_JVM      (NA,T,RH,JGAS)

      CASE (2)     ! INUC=2: VEHKAMAKI BINARY NUCLEATION SCHEME

        CALL NUCL_VEHKAMAKI(NA,T,RH,JGAS)

      CASE (3)     ! INUC=3: NAPARI TERNARY NUCLEATION SCHEME

        IF(MB.LT.0.1D+00) THEN
          CALL NUCL_VEHKAMAKI(NA,T,RH,JGAS)
        ELSE
          CALL NUCL_NAPARI(NA,MB,T,RH,JGAS)
        ENDIF
      
      CASE (4)     ! INUC=4: Lines from plots in Eisele and McMurry, 1997.

        CALL NUCL_EISELE_MCMURRY(NA,JGAS)

      END SELECT

      IF ( INCLUDE_ION_ION ) THEN
        CALL NUCL_TURCO (NA,Z,JION)    ! ION-ION RECOMBINATION SCHEME
      ELSE
        JION = 0.0D+00
      ENDIF

      IF( INUC.EQ.1 .OR. INUC.EQ.3 .OR. INUC.EQ.4 ) THEN
        JTOT = MIN( MAX( JGAS+JION, J_LOWER ), J_UPPER )
      ELSEIF( INUC.EQ.2 ) THEN   ! Higher upper limit for Vehkamaki et al. 2002
        JTOT = MIN( MAX( JGAS+JION, J_LOWER ), 1.0D+10 ) 
      ENDIF

      RETURN
      END SUBROUTINE NUCLEATION_RATE


      SUBROUTINE NUCL_NAPARI(NA,MB,T,RH,J)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8) :: NA   ! H2SO4 concentration [molecules/cm^3]
      REAL(8) :: MB   ! NH3 concentration [ppt]
      REAL(8) :: T    ! temperature [K]
      REAL(8) :: RH   ! relative humidity [%]

      ! Output arguments.

      REAL(8) :: J    ! nucleation rate [#/cm^3/s]

      ! Local variables.

      LOGICAL :: VALID_INPUT

      ! Parameters.

      REAL(8), PARAMETER :: PI = 3.141592653589793D+00
      REAL(8), PARAMETER :: AVO = 6.0221367D+23
      REAL(8), PARAMETER :: MWH2SO4 = 98.07948D+00
      REAL(8), PARAMETER :: MWNH3   = 17.03356D+00
      REAL(8), PARAMETER :: MWH2O   = 18.01528D+00

      MB = MIN ( MB, 1.00D+02 )   ! Cap at the maximum value for which the parameterization is valid.
      NA = MIN ( NA, 1.00D+09 )   ! Cap at the maximum value for which the parameterization is valid.

      ! Check the conditions of validity for input parameters.

      VALID_INPUT = .TRUE.
      IF     ( T  .LT. 240.00D+00  .OR.  T  .GT. 300.00D+00 ) THEN
        VALID_INPUT = .FALSE.
      ELSEIF ( RH .LT.   0.50D+00  .OR.  RH .GT.  95.00D+00 ) THEN
        VALID_INPUT = .FALSE.
      ELSEIF ( NA .LT.   1.00D+04  .OR.  NA .GT.   1.00D+09 ) THEN    ! Upper limit should have no effect.
        VALID_INPUT = .FALSE.
      ELSEIF ( MB .LT.   1.00D-01  .OR.  MB .GT.   1.00D+02 ) THEN    ! Upper limit should have no effect.
        VALID_INPUT = .FALSE.
      ENDIF

      IF ( .NOT. VALID_INPUT ) THEN
        J = J_LOWER
        RETURN
      ENDIF

      J = J_NAPARI(NA,MB,T,RH)
      IF ( J .LT. 1.0D-05 ) THEN      
        J = J_LOWER
        RETURN
      ELSEIF (  J .GT. 1.0D+06 ) THEN  
        J = J_UPPER
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE NUCL_NAPARI


      REAL(8) FUNCTION NH2SO4_NAPARI(J,T)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: J    ! nucleation rate [#/cm^3/s]
      REAL(8) :: T    ! K
  
      ! Scratch local variables.

      REAL(8) :: LNJ

      LNJ = LOG ( J )

      NH2SO4_NAPARI = 38.1645D+00 + 0.77410D+00*LNJ + 
     &                 0.00298879D+00*LNJ*LNJ-0.357605D+00*T - 
     &                 0.00366358D+00*T*LNJ + 0.0008553D+00*T*T
      NH2SO4_NAPARI = MAX (NH2SO4_NAPARI, 1.0D-30)

      RETURN
      END FUNCTION NH2SO4_NAPARI


      REAL(8) FUNCTION NNH3_NAPARI(J,T)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: J    ! nucleation rate [#/cm^3/s]
      REAL(8) :: T    ! K
  
      ! Scratch local variables.

      REAL(8) :: LNJ

      LNJ = LOG ( J )

      NNH3_NAPARI =   26.8982D+00 + 0.682905D+00*LNJ + 
     &                0.0035752D+00*LNJ*LNJ -0.265748D+00*T - 
     &                0.00341895D+00*T*LNJ + 0.000673454D+00*T*T
      NNH3_NAPARI = MAX (NNH3_NAPARI, 1.0D-30)

      RETURN
      END FUNCTION NNH3_NAPARI


      REAL(8) FUNCTION NTOT_NAPARI(J,T)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: T    ! [K]
      REAL(8) :: J    ! nucleation rate [#/cm^3/s]

      ! Scratch local variables.

      REAL(8),PARAMETER :: A = 79.3484D+00        ! Coefficients
      REAL(8),PARAMETER :: B = 1.7384D+00         ! Coefficients
      REAL(8),PARAMETER :: C = 0.00711403D+00     ! Coefficients
      REAL(8),PARAMETER :: D = -0.74493D+00       ! Coefficients
      REAL(8),PARAMETER :: E = -0.008202608D+00   ! Coefficients
      REAL(8),PARAMETER :: F = 0.0017855D+00      ! Coefficients
      REAL(8) :: LNJ

      LNJ  = LOG ( J )

      NTOT_NAPARI = A + B*LNJ + C*LNJ*LNJ + D*T + E*T*LNJ + F*T*T 
      NTOT_NAPARI = MAX (NTOT_NAPARI, 1.0D-30)

      RETURN
      END FUNCTION NTOT_NAPARI


      REAL(8) FUNCTION RSTAR_NAPARI(J,T)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: T    ! temperature [K]
      REAL(8) :: J    ! nucleation rate [#/cm^3/s]

      ! Scratch local variables.
      REAL(8) :: LNJ    

      LNJ = LOG(J)
      RSTAR_NAPARI = 0.141027D+00 - 0.00122625D+00*LNJ - 
     &               7.82211D-06*LNJ*LNJ - 0.00156727D+00*T - 
     &               0.00003076D+00*T*LNJ + 0.0000108375D+00*T*T

      RETURN
      END FUNCTION RSTAR_NAPARI


      REAL(8) FUNCTION J_NAPARI(NA,MB,T,RH)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: NA   ! [molecules/cm^3]
      REAL(8) :: T    ! [K]
      REAL(8) :: RH   ! [%]
      REAL(8) :: MB   ! [ppt]

      ! Scratch local variables.

      REAL(8) :: X,Y,EX,RY,Z,W,LNJ
      REAL(8), DIMENSION(20) :: F
      INTEGER :: J

      ! Parameters.

      REAL(8), DIMENSION(20,0:3) :: A

      DATA A( 1,0:3)/  -0.355297000D+00,  -0.338449000D+02,   0.345360000D+00,  -0.824007000D-03/
      DATA A( 2,0:3)/   0.313735000D+01,  -0.772861000D+00,   0.561204000D-02,  -0.974576000D-05/
      DATA A( 3,0:3)/   0.190359000D+02,  -0.170957000D+00,   0.479808000D-03,  -0.414699000D-06/
      DATA A( 4,0:3)/   0.107605000D+01,   0.148932000D+01,  -0.796052000D-02,   0.761229000D-05/
      DATA A( 5,0:3)/   0.609160000D+01,  -0.125378000D+01,   0.939836000D-02,  -0.174927000D-04/
      DATA A( 6,0:3)/   0.311760000D+00,   0.164009000D+01,  -0.343852000D-02,  -0.109753000D-04/
      DATA A( 7,0:3)/  -0.200738000D-01,  -0.752115000D+00,   0.525813000D-02,  -0.898038000D-05/
      DATA A( 8,0:3)/   0.165536000D+00,   0.326623000D+01,  -0.489703000D-01,   0.146967000D-03/
      DATA A( 9,0:3)/   0.652645000D+01,  -0.258002000D+00,   0.143456000D-02,  -0.202036000D-05/
      DATA A(10,0:3)/   0.368024000D+01,  -0.204098000D+00,   0.106259000D-02,  -0.126560000D-05/
      DATA A(11,0:3)/  -0.665140000D-01,  -0.782382000D+01,   0.122938000D-01,   0.618554000D-04/
      DATA A(12,0:3)/   0.658740000D+00,   0.190542000D+00,  -0.165718000D-02,   0.341744000D-05/
      DATA A(13,0:3)/   0.599321000D-01,   0.596475000D+01,  -0.362432000D-01,   0.493370000D-04/
      DATA A(14,0:3)/  -0.732731000D+00,  -0.184179000D-01,   0.147186000D-03,  -0.237711000D-06/
      DATA A(15,0:3)/   0.728429000D+00,   0.364736000D+01,  -0.274220000D-01,   0.493478000D-04/
      DATA A(16,0:3)/   0.413016000D+02,  -0.357520000D+00,   0.904383000D-03,  -0.573788000D-06/
      DATA A(17,0:3)/  -0.160336000D+00,   0.889881000D-02,  -0.539514000D-04,   0.839522000D-07/
      DATA A(18,0:3)/   0.857868000D+01,  -0.112358000D+00,   0.472626000D-03,  -0.648365000D-06/
      DATA A(19,0:3)/   0.530167000D-01,  -0.198815000D+01,   0.157827000D-01,  -0.293564000D-04/
      DATA A(20,0:3)/  -0.232736000D+01,   0.234646000D-01,  -0.765190000D-04,   0.804590000D-07/

      ! Statement function.

      REAL(8) :: Z0,Z1,Z2,Z3,ZT
      Z(Z0,Z1,Z2,Z3,ZT) = Z0 + Z1*ZT + Z2*ZT*ZT + Z3*ZT*ZT*ZT
    
      X  = LOG ( RH * 0.01D+00 )
      EX = RH * 0.01D+00
      Y  = LOG ( NA )
      RY = 1.00D+00/Y
      W  = LOG ( MB ) 

      DO J=1,20
        F(J) = Z(A(J,0),A(J,1),A(J,2),A(J,3),T)
      ENDDO
 
      LNJ= -84.7551D+00 + F(1)*RY + F(2)*Y + F(3)*Y*Y + F(4)*W 
     &                  + F(5)*W*W + F(6)*EX + F(7)*X + F(8)*W*RY
     &                  + F(9)*W*Y + F(10)*EX*Y + F(11)*EX*RY
     &                  + F(12)*EX*W + F(13)*X*RY + F(14)*X*W
     &                  + F(15)*W*W*RY + F(16)*Y*W*W + F(17)*Y*Y*W
     &                  + F(18)*EX*W*W + F(19)*EX*W*RY + F(20)*Y*Y*W*W

      J_NAPARI = EXP ( LNJ )

      RETURN
      END FUNCTION J_NAPARI


      SUBROUTINE NUCL_VEHKAMAKI(NA,T,RH,J)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!     DLW:062005: Created and checked for J value for 82 points in 
!                 Vehkamaki et al. 2002.
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8) :: NA   ! H2SO4 concentration [molecules/cm^3]
      REAL(8) :: T    ! temperature [K]
      REAL(8) :: RH   ! relative humidity [%]

      ! Output arguments.

      REAL(8) :: J    ! nucleation rate [#/cm^3/s]

      ! Local variables.

      REAL(8) :: X    ! mole fraction H2SO4
      REAL(8) :: N    ! number of molecules in the critical nucleus
      LOGICAL :: VALID_INPUT

      ! Parameters.

      REAL(8), PARAMETER :: PI = 3.141592653589793D+00
      REAL(8), PARAMETER :: AVO = 6.0221367D+23
      REAL(8), PARAMETER :: MWH2SO4 = 98.07948D+00
      REAL(8), PARAMETER :: MWH2O   = 18.01528D+00


      ! Check for conditions of validity for input parameters.
 
      NA = MIN ( NA, 1.00D+11 )  ! Cap at the maximum value valid for the parameterization. 

      VALID_INPUT = .TRUE.
      IF     ( T  .LT. 190.00D+00  .OR.  T  .GT. 305.15D+00 ) THEN
        VALID_INPUT = .FALSE.
      ELSEIF ( RH .LT.   0.01D+00  .OR.  RH .GT. 100.000001D+00 ) THEN
        VALID_INPUT = .FALSE.
      ELSEIF ( NA .LT.   1.00D+04  .OR.  NA .GT.   1.00D+11 ) THEN     ! Upper limit should have no effect. 
        VALID_INPUT = .FALSE.
      ENDIF
      IF ( .NOT. VALID_INPUT ) THEN
        J = J_LOWER
        RETURN
      ENDIF

      X = XSTAR_VEHKAMAKI(NA,T,RH)

      J = J_VEHKAMAKI(NA,T,RH,X)
      IF ( J .LT. 1.0D-07 ) THEN 
        J = J_LOWER
        RETURN
      ELSEIF ( J .GT. 1.0D+10 ) THEN 
        J = 1.0D+10                     ! Maximum value valid for this parameterization.
        RETURN
      ENDIF

      ! Properties of the critical nucleus.

      N = NTOT_VEHKAMAKI(NA,T,RH,X)
      IF ( N .LT. 4.0D+00 ) THEN        ! Check for condition on NTOT.
        J = J_LOWER
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE NUCL_VEHKAMAKI


      REAL(8) FUNCTION XSTAR_VEHKAMAKI(NA,T,RH)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: NA   ! molecules/cm^3
      REAL(8) :: T    ! K
      REAL(8) :: RH   ! %

      ! Scratch local variables.

      REAL(8) :: X

      X = LOG ( RH * 0.01D+00 )

      XSTAR_VEHKAMAKI =   0.740997D+00       - 0.00266379D+00   *T
     &       + ( 0.0000504022D+00*T - 0.00349998D+00      ) * LOG ( NA )
     &       + ( 0.00201048D+00     - 0.000183289D+00  *T ) * X
     &       + ( 0.00157407D+00     - 0.0000179059D+00 *T ) * X*X
     &       + ( 0.000184403D+00    - 1.50345D-06      *T ) * X**3

      RETURN
      END FUNCTION XSTAR_VEHKAMAKI


      REAL(8) FUNCTION NTOT_VEHKAMAKI(NA,T,RH,XS)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: NA   ! [molecules/cm^3]
      REAL(8) :: T    ! [K]
      REAL(8) :: RH   ! [%]
      REAL(8) :: XS   ! mole fraction H2SO4 in the critical nucleus

      ! Scratch local variables.

      REAL(8) :: A,B,C,D,E,F,G,H,I,J    ! Coefficients
      REAL(8) :: X,Y,RX

      ! Parameters.

      REAL(8), PARAMETER :: A1 = - 0.00295413D+00
      REAL(8), PARAMETER :: A2 = - 0.0976834D+00
      REAL(8), PARAMETER :: A3 = + 0.00102485D+00
      REAL(8), PARAMETER :: A4 = - 2.18646D-06
      REAL(8), PARAMETER :: A5 = - 0.101717D+00

      REAL(8), PARAMETER :: B1 = - 0.00205064D+00
      REAL(8), PARAMETER :: B2 = - 0.00758504D+00
      REAL(8), PARAMETER :: B3 = + 0.000192654D+00
      REAL(8), PARAMETER :: B4 = - 6.7043D-07 
      REAL(8), PARAMETER :: B5 = - 0.255774D+00

      REAL(8), PARAMETER :: C1 = + 0.00322308D+00
      REAL(8), PARAMETER :: C2 = + 0.000852637D+00
      REAL(8), PARAMETER :: C3 = - 0.0000154757D+00
      REAL(8), PARAMETER :: C4 = + 5.66661D-08
      REAL(8), PARAMETER :: C5 = + 0.0338444D+00

      REAL(8), PARAMETER :: D1 = + 0.0474323D+00 
      REAL(8), PARAMETER :: D2 = - 0.000625104D+00
      REAL(8), PARAMETER :: D3 = + 2.65066D-06 
      REAL(8), PARAMETER :: D4 = - 3.67471D-09
      REAL(8), PARAMETER :: D5 = - 0.000267251D+00

      REAL(8), PARAMETER :: E1 = - 0.0125211D+00
      REAL(8), PARAMETER :: E2 = + 0.00580655D+00
      REAL(8), PARAMETER :: E3 = - 0.000101674D+00
      REAL(8), PARAMETER :: E4 = + 2.88195D-07
      REAL(8), PARAMETER :: E5 = + 0.0942243D+00

      REAL(8), PARAMETER :: F1 = - 0.038546D+00
      REAL(8), PARAMETER :: F2 = - 0.000672316D+00
      REAL(8), PARAMETER :: F3 = + 2.60288D-06    
      REAL(8), PARAMETER :: F4 = + 1.19416D-08
      REAL(8), PARAMETER :: F5 = - 0.00851515D+00

      REAL(8), PARAMETER :: G1 = - 0.0183749D+00
      REAL(8), PARAMETER :: G2 = + 0.000172072D+00
      REAL(8), PARAMETER :: G3 = - 3.71766D-07    
      REAL(8), PARAMETER :: G4 = - 5.14875D-10
      REAL(8), PARAMETER :: G5 = + 0.00026866D+00

      REAL(8), PARAMETER :: H1 = - 0.0619974D+00
      REAL(8), PARAMETER :: H2 = + 0.000906958D+00
      REAL(8), PARAMETER :: H3 = - 9.11728D-07    
      REAL(8), PARAMETER :: H4 = - 5.36796D-09
      REAL(8), PARAMETER :: H5 = - 0.00774234D+00

      REAL(8), PARAMETER :: I1 = + 0.0121827D+00
      REAL(8), PARAMETER :: I2 = - 0.00010655D+00
      REAL(8), PARAMETER :: I3 = + 2.5346D-07      
      REAL(8), PARAMETER :: I4 = - 3.63519D-10 
      REAL(8), PARAMETER :: I5 = + 0.000610065D+00

      REAL(8), PARAMETER :: J1 = + 0.000320184D+00
      REAL(8), PARAMETER :: J2 = - 0.0000174762D+00 
      REAL(8), PARAMETER :: J3 = + 6.06504D-08      
      REAL(8), PARAMETER :: J4 = - 1.42177D-11
      REAL(8), PARAMETER :: J5 = + 0.000135751D+00

      ! Statement function.

      REAL(8) :: Z
      REAL(8) :: Z1,Z2,Z3,Z4,Z5,ZT,ZX
      Z(Z1,Z2,Z3,Z4,Z5,ZT,ZX) = Z1 + Z2*ZT + Z3*ZT*ZT + Z4*ZT**3 + Z5*ZX

      RX = 1.0D+00/XS
      X  = LOG ( RH * 0.01D+00 )
      Y  = LOG ( NA )

      A = Z(A1,A2,A3,A4,A5,T,RX)
      B = Z(B1,B2,B3,B4,B5,T,RX)
      C = Z(C1,C2,C3,C4,C5,T,RX)
      D = Z(D1,D2,D3,D4,D5,T,RX)
      E = Z(E1,E2,E3,E4,E5,T,RX)
      F = Z(F1,F2,F3,F4,F5,T,RX)
      G = Z(G1,G2,G3,G4,G5,T,RX)
      H = Z(H1,H2,H3,H4,H5,T,RX)
      I = Z(I1,I2,I3,I4,I5,T,RX)
      J = Z(J1,J2,J3,J4,J5,T,RX)

      NTOT_VEHKAMAKI = EXP ( A + B*X + C*X*X + D*X**3 + E*Y + F*X*Y + 
     &                       G*X*X*Y + H*Y*Y + I*X*Y*Y + J*Y**3 )

      RETURN
      END FUNCTION NTOT_VEHKAMAKI


      REAL(8) FUNCTION RSTAR_VEHKAMAKI(XS,NTOT)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: NTOT  ! # molecules in the critical nucleus
      REAL(8) :: XS    ! mole fraction H2SO4

      RSTAR_VEHKAMAKI = EXP( -1.6524245D+00 + 0.42316402D+00*XS
     &                                     + 0.3346648D+00*LOG( NTOT ) )

      RETURN
      END FUNCTION RSTAR_VEHKAMAKI


      REAL(8) FUNCTION J_VEHKAMAKI(NA,T,RH,XS)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: NA   ! [molecules/cm^3]
      REAL(8) :: T    ! [K]
      REAL(8) :: RH   ! [%]
      REAL(8) :: XS   ! mole fraction H2SO4 in the critical nucleus

      ! Scratch local variables.

      REAL(8) :: A,B,C,D,E,F,G,H,I,J    ! Coefficients
      REAL(8) :: X,Y,RX

      ! Parameters.

      REAL(8), PARAMETER :: A1 = + 0.14309D+00
      REAL(8), PARAMETER :: A2 = + 2.21956D+00
      REAL(8), PARAMETER :: A3 = - 0.0273911D+00
      REAL(8), PARAMETER :: A4 = + 0.0000722811D+00
      REAL(8), PARAMETER :: A5 = + 5.91822D+00

      REAL(8), PARAMETER :: B1 = + 0.117489D+00
      REAL(8), PARAMETER :: B2 = + 0.462532D+00
      REAL(8), PARAMETER :: B3 = - 0.0118059D+00
      REAL(8), PARAMETER :: B4 = + 0.0000404196D+00
      REAL(8), PARAMETER :: B5 = + 15.7963D+00

      REAL(8), PARAMETER :: C1 = - 0.215554D+00
      REAL(8), PARAMETER :: C2 = - 0.0810269D+00
      REAL(8), PARAMETER :: C3 = + 0.00143581D+00
      REAL(8), PARAMETER :: C4 = - 4.7758D-06
      REAL(8), PARAMETER :: C5 = - 2.91297D+00

      REAL(8), PARAMETER :: D1 = - 3.58856D+00 
      REAL(8), PARAMETER :: D2 = + 0.049508D+00
      REAL(8), PARAMETER :: D3 = - 0.00021382D+00
      REAL(8), PARAMETER :: D4 = + 3.10801D-07
      REAL(8), PARAMETER :: D5 = - 0.0293333D+00

      REAL(8), PARAMETER :: E1 = + 1.14598D+00
      REAL(8), PARAMETER :: E2 = - 0.600796D+00
      REAL(8), PARAMETER :: E3 = + 0.00864245D+00  
      REAL(8), PARAMETER :: E4 = - 0.0000228947D+00
      REAL(8), PARAMETER :: E5 = - 8.44985D+00  

      REAL(8), PARAMETER :: F1 = + 2.15855D+00
      REAL(8), PARAMETER :: F2 = + 0.0808121D+00
      REAL(8), PARAMETER :: F3 = - 0.000407382D+00
      REAL(8), PARAMETER :: F4 = - 4.01957D-07
      REAL(8), PARAMETER :: F5 = + 0.721326D+00

      REAL(8), PARAMETER :: G1 = + 1.6241D+00
      REAL(8), PARAMETER :: G2 = - 0.0160106D+00
      REAL(8), PARAMETER :: G3 = + 0.0000377124D+00
      REAL(8), PARAMETER :: G4 = + 3.21794D-08
      REAL(8), PARAMETER :: G5 = - 0.0113255D+00

      REAL(8), PARAMETER :: H1 = + 9.71682D+00
      REAL(8), PARAMETER :: H2 = - 0.115048D+00
      REAL(8), PARAMETER :: H3 = + 0.000157098D+00
      REAL(8), PARAMETER :: H4 = + 4.00914D-07 
      REAL(8), PARAMETER :: H5 = + 0.71186D+00

      REAL(8), PARAMETER :: I1 = - 1.05611D+00
      REAL(8), PARAMETER :: I2 = + 0.00903378D+00
      REAL(8), PARAMETER :: I3 = - 0.0000198417D+00
      REAL(8), PARAMETER :: I4 = + 2.46048D-08 
      REAL(8), PARAMETER :: I5 = - 0.0579087D+00

      REAL(8), PARAMETER :: J1 = - 0.148712D+00
      REAL(8), PARAMETER :: J2 = + 0.00283508D+00 
      REAL(8), PARAMETER :: J3 = - 9.24619D-06     
      REAL(8), PARAMETER :: J4 = + 5.00427D-09
      REAL(8), PARAMETER :: J5 = - 0.0127081D+00   

      ! Statement function.

      REAL(8) :: Z
      REAL(8) :: Z1,Z2,Z3,Z4,Z5,ZT,ZX
      Z(Z1,Z2,Z3,Z4,Z5,ZT,ZX) = Z1 + Z2*ZT + Z3*ZT*ZT + Z4*ZT**3 + Z5*ZX

      RX = 1.0D+00/XS
      X  = LOG ( RH * 0.01D+00 )
      Y  = LOG ( NA )

      A = Z(A1,A2,A3,A4,A5,T,RX)
      B = Z(B1,B2,B3,B4,B5,T,RX)
      C = Z(C1,C2,C3,C4,C5,T,RX)
      D = Z(D1,D2,D3,D4,D5,T,RX)
      E = Z(E1,E2,E3,E4,E5,T,RX)
      F = Z(F1,F2,F3,F4,F5,T,RX)
      G = Z(G1,G2,G3,G4,G5,T,RX)
      H = Z(H1,H2,H3,H4,H5,T,RX)
      I = Z(I1,I2,I3,I4,I5,T,RX)
      J = Z(J1,J2,J3,J4,J5,T,RX)

      J_VEHKAMAKI = EXP ( A + B*X + C*X*X + D*X**3 + E*Y + F*X*Y + 
     &                    G*X*X*Y + H*Y*Y + I*X*Y*Y + J*Y**3 )

      RETURN
      END FUNCTION J_VEHKAMAKI


      SUBROUTINE NUCL_JVM(NA,T,RH,J)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8) :: NA   ! H2SO4 concentration [molecules/cm^3]
      REAL(8) :: T    ! temperature [K]
      REAL(8) :: RH   ! relative humidity [%]

      ! Output arguments.

      REAL(8) :: J    ! nucleation rate [#/cm^3/s]

      ! Local variables.

      REAL(8), DIMENSION(0:2) :: B    ! coefficient of fitted curve function
      LOGICAL :: VALID_INPUT
      REAL(8) :: LOG10J, XX, XLOG10NA

      NA = MIN ( NA, 1.00D+14 ) 

      ! Check for conditions of validity for input parameters.

      CALL JVM_COEF (T,RH,B)

      XX = -B(1) / (2.0D+00 * B(2))
      XLOG10NA = LOG10(NA)

      VALID_INPUT = .TRUE.
      IF(XLOG10NA .GT. XX) THEN
        VALID_INPUT = .FALSE.
      ELSEIF ( NA .LT. 1.00D+04  .OR.  NA .GT. 1.00D+14 ) THEN     ! Upper limit should have no effect. 
        VALID_INPUT = .FALSE.
      ENDIF

      IF ( .NOT. VALID_INPUT ) THEN
        J = J_LOWER
        RETURN
      ENDIF

      ! Calculate the nucleation rate.

      LOG10J = B(0) + B(1)*XLOG10NA + B(2)*XLOG10NA*XLOG10NA

      J = 10.0**( LOG10J )

      IF ( J .LT. 1.0D-03 ) THEN 
        J = J_LOWER 
        RETURN
      ELSEIF (  J .GT. 1.0D+05 ) THEN 
        J = J_UPPER  
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE NUCL_JVM


      SUBROUTINE JVM_COEF(T,RH,B)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input Arguments.

      REAL(8) :: T    ! [K]
      REAL(8) :: RH   ! [%]

      ! Output Arguments.

      REAL(8), DIMENSION(0:2) :: B

      ! Scratch local variables.

      REAL(8) :: X,Y,Z
      REAL(8), DIMENSION(0:4) :: BB0,BB1,BB2

      ! Parameters.

      REAL(8), DIMENSION(0:4,0:4) :: C0,C1,C2 
      INTEGER :: J

      DATA C0(0,0:4)/  -0.493166930D+03,  -0.746060630D+03,   0.427001110D+04,   0.403205190D+04,  -0.103305460D+05/
      DATA C0(1,0:4)/   0.227147650D+04,   0.296173760D+04,  -0.345957800D+05,  -0.435346380D+05,   0.583845900D+05/
      DATA C0(2,0:4)/  -0.596070430D+04,  -0.921088590D+04,   0.919397290D+05,   0.131516870D+06,  -0.130410870D+06/
      DATA C0(3,0:4)/   0.673050820D+04,   0.112160080D+05,  -0.101584900D+06,  -0.151340760D+06,   0.133574610D+06/
      DATA C0(4,0:4)/  -0.261458200D+04,  -0.444305610D+04,   0.394068830D+05,   0.594771880D+05,  -0.496408160D+05/
      DATA C1(0,0:4)/   0.100665720D+03,   0.108313420D+03,  -0.113683880D+04,  -0.888430520D+03,   0.302914500D+04/
      DATA C1(1,0:4)/  -0.504046310D+03,  -0.474794750D+03,   0.904264540D+04,   0.909745930D+04,  -0.201185820D+05/
      DATA C1(2,0:4)/   0.137492990D+04,   0.161927870D+04,  -0.243558850D+05,  -0.269705410D+05,   0.498079760D+05/
      DATA C1(3,0:4)/  -0.156449310D+04,  -0.200427440D+04,   0.270531380D+05,   0.307133920D+05,  -0.531638600D+05/
      DATA C1(4,0:4)/   0.604879820D+03,   0.788620250D+03,  -0.104754690D+05,  -0.119805370D+05,   0.200426110D+05/
      DATA C2(0,0:4)/  -0.516074080D+01,  -0.398706400D+01,   0.720177880D+02,   0.509391250D+02,  -0.211700440D+03/
      DATA C2(1,0:4)/   0.284462970D+02,   0.212336580D+02,  -0.573120500D+03,  -0.497927540D+03,   0.152757040D+04/
      DATA C2(2,0:4)/  -0.795352710D+02,  -0.764306230D+02,   0.155886490D+04,   0.143625300D+04,  -0.392583270D+04/
      DATA C2(3,0:4)/   0.904137910D+02,   0.939499970D+02,  -0.173349330D+04,  -0.160777790D+04,   0.422689150D+04/
      DATA C2(4,0:4)/  -0.345602640D+02,  -0.362758330D+02,   0.667993290D+03,   0.619425480D+03,  -0.159032770D+04/

      ! Statement function.

      REAL(8) :: Z0,Z1,Z2,Z3,Z4,V
      Z(Z0,Z1,Z2,Z3,Z4,V) = Z0+Z1*V+Z2*V*V+Z3*V*V*V+Z4*V*V*V*V

      X  = ( T - 273.15D+00 ) * 0.01D+00
      Y  = RH * 0.01D+00

      DO J = 0, 4
        BB0(J) = Z(C0(J,0),C0(J,1),C0(J,2),C0(J,3),C0(J,4),X)
      ENDDO
      DO J = 0, 4
        BB1(J) = Z(C1(J,0),C1(J,1),C1(J,2),C1(J,3),C1(J,4),X)
      ENDDO
      DO J = 0, 4
        BB2(J) = Z(C2(J,0),C2(J,1),C2(J,2),C2(J,3),C2(J,4),X)
      ENDDO

      B(0) = Z(BB0(0),BB0(1),BB0(2),BB0(3),BB0(4),Y)
      B(1) = Z(BB1(0),BB1(1),BB1(2),BB1(3),BB1(4),Y)
      B(2) = Z(BB2(0),BB2(1),BB2(2),BB2(3),BB2(4),Y)

      RETURN
      END SUBROUTINE JVM_COEF


      SUBROUTINE NUCL_TURCO(NA,Z,J)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8) :: NA   ! H2SO4 concentration [molecules/cm^3]
      REAL(8) :: Z    ! geographical height [km]

      ! Output arguments.

      REAL(8) :: J    ! nucleation rate [#/cm^3/s]

      J = J_TURCO(NA,Z)

      RETURN
      END SUBROUTINE NUCL_TURCO


      REAL(8) FUNCTION J_TURCO(NA,Z)
!-----------------------------------------------------------------------------------------------------------------------
!     LSC/DLW 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: NA   ! H2SO4 concentration [molecules/cm^3]
      REAL(8) :: Z    ! geographical height [Km]

      ! Local variable
      REAL(8) :: Q    ! local ionization rate [/cm^3/s]

      ! Parameters.

      REAL(8), PARAMETER :: NA0 = 5.00D+06
      REAL(8), PARAMETER :: F0 = 1.00D-03
      INTEGER, PARAMETER :: NSTAR = 3

 
      IF(Z .LT. 2.25) THEN
        Q = 2.00D+00
      ELSEIF(Z .GE. 2.25D+00 .AND. Z .LT. 11.00D+00) THEN
        Q = 3.10D+00 * (Z - 2.25D+00) + 2.00D+00
      ELSEIF ( Z .GE. 11.00D+00) THEN
        Q = 29.00D+00 + (Z - 11.00D+00)
      ENDIF

      J_TURCO = Q * F0 * ( NA / NA0 )**NSTAR
      J_TURCO = MIN( Q, J_TURCO )

      RETURN
      END FUNCTION J_TURCO


      SUBROUTINE NUCL_EISELE_MCMURRY(NA,J)
!-----------------------------------------------------------------------------------------------------------------------
!     DLW: 7-27-06: These expressions were derived from Figure 7 of 
!          Eisele, F. L., and McMurry, P. H. (1997). 
!          Recent progress in understanding particle nucleation and growth,
!          Phil. Trans. R. Soc. Lond. B, 352, 191-201.
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8) :: NA   ! H2SO4 concentration [molecules/cm^3]

      ! Output arguments.

      REAL(8) :: J    ! nucleation rate [#/cm^3/s]

      ! Parameters.

      REAL(8), PARAMETER :: K1 = 5.8D-13
      REAL(8), PARAMETER :: K2 = 3.5D-15
      REAL(8), PARAMETER :: K3 = 3.7D-14
      LOGICAL, PARAMETER :: LOWER_CURVE = .FALSE.
      LOGICAL, PARAMETER :: UPPER_CURVE = .FALSE.
      LOGICAL, PARAMETER :: AVGUL_CURVE = .TRUE.

      ! Local variables.

      ! INTEGER :: I
      ! REAL(8) :: S

      IF(LOWER_CURVE) J = K1*NA
      IF(UPPER_CURVE) J = K2*NA*NA
      IF(AVGUL_CURVE) J = K3*NA**1.5
       
!-----------------------------------------------------------------------------------------------------------------------
!     Plot the expressions.
!-----------------------------------------------------------------------------------------------------------------------
!     DO I=1,5
!       S = 10.0**(4+I-1)
!       WRITE(78,91 )S, K1*S, K3*S**1.5, K2*S*S
!     ENDDO
!-----------------------------------------------------------------------------------------------------------------------

91    FORMAT(4D15.5)
      RETURN
      END SUBROUTINE NUCL_EISELE_MCMURRY


      SUBROUTINE F_KK02(PRS,TK,RHPERCENT,NAPERCM3,NH3PPT,KC,DNPF_NM,DNUC_NM,JP_TO_J)
!-----------------------------------------------------------------------------------------------------------------------
!     DLW, 8-1-06.
!     Routine to calculate the ratio of the new particle formation rate at
!     a user-specified diameter to the nucleation rate.
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8) :: PRS       ! ambient pressure               [Pa]
      REAL(8) :: TK        ! ambient temperature            [K]
      REAL(8) :: RHPERCENT ! relative humidity              [%]
      REAL(8) :: NAPERCM3  ! sulfuric acid concentration    [#/cm^3]
      REAL(8) :: NH3PPT    ! ammonia concentration          [ppt]         
      REAL(8) :: KC        ! condensational sink            [1/s]         
      REAL(8) :: DNPF_NM   ! user-selected diameter for new particles    [nm]         
      REAL(8) :: DNUC_NM   ! initial diameter before growth to size DNPF [nm]         

      ! Output arguments.

      REAL(8) :: JP_TO_J   ! ratio of new particle formation rate at diameter DNP
                           ! to the nucleation rate. [1]

      ! Local scratch variables.

      REAL(8) :: M_EFFECTIVE  ! effective molar mass of H2SO4 including the water
                              ! and ammonia that instantaneously equilibrate with
                              ! the particle [g/mol]
      REAL(8) :: GR_INORG     ! growth rate in the free-molecular regime due to 
                              ! inorganic species [nm/h]
      REAL(8) :: X_NH3        ! # of NH3 molecules condensing per H2SO4 molecule
                              ! - lies in the range 0-2 [1]
      REAL(8) :: X_H2O        ! # of H2O molecules condensing per H2SO4 molecule [1]
      REAL(8) :: GAMMA        ! Eq.(22) of KK2002. [nm^2 m^2 h^-1]
      REAL(8) :: DMEAN_NM     ! number mean diameter over all modes [nm]
      REAL(8) :: CS_PRIME     ! condensational sink espressed as  s [m^-2]
      REAL(8) :: ETA          ! ETA parameter of the KK2002 model   [nm]

      !-----------------------------------------------------------------------------------------------------------------
      ! The growth rate due to H2SO4+H2O+NH3 condensation can be expressed as
      !
      !     GR(nm/h) = GR_CONST * T^0.5 * Meff(g/mol) * C(#/cm^3) 
      !
      ! where Meff is the effective molar mass of H2SO4 with water and ammonia
      ! instantaneously equilibrating to the particle. See notes of 7-19-05.
      ! Free-molecular growth is assumed.
      ! 
      ! ALPHA is the mass accommodation coefficient[1]. See Eq.(20) of KK02:
      !   ALPHA should be unity here, even if less than that elsewhere.  
      ! 
      ! DSTAR_DENSITY is the density of a critical nucleus [g/cm^3].
      !
      ! The mean thermal speed of a H2SO4 molecule [m/s] is CMEAN_H2SO4 * T^0.5,
      !   with T the absolute temperature [K].
      !
      ! AVO is Avogadro's number [#/mol]. 
      !
      ! The factor 3.6D+12 converts [m/s] to [nm/h].
      !-----------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: ALPHA = 1.0D+00                           ! [1]
      REAL(8), PARAMETER :: DSTAR_DENSITY = 1.6D+00                   ! [g/cm^3]
      REAL(8), PARAMETER :: CMEAN_H2SO4 = 14.692171525317413D+00      ! [m/s/K^0.5]
      REAL(8), PARAMETER :: GR_CONST = 0.5D+00 * ALPHA * CMEAN_H2SO4 * 3.6D+12
     &                               / ( AVO * DSTAR_DENSITY )
      REAL(8), PARAMETER :: FOURPI = 4.0D+00 * PI 

      !-----------------------------------------------------------------------------------------------------------------
      ! The proportionality factor GAMMA for eq.22 of KK2002 have this 
      ! constant extracted.
      !-----------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: GAMMA_PRIME = 11.79348270D+00


      !-----------------------------------------------------------------------------------------------------------------
      ! Get the number mean diameter over all modes [nm].
      !-----------------------------------------------------------------------------------------------------------------
      DMEAN_NM = 1.0D+09 * AVG_DP_OF_AVG_MASS_METERS   ! Stored in aero_param.f.

      !-----------------------------------------------------------------------------------------------------------------
      ! Get GAMMA from Eq.(22) in KK2002 [nm^2 m^2 /h].
      !
      ! DSTAR_DENSITY is in [g/cm^3], with the conversion from [kg/m^3] in
      ! KK02 Eq.(22) folded into the parameter GAMMA_PRIME.
      !-----------------------------------------------------------------------------------------------------------------
      GAMMA = GAMMA_PRIME * DNUC_NM**0.2 * DNPF_NM**0.075 * DMEAN_NM**0.048 
     &                    * DSTAR_DENSITY**(-0.33)        * TK**(-0.75)

      !-----------------------------------------------------------------------------------------------------------------
      ! Get the effective molecular weight of condensing H2SO4 [g/mol].
      !-----------------------------------------------------------------------------------------------------------------
      CALL EFFECTIVE_MW(PRS,TK,RHPERCENT,NAPERCM3,NH3PPT,X_NH3,X_H2O)
      M_EFFECTIVE = MW_H2SO4 + X_NH3*MW_NH3 + X_H2O*MW_H2O

      !-----------------------------------------------------------------------------------------------------------------
      ! Get the growth rate due to condensation of inorganics [nm/h].
      !-----------------------------------------------------------------------------------------------------------------
      GR_INORG = GR_CONST * SQRT(TK) * M_EFFECTIVE * NAPERCM3 

      !-----------------------------------------------------------------------------------------------------------------
      ! Get the condensational sink (CS') in [m^-2].
      !-----------------------------------------------------------------------------------------------------------------
      CS_PRIME = KC / ( FOURPI * DIFFCOEF_M2S(ILAY) ) 

      !-----------------------------------------------------------------------------------------------------------------
      ! Get ETA of the KK2002 formulation [nm].
      !-----------------------------------------------------------------------------------------------------------------
      ETA = GAMMA * CS_PRIME / MAX( GR_INORG, 1.0D-30 )
      ! WRITE(34,*)'ETA,CS_PRIME,GR_INORG,NAPERCM3 = ',ETA,CS_PRIME,GR_INORG,NAPERCM3 

      IF(ETA .LT. 56.0D+00) THEN
        JP_TO_J = EXP( ETA * ( (1.0D+00/DNPF_NM) - (1.0D+00/DNUC_NM) ) )
      ELSE 
        JP_TO_J = 1.0D-17
      ENDIF
      ! WRITE(34,*)'JP_TO_J = ', JP_TO_J 

      !-----------------------------------------------------------------------------------------------------------------
      ! DLW, 8-1-06: All of these variables were fine.
      !-----------------------------------------------------------------------------------------------------------------
      IF( WRITE_F_KK02 ) THEN
        WRITE(AUNIT1,'(/A/)')   'IN SUBROUTINE F_KK02'
        WRITE(AUNIT1,*)         'DSTAR_DENSITY,TK,DMEAN_NM,GAMMA,DIFFCOEF_M2S(ILAY)'
        WRITE(AUNIT1,'(5D15.5)') DSTAR_DENSITY,TK,DMEAN_NM,GAMMA,DIFFCOEF_M2S(ILAY)
        WRITE(AUNIT1,*)         'X_H2O,X_NH3,M_EFFECTIVE,DNUC_NM,DNPF_NM'
        WRITE(AUNIT1,'(5D15.5)') X_H2O,X_NH3,M_EFFECTIVE,DNUC_NM,DNPF_NM 
        WRITE(AUNIT1,*)         'CS_PRIME,GR_INORG,ETA,JP_TO_J'
        WRITE(AUNIT1,'(5D15.5)') CS_PRIME,GR_INORG,ETA,JP_TO_J 
        WRITE(AUNIT1,'(A)')     '  '
      ENDIF
      !-----------------------------------------------------------------------------------------------------------------

      RETURN
      END SUBROUTINE F_KK02


      SUBROUTINE EFFECTIVE_MW(AIRPRS,TK,RHPERCENT,NAPERCM3,NH3PPT,X_NH3,X_H2O)
!-----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate X_NH3 and X_H2O as defined below.
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8) :: AIRPRS    ! ambient pressure               [Pa]
      REAL(8) :: TK        ! ambient temperature            [K]
      REAL(8) :: RHPERCENT ! relative humidity              [%]
      REAL(8) :: NAPERCM3  ! sulfuric acid concentration    [#/cm^3]
      REAL(8) :: NH3PPT    ! ammonia mixing ratio           [ppt]

      ! Output arguments.

      REAL(8) :: X_NH3     ! # of NH3 molecules condensing per H2SO4 [#]
      REAL(8) :: X_H2O     ! # of H2O molecules condensing per H2SO4 [#]

      ! Parameters.

!     REAL(8), PARAMETER :: DNU_NM   = 10.0D+00        ! diameter of new particles [nm]
      REAL(8), PARAMETER :: RRNU_NM  = 2.0D+00/DNU_NM  ! recipocal radius ... [nm^-1]
      REAL(8), PARAMETER :: SQRT_MWH2SO4_OVER_MWNH3 = 2.400D+00
      REAL(8), PARAMETER :: A0_H2SO4 =  0.300D+00
      REAL(8), PARAMETER :: A1_H2SO4 = -0.608D+00
      REAL(8), PARAMETER :: A2_H2SO4 =  0.701D+00
      REAL(8), PARAMETER :: A3_H2SO4 = -0.392D+00
      REAL(8), PARAMETER :: A0_BISO4 =  0.866D+00
      REAL(8), PARAMETER :: A1_BISO4 = -1.965D+00
      REAL(8), PARAMETER :: A2_BISO4 =  1.897D+00
      REAL(8), PARAMETER :: A3_BISO4 = -0.800D+00
      REAL(8), PARAMETER :: A0_AMSO4 =  0.951D+00
      REAL(8), PARAMETER :: A1_AMSO4 = -2.459D+00
      REAL(8), PARAMETER :: A2_AMSO4 =  2.650D+00
      REAL(8), PARAMETER :: A3_AMSO4 = -1.142D+00
!     LOGICAL, PARAMETER :: INCLUDE_KELVIN_EFFECT = .TRUE.
      LOGICAL, PARAMETER :: INCLUDE_KELVIN_EFFECT = .FALSE.

      ! Local scratch variables.

      REAL(8)       :: NH3PERCM3      ! ammonia number concentration [# cm^-3]
      REAL(8)       :: MFS            ! mole fraction sulfur in the particle:
                                      !   ( mol S ) / ( mol S + mol H2O )
      REAL(8)       :: RHF            ! (fractional RH) * EXP factor for Kelvin effect
      REAL(8)       :: A              ! scratch variable for Kelvin effect

      NAPERCM3 = MAX ( NAPERCM3, 1.0D-30 )
      NH3PPT   = MAX ( NH3PPT,   1.0D-30 )

      !-----------------------------------------------------------------------------------------------------------------
      ! Convert NH3 from ppt to molecules cm^-3.
      !-----------------------------------------------------------------------------------------------------------------
      NH3PERCM3 = 7.25D+04 * NH3PPT * AIRPRS / TK

      !-----------------------------------------------------------------------------------------------------------------
      ! See notes of 7-08-05 for the X_NH3 calculation.
      ! The square root arises from a ratio of the thermal velocities for 
      !   the two molecules. 
      !-----------------------------------------------------------------------------------------------------------------
      X_NH3 = MIN ( SQRT_MWH2SO4_OVER_MWNH3 * NH3PERCM3 / NAPERCM3, 2.0D+00 )

      !-----------------------------------------------------------------------------------------------------------------
      ! Now use the X_NH3 value to classify the new particles as either
      !
      !  (1) acidic            - treat as pure sulfuric acid
      !  (2) half-neutralized  - treat as pure ammonium bisulfate
      !  (3) fully-neutralized - treat as pure ammonium sulfate
      !
      ! Then compute the water uptake for the selected composition.
      !-----------------------------------------------------------------------------------------------------------------
      IF ( X_NH3 .LT. 0.5D+00 ) THEN

        ! Acidic case.

        NPFMASS_REGIME = 0                          ! second index to NPFMASS(:,:)
        IF ( INCLUDE_KELVIN_EFFECT ) THEN
          A = 1.2D+00 - 0.0072D+00*(TK-273.15D+00)
          RHF = 0.01D+00 * RHPERCENT * EXP ( -A * RRNU_NM )
        ELSE
          RHF = 0.01D+00 * RHPERCENT
        ENDIF
        MFS = A0_H2SO4 + A1_H2SO4*RHF + A2_H2SO4*RHF*RHF + A3_H2SO4*RHF**3

      ELSEIF ( X_NH3 .GE. 0.5D+00 .AND. X_NH3 .LE. 1.5D+00 ) THEN

        ! Half-neutralized case.

        NPFMASS_REGIME = 1                          ! second index to NPFMASS(:,:)
        IF ( INCLUDE_KELVIN_EFFECT ) THEN
          RHF = 0.01D+00 * RHPERCENT
          A = ( 1.2D+00 - 0.0072D+00*(TK-273.15D+00) )
     &      * ( 1.0D+00 - 0.17D+00*(1.0D+00-RHF))
     &      * ( 1.0D+00 + 0.95D+00*(1.0D+00-RHF))
          RHF = RHF * EXP ( -A * RRNU_NM )
        ELSE
          RHF = 0.01D+00 * RHPERCENT
        ENDIF
        MFS = A0_BISO4 + A1_BISO4*RHF + A2_BISO4*RHF*RHF + A3_BISO4*RHF**3

      ELSE

        ! Fully-neutralized case.

        NPFMASS_REGIME = 2                          ! second index to NPFMASS(:,:)
        IF ( INCLUDE_KELVIN_EFFECT ) THEN
          RHF = 0.01D+00 * RHPERCENT
          A = ( 1.2D+00 - 0.0072D+00*(TK-273.15D+00) )
     &      * ( 1.0D+00 - 0.17D+00*(1.0D+00-RHF))
     &      * ( 1.0D+00 + 0.95D+00*(1.0D+00-RHF))
          RHF = RHF * EXP ( -A * RRNU_NM )
        ELSE
          RHF = 0.01D+00 * RHPERCENT
        ENDIF
        MFS = A0_AMSO4 + A1_AMSO4*RHF + A2_AMSO4*RHF*RHF + A3_AMSO4*RHF**3

      ENDIF

      IF ( MFS .GT. 0.0D+00 ) THEN   ! The mole fraction S in the particle should exceed zero.
        X_H2O = (1.0D+00 - MFS ) / MFS
      ELSE                           ! This case should not occur.
        X_H2O = 10.0D+00
      ENDIF

      IF( WRITE_LOG ) WRITE(AUNIT1,'(/A,F6.2,I6/)') 'X_NH3, NPFMASS_REGIME = ', X_NH3, NPFMASS_REGIME
      RETURN
      END SUBROUTINE EFFECTIVE_MW


      SUBROUTINE STEADY_STATE_H2SO4(PRS,RH,TEMP,XH2SO4_SS,SO4RATE,XNH3,KC,DT,XH2SO4_SS_WNPF)
!------------------------------------------------------------------------------------------------------------------
!     101706, DLW: Routine to estimate the steady-state concentration of sulfuric acid
!                  including the consumption of H2SO4 by new particle formation during the current time step.
!------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Input arguments.

      REAL(8), INTENT(IN)   :: PRS              ! pressure [Pa]
      REAL(8), INTENT(IN)   :: RH               ! fractional relative humidity [1]
      REAL(8), INTENT(IN)   :: TEMP             ! ambient temperature [K]
      REAL(8), INTENT(IN)   :: XH2SO4_SS        ! initial steady-state [H2SO4] (as SO4) 
                                                !   neglecting new particle formation [ugSO4/m^3]
      REAL(8), INTENT(IN)   :: SO4RATE          ! gas-phase H2SO4 (as SO4) production rate [ugSO4/m^3 s]
      REAL(8), INTENT(IN)   :: XNH3             ! ammonia mixing ratio [ppmV]
      REAL(8), INTENT(IN)   :: KC               ! condensational sink due to pre-existing aerosol [1/s]
      REAL(8), INTENT(IN)   :: DT               ! model physics time step [s]

      ! Output arguments.

      REAL(8), INTENT(OUT)  :: XH2SO4_SS_WNPF   ! steady-state [H2SO4] including new particle formation [ugSO4/m^3]

      ! Scratch local variables.

      INTEGER :: I                              ! loop counter
      REAL(8) :: DNDT                           ! new particle formation rate [particles/m^3/s]
      REAL(8) :: DMDT_SO4                       ! npf mass production rates [ugSO4/m^3/s]
      REAL(8) :: FX                             ! steady-state equation is FX = 0. [ugSO4/m^3/s]
      INTEGER, PARAMETER :: ITMAX = 100         ! loop limit for development code
      INTEGER, PARAMETER :: ICALLNPFRATE = 1    ! =0 impose mass limitation, >0 do not impose mass limitation
      REAL(8), PARAMETER :: XH2SO4_THRES_NCM3 = 1.001D+04  ! If [H2SO4] is below this NPF can be neglected.[#/cm^3] 
      REAL(8), PARAMETER :: XH2SO4_THRES = XH2SO4_THRES_NCM3 * MW_SO4 * 1.0D+12 / AVO   ! converted to [ugSO4/m^3] 
      REAL(8), PARAMETER :: EPS_XH2SO4_NCM3   = 1.00D+00                                ! tiny [H2SO4] [#/cm^3] 
      REAL(8), PARAMETER :: EPS_XH2SO4 = EPS_XH2SO4_NCM3 * MW_SO4 * 1.0D+12 / AVO       ! convert to [ugSO4/m^3] 
      REAL(8), PARAMETER :: REDUCTION_FACTOR = 1.2D+00   ! factor by which [H2SO4]ss is reduced each iteration [1] 
      LOGICAL, PARAMETER :: EARLY_RETURN = .FALSE.       ! flag for no-operation early exit 

      IF( XH2SO4_SS.LT.XH2SO4_THRES .OR. INUC.NE.3 .OR. EARLY_RETURN ) THEN     ! INUC=3 is the Napari et al. 
        XH2SO4_SS_WNPF = XH2SO4_SS                                              ! [ugSO4/m^3]
        RETURN
      ENDIF
      XH2SO4_SS_WNPF = XH2SO4_SS + EPS_XH2SO4                                   ! [ugSO4/m^3]
      CALL NPFRATE(PRS,RH,TEMP,XH2SO4_SS_WNPF,SO4RATE,XNH3,KC,DNDT,DMDT_SO4,ICALLNPFRATE) 
      FX = SO4RATE - KC*XH2SO4_SS_WNPF - DMDT_SO4                               ! evaluate function [ugSO4/m^3/s]
      IF( FX .GT. 0.0D+00 ) THEN
        ! WRITE(34,*)'FX(XMAX) .GT. 0.0D+00 in STEADY_STATE_H2SO4: FX = ', FX
        RETURN
      ENDIF
      ! WRITE(34,'(A,I5,5D13.4)')'I,X_SS,X_SS_WNPF,P,DMDT_SO4,FX=',0,XH2SO4_SS,XH2SO4_SS_WNPF,SO4RATE,DMDT_SO4,FX
!------------------------------------------------------------------------------------------------------------------
!     Reduce the steady-state H2SO4 until FX changes sign from negative to positive. 
!     Then the current value of XH2SO4_SS_WNPF is within a factor of REDUCTION_FACTOR 
!     of the actual steady-state value. 
!------------------------------------------------------------------------------------------------------------------
      DO I=1, ITMAX
        XH2SO4_SS_WNPF = XH2SO4_SS_WNPF / REDUCTION_FACTOR                      ! [ugSO4/m^3]
        CALL NPFRATE(PRS,RH,TEMP,XH2SO4_SS_WNPF,SO4RATE,XNH3,KC,DNDT,DMDT_SO4,ICALLNPFRATE) 
        FX = SO4RATE - KC*XH2SO4_SS_WNPF - DMDT_SO4                             ! evaluate function [ugSO4/m^3/s]
        ! WRITE(34,'(A,I5,5D13.4)')'I,X_SS,X_SS_WNPF,P,DMDT_SO4,FX=',I,XH2SO4_SS,XH2SO4_SS_WNPF,SO4RATE,DMDT_SO4,FX
        IF( FX .GE. 0.0D+00 ) EXIT
        IF( XH2SO4_SS_WNPF .LT. EPS_XH2SO4 ) RETURN                             ! [ugSO4/m^3]
      ENDDO
      RETURN
      END SUBROUTINE STEADY_STATE_H2SO4


      END MODULE AERO_NPF


