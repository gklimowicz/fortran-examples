      MODULE AERO_SETUP
!-------------------------------------------------------------------------------
!
!@sum     This module contains various aerosol microphysical variables and routines.
!@auth    Susanne Bauer/Doug Wright
!-------------------------------------------------------------------------------
      USE AERO_PARAM
      USE AERO_CONFIG
      IMPLICIT NONE

      INTEGER, SAVE :: MODE_NUMB_AKK, MODE_NUMB_ACC, MODE_NUMB_DD1, MODE_NUMB_DD2
      INTEGER, SAVE :: MODE_NUMB_DS1, MODE_NUMB_DS2, MODE_NUMB_SSA, MODE_NUMB_SSC
      INTEGER, SAVE :: MODE_NUMB_SSS, MODE_NUMB_OCC, MODE_NUMB_BC1, MODE_NUMB_BC2
      INTEGER, SAVE :: MODE_NUMB_BC3, MODE_NUMB_DBC, MODE_NUMB_BOC, MODE_NUMB_BCS
      INTEGER, SAVE :: MODE_NUMB_OCS, MODE_NUMB_MXX
      !-------------------------------------------------------------------------
      ! NMODES_XXXX is the number of modes containing species XXXX.
      ! NMODES_SEAS is the number of modes containing sea salt, not the 
      !             number of sea salt modes.
      ! NUMBER_OF_SEASALT_MODES is the number of sea salt modes, either 1 (SSS)
      !                         or 2 (SSA and SSC).
      !-------------------------------------------------------------------------
      INTEGER, SAVE :: NMODES_SULF, NMODES_BCAR, NMODES_OCAR, NMODES_DUST, NMODES_SEAS 
      INTEGER, SAVE :: NUMBER_OF_SEASALT_MODES
      INTEGER, SAVE :: NUMB_MAP(NWEIGHTS)                           ! [1]
      INTEGER, SAVE :: MASS_MAP(NWEIGHTS,NMASS_SPCS)                ! [1]
      INTEGER, SAVE :: PROD_INDEX(NWEIGHTS,NMASS_SPCS)              ! [1]
      INTEGER, SAVE, ALLOCATABLE :: SULF_MAP(:)                     ! [1]
      INTEGER, SAVE, ALLOCATABLE :: BCAR_MAP(:)                     ! [1]
      INTEGER, SAVE, ALLOCATABLE :: OCAR_MAP(:)                     ! [1]
      INTEGER, SAVE, ALLOCATABLE :: DUST_MAP(:)                     ! [1]
      INTEGER, SAVE, ALLOCATABLE :: SEAS_MAP(:)                     ! [1]
      INTEGER, SAVE, ALLOCATABLE :: SEAS_MODE_MAP(:)                ! [1]
      INTEGER, SAVE, ALLOCATABLE :: SEAS_MODE_MASS_MAP(:)           ! [1]
      INTEGER, SAVE, ALLOCATABLE :: MODE_NUMB_SEAS(:)               ! [1]
      INTEGER, SAVE :: NM(NWEIGHTS)                                 ! [1]
      CHARACTER(LEN=4), SAVE :: NM_SPC_NAME(NWEIGHTS,NMASS_SPCS)

      type GIKLQ_type
        integer :: n
        integer, allocatable :: l(:), k(:)
        integer, allocatable :: qq(:)
      end type GIKLQ_type

      type DIKL_type
        integer :: i
        integer :: k
        integer :: l
      end type DIKL_type
      
      type (GIKLQ_type), allocatable, save :: GIKLQ_control(:)
      type (DIKL_type), allocatable, save :: DIKL_control(:)
      integer, save :: nDIKL

      INTEGER, SAVE :: GIKLQ(NWEIGHTS,NWEIGHTS,NWEIGHTS,NMASS_SPCS) ! [1]
      INTEGER, SAVE :: DIKL (NWEIGHTS,NWEIGHTS,NWEIGHTS)            ! [1]
      INTEGER, SAVE :: DIJ  (NWEIGHTS,NWEIGHTS)                     ! [1]
      REAL(8), SAVE :: xDIJ (NWEIGHTS,NWEIGHTS) ! [1]
      REAL(8), SAVE :: RECIP_PART_MASS(NEMIS_SPCS)                  ! [1/ug]
      REAL(8), SAVE :: KCI_COEF_DP     (NWEIGHTS,NLAYS)             ! [m^3/s]
      REAL(8), SAVE :: KCI_COEF_DP_AEQ1(NWEIGHTS,NLAYS)             ! [m^3/s]
      REAL(8), SAVE :: THETA_POLY (NWEIGHTS)                        ! [1]
      REAL(8), SAVE :: DP0(NWEIGHTS)      ! default mode diameters of average mass [m]
                                          !   calculated from the DGN0 and SIG0 values
      REAL(8), SAVE :: DP0_EMIS(NWEIGHTS) ! mode diameters of average mass [m] for emissions,
                                          !   calculated from the DGN0_EMIS and SIG0_EMIS values
      !-------------------------------------------------------------------------
      ! DIFFCOEF_M2S(I) is the diffusivity of H2SO4 in air [m^2/s] for layer L.
      !-------------------------------------------------------------------------
      REAL(8), SAVE :: DIFFCOEF_M2S(NLAYS)
      !-------------------------------------------------------------------------
      ! KAPPAI(I) is the activating fraction for mode I.
      !-------------------------------------------------------------------------
      REAL(8), SAVE :: KAPPAI(NWEIGHTS)
      !-------------------------------------------------------------------------
      ! DENSPI(I) is the default particle density for mode I.
      ! DENS_COMP(I) is the density of chemical component I.
      ! RECIP_DENS_COMP(I) is the reciprocal density of chemical component I.
      !-------------------------------------------------------------------------
      REAL(8), SAVE :: DENSPI(NWEIGHTS)           ! [g/cm^3]
      REAL(8), SAVE :: DENS_COMP(NWEIGHTS)        ! [g/cm^3] 
      REAL(8), SAVE :: RECIP_DENS_COMP(NWEIGHTS)  ! [cm^3/g] 
      !-------------------------------------------------------------------------
      ! Characteristic lognormal parameters for each mode: DGN0 [um], SIG0 [1].            
      !-------------------------------------------------------------------------
      REAL(8), SAVE :: DGN0(NWEIGHTS)
      REAL(8), SAVE :: SIG0(NWEIGHTS)
      REAL(8), SAVE :: LNSIG0(NWEIGHTS)  ! ln( SIG0 )
      !-------------------------------------------------------------------------
      ! Lognormal parameters for emissions into each mode: 
      !   DGN0_EMIS [um], SIG0_EMIS [1].            
      !   These are used to convert mass emission rates to number emission rates.
      !-------------------------------------------------------------------------
      REAL(8), SAVE :: DGN0_EMIS(NWEIGHTS)
      REAL(8), SAVE :: SIG0_EMIS(NWEIGHTS)
      !-------------------------------------------------------------------------
      ! CONV_DPAM_TO_DGN(I) converts the diameter of average mass to the
      ! geometric mean diameter of the number size distribution for mode I,
      ! based on an assumed standard deviation for each mode.
      !-------------------------------------------------------------------------
      REAL(8), SAVE :: CONV_DPAM_TO_DGN(NWEIGHTS)
      !-------------------------------------------------------------------------
      ! EMIS_MODE_MAP and EMIS_SPCS_MAP have elements corresponding to 
      !   the aerosol types (in this order): AKK(=1), ACC(=2), BCC(=8), OCC(=7),
      !               DD1(=3), SSA(=5), SSC(=6), BOC(BC=8), BOC(OC=9), DD2(=10).
      ! EMIS_MODE_MAP(J) is mode number receiving the emissions held 
      !                  in EMIS_MASS(J).
      ! EMIS_SPCS_MAP(J) is the chemical species number (1-5) of the chemical
      !                  species held in EMIS_MASS(J).
      !-------------------------------------------------------------------------
      INTEGER,                  SAVE :: EMIS_MODE_MAP(NEMIS_SPCS)
      INTEGER, DIMENSION(NEMIS_SPCS) :: EMIS_SPCS_MAP = (/1,1,2,3,4,5,5,2,3,4/)
      !-------------------------------------------------------------------------
      ! The dimensions of these arrays depends upon mechanism.
      ! SEAS_MAP(I) is the mean mass per particle for sea salt mode I.
      !-------------------------------------------------------------------------
      CHARACTER(LEN= 3), SAVE, ALLOCATABLE :: MODE_NAME(:)      
      INTEGER,           SAVE, ALLOCATABLE :: MODE_SPCS(:,:)
      CHARACTER(LEN=16), SAVE, ALLOCATABLE :: AERO_SPCS(:)
      INTEGER,           SAVE, ALLOCATABLE :: ICOND(:)
      REAL,              SAVE, ALLOCATABLE :: RECIP_SEAS_MPP(:)         ! [1/ug]
      CHARACTER(LEN=3),  SAVE, ALLOCATABLE :: CITABLE(:,:)
      LOGICAL,                        SAVE :: INTERMODAL_TRANSFER
      !-------------------------------------------------------------------------------------------------------------------
      ! Diameter of average mass, averaged over all modes, used in the KK02 parameterization in subr. NPFRATE.
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8) :: AVG_DP_OF_AVG_MASS_METERS = 150.0D-09    ! [m] initial value; updated in subr. MATRIX.
      !-------------------------------------------------------------------------------------------------------------------
      ! Variables for the lookup table for condensational growth.
      !-------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: N_DP_CONDTABLE = 1000              ! [1] number of tabulated particle diameters
      REAL(8), PARAMETER :: DP_CONDTABLE_MIN = DPMIN_GLOBAL    ! [m] minimum ambient particle diameter
      REAL(8), PARAMETER :: DP_CONDTABLE_MAX = DPMAX_GLOBAL    ! [m] maximum ambient particle diameter
      REAL(8) :: DP_CONDTABLE(N_DP_CONDTABLE)                  ! [m] tabulated particle ambient diameters.
      REAL(8), SAVE :: XLN_SCALE_DP                            ! [1] ln of table diameter ratio.
      REAL(8), SAVE :: KCI_DP_CONDTABLE     (N_DP_CONDTABLE,NLAYS)  ! [m^3/s] tabulated condensational growth factors.
      REAL(8), SAVE :: KCI_DP_CONDTABLE_AEQ1(N_DP_CONDTABLE,NLAYS)  ! [m^3/s] tabulated condensational growth factors
                                                                    !   for the mass accommodation coefficient 
                                                                    !   set to unity.

      CONTAINS


      SUBROUTINE SETUP_CONFIG
!-------------------------------------------------------------------------------
!     Routine to initialize variables that depend upon choice of mechanism.
!-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,J,K,INDEX,IDIM
      LOGICAL :: FOUND
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      IF ( FIRSTIME ) THEN
        FIRSTIME = .FALSE.
        !-----------------------------------------------------------------------
        ! Get the number of modes in the selected mechanism.
        !-----------------------------------------------------------------------
        IDIM = 0        
        IF ( MECH .EQ. 1 ) IDIM = NM1
        IF ( MECH .EQ. 2 ) IDIM = NM2
        IF ( MECH .EQ. 3 ) IDIM = NM3
        IF ( MECH .EQ. 4 ) IDIM = NM4
        IF ( MECH .EQ. 5 ) IDIM = NM5
        IF ( MECH .EQ. 6 ) IDIM = NM6
        IF ( MECH .EQ. 7 ) IDIM = NM7
        IF ( MECH .EQ. 8 ) IDIM = NM8
        IF ( IDIM .NE. NMODES ) THEN
          WRITE(*,*)'ERROR in mechanism number MECH: MECH = ', MECH
          STOP
        ENDIF
        ALLOCATE( MODE_NAME(IDIM) )      
        ALLOCATE( MODE_SPCS(NMASS_SPCS,IDIM) )
        ALLOCATE( AERO_SPCS(NAEROBOX) )      
        ALLOCATE( ICOND(IDIM) )      
        ALLOCATE( CITABLE(IDIM,IDIM) )
        MODE_NAME(:)   = '   '
        MODE_SPCS(:,:) = 0
        AERO_SPCS(:)   = '                '     
        ICOND(:)       = 0
        CITABLE(:,:)   = '   '
      ENDIF
      !-------------------------------------------------------------------------
      ! Initialize arrays for coagulation interactions, mode names, mode
      ! species, condensation flags, and whether mode BC3 is present in the
      ! selected mechanism.
      !-------------------------------------------------------------------------
      IF     ( MECH .EQ. 1 ) THEN
        CITABLE(:,:) = CITABLE1(:,:)
        MODE_NAME(:) = MNAME(MODES1(:))
        MODE_SPCS(:,:) = MSPCS(:,MODES1(:))
        ICOND(:) = ICOND1(:)
        INCLUDE_BC3 = .TRUE.
      ELSEIF ( MECH .EQ. 2 ) THEN
        CITABLE(:,:) = CITABLE2(:,:)
        MODE_NAME(:) = MNAME(MODES2(:))
        MODE_SPCS(:,:) = MSPCS(:,MODES2(:))
        ICOND(:) = ICOND2(:)
        INCLUDE_BC3 = .FALSE.
      ELSEIF ( MECH .EQ. 3 ) THEN
        CITABLE(:,:) = CITABLE3(:,:)
        MODE_NAME(:) = MNAME(MODES3(:))
        MODE_SPCS(:,:) = MSPCS(:,MODES3(:))
        ICOND(:) = ICOND3(:)
        INCLUDE_BC3 = .FALSE.
      ELSEIF ( MECH .EQ. 4 ) THEN
        CITABLE(:,:) = CITABLE4(:,:)
        MODE_NAME(:) = MNAME(MODES4(:))
        MODE_SPCS(:,:) = MSPCS(:,MODES4(:))
        ICOND(:) = ICOND4(:)
        INCLUDE_BC3 = .FALSE.
      ELSEIF ( MECH .EQ. 5 ) THEN
        CITABLE(:,:) = CITABLE5(:,:)
        MODE_NAME(:) = MNAME(MODES5(:))
        MODE_SPCS(:,:) = MSPCS(:,MODES5(:))
        ICOND(:) = ICOND5(:)
        INCLUDE_BC3 = .TRUE.  
      ELSEIF ( MECH .EQ. 6 ) THEN
        CITABLE(:,:) = CITABLE6(:,:)
        MODE_NAME(:) = MNAME(MODES6(:))
        MODE_SPCS(:,:) = MSPCS(:,MODES6(:))
        ICOND(:) = ICOND6(:)
        INCLUDE_BC3 = .FALSE. 
      ELSEIF ( MECH .EQ. 7 ) THEN
        CITABLE(:,:) = CITABLE7(:,:)
        MODE_NAME(:) = MNAME(MODES7(:))
        MODE_SPCS(:,:) = MSPCS(:,MODES7(:))
        ICOND(:) = ICOND7(:)
        INCLUDE_BC3 = .FALSE. 
      ELSEIF ( MECH .EQ. 8 ) THEN
        CITABLE(:,:) = CITABLE8(:,:)
        MODE_NAME(:) = MNAME(MODES8(:))
        MODE_SPCS(:,:) = MSPCS(:,MODES8(:))
        ICOND(:) = ICOND8(:)
        INCLUDE_BC3 = .FALSE. 
      ENDIF
      !-------------------------------------------------------------------------
      ! Set INTERMODAL_TRANSFER according to setting in aero_param.f.
      ! If mode AKK is not present in the current mechanism, set to .FALSE.
      !-------------------------------------------------------------------------
      INTERMODAL_TRANSFER = SET_INTERMODAL_TRANSFER
      FOUND = .FALSE.
      DO I=1, NMODES
        IF ( MODE_NAME(I) .EQ. 'AKK' ) FOUND = .TRUE.
      ENDDO
      IF ( .NOT. FOUND ) THEN
        IF ( WRITE_LOG ) WRITE(AUNIT1,'(/A/)')
     &    'INTERMODAL TRANSFER (AKK->ACC) TURNED OFF SINCE MODE AKK IS ABSENT.'
        INTERMODAL_TRANSFER = .FALSE.
      ENDIF
      !-------------------------------------------------------------------------
      ! Check that all receptor modes in the CITABLE array are defined modes
      !   in the present mechanism.
      ! Also check that CITABLE is a symmetric matrix.
      !-------------------------------------------------------------------------
      DO I=1, NMODES  
        DO J=1, NMODES
          IF( CITABLE(I,J) .NE. CITABLE(J,I) ) THEN
            WRITE(*,*) 'CITABLE(I,J) must be a symmetric matrix.'
            WRITE(*,*) 'The CITABLE(I,J) set in aero_config.f is asymmetric for I, J = ', I, J
            STOP
          ENDIF
          FOUND = .FALSE.
          DO K=1, NMODES
            IF( CITABLE(I,J) .EQ. MODE_NAME(K) ) FOUND = .TRUE.
          ENDDO
          IF( CITABLE(I,J) .EQ. 'OFF' ) THEN 
            FOUND = .TRUE.    ! I-J interaction has been turned off
            ! WRITE(36,'(A,2I5,A5)')'I,J,CITABLE(I,J) = ', I,J,CITABLE(I,J)
          ENDIF
          IF( .NOT. FOUND ) THEN
            WRITE(*,*)'INVALID RECEPTOR MODE NAME: I, J, CITABLE(I,J) = ',
     &                                             I, J, CITABLE(I,J)
            STOP
          ENDIF
        ENDDO
      ENDDO
      !-------------------------------------------------------------------------
      ! Setup the indices to the AERO array.
      !-------------------------------------------------------------------------
      AERO_SPCS(1) = 'MASS_NITRATE'
      AERO_SPCS(2) = 'MASS_AMMONIUM'
      AERO_SPCS(3) = 'MASS_WATER'
      INDEX = 3          ! The first three values of INDEX (1, 2, 3) are already 
                         ! assigned to NO3, NH4, and H2O.
      IF ( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/2A/)') 'MODE #   MODE_NAME   CHEM_SPC #   CHEM_SPC_NAME   location in AERO',
     &                         '         AERO_SPCS'
        WRITE(AUNIT1,90) 0,'NO3',0,'ANO3',1,'MASS_NO3        '
        WRITE(AUNIT1,90) 0,'NH4',0,'ANH4',2,'MASS_NH4        '
        WRITE(AUNIT1,90) 0,'H2O',0,'AH2O',3,'MASS_H2O        '
      ENDIF
      DO I=1, NMODES
        DO J=1, NMASS_SPCS
          IF ( MODE_SPCS(J,I) .GT. 0 ) THEN  ! This mode contains species J.
            INDEX = INDEX + 1
            CALL SETUP_INDICES(I,J,INDEX,0)
            AERO_SPCS(INDEX) = 'MASS_'//MODE_NAME(I)//'_'//CHEM_SPC_NAME(J)
            IF ( WRITE_LOG ) THEN
              WRITE(AUNIT1,90) I,MODE_NAME(I),J,CHEM_SPC_NAME(J),INDEX,AERO_SPCS(INDEX)
              ! WRITE(31,'(8X,A23)') AERO_SPCS(INDEX)//' = OMIT'
            ENDIF
          ENDIF
        ENDDO
        DO J=1, NPOINTS
          INDEX = INDEX + 1
          CALL SETUP_INDICES(I,J,INDEX,1)     ! Set number conc. indices.
          IF ( J .EQ. 1 )  AERO_SPCS(INDEX) = 'NUMB_'//MODE_NAME(I)//'_1'   ! first  quadrature point
          IF ( J .EQ. 2 )  AERO_SPCS(INDEX) = 'NUMB_'//MODE_NAME(I)//'_2'   ! second quadrature point
          IF ( WRITE_LOG ) THEN
            WRITE(AUNIT1,90) I,MODE_NAME(I),J,'NUMB',INDEX,AERO_SPCS(INDEX)
            ! WRITE(31,'(8X,A23)') AERO_SPCS(INDEX)//' = OMIT'
          ENDIF
        ENDDO
      ENDDO
      IF ( INDEX .NE. NAEROBOX ) THEN
        WRITE(*,*)'INDEX .NE. NAEROBOX', INDEX, NAEROBOX    ! size of the AERO and AERO_SPCS arrays
        STOP
      ENDIF
      !-------------------------------------------------------------------------------------------------------------------
      ! Setup maps etc. needed to convert sea salt mass concentrations to number concentrations
      !   using mean particle masses for the sea salt modes.
      !-------------------------------------------------------------------------------------------------------------------
      CALL SETUP_SEASALT_MAPS
!-------------------------------------------------------------------------------------------------------------------------
!     Indices of the AERO array.
!-------------------------------------------------------------------------------------------------------------------------
!     WRITE(*,*)       MASS_NO3,   MASS_NH4,   MASS_H2O
!     WRITE(*,*)       NUMB_AKK_1, NUMB_AKK_2, MASS_AKK_SULF,  
!    &                 NUMB_ACC_1, NUMB_ACC_2, MASS_ACC_SULF,
!    &                 NUMB_DD1_1, NUMB_DD1_2, MASS_DD1_SULF,                               MASS_DD1_DUST, 
!    &                 NUMB_DS1_1, NUMB_DS1_2, MASS_DS1_SULF,                               MASS_DS1_DUST, 
!    &                 NUMB_DD2_1, NUMB_DD2_2, MASS_DD2_SULF,                               MASS_DD2_DUST, 
!    &                 NUMB_DS2_1, NUMB_DS2_2, MASS_DS2_SULF,                               MASS_DS2_DUST, 
!    &                 NUMB_SSA_1, NUMB_SSA_2, MASS_SSA_SULF,                                              MASS_SSA_SEAS, 
!    &                 NUMB_SSC_1, NUMB_SSC_2, MASS_SSC_SULF,                                              MASS_SSC_SEAS,
!    &                 NUMB_SSS_1, NUMB_SSS_2, MASS_SSS_SULF,                                              MASS_SSS_SEAS,
!    &                 NUMB_OCC_1, NUMB_OCC_2, MASS_OCC_SULF,                MASS_OCC_OCAR,
!    &                 NUMB_BC1_1, NUMB_BC1_2, MASS_BC1_SULF, MASS_BC1_BCAR,
!    &                 NUMB_BC2_1, NUMB_BC2_2, MASS_BC2_SULF, MASS_BC2_BCAR,
!    &                 NUMB_BC3_1, NUMB_BC3_2, MASS_BC3_SULF, MASS_BC3_BCAR,
!    &                 NUMB_OCS_1, NUMB_OCS_2, MASS_OCS_SULF,                MASS_OCS_OCAR,
!    &                 NUMB_DBC_1, NUMB_DBC_2, MASS_DBC_SULF, MASS_DBC_BCAR,                MASS_DBC_DUST,
!    &                 NUMB_BOC_1, NUMB_BOC_2, MASS_BOC_SULF, MASS_BOC_BCAR, MASS_BOC_OCAR,
!    &                 NUMB_BCS_1, NUMB_BCS_2, MASS_BCS_SULF, MASS_BCS_BCAR,
!    &                 NUMB_MXX_1, NUMB_MXX_2, MASS_MXX_SULF, MASS_MXX_BCAR, MASS_MXX_OCAR, MASS_MXX_DUST, MASS_MXX_SEAS
!-------------------------------------------------------------------------------------------------------------------------
   90 FORMAT(I6,9X,A3,9X,I4,12X,A4,I19,5X,A16)
      RETURN
      END SUBROUTINE SETUP_CONFIG


      SUBROUTINE SETUP_INDICES(I,J,INDEX,IN)
!-------------------------------------------------------------------------------
!     Routine to initialize indices of the AERO array.
!-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I, J, INDEX, IN
      INTEGER, PARAMETER :: OMIT = 0
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      IF ( FIRSTIME ) THEN
        FIRSTIME = .FALSE.
        !-----------------------------------------------------------------------
        ! Set all AERO indices to default values indicating that the species 
        ! is not defined in the current mechanism.
        !-----------------------------------------------------------------------
        MASS_AKK_SULF = OMIT
        NUMB_AKK_1    = OMIT
        MASS_ACC_SULF = OMIT
        NUMB_ACC_1    = OMIT
        MASS_DD1_SULF = OMIT
        MASS_DD1_DUST = OMIT
        NUMB_DD1_1    = OMIT
        MASS_DS1_SULF = OMIT
        MASS_DS1_DUST = OMIT
        NUMB_DS1_1    = OMIT
        MASS_DD2_SULF = OMIT
        MASS_DD2_DUST = OMIT
        NUMB_DD2_1    = OMIT
        MASS_DS2_SULF = OMIT
        MASS_DS2_DUST = OMIT
        NUMB_DS2_1    = OMIT
        MASS_SSA_SULF = OMIT
        MASS_SSA_SEAS = OMIT
        NUMB_SSA_1    = OMIT
        MASS_SSC_SULF = OMIT
        MASS_SSC_SEAS = OMIT
        NUMB_SSC_1    = OMIT
        MASS_SSS_SULF = OMIT
        MASS_SSS_SEAS = OMIT
        NUMB_SSS_1    = OMIT
        MASS_OCC_SULF = OMIT
        MASS_OCC_OCAR = OMIT
        NUMB_OCC_1    = OMIT
        MASS_BC1_SULF = OMIT
        MASS_BC1_BCAR = OMIT
        NUMB_BC1_1    = OMIT
        MASS_BC2_SULF = OMIT
        MASS_BC2_BCAR = OMIT
        NUMB_BC2_1    = OMIT
        MASS_BC3_SULF = OMIT
        MASS_BC3_BCAR = OMIT
        NUMB_BC3_1    = OMIT
        MASS_DBC_SULF = OMIT
        MASS_DBC_BCAR = OMIT
        MASS_DBC_DUST = OMIT
        NUMB_DBC_1    = OMIT
        MASS_BOC_SULF = OMIT
        MASS_BOC_BCAR = OMIT
        MASS_BOC_OCAR = OMIT
        NUMB_BOC_1    = OMIT
        MASS_BCS_SULF = OMIT
        MASS_BCS_BCAR = OMIT
        NUMB_BCS_1    = OMIT
        MASS_OCS_SULF = OMIT
        MASS_OCS_OCAR = OMIT
        NUMB_OCS_1    = OMIT
        MASS_MXX_SULF = OMIT
        MASS_MXX_BCAR = OMIT
        MASS_MXX_OCAR = OMIT
        MASS_MXX_DUST = OMIT
        MASS_MXX_SEAS = OMIT
        NUMB_MXX_1    = OMIT
      ENDIF
      IF ( IN .EQ. 0 ) THEN       ! This is a mass concentration.
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'AKK_SULF' ) MASS_AKK_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'ACC_SULF' ) MASS_ACC_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DD1_SULF' ) MASS_DD1_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DD1_DUST' ) MASS_DD1_DUST = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DS1_SULF' ) MASS_DS1_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DS1_DUST' ) MASS_DS1_DUST = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DD2_SULF' ) MASS_DD2_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DD2_DUST' ) MASS_DD2_DUST = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DS2_SULF' ) MASS_DS2_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DS2_DUST' ) MASS_DS2_DUST = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'SSA_SULF' ) MASS_SSA_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'SSA_SEAS' ) MASS_SSA_SEAS = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'SSC_SULF' ) MASS_SSC_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'SSC_SEAS' ) MASS_SSC_SEAS = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'SSS_SULF' ) MASS_SSS_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'SSS_SEAS' ) MASS_SSS_SEAS = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'OCC_SULF' ) MASS_OCC_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'OCC_OCAR' ) MASS_OCC_OCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BC1_SULF' ) MASS_BC1_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BC1_BCAR' ) MASS_BC1_BCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BC2_SULF' ) MASS_BC2_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BC2_BCAR' ) MASS_BC2_BCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BC3_SULF' ) MASS_BC3_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BC3_BCAR' ) MASS_BC3_BCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'OCS_SULF' ) MASS_OCS_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'OCS_OCAR' ) MASS_OCS_OCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DBC_SULF' ) MASS_DBC_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DBC_BCAR' ) MASS_DBC_BCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'DBC_DUST' ) MASS_DBC_DUST = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BOC_SULF' ) MASS_BOC_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BOC_BCAR' ) MASS_BOC_BCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BOC_OCAR' ) MASS_BOC_OCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BCS_SULF' ) MASS_BCS_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'BCS_BCAR' ) MASS_BCS_BCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'MXX_SULF' ) MASS_MXX_SULF = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'MXX_BCAR' ) MASS_MXX_BCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'MXX_OCAR' ) MASS_MXX_OCAR = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'MXX_DUST' ) MASS_MXX_DUST = INDEX
        IF ( MODE_NAME(I)//'_'//CHEM_SPC_NAME(J) .EQ. 'MXX_SEAS' ) MASS_MXX_SEAS = INDEX
      ELSEIF ( IN .EQ. 1 ) THEN   ! This is a number concentration.
        IF ( MODE_NAME(I).EQ.'AKK' .AND. J.EQ.1 ) NUMB_AKK_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'AKK' .AND. J.EQ.2 ) NUMB_AKK_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'ACC' .AND. J.EQ.1 ) NUMB_ACC_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'ACC' .AND. J.EQ.2 ) NUMB_ACC_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'DD1' .AND. J.EQ.1 ) NUMB_DD1_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'DD1' .AND. J.EQ.2 ) NUMB_DD1_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'DS1' .AND. J.EQ.1 ) NUMB_DS1_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'DS1' .AND. J.EQ.2 ) NUMB_DS1_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'DD2' .AND. J.EQ.1 ) NUMB_DD2_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'DD2' .AND. J.EQ.2 ) NUMB_DD2_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'DS2' .AND. J.EQ.1 ) NUMB_DS2_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'DS2' .AND. J.EQ.2 ) NUMB_DS2_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'SSA' .AND. J.EQ.1 ) NUMB_SSA_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'SSA' .AND. J.EQ.2 ) NUMB_SSA_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'SSC' .AND. J.EQ.1 ) NUMB_SSC_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'SSC' .AND. J.EQ.2 ) NUMB_SSC_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'SSS' .AND. J.EQ.1 ) NUMB_SSS_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'SSS' .AND. J.EQ.2 ) NUMB_SSS_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'OCC' .AND. J.EQ.1 ) NUMB_OCC_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'OCC' .AND. J.EQ.2 ) NUMB_OCC_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'BC1' .AND. J.EQ.1 ) NUMB_BC1_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'BC1' .AND. J.EQ.2 ) NUMB_BC1_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'BC2' .AND. J.EQ.1 ) NUMB_BC2_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'BC2' .AND. J.EQ.2 ) NUMB_BC2_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'BC3' .AND. J.EQ.1 ) NUMB_BC3_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'BC3' .AND. J.EQ.2 ) NUMB_BC3_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'OCS' .AND. J.EQ.1 ) NUMB_OCS_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'OCS' .AND. J.EQ.2 ) NUMB_OCS_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'DBC' .AND. J.EQ.1 ) NUMB_DBC_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'DBC' .AND. J.EQ.2 ) NUMB_DBC_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'BOC' .AND. J.EQ.1 ) NUMB_BOC_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'BOC' .AND. J.EQ.2 ) NUMB_BOC_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'BCS' .AND. J.EQ.1 ) NUMB_BCS_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'BCS' .AND. J.EQ.2 ) NUMB_BCS_2 = INDEX     
        IF ( MODE_NAME(I).EQ.'MXX' .AND. J.EQ.1 ) NUMB_MXX_1 = INDEX     
        IF ( MODE_NAME(I).EQ.'MXX' .AND. J.EQ.2 ) NUMB_MXX_2 = INDEX     
      ENDIF
      RETURN
      END SUBROUTINE SETUP_INDICES


      SUBROUTINE SETUP_SPECIES_MAPS
!-------------------------------------------------------------------------------
!     Routine to ... assign a mode number to each mode.
!                ... setup mass maps for species SULF, BCAR, OCAR, DUST, SEAS.
!                ... setup the number of mass species for mode I, NM(I).
!                ... setup the mass species names for mode I, NM_SPC_NAME(I).
!-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,J
      INTEGER :: ISULF,IBCAR,IOCAR,IDUST,ISEAS,INM
      INTEGER, PARAMETER :: INACTIVE = 0
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      IF ( FIRSTIME ) THEN
        FIRSTIME = .FALSE.  
        !-----------------------------------------------------------------------
        ! Get the number of modes containing each species: SULF, BCAR, OCAR,
        !   DUST, SEAS.
        !-----------------------------------------------------------------------
        ISULF = 0
        IBCAR = 0
        IOCAR = 0
        IDUST = 0
        ISEAS = 0
        DO I=1, NMODES
          DO J=1, NAEROBOX
            IF ( AERO_SPCS(J)(1:13).EQ.'MASS_'//MODE_NAME(I)//'_SULF' ) ISULF=ISULF+1
            IF ( AERO_SPCS(J)(1:13).EQ.'MASS_'//MODE_NAME(I)//'_BCAR' ) IBCAR=IBCAR+1
            IF ( AERO_SPCS(J)(1:13).EQ.'MASS_'//MODE_NAME(I)//'_OCAR' ) IOCAR=IOCAR+1
            IF ( AERO_SPCS(J)(1:13).EQ.'MASS_'//MODE_NAME(I)//'_DUST' ) IDUST=IDUST+1
            IF ( AERO_SPCS(J)(1:13).EQ.'MASS_'//MODE_NAME(I)//'_SEAS' ) ISEAS=ISEAS+1
          ENDDO
        ENDDO
        IF( WRITE_LOG ) THEN
          WRITE(AUNIT1,'(/A,5I4,/)')'Number of modes containing each species: ISULF,IBCAR,IOCAR,IDUST,ISEAS=',
     &                                                                        ISULF,IBCAR,IOCAR,IDUST,ISEAS
        ENDIF
        NMODES_SULF = ISULF
        NMODES_BCAR = IBCAR
        NMODES_OCAR = IOCAR
        NMODES_DUST = IDUST
        NMODES_SEAS = ISEAS
        ALLOCATE ( SULF_MAP(NMODES_SULF) ) 
        ALLOCATE ( BCAR_MAP(NMODES_BCAR) ) 
        ALLOCATE ( OCAR_MAP(NMODES_OCAR) ) 
        ALLOCATE ( DUST_MAP(NMODES_DUST) ) 
        ALLOCATE ( SEAS_MAP(NMODES_SEAS) ) 
        ALLOCATE ( MODE_NUMB_SEAS(NMODES_SEAS) ) 
        SULF_MAP(:) = 0
        BCAR_MAP(:) = 0
        OCAR_MAP(:) = 0
        DUST_MAP(:) = 0
        SEAS_MAP(:) = 0
        MODE_NUMB_SEAS(:) = 0
      ENDIF
      !-------------------------------------------------------------------------
      ! Assign a mode number to each mode.
      !-------------------------------------------------------------------------
      MODE_NUMB_AKK = INACTIVE
      MODE_NUMB_ACC = INACTIVE
      MODE_NUMB_DD1 = INACTIVE
      MODE_NUMB_DD2 = INACTIVE
      MODE_NUMB_DS1 = INACTIVE
      MODE_NUMB_DS2 = INACTIVE
      MODE_NUMB_SSA = INACTIVE
      MODE_NUMB_SSC = INACTIVE
      MODE_NUMB_SSS = INACTIVE
      MODE_NUMB_OCC = INACTIVE
      MODE_NUMB_BC1 = INACTIVE
      MODE_NUMB_BC2 = INACTIVE
      MODE_NUMB_BC3 = INACTIVE
      MODE_NUMB_DBC = INACTIVE
      MODE_NUMB_BOC = INACTIVE
      MODE_NUMB_BCS = INACTIVE
      MODE_NUMB_OCS = INACTIVE
      MODE_NUMB_MXX = INACTIVE
      DO I=1, NMODES
        IF ( MODE_NAME(I) .EQ. 'AKK' ) MODE_NUMB_AKK = I
        IF ( MODE_NAME(I) .EQ. 'ACC' ) MODE_NUMB_ACC = I
        IF ( MODE_NAME(I) .EQ. 'DD1' ) MODE_NUMB_DD1 = I
        IF ( MODE_NAME(I) .EQ. 'DD2' ) MODE_NUMB_DD2 = I
        IF ( MODE_NAME(I) .EQ. 'DS1' ) MODE_NUMB_DS1 = I
        IF ( MODE_NAME(I) .EQ. 'DS2' ) MODE_NUMB_DS2 = I
        IF ( MODE_NAME(I) .EQ. 'SSA' ) MODE_NUMB_SSA = I
        IF ( MODE_NAME(I) .EQ. 'SSC' ) MODE_NUMB_SSC = I
        IF ( MODE_NAME(I) .EQ. 'SSS' ) MODE_NUMB_SSS = I
        IF ( MODE_NAME(I) .EQ. 'OCC' ) MODE_NUMB_OCC = I
        IF ( MODE_NAME(I) .EQ. 'BC1' ) MODE_NUMB_BC1 = I
        IF ( MODE_NAME(I) .EQ. 'BC2' ) MODE_NUMB_BC2 = I
        IF ( MODE_NAME(I) .EQ. 'BC3' ) MODE_NUMB_BC3 = I
        IF ( MODE_NAME(I) .EQ. 'DBC' ) MODE_NUMB_DBC = I
        IF ( MODE_NAME(I) .EQ. 'BOC' ) MODE_NUMB_BOC = I
        IF ( MODE_NAME(I) .EQ. 'BCS' ) MODE_NUMB_BCS = I
        IF ( MODE_NAME(I) .EQ. 'OCS' ) MODE_NUMB_OCS = I
        IF ( MODE_NAME(I) .EQ. 'MXX' ) MODE_NUMB_MXX = I
      ENDDO  
      IF ( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A/)')        'MODE_NAME( MODE_NUMB_XXX ), MODE_NUMB_XXX'
        IF ( MODE_NUMB_AKK .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_AKK ), MODE_NUMB_AKK
        ENDIF
        IF ( MODE_NUMB_ACC .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_ACC ), MODE_NUMB_ACC
        ENDIF
        IF ( MODE_NUMB_DD1 .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_DD1 ), MODE_NUMB_DD1
        ENDIF
        IF ( MODE_NUMB_DS1 .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_DS1 ), MODE_NUMB_DS1
        ENDIF
        IF ( MODE_NUMB_DD2 .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_DD2 ), MODE_NUMB_DD2
        ENDIF
        IF ( MODE_NUMB_DS2 .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_DS2 ), MODE_NUMB_DS2
        ENDIF
        IF ( MODE_NUMB_SSA .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_SSA ), MODE_NUMB_SSA
        ENDIF
        IF ( MODE_NUMB_SSC .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_SSC ), MODE_NUMB_SSC
        ENDIF
        IF ( MODE_NUMB_SSS .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_SSS ), MODE_NUMB_SSS
        ENDIF
        IF ( MODE_NUMB_OCC .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_OCC ), MODE_NUMB_OCC
        ENDIF
        IF ( MODE_NUMB_BC1 .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_BC1 ), MODE_NUMB_BC1
        ENDIF
        IF ( MODE_NUMB_BC2 .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_BC2 ), MODE_NUMB_BC2
        ENDIF
        IF ( MODE_NUMB_BC3 .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_BC3 ), MODE_NUMB_BC3
        ENDIF
        IF ( MODE_NUMB_DBC .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_DBC ), MODE_NUMB_DBC
        ENDIF
        IF ( MODE_NUMB_BOC .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_BOC ), MODE_NUMB_BOC
        ENDIF
        IF ( MODE_NUMB_BCS .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_BCS ), MODE_NUMB_BCS
        ENDIF
        IF ( MODE_NUMB_OCS .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_OCS ), MODE_NUMB_OCS
        ENDIF
        IF ( MODE_NUMB_MXX .GT. 0 ) THEN
          WRITE(AUNIT1,'(3X,A3,3X,I4)') MODE_NAME( MODE_NUMB_MXX ), MODE_NUMB_MXX
        ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      ! Setup NUMB_MAP(I): location of the Ith number  conc. in the AERO array
      ! Setup SULF_MAP(I): location of the Ith sulfate conc. in the AERO array
      ! Setup BCAR_MAP(I): location of the Ith BC      conc. in the AERO array
      ! Setup OCAR_MAP(I): location of the Ith OC      conc. in the AERO array
      ! Setup DUST_MAP(I): location of the Ith dust    conc. in the AERO array
      ! Setup SEAS_MAP(I): location of the Ith sea salt conc. in the AERO array
      ! Setup NM(I):       number of mass concs. defined for mode I      
      ! Setup NM_SPC_NAME(I,J): name of the Jth mass conc. in mode I
      !-------------------------------------------------------------------------
      IBCAR = 1
      IOCAR = 1
      IDUST = 1
      ISEAS = 1
      NM(:) = 0
      NM_SPC_NAME(:,:) = '    '     ! LEN=4 character variable.
      DO I=1, NMODES
        INM = 0
        DO J=1, NAEROBOX
          IF ( AERO_SPCS(J)(1:8)  .EQ. 'NUMB_'//MODE_NAME(I) )          NUMB_MAP(I) = J
          IF ( AERO_SPCS(J)(1:13) .EQ. 'MASS_'//MODE_NAME(I)//'_SULF' ) SULF_MAP(I) = J
          IF ( AERO_SPCS(J)(1:13) .EQ. 'MASS_'//MODE_NAME(I)//'_BCAR' ) THEN
            BCAR_MAP(IBCAR) = J          ! location of this mass conc. in the AERO array
            IBCAR = IBCAR + 1
          ENDIF
          IF ( AERO_SPCS(J)(1:13) .EQ. 'MASS_'//MODE_NAME(I)//'_OCAR' ) THEN
            OCAR_MAP(IOCAR) = J          ! location of this mass conc. in the AERO array
            IOCAR = IOCAR + 1
          ENDIF
          IF ( AERO_SPCS(J)(1:13) .EQ. 'MASS_'//MODE_NAME(I)//'_DUST' ) THEN
            DUST_MAP(IDUST) = J          ! location of this mass conc. in the AERO array
            IDUST = IDUST + 1
          ENDIF
          IF ( AERO_SPCS(J)(1:13) .EQ. 'MASS_'//MODE_NAME(I)//'_SEAS' ) THEN
            SEAS_MAP(ISEAS) = J          ! location of this mass conc. in the AERO array
            MODE_NUMB_SEAS(ISEAS) = I    ! mode number for this sea salt-containing mode
            ISEAS = ISEAS + 1
          ENDIF
          IF ( AERO_SPCS(J)(1:8) .EQ. 'MASS_'//MODE_NAME(I) ) THEN
            INM = INM + 1
            NM_SPC_NAME(I,INM) = AERO_SPCS(J)(10:13)
          ENDIF
        ENDDO
        NM(I) = INM
      ENDDO
      IF ( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A11,16I4)') 'NUMB_MAP = ', NUMB_MAP(:)
        WRITE(AUNIT1,'(/A11,16I4)') 'SULF_MAP = ', SULF_MAP(:)
        WRITE(AUNIT1,'(/A11,16I4)') 'BCAR_MAP = ', BCAR_MAP(:)
        WRITE(AUNIT1,'(/A11,16I4)') 'OCAR_MAP = ', OCAR_MAP(:)
        WRITE(AUNIT1,'(/A11,16I4)') 'DUST_MAP = ', DUST_MAP(:)
        WRITE(AUNIT1,'(/A11,16I4)') 'SEAS_MAP = ', SEAS_MAP(:)
        WRITE(AUNIT1,'(/A11,16I4)') 'NM(I)    = ', NM(:)
        DO I=1, NMODES
          WRITE(AUNIT1,'(/A40,A5,I4,5A5)') 'MODE_NAME(I), NM(I), NM_SPC_NAME(I,:) = ', MODE_NAME(I), NM(I), NM_SPC_NAME(I,:)
        ENDDO
      ENDIF
      !----------------------------------------------------------------------------------------------------
      ! Setup EMIS_MODE_MAP. There are presently 10 emitted species.
      !
      ! EMIS_MODE_MAP and EMIS_SPCS_MAP have elements corresponding to the aerosol types (in this order):
      !   AKK(=1), ACC(=2), BCC(=8), OCC(=7), DD1(=3), SSA(=5), SSC(=6), BOC(BC=8), BOC(OC=9), DD2(=10).
      !
      ! EMIS_MODE_MAP(J) is mode number receiving the emissions held in EMIS_MASS(J).
      ! EMIS_SPCS_MAP(J) is the chemical species number (1-5) of the chemical species held in EMIS_MASS(J).
      ! EMIS_SPCS_MAP = (/1,1,2,3,4,5,5,2,3,4/) is set at the top of the module.
      !----------------------------------------------------------------------------------------------------
      EMIS_MODE_MAP(:) = 0 
      EMIS_MODE_MAP(1) = 1  ! Aitken mode sulfate always goes in the first mode, whether it is AKK or ACC.
      DO I=1, NMODES        ! If no Aitken mode, then the accumulation mode is the first mode.
        IF( MODE_NAME(I) .EQ. 'ACC' ) EMIS_MODE_MAP(2) = I
        IF( MODE_NAME(I) .EQ. 'BC1' ) EMIS_MODE_MAP(3) = I
        IF( MODE_NAME(I) .EQ. 'OCC' ) EMIS_MODE_MAP(4) = I
        IF( MODE_NAME(I) .EQ. 'DD1' ) THEN
          EMIS_MODE_MAP(5) = I
          IF( MECH .GE. 5 .AND. MECH .LE. 8 ) EMIS_MODE_MAP(10) = I   ! emissions for both dust modes go into mode DD1
        ENDIF
        IF( MODE_NAME(I) .EQ. 'SSA' ) EMIS_MODE_MAP(6) = I
        IF( MODE_NAME(I) .EQ. 'SSC' ) EMIS_MODE_MAP(7) = I
        IF( MODE_NAME(I) .EQ. 'SSS' ) EMIS_MODE_MAP(6) = I  ! emissions for both sea salt modes go into mode SSS
        IF( MODE_NAME(I) .EQ. 'SSS' ) EMIS_MODE_MAP(7) = I  ! emissions for both sea salt modes go into mode SSS
        IF( MODE_NAME(I) .EQ. 'BOC' ) EMIS_MODE_MAP(8) = I
        IF( MODE_NAME(I) .EQ. 'BOC' ) EMIS_MODE_MAP(9) = I
        IF( MODE_NAME(I) .EQ. 'DD2' ) EMIS_MODE_MAP(10) = I
      ENDDO
      !-------------------------------------------------------------------------
      ! If the mechanism does not have the mode BOC to receive the 
      ! mixed BC-OC emissions, put these directly into the BC1 and OCC modes.
      !-------------------------------------------------------------------------
      IF ( EMIS_MODE_MAP(8) .EQ. 0 ) THEN   ! the BC in BC-OC emissions
        DO I=1, NMODES   
          IF( MODE_NAME(I) .EQ. 'BC1' ) THEN
            EMIS_MODE_MAP(8) = I
            IF ( WRITE_LOG ) WRITE(AUNIT1,'(/2A/)')'BC of BO-OC put into mode ',MODE_NAME(I)
          ENDIF
        ENDDO
      ENDIF
      IF ( EMIS_MODE_MAP(9) .EQ. 0 ) THEN   ! the OC in BC-OC emissions
        DO I=1, NMODES   
          IF( MODE_NAME(I) .EQ. 'OCC' ) THEN
            EMIS_MODE_MAP(9) = I
            IF ( WRITE_LOG ) WRITE(AUNIT1,'(/2A/)')'OC of BO-OC put into mode ',MODE_NAME(I)
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SETUP_SPECIES_MAPS


      SUBROUTINE SETUP_SEASALT_MAPS
!---------------------------------------------------------------------------------------------------------
!     Routine to setup maps and variables needed to derive sea salt 
!     number concentrations from sea salt mass concentrations and
!     characteristic mean masses per particle.
!
!     Three arrays are set up:
!
!       SEAS_MODE_MAP(II)      ! mode number of the IIth SS mode
!       SEAS_MODE_MASS_MAP(II) ! location in the AERO array of the SS mass for the IIth SS mode
!       RECIP_SEAS_MPP(II)     ! reciprocal of the mean mass per particle for the IIth SS mode
!
!     There are either two sea salt modes, SSA and SSC, or only one, SSS.
!---------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,IDIM,II,J
      REAL(8) :: DPS
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      IF ( FIRSTIME ) THEN
        FIRSTIME = .FALSE.
        IDIM = 0
        DO I=1, NMODES
          IF ( MODE_NAME(I)(1:2) .EQ. 'SS' ) IDIM = IDIM + 1
        ENDDO
        IF ( IDIM .LT. 1 .OR. IDIM .GT. 2 ) THEN
          WRITE(*,*)'SUBROUTIINE SETUP_SEASALT_MAPS:'
          WRITE(*,*)'NUMBER OF SEAS SALT MODES IS ZERO OR GREATER THAN TWO - INCORRECT'
          STOP
        ENDIF
        NUMBER_OF_SEASALT_MODES = IDIM
        ALLOCATE( SEAS_MODE_MAP( IDIM ) )     
        ALLOCATE( SEAS_MODE_MASS_MAP( IDIM ) )
        ALLOCATE( RECIP_SEAS_MPP( IDIM ) )    
        SEAS_MODE_MAP(:) = 0
        SEAS_MODE_MASS_MAP(:) = 0
        RECIP_SEAS_MPP(:) = 0.0D+00
        II = 0
        DO I=1, NMODES
          IF ( MODE_NAME(I)(1:2) .EQ. 'SS' ) THEN  ! This is a sea salt mode, either SSA, SSC, or SSS.
            II = II + 1
            SEAS_MODE_MAP(II) = I
            DO J=1, NAEROBOX
              IF( AERO_SPCS(J)(1:13) .EQ. 'MASS_'//MODE_NAME(I)//'_SEAS' ) THEN ! This is the sea salt mass in the mode.
                SEAS_MODE_MASS_MAP(II) = J
                ! WRITE(*,*) AERO_SPCS(J)(1:13)     ! Checked: the correct species are identified. 
              ENDIF
            ENDDO      
            !-------------------------------------------------------------------------------------------------
            ! Calculate the diameter of average mass for each sea salt mode for the emissions lognormals.
            !
            ! The emissions lognormals are used since it is the dry sea salt concentration that
            ! is divided by the average dry mass per particle (in subr. matrix) to obtain the current
            ! number concentration from the dry sea salt mass concentration for each sea salt mode.
            ! The dry NaCl mass per particle changes little over time for sea salt, given low coagulation rates.
            ! Accreted water, sulfate, nitrate, and ammonium are irrelevant here. The dry NaCl per emitted
            ! sea salt particle is the appropriate dry mass to divide into the current dry NaCl concentration
            ! to obtain the particle number concentration.
            ! 
            ! Dam [um] = ( diameter moment 3  /  diameter moment 0    )**(1/3)
            !          = ( dg**3 * sg**9                              )**(1/3)
            !          = ( dg**3 * [ exp( 0.5*(log(sigmag))**2 ) ]**9 )**(1/3)
            !          =   dg    * [ exp( 0.5*(log(sigmag))**2 ) ]**3
            !          =   dg    * [ exp( 1.5*(log(sigmag))**2 ) ]
            !-------------------------------------------------------------------------------------------------
            IF( MODE_NAME(I).EQ.'SSA') DPS = 1.0D-06 * DG_SSA_EMIS * EXP( 1.5D+00 * ( LOG(SG_SSA_EMIS) )**2 )
            IF( MODE_NAME(I).EQ.'SSC') DPS = 1.0D-06 * DG_SSC_EMIS * EXP( 1.5D+00 * ( LOG(SG_SSC_EMIS) )**2 )
            IF( MODE_NAME(I).EQ.'SSS') DPS = 1.0D-06 * DG_SSS_EMIS * EXP( 1.5D+00 * ( LOG(SG_SSS_EMIS) )**2 )
            IF( ACTIVATION_COMPARISON .AND. ( MECH.EQ.4 .OR. MECH.EQ.8 ) ) THEN   ! Activation test. 
              DPS = 1.0D-06 * 0.3308445477D+00 * EXP( 1.5D+00 * ( LOG( SG_SSS_EMIS) )**2 )
              WRITE(*,*)'Special value set for DPS and RECIP_SEAS_MPP for mode SSS in aero_setup.f.'
            ENDIF
            !-------------------------------------------------------------------------------------------------
            ! DPS is the diameter of average mass for mode J for the dry emitted sea salt, and when cubed
            ! and multiplied by pi/6 it yields the average dry sea salt particle volume in emissions mode J.
            ! Multiplication by the emitted dry particle sea salt density then yields the average
            ! dry sea salt mass per particle emitted into the mode. 
            !-------------------------------------------------------------------------------------------------
            IF( DISCRETE_EVAL_OPTION ) THEN
              RECIP_SEAS_MPP(II) = 1.0D+00 / ( 1.0D+12 * DENSP          * PI6 * DPS**3 )
            ELSE
              RECIP_SEAS_MPP(II) = 1.0D+00 / ( 1.0D+12 * EMIS_DENS_SEAS * PI6 * DPS**3 )
            ENDIF
            IF( WRITE_LOG ) THEN
              WRITE(AUNIT1,'(/A/)')'II,I,MODE_NAME(I),AERO_SPCS(SEAS_MODE_MASS_MAP(II)),RECIP_SEAS_MPP(II)'
              WRITE(AUNIT1,90000)   II,I,MODE_NAME(I),AERO_SPCS(SEAS_MODE_MASS_MAP(II)),RECIP_SEAS_MPP(II)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
90000 FORMAT(2I4,4X,A6,4X,A16,4X,D15.5)
      RETURN
      END SUBROUTINE SETUP_SEASALT_MAPS


      SUBROUTINE SETUP_KCI         
!-----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the coefficients that multiply the number
!     concentrations, or the number concentrations times the particle diameters,
!     to obtain the condensational sink for each mode or quadrature point.
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I, L     ! indices
      REAL    :: SIGMA    ! See subr. ATMOSPHERE below.
      REAL    :: DELTA    ! See subr. ATMOSPHERE below.
      REAL    :: THETA    ! See subr. ATMOSPHERE below.
      REAL(8) :: P        ! ambient pressure [Pa]
      REAL(8) :: T        ! ambient temperature [K]
      REAL(8) :: D        ! molecular diffusivity of H2SO4 in air [m^2/s]
      REAL(8) :: C        ! mean molecular speed of H2SO4 [m/s]
      REAL(8) :: LA       ! mean free path in air [m]
                          ! 6.6328D-08 is the sea level value given in Table I.2.8
                          ! on p.10 of U.S. Standard Atmosphere 1962
      REAL(8) :: LH       ! mean free path of H2SO4 in air [m]
      REAL(8) :: BETA     ! transition regime correction to the condensational flux [1]
      REAL(8) :: KN       ! Knudsen number for H2SO4 in air [1]
      REAL(8) :: THETAI   ! monodispersity correction factor [1] (Okuyama et al. 1988)
      REAL(8) :: SCALE_DP ! scale factor for defined table of ambient particle diameters [1]

      REAL(8), PARAMETER :: ALPHA    = 0.86D+00      ! mass accommodation coefficient [1] (Hansen, 2005)
      REAL(8), PARAMETER :: D_STDATM = 1.2250D+00    ! sea-level std. density [kg/m^3]
      REAL(8), PARAMETER :: P_STDATM = 101325.0D+00  ! sea-level std. pressure [Pa]
      REAL(8), PARAMETER :: T_STDATM = 288.15D+00    ! sea-level std. temperature [K]
      REAL(8), PARAMETER :: P0       = 101325.0D+00  ! reference pressure    [Pa] for D
      REAL(8), PARAMETER :: T0       = 273.16D+00    ! reference temperature [K]  for D
      REAL(8), PARAMETER :: D0       = 9.36D-06      ! diffusivity of H2SO4 in air [m^2/s]
                                                     ! calculated using Eqn 11-4.4 of Reid 
                                                     ! et al.(1987) at 273.16 K and 101325 Pa

      IF( WRITE_LOG ) WRITE(AUNIT1,90002)'I','L','ZHEIGHT','P','T','D','C','LA','LH',
     &                                   'DP0','KN','THETAI','BETA'
      KCI_COEF_DP          (:,:) = 0.0D+00   ! mass accommodation coefficient is arbitrary
      KCI_COEF_DP_AEQ1     (:,:) = 0.0D+00   ! mass accommodation coefficient is unity
      KCI_DP_CONDTABLE     (:,:) = 0.0D+00   ! mass accommodation coefficient is arbitrary
      KCI_DP_CONDTABLE_AEQ1(:,:) = 0.0D+00   ! mass accommodation coefficient is unity
      !----------------------------------------------------------------------------------------------------------------
      ! Setup table of ambient particle diameters.
      !----------------------------------------------------------------------------------------------------------------
      SCALE_DP = ( DP_CONDTABLE_MAX / DP_CONDTABLE_MIN )**(1.0D+00/REAL(N_DP_CONDTABLE-1))
      XLN_SCALE_DP = LOG( SCALE_DP )
      DO I=1, N_DP_CONDTABLE
        DP_CONDTABLE(I) = DP_CONDTABLE_MIN * SCALE_DP**(I-1)              ! [m]
      ENDDO

      DO L=1, NLAYS
        CALL ATMOSPHERE( REAL( ZHEIGHT(L) ), SIGMA, DELTA, THETA )
        P = DELTA * P_STDATM                                           ! [Pa}
        T = THETA * T_STDATM                                           ! [K]
        D = D0 * ( P0 / P ) * ( T / T0 )**1.75                         ! [m^2/s]
        DIFFCOEF_M2S(L) = D                                            ! [m^2/s]
        C = SQRT( 8.0D+00 * RGAS_SI * T / ( PI * MW_H2SO4*1.0D-03 ) )  ! [m/s]
        LA = 6.6332D-08 * ( P_STDATM / P ) * ( T / T_STDATM )          ! [m]
        LH = 3.0D+00 * D / C                                           ! [m]
        DO I=1, NWEIGHTS
          THETAI = EXP( - ( LOG(SIG0(I)) )**2 )                        ! [1]
          THETA_POLY(I) = THETAI                                       ! [1] polydispersity adjustment factor
          !------------------------------------------------------------------------------------------------------------
          ! For the condensation sink for general use, the mean free path is
          ! that for the condensing vapor (H2SO4), and the mass accommodation coefficient is adjustable. 
          !------------------------------------------------------------------------------------------------------------
          KN = 2.0D+00 * LH / DP0(I)                                   ! LH and DP0 in [m]
          BETA = ( 1.0D+00 + KN )                                      ! [1]
     &         / ( 1.0D+00 + 0.377D+00*KN + 1.33D+00*KN*(1.0D+00 + KN)/ALPHA ) 
          KCI_COEF_DP(I,L)      = 2.0D+00 * PI * THETAI * D * BETA * DP0(I) ! [m^3/s] may be updated in subr. matrix
          !------------------------------------------------------------------------------------------------------------
          ! For the condensation sink for use in Kerminen and Kulmala (2002), the mean free path is
          ! that for air, and the mass accommodation coefficient is set to unity. 
          !------------------------------------------------------------------------------------------------------------
          KN = 2.0D+00 * LA / DP0(I)                                   ! LA and DP0 in [m]
          BETA = ( 1.0D+00 + KN )                                      ! [1]
     &         / ( 1.0D+00 + 0.377D+00*KN + 1.33D+00*KN*(1.0D+00 + KN)/1.0D+00 ) 
          KCI_COEF_DP_AEQ1(I,L) = 2.0D+00 * PI * THETAI * D * BETA * DP0(I) ! [m^3/s] may be updated in subr. matrix
          ! IF( WRITE_LOG ) WRITE(AUNIT1,90000)I,L,ZHEIGHT(L),P,T,D,C,LA,LH,DP0(I)*1.0D+06,KN,THETAI,BETA
        ENDDO
        DO I=1, N_DP_CONDTABLE
          !------------------------------------------------------------------------------------------------------------
          ! For the condensation sink for general use, the mean free path is
          ! that for the condensing vapor (H2SO4), and the mass accommodation coefficient is adjustable. 
          ! The THETAI factor is included later in aero_matrix.f.
          !------------------------------------------------------------------------------------------------------------
          KN = 2.0D+00 * LH / DP_CONDTABLE(I)                          ! LH and DP in [m]
          BETA = ( 1.0D+00 + KN )                                      ! [1]
     &         / ( 1.0D+00 + 0.377D+00*KN + 1.33D+00*KN*(1.0D+00 + KN)/ALPHA ) 
          KCI_DP_CONDTABLE(I,L)      = 2.0D+00 * PI * D * BETA * DP_CONDTABLE(I)  ! [m^3/s]
          !------------------------------------------------------------------------------------------------------------
          ! For the condensation sink for use in Kerminen and Kulmala (2002), the mean free path is
          ! that for air, and the mass accommodation coefficient is set to unity. 
          ! The THETAI factor is included later in aero_matrix.f.
          !------------------------------------------------------------------------------------------------------------
          KN = 2.0D+00 * LA / DP_CONDTABLE(I)                          ! LA and DP in [m]
          BETA = ( 1.0D+00 + KN )                                      ! [1]
     &         / ( 1.0D+00 + 0.377D+00*KN + 1.33D+00*KN*(1.0D+00 + KN)/1.0D+00 ) 
          KCI_DP_CONDTABLE_AEQ1(I,L) = 2.0D+00 * PI * D * BETA * DP_CONDTABLE(I)  ! [m^3/s]
          ! WRITE(AUNIT1,90000) I, L, ZHEIGHT(L), P, T, D, C, LA, LH, DP_CONDTABLE(I)*1.0D+06, KN, THETAI, BETA
        ENDDO
      ENDDO
      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A/)')'  I  L     KCI_COEF_DP[m^3/s] KCI_COEF_DP_AEQ1[m^3/s]'    
        DO L=1, NLAYS
          DO I=1, NWEIGHTS
            WRITE(AUNIT1,90001) I, L, KCI_COEF_DP(I,L), KCI_COEF_DP_AEQ1(I,L)
          ENDDO
        ENDDO
      ENDIF
      DO I=1, NMODES
        IF ( ICOND(I) .EQ. 0 ) THEN 
          KCI_COEF_DP     (I,:) = 0.0D+00
          KCI_COEF_DP_AEQ1(I,:) = 0.0D+00
          THETA_POLY(I)         = 0.0D+00
          IF( WRITE_LOG ) WRITE(AUNIT1,'(/2A/)') 'Condensational growth turned off for mode ', MODE_NAME(I)
        ENDIF
      ENDDO

90000 FORMAT( 2I3,F8.4,F9.1,F7.2,D10.3,F6.1,2D10.3,4F8.4)
90002 FORMAT(/2A3,A8,  A9,  A7,  A10,  A6,  2A10,  4A8 /)
90001 FORMAT(2I3,4D20.4)
      RETURN
      END SUBROUTINE SETUP_KCI 


      SUBROUTINE SETUP_EMIS        
!-------------------------------------------------------------------------------
!     Routine to calculate the reciprocal of the average particle mass [ug]
!     for each mode for conversion of mode mass emission rates [ug/m^3/s] 
!     to mode number emission rates [#/m^3/s].
!
!     The factor 1.0D+12 converts [m^3] to [cm^3] and [g] to [ug].
!     DP0_EMIS(:) is in [m].
!
!     EMIS_MODE_MAP has elements corresponding to the aerosol types (in this order):
!       AKK(=1), ACC(=2), BCC(=8), OCC(=7),
!       DD1(=3), SSA(=5), SSC(=6), BOC(BC=8), BOC(OC=9), DD2(=10).
!     EMIS_MODE_MAP(J) is mode number receiving the emissions held in EMIS_MASS(J).
!-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I, J

      RECIP_PART_MASS(:) = 0.0D+00    ! [1/ug]

      IF ( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/2A5,A17,A32,A32/)') 'I','J','  DP0_EMIS(J)[um]',
     &               'RECIP_PART_MASS(I)[1/ug]','EMIS_DENS(I)[g/cm^3]'
      ENDIF
      DO I=1, NEMIS_SPCS       ! currently, NEMIS_SPCS = 10 emitted species
        J = EMIS_MODE_MAP(I)   ! J ranges over the number of modes.
        !-----------------------------------------------------------------------
        ! DP0_EMIS(J) is the diameter of average mass for mode J, and when cubed
        ! and multiplied by pi/6 it yields the average particle volume in mode J.
        ! Multiplication by the emitted particle density then yields the average
        ! mass per particle emitted into the mode.
        ! 
        ! EMIS_DENS(:)       ranges over the emission species only. 
        ! RECIP_PART_MASS(:) ranges over the emission species only. 
        ! DP0_EMIS(:)        ranges over all modes in the mechanism.
        !-----------------------------------------------------------------------
        RECIP_PART_MASS(I) = 1.0D+00 / ( 1.0D+12 * EMIS_DENS(I) * PI6 * DP0_EMIS(J)**3 )
        IF( WRITE_LOG) WRITE(AUNIT1,90000) I, J, DP0_EMIS(J)*1.0D+06, RECIP_PART_MASS(I), EMIS_DENS(I)
      ENDDO

90000 FORMAT(2I5,F17.6,D32.6,F32.4)
      RETURN
      END SUBROUTINE SETUP_EMIS


      SUBROUTINE SETUP_DP0         
!-------------------------------------------------------------------------------
!     Routine to calculate the diameter of average mass [m] for each mode.
!-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I

      !-------------------------------------------------------------------------
      ! Set lognormal parameters and activating fraction for each mode.
      !-------------------------------------------------------------------------
      DGN0(:)      = 0.0D+00
      SIG0(:)      = 0.0D+00
      LNSIG0(:)    = 0.0D+00
      DGN0_EMIS(:) = 0.0D+00
      SIG0_EMIS(:) = 0.0D+00
      KAPPAI(:)    = 0.0D+00
      DP0(:)       = 0.0D+00
      DP0_EMIS(:)  = 0.0D+00
      DO I=1, NMODES
        IF ( MODE_NAME(I) .EQ. 'AKK') DGN0(I) = DG_AKK     ! DGN0 for characteristic lognormal
        IF ( MODE_NAME(I) .EQ. 'ACC') THEN
          IF( NUMB_AKK_1 .EQ. 0 ) THEN
            IF ( ACTIVATION_COMPARISON .AND. ( MECH .EQ. 4 .OR. MECH .EQ. 8 ) ) THEN  ! No mode AKK 
              DGN0(I) = 0.0506736714D+00        ! For droplet activation test. 
              WRITE(*,*)'DGN0 for mode ACC set to special value of activation test.'
            ELSE
              DGN0(I) = SQRT( DG_AKK * DG_ACC ) ! No mode AKK; reduce DG_AKK
              WRITE(*,*)'DGN0 for mode ACC reduced in this mechanism.'
            ENDIF
          ELSE
            DGN0(I) = DG_ACC
          ENDIF
        ENDIF
        IF ( MODE_NAME(I) .EQ. 'DD1') THEN
          IF ( ACTIVATION_COMPARISON .AND. ( MECH .GE. 5 .AND. MECH .LE. 8 ) ) THEN  ! No mode DD2 or DS2 
            DGN0(I) = 0.4516304494D+00        ! For droplet activation test.
            WRITE(*,*)'DGN0 for mode DD1 set to special value of activation test.'
          ELSE
            DGN0(I) = DG_DD1
          ENDIF
        ENDIF
        IF ( MODE_NAME(I) .EQ. 'DD2') DGN0(I) = DG_DD2
        IF ( MODE_NAME(I) .EQ. 'DS1') DGN0(I) = DG_DS1
        IF ( MODE_NAME(I) .EQ. 'DS2') DGN0(I) = DG_DS2
        IF ( MODE_NAME(I) .EQ. 'SSA') DGN0(I) = DG_SSA
        IF ( MODE_NAME(I) .EQ. 'SSC') DGN0(I) = DG_SSC
        IF ( MODE_NAME(I) .EQ. 'SSS') THEN 
          IF ( ACTIVATION_COMPARISON .AND. ( MECH .EQ. 4 .OR. MECH .EQ. 8 ) ) THEN  ! No mode SSA or SSC 
            DGN0(I) = 0.3308445477D+00        ! For droplet activation test. 
            WRITE(*,*)'DGN0 for mode SSS set to special value for activation test.'
          ELSE
            DGN0(I) = DG_SSS
          ENDIF 
        ENDIF
        IF ( MODE_NAME(I) .EQ. 'OCC') DGN0(I) = DG_OCC
        IF ( MODE_NAME(I) .EQ. 'BC1') DGN0(I) = DG_BC1
        IF ( MODE_NAME(I) .EQ. 'BC2') DGN0(I) = DG_BC2
        IF ( MODE_NAME(I) .EQ. 'BC3') DGN0(I) = DG_BC3
        IF ( MODE_NAME(I) .EQ. 'DBC') DGN0(I) = DG_DBC
        IF ( MODE_NAME(I) .EQ. 'BOC') DGN0(I) = DG_BOC
        IF ( MODE_NAME(I) .EQ. 'BCS') DGN0(I) = DG_BCS
        IF ( MODE_NAME(I) .EQ. 'OCS') DGN0(I) = DG_OCS
        IF ( MODE_NAME(I) .EQ. 'MXX') DGN0(I) = DG_MXX
        IF ( MODE_NAME(I) .EQ. 'AKK') SIG0(I) = SG_AKK     ! SIG0 for characteristic lognormal
        IF ( MODE_NAME(I) .EQ. 'ACC') SIG0(I) = SG_ACC
        IF ( MODE_NAME(I) .EQ. 'DD1') SIG0(I) = SG_DD1
        IF ( MODE_NAME(I) .EQ. 'DD2') SIG0(I) = SG_DD2
        IF ( MODE_NAME(I) .EQ. 'DS1') SIG0(I) = SG_DS1
        IF ( MODE_NAME(I) .EQ. 'DS2') SIG0(I) = SG_DS2
        IF ( MODE_NAME(I) .EQ. 'SSA') SIG0(I) = SG_SSA
        IF ( MODE_NAME(I) .EQ. 'SSC') SIG0(I) = SG_SSC
        IF ( MODE_NAME(I) .EQ. 'SSS') SIG0(I) = SG_SSS
        IF ( MODE_NAME(I) .EQ. 'OCC') SIG0(I) = SG_OCC
        IF ( MODE_NAME(I) .EQ. 'BC1') SIG0(I) = SG_BC1
        IF ( MODE_NAME(I) .EQ. 'BC2') SIG0(I) = SG_BC2
        IF ( MODE_NAME(I) .EQ. 'BC3') SIG0(I) = SG_BC3
        IF ( MODE_NAME(I) .EQ. 'DBC') SIG0(I) = SG_DBC
        IF ( MODE_NAME(I) .EQ. 'BOC') SIG0(I) = SG_BOC
        IF ( MODE_NAME(I) .EQ. 'BCS') SIG0(I) = SG_BCS
        IF ( MODE_NAME(I) .EQ. 'OCS') SIG0(I) = SG_OCS
        IF ( MODE_NAME(I) .EQ. 'MXX') SIG0(I) = SG_MXX
        IF ( MODE_NAME(I) .EQ. 'AKK') DGN0_EMIS(I) = DG_AKK_EMIS     ! DGN0_EMIS for emissions lognormal
        IF ( MODE_NAME(I) .EQ. 'ACC') THEN
          DGN0_EMIS(I) = DG_ACC_EMIS
          IF( NUMB_AKK_1 .EQ. 0 ) THEN 
            DGN0_EMIS(I) = SQRT( DG_AKK_EMIS * DG_ACC_EMIS ) ! No mode AKK; reduce DG_AKK_EMIS
            WRITE(*,*)'DGN0_EMIS for mode ACC reduced in this mechanism.'
          ENDIF
        ENDIF
        IF ( MODE_NAME(I) .EQ. 'DD1') DGN0_EMIS(I) = DG_DD1_EMIS
        IF ( MODE_NAME(I) .EQ. 'DD2') DGN0_EMIS(I) = DG_DD2_EMIS
        IF ( MODE_NAME(I) .EQ. 'DS1') DGN0_EMIS(I) = DG_DS1_EMIS
        IF ( MODE_NAME(I) .EQ. 'DS2') DGN0_EMIS(I) = DG_DS2_EMIS
        IF ( MODE_NAME(I) .EQ. 'SSA') DGN0_EMIS(I) = DG_SSA_EMIS
        IF ( MODE_NAME(I) .EQ. 'SSC') DGN0_EMIS(I) = DG_SSC_EMIS
        IF ( MODE_NAME(I) .EQ. 'SSS') THEN 
          IF ( ACTIVATION_COMPARISON .AND. ( MECH .EQ. 4 .OR. MECH .EQ. 8 ) ) THEN  ! No mode SSA or SSC 
            DGN0_EMIS(I) = 0.3308445477D+00        ! For droplet activation test. 
            WRITE(*,*)'DGN0_EMIS for mode SSS set to special value for activation test.'
          ELSE
            DGN0_EMIS(I) = DG_SSS_EMIS
          ENDIF 
        ENDIF
        IF ( MODE_NAME(I) .EQ. 'OCC') DGN0_EMIS(I) = DG_OCC_EMIS
        IF ( MODE_NAME(I) .EQ. 'BC1') DGN0_EMIS(I) = DG_BC1_EMIS
        IF ( MODE_NAME(I) .EQ. 'BC2') DGN0_EMIS(I) = DG_BC2_EMIS
        IF ( MODE_NAME(I) .EQ. 'BC3') DGN0_EMIS(I) = DG_BC3_EMIS
        IF ( MODE_NAME(I) .EQ. 'DBC') DGN0_EMIS(I) = DG_DBC_EMIS
        IF ( MODE_NAME(I) .EQ. 'BOC') DGN0_EMIS(I) = DG_BOC_EMIS
        IF ( MODE_NAME(I) .EQ. 'BCS') DGN0_EMIS(I) = DG_BCS_EMIS
        IF ( MODE_NAME(I) .EQ. 'OCS') DGN0_EMIS(I) = DG_OCS_EMIS
        IF ( MODE_NAME(I) .EQ. 'MXX') DGN0_EMIS(I) = DG_MXX_EMIS
        IF ( MODE_NAME(I) .EQ. 'AKK') SIG0_EMIS(I) = SG_AKK_EMIS     ! SIG0_EMIS for emissions lognormal
        IF ( MODE_NAME(I) .EQ. 'ACC') SIG0_EMIS(I) = SG_ACC_EMIS
        IF ( MODE_NAME(I) .EQ. 'DD1') SIG0_EMIS(I) = SG_DD1_EMIS
        IF ( MODE_NAME(I) .EQ. 'DD2') SIG0_EMIS(I) = SG_DD2_EMIS
        IF ( MODE_NAME(I) .EQ. 'DS1') SIG0_EMIS(I) = SG_DS1_EMIS
        IF ( MODE_NAME(I) .EQ. 'DS2') SIG0_EMIS(I) = SG_DS2_EMIS
        IF ( MODE_NAME(I) .EQ. 'SSA') SIG0_EMIS(I) = SG_SSA_EMIS
        IF ( MODE_NAME(I) .EQ. 'SSC') SIG0_EMIS(I) = SG_SSC_EMIS
        IF ( MODE_NAME(I) .EQ. 'SSS') SIG0_EMIS(I) = SG_SSS_EMIS
        IF ( MODE_NAME(I) .EQ. 'OCC') SIG0_EMIS(I) = SG_OCC_EMIS
        IF ( MODE_NAME(I) .EQ. 'BC1') SIG0_EMIS(I) = SG_BC1_EMIS
        IF ( MODE_NAME(I) .EQ. 'BC2') SIG0_EMIS(I) = SG_BC2_EMIS
        IF ( MODE_NAME(I) .EQ. 'BC3') SIG0_EMIS(I) = SG_BC3_EMIS
        IF ( MODE_NAME(I) .EQ. 'DBC') SIG0_EMIS(I) = SG_DBC_EMIS
        IF ( MODE_NAME(I) .EQ. 'BOC') SIG0_EMIS(I) = SG_BOC_EMIS
        IF ( MODE_NAME(I) .EQ. 'BCS') SIG0_EMIS(I) = SG_BCS_EMIS
        IF ( MODE_NAME(I) .EQ. 'OCS') SIG0_EMIS(I) = SG_OCS_EMIS
        IF ( MODE_NAME(I) .EQ. 'MXX') SIG0_EMIS(I) = SG_MXX_EMIS
        IF ( MODE_NAME(I) .EQ. 'AKK') KAPPAI(I) = KAPPAI_AKK         ! KAPPAI - activating fraction
        IF ( MODE_NAME(I) .EQ. 'ACC') KAPPAI(I) = KAPPAI_ACC
        IF ( MODE_NAME(I) .EQ. 'DD1') KAPPAI(I) = KAPPAI_DD1
        IF ( MODE_NAME(I) .EQ. 'DD2') KAPPAI(I) = KAPPAI_DD2
        IF ( MODE_NAME(I) .EQ. 'DS1') KAPPAI(I) = KAPPAI_DS1
        IF ( MODE_NAME(I) .EQ. 'DS2') KAPPAI(I) = KAPPAI_DS2
        IF ( MODE_NAME(I) .EQ. 'SSA') KAPPAI(I) = KAPPAI_SSA
        IF ( MODE_NAME(I) .EQ. 'SSC') KAPPAI(I) = KAPPAI_SSC
        IF ( MODE_NAME(I) .EQ. 'SSS') KAPPAI(I) = KAPPAI_SSS
        IF ( MODE_NAME(I) .EQ. 'OCC') KAPPAI(I) = KAPPAI_OCC
        IF ( MODE_NAME(I) .EQ. 'BC1') KAPPAI(I) = KAPPAI_BC1
        IF ( MODE_NAME(I) .EQ. 'BC2') KAPPAI(I) = KAPPAI_BC2
        IF ( MODE_NAME(I) .EQ. 'BC3') KAPPAI(I) = KAPPAI_BC3
        IF ( MODE_NAME(I) .EQ. 'DBC') KAPPAI(I) = KAPPAI_DBC
        IF ( MODE_NAME(I) .EQ. 'BOC') KAPPAI(I) = KAPPAI_BOC
        IF ( MODE_NAME(I) .EQ. 'BCS') KAPPAI(I) = KAPPAI_BCS
        IF ( MODE_NAME(I) .EQ. 'OCS') KAPPAI(I) = KAPPAI_OCS
        IF ( MODE_NAME(I) .EQ. 'MXX') KAPPAI(I) = KAPPAI_MXX
        IF ( MODE_NAME(I) .EQ. 'AKK') DENSPI(I) = EMIS_DENS_SULF     ! DENSPI - default density for mode I
        IF ( MODE_NAME(I) .EQ. 'ACC') DENSPI(I) = EMIS_DENS_SULF
        IF ( MODE_NAME(I) .EQ. 'DD1') DENSPI(I) = EMIS_DENS_DUST
        IF ( MODE_NAME(I) .EQ. 'DD2') DENSPI(I) = EMIS_DENS_DUST
        IF ( MODE_NAME(I) .EQ. 'DS1') DENSPI(I) = EMIS_DENS_DUST
        IF ( MODE_NAME(I) .EQ. 'DS2') DENSPI(I) = EMIS_DENS_DUST
        IF ( MODE_NAME(I) .EQ. 'SSA') DENSPI(I) = EMIS_DENS_SEAS
        IF ( MODE_NAME(I) .EQ. 'SSC') DENSPI(I) = EMIS_DENS_SEAS
        IF ( MODE_NAME(I) .EQ. 'SSS') DENSPI(I) = EMIS_DENS_SEAS
        IF ( MODE_NAME(I) .EQ. 'OCC') DENSPI(I) = EMIS_DENS_OCAR
        IF ( MODE_NAME(I) .EQ. 'BC1') DENSPI(I) = EMIS_DENS_BCAR
        IF ( MODE_NAME(I) .EQ. 'BC2') DENSPI(I) = EMIS_DENS_BCAR
        IF ( MODE_NAME(I) .EQ. 'BC3') DENSPI(I) = EMIS_DENS_BCAR
        IF ( MODE_NAME(I) .EQ. 'DBC') DENSPI(I) = 0.5D+00 * ( EMIS_DENS_DUST + EMIS_DENS_BCAR )
        IF ( MODE_NAME(I) .EQ. 'BOC') DENSPI(I) = EMIS_DENS_BOCC
        IF ( MODE_NAME(I) .EQ. 'BCS') DENSPI(I) = 0.5D+00 * ( EMIS_DENS_BCAR + EMIS_DENS_SULF )
        IF ( MODE_NAME(I) .EQ. 'OCS') DENSPI(I) = 0.5D+00 * ( EMIS_DENS_OCAR + EMIS_DENS_SULF )
        IF ( MODE_NAME(I) .EQ. 'MXX') DENSPI(I) = 0.2D+00 * ( EMIS_DENS_SULF + EMIS_DENS_DUST 
     &                                                    +   EMIS_DENS_BCAR + EMIS_DENS_OCAR + EMIS_DENS_SEAS )
      ENDDO

      !---------------------------------------------------------------------------------------------------------------------
      ! Calculate the diameter of average mass for both the (default)
      ! characteristic lognormal and the emissions lognormal for each mode.
      ! 
      ! Dam [um] = ( diameter moment 3  /  diameter moment 0    )**(1/3)
      !          = ( dg**3 * sg**9                              )**(1/3)
      !          = ( dg**3 * [ exp( 0.5*(log(sigmag))**2 ) ]**9 )**(1/3)
      !          =   dg    * [ exp( 0.5*(log(sigmag))**2 ) ]**3
      !          =   dg    * [ exp( 1.5*(log(sigmag))**2 ) ]
      !
      ! Also store the natural logarithms of the geo. std. deviations. 
      !---------------------------------------------------------------------------------------------------------------------
      IF( WRITE_LOG ) WRITE(AUNIT1,'(/8A12/)') 'I','  MODE','DGN0 [um]','SIG0 [1]','DP0 [um]',
     &                                         'DGN0_E [um]','SIG0_E [1]','DP0_E [um]'
      DO I=1, NWEIGHTS
        DP0(I)      = 1.0D-06 * DGN0(I)      * EXP( 1.5D+00 * ( LOG(SIG0(I))      )**2 )  ! convert from [um] to [m]
        DP0_EMIS(I) = 1.0D-06 * DGN0_EMIS(I) * EXP( 1.5D+00 * ( LOG(SIG0_EMIS(I)) )**2 )  ! convert from [um] to [m]
        CONV_DPAM_TO_DGN(I) = 1.0D+00 / EXP( 1.5D+00 * ( LOG(SIG0_EMIS(I)) )**2 ) 
        LNSIG0(I) = LOG( SIG0(I) ) 
        IF( WRITE_LOG ) WRITE(AUNIT1,90000)I,MODE_NAME(I),DGN0(I),SIG0(I),DP0(I)*1.0D+06,
     &                                     DGN0_EMIS(I),SIG0_EMIS(I),DP0_EMIS(I)*1.0D+06
      ENDDO

      !---------------------------------------------------------------------------------------------------------------------
      ! Set densities and their reciprocals for each chemical component of any mode. 
      !---------------------------------------------------------------------------------------------------------------------
      DENS_COMP(1) = RHO_NH42SO4     ! [g/cm^3] sulfate 
      DENS_COMP(2) = EMIS_DENS_BCAR  ! [g/cm^3] BC
      DENS_COMP(3) = EMIS_DENS_OCAR  ! [g/cm^3] OC
      DENS_COMP(4) = EMIS_DENS_DUST  ! [g/cm^3] dust 
      DENS_COMP(5) = EMIS_DENS_SEAS  ! [g/cm^3] sea salt
      DENS_COMP(6) = RHO_NH42SO4     ! [g/cm^3] nitrate
      DENS_COMP(7) = RHO_NH42SO4     ! [g/cm^3] ammonium
      DENS_COMP(8) = RHO_H2O         ! [g/cm^3] water 
      DO I=1, 8
        RECIP_DENS_COMP(I) = 1.0D+00 / DENS_COMP(I)   ! [cm^3/g] sulfate 
        ! WRITE(*,'(I4,2F10.4)') I, DENS_COMP(I), RECIP_DENS_COMP(I)
      ENDDO

      !---------------------------------------------------------------------------------------------------------------------
      ! If doing comparison with the discrete pdf model, set all mode particle densities to the same value.
      !---------------------------------------------------------------------------------------------------------------------
      IF( DISCRETE_EVAL_OPTION ) THEN
        IF( WRITE_LOG ) WRITE(AUNIT1,'(/A,F7.3/)') 'Setting particle densities for all modes to (g/cm^3) ', DENSP
        DENSPI(:) = DENSP
        DENS_COMP(:) = DENSP 
        RECIP_DENS_COMP(:) = 1.0D+00 / DENSP 
      ENDIF

90000 FORMAT(I12,A12,6F12.6)
      RETURN
      END SUBROUTINE SETUP_DP0


      SUBROUTINE SETUP_COAG_TENSORS
!-------------------------------------------------------------------------------------------------------------------
!     Routine to define the g_ikl,q, the d_ikl, and the d_ij.
!
!     All elements GIKLQ(I,K,L,Q), DIKL(I,K,L), and DIJ(I,J) were checked through printouts available below. 
!
!     NM_SPC_NAME(I,:) contains the names of the mass species defined for mode I.
!-------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,J,K,L,Q,QQ, ntot, n
      LOGICAL, PARAMETER :: WRITE_TENSORS = .FALSE.

      IF ( WRITE_TENSORS ) THEN
        WRITE(AUNIT1,'(/A/)') 'CITABLE'
        DO I=1, NMODES
          WRITE(AUNIT1,90000) CITABLE(1:NMODES,I)
        ENDDO
        WRITE(AUNIT1,'(A)') '  '
      ENDIF

      GIKLQ(:,:,:,:) = 0
      DIKL(:,:,:) = 0
      DIJ(:,:) = 0

      !-------------------------------------------------------------------------------------------------------------
      ! The tensors g_ikl,q and d_ikl are symmetric in K and L.
      !
      ! GIKLQ is unity if coagulation of modes K and L produce mass of species Q
      !       in mode I, and zero otherwise.
      !
      ! DIKL is unity if coagulation of modes K and L produce particles 
      !      in mode I, and zero otherwise. 
      !      Neither mode K nor mode L can be mode I for a nonzero DIKL:
      !      all three modes I, K, L must be different modes. 
      !-------------------------------------------------------------------------------------------------------------
      DO I=1, NMODES
      DO K=1, NMODES
      DO L=K+1, NMODES                              ! Mode L is the same as mode K.
        IF ( CITABLE(K,L) .EQ. MODE_NAME(I) ) THEN  ! modes K and L produce mode I
          ! WRITE(36,*)'MODE_NAME(I) = ', MODE_NAME(I)
          IF ( I .NE. K  .AND. I .NE. L ) THEN      ! omit intramodal coagulation
            DIKL(I,K,L) = 1
            DIKL(I,L,K) = 1
          ENDIF
          DO Q=1, NM(I)                             ! loop over all mass species in mode I
            DO QQ=1, NMASS_SPCS                     ! loop over all principal mass species 
              !-----------------------------------------------------------------------------------------------------
              ! Compare the name of mass species Q in mode I with that of mass species QQ in mode K (or L).
              ! The inner loop is over all principal mass species since all species must be checked for 
              !   mode K (or L) for a potential match with species Q in mode I.
              !-----------------------------------------------------------------------------------------------------
              IF( NM_SPC_NAME(K,QQ) .EQ. NM_SPC_NAME(I,Q) ) THEN   ! mode K contains Q
                IF( I .NE. K ) GIKLQ(I,K,L,Q) = 1   ! I and K must be different modes
                IF( I .NE. K ) GIKLQ(I,L,K,Q) = 1   ! I and K must be different modes
              ENDIF
              IF( NM_SPC_NAME(L,QQ) .EQ. NM_SPC_NAME(I,Q) ) THEN   ! mode L contains Q
                IF( I .NE. L ) GIKLQ(I,K,L,Q) = 1   ! I and L must be different modes
                IF( I .NE. L ) GIKLQ(I,L,K,Q) = 1   ! I and L must be different modes
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      ENDDO
      ENDDO

      if (allocated(dikl_control)) then
        deallocate(dikl_control)
      end if
      allocate(dikl_control(count(dikl /= 0)))

      call initializeDiklControl(dikl_control, DIKL)

      if (allocated(giklq_control)) then
        do i = 1, nweights
          deallocate(giklq_control(i)%k)
          deallocate(giklq_control(i)%l)
          deallocate(giklq_control(i)%qq)
        end do
      else
        allocate(GIKLQ_control(NWEIGHTS))
      end if

      call initializeGiklqControl(GIKLQ_control, GIKLQ)

      !-------------------------------------------------------------------------------------------------------------
      ! The tensor d_ij is not symmetric in I,J.
      !
      ! DIJ(I,J) is unity if coagulation of mode I with mode J results
      !   in the removal of particles from mode I, and zero otherwise. 
      !-------------------------------------------------------------------------------------------------------------
      DO I=1, NMODES
      DO J=1, NMODES
        DO K=1, NMODES                               ! Find the product mode of the I-J coagulation.
          IF( I .EQ. J ) CYCLE                       ! Omit intramodal interactions: --> I .NE. J .
          IF( CITABLE(I,J) .EQ. MODE_NAME(K) ) THEN  ! I-particles and J-particles are lost; K-particles are formed.
            IF( I .NE. K ) DIJ(I,J) = 1              ! The K-particles are not I-particles (but may be J-particles),
          ENDIF                                      !   so I-particles are lost by this I-J interaction.
        ENDDO
      ENDDO
      ENDDO
      xDIJ = DIJ

      IF( .NOT. WRITE_TENSORS ) RETURN

      !-------------------------------------------------------------------------
      ! Write the g_ikl,q.
      !-------------------------------------------------------------------------
      DO I=1, NMODES
        WRITE(AUNIT1,'(/2A)') 'g_iklq for MODE ', MODE_NAME(I)
        DO Q=1, NM(I)
          WRITE(AUNIT1,'(/A,I3,3X,3A5/)') 'Q, NM_SPC_NAME(I,Q), MODE',
     &                        Q, NM_SPC_NAME(I,Q), '-->', MODE_NAME(I)
          IF ( SUM( GIKLQ(I,1:NMODES,1:NMODES,Q) ) .EQ. 0 ) CYCLE
          WRITE(AUNIT1,'(5X,16A5)') MODE_NAME(1:NMODES)
          DO K=1, NMODES
            WRITE(AUNIT1,'(A5,16I5)') MODE_NAME(K),GIKLQ(I,K,1:NMODES,Q)
          ENDDO
        ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      ! Write the d_ikl.
      !-------------------------------------------------------------------------
      DO I=1, NMODES
        WRITE(AUNIT1,'(/2A)') 'd_ikl for MODE ', MODE_NAME(I)
        IF ( SUM( DIKL(I,1:NMODES,1:NMODES) ) .EQ. 0 ) CYCLE
        WRITE(AUNIT1,'(5X,16A5)') MODE_NAME(1:NMODES)
        DO K=1, NMODES
          WRITE(AUNIT1,'(A5,16I5)') MODE_NAME(K),DIKL(I,K,1:NMODES)
        ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      ! Write the d_ij.
      !-------------------------------------------------------------------------
      WRITE(AUNIT1,'(/2A)') 'd_ij'
      WRITE(AUNIT1,'(5X,16A5)') MODE_NAME(1:NMODES)
      DO I=1, NMODES
        WRITE(AUNIT1,'(A5,16I5)') MODE_NAME(I),DIJ(I,1:NMODES)
      ENDDO

90000 FORMAT(14A4)
      RETURN

      contains

      subroutine initializeDiklControl(control, mask)
      type (dikl_type) :: control(:)
      integer, intent(in) :: mask(:,:,:)
      integer :: i, k, l, n

      n = 0
      do k = 1, NWEIGHTS
        do l = k+1, NWEIGHTS
          do i = 1, NWEIGHTS
            if (mask(i,k,l) /= 0) then
              n = n + 1
              control(n)%i = i
              control(n)%k = k
              control(n)%l = l
            end if          
          end do
        end do
      end do
      NDIKL = n
      end subroutine initializeDiklControl

      subroutine initializeGiklqControl(control, mask)
      type (GIKLQ_type) :: control(:)
      integer, intent(in) :: mask(:,:,:,:)

      integer :: i, q, k, l, n, nTotal

      do i = 1, NWEIGHTS
        ! 1) count contributing cases for mode i
        n = 0
        do q = 1, nm(i)
          do k = 1, nmodes
            do l = k+1, nmodes
              if (mask(i,k,l,q) /= 0) then
                if (i /= l) then
                  n = n + 1
                end if
                if (i /= k) then
                  n = n + 1
                end if
              end if
            end do
          end do
        end do
        nTotal = n
        ! 2) allocate nTotal entries
        control(i)%n = nTotal
        allocate(control(i)%k(nTotal))
        allocate(control(i)%l(nTotal))
        allocate(control(i)%qq(nTotal))

        ! 3) repeat sweep, but now assign k,l,qq
        n = 0
        do q = 1, nm(i)
          do k = 1, nmodes
            do l = k+1, nmodes
              if (mask(i,k,l,q) /= 0) then
                if (i /= l) then
                  n = n + 1
                  control(i)%k(n) = k
                  control(i)%l(n) = l
                  qq = prod_index(i,q)
                  control(i)%qq(n) = qq
                end if
                if (i /= k) then
                  n = n + 1
                  control(i)%k(n) = l
                  control(i)%l(n) = k
                  qq = prod_index(i,q)
                  control(i)%qq(n) = qq
                end if
              end if
            end do
          end do
        end do
        control(i)%n = n
      end do

      end subroutine initializeGiklqControl


      END SUBROUTINE SETUP_COAG_TENSORS


      SUBROUTINE SETUP_AERO_MASS_MAP
!-------------------------------------------------------------------------------
!     Defines a map giving the AERO locations of the NM(I) masses for each mode.
!
!     MASS_MAP(I,Q) is the location in AERO(:) of the Qth mass in mode I.
!
!     PROD_INDEX(I,Q) is the location in array PIQ(I,Q) of chemical species
!       CHEM_SPC_NAME(Q) for mode (quadrature weight) I.
!-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,Q,J

      MASS_MAP(:,:) = 0
      PROD_INDEX(:,:) = 0

      IF( WRITE_LOG ) WRITE(AUNIT1,'(/A/)')'I,J,Q,MODE_NAME(I),AERO_SPCS(J),PROD_INDEX(I,Q),MASS_MAP(I,Q)'

      DO I=1, NWEIGHTS                                       ! loop over modes (quadrature points)
        Q = 1
        DO J=1, NAEROBOX                                     ! loop over AERO species
          IF ( AERO_SPCS(J)(6:8) .NE. MODE_NAME(I) ) CYCLE   ! This AERO species is not for mode I.
          IF ( AERO_SPCS(J)(1:4) .EQ. 'MASS' ) THEN          ! This is a mass species for mode I. 
            IF ( AERO_SPCS(J)(10:13) .EQ. NM_SPC_NAME(I,Q) ) MASS_MAP(I,Q) = J  ! location of this species in AERO
            IF ( AERO_SPCS(J)(10:13) .EQ. 'SULF' ) THEN      ! This mass species is sulfate.
              PROD_INDEX(I,Q) = PROD_INDEX_SULF
            ENDIF
            IF ( AERO_SPCS(J)(10:13) .EQ. 'BCAR' ) THEN      ! This mass species is BC.
              PROD_INDEX(I,Q) = PROD_INDEX_BCAR
            ENDIF
            IF ( AERO_SPCS(J)(10:13) .EQ. 'OCAR' ) THEN      ! This mass species is OC.
              PROD_INDEX(I,Q) = PROD_INDEX_OCAR
            ENDIF
            IF ( AERO_SPCS(J)(10:13) .EQ. 'DUST' ) THEN      ! This mass species is dust.
              PROD_INDEX(I,Q) = PROD_INDEX_DUST
            ENDIF
            IF ( AERO_SPCS(J)(10:13) .EQ. 'SEAS' ) THEN      ! This mass species is sea salt.
              PROD_INDEX(I,Q) = PROD_INDEX_SEAS
            ENDIF
            IF( WRITE_LOG ) WRITE(AUNIT1,90000)I,J,Q,MODE_NAME(I),AERO_SPCS(J),PROD_INDEX(I,Q),MASS_MAP(I,Q)
            Q = Q + 1
            IF (Q .GT. NM(I) ) GOTO 10
          ENDIF
        ENDDO
10      CONTINUE
      ENDDO

90000 FORMAT(3I4,A8,4X,A16,5X,2I4)
      RETURN
      END SUBROUTINE SETUP_AERO_MASS_MAP


      SUBROUTINE ATMOSPHERE( ALT, SIGMA, DELTA, THETA )
!----------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
! AUTHOR  - Ralph Carmichael, Public Domain Aeronautical Software
! Reformatted for fixed-form Fortran 90 by D. Wright, 1-9-06.                 
! NOTE - If ALT > 86, the values returned will not be correct, but they will
!   not be too far removed from the correct values for density.
!   The reference document does not use the terms pressure and temperature
!   above 86 km.
!----------------------------------------------------------------------------
      USE CONSTANT, only : radius 
      IMPLICIT NONE
!============================================================================
!     A R G U M E N T S                                                     |
!============================================================================
      REAL,INTENT(IN)::  ALT    ! geometric ALTitude, km.
      REAL,INTENT(OUT):: SIGMA  ! density/sea-level standard density
      REAL,INTENT(OUT):: DELTA  ! pressure/sea-level standard pressure
      REAL,INTENT(OUT):: THETA  ! temperature/sea-level standard temperature
!============================================================================
!     L O C A L   C O N S T A N T S                                         |
!============================================================================
      REAL,PARAMETER:: REARTH = radius/1000.           ! radius of the Earth (km)
      REAL,PARAMETER:: GMR = 34.163195                 ! hydrostatic constant
      INTEGER,PARAMETER:: NTAB=8   ! number of entries in the defining tables
!============================================================================
!     L O C A L   V A R I A B L E S                                         |
!============================================================================
      INTEGER:: I,J,K                                              ! counters
      REAL:: H                                   ! geopotential ALTitude (km)
      REAL:: TGRAD, TBASE  ! temperature gradient and base temp of this layer
      REAL:: TLOCAL                                       ! local temperature
      REAL:: DELTAH                         ! height above base of this layer
!============================================================================
!     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
!============================================================================
      REAL,DIMENSION(NTAB),PARAMETER:: HTAB= 
     &  (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
      REAL,DIMENSION(NTAB),PARAMETER:: TTAB= 
     &    (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
      REAL,DIMENSION(NTAB),PARAMETER:: PTAB= 
     &     (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, 
     &                         6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
      REAL,DIMENSION(NTAB),PARAMETER:: GTAB= 
     &              (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)

      H=ALT*REARTH/(ALT+REARTH)  ! convert geometric to geopotential altitude

      I=1
      J=NTAB                                   ! setting up for binary search
      DO
        K=(I+J)/2                                          ! integer division
        IF (H < HTAB(K)) THEN
          J=K
        ELSE
          I=K
        END IF
        IF (J <= I+1) EXIT
      END DO

      TGRAD=GTAB(I)                                 ! I will be in 1...NTAB-1
      TBASE=TTAB(I)
      DELTAH=H-HTAB(I)
      TLOCAL=TBASE+TGRAD*DELTAH
      THETA=TLOCAL/TTAB(1)                                ! temperature ratio

      IF (TGRAD == 0.0) THEN                                 ! pressure ratio
        DELTA=PTAB(I)*EXP(-GMR*DELTAH/TBASE)
      ELSE
        DELTA=PTAB(I)*(TBASE/TLOCAL)**(GMR/TGRAD)
      END IF

      SIGMA=DELTA/THETA                                       ! density ratio

      RETURN
      END SUBROUTINE ATMOSPHERE 


      END MODULE AERO_SETUP
 
