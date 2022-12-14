      MODULE AERO_PARAM  
!-------------------------------------------------------------------------------------------------------------------------
!@sum     AEROSOL PARAMETERS AND VARIABLES THAT ARE INDEPENDENT OF CONFIGURATION. 
!@auth    Susanne Bauer/Doug Wright
!------------------------------------------------------------------------------------------------------------------------  
      IMPLICIT NONE
!-------------------------------------------------------------------------------------------------------------------------
!
!     GEOMETRIC AND SCIENTIFIC CONSTANTS; DERIVED CONVERSION FACTORS.
!
!-------------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: PI          = 3.141592653589793D+00
      REAL(8), PARAMETER :: PI6         = PI/6.0D+00
      REAL(8), PARAMETER :: AVO         = 6.0221367D+23            ! Avogadro's number     [#/mole]
      REAL(8), PARAMETER :: RGAS_SI     = 8.3145D+00               ! univers. gas constant [J/mol/K]
      REAL(8), PARAMETER :: MW_H2SO4    = 98.07948D+00             ! molar mass of H2SO4   [g/mole]
      REAL(8), PARAMETER :: MW_SO4      = 96.06360D+00             ! molar mass of SO4=    [g/mole]
      REAL(8), PARAMETER :: MW_NH3      = 17.03056D+00             ! molar mass of NH3     [g/mole]
      REAL(8), PARAMETER :: MW_H2O      = 18.01528D+00             ! molar mass of H2O     [g/mole]
      REAL(8), PARAMETER :: MW_NH42SO4  = 132.1406D+00             ! molar mass of NH42SO4 [g/mole]
      REAL(8), PARAMETER :: UGM3_NCM3   = 1.0D-12 * AVO / MW_SO4   ! [ugSO4/m^3] to [#/cm^3]
      REAL(8), PARAMETER :: CONVNH3     = RGAS_SI / MW_NH3         ! used in [ug NH3/m^3] to [ppmV]
      REAL(8), PARAMETER :: RHO_NH42SO4 = 1.77D+00                 ! density of dry (NH4)2SO4 [g/cm^3]
      REAL(8), PARAMETER :: RHO_H2SO4   = 1.84D+00                 ! density of pure H2SO4    [g/cm^3] - CRC
      REAL(8), PARAMETER :: RHO_H2O     = 1.00D+00                 ! density of pure H2O      [g/cm^3] - CRC
      REAL(8), PARAMETER :: DENSP       = 1.40D+00                 ! default ambient particle density [g/cm^3]
      !-------------------------------------------------------------------------------------------------------------------

      ! CONV_DP_TO_MASS converts Dp^3 [m^3] to particle mass [ug].
      ! CONV_MASS_TO_DP converts particle mass [ug] to Dp^3 [m^3].
      ! CONV_VOL_TO_DP_FAC converts particle volume [ug] to Dp^3 [m^3].
      ! The factor 1.0D+12 converts [m^3] to [cm^3] and [g] to [ug].
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: CONV_DP_TO_MASS    = 1.0D+12 * PI6 * DENSP
      REAL(8), PARAMETER :: CONV_MASS_TO_DP    = 1.0D+00 / CONV_DP_TO_MASS
      REAL(8), PARAMETER :: CONV_VOL_TO_DP_FAC = 1.0D-12 / PI6
      !-------------------------------------------------------------------------------------------------------------------
      ! Miniumum and maximum values of average mode diameters. [m]
      ! These are needed when number and/or mass concentrations are very small.
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DPMIN_GLOBAL =  0.001D-06   ! [m] -  1 nm
      REAL(8), PARAMETER :: DPMAX_GLOBAL = 20.000D-06   ! [m] - 20 um
!-------------------------------------------------------------------------------------------------------------------------
!
!     MODEL PARAMETERS AND VARIABLES THAT MAY NEED TO BE SET BY THE USER.
!
!-------------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: AUNIT1            = 90 ! logical unit # - log file of module
      INTEGER, PARAMETER :: AUNIT2            = 91 ! logical unit # - test of coag. coef.
      INTEGER, PARAMETER :: NEMIS_SPCS        = 10 ! number of emissions variables
      INTEGER, PARAMETER :: NDIAG_AERO        = 15 ! number of aerosol diagnostics collected
      INTEGER, PARAMETER :: KIJ_NDGS_SET      = 31 ! default value=81; if NO_MICROPHYSICS=.TRUE., set to 3 to save storage
      INTEGER, PARAMETER :: IMTR_METHOD       =  1 ! =1 no cut of pdf, =2 fixed-Dp cut, =3 variable-Dp cut as in CMAQ
      INTEGER, PARAMETER :: ACTIVATION_SCHEME =  2 ! =1 uses typical solubility only, =2 detailed multimodal activation
      INTEGER, PARAMETER :: UPDATE_KIJ        =  1 ! =0 use time-independent coagulation coefficients from lookup tables
                                                   ! =1 use time-  dependent coagulation coefficients from lookup tables
      LOGICAL, PARAMETER :: WRITE_LOG   = .FALSE.  ! WRITE MATRIX log to unit AUNIT1: default setting is .FALSE.
      LOGICAL, PARAMETER :: MASS_ADJ    = .TRUE.   ! enforce precise mass conservation: default setting is .TRUE.
      LOGICAL, PARAMETER :: CPU_STATS   = .FALSE.  ! timer for sections of the MATRIX code
      LOGICAL, PARAMETER :: UPDATE_DP   = .TRUE.   ! update particle diameters at each time step: default is .TRUE.
      LOGICAL, PARAMETER :: UPDATE_VDEP = .FALSE.  ! update particle diameters at each time step in global array
      LOGICAL, PARAMETER :: SET_INTERMODAL_TRANSFER  = .TRUE.  ! do AKK -> ACC transfer if mode AKK is defined
      LOGICAL, PARAMETER :: NO_MICROPHYSICS          = .FALSE. ! no-microphysics option
      LOGICAL, PARAMETER :: NO_MICROPHYSICS_W_THERMO = .TRUE.  ! do gas-particle partitioning when NO_MICROPHYSICS=.TRUE.
      LOGICAL, PARAMETER :: DO_NPF                   = .TRUE.  ! include secondary particle formation
      LOGICAL, PARAMETER :: DISCRETE_EVAL_OPTION     = .FALSE. ! for evaluation with results from the discrete pdf model
      LOGICAL, PARAMETER :: ACTIVATION_COMPARISON    = .FALSE. ! for comparison of aerosol activation w/ all 8 mechanisms
      !-------------------------------------------------------------------------------------------------------------------
      ! AQSO4RATE_MIN is the min, aqueous SO4 production rate for call of the activation routine. 
      ! The default value is 4.43D-11 [ug/m^3/s], equivalent to 1000 [molecules/cm3/h] = 1.6D-07 [ugSO4/m^3/h].
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: AQSO4RATE_MIN  = 4.43D-11 
      !-------------------------------------------------------------------------------------------------------------------
      ! The Maximum Inorganic Volume Fraction (MIVF) in modes DD1, DD2, BC1, and BC2.
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: MIVF_DDD = 0.05D+00   ! 
      REAL(8), PARAMETER :: MIVF_BC1 = 0.05D+00   !- These two are from MZJ 2002, "Analysis ..."
      REAL(8), PARAMETER :: MIVF_BC2 = 0.20D+00   !/
      !-------------------------------------------------------------------------------------------------------------------
      ! Scale factor for the geometric mean diameters of the emissions lognormals.
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: SCALE_EMIS_DIAM = 1.0D+00 
      !-------------------------------------------------------------------------------------------------------------------
      ! The minimum value of a number concentration or mass concentration leaving MATRIX.
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: MINCONC   = 1.0D-15   ! [ug/m^3] and [#/m^3]
      REAL(8), PARAMETER :: TINYNUMER = 1.0D-30   ! 
      REAL(8), PARAMETER :: TINYDENOM = 1.0D-30   !  
      !-------------------------------------------------------------------------------------------------------------------
      ! ZHEIGHT(I) is the mid-level height of model vertical layer I in [km].
      ! It is used in the calculation of ionization rates in the ion-ion
      ! recombination nucleation scheme, and in computing the pre-calculated
      ! factor in the condensational sink.
      !
      ! Set ZHEIGHT to global-average values typical of the vertical structure of the host GCM.
      !-------------------------------------------------------------------------------------------------------------------
       INTEGER, PARAMETER :: NLAYS = 42          ! number of model vertical layers
c        INTEGER, PARAMETER :: NLAYS = 103          ! number of model vertical layers
        REAL(8) :: ZHEIGHT(NLAYS)                 ! typical (CMAQ) mid-layer heights [km]
c M 20 model
c      DATA ZHEIGHT/ 0.007D+00,  0.024D+00,  0.052D+00,  0.100D+00,  0.210D+00,
c     &              0.390D+00,  0.640D+00,  0.950D+00,  1.300D+00,  1.740D+00,
c     &              2.260D+00,  2.810D+00,  3.390D+00,  4.000D+00,  4.700D+00,
c     &              5.400D+00,  6.200D+00,  7.200D+00,  8.400D+00,  10.00D+00,
c     &              12.40D+00 /
c F 40 model
      DATA ZHEIGHT/0.01D+00,0.02D+00,0.04D+00,0.06D+00,0.09D+00,1.2D+00,1.5D+00,2.D+00,2.4D+00,3.D+00
     &             ,3.5D+00,4.D+00,6.7D+00,7.4D+00,8.1D+00,8.5D+00,9.D+00,10.D+00,11.D+00,12.D+00
     &             ,13.D+00,14.D+00,15.D+00,16.D+00,18.D+00,19.D+00,21.D+00,24.D+00,28.D+00,32.D+00
     &             ,36.D+00,40.D+00,44.D+00,48.D+00,53.D+00,58.D+00,62.D+00,66.D+00,72.D+00,80.D+00
     &             ,85.D+00,90.D+00/
c F 102 model
c      DATA ZHEIGHT/0.0065, 0.19, 0.33, 0.46, 0.60, 0.74, 0.88, 1.03, 1.17, 1.32,
c     &               1.47, 1.62, 1.79, 1.93, 2.09, 2.25, 2.42, 2.59, 2.77, 2.96,
c     &               3.19, 3.46, 3.73, 4.02, 4.33, 4.64, 4.96, 5.29, 5.63, 5.96,
c     &               6.30, 6.64, 6.97, 7.30, 7.63, 7.96, 8.27, 8.59, 8.91, 9.23,
c     &               9.55, 9.87,10.19,10.52,10.85,11.19,11.53,11.89,12.26,12.64,
c     &               13.02,13.43,13.85,14.26,14.68,15.13,15.62,16.11,16.61,17.11,
c     &               17.62,18.18,18.75,19.32,19.94,20.58,21.22,21.95,22.69,23.44,
c     &               24.30,25.17,26.03,27.03,28.22,29.41,30.56,31.73,32.92,34.14,
c     &               35.37,36.63,37.88,39.15,40.43,41.73,43.03,44.41,45.81,47.35, 
c     &               48.99,50.76,52.80,55.33,58.61,62.42,66.16,69.79,73.34,76.94, 
c     &               80.81, 85.58, 90./
!-------------------------------------------------------------------------------------------------------------------------
!    Default values for the box model - do not change these. 
!-------------------------------------------------------------------------------------------------------------------------
!     INTEGER, PARAMETER :: NLAYS = 21          ! number of model vertical layers
!     REAL(8) :: ZHEIGHT(NLAYS)                 ! typical (CMAQ) mid-layer heights [km]
!     DATA ZHEIGHT/ 0.007D+00,  0.024D+00,  0.052D+00,  0.100D+00,  0.210D+00,
!    &              0.390D+00,  0.640D+00,  0.950D+00,  1.300D+00,  1.740D+00,
!    &              2.260D+00,  2.810D+00,  3.390D+00,  4.000D+00,  4.700D+00,
!    &              5.400D+00,  6.200D+00,  7.200D+00,  8.400D+00,  10.00D+00,
!    &              12.40D+00 /
!-------------------------------------------------------------------------------------------------------------------------
      ! Characteristic lognormal parameters for each mode: DG_XXX[um],SG_XXX[1]. 
      ! 
      ! The Dg values from Easter et al., 2004 (E04) are the "diagnosed number" values, not the "emitted" values. 
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DG_AKK = 0.026D+00      ! E04, Table 2, Aitken       mode
      REAL(8), PARAMETER :: DG_ACC = 0.110D+00      ! E04, Table 2, accumulation mode         
      !-------------------------------------------------------------------------------------------------------------------
      ! DLW: 021507: These dust and sea salt Dg values were calculated approximate emissions sizes used at GISS. 
      !
      ! The number mean diameters for the four GISS dust size classes were 0.46, 2.94, 5.88, and 11.77 micrometers.
      ! Size classes 1 and 2 were averaged with a 10:1 ratio, and the average converted to a lognormal geometric
      ! mean diameter for an assumed geometric standard deviation of 1.8. 
      ! Likewise, size classes 3 and 4 were averaged with a 10:1 ratio, and the average converted to a lognormal 
      ! geometric mean diameter for an assumed geometric standard deviation of 1.8.  
      !
      ! The number mean diameters for the two GISS sea salt size classes were 0.44 and 5.0 micrometers, and were
      ! converted to lognormal geometric mean diameters for an assumed geometric standard deviation of 1.8 for the
      ! smaller (accumulation) size class and a standard deviation of 2.0 for the larger (coarse) size class. 
       !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DG_DD1 = 0.580D+00 *2.     ! set to match GISS dust emissions for average of sizes 1 & 2         
      REAL(8), PARAMETER :: DG_DD2 = 1.000D+00 *2.     ! set to match GISS dust emissions for average of sizes 3 & 4         
      REAL(8), PARAMETER :: DG_DS1 = 0.580D+00 *2.     ! set to match GISS dust emissions for average of sizes 1 & 2         
      REAL(8), PARAMETER :: DG_DS2 = 1.00D+00 *2.     ! set to match GISS dust emissions for average of sizes 3 & 4         
      REAL(8), PARAMETER :: DG_SSA = 0.440D+00 *2.     ! set to match GISS sea salt emissions         
      REAL(8), PARAMETER :: DG_SSC = 1.D+00 *2.     ! set to match GISS sea salt emissions         
      REAL(8), PARAMETER :: DG_SSS = 0.690D+00 *2.     ! 10:1 average of modes SSA and SSC
      !-------------------------------------------------------------------------------------------------------------------
      ! DLW: 021507: End of emissions sizes used at GISS. 
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DG_OCC = 0.050D+00      ! geo. avg. of E04, Table 2, Aitken and accumulation modes         
      REAL(8), PARAMETER :: DG_BC1 = 0.050D+00      ! geo. avg. of E04, Table 2, Aitken and accumulation modes         
      REAL(8), PARAMETER :: DG_BC2 = 0.100D+00      ! 1.0164 * geo. avg. of E04, Table 2, AKK and ACC modes, w/  5% shell         
      REAL(8), PARAMETER :: DG_BC3 = 0.100D+00      ! 1.0627 * geo. avg. of E04, Table 2, AKK and ACC modes, w/ 20% shell       
      REAL(8), PARAMETER :: DG_DBC = 0.330D+00      ! geo. avg. of E04, Table 2, accumulation and coarse dust modes
      REAL(8), PARAMETER :: DG_BOC = 0.100D+00      ! assuming additive volumes for BC1 and OCC
      REAL(8), PARAMETER :: DG_BCS = 0.070D+00      ! assuming add. vol. for BC1 and AKK (ACC) and greater weight for AKK
      REAL(8), PARAMETER :: DG_OCS = 0.070D+00      ! assuming add. vol. for BC1 and AKK (ACC) and greater weight for AKK
      REAL(8), PARAMETER :: DG_MXX = 0.300D+00      ! value is midrange considering all modes
      REAL(8), PARAMETER :: SG_AKK = 1.600D+00      ! E04, Table 2, Aitken       mode
      REAL(8), PARAMETER :: SG_ACC = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_DD1 = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_DD2 = 1.800D+00      ! E04, Table 2, coarse       mode         
      REAL(8), PARAMETER :: SG_DS1 = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_DS2 = 1.800D+00      ! E04, Table 2, coarse       mode         
      REAL(8), PARAMETER :: SG_SSA = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_SSC = 2.000D+00      ! E04, Table 2, coarse       mode         
      REAL(8), PARAMETER :: SG_SSS = 2.000D+00      ! same as SSC
      REAL(8), PARAMETER :: SG_OCC = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_BC1 = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_BC2 = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_BC3 = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_DBC = 1.800D+00      ! same as parent modes
      REAL(8), PARAMETER :: SG_BOC = 1.800D+00      ! same as parent modes
      REAL(8), PARAMETER :: SG_BCS = 1.800D+00      ! same as parent modes
      REAL(8), PARAMETER :: SG_OCS = 1.800D+00      ! same as parent modes
      REAL(8), PARAMETER :: SG_MXX = 2.000D+00      ! likely a broad mode
      !-------------------------------------------------------------------------------------------------------------------
      ! Lognormal parameters for emissions into each mode: DG_XXX_EMIS[um],SG_XXX_EMIS[1]. 
      ! 
      ! These are used to convert mass emission rates to number emission rates.
      ! The Dg values from Easter et al., 2004 (E04) are the "emitted" values, not the "diagnosed number" values. 
      ! All modes are assigned a value, even if they do not receive primary particles. 
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DG_AKK_EMIS = 0.013D+00      ! E04, Table 2, Aitken       mode
      REAL(8), PARAMETER :: DG_ACC_EMIS = 0.068D+00      ! E04, Table 2, accumulation mode         
      !-------------------------------------------------------------------------------------------------------------------
      ! DLW: 021507: These dust and sea salt Dg values were calculated approximate emissions sizes used at GISS. 
      !
      ! The number mean diameters for the four GISS dust size classes were 0.46, 2.94, 5.88, and 11.77 micrometers.
      ! Size classes 1 and 2 were averaged with a 10:1 ratio, and the average converted to a lognormal geometric
      ! mean diameter for an assumed geometric standard deviation of 1.8. 
      ! Likewise, size classes 3 and 4 were averaged with a 10:1 ratio, and the average converted to a lognormal 
      ! geometric mean diameter for an assumed geometric standard deviation of 1.8.  
      !
      ! The number mean diameters for the two GISS sea salt size classes were 0.44 and 5.0 micrometers, and were
      ! converted to lognormal geometric mean diameters for an assumed geometric standard deviation of 1.8 for the
      ! smaller (accumulation) size class and a standard deviation of 2.0 for the larger (coarse) size class. 
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DG_DD1_EMIS = 0.580D+00 *2.     ! set to match GISS dust emissions for average of sizes 1 & 2         
      REAL(8), PARAMETER :: DG_DD2_EMIS = 1.000D+00 *2.     ! set to match GISS dust emissions for average of sizes 3 & 4         
      REAL(8), PARAMETER :: DG_DS1_EMIS = 0.580D+00 *2.     ! set to match GISS dust emissions for average of sizes 1 & 2         
      REAL(8), PARAMETER :: DG_DS2_EMIS = 1.00D+00 *2.     ! set to match GISS dust emissions for average of sizes 3 & 4         
      REAL(8), PARAMETER :: DG_SSA_EMIS = 0.440D+00 *2     ! set to match GISS sea salt emissions         
      REAL(8), PARAMETER :: DG_SSC_EMIS = 1.000D+00 *2.     ! set to match GISS sea salt emissions         
      REAL(8), PARAMETER :: DG_SSS_EMIS = 0.690D+00 *2.     ! 10:1 average of modes SSA and SSC
      !-------------------------------------------------------------------------------------------------------------------
      ! DLW: 021507: End of emissions sizes used at GISS. 
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DG_OCC_EMIS = 0.050D+00      ! geometric average of the AKK and ACC values in E04 
      REAL(8), PARAMETER :: DG_BC1_EMIS = 0.050D+00      ! geometric average of the AKK and ACC values in E04 
      REAL(8), PARAMETER :: DG_BC2_EMIS = 0.100D+00      !   currently no emissions into this mode 
      REAL(8), PARAMETER :: DG_BC3_EMIS = 0.100D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: DG_DBC_EMIS = 0.300D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: DG_BOC_EMIS = 0.100D+00      ! assuming additive volumes for BC1 and OCC
      REAL(8), PARAMETER :: DG_BCS_EMIS = 0.140D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: DG_OCS_EMIS = 0.140D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: DG_MXX_EMIS = 0.500D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: SG_AKK_EMIS = 1.600D+00      ! E04, Table 2, Aitken       mode
      REAL(8), PARAMETER :: SG_ACC_EMIS = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_DD1_EMIS = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_DD2_EMIS = 1.800D+00      ! E04, Table 2, coarse       mode         
      REAL(8), PARAMETER :: SG_DS1_EMIS = 1.800D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: SG_DS2_EMIS = 1.800D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: SG_SSA_EMIS = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_SSC_EMIS = 2.000D+00      ! E04, Table 2, coarse       mode         
      REAL(8), PARAMETER :: SG_SSS_EMIS = 2.000D+00      ! same as SSC
      REAL(8), PARAMETER :: SG_OCC_EMIS = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_BC1_EMIS = 1.800D+00      ! E04, Table 2, accumulation mode         
      REAL(8), PARAMETER :: SG_BC2_EMIS = 1.800D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: SG_BC3_EMIS = 1.800D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: SG_DBC_EMIS = 1.800D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: SG_BOC_EMIS = 1.800D+00      ! same as BC1 and OCC modes 
      REAL(8), PARAMETER :: SG_BCS_EMIS = 1.800D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: SG_OCS_EMIS = 1.800D+00      !   currently no emissions into this mode
      REAL(8), PARAMETER :: SG_MXX_EMIS = 2.000D+00      !   currently no emissions into this mode
      !-------------------------------------------------------------------------------------------------------------------
      ! KAPPAI_XXX is the activating fraction for mode XXX.
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: KAPPAI_AKK = 0.0D+00
      REAL(8), PARAMETER :: KAPPAI_ACC = 1.0D+00  
      REAL(8), PARAMETER :: KAPPAI_DD1 = 0.0D+00
      REAL(8), PARAMETER :: KAPPAI_DD2 = 0.0D+00
      REAL(8), PARAMETER :: KAPPAI_DS1 = 1.0D+00
      REAL(8), PARAMETER :: KAPPAI_DS2 = 1.0D+00
      REAL(8), PARAMETER :: KAPPAI_SSA = 1.0D+00
      REAL(8), PARAMETER :: KAPPAI_SSC = 1.0D+00
      REAL(8), PARAMETER :: KAPPAI_SSS = 1.0D+00
      REAL(8), PARAMETER :: KAPPAI_OCC = 0.7D+00
      REAL(8), PARAMETER :: KAPPAI_BC1 = 0.0D+00
      REAL(8), PARAMETER :: KAPPAI_BC2 = 1.0D+00
      REAL(8), PARAMETER :: KAPPAI_BC3 = 1.0D+00
      REAL(8), PARAMETER :: KAPPAI_DBC = 0.0D+00
      REAL(8), PARAMETER :: KAPPAI_BOC = 0.5D+00
      REAL(8), PARAMETER :: KAPPAI_BCS = 1.0D+00
      REAL(8), PARAMETER :: KAPPAI_OCS = 1.0D+00
      REAL(8), PARAMETER :: KAPPAI_MXX = 1.0D+00
!-------------------------------------------------------------------------------------------------------------------
! Solubility per mode
!-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: SOLU_AKK = 1.0D+00
      REAL(8), PARAMETER :: SOLU_ACC = 1.0D+00  
      REAL(8), PARAMETER :: SOLU_DD1 = 0.5D+00
      REAL(8), PARAMETER :: SOLU_DD2 = 0.5D+00
      REAL(8), PARAMETER :: SOLU_DS1 = 1.0D+00
      REAL(8), PARAMETER :: SOLU_DS2 = 1.0D+00
      REAL(8), PARAMETER :: SOLU_SSA = 1.0D+00
      REAL(8), PARAMETER :: SOLU_SSC = 1.0D+00
      REAL(8), PARAMETER :: SOLU_SSS = 1.0D+00
      REAL(8), PARAMETER :: SOLU_OCC = 0.4D+00
      REAL(8), PARAMETER :: SOLU_BC1 = 0.4D+00
      REAL(8), PARAMETER :: SOLU_BC2 = 0.8D+00
      REAL(8), PARAMETER :: SOLU_BC3 = 1.0D+00
      REAL(8), PARAMETER :: SOLU_DBC = 0.0D+00
      REAL(8), PARAMETER :: SOLU_BOC = 0.6D+00
      REAL(8), PARAMETER :: SOLU_BCS = 1.0D+00
      REAL(8), PARAMETER :: SOLU_OCS = 1.0D+00
      REAL(8), PARAMETER :: SOLU_MXX = 1.0D+00
!-------------------------------------------------------------------------------------------------------------------------
!
!     MODEL PARAMETERS AND VARIABLES THAT PROBABLY DO NOT NEED TO BE CHANGED.
!
!-------------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: NGASES     = 3      ! number of gas-phase species
      INTEGER, PARAMETER :: NMASS_SPCS = 5      ! total number of mass species
      INTEGER, PARAMETER :: GAS_H2SO4  = 1      !\.
      INTEGER, PARAMETER :: GAS_HNO3   = 2      !-indices in the GAS array
      INTEGER, PARAMETER :: GAS_NH3    = 3      !/
      INTEGER, PARAMETER :: PROD_INDEX_SULF = 1 ! SULF index in PROD_INDEX(:,:)
      INTEGER, PARAMETER :: PROD_INDEX_BCAR = 2 ! BCAR index in PROD_INDEX(:,:)
      INTEGER, PARAMETER :: PROD_INDEX_OCAR = 3 ! OCAR index in PROD_INDEX(:,:)
      INTEGER, PARAMETER :: PROD_INDEX_DUST = 4 ! DUST index in PROD_INDEX(:,:)
      INTEGER, PARAMETER :: PROD_INDEX_SEAS = 5 ! SEAS index in PROD_INDEX(:,:)
      !-------------------------------------------------------------------------------------------------------------------
      ! EMIS_DENS_XXXX is the dry particle density of emitted species XXXX.
      !
      ! These are used only for deriving number emission rates from mass emission rates.
      !
      ! For sulfate, the emissions are treated as pure dry ammonium sulfate.
      ! For sea salt, the emissions are treated as pure dry sodium chloride.
      !
      ! For the sulfate density, volume is converted to ammonium sulfate mass using
      ! the density of ammonium sulfate, which is then converted to sulfate (only) mass.
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: EMIS_DENS_SULF = 1.77D+00   ! [gSO4/cm^3] - NH42SO4
     &                                     * MW_SO4 / MW_NH42SO4
      REAL(8), PARAMETER :: EMIS_DENS_BCAR = 1.70D+00   ! [g/cm^3] - Ghan et al. (2001) - MIRAGE
      REAL(8), PARAMETER :: EMIS_DENS_OCAR = 1.00D+00   ! [g/cm^3] - Ghan et al. (2001) - MIRAGE
      REAL(8), PARAMETER :: EMIS_DENS_DUST = 2.60D+00   ! [g/cm^3] - Ghan et al. (2001) - MIRAGE
      REAL(8), PARAMETER :: EMIS_DENS_SEAS = 2.165D+00  ! [g/cm^3] - NaCl
      REAL(8), PARAMETER :: EMIS_DENS_BOCC = 0.50D+00   ! [g/cm^3] - average
     &                                     * ( EMIS_DENS_BCAR + EMIS_DENS_OCAR ) 
      REAL, DIMENSION(NEMIS_SPCS) :: EMIS_DENS = (/  EMIS_DENS_SULF,
     &               EMIS_DENS_SULF, EMIS_DENS_BCAR, EMIS_DENS_OCAR,
     &               EMIS_DENS_DUST, EMIS_DENS_SEAS, EMIS_DENS_SEAS,
     &               EMIS_DENS_BOCC, EMIS_DENS_BOCC, EMIS_DENS_DUST /)
      !-------------------------------------------------------------------------------------------------------------------
      ! The aerosol chemical species are SO4, BC, OC, mineral dust, and sea salt.
      ! Nitrate, ammonium and water are not included here.
      !-------------------------------------------------------------------------------------------------------------------
      CHARACTER(LEN=4) :: CHEM_SPC_NAME(NMASS_SPCS)
     &                 = (/'SULF','BCAR','OCAR','DUST','SEAS'/)
      !-------------------------------------------------------------------------------------------------------------------
      ! The Maximum Inorganic Mass Ratio (MIMR) in modes DD1, DD2, BC1, and BC2.
      ! 
      ! The above MIVF values are converted to MIMR values for computational efficiency.
      ! The volume of inorganic coating is converted to mass using the default ambient aerosol density.
      !-------------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: MIMR_DDD = ( DENSP / EMIS_DENS_DUST )
     &                               * MIVF_DDD / ( 1.0D+00 - MIVF_DDD )
      REAL(8), PARAMETER :: MIMR_BC1 = ( DENSP / EMIS_DENS_BCAR )
     &                               * MIVF_BC1 / ( 1.0D+00 - MIVF_BC1 )
      REAL(8), PARAMETER :: MIMR_BC2 = ( DENSP / EMIS_DENS_BCAR )
     &                               * MIVF_BC2 / ( 1.0D+00 - MIVF_BC2 )
!-------------------------------------------------------------------------------------------------------------------------
!
!     VARIABLES HELD IN THE MANNER OF A COMMON BLOCK.
!
!-------------------------------------------------------------------------------------------------------------------------
      INTEGER, SAVE :: IXXX, IYYY, ILAY    ! current grid cell indices 
      LOGICAL, SAVE :: INCLUDE_BC3         ! true if mechanism includes mode BC3; false otherwise
!-------------------------------------------------------------------------------------------------------------------------
!     Aerosol mode names (and numbers) that might appear in one or more mechanisms.
!     These mode number only pertain to this set of all possible modes, and are not the mode numbers
!     used for any specific mechanism.
!-------------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: NMODES_MAX=18 
      CHARACTER(LEN=3) :: MNAME(NMODES_MAX)
      ! Mode #     1     2     3     4     5     6     7     8     9    ! # to identify the mode in MODES1, etc. below. 
      DATA MNAME/'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','SSS',
     &           'OCC','BC1','BC2','BC3','OCS','DBC','BOC','BCS','MXX'/
      ! Mode #     10    11    12    13    14    15    16    17    18   ! # to identify the mode in MODES1, etc. below. 
!-------------------------------------------------------------------------------------------------------------------------
!     Aerosol species defined for each mode in each mechanism.
!-------------------------------------------------------------------------------------------------------------------------
      INTEGER, SAVE :: MSPCS(NMASS_SPCS,NMODES_MAX) 
      DATA MSPCS(1,1:NMODES_MAX)/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/   ! SULF: =0 no sulfate, =1 has sulfate
      DATA MSPCS(2,1:NMODES_MAX)/0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1/   ! BCAR: =0 no BC     , =1 has BC
      DATA MSPCS(3,1:NMODES_MAX)/0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1/   ! OCAR: =0 no OC     , =1 has OC
      DATA MSPCS(4,1:NMODES_MAX)/0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,1/   ! DUST: =0 no dust   , =1 has dust     
      DATA MSPCS(5,1:NMODES_MAX)/0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1/   ! SEAS: =0 no seasalt, =1 has seasalt
!-------------------------------------------------------------------------------------------------------------------------
!     Aerosol modes used for each mechanism.
!-------------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: NM1=16,NM2=16,NM3=13,NM4=10
      INTEGER, PARAMETER :: NM5=14,NM6=14,NM7=11,NM8=8  
      INTEGER :: MODES1(NM1),MODES2(NM2),MODES3(NM3),MODES4(NM4)
      INTEGER :: MODES5(NM5),MODES6(NM6),MODES7(NM7),MODES8(NM8)
      DATA MODES1/ 1, 2, 3, 4, 5, 6, 7, 8,10,11,12,13,15,16,17,18/
      DATA MODES2/ 1, 2, 3, 4, 5, 6, 7, 8,10,11,12,14,15,16,17,18/
      DATA MODES3/ 1, 2, 3, 4, 5, 6, 7, 8,10,11,12,16,18/
      DATA MODES4/ 2, 3, 4, 5, 6, 9,10,11,12,18/
      DATA MODES5/ 1, 2, 3, 4, 7, 8,10,11,12,13,15,16,17,18/
      DATA MODES6/ 1, 2, 3, 4, 7, 8,10,11,12,14,15,16,17,18/
      DATA MODES7/ 1, 2, 3, 4, 7, 8,10,11,12,16,18/
      DATA MODES8/ 2, 3, 4, 9,10,11,12,18/
!-------------------------------------------------------------------------------------------------------------------------
!     Indices of the AERO array. There are 78 possible indices.
!-------------------------------------------------------------------------------------------------------------------------
      INTEGER       :: MASS_NO3=1, MASS_NH4=2, MASS_H2O=3 
      INTEGER, SAVE :: NUMB_AKK_1, NUMB_AKK_2, MASS_AKK_SULF,  
     &                 NUMB_ACC_1, NUMB_ACC_2, MASS_ACC_SULF,
     &                 NUMB_DD1_1, NUMB_DD1_2, MASS_DD1_SULF,                               MASS_DD1_DUST, 
     &                 NUMB_DS1_1, NUMB_DS1_2, MASS_DS1_SULF,                               MASS_DS1_DUST, 
     &                 NUMB_DD2_1, NUMB_DD2_2, MASS_DD2_SULF,                               MASS_DD2_DUST, 
     &                 NUMB_DS2_1, NUMB_DS2_2, MASS_DS2_SULF,                               MASS_DS2_DUST, 
     &                 NUMB_SSA_1, NUMB_SSA_2, MASS_SSA_SULF,                                              MASS_SSA_SEAS, 
     &                 NUMB_SSC_1, NUMB_SSC_2, MASS_SSC_SULF,                                              MASS_SSC_SEAS,
     &                 NUMB_SSS_1, NUMB_SSS_2, MASS_SSS_SULF,                                              MASS_SSS_SEAS,
     &                 NUMB_OCC_1, NUMB_OCC_2, MASS_OCC_SULF,                MASS_OCC_OCAR,
     &                 NUMB_BC1_1, NUMB_BC1_2, MASS_BC1_SULF, MASS_BC1_BCAR,
     &                 NUMB_BC2_1, NUMB_BC2_2, MASS_BC2_SULF, MASS_BC2_BCAR,
     &                 NUMB_BC3_1, NUMB_BC3_2, MASS_BC3_SULF, MASS_BC3_BCAR,
     &                 NUMB_OCS_1, NUMB_OCS_2, MASS_OCS_SULF,                MASS_OCS_OCAR,
     &                 NUMB_DBC_1, NUMB_DBC_2, MASS_DBC_SULF, MASS_DBC_BCAR,                MASS_DBC_DUST,
     &                 NUMB_BOC_1, NUMB_BOC_2, MASS_BOC_SULF, MASS_BOC_BCAR, MASS_BOC_OCAR,
     &                 NUMB_BCS_1, NUMB_BCS_2, MASS_BCS_SULF, MASS_BCS_BCAR,
     &                 NUMB_MXX_1, NUMB_MXX_2, MASS_MXX_SULF, MASS_MXX_BCAR, MASS_MXX_OCAR, MASS_MXX_DUST, MASS_MXX_SEAS

      END MODULE AERO_PARAM  
