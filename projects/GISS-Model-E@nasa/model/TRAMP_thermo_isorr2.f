       SUBROUTINE AERO_THERMO(ASO4,ANO3,ANH4,AH2O,GNH3,GHNO3,TOT_DUST,
     &                       SEAS,SSH2O,TK,RH,PRES,RHD,RHC)
!@sum
!@+     This routine sets up for and calls the thermodynamic module for aerosol
!@+     gas-particle partitioning.
!@+
!@+     This version of AERO_THERMO is for use with the ISORROPIA thermodynamic module.
!@auth Susanne Bauer/Doug Wright

      USE AERO_PARAM, ONLY: WRITE_LOG, AUNIT1
      USE AERO_PARAM, only: IXXX, IYYY, ILAY
      IMPLICIT NONE

      !------------------------------------------------------------------------------------------------------
      ! Arguments.
      !------------------------------------------------------------------------------------------------------
      REAL(8), INTENT(INOUT) :: ASO4      ! aerosol sulfate       [ug/m^3]
      REAL(8), INTENT(INOUT) :: ANO3      ! aerosol nitrate       [ug/m^3]
      REAL(8), INTENT(INOUT) :: ANH4      ! aerosol ammonium      [ug/m^3]
      REAL(8), INTENT(INOUT) :: AH2O      ! aerosol water         [ug/m^3]
      REAL(8), INTENT(INOUT) :: GNH3      ! gas-phase ammonia     [ugNH4/m^3] as ammonium (MW)
      REAL(8), INTENT(INOUT) :: GHNO3     ! gas-phase nitric acid [ugNO3/m^3] as nitrate  (MW)
      REAL(8), INTENT(IN)    :: TOT_DUST  ! total dust(sol+insol) [ug/m^3]
      REAL(8), INTENT(IN)    :: SEAS      ! sea salt (NaCl)       [ug/m^3]
      REAL(8), INTENT(OUT)   :: SSH2O     ! sea salt assoc. H2O   [ug/m^3]
      REAL(8), INTENT(IN)    :: TK        ! absolute temperature  [K]          
      REAL(8), INTENT(IN)    :: RH        ! relative humidity     [0-1]
      REAL(8), INTENT(IN)    :: PRES      ! ambient pressure      [Pa]  
      REAL(8), INTENT(OUT)   :: RHD       ! RH of deliquescence   [0-1]
      REAL(8), INTENT(OUT)   :: RHC       ! RH of crystallization [0-1]

      !------------------------------------------------------------------------------------------------------
      ! Input to ISOROPIA.
      !------------------------------------------------------------------------------------------------------
      REAL(8) :: WI(8)        ! [moles/m^3]
      REAL(8) :: RHI          ! [0.0-1.0]
      REAL(8) :: TEMPI        ! [K]
      REAL(8) :: CNTRL(2)     ! [1] control variables

      !------------------------------------------------------------------------------------------------------
      ! Output from ISOROPIA.
      !------------------------------------------------------------------------------------------------------
      REAL(8) :: WT(8)        ! [moles/m^3]
      REAL(8) :: GAS(3)       ! [moles/m^3]
      REAL(8) :: AERLIQ(15)   ! [moles/m^3]
      REAL(8) :: AERSLD(19)   ! [moles/m^3]
      REAL(8) :: OTHER(6)     ! 
      CHARACTER(LEN=15) :: SCASI = '               '
   
      !------------------------------------------------------------------------------------------------------
      ! Parameters. Double-precision molecular weights [g/mol] and their reciprocals.
      !------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: MW_ANH4   = 18.03850D+00  ! [g/mol]
      REAL(8), PARAMETER :: MW_GNH3   = MW_ANH4       ! [g/mol] NH3  is passed as equivalent conc. of NH4+
      REAL(8), PARAMETER :: MW_ANO3   = 62.00494D+00  ! [g/mol]
      REAL(8), PARAMETER :: MW_GHNO3  = MW_ANO3       ! [g/mol] HNO3 is passed as equivalent conc. of NO3-
      REAL(8), PARAMETER :: MW_ASO4   = 96.0636D+00   ! [g/mol]
      REAL(8), PARAMETER :: MW_NA     = 22.989768D+00 ! [g/mol]
      REAL(8), PARAMETER :: MW_CL     = 35.4527D+00   ! [g/mol]
      REAL(8), PARAMETER :: MW_NACL   = 58.442468D+00 ! [g/mol]
      REAL(8), PARAMETER :: MW_H2O    = 18.01528D+00  ! [g/mol]
      REAL(4), PARAMETER :: MW_K      = 39.0983       ! [g/mol]
      REAL(4), PARAMETER :: MW_CA     = 40.078        ! [g/mol]
      REAL(4), PARAMETER :: MW_MG     = 24.3050       ! [g/mol]  
      REAL(8), PARAMETER :: RMW_NA    = 1.0D-06 / MW_NA    ! [mol/g]
      REAL(8), PARAMETER :: RMW_ASO4  = 1.0D-06 / MW_ASO4  ! [mol/g]
      REAL(8), PARAMETER :: RMW_ANH4  = 1.0D-06 / MW_ANH4  ! [mol/g]
      REAL(8), PARAMETER :: RMW_GNH3  = 1.0D-06 / MW_GNH3  ! [mol/g]
      REAL(8), PARAMETER :: RMW_ANO3  = 1.0D-06 / MW_ANO3  ! [mol/g]
      REAL(8), PARAMETER :: RMW_GHNO3 = 1.0D-06 / MW_GHNO3 ! [mol/g]
      REAL(8), PARAMETER :: RMW_CL    = 1.0D-06 / MW_CL    ! [mol/g]
      REAL(8), PARAMETER :: CMW_NA    = 1.0D+06 * MW_NA    ! [ug/mol]
      REAL(8), PARAMETER :: CMW_ASO4  = 1.0D+06 * MW_ASO4  ! [ug/mol]
      REAL(8), PARAMETER :: CMW_ANH4  = 1.0D+06 * MW_ANH4  ! [ug/mol]
      REAL(8), PARAMETER :: CMW_GNH3  = 1.0D+06 * MW_GNH3  ! [ug/mol]
      REAL(8), PARAMETER :: CMW_ANO3  = 1.0D+06 * MW_ANO3  ! [ug/mol]
      REAL(8), PARAMETER :: CMW_GHNO3 = 1.0D+06 * MW_GHNO3 ! [ug/mol]
      REAL(8), PARAMETER :: CMW_CL    = 1.0D+06 * MW_CL    ! [ug/mol]
      REAL(8), PARAMETER :: CMW_H2O   = 1.0D+06 * MW_H2O   ! [g/mol]

      !------------------------------------------------------------------------------------------------------
      ! Fraction of sea salt (NaCl) mass that is Na, and is Cl.
      !------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: RAT_NA = MW_NA / ( MW_NA + MW_CL )  ! [1] 
      REAL(8), PARAMETER :: RAT_CL = MW_CL / ( MW_NA + MW_CL )  ! [1] 
      !------------------------------------------------------------------------------------------------------
      ! Fraction of dust mass that is K, Mg, Cl-, and Ca
      !------------------------------------------------------------------------------------------------------
      REAL(4), PARAMETER :: MASS_FRAC_K  = 0.0028  ! From Ghan et al. (2001).
      REAL(4), PARAMETER :: MASS_FRAC_CA = 0.024   !   JGR, Vol. 106, p. 5295-5316.
      REAL(4), PARAMETER :: MASS_FRAC_MG = 0.0038  !   on p. 5296
      REAL(4), PARAMETER :: MASS_FRAC_NA = 0.014   !   "water sol. mass frac. in soil dust"
      REAL(4), PARAMETER :: FRAC_DUST  = 0.d0      ! [1] fraction of dust conc. passed to Thermodynamics         
      REAL(4), PARAMETER :: FRAC_SALT  = 0.d0      ! [1] fraction of salt conc. passed to Thermodynamics
      REAL(4), PARAMETER :: CONV_KION  = FRAC_DUST * MASS_FRAC_K  / MW_K  ! [mol/g]
      REAL(4), PARAMETER :: CONV_CAION = FRAC_DUST * MASS_FRAC_CA / MW_CA ! [mol/g]
      REAL(4), PARAMETER :: CONV_MGION = FRAC_DUST * MASS_FRAC_MG / MW_MG ! [mol/g]
      REAL(4), PARAMETER :: CONV_NAION = FRAC_DUST * MASS_FRAC_NA / MW_NA ! [mol/g]

      !------------------------------------------------------------------------------------------------------
      ! Other parameters.
      !------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DH2O   = 1.000D+00   ! density of water [g/cm^3]
      REAL(8), PARAMETER :: DNACL  = 2.165D+00   ! density of NaCl  [g/cm^3]
      REAL(8), PARAMETER :: CSS    = 1.08D+00    ! for sea salt ...
      REAL(8), PARAMETER :: BSS    = 1.2D+00     ! for sea salt ...              
      REAL(8), PARAMETER :: SSH2OA = (CSS*CSS*CSS*BSS-1.0D+00)*DH2O/DNACL
      REAL(8), PARAMETER :: SSH2OB = (CSS*CSS*CSS            )*DH2O/DNACL
      REAL(8), PARAMETER :: RHMAX  = 0.995D+00   ! [0-1]
      REAL(8), PARAMETER :: RHMIN  = 0.010D+00   ! [0-1]   
      REAL(8)            :: H                    ! local RH, with RHMIN < H < RHMAX

      !------------------------------------------------------------------------------------------------------
      ! Call for the bulk non-sea salt inorganic aerosol.
      !------------------------------------------------------------------------------------------------------

      H = MAX( MIN( RH, RHMAX ), RHMIN )
      WI(1) = RAT_NA*SEAS*RMW_NA*FRAC_SALT             ! sodium from [ug/m^3] to [mol/m^3]
      WI(2) =        ASO4*RMW_ASO4                     ! sulfate from [ug/m^3] to [mol/m^3]
      WI(3) =        ANH4*RMW_ANH4 +  GNH3*RMW_GNH3    ! ammonium from [ug/m^3] to [mol/m^3]
      WI(4) =        ANO3*RMW_ANO3 + GHNO3*RMW_GHNO3   ! nitrate from [ug/m^3] to [mol/m^3]
      WI(5) = RAT_CL*SEAS*RMW_CL*FRAC_SALT             ! chloride from [ug/m^3] to [mol/m^3]
      WI(6) = TOT_DUST*CONV_CAION*1.0D-06              ! calcium
      WI(7) = TOT_DUST*CONV_KION *1.0D-06              ! potassium
      WI(8) = TOT_DUST*CONV_MGION*1.0D-06              ! magnesium

      CNTRL(1) = 0.0D+00  ! Forward problem: WI contains the gas+aerosol concentrations
      CNTRL(2) = 0.0D+00  ! 0 (solid & liquid phases), 1 (liquid only, metastable)

      WT(:)     = 0.0D+00
      GAS(:)    = 0.0D+00
      AERLIQ(:) = 0.0D+00
      AERSLD(:) = 0.0D+00
      OTHER(:)  = 0.0D+00
 

      CALL ISOROPIA ( WI, H, TK, CNTRL, WT, GAS, AERLIQ, AERSLD, SCASI, OTHER )


      GNH3  = MAX( GAS(1)*CMW_GNH3,  0.0D+00 )    ! from [mol/m^3] to [ug/m^3]
      GHNO3 = MAX( GAS(2)*CMW_GHNO3, 0.0D+00 )    ! from [mol/m^3] to [ug/m^3]
      ASO4  = WT(2)*CMW_ASO4                      ! from [mol/m^3] to [ug/m^3]
      ANH4  = WT(3)*CMW_ANH4 - GNH3               ! from [mol/m^3] to [ug/m^3]
      ANO3  = WT(4)*CMW_ANO3 - GHNO3              ! from [mol/m^3] to [ug/m^3]
      AH2O  = AERLIQ(8)*CMW_H2O                   ! from [mol/m^3] to [ug/m^3]
      ANH4  = MAX( ANH4, 0.0D+00 )                ! [ug/m^3]
      ANO3  = MAX( ANO3, 0.0D+00 )                ! [ug/m^3]

      RHD   = 0.80D+00                            ! RHD = 0.80 for ammonium sulfate (Ghan et al., 2001).
      RHC   = 0.35D+00                            ! RHC = 0.35 for ammonium sulfate (Ghan et al., 2001).

 
      !-------------------------------------------------------------------------
      ! Get the sea salt-associated water (only).
      !
      ! A simple parameterization provided by E. Lewis is used.
      !-------------------------------------------------------------------------

      IF ( H .GT. 0.45D+00 ) THEN     ! ... then we are above the crystallization RH of NaCl
        SSH2O = SEAS * ( SSH2OA + SSH2OB / ( 1.0D+00 - H ) )
      ELSE
        SSH2O = 0.0D+00
      ENDIF

      RETURN
      END SUBROUTINE AERO_THERMO


