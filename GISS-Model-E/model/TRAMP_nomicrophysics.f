      SUBROUTINE AERO_NOMICROPHYSICS(AERO,GAS,EMIS_MASS,TSTEP,TK,RH,PRES,AQSO4RATE)
!-------------------------------------------------------------------------------------------------------------
!     DLW, 092106: Routine for the no-microphysics option.
!
!     When this option is chosen, mass emissions are added in, all gas-phase H2SO4 and sulfate from
!     in-cloud oxidation is added to the sulfate accumulation mode, and the routine for gas-particle
!     partitioning of NH3/NH4+, HNO3/NO3-, and aerosol H2O is called. All sulfate emissions are here 
!     put into the accumulation mode and none into the Aitken mode.
!-------------------------------------------------------------------------------------------------------------
      USE AERO_CONFIG
      USE AERO_SETUP 
      USE AERO_SUBS
      IMPLICIT NONE

      ! Arguments.
 
      REAL(8), INTENT(INOUT) :: AERO(NAEROBOX)        ! aerosol conc. [ug/m^3] or [#/m^3]
      REAL(8), INTENT(INOUT) :: GAS(NGASES)           ! gas-phase conc. [ug/m^3]
      REAL(8), INTENT(IN)    :: EMIS_MASS(NEMIS_SPCS) ! mass emission rates [ug/m^3/s]
      REAL(8), INTENT(IN)    :: TSTEP                 ! model physics time step [s]
      REAL(8), INTENT(IN)    :: TK                    ! absolute temperature [K]
      REAL(8), INTENT(IN)    :: RH                    ! relative humidity [0-1]
      REAL(8), INTENT(IN)    :: PRES                  ! ambient pressure [Pa]  
      REAL(8), INTENT(IN)    :: AQSO4RATE             ! in-cloud SO4 production rate [ug/m^3/s]

      ! Local variables.

      INTEGER :: I,J,Q                     ! indices
      REAL(8) :: PIQ(NWEIGHTS,NMASS_SPCS)  ! production terms for mass conc. [ug/m^3/s]
      REAL(8) :: TOT_SULF                  ! total sulfate conc.  [ug/m^3]
      REAL(8) :: TOT_DUST                  ! total dust conc.     [ug/m^3]
      REAL(8) :: TOT_SEAS                  ! total sea salt conc. [ug/m^3]
      REAL(8) :: SSH2O                     ! total sea salt H2O   [ug/m^3]
      REAL(8) :: AERO_WATER_ACTUAL         ! actual aerosol H2O conc. [ug/m^3]
      REAL(8) :: AERO_WATER_WET            ! wet    aerosol H2O conc. [ug/m^3]
      REAL(8) :: RHD                       ! deliquescence   RH [0-1]
      REAL(8) :: RHC                       ! crystallization RH [0-1]
      REAL(8) :: SPCMASS1(NMASS_SPCS+2)    ! initial total mass conc. of each model species [ug/m^3]
      REAL(8) :: SPCMASS2(NMASS_SPCS+2)    ! final   total mass conc. of each model species [ug/m^3]
      INTEGER, SAVE :: INDEX_ACC           ! index of sulfate accumulation mode
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      IF( FIRSTIME ) THEN
        FIRSTIME = .FALSE.
        !-----------------------------------------------------------------------------------------------------
        ! Get the mode number of the sulfate accumulation mode.
        !-----------------------------------------------------------------------------------------------------
        DO I=1, NMODES
          IF( MODE_NAME(I) .EQ. 'ACC' ) THEN
            INDEX_ACC = I
            GOTO 100
          ENDIF
        ENDDO
  100   CONTINUE
      ENDIF

      !-------------------------------------------------------------------------------------------------------
      ! Get the initial total mass concentration for each model species [ug/m^3].
      !-------------------------------------------------------------------------------------------------------
      IF ( MASS_ADJ ) CALL SPCMASSES(AERO,GAS,SPCMASS1)

      !-------------------------------------------------------------------------------------------------------
      ! Add the particle emissions.
      !-------------------------------------------------------------------------------------------------------
      ! The production terms PIQ(NWEIGHTS,NMASS_SPCS) are in [ug/m^3/s].
      !   The first  index runs over all modes or quadrature points.
      !   The second index runs over the NMASS_SPCS (=5) principal mass
      !     species: sulfate, BC, OC, dust, sea salt. 
      !
      ! Get the emissions production terms (Pemis_i,q in the manuscript) in [ug/m^3/s] and put them 
      !   in the total production rate array PIQ (the P_i,q in the manuscript).
      !
      ! EMIS_MODE_MAP and EMIS_SPCS_MAP have elements corresponding to 
      !   the aerosol types (in this order): AKK(=1), ACC(=2), BCC(=8), OCC(=7),
      !               DD1(=3), SSA(=5), SSC(=6), BOC(BC=8), BOC(OC=9), DD2(=10).
      ! EMIS_MODE_MAP(J) is mode number receiving the emissions held 
      !                  in EMIS_MASS(J).
      ! EMIS_SPCS_MAP(J) is the chemical species number (1-5) of the species
      !                  held in EMIS_MASS(J).
      ! Currently, EMIS_SPCS_MAP = (/1,1,2,3,4,5,5,2,3,4/)
      !-------------------------------------------------------------------------------------------------------
      PIQ(:,:) = 0.0D+00                                       ! Zero all production terms for this time step.
      DO J=1, NEMIS_SPCS                                       ! Loop over the emitted species. 
        PIQ( EMIS_MODE_MAP(J), EMIS_SPCS_MAP(J) ) = 
     &  PIQ( EMIS_MODE_MAP(J), EMIS_SPCS_MAP(J) ) + EMIS_MASS(J)
      ENDDO

      !-------------------------------------------------------------------------------------------------------
      ! Add all SO4(aq) and H2SO4(g) to mode ACC. Zero the H2SO4(g).
      !-------------------------------------------------------------------------------------------------------
      PIQ(INDEX_ACC,PROD_INDEX_SULF) = PIQ(INDEX_ACC,PROD_INDEX_SULF) + AQSO4RATE + GAS(GAS_H2SO4)/TSTEP
      GAS(GAS_H2SO4) = 0.0D+00

      !-------------------------------------------------------------------------------------------------------
      ! Update the sulfate, BC, OC, dust, and sea salt concentrations.
      !-------------------------------------------------------------------------------------------------------
      ! MASS_MAP(I,Q) is the location in AERO(:) of the Qth mass in mode I.
      ! Mode I has NM(I) mass species defined for it, and NM(I) varies between
      !   1 and NMASS_SPCS (=5).
      !
      ! The second index of PROD_INDEX has NMASS_SPCS (=5) values:
      !   1=sulfate, 2=BC, 3=OC, 4=dust, 5=sea salt.
      !   PROD_INDEX(I,Q) is the location in array PIQ(I,Q) of chemical species
      !   CHEM_SPC_NAME(Q) for mode (quadrature weight) I.
      !-------------------------------------------------------------------------------------------------------
      DO I=1, NWEIGHTS
        DO Q=1, NM(I)
          AERO(MASS_MAP(I,Q)) = AERO(MASS_MAP(I,Q)) + PIQ(I,PROD_INDEX(I,Q)) * TSTEP
        ENDDO
      ENDDO

      IF( NO_MICROPHYSICS_W_THERMO ) THEN
  
        !-------------------------------------------------------------------------------------------------------
        ! Get the total sulfate, total mineral dust, and total sea salt (NaCl) mass concentrations summed 
        ! over all modes [ug/m^3] for use in the subsequent gas-particle partitioning calculation.
        !-------------------------------------------------------------------------------------------------------
        TOT_SULF = SUM( AERO(SULF_MAP(:)) )
        TOT_DUST = SUM( AERO(DUST_MAP(:)) )
        TOT_SEAS = SUM( AERO(SEAS_MAP(:)) )

        !-------------------------------------------------------------------------------------------------------
        ! Do the gas-particle mass transfer after adding mass emissions, SO4(aq), and H2SO4(g).
        !-------------------------------------------------------------------------------------------------------
        ! Call the aerosol thermodynamic module to determine the bulk gas-particle partitioning of the 
        ! inorganic species and the water content associated with those species.
        ! Also determine the water content associated with the NaCl of sea salt.
        !
        ! Note that AERO(MASS_H2O) includes the sea-salt associated water when it is passed to this routine
        ! and passed to AERO_THERMO. Upon return from AERO_THERMO, AERO(MASS_H2O) contains only the 
        ! non-sea salt-associated water. The sea salt-associated water is in SSH2O.
        !-------------------------------------------------------------------------------------------------------
        AERO_WATER_ACTUAL = AERO(MASS_H2O)                                ! actual tracked aerosol water conc.
        CALL AERO_THERMO(TOT_SULF,AERO(MASS_NO3),AERO(MASS_NH4),
     &                   AERO(MASS_H2O),GAS(GAS_NH3),GAS(GAS_HNO3),
     &                   TOT_DUST,TOT_SEAS,SSH2O,TK,RH,PRES,RHD,RHC)
        AERO_WATER_WET = AERO(MASS_H2O) + SSH2O                           ! total metastable aerosol water conc.
  
        !-------------------------------------------------------------------------------------------------------
        ! Adjust the aerosol water concentration for hysteresis.
        !-------------------------------------------------------------------------------------------------------
        ! If the aerosol water concentration is less than half its metastable (wet)
        ! concentration, treat the aerosol as dry. Otherwise, the values
        ! in AERO(MASS_H2O) and SSH2O remain at their wet values. Since the water
        ! associated with sea salt is not tracked separately, this treatment of
        ! hysteresis is based on the total aerosol water including that of sea
        ! salt, although the RHD governing this is that obtained (or set) for the
        ! non-sea salt inorganics.
        !
        ! This is done only for RH between the crystallization and deliquescence
        ! RHs of the non-sea salt inorganics. 
        !
        ! On 8-7-06: RHC = 80%   as for ammonium sulfate (Ghan et al. 2001)
        !            RHD = 35%   as for ammonium sulfate (Ghan et al. 2001)
        !-------------------------------------------------------------------------------------------------------
        IF ( RH .GT. RHC  .AND.  RH .LT. RHD ) THEN
          IF ( AERO_WATER_WET .GT. 0.0D+00 ) THEN
            IF ( AERO_WATER_ACTUAL/AERO_WATER_WET .LT. 0.5D+00 ) THEN  ! dry aerosol
              AERO(MASS_H2O) = 0.0D+00                                 ! Zero the non-sea salt-associated water.
              SSH2O = 0.0D+00                                          ! Zero the     sea salt-associated water.
            ELSE                                                       ! wet (metastable) aerosol
                                                                       ! Leave AERO(MASS_H2O) and SSH2O 
            ENDIF                                                      !   at their metastable (wet) values. 
          ENDIF
        ELSEIF ( RH .LE. RHC ) THEN                                    ! Insure that the aerosol is dry below the RHC.
          AERO(MASS_H2O) = 0.0D+00
          SSH2O = 0.0D+00                                              ! Since the RHC of NaCl is set to 45%, and RHC
        ENDIF                                                          ! is currently set to 35% (ammonium sulfate),
                                                                       ! SSH2O should already be zero here. 
        !-------------------------------------------------------------------------------------------------------
        ! The total aerosol water must exit the routine in AERO(MASS_H2O).
        !-------------------------------------------------------------------------------------------------------
        AERO(MASS_H2O) = AERO(MASS_H2O) + SSH2O                        ! Total aerosol water now.

      ENDIF

      !-------------------------------------------------------------------------------------------------------
      ! Get final total mass concentration for each model species [ug/m^3].
      ! Rescale final aerosol and gas-phase species to enforce mass conservation to machine precision.
      !-------------------------------------------------------------------------------------------------------
      IF ( MASS_ADJ ) THEN
        CALL SPCMASSES(AERO,GAS,SPCMASS2)
        CALL MASSADJ(AERO,GAS,SPCMASS1,SPCMASS2,EMIS_MASS,AQSO4RATE,TSTEP)
      ENDIF

      !-------------------------------------------------------------------------------------------------------
      ! Limit low mass or number concentrations.
      !-------------------------------------------------------------------------------------------------------
      AERO(:) = MAX( AERO(:), MINCONC )

      RETURN
      END SUBROUTINE AERO_NOMICROPHYSICS

