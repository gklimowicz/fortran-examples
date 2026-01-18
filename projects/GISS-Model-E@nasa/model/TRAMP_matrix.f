      SUBROUTINE MATRIX(AERO,GAS,EMIS_MASS,TSTEP,TK,RH,PRES,AQSO4RATE,WUPDRAFT,DIAG,PM)
!-----------------------------------------------------------------------------------------------------------------------
!
!@sum     This is the top-level routine of the MATRIX aerosol microphysical module.
!@auth    Susanne Bauer/Doug Wright
!
!     Here are some items to be familiar with in setting up the call to MATRIX.
!
!     1. The GAS array: 
!
!        GAS(1) contains the current H2SO4(g) concentration expressed as
!        ug SO4=/m^3, including the H2SO4(g) produced by gas-phase chemistry
!        during the current time step.
!     
!        GAS(2) contains the current HNO3(g) concentration expressed as
!        ug NO3-/m^3, including the HNO3(g) produced by gas-phase chemistry
!        during the current time step.
!     
!        GAS(3) contains the current NH3(g) concentration expressed as
!        ug NH4+/m^3. Any sources of NH3(g) during the time step should  
!        probably be added in before the call to MATRIX.
!
!     2. H2SO4 is always represented as SO4= (MW=96g/mol) rather than 
!        H2SO4 (MW=98g/mol) (except in the calculation of the mean thermal
!        velocity of an H2SO4 molecule). Likewise, HNO3 is always 
!        represented as NO3- and NH3 as NH4+.
!
!     3. In this model, sea salt contains NaCl and potentially non-sea salt 
!        sulfate, nitrate, ammonium, and water, depending upon conditions.
!        The model sea salt species MASS_SSA_SEAS (accumulation mode),
!        MASS_SSC_SEAS (coarse mode), and MASS_MXX_SEAS (sea salt in 
!        a mixed mode) do not include the non-sea salt sulfate, nitrate,
!        ammonium, and water, are treated as pure NaCl, and the density
!        and MW of NaCl are used with these species.
!
!     4. EMIS_MASS(J) contains the mass emissions rates for species J. 
!        The species and units are:
!
!          J               SPECIES                      UNITS
!        -----   -----------------------------    ------------------
!          1       Aitken-mode sulfate              [ug SO4 /m^3/s]
!          2       accumulation-mode sulfate        [ug SO4 /m^3/s]
!          3       black carbon                     [ug BC  /m^3/s]
!          4       organic carbon                   [ug OC  /m^3/s]
!          5       mineral dust - size 1            [ug DUST/m^3/s]
!          6       accumulation-mode sea salt       [ug NaCl/m^3/s]
!          7       coarse-mode sea salt             [ug NaCl/m^3/s]
!          8       BC in the mixed BC-OC mode       [ug BC  /m^3/s]
!          9       OC in the mixed BC-OC mode       [ug OC  /m^3/s]
!         10       mineral dust - size 2            [ug DUST/m^3/s]
!
!-----------------------------------------------------------------------------------------------------------------------
      USE AERO_CONFIG   
      USE AERO_SETUP     ! Module AERO_SETUP uses module AERO_PARAM, so the statement USE AERO_PARAM is not needed here.  
      USE AERO_SUBS  
      USE AERO_COAG,   ONLY: SETUP_KIJ_DIAMETERS, SETUP_KIJ_TABLES, GET_KBARNIJ
      USE AERO_NPF,    ONLY: DNU, NPFRATE, SETUP_NPFMASS, STEADY_STATE_H2SO4   
!      USE AERO_DIAM,   ONLY: DIAM !, DIAM_HISTOGRAM
      USE AERO_ACTV,   ONLY: GETACTFRAC 
      USE AMP_AEROSOL, ONLY: NACTV , DIAM_DRY, DIAM  
      USE AERO_DEPV,   ONLY: GET_AERO_DEPV, VDDEP_AERO  
      IMPLICIT NONE

      ! Arguments.
 
      REAL(8), INTENT(INOUT) :: AERO(NAEROBOX)             ! aerosol conc. [ug/m^3] or [#/m^3]
      REAL(8), INTENT(INOUT) :: GAS(NGASES)                ! gas-phase conc. [ug/m^3]
      REAL(8), INTENT(INOUT) :: EMIS_MASS(NEMIS_SPCS)      ! mass emission rates [ug/m^3/s]
      REAL(8), INTENT(IN)    :: TSTEP                      ! model physics time step [s]
      REAL(8), INTENT(IN)    :: TK                         ! absolute temperature [K]
      REAL(8), INTENT(IN)    :: RH                         ! relative humidity [0-1]
      REAL(8), INTENT(IN)    :: PRES                       ! ambient pressure [Pa]  
      REAL(8), INTENT(IN)    :: AQSO4RATE                  ! in-cloud SO4 production rate [ug/m^3/s]
      REAL(8), INTENT(IN)    :: WUPDRAFT                   ! cloud updraft velocity [m/s]
      REAL(8), INTENT(INOUT) :: DIAG(NDIAG_AERO,NAEROBOX)  ! budget or tendency diagnostics [ug/m^3/s] or [#/m^3/s]      
      REAL(8), INTENT(INOUT) :: PM(3)                      ! Particulate Matter, PM1, PM2,5, PM10 at ambient humidity
      ! Local variables.

      INTEGER :: I,J,K,L,Q,QQ              ! indices
      INTEGER :: INDEX_DP!, INDEX_DP_DRY    ! index for condensation factor lookup table
      INTEGER :: IBRANCH                   ! scratch debugging variable [1]
      REAL(8) :: BI(NWEIGHTS)              ! number conc. coefficients [1/s]
      REAL(8) :: CI(NWEIGHTS)              ! number conc. coefficients [#/m^3/s]
      REAL(8) :: RI(NWEIGHTS)              ! number conc. coefficients [#/m^3/s]
      REAL(8) :: NI(NWEIGHTS)              ! number concentrations [#/m^3]
      REAL(8) :: MJQ(NWEIGHTS,NMASS_SPCS)  ! MJQ(J,Q) is the avg. mass/particle of species Q (=1-5) for mode J
      REAL(8) :: PIQ(NWEIGHTS,NMASS_SPCS)  ! production terms for mass conc. [ug/m^3/s]
      REAL(8) :: FI(NWEIGHTS)              ! loss coefficients for mass conc. [1/s]
      REAL(8) :: TOT_SULF                  ! total sulfate conc.  [ug/m^3]
      REAL(8) :: TOT_DUST                  ! total dust conc.     [ug/m^3]
      REAL(8) :: TOT_SEAS                  ! total sea salt conc. [ug/m^3]
      REAL(8) :: SSH2O                     ! total sea salt H2O   [ug/m^3]
      REAL(8) :: SSH2O_PER_SSMASS          ! total sea salt H2O / total sea-salt dry mass
      REAL(8) :: FTMP, VOLTMP, VOLTMP_DRY  ! scratch variables for computing mode mean diameters 
      REAL(8) :: TOT_MASS    (NMODES)      ! total ambient mass conc. for each mode [ug/m^3]
      REAL(8) :: TOT_MASS_DRY(NMODES)      ! total dry     mass conc. for each mode [ug/m^3]
      REAL(8) :: H2O_MASS(NMODES)          ! water mass mass conc. for each mode [ug/m^3]
      REAL(8) :: MASS_COMP(NMODES,8)       ! mass conc. of each component for each mode [ug/m^3]
      REAL(8) :: DENS_MODE    (NMODES)     ! average mode density calculated from component concentrations [g/cm^3]
      REAL(8) :: DENS_MODE_DRY(NMODES)     ! average mode density calculated from component concentrations [g/cm^3]
      REAL(8) :: OPTOT_NO3NH4H2O_TO_SULF   ! 1 + ( total NO3+NH4+H2O / total SO4 )
      REAL(8) :: AERO_WATER_ACTUAL         ! actual aerosol H2O conc. [ug/m^3]
      REAL(8) :: AERO_WATER_WET            ! wet    aerosol H2O conc. [ug/m^3]
      REAL(8) :: RHD                       ! deliquescence   RH [0-1]
      REAL(8) :: RHC                       ! crystallization RH [0-1]
      REAL(8) :: DGN(NWEIGHTS)             ! ambient geometric mean diameter of the number distribution for each mode [m]
      REAL(8) :: DGN_DRY(NWEIGHTS)         ! geometric mean dry diameter of the number distribution for each mode [m]
      REAL(8) :: DP(NWEIGHTS)              ! ambient diameter of average mass of the number distribution for each mode [m]
      REAL(8) :: DP_DRY(NWEIGHTS)          ! dry diameter of average mass of the number distribution for each mode [m]
      REAL(8) :: P_EMIS_NUMB(NMODES)       ! number emission rates [#/m^3/s]       
      REAL(8) :: SPCMASS1(NMASS_SPCS+2)    ! initial total mass conc. of each model species [ug/m^3]
      REAL(8) :: SPCMASS2(NMASS_SPCS+2)    ! final   total mass conc. of each model species [ug/m^3]
      REAL(8) :: FSEASSULF                 ! fraction of sulfate assigned to sea salt in mode SSC (coarse mode)
      REAL(8) :: FACTOR                    ! scratch variable in mass conc. equation solver
      REAL(8) :: Y0,Y,A,B,C,DELTA,R1,R2    ! scratch variables in the number concentration analytic solution
      REAL(8) :: GAMMA,GEXPDT,EXPDT        ! scratch variables in the number concentration analytic solution
      REAL(8), PARAMETER :: PIQ_THRESH = 1.0D-08  ! [1] threshold in number/mass conc. solver

      ! For the condensational sink, condensational growth, and cou. 

      REAL(8) :: KC                         ! total condensational sink with arbitrary mass accommodation coef. [1/s]
      REAL(8) :: KC_AEQ1                    ! total condensational sink with unity     mass accommodation coef. [1/s]
      REAL(8) :: PQ_GROWTH                  ! avg. cond. rate of H2SO4 (as SO4) over the time step [ugSO4/m^3/s]
      REAL(8) :: TOT_H2SO4_LOSS             ! net loss of H2SO4 (as SO4) [ugSO4/m^3]
      !----------------------------------------------------------------------------------------------------------------
      ! This minimum value for the condensational sink, 1.0D-08 [1/s], is calculated as 
      ! 
      !   KC = 2 * PI * D * Beta * Dpbar * N = 2 * 3.14 * (1.0D-05 m^2/s) * (1) * (1.0D-07 m) * (1.0D+04 #/m^3)
      !      = 6.28D-08 1/s --> 1.0D-08 1/s  
      !----------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: KCMIN = 1.0D-08 ! [1/s] minimum condensational sink - see notes of 10-18-06

      ! Mode-average coagulation coefficients. 

      REAL(8), SAVE :: KBAR0IJ(NWEIGHTS,NWEIGHTS)     ! mode-average coagulation coefficients [m^3/s] for number
      REAL(8), SAVE :: KBAR3IJ(NWEIGHTS,NWEIGHTS)     ! mode-average coagulation coefficients [m^3/s] for mass 

      ! These variables are used in intermodal transfer.

      INTEGER, PARAMETER :: IMTR_EXP = 4              ! exponent for AKK --> ACC intermodal transfer
      REAL(8)            :: DPAKK                     ! diameter of average mass for mode AKK [m]
      REAL(8)            :: DPACC                     ! diameter of average mass for mode ACC [m] 
      REAL(8)            :: FNUM                      ! fraction of AKK number transferred over the time step [1]
      REAL(8)            :: F3                        ! fraction of AKK mass   transferred over the time step [1]
      REAL(8)            :: XNUM, X3                  ! error function complement arguments [1]
      REAL(8)            :: DGN_AKK_IMTR              ! geo. mean diam. of the AKK mode number distribution [m]
      REAL(8)            :: DGN_ACC_IMTR              ! geo. mean diam. of the ACC mode number distribution [m]
      REAL(8)            :: DEL_NUMB                  ! number conc. transferred from AKK to ACC [ #/m^3]
      REAL(8)            :: DEL_MASS                  ! mass   conc. transferred from AKK to ACC [ug/m^3]
      REAL(8), SAVE      :: DPCUT_IMTR      = 0.0D+00 ! fixed diameter for intermodal transfer [um]
      REAL(8), SAVE      :: XNUM_FACTOR     = 0.0D+00 ! factor in XNUM expression [1]
      REAL(8), SAVE      :: X3_TERM         = 0.0D+00 ! term   in X3   expression [1]
      REAL(8), SAVE      :: LNSG_AKK        = 0.0D+00 ! ln(SG_AKK) [1]
      REAL(8), SAVE      :: LNSG_ACC        = 0.0D+00 ! ln(SG_ACC) [1]
      REAL(8), PARAMETER :: DPAKK0          = DNU     ! min. diameter of average mass for mode AKK [m]
      REAL(8), PARAMETER :: FNUM_MAX        = 0.5D+00 ! max. value of FNUM [1]
      REAL(8), PARAMETER :: AKK_MINNUM_IMTR = 1.0D+06 ! min. AKK number conc. to enable IMTR [#/m^3]

      ! These variables and parameters are used in nucleation and new particle formation.

      integer :: klq
      integer :: ikl
      INTEGER :: IH2SO4_PATH               ! index for branch in the [H2SO4] calc. for nucl and GR [1]
      REAL(8) :: XH2SO4_SS                 ! steady-state H2SO4 (as SO4) conc. [ugSO4/m^3]
      REAL(8) :: XH2SO4_SS_WNPF            ! steady-state H2SO4 (as SO4) conc. including NPF [ugSO4/m^3]
      REAL(8) :: XH2SO4_NUCL               ! H2SO4 (as SO4) conc. used in nucleation and GR calculation [ugSO4/m^3]
      REAL(8) :: XH2SO4_INIT               ! H2SO4 at the top of the time step [ugSO4/m^3]
      REAL(8) :: SO4RATE                   ! H2SO4(g) production rate in/as [ugSO4/m^3/s]
      REAL(8) :: XNH3                      ! ammonia concentration [ppmV]
      REAL(8) :: DNDT                      ! number production rate [1/m^3/s]     from npf
      REAL(8) :: DMDT_SO4                  ! SO4 mass product. rate [ugSO4/m^3/s] from npf
      REAL(8), PARAMETER :: XH2SO4_NUCL_MIN_NCM3 = 1.00D+03 ! min. [H2SO4] to enter nucleation calculations [#/cm^3]  
      REAL(8), PARAMETER :: XH2SO4_NUCL_MIN = XH2SO4_NUCL_MIN_NCM3 * MW_SO4 * 1.0D+12 / AVO  ! convert to [ugSO4/m^3] 
      REAL(8), PARAMETER :: XNTAU = 2.0D+00                 ! number of time constants in the current time step  
                                                            ! required to invoke the steady-state approx. [1] 
      ! These variables are used in aerosol activation.

      REAL(8) :: RSUM_ACTIV        ! used in Pcloud_i,q terms
      REAL(8) :: NSOL    (NMODES)  ! number concentration of soluble particles for each mode [#/m^3] 
      REAL(8) :: AC      (NMODES)  ! minimum dry radius to activate for each mode [um].
      REAL(8) :: FRACACTN(NMODES)  ! activated fraction of number concentration for each mode [1]
      REAL(8) :: FRACACTM(NMODES)  ! activated fraction of mass   concentration for each mode [1]
      REAL(8) :: NACT    (NMODES)  ! activated number concentration for each mode [#/m^3] 
      REAL(8) :: MACT    (NMODES)  ! activated mass   concentration for each mode [ug/m^3] 
      REAL(8) :: MI5     (NMODES,NMASS_SPCS)  ! mass conc. of each species for each mode [ug/m^3] 

      ! These variables are used as default values in computing dry deposition velocities.

      REAL(8), PARAMETER :: TEMP_DDEP = 288.15D+00      ! temperature [ K ]
      REAL(8), PARAMETER :: RHOA_DDEP =  1.225D+00      ! air density [ kg/m^3 ]
      REAL(8), PARAMETER :: LAMB_DDEP = 6.6332D-08      ! mean free path of air [ m ]
      REAL(8), PARAMETER :: DVIS_DDEP = 1.7894D-05      ! dynamic viscosity of air [ kg/(m s) ]
      REAL(8), PARAMETER :: WSTR_DDEP =    1.0D+00      ! convective velocity scale [ m/s ]
      REAL(8), PARAMETER :: USTR_DDEP =    0.5D+00      ! friction velocity [ m/s ]
      REAL(8), PARAMETER :: RAER_DDEP =    5.0D+00      ! aerodynamic resistance [ s/m ]
      LOGICAL, PARAMETER :: GET_DEP_VEL_ONLY = .FALSE.  ! For early RETURN after getting dep. velocities

      ! These variables are used in computing the DIAG or DIAM_HISTOGRAM diagnostic arrays. 

      REAL(8) :: DIAGTMP1(NDIAG_AERO,NAEROBOX)          ! scratch work array [ug/m^3/s] or [#/m^3/s]
      REAL(8) :: AEROTMP1(NAEROBOX)                     ! input   aerosol concentrations [ug/m^3] or [#/m^3]
      REAL(8) :: AEROTMP2(NAEROBOX)                     ! scratch aerosol concentrations [ug/m^3] or [#/m^3]
      REAL(8) :: PIQTMP(NWEIGHTS,NMASS_SPCS)            ! scratch work array for mass production terms [ug/m^3/s]
!      REAL(8), PARAMETER :: N_MIN_DIAM_HISTOGRAM = 1.0D+04  ! min. # conc. for count in DIAM_HISTOGRAM [#/m^3] 
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      !----------------------------------------------------------------------------------------------------------------
      ! Error function statement function derived from code in CMAQ v4.4. 
      ! 
      ! The CMAQ source is Meng, Z., and J.H.Seinfeld (1994) On the source of the submicrometer droplet mode of 
      ! urban and regional aerosols. Aerosol Sci. and Technology, 20:253-265, who cite
      ! Reasearch & Education Association (REA), (1991) Handbook of Mathematical, Scientific, 
      ! and Engineering Formulas, Tables, Functions, Graphs, Transforms: REA, Piscataway, NJ. p. 493.
      !----------------------------------------------------------------------------------------------------------------
       REAL(8), PARAMETER :: ERFCONST = -4.0D+00 / PI
       REAL(8)            :: XX, ERF1, ERFC, ERF                          ! Error Function, Error Function Complement
       ERF1(XX) = SQRT(1.0D+00 - EXP( ERFCONST * XX * XX ) )
       ERFC(XX) = 1.0D+00 - ERF1(XX)

!----------------------------------------------------------------------------------------------------------------------
!     Begin execution.
!----------------------------------------------------------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
        FIRSTIME = .FALSE.
        CALL SETUP_CONFIG           
        CALL SETUP_SPECIES_MAPS
        CALL SETUP_AERO_MASS_MAP 
        IF ( .NOT. NO_MICROPHYSICS ) THEN
          CALL SETUP_COAG_TENSORS
          CALL SETUP_DP0     
          CALL SETUP_KIJ_DIAMETERS
          CALL SETUP_KIJ_TABLES
          IF( UPDATE_KIJ .EQ. 0 ) THEN
            DGN(:)     = DGN0(:)
            DGN_DRY(:) = DGN0(:)
            CALL GET_KBARNIJ(0,TK,PRES,DGN,KBAR0IJ,KBAR3IJ)
          ENDIF
        ENDIF
        CALL SETUP_EMIS  
        CALL SETUP_KCI   
        CALL SETUP_NPFMASS
        DP    (:) = DP0(:)  ! Set particle diameters to default values [m] - optionally overwritten below.
        DP_DRY(:) = DP0(:)  ! Set particle diameters to default values [m] - optionally overwritten below.
        
        ! These variables are used in intermodal transfer.
        
        IF( INTERMODAL_TRANSFER .AND. NUMB_AKK_1 .EQ. 0 ) THEN
          WRITE(*,*)'INTERMODAL_TRANSFER must be set to .FALSE. for this mechanism.'
          STOP
        ENDIF
        DPCUT_IMTR  = SQRT( DG_AKK * DG_ACC )   ! This is the formula of Easter et al. 2004. (= 0.053 um)
        LNSG_AKK    = LOG( SG_AKK )
        LNSG_ACC    = LOG( SG_ACC )
        XNUM_FACTOR = 1.0D+00 / ( SQRT( 2.0D+00 ) * LNSG_AKK )
        X3_TERM     = 3.0D+00 * LNSG_AKK / SQRT( 2.0D+00 )
        IF( INTERMODAL_TRANSFER ) THEN
          IF( WRITE_LOG ) WRITE(AUNIT1,'(/A/)') 'INTERMODAL TRANSFER (AKK->ACC) IS TURNED ON.'
          IF( WRITE_LOG ) WRITE(AUNIT1,'(A,5F12.6)')'DPCUT_IMTR, XNUM_FACTOR, X3_TERM, LNSG_AKK, LNSG_ACC = ', 
     &                                               DPCUT_IMTR, XNUM_FACTOR, X3_TERM, LNSG_AKK, LNSG_ACC
        ELSE
          IF( WRITE_LOG ) WRITE(AUNIT1,'(/A/)') 'INTERMODAL TRANSFER (AKK->ACC) IS TURNED OFF.'
        ENDIF


        ! Calculate and store dry deposition velocities for the entire X-Y grid. 

        CALL GET_AERO_DEPV(NMODES,TEMP_DDEP,RHOA_DDEP,LAMB_DDEP,DVIS_DDEP,
     &                     WSTR_DDEP,USTR_DDEP,RAER_DDEP,DGN0,LNSIG0,DENSPI)
        DO I=1, NMODES   ! Initialize entire X-Y grid. 
          VDDEP_AERO(:,:,I,1) = VDDEP_AERO(IXXX,IYYY,I,1)  ! For deposition of number concentrations; [m/s]. 
          VDDEP_AERO(:,:,I,2) = VDDEP_AERO(IXXX,IYYY,I,2)  ! For deposition of mass   concentrations; [m/s].
        ENDDO
      ENDIF                           

      !----------------------------------------------------------------------------------------------------------------
      ! The mass-only, no-microphysics option.
      ! Add in emissions, H2SO4(g), SO4(aq), and do gas-particle partitioning only.
      !---------------------------------------------------------------------------------------------------------------- 
      IF( NO_MICROPHYSICS ) THEN
         EMIS_MASS(:) = 0.0D+00 
        CALL AERO_NOMICROPHYSICS(AERO,GAS,EMIS_MASS,TSTEP,TK,RH,PRES,AQSO4RATE)
        RETURN
      ENDIF

      !----------------------------------------------------------------------------------------------------------------
      ! Zero entire diagnostics array. 
      !----------------------------------------------------------------------------------------------------------------
      DIAGTMP1(:,:) = 0.0D+00 
      AEROTMP1(:)   = AERO(:)   ! Input concentrations
      
      !----------------------------------------------------------------------------------------------------------------
      ! Get initial total mass conc. for each model species [ug/m^3].
      ! These are (optionally) used to enforce precise mass conservation.
      !----------------------------------------------------------------------------------------------------------------
      IF( MASS_ADJ ) CALL SPCMASSES(AERO,GAS,SPCMASS1)

      !----------------------------------------------------------------------------------------------------------------
      ! Load the number concentrations [#/m^3] into the local work array.
      !----------------------------------------------------------------------------------------------------------------
      NI(:) = AERO( NUMB_MAP(:) ) + TINYDENOM

      !----------------------------------------------------------------------------------------------------------------
      ! IF( WRITE_LOG ) THEN
      !   WRITE(AUNIT1,'(/5A20/)') 'I','NI [#/m^3] - A'
      !   DO I=1, NWEIGHTS
      !     WRITE(AUNIT1,90004) I, NI(I)
      !   ENDDO
      ! ENDIF
      !----------------------------------------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------------------------
      ! For the sea salt modes, use characteristic mean particle masses
      ! to derive the number concentration from the mass concentration.
      !
      ! There may be two sea salt modes, SSA and SSC, or just one, SSS.
      !
      ! SEAS_MODE_MAP(I)      is the mode number of sea salt mode I.
      ! SEAS_MODE_MASS_MAP(I) is the location in AERO of the sea salt (NaCl)
      !                         mass concentration for sea salt mode I.
      ! RECIP_SEAS_MPP(I)     is the reciprocal of the mean NaCl mass 
      !                         per sea salt particle for sea salt mode I.
      !----------------------------------------------------------------------------------------------------------------
      NI(SEAS_MODE_MAP(:)) = AERO( SEAS_MODE_MASS_MAP(:) ) * RECIP_SEAS_MPP(:)
      ! WRITE(*,*) NI(SEAS_MODE_MAP(:)) 

      !----------------------------------------------------------------------------------------------------------------
      ! IF( WRITE_LOG ) THEN
      !   WRITE(AUNIT1,'(/5A20/)') 'I','NI [#/m^3] - B'
      !   DO I=1, NWEIGHTS
      !     WRITE(AUNIT1,90004) I, NI(I)
      !   ENDDO
      ! ENDIF
      !----------------------------------------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------------------------
      ! Partition the sea salt sulfate between the SSA and SSC modes, if both
      ! of these modes are present. In this case, all sea salt sulfate is
      ! passed to this routine in mass species MASS_SSA_SULF.
      !----------------------------------------------------------------------------------------------------------------
      ! WRITE(*,*)'NUMBER_OF_SEASALT_MODES = ', NUMBER_OF_SEASALT_MODES
      IF( NUMBER_OF_SEASALT_MODES .EQ. 2 ) THEN
        FSEASSULF = AERO( MASS_SSC_SEAS ) / ( AERO( MASS_SSC_SEAS ) + AERO( MASS_SSA_SEAS ) + TINYDENOM )
        AERO( MASS_SSC_SULF ) = AERO( MASS_SSA_SULF ) * FSEASSULF
        AERO( MASS_SSA_SULF ) = AERO( MASS_SSA_SULF ) * ( 1.0D+00 - FSEASSULF )
      ENDIF

      !----------------------------------------------------------------------------------------------------------------
      ! Get the total sulfate, total mineral dust, and total sea salt (NaCl)
      ! mass concentrations summed over all quadrature points/modes. [ug/m^3]
      !----------------------------------------------------------------------------------------------------------------
      TOT_SULF = SUM( AERO(SULF_MAP(:)) ) + TINYNUMER   ! insure nonzero value
      TOT_DUST = SUM( AERO(DUST_MAP(:)) )
      TOT_SEAS = SUM( AERO(SEAS_MAP(:)) ) + TINYNUMER   ! insure nonzero value

      IF( WRITE_LOG ) WRITE(AUNIT1,'(/A,F12.4/)') 'TOT_SULF = ', TOT_SULF


      !----------------------------------------------------------------------------------------------------------------
      ! Call the aerosol thermodynamic module to determine the bulk gas-particle 
      ! partitioning of the inorganic species and the water content
      ! associated with those species. Also determine the water content
      ! associated with the NaCl of sea salt.
      !
      ! Note that AERO(MASS_H2O) includes the sea-salt associated water 
      ! when it is passed to this routine and passed to AERO_THERMO.
      ! Upon return from AERO_THERMO, AERO(MASS_H2O) contains only the 
      ! non-sea salt-associated water. The sea salt-associated water is in SSH2O.
      !----------------------------------------------------------------------------------------------------------------
      AERO_WATER_ACTUAL = AERO(MASS_H2O)       ! actual tracked aerosol water conc.
      CALL AERO_THERMO(TOT_SULF,AERO(MASS_NO3),AERO(MASS_NH4),AERO(MASS_H2O),GAS(GAS_NH3),
     &                 GAS(GAS_HNO3),TOT_DUST,TOT_SEAS,SSH2O,TK,RH,PRES,RHD,RHC)
      AERO_WATER_WET = AERO(MASS_H2O) + SSH2O  ! total metastable aerosol water conc.

      !----------------------------------------------------------------------------------------------------------------
      ! Adjust the aerosol water concentration for hysteresis.
      !
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
      ! Currently, RHC = 80%, as for ammonium sulfate (Ghan et al. 2001), and 
      !            RHD = 35%, as for ammonium sulfate (Ghan et al. 2001), both set in subr. AERO_THERMO. 
      !----------------------------------------------------------------------------------------------------------------
      IF ( RH .GT. RHC  .AND.  RH .LT. RHD ) THEN
        IF ( AERO_WATER_WET .GT. 0.0D+00 ) THEN
          IF ( AERO_WATER_ACTUAL/AERO_WATER_WET .LT. 0.5D+00 ) THEN  ! dry aerosol
            AERO(MASS_H2O) = 0.0D+00  ! Zero the non-sea salt-associated water.
            SSH2O = 0.0D+00           ! Zero the     sea salt-associated water.
          ELSE                                                       ! wet (metastable) aerosol
            ! Leave AERO(MASS_H2O) and SSH2O at their metastable (wet) values. 
          ENDIF
        ENDIF
      ELSEIF ( RH .LE. RHC ) THEN     ! Insure that the aerosol is dry below the RHC.
        AERO(MASS_H2O) = 0.0D+00
        SSH2O = 0.0D+00               ! Since the RHC of NaCl is set to 45%, and RHC
      ENDIF                           ! is currently set to 35% (ammonium sulfate),
                                      ! SSH2O should already be zero here. 


      !----------------------------------------------------------------------------------------------------------------
      ! Get the mass of species Q in mode (or quadrature point) I for the
      !   principal aerosol-phase chemical species.
      !   These species are sulfate, BC, OC, dust, and sea salt.
      !
      ! MASS_MAP(I,Q) is the location in AERO(:) of the Qth mass in mode I.
      ! Mode I has NM(I) mass species defined for it, and NM(I) varies between 1 and NMASS_SPCS (=5).
      ! The second index of PROD_INDEX has NMASS_SPCS (=5) values:
      !   1=sulfate, 2=BC, 3=OC, 4=dust, 5=sea salt.
      ! The second index of MJQ also has these same NMASS_SPCS (=5) values.
      !----------------------------------------------------------------------------------------------------------------
      IF( WRITE_LOG ) WRITE(AUNIT1,'(/A/)')'I,Q,MJQ(I,QQ),AERO(MASS_MAP(I,Q)),NI(I),MASS_MAP(I,Q)'
      DO I=1, NWEIGHTS                               ! loop over modes (quadrature points)
        IF ( NI(I) .GT. 1.0D-10 ) THEN        
          MJQ(I,:) = 1.0D-15                         ! [ug/particle] - all species in mode I
          DO Q=1, NM(I)                              ! loop over species defined for mode I
            QQ = PROD_INDEX(I,Q)                     ! QQ is in the range 1 to 5.
            MJQ(I,QQ) = AERO(MASS_MAP(I,Q)) / NI(I)  ! [ug/particle]
            MJQ(I,QQ) = MAX( MJQ(I,QQ), 1.0D-15)     ! [ug/particle] - all species in mode I
            ! WRITE(30,'(2I4,3D15.5,I4)') I,Q,MJQ(I,QQ),AERO(MASS_MAP(I,Q)),NI(I),MASS_MAP(I,Q)
            IF( WRITE_LOG ) WRITE(AUNIT1,'(2I4,3D15.5,I4)') I,Q,MJQ(I,QQ),AERO(MASS_MAP(I,Q)),NI(I),MASS_MAP(I,Q)
          ENDDO
        ELSE
          MJQ(I,:) = 1.0D-15                         ! [ug/particle] - all species in mode I
        ENDIF
      ENDDO

      !----------------------------------------------------------------------------------------------------------------
      ! Get the ambient diameter of average mass for each mode.
      !----------------------------------------------------------------------------------------------------------------

      IF ( UPDATE_DP ) THEN
        TOT_SULF = SUM( AERO(SULF_MAP(:)) ) + TINYNUMER   ! insure nonzero value
        MASS_COMP(:,:)  = 0.0D+00    ! [ug/m^3] mass of each component in each mode
        TOT_MASS    (:) = 0.0D+00    ! [ug/m^3] total mass in each mode
        TOT_MASS_DRY(:) = 0.0D+00    ! [ug/m^3] total mass in each mode
        DO I=1, NWEIGHTS
          DO Q=1, NM(I)
            QQ = PROD_INDEX(I,Q)    ! QQ is in the range 1 to 5.
            !----------------------------------------------------------------------------------------------------------
            ! Sulfate, BC, OC, dust, and sea salt concentrations.
            !----------------------------------------------------------------------------------------------------------
            MASS_COMP(I,QQ) = AERO(MASS_MAP(I,Q))   ! [ug/m^3]
          ENDDO
          !------------------------------------------------------------------------------------------------------------
          ! Nitrate, ammonium, and (non-sea salt) water concentrations.
          !------------------------------------------------------------------------------------------------------------
          FTMP = AERO(SULF_MAP(I)) / TOT_SULF
          MASS_COMP(I,6) = FTMP * AERO(1)           ! [ug/m^3] nitrate
          MASS_COMP(I,7) = FTMP * AERO(2)           ! [ug/m^3] ammonium
          MASS_COMP(I,8) = FTMP * AERO(3)           ! [ug/m^3] water
          ! WRITE(*,'(I5,6D15.5)')I,FTMP,AERO(1:3), TOT_SULF, AERO(SULF_MAP(I))
        ENDDO
        !--------------------------------------------------------------------------------------------------------------
        ! Distribute the sea-salt associated water over modes in proportion to the sea salt mass present.
        !--------------------------------------------------------------------------------------------------------------

        SSH2O_PER_SSMASS = SSH2O / TOT_SEAS                                  ! [ugH2O/ugNaCl]
        DO J=1, NMODES_SEAS    ! loop over all modes containing sea salt
          MASS_COMP(MODE_NUMB_SEAS(J),8) = MASS_COMP(MODE_NUMB_SEAS(J),8) 
     &                                   + SSH2O_PER_SSMASS * AERO(SEAS_MAP(J)) ! [ug/m^3]
        ENDDO
        DO I=1, NWEIGHTS
          DO J=1, 7                                                  ! all components except water 
            TOT_MASS_DRY(I) = TOT_MASS_DRY(I) + MASS_COMP(I,J)       ! [ug/m^3]
          ENDDO
          TOT_MASS(I) = TOT_MASS_DRY(I) + MASS_COMP(I,8)             ! add in water [ug/m^3]
          !------------------------------------------------------------------------------------------------------------
          ! IF(I.EQ. 3) WRITE(*,'(I4,A6,9F12.6 )') I, MODE_NAME(I), TOT_MASS(I), MASS_COMP(I,:)
          ! IF(I.EQ. 8) WRITE(*,'(I4,A6,9F12.6 )') I, MODE_NAME(I), TOT_MASS(I), MASS_COMP(I,:)
          ! IF(I.EQ.16) WRITE(*,'(I4,A6,9F12.6/)') I, MODE_NAME(I), TOT_MASS(I), MASS_COMP(I,:)
          !------------------------------------------------------------------------------------------------------------
        ENDDO      

        DO I=1, NWEIGHTS
          !------------------------------------------------------------------------------------------------------------
          ! Get the ambient and dry diameter of average mass for each mode (quadrature point).
          !------------------------------------------------------------------------------------------------------------
          VOLTMP_DRY = 1.0D-30
          DO J=1, 7                                                        ! all components except water 
            VOLTMP_DRY = VOLTMP_DRY + MASS_COMP(I,J) * RECIP_DENS_COMP(J)  ! mode dry volume conc. [10^6 cm^3/m^3]
          ENDDO
          VOLTMP = VOLTMP_DRY + MASS_COMP(I,8) * RECIP_DENS_COMP(8)        ! mode ambient volume conc. [10^6 cm^3/m^3]
          DENS_MODE    (I) = TOT_MASS    (I) / VOLTMP                      ! mode ambient density [g/cm^3] 
          DENS_MODE_DRY(I) = TOT_MASS_DRY(I) / VOLTMP_DRY                  ! mode dry     density [g/cm^3] 
          DP    (I) = ( CONV_VOL_TO_DP_FAC * VOLTMP     / NI(I) )**0.333333333333333  ! [m]
          DP_DRY(I) = ( CONV_VOL_TO_DP_FAC * VOLTMP_DRY / NI(I) )**0.333333333333333  ! [m]
          IF (NI(I) .LT. 1.) THEN
             DP(I)     = DP0(I)
             DP_DRY(I) = DP0(I)
          ENDIF   
          ! IF( DP(I) .GT. DPMAX_GLOBAL ) WRITE(*,'(I4,3D15.5)') I,DP(I),NI(I),TOT_MASS(I)
          DP    (I) = MIN( MAX( DP    (I), DPMIN_GLOBAL ), DPMAX_GLOBAL )
          DP_DRY(I) = MIN( MAX( DP_DRY(I), DPMIN_GLOBAL ), DPMAX_GLOBAL )
         DIAM(IXXX,IYYY,ILAY,I)     = DP(I)   ! [m] - Store for use outside this routine.
         DIAM_dry(IXXX,IYYY,ILAY,I) = DP_DRY(I)   ! [m] - Store for use outside this routine.
          !------------------------------------------------------------------------------------------------------------
          ! Update values of KCI_COEF_DP(I) for the current diameter of average mass for each mode.
          ! THETA_POLY(I) prevents excessive condensation due to treating the mode as monodisperse.
          !   See Okuyama et al., Studies in binary nucleation: The dibutylphthalate/dioctylphthalate system,
          !   J. Chem. Phys. 89, p. 6442, 1988. 
          !------------------------------------------------------------------------------------------------------------
          INDEX_DP = NINT( LOG( DP(I) / DP_CONDTABLE_MIN ) / XLN_SCALE_DP ) + 1
          INDEX_DP = MAX( MIN( INDEX_DP, N_DP_CONDTABLE ), 1 )
!          IF( NI(I) .GT. N_MIN_DIAM_HISTOGRAM ) THEN
!            INDEX_DP_DRY = NINT( LOG( DP_DRY(I) / DP_CONDTABLE_MIN ) / XLN_SCALE_DP ) + 1
!            INDEX_DP_DRY = MAX( MIN( INDEX_DP_DRY, N_DP_CONDTABLE ), 1 )
!            DIAM_HISTOGRAM(I,INDEX_DP,    1) = DIAM_HISTOGRAM(I,INDEX_DP,    1) + 1.0D+00 
!            DIAM_HISTOGRAM(I,INDEX_DP_DRY,2) = DIAM_HISTOGRAM(I,INDEX_DP_DRY,2) + 1.0D+00 
!            !----------------------------------------------------------------------------------------------------------
!            ! IF (I.EQ.5.OR.I.EQ.6) WRITE(*,'(8F13.5)') 1D-6*NI(5),TOT_MASS(5),1D6*DP(5),1D6*DP_DRY(5),
!            ! &                                                1D-6*NI(6),TOT_MASS(6),1D6*DP(6),1D6*DP_DRY(6)
!            !----------------------------------------------------------------------------------------------------------
!          ENDIF
          KCI_COEF_DP     (I,ILAY) = THETA_POLY(I) * KCI_DP_CONDTABLE     (INDEX_DP,ILAY)
          KCI_COEF_DP_AEQ1(I,ILAY) = THETA_POLY(I) * KCI_DP_CONDTABLE_AEQ1(INDEX_DP,ILAY)
          !------------------------------------------------------------------------------------------------------------
          ! WRITE(*,90008)I,INDEX_DP,DP(I),DP_CONDTABLE(INDEX_DP-1),DP_CONDTABLE(INDEX_DP),DP_CONDTABLE(INDEX_DP+1) 
          ! WRITE(*,'(I4,D14.5,2F15.8)') I,VOLTMP,DENS_MODE(I),TOT_MASS(I) 
          ! WRITE(*,'(I4,8F12.6)') I,MASS_COMP(I,:) 
          !------------------------------------------------------------------------------------------------------------
        ENDDO
        IF( WRITE_LOG ) THEN
          WRITE(AUNIT1,'(/A8,4A25/)') 'I','DP0(I) [um]','DP [um]','DENS_MODE [g/cm^3]','DENS_MODE_DRY [g/cm^3]'
          DO I=1, NWEIGHTS
            WRITE(AUNIT1,'(I8,4D25.6)') I, DP0(I)*1.0D+06, DP(I)*1.0D+06, DENS_MODE(I), DENS_MODE_DRY(I)
          ENDDO
        ENDIF
        AVG_DP_OF_AVG_MASS_METERS = SUM( DP(:)*NI(:) ) / SUM( NI(:) )  ! [m] used in AERO_NPF, KK02 gamma expression 
      ELSE
         DO I=1,NWEIGHTS
          DIAM(IXXX,IYYY,ILAY,I) = DP(I)   ! [m] - Store for use outside this routine.
         ENDDO
      ENDIF
      
      !----------------------------------------------------------------------------------------------------------------
      ! Update the ambient and dry geometric mean diameters. 
      !----------------------------------------------------------------------------------------------------------------
      DGN(:)     = 1.0D+06 * DP    (:) * CONV_DPAM_TO_DGN(:)   ! [um]
      DGN_DRY(:) = 1.0D+06 * DP_DRY(:) * CONV_DPAM_TO_DGN(:)   ! [um]
      !----------------------------------------------------------------------------------------------------------------
      ! DO I=1, NWEIGHTS
      !   WRITE(*,'(I6,3F18.8)') I, 1.0D+06*DP_DRY(I), CONV_DPAM_TO_DGN(I), DGN_DRY(I)
      ! ENDDO
      ! WRITE(*,'(A)') 
      !----------------------------------------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------------------------
      ! Update the dry deposition velocities if desired. 
      !----------------------------------------------------------------------------------------------------------------
      IF( UPDATE_VDEP .AND. ILAY .EQ. 1 ) THEN
        CALL GET_AERO_DEPV(NMODES,TEMP_DDEP,RHOA_DDEP,LAMB_DDEP,DVIS_DDEP,
     &                     WSTR_DDEP,USTR_DDEP,RAER_DDEP,DGN,LNSIG0,DENSPI)
        IF( GET_DEP_VEL_ONLY ) RETURN
      ENDIF


      !----------------------------------------------------------------------------------------------------------------
      ! Update the coagulation coefficients, if desired. [m^3/s]
      !----------------------------------------------------------------------------------------------------------------
      SELECT CASE( UPDATE_KIJ )
      CASE ( 0 )
      CASE ( 1 )
        CALL GET_KBARNIJ(1,TK,PRES,DGN,KBAR0IJ,KBAR3IJ)
        !--------------------------------------------------------------------------------------------------------------
        ! KBAR0IJ(:,:) = 1.0D-14
        ! KBAR3IJ(:,:) = 1.0D-14
        ! KBAR3IJ(:,:) = KBAR0IJ(:,:)
        !--------------------------------------------------------------------------------------------------------------
      END SELECT


      !----------------------------------------------------------------------------------------------------------------
      !
      ! GET THE PRODUCTION AND LOSS TERMS FOR THE NUMBER CONCENTRATIONS.
      !
      !----------------------------------------------------------------------------------------------------------------
      ! Emissions terms. 
      !
      ! P_EMIS_NUMB is in [#/m^3/s].
      ! EMIS_MODE_MAP(J) is the mode number receiving the emissions species J.
      ! EMIS_MASS is the mass emissions rate [ug/m^3/s].
      ! RECIP_PART_MASS is the reciprocal of the mean mass
      !   of an emitted particle for this emissions species [1/ug].
      !----------------------------------------------------------------------------------------------------------------
      P_EMIS_NUMB(:) = 0.0D+00
      DO J=1, NEMIS_SPCS
        P_EMIS_NUMB(EMIS_MODE_MAP(J)) = P_EMIS_NUMB(EMIS_MODE_MAP(J)) 
     &                                + EMIS_MASS(J)*RECIP_PART_MASS(J)
       ENDDO 

      DIAGTMP1(1,NUMB_MAP(:)) = P_EMIS_NUMB(:)  
! Mass emissions are already treated before
      EMIS_MASS(:) = 0.0D+00 

      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A4,A30/)')'I','P_EMIS_NUMB(I) [ #/m^3/s]'
        DO I=1, NMODES
          ! WRITE(AUNIT1,'(I4,2D30.6)') I, P_EMIS_NUMB(I)
        ENDDO
      ENDIF

      !----------------------------------------------------------------------------------------------------------------
      ! Get the secondary-particle formation rate DNDT [#/m^3/s] and the 
      !   secondary-particle sulfate mass production rate DMDT_SO4 [ugSO4/m^3/s].
      !   The new particle formation rate DNDT is smaller than the nucleation 
      !   rate J and is (for most parameterizations) derived from the nucleation rate. 
      !
      ! For the calculation of DNDT and DMDT_SO4, it is assumed that all of the
      !   current H2SO4 concentration was produced by gas-phase oxidation of 
      !   SO2 during the current time step. SO4RATE is the average production rate 
      !   of H2SO4 (represented as SO4, MW=96g/mol) over the time step based on this
      !   assumption. DMDT_SO4 and DNDT are limited by the available H2SO4. 
      !
      ! When appropriate, the nucleation rate is calculated using a steady-state 
      !   H2SO4 concentration. This is done to avoid the spuriously high nucleation 
      !   rates that would result if the current H2SO4 concentration, which has 
      !   accumulated without loss over the time step, were used. When the  
      !   condensational sink (KC) is small enough such that steady-state will not
      !   be reached during the time step, an estimate of the H2SO4 concentration
      !   at the mid-point of the time step is used in nucleation and condensation
      !   calculations. 
      !
      ! The condensational sink KC is the first-order rate constant for the
      !   loss of H2SO4 due to condensation. This is summed over all modes.
      !   For the Kerminen and Kulmala (2002) parameterization for the conversion
      !   of the nucleation rate to the new particle formation rate, the 
      !   condensational sink KC_AEQ1 should obtained with the mass accommodation 
      !   coefficient set of unity. 
      !----------------------------------------------------------------------------------------------------------------
      XH2SO4_INIT = GAS( GAS_H2SO4 ) + TINYNUMER             ! input H2SO4 concentration [ugSO4/m^3]
      XH2SO4_NUCL =  XH2SO4_NUCL_MIN !TINYNUMER ! XH2SO4_NUCL_MIN              ! for the case  XH2SO4_INIT .LT. XH2SO4_NUCL_MIN
      KC = SUM( KCI_COEF_DP(:,ILAY)*NI(:) )                  ! total condensational sink [1/s] for any value of the
      KC = MAX( KC, KCMIN )                                  !   mass accommodation coefficient [1/s]
      IF( DO_NPF .AND. XH2SO4_INIT .GT. XH2SO4_NUCL_MIN ) THEN
        KC_AEQ1 = SUM( KCI_COEF_DP_AEQ1(:,ILAY)*NI(:) )      ! total condensational sink for the
        KC_AEQ1 = MAX( KC_AEQ1, KCMIN )                      !   mass accommodation coefficient set to unity  [1/s]
        XNH3 = CONVNH3 * GAS( GAS_NH3 ) * TK / PRES          ! NH3 concentration; from [ug/m^3] to [ppmV]
        SO4RATE = XH2SO4_INIT / TSTEP                        ! average H2SO4 production rate [ugSO4/m^3/s]
        IF(KC*TSTEP .GE. XNTAU ) THEN                        ! invoke steady-state assumption
          IH2SO4_PATH = 1
          XH2SO4_SS = MIN( SO4RATE/KC, XH2SO4_INIT )         ! steady-state H2SO4 [ugSO4/m^3]
          CALL STEADY_STATE_H2SO4(PRES,RH,TK,XH2SO4_SS,SO4RATE,XNH3,KC,TSTEP,XH2SO4_SS_WNPF)
          XH2SO4_NUCL = XH2SO4_SS_WNPF                       ! [H2SO4] for nucl., GR, and cond. calculation [ugSO4/m^3]
        ELSE
          IH2SO4_PATH = 2
          XH2SO4_NUCL = SO4RATE / ( (2.0D+00/TSTEP) + KC )   ! use [H2SO4] at mid-time step [ugSO4/m^3]
        ENDIF
        CALL NPFRATE(PRES,RH,TK,XH2SO4_NUCL,SO4RATE,XNH3,KC_AEQ1,DNDT,DMDT_SO4,0)
        ! WRITE(34,90010) TK,RH,UGM3_NCM3*XH2SO4_INIT,UGM3_NCM3*XH2SO4_NUCL,KC,XNH3,1.0D-06*DNDT,IH2SO4_PATH
        IF( WRITE_LOG ) THEN
          WRITE(AUNIT1,'(/A)')'NEW PARTICLE FORMATION'                          
          WRITE(AUNIT1,'(/A)')'PRES,RH,TK,XH2SO4_INIT(ug/m3),XH2SO4_NUCL(#/cm3),SO4RATE,XNH3,KC,DNDT,DMDT_SO4'
          WRITE(AUNIT1,90003)  PRES,RH,TK,XH2SO4_INIT,       XH2SO4_NUCL*UGM3_NCM3,SO4RATE,XNH3,KC,DNDT,DMDT_SO4
          WRITE(AUNIT1,'(/A)')'TK,RH,H2SO4(#.cm^3),H2SO4_NUCL(#/cm3),KC,NH3,DNDT(#/cm3/s),IH2SO4_PATH'      
          WRITE(AUNIT1,90010)  TK,RH,UGM3_NCM3*XH2SO4_INIT,UGM3_NCM3*XH2SO4_NUCL,KC,XNH3,1.0D-06*DNDT,IH2SO4_PATH
        ENDIF
      ELSE
        DNDT     = 0.0D+00
        DMDT_SO4 = 0.0D+00
      ENDIF


      !----------------------------------------------------------------------------------------------------------------
      ! Get the B_i loss       terms due to intermodal coagulation. [1/s]
      ! Get the R_i production terms due to intermodal coagulation. [#/m^3/s]
      ! Get the C_i terms, which include all source terms.          [#/m^3/s]
      ! For the C_i terms, the secondary particle formation term DNDT must be 
      !   added in after coupling to condensation below.
      ! The A_i terms for intramodal coagulation are directly computed 
      !   from the coagulation coefficients when the number equations 
      !   are integrated.
      ! If DIKL(I,K,L) = 0, then modes K and L to not coagulate to form mode I.
      ! DIJ(I,J) is unity if coagulation of mode I with mode J results
      !   in the removal of particles from mode I, and zero otherwise. 
      !----------------------------------------------------------------------------------------------------------------
      BI(:) = 0.0D+00
      RI(:) = 0.0D+00
      do ikl = 1, nDIKL
        i = dikl_control(ikl)%i
        k = dikl_control(ikl)%k
        l = dikl_control(ikl)%l
        RI(i) = RI(i) + KBAR0IJ(k,l) * ni(k) * ni(l)
      end do
      do k=1, NWEIGHTS
        do i = 1, NWEIGHTS
          BI(i) = BI(i) + KBAR0IJ(i,k)*NI(k)*xDIJ(i,k) ! Coagulation of modes I and J removes
        end do
      end do
      CI(:) = RI(:) + P_EMIS_NUMB(:)                              ! all source terms for number concentration 
      DIAGTMP1(3,NUMB_MAP(:)) = RI(:)  

      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/6A20/)') 'I','BI [1/s]','CI [#/m^3/s]','RI [#/m^3/s]','NI [#/m^3]','KII[m^3/s]'
        DO I=1, NWEIGHTS
          WRITE(AUNIT1,90004) I, BI(I), CI(I), RI(I), NI(I), KBAR0IJ(I,I)
        ENDDO
      ENDIF


      !----------------------------------------------------------------------------------------------------------------
      !
      ! GET THE PRODUCTION AND LOSS TERMS FOR THE MASS CONCENTRATIONS.
      !
      !----------------------------------------------------------------------------------------------------------------
      ! The production terms PIQ(NWEIGHTS,NMASS_SPCS) are in [ug/m^3/s].
      !   The first  index runs over all modes or quadrature points.
      !   The second index runs over the NMASS_SPCS (=5) principal mass
      !     species: sulfate, BC, OC, dust, sea salt. 
      !
      ! Get the emissions production terms (Pemis_i,q in the manuscript)
      !   in [ug/m^3/s] and put them in the total production rate array PIQ
      !   (the P_i,q in the manuscript).
      !
      ! EMIS_MODE_MAP and EMIS_SPCS_MAP have elements corresponding to 
      !   the aerosol types (in this order): AKK(=1), ACC(=2), BCC(=8), OCC(=7),
      !               DD1(=3), SSA(=5), SSC(=6), BOC(BC=8), BOC(OC=9), DD2(=10).
      ! EMIS_MODE_MAP(J) is mode number receiving the emissions held 
      !                  in EMIS_MASS(J).
      ! EMIS_SPCS_MAP(J) is the chemical species number (1-5) of the species
      !                  held in EMIS_MASS(J).
      ! Currently, EMIS_SPCS_MAP = (/1,1,2,3,4,5,5,2,3,4/)
      !----------------------------------------------------------------------------------------------------------------
      PIQ(:,:) = 0.0D+00          ! Zero all production terms for this time step.
      DO J=1, NEMIS_SPCS          ! Loop over the emitted species. 
        PIQ( EMIS_MODE_MAP(J), EMIS_SPCS_MAP(J) ) = 
     &  PIQ( EMIS_MODE_MAP(J), EMIS_SPCS_MAP(J) ) + EMIS_MASS(J)
      ENDDO
      DO I=1, NMODES
        DO Q=1, NM(I)
          DIAGTMP1(8,MASS_MAP(I,Q)) = PIQ(I,PROD_INDEX(I,Q))  
        ENDDO
      ENDDO      


      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A/)') 'EMISSIONS ONLY: I, PIQ(I,:)'
        DO I=1, NMODES
          ! WRITE(AUNIT1,'(I5,5D14.4)') I, PIQ(I,:)
        ENDDO
      ENDIF


      !-------------------------------------------------------------------------------------------------------------------
      ! Get the Pcoag_i,q terms for the production of mass species Q in mode or quadrature point I 
      !   due to coagulation between modes or quadrature points K and L.
      ! These are at the same time mapped to the PIQ array.
      !
      ! QQ = PROD_INDEX(I,Q) is the principal species index that is the second index of the PIQ and MJQ arrays.
      ! It points to chemical species CHEM_SPC_NAME(Q) where
      !   CHEM_SPC_NAME(NMASS_SPCS) = (/'SULF','BCAR','OCAR','DUST','SEAS'/)
      !-------------------------------------------------------------------------------------------------------------------
      ! WRITE(36,'(A)')'I,Q,K,L,PIQ(I,QQ),KBAR3IJ(K,L)*NI(K)*NI(L)*(MJQ(K,Q)+MJQ(L,Q))'
      ! WRITE(36,'(A)')'I,Q,K,L,KBAR3IJ(K,L),NI(K),NI(L),MJQ(K,Q),MJQ(L,Q)'
      !-------------------------------------------------------------------------------------------------------------------

      PIQTMP(:,:) = 0.0D+00

        !-----------------------------------------------------------------------------------------------------------------
        ! Sum over all K-L interactions and identify which interactions produce species Q in mode I.
        !
        ! There are three cases where mass of Q is transferred to mode I, if the Q concentration 
        ! in the donor mode(s) K and/or L is non-zero.
        !
        !   1. Mode I is the same as mode K but different from mode L, and Q-mass from mode L is 
        !      transferred to mode I.
        !   2. Mode I is the same as mode L but different from mode K, and Q-mass from mode K is
        !      transferred to mode I.
        !   3. Mode I is different from either mode K or mode L, and Q-mass from both mode K and mode L is 
        !      transferred to mode I.
        !
        ! The K-L double sum are symmetric in K-L, so we can sum over either the 
        ! superdiagonal or subdiagonal part of the 'matrix', but not both. Note that 
        ! KBAR3IJ(K,L) is not symmetric in K-L, so KBAR3IJ(K,L) and KBAR3IJ(L,K) are not interchangeable.
        !---------------------------------------------------------------------------------------------------------------

      do I=1, NWEIGHTS          ! loop over all modes or quadrature points receiving mass species Q
        do klq = 1, giklq_control(i)%n
          k = giklq_control(i)%k(klq)
          l = giklq_control(i)%l(klq)
          qq = giklq_control(i)%qq(klq)

          PIQTMP(i,qq) = PIQTMP(i,qq) + 
     &         ni(k)*ni(l) *KBAR3IJ(l,k) * MJQ(l,qq)
        end do
      end do

          !-------------------------------------------------------------------------------------------------------------
          ! IF( I.EQ.15 .AND. K.LE.2 .AND. L.EQ.10 ) THEN
          !   WRITE(36,'(4I3,3D15.6,I6)') I,Q,K,L,PIQTMP(I,QQ),NI(K),NI(L),IBRANCH
          !   WRITE(36,'(12X,5D15.6)')             KBAR3IJ(K,L),KBAR3IJ(L,K),MJQ(K,QQ),MJQ(L,QQ)
          !   WRITE(36,'(12X,5D15.6)')             KBAR3IJ(K,L)*MJQ(K,QQ),KBAR3IJ(L,K)*MJQ(L,QQ)
          !   WRITE(36,'(12X,5D15.6)')             KBAR3IJ(K,L)*MJQ(K,QQ)+KBAR3IJ(L,K)*MJQ(L,QQ)
          !   WRITE(36,'(12X,5D15.6)')NI(K)*NI(L)*(KBAR3IJ(K,L)*MJQ(K,QQ)+KBAR3IJ(L,K)*MJQ(L,QQ))
          ! ENDIF
          !-------------------------------------------------------------------------------------------------------------
      DO I=1, NMODES
        DO Q=1, NM(I)
          DIAGTMP1( 8,MASS_MAP(I,Q)) = PIQ   (I,PROD_INDEX(I,Q))  ! Due to mass emissions
          DIAGTMP1(10,MASS_MAP(I,Q)) = PIQTMP(I,PROD_INDEX(I,Q))  ! Due to coagulation
        ENDDO
      ENDDO      
      PIQ(:,:) = PIQ(:,:) + PIQTMP(:,:)

      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A/)') 'I, Q, QQ, PIQ(I,QQ) [ug/m^3/s] - after coagulation mass terms'
        DO I=1, NWEIGHTS
        DO Q=1, NM(I)
          QQ = PROD_INDEX(I,Q)
          ! WRITE(AUNIT1,90000) I, Q, QQ, PIQ(I,QQ)
        ENDDO
        ENDDO
      ENDIF

      !----------------------------------------------------------------------------------------------------------------
      ! Get the F_i terms for the loss of species Q in quadrature point I due
      !   to intermodal coagulation of I with other quadrature points. 
      ! As all species in quadrature point I are lost through this coagulation,
      !   the species index Q is not needed.
      ! DIJ(I,J) is unity if coagulation of mode I with mode J results
      !   in the removal of particles from mode I, and zero otherwise. 
      !   DIJ(I,J) is not symmetric in I-J.
      ! KBAR3IJ(I,J) is not symmetric in I-J, and the first index I refers
      !   to the 'donor' mode in the coagulation interaction. Since mode I
      !   is here losing mass due to coagulation, mode I is the 'donor' mode.
      !----------------------------------------------------------------------------------------------------------------
      ! WRITE(35,'(/A/)') 'DIJ(I,J), I, J, KBAR3IJ(I,J), NI(J), FI(I)'
      !----------------------------------------------------------------------------------------------------------------
      FI(:) = 0.0D+00
      DO I=1, NWEIGHTS
        DO J=1, NWEIGHTS    ! For each mode I, we must sum over all modes J.
          IF( DIJ(I,J).GT.0 ) FI(I) = FI(I) + KBAR3IJ(I,J)*NI(J)
          !------------------------------------------------------------------------------------------------------------
          ! WRITE(35,'(3I4,3D15.7)') DIJ(I,J), I, J, KBAR3IJ(I,J), NI(J), FI(I)
          !------------------------------------------------------------------------------------------------------------
        ENDDO
      ENDDO
      DO I=1, NMODES
        DO Q=1, NM(I)
          DIAGTMP1(13,MASS_MAP(I,Q)) = - FI(I) * AERO( MASS_MAP(I,Q) )  
        ENDDO
      ENDDO

      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A/)') 'I, FI [1/s], NI [#/m^3]'
        DO I=1, NWEIGHTS
          ! WRITE(AUNIT1,90001) I, FI(I), NI(I)
        ENDDO
      ENDIF


      !----------------------------------------------------------------------------------------------------------------
      ! Get the Pcloud_i,q terms due to the in-cloud production of sulfate
      ! for each mode or quadrature point in [ugSO4/m^3/s].
      !
      ! AQSO4RATE is the in-cloud sulfate production rate [ug/m^3/s].
      ! 
      ! KAPPAI(I) is the soluble (or activating) fraction for mode I, as specified in aero_param.f.
      ! 
      ! If all solule particles in mode I are assumed to activate, then
      ! the product KAPPAI(I) * NI(I) * RSUM_ACTIV is the fraction of the
      ! in-cloud sulfate going into mode I.
      ! 
      ! If a droplet activation calculation is done for each mode, then 
      ! the product NACT(I) * RSUM_ACTIV is the fraction of the
      ! in-cloud sulfate going into mode I. The geometric mean radii passed to
      ! GETACTFRAC is for the dry aerosol. 
      !----------------------------------------------------------------------------------------------------------------
      PIQTMP(:,:) = 0.0D+00
      IF( AQSO4RATE .GT. AQSO4RATE_MIN ) THEN  
        IF(     ACTIVATION_SCHEME .EQ. 1 ) THEN  
          NSOL(:) = KAPPAI(:)*NI(:)
          RSUM_ACTIV = 1.0D+00 / SUM( NSOL(:) + TINYDENOM ) 
          PIQTMP(:,PROD_INDEX_SULF) = ( KAPPAI(:)*NI(:)*RSUM_ACTIV ) * AQSO4RATE
          NACTV(IXXX,IYYY,ILAY,:) = NSOL(:)       ! [#/m] - Store for use outside this routine.
          !------------------------------------------------------------------------------------------------------------
          ! WRITE(40,'(/A,F15.6/)')'Total number soluble (#/cm^3) = ', 1.0D-06/RSUM_ACTIV 
          ! DO I=1, NMODES
          !   WRITE(40,'(I3,F5.2,3F12.4,A5)')I,KAPPAI(I),1.0D-06*NI(I),1.0D-06*NSOL(I),NSOL(I)*RSUM_ACTIV,MODE_NAME(I)
          ! ENDDO  
          !------------------------------------------------------------------------------------------------------------
        ELSEIF( ACTIVATION_SCHEME .EQ. 2 ) THEN
          DO I=1, NMODES
            MI5(I,:) = MJQ(I,:) * NI(I)     ! mass conc. of each component in mode I [ug/m^3]
          ENDDO  
          CALL GETACTFRAC(NMODES,NI,MI5,0.5D+00*DGN_DRY,SIG0,TK,PRES,WUPDRAFT,AC,FRACACTN,FRACACTM,NACT,MACT)
          RSUM_ACTIV = 1.0D+00 / SUM ( NACT(:) + TINYDENOM ) 
          PIQTMP(:,PROD_INDEX_SULF) = ( NACT(:)*RSUM_ACTIV ) * AQSO4RATE
          NACTV(IXXX,IYYY,ILAY,:) = NACT(:)       ! [#/m] - Store for use outside this routine.
          !------------------------------------------------------------------------------------------------------------
          ! WRITE(40,'(/A,F15.6/)')'Total number activated (#/cm^3) = ', 1.0D-06/RSUM_ACTIV 
          ! DO I=1, NMODES
          !   WRITE(40,90009)I,FRACACTN(I),1.0D-06*NI(I),1.0D-06*NACT(I),1.0D-06*NACT(I),NACT(I)*RSUM_ACTIV,
          !&                       MODE_NAME(I),DGN(I),DGN_DRY(I)
          ! ENDDO  
          !------------------------------------------------------------------------------------------------------------
        ENDIF
      ELSE                                       ! Put all in-cloud sulfate into the accumulation mode.
        IF( NUMB_AKK_1 .NE. 0 ) THEN             ! Mode AKK exists and mode ACC has mode number = 2.
          PIQTMP(2,PROD_INDEX_SULF) = AQSO4RATE  ! [ugSO4/m^3/s]
        ELSE                                     ! Mode AKK does not exist and mode ACC has mode number = 1.
          PIQTMP(1,PROD_INDEX_SULF) = AQSO4RATE  ! [ugSO4/m^3/s]
        ENDIF      
      ENDIF
      PIQ(:,PROD_INDEX_SULF) = PIQ(:,PROD_INDEX_SULF) + PIQTMP(:,PROD_INDEX_SULF) 
      DIAGTMP1(12,MASS_MAP(:,PROD_INDEX_SULF)) = PIQTMP(:,PROD_INDEX_SULF)  

      !----------------------------------------------------------------------------------------------------------------
      ! Get the Pgrowth_i,q terms due to condensation and gas-particle 
      !   mass transfer for each mode (or quadrature point) and add them to the PIQ array. 
      !
      ! The net loss of H2SO4 due to both secondary particle formation
      !   and condensation should not exceed the current H2SO4 concentration.
      !   This is enforced by rescaling the two H2SO4 consumption rates in a way 
      !   that preserves the relative magnitudes of these two loss processes. 
      !   The net condensation rate PQ_GROWTH is calculated from XH2SO4_NUCL, not
      !   the total accumulated H2SO4 in XH2SO4_INIT, for balance with the 
      !   treatment of new particle formation above. 
      !
      ! The expression in parentheses on the rhs of PIQ is that for h_i, 
      !
      !   h_i = KCI_COEF_DP(i,ILAY) * NI(i) / KC
      !
      !   the ratio of the condensational sink of mode or quadrature point I 
      !   to the total condensational sink.
      !
      ! At this point, XH2SO4_INIT = GAS( GAS_H2SO4 ) + TINYNUMER.
      !
      ! The H2SO4 concentration GAS( GAS_H2SO4 ) is updated here.
      !----------------------------------------------------------------------------------------------------------------
      PQ_GROWTH = XH2SO4_NUCL * ( 1.0D+00 - EXP(-KC*TSTEP) ) / TSTEP            ! [ugSO4/m^3/s]
      TOT_H2SO4_LOSS = ( DMDT_SO4 + PQ_GROWTH ) * TSTEP                         ! [ugSO4/m^3]
      IF ( TOT_H2SO4_LOSS .GT. XH2SO4_INIT ) THEN                               ! XH2SO4_INIT=GAS(GAS_H2SO4)+TINYNUMER
        DMDT_SO4  = DMDT_SO4  * ( XH2SO4_INIT / ( TOT_H2SO4_LOSS + TINYDENOM ) )! [ugSO4/m^3/s]
        DNDT      = DNDT      * ( XH2SO4_INIT / ( TOT_H2SO4_LOSS + TINYDENOM ) )! [  #  /m^3/s]
        PQ_GROWTH = PQ_GROWTH * ( XH2SO4_INIT / ( TOT_H2SO4_LOSS + TINYDENOM ) )! [ugSO4/m^3/s]
        GAS( GAS_H2SO4 ) = TINYNUMER
      ELSE
        GAS( GAS_H2SO4 ) = GAS( GAS_H2SO4 ) - TOT_H2SO4_LOSS + TINYNUMER
      ENDIF
      DIAGTMP1(2,NUMB_MAP(1))               = DNDT  
      DIAGTMP1(9,SULF_MAP(PROD_INDEX_SULF)) = DMDT_SO4  
      
      ! WRITE(35,'(8D13.5)')PIQ(:,PROD_INDEX_SULF), KCI_COEF_DP(:,ILAY),NI(:),KC
      
      PIQTMP(:,PROD_INDEX_SULF) = ( KCI_COEF_DP(:,ILAY)*NI(:)/KC ) * PQ_GROWTH
      PIQ   (:,PROD_INDEX_SULF) = PIQ(:,PROD_INDEX_SULF) + PIQTMP(:,PROD_INDEX_SULF) 

      DIAGTMP1(11,MASS_MAP(:,PROD_INDEX_SULF)) = PIQTMP(:,PROD_INDEX_SULF)  
      !-----------------------------------------------------------------------------------------------------------------
      ! Add the secondary particle formation term DNDT [#/m^3/s] for the number concentration term.
      ! Add the secondary particle formation term DMDT_SO4 [ug/m^3/s], calculated above, 
      !   to the total mass production rate array PIQ for the AKK mode.
      ! If the AKK mode is absent, the secondary particle formation terms still go into mode 1.
      !-----------------------------------------------------------------------------------------------------------------
      CI(1) = CI(1) + DNDT                                        ! add secondary particle formation number term
      PIQ(1,PROD_INDEX_SULF) = PIQ(1,PROD_INDEX_SULF) + DMDT_SO4  ! add secondary particle formation mass   term

      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A,5X,3D15.8)')'XH2SO4_INIT, XH2SO4_NUCL, PQ_GROWTH = ', XH2SO4_INIT, XH2SO4_NUCL, PQ_GROWTH
        WRITE(AUNIT1,*)'PIQ(1,PROD_INDEX_SULF) = ', PIQ(1,PROD_INDEX_SULF)
      ENDIF


      !----------------------------------------------------------------------------------------------------------------
      !
      ! SOLVE FOR THE UPDATED NUMBER CONCENTRATIONS (QUADRATURE WEIGHTS). 
      !
      !----------------------------------------------------------------------------------------------------------------
      DO I=1, NWEIGHTS
        Y0 = NI(I)                                  ! Initial number concentration [#/m^3/s].
        A = 0.5D+00*KBAR0IJ(I,I)
        B = BI(I)   
        C = CI(I)   
        IF( C .GT. 1.0D-30 ) THEN
          DELTA = SQRT( B * B + 4.0D+00 * A * C )    
          R1 = 2.0D+00 * A * C / ( B + DELTA )       
          R2 = - 0.5D0 * ( B + DELTA )             
          GAMMA =  - ( R1 - A * Y0 ) / ( R2 - A * Y0 )
          GEXPDT = GAMMA * EXP( - DELTA * TSTEP )
          Y = (  R1 + R2 * GEXPDT ) / ( A * ( 1.0D+00 + GEXPDT ) )
        ELSE                                        ! When C = 0.0D+00, as we assume C is not negative.
          EXPDT = EXP( - B * TSTEP )
          IF( 1.0D+00-EXPDT .GT. PIQ_THRESH ) THEN                     ! IF( EXPDT .LT. 1.0D+00 ) THEN
            Y = B * Y0 * EXPDT / ( B + A * Y0 * ( 1.0D+00 - EXPDT ) )
          ELSE
            Y = Y0 / ( 1.0D+00 + A * Y0 * TSTEP )
          ENDIF
        ENDIF
        AERO( NUMB_MAP(I) ) = MAX ( Y, MINCONC )   ! update the output array, not the work array NI(I)
        DIAGTMP1(4,NUMB_MAP(I)) = -B*Y0  
        DIAGTMP1(5,NUMB_MAP(I)) = -A*Y0*Y0
        ! WRITE(*,'(4D15.5)') A,B,C,Y0  
      ENDDO

      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A6,2A30/)') 'I','NI(I) [#/m^3]','AERO( NUMB_MAP(I) ) [#/m^3]'
        DO I=1, NWEIGHTS
          ! WRITE(AUNIT1,90006) I, NI(I), AERO ( NUMB_MAP(I) )
        ENDDO
      ENDIF


      !----------------------------------------------------------------------------------------------------------------
      !
      ! SOLVE FOR THE UPDATED MASS CONCENTRATIONS. 
      !
      !----------------------------------------------------------------------------------------------------------------
      ! Update the sulfate, BC, OC, dust, and sea salt concentrations.
      !
      ! MASS_MAP(I,Q) is the location in AERO(:) of the Qth mass in mode I.
      ! Mode I has NM(I) mass species defined for it, and NM(I) varies between
      !   1 and NMASS_SPCS (=5).
      !
      ! The second index of PROD_INDEX has NMASS_SPCS (=5) values:
      !   1=sulfate, 2=BC, 3=OC, 4=dust, 5=sea salt.
      !   PROD_INDEX(I,Q) is the location in array PIQ(I,Q) of chemical species
      !   CHEM_SPC_NAME(Q) for mode (quadrature weight) I.
      !
      ! IF ( 1.0D+00-EXPDT .GT. PIQ_THRESH ) --> The first (loss) term in AERO update may be significant.
      !----------------------------------------------------------------------------------------------------------------
      DO I=1, NWEIGHTS
        EXPDT = EXP( - FI(I) * TSTEP ) 
        IF ( 1.0D+00-EXPDT .GT. PIQ_THRESH ) THEN 
          FACTOR = ( 1.0D+00 - EXPDT ) / FI(I)
          DO Q=1, NM(I)
            AERO(MASS_MAP(I,Q)) = AERO(MASS_MAP(I,Q)) * EXPDT + PIQ(I,PROD_INDEX(I,Q)) * FACTOR
          ENDDO
        ELSE
          DO Q=1, NM(I)
            AERO(MASS_MAP(I,Q)) = AERO(MASS_MAP(I,Q)) + PIQ(I,PROD_INDEX(I,Q)) * TSTEP
          ENDDO
        ENDIF
      ENDDO


      !----------------------------------------------------------------------------------------------------------------
      ! Update the nitrate, ammonium and aerosol water concentrations.
      !
      ! Call the thermodynamic module again to determine the bulk gas-particle 
      ! partitioning of the inorganic species and the water content
      ! associated with those species. Also determine the water content
      ! associated with the NaCl of sea salt.
      !
      ! Note that AERO(MASS_H2O) must include the sea-salt associated water 
      !   upon exit from this routine. 
      ! Upon return from AERO_THERMO, AERO(MASS_H2O) contains only the 
      !   non-sea salt-associated water. The sea salt-associated water is in SSH2O.
      !----------------------------------------------------------------------------------------------------------------
      TOT_SULF = SUM( AERO(SULF_MAP(:)) )
      TOT_DUST = SUM( AERO(DUST_MAP(:)) )
      TOT_SEAS = SUM( AERO(SEAS_MAP(:)) )
      IF( WRITE_LOG ) WRITE(AUNIT1,'(/A,F12.4/)') 'TOT_SULF = ', TOT_SULF
      AERO_WATER_ACTUAL = AERO(MASS_H2O) + SSH2O
      CALL AERO_THERMO(TOT_SULF,AERO(MASS_NO3),AERO(MASS_NH4),AERO(MASS_H2O),GAS(GAS_NH3),
     &                 GAS(GAS_HNO3),TOT_DUST,TOT_SEAS,SSH2O,TK,RH,PRES,RHD,RHC)
      AERO_WATER_WET = AERO(MASS_H2O) + SSH2O

      !----------------------------------------------------------------------------------------------------------------
      ! Again, adjust the water concentration for hysteresis.
      !----------------------------------------------------------------------------------------------------------------
      IF ( RH .GT. RHC  .AND.  RH .LT. RHD ) THEN
        IF ( AERO_WATER_WET .GT. 0.0D+00 ) THEN
          IF ( AERO_WATER_ACTUAL/AERO_WATER_WET .LT. 0.5D+00 ) THEN
            AERO(MASS_H2O) = 0.0D+00  ! Zero the non-sea salt-associated water.
            SSH2O = 0.0D+00           ! Zero the     sea salt-associated water.
          ENDIF
        ENDIF
      ELSEIF ( RH .LE. RHC ) THEN  
        AERO(MASS_H2O) = 0.0D+00
        SSH2O = 0.0D+00
      ENDIF
      !----------------------------------------------------------------------------------------------------------------
      ! The total aerosol water must exit the routine in AERO(MASS_H2O).
      !----------------------------------------------------------------------------------------------------------------
      AERO(MASS_H2O) = AERO(MASS_H2O) + SSH2O         ! Total aerosol water now.


      !----------------------------------------------------------------------------------------------------------------
      ! Transfer mass and number concentrations from DD1 to DS1, DD2 to DS2,
      !   BC1 to BC2, and BC2 to BC3, if the appropriate modes are defined 
      !   in the present configuration.
      !
      ! This transfer is based upon the volume fraction of inorganic constituents.
      !   For computational efficiency, the maximum inorganic volume fraction (MIVF)
      !   for a mode is transformed into the maximum inorganic mass ratio (MIMR), 
      !   and the MIMR are precomputed parameters.
      !
      ! In the following IF statements, the first (and lengthy) quantity is the current
      !   inorganic-to-dust or inorganic-to-BC mass ratio. For example, in the second IF statement below,  
      !
      !     AERO(MASS_DD1_SULF)*OPTOT_NO3NH4H2O_TO_SULF 
      !
      !   is the total mass concentration of inorganic coating 
      !   (sulfate + nitrate + ammonium + water) in mode DD1, and 
      !
      !     AERO(MASS_DD1_SULF)*OPTOT_NO3NH4H2O_TO_SULF/AERO(MASS_DD1_DUST) is the 
      !
      !   inorganic-to-dust mass ratio for mode DD1.
      !
      ! There is no need to explicitly transfer nitrate, ammonium, or aerosol water. 
      !
      ! OPTOT_NO3NH4H2O_TO_SULF is the total NO3+NH4+H2O mass per unit mass SO4, plus 1.
      !----------------------------------------------------------------------------------------------------------------
      AEROTMP2(:) = AERO(:) 
      AERO(MASS_DD1_DUST) = MAX( AERO(MASS_DD1_DUST), TINYNUMER )
      AERO(MASS_BC1_BCAR) = MAX( AERO(MASS_BC1_BCAR), TINYNUMER )
      AERO(MASS_BC2_BCAR) = MAX( AERO(MASS_BC2_BCAR), TINYNUMER )
      IF( MASS_DD2_DUST .GT. 0 ) AERO(MASS_DD2_DUST) = MAX( AERO(MASS_DD2_DUST), TINYNUMER )   ! If mode DD2 exists. 

      TOT_SULF = SUM( AERO(SULF_MAP(:)) ) + TINYNUMER
      OPTOT_NO3NH4H2O_TO_SULF = 1.0D+00 + SUM(AERO(1:3)) / TOT_SULF

      IF( AERO(MASS_DD1_SULF)*OPTOT_NO3NH4H2O_TO_SULF/AERO(MASS_DD1_DUST) .GT. MIMR_DDD ) THEN
        !--------------------------------------------------------------------------------------------------------------
        ! Transfer mode DD1 to mode DS1.
        !--------------------------------------------------------------------------------------------------------------
        AERO(MASS_DS1_SULF) = AERO(MASS_DS1_SULF) + AERO(MASS_DD1_SULF)
        AERO(MASS_DS1_DUST) = AERO(MASS_DS1_DUST) + AERO(MASS_DD1_DUST)
        AERO(NUMB_DS1_1   ) = AERO(NUMB_DS1_1   ) + AERO(NUMB_DD1_1   )
        AERO(MASS_DD1_SULF) = TINYNUMER
        AERO(MASS_DD1_DUST) = TINYNUMER
        AERO(NUMB_DD1_1   ) = TINYNUMER
      ENDIF
  
      IF( MASS_DD2_DUST .GT. 0.0D+00 ) THEN
        IF( AERO(MASS_DD2_SULF)*OPTOT_NO3NH4H2O_TO_SULF/AERO(MASS_DD2_DUST) .GT. MIMR_DDD ) THEN
          !------------------------------------------------------------------------------------------------------------
          ! Transfer mode DD2 to mode DS2.
          !------------------------------------------------------------------------------------------------------------
          AERO(MASS_DS2_SULF) = AERO(MASS_DS2_SULF) + AERO(MASS_DD2_SULF)
          AERO(MASS_DS2_DUST) = AERO(MASS_DS2_DUST) + AERO(MASS_DD2_DUST)
          AERO(NUMB_DS2_1   ) = AERO(NUMB_DS2_1   ) + AERO(NUMB_DD2_1   )
          AERO(MASS_DD2_SULF) = TINYNUMER
          AERO(MASS_DD2_DUST) = TINYNUMER
          AERO(NUMB_DD2_1   ) = TINYNUMER
        ENDIF
      ENDIF

      IF( AERO(MASS_BC1_SULF)*OPTOT_NO3NH4H2O_TO_SULF/AERO(MASS_BC1_BCAR) .GT. MIMR_BC1 ) THEN
        !--------------------------------------------------------------------------------------------------------------
        ! Transfer mode BC1 to mode BC2.
        !--------------------------------------------------------------------------------------------------------------
        AERO(MASS_BC2_SULF) = AERO(MASS_BC2_SULF) + AERO(MASS_BC1_SULF)
        AERO(MASS_BC2_BCAR) = AERO(MASS_BC2_BCAR) + AERO(MASS_BC1_BCAR)
        AERO(NUMB_BC2_1   ) = AERO(NUMB_BC2_1   ) + AERO(NUMB_BC1_1   )
        AERO(MASS_BC1_SULF) = TINYNUMER
        AERO(MASS_BC1_BCAR) = TINYNUMER
        AERO(NUMB_BC1_1   ) = TINYNUMER
      ENDIF

      IF( AERO(MASS_BC2_SULF)*OPTOT_NO3NH4H2O_TO_SULF/AERO(MASS_BC2_BCAR) .GT. MIMR_BC2 ) THEN
        !--------------------------------------------------------------------------------------------------------------
        ! Transfer mode BC2 to mode BC3.
        !--------------------------------------------------------------------------------------------------------------
        IF( INCLUDE_BC3 ) THEN
          AERO(MASS_BC3_SULF) = AERO(MASS_BC3_SULF) + AERO(MASS_BC2_SULF)
          AERO(MASS_BC3_BCAR) = AERO(MASS_BC3_BCAR) + AERO(MASS_BC2_BCAR)
          AERO(NUMB_BC3_1   ) = AERO(NUMB_BC3_1   ) + AERO(NUMB_BC2_1   )
          AERO(MASS_BC2_SULF) = TINYNUMER
          AERO(MASS_BC2_BCAR) = TINYNUMER
          AERO(NUMB_BC2_1   ) = TINYNUMER
        ENDIF
      ENDIF

      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A,2F15.6)')'MIVF_DDD, MIMR_DDD = ', MIVF_DDD, MIMR_DDD 
        WRITE(AUNIT1,'( A,2F15.6)')'MIVF_BC1, MIMR_BC1 = ', MIVF_BC1, MIMR_BC1 
        WRITE(AUNIT1,'( A,2F15.6)')'MIVF_BC2, MIMR_BC2 = ', MIVF_BC2, MIMR_BC2 
        WRITE(AUNIT1,'(/A,2D15.6)')'TOT_SULF, OPTOT_NO3NH4H2O_TO_SULF = ', TOT_SULF, OPTOT_NO3NH4H2O_TO_SULF
      ENDIF

      !----------------------------------------------------------------------------------------------------------------
      ! Intermodal transfer from the Aitken (AKK) mode to the accumulation (ACC) mode.
      !
      ! This is not a physical process, but a reclassification of particles. See Binkowski and Roselle (2003).
      !
      ! AERO(MASS_AXX_SULF)*OPTOT_NO3NH4H2O_TO_SULF is the total mass in mode AXX (X=K, C) [ug/m^3]. 
      !   Division by AERO(NUMB_AXX_1) converts to mean mass per particle [ug].
      !   Multiplication by CONV_MASS_TO_DP converts to Dp^3 [m^3].
      !   Taking the cube root yields the diameter of average mass for the mode [m].
      !
      ! FNUM is the fraction of the AKK mode mass and number concentrations transferred to the ACC mode
      !      during this time step. Binkowski and Roselle (2003) limit FNUM to a maximum of 0.5 for
      !      numerical stability, and the same is done here.
      !----------------------------------------------------------------------------------------------------------------
      IF( INTERMODAL_TRANSFER .AND. AERO( NUMB_MAP(1) ) .GT. AKK_MINNUM_IMTR ) THEN
        DPAKK = ( CONV_MASS_TO_DP * AERO( MASS_AKK_SULF )                ! diameter of average mass for AKK [m]
     &          * OPTOT_NO3NH4H2O_TO_SULF / AERO( NUMB_AKK_1 ) )**0.33333333333
        DPACC = ( CONV_MASS_TO_DP * AERO( MASS_ACC_SULF )                ! diameter of average mass for ACC [m]
     &          * OPTOT_NO3NH4H2O_TO_SULF / AERO( NUMB_ACC_1 ) )**0.33333333333
        IF(     IMTR_METHOD .EQ. 1 ) THEN
          !------------------------------------------------------------------------------------------------------------
          ! Calculate the fraction transferred based on the relative difference 
          !   in mass mean diameters of the AKK and ACC modes. 
          !------------------------------------------------------------------------------------------------------------
          IF( DPAKK .GE. DPAKK0 ) THEN                                     ! [m]
            FNUM = ( ( DPAKK - DPAKK0 ) / ( DPACC - DPAKK0 ) )**IMTR_EXP   ! fraction transferred from AKK to ACC
            FNUM = MAX( MIN( FNUM, FNUM_MAX ), 0.0D+00 )                   ! limit transfer in a single transfer
          ELSE
            FNUM = 0.0D+00
          ENDIF
          F3 = FNUM
          ! WRITE(34,'(7D12.4)')FNUM,F3
        ELSEIF( IMTR_METHOD .EQ. 2 ) THEN
          !------------------------------------------------------------------------------------------------------------
          ! Calculate the fraction transferred based on a fixed 
          !   threshold diameter DPCUT_IMTR. 
          !------------------------------------------------------------------------------------------------------------
          DGN_AKK_IMTR = 1.0D+06 * DPAKK * CONV_DPAM_TO_DGN(1)           ! [um]
          XNUM = XNUM_FACTOR * LOG( DPCUT_IMTR / DGN_AKK_IMTR )          ! [1]
          XNUM = MAX( XNUM, X3_TERM )                                    ! limit for stability as in BS2003
          X3 = XNUM - X3_TERM                                            ! [1]
          FNUM = 0.5D+00 * ERFC( XNUM )                                  ! number fraction transferred from AKK to ACC
          F3   = 0.5D+00 * ERFC( X3   )                                  ! mass   fraction transferred from AKK to ACC
          ! WRITE(34,'(9D12.4)')DGN_AKK_IMTR,DPCUT_IMTR,DPAKK*1.0D+06,AERO(NUMB_AKK_1),AERO(MASS_AKK_SULF),FNUM,F3
        ELSEIF( IMTR_METHOD .EQ. 3 ) THEN
          !------------------------------------------------------------------------------------------------------------
          ! Calculate the fraction transferred based on the 
          !   diameter of intersection of the AKK and ACC modes. 
          !------------------------------------------------------------------------------------------------------------
          DGN_AKK_IMTR = 1.0D+06 * DPAKK * CONV_DPAM_TO_DGN(1)           ! [um]
          DGN_ACC_IMTR = 1.0D+06 * DPACC * CONV_DPAM_TO_DGN(2)           ! [um]
          IF( AERO( NUMB_ACC_1 ) .GT. 1.0D+06 ) THEN                     ! mode ACC not essentially empty
            XNUM = GETXNUM(AERO(NUMB_AKK_1),AERO(NUMB_ACC_1),DGN_AKK_IMTR,DGN_ACC_IMTR,LNSG_AKK,LNSG_ACC)   ! [1]
          ELSE                                                           ! mode ACC essentially empty - use Method 2
            XNUM = XNUM_FACTOR * LOG( DPCUT_IMTR / DGN_AKK_IMTR )        ! [1]
          ENDIF
          XNUM = MAX( XNUM, X3_TERM )                                    ! limit for stability as in BS2003
          X3 = XNUM - X3_TERM                                            ! [1]
          FNUM = 0.5D+00 * ERFC( XNUM )                                  ! number fraction transferred from AKK to ACC
          F3   = 0.5D+00 * ERFC( X3   )                                  ! mass   fraction transferred from AKK to ACC
          ! WRITE(34,'(9D12.4)')DGN_AKK_IMTR,DGN_ACC_IMTR,AERO(NUMB_AKK_1),AERO(NUMB_ACC_1),FNUM,F3
        ENDIF
        DEL_MASS = AERO( MASS_AKK_SULF ) * F3                            ! mass   concentration transferred [ug/m^3]
        DEL_NUMB = AERO( NUMB_AKK_1    ) * FNUM                          ! number concentration transferred [# /m^3]
        AERO( MASS_AKK_SULF ) = AERO( MASS_AKK_SULF ) - DEL_MASS         ! update AKK mass   concentration  [ug/m^3]
        AERO( NUMB_AKK_1    ) = AERO( NUMB_AKK_1    ) - DEL_NUMB         ! update AKK number concentration  [# /m^3]
        AERO( MASS_ACC_SULF ) = AERO( MASS_ACC_SULF ) + DEL_MASS         ! update ACC mass   concentration  [ug/m^3]
        AERO( NUMB_ACC_1    ) = AERO( NUMB_ACC_1    ) + DEL_NUMB         ! update ACC number concentration  [# /m^3]
        IF( WRITE_LOG ) THEN
          WRITE(AUNIT1,'(/A5,9D11.3)') 'IMTR:',DPAKK,DPACC,DPAKK0,FNUM,F3,
     &                                         DEL_MASS,DEL_NUMB,AERO( NUMB_AKK_1 ), AERO( NUMB_ACC_1 )
        ENDIF
      ENDIF

      DIAGTMP1( 6,:) = ( AERO(:) - AEROTMP2(:) ) / TSTEP
      DIAGTMP1(14,:) = ( AERO(:) - AEROTMP2(:) ) / TSTEP
       
      !----------------------------------------------------------------------------------------------------------------
      ! Put all sea salt sulfate back into the accumulation mode (SSA) if the
      !   mechanism uses both the SSA and SSC modes.
      !----------------------------------------------------------------------------------------------------------------
      IF ( NUMBER_OF_SEASALT_MODES .EQ. 2 ) THEN
        AERO( MASS_SSA_SULF ) = AERO( MASS_SSA_SULF ) +  AERO( MASS_SSC_SULF )
        AERO( MASS_SSC_SULF ) = TINYNUMER
      ENDIF

      !----------------------------------------------------------------------------------------------------------------
      ! Get final total mass concentration for each model species [ug/m^3].
      ! Adjust final aerosol and gas-phase species by rescaling to enforce mass 
      !   conservation to machine precision.
      ! Precise mass conservation to machine precision is not necessarily
      !   conserved otherwise due to formulation of the model equations
      !   in terms of production and loss terms that may not precisely
      !   cancel on some occasions that they should due to their distinct
      !   treatment.
      !--------------------- -------------------------------------------------------------------------------------------
      ! WRITE(*,*)'TOT_SULF 1 = ',SUM( AERO( SULF_MAP(:)) )
      IF( MASS_ADJ ) THEN
        CALL SPCMASSES(AERO,GAS,SPCMASS2)
        CALL MASSADJ(AERO,GAS,SPCMASS1,SPCMASS2,EMIS_MASS,AQSO4RATE,TSTEP)
      ENDIF
      ! WRITE(*,*)'TOT_SULF 2 = ',SUM( AERO( SULF_MAP(:)) )

      !----------------------------------------------------------------------------------------------------------------
      ! Limit low mass or number concentrations.
      !----------------------------------------------------------------------------------------------------------------
      AERO(:) = MAX( AERO(:), MINCONC )

      !----------------------------------------------------------------------------------------------------------------
      ! Budget diagnostics.
      !----------------------------------------------------------------------------------------------------------------
      DIAGTMP1(7,:)          = ( AERO(:) - AEROTMP1(:) ) / TSTEP   ! Actual total differences for both mass and number.
      DIAG(:,:)              = 0.0D+00 
      DIAG(1: 7,NUMB_MAP(:)) = DIAGTMP1(1: 7,NUMB_MAP(:))  ! Save diagnostics 1-7  for the number concentrations only.  
      DIAG(8:14,:)           = DIAGTMP1(8:14,:)            ! Save diagnostics 8-15 for the mass   concentrations only.  
      DIAG(15,:)             = DIAGTMP1(7,:)               ! ... cont'd ...                      
      DIAG(8:15,NUMB_MAP(:)) = 0.0D+00                     ! Zero diagnostics 8-15 for the number concentrations.
      !----------------------------------------------------------------------------------------------------------------
      ! Rescale each term of the budget for species I to get the correct sum. 
      !----------------------------------------------------------------------------------------------------------------
c      DO I=1, NAEROBOX
c        AEROTMP2(I) = 0.0D+00
c        DO K=1, 6 
c         AEROTMP2(I) = AEROTMP2(I) + DIAG(K,I) + DIAG(K+7,I)
c        ENDDO 
c        AEROTMP2(I) = AEROTMP2(I) + DIAG(14,I)
c         DO J=1, 6
c          DIAG(J,I) = DIAG(J,I) * ( DIAGTMP1(7,I) / ( AEROTMP2(I) + TINYDENOM ) )
c        ENDDO 
c        DO J=8, NDIAG_AERO-1
c          DIAG(J,I) = DIAG(J,I) * ( DIAGTMP1(7,I) / ( AEROTMP2(I) + TINYDENOM ) )
c        ENDDO 
c      ENDDO 

c Diagnostics for PM1 PM(1), PM2.5 PM(2) and PM10 PM(3)
        IF ( ILAY.eq.1) THEN 
        PM(:) = 0.d0

      DO I=1, NWEIGHTS                               ! loop over modes (quadrature points) 
          DO Q=1, NM(I)                              ! loop over species defined for mode I 
           
               PM(1) = PM(1) + AERO(MASS_MAP(I,Q))*
     &          0.5d0*(erf(log(1.d0/DGN(I))/(sqrt(2.d0)*log(sig0(I)))) 
     &              -erf(log(0.001d0/DGN(I))/(sqrt(2.d0)*log(sig0(I)))))

               PM(2) = PM(2) + AERO(MASS_MAP(I,Q))*
     &          0.5d0*(erf(log(2.5d0/DGN(I))/(sqrt(2.d0)*log(sig0(I)))) 
     &              -erf(log(0.001d0/DGN(I))/(sqrt(2.d0)*log(sig0(I)))))

               PM(3) = PM(3) + AERO(MASS_MAP(I,Q))*
     &          0.5d0*(erf(log(10.d0/DGN(I))/(sqrt(2.d0)*log(sig0(I)))) 
     &              -erf(log(0.001d0/DGN(I))/(sqrt(2.d0)*log(sig0(I)))))

         ENDDO
         ENDDO
        ENDIF



90000 FORMAT(3I6,D15.5)
90001 FORMAT(1I6,2D15.5)
90002 FORMAT(1I6,D15.5)
90003 FORMAT(F9.1,F6.3,F7.2,7D13.3)
90004 FORMAT(I20,5D20.6)
90005 FORMAT(I6,3D30.6)
90006 FORMAT(I6,2D30.6)
90007 FORMAT(F7.2,F5.2,7D13.5,F8.5)
90008 FORMAT(2I5,4D15.5)
90009 FORMAT(I3,F11.7,2F13.3,D11.3,F8.3,A5,2F8.3)
90010 FORMAT(F7.2,F5.2,4D13.5,F16.8,I3)
      END SUBROUTINE MATRIX

