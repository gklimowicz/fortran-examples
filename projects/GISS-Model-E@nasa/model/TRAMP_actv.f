      MODULE AERO_ACTV
!      USE AMP_AEROSOL, ONLY: NACTV
      USE AERO_PARAM,  ONLY: NLAYS, AUNIT1
      USE AERO_CONFIG, ONLY: NMODES
!-------------------------------------------------------------------------------------------------------------------------
!@sum     The array NACTV(X,Y,Z,I) contains current values of the number of aerosol particles 
!@+       activated in clouds for each mode I for use outside of the MATRIX microphysical module.
!@+       Values in NACTV are saved in subr. MATRIX at each time step. 
!
!@auth    Susanne Bauer/Doug Wright
!
!-------------------------------------------------------------------------------------------------------------------------
      
      REAL(8), PARAMETER :: DENS_SULF = 1.77D+03    ! [kg/m^3] NH42SO4
      REAL(8), PARAMETER :: DENS_BCAR = 1.70D+03    ! [kg/m^3] Ghan et al. (2001) - MIRAGE
      REAL(8), PARAMETER :: DENS_OCAR = 1.00D+03    ! [kg/m^3] Ghan et al. (2001) - MIRAGE
      REAL(8), PARAMETER :: DENS_DUST = 2.60D+03    ! [kg/m^3] Ghan et al. (2001) - MIRAGE
      REAL(8), PARAMETER :: DENS_SEAS = 2.165D+03   ! [kg/m^3] NaCl, Ghan et al. (2001) used 1.90D+03

      CONTAINS


      SUBROUTINE GETACTFRAC(NMODEX,XNAP,XMAP5,RG,SIGMAG,TKELVIN,PTOT,WUPDRAFT,
     &                      AC,FRACACTN,FRACACTM,NACT,MACT)
!----------------------------------------------------------------------------------------------------------------------
!     12-12-06, DLW: Routine to set up the call to subr. ACTFRAC_MAT to calculate the 
!                    activated fraction of the number and mass concentrations, 
!                    as well as the number and mass concentrations activated 
!                    for each of NMODEX modes. The minimum dry radius for activation 
!                    for each mode is also returned. 
!
!     Each mode is assumed to potentially contains 5 chemical species:
!         (1) sulfate 
!         (2) BC 
!         (3) OC
!         (4) mineral dust
!         (5) sea salt 
!
!     The aerosol activation parameterizations are described in 
!
!         1. Abdul-Razzak et al.   1998, JGR, vol.103, p.6123-6131.
!         2. Abdul-Razzak and Ghan 2000, JGR, vol.105, p.6837-6844. 
!
!     and values for many of the required parameters were taken from 
!
!         3. Ghan et al. 2001, JGR vol 106, p.5295-5316.
!
!     With the density of sea salt set to the value used in ref. 3 (1900 kg/m^3), this routine 
!     yields values for the hygroscopicity parameters Bi in agreement with ref. 3. 
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, PARAMETER :: NCOMPS = 5

      ! Arguments.
      
      INTEGER :: NMODEX               ! number of modes [1]      
      REAL(8) :: XNAP(NMODEX)         ! number concentration for each mode [#/m^3]
      REAL(8) :: XMAP5(NMODEX,NCOMPS) ! mass concentration of each of the 5 species for each mode [ug/m^3]
      REAL(8) :: RG(NMODEX)           ! geometric mean dry radius for each mode [um]
      REAL(8) :: SIGMAG(NMODEX)       ! geometric standard deviation for each mode [um]
      REAL(8) :: TKELVIN              ! absolute temperature [K]
      REAL(8) :: PTOT                 ! ambient pressure [Pa]
      REAL(8) :: WUPDRAFT             ! updraft velocity [m/s]
      REAL(8) :: AC(NMODEX)           ! minimum dry radius for activation for each mode [um]
      REAL(8) :: FRACACTN(NMODEX)     ! activating fraction of number conc. for each mode [1]
      REAL(8) :: FRACACTM(NMODEX)     ! activating fraction of mass   conc. for each mode [1]
      REAL(8) :: NACT(NMODEX)         ! activating number concentration for each mode [#/m^3]
      REAL(8) :: MACT(NMODEX)         ! activating mass   concentration for each mode [ug/m^3]

      ! Local variables. 

      INTEGER :: I, J                 ! loop counters 
      REAL(8) :: XMAP(NMODEX)         ! total mass concentration for each mode [ug/m^3]
      REAL(8) :: BIBAR(NMODEX)        ! hygroscopicity parameter for each mode [1]

      ! Variables for mode-average hygroscopicity parameters. 
      
      REAL(8)       :: XR  (NMODEX,NCOMPS)  ! mass fraction for component J in mode I [1]   
      REAL(8), SAVE :: XNU (NCOMPS)         ! # of ions formed per formula unit solute for component J in mode I [1]
      REAL(8), SAVE :: XPHI(NCOMPS)         ! osmotic coefficient for component J in mode I [1]
      REAL(8), SAVE :: XMW (NCOMPS)         ! molecular weight for component J in mode I [kg/mol]
      REAL(8), SAVE :: XRHO(NCOMPS)         ! density of component J in mode I [kg/m^3]
      REAL(8), SAVE :: XEPS(NCOMPS)         ! soluble fraction of component J in mode I [1]

      REAL(8) :: SUMNUMER, SUMDENOM         ! scratch variables 

      REAL(8), PARAMETER :: NION_SULF = 3.00D+00    ! [1]
      REAL(8), PARAMETER :: NION_BCAR = 1.00D+00    ! [1]
      REAL(8), PARAMETER :: NION_OCAR = 1.00D+00    ! [1]
      REAL(8), PARAMETER :: NION_DUST = 2.30D+00    ! [1]
      REAL(8), PARAMETER :: NION_SEAS = 2.00D+00    ! [1] NaCl

      REAL(8), PARAMETER :: XPHI_SULF = 0.70D+00    ! [1]
      REAL(8), PARAMETER :: XPHI_BCAR = 1.00D+00    ! [1]
      REAL(8), PARAMETER :: XPHI_OCAR = 1.00D+00    ! [1]
      REAL(8), PARAMETER :: XPHI_DUST = 1.00D+00    ! [1]
      REAL(8), PARAMETER :: XPHI_SEAS = 1.00D+00    ! [1] NaCl

      REAL(8), PARAMETER :: MOLW_SULF = 132.0D-03   ! [kg/mol]
      REAL(8), PARAMETER :: MOLW_BCAR = 100.0D-03   ! [kg/mol]
      REAL(8), PARAMETER :: MOLW_OCAR = 100.0D-03   ! [kg/mol]
      REAL(8), PARAMETER :: MOLW_DUST = 100.0D-03   ! [kg/mol]
      REAL(8), PARAMETER :: MOLW_SEAS = 58.44D-03   ! [kg/m^3] NaCl

      REAL(8), PARAMETER :: XEPS_SULF = 1.00D+00    ! [1]
      REAL(8), PARAMETER :: XEPS_BCAR = 1.67D-06    ! [1]
      REAL(8), PARAMETER :: XEPS_OCAR = 0.78D+00    ! [1]
      REAL(8), PARAMETER :: XEPS_DUST = 0.13D+00    ! [1]
      REAL(8), PARAMETER :: XEPS_SEAS = 1.00D+00    ! [1] NaCl

      REAL(8), PARAMETER :: WMOLMASS = 18.01528D-03 ! molar mass of H2O     [kg/mol]
      REAL(8), PARAMETER :: DENH2O   =  1.00D+03    ! density of water [kg/m^3]

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      
      IF( FIRSTIME ) THEN
        FIRSTIME = .FALSE.
        XNU (1) = NION_SULF
        XNU (2) = NION_BCAR
        XNU (3) = NION_OCAR
        XNU (4) = NION_DUST
        XNU (5) = NION_SEAS
        XPHI(1) = XPHI_SULF
        XPHI(2) = XPHI_BCAR
        XPHI(3) = XPHI_OCAR
        XPHI(4) = XPHI_DUST
        XPHI(5) = XPHI_SEAS
        XMW (1) = MOLW_SULF
        XMW (2) = MOLW_BCAR
        XMW (3) = MOLW_OCAR
        XMW (4) = MOLW_DUST
        XMW (5) = MOLW_SEAS
        XRHO(1) = DENS_SULF
        XRHO(2) = DENS_BCAR
        XRHO(3) = DENS_OCAR
        XRHO(4) = DENS_DUST
        XRHO(5) = DENS_SEAS
        XEPS(1) = XEPS_SULF
        XEPS(2) = XEPS_BCAR
        XEPS(3) = XEPS_OCAR
        XEPS(4) = XEPS_DUST
        XEPS(5) = XEPS_SEAS
      ENDIF

      !--------------------------------------------------------------------------------------------------------------
      ! Calculate the mass fraction component J for each mode I. 
      !--------------------------------------------------------------------------------------------------------------
      DO I=1, NMODEX
        XMAP(I) = 0.0D+00
        DO J=1, NCOMPS
          XMAP(I) = XMAP(I) + XMAP5(I,J)
        ENDDO
        XR(I,:) = XMAP5(I,:) / MAX( XMAP(I), 1.0D-30 )   
        ! WRITE(*,'(I4,5F12.6)') I,XR(I,:)
      ENDDO

      !--------------------------------------------------------------------------------------------------------------
      ! Calculate the hygroscopicity parameter for each mode. 
      !--------------------------------------------------------------------------------------------------------------
      DO I=1, NMODEX
        SUMNUMER = 0.0D+00
        SUMDENOM = 0.0D+00
        DO J=1, NCOMPS
          SUMNUMER = SUMNUMER + XR(I,J)*XNU(J)*XPHI(J)*XEPS(J)/XMW(J)     ! [mol/kg] 
          SUMDENOM = SUMDENOM + XR(I,J)/XRHO(J)                           ! [m^3/kg] 
        ENDDO
        BIBAR(I) = ( WMOLMASS*SUMNUMER ) / ( DENH2O*SUMDENOM )            ! [1] 
      ENDDO
      ! WRITE(*,'(8D15.6)') BIBAR(:)

      !--------------------------------------------------------------------------------------------------------------
      ! Calculate the droplet activation parameters for each mode. 
      !--------------------------------------------------------------------------------------------------------------
      CALL ACTFRAC_MAT(NMODEX,XNAP,XMAP,RG,SIGMAG,BIBAR,TKELVIN,PTOT,WUPDRAFT,
     &                 AC,FRACACTN,FRACACTM,NACT,MACT)

      DO I=1, NMODEX
        IF(XNAP(I) .LT. 1.0D-06 ) FRACACTN(I) = 1.0D-30
      ENDDO

      RETURN 
      END SUBROUTINE GETACTFRAC


      SUBROUTINE ACTFRAC_MAT(NMODEX,XNAP,XMAP,RG,SIGMAG,BIBAR,TKELVIN,PTOT,WUPDRAFT,
     &                       AC,FRACACTN,FRACACTM,NACT,MACT)
!----------------------------------------------------------------------------------------------------------------------
!     12-12-06, DLW: Routine to calculate the activated fraction of the number 
!                    and mass concentrations, as well as the number and mass 
!                    concentrations activated for each of NMODEX modes. The 
!                    minimum dry radius for activation for each mode is also returned. 
!
!     The aerosol activation parameterizations are described in 
!
!         1. Abdul-Razzak et al.   1998, JGR, vol.103, p.6123-6131.
!         2. Abdul-Razzak and Ghan 2000, JGR, vol.105, p.6837-6844. 
! 
!     This routine is for the multiple-aerosol type parameterization. 
!----------------------------------------------------------------------------------------------------------------------
      USE DOMAIN_DECOMP_ATM,only: am_i_root
      IMPLICIT NONE

      ! Arguments.
      
      INTEGER :: NMODEX            ! number of modes [1]      
      REAL(8) :: XNAP(NMODEX)      ! number concentration for each mode [#/m^3]
      REAL(8) :: XMAP(NMODEX)      ! mass   concentration for each mode [ug/m^3]
      REAL(8) :: RG(NMODEX)        ! geometric mean radius for each mode [um]
      REAL(8) :: SIGMAG(NMODEX)    ! geometric standard deviation for each mode [um]
      REAL(8) :: BIBAR(NMODEX)     ! hygroscopicity parameter for each mode [1]
      REAL(8) :: TKELVIN           ! absolute temperature [K]
      REAL(8) :: PTOT              ! ambient pressure [Pa]
      REAL(8) :: WUPDRAFT          ! updraft velocity [m/s]
      REAL(8) :: AC(NMODEX)        ! minimum dry radius for activation for each mode [um]
      REAL(8) :: AC_2(NMODEX)        ! minimum dry radius for activation for each mode [um]
      REAL(8) :: AC_3(NMODEX)        ! minimum dry radius for activation for each mode [um]
      REAL(8) :: AC_5(NMODEX)        ! minimum dry radius for activation for each mode [um]
      REAL(8) :: FRACACTN(NMODEX)  ! activating fraction of number conc. for each mode [1]
      REAL(8) :: FRACACTN_2(NMODEX)  ! activating fraction of number conc. for each mode [1]
      REAL(8) :: FRACACTN_3(NMODEX)  ! activating fraction of number conc. for each mode [1]
      REAL(8) :: FRACACTN_5(NMODEX)  ! activating fraction of number conc. for each mode [1]
      REAL(8) :: FRACACTM(NMODEX)  ! activating fraction of mass   conc. for each mode [1]
      REAL(8) :: NACT(NMODEX)      ! activating number concentration for each mode [#/m^3]
      REAL(8) :: MACT(NMODEX)      ! activating mass   concentration for each mode [ug/m^3]

      ! Parameters.
      
      REAL(8), PARAMETER :: PI            = 3.141592653589793D+00
      REAL(8), PARAMETER :: TWOPI         = 2.0D+00 * PI
      REAL(8), PARAMETER :: SQRT2         = 1.414213562D+00
      REAL(8), PARAMETER :: THREESQRT2BY2 = 1.5D+00 * SQRT2

      REAL(8), PARAMETER :: AVGNUM   = 6.0221367D+23       ! [1/mol]
      REAL(8), PARAMETER :: RGASJMOL = 8.31451D+00         ! [J/mol/K]
      REAL(8), PARAMETER :: WMOLMASS = 18.01528D-03        ! molar mass of H2O     [kg/mol]
      REAL(8), PARAMETER :: AMOLMASS = 28.966D-03          ! molar mass of air     [kg/mol]
      REAL(8), PARAMETER :: ASMOLMSS = 132.1406D-03        ! molar mass of NH42SO4 [kg/mol]
      REAL(8), PARAMETER :: DENH2O   = 1.00D+03            ! density of water [kg/m^3]
      REAL(8), PARAMETER :: DENAMSUL = 1.77D+03            ! density of pure ammonium sulfate [kg/m^3]
      REAL(8), PARAMETER :: XNUAMSUL = 3.00D+00            ! # of ions formed when the salt is dissolved in water [1]
      REAL(8), PARAMETER :: PHIAMSUL = 1.000D+00           ! osmotic coefficient value in A-R 1998. [1] 
      REAL(8), PARAMETER :: GRAVITY  = 9.81D+00            ! grav. accel. at the Earth's surface [m/s/s] 
      REAL(8), PARAMETER :: HEATVAP  = 40.66D+03/WMOLMASS  ! latent heat of vap. for water and Tnbp [J/kg] 
      REAL(8), PARAMETER :: CPAIR    = 1006.0D+00          ! heat capacity of air [J/kg/K] 
      REAL(8), PARAMETER :: T0DIJ    = 273.15D+00          ! reference temp. for DV [K] 
      REAL(8), PARAMETER :: P0DIJ    = 101325.0D+00        ! reference pressure for DV [Pa] 
      REAL(8), PARAMETER :: DIJH2O0  = 0.211D-04           ! reference value of DV [m^2/s] (P&K,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------    
      ! REAL(8), PARAMETER :: T0DIJ    = 283.15D+00          ! reference temp. for DV [K] 
      ! REAL(8), PARAMETER :: P0DIJ    = 80000.0D+00         ! reference pressure for DV [Pa] 
      ! REAL(8), PARAMETER :: DIJH2O0  = 0.300D-04           ! reference value of DV [m^2/s] (P&K,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: DELTAV   = 1.096D-07           ! vapor jump length [m]  
      REAL(8), PARAMETER :: DELTAT   = 2.160D-07           ! thermal jump length [m]  
      REAL(8), PARAMETER :: ALPHAC   = 1.000D+00           ! condensation mass accommodation coefficient [1]  
      REAL(8), PARAMETER :: ALPHAT   = 0.960D+00           ! thermal accommodation coefficient [1]  

      ! Local variables. 

      INTEGER            :: I                              ! loop counter 
      REAL(8)            :: DV                             ! diffusion coefficient for water [m^2/s] 
      REAL(8)            :: DVPRIME                        ! modified diffusion coefficient for water [m^2/s] 
      REAL(8)            :: DUMW, DUMA                     ! scratch variables [s/m] 
      REAL(8)            :: WPE                            ! saturation vapor pressure of water [Pa]  
      REAL(8)            :: SURTEN                         ! surface tension of air-water interface [J/m^2] 
      REAL(8)            :: XKA                            ! thermal conductivity of air [J/m/s/K]  
      REAL(8)            :: XKAPRIME                       ! modified thermal conductivity of air [J/m/s/K]  
      REAL(8)            :: ETA(NMODEX)                    ! model parameter [1]  
      REAL(8)            :: ZETA                           ! model parameter [1]  
      REAL(8)            :: XLOGSIGM(NMODEX)               ! ln(sigmag) [1]   
      REAL(8)            :: A                              ! [m]
      REAL(8)            :: G                              ! [m^2/s]   
      REAL(8)            :: RDRP                           ! [m]   
      REAL(8)            :: F1                             ! [1]   
      REAL(8)            :: F2                             ! [1]
      REAL(8)            :: ALPHA                          ! [1/m]
      REAL(8)            :: GAMMA                          ! [m^3/kg]   
      REAL(8)            :: SM(NMODEX)                     ! [1]   
      REAL(8)            :: DUM                            ! [1/m]    
      REAL(8)            :: U                              ! argument to error function [1]
      REAL(8)            :: ERF                            ! error function [1], but not declared in an f90 module 
      REAL(8)            :: SMAX                           ! maximum supersaturation [1]
      
!----------------------------------------------------------------------------------------------------------------------
!     RDRP is the radius value used in Eqs.(17) & (18) and was adjusted to yield eta and zeta 
!     values close to those given in A-Z et al. 1998 Figure 5. 
!----------------------------------------------------------------------------------------------------------------------
      RDRP = 0.105D-06   ! [m] Tuned to approximate the results in Figures 1-5 in A-Z et al. 1998.  
!----------------------------------------------------------------------------------------------------------------------
!     These variables are common to all modes and need only be computed once. 
!----------------------------------------------------------------------------------------------------------------------
      DV = DIJH2O0*(P0DIJ/PTOT)*(TKELVIN/T0DIJ)**1.94D+00                 ! [m^2/s] (P&K,2nd ed., p.503)
      SURTEN = 76.10D-03 - 0.155D-03 * (TKELVIN-273.15D+00)               ! [J/m^2] 
      WPE = EXP( 77.34491296D+00 - 7235.424651D+00/TKELVIN - 8.2D+00*LOG(TKELVIN) + TKELVIN*5.7113D-03 )  ! [Pa] 
      DUMW = SQRT(TWOPI*WMOLMASS/RGASJMOL/TKELVIN)                        ! [s/m] 
      DVPRIME = DV / ( (RDRP/(RDRP+DELTAV)) + (DV*DUMW/(RDRP*ALPHAC)) )   ! [m^2/s] - Eq. (17) 
      XKA = (5.69D+00+0.017D+00*(TKELVIN-273.15D+00))*418.4D-05           ! [J/m/s/K] (0.0238 J/m/s/K at 273.15 K)
      DUMA = SQRT(TWOPI*AMOLMASS/RGASJMOL/TKELVIN)                        ! [s/m]
      XKAPRIME = XKA / ( ( RDRP/(RDRP+DELTAT) ) + ( XKA*DUMA/(RDRP*ALPHAT*DENH2O*CPAIR) ) )   ! [J/m/s/K]
      G = 1.0D+00 / ( (DENH2O*RGASJMOL*TKELVIN) / (WPE*DVPRIME*WMOLMASS)
     :                + ( (HEATVAP*DENH2O) / (XKAPRIME*TKELVIN) )
     :                * ( (HEATVAP*WMOLMASS) / (RGASJMOL*TKELVIN) - 1.0D+00 ) )               ! [m^2/s]
      A = (2.0D+00*SURTEN*WMOLMASS)/(DENH2O*RGASJMOL*TKELVIN)                                 ! [m] 
      ALPHA = (GRAVITY/(RGASJMOL*TKELVIN))*((WMOLMASS*HEATVAP)/(CPAIR*TKELVIN) - AMOLMASS)    ! [1/m] 
      GAMMA = (RGASJMOL*TKELVIN)/(WPE*WMOLMASS) 
     &      + (WMOLMASS*HEATVAP*HEATVAP)/(CPAIR*PTOT*AMOLMASS*TKELVIN)                        ! [m^3/kg]
      DUM = SQRT(ALPHA*WUPDRAFT/G)                  ! [1/m] 
      ZETA = 2.D+00*A*DUM/3.D+00                    ! [1] 
      !----------------------------------------------------------------------------------------------------------------
      ! WRITE(1,'(A27,4D15.5)')'SURTEN,WPE,A            =',SURTEN,WPE,A
      ! WRITE(1,'(A27,4D15.5)')'XKA,XKAPRIME,DV,DVPRIME =',XKA,XKAPRIME,DV,DVPRIME
      ! WRITE(1,'(A27,4D15.5)')'ALPHA,GAMMA,G, ZETA     =',ALPHA,GAMMA,G,ZETA
!----------------------------------------------------------------------------------------------------------------------
!     These variables must be computed for each mode. 
!----------------------------------------------------------------------------------------------------------------------
      XLOGSIGM(:) = LOG(SIGMAG(:))                                                    ! [1] 
      SMAX = 0.0D+00                                                                  ! [1]
      DO I=1, NMODEX
        SM(I) = ( 2.0D+00/SQRT(BIBAR(I)) ) * ( A/(3.0D-06*RG(I)) )**1.5D+00           ! [1] 
        ETA(I) = DUM**3 / (TWOPI*DENH2O*GAMMA*XNAP(I))                                ! [1] 
        !--------------------------------------------------------------------------------------------------------------
        ! WRITE(1,'(A27,I4,4D15.5)')'I,ETA(I),SM(I) =',I,ETA(I),SM(I)
        !--------------------------------------------------------------------------------------------------------------
        F1 = 0.5D+00 * EXP(2.50D+00 * XLOGSIGM(I)**2)                                 ! [1] 
        F2 = 1.0D+00 +     0.25D+00 * XLOGSIGM(I)                                     ! [1] 
        SMAX = SMAX + (   F1*(  ZETA  / ETA(I)              )**1.50D+00 
     &                  + F2*(SM(I)**2/(ETA(I)+3.0D+00*ZETA))**0.75D+00 ) / SM(I)**2  ! [1] - Eq. (6)
      ENDDO 
      SMAX = 1.0D+00 / SQRT(SMAX)                                                     ! [1]
      DO I=1, NMODEX
        AC(I)       = RG(I) * ( SM(I) / SMAX )**0.66666666666666667D+00               ! [um]
        U           = LOG(AC(I)/RG(I)) / ( SQRT2 * XLOGSIGM(I) )                      ! [1]
        FRACACTN(I) = 0.5D+00 * (1.0D+00 - ERF(U))                                    ! [1]
        FRACACTM(I) = 0.5D+00 * (1.0D+00 - ERF(U - THREESQRT2BY2*XLOGSIGM(I) ) )      ! [1]
        NACT(I)     = FRACACTN(I) * XNAP(I)                                           ! [#/m^3]
        MACT(I)     = FRACACTM(I) * XMAP(I)                                           ! [ug/m^3]
        !--------------------------------------------------------------------------------------------------------------
      ENDDO 

      RETURN
      END SUBROUTINE ACTFRAC_MAT


      SUBROUTINE GCF(GAMMCF,A,X,GLN)

      IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!-----------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: ITMAX=10000
      REAL(8), PARAMETER :: EPS=3.0D-07
      REAL(8), PARAMETER :: FPMIN=1.0D-30
      REAL(8) :: A,GAMMCF,GLN,X
      INTEGER :: I
      REAL(8) :: AN,B,C,D,DEL,H
      GLN=GAMMLN(A)
      B=X+1.0D+00-A
      C=1.0D+00/FPMIN
      D=1.0D+00/B
      H=D
      DO I=1,ITMAX
        AN=-I*(I-A)
        B=B+2.0D+00
        D=AN*D+B
        IF(ABS(D).LT.FPMIN)D=FPMIN
        C=B+AN/C
        IF(ABS(C).LT.FPMIN)C=FPMIN
        D=1.0D+00/D
        DEL=D*C
        H=H*DEL
        IF(ABS(DEL-1.0D+00).LT.EPS)GOTO 1
      ENDDO
      WRITE(*,*)'AERO_ACTV: SUBROUTINE GCF: A TOO LARGE, ITMAX TOO SMALL', GAMMCF,A,X,GLN
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*H
      RETURN
      END SUBROUTINE GCF


      SUBROUTINE GSER(GAMSER,A,X,GLN)

      IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!-----------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: ITMAX=10000  ! was ITMAX=100   in Press et al. 
      REAL(8), PARAMETER :: EPS=3.0D-09  ! was EPS=3.0D-07 in Press et al.
      REAL(8) :: A,GAMSER,GLN,X
      INTEGER :: N
      REAL(8) :: AP,DEL,SUM
      GLN=GAMMLN(A)
      IF(X.LE.0.D+00)THEN
        IF(X.LT.0.)STOP 'AERO_ACTV: SUBROUTINE GSER: X < 0 IN GSER'
        GAMSER=0.D+00
        RETURN
      ENDIF
      AP=A
      SUM=1.D+00/A
      DEL=SUM
      DO N=1,ITMAX
        AP=AP+1.D+00
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GOTO 1
      ENDDO
      WRITE(*,*)'AERO_ACTV: SUBROUTINE GSER: A TOO LARGE, ITMAX TOO SMALL'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END SUBROUTINE GSER


      DOUBLE PRECISION FUNCTION GAMMLN(XX)

      IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!-----------------------------------------------------------------------------------------------------------------------
      REAL(8) :: XX
      INTEGER J
      DOUBLE PRECISION SER,STP,TMP,X,Y,COF(6)
      SAVE COF,STP
      DATA COF,STP/76.18009172947146D0,-86.50532032941677D0,
     &24.01409824083091D0,-1.231739572450155D0,.1208650973866179D-2,
     &-.5395239384953D-5,2.5066282746310005D0/
      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO J=1,6
        Y=Y+1.D0
        SER=SER+COF(J)/Y
      ENDDO
      GAMMLN=TMP+LOG(STP*SER/X)
      RETURN
      END FUNCTION GAMMLN 


      DOUBLE PRECISION FUNCTION ERF(X)
      IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!-----------------------------------------------------------------------------------------------------------------------
      REAL(8) :: X
!U    USES GAMMP
      ERF = 0.d0
      IF(X.LT.0.0D+00)THEN
        ERF=-GAMMP(0.5D0,X**2)
      ELSE
        ERF= GAMMP(0.5D0,X**2)
      ENDIF
      RETURN
      END FUNCTION ERF


      DOUBLE PRECISION FUNCTION GAMMP(A,X)
      IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!-----------------------------------------------------------------------------------------------------------------------
      REAL(8) :: A,X
      REAL(8) :: GAMMCF,GAMSER,GLN
      IF(X.LT.0.0D+00.OR.A.LE.0.0D+00)THEN
        WRITE(*,*)'AERO_ACTV: FUNCTION GAMMP: BAD ARGUMENTS'
      ENDIF
      IF(X.LT.A+1.0D+00)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.0D+00-GAMMCF
      ENDIF
      RETURN
      END FUNCTION GAMMP


      END MODULE AERO_ACTV
      
