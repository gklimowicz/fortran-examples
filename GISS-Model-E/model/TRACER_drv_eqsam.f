#include "rundeck_opts.h"
      SUBROUTINE NITRATE_THERMO_DRV(LTOP)
!@sum
!@+     This routine sets up for and calls the thermodynamic module for aerosol
!@+     gas-particle partitioning.
!@+
!@+      A version of EQSAM (eqsam_v03d) is the current thermodynamic model. 
!@auth Susanne Bauer


!----------------------------------------------------------------------------------------------------------------------
!     This routine sets up for and calls the thermodynamic module for aerosol
!     gas-particle partitioning.
!
!     A version of EQSAM (eqsam_v03d) is the current thermodynamic model. 
!
!     EQSAM is called with control variable IOPT=1. 
!
!     Although EQSAM takes as input the total S(VI) (H2SO4+SO4=), since the
!     aerosol model does not necessarily transfer all H2SO4 to the aerosol
!     phase (depending on configuration), we pass only the particulate SO4
!     as the total sulfate to EQSAM.
!
!     Also, this version of EQSAM takes as input the mineral cation 
!     concentrations K+, Ca++, Mg++, Na+. Given the 'well-mixed' treatment
!     of inorganic aerosol constituents, these cations are included.
!----------------------------------------------------------------------------------------------------------------------
      USE TRACER_COM, only: ntm, trm
      use TRACER_COM, only: n_Clay, n_HNO3, n_NH3, n_NH4, n_NO3p
#ifdef TRACERS_AEROSOLS_SEASALT
      use TRACER_COM, only: n_seasalt1, n_seasalt2
#endif  /* TRACERS_AEROSOLS_SEASALT */
      use TRACER_COM, only: n_Silt1, n_Silt2, n_Silt3
      use TRACER_COM, only: n_SO4, n_SO4_d1, n_SO4_d2, n_SO4_d3
      USE AEROSOL_SOURCES, only: off_HNO3

      USE RESOLUTION, only : im,jm,lm     ! dimensions
      USE ATM_COM, only :   t            ! potential temperature (C)
     $                     ,q            ! saturatered pressure
      USE MODEL_COM, only : dtsrc
      USE GEOM, only: axyp,BYAXYP
      USE CONSTANT,   only: mair,gasc,lhe
      USE FLUXES, only: tr3Dsource
      USE ATM_COM,   only: pmid,pk,MA   ! midpoint pressure in hPa (mb)
!                                             and pk is t mess up factor
      USE DOMAIN_DECOMP_ATM,only: GRID, getDomainBounds
      use TRDIAG_COM, only: taijls=>taijls_loc,ijlt_aH2O,ijlt_apH

      IMPLICIT NONE

      INTEGER:: j,l,i,J_0, J_1,n,I_0,I_1,LTOP
      INTEGER, PARAMETER :: NCA  = 11    ! fixed number of input variables
      INTEGER, PARAMETER :: NCO  = 37    ! fixed number of output variables
      INTEGER, PARAMETER :: IOPT =  1    ! =1 selects the metastable (wet) state and history
!     INTEGER, PARAMETER :: IOPT =  2    ! =2 selects the solid      (dry) state and history
      INTEGER, PARAMETER :: LOOP =  1    ! only a single time step done
      INTEGER, PARAMETER :: IMAX =  1    ! only a single time step done
      INTEGER, PARAMETER :: AUNIT1 =  66

      ! Functions
      REAL*8 :: QSAT, AVOL

      ! Variables
    
      REAL(8) :: ASO4      ! aerosol sulfate       [ug/m^3]
      REAL(8) :: ANO3      ! aerosol nitrate       [ug/m^3]
      REAL(8) :: ANH4      ! aerosol ammonium      [ug/m^3]
      REAL(8) :: AH2O      ! aerosol water         [ug/m^3]
      REAL(8) :: GNH3      ! gas-phase ammonia     [ugNH4/m^3] as ammonium (MW)
      REAL(8) :: GHNO3     ! gas-phase nitric acid [ugNO3/m^3] as nitrate  (MW)
      REAL(8) :: DUST      ! fine dust(sol+insol) [ug/m^3]
      REAL(8) :: SALT      ! fine salt(sol+insol) [ug/m^3]
      REAL(8) :: TK        ! absolute temperature  [K]          
      REAL(8) :: RH        ! relative humidity     [0-1] w/r/t liquid water
      REAL(8) :: RHD       ! RH of deliquescence   [0-1]
      REAL(8) :: RHC       ! RH of crystallization [0-1]
       
      !------------------------------------------------------------------------------------------------------
      ! Input/Output to/from EQSAM
      !------------------------------------------------------------------------------------------------------
      REAL(4) :: YI(IMAX,NCA)            ! [umol/m^3] for chemical species - input
      REAL(4) :: YO(IMAX,NCO)            ! [umol/m^3] for chemical species - output
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
  

      REAL(4), PARAMETER :: MASS_FRAC_K  = 0.0028  ! From Ghan et al. (2001).
      REAL(4), PARAMETER :: MASS_FRAC_CA = 0.024   !   JGR, Vol. 106, p. 5295-5316.
      REAL(4), PARAMETER :: MASS_FRAC_MG = 0.0038  !   on p. 5296
      REAL(4), PARAMETER :: MASS_FRAC_NA = 0.014   !   "water sol. mass frac. in soil dust"

      REAL(4), PARAMETER :: FRAC_DUST  = 1.d0                             ! [1] fraction of dust conc. passed to EQSAM         
      REAL(4), PARAMETER :: FRAC_SALT  = 0.001                            ! [1] fraction of salt conc. passed to EQSAM         
      REAL(4), PARAMETER :: CONV_KION  = FRAC_DUST * MASS_FRAC_K  / MW_K  ! [mol/g]
      REAL(4), PARAMETER :: CONV_CAION = FRAC_DUST * MASS_FRAC_CA / MW_CA ! [mol/g]
      REAL(4), PARAMETER :: CONV_MGION = FRAC_DUST * MASS_FRAC_MG / MW_MG ! [mol/g]
      REAL(4), PARAMETER :: CONV_NAION = FRAC_DUST * MASS_FRAC_NA / MW_NA ! [mol/g]
      !------------------------------------------------------------------------------------------------------
      ! Fraction of sea salt (NaCl) mass that is Na, and is Cl.
      !------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: RAT_NA = MW_NA / ( MW_NA + MW_CL ) ! [1] 
      REAL(8), PARAMETER :: RAT_CL = MW_CL / ( MW_NA + MW_CL ) ! [1] 
      REAL(4), PARAMETER :: RMW_GNH3  = 1.d0 / MW_GNH3          ! [mol/g]
      REAL(4), PARAMETER :: RMW_ANH4  = 1.d0 / MW_ANH4          ! [mol/g]
      REAL(4), PARAMETER :: RMW_GHNO3 = 1.d0 / MW_GHNO3         ! [mol/g]
      REAL(4), PARAMETER :: RMW_ANO3  = 1.d0 / MW_ANO3          ! [mol/g]
      REAL(4), PARAMETER :: RMW_ASO4  = 1.d0 / MW_ASO4          ! [mol/g]
      REAL(4), PARAMETER :: RMW_NACL  = 1.d0 / MW_NACL          ! [mol/g]
      REAL(4), PARAMETER :: RMW_NA    = 1.d0 / MW_NA            ! [mol/g]
      REAL(8), PARAMETER :: RMW_CL    = 1.d0 / MW_CL            ! [mol/g]
      REAL(8), PARAMETER :: RHMAX  = 0.995D+00                 ! [0-1]
      REAL(8), PARAMETER :: RHMIN  = 0.010D+00                 ! [0-1]   
      REAL(8), PARAMETER :: SMALL_SO4 = 1.0D-05                ! [umol SO4/m^3] EQSAM has crashed at low RH and low sulfate conc.

      REAL(8) :: H   ! local RH, with RHMIN < H < RHMAX

      call getDomainBounds(grid, J_STRT =J_0, J_STOP =J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
#ifndef  TRACERS_SPECIAL_Shindell
      CALL READ_OFFHNO3(OFF_HNO3)
#endif
      YI(1,:) = 0.d0
      YO(1,:) = 0.d0

      DO L=1,LTOP                        
      DO J=J_0,J_1                               
      DO I=I_0,I_1
! meteo
      TK = pk(l,i,j)*t(i,j,l)           ! in [K]
      RH = q(i,j,l)/QSAT(pk(l,i,j)*t(i,j,l),lhe,pmid(l,i,j)) ! rH [0-1]
c avol [m3/gb] mass of air pro m3     
      AVOL = MA(l,i,j)*axyp(i,j)/mair*1000.d0*gasc*tk/(pmid(l,i,j)*100.d0) 
! gas and aerosol trm [kg/gb] -> [ug/m^3]
      GNH3 = trm(i,j,l,n_NH3)       *1.d9 /AVOL
      ANH4 = trm(i,j,l,n_NH4)       *1.d9 /AVOL
      ASO4 = trm(i,j,l,n_SO4)       *1.d9 /AVOL
      ANO3 = trm(i,j,l,n_NO3p)      *1.d9 /AVOL
      GHNO3= trm(i,j,l,n_HNO3)      *1.d9 /AVOL
      DUST = trm(i,j,l,n_Clay)      *1.d9 /AVOL
      SALT = trm(i,j,l,n_seasalt1)  *1.d9 /AVOL

      H = MAX( MIN( RH, RHMAX ), RHMIN )
      
      YI(1,1)  = TK                               ! [K]
      YI(1,2)  = H                                ! [0-1]
      YI(1,3)  = GNH3*RMW_GNH3   + ANH4*RMW_ANH4  ! from [ug/m^3] to [umol/m^3]
      YI(1,4)  =                   ASO4*RMW_ASO4  ! from [ug/m^3] to [umol/m^3]
      YI(1,5)  = GHNO3*RMW_GHNO3 + ANO3*RMW_ANO3  ! from [ug/m^3] to [umol/m^3]
      YI(1,6)  = RAT_NA*SALT*RMW_NA * FRAC_SALT   ! Sodium    from [ug dust/m^3] to [umol Na+/m^3]
      YI(1,7)  = RAT_CL*SALT*RMW_CL * FRAC_SALT   ! (HCl + Cl-)
      YI(1,8)  = DUST*CONV_KION                   ! Potassium from [ug dust/m^3] to [umol K+ /m^3]
      YI(1,9)  = DUST*CONV_CAION                  ! Calcium   from [ug dust/m^3] to [umol Ca+/m^3]
      YI(1,10) = DUST*CONV_MGION                  ! Magnesium from [ug dust/m^3] to [umol Mg+/m^3]
      YI(1,11) = pmid(l,i,j)                      ! [hPa]
      YI(1, :) = MAX( YI(1,:), 0.0E-10 )          ! Lower limit was 1.0E-10 before 102406.
      YI(1,4)  = YI(1,4) + SMALL_SO4              ! EQSAM has crashed at low RH and low sulfate conc.


      CALL EQSAM_V03D(YI,YO,NCA,NCO,IOPT,LOOP,IMAX,AUNIT1)


      GHNO3 = MAX(YO(1, 9) * MW_GHNO3, 0.d0)  ! from [umol/m^3] to [ug/m^3]
      GNH3  = MAX(YO(1,10) * MW_GNH3 , 0.d0)  ! from [umol/m^3] to [ug/m^3]
      AH2O  = MAX(YO(1,12)           , 0.d0)  ! already in [ugH2O/m^3]
      ANH4  = MAX(YO(1,19) * MW_ANH4 , 0.d0)  ! from [umol/m^3] to [ug/m^3]
      ANO3  = MAX(YO(1,20) * MW_ANO3 , 0.d0)  ! from [umol/m^3] to [ug/m^3]
      ASO4  = ( YO(1,21) - SMALL_SO4 ) * MW_ASO4  ! from [umol/m^3] to [ug/m^3]
      ASO4  = MAX( ASO4, 0.d0 )              ! 

      RHD   = 0.80D+00                            ! RHD = 0.80 for ammonium sulfate (Ghan et al., 2001).
      RHC   = 0.35D+00                            ! RHC = 0.35 for ammonium sulfate (Ghan et al., 2001).
 
! save aerosol water (ug/m3) and aerosol pH (dimensionless)
      taijls(I,J,L,ijlt_aH2O)=taijls(I,J,L,ijlt_aH2O)+AH2O
      taijls(I,J,L,ijlt_apH)=taijls(I,J,L,ijlt_apH)+(-log10(YO(1,37)+tiny(1.e0)))

! Nitrate production   from [ug/m^3] -> trm [kg/gb]
      tr3Dsource(i,j,l,1,n_NO3p)= ((ANO3 * 1.d-9 *AVOL) -trm(i,j,l,n_NO3p)) /dtsrc
! Ammonia residual
      tr3Dsource(i,j,l,1,n_NH3)= ((GNH3 * 1.d-9 *AVOL) -trm(i,j,l,n_NH3)) /dtsrc
! Ammonium production
      tr3Dsource(i,j,l,1,n_NH4)= ((ANH4 * 1.d-9 *AVOL) -trm(i,j,l,n_NH4)) /dtsrc
! Nitric Acid residual
      tr3Dsource(i,j,l,3,n_HNO3)= ((GHNO3 * 1.d-9 *AVOL) -trm(i,j,l,n_HNO3)) /dtsrc


      ENDDO
      ENDDO
      ENDDO

 ! Set source to zero above the chemistry:
      DO L=LTOP+1,LM
      DO J=J_0,J_1
      DO I=I_0,I_1
      tr3Dsource(i,j,l,1,n_NO3p)= 0.d0
      tr3Dsource(i,j,l,3,n_HNO3)= 0.d0
      tr3Dsource(i,j,l,1,n_NH3)= 0.d0
      tr3Dsource(i,j,l,1,n_NH4)= 0.d0
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE NITRATE_THERMO_DRV

   
