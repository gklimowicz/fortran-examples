#include "rundeck_opts.h"

      MODULE TRCHEM_Shindell_COM
!@sum  TRCHEM_Shindell_COM declares variables for tracer chemistry
!@+    and sources.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
c
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only  : dtsrc,Itime,ItimeI
      USE CONSTANT, only   : pi,mair,mwat,radian,byavog,avog,gasc
      USE ATM_COM, only    : MA, byMA, PMID, PK
      USE TRACER_COM, only : trm, ntm_chem
      use OldTracer_mod, only: TR_MM

      IMPLICIT NONE
      SAVE

      type rrmono_index
! monomolecular break ups
        integer :: HO2NO2_M__HO2_NO2=0
        integer :: N2O5_M__NO3_NO2=0
        integer :: Cl2O2_M__ClO_ClO=0
      end type rrmono_index

      type rrbi_index
! bimolecular reactions
        integer :: OH_HO2__H2O_O2=0
        integer :: OH_O3__HO2_O2=0
        integer :: OH_OH__H2O_O=0
        integer :: HO2_O3__OH_O2=0
        integer :: O3_NO__NO2_O2=0
        integer :: HO2_NO__OH_NO2=0
        integer :: NO2_O3__NO3_O2=0
        integer :: O1D_O2__O_O2=0
        integer :: O1D_M__O_M=0
        integer :: O1D_H2O__OH_OH=0
        integer :: O1D_CH4__OH_CH3O2=0
        integer :: CH4_OH__H2O_CH3O2=0
        integer :: CO_OH__HO2_O2=0
        integer :: OH_H2O2__H2O_HO2=0
        integer :: HO2_HO2__H2O2_O2=0
        integer :: OH_HNO3__H2O_NO3=0
        integer :: NO3_NO__NO2_NO2=0
        integer :: OH_HO2NO2__H2O_NO2=0
        integer :: H2_OH__HO2_H2O=0
        integer :: CH3O2_NO__HCHO_NO2=0
        integer :: HCHO_OH__HO2_CO=0
        integer :: CH3O2_HO2__CH3OOH_O2=0
        integer :: CH3OOH_OH__CH3O2_H2O=0
        integer :: NO2_NO3__NO_NO2=0
        integer :: NO3_NO3__NO2_NO2=0
        integer :: O_NO2__NO_O2=0
        integer :: CH3O2_CH3O2__HCHO_HCHO=0
        integer :: NO3_HCHO__HNO3_CO=0
        integer :: PAN_M__C2O3_NO2=0
        integer :: Isoprene_OH__HCHO_Alkenes=0
        integer :: Isoprene_O3__HCHO_Alkenes=0
        integer :: Isoprene_NO3__HO2_Alkenes=0
        integer :: AlkylNit_OH__NO2_XO2=0
        integer :: Alkenes_OH__HCHO_HO2=0
        integer :: Alkenes_O3__HCHO_CO=0
        integer :: Alkenes_NO3__HCHO_NO2=0
        integer :: Paraffin_OH__HO2_M=0
        integer :: Paraffin_RXPAR__M_M=0
        integer :: Aldehyde_OH__C2O3_M=0
        integer :: C2O3_NO__HCHO_NO2=0
        integer :: C2O3_C2O3__HCHO_HCHO=0
        integer :: C2O3_HO2__HCHO_HO2=0
        integer :: ROR_M__Aldehyde_HO2=0
        integer :: ROR_M__HO2_M=0
        integer :: XO2_NO__NO2_M=0
        integer :: XO2_HO2__CH3OOH_O2=0
        integer :: XO2_XO2__M_M=0
        integer :: XO2N_NO__AlkylNit_M=0
        integer :: XO2N_HO2__CH3OOH_O2=0
        integer :: Cl_O3__ClO_O2=0
        integer :: ClO_O__Cl_O2=0
        integer :: Cl_OClO__ClO_ClO=0
        integer :: ClO_O3__Cl_O2=0
        integer :: ClO_O3__OClO_O2=0
        integer :: O_OClO__ClO_O2=0
        integer :: OH_Cl2__HOCl_Cl=0
        integer :: OH_HCl__H2O_Cl=0
        integer :: OH_HOCl__H2O_ClO=0
        integer :: O_HCl__OH_Cl=0
        integer :: O_HOCl__OH_ClO=0
        integer :: OClO_OH__HOCl_O2=0
        integer :: Cl_HOCl__Cl2_OH=0
        integer :: Cl_H2O2__HCl_HO2=0
        integer :: Cl_HO2__HCl_O2=0
        integer :: Cl_HO2__OH_ClO=0
        integer :: ClO_OH__HO2_Cl=0
        integer :: ClO_OH__HCl_O2=0
        integer :: ClO_HO2__HOCl_O2=0
        integer :: ClO_NO__NO2_Cl=0
        integer :: ClONO2_O__ClO_NO3=0
        integer :: HCl_O1D__Cl_OH=0
        integer :: NO_OClO__NO2_ClO=0
        integer :: HBr_OH__H2O_Br=0
        integer :: BrO_O__Br_O2=0
        integer :: Br_O3__BrO_O2=0
        integer :: BrO_NO__Br_NO2=0
        integer :: Br_HO2__HBr_O2=0
        integer :: BrO_HO2__HOBr_O2=0
        integer :: Br_OClO__BrO_ClO=0
        integer :: BrO_ClO__OClO_Br=0
        integer :: BrO_ClO__Br_Cl=0
        integer :: O1D_CH4__HCHO_H2=0
        integer :: BrO_BrO__Br_Br=0
        integer :: Br_H2O2__HBr_HO2=0
        integer :: BrO_OH__Br_HO2=0
        integer :: BrO_OH__HBr_O2=0
        integer :: Cl_CH4__HCl_CH3O2=0
        integer :: Cl_H2__HCl_HO2=0
        integer :: O_HBr__OH_Br=0
        integer :: ClO_CH3O2__Cl_HCHO=0
        integer :: N2O_O1D__N2_O2=0
        integer :: N2O_O1D__NO_NO=0
        integer :: O_O3__O2_O2=0
        integer :: O_OH__O2_H=0
        integer :: O_HO2__OH_O2=0
        integer :: HO2_NO__HNO3_M=0
#ifdef TRACERS_TERP
        integer :: Terpenes_OH__HCHO_Alkenes=0
        integer :: Terpenes_O3__HCHO_Alkenes=0
        integer :: Terpenes_NO3__HO2_Alkenes=0
#endif  /* TRACERS_TERP */
#ifdef TRACERS_dCO
        integer :: O1D_CH4__OH_dCH317O2=0
        integer :: CH4_OH__H2O_dCH317O2=0
        integer :: dC17O_OH__HO2_O2=0
        integer :: dCH317O2_NO__dHCH17O_NO2=0
        integer :: dHCH17O_OH__HO2_dC17O=0
        integer :: dCH317O2_HO2__dMe17OOH_O2=0
        integer :: dMe17OOH_OH__dCH317O2_H2O=0
        integer :: dCH317O2_CH3O2__dHCH17O_HCHO=0
        integer :: CH3O2_dCH317O2__HCHO_dHCH17O=0
        integer :: NO3_dHCH17O__HNO3_dC17O=0
        integer :: d17OPAN_M__dC217O3_NO2=0
        integer :: Isoprene_OH__dHCH17O_Alkenes=0
        integer :: Isoprene_O3__dHCH17O_Alkenes=0
        integer :: Alkenes_OH__dHCH17O_HO2=0
        integer :: Alkenes_O3__dHCH17O_CO=0
        integer :: Alkenes_O3__HCHO_dC17O=0
        integer :: Alkenes_NO3__dHCH17O_NO2=0
        integer :: d17Oald_OH__dC217O3_M=0
        integer :: dC217O3_NO__dHCH17O_NO2=0
        integer :: dC217O3_dC217O3__dHCH17O_dHCH17O=0
        integer :: dC217O3_HO2__dHCH17O_HO2=0
        integer :: d17OROR_M__d17Oald_HO2=0
        integer :: d17OROR_M__HO2_M=0
        integer :: O1D_CH4__dHCH17O_H2=0
        integer :: Cl_CH4__HCl_dCH317O2=0
        integer :: ClO_dCH317O2__Cl_dHCH17O=0
        integer :: Terpenes_OH__dHCH17O_Alkenes=0
        integer :: Terpenes_O3__dHCH17O_Alkenes=0

        integer :: O1D_CH4__OH_dCH318O2=0
        integer :: CH4_OH__H2O_dCH318O2=0
        integer :: dC18O_OH__HO2_O2=0
        integer :: dCH318O2_NO__dHCH18O_NO2=0
        integer :: dHCH18O_OH__HO2_dC18O=0
        integer :: dCH318O2_HO2__dMe18OOH_O2=0
        integer :: dMe18OOH_OH__dCH318O2_H2O=0
        integer :: dCH318O2_CH3O2__dHCH18O_HCHO=0
        integer :: CH3O2_dCH318O2__HCHO_dHCH18O=0
        integer :: NO3_dHCH18O__HNO3_dC18O=0
        integer :: d18OPAN_M__dC218O3_NO2=0
        integer :: Isoprene_OH__dHCH18O_Alkenes=0
        integer :: Isoprene_O3__dHCH18O_Alkenes=0
        integer :: Alkenes_OH__dHCH18O_HO2=0
        integer :: Alkenes_O3__dHCH18O_CO=0
        integer :: Alkenes_O3__HCHO_dC18O=0
        integer :: Alkenes_NO3__dHCH18O_NO2=0
        integer :: d18Oald_OH__dC218O3_M=0
        integer :: dC218O3_NO__dHCH18O_NO2=0
        integer :: dC218O3_dC218O3__dHCH18O_dHCH18O=0
        integer :: dC218O3_HO2__dHCH18O_HO2=0
        integer :: d18OROR_M__d18Oald_HO2=0
        integer :: d18OROR_M__HO2_M=0
        integer :: O1D_CH4__dHCH18O_H2=0
        integer :: Cl_CH4__HCl_dCH318O2=0
        integer :: ClO_dCH318O2__Cl_dHCH18O=0
        integer :: Terpenes_OH__dHCH18O_Alkenes=0
        integer :: Terpenes_O3__dHCH18O_Alkenes=0

        integer :: O1D_CH4__OH_d13CH3O2=0
        integer :: CH4_OH__H2O_d13CH3O2=0
        integer :: d13CO_OH__HO2_O2=0
        integer :: d13CH3O2_NO__dH13CHO_NO2=0
        integer :: dH13CHO_OH__HO2_d13CO=0
        integer :: d13CH3O2_HO2__d13MeOOH_O2=0
        integer :: d13MeOOH_OH__d13CH3O2_H2O=0
        integer :: d13CH3O2_CH3O2__dH13CHO_HCHO=0
        integer :: CH3O2_d13CH3O2__HCHO_dH13CHO=0
        integer :: NO3_dH13CHO__HNO3_d13CO=0
        integer :: d13CPAN_M__d13C2O3_NO2=0
        integer :: Isoprene_OH__dH13CHO_d13Calke=0
        integer :: Isoprene_O3__dH13CHO_d13Calke=0
        integer :: Isoprene_NO3__HO2_d13Calke=0
        integer :: d13Calke_OH__dH13CHO_HO2=0
        integer :: d13Calke_O3__dH13CHO_d13CO=0
        integer :: d13Calke_NO3__dH13CHO_NO2=0
        integer :: d13CPAR_OH__HO2_M=0
        integer :: d13CPAR_d13CXPAR__M_M=0
        integer :: d13Cald_OH__d13C2O3_M=0
        integer :: d13C2O3_NO__dH13CHO_NO2=0
        integer :: d13C2O3_d13C2O3__dH13CHO_dH13CHO=0
        integer :: d13C2O3_HO2__dH13CHO_HO2=0
        integer :: d13CROR_M__d13Cald_HO2=0
        integer :: d13CROR_M__HO2_M=0
        integer :: O1D_CH4__dH13CHO_H2=0
        integer :: Cl_CH4__HCl_d13CH3O2=0
        integer :: ClO_d13CH3O2__Cl_dH13CHO=0
        integer :: Terpenes_OH__dH13CHO_d13Calke=0
        integer :: Terpenes_O3__dH13CHO_d13Calke=0
        integer :: Terpenes_NO3__HO2_d13Calke=0
#endif  /* TRACERS_dCO */
      end type rrbi_index

      type rrtri_index
! trimolecular reactions
        integer :: O_O2__O3_M=0
        integer :: NO_O__NO2_M=0
        integer :: OH_OH__H2O2_M=0
        integer :: OH_NO2__HNO3_M=0
        integer :: HO2_NO2__HO2NO2_M=0
        integer :: NO3_NO2__N2O5_M=0
        integer :: OH_NO__HONO_M=0
        integer :: C2O3_NO2__PAN_M=0
        integer :: ClO_ClO__Cl2O2_M=0
        integer :: ClO_NO2__ClONO2_M=0
        integer :: BrO_NO2__BrONO2_M=0
#ifdef TRACERS_dCO
        integer :: dC217O3_NO2__d17OPAN_M=0
        integer :: dC218O3_NO2__d18OPAN_M=0
        integer :: d13C2O3_NO2__d13CPAN_M=0
#endif  /* TRACERS_dCO */
      end type rrtri_index

      type rrhet_index
! heterogeneous reactions
        integer :: N2O5_H2O__HNO3_HNO3=0
        integer :: ClONO2_H2O__HOCl_HNO3=0
        integer :: ClONO2_HCl__Cl_HNO3=0
        integer :: HOCl_HCl__Cl_H2O=0
        integer :: N2O5_HCl__Cl_HNO3=0
      end type rrhet_index

      type(rrmono_index) :: rrmono
      type(rrbi_index) :: rrbi
      type(rrtri_index) :: rrtri
      type(rrhet_index) :: rrhet

C**************  P  A  R  A  M  E  T  E  R  S  *******************
!@param p_1 number of reactants or products per reaction
!@param n_rx maximum number of chemical reactions in JPLRX
!@param n_bi maximum number of bimolecular reactions in JPLRX
!@param n_tri maximum number of trimolecular reactions in JPLRX
!@param n_nst maximum number of monomolecular decompositions in JPLRX
!@param n_het maximum number of heterogeneous reactions in JPLRX
!@param numfam number of chemical families in JPLRX
!@param n_rj number of photolysis reactions in JPLPH
!@param luselb Use reflective photolysis boundary treatment
!@param zlbatm Optical depth above which to set lower boundary
!@param CMEQ1 ?
!@param nc total number of molecules included (incl. O2 and N2)
!@param ny number of chemically calculated gases (no O2 or N2)
!@param O3MULT =2.14d-2 This is the conversion from (atm*cm) units
!@+     (i.e. 1000 Dobson Units) to KG/m2. It is: 
!@+     1.E4*2.69E19*48./6.02E26 where 1.E4 is cm2/m2, 2.69E19 is 
!@+     molecules/cm3 at 1 atm pressure, 48. is molecular wt of O3,
!@+     and 6.02E26 is Avogadro's number in molecules/Kmol.
!@param cpd conversion from molecules/cm3 to mole/m3
!@param BYO3MULT = 1/O3MULT
!@param pfix_H2 fixed ratio of H2/M
!@param pfix_Aldehyde fixed ratio of Aldehyde/M for initial conditions
!@param MWabyMWw ratio of molecular weights of air/water
!@param RKBYPIM=8.*RBOLTZ/pi/MASSN2O55=8.*1.38062D-23/3.14159/1.793D-25
!@param cboltz Boltzman's Constant = 1.3806d-19
!@param byradian 1/radian = conversion from radians to degrees
!@param LCOalt number of levels in the several tracer IC arrays
!@param LCH4alt number of levels in the CH4altIN array
!@param PCOalt pressures at LCOalt levels
!@param PCH4alt pressures at LCH4alt levels
!@param T_thresh threshold temperature used in master chem
!@param n2o_pppv default N2O L=1 overwriting in pppv
!@param cfc_pppv default CFC L=1 overwriting in pppv
!@param cfc_rad95 the average L=1 radiation code CFC11+CFC12 value
!@+     for 1995 (pppv). 
!@param fact_cfc ratio of our default CFC L=1 overwriting to the 
!@+     radiation's 1995 L=1 CFC11+CFC12 value.
!@param PSClatS SH latitude limit for PSCs
!@param PSClatN NH latitude limit for PSCs
!@param minKG minimum kg for trm before we set to this after change
!@param kbolt Boltzmann's constant in erg K-1 or (1d-7 J) K-1
!@param boltAvog8byPi kbolt*avog*8/pi derived quantity used in chemistry
      INTEGER, PARAMETER ::
     & LCOalt =   23,
     & LCH4alt=    6,
#ifdef TRACERS_TERP
     & n_bi_terp = 3, ! number of terpenes bimolecular reactions
#else
     & n_bi_terp = 0,
#endif  /* TRACERS_TERP */
     & ntm_shindell_nontransp = 26, ! number of non-transported Shindell tracers
#ifdef TRACERS_dCO
     & ntm_dCO_nontransp = 13, ! number of non-transported dCO tracers
     & n_bi_dCO = 87, ! number of dCO bimolecular reactions
     & n_tri_dCO = 3, ! number of dCO trimolecular reactions
     & n_rj_dCO = 21, ! number of dCO photochemical reactions
#else
     & ntm_dCO_nontransp = 0,
     & n_bi_dCO = 0,
     & n_tri_dCO = 0,
     & n_rj_dCO = 0,
#endif  /* TRACERS_dCO */
     & n_bi  =    96+n_bi_terp+n_bi_dCO,
     & n_nst =     3,
     & n_tri =    11+n_tri_dCO,
     & n_het =     5,
     & n_rx  = n_bi+n_nst+n_tri+n_het,
     & ntm_chem_nontransp=ntm_shindell_nontransp+ntm_dCO_nontransp,
     & ntm_chem_extra=2,
     & ny     = ntm_chem+ntm_chem_nontransp,
     & nc     = ny+ntm_chem_extra,
     & numfam =    4,
#ifdef TRACERS_dCO
! define below ntm_dCO_nontransp tracers
     & ndC217O3=   1+ntm_chem,
     & ndC218O3=   2+ntm_chem,
     & nd13C2O3=   3+ntm_chem,
     & nd13CXPAR=  4+ntm_chem,
     & nd17OROR =  5+ntm_chem,
     & nd18OROR =  6+ntm_chem,
     & nd13CROR =  7+ntm_chem,
     & nd17Oald =  8+ntm_chem,
     & nd18Oald =  9+ntm_chem,
     & nd13Cald = 10+ntm_chem,
     & ndCH317O2= 11+ntm_chem,
     & ndCH318O2= 12+ntm_chem,
     & nd13CH3O2= 13+ntm_chem,
#endif  /* TRACERS_dCO */
! define below ntm_shindell_nontransp tracers
     & nC2O3=      1+ntm_chem+ntm_dCO_nontransp,
     & nXO2=       2+ntm_chem+ntm_dCO_nontransp,
     & nXO2N=      3+ntm_chem+ntm_dCO_nontransp,
     & nRXPAR=     4+ntm_chem+ntm_dCO_nontransp,
     & nROR=       5+ntm_chem+ntm_dCO_nontransp,
     & nAldehyde=  6+ntm_chem+ntm_dCO_nontransp,
     & nH2O=       7+ntm_chem+ntm_dCO_nontransp,
     & nCH3O2=     8+ntm_chem+ntm_dCO_nontransp,
     & nH2=        9+ntm_chem+ntm_dCO_nontransp,
     & nOH=       10+ntm_chem+ntm_dCO_nontransp,
     & nHO2=      11+ntm_chem+ntm_dCO_nontransp,
     & nO3=       12+ntm_chem+ntm_dCO_nontransp,
     & nO=        13+ntm_chem+ntm_dCO_nontransp,
     & nO1D=      14+ntm_chem+ntm_dCO_nontransp,
     & nNO=       15+ntm_chem+ntm_dCO_nontransp,
     & nNO2=      16+ntm_chem+ntm_dCO_nontransp,
     & nNO3=      17+ntm_chem+ntm_dCO_nontransp,
     & nHONO=     18+ntm_chem+ntm_dCO_nontransp,
     & nCl2O2=    19+ntm_chem+ntm_dCO_nontransp,
     & nClO=      20+ntm_chem+ntm_dCO_nontransp,
     & nOClO=     21+ntm_chem+ntm_dCO_nontransp,
     & nCl2=      22+ntm_chem+ntm_dCO_nontransp,
     & nCl=       23+ntm_chem+ntm_dCO_nontransp,
     & nBrCl=     24+ntm_chem+ntm_dCO_nontransp,
     & nBrO=      25+ntm_chem+ntm_dCO_nontransp,
     & nBr=       26+ntm_chem+ntm_dCO_nontransp,
! define below ntm_chem_extra tracers
     & nO2=        1+ntm_chem+ntm_chem_nontransp,
     & nM=         2+ntm_chem+ntm_chem_nontransp, !you must always put nM last (highest number)
     & n_rj  =    28+n_rj_dCO,
     & p_1   =     2
C ----------------------------------------------     
c     & n_Ox=        1,    ! note, these
c     & n_NOx=       2,    ! first 15 species are
c     & n_N2O5=      3,    ! tracers, and therefore
c     & n_HNO3=      4,    ! these parameters are
c     & n_H2O2=      5,    ! to be defined in 
c     & n_CH3OOH=    6,    ! TRACER_COM.f.
c     & n_HCHO=      7,    ! Note the UNDERSCORE!
c     & n_HO2NO2=    8,    !  T
c     & n_CO=        9,    !  R
c     & n_CH4=      10,    !  A
c     & n_PAN=      11,    !  C
c     & n_Isoprene= 12,    !  E
c     & n_AlkylNit= 13,    !  R
c     & n_Alkenes=  14,    !  S
c     & n_Paraffin= 15,    !
c     & n_Terpenes= 16,    ! ---------------
C ----------------------------------------------   
     
      REAL*8, PARAMETER ::  O3MULT       = 2.14d-2,
     &                      BYO3MULT     = 1./O3MULT,
     &                      T_thresh     = 200.d0,
     &                      pfix_H2      = 560.d-9,
     &                      pfix_Aldehyde= 2.d-9,
     &                      MWabyMWw     = mair/mwat,
     &                      RKBYPIM      = 1.961d2,
     &                      cboltz       = 1.3806d-19,
     &                      zlbatm       = 4.d0,
     &                      CMEQ1        = 0.25d0,
     &                      byradian     = 1.d0/radian,
     &                      cpd          = 1.d6*byavog,
     &                      minKG        = 0.d0,
     &                      cfc_pppv     = 1722.d-12,
     &                      n2o_pppv     = 316.3d-9,
     &                      cfc_rad95    = 794.d-12,
     &                      fact_cfc     = cfc_pppv/cfc_rad95,
     &                      kbolt        = 1.d7*gasc/avog,
     &                      boltAvog8byPi= kbolt*avog*8.d0/pi

C Please note: since PCOalt is essentially the nominal 
C pressures for the 23-level GCM, I'm going to use it
C to define BrOx,ClOx,ClONOs,HCL,COIC,OxIC,CFCIC,N2OICX,CH4ICX too:
      REAL*8, PARAMETER, DIMENSION(LCOalt) :: PCOalt = (/
     & 0.9720D+03,0.9445D+03,0.9065D+03,
     & 0.8515D+03,0.7645D+03,0.6400D+03,0.4975D+03,0.3695D+03,
     & 0.2795D+03,0.2185D+03,0.1710D+03,0.1335D+03,0.1016D+03,
     & 0.7120D+02,0.4390D+02,0.2470D+02,0.1390D+02,0.7315D+01,
     & 0.3045D+01,0.9605D+00,0.3030D+00,0.8810D-01,0.1663D-01/)
      REAL*8, PARAMETER, DIMENSION(LCOalt) ::  
     &     BrOxaltIN = (/1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,
     &     1.d-2,1.d-2,1.d-2,1.d-2,0.12d0,0.12d0,0.12d0,0.12d0,0.06d0,
     &     0.06d0,0.06d0,0.06d0,0.06d0,0.06d0,0.06d0,0.06d0/)
     &     ,ClOxaltIN = (/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     &     1.d0,1.d0,8.d0,8.d0,8.d0,8.d0,8.d1,8.d1,8.d1,8.d1,8.d0,8.d0,
     &     8.d0,8.d0/)
     &     ,ClONO2altIN = (/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     &     1.d0,1.d0,1.d0,1.d0,1.d0,5.d1,5.d1,5.d1,5.d1,5.d1,5.d1,5.d1,
     &     5.d1,5.d1,5.d1/)
     &     ,HClaltIN = (/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     &     1.d0,1.d0,2.5d1,4.0d1,9.0d1,1.7d2,1.9d2,2.5d2,2.5d2,2.5d2,
     &     2.5d2,2.5d2,2.5d2,2.5d2/)
      REAL*8, PARAMETER, DIMENSION(LCH4alt) :: PCH4alt = 
     &                     (/569d0, 150d0, 100d0, 32d0, 3.2d0, 0.23d0/)
      REAL*8, PARAMETER, DIMENSION(LCH4alt) ::   
     &   CH4altINT =(/1.79d0, 1.75d0, 1.620d0,1.460d0,0.812d0,0.230d0/),
     &   CH4altINX =(/1.79d0, 1.75d0, 1.440d0,1.130d0,0.473d0,0.202d0/)    
     
!@dbparam Tpsc_offset_N NH offset for the above T_thresh
!@dbparam Tpsc_offset_S SH offset for the above T_thresh
!@dbparam reg1Power_SpherO2andN2Ocorr first from surface region power of 
!@+ cos(sza)^x of spherical correction to ss(rj%O2__O_O) and ss(rj%N2O__M_O1D)
!@dbparam reg2Power_SpherO2andN2Ocorr second from surface region power of 
!@+ cos(sza)^x of spherical correction to ss(rj%O2__O_O) and ss(rj%N2O__M_O1D)
!@dbparam reg3Power_SpherO2andN2Ocorr third from surface region power of 
!@+ cos(sza)^x of spherical correction to ss(rj%O2__O_O) and ss(rj%N2O__M_O1D)
!@dbparam reg4Power_SpherO2andN2Ocorr fourth and last from surface region power of 
!@+ cos(sza)^x of spherical correction to ss(rj%O2__O_O) and ss(rj%N2O__M_O1D)
!@dbparam reg1TopPres_SpherO2andN2Ocorr pressure at top of first from surface region
!@+ for spherical correction to ss(rj%O2__O_O) and ss(rj%N2O__M_O1D) (hPa)
!@dbparam reg2TopPres_SpherO2andN2Ocorr pressure at top of second from surface region
!@+ for spherical correction to ss(rj%O2__O_O) and ss(rj%N2O__M_O1D) (hPa)
!@dbparam reg3TopPres_SpherO2andN2Ocorr pressure at top of third from surface region
!@+ for spherical correction to ss(rj%O2__O_O) and ss(rj%N2O__M_O1D) (hPa)
! (fourth = top region needs no upper pressure)
!@dbparam windowO2corr linear correction to ss(rj%O2__O_O) O2 in window region (in addition to spherical)
!@dbparam windowN2Ocorr linear correction to ss(rj%N2O__M_O1D) N2O in window region (in addition to spherical)
!@dbparam ch4_init_sh,ch4_init_nh initial methane conc. (ppmv) 
!@+       defaults are for 1990
!@dbparam allowSomeChemReinit (1=YES) to allow some chemistry variables
!@+       to cold-start even if tracers don't. Warning: this includes
!@+       model strat Q( ) spec. hum. reinitialization, and default =1!
!@dbparam fix_CH4_chemistry (1=YES 0=NO) whether or not to used a fixed
!@+       value for methane in the chemistry code. USE -1 for initial
!@+       conditions from file CH4_IC (but L>LS1-1 only now!)
!@+       but use_rad_CH4=1 overrides this.
!@dbparam scale_ch4_IC_file multiplicative factor of CH4 IC if 
!@+       fix_CH4_chemistry=-1 (but only above LS1-1 !)
!@dbparam use_rad_ch4 =1 replaces CH4 surface sources with L=1
!+        overwriting with radiation code values.
!@dbparam use_rad_n2o =1 as ch4 case above
!@dbparam use_rad_cfc =1 as ch4 case above
!@dbparam Lmax_rad_O3 model levels to use tracer Ox in rad code (if on)
!@dbparam Lmax_rad_CH4 model levels to use tracer CH4 in rad code(if on)
!@dbparam which_trop 1=ls1-1 is tropopause, 0=LTROPO(I,J) is tropopause
!@dbparam PI_run used to turn on (1) and off (0) use of PI_ratio*
!@dbparam PIratio_N to scale NOx, HNO3, N2O5, HO2NO2
!@+       initial conditions and stratospheric overwriting.
!@dbparam PIratio_CO_T to scale tropospheric CO IC and overwrite
!@dbparam PIratio_CO_S to scale stratospheric CO IC and overwrite
!@dbparam PIratio_other to scale PAN,Isoprene,AlkyNit,Alkenes,Paraffin
!@+       ,Terpenes
!@+       initial conditions and stratospheric overwriting.
!@dbparam PIratio_N2O preindustrial ratio for N2O ICs and L=1 overwrite
!@dbparam PIratio_CFC preindustrial ratio for CFC ICs and L=1 overwrite
!@+       with model time (JYEAR, JMON, JDAY) 
!@dbparam PltOx for pres<PltOx Ox, NOx, ClOx, and BrOx get overwritten

      INTEGER ::        fix_CH4_chemistry = 0
     &                 ,which_trop        = 0
     &                 ,PI_run            = 0
     &                 ,use_rad_ch4       = 0
     &                 ,use_rad_n2o       = 0
     &                 ,use_rad_cfc       = 0
     &                 ,Lmax_rad_O3       = LM ! not topLevelOfChemistry
     &                 ,Lmax_rad_CH4      = LM ! not topLevelOfChemistry
     &                 ,allowSomeChemReinit = 1
      REAL*8 ::             ch4_init_sh   = 1.750d0
     &                     ,ch4_init_nh   = 1.855d0
     &                     ,scale_ch4_IC_file= 1.d0 
     &                     ,PIratio_N     = 0.667d0
     &                     ,PIratio_CO_T  = 0.667d0
     &                     ,PIratio_CO_S  = 0.500d0
     &                     ,PIratio_other = 0.500d0
     &                     ,PIratio_N2O   = 0.896d0
     &                     ,PIratio_CFC   = 0.000d0
     &                     ,PltOx         = 0.000d0
     &                     ,Tpsc_offset_N = -10.d0
     &                     ,Tpsc_offset_S = -10.d0
     &                     ,reg1Power_SpherO2andN2Ocorr = 2.0d0
     &                     ,reg2Power_SpherO2andN2Ocorr = 2.0d0
     &                     ,reg3Power_SpherO2andN2Ocorr = 1.0d0
     &                     ,reg4Power_SpherO2andN2Ocorr = 0.5d0
     &                     ,reg1TopPres_SpherO2andN2Ocorr = 50.d0
     &                     ,reg2TopPres_SpherO2andN2Ocorr = 10.d0
     &                     ,reg3TopPres_SpherO2andN2Ocorr = 5.d0
     &                     ,windowN2Ocorr = 0.8d0
     &                     ,windowO2corr  = 0.8d0
     &                     ,PSClatS       = -50.d0
     &                     ,PSClatN       =  50.d0

      LOGICAL, PARAMETER :: luselb            = .false.

C**************  V  A  R  I  A  B  L  E  S *******************  
!@var topLevelOfChemistry the model level above which no chemistry is done
!@var nn name of species that reacts, as defined in the MOLEC file. The
!@+      first index denotes the reactant 1 or 2, and the second the reaction
!@+      number, as defined in the JPLRX file
!@var nnr same as nn, for products
!@var nps reaction index for production as defined in JPLPH, given the
!@+       accumulated ireac (one ireac element per reaction per unique
!@+       reactant)
!@var nds same as nps for destruction
!@var npnr same as nps for thermal reactions in JPLRX
!@var ndnr same as npnr for destruction
!@var kps index of JPLPH reaction (production) per photodissociating
!@+       species found in MOLEC
!@var kds same as kps for destruction
!@var kpnr same as kps for thermal reactions in JPLRX
!@var kdnr same as kpnr for destruction
!@var nst reverse reaction number for dissociation reactions
!@var lprn,jprn,iprn l, j, and i point for chemistry debugging
!@var ay name of gas being considered, as defined in MOLEC
!@var y concentration of gas, 1st index=gas number, 2nd=verticle level
!@var rr rate constant of chemical reaction, first index - reaction
!@+   number, 2nd is verticle level
!@var ss photodissociation coefficient, indicies; rxn #,L,I,J
!@var pe rate constant for bimolecular chemical reaction
!@var ea activation energy constant for bimolecular chemical reactions
!@var ro,r1,sn,sb rate parameters for trimolecular reactions
!@var conc concentration of optically important gases (O2 & O3), first
!@+   vertivle level, second=gas number (1=O2,2=O3)
!@var TXL temperature profile
!@var prnrts logical: print rate of each chemical reaction?
!@var prnchg logical: print chemical changes?
!@var prnls logical: print reaction lists by species?
!@var yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR?
!@var mNO2 3D vol mixing ratio of NO2 saved for subdaily diagnostics
!@var yCl2,yCl2O2 3D arrays to remember some non-tracer species...
!@var NCFASTJ number of levels in the fastj atmosphere
!@var MIEDX Type of aerosol scattering, currently 6 set up:
!@+   1=Rayly 2=iso 3=iso-equiv 4=bkgrd-sulf,5=volc-sulf,6=liq water
!@var PFASTJ pressure sent to FASTJ
!@var   Rayleigh parameters (effective cross-section) (cm2)
!@var odtmp Optical depth (temporary array)
!@var XLTAU    TTAU along the slant path
!@var XL      Slant path between points
!@var nfam number of beginning molecule of each chemical family
!@var SALBFJ surface albedo parameter from radiation to fastj
!@var OxICIN Ox initial conditions (unit=PPPM,LCOalt levels)
!@var OxICINL column version of OxICIN
!@var COICIN CO initial conditions (unit=PPPM,LCOalt levels)
!@var COICINL column version of OxICIN
!@var N2OICIN N2O initial conditions (unit=PPPM,LCOalt levels)
!@var N2OICINL column version of N2OICIN
!@var CH4ICIN CH4 initial conditions (unit=PPPM,LCOalt levels)
!@var CH4ICINL column version of CH4ICIN
!@var CFCICIN CFC initial conditions (unit=PPPM,LCOalt levels)
!@var CFCICINL column version of CFCICIN
!@var BrOxaltIN altitude dependence BrOx (unitless,LCOalt levels)
!@var ClOxaltIN altitude dependence ClOx (unitless,LCOalt levels)
!@var ClONO2altIN altitude dependence ClONO2 (unitless,LCOalt levels)
!@var HClaltIN altitude dependence HCl (unitless,LCOalt levels)
!@var CH4altINT tropical strat adjustments to CH4 (LCH4alt levels)
!@var CH4altINX xtra-tropical strat adjustments to CH4 LCH4alt levels)
!@var OxIC Ox initial conditions (unit=KG,LM levels)
!@var OxICL column version of OxIC
!@var COIC CO initial conditions (unit=KG,LM levels)
!@var COICL column version of COIC
!@var N2OICX N2O initial conditions (unit=KG,LM levels) X=not Jean's
!@var N2OICL column version of N2OICX
!@var CH4ICX CH4 initial conditions (unit=KG,LM levels) X=not Jean's
!@var CH4ICL column version of CH4ICX
!@var CFCIC CFC initial conditions (unit=KG,LM levels)
!@var CFCICL column version of CFCIC
!@var BrOxalt altitude dependence BrOx (unitless,LM levels)
!@var ClOxalt altitude dependence ClOx (unitless,LM levels)
!@var ClONO2alt altitude dependence ClONO2 (unitless,LM levels)
!@var HClalt altitude dependence HCl (unitless,LM levels)
!@var CH4altT tropical strat adjustments to CH4 (unitless, LM levels)
!@var CH4altX xtra-tropical strat adjustments to CH4 (LM levels)
!@var BYFJM = 1/JM
!@var TX temperature variable for master chem
!@var ta local array to hold temperature
!@var rh local array to hold relative humidity
!@var FASTJLAT,FASTJLON latitude & LONGITUDE (degrees) for use in fastj
!@var sulfate N2O5 sulfate sink (formerly SRC(I,J,L,20) variable)   
!@var dms_offline DMS concentration for HOx sink reactions
!@var so2_offline SO2 concentration for HOx conversion reactions
!@var prod_sulfate  N2O5 change by sulfate reactions in mass units
!@var wprod_sulf N2O5 change by sulfate reactions in molecules/cm3/s
!@var DT2 variable chemical time step, set in masterchem
!@var ratioNs,ratioN2,rNO2frac,rNOfrac,rNOdenom variables for nitrogen
!@+   conservation (strat)
!@var chemrate reaction rate per layer
!@var photrate photolysis rate per layer
!@var L75P first model level above nominal 75 hPa
!@var L75M first model level below nominal 75 hPa
!@var F75P interpolation coeff. of higher altitude value (units ln(P))
!@var F75M interpolation coeff. of lower altitude value (units ln(P))
!@var L569P first model level above nominal 569 hPa
!@var L569M first model level below nominal 569 hPa
!@var F569P interpolation coeff. of higher altitude value (units ln(P))
!@var F569M interpolation coeff. of lower altitude value (units ln(P))
!@var DU_O3 total column ozone in latitude band
!@var SF3 is H2O photolysis in Schumann-Runge Bands
!@var SF2 is NO photolysis in Schumann-Runge Bands
!@var Jacet photolysis rate for acetone (not done through fastj)
!@var acetone acetone column mixing ratio for the curren I,J (static for now)
!@var pscX column logical for the existance of polar strat clouds(PSCs)
!@var save_NO2column instantaneous NO2 column (for SUBDD exporting)
!@var RGAMMASULF N2O5-->HNO3 conversion on aerosols?
!@var changeL 2D array holds the local change due to chem until
!@+   adding to tr3Dsource
!@var bythick recipricol thickness of each layer (1/m) saved on
!@+ model layers.
!@var ClOx_old total ClOx at start of chemical timestep
!@var aero yes(1) or no(0) tag of non-zero rkext from Crates
!@var mostRecentNonZeroAlbedo remembers last time that ALB(I,J,1) was non-zer
!@+ for given I,J point (saved in rsf for reproducibilty purposes)
      INTEGER :: L75P,L75M,L569P,L569M,
     &lprn,jprn,iprn,MIEDX,NCFASTJ,topLevelOfChemistry
      INTEGER, DIMENSION(numfam+1)     :: nfam = (/0,0,0,0,ny+1/)
      INTEGER, DIMENSION(p_1,n_rx)     :: nn, nnr
      INTEGER, DIMENSION(p_1*n_rx)     :: npnr, ndnr
      INTEGER, DIMENSION(p_1*n_rj)     :: nps, nds
      INTEGER, DIMENSION(nc)           :: kps, kds, kpnr, kdnr
      INTEGER, DIMENSION(n_nst)        :: nst
      INTEGER, ALLOCATABLE, DIMENSION(:) :: aero

C**************  Latitude-Dependant (allocatable) *******************
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: DU_O3
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: acetone
#ifdef TRACERS_dCO
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: d17Oacetone
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: d18Oacetone
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: d13Cacetone
#endif  /* TRACERS_dCO */
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: ss
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: yNO3,pHOx,pNOx,pOx,
     & yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR,TX,sulfate,OxIC,
#ifdef TRACERS_dCO
     & ydC217O3,ydC218O3,yd13C2O3,
     & yd13CXPAR,
     & yd17OROR,yd18OROR,yd13CROR,
     & yd17Oald,yd18Oald,yd13Cald,
     & ydCH317O2,ydCH318O2,yd13CH3O2,
#endif  /* TRACERS_dCO */
     & CH4ICX,dms_offline,so2_offline,yso2,ydms,mNO2,COIC,pNO3
     & ,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2,N2OICX,CFCIC,SF3,SF2
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: COICIN,OxICIN,CH4ICIN
     &                                       ,N2OICIN,CFCICIN
      REAL*8, ALLOCATABLE, DIMENSION(:,:):: save_NO2column
      REAL*8, ALLOCATABLE, DIMENSION(:,:):: mostRecentNonZeroAlbedo

C**************  Not Latitude-Dependant ****************************      
      REAL*8 :: XLTAU,BYFJM,
     & FASTJLAT,FASTJLON,DT2,F75P,F75M,F569P,F569M,RGAMMASULF
     & ,ratioNs,ratioN2,rNO2frac,rNOfrac,rNOdenom
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: y
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rr
      REAL*8, ALLOCATABLE, DIMENSION(:)   :: odtmp,ta,Jacet,rh,bythick
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: chemrate, photrate
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: dest, prod
      REAL*8, ALLOCATABLE, DIMENSION(:)   :: OxlossbyH, ClOx_old
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: changeL
      REAL*8, DIMENSION(n_bi+n_nst)       :: pe, ea
      REAL*8, DIMENSION(n_tri)            :: ro, r1, sn, sb
      REAL*8, DIMENSION(LCOalt)           :: COICINL,OxICINL,CH4ICINL
     &                                       ,N2OICINL,CFCICINL
      REAL*8, DIMENSION(LM)  :: CH4altT,CH4altX,COICL,OxICL,CH4ICL ! stays LM
     &                        ,BrOxalt,ClOxalt,ClONO2alt,HClalt
     &                        ,N2OICL,CFCICL  

      LOGICAL                             :: prnrts,prnchg,prnls
      LOGICAL, ALLOCATABLE, DIMENSION(:)  :: pscX

      CHARACTER*8, DIMENSION(nc)          :: ay
      
!@dbparam tune_NOx Factor to multiply NOx emissions with
      real*8 :: tune_NOx = 1.d0
!@dbparam tune_BVOC Factor to multiply biogenic VOC emissions with
      real*8 :: tune_BVOC = 1.d0
!@dbparam NOx_yr Year of NOx emissions to use
      integer :: NOx_yr=0
!@dbparam CO_yr Year of CO emissions to use
      integer :: CO_yr=0
!@dbparam VOC_yr Year of anthropogenic VOC emissions to use
      integer :: VOC_yr=0

      END MODULE TRCHEM_Shindell_COM
      
      
      
      subroutine alloc_trchem_shindell_com(grid)
!@SUM  To allocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth G.Faluvegi
      use Dictionary_mod, only : get_param
      use domain_decomp_atm, only: dist_grid, getDomainBounds
      use resolution, only: im,lm,Plbot
      use tracer_com, only: ntm
      use TRCHEM_Shindell_COM, only: DU_O3,ss,yNO3,
     & pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR,
     & TX,sulfate,COIC,OxIC,CH4ICX,dms_offline,so2_offline,yso2,ydms,
#ifdef TRACERS_dCO
     & ydC217O3,ydC218O3,yd13C2O3,
     & yd13CXPAR,
     & yd17OROR,yd18OROR,yd13CROR,
     & yd17Oald,yd18Oald,yd13Cald,
     & ydCH317O2,ydCH318O2,yd13CH3O2,
     & d17Oacetone,d18Oacetone,d13Cacetone,
#endif  /* TRACERS_dCO */
     & COICIN,OxICIN,CH4ICIN,n_rj,LCOalt,acetone,mNO2,
     & save_NO2column,pNO3
     & ,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2,N2OICX,CFCIC,SF3,SF2,
     & N2OICIN,CFCICIN,y,rr,odtmp,ta,Jacet,chemrate,photrate,dest,prod,
     & OxlossbyH,pscX,nc,n_rx,ny,changeL,rh,bythick,ClOx_old,aero

      use TRCHEM_Shindell_COM, only: topLevelOfChemistry ! define here
      use TRCHEM_Shindell_COM, only: mostRecentNonZeroAlbedo

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H, L
      logical :: init = .false.
      real*8  :: x

      if(init)return
      init=.true.
    
      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO

      ! First, determine the number of layers where we'll do
      ! chemistry. (Currently, this includes photolysis, as
      ! topLevelOfChemistry sets NLGCM). The value here is
      ! set based on a pressure criterion. Note that this is
      ! with respect to PLbot array only, as sigmas were not
      ! yet defined at this point in the model. The 0.1 is 
      ! chosen below such that chemistry will use all levels of
      ! e.g. the traditional AR5 40-layer model:
  
      topLevelOfChemistry=0
      do L=1,LM
        ! next line tries to prevent picking wrong number of
        ! layers when e.g. the model top is 0.01=0.0099999993423432
        ! because of representation errror: 
        x=nint(PLbot(L+1) * 1.E3)*1.E-3 
        if(x>=0.1d0)topLevelOfChemistry=L
      end do

      ! Allow user to override this level from the rundeck:
      call get_param('override_LM_chem',topLevelOfChemistry,
     &               default=topLevelOfChemistry)

      if(topLevelOfChemistry == 0 .or. topLevelOfChemistry>LM)
     & call stop_model(
     & 'topLevelOfChemistry not determined.',255)
      ! Think about whether you want to also set
      ! Lmax_rad_O3 and Lmax_rad_CH4 to topLevelOfChemistry
      ! here. At the moment, I think no.

      ! Things allocatable *because* of above-chemistry model layers:
      allocate(        y(nc,   topLevelOfChemistry) )
      allocate(       rr(n_rx, topLevelOfChemistry) )
      allocate(    odtmp(      topLevelOfChemistry) )
      allocate(       ta(      topLevelOfChemistry) )
      allocate(       rh(      topLevelOfChemistry) )
      allocate(  bythick(      topLevelOfChemistry) )
      allocate( ClOx_old(      topLevelOfChemistry) )
      allocate(    Jacet(      topLevelOfChemistry) )
      allocate(     aero(      topLevelOfChemistry) )
      allocate(     pscX(      topLevelOfChemistry) )
      allocate( chemrate(n_rx, topLevelOfChemistry) )
      allocate( photrate(n_rj, topLevelOfChemistry) )
      allocate(     dest(ny,   topLevelOfChemistry) )
      allocate(     prod(ny,   topLevelOfChemistry) )
      allocate(OxlossbyH(      topLevelOfChemistry) )
      allocate(  changeL(      topLevelOfChemistry, ntm) )

      ! Normally allocated things:
      allocate(save_NO2column(I_0H:I_1H,J_0H:J_1H) )
      allocate(         DU_O3(          J_0H:J_1H) )
      allocate(ss(n_rj, topLevelOfChemistry,
     &                      I_0H:I_1H,J_0H:J_1H) )
      allocate(     acetone(topLevelOfChemistry) )
#ifdef TRACERS_dCO
      allocate( d17Oacetone(topLevelOfChemistry) )
      allocate( d18Oacetone(topLevelOfChemistry) )
      allocate( d13Cacetone(topLevelOfChemistry) )
#endif  /* TRACERS_dCO */
      allocate(        yNO3(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(        pHOx(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(        pNOx(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(        pNO3(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(         pOx(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(      yCH3O2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#ifdef TRACERS_dCO
      allocate(   ydCH317O2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(   ydCH318O2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(   yd13CH3O2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#endif  /* TRACERS_dCO */
      allocate(       yC2O3(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#ifdef TRACERS_dCO
      allocate(    ydC217O3(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(    ydC218O3(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(    yd13C2O3(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#endif  /* TRACERS_dCO */
      allocate(        yROR(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#ifdef TRACERS_dCO
      allocate(    yd17OROR(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(    yd18OROR(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(    yd13CROR(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#endif  /* TRACERS_dCO */
      allocate(        yXO2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(   yAldehyde(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#ifdef TRACERS_dCO
      allocate(    yd17Oald(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(    yd18Oald(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(    yd13Cald(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#endif  /* TRACERS_dCO */
      allocate(       yXO2N(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(      yRXPAR(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#ifdef TRACERS_dCO
      allocate(   yd13CXPAR(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
#endif  /* TRACERS_dCO */
      allocate(        yso2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(        ydms(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(       pClOx(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(        pClx(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(      pOClOx(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) ) 
      allocate(       pBrOx(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(        yCl2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) ) 
      allocate(      yCl2O2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(         SF3(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(         SF2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) )
      allocate(        mNO2(I_0H:I_1H,J_0H:J_1H,topLevelOfChemistry) ) ! set to undef above chem in DIAG.f
      allocate(      OxICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
      allocate(      COICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
      allocate(     N2OICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
      allocate(     CFCICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
      allocate(     CH4ICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
      allocate(        OxIC(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        COIC(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(       CFCIC(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(      CH4ICX(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(      N2OICX(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(          TX(I_0H:I_1H,J_0H:J_1H,LM)      ) 
      allocate( dms_offline(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate( so2_offline(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(     sulfate(I_0H:I_1H,J_0H:J_1H,LM)      ) ! could be read from 3D file

      allocate( mostRecentNonZeroAlbedo(I_0H:I_1H,J_0H:J_1H))
      mostRecentNonZeroAlbedo=0.d0
      
      return
      end subroutine alloc_trchem_shindell_com
      
      subroutine set_rrate_index(irr, ate)
!@sum Dynamically assign reaction rate indices to variables, for use in
!@+   chemistry. The reactions are denoted REAC1_REAC2__PROD1_PROD2, and
!@+   the reaction types are rrmono (monomolecular breakups), rrbi
!@+   bimolecular), rrtri (trimolecular), and rrhet (heterogeneous).
!@auth Kostas Tsigaridis
      use trchem_shindell_com, only: rrmono, rrbi,rrtri,rrhet
      implicit none

      integer, intent(in) :: irr
      character(len=8), dimension(4), intent(in) :: ate
      character(len=36) :: reaction

      reaction = trim(ate(1))//'_'//trim(ate(2))//'__'//
     &           trim(ate(3))//'_'//trim(ate(4))

      select case(reaction)

! monomolecular break ups
        case('HO2NO2_M__HO2_NO2')
          rrmono%HO2NO2_M__HO2_NO2=irr
        case('N2O5_M__NO3_NO2')
          rrmono%N2O5_M__NO3_NO2=irr
        case('Cl2O2_M__ClO_ClO')
          rrmono%Cl2O2_M__ClO_ClO=irr

! bimolecular reactions
        case('OH_HO2__H2O_O2')
          rrbi%OH_HO2__H2O_O2=irr
        case('OH_O3__HO2_O2')
          rrbi%OH_O3__HO2_O2=irr
        case('OH_OH__H2O_O')
          rrbi%OH_OH__H2O_O=irr
        case('HO2_O3__OH_O2')
          rrbi%HO2_O3__OH_O2=irr
        case('O3_NO__NO2_O2')
          rrbi%O3_NO__NO2_O2=irr
        case('HO2_NO__OH_NO2')
          rrbi%HO2_NO__OH_NO2=irr
        case('NO2_O3__NO3_O2')
          rrbi%NO2_O3__NO3_O2=irr
        case('O(1D)_O2__O_O2')
          rrbi%O1D_O2__O_O2=irr
        case('O(1D)_M__O_M')
          rrbi%O1D_M__O_M=irr
        case('O(1D)_H2O__OH_OH')
          rrbi%O1D_H2O__OH_OH=irr
        case('O(1D)_CH4__OH_CH3O2')
          rrbi%O1D_CH4__OH_CH3O2=irr
        case('CH4_OH__H2O_CH3O2')
          rrbi%CH4_OH__H2O_CH3O2=irr
        case('CO_OH__HO2_O2')
          rrbi%CO_OH__HO2_O2=irr
        case('OH_H2O2__H2O_HO2')
          rrbi%OH_H2O2__H2O_HO2=irr
        case('HO2_HO2__H2O2_O2')
          rrbi%HO2_HO2__H2O2_O2=irr
        case('OH_HNO3__H2O_NO3')
          rrbi%OH_HNO3__H2O_NO3=irr
        case('NO3_NO__NO2_NO2')
          rrbi%NO3_NO__NO2_NO2=irr
        case('OH_HO2NO2__H2O_NO2')
          rrbi%OH_HO2NO2__H2O_NO2=irr
        case('H2_OH__HO2_H2O')
          rrbi%H2_OH__HO2_H2O=irr
        case('CH3O2_NO__HCHO_NO2')
          rrbi%CH3O2_NO__HCHO_NO2=irr
        case('HCHO_OH__HO2_CO')
          rrbi%HCHO_OH__HO2_CO=irr
        case('CH3O2_HO2__CH3OOH_O2')
          rrbi%CH3O2_HO2__CH3OOH_O2=irr
        case('CH3OOH_OH__CH3O2_H2O')
          rrbi%CH3OOH_OH__CH3O2_H2O=irr
        case('NO2_NO3__NO_NO2')
          rrbi%NO2_NO3__NO_NO2=irr
        case('NO3_NO3__NO2_NO2')
          rrbi%NO3_NO3__NO2_NO2=irr
        case('O_NO2__NO_O2')
          rrbi%O_NO2__NO_O2=irr
        case('CH3O2_CH3O2__HCHO_HCHO')
          rrbi%CH3O2_CH3O2__HCHO_HCHO=irr
        case('NO3_HCHO__HNO3_CO')
          rrbi%NO3_HCHO__HNO3_CO=irr
        case('PAN_M__C2O3_NO2')
          rrbi%PAN_M__C2O3_NO2=irr
        case('Isoprene_OH__HCHO_Alkenes')
          rrbi%Isoprene_OH__HCHO_Alkenes=irr
        case('Isoprene_O3__HCHO_Alkenes')
          rrbi%Isoprene_O3__HCHO_Alkenes=irr
        case('Isoprene_NO3__HO2_Alkenes')
          rrbi%Isoprene_NO3__HO2_Alkenes=irr
        case('AlkylNit_OH__NO2_XO2')
          rrbi%AlkylNit_OH__NO2_XO2=irr
        case('Alkenes_OH__HCHO_HO2')
          rrbi%Alkenes_OH__HCHO_HO2=irr
        case('Alkenes_O3__HCHO_CO')
          rrbi%Alkenes_O3__HCHO_CO=irr
        case('Alkenes_NO3__HCHO_NO2')
          rrbi%Alkenes_NO3__HCHO_NO2=irr
        case('Paraffin_OH__HO2_M')
          rrbi%Paraffin_OH__HO2_M=irr
        case('Paraffin_RXPAR__M_M')
          rrbi%Paraffin_RXPAR__M_M=irr
        case('Aldehyde_OH__C2O3_M')
          rrbi%Aldehyde_OH__C2O3_M=irr
        case('C2O3_NO__HCHO_NO2')
          rrbi%C2O3_NO__HCHO_NO2=irr
        case('C2O3_C2O3__HCHO_HCHO')
          rrbi%C2O3_C2O3__HCHO_HCHO=irr
        case('C2O3_HO2__HCHO_HO2')
          rrbi%C2O3_HO2__HCHO_HO2=irr
        case('ROR_M__Aldehyde_HO2')
          rrbi%ROR_M__Aldehyde_HO2=irr
        case('ROR_M__HO2_M')
          rrbi%ROR_M__HO2_M=irr
        case('XO2_NO__NO2_M')
          rrbi%XO2_NO__NO2_M=irr
        case('XO2_HO2__CH3OOH_O2')
          rrbi%XO2_HO2__CH3OOH_O2=irr
        case('XO2_XO2__M_M')
          rrbi%XO2_XO2__M_M=irr
        case('XO2N_NO__AlkylNit_M')
          rrbi%XO2N_NO__AlkylNit_M=irr
        case('XO2N_HO2__CH3OOH_O2')
          rrbi%XO2N_HO2__CH3OOH_O2=irr
        case('Cl_O3__ClO_O2')
          rrbi%Cl_O3__ClO_O2=irr
        case('ClO_O__Cl_O2')
          rrbi%ClO_O__Cl_O2=irr
        case('Cl_OClO__ClO_ClO')
          rrbi%Cl_OClO__ClO_ClO=irr
        case('ClO_O3__Cl_O2')
          rrbi%ClO_O3__Cl_O2=irr
        case('ClO_O3__OClO_O2')
          rrbi%ClO_O3__OClO_O2=irr
        case('O_OClO__ClO_O2')
          rrbi%O_OClO__ClO_O2=irr
        case('OH_Cl2__HOCl_Cl')
          rrbi%OH_Cl2__HOCl_Cl=irr
        case('OH_HCl__H2O_Cl')
          rrbi%OH_HCl__H2O_Cl=irr
        case('OH_HOCl__H2O_ClO')
          rrbi%OH_HOCl__H2O_ClO=irr
        case('O_HCl__OH_Cl')
          rrbi%O_HCl__OH_Cl=irr
        case('O_HOCl__OH_ClO')
          rrbi%O_HOCl__OH_ClO=irr
        case('OClO_OH__HOCl_O2')
          rrbi%OClO_OH__HOCl_O2=irr
        case('Cl_HOCl__Cl2_OH')
          rrbi%Cl_HOCl__Cl2_OH=irr
        case('Cl_H2O2__HCl_HO2')
          rrbi%Cl_H2O2__HCl_HO2=irr
        case('Cl_HO2__HCl_O2')
          rrbi%Cl_HO2__HCl_O2=irr
        case('Cl_HO2__OH_ClO')
          rrbi%Cl_HO2__OH_ClO=irr
        case('ClO_OH__HO2_Cl')
          rrbi%ClO_OH__HO2_Cl=irr
        case('ClO_OH__HCl_O2')
          rrbi%ClO_OH__HCl_O2=irr
        case('ClO_HO2__HOCl_O2')
          rrbi%ClO_HO2__HOCl_O2=irr
        case('ClO_NO__NO2_Cl')
          rrbi%ClO_NO__NO2_Cl=irr
        case('ClONO2_O__ClO_NO3')
          rrbi%ClONO2_O__ClO_NO3=irr
        case('HCl_O(1D)__Cl_OH')
          rrbi%HCl_O1D__Cl_OH=irr
        case('NO_OClO__NO2_ClO')
          rrbi%NO_OClO__NO2_ClO=irr
        case('HBr_OH__H2O_Br')
          rrbi%HBr_OH__H2O_Br=irr
        case('BrO_O__Br_O2')
          rrbi%BrO_O__Br_O2=irr
        case('Br_O3__BrO_O2')
          rrbi%Br_O3__BrO_O2=irr
        case('BrO_NO__Br_NO2')
          rrbi%BrO_NO__Br_NO2=irr
        case('Br_HO2__HBr_O2')
          rrbi%Br_HO2__HBr_O2=irr
        case('BrO_HO2__HOBr_O2')
          rrbi%BrO_HO2__HOBr_O2=irr
        case('Br_OClO__BrO_ClO')
          rrbi%Br_OClO__BrO_ClO=irr
        case('BrO_ClO__OClO_Br')
          rrbi%BrO_ClO__OClO_Br=irr
        case('BrO_ClO__Br_Cl')
          rrbi%BrO_ClO__Br_Cl=irr
        case('O(1D)_CH4__HCHO_H2')
          rrbi%O1D_CH4__HCHO_H2=irr
        case('BrO_BrO__Br_Br')
          rrbi%BrO_BrO__Br_Br=irr
        case('Br_H2O2__HBr_HO2')
          rrbi%Br_H2O2__HBr_HO2=irr
        case('BrO_OH__Br_HO2')
          rrbi%BrO_OH__Br_HO2=irr
        case('BrO_OH__HBr_O2')
          rrbi%BrO_OH__HBr_O2=irr
        case('Cl_CH4__HCl_CH3O2')
          rrbi%Cl_CH4__HCl_CH3O2=irr
        case('Cl_H2__HCl_HO2')
          rrbi%Cl_H2__HCl_HO2=irr
        case('O_HBr__OH_Br')
          rrbi%O_HBr__OH_Br=irr
        case('ClO_CH3O2__Cl_HCHO')
          rrbi%ClO_CH3O2__Cl_HCHO=irr
        case('N2O_O(1D)__N2_O2')
          rrbi%N2O_O1D__N2_O2=irr
        case('N2O_O(1D)__NO_NO')
          rrbi%N2O_O1D__NO_NO=irr
        case('O_O3__O2_O2')
          rrbi%O_O3__O2_O2=irr
        case('O_OH__O2_H')
          rrbi%O_OH__O2_H=irr
        case('O_HO2__OH_O2')
          rrbi%O_HO2__OH_O2=irr
        case('HO2_NO__HNO3_M')
          rrbi%HO2_NO__HNO3_M=irr
#ifdef TRACERS_TERP
        case('Terpenes_OH__HCHO_Alkenes')
          rrbi%Terpenes_OH__HCHO_Alkenes=irr
        case('Terpenes_O3__HCHO_Alkenes')
          rrbi%Terpenes_O3__HCHO_Alkenes=irr
        case('Terpenes_NO3__HO2_Alkenes')
          rrbi%Terpenes_NO3__HO2_Alkenes=irr
#endif  /* TRACERS_TERP */
#ifdef TRACERS_dCO
        case('O(1D)_CH4__OH_dCH317O2')
          rrbi%O1D_CH4__OH_dCH317O2=irr
        case('CH4_OH__H2O_dCH317O2')
          rrbi%CH4_OH__H2O_dCH317O2=irr
        case('dC17O_OH__HO2_O2')
          rrbi%dC17O_OH__HO2_O2=irr
        case('dCH317O2_NO__dHCH17O_NO2')
          rrbi%dCH317O2_NO__dHCH17O_NO2=irr
        case('dHCH17O_OH__HO2_dC17O')
          rrbi%dHCH17O_OH__HO2_dC17O=irr
        case('dCH317O2_HO2__dMe17OOH_O2')
          rrbi%dCH317O2_HO2__dMe17OOH_O2=irr
        case('dMe17OOH_OH__dCH317O2_H2O')
          rrbi%dMe17OOH_OH__dCH317O2_H2O=irr
        case('dCH317O2_CH3O2__dHCH17O_HCHO')
          rrbi%dCH317O2_CH3O2__dHCH17O_HCHO=irr
        case('CH3O2_dCH317O2__HCHO_dHCH17O')
          rrbi%CH3O2_dCH317O2__HCHO_dHCH17O=irr
        case('NO3_dHCH17O__HNO3_dC17O')
          rrbi%NO3_dHCH17O__HNO3_dC17O=irr
        case('d17OPAN_M__dC217O3_NO2')
          rrbi%d17OPAN_M__dC217O3_NO2=irr
        case('Isoprene_OH__dHCH17O_Alkenes')
          rrbi%Isoprene_OH__dHCH17O_Alkenes=irr
        case('Isoprene_O3__dHCH17O_Alkenes')
          rrbi%Isoprene_O3__dHCH17O_Alkenes=irr
        case('Alkenes_OH__dHCH17O_HO2')
          rrbi%Alkenes_OH__dHCH17O_HO2=irr
        case('Alkenes_O3__dHCH17O_CO')
          rrbi%Alkenes_O3__dHCH17O_CO=irr
        case('Alkenes_O3__HCHO_dC17O')
          rrbi%Alkenes_O3__HCHO_dC17O=irr
        case('Alkenes_NO3__dHCH17O_NO2')
          rrbi%Alkenes_NO3__dHCH17O_NO2=irr
        case('d17Oald_OH__dC217O3_M')
          rrbi%d17Oald_OH__dC217O3_M=irr
        case('dC217O3_NO__dHCH17O_NO2')
          rrbi%dC217O3_NO__dHCH17O_NO2=irr
        case('dC217O3_dC217O3__dHCH17O_dHCH17O')
          rrbi%dC217O3_dC217O3__dHCH17O_dHCH17O=irr
        case('dC217O3_HO2__dHCH17O_HO2')
          rrbi%dC217O3_HO2__dHCH17O_HO2=irr
        case('d17OROR_M__d17Oald_HO2')
          rrbi%d17OROR_M__d17Oald_HO2=irr
        case('d17OROR_M__HO2_M')
          rrbi%d17OROR_M__HO2_M=irr
        case('O(1D)_CH4__dHCH17O_H2')
          rrbi%O1D_CH4__dHCH17O_H2=irr
        case('Cl_CH4__HCl_dCH317O2')
          rrbi%Cl_CH4__HCl_dCH317O2=irr
        case('ClO_dCH317O2__Cl_dHCH17O')
          rrbi%ClO_dCH317O2__Cl_dHCH17O=irr
        case('Terpenes_OH__dHCH17O_Alkenes')
          rrbi%Terpenes_OH__dHCH17O_Alkenes=irr
        case('Terpenes_O3__dHCH17O_Alkenes')
          rrbi%Terpenes_O3__dHCH17O_Alkenes=irr

        case('O(1D)_CH4__OH_dCH318O2')
          rrbi%O1D_CH4__OH_dCH318O2=irr
        case('CH4_OH__H2O_dCH318O2')
          rrbi%CH4_OH__H2O_dCH318O2=irr
        case('dC18O_OH__HO2_O2')
          rrbi%dC18O_OH__HO2_O2=irr
        case('dCH318O2_NO__dHCH18O_NO2')
          rrbi%dCH318O2_NO__dHCH18O_NO2=irr
        case('dHCH18O_OH__HO2_dC18O')
          rrbi%dHCH18O_OH__HO2_dC18O=irr
        case('dCH318O2_HO2__dMe18OOH_O2')
          rrbi%dCH318O2_HO2__dMe18OOH_O2=irr
        case('dMe18OOH_OH__dCH318O2_H2O')
          rrbi%dMe18OOH_OH__dCH318O2_H2O=irr
        case('dCH318O2_CH3O2__dHCH18O_HCHO')
          rrbi%dCH318O2_CH3O2__dHCH18O_HCHO=irr
        case('CH3O2_dCH318O2__HCHO_dHCH18O')
          rrbi%CH3O2_dCH318O2__HCHO_dHCH18O=irr
        case('NO3_dHCH18O__HNO3_dC18O')
          rrbi%NO3_dHCH18O__HNO3_dC18O=irr
        case('d18OPAN_M__dC218O3_NO2')
          rrbi%d18OPAN_M__dC218O3_NO2=irr
        case('Isoprene_OH__dHCH18O_Alkenes')
          rrbi%Isoprene_OH__dHCH18O_Alkenes=irr
        case('Isoprene_O3__dHCH18O_Alkenes')
          rrbi%Isoprene_O3__dHCH18O_Alkenes=irr
        case('Alkenes_OH__dHCH18O_HO2')
          rrbi%Alkenes_OH__dHCH18O_HO2=irr
        case('Alkenes_O3__dHCH18O_CO')
          rrbi%Alkenes_O3__dHCH18O_CO=irr
        case('Alkenes_O3__HCHO_dC18O')
          rrbi%Alkenes_O3__HCHO_dC18O=irr
        case('Alkenes_NO3__dHCH18O_NO2')
          rrbi%Alkenes_NO3__dHCH18O_NO2=irr
        case('d18Oald_OH__dC218O3_M')
          rrbi%d18Oald_OH__dC218O3_M=irr
        case('dC218O3_NO__dHCH18O_NO2')
          rrbi%dC218O3_NO__dHCH18O_NO2=irr
        case('dC218O3_dC218O3__dHCH18O_dHCH18O')
          rrbi%dC218O3_dC218O3__dHCH18O_dHCH18O=irr
        case('dC218O3_HO2__dHCH18O_HO2')
          rrbi%dC218O3_HO2__dHCH18O_HO2=irr
        case('d18OROR_M__d18Oald_HO2')
          rrbi%d18OROR_M__d18Oald_HO2=irr
        case('d18OROR_M__HO2_M')
          rrbi%d18OROR_M__HO2_M=irr
        case('O(1D)_CH4__dHCH18O_H2')
          rrbi%O1D_CH4__dHCH18O_H2=irr
        case('Cl_CH4__HCl_dCH318O2')
          rrbi%Cl_CH4__HCl_dCH318O2=irr
        case('ClO_dCH318O2__Cl_dHCH18O')
          rrbi%ClO_dCH318O2__Cl_dHCH18O=irr
        case('Terpenes_OH__dHCH18O_Alkenes')
          rrbi%Terpenes_OH__dHCH18O_Alkenes=irr
        case('Terpenes_O3__dHCH18O_Alkenes')
          rrbi%Terpenes_O3__dHCH18O_Alkenes=irr

        case('O(1D)_CH4__OH_d13CH3O2')
          rrbi%O1D_CH4__OH_d13CH3O2=irr
        case('CH4_OH__H2O_d13CH3O2')
          rrbi%CH4_OH__H2O_d13CH3O2=irr
        case('d13CO_OH__HO2_O2')
          rrbi%d13CO_OH__HO2_O2=irr
        case('d13CH3O2_NO__dH13CHO_NO2')
          rrbi%d13CH3O2_NO__dH13CHO_NO2=irr
        case('dH13CHO_OH__HO2_d13CO')
          rrbi%dH13CHO_OH__HO2_d13CO=irr
        case('d13CH3O2_HO2__d13MeOOH_O2')
          rrbi%d13CH3O2_HO2__d13MeOOH_O2=irr
        case('d13MeOOH_OH__d13CH3O2_H2O')
          rrbi%d13MeOOH_OH__d13CH3O2_H2O=irr
        case('d13CH3O2_CH3O2__dH13CHO_HCHO')
          rrbi%d13CH3O2_CH3O2__dH13CHO_HCHO=irr
        case('CH3O2_d13CH3O2__HCHO_dH13CHO')
          rrbi%CH3O2_d13CH3O2__HCHO_dH13CHO=irr
        case('NO3_dH13CHO__HNO3_d13CO')
          rrbi%NO3_dH13CHO__HNO3_d13CO=irr
        case('d13CPAN_M__d13C2O3_NO2')
          rrbi%d13CPAN_M__d13C2O3_NO2=irr
        case('Isoprene_OH__dH13CHO_d13Calke')
          rrbi%Isoprene_OH__dH13CHO_d13Calke=irr
        case('Isoprene_O3__dH13CHO_d13Calke')
          rrbi%Isoprene_O3__dH13CHO_d13Calke=irr
        case('Isoprene_NO3__HO2_d13Calke')
          rrbi%Isoprene_NO3__HO2_d13Calke=irr
        case('d13Calke_OH__dH13CHO_HO2')
          rrbi%d13Calke_OH__dH13CHO_HO2=irr
        case('d13Calke_O3__dH13CHO_d13CO')
          rrbi%d13Calke_O3__dH13CHO_d13CO=irr
        case('d13Calke_NO3__dH13CHO_NO2')
          rrbi%d13Calke_NO3__dH13CHO_NO2=irr
        case('d13CPAR_OH__HO2_M')
          rrbi%d13CPAR_OH__HO2_M=irr
        case('d13CPAR_d13CXPAR__M_M')
          rrbi%d13CPAR_d13CXPAR__M_M=irr
        case('d13Cald_OH__d13C2O3_M')
          rrbi%d13Cald_OH__d13C2O3_M=irr
        case('d13C2O3_NO__dH13CHO_NO2')
          rrbi%d13C2O3_NO__dH13CHO_NO2=irr
        case('d13C2O3_d13C2O3__dH13CHO_dH13CHO')
          rrbi%d13C2O3_d13C2O3__dH13CHO_dH13CHO=irr
        case('d13C2O3_HO2__dH13CHO_HO2')
          rrbi%d13C2O3_HO2__dH13CHO_HO2=irr
        case('d13CROR_M__d13Cald_HO2')
          rrbi%d13CROR_M__d13Cald_HO2=irr
        case('d13CROR_M__HO2_M')
          rrbi%d13CROR_M__HO2_M=irr
        case('O(1D)_CH4__dH13CHO_H2')
          rrbi%O1D_CH4__dH13CHO_H2=irr
        case('Cl_CH4__HCl_d13CH3O2')
          rrbi%Cl_CH4__HCl_d13CH3O2=irr
        case('ClO_d13CH3O2__Cl_dH13CHO')
          rrbi%ClO_d13CH3O2__Cl_dH13CHO=irr
        case('Terpenes_OH__dH13CHO_d13Calke')
          rrbi%Terpenes_OH__dH13CHO_d13Calke=irr
        case('Terpenes_O3__dH13CHO_d13Calke')
          rrbi%Terpenes_O3__dH13CHO_d13Calke=irr
        case('Terpenes_NO3__HO2_d13Calke')
          rrbi%Terpenes_NO3__HO2_d13Calke=irr
#endif  /* TRACERS_dCO */

! trimolecular reactions
        case('O_O2__O3_M')
          rrtri%O_O2__O3_M=irr
        case('NO_O__NO2_M')
          rrtri%NO_O__NO2_M=irr
        case('OH_OH__H2O2_M')
          rrtri%OH_OH__H2O2_M=irr
        case('OH_NO2__HNO3_M')
          rrtri%OH_NO2__HNO3_M=irr
        case('HO2_NO2__HO2NO2_M')
          rrtri%HO2_NO2__HO2NO2_M=irr
        case('NO3_NO2__N2O5_M')
          rrtri%NO3_NO2__N2O5_M=irr
        case('OH_NO__HONO_M')
          rrtri%OH_NO__HONO_M=irr
        case('C2O3_NO2__PAN_M')
          rrtri%C2O3_NO2__PAN_M=irr
        case('ClO_ClO__Cl2O2_M')
          rrtri%ClO_ClO__Cl2O2_M=irr
        case('ClO_NO2__ClONO2_M')
          rrtri%ClO_NO2__ClONO2_M=irr
        case('BrO_NO2__BrONO2_M')
          rrtri%BrO_NO2__BrONO2_M=irr
#ifdef TRACERS_dCO
        case('dC217O3_NO2__d17OPAN_M')
          rrtri%dC217O3_NO2__d17OPAN_M=irr
        case('dC218O3_NO2__d18OPAN_M')
          rrtri%dC218O3_NO2__d18OPAN_M=irr
        case('d13C2O3_NO2__d13CPAN_M')
          rrtri%d13C2O3_NO2__d13CPAN_M=irr
#endif  /* TRACERS_dCO */

! heterogeneous reactions
        case('N2O5_H2O__HNO3_HNO3')
          rrhet%N2O5_H2O__HNO3_HNO3=irr
        case('ClONO2_H2O__HOCl_HNO3')
          rrhet%ClONO2_H2O__HOCl_HNO3=irr
        case('ClONO2_HCl__Cl_HNO3')
          rrhet%ClONO2_HCl__Cl_HNO3=irr
        case('HOCl_HCl__Cl_H2O')
          rrhet%HOCl_HCl__Cl_H2O=irr
        case('N2O5_HCl__Cl_HNO3')
          rrhet%N2O5_HCl__Cl_HNO3=irr
        case default
          call stop_model('Index for '//trim(reaction)//' missing',255)
      end select

      end subroutine set_rrate_index
