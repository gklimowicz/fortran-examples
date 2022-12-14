#include "rundeck_opts.h"
module TRACERS_SOA
!-------------------------------------------------------------------------------
! SOA formation - created by Kostas Tsigaridis
!-------------------------------------------------------------------------------
!@sum module for the calculation of secondary organic aerosols (SOA). Requires
!@+   tropospheric chemistry and aerosols to be activated.
!@auth Kostas Tsigaridis (ktsigaridis@giss.nasa.gov)
use RESOLUTION, only: LM
use DOMAIN_DECOMP_ATM,only: write_parallel, am_i_root
use OldTracer_mod, only: tr_mm
use TRACER_COM, only: NTM,nsoa,&
                      n_Isoprene,&

#ifdef TRACERS_TERP
                      n_Terpenes,&
                      n_apinp1g,n_apinp1a,n_apinp2g,n_apinp2a,&
#endif  /* TRACERS_TERP */
                      n_isopp1g,n_isopp1a,n_isopp2g,n_isopp2a
implicit none

!@var n_soa_i the first SOA-related species
!@var n_soa_e the last SOA-related species
integer                    :: n_soa_i,n_soa_e
!@var mw the molecular weight of all tracers, in units of tracer mass per mole
!@+      in order to replace the tm_mm which in some exceptional cases is in units of
!@+      a certain atom (typically C or S) per mole.
real*8, allocatable, dimension(:)     :: mw
!@var apartmass the mass-based yield of semivolatile species from chemistry
!@var apartmolar the molar-based yield of semivolatile species from chemistry
real*8, dimension(LM,nsoa) :: apartmass,apartmolar
!@var voc2nox factor deciding how much SOA is produced from the high-NOx and how much from
!@+           the low-NOx pathway, depending on the RO2+NO/RO2+HO2/RO2+RO2 branching ratio
!@+           as described by Lane et al., 2008 (ratio B). A value close to one means that
!@+           the high-NOx pathway dominates, while a value close to zero means that the
!@+           low-NOx pathway prevails. CH3O2 radical is not taken into account.
real*8, dimension(LM)      :: voc2nox
!@var apartmass_ref the low-NOx mass-based yield of semivolatile species from chemistry
!@var apartmass_nox_ref the high-NOx mass-based yield of semivolatile species from chemistry
real*8, dimension(nsoa)    :: apartmass_ref,apartmass_nox_ref
!@var  molec2ug converts molec/cm3 to ug/m3. 1.d0/molec2ug converts ug/m3 to molec/cm3
real*8, allocatable, dimension(:)     :: molec2ug
!@param LM_soa the uppermost level where chemical production of semivolatile gases is allowed
!WRONG ON PURPOSE
! SOA production from chemistry should be allowed everywhere, but
! due to convection Isoprene and Terpenes have a local maximum in the
! upper layers. This is unlikely to be the case in the real atmosphere,
! thus chemical (but not partitioning) production is disabled above level 7
integer, parameter :: LM_soa=LM!floor(float(LM)/3.)

!
! soa semivolatile products in aerosol phase
!
!@var whichsoa converts tracer index to soa index
integer, allocatable, dimension(:)     :: whichsoa
!@var issoa converts soa index to tracer index
integer, dimension(nsoa) :: issoa

!@param soacomp number of different cases in Lambda calculations
integer, parameter                 :: soacomp=11
!@var Lambda empirical factor describing the affinity of species i with species j
!@+          for the activity coefficient calculations. A value close to unity means
!@+          high chemical similarity of the species. A value of -1.0 declares
!@+          symmetric treatment: Lambda(i,j)=Lambda(j,i).
real*8, dimension(soacomp,soacomp) :: Lambda
integer, parameter                 :: imfapin=1
integer, parameter                 :: imfaro=2
integer, parameter                 :: imfisop=3
integer, parameter                 :: imfocii=4
integer, parameter                 :: imfocia=5
integer, parameter                 :: imfocb=6
integer, parameter                 :: imfococean=7
integer, parameter                 :: imfbcii=8
integer, parameter                 :: imfbcia=9
integer, parameter                 :: imfbcb=10
integer, parameter                 :: imfinorg=11
!
!                                     ter     aro    isop    ocii    ocia     ocb   ococean   bcii    bcia    bcb    inorg
data Lambda(imfapin,1:soacomp)    / 1.00d0, 1.00d0, 1.00d0, 0.80d0, 0.90d0, 0.85d0, 0.80d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfaro,1:soacomp)     /-1.00d0, 1.00d0, 1.00d0, 0.80d0, 0.90d0, 0.85d0, 0.80d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfisop,1:soacomp)    /-1.00d0,-1.00d0, 1.00d0, 0.80d0, 0.90d0, 0.85d0, 0.80d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfocii,1:soacomp)    /-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.90d0, 0.85d0, 0.80d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfocia,1:soacomp)    /-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.85d0, 0.80d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfocb,1:soacomp)     /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.80d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfococean,1:soacomp) /-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.90d0,-1.00d0, 1.00d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfbcii,1:soacomp)    /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfbcia,1:soacomp)    /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.75d0, 0.70d0/
data Lambda(imfbcb,1:soacomp)     /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.70d0/
data Lambda(imfinorg,1:soacomp)   /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0/
!
! enthalpies of vaporization (dH, in kJ/mol):
!    0.     No temperature dependence on vapor pressure
!   42.     Chung and Seinfeld, 2002; Kleindienst et al., GRL, 2007 (for isoprene products)
!   72.4    Pun et al., 2003
!   79.     Andersson-Skold and Simpson, 2001
!  109.     for pinic acid; Bilde and Pandis, EST, 2001
!  156.     Strader et al., 1999
!
!@param dH_isoprene enthalpy of vaporization for the isoprene-produced SOA species (KJ/mol)
real*8, parameter :: dH_isoprene=42.d0 ! Chung and Seinfeld, JGR, 2002; Henze and Seinfeld, GRL, 2006; Kleindienst et al., GRL, 2007
real*8, parameter :: dH_apinene=72.9d0
!real*8, parameter :: dH_terpenes_pinic=109.d0
!real*8, parameter :: dH_aromatics=72.9d0
!real*8, parameter :: dH_terpenes=42.d0
!real*8, parameter :: dH_terpenes_pinic=42.d0
!real*8, parameter :: dH_aromatics=15.d0

!@var kpart partitioning coefficient of SOA species (m3/ug)
real*8, dimension(LM,nsoa)         :: kpart
!@param kpart_ref Low-NOx partitioning coefficient of SOA species (m3/ug) at the reference temperature kpart_temp_ref
real*8, dimension(nsoa), parameter :: kpart_ref=(/&
! isoprene + OH low NOx SOAb formation, Henze and Seinfeld, 2006
1.62d0, 0.00862d0         & ! iisopp1a, iisopp2a ! WARNING!!! Species indices are inverted compared to the paper, for consistency:
                                                !            less volatile first, more volatile second
#ifdef TRACERS_TERP
! a-pinene + O3 low NOx SOAb formation, Presto et al., 2005
,1.d0/15.7d0, 1.d0/385d0  & ! iapinp1a, iapinp2a
#endif  /* TRACERS_TERP */
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!0.195,   0.003, & ! ibpinp1a_nox, ibpinp2a_nox: Griffin et al., JGR, 1999
!modelE!#  endif
!modelE!0.195,   0.003, & ! ibpinp1a_hox, ibpinp2a_hox: Griffin et al., JGR, 1999
!modelE!#endif
!modelE!#ifdef SOA_FULL
!modelE!0.053,   0.0019,& ! itolp1a_nox,  itolp2a_nox: Odum et al., Science, 1997
!modelE!#endif
!modelE!0.053,   0.0019 & ! itolp1a_hox,  itolp2a_hox: Odum et al., Science, 1997
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!,0.301,   0.008  & ! ixylp1a_nox,  ixylp2a_nox: Song et al., EST, 2005 !!!!!!!!!!!!!!CHANGE ALSO HARDCODED PART IN sources_sinks.f90!!!!!!!!!!!!!!!
!modelE!#  endif
!modelE!,0.229,   0.004  & ! ixylp1a_hox,  ixylp2a_hox: Song et al., EST, 2005 !!!!!!!!!!!!!!CHANGE ALSO HARDCODED PART IN sources_sinks.f90!!!!!!!!!!!!!!!
!modelE!#endif
!modelE!!#ifdef SOA_FULL
!modelE!!0.430,   0.047, & ! itolp1a_nox,  itolp2a_nox: Ng et al., ACP, 2007
!modelE!!#endif
!modelE!!1.e10,   1.e10  & ! itolp1a_hox,  itolp2a_hox: Ng et al., ACP, 2007
!modelE!!#ifndef SOA_MINIMUM
!modelE!!#  ifdef SOA_FULL
!modelE!!,0.761,   0.029  & ! ixylp1a_nox,  ixylp2a_nox: Ng et al., ACP, 2007
!modelE!!#  endif
!modelE!!,1.e10,   1.e10  & ! ixylp1a_hox,  ixylp2a_hox: Ng et al., ACP, 2007
!modelE!!#endif
/)

!@param kpart_nox_ref High-NOx partitioning coefficient of SOA species (m3/ug) at the reference temperature kpart_temp_nox_ref
real, dimension(nsoa), parameter :: kpart_nox_ref=(/&
1.62d0, 0.00862d0         & ! iisopp1a, iisopp2a ! WARNING!!! Species indices are inverted compared to the paper, for consistency:
                                                !            less volatile first, more volatile second
#ifdef TRACERS_TERP
,1.d0/15.7d0, 1.d0/385.d0 & ! iapinp1a, iapinp2a: Presto et al., EST, 2005
#endif  /* TRACERS_TERP */
!modelE!0.195,   0.003, & ! ibpinp1a_nox, ibpinp2a_nox: Griffin et al., JGR, 1999
!modelE!0.053,   0.0019,& ! itolp1a_nox,  itolp2a_nox: Odum et al., Science, 1997
!modelE!0.301,   0.008  & ! ixylp1a_nox,  ixylp2a_nox: Song et al., EST, 2005
!modelE!!1.e10,   1.e10, & ! itolp1a_nox,  itolp2a_nox: Ng et al., ACP, 2007
!modelE!!1.e10,   1.e10  & ! ixylp1a_nox,  ixylp2a_nox: Ng et al., ACP, 2007
/)

!@param kpart_temp_ref reference temperature where kpart_ref was derived
real*8, dimension(nsoa), parameter :: kpart_temp_ref=(/&
! isoprene + OH low NOx SOAb formation, Henze and Seinfeld, GRL 2006
295.d0,295.d0   & ! iisopp1a_hox, iisopp2a_hox
#ifdef TRACERS_TERP
,295.d0,295.d0  & ! iapinp1a_hox, iapinp2a_hox: Presto et al., EST, 2005
#endif  /* TRACERS_TERP */
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!298.,298., & ! ibpinp1a_nox, ibpinp2a_nox: Griffin et al., JGR, 1999
!modelE!#  endif
!modelE!298.,298., & ! ibpinp1a_hox, ibpinp2a_hox: Griffin et al., JGR, 1999
!modelE!#endif
!modelE!#ifdef SOA_FULL
!modelE!298.,298., & ! itolp1a_nox,  itolp2a_nox: Odum et al., Science, 1997
!modelE!#endif
!modelE!298.,298.  & ! itolp1a_hox,  itolp2a_hox: Odum et al., Science, 1997
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!,300.,300.  & ! ixylp1a_nox,  ixylp2a_nox: Song et al., EST, 2005
!modelE!#  endif
!modelE!,300.,300.  & ! ixylp1a_hox,  ixylp2a_hox: Song et al., EST, 2005
!modelE!#endif
!modelE!!#ifdef SOA_FULL
!modelE!!295.,295., & ! itolp1a_nox,  itolp2a_nox: Ng et al., ACP, 2007
!modelE!!#endif
!modelE!!295.,295.  & ! itolp1a_hox,  itolp2a_hox: Ng et al., ACP, 2007
!modelE!!#ifndef SOA_MINIMUM
!modelE!!#  ifdef SOA_FULL
!modelE!!,295.,295.  & ! ixylp1a_nox,  ixylp2a_nox: Ng et al., ACP, 2007
!modelE!!#  endif
!modelE!!,295.,295.  & ! ixylp1a_hox,  ixylp2a_hox: Ng et al., ACP, 2007
!modelE!!#endif
/)

!@param kpart_temp_nox_ref reference temperature where kpart_nox_ref was derived
real, dimension(nsoa), parameter :: kpart_temp_nox_ref=(/&
295.d0,295.d0   & ! iisopp1a_nox, iisopp2a_nox
#ifdef TRACERS_TERP
,295.d0,295.d0  & ! iapinp1a_nox, iapinp2a_nox: Presto et al., EST, 2005
#endif  /* TRACERS_TERP */
!modelE!298.,298., & ! ibpinp1a_nox, ibpinp2a_nox: Griffin et al., JGR, 1999
!modelE!298.,298., & ! itolp1a_nox,  itolp2a_nox: Odum et al., Science, 1997
!modelE!300.,300.  & ! ixylp1a_nox,  ixylp2a_nox: Song et al., EST, 2005
!modelE!!295.,295., & ! itolp1a_nox,  itolp2a_nox: Ng et al., ACP, 2007
!modelE!!295.,295.  & ! ixylp1a_nox,  ixylp2a_nox: Ng et al., ACP, 2007
/)

character(len=300) :: out_line


contains

subroutine soa_init

!use TRACER_COM, only: n_bcii,n_bcia,n_bcb,n_ocii,n_ocia,n_ocb,n_ococean
use CONSTANT, only: byavog
implicit none

integer   :: i,j
!@param mw_c atomic weight of carbon
real, parameter :: mw_c=12.01078d0
!@param mw_c atomic weight of hydrogen
real, parameter :: mw_h=1.007947d0
!@param mw_c atomic weight of oxygen
real, parameter :: mw_o=15.99943d0
!modelE!real*8    :: apartmolar_tot

!
! define soa species
!
issoa(1)=n_isopp1a
issoa(2)=n_isopp2a
#ifdef TRACERS_TERP
issoa(3)=n_apinp1a
issoa(4)=n_apinp2a
#endif  /* TRACERS_TERP */

allocate(mw(ntm))
allocate(molec2ug(ntm))
!
! create whichsoa from issoa, in order to correlate the two variables
!
allocate(whichsoa(ntm))
whichsoa=0
do i=1,nsoa
  do j=1,ntm
    if (issoa(i)==j) then
      whichsoa(j)=i
    endif
  enddo
enddo

! Correct tr_mm in order to use the real MW of species in meanmw and activity coefficient calculations only
mw=tr_mm()
mw(n_Isoprene)=5.d0*mw_c+8.d0*mw_h ! C5H8
#ifdef TRACERS_TERP
mw(n_Terpenes)=10.d0*mw_c+16.d0*mw_h ! C10H16
#endif  /* TRACERS_TERP */
!mw(n_bcii)=170.d0
!mw(n_bcia)=170.d0
!mw(n_bcb)=170.d0
!mw(n_ocii)=170.d0
!mw(n_ocia)=170.d0
!mw(n_ocb)=170.d0
!mw(n_ococean)=170.d0

!
! High and low-NOx mass based stoicheiometric coefficients
!
apartmass_ref=0.d0
apartmass_nox_ref=0.d0
#ifdef TRACERS_TERP
!!! WARNING !!! Change also isoprene values for the case where TRACERS_TERP are off
! a-pinene SOAb formation, Presto et al., EST, 2005
apartmass_ref(whichsoa(n_apinp1a))=0.192d0
apartmass_ref(whichsoa(n_apinp2a))=0.215d0
apartmass_nox_ref(whichsoa(n_apinp1a))=0.0138d0
apartmass_nox_ref(whichsoa(n_apinp2a))=0.461d0
#endif  /* TRACERS_TERP */
!modelE!! xylene SOAa formation, Song et al., EST, 2005 (first two lines) or Ng et al., ACP, 2007 (latter two lines)
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!apartmass_ref(whichsoa(ixylp1a_nox))=0.049
!modelE!apartmass_ref(whichsoa(ixylp2a_nox))=0.178
!modelE!!apartmass_ref(whichsoa(ixylp1a_nox))=0.031
!modelE!!apartmass_ref(whichsoa(ixylp2a_nox))=0.090
!modelE!#  endif
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmass_nox_ref(whichsoa(ixylp1a_hox))=0.049
!modelE!apartmass_nox_ref(whichsoa(ixylp2a_hox))=0.178
!modelE!!apartmass_nox_ref(whichsoa(ixylp1a_hox))=0.031
!modelE!!apartmass_nox_ref(whichsoa(ixylp2a_hox))=0.090
!modelE!#  endif
!modelE!apartmass_ref(whichsoa(ixylp1a_hox))=0.024
!modelE!apartmass_ref(whichsoa(ixylp2a_hox))=0.152
!modelE!!apartmass_ref(whichsoa(ixylp1a_hox))=0.30
!modelE!!apartmass_ref(whichsoa(ixylp2a_hox))=0.
!modelE!#endif
! isoprene + OH low NOx SOAb formation, Henze and Seinfeld, GRL 2006; scaled, based on a-pinene
! WARNING!!! Species indices are inverted compared to the paper, for consistency:
!            less volatile first, more volatile second
apartmass_ref(whichsoa(n_isopp1a))=0.0288d0
apartmass_ref(whichsoa(n_isopp2a))=0.232d0
#ifdef TRACERS_TERP
apartmass_nox_ref(whichsoa(n_isopp1a))=apartmass_ref(whichsoa(n_isopp1a))*apartmass_nox_ref(whichsoa(n_apinp1a))/&
                                                                          apartmass_ref(whichsoa(n_apinp1a))
apartmass_nox_ref(whichsoa(n_isopp2a))=apartmass_ref(whichsoa(n_isopp2a))*apartmass_nox_ref(whichsoa(n_apinp2a))/&
                                                                          apartmass_ref(whichsoa(n_apinp2a))
#else
apartmass_nox_ref(whichsoa(n_isopp1a))=apartmass_ref(whichsoa(n_isopp1a))*0.0138d0/0.192d0 !!! HARDCODED !!!
apartmass_nox_ref(whichsoa(n_isopp2a))=apartmass_ref(whichsoa(n_isopp2a))*0.461d0/0.215d0  !!! HARDCODED !!!
#endif  /* TRACERS_TERP */
!modelE!! b-pinene SOAb formation, Griffin et al., JGR, 1999; scaled, based on a-pinene
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!apartmass_ref(whichsoa(ibpinp1a_nox))=0.026/0.125*apartmass_ref(whichsoa(iapinp1a_nox))
!modelE!apartmass_ref(whichsoa(ibpinp2a_nox))=0.485/0.102*apartmass_ref(whichsoa(iapinp2a_nox))
!modelE!#  endif
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmass_nox_ref(whichsoa(ibpinp1a_hox))=0.026/0.125*apartmass_nox_ref(whichsoa(iapinp1a_hox))
!modelE!apartmass_nox_ref(whichsoa(ibpinp2a_hox))=0.485/0.102*apartmass_nox_ref(whichsoa(iapinp2a_hox))
!modelE!#  endif
!modelE!apartmass_ref(whichsoa(ibpinp1a_hox))=0.026/0.125*apartmass_ref(whichsoa(iapinp1a_hox))
!modelE!apartmass_ref(whichsoa(ibpinp2a_hox))=0.485/0.102*apartmass_ref(whichsoa(iapinp2a_hox))
!modelE!#endif
!modelE!! toluene SOAa formation, Odum et al., Science, 1997; scaled, based on xylene (first two lines) or Ng et al., ACP, 2007, no scaling (latter two lines)
!modelE!#ifdef SOA_FULL
!modelE!apartmass_ref(whichsoa(itolp1a_nox))=0.071/0.038*apartmass_ref(whichsoa(ixylp1a_nox))
!modelE!apartmass_ref(whichsoa(itolp2a_nox))=0.138/0.167*apartmass_ref(whichsoa(ixylp2a_nox))
!modelE!!apartmass_ref(whichsoa(itolp1a_nox))=0.058
!modelE!!apartmass_ref(whichsoa(itolp2a_nox))=0.113
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmass_nox_ref(whichsoa(itolp1a_hox))=0.071/0.038*0.049!apartmass_nox_ref(whichsoa(ixylp1a_hox)) ! HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!modelE!apartmass_nox_ref(whichsoa(itolp2a_hox))=0.138/0.167*0.178!apartmass_nox_ref(whichsoa(ixylp2a_hox)) ! HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!modelE!!apartmass_nox_ref(whichsoa(itolp1a_hox))=0.058
!modelE!!apartmass_nox_ref(whichsoa(itolp2a_hox))=0.113
!modelE!#endif
!modelE!apartmass_ref(whichsoa(itolp1a_hox))=0.071/0.038*0.024!apartmass_ref(whichsoa(ixylp1a_hox)) ! HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!modelE!apartmass_ref(whichsoa(itolp2a_hox))=0.138/0.167*0.152!apartmass_ref(whichsoa(ixylp2a_hox)) ! HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!modelE!!apartmass_ref(whichsoa(itolp1a_hox))=0.36
!modelE!!apartmass_ref(whichsoa(itolp2a_hox))=0.

do i=1,soacomp
  do j=1,soacomp
    if (Lambda(i,j).eq.-1.0d0) Lambda(i,j)=Lambda(j,i) ! symmetric interactions
  enddo
enddo

do i=1,ntm
  molec2ug(i)=tr_mm(i)*1.d12*byavog
enddo

if (am_i_root()) then
  out_line="Initialization of SOA formation completed"
  call write_parallel(trim(out_line),crit=.true.)
endif

end subroutine soa_init


#ifdef SOA_DIAGS
subroutine soa_apart(III,JJJ)
#else
subroutine soa_apart
#endif  /* SOA_DIAGS */

#ifdef SOA_DIAGS
use TRDIAG_COM, only: taijls=>taijls_loc,&
                      ijlt_soa_voc2nox,ijlt_soa_apartmass
#endif  /* SOA_DIAGS */
implicit none

#ifdef SOA_DIAGS
integer, intent(in) :: III,JJJ
#endif  /* SOA_DIAGS */
integer             :: jl,i

apartmass=0.d0
apartmolar=0.d0

!
! mass based stoicheiometric coefficients
!
do jl=1,LM
#ifdef SOA_DIAGS
  taijls(III,JJJ,jl,ijlt_soa_voc2nox)=taijls(III,JJJ,jl,ijlt_soa_voc2nox)+voc2nox(jl)
#endif  /* SOA_DIAGS */
  do i=1,nsoa
    apartmass(jl,i)=voc2nox(jl)*apartmass_nox_ref(i)+(1.d0-voc2nox(jl))*apartmass_ref(i)
#ifdef SOA_DIAGS
    taijls(III,JJJ,jl,ijlt_soa_apartmass(i))=taijls(III,JJJ,jl,ijlt_soa_apartmass(i))+apartmass(jl,i)
#endif  /* SOA_DIAGS */
  enddo
enddo

!
! molar based stoicheiometric coefficients
!
apartmolar(:,whichsoa(n_isopp1a))=apartmass(:,whichsoa(n_isopp1a))*mw(n_Isoprene)/tr_mm(n_isopp1a)
apartmolar(:,whichsoa(n_isopp2a))=apartmass(:,whichsoa(n_isopp2a))*mw(n_Isoprene)/tr_mm(n_isopp2a)
#ifdef TRACERS_TERP
apartmolar(:,whichsoa(n_apinp1a))=apartmass(:,whichsoa(n_apinp1a))*mw(n_Terpenes)/tr_mm(n_apinp1a)
apartmolar(:,whichsoa(n_apinp2a))=apartmass(:,whichsoa(n_apinp2a))*mw(n_Terpenes)/tr_mm(n_apinp2a)
#endif  /* TRACERS_TERP */

end subroutine soa_apart


subroutine soa_aerosolphase(III,JJJ,L,changeL,bypfactor)
use OldTracer_mod, only: mass2vol
use TRACER_COM, only: trm,n_bcii,n_bcia,n_bcb,n_ocii,n_ocia,n_ocb,n_ococean,&
#ifdef TRACERS_NITRATE
                      n_nh4,n_no3p,&
#endif
                      n_msa,n_so4
#ifdef TRACERS_AEROSOLS_VBS
use TRACERS_VBS, only: vbs_tr
#endif
#ifdef SOA_DIAGS
use TRDIAG_COM, only: taijls=>taijls_loc,&
                      ijlt_soa_changeL_isoprene,&
#ifdef TRACERS_TERP
                      ijlt_soa_changeL_terpenes,&
#endif  /* TRACERS_TERP */
                      ijlt_soa_pcp,ijlt_soa_aerotot,ijlt_soa_aerotot_gas,&
                      ijlt_soa_xmf_isop,ijlt_soa_xmf_apin,&
                      ijlt_soa_zcoef_isop,ijlt_soa_zcoef_apin,ijlt_soa_meanmw,&
                      ijlt_soa_iternum,ijlt_soa_m0,&
                      ijlt_soa_y0_ug_g,ijlt_soa_y0_ug_a,&
                      ijlt_soa_y_ug_g,ijlt_soa_y_ug_a,&
                      ijlt_soa_changeL_g_before,ijlt_soa_changeL_a_before,&
                      ijlt_soa_changeL_g_after,ijlt_soa_changeL_a_after,&
                      ijlt_soa_kpart,ijlt_soa_kp,ijlt_soa_soamass,ijlt_soa_partfact,&
                      ijlt_soa_evap,ijlt_soa_cond,ijlt_soa_chem
#endif  /* SOA_DIAGS */
implicit none

!@var y0_ug concentration of species (in molecules/cm3) before chemistry and partitioning
!@var y_ug concentration of species (in molecules/cm3) after chemistry and partitioning
!@var y_mw =y/tr_mm*mw
real*8,dimension(ntm)                    :: y0_ug,y_ug,y_mw
real*8, intent(in)                       :: bypfactor
real*8, intent(inout) :: changeL(:,:) ! automatic
!@var jl looping index
integer                                  :: jl
integer, intent(in)                      :: III,JJJ,L
!@var iterbackward logical that imposes the iterative solution direction
logical                                  :: iterbackward
!
! For correction of Kp due to composition
!
!@var AEROtot the total aerosol concentration, which includes both the semivolatile aerosol phase and
!@+           the non-volatile aerosol phase that affects condensation, BEFORE new condensation
!@var AEROtot_gas same as AEROtot but for the gas phase (used as a replacement when AEROtot=0.d0)
!@var sigma-ij,sigma_ik,sigma_jk used for the activity coefficient calculations
!@var zc the activity coefficient BEFORE new condensation and BEFORE the validity check (needed only for diagnistic output)
!@var meanmw the mean molecular weight of the total condensation affected aerosol phase, BEFORE new condensation
real*8                      :: AEROtot,AEROtot_gas,sigma_ij,sigma_ik,sigma_kj,zc,meanmw
!@var imf,jmf,kmf indices for the activity coefficient calculation
integer                     :: imf,jmf,kmf
!@var xmf the mole fraction of species in the aerosol phase BEFORE new condensation
!@var zcoef final value of activity coefficient BEFORE new condensation (this is the value that is
!@+         going to be used). The allowed range is 0.3 < zcoef < 5.0
real*8, dimension(soacomp)  :: xmf,zcoef
!@var kp final value of partitioning coefficient after applying all relevant parameters to the reference case kpart_ref
real*8, dimension(nsoa)     :: kp
! Other SOA parameters
!@var M0err maximum error requested by the iterative solution. THIS IS COMPUTER ARCHITECTURE DEPENDENT!
!@+   For single precision computers, 1.e-6 is the lowest it can go. For double precision this number
!@+   should be able to go down to 1.d-12, but it is not tested. 1.d-10 is more than sufficient.
real*8, parameter           :: M0err=1.d-10
!@var M0 the total aerosol concentration. After the end of the iterations, this is the solution.
!@var PCP the total aerosol concentration that cannot evaporate. The abbreviation comes from
!@+       Primary Carbonaceous Particles but it can also contain inorganic species or non-volatile SOA
!@var M0temp the M0 before an iterative solution is found
!@var M0a lower limit of M0 during iteration
!@var M0b upper limit of M0 during iteration
!@var M0err_curr the M0 error for the current iteration, based on M0a, M0b, M0err and computer architecture
real*8                      :: M0,PCP,M0temp,M0a,M0b,M0err_curr
!@var soamass the total concentration that can partition (ug m-3)
!@var partfact,partfact1,partfact2 help variable for the M0 calculation
real*8, dimension(nsoa)     :: soamass,partfact,partfact1,partfact2
!@var iternum counter for the number of iterations per box (not saved)
integer                     :: i,iternum
!@var x1,x2,y1,y2,a,b iteration help parameter
real*8                      :: x1,x2,y1,y2,a,b

!@var SO4part set to true if partitioning can occur in sulfur-containing aerosols
logical, parameter :: SO4part=.false.
#ifdef TRACERS_NITRATE
!@var NH4part set to true if partitioning can occur in nitrogen-containing aerosols
logical, parameter :: NH4part=.false.
#endif
!@var SOAevap set to true if SOA can evaporate after condensation
logical, parameter :: SOAevap=.true.

!DO JL=1,LM ! if this gets put back in, use size( ) to get changeL L-dim instead of LM!
DO JL=L,L

!
! Change concentrations to ug/m3
! Note that we want the concentrations AFTER chemistry, 
! thus we apply changeL artificially (but not hardcoded)
!
#ifdef SOA_DIAGS
  taijls(III,JJJ,jl,ijlt_soa_changeL_isoprene)=taijls(III,JJJ,jl,ijlt_soa_changeL_isoprene)+&
                                               changeL(jl,n_isoprene)*bypfactor*mass2vol(n_isoprene)*molec2ug(n_isoprene)
#ifdef TRACERS_TERP
  taijls(III,JJJ,jl,ijlt_soa_changeL_terpenes)=taijls(III,JJJ,jl,ijlt_soa_changeL_terpenes)+&
                                               changeL(jl,n_terpenes)*bypfactor*mass2vol(n_terpenes)*molec2ug(n_terpenes)
#endif  /* TRACERS_TERP */
#endif  /* SOA_DIAGS */
  do i=1,ntm
    y0_ug(i)=trm(III,JJJ,jl,i)*bypfactor*mass2vol(i)*molec2ug(i)
    y_ug(i)=(trm(III,JJJ,jl,i)+changeL(jl,i))*bypfactor*mass2vol(i)*molec2ug(i)
  enddo
#ifdef SOA_DIAGS
  do i=1,nsoa
    taijls(III,JJJ,jl,ijlt_soa_y0_ug_g(i))=taijls(III,JJJ,jl,ijlt_soa_y0_ug_g(i))+y0_ug(issoa(i)-1)
    taijls(III,JJJ,jl,ijlt_soa_y_ug_g(i))=taijls(III,JJJ,jl,ijlt_soa_y_ug_g(i))+y_ug(issoa(i)-1)
    taijls(III,JJJ,jl,ijlt_soa_changeL_g_before(i))=taijls(III,JJJ,jl,ijlt_soa_changeL_g_before(i))+&
                                                    changeL(jl,issoa(i)-1)*bypfactor*mass2vol(issoa(i)-1)*molec2ug(issoa(i)-1)
    taijls(III,JJJ,jl,ijlt_soa_y0_ug_a(i))=taijls(III,JJJ,jl,ijlt_soa_y0_ug_a(i))+y0_ug(issoa(i))
    taijls(III,JJJ,jl,ijlt_soa_y_ug_a(i))=taijls(III,JJJ,jl,ijlt_soa_y_ug_a(i))+y_ug(issoa(i))
    taijls(III,JJJ,jl,ijlt_soa_changeL_a_before(i))=taijls(III,JJJ,jl,ijlt_soa_changeL_a_before(i))+&
                                                    changeL(jl,issoa(i))*bypfactor*mass2vol(issoa(i))*molec2ug(issoa(i))
  enddo
#endif  /* SOA_DIAGS */

!
! Calculate the primary aerosol able to absorb semivolatile compounds
! Partitioning in both OC and BC always happens
! Partitioning in (SO4+MSA) and/or NH4/NO3 is also an option
!
  PCP=y_ug(n_bcii)+y_ug(n_bcia)+y_ug(n_bcb)
#ifdef TRACERS_AEROSOLS_VBS
  PCP=PCP+sum(vbs_tr%aer(:))
#else
  PCP=PCP+y_ug(n_ocii)+y_ug(n_ocia)+y_ug(n_ocb)
#endif /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_OCEAN
  PCP=PCP+y_ug(n_ococean)
#endif  /* TRACERS_AEROSOLS_OCEAN */
  if(SO4part) PCP=PCP+y_ug(n_msa)+y_ug(n_so4)
#ifdef TRACERS_NITRATE
  if(NH4part) PCP=PCP+y_ug(n_nh4)+y_ug(n_no3p)
#endif
#ifdef SOA_DIAGS
  taijls(III,JJJ,jl,ijlt_soa_pcp)=taijls(III,JJJ,jl,ijlt_soa_pcp)+PCP
#endif  /* SOA_DIAGS */
!
! Correct the Kp to take into account the change of activity coefficient due to change in composition
! using the Wilson equation. Method described in Bowman and Karamalegos, EST, 2002, 36, 2701-2707.
!
  do i=1,ntm
    y_mw(i)=y_ug(i)/mw(i)
  enddo
  AEROtot=y_mw(n_bcii)+y_mw(n_bcia)+y_mw(n_bcb)
#ifdef TRACERS_AEROSOLS_VBS
  AEROtot=AEROtot+sum(y_mw(vbs_tr%iaer))
#else
  AEROtot=AEROtot+y_mw(n_ocii)+y_mw(n_ocia)+y_mw(n_ocb)
#endif /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_OCEAN
  AEROtot=AEROtot+y_mw(n_ococean)
#endif  /* TRACERS_AEROSOLS_OCEAN */
  if(SO4part) AEROtot=AEROtot+y_mw(n_msa)+y_mw(n_so4)
#ifdef TRACERS_NITRATE
  if(NH4part) AEROtot=AEROtot+y_mw(n_nh4)+y_mw(n_no3p)
#endif
  do i=1,nsoa
    AEROtot=AEROtot+y_mw(issoa(i))
  enddo
#ifdef SOA_DIAGS
  taijls(III,JJJ,jl,ijlt_soa_aerotot)=taijls(III,JJJ,jl,ijlt_soa_aerotot)+AEROtot
#endif  /* SOA_DIAGS */
  if (AEROtot <= 0.d0) then ! also grabs cases with negative aerosol concentrations, which should never occur
    AEROtot_gas=0.d0
    do i=1,nsoa
      AEROtot_gas=AEROtot_gas+y_mw(issoa(i)-1)
    enddo
    if (AEROtot_gas <= 0.d0) goto 60 ! Neither aerosols nor semivolatile gases exist, do nothing.
#ifdef SOA_DIAGS
    taijls(III,JJJ,jl,ijlt_soa_aerotot_gas)=taijls(III,JJJ,jl,ijlt_soa_aerotot_gas)+AEROtot_gas
#endif  /* SOA_DIAGS */
  endif
!ktt  if (AEROtot.le.0.d0) goto 60 ! This part will be important on nucleation in the future versions
!
! Calculate mole fraction of individual species.
!
  xmf=0.d0
  if (AEROtot > 0.d0) then
    xmf(imfisop)=(y_mw(n_isopp1a)+y_mw(n_isopp2a))/AEROtot
#ifdef TRACERS_TERP
    xmf(imfapin)=(y_mw(n_apinp1a)+y_mw(n_apinp2a))/AEROtot
#endif  /* TRACERS_TERP */
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!              +y(jl,ibpinp1a_nox)/tr_mm(ibpinp1a_nox)*mw(ibpinp1a_nox)+y(jl,ibpinp2a_nox)/tr_mm(ibpinp2a_nox)*mw(ibpinp2a_nox)&
!modelE!#  endif
!modelE!              +y(jl,ibpinp1a_hox)/tr_mm(ibpinp1a_hox)*mw(ibpinp1a_hox)+y(jl,ibpinp2a_hox)/tr_mm(ibpinp2a_hox)*mw(ibpinp2a_hox)&
!modelE!#endif
!modelE!              )/AEROtot
!modelE!  xmf(imfaro)=(&
!modelE!#ifdef SOA_FULL
!modelE!               y(jl,itolp1a_nox)/tr_mm(itolp1a_nox)*mw(itolp1a_nox)+y(jl,itolp2a_nox)/tr_mm(itolp2a_nox)*mw(itolp2a_nox)&
!modelE!#endif
!modelE!              +y(jl,itolp1a_hox)/tr_mm(itolp1a_hox)*mw(itolp1a_hox)+y(jl,itolp2a_hox)/tr_mm(itolp2a_hox)*mw(itolp2a_hox)&
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!              +y(jl,ixylp1a_nox)/tr_mm(ixylp1a_nox)*mw(ixylp1a_nox)+y(jl,ixylp2a_nox)/tr_mm(ixylp2a_nox)*mw(ixylp2a_nox)&
!modelE!#  endif
!modelE!              +y(jl,ixylp1a_hox)/tr_mm(ixylp1a_hox)*mw(ixylp1a_hox)+y(jl,ixylp2a_hox)/tr_mm(ixylp2a_hox)*mw(ixylp2a_hox)&
!modelE!#endif
!modelE!              )/AEROtot
    xmf(imfbcii)=y_mw(n_bcii)/AEROtot
    xmf(imfbcia)=y_mw(n_bcia)/AEROtot
    xmf(imfbcb)=y_mw(n_bcb)/AEROtot
#ifdef TRACERS_AEROSOLS_VBS
    xmf(imfocii)=0.d0
    xmf(imfocia)=sum(y_mw(vbs_tr%iaer))/AEROtot
    xmf(imfocb)=0.d0
#else
    xmf(imfocii)=y_mw(n_ocii)/AEROtot
    xmf(imfocia)=y_mw(n_ocia)/AEROtot
    xmf(imfocb)=y_mw(n_ocb)/AEROtot
#endif /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_OCEAN
    xmf(imfococean)=y_mw(n_ococean)/AEROtot
#endif  /* TRACERS_AEROSOLS_OCEAN */
!
! If partitioning does not occur on some aerosol species, it's mole fraction equals to zero, no matter the real
! concentration is, in order not to affect the calculation of the activity coefficient.
!
    if(SO4part) then
      xmf(imfinorg)=xmf(imfinorg)+y_mw(n_msa)/AEROtot
      xmf(imfinorg)=xmf(imfinorg)+y_mw(n_so4)/AEROtot
    endif
#ifdef TRACERS_NITRATE
    if(NH4part) then
      xmf(imfinorg)=xmf(imfinorg)+y_mw(n_nh4)/AEROtot
      xmf(imfinorg)=xmf(imfinorg)+y_mw(n_no3p)/AEROtot
    endif
#endif
  else
    xmf(imfisop)=(y_mw(n_isopp1g)+y_mw(n_isopp2g))/AEROtot_gas
#ifdef TRACERS_TERP
    xmf(imfapin)=(y_mw(n_apinp1g)+y_mw(n_apinp2g))/AEROtot_gas
#endif  /* TRACERS_TERP */
  endif ! AEROtot > 0.d0
#ifdef SOA_DIAGS
  taijls(III,JJJ,jl,ijlt_soa_xmf_isop)=taijls(III,JJJ,jl,ijlt_soa_xmf_isop)+xmf(imfisop)
  taijls(III,JJJ,jl,ijlt_soa_xmf_apin)=taijls(III,JJJ,jl,ijlt_soa_xmf_apin)+xmf(imfapin)
#endif  /* SOA_DIAGS */
!
! Calculate activity coefficient zcoef
!
  zcoef=1.d0
!  do kmf=1,soacomp
  do kmf=imfapin,imfisop
    sigma_kj=0.d0
    do jmf=1,soacomp
      sigma_kj=sigma_kj+xmf(jmf)*Lambda(kmf,jmf)
    enddo
    sigma_ik=0.d0
    do imf=1,soacomp
      sigma_ij=0.d0
      do jmf=1,soacomp
        sigma_ij=sigma_ij+xmf(jmf)*Lambda(imf,jmf)
      enddo
      sigma_ik=sigma_ik+xmf(imf)*Lambda(imf,kmf)/sigma_ij
    enddo
    zc=exp(1.d0-log(sigma_kj)-sigma_ik)
    zcoef(kmf)=max(0.3d0,min(5.0d0,zc))
    if((zc.lt.0.3d0).or.(zc.gt.5.0d0)) then
      write(out_line,*)'WARNING: zcoef set to',zcoef(kmf),' (was',zc,')'
      call write_parallel(trim(out_line),crit=.true.)
    endif
  enddo
#ifdef SOA_DIAGS
  taijls(III,JJJ,jl,ijlt_soa_zcoef_isop)=taijls(III,JJJ,jl,ijlt_soa_zcoef_isop)+zcoef(imfisop)
  taijls(III,JJJ,jl,ijlt_soa_zcoef_apin)=taijls(III,JJJ,jl,ijlt_soa_zcoef_apin)+zcoef(imfapin)
#endif  /* SOA_DIAGS */
!
! Calculate mean molecular weight
!
  if (AEROtot > 0.d0) then
    meanmw=0.d0
#ifdef TRACERS_AEROSOLS_VBS
    do i=1,vbs_tr%nbins
      meanmw=meanmw+y_mw(vbs_tr%iaer(i))*mw(vbs_tr%iaer(i))/AEROtot
    enddo
#else
    meanmw=meanmw+y_mw(n_ocii)*mw(n_ocii)/AEROtot
    meanmw=meanmw+y_mw(n_ocia)*mw(n_ocia)/AEROtot
    meanmw=meanmw+y_mw(n_ocb)*mw(n_ocb)/AEROtot
#endif /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_OCEAN
    meanmw=meanmw+y_mw(n_ococean)*mw(n_ococean)/AEROtot
#endif  /* TRACERS_AEROSOLS_OCEAN */
    meanmw=meanmw+y_mw(n_bcii)*mw(n_bcii)/AEROtot
    meanmw=meanmw+y_mw(n_bcia)*mw(n_bcia)/AEROtot
    meanmw=meanmw+y_mw(n_bcb)*mw(n_bcb)/AEROtot
    if(SO4part) then
      meanmw=meanmw+y_mw(n_msa)*mw(n_msa)/AEROtot
      meanmw=meanmw+y_mw(n_so4)*mw(n_so4)/AEROtot
    endif
#ifdef TRACERS_NITRATE
    if(NH4part) then
      meanmw=meanmw+y_mw(n_nh4)*mw(n_nh4)/AEROtot
      meanmw=meanmw+y_mw(n_no3p)*mw(n_no3p)/AEROtot
    endif
#endif
    do i=1,nsoa
      meanmw=meanmw+y_mw(issoa(i))*mw(issoa(i))/AEROtot
    enddo
  else
    do i=1,nsoa
! This is based on the gas-phase concentrations but with the aerosol phase mw (due to the mw array size), since AEROtot=0.d0
      meanmw=meanmw+y_mw(issoa(i)-1)*mw(issoa(i))/AEROtot_gas
    enddo
  endif ! AEROtot > 0.d0
#ifdef SOA_DIAGS
  taijls(III,JJJ,jl,ijlt_soa_meanmw)=taijls(III,JJJ,jl,ijlt_soa_meanmw)+meanmw
#endif  /* SOA_DIAGS */
!
! Calculate final value of partitioning coefficient
!
  kp(whichsoa(n_isopp1a))=kpart(jl,whichsoa(n_isopp1a))/zcoef(imfisop)*mw(n_isopp1a)/meanmw
  kp(whichsoa(n_isopp2a))=kpart(jl,whichsoa(n_isopp2a))/zcoef(imfisop)*mw(n_isopp2a)/meanmw
#ifdef TRACERS_TERP
  kp(whichsoa(n_apinp1a))=kpart(jl,whichsoa(n_apinp1a))/zcoef(imfapin)*mw(n_apinp1a)/meanmw
  kp(whichsoa(n_apinp2a))=kpart(jl,whichsoa(n_apinp2a))/zcoef(imfapin)*mw(n_apinp2a)/meanmw
#endif  /* TRACERS_TERP */
#ifdef SOA_DIAGS
  do i=1,nsoa
    taijls(III,JJJ,jl,ijlt_soa_kpart(i))=taijls(III,JJJ,jl,ijlt_soa_kpart(i))+kpart(jl,i)
    taijls(III,JJJ,jl,ijlt_soa_kp(i))=taijls(III,JJJ,jl,ijlt_soa_kp(i))+kp(i)
  enddo
#endif  /* SOA_DIAGS */
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!  kp(whichsoa(ibpinp1a_nox))=kpart(jl,whichsoa(ibpinp1a_nox))/zcoef(imfapin)*mw(ibpinp1a_nox)/meanmw
!modelE!  kp(whichsoa(ibpinp2a_nox))=kpart(jl,whichsoa(ibpinp2a_nox))/zcoef(imfapin)*mw(ibpinp2a_nox)/meanmw
!modelE!#  endif
!modelE!  kp(whichsoa(ibpinp1a_hox))=kpart(jl,whichsoa(ibpinp1a_hox))/zcoef(imfapin)*mw(ibpinp1a_hox)/meanmw
!modelE!  kp(whichsoa(ibpinp2a_hox))=kpart(jl,whichsoa(ibpinp2a_hox))/zcoef(imfapin)*mw(ibpinp2a_hox)/meanmw
!modelE!#endif
!modelE!#ifdef SOA_FULL
!modelE!  kp(whichsoa(itolp1a_nox))=kpart(jl,whichsoa(itolp1a_nox))/zcoef(imfaro)*mw(itolp1a_nox)/meanmw
!modelE!  kp(whichsoa(itolp2a_nox))=kpart(jl,whichsoa(itolp2a_nox))/zcoef(imfaro)*mw(itolp2a_nox)/meanmw
!modelE!#endif
!modelE!  kp(whichsoa(itolp1a_hox))=kpart(jl,whichsoa(itolp1a_hox))/zcoef(imfaro)*mw(itolp1a_hox)/meanmw
!modelE!  kp(whichsoa(itolp2a_hox))=kpart(jl,whichsoa(itolp2a_hox))/zcoef(imfaro)*mw(itolp2a_hox)/meanmw
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!  kp(whichsoa(ixylp1a_nox))=kpart(jl,whichsoa(ixylp1a_nox))/zcoef(imfaro)*mw(ixylp1a_nox)/meanmw
!modelE!  kp(whichsoa(ixylp2a_nox))=kpart(jl,whichsoa(ixylp2a_nox))/zcoef(imfaro)*mw(ixylp2a_nox)/meanmw
!modelE!#  endif
!modelE!  kp(whichsoa(ixylp1a_hox))=kpart(jl,whichsoa(ixylp1a_hox))/zcoef(imfaro)*mw(ixylp1a_hox)/meanmw
!modelE!  kp(whichsoa(ixylp2a_hox))=kpart(jl,whichsoa(ixylp2a_hox))/zcoef(imfaro)*mw(ixylp2a_hox)/meanmw
!modelE!#endif
!
! Add soa species gas and particulate phase in variable soamass (ug m-3)
! If no evaporation, add only gas phase
!
  do i=1,nsoa
    soamass(i)=y_ug(issoa(i)-1)
    if (SOAevap) soamass(i)=soamass(i)+y_ug(issoa(i))
#ifdef SOA_DIAGS
    taijls(III,JJJ,jl,ijlt_soa_soamass(i))=taijls(III,JJJ,jl,ijlt_soa_soamass(i))+soamass(i)
#endif  /* SOA_DIAGS */
  enddo
  if (sum(soamass)==0.d0) goto 60 ! no semivolatiles, nothing to do
  M0a=PCP
  M0b=PCP
  do i=1,nsoa
    M0b=M0b+soamass(i)
    if(.not.SOAevap) then
      M0a=M0a+y_ug(issoa(i))
      M0b=M0b+y_ug(issoa(i))
    endif
  enddo
  iternum=0
!
! Iteration - chord method
!
  M0=M0a
  iterbackward=.false.
11  continue ! Try with another M0
  if (iternum==0) then
    x2=M0
    M0err_curr=max(tiny(x2),x2*M0err)
    if (x2 <= 0.d0) M0err_curr=max(tiny(x2),(x2+sum(soamass))*M0err)
    x1=x2-M0err_curr
  else if (iterbackward) then
    x2=M0
    M0err_curr=max(tiny(x2),x2*M0err)
    if (x2 <= 0.d0) M0err_curr=max(tiny(x2),(x2+sum(soamass))*M0err)
    x1=x2-M0err_curr
  else
    x1=M0
    M0err_curr=max(tiny(x1),x1*M0err)
    x2=x1+M0err_curr
  endif
  y1=PCP-x1
  y2=PCP-x2
  do i=1,nsoa
    partfact2(i)=kp(i)*x2/(1.d0+kp(i)*x2)
    y2=y2+soamass(i)*partfact2(i)
    partfact1(i)=kp(i)*x1/(1.d0+kp(i)*x1)
    y1=y1+soamass(i)*partfact1(i)
  enddo
  if(abs(y1).lt.M0err) then
    if (abs(y2) < abs(y1)) then
      partfact(:)=partfact2(:)
      M0=x2
    else
      partfact(:)=partfact1(:)
      M0=x1
    endif
    goto 20
  endif
  if(abs(y2).lt.M0err) then
    partfact(:)=partfact2(:)
    M0=x2
    goto 20
  endif
!
! Stop iteration after 100 iterations (no solution found)
! This number should typically be below 5
!
  iternum=iternum+1
  if(iternum.ge.100) goto 30
!
! If -M0err/2<M0temp<M0err/2, M0 is the solution (goto 20)
! Otherwise change the limits of M0a and M0b and restart (goto 11)
!
  a=(y2-y1)/(x2-x1)
  if (a.ge.0.d0) then
    M0=M0+sum(soamass)*0.1d0
    goto 11
  endif
  b=y1-a*x1
  M0=-b/a
  if (((iternum==1).and.(M0==x2)).or.((iternum>1).and.(M0==x1))) goto 20
  if (M0<M0a .or. M0>M0b) then
    iterbackward=.true.
    M0=M0b
  endif
  goto 11
!
! No solution found, aerosol phase set to previous values, chemistry still applies to gas-phase species
!
30 continue
  write(out_line,"('WARNING: Too many iteration steps. SOA forced to previous values. M0=',1pe9.2,' PCP=',1pe9.2)") M0,PCP
  call write_parallel(trim(out_line),crit=.true.)
  ! During some testing, when the model got to this section of code (too many
  ! interations), NCCS compute nodes started to run out of memory, so there
  ! is a model stop for now:
  call stop_model('soa_aerosolphase: stop to prevent memory issue',255)
  goto 60
!
! Found solution for M0
!
20 continue
#ifdef SOA_DIAGS
  taijls(III,JJJ,jl,ijlt_soa_iternum)=taijls(III,JJJ,jl,ijlt_soa_iternum)+float(iternum)
  do i=1,nsoa
    taijls(III,JJJ,jl,ijlt_soa_partfact(i))=taijls(III,JJJ,jl,ijlt_soa_partfact(i))+partfact(i)
  enddo
  taijls(III,JJJ,jl,ijlt_soa_m0)=taijls(III,JJJ,jl,ijlt_soa_m0)+M0
#endif  /* SOA_DIAGS */
  if (M0.ne.0.d0) then
    do i=1,nsoa
      if(SOAevap) then
        y_ug(issoa(i))=soamass(i)*partfact(i)      ! Calculate new aerosol-phase concentration
        y_ug(issoa(i)-1)=y_ug(issoa(i))/(kp(i)*M0) ! Calculate new gas-phase concentration
      else
        y_ug(issoa(i))=y_ug(issoa(i))+soamass(i)*partfact(i)
        y_ug(issoa(i)-1)=max(0.d0,y_ug(issoa(i)-1)-soamass(i)*partfact(i))
      endif
#ifdef SOA_DIAGS
    taijls(III,JJJ,jl,ijlt_soa_evap(i))=taijls(III,JJJ,jl,ijlt_soa_evap(i))+&
                                        y0_ug(issoa(i))*(1.d0-partfact(i))
    taijls(III,JJJ,jl,ijlt_soa_cond(i))=taijls(III,JJJ,jl,ijlt_soa_cond(i))+&
                                        y0_ug(issoa(i)-1)*partfact(i)
    taijls(III,JJJ,jl,ijlt_soa_chem(i))=taijls(III,JJJ,jl,ijlt_soa_chem(i))+&
                                        (y_ug(issoa(i)-1)-y0_ug(issoa(i)-1))*partfact(i)
#endif  /* SOA_DIAGS */
    enddo
  endif

60 continue ! everytime aerosols are zero, goto 60 is called
!
! Save results to changeL (kg)
!
  do i=1,nsoa
    changeL(jl,issoa(i)-1)=(y_ug(issoa(i)-1)-y0_ug(issoa(i)-1))/&
                           (molec2ug(issoa(i)-1)*bypfactor*mass2vol(issoa(i)-1))
    changeL(jl,issoa(i))=(y_ug(issoa(i))-y0_ug(issoa(i)))/&
                         (molec2ug(issoa(i))*bypfactor*mass2vol(issoa(i)))
#ifdef SOA_DIAGS
    taijls(III,JJJ,jl,ijlt_soa_changeL_g_after(i))=taijls(III,JJJ,jl,ijlt_soa_changeL_g_after(i))+&
                                                   changeL(jl,issoa(i)-1)*bypfactor*mass2vol(issoa(i)-1)*molec2ug(issoa(i)-1)
    taijls(III,JJJ,jl,ijlt_soa_changeL_a_after(i))=taijls(III,JJJ,jl,ijlt_soa_changeL_a_after(i))+&
                                                   changeL(jl,issoa(i))*bypfactor*mass2vol(issoa(i))*molec2ug(issoa(i))
#endif  /* SOA_DIAGS */
  enddo

ENDDO !JL

end subroutine soa_aerosolphase


real*8 function KpCALC(dH,Ksc,temp,Tsc)
!-------------------------------------------------------------------------------
!KpCALC calculates the temperature dependence of the partitioning coefficients
!-------------------------------------------------------------------------------
use CONSTANT, only: gasc
implicit none
!@var dH enthalpy of vaporization
real*8, intent(in):: dH         ! KJ/mol
!@var Ksc reference partitioning coefficient (from champber experiments at temperature Tsc)
!@var temp temperature that the partitioning coefficient will be calculated
!@var Tsc reference temperature (from champber experiments with partitioning coefficient Ksc)
real*8, intent(in):: Ksc,temp,Tsc

KpCALC=Ksc*(temp/Tsc)*exp(dH*1.d3/gasc*(1.d0/temp-1.d0/Tsc))

end function KpCALC


end module TRACERS_SOA
