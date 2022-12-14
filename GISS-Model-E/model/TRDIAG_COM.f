#include "rundeck_opts.h"

      MODULE TRDIAG_COM
!@sum Tracer diagnostic arrays
!@+    Mostly tracer independent, but this may depend on applications
!@auth Jean Lerner
!ver   1.0
      USE RESOLUTION, only: im,jm,lm
      USE DIAG_COM, only: npts !npts are conservation quantities
     &     ,jm_budg
      USE MDIAG_COM, only : sname_strlen,units_strlen,lname_strlen
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: NTM
#ifdef TRACERS_AEROSOLS_SOA
     &                     ,nsoa
#endif  /* TRACERS_AEROSOLS_SOA */
      use Tracer_mod, only: ntsurfsrcmax, nt3Dsrcmax
#ifdef TRACERS_AMP
      USE AERO_CONFIG,only: nbins
#endif
#endif   /* TRACERS_ON or OCEAN */
      use cdl_mod
      use OldTracer_mod, only: to_conc, set_to_conc
      use OldTracer_mod, only: to_volume_MixRat, set_to_volume_MixRat
      IMPLICIT NONE
      SAVE

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>

#ifdef TRACERS_ON
#endif  /* TRACERS_ON */
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@dbparam to_per_mil For printout of tracer concentration in permil
      INTEGER, ALLOCATABLE, DIMENSION(:) :: to_per_mil
#endif
#ifdef TRACERS_ON
!@var MMR_to_VMR: converts tracer mass mixing ratio to volume mr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: MMR_to_VMR

!!! WARNING: if new diagnostics are added, keep io_trdiag up-to-date !!!
C**** TAIJLN
!@var TAIJLN 3D tracer diagnostics (all tracers)
      real*8, allocatable, dimension(:,:,:,:) :: taijln
      real*8, allocatable, dimension(:,:,:,:) :: taijln_loc
!@var TSCF3D 3D tracer advection flxes (all tracers)
      real*8, allocatable, dimension(:,:,:,:) :: TSCF3D
      real*8, allocatable, dimension(:,:,:,:) :: TSCF3D_loc
!@var SNAME_IJT, UNITS_IJT: Names and units of lat-sigma tracer IJ diags
      character(len=sname_strlen), allocatable,dimension(:) :: sname_ijt
      character(len=units_strlen), allocatable,dimension(:) :: units_ijt
!@var LNAME_IJT: descriptions of tracer IJ diagnostics
      character(len=lname_strlen), allocatable, dimension(:) ::lname_ijt
!@var SCALE_IJT: printout scaling factor for tracer IJ diagnostics
      REAL*8, allocatable, dimension(:) :: scale_ijt
!@var IJTC_POWER: power of 10 used for tracer IJ concentration diags
      integer, allocatable, dimension(:) :: ijtc_power
!@var IJTM_POWER: power of 10 used for tracer IJ mass unit diags
      integer, allocatable, dimension(:) :: ijtm_power

C**** TAIJN
!@param KTAIJ number of 2D diags describing surface and column load along 
!@+   with wet and dry deposition
!@+   please just increase this if needed - do not bother with pp options
      integer, parameter :: ktaij=23

!@var IJT_XX names for taijn diagnostics
      integer tij_conc,tij_surf,tij_surfbv,tij_mass,tij_strop
!@var IJT_XX names for water-based taijn diagnostics
      integer tij_rvr,tij_prec,tij_evap,tij_grnd,tij_lk1
     *     ,tij_lk2,tij_soil,tij_snow,tij_uflx,tij_vflx
     *     ,tij_rvro,tij_icb
      integer tij_drydep,tij_gsdep ! TRACERS_DRYDEP

#ifdef TRACERS_SPECIAL_O18
      integer tij_owiso  !Sea surface water isotope ratios
#endif

!@var TAIJN lat/lon tracer diagnostics (all tracers)
      real*8, allocatable, dimension(:,:,:,:)      :: taijn
      real*8, allocatable, dimension(:,:,:,:) :: taijn_loc

!@var SCALE_TIJ: printout scaling factor for tracer IJK diagnostics
      REAL*8, allocatable :: scale_tij(:,:)
!@var SNAME_TIJ,UNITS_TIJ: Names and units of lat-sigma tracer diags
      character(len=sname_strlen), allocatable :: sname_tij(:,:)
      character(len=units_strlen), allocatable :: units_tij(:,:)
!@var DNAME_TIJ, DENOM_TIJ: Short names, indices of tij denominators.
!@+   Currently, denom is specified along with the standard metadata,
!@+   and the dnames are looked up afterward.  Exception: denominators
!@+   that are not part of the taijn array must be specified using dname.
!@+   For now the denom index indicates which _tracer_ index is to be
!@+   used as a denominator.  Ratios of 2 accumulation types are
!@+   not yet needed.
      character(len=sname_strlen),allocatable,dimension(:,:)::dname_tij
      integer, allocatable, dimension(:,:) :: denom_tij
!@var LNAME_TIJ: descriptions of tracer IJK diags
      character(len=lname_strlen),allocatable, dimension(:,:)::lname_tij

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
!@var ijs_XXX index for diags not specific to a certain tracer
      INTEGER :: ijs_isoprene,ijs_NO2_1030,ijs_NO2_1030c,
     &ijs_NO2_1330,ijs_NO2_1330c,ijts_Sdrydep,ijs_O3mass,
     &ijts_clrsky=0,ijts_pocean=0,ijts_sunlit_snow=0

!@param KTAIJS number of special lat/lon tracer diagnostics
!@+   please just increase this if needed - don't bother with pp options
      INTEGER,PARAMETER :: ktaijs=6500
!@param MaxSubCl Maximum number of sub classes of tracers for rad. diagnostics
      INTEGER,PARAMETER :: MaxSubCl=4
!@param MaxDMc Maximum number of special wet depo diags for MC clouds
!@param MaxDLs Maximum number of special wet depo diags for LS clouds
      INTEGER,PARAMETER :: MaxDMc=6,MaxDLs=6
!@param MaxSpec Maximum number special diagnostics not associated with specific
!@+ tracer
      INTEGER,PARAMETER :: MaxSpec=3
!@dbparam diag_rad switches on/off comprehensive radiative diags for tracers
      INTEGER :: diag_rad=0 ! =off (default)
!@dbparam diag_aod_3d outputs 3d aod and aaod (band6) for all-sky (=1),
!@+       clear-sky(=2), dry aerosol(=4), or all (=3). Notice weird order!
      INTEGER :: diag_aod_3d=0 ! =off (default)
!@dbparam save_dry_aod outputs 2d dry aod for all bands.
      INTEGER :: save_dry_aod=0 ! =off (default)
!@var TAIJS  lat/lon special tracer diagnostics; sources, sinks, etc.
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAIJS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAIJS_loc
!@var ijts_source tracer independent array for TAIJS surface src. diags
      INTEGER, ALLOCATABLE ::ijts_source(:,:)
!@var ijts_isrc tracer independent array for TAIJS interactive srf. src.
      INTEGER, ALLOCATABLE ::ijts_isrc(:,:)
!@var ijts_aq tracer independent array for TAIJS aqueous change
      INTEGER, ALLOCATABLE ::ijts_aq(:)
!@var ijts_alb BC impact on snow/ice albedo, grain size,sw and lw radiation
      INTEGER ijts_alb(2)
!@var ijts_tau tracer independent array for TAIJS hydrated opt. thick.
      INTEGER, ALLOCATABLE ::ijts_tau(:,:)
!@var ijts_tausub index for TAIJS opt. thick. for tracer sub classes
      INTEGER, ALLOCATABLE ::ijts_tausub(:,:,:)
!@var ijts_sqex index for TAIJS total extinction for 6 radiation bands
      INTEGER, ALLOCATABLE ::ijts_sqex(:,:,:)
!@var ijts_sqexsub index for TAIJS total extinction for 6 radiation bands for
!@+   tracer sub classes
      INTEGER, ALLOCATABLE ::ijts_sqexsub(:,:,:,:)
!@var ijts_sqsc index for TAIJS scattering extinction for 6 radiation bands
      INTEGER, ALLOCATABLE ::ijts_sqsc(:,:,:)
!@var ijts_sqscsub index for TAIJS scattering extinction for 6 radiation
!@+   bands for tracer sub classes
      INTEGER, ALLOCATABLE ::ijts_sqscsub(:,:,:,:)
!@var ijts_sqcb index for TAIJS sct asymmetry factor for 6 radiation bands
      INTEGER, ALLOCATABLE ::ijts_sqcb(:,:,:)
!@var ijts_sqcbsub index for TAIJS sct asymmetry factor for 6 radiation bands
!@+   for tracer sub classes
      INTEGER, ALLOCATABLE ::ijts_sqcbsub(:,:,:,:)
!@var ijts_fc tracer independent array for TAIJS SW/LW rad. forcings
      INTEGER, ALLOCATABLE ::ijts_fc(:,:)
#ifdef AUXILIARY_OX_RADF
!@var ijts_auxfc auxiliary Ox array for TAIJS SW/LW rad. forcings
      INTEGER ijts_auxfc(4)
#endif
!@var ijts_fcsub index for TAIJS SW/LW rad. forc. for tracer sub classes
      INTEGER, ALLOCATABLE ::ijts_fcsub(:,:,:)
!@var ijts_spec index for TAIJS diags. not associated with single tracer
      INTEGER :: ijts_spec(MaxSpec)
!@var ijts_gasex index for TAIJS associated with atm/oc gas exchange
      INTEGER, ALLOCATABLE :: ijts_gasex(:,:)
#ifdef TRACERS_AMP
!@var ijts_AMPpdf special diagnostic for not-transported tracers
      INTEGER ijts_AMPpdf(1,nbins)
#endif
!@var ijts_AMPe tracer independent array for extra diadnostics
      INTEGER, ALLOCATABLE ::ijts_AMPe(:)
!@var ijts_AMPp tracer independent array for AMP processes
      INTEGER, ALLOCATABLE ::ijts_AMPp(:,:)
#ifdef TRACERS_TOMAS
!@var ijts_TOMAS tracer independent array for TOMAS processes
      INTEGER, ALLOCATABLE ::ijts_TOMAS(:,:)
!@var ijts_TOMAS tracer independent array for TOMAS subgrid coagulation
      INTEGER, ALLOCATABLE ::ijts_subcoag(:)
#endif
!@var ijts_trdpmc indices of taijs special wet depo diags for MC clouds
      INTEGER, ALLOCATABLE :: ijts_trdpmc(:,:)
!@var ijts_trdpls indices of taijs special wet depo diags for LS clouds
      INTEGER, ALLOCATABLE :: ijts_trdpls(:,:)
!@var ijts_wet tracer independent array for TAIJS wet depo diagnostics
      INTEGER, ALLOCATABLE :: ijts_wet(:)
!@var ijts_3Dsource tracer independent array for TAIJS 3D src. diags
      INTEGER, ALLOCATABLE ::ijts_3Dsource(:,:)
!@var SNAME_IJTS, UNITS_IJTS: Names & units of lat-sigma tracer diags
      character(len=sname_strlen), dimension(ktaijs) :: sname_ijts
      character(len=units_strlen), dimension(ktaijs) :: units_ijts
!@var LNAME_IJTS: descriptions of tracer IJTS diags
      character(len=lname_strlen), dimension(ktaijs) ::
     &     lname_ijts = 'unused'
!@var DNAME_IJTS, DENOM_IJTS: Short names, indices of ijts denominators.
!@+   Currently, dname is specified along with the standard metadata and
!@+   the denom indices are looked up afterward.
      character(len=sname_strlen), dimension(ktaijs) :: dname_ijts=''
      integer, dimension(ktaijs) :: denom_ijts=0
!@var SCALE_IJTS: printout scaling factor for tracer IJTS diagnostics
      REAL*8, dimension(ktaijs) :: scale_ijts
!@var IA_IJTS: idacc-number for tracer source/sink IJ diags
      integer ia_ijts(ktaijs)
!@var ijts_power: power of 10 used for tracer IJ source/sink diags
      INTEGER, DIMENSION(ktaijs) :: ijts_power
!@var ijts_HasArea: does accumulation need to be divided by grid area?
      LOGICAL, DIMENSION(ktaijs) :: ijts_HasArea

C**** TAIJLS 3D special tracer diagnostics

!@param ktaijl number of TAIJLS tracer diagnostics;
      INTEGER, PARAMETER :: ktaijl=116
#ifdef ACCMIP_LIKE_DIAGS 
     &                            + 17
#endif
#ifdef SOA_DIAGS
     &                            + 12
#ifdef TRACERS_TERP
     &                            + 1
#endif  /* TRACERS_TERP */
     &                            + 16*nsoa
#endif  /* SOA_DIAGS */
!@var TAIJLS  3D tracer diagnostics (tracer dependent)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TAIJLS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TAIJLS_loc
!@var SNAME_IJLT: Names of 3D tracer IJL diagnostics
      character(len=sname_strlen), dimension(ktaijl) :: sname_ijlt
!@var DNAME_IJLT, DENOM_IJLT: Short names, indices of taijls denominators.
!@+   Currently, dname is specified along with the standard metadata and
!@+   the denom indices are looked up afterward.
      character(len=sname_strlen), dimension(ktaijl) :: dname_ijlt=''
      integer, dimension(ktaijl) :: denom_ijlt=0
!@var LNAME_IJLT,UNITS_IJLT: descriptions/units of 3D tracer diagnostics
      character(len=lname_strlen), dimension(ktaijl) ::
     &     lname_ijlt = 'unused'
      character(len=units_strlen), dimension(ktaijl) :: units_ijlt
!@var SCALE_IJLT: printout scaling factor for 3D tracer diagnostics
      REAL*8, dimension(ktaijl) :: scale_ijlt
!@var IR_IJLT: range index of IJL diagnostics
      integer, dimension(ktaijl) :: ir_ijlt
!@var IA_IJLT: accumulation index for IJL diagnostics
      integer, dimension(ktaijl) :: ia_ijlt
!@var ijlt_power: power of 10 used for tracer IJL 3D diags
      INTEGER, DIMENSION(ktaijs) :: ijlt_power
!@var ijlt_XXX diag names associated with 3D tracer special diags
      INTEGER :: ijlt_OHvmr,ijlt_OHconc,
     & ijlt_NO3,ijlt_HO2,ijlt_COp,ijlt_COd,
     & ijlt_Oxp,ijlt_Oxd,ijlt_CH4d,ijlt_OxpHO2,ijlt_OxpCH3O2,ijlt_OxpRO2
     & ,ijlt_OxlOH,ijlt_OxlHO2,ijlt_OxlALK,ijlt_phO1D,ijlt_pO1D,ijlt_pOH
     & ,ijlt_NOxLgt,ijlt_NOvmr,ijlt_NO2vmr,ijlt_JO1D,ijlt_JNO2
     & ,ijlt_JH2O2,ijlt_prodSO4aq,ijlt_prodSO4gs,ijlt_O3ppbv
     & ,ijlt_O3cmatm
     & ,ijlt_clrsky2d=0
     & ,ijlt_airmass=0
!@var ijlt_aH2O aerosol H2O from thermodynamics (ug/m3)
!@var ijlt_apH aerosol pH from thermodynamics (dimensionless)
      integer :: ijlt_aH2O,ijlt_apH
#ifdef SOA_DIAGS
!@var ijlt_soa_changeL_isoprene gas-phase changeL of isoprene SOA (ug/m3)
#ifdef TRACERS_TERP
!@var ijlt_soa_changeL_terpenes gas-phase changeL of terpenes (ug/m3)
#endif  /* TRACERS_TERP */
!@var ijlt_soa_voc2nox VOC/NOx ratio (ppbC/ppb)
!@var ijlt_soa_pcp Total non-volatile SOA-absorbing mass (ug/m3)
!@var ijlt_soa_aerotot PCP plus SOA (g+a) (ug/m3)
!@var ijlt_soa_aerotot_gas Gas-phase semivolatile potential SOA (g+a) (ug/m3)
!@var ijlt_soa_xmf_isop Molar fraction of isoprene SOA
!@var ijlt_soa_xmf_apin Molar fraction of a-pinene SOA
!@var ijlt_soa_zcoef_isop Activity coefficient for isoprene SOA
!@var ijlt_soa_zcoef_apin Activity coefficient for a-pinene SOA
!@var ijlt_soa_meanmw Mean organic aerosol molecular weight (g/mol)
!@var ijlt_soa_iternum Total iterations for SOA calculations (count)
!@var ijlt_soa_M0 Final M0 value
!@var ijlt_soa_y0_ug_g gas-phase y0_ug (ug/m3)
!@var ijlt_soa_y0_ug_a aerosol-phase y0_ug (ug/m3)
!@var ijlt_soa_y_ug_g gas-phase y_ug (ug/m3)
!@var ijlt_soa_y_ug_a aerosol-phase y_ug (ug/m3)
!@var ijlt_soa_changeL_g_before gas-phase changeL before SOA (ug/m3)
!@var ijlt_soa_changeL_a_before aerosol-phase changeL before SOA (ug/m3)
!@var ijlt_soa_changeL_g_after gas-phase changeL after SOA (ug/m3)
!@var ijlt_soa_changeL_a_after aerosol-phase changeL after SOA (ug/m3)
!@var ijlt_soa_apartmass Effective apartmass
!@var ijlt_soa_kpart Partitioning coefficient (m3/ug)
!@var ijlt_soa_kp Final partitioning coefficient (m3/ug)
!@var ijlt_soa_soamass Final soamass value
!@var ijlt_soa_partfact Final partfact value
!@var ijlt_soa_evap Previous time-step soa that survived evaporation (ug/m3)
!@var ijlt_soa_cond Previous time-step soa precursor that condensed (ug/m3)
!@var ijlt_soa_chem Condensed soa precursor that was produced chemically at this time-step (ug/m3)
      integer :: ijlt_soa_changeL_isoprene
#ifdef TRACERS_TERP
      integer :: ijlt_soa_changeL_terpenes
#endif  /* TRACERS_TERP */
      integer :: ijlt_soa_voc2nox, ijlt_soa_pcp, ijlt_soa_aerotot,
     $     ijlt_soa_aerotot_gas, ijlt_soa_xmf_isop, ijlt_soa_xmf_apin,
     $     ijlt_soa_zcoef_isop, ijlt_soa_zcoef_apin, ijlt_soa_meanmw,
     $     ijlt_soa_iternum,  ijlt_soa_m0
      integer, dimension(nsoa) :: ijlt_soa_y0_ug_g,  ijlt_soa_y0_ug_a,
     $     ijlt_soa_y_ug_g, ijlt_soa_y_ug_a, ijlt_soa_changeL_g_before,
     $     ijlt_soa_changeL_a_before, ijlt_soa_changeL_g_after,
     $     ijlt_soa_changeL_a_after, ijlt_soa_apartmass, ijlt_soa_kpart,
     $     ijlt_soa_kp, ijlt_soa_soamass, ijlt_soa_partfact,
     $     ijlt_soa_evap, ijlt_soa_cond, ijlt_soa_chem
#endif  /* SOA_DIAGS */
#ifdef TRACERS_AMP
!@var ijlt_AMPm tracer independent array for AMP modes
      integer, allocatable :: ijlt_AMPm(:,:)
#endif 
!@var ijlt_3Dtau 3D tracer independent array for all-sky hydrated opt. thick.
      integer, allocatable :: ijlt_3Dtau(:)
!@var ijlt_3DtauCS 3D tracer independent array for clear-sky hydrated opt. thick.
      integer, allocatable :: ijlt_3DtauCS(:)
!@var ijlt_3DtauDRY 3D tracer independent array for dry opt. thick.
      integer, allocatable :: ijlt_3DtauDRY(:)
!@var ijlt_3Daaod 3D tracer independent array for all-sky hydrated absorption
      INTEGER, allocatable :: ijlt_3Daaod(:)
!@var ijlt_3DaaodCS 3D tracer independent array for clear-sky hydrated absorption
      INTEGER, allocatable :: ijlt_3DaaodCS(:)
!@var ijlt_3DaaodDRY 3D tracer independent array for dry absorption
      INTEGER, allocatable :: ijlt_3DaaodDRY(:)
#ifdef SAVE_AEROSOL_3DMASS_FOR_NINT
!@var ijlt_3Dmass 3D tracer independent array for layer MASS (or load)
      INTEGER, allocatable :: ijlt_3Dmass(:)
#endif
#ifdef TRACERS_TOMAS
!@var ijlt_ccn_01-ccn_03 CCN diagnostic
      INTEGER :: ijlt_ccn_01,ijlt_ccn_02,ijlt_ccn_03
#endif
C**** TAJLN
!@param ktajl,ktajlx number of TAJL tracer diagnostics;
!@+          ktajlx includes composites
      INTEGER, PARAMETER :: ktajl=10, ktajlx=ktajl+2
!@var TAJLN  vertical tracer diagnostics (all tracers)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: TAJLN
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TAJLN_loc
!@var jlnt_xx Names for TAJLN diagnostics
      INTEGER jlnt_conc,jlnt_mass,jlnt_nt_tot,jlnt_nt_mm,jlnt_vt_tot,
     &  jlnt_vt_mm,jlnt_mc,jlnt_turb,jlnt_lscond, jlnt_bebe,
     &  jlnt_bepb,jlnt_cldh2o

!@var SNAME_JLN: Names of lat-sigma tracer JL diagnostics
      character(len=sname_strlen), allocatable ::sname_jln(:,:)
!@var LNAME_JLN,UNITS_JLN: descriptions/units of tracer JL diagnostics
      character(len=lname_strlen), allocatable :: lname_jln(:,:)
      character(len=units_strlen), allocatable :: units_jln(:,:)
!@var SCALE_JLQ: printout scaling factor for tracer JL diagnostics
      REAL*8, dimension(ktajlx) :: scale_jlq
!@var IA_JLQ,JGRID_JLQ: idacc-numbers,gridtypes for tracer JL diags
      integer, dimension(ktajlx) :: ia_jlq,jgrid_jlq
!@var JLQ_POWER: power of 10 used for tracer JL diagnostics
!@+        It is associated with a specific physical process
      integer jlq_power(ktajlx)
!@var scale_jln: Scale for jl maps
      REAL*8, allocatable, DIMENSION(:) :: scale_jln

C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>
!@param ktajls number of source/sink TAJLS tracer diagnostics;
!@+   please just increase this if needed - do not bother with pp options
#ifndef TRACERS_TOMAS
      INTEGER,PARAMETER :: ktajls=1350
#else
      INTEGER,PARAMETER :: ktajls=3400
#endif
!@var jls_XXX index for non-tracer specific or special diags
      INTEGER jls_OHconk,jls_HO2con,jls_NO3,jls_O3vmr
     *     ,jls_phot,jls_OHcon,jls_H2Omr
     *     ,jls_N2O5sulf,jls_day,jls_COd,jls_COp,jls_Oxd,jls_Oxp
     *     ,jls_ClOcon,jls_H2Ocon,jls_H2Ochem,jls_OxdT,jls_OxpT
      integer, allocatable :: jls_incloud(:,:)
!@var TAJLS  JL special tracer diagnostics for sources, sinks, etc
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAJLS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAJLS_loc
!@var jls_source tracer independent array for TAJLS surface src. diags
      INTEGER, ALLOCATABLE :: jls_source(:,:)
!@var jls_isrc tracer indep. array for TAJLS interactive surface src. diags
      INTEGER, ALLOCATABLE :: jls_isrc(:,:)
!@var jls_3Dsource tracer independent array for TAJLS 3D source diags
      INTEGER, ALLOCATABLE :: jls_3Dsource(:,:)
!@var jls_decay tracer independent array for radioactive sinks
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jls_decay
!@var jls_grav tracer independent array for grav. settling sink
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jls_grav
!@var jls_prec tracer independent array for precipitation/wet dep
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jls_prec
!@var jls_trdpmc14 indices of tajls special wet depo diags for MC clouds
      INTEGER, ALLOCATABLE:: jls_trdpmc(:,:)
!@var jls_trdpls indices of tajls special wet depo diags for LS clouds
      INTEGER, ALLOCATABLE:: jls_trdpls(:,:)
!@var jls_wet tracer independent array for wet deposition (for old dust)
      INTEGER, allocatable, DIMENSION(:) :: jls_wet

!@var jls_spec index for TAJLS for special diagnostics not associated with
!@var jls_spec single tracer
      INTEGER :: jls_spec(MaxSpec)
!@var jwt_jls: Weighting index for jls diags 1=simple average, 2=by area
      integer, dimension(ktajls) :: jwt_jls
!@var SNAME_JLS: Names of lat-sigma tracer JL sources/sinks
      character(len=sname_strlen), dimension(ktajls) :: sname_jls
!@var LNAME_JLS,UNITS_JLS: descriptions/units of tracer JLS diags
      character(len=lname_strlen), dimension(ktajls) ::
     &     lname_jls = 'unused'
      character(len=units_strlen), dimension(ktajls) :: units_jls
!@var SCALE_JLS: printout scaling factor for tracer JLS diagnostics
      REAL*8, dimension(ktajls) :: scale_jls
!@var IA_JLS,JGRID_JLS: idacc-numbers,gridtypes for tracer JL diags
      integer, dimension(ktajls) :: ia_jls,jgrid_jls
!@var JLS_POWER: power of 10 used for tracer JLS diagnostics
      integer jls_power(ktajls)
!@var jls_ltop: Top layer for this diagnostic
      integer jls_ltop(ktajls)
#endif  /* TRACERS_ON */

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** TCONSRV
!@param NTCONS Maximum Number of special tracer conservation points
      INTEGER, PARAMETER :: ntcons=20
#ifdef TRACERS_SPECIAL_Lerner
     &                             +10
#endif
#ifdef TRACERS_AMP
     &                             +4
#endif
#ifdef TRACERS_TOMAS
     &                             +6+15
#endif
!@param npts_common total number of conservation diagnostics outside
!@+                 those defined for tracers
      integer, parameter :: npts_common=npts+1
!@param KTCON total number of conservation diagnostics for tracers
      INTEGER, PARAMETER :: KTCON=npts_common+ntcons+1
!@param ntmxcon total number of conservation quantities
      integer :: ntmxcon
      logical :: qcon(KTCON-1), qsum(KTCON-1)

!@var TCONSRV conservation diagnostics for tracers
      REAL*8, allocatable, DIMENSION(:,:,:) :: TCONSRV,TCONSRV_loc 
!@var SCALE_TCON scales for tracer conservation diagnostics
      REAL*8, allocatable, DIMENSION(:,:) :: SCALE_TCON
!@var TITLE_TCON titles for tracer conservation diagnostics
      CHARACTER*38, allocatable, DIMENSION(:,:) :: TITLE_TCON
!@var IA_TCON IDACC numbers for tracer conservation diagnostics
      INTEGER, allocatable, DIMENSION(:,:) ::  IA_TCON
!@var NSUM_TCON indices for summation of conservation diagnostics
      INTEGER, allocatable, DIMENSION(:,:) :: NSUM_TCON
!@var NOFMT indices for TCONSRV array
      INTEGER, allocatable, DIMENSION(:,:) :: NOFMT
!@var CONPTS names of special processes for tracer conservation diags
      CHARACTER*16, DIMENSION(ntcons) :: CONPTS=''
!@var kt_power_inst,kt_power_change: Exponents for tracer conservation
      INTEGER, allocatable, DIMENSION(:):: kt_power_inst,kt_power_change
!@var name_tconsrv,lname_tconsrv,units_tconsrv: for tracer conservation
      character(len=sname_strlen), allocatable, dimension(:,:) :: 
     &     name_tconsrv
      character(len=units_strlen), allocatable, dimension(:,:) ::
     &     units_tconsrv
      character(len=lname_strlen), allocatable, dimension(:,:) ::
     &     lname_tconsrv
!@var SCALE_INST,SCALE_CHANGE: Scale factors for tracer conservation
      REAL*8, allocatable, dimension(:) :: SCALE_INST,SCALE_CHANGE
!@var itcon_surf Index array for surface source/sink conservation diags
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itcon_surf
!@var itcon_3Dsrc Index array for 3D source/sink conservation diags
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itcon_3Dsrc
!@var itcon_decay Index array for decay conservation diags
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itcon_decay
!@var itcon_mc Index array for moist convection conserv. diags
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itcon_mc
!@var itcon_ss Index array for large-scale condensation conserv. diags
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itcon_ss
!@var itcon_amp Index array for microphysical processes diags
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itcon_amp
!@var itcon_amp Index array for microphysical processes diags
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itcon_ampe
!@var itcon_amp Index array for microphysical processes diags
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itcon_ampm
!@var itcon_dd Index array for dry deposition conserv. diags
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itcon_dd
!@var itcon_wt Index array for dust/mineral dust deposition conserv. diags
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itcon_wt
!@var natmtrcons, nocntrcons number of atmospheric/ocean tcon diags
      INTEGER :: natmtrcons=0, nocntrcons=0
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: maxntmocn=3
#else 
#ifdef TRACERS_OCEAN
      integer, parameter :: maxntmocn=25
#else
      integer, parameter :: maxntmocn=0
#endif
#endif
#ifdef TRACERS_TOMAS
!@var itcon_TOMAS Index array for microphysical processes diags
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itcon_TOMAS
!@var itcon_TOMAS Index array for subgrid coagulation diags
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itcon_subcoag
#endif
#endif  /* TRACERS_ON  or  TRACERS_OCEAN */

!@var PDSIGJL temporary storage for mean pressures for jl diags
      REAL*8, DIMENSION(JM,LM)            :: PDSIGJL

#ifdef TRACERS_WATER
!@var TRP_acc, TRE_acc accumulation arrays for some SUBDD diags
!!    REAL*8 TRP_acc(ntm,IM,JM), TRE_acc(ntm,IM,JM)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: TRP_acc,TRE_acc
#endif

!@var trcsurf global array of tracer mixing ratio at surface [kg/kg]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: trcsurf
!@var trcSurfByVol global array of tracer concentration at surface [kg/m^3]
      real(kind=8),allocatable,dimension(:,:,:) :: trcSurfByVol

!@var trcSurfMixR_acc global array of tracers to accumulate mxixing ratio at
!@var                 surface for subdd diagnostics [kg/kg]
!@var trcSurfByVol_acc global array of tracers to accumulate concentration at
!@var                  surface for subdd diagnostics [kg/m^3]
      real(kind=8),allocatable,dimension(:,:,:) :: trcSurfMixR_acc
     &     ,trcSurfByVol_acc

#ifdef TRACERS_ON
! This section declares arrays used to write tracer
! diagnostics accumulations in a format suitable for offline
! postprocessing.  The size of these arrays cannot be known
! precisely a priori due to the large number of outputs not
! declared in advance.  Thus, metadata arrays are declared with
! a larger-than-needed size, and acc arrays are re-allocated
! to the correct size after the size is determined.
! The "n" and "s" instances of diagnostics classes are
! merged as follows:
!@var taijl_out combines taijln and taijls.  Distributed.
!@var taij_out  combines taijn  and taijs.   Distributed.
!@var tajl_out  combines tajln  and tajls.
!@var tconsrv_out combines the ktcon/ntmxcon dims of tconsrv
!@+               and omits its unused elements.
! All contents of txxx_out arrays are PER UNIT AREA.
! Denominator information is stored in denom_xxx.
! Metadata for scaled outputs is consolidated in cdl_xxx.
! 
      ! This was NTM, but NTM is now dynamic. Introducing MAXNTM instead.
      integer, parameter :: MAXNTM = 1000 ! rather than making all ktaij_ arrays allocatable - TLC
      integer, parameter :: ktaij_ = (ktaij*MAXNTM+ktaijs
#ifdef TRACERS_SPECIAL_O18
     &     + 4                  ! include dexcess + D17O diags
#endif
#ifdef TRACERS_DRYDEP
     &     + MAXNTM             ! include dry dep % diags
#endif
     &     )*3/2  ! make 50% larger for denoms and extra specials
      integer :: ktaij_out ! actual number of qtys in taij_out
      real*8, dimension(:,:,:), allocatable :: taij_out
      integer, dimension(ktaij_) :: ir_taij,ia_taij,denom_taij
      character(len=lname_strlen), dimension(ktaij_) :: lname_taij
      character(len=sname_strlen), dimension(ktaij_) :: sname_taij
      character(len=units_strlen), dimension(ktaij_) :: units_taij
      real*8, dimension(ktaij_) :: scale_taij
      type(cdl_type) :: cdl_taij,cdl_taij_latlon
      real*8, dimension(:,:,:), allocatable :: hemis_taij

      integer :: ktaijl_
      integer :: ktaijl_out ! actual number of qtys in taijl_out
      real*8, dimension(:,:,:,:), allocatable :: taijl_out
      integer, allocatable, dimension(:) ::ir_taijl,ia_taijl,denom_taijl
      character(len=lname_strlen),allocatable,dimension(:) ::lname_taijl
      character(len=sname_strlen),allocatable,dimension(:) ::sname_taijl
      character(len=units_strlen),allocatable,dimension(:) ::units_taijl
      real*8, allocatable, dimension(:) :: scale_taijl
      type(cdl_type) :: cdl_taijl,cdl_taijl_latlon

      integer :: ktajl_
      integer :: ktajl_out ! actual number of qtys in tajl_out
      real*8, dimension(:,:,:), allocatable :: tajl_out
      integer, allocatable, dimension(:) :: pow_tajl,ia_tajl,denom_tajl,
     &     jgrid_tajl,lgrid_tajl,ltop_tajl
      character(len=lname_strlen),allocatable,dimension(:) :: lname_tajl
      character(len=sname_strlen),allocatable,dimension(:) :: sname_tajl
      character(len=units_strlen),allocatable,dimension(:) :: units_tajl
      real*8, allocatable, dimension(:) :: scale_tajl
      type(cdl_type) :: cdl_tajl
      real*8, dimension(:,:,:), allocatable :: hemis_tajl,vmean_tajl
#endif /* TRACERS_ON */

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      integer :: ktcon_out ! actual number of qtys in tconsrv_out
      real*8, dimension(:,:), allocatable :: tconsrv_out,hemis_tconsrv
      real*8, dimension(:), allocatable :: scale_tcon_out
      integer, dimension(:), allocatable ::  ia_tcon_out
      character(len=sname_strlen), dimension(:), allocatable ::
     &     sname_tconsrv_out
      type(cdl_type) :: cdl_tconsrv
#endif

      target :: taijn_loc,tconsrv_loc,nofmt

      END MODULE TRDIAG_COM

#ifdef TRACERS_ON
      SUBROUTINE SET_TCON(QCON,NAME_CON,QSUM,INST_UNIT,SUM_UNIT
     *     ,INST_SC,CHNG_SC, itr,CONPTs,CONPT0)
!@sum  SET_TCON assigns conservation diagnostic array indices
!@auth Gavin Schmidt
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE MODEL_COM, only: dtsrc
      USE DYNAMICS, only : nfiltr
      USE DIAG_COM, only: npts,ia_d5d,ia_d5s,ia_filt,ia_12hr,ia_src
      USE TRDIAG_COM, only: ktcon,title_tcon,scale_tcon,nsum_tcon
     *     ,nofmt,ia_tcon,name_tconsrv,lname_tconsrv,units_tconsrv
     *     ,ntcons,npts_common
      IMPLICIT NONE
!@var QCON denotes at which points conservation diags are saved
      LOGICAL, INTENT(IN),DIMENSION(ktcon-1) :: QCON
!@var QSUM sets whether each diag is included in final sum
!@+   should be zero for diags set using DIAGTCB (i.e. using difference)
      LOGICAL, INTENT(IN),DIMENSION(ktcon-1) :: QSUM
      LOGICAL, DIMENSION(ktcon-1) :: QSUM_CON   ! local version
!@var INST_SC scale for instantaneous value
      REAL*8, INTENT(IN) :: INST_SC
!@var CHNG_SC scale for changes
      REAL*8, INTENT(IN) :: CHNG_SC
!@var NAME_CON name of conservation quantity
      CHARACTER(len=*), INTENT(IN) :: NAME_CON
!@var sname name of conservation quantity (no spaces)
      CHARACTER*8 :: sname
!@var INST_UNIT string for unit for instant. values
      CHARACTER*20, INTENT(IN) :: INST_UNIT
!@var SUM_UNIT string for unit for summed changes
      CHARACTER*20, INTENT(IN) :: SUM_UNIT
!@var ITR index for the tracer
      INTEGER, INTENT(IN) :: ITR
!@var CONPTS conservation diag points for special tracer diags
      CHARACTER*16, INTENT(IN), DIMENSION(ntcons) :: CONPTS
!@var CONPTS_sname like CONPTS but without spaces
      CHARACTER*16, DIMENSION(ntcons) :: CONPTS_sname
!@var CONPT0 conservation diag points for common points
      CHARACTER*10, INTENT(IN), DIMENSION(npts) :: CONPT0
!@var CONPT0_sname like CONPT0 but without spaces
      CHARACTER*10, DIMENSION(npts) :: CONPT0_sname
      CHARACTER*11 CHGSTR
      INTEGER NI,NM,NS,N,k
      CHARACTER*40 clean_str

C**** make nice netcdf names
      sname=trim(clean_str(name_con))
      do n=1,ntcons
         conpts_sname(n) = trim(clean_str(conpts(n)))
      enddo
      do n=1,npts
         conpt0_sname(n) = trim(clean_str(conpt0(n)))
      enddo
C****
      NI=1
      NOFMT(1,itr) = NI
      TITLE_TCON(NI,itr) =
     *     "0INSTANT "//TRIM(NAME_CON)//" "//TRIM(INST_UNIT)
      SCALE_TCON(NI,itr) = INST_SC
      name_tconsrv(NI,itr) ="inst_"//sname
      lname_tconsrv(NI,itr) = "INSTANT "//TRIM(NAME_CON)
      units_tconsrv(NI,itr) = INST_UNIT
      NSUM_TCON(NI,itr) = -1
      IA_TCON(NI,itr) = 12
      NM=NI
      DO N=2,ktcon-1
        IF (QCON(N)) THEN
          NM = NM + 1
          NOFMT(N,itr) = NM
          QSUM_CON(NM)=.FALSE.
          IF (QSUM(N)) QSUM_CON(NM)=.TRUE.
          CHGSTR=" CHANGE OF "
          if (n.le.npts_common) then
            TITLE_TCON(NM,itr) = CHGSTR//TRIM(NAME_CON)//" BY "//
     *         CONPT0(N-1)
            name_tconsrv(NM,itr) =
     *           "chg_"//trim(sname)//"_"//TRIM(CONPT0_sname(N-1))
          else
            IF (.not. QSUM(N)) CHGSTR="     DELTA "
            TITLE_TCON(NM,itr) = CHGSTR//TRIM(NAME_CON)//" BY "//
     *           CONPTs(N-npts_common)
            name_tconsrv(NM,itr) =
     *           "chg_"//trim(sname)//"_"//
     *           TRIM(CONPTs_sname(N-npts_common))
          end if
          lname_tconsrv(NM,itr) = TITLE_TCON(NM,itr)
          units_tconsrv(NM,itr) = SUM_UNIT
          SELECT CASE (N)
          CASE (2)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_d5d
          CASE (3,4,5,6,7,9,11,12)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_d5s
          CASE (8)
            SCALE_TCON(NM,itr) = CHNG_SC/(NFILTR*DTSRC)
            IA_TCON(NM,itr) = ia_filt
          CASE (10)
            SCALE_TCON(NM,itr) = CHNG_SC*2./SECONDS_PER_DAY
            IA_TCON(NM,itr) = ia_12hr
          CASE (13:)   ! special tracer sources
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_src
          END SELECT
        ELSE
          NOFMT(N,itr) = 0
        END IF
      END DO
      NS=NM+1
      IF (NS.gt.KTCON) THEN
        WRITE(6,*) "KTCON not large enough for extra conserv diags",
     *       KTCON,NI,NM,NS,NAME_CON
        call stop_model(
     &       "Change KTCON in tracer diagnostic common block",255)
      END IF
      DO NM=NI+1,NS-1
        NSUM_TCON(NM,itr) = -1
        IF (QSUM_CON(NM)) NSUM_TCON(NM,itr) = NS
      END DO
      TITLE_TCON(NS,itr) = " SUM OF CHANGES "//TRIM(SUM_UNIT)
      name_Tconsrv(NS,itr) ="sum_chg_"//trim(sname)
      lname_Tconsrv(NS,itr) = " SUM OF CHANGES OF "//TRIM(NAME_CON)
      units_Tconsrv(NS,itr) = SUM_UNIT
      SCALE_TCON(NS,itr) = 1.
      IA_TCON(NS,itr) = 12
      NSUM_TCON(NS,itr) = 0
      RETURN
      END SUBROUTINE set_tcon

#ifdef TRACERS_OCEAN
      SUBROUTINE SET_TCONO(NAME_CON,INST_UNIT,SUM_UNIT,
     &    INST_SC,CHNG_SC, itr0, extra_pt_n, extra_pt_idx, extra_pt_str)
!@sum  SET_TCONO assigns ocean conservation diagnostic array indices
!@auth Gavin Schmidt
      USE TimeConstants_mod, only: SECONDS_PER_DAY
      USE MODEL_COM, only: dtsrc
      USE DIAG_COM, only: npts,ia_d5s,ia_12hr,ia_src,conpt0
      USE TRDIAG_COM, only: ktcon,title_tcon,scale_tcon,nsum_tcon
     *     ,nofmt,ia_tcon,name_tconsrv,lname_tconsrv,units_tconsrv
     *     ,natmtrcons,nocntrcons,maxntmocn
      IMPLICIT NONE
!@var NAME_CON name of conservation quantity
      CHARACTER*8, INTENT(IN) :: NAME_CON
!@var INST_UNIT string for unit for instant. values
      CHARACTER*20, INTENT(IN) :: INST_UNIT
!@var SUM_UNIT string for unit for summed changes
      CHARACTER*20, INTENT(IN) :: SUM_UNIT
!@var INST_SC scale for instantaneous value
      REAL*8, INTENT(IN) :: INST_SC
!@var CHNG_SC scale for changes
      REAL*8, INTENT(IN) :: CHNG_SC
!@var ITR index for the tracer
      INTEGER, INTENT(IN) :: ITR0
!@var extra_pt_n number of extra points
      integer, intent(in) :: extra_pt_n
!@var extra_pt_n extra point indices
      integer, dimension(extra_pt_n), intent(in) :: extra_pt_idx
!@var extra_pt_str extra point labels
      character*10, dimension(extra_pt_n), intent(in) :: extra_pt_str
!@var QCON denotes at which points conservation diags are saved
      LOGICAL, DIMENSION(ktcon) :: QCON
!@var QSUM sets whether each diag is included in final sum
!@+   should be zero for diags set using DIAGTCB (i.e. using difference)
      LOGICAL, DIMENSION(ktcon) :: QSUM
      LOGICAL, DIMENSION(ktcon) :: QSUM_CON   ! local version
!@var sname name of conservation quantity (no spaces)
      CHARACTER*8 :: sname
!@var CONPT0_sname like CONPT0 but without spaces
      CHARACTER*10, DIMENSION(ktcon) :: CONPT0_sname, CONPT
      CHARACTER*11 CHGSTR
      CHARACTER*40 clean_str
      INTEGER NI,NM,NS,N,k,itr
      LOGICAL, PARAMETER :: T=.true., F=.false.

      if (itr0>maxntmocn) call
     &     stop_model('trdiag_com: increase maxntmocn', 255)
      nocntrcons=max(nocntrcons,itr0)
      CONPT(1:npts)=CONPT0
      CONPT(8)="OCN PHYS"
      QCON=F
      QCON(1:npts)=(/ F, F, F, T, F, F, F, T, T, T, T/)
      QSUM(:) = T

C**** make nice netcdf names
      sname=trim(clean_str(name_con))
      do n=1,npts
         conpt0_sname(n) = trim(clean_str(conpt(n)))
      enddo
      do k=1,extra_pt_n
         n=extra_pt_idx(k)
         conpt(n)=extra_pt_str(k)
         conpt0_sname(n)=trim(clean_str(conpt(n)))
         qcon(n)=T
         qsum(n)=T
      end do
C****
      NI=1
      itr=itr0+natmtrcons
      NOFMT(1,itr) = NI
      TITLE_TCON(NI,itr) =
     *     "0INSTANT "//TRIM(NAME_CON)//" "//TRIM(INST_UNIT)
      SCALE_TCON(NI,itr) = INST_SC
      name_tconsrv(NI,itr) ="inst_oc_"//sname
      lname_tconsrv(NI,itr) = "INSTANT "//TRIM(NAME_CON)
      units_tconsrv(NI,itr) = INST_UNIT
      NSUM_TCON(NI,itr) = -1
      IA_TCON(NI,itr) = 12
      NM=NI
      DO N=2,ktcon
        IF (QCON(N-1)) THEN
          NM = NM + 1
          NOFMT(N,itr) = NM
          QSUM_CON(NM)=.FALSE.
          IF (QSUM(N-1)) QSUM_CON(NM)=.TRUE.
          CHGSTR=" CHANGE OF "
            TITLE_TCON(NM,itr) = CHGSTR//TRIM(NAME_CON)//" BY "//
     *         CONPT(N-1)
            name_tconsrv(NM,itr) =
     *           "chg_oc_"//trim(sname)//"_"//TRIM(CONPT0_sname(N-1))
          lname_tconsrv(NM,itr) = TITLE_TCON(NM,itr)
          units_tconsrv(NM,itr) = SUM_UNIT
          SELECT CASE (N)
          CASE (2,8) ! ocean not affected by certain atm processes
            call stop_model('nonsensical choice in set_tcono',255)
          CASE (3,4,5,6,7,9,11,12)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_d5s
          CASE (10)
            SCALE_TCON(NM,itr) = CHNG_SC*2./SECONDS_PER_DAY
            IA_TCON(NM,itr) = ia_12hr
          CASE (13:)   ! special tracer sources
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_src
          END SELECT
        ELSE
          NOFMT(N,itr) = 0
        END IF
      END DO
      NS=NM+1
      IF (NS.gt.KTCON) THEN
        WRITE(6,*) "KTCON not large enough for extra conserv diags",
     *       KTCON,NI,NM,NS,NAME_CON
        call stop_model(
     &       "Change KTCON in tracer diagnostic common block",255)
      END IF
      DO NM=NI+1,NS-1
        NSUM_TCON(NM,itr) = -1
        IF (QSUM_CON(NM)) NSUM_TCON(NM,itr) = NS
      END DO
      TITLE_TCON(NS,itr) = " SUM OF CHANGES "//TRIM(SUM_UNIT)
      name_Tconsrv(NS,itr) ="sum_chg_oc_"//trim(sname)
      lname_Tconsrv(NS,itr) = " SUM OF CHANGES OF "//TRIM(NAME_CON)
      units_Tconsrv(NS,itr) = SUM_UNIT
      SCALE_TCON(NS,itr) = 1.
      IA_TCON(NS,itr) = 12
      NSUM_TCON(NS,itr) = 0
      RETURN
      END SUBROUTINE set_tcono
#endif /* TRACERS_OCEAN */

#endif

      subroutine write_src_dist_data(fid, def)
      use pario, only: defvar, write_dist_data, write_data
      use domain_decomp_atm, only : grid
      use trdiag_com, only: taijln=>taijln_loc,
     &                                   taijn=>taijn_loc, tij_prec
      use oldtracer_mod, only: src_dist_base, src_dist_index
      use tracer_com, only: ntm, xyztr, ntm_sph, ntm_reg
      implicit none
      integer, intent(in) :: fid
      logical, intent(in) :: def
      real*8, dimension(:, :, :, :), allocatable, save :: tr2
      real*8, dimension(:, :, :, :, :), allocatable, save :: tr3
      integer, dimension(ntm) :: lst
      integer :: i, n, ndist, nindex

      if (def) then
        ndist=0
        nindex=0
        do n=1, ntm
          if (src_dist_index(n)==1) then
            ndist=ndist+1
            lst(n)=ndist
          end if
          nindex=max(nindex, src_dist_index(n))
        end do
        if (ndist==0) return
        allocate(tr2(size(taijn, 1), size(taijn, 2), ndist, nindex))
        allocate(tr3(size(taijln, 1), size(taijln, 2), size(taijln, 3),
     &    ndist, nindex))
        do n=1, ntm
          if (src_dist_index(n)==1) then
            do i=n, ntm
              if (src_dist_base(i)==src_dist_base(n)) then
                tr2(:, :, lst(n), src_dist_index(i))=
     &                                 taijn(:, :, tij_prec, i)
                tr3(:, :, :, lst(n), src_dist_index(i))=
     &                                 taijln(:, :, :, i)
              end if
            end do
          end if
        end do
        call defvar(grid,fid,tr2,
     &            'src_dist2(dist_im,dist_jm,ndist,nbasis)')
        call defvar(grid,fid,tr3,
     &            'src_dist3(dist_im,dist_jm,lm,ndist,nbasis)')
        call defvar(grid,fid,xyztr,
     &            'src_dist_basis(nbasis,dist_im,dist_jm)')
        call defvar(grid,fid,ntm_sph,'ntm_sph')
        call defvar(grid,fid,ntm_reg,'ntm_reg')
      else
        if (allocated(tr2).and.allocated(tr3)) then
          call write_dist_data(grid,fid,'src_dist2',tr2)
          call write_dist_data(grid,fid,'src_dist3',tr3)
          call write_data(grid,fid,'ntm_sph',ntm_sph)
          call write_data(grid,fid,'ntm_reg',ntm_reg)
          deallocate(tr2)
          deallocate(tr3)
        end if
        if (allocated(xyztr))
     &    call write_dist_data(grid,fid,'src_dist_basis',xyztr,jdim=3)
      endif
      end subroutine write_src_dist_data

#ifdef TRACERS_ON
      subroutine def_rsf_trdiag(fid,r4_on_disk)
!@sum  def_rsf_trdiag defines tracer diag array structure in restart+acc files
!@auth M. Kelley
!@ver  beta
      use trdiag_com, only :
     &     TAIJLN=>TAIJLN_loc,
     &     TAIJLS=>TAIJLS_loc,
     &     TAIJN=>TAIJN_loc,
     &     TAIJS=>TAIJS_loc,
     &     TAJLN,
     &     TAJLS,
     &     TAIJL=>TAIJL_out,
     &     TAIJ=>TAIJ_out,
     &     TAJL=>TAJL_out
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid           !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
      if(r4_on_disk) then ! acc file
        call defvar(grid,fid,taijl,
     &       'taijl(dist_im,dist_jm,lm,ktaijl)',r4_on_disk=.true.)
        call defvar(grid,fid,taij,
     &       'taij(dist_im,dist_jm,ktaij)',r4_on_disk=.true.)
        call defvar(grid,fid,tajl,
     &       'tajl(jm_budg,lm,ktajl)',r4_on_disk=.true.)
        call write_src_dist_data(fid, .true.)
      else
        call defvar(grid,fid,taijln,'taijln(dist_im,dist_jm,lm,ntm)')
        call defvar(grid,fid,taijls,'taijls(dist_im,dist_jm,lm,ktaijl)')
        call defvar(grid,fid,taijs,'taijs(dist_im,dist_jm,ktaijs)')
        call defvar(grid,fid,taijn,'taijn(dist_im,dist_jm,ktaij,ntm)')
        call defvar(grid,fid,tajln,'tajln(jm_budg,lm,ktajlx,ntm)')
        call defvar(grid,fid,tajls,'tajls(jm_budg,lm,ktajls)')
      endif

      call def_rsf_tcons(fid,r4_on_disk)

      return
      end subroutine def_rsf_trdiag

      subroutine new_io_trdiag(fid,iaction)
!@sum  new_io_trdiag read/write tracer acc arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite,iowrite_single
      use trdiag_com, only :
     &     TAIJLN=>TAIJLN_loc,
     &     TAIJLS=>TAIJLS_loc,
     &     TAIJN=>TAIJN_loc,
     &     TAIJS=>TAIJS_loc,
     &     TAJLN,
     &     TAJLS,
     &     TAIJL=>TAIJL_out,
     &     TAIJ=>TAIJ_out,
     &     TAJL=>TAJL_out
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite_single)     ! output to acc file
        call write_dist_data(grid,fid,'taijl',taijl)
        call write_dist_data(grid,fid,'taij',taij)
        call write_data(grid,fid,'tajl',tajl)
        call write_src_dist_data(fid, .false.)
      case (iowrite)            ! output to restart file
        call gather_zonal_trdiag
        call write_dist_data(grid,fid,'taijln',taijln)
        call write_dist_data(grid,fid,'taijls',taijls)
        call write_dist_data(grid,fid,'taijs',taijs)
        call write_dist_data(grid,fid,'taijn',taijn)
        call write_data(grid,fid,'tajln',tajln)
        call write_data(grid,fid,'tajls',tajls)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'taijln',taijln)
        call read_dist_data(grid,fid,'taijls',taijls)
        call read_dist_data(grid,fid,'taijs',taijs)
        call read_dist_data(grid,fid,'taijn',taijn)
        call read_data(grid,fid,'tajln',tajln)
        call read_data(grid,fid,'tajls',tajls)
        call scatter_zonal_trdiag
      end select

      call new_io_tcons(fid,iaction)

      return
      end subroutine new_io_trdiag

      subroutine def_meta_trdiag(fid)
!@sum  def_meta_trdiag defines tracer metadata in acc files
!@auth M. Kelley
!@ver  beta
      use trdiag_com
      use pario, only : defvar,write_attr
      use domain_decomp_atm, only : grid
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_attr(grid,fid,'taij','reduction','sum')
      call write_attr(grid,fid,'taij','split_dim',3)
      call defvar(grid,fid,ia_taij(1:ktaij_out),
     &     'ia_taij(ktaij)')
      call defvar(grid,fid,denom_taij(1:ktaij_out),
     &     'denom_taij(ktaij)')
      call defvar(grid,fid,scale_taij(1:ktaij_out),
     &     'scale_taij(ktaij)')
      call defvar(grid,fid,sname_taij(1:ktaij_out),
     &     'sname_taij(sname_strlen,ktaij)')
      call defvar(grid,fid,hemis_taij,'hemis_taij(one,shnhgm,ktaij)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_taij','reduction','sum')
      call defvar_cdl(grid,fid,cdl_taij,
     &     'cdl_taij(cdl_strlen,kcdl_taij)')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_taij_latlon,
     &     'cdl_taij_latlon(cdl_strlen,kcdl_taij_latlon)')
#endif

      call write_attr(grid,fid,'taijl','reduction','sum')
      call write_attr(grid,fid,'taijl','split_dim',4)
      call defvar(grid,fid,ia_taijl(1:ktaijl_out),
     &     'ia_taijl(ktaijl)')
      call defvar(grid,fid,denom_taijl(1:ktaijl_out),
     &     'denom_taijl(ktaijl)')
      call defvar(grid,fid,scale_taijl(1:ktaijl_out),
     &     'scale_taijl(ktaijl)')
      call defvar(grid,fid,sname_taijl(1:ktaijl_out),
     &     'sname_taijl(sname_strlen,ktaijl)')
      call defvar_cdl(grid,fid,cdl_taijl,
     &     'cdl_taijl(cdl_strlen,kcdl_taijl)')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_taijl_latlon,
     &     'cdl_taijl_latlon(cdl_strlen,kcdl_taijl_latlon)')
#endif

      call write_attr(grid,fid,'tajl','reduction','sum')
      call write_attr(grid,fid,'tajl','split_dim',3)
      call defvar(grid,fid,ia_tajl(1:ktajl_out),
     &     'ia_tajl(ktajl)')
      call defvar(grid,fid,denom_tajl(1:ktajl_out),
     &     'denom_tajl(ktajl)')
      call defvar(grid,fid,scale_tajl(1:ktajl_out),
     &     'scale_tajl(ktajl)')
      call defvar(grid,fid,sname_tajl(1:ktajl_out),
     &     'sname_tajl(sname_strlen,ktajl)')
      call defvar_cdl(grid,fid,cdl_tajl,
     &     'cdl_tajl(cdl_strlen,kcdl_tajl)')
      call defvar(grid,fid,hemis_tajl,'hemis_tajl(shnhgm,lm,ktajl)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_tajl','reduction','sum')
      call defvar(grid,fid,vmean_tajl,
     &     'vmean_tajl(jm_budg_plus3,one,ktajl)',r4_on_disk=.true.)
      call write_attr(grid,fid,'vmean_tajl','reduction','sum')

      call def_meta_tcons(fid)

      return
      end subroutine def_meta_trdiag

      subroutine write_meta_trdiag(fid)
!@sum  write_meta_trdiag write tracer accumulation metadata to file
!@auth M. Kelley
      use trdiag_com
      use pario, only : write_dist_data,write_data
      use domain_decomp_atm, only : grid
      use cdl_mod, only : write_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_data(grid,fid,'hemis_taij',hemis_taij)
      call write_data(grid,fid,'ia_taij',ia_taij(1:ktaij_out))
      call write_data(grid,fid,'denom_taij',denom_taij(1:ktaij_out))
      call write_data(grid,fid,'scale_taij',scale_taij(1:ktaij_out))
      call write_data(grid,fid,'sname_taij',sname_taij(1:ktaij_out))
      call write_cdl(grid,fid,'cdl_taij',cdl_taij)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_taij_latlon',cdl_taij_latlon)
#endif

      call write_data(grid,fid,'ia_taijl',ia_taijl(1:ktaijl_out))
      call write_data(grid,fid,'denom_taijl',denom_taijl(1:ktaijl_out))
      call write_data(grid,fid,'scale_taijl',scale_taijl(1:ktaijl_out))
      call write_data(grid,fid,'sname_taijl',sname_taijl(1:ktaijl_out))
      call write_cdl(grid,fid,'cdl_taijl',cdl_taijl)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_taijl_latlon',cdl_taijl_latlon)
#endif

      call write_data(grid,fid,'hemis_tajl',hemis_tajl)
      call write_data(grid,fid,'vmean_tajl',vmean_tajl)
      call write_data(grid,fid,'ia_tajl',ia_tajl(1:ktajl_out))
      call write_data(grid,fid,'denom_tajl',denom_tajl(1:ktajl_out))
      call write_data(grid,fid,'scale_tajl',scale_tajl(1:ktajl_out))
      call write_data(grid,fid,'sname_tajl',sname_tajl(1:ktajl_out))
      call write_cdl(grid,fid,'cdl_tajl',cdl_tajl)

      call write_meta_tcons(fid)

      return
      end subroutine write_meta_trdiag

#endif /* TRACERS_ON */

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)

      subroutine def_rsf_tcons(fid,r4_on_disk)
!@sum  def_rsf_tcons defines tracer diag array structure in restart+acc files
!@auth M. Kelley
!@ver  beta
      use trdiag_com, only :
     &     TCONSRV,
     &     TCONSRV_out
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid           !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
      if(r4_on_disk) then ! acc file
        call defvar(grid,fid,tconsrv_out,
     &       'tconsrv(jm_budg,ktcon)',r4_on_disk=.true.)
      else
        call defvar(grid,fid,tconsrv,
     &       'tconsrv(jm_budg,ktcon,ntmxcon)',r4_on_disk=r4_on_disk)
      endif
      return
      end subroutine def_rsf_tcons

      subroutine new_io_tcons(fid,iaction)
!@sum  new_io_tcons read/write tconsrv arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite,iowrite_single
      use trdiag_com, only :
     &     TCONSRV,
     &     TCONSRV_out
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite_single)     ! output to acc file
        call write_data(grid,fid,'tconsrv',tconsrv_out)
      case (iowrite)            ! output to restart file
        call gather_zonal_tcons
        call write_data(grid,fid,'tconsrv',tconsrv)
      case (ioread)            ! input from restart file
        call read_data(grid,fid,'tconsrv',tconsrv)
        call scatter_zonal_tcons
      end select
      return
      end subroutine new_io_tcons

      subroutine def_meta_tcons(fid)
!@sum  def_meta_trdiag defines tconsrv metadata in acc files
!@auth M. Kelley
!@ver  beta
      use trdiag_com
      use pario, only : defvar,write_attr
      use domain_decomp_atm, only : grid
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_attr(grid,fid,'tconsrv','reduction','sum')
      call write_attr(grid,fid,'tconsrv','split_dim',2)
      call defvar(grid,fid,hemis_tconsrv,'hemis_tconsrv(shnhgm,ktcon)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_tconsrv','reduction','sum')
      call defvar(grid,fid,ia_tcon_out,'ia_tconsrv(ktcon)')
      call defvar(grid,fid,scale_tcon_out,'scale_tconsrv(ktcon)')
      call defvar(grid,fid,sname_tconsrv_out,
     &     'sname_tconsrv(sname_strlen,ktcon)')
      call defvar_cdl(grid,fid,cdl_tconsrv,
     &     'cdl_tconsrv(cdl_strlen,kcdl_tconsrv)')

      return
      end subroutine def_meta_tcons

      subroutine write_meta_tcons(fid)
!@sum  write_meta_tcons write tconsrv accumulation metadata to file
!@auth M. Kelley
      use trdiag_com
      use pario, only : write_dist_data,write_data
      use domain_decomp_atm, only : grid
      use cdl_mod, only : write_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_data(grid,fid,'hemis_tconsrv',hemis_tconsrv)
      call write_data(grid,fid,'ia_tconsrv',ia_tcon_out)
      call write_data(grid,fid,'scale_tconsrv',scale_tcon_out)
      call write_data(grid,fid,'sname_tconsrv',sname_tconsrv_out)
      call write_cdl(grid,fid,'cdl_tconsrv',cdl_tconsrv)

      return
      end subroutine write_meta_tcons

#endif /* TRACERS_ON or TRACERS_OCEAN */

      SUBROUTINE ALLOC_TRDIAG_COM
      USE DIAG_COM, only : jm_budg
      USE TRDIAG_COM
      USE DOMAIN_DECOMP_ATM, only : getDomainBounds, AM_I_ROOT, GRID
      use diag_zonal, only : get_alloc_bounds
      use fluxes, only : atmocn
      implicit none
      INTEGER :: J_0H,J_1H, I_0H,I_1H
      INTEGER :: status
      integer :: j_0budg,j_1budg
      integer :: img, jmg
      
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO

      call get_alloc_bounds(grid,
     &     j_strt_budg=j_0budg,j_stop_budg=j_1budg)

#ifdef TRACERS_WATER
      ALLOCATE ( TRP_acc(ntm,I_0H:I_1H,J_0H:J_1H),stat=status)
      ALLOCATE ( TRE_acc(ntm,I_0H:I_1H,J_0H:J_1H),stat=status)
#endif
#ifdef TRACERS_ON 
      ALLOCATE(trcsurf(I_0H:I_1H,J_0H:J_1H,Ntm),stat=status)
      ALLOCATE(trcSurfByVol(I_0H:I_1H,J_0H:J_1H,Ntm),stat=status)
      ALLOCATE(trcSurfMixR_acc(I_0H:I_1H,J_0H:J_1H,Ntm),stat=status)
      ALLOCATE(trcSurfByVol_acc(I_0H:I_1H,J_0H:J_1H,Ntm),stat=status)
      ALLOCATE ( TAIJLN_loc(I_0H:I_1H,J_0H:J_1H,LM,ntm), stat=status )
      ALLOCATE ( TAIJLS_loc(I_0H:I_1H,J_0H:J_1H,LM,ktaijl), stat=status)
      ALLOCATE ( TAIJN_loc( I_0H:I_1H,J_0H:J_1H,ktaij,ntm),stat=status )
      ALLOCATE ( TAIJS_loc( I_0H:I_1H,J_0H:J_1H,ktaijs   ),stat=status )
      ALLOCATE ( TAJLN_loc(  J_0BUDG:J_1BUDG,LM,ktajlx,ntm),stat=status)
      ALLOCATE ( TAJLS_loc(  J_0BUDG:J_1BUDG,LM,ktajls    ),stat=status)
      ALLOCATE ( TSCF3d_loc( I_0H:I_1H,J_0H:J_1H,lm,ntm),stat=status )
      if(am_i_root()) then
         img = IM
         jmg = JM
      else
         img = 1
         jmg = 1
      end if
      ALLOCATE ( TAIJLN(img,jmg,LM,ntm), stat=status )
      ALLOCATE ( TSCF3d(img,jmg,LM,ntm), stat=status )
      ALLOCATE ( TAIJLS(img,jmg,LM,ktaijl), stat=status )
      ALLOCATE ( TAIJN( img,jmg,ktaij,ntm), stat=status )
      ALLOCATE ( TAIJS( img,jmg,ktaijs   ), stat=status )

      ALLOCATE ( TAJLN(JM_BUDG,LM,ktajlx,ntm), stat=status )
      ALLOCATE ( TAJLS(JM_BUDG,LM,ktajls    ), stat=status )

      allocate(jls_incloud(2,ntm))
      jls_incloud = 0
      allocate(jls_source(ntsurfsrcmax,NTM))
      jls_source = 0
      allocate(jls_isrc(ntsurfsrcmax,NTM))
      jls_isrc = 0
      allocate(jls_3Dsource(nt3Dsrcmax,NTM))
      jls_3Dsource = 0
      allocate(jls_decay(NTM))
      jls_decay = 0
      allocate(jls_grav(NTM))
      jls_grav = 0
      allocate(jls_prec(2,NTM))
      jls_prec = 0
      allocate(jls_trdpmc(MaxDMc,NTM))
      jls_trdpmc = 0
      allocate(jls_trdpls(MaxDLs,NTM))
      jls_trdpls = 0
      allocate(jls_wet(NTM))
      jls_wet = 0


#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      allocate(to_per_mil(NTM))
      to_per_mil = 0
#endif
      allocate(MMR_to_VMR(NTM))

      allocate(scale_tij(ktaij,ntm))
      allocate(sname_tij(ktaij,ntm))
      allocate(units_tij(ktaij,ntm))
      allocate(dname_tij(ktaij,NTM))
      dname_tij=''
      allocate(denom_tij(ktaij,NTM))
      denom_tij = 0
      allocate(lname_tij(ktaij,NTM))
      lname_tij = 'unused'

      allocate(sname_ijt(NTM))
      allocate(units_ijt(NTM))
      allocate(lname_ijt(NTM))
      lname_ijt(:) = 'unused'

      allocate(scale_ijt(NTM))
      allocate(ijtc_power(NTM))
      ijtc_power = 0
      allocate(ijtm_power(NTM))
      ijtm_power = 0
      allocate(ijts_source(ntsurfsrcmax,ntm))
      ijts_source = 0
      allocate(ijts_isrc(ntsurfsrcmax,ntm))
      ijts_isrc = 0
      allocate(ijts_aq(ntm))
      ijts_aq = 0
      allocate(ijts_tau(3,ntm))
      ijts_tau = 0
      allocate(ijts_tausub(3,Ntm,MaxSubCl))
      ijts_tausub = 0
      allocate(ijts_sqex(3,6,Ntm))
      ijts_sqex = 0
      allocate(ijts_sqexsub(3,6,Ntm,MaxSubCl))
      ijts_sqexsub = 0
      allocate(ijts_sqsc(3,6,Ntm))
      ijts_sqsc = 0
      allocate(ijts_sqscsub(3,6,Ntm,MaxSubCl))
      ijts_sqscsub = 0
      allocate(ijts_sqcb(3,6,Ntm))
      ijts_sqcb = 0
      allocate(ijts_sqcbsub(3,6,Ntm,MaxSubCl))
      ijts_sqcbsub = 0
      allocate(ijts_fc(8,ntm))
      ijts_fc = 0
      allocate(ijts_fcsub(8,Ntm,MaxSubCl))
      ijts_fcsub = 0
      allocate(ijts_gasex(3,Ntm))
      ijts_gasex = 0
      allocate(ijts_AMPe(Ntm))
      ijts_AMPe = 0
      allocate(ijts_AMPp(7,Ntm))
      ijts_AMPp = 0
#ifdef TRACERS_TOMAS
      allocate(ijts_TOMAS(7,Ntm))
      ijts_TOMAS = 0
      allocate(ijts_subcoag(Ntm))
      ijts_subcoag = 0
#endif
      allocate(ijts_trdpmc(MaxDMc,Ntm))
      ijts_trdpmc = 0
      allocate(ijts_trdpls(MaxDLs,Ntm))
      ijts_trdpls = 0
      allocate(ijts_wet(Ntm))
      ijts_wet = 0
      allocate(ijts_3Dsource(nt3Dsrcmax,ntm))
      ijts_3Dsource = 0
#endif  /* TRACERS_ON */

#ifdef TRACERS_AMP
      allocate(ijlt_AMPm(3,ntm))
      ijlt_AMPm = 0
#endif 
#ifdef SAVE_AEROSOL_3DMASS_FOR_NINT
      allocate(ijlt_3Dmass(ntm))
      ijlt_3Dmass = 0
#endif
      allocate(sname_jln(ktajlx,ntm))
      allocate(lname_jln(ktajlx,ntm))
      lname_jln = 'unused'
      allocate(units_jln(ktajlx,ntm))
      allocate(scale_jln(ntm))


#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)

      ntmxcon = ntm + maxntmocn
      allocate(TCONSRV(JM_BUDG,ktcon,ntmxcon))
      allocate(TCONSRV_loc(JM_BUDG,ktcon,ntmxcon))

      allocate(SCALE_TCON(ktcon,ntmxcon))
      allocate(TITLE_TCON(ktcon,ntmxcon))
      allocate(IA_TCON(ktcon,ntmxcon))
      IA_TCON = 1
      allocate(NSUM_TCON(ktcon,ntmxcon))
      NSUM_TCON = 0
      allocate(NOFMT(ktcon,ntmxcon))
      NOFMT = 0

#ifdef TRACERS_OCEAN
      atmocn%tconsrv => tconsrv_loc
      atmocn%nofmt => nofmt
#endif

      allocate(kt_power_inst(ntm))
      kt_power_inst = 0
      allocate(kt_power_change(ntm))
      kt_power_change = 0
      allocate(name_tconsrv(ktcon,ntmxcon))
      name_tconsrv='unused'
      allocate(units_tconsrv(ktcon,ntmxcon))
      allocate(lname_tconsrv(ktcon,ntmxcon))

      allocate(SCALE_INST(ntmxcon))
      allocate(SCALE_CHANGE(ntmxcon))
      allocate(itcon_surf(ntsurfsrcmax,ntmxcon))
      itcon_surf = 0
      allocate(itcon_3Dsrc(nt3Dsrcmax,ntmxcon))
      itcon_3Dsrc = 0
      allocate(itcon_decay(ntmxcon))
      itcon_decay = 0
      allocate(itcon_mc(ntmxcon))
      itcon_mc = 0
      allocate(itcon_ss(ntmxcon))
      itcon_ss = 0
      allocate(itcon_amp(7,ntmxcon))
      itcon_amp = 0
      allocate(itcon_ampm(3,ntmxcon))
      itcon_ampm = 0
      allocate(itcon_ampe(ntmxcon))
      itcon_ampe = 0
      allocate(itcon_dd(ntmxcon,2))
      itcon_dd = 0
      allocate(itcon_wt(ntmxcon))
      itcon_wt = 0

#ifdef TRACERS_TOMAS
      allocate(itcon_TOMAS(7,ntmxcon))
      itcon_TOMAS = 0
      allocate(itcon_subcoag(ntmxcon))
      itcon_subcoag = 0
#endif
#endif  /* TRACERS_ON  or  TRACERS_OCEAN */

      ktaijl_ = (ntm + ktaijl
#ifdef TRACERS_SPECIAL_O18
     *     + 2        ! include dexcess + D17O diags
#endif
     &     )*3/2  ! make 50% larger for denoms and extra specials


      allocate(ir_taijl(ktaijl_))
      ir_taijl = 0
      allocate(ia_taijl(ktaijl_))
      ia_taijl = 0
      allocate(denom_taijl(ktaijl_))
      denom_taijl = 0
      allocate(lname_taijl(ktaijl_))
      allocate(sname_taijl(ktaijl_))
      allocate(units_taijl(ktaijl_))
      allocate(scale_taijl(ktaijl_))

      ktajl_ = (ktajlx*ntm+ktajls
#ifdef TRACERS_SPECIAL_O18
     &     + 4                  ! include dexcess + D17O diags
#endif
#ifdef TRACERS_COSMO
     &     + 2                  ! beryllium ratios
#endif
     &     )*3/2  ! make 50% larger for denoms and extra specials
      allocate(pow_tajl(ktajl_))
      pow_tajl = 0
      allocate(ia_tajl(ktajl_))
      ia_tajl = 0
      allocate(denom_tajl(ktajl_))
      denom_tajl = 0
      allocate(jgrid_tajl(ktajl_))
      jgrid_tajl = 0 
      allocate(lgrid_tajl(ktajl_))
      lgrid_tajl = 0
      allocate(ltop_tajl(ktajl_))
      ltop_tajl = 0
      allocate(lname_tajl(ktajl_))
      allocate(sname_tajl(ktajl_))
      allocate(units_tajl(ktajl_))
      allocate(scale_tajl(ktajl_))

! ktaij_, ktaijl_, ktajl_ are larger than necessary.  These arrays
! will be reallocated to the proper sizes later.
      ALLOCATE ( TAIJ_out( I_0H:I_1H,J_0H:J_1H,ktaij_),stat=status )
      ALLOCATE ( TAIJL_out( I_0H:I_1H,J_0H:J_1H,LM,ktaijl_),stat=status)
      if(am_i_root()) then
        ALLOCATE ( TAJL_out(JM_BUDG,LM,ktajl_), stat=status )
      endif

      RETURN
      END SUBROUTINE ALLOC_TRDIAG_COM

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** reset/gather/scatter routines for tconsrv special case
C**** can be called from both ocean tracer or atm tracer code
      subroutine reset_tcons
      USE TRDIAG_COM, only: TCONSRV_loc, TCONSRV
      USE DOMAIN_DECOMP_ATM, only : am_i_root
      implicit none
      TCONSRV_loc=0.
      if(am_i_root()) TCONSRV=0.
      return
      end subroutine reset_tcons

      subroutine gather_zonal_tcons
      USE DIAG_COM, only : ia_inst,jm_budg
      USE TRDIAG_COM, only: TCONSRV_loc, TCONSRV, ia_tcon,ntmxcon,ktcon
      USE DOMAIN_DECOMP_ATM, only : GRID,am_i_root,sumxpe
      implicit none
      real*8, dimension(:,:,:), allocatable :: tconsrv_sv
      integer :: n,k

c TCONSRV is a mixture of accumulations and instantaneous values, hence the
c complicated logic
! set inst qtys to zero in the "global" array before summing over PEs
      if(am_i_root()) then
        allocate(tconsrv_sv(jm_budg,ktcon,ntmxcon)); tconsrv_sv=tconsrv
        do n=1,ntmxcon
          do k=1,ktcon
            if(ia_tcon(k,n).eq.ia_inst) tconsrv(:,k,n)=0.
          enddo
        enddo
      endif

      call sumxpe (TCONSRV_loc,TCONSRV, increment=.true. )

      do n=1,ntmxcon ! keep instantaneous values
        do k=1,ktcon
          if(ia_tcon(k,n).ne.ia_inst) tconsrv_loc(:,k,n)=0. 
        enddo
      enddo

      if(am_i_root()) then
        do n=1,ntmxcon
          do k=1,ktcon
            if(ia_tcon(k,n).eq.ia_inst .and. all(tconsrv(:,k,n)==0.))
     &         tconsrv(:,k,n) = tconsrv_sv(:,k,n)
          enddo
        enddo
        deallocate(tconsrv_sv)
      endif

      return
      end subroutine gather_zonal_tcons

      subroutine scatter_zonal_tcons
      USE TRDIAG_COM, only: TCONSRV_loc, TCONSRV
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DIAG_ZONAL, only : unpack_lc
      implicit none
c      call unpack_lc   (grid, TCONSRV, TCONSRV_loc)
      tconsrv_loc = 0
      return
      end subroutine scatter_zonal_tcons
#endif

#ifdef TRACERS_ON
      subroutine reset_trdiag
      USE TRDIAG_COM, only: TAIJLN_loc, TAIJLS_loc, TAIJN_loc, 
     *     TAIJS_loc, TAJLN_loc, TAJLS_loc, TAJLN, TAJLS
      USE DOMAIN_DECOMP_ATM, only : am_i_root
      implicit none

       TAJLN_loc=0. ; TAJLS_loc=0. 
      TAIJLN_loc=0. ; TAIJLS_loc=0. ; TAIJN_loc=0. ; TAIJS_loc=0.
      if(am_i_root()) then
        TAJLN=0. ; TAJLS=0. 
      endif
      call reset_tcons

      return
      end subroutine reset_trdiag

      subroutine gather_trdiag
      USE TRDIAG_COM, only : TAIJLN, TAIJLN_loc, TAIJLS, TAIJLS_loc,
     *     TAIJN, TAIJN_loc,TAIJS, TAIJS_loc
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, ONLY : PACK_DATA
      implicit none

      CALL PACK_DATA (GRID, TAIJLN_loc, TAIJLN)
      CALL PACK_DATA (GRID, TAIJLS_loc, TAIJLS)
      CALL PACK_DATA (GRID, TAIJN_loc , TAIJN)
      CALL PACK_DATA (GRID, TAIJS_loc , TAIJS)
      call gather_zonal_trdiag
      call gather_zonal_tcons

      return
      end subroutine gather_trdiag

      subroutine gather_zonal_trdiag
      USE TRDIAG_COM, only : TAJLN , TAJLN_loc, TAJLS, TAJLS_loc
      USE DOMAIN_DECOMP_ATM, ONLY : GRID
      USE DIAG_ZONAL, only : pack_lc
      implicit none
      call pack_lc   (grid, TAJLN_loc,  TAJLN )
      call pack_lc   (grid, TAJLS_loc,  TAJLS )
      return
      end subroutine gather_zonal_trdiag

      subroutine scatter_trdiag
      USE TRDIAG_COM, only : TAIJLN, TAIJLN_loc, TAIJLS, TAIJLS_loc,
     *     TAIJN, TAIJN_loc,TAIJS, TAIJS_loc
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, ONLY : UNPACK_DATA
      implicit none
      CALL UNPACK_DATA (GRID, TAIJLN, TAIJLN_loc)
      CALL UNPACK_DATA (GRID, TAIJLS, TAIJLS_loc)
      CALL UNPACK_DATA (GRID, TAIJN , TAIJN_loc)
      CALL UNPACK_DATA (GRID, TAIJS , TAIJS_loc)
      call scatter_zonal_trdiag
      call scatter_zonal_tcons
      return
      end subroutine scatter_trdiag

      subroutine scatter_zonal_trdiag
      USE TRDIAG_COM, only : TAJLN , TAJLN_loc, TAJLS, TAJLS_loc
      USE DOMAIN_DECOMP_ATM, ONLY : GRID
      USE DIAG_ZONAL, only : unpack_lc
      implicit none
      call unpack_lc (grid, TAJLN,  TAJLN_loc )
      call unpack_lc (grid, TAJLS,  TAJLS_loc )
      return
      end subroutine scatter_zonal_trdiag

C**** routines for accumulating zonal mean diags (lat/lon grid)

      SUBROUTINE INC_TAJLS(I,J,L,TJL_INDEX,ACC)
!@sum inc_tajl adds ACC located at atmospheric gridpoint I,J,L
!@+   to the latitude-height zonal sum TAJLS(TJL_INDEX).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajls=>tajls_loc
      USE DIAG_COM, only : wtbudg,j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC
C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      TAJLS(J_BUDG(I,J),L,TJL_INDEX) = TAJLS(J_BUDG(I,J),L,TJL_INDEX) +
     *     ACC!*wtbudg(I,J) ! cannot use wtbudg!=1 b/c kg units for tracers

      RETURN
      END SUBROUTINE INC_TAJLS

      SUBROUTINE INC_TAJLS2(I,J,L,TJL_INDEX,ACC)
!@sum inc_tajl adds ACC located at atmospheric gridpoint I,J,L
!@+   to the latitude-height zonal sum TAJLS(TJL_INDEX).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajls=>tajls_loc
      USE DIAG_COM, only : wtbudg2,j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC
C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      TAJLS(J_BUDG(I,J),L,TJL_INDEX) = TAJLS(J_BUDG(I,J),L,TJL_INDEX) +
     *     ACC*wtbudg2(I,J)

      RETURN
      END SUBROUTINE INC_TAJLS2

      SUBROUTINE INC_TAJLS_COLUMN(I,J,L1,L2,NL,TJL_INDEX,ACC)
!@sum inc_tajl_column adds ACC(L1:L2) located at atmospheric gridpoint I,J
!@+   to the latitude-height zonal sums TAJLS(:,L1:L2,TJL_INDEX).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajls=>tajls_loc
      USE DIAG_COM, only : wtbudg,j_budg
      IMPLICIT NONE
!@var I,J,L1,L2 atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L1,L2,NL
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC(NL)
      INTEGER :: JB,L
C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      JB = J_BUDG(I,J)
      DO L=L1,L2
        TAJLS(JB,L,TJL_INDEX) = TAJLS(JB,L,TJL_INDEX) +
     *     ACC(L)!*wtbudg(I,J) ! cannot use wtbudg!=1 b/c kg units for tracers
      ENDDO
      RETURN
      END SUBROUTINE INC_TAJLS_COLUMN

      SUBROUTINE INC_TAJLS2_COLUMN(I,J,L1,L2,NL,TJL_INDEX,ACC)
!@sum inc_tajl_column adds ACC(L1:L2) located at atmospheric gridpoint I,J
!@+   to the latitude-height zonal sums TAJLS(:,L1:L2,TJL_INDEX).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajls=>tajls_loc
      USE DIAG_COM, only : wtbudg2,j_budg
      IMPLICIT NONE
!@var I,J,L1,L2 atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L1,L2,NL
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC(NL)
      INTEGER :: JB,L
C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      JB = J_BUDG(I,J)
      DO L=L1,L2
        TAJLS(JB,L,TJL_INDEX) = TAJLS(JB,L,TJL_INDEX) +
     *     ACC(L)*wtbudg2(I,J)
      ENDDO
      RETURN
      END SUBROUTINE INC_TAJLS2_COLUMN

      SUBROUTINE INC_TAJLN(I,J,L,TJL_INDEX,N,ACC)
!@sum inc_tajln adds ACC located at atmospheric gridpoint I,J,L
!@+   and tracer n to the latitude-height zonal sum TAJLN(TJL_INDEX,N).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajln=>tajln_loc
      USE DIAG_COM, only : wtbudg,j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices, N tracer no. for the accumulation
      INTEGER, INTENT(IN) :: I,J,L,N
!@var TJL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC
      TAJLN(J_BUDG(I,J),L,TJL_INDEX,N)=TAJLN(J_BUDG(I,J),L,TJL_INDEX,N)
     *     + ACC!*wtbudg(I,J) ! cannot use wtbudg!=1 b/c kg units for tracers

      RETURN
      END SUBROUTINE INC_TAJLN

      SUBROUTINE INC_TAJLN_COLUMN(I,J,L1,L2,NL,TJL_INDEX,N,ACC)
!@sum inc_tajln adds ACC located at atmospheric gridpoint I,J,L
!@+   and tracer n to the latitude-height zonal sum TAJLN(TJL_INDEX,N).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajln=>tajln_loc
      USE DIAG_COM, only : wtbudg,j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices, N tracer no. for the accumulation
      INTEGER, INTENT(IN) :: I,J,L1,L2,NL,N
!@var TJL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC(NL)
      INTEGER :: JB,L
      JB = J_BUDG(I,J)
      DO L=L1,L2
        TAJLN(JB,L,TJL_INDEX,N)=TAJLN(JB,L,TJL_INDEX,N) +
     *    ACC(L) !*wtbudg(I,J) ! cannot use wtbudg!=1 b/c kg units for tracers
      ENDDO
      RETURN
      END SUBROUTINE INC_TAJLN_COLUMN
#endif   /* TRACERS_ON */
