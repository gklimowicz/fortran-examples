#include "rundeck_opts.h"

! Remove the following line and directive when the problems are fixed.
! This directive is simply to permit running a few timesteps.

!@sum  TRACERS_DRV: tracer-dependent routines for air/water mass
!@+    and ocean tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers:
!@+        Diagnostic specs: init_tracer
!@+        Tracer initialisation + sources: tracer_ic, set_tracer_source
!@+        Entry points: daily_tracer
!@auth Jean Lerner/Gavin Schmidt
!=======================================================================
      integer function get_src_index(n)
!@var get_src_index If an emission file contains information for more than one
!@+                 tracer, first read tracer n_XXX, then set src_index=n_XXX.
!@+                 Note the order! Notable case is SO2/SO4
!@auth Kostas Tsigaridis
      use OldTracer_mod, only: trname
      use TRACER_COM, only: n_SO2
      implicit none
!@var n index of current tracer whose emissions index is being seeked
      integer, intent(in) :: n

      select case (trname(n))
        case ('SO4', 'M_ACC_SU', 'ASO4__01')
          get_src_index=n_SO2
        case ('M_AKK_SU')
          get_src_index=n_SO2
        case default
          get_src_index=n
      end select

      end function get_src_index
!=======================================================================
      real*8 function get_src_fact(n,OA_not_OC)
!@var src_fact Factor to multiply aerosol emissions with. Default is 1.
!@+            Notable exceptions are SO2/SO4, where one file is being read
!@+            and distributed to both tracers,and organics, where emissions
!@+            of C are multiplied with OM/OC, and VBS tracers.
!@auth Kostas Tsigaridis
      use OldTracer_mod, only: trname
      use OldTracer_mod, only: tr_mm
      use OldTracer_mod, only: om2oc
      use TRACER_COM, only: n_M_AKK_SU
#ifdef TRACERS_AEROSOLS_VBS
      use aerosol_sources, only: VBSemifact
      use tracers_vbs, only: vbs_tr
#endif  /* TRACERS_AEROSOLS_VBS */
      implicit none
!@var n index of current tracer whose emissions factor is being seeked
!@var OA_not_OC true if the tracer has emissions in OA units, not OC
!@var so4_fraction mole fraction of so2 to be emitted as so4
      integer, intent(in) :: n
      logical, intent(in), optional :: OA_not_OC
      real*8, parameter :: so4_fraction=0.025d0
      real*8 :: akk_fraction
      logical OA_not_OC_local
      integer get_src_index

      if (n_M_AKK_SU>0) then
        akk_fraction=0.01d0
      else
        akk_fraction=0.d0
      endif

      OA_not_OC_local=.false.
      if (present(OA_not_OC)) OA_not_OC_local=OA_not_OC

      select case (trname(n))
        case ('SO2')
          get_src_fact=1.d0-so4_fraction
        case ('SO4', 'M_ACC_SU', 'ASO4__01')
          get_src_fact=so4_fraction*tr_mm(n)/tr_mm(get_src_index(n))*
     &                 (1.d0-akk_fraction)
        case ('M_AKK_SU')
          get_src_fact=so4_fraction*tr_mm(n)/tr_mm(get_src_index(n))*
     &                 akk_fraction
        case ('OCII', 'OCIA', 'OCB', 'M_OCC_OC', 'M_BOC_OC', 'AOCOB_01')
          get_src_fact=1.d0
          if (.not.OA_not_OC_local) get_src_fact=get_src_fact*om2oc(n)
#ifdef TRACERS_AEROSOLS_VBS
        case ('vbsAm2', 'vbsAm1', 'vbsAz', 'vbsAp1', 'vbsAp2',
     &        'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
          get_src_fact=VBSemifact(vbs_tr%iaerinv(n))
          if (.not.OA_not_OC_local) get_src_fact=get_src_fact*om2oc(n)
#endif  /* TRACERS_AEROSOLS_VBS */
        case default
          get_src_fact=1.d0
      end select

      end function get_src_fact
!=======================================================================
      integer function tr_con_diag(vconpts, vqcon, vqsum)
!@sum tr_con_diag populate tracer conservation diagnostics
!@auth Kostas Tsigaridis
      use TRDIAG_COM, only: ntcons,conpts,npts_common,qcon,qsum
      implicit none
!@var vconpts local value of conpts
!@var vqcon local value of qcon
!@var vqsum local value of qsum
!@var g index to be assigned to the current diagnostic
!@var i local loop index
      character(len=*), intent(in) :: vconpts
      logical, intent(in), optional :: vqcon, vqsum
      integer :: g,i

      g=0
      do i=1,ntcons ! brute force, but only happens during initialization
        if (trim(conpts(i))=='') then
          g=npts_common+i
          exit
        endif
      enddo
      if (g==0) call stop_model('ntcons too small',255)

      tr_con_diag=g
      conpts(g-npts_common)=trim(vconpts)
      if (present(vqcon)) qcon(g)=vqcon
      if (present(vqsum)) qsum(g)=vqsum

      end function tr_con_diag
!=======================================================================
      subroutine init_tracer_cons_diag
!@sum init_tracer_cons_diag Initialize tracer conservation diagnostics
!@auth Gavin Schmidt
      use AbstractAttribute_mod, only: AbstractAttribute
      use Attributes_mod, only: assignment(=)
      use Attributes_mod, only: toPointer
      use Tracer_mod
      use TracerBundle_mod
      use TracerHashMap_mod
      use AttributeDictionary_mod, only: assignment(=)
      use TracerSurfaceSource_mod, only: TracerSurfaceSource
      USE TRACER_COM, only: ntm
      USE TRACER_COM, only: noverwrite
      USE TRACER_COM, only: nvolcanic
      USE TRACER_COM, only: nother
      use OldTracer_mod, only: ntm_power, dowetdep, dodrydep
      use OldTracer_mod, only: tr_wd_type, nPart
      use OldTracer_mod, only: nBBsources,trname,do_fire
      use TRACER_COM, only: nchemloss
      use TRACER_COM, only: nchemistry
      use TRACER_COM, only: nbiomass
      use TRACER_COM, only: naircraft, do_aircraft
      use TRACER_COM, only: ntsurfsrc
      use TRACER_COM, only: tracers
      use TRACER_COM, only: n_SO2
      use Tracer_mod, only: Tracer
      use DIAG_COM, only: conpt0,npts
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: n_AH2O, n_AECOB, n_AOCOB, n_ANUM
      use TRACER_COM, only: N_AECOB, n_ASO4, xk
#endif
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif /* TRACERS_ON */
      USE FLUXES, only : atmocn
      implicit none
      character*20 sum_unit(NTM),inst_unit(NTM)   ! for conservation
      character*50 :: unit_string
      character*10 :: conpt(npts)
#ifdef TRACERS_ON
      logical :: T=.TRUE. , F=.FALSE.
      logical :: Qf
      integer n,n_src,kk
      integer, pointer :: index=> null()
      class (AbstractAttribute), pointer :: pa
      class (Tracer), pointer :: pTracer,pTracerSrc
      type (TracerSurfaceSource), pointer :: sources(:)
#endif
      type (TracerIterator) :: iter
      interface
        integer function tr_con_diag(vconpts, vqcon, vqsum)
          character(len=*), intent(in) :: vconpts
          logical, intent(in), optional :: vqcon, vqsum
        end function tr_con_diag
      end interface

#ifdef TRACERS_ON

C**** To add a new conservation diagnostic:
C****       Set up a QCON, and call SET_TCON to allocate array numbers,
C****       set up scales, titles, etc.
C**** QCON denotes when the conservation diags should be accumulated
C**** QSUM says whether that diag is to be used in summation (if the
C****      routine DIAGTCB is used, this must be false).
C**** 1:NPTS+1 ==> INST,  DYN,   COND,   RAD,   PREC,   LAND,  SURF,
C****            FILTER,STRDG/OCEAN, DAILY, OCEAN1, OCEAN2,
C**** First 12 (npts_common) are standard for all tracers and GCM
C**** Later indices are configurable - you provide title and itcon
C**** index (which is used wherever you want to check point)
C**** For example, separate Moist convection/Large scale condensation
!      itcon_mc(n)=tr_con_diag('MOIST CONV',T)
!      itcon_ss(n)=tr_con_diag('LS COND',T)

#ifdef CUBED_SPHERE
      Qf = .false.  ! no SLP filter
#else
      Qf = .true.   ! SLP filter on
#endif

      qcon(1:npts_common)=(/T,                                !instant. (1)
     *                      T, T, F, F, T, T,Qf, T, F, F, F/) !2-12 (npts)
      qcon(npts_common+1:npts_common+ntcons)=F                !13-ktcon-1
      qsum(1:npts_common)=(/F,                                !instant. (1)
     *                      T, T, F, F, T, T,Qf, T, F, F, F/) !2-12 (npts)
      qsum(npts_common+1:npts_common+ntcons)=F                !13-ktcon-1
C**** this allows you to configure the common check points names.
      conpt=conpt0

      do n=1,NTM
        kt_power_inst(n)   = ntm_power(n)+2
        kt_power_change(n) = ntm_power(n)-4
      end do

C**** set some defaults
      itcon_AMP(:,:)=0
      itcon_AMPe(:)=0
      itcon_AMPm(:,:)=0
      itcon_surf(:,:)=0
      itcon_3Dsrc(:,:)=0
      itcon_decay(:)=0
      itcon_wt(:)=0
#ifdef TRACERS_WATER
      itcon_mc(:)=0
      itcon_ss(:)=0
#endif
#ifdef TRACERS_DRYDEP
      itcon_dd(:,:)=0
#endif
#ifdef TRACERS_TOMAS
      itcon_TOMAS(:,:)=0
      itcon_subcoag(:)=0
#endif

      iter = tracers%begin()
      do while (iter /= tracers%last())
        pTracer => iter%value()
        index => toPointer(pTracer%getReference('index'), index)
        n = index

! handle exceptions first (e.g. SO4 emissions are listed under SO2)
        select case (trim(pTracer%getName()))
        case ('SO4',
     &        'M_AKK_SU','M_ACC_SU',
     &        'ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     &        'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     &        'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15')
          n_src = n_SO2
#ifdef TRACERS_TOMAS
        case ('AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05',
     &        'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10',
     &        'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15',
     &        'AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05',
     &        'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10',
     &        'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15')
          n_src = n_AECOB(1)
        case ('AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05',
     &        'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10',
     &        'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15',
     &        'AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05',
     &        'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10',
     &        'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15')
          n_src = n_AOCOB(1)
#endif  /* TRACERS_TOMAS */
        case default
          n_src = n
        end select
        pTracerSrc => tracers%getReference(trname(n_src))
        sources => pTracerSrc%surfaceSources

!-----
! diagnostics for all tracers, if they meet certain conditions
!-----
        conpt(8)="SRCS+SNKS"

#ifdef TRACERS_WATER
        if(dowetdep(n)) then
          itcon_mc(n)=tr_con_diag('MOIST CONV',T)
          itcon_ss(n)=tr_con_diag('LS COND',T)
        endif
#endif
#ifdef TRACERS_DRYDEP
        if(dodrydep(n)) then
          itcon_dd(n,1)=tr_con_diag('TURB DEP',T)
          if (tr_wd_type(n)==nPart) then
            itcon_dd(n,2)=tr_con_diag('GRAV SET',T)
          endif
        end if
#endif
        if(do_aircraft(n_src))then
          itcon_3Dsrc(nAircraft,n)=tr_con_diag('Aircraft src',T,T)
        endif
        if (nBBsources(n_src)>0 .or. do_fire(n_src)) then
          itcon_3Dsrc(nBiomass,n)=tr_con_diag('Biomass src',T,T)
        endif
        do kk=1,ntsurfsrc(n_src)
          itcon_surf(kk,n)=tr_con_diag(trim(sources(kk)%sourceLname),T)
        enddo

!-----
!     per-tracer diagnostics
!     Note: Do not include already registered surface sources
!-----
        select case (trim(pTracer%getName()))

        case ('CO2n')
          qcon(10) = .true.
          qsum(10) = .true.

        case ('Rn222')
          itcon_decay(n)=tr_con_diag('DECAY',T,T)

        case ('CO2')
c          itcon_surf(1,N)=tr_con_diag('FossilFuel',T)
c          itcon_surf(2,N)=tr_con_diag('Fertilization',T)
c          itcon_surf(3,N)=tr_con_diag('Forest Regrowth',T)
c          itcon_surf(4,N)=tr_con_diag('Land Use',T)
c          itcon_surf(5,N)=tr_con_diag('Ecosystem Exch',T)
c          itcon_surf(6,N)=tr_con_diag('Ocean Exch',T)

        case ('N2O')   ! two versions dependent on configuration
#ifdef TRACERS_SPECIAL_Lerner
c          itcon_surf(1,N)=tr_con_diag('Reset in L1',T)
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Strat. Chem.',T,T)
#endif
#ifdef TRACERS_SPECIAL_Shindell
          kt_power_change(n) = -14
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Chemistry',T,T)
          itcon_3Dsrc(nOverwrite,n)=tr_con_diag('Overwrite',T,T)
#endif

        case ('CFC11')
c          itcon_surf(1,N)=tr_con_diag('L1 Source',T)
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Strat. Chem.',T,T)

        case ('14CO2')
          itcon_surf(1,N)=tr_con_diag('Bombs and drift',T)


        case ('nh5')
          kt_power_change(n) = -13
          itcon_decay(n)=tr_con_diag('DECAY',T,T)

        case ('nh15')
          kt_power_change(n) = -13
          itcon_decay(n)=tr_con_diag('DECAY',T,T)

        case ('nh50')
          kt_power_change(n) = -13
          itcon_decay(n)=tr_con_diag('DECAY',T,T)

        case ('e90')
          kt_power_change(n) = -13
          itcon_decay(n)=tr_con_diag('DECAY',T,T)

        case ('aoa')
          kt_power_change(n) = -17
          itcon_3Dsrc(1,n)=tr_con_diag('L1 overwriting',T,T)

        case ('aoanh')
          kt_power_change(n) = -17
          itcon_3Dsrc(1,n)=tr_con_diag('L1 overwriting',T,T)


        case ('CH4')            ! two versions
#ifdef TRACERS_SPECIAL_Shindell
          kt_power_change(n) = -13
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Chemistry',T,T)
          itcon_3Dsrc(nOverwrite,n)=tr_con_diag('Overwrite',T,T)
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_SPECIAL_Lerner
c          itcon_surf(1,N)=tr_con_diag('Animal source',T)
c          itcon_surf(2,N)=tr_con_diag('Coal Mine source',T)
c          itcon_surf(3,N)=tr_con_diag('Gas Leak source',T)
c          itcon_surf(4,N)=tr_con_diag('Gas Vent source',T)
c          itcon_surf(5,N)=tr_con_diag('City Dump source',T)
c          itcon_surf(6,N)=tr_con_diag('Soil sink',T)
c          itcon_surf(7,N)=tr_con_diag('Termite Source',T)
c          itcon_surf(8,N)=tr_con_diag('Coal Combustion',T)
c          itcon_surf(9,N)=tr_con_diag('Ocean source',T)
c          itcon_surf(10,N)=tr_con_diag('Lake source',T)
c          itcon_surf(11,N)=tr_con_diag('Misc. Ground source',T)
c          itcon_surf(12,N)=tr_con_diag('Biomass Burning',T)
c          itcon_surf(13,N)=tr_con_diag('Rice source',T)
c          itcon_surf(14,N)=tr_con_diag('Wetlands+Tundra',T)
          itcon_3Dsrc(1,n)=tr_con_diag('Tropos. Chem.',T,T)
          itcon_3Dsrc(2,n)=tr_con_diag('Stratos. Chem.',T,T)
#endif /* TRACERS_SPECIAL_Lerner */

        case ('O3')
c          itcon_surf(1,N)=tr_con_diag('Deposition',T)
          itcon_3Dsrc(1,n)=tr_con_diag('Stratos. Chem.',T,T)
          itcon_3Dsrc(2,n)=tr_con_diag('Trop. Chem. Prod.',T,T)
          itcon_3Dsrc(3,n)=tr_con_diag('Trop. Chem. Loss',T,T)

        case ('Ox','N2O5','HNO3','H2O2','CH3OOH','HCHO','HO2NO2','PAN'
     *       ,'AlkylNit','ClOx','BrOx','HCl','HOCl','ClONO2','HBr'
#ifdef TRACERS_dCO
     *       ,'d13Calke','d13CPAR'
     *       ,'d17OPAN', 'd18OPAN', 'd13CPAN'
     *       ,'dMe17OOH', 'dMe18OOH', 'd13MeOOH'
     *       ,'dHCH17O', 'dHCH18O', 'dH13CHO'
     *       ,'dC17O', 'dC18O', 'd13CO'
#endif  /* TRACERS_dCO */
     *       ,'HOBr','BrONO2','CFC','NOx','CO','Isoprene','Alkenes'
     *       ,'Paraffin','stratOx','Terpenes') ! N2O done above
          select case (trim(pTracer%getName()))
            case ('N2O5','CH3OOH','HCHO','HO2NO2','PAN','AlkylNit','CFC'
#ifdef TRACERS_dCO
     *           ,'d13CPAR'
     *           ,'d17OPAN', 'd18OPAN', 'd13CPAN'
     *           ,'dMe17OOH', 'dMe18OOH', 'd13MeOOH'
     *           ,'dHCH17O', 'dHCH18O', 'dH13CHO'
#endif  /* TRACERS_dCO */
     *           ,'ClOx','BrOx','HCl','HOCl','ClONO2','HBr','HOBr'
     *           ,'BrONO2','NOx')
              kt_power_change(n) = -14
            case ('HNO3','H2O2','CO','Isoprene','Alkenes','Paraffin'
#ifdef TRACERS_dCO
     *           ,'d13Calke'
     *           ,'dC17O', 'dC18O', 'd13CO'
#endif  /* TRACERS_dCO */
     *           ,'Terpenes')
              kt_power_change(n) = -13
            case default
              kt_power_change(n) = -12
          end select

          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Chemistry',T,T)
          itcon_3Dsrc(nOverwrite,n)=tr_con_diag('Overwrite',T,T)
          select case(trim(pTracer%getName()))
          case ('NOx')
            itcon_3Dsrc(nOther,n)=tr_con_diag('Lightning',T,T)
          end select
#ifdef TRACERS_NITRATE
          select case (trim(pTracer%getName()))
            case ('HNO3')
              itcon_3Dsrc(nOther,n)=tr_con_diag('Nitrate Chemistry',T,T)
          end select
#endif

        case ('codirect')
          kt_power_change(n) = -13
          itcon_decay(n)=tr_con_diag('DECAY',T,T)

#ifdef TRACERS_AEROSOLS_SOA
        case ('isopp1g','isopp1a','isopp2g','isopp2a',
     &        'apinp1g','apinp1a','apinp2g','apinp2a')
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Chemistry',T,T)
#endif  /* TRACERS_AEROSOLS_SOA */

        case ('GLT')
          kt_power_change(n) = -17
          itcon_3Dsrc(1,n)=tr_con_diag('L1 overwriting',T,T)

        case ('HTO')
          itcon_decay(n)=tr_con_diag('DECAY',T,T)

        case ('DMS')
          itcon_surf(1,n)=tr_con_diag('Ocean src',T)
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Chemistry',T,T)

        case ('MSA')
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Chemistry',T,T)

        case ('SO2')
          itcon_3Dsrc(nVolcanic,n)=tr_con_diag('Volcanic src',T,T)
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Chem. src',T,T)
          itcon_3Dsrc(nChemLoss,n)=tr_con_diag('Chem. sink',T,T)

        case ('SO4')
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Gas phase src',T,T)
          itcon_3Dsrc(nVolcanic,n)=tr_con_diag('Volcanic src',T,T)

        case ('BCII', 'BCIA', 'BCB', 'OCII', 'OCIA', 'OCB',
     &        'vbsGm2', 'vbsGm1', 'vbsGz',  'vbsGp1', 'vbsGp2',
     &        'vbsGp3', 'vbsGp4', 'vbsGp5', 'vbsGp6',
     &        'vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &        'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
          select case(trim(pTracer%getName()))
          case ('vbsGm2', 'vbsGm1', 'vbsGz',  'vbsGp1', 'vbsGp2',
     &          'vbsGp3', 'vbsGp4', 'vbsGp5', 'vbsGp6')
            itcon_3Dsrc(nChemistry,n)=tr_con_diag('Aging source',T,T)
            itcon_3Dsrc(nChemLoss,n)=tr_con_diag('Aging loss',T,T)
            itcon_3Dsrc(nOther,n)=tr_con_diag('Part. loss',T,T)
          case ('vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &          'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
            itcon_3Dsrc(nChemistry,n)=tr_con_diag('Part. source',T,T)
          case ('BCII', 'OCII')
            itcon_3Dsrc(nChemistry,n)=tr_con_diag('Aging loss',T,T)
          case ('BCIA', 'OCIA')
            itcon_3Dsrc(nChemistry,n)=tr_con_diag('Aging source',T,T)
          end select

        case ('SO4_d1', 'SO4_d2','SO4_d3','N_d1','N_d2','N_d3')
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Gas phase change',T,T)

        case ('NH3','H2SO4')
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Gas phase change',T,T)
#ifdef TRACERS_TOMAS
          select case (trim(pTracer%getName()))
          case ('H2SO4')
            itcon_3Dsrc(nOther,n)=tr_con_diag('Microphysics change',T,T)
          end select
#endif

        case ('NH4', 'NO3p')
#ifndef TRACERS_TOMAS
          itcon_3Dsrc(1,n)=tr_con_diag('Gas phase change',T,T)
#else
          itcon_3Dsrc(1,n)=tr_con_diag('Microphysics change',T,T)
#endif

        case ('Be7', 'Be10')
          itcon_3Dsrc(1,n)=tr_con_diag('COSMO SRC',T,T)
          if (trim(pTracer%getName()).eq."Be7") then
            itcon_decay(n)=tr_con_diag('DECAY',T,T)
          end if

        case ('Pb210')
          itcon_3Dsrc(1,n)=tr_con_diag('RADIO SRC',T,T)
          itcon_decay(n)=tr_con_diag('DECAY',T,T)

        case ('H2O2_s')
          itcon_3Dsrc(1,n)=tr_con_diag('Gas phase src',T,T)
          itcon_3Dsrc(2,n)=tr_con_diag('Gas phase sink',T,T)

#ifndef TRACERS_WATER
        case ('seasalt1','seasalt2','OCocean'
     &         ,'Clay','Silt1','Silt2','Silt3','Silt4','Silt5'
     &         ,'ClayIlli' ,'ClayKaol','ClaySmec','ClayCalc','ClayQuar'
     &         ,'ClayFeld' ,'ClayHema','ClayGyps','ClayIlHe','ClayKaHe'
     &         ,'ClaySmHe' ,'ClayCaHe','ClayQuHe','ClayFeHe','ClayGyHe'
     &         ,'Sil1Quar' ,'Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps'
     &         ,'Sil1Illi' ,'Sil1Kaol','Sil1Smec','Sil1QuHe','Sil1FeHe'
     &         ,'Sil1CaHe' ,'Sil1GyHe','Sil1IlHe','Sil1KaHe','Sil1SmHe'
     &         ,'Sil2Quar' ,'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps'
     &         ,'Sil2Illi' ,'Sil2Kaol','Sil2Smec','Sil2QuHe','Sil2FeHe'
     &         ,'Sil2CaHe' ,'Sil2GyHe','Sil2IlHe','Sil2KaHe','Sil2SmHe'
     &         ,'Sil3Quar' ,'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps'
     &         ,'Sil3Illi' ,'Sil3Kaol','Sil3Smec','Sil3QuHe','Sil3FeHe'
     &         ,'Sil3CaHe' ,'Sil3GyHe','Sil3IlHe','Sil3KaHe','Sil3SmHe'
     &         ,'Sil4Quar' ,'Sil4Feld','Sil4Calc','Sil4Hema','Sil4Gyps'
     &         ,'Sil4Illi' ,'Sil4Kaol','Sil4Smec','Sil4QuHe','Sil4FeHe'
     &         ,'Sil4CaHe' ,'Sil4GyHe','Sil4IlHe','Sil4KaHe','Sil4SmHe'
     &         ,'Sil5Quar' ,'Sil5Feld','Sil5Calc','Sil5Hema','Sil5Gyps'
     &         ,'Sil5Illi' ,'Sil5Kaol','Sil5Smec','Sil5QuHe','Sil5FeHe'
     &         ,'Sil5CaHe' ,'Sil5GyHe','Sil5IlHe','Sil5KaHe','Sil5SmHe')
          itcon_wt(n)=tr_con_diag('WET DEP',T)
#endif  /* not TRACERS_WATER */

c- Species including AMP  emissions - 2D sources and 3D sources
        case('M_AKK_SU','M_ACC_SU','M_OCC_OC','M_BC1_BC',
     *       'M_SSA_SS','M_SSC_SS','M_SSS_SS','M_DD1_DU','M_DD2_DU',
     *       'M_BOC_BC','M_BOC_OC',
     *       'M_NO3   ','M_NH4   ','M_H2O   ','M_DD1_SU',
     *       'M_DS1_SU','M_DS1_DU','M_DD2_SU',
     *       'M_DS2_SU','M_DS2_DU','M_SSA_SU',
     *       'M_OCC_SU','M_BC1_SU',
     *       'M_BC2_SU','M_BC2_BC','M_BC3_SU',
     *       'M_BC3_BC','M_DBC_SU','M_DBC_BC','M_DBC_DU',
     *       'M_BOC_SU',
     *       'M_BCS_SU','M_BCS_BC','M_MXX_SU','M_MXX_BC',
     *       'M_MXX_OC','M_MXX_DU','M_MXX_SS','M_OCS_SU',
     *       'M_OCS_OC','M_SSS_SU')
          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Gas phase change',T,T)
          itcon_3Dsrc(nOther,n)=tr_con_diag('AMP source',T,T)
          select case (trim(pTracer%getName()))
            case ('M_SSA_SS','M_SSC_SS','M_SSS_SS','M_DD1_DU'
     *           ,'M_DD2_DU')
              itcon_surf(1,n)=tr_con_diag('Emission 2D AMP',T,T)
            case ('M_AKK_SU','M_ACC_SU',
     &            'M_BC1_BC','M_OCC_OC','M_BOC_BC','M_BOC_OC')
              select case (trim(pTracer%getName()))
              case ('M_AKK_SU','M_ACC_SU')
                itcon_3Dsrc(nVolcanic,n)=tr_con_diag('Volcanic src',T,T)
              end select
          end select
c Processes AMP Budget
          itcon_AMP(1,n)=tr_con_diag('P1 Nucleation',T,T)
          itcon_AMP(2,n)=tr_con_diag('P2 Coagulation',T,T)
          itcon_AMP(3,n)=tr_con_diag('P3 Condensation',T,T)
          itcon_AMP(4,n)=tr_con_diag('P4 Incloud',T,T)
          itcon_AMP(5,n)=tr_con_diag('P5 Intermode Loss',T,T)
          itcon_AMP(6,n)=tr_con_diag('P6 Mode Transf',T,T)
          itcon_AMP(7,n)=tr_con_diag('P7 AMP Budget',T,T)

        case ('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 '
     *       ,'N_DS2_1 ','N_SSA_1 ','N_SSC_1 ','N_OCC_1 ','N_BC1_1 '
     *       ,'N_BC2_1 ','N_BC3_1 ','N_DBC_1 ','N_BOC_1 ','N_BCS_1 '
     *       ,'N_MXX_1 ','N_OCS_1 ')

          kt_power_change(n) = 5
          kt_power_inst(n) = 3

          itcon_3Dsrc(nChemistry,n)=tr_con_diag('Gas phase change',T,T)
          itcon_AMPm(1,n)=tr_con_diag('Wet Diameter',T)
          itcon_AMPm(2,n)=tr_con_diag('Mode AktivPart',T)
          itcon_AMPm(3,n)=tr_con_diag('Dry Diameter',T)
c     Processes AMP Budget
          itcon_AMP(1,n)=tr_con_diag('P1 Nucleation',T,T)
          itcon_AMP(2,n)=tr_con_diag('P2 Coagulation',T,T)
          itcon_AMP(3,n)=tr_con_diag('P3 NOTHING',T,T)
          itcon_AMP(4,n)=tr_con_diag('P4 Intermode Coag',T,T)
          itcon_AMP(5,n)=tr_con_diag('P5 Intramode Tr',T,T)
          itcon_AMP(6,n)=tr_con_diag('P6 Mode Transf',T,T)
          itcon_AMP(7,n)=tr_con_diag('P7 AMP Budget',T,T)

#ifdef TRACERS_TOMAS

        case ('SOAgas')
!TOMAS - here needs lots of work~!
          itcon_3Dsrc(1,n)=tr_con_diag('Microphysics change',T,T)
!          itcon_surf(1,n)=tr_con_diag('Terpene_source',T)

       case('ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     *    'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     *    'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15',
     *    'ANACL_01','ANACL_02','ANACL_03','ANACL_04','ANACL_05',
     *    'ANACL_06','ANACL_07','ANACL_08','ANACL_09','ANACL_10',
     *    'ANACL_11','ANACL_12','ANACL_13','ANACL_14','ANACL_15',
     *    'AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05',
     *    'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10',
     *    'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15',
     *    'AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05',
     *    'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10',
     *    'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15',
     *    'AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05',
     *    'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10',
     *    'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15',
     *    'AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05',
     *    'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10',
     *    'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15',
     *    'ADUST_01','ADUST_02','ADUST_03','ADUST_04','ADUST_05',
     *    'ADUST_06','ADUST_07','ADUST_08','ADUST_09','ADUST_10',
     *    'ADUST_11','ADUST_12','ADUST_13','ADUST_14','ADUST_15',
     *    'ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     *    'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     *    'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15')

          itcon_3Dsrc(nOther,n)=tr_con_diag('Microphysics',T,T)
c     Processes TOMAS Budget
          itcon_TOMAS(1,n)=tr_con_diag('Condensation',T)
          itcon_TOMAS(2,n)=tr_con_diag('Coagulation',T)
          itcon_TOMAS(3,n)=tr_con_diag('Nucleation',T)
          itcon_TOMAS(4,n)=tr_con_diag('Aqoxid SO4 MCV',T)
          itcon_TOMAS(5,n)=tr_con_diag('Aqoxid SO4 LGS',T)
          itcon_TOMAS(6,n)=tr_con_diag('Mk_Nk Fix',T)
          itcon_TOMAS(7,n)=tr_con_diag('Aeroupdate',T)
          itcon_subcoag(n)=tr_con_diag('subgrid coag',T)

       select case (trim(pTracer%getName()))

         case ('ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     *        'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     *        'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15')

         itcon_3Dsrc(nVolcanic,n)=tr_con_diag('Volcanic src',T,T)

         case ('AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05',
     *        'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10',
     *        'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15',
     *        'AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05',
     *        'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10',
     *        'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15')

          itcon_3Dsrc(nChemistry,n)=tr_con_diag('ECOB Aging',T,T)

         case ('AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05',
     *        'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10',
     *        'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15',
     *        'AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05',
     *        'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10',
     *        'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15')

          itcon_3Dsrc(nChemistry,n)=tr_con_diag('OCOB Aging',T,T)

c     - Species including TOMAS  emissions - 2D sources and 3D sources
         case('ANACL_01','ANACL_02','ANACL_03','ANACL_04','ANACL_05',
     *        'ANACL_06','ANACL_07','ANACL_08','ANACL_09','ANACL_10',
     *        'ANACL_11','ANACL_12','ANACL_13','ANACL_14','ANACL_15')

          itcon_surf(1,n)=tr_con_diag('2D src',T)

         case('ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     *        'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     *        'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15')

          itcon_3Dsrc(1,n)=tr_con_diag('SO4 3D src',T,T)
          itcon_3Dsrc(2,n)=tr_con_diag('EC 3D src',T,T)
          itcon_3Dsrc(4,n)=tr_con_diag('OC 3D src',T,T) ! why 4 and not 3?
          itcon_surf(1,n)=tr_con_diag('2D src by SO4',T)
          itcon_surf(2,n)=tr_con_diag('2D src by EC',T)
          itcon_surf(3,n)=tr_con_diag('2D src by OC',T)
          itcon_surf(4,n)=tr_con_diag('2D src by SS',T)
          itcon_surf(5,n)=tr_con_diag('2D src by DU',T)

         case('ADUST_01','ADUST_02','ADUST_03','ADUST_04','ADUST_05',
     *        'ADUST_06','ADUST_07','ADUST_08','ADUST_09','ADUST_10',
     *        'ADUST_11','ADUST_12','ADUST_13','ADUST_14','ADUST_15')
          itcon_surf(1,n)=tr_con_diag('2D src',T)

       end select
#endif /* TRACERS_TOMAS */
        end select

        scale_inst(n)   = 10d0**(-kt_power_inst(n))
        scale_change(n) = 10d0**(-kt_power_change(n))
#ifdef TRACERS_TOMAS
        if(n.ge.n_ANUM(1).and.n.lt.n_AH2O(1))THEN
        inst_unit(n) = unit_string(kt_power_inst(n),  '# m-2)')
        sum_unit(n)  = unit_string(kt_power_change(n),'# m-2 s-1)')
        else
        inst_unit(n) = unit_string(kt_power_inst(n),  'kg m-2)')
        sum_unit(n)  = unit_string(kt_power_change(n),'kg m-2 s-1)')
        endif
#else
        inst_unit(n) = unit_string(kt_power_inst(n),  'kg m-2)')
        sum_unit(n)  = unit_string(kt_power_change(n),'kg m-2 s-1)')
#endif

        CALL SET_TCON(QCON,pTracer%getName(),QSUM,inst_unit(n),
     *       sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs,CONPT)
        qcon(npts_common+1:) = .false. ! reset to defaults for next tracer
        qsum(npts_common+1:) = .false. ! reset to defaults for next tracer
        qcon(10)  = .false.     ! reset to defaults for next tracer
        qsum(10)  = .false.     ! reset to defaults for next tracer
        conpts=''
        conpt=conpt0

        call iter%next()
      end do

      natmtrcons = tracers%size()

#ifdef TRACERS_OCEAN
      atmocn%natmtrcons = natmtrcons
#endif

#endif /* TRACERS_ON */

      return
      end subroutine init_tracer_cons_diag

      subroutine init_jls_diag
!@sum init_jls_diag Initialise zonal mean/height tracer diags
!@auth Gavin Schmidt
      use Tracer_mod, only: Tracer
      use TracerSurfaceSource_mod, only: TracerSurfaceSource
      use TRACER_COM, only: ntm
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE MODEL_COM, only: dtsrc
      use TRACER_COM, only: n_SO2, naircraft, nbiomass, nchemistry
      use TRACER_COM, only: nOther, nOverwrite, nVolcanic, nChemloss
      use TRACER_COM, only: ntsurfsrc, tracers, do_aircraft, aqchem_list
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: n_ANUM, n_AECOB, n_AOCOB
#endif
      USE DIAG_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      use trdust_mod, only: nDustEmjl, nDustEm2jl, nDustEv1jl,
     &   nDustEv2jl, nDustWthjl, imDust
#endif
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
      USE CLOUDS, ONLY : diag_wetdep
#endif
#endif /* TRACERS_ON */
      use OldTracer_mod, only: trname, ntm_power, src_dist_index,
     &                         nBBsources,do_fire
      implicit none
      integer k,n,kk,ltop,n_src
      character*50 :: unit_string
      class (Tracer), pointer :: pTracer
      type (TracerSurfaceSource), pointer :: sources(:)
      type (TracerSurfaceSource), pointer :: SO2sources(:)
      type (TracerSurfaceSource), pointer :: AECOB01sources(:)
      type (TracerSurfaceSource), pointer :: AOCOB01sources(:)

C**** Please note that short names for diags i.e. sname_jls are used
C**** in special ways and MUST NOT contain spaces, commas or % signs.
C**** Underscores and minus signs are allowed.

C**** Define a max layer for some optionally trop/strat tracers
      LTOP = LM

#ifdef TRACERS_ON
C**** Tracer sources and sinks
C**** Defaults for jls (sources, sinks, etc.)
C**** These need to be 'hand coded' depending on circumstances
      do k=1,ktajls             ! max number of sources and sinks
        jgrid_jls(k) = 1
        jwt_jls(k) = 1
        ia_jls(k) = ia_src
        scale_jls(k) = 1./DTsrc
      end do
      jls_grav=0
#ifdef TRACERS_WATER
C**** set defaults for some precip/wet-dep related diags
      jls_prec(:,:)=0
#endif

      k = 0
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      pTracer => tracers%getReference('SO2')
      SO2sources => pTracer%surfaceSources
#endif
#ifdef TRACERS_TOMAS
      pTracer => tracers%getReference('AECOB_01')
      AECOB01sources => pTracer%surfaceSources
      pTracer => tracers%getReference('AOCOB_01')
      AOCOB01sources => pTracer%surfaceSources
#endif
      do n=1,NTM
        if (src_dist_index(n)/=0) cycle

!=============================================!
! emissions for all tracers, if they have any !
!=============================================!

! handle exceptions first (e.g. SO4 emissions are listed under SO2)
      select case (trname(n))
      case ('SO4',
     &      'M_AKK_SU','M_ACC_SU',
     &      'ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     &      'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     &      'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15')
        n_src = n_SO2
#ifdef TRACERS_TOMAS
      case ('AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05',
     &      'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10',
     &      'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15',
     &      'AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05',
     &      'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10',
     &      'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15')
        n_src = n_AECOB(1)
      case ('AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05',
     &      'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10',
     &      'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15',
     &      'AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05',
     &      'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10',
     &      'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15')
        n_src = n_AOCOB(1)
#endif  /* TRACERS_TOMAS */
      case default
        n_src = n
      end select
      pTracer => tracers%getReference(trname(n_src))
      sources => pTracer%surfaceSources

! aqueous chemistry sources and sinks
      if (allocated(aqchem_list)) then
      if (any(n.eq.aqchem_list)) then
        k = k + 1
        jls_incloud(1,n) = k
        sname_jls(k) = trim(trname(n))//'_mc_cloud_aqchem'
        lname_jls(k) = trim(trname(n))//' mc cloud aqchem'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

        k = k + 1
        jls_incloud(2,n) = k
        sname_jls(k) = trim(trname(n))//'_ss_cloud_aqchem'
        lname_jls(k) = trim(trname(n))//' ss cloud aqchem'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
      endif
      endif

! surface emissions
      do kk=1,ntsurfsrc(n_src)
        k = k + 1
        jls_source(kk,n) = k
        sname_jls(k) = trim(trname(n))//'_'//
     &                 trim(sources(kk)%sourceName)
        lname_jls(k) = trim(trname(n))//' '//
     &                 trim(sources(kk)%sourceLname)
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+11
        select case(trname(n))
        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     *    'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     *    'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15')
          units_jls(k) = unit_string(jls_power(k),'# s-1')
        case default
          units_jls(k) = unit_string(jls_power(k),'kg s-1')
        end select
      end do

! aircraft emissions
      if(do_aircraft(n_src)) then
        k = k + 1
        jls_3Dsource(nAircraft,n) = k
        sname_jls(k) = trim(trname(n))//'_aircraft_src'
        lname_jls(k) = trim(trname(n))//' aircraft source'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
      end if

! biomass burning emissions
      if (nBBsources(n_src) .gt. 0) then
        k = k + 1
        jls_3Dsource(nBiomass,n) = k
        sname_jls(k) = trim(trname(n))//'_biomass_src'
        lname_jls(k) = trim(trname(n))//' biomass source'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
      endif

!=============================!
! Tracer-specific diagnostics !
!=============================!
      select case (trname(n))

c      case ('SF6','SF6_c','nh5','nh50','e90','st8025','tape_rec','aoa','aoanh','nh15')
c        call layer1_init_jls(k,n,trname(n))
c      case ('CFCn')
c        call layer1_init_jls(k,n,trname(n))
      case ('CO2n')
        call CO2n_init_jls(k,n,'CO2n')
      case ('Rn222')
        call Rn222_init_jls(k,n,'Rn222')
! keep AIJ and AJL CO2 sources in same order !!
c      case ('CO2')
c        call CO2_init_jls(k,n,'CO2')
      case ('N2O')
        call N2O_init_jls(k,n,'N2O')
      case ('CFC11')   !!! should start April 1
c        k = k + 1
c        jls_source(1,n) = k
c        sname_jls(k) = 'L1_sink_'//trim(trname(n))
c        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//' BY SOURCE, L1'
c        jls_ltop(k) = 1
c        jls_power(k) = -1
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Stratos_chem_change_'//trim(trname(n))
        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//
     &                 ' BY CHEMISTRY IN STRATOS'
        jls_ltop(k) = lm
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

c      case ('14CO2')   !!! should start 10/16
c        k = k + 1
c        jls_source(1,n) = k
c        sname_jls(k) = 'L1_sink_'//trim(trname(n))
c        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//' by SINK, L1'
c        jls_ltop(k) = 1
c        jls_power(k) = -4
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of_'//trim(trname(n))
        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//' BY CHEMISTRY'
        jls_ltop(k) = LTOP
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_3Dsource(nOverwrite,n) = k
        sname_jls(k) = 'overwrite_source_of_'//trim(trname(n))
        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//' BY OVERWRITE'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
#else
c        k = k + 1
c        jls_source(6,n) = k
c        sname_jls(k) = 'Soil_sink_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' sink due to soil absorption'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(7,n) = k
c        sname_jls(k) = 'Termite_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Termite source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(9,n) = k
c        sname_jls(k) = 'Ocean_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Ocean source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(10,n) = k
c        sname_jls(k) = 'Fresh_Water_lake_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Fresh Water lake source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(11,n) = k
c        sname_jls(k) = 'Misc_Ground_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Misc_Ground source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(14,n) = k
c        sname_jls(k) = 'Wetlands+Tundra_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Wetlands+Tundra source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(1,n) = k
c        sname_jls(k) = 'Animal_source_of_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Animal source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(2,n) = k
c        sname_jls(k) = 'Coal_Mine_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Coal Mine source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(3,n) = k
c        sname_jls(k) = 'Gas_Leak_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Gas Leak source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(4,n) = k
c        sname_jls(k) = 'Gas_Venting_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Gas Venting source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(5,n) = k
c        sname_jls(k) = 'Municipal_solid_waste_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Municipal solid waste source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(8,n) = k
c        sname_jls(k) = 'Coal_combustion_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Coal combustion source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(12,n) = k
c        sname_jls(k) = 'Biomass_burning_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Biomass burning source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(13,n) = k
c        sname_jls(k) = 'Rice_Cultivation_source_'//trim(trname(n))
c        lname_jls(k) = trim(trname(n))//' Rice Cultivation source'
c        jls_ltop(k) = 1
c        jls_power(k) = 0
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Tropos_Chem_change_'//trim(trname(n))
        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//
     &                 ' BY CHEMISTRY IN TROPOSPHERE'
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Stratos_Chem_change_'//trim(trname(n))
        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//
     &                 ' BY CHEMISTRY IN STRATOS'
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
#endif

      case ('O3')
c       k = k + 1
c        jls_source(1,n) = k
c        sname_jls(k) = 'Deposition_L1_'//trim(trname(n))
c        lname_jls(k) = 'Change of O3 by Deposition in Layer 1'
c        jls_ltop(k) = 1
c        jls_power(k) = 1
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
       k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Strat_Chem_change_'//trim(trname(n))
        lname_jls(k) = 'Change of O3 by Chemistry in Stratos'
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
       k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Trop_Chem_Prod_change_'//trim(trname(n))
        lname_jls(k) = 'Change of O3 by Chem Prod. in Troposphere'
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
       k = k + 1
        jls_3Dsource(3,n) = k
        sname_jls(k) = 'Trop_Chem_Loss_change_'//trim(trname(n))
        lname_jls(k) = 'Change of O3 by Chem Loss in Troposphere'
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

#ifdef TRACERS_PASSIVE

      case ('nh5')
        k = k + 1
        jls_decay(n) = k   ! decay loss
        sname_jls(k) = 'Decay_of_'//trim(trname(n))
        lname_jls(k) = 'LOSS OF '//trim(trname(n))//' BY DECAY'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      case ('nh50')
        k = k + 1
        jls_decay(n) = k   ! decay loss
        sname_jls(k) = 'Decay_of_'//trim(trname(n))
        lname_jls(k) = 'LOSS OF '//trim(trname(n))//' BY DECAY'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      case ('e90')
        k = k + 1
        jls_decay(n) = k   ! decay loss
        sname_jls(k) = 'Decay_of_'//trim(trname(n))
        lname_jls(k) = 'LOSS OF '//trim(trname(n))//' BY DECAY'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      case ('nh15')
        k = k + 1
        jls_decay(n) = k   ! decay loss
        sname_jls(k) = 'Decay_of_'//trim(trname(n))
        lname_jls(k) = 'LOSS OF '//trim(trname(n))//' BY DECAY'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

#endif

#ifdef TRACERS_WATER
C**** generic ones for many water tracers
      case ('Water', 'H2O18', 'HDO', 'HTO', 'H2O17' )
       k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Evap_'//trim(trname(n))
        lname_jls(k) = 'EVAPORATION OF '//trim(trname(n))
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SECONDS_PER_DAY/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm day-1')
       k = k + 1
        jls_isrc(2,n) = k
        sname_jls(k) = 'Ocn_Evap_'//trim(trname(n))
        lname_jls(k) = 'OCEAN EVAP OF '//trim(trname(n))
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SECONDS_PER_DAY/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm day-1')
       k = k + 1
        jls_prec(1,n)=k
        sname_jls(k) = 'Precip_'//trim(trname(n))
        lname_jls(k) = 'PRECIPITATION OF '//trim(trname(n))
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SECONDS_PER_DAY/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm day-1')
       k = k + 1
        jls_prec(2,n)=k
        sname_jls(k) = 'Ocn_Precip_'//trim(trname(n))
        lname_jls(k) = 'OCEAN PRECIP OF '//trim(trname(n))
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SECONDS_PER_DAY/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm day-1')

C**** special one unique to HTO
      if (trname(n).eq."HTO") then
       k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trim(trname(n))
        lname_jls(k) = 'LOSS OF '//TRIM(trname(n))//' BY DECAY'
        jls_ltop(k) = lm
        jls_power(k) = ntm_power(n)+8
        scale_jls(k) = 1./DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
      end if
#endif

!#ifdef TRACERS_NITRATE
!       case ('HNO3')
!        k = k + 1
!        jls_3Dsource(nChemistry,n) = k
!        sname_jls(k) = 'chemistry_nitrate_of_'//trim(trname(n))
!        lname_jls(k) = 'CHANGE OF HNO3 BY NITRATE CHEM'
!        jls_ltop(k) = LTOP
!        jls_power(k) = 0
!        units_jls(k) = unit_string(jls_power(k),'kg s-1')
!#endif

#ifdef SHINDELL_STRAT_EXTRA
      case ('GLT')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'L1_overwrite_soure_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' L1 overwrite source'
        jls_ltop(k) = 1
        jls_power(k) = -5
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
#endif

      case ('codirect')
        k = k + 1
        jls_decay(n) = k   ! decay loss
        sname_jls(k) = 'Decay_of_'//trim(trname(n))
        lname_jls(k) = 'LOSS OF '//trim(trname(n))//' BY DECAY'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      case ('HCl','HOCl','ClONO2','HBr','HOBr','BrONO2','CFC',
     &      'BrOx','ClOx','Alkenes','Paraffin','Isoprene','CO',
#ifdef TRACERS_dCO
     *      'd13Calke','d13CPAR',
     *      'd17OPAN', 'd18OPAN', 'd13CPAN',
     *      'dMe17OOH', 'dMe18OOH', 'd13MeOOH',
     *      'dHCH17O', 'dHCH18O', 'dH13CHO',
     *      'dC17O', 'dC18O', 'd13CO',
#endif  /* TRACERS_dCO */
     &      'N2O5','HNO3','H2O2','CH3OOH','HCHO','HO2NO2','PAN',
     &      'AlkylNit','Ox','NOx','stratOx','Terpenes')
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of_'//trim(trname(n))
        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//' BY CHEMISTRY'
        jls_ltop(k) = LM
        select case(trname(n))
        case ('Ox','stratOx')
          jls_power(k) = 1
        case default
          jls_power(k) = -1
        end select
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        select case(trname(n))
        case ('Alkenes','Paraffin','Isoprene','CO','N2O5','HNO3',
#ifdef TRACERS_dCO
     *      'd13Calke','d13CPAR',
     *      'd17OPAN', 'd18OPAN', 'd13CPAN',
     *      'dMe17OOH', 'dMe18OOH', 'd13MeOOH',
     *      'dHCH17O', 'dHCH18O', 'dH13CHO',
     *      'dC17O', 'dC18O', 'd13CO',
#endif  /* TRACERS_dCO */
     &  'H2O2','CH3OOH','HCHO','HO2NO2','PAN','AlkylNit','Ox',
     &  'Terpenes','NOx','stratOx','BrOx','ClOx')
          k = k + 1
          jls_3Dsource(nOverwrite,n) = k
          sname_jls(k) = 'overwrite_source_of_'//trim(trname(n))
          lname_jls(k) =
     &    'CHANGE OF '//trim(trname(n))//' BY OVERWRITE'
          jls_ltop(k) = LM
          jls_power(k) = -1
          units_jls(k) = unit_string(jls_power(k),'kg s-1')
        case ('CFC')  ! L=1 overwrite only.
          k = k + 1
          jls_3Dsource(nOverwrite,n) = k
          sname_jls(k) = 'overwrite_source_of_'//trim(trname(n))
          lname_jls(k) =
     &    'CHANGE OF '//trname(n)//' BY OVERWRITE'
          jls_ltop(k) = 1 ! L=1 overwrite only
          jls_power(k) = -1
          units_jls(k) = unit_string(jls_power(k),'kg s-1')
        end select
        select case(trname(n))
        case('NOx')
          k = k + 1
          jls_3Dsource(nOther,n) = k
          sname_jls(k) = 'lightning_source_of_'//trim(trname(n))
          lname_jls(k) = 'CHANGE OF '//trim(trname(n))//' BY LIGHTNING'
          jls_ltop(k) = LM
          jls_power(k) = -2
          units_jls(k) = unit_string(jls_power(k),'kg s-1')
        end select

#ifdef TRACERS_AEROSOLS_SOA
      case ('isopp1g','isopp2g','apinp1g','apinp2g')
c put in chemical production
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of_'//trim(trname(n))
        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//' BY CHEMISTRY'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      case ('isopp1a','isopp2a','apinp1a','apinp2a')
c put in chemical production
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of_'//trim(trname(n))
        lname_jls(k) = 'CHANGE OF '//trim(trname(n))//' BY CHEMISTRY'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c gravitational settling of SOA
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
        lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
#endif  /* TRACERS_AEROSOLS_SOA*/

      case ('DMS')
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Ocean_source_of_'//trim(trname(n))
        lname_jls(k) = 'DMS ocean source'
        jls_ltop(k) = 1
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
C
        k = k + 1
        jls_isrc(2,n) = k
        sname_jls(k) = 'TKE_Contribution_'//trim(trname(n))
        lname_jls(k) = 'SGSWSP TKE'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_isrc(3,n) = k
        sname_jls(k) = 'Wet_Conv_Contr_'//trim(trname(n))
        lname_jls(k) = 'SGSWSP Wet Conv'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_isrc(4,n) = k
        sname_jls(k) = 'Dry_Conv_Contr_'//trim(trname(n))
        lname_jls(k) = 'SGSWSP Dry Conv'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_isrc(5,n) = k
        sname_jls(k) = 'SGSWSP-old_'//trim(trname(n))
        lname_jls(k) = 'DMS SGSWP-old/old'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Chemical_sink_of_'//trim(trname(n))
        lname_jls(k) = 'DMS chemical loss'
        jls_ltop(k) =LM
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

       case ('MSA')
c put in chemical production of MSA
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'chemistry_source_of_'//trim(trname(n))
        lname_jls(k) = 'Chemical production of MSA'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c gravitational settling of MSA
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
        lname_jls(k) = 'Gravitational Settling of MSA'
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

       case ('SO2')
c volcanic production of SO2
        k = k + 1
        jls_3Dsource(nVolcanic,n) = k
        sname_jls(k) = trim(trname(n))//'_volcanic_src'
        lname_jls(k) = trim(trname(n))//' volcanic source'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c put in chemical production of SO2
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'dms_source_of_'//trim(trname(n))
        lname_jls(k) = 'production of SO2 from DMS'
        jls_ltop(k) = LM
        jls_power(k) =  1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c put in chemical sink of SO2
        k = k + 1
        jls_3Dsource(nChemloss,n) = k
        sname_jls(k) = 'chem_sink_of_'//trim(trname(n))
        lname_jls(k) = 'chemical sink of SO2'
        jls_ltop(k) = LM
        jls_power(k) =  1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        case ('SO4')
c gas phase source of SO4
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'gas_phase_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c volcanic source of SO4
        k = k + 1
        jls_3Dsource(nVolcanic,n) = k
        sname_jls(k) = trim(trname(n))//'_volcanic_src'
        lname_jls(k) = trim(trname(n))//' volcanic source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c gravitational settling of SO4
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
        lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

        case ('SO4_d1', 'SO4_d2', 'SO4_d3')
c gas phase source
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c gravitational settling
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
        lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

        case ('Be7')
c cosmogenic source from file
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Cosmogenic_src_of_'//trim(trname(n))
        lname_jls(k) = 'Be7 cosmogenic src'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg mb-1 m-2 s-1')
        jwt_jls(k) = 3
c radioactive decay
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trim(trname(n))
        lname_jls(k) = 'Loss of Be7 by decay'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg mb-1 m-2 s-1')
        jwt_jls(k) = 3
c gravitational settling
        k = k + 1
        jls_grav(n) = k   ! special array grav. settling sinks
        sname_jls(k) = 'Grav_Settle_of_'//trim(trname(n))
        lname_jls(k) = 'Loss of Be7 by grav settling'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg mb-1 m-2 s-1')
        jwt_jls(k) = 3

        case ('Be10')
c cosmogenic source from file/same as Be7
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Cosmogenic_src_of_'//trim(trname(n))
        lname_jls(k) = 'Be10 cosmogenic src'
        jls_ltop(k) = lm
        jls_power(k) = -28  !may need changing around
        units_jls(k) = unit_string(jls_power(k),'kg mb-1 m-2 s-1')
        jwt_jls(k) = 3
c gravitational settling
        k = k + 1
        jls_grav(n) = k   ! special array grav. settling sinks
        sname_jls(k) = 'Grav_Settle_of_'//trim(trname(n))
        lname_jls(k) = 'Loss of Be10 by grav settling'
        jls_ltop(k) = lm
        jls_power(k) = -28  !may need changing around
        units_jls(k) = unit_string(jls_power(k),'kg mb-1 m-2 s-1')
        jwt_jls(k) = 3

        case ('Pb210')
c source of Pb210 from Rn222 decay
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Radioactive_src_of_'//trim(trname(n))
        lname_jls(k) = 'Pb210 radioactive src'
        jls_ltop(k) = lm
        jls_power(k) =-26   ! -10  !may need to be changed
        units_jls(k) = unit_string(jls_power(k),'kg mb-1 m-2 s-1')
        jwt_jls(k) = 3
c radioactive decay
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trim(trname(n))
        lname_jls(k) = 'Loss of Pb210 by decay'
        jls_ltop(k) = lm
        jls_power(k) =-26   ! -10  !may need to be changed
        units_jls(k) = unit_string(jls_power(k),'kg mb-1 m-2 s-1')
        jwt_jls(k) = 3
c gravitational settling
        k = k + 1
        jls_grav(n) = k   ! special array grav. settling sinks
        sname_jls(k) = 'Grav_Settle_of_'//trim(trname(n))
        lname_jls(k) = 'Loss of Pb210 by grav settling'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg mb-1 m-2 s-1')
        jwt_jls(k) = 3

        case ('H2O2_s')
c gas phase source and sink of H2O2
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 2
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'gas_phase_sink_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' gas phase sink'
        jls_ltop(k) = LM
        jls_power(k) = 2
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c photolysis rate
        k = k + 1
        jls_phot = k
        sname_jls(k) = 'photolysis_rate_of_'//trim(trname(n))
        lname_jls(k) = 'photolysis rate of '//trim(trname(n))
        jls_ltop(k) =LM
        jls_power(k) =-9
        units_jls(k) = unit_string(jls_power(k),'/s')
      case ('vbsGm2', 'vbsGm1', 'vbsGz',  'vbsGp1', 'vbsGp2',
     &      'vbsGp3', 'vbsGp4', 'vbsGp5', 'vbsGp6')
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = trim(trname(n))//'_aging_source'
        lname_jls(k) = trim(trname(n))//' aging source'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_3Dsource(nChemloss,n) = k
        sname_jls(k) = trim(trname(n))//'_aging_loss'
        lname_jls(k) = trim(trname(n))//' aging loss'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_3Dsource(nOther,n) = k
        sname_jls(k) = trim(trname(n))//'_partitioning'
        lname_jls(k) = trim(trname(n))//' partitioning'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      case ('BCII', 'BCIA', 'BCB', 'OCII', 'OCIA', 'OCB',
     &      'vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &      'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
        select case(trname(n))
        case ('BCII', 'BCB', 'OCII', 'OCB',
     &        'vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &        'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
          select case(trname(n))
          case ('vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &          'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
            k = k + 1
            jls_3Dsource(nChemloss,n) = k
            sname_jls(k) = trim(trname(n))//'_partitioning'
            lname_jls(k) = trim(trname(n))//' partitioning'
            jls_ltop(k) = LM
            jls_power(k) = -1
            units_jls(k) = unit_string(jls_power(k),'kg s-1')
          end select
        case ('BCIA', 'OCIA')
          k = k + 1
          jls_3Dsource(nChemistry,n) = k
          sname_jls(k) = 'Aging_source_of_'//trim(trname(n))
          lname_jls(k) = trim(trname(n))//' aging source'
          jls_ltop(k) = LM
          jls_power(k) = -1
          units_jls(k) = unit_string(jls_power(k),'kg s-1')
        end select
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
        lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

#ifdef TRACERS_TOMAS
       case('ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     *    'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     *    'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15',
     *    'ANACL_01','ANACL_02','ANACL_03','ANACL_04','ANACL_05',
     *    'ANACL_06','ANACL_07','ANACL_08','ANACL_09','ANACL_10',
     *    'ANACL_11','ANACL_12','ANACL_13','ANACL_14','ANACL_15',
     *    'AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05',
     *    'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10',
     *    'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15',
     *    'AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05',
     *    'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10',
     *    'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15',
     *    'AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05',
     *    'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10',
     *    'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15',
     *    'AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05',
     *    'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10',
     *    'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15',
     *    'ADUST_01','ADUST_02','ADUST_03','ADUST_04','ADUST_05',
     *    'ADUST_06','ADUST_07','ADUST_08','ADUST_09','ADUST_10',
     *    'ADUST_11','ADUST_12','ADUST_13','ADUST_14','ADUST_15')
        k = k + 1
        jls_3Dsource(nOther,n) = k
        sname_jls(k) = 'Microphysics_src_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//'Microphysics src'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
        lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

        select case (trname(n))

        case ('ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     *       'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     *       'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15')

c volcanic source of SO4
        k = k + 1
        jls_3Dsource(nVolcanic,n) = k
        sname_jls(k) = trim(trname(n))//'_volcanic_src'
        lname_jls(k) = trim(trname(n))//' volcanic source'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c industrial source
        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     *    'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     *    'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15')
c SO4
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'SO4_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' SO4 source'
        jls_ltop(k) = LM
        jls_power(k) = 10
        units_jls(k) = unit_string(jls_power(k),'# s-1')
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'EC_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//'EC source'
        jls_ltop(k) = LM
        jls_power(k) = 10
        units_jls(k) = unit_string(jls_power(k),'# s-1')
        k = k + 1
        jls_3Dsource(4,n) = k
        sname_jls(k) = 'OC_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//'OC source'
        jls_ltop(k) = LM
        jls_power(k) = 10
        units_jls(k) = unit_string(jls_power(k),'# s-1')
        k = k + 1
        jls_3Dsource(nOther,n) = k
        sname_jls(k) = 'Microphysics_src_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//'Microphysics src'
        jls_ltop(k) = LM
        jls_power(k) = 10
        units_jls(k) = unit_string(jls_power(k),'# s-1')
c industrial source
        do kk=1,ntsurfsrc(n_ANUM(1))
          k = k + 1
          jls_source(kk,n) = k
          sname_jls(k) = trim(trname(n))//'_'//
     &                   trim(sources(kk)%sourceName)
          lname_jls(k) = trim(trname(n))//' '//
     &                   trim(sources(kk)%sourceLname)
          jls_ltop(k) = 1
          jls_power(k) =10
          units_jls(k) = unit_string(jls_power(k),'# s-1')
        enddo
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'NACL_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//'ANACL source'
        jls_ltop(k) = 1
        jls_power(k) =10
        units_jls(k) = unit_string(jls_power(k),'# s-1')
        k = k + 1
        jls_isrc(2,n) = k
        sname_jls(k) = 'Dust_source_of_'//trim(trname(n))
        lname_jls(k) =  trim(trname(n))//'ADUST source'
        jls_ltop(k) = 1
        jls_power(k) =1
        units_jls(k) = unit_string(jls_power(k),'# s-1')
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
        lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = 10
        units_jls(k) = unit_string(jls_power(k),'# s-1')

      case ('ANACL_01','ANACL_02','ANACL_03','ANACL_04','ANACL_05',
     *    'ANACL_06','ANACL_07','ANACL_08','ANACL_09','ANACL_10',
     *    'ANACL_11','ANACL_12','ANACL_13','ANACL_14','ANACL_15')
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Ocean_source_of_'//trim(trname(n))
        lname_jls(k) = 'Ocean source of '//trim(trname(n))
        jls_ltop(k) = 1
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      case ('AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05',
     *    'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10',
     *    'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15',
     *    'AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05',
     *    'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10',
     *    'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Aging_loss_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' aging loss'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      case ('AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05',
     *    'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10',
     *    'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15',
     *    'AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05',
     *    'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10',
     *    'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Aging_loss_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' aging loss'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

! TOMAS  : should I exclude aerosol water??

        case('ADUST_01','ADUST_02','ADUST_03','ADUST_04','ADUST_05',
     *    'ADUST_06','ADUST_07','ADUST_08','ADUST_09','ADUST_10',
     *    'ADUST_11','ADUST_12','ADUST_13','ADUST_14','ADUST_15')

        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Dust_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' dust source'
        jls_ltop(k) = 1
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

        end select

#endif /* TRACERS_TOMAS*/

      case ('seasalt1', 'seasalt2', 'OCocean')
c ocean source
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = trim(trname(n))//'_ocean_src'
        lname_jls(k) = trim(trname(n))//' ocean source'
        jls_ltop(k) = 1
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c gravitational settling
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = trim(trname(n))//'_grav_sett'
        lname_jls(k) = trim(trname(n))//' gravitational settling'
        jls_ltop(k) = LM
        select case (trname(n))
        case ('seasalt1', 'OCocean')
          jls_power(k) = -2
        case ('seasalt2')
          jls_power(k) =0
        end select
        units_jls(k) = unit_string(jls_power(k),'kg s-1')

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
        CASE('Clay','Silt1','Silt2','Silt3','Silt4','Silt5','ClayIlli'
     &         ,'ClayKaol','ClaySmec','ClayCalc','ClayQuar','ClayFeld'
     &         ,'ClayHema','ClayGyps','ClayIlHe','ClayKaHe','ClaySmHe'
     &         ,'ClayCaHe','ClayQuHe','ClayFeHe','ClayGyHe','Sil1Quar'
     &         ,'Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps','Sil1Illi'
     &         ,'Sil1Kaol','Sil1Smec','Sil1QuHe','Sil1FeHe','Sil1CaHe'
     &         ,'Sil1GyHe','Sil1IlHe','Sil1KaHe','Sil1SmHe','Sil2Quar'
     &         ,'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps','Sil2Illi'
     &         ,'Sil2Kaol','Sil2Smec','Sil2QuHe','Sil2FeHe','Sil2CaHe'
     &         ,'Sil2GyHe','Sil2IlHe','Sil2KaHe','Sil2SmHe','Sil3Quar'
     &         ,'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps','Sil3Illi'
     &         ,'Sil3Kaol','Sil3Smec','Sil3QuHe','Sil3FeHe','Sil3CaHe'
     &         ,'Sil3GyHe','Sil3IlHe','Sil3KaHe','Sil3SmHe','Sil4Quar'
     &         ,'Sil4Feld','Sil4Calc','Sil4Hema','Sil4Gyps','Sil4Illi'
     &         ,'Sil4Kaol','Sil4Smec','Sil4QuHe','Sil4FeHe','Sil4CaHe'
     &         ,'Sil4GyHe','Sil4IlHe','Sil4KaHe','Sil4SmHe','Sil5Quar'
     &         ,'Sil5Feld','Sil5Calc','Sil5Hema','Sil5Gyps','Sil5Illi'
     &         ,'Sil5Kaol','Sil5Smec','Sil5QuHe','Sil5FeHe','Sil5CaHe'
     &         ,'Sil5GyHe','Sil5IlHe','Sil5KaHe','Sil5SmHe')

        k=k+1
          jls_isrc(nDustEmjl,n)=k
          lname_jls(k)='Emission of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_emission'
          jls_ltop(k)=1
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg s-1')
        IF ( imDust == 0 .or. imDust >= 3 ) THEN
          k=k+1
          jls_isrc(nDustEm2jl,n)=k
          lname_jls(k)='Cubic emission of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_emission2'
          jls_ltop(k)=1
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg s-1')
        END IF
#ifndef TRACERS_DRYDEP
        k=k+1
          jls_isrc(nDustTurbjl,n)=k
          lname_jls(k)='Turbulent deposition of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_turb_depo'
          jls_ltop(k)=1
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg s-1')
#endif
        k=k+1
          jls_grav(n)=k
          lname_jls(k)='Gain by gravitational settling of '
     &         //TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_grav_sett'
          jls_ltop(k)=Lm
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg s-1')
#ifndef TRACERS_WATER
        k=k+1
          jls_wet(n)=k
          lname_jls(k)='Loss by wet deposition of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_wet_depo'
          jls_ltop(k)=Lm
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg s-1')
#endif
#endif /* TRACERS_DUST || TRACERS_MINERALS */

C**** Here are some more examples of generalised diag. configuration
c      n = n_dust
c        k = k + 1
c        jls_grav(n) = k   ! special array grav. settling sinks
c        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
c        lname_jls(k) = 'LOSS OF DUST BY SETTLING'
c        jls_ltop(k) = lm
c        jls_power(k) = -11
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')

      end select

#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
c**** additional wet deposition diagnostics
      IF (diag_wetdep == 1) THEN
        k=k+1
        jls_trdpmc(1,n)=k
        lname_jls(k)='MC Condensation of '//TRIM(trname(n))
        sname_jls(k)=TRIM(trname(n))//'_cond_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpmc(2,n)=k
        lname_jls(k)='Evaporated '//TRIM(trname(n))//' in MC Downdrafts'
        sname_jls(k)=TRIM(trname(n))//'_downeva_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpmc(3,n)=k
        lname_jls(k)='Condensed '//TRIM(trname(n))//' in MC CLW'
        sname_jls(k)=TRIM(trname(n))//'_conclw_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpmc(4,n)=k
        lname_jls(k)='Precipitated '//TRIM(trname(n))//' by MC'
        sname_jls(k)=TRIM(trname(n))//'_precip_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpmc(5,n)=k
        lname_jls(k)='Reevaporated '//TRIM(trname(n))//' from MC Precip'
        sname_jls(k)=TRIM(trname(n))//'_reevap_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpmc(6,n)=k
        lname_jls(k)='MC Washout of '//TRIM(trname(n))
        sname_jls(k)=TRIM(trname(n))//'_washout_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpls(1,n)=k
        lname_jls(k)='LS Washout of '//TRIM(trname(n))
        sname_jls(k)=TRIM(trname(n))//'_washout_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpls(2,n)=k
        lname_jls(k)='Precipitated '//TRIM(trname(n))//' by LS'
        sname_jls(k)=TRIM(trname(n))//'_precip_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpls(3,n)=k
        lname_jls(k)='Condensed '//TRIM(trname(n))// ' in LS CLW'
        sname_jls(k)=TRIM(trname(n))//'_conclw_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpls(4,n)=k
        lname_jls(k)='Reevaporated '//TRIM(trname(n))//' from LS Precip'
        sname_jls(k)=TRIM(trname(n))//'_reevap_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpls(5,n)=k
        lname_jls(k)='Evaporated '//TRIM(trname(n))//' from LS CLW'
        sname_jls(k)=TRIM(trname(n))//'_clwevap_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
        k=k+1
        jls_trdpls(6,n)=k
        lname_jls(k)='LS Condensation of '//TRIM(trname(n))
        sname_jls(k)=TRIM(trname(n))//'_cond_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg s-1')
      END IF
#endif

c
C**** Checks
      if (ntsurfsrc(n).gt.ntsurfsrcmax) then
!       write(6,*) ' ntsurfsrc too large for ',trname(n)
        if (am_i_root())
     &      write(6,*) ' Increase ntsurfsrcmax to at least',ntsurfsrc(n)
        call stop_model(
     &       ' Ntsurfsrc too large.  Increase ntsurfsrcmax',255)
      end if

      end do

C**** Additional Special JL diagnostics
C**** (not necessary associated with a particular tracer)
#ifdef TRACERS_SPECIAL_Shindell
        k = k + 1
        jls_ClOcon=k
        sname_jls(k) = 'ClO_conc'
        lname_jls(k) = 'ClO concentration'
        jls_ltop(k)  = LTOP
        jls_power(k) = -11
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'V/V air')
        k = k + 1
        jls_H2Ocon=k
        sname_jls(k) = 'H2O_conc'
        lname_jls(k) = 'H2O concentration'
        jls_ltop(k)  = LTOP
        jls_power(k) = -7
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'V/V air')
        k = k + 1
        jls_H2Ochem=k
        sname_jls(k) = 'H2O_chem'
        lname_jls(k) = 'H2O change due to chemistry'
        jls_ltop(k)  = LTOP
        jls_power(k) = -4
        scale_jls(k) = 1./DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_Oxp=k
        sname_jls(k) = 'Ox_chem_prod'
        lname_jls(k) = 'Ox production due to chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 2
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_Oxd=k
        sname_jls(k) = 'Ox_chem_dest'
        lname_jls(k) = 'Ox destruction due to chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 2
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_OxpT=k
        sname_jls(k) = 'trop_Ox_chem_prod'
        lname_jls(k) = 'Troposphere Ox prod by chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 2
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_OxdT=k
        sname_jls(k) = 'trop_Ox_chem_dest'
        lname_jls(k) = 'Troposphere Ox dest by chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 2
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_COp=k
        sname_jls(k) = 'CO_chem_prod'
        lname_jls(k) = 'CO production due to chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 1
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_COd=k
        sname_jls(k) = 'CO_chem_dest'
        lname_jls(k) = 'CO destruction due to chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 1
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
        k = k + 1
        jls_OHcon=k
        sname_jls(k) = 'OH_conc'
        lname_jls(k) = 'OH concentration'
        jls_ltop(k)  = LTOP
        jls_power(k) = 5
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'molecules cm-3')
c
        k = k + 1
        jls_H2Omr=k
        sname_jls(k) = 'H2O_mr'
        lname_jls(k) = 'H2O mixing ratio (weighted by daylight)'
        jls_ltop(k)  = LTOP
        jls_power(k) = -4
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'parts/vol')
c
        k = k + 1
        jls_day=k
        sname_jls(k) = 'daylight'   ! not output
        lname_jls(k) = 'Daylight weighting'
        jls_ltop(k)  = 1
        jls_power(k) = 0
        scale_jls(k) = 100.
        units_jls(k) = unit_string(jls_power(k),'%')
c
        k = k + 1
        jls_N2O5sulf=k
        sname_jls(k) = 'N2O5_sulf'
        lname_jls(k) = 'N2O5 sulfate sink'
        jls_ltop(k)  = LTOP
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c
        k = k + 1
        jls_O3vmr=k
        sname_jls(k) = 'O3_VMR'
        lname_jls(k) = 'O3 volume mixing ratio'
        jls_ltop(k)  = LTOP
        jls_power(k) = -8 ! simply to match Ox_CONCENTRATION diag
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'V/V air')

#endif  /* TRACERS_SPECIAL_Shindell */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
c Oxidants
#ifndef TRACERS_SPECIAL_Shindell
        k = k + 1
        jls_OHconk = k
        sname_jls(k) = 'OH_conc'
        lname_jls(k) = 'OH Concentration'
        jls_ltop(k) = LM
        jls_power(k) =5
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'molecules cm-3')
#endif
        k = k + 1
        jls_HO2con = k
        sname_jls(k) = 'HO2_conc'
        lname_jls(k) = 'HO2 Concentration'
        jls_ltop(k) =LM
        jls_power(k) =7
        scale_jls(k) =1.
        units_jls(k) = unit_string(jls_power(k),'molecules cm-3')

        k = k + 1
        jls_NO3 = k
        sname_jls(k) = 'NO3_conc'
        lname_jls(k) = 'NO3 Concentration'
        jls_ltop(k) =LM
        jls_power(k) =5
        scale_jls(k) =1.
        units_jls(k) = unit_string(jls_power(k),'molecules cm-3')
#endif  /* TRACERS_AEROSOLS_Koch || TRACERS_AMP || TRACERS_TOMAS */

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      k = k + 1
      jls_spec(nDustEv1jl)=k
      lname_jls(k)='No. dust events'
      sname_jls(k)='no_dust_ev1'
      jls_ltop(k)=1
      scale_jls(k)=SECONDS_PER_DAY/Dtsrc
      units_jls(k)='1/d'
      k = k + 1
      jls_spec(nDustEv2jl)=k
      lname_jls(k)='No. dust events above threshold wind'
      sname_jls(k)='no_dust_ev2'
      jls_ltop(k)=1
      scale_jls(k)=SECONDS_PER_DAY/Dtsrc
      units_jls(k)='1/d'
      k = k + 1
      jls_spec(nDustWthjl)=k
      lname_jls(k)='Threshold velocity for dust emission'
      sname_jls(k)='wtrsh'
      jls_ltop(k)=1
      scale_jls(k)=1.
      units_jls(k)='m s-1'
#endif

      if (k.gt. ktajls) then
        if (AM_I_ROOT()) write (6,*)
     &   'tjl_defs: Increase ktajls=',ktajls,' to at least ',k
         print *,'k should be ',k
        call stop_model('ktajls too small',255)
      end if
#endif /* TRACERS_ON */

      return

      contains

c      subroutine layer1_init_jls(k,n, name)
c      integer, intent(inout) :: k
c      integer, intent(in) :: n
c      character(len=*), intent(in) :: name
c      k = k + 1
c      jls_source(1,n) = k
c      sname_jls(k) = 'Layer_1_source_of_'//trim(trname(n))
c      lname_jls(k) = trim(trname(n))//' GRID SOURCE, LAYER 1'
c      jls_ltop(k) = 1
c      jls_power(k) = -3
c      units_jls(k) = unit_string(jls_power(k),'kg s-1')
c      end subroutine layer1_init_jls

      subroutine CO2n_init_jls(k,n,name)
      integer, intent(inout) :: k
      integer, intent(in) :: n
      character(len=*), intent(in) :: name
      k = k + 1
      jls_isrc(1,n) = k
      sname_jls(k) = 'Ocean_Gas_Exchange_'//trim(trname(n))
      lname_jls(k) = trim(trname(n))//' Ocean/Atmos. Gas Exchange'
      jls_ltop(k) = 1
      jls_power(k) = 3
      units_jls(k) = unit_string(jls_power(k),'kg s-1')
      end subroutine CO2n_init_jls

      subroutine Rn222_init_jls(k,n,name)
      integer, intent(inout) :: k
      integer, intent(in) :: n
      character(len=*), intent(in) :: name
      k = k + 1
      jls_decay(n) = k          ! special array for all radioactive sinks
      sname_jls(k) = 'Decay_of_'//trim(trname(n))
      lname_jls(k) = 'LOSS OF '//trim(trname(n))//' BY DECAY'
      jls_ltop(k) = lm
      jls_power(k) = -26
      units_jls(k) = unit_string(jls_power(k),'kg mb-1 m-2 s-1')
      jwt_jls(k)=3

c      k = k + 1
c      jls_source(1,n) = k
c      sname_jls(k) = 'Ground_Source_of_'//trim(trname(n))
c      lname_jls(k) = 'RADON-222 SOURCE, LAYER 1'
c      jls_ltop(k) = 1
c      jls_power(k) = -10
c      units_jls(k) = unit_string(jls_power(k),'kg s-1')
      end subroutine Rn222_init_jls

c      subroutine CO2_init_jls(k,n,name)
c      integer, intent(inout) :: k
c      integer, intent(in) :: n
c      character(len=*), intent(in) :: name
c        k = k + 1
c        jls_source(1,n) = k
c        sname_jls(k) = 'Fossil_fuel_source_'//trim(trname(n))
c        lname_jls(k) = 'CO2 Fossil fuel source (Marland)'
c        jls_ltop(k) = 1
c        jls_power(k) = 3
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(2,n) = k
c        sname_jls(k) = 'fertilization_sink_'//trim(trname(n))
c        lname_jls(k) = 'CO2 fertilization sink (Friedlingstein)'
c        jls_ltop(k) = 1
c        jls_power(k) = 3
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(3,n) = k
c        sname_jls(k) = 'Northern_forest_regrowth_'//trim(trname(n))
c        lname_jls(k) = 'CO2 Northern forest regrowth sink'
c        jls_ltop(k) = 1
c        jls_power(k) = 3
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(4,n) = k
c        sname_jls(k) = 'Land_Use_Modification_'//trim(trname(n))
c        lname_jls(k) = 'CO2 from Land use modification (Houton)'
c        jls_ltop(k) = 1
c        jls_power(k) = 3
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(5,n) = k
c        sname_jls(k) = 'Ecosystem_exchange_'//trim(trname(n))
c        lname_jls(k) = 'CO2 Ecosystem exchange (Matthews)'
c        jls_ltop(k) = 1
c        jls_power(k) = 3
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c        k = k + 1
c        jls_source(6,n) = k
c        sname_jls(k) = 'Ocean_exchange_'//trim(trname(n))
c        lname_jls(k) = 'CO2 Ocean exchange'
c        jls_ltop(k) = 1
c        jls_power(k) = 3
c        units_jls(k) = unit_string(jls_power(k),'kg s-1')
c
c      end subroutine CO2_init_jls

      subroutine N2O_init_jls(k,n,name)
      integer, intent(inout) :: k
      integer, intent(in) :: n
      character(len=*), intent(in) :: name
#ifdef TRACERS_SPECIAL_Shindell
      k = k + 1
      jls_3Dsource(nChemistry,n) = k
      sname_jls(k) = 'chemistry_source_of_'//trim(trname(n))
      lname_jls(k) = 'CHANGE OF '//trim(trname(n))//' BY CHEMISTRY'
      jls_ltop(k) = LM
      jls_power(k) = -1
      units_jls(k) = unit_string(jls_power(k),'kg s-1')
      k = k + 1
      jls_3Dsource(nOverwrite,n) = k
      sname_jls(k) = 'overwrite_source_of_'//trim(trname(n))
      lname_jls(k) =
     &     'CHANGE OF '//trim(trname(n))//' BY OVERWRITE'
      jls_ltop(k) = 1           ! really L=1 overwrite only
      jls_power(k) = -1
      units_jls(k) = unit_string(jls_power(k),'kg s-1')
#endif
#ifdef TRACERS_SPECIAL_Lerner
c      k = k + 1
c      jls_source(1,n) = k
c      sname_jls(k) = 'L1_sink_'//trim(trname(n))
c      lname_jls(k) = 'CHANGE OF '//trim(trname(n))//
c     &               ' BY RESETTING TO 462.2d-9, L1'
c      jls_ltop(k) = 1
c      jls_power(k) = 0
c      units_jls(k) = unit_string(jls_power(k),'kg s-1')
      k = k + 1
      jls_3Dsource(1,n) = k
      sname_jls(k) = 'Stratos_chem_change_'//trim(trname(n))
      lname_jls(k) = 'CHANGE OF '//trim(trname(n))//
     &               ' BY CHEMISTRY IN STRATOS'
      jls_ltop(k) = lm
      jls_power(k) = -1
      units_jls(k) = unit_string(jls_power(k),'kg s-1')
#endif
      end subroutine N2O_init_jls

      end subroutine init_jls_diag

      subroutine init_ijts_diag
!@sum init_ijts_diag Initialise lat/lon tracer diags
!@auth Gavin Schmidt
      use Tracer_mod, only: Tracer
      use TracerSurfaceSource_mod, only: TracerSurfaceSource
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE MODEL_COM, only: dtsrc
      use TRACER_COM, only: ntm, n_SO2, naircraft, nbiomass, nchemistry
      use TRACER_COM, only: nOther, nOverwrite, nVolcanic, nChemloss
      use TRACER_COM, only: ntsurfsrc, tracers, do_aircraft, aqchem_list
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: n_AOCOB, n_ANUM, n_AECOB
#endif
      USE DIAG_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      use trdust_mod, only: nDustEmij, nDustEm2ij, nDustEv1ij
     &   ,nDustEv2ij, nDustWthij, imDust, nSubClays
#endif
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
      USE CLOUDS, ONLY : diag_wetdep
#endif
      use RAD_COM, only: diag_fc
#endif /* TRACERS_ON */
#ifdef TRACERS_AMP
      use tracer_com, only: n_N_AKK_1
#endif
      use OldTracer_mod, only: has_chemistry
      use OldTracer_mod, only: trname, ntm_power, dodrydep,
     &          src_dist_index,nBBsources,do_fire
      use rad_com, only: nradfrc
      implicit none

      interface
        subroutine set_diag_aod(n,k,n_subclasses)
        integer, intent(inout) :: k
        integer, intent(in) :: n
        integer, optional, intent(in) :: n_subclasses
        end subroutine set_diag_aod
      end interface
      interface
        subroutine set_diag_rf(n,k,n_subclasses)
        integer, intent(inout) :: k
        integer, intent(in) :: n
        integer, optional, intent(in) :: n_subclasses
        end subroutine set_diag_rf
      end interface

      integer k,n,n1,kr,ktaijs_out,n_src
      character*50 :: unit_string
      CHARACTER*17 :: cform
      class (Tracer), pointer :: pTracer
      type (TracerSurfaceSource), pointer :: sources(:)

#ifdef TRACERS_ON
C**** Defaults for ijts (sources, sinks, etc.)
      ijts_fc(:,:)=0
      ijts_3Dsource(:,:)=0
      ijts_aq(:)=0
      ijts_isrc(:,:)=0
      ijts_gasex(:,:)=0
      ijts_HasArea(:) = .true. ! default applies to >50% of cases ???
      denom_ijts(:) = 0
#ifdef TRACERS_AMP
      ijts_AMPe(:)=0
      ijts_AMPp(:,:)=0
      ijts_AMPpdf(:,:)=0
#endif
#ifdef TRACERS_TOMAS
      ijts_TOMAS(:,:)=0
      ijts_subcoag(:)=0
#endif
C**** This needs to be 'hand coded' depending on circumstances
      k = 0
      do n=1,NTM
        if (src_dist_index(n)/=0) cycle

!=============================================!
! emissions for all tracers, if they have any !
!=============================================!

! handle exceptions first (e.g. SO4 emissions are listed under SO2)
      select case (trname(n))
      case ('SO4',
     &      'M_AKK_SU','M_ACC_SU',
     &      'ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     &      'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     &      'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15')
        n_src = n_SO2
#ifdef TRACERS_TOMAS
      case ('AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05',
     &      'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10',
     &      'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15',
     &      'AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05',
     &      'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10',
     &      'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15')
        n_src = n_AECOB(1)
      case ('AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05',
     &      'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10',
     &      'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15',
     &      'AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05',
     &      'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10',
     &      'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15')
        n_src = n_AOCOB(1)
#endif  /* TRACERS_TOMAS */
      case default
        n_src = n
      end select
      pTracer => tracers%getReference(trname(n_src))
      sources => pTracer%surfaceSources

! aqueous chemistry sources and sinks
      if (allocated(aqchem_list)) then
      if (any(n.eq.aqchem_list)) then
        k = k + 1
        ijts_aq(n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' aqueous chemistry change'
        sname_ijts(k) = trim(trname(n))//'_aqchem'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      endif
      endif

! surface emissions
      do kr=1,ntsurfsrc(n_src)
        k = k+1
        ijts_source(kr,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = trim(trname(n))//'_'//
     &                  trim(sources(kr)%sourceName)
        lname_ijts(k) = trim(trname(n))//' '//
     &                  trim(sources(kr)%sourceLname)
        ijts_power(k) = -15
        select case(trname(n))
        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     *    'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     *    'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15')
          units_ijts(k) = unit_string(ijts_power(k),'# m-2 s-1')
        case default
          units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        end select
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      end do

! aircraft emissions
      if(do_aircraft(n_src))then
        k = k + 1
        ijts_3Dsource(nAircraft,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = trim(trname(n))//'_aircraft_src'
        lname_ijts(k) = trim(trname(n))//' aircraft source'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      end if

! biomass burning emissions
      if (nBBsources(n_src) .gt. 0 .or. do_fire(n_src)) then
        k = k + 1
        ijts_3Dsource(nBiomass,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = trim(trname(n))//'_biomass_src'
        lname_ijts(k) = trim(trname(n))//' biomass source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      endif

!============================================!
! Chemical source (+) or sink (-) of tracers !
!============================================!
      if (has_chemistry(n)) then
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Chemistry'
        sname_ijts(k) = trim(trname(n))//'_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      endif

!=============================!
! Tracer-specific diagnostics !
!=============================!
      select case (trname(n))

      case ('CFCn')
      k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Ocean source'
        sname_ijts(k) = trim(trname(n))//'_Ocean_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      k = k+1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Layer 1 SOURCE'
        sname_ijts(k) = trim(trname(n))//'_CFC-GRID_SOURCE_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      k = k+1  ! Gas Exchange Coefficient (piston velocity) (open ocean only)
        ijts_gasex(1,n)  = k
        ia_ijts(k) = ia_srf
        sname_ijts(k) = 'Piston_Veloc_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Piston Velocity'
        ijtc_power(n) = -5
        units_ijts(k) = unit_string(ijtc_power(n),'m s-1')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      k = k+1  ! Gas Exchange Solubility coefficient
        ijts_gasex(2,n)  = k
        ia_ijts(k) = ia_srf
        sname_ijts(k) = 'Solubility_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Solubility'
        ijtc_power(n) = -5
        units_ijts(k) = unit_string(ijtc_power(n),'mol/m3/uatm')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      k = k+1  ! Gas exchange
        ijts_gasex(3,n)  = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Gas_Exchange_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Gas Exchange'
        ijtc_power(n) = 0
        units_ijts(k) = unit_string(ijtc_power(n),'molCFC/m2/yr')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      case ('CO2n')
      k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'CO2_GASX'
        lname_ijts(k) = 'AO GASEX CO2'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg,CO2/m2/s')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      k = k+1  ! Gas Exchange Coefficient (piston velocity) (open ocean only)
        ijts_gasex(1,n)  = k
        ia_ijts(k) = ia_srf
        sname_ijts(k) = 'Piston_Veloc_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Piston Velocity'
        ijtc_power(n) = -5
        units_ijts(k) = unit_string(ijtc_power(n),'m s-1')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      k = k+1  ! Gas Exchange Solubility coefficient
        ijts_gasex(2,n)  = k
        ia_ijts(k) = ia_srf
        sname_ijts(k) = 'Solubility_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Solubility'
        ijtc_power(n) = -5
        units_ijts(k) = unit_string(ijtc_power(n),'mol/m3/uatm')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      k = k+1  ! Gas exchange
        ijts_gasex(3,n)  = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Gas_Exchange_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Gas Exchange'
        ijtc_power(n) = 0
        units_ijts(k) = unit_string(ijtc_power(n),'molCO2/m2/yr')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      case ('SF6','SF6_c')
      k = k+1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Layer 1 SOURCE'
        sname_ijts(k) = trim(trname(n))//'_CFC-GRID_SOURCE_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('nh5','nh50','nh15')

      k = k+1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Layer 1 SOURCE'
        sname_ijts(k) = trim(trname(n))//'_NH_Mid_SOURCE_LAYER_1'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('e90')

      k = k+1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Layer 1 SOURCE'
        sname_ijts(k) = trim(trname(n))//'_SOURCE_LAYER_1'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('tape_rec')

      k = k+1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' UTLS_source'
        sname_ijts(k) = trim(trname(n))//'_SOURCE_UTLS'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('aoanh','aoa')

      k = k+1
        ijts_3Dsource(nOverwrite,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' overwrite'
        sname_ijts(k) = trim(trname(n))//'_overw'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'days')
        scale_ijts(k) = 10.**(-ijts_power(k))

      case ('Rn222')
      k = k+1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' L 1 SOURCE'
        sname_ijts(k) = trim(trname(n))//'_SOURCE_Layer_1'
        ijts_power(k) = -21
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('CO2')
! keep AIJ and AJL CO2 sources in same order !!
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Fossil_fuel_source_'//trim(trname(n))
        lname_ijts(k) = trim(trname(n))//' Fossil fuel src'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(2,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'fertilization_sink_'//trim(trname(n))
        lname_ijts(k) = trim(trname(n))//' fertilization'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(3,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Northern_forest_regrowth_'//trim(trname(n))
        lname_ijts(k) = trim(trname(n))//' North forest regrowth'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(4,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Land_Use_Modification_'//trim(trname(n))
        lname_ijts(k) = trim(trname(n))//' from Land use mods'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(5,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Ecosystem_exchange_'//trim(trname(n))
        lname_ijts(k) = trim(trname(n))//' Ecosystem exch'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(6,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Ocean_exchange_'//trim(trname(n))
        lname_ijts(k) = trim(trname(n))//' Ocean exchange'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('N2O')
#ifdef TRACERS_SPECIAL_Shindell
      k = k + 1
        ijts_3Dsource(nOverwrite,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Overwrite'
        sname_ijts(k) = trim(trname(n))//'_overw'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
#ifdef TRACERS_SPECIAL_Lerner
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' CHANGE IN L 1'
        sname_ijts(k) = trim(trname(n))//'_CHANGE_IN_L_1'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('CFC11')
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' L 1 SOURCE'
        sname_ijts(k) = trim(trname(n))//'_SOURCE_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('14CO2')
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' L 1 Sink'
        sname_ijts(k) = trim(trname(n))//'_L1_Sink'
        ijts_power(k) = -21
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('NOx','CO','Isoprene','Alkenes','Paraffin',
#ifdef TRACERS_dCO
     *'d13Calke','d13CPAR',
     *'d17OPAN', 'd18OPAN', 'd13CPAN',
     *'dMe17OOH', 'dMe18OOH', 'd13MeOOH',
     *'dHCH17O', 'dHCH18O', 'dH13CHO',
     *'dC17O', 'dC18O', 'd13CO',
#endif  /* TRACERS_dCO */
     &'ClOx','BrOx','HCl','HOCl','ClONO2','HBr','HOBr','BrONO2',
     &'CFC','H2O2','CH3OOH','Ox','N2O5','HNO3','HCHO','Terpenes',
     &'HO2NO2','PAN','AlkylNit','stratOx')

        select case(trname(n))
        case('NOx','CO','Isoprene','Alkenes','Paraffin',
#ifdef TRACERS_dCO
     *  'd13Calke','d13CPAR',
     *  'd17OPAN', 'd18OPAN', 'd13CPAN',
     *  'dMe17OOH', 'dMe18OOH', 'd13MeOOH',
     *  'dHCH17O', 'dHCH18O', 'dH13CHO',
     *  'dC17O', 'dC18O', 'd13CO',
#endif  /* TRACERS_dCO */
     &  'CFC','H2O2','CH3OOH','Ox','N2O5','HNO3','HCHO',
     &  'Terpenes','HO2NO2','PAN','AlkylNit','stratOx')
          k = k + 1
          ijts_3Dsource(nOverwrite,n) = k
          ia_ijts(k) = ia_src
          lname_ijts(k) = trim(trname(n))//' Overwrite'
          sname_ijts(k) = trim(trname(n))//'_overw'
          ijts_power(k) = -12
          units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        select case(trname(n))
        case('NOx')
          k = k + 1
          ijts_3Dsource(nOther,n) = k
          ia_ijts(k) = ia_src
          lname_ijts(k) = trim(trname(n))//' Lightning Source'
          sname_ijts(k) = trim(trname(n))//'_lightning'
          ijts_power(k) = -12
          units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        case('Ox','stratOx')
          if (nradfrc>0) then
            k = k + 1
            ijts_fc(1,n) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//' tropopause SW rad forc'
            sname_ijts(k) = 'swf_tp_'//trim(trname(n))
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W m-2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
            k = k + 1
            ijts_fc(2,n) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//' tropopause LW rad forc'
            sname_ijts(k) = 'lwf_tp_'//trim(trname(n))
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W m-2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
            k = k + 1
            ijts_fc(3,n) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//' TOA SW rad forc'
            sname_ijts(k) = 'swf_toa_'//trim(trname(n))
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W m-2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
            k = k + 1
            ijts_fc(4,n) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//' TOA LW rad forc'
            sname_ijts(k) = 'lwf_toa_'//trim(trname(n))
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W m-2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
#ifdef AUX_OX_RADF_TROP
#ifndef AUXILIARY_OX_RADF
            call stop_model
     &      ('AUX_OX_RADF_TROP needs AUXILIARY_OX_RADF',255)
#endif
#endif
#ifdef AUXILIARY_OX_RADF
            if(trname(n)=='Ox')then
              k = k + 1
              ijts_auxfc(1) = k
              ia_ijts(k) = ia_rad_frc
              lname_ijts(k) = trim(trname(n))//' AUX tropp SW rad forc'
              sname_ijts(k) = 'swfauxtp_'//trim(trname(n))
              ijts_power(k) = -2
              units_ijts(k) = unit_string(ijts_power(k),'W m-2')
              scale_ijts(k) = 10.**(-ijts_power(k))
              ijts_HasArea(k) = .false.
              k = k + 1
              ijts_auxfc(2) = k
              ia_ijts(k) = ia_rad_frc
              lname_ijts(k) = trim(trname(n))//' AUX tropp LW rad forc'
              sname_ijts(k) = 'lwfauxtp_'//trim(trname(n))
              ijts_power(k) = -2
              units_ijts(k) = unit_string(ijts_power(k),'W m-2')
              scale_ijts(k) = 10.**(-ijts_power(k))
              ijts_HasArea(k) = .false.
              k = k + 1
              ijts_auxfc(3) = k
              ia_ijts(k) = ia_rad_frc
              lname_ijts(k) = trim(trname(n))//' AUX TOA SW rad forc'
              sname_ijts(k) = 'swfauxtoa_'//trim(trname(n))
              ijts_power(k) = -2
              units_ijts(k) = unit_string(ijts_power(k),'W m-2')
              scale_ijts(k) = 10.**(-ijts_power(k))
              ijts_HasArea(k) = .false.
              k = k + 1
              ijts_auxfc(4) = k
              ia_ijts(k) = ia_rad_frc
              lname_ijts(k) = trim(trname(n))//' AUX TOA LW rad forc'
              sname_ijts(k) = 'lwfauxtoa_'//trim(trname(n))
              ijts_power(k) = -2
              units_ijts(k) = unit_string(ijts_power(k),'W m-2')
              scale_ijts(k) = 10.**(-ijts_power(k))
              ijts_HasArea(k) = .false.
            endif
#endif /* AUXILIARY_OX_RADF */
          endif
#ifdef ACCMIP_LIKE_DIAGS
          if(trname(n)=='Ox' .and. dodrydep(n))then
            k = k+1
            ijts_Sdrydep = k
            ia_ijts(k) = ia_src
            lname_ijts(k) = trim(trname(n))//' stomatal drydep flux'
            sname_ijts(k) = 'stomatal_'//trim(trname(n))
            ijts_power(k) = ntm_power(n)-4
            units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
            scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
            ijts_HasArea(k) = .false.
          end if
#endif /* ACCMIP_LIKE_DIAGS */
        end select
      end select

      case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
        k = k + 1
        ijts_3Dsource(nOverwrite,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Overwrite'
        sname_ijts(k) = trim(trname(n))//'_overw'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#else
      k = k + 1
        ijts_source(6,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 sink due to soil absorp.'
        sname_ijts(k) = 'CH4_soil_sink.'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(7,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Termite source'
        sname_ijts(k) = 'CH4_Termite_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(9,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Ocean source'
        sname_ijts(k) = 'CH4_Ocean_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(10,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Fresh Water lake source'
        sname_ijts(k) = 'CH4_lake_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(11,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Misc. Ground source'
        sname_ijts(k) = 'CH4_Misc._Ground_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(14,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Wetlands+Tundra source'
        sname_ijts(k) = 'CH4_Wetlands+Tundra_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Animal source'
        sname_ijts(k) = 'CH4_Animal_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal Mine source'
        sname_ijts(k) = 'CH4_Coal_Mine_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(3,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Leak source'
        sname_ijts(k) = 'CH4_Gas_Leak_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(4,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Venting source'
        sname_ijts(k) = 'CH4_Gas_Venting_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(5,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Municipal solid waste src'
        sname_ijts(k) = 'CH4_MSW_src'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(8,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal combustion source'
        sname_ijts(k) = 'CH4_Coal_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(12,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Biomass burning source'
        sname_ijts(k) = 'CH4_Biomass_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(13,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Rice cultivation source'
        sname_ijts(k) = 'CH4_Rice_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Tropospheric Chemistry'
        sname_ijts(k) = 'CH4_trop_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Stratospheric Chemistry'
        sname_ijts(k) = 'CH4_strat_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('O3')
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 deposition, layer 1'
        sname_ijts(k) = 'O3_deposition_L1'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 Stratospheric Chem.'
        sname_ijts(k) = 'O3_strat_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 Tropo. Chem. Production'
        sname_ijts(k) = 'O3_trop_chem_prod'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(3,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 Tropo. Chemistry Loss'
        sname_ijts(k) = 'O3_trop_chem_loss'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

#ifdef TRACERS_WATER
      case ('Water', 'H2O18', 'H2O17', 'HDO', 'HTO' )
          ! nothing I can think of....
#endif

#ifdef SHINDELL_STRAT_EXTRA
      case ('GLT')
      k = k+1
        ijts_3Dsource(nOverwrite,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' overwrite'
        sname_ijts(k) = trim(trname(n))//'_overw'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('BCB', 'OCB', 'BCIA', 'OCIA', 'NO3p', 'isopp1a')
        call set_diag_aod(n,k)
        if (diag_fc==2) call set_diag_rf(n,k)

      case ('SO2')
c production from volcanic emissions
        k = k + 1
        ijts_3Dsource(nVolcanic,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = trim(trname(n))//'_volcanic_src'
        lname_ijts(k) = trim(trname(n))//' volcanic source'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c production from DMS
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' source from DMS'
        sname_ijts(k) = trim(trname(n))//'_source_from_DMS'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c chemical loss
        k = k + 1
        ijts_3Dsource(nChemloss,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Chemical sink'
        sname_ijts(k) = trim(trname(n))//'_chem_sink'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &      'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' partitioning'
        sname_ijts(k) = trim(trname(n))//'_partitioning'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        select case(trname(n))
        case ('vbsAm2')
          call set_diag_aod(n,k)
          if (diag_fc==2) call set_diag_rf(n,k)
        end select

      case ('DMS')
        k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = trim(trname(n))//'_ocean_src'
        lname_ijts(k) = trim(trname(n))//' ocean source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('SO4')
c put in production of SO4 from gas phase
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' gas phase source'
        sname_ijts(k) = trim(trname(n))//'_gas_phase_source'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(nVolcanic,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = trim(trname(n))//'_volcanic_src'
        lname_ijts(k) = trim(trname(n))//' volcanic source'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        call set_diag_aod(n,k)
        if (diag_fc==2) call set_diag_rf(n,k)

      case ('vbsGm2', 'vbsGm1', 'vbsGz',  'vbsGp1', 'vbsGp2',
     &      'vbsGp3', 'vbsGp4', 'vbsGp5', 'vbsGp6')
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' aging source'
        sname_ijts(k) = trim(trname(n))//'_aging_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(nChemloss,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' aging loss'
        sname_ijts(k) = trim(trname(n))//'_aging_loss'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(nOther,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' partitioning'
        sname_ijts(k) = trim(trname(n))//'_partitioning'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

#ifdef TRACERS_AMP
      case ('M_NO3   ','M_NH4   ','M_H2O   ','M_AKK_SU','N_AKK_1 ',!AKK
     *    'M_ACC_SU','N_ACC_1 ','M_DD1_SU','M_DD1_DU','N_DD1_1 ',!ACC,DD1
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ','M_DD2_SU','M_DD2_DU',!DS1,DD2
     *    'N_DD2_1 ','M_DS2_SU','M_DS2_DU','N_DS2_1 ','M_SSA_SU',!DD2,DS2,SSA
     *    'M_SSA_SS','N_SSA_1 ','M_SSC_SS','N_SSC_1',            !SSA,SSC
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',!OCC,BC1
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_BC3_SU',!BC1,BC2,BC3
     *    'M_BC3_BC','N_BC3_1 ','M_DBC_SU','M_DBC_BC','M_DBC_DU',!BC3,DBC
     *    'N_DBC_1 ','M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',!DBC,BOC
     *    'M_BCS_SU','M_BCS_BC','N_BCS_1 ','M_MXX_SU','M_MXX_BC',!BCS,MXX
     *    'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','M_OCS_SU',
     *    'M_OCS_OC','N_OCS_1 ','M_SSS_SS','M_SSS_SU')
       k = k + 1
         ijts_3Dsource(nChemistry,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'AMP_src_'//trim(trname(n))
         sname_ijts(k) = 'AMP_src_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
       k = k + 1
         ijts_AMPp(1,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P1_Nucl_'//trim(trname(n))
         sname_ijts(k) = 'P1_Nucl_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(2,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P2_Coag_'//trim(trname(n))
         sname_ijts(k) = 'P2_Coag_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(3,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P3_Cond_'//trim(trname(n))
         sname_ijts(k) = 'P3_Cond_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(4,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P4_Incld_NIMC_'//trim(trname(n))
         sname_ijts(k) = 'P4_Incld_NIMC_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(5,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P5_IMLoss_NIAC_'//trim(trname(n))
         sname_ijts(k) = 'P5_IMLoss_NIAC_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(6,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P6_Mode_Trans_'//trim(trname(n))
         sname_ijts(k) = 'P6_Mode_Trans_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(7,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P7_Total_Change_'//trim(trname(n))
         sname_ijts(k) = 'P7_Total_Change_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
#endif
#ifdef TRACERS_TOMAS

      case ('SOAgas')

c put in production of SO4 from gas phase
        do kr=1,ntsurfsrc(n)
          k = k + 1
          ijts_source(kr,n) = k
          ia_ijts(k) = ia_src
          sname_ijts(k) = trim(trname(n))//'_terpenes_src'
          lname_ijts(k) = trim(trname(n))//' terpenes source'
          ijts_power(k) = -15
          units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        enddo

      case ('H2SO4')

c put in production of SO4 from gas phase
        k = k + 1
        ijts_3Dsource(nOther,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Microphysics change '//trim(trname(n))
        sname_ijts(k) = 'Microphysics_chg_'//trim(trname(n))
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case('ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     *    'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     *    'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15',
     *    'ANACL_01','ANACL_02','ANACL_03','ANACL_04','ANACL_05',
     *    'ANACL_06','ANACL_07','ANACL_08','ANACL_09','ANACL_10',
     *    'ANACL_11','ANACL_12','ANACL_13','ANACL_14','ANACL_15',
     *    'AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05',
     *    'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10',
     *    'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15',
     *    'AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05',
     *    'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10',
     *    'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15',
     *    'AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05',
     *    'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10',
     *    'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15',
     *    'AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05',
     *    'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10',
     *    'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15',
     *    'ADUST_01','ADUST_02','ADUST_03','ADUST_04','ADUST_05',
     *    'ADUST_06','ADUST_07','ADUST_08','ADUST_09','ADUST_10',
     *    'ADUST_11','ADUST_12','ADUST_13','ADUST_14','ADUST_15',
     *    'ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     *    'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     *    'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15')

       k = k + 1
         ijts_3Dsource(nOther,n) = k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'Microphysics change '//trim(trname(n))
         sname_ijts(k) = 'Microphysics_chg_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc


       k = k + 1
         ijts_TOMAS(1,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'MP1_Cond_'//trim(trname(n))
         sname_ijts(k) = 'MP1_Cond_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
       k = k + 1
         ijts_TOMAS(2,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'MP2_Coag_'//trim(trname(n))
         sname_ijts(k) = 'MP2_Coag_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
       k = k + 1
         ijts_TOMAS(3,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'MP3_Nucl_'//trim(trname(n))
         sname_ijts(k) = 'MP3_Nucl_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
       k = k + 1
         ijts_TOMAS(4,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'MP4_Aqoxid_MC_'//trim(trname(n))
         sname_ijts(k) = 'MP4_Aqoxid_MC_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
       k = k + 1
         ijts_TOMAS(5,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'MP5_Aqoxid_LS_'//trim(trname(n))
         sname_ijts(k) = 'MP5_Aqoxid_LS_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
       k = k + 1
         ijts_TOMAS(6,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'MP6_Mk_Nk_Fix_'//trim(trname(n))
         sname_ijts(k) = 'MP6_Mk_Nk_Fix_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
       k = k + 1
         ijts_TOMAS(7,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'MP7_Aeroupdate_'//trim(trname(n))
         sname_ijts(k) = 'MP7_Aeroupdate_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        k = k + 1
        ijts_subcoag(n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Subgrid_coag_'//trim(trname(n))
        sname_ijts(k) = 'Subgrid_coag_'//trim(trname(n))
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc


         select case(trname(n))
        case ('ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     *       'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     *       'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15')

        k = k + 1
        ijts_3Dsource(nVolcanic,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = trim(trname(n))//'_volcanic_src'
        lname_ijts(k) = trim(trname(n))//' volcanic source'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     *    'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     *    'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15')

        k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4 source '//trim(trname(n))
        sname_ijts(k) = 'SO4_src_'//trim(trname(n))
        ijts_power(k) = 10
        units_ijts(k) = unit_string(ijts_power(k),'# m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        k = k + 1
        ijts_3Dsource(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'EC source '//trim(trname(n))
        sname_ijts(k) = 'EC_src_'//trim(trname(n))
        ijts_power(k) = 10
        units_ijts(k) = unit_string(ijts_power(k),'# m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
!TOMAS!#endif
        k = k + 1
        ijts_3Dsource(4,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'OC source '//trim(trname(n))
        sname_ijts(k) = 'OC_src_'//trim(trname(n))
        ijts_power(k) = 10
        units_ijts(k) = unit_string(ijts_power(k),'# m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NACL source '//trim(trname(n))
        sname_ijts(k) = 'NACL_src_'//trim(trname(n))
        ijts_power(k) = 10
        units_ijts(k) = unit_string(ijts_power(k),'# m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        k = k + 1
        ijts_isrc(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'DUST source '//trim(trname(n))
        sname_ijts(k) = 'DUST_src_'//trim(trname(n))
        ijts_power(k) = 10
        units_ijts(k) = unit_string(ijts_power(k),'# m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('ANACL_01','ANACL_02','ANACL_03','ANACL_04','ANACL_05',
     *    'ANACL_06','ANACL_07','ANACL_08','ANACL_09','ANACL_10',
     *    'ANACL_11','ANACL_12','ANACL_13','ANACL_14','ANACL_15')

        k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = trim(trname(n))//'_emission'
        lname_ijts(k) = trim(trname(n))//' Ocean source'

        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        case('ADUST_01','ADUST_02','ADUST_03','ADUST_04','ADUST_05',
     *    'ADUST_06','ADUST_07','ADUST_08','ADUST_09','ADUST_10',
     *    'ADUST_11','ADUST_12','ADUST_13','ADUST_14','ADUST_15')

        k=k+1
        ijts_isrc(1,n)=k
        lname_ijts(k)=trim(trname(n))//' Emission'
        sname_ijts(k)=trim(trname(n))//'_emission'
        ia_ijts(k)=ia_src
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        end select

        select case(trname(n))

        case('ASO4__01','ANACL_01','AECOB_01','AECIL_01',
     &       'AOCOB_01','AOCIL_01','ADUST_01')

        call set_diag_aod(n,k)
        if (diag_fc==2) then
          call set_diag_rf(n,k)
        else if (diag_fc==1) then
          select case (trname(n))
            case ('ASO4__01')
              call set_diag_rf(n,k)
          end select
        endif

      end select

#endif
      case ('H2O2_s')
c put in production of H2O2 from gas phase
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' gas phase source'
        sname_ijts(k) = trim(trname(n))//'_gas_phase_source'
        ijts_power(k) = -10
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in production of H2O2 from gas phase
        k = k + 1
        ijts_3Dsource(nChemLoss,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' gas phase sink'
        sname_ijts(k) = trim(trname(n))//'_gas_phase_sink'
        ijts_power(k) = -10
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('seasalt1', 'seasalt2', 'OCocean')
        k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = trim(trname(n))//'_ocean_src'
        lname_ijts(k) = trim(trname(n))//' ocean source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_SEASALT
        select case (trname(n))
        case ('seasalt1')
          call set_diag_aod(n,k)
          if (diag_fc>0) call set_diag_rf(n,k) ! this handles [sl]wf_OMA if =1
        case ('seasalt2')
          call set_diag_aod(n,k)
        end select

#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      CASE('Clay','Silt1','Silt2','Silt3','Silt4','Silt5','ClayIlli'
     &      ,'ClayKaol','ClaySmec','ClayCalc','ClayQuar','ClayFeld'
     &      ,'ClayHema','ClayGyps','ClayIlHe','ClayKaHe','ClaySmHe'
     &      ,'ClayCaHe','ClayQuHe','ClayFeHe','ClayGyHe','Sil1Quar'
     &      ,'Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps','Sil1Illi'
     &      ,'Sil1Kaol','Sil1Smec','Sil1QuHe','Sil1FeHe','Sil1CaHe'
     &      ,'Sil1GyHe','Sil1IlHe','Sil1KaHe','Sil1SmHe','Sil2Quar'
     &      ,'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps','Sil2Illi'
     &      ,'Sil2Kaol','Sil2Smec','Sil2QuHe','Sil2FeHe','Sil2CaHe'
     &      ,'Sil2GyHe','Sil2IlHe','Sil2KaHe','Sil2SmHe','Sil3Quar'
     &      ,'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps','Sil3Illi'
     &      ,'Sil3Kaol','Sil3Smec','Sil3QuHe','Sil3FeHe','Sil3CaHe'
     &      ,'Sil3GyHe','Sil3IlHe','Sil3KaHe','Sil3SmHe','Sil4Quar'
     &      ,'Sil4Feld','Sil4Calc','Sil4Hema','Sil4Gyps','Sil4Illi'
     &      ,'Sil4Kaol','Sil4Smec','Sil4QuHe','Sil4FeHe','Sil4CaHe'
     &      ,'Sil4GyHe','Sil4IlHe','Sil4KaHe','Sil4SmHe','Sil5Quar'
     &      ,'Sil5Feld','Sil5Calc','Sil5Hema','Sil5Gyps','Sil5Illi'
     &      ,'Sil5Kaol','Sil5Smec','Sil5QuHe','Sil5FeHe','Sil5CaHe'
     &      ,'Sil5GyHe','Sil5IlHe','Sil5KaHe','Sil5SmHe')
        k=k+1
        ijts_isrc(nDustEmij,n)=k
        lname_ijts(k)='Emission of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_emission'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        IF ( imDust == 0 .or. imDust >= 3 ) THEN
        k=k+1
        ijts_isrc(nDustEm2ij,n)=k
        lname_ijts(k)='Cubic emission of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_emission2'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        END IF
#ifndef TRACERS_DRYDEP
      k=k+1
        ijts_isrc(nDustTurbij,n)=k
        lname_ijts(k)='Turbulent Deposition of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_turb_depo'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
#ifndef TRACERS_WATER
      k=k+1
        ijts_wet(n)=k
        lname_ijts(k)='Wet deposition of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_wet_depo'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
        SELECT CASE (trname(n))
        CASE ('Clay','ClayIlli','ClayKaol','ClaySmec','ClayCalc'
     &         ,'ClayQuar','ClayFeld' ,'ClayHema','ClayGyps','ClayIlHe'
     &         ,'ClayKaHe','ClaySmHe' ,'ClayCaHe','ClayQuHe','ClayFeHe'
     &         ,'ClayGyHe')

        call set_diag_aod(n,k,nSubClays)
        if (diag_fc==2) call set_diag_rf(n,k,nSubClays)

        CASE('Silt1','Silt2','Silt3','Silt4','Silt5','Sil1Quar'
     &         ,'Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps','Sil1Illi'
     &         ,'Sil1Kaol','Sil1Smec','Sil1QuHe','Sil1FeHe','Sil1CaHe'
     &         ,'Sil1GyHe','Sil1IlHe','Sil1KaHe','Sil1SmHe','Sil2Quar'
     &         ,'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps','Sil2Illi'
     &         ,'Sil2Kaol','Sil2Smec','Sil2QuHe','Sil2FeHe','Sil2CaHe'
     &         ,'Sil2GyHe','Sil2IlHe','Sil2KaHe','Sil2SmHe','Sil3Quar'
     &         ,'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps','Sil3Illi'
     &         ,'Sil3Kaol','Sil3Smec','Sil3QuHe','Sil3FeHe','Sil3CaHe'
     &         ,'Sil3GyHe','Sil3IlHe','Sil3KaHe','Sil3SmHe','Sil4Quar'
     &         ,'Sil4Feld','Sil4Calc','Sil4Hema','Sil4Gyps','Sil4Illi'
     &         ,'Sil4Kaol','Sil4Smec','Sil4QuHe','Sil4FeHe','Sil4CaHe'
     &         ,'Sil4GyHe','Sil4IlHe','Sil4KaHe','Sil4SmHe','Sil5Quar'
     &         ,'Sil5Feld','Sil5Calc','Sil5Hema','Sil5Gyps','Sil5Illi'
     &         ,'Sil5Kaol','Sil5Smec','Sil5QuHe','Sil5FeHe','Sil5CaHe'
     &         ,'Sil5GyHe','Sil5IlHe','Sil5KaHe','Sil5SmHe')

          call set_diag_aod(n,k)
          if (diag_fc==2) call set_diag_rf(n,k)

        END SELECT
#endif  /* TRACERS_DUST || TRACERS_MINERALS */

      end select

#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
c**** additional wet deposition diagnostics
      IF (diag_wetdep == 1) THEN
        k=k+1
        ijts_trdpmc(1,n)=k
        lname_ijts(k)='MC Condensation of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_cond_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(2,n)=k
        lname_ijts(k)='Evaporated '//TRIM(trname(n))
     &       //' in MC Downdrafts'
        sname_ijts(k)=TRIM(trname(n))//'_downeva_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(3,n)=k
        lname_ijts(k)='Condensed '//TRIM(trname(n))//' in MC CLW'
        sname_ijts(k)=TRIM(trname(n))//'_conclw_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(4,n)=k
        lname_ijts(k)='Precipitated '//TRIM(trname(n))//' by MC'
        sname_ijts(k)=TRIM(trname(n))//'_precip_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(5,n)=k
        lname_ijts(k)='Reevaporated '//TRIM(trname(n))
     &       //' from MC Precip'
        sname_ijts(k)=TRIM(trname(n))//'_reevap_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(6,n)=k
        lname_ijts(k)='MC Washout of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_washout_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(1,n)=k
        lname_ijts(k)='LS Washout of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_washout_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(2,n)=k
        lname_ijts(k)='Precipitated '//TRIM(trname(n))//' by LS'
        sname_ijts(k)=TRIM(trname(n))//'_precip_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(3,n)=k
        lname_ijts(k)='Condensed '//TRIM(trname(n))// ' in LS CLW'
        sname_ijts(k)=TRIM(trname(n))//'_conclw_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(4,n)=k
        lname_ijts(k)='Reevaporated '//TRIM(trname(n))
     &       //' from LS Precip'
        sname_ijts(k)=TRIM(trname(n))//'_reevap_ls'

        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(5,n)=k
        lname_ijts(k)='Evaporated '//TRIM(trname(n))//' from LS CLW'
        sname_ijts(k)=TRIM(trname(n))//'_clwevap_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(6,n)=k
        lname_ijts(k)='LS Condensation of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_cond_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      END IF
#endif
      end do

C**** Additional Special IJ diagnostics
C**** (not necessary associated with a particular tracer)
#ifdef BC_ALB
c BC impact on grain size
c         k = k + 1
c         ijts_alb(2,n) = k
c         ia_ijts(k) = ia_rad????
c         lname_ijts(k) = 'BC impact on grain size'
c         sname_ijts(k) = 'grain_BC'
c         ijts_power(k) = -9
c         units_ijts(k) = unit_string(ijts_power(k),' ')
c         scale_ijts(k) = 10.**(-ijts_power(k))
c BC impact on albedo
        if (nradfrc>0) then
          k = k + 1
          ijts_alb(1) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'BC impact on snow albedo of land/seaice'
          sname_ijts(k) = 'alb_BC'
          ijts_power(k) = 0
          units_ijts(k) = unit_string(ijts_power(k),'%')
          scale_ijts(k) = 100.
          ijts_HasArea(k) = .false.
          dname_ijts(k) = 'sunlit_snow_freq'

c SW forcing from albedo change
          k = k + 1
          ijts_alb(2) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'BCalb SW radiative forcing'
          sname_ijts(k) = 'swf_BCALB'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W m-2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
        endif

#endif
#ifdef TRACERS_SPECIAL_Shindell
#ifdef BIOGENIC_EMISSIONS
      k = k+1
        ijs_isoprene=k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Interactive isoprene source'
        sname_ijts(k) = 'Int_isop'
        ijts_power(k) = -10
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
      k = k + 1
        ijs_NO2_1030=k
        ia_ijts(k) = ia_src
        write(lname_ijts(k),'(a18)')'NO2 10:30 trop col'
        sname_ijts(k) = 'NO2_1030'
        dname_ijts(k) = 'NO2_1030c'
        ijts_power(k) = 15
        units_ijts(k) = unit_string(ijts_power(k),'molecules cm-2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      k = k + 1
        ijs_NO2_1030c=k
        ia_ijts(k) = ia_src ! overridden in TRACER_PRT...
        write(lname_ijts(k),'(a24)')'count NO2 10:30 trop col'
        write(sname_ijts(k),'(a9)')'NO2_1030c'
        ijts_power(k) = 0
        units_ijts(k) = unit_string(ijts_power(k),'number of accum')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      k = k + 1
        ijs_NO2_1330=k
        ia_ijts(k) = ia_src
        write(lname_ijts(k),'(a18)')'NO2 13:30 trop col'
        sname_ijts(k) = 'NO2_1330'
        dname_ijts(k) = 'NO2_1330c'
        ijts_power(k) = 15
        units_ijts(k) = unit_string(ijts_power(k),'molecules cm-2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      k = k + 1
        ijs_NO2_1330c=k
        ia_ijts(k) = ia_src ! overridden in TRACER_PRT...
        write(lname_ijts(k),'(a24)')'count NO2 13:30 trop col'
        write(sname_ijts(k),'(a9)')'NO2_1330c'
        ijts_power(k) = 0
        units_ijts(k) = unit_string(ijts_power(k),'number of accum')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      k = k + 1
        ijs_O3mass=k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Total Column Ozone (not Ox) Mass'
        sname_ijts(k) = 'O3_Total_Mass'
        ijts_power(k) = -4
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2') ! to match tracers
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
#endif  /* TRACERS_SPECIAL_Shindell */
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      k = k + 1
      ijts_spec(nDustEv1ij)=k
      lname_ijts(k)='No. dust events'
      sname_ijts(k)='no_dust_ev1'
      ia_ijts(k)=ia_src
      scale_ijts(k)=SECONDS_PER_DAY/Dtsrc
      units_ijts(k)='1/d'
      ijts_HasArea(k) = .false.
      k = k + 1
      ijts_spec(nDustEv2ij)=k
      lname_ijts(k)='No. dust events above threshold wind'
      sname_ijts(k)='no_dust_ev2'
      ia_ijts(k)=ia_src
      scale_ijts(k)=SECONDS_PER_DAY/Dtsrc
      units_ijts(k)='1/d'
      ijts_HasArea(k) = .false.
      k = k + 1
      ijts_spec(nDustWthij)=k
      lname_ijts(k)='Threshold velocity for dust emission'
      sname_ijts(k)='wtrsh'
      ia_ijts(k)=ia_src
      scale_ijts(k)=1.
      units_ijts(k)='m s-1'
      ijts_HasArea(k) = .false.
#endif

#ifdef TRACERS_AMP


      do n=1,NTM
        pTracer => tracers%getReference(trname(n))
        sources => pTracer%surfaceSources
      select case(trname(n))
        case('M_AKK_SU','M_ACC_SU','Water')
        k = k + 1
          ijts_3Dsource(nVolcanic,n)=k
          ia_ijts(k) = ia_src
          sname_ijts(k) = trim(trname(n))//'_volcanic_src'
          lname_ijts(k) = trim(trname(n))//' volcanic source'
          ijts_power(k) = -15
          units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c- interactive sources diagnostic
      CASE('M_DD1_DU','M_SSA_SS','M_SSC_SS','M_DD2_DU','M_SSS_SS')
      k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Emission_'//trim(trname(n))
        sname_ijts(k) = 'Emission_'//trim(trname(n))
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg m-2 s-1')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      CASE('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 ',
     *     'N_DS2_1 ','N_SSA_1 ','N_SSC_1 ','N_OCC_1 ','N_BC1_1 ',
     *     'N_BC2_1 ','N_BC3_1 ','N_DBC_1 ','N_BOC_1 ','N_BCS_1 ',
     *     'N_MXX_1 ','N_OCS_1 ')

        call set_diag_aod(n,k)
        if (diag_fc==2) then
          call set_diag_rf(n,k)
        else if (diag_fc==1) then
          select case (trname(n))
            case ('N_AKK_1')
              call set_diag_rf(n,k)
          end select
        endif

      end select
      end do

c - Tracer independent Diagnostic (stays here if 2D, moves to ijlt if 3D)
       k = k + 1
         ijts_AMPe(1)=k
         ia_ijts(k) = ia_src
         sname_ijts(k) = 'PM1'
         lname_ijts(k) = 'PM1 Mixing ratio'
         ijts_power(k) =  -9
         units_ijts(k) = unit_string(ijts_power(k),'kg kg-1')
         scale_ijts(k) = 10.**(-ijts_power(k))
         ijts_HasArea(k) = .false.

       k = k + 1
         ijts_AMPe(2)=k
         ia_ijts(k) = ia_src
         sname_ijts(k) = 'PM2.5'
         lname_ijts(k) = 'PM2.5 Mixing ratio'
         ijts_power(k) =  -9
         units_ijts(k) = unit_string(ijts_power(k),'kg kg-1')
         scale_ijts(k) = 10.**(-ijts_power(k))
         ijts_HasArea(k) = .false.

       k = k + 1
         ijts_AMPe(3)=k
         ia_ijts(k) = ia_src
         sname_ijts(k) = 'PM10'
         lname_ijts(k) = 'PM10 Mixing ratio'
         ijts_power(k) =  -9
         units_ijts(k) = unit_string(ijts_power(k),'kg kg-1')
         scale_ijts(k) = 10.**(-ijts_power(k))
         ijts_HasArea(k) = .false.

c      do L=1,1    !LTOP
c      do m=1,NBINS
c        k = k + 1
c         ijts_AMPpdf(l,m)=k
c         ia_ijts(k) = ia_src
c         write(lname_ijts(k),'(a15,i2.2,i2.2)') 'NUMB_PDF BIN L=',L,M
c         write(sname_ijts(k),'(a9,i2.2,i2.2)') 'N_PDF_BIN',L,M
c         ijts_power(k) = -2
c         units_ijts(k) = unit_string(ijts_power(k),'#')
c         scale_ijts(k) = 10.**(-ijts_power(k))
c      end do
c      end do
#endif  /* TRACERS_AMP */

c
c Append some denominator fields if necessary
c
      if(any(dname_ijts(1:k).eq.'clrsky')) then
        k = k + 1
        ijts_clrsky = k
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'CLEAR SKY FRACTION'
        sname_ijts(k) = 'clrsky'
        units_ijts(k) = '%'
        scale_ijts(k) = 100.
        ijts_HasArea(k) = .false.
      endif

      if(any(dname_ijts(1:k).eq.'ocnfr')) then
        k = k + 1
        ijts_pocean = k
        lname_ijts(k) = 'OCEAN FRACTION'
        units_ijts(k) = '%'
        sname_ijts(k) = 'ocnfr'
        ia_ijts(k) = ia_src     ! ia_ij(ij_pocean) is not initialized yet :(
        scale_ijts(k) = 100.
        ijts_HasArea(k) = .false.
      endif

      if(any(dname_ijts(1:k).eq.'sunlit_snow_freq')) then ! snow albedo weight
        k = k + 1
        ijts_sunlit_snow = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SUNLIT SNOW FREQUENCY'
        sname_ijts(k) = 'sunlit_snow_freq'
        units_ijts(k) = '%'
        scale_ijts(k) = 100.
        ijts_HasArea(k) = .false.
      endif

      ktaijs_out = k
      if (ktaijs_out .gt. ktaijs) then
        if (AM_I_ROOT())
     *       write (6,*)'ijt_defs: Increase ktaijs=',ktaijs
     *       ,' to at least ',ktaijs_out
        call stop_model('ktaijs too small',255)
      end if

c find indices of denominators
      call FindStrings(dname_ijts,sname_ijts,denom_ijts,ktaijs_out)
c      do k=1,ktaijs_out
c        if(len_trim(dname_ijts(k)).gt.0) then
c          do kk=ktaijs_out,1,-1
c            if(trim(sname_ijts(kk)).eq.trim(dname_ijts(k))) then
c              denom_ijts(k) = kk
c              exit
c            endif
c          enddo
c          if(denom_ijts(k).eq.0) then
c            if(am_i_root()) then
c              write(6,*)
c     &             'init_ijts_diag: denominator '//trim(dname_ijts(k))
c              write(6,*) 'not found for field '//trim(sname_ijts(k))
c            endif
c            call stop_model('init_ijts_diag: denominator not found',255)
c          endif
c        endif
c      enddo

#endif /* TRACERS_ON */

      return
      end subroutine init_ijts_diag

#ifdef TRACERS_ON
      subroutine FindStrings(StringsToFind,ListOfStrings,Indices,n)
!@sum FindStrings finds the positions of a list of strings in a 2nd list.
!     Needs optimization.
      use mdiag_com, only : sname_strlen
      implicit none
      integer :: n
      character(len=sname_strlen), dimension(n) ::
     &     StringsToFind,ListOfStrings
      integer, dimension(n) :: Indices
      integer :: k,kk
      logical :: found
      do k=1,n
        if(len_trim(StringsToFind(k)).gt.0) then
          found = .false.
          do kk=1,n
            if(trim(ListOfStrings(kk)).eq.trim(StringsToFind(k))) then
              Indices(k) = kk
              found = .true.
              exit
            endif
          enddo
          if(.not.found) then
            write(6,*) 'FindStrings: string '//
     &           trim(StringsToFind(k))//' not found'
            call stop_model('FindStrings: string not found',255)
          endif
        endif
      enddo
      end subroutine FindStrings
#endif

#ifdef TRACERS_ON
      subroutine set_diag_aod(n,k,n_subclasses)
!@sum set_diag_aod saves extinction, scattering and asymmetry parameter diags
!@auth Dorothy Koch, modified by Kostas Tsigaridis
      use OldTracer_mod, only: trname
      use mdiag_com, only : sname_strlen
      USE TRDIAG_COM, only: diag_rad,ijts_tau,ijts_sqex,ijts_sqsc
     &                     ,ijts_sqcb,ia_ijts,sname_ijts
     &                     ,lname_ijts,dname_ijts,ijts_power
     &                     ,units_ijts,scale_ijts,ijts_HasArea
     &                     ,ijts_tausub,ijts_sqexsub,ijts_sqscsub
     &                     ,ijts_sqcbsub,save_dry_aod
      USE DIAG_COM, only: ia_rad
      implicit none

      integer, intent(inout) :: k
      integer, intent(in) :: n
!@var n_subclasses optional argument for the number of sub classes of a given
!@+  tracer (>= 1)
      integer, optional, intent(in) :: n_subclasses
      character*50 :: unit_string
!@param sascs short name of all-sky/clear-sky selector
!@param lascs long name of all-sky/clear-sky selector
!@var s index of sascs and lascs
!@var kr index of solar bands
!@var skr value of kr as a string
!@var sn1 value of n1 as a string
      character(len=sname_strlen), parameter :: dname='clrsky'
      character(len=10), parameter, dimension(3) ::
     &  sascs=(/'    ','CS_ ','DRY_'/),
     &  lascs=(/'         ','clear sky','dry aeros'/)
      integer :: kr,s,n1,n_sub
      character(len=1) :: skr,sn1

      n_sub=1
      if (present(n_subclasses)) n_sub=n_subclasses

! aerosol optical depth and related diagnostics

      do s=1,size(sascs)
        if (trim(sascs(s)).eq.'DRY_' .and. save_dry_aod.eq.0) cycle
        IF (diag_rad /= 1) THEN
! aerosol optical depth for band6
          do n1=1,n_sub
            k = k + 1
            if (n_sub == 1) then
              ijts_tau(s,n) = k
              sn1=' '
            else
              ijts_tausub(s,n,n1) = k
              sn1=char(48+n1)
            end if
            ia_ijts(k)=ia_rad
            sname_ijts(k) = 'tau_'//trim(sascs(s))//trim(trname(n))//
     &           trim(sn1)
            lname_ijts(k) = trim(trname(n))//trim(sn1)//' '//
     &           trim(lascs(s))//' aerosol optical depth'
            if (trim(sascs(s))=='CS_') dname_ijts(k) = trim(dname)
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),' ')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do                ! n1
        ELSE
          DO kr=1,6
            write (skr,'(i1)') kr
! extinction aerosol optical depth in six solar bands
            do n1=1,n_sub
              k=k+1
              if (n_sub == 1) then
                ijts_sqex(s,kr,n)=k
                sn1=' '
              else
                ijts_sqexsub(s,kr,n,n1)=k
                sn1=char(48+n1)
              end if
              ia_ijts(k)=ia_rad
              sname_ijts(k)='ext_'//trim(sascs(s))//'band'//skr//'_'//
     &                    trim(trname(n))//trim(sn1)
              lname_ijts(k)=trim(trname(n))//trim(sn1)//' '//
     &             trim(lascs(s))//' SW extinction band '//skr
              if (trim(sascs(s))=='CS_') dname_ijts(k) = trim(dname)
              ijts_power(k) = -4
              units_ijts(k) = unit_string(ijts_power(k),' ')
              scale_ijts(k) = 10.**(-ijts_power(k))
              ijts_HasArea(k) = .false.
! scattering aerosol optical depth in six solar bands
              k=k+1
              if (n_sub == 1) then
                ijts_sqsc(s,kr,n)=k
                sn1=' '
              else
                ijts_sqscsub(s,kr,n,n1)=k
                sn1=char(48+n1)
              end if
              ia_ijts(k)=ia_rad
              sname_ijts(k)='sct_'//trim(sascs(s))//'band'//skr//'_'//
     &             trim(trname(n))//trim(sn1)
              lname_ijts(k)=trim(trname(n))//trim(sn1)//' '//
     &             trim(lascs(s))//' SW scattering band '//skr
              if (trim(sascs(s))=='CS_') dname_ijts(k) = trim(dname)
              ijts_power(k) = -4
              units_ijts(k) = unit_string(ijts_power(k),' ')
              scale_ijts(k) = 10.**(-ijts_power(k))
              ijts_HasArea(k) = .false.
! scattering asymmetry factor in six solar bands
              k=k+1
              if (n_sub == 1) then
                ijts_sqcb(s,kr,n)=k
                sn1=' '
              else
                ijts_sqcbsub(s,kr,n,n1)=k
                sn1=char(48+n1)
              end if
              ia_ijts(k)=ia_rad
              sname_ijts(k)='asf_'//trim(sascs(s))//'band'//skr//'_'//
     &             trim(trname(n))//(trim(sn1))
              lname_ijts(k)=trim(trname(n))//trim(sn1)//' '//
     &             trim(lascs(s))//' SW assymetry factor band '//skr
              if (trim(sascs(s))=='CS_') dname_ijts(k) = trim(dname)
              ijts_power(k) = -2
              units_ijts(k) = unit_string(ijts_power(k),' ')
              scale_ijts(k) = 10.**(-ijts_power(k))
              ijts_HasArea(k) = .false.
            end do              ! n1
          END DO                ! kr
        END IF
      enddo                     ! s

      return
      end subroutine set_diag_aod


      subroutine set_diag_rf(n,k,n_subclasses)
!@sum set_diag_rf saves shortwave and longwave forcing, for all-sky and
!@+               clear-sky, at surface and TOA
!@auth Kostas Tsigaridis
      use OldTracer_mod, only: trname
      use RunTimeControls_mod, only: tracers_amp, tracers_tomas
      use mdiag_com, only : sname_strlen
      USE TRDIAG_COM, only: ia_ijts,sname_ijts
     &                     ,lname_ijts,dname_ijts,ijts_power
     &                     ,units_ijts,scale_ijts,ijts_HasArea
     &                     ,ijts_fc,ijts_fcsub
      USE DIAG_COM, only: ia_rad_frc
      use RAD_COM, only: nradfrc,diag_fc
      implicit none

      integer, intent(inout) :: k
      integer, intent(in) :: n
!@var n_subclasses optional argument for the number of sub classes of a given
!@+  tracer (>= 1)
      integer, optional, intent(in) :: n_subclasses
      character*50 :: unit_string
!@param sascs short name of all-sky/clear-sky selector
!@param lascs long name of all-sky/clear-sky selector
!@param stoasrf short name of toa/surf selector
!@param ltoasrf long name of toa/surf selector
!@param sswlw short name of swf/lwf selector
!@param lswlw long name of swf/lwf selector
!@var s index of sascs and lascs
!@var l index of stoasrf and ltoasrf
!@var f index of sswlw and lswlw
!@var i combined index of s,l,f
!@var kr index of solar bands
!@var skr value of kr as a string
!@var sn1 value of n1 as a string
      character(len=sname_strlen), parameter :: dname='clrsky'
      character(len=10), parameter, dimension(2) ::
     &  sascs=(/'   ','CS_'/),lascs=(/'         ','clear sky'/),
     &  stoasrf=(/'     ','surf_'/),ltoasrf=(/'TOA    ','surface'/),
     &  sswlw=(/'swf_','lwf_'/),lswlw=(/'shortwave','longwave '/)
      integer :: kr,s,l,f,i,n1,n_sub
      character(len=1) :: skr,sn1
      character(len=20) :: spcname ! following MAX_LEN_NAME=20 for trname

      n_sub=1
      if (present(n_subclasses)) n_sub=n_subclasses

! radiative forcing and related diagnostics

      if (diag_fc==2) then
        spcname=trim(trname(n))
      else if (diag_fc==1) then
        if (tracers_amp) then
          spcname='AMP'
        elseif (tracers_tomas) then
          spcname='TOMAS'
        else
          spcname='OMA'
        endif
      endif

      if (nradfrc>0) then
        do s=1,size(sascs)
        do l=1,size(stoasrf)
        do f=1,size(sswlw)
          i=(s-1)*size(stoasrf)*size(sswlw)+(l-1)*size(sswlw)+f
          do n1=1,n_sub
            k = k + 1
            if (n_sub == 1) then
              ijts_fc(i,n) = k
              sn1=' '
            else
              ijts_fcsub(i,n,n1) = k
              sn1=char(48+n1)
            end if
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(spcname)//trim(sn1)//' '//
     &                      trim(lswlw(f))
            sname_ijts(k) = trim(sswlw(f))
            if (trim(sascs(s))=='CS_') then
              lname_ijts(k) = trim(lname_ijts(k))//' '//trim(lascs(s))
              sname_ijts(k) = trim(sname_ijts(k))//trim(sascs(s))
              dname_ijts(k) = trim(dname)
            endif
            if (trim(stoasrf(l))=='surf_') then
              lname_ijts(k) = trim(lname_ijts(k))//' '//trim(ltoasrf(l))
              sname_ijts(k) = trim(sname_ijts(k))//trim(stoasrf(l))
            endif
            lname_ijts(k) = trim(lname_ijts(k))//' radiative forcing'
            sname_ijts(k) = trim(sname_ijts(k))//trim(spcname)//
     &           trim(sn1)
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W m-2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do                ! n1
        enddo ! f
        enddo ! l
        enddo ! s
      endif

      return
      end subroutine set_diag_rf
#endif  /* TRACERS_ON */

      subroutine init_ijlts_diag
!@sum init_ijlts_diag Initialise lat/lon/height tracer diags
!@auth Gavin Schmidt
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use OldTracer_mod, only: trname,max_len_name
      USE TRACER_COM, only: ntm
      use rad_com, only: nraero_aod,ntrix_aod
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif /* TRACERS_ON */
#ifdef ACCMIP_LIKE_DIAGS
      USE MODEL_COM, only: dtsrc
#endif
      USE DIAG_COM
#ifdef SOA_DIAGS
      use tracers_soa, only: issoa
#endif  /* SOA_DIAGS */
      implicit none
      integer k,n,i
      character*50 :: unit_string
      character(len=max_len_name) :: trname_curr
      character*1 :: clay_num
      integer :: iclay
      integer :: ktaijlt_out

#ifdef TRACERS_ON
      ir_ijlt = ir_log2  ! default
      ia_ijlt = ia_src   ! default
      denom_ijlt(:) = 0
#ifdef TRACERS_AMP
      ijlt_AMPm(:,:)=0
#endif

      k=0
C**** use this routine to set 3D tracer-related diagnostics.

C**** some tracer specific 3D arrays
      if (diag_aod_3d>0 .and. diag_aod_3d<5) then ! valid values are 1-4
        allocate(ijlt_3Dtau(nraero_aod))    ; ijlt_3Dtau = 0
        allocate(ijlt_3DtauCS(nraero_aod))  ; ijlt_3DtauCS = 0
        allocate(ijlt_3DtauDRY(nraero_aod))  ; ijlt_3DtauDRY = 0
        allocate(ijlt_3Daaod(nraero_aod))   ; ijlt_3Daaod = 0
        allocate(ijlt_3DaaodCS(nraero_aod)) ; ijlt_3DaaodCS = 0
        allocate(ijlt_3DaaodDRY(nraero_aod)) ; ijlt_3DaaodDRY = 0

        iclay=0
        do n=1,nraero_aod
          trname_curr=trim(trname(ntrix_aod(n)))
          select case (trname(ntrix_aod(n)))
            case ('Clay','ClayIlli','ClayKaol','ClaySmec','ClayCalc'
     &           ,'ClayQuar','ClayFeld','ClayHema','ClayGyps'
     &           ,'ClayIlHe','ClayKaHe','ClaySmHe','ClayCaHe'
     &           ,'ClayQuHe','ClayFeHe','ClayGyHe')
              iclay=iclay+1
              write(clay_num, '(i1)') iclay
              trname_curr=trim(trname_curr)//trim(clay_num)
          end select

          if (diag_aod_3d==1 .or. diag_aod_3d==3) then
            k = k + 1
            ijlt_3Dtau(n)=k
            ia_ijlt(k) = ia_rad
            lname_ijlt(k) = trim(trname_curr)//' tau'
            sname_ijlt(k) = 'tau_3D_'//trim(trname_curr)
            ijlt_power(k) = -2
            units_ijlt(k) = unit_string(ijlt_power(k),' ')
            scale_ijlt(k) = 10.**(-ijlt_power(k))
            k = k + 1
            ijlt_3Daaod(n)=k
            ia_ijlt(k) = ia_rad
            lname_ijlt(k) = trim(trname_curr)//' aaod'
            sname_ijlt(k) = 'aaod_3D_'//trim(trname_curr)
            ijlt_power(k) = -2
            units_ijlt(k) = unit_string(ijlt_power(k),' ')
            scale_ijlt(k) = 10.**(-ijlt_power(k))
          endif ! diag_aod_3d = 1 or 3

          if (diag_aod_3d==2 .or. diag_aod_3d==3) then
            k = k + 1
            ijlt_3DtauCS(n)=k
            ia_ijlt(k) = ia_rad
            lname_ijlt(k) = trim(trname_curr)//' CS tau'
            sname_ijlt(k) = 'tau_3D_CS_'//trim(trname_curr)
            dname_ijlt(k) = 'clrsky2d'
            ijlt_power(k) = -2
            units_ijlt(k) = unit_string(ijlt_power(k),' ')
            scale_ijlt(k) = 10.**(-ijlt_power(k))
            k = k + 1
            ijlt_3DaaodCS(n)=k
            ia_ijlt(k) = ia_rad
            lname_ijlt(k) = trim(trname_curr)//' CS aaod'
            sname_ijlt(k) = 'aaod_3D_CS_'//trim(trname_curr)
            dname_ijlt(k) = 'clrsky2d'
            ijlt_power(k) = -2
            units_ijlt(k) = unit_string(ijlt_power(k),' ')
            scale_ijlt(k) = 10.**(-ijlt_power(k))
          endif ! diag_aod_3d = 2 or 3

          if (diag_aod_3d==4 .or. diag_aod_3d==3) then
          if (save_dry_aod>0) then
            k = k + 1
            ijlt_3DtauDRY(n)=k
            ia_ijlt(k) = ia_rad
            lname_ijlt(k) = trim(trname_curr)//' DRY tau'
            sname_ijlt(k) = 'tau_3D_DRY_'//trim(trname_curr)
            ijlt_power(k) = -2
            units_ijlt(k) = unit_string(ijlt_power(k),' ')
            scale_ijlt(k) = 10.**(-ijlt_power(k))
            k = k + 1
            ijlt_3DaaodDRY(n)=k
            ia_ijlt(k) = ia_rad
            lname_ijlt(k) = trim(trname_curr)//' DRY aaod'
            sname_ijlt(k) = 'aaod_3D_DRY_'//trim(trname_curr)
            ijlt_power(k) = -2
            units_ijlt(k) = unit_string(ijlt_power(k),' ')
            scale_ijlt(k) = 10.**(-ijlt_power(k))
          endif ! save_dry_aod>0
          endif ! diag_aod_3d = 4 or 3

        enddo ! nraero_aod
      else if (diag_aod_3d<0 .and. diag_aod_3d>-5) then ! if negative, save total
        allocate(ijlt_3Dtau(1))    ; ijlt_3Dtau = 0
        allocate(ijlt_3DtauCS(1))  ; ijlt_3DtauCS = 0
        allocate(ijlt_3DtauDRY(1))  ; ijlt_3DtauDRY = 0
        allocate(ijlt_3Daaod(1))   ; ijlt_3Daaod = 0
        allocate(ijlt_3DaaodCS(1)) ; ijlt_3DaaodCS = 0
        allocate(ijlt_3DaaodDRY(1)) ; ijlt_3DaaodDRY = 0

        if (diag_aod_3d==-1 .or. diag_aod_3d==-3) then
          k = k + 1
          ijlt_3Dtau(1)=k
          ia_ijlt(k) = ia_rad
          lname_ijlt(k) = 'tau'
          sname_ijlt(k) = 'tau_3D'
          ijlt_power(k) = -2
          units_ijlt(k) = unit_string(ijlt_power(k),' ')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
          k = k + 1
          ijlt_3Daaod(1)=k
          ia_ijlt(k) = ia_rad
          lname_ijlt(k) = 'aaod'
          sname_ijlt(k) = 'aaod_3D'
          ijlt_power(k) = -2
          units_ijlt(k) = unit_string(ijlt_power(k),' ')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        endif ! diag_aod_3d = -1 or -3

        if (diag_aod_3d==-2 .or. diag_aod_3d==-3) then
          k = k + 1
          ijlt_3DtauCS(1)=k
          ia_ijlt(k) = ia_rad
          lname_ijlt(k) = 'CS tau'
          sname_ijlt(k) = 'tau_3D_CS'
          dname_ijlt(k) = 'clrsky2d'
          ijlt_power(k) = -2
          units_ijlt(k) = unit_string(ijlt_power(k),' ')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
          k = k + 1
          ijlt_3DaaodCS(1)=k
          ia_ijlt(k) = ia_rad
          lname_ijlt(k) = 'CS aaod'
          sname_ijlt(k) = 'aaod_3D_CS'
          dname_ijlt(k) = 'clrsky2d'
          ijlt_power(k) = -2
          units_ijlt(k) = unit_string(ijlt_power(k),' ')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        endif ! diag_aod_3d = -2 or -3

        if (diag_aod_3d==-4 .or. diag_aod_3d==-3) then
        if (save_dry_aod>0) then
          k = k + 1
          ijlt_3DtauDRY(1)=k
          ia_ijlt(k) = ia_rad
          lname_ijlt(k) = 'DRY tau'
          sname_ijlt(k) = 'tau_3D_DRY'
          ijlt_power(k) = -2
          units_ijlt(k) = unit_string(ijlt_power(k),' ')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
          k = k + 1
          ijlt_3DaaodDRY(1)=k
          ia_ijlt(k) = ia_rad
          lname_ijlt(k) = 'DRY aaod'
          sname_ijlt(k) = 'aaod_3D_DRY'
          ijlt_power(k) = -2
          units_ijlt(k) = unit_string(ijlt_power(k),' ')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        endif ! save_dry_aod>0
        endif ! diag_aod_3d = -4 or -3
      endif ! 0<diag_aod_3d<5 or 0>diag_aod_3d>-5

      do n=1,NTM
        select case(trname(n))

        CASE('Clay','Silt1','Silt2','Silt3','Silt4','Silt5', 'isopp1a'
     $       ,'isopp2a','apinp1a','apinp2a','OCB','OCII','OCIA','BCB'
     $       ,'BCII' ,'BCIA', 'SO4','MSA','NO3p','NH4','seasalt1'
     $       ,'seasalt2','SO4_d1','SO4_d2','SO4_d3','N_d1','N_d2'
     $       ,'N_d3')

#ifdef SAVE_AEROSOL_3DMASS_FOR_NINT
        k = k + 1
        ijlt_3Dmass(n)=k
        ia_ijlt(k) = ia_src     ! just guessing
        lname_ijlt(k) = trim(trname(n))//' Mass'
        sname_ijlt(k) = 'Mass_3D_'//trname(n)
        ijlt_power(k) = -5
        units_ijlt(k) = unit_string(ijlt_power(k),' kg m-2')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
#endif /* define io parameters for 3Dmass diagnostic (Ron) */

#ifdef TRACERS_AMP
c- 3D diagnostic per mode
      CASE('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 ',
     *     'N_DS2_1 ','N_SSA_1 ','N_SSC_1 ','N_OCC_1 ','N_BC1_1 ',
     *     'N_BC2_1 ','N_BC3_1 ','N_DBC_1 ','N_BOC_1 ','N_BCS_1 ',
     *     'N_MXX_1 ','N_OCS_1 ')
        k = k + 1
         ijlt_AMPm(1,n)=k
         lname_ijlt(k) = TRIM(trname(n))//' DIAM'
         sname_ijlt(k) = 'DIAM_'//TRIM(trname(n))
         ijlt_power(k) = -2.
         units_ijlt(k) = unit_string(ijlt_power(k),'m')
         scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
         ijlt_AMPm(3,n)=k
         lname_ijlt(k) = TRIM(trname(n))//' DIAM_DRY'
         sname_ijlt(k) = 'DIAM_DRY_'//TRIM(trname(n))
         ijlt_power(k) = -2.
         units_ijlt(k) = unit_string(ijlt_power(k),'m')
         scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
         ijlt_AMPm(2,n)=k
         lname_ijlt(k) = TRIM(trname(n))//' ACTI'
         sname_ijlt(k) = 'ACTI3D_'//TRIM(trname(n))
         ijlt_power(k) = -2.
         units_ijlt(k) = unit_string(ijlt_power(k),'#')
         scale_ijlt(k) = 10.**(-ijlt_power(k))
#endif
      end select
      end do

C**** 3D tracer-related arrays but not attached to any one tracer

#ifdef TRACERS_SPECIAL_Shindell
#ifdef ACCMIP_LIKE_DIAGS
      k = k + 1
        ijlt_OHvmr=k
        lname_ijlt(k) = 'OH mixing ratio'
        sname_ijlt(k) = 'OH_vmr'
        ijlt_power(k) = -10
        units_ijlt(k) = unit_string(ijlt_power(k),'V/V air')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
#endif
      k = k + 1
        ijlt_OHconc=k
        lname_ijlt(k) = 'OH concentration'
        sname_ijlt(k) = 'OH_conc'
        ijlt_power(k) = 5
        units_ijlt(k) = unit_string(ijlt_power(k),'molecules cm-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_NO3=k
        lname_ijlt(k) = 'NO3 concentration'
        sname_ijlt(k) = 'NO3_conc'
        ijlt_power(k) = 5
        units_ijlt(k) = unit_string(ijlt_power(k),'molecules cm-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_HO2=k
        lname_ijlt(k) = 'HO2 concentration'
        sname_ijlt(k) = 'HO2_conc'
        ijlt_power(k) = 7
        units_ijlt(k) = unit_string(ijlt_power(k),'molecules cm-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_JO1D=k
        lname_ijlt(k) = 'Ox to O1D photolysis rate'
        sname_ijlt(k) = 'JO1D'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_JNO2=k
        lname_ijlt(k) = 'NO2 photolysis rate'
        sname_ijlt(k) = 'JNO2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_JH2O2=k
        lname_ijlt(k) = 'H2O2 photolysis rate'
        sname_ijlt(k) = 'JH2O2'
        ijlt_power(k) = 2
        units_ijlt(k) = unit_string(ijlt_power(k),'s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_O3ppbv=k
        lname_ijlt(k) = 'O3 not Ox volume mixing ratio'
        sname_ijlt(k) = 'O3_vmr'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ppbv')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_O3cmatm=k
        lname_ijlt(k) = 'O3 not Ox in cm-atm units'
        sname_ijlt(k) = 'O3_cm_atm'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'cm-atm')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
#ifdef ACCMIP_LIKE_DIAGS
      k = k + 1
        ijlt_COp=k
        lname_ijlt(k) = 'CO production rate'
        sname_ijlt(k) = 'COprod'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_COd=k
        lname_ijlt(k) = 'CO destruction rate'
        sname_ijlt(k) = 'COdest'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_Oxp=k
        lname_ijlt(k) = 'Ox production rate'
        sname_ijlt(k) = 'Oxprod'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_Oxd=k
        lname_ijlt(k) = 'Ox destruction rate'
        sname_ijlt(k) = 'Oxdest'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_CH4d=k
        lname_ijlt(k) = 'CH4 destruction rate'
        sname_ijlt(k) = 'CH4dest'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_OxpHO2=k
        lname_ijlt(k) = 'Ox prod rate via HO2+NO'
        sname_ijlt(k) = 'OxpHO2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxpCH3O2=k
        lname_ijlt(k) = 'Ox prod rate via CH3O2+NO'
        sname_ijlt(k) = 'OxpCH3O2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxpRO2=k
        lname_ijlt(k) = 'Ox prod rate via RO2+NO'
        sname_ijlt(k) = 'OxpRO2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxlOH=k
        lname_ijlt(k) = 'Ox loss rate via OH'
        sname_ijlt(k) = 'OxlOH'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxlHO2=k
        lname_ijlt(k) = 'Ox loss rate via HO2'
        sname_ijlt(k) = 'OxlHO2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxlALK=k
        lname_ijlt(k) = 'Ox loss rate via Alkenes'
        sname_ijlt(k) = 'OxlALK'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_pO1D=k
        lname_ijlt(k) = 'O1D production from ozone'
        sname_ijlt(k) = 'pO1d'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_pOH=k
        lname_ijlt(k) = 'OH production from O1D+H2O'
        sname_ijlt(k) = 'pOH'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_NOxLgt=k
        lname_ijlt(k) = 'NOx production from Lightning'
        sname_ijlt(k) = 'NOx_Lightn'
        ijlt_power(k) = -15
        units_ijlt(k) = unit_string(ijlt_power(k),'kg(N) m-2 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_NOvmr=k
        lname_ijlt(k) = 'NO mixing ratio'
        sname_ijlt(k) = 'NO_vmr'
        ijlt_power(k) = -10 ! to match NOx
        units_ijlt(k) = unit_string(ijlt_power(k),'V/V air')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_NO2vmr=k
        lname_ijlt(k) = 'NO2 mixing ratio'
        sname_ijlt(k) = 'NO2_vmr'
        ijlt_power(k) = -10 ! to match NOx
        units_ijlt(k) = unit_string(ijlt_power(k),'V/V air')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_prodSO4aq=k
        lname_ijlt(k) = 'SO4 aqueous chem source 3D'
        sname_ijlt(k) = 'SO4aqSrc3D'
        ijlt_power(k) = -15 ! to match ijts 2D
        units_ijlt(k) = unit_string(ijlt_power(k),'kg m-2 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_prodSO4gs=k
        lname_ijlt(k) = 'SO4 gas phase source 3D'
        sname_ijlt(k) = 'SO4gasSrc3D'
        ijlt_power(k) = -15 ! to match ijts 2D
        units_ijlt(k) = unit_string(ijlt_power(k),'kg m-2 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
#endif /* ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */

#ifdef TRACERS_NITRATE
      k = k + 1
        ijlt_aH2O=k
        lname_ijlt(k) = 'aerosol H2O'
        sname_ijlt(k) = 'aerosol_H2O'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_apH=k
        lname_ijlt(k) = 'aerosol pH'
        sname_ijlt(k) = 'aerosol_pH'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
#endif  /* TRACERS_NITRATE */

#ifdef SOA_DIAGS
      k = k + 1
        ijlt_soa_changeL_isoprene=k
        lname_ijlt(k) = 'changeL of isoprene'
        sname_ijlt(k) = 'SOA_changeL_isoprene'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_changeL_terpenes=k
        lname_ijlt(k) = 'changeL of terpenes'
        sname_ijlt(k) = 'SOA_changeL_terpenes'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_voc2nox=k
        lname_ijlt(k) = 'VOC/NOx ratio'
        sname_ijlt(k) = 'SOA_voc2nox'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ppbC/ppb')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_pcp=k
        lname_ijlt(k) = 'Total non-volatile SOA-absorbing mass'
        sname_ijlt(k) = 'SOA_pcp'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_aerotot=k
        lname_ijlt(k) = 'PCP plus SOA'
        sname_ijlt(k) = 'SOA_aerotot'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3 per MW')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_aerotot_gas=k
        lname_ijlt(k) = 'Gas-phase semivolatile potential SOA'
        sname_ijlt(k) = 'SOA_aerotot_gas'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3 per MW')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_xmf_isop=k
        lname_ijlt(k) = 'Molar fraction of isoprene SOA'
        sname_ijlt(k) = 'SOA_xmf_isop'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'fraction')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_xmf_apin=k
        lname_ijlt(k) = 'Molar fraction of a-pinene SOA'
        sname_ijlt(k) = 'SOA_xmf_apin'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'fraction')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_zcoef_isop=k
        lname_ijlt(k) = 'Activity coefficient for isoprene SOA'
        sname_ijlt(k) = 'SOA_zcoef_isop'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'dimensionless')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_zcoef_apin=k
        lname_ijlt(k) = 'Activity coefficient for a-pinene SOA'
        sname_ijlt(k) = 'SOA_zcoef_apin'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'dimensionless')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_meanmw=k
        lname_ijlt(k) = 'Mean organic aerosol molecular weight'
        sname_ijlt(k) = 'SOA_meanmw'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'g mol-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_iternum=k
        lname_ijlt(k) = 'Total iterations for SOA calculations'
        sname_ijlt(k) = 'SOA_iternum'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'count')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_m0=k
        lname_ijlt(k) = 'Final M0 value'
        sname_ijlt(k) = 'SOA_M0'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      do i=1,nsoa
        k = k + 1
          ijlt_soa_y0_ug_g(i)=k
          lname_ijlt(k) = 'y0_ug of '//trim(trname(issoa(i)-1))
          sname_ijlt(k) = 'SOA_y0_ug_'//trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_y0_ug_a(i)=k
          lname_ijlt(k) = 'y0_ug of '//trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_y0_ug_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_y_ug_g(i)=k
          lname_ijlt(k) = 'y_ug of '//trim(trname(issoa(i)-1))
          sname_ijlt(k) = 'SOA_y_ug_'//trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_y_ug_a(i)=k
          lname_ijlt(k) = 'y_ug of '//trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_y_ug_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_changeL_g_before(i)=k
          lname_ijlt(k) = 'changeL of '//trim(trname(issoa(i)-1))//
     &                    ' before SOA'
          sname_ijlt(k) = 'SOA_changeL_before_'//
     &                    trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_changeL_a_before(i)=k
          lname_ijlt(k) = 'changeL of '//trim(trname(issoa(i)))//
     &                    ' before SOA'
          sname_ijlt(k) = 'SOA_changeL_before_'//
     &                    trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_changeL_g_after(i)=k
          lname_ijlt(k) = 'changeL of '//trim(trname(issoa(i)-1))//
     &                    ' after SOA'
          sname_ijlt(k) = 'SOA_changeL_after_'//
     &                    trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_changeL_a_after(i)=k
          lname_ijlt(k) = 'changeL of '//trim(trname(issoa(i)))//
     &                    ' after SOA'
          sname_ijlt(k) = 'SOA_changeL_after_'//
     &                    trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_apartmass(i)=k
          lname_ijlt(k) = 'Effective apartmass of '//
     &                     trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_apartmass_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'dimensionless')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_kpart(i)=k
          lname_ijlt(k) = 'Partitioning coefficient of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_kpart_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'m3 ug-1')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_kp(i)=k
          lname_ijlt(k) = 'Final partitioning coefficient of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_kp_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'m3 ug-1')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_soamass(i)=k
          lname_ijlt(k) = 'Potential SOA mass of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_soamass_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_partfact(i)=k
          lname_ijlt(k) = 'Final partfact value of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_partfact_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'dimensionless')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_evap(i)=k
          lname_ijlt(k) = 'Evaporation of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_evap_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_cond(i)=k
          lname_ijlt(k) = 'Condensation of pre-existing '//
     &                    trim(trname(issoa(i)-1))
          sname_ijlt(k) = 'SOA_cond_'//trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_chem(i)=k
          lname_ijlt(k) = 'Condensation of same-step produced '//
     &                    trim(trname(issoa(i)-1))
          sname_ijlt(k) = 'SOA_chem_'//trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
      enddo
#endif  /* SOA_DIAGS */

#ifdef TRACERS_TOMAS

      k = k + 1
        ijlt_ccn_01=k
        lname_ijlt(k) = 'CCN 0.1% '
        sname_ijlt(k) = 'CCN_01_SS'
        ijlt_power(k) = 0 ! to match ijts 2D
        units_ijlt(k) = unit_string(ijlt_power(k),'cm-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_ccn_02=k
        lname_ijlt(k) = 'CCN 0.2%'
        sname_ijlt(k) = 'CCN_02_SS'
        ijlt_power(k) = 0 ! to match ijts 2D
        units_ijlt(k) = unit_string(ijlt_power(k),'cm-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_ccn_03=k
        lname_ijlt(k) = 'CCN 0.3%'
        sname_ijlt(k) = 'CCN_03_SS'
        ijlt_power(k) = 0 ! to match ijts 2D
        units_ijlt(k) = unit_string(ijlt_power(k),'cm-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))

#endif /* TRACERS_TOMAS */

c
c Append some denominator fields if necessary
c

c nothing is using this as a denominator yet, but some fields _could_.
c      if(any(dname_ijlt(1:k).eq.'airmass')) then
        k = k + 1
        ijlt_airmass = k
        lname_ijlt(k) = 'Air Mass'
        sname_ijlt(k) = 'airmass'
        units_ijlt(k) = 'kg/m2/layer'
        scale_ijlt(k) = 1.
c      endif

      if(any(dname_ijlt(1:k).eq.'clrsky2d')) then
        k = k + 1
        ijlt_clrsky2d = k
        ia_ijlt(k) = ia_rad
        lname_ijlt(k) = 'CLEAR SKY FRACTION'
        sname_ijlt(k) = 'clrsky2d'
        units_ijlt(k) = '%'
        scale_ijlt(k) = 100.
      endif

      ktaijlt_out = k
      if (k .gt. ktaijl) then
       if (AM_I_ROOT())
     *       write (6,*)'ijlt_defs: Increase ktaijl=',ktaijl
     *       ,' to at least ',k
        call stop_model('ktaijl too small',255)
      end if

c find indices of denominators
      call FindStrings(dname_ijlt,sname_ijlt,denom_ijlt,ktaijlt_out)

#endif /* TRACERS_ON */

      return
      end subroutine init_ijlts_diag

      function get_atmco2()
      USE RADPAR, only : xnow
      use runtimecontrols_mod, only: constco2
      use dictionary_mod, only: sync_param
      use model_com, only: modelEclock
      use ghg_mod, only : CO2_trend
      implicit none
      real*8 :: get_atmco2
      real*8, save :: atmco2=-1.
      integer :: year, dayOfYear
      real*8 :: tnow

      if (constco2) then
        if (atmco2<0.) then    ! uninitialized
          atmco2=280.
          call sync_param("atmCO2", atmco2)
        endif
        if ( atmco2 > 0.d0 ) then
          get_atmco2=atmco2
        else
          call modelEclock%get(year=year, dayOfYear=dayOfYear)
          tnow = year + (dayOfYear-0.999d0)/366.d0
          call CO2_trend(get_atmco2, tnow)
        endif
      else
        get_atmco2=xnow(1)
      endif
      return
      end function get_atmco2

      SUBROUTINE tracer_IC
!@sum tracer_IC initializes tracers when they are first switched on
!@vers 2013/03/27
!@auth Jean Lerner
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT,readt_parallel,
     &     readt8_column, skip_parallel
      USE Dictionary_mod, only : get_param, is_set_param
#ifdef TRACERS_ON
      USE FLUXES, only : atmocn,atmice,atmglas,atmlnd,atmsrf,asflx
#ifdef GLINT2
      USE FLUXES, only : atmglas_hp
#endif
      USE CONSTANT, only: mair,rhow,grav,tf,avog,rgas
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE resolution,ONLY : Im,Jm,Lm,Ls1=>ls1_nominal
      USE ATM_COM, only : q,qcl,qci
      use model_com, only: modelEclock
      USE MODEL_COM, only: itime,dtsrc,itimeI
      USE ATM_COM, only: pmidl00
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds,write_parallel
      USE SOMTQ_COM, only : qmom,mz,mzz
      use OldTracer_mod, only: trname, itime_tr0, vol2mass, tr_mm
      use OldTracer_mod, only: trsi0, needtrs
      use OldTracer_mod, only: set_itime_tr0
      USE TRACER_COM, only: NTM, trm, trmom, rnsrc, tracers
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only:
     *     n_ASO4,n_AOCOB,n_ASO4,n_ANUM,xk,nbins
#endif
#ifdef TRACERS_WATER
      use OldTracer_mod, only: trw0, tr_wd_type, nWATER
      USE TRACER_COM,only: trwm,n_HDO,n_H2O18, n_OCII
      USE LANDICE, only : ace1li,ace2li
#ifndef TRACERS_ATM_ONLY
      USE LANDICE_COM, only : trsnowli,trlndi
      USE LAKES_COM, only : trlake
      USE GHY_COM, only : tr_w_ij,tr_wsn_ij
#endif
      USE LANDICE_COM, only : snowli
      USE SEAICE_COM, only : si_atm,si_ocn
      USE LAKES_COM, only : mwl,mldlk,flake
      USE GHY_COM, only : w_ij,wsn_ij,nsn_ij,fr_snow_ij,fearth
      USE FLUXES, only : flice,focean
#endif
      USE GEOM, only: axyp,byaxyp,lat2d_dg,lonlat_to_ij
      USE ATM_COM, only: MA,byMA  ! Air mass of each box (kg m-2)
      USE PBLCOM, only: npbl
#ifdef TRACERS_SPECIAL_Lerner
      USE LINOZ_CHEM_COM, only: tlt0m,tltzm, tltzzm
      USE PRATHER_CHEM_COM, only: nstrtc
#endif
      USE FILEMANAGER, only: openunit,closeunit,nameunit
#ifdef TRACERS_SPECIAL_Shindell
      USE RAD_COM, only : chem_tracer_save,rad_to_file,ghg_yr
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      USE RAD_COM, only: stratO3_tracer_save
#endif
#ifdef TRACERS_dCO
      use tracers_dCO, only: dalke_IC_fact
      use tracers_dCO, only: dPAR_IC_fact
      use tracers_dCO, only: dPAN_IC_fact
      use tracers_dCO, only: dMeOOH_IC_fact
      use tracers_dCO, only: dHCHO_IC_fact
      use tracers_dCO, only: dCO_IC_fact
#endif  /* TRACERS_dCO */
      USE TRCHEM_Shindell_COM,only:O3MULT,ch4icx,
     &  OxIC,COIC,byO3MULT,PI_run,fix_CH4_chemistry,
     &  PIratio_N,PIratio_CO_T,PIratio_CO_S,PIratio_other
     &  ,use_rad_n2o,use_rad_cfc,use_rad_ch4
     &  ,ClOxalt,BrOxalt,ClONO2alt,HClalt,N2OICX,CFCIC
     &  ,PIratio_N2O,PIratio_CFC
#ifdef INTERACTIVE_WETLANDS_CH4
      USE TRACER_SOURCES, only:first_mod,first_ncep,avg_model,avg_ncep,
     & PRS_ch4,sum_ncep
#endif
#ifdef SHINDELL_STRAT_EXTRA
      USE TRACER_SOURCES, only:GLTic
#endif
#endif /* TRACERS_SPECIAL_Shindell */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      use OldTracer_mod, only: om2oc
      USE AEROSOL_SOURCES, only: DMSinput
#ifndef TRACERS_AEROSOLS_SOA
      USE AEROSOL_SOURCES, only: OCT_src
#endif  /* TRACERS_AEROSOLS_SOA */
      USE AEROSOL_SOURCES, only: SO2_src_3D,iso2volcano
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      USE trdust_mod,ONLY : hbaij,ricntd
      use trdust_drv, only: tracer_ic_soildust
#endif
      use OldTracer_mod, only: trli0
#ifdef TRACERS_AMP
      use TRACER_COM, only:n_M_OCC_OC
#else
      use TRACER_COM, only:n_OCII
#endif
#endif /* TRACERS_ON */
      use oldtracer_mod, only: src_dist_base, src_dist_index
      use tracer_com, only: xyztr

      IMPLICIT NONE
      real*8,parameter :: d18oT_slope=0.45,tracerT0=25
      INTEGER i,n,l,j,iu_data,ipbl,it,lr,m,ls,lt,ipatch
      CHARACTER*80 title
      CHARACTER*300 out_line
      REAL*8 CFC11ic,conv
      REAL*8 :: trinit =1., tmominit=0.
      real*8 tracerTs

#ifdef TRACERS_ON
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     *                                CO2ic
      REAL*8, DIMENSION(GRID%I_STRT:GRID%I_STOP,
     &                  GRID%J_STRT:GRID%J_STOP,lm) ::
     *                                ic14CO2
      REAL*4, DIMENSION(jm,lm)    ::  N2Oic   !each proc. reads global array
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     *                                                      CH4ic
#ifdef TRACERS_SPECIAL_Lerner
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: icCFC
      REAL*8 stratm,xlat,pdn,pup
#endif
      real*8 :: get_atmco2

!@param bymair 1/molecular wt. of air = 1/mair
!@param byjm 1./JM
      REAL*8, PARAMETER :: bymair = 1.d0/mair, byjm =1.d0/JM
#ifdef TRACERS_SPECIAL_Shindell
      character*4 ghg_name
      character*80 ghg_file
      real*8, dimension(LM,GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                     GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: ghg_in
!@var imonth dummy index for choosing the right month
!@var ICfactor varying factor for altering initial conditions
!@var dICfactor varying factor for altering initial conditions of dCO tracers
      INTEGER imonth, J2
      REAL*8 ICfactor,dICfactor
!@var PRES local nominal pressure for vertical interpolations
      REAL*8, DIMENSION(LM) :: PRES
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)
      include 'netcdf.inc'
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)
      integer start(3),count(3),status,ncidu,id1
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS)
      INTEGER ii,jj,ir,mm,iuc,mmm,ll,iudms
      INTEGER iuc2,lmax
#endif
#endif
#if defined (TRACERS_AEROSOLS_Koch) || defined (TRACERS_AMP) ||\
    defined (TRACERS_TOMAS)
#ifdef CUBED_SPHERE /* 1-deg volc. emiss */
      real*8 :: volc_lons(360),volc_lats(180),
     &     volc_pup(360,180),volc_emiss(360,180)
#else /* volc. emiss on model grid, 1 extra lat at SP */
      real*8 :: volc_lons(Im)
      real*8, allocatable, dimension(:)   ::volc_lats
      real*8, allocatable, dimension(:,:) ::volc_pup, volc_emiss
#endif
      real*8 :: x1d(lm),amref(lm),pednref(lm+1),amsum
      real*8, allocatable, dimension(:,:) :: psref
      integer :: iu_ps,file_id,vid,ilon,jlat,volc_ij(2)
#endif
#ifdef TRACERS_TOMAS
      integer k
#endif
      INTEGER J_0, J_1, I_0, I_1
      INTEGER J_0H, J_1H
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      integer :: initial_GHG_setup
      integer :: lat_val
#endif /* TRACERS_ON */
      character(len=:), allocatable :: name

#ifdef TRACERS_ON
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0,       J_STOP=J_1,
     *               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     *               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP


#ifdef TRACERS_SPECIAL_Shindell
      PRES(1:LM)=PMIDL00(1:LM)
#endif
      do n=1,ntm
      if (itime.eq.itime_tr0(n)) then

C**** set some defaults for air mass tracers
      trm(:,J_0:J_1,:,n) = 0.
      trmom(:,:,J_0:J_1,:,n) = 0.

#ifdef TRACERS_WATER
C**** set some defaults for water tracers
      trwm(:,J_0:J_1,:,n)=0. ! cloud liquid water
#ifndef TRACERS_ATM_ONLY
      trlake(n,:,:,J_0:J_1)=0.
      si_atm%trsi(n,:,:,J_0:J_1)=0.
#endif
      if(si_ocn%grid%im_world .ne. im) then
        call stop_model(
     &       'TRACER_IC: tracers in sea ice are no longer on the '//
     &       'atm. grid - please move the si_ocn references',255)
      endif
#ifndef TRACERS_ATM_ONLY
      si_ocn%trsi(n,:,:,J_0:J_1)=0.
      trlndi(n,:,J_0:J_1,:)=0.
      trsnowli(n,:,J_0:J_1,:)=0.
      tr_w_ij(n,:,:,:,J_0:J_1)=0.
      tr_wsn_ij(n,:,:,:,J_0:J_1)=0.
#endif
#endif
      name=trname(n)
      if (src_dist_index(n)/=0) name=trname(src_dist_base(n))
      select case (name)

        case default
          write(6,*) 'In TRACER_IC:',name,' does not exist '
          call stop_model("TRACER_IC",255)

        case ('Air')
          do l=1,lm
          do j=J_0,J_1
            trm(:,j,l,n) = MA(l,:,j)*axyp(:,j)
          end do; enddo

        case ('SF6','SF6_c','nh5','nh50','e90')

         ! defaults ok

        case ('aoa','aoanh')

         ! defaults ok

        case ('st8025','tape_rec','nh15')

         ! defaults ok

        case ('Be7', 'Be10', 'Pb210', 'Rn222')
          ! defaults ok

        case ('CO2')
          call openunit('CO2_IC',iu_data,.true.,.true.)
          CALL READT_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),CO2IC,0)
          call closeunit(iu_data)
          do l=1,lm         !ppmv==>ppmm
          do j=J_0,J_1
            trm(:,j,l,n) = co2ic(:,j,l)*MA(l,:,j)*axyp(:,j)*1.54d-6
          enddo; enddo

        case ('N2O')
#ifdef TRACERS_SPECIAL_Lerner
          call openunit('N2O_IC',iu_data,.true.,.true.)
C**** ESMF: Each processor reads the global array: N2Oic
          read (iu_data) title,N2Oic     ! unit is PPMM/(M*AXYP)
          call closeunit(iu_data)
          if (AM_I_ROOT()) write(6,*) title,' read from N2O_IC'
          do l=1,lm         !ppmv==>ppmm
          do j=J_0,J_1
            trm(:,j,l,n) = MA(l,:,j)*axyp(:,j)*N2Oic(j,l)
          enddo; enddo
#endif


#ifdef TRACERS_SPECIAL_Shindell
         if(use_rad_n2o <= 0)then
           select case(PI_run)
           case(1)     ; ICfactor=PIratio_N2O
           case default; ICfactor=1.d0
           end select
           do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
             trm(i,j,l,n) = N2OICX(i,j,l)*ICfactor
           end do   ; end do   ; end do
         else
           if (is_set_param('initial_GHG_setup')) then
             call get_param('initial_GHG_setup', initial_GHG_setup)
             if (initial_GHG_setup == 1 .and. itime == itimeI) then
               select case(PI_run)
               case(1)     ; ICfactor=PIratio_N2O
               case default; ICfactor=1.d0
               end select
               do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
                 trm(i,j,l,n) = N2OICX(i,j,l)*ICfactor
               end do   ; end do   ; end do
             else
               if(ghg_yr/=0) then
                 write(ghg_name,'(I4.4)') ghg_yr
               else
                 write(ghg_name,'(I4.4)') modelEclock%getYear()
               endif
               ghg_file='GHG_IC_'//ghg_name
               call openunit(ghg_file,iu_data,.true.,.true.)
               do m=1,3
                 CALL READT8_COLUMN
     &           (grid,iu_data,NAMEUNIT(iu_data),GHG_IN,0)
                 rad_to_file(m,:,I_0:I_1,J_0:J_1)=
     &           ghg_in(:,I_0:I_1,J_0:J_1)
               enddo
               call closeunit(iu_data)
               do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
                 trm(I,J,L,n) = rad_to_file(3,l,i,j)
               end do   ; end do   ; end do
             end if
           endif
         end if
#endif

        case ('CFC11')   !!! should start April 1
          CFC11ic = 268.D-12*136.5/29.029    !268 PPTV
          do l=1,lm
          do j=J_0,J_1
            trm(:,j,l,n) = MA(l,:,j)*axyp(:,j)*CFC11ic
          enddo; enddo
#ifdef TRACERS_SPECIAL_Lerner
C****
C**** Read in first layer distribution; This is used up to about 100 mb
C****
      call openunit('CFCic_Lerner',iu_data,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),icCFC,0)
      call closeunit(iu_data)
C**** Fill in the tracer; above 100 mb interpolate linearly with P to 0 at top
      stratm = 101.9368
      DO J=J_0,J_1
      DO I=I_0,I_1
        trm(i,j,:,n) = 0.
        trm(i,j,:,n) = 0.
        PUP = STRATM*GRAV
        DO LS=LM,1,-1
          PDN = PUP + MA(ls,I,J)*GRAV
          IF(PDN.GT.10000.d0)  GO TO 450
          trm(I,J,LS,N) =
     *      MA(ls,I,J)*AXYP(I,J)*icCFC(i,j)*.5*(PUP+PDN)/10000.d0
          PUP = PDN
        enddo
  450   CONTINUE
        trm(I,J,LS,N) = MA(ls,I,J)*AXYP(I,J)*icCFC(i,j)*
     *    (1.-.5*(10000.-PUP)*(10000.-PUP)/(10000.*(PDN-PUP)))
        DO LT=1,LS-1
          trm(I,J,LT,N) = MA(lt,I,J)*AXYP(I,J)*icCFC(i,j)
        enddo
      enddo; enddo
#endif


        case ('14CO2')   !!! this tracer is supposed to start 10/16
#ifdef TRACERS_SPECIAL_Lerner
          call get_14CO2_IC(ic14CO2)
          do l=1,lm         !ppmv==>ppmm
          do j=J_0,J_1
            trm(:,j,l,n) = MA(l,:,j)*axyp(:,j)*ic14CO2(:,j,l)*1.d-18
          enddo; enddo
#endif

        case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
         if(use_rad_ch4 <= 0)then
          select case (fix_CH4_chemistry)
          case default
            call get_CH4_IC(0) ! defines trm(:,:,:,n_CH4) within
          case(-1) ! ICs from file...
            call get_CH4_IC(0) ! defines trm(:,:,:,n_CH4) within
            do l=ls1,lm; do j=J_0,J_1; do i=I_0,I_1
              trm(I,J,L,n) = CH4ICX(I,J,L)
            end do   ; end do   ; end do
          end select
#ifdef INTERACTIVE_WETLANDS_CH4
          first_mod(:,:,:)=1
          first_ncep(:)=1
          avg_model(:,:,:)=0.d0
          avg_ncep(:,:,:)=0.d0
          PRS_ch4(:,:,:)=0.d0
          sum_ncep(:,:,:)=0.d0
#endif
         else
           if (is_set_param('initial_GHG_setup')) then
             call get_param('initial_GHG_setup', initial_GHG_setup)
             if (initial_GHG_setup == 1 .and. itime == itimeI) then
               select case (fix_CH4_chemistry)
               case default
                 call get_CH4_IC(0) ! defines trm(:,:,:,n_CH4) within
               case(-1)         ! ICs from file...
                 call get_CH4_IC(0) ! defines trm(:,:,:,n_CH4) within
                 do l=ls1,lm; do j=J_0,J_1; do i=I_0,I_1
                   trm(I,J,L,n) = CH4ICX(I,J,L)
                 end do   ; end do   ; end do
               end select
             else
               if(ghg_yr/=0) then
                 write(ghg_name,'(I4.4)') ghg_yr
               else
                 write(ghg_name,'(I4.4)') modelEclock%getYear()
               endif
               ghg_file='GHG_IC_'//ghg_name
               call openunit(ghg_file,iu_data,.true.,.true.)
               do m=1,4
                 CALL READT8_COLUMN(grid,iu_data,NAMEUNIT(iu_data),
     &                GHG_IN,0)
                 rad_to_file(m,:,I_0:I_1,J_0:J_1)=
     &                ghg_in(:,I_0:I_1,J_0:J_1)
               enddo
               call closeunit(iu_data)
               do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
                 trm(I,J,L,n) = rad_to_file(4,l,i,j)
               end do   ; end do   ; end do
             end if
           end if
         end if
         do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
           chem_tracer_save(2,L,I,J)=trm(I,J,L,n)
     &          *byaxyp(i,j)*avog/(tr_mm(n)*2.69e20) ! to atm*cm
         end do   ; end do   ; end do
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_SPECIAL_Lerner
          call get_wofsy_gas_IC(name,CH4ic)
          do l=1,lm         !ppbv==>ppbm
          do j=J_0,J_1
            trm(:,j,l,n) = MA(l,:,j)*axyp(:,j)*CH4ic(j,l)*0.552d-9
          enddo; enddo
#endif

        case ('O3')
          do l=1,lm
          do j=J_0,J_1
            trm(:,j,l,n) = MA(l,:,j)*axyp(:,j)*20.d-9*vol2mass(n)
          enddo; enddo
#ifdef TRACERS_SPECIAL_Lerner
          do l=lm,lm+1-nstrtc,-1
          lr = lm+1-l
            do j=J_0,J_1
            if (tlt0m(j,lr,5) /= 0.) then
            trm(:,j,l,n) =
     *          tlt0m(j,lr,1)*MA(l,:,j)*axyp(:,j)*vol2mass(n)
            trmom(mz,:,j,l,n)  =
     *          tltzm(j,lr,1)*MA(l,:,j)*axyp(:,j)*vol2mass(n)
            trmom(mzz,:,j,l,n)  =
     *         tltzzm(j,lr,1)*MA(l,:,j)*axyp(:,j)*vol2mass(n)
            end if
            end do
          end do
#endif

#ifdef TRACERS_WATER
      case ('Water', 'H2O18', 'HDO', 'HTO', 'H2O17')

C**** initial atmospheric conc. needs to be defined for each tracer
        select case (name)
        case ('Water')
          trinit=1.
C**** for gradients defined on air mass
          tmominit = 1.
C**** for gradients defined on water mass (should be an option?)
c     tmominit = 0.
        case ('H2O18')        ! d18O=-80
          trinit=0.92d0*trw0(n)
          tmominit = trinit
        case ('H2O17')        ! d18O=-43.15  (D17O=0)
          trinit=.95685d0*trw0(n)
          tmominit = trinit
        case ('HDO')   ! dD=-630
          trinit=0.37d0*trw0(n)
          tmominit = trinit
        case ('HTO')
          trinit=0.
          tmominit = trinit
        end select

        do l=1,lm
        do j=J_0,J_1
          do i=I_0,I_1
            trm(i,j,l,n) =  q(i,j,l)*MA(l,i,j)*axyp(i,j)*trinit
            trwm(i,j,l,n)= (qcl(i,j,l)+qci(i,j,l))*MA(l,i,j)
     &           *axyp(i,j)*trinit
            trmom(:,i,j,l,n) = qmom(:,i,j,l)*MA(l,i,j)*axyp(i,j)
     *           *tmominit
            if (src_dist_index(n)/=0) then
              trm(i, j, l, n)=trm(i, j, l, n)*
     &                                 xyztr(src_dist_index(n), i, j)
              trwm(i, j, l, n)=trwm(i, j, l, n)*
     &                                 xyztr(src_dist_index(n), i, j)
              trmom(:, i, j, l, n)=trmom(:, i, j, l, n)*
     &                                 xyztr(src_dist_index(n), i, j)
            endif
          end do
        end do
        end do
        if (HAVE_SOUTH_POLE) then
           do i=2,im
              trm(i,1,:,n) =  trm(1,1,:,n) !poles
              trwm(i, 1,:,n)= trwm(1, 1,:,n) !poles
              trmom(:,i, 1,:,n)=0.
           enddo
        endif
        if (HAVE_NORTH_POLE) then
           do i=2,im
              trm(i,jm,:,n) = trm(1,jm,:,n) !poles
              trwm(i,jm,:,n)= trwm(1,jm,:,n) !poles
              trmom(:,i,jm,:,n)=0.
          enddo
        endif
        if (name.eq."HTO") then ! initialise bomb source
          do l=ls1-1,ls1+1      ! strat. source lat 44 N - 56 N
          do j=J_0,J_1
          do i=I_0,I_1
            if(nint(lat2d_dg(i,j)).ge.44.and.nint(lat2d_dg(i,j)).le.56)
     *           trm(i,j,l,n)= q(i,j,l)*MA(l,i,j)*axyp(i,j)*1d10*1d-18
          end do
          end do
          end do
        end if

#ifndef TRACERS_ATM_ONLY
        call init_single_seaice_tracer(si_atm,n,trsi0(n))
        call init_single_seaice_tracer(si_ocn,n,trsi0(n))

        do j=J_0,J_1
          do i=I_0,I_1
            tracerTs=trw0(n)
#ifdef TRACERS_SPECIAL_O18
c Define a simple d18O based on Tsurf for GIC, put dD on meteoric water line
            if(name.eq."H2O18") tracerTs=TRW0(n_H2O18)*(1.+1d-3*
     *           ((atmsrf%tsavg(i,j)-(tf+tracerT0))*d18oT_slope))
            if(name.eq."HDO") tracerTs=TRW0(n_HDO)*(1.+(1d-3*
     *         (((atmsrf%tsavg(i,j)-(tf+tracerT0))*d18oT_slope)*8+1d1)))
#endif
C**** lakes
            if (flake(i,j).gt.0) then
              trlake(n,1,i,j)=tracerTs*mldlk(i,j)*rhow*flake(i,j)
     *             *axyp(i,j)
              if (mwl(i,j)-mldlk(i,j)*rhow*flake(i,j)*axyp(i,j).gt.1d-10
     *             *mwl(i,j)) then
                trlake(n,2,i,j)=tracerTs*mwl(i,j)-trlake(n,1,i,j)
              else
                trlake(n,2,i,j)=0.
              end if
              atmocn%gtracer(n,i,j)=trw0(n)
            else !if (focean(i,j).eq.0) then
              trlake(n,1,i,j)=trw0(n)*mwl(i,j)
              trlake(n,2,i,j)=0.
c            else
c              trlake(n,1:2,i,j)=0.
            end if
c**** ice
            if (si_atm%msi(i,j).gt.0) then
              atmice%gtracer(n,i,j)=trsi0(n)
            end if
c**** landice
            if (flice(i,j).gt.0) then
              trlndi(n,i,j,:)=trli0(n)*(ace1li+ace2li)	! calls trli0_s()
              trsnowli(n,i,j,:)=trli0(n)*snowli(i,j,:)
              do ipatch=1,ubound(atmglas,1)
#ifdef GLINT2
                atmglas_hp(ipatch)%gtracer(n,i,j)=trli0(n)
#endif
                atmglas(ipatch)%gtracer(n,i,j)=trli0(n)
              enddo
            else
              trlndi(n,i,j,:)=0.
              trsnowli(n,i,j,:)=0.
              do ipatch=1,ubound(atmglas,1)
#ifdef GLINT2
                atmglas_hp(ipatch)%gtracer(n,i,j)=0.
#endif
                atmglas(ipatch)%gtracer(n,i,j)=0.
              enddo
            end if
c**** earth
            !!!if (fearth(i,j).gt.0) then
            if (focean(i,j) < 1.d0) then
              conv=rhow         ! convert from m to kg m-2
              tr_w_ij  (n,:,:,i,j)=tracerTs*w_ij (:,:,i,j)*conv
              tr_wsn_ij(n,1:nsn_ij(1,i,j),1,i,j)=
     &             tracerTs*wsn_ij(1:nsn_ij(1,i,j),1,i,j)
     &             *fr_snow_ij(1,i,j)*conv
              tr_wsn_ij(n,1:nsn_ij(2,i,j),2,i,j)=
     &             tracerTs*wsn_ij(1:nsn_ij(2,i,j),2,i,j)
     &             *fr_snow_ij(2,i,j)*conv
              !trsnowbv(n,2,i,j)=trw0(n)*snowbv(2,i,j)*conv
              atmlnd%gtracer (n,i,j)=trw0(n)
            else
              tr_w_ij  (n,:,:,i,j)=0.
              tr_wsn_ij(n,:,:,i,j)=0.
              !trsnowbv(n,1,i,j)=0.
              !trsnowbv(n,2,i,j)=0.
              atmlnd%gtracer(n,i,j)=0.
            end if
          end do
          end do
#endif
#ifdef TRACERS_SPECIAL_O18
          if (AM_I_ROOT()) then
            if(name.eq."H2O18") write(6,'(A52,f6.2,A15,f8.4,A18)')
     *            "Initialized trlake tr_w_ij tr_wsn_ij using Tsurf at"
     *           ,tracerT0,"degC, 0 permil",d18oT_slope
     *           ,"permil d18O/degC"
          endif
#endif

#endif /* TRACERS_WATER */

#ifdef TRACERS_SPECIAL_Shindell
        case ('Ox')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(I,J,L,n) = OxIC(I,J,L)
            chem_tracer_save(1,L,I,J)=OxIC(I,J,L)*byO3MULT*byaxyp(i,j)
          end do   ; end do   ; end do
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        case ('stratOx')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(I,J,L,n) = OxIC(I,J,L)
            stratO3_tracer_save(L,I,J)=OxIC(I,J,L)*byO3MULT*byaxyp(i,j)
          end do   ; end do   ; end do
#endif

        case ('NOx')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*1.d-11*ICfactor
            if(PRES(L).lt.10.)trm(i,j,l,n)=trm(i,j,l,n)*3.d2
          end do; end do; end do

        case ('ClOx')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11*ClOxalt(l)
          end do; end do; end do

        case ('BrOx')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11*BrOxalt(l)
          end do; end do; end do

        case ('HCl')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11*HClalt(l)
          end do; end do; end do

        case ('ClONO2')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11*ClONO2alt(l)
          end do; end do; end do

        case ('N2O5')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*1.d-12*ICfactor
          end do; end do; end do

        case ('HNO3')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*1.d-10*ICfactor
            if(PRES(L).lt.50.and.PRES(L).gt.10.)
     &      trm(i,j,l,n)=trm(i,j,l,n)*1.d2
          end do; end do; end do
#endif /* TRACERS_SPECIAL_Shindell */

        case ('H2O2')
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*5.d-10
          end do; end do; end do

#ifdef SHINDELL_STRAT_EXTRA
        case ('GLT')
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = GLTic*vol2mass(n)*MA(l,i,j)*axyp(i,j)
          end do; end do; end do
#endif

#ifdef TRACERS_SPECIAL_Shindell
        case ('CH3OOH',
#ifdef TRACERS_dCO
     *        'dMe17OOH', 'dMe18OOH', 'd13MeOOH',
     *        'dHCH17O', 'dHCH18O', 'dH13CHO',
#endif  /* TRACERS_dCO */
     &        'HCHO')
          select case (trname(n))
#ifdef TRACERS_dCO
            case ('dMe17OOH', 'dMe18OOH', 'd13MeOOH')
              dICfactor=dMeOOH_IC_fact
            case ('dHCH17O', 'dHCH18O', 'dH13CHO')
              dICfactor=dHCHO_IC_fact
#endif  /* TRACERS_dCO */
            case default
              dICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*1.d-11*dICfactor
          end do; end do; end do

        case ('HO2NO2')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*1.d-12*ICfactor
          end do; end do; end do

        case ('CO'
#ifdef TRACERS_dCO
     *       ,'dC17O','dC18O','d13CO'
#endif  /* TRACERS_dCO */
     *       )
          select case (trname(n))
#ifdef TRACERS_dCO
            case ('dC17O','dC18O','d13CO')
              dICfactor=dCO_IC_fact
#endif  /* TRACERS_dCO */
            case default
              dICfactor=1.d0
          end select
          do l=1,lm
            select case(PI_run)
            case(1) ! ise scaling
              if(L.le.LS1-1) then
                ICfactor=PIratio_CO_T ! troposphere
              else
                ICfactor=PIratio_CO_S ! stratosphere
              end if
            case default; ICfactor=1.d0
            end select
            do j=J_0,J_1; do i=I_0,I_1
              trm(I,J,L,n) = COIC(I,J,L)*ICfactor*dICfactor
            end do   ; end do
          end do

        case ('codirect')
          ! supposed to start from zero conc for the HTAP TP anyway...
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) = 0.d0
          end do; end do; end do

        case ('PAN'
#ifdef TRACERS_dCO
     *       ,'d17OPAN','d18OPAN','d13CPAN'
#endif  /* TRACERS_dCO */
     *       )
          select case (trname(n))
#ifdef TRACERS_dCO
            case ('d17OPAN','d18OPAN','d13CPAN')
              dICfactor=dPAN_IC_fact
#endif  /* TRACERS_dCO */
            case default
              dICfactor=1.d0
          end select
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*4.d-11*ICfactor*dICfactor
          end do; end do; end do

        case ('Isoprene')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*0.d-11*ICfactor
          end do; end do; end do

        case ('AlkylNit')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*2.d-10*ICfactor
          end do; end do; end do

        case('Alkenes'
#ifdef TRACERS_dCO
     *      ,'d13Calke'
#endif  /* TRACERS_dCO */
     *       )
          select case (trname(n))
#ifdef TRACERS_dCO
            case ('d13Calke')
              dICfactor=dalke_IC_fact
#endif  /* TRACERS_dCO */
            case default
              dICfactor=1.d0
          end select
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*4.d-10*ICfactor*dICfactor
          end do; end do; end do

        case('Paraffin'
#ifdef TRACERS_dCO
     *      ,'d13CPAR'
#endif  /* TRACERS_dCO */
     *       )
          select case (trname(n))
#ifdef TRACERS_dCO
            case ('d13CPAR')
              dICfactor=dPAR_IC_fact
#endif  /* TRACERS_dCO */
            case default
              dICfactor=1.d0
          end select
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*5.d-10*ICfactor*dICfactor
          end do; end do; end do

        case('Terpenes'
#ifdef TRACERS_AEROSOLS_SOA
     &      ,'isopp1g','isopp1a','isopp2g','isopp2a'
     &      ,'apinp1g','apinp1a','apinp2g','apinp2a'
#endif
#ifdef TRACERS_AEROSOLS_VBS
     *      ,'vbsGm2', 'vbsGm1', 'vbsGz',  'vbsGp1', 'vbsGp2'
     *      ,'vbsGp3', 'vbsGp4', 'vbsGp5', 'vbsGp6'
     *      ,'vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2'
     *      ,'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6'
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_OCEAN
     &      ,'OCocean'
#endif
     &      )
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      MA(l,i,j)*axyp(i,j)*vol2mass(n)*0.d0*ICfactor*5.d-14
          end do; end do; end do
#endif /* TRACERS_SPECIAL_Shindell */

        case ('CO2n')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
             !units: [am]=kg_air/m2, [axyp]=m2, [tr_mm]=kg_CO2,
             !       [bymair]=1/kg_air, [atmCO2]=ppmv=10^(-6)kg_CO2/kg_air
             !       [vol2mass]=(gr,CO2/moleCO2)/(gr,air/mole air)
             trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*vol2mass(n)
     .                    * get_atmCO2()*1.d-6
             atmocn%gtracer(n,i,j) = vol2mass(n)
     .                    * get_atmCO2() * 1.d-6      !initialize gtracer
          end do; end do; end do
        case ('CFCn')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            if(l.ge.LS1) then
              trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*vol2mass(n)*2.d-13
            else
              trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-13
            end if
          end do; end do; end do

#ifdef TRACERS_SPECIAL_Shindell
        case ('CFC')
         if(use_rad_cfc.le.0)then
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_CFC
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(I,J,L,n) = CFCIC(I,J,L)*ICfactor
          end do   ; end do   ; end do
         else
           if (is_set_param('initial_GHG_setup')) then
             call get_param('initial_GHG_setup', initial_GHG_setup)
             if (initial_GHG_setup == 1 .and. itime == itimeI) then
               select case(PI_run)
             case(1)     ; ICfactor=PIratio_CFC
               case default; ICfactor=1.d0
             end select
             do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
               trm(I,J,L,n) = CFCIC(I,J,L)*ICfactor
             end do   ; end do   ; end do
           else
             if(ghg_yr/=0)then; write(ghg_name,'(I4.4)') ghg_yr
             else; write(ghg_name,'(I4.4)') modelEclock%getYear(); endif
             ghg_file='GHG_IC_'//ghg_name
             call openunit(ghg_file,iu_data,.true.,.true.)
             do m=1,5
               CALL READT8_COLUMN(grid,iu_data,NAMEUNIT(iu_data),GHG_IN,
     &              0)
               rad_to_file(m,:,I_0:I_1,J_0:J_1)=
     &              ghg_in(:,I_0:I_1,J_0:J_1)
             enddo
             call closeunit(iu_data)
             do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
               trm(I,J,L,n) = rad_to_file(5,l,i,j)
             end do   ; end do   ; end do
           end if
         endif
       end if
#endif /* TRACERS_SPECIAL_Shindell */

        case ('BrONO2','HBr','HOBr')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            if(l.ge.LS1) then
              trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*vol2mass(n)*2.d-13
            else
              trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-13
            end if
          end do; end do; end do

        case ('HOCl')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            if(l.ge.LS1) then
              trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*vol2mass(n)*5.d-11
            else
              trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11
            end if
          end do; end do; end do

        case('DMS')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*vol2mass(n)*5.d-13
          end do; end do; end do

#ifndef TRACERS_TOMAS
        case('MSA', 'SO2', 'SO4', 'SO4_d1', 'SO4_d2', 'SO4_d3',
     *         'N_d1','N_d2','N_d3','NH3','NH4','NO3p',
     *         'BCII', 'BCIA', 'BCB', 'OCII', 'OCIA', 'OCB', 'H2O2_s',
     *         'seasalt1', 'seasalt2',
     *         'M_NO3   ','M_NH4   ','M_H2O   ','N_AKK_1 ',
     *         'N_ACC_1 ','M_DD1_SU','N_DD1_1 ',
     *         'M_DS1_SU','M_DS1_DU','N_DS1_1 ','M_DD2_SU','M_DD2_DU',
     *         'N_DD2_1 ','M_DS2_SU','M_DS2_DU','N_DS2_1 ','M_SSA_SU',
     *         'M_OCC_SU','N_OCC_1 ','M_BC1_SU','N_SSA_1 ','N_SSC_1 ',
     *         'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_BC3_SU',
     *         'M_BC3_BC','N_BC3_1 ','M_DBC_SU','M_DBC_BC','M_DBC_DU',
     *         'N_DBC_1 ','M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',
     *         'M_BCS_SU','M_BCS_BC','N_BCS_1 ','M_MXX_SU','M_MXX_BC',
     *         'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','M_OCS_SU',
     *         'M_OCS_OC','N_OCS_1 ','H2SO4',
     *         'M_AKK_SU','M_ACC_SU','M_DD1_DU',
     *         'M_SSA_SS','M_SSC_SS','M_BC1_BC','M_OCC_OC',
     *         'M_SSS_SS','M_SSS_SU')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) = MA(l,i,j)*axyp(i,j)*vol2mass(n)*5.d-14
          end do; end do; end do
#endif
#ifdef TRACERS_TOMAS
        case('SO2','NH3','NH4','H2SO4','SOAgas','H2O2_s')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =MA(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-30
          end do; end do; end do

       case('ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     *    'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     *    'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15',
     *    'ANACL_01','ANACL_02','ANACL_03','ANACL_04','ANACL_05',
     *    'ANACL_06','ANACL_07','ANACL_08','ANACL_09','ANACL_10',
     *    'ANACL_11','ANACL_12','ANACL_13','ANACL_14','ANACL_15',
     *    'AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05',
     *    'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10',
     *    'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15',
     *    'AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05',
     *    'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10',
     *    'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15',
     *    'AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05',
     *    'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10',
     *    'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15',
     *    'AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05',
     *    'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10',
     *    'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15',
     *    'ADUST_01','ADUST_02','ADUST_03','ADUST_04','ADUST_05',
     *    'ADUST_06','ADUST_07','ADUST_08','ADUST_09','ADUST_10',
     *    'ADUST_11','ADUST_12','ADUST_13','ADUST_14','ADUST_15',
     *    'AH2O__01','AH2O__02','AH2O__03','AH2O__04','AH2O__05',
     *    'AH2O__06','AH2O__07','AH2O__08','AH2O__09','AH2O__10',
     *    'AH2O__11','AH2O__12','AH2O__13','AH2O__14','AH2O__15')

          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
                trm(i,j,l,n) =MA(l,i,j)*axyp(i,j)*1.d-20
          end do; end do; end do

      CASE('ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     *    'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     *    'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15')

           k=n-n_ANUM(1)+1

          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
                trm(i,j,l,n) =MA(l,i,j)*axyp(i,j)*7.d-20
     &               /(sqrt(xk(k+1)*xk(k))) !MA(l,i,j)*axyp(i,j)*vol2mass(n)*5.d-14
          end do; end do; end do
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
        CASE('Clay','Silt1','Silt2','Silt3','Silt4','Silt5','ClayIlli'
     &         ,'ClayKaol','ClaySmec','ClayCalc','ClayQuar','ClayFeld'
     &         ,'ClayHema','ClayGyps','ClayIlHe','ClayKaHe','ClaySmHe'
     &         ,'ClayCaHe','ClayQuHe','ClayFeHe','ClayGyHe','Sil1Quar'
     &         ,'Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps','Sil1Illi'
     &         ,'Sil1Kaol','Sil1Smec','Sil1QuHe','Sil1FeHe','Sil1CaHe'
     &         ,'Sil1GyHe','Sil1IlHe','Sil1KaHe','Sil1SmHe','Sil2Quar'
     &         ,'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps','Sil2Illi'
     &         ,'Sil2Kaol','Sil2Smec','Sil2QuHe','Sil2FeHe','Sil2CaHe'
     &         ,'Sil2GyHe','Sil2IlHe','Sil2KaHe','Sil2SmHe','Sil3Quar'
     &         ,'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps','Sil3Illi'
     &         ,'Sil3Kaol','Sil3Smec','Sil3QuHe','Sil3FeHe','Sil3CaHe'
     &         ,'Sil3GyHe','Sil3IlHe','Sil3KaHe','Sil3SmHe','Sil4Quar'
     &         ,'Sil4Feld','Sil4Calc','Sil4Hema','Sil4Gyps','Sil4Illi'
     &         ,'Sil4Kaol','Sil4Smec','Sil4QuHe','Sil4FeHe','Sil4CaHe'
     &         ,'Sil4GyHe','Sil4IlHe','Sil4KaHe','Sil4SmHe','Sil5Quar'
     &         ,'Sil5Feld','Sil5Calc','Sil5Hema','Sil5Gyps','Sil5Illi'
     &         ,'Sil5Kaol','Sil5Smec','Sil5QuHe','Sil5FeHe','Sil5CaHe'
     &         ,'Sil5GyHe','Sil5IlHe','Sil5KaHe','Sil5SmHe')
          ! defaults ok
          hbaij=0D0
          ricntd=0D0
#endif

      end select

C**** Initialise pbl profile if necessary
      if (needtrs(n)) then
        do ipatch=1,size(asflx)
        do j=J_0,J_1
        do ipbl=1,npbl
#ifdef TRACERS_WATER
          if(tr_wd_type(n).eq.nWATER)THEN
            asflx(ipatch)%trabl(ipbl,n,:,j) =
     &           trinit*asflx(ipatch)%qabl(ipbl,:,j)
            if (src_dist_index(n)/=0) asflx(ipatch)%trabl(ipbl,n,:,j) =
     &              asflx(ipatch)%trabl(ipbl,n,:,j)*
     &                                 xyztr(src_dist_index(n), :, j)
          ELSE
            asflx(ipatch)%trabl(ipbl,n,:,j) =
     &           trm(:,j,1,n)*byMA(1,:,j)*byaxyp(:,j)
          END IF
#else
            asflx(ipatch)%trabl(ipbl,n,:,j) =
     &         trm(:,j,1,n)*byMA(1,:,j)*byaxyp(:,j)
#endif
        end do
        end do
        end do
      end if

      write(out_line,*) ' Tracer ',trname(n),' initialized at itime='
     *     ,itime
      call write_parallel(trim(out_line))

      end if
      end do
#endif /* TRACERS_ON */
C****
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
c read in DMS source
      call openunit('DMS_SEA',iudms,.true.,.true.)
      DMSinput(:,:,:)= 0.d0
      do mm=1,12
        call readt_parallel(grid,iudms,nameunit(iudms),
     *                      DMSinput(:,:,mm),0)
      end do
      call closeunit(iudms)
 901  FORMAT(3X,3(I4),E11.3)

c read in SO2 emissions
c volcano - continuous
C    Initialize:
      if (iso2volcano>0) then
      so2_src_3D(:,:,:,iso2volcano)= 0.d0
c read lat-lon netcdf file and convert lat,lon,pres to i,j,l.
c NOTE: the input file specifies integrals over its gridboxes.
      ALLOCATE(  psref(grid%i_strt_halo:grid%i_stop_halo,
     &                 grid%j_strt_halo:grid%j_stop_halo) )
      call openunit('PSREF',iu_ps,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_ps,'PSREF',psref,1)
      call closeunit(iu_ps)
      status = nf_open('SO2_VOLCANO',nf_nowrite,file_id)
      status = nf_inq_varid(file_id,'lon',vid)
      status = nf_get_var_double(file_id,vid,volc_lons)
      status = nf_inq_varid(file_id,'lat',vid)
#ifndef CUBED_SPHERE
      status = nf_inq_dimlen(file_id,vid,lat_val)
      allocate(volc_lats(lat_val))
      allocate(volc_pup(im,lat_val))
      allocate(volc_emiss(im,lat_val))
#endif
      status = nf_get_var_double(file_id,vid,volc_lats)
      status = nf_inq_varid(file_id,'Pres_CONTmax',vid)
      status = nf_get_var_double(file_id,vid,volc_pup)
      status = nf_inq_varid(file_id,'VOLC_CONT',vid)
      status = nf_get_var_double(file_id,vid,volc_emiss)
      status = nf_close(file_id)
      do jlat=1,ubound(volc_lats,1)
        do ilon=1,ubound(volc_lons,1)
          if(volc_emiss(ilon,jlat) <= 0.) cycle
          call lonlat_to_ij(
     &         (/volc_lons(ilon),volc_lats(jlat)/),volc_ij)
          ii = volc_ij(1); jj = volc_ij(2)
          if(jj<j_0 .or. jj>j_1) cycle
          if(ii<i_0 .or. ii>i_1) cycle
          Call CALC_VERT_AMP (psref(ii,jj),lm, amref,x1d,pednref,x1d)
          lmax = 1
          do while(pednref(lmax) > volc_pup(ilon,jlat))
            lmax = lmax + 1
          enddo
          amsum = sum(amref(1:lmax))
          do ll=1,lmax ! add source between surf and max height
            so2_src_3d(ii,jj,ll,iso2volcano)=
     &        so2_src_3d(ii,jj,ll,iso2volcano)
     &          +(amref(ll)/amsum)*
     &           volc_emiss(ilon,jlat)/(SECONDS_PER_DAY*30.4d0)/12.d0
          enddo
        enddo
      enddo
#ifndef CUBED_SPHERE
      deallocate(volc_lats, volc_pup, volc_emiss)
#endif
      deallocate(psref)
      endif ! iso2volcano>0
#endif
! ---------------------------------------------------
#ifndef TRACERS_AEROSOLS_SOA
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
c Terpenes
      OCT_src(:,:,:)=0.d0
      call openunit('Terpenes_01',iuc,.true.,.true.)
      call skip_parallel(iuc)
      do mm=1,12
        call readt_parallel(grid,iuc,nameunit(iuc),OCT_src(:,:,mm),0)
      end do
      call closeunit(iuc)
c units are mg Terpene/m2/month
      do i=I_0,I_1; do j=J_0,J_1; do mm=1,12
! 10% of terpenes end up being SOA
#ifdef TRACERS_TOMAS
        OCT_src(i,j,mm)=OCT_src(i,j,mm)*axyp(i,j)*0.1d0
     +                  *om2oc(n_AOCOB(1))
#else
#ifdef TRACERS_AMP
        OCT_src(i,j,mm)=OCT_src(i,j,mm)*axyp(i,j)*0.1d0
     +                  *om2oc(n_M_OCC_OC)
#else
        OCT_src(i,j,mm)=OCT_src(i,j,mm)*axyp(i,j)*0.1d0
     +                  *om2oc(n_OCII)
#endif
#endif
      end do; end do; end do
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
! ---------------------------------------------------
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
c **** reads in files for dust/mineral tracers
      call tracer_ic_soildust
#endif

      end subroutine tracer_IC


      subroutine daily_tracer(end_of_day)
!@sum daily_tracer is called once a day for tracers
!@+   SUBROUTINE tracer_IC is called from daily_tracer to allow for
!@+     tracers that 'turn on' on different dates.
!@auth Jean Lerner
C**** Note this routine must always exist (but can be a dummy routine)
      USE RESOLUTION, only : lm
      Use ATM_COM,    Only: MA,T
      use model_com, only: modelEclock
      USE MODEL_COM, only:itime
      USE FLUXES, only : fearth0,focean,flake0
      USE SOMTQ_COM, only : tmom,mz
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds,
     & write_parallel
      USE RAD_COM, only: o3_yr
#ifdef TRACERS_VOLCEXP
      USE AEROSOL_SOURCES, only: so2_src_3d,iso2volcanoexpl
      USE timestream_mod, only: init_stream,read_stream
      USE tracer_com, only: SO2_volc_stream,SO2_vphe_stream
#endif
      use GEOM, only: lat_to_j
      use GEOM, only: lon_to_i
      USE ATM_COM, only: byMA
      USE GEOM, only: byaxyp
#ifdef TRACERS_SPECIAL_Shindell
      use TRCHEM_Shindell_COM, only: NOx_yr
      use TRCHEM_Shindell_COM, only: CO_yr
      use TRCHEM_Shindell_COM, only: VOC_yr
#endif  /* TRACERS_SPECIAL_Shindell */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      use TRACER_COM, only: aer_int_yr
      use TRACER_COM, only: SO2_int_yr
      use TRACER_COM, only: NH3_int_yr
      use TRACER_COM, only: BC_int_yr
      use TRACER_COM, only: OC_int_yr
      use TRACER_COM, only: ex_volc_num
      use TRACER_COM, only: ex_volc_jday
      use TRACER_COM, only: ex_volc_year
      use TRACER_COM, only: ex_volc_lat
      use TRACER_COM, only: ex_volc_lon
      use TRACER_COM, only: ex_volc_bot
      use TRACER_COM, only: ex_volc_top
      use TRACER_COM, only: ex_volc_SO2
      use TRACER_COM, only: ex_volc_H2O
      USE AEROSOL_SOURCES, only: so2_src_3d,iso2exvolc,H2O_src_3d
#endif
#ifdef CUBED_SPHERE
      USE tracer_com, only: AIRCstreams
#endif

#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : variable_phi
#endif
      USE CONSTANT, only: grav
      use RunTimeControls_mod, only: tracers_amp
      use RunTimeControls_mod, only: tracers_tomas
      use RunTimeControls_mod, only: tracers_aerosols_soa
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use OldTracer_mod, only: trname, itime_tr0, MAX_LEN_NAME
      use OldTracer_mod, only: nBBsources,do_fire,vol2mass,do_aircraft
      use TRACER_COM, only: tracers, set_ntsurfsrc
      USE TRACER_COM, only: coupled_chem,daily_z
      USE TRACER_COM, only: n_CO2n
      USE TRACER_COM, only: NTM,nOther,nAircraft,
     & n_CH4,n_Isoprene,n_codirect,sfc_src,ntsurfsrc,
     & trans_emis_overr_yr,trans_emis_overr_day
      use TRACER_COM, only: ntm_chem_beg,ntm_chem_end
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      use TRACER_COM, only:
     *  n_NH3,n_SO2,n_SO4,n_BCII,n_BCB,n_OCII,n_OCB
     * ,n_M_ACC_SU,n_M_AKK_SU,n_M_BC1_BC,n_M_OCC_OC,n_M_BOC_BC
     * ,n_M_BOC_OC
#ifdef TRACERS_TOMAS
      use TRACER_COM, only:
     * n_AECOB,n_AOCOB,n_ASO4,nbins,n_AECOB,n_ADUST
#endif
#endif
#ifdef TRACERS_SPECIAL_Lerner
      use tracer_com, only: n_O3,n_CO2,n_CH4
      USE TRACERS_MPchem_COM, only: STRATCHEM_SETUP
      USE LINOZ_CHEM_COM, only: LINOZ_SETUP
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE FLUXES, only: tr3Dsource
      USE TRCHEM_Shindell_COM,only:
     & dms_offline,so2_offline,sulfate,fix_CH4_chemistry
      USE TRCHEM_Shindell_COM, only: tune_NOx
      USE TRCHEM_Shindell_COM, only: tune_BVOC
      use photolysis, only: rad_FL,read_FL
#endif
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : variable_phi
#endif
      use Tracer_mod, only: Tracer, readSurfaceSources
      IMPLICIT NONE
      INTEGER n,last_month,kk,nread,xday,xyear,ns
      logical :: checkname
      LOGICAL, INTENT(IN) :: end_of_day
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
     &     :: daily_gz
      data last_month/-1/
      INTEGER J_0, J_1, I_0, I_1,I,J,ll,lmax,lmin
#ifdef TRACERS_TOMAS
      integer km, najl_num,naij_num,k
      real*8 :: scalesize(nbins+nbins) !temporal emission mass fraction
#endif
#ifdef TRACERS_VOLCEXP
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO)
     &     :: SO2_volc_emis_expl, Plume_hei_volc_emis_expl !  volc emiss
#endif
      integer :: ex
      real*8 :: dz
      class (Tracer), pointer :: pTracer
C****
      integer :: year, month, dayOfYear
      character(len=MAX_LEN_NAME) :: tmpString
      logical :: isChemTracer

      call modelEclock%get(year=year, month=month,
     *     dayOfYear=dayOfYear)
CC****
C**** Extract useful local domain parameters from "grid"
C****

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      if(end_of_day) then
         Call COMPUTE_GZ (MA,T,TMOM(MZ,:,:,:), DAILY_Z)
        daily_z = daily_z/grav
      endif
      daily_gz = grav*daily_z

#ifdef TRACERS_VOLCEXP
! Reading explosive volcano emissions for SO2
      if(.not. end_of_day) then ! synonym for model init phase
        ! initialize the file handle

        call init_stream(grid,SO2_volc_stream,'SO2_VOLCANO_EXPL','SO2'
     &       ,0d0,1d30,'linm2m',year,dayofyear)

        call init_stream(grid,SO2_vphe_stream,'SO2_VOLCANO_EXPL',
     &       'Plume_height',0d0,1d30,'linm2m',year,dayofyear)
      endif

      call read_stream(grid,SO2_volc_stream,year,dayofyear,
     &                 SO2_volc_emis_expl)
      call read_stream(grid,SO2_vphe_stream,year,dayofyear,
     &                 Plume_hei_volc_emis_expl)

      so2_src_3d(:,:,:,iso2volcanoexpl) = 0.d0

      DO J=J_0,J_1
      DO I=I_0,I_1

        if(so2_volc_emis_expl(i,j) <= 0.d0) cycle
          lmax = 1
          do while(daily_z(i,j,lmax) < Plume_hei_volc_emis_expl(i,j))
            lmax = lmax + 1
          enddo
            lmax = lmax + 1 ! adding one layer as to not have plume height identical to mixing height
            lmin=max(1,lmax - lmax/3)
          do ll=lmin,lmax ! add source into the upper 1/3 of the plume
                          ! conversion kt d-1 into kg s-1
          if (lmax <= 2) then
            so2_src_3d(i,j,1,iso2volcanoexpl)=
     &        so2_src_3d(i,j,1,iso2volcanoexpl)
     &                + so2_volc_emis_expl(i,j)/SECONDS_PER_DAY*1.d6
          else
            so2_src_3d(i,j,ll,iso2volcanoexpl)=
     &        so2_src_3d(i,j,ll,iso2volcanoexpl)
     &                + (1./(float(lmax-lmin)+1))
     &                * so2_volc_emis_expl(i,j)/SECONDS_PER_DAY*1.d6
          endif
          enddo
        enddo
      enddo
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
! explosive volcano injections based on rundeck parameters
      if (iso2exvolc>0) then
        so2_src_3d(:,:,:,iso2exvolc)=0.d0
        H2O_src_3d(:,:,:)=0.d0
        do ex=1,ex_volc_num

! check data
          if (ex_volc_jday(ex)/=dayOfYear) cycle
          if (ex_volc_year(ex)/=year) cycle

! check location
          j=lat_to_j(ex_volc_lat(ex))
          if (j<j_0 .or. j>j_1) cycle
          i=lon_to_i(ex_volc_lon(ex))
          if (i<i_0 .or. i>i_1) cycle

! find layers range
          do lmin=2,lm ! start from 2, to avoid mixing layer
            if (daily_z(i,j,lmin)>=ex_volc_bot(ex)) exit
          enddo
          do lmax=lmin,lm
            if (daily_z(i,j,lmax)>=ex_volc_top(ex)) exit
          enddo

! set emissions
          dz=sum(daily_z(i,j,lmin:lmax))
          do ll=lmin,lmax
            so2_src_3d(i,j,ll,iso2exvolc)=so2_src_3d(i,j,ll,iso2exvolc)+
     &        daily_z(i,j,ll)/dz*ex_volc_SO2(ex)/SECONDS_PER_DAY*1.d9 ! kg s-1
! Water not implemented yet, and it might never be in this branch.
!            H2O_src_3d(i,j,ll)=H2O_src_3d(i,j,ll)+
!     &        daily_z(i,j,ll)/dz*ex_volc_H2O(ex)/SECONDS_PER_DAY*1.d9*
!     &        byMA(ll,i,j) ! kg kg-1 s-1
          enddo

        enddo ! ex
      endif
#endif

#ifdef TRACERS_SPECIAL_Lerner
      if (.not. end_of_day) then
C**** Initialize tables for linoz
        if (itime.ge.itime_tr0(n_O3)) call linoz_setup(n_O3)

C**** Initialize tables for Prather StratChem tracers
        call stratchem_setup
      end if  ! not end of day

C**** Prather StratChem tracers and linoz tables change each month
      IF (modelEclock%getMonth().NE.last_month) THEN
        CALL STRTL  ! one call does all based on n_MPtable_max
        if (itime.ge.itime_tr0(n_O3)) CALL linoz_STRATL
        last_month = modelEclock%getMonth()
      END IF

C**** Tracer specific call for CO2
      call read_CO2_sources(n_CO2)

C**** Tracer specific call for CH4
      call read_CH4_sources(n_CH4)
#endif

#ifdef TRACERS_COSMO
      if (variable_phi .eq. 0) then
         call read_Be_source_noAlpha
         print*, "called old version of Be source"
      end if

      if (variable_phi .eq. 1) then
         call read_Be_source
         print*, "called new version of Be source"
      end if

      if (variable_phi .eq. 2) then
         if ((dayOfYear .eq. 1) .or. (.not. end_of_day)) then
            call update_annual_phi
            print*, "called update_annual_phi"
         end if
      end if

      if (variable_phi .eq. 3) then
         call update_daily_phi
         print*, "called update_daily_phi"
      end if
#endif

!===============================================================================
! Chemistry/OMA/MATRIX/TOMAS case, where surface emissions are of type
! TRACERNAME_XX
!===============================================================================
#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_GASEXCH_GCC)

      !! xday is used by multiple sources below
      xday=dayOfYear

#ifdef TRACERS_SPECIAL_Shindell
C**** Next line for fastj photon fluxes to vary with time:
      if(rad_FL.gt.0) call READ_FL(end_of_day)
C**** Daily tracer-specific calls to read 2D and 3D sources:
      if (COUPLED_CHEM.ne.1) then
        call read_aero(dms_offline,'DMS_FIELD') !not applied directly to tracer
        call read_aero(so2_offline,'SO2_FIELD') !not applied directly to tracer
      endif
#endif /* TRACERS_SPECIAL_Shindell */

#ifdef CUBED_SPHERE
      ! Currently, for the cubed sphere case, the aircraft sources
      ! are applied each time step but only updated daily here:
      ! Thus the logical argument .true. to read from disk and distribute:
      ! All Tracers! (formerly just hardcoded set allowed)
      do n=1,ntm
        if(do_aircraft(n)) then
          call get_aircraft_tracer
     &    (n,trim(trname(n))//'_AIRC',year,xday,daily_gz,.true.,
     &     AIRCstreams(n))
          ! for TOMAS, is trname(n_AECOB(1))=='AECOB_01' ?
        end if
      end do
#endif /* CUBED_SPHERE */

#if defined DYNAMIC_BIOMASS_BURNING && defined ANTHROPOGENIC_FIRE_MODEL
      trans_emis_overr_yr=ABS(o3_yr) ! note: for now, ignores aer_int_yr
      if(trans_emis_overr_yr > 0)then
        xyear=trans_emis_overr_yr
      else
        xyear=year
      endif
      call readflamPopDens(xyear,xday)
#endif

!-------------------------------------------------------------------------------
! tracers loop
!-------------------------------------------------------------------------------
      do n=1,ntm
        pTracer => tracers%getReference(trname(n))

        if ((n>=ntm_chem_beg).and.(n<=ntm_chem_end)) then
          isChemTracer=.true. ! careful: only for this n loop
        else
          isChemTracer=.false.
        end if
!**** Allow overriding of ozone precursor transient emissions date:
! for now, tying this to O3_yr becasue Gavin
! didn't want a new parameter, also not allowing
! day overriding yet, because of that.
#ifdef TRACERS_SPECIAL_Shindell
        if (isChemTracer) then
          trans_emis_overr_yr=ABS(o3_yr)
          if(trans_emis_overr_yr > 0)then
            xyear=trans_emis_overr_yr
          else
            xyear=year
          endif

          select case (trname(n))
          case ('NOx')
            if (NOx_yr > 0) xyear=NOx_yr
          case ('CO')
            if (CO_yr > 0) xyear=CO_yr
          case ('Alkenes', 'Paraffin')
            if (VOC_yr > 0) xyear=VOC_yr
          end select
        else
#endif

! allow overriding of transient aerosol emissions date
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
          xyear=year
          if(aer_int_yr > 0) xyear=aer_int_yr
          select case (trname(n))
          case ('SO2', 'SO4', 'M_ACC_SU', 'M_AKK_SU', 'ASO4__01')
            if (SO2_int_yr > 0) xyear=SO2_int_yr
          case ('NH3')
            if (NH3_int_yr > 0) xyear=NH3_int_yr
          case ('BCII', 'BCB', 'M_BC1_BC', 'M_BOC_BC', 'AECOB_01')
            if (BC_int_yr > 0) xyear=BC_int_yr
          case ('OCII', 'OCB', 'M_OCC_OC', 'M_BOC_OC', 'AOCOB_01',
     &          'vbsAm2', 'vbsAm1', 'vbsAz', 'vbsAp1', 'vbsAp2',
     &          'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
            if (OC_int_yr > 0) xyear=OC_int_yr
          end select
#endif
#ifdef TRACERS_SPECIAL_Shindell
        end if
#endif
#ifdef TRACERS_GASEXCH_GCC
        xyear=year
#endif

! define nread and checkname per tracer
        nread=ntsurfsrc(n)+nBBsources(n)
        if (.not.tracers_amp .and. .not.tracers_tomas) then
          checkname=.false.
        else
          checkname=.true.
        endif

        select case (trname(n))
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        case ('codirect')
          checkname=.false.
          isChemTracer=.true.
          if(trans_emis_overr_yr > 0)then
            xyear=trans_emis_overr_yr
          else
            xyear=year
          endif
#endif

        case ('OCII','M_OCC_OC','SOAgas') ! Koch/AMP/TOMAS cases
          if (.not.tracers_aerosols_soa) nread=nread-1
        case ('SO4','M_AKK_SU','M_ACC_SU',
     &        'ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     &        'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     &        'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15',
     &        'ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05',
     &        'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10',
     &        'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15')
          nread=0
        case ('vbsAm2', 'vbsAm1', 'vbsAz', 'vbsAp1', 'vbsAp2',
     &        'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
          checkname=.false.
        case ('SF6', 'SF6_c', 'nh5', 'nh50', 'e90',
     &        'st8025', 'aoa', 'aoanh', 'tape_rec', 'nh15')
          nread=0 ! regional sources calculated in the code, not via a file
        end select

!-------------------------------------------------------------------------------
! read surface sources of all tracers
!-------------------------------------------------------------------------------
        call readSurfaceSources(pTracer,n,nread,xyear,xday,checkname,
     &                          itime,itime_tr0(n),sfc_src,isChemTracer)
!-------------------------------------------------------------------------------

! post-read calculations
        select case (trname(n))
#ifdef TRACERS_SPECIAL_Shindell
        case ('CH4')
#ifdef WATER_MISC_GRND_CH4_SRC
          do ns=1,ntsurfsrc(n)
            if(pTracer%surfaceSources(ns)%sourceName==
     &         'gsfMGOLjal_src') then
              sfc_src(I_0:I_1,J_0:J_1,n,ns)=
     &          1.698d-12*fearth0(I_0:I_1,J_0:J_1) + ! 5.3558e-5 Jean
     &          5.495d-11*flake0(I_0:I_1,J_0:J_1)  + ! 17.330e-4 Jean
     &          1.141d-12*focean(I_0:I_1,J_0:J_1)    ! 3.5997e-5 Jean
            endif
          end do
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
          if(nread>0) call read_ncep_for_wetlands(end_of_day)
#endif

        case ('N2O5')
          tr3Dsource(I_0:I_1,J_0:J_1,:,:,n) = 0.
          if (COUPLED_CHEM.ne.1)
     &      call read_aero(sulfate,'SULFATE_SA') !not applied directly

       case ('NOx') ! use : for sources, to include BB
         sfc_src(:,J_0:J_1,n,:)=tune_NOx*sfc_src(:,J_0:J_1,n,:)

       case ('Isoprene', 'Terpenes')
         sfc_src(:,J_0:J_1,n,1:ntsurfsrc(n))=
     &     tune_BVOC*sfc_src(:,J_0:J_1,n,1:ntsurfsrc(n))
#endif  /* TRACERS_SPECIAL_Shindell */

        case ('M_OCC_OC', 'OCII')
          if (.not.tracers_aerosols_soa) then
            sfc_src(:,J_0:J_1,n,ntsurfsrc(n):
     &                          ntsurfsrc(n)+nBBsources(n))=
     &      sfc_src(:,J_0:J_1,n,ntsurfsrc(n)-1:
     &                          ntsurfsrc(n)+nBBsources(n)-1)
            sfc_src(:,J_0:J_1,n,ntsurfsrc(n))=0.d0 ! this will become terpene sources
          endif
        end select

      end do ! ntm
!-------------------------------------------------------------------------------
! end tracers loop
!-------------------------------------------------------------------------------

! define ntsurfsrc for sulfate, AFTER reading and AFTER the tracers loop
#ifdef TRACERS_AEROSOLS_Koch
        call set_ntsurfsrc(n_SO4,ntsurfsrc(n_SO2))
#endif
#ifdef TRACERS_AMP
        call set_ntsurfsrc(n_M_ACC_SU, ntsurfsrc(n_SO2))
#ifndef TRACERS_AMP_M4
        call set_ntsurfsrc(n_M_AKK_SU, ntsurfsrc(n_SO2))
#endif
#endif
#ifdef TRACERS_TOMAS
        call set_ntsurfsrc(n_ASO4(1),ntsurfsrc(n_SO2))
#endif

#endif /* TRACERS_SPECIAL_Shindell || TRACERS_AEROSOLS_Koch || TRACERS_AMP || TRACERS_TOMAS */
!===============================================================================
! End of Chemistry/OMA/MATRIX/TOMAS case
!===============================================================================

C**** Initialize tracers here to allow for tracers that 'turn on'
C**** at the start of any day
      call tracer_IC

#ifdef TRACERS_GASEXCH_GCC
      return
#endif
      if ((n_co2n>0).and.end_of_day) call adjust_co2

      return
      end subroutine daily_tracer

!gas exchange CO2 case reset trm here
!for the constCO2 case just reset to atmCO2 which is defined in the rundeck
!for the variable case (presently default) reset to the value scaled by
!the xnow value.
      subroutine adjust_co2
      use resolution, only : lm, psf
      use domain_decomp_atm, only : grid, globalsum, am_i_root
      use tracer_com, only: trm, trmom, n_co2n
      use oldtracer_mod, only: vol2mass
      use geom, only: axyp
      use constant, only: grav
      use model_com, only : nstep=>itime
      implicit none

      integer i,j,l,i_0, i_1, j_0, j_1, n
      real*8 :: sarea,
     &          trm_vert(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                   GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     &          trm_glbavg,factor,atm_glbavg
      real*8 :: get_atmco2

      n=n_CO2n
      i_0=grid%i_strt
      i_1=grid%i_stop
      j_0=grid%j_strt
      j_1=grid%j_stop

         !area weighted tracer global average
      do j=J_0,J_1 ; do i=I_0,I_1
        trm_vert(i,j) = sum(trm(i,j,1:lm,n))
      enddo; enddo

      CALL GLOBALSUM(grid,axyp,    sarea,     all=.true.)
      CALL GLOBALSUM(grid,trm_vert,trm_glbavg,all=.true.)

         !total atm mass
      atm_glbavg = PSF*sarea*100.d0/grav

         !current concentration to new concentration
      factor = get_atmCO2()*atm_glbavg/trm_glbavg *vol2mass(n)*1.d-6
      if(AM_I_ROOT( ))then
        write(*,'(a,i5,8e12.4)')
     .           "TRACER_DRV, factor", nstep,factor,
     .            atm_glbavg,vol2mass(n),
     .            get_atmCO2(),sarea,
     .            PSF,grav,
     .            trm_glbavg/(atm_glbavg*vol2mass(n)*1.d-6)
      endif
#ifdef TRACERS_GASEXCH_CO2_DO_NOT_ADJUST
      return
#endif
      do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
        trm(i,j,l,n) = factor*trm(i,j,l,n)
      enddo; enddo; enddo

      if (factor .lt. 1.d0) then ! adjust moments
        do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
          trmom(:,i,j,l,n)=factor*trmom(:,i,j,l,n)
             !?? do we need trmom=0 for the atmco2 case?
        enddo; end do; end do
      end if

      return
      end subroutine adjust_co2


#ifdef TRACERS_ON
      SUBROUTINE set_tracer_2Dsource
!@sum tracer_source calculates non-interactive sources for tracers
!@vers 2013/03/26
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only: itime,dtsrc,nday
      use SpecialIO_mod, only: write_parallel
      use TracerSurfaceSource_mod, only: TracerSurfaceSource
      use Attributes_mod
      use AttributeDictionary_mod
      use TracerBundle_mod
      use Tracer_mod, only: Tracer
      use model_com, only: modelEclock
      use RESOLUTION, only: im
      USE DOMAIN_DECOMP_ATM, only : GRID, GLOBALSUM,AM_I_ROOT
     *   ,globalmax, getDomainBounds
      USE FLUXES, only: fland,flice,focean,atmsrf
      USE SEAICE_COM, only : si_atm
      USE GHY_COM, only : fearth
      USE CONSTANT, only: tf,bygrav,mair,pi,teeny
      USE LAKES_COM, only : flake
      use OldTracer_mod, only: vol2mass
      use OldTracer_mod, only: trname
      use OldTracer_mod, only: itime_tr0
      use OldTracer_mod, only: do_fire
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR,
     &           SECONDS_PER_HOUR, HOURS_PER_DAY, INT_MONTHS_PER_YEAR
      USE ATM_COM, only: MA  ! Air mass of each box (kg m-2)
      USE TRACER_COM, only: ntm
#ifndef SKIP_TRACER_SRCS
      USE FLUXES, only: trsource
#endif
      use TRACER_COM, only: tracers
      use EmissionRegion_mod, only: numRegions,regions
      use TRACER_COM, only: ef_FACT
      USE RESOLUTION, only : pmtop,psf
      USE GEOM, only: axyp,areag,lat2d_dg,lon2d_dg,imaxj,lat2d
      USE QUSDEF
      USE TRACER_COM, only: sfc_src
      USE TRACER_COM, only: alter_sources
      use TRACER_COM, only: n_isoprene, n_SO2, no_emis_over_ice
      use TRACER_COM, only: trm, ntsurfsrc, rnsrc
      use TRACER_COM, only: tracers, ef_FACT
      use OldTracer_mod, only: itime_tr0, vol2mass, trname
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: n_AH2O, n_AECOB
      use TRACER_COM, only: n_ANUM, n_AECIL, n_AOCIL, n_AOCOB
      use TRACER_COM, only: nbins, n_ASO4, xk
#endif
#if (defined INTERACTIVE_WETLANDS_CH4) && (defined TRACERS_SPECIAL_Shindell)
      USE TRACER_SOURCES, only: ns_wet,add_wet_src
#endif
#ifdef TRACERS_SPECIAL_Lerner
      USE CO2_SOURCES, only: co2_src
      USE CH4_SOURCES, only: ch4_src
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
#ifndef TRACERS_AEROSOLS_SOA
      USE AEROSOL_SOURCES, only: OCT_src
#endif  /* TRACERS_AEROSOLS_SOA */
#endif
#ifdef TRACERS_RADON
      USE AEROSOL_SOURCES, only: rn_src
#endif
#if (defined TRACERS_NITRATE) || (defined TRACERS_AMP) || \
    (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_TOMAS)
      USE apply3d, only : apply_tracer_3Dsource
      USE RAD_COM,  only : cosz1,cosz_day
      use tracer_com, only: seasonalNH3src
#endif
#ifdef TRACERS_AMP
      USE AERO_SETUP, only : RECIP_PART_MASS
      USE TRDIAG_COM, only : taijs=>taijs_loc,ijts_AMPe
#endif
#ifdef TRACERS_TOMAS
      USE TOMAS_EMIS, only : scalesizeSO4,scalesizeCARBO30
#endif
      use TracerHashMap_mod, only:
     &     TracerIterator, operator(/=)
      use Attributes_mod
      use AbstractAttribute_mod
      USE FILEMANAGER, only: openunit,closeunit
      USE Dictionary_mod, only: sync_param
      implicit none
      integer :: i,j,ns,ns_isop,l,ky,n,nsect,kreg
      REAL*8 :: sarea,steppy,base,steppd,x,airm,anngas,
     *  tmon,bydt,tnew,fice
      REAL*8 :: sarea_prt(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                    GRID%J_STRT_HALO:GRID%J_STOP_HALO)
#ifdef TRACERS_SPECIAL_Shindell
c      real*8 :: factj(GRID%J_STRT_HALO:GRID%J_STOP_HALO)
c      real*8 :: nlight, max_COSZ1, fact0
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
!@var src_index source index for the current tracer
!@var src_fact source factor for the current tracer
      integer :: src_index,get_src_index
      real*8 :: src_fact
      interface
        real*8 function get_src_fact(n,OA_not_OC)
          integer, intent(in) :: n
          logical, intent(in), optional :: OA_not_OC
        end function get_src_fact
      end interface
#endif

#ifdef TRACERS_TERP
!@param orvoc_fact Fraction of ORVOC added to Terpenes, for SOA production (Griffin et al., 1999)
      real*8, parameter :: orvoc_fact=0.32d0
      real*8 :: max_isop_flux
#endif  /* TRACERS_TERP */
      integer :: i_ocmip, iu_data
      real*8  :: factor
      real*8  :: trsource_prt(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                        GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      real*8, dimension(NTM) :: trsource_glbavg
!@var ocmip_cfc: CFC-11 emissions estimated from OCMIP surf.conc.
      !60years (1939--1998) OCMIP surfc. concentr. converted to
      !global averaged emission rates
      !each value corresponds to the annual value
!      REAL*8, DIMENSION(:), allocatable, save :: ocmip_cfc
      INTEGER I_0, I_1, J_0, J_1
      class (Tracer), pointer :: pTracer
      integer :: index
      type (TracerSurfaceSource), pointer :: sources(:)
#ifdef TRACERS_TOMAS
      integer :: k, kn
      real*8 :: tot_emis(GRID%I_STRT:GRID%I_STOP,
     &     GRID%J_STRT:GRID%J_STOP)
#endif
      integer :: year, month, dayOfYear,hour,localTimeIndex

      type (TracerIterator) :: iter
      class (AbstractAttribute), pointer :: pa

      call modelEclock%get(year=year, month=month,
     *     dayOfYear=dayOfYear, hour=hour)
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     &     I_STRT=I_0, I_STOP=I_1)

      bydt = 1./DTsrc
#ifdef TRACERS_TOMAS
#ifndef SKIP_TRACER_SRCS
        do k=1,nbins
           trsource(:,J_0:J_1,1,n_ANUM(1)+k-1)=0.
           trsource(:,J_0:J_1,2,n_ANUM(1)+k-1)=0.
           trsource(:,J_0:J_1,3,n_ANUM(1)+k-1)=0.
        enddo
#endif
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
      call calculate_fire_count
#endif
C**** All sources are saved as kg s-1
      iter = tracers%begin()
      do while (iter /= tracers%last())
        pTracer => iter%value()
        pa => pTracer%getReference('index')
        index = pa
        n = index
        pTracer => tracers%getReference(trname(n))
        sources => pTracer%surfaceSources

        if (itime.lt.itime_tr0(n)) then
          call iter%next()
          cycle
        end if

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      src_index=get_src_index(n)
      src_fact=get_src_fact(n)
#endif

      select case (trim(pTracer%getName()))

      case default
!     write(6,*) ' Sources for ',trim(pTracer%getName()),' are not in this routine'
C****
C**** Surface Sources of SF6 and CFCn (Same grid as CFC11)
C****
      case ('SF6','CFC11','CFCn','SF6_c')
#ifndef SKIP_TRACER_SRCS
        trsource(:,:,:,n)=0
#endif
C**** SF6 source increases each year by .3pptv/year
C**** SF6_c source is constant, same as first year SF6, but always
C**** CFCn source increases each year so that the glbavg is from obs
C**** CFC source is the same each year
C**** Distribute source over ice-free land
        steppy = 1./(SECONDS_PER_DAY*INT_DAYS_PER_YEAR)
        if (trim(pTracer%getName()).eq.'SF6' .or.
     *      trim(pTracer%getName()).eq.'CFCn' .or.
     *      trim(pTracer%getName()).eq.'SF6_c') then
C         Make sure index KY=1 in year that tracer turns on
          ky = 1 + (itime-itime_tr0(n))/(nday*INT_DAYS_PER_YEAR)
          if (trim(pTracer%getName()).eq.'SF6_c') ky = 1
          base = (0.3d-12)*vol2mass(n) !pptm
          x = base*ky
          airm = (psf-pmtop)*100.*bygrav*AREAG !(kg/m**2 X m**2 = kg)
          anngas = x*airm
        else if (trim(pTracer%getName()).eq.'CFC11') then
          anngas = 310.d6
        endif

c Could the masks for latlon rectangles be precomputed at
c initialization and their areas saved? (Or do whenever
c fearth changes.)

C**** Source over United States and Canada
        call regional_src(n, .37d0*anngas*steppy,
     &                    -70.d0, -125.d0, 50.d0, 30.d0)
C**** Source over Europe and Russia
        call regional_src(n, .37d0*anngas*steppy,
     &                    45.d0, -10.d0, 65.d0, 36.1d0) ! 0.1 deg offset avoids overlap with Middle East
C**** Source over Far East
        call regional_src(n, .13d0*anngas*steppy,
     &                    150.d0, 120.d0, 45.d0, 20.d0)
C**** Source over Middle East
        call regional_src(n, .05d0*anngas*steppy,
     &                    75.d0, 30.d0, 35.9d0, 15.d0) ! 0.1 deg offset avoids overlap with Europe
C**** Source over South America
        call regional_src(n, .04d0*anngas*steppy,
     &                    -40.d0, -50.d0, -22.5d0, -23.5d0)
C**** Source over South Africa
        call regional_src(n, .02d0*anngas*steppy,
     &                    30.d0, 25.d0, -24.d0, -28.d0)
C**** Source over Australia and New Zealand
        call regional_src(n, .02d0*anngas*steppy,
     &                    150.5d0, 149.5d0, -33.5d0, -34.5d0)

        if (trim(pTracer%getName()).eq.'CFCn') then
          !print out global average for each time step before weighing
          !in the OCMIP values
          sarea  = 0.
          trsource_glbavg(n)=0.
          sarea_prt(:,:)  = 0.
          trsource_prt(:,:) = 0.
          do j=J_0,J_1
            do i=I_0,I_1
             factor = axyp(i,j)*fearth(i,j)
             sarea_prt(i,j)= FACTOR
#ifndef SKIP_TRACER_SRCS
             trsource_prt(i,j) = trsource(i,j,1,n)*FACTOR
#endif
            enddo
          enddo

          CALL GLOBALSUM(grid, sarea_prt,    sarea, all=.true.)
          CALL GLOBALSUM(grid,trsource_prt,trsource_glbavg(n),
     &                                                all=.true.)

          trsource_glbavg(n)=trsource_glbavg(n)/sarea

          !weight trsource by ocmip_cfc global average
          !number of steps/year=INT_DAYS_PER_YEAR*SECONDS_PER_DAY/dtsrc
          !                    =365*86400/1800 =17520
!         if (.not.allocated(ocmip_cfc)) then
!           !read in OCMIP based CFC-11 global emissions
!           !=sum(dC/dt) for each hemisphere
!           !these are *annual global averages* and need to be
!           !converted to our timestep value
!           allocate(ocmip_cfc(67))
!           print*, 'opening file=OCMIP_cfc.dat'
!           call openunit('OCMIP_cfc',iu_data,.false.,.true.)
!           do i=1,67
!             read(iu_data,'(5x,e12.4)')ocmip_cfc(i)
!           enddo
!           call closeunit(iu_data)
!         endif
!         i_ocmip=(itime-itime_tr0(n))/INT_DAYS_PER_YEAR/
!    &            int(SECONDS_PER_DAY/dtsrc)+1
!         if (mod(itime,INT_DAYS_PER_YEAR*int(SECONDS_PER_DAY/dtsrc))
!    &        .eq. 0.) then
!           write(6,'(a,2i5)'),'TRACERS_DRV, new year: itime, i_ocmip=',
!    &                         itime,i_ocmip
!         endif
!#ifndef SKIP_TRACER_SRCS
!          do j=J_0,J_1 ! TNL
!            do i=1,72
!               trsource(i,j,1,n) = trsource(i,j,1,n)*
!     &           (ocmip_cfc(i_ocmip)/(INT_DAYS_PER_YEAR*
!     &           SECONDS_PER_DAY/dtsrc)) / trsource_glbavg(n)
!            enddo
!          enddo
!#endif

          !recompute global average after weighting in OCMIP
          sarea  = 0.
          trsource_glbavg(n)=0.
          sarea_prt(:,:)  = 0.
          trsource_prt(:,:) = 0.
          do j=J_0,J_1
            do i=I_0,I_1
             factor = axyp(i,j)*fearth(i,j)
             sarea_prt(i,j)= FACTOR
#ifndef SKIP_TRACER_SRCS
             trsource_prt(i,j) = trsource(i,j,1,n)*FACTOR
#endif
            enddo
          enddo

          CALL GLOBALSUM(grid, sarea_prt,    sarea, all=.true.)
          CALL GLOBALSUM(grid,trsource_prt,trsource_glbavg(n),
     &                                                 all=.true.)

          trsource_glbavg(n)=trsource_glbavg(n)/sarea

        endif  ! pTracer==CFCn

#ifdef TRACERS_PASSIVE

       case ('nh5','nh50','nh15')

        trsource(:,:,:,n)=0

       case ('e90')

        trsource(:,:,:,n)=0

#endif

C****
C**** Surface Sources for Radon-222
C****
      case ('Rn222')
#ifndef SKIP_TRACER_SRCS
        trsource(:,J_0:J_1,:,n)=0
C**** ground source
        steppd = 1./SECONDS_PER_DAY
        do j=J_0,J_1
          do i=I_0,I_1
          if (rnsrc.eq.0) then !standard source
C**** source from ice-free land
            if(atmsrf%tsavg(i,j).lt.tf) then !composite surface air temperature
              trsource(i,j,1,n) = 1.0d-16*steppd*axyp(i,j)*fearth(i,j)
            else  ! 1 atom/cm^2/s
              trsource(i,j,1,n) = 3.2d-16*steppd*axyp(i,j)*fearth(i,j)
            end if
         else if (rnsrc.eq.1) then !Conen and Robertson
              trsource(i,j,1,n) = 3.2d-16*steppd*axyp(i,j)*fearth(i,j)
c add code to implement Conen and Robertson - linear decrease in Rn222
c   emission from 1 at 30N to 0.2 at 70N and 0.2 north of 70N
           if (nint(lat2d_dg(i,j)).gt.30 .and.
     &         nint(lat2d_dg(i,j)).lt.70) then
             trsource(i,j,1,n)=trsource(i,j,1,n)*
     &            (1.d0-(lat2d_dg(i,j)-30.d0)/40.d0*0.8d0)
           else if (nint(lat2d_dg(i,j)).ge.70) then
             trsource(i,j,1,n)=0.2*trsource(i,j,1,n)
           endif
          else if (rnsrc.eq.2) then !Schery and Wasiolek
#ifdef TRACERS_RADON
c Schery source
          trsource(i,j,1,n)=rn_src(i,j,month)
#endif
          endif
          if (rnsrc.le.1) then
C**** source from ice-free ocean
            trsource(i,j,1,n) =trsource(i,j,1,n)+ 1.6d-18*steppd*axyp(i
     *           ,j)*(1.-fland(i,j))*(1.-si_atm%rsi(i,j))
          endif
          enddo                 !i
        enddo                   !j
#endif

#ifdef TRACERS_SPECIAL_Lerner
C****
C**** Sources and sinks for CO2 (kg s-1)
C****
      case ('CO2')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
            trsource(:,j,ns,n) = co2_src(:,j,ns)*axyp(:,j)
          end do
        end do

C****
C**** Sources and sinks for CH4 (kg s-1)
C****
      case ('CH4')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
            trsource(:,j,ns,n) = ch4_src(:,j,ns)*axyp(:,j)
          end do
        end do
C****
C**** Sources and sinks for N2O:
C**** First layer is set to a constant 462.2 ppbm. (300 PPB V)
C****
      case ('N2O')
      do j=J_0,J_1
        trsource(:,j,1,n) = (MA(1,:,j)*axyp(:,j)*462.2d-9
     *   -trm(:,j,1,n))*bydt
      end do
C****
C**** Linoz Deposition from layer 1
C****
      case ('O3')
      call linoz_depo(1,n)
#endif
C****
C**** Sources and sinks for 14CO2
C**** NOTE: This tracer is supposed to start on 10/16
C**** Decay is a function of the number of months since itime_tr0
C**** The tracer is reset to specific values in layer 1 only if
C****   this results in a sink
C****
      case ('14CO2')
#ifndef SKIP_TRACER_SRCS
      tmon = (itime-itime_tr0(n))*INT_MONTHS_PER_YEAR/
     &       (nday*INT_DAYS_PER_YEAR)
      trsource(:,J_0:J_1,1,n) = 0.
      do j=J_0,J_1
      do i=I_0,I_1
         if (lat2d(i,j).lt.0.) then
               tnew = MA(1,i,j)*axyp(i,j)*(4.82d-18*46./mair)*
     *          (44.5 + tmon*(1.02535d0 - tmon*
     *                  (2.13565d-2 - tmon*8.61853d-5)))
               if (tnew.lt.trm(i,j,1,n))
     *             trsource(i,j,1,n) = (tnew-trm(i,j,1,n))*bydt
         else
               tnew = MA(1,i,j)*axyp(i,j)*(4.82d-18*46./mair)*
     *          (73.0 - tmon*(0.27823d0 + tmon*
     *                  (3.45648d-3 - tmon*4.21159d-5)))
               if (tnew.lt.trm(i,j,1,n))
     *             trsource(i,j,1,n) = (tnew-trm(i,j,1,n))*bydt
         endif
      end do
      end do
#endif

C****
C**** No non-interactive surface sources of Water
C****
      case ('Water')
#ifndef SKIP_TRACER_SRCS
        trsource(:,J_0:J_1,:,n)=0.d0
#endif

#ifdef TRACERS_SPECIAL_Shindell
      case ('Ox','NOx','ClOx','BrOx','N2O5','HNO3','H2O2','CH3OOH',
     &      'HCHO','HO2NO2','CO','PAN','AlkylNit','Alkenes','Paraffin',
#ifdef TRACERS_dCO
     *      'd13Calke','d13CPAR',
     *      'd17OPAN','d18OPAN','d13CPAN',
     *      'dMe17OOH', 'dMe18OOH', 'd13MeOOH',
     *      'dHCH17O', 'dHCH18O', 'dH13CHO',
     *      'dC17O', 'dC18O', 'd13CO',
#endif  /* TRACERS_dCO */
     &      'HCl','HOCl','ClONO2','HBr','HOBr','BrONO2','N2O','CFC',
     &      'stratOx','codirect')
#ifdef DYNAMIC_BIOMASS_BURNING
        if(do_fire(n))call dynamic_biomass_burning(n,ntsurfsrc(n)+1)
#endif
        do ns=1,ntsurfsrc(n); do j=J_0,J_1
          trsource(I_0:I_1,j,ns,n)=
     &    sfc_src(I_0:I_1,j,n,ns)*axyp(I_0:I_1,j)
        end do ; end do
#ifdef TRACERS_TERP
      case ('Terpenes')
        do ns=1,ntsurfsrc(n)
          if(ns==1) then
            do j=J_0,J_1
              trsource(I_0:I_1,j,ns,n)=
     &        sfc_src(I_0:I_1,j,n,ns)*axyp(I_0:I_1,j)
            end do
! If no orvoc file provided, scale up the terpenes one instead.
! 0.4371 is the ratio of orvoc/isoprene emissions in the Lathiere et al. (2005) results
            if(ntsurfsrc(n)==1) then ! no orvoc file exists
              call globalmax(grid,
     &             maxval(sfc_src(I_0:I_1,J_0:J_1,n_Isoprene,
     &                            1:ntsurfsrc(n_Isoprene))),
     &             max_isop_flux)
              if (max_isop_flux <= 0.d0) call stop_model(
     &          'Offline isoprene sources are needed', 255)
              do ns_isop=1,ntsurfsrc(n_Isoprene) ! use all Isoprene sources for orvoc scaling
                do j=J_0,J_1
                  trsource(I_0:I_1,j,ns,n)=trsource(I_0:I_1,j,ns,n)+
     &            orvoc_fact*0.4371*axyp(I_0:I_1,j)*
     &            sfc_src(I_0:I_1,j,n_Isoprene,ns_isop)
                end do
              end do
            end if
          else ! use the orvoc file
            do j=J_0,J_1
              trsource(I_0:I_1,j,ns,n)=orvoc_fact*
     &        sfc_src(I_0:I_1,j,n,ns)*axyp(I_0:I_1,j)
            end do
          endif
        end do
#endif  /* TRACERS_TERP */
      case ('CH4')
#ifdef DYNAMIC_BIOMASS_BURNING
        if(do_fire(n))call dynamic_biomass_burning(n,ntsurfsrc(n)+1)
#endif
        do ns=1,ntsurfsrc(n); do j=J_0,J_1
          trsource(I_0:I_1,j,ns,n)=
     &    sfc_src(I_0:I_1,j,n,ns)*axyp(I_0:I_1,j)
        end do ; end do
#ifdef INTERACTIVE_WETLANDS_CH4
        if(ntsurfsrc(n) > 0) then
          call alter_wetlands_source(n,ns_wet)
          do j=J_0,J_1
            trsource(I_0:I_1,j,ns_wet,n)=trsource(I_0:I_1,j,ns_wet,n)+
     &      add_wet_src(I_0:I_1,j)*axyp(I_0:I_1,j)
          enddo
        endif
#endif
#if !defined(PS_BVOC) && !defined(BIOGENIC_EMISSIONS)
      case ('Isoprene')
! Isoprene sources to be emitted only during sunlight, and
! weighted by cos of solar zenith angle:
        do ns=1,ntsurfsrc(n); do j=J_0,J_1; do i=I_0,I_1
          if(COSZ1(i,j)>0.)then
            trsource(i,j,ns,n)=(COSZ1(i,j)/(COSZ_day(i,j)+teeny))*
     &      sfc_src(i,j,n,ns)*axyp(i,j)
          else
            trsource(i,j,ns,n)=0.d0
          endif
        end do ; end do; enddo
#endif
#endif /* TRACERS_SPECIAL_Shindell */

#ifdef TRACERS_AEROSOLS_OCEAN
      case ('OCocean')
        call read_seawifs_chla(month) ! CHECK this has to be called once per month, not every timestep
#endif  /* TRACERS_AEROSOLS_OCEAN */

#ifndef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_TOMAS
        case ('SOAgas')
!OCT_src is kg/month? or kg/sec??
        do j=J_0,J_1; do i=I_0,I_1
           trsource(i,j,ntsurfsrc(n),n)=OCT_src(i,j,month)
         end do; enddo
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      case ('SO2', 'SO4', 'M_ACC_SU', 'M_AKK_SU',
     &      'BCII', 'BCB', 'OCII', 'OCB',
     &      'vbsAm2', 'vbsAm1', 'vbsAz', 'vbsAp1', 'vbsAp2',
     &      'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6',
     &      'M_BC1_BC', 'M_OCC_OC', 'M_BOC_BC', 'M_BOC_OC',
     &      'ASO4__01','AOCOB_01','AECOB_01')

#ifdef DYNAMIC_BIOMASS_BURNING
        if(do_fire(n))call dynamic_biomass_burning(n,ntsurfsrc(n)+1)
#endif

#ifndef TRACERS_AEROSOLS_SOA
        select case (trim(pTracer%getName()))
        case ('OCII', 'M_OCC_OC')
          sfc_src(:,J_0:J_1,src_index,ntsurfsrc(n))=
     &      OCT_src(:,J_0:J_1,month)/axyp(:,J_0:J_1)/src_fact
        end select
#endif  /* TRACERS_AEROSOLS_SOA */

        do ns=1,ntsurfsrc(src_index)
          trsource(:,J_0:J_1,ns,n)=
     &      sfc_src(:,J_0:J_1,src_index,ns)
     &      *axyp(:,J_0:J_1)*src_fact

#ifdef TRACERS_TOMAS
!ntsurfsrc(3) for number

!     ns=1 : SO4 number
!     ns=2 : EC number
!     ns=3 : OC number
        tot_emis(:,J_0:J_1)=0.0
          if(n.eq.n_ASO4(1))then
             tot_emis(:,J_0:J_1)= trsource(:,J_0:J_1,ns,n_ASO4(1))

             do k=1,nbins
                trsource(:,J_0:J_1,ns,n_ASO4(1)+k-1)=
     &              tot_emis(:,J_0:J_1)*scalesizeSO4(k)

                trsource(:,J_0:J_1,1,n_ANUM(1)+k-1)=
     &           trsource(:,J_0:J_1,1,n_ANUM(1)+k-1) +
     &               trsource(:,J_0:J_1,ns,n_ASO4(1)+k-1)
     &               /sqrt(xk(k)*xk(k+1))
              enddo


          elseif(n.eq.n_AECOB(1))then

             tot_emis(:,J_0:J_1)= trsource(:,J_0:J_1,ns,n_AECOB(1))

             do k=1,nbins
                trsource(:,J_0:J_1,ns,n_AECOB(1)+k-1)=
     &               tot_emis(:,J_0:J_1)*scalesizeCARBO30(k)*0.8

                trsource(:,J_0:J_1,ns,n_AECIL(1)+k-1)=
     &               tot_emis(:,J_0:J_1)*scalesizeCARBO30(k)*0.2

                trsource(:,J_0:J_1,2,n_ANUM(1)+k-1)=
     &           trsource(:,J_0:J_1,2,n_ANUM(1)+k-1) +
     &             ( trsource(:,J_0:J_1,ns,n_AECOB(1)+k-1)+
     &                 trsource(:,J_0:J_1,ns,n_AECIL(1)+k-1))
     &             /sqrt(xk(k)*xk(k+1))
             enddo
          elseif(n.eq.n_AOCOB(1))then

             tot_emis(:,J_0:J_1)= trsource(:,J_0:J_1,ns,n_AOCOB(1))

             do k=1,nbins
                trsource(:,J_0:J_1,ns,n_AOCOB(1)+k-1)=
     &               tot_emis(:,J_0:J_1)*scalesizeCARBO30(k)*0.5

                trsource(:,J_0:J_1,ns,n_AOCIL(1)+k-1)=
     &               tot_emis(:,J_0:J_1)*scalesizeCARBO30(k)*0.5

                trsource(:,J_0:J_1,3,n_ANUM(1)+k-1)=
     &           trsource(:,J_0:J_1,3,n_ANUM(1)+k-1) +
     &              ( trsource(:,J_0:J_1,ns,n_AOCOB(1)+k-1)+
     &                trsource(:,J_0:J_1,ns,n_AOCIL(1)+k-1))
     &            /sqrt(xk(k)*xk(k+1))
             enddo
          endif

#endif
        enddo ! ns
#if (defined TRACERS_NITRATE) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      case ('NH3')
#ifdef DYNAMIC_BIOMASS_BURNING
        if(do_fire(n))call dynamic_biomass_burning(n,ntsurfsrc(n)+1)
#endif
        do ns=1,ntsurfsrc(n)
          if (ns == seasonalNH3src) then
! add annual cycle to agricultural emissions
            do j=J_0,J_1; do i=I_0,I_1
              if (cosz1(i,j) > 0.) then
                trsource(i,j,ns,n)=sfc_src(i,j,n,ns)
     &            *axyp(i,j)*cosz1(i,j)*4.d0
              else
                trsource(i,j,ns,n)=0.d0
              endif
            enddo; enddo
          else
            trsource(:,J_0:J_1,ns,n)=sfc_src(:,J_0:J_1,n,ns)
     &        *axyp(:,J_0:J_1)
          endif
        enddo

#endif /* TRACERS_NITRATE || TRACERS_AMP || TRACERS_TOMAS */
#endif /* (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) || (defined TRACERS_TOMAS) || (defined TRACERS_GASEXCH_GCC) */

! Addition of surface emission sources for GCC
#ifdef TRACERS_GASEXCH_GCC
      case ('CO2n')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
#ifdef GCC_ZERO_EMISSION
            trsource(I_0:I_1,j,ns,n)=sfc_src(I_0:I_1,j,n,ns)
     &                              *axyp(I_0:I_1,j) * 0.0
#else
            trsource(I_0:I_1,j,ns,n)=sfc_src(I_0:I_1,j,n,ns)
     &                              *axyp(I_0:I_1,j)
#endif
          enddo
        enddo
#endif
      end select

! please keep at end of tracer loop :
! TODO: should be able to uncomment
!       this and delete the subsequent block when F2003 compilers are ready.
! Note 1: Whoever implements, please check that the "nn" index in the commented call
!         is correct.
! Note 2: When the diurnal cycle functionality was added the below-commented
!         routine was not updated.
c$$$      do ns = 1, size(sources)      ! loop over source
c$$$        call emissionScenario%scaleSource(trsource(:,:,ns,n),
c$$$     &       sources(ns)%trsect_index(1:sources(nn)%num_tr_sectors))
c$$$      end do

#ifndef SKIP_TRACER_SRCS
      ! First regional sector alterations:
      if(alter_sources)then                     ! if any region/sector altering requested
        do ns=1,ntsurfsrc(n)                    ! loop over sources
          do nsect=1,sources(ns)%num_tr_sectors ! and sectors for that source
            do j=J_0,J_1                        ! and horizonal space
              do i=I_0,imaxj(j)
                do kreg=1,numRegions            ! loop defined regions
                  if(
     &            lat2d_dg(i,j)>=regions(kreg)%southernEdge .and.  ! check if
     &            lat2d_dg(i,j)<=regions(kreg)%northernEdge .and.  ! currently
     &            lon2d_dg(i,j)>=regions(kreg)%westernEdge .and.   ! in region
     &            lon2d_dg(i,j)< regions(kreg)%easternEdge) then ! change to <= after thinking about it.
       if(ef_fact(sources(ns)%tr_sect_index(nsect),kreg) > -1.e20)then
         trsource(i,j,ns,n)=trsource(i,j,ns,n)*
     &   ef_fact(sources(ns)%tr_sect_index(nsect),kreg)
       end if
                  end if ! in-region
                end do   ! regions
              end do     ! i
            end do       ! j
          end do         ! sector
        end do           ! sources
      end if             ! any region/sector altering of sources requested

      ! Then diurnal cycle application:
      do ns=1,ntsurfsrc(n)                    ! loop over sources
        if(sources(ns)%applyDiurnalCycle)then ! does source have a diurnal cycle defined?
          ! this might not work if there aren't an equal number of timesteps each hour:
          if(MOD(SECONDS_PER_HOUR,dtsrc).ne.0.d0)then
            call write_parallel('Diurnal emissions steps/hr problem')
            call stop_model('Problem w/ emissions diurnal cycle 2',255)
          end if
          do j=J_0,J_1                        ! loop horizontal space
            do i=I_0,imaxj(j)
              ! intendinf here for localTimeIndex an integer index ranging from 1 to INT_HOURS_PER_DAY
              localTimeIndex=(hour+1)
     &            +NINT((i-(IM+1)/2.)*HOURS_PER_DAY/float(IM))
              if(localTimeIndex>HOURS_PER_DAY)
     &            localTimeIndex=localTimeIndex-HOURS_PER_DAY
              if(localTimeIndex<1)
     &            localTimeIndex=localTimeIndex+HOURS_PER_DAY
              trsource(i,j,ns,n)=trsource(i,j,ns,n)*
     &            sources(ns)%diurnalCycle(localTimeIndex)
            end do     ! i
          end do       ! j
        end if         ! this source has a diurnal cycle
      end do           ! sources

      ! Optionally set sources to zero over (>90%) ice:
      if(no_emis_over_ice > 0)then
        do j=J_0,J_1
          do i=I_0,imaxj(j)
            fice=flice(i,j)+si_atm%rsi(i,j)*(focean(i,j)+flake(i,j))
            if(fice > 0.9d0) trsource(i,j,:,n)=0.d0
          end do
        end do
      end if
#endif

      call iter%next()
      end do ! n - main tracer loop

#if defined(DYNAMIC_BIOMASS_BURNING) && (defined DETAILED_FIRE_OUTPUT)
      call accumulateVegTypesDiag
#endif

      END SUBROUTINE set_tracer_2Dsource

      subroutine regional_src(n,source,lon_e,lon_w,lat_n,lat_s)
!@sum Assign regional 2d sources
!@auth Kostas Tsigaridis, based on old Lerner code
        use DOMAIN_DECOMP_ATM, only : globalsum,grid,getDomainBounds
        use GEOM, only: axyp
        use GHY_COM, only : fearth
#ifndef SKIP_TRACER_SRCS
        use FLUXES, only: trsource
#endif
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: source, lon_e, lon_w, lat_n, lat_s
        real*8 :: sarea_prt(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                      grid%J_STRT_HALO:grid%J_STOP_HALO)
        real*8 :: sarea
        integer :: i_0,i_1,j_0,j_1
        integer :: i,j

        call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1,
     &                             J_STRT=J_0, J_STOP=J_1)

        call get_latlon_mask(lon_w,lon_e,lat_s,lat_n,sarea_prt)
        do j=j_0,j_1; do i=i_0,i_1
          sarea_prt(i,j) = sarea_prt(i,j)*axyp(i,j)*fearth(i,j)
        enddo; enddo
        call globalsum(grid, sarea_prt, sarea, all=.true.)
#ifndef SKIP_TRACER_SRCS
        do j=j_0,j_1; do i=i_0,i_1
            trsource(i,j,1,n) = trsource(i,j,1,n) +
     &         source*sarea_prt(i,j)/sarea
        enddo; enddo
#endif
      end subroutine regional_src

      subroutine get_latlon_mask(lon_w,lon_e,lat_s,lat_n,latlon_mask)
!@sum Set mask array to 1 for all cells overlapping a lat-lon rectangle
!@auth Kelley
      use domain_decomp_atm, only : getDomainBounds,grid
#ifdef CUBED_SPHERE
      use geom, only : lon2d_dg,lat2d_dg
#else
      use geom, only : lon_to_i,lat_to_j
#endif
      implicit none
      real*8 :: lon_w,lon_e,lat_s,lat_n
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     latlon_mask
      integer :: i,j, i_0,i_1,j_0,j_1
      integer :: ie,iw,js,jn

      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1,
     &     J_STRT=J_0, J_STOP=J_1)

      latlon_mask(:,:) = 0d0

#ifdef CUBED_SPHERE
      do j=j_0,j_1
      do i=i_0,i_1
        if(lon2d_dg(i,j) >= lon_w .and.
     &     lon2d_dg(i,j) <= lon_e .and.
     &     lat2d_dg(i,j) >= lat_s .and.
     &     lat2d_dg(i,j) <= lat_n) then
          latlon_mask(i,j) = 1d0
        endif
      enddo
      enddo
c The above approach to defining masks does not work for latlon
c rectangles whose lon and/or lat lengths are less than the
c gridsize.  Solution: if a rectangle is narrow in either
c dimension, trace a discretized path from its SW to its NE
c corner and call lonlat_to_ij for each latlon point along the
c path, setting the mask to 1 for the returned i,j that are in
c the local domain.  Or, trace all 4 edges of the rectangle.
c      if(narrow rectangle) then
c        lon = lon_w
c        lat = lat_s
c        do point=1,npoints_traverse
c          lon = lon + (lon_e-lon_w)/npoints_traverse
c          lat = lat + (lat_n-lat_s)/npoints_traverse
c          call lonlat_to_ij((/lon,lat/),ij)
c          if(i,j in local domain) then
c            latlon_mask(i,j) = 1d0
c          endif
c        enddo
c      endif
#else
c latlon grid
      ie = lon_to_i(lon_e)
      iw = lon_to_i(lon_w)
      jn = lat_to_j(lat_n)
      js = lat_to_j(lat_s)
      do j=max(js,j_0),min(jn,j_1)
        latlon_mask(iw:ie,j) = 1d0
      enddo
#endif
      return
      end subroutine get_latlon_mask

      SUBROUTINE tracer_3Dsource
!@sum tracer_3Dsource calculates interactive sources for tracers
!@+   Please note that if the generic routine 'apply_tracer_3Dsource'
!@+   is used, all diagnostics and moments are updated automatically.
!@vers 2013/03/27
!@auth Jean Lerner/Greg Faluvegi
!@calls DIAGTCA, masterchem, apply_tracer_3Dsource
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds,
     & write_parallel,AM_I_ROOT
      use RESOLUTION, only: LM
c$$$      use OldTracer_mod, only: itime_tr0, do_fire, trname, do_aircraft
c$$$      use OldTracer_mod, only: tr_mm, nBBsources, mass2vol
      use OldTracer_mod
#ifdef TRACERS_AEROSOLS_VBS
      USE TRACERS_VBS, only: vbs_tr
#endif  /* TRACERS_AEROSOLS_VBS */
      USE TRACER_COM, only: ntm, sfc_src, trm
      use TRACER_COM, only: mchem, mtrace, n_BCIA, n_BCII, n_CFC, n_CH4
      use TRACER_COM, only: n_DMS, n_H2O2_s, n_HNO3, n_MSA, N_N2O
      use TRACER_COM, only: n_N_d1, n_N_d2, n_N_d3, n_NH3, n_NH4
      use TRACER_COM, only: n_NOx, n_NO3p, n_OCIA, n_OCII
      use TRACER_COM, only: n_SO4, n_SO4_d1, n_SO4_d2, n_SO4_d3
      use TRACER_COM, only: n_SO2
      use TRACER_COM, only: ntsurfsrc
      use TRACER_COM, only: ntm_chem_beg, ntm_chem_end
      use TRACER_COM, only: n_NOx, naircraft, nBiomass, nChemistry
      use TRACER_COM, only: nVolcanic, nOverwrite, nChemloss, nOther
      use TRACER_COM, only:  trans_emis_overr_day, trans_emis_overr_yr
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: nbins, n_AH2O
      use TRACER_COM, only: n_AOCIL, n_ANUM, n_ANACL, n_ADUST
      use TRACER_COM, only: n_AECOB, n_AOCOB, n_ASO4, n_H2SO4, n_SOAGAS
      use TRACER_COM, only: nChemistry, n_AECIL, ntm_tomas
#endif
#ifdef TRACERS_AMP
      use TRACER_COM, only: n_H2SO4,n_M_BC1_BC,n_M_NO3,n_M_NH4
      use TRACER_COM, only: ntmAMPi, ntmAMPe
#endif
#ifdef SHINDELL_STRAT_EXTRA
      use TRACER_COM, only: n_GLT, n_stratOx
#endif
      use TRACER_COM, only: tune_BBsources
      use TRACER_COM, only: n_aoa, n_aoanh
      USE CONSTANT, only : mair, byavog, pi
#ifndef SKIP_TRACER_SRCS
      USE FLUXES, only: tr3Dsource
#endif
      use model_com, only: modelEclock
      USE MODEL_COM, only: itime,dtsrc,itimeI
      USE ATM_COM, only: MA,byMA ! Air mass of each box (kg m-2)
      use ATM_COM, only: phi
      USE ATM_COM, only: pmid,pmidl00
      USE apply3d, only : apply_tracer_3Dsource
      USE GEOM, only : byaxyp,axyp
      USE RAD_COM, only: o3_yr
      USE Dictionary_mod, only : get_param, is_set_param
      use trdiag_com, only : taijls=>taijls_loc,ijlt_prodSO4gs
#if (defined TRACERS_COSMO)
      USE COSMO_SOURCES, only: be7_src_3d, be10_src_3d, be7_src_param
#endif
#ifdef TRACERS_AEROSOLS_SOA
      USE TRACERS_SOA, only: n_soa_i,n_soa_e
#endif  /* TRACERS_AEROSOLS_SOA */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      USE AEROSOL_SOURCES, only: so2_src_3d
#endif
      USE PBLCOM, only: dclev
#ifdef TRACERS_AMP
      USE AERO_SETUP, only : RECIP_PART_MASS
      USE TRDIAG_COM, only : itcon_AMP, itcon_AMPe,itcon_AMPm
      USE TRDIAG_COM, only : taijs=>taijs_loc,ijts_AMPe
#endif
#ifdef TRACERS_TOMAS
      USE CONSTANT, only : pi
      USE TRDIAG_COM, only : itcon_TOMAS,itcon_subcoag
      USE TRDIAG_COM, only : taijs=>taijs_loc,ijts_TOMAS
      USE TOMAS_AEROSOL, only : TRM_EMIS,xk,icomp,idiag
      USE TOMAS_EMIS, only : scalesizeCARBO100,
     &     scalesizeCARBO30, scalesizeSO4,
     &     scalesizeSO4_vol,scalesizeSO4_bio
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: fix_CH4_chemistry
#endif

#ifdef TRACERS_SPECIAL_Shindell
      use RAD_COM, only: rad_to_chem
      use TRCHEM_Shindell_COM, only: fact_cfc,
     &     use_rad_n2o, use_rad_ch4, use_rad_cfc, topLevelOfChemistry
#endif
#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS) ||\
    (defined TRACERS_GASEXCH_GCC)

      use TRACER_COM, only: AIRCstreams
#endif
      USE GEOM, only: lat2d_dg

      implicit none
      INTEGER n,ns,najl,i,j,l,blay,xday   ; real*8 now
      REAL*8 factor
      INTEGER J_0, J_1, I_0, I_1
!@var src_index source index for the current tracer
!@var src_fact source factor for the current tracer
      integer :: src_index,get_src_index,bb_i,bb_e
      real*8 :: src_fact
      real*8 :: bydt

      interface
        real*8 function get_src_fact(n,OA_not_OC)
          integer, intent(in) :: n
          logical, intent(in), optional :: OA_not_OC
        end function get_src_fact
      end interface
      integer :: initial_ghg_setup
!@var blsrc (m2/s) tr3Dsource (kg s-1) in boundary layer,
!@+                per unit of air mass (kg/m2)
      real*8 :: blsrc
#ifdef CUBED_SPHERE
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
     &     :: dummy3d
#endif
#ifdef TRACERS_TOMAS
      integer :: k, kk,kn,jc,tracnum
      real*8, dimension (GRID%I_STRT:GRID%I_STOP,
     &     GRID%J_STRT:GRID%J_STOP,LM,NBINS) :: TOMAS_bio,TOMAS_air
#endif
      integer :: year, dayOfYear

      call modelEclock%get(year=year, dayOfYear=dayOfYear)
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      bydt = 1./DTsrc

C**** All sources are saved as kg s-1
      do n=1,NTM
      if (itime.lt.itime_tr0(n)) cycle
      src_index=get_src_index(n)
      src_fact=get_src_fact(n)

      select case (trname(n))

      case default
#ifdef TRACERS_SPECIAL_Lerner
C****
      case ('CH4')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Trop_chem_CH4(1,n)
      call apply_tracer_3Dsource(1,n)
      call Strat_chem_Prather(2,n)
      call apply_tracer_3Dsource(2,n,.false.)
C****
      case ('O3')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Trop_chem_O3(2,3,n)
        call apply_tracer_3Dsource(2,n,.false.)
        call apply_tracer_3Dsource(3,n,.false.)
      call Strat_chem_O3(1,n)
        call apply_tracer_3Dsource(1,n,.false.)
C****
      case ('N2O','CFC11')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Strat_chem_Prather(nChemistry,n)
      call apply_tracer_3Dsource(nChemistry,n,.FALSE.)

#endif

#ifdef TRACERS_PASSIVE

C**** NH5, NH50 and NH15: idealized loss tracers set over NH midlatitudes (30N-50N)

      case ('nh5','nh50','nh15')

      tr3Dsource(:,:,:,:,n) = 0.

       do l=1,lm
        do j=J_0,J_1
          do i=I_0,I_1
           if (l==1) then !enforce concentrations in surface layer
            if (nint(lat2d_dg(i,j)).gt.30 .and.
     &         nint(lat2d_dg(i,j)).lt.50) then
      tr3Dsource(i,j,l,:,n) = (MA(l,i,j)*axyp(i,j)*100.0d-9
     *        -trm(i,j,l,n))*bydt
            endif
           else
              tr3Dsource(i,j,l,:,n) = 0.
           endif
          enddo
         enddo
        enddo

      call apply_tracer_3Dsource(nChemistry,n,.FALSE.)

C**** E90: idealized loss tracer set over entire surface layer

      case ('e90')

      tr3Dsource(:,:,:,:,n) = 0.

       do l=1,lm
        do j=J_0,J_1
          do i=I_0,I_1
           if (l==1) then
      tr3Dsource(i,j,l,:,n) = (MA(l,i,j)*axyp(i,j)*100.0d-9
     *        -trm(i,j,l,n))*bydt
           else
              tr3Dsource(i,j,l,:,n) = 0.
           endif
          enddo
         enddo
        enddo

      call apply_tracer_3Dsource(nChemistry,n,.FALSE.)

C****AOANH and AOA: Two mean age tracers are defined, one with respect to the NH midlatitude surface and one with respect to the Earth's surface.  In lieu of the "clock-tracer" implementation (see GLT tracer), here we solve for the mean age as the solution to d(G)dt=1, where G is the mean age and d/dt is the advective derivative. Zero boundary conditions are enforced over the source region (either NH midlatitude surface layer (30N-50N) or the entire Earth's surface). Units are in days.

      case ('aoanh')

       tr3Dsource(:,:,:,:,n) = 0.

       factor = DTsrc/(60.d0*60.d0*24.d0) !ensure units of days

       do l=1,lm
        do j=J_0,J_1
          do i=I_0,I_1
          tr3Dsource(i,j,l,:,n) = (MA(l,i,j)*axyp(i,j))*factor/DTsrc
          enddo
         enddo
        enddo

      call apply_tracer_3Dsource(nChemistry,n,.FALSE.)

        do j=J_0,J_1
          do i=I_0,I_1
           if (nint(lat2d_dg(i,j)).ge.30 .and.
     &         nint(lat2d_dg(i,j)).le.50) then
              trm(i,j,1,n) = 0.d0
           endif
          enddo
        enddo

      case ('aoa')

       tr3Dsource(:,:,:,:,n) = 0.

       factor = DTsrc/(60.d0*60.d0*24.d0)

       do l=1,lm
        do j=J_0,J_1
          do i=I_0,I_1
          tr3Dsource(i,j,l,:,n) = (MA(l,i,j)*axyp(i,j))*factor/DTsrc
          enddo
         enddo
        enddo

      call apply_tracer_3Dsource(nChemistry,n,.FALSE.)

      trm(:,:,1,n) = 0.d0

C****ST8025: An idealized loss tracer with a stratospheric source (fixed concentration above 80 mb) and idealized exponential decay.  Can be used to evaluate stratosphere-troposphere-exchange.

      case ('st8025')

      tr3Dsource(:,:,:,:,n) = 0.

      call exp_loss_trop(nChemistry,n)

      tr3Dsource(:,:,:,:,n) = 0.

       do l=1,lm
        if (nint(pmidl00(l)).le.82.) then
           do j=J_0,J_1
             do i=I_0,I_1
               tr3Dsource(i,j,l,:,n) = (MA(l,i,j)*axyp(i,j)*200.0d-9
     *        -trm(i,j,l,n))*bydt
             enddo
           enddo
        else
            tr3Dsource(:,:,l,:,n) = 0.
        endif
      enddo

      call apply_tracer_3Dsource(nChemistry,n,.FALSE.)


C****TAPE_REC: An idealized oscillating tracer in the tropical lower stratosphere whose period mimics the seasonal cycle of water vapor.  Can be used to evaluate vertical ascent in the lower stratosphere.

      case ('tape_rec')

      tr3Dsource(:,:,:,:,n) = 0.

       do l=1,lm
           if (nint(pmidl00(l)).ge.99. .and.
     &         nint(pmidl00(l)).lt.110.) then
           do j=J_0,J_1
             do i=I_0,I_1
            if (nint(lat2d_dg(i,j)).gt.-10 .and.
     &         nint(lat2d_dg(i,j)).lt.10) then
               xday=dayOfYear
               factor=(1.d0 + sin(2*(xday*pi/365.d0)-pi))*10.0d-9
               tr3Dsource(i,j,l,:,n) = (MA(l,i,j)*axyp(i,j)*factor
     *        -trm(i,j,l,n))*bydt
            else
              tr3Dsource(i,j,l,:,n) = 0.
            endif
             enddo
           enddo
           else
              tr3Dsource(:,:,l,:,n) = 0.
           endif
       enddo

      call apply_tracer_3Dsource(nChemistry,n,.FALSE.)

#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_TOMAS)
      case ('Alkenes', 'CO', 'NOx', 'Paraffin','CH4','codirect',
#ifdef TRACERS_dCO
     *      'd13Calke','d13CPAR',
     *      'dC17O', 'dC18O', 'd13CO',
#endif  /* TRACERS_dCO */
     &      'NH3', 'SO2', 'SO4', 'BCII', 'BCB', 'OCII', 'OCB',
     &      'vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &      'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6',
     &      'M_ACC_SU', 'M_AKK_SU',
     &      'M_BC1_BC', 'M_OCC_OC', 'M_BOC_BC', 'M_BOC_OC'
     &      ,'ASO4__01','AECOB_01','AOCOB_01')

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
C**** 3D volcanic source
        select case (trname(n))
        case ('SO2', 'SO4', 'M_ACC_SU', 'M_AKK_SU')
          tr3Dsource(:,J_0:J_1,:,nVolcanic,n)=
     &      sum(so2_src_3d(:,J_0:J_1,:,:),4)*src_fact
          call apply_tracer_3Dsource(nVolcanic,n)
! Water not implemented yet, and it might never be in this branch.
!          select case(trname(n))
!          case ('SO2') ! apply changes to q. do not update qmom,
!                       ! since the concentration always increases.
!            q(i,j,:)=XXX
!#ifdef TRACERS_WATER
!C**** Add water to relevant tracers as well
!            do nn=1,ntm
!              select case (tr_wd_type(n))
!              case (nWater)       ! water: initialise tracers
!                trm(i,j,:,n)=trm(i,j,:,n)*
!     &            XXX
!              end select
!            end do
!#endif
!          end select
        end select
#endif
C**** 3D biomass source
        tr3Dsource(:,J_0:J_1,:,nBiomass,n) = 0.
        if(do_fire(src_index) .or. nBBsources(src_index) > 0) then
          bb_i=ntsurfsrc(src_index)+1 ! index of first BB source
          if(do_fire(src_index))then
            bb_e=bb_i ! index last BB source
          else
            bb_e=ntsurfsrc(src_index)+nBBsources(src_index) ! index last BB source
          end if
          do j=J_0,J_1; do i=I_0,I_1
            blay=int(dclev(i,j)+0.5d0)
            blsrc = axyp(i,j)*src_fact*
     &       sum(sfc_src(i,j,src_index,bb_i:bb_e))/sum(MA(1:blay,i,j))
            do l=1,blay
             tr3Dsource(i,j,l,nBiomass,n)=tune_BBsources*blsrc*MA(l,i,j)
            end do
          end do; end do
        end if
#ifndef TRACERS_TOMAS
        call apply_tracer_3Dsource(nBiomass,n)
#endif

#ifdef TRACERS_TOMAS
        if(n<n_ASO4(1)) call apply_tracer_3Dsource(nBiomass,n)

!Initialize
       TOMAS_bio(:,J_0:J_1,:,:)=0.0
       TOMAS_air(:,J_0:J_1,:,:)=0.0


        select case (trname(n))
        case ('ASO4__01')

       do kk=1,nbins
         TOMAS_bio(:,J_0:J_1,:,kk)=
     &        tr3Dsource(:,J_0:J_1,:,nBiomass,n_ASO4(1))
     &        *scalesizeSO4_bio(kk)
       enddo

       do k=1,nbins
         tr3Dsource(:,J_0:J_1,:,nVolcanic,n_ASO4(1)+k-1)=
     &     sum(so2_src_3d(:,J_0:J_1,:,:),4)*scalesizeSO4_vol(k)*src_fact

         tr3Dsource(:,J_0:J_1,:,nBiomass,n_ASO4(1)+k-1)=
     *        TOMAS_bio(:,J_0:J_1,:,k)

         tr3Dsource(:,J_0:J_1,:,1,n_ANUM(1)+k-1)=
     &        (tr3Dsource(:,J_0:J_1,:,nVolcanic,n_ASO4(1)+k-1)
     &        +tr3Dsource(:,J_0:J_1,:,nBiomass,n_ASO4(1)+k-1))
     &        /(sqrt(xk(k)*xk(k+1)))

       enddo
       end select
#endif

#endif /* TRACERS_AEROSOLS_Koch || TRACERS_AMP || TRACERS_SPECIAL_Shindell || TRACERS_TOMAS*/

#if (defined TRACERS_COSMO)
C****
      case ('Be7')
c cosmogenic src
        do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
          tr3Dsource(i,j,l,nChemistry,n) = MA(l,i,j)*be7_src_3d(i,j,l)
        end do; end do; end do

        call apply_tracer_3Dsource(nChemistry,n)
C****
      case ('Be10')
c cosmogenic src
        do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
           tr3Dsource(i,j,l,nChemistry,n) = MA(l,i,j)*be10_src_3d(i,j,l)
        end do; end do; end do

        call apply_tracer_3Dsource(nChemistry,n)
C****
#endif
      case('Pb210')
        call apply_tracer_3Dsource(nChemistry,n) !radioactive decay of Rn222

      end select

      end do

#if (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
c Calculation of gas phase reaction rates for sulfur chemistry
      CALL GET_SULF_GAS_RATES
#endif

#ifdef TRACERS_SPECIAL_Shindell
C Apply non-chemistry 3D sources, so they can be "seen" by chemistry:
C (Note: using this method, tracer moments are changed just like they
C are done for chemistry.  It might be better to do it like surface
C sources are done? -- GSF 11/26/02)
c
      CALL TIMER (NOW,MTRACE)

#ifdef SHINDELL_STRAT_EXTRA
      tr3Dsource(I_0:I_1,J_0:J_1,:,nOverwrite,n_GLT) = 0.d0
      call overwrite_GLT
      call apply_tracer_3Dsource(nOverwrite,n_GLT)
#endif
#endif /* TRACERS_SPECIAL_Shindell */

#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS) ||\
    (defined TRACERS_GASEXCH_GCC)
      !  Aircraft Sources Here: All Tracers! (formerly just hardcoded set allowed)
      do n=1,ntm
        src_index=get_src_index(n)
        if(do_aircraft(src_index)) then
          xday=dayOfYear
          tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n)  = 0.d0
#ifdef CUBED_SPHERE
          call get_aircraft_tracer ! logical read from disk
     &     (n,trim(trname(src_index))//'_AIRC',year,xday,
     &      dummy3d,.false.,AIRCstreams(n))
#else
          call get_aircraft_tracer
     &     (n,trim(trname(src_index))//'_AIRC',year,xday,
     &      phi,.true.,AIRCstreams(n))
#endif
#ifdef TRACERS_TOMAS
          ! TOMAS has to apply this among tracers in its own section below.
          if(n /= n_AECOB(1))
     &    call apply_tracer_3Dsource(nAircraft,n)
#elif defined GCC_ZERO_EMISSION
#else
          call apply_tracer_3Dsource(nAircraft,n)
#endif
        end if
      end do
#endif /* Shindell or Koch or AMP or TOMAS or GASEXCH_GCC */

#ifdef TRACERS_SPECIAL_Shindell
      tr3Dsource(I_0:I_1,J_0:J_1,:,nOther,n_NOx) = 0.d0
      call get_lightning_NOx
      call apply_tracer_3Dsource(nOther,n_NOx)

C**** Make sure that these 3D sources for all chem tracers start at 0.:
      ! I think this zeroing is more important, now that the chemistry
      ! may not reach the top model layers:
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,ntm_chem_beg:ntm_chem_end)
     &  = 0.d0
      tr3Dsource(I_0:I_1,J_0:J_1,:,nOverwrite,ntm_chem_beg:ntm_chem_end)
     &  = 0.d0
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_stratOx)  = 0.d0
      tr3Dsource(I_0:I_1,J_0:J_1,:,nOverwrite,n_stratOx) = 0.d0
#endif
#if (defined TRACERS_HETCHEM) && (defined TRACERS_NITRATE)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_N_d1)  = 0.d0
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_N_d2)  = 0.d0
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_N_d3)  = 0.d0
#endif
#ifdef TRACERS_AEROSOLS_SOA
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_soa_i:n_soa_e)  = 0.d0
#endif  /* TRACERS_AEROSOLS_SOA */

      if (is_set_param('initial_ghg_setup')) then
        call get_param('initial_GHG_setup', initial_GHG_setup)
        if (initial_GHG_setup == 1 .and. itime == itimeI) then

          if (use_rad_n2o > 0) call applyRadChem(3, n_N2O, 1.d+0)
          if (use_rad_ch4 > 0) call applyRadChem(4, n_CH4, 1.d+0)
          if (use_rad_cfc > 0) call applyRadChem(5, n_CFC, fact_CFC)

        end if
      end if

C**** Call the model CHEMISTRY and OVERWRITEs:
      call masterchem ! does chemistry and over-writing.
                      ! tr3Dsource defined within, for both processes

C**** Apply chemistry and overwrite changes:
      do n=ntm_chem_beg, ntm_chem_end
        call apply_tracer_3Dsource(nChemistry,n)
        call apply_tracer_3Dsource(nOverwrite,n)
      end do
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      call apply_tracer_3Dsource(nChemistry,n_stratOx)
      call apply_tracer_3Dsource(nOverwrite,n_stratOx)
#endif
#if (defined TRACERS_HETCHEM) && (defined TRACERS_NITRATE)
       call apply_tracer_3Dsource(nChemistry,n_N_d1) ! NO3 chem prod on dust
       call apply_tracer_3Dsource(nChemistry,n_N_d2) ! NO3 chem prod on dust
       call apply_tracer_3Dsource(nChemistry,n_N_d3) ! NO3 chem prod on dust
#endif
      CALL TIMER (NOW,MCHEM)
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_NITRATE
#ifdef TRACERS_SPECIAL_Shindell
       tr3Dsource(I_0:I_1,J_0:J_1,:,3,n_HNO3) = 0.d0
#endif
       tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_NO3p) = 0.d0
       tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_NH4)  = 0.d0
       tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_NH3)  = 0.d0

#ifdef TRACERS_SPECIAL_Shindell
       call NITRATE_THERMO_DRV(topLevelOfChemistry)
#else
       call NITRATE_THERMO_DRV(LM)
#endif

#ifdef TRACERS_SPECIAL_Shindell
       call apply_tracer_3Dsource(3,n_HNO3) ! NO3 chem prod
#endif
       call apply_tracer_3Dsource(nChemistry,n_NO3p) ! NO3 chem prod
       call apply_tracer_3Dsource(nChemistry,n_NH4)  ! NO3 chem prod
       call apply_tracer_3Dsource(nChemistry,n_NH3)  ! NH3

#endif /* TRACERS_NITRATE */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) || \
    (defined TRACERS_TOMAS)
       call aerosol_gas_chem
#endif

#ifdef TRACERS_TOMAS
!H2SO4 chem prod is zero for TOMAS (H2SO4 will use directly in TOMAS_DRV)
!But it still calls to save the diagnostics.
       call apply_tracer_3Dsource(nChemistry,n_H2SO4) ! H2SO4 chem prod
       call apply_tracer_3Dsource(nChemistry,n_DMS)  ! DMS chem sink
       call apply_tracer_3Dsource(nChemistry,n_MSA)  ! MSA chem source
       call apply_tracer_3Dsource(nChemistry,n_SO2)  ! SO2 chem source
       call apply_tracer_3Dsource(nChemloss,n_SO2)  ! SO2 chem sink
       call apply_tracer_3Dsource(nChemistry,n_H2O2_s) ! H2O2 chem source
       call apply_tracer_3Dsource(2,n_H2O2_s) ! H2O2 chem sink
! EC/OC aging

       do k=1,nbins
          call apply_tracer_3Dsource(nChemistry,n_AECOB(k))
          call apply_tracer_3Dsource(nChemistry,n_AECIL(k))
          call apply_tracer_3Dsource(nChemistry,n_AOCOB(k))
          call apply_tracer_3Dsource(nChemistry,n_AOCIL(k))
       enddo

       do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
         trm_emis(i,j,l,:)=trm(i,j,l,:)
       end do; end do; end do

       TOMAS_bio(:,J_0:J_1,:,:)=0.0
       TOMAS_air(:,J_0:J_1,:,:)=0.0

       do kk=1,nbins
         TOMAS_bio(:,J_0:J_1,:,kk)=
     &        tr3Dsource(:,J_0:J_1,:,nBiomass,n_AECOB(1))
     &        *scalesizeCARBO100(kk)
c$$$
       enddo
       if(do_aircraft(n_AECOB(1)))then
         do kk=1,nbins
           TOMAS_air(:,J_0:J_1,:,kk)=
     &        tr3Dsource(:,J_0:J_1,:,nAircraft,n_AECOB(1))
     &        *scalesizeCARBO30(kk)

         enddo
       endif

       !TODO: once reproducibility is determined, pull these
       ! if's out of the k loop and do a second conditional
       ! k-loop instead:
       do k=1,nbins

         tr3Dsource(:,J_0:J_1,:,nBiomass,n_AECOB(1)+k-1)=
     *        TOMAS_bio(:,J_0:J_1,:,k)*0.8
         tr3Dsource(:,J_0:J_1,:,nBiomass,n_AECIL(1)+k-1)=
     *        TOMAS_bio(:,J_0:J_1,:,k)*0.2

         if(do_aircraft(n_AECOB(1))) then
           tr3Dsource(:,J_0:J_1,:,nAircraft,n_AECOB(1)+k-1)=
     *        TOMAS_air(:,J_0:J_1,:,k)*0.8
           tr3Dsource(:,J_0:J_1,:,nAircraft,n_AECIL(1)+k-1)=
     *        TOMAS_air(:,J_0:J_1,:,k)*0.2
         end if

         ! Here TOMAS_air() would be 0 when do_aircraft(n_AECOB(1)) is false,
         ! so leaving it unconditional:
         tr3Dsource(:,J_0:J_1,:,2,n_ANUM(1)+k-1)=
     &      (TOMAS_bio(:,J_0:J_1,:,k)+TOMAS_air(:,J_0:J_1,:,k))
     &      /(sqrt(xk(k)*xk(k+1)))

         call apply_tracer_3Dsource(nBiomass, n_AECOB(1)+k-1)
         if(do_aircraft(n_AECOB(1)))
     &    call apply_tracer_3Dsource(nAircraft,n_AECOB(1)+k-1)
         call apply_tracer_3Dsource(nBiomass, n_AECIL(1)+k-1)
         if(do_aircraft(n_AECOB(1)))
     &    call apply_tracer_3Dsource(nAircraft,n_AECIL(1)+k-1)
         call apply_tracer_3Dsource(2,       n_ANUM(1)+k-1)

         call apply_tracer_3Dsource(nVolcanic,n_ASO4(1)+k-1)
         call apply_tracer_3Dsource(nBiomass, n_ASO4(1)+k-1)
         call apply_tracer_3Dsource(1,       n_ANUM(1)+k-1)

       enddo

       TOMAS_bio(:,J_0:J_1,:,:)=0.0
       TOMAS_air(:,J_0:J_1,:,:)=0.0

       do kk=1,nbins
         TOMAS_bio(:,J_0:J_1,:,kk)=
     &        tr3Dsource(:,J_0:J_1,:,nBiomass,n_AOCOB(1))
     &        *scalesizeCARBO100(kk)
       enddo

       do k=1,nbins

         tr3Dsource(:,J_0:J_1,:,nBiomass,n_AOCOB(1)+k-1)=
     *        TOMAS_bio(:,J_0:J_1,:,k)*0.5
         tr3Dsource(:,J_0:J_1,:,nBiomass,n_AOCIL(1)+k-1)=
     *        TOMAS_bio(:,J_0:J_1,:,k)*0.5

         tr3Dsource(:,J_0:J_1,:,4,n_ANUM(1)+k-1)=
     &        (TOMAS_bio(:,J_0:J_1,:,k)
     &        )/(sqrt(xk(k)*xk(k+1)))


         call apply_tracer_3Dsource(nBiomass, n_AOCOB(1)+k-1)
         call apply_tracer_3Dsource(nBiomass, n_AOCIL(1)+k-1)
!     ntsurfsrc(n=3) is used for microphysics, so it is 4.
         call apply_tracer_3Dsource(4,       n_ANUM(1)+k-1)

       enddo

!for debugging!
!       do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
       call subgridcoag_drv(dtsrc)
!       end do; end do; end do

      DO n=1,ntm_TOMAS
!         if(am_i_root()) print*,'tr3dsource',trname(n_ASO4(1)+n-1)
        tr3Dsource(I_0:I_1,J_0:J_1,:,nOther,n_ASO4(1)+n-1) = 0.d0! Aerosol Mirophysics
      ENDDO
        tr3Dsource(I_0:I_1,J_0:J_1,:,nOther,n_H2SO4) = 0.d0! Aerosol Mirophysics
        tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_NH3) = 0.d0! Aerosol Mirophysics
        tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_NH4) = 0.d0! Aerosol Mirophysics
        tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_SOAgas) = 0.d0! Aerosol Mirophysics
c$$$#ifdef  TRACERS_SPECIAL_Shindell
c$$$        tr3Dsource(:,J_0:J_1,:,3,n_HNO3)  = 0.d0! Aerosol Mirophysics
c$$$#endif

        call TOMAS_DRV
        if(am_i_root()) print*,'exit TOMAS DRV'

      DO n=1,ntm_TOMAS-nbins ! exclude h2o
        call apply_tracer_3Dsource(nOther,n_ASO4(1)+n-1)! Aerosol Mirophysics
      ENDDO

      call apply_tracer_3Dsource(nChemistry,n_NH3)  !simple equilibrium model in TOMAS
      call apply_tracer_3Dsource(nChemistry,n_NH4) !simple equilibrium model in TOMAS
      call apply_tracer_3Dsource(nOther,n_H2SO4) ! H2SO4 chem prod
      call apply_tracer_3Dsource(nChemistry,n_SOAgas) ! H2SO4 chem prod
c$$$#ifdef  TRACERS_SPECIAL_Shindell
c$$$       call apply_tracer_3Dsource(3,n_HNO3) ! H2SO4 chem prod
c$$$#endif
C       stop

#endif /* TRACERS_TOMAS */

#ifdef TRACERS_AEROSOLS_Koch
       call apply_tracer_3Dsource(nChemistry,n_DMS)  ! DMS chem sink
       call apply_tracer_3Dsource(nChemistry,n_MSA)  ! MSA chem source
       call apply_tracer_3Dsource(nChemistry,n_SO2)  ! SO2 chem source
       call apply_tracer_3Dsource(nChemloss,n_SO2)  ! SO2 chem sink
#ifdef ACCMIP_LIKE_DIAGS
       do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
         taijls(i,j,l,ijlt_prodSO4gs)=taijls(i,j,l,ijlt_prodSO4gs)+
     &   tr3Dsource(i,j,l,nChemistry,n_SO4)*dtsrc*byaxyp(i,j)
       end do; end do; end do
#endif
       call apply_tracer_3Dsource(nChemistry,n_SO4)  ! SO4 chem source
       call apply_tracer_3Dsource(1,n_H2O2_s) ! H2O2 chem source
       call apply_tracer_3Dsource(2,n_H2O2_s) ! H2O2 chem sink
       call apply_tracer_3Dsource(nChemistry,n_BCII) ! BCII aging sink
       call apply_tracer_3Dsource(nChemistry,n_BCIA) ! BCIA aging source
#ifdef TRACERS_AEROSOLS_VBS
       do i=1,vbs_tr%nbins
         call apply_tracer_3Dsource(nChemistry,vbs_tr%igas(i)) ! aging source
         call apply_tracer_3Dsource(nChemloss,vbs_tr%igas(i))  ! aging loss
         call apply_tracer_3Dsource(nOther,vbs_tr%igas(i))     ! partitioning
         call apply_tracer_3Dsource(nChemistry,vbs_tr%iaer(i)) ! partitioning
       enddo
#else
       call apply_tracer_3Dsource(nChemistry,n_OCII) ! OCII aging sink
       call apply_tracer_3Dsource(nChemistry,n_OCIA) ! OCIA aging source
#endif

#ifdef TRACERS_HETCHEM
       call apply_tracer_3Dsource(nChemistry,n_SO4_d1) ! SO4 chem prod on dust
       call apply_tracer_3Dsource(nChemistry,n_SO4_d2) ! SO4 chem prod on dust
       call apply_tracer_3Dsource(nChemistry,n_SO4_d3) ! SO4 chem prod on dust
#endif
#endif  /* TRACERS_AEROSOLS_Koch */

#ifdef TRACERS_AMP
       call apply_tracer_3Dsource(nChemistry,n_H2SO4) ! H2SO4 chem prod
       call apply_tracer_3Dsource(nChemistry,n_DMS)  ! DMS chem sink
       call apply_tracer_3Dsource(nChemistry,n_SO2)  ! SO2 chem source
       call apply_tracer_3Dsource(nChemloss,n_SO2)  ! SO2 chem sink
       call apply_tracer_3Dsource(nChemistry,n_H2O2_s)! H2O2 chem source
       call apply_tracer_3Dsource(2,n_H2O2_s)! H2O2 chem sink
      DO n=ntmAMPi,ntmAMPe
        tr3Dsource(:,J_0:J_1,:,nChemistry,n)  = 0.d0! Aerosol Mirophysics !kt is the index correct?
      ENDDO
        tr3Dsource(:,J_0:J_1,:,nOther,n_H2SO4)  = 0.d0! Aerosol Mirophysics
        tr3Dsource(:,J_0:J_1,:,nChemistry,n_NH3)  = 0.d0! Aerosol Mirophysics !kt is the index correct?
#ifdef  TRACERS_SPECIAL_Shindell
        tr3Dsource(:,J_0:J_1,:,3,n_HNO3)  = 0.d0! Aerosol Mirophysics
#endif
        call MATRIX_DRV

      DO n=ntmAMPi,ntmAMPe
       call apply_tracer_3Dsource(nChemistry,n) ! Aerosol Mirophysics !kt is the index correct?
      ENDDO

       call apply_tracer_3Dsource(nOther,n_H2SO4) ! H2SO4 chem prod
       call apply_tracer_3Dsource(nChemistry,n_NH3)  ! NH3 !kt is the index correct?
#ifdef  TRACERS_SPECIAL_Shindell
       call apply_tracer_3Dsource(nOther,n_HNO3) ! HNO3 chem prod
#endif
#endif /* TRACERS_AMP */

#ifdef CACHED_SUBDD
      ! Accumulate the tracer-related subdaily diagnostics
      ! (seems like a reasonable place to put this, as tracer_3Dsource
      ! is called each time step, and here chemistry has been done, etc.,
      ! but we could move it):
      call accumCachedTracerSUBDDs
#endif

      return

#ifdef TRACERS_SPECIAL_Shindell
      contains

      subroutine applyRadChem(index, n, factor)
      integer, intent(in) :: index
      integer, intent(in) :: n
      real*8, intent(in) :: factor

      integer :: L
      do L = 1, LM
        tr3Dsource(I_0:I_1,J_0:J_1,L,nOverwrite,n) =
     &       (rad_to_chem(index,L,I_0:I_1,J_0:J_1)*2.69e20*byavog*
     &       axyp(I_0:I_1,J_0:J_1)*tr_mm(n) * factor -
     &       trm(I_0:I_1,J_0:J_1,L,n)) / dtsrc
      end do
      call apply_tracer_3Dsource(nOverwrite,n)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nOverwrite,n) = 0.d0

      end subroutine applyRadChem
#endif

      END SUBROUTINE tracer_3Dsource
#endif /* TRACERS_ON */

#ifdef TRACERS_WATER
C---SUBROUTINES FOR TRACER WET DEPOSITION-------------------------------

      SUBROUTINE GET_COND_FACTOR(
     &     NTX,WMXTR,TEMP,TEMP0,LHX,FCLOUD,
     &     FQ0,fq,TR_CONV,TRWML,TM,THLAW,TR_LEF,pl,ntix,CLDSAVT)
!@sum  GET_COND_FACTOR calculation of condensate fraction for tracers
!@+    within or below convective or large-scale clouds. Gas
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
! NOTE: THLAW is only computed for the tracers in hlaw_list!
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only: BYGASC, MAIR,teeny,LHE,tf,by3

      USE TRACER_COM, only :
     &     aero_count,water_count,hlaw_count,
! NB: these lists are often used for implicit loops
     &     aero_list,water_list,hlaw_list

      use OldTracer_mod, only: tr_RKD, tr_DHD, tr_wd_type
      use OldTracer_mod, only: nWater, ngas,nPART
      use OldTracer_mod, only: trname, t_qlimit, fq_aer, trpdens
      USE TRACER_COM, only:
     *     NTM,n_SO2,n_H2O2,n_H2O2_s
#ifdef TRACERS_SPECIAL_O18
      USE TRACER_COM, only: supsatfac
#endif
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only: NBS,NBINS
     &     ,n_ANUM,n_ASO4,n_ANACL,n_AECIL,n_AECOB
     &     ,n_AOCIL,n_AOCOB,n_ADUST,n_AH2O
#endif
#ifdef TRACERS_HETCHEM
      USE TRACER_COM, only: trm ,n_SO4_d1, n_SO4_d2, n_SO4_d3,n_SO4
     *     ,n_N_d1,n_N_d2,n_N_d3,n_NO3p, n_Clay,n_Silt1,n_Silt2
      USE MODEL_COM, only  : dtsrc
#endif
      use OldTracer_mod, only: set_fq_aer
      IMPLICIT NONE
C**** Local parameters and variables and arguments:
!@param BY298K unknown meaning for now (assumed= 1./298K)
!@var Ppas pressure at current altitude (in Pascal=kg/s2/m)
!@var TFAC exponential coeffiecient of tracer condensation temperature
!@+   dependence (mole/joule)
!@var FCLOUD fraction of cloud available for tracer condensation
!@var SSFAC dummy variable (assumed units= kg water?)
!@var FQ            fraction of tracer that goes into condensate
!@var FQ0 default fraction of water tracer that goes into condensate
!@var L index for altitude loop
!@var N index for tracer number loop
!@var WMXTR mixing ratio of water available for tracer condensation?
!@var SUPSAT super-saturation ratio for cloud droplets
!@var LHX latent heat flag for whether condensation is to ice or water
!@var RKD dummy variable (= tr_RKD*EXP[ ])
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, RKD,CLDINC,trlef
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvs,fracvl,kin_cond_ice
#endif
      REAL*8,  INTENT(IN) :: fq0, FCLOUD, WMXTR, TEMP, TEMP0,LHX
     &     , TR_LEF(NTM), pl,CLDSAVT
      REAL*8,  INTENT(IN), DIMENSION(NTM) :: trwml
      REAL*8,  INTENT(IN), DIMENSION(NTM) :: TM
      REAL*8,  INTENT(OUT):: fq(NTM),thlaw(NTM)
      INTEGER, INTENT(IN) :: NTX, ntix(NTM)
      LOGICAL TR_CONV
      REAL*8 :: FQ0FAC,SUPSAT,SSFAC(NTM),SSFAC0
      INTEGER :: N,IHLAW,IAERO,IWAT
#ifdef TRACERS_TOMAS
      integer :: k
      real*8,dimension(nbins):: fraction !where to read fraction?
#endif
c      thlaw(:) = 0. ! default

#if (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_SPECIAL_Shindell) ||\
    (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
c
c gases with a henry law constant
c
c     cldinc=max(0.,cldsavt-fcloud)
      if(lhx.eq.lhe .and. fcloud.ge.1d-16) then
        Ppas = PL*1.D2          ! pressure to pascals
        tfac = (1.D0/TEMP - BY298K)*BYGASC
        ssfac0 = WMXTR*MAIR*1.D-6*Ppas/(CLDSAVT+teeny)
        ssfac(hlaw_list) = ssfac0*tr_RKD(hlaw_list)
     &    *exp(-tr_DHD(hlaw_list)*tfac)
        if(tr_conv) then ! convective cloud
          fq0fac = 1.
          if (fq0.eq.0.) fq0fac=0.d0
          do ihlaw=1,hlaw_count
            n = hlaw_list(ihlaw)
            fq(n) = fq0fac*ssfac(n) / (1d0 + ssfac(n))
            thlaw(n) = 0.
          enddo
        else             ! stratiform cloud
          do ihlaw=1,hlaw_count
            n = hlaw_list(ihlaw)
            fq(n) = 0.
c limit gas dissolution to incremental cloud change after cloud forms
c   only apply to non-aqueous sulfur species since this is already
c   done in GET_SULFATE
c but H2O2 should be limited if not coupled with sulfate, have not done this
c           if (n.ne.n_h2O2.and.n.ne.n_so2.and.n.ne.n_h2O2_s) then
c           if (FCLOUD.ne.0.) tr_lef(n)=cldinc
c           endif
            trlef=min(tr_lef(n),cldsavt)
            thlaw(n) = min(tm(n),max(0d0,
     &           (ssfac(n)*trlef*tm(n)-TRWML(n))
     &           /(1.D0+ssfac(n)) ))
          enddo
        endif
      else
        fq(hlaw_list) = 0.
        thlaw(hlaw_list) = 0.
      endif
#endif /* dissolved gases with a henry law constant */

c
c loop over water species
c
      do iwat=1,water_count
        n = water_list(iwat)
        fq(n) = fq0
#ifdef TRACERS_SPECIAL_O18
          if (fq0.gt.0. .and. fq0.lt.1.) then
C**** If process occurs at constant temperature, calculate condensate
C**** in equilibrium with source vapour. Otherwise, use mid-point
C**** temperature and estimate instantaneous fractionation. This gives
C**** a very good estimate to complete integral
C****
            if (abs(temp-temp0).gt.1d-14) then  ! use instantaneous frac
              tdegc=0.5*(temp0 + temp) -tf
C**** Calculate alpha (fractionation coefficient)
                if (LHX.eq.LHE) then ! cond to water
                  alph=1./fracvl(tdegc,ntix(n))
                else            ! cond to ice
                  alph=1./fracvs(tdegc,ntix(n))
C**** kinetic fractionation can occur as a function of supersaturation
C**** this is a parameterisation from Georg Hoffmann
                  supsat=1d0-supsatfac*tdegc
                  if (supsat .gt. 1.) alph=kin_cond_ice(alph,supsat
     *                 ,ntix(n))
                end if
                fq(n) = 1.- (1.-fq0)**alph
            else
C**** assume condensate in equilibrium with vapour at temp
              tdegc=temp -tf
              if (LHX.eq.LHE) then ! cond to water
                alph=1./fracvl(tdegc,ntix(n))
              else              ! cond to ice
                alph=1./fracvs(tdegc,ntix(n))
C**** kinetic fractionation can occur as a function of supersaturation
C**** this is a parameterisation from Georg Hoffmann
                supsat=1d0-supsatfac*tdegc
                if (supsat .gt. 1.) alph=kin_cond_ice(alph,supsat
     *               ,ntix(n))
              end if
              fq(n) = alph * fq0/(1.+(alph-1.)*fq0)
            end if
          else
            fq(n) = fq0
          end if
#endif
      enddo ! end loop over water species

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_COSMO) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_RADON) ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)

c
c aerosols
c

c     if (FCLOUD.lt.1.D-16 .or. fq0.eq.0.) then
      if (CLDSAVT.lt.1.D-16 .or. fq0.eq.0.) then
        fq(aero_list) = 0.      ! defaults to zero
      else

#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_DUST) &&\
    (defined TRACERS_HETCHEM)

      n = n_Clay
      if ( ( TM(ntix(n_SO4_d1)) /trpdens(n_SO4)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        call set_fq_aer(NTIX(N), 1.d0)
      else
        call set_fq_aer(NTIX(N), 0.d0)
      endif
#ifdef TRACERS_NITRATE
      if ( ( TM(ntix(n_N_d1)) /trpdens(n_NO3p)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        call set_fq_aer(NTIX(N), 1.d0)
      endif
#endif

      n = n_Silt1
      if ( ( TM(ntix(n_SO4_d2)) /trpdens(n_SO4)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        call set_fq_aer(NTIX(N), 1.d0)
      else
        call set_fq_aer(NTIX(N), 0.d0)
      endif
#ifdef TRACERS_NITRATE
      if ( ( TM(ntix(n_N_d2)) /trpdens(n_NO3p)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        call set_fq_aer(NTIX(N), 1.d0)
      endif
#endif

      n = n_Silt2
      if ( ( TM(ntix(n_SO4_d3)) /trpdens(n_SO4)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        call set_fq_aer(NTIX(N), 1.d0)
      else
        call set_fq_aer(NTIX(N), 0.d0)
      endif
#ifdef TRACERS_NITRATE
      if ( ( TM(ntix(n_N_d3)) /trpdens(n_NO3p)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        call set_fq_aer(NTIX(N), 1.d0)
      endif
#endif

#endif
#ifdef TRACERS_TOMAS

      if(tr_conv)then

         CALL getfraction (.true.,TM,FRACTION) !1% supersaturation assumption

      else                      ! large-scale clouds
         CALL getfraction (.false.,TM,FRACTION) !0.2% supersaturation assumption

      endif

      do k=1,nbins

        call set_fq_aer(ntix(n_ANUM(1)+k-1),fraction(k))
        call set_fq_aer(ntix(n_ASO4(1)+k-1), fraction(k))
        call set_fq_aer(ntix(n_ANACL(1)+k-1),  fraction(k))
        call set_fq_aer(ntix(n_AECIL(1)+k-1),fraction(k))
        call set_fq_aer(ntix(n_AECOB(1)+k-1),fraction(k))
        call set_fq_aer(ntix(n_AOCIL(1)+k-1),fraction(k))
        call set_fq_aer(ntix(n_AOCOB(1)+k-1),fraction(k))
        call set_fq_aer(ntix(n_ADUST(1)+k-1),fraction(k))
        call set_fq_aer(ntix(n_AH2O(1)+k-1), fraction(k))

         if (fraction(k).gt.1.or.fraction(k).lt.0) then
            print*,'fraction>1 or fraction<0'
            call stop_model('wrong fraction',255)
         endif

      enddo
#endif

      cldinc=cldsavt-fcloud
      if(tr_conv) then          ! convective cloud
c complete dissolution in convective clouds
c with double dissolution if partially soluble
        if(lhx.eq.lhe) then
          fq(aero_list) = fq_aer(aero_list)
        else
          fq(aero_list) = fq_aer(aero_list)*0.12d0
        endif
c this should not work because cldinc should be fcld for
c    when cloud first forms
      elseif(fq0.gt.0 .and. cldinc.gt.0.) then ! growing stratiform cloud
        if(lhx.eq.lhe) then
          fq(aero_list) = fq_aer(aero_list)*cldinc
        else
          fq(aero_list) = fq_aer(aero_list)*cldinc*0.12d0
        endif
      else
        fq(aero_list) = 0.
      endif
      where(fq(aero_list).ge.1.d0) fq(aero_list)=0.9999
c
c use this code in place of the above if the commented-out formulas
c for (dust?) fq are reinstated
c

c      do iaero=1,aero_count ! loop over aerosols
c        n = aero_list(iaero)
cc complete dissolution in convective clouds
cc with double dissolution if partially soluble
c          if (TR_CONV) then ! convective cloud
c            if (LHX.EQ.LHE) then !liquid cloud
ccdust #if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
ccdust           IF (fq_aer(ntix(n)) > 0.)
ccdust #endif
c              fq(n)=fq_aer(ntix(n))
ccdust?              fq(n)=(1.d0+fq_aer(ntix(n)))/2.d0
ccdust?              fq(n)=(1.d0+3.d0*fq_aer(ntix(n)))/4.d0
c            else
ccdust #if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
ccdust           IF (fq_aer(ntix(n)) > 0.)
ccdust #endif
c              fq(n)=fq_aer(ntix(n))*0.12d0
ccdust?              fq(n)=(1.d0+fq_aer(ntix(n)))/2.d0*0.05d0
ccdust?              fq(n)=(1.d0+3.d0*fq_aer(ntix(n)))/4.d0*0.05d0
c
c            endif
c          elseif (fq0.gt.0.and.CLDINC.gt.0.) then ! stratiform cloud.
cc only dissolve if the cloud has grown
c            if(LHX.EQ.LHE) then !liquid cloud
c              fq(n) = fq_aer(ntix(n))*CLDINC
c            else                ! ice cloud - small dissolution
c              fq(n) = fq_aer(ntix(n))*CLDINC*0.12d0
c            endif
c          endif
c          if (fq(n).ge.1.d0) fq(n)=0.9999
c      enddo ! end loop over aerosols

      endif ! fcloud>0 and fq0.ne.0

#endif /* aerosols */

      RETURN
      END SUBROUTINE GET_COND_FACTOR


      SUBROUTINE GET_WASH_FACTOR(NTX,b_beta_DT,PREC,fq
     * ,TEMP,LHX,WMXTR,FCLOUD,TM,TRPR,THLAW,pl,ntix,BELOW_CLOUD
#ifdef TRACERS_TOMAS
     * ,I,J,L
#endif
     *)
!@sum  GET_WASH_FACTOR calculation of the fraction of tracer
!@+    scavanged by precipitation below convective clouds ("washout").
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
! NOTE: THLAW is only computed for the tracers in hlaw_list!
! NOTE: FQ is only computed for the tracers in aero_list!
c
C**** GLOBAL parameters and variables:
      use OldTracer_mod, only: nWATER, ngas, nPART, tr_wd_type
      use OldTracer_mod, only: tr_RKD, tr_DHD, rc_washt, trname
      USE TRACER_COM, only:
     * NTM
#ifdef TRACERS_AEROSOLS_SEASALT
c     * n_seasalt1,n_seasalt2
c     USE PBLCOM, only: wsavg
#endif

      USE TRACER_COM, only :
     &     aero_count,water_count,hlaw_count,
! NB: these lists are often used for implicit loops
     &     aero_list,water_list,hlaw_list
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only :
     &     NBS,NBINS,n_ANUM,n_ASO4,n_ANACL,xk
     &    ,n_AOCOB,n_AECIL,n_AECOB,n_AOCIL,n_ADUST,n_AH2O
      use OldTracer_mod, only: set_rc_washt
#endif
      USE CONSTANT, only: BYGASC,LHE,MAIR,teeny,pi

      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var FQ fraction of tracer scavenged by below-cloud precipitation
!@param rc_wash aerosol washout rate constant (mm-1)
!@var PREC precipitation amount from layer above for washout (mm)
!@var b_beta_DT precipitating grid box fraction from lowest
!@+   percipitating layer.
!@+   The name was chosen to correspond to Koch et al. p. 23,802.
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: NTX,ntix(NTM)
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: THLAW
      REAL*8, INTENT(INOUT), DIMENSION(NTM) :: FQ
      REAL*8, INTENT(IN) :: PREC,b_beta_DT,TEMP,LHX,WMXTR,FCLOUD,
     *  TM(NTM),pl, TRPR(NTM)
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac0, ssfac(NTM), bb_tmp
      INTEGER :: N,IHLAW,IAERO
      LOGICAL BELOW_CLOUD
C
#ifdef TRACERS_TOMAS
      integer k,i,j,l
!@var scavr/stratsav : below-cloud scavenging coefficient (per mm rain)
      real*8 scavr
      real stratscav
!@var dpaero : aerosol diameter [m]
      real*8 dpaero,mtot
      real*8,dimension(nbins) ::  getdp,density
#endif
c      thlaw(:)=0.

c
c gases with a henry law constant
c
c      fq(hlaw_list) = 0.D0
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      if(      LHX.EQ.LHE ! if not frozen
     &   .AND. FCLOUD.GE.1D-16 .AND. WMXTR.GT.0. AND . BELOW_CLOUD) THEN
        bb_tmp = max(b_beta_DT,0.d0) ! necessary check?
        Ppas = PL*1.D2          ! pressure to pascals
        tfac = (1.D0/TEMP - BY298K)*BYGASC
        ssfac0 = WMXTR*MAIR*1.D-6*Ppas/(FCLOUD+teeny)
        ssfac(hlaw_list) = ssfac0*tr_RKD(hlaw_list)
     &    *exp(-tr_DHD(hlaw_list)*tfac)
        do ihlaw=1,hlaw_count
          n = hlaw_list(ihlaw)
          thlaw(n) = min( tm(n),max( 0d0,(FCLOUD*
     &               ssfac(n)*tm(n)-TRPR(n))/(1.D0+ssfac(n)) ))
        enddo
      else
        thlaw(hlaw_list) = 0.
      endif
#endif
#ifdef TRACERS_TOMAS

      if(FCLOUD.GE.1.D-16 .and. prec.gt.0.) then

         call dep_getdp(i,j,l,getdp,density) !1 for tempk, dummy=vs
         do k=1,nbins
            dpaero=getdp(k)
            scavr=stratscav(dpaero)
            call set_rc_washt(ntix(n_ASO4(1)+k-1), scavr)
            call set_rc_washt(ntix(n_ANACL(1)+k-1),  scavr)
            call set_rc_washt(ntix(n_AECOB(1)+k-1),scavr)
            call set_rc_washt(ntix(n_AECIL(1)+k-1),scavr)
            call set_rc_washt(ntix(n_AOCOB(1)+k-1),scavr)
            call set_rc_washt(ntix(n_AOCIL(1)+k-1),scavr)
            call set_rc_washt(ntix(n_ADUST(1)+k-1),scavr)
            call set_rc_washt(ntix(n_AH2O(1)+k-1), scavr)
            call set_rc_washt(ntix(n_ANUM(1)+k-1),scavr)
         enddo

      endif

#endif
c
c water species
c
c      fq(water_list) = 0d0

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_COSMO) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AEROSOLS_SEASALT) ||\
    (defined SHINDELL_STRAT_EXTRA) || (defined TRACERS_AMP) ||\
    (defined TRACERS_RADON) || (defined TRACERS_TOMAS)
c
c aerosols
c
      if(FCLOUD.GE.1.D-16 .and. prec.gt.0.) then
        bb_tmp = max(b_beta_DT,0.d0) ! necessary check?
        do iaero=1,aero_count
          n = aero_list(iaero)
          fq(n) = bb_tmp*(1d0-exp(-prec*rc_washt(n)))
        enddo
      else
        fq(aero_list) = 0.
      endif
#endif

      RETURN
      END SUBROUTINE GET_WASH_FACTOR

      SUBROUTINE GET_EVAP_FACTOR(
     &     NTX,TEMP,LHX,HEFF,FQ0,fq,ntix)
!@sum  GET_EVAP_FACTOR calculation of the evaporation fraction
!@+    for tracers.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only : tf,lhe
      use OldTracer_mod, only: tr_wd_type,nwater,trname
      USE TRACER_COM, only: NTM, tr_evap_fact, water_count,water_list
c      USE CLOUDS, only: NTIX
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var FQ            fraction of tracer evaporated
!@var FQ0 [default] fraction of tracer evaporated
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: NTX,ntix(NTM)
      REAL*8,  INTENT(OUT):: FQ(NTM)
      REAL*8,  INTENT(IN) :: FQ0,TEMP,LHX
!@var HEFF effective relative humidity for evap occuring below cloud
      REAL*8, INTENT(IN) :: HEFF
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvl,fracvs,kin_evap_prec
      integer :: iwat
#endif
      integer :: n
c

      if(fq0.ge.1d0) then
        fq(1:ntx) = 1d0 ! total evaporation
      else
        do n = 1, ntx
          fq(n) = fq0*tr_evap_fact(tr_wd_type(ntix(n)))
        end do
      endif

#ifdef TRACERS_SPECIAL_O18
c overwrite fq for water isotopes
      tdegc=temp-tf
      do iwat=1,water_count
        n = water_list(iwat)
        if (lhx.eq.lhe) then
          alph=fracvl(tdegc,ntix(n))
C**** kinetic effects with evap into unsaturated air
          if (heff.lt.1.)
     &         alph=kin_evap_prec(alph,heff,ntix(n))
        else
C**** no fractionation for ice evap
          alph=1.
        end if
        if (fq0.ne.1.) then
          fq(n) = 1. - (1.-fq0)**alph
        else
          fq(n) = fq0
        end if
      enddo
#endif

      RETURN
      END SUBROUTINE GET_EVAP_FACTOR

#endif


#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      SUBROUTINE GET_SULF_GAS_RATES
!@sum  GET_SULF_GAS_RATES calculation of rate coefficients for
!@+    gas phase sulfur oxidation chemistry
!@vers 2013/03/26
!@auth Bell
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only: t
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds,write_parallel
      USE ATM_COM, only: pmid,MA,pk
      USE GEOM, only: axyp,imaxj
      USE TRACER_COM, only: rsulf1,rsulf2,rsulf3,rsulf4
      INTEGER J_0, J_1, I_0, I_1
      real*8 ppres,te,tt,mm,dmm,rk4,ek4,f
C Greg: certain things now done outside the loops for speed:
      real*8, parameter ::  a= 73.41463d20, ! 6.02d20/.082d0
     *     aa=1.d-20,
     *     b= 0.357d-22,         ! 1.7d-22*0.21d0*1.d-20/aa
     *     c= 1.155d-11,         ! 5.5d-20*0.21d0*1.d-11/aa
     *     d= 4.0d-11            ! 4.0d-20*1.d-11/aa

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C Initialise
      rsulf1(:,j_0:j_1,:)=0.d0
      rsulf2(:,j_0:j_1,:)=0.d0
      rsulf3(:,j_0:j_1,:)=0.d0
      rsulf4(:,j_0:j_1,:)=0.d0

C Reactions
C***1.DMS + OH -> 0.75SO2 + 0.25MSA
C***2.DMS + OH -> SO2
C***3.DMS + NO3 -> HNO3 + SO2
C***4.SO2 + OH -> SO4 + HO2


      do l=1,LM
      do j=J_0,J_1
      do i=I_0,imaxj(j)

C Calculate effective temperature

        ppres=pmid(l,i,j)*9.869d-4 !in atm
        te=pk(l,i,j)*t(i,j,l)
        mm=MA(l,i,j)*axyp(i,j)
        tt = 1.d0/te

c DMM is number density of air in molecules cm-3

        dmm=ppres*tt*a
        rsulf1(i,j,l) =
     & b*dmm*exp(7810.d0*tt)*aa/(1.d0+c*exp(7460.d0*tt)*dmm*aa)

        rsulf2(i,j,l) = 9.6d-12*exp(-234.d0*tt)

        rsulf3(i,j,l) = 1.9d-13*exp(520.d0*tt)

        rk4 = aa*((tt*300.d0)**(3.3d0))*dmm*d
        f=log10(0.5d12*rk4)
        ek4 = 1.d0/(1.d0 + (f*f))

        rsulf4(i,j,l) = (rk4/(1.d0 + 0.5d12*rk4  ))*(0.45d0**ek4)

      end do
      end do
      end do

      END SUBROUTINE GET_SULF_GAS_RATES
#endif

#ifdef TRACERS_TOMAS

!    **************************************************
!@sum  initbounds
!    **************************************************
!@+    This subroutine initializes the array, xk, which describes the
!@+    boundaries between the aerosol size bins.  xk is in terms of dry
!@+    single-particle mass (kg).  The aerosol microphysics algorithm
!@+    used here assumes mass doubling such that each xk is twice the
!@+    previous.

!@auth Peter Adams, November 1999 (Modified by Yunha Lee)

      SUBROUTINE initbounds()



C-----INCLUDE FILES-----------------------------------------------------

      USE TRACER_COM,only : xk,nbins

C-----VARIABLE DECLARATIONS---------------------------------------------
      IMPLICIT NONE

      integer k
!@var Mo : lower mass bound for first size bin (kg)
      real*8 Mo

C-----ADJUSTABLE PARAMETERS---------------------------------------------

#if (defined TOMAS_12_3NM) || (defined TOMAS_15_3NM)
      parameter(Mo=1.5625d-23) ! 3nm
#else
      parameter(Mo=1.0d-21)    ! 10nm
#endif


C-----CODE--------------------------------------------------------------

      do k=1,nbins+1
!YUNHA LEE - working on adding more version of TOMAS (Aug, 2012)
#if (defined TOMAS_12_10NM) || (defined TOMAS_12_3NM)
         if(k.lt.nbins)then
            xk(k)=Mo*4.d0**(k-1)
         else
            xk(k)=xk(k-1)*32.d0
         endif
#elif (defined TOMAS_15_10NM) || (defined TOMAS_15_3NM)
           xk(k)=Mo*4.d0**(k-1)
#elif (defined TOMAS_30_10NM) || (defined TOMAS_30_3NM)
           xk(k)=Mo*2.d0**(k-1)
#endif
      enddo

      RETURN
      END


!    **************************************************
!@sum   momentfix
!    **************************************************
!@+    This routine changes the first and second order moments of a
!@+    given tracer's distribution such that they match the shape of
!@+    another specified tracer.  Since the zeroth order moment is
!@+    unchanged, the routine conserves mass.

!@auth Peter Adams

C-----INPUTS------------------------------------------------------------

!@var     pn - the number of the tracer that will serve as the pattern
!@var     fn - the number of the tracer whose moments will be fixed

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE momentfix(pn,fn)

C-----INCLUDE FILES-----------------------------------------------------

      USE QUSDEF, only : nmom
      USE RESOLUTION, ONLY : IM,JM,LM
      USE TRACER_COM, only : trm, trmom
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE GEOM, only : imaxj

      integer pn, fn

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,l,n
      INTEGER I_0, I_1, J_0, J_1
      real*8 ratio

C-----CODE--------------------------------------------------------------
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     $     I_STRT=I_0, I_STOP=I_1)

      do l=1,lm; do j=J_0,J_1; do i=I_0,imaxj(j)
         if (TRM(i,j,l,pn) .ge. 1.e-10) then
            ratio=TRM(i,j,l,fn)/TRM(i,j,l,pn)
            do n=1,nmom
               trmom(n,i,j,l,fn)=ratio*trmom(n,i,j,l,pn)
            enddo
         else

            do n=1,nmom
               trmom(n,i,j,l,fn)=0.0
               trmom(n,i,j,l,pn)=0.0
            enddo
         endif
      end do; end do; end do

      RETURN
      END



!    **************************************************
!@sum  stratscav
!@    **************************************************
!@+    This function is basically a lookup table to get the below-cloud
!@+    scavenging rate (per mm of rainfall) as a function of particle
!@+    diameter.  The data are taken from Dana, M. T., and
!@+    J. M. Hales, Statistical Aspects of the Washout of Polydisperse
!@+    Aerosols, Atmos. Environ., 10, 45-50, 1976.  I am using the
!@+    monodisperse aerosol curve from Figure 2 which assumes a
!@+    lognormal distribution of rain drops with Rg=0.02 cm and a
!@+    sigma of 1.86, values typical of a frontal rain spectrum
!@+    (stratiform clouds).

!@auth  Peter Adams, January 2001

      real FUNCTION stratscav(dp)

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS------------------------------------------
!@var dp : particle diameter [m]
      real*8 dp

C-----VARIABLE DECLARATIONS------------------------------------------
!@param numpts : number of points in lookup table
!@var dpdat : particle diameter in lookup table [m]
!@var scdat : scavenging rate in lookup table [mm-1]
!@var n1/n2 : indices of nearest data points
      integer numpts
      real dpdat
      real scdat
      integer n1, n2

C-----VARIABLE COMMENTS----------------------------------------------

C-----ADJUSTABLE PARAMETERS------------------------------------------

      parameter(numpts=37)
      dimension dpdat(numpts), scdat(numpts)

      data dpdat/ 2.0E-09, 4.0E-09, 6.0E-09, 8.0E-09, 1.0E-08,
     &            1.2E-08, 1.4E-08, 1.6E-08, 1.8E-08, 2.0E-08,
     &            4.0E-08, 6.0E-08, 8.0E-08, 1.0E-07, 1.2E-07,
     &            1.4E-07, 1.6E-07, 1.8E-07, 2.0E-07, 4.0E-07,
     &            6.0E-07, 8.0E-07, 1.0E-06, 1.2E-06, 1.4E-06,
     &            1.6E-06, 1.8E-06, 2.0E-06, 4.0E-06, 6.0E-06,
     &            8.0E-06, 1.0E-05, 1.2E-05, 1.4E-05, 1.6E-05,
     &            1.8E-05, 2.0E-05/

      data scdat/ 6.99E-02, 2.61E-02, 1.46E-02, 9.67E-03, 7.07E-03,
     &            5.52E-03, 4.53E-03, 3.87E-03, 3.42E-03, 3.10E-03,
     &            1.46E-03, 1.08E-03, 9.75E-04, 9.77E-04, 1.03E-03,
     &            1.11E-03, 1.21E-03, 1.33E-03, 1.45E-03, 3.09E-03,
     &            4.86E-03, 7.24E-03, 1.02E-02, 1.36E-02, 1.76E-02,
     &            2.21E-02, 2.70E-02, 3.24E-02, 4.86E-01, 8.36E-01,
     &            1.14E+00, 1.39E+00, 1.59E+00, 1.75E+00, 1.85E+00,
     &            1.91E+00, 1.91E+00/

C-----CODE-----------------------------------------------------------

C If particle diameter is in bounds, interpolate to find value
      if ((dp .gt. dpdat(1)) .and. (dp .lt. dpdat(numpts))) then
         !loop over lookup table points to find nearest values
         n1=1
         do while (dp .gt. dpdat(n1+1))
            n1=n1+1
         enddo
         n2=n1+1
         stratscav=scdat(n1)+(scdat(n2)-scdat(n1))
     &             *(dp-dpdat(n1))/(dpdat(n2)-dpdat(n1))
      endif

C If particle diameter is out of bounds, return reasonable value
!YUNHA - (TOMAS bug) I changed the condition that gt ==> ge.  lt ==> le
!YUNHA - Because stratscav has no value when dp=dpdat(numpts) and dpdat(1).

      if (dp .ge. dpdat(numpts)) stratscav=2.0
      if (dp .le. dpdat(1))      stratscav=7.0e-2

      RETURN
      END FUNCTION stratscav


#endif

!--------------------------------------------------------------------
! sph_mod -- generate spherical harmonics
!--------------------------------------------------------------------
      module sph_mod
      implicit none

      private

      type, public :: sph
        private
        real*8 :: base
        real*8, dimension(:), allocatable :: coeff
        integer :: l, m
      contains
        procedure :: build
        procedure :: disp
        procedure :: value
      end type sph

      contains

      function binomial(n, k)
      implicit none
      integer, intent(in) :: n, k
      integer :: binomial
      integer :: i, kk

      kk=merge(k, n-k, (k<n-k).and.(k<=n))
      binomial=1
      do i=1, kk
        binomial=binomial*(n+1-i)/i
      end do
      return
      end function binomial

      subroutine deriv(a)
      implicit none
      real*8, dimension(:), intent(inout) :: a
      integer :: i

      do i=1, size(a)-1
        a(i)=a(i+1)*i
      end do
      a(size(a))=0
      return
      end subroutine deriv

      subroutine mult(a, b, res)
      implicit none
      real*8, dimension(:), intent(in) :: a, b
      real*8, dimension(:), allocatable, intent(out) :: res
      integer :: i, j, sz, r0, r1
      real*8 :: val

      sz=size(a)+size(b)-1
      allocate(res(sz))
      do i=0, sz-1
        r0=max(0, i+1-size(b))
        r1=min(i, size(a)-1)
        val=0
        do j=r0, r1
          val=val+a(j+1)*b(i-j+1)
        end do
        res(i+1)=val
      end do
      return
      end subroutine mult

      subroutine legendre_poly(n, coeff)
      implicit none
      integer, intent(in) :: n
      real*8, dimension(:), allocatable, intent(out) :: coeff
      integer :: i, j, k
      real*8 :: val

      allocate(coeff(n+1))
      coeff=0
      do i=mod(n, 2), n, 2
        val=1
        k=merge(i, n-i, i<n-i)
        do j=1, k
          val=val*(n+1-j)/j
        end do
        do j=1, n
          val=val*(n+i+1-j*2)/j
        end do
        coeff(i+1)=val
      end do
      return
      end subroutine legendre_poly

!----------------------
! build(l, m) builds the spherical harmonic function l, m
!----------------------
      subroutine build(this, l, m)
      implicit none
      class(sph), intent(inout) :: this
      integer, intent(in) :: l, m
      integer :: i, absm
      real*8, parameter :: pi=4.d0*datan(1.d0)
      real*8, dimension(:), allocatable :: coeff
      real*8 :: term

      this%l=l
      this%m=m
      absm=abs(m)
      term=1.
      do i=l-absm+1, l+absm
        term=term*i
      end do
      this%base=sqrt((.5*l+.25)/term/pi)
      if (m/=0) this%base=this%base*sqrt(2.)
      call legendre_poly(l, coeff)
      do i=1, absm
        call deriv(coeff)
      end do
      allocate(this%coeff(l+1-absm))
      this%coeff(:)=coeff(1:(l+1-absm))
      return
      end subroutine build

      subroutine disp(this)
      implicit none
      class(sph), intent(inout) :: this
      integer :: i

      write(*,*)'Y',this%l,this%m,'= ...'
      write(*,*)this%base,'*'
      if (this%m/=0) then
        write(*,*)'sin'
        write(*,*)'^',abs(this%m)
        write(*,*)'(t)*'
      end if
      write(*,*)'(',this%coeff(1)
      do i=2, size(this%coeff)
        write(*,*)'+',this%coeff(i),'*cos^',i-1,'(t)'
      end do
      write(*,*)')'
      if (this%m<0) write(*,*)'*sin(',abs(this%m),'*p)'
      if (this%m>0) write(*,*)'*cos(',abs(this%m),'*p)'
      return
      end subroutine disp

!----------------------
! value(t, p) value at a given point (in spherical coordinates theta, phi)
!----------------------
      function value(this, t, p)
      implicit none
      class(sph), intent(in) :: this
      real*8, intent(in) :: t, p
      real*8 :: value
      integer :: i

      value=0.
      do i=1, size(this%coeff)
        value=value+this%coeff(i)*cos(t)**(i-1)
      end do
      value=value*this%base
      if (this%m/=0) then
        value=value*sin(t)**abs(this%m)
      end if
      if (this%m<0) value=value*sin(-this%m*p)
      if (this%m>0) value=value*cos(this%m*p)
      return
      end function value

      end module sph_mod

!--------------------------------------------------------------------
!--------------------------------------------------------------------


      subroutine src_dist_config
      use resolution, only : im,jm
      USE DOMAIN_DECOMP_atm, ONLY : GRID, getdomainbounds, am_i_root
      use tracer_com, only: xyztr, ntm_sph, ntm_reg
      use dictionary_mod, only: get_param, is_set_param
      implicit none
      include 'netcdf.inc'
      interface
        subroutine src_dist_config_sph(l, m, n, arr)
        integer, intent(in) :: l, m, n
        real*8, dimension(:, :, :), allocatable, intent(inout) :: arr
        end subroutine src_dist_config_sph
      end interface

      integer :: status,fid,vid,did,srt(3),cnt(3)
c
      character(len=80) :: xyzfile='src_dist_cfg'
      integer, dimension(im,jm) :: regions
      real*8, dimension(:,:,:), allocatable :: xyztr_sph
      integer :: i,j,n,l,m,ntm,n_sph
      INTEGER :: J_0, J_1
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL getdomainbounds(grid, J_STRT=J_0, J_STOP=J_1)
      ntm_sph = 1  ! always keep trname(1)='Water   '
      ntm_reg = 0
      status = nf_open(trim(xyzfile),nf_nowrite,fid)
      if(status.ne.nf_noerr) then
        if(am_i_root())
     &       write(6,*) 'input file ',trim(xyzfile),' not found'
        call stop_model('missing input file in src_dist_config',255)
      endif
      status = nf_inq_dimid(fid,'tracer',did)
      if(status.eq.nf_noerr) then
        status = nf_inq_dimlen(fid,did,ntm_sph)
      endif
      status = nf_inq_varid(fid,'regions',vid)
      if(status.eq.nf_noerr) then
        status = nf_get_var_int(fid,vid,regions)
        ntm_reg = maxval(regions)
        do n=1,ntm_reg
          if(.not.any(regions.eq.n)) then
            if(am_i_root()) write(6,*) 'missing region ',n
            call stop_model('missing region in src_dist_config',255)
          endif
        enddo
      endif
      n_sph=0
      if (is_set_param('src_dist_sph'))
     &   call get_param('src_dist_sph', n_sph)
      ntm=ntm_sph+ntm_reg+n_sph**2
      allocate(xyztr(ntm,IM,J_0:J_1))

c read spherical harmonic basis functions
      status = nf_inq_varid(fid,'xyztr',vid)
      if(status.eq.nf_noerr .and. ntm_sph.gt.1) then
        srt(1:3) = (/ 1, 1, j_0 /)
        cnt(1:3) = (/ ntm_sph, im, 1+j_1-j_0 /)
        if(ntm_reg.gt.0) then
          allocate(xyztr_sph(ntm_sph,IM,J_0:J_1))
          status = nf_get_vara_double(fid,vid,srt,cnt,xyztr_sph)
        else
          status = nf_get_vara_double(fid,vid,srt,cnt,xyztr)
        endif
      else  ! if no xyztr in file, retain normal water as a tracer
        xyztr(1,:,:) = 1d0
      endif
      status = nf_close(fid)
c copy spherical funcs and region tracers into xyztr
      do j=j_0,j_1
        do i=1,im
          if(ntm_sph.gt.1 .and. ntm_reg.gt.0) then
            xyztr(1:ntm_sph,i,j)=xyztr_sph(1:ntm_sph,i,j)
          endif
          if(ntm_reg.gt.0) then
            xyztr(ntm_sph+1:ntm,i,j) = 0d0
            xyztr(ntm_sph+regions(i,j),i,j) = 1d0
          endif
        enddo
      enddo
      do l=0, n_sph-1
        call src_dist_config_sph(l, 0,
     &                  ntm_sph+ntm_reg+l**2+1, xyztr)
        do m=1, l
          call src_dist_config_sph(l, m,
     &                  ntm_sph+ntm_reg+l**2+m*2, xyztr)
          call src_dist_config_sph(l, -m,
     &                  ntm_sph+ntm_reg+l**2+m*2+1, xyztr)
        end do
      end do
      return
      end subroutine src_dist_config

      subroutine src_dist_config_sph(l, m, n, arr)
      use sph_mod, only: sph
      use geom, only: lat_dg,lon_dg
      use constant, only: pi
      use domain_decomp_atm, only : grid, getdomainbounds
      implicit none
      integer, intent(in) :: l, m, n
      real*8, dimension(:, :, :), allocatable, intent(inout) :: arr
      integer :: i, j, i0, i1, j0, j1
      real*8, parameter :: factor=sqrt(4.*pi)
      type(sph) :: gen

      call getdomainbounds(grid,
     &                 i_strt=i0, i_stop=i1, j_strt=j0, j_stop=j1)
      call gen%build(l, m)
      do i=i0, i1
        do j=j0, j1
          arr(n,i,j)=gen%value((90.-lat_dg(j,1))*pi/180.,
     &                                 lon_dg(i,1)*pi/180.)*factor
        end do
      end do

      return
      end subroutine src_dist_config_sph

      subroutine init_src_dist
      use resolution, only : im
      use tracer_com, only: ntm, xyztr, ntm_sph, ntm_reg
      use domain_decomp_atm, only : grid, getdomainbounds
      use fluxes, only: focean, asflx
      use oldtracer_mod, only: src_dist_index
      implicit none
      integer :: i, j, n, j_0, j_1

      if (allocated(xyztr)) then
        call getdomainbounds(grid, j_strt=j_0, j_stop=j_1)
        do j=j_0, j_1
          do i=1, im
            do n=1, ntm_reg
              xyztr(n+ntm_sph, i, j)=xyztr(n+ntm_sph, i, j)*focean(i, j)
            enddo
          enddo
        enddo
        do i=1, size(asflx)
          do n=1, ntm
            if (src_dist_index(n)>0) asflx(i)%gtracer(n, :, j_0:j_1)=
     &                             xyztr(src_dist_index(n), :, j_0:j_1)
          enddo
        enddo
      endif
      end subroutine init_src_dist
