#include "rundeck_opts.h"
!------------------------------------------------------------------------------
      subroutine setDefaultSpec(n, pTracer)
!------------------------------------------------------------------------------
      use Dictionary_mod, only: sync_param
      use RunTimeControls_mod, only: tracers_amp
      use RunTimeControls_mod, only: tracers_tomas
      use RunTimeControls_mod, only: tracers_aerosols_vbs
      use OldTracer_mod, only: trName, do_fire, do_aircraft
      use OldTracer_mod, only: set_do_fire, set_do_aircraft
      use OldTracer_mod, only: set_first_aircraft, first_aircraft
      use OldTracer_mod, only: nBBsources, set_nBBsources
      use DOMAIN_DECOMP_ATM, only: am_i_root
      use TRACER_COM, only: tracers
      use TRACER_COM, only: sect_name
      use TRACER_COM, only: set_ntsurfsrc, ntsurfsrc
      use Tracer_mod, only: ntsurfsrcmax
      use TRACER_COM, only: num_sectors
      use Tracer_mod, only: Tracer
      use Tracer_mod, only: findSurfaceSources
      use Tracer_mod, only: addSurfaceSource
      use SystemTools, only: stLinkStatus
#ifdef TRACERS_SPECIAL_Shindell
      use TRCHEM_Shindell_COM, only: use_rad_ch4
#endif
      implicit none

      integer, intent(in) :: n
      class (Tracer), pointer :: pTracer

      logical :: checkSourceName
      integer :: val, linkstatus

      call pTracer%insert('ntSurfSrc', 0)

!     The following section will check for rundeck file of
!     the form: trname_01, trname_02... and thereby define
!     the ntsurfsrc(n). If those files exist it reads an
!     80 char header to get information including the
!     source name (ssame-->{sname,lname,etc.}. ntsurfsrc(n)
!     get set to zero if those files aren't found:
!     (I can enclose this in an ifdef if it causes problems
!     for people). findSurfaceSources routine also assigns
!     sources to sectors, if desired:
!     general case:

      if (tracers_amp .or. tracers_tomas .or. tracers_aerosols_vbs) then
         checkSourceName = .false.
      else if (trname(n) == 'codirect') then 
         checkSourceName = .false.
      else
         checkSourceName = .true.
      end if

      call findSurfaceSources(pTracer, checkSourceName, 
     &     sect_name(1:num_sectors))

!     Next, check whether tracers have 3D aircraft source files/dirs:
      call stLinkStatus(trim(trname(n)//'_AIRC'),linkstatus)
      select case(linkstatus)
      case(1,2) ! TODO: no hardcoded integers
        call set_do_aircraft(n, .true.)
        call set_first_aircraft(n, .true.)
      end select

#ifdef DYNAMIC_BIOMASS_BURNING
!-------------------------------------------------------------------------------
!     allow some tracers to have biomass burning based on fire model:
!-------------------------------------------------------------------------------
        select case (trname(n))
          case('NOx','CO','Alkenes','Paraffin','BCB','OCB','NH3','SO2',
#ifdef TRACERS_AMP
     &         'M_BC1_BC','M_OCC_OC','M_ACC_SU','M_AKK_SU',
#endif
#ifdef TRACERS_dCO
     &         'd13Calke', 'd13CPAR',
     &         'dC17O', 'dC18O', 'd13CO',
#endif  /* TRACERS_dCO */
     &         'vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &         'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6'
#ifdef TRACERS_TOMAS
     &         ,'AECOB_01','AOCOB_01' !BCB and OCB hygroscopities? Need to put emission into OB and IL.
#endif
     &         )
            call set_do_fire(n, .true.)
#ifdef TRACERS_SPECIAL_Shindell
          case('CH4') ! in here to avoid potential Lerner tracers conflict
            if(use_rad_ch4==0) call set_do_fire(n, .true.)
#endif
        end select
#endif /* DYNAMIC_BIOMASS_BURNING */

!-------------------------------------------------------------------------------
!     allow some tracers to have biomass burning sources that mix over
!     PBL layers (these become 3D sources no longer within ntsurfsrc(n)):
!-------------------------------------------------------------------------------
        val = nBBsources(n)
        call sync_param(trim(trname(n))//"_nBBsources",val)
        call set_nBBsources(n, val)
        if(nBBsources(n)>0)then
          if(do_fire(n))then
            if(am_i_root())write(6,*)
     &           'nBBsource>0 for ',trim(trname(n)),' do_fire=t'
            call stop_model('nBBsource do_fire conflict',13)
          else
            call set_ntsurfsrc(n, ntsurfsrc(n)-nBBsources(n))
          end if
        end if
        if(do_fire(n) .and.  (ntsurfsrc(n)+1 > ntsurfsrcmax))then
          write(6,*)trname(n),'ntsurfsrc+1 > max of ',ntsurfsrcmax
          call stop_model('do_fire+ntsurfsrc too large',13)
        end if
        if(ntsurfsrc(n)+nBBsources(n) > ntsurfsrcmax)then
          write(6,*)trname(n),'ntsurfsrc+nBBsources > max of ',
     &         ntsurfsrcmax
          call stop_model('ntsurfsrc+nBBsources too large',13)
        end if

!-------------------------------------------------------------------------------
! add surface sources to tracers that don't follow the TRACERNAME_XX convention
!-------------------------------------------------------------------------------
        select case (trname(n))

#ifndef TRACERS_AEROSOLS_SOA
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
        (defined TRACERS_TOMAS)
        case ('OCII', 'M_OCC_OC', 'SOAgas') ! this handles OCT_src (terpene source)
          call addSurfaceSource(pTracer, "Terpenes_src",
     &                                   "Terpenes source")
#endif
#endif  /* TRACERS_AEROSOLS_SOA */

#ifdef TRACERS_SPECIAL_Shindell
        case('CH4') ! in case use_rad_ch4/=0 but the rundeck lists CH4_XX files
          if (use_rad_ch4/=0) call set_ntsurfsrc(n,0)
#endif

#ifdef TRACERS_SPECIAL_Lerner
        case ('N2O')
          call addSurfaceSource(pTracer, "overwrite_at_surface",
     *          "Overwrite")
        case ('CFC11', 'Rn222')
          call addSurfaceSource(pTracer, "surface_src", "Surface Src")
        case ('CH4')
          call addSurfaceSource(pTracer, "animal_src", "Animals")
          call addSurfaceSource(pTracer, "coal_mine_src", "Coal Mines")
          call addSurfaceSource(pTracer, "gas_leak_src", "Gas Leaks")
          call addSurfaceSource(pTracer, "gas_vent_src", "Gas Venting")
          call addSurfaceSource(pTracer, "city_dump_src", "City Dumps")
          call addSurfaceSource(pTracer, "soil_sink", "Soil Sink")
          call addSurfaceSource(pTracer, "termite_src", "Termites")
          call addSurfaceSource(pTracer, "coal_combustion_src", 
     *         "Coal Combustion")
          call addSurfaceSource(pTracer, "ocean_src", "Ocean Src")
          call addSurfaceSource(pTracer, "lake_src", "Lake Src")
          call addSurfaceSource(pTracer, "misc_ground_src", "Misc Src")
          call addSurfaceSource(pTracer, "biomass_src", "Biomass")
          call addSurfaceSource(pTracer, "rice_src", "Rice")
          call addSurfaceSource(pTracer, "wetlands_tundra_src",
     *         "WetlandsTundra")
        case ('O3')
          call addSurfaceSource(pTracer, "deposition_sink",
     *          "Deposition")
        case ('CO2')
          call addSurfaceSource(pTracer, "fossil_fuel_src",
     *          "Fossil Fuel")
          call addSurfaceSource(pTracer, "fertilization_sink",
     *         "Fertilization")
          call addSurfaceSource(pTracer, "north_forest_regrowth_src",
     *         "Forest Regrow")
          call addSurfaceSource(pTracer, "land_use_modification",
     *         "Land Use")
          call addSurfaceSource(pTracer, "ecosystem_exchange",
     *         "Ecosystem Exch")
          call addSurfaceSource(pTracer, "ocean_exchange", "Ocean Exch")
        case ('14CO2')
          call addSurfaceSource(pTracer, "surface_sink", "Surface Sink")
#endif  /* TRACERS_SPECIAL_Lerner */

#ifdef TRACERS_PASSIVE
      case ('SF6', 'SF6_c', 'nh5', 'nh50', 'aoanh', 'aoa',
     &     'e90', 'st8025', 'tape_rec', 'nh15')

          call addSurfaceSource(pTracer, "surface_src", "Surface Src")
#endif  /* TRACERS_PASSIVE */

#ifdef TRACERS_TOMAS
        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05',
     &        'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10',
     &        'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15')
          call addSurfaceSource(pTracer, "SO4_src", "SO4 source")
          call addSurfaceSource(pTracer, "EC_src", "EC source")
          call addSurfaceSource(pTracer, "OC_src", "OC source")
#endif  /* TRACERS_TOMAS */

        end select

      end subroutine setDefaultSpec

!------------------------------------------------------------------------------
      subroutine initTracerMetadata()
!------------------------------------------------------------------------------
      use TRACER_COM, only: COUPLED_CHEM
      use TRACER_COM, only: ex_volc_num
      use TRACER_COM, only: ex_volc_jday
      use TRACER_COM, only: ex_volc_year
      use TRACER_COM, only: ex_volc_lat
      use TRACER_COM, only: ex_volc_lon
      use TRACER_COM, only: ex_volc_bot
      use TRACER_COM, only: ex_volc_top
      use TRACER_COM, only: ex_volc_SO2
      use TRACER_COM, only: ex_volc_H2O
      use Dictionary_mod, only: set_param, sync_param
      use RunTimeControls_mod, only: tracers_special_shindell
      use RunTimeControls_mod, only: tracers_terp
      use RunTimeControls_mod, only: shindell_strat_extra
      use RunTimeControls_mod, only: accmip_like_diags
      use RunTimeControls_mod, only: tracers_drydep
      use RunTimeControls_mod, only: tracers_tomas
      use RunTimeControls_mod, only: tracers_water
      use RunTimeControls_mod, only: tracers_special_o18
      use RunTimeControls_mod, only: tracers_gasexch_ocean_cfc
      use RunTimeControls_mod, only: tracers_gasexch_ocean_co2
      use RunTimeControls_mod, only: tracers_gasexch_land_co2
      use RunTimeControls_mod, only: tracers_gasexch_gcc
      use RunTimeControls_mod, only: tracers_special_lerner
      use RunTimeControls_mod, only: tracers_passive
      use RunTimeControls_mod, only: tracers_aerosols_koch
      use RunTimeControls_mod, only: tracers_aerosols_seasalt
      use RunTimeControls_mod, only: tracers_aerosols_ocean
      use RunTimeControls_mod, only: tracers_nitrate
      use RunTimeControls_mod, only: tracers_dust
      use RunTimeControls_mod, only: tracers_dust_silt4
      use RunTimeControls_mod, only: tracers_dust_silt5
      use RunTimeControls_mod, only: tracers_hetchem
      use RunTimeControls_mod, only: tracers_cosmo
      use RunTimeControls_mod, only: tracers_radon
      use RunTimeControls_mod, only: tracers_minerals
      use RunTimeControls_mod, only: tracers_on
      use RunTimeControls_mod, only: tracers_air
      use RunTimeControls_mod, only: tracers_amp
      use OldTracer_mod, only: HSTAR
      use OldTracer_mod, only: F0
      use OldTracer_mod, only: nGas
      use OldTracer_mod, only: nPart
      use OldTracer_mod, only: nWater
      use OldTracer_mod, only: tr_wd_type
      use OldTracer_mod, only: tr_mm
      use OldTracer_mod, only: tr_rkd
      use OldTracer_mod, only: trname
      use OldTracer_mod, only: initializeOldTracers
      use OldTracer_mod, only: set_needtrs
      use OldTracer_mod, only: set_mass2vol 
      use OldTracer_mod, only: set_vol2mass
      use OldTracer_mod, only: set_to_conc
      use OldTracer_mod, only: set_dowetdep
      use OldTracer_mod, only: set_dodrydep
      use OldTracer_mod, only: set_to_volume_MixRat
      use Tracer_mod, only: Tracer
#ifdef TRACERS_SPECIAL_Lerner
      use LernerTracersMetadata_mod
#endif
#ifdef TRACERS_PASSIVE
      use PassivetracersMetadata_mod
#endif
#ifdef TRACERS_SPECIAL_Shindell
      use ShindellTracersMetadata_mod
      use TRCHEM_Shindell_COM, only: use_rad_ch4
#endif   
#ifdef TRACERS_TOMAS
      use TomasTracersMetadata_mod
#endif    
#ifdef TRACERS_AMP
      use AmpTracersMetadata_mod
#endif   
#ifdef TRACERS_AEROSOLS_Koch
      use KochTracersMetadata_mod
#endif   
#ifdef TRACERS_AEROSOLS_SEASALT
      use SeasaltTracersMetadata_mod
#endif   
#ifdef TRACERS_NITRATE
      use sharedTracersMetadata_mod, only:
     &  NH3_setSpec, NH4_setSpec
#endif
#ifdef TRACERS_RADON
      use sharedTracersMetadata_mod, only: Rn222_setSpec
#endif
#ifdef TRACERS_MINERALS
      use MineralsTracersMetadata_mod
#endif
      use MiscTracersMetadata_mod
      USE CONSTANT, only: mair
      USE TRACER_COM, only: ntm
      USE TRACER_COM, only: tracers

      implicit none
      class (Tracer), pointer :: pTracer
      external setDefaultSpec
      integer :: i
      integer :: ex

      call sync_param( "COUPLED_CHEM", COUPLED_CHEM )
#ifdef TRACERS_SPECIAL_Shindell
      call sync_param( "use_rad_ch4", use_rad_ch4 )
#endif

! call routine to read/set up sectors for emissions:
      call setup_emis_sectors()

! explosive volcano injections based on rundeck parameters
      call sync_param("ex_volc_num", ex_volc_num)
      if (ex_volc_num>0) then
        allocate(ex_volc_jday(ex_volc_num))
        allocate(ex_volc_year(ex_volc_num))
        allocate(ex_volc_lat(ex_volc_num))
        allocate(ex_volc_lon(ex_volc_num))
        allocate(ex_volc_bot(ex_volc_num))
        allocate(ex_volc_top(ex_volc_num))
        allocate(ex_volc_SO2(ex_volc_num))
        allocate(ex_volc_H2O(ex_volc_num))

! set default emissions to zero
        ex_volc_SO2(:)=0.d0
        ex_volc_H2O(:)=0.d0

        call sync_param("ex_volc_jday", ex_volc_jday, ex_volc_num)
        call sync_param("ex_volc_year", ex_volc_year, ex_volc_num)
        call sync_param("ex_volc_lat", ex_volc_lat, ex_volc_num)
        call sync_param("ex_volc_lon", ex_volc_lon, ex_volc_num)
        call sync_param("ex_volc_bot", ex_volc_bot, ex_volc_num)
        call sync_param("ex_volc_top", ex_volc_top, ex_volc_num)
        call sync_param("ex_volc_SO2", ex_volc_SO2, ex_volc_num)
        call sync_param("ex_volc_H2O", ex_volc_H2O, ex_volc_num)

        do ex=1,ex_volc_num
          if (ex_volc_jday(ex)<=0 .or. ex_volc_jday(ex)>366)
     &      call stop_model('ex_volc_jday(ex) out of bounds', 255)
          if (ex_volc_year(ex)<=0)
     &      call stop_model('ex_volc_year(ex)<=0', 255)
          if (ex_volc_lat(ex)<=-90.d0 .or. ex_volc_lat(ex)>90.d0)
     &      call stop_model('ex_volc_lat(ex) out of bounds', 255)
          if (ex_volc_lon(ex)<=-180.d0 .or. ex_volc_lon(ex)>180.d0)
     &      call stop_model('ex_volc_lon(ex) out of bounds', 255)
          if (ex_volc_bot(ex)>=ex_volc_top(ex))
     &      call stop_model('ex_volc_bot(ex)>=ex_volc_top(ex)', 255)
          if (ex_volc_SO2(ex)<0.d0)
     &      call stop_model('ex_volc_SO2(ex)<0.d0', 255)
          if (ex_volc_H2O(ex)<0.d0)
     &      call stop_model('ex_volc_H2O(ex)<0.d0', 255)
        enddo
      endif

      call initializeOldTracers(tracers, setDefaultSpec)

! ***  BEGIN TRACER METADATA INITIALIZATION

#ifdef TRACERS_SPECIAL_Shindell
        if (tracers_special_shindell) then
          call SHINDELL_InitMetadata(pTracer)
        end if
#endif

#ifdef TRACERS_SPECIAL_Lerner
      if (tracers_special_lerner) then
        call Lerner_InitMetadata(pTracer)
        if (tracers_special_shindell) 
     &    call stop_model('contradictory tracer specs')
      end if
#endif

#ifdef TRACERS_PASSIVE
      if (tracers_passive) then
        call Passive_InitMetadata(pTracer)
      end if
#endif

      if ((.not. tracers_amp) .and. tracers_water) then
        call  Water_setSpec('Water')
      end if
     
#ifdef TRACERS_SPECIAL_O18
        if (tracers_special_o18) then
          call H2O18_setSpec('H2O18')
          call HDO_setSpec('HDO')
#ifdef TRACERS_WISO_O17
          call H2O17_setSpec('H2O17')
#endif
        end if
#endif

      if (tracers_gasexch_ocean_co2 .or. tracers_gasexch_land_co2
     &    .or. tracers_gasexch_gcc) then
        call  CO2n_setSpec('CO2n')
      end if
      
      if (tracers_gasexch_ocean_cfc) then
        call  CFCn_setSpec('CFCn')
      end if

#ifdef TRACERS_AEROSOLS_SEASALT
      if (tracers_aerosols_seasalt) then
        call Seasalt_InitMetadata(pTracer)
      end if
#endif

#ifdef TRACERS_AEROSOLS_Koch
      if (tracers_aerosols_koch) then
        call KOCH_InitMetadata(pTracer)
      end if
#endif

      if (tracers_aerosols_ocean) then
        call  OCocean_setSpec('OCocean') !Insoluble oceanic organic mass
      end if

      if (tracers_dust) then
        call  clay_setSpec('Clay')
        call  Silt1_setSpec('Silt1')
        call  Silt2_setSpec('Silt2')
        call  Silt3_setSpec('Silt3')
        if (tracers_dust_Silt4) call  Silt4_setSpec('Silt4')
        if (tracers_dust_Silt5) call  Silt5_setSpec('Silt5')
      end if

#ifdef TRACERS_NITRATE
      if (tracers_nitrate) then
        call  NH3_setSpec('NH3')
        call  NH4_setSpec('NH4')
        call  NO3p_setSpec('NO3p')
      end if
#endif

      if (tracers_hetchem) then
        call  SO4_d1_setSpec('SO4_d1')
        call  SO4_d2_setSpec('SO4_d2')
        call  SO4_d3_setSpec('SO4_d3')
        if (tracers_nitrate) then
          call  N_d1_setSpec('N_d1')
          call  N_d2_setSpec('N_d2')
          call  N_d3_setSpec('N_d3')
        end if
      end if

      if (tracers_cosmo) then
#ifdef TRACERS_RADON
        if (tracers_radon) then
          call  Pb210_setSpec('Pb210')
        end if
#endif
        call  Be7_setSpec('Be7')
        call  Be10_setSpec('Be10')
#ifdef TRACERS_RADON
        if (tracers_radon) then
          call  Rn222_setSpec('Rn222') ! duplicate with Lerner
        end if
#endif
        if (tracers_special_lerner) 
     &       call stop_model('contradictory tracer specs')
      end if

#ifdef TRACERS_MINERALS
       if (tracers_minerals) then
         call Minerals_InitMetadata(pTracer)
       end if
#endif

      if (tracers_air .or. accmip_like_diags) then
        call  air_setSpec('Air')
      end if

#ifdef TRACERS_AMP
      if (tracers_amp) then
        call AMP_InitMetadata(pTracer)
      end if
#endif

#ifdef TRACERS_TOMAS
      if (tracers_tomas) then
        call TOMAS_InitMetadata(pTracer)
      end if
#endif
      call init_source_distrib

! ***  END TRACER METADATA INITIALIZATION

      ! Generic tracer work
      ! All tracers must have been declared before reaching this point!!!
      ntm = tracers%size()

      call set_param("NTM",NTM,'o')
      call set_param("TRNAME",trName(),ntm,'o')
      call printTracerNames(trName())

      ! Generic tracer work
      do i = 1, ntm
        if (tracers_water) then
!**** Tracers that are soluble or are scavenged or are water => wet dep
          if (tr_wd_type(i).eq.nWater.or.tr_wd_type(i) .EQ. nPART .or.
     &      tr_RKD(i).gt.0) then
            call set_dowetdep(i, .true.)
          end if
        end if
        if (tracers_drydep) then
!**** If tracers are particles or have non-zero HSTAR or F0 do dry dep:
!**** Any tracers that dry deposits needs the surface concentration:
          if(HSTAR(i).GT.0..OR.F0(i).GT.0..OR.tr_wd_type(i).eq.nPART) 
     &      then
            call set_dodrydep(i, .true.)
            call set_needtrs(i, .true.)
            if (tracers_water) then
              if (tr_wd_type(i).eq.nWATER) call stop_model
     &         ('A water tracer should not undergo dry deposition.',255)
            end if
          end if
        end if

        if (tracers_on) then
!**** Define the conversion from mass to volume units here
          call set_mass2vol(i, mair/tr_mm(i))
          call set_vol2mass(i, tr_mm(i)/mair)
          call set_to_conc(i, 0)
          if (tracers_special_shindell) then
!**** Aerosol tracer output should be mass mixing ratio
            select case (tr_wd_TYPE(i))
            case (nGAS)
              call set_to_volume_MixRat(i, 1) !gas output to volume mixing ratio
            case (nPART, nWATER)
              call set_to_volume_MixRat(i, 0) ! aerosol/water output to mass mixing ratio
            case default
              call set_to_volume_MixRat(i, 0) !default output to mass mixing ratio
            end select
          end if
          if (tracers_gasexch_ocean_co2 .or. tracers_gasexch_land_co2) 
     &      then
            call set_to_volume_MixRat(i, 1) !gas output to volume mixing ratio
          end if
       end if
      end do

      contains

      subroutine printTracerNames(tracerNames)
      use domain_decomp_atm, only: am_i_root

      character(len=*) :: tracerNames(:)
      integer :: i
      
      if (am_i_root()) then
        do i = 1, size(tracerNames)
          write(6,*) 'TRACER',i,trim(tracerNames(i))
        end do
      end if
      
      end subroutine printTracerNames

      end subroutine initTracerMetadata


!------------------------------------------------------------------------------
      subroutine init_source_distrib
      use dictionary_mod, only : is_set_param, get_param
      use oldtracer_mod, only: oldaddtracer, set_t_qlimit, findtracer,
     &  set_src_dist_base, set_src_dist_index
      use tracer_com, only: xyztr
      implicit none
      character(len=1024) :: list
      integer :: str_pos, i, nt, nt_orig, ndigits
      character(len=10) :: name, basename
      character(len=8) :: fm

      if (is_set_param('src_dist_tr')) then
        call src_dist_config
        call get_param('src_dist_tr', list)
        list=adjustl(list)
        do while (len_trim(list).gt.0)
          str_pos=index(list, ' ')
          name=list(1:str_pos-1)
          basename=name
          nt_orig=findtracer(basename)
          call set_t_qlimit(nt_orig, .false.)
          call set_src_dist_base(nt_orig, nt_orig)
          call set_src_dist_index(nt_orig, 1)
          ndigits=int(log10(size(xyztr, 1)*1d0))+1
          if (ndigits>7) call stop_model('too many digits',255)
          write(fm,'(a,i1.1,a,i1.1,a)') '(a,i',ndigits,'.',ndigits,')'
          do i=2, size(xyztr, 1)
            write(name, fm) trim(basename(1:8-ndigits)), i
            nt=oldaddtracer(name, basename)
            call set_src_dist_index(nt, i)
          end do
          list=adjustl(list(str_pos:))
        end do
      endif

      return
      end subroutine init_source_distrib
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      subroutine laterInitTracerMetadata()
!------------------------------------------------------------------------------
      USE MODEL_COM, only: itime,master_yr
      use OldTracer_mod, only: itime_tr0
      use OldTracer_mod, only: set_itime_tr0
      USE TRACER_COM, only: NTM, tracers, syncProperty
      use TRACER_COM, only: nc_emis_use_ppm_interp
      use Dictionary_mod, only: sync_param,is_set_param,get_param
      use RAD_COM, only: diag_fc
#ifdef TRACERS_SPECIAL_O18
      use tracer_com, only: supsatfac
#endif
#ifdef TRACERS_WATER
      use TRDIAG_com, only: to_per_mil
#endif
#ifdef TRACERS_SPECIAL_Shindell
      use TRCHEM_Shindell_COM, only: tune_NOx
      use TRCHEM_Shindell_COM, only: tune_BVOC
#endif  /* TRACERS_SPECIAL_Shindell */
      use TRACER_COM, only: tune_BBsources
#ifdef TRACERS_AEROSOLS_Koch
      use aerosol_sources, only: tune_DMS
#endif  /* TRACERS_AEROSOLS_Koch */
#ifdef TRACERS_AEROSOLS_SEASALT
      use tracers_seasalt, only: tune_ss1, tune_ss2
#endif  /* TRACERS_AEROSOLS_SEASALT */
      use TRDIAG_COM, only: diag_rad,diag_aod_3d,save_dry_aod
      use TRACER_COM, only: ntm ! should be available by this procedure call
#ifdef TRACERS_WATER
#ifdef TRDIAG_WETDEPO
      USE CLOUDS, ONLY : diag_wetdep
#endif
#endif /* TRACERS_WATER */
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM,only:LCOalt,PCOalt,
     &     CH4altINT,CH4altINX,LCH4alt,PCH4alt,
     &     CH4altX,CH4altT,scale_ch4_IC_file,
     &     OxICIN,OxIC,OxICINL,OxICL,
     &     fix_CH4_chemistry,which_trop,PI_run,PIratio_N,PIratio_CO_T,
     &     PIratio_CO_S,PIratio_other,allowSomeChemReinit,
     &     CH4ICIN,CH4ICX,CH4ICINL,CH4ICL,
     &     COICIN,COIC,COICINL,COICL,Lmax_rad_O3,Lmax_rad_CH4
     &     ,BrOxaltIN,ClOxaltIN,ClONO2altIN,HClaltIN,BrOxalt,
     &     ClOxalt,ClONO2alt,HClalt,N2OICIN,N2OICX,N2OICINL,N2OICL,
     &     CFCICIN,CFCIC,CFCICINL,CFCICL,PIratio_N2O,PIratio_CFC,
     &     use_rad_n2o,use_rad_cfc,cfc_rad95,PltOx,Tpsc_offset_N,
     &     Tpsc_offset_S,windowN2Ocorr,windowO2corr,
     &     reg1Power_SpherO2andN2Ocorr,reg1TopPres_SpherO2andN2Ocorr,
     &     reg2Power_SpherO2andN2Ocorr,reg2TopPres_SpherO2andN2Ocorr,
     &     reg3Power_SpherO2andN2Ocorr,reg3TopPres_SpherO2andN2Ocorr,
     &     reg4Power_SpherO2andN2Ocorr
      use photolysis, only: rad_FL
#ifdef INTERACTIVE_WETLANDS_CH4
      USE TRACER_SOURCES, only:int_wet_dist,topo_lim,sat_lim,gw_ulim,
     &  gw_llim,sw_lim,exclude_us_eu,nn_or_zon,ice_age,nday_ch4,
     &  max_days,ns_wet,nra_ch4
#endif
#ifdef BIOGENIC_EMISSIONS
      use biogenic_emis, only: base_isopreneX
#endif
      use RAD_COM, only: O3_yr
      use TRCHEM_Shindell_COM, only: NOx_yr
      use TRCHEM_Shindell_COM, only: CO_yr
      use TRCHEM_Shindell_COM, only: VOC_yr
#endif /* TRACERS_SPECIAL_Shindell */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)  || (defined TRACERS_AEROSOLS_SEASALT)
      use TRACER_COM, only: aer_int_yr
      use TRACER_COM, only: SO2_int_yr
      use TRACER_COM, only: NH3_int_yr
      use TRACER_COM, only: BC_int_yr
      use TRACER_COM, only: OC_int_yr
#endif
#ifdef TRACERS_AMP
      USE AMP_AEROSOL, only: AMP_RAD_KEY
#endif
#if (defined TRACERS_COSMO)
      USE COSMO_SOURCES, only: be7_src_param
#endif
#ifdef TRACERS_AEROSOLS_VBS
      USE AEROSOL_SOURCES, only: VBSemifact
      USE TRACERS_VBS, only: vbs_tr
#endif  /* TRACERS_AEROSOLS_VBS */
      USE TRACER_COM, only: no_emis_over_ice
#ifdef TRACERS_MINERALS
      use trdust_mod, only: frIronOxideInAggregate,
     &     noAggregateByTotalFeox
#endif
      use Model_com, only: itime
      implicit none
      integer :: n

! call routine to read/set up regions for emissions:
      call setup_emis_sectors_regions()

C**** 
C**** Set some documentary parameters in the database
C**** 
      do n = 1, NTM
        call set_itime_tr0(n, itime)
      end do
      call syncProperty(tracers, "itime_tr0", set_itime_tr0,itime_tr0())

C**** Synchronise tracer related parameters from rundeck

      call sync_param("nc_emis_use_ppm_interp",nc_emis_use_ppm_interp)
 
#ifdef TRACERS_WATER
C**** Decide on water tracer conc. units from rundeck if it exists
      call sync_param("to_per_mil",to_per_mil,ntm)
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call sync_param("tune_NOx",tune_NOx)
      call sync_param("tune_BVOC",tune_BVOC)
      call get_param("O3_yr", O3_yr, default=master_yr) ! duplicate of RAD_DRV
      if (is_set_param("NOx_yr")) then
        call get_param("NOx_yr",NOx_yr)
      else
        NOx_yr=O3_yr
      endif
      if (is_set_param("CO_yr")) then
        call get_param("CO_yr",CO_yr)
      else
        CO_yr=O3_yr
      endif
      if (is_set_param("VOC_yr")) then
        call get_param("VOC_yr",VOC_yr)
      else
        VOC_yr=O3_yr
      endif
#endif  /* TRACERS_SPECIAL_Shindell */
      call sync_param("tune_BBsources",tune_BBsources)
#ifdef TRACERS_AEROSOLS_Koch
      call sync_param("tune_DMS",tune_DMS)
#endif  /* TRACERS_AEROSOLS_Koch */
#ifdef TRACERS_AEROSOLS_SEASALT
      call sync_param("tune_ss1",tune_ss1)
      call sync_param("tune_ss2",tune_ss2)
#endif  /* TRACERS_AEROSOLS_SEASALT */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
      (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)
C**** determine year of emissions
      if (is_set_param("aer_int_yr")) then
        call get_param("aer_int_yr",aer_int_yr)
      else
        aer_int_yr=master_yr
      endif
      if (is_set_param("SO2_int_yr")) then
        call get_param("SO2_int_yr",SO2_int_yr)
      else
        SO2_int_yr=aer_int_yr
      endif
      if (is_set_param("NH3_int_yr")) then
        call get_param("NH3_int_yr",NH3_int_yr)
      else
        NH3_int_yr=aer_int_yr
      endif
      if (is_set_param("BC_int_yr")) then
        call get_param("BC_int_yr",BC_int_yr)
      else
        BC_int_yr=aer_int_yr
      endif
      if (is_set_param("OC_int_yr")) then
        call get_param("OC_int_yr",OC_int_yr)
      else
        OC_int_yr=aer_int_yr
      endif
#endif
#ifdef TRACERS_AEROSOLS_VBS
      call sync_param("VBSemifact",VBSemifact,vbs_tr%nbins)
#endif
#ifdef TRACERS_SPECIAL_O18
C**** set super saturation parameter for isotopes if needed
      call sync_param("supsatfac",supsatfac)
#endif
#ifdef TRACERS_ON
      CALL sync_param("diag_rad",diag_rad)
      CALL sync_param("diag_aod_3d",diag_aod_3d)
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
      CALL sync_param("diag_wetdep",diag_wetdep)
#endif
!     not params call sync_param("trans_emis_overr_day",trans_emis_overr_day)
!     not params call sync_param("trans_emis_overr_yr", trans_emis_overr_yr )
#endif /* TRACERS_ON */
#ifdef TRACERS_SPECIAL_Shindell
      call sync_param("allowSomeChemReinit",allowSomeChemReinit)
      call sync_param("which_trop",which_trop)
      if (is_set_param("PI_run")) then
        call get_param("PI_run",PI_run)
      else
        if (master_yr == 1850) then
          PI_run=1
        else
          PI_run=0
        endif
      endif
      call sync_param("PIratio_N",PIratio_N)
      call sync_param("PIratio_CO_T",PIratio_CO_T)
      call sync_param("PIratio_CO_S",PIratio_CO_S)
      call sync_param("PIratio_other",PIratio_other)
      call sync_param("rad_FL",rad_fl)
      call sync_param("Lmax_rad_O3",Lmax_rad_O3)
      call sync_param("Lmax_rad_CH4",Lmax_rad_CH4)
      call sync_param("use_rad_n2o",use_rad_n2o)
      call sync_param("use_rad_cfc",use_rad_cfc)
      call sync_param("PIratio_N2O",PIratio_N2O)
      call sync_param("PIratio_CFC",PIratio_CFC)
      call sync_param("PltOx",PltOx)
      call sync_param("Tpsc_offset_N",Tpsc_offset_N)
      call sync_param("Tpsc_offset_S",Tpsc_offset_S)
      call sync_param("reg1Power_SpherO2andN2Ocorr",
     &                 reg1Power_SpherO2andN2Ocorr)
      call sync_param("reg2Power_SpherO2andN2Ocorr",
     &                 reg2Power_SpherO2andN2Ocorr)
      call sync_param("reg3Power_SpherO2andN2Ocorr",
     &                 reg3Power_SpherO2andN2Ocorr)
      call sync_param("reg4Power_SpherO2andN2Ocorr",
     &                 reg4Power_SpherO2andN2Ocorr)
      call sync_param("reg1TopPres_SpherO2andN2Ocorr",
     &                 reg1TopPres_SpherO2andN2Ocorr)
      call sync_param("reg2TopPres_SpherO2andN2Ocorr",
     &                 reg2TopPres_SpherO2andN2Ocorr)
      call sync_param("reg3TopPres_SpherO2andN2Ocorr",
     &                 reg3TopPres_SpherO2andN2Ocorr)
      call sync_param("windowN2Ocorr",windowN2Ocorr)
      call sync_param("windowO2corr",windowO2corr)
#ifdef BIOGENIC_EMISSIONS
      call sync_param("base_isopreneX",base_isopreneX)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      call sync_param("ice_age",ice_age)
      call sync_param("ns_wet",ns_wet)
      call sync_param("int_wet_dist",int_wet_dist)
      call sync_param("topo_lim",topo_lim)
      call sync_param("sat_lim",sat_lim)
      call sync_param("gw_ulim",gw_ulim)
      call sync_param("gw_llim",gw_llim)
      call sync_param("sw_lim",sw_lim)
      call sync_param("exclude_us_eu",exclude_us_eu)
      call sync_param("nn_or_zon",nn_or_zon)
      do n=1,nra_ch4
        if(nday_ch4(n) > max_days .or. nday_ch4(n) < 1)
     &       call stop_model('nday_ch4 out of range',255)
      end do
#endif

#endif /* TRACERS_SPECIAL_Shindell */
      call sync_param("diag_fc",diag_fc)

#if (defined TRACERS_AMP)
C**** Decide Radiative Mixing Rules - Volume - Core Shell - Maxwell Garnett, default Volume
      call sync_param("AMP_RAD_KEY",AMP_RAD_KEY)
#endif
#ifdef TRACERS_MINERALS
      call sync_param('frIronOxideInAggregate', frIronOxideInAggregate)
      call sync_param('noAggregateByTotalFeox', noAggregateByTotalFeox)
#endif

#if (defined TRACERS_COSMO)
C**** get rundeck parameter for cosmogenic source factor
      call sync_param("be7_src_param", be7_src_param)
#endif
      call sync_param("no_emis_over_ice",no_emis_over_ice)

      end subroutine laterInitTracerMetadata

!------------------------------------------------------------------------------
      subroutine InitTracerMetadataAtmOcnCpler()
!------------------------------------------------------------------------------
      use Dictionary_mod, only: sync_param
#if (defined TRACERS_OCEAN) && !defined(TRACERS_OCEAN_INDEP)
! atmosphere copies atmosphere-declared tracer info to ocean
! so that the ocean can "inherit" it without referencing atm. code
      use ocn_tracer_com, only : tracerlist, ocn_tracer_entry,
     &     n_Water_ocn      => n_Water
      use oldtracer_mod, only: itime_tr0, ntrocn, t_qlimit,
     &            conc_from_fw, trdecay, vol2mass
      use tracer_com, only: n_water
      use trdiag_com, only: to_per_mil
#endif
      USE FLUXES, only : atmocn
      use OldTracer_mod, only: vol2mass
      USE TRACER_COM, only: ntm, gasex_index
      use OldTracer_mod, only: trw0
      implicit none
      integer :: n
#if (defined TRACERS_OCEAN) && !defined(TRACERS_OCEAN_INDEP)
      type(ocn_tracer_entry), pointer :: entry
#endif

#if (defined TRACERS_OCEAN) && !defined(TRACERS_OCEAN_INDEP)
! atmosphere copies atmosphere-declared tracer info to ocean module
! so that the ocean can "inherit" it without referencing atm. code
      n_Water_ocn = n_Water
      do n=1,ntm
        entry=>tracerlist%at(n)
        entry%itime_tr0    = itime_tr0(n)
        entry%ntrocn       = ntrocn(n)
        entry%to_per_mil   = to_per_mil(n)
        entry%t_qlimit     = t_qlimit(n)
        entry%conc_from_fw = conc_from_fw(n) 
        entry%trdecay      = trdecay(n)
        entry%trw0         = trw0(n)
#ifdef TRACERS_WATER
        entry%need_ic      = .true.
#endif
      enddo
#endif

! copy atmosphere-declared tracer info to atm-ocean coupler data
! structure for uses within ocean codes
      allocate(atmocn%trw0(ntm))
      do n=1,ntm
        atmocn%trw0(n) = trw0(n)
      enddo
      allocate(atmocn%vol2mass(gasex_index%getsize()))
      do n=1,gasex_index%getsize()
        atmocn%vol2mass(n) = vol2mass(gasex_index%at(n))
      enddo

      end subroutine InitTracerMetadataAtmOcnCpler

!------------------------------------------------------------------------------
      subroutine InitTracerDiagMetadata()
!------------------------------------------------------------------------------
      use TRACER_COM, only: remake_tracer_lists
      implicit none

C**** create tracer lists, needed for some diagnostic output decisions
      call remake_tracer_lists()

C**** Set some diags that are the same regardless
      call set_generic_tracer_diags

C**** Zonal mean/height diags
      call init_jls_diag

C**** lat/lon tracer sources, sinks and specials
      call init_ijts_diag

C**** lat/lon/height tracer specials
      call init_ijlts_diag

C**** Initialize conservation diagnostics
      call init_tracer_cons_diag

      end subroutine InitTracerDiagMetadata
