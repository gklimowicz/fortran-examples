#include "rundeck_opts.h"

#define CHECK_CARBON_CONSERVATION

      module ent_drv
!@sum ent_drv contains variables and routines for vegetation driver
!@auth I. Aleinov, N. Kiang, Y. Kim

      use resolution, only : im,jm
#ifdef HEALY_LM_DIAGS
      use diag_com, only : CROPS_DIAG
#endif
      use ent_mod
      use timestream_mod, only : timestream
      implicit none
      private
      save

      public init_module_ent, update_vegetation_data
      public map_ent2giss !YKIM- temporary hack to use Ent pfts in modelE

      logical :: initialized = .false.
      integer :: crops_yr = 0
      integer :: do_soilresp
      integer :: do_phenology_activegrowth,do_structuralgrowth
      integer :: do_frost_hardiness,do_patchdynamics, do_init_geo
      logical :: do_modis_lai = .false. 

      type(timestream), allocatable :: LAIstream(:)

      contains



      subroutine init_module_ent(iniENT_in, Jday, Jyear, istart)
!@sum initializes vegetation
      use Dictionary_mod
      use ent_com, only : entcells,Cint,Qfol,cnc_ij,excess_C
      !use ent_prescr_veg, only : prescr_calcconst
      use ent_prescribed_drv, only : ent_init_params
      use fluxes, only : focean, FLICE
      use MODEL_COM, only: master_yr
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      use timestream_mod, only : init_stream
      use filemanager, only : file_exists
      integer, intent(in) :: Jday, Jyear
      logical, intent(in) :: iniENT_in
      !---
      integer I_0, I_1, J_0, J_1, i, j
      ! the following are rundeck parameters which need to be
      ! passed to ent (and used there)
      integer year
      integer :: force_init_ent=0
      integer :: ent_io_plain_array=1
      logical iniENT
      integer pft
      integer istart

      if ( initialized ) return
      initialized = .true.

      iniENT = iniENT_in

      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               I_STRT     =I_0,    I_STOP     =I_1)

!!! hack
      call sync_param( "init_ent", force_init_ent)
      if ( force_init_ent == 1 .and. istart <= 8 ) then
        iniENT = .true.
        print *,"!!! Will reinitialize Ent !!!"
      endif

      call sync_param( "ent_io_plain_array", ent_io_plain_array)

      !YKIM - default parameters
      do_soilresp = 1 !true
      do_phenology_activegrowth = 0 !false
      do_structuralgrowth = 0 !false
      do_frost_hardiness = 1 !true
      do_patchdynamics = 0 !false
      do_init_geo = 0 !false

      !--- read rundeck parameters
      call get_param( "crops_yr", crops_yr, default=master_yr )
      !YKIM - add new options for modelE
      call sync_param( "do_soilresp",do_soilresp)
      call sync_param( "do_phenology_activegrowth",
     &                 do_phenology_activegrowth)
      call sync_param( "do_structuralgrowth",do_structuralgrowth)
      call sync_param( "do_frost_hardiness",do_frost_hardiness)
      call sync_param( "do_patchdynamics",do_patchdynamics)
      call sync_param( "do_init_geo",do_init_geo)

      if ( crops_yr .ne. 0 ) then
        year = crops_yr
      else
        year = Jyear
      endif

      if ( file_exists('LAI') ) then
        if ( .not. file_exists('LAIMAX') ) then
          call stop_model("Prescribed LAI requires LAIMAX file",255)
        endif
        if ( .not. file_exists('HITEent') ) then
          call stop_model("Prescribed LAI requires HITEent file",255)
        endif
        do_modis_lai = .true.
        allocate( LAIstream(N_PFT) )
        do pft=1,N_PFT
          call init_stream(grid,LAIstream(pft)
     &         ,'LAI'           ! name of file in rundeck
     &         ,trim(ent_cover_names(pft)) ! netcdf name of the variable to be read
     &         ,0d0,1d30        ! min/max valid data values (irrelevant for linm2m method)
     &         ,'linm2m'        ! time interp method.  For monthly data, daily data created
                                ! by linearly interpolating between month midpoints.
     &         ,Jyear,Jday)
        enddo
      endif

      ! maybe call "ent_initialize" ? , i.e.
      ! call ent_initialize(cond_scheme,vegCO2X_off,crops_yr,nl_soil, etc)
      ! ask Max if OK

      call ent_init_config(
     &     do_soilresp=(do_soilresp==1) !.true.
     &     ,do_phenology_activegrowth=(do_phenology_activegrowth==1) !.true.
     &     ,do_structuralgrowth=(do_structuralgrowth==1) 
     &     ,do_frost_hardiness=(do_frost_hardiness==1) !.true. !.false.
     &     ,do_patchdynamics=(do_patchdynamics==1) !.false.
     &     ,do_init_geo=(do_init_geo==1) !default false
     &     )

      !call prescr_calcconst() ! moved above
      call ent_init_params()

      if (iniENT ) then
      ! initialize ent cells to something meaningful
        do j=J_0,J_1
          do i=I_0,I_1
            !!if (focean(i,j) <= 0) then
            if ( focean(i,j) < 1.d0 ) then
              !print *,"EARTH",i,j,focean(i,j)
              !print *,"ent_construct",i,j
              call ent_cell_construct( entcells(i,j) )
            else
              !print *,"OCEAN",i,j,focean(i,j)
            end if
          enddo
        enddo

        !print *,'init_module_ent printing cells 1'
        !call ent_cell_print(6, entcells)
        !print *,'init_module_ent end printing cells 1'

        call set_vegetation_data( entcells, ! (I_0:I_1,J_0:J_1),
     &       IM, JM, I_0, I_1, J_0, J_1, jday, year, .true. ,
     &       prog_veg = (do_phenology_activegrowth==1))

        ! probably it is ok to initialize these here since 
        ! if .not. iniENT we should read everything from restart file
        Cint(:,:)=0.0127D0      ! internal CO2
        Qfol(:,:)=3.D-6         ! surface mixing ratio
        cnc_ij(:,:) = 0.d0
        excess_C(:,:) = 0.d0

      else ! i.e. *not* iniENT

        ! just in case, do nothing, just set heat capacities
        call ent_prescribe_vegupdateB(entcells)        

      endif
      end subroutine init_module_ent


      subroutine set_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, year, reinitialize, 
     &     prog_veg)
!@sum read standard GISS vegetation BC''s and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
!@+   For prog_veg, added option to calculate vegetation structure given
!@+   geographic input files.
      use DOMAIN_DECOMP_ATM, only : GRID
      use geom, only : lat2d
      use ent_prescribed_drv, only : init_canopy_physical!,prescr_vegdata

      use ent_prescribed_drv, only:
     &     prescr_get_laidata,prescr_veg_albedodata,
     &     prescr_get_hdata,
     &     prescr_calc_canopy_geometry,
!     &     prescr_get_woodydiameter,prescr_get_pop,
!     &     prescr_get_crownrad,
     &     prescr_get_carbonplant,
     &     prescr_get_pft_vars,
!     &     prescr_get_initnm,
!     &     prescr_get_rootprof,prescr_get_soilcolor,
     &     prescr_get_soilpools,prescr_get_soil_C_total,
     &     init_ent_laidata, init_ent_hdata,  prescr_get_ent_plant,
     &     prescr_get_laimaxdata

      use ent_prescribed_drv_geo, only : init_entvegdata_geo
     &     ,init_ent_laimax_geo

      use ghy_com, only : q_ij, qk_ij, dz_ij
      use ghy_com, only : fearth
      use timestream_mod, only : read_stream
      use veg_drv, only : get_vdata, get_soil_c_total, get_cropdata
      type(entcelltype_public), intent(inout) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      logical :: reinitialize
      logical, intent(in), optional :: prog_veg
      !Local variables
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laimaxdata  !cohort
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: dbhdata !Diameter at breast height for woody veg.(cm)
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: craddata !Crown radius (m)
      real*8, dimension(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1) :: cpooldata !Carbon pools in individuals
      integer, dimension(N_COVERTYPES) :: soil_color ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      !real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture1
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpool_ini
      !-----Local---------
      integer i,j,k
      real*8 heat_capacity
      real*8 :: vdata_H(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO,N_COVERTYPES)
      real*8 :: cropdata_H(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)
!      real*8 :: laimaxdata_H(grid%I_STRT_HALO:grid%I_STOP_HALO,
!     &     grid%J_STRT_HALO:grid%J_STOP_HALO,N_COVERTYPES)
      REAL*8 :: soil_C_total_H(N_CASA_LAYERS,
     &     grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)
      integer hemi(I0:I1,J0:J1)

      real*8 :: laidata_h(N_PFT,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)
      integer :: pft
      character*80 :: tilelai, ttt

      ttt = "testing"

      !* Set hemisphere flags.
      where(lat2d(I0:I1,J0:J1) <= 0.)
        hemi(I0:I1,J0:J1) = -1  ! S
      elsewhere
        hemi(I0:I1,J0:J1) = +1  ! N
      end where

      call init_canopy_physical(I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

      if ( do_modis_lai ) then
        do pft=1,N_PFT
          call read_stream(grid,LAIstream(pft),year,jday,
     &         laidata_h(pft,:,:))
        enddo
      endif

      call get_vdata(year, vdata_H, ent_cover_names)
      if ( year == -1 ) then
        do k=1,N_COVERTYPES
          vegdata(k,I0:I1,J0:J1) = vdata_H(I0:I1,J0:J1,k)
        enddo        
      else
        call get_cropdata(year, cropdata_H)
        do k=1,N_COVERTYPES
          vegdata(k,I0:I1,J0:J1) = vdata_H(I0:I1,J0:J1,k)
     &         *(1.d0-cropdata_H(I0:I1,J0:J1))
        enddo
        vegdata(CROPS+COVEROFFSET,I0:I1,J0:J1) = cropdata_H(I0:I1,J0:J1)
#ifdef HEALY_LM_DIAGS
        CROPS_DIAG(I0:I1,J0:J1) = cropdata_H(I0:I1,J0:J1)
#endif
      endif


      !* Canopy structure
      if (.not. present(prog_veg) .or. .not.prog_veg) then
         !* Variables that may be from files or prescribed model values.
        if ( do_modis_lai ) then
          laidata(:,:,:) = 0.d0
          laidata(1:N_PFT,I0:I1,J0:J1)=laidata_h(:,I0:I1,J0:J1)
          call read_height(I0,I1,J0,J1, hdata)
          call read_laimax(I0,I1,J0,J1,laimaxdata)
        else
          call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata)
          call prescr_get_hdata(I0,I1,J0,J1,hdata) !height
          call prescr_get_laimaxdata(I0,I1,J0,J1,laimaxdata)
        endif
         !do k=1,N_COVERTYPES
         !  write(922) ttt, real(hdata(k,:,:),kind=4)
         !enddo
         !do k=1,N_COVERTYPES
         !  write(923) ttt, real(laimaxdata(k,:,:),kind=4)
         !enddo
         
         !* Variables derived from internal allometry. Subject to change.
         call prescr_calc_canopy_geometry(I0,I1,J0,J1
     &        ,hdata,laimaxdata,dbhdata,popdata,craddata)
!call prescr_calc_canopy_geometry_cpools replaces loop below.
!prescr_get_hdata is separated for case of reading heights from file.
!        do j=J0,J1
!          do i=I0,I1
!            call prescr_get_hdata(hdata(:,i,j)) !height
!            call prescr_get_woodydiameter(hdata(:,i,j), dbhdata(:,i,j))
!            call prescr_get_pop(dbhdata(:,i,j), popdata(:,i,j))
!            call prescr_get_crownrad(popdata(:,i,j), craddata(:,i,j))
!          enddo
!        enddo
         call prescr_get_carbonplant(I0,I1,J0,J1,
     &        laidata,hdata,dbhdata,popdata
     &        ,cpooldata)
      else                      !if prog_veg=true
         call init_ent_laimax_geo(IM,JM,I0,I1,J0,J1,laimaxdata) !Off-line routine
         !call get_laimaxdata(laimaxdata_H) !VEG_DRV version of init_ent_laimax_geo
         call init_ent_hdata(IM,JM,I0,I1,J0,J1,hdata) !height
         if (do_init_geo==0) then !This is hack if no LAI init file is avail.
	    call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata) !lai
c            call prescr_get_ent_plant(I0,I1,J0,J1, 
c     &           laidata,hdata,laimaxdata
c     &           ,dbhdata,popdata,craddata,cpooldata)
         else  !Geographic LAI
            call init_ent_laidata(IM,JM,I0,I1,J0,J1,laidata) !lai
         endif
         call init_entvegdata_geo( IM,JM,I0,I1,J0,J1
     i        ,laidata,hdata,laimaxdata
     o        ,popdata ,dbhdata,craddata,cpooldata) 
      endif

      !* Canopy albedo - moved to after specification of veg structure.
      call prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)

      !* Other plant characteristics
      call prescr_get_pft_vars(nmdata,rootprofdata,soil_color)
      !call prescr_get_initnm(nmdata) !nm ! mean canopy nitrogen
      !call prescr_get_rootprof(rootprofdata)
      !call prescr_get_soilcolor(soil_color)

      !* Soil
#ifdef SET_SOILCARBON_GLOBAL_TO_ZERO
      Tpool_ini = 0.d0
#else
      !call prescr_get_soil_C_total(IM,JM,I0,I1,J0,J1,soil_C_total)
      call get_soil_C_total(N_CASA_LAYERS, soil_C_total_H)
      call prescr_get_soilpools(I0,I1,J0,J1,
     &      soil_C_total_H(:,I0:I1,J0:J1), Tpool_ini)
#endif

      !Compute soil textures for upper 30cm
      do j=j0,j1
        do i=i0,i1
          !if ( fearth(i,j) < .01 ) cycle
          soil_texture(:,i,j) = (
     &            q_ij(i,j,:,1)*dz_ij(i,j,1)
     &         +  q_ij(i,j,:,2)*dz_ij(i,j,2)
     &         +  q_ij(i,j,:,3)*(.3d0 - dz_ij(i,j,1) - dz_ij(i,j,2))
     &         ) / .3d0

        enddo
      enddo

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
      call ent_cell_set(entcells, vegdata, popdata, laidata,
     &     hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata, 
     &     soil_color, albedodata, soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini,
     &     reinitialize)

      ! just in case, do nothing, just set heat capacities
      call ent_prescribe_vegupdateB(entcells)

      end subroutine set_vegetation_data



      subroutine update_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, jyear )
!@sum read standard GISS vegetation BC''s and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      use DOMAIN_DECOMP_ATM, only : GRID
      use geom, only : lat2d
      use ent_prescribed_drv, only:
     &     prescr_get_laidata,prescr_veg_albedodata,prescr_get_cropdata
      !use ent_prescr_veg, only: prescr_get_laidata,prescr_veg_albedodata
      use ent_prescr_veg, only : prescr_calc_lai
      use timestream_mod, only : read_stream
      type(entcelltype_public), intent(inout) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, jyear
      !Local variables
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      !real*8, dimension(I0:I1,J0:J1) :: cropsdata
      !-----Local---------
      integer hemi(I0:I1,J0:J1)
      integer i,j
      integer year
      integer, save :: year_old = -1
      real*8 :: cropdata_H(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)

      real*8 :: laidata_h(N_PFT,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)
      integer :: pft
      character*80 :: tilelai

#ifdef CHECK_CARBON_CONSERVATION
      real*8, dimension(I0:I1,J0:J1) :: C_entcell_start, C_entcell
      real*8 :: dC

      call ent_get_exports( entcells,
     &     C_entcell=C_entcell_start
     &     )
#endif
! update crops here 
!          year1 = 1965
      if ( crops_yr .ne. 0 ) then
        year = crops_yr
      else
        year = Jyear
      endif

      if( year .ne. year_old ) then
        !call prescr_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropsdata)
cddd        call get_cropdata(year, cropdata_H)
cddd        call ent_prescribe_vegupdate(entcells,
cddd     &       do_giss_lai=.false.,
cddd     &       cropsdata=cropdata_H(I0:I1,J0:J1) )
        call set_vegetation_data( entcells,
     &       im, jm, i0, i1, j0, j1, jday, year, .false. )
        year_old = year
      endif

            !* Set hemisphere flags.
      where(lat2d(I0:I1,J0:J1) <= 0.)
        hemi(I0:I1,J0:J1) = -1  ! S
      elsewhere
        hemi(I0:I1,J0:J1) = +1  ! N
      end where


      if ( do_modis_lai ) then
        laidata_h(:,:,:) = 0.d0
        do pft=1,N_PFT
          call read_stream(grid,LAIstream(pft),year,jday,
     &         laidata_h(pft,:,:))
        enddo
        
        call ent_prescribe_vegupdate(entcells,hemi,jday,year,
     &       do_giss_phenology=(do_phenology_activegrowth==0), !.false.,
     &       do_giss_albedo= .true.,
     &       do_giss_lai=.false., !.false.,
     &       update_crops=.false.,
     &       laidata=laidata_h(:,I0:I1,J0:J1)
     &       )
      else
        call ent_prescribe_vegupdateC(entcells,hemi,jday,year,
     &       do_giss_phenology=(do_phenology_activegrowth==0), !.false.,
     &       do_giss_albedo= .true.,
     &       do_giss_lai=(do_phenology_activegrowth==0), !.false.,
     &       update_crops=.false. )
      endif
      ! hack to avoid descrepancy with ent_standalone setup
      ! but really should do update as below

#ifdef CHECK_CARBON_CONSERVATION
      call ent_get_exports( entcells,
     &     C_entcell=C_entcell
     &     )
      do j=j0,j1
        do i=i0,i1
          dC = C_entcell(i,j) - C_entcell_start(i,j)
          if( abs(dC) > 1.d-13 ) then
            write(700+j0,*) jday,i,j, dC
          endif
        enddo
      enddo
#endif

      return

      call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata)
      call prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)

      call ent_prescribe_vegupdateD(entcells
     &     ,laidata=laidata(2:2+N_PFT-1,:,:)
     &     ,albedodata=albedodata(:,2:2+N_PFT-1,:,:)
     &     )

      end subroutine update_vegetation_data

      subroutine map_ent2giss(v_ent,h_ent,v_giss)
      implicit none
      real*8, dimension(:), intent(in) :: v_ent
      real*8, dimension(:), intent(in) :: h_ent
      real*8, dimension(:), intent(out) :: v_giss

      if (N_COVERTYPES == 12) then
        v_giss(1:12) = v_ent(1:12)
      else if (N_COVERTYPES == 18) then
        !18->12
        v_giss(1) = v_ent(17)   !sand
        v_giss(2) = 0.d0 ! v_ent(9)    !tundra
        v_giss(3) = v_ent(11) + v_ent(12) 
     &          + v_ent(13) + v_ent(14) !grass
        v_giss(4) = v_ent(10) !!+ v_ent(9)  !shrub
        v_giss(5) = 0.d0        !tress
        v_giss(6) = v_ent(5) + v_ent(6) 
     &          + v_ent(7) + v_ent(8)
     &          + v_ent(16)             !deci
        v_giss(7) = v_ent(3) + v_ent(4) !evergr
        v_giss(8) = v_ent(1) + v_ent(2) !rainf
        v_giss(9) = v_ent(15)   !crops
        v_giss(10)= v_ent(18)   !dirt
        v_giss(11)= 0.d0
        v_giss(12)= 0.d0
        if ( h_ent(9) > .5d0 ) then
          v_giss(4) = v_giss(4) + v_ent(9)
        else
          v_giss(2) = v_giss(2) + v_ent(9)
        endif
      else
        call stop_model("map_ent2giss: unsupported N_COVERTYPES",255)
      endif


      end subroutine map_ent2giss


      subroutine read_laimax(I0,I1,J0,J1, laimax)
!@sum read maximum LAI from a file, 
!@+   return array with or without halo depending on I0,I1,J0,J1
      use domain_decomp_atm, only : grid
      use pario, only : par_open,par_close,read_dist_data
      integer, intent(in) :: I0,I1,J0,J1
      real*8, intent(out) :: laimax(:, I0:, J0:)
      !---
      real*8 :: buf( grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO )
      character*32 :: name
      character*8 :: i_char
      integer :: k, fid


cddd      call get(grid, J_STRT     =J_0,    J_STOP     =J_1,
cddd     &               I_STRT     =I_0,    I_STOP     =I_1)

      ! make sure it is 0 for non-veg cover types
      laimax(:,:,:) = 0.d0

      fid = par_open(grid,'LAIMAX','read')
      do k=1, N_PFT
        write(i_char,"(i8)") k
        !name = 'var_'//adjustl(i_char)
        name = ent_cover_names(k)
!        print *,"read_laimax: reading ", name
        buf = 0
        call read_dist_data(grid, fid, trim(name), buf)
!        print *, sum(buf)
!!! hack : make sure that laimax is always >= .5 to avoid zero-divisions
        where ( buf(I0:I1, J0:J1) > .1d0 )
          laimax(k, I0:I1, J0:J1) = buf(I0:I1, J0:J1)
        elsewhere
          laimax(k, I0:I1, J0:J1) = .1d0
        end where
      enddo
      call par_close(grid,fid)

      end subroutine read_laimax

      subroutine read_height(I0,I1,J0,J1, height)
!@sum read height of vegetation  from a file, 
!@+   return array with or without halo depending on I0,I1,J0,J1
      use domain_decomp_atm, only : grid
      use pario, only : par_open,par_close,read_dist_data
      integer, intent(in) :: I0,I1,J0,J1
      real*8, intent(out) :: height(:, I0:, J0:)
      !---
      real*8 :: buf( grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO )
      character*32 :: name
      character*8 :: i_char
      integer :: k, fid


cddd      call get(grid, J_STRT     =J_0,    J_STOP     =J_1,
cddd     &               I_STRT     =I_0,    I_STOP     =I_1)

      ! make sure it is 0 for non-veg cover types
      height(:,:,:) = 0.d0

      fid = par_open(grid,'HITEent','read')
      do k=1, N_PFT
        write(i_char,"(i8)") k
        !name = 'var_'//adjustl(i_char)
        name = "hgt_"//ent_cover_names(k)
!        print *,"read_height: reading ", name
        call read_dist_data(grid, fid, trim(name), buf)
!!! hack : make sure that height is always >= .5 to avoid zero-divisions
        where ( buf(I0:I1, J0:J1) > .5d0 )
          height(k, I0:I1, J0:J1) = buf(I0:I1, J0:J1)
        elsewhere
          height(k, I0:I1, J0:J1) = .5d0
        end where
      enddo
      call par_close(grid,fid)

      end subroutine read_height

      end module ent_drv
