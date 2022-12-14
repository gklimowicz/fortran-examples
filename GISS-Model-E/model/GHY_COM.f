#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif

      MODULE GHY_COM
!@sum  GHY_COM contains the areas used by the TerraE Global Land Model; 
!@+    code for TerraE is in 'model/giss_LSM' 
!@auth Frank Abramopolus/Igor Aleinov
      USE RESOLUTION, only : im,jm
!!!      USE SLE001, only : ngm,imt,nlsn
#ifdef TRACERS_WATER
      USE TRACER_COM, only : NTM
#endif


      IMPLICIT NONE
      SAVE

ccc CONSTANTS

ccc dimensions of the GHY arrays
!@var ngm number of soil layers
!@var imt number of soil textures
!@var nlsn max number of snow layers
      integer, parameter, public :: ngm=6, imt=5, nlsn=3
!@dbparam WSN_MAX snow amount limit (m, water equivalent)
      real*8, public :: WSN_MAX = 2.d0

!@var LS_NFRAC number of land surface fractions
      integer, parameter, public :: LS_NFRAC=3

!@var shc_soil_texture specific heat capacity of soil texture (J/K/M^3)
      real*8, parameter,public :: shc_soil_texture(imt)
     &     = (/2d6,2d6,2d6,2.5d6,2.4d6/)


ccc variable earth fraction
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FEARTH

!@var WFCS water field capacity of first ground layer (kg/m2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: WFCS

ccc bare/veg not in merged array because WBARE does not contain
ccc 0 index for legacy reasons
      !REAL*8, POINTER, DIMENSION(:,:,:) :: WBARE
      !REAL*8, POINTER, DIMENSION(:,:,:) :: WVEGE
      !REAL*8, POINTER, DIMENSION(:,:,:) :: HTBARE
      !REAL*8, POINTER, DIMENSION(:,:,:) :: HTVEGE
      REAL*8, ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: W_IJ
      REAL*8, ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: HT_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SNOWBV

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: DZ_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: QK_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:)     :: SL_IJ

ccc the following arrays contain prognostic variables for the snow model
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:)     :: NSN_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: DZSN_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: WSN_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: HSN_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)      :: FR_SNOW_IJ
C**** replacements for GDATA
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: TEARTH
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: WEARTH
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: AIEARTH

C**** arrays needed to restart from a SOILIC file with intensive units.
C**** Normally these are only needed during cold starts.
C**** If the SOILIC file is present, they are allocated by read_landice_ic,
C**** used in init_land_surface, and deallocated after use.
!@var earth_sat saturation of each soil layer (1 - completely saturated)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: earth_sat
!@var earth_ice fraction of frozen water in the layer
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: earth_ice
!@var earth_tp temperature of layer (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: earth_tp
!@var tsnowtop skin temperature of snow (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tsnowtop

!@var GDEEP keeps average (2:n) values of temperature, water and ice
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GDEEP
!@var GSAVEL indiv layers temp,water,ice (for diag exporting)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: GSAVEL

ccc topmodel input data and standard deviation of the elevation
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: TOP_INDEX_IJ, top_dev_ij

ccc evaporation limits from previous time step
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: evap_max_ij, fr_sat_ij,
     &     qg_ij

!@var tsns_ij surface temperature corresponding to sensible heat flux (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tsns_ij

!@var SNOWD snow depth (m) over bare and vegetated soil
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SNOWD

#ifdef TRACERS_WATER_OLD
!@var TRBARE,TRVEGE tracers in bare and veg. soil fraction (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRBARE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRVEGE
C**** What is the prognostic variable for snow here?
!@var TRSNOWBV tracer amount in snow over bare and veg. soil (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRSNOWBV
#endif
#ifdef TRACERS_WATER
ccc new tracers
      !integer, parameter :: NTM = 3
!!!@var TR_WBARE tracers in bare soil fraction (kg/m^2)
!!!@var TR_WVEGE tracers in vegetated soil fraction (kg/m^2)
!@var TR_W_IJ water tracers in soil (kg/m^2)
!@var TR_WSN_IJ tracer amount in snow (multiplied by fr_snow) (kg/m^2)
      !REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TR_WBARE
      !REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TR_WVEGE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TR_W_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TR_WSN_IJ
ccc TRSNOWBV is not used
!@var ntixw index array for tracers (shared by OMP threads)
!!      integer ntixw(ntm)
#endif
ccc stuff that got back from VEG_COM, maybe should be relocated to Ent
!@var Cint Internal foliage CO2 concentration (mol/m3)
      real*8, ALLOCATABLE, dimension(:,:) :: Ci_ij
!@var Qfol Foliage surface mixing ratio (kg/kg)
      real*8, ALLOCATABLE, dimension(:,:) :: Qf_ij
!@var cnc_ij canopy conductance
      real*8, ALLOCATABLE, dimension(:,:) :: cnc_ij

!@var aalbveg vegetation albedo, eventually should be moved to a
!@+   better place DO WE NEED THIS ???
      real*8, ALLOCATABLE, dimension(:,:) :: aalbveg

!@var soil_surf_moist near surf soil moisture (kg/m^3) for subdd
      real*8, ALLOCATABLE, dimension(:,:) :: soil_surf_moist
      END MODULE GHY_COM

      SUBROUTINE ALLOC_GHY_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE GHY_COM
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, getDomainBounds
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE(     FEARTH(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(WFCS(I_0H:I_1H,J_0H:J_1H), STAT = IER)

cddd      ALLOCATE(      WBARE(  NGM,I_0H:I_1H,J_0H:J_1H),
cddd     *               WVEGE(0:NGM,I_0H:I_1H,J_0H:J_1H),
cddd     *              HTBARE(0:NGM,I_0H:I_1H,J_0H:J_1H),
cddd     *              HTVEGE(0:NGM,I_0H:I_1H,J_0H:J_1H),
cddd     *              SNOWBV(    2,I_0H:I_1H,J_0H:J_1H),
cddd     *         STAT=IER)

      ALLOCATE(      W_IJ(0:NGM,LS_NFRAC,I_0H:I_1H,J_0H:J_1H),
     *              HT_IJ(0:NGM,LS_NFRAC,I_0H:I_1H,J_0H:J_1H),
     *              SNOWBV(     LS_NFRAC,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      !WBARE => W_IJ(1:,1,1:,1:)
      !WVEGE => W_IJ(1:,1,1:,1:)
      !HTBARE => HT_IJ(1:,1,1:,1:)
      !HTVEGE => HT_IJ(1:,1,1:,1:)

      ALLOCATE(     DZ_IJ(I_0H:I_1H,J_0H:J_1H,NGM),
     *               Q_IJ(I_0H:I_1H,J_0H:J_1H,IMT,NGM),
     *              QK_IJ(I_0H:I_1H,J_0H:J_1H,IMT,NGM),
     *              SL_IJ(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(      NSN_IJ(     2,I_0H:I_1H,J_0H:J_1H),
     *              DZSN_IJ(NLSN,2,I_0H:I_1H,J_0H:J_1H),
     *               WSN_IJ(NLSN,2,I_0H:I_1H,J_0H:J_1H),
     *               HSN_IJ(NLSN,2,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(         FR_SNOW_IJ(2,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(
     *              TEARTH(  I_0H:I_1H,J_0H:J_1H),
     *              WEARTH(  I_0H:I_1H,J_0H:J_1H),
     *             AIEARTH(  I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(     GDEEP(I_0H:I_1H,J_0H:J_1H,3),
     *         STAT=IER)
      ALLOCATE(    GSAVEL(I_0H:I_1H,J_0H:J_1H,ngm,3),
     *         STAT=IER)

      ALLOCATE(      TOP_INDEX_IJ(I_0H:I_1H,J_0H:J_1H),
     *                 top_dev_ij(I_0H:I_1H,J_0H:J_1H),
     *                evap_max_ij(I_0H:I_1H,J_0H:J_1H),
     *                  fr_sat_ij(I_0H:I_1H,J_0H:J_1H),
     *                      qg_ij(I_0H:I_1H,J_0H:J_1H),
     *           STAT=IER)
      ALLOCATE(     tsns_ij(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(      SNOWD(2,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(     aalbveg(I_0H:I_1H,J_0H:J_1H),
     *              Ci_ij(I_0H:I_1H,J_0H:J_1H),
     *              Qf_ij(I_0H:I_1H,J_0H:J_1H),
     *              cnc_ij(I_0H:I_1H,J_0H:J_1H),
     *           STAT=IER)

C**** Initialize evaporation limits
      evap_max_ij(:,J_0H:J_1H)=-1.
      fr_sat_ij(:,J_0H:J_1H)=1.
      qg_ij(:,J_0H:J_1H)=0.

ccc init snow arrays to prevent addressing uninitialized vars
      nsn_ij    (:,:,J_0H:J_1H)=0
      fr_snow_ij(:,:,J_0H:J_1H)=0.
      dzsn_ij (:,:,:,J_0H:J_1H)=0.
      hsn_ij  (:,:,:,J_0H:J_1H)=0.
      wsn_ij  (:,:,:,J_0H:J_1H)=0.

ccc make sure SNOWD is initialized
      SNOWD(:,:,:) = 0.d0

#ifdef TRACERS_WATER_OLD
      ALLOCATE(     TRBARE(NTM,  NGM,I_0H:I_1H,J_0H:J_1H),
     *              TRVEGE(NTM,0:NGM,I_0H:I_1H,J_0H:J_1H),
     *            TRSNOWBV(NTM,    2,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)
#endif
#ifdef TRACERS_WATER
cddd      ALLOCATE(     TR_WBARE(NTM,  NGM,I_0H:I_1H,J_0H:J_1H),
cddd     *              TR_WVEGE(NTM,0:NGM,I_0H:I_1H,J_0H:J_1H),
      ALLOCATE(      TR_W_IJ(NTM,0:NGM,LS_NFRAC,I_0H:I_1H,J_0H:J_1H),
     *             TR_WSN_IJ(NTM,NLSN,        2,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)
C**** Initialize to zero
      !TR_WBARE(:,:,:,J_0H:J_1H)=0.d0
      !TR_WVEGE(:,:,:,J_0H:J_1H)=0.d0
      TR_W_IJ(:,:,:,:,J_0H:J_1H)=0.d0
      TR_WSN_IJ(:,:,:,:,J_0H:J_1H)=0.d0
#endif

      ALLOCATE( soil_surf_moist(I_0H:I_1H,J_0H:J_1H) )
      soil_surf_moist(:,:) = 0.d0

      END SUBROUTINE ALLOC_GHY_COM

      subroutine read_landsurf_ic
!@sum   read_landsurf_ic read land surface initial conditions file.
      use model_com, only : ioread
      use domain_decomp_atm, only: grid,getdomainbounds
      use pario, only : par_open,par_close,read_dist_data
      use filemanager, only : file_exists
      use ghy_com
      implicit none
      integer :: fid
      integer :: ifrac
      real*8, dimension(:,:,:), allocatable :: temp,wetness
      real*8, dimension(:,:), allocatable :: snowdp
      integer :: i_0,i_1,j_0,j_1,i,j

      if(file_exists('SOILIC')) then
        fid = par_open(grid,'SOILIC','read')
        ! read IC having intensive units.  Bare and vegetated
        ! fractions start with the same temp/saturation/snow.
        call getdomainbounds(grid,
     &       i_strt=i_0,i_stop=i_1, j_strt=j_0,j_stop=j_1)
        allocate(temp(ngm,i_0:i_1,j_0:j_1),
     &        wetness(ngm,i_0:i_1,j_0:j_1),
     &         snowdp(i_0:i_1,j_0:j_1),
     &       tsnowtop(i_0:i_1,j_0:j_1),
     &        earth_tp(0:ngm,ls_nfrac,i_0:i_1,j_0:j_1),
     &       earth_sat(0:ngm,ls_nfrac,i_0:i_1,j_0:j_1),
     &       earth_ice(0:ngm,ls_nfrac,i_0:i_1,j_0:j_1))
        call read_dist_data(grid,fid,'temp',temp,jdim=3)
        wetness = .8d0 ! default in case wetness not present
        call read_dist_data(grid,fid,'wetness',wetness,jdim=3)
        call read_dist_data(grid,fid,'snow',snowdp)
        tsnowtop = temp(1,:,:) ! default in case tsnowtop not present
        call read_dist_data(grid,fid,'tsnowtop',tsnowtop)
        do ifrac=1,ls_nfrac
          do j=j_0,j_1
          do i=i_0,i_1
            earth_tp(1:ngm,ifrac,i,j) = temp(:,i,j)
            earth_tp(0,ifrac,i,j) = temp(1,i,j)
            earth_sat(1:ngm,ifrac,i,j) = wetness(:,i,j)
            earth_sat(0,ifrac,i,j) = 0.
            earth_ice(:,ifrac,i,j) = 0.
            where(temp(:,i,j).lt.0.) earth_ice(1:ngm,ifrac,i,j) = 1.
            snowbv(ifrac,i,j) = snowdp(i,j)
            if(temp(1,i,j).gt.0.) snowbv(ifrac,i,j) = 0.
          enddo
          enddo
        enddo
        call par_close(grid,fid)
      elseif(file_exists('SOIL')) then ! check SOIL presence, not GIC
        ! IC file with extensive units (per-layer heat and water amounts)
        fid = par_open(grid,'GIC','read')
        call new_io_earth  (fid,ioread)
        call new_io_soils  (fid,ioread)
        call par_close(grid,fid)
      endif
      return
      end subroutine read_landsurf_ic

      subroutine def_rsf_earth(fid)
!@sum  def_rsf_earth defines earth array structure in restart files
!@auth M. Kelley
!@ver  beta
      use ghy_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use fluxes, only : atmlnd
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,tearth,'tearth(dist_im,dist_jm)')
      call defvar(grid,fid,atmlnd%snowe,'snowe(dist_im,dist_jm)')
      call defvar(grid,fid,wearth,'wearth(dist_im,dist_jm)')
      call defvar(grid,fid,aiearth,'aiearth(dist_im,dist_jm)')
      call defvar(grid,fid,evap_max_ij,
     &     'evap_max_ij(dist_im,dist_jm)')
      call defvar(grid,fid,fr_sat_ij,
     &     'fr_sat_ij(dist_im,dist_jm)')
      call defvar(grid,fid,qg_ij,'qg_ij(dist_im,dist_jm)')
      call defvar(grid,fid,tsns_ij,'tsns_ij(dist_im,dist_jm)')
      call defvar(grid,fid,snowd,'snowd(bv,dist_im,dist_jm)')
      return
      end subroutine def_rsf_earth

      subroutine new_io_earth(fid,iaction)
!@sum  new_io_earth read/write earth arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use ghy_com
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use fluxes, only : atmlnd
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'tearth',tearth)
        call write_dist_data(grid,fid,'snowe',atmlnd%snowe)
        call write_dist_data(grid,fid,'wearth',wearth)
        call write_dist_data(grid,fid,'aiearth',aiearth)
        call write_dist_data(grid,fid,'evap_max_ij',evap_max_ij)
        call write_dist_data(grid,fid,'fr_sat_ij',fr_sat_ij)
        call write_dist_data(grid,fid,'qg_ij',qg_ij)
        call write_dist_data(grid,fid,'tsns_ij',tsns_ij)
        call write_dist_data(grid,fid,'snowd',snowd,jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'tearth',tearth)
        call read_dist_data(grid,fid,'snowe',atmlnd%snowe)
        call read_dist_data(grid,fid,'wearth',wearth)
        call read_dist_data(grid,fid,'aiearth',aiearth)
        call read_dist_data(grid,fid,'evap_max_ij',evap_max_ij)
        call read_dist_data(grid,fid,'fr_sat_ij',fr_sat_ij)
        call read_dist_data(grid,fid,'qg_ij',qg_ij)
        tsns_ij(:,:) = tearth(:,:) ! default if not in input file
        call read_dist_data(grid,fid,'tsns_ij',tsns_ij)
        call read_dist_data(grid,fid,'snowd',snowd,jdim=3)
      end select
      return
      end subroutine new_io_earth

      subroutine def_rsf_soils(fid)
!@sum  def_rsf_soils defines soils array structure in restart files
!@auth M. Kelley
!@ver  beta
      use ghy_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use conserv_diags
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,w_ij,
     &     'w_ij(zero_to_ngm,ls_nfrac,dist_im,dist_jm)')
      call defvar(grid,fid,ht_ij,
     &     'ht_ij(zero_to_ngm,ls_nfrac,dist_im,dist_jm)')
      call defvar(grid,fid,snowbv,
     &     'snowbv(ls_nfrac,dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid,fid,tr_w_ij,
     &     'tr_w_ij(ntm,zero_to_ngm,ls_nfrac,dist_im,dist_jm)')
#endif
      call declare_conserv_diags( grid, fid, 'wgrnd(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'egrnd(dist_im,dist_jm)' )
      return
      end subroutine def_rsf_soils

      subroutine new_io_soils(fid,iaction)
!@sum  new_io_soils read/write soil arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only: grid
      use pario, only : write_dist_data,read_dist_data
      use ghy_com
      use conserv_diags
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      external conserv_WTG, conserv_HTG
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'w_ij',w_ij,jdim=4)
        call write_dist_data(grid,fid,'ht_ij',ht_ij,jdim=4)
        call write_dist_data(grid,fid,'snowbv',snowbv,jdim=3)
#ifdef TRACERS_WATER
        call write_dist_data(grid,fid,'tr_w_ij',tr_w_ij,jdim=5)
#endif
        call dump_conserv_diags( grid, fid, 'wgrnd', conserv_WTG )
        call dump_conserv_diags( grid, fid, 'egrnd', conserv_HTG )
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'w_ij',w_ij,jdim=4)
        call read_dist_data(grid,fid,'ht_ij',ht_ij,jdim=4)
        call read_dist_data(grid,fid,'snowbv',snowbv,jdim=3)
#ifdef TRACERS_WATER
        call read_dist_data(grid,fid,'tr_w_ij',tr_w_ij,jdim=5)
#endif
      end select
      return
      end subroutine new_io_soils

      subroutine def_rsf_snow(fid)
!@sum  def_rsf_snow defines snow array structure in restart files
!@auth M. Kelley
!@ver  beta
      use ghy_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,dzsn_ij,
     &     'dzsn_ij(nlsn,bv,dist_im,dist_jm)')
      call defvar(grid,fid,wsn_ij,
     &     'wsn_ij(nlsn,bv,dist_im,dist_jm)')
      call defvar(grid,fid,hsn_ij,
     &     'hsn_ij(nlsn,bv,dist_im,dist_jm)')
      call defvar(grid,fid,nsn_ij,
     &     'nsn_ij(bv,dist_im,dist_jm)')
      call defvar(grid,fid,fr_snow_ij,
     &     'fr_snow_ij(bv,dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid,fid,tr_wsn_ij,
     &     'tr_wsn_ij(ntm,nlsn,bv,dist_im,dist_jm)')
#endif
      return
      end subroutine def_rsf_snow

      subroutine new_io_snow(fid,iaction)
!@sum  new_io_snow read/write snow arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use ghy_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to restart file
        call write_dist_data(grid, fid, 'dzsn_ij', dzsn_ij,jdim=4)
        call write_dist_data(grid, fid, 'wsn_ij', wsn_ij,jdim=4)
        call write_dist_data(grid, fid, 'hsn_ij', hsn_ij,jdim=4)
        call write_dist_data(grid, fid, 'nsn_ij', nsn_ij, jdim=3)
        call write_dist_data(grid, fid, 'fr_snow_ij',
     &       fr_snow_ij,jdim=3)
#ifdef TRACERS_WATER
        call write_dist_data(grid, fid, 'tr_wsn_ij',
     &       tr_wsn_ij,jdim=5)
#endif
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'dzsn_ij', dzsn_ij,jdim=4)
        call read_dist_data(grid, fid, 'wsn_ij', wsn_ij,jdim=4)
        call read_dist_data(grid, fid, 'hsn_ij', hsn_ij,jdim=4)
        call read_dist_data(grid, fid, 'nsn_ij', nsn_ij, jdim=3)
        call read_dist_data(grid, fid, 'fr_snow_ij',
     &       fr_snow_ij,jdim=3)
#ifdef TRACERS_WATER
        call read_dist_data(grid, fid, 'tr_wsn_ij',
     &       tr_wsn_ij,jdim=5)
#endif
      end select
      return      
      end subroutine new_io_snow

      subroutine def_rsf_veg_related(fid)
!@sum  def_rsf_snow defines veg_related array structure in restart files
!@auth M. Kelley
!@ver  beta
      use ghy_com, only : Ci_ij, Qf_ij, cnc_ij
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,ci_ij,'ci_ij(dist_im,dist_jm)')
      call defvar(grid,fid,qf_ij,'qf_ij(dist_im,dist_jm)')
      call defvar(grid,fid,cnc_ij,'cnc_ij_veg_related(dist_im,dist_jm)')
      return
      end subroutine def_rsf_veg_related

      subroutine new_io_veg_related(fid,iaction)
!@sum  new_io_snow read/write veg_related arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use ghy_com, only : Ci_ij, Qf_ij, cnc_ij
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to restart file
        call write_dist_data(grid, fid, 'ci_ij', ci_ij)
        call write_dist_data(grid, fid, 'qf_ij', qf_ij)
        call write_dist_data(grid, fid,
     &       'cnc_ij_veg_related', cnc_ij)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'ci_ij', ci_ij)
        call read_dist_data(grid, fid, 'qf_ij', qf_ij)
        call read_dist_data(grid, fid,
     &       'cnc_ij_veg_related', cnc_ij)
      end select
      return      
      end subroutine new_io_veg_related
