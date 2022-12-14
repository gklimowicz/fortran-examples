#include "rundeck_opts.h"

module cosp_drv
!@sum   Provides the ModelE interface with COSP
!@+     (CFMIP Observation Simulator Package).
!@auth  Mike Bauer
!@cont  init_cosp,init_cosp_gbx,run_cosp_sims,free_cosp_sims,
!@+     remap_cosp_2d,remap_cosp_3d
!@calls mod_cosp_constants,mod_cosp_modis_simulator,mod_cosp_types,
!@+     filemanager

  !----------------------------------------------------------------------------
  ! Purpose : Provides the ModelE interface with COSP.
  ! Platform: Linux
  ! Language: Fortran 90
  ! ---------------------------------------------------------------------------
  ! Description: Basic glue between the code-bases of COSP and ModelE.
  !
  ! Author: Mike Bauer
  !
  ! ToDo:
  !
  ! History:
  !   Jun 2011 - M. Bauer - Tweak to get COSP to run in-line in GISS modelE
  !   Nov 2011 - M. Bauer - Tweaked and added to ModelE/AR5 code base.
  !   Jul 2017 - M. Bauer - Tweaked and added to ModelE/AR5 code base for v1.4.1
  !   Dec 2017 - M. Bauer - Adapted to ModelE/CloudSandbox code base for v1.4.2
  !   May 2019 - M. Bauer - Adapted to ModelE/E21 code base for v1.4.3
  ! ---------------------------------------------------------------------------

  ! MODULE Imports
  ! ===========================================================================

  ! From COSP
  ! ---------
  ! Imported Parameters:
  ! use mod_cosp_constants, only : dbze_bins,i_cvcice,i_cvcliq,i_cvrain,i_cvsnow,&
  !   i_lscice,i_lscliq,i_lsgrpl,i_lsrain,i_lssnow,lidar_ntemp,lidar_ntemp,      &
  !   n_out_list,parasol_nrefl,rttov_max_channels,sr_bins
  use mod_cosp_constants, only : dbze_bins,n_out_list,parasol_nrefl,           &
    rttov_max_channels,sr_bins,misr_n_cth
  use mod_cosp_constants, only : i_lscliq,i_lscice,i_lsrain,i_lssnow,i_cvcliq, &
    i_cvcice,i_cvrain,i_cvsnow,i_lsgrpl
#ifdef COSP_DEBUG
    use mod_cosp_constants, only : isccp_tau,isccp_pc
#endif
  use mod_modis_sim, only : nummodisrefficebins,nummodisreffliqbins

  ! Imported Type Definitions:
  use mod_cosp_types
  use mod_cosp_modis_simulator, only : cosp_modis

  ! From ModelE
  ! -----------
  ! COSP_Note_1: As of Nov 2010 ModelE does not allow 'components' such as
  ! COSP to access main modelE modules such as MODEL_COM. Instead these
  ! values must be passed as arguments to subroutines etc. This might
  ! change at some point in which case USE statements are preferred.
  !
  ! Exceptions
  !   /Constants_mod.F90:constant
  !   /shared/FileManager_mod.F90:filemanager
  use filemanager, only: closeunit,openunit

  ! All variables, parameters, and functions must be declared
  ! ===========================================================================
  implicit none

  ! Limit Access, unless explicitly declared PUBLIC
  ! ===========================================================================
  private

  ! Module Parameters/Constants
  ! ===========================================================================
  ! CALIPSO array index for cloud layer phase types
  !   Index = Hydrometeor
  !   1     = ice
  !   2     = liquid
  !   3     = undefined
  !   4     = false ice
  !   5     = false liquid
  !   6     = Percent of ice
  integer, parameter :: nphase = 6
  ! Used for temperature phase arrays
  integer, parameter :: nphaset = 5
  ! Used for ISCCP histograms
  integer,parameter :: cosp_ntau = 7
  integer,parameter :: cosp_npres = 7
  ! Scaling factors
  real*8, parameter :: cosp_bywc = 1./2.56d0
  real*8, parameter :: cosp_byic = 1./2.13d0

  ! Missing value
  !  The COSP and ModelE missing value are the same
  !   In /shared/Constants_mod.F90
  !     ! Missing value
  !     real*8,parameter :: undef = -1.d30
  !     ! Small positive value used in num/(den+teeny) to avoid 0/0
  !     real*8,parameter :: teeny = 1.d-30
  !   In cosp_constants.F90
  !     real,parameter :: r_undef = -1.0E30
  real*8,parameter :: cosp_undef = -1.d30
  real*8,parameter :: cosp_teeny = 1.d-30
  !
  ! WARNING of a major gotcha!
  !
  ! Be aware that if COSP is told to use fixed pressure levels (USE_VGRID = True)
  ! then COSP output fields will contain cosp_undef where the topography
  ! intersects a pressure level. However, if USE_VGRID = False, COSP output
  ! fields are on model (sigma) levels and cosp_undef only occur for
  ! instances like night-time for ISCCP.
  !
  ! Not taking this into account clearly will corrupt your statistics,
  ! weather state checks and the like.

  ! Module Scalers/Variables
  ! ===========================================================================
  integer :: channel_axid,column_axid,dbze_axid,grid_id,height_axid,           &
    height_mlev_axid,lat_axid,latvar_id,lon_axid,lonvar_id,misr_cth_axid,      &
    pressure2_axid,reffice_axid,reffliq_axid,sratio_axid,sza_axid,tau_axid,    &
    temp_axid,time_axid
  ! Special trigger for when want to collect ISCCP Cloud Area Fraction CFAD
  ! every cloud-time step
  logical :: isccp_always = .false.
  logical :: isccp_ave_now = .false.

  ! COSP_DEBUG Triggers very local, debugging etc., directed to UNIT=0 == stderr
  !  To redirect this to file as follows: $ command 2> errors.txt
  integer :: cosp_out_unit

#ifdef COSP_DEBUG
  ! For deriving statistics of the model values input into COSP
  real*8 :: min_dtau_c,min_dem_c_cvcliq,min_dem_c_cvcice,min_reff_cvcliq,      &
    min_mrhy_cvcliq,min_reff_cvcice,min_mrhy_cvcice,min_reff_cvrain,           &
    min_mrhy_cvrain,min_reff_cvsnow,min_mrhy_cvsnow
  real*8 :: max_dtau_c,max_dem_c_cvcliq,max_dem_c_cvcice,max_reff_cvcliq,      &
    max_mrhy_cvcliq,max_reff_cvcice,max_mrhy_cvcice,max_reff_cvrain,           &
    max_mrhy_cvrain,max_reff_cvsnow,max_mrhy_cvsnow
  real*8 :: min_dtau_s,min_dem_s_lscliq,min_dem_s_lscice,min_reff_lscliq,      &
    min_mrhy_lscliq,min_reff_lscice,min_mrhy_lscice,min_reff_lsrain,           &
    min_mrhy_lsrain,min_reff_lssnow,min_mrhy_lssnow
  real*8 :: max_dtau_s,max_dem_s_lscliq,max_dem_s_lscice,max_reff_lscliq,      &
    max_mrhy_lscliq,max_reff_lscice,max_mrhy_lscice, max_reff_lsrain,          &
    max_mrhy_lsrain,max_reff_lssnow,max_mrhy_lssnow
#endif

  ! Module Derived Type Definitions
  ! ===========================================================================
  type cosp_nmlst
    ! Holds configuration options (input, normally from cosp_input_nl.txt)
    !   See cosp_input_nl.txt and COSP_INPUT and in init_cosp().
    !
    ! Similar to type(cosp_config) :: cfg provided by COSP for working
    !   with cosp_output_nl.txt and COSP_OUTPUT in init_cosp().

    ! Type Scalers/Variables
    integer :: do_ray,instrument,isccp_topheight,isccp_topheight_direction,    &
      lidar_ice_type,melt_lay,naero,nchannels,ncolumns,nlevels,nlr,npoints,    &
      npoints_it,nprmts_max_aero,nprmts_max_hydro,overlap,platform,            &
      satellite,surface_radar,use_gas_abs,use_mie_tables
    logical :: csat_vgrid,use_precipitation_fluxes,use_reff,use_vgrid
    real :: ch4,co,co2,k2,n2o,radar_freq,zenang
    character(len=512) :: finput,cmor_nl
    character(len=600) :: dinput

    ! Type Arrays
    integer,dimension(rttov_max_channels) :: channels
    real,dimension(rttov_max_channels) :: surfem
  end type cosp_nmlst

  type subdd_info
    ! Holds various scalers/variables for the model-COSP interface.
    logical           :: cflag     ! COSP output flag
    character(len=16) :: sname     ! COSP variable short-name
    character(len=52) :: lname     ! COSP variable long-name/title
    character(len=6)  :: units     ! COSP variable units
    integer           :: dim1      ! Length of extra-dimension 1
    integer           :: dim2      ! Length of extra-dimension 2
    logical           :: is_3d     ! Call remap_cosp_3d
  end type subdd_info

  ! Module Derived Type Declarations
  ! ===========================================================================
  ! Define COSP data structures
  type(cosp_config)     :: cfg
  type(cosp_gridbox)    :: gbx
  type(cosp_isccp)      :: isccp
  type(cosp_lidarstats) :: stlidar
  type(cosp_misr)       :: misr
  type(cosp_modis)      :: modis
  type(cosp_nmlst)      :: cosp_inputs
  type(cosp_radarstats) :: stradar
  type(cosp_sglidar)    :: sglidar
  type(cosp_sgradar)    :: sgradar
  type(cosp_subgrid)    :: sgx
  type(cosp_vgrid)      :: vgrid
  ! RTTOV not used for CFMIP
#ifdef RTTOV
  type(cosp_rttov)      :: rttov
#endif
  type(subdd_info),dimension(n_out_list) :: cvar

  ! ^^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^^
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  ! Manage Access, by explicitly declaring PUBLIC
  ! ===========================================================================

  ! Public Module Parameters:
  public :: cosp_byic,cosp_bywc,dbze_bins,lidar_ntemp,nphase,nphaset,          &
    parasol_nrefl,sr_bins,cosp_ntau,cosp_npres,n_out_list,misr_n_cth,          &
    nummodisrefficebins,nummodisreffliqbins,isccp_always,isccp_ave_now
  public :: i_lscliq,i_lscice,i_lsrain,i_lssnow,i_cvcliq,i_cvcice,i_cvrain,    &
    i_cvsnow,i_lsgrpl
  !public : cosp_snames, cosp_units, cosp_lnames

  ! Public Module Scalar Variables:
#ifdef COSP_DEBUG
  public :: min_dtau_c,min_dem_c_cvcliq,min_dem_c_cvcice,min_reff_cvcliq,      &
    min_mrhy_cvcliq,min_reff_cvcice,min_mrhy_cvcice,min_reff_cvrain,           &
    min_mrhy_cvrain,min_reff_cvsnow,min_mrhy_cvsnow,max_dtau_c,                &
    max_dem_c_cvcliq,max_dem_c_cvcice,max_reff_cvcliq,max_mrhy_cvcliq,         &
    max_reff_cvcice,max_mrhy_cvcice,max_reff_cvrain,max_mrhy_cvrain,           &
    max_reff_cvsnow,max_mrhy_cvsnow,min_dtau_s,min_dem_s_lscliq,               &
    min_dem_s_lscice,min_reff_lscliq,min_mrhy_lscliq,min_reff_lscice,          &
    min_mrhy_lscice,min_reff_lsrain,min_mrhy_lsrain,min_reff_lssnow,           &
    min_mrhy_lssnow,max_dtau_s,max_dem_s_lscliq,max_dem_s_lscice,              &
    max_reff_lscliq,max_mrhy_lscliq,max_reff_lscice,max_mrhy_lscice,           &
    max_reff_lsrain,max_mrhy_lsrain,max_reff_lssnow,max_mrhy_lssnow
#endif

  ! Public Module Type Definitions:
  public :: cosp_nmlst

  ! Public Module Type Declarations:
  public :: cfg,gbx,isccp,stlidar,misr,modis,stradar,cosp_inputs,cvar

  ! Public Module Routines:
  public :: init_cosp,init_cosp_gbx,run_cosp_sims,free_cosp_sims,              &
    remap_cosp_2d,remap_cosp_3d
#ifdef COSP_DEBUG
  public :: cosp_print_pctau,print_cosp_gbx,cosp_input_check
#endif

  ! Internal Procedures
  ! ===========================================================================
  contains

  ! ***************************************************************************
  ! PUBLIC SUBROUTINE: init_cosp
  ! ****************************
  subroutine init_cosp(model_npoints,model_im,model_lm,flag_pfluxes,flag_re)
!@sum   Read COSP namelist options and configure COSP input/output options.
!@auth  Mike Bauer
!@usage Only needs to be called once early in run. Depends on namelist files.
!@calls openunit,closeunit
    ! -------------------------------------------------------------------------
    ! Inputs:
    !   model_npoints: Number of model grids/Columns being passed to COSP.
    !   model_im     : Number of model longitude grids.
    !   model_lm     : Number of model layers
    !   flag_pfluxes : Reflects the rundeck pre-processor option COSP_PFLUX.
    !                  This must be consistent with use_precipitation_fluxes
    !                  from COSP_INPUT.
    !   flag_re      : Reflects the rundeck pre-processor option COSP_USERE.
    !                  This must be consistent with use_reff from COSP_INPUT.
    ! Outputs:
    !   Populates the COSP data structures cfg and cosp_inputs.
    !
    ! Reads COSP namelists (COSPIN and COSPOUT in RunDeck) which
    !   * Sets model specific parameters needed by COSP.
    !   * Sets which COSP simulators and output variables are being requested.
    !
    ! ToDo:
    !
    ! Language: Fortran 90
    ! -------------------------------------------------------------------------

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! Subroutine Arguments, in order of appearance
    ! =========================================================================
    integer,intent(in) :: model_npoints,model_im,model_lm
    logical,intent(in) :: flag_pfluxes,flag_re

    ! Subroutine Scalers/Variables
    ! =========================================================================
    integer :: iu_cospin,iu_cospout
    character(len=512) :: finput,cmor_nl
    character(len=600) :: dinput
    character(len=16) :: search_str
    integer :: do_ray,i,instrument,nchannels,isccp_topheight,                  &
        isccp_topheight_direction,lidar_ice_type,melt_lay,naero,ncolumns,      &
        nlevels,nlr,npoints,npoints_it,nprmts_max_aero,nprmts_max_hydro,       &
        overlap,platform,satellite,surface_radar,use_gas_abs,use_mie_tables,   &
        nlev
    logical :: csat_vgrid,use_precipitation_fluxes,use_reff,use_vgrid
    logical :: lalbisccp,latb532,lboxptopisccp,lboxtauisccp,lcfaddbze94,       &
      lcfadlidarsr532,lclcalipso,lclcalipso2,lclcalipsoice,lclcalipsoliq,      &
      lclcalipsotmp,lclcalipsotmpice,lclcalipsotmpliq,lclcalipsotmpun,         &
      lclcalipsoun,lclhcalipso,lclhcalipsoice,lclhcalipsoliq,lclhcalipsoun,    &
      lclhmodis,lclimodis,lclisccp,lcllcalipso,lcllcalipsoice,lcllcalipsoliq,  &
      lcllcalipsoun,lcllmodis,lclmcalipso,lclmcalipsoice,lclmcalipsoliq,       &
      lclmcalipsoun,lclmisr,lclmmodis,lclmodis,lcltcalipso,lcltcalipsoice,     &
      lcltcalipsoliq,lcltcalipsoun,lcltisccp,lcltlidarradar,lcltmodis,         &
      lclwmodis,lcrimodis,lcrlmodis,ldbze94,lfracout,lisccp_sim,liwpmodis,     &
      llidar_sim,llidarbetamol532,llwpmodis,lmeantbclrisccp,lmeantbisccp,      &
      lmisr_sim,lmodis_sim,lparasolrefl,lpctisccp,lpctmodis,lradar_sim,        &
      lreffclimodis,lreffclwmodis,lrttov_sim,ltauilogmodis,ltauimodis,         &
      ltauisccp,ltautlogmodis,ltautmodis,ltauwlogmodis,ltauwmodis,ltbrttov,    &
      ltoffset
    real :: ch4,co,co2,k2,n2o,radar_freq,zenang
    character(len=*), parameter :: twod_fmt = '("i = ",I2,/,5x,"cfg%",A,       &
      " = ",L2,/,5x,"sname = ",A,/,5x,"lname = ",A,/,5x,"units = ",A)'
    character(len=*), parameter :: threed_fmt = '("i = ",I2,/,5x,"cfg%",A,     &
      " = ",L2,/,5x,"sname = ",A,/,5x,"lname = ",A,/,5x,"units = ",A,/,5x,     &
      "DIM1  = ",I4,/,5x,"DIM2  = ",I4)'

    ! Subroutine Arrays
    ! =========================================================================
    integer,dimension(rttov_max_channels) :: channels
    logical,dimension(n_out_list) :: out_list_flags
    real,dimension(rttov_max_channels) :: surfem
    character(len=16),dimension(n_out_list) :: out_list_labels
    character(len=16),dimension(n_out_list) :: cosp_snames = (/                &
      'Lalbisccp       ', 'Latb532         ', 'Lboxptopisccp   ',              &
      'Lboxtauisccp    ', 'Lcfaddbze94     ', 'Lcfadlidarsr532 ',              &
      'Lclcalipso2     ', 'Lclcalipso      ', 'Lclhcalipso     ',              &
      'Lclisccp        ', 'Lcllcalipso     ', 'Lclmcalipso     ',              &
      'Lcltcalipso     ', 'Lcllcalipsoice  ', 'Lclmcalipsoice  ',              &
      'Lclhcalipsoice  ', 'Lcltcalipsoice  ', 'Lcllcalipsoliq  ',              &
      'Lclmcalipsoliq  ', 'Lclhcalipsoliq  ', 'Lcltcalipsoliq  ',              &
      'Lcllcalipsoun   ', 'Lclmcalipsoun   ', 'Lclhcalipsoun   ',              &
      'Lcltcalipsoun   ', 'Lclcalipsoice   ', 'Lclcalipsoliq   ',              &
      'Lclcalipsoun    ', 'Lclcalipsotmp   ', 'Lclcalipsotmpice',              &
      'Lclcalipsotmpliq', 'Lclcalipsotmpun ', 'Lcltlidarradar  ',              &
      'Lpctisccp       ', 'Ldbze94         ', 'Ltauisccp       ',              &
      'Lcltisccp       ', 'Ltoffset        ', 'Lparasolrefl    ',              &
      'Lclmisr         ', 'Lmeantbisccp    ', 'Lmeantbclrisccp ',              &
      'Lfracout        ', 'Llidarbetamol532', 'Ltbrttov        ',              &
      'Lcltmodis       ', 'Lclwmodis       ', 'Lclimodis       ',              &
      'Lclhmodis       ', 'Lclmmodis       ', 'Lcllmodis       ',              &
      'Ltautmodis      ', 'Ltauwmodis      ', 'Ltauimodis      ',              &
      'Ltautlogmodis   ', 'Ltauwlogmodis   ', 'Ltauilogmodis   ',              &
      'Lreffclwmodis   ', 'Lreffclimodis   ', 'Lpctmodis       ',              &
      'Llwpmodis       ', 'Liwpmodis       ', 'Lclmodis        ',              &
      'Lcrimodis       ', 'Lcrlmodis       ' /)
    character(len=16),dimension(n_out_list) :: cosp_units = (/                 &
      '1     ','      ','      ','      ','1     ','1     ','%     ','%     ', &
      '%     ','%     ','%     ','%     ','%     ','%     ','%     ','%     ', &
      '%     ','%     ','%     ','%     ','%     ','%     ','%     ','%     ', &
      '%     ','%     ','%     ','%     ','C     ','%     ','%     ','%     ', &
      '%     ','Pa    ','      ','1     ','%     ','      ','1     ','%     ', &
      'K     ','K     ','      ','      ','      ','%     ','%     ','%     ', &
      '%     ','%     ','%     ','1     ','1     ','1     ','1     ','1     ', &
      '1     ','m     ','m     ','Pa    ','kg m-2','kg m-2','1     ','%     ', &
      '%     ' /)
    character(len=50),dimension(n_out_list) :: cosp_lnames = (/                &
      'ISCCP Mean Cloud Albedo                             ',                  &
      '                                                    ',                  &
      '                                                    ',                  &
      '                                                    ',                  &
      'CloudSat Radar Reflectivity CFAD                    ',                  &
      'CALIPSO Scattering Ratio CFAD                       ',                  &
      'CALIPSO Cloud Fraction Undetected by CloudSat       ',                  &
      'CALIPSO Cloud Fraction                              ',                  &
      'CALIPSO High Level Cloud Fraction                   ',                  &
      'ISCCP Cloud Area Fraction CFAD                      ',                  &
      'CALIPSO Low Level Cloud Fraction                    ',                  &
      'CALIPSO Mid Level Cloud Fraction                    ',                  &
      'CALIPSO Total Cloud Area Fraction                   ',                  &
      'CALIPSO Low Level Ice Cloud Fraction                ',                  &
      'CALIPSO Mid Level Ice Cloud Fraction                ',                  &
      'CALIPSO High Level Ice Cloud Fraction               ',                  &
      'CALIPSO Ice Total Cloud Fraction                    ',                  &
      'CALIPSO Low Level Liquid Cloud Fraction             ',                  &
      'CALIPSO Mid Level Liquid Cloud Fraction             ',                  &
      'CALIPSO High Level Liquid Cloud Fraction            ',                  &
      'CALIPSO Liquid Total Cloud Fraction                 ',                  &
      'CALIPSO Low Undefined-Phase Cloud Fraction          ',                  &
      'CALIPSO Mid Undefined-Phase Cloud Fraction          ',                  &
      'CALIPSO High Undefined-Phase Cloud Fraction         ',                  &
      'CALIPSO Undefined-Phase Total Cloud Fraction        ',                  &
      'CALIPSO 3D Ice Cloud Fraction                       ',                  &
      'CALIPSO 3D Liquid Cloud Fraction                    ',                  &
      'CALIPSO 3D Undefined-Phase Cloud Fraction           ',                  &
      'CALIPSO 3D Cloud Temperature                        ',                  &
      'CALIPSO Ice Cloud Fraction Undetected by CloudSat   ',                  &
      'CALIPSO Liquid Cloud Fraction Undetected by CloudSat',                  &
      'CALIPSO 3D Undefined-Phase Cloud Temperature        ',                  &
      'Lidar and Radar Total Cloud Fraction                ',                  &
      'ISCCP Mean Cloud Top Pressure                       ',                  &
      '                                                    ',                  &
      'ISCCP Mean Optical Depth                            ',                  &
      'ISCCP Total Cloud Fraction                          ',                  &
      '                                                    ',                  &
      'PARASOL Reflectance                                 ',                  &
      'MISR Cloud Fraction CFAD                            ',                  &
      'ISCCP Mean all-sky 10.5um brightness temp           ',                  &
      'ISCCP Mean clear-sky 10.5um brightness temp         ',                  &
      '                                                    ',                  &
      '                                                    ',                  &
      '                                                    ',                  &
      'MODIS Total Cloud Fraction                          ',                  &
      'MODIS Liquid Cloud Fraction                         ',                  &
      'MODIS Ice Cloud Fraction                            ',                  &
      'MODIS High Level Cloud Fraction                     ',                  &
      'MODIS Mid Level Cloud Fraction                      ',                  &
      'MODIS Low Level Cloud Fraction                      ',                  &
      'MODIS Total Cloud Optical Thickness                 ',                  &
      'MODIS Liquid Cloud Optical Thickness                ',                  &
      'MODIS Ice Cloud Optical Thickness                   ',                  &
      'MODIS Total Cloud Optical Thickness (Log10 Mean)    ',                  &
      'MODIS Liquid Cloud Optical Thickness (Log10 Mean)   ',                  &
      'MODIS Ice Cloud Optical Thickness (Log10 Mean)      ',                  &
      'MODIS Liquid Cloud Particle Size                    ',                  &
      'MODIS Ice Cloud Particle Size                       ',                  &
      'MODIS Cloud Top Pressure                            ',                  &
      'MODIS Cloud Liquid Water Path                       ',                  &
      'MODIS Cloud Ice Water Path                          ',                  &
      'MODIS PC-Tau Histogram                              ',                  &
      'MODIS Optical_Thickness_vs_ReffIce Histogram        ',                  &
      'MODIS Optical_Thickness_vs_ReffLiq Histogram        ' /)

    integer,dimension(n_out_list) :: dim1 = (/ 0,0,0,0,dbze_bins,sr_bins,      &
      1,1,0,cosp_ntau,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,               &
      1,0,0,0,0,0,0,1,cosp_ntau,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,   &
      cosp_ntau,cosp_ntau,cosp_ntau/)
    integer,dimension(n_out_list) :: dim2

    NAMELIST/COSP_INPUT/cmor_nl,npoints,npoints_it,ncolumns,nlevels,use_vgrid, &
      nlr,csat_vgrid,dinput,finput,radar_freq,surface_radar,use_mie_tables,    &
      use_gas_abs,do_ray,melt_lay,k2,use_reff,use_precipitation_fluxes,        &
      nprmts_max_hydro,naero,nprmts_max_aero,lidar_ice_type,overlap,           &
      isccp_topheight,isccp_topheight_direction,platform,satellite,instrument, &
      nchannels,channels,surfem,zenang,co2,ch4,n2o,co

    NAMELIST/COSP_OUTPUT/lradar_sim,llidar_sim,lisccp_sim,lmodis_sim,lmisr_sim,&
      lrttov_sim,lalbisccp,latb532,lboxptopisccp,lboxtauisccp,lcfaddbze94,     &
      lcfadlidarsr532,lclcalipso2,lclcalipso,lclhcalipso,lclisccp,lcllcalipso, &
      lclmcalipso,lcltcalipso,lcltlidarradar,lpctisccp,ldbze94,ltauisccp,      &
      lclcalipsoliq,lclcalipsoice,lclcalipsoun,lclcalipsotmp,lclcalipsotmpliq, &
      lclcalipsotmpice,lclcalipsotmpun,lcltcalipsoliq,lcltcalipsoice,          &
      lcltcalipsoun,lclhcalipsoliq,lclhcalipsoice,lclhcalipsoun,lclmcalipsoliq,&
      lclmcalipsoice,lclmcalipsoun,lcllcalipsoliq,lcllcalipsoice,lcllcalipsoun,&
      lcltisccp,ltoffset,lparasolrefl,lclmisr,lmeantbclrisccp,lmeantbisccp,    &
      lfracout,llidarbetamol532,ltbrttov,lcltmodis,lclwmodis,lclimodis,        &
      lclhmodis,lclmmodis,lcllmodis,ltautmodis,ltauwmodis,ltauimodis,          &
      ltautlogmodis,ltauwlogmodis,ltauilogmodis,lreffclwmodis,lreffclimodis,   &
      lpctmodis,llwpmodis,liwpmodis,lclmodis,lcrimodis,lcrlmodis
    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! Manage Access, by explicitly declaring PUBLIC
    ! ===========================================================================

#ifdef COSP_DEBUG
    write(0,fmt='(25x,"model_npoints      = ",I6)') model_npoints
    write(0,fmt='(25x,"model_im           = ",i5)') model_im
    write(0,fmt='(25x,"model_lm           = ",i5)') model_lm
    write(0,fmt='(25x,"flag_pfluxes       = ",l5)') flag_pfluxes
    write(0,fmt='(25x,"flag_re            = ",l5)') flag_re
#endif

    ! Read COSP namelist options
    ! -------------------------------------------------------------------------
    call openunit('COSPIN',iu_cospin,.false.,.true.)
    read(iu_cospin,nml=cosp_input)
    call closeunit(iu_cospin)

    call openunit('COSPOUT',iu_cospout,.false.,.true.)
    read(iu_cospout,nml=cosp_output)
    call closeunit(iu_cospout)

    ! Initialize COSP data-structure (configure & allocate memory)
    ! -------------------------------------------------------------------------
    cosp_inputs%finput                           = finput
    cosp_inputs%dinput                           = dinput
    cosp_inputs%cmor_nl                          = cmor_nl
    cosp_inputs%overlap                          = overlap
    ! Because ModelE has a polar cap... and because of domain decomposition
    !  npoints might not be the full model grid resolution.
    cosp_inputs%npoints                          = model_npoints !npoints
    ! Tell COSP to iterate over whole latitude circles (more memory efficient)
    cosp_inputs%npoints_it                       = model_im
    !cosp_inputs%npoints_it                       = model_npoints
    cosp_inputs%ncolumns                         = ncolumns
    cosp_inputs%nlevels                          = model_lm !nlevels
    cosp_inputs%nlr                              = nlr
    cosp_inputs%surface_radar                    = surface_radar
    cosp_inputs%use_mie_tables                   = use_mie_tables
    cosp_inputs%use_gas_abs                      = use_gas_abs
    cosp_inputs%do_ray                           = do_ray
    cosp_inputs%melt_lay                         = melt_lay
    cosp_inputs%k2                               = k2
    cosp_inputs%nprmts_max_hydro                 = nprmts_max_hydro
    cosp_inputs%naero                            = naero
    cosp_inputs%nprmts_max_aero                  = nprmts_max_aero
    cosp_inputs%lidar_ice_type                   = lidar_ice_type
    cosp_inputs%isccp_topheight                  = isccp_topheight
    cosp_inputs%isccp_topheight_direction        = isccp_topheight_direction
    cosp_inputs%platform                         = platform
    cosp_inputs%satellite                        = satellite
    cosp_inputs%instrument                       = instrument
    cosp_inputs%nchannels                        = nchannels
    cosp_inputs%use_vgrid                        = use_vgrid
    cosp_inputs%csat_vgrid                       = csat_vgrid
    cosp_inputs%use_reff                         = use_reff
    cosp_inputs%use_precipitation_fluxes         = use_precipitation_fluxes
    cosp_inputs%radar_freq                       = radar_freq
    cosp_inputs%zenang                           = zenang
    cosp_inputs%co2                              = co2
    cosp_inputs%ch4                              = ch4
    cosp_inputs%n2o                              = n2o
    cosp_inputs%co                               = co
    cosp_inputs%channels(:cosp_inputs%nchannels) = channels
    cosp_inputs%surfem(:cosp_inputs%nchannels)   = surfem

    ! Check consistency of rundeck options and inputs
    ! -------------------------------------------------------------------------
    if (model_lm.ne.nlevels) then
      write(0,*) "RunDeck nlevels conflicts model_lm"
      call stop_model("nlevels conflicts model_lm",255)
    endif
    if (flag_pfluxes.ne.cosp_inputs%use_precipitation_fluxes) then
      write(0,*) "RunDeck COSP_PFLUX conflicts use_precipitation_fluxes"
      call stop_model("COSP_PFLUX conflicts use_precipitation_fluxes",255)
    endif
    if (flag_re.ne.cosp_inputs%use_reff) then
      write(0,*) "RunDeck COSP_USERE conflicts Input use_reff"
      call stop_model("RunDeck COSP_USERE conflicts Input use_reff",255)
    endif

#ifdef COSP_DEBUG
    call cosp_input_print(cosp_inputs)
#endif

    ! Deal with dependencies in requested simulators and outputs
    ! -------------------------------------------------------------------------
    if (.not.lradar_sim) then
      lcfaddbze94      = .false.; lclcalipso2      = .false.
      lclcalipso2      = .false.; lcltlidarradar   = .false.
      ldbze94          = .false.
    endif
    if (.not.llidar_sim) then
      latb532          = .false.; lcfadlidarsr532  = .false.
      lclcalipso       = .false.; lclcalipso2      = .false.
      lclhcalipso      = .false.; lcllcalipso      = .false.
      lclmcalipso      = .false.; lcltcalipso      = .false.
      lcltlidarradar   = .false.; lcltlidarradar   = .false.
      llidarbetamol532 = .false.; lparasolrefl     = .false.
      lclcalipsoice    = .false.; lclcalipsoliq    = .false.
      lclcalipsotmp    = .false.; lclcalipsotmpice = .false.
      lclcalipsotmpliq = .false.; lclcalipsotmpun  = .false.
      lclcalipsoun     = .false.; lclhcalipsoice   = .false.
      lclhcalipsoliq   = .false.; lclhcalipsoun    = .false.
      lcllcalipsoice   = .false.; lcllcalipsoliq   = .false.
      lcllcalipsoun    = .false.; lclmcalipsoice   = .false.
      lclmcalipsoliq   = .false.; lclmcalipsoun    = .false.
      lcltcalipsoice   = .false.; lcltcalipsoliq   = .false.
      lcltcalipsoun    = .false.
    endif
    if (.not.lisccp_sim) then
      lalbisccp        = .false.; lboxptopisccp    = .false.
      lboxtauisccp     = .false.; lclisccp         = .false.
      lcltisccp        = .false.; lmeantbclrisccp  = .false.
      lmeantbisccp     = .false.; lpctisccp        = .false.
      ltauisccp        = .false.
    endif
    if (.not.lmisr_sim) then
      lclmisr          = .false.
    endif
    if (.not.lrttov_sim) then
      ltbrttov         = .false.
    endif
    if ((.not.lradar_sim).and.(.not.llidar_sim).and.                           &
      (.not.lisccp_sim).and.(.not.lmisr_sim)) then
        lfracout         = .false.
    endif
    if (.not.lmodis_sim) then
      lcltmodis        = .false.; lclwmodis        = .false.
      lclimodis        = .false.; lclhmodis        = .false.
      lclmmodis        = .false.; lcllmodis        = .false.
      ltautmodis       = .false.; ltauwmodis       = .false.
      ltauimodis       = .false.; ltautlogmodis    = .false.
      ltauwlogmodis    = .false.; ltauilogmodis    = .false.
      lreffclwmodis    = .false.; lreffclimodis    = .false.
      lpctmodis        = .false.; llwpmodis        = .false.
      liwpmodis        = .false.; lclmodis         = .false.
      lcrimodis        = .false.; lcrlmodis        = .false.
    endif
    if (lmodis_sim) lisccp_sim = .true.
    cfg%lstats           = .false.
    if ((lradar_sim).or.(llidar_sim).or.(lisccp_sim)) cfg%lstats = .true.

    ! Warning the order of this matters, so ensure correct
    out_list_flags = (/ lalbisccp,latb532,lboxptopisccp,lboxtauisccp,          &
      lcfaddbze94,lcfadlidarsr532,lclcalipso2,lclcalipso,lclhcalipso,lclisccp, &
      lcllcalipso,lclmcalipso,lcltcalipso,lcllcalipsoice,lclmcalipsoice,       &
      lclhcalipsoice,lcltcalipsoice,lcllcalipsoliq,lclmcalipsoliq,             &
      lclhcalipsoliq,lcltcalipsoliq,lcllcalipsoun,lclmcalipsoun,lclhcalipsoun, &
      lcltcalipsoun,lclcalipsoice,lclcalipsoliq,lclcalipsoun,lclcalipsotmp,    &
      lclcalipsotmpice,lclcalipsotmpliq,lclcalipsotmpun,lcltlidarradar,        &
      lpctisccp,ldbze94,ltauisccp,lcltisccp,ltoffset,lparasolrefl,lclmisr,     &
      lmeantbisccp,lmeantbclrisccp,lfracout,llidarbetamol532,ltbrttov,         &
      lcltmodis,lclwmodis,lclimodis,lclhmodis,lclmmodis,lcllmodis,ltautmodis,  &
      ltauwmodis,ltauimodis,ltautlogmodis,ltauwlogmodis,ltauilogmodis,         &
      lreffclwmodis,lreffclimodis,lpctmodis,llwpmodis,liwpmodis,lclmodis,      &
      lcrimodis,lcrlmodis /)

    ! Warning the order of this matters, so ensure correct
    out_list_labels = (/ 'albisccp','atb532','boxptopisccp','boxtauisccp',     &
      'cfaddbze94','cfadlidarsr532','clcalipso2','clcalipso','clhcalipso',     &
      'clisccp','cllcalipso','clmcalipso','cltcalipso','cllcalipsoice',        &
      'clmcalipsoice','clhcalipsoice','cltcalipsoice','cllcalipsoliq',         &
      'clmcalipsoliq','clhcalipsoliq','cltcalipsoliq','cllcalipsoun',          &
      'clmcalipsoun','clhcalipsoun','cltcalipsoun','clcalipsoice',             &
      'clcalipsoliq','clcalipsoun','clcalipsotmp','clcalipsotmpice',           &
      'clcalipsotmpliq','clcalipsotmpun','cltlidarradar','pctisccp','dbze94',  &
      'tauisccp','cltisccp','toffset','parasolrefl','clmisr','meantbisccp',    &
      'meantbclrisccp','fracout','lidarbetamol532','tbrttov','cltmodis',       &
      'clwmodis','climodis','clhmodis','clmmodis','cllmodis','tautmodis',      &
      'tauwmodis','tauimodis','tautlogmodis','tauwlogmodis','tauilogmodis',    &
      'reffclwmodis','reffclimodis','pctmodis','lwpmodis','iwpmodis',          &
      'clmodis','crimodis','crlmodis'/)

    ! Initialize dim2, check requested vertical grid
    if (cosp_inputs%use_vgrid) then
      nlev = cosp_inputs%nlr
    else
      nlev = cosp_inputs%nlevels
    endif

    dim2 = (/ 0,0,0,0,nlev,nlev,nlev,nlev,0,cosp_npres,0,0,0,0,0,0,0,0,0,0,0,0,&
      0,0,0,nlev,nlev,nlev,lidar_ntemp,lidar_ntemp,lidar_ntemp,lidar_ntemp,0,  &
      0,0,0,0,0,parasol_nrefl,misr_n_cth,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  &
      0,0,0,0,cosp_npres,nummodisrefficebins,nummodisreffliqbins /)

    ! ! Finds the index based on a string search
    ! !   use to test out_list_flags, out_list_labels
    ! search_str = 'albisccp'
    ! print *, "The Index for '",TRIM(search_str),"'' is ",                      &
    !   minval(pack([(ii,ii=1,n_out_list)],out_list_labels.eq.search_str))

    ! Pass requested simulators to cfg
    ! -------------------------------------------------------------------------
    cfg%lradar_sim = lradar_sim
    cfg%llidar_sim = llidar_sim
    cfg%lisccp_sim = lisccp_sim
    cfg%lmodis_sim = lmodis_sim
    cfg%lmisr_sim  = lmisr_sim
    cfg%lrttov_sim = lrttov_sim

    ! Flag to control output to file
    ! -------------------------------------------------------------------------
    cfg%lwrite_output    = .false.
    if (cfg%lstats.or.cfg%lmisr_sim.or.cfg%lrttov_sim) then
        cfg%lwrite_output = .true.
    endif

    ! Pass requested outputs/diagnostics to cfg and cvar
    ! -------------------------------------------------------------------------
    do i = 1, n_out_list
      if (out_list_flags(i)) cfg%out_list(i) = out_list_labels(i)
#ifdef COSP_DEBUG
      if (.not. out_list_flags(i)) then
        write(0,fmt='("i = ",i2,1x,A," = ",L2," to ",A)')                      &
          i,trim(out_list_labels(i)(1:13)),out_list_flags(i),                  &
          trim(out_list_labels(i))
      endif
#endif
      if (cosp_units(i)=='      ') CYCLE
      cvar(i)%cflag = out_list_flags(i)
      cvar(i)%sname = cosp_snames(i)
      cvar(i)%lname = cosp_lnames(i)
      cvar(i)%units = cosp_units(i)
      cvar(i)%dim1  = dim1(i)
      cvar(i)%dim2  = dim2(i)
      if (dim1(i)+dim2(i).gt.0) then
        cvar(i)%is_3d = .true.
      else
        cvar(i)%is_3d = .false.
      endif
#ifdef COSP_DEBUG
      if (.not. cvar(i)%cflag) CYCLE
      if (cvar(i)%is_3d) then
        write(0,threed_fmt) i,TRIM(cvar(i)%sname),cvar(i)%cflag,cvar(i)%sname, &
          cvar(i)%lname,cvar(i)%units,cvar(i)%dim1,cvar(i)%dim2
      else
        write(0,twod_fmt) i,TRIM(cvar(i)%sname),cvar(i)%cflag,cvar(i)%sname,   &
          cvar(i)%lname,cvar(i)%units
      endif
#endif
    enddo

    ! i   |sname              | lname                                                     | units     | dim1          | dim2               | is_3d
    ! ---------------------------------------------------------------------------------------------------------------------------------------------
    ! 01  'Lalbisccp       ', 'ISCCP Mean Cloud Albedo                             '      '1     '    0               0                     .false.
    ! 02  'Latb532         ', '                                                    '      '      '    0               0                     .false.
    ! 03  'Lboxptopisccp   ', '                                                    '      '      '    0               0                     .false.
    ! 04  'Lboxtauisccp    ', '                                                    '      '      '    0               0                     .false.
    ! 05  'Lcfaddbze94     ', 'CloudSat Radar Reflectivity CFAD                    '      '1     '    dbze_bins       nlev                  .true.
    ! 06  'Lcfadlidarsr532 ', 'CALIPSO Scattering Ratio CFAD                       '      '1     '    sr_bins         nlev                  .true.
    ! 07  'Lclcalipso2     ', 'CALIPSO Cloud Fraction Undetected by CloudSat       '      '%     '    1               nlev                  .true.
    ! 08  'Lclcalipso      ', 'CALIPSO Cloud Fraction                              '      '%     '    1               nlev                  .true.
    ! 09  'Lclhcalipso     ', 'CALIPSO High Level Cloud Fraction                   '      '%     '    0               0                     .false.
    ! 10  'Lclisccp        ', 'ISCCP Cloud Area Fraction CFAD                      '      '%     '    isccp_ntau      isccp_npres           .true.
    ! 11  'Lcllcalipso     ', 'CALIPSO Low Level Cloud Fraction                    '      '%     '    0               0                     .false.
    ! 12  'Lclmcalipso     ', 'CALIPSO Mid Level Cloud Fraction                    '      '%     '    0               0                     .false.
    ! 13  'Lcltcalipso     ', 'CALIPSO Total Cloud Area Fraction                   '      '%     '    0               0                     .false.
    ! 14  'Lcllcalipsoice  ', 'CALIPSO Low Level Ice Cloud Fraction                '      '%     '    0               0                     .false.
    ! 15  'Lclmcalipsoice  ', 'CALIPSO Mid Level Ice Cloud Fraction                '      '%     '    0               0                     .false.
    ! 16  'Lclhcalipsoice  ', 'CALIPSO High Level Ice Cloud Fraction               '      '%     '    0               0                     .false.
    ! 17  'Lcltcalipsoice  ', 'CALIPSO Ice Total Cloud Fraction                    '      '%     '    0               0                     .false.
    ! 18  'Lcllcalipsoliq  ', 'CALIPSO Low Level Liquid Cloud Fraction             '      '%     '    0               0                     .false.
    ! 19  'Lclmcalipsoliq  ', 'CALIPSO Mid Level Liquid Cloud Fraction             '      '%     '    0               0                     .false.
    ! 20  'Lclhcalipsoliq  ', 'CALIPSO High Level Liquid Cloud Fraction            '      '%     '    0               0                     .false.
    ! 21  'Lcltcalipsoliq  ', 'CALIPSO Liquid Total Cloud Fraction                 '      '%     '    0               0                     .false.
    ! 22  'Lcllcalipsoun   ', 'CALIPSO Low Undefined-Phase Cloud Fraction          '      '%     '    0               0                     .false.
    ! 23  'Lclmcalipsoun   ', 'CALIPSO Mid Undefined-Phase Cloud Fraction          '      '%     '    0               0                     .false.
    ! 24  'Lclhcalipsoun   ', 'CALIPSO High Undefined-Phase Cloud Fraction         '      '%     '    0               0                     .false.
    ! 25  'Lcltcalipsoun   ', 'CALIPSO Undefined-Phase Total Cloud Fraction        '      '%     '    0               0                     .false.
    ! 26  'Lclcalipsoice   ', 'CALIPSO 3D Ice Cloud Fraction                       '      '%     '    1               nlev                  .true.
    ! 27  'Lclcalipsoliq   ', 'CALIPSO 3D Liquid Cloud Fraction                    '      '%     '    1               nlev                  .true.
    ! 28  'Lclcalipsoun    ', 'CALIPSO 3D Undefined-Phase Cloud Fraction           '      '%     '    1               nlev                  .true.
    ! 29  'Lclcalipsotmp   ', 'CALIPSO 3D Cloud Temperature                        '      'C     '    1               lidar_ntemp           .true.
    ! 30  'Lclcalipsotmpice', 'CALIPSO Ice Cloud Fraction Undetected by CloudSat   '      '%     '    1               lidar_ntemp           .true.
    ! 31  'Lclcalipsotmpliq', 'CALIPSO Liquid Cloud Fraction Undetected by CloudSat'      '%     '    1               lidar_ntemp           .true.
    ! 32  'Lclcalipsotmpun ', 'CALIPSO 3D Undefined-Phase Cloud Temperature        '      '%     '    1               lidar_ntemp           .true.
    ! 33  'Lcltlidarradar  ', 'Lidar and Radar Total Cloud Fraction                '      '%     '    0               0                     .false.
    ! 34  'Lpctisccp       ', 'ISCCP Mean Cloud Top Pressure                       '      'Pa    '    0               0                     .false.
    ! 35  'Ldbze94         ', '                                                    '      '      '    0               0                     .false.
    ! 36  'Ltauisccp       ', 'ISCCP Mean Optical Depth                            '      '1     '    0               0                     .false.
    ! 37  'Lcltisccp       ', 'ISCCP Total Cloud Fraction                          '      '%     '    0               0                     .false.
    ! 38  'Ltoffset        ', '                                                    '      '      '    0               0                     .false.
    ! 39  'Lparasolrefl    ', 'PARASOL Reflectance                                 '      '1     '    1               parasol_nrefl         .true.
    ! 40  'Lclmisr         ', 'MISR Cloud Fraction CFAD                            '      '%     '    isccp_ntau      misr_n_cth            .true.
    ! 41  'Lmeantbisccp    ', 'ISCCP Mean all-sky 10.5um brightness temp           '      'K     '    0               0                     .false.
    ! 42  'Lmeantbclrisccp ', 'ISCCP Mean clear-sky 10.5um brightness temp         '      'K     '    0               0                     .false.
    ! 43  'Lfracout        ', '                                                    '      '      '    0               0                     .false.
    ! 44  'Llidarbetamol532', '                                                    '      '      '    1               0                     .false.
    ! 45  'Ltbrttov        ', '                                                    '      '      '    0               0                     .false.
    ! 46  'Lcltmodis       ', 'MODIS Total Cloud Fraction                          '      '%     '    0               0                     .false.
    ! 47  'Lclwmodis       ', 'MODIS Liquid Cloud Fraction                         '      '%     '    0               0                     .false.
    ! 48  'Lclimodis       ', 'MODIS Ice Cloud Fraction                            '      '%     '    0               0                     .false.
    ! 49  'Lclhmodis       ', 'MODIS High Level Cloud Fraction                     '      '%     '    0               0                     .false.
    ! 50  'Lclmmodis       ', 'MODIS Mid Level Cloud Fraction                      '      '%     '    0               0                     .false.
    ! 51  'Lcllmodis       ', 'MODIS Low Level Cloud Fraction                      '      '%     '    0               0                     .false.
    ! 52  'Ltautmodis      ', 'MODIS Total Cloud Optical Thickness                 '      '1     '    0               0                     .false.
    ! 53  'Ltauwmodis      ', 'MODIS Liquid Cloud Optical Thickness                '      '1     '    0               0                     .false.
    ! 54  'Ltauimodis      ', 'MODIS Ice Cloud Optical Thickness                   '      '1     '    0               0                     .false.
    ! 55  'Ltautlogmodis   ', 'MODIS Total Cloud Optical Thickness (Log10 Mean)    '      '1     '    0               0                     .false.
    ! 56  'Ltauwlogmodis   ', 'MODIS Liquid Cloud Optical Thickness (Log10 Mean)   '      '1     '    0               0                     .false.
    ! 57  'Ltauilogmodis   ', 'MODIS Ice Cloud Optical Thickness (Log10 Mean)      '      '1     '    0               0                     .false.
    ! 58  'Lreffclwmodis   ', 'MODIS Liquid Cloud Particle Size                    '      'm     '    0               0                     .false.
    ! 59  'Lreffclimodis   ', 'MODIS Ice Cloud Particle Size                       '      'm     '    0               0                     .false.
    ! 60  'Lpctmodis       ', 'MODIS Cloud Top Pressure                            '      'Pa    '    0               0                     .false.
    ! 61  'Llwpmodis       ', 'MODIS Cloud Liquid Water Path                       '      'kg m-2'    0               0                     .false.
    ! 62  'Liwpmodis       ', 'MODIS Cloud Ice Water Path                          '      'kg m-2'    0               0                     .false.
    ! 63  'Lclmodis        ', 'MODIS PC-Tau Histogram                              '      '%     '    isccp_ntau      isccp_npres           .true.
    ! 64  'Lcrimodis       ', 'MODIS Optical_Thickness_vs_ReffIce Histogram        '      '%     '    isccp_ntau      nummodisrefficebins   .true.
    ! 65  'Lcrlmodis       ', 'MODIS Optical_Thickness_vs_ReffLiq Histogram        '      '%     '    isccp_ntau      nummodisreffliqbins   .true.

    ! Note the variable crimodis is sometimes referred to as jpdftaureicemodis.
    ! Likewise for crlmodis (jpdftaureliqmodis),

    ! Sync requested outputs/diagnostics to cfg
    ! -------------------------------------------------------------------------
    ! ISCCP simulator
    cfg%latb532          = latb532
    cfg%lboxptopisccp    = lboxptopisccp
    cfg%lboxtauisccp     = lboxtauisccp
    cfg%lmeantbisccp     = lmeantbisccp
    cfg%lmeantbclrisccp  = lmeantbclrisccp
    cfg%lclisccp         = lclisccp
    cfg%lpctisccp        = lpctisccp
    cfg%ltauisccp        = ltauisccp
    cfg%lcltisccp        = lcltisccp

    ! CloudSAT simulator
    cfg%ldbze94          = ldbze94
    cfg%lcfaddbze94      = lcfaddbze94

    ! CALIPSO/PARASOL simulator
    cfg%lcfadlidarsr532  = lcfadlidarsr532
    cfg%lclcalipso2      = lclcalipso2
    cfg%lclcalipso       = lclcalipso
    cfg%lclhcalipso      = lclhcalipso
    cfg%lcllcalipso      = lcllcalipso
    cfg%lclmcalipso      = lclmcalipso
    cfg%lcltcalipso      = lcltcalipso
    cfg%lclhcalipsoice   = lclhcalipsoice
    cfg%lcllcalipsoice   = lcllcalipsoice
    cfg%lclmcalipsoice   = lclmcalipsoice
    cfg%lcltcalipsoice   = lcltcalipsoice
    cfg%lclhcalipsoliq   = lclhcalipsoliq
    cfg%lcllcalipsoliq   = lcllcalipsoliq
    cfg%lclmcalipsoliq   = lclmcalipsoliq
    cfg%lcltcalipsoliq   = lcltcalipsoliq
    cfg%lclhcalipsoun    = lclhcalipsoun
    cfg%lcllcalipsoun    = lcllcalipsoun
    cfg%lclmcalipsoun    = lclmcalipsoun
    cfg%lcltcalipsoun    = lcltcalipsoun
    cfg%lclcalipsoice    = lclcalipsoice
    cfg%lclcalipsoliq    = lclcalipsoliq
    cfg%lclcalipsoun     = lclcalipsoun
    cfg%lclcalipsotmp    = lclcalipsotmp
    cfg%lclcalipsotmpice = lclcalipsotmpice
    cfg%lclcalipsotmpliq = lclcalipsotmpliq
    cfg%lclcalipsotmpun  = lclcalipsotmpun
    cfg%lcltlidarradar   = lcltlidarradar
    cfg%lparasolrefl     = lparasolrefl

    ! MISR simulator
    cfg%lclmisr          = lclmisr

    ! OTHER
    cfg%ltoffset         = ltoffset
    cfg%lfracout         = lfracout
    cfg%llidarbetamol532 = llidarbetamol532

    ! RTTOV
    cfg%ltbrttov         = ltbrttov

    ! MODIS simulator
    cfg%lcltmodis        = lcltmodis
    cfg%lclwmodis        = lclwmodis
    cfg%lclimodis        = lclimodis
    cfg%lclhmodis        = lclhmodis
    cfg%lclmmodis        = lclmmodis
    cfg%lcllmodis        = lcllmodis
    cfg%ltautmodis       = ltautmodis
    cfg%ltauwmodis       = ltauwmodis
    cfg%ltauimodis       = ltauimodis
    cfg%ltautlogmodis    = ltautlogmodis
    cfg%ltauwlogmodis    = ltauwlogmodis
    cfg%ltauilogmodis    = ltauilogmodis
    cfg%lreffclwmodis    = lreffclwmodis
    cfg%lreffclimodis    = lreffclimodis
    cfg%lpctmodis        = lpctmodis
    cfg%llwpmodis        = llwpmodis
    cfg%liwpmodis        = liwpmodis
    cfg%lclmodis         = lclmodis
    cfg%lcrimodis        = lcrimodis
    cfg%lcrlmodis        = lcrlmodis

#ifdef RTTOV
    call construct_cosp_rttov(cosp_inputs%npoints,cosp_inputs%nchannels,rttov)
#ifdef COSP_DEBUG
    call cosp_rttpv_print(rttov)
#endif
#endif

#ifdef COSP_DEBUG
    call cosp_output_print(cfg)
    write(0,fmt='(15x,a,/)') "end init_cosp()"

    ! Set min to large number
    min_dtau_c       = 1e30
    min_dem_c_cvcliq = 1e30
    min_dem_c_cvcice = 1e30
    min_reff_cvcliq  = 1e30
    min_mrhy_cvcliq  = 1e30
    min_reff_cvcice  = 1e30
    min_mrhy_cvcice  = 1e30
    min_reff_cvrain  = 1e30
    min_mrhy_cvrain  = 1e30
    min_reff_cvsnow  = 1e30
    min_mrhy_cvsnow  = 1e30
    min_dtau_s       = 1e30
    min_dem_s_lscliq = 1e30
    min_dem_s_lscice = 1e30
    min_reff_lscliq  = 1e30
    min_mrhy_lscliq  = 1e30
    min_reff_lscice  = 1e30
    min_mrhy_lscice  = 1e30
    min_reff_lsrain  = 1e30
    min_mrhy_lsrain  = 1e30
    min_reff_lssnow  = 1e30
    min_mrhy_lssnow  = 1e30

    ! Set Max to small number
    max_dtau_c       = -1e30
    max_dem_c_cvcliq = -1e30
    max_dem_c_cvcice = -1e30
    max_reff_cvcliq  = -1e30
    max_mrhy_cvcliq  = -1e30
    max_reff_cvcice  = -1e30
    max_mrhy_cvcice  = -1e30
    max_reff_cvrain  = -1e30
    max_mrhy_cvrain  = -1e30
    max_reff_cvsnow  = -1e30
    max_mrhy_cvsnow  = -1e30
    max_dtau_s       = -1e30
    max_dem_s_lscliq = -1e30
    max_dem_s_lscice = -1e30
    max_reff_lscliq  = -1e30
    max_mrhy_lscliq  = -1e30
    max_reff_lscice  = -1e30
    max_mrhy_lscice  = -1e30
    max_reff_lsrain  = -1e30
    max_mrhy_lsrain  = -1e30
    max_reff_lssnow  = -1e30
    max_mrhy_lssnow  = -1e30
#endif
    return
  end subroutine init_cosp

  ! ***************************************************************************
  ! PUBLIC SUBROUTINE: init_cosp_gbx
  ! ********************************
  subroutine init_cosp_gbx()
!@sum   Define and initialize COSP gridbox data structure for passing modelE
!@+     data to COSP.
!@auth  Mike Bauer
!@usage Run once at beginning of every time-step where the COSP simulators
!@+       will be called and before the collection of column/grid data.
!@      For example, near the bottom of CLOUDS_DRV.F90:CONDSE_column(i,j)
!0        inside the column physics i,j loop of ATM_DRV.f:atm_colphys().
!@calls construct_cosp_gridbox
    ! -------------------------------------------------------------------------
    ! Inputs:
    !
    ! Outputs:
    !   Initializes the COSP data structures gbx.
    !
    ! ToDo:
    !
    ! Language: Fortran 90
    ! -------------------------------------------------------------------------
    !
    ! Notes General:
    !   * ModelE does not allow 'components' such as COSP to access main
    !     modelE modules such as MODEL_COM. This results in the large
    !     argument list of this routine. This might change at some point
    !     in which case USE statements are preferred.
    !   * Layer data must be terrain following. That is, height is above
    !     sea level, not above the surface.
    !   * The rundeck pre-processor option COSP_PFLUX (sets flag_pfluxes), and
    !     the namelist parameter use_precipitation_fluxes, tell
    !     COSP to expect precipitation data as a flux, otherwise it expects
    !     a hydrometeor mixing ratio.
    !   * The rundeck pre-processor option COSP_USERE (sets flag_re), and the
    !     namelist parameter use_reff, tell COSP to expect model derived
    !     hydrometeor effective radii, otherwise COSP calculates them itself
    !     using distributions set in cosp_constants.F90 HCLASS and N_HYDRO.
    !
    ! Notes from icarus/isccp simulator:
    !   * Reversed vertical order (nLevels:1:-1) from gbx
    !       - p,ph,sh,tca,cca,dtau_c,dtau_s,t,dem_s,dem_c
    !   * pfull
    !       - index 1 is the top of the model
    !       - index nlev is the bottom of the model
    !   * phalf
    !       - index 1 is the top of model
    !       - index nlev+1 is the surface pressure
    !   * dtau_s, dtau_c, dem_s, dem_c
    !       - This the cloud optical depth of only the cloudy part of
    !         the grid box, it is not weighted with the 0 cloud optical
    !         depth of the clear part of the grid box
    !   * The following are used only if top_height = 1 or top_height = 3
    !       - skt, emsfc_lw, at, dem_s, dem_c, frac_out
    !       - The COSP default is top_height = 1
    !
    ! Notes from quickbeam/radar simulator:
    !   * These arrays must be in order from closest to farthest from the
    !     radar, for CFMIP we only use space-borne radar so index 1 is TOA.
    !       - p_matrix,hgt_matrix,t_matrix,rh_matrix,re_matrix
    !         From cosp_radar.f90
    !           p_matrix   = gbx%p/100.0     ! From Pa to hPa
    !           hgt_matrix = gbx%zlev/1000.0 ! From m to km
    !           t_matrix   = gbx%T-273.15    ! From K to C
    !           rh_matrix  = gbx%q
    !           re_matrix  = 0.0
    !
    ! Notes from MODIS simulator:
    !   * Reversed vertical order (nLevels:1:-1) from gbx
    !       - T,ph,ph,mr_hydro,dtau_c,reff
    !
    ! Notes from the lidar/calipso simulator:
    !   * presf: pressure half levels (as with gbx)
    !       - index 1 is sfc pressure
    !       - index nlev+1 TOA pressure
    !       inputs:
    !           real pres(npoints,nlev)    ! pressure full levels
    !           real presf(npoints,nlev+1) ! pressure half levels
    !       passed from cosp_lidar.f90
    !           presf(:,1:sgx%Nlevels) = gbx%ph
    !           presf(:,sgx%Nlevels + 1) = 0.0
    !
    ! Notes from the MISR simulator:
    !   * Reversed vertical order (nLevels:1:-1) from gbx
    !       - dtau_c,dtau_s,t,zlev
    !   * zfull
    !       - index 1 is top of model
    !       - index nlev is the bottom level of model
    !
    ! =========================================================================
    !         Details of the mapping of ModelE variables to COSP gbx
    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    ! Assume called from CLOUDS_DRV.F90:CONDSE_column(i,j) inside the column
    !   physics (i,j) loop of ATM_DRV.f:atm_colphys.
    !
    ! i       = grid/column longitude index
    ! j       = grid/column latitude index
    ! nlev    = number of model layers
    ! np_cosp = grid/column index being passed to COSP
    !
    !       Geo-location Information for each grid being processed by COSP
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! Longitude [0:360, degrees east, must be + signed]
    !   gbx%longitude(np_cosp)
    ! Latitude [-90:90, degrees north]
    !   gbx%latitude(np_cosp)
    !
    !           Surface Information for each grid being processed by COSP
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! Landmask [0 - ocean, 1 - land]
    !   gbx%land(np_cosp)
    ! Skin temperature (radiative temperature) [K]
    !   gbx%skt(np_cosp)
    ! Surface pressure [Pa]
    !   gbx%psfc(np_cosp)
    ! Sunlit [1 for day, 0 for night]
    !   gbx%sunlit(np_cosp)
    !
    !        Column Based Information for each grid being processed by COSP
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! Note: For COSP Levels are such at the 1st level is nearest the SFC.
    !       Also, half levels refer to lower layer edges (bottoms).
    ! Pressure at model levels [Pa, by layer]
    !   gbx%p(np_cosp,nlev)
    ! Pressure at half model levels (upper) [Pa]
    !   gbx%ph(np_cosp,nlev)
    ! Temperature at model levels [K]
    !   gbx%t(np_cosp,nlev)
    ! Relative humidity at model levels (%, to water)
    !   gbx%q(np_cosp,nlev)
    ! Specific humidity [kg/kg, to water]
    !   gbx%sh(np_cosp,nlev)
    ! Height of model levels [m]
    !   gbx%zlev(np_cosp,nlev)
    ! Height at half model levels [m]
    !   gbx%zlev_half(np_cosp,nlev)
    ! Total cloud fraction at model levels [unitless]
    !   gbx%tca(np_cosp,nlev)
    !
    !             Stratiform/Large-Scale Cloud Info processed by COSP
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !
    ! Large-Scale Cloud optical depth (mean 0.67 micron) [unitless]
    !   gbx%dtau_s(np_cosp,nlev)
    ! Large-Scale Cloud longwave emissivity (10.5 micron) [unitless]
    !   gbx%dem_s(np_cosp,nlev)
    ! Mixing ratio of large-scale cloud liquid [kg/kg]
    !   gbx%mr_hydro(np_cosp,nlev,I_LSCLIQ)
    ! Effective radius of large-scale cloud liquid [m]
    !   gbx%reff(np_cosp,nlev,I_LSCLIQ), only if COSP_USERE
    ! Mixing ratio of large-scale cloud ice [kg/kg]
    !   gbx%mr_hydro(np_cosp,nlev,I_LSCICE)
    ! Effective radius of large-scale cloud ice [m]
    !   gbx%reff(np_cosp,nlev,I_LSCICE), only if COSP_USERE
    ! Precipitation in ModleE is defined at layer edges,
    !   whereas COSP needs at layer levels (average of bounding edges).
    ! Large-scale rain
    !   if COSP_PFLUX
    !     Flux of large-scale rain [kg/m2.s]
    !       gbx%rain_ls(np_cosp,nlev)
    !   else
    !     Mixing ratio of large-scale rain [kg/kg]
    !       gbx%mr_hydro(np_cosp,nlev,I_LSRAIN)
    ! Effective radius of large-scale rain [m]
    !   gbx%reff(np_cosp,nlev,I_LSRAIN), only if COSP_USERE
    ! Large-scale ice
    !   if COSP_PFLUX
    !     Flux of large-scale snow [kg/m2.s]
    !       gbx%snow_ls(np_cosp,nlev)
    !   else
    !     Mixing ratio of large-scale snow [kg/kg]
    !       gbx%mr_hydro(np_cosp,nlev,I_LSSNOW)
    ! Effective radius of large-scale snow [m]
    !   gbx%reff(np_cosp,nlev,I_LSSNOW), only if COSP_USERE
    ! Large-scale graupel
    !   Not defined in ModelE
    !   if COSP_PFLUX
    !     Flux of large-scale graupel [kg/m2.s]
    !       gbx%grpl_ls(np_cosp,nlev)
    !   else
    !     Mixing ratio of large-scale graupel [kg/kg]
    !     gbx%mr_hydro(np_cosp,nlev,I_LSGRPL)
    ! Effective radius of large-scale graupel [m]
    !   gbx%reff(np_cosp,nlev,I_LSGRPL)
    !
    !                 Convective Cloud Info processed by COSP
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !
    ! Convective cloud fraction at model levels [unitless]
    !   gbx%cca(np_cosp,nlev)
    ! Convective Cloud optical depth (mean 0.67 micron) [unitless]
    !   gbx%dtau_c(np_cosp,nlev)
    ! Convective Cloud longwave emissivity (10.5 micron) [unitless]
    !   gbx%dem_c(np_cosp,nlev)
    ! Mixing ratio of convective cloud liquid [kg/kg]
    !   gbx%mr_hydro(np_cosp,nlev,I_CVCLIQ)
    ! Effective radius of convective cloud liquid [m]
    !   gbx%reff(np_cosp,nlev,I_CVCLIQ), only if COSP_USERE
    ! Mixing ratio of convective cloud ice [kg/kg]
    !   gbx%mr_hydro(np_cosp,nlev,I_CVCICE)
    ! Effective radius of convective cloud ice [m]
    !   gbx%reff(np_cosp,nlev,I_CVCICE), only if COSP_USERE
    ! Precipitation in ModleE is defined at layer edges,
    !   whereas COSP needs at layer levels (average of bounding edges).
    ! Convective rain
    !   if COSP_PFLUX
    !     Flux of convective rain (average crossing layer edges) [kg/m2.s]
    !       gbx%rain_cv(np_cosp,nlev)
    !   else
    !      Mixing ratio of convective rain (average crossing layer edges) [kg/kg]
    !       gbx%mr_hydro(np_cosp,nlev,I_CVRAIN)
    ! Effective radius of convective rain [m]
    !   gbx%reff(np_cosp,nlev,I_CVRAIN), only if COSP_USERE
    ! Convective snow
    !   if COSP_PFLUX
    !     Mixing ratio of convective snow (average crossing layer edges) [kg/kg]
    !       gbx%snow_cv(np_cosp,nlev)
    !   else
    !     Mixing ratio of convective snow [kg/kg]
    !       gbx%mr_hydro(np_cosp,nlev,I_CVSNOW)
    ! Effective radius of convective snow [m]
    !   gbx%reff(np_cosp,nlev,I_CVSNOW), only if COSP_USERE
    !--------------------------------------------------------------------------

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! Subroutine Scalers/Variables
    ! =========================================================================
    ! 10.5 micron emissivity of surface (fraction)
    ! Same value that the one originally used in the ISCCP simulator
    real :: isccp_emsfc_lw = 0.99
    double precision :: cosp_time

    ! Subroutine Arrays
    ! =========================================================================
    double precision :: cosp_time_bnds(2)
    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! Allocate, and check consistency, for COSP gridbox type (gbx)
    ! -----------------------------------------------------------------------
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_gridbox()"
#endif
    call construct_cosp_gridbox(cosp_time,cosp_time_bnds,                      &
      cosp_inputs%radar_freq,cosp_inputs%surface_radar,                        &
      cosp_inputs%use_mie_tables,cosp_inputs%use_gas_abs,cosp_inputs%do_ray,   &
      cosp_inputs%melt_lay,cosp_inputs%k2,cosp_inputs%npoints,                 &
      cosp_inputs%nlevels,cosp_inputs%ncolumns,n_hydro,                        &
      cosp_inputs%nprmts_max_hydro,cosp_inputs%naero,                          &
      cosp_inputs%nprmts_max_aero,cosp_inputs%npoints_it,                      &
      cosp_inputs%lidar_ice_type,cosp_inputs%isccp_topheight,                  &
      cosp_inputs%isccp_topheight_direction,cosp_inputs%overlap,               &
      isccp_emsfc_lw,cosp_inputs%use_precipitation_fluxes,                     &
      cosp_inputs%use_reff,cosp_inputs%platform,cosp_inputs%satellite,         &
      cosp_inputs%instrument,cosp_inputs%nchannels,cosp_inputs%zenang,         &
      cosp_inputs%channels(1:cosp_inputs%nchannels),                           &
      cosp_inputs%surfem(1:cosp_inputs%nchannels),cosp_inputs%co2,             &
      cosp_inputs%ch4,cosp_inputs%n2o,cosp_inputs%co,gbx)

    return
  end subroutine init_cosp_gbx

  ! ***************************************************************************
  ! PUBLIC SUBROUTINE: run_cosp_sims
  ! ********************************
  subroutine run_cosp_sims()
!@sum   Initialize the COSP simulators, pass in the model data, and run.
!@auth  Mike Bauer
!@usage Run once for every time-step where the COSP simulators are asked for,
!@+       after the looping over all column/grids for that domain.
!@      For example, after the column physics i,j loop in
!@+       ATM_DRV.f:atm_colphys().
!@calls cosp,construct_cosp_vgrid,construct_cosp_subgrid,construct_cosp_sgradar,
!@+     construct_cosp_radarstats,construct_cosp_sglidar,construct_cosp_isccp,
!@+     construct_cosp_misr,construct_cosp_modis
    ! -------------------------------------------------------------------------
    ! Inputs:
    !
    ! Outputs:
    !   Results of the COSP simulators (isccp,stlidar,misr,modis,stradar)
    !
    ! ToDo:
    !
    ! Language: Fortran 90
    ! -------------------------------------------------------------------------
    !
    ! Extra initialization needed for the COSP simulators.
    !  * Defines a new vertical grid to report COSP output (vgrid).
    !   * See use_vgrid in COSPIN.
    !   * Requires knowledge of model vertical grid (gbx%zlev, gbx%zlev_half).
    !     ! Height of model levels and mid-levels [m]
    !     gbx%zlev(np_cosp,:) = bygrav*gz(I,J,:)
    !     gbx%zlev_half(np_cosp,L)=0.5*bygrav*(gz(I,J,L)+gz(I,J,L-1))
    !  * Allocates and initialize the COSP subgrid (sgx)
    !  * Allocate and initialize simulator data structures
    !     (sgradar, sglidar, isccp, misr, modis)
    !  * Allocate and initialize statistics data structures
    !     (stradar, stlidar)
    !--------------------------------------------------------------------------

    ! MODULE Imports
    ! =========================================================================

    ! Imported Routines:
    use mod_cosp, only : cosp
    use mod_cosp_modis_simulator, only : construct_cosp_modis

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#ifdef COSP_DEBUG
    ! write(0,fmt='(/,15x,A)') "CV Cloud"

    ! write(0,fmt='(20x,A)') "dtau_c: Optical Depth"
    ! write(0,fmt='(20x,"Min = ",E10.3," Max = ",E10.3)') min_dtau_c,max_dtau_c

    ! write(0,fmt='(/,20x,A)') "CV Cloud Water"
    ! write(0,fmt='(25x,A)') "Mixing Ratio"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_mrhy_cvcliq,max_mrhy_cvcliq
    ! write(0,fmt='(25x,A)') "Longwave emissivity (dem_c)"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_dem_c_cvcliq,max_dem_c_cvcliq
    ! write(0,fmt='(25x,A)') "Effective Radius"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_reff_cvcliq,max_reff_cvcliq

    ! write(0,fmt='(/,20x,A)') "CV Cloud Ice"
    ! write(0,fmt='(25x,A)') "Mixing Ratio"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_mrhy_cvcice,max_mrhy_cvcice
    ! write(0,fmt='(25x,A)') "Longwave emissivity (dem_c)"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_dem_c_cvcice,max_dem_c_cvcice
    ! write(0,fmt='(25x,A)') "Effective Radius"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_reff_cvcice,max_reff_cvcice

    ! write(0,fmt='(/,20x,A)') "CV Rain"
    ! write(0,fmt='(25x,A)') "Mixing Ratio"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_mrhy_cvrain,max_mrhy_cvrain
    ! write(0,fmt='(25x,A)') "Effective Radius"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_reff_cvrain,max_reff_cvrain

    ! write(0,fmt='(/,20x,A)') "CV Snow"
    ! write(0,fmt='(25x,A)') "Mixing Ratio"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_mrhy_cvsnow,max_mrhy_cvsnow
    ! write(0,fmt='(25x,A)') "Effective Radius"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_reff_cvsnow,max_reff_cvsnow

    ! write(0,fmt='(/,15x,A)') "LS Cloud"

    ! write(0,fmt='(20x,A)') "dtau_s: Optical Depth"
    ! write(0,fmt='(20x,"Min = ",E10.3,", Max = ",E10.3)') min_dtau_s,max_dtau_s

    ! write(0,fmt='(/,20x,A)') "LS Cloud Water"
    ! write(0,fmt='(25x,A)') "Mixing Ratio"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_mrhy_lscliq,max_mrhy_lscliq
    ! write(0,fmt='(25x,A)') "Longwave emissivity (dem_s)"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_dem_s_lscliq,max_dem_s_lscliq
    ! write(0,fmt='(25x,A)') "Effective Radius"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_reff_lscliq,max_reff_lscliq

    ! write(0,fmt='(/,20x,A)') "LS Cloud Ice"
    ! write(0,fmt='(25x,A)') "Mixing Ratio"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_mrhy_lscice,max_mrhy_lscice
    ! write(0,fmt='(25x,A)') "Longwave emissivity (dem_s)"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_dem_s_lscice,max_dem_s_lscice
    ! write(0,fmt='(25x,A)') "Effective Radius"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_reff_lscice,max_reff_lscice

    ! write(0,fmt='(/,20x,A)') "LS Rain"
    ! write(0,fmt='(25x,A)') "Mixing Ratio"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_mrhy_lsrain,max_mrhy_lsrain
    ! write(0,fmt='(25x,A)') "Effective Radius"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_reff_lsrain,max_reff_lsrain

    ! write(0,fmt='(/,20x,A)') "LS Snow"
    ! write(0,fmt='(25x,A)') "Mixing Ratio"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_mrhy_lssnow,max_mrhy_lssnow
    ! write(0,fmt='(25x,A)') "Effective Radius"
    ! write(0,fmt='(25x,"Min = ",E10.3,", Max = ",E10.3)') min_reff_lssnow,max_reff_lssnow
#endif

    ! Define and initialize COSP vertical grid (vgrid)
    ! ------------------------------------------------------------------------
    call construct_cosp_vgrid(gbx,cosp_inputs%nlr,cosp_inputs%use_vgrid,      &
      cosp_inputs%csat_vgrid,vgrid)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_vgrid()"
    !call cosp_vgrid_print(vgrid)
#endif

    ! Allocates and initialize the COSP subgrid (sgx)
    ! ------------------------------------------------------------------------
    call construct_cosp_subgrid(cosp_inputs%npoints,cosp_inputs%Ncolumns,    &
      cosp_inputs%nlevels,sgx)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_subgrid()"
    !call cosp_subgrid_print(sgx)
#endif

    ! Allocate and initialize ISCCP output data structures (isccp)
    ! ------------------------------------------------------------------------
    call construct_cosp_isccp(cfg,cosp_inputs%npoints,cosp_inputs%ncolumns,  &
      cosp_inputs%nlevels,isccp)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_isccp()"
    !call cosp_isccp_print(isccp)
#endif

    ! Allocate and initialize radar output data structures (sgradar)
    ! ------------------------------------------------------------------------
    call construct_cosp_sgradar(cfg,cosp_inputs%npoints,cosp_inputs%ncolumns, &
      cosp_inputs%nlevels,n_hydro,sgradar)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_sgradar()"
    !call cosp_sgradar_print(sgradar)
#endif

    ! Allocate and initialize radar statistics data structures (stradar)
    ! ------------------------------------------------------------------------
    call construct_cosp_radarstats(cfg,cosp_inputs%npoints,                  &
      cosp_inputs%ncolumns,vgrid%nlvgrid,n_hydro,stradar)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_radarstats()"
    !call cosp_radarstats_print(stradar)
#endif

    ! Allocate and initialize lidar output data structures (sglidar)
    ! ------------------------------------------------------------------------
    call construct_cosp_sglidar(cfg,cosp_inputs%npoints,cosp_inputs%ncolumns, &
      cosp_inputs%nlevels,n_hydro,parasol_nrefl,sglidar)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_sglidar()"
    !call cosp_sglidar_print(sglidar)
#endif

    ! Allocate and initialize lidar statistics data structures (stlidar)
    ! ------------------------------------------------------------------------
    call construct_cosp_lidarstats(cfg,cosp_inputs%npoints,                  &
      cosp_inputs%ncolumns,vgrid%nlvgrid,n_hydro,parasol_nrefl,stlidar)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_lidarstats()"
    !call cosp_lidarstats_print(stlidar)
#endif

    ! Allocate and initialize MISR output data structures (misr)
    ! ------------------------------------------------------------------------
    call construct_cosp_misr(cfg,cosp_inputs%npoints,misr)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_misr()"
    !call cosp_misr_print(misr)
#endif

    ! Allocate and initialize MODIS output data structures (modis)
    ! ------------------------------------------------------------------------
    call construct_cosp_modis(cfg,cosp_inputs%npoints,modis)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_modis()"
    !!Not defined in COSP: call cosp_modis_print(modis)
#endif

    ! Allocate and initialize RTTOV data structures
    ! ------------------------------------------------------------------------
#ifdef RTTOV
    !call construct_cosp_rttov(cosp_inputs%npoints,cosp_inputs%nchannels,rttov)
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call construct_cosp_rttov()"
    !!Not defined in COSP: call cosp_rttpv_print(rttov)
#endif
#endif

#ifdef COSP_DEBUG
    !call cosp_gridbox_print(gbx)
    !write(0,fmt='(15x,A)') "call cosp()"
#endif

    ! Execute the COSP simulators; after COSP data structures gather data
    ! -------------------------------------------------------------------------
    call cosp(cosp_inputs%overlap,cosp_inputs%ncolumns,cfg,vgrid,gbx,sgx,      &
      sgradar,sglidar,isccp,misr,modis,stradar,stlidar)

    return
  end subroutine run_cosp_sims

  ! ***************************************************************************
  ! PUBLIC SUBROUTINE: free_cosp_sims
  ! *********************************
  subroutine free_cosp_sims()
!@sum   Free memory allocated by init_cosp_gbx() and run_cosp_sims().
!@auth  Mike Bauer
!@usage Run once for every time-step where the COSP simulators are asked for,
!@+       after run_cosp_sims() AND after saving the cosp data save_cosp().
!@      For example, after the column physics i,j loop in
!@+       ATM_DRV.f:atm_colphys().
!------------------------------------------------------------------------------

    ! MODULE Imports
    ! =========================================================================

    ! Imported Routines:
    use mod_cosp_modis_simulator, only : free_cosp_modis

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! Deallocate memory in derived types
    ! -------------------------------------------------------------------------
#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_gridbox()"
#endif
    call free_cosp_gridbox(gbx)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_subgrid()"
#endif
    call free_cosp_subgrid(sgx)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_sgradar()"
#endif
    call free_cosp_sgradar(sgradar)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_radarstats()"
#endif
    call free_cosp_radarstats(stradar)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_sglidar()"
#endif
    call free_cosp_sglidar(sglidar)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_lidarstats()"
#endif
    call free_cosp_lidarstats(stlidar)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_isccp()"
#endif
    call free_cosp_isccp(isccp)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_misr()"
#endif
    call free_cosp_misr(misr)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_modis()"
#endif
    call free_cosp_modis(modis)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "call free_cosp_vgrid()"
#endif
    call free_cosp_vgrid(vgrid)

#ifdef COSP_DEBUG
    !write(0,fmt='(15x,A)') "END free_cosp_sims()"
#endif
    return
  end subroutine free_cosp_sims

  ! ***************************************************************************
  ! PUBLIC SUBROUTINE: remap_cosp_2d
  ! ********************************
  subroutine remap_cosp_2d(vname,imaxj,geom_2d,blob)
!@sum  Maps native 2D COSP variables so suitable for the SUBDD framework.
!@auth Mike Bauer
    ! -------------------------------------------------------------------------
    ! Inputs:
    !   vname   : Name of COSP variable to be processed
    !   imaxj   : Varying number of used longitudes (See GEOM_B.f)
    !   geom_2d : Domain decomposition parameters (See Below)
    !             (I_0, I_1, J_0, J_1, I_0H, I_1H, J_0H, J_1H)
    !              1    2    3    4    5     6     7     8
    !   blob    : Allocatable ModelE data array (2D)
    !
    ! Outputs:
    !   blob  : 2D ModelE data array (I_0H:I_1H, J_0H:J_1H)
    ! Terminology:
    !
    ! Incoming COSP data are generally of size (nsubdd_cosp), which does not
    ! include the HALO and only has a single longitude for each
    ! of the polar latitudes. That is,
    !
    !   nsubdd_cosp = (1 + I_1 - I_0) * (1 + J_1 - J_0)
    !
    ! or if the domain contains a pole then
    !
    !   nsubdd_cosp = nsubdd_cosp - IM + 1
    !   nsubdd_cosp = (IM * JM) - (2 * IM) + 2
    !
    !   if (cosp_i .eq. 1) then
    !     ! Polar Cap (SH)
    !     j = 1
    !     i = 1
    !   else if (cosp_i .eq. nsubdd_cosp) then
    !       ! Polar Cap (NH)
    !       j = jm
    !       i = 1
    !   else
    !       j = floor(real(cosp_i - 2) / real(im)) + 2
    !       i = modulo(cosp_i - 2, im) + 1
    !   endif
    !
    ! This routine takes the COSP 2D array and remaps it to fill the ModelE
    ! domain decomposition outlined by geom_2d. This mapping basically returns
    ! the HALO and fills the polar rows with the single valid value from that
    ! latitude from COSP.
    !
    !     Domain Decomposition of the Earth System Modeling Framework (ESMF)
    !     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! see model/MPI_Support/DomainDecomposition_mod.F90
    ! -------------------------------------------------------------------------
    !   Parameters for Global domain
    !     IM_WORLD : Number of Longitudes
    !     JM_WORLD : Number of latitudes
    !   Parameters for local domain
    !     I_0  = I_STRT    : Start of the local domain longitude index
    !     I_1  = I_STOP    : End of the local domain longitude index
    !     J_0  = J_STRT    : Start of the of local domain latitude index
    !     J_1  = J_STOP    : End of the local domain latitude index
    !     J_0S = J_STRT_SKP: Start of the of local domain exclusive of S pole
    !     J_1S = J_STOP_SKP: End of the local domain exclusive of N pole
    !   Parameters for halo of local domain
    !     I_0H = I_STRT_HALO : Start of the halo longitude index
    !     I_1H = I_STOP_HALO : End of the halo longitude index
    !     J_0H = J_STRT_HALO : Start of the halo latitude index
    !     J_1H = J_STOP_HALO : End of the halo latitude index
    !
    ! Example (JM = 90, IM = 144, ncosp = 12,674)
    ! ncosp =     1 (i =   1, j =  1)  (-178.75,-90.00)      ...
    ! ncosp =     2 (i =   1, j =  2)  (-178.75,-87.00)      ...    ncosp =   145 (i = 144, j =  2)  ( 178.75,-87.00)
    ! ncosp =   146 (i =   1, j =  3)  (-178.75,-85.00)      ...    ncosp =   289 (i = 144, j =  3)  ( 178.75,-85.00)
    ! ncosp =   290 (i =   1, j =  4)  (-178.75,-83.00)      ...    ncosp =   433 (i = 144, j =  4)  ( 178.75,-83.00)
    ! ncosp =   434 (i =   1, j =  5)  (-178.75,-81.00)      ...    ncosp =   577 (i = 144, j =  5)  ( 178.75,-81.00)
    ! ncosp =   578 (i =   1, j =  6)  (-178.75,-79.00)      ...    ncosp =   721 (i = 144, j =  6)  ( 178.75,-79.00)
    ! ncosp =   722 (i =   1, j =  7)  (-178.75,-77.00)      ...    ncosp =   865 (i = 144, j =  7)  ( 178.75,-77.00)
    ! ncosp =   866 (i =   1, j =  8)  (-178.75,-75.00)      ...    ncosp =  1009 (i = 144, j =  8)  ( 178.75,-75.00)
    ! ncosp =  1010 (i =   1, j =  9)  (-178.75,-73.00)      ...    ncosp =  1153 (i = 144, j =  9)  ( 178.75,-73.00)
    ! ncosp =  1154 (i =   1, j = 10)  (-178.75,-71.00)      ...    ncosp =  1297 (i = 144, j = 10)  ( 178.75,-71.00)
    ! ncosp =  1298 (i =   1, j = 11)  (-178.75,-69.00)      ...    ncosp =  1441 (i = 144, j = 11)  ( 178.75,-69.00)
    ! ncosp =  1442 (i =   1, j = 12)  (-178.75,-67.00)      ...    ncosp =  1585 (i = 144, j = 12)  ( 178.75,-67.00)
    ! ncosp =  1586 (i =   1, j = 13)  (-178.75,-65.00)      ...    ncosp =  1729 (i = 144, j = 13)  ( 178.75,-65.00)
    ! ncosp =  1730 (i =   1, j = 14)  (-178.75,-63.00)      ...    ncosp =  1873 (i = 144, j = 14)  ( 178.75,-63.00)
    ! ncosp =  1874 (i =   1, j = 15)  (-178.75,-61.00)      ...    ncosp =  2017 (i = 144, j = 15)  ( 178.75,-61.00)
    ! ncosp =  2018 (i =   1, j = 16)  (-178.75,-59.00)      ...    ncosp =  2161 (i = 144, j = 16)  ( 178.75,-59.00)
    ! ncosp =  2162 (i =   1, j = 17)  (-178.75,-57.00)      ...    ncosp =  2305 (i = 144, j = 17)  ( 178.75,-57.00)
    ! ncosp =  2306 (i =   1, j = 18)  (-178.75,-55.00)      ...    ncosp =  2449 (i = 144, j = 18)  ( 178.75,-55.00)
    ! ncosp =  2450 (i =   1, j = 19)  (-178.75,-53.00)      ...    ncosp =  2593 (i = 144, j = 19)  ( 178.75,-53.00)
    ! ncosp =  2594 (i =   1, j = 20)  (-178.75,-51.00)      ...    ncosp =  2737 (i = 144, j = 20)  ( 178.75,-51.00)
    ! ncosp =  2738 (i =   1, j = 21)  (-178.75,-49.00)      ...    ncosp =  2881 (i = 144, j = 21)  ( 178.75,-49.00)
    ! ncosp =  2882 (i =   1, j = 22)  (-178.75,-47.00)      ...    ncosp =  3025 (i = 144, j = 22)  ( 178.75,-47.00)
    ! ncosp =  3026 (i =   1, j = 23)  (-178.75,-45.00)      ...    ncosp =  3169 (i = 144, j = 23)  ( 178.75,-45.00)
    ! ncosp =  3170 (i =   1, j = 24)  (-178.75,-43.00)      ...    ncosp =  3313 (i = 144, j = 24)  ( 178.75,-43.00)
    ! ncosp =  3314 (i =   1, j = 25)  (-178.75,-41.00)      ...    ncosp =  3457 (i = 144, j = 25)  ( 178.75,-41.00)
    ! ncosp =  3458 (i =   1, j = 26)  (-178.75,-39.00)      ...    ncosp =  3601 (i = 144, j = 26)  ( 178.75,-39.00)
    ! ncosp =  3602 (i =   1, j = 27)  (-178.75,-37.00)      ...    ncosp =  3745 (i = 144, j = 27)  ( 178.75,-37.00)
    ! ncosp =  3746 (i =   1, j = 28)  (-178.75,-35.00)      ...    ncosp =  3889 (i = 144, j = 28)  ( 178.75,-35.00)
    ! ncosp =  3890 (i =   1, j = 29)  (-178.75,-33.00)      ...    ncosp =  4033 (i = 144, j = 29)  ( 178.75,-33.00)
    ! ncosp =  4034 (i =   1, j = 30)  (-178.75,-31.00)      ...    ncosp =  4177 (i = 144, j = 30)  ( 178.75,-31.00)
    ! ncosp =  4178 (i =   1, j = 31)  (-178.75,-29.00)      ...    ncosp =  4321 (i = 144, j = 31)  ( 178.75,-29.00)
    ! ncosp =  4322 (i =   1, j = 32)  (-178.75,-27.00)      ...    ncosp =  4465 (i = 144, j = 32)  ( 178.75,-27.00)
    ! ncosp =  4466 (i =   1, j = 33)  (-178.75,-25.00)      ...    ncosp =  4609 (i = 144, j = 33)  ( 178.75,-25.00)
    ! ncosp =  4610 (i =   1, j = 34)  (-178.75,-23.00)      ...    ncosp =  4753 (i = 144, j = 34)  ( 178.75,-23.00)
    ! ncosp =  4754 (i =   1, j = 35)  (-178.75,-21.00)      ...    ncosp =  4897 (i = 144, j = 35)  ( 178.75,-21.00)
    ! ncosp =  4898 (i =   1, j = 36)  (-178.75,-19.00)      ...    ncosp =  5041 (i = 144, j = 36)  ( 178.75,-19.00)
    ! ncosp =  5042 (i =   1, j = 37)  (-178.75,-17.00)      ...    ncosp =  5185 (i = 144, j = 37)  ( 178.75,-17.00)
    ! ncosp =  5186 (i =   1, j = 38)  (-178.75,-15.00)      ...    ncosp =  5329 (i = 144, j = 38)  ( 178.75,-15.00)
    ! ncosp =  5330 (i =   1, j = 39)  (-178.75,-13.00)      ...    ncosp =  5473 (i = 144, j = 39)  ( 178.75,-13.00)
    ! ncosp =  5474 (i =   1, j = 40)  (-178.75,-11.00)      ...    ncosp =  5617 (i = 144, j = 40)  ( 178.75,-11.00)
    ! ncosp =  5618 (i =   1, j = 41)  (-178.75, -9.00)      ...    ncosp =  5761 (i = 144, j = 41)  ( 178.75, -9.00)
    ! ncosp =  5762 (i =   1, j = 42)  (-178.75, -7.00)      ...    ncosp =  5905 (i = 144, j = 42)  ( 178.75, -7.00)
    ! ncosp =  5906 (i =   1, j = 43)  (-178.75, -5.00)      ...    ncosp =  6049 (i = 144, j = 43)  ( 178.75, -5.00)
    ! ncosp =  6050 (i =   1, j = 44)  (-178.75, -3.00)      ...    ncosp =  6193 (i = 144, j = 44)  ( 178.75, -3.00)
    ! ncosp =  6194 (i =   1, j = 45)  (-178.75, -1.00)      ...    ncosp =  6337 (i = 144, j = 45)  ( 178.75, -1.00)
    ! ncosp =  6338 (i =   1, j = 46)  (-178.75,  1.00)      ...    ncosp =  6481 (i = 144, j = 46)  ( 178.75,  1.00)
    ! ncosp =  6482 (i =   1, j = 47)  (-178.75,  3.00)      ...    ncosp =  6625 (i = 144, j = 47)  ( 178.75,  3.00)
    ! ncosp =  6626 (i =   1, j = 48)  (-178.75,  5.00)      ...    ncosp =  6769 (i = 144, j = 48)  ( 178.75,  5.00)
    ! ncosp =  6770 (i =   1, j = 49)  (-178.75,  7.00)      ...    ncosp =  6913 (i = 144, j = 49)  ( 178.75,  7.00)
    ! ncosp =  6914 (i =   1, j = 50)  (-178.75,  9.00)      ...    ncosp =  7057 (i = 144, j = 50)  ( 178.75,  9.00)
    ! ncosp =  7058 (i =   1, j = 51)  (-178.75, 11.00)      ...    ncosp =  7201 (i = 144, j = 51)  ( 178.75, 11.00)
    ! ncosp =  7202 (i =   1, j = 52)  (-178.75, 13.00)      ...    ncosp =  7345 (i = 144, j = 52)  ( 178.75, 13.00)
    ! ncosp =  7346 (i =   1, j = 53)  (-178.75, 15.00)      ...    ncosp =  7489 (i = 144, j = 53)  ( 178.75, 15.00)
    ! ncosp =  7490 (i =   1, j = 54)  (-178.75, 17.00)      ...    ncosp =  7633 (i = 144, j = 54)  ( 178.75, 17.00)
    ! ncosp =  7634 (i =   1, j = 55)  (-178.75, 19.00)      ...    ncosp =  7777 (i = 144, j = 55)  ( 178.75, 19.00)
    ! ncosp =  7778 (i =   1, j = 56)  (-178.75, 21.00)      ...    ncosp =  7921 (i = 144, j = 56)  ( 178.75, 21.00)
    ! ncosp =  7922 (i =   1, j = 57)  (-178.75, 23.00)      ...    ncosp =  8065 (i = 144, j = 57)  ( 178.75, 23.00)
    ! ncosp =  8066 (i =   1, j = 58)  (-178.75, 25.00)      ...    ncosp =  8209 (i = 144, j = 58)  ( 178.75, 25.00)
    ! ncosp =  8210 (i =   1, j = 59)  (-178.75, 27.00)      ...    ncosp =  8353 (i = 144, j = 59)  ( 178.75, 27.00)
    ! ncosp =  8354 (i =   1, j = 60)  (-178.75, 29.00)      ...    ncosp =  8497 (i = 144, j = 60)  ( 178.75, 29.00)
    ! ncosp =  8498 (i =   1, j = 61)  (-178.75, 31.00)      ...    ncosp =  8641 (i = 144, j = 61)  ( 178.75, 31.00)
    ! ncosp =  8642 (i =   1, j = 62)  (-178.75, 33.00)      ...    ncosp =  8785 (i = 144, j = 62)  ( 178.75, 33.00)
    ! ncosp =  8786 (i =   1, j = 63)  (-178.75, 35.00)      ...    ncosp =  8929 (i = 144, j = 63)  ( 178.75, 35.00)
    ! ncosp =  8930 (i =   1, j = 64)  (-178.75, 37.00)      ...    ncosp =  9073 (i = 144, j = 64)  ( 178.75, 37.00)
    ! ncosp =  9074 (i =   1, j = 65)  (-178.75, 39.00)      ...    ncosp =  9217 (i = 144, j = 65)  ( 178.75, 39.00)
    ! ncosp =  9218 (i =   1, j = 66)  (-178.75, 41.00)      ...    ncosp =  9361 (i = 144, j = 66)  ( 178.75, 41.00)
    ! ncosp =  9362 (i =   1, j = 67)  (-178.75, 43.00)      ...    ncosp =  9505 (i = 144, j = 67)  ( 178.75, 43.00)
    ! ncosp =  9506 (i =   1, j = 68)  (-178.75, 45.00)      ...    ncosp =  9649 (i = 144, j = 68)  ( 178.75, 45.00)
    ! ncosp =  9650 (i =   1, j = 69)  (-178.75, 47.00)      ...    ncosp =  9793 (i = 144, j = 69)  ( 178.75, 47.00)
    ! ncosp =  9794 (i =   1, j = 70)  (-178.75, 49.00)      ...    ncosp =  9937 (i = 144, j = 70)  ( 178.75, 49.00)
    ! ncosp =  9938 (i =   1, j = 71)  (-178.75, 51.00)      ...    ncosp = 10081 (i = 144, j = 71)  ( 178.75, 51.00)
    ! ncosp = 10082 (i =   1, j = 72)  (-178.75, 53.00)      ...    ncosp = 10225 (i = 144, j = 72)  ( 178.75, 53.00)
    ! ncosp = 10226 (i =   1, j = 73)  (-178.75, 55.00)      ...    ncosp = 10369 (i = 144, j = 73)  ( 178.75, 55.00)
    ! ncosp = 10370 (i =   1, j = 74)  (-178.75, 57.00)      ...    ncosp = 10513 (i = 144, j = 74)  ( 178.75, 57.00)
    ! ncosp = 10514 (i =   1, j = 75)  (-178.75, 59.00)      ...    ncosp = 10657 (i = 144, j = 75)  ( 178.75, 59.00)
    ! ncosp = 10658 (i =   1, j = 76)  (-178.75, 61.00)      ...    ncosp = 10801 (i = 144, j = 76)  ( 178.75, 61.00)
    ! ncosp = 10802 (i =   1, j = 77)  (-178.75, 63.00)      ...    ncosp = 10945 (i = 144, j = 77)  ( 178.75, 63.00)
    ! ncosp = 10946 (i =   1, j = 78)  (-178.75, 65.00)      ...    ncosp = 11089 (i = 144, j = 78)  ( 178.75, 65.00)
    ! ncosp = 11090 (i =   1, j = 79)  (-178.75, 67.00)      ...    ncosp = 11233 (i = 144, j = 79)  ( 178.75, 67.00)
    ! ncosp = 11234 (i =   1, j = 80)  (-178.75, 69.00)      ...    ncosp = 11377 (i = 144, j = 80)  ( 178.75, 69.00)
    ! ncosp = 11378 (i =   1, j = 81)  (-178.75, 71.00)      ...    ncosp = 11521 (i = 144, j = 81)  ( 178.75, 71.00)
    ! ncosp = 11522 (i =   1, j = 82)  (-178.75, 73.00)      ...    ncosp = 11665 (i = 144, j = 82)  ( 178.75, 73.00)
    ! ncosp = 11666 (i =   1, j = 83)  (-178.75, 75.00)      ...    ncosp = 11809 (i = 144, j = 83)  ( 178.75, 75.00)
    ! ncosp = 11810 (i =   1, j = 84)  (-178.75, 77.00)      ...    ncosp = 11953 (i = 144, j = 84)  ( 178.75, 77.00)
    ! ncosp = 11954 (i =   1, j = 85)  (-178.75, 79.00)      ...    ncosp = 12097 (i = 144, j = 85)  ( 178.75, 79.00)
    ! ncosp = 12098 (i =   1, j = 86)  (-178.75, 81.00)      ...    ncosp = 12241 (i = 144, j = 86)  ( 178.75, 81.00)
    ! ncosp = 12242 (i =   1, j = 87)  (-178.75, 83.00)      ...    ncosp = 12385 (i = 144, j = 87)  ( 178.75, 83.00)
    ! ncosp = 12386 (i =   1, j = 88)  (-178.75, 85.00)      ...    ncosp = 12529 (i = 144, j = 88)  ( 178.75, 85.00)
    ! ncosp = 12530 (i =   1, j = 89)  (-178.75, 87.00)      ...    ncosp = 12673 (i = 144, j = 89)  ( 178.75, 87.00)
    ! ncosp = 12674 (i =   1, j = 90)  (-178.75, 90.00)      ...
    !
    ! ToDo:
    !
    ! Language: Fortran 90
    ! -------------------------------------------------------------------------

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! Subroutine Arguments, in order of appearance
    ! =========================================================================
    character(len=*),intent(in) :: vname
    integer,intent(in),dimension(:) :: imaxj
    integer,intent(in),dimension(7) :: geom_2d
    ! Note:
    !   Allocatable arrays with intent(out) are automatically deallocated when
    !   the program leaves the scope where they are defined but not allocated,
    !   that is, the routine calling remap_cosp_2d(), which is save_cosp().
    real*8,intent(out),dimension(:,:),allocatable :: blob

    ! Subroutine Scalers/Variables
    ! =========================================================================
    integer :: m_cosp,i,j
#ifdef COSP_DEBUG
    logical,dimension(:,:),allocatable :: mask
    logical,dimension(:),allocatable :: mask2
#endif
    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! Allocate domain decomposition aware working array for saving COSP
    ! output via SUBDD.
    allocate(blob(geom_2d(5):geom_2d(6),geom_2d(7):geom_2d(8)))
    ! Initialize as missing value
    !blob = cosp_undef
    blob = 0.0

#ifdef COSP_DEBUG
    ! allocate(mask(geom_2d(5):geom_2d(6),geom_2d(7):geom_2d(8)))
    ! allocate(mask2(SIZE(isccp%totalcldarea)))
    ! mask2 = isccp%totalcldarea.ge.0

    write(0,fmt='(/,25x,"START remap_cosp_2d() for ",A)') trim(vname)
    write(0,fmt='(30x,"Size Blob",2I10)') size(blob,1),size(blob,2)
    select case (trim(vname))
      case ("Lclhcalipso")
        write(0,fmt='(30x,"SIZE clhcalipso",2I10)') &
          SIZE(stlidar%cldlayer,1),SIZE(stlidar%cldlayer,2)
      case ("Lclmcalipso")
        write(0,fmt='(30x,"SIZE clmcalipso",2I10)') &
          SIZE(stlidar%cldlayer,1),SIZE(stlidar%cldlayer,2)
      case ("Lcllcalipso")
        write(0,fmt='(30x,"SIZE cllcalipso",2I10)') &
          SIZE(stlidar%cldlayer,1),SIZE(stlidar%cldlayer,2)
      case ("Lcltcalipso")
        write(0,fmt='(30x,"SIZE cltcalipso",2I10)') &
          SIZE(stlidar%cldlayer,1),SIZE(stlidar%cldlayer,2)
      case ("Lpctisccp")
        write(0,fmt='(30x,"SIZE pctisccp",I10)') SIZE(isccp%meanptop)
      case ("Lcltisccp")
        write(0,fmt='(30x,"SIZE cltisccp",I10)') SIZE(isccp%totalcldarea)
      case ("Lcltmodis")
        write(0,fmt='(30x,"SIZE cltmodis",I10)') SIZE(modis%cloud_fraction_total_mean)
      case ("Lclwmodis")
        write(0,fmt='(30x,"SIZE clwmodis",I10)') SIZE(modis%cloud_fraction_water_mean)
      case ("Lclimodis")
        write(0,fmt='(30x,"SIZE climodis",I10)') SIZE(modis%cloud_fraction_ice_mean)
      case ("Lalbisccp")
        write(0,fmt='(30x,"SIZE albisccp",I10)') SIZE(isccp%meanalbedocld)
      case default
        write(0,*) trim(vname)," unknown to remap_cosp_2d"
      end select

    ! write(0,fmt='(/,30x,"isccp%totalcldarea non-missing count",2I10)') COUNT(mask2)
    ! WRITE(0,FMT='(30x,"MINVAL isccp%totalcldarea",E10.2)') MINVAL(isccp%totalcldarea,MASk=mask2)
    ! WRITE(0,FMT='(30x,"MAXVAL isccp%totalcldarea",F8.2)') MAXVAL(isccp%totalcldarea,MASk=mask2)
#endif
    m_cosp = 0
    do j=geom_2d(3),geom_2d(4)
      do i=geom_2d(1),imaxj(j)
        m_cosp = m_cosp + 1
        select case (trim(vname))
          case ("Lalbisccp")
            blob(i,j) = isccp%meanalbedocld(m_cosp)
          case ("Lpctisccp")
            blob(i,j) = isccp%meanptop(m_cosp)
          case ("Ltauisccp")
            blob(i,j) = isccp%meantaucld(m_cosp)
          case ("Lcltisccp")
            blob(i,j) = isccp%totalcldarea(m_cosp)
          case ("Lmeantbisccp")
            blob(i,j) = isccp%meantb(m_cosp)
          case ("Lmeantbclrisccp")
            blob(i,j) = isccp%meantbclr(m_cosp)
          case ("Lclhcalipso")
            blob(i,j) = stlidar%cldlayer(m_cosp,3)
          case ("Lclmcalipso")
            blob(i,j) = stlidar%cldlayer(m_cosp,2)
          case ("Lcllcalipso")
            blob(i,j) = stlidar%cldlayer(m_cosp,1)
          case ("Lcltcalipso")
            blob(i,j) = stlidar%cldlayer(m_cosp,4)
          case ("Lcltcalipsoice")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,4,1)
          case ("Lcltcalipsoliq")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,4,2)
          case ("Lcltcalipsoun")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,4,3)
          case ("Lcllcalipsoice")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,1,1)
          case ("Lclmcalipsoice")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,2,1)
          case ("Lclhcalipsoice")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,3,1)
          case ("Lcllcalipsoliq")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,1,2)
          case ("Lclmcalipsoliq")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,2,2)
          case ("Lclhcalipsoliq")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,3,2)
          case ("Lcllcalipsoun")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,1,3)
          case ("Lclmcalipsoun")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,2,3)
          case ("Lclhcalipsoun")
            blob(i,j) = stlidar%cldlayerphase(m_cosp,3,3)
          case ("Lcltlidarradar")
            blob(i,j) = stradar%radar_lidar_tcc(m_cosp)
          case ("Lcltmodis")
            blob(i,j) = modis%cloud_fraction_total_mean(m_cosp)
          case ("Lclwmodis")
            blob(i,j) = modis%cloud_fraction_water_mean(m_cosp)
          case ("Lclimodis")
            blob(i,j) = modis%cloud_fraction_ice_mean(m_cosp)
          case ("Lclhmodis")
            blob(i,j) = modis%cloud_fraction_high_mean(m_cosp)
          case ("Lclmmodis")
            blob(i,j) = modis%cloud_fraction_mid_mean(m_cosp)
          case ("Lcllmodis")
            blob(i,j) = modis%cloud_fraction_low_mean(m_cosp)
          case ("Ltautmodis")
           blob(i,j) = modis%optical_thickness_total_mean(m_cosp)
          case ("Ltauwmodis")
           blob(i,j) = modis%optical_thickness_water_mean(m_cosp)
          case ("Ltauimodis")
           blob(i,j) = modis%optical_thickness_ice_mean(m_cosp)
          case ("Ltautlogmodis")
            blob(i,j)=modis%optical_thickness_total_logmean(m_cosp)
          case ("Ltauwlogmodis")
            blob(i,j)=modis%optical_thickness_water_logmean(m_cosp)
          case ("Ltauilogmodis")
            blob(i,j) = modis%optical_thickness_ice_logmean(m_cosp)
          case ("Lreffclwmodis")
            blob(i,j)=modis%cloud_particle_size_water_mean(m_cosp)
          case ("Lreffclimodis")
            blob(i,j)=modis%cloud_particle_size_ice_mean(m_cosp)
          case ("Lpctmodis")
            blob(i,j)=modis%cloud_top_pressure_total_mean(m_cosp)
          case ("Llwpmodis")
            blob(i,j)=modis%liquid_water_path_mean(m_cosp)
          case ("Liwpmodis")
            blob(i,j)=modis%ice_water_path_mean(m_cosp)
          case default
            write(0,*) trim(vname)," unknown to remap_cosp_2d"
            call stop_model("unknown case value in remap_cosp_2d",255)
        end select
      end do
    end do

#ifdef COSP_DEBUG
    ! write(0,fmt='(/,30x,"MINVAL blob",E10.2)') minval(blob)
    ! write(0,fmt='(30x,"MAXVAL blob",F8.2)') maxval(blob)
    ! mask = blob.ge.0
    ! write(0,fmt='(/,30x,"blob non-missing count",2I10)') COUNT(mask)
    ! write(0,fmt='(30x,"minval blob",F8.2)') MINVAL(blob,mask=mask)
    ! write(0,fmt='(30x,"maxval blob",F8.2)') MAXVAL(blob,mask=mask)
    ! deallocate(mask)
    ! deallocate(mask2)
    WRITE(0,FMT='(25x,"END remap_cosp_2d()")')
#endif
    return
end subroutine remap_cosp_2d

  ! ***************************************************************************
  ! PUBLIC SUBROUTINE: remap_cosp_3d
  ! *********************************
  subroutine remap_cosp_3d(vname,imaxj,geom_2d,dim1,dim2,blob)
!@sum  Maps native 3D+ COSP variables so suitable for the SUBDD framework.
!@auth Mike Bauer
    ! -------------------------------------------------------------------------
    ! Inputs:
    !   vname   : Name of COSP variable to be processed
    !   imaxj   : Varying number of used longitudes (See GEOM_B.f)
    !   geom_2d : Domain decomposition parameters (See Below)
    !             (I_0, I_1, J_0, J_1, I_0H, I_1H, J_0H, J_1H)
    !              1    2    3    4    5     6     7     8
    !   dim1    : Length of extra-dimension 1
    !   dim2    : Length of extra-dimension 2
    !   blob    : Allocatable ModelE data array (3D)
    !
    ! Outputs:
    !   blob  : 3D ModelE data array (3D) (I_0H:I_1H, J_0H:J_1H), dim1*dim2)
    !
    ! The same as remap_cosp_2d, except that it deals with higher order COSP
    ! data; compressing the extra dimensions so that they can be stored with
    ! the 3D API of SUBDD. As a result the stored netcdf files for these
    ! variables will require an extra post-processing decompression.
    !
    ! Note that dim1 is generally the number of levels (LM) being stored and
    ! dim2 is a number of fixed collection bins.
    !
    ! Note in cases where the COSP data is 3D, dim1 is set to 1.
    !
    ! Due to limits with the SUBDD sub-daily diagnostics routine some COSP
    ! output variables are not stored in its native shape. This is because
    ! ModelE’s SUBDD routine only allows for 2 and 3D arrays, whereas some COSP
    ! arrays are 4 or 5D.
    !
    ! For these variables we simply collapse the extra dimensions down to 3D,
    ! which then requires post-processing with tool provided by COSP_GISS_TOOLS.
    !
    ! In addition, netcdf file size limitations generally necessitates that COSP
    ! data is stored on a daily basis, rather than monthly. This of course
    ! depends on the model resolution and the time frequency at which COSP is
    ! called.
    !
    ! The primary COSP output variables that have these problems are
    ! LparasolRefl, LcfadLidarsr532 and especially Lcfaddbze94. If these
    ! variables are not requested SUBDD can be safely directed to store
    ! monthly output (see write_daily_files in the rundeck).
    !
    ! ToDo:
    !
    ! Language: Fortran 90
    ! -------------------------------------------------------------------------

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! Subroutine Arguments, in order of appearance
    ! =========================================================================
    character(len=*),intent(in) :: vname
    integer,intent(in),dimension(:) :: imaxj
    integer,intent(in),dimension(7) :: geom_2d
    integer,intent(in) :: dim1,dim2
    real*8,intent(out),dimension(:,:,:),allocatable :: blob
#ifdef COSP_DEBUG
    logical,dimension(:,:,:),allocatable :: mask
    !logical,dimension(:,:,:),allocatable :: mask2
#endif
    ! Subroutine Scalers/Variables
    ! =========================================================================
    integer :: m_cosp,i,j,n_cosp,k,l
    integer :: itau, ipc
    real*8 :: icell
    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! Allocate domain decomposition aware working array for saving COSP
    ! output via SUBDD.
    allocate(blob(geom_2d(5):geom_2d(6),geom_2d(7):geom_2d(8),dim1*dim2))

    ! Initialize as missing value
    !blob = cosp_undef
    blob = 0.0

    ! Test isccp%fq_isccp
    !
    ! write(0,fmt='(30x,"Size clisccp",I10,I10,I10)') &
    !   size(isccp%fq_isccp,1),size(isccp%fq_isccp,2),size(isccp%fq_isccp,3)

    ! Set CF values 1 to 49 by their fq_isccp cell.
    !                (0,-90)->(360,90), thin->thick, low->high
    ! isccp%fq_isccp(npoint,            tau_axid,    pressure2_axid)
    !
    ! Test All grids same PCTau histogram
    !
    ! PC | 1   2   3   4   5   6   7 Tau
    ! -------------------------------
    ! 1  | 1   8   15  22  29  36  43
    ! 2  | 2   9   16  23  30  37  44
    ! 3  | 3   10  17  24  31  38  45
    ! 4  | 4   11  18  25  32  39  46
    ! 5  | 5   12  19  26  33  40  47
    ! 6  | 6   13  20  27  34  41  48
    ! 7  | 7   14  21  28  35  42  49
    ! icell = 0.0
    ! do itau=1,7
    !   ! Increment Tau
    !   do ipc=1,7
    !     ! Increment PC
    !     icell = icell + 1.0
    !     isccp%fq_isccp(:, itau, ipc) = icell
    !   enddo
    ! enddo

#ifdef COSP_DEBUG
    ! allocate(mask(geom_2d(5):geom_2d(6),geom_2d(7):geom_2d(8),dim1*dim2))
    ! !allocate(mask2((size(isccp%fq_isccp,dim=1),dim1,dim2))
    ! !mask2 = isccp%fq_isccp.ge.0
    write(0,fmt='(/,25x,"start remap_cosp_3d()")')

    ! !write(0,fmt='(/,30x,"isccp%fq_isccp non-missing count",2i10)') count(mask2)
    ! !write(0,fmt='(/,30x,"MINVAL isccp%fq_isccp",f8.2)') minval(isccp%fq_isccp,mask=mask2)
    ! !write(0,fmt='(30x,"MAXVAL isccp%fq_isccp",f8.2)') maxval(isccp%fq_isccp,mask=mask2)
    ! write(0,fmt='(/,30x,"MINVAL isccp%fq_isccp",E10.2)') minval(isccp%fq_isccp)
    ! write(0,fmt='(30x,"MAXVAL isccp%fq_isccp",f8.2)') maxval(isccp%fq_isccp)
    select case (trim(vname))
      case ("Lcfaddbze94")
        write(0,fmt='(30x,"Size cfaddbze94"3I10)') &
          size(stradar%cfad_ze,1),size(stradar%cfad_ze,2),size(stradar%cfad_ze,3)
      case ("Lcfadlidarsr532")
         write(0,fmt='(30x,"Size cfadlidarsr532",3I10)') &
          size(stlidar%cfad_sr,1),size(stlidar%cfad_sr,2),size(stlidar%cfad_sr,3)
      case ("Lclcalipso2")
         write(0,fmt='(30x,"Size clcalipso2",2I10)') &
          size(stradar%lidar_only_freq_cloud,1),size(stradar%lidar_only_freq_cloud,2)
      case ("Lclcalipso")
         write(0,fmt='(30x,"Size clcalipso",2I10)') &
          size(stlidar%lidarcld,1),size(stlidar%lidarcld,2)
      case ("Lclisccp")
         write(0,fmt='(30x,"Size clisccp",3I10)') &
          size(isccp%fq_isccp,1),size(isccp%fq_isccp,2),size(isccp%fq_isccp,3)
      case ("Lclcalipsoice")
         write(0,fmt='(30x,"Size clcalipsoice",2I10)') &
          size(stlidar%lidarcldphase,1),size(stlidar%lidarcldphase,2)
      case ("Lclcalipsoliq")
         write(0,fmt='(30x,"Size clcalipsoliq",2I10)') &
          size(stlidar%lidarcldphase,1),size(stlidar%lidarcldphase,2)
      case ("Lparasolrefl")
         write(0,fmt='(30x,"Size parasolrefl",2I10)') &
          size(stlidar%parasolrefl,1),size(stlidar%parasolrefl,2)
      case ("Lclmisr")
         write(0,fmt='(30x,"Size clmisr",3I10)') &
          size(misr%fq_MISR,1),size(misr%fq_MISR,2),size(misr%fq_MISR,3)
      case ("Lcrimodis")
        write(0,fmt='(30x,"SIZE crimodis",3I10)') &
          size(modis%Optical_Thickness_vs_ReffIce,3)
      case ("Lcrlmodis")
        write(0,fmt='(30x,"SIZE crlmodis",3I10)') &
          size(modis%Optical_Thickness_vs_ReffLiq,3)
      case ("Lclmodis")
         write(0,fmt='(30x,"Size clmodis",3I10)') &
          size(modis%Optical_Thickness_vs_Cloud_Top_Pressure,1),size(modis%Optical_Thickness_vs_Cloud_Top_Pressure,2),size(modis%Optical_Thickness_vs_Cloud_Top_Pressure,3)          
    end select
#endif
    ! write(0,fmt='(/,25x,"START remap_cosp_3d() for ",A)') trim(vname)
    ! write(0,fmt='(30x,"Size Blob",I10,I10,I10)') size(blob,1),size(blob,2),size(blob,3)
    ! write(0,fmt='(30x,"Size clisccp",I10,I10,I10)') size(isccp%fq_isccp,1),size(isccp%fq_isccp,2),size(isccp%fq_isccp,3)

    m_cosp=0
    do j=geom_2d(3),geom_2d(4)  ! Loop over Latitude  (J_0 -> J_1)
      do i=geom_2d(1),imaxj(j)  ! Loop over Longitude (I_0 -> imaxj(j))
        m_cosp = m_cosp + 1     ! Dimension Folding (lon, lat) => m_cosp
        n_cosp = 0
        do k=1,dim1             ! Loop over dim1
          do l=1,dim2           ! Loop over dim2
            n_cosp = n_cosp + 1 ! Dimension Folding (dim1, dim2) => n_cosp
            select case (trim(vname))
              case ("Lclisccp")
                !write(0, *) "Lclisccp","blob(i,j,n_cosp)",i,j,n_cosp,"isccp%fq_isccp(m_cosp,k,l)",m_cosp,k,l
                blob(i,j,n_cosp)=isccp%fq_isccp(m_cosp,k,l)
              case ("Lcfaddbze94")
                blob(i,j,n_cosp)=stradar%cfad_ze(m_cosp,k,l)
              case ("Lclmisr")
                blob(i,j,n_cosp)=misr%fq_misr(m_cosp,k,l)
              case ("Lcfadlidarsr532")
                blob(i,j,n_cosp)=stlidar%cfad_sr(m_cosp,k,l)
              case ("Lclcalipso")
                blob(i,j,n_cosp)=stlidar%lidarcld(m_cosp,l)
              case ("Lparasolrefl")
                blob(i,j,n_cosp)=stlidar%parasolrefl(m_cosp,l)
              case ("Lclcalipso2")
                blob(i,j,n_cosp)=stradar%lidar_only_freq_cloud(m_cosp,l)
              case ("Lclmodis")
                blob(i,j,n_cosp)=                                              &
                  modis%Optical_Thickness_vs_Cloud_Top_Pressure(m_cosp,k,l)
              case ("Lclcalipsoice")
                blob(i,j,n_cosp)=stlidar%lidarcldphase(m_cosp,l,1)
              case ("Lclcalipsoliq")
                blob(i,j,n_cosp)=stlidar%lidarcldphase(m_cosp,l,2)
              case ("Lclcalipsoun")
                blob(i,j,n_cosp)=stlidar%lidarcldphase(m_cosp,l,3)
              case ("Lclcalipsotmp")
                blob(i,j,n_cosp)=stlidar%lidarcldtmp(m_cosp,l,1)
              case ("Lclcalipsotmpice")
                blob(i,j,n_cosp)=stlidar%lidarcldtmp(m_cosp,l,2)
              case ("Lclcalipsotmpliq")
                blob(i,j,n_cosp)=stlidar%lidarcldtmp(m_cosp,l,3)
              case ("Lclcalipsotmpun")
                blob(i,j,n_cosp)=stlidar%lidarcldtmp(m_cosp,l,4)
              case ("Lcrimodis")
                blob(i,j,n_cosp)=modis%Optical_Thickness_vs_ReffIce(m_cosp,k,l)
              case ("Lcrlmodis")
                blob(i,j,n_cosp)=modis%Optical_Thickness_vs_ReffLiq(m_cosp,k,l)
              case default
                write(0,*) trim(vname)," unknown to remap_cosp_3d"
                call stop_model("unknown case value in remap_cosp_3d",255)
            end select
          end do
        end do
      end do
    end do

#ifdef COSP_DEBUG
    ! write(0,fmt='(30x,"minval blob",E10.2)') MINVAL(blob)
    ! write(0,fmt='(30x,"maxval blob",F8.2)') MAXVAL(blob)
    ! mask = blob.ge.0
    ! write(0,fmt='(/,30x,"blob non-missing count",3I10)') COUNT(mask)
    ! write(0,fmt='(30x,"minval blob",F8.2)') MINVAL(blob,mask=mask)
    ! write(0,fmt='(30x,"maxval blob",F8.2)') MAXVAL(blob,mask=mask)
    ! deallocate(mask)
    ! !deallocate(mask2)
    write(0,fmt='(25x,"END remap_cosp_3d()")')
#endif
end subroutine remap_cosp_3d

#ifdef COSP_DEBUG
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !---------------------------- PRINT SUBROUTINES -----------------------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! ***************************************************************************
  ! PRIVATE SUBROUTINE: cosp_print_pctau
  ! ************************************
  subroutine cosp_print_pctau(in_pctau,in_cospn,in_im,in_jm)
!@sum  Write ISCCP PC-Tau histograms to a text file.
!@+    Use in save_cosp().
!@auth Mike Bauer


    !  Write individual ISCCP PC-Tau histogrrams to a text file.

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! Subroutine Arguments, in order of appearance
    ! =========================================================================
    real*8,dimension(:,:),intent(in) :: in_pctau
    integer,intent(in) :: in_cospn,in_im,in_jm

    ! Subroutine Parameters/Constants
    ! =========================================================================
    integer,parameter :: ntau = 7
    integer,parameter :: npres = 7
    real*8,parameter :: undef = -1.d30

    ! Subroutine Scalers/Variables
    ! =========================================================================
    integer :: ii,jj,fakei,fakej,ncosp
    real*8 :: total_cf

    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ncosp = (in_im*in_jm) - (2*in_im) + 2
    if (in_cospn .eq. 1) then
        ! Polar Cap (SH)
        fakej = 1
        fakei = 1
    else if (in_cospn .eq. ncosp) then
        ! Polar Cap (NH)
        fakej = in_jm
        fakei = 1
    else
        fakej = floor(real(in_cospn-2)/real(in_im)) + 2
        fakei = modulo(in_cospn-2,in_im) + 1
    endif

    !fakej = floor(real(in_cospn-1)/real(144))+1
    !fakei = modulo((in_cospn-1),144)+1
    total_cf = sum(in_pctau)

    write(0,fmt='(8x,"ncosp = ",i5," (i = ",i3,", j = ",i2,")",2x,"CF = ",f8.4)') &
      in_cospn,fakei,fakej,total_cf
    write(0,fmt='(a62)') '------------------------------------------------------'
    ! do ii=npres,1,-1
    !   write(0,fmt='(f7.1," | ",e10.1,6f7.2)') isccp_pc(ii)*0.01,in_pctau(ii,1),&
    !     (in_pctau(ii,jj), jj=2,ntau)
    ! enddo
    do ii=npres,1,-1
      write(0,fmt='(f7.1," | ",e10.1,6f7.2)') isccp_pc(ii)*0.01,in_pctau(1,ii),&
        (in_pctau(jj,ii), jj=2,ntau)
    enddo
    write(0,fmt='(a62)') '------------------------------------------------------'
    write(0,fmt='(10x,f10.2,6f7.2)') isccp_tau(1), (isccp_tau(jj), jj=2,ntau)
    return

  end subroutine cosp_print_pctau

  ! ***************************************************************************
  ! PRIVATE SUBROUTINE: cosp_input_print
  ! ************************************
  subroutine cosp_input_print(x)
!@sum  Write COSP namelist for inputs (i.e., cosp_inputs) to text file.
!@+    Use in init_cosp().
!@auth Mike Bauer
    ! -------------------------------------------------------------------------
    ! Inputs:
    !   x: COSP namelist for inputs type(cosp_nmlst), i.e., cosp_inputs
    !
    ! Outputs:
    !   Writes to unit=0
    !
    ! ToDo:
    !
    ! Language: Fortran 90
    ! -------------------------------------------------------------------------

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! Subroutine Arguments, in order of appearance
    ! =========================================================================
    type(cosp_nmlst),intent(in) :: x

    ! Subroutine Scalers/Variables
    ! =========================================================================
    integer :: i
    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    write(0,fmt='(/,30x,"%%%%-- Information on COSP Input (cosp_inputs)")')
    write(0,fmt='(30x,"CMOR_NL                   : ",A)') trim(x%cmor_nl)
    write(0,fmt='(30x,"NPOINTS                   = ",I6)') x%npoints
    write(0,fmt='(30x,"NPOINTS_IT                = ",I6)') x%npoints_it
    write(0,fmt='(30x,"NCOLUMNS                  = ",I5)') x%ncolumns
    write(0,fmt='(30x,"NLEVELS                   = ",I3)') x%nlevels
    write(0,fmt='(30x,"USE_VGRID                 = ",L2)') x%use_vgrid
    write(0,fmt='(30x,"NLR                       = ",I3)') x%nlr
    write(0,fmt='(30x,"CSAT_VGRID                = ",L2)') x%csat_vgrid
    write(0,fmt='(30x,"DINPUT                    : ",A)') trim(x%dinput)
    write(0,fmt='(30x,"FINPUT                    : ",A)') trim(x%finput)
    write(0,fmt='(30x,"RADAR_FREQ                = ",F6.2)') x%radar_freq
    write(0,fmt='(30x,"SURFACE_RADAR             = ",I2)') x%surface_radar
    write(0,fmt='(30x,"USE_MIE_TABLES            = ",I2)') x%use_mie_tables
    write(0,fmt='(30x,"USE_GAS_ABS               = ",I2)') x%use_gas_abs
    write(0,fmt='(30x,"DO_RAY                    = ",I2)') x%do_ray
    write(0,fmt='(30x,"MELT_LAY                  = ",I2)') x%melt_lay
    write(0,fmt='(30x,"K2                        = ",F6.2)') x%k2
    write(0,fmt='(30x,"USE_REFF                  = ",L2)') x%use_reff
    write(0,fmt='(30x,"USE_PRECIPITATION_FLUXES  = ",L2)') x%use_precipitation_fluxes
    write(0,fmt='(30x,"NPRMTS_MAX_HYDRO          = ",I3)') x%nprmts_max_hydro
    write(0,fmt='(30x,"NAERO                     = ",I3)') x%naero
    write(0,fmt='(30x,"NPRMTS_MAX_AERO           = ",I3)') x%nprmts_max_aero
    write(0,fmt='(30x,"LIDAR_ICE_TYPE            = ",I3)') x%lidar_ice_type
    write(0,fmt='(30x,"OVERLAP                   = ",I3)') x%overlap
    write(0,fmt='(30x,"ISCCP_TOPHEIGHT           = ",I3)') x%isccp_topheight
    write(0,fmt='(30x,"ISCCP_TOPHEIGHT_DIRECTION = ",I3)') x%isccp_topheight_direction
    write(0,fmt='(30x,"PLATFORM                  = ",I3)') x%platform
    write(0,fmt='(30x,"SATELLITE                 = ",I3)') x%satellite
    write(0,fmt='(30x,"INSTRUMENT                = ",I3)') x%instrument
    write(0,fmt='(30x,"NCHANNELS                 = ",I3)') x%nchannels
    write(0,fmt='(30x,"CHANNELS                 = ",50I3)') (x%channels(i), i=1,rttov_max_channels)
    write(0,fmt='(30x,"SURFEM                   = ",50F6.2)') (x%surfem(i), i=1,rttov_max_channels)
    write(0,fmt='(30x,"ZENANG                    = ",F6.2)') x%zenang
    write(0,fmt='(30x,"CO2                       = ",E10.3E3)') x%co2
    write(0,fmt='(30x,"CH4                       = ",E10.3E3)') x%ch4
    write(0,fmt='(30x,"N2O                       = ",E10.3E3)') x%n2o
    write(0,fmt='(30x,"CO                        = ",E10.3E3)') x%co
    write(0,fmt='(30x,"%%%%-- END (cosp_inputs)")')
    return
  end subroutine cosp_input_print

  ! ***************************************************************************
  ! PRIVATE SUBROUTINE: cosp_output_print
  ! *************************************
  subroutine cosp_output_print(x)
!@sum  Write COSP namelist for outputs (i.e., cfg) to text file.
!@+    Use in init_cosp().
!@auth Mike Bauer
    ! -------------------------------------------------------------------------
    ! Inputs:
    !   x: COSP namelist for outputs type(cosp_config), i.e., cfg
    !
    ! Outputs:
    !   Writes to unit=0
    !
    ! ToDo:
    !
    ! Language: Fortran 90
    ! -------------------------------------------------------------------------

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! Subroutine Arguments, in order of appearance
    ! =========================================================================
    type(cosp_config),intent(in) :: x

    ! Subroutine Scalers/Variables
    ! =========================================================================
    integer :: i
    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    write(0,fmt='(/,30x,"%%%%-- Information on COSP Output (cosp_config)")')
    write(0,fmt='(30x,"Lradar_sim       = ",L2)') x%lradar_sim
    write(0,fmt='(30x,"Llidar_sim       = ",L2)') x%llidar_sim
    write(0,fmt='(30x,"Lisccp_sim       = ",L2)') x%lisccp_sim
    write(0,fmt='(30x,"Lmodis_sim       = ",L2)') x%lmodis_sim
    write(0,fmt='(30x,"Lmisr_sim        = ",L2)') x%lmisr_sim
    write(0,fmt='(30x,"Lrttov_sim       = ",L2)') x%lrttov_sim
    write(0,fmt='(30x,"Ltoffset         = ",L2)') x%ltoffset
    write(0,fmt='(30x,"LcfadDbze94      = ",L2)') x%lcfaddbze94
    write(0,fmt='(30x,"LcfadLidarsr532  = ",L2)') x%lcfadlidarsr532
    write(0,fmt='(30x,"Lclcalipso       = ",L2)') x%lclcalipso
    write(0,fmt='(30x,"Lclhcalipso      = ",L2)') x%lclhcalipso
    write(0,fmt='(30x,"Lcllcalipso      = ",L2)') x%lcllcalipso
    write(0,fmt='(30x,"Lclmcalipso      = ",L2)') x%lclmcalipso
    write(0,fmt='(30x,"Lcltcalipso      = ",L2)') x%lcltcalipso
    write(0,fmt='(30x,"Lclcalipsoliq    = ",L2)') x%lclcalipsoliq
    write(0,fmt='(30x,"Lclcalipsoice    = ",L2)') x%lclcalipsoice
    write(0,fmt='(30x,"Lclcalipsoun     = ",L2)') x%lclcalipsoun
    write(0,fmt='(30x,"Lclcalipsotmp    = ",L2)') x%lclcalipsotmp
    write(0,fmt='(30x,"Lclcalipsotmpun  = ",L2)') x%lclcalipsotmpun
    write(0,fmt='(30x,"Lclcalipsotmpliq = ",L2)') x%lclcalipsotmpliq
    write(0,fmt='(30x,"Lclcalipsotmpice = ",L2)') x%lclcalipsotmpice
    write(0,fmt='(30x,"Lclhcalipsoliq   = ",L2)') x%lclhcalipsoliq
    write(0,fmt='(30x,"Lcllcalipsoliq   = ",L2)') x%lcllcalipsoliq
    write(0,fmt='(30x,"Lclmcalipsoliq   = ",L2)') x%lclmcalipsoliq
    write(0,fmt='(30x,"Lcltcalipsoliq   = ",L2)') x%lcltcalipsoliq
    write(0,fmt='(30x,"Lclhcalipsoice   = ",L2)') x%lclhcalipsoice
    write(0,fmt='(30x,"Lcllcalipsoice   = ",L2)') x%lcllcalipsoice
    write(0,fmt='(30x,"Lclmcalipsoice   = ",L2)') x%lclmcalipsoice
    write(0,fmt='(30x,"Lcltcalipsoice   = ",L2)') x%lcltcalipsoice
    write(0,fmt='(30x,"Lclhcalipsoun    = ",L2)') x%lclhcalipsoun
    write(0,fmt='(30x,"Lcllcalipsoun    = ",L2)') x%lcllcalipsoun
    write(0,fmt='(30x,"Lclmcalipsoun    = ",L2)') x%lclmcalipsoun
    write(0,fmt='(30x,"Lcltcalipsoun    = ",L2)') x%lcltcalipsoun
    write(0,fmt='(30x,"LparasolRefl     = ",L2)') x%lparasolrefl
    write(0,fmt='(30x,"Lclcalipso2      = ",L2)') x%lclcalipso2
    write(0,fmt='(30x,"Lcltlidarradar   = ",L2)') x%lcltlidarradar
    write(0,fmt='(30x,"Lalbisccp        = ",L2)') x%lalbisccp
    write(0,fmt='(30x,"Lalbisccp        = ",L2)') x%lalbisccp
    write(0,fmt='(30x,"Lpctisccp        = ",L2)') x%lpctisccp
    write(0,fmt='(30x,"Lcltisccp        = ",L2)') x%lcltisccp
    write(0,fmt='(30x,"Lclisccp         = ",L2)') x%lclisccp
    write(0,fmt='(30x,"Ltauisccp        = ",L2)') x%ltauisccp
    write(0,fmt='(30x,"Lmeantbisccp     = ",L2)') x%lmeantbisccp
    write(0,fmt='(30x,"Lmeantbclrisccp  = ",L2)') x%lmeantbclrisccp
    write(0,fmt='(30x,"LclMISR          = ",L2)') x%lclmisr
    write(0,fmt='(30x,"Lcltmodis        = ",L2)') x%lcltmodis
    write(0,fmt='(30x,"Lclwmodis        = ",L2)') x%lclwmodis
    write(0,fmt='(30x,"Lclimodis        = ",L2)') x%lclimodis
    write(0,fmt='(30x,"Lclhmodis        = ",L2)') x%lclhmodis
    write(0,fmt='(30x,"Lclmmodis        = ",L2)') x%lclmmodis
    write(0,fmt='(30x,"Lcllmodis        = ",L2)') x%lcllmodis
    write(0,fmt='(30x,"Ltautmodis       = ",L2)') x%ltautmodis
    write(0,fmt='(30x,"Ltauwmodis       = ",L2)') x%ltauwmodis
    write(0,fmt='(30x,"Ltauimodis       = ",L2)') x%ltauimodis
    write(0,fmt='(30x,"Ltautlogmodis    = ",L2)') x%ltautlogmodis
    write(0,fmt='(30x,"Ltauwlogmodis    = ",L2)') x%ltauwlogmodis
    write(0,fmt='(30x,"Ltauilogmodis    = ",L2)') x%ltauilogmodis
    write(0,fmt='(30x,"Lreffclwmodis    = ",L2)') x%lreffclwmodis
    write(0,fmt='(30x,"Lreffclimodis    = ",L2)') x%lreffclimodis
    write(0,fmt='(30x,"Lpctmodis        = ",L2)') x%lpctmodis
    write(0,fmt='(30x,"Llwpmodis        = ",L2)') x%llwpmodis
    write(0,fmt='(30x,"Liwpmodis        = ",L2)') x%liwpmodis
    write(0,fmt='(30x,"Lclmodis         = ",L2)') x%lclmodis
    write(0,fmt='(30x,"Lfracout         = ",L2)') x%lfracout
    write(0,fmt='(30x,"Lcrimodis        = ",L2)') x%lcrimodis
    write(0,fmt='(30x,"Lcrlmodis        = ",L2)') x%lcrlmodis
    write(0,fmt='(30x,"LlidarBetaMol532 = ",L2)') x%llidarbetamol532
    write(0,fmt='(30x,"Ldbze94          = ",L2)') x%ldbze94
    write(0,fmt='(30x,"Lboxptopisccp    = ",L2)') x%lboxptopisccp
    write(0,fmt='(30x,"Lboxtauisccp     = ",L2)') x%lboxtauisccp
    write(0,fmt='(30x,"Latb532          = ",L2)') x%latb532
    write(0,fmt='(30x,"Ltbrttov         = ",L2)') x%ltbrttov
    write(0,fmt='(30x,"Lstats           = ",L2)') x%lstats
    write(0,fmt='(30x,"Lwrite_output    = ",L2)') x%lwrite_output
    !WRITE(0,*) "cfg%out_list"
    !WRITE(0,"(45A)") (x%out_list(i),i=1,N_OUT_LIST)
    write(0,fmt='(30x,"%%%%-- END (cosp_config)")')
    return
  end subroutine cosp_output_print

  ! ***************************************************************************
  ! PRIVATE SUBROUTINE: cosp_input_check
  ! ************************************
  subroutine cosp_input_check(np_cosp,itime,lm)
!@sum  Write statistics MIN/MEAN/MAX of all modelE data inputs going into
!@+    COSP for a single time-step (i.e., gbx) to text file.
!@+    Use in ModelE just before calling run_cosp_sims().
!@+    Use in ModelE just after filling gbx for a individual model grid.
!@auth Mike Bauer
    !
    ! Prints out column view of cosp input data (min/mean/max/cnt) for
    !   checking validity of input values
    !
    implicit none
    integer,intent(in) :: np_cosp,lm
    integer,intent(in) :: itime
    integer :: cosp_i,cosp_cnter,cosp_cnter1,l,hydo_cnter
    integer, dimension(1:9), parameter :: non_prec = (/ 1,1,0,0,1,1,0,0,0 /)
    real*8 :: cosp_maxval,cosp_minval,cosp_sumval,cosp_cntval,cosp_meanval
    real*8 :: cosp_tmp_value,cosp_tmp2_value,cosp_tmp3_value,                  &
            cosp_tmp4_value,cosp_tmp5_value,cosp_tmp6_value,                   &
            cosp_tmp7_value,cosp_tmp8_value,cosp_tmp9_value,                   &
            cosp_tmp10_value
    logical inconsistent
    character(len=8), parameter :: hydroname(9) = (/"I_LSCLIQ","I_LSCICE",     &
      "I_LSRAIN","I_LSSNOW","I_CVCLIQ","I_CVCICE","I_CVRAIN","I_CVSNOW",       &
      "I_LSGRPL"/)
    character(len=10), parameter :: cldtype(9) = (/"Stratiform","Stratiform",  &
      "Stratiform","Stratiform","Convective","Convective","Convective",        &
      "Convective","Stratiform"/)
    character(len=12), parameter :: hydtype(9) = (/"Cloud Liquid","Cloud Ice", &
      "Rain","Snow","Cloud Liquid","Cloud Ice","Rain","Snow","Graupel"/)
    character(len=20) :: hydro_label
    !---------------- End of declaration of variables --------------
    if (cosp_inputs%use_precipitation_fluxes) then
      hydro_label = "Flux [kg/m2.s,"
    else
      hydro_label = "Mixing Ratio [kg/kg,"
    endif

    201 format("Layer (",I2,")",5x,"Pressure ",F8.3," gbx%reff(",A8,           &
      ") = MIN/MEAN/MAX/CNT ",3e12.3, F8.0)
    202 format("Layer (",I2,")",5x,"Ave. Pressure ",F8.3," gbx%mr_hydro(",A8,  &
      ") = MIN/MEAN/MAX/CNT ",4e12.3)
    203 format("Layer (",I2,")",5x,"Ave. Pressure ",F8.3,"|      MIN = ",F8.2, &
      ",   MEAN = ",F8.2,",   MAX = ",F8.2,",   Valid CNT =",F8.0)
    204 format("Mean Temperature [K, gbx%t]")
    206 format("Mean Relative Humidity [%, to H2O, gbx%q]")
    207 format("Layer (",I2,")",5x,"Ave. Pressure ",F8.3,"|      MIN = ",e12.3,&
      ",   MEAN = ",e12.3,",   MAX = ",e12.3,",   Valid CNT =",F8.0)
    208 format("Mean Specific Humidity [kg/kg, to water, gbx%sh]")
    210 format("Mean Height above Sea-Level [m, terrain following, gbx%zlev]")
    211 format("Layer (",I2,")",5x,"Ave. Pressure ",F8.3,"|      MIN = ",F8.4, &
      ",   MEAN = ",F8.4,",   MAX = ",F8.4,",   Valid CNT =",F8.0)
    212 format("Mean Total Cloud Fraction [-, gbx%tca]")
    214 format("Mean Stratiform Cloud Optical Depth [-, gbx%dtau_s]")
    216 format("Mean Stratiform Cloud Longwave Emissivity [-, gbx%dem_s]")
    218 format("Mean Convective Cloud Optical Depth [-, gbx%dtau_c]")
    220 format("Mean Convective Cloud Longwave Emissivity [-, gbx%dem_c]")

    300 format("Mean ",A10,1x,A12,1x,"Effective Radius [m, gbx%reff(",A8,")]")
    301 format("Mean ",A10,1x,A12,1x,"Mixing Ratio [kg/kg, gbx%mr_hydro(",A8,  &
      ")]")
    302 format("Mean ",A10,1x,A12,1x,A20,1x,"gbx%mr_hydro(",A8,")]")
    303 format("-------------------------------------------------------------")
    304 format(" Time (",I10,")")
    305 format("Layer (",I2,")",5x,"No Qualifying Grids")
    307 format("Mean Pressure [hPa, gbx%p]")

    ! Dump various non=hydrometeor info
    write(0, *)
    write(0,304) itime
    cospvars: do cosp_cnter=0,9
      write(0,303)
      if(cosp_cnter.eq.0) then
        write(0,307)
      else if(cosp_cnter.eq.1) then
        write(0,204)
      else if(cosp_cnter.eq.2) then
        write(0,206)
      else if(cosp_cnter.eq.3) then
        write(0,208)
      else if(cosp_cnter.eq.4) then
        write(0,210)
      else if(cosp_cnter.eq.5) then
        write(0,212)
      else if(cosp_cnter.eq.6) then
        write(0,214)
      else if(cosp_cnter.eq.7) then
        write(0,216)
      else if(cosp_cnter.eq.8) then
        write(0,218)
      else if(cosp_cnter.eq.9) then
        write(0,220)
      endif
      write(0,303)
      ! Loop over Layers
      layers: do l=1,lm
        ! Conditional Average over all non-zero data grids

        ! Conditional non-zero count and sum
        cosp_tmp2_value = 0
        cosp_tmp_value = 0
        ! Conditional non-zero Minimum over all grids
        !cosp_tmp3_value = 1
        cosp_tmp3_value = cosp_teeny
        ! Conditional non-zero Maximum over all grids
        cosp_tmp4_value = -1
        cosp_tmp5_value = 0
        cosp_tmp6_value = 0

        !WRITE(0, *) "Layer",L
        ! Loop over all Grids
        grids: do cosp_i = 1, cosp_inputs%npoints

          ! Screen grids by some criteria (Optional)
          ! ----------------------------------------

          ! ! Screen by Sunlit [1 for day, 0 for night]
          ! IF(gbx%sunlit(cosp_i).EQ.0) CYCLE GRIDS

          ! ! Screen by Total Clear-sky
          ! IF(SUM(gbx%tca(cosp_i,:)).GT.0.000001) CYCLE GRIDS

          ! ! Screen by Layer Clear-sky
          ! IF(gbx%tca(cosp_i,LM+1-L).GT.0.000001) CYCLE GRIDS

          ! ! Screen by non-zero Layer Cloud
          ! IF(gbx%tca(cosp_i,LM+1-L).LE.0.000001) CYCLE GRIDS

          ! ! Screen by non-zero Total Column Cloud
          ! IF(SUM(gbx%tca(cosp_i,:)).LE.0.000001) CYCLE GRIDS
          ! IF(SUM(gbx%tca(cosp_i,:)).LE.0.01) CYCLE GRIDS

          ! ! Screen by inconsistency
          ! inconsistent=.false.
          ! DO cosp_cnter=1,8
          !     ! ! Non-Zero hydro but zero radius
          !     ! IF ((gbx%mr_hydro(cosp_i,LM+1-L,cosp_cnter)>0.0).and.( &
          !     !     gbx%Reff(cosp_i,LM+1-L,cosp_cnter)<=0.0)) THEN
          !     !         cosp_tmp10_value = 1
          !     !         inconsistent=.true.
          !     !         ENDIF
          !     ! ! Zero hydro but non-zero radius
          !     ! IF ((gbx%mr_hydro(cosp_i,LM+1-L,cosp_cnter)<=0.0).and.( &
          !     !     gbx%Reff(cosp_i,LM+1-L,cosp_cnter)>0.0)) THEN
          !     !         inconsistent=.true.
          !     ! ENDIF
          !     ! Non-Zero hydro and non-zero radius
          !     IF ((gbx%mr_hydro(cosp_i,LM+1-L,cosp_cnter).GT.0.0).AND.(        &
          !       gbx%Reff(cosp_i,LM+1-L,cosp_cnter).GT.0.0)) THEN
          !         inconsistent=.true.
          !     ENDIF
          !     ! ! Zero hydro and zero radius
          !     ! IF ((gbx%mr_hydro(cosp_i,LM+1-L,cosp_cnter).GE.0.0).and.( &
          !     !     gbx%Reff(cosp_i,LM+1-L,cosp_cnter).GE.0.0)) THEN
          !     !         inconsistent=.true.
          !     ! ENDIF
          ! ENDDO
          ! !IF(inconsistent) CYCLE GRIDS
          ! IF(.NOT.inconsistent) CYCLE GRIDS

          ! Conditional on Screen Data Collection
          ! -------------------------------------
          ! Layer Pressure for this Grid
          !cosp_tmp8_value = gbx%p(cosp_i,LM+1-L)*0.01

          ! Conditional on Screen sum of grid layer pressures
          cosp_tmp5_value = cosp_tmp5_value+gbx%p(cosp_i,LM+1-L)*0.01
          ! Conditional on Screen grid count (= np_cosp if no screens)
          cosp_tmp6_value = cosp_tmp6_value+1
          !cosp_tmp6_value = cosp_teeny

          if(cosp_cnter.eq.0) then
            ! pressure at model levels [hpa]
            cosp_tmp9_value = gbx%p(cosp_i,lm+1-l)*0.01
          else if(cosp_cnter.eq.1) then
            ! temperature at model levels [k]
            cosp_tmp9_value = gbx%t(cosp_i,lm+1-l)
          else if(cosp_cnter.eq.2) then
            !  relative humidity at model levels (%, to water)
            cosp_tmp9_value = gbx%q(cosp_i,lm+1-l)
          else if(cosp_cnter.eq.3) then
            ! specific humidity [kg/kg, to water]
            cosp_tmp9_value = gbx%sh(cosp_i,lm+1-l)
          else if(cosp_cnter.eq.4) then
            ! height of model levels [m]
            cosp_tmp9_value = gbx%zlev(cosp_i,lm+1-l)
          else if(cosp_cnter.eq.5) then
            ! total cloud fraction at model levels [unitless]
            cosp_tmp9_value = gbx%tca(cosp_i,lm+1-l)
          else if(cosp_cnter.eq.6) then
            ! mean 0.67 micron optical depth of stratiform clouds at model levels
            cosp_tmp9_value = gbx%dtau_s(cosp_i,lm+1-l)
          else if(cosp_cnter.eq.7) then
            ! 10.5 micron longwave emissivity of stratiform clouds at model levels
            cosp_tmp9_value = gbx%dem_s(cosp_i,lm+1-l)
          else if(cosp_cnter.eq.8) then
            ! mean 0.67 micron optical depth of convective clouds at model levels
            cosp_tmp9_value = gbx%dtau_c(cosp_i,lm+1-l)
          else if(cosp_cnter.eq.9) then
            ! 10.5 micron longwave emissivity of convective clouds at model levels [unitless]
            cosp_tmp9_value = gbx%dem_c(cosp_i,lm+1-l)
          endif

          ! conditional non-zero grid-value collection
          ! ------------------------------------------
          !if(cosp_tmp9_value.ge.0.0) then
          if(cosp_tmp9_value.gt.0.0) then
            ! accumulate sum and cnt
            cosp_tmp_value=cosp_tmp_value+cosp_tmp9_value
            cosp_tmp2_value=cosp_tmp2_value+1
            if(cosp_tmp9_value > cosp_tmp4_value) then
              cosp_tmp4_value = cosp_tmp9_value
            endif
            if(cosp_tmp9_value < cosp_tmp3_value) then
              cosp_tmp3_value = cosp_tmp9_value
            endif
          endif
        end do grids

        ! Find Conditional non-zero Average
        if(cosp_tmp2_value.gt.0.0) then
          cosp_tmp8_value = cosp_tmp_value/cosp_tmp2_value
        else
          ! no valid grid values
          cosp_tmp8_value =  -999
          cosp_tmp4_value =  -999
          cosp_tmp3_value =  -999
          write(0,305) lm+1-l
          cycle layers
        endif

        ! Conditional on Screen Average layer pressure
        cosp_tmp7_value = cosp_tmp5_value/cosp_tmp6_value

        if(cosp_cnter.eq.0) then
          write(0,203) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        else if(cosp_cnter.eq.1) then
          write(0,203) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        else if(cosp_cnter.eq.2) then
          write(0,203) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        else if(cosp_cnter.eq.3) then
          write(0,207) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        else if(cosp_cnter.eq.4) then
          write(0,203) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        else if(cosp_cnter.eq.5) then
          write(0,211) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        else if(cosp_cnter.eq.6) then
          write(0,207) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        else if(cosp_cnter.eq.7) then
          write(0,207) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        else if(cosp_cnter.eq.8) then
          write(0,207) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        else if(cosp_cnter.eq.9) then
          write(0,207) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,                 &
            cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
        endif

      enddo layers
    enddo cospvars

    ! print hydrometero data
    do hydo_cnter=1,2
      do cosp_cnter=1,8
        write(0,303)
        if(hydo_cnter.eq.1) then
            ! effective radius
            write(0,300) cldtype(cosp_cnter), hydtype(cosp_cnter), hydroname(cosp_cnter)
        else
            ! mixing ratio or flux
            if (non_prec(cosp_cnter).eq.1) then
                ! non-precipitation types
                write(0,301) cldtype(cosp_cnter), hydtype(cosp_cnter), hydroname(cosp_cnter)
            else
                ! precipitation
                write(0,302) cldtype(cosp_cnter), hydtype(cosp_cnter), hydro_label, hydroname(cosp_cnter)
            endif
        endif
        write(0,303)
        ! loop over layers
        layers1: do l=1,lm
          ! conditional average over all non-zero data grids

          ! conditional non-zero count and sum
          cosp_tmp2_value = 0
          cosp_tmp_value = 0

          ! conditional non-zero minimum over all grids
          cosp_tmp3_value = 1e6

          ! conditional non-zero maximum over all grids
          cosp_tmp4_value = -1e6

          cosp_tmp5_value = 0
          cosp_tmp6_value = 0

          ! loop over all grids
          grids1: do cosp_i = 1, cosp_inputs%npoints

            ! screen grids by some criteria (optional)
            ! ----------------------------------------

            ! screen by sunlit [1 for day, 0 for night]
            !if(gbx%sunlit(cosp_i).eq.0) cycle

            ! screen by total column clear-sky
            !if(sum(gbx%tca(cosp_i,:)).gt.0.000001) cycle grids1

            ! screen by layer clear-sky
            !if(gbx%tca(cosp_i,lm+1-l).gt.0.000001) cycle grids1

            ! screen by layer cloud
            !if(gbx%tca(cosp_i,lm+1-l).le.0.000001) cycle grids1

            !screen by non-zero total column cloud
            !if(sum(gbx%tca(cosp_i,:)).le.0.000001) cycle grids1
            !if(sum(gbx%tca(cosp_i,:)).le.0.01) cycle grids1

            ! screen by inconsistency
            ! inconsistent=.false.
            if(hydo_cnter.eq.1) then
              ! effective radius
              cosp_tmp9_value = gbx%reff(cosp_i,lm+1-l,cosp_cnter)
              ! screen by inconsistency
              ! ! non-zero hydro but zero radius
              ! if ((gbx%mr_hydro(cosp_i,lm+1-l,cosp_cnter)>0.0).and.(&
              !     gbx%reff(cosp_i,lm+1-l,cosp_cnter)<=0.0)) then
              !         inconsistent=.true.
              ! endif
              ! ! zero hydro but non-zero radius
              ! if ((gbx%mr_hydro(cosp_i,lm+1-l,cosp_cnter)<=0.0).and.(&
              !     gbx%reff(cosp_i,lm+1-l,cosp_cnter)>0.0)) then
              !         inconsistent=.true.
              ! endif
              ! non-zero hydro and non-zero radius
              ! if ((gbx%mr_hydro(cosp_i,lm+1-l,cosp_cnter)>0.0).and.(&
              !     gbx%reff(cosp_i,lm+1-l,cosp_cnter)>0.0)) then
              !         inconsistent=.true.
              ! endif
              ! !  zero hydro and  zero radius
              ! if ((gbx%mr_hydro(cosp_i,lm+1-l,cosp_cnter).eq.0.0).and.(
              !     gbx%reff(cosp_i,lm+1-l,cosp_cnter).eq.0.0)) then
              !         inconsistent=.true.
              ! endif
            else
              ! mixing ratio or flux
              cosp_tmp9_value = gbx%mr_hydro(cosp_i,lm+1-l,cosp_cnter)
              ! screen by inconsistency
              ! ! non-zero hydro but zero radius
              ! if ((gbx%mr_hydro(cosp_i,lm+1-l,cosp_cnter)>0.0).and.(&
              !     gbx%reff(cosp_i,lm+1-l,cosp_cnter)<=0.0)) then
              !     inconsistent=.true.
              ! endif
              ! ! zero hydro but non-zero radius
              ! if ((gbx%mr_hydro(cosp_i,lm+1-l,cosp_cnter)<=0.0).and.(&
              !     gbx%reff(cosp_i,lm+1-l,cosp_cnter)>0.0)) then
              !     inconsistent=.true.
              ! endif
              ! non-zero hydro and non-zero radius
              ! if ((gbx%mr_hydro(cosp_i,lm+1-l,cosp_cnter)>0.0).and.(&
              !     gbx%reff(cosp_i,lm+1-l,cosp_cnter)>0.0)) then
              !     inconsistent=.true.
              ! endif
              ! ! zero hydro and  zero radius
              ! if ((gbx%mr_hydro(cosp_i,lm+1-l,cosp_cnter).eq.0.0).and.(&
              !     gbx%reff(cosp_i,lm+1-l,cosp_cnter).eq.0.0)) then
              !     inconsistent=.true.
              ! endif
            endif
            ! !if(inconsistent) cycle grids1
            ! if(.not.inconsistent) cycle grids1

            ! conditional on screen data collection
            ! -------------------------------------
            ! layer pressure for this grid
            !cosp_tmp8_value = gbx%p(cosp_i,lm+1-l)*0.01

            ! conditional on screen sum of grid layer pressures
            cosp_tmp5_value = cosp_tmp5_value+gbx%p(cosp_i,lm+1-l)*0.01
            ! conditional on screen grid count (= np_cosp if no screens)
            cosp_tmp6_value = cosp_tmp6_value + 1

            ! conditional non-zero grid-value collection
            ! ------------------------------------------
            !if(cosp_tmp9_value.ge.0.0) then
            if(cosp_tmp9_value.gt.0.0) then
              ! accumulate sum and cnt
              cosp_tmp_value=cosp_tmp_value+cosp_tmp9_value
              cosp_tmp2_value=cosp_tmp2_value+1
              if(cosp_tmp9_value > cosp_tmp4_value) then
                cosp_tmp4_value = cosp_tmp9_value
              endif
              if(cosp_tmp9_value < cosp_tmp3_value) then
                cosp_tmp3_value = cosp_tmp9_value
              endif
            endif
          end do grids1

          ! find conditional non-zero average
          if(cosp_tmp2_value.gt.0.0) then
            cosp_tmp8_value = cosp_tmp_value/cosp_tmp2_value
          else
            ! no valid grid values
            cosp_tmp8_value = -999
            cosp_tmp4_value = -999
            cosp_tmp3_value = -999
            write(0,305) lm+1-l
            cycle layers1
          endif

          ! conditional on screen average layer pressure
          cosp_tmp7_value = cosp_tmp5_value/cosp_tmp6_value

          if(hydo_cnter.eq.1) then
            ! effective radius
            write(0,207) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,             &
              cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
          else
            ! mixing ratio or flux
            write(0,207) lm+1-l,cosp_tmp7_value,cosp_tmp3_value,             &
              cosp_tmp8_value,cosp_tmp4_value,cosp_tmp2_value
          endif
        enddo layers1
      enddo
    enddo
    return
  end subroutine cosp_input_check

  ! ***************************************************************************
  ! PRIVATE SUBROUTINE: print_cosp_gbx
  ! **********************************
  subroutine print_cosp_gbx(np_cosp, pfluxes, reff, one_gbx, itime)
!@sum  Write details of modelE data inputs going into individual COSP
!@+    model grid (i.e., gbx) to text file. Use after gbx filled for each grid.
!@auth Mike Bauer

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! Subroutine Arguments, in order of appearance
    ! =========================================================================
    integer,intent(in) :: np_cosp
    logical,intent(in) :: pfluxes
    logical,intent(in) :: reff
    type(cosp_gridbox), intent(in) :: one_gbx
    integer,intent(in) :: itime

    ! Subroutine Parameters/Constants
    ! =========================================================================
    character (len=8), parameter :: hydroname(9) = (/"I_LSCLIQ","I_LSCICE",    &
      "I_LSRAIN","I_LSSNOW","I_CVCLIQ","I_CVCICE","I_CVRAIN","I_CVSNOW",       &
      "I_LSGRPL"/)
    character (len=10), parameter :: cldtype(9) = (/"Stratiform","Stratiform", &
      "Stratiform","Stratiform","Convective","Convective","Convective",        &
      "Convective","Stratiform"/)
    character (len=12), parameter :: hydtype(9) = (/"Cloud Liquid","Cloud Ice",&
      "Rain","Snow","Cloud Liquid","Cloud Ice","Rain","Snow","Graupel"/)
    character(len=*),parameter :: fmt0 =                                       &
      '("itime",i0,/,"np_cosp: ",i0,/,2x,"Lat  :",f0.2,/,2x,"Lon  :",f0.2,/)'
    ! The shape of an array is a 1-dimensional array of the number of elements in each dimension
    character(len=*),parameter :: fmt1 = '(4x,"SHAPE(gbx%",A16,"): (",i0,")")'
    character(len=*),parameter :: fmt2 = '(4x,"SHAPE(gbx%",A16,"): (",i0,1x,i0,")")'
    character(len=*),parameter :: fmt3 = '(4x,"SHAPE(gbx%",A16,"): (",i0,1x,i0,1x,i0")")'
    character(len=*),parameter :: fmt4 = '(4x,"SHAPE(gbx%",A16,"): (",i0,1x,i0,1x,i0,1x,i0")")'

    character(len=20) :: hydro_label
    integer,dimension(1:9), parameter :: non_prec = (/ 1,1,0,0,1,1,0,0,0 /)

    ! Subroutine Scalers/Variables
    ! =========================================================================
    integer :: ii,jj

    ! Subroutine Arrays
    ! =========================================================================
    integer,dimension(2) :: shape2
    integer,dimension(3) :: shape3
    integer,dimension(4) :: shape4

    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    if (cosp_inputs%use_precipitation_fluxes) then
      hydro_label = "Flux [kg/m2.s,"
    else
      hydro_label = "Mixing Ratio [kg/kg,"
    endif
    ! Example:
    !   --------------------------------------------------------------------------------
    !   p_cosp:        1
    !    Lat  :  181.25
    !    Lon  :  181.25
    !
    !    Dimensions:
    !
    !      SHAPE(gbx%       longitude): (12960)
    !      SHAPE(gbx%        latitude): (12960)
    !      SHAPE(gbx%            land): (12960)
    !      SHAPE(gbx%            psfc): (12960)
    !      SHAPE(gbx%          sunlit): (12960)
    !      SHAPE(gbx%             skt): (12960)
    !      SHAPE(gbx%            zlev): (12960 40)
    !      SHAPE(gbx%       zlev_half): (12960 40)
    !      SHAPE(gbx%            dlev): (12960 40)
    !      SHAPE(gbx%               p): (12960 40)
    !      SHAPE(gbx%              ph): (12960 40)
    !      SHAPE(gbx%               t): (12960 40)
    !      SHAPE(gbx%               q): (12960 40)
    !      SHAPE(gbx%              sh): (12960 40)
    !      SHAPE(gbx%          dtau_s): (12960 40)
    !      SHAPE(gbx%          dtau_c): (12960 40)
    !      SHAPE(gbx%           dem_s): (12960 40)
    !      SHAPE(gbx%           dem_c): (12960 40)
    !      SHAPE(gbx%        mr_ozone): (12960 40)
    !      SHAPE(gbx%             tca): (12960 40)
    !      SHAPE(gbx%             cca): (12960 40)
    !      SHAPE(gbx%         rain_ls): (12960 40)
    !      SHAPE(gbx%         rain_cv): (12960 40)
    !      SHAPE(gbx%         snow_ls): (12960 40)
    !      SHAPE(gbx%         snow_cv): (12960 40)
    !      SHAPE(gbx%         grpl_ls): (12960 40)
    !      SHAPE(gbx%        mr_hydro): (12960 40)
    !      SHAPE(gbx%

    write(0,'(/,A80)') repeat("-", 80)
    write(0,fmt0) itime, np_cosp,one_gbx%latitude(np_cosp),one_gbx%longitude(np_cosp)

    write(0,'(2x,"Dimensions:")')
    ! point information (npoints)
    write(0,fmt1) 'longitude',shape(one_gbx%longitude)
    write(0,fmt1) 'latitude',shape(one_gbx%latitude)
    write(0,fmt1) 'land',shape(one_gbx%land)
    write(0,fmt1) 'psfc',shape(one_gbx%psfc)
    write(0,fmt1) 'sunlit',shape(one_gbx%sunlit)
    write(0,fmt1) 'skt',shape(one_gbx%skt)
    ! gridbox information (npoints,nlevels)
    shape2 = shape(one_gbx%zlev)
    write(0,fmt2) 'zlev',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%zlev_half)
    write(0,fmt2) 'zlev_half',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%dlev)
    write(0,fmt2) 'dlev',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%p)
    write(0,fmt2) 'p',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%ph)
    write(0,fmt2) 'ph',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%t)
    write(0,fmt2) 't',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%q)
    write(0,fmt2) 'q',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%sh)
    write(0,fmt2) 'sh',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%dtau_s)
    write(0,fmt2) 'dtau_s',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%dtau_c)
    write(0,fmt2) 'dtau_c',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%dem_s)
    write(0,fmt2) 'dem_s',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%dem_c)
    write(0,fmt2) 'dem_c',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%mr_ozone)
    write(0,fmt2) 'mr_ozone',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%tca)
    write(0,fmt2) 'tca',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%cca)
    write(0,fmt2) 'cca',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%rain_ls)
    write(0,fmt2) 'rain_ls',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%rain_cv)
    write(0,fmt2) 'rain_cv',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%snow_ls)
    write(0,fmt2) 'snow_ls',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%snow_cv)
    write(0,fmt2) 'snow_cv',(shape2(ii),ii=1,2)
    shape2 = shape(one_gbx%grpl_ls)
    write(0,fmt2) 'grpl_ls',(shape2(ii),ii=1,2)
    ! hydrometeors concentration and distribution parameters (npoints,nlevels,nhydro)
    shape3 = shape(one_gbx%mr_hydro)
    write(0,fmt3) 'mr_hydro',(shape3(ii),ii=1,3)
    shape3 = shape(one_gbx%reff)
    write(0,fmt3) 'reff',(shape3(ii),ii=1,3)
    shape3 = shape(one_gbx%np)
    write(0,fmt3) 'np',(shape3(ii),ii=1,3)
    ! hydrometeors concentration and distribution parameters (nprmts_max_hydro,nhydro)
    shape2 = shape(one_gbx%dist_prmts_hydro)
    write(0,fmt2) 'dist_prmts_hydro',(shape2(ii),ii=1,2)
    ! aerosol concentration and distribution parameters (npoints,nlevels,naero)
    shape3 = shape(one_gbx%conc_aero)
    write(0,fmt3) 'conc_aero',(shape3(ii),ii=1,3)
    ! particle size distribution type for each aerosol species (naero)
    write(0,fmt1) 'dist_type_aero',shape(one_gbx%dist_type_aero)
    ! distributional parameters for aerosols (npoints,nlevels,nprmts_max_aero,naero)
    shape4 = shape(one_gbx%dist_prmts_aero)
    write(0,fmt4) 'dist_prmts_aero',(shape4(ii),ii=1,4)

    ii = np_cosp
    ! point information (npoints)
    write(0,fmt='(/,a)') "Point Data"
    write(0,fmt='(5x,a35,a1,1x,f0.2)') "Landmask [0=ocean, 1=land]",":",one_gbx%land(ii)
    write(0,fmt='(5x,a35,a1,1x,f10.2)') "Surface Pressure [Pa]",":",one_gbx%psfc(ii)
    write(0,fmt='(5x,a35,a1,1x,f0.2)') "Skin Temperature [K]",":",one_gbx%skt(ii)
    write(0,fmt='(5x,a35,a1,1x,f0.2)') "Sunlit [1 For Day, 0 For Night]",":",one_gbx%sunlit(ii)
    ! Gridbox Information (npoints,nlevels)
    write(0,fmt='(/,5x,a)') "Gridbox Data"
    do jj=one_gbx%nlevels,1,-1
      write(0,fmt='(/,10x,a35,a1,1x,i0)') "Level",":",jj
      write(0,fmt='(15x,a35,a1,1x,f0.2)') "Pressure [Pa]",":",one_gbx%p(ii,jj)
      write(0,fmt='(15x,a35,a1,1x,f0.2)') "Pressure at Lower Edge [Pa]",":",one_gbx%ph(ii,jj)
      write(0,fmt='(15x,a35,a1,1x,f0.2)') "Temperature [K]",":",one_gbx%t(ii,jj)
      write(0,fmt='(15x,a35,a1,1x,f0.2)') "Height [m]",":",one_gbx%zlev(ii,jj)
      write(0,fmt='(15x,a35,a1,1x,f0.2)') "Height at Lower Edge [m]",":",one_gbx%zlev_half(ii,jj)
      write(0,fmt='(15x,a35,a1,1x,f0.2)') "Depth of Model Layer [m]",":",one_gbx%dlev(ii,jj)
      write(0,fmt='(15x,a35,a1,1x,f0.2)') "Relative Humidity [%]",":",one_gbx%q(ii,jj)
      write(0,fmt='(15x,a35,a1,1x,e10.3)') "Specific Humidity [kg/kg]",":",one_gbx%sh(ii,jj)
      write(0,fmt='(15x,a35,a1,1x,e10.3)') "Total Cloud Fraction [unitless]",":",one_gbx%tca(ii,jj)
      write(0,fmt='(15x,a)') "Stratiform/Large-Scale Cloud Info"
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "0.67 um Optical Depth [unitless]",":",one_gbx%dtau_s(ii,jj)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "10.5 um Emissivity [unitless]",":",one_gbx%dem_s(ii,jj)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Mixing Ratio of Cloud Liq. [kg/kg]",":",one_gbx%mr_hydro(ii,jj,i_lscliq)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Effective Radius of Cloud Liq. [m]",":",one_gbx%reff(ii,jj,i_lscliq)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Mixing Ratio of Cloud ice [kg/kg]",":",one_gbx%mr_hydro(ii,jj,i_lscice)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Effective Radius of Cloud Ice [m]",":", one_gbx%reff(ii,jj,i_lscice)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Flux of Rain [kg/m2.s]",":",one_gbx%mr_hydro(ii,jj,i_lsrain)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Effective Radius of Rain [m]",":",one_gbx%reff(ii,jj,i_lsrain)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Flux of Snow [kg/m2.s]",":",one_gbx%mr_hydro(ii,jj,i_lssnow)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Effective Radius of Snow [m]",":",one_gbx%reff(ii,jj,i_lssnow)
      !write(0,fmt='(20x,a35,a1,1x,a)') "Flux of Graupel [kg/m2.s]",":","n/a"
      !write(0,fmt='(20x,a35,a1,1x,a)') "Effective Radius of Graupel [m]",":","n/a"
      write(0,fmt='(15x,a)') "Column Based Convective Cloud Info"
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "0.67 um Optical Depth [unitless]",":",one_gbx%dtau_c(ii,jj)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "10.5 um Emissivity [unitless]",":",one_gbx%dem_c(ii,jj)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Mixing Ratio of Cloud Liq. [kg/kg]",":",one_gbx%mr_hydro(ii,jj,i_cvcliq)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Effective Radius of Cloud Liq [m]",":",one_gbx%reff(ii,jj,i_cvcliq)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Mixing Ratio of Cloud Ice [kg/kg]",":",one_gbx%mr_hydro(ii,jj,i_cvcice)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Effective Radius of Cloud Ice [m]",":",one_gbx%reff(ii,jj,i_cvcice)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Flux of Rain [kg/m2.s]",":", one_gbx%mr_hydro(ii,jj,i_cvrain)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Effective Radius of Rain [m]",":",one_gbx%reff(ii,jj,i_cvrain)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Flux of Snow [kg/m2.s]",":",one_gbx%mr_hydro(ii,jj,i_cvsnow)
      write(0,fmt='(20x,a35,a1,1x,e10.3)') "Effective Radius of Snow [m]",":",one_gbx%reff(ii,jj,i_cvsnow)
    end do ! nlevels
    write(0,'(A55)') '-------------------------------------------------------'
    !call stop_model("nlevels conflicts model_lm",255)
    return
  end subroutine print_cosp_gbx
#endif

end module cosp_drv
! >>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<
! >>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<
! >>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<
