#include "rundeck_opts.h"

      module flammability_com    
!@sum for routines to calculate flammability potential of surface 
!@+   vegetation. Optionally also altering tracer biomass sources.
!@auth Greg Faluvegi based on direction from Olga Pechony including
!@+ her document Flammability.doc
!@ Modified by Keren Mezuman
!@var ij_flamV indicies for aij output 
      use ent_const, only : N_COVERTYPES

      implicit none
      save

! main variables:
!@var flammability unitless flammability coefficient
!@var vegetation density for purpose of flammability calculation
      real*8, allocatable, dimension(:,:) :: flammability,veg_density
!@param mfcc MODIS fire count calibration (goes with EPFC by veg type
!@+ params.) Units=fires/m2/yr when multiplied by the unitless flammability
      real*8, parameter :: mfcc=2.2d-5
      real*8, allocatable, dimension(:,:) :: EPFCByVegType 
      integer :: ij_flamV 
! rest is for the running average:
!@dbparam allowFlammabilityReinit (default 1=YES) allows the
!@+ flammability to initialize to undef when Itime=ItimeI and
!@+ the averaging period to start from the beginning (thus no emissions
!@+ for nday_prec days.) Override this with allowFlammabilityReinit=0
!@+ in the rundeck, and values from the AIC file should be used (if 
!@+ present.) 
!@var maxHR_prec maximum number of sub-daily accumulations
!@param nday_prec number of days in running average for prec
!@param nday_lai number of days in running average for Ent LAI
!@var DRAfl daily running average of model prec for flammability
!@var ravg_prec period running average of model prec for flammability
!@var PRSfl period running sum of model prec for flammability
!@var HRAfl hourly running average of model prec for flammability
!@var iHfl "hourly" index for averages of model prec
!@var iDfl "daily"  index for averages of model prec
!@var i0fl ponter to current index in running sum of model prec
!@var first_prec whether in the first model averaging per. for prec
!@var raP_acc accumulate running avg precip for SUBDD output
      integer, parameter :: nday_prec=30, nday_lai=365 ! unhappy with this making a large I,J array
      integer :: allowFlammabilityReinit = 1
      real*8, allocatable, dimension(:,:):: raP_acc

      real*8, allocatable, dimension(:,:,:):: DRAfl,burnt_area
      real*8, allocatable, dimension(:,:)  :: ravg_prec,PRSfl,iHfl,iDfl,
     &                                        i0fl,first_prec
      real*8, allocatable, dimension(:,:,:):: HRAfl
      integer :: maxHR_prec
      real*8, allocatable, dimension(:,:,:):: DRAlai
      real*8, allocatable, dimension(:,:):: ravg_lai,PRSlai,iHlai,iDlai,
     &                                      i0lai,first_lai
      real*8, allocatable, dimension(:,:,:):: HRAlai
      integer :: maxHR_lai
#ifdef ANTHROPOGENIC_FIRE_MODEL
      real*8, allocatable, dimension(:,:) :: populationDensity,flamPopA,
     &                                       flamPopB,nonSuppressFrac
      logical :: firstFlamPop = .true. 
      integer :: flamPopYearStart=2000, flamPopYearEnd=2000, 
     & flamPopDelYear=10
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
      real*8, allocatable, dimension(:,:) :: saveFireCount
#endif

      end module flammability_com


      subroutine alloc_flammability(grid)
!@SUM  alllocates arrays whose sizes need to be determined
!@+    at run-time
!@auth Greg Faluvegi
      use TRACER_COM, only: ntm
      use ent_const, only : N_COVERTYPES
      use domain_decomp_atm, only: dist_grid, getDomainBounds
      use dictionary_mod, only : get_param, is_set_param
      use model_com, only: dtsrc
      use TimeConstants_mod, only: SECONDS_PER_HOUR, HOURS_PER_DAY
      use flammability_com, only: flammability,veg_density,
     & first_prec,iHfl,iDfl,i0fl,DRAfl,ravg_prec,PRSfl,HRAfl,
     & nday_prec,maxHR_prec,raP_acc,burnt_area,EPFCByVegType
      use flammability_com, only: 
     & first_lai,iHlai,iDlai,i0lai,DRAlai,ravg_lai,PRSlai,HRAlai,
     & nday_lai,maxHR_lai
#ifdef ANTHROPOGENIC_FIRE_MODEL
      use flammability_com, only: populationDensity,flamPopA,flamPopB,
     & nonSuppressFrac
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
      use flammability_com, only: saveFireCount
#endif

      implicit none

      real*8 :: DTsrc_LOCAL
      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H

      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO 

      ! maxHR_prec is to be defined as the number of precip
      ! running-average steps that define one day. (The 1.d0
      ! in the line is left in to represent the number of calls 
      ! per DTsrc timestep. Same for Ent LAI:
      ! I believe the real DTsrc is not available yet so I use
      ! a local copy from the database:

      DTsrc_LOCAL = DTsrc
      if(is_set_param("DTsrc"))call get_param("DTsrc",DTsrc_LOCAL)
      maxHR_prec = NINT(HOURS_PER_DAY*1.d0*SECONDS_PER_HOUR/DTsrc_LOCAL)
      maxHR_lai  = NINT(HOURS_PER_DAY*1.d0*SECONDS_PER_HOUR/DTsrc_LOCAL)

      allocate( flammability(I_0H:I_1H,J_0H:J_1H) )
      allocate( veg_density (I_0H:I_1H,J_0H:J_1H) )
      allocate( first_prec  (I_0H:I_1H,J_0H:J_1H) )
      allocate( iHfl        (I_0H:I_1H,J_0H:J_1H) )
      allocate( iDfl        (I_0H:I_1H,J_0H:J_1H) )
      allocate( i0fl        (I_0H:I_1H,J_0H:J_1H) )
      allocate( DRAfl       (nday_prec,I_0H:I_1H,J_0H:J_1H) )
      allocate( ravg_prec   (I_0H:I_1H,J_0H:J_1H) )
      allocate( EPFCByVegType  (N_COVERTYPES,ntm) )
      allocate( burnt_area  (N_COVERTYPES,I_0H:I_1H,J_0H:J_1H) )
      allocate( PRSfl       (I_0H:I_1H,J_0H:J_1H) )
      allocate( raP_acc     (I_0H:I_1H,J_0H:J_1H) )
      allocate( HRAfl       (maxHR_prec,I_0H:I_1H,J_0H:J_1H) )
      allocate( first_lai   (I_0H:I_1H,J_0H:J_1H) )
      allocate( iHlai       (I_0H:I_1H,J_0H:J_1H) )
      allocate( iDlai       (I_0H:I_1H,J_0H:J_1H) )
      allocate( i0lai       (I_0H:I_1H,J_0H:J_1H) )
      allocate( DRAlai      (nday_lai,I_0H:I_1H,J_0H:J_1H) )
      allocate( ravg_lai    (I_0H:I_1H,J_0H:J_1H) )
      allocate( PRSlai      (I_0H:I_1H,J_0H:J_1H) )
      allocate( HRAlai      (maxHR_lai,I_0H:I_1H,J_0H:J_1H) )
#ifdef ANTHROPOGENIC_FIRE_MODEL
      allocate( nonSuppressFrac(I_0H:I_1H,J_0H:J_1H) )
      allocate( populationDensity(I_0H:I_1H,J_0H:J_1H) )
      allocate( flamPopA         (I_0H:I_1H,J_0H:J_1H) )
      allocate( flamPopB         (I_0H:I_1H,J_0H:J_1H) )
#endif /* ANTHROPOGENIC_FIRE_MODEL */
#ifdef DYNAMIC_BIOMASS_BURNING
      allocate( saveFireCount    (I_0H:I_1H,J_0H:J_1H) )
#endif

      return
      end subroutine alloc_flammability


      subroutine init_flammability
!@sum initialize flamability, veg density, etc. for fire model
!@auth Greg Faluvegi based on direction from Olga Pechony
      use model_com, only: Itime,ItimeI
      use constant, only: undef
      use OldTracer_mod, only: trname
      use TRACER_COM, only: whichEPFCs,ntm
      use ent_const, only: N_COVERTYPES
      use ent_pfts, only: ent_cover_names
      use dictionary_mod, only: sync_param
      use flammability_com, only: flammability, veg_density, first_prec
     & ,allowFlammabilityReinit,DRAfl,hrafl,prsfl,i0fl,iDfl,iHfl
     & ,ravg_prec,burnt_area,EPFCByVegType
      use flammability_com, only: first_lai,DRAlai,HRAlai,PRSlai,i0lai
     & ,iDlai,iHlai,ravg_lai
      use domain_decomp_atm,only: grid, getDomainBounds, 
     & am_i_root, readt_parallel
      use filemanager, only: openunit, closeunit, nameunit

      implicit none
      character*80 :: title,fname

      integer :: I_1H, I_0H, J_1H, J_0H, iu_data, n, v, tracerIndex
!!    EPFCByVegType tracer emissions per fire count as a
!!    function of ENT vegetation types

      call getDomainBounds(grid,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)
      call getDomainBounds(grid,I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

!! #ifdef FLAM_USE_OFFLINE_VEG_DENS
!!    Example of reading in one's own vegetation density and not using Ent
!!    call openunit('VEG_DENSE',iu_data,.true.,.true.)
!!    call readt_parallel(grid,iu_data,nameunit(iu_data),veg_density,1)
!!    call closeunit(iu_data)
!! #else
      veg_density(:,:)=undef ! defined later
!! #endif

      call sync_param
     &("allowFlammabilityReinit",allowFlammabilityReinit)

      if( (Itime==ItimeI .and. allowFlammabilityReinit==1) .or. 
     &  allowFlammabilityReinit == -1 )then
        flammability(:,:)=undef
        first_prec(:,:)=1.d0
        burnt_area(:,:,:)=0.d0
        DRAfl=0.d0
        hrafl=0.d0
        prsfl=0.d0
        i0fl=0.d0
        iDfl=0.d0
        iHfl=0.d0
        ravg_prec=0.d0
        first_lai(:,:)=1.d0
        DRAlai=0.d0
        HRAlai=0.d0
        PRSlai=0.d0
        i0lai=0.d0
        iDlai=0.d0
        iHlai=0.d0
        ravg_lai=0.d0
      end if
      !loop over tracers and have a select case for every tracer
      ! have multiple tracers at once e.g. BCB, M_BC1_BC, tomas one..
      ! make sure all the variables below are well defined at the
      ! beginning of the subroutine
      call sync_param("whichEPFCs",whichEPFCs)
      EPFCByVegType= 0.d0
      select case(whichEPFCs)!kgC/#fire/(vegfrac in gridbox)
        case(1) ! AR6 2001-2009
          do n=1,ntm
          select case(trname(n))
            case('M_BC1_BC','BCB')
              do v=1,N_COVERTYPES
              select case(ent_cover_names(v))
                case('arid_shrub')
                  EPFCByVegType(v,n)=238.d0
                case('c3_grass_ann')
                  EPFCByVegType(v,n)=173.d0
                case('c3_grass_arct')
                  EPFCByVegType(v,n)=1159.d0
                case('c3_grass_per')
                  EPFCByVegType(v,n)=257.d0
                case('c4_grass')
                  EPFCByVegType(v,n)=726.d0
                case('cold_br_late')
                  EPFCByVegType(v,n)=767.d0
                case('cold_shrub')
                  EPFCByVegType(v,n)=357.d0
                case('decid_nd')
                  EPFCByVegType(v,n)=1844.d0
                case('drought_br')
                  EPFCByVegType(v,n)=1382.d0
                case('ever_br_late')
                  EPFCByVegType(v,n)=1434.d0
                case('ever_nd_late')
                  EPFCByVegType(v,n)=821.d0
                case default
                  EPFCByVegType(v,n)=0.d0
              end select
              end do
            case('M_OCC_OC','OCB')
              do v=1,N_COVERTYPES
              select case(ent_cover_names(v))
                case('arid_shrub')
                  EPFCByVegType(v,n)=1479.d0
                case('c3_grass_ann')
                  EPFCByVegType(v,n)=728.d0
                case('c3_grass_arct')
                  EPFCByVegType(v,n)=15551.d0
                case('c3_grass_per')
                  EPFCByVegType(v,n)=1504.d0
                case('c4_grass')
                  EPFCByVegType(v,n)=4339.d0
                case('cold_br_late')
                  EPFCByVegType(v,n)=3437.d0
                case('cold_shrub')
                  EPFCByVegType(v,n)=6562.d0
                case('decid_nd')
                  EPFCByVegType(v,n)=36753.d0
                case('drought_br')
                  EPFCByVegType(v,n)=10667.d0
                case('ever_br_late')
                  EPFCByVegType(v,n)=10941.d0
                case('ever_nd_late')
                  EPFCByVegType(v,n)=6537.d0
                case default
                  EPFCByVegType(v,n)=0.d0
              end select
              end do
            case('NOx') 
              do v=1,N_COVERTYPES
              select case(ent_cover_names(v))
                case('arid_shrub')
                  EPFCByVegType(v,n)=1009.d0
                case('c3_grass_ann')
                  EPFCByVegType(v,n)=690.d0
                case('c3_grass_arct')
                  EPFCByVegType(v,n)=1094.d0
                case('c3_grass_per')
                  EPFCByVegType(v,n)=908.d0
                case('c4_grass')
                  EPFCByVegType(v,n)=3152.d0
                case('cold_br_late')
                  EPFCByVegType(v,n)=1529.d0
                case('cold_shrub')
                  EPFCByVegType(v,n)=241.d0
                case('decid_nd')
                  EPFCByVegType(v,n)=1559.d0
                case('drought_br')
                  EPFCByVegType(v,n)=4835.d0
                case('ever_br_late')
                  EPFCByVegType(v,n)=4905.d0
                case('ever_nd_late')
                  EPFCByVegType(v,n)=1197.d0
                case default
                  EPFCByVegType(v,n)=0.d0
              end select
              end do
            case('CO')
              do v=1,N_COVERTYPES
              select case(ent_cover_names(v))
                case('arid_shrub')
                  EPFCByVegType(v,n)=39268.d0
                case('c3_grass_ann')
                  EPFCByVegType(v,n)=26761.d0
                case('c3_grass_arct')
                  EPFCByVegType(v,n)=251702.d0
                case('c3_grass_per')
                  EPFCByVegType(v,n)=41043.d0
                case('c4_grass')
                  EPFCByVegType(v,n)=117577.d0
                case('cold_br_late')
                  EPFCByVegType(v,n)=113392.d0
                case('cold_shrub')
                  EPFCByVegType(v,n)=105936.d0
                case('decid_nd')
                  EPFCByVegType(v,n)=481485.d0
                case('drought_br')
                  EPFCByVegType(v,n)=230829.d0
                case('ever_br_late')
                  EPFCByVegType(v,n)=249906.d0
                case('ever_nd_late')
                  EPFCByVegType(v,n)=146622.d0
                case default
                  EPFCByVegType(v,n)=0.d0
              end select
              end do
            case('Alkenes')
              do v=1,N_COVERTYPES
              select case(ent_cover_names(v))
                case('arid_shrub')
                  EPFCByVegType(v,n)=36.6d0
                case('c3_grass_ann')
                  EPFCByVegType(v,n)=25.1d0
                case('c3_grass_arct')
                  EPFCByVegType(v,n)=489.d0
                case('c3_grass_per')
                  EPFCByVegType(v,n)=38.8d0
                case('c4_grass')
                  EPFCByVegType(v,n)=110.d0
                case('cold_br_late')
                  EPFCByVegType(v,n)=106.d0
                case('cold_shrub')
                  EPFCByVegType(v,n)=104.d0
                case('decid_nd')
                  EPFCByVegType(v,n)=422.d0
                case('drought_br')
                  EPFCByVegType(v,n)=214.d0
                case('ever_br_late')
                  EPFCByVegType(v,n)=220.d0
                case('ever_nd_late')
                  EPFCByVegType(v,n)=137.d0
                case default
                  EPFCByVegType(v,n)=0.d0
              end select
              end do
            case('Paraffin')
              do v=1,N_COVERTYPES
              select case(ent_cover_names(v))
                case('arid_shrub')
                  EPFCByVegType(v,n)=18.5d0
                case('c3_grass_ann')
                  EPFCByVegType(v,n)=13.9d0
                case('c3_grass_arct')
                  EPFCByVegType(v,n)=226.d0
                case('c3_grass_per')
                  EPFCByVegType(v,n)=20.7d0
                case('c4_grass')
                  EPFCByVegType(v,n)=57.0d0
                case('cold_br_late')
                  EPFCByVegType(v,n)=69.8d0
                case('cold_shrub')
                  EPFCByVegType(v,n)=72.1d0
                case('decid_nd')
                  EPFCByVegType(v,n)=373.d0
                case('drought_br')
                  EPFCByVegType(v,n)=108.d0
                case('ever_br_late')
                  EPFCByVegType(v,n)=102.d0
                case('ever_nd_late')
                  EPFCByVegType(v,n)=89.1d0
                case default
                  EPFCByVegType(v,n)=0.d0
              end select
              end do
            case('SO2')
              do v=1,N_COVERTYPES
              select case(ent_cover_names(v))
                case('arid_shrub')
                  EPFCByVegType(v,n)=262.d0
                case('c3_grass_ann')
                  EPFCByVegType(v,n)=147.d0
                case('c3_grass_arct')
                  EPFCByVegType(v,n)=2315.d0
                case('c3_grass_per')
                  EPFCByVegType(v,n)=270.d0
                case('c4_grass')
                  EPFCByVegType(v,n)=795.d0
                case('cold_br_late')
                  EPFCByVegType(v,n)=555.d0
                case('cold_shrub')
                  EPFCByVegType(v,n)=878.d0
                case('decid_nd')
                  EPFCByVegType(v,n)=4168.d0
                case('drought_br')
                  EPFCByVegType(v,n)=1687.d0
                case('ever_br_late')
                  EPFCByVegType(v,n)=1438.d0
                case('ever_nd_late')
                  EPFCByVegType(v,n)=972.d0
                case default
                  EPFCByVegType(v,n)=0.d0
              end select
              end do
            case('NH3')
              do v=1,N_COVERTYPES
              select case(ent_cover_names(v))
                case('arid_shrub')
                  EPFCByVegType(v,n)=378.d0
                case('c3_grass_ann')
                  EPFCByVegType(v,n)=313.d0
                case('c3_grass_arct')
                  EPFCByVegType(v,n)=5065.d0
                case('c3_grass_per')
                  EPFCByVegType(v,n)=438.d0
                case('c4_grass')
                  EPFCByVegType(v,n)=1196.d0
                case('cold_br_late')
                  EPFCByVegType(v,n)=2101.d0
                case('cold_shrub')
                  EPFCByVegType(v,n)=2006.d0
                case('decid_nd')
                  EPFCByVegType(v,n)=10722.d0
                case('drought_br')
                  EPFCByVegType(v,n)=2340.d0
                case('ever_br_late')
                  EPFCByVegType(v,n)=2847.d0
                case('ever_nd_late')
                  EPFCByVegType(v,n)=2277.d0
                case default
                  EPFCByVegType(v,n)=0.d0
              end select
              end do
            case('CH4')
              do v=1,N_COVERTYPES
              select case(ent_cover_names(v))
                case('arid_shrub')
                  EPFCByVegType(v,n)=1351.d0
                case('c3_grass_ann')
                  EPFCByVegType(v,n)=1013.d0
                case('c3_grass_arct')
                  EPFCByVegType(v,n)=11574.d0
                case('c3_grass_per')
                  EPFCByVegType(v,n)=1494.d0
                case('c4_grass')
                  EPFCByVegType(v,n)=4113.d0
                case('cold_br_late')
                  EPFCByVegType(v,n)=5915.d0
                case('cold_shrub')
                  EPFCByVegType(v,n)=4505.d0
                case('decid_nd')
                  EPFCByVegType(v,n)=22977.d0
                case('drought_br')
                  EPFCByVegType(v,n)=8461.d0
                case('ever_br_late')
                  EPFCByVegType(v,n)=10477.d0
                case('ever_nd_late')
                  EPFCByVegType(v,n)=6863.d0
                case default
                  EPFCByVegType(v,n)=0.d0
              end select
              end do
          end select
          end do
        case default
          call stop_model('whichEPFCs unknown',255)
      end select

      return
      end subroutine init_flammability


      subroutine def_rsf_flammability(fid)
!@sum  def_rsf_flammability defines flammability array structure in 
!@+    restart files
!@auth Greg Faluvegi (directly from M. Kelley's def_rsf_lakes)
!@ver  beta
      use flammability_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id

      call defvar(grid,fid,drafl,'drafl(nday_prec,dist_im,dist_jm)')
      call defvar(grid,fid,hrafl,'hrafl(maxHR_prec,dist_im,dist_jm)')
      call defvar(grid,fid,prsfl,'prsfl(dist_im,dist_jm)')
      call defvar(grid,fid,i0fl,'i0fl(dist_im,dist_jm)') ! real
      call defvar(grid,fid,iDfl,'iDfl(dist_im,dist_jm)') ! real
      call defvar(grid,fid,iHfl,'iHfl(dist_im,dist_jm)') ! real
      call defvar(grid,fid,first_prec,'first_prec(dist_im,dist_jm)')
      call defvar(grid,fid,ravg_prec,'ravg_prec(dist_im,dist_jm)')
      call defvar(grid,fid,burnt_area,'burnt_area(N_COVERTYPES,dist_im,
     &     dist_jm)')
      call defvar(grid,fid,flammability,'flammability(dist_im,dist_jm)')
      call defvar(grid,fid,raP_acc,'raP_acc(dist_im,dist_jm)')
      call defvar(grid,fid,dralai,'dralai(nday_lai,dist_im,dist_jm)')
      call defvar(grid,fid,hralai,'hralai(maxHR_lai,dist_im,dist_jm)')
      call defvar(grid,fid,prslai,'prslai(dist_im,dist_jm)')
      call defvar(grid,fid,i0lai,'i0lai(dist_im,dist_jm)') ! real
      call defvar(grid,fid,iDlai,'iDlai(dist_im,dist_jm)') ! real
      call defvar(grid,fid,iHlai,'iHlai(dist_im,dist_jm)') ! real
      call defvar(grid,fid,first_lai,'first_lai(dist_im,dist_jm)')
      call defvar(grid,fid,ravg_lai,'ravg_lai(dist_im,dist_jm)')

      return
      end subroutine def_rsf_flammability

      subroutine new_io_flammability(fid,iaction)
!@sum  new_io_flammability read/write arrays from/to restart files
!@auth Greg Faluvegi (directly from M. Kelley's new_io_lakes)
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use flammability_com
      implicit none
!@var fid unit number of read/write
      integer fid   
!@var iaction flag for reading or writing to file
      integer iaction 
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'drafl', drafl, jdim=3 )
        call write_dist_data(grid, fid, 'hrafl', hrafl, jdim=3 )
        call write_dist_data(grid, fid, 'prsfl', prsfl )
        call write_dist_data(grid, fid, 'i0fl', i0fl )
        call write_dist_data(grid, fid, 'iDfl', iDfl )
        call write_dist_data(grid, fid, 'iHfl', iHfl )
        call write_dist_data(grid, fid, 'first_prec', first_prec )
        call write_dist_data(grid, fid, 'ravg_prec', ravg_prec )
        call write_dist_data(grid, fid, 'burnt_area',burnt_area,jdim=3 )
        call write_dist_data(grid, fid, 'flammability', flammability )
        call write_dist_data(grid, fid, 'raP_acc', raP_acc )
        call write_dist_data(grid, fid, 'dralai', dralai, jdim=3 )
        call write_dist_data(grid, fid, 'hralai', hralai, jdim=3 )
        call write_dist_data(grid, fid, 'prslai', prslai )
        call write_dist_data(grid, fid, 'i0lai', i0lai )
        call write_dist_data(grid, fid, 'iDlai', iDlai )
        call write_dist_data(grid, fid, 'iHlai', iHlai )
        call write_dist_data(grid, fid, 'first_lai', first_lai )
        call write_dist_data(grid, fid, 'ravg_lai', ravg_lai )

      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'drafl', drafl, jdim=3 )
        call read_dist_data(grid, fid, 'hrafl', hrafl, jdim=3 )
        call read_dist_data(grid, fid, 'prsfl', prsfl )
        call read_dist_data(grid, fid, 'i0fl', i0fl )
        call read_dist_data(grid, fid, 'iDfl', iDfl )
        call read_dist_data(grid, fid, 'iHfl', iHfl )
        call read_dist_data(grid, fid, 'first_prec', first_prec )
        call read_dist_data(grid, fid, 'ravg_prec', ravg_prec )
        call read_dist_data(grid, fid, 'burnt_area',burnt_area,jdim=3 )
        call read_dist_data(grid, fid, 'flammability', flammability )
        call read_dist_data(grid, fid, 'raP_acc', raP_acc )
        call read_dist_data(grid, fid, 'dralai', dralai, jdim=3 )
        call read_dist_data(grid, fid, 'hralai', hralai, jdim=3 )
        call read_dist_data(grid, fid, 'prslai', prslai )
        call read_dist_data(grid, fid, 'i0lai', i0lai )
        call read_dist_data(grid, fid, 'iDlai', iDlai )
        call read_dist_data(grid, fid, 'iHlai', iHlai )
        call read_dist_data(grid, fid, 'first_lai', first_lai )
        call read_dist_data(grid, fid, 'ravg_lai', ravg_lai )
      end select
      return
      end subroutine new_io_flammability

      subroutine flammability_drv
!@sum driver routine for flammability potential of surface
!@+   vegetation calculation.
!@auth Greg Faluvegi based on direction from Olga Pechony
!@+ Later modified by Keren Mezuman
!@ver  1.0 
      use geom, only: axyp
      use model_com, only: dtsrc
      use resolution, only : jm
      use atm_com, only : pedn,Q,PMID,pk,t
      use domain_decomp_atm,only: grid, getDomainBounds
      use flammability_com, only: flammability,veg_density,ravg_prec,
     & iHfl,iDfl,i0fl,first_prec,HRAfl,DRAfl,PRSfl,raP_acc,
     & saveFireCount,burnt_area

      use fluxes, only: prec,atmsrf
      use constant, only: lhe, undef
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use diag_com, only: ij_flam,ij_flam_rh,ij_flam_prec,ij_flam_tsurf,
     &  ij_barh1,ij_bawsurf,aij=>aij_loc
      use flammability_com, only: ravg_lai,iHlai,iDlai,i0lai,
     & first_lai,HRAlai,DRAlai,PRSlai
      use ent_const, only : N_COVERTYPES
      use ent_pfts, only: ent_cover_names
      use ghy_com, only: fearth
      use diag_com, only: ij_fvden,ij_ba_tree,ij_ba_shrub,ij_ba_grass
      use ent_com, only: entcells
      use ent_mod, only: ent_get_exports
#ifdef CACHED_SUBDD
      use subdd_mod, only : subdd_groups,subdd_ngroups,subdd_type
     &     ,inc_subdd,find_groups
#endif

      implicit none

#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups),k
      type(subdd_type), pointer :: subdd
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: sddarr2d
#endif
      integer :: J_0S, J_1S, I_0H, I_1H, i, j, v
      logical :: have_south_pole, have_north_pole     
      real*8 :: qsat ! this is a function in UTILDBL.f
      real*8 :: tsurf,qsurf,RH1,wsurf,fearth_axyp,TK
      ! the 7.9 here was from running a year or two under 2005 conditions
      ! and seeing what was the maximum LAI returned by Ent. Therefore,
      ! under other climate conditions, the vegetation density may reach > 1.0. 
      ! Olga thought this would not cause any mathematical probems. 
      ! I am however limiting the flammability to 1.0, but normally its
      ! values seem to be much much lower than that anyway...
      ! This seemed to result in too much emissions, so trying 10.0, which
      ! is the MODIS-based number Olga used in her original offline VEG_DENS
      ! file: 
      real*8, parameter :: byLaiMax=1.d0/10.0d0 !! 7.9d0
      real*8 :: lai
!@var pvt fraction vegetation type for N_COVERTYPES (fraction)
!@var N_COVERTYPES = N_PFT + N_SOILCOV + N_OTHER = 16 + 2 + 0
      real*8, dimension(N_COVERTYPES):: pvt
      real*8 :: fracVegNonCrops, fracBare
      real*8, parameter :: critFracBare = 0.8d0 ! 80% of box is bare soils

      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE = have_south_pole,
     &               HAVE_NORTH_POLE = have_north_pole)
      call getDomainBounds(grid,I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

      if(have_north_pole)flammability(I_0H:I_1H,JM)=undef
      if(have_south_pole)flammability(I_0H:I_1H,1) =undef

      do j=J_0S,J_1S
        do i=I_0H,I_1H

          ! update the precipitation running average:
          call prec_running_average(prec(i,j),ravg_prec(i,j), 
     &    iHfl(i,j),iDfl(i,j),i0fl(i,j),first_prec(i,j),HRAfl(:,i,j),
     &    DRAfl(:,i,j),PRSfl(i,j))
          ! for sub-daily diag purposes, accumulate the running avg:
          raP_acc(i,j)=raP_acc(i,j)+ravg_prec(i,j)
          ! if the first period has elapsed, calculate the flammability
          if(first_prec(i,j)/=0.) cycle

          ! and the LAI running average from Ent:
          if(fearth(i,j)>0.d0) then
            call ent_get_exports( entcells(i,j),leaf_area_index=lai)
            ! I guess that is the lai from the last surface timestep only?
            ! (But Igor says this is OK as LAI is only computed once per day.)
          else
            lai=0.d0
          end if
          call lai_running_average(lai,ravg_lai(i,j), 
     &    iHlai(i,j),iDlai(i,j),i0lai(i,j),first_lai(i,j),HRAlai(:,i,j),
     &    DRAlai(:,i,j),PRSlai(i,j))


!! #ifndef FLAM_USE_OFFLINE_VEG_DENS /* NOT */
!! I.e. do not define/limit the veg_density here if it is prescribed...
          tsurf = atmsrf%tsavg(i,j)![K]
          qsurf = atmsrf%qsavg(i,j)
          if(fearth(i,j)>0.d0) then
            if(first_lai(i,j)==0.) then
              veg_density(i,j) = ravg_lai(i,j)*byLaiMax*fearth(i,j)
            else
              ! for the first year of a run (since we don't have an annual 
              ! average yet, use the concurrent LAI):
              veg_density(i,j) =  lai*byLaiMax*fearth(i,j)
            end if
#ifdef LIMIT_BARREN_FLAMMABILITY
            ! Because of unrealistic LAI (therefore veg density) in deserts
            ! due to crops+pasture cover, we need to set the veg densitry to zero
            ! when either a box has 80% or more bare soil (light+dark) or the box 
            ! had zero non-crops vegetation:
            fracVegNonCrops=0.d0
            fracBare=0.d0
            call ent_get_exports(entcells(i,j),
     &         vegetation_fractions=pvt)
            do v=1,N_COVERTYPES
            select case(ent_cover_names(v))
              case('bare_bright')
                fracBare = fracBare + pvt(v) * fearth(i,j) 
              case('bare_dark')
                fracBare = fracBare + pvt(v) * fearth(i,j) 
              case('c3_grass_ann')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('c3_grass_arct')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('c3_grass_per')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('c4_grass')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('arid_shrub')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('cold_shrub')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('cold_br_late')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('decid_nd')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('drought_br')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('ever_br_late')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
              case('ever_nd_late')
                fracVegNonCrops = fracVegNonCrops + pvt(v) *
     &          fearth(i,j)
            end select
            end do

            if (fracBare + fracVegNonCrops > fearth(i,j)+0.00001) then
              call stop_model('cover types sum greater than'// 
     &                         'fearth',255) 
            endif
            if((fracBare >= critFracBare).or.(fracVegNonCrops == 0.)) 
     &        veg_density(i,j) = 0.d0
!@var atmsrf contains atm-surf interaction
!@+ quantities averaged over all surface types
!@var WSAVG SURFACE WIND MAGNITUDE (M/S)
!@var wsurf surface wind velocity (m/s)
            wsurf=atmsrf%wsavg(i,j)
!@var SECONDS_PER_DAY (s/day)
!@var DTSRC source time step (s) 
!@var TK ?
!@var lhe ?
!@var  PMID Pressure at mid point of box (mb)
!@var QSAT saturation specific humidity (?)
!@var Q specific humidity (kg water vapor/kg air)
!@var RH1 relative humidity in layer 1 (fraction)
!@var RH1 relative humidity at surface (fraction)
            TK = pk(1,i,j)*t(i,j,1)           !should be in [K]
            RH1=min(1.d0,qsurf/qsat(tsurf,lhe,pedn(1,i,j)))![fraction]
!@var fearth soil covered land fraction (fraction)
!@var  axyp,byaxyp area of grid box (+inverse) (m^2)
            fearth_axyp=fearth(i,j)*axyp(i,j)
            ! update the burnt area  average:
!@var burnt_area (m^2)
!@var saveFireCount fire count rate (fire/m2/s)
            aij(i,j,ij_barh1)=aij(i,j,ij_barh1)+RH1
            aij(i,j,ij_bawsurf)=aij(i,j,ij_bawsurf)+wsurf
            call step_ba(burnt_area(:,i,j),RH1,wsurf,
     &                   saveFireCount(i,j),pvt,fearth_axyp,i,j)
#endif /* LIMIT_BARREN_FLAMMABILITY */
          else
            veg_density(i,j) = 0.d0
          end if
!! #endif /* FLAM_USE_OFFLINE_VEG_DENS NOT DEFINED */

          call calc_flammability(tsurf,SECONDS_PER_DAY
     &     *ravg_prec(i,j)/dtsrc,min(1.d0,qsurf/
     &     qsat(tsurf,lhe,pedn(1,i,j))),veg_density(i,j),
     &     flammability(i,j),sum(burnt_area(:,i,j)),fearth_axyp)
          ! update diagnostic
          aij(i,j,ij_flam)=aij(i,j,ij_flam)+flammability(i,j)
          aij(i,j,ij_flam_tsurf)=aij(i,j,ij_flam_tsurf)+tsurf
          aij(i,j,ij_flam_prec)=aij(i,j,ij_flam_prec)+
     &      SECONDS_PER_DAY*ravg_prec(i,j)/dtsrc
          aij(i,j,ij_flam_rh)=aij(i,j,ij_flam_rh)+
     &      min(1.d0,qsurf/qsat(tsurf,lhe,pedn(1,i,j)))
          aij(i,j,ij_fvden)=aij(i,j,ij_fvden)+veg_density(i,j)
          do v=1,N_COVERTYPES
          select case(ent_cover_names(v))
            case('arid_shrub')
              aij(i,j,ij_ba_shrub) = aij(i,j,ij_ba_shrub) +
     &        burnt_area(v,i,j)
            case('c3_grass_ann')
              aij(i,j,ij_ba_grass) = aij(i,j,ij_ba_grass) +
     &        burnt_area(v,i,j)
            case('c3_grass_arct')
              aij(i,j,ij_ba_grass) = aij(i,j,ij_ba_grass) +
     &        burnt_area(v,i,j)
            case('c3_grass_per')
              aij(i,j,ij_ba_grass) = aij(i,j,ij_ba_grass) +
     &        burnt_area(v,i,j)
            case('c4_grass')
              aij(i,j,ij_ba_grass) = aij(i,j,ij_ba_grass) +
     &        burnt_area(v,i,j)
            case('cold_br_late')
              aij(i,j,ij_ba_tree) =  aij(i,j,ij_ba_tree) +
     &        burnt_area(v,i,j)
            case('cold_shrub')
              aij(i,j,ij_ba_shrub) = aij(i,j,ij_ba_shrub) + 
     &        burnt_area(v,i,j)
            case('decid_nd')
              aij(i,j,ij_ba_tree) = aij(i,j,ij_ba_tree) + 
     &        burnt_area(v,i,j)
            case('drought_br')
              aij(i,j,ij_ba_tree) = aij(i,j,ij_ba_tree) + 
     &        burnt_area(v,i,j)
            case('ever_br_late')
              aij(i,j,ij_ba_tree) = aij(i,j,ij_ba_tree) + 
     &        burnt_area(v,i,j)
            case('ever_nd_late')
              aij(i,j,ij_ba_tree) = aij(i,j,ij_ba_tree) + 
     &        burnt_area(v,i,j)
          end select
          end do

        end do
      end do
  
#ifdef CACHED_SUBDD
****
**** Collect some high-frequency outputs
      call find_groups('aijh',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          select case (subdd%name(k))
          case ('FLAMM')
            sddarr2d = flammability(:,:)
            call inc_subdd(subdd,k,sddarr2d)
          case ('FLAMM_prec')
            sddarr2d = SECONDS_PER_DAY*ravg_prec(:,:)/dtsrc
            call inc_subdd(subdd,k,sddarr2d)
          case ('FLAMM_rh')
            sddarr2d = min(1.d0,qsurf/qsat(tsurf,lhe,pedn(1,:,:)))
            call inc_subdd(subdd,k,sddarr2d)
          case ('FVDEN')
            sddarr2d = veg_density(:,:)
            call inc_subdd(subdd,k,sddarr2d)
          end select
        enddo
      enddo
#endif
      return
      end subroutine flammability_drv

