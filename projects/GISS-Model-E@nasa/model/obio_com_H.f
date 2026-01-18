#include "rundeck_opts.h"
      MODULE obio_com
!@sum  obio_com contains the parameters, arrays and definitions
!@+    necessary for the OceanBiology routines
!@auth NR

      USE obio_dim
      use ocalbedo_mod, only: nlt
!     USE Constant, only: sday     ! sday=86400.0    !seconds per day

      USE hycom_dim_glob, only: kdm
      use hycom_dim, only: idm,jdm,ntrcr

      implicit none

c --- dobio       activate Watson Gregg's ocean biology code
      logical dobio
      data dobio/.true./
c

      real, ALLOCATABLE, DIMENSION(:,:)    :: tzoo2d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tfac3d,wshc3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: Fescav3d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: rmuplsr3d,rikd3d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: acdom3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: gcmax         !cocco max growth rate
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2tot_day    !net pp total per day
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2diat_day
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2chlo_day
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2cyan_day
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2cocc_day
      real, ALLOCATABLE, DIMENSION(:,:)    :: tot_chlo      !tot chlorophyl at surf. layer
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: rhs_obio      !rhs matrix
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: chng_by       !integr tendency for total C

      real, ALLOCATABLE, DIMENSION(:,:) :: pCO2av, pCO2av_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2tot_dayav
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2tot_dayav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: cexpav, cexpav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: caexpav, caexpav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: ao_co2fluxav,ao_co2fluxav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2diat_dayav
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2diat_dayav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2chlo_dayav
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2chlo_dayav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2cyan_dayav
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2cyan_dayav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2cocc_dayav
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2cocc_dayav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pHav,pHav_loc
      
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: tracer
      real, ALLOCATABLE, DIMENSION(:,:)::  Edz,Euz,Esz
      real, ALLOCATABLE, DIMENSION(:,:)::  Kd       !absorption+scattering in seawater due to chl
      real, ALLOCATABLE, DIMENSION(:)::  Kpar     !kpar from NBOM
      real, ALLOCATABLE, DIMENSION(:)::  delta_temp1d  !change in T due to kpar

      integer :: nstep0=0

      integer:: num_tracers

      !test point
!!    integer, parameter :: itest=16, jtest=45    !equatorial Pacific                  2deg ocean
!!    integer, parameter :: itest=32, jtest=20    !southern ocean; Pacific          
!     integer, parameter :: itest=(220,320) equator Atlant; (245,275) 0.6S;274.5E Nino3
!     integer, parameter :: itest=316, jtest=258    !257.5E;-50.7S
      integer, parameter :: itest=243, jtest=1      !equator,dateline
      real :: diag_counter


      integer, parameter :: EUZ_DEFINED=1

      real, parameter :: sday=86400.0
      real, parameter :: rlamz=1.0            !Ivlev constant
      real, parameter :: greff=0.25           !grazing efficiency     
!     real, parameter :: drate=0.05/24.0      !phytoplankton death rate/hr
      real, parameter :: drate=0.05/sday      !phytoplankton death rate/s    !July 2016
!     real, parameter :: dratez1=0.1/24.0     !zooplankton death rate/hr
      real, parameter :: dratez1=0.1/sday     !zooplankton death rate/s      !July 2016
!     real, parameter :: dratez2=0.5/24.0     !zooplankton death rate/hr
      real, parameter :: dratez2=0.5/sday     !zooplankton death rate/s      !July 2016
      real, parameter :: regen=0.25           !regeneration fraction

      real ::  obio_deltath,obio_deltat       !time steps in s  !July 2016
      real ::  co2mon(26,12)        !26 years 1979-2004, 12 months

      integer npst,npnd   !starting and ending array index for PAR
      data npst,npnd /3,17/

      real WtoQ(nlt)           !Watts/m2 to quanta/m2/s conversion

! reduced rank arrays for obio_model calculations
      real cexp, caexp
      real temp1d(kdm),dp1d(kdm),obio_P(kdm,ntyp)
     .                 ,det(kdm,ndet),car(kdm,ncar),avgq1d(kdm)
     .                 ,gcmax1d(kdm),saln1d(kdm),p1d(kdm+1)
     .                 ,alk1d(kdm),flimit(kdm,nchl,5)

      real atmFe_ij,covice_ij

      integer inwst,inwnd,jnwst,jnwnd     !starting and ending indices 
                                          !for daylight 
                                          !in i and j directions
      real acdom(kdm,nlt)                 !absorption coefficient of CDOM
      real P_tend(kdm,ntyp)               !bio tendency (dP/dt)

#ifdef TRACERS_Alkalinity
      real A_tend(kdm), co3_conc
      real ca_det_calc1d(kdm),Ca_tend(kdm)
#endif
      real rmuplsr(kdm,nchl)                  !growth+resp 
      real D_tend(kdm,ndet)                   !detrtial tendency
      real obio_ws(kdm+1,nchl)                !phyto sinking rate
      real tfac(kdm)                          !phyto T-dependence
      real pnoice(kdm)                        !pct ice-free
      real wsdet(kdm+1,ndet)                  !detrital sinking rate
      real rikd(kdm,nchl)                     !photoadaption state
      real tzoo                               !herbivore T-dependence
      real Fescav(kdm)                        !iron scavenging rate

C if NCHL_DEFINED > 3
      real wshc(kdm)                          !cocco sinking rate
C endif

      real :: C_tend(kdm,ncar)                !carbon tendency
      real :: pCO2_ij,pHsfc                   !partial pressure of CO2, pH
      real :: gro(kdm,nchl)                   !realized growth rate
      integer :: day_of_month, hour_of_day

      real :: rhs(kdm,ntrac,17)         !secord arg-refers to tracer definition 

      real :: pp2_1d(kdm,nchl)          !net primary production

      real*8 :: co2flux
      integer kzc
      real*8 :: carb_old,iron_old    !prev timesetep total carbon inventory

#ifdef restoreIRON
!this is an AR5 preprocessor option
!per year change dI/I, + for sink, - for source
      !!!real*8 :: Iron_BC = 0.002
      real*8 :: Iron_BC = -0.005
#endif

      real*8, dimension(:, :, :), allocatable :: ze

      character(len=50) :: arg2d, arg3d


      contains

      subroutine build_ze
      use hycom_dim, only: ogrid, kdm
      USE hycom_arrays, only : dpinit
      USE hycom_scalars, only: onem
      
      implicit none
      integer :: k
      real :: pres,g,s
      integer :: i,j
      real*8 :: volgsp

      if (.not.allocated(ze)) allocate(ze(ogrid%i_strt:ogrid%i_stop,
     &                                ogrid%j_strt:ogrid%j_stop, 0:kdm))
      ze=0
      
      do k=1, kdm
        ze(:, :, k)=ze(:, :, k-1)+dpinit(ogrid%i_strt:ogrid%i_stop,
     &                 ogrid%j_strt:ogrid%j_stop, k)/onem
      end do
      
      end subroutine build_ze

      END MODULE obio_com



!------------------------------------------------------------------------------
      module vector_str30_mod
#define _entry character(30)
#include "containers/vector.fh"
      end module vector_str30_mod

      module vector_str80_mod
#define _entry character(80)
#include "containers/vector.fh"
      end module vector_str80_mod


      module obio_diag
      use vector_str30_mod, only: vector_str30=>vector
      use vector_str80_mod, only: vector_str80=>vector
      use vector_integer_mod, only: vector_integer=>vector
      use vector_real8_mod, only: vector_real8=>vector
      use cdl_mod, only: cdl_type
      use obio_dim, only: ntrac
      implicit none
      private

      public :: init_obio_diag, new_io_obio_diag, def_meta_obio_diag,
     &   write_meta_obio_diag, def_rsf_obio_diag, reset_obio_diag,
     &   add_diag

      real*8, dimension(:, :, :), allocatable, public :: obio_ij
      real*8, dimension(:, :, :, :), allocatable, public :: obio_ijl
      integer, public :: ij_solz, ij_sunz, ij_dayl, ij_ed, ij_es,
     &   ij_nitr, ij_amm, ij_sil, ij_iron, ij_diat, ij_chlo, ij_cyan,
     &   ij_cocc, ij_herb, ij_doc, ij_dic, ij_pco2, ij_alk, ij_flux,
     &   ij_cexp, ij_ndet, ij_setl, ij_sink, ij_xchl, ij_fca, 
     &   ij_rnitrmflo,
     &   ij_rnitrconc, ij_rdicconc, ij_rdocconc, ij_rsiliconc,
     &   ij_rironconc, ij_rpocconc, ij_ralkconc, ij_pp, ij_lim(4, 5),
     &   ij_rhs(ntrac, 17), ij_pp1, ij_pp2, ij_pp3, ij_pp4, ij_co3,
     &   ij_ph
      integer, public :: ijl_avgq, ijl_kpar, ijl_dtemp
      type(vector_str30) :: sname_ij, units_ij
      type(vector_str30) :: sname_ijl, units_ijl
      type(vector_str80) :: lname_ij, lname_ijl
      type(vector_integer) :: ia_ij, ia_ijl
      type(vector_real8) :: scale_ij, scale_ijl
      type(cdl_type) :: cdl_ij, cdl_ijl
      type(cdl_type), pointer :: cdl_lons, cdl_lats, cdl_depths
      type(cdl_type), allocatable, target ::
     &                         hycom_lons, hycom_lats, hycom_depths
      
      contains


      subroutine add_diag(lname, sname, units, dim3, idx)
      use mdiag_com, only : ia_cpl
      implicit none

      character(len=*), target, intent(in) :: lname, sname, units
      logical, intent(in) :: dim3
      integer, intent(out) :: idx
      character(len=30) :: sname1, units1
      character(len=80) :: lname1

      lname1=lname
      sname1=sname
      units1=units
      if (dim3) then
        call lname_ijl%push_back(lname1)
        call sname_ijl%push_back(sname1)
        call units_ijl%push_back(units1)
        call scale_ijl%push_back(1.d0)
        call ia_ijl%push_back(ia_cpl)
        idx=lname_ijl%getsize()
      else
        call lname_ij%push_back(lname1)
        call sname_ij%push_back(sname1)
        call units_ij%push_back(units1)
        call scale_ij%push_back(1.d0)
        call ia_ij%push_back(ia_cpl)
        idx=lname_ij%getsize()
      endif
      end subroutine add_diag


      subroutine init_obio_diag
      
      USE hycom_dim, only: kdm, ogrid
      use cdl_mod, only: add_coord, init_cdl_type
      use domain_decomp_1d, only: am_i_root
      use obio_com, only: arg2d, arg3d

      implicit none

      arg2d='idm,dist_jdm'
      arg3d='idm,dist_jdm,kdm'
      allocate(hycom_lons, hycom_lats, hycom_depths)
      cdl_lons=>hycom_lons
      cdl_lats=>hycom_lats
      cdl_depths=>hycom_depths
      if (am_i_root()) then
            ! cdl not defined for hycom yet, this is for the time being:
        call init_cdl_type('cdl_obio_lons', hycom_lons)
        call add_coord(hycom_lons, 'lono', ogrid%im_world,
     &      units='degrees_east')
        call init_cdl_type('cdl_obio_lats', hycom_lats)
        call add_coord(hycom_lats, 'lato', ogrid%jm_world,
     &      units='degrees_north')
        call init_cdl_type('cdl_obio_depths', hycom_depths)
        call add_coord(hycom_depths, 'zoc', kdm, units='m')
      endif
      
      allocate(obio_ij(ogrid%i_strt:ogrid%i_stop,
     &         ogrid%j_strt:ogrid%j_stop, lname_ij%getsize()))
      allocate(obio_ijl(ogrid%i_strt:ogrid%i_stop,
     &         ogrid%j_strt:ogrid%j_stop, kdm, lname_ijl%getsize()))

      end subroutine init_obio_diag


      subroutine def_rsf_obio_diag(fid, r4_on_disk)
      use pario, only: defvar
      use obio_com, only: arg2d, arg3d
      USE hycom_dim, only: ogrid
      
      implicit none

      integer, intent(in) :: fid
      logical, intent(in) :: r4_on_disk

      call defvar(ogrid, fid, obio_ij,
     &   'obio_ij('//trim(arg2d)//',kobio_ij)', r4_on_disk=r4_on_disk)
      call defvar(ogrid, fid, obio_ijl,
     &   'obio_ijl('//trim(arg3d)//',kobio_ijl)', r4_on_disk=r4_on_disk)

      end subroutine def_rsf_obio_diag


      subroutine new_io_obio_diag(fid, iaction)
      use model_com, only: ioread
      use pario, only: write_dist_data, read_dist_data
      USE hycom_dim, only: ogrid
      
      implicit none

      integer, intent(in) :: fid
      integer, intent(in) :: iaction

      if (iaction.eq.ioread) then
        call read_dist_data(ogrid, fid, 'obio_ij', obio_ij)
        call read_dist_data(ogrid, fid, 'obio_ijl', obio_ijl)
      else
        call write_dist_data(ogrid, fid, 'obio_ij', obio_ij)
        call write_dist_data(ogrid, fid, 'obio_ijl', obio_ijl)
      end if

      end subroutine new_io_obio_diag


      subroutine def_meta_obio_diag(fid)
      use pario, only: defvar, write_attr
      use cdl_mod, only: defvar_cdl, merge_cdl, add_var
      use domain_decomp_1d, only: am_i_root
      USE hycom_dim, only: ogrid, kdm
      use hycom_arrays, only : lonij, latij
      use obio_com, only: ze, build_ze, arg2d, arg3d
      
      implicit none

      integer, intent(in) :: fid
      integer :: k

      if (associated(cdl_lons)) then
        if (am_i_root()) then
          call merge_cdl(cdl_lons, cdl_lats, cdl_ij)
          call add_var(cdl_ij, 'float latij(lato,lono) ;',
     &         long_name='gridbox latitude', units='degrees')
          call add_var(cdl_ij, 'float lonij(lato,lono) ;',
     &         long_name='gridbox longitude', units='degrees')
          call merge_cdl(cdl_ij, cdl_depths, cdl_ijl)
          call add_var(cdl_ijl, 'float depths(lato,lono,zoc) ;',
     &         long_name='gridbox depth', units='m')
          
       do k=1, sname_ij%getsize()
            call add_var(cdl_ij,
     &         'float '//trim(sname_ij%at(k))//'(lato,lono) ;',
     &         long_name=trim(lname_ij%at(k)),
     &         units=trim(units_ij%at(k)) )
          enddo
          do k=1, sname_ijl%getsize()
            call add_var(cdl_ijl,
     &         'float '//trim(sname_ijl%at(k))//'(zoc,lato,lono) ;',
     &         long_name=trim(lname_ijl%at(k)),
     &         units=trim(units_ijl%at(k)),
     &         set_miss=.true.)
          enddo
        endif
        call defvar_cdl(ogrid, fid, cdl_ij,
     &                  'cdl_obio_ij(cdl_strlen,kcdl_obio_ij)')
        call defvar_cdl(ogrid, fid, cdl_ijl,
     &       'cdl_obio_ijl(cdl_strlen,kcdl_obio_ijl)')
      endif

      call write_attr(ogrid, fid, 'obio_ij', 'reduction', 'sum')
      call write_attr(ogrid, fid, 'obio_ij', 'split_dim', 3)
      call defvar(ogrid, fid, ia_ij%getdata(), 'ia_obio_ij(kobio_ij)')
      call defvar(ogrid, fid, scale_ij%getdata(),
     &            'scale_obio_ij(kobio_ij)')
      call defvar(ogrid, fid, sname_ij%getdata(),
     &            'sname_obio_ij(sname_strlen,kobio_ij)')


      call write_attr(ogrid, fid, 'obio_ijl', 'reduction', 'sum')
      call write_attr(ogrid, fid, 'obio_ijl', 'split_dim', 4)
      call defvar(ogrid, fid, ia_ijl%getdata(),
     &            'ia_obio_ijl(kobio_ijl)')
      call defvar(ogrid, fid, scale_ijl%getdata(),
     &            'scale_obio_ijl(kobio_ijl)')
      call defvar(ogrid, fid, sname_ijl%getdata(),
     &            'sname_obio_ijl(sname_strlen,kobio_ijl)')
      call defvar(ogrid, fid,latij(:,:,3), 'latij('//trim(arg2d)//')')
      call defvar(ogrid, fid,lonij(:,:,3), 'lonij('//trim(arg2d)//')')
      call build_ze
      call defvar(ogrid,fid,ze(:,:,1:), 'depths('//trim(arg3d)//')')

      end subroutine def_meta_obio_diag


      subroutine write_meta_obio_diag(fid)
      use pario, only: write_data, write_dist_data
      use cdl_mod, only: write_cdl
      USE hycom_dim, only: ogrid
      use hycom_arrays, only : lonij, latij
      use obio_com, only: ze
      
      implicit none

      integer, intent(in) :: fid

      call write_dist_data(ogrid, fid, 'latij', latij(:, :, 3))
      call write_dist_data(ogrid, fid, 'lonij', lonij(:, :, 3))
      call write_dist_data(ogrid, fid, 'depths', ze)
      call write_data(ogrid, fid, 'ia_obio_ij', ia_ij%getdata())
      call write_data(ogrid, fid, 'scale_obio_ij', scale_ij%getdata())
      call write_data(ogrid, fid, 'sname_obio_ij', sname_ij%getdata())
      if (associated(cdl_lons))
     &      call write_cdl(ogrid, fid, 'cdl_obio_ij', cdl_ij)

      call write_data(ogrid, fid, 'ia_obio_ijl', ia_ijl%getdata())
      call write_data(ogrid, fid, 'scale_obio_ijl', scale_ijl%getdata())
      call write_data(ogrid, fid, 'sname_obio_ijl', sname_ijl%getdata())
      if (associated(cdl_lons))
     &      call write_cdl(ogrid, fid, 'cdl_obio_ijl', cdl_ijl)

      end subroutine write_meta_obio_diag


      subroutine reset_obio_diag
      implicit none

      obio_ij=0.
      obio_ijl=0.

      end subroutine reset_obio_diag


      end module obio_diag



!------------------------------------------------------------------------------
      subroutine alloc_obio_com

      USE obio_com
      USE obio_dim
      use obio_diag, only: init_obio_diag

      USE hycom_dim, only: kdm,ogrid

      implicit none

c**** Extract domain decomposition info
      INTEGER :: j_0,j_1,i_0,i_1

      I_0 = ogrid%I_STRT
      I_1 = ogrid%I_STOP
      J_0 = ogrid%J_STRT
      J_1 = ogrid%J_STOP


      ALLOCATE(tracer(i_0:i_1,j_0:j_1,kdm,ntrac))

!NOT FOR HYCOM:idm and jdm were not originally passes to hycom
      call alloc_obio_forc(kdm,ogrid,idm,jdm)

      ALLOCATE(tzoo2d(i_0:i_1,j_0:j_1))
      ALLOCATE(wshc3d(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(Fescav3d(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(rmuplsr3d(i_0:i_1,j_0:j_1,kdm,nchl),
     &            rikd3d(i_0:i_1,j_0:j_1,kdm,nchl))
      ALLOCATE(acdom3d(i_0:i_1,j_0:j_1,kdm,nlt))
      ALLOCATE(tfac3d(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(gcmax(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(pp2tot_day(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2diat_day(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2chlo_day(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2cyan_day(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2cocc_day(i_0:i_1,j_0:j_1))
      ALLOCATE(tot_chlo(i_0:i_1,j_0:j_1))
      ALLOCATE(rhs_obio(i_0:i_1,j_0:j_1,ntrac,17))
      ALLOCATE(chng_by(i_0:i_1,j_0:j_1,14))

      ALLOCATE(pCO2av(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(pp2tot_dayav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(cexpav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(caexpav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(pp2diat_dayav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(pp2chlo_dayav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(pp2cyan_dayav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(pp2cocc_dayav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(pHav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(ao_co2fluxav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(pCO2av_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2tot_dayav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(cexpav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(caexpav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(ao_co2fluxav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2diat_dayav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2chlo_dayav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2cyan_dayav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2cocc_dayav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(pHav_loc(i_0:i_1,j_0:j_1))

      ALLOCATE(Edz(nlt,kdm))
      ALLOCATE(Esz(nlt,kdm))
      ALLOCATE(Euz(nlt,kdm))
      ALLOCATE(Kd(nlt,kdm))
      ALLOCATE(Kpar(kdm))
      ALLOCATE(delta_temp1d(kdm))

      call init_obio_diag

      end subroutine alloc_obio_com

!------------------------------------------------------------------------------

      subroutine def_rsf_obio(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      USE obio_forc, only : avgq,tirrq3d
      USE obio_com, only : gcmax,nstep0,pp2tot_day, arg2d, arg3d
      USE HYCOM_DIM, only : grid=>ogrid
      use pario, only : defvar
      USE obio_com, only : pCO2av=>pCO2av_loc,
     &    ao_co2fluxav=>ao_co2fluxav_loc,
     &    pp2tot_dayav=>pp2tot_dayav_loc,
     &    cexpav=>cexpav_loc, diag_counter,
     .    pp2diat_dayav_loc,pp2chlo_dayav_loc,pp2cyan_dayav_loc,
     .    pp2cocc_dayav_loc,pHav_loc
      
      implicit none
      integer, intent(in) :: fid   !@var fid file id

      call defvar(grid,fid,nstep0,'obio_nstep0')
      call defvar(grid,fid,avgq,'avgq('//trim(arg3d)//')')
      call defvar(grid,fid,gcmax,'gcmax('//trim(arg3d)//')')
      call defvar(grid,fid,tirrq3d,'tirrq3d('//trim(arg3d)//')')
      call defvar(grid,fid,pp2tot_day,'pp2tot_day('//trim(arg2d)//')')
      call defvar(grid,fid,diag_counter,'obio_diag_counter')
      call defvar(grid,fid,pCO2av,'pCO2av('//trim(arg2d)//')')
      call defvar(grid,fid,pp2tot_dayav,'pp2tot_dayav('//
     &                                                 trim(arg2d)//')')
      call defvar(grid,fid,ao_co2fluxav,'ao_co2fluxav('//
     &                                                 trim(arg2d)//')')
      call defvar(grid,fid,cexpav,'cexpav('//trim(arg2d)//')')
      call defvar(grid,fid,pp2diat_dayav_loc,'pp2diat_day('//
     &                                                 trim(arg2d)//')')
      call defvar(grid,fid,pp2chlo_dayav_loc,'pp2chlo_day('//
     &                                                 trim(arg2d)//')')
      call defvar(grid,fid,pp2cyan_dayav_loc,'pp2cyan_day('//
     &                                                 trim(arg2d)//')')
      call defvar(grid,fid,pp2cocc_dayav_loc,'pp2cocc_day('//
     &                                                 trim(arg2d)//')')
      call defvar(grid,fid,phav_loc,'pHav('//trim(arg2d)//')')
      
      return
      end subroutine def_rsf_obio


      subroutine new_io_obio(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      USE obio_forc, only : avgq,tirrq3d
      USE obio_com, only : gcmax,nstep0
     &     ,pp2tot_day
      USE HYCOM_DIM, only : grid=>ogrid
      USE obio_com, only : pCO2av=>pCO2av_loc,
     &     ao_co2fluxav=>ao_co2fluxav_loc,
     &     pp2tot_dayav=>pp2tot_dayav_loc,
     &     cexpav=>cexpav_loc, diag_counter
     .    ,pp2diat_dayav_loc,pp2chlo_dayav_loc,pp2cyan_dayav_loc
     .    ,pp2cocc_dayav_loc,pHav_loc
      
      implicit none

      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_data(grid,fid,'obio_nstep0',nstep0)
        call write_dist_data(grid,fid,'avgq',avgq)
        call write_dist_data(grid,fid,'gcmax',gcmax)
        call write_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call write_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
        call write_data(grid,fid,'obio_diag_counter',diag_counter)
        call write_dist_data(grid,fid,'pCO2av',pCO2av)
        call write_dist_data(grid,fid,'pp2tot_dayav',pp2tot_dayav)
        call write_dist_data(grid,fid,'ao_co2fluxav',ao_co2fluxav)
        call write_dist_data(grid,fid,'cexpav',cexpav)
        call write_dist_data(grid,fid,'pp2diat_day',pp2diat_dayav_loc)
        call write_dist_data(grid,fid,'pp2chlo_day',pp2chlo_dayav_loc)
        call write_dist_data(grid,fid,'pp2cyan_day',pp2cyan_dayav_loc)
        call write_dist_data(grid,fid,'pp2cocc_day',pp2cocc_dayav_loc)
        call write_dist_data(grid,fid,'pHav',phav_loc)
      case (ioread)            ! input from restart file
        call read_data(grid,fid,'obio_nstep0',nstep0,
     &       bcast_all=.true.)
        call read_dist_data(grid,fid,'avgq',avgq)
        call read_dist_data(grid,fid,'gcmax',gcmax)
        call read_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call read_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
        call read_data(grid,fid,'obio_diag_counter',diag_counter)
        call read_dist_data(grid,fid,'pCO2av',pCO2av)
        call read_dist_data(grid,fid,'pp2tot_dayav',pp2tot_dayav)
        call read_dist_data(grid,fid,'ao_co2fluxav',ao_co2fluxav)
        call read_dist_data(grid,fid,'cexpav',cexpav)
        call read_dist_data(grid,fid,'pp2diat_day',pp2diat_dayav_loc)
        call read_dist_data(grid,fid,'pp2chlo_day',pp2chlo_dayav_loc)
        call read_dist_data(grid,fid,'pp2cyan_day',pp2cyan_dayav_loc)
        call read_dist_data(grid,fid,'pp2cocc_day',pp2cocc_dayav_loc)
        call read_dist_data(grid,fid,'pHav',phav_loc)
      end select
      return
      end subroutine new_io_obio


      subroutine obio_set_data_after_archiv
      USE obio_com, only:
     .     diag_counter
     .    ,ao_co2fluxav_loc, pp2tot_dayav_loc
     .    ,pCO2av_loc, cexpav_loc
     .    ,pp2diat_dayav_loc,pp2chlo_dayav_loc,pp2cyan_dayav_loc
     .    ,pp2cocc_dayav_loc,pHav_loc
#ifdef TRACERS_Alkalinity
     .    ,caexpav_loc
#endif
      implicit none
      diag_counter=0
      ao_co2fluxav_loc=0
      pCO2av_loc = 0
      pp2tot_dayav_loc = 0
      cexpav_loc = 0
#ifdef TRACERS_Alkalinity
      caexpav_loc = 0
#endif
      pp2diat_dayav_loc=0
      pp2chlo_dayav_loc=0
      pp2cyan_dayav_loc=0
      pp2cocc_dayav_loc=0
      pHav_loc=0
      end subroutine obio_set_data_after_archiv

      subroutine obio_gather_before_archive
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
      USE obio_com, only:
     .     ao_co2fluxav_loc, ao_co2fluxav
     .    ,pCO2av_loc, pCO2av
     .    ,pp2tot_dayav_loc, pp2tot_dayav
     .    ,cexpav_loc, cexpav
     .    ,pp2diat_dayav_loc,pp2chlo_dayav_loc,pp2cyan_dayav_loc
     .    ,pp2cocc_dayav_loc,pHav_loc
     .    ,pp2diat_dayav,pp2chlo_dayav,pp2cyan_dayav
     .    ,pp2cocc_dayav,pHav
#ifdef TRACERS_Alkalinity
     .    ,caexpav_loc, caexpav
#endif
      implicit none 

      call pack_data(ogrid, ao_co2fluxav_loc, ao_co2fluxav)
      call pack_data(ogrid, pCO2av_loc, pCO2av)
      call pack_data(ogrid, pp2tot_dayav_loc, pp2tot_dayav)
      call pack_data(ogrid, cexpav_loc, cexpav)
#ifdef TRACERS_Alkalinity
      call pack_data(ogrid, caexpav_loc, caexpav)
#endif
      call pack_data(ogrid, pp2diat_dayav_loc, pp2diat_dayav)
      call pack_data(ogrid, pp2chlo_dayav_loc, pp2chlo_dayav)
      call pack_data(ogrid, pp2cyan_dayav_loc, pp2cyan_dayav)
      call pack_data(ogrid, pp2cocc_dayav_loc, pp2cocc_dayav)
      call pack_data(ogrid, pHav_loc, pHav)
      return
      end subroutine obio_gather_before_archive

