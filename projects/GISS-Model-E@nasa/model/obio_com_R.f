#include "rundeck_opts.h"
      MODULE obio_com
!@sum  obio_com contains the parameters, arrays and definitions
!@+    necessary for the OceanBiology routines
!@auth NR

      USE obio_dim
      use ocalbedo_mod, only: nlt
!     USE Constant, only: sday     ! sday=86400.0    !seconds per day

      USE OCEANRES, only : kdm=>lmo
      use ocean, only : jm,use_qus

      implicit none

c --- dobio       activate Watson Gregg's ocean biology code
      logical dobio
      data dobio/.true./
c

      real, ALLOCATABLE, DIMENSION(:,:)    :: tzoo2d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tfac3d,wshc3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: Fescav3d

#ifdef restart_add_o2
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: o2rst
#endif
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

      real, ALLOCATABLE, DIMENSION(:,:,:,:):: tracer
      real, ALLOCATABLE, DIMENSION(:,:)::  Edz,Euz,Esz
      real, ALLOCATABLE, DIMENSION(:,:)::  Kd,Kd_em2d !absorption+scattering in seawater due to chl
      real, ALLOCATABLE, DIMENSION(:)::  Kpar     !kpar from NBOM
      real, ALLOCATABLE, DIMENSION(:)::  Kpar_em2d !kpar from NBOM in Einstein/m2/day
      real, ALLOCATABLE, DIMENSION(:)::  delta_temp1d  !change in T due to kpar
      real, ALLOCATABLE, DIMENSION(:)::  obio_lambdas  !wavelengths in water column

      integer :: nstep0=0

      integer:: num_tracers
      integer:: errchk1=0
      integer:: errchk2=0    !@PL check if there are NaNs in O2     

      !test point
!!    integer, parameter :: itest=16, jtest=45    !equatorial Pacific                  2deg ocean
!!    integer, parameter :: itest=32, jtest=20    !southern ocean; Pacific          
!!    integer, parameter :: itest=1,  jtest=jm/2     !equator Pacific
#ifdef OBIO_QUIET_MODE
      integer, parameter :: itest=2,  jtest=1     ! point with no ocean
#else
      integer, parameter :: itest=272,  jtest=22     !regr test point
#endif

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
#ifdef TRACERS_degC
      real cexpdeg
#endif
      real temp1d(kdm),dp1d(kdm),obio_P(kdm,ntyp)
     .                 ,det(kdm,ndet),car(kdm,ncar),avgq1d(kdm)
     .                 ,gcmax1d(kdm),saln1d(kdm),p1d(kdm+1)
     .                 ,alk1d(kdm),flimit(kdm,nchl,5),rho1d(kdm)
!PLdbg
#ifdef TRACERS_Ocean_O2
     .                 ,o21d(kdm)   ! oxygen 1d array
     .                 ,abo21d(kdm) ! abiotic oxygen 1d array
#endif
#ifdef TRACERS_degC
     .                 ,ndegC1d(kdm) !@PL nondegradable carbon 1d array
#endif
      real rho_water
      real atmFe_ij,covice_ij

      real, allocatable, DIMENSION(:,:) :: daily_atmFe
      integer inwst,inwnd,jnwst,jnwnd     !starting and ending indices 
                                          !for daylight 
                                          !in i and j directions
      real acdom(kdm,nlt)                 !absorption coefficient of CDOM
      real P_tend(kdm,ntyp)               !bio tendency (dP/dt)

#ifdef TRACERS_Alkalinity
      real A_tend(kdm), co3_conc
      real ca_det_calc1d(kdm),Ca_tend(kdm)
#endif
!PLdbg
#ifdef TRACERS_Ocean_O2
      real O_tend(kdm),Abo_tend(kdm) ! oxygen tendency terms
#endif
      real rmuplsr(kdm,nchl)                  !growth+resp 
      real D_tend(kdm,ndet)                   !detrtial tendency
#ifdef TRACERS_degC
      real Ndeg_tend(kdm)                     !nondegradable C tendency
#endif 
      real obio_ws(kdm+1,nchl)                !phyto sinking rate
      real tfac(kdm)                          !phyto T-dependence
      real pnoice(kdm)                        !pct ice-free
      real wsdet(kdm+1,ndet)                  !detrital sinking rate
      real rikd(kdm,nchl)                     !photoadaption state
      real tzoo                               !herbivore T-dependence
      real docbac(kdm)                             !bacterial loss of DOC
      real rmu3(kdm,nchl)                      !grwoth on nitrate
      real rmu4(kdm,nchl)                      !growth on ammonium
      real gronfix(kdm)                       !growth from nitrogen fixation
      real dicresp(nchl)                      !phyto respiration 
      real Fescav(kdm)                        !iron scavenging rate

C if NCHL_DEFINED > 3
      real wshc(kdm)                          !cocco sinking rate
C endif

      real :: C_tend(kdm,ncar)                !carbon tendency
      real :: pCO2_ij,pHsfc                   !partial pressure of CO2, pH
!PLdbg
#ifdef TRACERS_Ocean_O2
      real :: pO2_ij,pabO2_ij                          !partial presure O2
#endif
      real :: gro(kdm,nchl)                   !realized growth rate
      integer :: day_of_month, hour_of_day

      real :: rhs(kdm,ntrac,17)         !secord arg-refers to tracer definition 

      real :: pp2_1d(kdm,nchl)          !net primary production

      real*8 :: co2flux
!PLdbg
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
     &          ,o2flux
#endif
#ifdef TRACERS_abio_O2 
     &          ,abo2flux ! Air-sea o2 fluxes
#endif
#endif
      integer kzc
      real*8 :: carb_old,iron_old    !prev timesetep total carbon inventory

      
      real*8 trmo_unit_factor(kdm,ntrac)

#ifdef restoreIRON
!this is an AR5 preprocessor option
!per year change dI/I, + for sink, - for source
      !!!real*8 :: Iron_BC = 0.002
      real*8 :: Iron_BC = -0.005
#endif

      real*8, dimension(:, :, :), allocatable :: ze

#ifdef OBIO_RUNOFF
!     real, ALLOCATABLE, DIMENSION(:,:)    :: rnitrmflo_loc      ! riverine nitrate mass flow rate (kg/s)
      real, ALLOCATABLE, DIMENSION(:,:)    :: rnitrconc_loc      ! riverine nitrate concentration (kg/kg)
      real, ALLOCATABLE, DIMENSION(:,:)    :: rdicconc_loc       ! riverine dic concentration (kg/kg)
      real, ALLOCATABLE, DIMENSION(:,:)    :: rdocconc_loc       ! riverine doc concentration (kg/kg)
      real, ALLOCATABLE, DIMENSION(:,:)    :: rsiliconc_loc      ! riverine silica concentration (kg/kg)
      real, ALLOCATABLE, DIMENSION(:,:)    :: rironconc_loc      ! riverine iron concentration (kg/kg)
      real, ALLOCATABLE, DIMENSION(:,:)    :: rpocconc_loc       ! riverine poc concentration (kg/kg)
      real, ALLOCATABLE, DIMENSION(:,:)    :: ralkconc_loc       ! riverine alkalinity concentration (mol/kg)

      real rnitrconc
!   .    , rnitrmflo
      real rdicconc
      real rdocconc
      real rsiliconc
      real rironconc
      real rpocconc
      real ralkconc
#endif  /* obio_runoff */



      character(len=50) :: arg2d, arg3d

      contains


      subroutine build_ze
      
      use oceanr_dim, only: ogrid
      use oceanres, only: kdm=>lmo
      use ofluxes, only: oapress
      use ocean, only: g0m,s0m,mo,dxypo,lmm
      use constant, only: grav
      
      implicit none
      integer :: k
      real :: pres,g,s
      integer :: i,j
      real*8 :: volgsp

      if (.not.allocated(ze)) allocate(ze(ogrid%i_strt:ogrid%i_stop,
     &                                ogrid%j_strt:ogrid%j_stop, 0:kdm))
      ze=0
      
      do i=ogrid%i_strt,ogrid%i_stop
        do j=ogrid%j_strt,ogrid%j_stop
          pres=oapress(i,j) 
          do k=1, lmm(i,j)
            pres=pres+mo(i,j,k)*grav*.5
            g=g0m(i,j,k)/(mo(i,j,k)*dxypo(j))
            s=s0m(i,j,k)/(mo(i,j,k)*dxypo(j))
            ze(i,j,k)=ze(i,j,k-1)+mo(i,j,k)*volgsp(g,s,pres)
            pres=pres+mo(i,j,k)*grav*.5
          end do
        end do
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
      use obio_dim, only: ntrac,nchl,ndet
      implicit none
      private

      public :: init_obio_diag, new_io_obio_diag, def_meta_obio_diag,
     &   write_meta_obio_diag, def_rsf_obio_diag, reset_obio_diag,
     &   add_diag

      real*8, dimension(:, :, :), allocatable, public :: obio_ij
!      real*8, dimension(:, :, :), allocatable, public :: hemis_obio_ij
      real*8, dimension(:, :, :, :), allocatable, public :: obio_ijl
#ifdef obio_rhsdiags
      real*8, dimension(:, :, :, :), allocatable, public :: rhs_ijl
#endif
      integer, public :: ij_solz, ij_sunz, ij_dayl, ij_ed, ij_es,
     &   ij_nitr, ij_amm, ij_sil, ij_iron, ij_diat, ij_chlo, ij_cyan,
     &   ij_cocc, ij_herb, ij_doc, ij_dic, ij_pco2, ij_alk, ij_flux,
     &   ij_cexp, ij_ndet, ij_setl, ij_sink, ij_xchl, ij_fca, 
     &   ij_rnitrmflo,
     &   ij_rnitrconc, ij_rdicconc, ij_rdocconc, ij_rsiliconc,
     &   ij_rironconc, ij_rpocconc, ij_ralkconc, ij_pp, ij_lim(4,5),
     &   ij_rhs(ntrac,17),ij_pp1, ij_pp2, ij_pp3, ij_pp4, ij_co3,
     &   ij_ph,kobio_ij
!PLdbg
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
     &  ,ij_o2,ij_oflx,ij_po2
#endif
#ifdef TRACERS_abio_O2
     &  ,ij_abo2,ij_aboflx,ij_pabo2
#endif
#endif
#ifdef TRACERS_degC
     &  ,ij_degC,ij_cexpdeg
#endif
      integer, public :: ijl_avgq, ijl_kpar,ijl_kpar_em2d,ijl_dtemp
     .                  ,ijl_wss(nchl),ijl_wsdet(ndet)
     .                  ,ijl_pp(nchl)
     .                  ,ijl_lim1(nchl)
     .                  ,ijl_lim2(nchl)
     .                  ,ijl_lim3(nchl)
     .                  ,ijl_lim4(nchl)
     .                  ,ijl_lim5(nchl)
     .                  ,ijl_rhs3(ntrac,17)
!PLdbg
#ifdef TRACERS_bio_O2
     .                  ,ijl_cprod
     .                  ,ijl_cdet
     .                  ,ijl_cdoc
     .                  ,ijl_caresp
     .                  ,ijl_chresp
     .                  ,ijl_oprodam
     .                  ,ijl_oprodnit
     .                  ,ijl_odet
     .                  ,ijl_odoc
     .                  ,ijl_oaresp
     .                  ,ijl_ohresp
#endif
      type(vector_str30) :: sname_ij, units_ij
      type(vector_str30) :: sname_ijl, units_ijl
      type(vector_str80) :: lname_ij, lname_ijl
      type(vector_integer) :: ia_ij, ia_ijl
      type(vector_integer) :: denom_ij
      type(vector_real8) :: scale_ij, scale_ijl
      type(cdl_type) :: cdl_ij, cdl_ijl
#ifdef obio_rhsdiags
      type(vector_str30) :: sname_rhs_ijl, units_rhs_ijl
      type(vector_str80) :: lname_rhs_ijl
      type(vector_integer) :: ia_rhs_ijl
      type(vector_real8) :: scale_rhs_ijl
      type(cdl_type) :: cdl_rhs_ijl
#endif
      type(cdl_type), pointer :: cdl_lons, cdl_lats, cdl_depths
      
      contains

      subroutine add_diag(lname, sname, units, dim3, idx)

      use mdiag_com, only : ia_cpl
      use DIAG_COM, only : IJ_POCEAN 

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
#ifdef obio_rhsdiags
      if (sname1(5:7).eq.'rhs') then
        call lname_rhs_ijl%push_back(lname1)
        call sname_rhs_ijl%push_back(sname1)
        call units_rhs_ijl%push_back(units1)
        call scale_rhs_ijl%push_back(1.d0)
        call ia_rhs_ijl%push_back(ia_cpl)
        idx=lname_rhs_ijl%getsize()
      else
#endif
        call lname_ijl%push_back(lname1)
        call sname_ijl%push_back(sname1)
        call units_ijl%push_back(units1)
        call scale_ijl%push_back(1.d0)
        call ia_ijl%push_back(ia_cpl)
        idx=lname_ijl%getsize()
#ifdef obio_rhsdiags
      endif
#endif
      else
        call lname_ij%push_back(lname1)
        call sname_ij%push_back(sname1)
        call units_ij%push_back(units1)
        call scale_ij%push_back(1.d0)
        call denom_ij%push_back(IJ_POCEAN)
        call ia_ij%push_back(ia_cpl)
        idx=lname_ij%getsize()
      endif
      end subroutine add_diag


      subroutine init_obio_diag
      
      USE OCEANR_DIM, only: ogrid
      USE OCEANRES, only: kdm=>lmo
      use odiag, only: cdl_olons, cdl_olats, cdl_odepths
      use obio_com, only: arg2d, arg3d
      implicit none

      cdl_lons=>cdl_olons
      cdl_lats=>cdl_olats
      cdl_depths=>cdl_odepths
      arg2d='dist_imo,dist_jmo'
      arg3d='dist_imo,dist_jmo,lmo'
      
      allocate(obio_ij(ogrid%i_strt:ogrid%i_stop,
     &         ogrid%j_strt:ogrid%j_stop, lname_ij%getsize()))
      allocate(obio_ijl(ogrid%i_strt:ogrid%i_stop,
     &         ogrid%j_strt:ogrid%j_stop, kdm, lname_ijl%getsize()))
!@PL
!      allocate(hemis_obio_ij(1,3,lname_ij%getsize()))
!@PL

#ifdef obio_rhsdiags
      allocate(rhs_ijl(ogrid%i_strt:ogrid%i_stop,
     &         ogrid%j_strt:ogrid%j_stop, kdm, lname_rhs_ijl%getsize()))
#endif

      end subroutine init_obio_diag


      subroutine def_rsf_obio_diag(fid, r4_on_disk)
      
      use pario, only: defvar
      use obio_com, only: arg2d, arg3d
      USE OCEANR_DIM, only: ogrid
      
      implicit none

      integer, intent(in) :: fid
      logical, intent(in) :: r4_on_disk

      call defvar(ogrid, fid, obio_ij,
     &   'obio_ij('//trim(arg2d)//',kobio_ij)', r4_on_disk=r4_on_disk)
      call defvar(ogrid, fid, obio_ijl,
     &   'obio_ijl('//trim(arg3d)//',kobio_ijl)', r4_on_disk=r4_on_disk)
!@PL
!      call defvar(ogrid,fid,hemis_obio_ij,
!     &      'hemis_obio_ij(one,shnhgm,kobio_ij)', r4_on_disk=r4_on_disk)
!@PL
#ifdef obio_rhsdiags
      call defvar(ogrid, fid, rhs_ijl,
     &   'rhs_ijl('//trim(arg3d)//',krhs_ijl)', r4_on_disk=r4_on_disk)
#endif

      end subroutine def_rsf_obio_diag


      subroutine new_io_obio_diag(fid, iaction)
     
      use model_com, only: ioread
      use pario, only: write_dist_data, read_dist_data
      USE OCEANR_DIM, only: ogrid
      
      implicit none

      integer, intent(in) :: fid
      integer, intent(in) :: iaction

      if (iaction.eq.ioread) then
        call read_dist_data(ogrid, fid, 'obio_ij', obio_ij)
        call read_dist_data(ogrid, fid, 'obio_ijl', obio_ijl)
#ifdef obio_rhsdiags
        call read_dist_data(ogrid, fid, 'rhs_ijl', rhs_ijl)
#endif
      else
        call write_dist_data(ogrid, fid, 'obio_ij', obio_ij)
        call write_dist_data(ogrid, fid, 'obio_ijl', obio_ijl)
#ifdef obio_rhsdiags
        call write_dist_data(ogrid, fid, 'rhs_ijl', rhs_ijl)
#endif
      end if

      end subroutine new_io_obio_diag


      subroutine def_meta_obio_diag(fid)
     
      use pario, only: defvar, write_attr
      use cdl_mod, only: defvar_cdl, merge_cdl, add_var
      use domain_decomp_1d, only: am_i_root
      USE OCEANR_DIM, only: ogrid
      
      implicit none

      integer, intent(in) :: fid
      integer :: k

      if (associated(cdl_lons)) then
        if (am_i_root()) then
          call merge_cdl(cdl_lons, cdl_lats, cdl_ij)
          call merge_cdl(cdl_ij, cdl_depths, cdl_ijl)
#ifdef obio_rhsdiags
          call merge_cdl(cdl_ij, cdl_depths, cdl_rhs_ijl)
#endif
          
        do k=1, sname_ij%getsize()
            call add_var(cdl_ij,
     &         'float '//trim(sname_ij%at(k))//'(lato,lono) ;',
!     &       auxvar_string=                                            !@PL aux var for hemis
!     &           'float '//trim(sname_ij%at(k))//'_hemis(shnhgm);',    !@PL
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
#ifdef obio_rhsdiags
          do k=1, sname_rhs_ijl%getsize()
            call add_var(cdl_rhs_ijl,
     &         'float '//trim(sname_rhs_ijl%at(k))//'(zoc,lato,lono) ;',
     &         long_name=trim(lname_rhs_ijl%at(k)),
     &         units=trim(units_rhs_ijl%at(k)),
     &         set_miss=.true.)
          enddo
#endif
        endif
        call defvar_cdl(ogrid, fid, cdl_ij,
     &                  'cdl_obio_ij(cdl_strlen,kcdl_obio_ij)')
        call defvar_cdl(ogrid, fid, cdl_ijl,
     &       'cdl_obio_ijl(cdl_strlen,kcdl_obio_ijl)')
#ifdef obio_rhsdiags
        call defvar_cdl(ogrid, fid, cdl_rhs_ijl,
     &       'cdl_rhs_ijl(cdl_strlen,kcdl_rhs_ijl)')
#endif
      endif

      call write_attr(ogrid, fid, 'obio_ij', 'reduction', 'sum')
      call write_attr(ogrid, fid, 'obio_ij', 'split_dim', 3)

!@PL
!      call write_attr(ogrid,fid,'hemis_obio_ij','reduction','sum')
!@PL
      call defvar(ogrid, fid, ia_ij%getdata(), 'ia_obio_ij(kobio_ij)')
      call defvar(ogrid, fid, scale_ij%getdata(),
     &            'scale_obio_ij(kobio_ij)')
      call defvar(ogrid, fid, sname_ij%getdata(),
     &            'sname_obio_ij(sname_strlen,kobio_ij)')


!@PL
      call defvar(ogrid, fid, denom_ij%getdata(),
     &            'denom_obio_ij(kobio_ij)')
!@PL

      call write_attr(ogrid, fid, 'obio_ijl', 'reduction', 'sum')
      call write_attr(ogrid, fid, 'obio_ijl', 'split_dim', 4)
      call defvar(ogrid, fid, ia_ijl%getdata(),
     &            'ia_obio_ijl(kobio_ijl)')
      call defvar(ogrid, fid, scale_ijl%getdata(),
     &            'scale_obio_ijl(kobio_ijl)')
      call defvar(ogrid, fid, sname_ijl%getdata(),
     &            'sname_obio_ijl(sname_strlen,kobio_ijl)')

#ifdef obio_rhsdiags
      call write_attr(ogrid, fid, 'rhs_ijl', 'reduction', 'sum')
      call write_attr(ogrid, fid, 'rhs_ijl', 'split_dim', 4)
      call defvar(ogrid, fid, ia_rhs_ijl%getdata(),
     &            'ia_rhs_ijl(krhs_ijl)')
      call defvar(ogrid, fid, scale_rhs_ijl%getdata(),
     &            'scale_rhs_ijl(krhs_ijl)')
      call defvar(ogrid, fid, sname_rhs_ijl%getdata(),
     &            'sname_rhs_ijl(sname_strlen,krhs_ijl)')
#endif

      !@PL size of obio_ij
      kobio_ij = sname_ij%getsize()
      !@PL

      end subroutine def_meta_obio_diag


      subroutine write_meta_obio_diag(fid)
      use pario, only: write_data, write_dist_data
      use cdl_mod, only: write_cdl
      USE OCEANR_DIM, only: ogrid
      implicit none

      integer, intent(in) :: fid

      call write_data(ogrid, fid, 'ia_obio_ij', ia_ij%getdata())
      call write_data(ogrid, fid, 'scale_obio_ij', scale_ij%getdata())
      call write_data(ogrid, fid, 'sname_obio_ij', sname_ij%getdata())
      !@PL
      call write_data(ogrid, fid, 'denom_obio_ij', denom_ij%getdata())
!      call write_data(ogrid,fid,'hemis_obio_ij',hemis_obio_ij)

      !@PL
      if (associated(cdl_lons))
     &      call write_cdl(ogrid, fid, 'cdl_obio_ij', cdl_ij)

      call write_data(ogrid, fid, 'ia_obio_ijl', ia_ijl%getdata())
      call write_data(ogrid, fid, 'scale_obio_ijl', scale_ijl%getdata())
      call write_data(ogrid, fid, 'sname_obio_ijl', sname_ijl%getdata())
      if (associated(cdl_lons))
     &      call write_cdl(ogrid, fid, 'cdl_obio_ijl', cdl_ijl)

#ifdef obio_rhsdiags
      call write_data(ogrid, fid,
     .                       'ia_rhs_ijl', ia_rhs_ijl%getdata())
      call write_data(ogrid, fid,
     .                       'scale_rhs_ijl', scale_rhs_ijl%getdata())
      call write_data(ogrid, fid,
     .                       'sname_rhs_ijl', sname_rhs_ijl%getdata())
      if (associated(cdl_lons))
     &      call write_cdl(ogrid, fid, 'cdl_rhs_ijl', cdl_rhs_ijl)
#endif

      end subroutine write_meta_obio_diag


      subroutine reset_obio_diag
      implicit none

      obio_ij=0.
      obio_ijl=0.
#ifdef obio_rhsdiags
      rhs_ijl=0.
#endif

      end subroutine reset_obio_diag


      end module obio_diag



!------------------------------------------------------------------------------
      subroutine alloc_obio_com

      USE obio_com
      USE obio_dim
      use obio_diag, only: init_obio_diag

      USE OCEANR_DIM, only : ogrid
      USE OCEANRES, only :kdm=>lmo, idm=>imo,jdm=>jmo

      implicit none

c**** Extract domain decomposition info
      INTEGER :: j_0,j_1,i_0,i_1

      I_0 = ogrid%I_STRT
      I_1 = ogrid%I_STOP
      J_0 = ogrid%J_STRT
      J_1 = ogrid%J_STOP


      ALLOCATE(tracer(i_0:i_1,j_0:j_1,kdm,ntrac))

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

#ifdef restart_add_o2
      ALLOCATE(o2rst(i_0:i_1,j_0:j_1,kdm))
#endif

      ALLOCATE(Edz(nlt,kdm))
      ALLOCATE(Esz(nlt,kdm))
      ALLOCATE(Euz(nlt,kdm))
      ALLOCATE(Kd(nlt,kdm))
      ALLOCATE(Kd_em2d(nlt,kdm))
      ALLOCATE(Kpar(kdm))
      ALLOCATE(Kpar_em2d(kdm))
      ALLOCATE(delta_temp1d(kdm))

#ifdef OBIO_RUNOFF
!     ALLOCATE(rnitrmflo_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(rnitrconc_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(rdicconc_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(rdocconc_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(rsiliconc_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(rironconc_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(rpocconc_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(ralkconc_loc(i_0:i_1,j_0:j_1))
#endif

      call init_obio_diag

      end subroutine alloc_obio_com

!------------------------------------------------------------------------------

      subroutine def_rsf_obio(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      USE obio_forc, only : avgq,tirrq3d
      USE obio_com, only : gcmax,nstep0,pp2tot_day, arg2d, arg3d
      USE OCEANR_DIM, only : grid=>ogrid
      use pario, only : defvar
      implicit none
      integer, intent(in) :: fid   !@var fid file id

      call defvar(grid,fid,nstep0,'obio_nstep0')
      call defvar(grid,fid,avgq,'avgq('//trim(arg3d)//')')
      call defvar(grid,fid,gcmax,'gcmax('//trim(arg3d)//')')
      call defvar(grid,fid,tirrq3d,'tirrq3d('//trim(arg3d)//')')
      call defvar(grid,fid,pp2tot_day,'pp2tot_day('//trim(arg2d)//')')
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
      USE OCEANR_DIM, only : grid=>ogrid
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
      case (ioread)            ! input from restart file
        call read_data(grid,fid,'obio_nstep0',nstep0,
     &       bcast_all=.true.)
        call read_dist_data(grid,fid,'avgq',avgq)
        call read_dist_data(grid,fid,'gcmax',gcmax)
        call read_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call read_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
      end select
      return
      end subroutine new_io_obio

!---------------------------------------------------------------------------

      subroutine setup_obio

      use ocn_tracer_com, only: add_ocn_tracer
      use runtimecontrols_mod, only: tracers_alkalinity
      use obio_dim, only: ntrac,nchl,ndet,nnut,ntyp,ndimc
!PLdbg
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
     .                   ,no2,ndimo2
#endif
#ifdef TRACERS_abio_O2
     .                   ,nabo2,ndimabo2
#endif
#endif
      use obio_diag

       implicit none

      integer, dimension(1) :: con_idx
      character(len=10), dimension(1) :: con_str
      integer :: nt, ilim, ll
      character :: str1*5,str2*9,str3*10,str4*6,str5*7,str6*9
      character :: unit_str*8
      character(len=1), parameter :: lim_sym(4)=(/'d', 'h', 'b', 'c'/)
! diatoms, chloroph, cyanobact, coccoliths
!@PL ifdef statements for O2 and alkalinity options
      character(len=4), parameter :: rhs_sym(ntrac)=(/ 'nitr', 'ammo',
     &     'sili', 'iron', 'diat', 'chlo', 'cyan', 'cocc', 'herb',
     &     'ndet', 'sdet', 'idet', 'doc_', 'dic_'
#ifdef TRACERS_Alkalinity
     &      ,'alk_'
#endif
!PLdbg
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
     &      ,'o2__'
#endif
#ifdef TRACERS_abio_O2
     &      ,'abo2'
#endif
#ifdef TRACERS_degC   
     &      ,'ndeg' !@PL nondegradable carbon (currently labeled as nitrogen, but is actually carbon)
#endif
#endif
     &           /)



!@PL
      con_idx=[12]
      con_str=['OCN BIOL']

      call add_ocn_tracer('Nitr      ', i_ntrocn=-4, i_ntrocn_delta=-12,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Ammo      ', i_ntrocn=-6, i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Sili      ', i_ntrocn=-4, i_ntrocn_delta=-12,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Iron      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Diat      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Chlo      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Cyan      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Cocc      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Herb      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('N_det     ', i_ntrocn=-6, i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('S_det     ', i_ntrocn=-6, i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
       call add_ocn_tracer('I_det     ',i_ntrocn=-10,i_ntrocn_delta=-18,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('DOC       ', i_ntrocn=-6, i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('DIC       ', i_ntrocn=-3, i_ntrocn_delta=-11,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      if (tracers_alkalinity)
     &  call add_ocn_tracer('Alk       ',i_ntrocn=-6,i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
!PLdbg
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
        call add_ocn_tracer('O2        ',i_ntrocn=-3,i_ntrocn_delta=-11,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
#endif
#ifdef TRACERS_abio_O2
        call add_ocn_tracer('abO2      ',i_ntrocn=-3,i_ntrocn_delta=-11,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
#endif
#endif
#ifdef TRACERS_degC
      call add_ocn_tracer('Nndeg     ', i_ntrocn=-6, i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
#endif
#ifdef TOPAZ_params
      call add_diag("co3 ", "oij_co3",
     &              "????????", .false., IJ_co3)
#endif
      call add_diag("ocean surface pH", "oij_pH",
     &              "pH units", .false., IJ_pH)
      call add_diag("Cos Solar Zenith Angle", "oij_solz",
     &              "NaN", .false., IJ_solz)
      call add_diag("Solar Zenith Angle", "oij_sunz",
     &              "degrees", .false., IJ_sunz)
      call add_diag("Daylight length", "oij_dayl",
     &              "timesteps", .false., IJ_dayl)
       call add_diag("Surface Ocean Direct Sunlight",
     &              "oij_Ed", "Wm-2", .false., IJ_Ed)
      call add_diag("Surface Ocean Diffuse Sunlight",
     &              "oij_Es", "Wm-2", .false., IJ_Es)
      call add_diag("Surface ocean Nitrates", "oij_nitr",
     &              "uM", .false., IJ_nitr)
      call add_diag("Surface ocean Ammonium", "oij_amm",
     &              "uM", .false., IJ_amm)
      call add_diag("Surface ocean Silicate", "oij_sil",
     &              "uM", .false., IJ_sil)
      call add_diag("Surface ocean Iron", "oij_iron",
     &              "nM", .false., IJ_iron)
      call add_diag("Surface ocean Diatoms", "oij_diat",
     &              "mg/m3", .false., IJ_diat)
      call add_diag("Surface ocean Chlorophytes", "oij_chlo",
     &              "mg/m3", .false., IJ_chlo)
      call add_diag("Surface ocean Cyanobacteria", "oij_cyan",
     &              "mg/m3", .false., IJ_cyan)
      call add_diag("Surface ocean Coccolithophores", "oij_cocc",
     &              "mg/m3", .false., IJ_cocc)
      call add_diag("Surface ocean Herbivores", "oij_herb",
     &              "mg/m3", .false., IJ_herb)
       call add_diag("Surface ocean DOC", "oij_doc",
     &              "uM", .false., IJ_doc)
      call add_diag("Surface ocean DIC", "oij_dic",
     &              "uM", .false., IJ_dic)
!PLdbg
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
      call add_diag("Surface ocean O2", "oij_o2",
     &              "mmol/kg", .false., IJ_o2)
      call add_diag("AO Flux O2",  "oij_oflx",
     &              "molO2/m2/yr", .false., IJ_oflx)
      call add_diag("Surface ocean partial O2 pressure",
     &              "oij_pO2","atm",.false.,IJ_po2)
#endif
#ifdef TRACERS_abio_O2
      call add_diag("Surface ocean abiotic O2", "oij_abo2",
     &              "mmol/kg", .false., IJ_abo2)
      call add_diag("AO Flux abO2",  "oij_aboflx",
     &              "molO2/m2/yr", .false., IJ_aboflx)
      call add_diag("Surface ocean abiotic partial O2 pressure",
     &              "oij_pabO2","atm",.false.,IJ_pabo2)
#endif
#endif
      call add_diag("Surface ocean partial CO2 pressure",
     &              "oij_pCO2", "uatm", .false., IJ_pCO2)
      call add_diag("Surface ocean alkalinity", "oij_alk",
     &              "umol/kg", .false., IJ_alk)
      call add_diag("AO Flux CO2 (gr,CO2 or mol,CO2/m2/yr)", "oij_flux",
     &              "depends if on atm/ocean grid", .false., IJ_flux)
      call add_diag("C export flux at compensation depth", "oij_cexp",
     &              "PgC/yr", .false., IJ_cexp)
      call add_diag("N/C detritus at 74m", "oij_ndet",
     &              "ugC/l", .false., IJ_ndet)
#ifdef TRACERS_degC
      call add_diag("deg C export flux at comp. depth", "oij_cexd", !degradable carbon flux
     &              "PgC/yr", .false., IJ_cexpdeg)
      call add_diag("N/C deg detritus at 74m", "oij_degC", ! degradable cabron at  74m
     &              "ugC/l", .false., IJ_degC)
#endif
      call add_diag("settlvel n/cdet at 74m", "oij_setl",
     &              "m/s", .false., IJ_setl)
      call add_diag("sink vel phytopl at 74m", "oij_sink",
     &              "m/s", .false., IJ_sink)
      call add_diag("C export due to chloroph", "oij_xchl",
     &              "kg,C*m/s", .false., IJ_xchl)
      if (tracers_alkalinity) then
        call add_diag("CaCO3 export flux at compensation depth",
     &                "oij_fca", "mili-g,C/m2/s", .false., IJ_fca) 
      endif

      call add_diag("Depth integrated PP", "oij_pp",
     &              "mg,C/m2/day", .false., IJ_pp)
      call add_diag("PP-diat", "oij_pp1",
     &              "mg,C/m2/day", .false., IJ_pp1)
      call add_diag("PP-chlor", "oij_pp2",
     &              "mg,C/m2/day", .false., IJ_pp2)
      call add_diag("PP-cyan", "oij_pp3",
     &              "mg,C/m2/day", .false., IJ_pp3)
      call add_diag("PP-cocc", "oij_pp4",
     &              "mg,C/m2/day", .false., IJ_pp4)
      do nt=1, 4
        do ilim=1, 5
          write(str1, '(A1,A3,I1)') lim_sym(nt), 'lim', ilim
          call add_diag(str1, str1, "?", .false., ij_lim(nt, ilim))
        end do
      end do
      do nt=1, ntrac
        do ll=1, 17
          write(str2, '(A4,A3,I2.2)') rhs_sym(nt), 'rhs', ll
          call add_diag(str2, str2, "?", .false., ij_rhs(nt, ll))
        end do
      end do
#ifdef OBIO_RUNOFF
!      call add_diag("Nitrate mass flow from rivers", "oij_rnitrmflo",
!     &               "kg/s", IJ_rnitrmflo)
      call add_diag("Nitrate conc in runoff", "oij_rnitrconc",
     &              "kg/kg", .false., IJ_rnitrconc)
      call add_diag("DIC conc in runoff", "oij_rdicconc",
     &              "kg/kg", .false., IJ_rdicconc)
      call add_diag("DOC conc in runoff", "oij_rdocconc",
     &              "kg/kg", .false., IJ_rdocconc)
      call add_diag("silica conc in runoff", "oij_rsiliconc",
     &              "kg/kg", .false., IJ_rsiliconc)
      call add_diag("iron conc in runoff", "oij_rironconc",
     &              "kg/kg", .false., IJ_rironconc)
      call add_diag("poc conc in runoff", "oij_rpocconc",
     &              "kg/kg", .false., IJ_rpocconc)
      call add_diag("alkalinity conc in runoff", "oij_ralkconc",
     &              "mol/kg", .false., IJ_ralkconc)
#endif


      call add_diag("Mean daily irradiance", "avgq",
     &              "quanta/m2/s", .true., IJL_avgq)

#ifdef KPAR_2_OCEAN
      call add_diag("KPAR", "kpar",
     &              "w/m2", .true., IJL_kpar)
      call add_diag("KPAR_EM2D", "kpar_em2d",
     &              "Einstein/m2 day", .true., IJL_kpar_em2d)
      call add_diag("dtemp due to kpar", "dtemp_par",
     &              "C", .true., IJL_dtemp)
#endif

      do nt=1, nchl
        write(str5,'(A4,A3)') rhs_sym(nnut+nt), 'wss'    !chl records 5:8
        call add_diag(str5, str5,"m/s", .true., IJL_wss(nt))
      enddo
      do nt=1, ndet
        write(str2,'(A4,A5)') rhs_sym(ntyp+nt), 'wsdet'    !detr records 10:12
        call add_diag(str2, str2,"m/s", .true., IJL_wsdet(nt))
      enddo
      do nt=1, nchl
        write(str4,'(A4,A2)') rhs_sym(nnut+nt), 'pp'
        call add_diag(str4, str4,"mg,C/m2/day", .true., IJL_pp(nt))
      enddo

!@PL rhs diagnostics when rhsobio is undefined
!PLdbg
#ifdef TRACERS_bio_O2
      call add_diag("O2 Nitrate production", "Oprod",
     &              "kg/d", .true., IJL_oprodnit)
      call add_diag("O2 auto respiration", "Oaresp",
     &              "kg/d", .true., IJL_oaresp)
      call add_diag("O2 heter respiration", "Ohresp",
     &              "kg/d", .true., IJL_ohresp)
      call add_diag("O2 detritial degradation", "Odet",
     &              "kg/d", .true., IJL_odet)
      call add_diag("O2 doc degradation", "Odoc",
     &              "kg/d", .true., IJL_odoc)
      call add_diag("O2 Ammonium production", "Oprodam",
     &              "kg/d", .true., IJL_oprodam)
      call add_diag("DIC Nitrate production", "Cprod",
     &              "kg/d", .true., IJL_cprod)
      call add_diag("DIC auto respiration", "Caresp",
     &              "kg/d", .true., IJL_caresp)
      call add_diag("DIC heter respiration", "Chresp",
     &              "kg/d", .true., IJL_chresp)
      call add_diag("DIC detritial degradation", "Cdet",
     &              "kg/d", .true., IJL_cdet)
      call add_diag("DIC doc degradation", "Cdoc",
     &              "kg/d", .true., IJL_cdoc)

#endif

      do nt=1,nchl
        write(str6,'(A4,A5)')rhs_sym(nnut+nt),'_lim1'
        call add_diag(str6,str6,"?",.true.,IJL_lim1(nt))

        write(str6,'(A4,A5)')rhs_sym(nnut+nt),'_lim2'
        call add_diag(str6,str6,"?",.true.,IJL_lim2(nt))

        write(str6,'(A4,A5)')rhs_sym(nnut+nt),'_lim3'
        call add_diag(str6,str6,"?",.true.,IJL_lim3(nt))

        write(str6,'(A4,A5)')rhs_sym(nnut+nt),'_lim4'
        call add_diag(str6,str6,"?",.true.,IJL_lim4(nt))

        write(str6,'(A4,A5)')rhs_sym(nnut+nt),'_lim5'
        call add_diag(str6,str6,"?",.true.,IJL_lim5(nt))
      enddo
!@PL currently including all tracers in rhs diagnostics causes an error
! in the checkpoint file. So only DIC and O2 are included for now when O2 
! tracer is included
#ifdef obio_rhsdiags
#ifdef TRACERS_Ocean_O2
      do nt=ndimc, ntrac
       if (nt.eq.ndimc
#ifdef TRACERS_bio_O2
     & .or. nt.eq.ndimo2
#endif
#ifdef TRACERS_abio_O2
     & .or. nt.eq.ndimabo2
#endif
     &           ) then
      do ll=1, 17
        write(str3, '(A4,A4,I2.2)') rhs_sym(nt), 'rhs3', ll
        if (nt.eq.ndimc) unit_str='kg,C/s'
!@PL
#ifdef TRACERS_bio_O2
        if (nt.eq.ndimo2) unit_str='kg,O/s'
#endif
#ifdef TRACERS_abio_O2
        if (nt.eq.ndimabo2) unit_str='kg,O/s'
#endif
        call add_diag(str3,str3,unit_str,.true.,IJL_rhs3(nt,ll))
!@PL
      enddo
      endif
      enddo
#else

      do nt=1, ntrac
      do ll=1, 17
        write(str3, '(A4,A4,I2.2)') rhs_sym(nt), 'rhs3', ll
        if (nt.eq.1) unit_str='kg,N/s'
        if (nt.eq.2) unit_str='kg,A/s'
        if (nt.eq.3) unit_str='kg,S/s'
        if (nt.eq.4) unit_str='kg,I/s'
        if (nt.ge.5.and.nt.le.9) unit_str='kg,C/s'
        if (nt.eq.10) unit_str='kg,C/s'
        if (nt.eq.11) unit_str='kg,S/s'
        if (nt.eq.12) unit_str='kg,I/s'
        if (nt.eq.13) unit_str='kg,C/s'
        if (nt.eq.ndimc-1) unit_str='kg,C/s'
        if (nt.eq.ndimc) unit_str='kg,C/s'
        call add_diag(str3,str3,unit_str,.true.,IJL_rhs3(nt,ll))
      enddo
      enddo
#endif
#endif

      return
      end subroutine setup_obio
     
      

