#include "rundeck_opts.h"
module CLOUDS_COM
!@sum  CLOUDS_COM model variables for moist convction and
!@+    large-scale condensation
!@+    cloud droplet number added to list of saves for rsf files
!@auth M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
  implicit none
  save
!@var TTOLD,QTOLD previous potential temperature, humidity
  real*8, allocatable, dimension(:,:,:) :: TTOLD,QTOLD
!@var SVLHX,SVLAT previous latent heat of evaporation
  real*8, allocatable, dimension(:,:,:) :: SVLHX,SVLAT
!@var RHSAV previous relative humidity
  real*8, allocatable, dimension(:,:,:) :: RHSAV
  !**** some arrays here for compatility with new clouds
!@var CLDSAV, CLDSAV1 previous cloud cover area (percent)
  real*8, allocatable, dimension(:,:,:) :: CLDSAV,CLDSAV1
!@var ULS,VLS,UMC,VMC velocity work arrays
  real*8, allocatable, dimension(:,:,:) :: ULS,VLS,UMC,VMC
!@var TLS,QLS,TMC,QMC temperature and humidity work arrays
  real*8, allocatable, dimension(:,:,:) :: TLS,QLS,TMC,QMC
!@var FSS grid fraction for large-scale clouds
!@+   initialised as 1. for compatibility with previous clouds
  real*8, allocatable, dimension(:,:,:) :: FSS
#ifdef mjo_subdd
!@var Source/sink from total, shallow,deep convection and large-scale condensation
  real*8, allocatable, dimension(:,:,:) :: TMCDRY,SMCDRY,DMCDRY, &
       LSCDRY
#endif
#ifdef etc_subdd
!@var liquid water path
  real*8, allocatable, dimension(:,:)   :: LWP2D
!@var ice water path
  real*8, allocatable, dimension(:,:)   :: IWP2D
#endif
#if (defined mjo_subdd) || (defined etc_subdd)
!@var Cloud liquid and ice water content for subdaily output
  real*8, allocatable, dimension(:,:,:) :: CLWC3D,CIWC3D
!@var Total, shallow, deep and large-scale latent heating for subdaily output
  real*8, allocatable, dimension(:,:,:) :: TLH3D,SLH3D,DLH3D,LLH3D
#endif
#ifdef CLD_AER_CDNC
!@var NCL old CDNC,NCI old ice crystal
  real*8, allocatable, dimension(:,:,:) :: NCL,NCI
!@var N, Re, LWP for 3 hrly diag save
  real*8, allocatable, dimension(:,:,:) :: CDN3D,CRE3D
  real*8, allocatable, dimension(:,:)   :: CLWP
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
!@var LWC,Cld depth, cld tem for 3 hrly diag save
  real*8, allocatable, dimension(:,:,:) :: CL3D,CI3D,CD3D,CTEM
#endif

  !**** variables saved for radiation calculations
!@var TAUSS optical depth from super-saturated clouds
!@var TAUSSIP counterpart to TAUSS for ice precip in stratiform (ss) liquid clouds
  real*8, allocatable, dimension(:,:,:) :: TAUSS,TAUSSIP
!@var TAUMC optical depth from moist-convective clouds
  real*8, allocatable, dimension(:,:,:) :: TAUMC
!@var CLDSS super-saturated cloud cover area (percent)
  real*8, allocatable, dimension(:,:,:) :: CLDSS
!@var CLDMC moist convective cloud cover area (percent)
  real*8, allocatable, dimension(:,:,:) :: CLDMC
!@var CSIZMC,CSIZSS mc,ss effective cloud droplet radius (microns)
!@var CSIZSSIP counterpart to CSIZSS for ice precip in stratiform liquid clouds
  real*8, allocatable, dimension(:,:,:) :: CSIZMC,CSIZSS,CSIZSSIP
!@var QLss,QIss stratiform liquid, ice water available to radiation (kg/kg)
  real*8, allocatable, dimension(:,:,:) :: QLss,QIss
!@var QLmc,QImc convective liquid, ice water available to radiation (kg/kg)
  real*8, allocatable, dimension(:,:,:) :: QLmc,QImc

  !**** variables saved for surface wind spectrum calculations
!@var DDM1 downdraft mass flux / rho at lowest level (m/s)
!@var DDML lowest level of downdraft
!@var DDMS downdraft mass flux at level 1 (kg/s/m**2)
!@var TDN1 downdraft temperature (K)
!@var QDN1 downdraft humidity (kg/kg)
  real*8, allocatable, dimension(:,:) :: DDM1,DDMS,TDN1,QDN1
  integer, allocatable, dimension(:,:) :: DDML

  !**** variables used (and saved) for gravity wave drag calculations
!@var AIRX, AIRMX*AREA convective mass flux (kg/s)
  real*8, allocatable, dimension(:,:) :: AIRX
!@var LMC max layer of mc convective mass flux.
  integer, allocatable, dimension(:,:,:) :: LMC

!@var LLOW,LMID,LHI max levels for low, mid and high clouds
  integer LLOW,LMID,LHI

!@var ISCCP_REG2D distributed version of ISCCP_REG
  integer, dimension(:,:), allocatable :: ISCCP_REG2D

!@var UKM,VKM arrays for vertical momentum mixing
  real*8, dimension(:,:,:,:), allocatable :: UKM,VKM
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
  real*8, allocatable, dimension(:,:)  :: NACTC      ! = 1.0D-30  ![#/m3](l,nmodes)
  real*8, allocatable, dimension(:,:)  :: NAERC      ! = 1.0D-30  ![#/m3](l,nmodes)
#endif
#endif

  !**** ISCCP diagnostics related parameter
#ifdef SCM
  integer,parameter :: ncol = 100    !@var ncol number of subcolumns
#else
  integer,parameter :: ncol = 20    !@var ncol number of subcolumns
#endif
#ifdef CACHED_SUBDD
#ifdef COSP_SIM
!@dbparam nsubdd_cosp: sub-daily diag freq for COSP simulator (set in rundeck)
!@var npoints_cosp: number of grids in a given spatial sub-domain for the COSP simulator
!@var flag_pfluxes_cosp trigger reflecting rundeck pre-processor flag COSP_PFLUX
!@var flag_re_cosp trigger reflecting rundeck pre-processor flag COSP_USERE
  integer :: nsubdd_cosp,npoints_cosp
  logical :: flag_pfluxes_cosp,flag_re_cosp
#endif
#endif

contains
  subroutine get_cld_overlap (lmax, cldssl, cldmcl, CldTot, CldSS, CldMC, RandSS)
!@sum Defines the cloud overlap type and either adjusts the random numbers
!@+    accordingly and/or finds the probability for the total cloud cover
!@auth R.Ruedy
  integer,          intent(in)    :: lmax
  real*8,           intent(in)    :: cldssl(lmax)
  real*8, optional, intent(in)    :: cldmcl(lmax)
  real*8, optional, intent(out)   :: CldTot
!@var CldSS stratiform cloud cover
  real*8, optional, intent(out)   :: CldSS
!@var CldMC convective cloud cover
  real*8, optional, intent(out)   :: CldMC
  real*8, optional, intent(inout) :: RandSS(lmax)
  integer L
  logical same_cloud
  real*8 clearmc, clearss, clearss_part

!  Current scheme: - semi-random overlap for supersaturation (SS) clouds
!                  - full overlap for moist convective (MC) clouds
!                  - random overlap of SS and MC clouds

!  SS Clouds: Clouds in consecutive layers are treated as part of the
!             same cloud with maximal overlap (single random number)
!             Random overlap is used for clouds separated by clear layers

!  Note 1: The first block is a slimmed-down version of the remaining part
!          (with the RandSS-line un-commented); the current MC and MC/SS
!          overlap schemes are already realized by having picked a single MC
!          random number for the whole column in addition to the ones for SS
!          (That block was added to avoid unneeded computations in long runs)
!  Note 2: Reverse looping over L was only kept for bit-wise consistency
!          with the previous version of the code.

   if(present(RandSS)) then
     same_cloud = .false.
     do L=lmax,1,-1       ! better:  1,lmax  (and replace L+1 by L-1 below)
        if( cldssl(L) > 0 ) then
          if(same_cloud) RandSS(L) = RandSS(L+1)        ! use same random #
          same_cloud = .true.
        else
          same_cloud = .false.
        end if
     end do
   end if
   if (.not.(present(CldTot).or.present(CldMC).or.present(CldSS))) return

   same_cloud = .false. ; clearss=1. ; clearss_part=1. ; clearmc=1.
   do L=1,lmax
      if( cldssl(L) > 0 ) then
!!      if(same_cloud) RandSS(L) = RandSS(L-1)            ! use same random #
        clearss_part = min( clearss_part, 1-cldssl(L) )   ! total overlap
        same_cloud = .true.
      else                                                ! clear sky layer
        same_cloud = .false.
        clearss = clearss * clearss_part                  ! random overlap
        clearss_part = 1.                                 ! reset for next cloud
      end if
   end do
   clearss = clearss * clearss_part
   if( present(CldSS) ) CldSS = 1.-clearss

!! Treat convective clouds as a single cloud with max. overlap
!! (using the same random number for the whole column)

   if( present(cldmcl) )then
     do l=1,lmax
        clearmc = min (clearmc, 1. - cldmcl(l))
     end do
     if (clearmc < 0) clearmc = 0.
     if( present(CldMC) ) CldMC = 1.-clearmc
   endif

!! Use random overlap for MC and SS cloud
   if( present(CldTot) ) CldTot = 1 - clearss*clearmc
   return

   end subroutine get_cld_overlap

end module CLOUDS_COM

subroutine ALLOC_CLOUDS_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
  use DOMAIN_DECOMP_ATM, only : DIST_GRID
  use RESOLUTION, only : LM
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
  use AERO_CONFIG, only: NMODES
#endif
#endif
  use CLOUDS_COM, only : TTOLD,QTOLD,SVLHX,SVLAT,RHSAV,CLDSAV,CLDSAV1,FSS
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
  use CLOUDS_COM, only : CL3D,CI3D,CD3D,CTEM
#endif
#ifdef CLD_AER_CDNC
  use CLOUDS_COM, only :  NCL,NCI, CDN3D,CRE3D,CLWP
#endif
  use CLOUDS_COM, only : TAUSS,TAUMC, CLDSS,CLDMC,CSIZMC,CSIZSS, &
       ULS,VLS,UMC,VMC,TLS,QLS,TAUSSIP,CSIZSSIP, &
       QLss,QIss,QLmc,QImc, &
       TMC,QMC,DDM1,AIRX,LMC,DDMS,TDN1,QDN1,DDML
#if (defined mjo_subdd) || (defined etc_subdd)
  use CLOUDS_COM, only : CLWC3D,CIWC3D,TLH3D,SLH3D,DLH3D,LLH3D
#endif
#ifdef etc_subdd
  use CLOUDS_COM, only : LWP2D,IWP2D
#endif
#ifdef mjo_subdd
  use CLOUDS_COM, only : TMCDRY,SMCDRY,DMCDRY,LSCDRY
#endif
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
  use CLOUDS_COM, only : NACTC,NAERC
#endif
#endif
#ifdef CACHED_SUBDD
#ifdef COSP_SIM
  ! use dictionary_mod, only : get_param, is_set_param, sync_param
  use dictionary_mod, only : get_param, is_set_param
  use clouds_com, only : nsubdd_cosp,flag_pfluxes_cosp,flag_re_cosp
  ! use domain_decomp_atm, only : am_i_root
#endif
#endif
  implicit none
  type (DIST_GRID), intent(IN) :: grid

  integer :: I_1H, I_0H, J_1H, J_0H
  integer :: IER

  I_0H = grid%I_STRT_HALO
  I_1H = grid%I_STOP_HALO
  J_0H = grid%J_STRT_HALO
  J_1H = grid%J_STOP_HALO

  allocate( &
       TTOLD(LM,I_0H:I_1H,J_0H:J_1H), &
       QTOLD(LM,I_0H:I_1H,J_0H:J_1H), &
       SVLHX(LM,I_0H:I_1H,J_0H:J_1H), &
       SVLAT(LM,I_0H:I_1H,J_0H:J_1H), &
       RHSAV(LM,I_0H:I_1H,J_0H:J_1H), &
       CLDSAV(LM,I_0H:I_1H,J_0H:J_1H), &
       CLDSAV1(LM,I_0H:I_1H,J_0H:J_1H), &
       FSS(LM,I_0H:I_1H,J_0H:J_1H), &
       STAT=IER)
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
  allocate( &
       CTEM(LM,I_0H:I_1H,J_0H:J_1H), &
       CD3D(LM,I_0H:I_1H,J_0H:J_1H), &
       CL3D(LM,I_0H:I_1H,J_0H:J_1H), &
       CI3D(LM,I_0H:I_1H,J_0H:J_1H), &
       STAT=IER)
#endif
#ifdef CLD_AER_CDNC
  allocate( &
       NCL(LM,I_0H:I_1H,J_0H:J_1H), &
       NCI(LM,I_0H:I_1H,J_0H:J_1H), &
       CDN3D(LM,I_0H:I_1H,J_0H:J_1H), &
       CRE3D(LM,I_0H:I_1H,J_0H:J_1H), &
       CLWP(I_0H:I_1H,J_0H:J_1H), &
       STAT=IER)
#endif
  allocate( &
       TAUSS(LM,I_0H:I_1H,J_0H:J_1H), &
       TAUSSIP(LM,I_0H:I_1H,J_0H:J_1H), &
       TAUMC(LM,I_0H:I_1H,J_0H:J_1H), &
       CLDSS(LM,I_0H:I_1H,J_0H:J_1H), &
       CLDMC(LM,I_0H:I_1H,J_0H:J_1H), &
       CSIZMC(LM,I_0H:I_1H,J_0H:J_1H), &
       CSIZSS(LM,I_0H:I_1H,J_0H:J_1H), &
       CSIZSSIP(LM,I_0H:I_1H,J_0H:J_1H), &
       QLss(LM,I_0H:I_1H,J_0H:J_1H), &
       QIss(LM,I_0H:I_1H,J_0H:J_1H), &
       QLmc(LM,I_0H:I_1H,J_0H:J_1H), &
       QImc(LM,I_0H:I_1H,J_0H:J_1H), &
       STAT=IER)
#ifdef mjo_subdd
  allocate( &
       TMCDRY(LM,I_0H:I_1H,J_0H:J_1H), &
       SMCDRY(LM,I_0H:I_1H,J_0H:J_1H), &
       DMCDRY(LM,I_0H:I_1H,J_0H:J_1H), &
       LSCDRY(LM,I_0H:I_1H,J_0H:J_1H), &
       STAT=IER)
#endif
#ifdef etc_subdd
  allocate( &
       LWP2D(I_0H:I_1H,J_0H:J_1H), &
       IWP2D(I_0H:I_1H,J_0H:J_1H), &
       STAT=IER)
#endif
#if (defined mjo_subdd) || (defined etc_subdd)
  allocate( &
       CLWC3D(LM,I_0H:I_1H,J_0H:J_1H), &
       CIWC3D(LM,I_0H:I_1H,J_0H:J_1H), &
       TLH3D(LM,I_0H:I_1H,J_0H:J_1H), &
       SLH3D(LM,I_0H:I_1H,J_0H:J_1H), &
       DLH3D(LM,I_0H:I_1H,J_0H:J_1H), &
       LLH3D(LM,I_0H:I_1H,J_0H:J_1H), &
       STAT=IER)
#endif

!**** Default values for certain arrays
  RHSAV (:,:,:)=.85d0
  CLDSAV(:,:,:)=0.
  SVLHX (:,:,:)=0.

  allocate(     ULS(I_0H:I_1H,J_0H:J_1H,LM), &
       VLS(I_0H:I_1H,J_0H:J_1H,LM), &
       UMC(I_0H:I_1H,J_0H:J_1H,LM), &
       VMC(I_0H:I_1H,J_0H:J_1H,LM), &
       TLS(I_0H:I_1H,J_0H:J_1H,LM), &
       QLS(I_0H:I_1H,J_0H:J_1H,LM), &
       TMC(I_0H:I_1H,J_0H:J_1H,LM), &
       QMC(I_0H:I_1H,J_0H:J_1H,LM), &
       STAT=IER)


!@var FSS initialized to 1.
  FSS = 1.
  TAUSS = 0.
  TAUSSIP = 0.
#ifdef CLD_AER_CDNC
!@var NCL is initialized to 10.0 cm-3
!@var NCI is initialised to 0.1 l^-1 or 10^-4 cm-3
  NCL = 10.
  NCI = 1.d-4
#endif

  allocate(     DDM1(I_0H:I_1H,J_0H:J_1H), &
       AIRX(I_0H:I_1H,J_0H:J_1H), &
       DDMS(I_0H:I_1H,J_0H:J_1H), &
       TDN1(I_0H:I_1H,J_0H:J_1H), &
       QDN1(I_0H:I_1H,J_0H:J_1H), &
       DDML(I_0H:I_1H,J_0H:J_1H), &
       STAT=IER)

  allocate(     LMC(2,I_0H:I_1H,J_0H:J_1H), &
       STAT=IER)

  !**** Initialise some output used in dynamics
  LMC(:,:,J_0H:J_1H)=0
  AIRX(:,J_0H:J_1H)=0.

#ifdef BLK_2MOM
#ifdef TRACERS_AMP
  allocate(  NACTC(LM,nmodes) )
  NACTC      = 1.0D-30
  allocate(  NAERC(LM,nmodes) )
  NAERC      = 1.0D-30
#endif
#endif
#ifdef CACHED_SUBDD
#ifdef COSP_SIM /* COSP_SIM */
  ! Trigger SUBDD COSP sampling at nsubdd_cosp
  !call sync_param('nsubdd_cosp', nsubdd_cosp, default=0)
  !call sync_param('nsubdd_cosp', nsubdd_cosp)
  !if( am_i_root() .and. nsubdd_cosp <= 0 )then
  !  call stop_model('nsubdd_cosp must be +ve when COSP_SIM defined',255)
  !endif
  if(is_set_param("nsubdd_cosp")) call get_param("nsubdd_cosp",nsubdd_cosp)
#ifdef COSP_PFLUX
  flag_pfluxes_cosp = .true.
#else
  flag_pfluxes_cosp = .false.
#endif
#ifdef COSP_USERE
  flag_re_cosp = .true.
#else
  flag_re_cosp = .false.
#endif
#endif /* COSP_SIM */
#endif

end subroutine ALLOC_CLOUDS_COM

subroutine def_rsf_clouds(fid)
!@sum  def_rsf_clouds defines cloud array structure in restart files
!@auth M. Kelley
!@ver  beta
  use clouds_com
  use domain_decomp_atm, only : grid
  use pario, only : defvar
  implicit none
  integer fid   !@var fid file id
  character(len=20) :: lijstr
  lijstr='(lm,dist_im,dist_jm)'
  call defvar(grid,fid,ttold,'ttold'//lijstr)
  call defvar(grid,fid,qtold,'qtold'//lijstr)
  call defvar(grid,fid,svlhx,'svlhx'//lijstr)
  call defvar(grid,fid,rhsav,'rhsav'//lijstr)
  call defvar(grid,fid,cldsav,'cldsav'//lijstr)
#ifdef CLD_AER_CDNC
  call defvar(grid,fid,ncl,'ncl'//lijstr)
  call defvar(grid,fid,nci,'nci'//lijstr)
#endif
  call defvar(grid,fid,airx,'airx(dist_im,dist_jm)')
  call defvar(grid,fid,lmc,'lmc(two,dist_im,dist_jm)')
  return
end subroutine def_rsf_clouds

subroutine new_io_clouds(fid,iaction)
!@sum  new_io_clouds read/write cloud arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
  use model_com, only : ioread,iowrite
  use domain_decomp_atm, only : grid
  use pario, only : write_dist_data,read_dist_data
  use clouds_com
  implicit none
  integer fid   !@var fid unit number of read/write
  integer iaction !@var iaction flag for reading or writing to file
  select case (iaction)
  case (iowrite)           ! output to restart file
    call write_dist_data(grid, fid, 'ttold', ttold, jdim=3)
    call write_dist_data(grid, fid, 'qtold', qtold, jdim=3)
    call write_dist_data(grid, fid, 'svlhx', svlhx, jdim=3)
    call write_dist_data(grid, fid, 'rhsav', rhsav, jdim=3)
    call write_dist_data(grid, fid, 'cldsav', cldsav, jdim=3)
#ifdef CLD_AER_CDNC
    call write_dist_data(grid, fid, 'ncl', ncl, jdim=3)
    call write_dist_data(grid, fid, 'nci', nci, jdim=3)
#endif
    call write_dist_data(grid, fid, 'airx', airx)
    call write_dist_data(grid, fid, 'lmc', lmc, jdim=3)
  case (ioread)            ! input from restart file
    call read_dist_data(grid, fid, 'ttold', ttold, jdim=3)
    call read_dist_data(grid, fid, 'qtold', qtold, jdim=3)
    call read_dist_data(grid, fid, 'svlhx', svlhx, jdim=3)
    call read_dist_data(grid, fid, 'rhsav', rhsav, jdim=3)
    call read_dist_data(grid, fid, 'cldsav', cldsav, jdim=3)
#ifdef CLD_AER_CDNC
    call read_dist_data(grid, fid, 'ncl', ncl, jdim=3)
    call read_dist_data(grid, fid, 'nci', nci, jdim=3)
#endif
    call read_dist_data(grid, fid, 'airx', airx)
    call read_dist_data(grid, fid, 'lmc', lmc, jdim=3)
  end select
  return
end subroutine new_io_clouds

#ifdef COSP_SIM /* COSP_SIM */
  subroutine save_cosp(grid)
  !------------------------------------------------------------------------------
  !@sum  Pass COSP output data structures to modelE SUBDD for storage.
  !@auth Mike Bauer
  !@usage Run once for time-steps where the COSP simulators are asked for after
  !@+       run_cosp_sims() is complete and before free_cosp_sims() is executed.
  !@calls remap_cosp_2d,remap_cosp_3d

    ! MODULE Imports
    ! =========================================================================

    ! Imported Type Definitions:
    use domain_decomp_atm, only : dist_grid
    !use cosp_drv, only : cfg,gbx
    ! Imported Parameters:
    use clouds_com, only : nsubdd_cosp
    use cosp_drv, only : n_out_list

    ! Imported Variables:
    use cosp_drv, only : cvar
    use geom, only : imaxj

    ! Imported Routines:
    use cosp_drv, only : remap_cosp_2d,remap_cosp_3d
    use subdd_mod, only : inc_subdd

    ! All variables, parameters, and functions must be declared
    ! =========================================================================
    implicit none

    ! Subroutine Arguments, in order of appearance
    ! =========================================================================
    type (dist_grid),intent(in) :: grid

    ! Subroutine Scalers/Variables
    ! =========================================================================
    integer :: i_cosp!np_cosp,i_cosp,j_cosp
    real*8,dimension(:,:),allocatable :: blob2d
    real*8,dimension(:,:,:),allocatable :: blob3d
    ! real*8,dimension(:,:,:,:),allocatable :: blob4d
    integer,dimension(8) :: geom_2d
    ! ^^^^^^^^^^^^^^^^^^^^^ end of declaration of variables ^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! Deal with domain decomposition
    geom_2d = (/                                                                &
    ! 1=I_0       2=I_1       3=J_0       4=J_1
      grid%i_strt,grid%i_stop,grid%j_strt,grid%j_stop,                          &
    ! 5=I_0H           6=5=I_1H         7=J_0H           8=J_1H
      grid%i_strt_halo,grid%i_stop_halo,grid%j_strt_halo,grid%j_stop_halo /)

    ! General COSP Variables (e.g, axis data like IJ layers)
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! call inc_subdd(
    !   'some_var',               ! name of this field
    !   real*8  some_array,       ! data to be saved/accumulated
    !   integer Nsubdd,           ! output frequency for this field
    !   logical instant,          ! T for snapshots, F for averages
    !   optional units='xxxx',    ! units specification
    !   optional long_name='xxxx' ! description
    !   dim3name='isccp_ntau',
    !   dim4name='isccp_npres'
    !   )
    ! Generally ends up calling SUBDD.f:inc_subdd_solo_2d()
    !   inc_subdd_solo_2d(vname,arr,nsubdd,inst,units,long_name)

    do i_cosp=1,n_out_list
#ifdef COSP_DEBUG
        write(0, *) "i_cosp ",i_cosp
#endif
      if (cvar(i_cosp)%cflag) then
#ifdef COSP_DEBUG
        write(0, *) "cvar(i_cosp)%cflag ",cvar(i_cosp)%cflag
#endif
        ! Requested COSP Output
        ! ---------------------
        if (cvar(i_cosp)%is_3d) then
#ifdef COSP_DEBUG
         write(0, *) "3d ",trim(cvar(i_cosp)%sname)," ", &
          cvar(i_cosp)%dim1," ",cvar(i_cosp)%dim2
#endif
          ! 3D Output
          ! ----------
          call remap_cosp_3d(trim(cvar(i_cosp)%sname),imaxj,geom_2d,           &
            cvar(i_cosp)%dim1,cvar(i_cosp)%dim2,blob3d)

          call inc_subdd(trim(cvar(i_cosp)%sname(2:16)),blob3d,nsubdd_cosp,    &
            .true.,units=trim(cvar(i_cosp)%units),                             &
            long_name=trim(cvar(i_cosp)%lname))
        else
          ! 2D Output
          ! ----------
#ifdef COSP_DEBUG
          write(0, *) "2d ",trim(cvar(i_cosp)%sname)
#endif
          call remap_cosp_2d(trim(cvar(i_cosp)%sname),imaxj,geom_2d,blob2d)
          call inc_subdd(trim(cvar(i_cosp)%sname(2:16)),blob2d,nsubdd_cosp,    &
            .true.,units=trim(cvar(i_cosp)%units),                             &
            long_name=trim(cvar(i_cosp)%lname))
        endif ! is_3d
      endif ! cflag
    enddo
  end subroutine save_cosp
#endif /* COSP_SIM */
