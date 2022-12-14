#include "rundeck_opts.h"
c-----------------------------------------------------------------------------
      module hycom_dim
      USE DOMAIN_DECOMP_1D, only : DIST_GRID

      implicit none

      private

      public alloc_hycom_dim,init_hycom_grid

      public ii,jj,kk,ii1
      public I_0,  I_1,  J_0,  J_1, I_0H, I_1H, J_0H, J_1H
      public aI_0, aI_1, aJ_0, aJ_1,aI_0H,aI_1H,aJ_0H,aJ_1H
      public ogrid

      public ip,iu,iv,iq,
     .ifp,ilp,isp,jfp,jlp,jsp,
     .ifq,ilq,isq,jfq,jlq,jsq,
     .ifu,ilu,isu,jfu,jlu,jsu,
     .ifv,ilv,isv,jfv,jlv,jsv,
     .msk

#ifdef HYCOM2deg
      integer, public, parameter :: idm=195,jdm=180,kdm=26,ms=15
     .                             ,iold=181
#endif
#ifdef HYCOM1degRefined
      integer, public, parameter :: idm=387,jdm=360,ms=20,iold=359
#endif
#ifdef HYCOM1degUnrefined
      integer, public, parameter :: idm=359,jdm=360,ms=20,iold=359
#endif
#ifdef HYCOM26layers
      integer, public, parameter:: kdm=26
#endif
#ifdef HYCOM32layers
      integer, public, parameter:: kdm=32
#endif
#ifdef ATM4x5
      integer, public, parameter :: iia=72,jja=46
#endif
#ifdef ATM2x2h
      integer, public, parameter :: iia=144,jja=90
#endif
      integer, public, parameter :: iio=idm,jjo=jdm

#if defined(TRACERS_HYCOM_Ventilation) || defined(TRACERS_AGE_OCEAN) \
    || defined(TRACERS_OCEAN_WATER_MASSES) || defined(TRACERS_ZEBRA)
#ifdef TRACERS_ZEBRA
      integer, public, parameter :: ntrcr_zebra = 26
#else
      integer, public, parameter :: ntrcr_zebra = 0

#endif
#ifdef TRACERS_HYCOM_Ventilation
      integer, public, parameter :: ntrcr_vent = 1
#else
      integer, public, parameter :: ntrcr_vent = 0
#endif
#ifdef TRACERS_AGE_OCEAN
      integer, public, parameter :: ntrcr_age = 1
#else
      integer, public, parameter :: ntrcr_age = 0
#endif
#ifdef TRACERS_OCEAN_WATER_MASSES 
      integer, public, parameter :: ntrcr_wm = 2
#else
      integer, public, parameter :: ntrcr_wm = 0
#endif
#else
      integer, public, parameter :: ntrcr_zebra = 0
      integer, public, parameter :: ntrcr_vent = 0
      integer, public, parameter :: ntrcr_age = 0
      integer, public, parameter :: ntrcr_wm = 0
#endif

#ifdef TRACERS_OceanBiology
#ifdef TRACERS_Alkalinity
      integer, public, parameter :: ntrcr_obio = 16
#else
      integer, public, parameter :: ntrcr_obio = 15
#endif
#else
      integer, public, parameter :: ntrcr_obio = 0
#endif  /*oceanbiology*/

      integer, public ::  ntrcr = max(1,ntrcr_zebra + ntrcr_vent 
     .                          + ntrcr_age + ntrcr_wm + ntrcr_obio)

c
c --- ms-1  = max. number of interruptions of any grid row or column by land
c
      integer ii,jj,kk,ii1,nlato,nlongo,nlatn,nlongn
      parameter (ii=idm,jj=jdm,kk=kdm,ii1=ii-1)
      parameter (nlato=1,nlongo=1,nlatn=idm,nlongn=jdm)
c
c --- information in common block gindex keeps do loops from running into land
      integer, allocatable :: ip(:,:),iu(:,:),iv(:,:),iq(:,:),
     .ifp(:,:),ilp(:,:),isp(:),jfp(:,:),jlp(:,:),jsp(:),
     .ifq(:,:),ilq(:,:),isq(:),jfq(:,:),jlq(:,:),jsq(:),
     .ifu(:,:),ilu(:,:),isu(:),jfu(:,:),jlu(:,:),jsu(:),
     .ifv(:,:),ilv(:,:),isv(:),jfv(:,:),jlv(:,:),jsv(:),
     .msk(:,:)
c
!!      integer ip,iu,iv,iq,
!!     .        ifp,ilp,isp,jfp,jlp,jsp,ifq,ilq,isq,jfq,jlq,jsq,
!!     .        ifu,ilu,isu,jfu,jlu,jsu,ifv,ilv,isv,jfv,jlv,jsv,msk

      ! MPI decomposition stuff (maybe should move to a separate module?)
      TYPE(DIST_GRID), target :: ogrid   ! ocean (hycom) grid
      ! domain bounds
      integer :: I_0,  I_1,  J_0,  J_1
      integer ::aI_0, aI_1, aJ_0, aJ_1
      ! domain bounds with halos
      integer :: I_0H, I_1H, J_0H, J_1H
      integer ::aI_0H,aI_1H,aJ_0H,aJ_1H

      ! openmp decomposition parameter
      integer, public :: jchunk

      contains

      subroutine init_hycom_grid
      USE DOMAIN_DECOMP_1D, only : init_grid, getDomainBounds

      call init_grid( ogrid, idm, jdm, 1, bc_periodic=.true. )
      call getDomainBounds(ogrid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H ,
     &               I_STRT     =I_0,    I_STOP     =I_1,
     &               I_STRT_HALO=I_0H,   I_STOP_HALO=I_1H )
 
      ! hack to work with older code on 1 cpu
      print *,"init_hycom_grid: i0,i1,j0,j1", I_0,I_1,J_0,J_1
      print *,"init_hycom_grid: halos", I_0H,I_1H,J_0H,J_1H
!      J_0H = J_0
!      J_1H = J_1
!      I_0H = I_0
!      I_1H = I_1

!      I_0H = 1
!      I_1H = idm
!      J_0H = 1
!      J_1H = jdm


      end subroutine init_hycom_grid

      subroutine alloc_hycom_dim

      ! leave these global for time being

      allocate(
     . ip(I_0H:I_1H,J_0H:J_1H),iu(I_0H:I_1H,J_0H:J_1H),
     . iv(I_0H:I_1H,J_0H:J_1H),iq(I_0H:I_1H,J_0H:J_1H),
     . ifp(J_0H:J_1H,ms),ilp(J_0H:J_1H,ms),isp(J_0H:J_1H),
     . jfp(I_0H:I_1H,ms),jlp(I_0H:I_1H,ms),jsp(I_0H:I_1H),
     . ifq(J_0H:J_1H,ms),ilq(J_0H:J_1H,ms),isq(J_0H:J_1H),
     . jfq(I_0H:I_1H,ms),jlq(I_0H:I_1H,ms),jsq(I_0H:I_1H),
     . ifu(J_0H:J_1H,ms),ilu(J_0H:J_1H,ms),isu(J_0H:J_1H),
     . jfu(I_0H:I_1H,ms),jlu(I_0H:I_1H,ms),jsu(I_0H:I_1H),
     . ifv(J_0H:J_1H,ms),ilv(J_0H:J_1H,ms),isv(J_0H:J_1H),
     . jfv(I_0H:I_1H,ms),jlv(I_0H:I_1H,ms),jsv(I_0H:I_1H),
     . msk(I_0H:I_1H,J_0H:J_1H) )

cddd      allocate(
cddd     . ip(idm,jdm),iu(idm,jdm),iv(idm,jdm),iq(idm,jdm),
cddd     .ifp(jdm,ms),ilp(jdm,ms),isp(jdm),jfp(idm,ms),jlp(idm,ms),jsp(idm),
cddd     .ifq(jdm,ms),ilq(jdm,ms),isq(jdm),jfq(idm,ms),jlq(idm,ms),jsq(idm),
cddd     .ifu(jdm,ms),ilu(jdm,ms),isu(jdm),jfu(idm,ms),jlu(idm,ms),jsu(idm),
cddd     .ifv(jdm,ms),ilv(jdm,ms),isv(jdm),jfv(idm,ms),jlv(idm,ms),jsv(idm),
cddd     .msk(idm,jdm) )

      end subroutine alloc_hycom_dim

      end module hycom_dim
c-----------------------------------------------------------------------------
