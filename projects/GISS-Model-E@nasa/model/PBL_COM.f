#include "rundeck_opts.h"

      MODULE PBLCOM
!@sum  PBLCOM contains the arrays used by the Boundary Layer code
!@auth Greg Hartke/Ye Cheng
      USE SOCPBL, only : npbl=>n
      IMPLICIT NONE
      SAVE

!@var ROUGHL log10(zgs/roughness length), prescribed with zgs=30 m.
      REAL*8, allocatable, dimension(:,:) :: roughl

!@var DCLEV     LAYER TO WHICH DRY CONVECTION MIXES (1)
!@var ugeo,vgeo components of geostrophic wind at the top of the BL
!@var bldep     boundary layer depth (m)
      REAL*8, allocatable, dimension(:,:) ::
     &     dclev,pblht,pblptop,ugeo,vgeo,bldep

!@var [tuv]1_after_aturb first-layer temp/winds after ATURB completes
!@+   (used to compute tendencies seen by the PBL code)
      REAL*8, allocatable, dimension(:,:) ::
     &     t1_after_aturb,u1_after_aturb,v1_after_aturb

!@var egcm  3-d turbulent kinetic energy in the whole atmosphere
!@var w2gcm vertical component of egcm
!@var t2gcm 3-d turbulent temperature variance in the whole atmosphere
      real*8, allocatable, dimension(:,:,:) :: egcm,w2gcm,t2gcm

!@ egcm_init_max maximum initial vaule of egcm
      real*8, parameter :: egcm_init_max=0.5d0

      END MODULE PBLCOM

      subroutine def_rsf_pbl(fid)
!@sum  def_rsf_pbl defines pbl array structures in restart files
!@auth M. Kelley
!@ver  beta
      use pblcom
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      character(len=29) :: dimstr
      dimstr='(dist_im,dist_jm)'
      call defvar(grid,fid,t1_after_aturb,'t1_after_aturb'//dimstr)
      call defvar(grid,fid,u1_after_aturb,'u1_after_aturb'//dimstr)
      call defvar(grid,fid,v1_after_aturb,'v1_after_aturb'//dimstr)
      return
      end subroutine def_rsf_pbl

      subroutine new_io_pbl(fid,iaction)
!@sum  new_io_pbl read/write pbl arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use pblcom
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 't1_after_aturb',t1_after_aturb)
        call write_dist_data(grid, fid, 'u1_after_aturb',u1_after_aturb)
        call write_dist_data(grid, fid, 'v1_after_aturb',v1_after_aturb)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 't1_after_aturb',t1_after_aturb)
        call read_dist_data(grid, fid, 'u1_after_aturb',u1_after_aturb)
        call read_dist_data(grid, fid, 'v1_after_aturb',v1_after_aturb)
      end select
      return
      end subroutine new_io_pbl

      subroutine def_rsf_bldat(fid)
!@sum  def_rsf_bldat defines bldat array structure in restart files
!@auth M. Kelley
!@ver  beta
      use pblcom
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,dclev,'dclev(dist_im,dist_jm)')
      call defvar(grid,fid,pblht,'pblht(dist_im,dist_jm)')
      call defvar(grid,fid,pblptop,'pblptop(dist_im,dist_jm)')
      call defvar(grid,fid,egcm,'egcm(lm,dist_im,dist_jm)')
      call defvar(grid,fid,w2gcm,'w2gcm(lm,dist_im,dist_jm)')
      return
      end subroutine def_rsf_bldat

      subroutine new_io_bldat(fid,iaction)
!@sum  new_io_bldat read/write bldat arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use pblcom
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'dclev',dclev)
        call write_dist_data(grid,fid,'pblht',pblht)
        call write_dist_data(grid,fid,'pblptop',pblptop)
        call write_dist_data(grid,fid,'egcm',egcm, jdim=3)
        call write_dist_data(grid,fid,'w2gcm',w2gcm, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'dclev',dclev)
        call read_dist_data(grid,fid,'pblht',pblht)
        call read_dist_data(grid,fid,'pblptop',pblptop)
        call read_dist_data(grid,fid,'egcm',egcm, jdim=3)
        call read_dist_data(grid,fid,'w2gcm',w2gcm, jdim=3)
      end select
      return
      end subroutine new_io_bldat

      SUBROUTINE ALLOC_PBL_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE RESOLUTION, only : lm
      USE CONSTANT, only : by3
      USE PBLCOM
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, getDomainBounds
      USE FLUXES, only : atmocns,atmices,atmglas,atmlnds,asflx
#ifdef GLINT2
      USE FLUXES, only : atmglas_hp
#endif

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER, IP3, K, L

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE(    roughl(I_0H:I_1H,J_0H:J_1H),
     *              dclev(I_0H:I_1H,J_0H:J_1H),
     *              pblht(I_0H:I_1H,J_0H:J_1H),
     *              pblptop(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

C**** SET LAYER THROUGH WHICH DRY CONVECTION MIXES TO 1
      DCLEV(:,:)=1.

      ALLOCATE(    egcm(lm,I_0H:I_1H,J_0H:J_1H),
     *            w2gcm(lm,I_0H:I_1H,J_0H:J_1H),
     *            t2gcm(lm,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(
     &         ugeo(I_0H:I_1H,J_0H:J_1H),
     &         vgeo(I_0H:I_1H,J_0H:J_1H),
     &         bldep(I_0H:I_1H,J_0H:J_1H)
     &     )

      ip3 = 0
      do k=1,size(asflx)
        select case(asflx(k)%itype4)
        case (1)
          call alloc_atmsrf_pbl_vars(grid,
     &         atmocns(1)%atmsrf_xchng_vars,asflx(k))
        case (2)
          call alloc_atmsrf_pbl_vars(grid,
     &         atmices(1)%atmsrf_xchng_vars,asflx(k))
        case (3)
          ip3 = ip3 + 1
#ifdef GLINT2
          call alloc_atmsrf_pbl_vars(grid,
     &         atmglas_hp(ip3)%atmsrf_xchng_vars,asflx(k))
#endif
      ! dest in alloc_atmsrf_pbl_vars() should point to atmglas
      ! in the end, not atmglas_hp
          call alloc_atmsrf_pbl_vars(grid,
     &         atmglas(ip3)%atmsrf_xchng_vars,asflx(k))
        case (4)
          call alloc_atmsrf_pbl_vars(grid,
     &         atmlnds(1)%atmsrf_xchng_vars,asflx(k))
        end select
      enddo

      ALLOCATE(t1_after_aturb(I_0H:I_1H,J_0H:J_1H),
     &         u1_after_aturb(I_0H:I_1H,J_0H:J_1H),
     &         v1_after_aturb(I_0H:I_1H,J_0H:J_1H))
      t1_after_aturb(:,J_0H:J_1H) = 0.
      u1_after_aturb(:,J_0H:J_1H) = 0.
      v1_after_aturb(:,J_0H:J_1H) = 0.

C**** initialize egcm to be used in ATURB.f
      DO L=1,LM
        egcm(l,:,:)=egcm_init_max/(float(l)**2)
        w2gcm(l,:,:)=egcm(l,:,:)*2.*by3
      END DO

      END SUBROUTINE ALLOC_PBL_COM

      subroutine alloc_atmsrf_pbl_vars(grid,this,dest)
      use pblcom
      use exchange_types, only : atmsrf_xchng_vars
      use domain_decomp_1d, only : dist_grid, getDomainBounds
      implicit none
      type (dist_grid), intent(in) :: grid
      type(atmsrf_xchng_vars) :: this,dest
      integer :: i_0h, i_1h, j_1h, j_0h

      call getDomainBounds(grid, j_strt_halo=j_0h, j_stop_halo=j_1h)
      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo

      allocate(
     &     this % uabl(npbl,i_0h:i_1h,j_0h:j_1h)
     &    ,this % vabl(npbl,i_0h:i_1h,j_0h:j_1h)
     &    ,this % tabl(npbl,i_0h:i_1h,j_0h:j_1h)
     &    ,this % qabl(npbl,i_0h:i_1h,j_0h:j_1h)
     &    ,this % eabl(npbl,i_0h:i_1h,j_0h:j_1h)
#ifdef TRACERS_ON
     &    ,this % trabl(npbl,this%ntm,i_0h:i_1h,j_0h:j_1h)
#endif
     &     )
      this % qabl = 0.  ! initialise to make life easier

      dest%uabl => this%uabl
      dest%vabl => this%vabl
      dest%tabl => this%tabl
      dest%qabl => this%qabl
      dest%eabl => this%eabl
#ifdef TRACERS_ON
      dest%trabl => this%trabl
#endif
      return
      end subroutine alloc_atmsrf_pbl_vars
