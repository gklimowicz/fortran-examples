#include "rundeck_opts.h"

      module veg_com
!@sum  GHY_COM contains the areas used by the Ground Hydrology routines
!@auth N. Kiang, I. Aleinov
      use ghy_com, only : ngm
      implicit none
      save

!---  boundary conditions (read from file)
!@var vdata(:,:,k)  fraction of gridbox of veg.type k=1-12
      real*8, ALLOCATABLE, dimension(:,:,:) :: vdata

!---  prognostic variables (saved to restart file)
!@var Cint Internal foliage CO2 concentration (mol/m3)
      real*8, ALLOCATABLE, dimension(:,:) :: Cint
!@var Qfol Foliage surface mixing ratio (kg/kg)
      real*8, ALLOCATABLE, dimension(:,:) :: Qfol
!@var cnc_ij canopy conductance
      real*8, ALLOCATABLE, dimension(:,:) :: cnc_ij

!---  work arrays (recomputed for each restart)
      real*8, ALLOCATABLE, dimension(:,:,:) :: afr
      real*8, ALLOCATABLE, dimension(:,:,:) :: ala,acs,almass !nyk almass=leaf mass
      real*8, ALLOCATABLE, dimension(:,:,:,:) :: alaf     !nyk lai components by vegtype
      real*8, ALLOCATABLE, dimension(:,:,:) :: alaif   !nyk lai by vegtype (not * vdata!)
      real*8, ALLOCATABLE, dimension(:,:) :: afb,avh,aalbveg !nyk aalbveg
      real*8, ALLOCATABLE, dimension(:,:) :: can_w_capacity
      real*8, ALLOCATABLE, dimension(:,:) :: anm,anf ! adf

      end module veg_com

      subroutine def_rsf_vegetation(fid)
!@sum  def_rsf_vegetation defines vegetation array structure in restart files
!@auth M. Kelley
!@ver  beta
      use veg_com, only : Cint, Qfol, cnc_ij
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,cint,'cint(dist_im,dist_jm)')
      call defvar(grid,fid,qfol,'qfol(dist_im,dist_jm)')
      call defvar(grid,fid,cnc_ij,'cnc_ij(dist_im,dist_jm)')
      return
      end subroutine def_rsf_vegetation

      subroutine new_io_vegetation(fid,iaction)
!@sum  new_io_vegetation read/write vegetation arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use veg_com, only : Cint, Qfol, cnc_ij
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid, fid, 'cint', cint)
        call write_dist_data(grid, fid, 'qfol', qfol)
        call write_dist_data(grid, fid, 'cnc_ij', cnc_ij)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'cint', cint)
        call read_dist_data(grid, fid, 'qfol', qfol)
        call read_dist_data(grid, fid, 'cnc_ij', cnc_ij)
      end select
      return
      end subroutine new_io_vegetation

      SUBROUTINE ALLOC_VEG_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE VEG_COM
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

      ALLOCATE(    vdata(I_0H:I_1H,J_0H:J_1H,12),
     *              Cint(I_0H:I_1H,J_0H:J_1H),
     *              Qfol(I_0H:I_1H,J_0H:J_1H),
     *            cnc_ij(I_0H:I_1H,J_0H:J_1H),
     *               afr(ngm,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(    ala(3,  I_0H:I_1H,J_0H:J_1H),
     *             acs(3,  I_0H:I_1H,J_0H:J_1H),
     *          almass(3,  I_0H:I_1H,J_0H:J_1H),
     *            alaf(3,11,I_0H:I_1H,J_0H:J_1H),
     *           alaif(11,  I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(   afb(I_0H:I_1H,J_0H:J_1H),
     *            avh(I_0H:I_1H,J_0H:J_1H),
     *            aalbveg(I_0H:I_1H,J_0H:J_1H),
     *            can_w_capacity(I_0H:I_1H,J_0H:J_1H),
     *            anm(I_0H:I_1H,J_0H:J_1H),
     *            anf(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      END SUBROUTINE ALLOC_VEG_COM

