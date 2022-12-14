#include "rundeck_opts.h"
      MODULE SOMTQ_COM
!@sum  SOMTQ_COM contains the arrays containing second order moments
!@auth Gary Russell
      USE QUSDEF
      USE RESOLUTION, only : im,jm,lm
      IMPLICIT NONE
      SAVE
!     REAL*8, DIMENSION(NMOM,IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO ,LM)  
!    &        :: TMOM,QMOM
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TMOM,QMOM

      END MODULE SOMTQ_COM

      SUBROUTINE ALLOC_SMOMTQ(grid)
!@sum  init_smomtq allocates the arrays in this module which
!@+    must now be dynamic for the distributed memory implementation.
!@auth Rosalinda de Fainchtein
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE QUSDEF, ONLY : NMOM
      USE RESOLUTION, ONLY : LM
      USE SOMTQ_COM, ONLY : TMOM,QMOM
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE ( TMOM(NMOM , I_0H:I_1H , J_0H:J_1H , LM),
     &           QMOM(NMOM , I_0H:I_1H , J_0H:J_1H , LM),
     &   STAT=IER )

      TMOM = 0.
      QMOM = 0.

      END SUBROUTINE ALLOC_SMOMTQ

      subroutine def_rsf_somtq(fid)
!@sum  def_rsf_somtq defines QUS T/Q array structure in restart files
!@auth M. Kelley
!@ver  beta
      use somtq_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,tmom,'tmom(nmom,dist_im,dist_jm,lm)')
      call defvar(grid,fid,qmom,'qmom(nmom,dist_im,dist_jm,lm)')
      return
      end subroutine def_rsf_somtq

      subroutine new_io_somtq(fid,iaction)
!@sum  new_io_somtq read/write QUS T/Q arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use somtq_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to restart file
        call write_dist_data(grid, fid, 'tmom', tmom, jdim=3)
        call write_dist_data(grid, fid, 'qmom', qmom, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'tmom', tmom, jdim=3)
        call read_dist_data(grid, fid, 'qmom', qmom, jdim=3)
      end select
      return
      end subroutine new_io_somtq

      subroutine tq_zmom_init(t,q,pmid,pedn)
      USE DOMAIN_DECOMP_ATM, ONLY: grid
      USE SOMTQ_COM
      implicit none
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) :: t,q
      REAL*8 :: pmid(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo)
      REAL*8 :: pedn(lm+1,grid%i_strt_halo:grid%i_stop_halo,
     &                    grid%j_strt_halo:grid%j_stop_halo)
      integer :: i,j,l
      REAL*8 :: rdsig

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S
C****
C**** Extract useful local domain parameters from "grid"
C****
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP

C**** INITIALIZES VERTICAL SLOPES OF T,Q
      DO J=J_0,J_1
        DO I=I_0,I_1
          RDSIG=(PMID(1,I,J)-PEDN(2,I,J))/(PMID(1,I,J)-PMID(2,I,J))
          TMOM(MZ,I,J,1)=(T(I,J,2)-T(I,J,1))*RDSIG
          QMOM(MZ,I,J,1)=(Q(I,J,2)-Q(I,J,1))*RDSIG
          IF(Q(I,J,1)+QMOM(MZ,I,J,1).LT.0.) QMOM(MZ,I,J,1)=-Q(I,J,1)
          DO L=2,LM-1
            RDSIG=(PMID(L,I,J)-PEDN(L+1,I,J))/
     *           (PMID(L-1,I,J)-PMID(L+1,I,J))
            TMOM(MZ,I,J,L)=(T(I,J,L+1)-T(I,J,L-1))*RDSIG
            QMOM(MZ,I,J,L)=(Q(I,J,L+1)-Q(I,J,L-1))*RDSIG
            IF(Q(I,J,L)+QMOM(MZ,I,J,L).LT.0.) QMOM(MZ,I,J,L)=-Q(I,J,L)
          END DO
          RDSIG=(PMID(LM,I,J)-PEDN(LM+1,I,J))/
     *         (PMID(LM-1,I,J)-PMID(LM,I,J))
          TMOM(MZ,I,J,LM)=(T(I,J,LM)-T(I,J,LM-1))*RDSIG
          QMOM(MZ,I,J,LM)=(Q(I,J,LM)-Q(I,J,LM-1))*RDSIG
          IF(Q(I,J,LM)+QMOM(MZ,I,J,LM).LT.0.) QMOM(MZ,I,J,LM)=-Q(I,J,LM)
        END DO
      END DO
      return
      end subroutine tq_zmom_init
