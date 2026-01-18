#include "rundeck_opts.h"

      MODULE INT_AG2OG_MOD

!@sum INT_AG2OG_MOD contains subroutines for conversion 2D, 3D, etc. 
!!    arrays from atm. to the ocean grid 
!@auth Larissa Nazarenko
      use hntrp_mod
      IMPLICIT NONE
      SAVE
      PRIVATE
      PUBLIC INT_AG2OG

      Interface INT_AG2OG
c _par are parallelized versions
      Module Procedure INT_AG2OG_2Da_par
      Module Procedure INT_AG2OG_3Da_par
      Module Procedure INT_AG2OG_Vector1_par
!      Module Procedure INT_AG2OG_3Dc_par
!      Module Procedure INT_AG2OG_Vector2
!      Module Procedure INT_AG2OG_2Db ! need _par for tracers
!      Module Procedure INT_AG2OG_3Db_par
!      Module Procedure INT_AG2OG_4Da ! trgasex
!      Module Procedure INT_AG2OG_4Db ! need _par for tracers
      End Interface

      contains

      SUBROUTINE INT_AG2OG_2Da_par(aA,oA,aWEIGHT,copypole)

!@sum INT_AG2OG_2Da_par parallel version of INT_AG2OG_2Da
!@auth M. Kelley

      USE OCEAN,      only : oIM=>im,oJM=>jm
      Use OCEANR_DIM,       only : oGRID
      USE OCEAN, only : hntrp_a2o => remap_a2o
      IMPLICIT NONE

      real*8, dimension(hntrp_a2o%ima,
     &     hntrp_a2o%bpack%J_STRT_HALO:hntrp_a2o%bpack%J_STOP_HALO),
     &     intent(in)  :: aWEIGHT, aA
      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(out)  :: oA
      integer :: jmin,jmax
      real*8, dimension(:,:), allocatable :: aA_band,aWEIGHT_band
      logical, intent(in), optional :: copypole

      logical :: copypole_
      integer :: aIM,aJM

      aIM = hntrp_a2o%ima
      aJM = hntrp_a2o%jma_full

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(:aIM,:) = aA(:,:)
        return
      endif

C***  Gather the requisite atm latitude bands
      jmin = hntrp_a2o%bpack%jband_strt
      jmax = hntrp_a2o%bpack%jband_stop
      ALLOCATE(aA_band(aIM,jmin:jmax),aWEIGHT_band(aIM,jmin:jmax))
      CALL BAND_PACK (hntrp_a2o%bpack, aA, aA_band)
      CALL BAND_PACK (hntrp_a2o%bpack, aWEIGHT, aWEIGHT_band)

C***  Interpolate aA from atmospheric grid to ocean grid 
      copypole_ = .true.
      if(present(copypole)) copypole_ = copypole
      !if(copypole_ .and. hasNorthPole(aGRID)) then
      if(copypole_ .and. jmax==aJM) then
        aA_band(2:aIM,aJM) = aA_band(1,aJM)
      endif
      if(copypole_) then
        call HNTR8P_band(aWEIGHT_band, aA_band, hntrp_a2o, oA)
      else
        call HNTR8_band (aWEIGHT_band, aA_band, hntrp_a2o, oA)
      endif
      DEALLOCATE(aA_band,aWEIGHT_band)

      RETURN
      END SUBROUTINE INT_AG2OG_2Da_par

      SUBROUTINE INT_AG2OG_3Da_par(aA,oA,aWEIGHT,NT)
!@sum INT_OG2AG_3Da_par parallel version of INT_OG2AG_3Da
!@auth M. Kelley
      USE OCEAN,      only : oIM=>im,oJM=>jm
      Use OCEANR_DIM,       only : oGRID
      USE OCEAN, only : hntrp_a2o => remap_a2o
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NT

      real*8, dimension(NT,hntrp_a2o%ima,
     &     hntrp_a2o%bpack%J_STRT_HALO:hntrp_a2o%bpack%J_STOP_HALO),
     &     intent(in) :: aA
      real*8, dimension(hntrp_a2o%ima,
     &     hntrp_a2o%bpack%J_STRT_HALO:hntrp_a2o%bpack%J_STOP_HALO),
     &     intent(in) :: aWEIGHT
      real*8, dimension(NT,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(out)  :: oA

      real*8, dimension(:,:,:), allocatable :: aA_band
      real*8, dimension(:,:), allocatable :: aWEIGHT_band,oA2d,aA2d
      integer :: i,j,n, jmin,jmax
      integer :: aIM,aJM

      aIM = hntrp_a2o%ima
      aJM = hntrp_a2o%jma_full

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(:,:aIM,:) = aA(:,:,:)
        return
      endif

C***  Gather the requisite atm latitude bands
      jmin = hntrp_a2o%bpack%jband_strt
      jmax = hntrp_a2o%bpack%jband_stop
      ALLOCATE(aA_band(NT,aIM,jmin:jmax),aWEIGHT_band(aIM,jmin:jmax))
      CALL BAND_PACK_COLUMN (hntrp_a2o%bpack, aA, aA_band)
      CALL BAND_PACK (hntrp_a2o%bpack, aWEIGHT, aWEIGHT_band)

      allocate(aA2d(aIM,jmin:jmax))
      allocate(oA2d(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

C***  Interpolate aA from atmospheric grid to ocean grid 
      do n=1,NT
        aA2d(:,:) = aA_band(n,:,:)
        !if(hasNorthPole(aGRID)) then
        if(jmax==aJM) then
          aA2D(2:aIM,aJM) = aA2D(1,aJM)
        endif
        call HNTR8P_band(aWEIGHT_band, aA2D, hntrp_a2o, oA2d)
        oA(n,:,:) = oA2d(:,:)
      enddo

      DEALLOCATE(aA_band, aWEIGHT_band, oA2d,aA2d)

      RETURN
      END SUBROUTINE INT_AG2OG_3Da_par

      SUBROUTINE INT_AG2OG_Vector1_par(aU,aV,oU,oV,aWEIGHT,aFOCEAN,
     &     aSINI,aCOSI)
!@sum INT_AG2OG_Vector1_par parallel version of INT_AG2OG_Vector1
!@auth M. Kelley
      USE OCEAN,      only : oIM=>im,oJM=>jm
      USE OCEAN, only : oSINI=>SINIC,oCOSI=>COSIC
      USE DOMAIN_DECOMP_1D, only : hasNorthPole, hasSouthPole
      USE OCEANR_DIM,       only : oGRID
      USE OCEAN, only : hntrp_a2o => remap_a2o
      IMPLICIT NONE

      real*8, dimension(hntrp_a2o%ima,
     &     hntrp_a2o%bpack%J_STRT_HALO:hntrp_a2o%bpack%J_STOP_HALO),
     &     intent(in) :: aFOCEAN, aWEIGHT
      real*8, dimension(hntrp_a2o%ima) :: aSINI,aCOSI
      real*8, dimension(hntrp_a2o%ima,
     &     hntrp_a2o%bpack%J_STRT_HALO:hntrp_a2o%bpack%J_STOP_HALO),
     &     intent(in) :: aU, aV
      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(out) :: oU, oV

C***  Local vars
      integer :: i,j
      real*8 :: aUnp,aVnp, oUnp,oVnp
      real*8, dimension(:,:), allocatable :: aUtmp,aVtmp
      integer :: aIM,aJM,aJ_0,aJ_1,aJ_0H,aJ_1H

      aIM = hntrp_a2o%ima
      aJM = hntrp_a2o%jma_full
      aJ_0 = hntrp_a2o%bpack%J_STRT
      aJ_1 = hntrp_a2o%bpack%J_STOP
      aJ_0H = hntrp_a2o%bpack%J_STRT_HALO
      aJ_1H = hntrp_a2o%bpack%J_STOP_HALO

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oU(:aIM,:) = aU(:,:)
        oV(:aIM,:) = aV(:,:)
        return
      endif

      allocate(aUtmp(aIM,aJ_0H:aJ_1H),aVtmp(aIM,aJ_0H:aJ_1H))

      do j=max(2,aJ_0),min(aJ_1,aJM-1)
        do i=1,aIM
          if(aFOCEAN(i,j).gt.0.) then
            aUtmp(i,j) = aU(i,j)
            aVtmp(i,j) = aV(i,j)
          else
            aUtmp(i,j) = 0.
            aVtmp(i,j) = 0.
          endif
        enddo
      enddo
      if(aJ_0==1) then
        aUtmp(:,1) = 0.
        aVtmp(:,1) = 0.
      endif
      if(aJ_1==aJM) then
        aUnp = aU(1,aJM)
        aVnp = aV(1,aJM)
        aUtmp(:,aJM) = aUnp*aCOSI(:) + aVnp*aSINI(:)
        aVtmp(:,aJM) = aVnp*aCOSI(:) - aUnp*aSINI(:)
      endif
c      call INT_AG2OG(aUtmp,oU(:,:),aWEIGHT,copypole=.false.)
c      call INT_AG2OG(aVtmp,oV(:,:),aWEIGHT,copypole=.false.)
      call INT_AG2OG_2Da_par(aUtmp,oU,aWEIGHT,copypole=.false.)
      call INT_AG2OG_2Da_par(aVtmp,oV,aWEIGHT,copypole=.false.)
      if(hasNorthPole(oGRID)) then
        oUnp = SUM(oU(:,oJM)*oCOSI(:))*2/oIM
        oVnp = SUM(oU(:,oJM)*oSINI(:))*2/oIM
        oU(1,oJM) = oUnp
        oV(1,oJM) = oVnp
      endif

      DEALLOCATE(aUtmp,aVtmp)

      RETURN
      END SUBROUTINE INT_AG2OG_Vector1_par

!      SUBROUTINE INT_AG2OG_3Db_par(aA,oA,aWEIGHT, aN,oN)
!!@sum INT_AG2OG_3Db_par parallel version of INT_AG2OG_3Db
!!@auth M. Kelley
!      USE RESOLUTION, only : aIM=>im,aJM=>jm
!      USE OCEAN,      only : oIM=>im,oJM=>jm
!      USE DOMAIN_DECOMP_ATM, only : aGRID=>GRID
!      Use OCEANR_DIM,       only : oGRID
!      IMPLICIT NONE
!      integer, intent(in) :: aN,oN
!      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN),
!     &     intent(in)  :: aA
!      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
!     &     intent(in)  :: aWEIGHT
!      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN),
!     &     intent(out)  :: oA
!      call INT_AG2OG(aA(:,:,oN),oA(:,:,oN),aWEIGHT)
!      RETURN
!      END SUBROUTINE INT_AG2OG_3Db_par
!
!      SUBROUTINE INT_AG2OG_3Dc_par(aA,oA,aWEIGHT, aN,oN,oNin)
!!@sum INT_AG2OG_3Dc_par parallel version of INT_AG2OG_3Dc
!!@auth M. Kelley
!      USE RESOLUTION, only : aIM=>im,aJM=>jm
!      USE OCEAN,      only : oIM=>im,oJM=>jm
!      USE DOMAIN_DECOMP_ATM, only : aGRID=>GRID
!      Use OCEANR_DIM,       only : oGRID
!      IMPLICIT NONE
!      integer, intent(in) :: aN,oN,oNin
!      real*8, dimension(aN,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
!     &     intent(in)  :: aA
!      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
!     &     intent(in)  :: aWEIGHT
!      real*8, dimension(oN,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
!     &     intent(out)  :: oA
!      real*8, dimension(:,:), allocatable :: aA2d,oA2d
!
!      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
!        oA(oNin,:aIM,:) = aA(oNin,:,:)
!      else
!        allocate(aA2d(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
!        allocate(oA2d(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
!        aA2d = aA(oNin,:,:)
!        call INT_AG2OG(aA2d,oA2d,aWEIGHT)
!        oA(oNin,:,:) = oA2d
!        deallocate(aA2d,oA2d)
!      endif
!      RETURN
!      END SUBROUTINE INT_AG2OG_3Dc_par

      END MODULE INT_AG2OG_MOD

      MODULE INT_OG2AG_MOD

!@sum INT_OG2AG_MOD contains subroutines for conversion 2D, 3D, etc. 
!!    arrays from ocean to the atm. grid 
!@auth Larissa Nazarenko
      use hntrp_mod
      IMPLICIT NONE
      SAVE
      PRIVATE
      PUBLIC INT_OG2AG

      Interface INT_OG2AG
c _par are parallelized versions
      Module Procedure INT_OG2AG_2Da_par
      Module Procedure INT_OG2AG_3Da_par
      !Module Procedure INT_OG2AG_3Db_par

      End Interface

      contains

      SUBROUTINE INT_OG2AG_2Da_par(oA,aA,oWEIGHT, CopyPole, AvgPole)
!@sum INT_OG2AG_2Da_par parallel version of INT_OG2AG_2Da
!@auth M. Kelley

      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE DOMAIN_DECOMP_1D, only : hasNorthPole, hasSouthPole
      Use OCEANR_DIM,       only : oGRID
      USE OCEAN, only : hntrp_o2a => remap_o2a
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: CopyPole
      LOGICAL, INTENT(IN), OPTIONAL :: AvgPole

      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(in)  :: oWEIGHT, oA
      real*8, dimension(hntrp_o2a%imb,
     &     hntrp_o2a%J1B_HALO:hntrp_o2a%JNB_HALO),
     &     intent(out)  :: aA

      integer :: jmin,jmax
      real*8, dimension(:,:), allocatable :: oA_band,oWEIGHT_band
      logical :: AvgPole_
      integer :: aIM,aJM

      aIM = hntrp_o2a%imb
      aJM = hntrp_o2a%jmb_full

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        aA(:,:) = oA(:aIM,:)
        return
      endif

      AvgPole_ = .true.
      if(present(AvgPole)) AvgPole_ = AvgPole

      jmin = hntrp_o2a%bpack%jband_strt
      jmax = hntrp_o2a%bpack%jband_stop
      ALLOCATE(oA_band(oIM,jmin:jmax),oWEIGHT_band(oIM,jmin:jmax))

C***  Gather the requisite ocean latitude bands
      CALL BAND_PACK (hntrp_o2a%bpack, oA, oA_band)
      CALL BAND_PACK (hntrp_o2a%bpack, oWEIGHT, oWEIGHT_band)

C***  Interpolate oA_band from ocean grid to atmospheric grid 

      if(hasNorthPole(oGRID)) then
        if(CopyPole) oWEIGHT_band(2:oIM,oJM) = oWEIGHT_band(1,oJM)
        if(AvgPole_) oA_band(2:oIM,oJM) = oA_band(1,oJM)
      endif
      if(AvgPole_) then
        call HNTR8P_band(oWEIGHT_band, oA_band, hntrp_o2a, aA)
      else
        call HNTR8_band(oWEIGHT_band, oA_band, hntrp_o2a, aA)
      endif
      DEALLOCATE(oA_band, oWEIGHT_band)

      RETURN
      END SUBROUTINE INT_OG2AG_2Da_par

      SUBROUTINE INT_OG2AG_3Da_par(oA,aA,oWEIGHT,NT,aFOCEAN)
!@sum INT_OG2AG_3Da_par parallel version of INT_OG2AG_3Da
!@auth M. Kelley
      USE OCEAN,      only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_1D, only : hasNorthPole, hasSouthPole
      Use OCEANR_DIM,       only : oGRID
      USE OCEAN, only : hntrp_o2a => remap_o2a
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NT

      real*8, dimension(NT,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(in)  :: oA
      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(in)  :: oWEIGHT
      real*8, dimension(NT,hntrp_o2a%imb,
     &     hntrp_o2a%J1B_HALO:hntrp_o2a%JNB_HALO),
     &     intent(out) :: aA
      real*8, dimension(hntrp_o2a%imb,
     &     hntrp_o2a%J1B_HALO:hntrp_o2a%JNB_HALO),
     &     intent(in) :: aFOCEAN

      real*8, dimension(:,:,:), allocatable :: oA_band
      real*8, dimension(:,:), allocatable :: oWEIGHT_band,oA2d,aA2d
      integer :: i,j,n, jmin,jmax
      integer :: aIM,aJM,imaxj

      aIM = hntrp_o2a%imb
      aJM = hntrp_o2a%jmb_full

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        DO J=hntrp_o2a%J1B,hntrp_o2a%JNB
          if(j==1 .or. j==aJM) then
            imaxj = 1
          else
            imaxj = aIM
          endif
          DO I=1,imaxj
            IF (aFOCEAN(I,J).gt.0.) THEN
              DO N=1,NT
                aA(N,I,J) = oA(N,I,J)
              END DO
            END IF
          END DO
        END DO
        return
      endif

      allocate(aA2d(aIM,hntrp_o2a%J1B_HALO:hntrp_o2a%JNB_HALO))

      jmin = hntrp_o2a%bpack%jband_strt
      jmax = hntrp_o2a%bpack%jband_stop
      ALLOCATE(oA_band(NT,oIM,jmin:jmax),oWEIGHT_band(oIM,jmin:jmax))
      allocate(oA2d(oIM,jmin:jmax))

C***  Gather the requisite ocean latitude bands
      CALL BAND_PACK_COLUMN (hntrp_o2a%bpack, oA, oA_band)
      CALL BAND_PACK (hntrp_o2a%bpack, oWEIGHT, oWEIGHT_band)

C***  Interpolate oA from ocean grid to atmospheric grid 
      do n=1,NT
        oA2d(:,:) = oA_band(n,:,:)
        if(hasNorthPole(oGRID)) then
          oA2D(2:oIM,oJM) = oA2D(1,oJM)
        endif
        call HNTR8P_band(oWEIGHT_band, oA2D, hntrp_o2a, aA2d)
        do j=hntrp_o2a%J1B,hntrp_o2a%JNB
          if(j==1 .or. j==aJM) then
            imaxj = 1
          else
            imaxj = aIM
          endif
          do i=1,imaxj
            if (afocean(i,j).gt.0.) then
              aA(n,i,j) = aA2d(i,j)
            endif
          enddo
        enddo
      enddo

      DEALLOCATE(oA_band, oWEIGHT_band, oA2d,aA2d)

      RETURN
      END SUBROUTINE INT_OG2AG_3Da_par

!      SUBROUTINE INT_OG2AG_3Db_par(oA,aA,oWEIGHT,  oN,aN, CopyPole)
!!@sum INT_OG2AG_3Db_par parallel version of INT_OG2AG_3Db
!!@auth M. Kelley
!      USE OCEAN,      only : oIM=>im,oJM=>jm
!
!      USE DOMAIN_DECOMP_1D, only : hasNorthPole, hasSouthPole
!      Use OCEANR_DIM,       only : oGRID
!      USE OCEAN, only : hntrp_o2a => remap_o2a
!      IMPLICIT NONE
!
!      INTEGER, INTENT(IN) :: oN,aN
!      LOGICAL, INTENT(IN) :: CopyPole
!
!      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN),
!     &     intent(in)  :: oA
!      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
!     &     intent(in)  :: oWEIGHT
!      real*8, dimension(hntrp_o2a%imb,
!     &     hntrp_o2a%J1B_HALO:hntrp_o2a%JNB_HALO,aN),
!     &     intent(out)  :: aA
!
!      real*8, dimension(:,:,:), allocatable :: oA_band
!      real*8, dimension(:,:), allocatable :: oWEIGHT_band
!      integer :: n, jmin,jmax
!      integer :: aIM,aJM
!
!      aIM = hntrp_o2a%imb
!      aJM = hntrp_o2a%jmb_full
!
!      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
!        aA(:,:,1:aN) = oA(:aIM,:,1:aN)
!        return
!      endif
!
!      jmin = hntrp_o2a%bpack%jband_strt
!      jmax = hntrp_o2a%bpack%jband_stop
!      ALLOCATE(oA_band(oIM,jmin:jmax,aN),oWEIGHT_band(oIM,jmin:jmax))
!
!C***  Gather the requisite ocean latitude bands
!      CALL BAND_PACK (hntrp_o2a%bpack, oA(:,:,1:aN), oA_band)
!      CALL BAND_PACK (hntrp_o2a%bpack, oWEIGHT, oWEIGHT_band)
!      if(CopyPole .and. hasNorthPole(oGRID)) then
!        oWEIGHT_band(2:oIM,oJM) = oWEIGHT_band(1,oJM)
!      endif
!
!C***  Interpolate oA from ocean grid to atmospheric grid 
!      do n=1,aN
!        if(hasNorthPole(oGRID)) then
!          oA_band(2:oIM,oJM,n) = oA_band(1,oJM,n)
!        endif
!        call HNTR8P_band(oWEIGHT_band,oA_band(:,:,n),hntrp_o2a,
!     &       aA(:,:,n))
!      enddo
!
!      DEALLOCATE(oA_band, oWEIGHT_band)
!
!      RETURN
!      END SUBROUTINE INT_OG2AG_3Db_par

      END MODULE INT_OG2AG_MOD

      Subroutine HNTR80 (IMA,JMA,OFFIA,DLATA,
     *                   IMB,JMB,OFFIB,DLATB, DATMIS)
C****
C**** HNTR80 fills in the common block HNTRCB with coordinate
C**** parameters that will be used by subsequent calls to HNTR8.
C**** The 5 Real input values are expected to be Real*8.
C****
C**** Input: IMA = number of cells in east-west direction of grid A
C****        JMA = number of cells in north-south direction of grid A
C****      OFFIA = number of cells of grid A in east-west direction
C****              from IDL (180) to western edge of cell IA=1
C****      DLATA = minutes of latitude for non-polar cells on grid A
C****        IMB = number of cells in east-west direction of grid B
C****        JMB = number of cells in north-south direction of grid B
C****      OFFIB = number of cells of grid B in east-west direction
C****              from IDL (180) to western edge of cell IB=1
C****      DLATB = minutes of latitude for non-polar cells on grid B
C****     DATMIS = missing data value inserted in output array B when
C****              cell (IB,JB) has integrated value 0 of WTA
C****
C**** Output: common block /HNTRCB/
C**** SINA(JA) = sine of latitude of northern edge of cell JA on grid A
C**** SINB(JB) = sine of latitude of northern edge of cell JB on grid B
C**** FMIN(IB) = fraction of cell IMIN(IB) on grid A west of cell IB
C**** FMAX(IB) = fraction of cell IMAX(IB) on grid A east of cell IB
C**** GMIN(JB) = fraction of cell JMIN(JB) on grid A south of cell JB
C**** GMAX(JB) = fraction of cell JMAX(JB) on grid A north of cell JB
C**** IMIN(IB) = western most cell of grid A that intersects cell IB
C**** IMAX(IB) = eastern most cell of grid A that intersects cell IB
C**** JMIN(JB) = southern most cell of grid A that intersects cell JB
C**** JMAX(JB) = northern most cell of grid A that intersects cell JB
C****
      Implicit None
      Integer,Intent(In) :: IMA,JMA,IMB,JMB
      Real*8, Intent(In) :: OFFIA,DLATA,OFFIB,DLATB,DATMIS
      Real*8,Parameter :: TWOPI=6.283185307179586477d0
      Real*8  :: FMIN,FMAX,GMIN,GMAX, SINA,SINB, DATMCB
      Integer :: IMIN,IMAX,JMIN,JMAX, INA,JNA,INB,JNB
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMCB, INA,JNA, INB,JNB
      Integer :: IA,IB,JA,JB, IBp1
      Real*8  :: DIA,DIB, RIA,RIB,RJA,RJB, FJEQA,FJEQB
C****
      INA = IMA  ;  JNA = JMA
      INB = IMB  ;  JNB = JMB
      DATMCB = DATMIS
      If (IMA<1 .or. IMA>10800 .or. JMA<1 .or. JMA>5401 .or.
     *    IMB<1 .or. IMB>10800 .or. JMB<1 .or. JMB>5401)  GoTo 400
C****
C**** Partitions in east-west (I) direction
C**** Domain, around the globe, is scaled to fit from 0 to IMA*IMB
C****
      DIA = IMB  !  width of single A grid cell in scaled domain
      DIB = IMA  !  width of single B grid cell in scaled domain
      IA  = 1
      RIA = (IA+OFFIA - IMA)*IMB  !  scaled longitude of eastern edge
      IB  = IMB
      Do 150 IBp1=1,IMB
      RIB = (IBp1-1+OFFIB)*IMA    !  scaled longitude of eastern edge
  110 If (RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GoTo 110
C**** Eastern edges of cells IA of grid A and IB of grid B coincide
  130 IMAX(IB) = IA
      FMAX(IB) = 0
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 0
      GoTo 150
C**** Cell IA of grid A contains western edge of cell IB of grid B
  140 IMAX(IB) = IA
      FMAX(IB) = (RIA-RIB)/DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 1-FMAX(IB)
  150 IB = IBp1
      IMAX(IMB) = IMAX(IMB) + IMA
C       WRITE (0,915) 'IMIN=',IMIN(1:IMB)
C       WRITE (0,915) 'IMAX=',IMAX(1:IMB)
C       WRITE (0,916) 'FMIN=',FMIN(1:IMB)
C       WRITE (0,916) 'FMAX=',FMAX(1:IMB)
C****
C**** Partitions in the north-south (J) direction
C**** Domain is measured in minutes (1/60-th of a degree)
C****
      FJEQA = .5*(1+JMA)
      Do 210 JA=1,JMA-1
      RJA = (JA+.5-FJEQA)*DLATA  !  latitude in minutes of northern edge
  210 SINA(JA) = Sin (RJA*TWOPI/(360*60))
      SINA(0)  = -1
      SINA(JMA)=  1
C****
      FJEQB = .5*(1+JMB)
      Do 220 JB=1,JMB-1
      RJB = (JB+.5-FJEQB)*DLATB  !  latitude in minutes of northern edge
  220 SINB(JB) = Sin (RJB*TWOPI/(360*60))
      SINB(0)  = -1
      SINB(JMB)=  1
C****
      JMIN(1) = 1
      GMIN(1) = 0
      JA = 1
      Do 350 JB=1,JMB-1
  310 If (SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GoTo 310
C**** Northern edges of cells JA of grid A and JB of grid B coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0
      GoTo 350
C**** Cell JA of grid A contains northern edge of cell JB of grid B
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA) - SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB) - SINA(JA-1)
  350 Continue
      JMAX(JMB) = JMA
      GMAX(JMB) = 0
C       WRITE (0,915) 'JMIN=',JMIN(1:JMB)
C       WRITE (0,915) 'JMAX=',JMAX(1:JMB)
C       WRITE (0,916) 'GMIN=',GMIN(1:JMB)
C       WRITE (0,916) 'GMAX=',GMAX(1:JMB)
      Return
C****
C**** Invalid parameters or dimensions out of range
C****
  400 Write (0,940) IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,DLATB, DATMIS
      Stop 400
C****
C 915 Format (/ 1X,A5 / (20I6))
C 916 Format (/ 1X,A5 / (20F6.2))
  940 Format ('0Arguments received by HNTRP0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATA   = minutes of latitude for interior grid cell'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATB   = minute of latitude for interior grid cell'/
     *  E24.8,' = DATMIS  = missing data value to be put in B array',
     *                    ' when integrated WTA = 0'/
     *  '0These arguments are invalid or out of range.')
      End Subroutine HNTR80


      Subroutine HNTR8 (WTA,A,B)
C****
C**** HNTR8 performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value DATMIS.
C**** The area weighted integral of the quantity is conserved.
C**** The 3 Real input values are expected to be Real*8.
C****
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
      Implicit None
      Real*8,Intent(In)  :: A(*),WTA(*)
      Real*8,Intent(Out) :: B(*)
      Real*8  :: FMIN,FMAX,GMIN,GMAX, SINA,SINB, DATMIS
      Integer :: IMIN,IMAX,JMIN,JMAX, IMA,JMA,IMB,JMB
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
      Integer :: IA,IB,JA,JB, IJA,IJB, IAREV, IAMIN,IAMAX,JAMIN,JAMAX
      Real*8  :: WEIGHT,VALUE, F,G
C****
C**** Interpolate the A grid onto the B grid
C****
      Do 20 JB=1,JMB
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      Do 20 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      WEIGHT= 0
      VALUE = 0
      IAMIN = IMIN(IB)
      IAMAX = IMAX(IB)
      Do 10 JA=JAMIN,JAMAX
      G = SINA(JA)-SINA(JA-1)
      If (JA==JAMIN)  G = G - GMIN(JB)
      If (JA==JAMAX)  G = G - GMAX(JB)
      Do 10 IAREV=IAMIN,IAMAX
      IA  = 1 + Mod(IAREV-1,IMA)
      IJA = IA + IMA*(JA-1)
      F   = 1
      If (IAREV==IAMIN)  F = F - FMIN(IB)
      If (IAREV==IAMAX)  F = F - FMAX(IB)
      WEIGHT = WEIGHT + F*G*WTA(IJA)
   10 VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = DATMIS
      If (WEIGHT.ne.0)  B(IJB) = VALUE/WEIGHT
   20 Continue
      Return
      End Subroutine HNTR8


      Subroutine HNTR8P (WTA,A,B)
C****
C**** HNTR8P is similar to HNTR8 but polar values are replaced by
C**** their longitudinal mean.
C**** The 3 Real input values are expected to be Real*8.
C****
      Implicit None
      Real*8,Intent(In)  :: A(*),WTA(*)
      Real*8,Intent(Out) :: B(*)
      Real*8  :: FMIN,FMAX,GMIN,GMAX, SINA,SINB, DATMIS
      Integer :: IMIN,IMAX,JMIN,JMAX, IMA,JMA,IMB,JMB
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
      Integer :: IB,JB, IJB
      Real*8  :: BMEAN,WEIGHT,VALUE
C****
      Call HNTR8 (WTA,A,B)
C****
C**** Replace individual values near the poles by longitudinal mean
C****
      Do 40 JB=1,JMB,JMB-1
      BMEAN  = DATMIS
      WEIGHT = 0
      VALUE  = 0
      Do 10 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      If (B(IJB) == DATMIS)  GoTo 10
      WEIGHT = WEIGHT + 1
      VALUE  = VALUE  + B(IJB)
   10 Continue
      If (WEIGHT.ne.0)  BMEAN = VALUE/WEIGHT
   20 Do 30 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
   30 B(IJB) = BMEAN
   40 Continue
      Return
      End Subroutine HNTR8P

