#include "rundeck_opts.h"

      MODULE ICEDYN
      USE OCEAN,  only     : oIM=>im,oJM=>jm
      use DOMAIN_DECOMP_1D, only : DIST_GRID
      IMPLICIT NONE 
      SAVE
      INTEGER, parameter :: IMICDYN = oIM, JMICDYN = oJM !     resolution of icedyn grid set to ocean resolution
      TYPE(DIST_GRID) :: grid_ICDYN
      REAL*8  :: DLATM
      END MODULE ICEDYN

      MODULE ICEDYN_COM
      IMPLICIT NONE
      SAVE

      integer :: imic=0,kticij=0
      real*8,dimension(2) :: rsix=0,rsiy=0,usi=0,vsi=0
      real*8,dimension(2) :: USIDT=0,VSIDT=0,RSISAVE=0
      real*8,dimension(2) :: icij=0,ticij=0
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: UOSURF,VOSURF,DMUI,DMVI
      END MODULE ICEDYN_COM

      subroutine alloc_icedyn
      USE DOMAIN_DECOMP_1D, ONLY : INIT_GRID
      USE ICEDYN, only : grid_ICDYN,IMICDYN,JMICDYN
      LOGICAL, SAVE :: init = .false.

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

      CALL INIT_GRID(grid_ICDYN,IMICDYN,JMICDYN,1)

      return
      end subroutine alloc_icedyn

      subroutine GEOMICDYN
      USE ICEDYN, only : DLATM,IMICDYN,JMICDYN
      IMPLICIT NONE
      real*8 :: DLAT_DG

      DLAT_DG=180./REAL(JMICDYN)                   ! even spacing (default)
      IF (JMICDYN.eq.46) DLAT_DG=180./REAL(JMICDYN-1)   ! 1/2 box at pole for 4x5
      IF (JMICDYN.eq.24) DLAT_DG=180./REAL(JMICDYN-1.5) ! 1/4 box at pole, 'real' 8x10
      DLATM=60.*DLAT_DG
      end subroutine GEOMICDYN

      subroutine alloc_icedyn_com(grid)
      use DOMAIN_DECOMP_ATM, only : DIST_GRID
      use DOMAIN_DECOMP_1D, only : GET
      USE ICEDYN, only : grid_icdyn
      USE ICEDYN_COM, only : UOSURF,VOSURF,DMUI,DMVI
      IMPLICIT NONE
      LOGICAL, SAVE :: init=.false.
      INTEGER :: I_1H    , I_0H
      INTEGER :: J_1H    , J_0H
      TYPE (DIST_GRID), INTENT(IN) :: grid

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

      call getDomainBounds(grid_icdyn, 
     &     I_STRT_HALO=I_0H    , I_STOP_HALO=I_1H    ,
     &     J_STRT_HALO=J_0H    , J_STOP_HALO=J_1H    )

      ALLOCATE(DMUI( I_0H:I_1H , J_0H:J_1H ),
     &         DMVI( I_0H:I_1H , J_0H:J_1H ))

      ALLOCATE(UOSURF( I_0H:I_1H , J_0H:J_1H ),
     &         VOSURF( I_0H:I_1H , J_0H:J_1H ))
      UOSURF = 0. ! in case there is no dynamic ocean
      VOSURF = 0. ! ""

      return
      end subroutine alloc_icedyn_com

      SUBROUTINE ICEDYN_DUM
!@sum ICEDYN_DUM dummy routines to replace ice dynamics
c      ENTRY alloc_icedyn
c      ENTRY alloc_icedyn_com
      ENTRY gather_icdiags
      ENTRY io_icedyn
      ENTRY io_icdiag
      ENTRY reset_icdiag
      ENTRY ADVSI
      ENTRY init_icedyn
      ENTRY diag_ICEDYN
      entry def_rsf_icedyn
      entry new_io_icedyn
      entry def_rsf_icdiag
      entry new_io_icdiag
      entry set_ioptrs_iceacc_default
      entry def_meta_icdiag
      entry write_meta_icdiag
      RETURN
      END SUBROUTINE ICEDYN_DUM


      SUBROUTINE DYNSI
!@sum DYNSI simple coding to estimate ice-ocean friction velocity
!@auth Gavin Schmidt
      USE CONSTANT, only : rhows
      USE MODEL_COM, only : im,jm,kocean,focean,dtsrc
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds
      USE GEOM, only : imaxj
      USE SEAICE, only : oi_ustar0
      USE SEAICE_COM, only : rsi
      USE FLUXES, only : UI2rho,dmua,dmva

      IMPLICIT NONE
      INTEGER I,J
      REAL*8 ustar1
      INTEGER :: I_0,I_1, J_0,J_1


      IF (KOCEAN.eq.1) THEN

        call getDomainBounds(grid, J_STRT=J_0,   J_STOP=J_1 )
        I_0 = GRID%I_STRT
        I_1 = GRID%I_STOP
        DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
c          UI2rho(I,J) = rhows*(oi_ustar0)**2  ! default
C**** with wind stress dependence
          if (rsi(i,j)*focean(i,j).gt.0) then
            ustar1= SQRT(SQRT(DMUA(I,J,2)**2+DMVA(I,J,2)**2)/
     &           (DTSRC*RHOWS))
            UI2rho(I,J)=rhows*(oi_ustar0*max(1d0,1d3*ustar1))**2
          end if
        END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE DYNSI

