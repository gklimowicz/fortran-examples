
C PLEASE KEEP THIS NOTE OF MODEL-DEVELOPMENT HISTORY
C Matrix solve uses Thomas algorithm, 10/1991, Jinlun Zhang
C Spherical coordinate system, 10/27/93, Jinlun Zhang
C Latest finite differencing scheme for treatment of NP,
C 9/9/1996,Jinlun Zhang
C Alternating direction implicit (ADI) method is used, 10/1998,
C Jinlun Zhang
C For details about ADI dynamics model, see Zhang and Rothrock,
C "Modeling
C Arctic Sea ice with an efficient plastic solution",
C submitted to JGR, 1999
C Adapted for GISS coupled model Jiping Liu/Gavin Schmidt 2000
C - Further modularised May 2001
C*************************************************************
#include "rundeck_opts.h"
      MODULE ICEDYN
!@sum  ICEDYN holds local variables for dynamic sea ice
!@auth Gavin Schmidt (based on code from Jinlun Zhang)
      USE CONSTANT, only : radian,radius,TWOPI
      USE DOMAIN_DECOMP_1D, only : DIST_GRID
      USE SEAICE, only : osurf_tilt
      IMPLICIT NONE
      SAVE

C**** Dimensions of ice dynamics grids grid_NXY and grid_ICDYN
C**** Edit the definition of nx1,ny1 to change the grid for the
C**** rheology calculations without changing ADVSI grid.
!@var nx1 number of grid points in the longitudinal direction
!@+  (calculated points from 2 through nx1-1. End points are boundaries)
!@var ny1 number of grid points in the latitudinal direction
!@+  (calculated points from 2 through ny1-1. End points are boundaries)
!#ifdef CUBED_SPHERE
!      INTEGER, parameter :: IMICDYN = 2*IM, JMICDYN = 2*JM
!#else
!      INTEGER, parameter :: IMICDYN = IM, JMICDYN = JM
!#endif
      INTEGER :: IMICDYN, JMICDYN
!     dimensions of grid_NXY:
      !INTEGER, parameter :: NX1=IMICDYN+2, NY1=JMICDYN
      INTEGER :: NX1, NY1
      !INTEGER, parameter :: NYPOLE=NY1-1,NXLCYC=NX1-1
      INTEGER :: NYPOLE,NXLCYC
      integer :: NPOL=1,LCYC=1

      TYPE(DIST_GRID) :: grid_NXY
      TYPE(DIST_GRID), target :: grid_ICDYN

!@var FOCEAN land/ocean mask on ice dynamic grid
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: FOCEAN

C**** input
!@var HEFFM ice mass mask (1/0)
!@var UVM ice velocity mask (1/0)
!@var COR coriolis term for ice dynamic equation
!@var GAIRX,GAIRY atmosphere-ice stress (B grid)
!@var GWATX,GWATY ocean velocities (B grid) (m/s)
!@var PGFUB,PGFVB pressure accelaration force
!@var AMASS ice mass (kg/m^2)
!@var UICEC,VICEC velocity arrays  (m/s)
!@var UICE,VICE velocity arrays  (m/s)
!@var HEFF ice thickness (mean over box) (m)
!@var AREA ice area (frac)
C**** internal variables
!@var PRESS ice internal pressure (Pa)
!@var FORCEX,FORCEY external force
!@var DRAGS,DRAGA symmetric/anti-symmetric drag terms
!@var ZMAX,ZMIN max,min values of ZETA
!@var ETA,ZETA viscosities
C**** output
!@var DWATN non-linear water drag term
!@var DMU,DMV ice-ocean stress
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: PRESS,HEFFM,UVM,DWATN,COR
     *     ,ZMAX,ZMIN,ETA,ZETA,DRAGS,DRAGA,GAIRX,GAIRY
     *     ,GWATX,GWATY,PGFUB,PGFVB,FORCEX,FORCEY,AMASS,UICEC,VICEC
     *     ,DMU,DMV,HEFF,AREA,USI,VSI
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UICE,VICE


C**** Geometry
!@var SINEN sin(phi)
!@var BYDXDY
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: SINEN,BYDXDY
!@var DXT,DXU x-direction distances on tracer and velocity grid
!@var DYT,DYU y-direction distances on tracer and velocity grid
!      REAL*8, DIMENSION(NX1) :: DXT,DXU,BYDX2,BYDXR
!      REAL*8, DIMENSION(NY1) :: DYT,DYU,BYDY2,BYDYR,CST,
!     &                          CSU,TNGT,TNG,BYCSU

      REAL*8, DIMENSION(:), ALLOCATABLE ::
     &     DXT,DXU,BYDX2,BYDXR,
     &     DYT,DYU,BYDY2,BYDYR,CST,
     &     CSU,TNGT,TNG,BYCSU

!@var OIPHI ice-ocean turning angle (25 degrees)
!@var ECCEN value of eccentricity for yield curve ellipse
      REAL*8, PARAMETER :: ECCEN=2.0, OIPHI=25d0*radian
!@var SINWAT,COSWAT sin and cos of ice-ocean turning angle
      REAL*8 SINWAT,COSWAT

C**** Geometry inherited from geomb
      REAL*8, DIMENSION(:), ALLOCATABLE ::
     &     LAT,LATB,LON,LONB,DXYP,BYDXYP,DXYV,DXP,DYP,DXV,DYV,DXYS,DXYN,
     &     SINIU,COSIU
      REAL*8, DIMENSION(:,:), ALLOCATABLE ::
     &     LON_DG,LAT_DG

!@param  DLON grid spacing in longitude (deg)
      !REAL*8, PARAMETER :: FIM=IMICDYN, BYIM=1./FIM
      !REAL*8, PARAMETER :: DLON = TWOPI*BYIM
!@var DLAT,DLAT_DG,DLATM grid spacing in latitude (rad,deg,minutes)
      REAL*8  :: DLON,DLAT,DLON_DG,DLAT_DG,DLATM
      !REAL*8 :: FIM,BYIM

!!@param FJEQ equatorial value of J
!      !REAL*8, PARAMETER :: FJEQ=.5*(1+JMICDYN)
!!@var  LAT latitude of mid point of primary grid box (radians)
!      REAL*8, DIMENSION(JMICDYN) :: LAT,LATB
!!@var  LAT_DG latitude of mid points of primary and sec. grid boxs (deg)
!      REAL*8, DIMENSION(JMICDYN,2) :: LAT_DG  ! for diags
!!@var  LON longitude of mid points of primary grid box (radians)
!      REAL*8, DIMENSION(IMICDYN) :: LON,LONB
!!@var  LON_DG longitude of mid points of prim. and sec. grid boxes (deg)
!      REAL*8, DIMENSION(IMICDYN,2) :: LON_DG  ! for diags
!!@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
!C**** Note that this is not the exact area, but is what is required for
!C**** some B-grid conservation quantities
!      REAL*8, DIMENSION(JMICDYN) :: DXYP,BYDXYP
!      !REAL*8, DIMENSION(IMICDYN,JMICDYN) :: aDXYP
!!@var DXYV,BYDXYV area of grid box around velocity point (recip.)(m^2)
!      REAL*8, DIMENSION(JMICDYN) :: DXYV !,BYDXYV
!!@var  DXP,DYP,BYDXP,BYDYP distance between points on primary grid
!!@+     (+inverse)
!      REAL*8, DIMENSION(JMICDYN) :: DXP,DYP !,BYDXP,BYDYP
!!@var  DXV,DYV distance between velocity points (secondary grid)
!      REAL*8, DIMENSION(JMICDYN) :: DXV,DYV
!!@var  DXYN,DXYS half box areas to the North,South of primary grid point
!      REAL*8, DIMENSION(JMICDYN) :: DXYS,DXYN
!!@var  SINP sin of latitude at primary grid points
!      !REAL*8, DIMENSION(JMICDYN) :: SINP
!!@var  COSP, COSV cos of latitude at primary, secondary latitudes
!      !REAL*8, DIMENSION(JMICDYN) :: COSP,COSV
!!@var  IMAXJ varying number of used longitudes
!      !INTEGER, DIMENSION(JMICDYN) :: IMAXJ

!@var PSTAR maximum sea ice pressure per unit thickness (Pa/m)
      REAL*8, PARAMETER :: PSTAR=2.75d4

!@var BYDTS reciprocal of timestep in ice dynamics code
      REAL*8 :: BYDTS

      CONTAINS

      SUBROUTINE FORM
!@sum  FORM calculates ice dynamics input parameters for relaxation
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, NORTH,SOUTH
      USE DOMAIN_DECOMP_1D, ONLY : HALO_UPDATE,
     &     hasSouthPole, hasNorthPole
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 AAA
      INTEGER :: J_0,J_1,J_0S,J_1S

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid_ICDYN, J_STRT=J_0, J_STOP=J_1,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S )

C****
C**** Set up non linear water drag
C****
      DO J=J_0,J_1
      DO I=1,NX1-1
        DWATN(I,J)=5.5*SQRT((UICE(I,J,1)-GWATX(I,J))**2
     1       +(VICE(I,J,1)-GWATY(I,J))**2)
      END DO
      END DO
C NOW SET UP SYMMETRIC DRAG
      DO J=J_0,J_1S
      DO I=1,NX1-1
        DRAGS(I,J)=DWATN(I,J)*COSWAT
      END DO
      END DO
C NOW SET UP ANTI SYMMETRIC DRAG PLUS CORIOLIS
      DO J=J_0,J_1S
      DO I=1,NX1-1
        IF(J.GT.NY1/2) THEN
          DRAGA(I,J)=DWATN(I,J)*SINWAT+COR(I,J)
        ELSE
          DRAGA(I,J)=DWATN(I,J)*(-SINWAT)+COR(I,J)
        END IF
      END DO
      END DO
C NOW SET UP FORCING FIELD
      DO J=J_0,J_1S
      DO I=1,NX1-1

C FIRST DO WIND
        FORCEX(I,J)=GAIRX(i,j)
        FORCEY(I,J)=GAIRY(i,j)

C NOW ADD IN CURRENT FORCE
        IF(J.GT.NY1/2) THEN
          FORCEX(I,J)=FORCEX(I,J)+DWATN(I,J)*(COSWAT*GWATX(I,J)
     1         -SINWAT*GWATY(I,J))
          FORCEY(I,J)=FORCEY(I,J)+DWATN(I,J)*(SINWAT*GWATX(I,J)
     1         +COSWAT*GWATY(I,J))
        ELSE
          FORCEX(I,J)=FORCEX(I,J)+DWATN(I,J)*(COSWAT*GWATX(I,J)
     1         +SINWAT*GWATY(I,J))
          FORCEY(I,J)=FORCEY(I,J)+DWATN(I,J)*(-SINWAT*GWATX(I,J)
     1         +COSWAT*GWATY(I,J))
        END IF

C     NOW ADD IN TILT
        if (osurf_tilt.eq.1) then
C**** This assumes explicit knowledge of sea surface tilt
          FORCEX(I,J)=FORCEX(I,J)+AMASS(I,J)*PGFUB(I,J)
          FORCEY(I,J)=FORCEY(I,J)+AMASS(I,J)*PGFVB(I,J)
        else
C**** Otherwise estimate tilt using geostrophy
          FORCEX(I,J)=FORCEX(I,J)-COR(I,J)*GWATY(I,J)
          FORCEY(I,J)=FORCEY(I,J)+COR(I,J)*GWATX(I,J)
        end if
      END DO
      END DO

C NOW SET UP ICE PRESSURE AND VISCOSITIES
      DO J=J_0,J_1
      DO I=1,NX1
        PRESS(I,J)=PSTAR*HEFF(I,J)*EXP(-20.0*(1.0-AREA(I,J)))
        ZMAX(I,J)=(5d12/2d4)*PRESS(I,J)
c       ZMIN(I,J)=0.0D+00
        ZMIN(I,J)=4d8
      END DO
      END DO

      CALL PLAST

      if (hasNorthPole(grid_ICDYN)) then
       AAA=0.0
       DO I=2,NX1-1
         AAA=AAA+PRESS(I,NY1-1)
       END DO
       AAA=AAA/FLOAT(NX1-2)
       DO I=1,NX1
         PRESS(I,NY1)=AAA
       END DO
      end if

 8481 CONTINUE

      DO J=J_0,J_1
        PRESS(1,J)=PRESS(NX1-1,J)
        PRESS(NX1,J)=PRESS(2,J)
      END DO

C NOW SET VISCOSITIES AND PRESSURE EQUAL TO ZERO AT OUTFLOW PTS

      DO J=J_0,J_1
      DO I=1,NX1
        PRESS(I,J)=PRESS(I,J)*HEFFM(I,J)
        ETA(I,J)=ETA(I,J)*HEFFM(I,J)
        ZETA(I,J)=ZETA(I,J)*HEFFM(I,J)
      END DO
      END DO

C NOW CALCULATE PRESSURE FORCE AND ADD TO EXTERNAL FORCE
C**** Update halo of PRESS for distributed memory implementation
      CALL HALO_UPDATE(grid_ICDYN, PRESS, FROM=NORTH)
      DO J=J_0,J_1S
        DO I=1,NX1-1
          FORCEX(I,J)=FORCEX(I,J)-(0.25/(DXU(I)*CSU(J)))
     1         *(PRESS(I+1,J)+PRESS(I+1,J+1)-PRESS(I,J)-PRESS(I,J+1))
          FORCEY(I,J)=FORCEY(I,J)-0.25/DYU(J)
     1         *(PRESS(I,J+1)+PRESS(I+1,J+1)-PRESS(I,J)-PRESS(I+1,J))
C NOW PUT IN MINIMAL MASS FOR TIME STEPPING CALCULATIONS
        END DO
      END DO

      DO J=J_0,J_1
        FORCEX(1,J)=FORCEX(NX1-1,J)
        FORCEY(1,J)=FORCEY(NX1-1,J)
        FORCEX(NX1,J)=FORCEX(2,J)
        FORCEY(NX1,J)=FORCEY(2,J)
        DWATN(1,J)=DWATN(NX1-1,J)
        DWATN(NX1,J)=DWATN(2,J)
      END DO

      RETURN
      END SUBROUTINE FORM

      SUBROUTINE PLAST
!@sum  PLAST Calculates strain rates and viscosity for dynamic ice
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, NORTH,SOUTH
      USE DOMAIN_DECOMP_1D, ONLY : HALO_UPDATE,
     &     hasSouthPole, hasNorthPole
      IMPLICIT NONE
      REAL*8, DIMENSION(NX1,grid_ICDYN%j_strt_halo:
     &     grid_ICDYN%j_stop_halo) :: E11,E22,E12

      REAL*8, PARAMETER :: ECM2 = 1.0/(ECCEN**2),GMIN=1d-20
      REAL*8 DELT,DELT1,AAA
      INTEGER I,J

      INTEGER :: J_0,J_1,J_0S,J_1S

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid_ICDYN, J_STRT=J_0, J_STOP=J_1,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S )


C EVALUATE STRAIN RATES
      CALL HALO_UPDATE(grid_ICDYN, UICE, FROM=SOUTH)
      CALL HALO_UPDATE(grid_ICDYN, VICE, FROM=SOUTH)
      DO J=J_0S,J_1S
        DO I=2,NX1-1
          E11(I,J)=0.5/(DXT(I)*CST(J))*(UICE(I,J,1)+UICE(I,J-1,1)
     *         -UICE(I-1,J,1)-UICE(I-1,J-1,1))-0.25*(VICE(I,J,1)+
     *         VICE(I-1,J,1)+VICE(I-1,J-1,1)+VICE(I,J-1,1))*TNGT(J)
     *         /RADIUS
          E22(I,J)=0.5/DYT(J)*(VICE(I,J,1)+VICE(I-1,J,1)
     *         -VICE(I,J-1,1)-VICE(I-1,J-1,1))
          E12(I,J)=0.5*(0.5/DYT(J)*(UICE(I,J,1)+UICE(I-1,J,1)-
     *         UICE(I,J-1,1)-UICE(I-1,J-1,1))+0.5/(DXT(I)*CST(J))*
     *         (VICE(I,J,1)+VICE(I,J-1,1)-VICE(I-1,J,1)-VICE(I-1,J-1,1))
     *         +0.25*(UICE(I,J,1)+UICE(I-1,J,1)+UICE(I-1,J-1,1)+UICE(I,J
     *         -1,1))*TNGT(J)/RADIUS)
C NOW EVALUATE VISCOSITIES
          DELT=(E11(I,J)**2+E22(I,J)**2)*(1.0+ECM2)+4.0*ECM2*E12(I,J)**2
     *         +2.0*E11(I,J)*E22(I,J)*(1.0-ECM2)
          DELT1=SQRT(DELT)
          DELT1=MAX(GMIN,DELT1)
          ZETA(I,J)=0.5*PRESS(I,J)/DELT1
        END DO
      END DO

C NOW PUT MIN AND MAX VISCOSITIES IN
      DO J=J_0S,J_1S
        DO I=2,NX1-1
          ZETA(I,J)=MIN(ZMAX(I,J),ZETA(I,J))
          ZETA(I,J)=MAX(ZMIN(I,J),ZETA(I,J))
        END DO
      END DO

      if (hasNorthPole(grid_ICDYN)) then
        AAA=0.0
        DO I=2,NX1-1
          AAA=AAA+ZETA(I,NY1-1)
        END DO
        AAA=AAA/FLOAT(NX1-2)
        DO I=1,NX1
          ZETA(I,NY1)=AAA
        END DO
      end if

      if (hasSouthPole(grid_ICDYN)) then
        AAA=0.0
        DO I=2,NX1-1
          AAA=AAA+ZETA(I,2)
        END DO
        AAA=AAA/FLOAT(NX1-2)
        DO I=1,NX1
          ZETA(I,1)=AAA
        END DO
      end if

      DO J=J_0,J_1
        ZETA(1,J)=ZETA(NX1-1,J)
        ZETA(NX1,J)=ZETA(2,J)
      END DO

      DO J=J_0,J_1
        DO I=1,NX1
          ETA(I,J)=ECM2*ZETA(I,J)
c         E11(I,J)=E11(I,J)*HEFFM(I,J)
c         E22(I,J)=E22(I,J)*HEFFM(I,J)
c         E12(I,J)=E12(I,J)*HEFFM(I,J)
c         SS11=(ZETA(I,J)-ETA(I,J))*(E11(I,J)+E22(I,J))-PRESS(I,J)*0.5
        END DO
      END DO
      RETURN
      END SUBROUTINE PLAST

      SUBROUTINE RELAX
!@sum  RELAX calculates ice dynamics relaxation method
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, NORTH,SOUTH
      USE DOMAIN_DECOMP_1D, ONLY : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, ONLY : PACK_DATA, UNPACK_DATA, AM_I_ROOT,
     &     hasSouthPole, hasNorthPole
      USE TRIDIAG_MOD, only : TRIDIAG_new, TRIDIAG_cyclic
      IMPLICIT NONE

      REAL*8, DIMENSION(NX1,grid_ICDYN%J_STRT_HALO:
     &     grid_ICDYN%J_STOP_HALO) ::
     &         AU,BU,CU,FXY,FXY1
      REAL*8, DIMENSION(2:NX1-1,grid_ICDYN%J_STRT_HALO:
     &     grid_ICDYN%J_STOP_HALO) ::
     &         AV,BV,CV,VRT,U_tmp
      REAL*8, DIMENSION(NX1,grid_ICDYN%J_STRT_HALO:
     &     grid_ICDYN%J_STOP_HALO) ::
     &         FXYa,FXY1a
      REAL*8, DIMENSION(NX1) :: URT
      REAL*8, PARAMETER :: BYRAD2 = 1./(RADIUS*RADIUS)
      INTEGER I,J
      REAL*8 DELXY,DELXR,DELX2,DELY2,DELYR,ETAMEAN,ZETAMEAN,AA1,AA2
     *     ,AA3,AA4,AA5,AA6,AA9

      INTEGER :: J_0,J_1, J_0S,J_1S, J_0H,J_1H

C**** Replaces NYPOLE in loops.
      INTEGER :: J_NYP

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid_ICDYN, J_STRT=J_0,    J_STOP=J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S )

C****
C**** Modify if NYPOLE definition is modified.
C****
      J_NYP=J_1S

       DO J=J_0,J_1S
        DO I=1,NX1
          FORCEX(I,J)=FORCEX(I,J)*UVM(I,J)
          FORCEY(I,J)=FORCEY(I,J)*UVM(I,J)
        END DO
      END DO
C MUST UPDATE HEFF BEFORE CALLING RELAX
C FIRST SET U(2)=U(1)
       DO J=J_0,J_1S
        DO I=1,NX1
C NOW MAKE SURE BDRY PTS ARE EQUAL TO ZERO
          UICE(I,J,2)=UICE(I,J,1)
          VICE(I,J,2)=VICE(I,J,1)
          UICE(I,J,1)=UICE(I,J,3)*UVM(I,J)
          VICE(I,J,1)=VICE(I,J,3)*UVM(I,J)
        END DO
      END DO

      if (hasNorthPole(grid_ICDYN)) then
        DO I=1,NX1/2
         UICE(I,NY1,1)=-UICEC(I+(NX1-2)/2,NY1-1)
         VICE(I,NY1,1)=-VICEC(I+(NX1-2)/2,NY1-1)
         UICE(I,NY1,3)=-UICEC(I+(NX1-2)/2,NY1-1)
         VICE(I,NY1,3)=-VICEC(I+(NX1-2)/2,NY1-1)

         UICEC(I,NY1)=-UICEC(I+(NX1-2)/2,NY1-1)
         VICEC(I,NY1)=-VICEC(I+(NX1-2)/2,NY1-1)
        END DO

        DO I=NX1/2+1,NX1-1
         UICE(I,NY1,1)=-UICEC(I-(NX1-2)/2,NY1-1)
         VICE(I,NY1,1)=-VICEC(I-(NX1-2)/2,NY1-1)
         UICE(I,NY1,3)=-UICEC(I-(NX1-2)/2,NY1-1)
         VICE(I,NY1,3)=-VICEC(I-(NX1-2)/2,NY1-1)

         UICEC(I,NY1)=-UICEC(I-(NX1-2)/2,NY1-1)
         VICEC(I,NY1)=-VICEC(I-(NX1-2)/2,NY1-1)
        END DO
      end if

      DO J=J_0,J_1S
        UICE(1,J,1)=UICEC(NX1-1,J)
        VICE(1,J,1)=VICEC(NX1-1,J)
        UICE(NX1,J,1)=UICEC(2,J)
        VICE(NX1,J,1)=VICEC(2,J)
        UICE(1,J,3)=UICEC(NX1-1,J)
        VICE(1,J,3)=VICEC(NX1-1,J)
        UICE(NX1,J,3)=UICEC(2,J)
        VICE(NX1,J,3)=VICEC(2,J)

        UICEC(1,J)=UICEC(NX1-1,J)
        VICEC(1,J)=VICEC(NX1-1,J)
        UICEC(NX1,J)=UICEC(2,J)
        VICEC(NX1,J)=VICEC(2,J)
      END DO

C FIRST DO UICE
C THE FIRST HALF

C**Update halos for arrays eta,zeta,vicec,bycsu as needed in the next loop
        CALL HALO_UPDATE(grid_ICDYN, ETA, FROM=NORTH)
        CALL HALO_UPDATE(grid_ICDYN, ZETA, FROM=NORTH)
        CALL HALO_UPDATE(grid_ICDYN, VICEC, FROM=NORTH+SOUTH)
        CALL HALO_UPDATE(grid_ICDYN, BYCSU, FROM=NORTH+SOUTH)

      DO J=J_0S,J_NYP
      DO I=2,NXLCYC
      DELXY=BYDXDY(I,J)    ! 0.5/(DXU(I)*DYU(J))
      DELXR=BYDXR(I)       ! 0.5/(DXU(I)*RADIUS)
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      FXY(I,J)=DRAGA(I,J)*VICEC(I,J)+FORCEX(I,J)
     3+0.5*(ZETA(I+1,J+1)*(VICEC(I+1,J+1)+VICEC(I,J+1)
     3-VICEC(I+1,J)-VICEC(I,J))+ZETA(I+1,J)*(VICEC(I+1,J)
     3+VICEC(I,J)-VICEC(I+1,J-1)-VICEC(I,J-1))+ZETA(I,J+1)
     3*(VICEC(I,J)+VICEC(I-1,J)-VICEC(I,J+1)-VICEC(I-1,J+1))
     3+ZETA(I,J)*(VICEC(I,J-1)+VICEC(I-1,J-1)-VICEC(I,J)
     3-VICEC(I-1,J)))*DELXY*BYCSU(J)
     3
     4-0.5*(ETA(I+1,J+1)*(VICEC(I+1,J+1)+VICEC(I,J+1)
     4-VICEC(I+1,J)-VICEC(I,J))+ETA(I+1,J)*(VICEC(I+1,J)
     4+VICEC(I,J)-VICEC(I+1,J-1)-VICEC(I,J-1))+ETA(I,J+1)
     4*(VICEC(I,J)+VICEC(I-1,J)-VICEC(I,J+1)-VICEC(I-1,J+1))
     4+ETA(I,J)*(VICEC(I,J-1)+VICEC(I-1,J-1)-VICEC(I,J)
     4-VICEC(I-1,J)))*DELXY*BYCSU(J)
     4
     5+0.5*(VICEC(I+1,J)-VICEC(I-1,J))*(ETA(I,J+1)+ETA(I+1,J+1)
     5-ETA(I,J)-ETA(I+1,J))*DELXY*BYCSU(J)+0.5*ETAMEAN*((VICEC(I+1,J+1)
     5-VICEC(I-1,J+1))*BYCSU(J+1)-(VICEC(I+1,J-1)-VICEC(I-1,J-1))
     5*BYCSU(J-1))*DELXY
     5
     6-((ZETA(I+1,J+1)+ZETA(I+1,J)-ZETA(I,J)-ZETA(I,J+1))
     6+(ETA(I+1,J+1)+ETA(I+1,J)-ETA(I,J)-ETA(I,J+1)))
     6*TNG(J)*VICEC(I,J)*DELXR*BYCSU(J)
     6-(ETAMEAN+ZETAMEAN)*TNG(J)*(VICEC(I+1,J)-VICEC(I-1,J))
     6*DELXR*BYCSU(J)
     6
     7-ETAMEAN*2.0*TNG(J)*(VICEC(I+1,J)-VICEC(I-1,J))*DELXR*BYCSU(J)

      END DO
      END DO

      DO J=J_0S,J_NYP
      DO I=2,NXLCYC
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      AA1=((ETA(I+1,J)  +ZETA(I+1,J)  )*BYCSU(J)+
     *     (ETA(I+1,J+1)+ZETA(I+1,J+1))*BYCSU(J))*BYCSU(J)
      AA2=((ETA(I,J)+ZETA(I,J))*BYCSU(J)+(ETA(I,J+1)+ZETA(I,J+1))
     &     *BYCSU(J))*BYCSU(J)
      AA3=ETA(I,J+1)+ETA(I+1,J+1)
      AA4=ETA(I,J)+ETA(I+1,J)
      AA5=-(ETA(I,J+1)+ETA(I+1,J+1)-ETA(I,J)-ETA(I+1,J))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)
      AU(I,J)=-AA2*DELX2*UVM(I,J)
      BU(I,J)=((AA1+AA2)*DELX2+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CU(I,J)=-AA1*DELX2*UVM(I,J)
      END DO
      END DO

C**Update halos for UICE and TNG as needed for loop 1200
C**(ETA and ZETA were updted above)
      CALL HALO_UPDATE(grid_ICDYN, UICE, FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid_ICDYN, TNG, FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid_ICDYN, UICEC, FROM=SOUTH+NORTH)

      DO 1200 J=J_0S,J_NYP
      DO I=2,NXLCYC

      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))

      AA1=((ETA(I+1,J)  +ZETA(I+1,J)  )*BYCSU(J)+
     *     (ETA(I+1,J+1)+ZETA(I+1,J+1))*BYCSU(J))*BYCSU(J)
      AA2=((ETA(I,J)+ZETA(I,J))*BYCSU(J)+(ETA(I,J+1)+ZETA(I,J+1))
     &     *BYCSU(J))*BYCSU(J)
      AA3=ETA(I,J+1)+ETA(I+1,J+1)
      AA4=ETA(I,J)+ETA(I+1,J)
      AA5=-(ETA(I,J+1)+ETA(I+1,J+1)-ETA(I,J)-ETA(I+1,J))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      URT(I)=FXY(I,J)-AA5*DELYR*UICEC(I,J)
     1-(AA3+AA4)*DELY2*UICEC(I,J)
     1+(ETA(I,J+1)+ETA(I+1,J+1))*UICEC(I,J+1)*DELY2
     2+(ETA(I,J)+ETA(I+1,J))*UICEC(I,J-1)*DELY2
     3+ETAMEAN*DELYR*(UICEC(I,J+1)*TNG(J+1)-UICEC(I,J-1)*TNG(J-1))
     4-ETAMEAN*DELYR*2.0*TNG(J)*(UICEC(I,J+1)-UICEC(I,J-1))

      URT(I)=(URT(I)+AMASS(I,J)*BYDTS*UICE(I,J,2)*2.0)*UVM(I,J)

      END DO

      CALL TRIDIAG_cyclic(AU(2,J),BU(2,J),CU(2,J),URT(2),UICE(2,J,1),
     &               NXLCYC-1)
      UICE(1,J,1)   = UICE(NX1-1,J,1)
      UICE(NX1,J,1) = UICE(2,J,1)

 1200 CONTINUE

      DO I=2,NXLCYC
      DO J=J_0S,J_NYP
      UICE(I,J,3)=UICE(I,J,1)
      END DO
      END DO

C**ETA, ZETA, and TNG halos already updated above**
C NOW THE SECOND HALF
      DO I=2,NXLCYC
      DO J=J_0S,J_NYP
       DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
       DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
       DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
       ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
       ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

       AA1=ETA(I,J+1)+ETA(I+1,J+1)
       AA2=ETA(I,J)+ETA(I+1,J)
       AA5=-(ETA(I,J+1)+ETA(I+1,J+1)-ETA(I,J)-ETA(I+1,J))*TNG(J)
       AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

       AV(I,J)=(-AA2*DELY2+ETAMEAN*DELYR*(TNG(J-1)-2.0*TNG(J)))*UVM(I,J)
       BV(I,J)=((AA1+AA2)*DELY2+AA5*DELYR+AA6*BYRAD2
     & +AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
       CV(I,J)=(-AA1*DELY2-ETAMEAN*DELYR*(TNG(J+1)-2.0*TNG(J)))*UVM(I,J)
      END DO
      END DO

      if (hasSouthPole(grid_ICDYN)) then
        DO I=2,NXLCYC
          AV(I,2)=0.0
c         CV(2,I)=CV(2,I)/BV(2,I)  ! absorbed into TRIDIAG
        END DO
      end if

      if (hasNorthPole(grid_ICDYN)) then
        DO I=2,NXLCYC
          CV(I,NYPOLE)=0.0
        END DO
      end if

      DO I=2,NXLCYC
      DO J=J_0S,J_NYP
        DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
        DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
        DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
        ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))

        AA1=((ETA(I+1,J)  +ZETA(I+1,J)  )*BYCSU(J)+
     *      (ETA(I+1,J+1)+ZETA(I+1,J+1))*BYCSU(J))*BYCSU(J)
        AA2=((ETA(I,J)+ZETA(I,J))*BYCSU(J)+(ETA(I,J+1)+ZETA(I,J+1))
     &      *BYCSU(J))*BYCSU(J)

        IF(J.EQ.NYPOLE) THEN
          AA9=( (ETA(I,J+1)+ETA(I+1,J+1))*DELY2*UICEC(I,J+1)
     &    +ETAMEAN*DELYR*(TNG(J+1)-2.0*TNG(J))*UICEC(I,J+1) )*UVM(I,J)
     &       *FLOAT(NPOL-0)
        ELSE
          AA9=0.0
        END IF

        FXY1a(I,J)=AA9+AMASS(I,J)*BYDTS*UICE(I,J,1)*2.0
     5  -(AA1+AA2)*DELX2*UICE(I,J,1)
     6  +((ETA(I+1,J)+ZETA(I+1,J)+ETA(I+1,J+1)+ZETA(I+1,J+1))
     6  *UICE(I+1,J,1)
     6  +(ETA(I,J)+ZETA(I,J)+ETA(I,J+1)+ZETA(I,J+1))*UICE(I-1,J,1))
     6  *DELX2*BYCSU(J)*BYCSU(J)

      END DO
      END DO

      DO I=2,NXLCYC
      DO J=J_0S,J_NYP
      VRT(I,J)=FXY(I,J)+FXY1a(I,J)
      VRT(I,J)=VRT(I,J)*UVM(I,J)
      END DO
      END DO

      CALL TRIDIAG_new(AV, BV, CV, VRT, U_tmp, grid_ICDYN, 2, NYPOLE)

      DO I=2,NXLCYC
      DO J=J_0S,J_NYP
        UICE(I,J,1)=U_tmp(I,J)     ! VRT(J)   !
      END DO
      END DO

C NOW DO VICE
C THE FIRST HALF

      CALL HALO_UPDATE(GRID_ICDYN,  ETA,FROM=NORTH)
      CALL HALO_UPDATE(GRID_ICDYN, ZETA,FROM=NORTH)

      DO I=2,NXLCYC
      DO J=J_0S,J_NYP
      DELXY=BYDXDY(I,J)    ! 0.5/(DXU(I)*DYU(J))
      DELXR=BYDXR(I)       ! 0.5/(DXU(I)*RADIUS)
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      FXYa(I,J)=-DRAGA(I,J)*UICEC(I,J)+FORCEY(I,J)
     3+(0.5*(UICEC(I+1,J)-UICEC(I-1,J))*(ZETA(I,J+1)+ZETA(I+1,J+1)
     3-ZETA(I,J)-ZETA(I+1,J))*DELXY*BYCSU(J)+0.5*ZETAMEAN*
     3((UICEC(I+1,J+1)
     3-UICEC(I-1,J+1))*BYCSU(J+1)-(UICEC(I+1,J-1)-UICEC(I-1,J-1))
     3*BYCSU(J-1))*DELXY)
     3
     4-(0.5*(UICEC(I+1,J)-UICEC(I-1,J))*(ETA(I,J+1)+ETA(I+1,J+1)
     4-ETA(I,J)-ETA(I+1,J))*DELXY*BYCSU(J)+0.5*ETAMEAN*((UICEC(I+1,J+1)
     4-UICEC(I-1,J+1))*BYCSU(J+1)-(UICEC(I+1,J-1)-UICEC(I-1,J-1))
     4*BYCSU(J-1))*DELXY)
     4
     5+0.5*(ETA(I+1,J+1)*(UICEC(I+1,J+1)+UICEC(I,J+1)
     5-UICEC(I+1,J)-UICEC(I,J))+ETA(I+1,J)*(UICEC(I+1,J)
     5+UICEC(I,J)-UICEC(I+1,J-1)-UICEC(I,J-1))+ETA(I,J+1)
     5*(UICEC(I,J)+UICEC(I-1,J)-UICEC(I,J+1)-UICEC(I-1,J+1))
     5+ETA(I,J)*(UICEC(I,J-1)+UICEC(I-1,J-1)-UICEC(I,J)
     5-UICEC(I-1,J)))*DELXY*BYCSU(J)
     5
     6+(ETA(I+1,J+1)+ETA(I+1,J)-ETA(I,J)-ETA(I,J+1))
     6*TNG(J)*UICEC(I,J)*DELXR*BYCSU(J)
     6+ETAMEAN*TNG(J)*(UICEC(I+1,J)-UICEC(I-1,J))*DELXR*BYCSU(J)
     6
     7+ETAMEAN*2.0*TNG(J)*(UICEC(I+1,J)-UICEC(I-1,J))*DELXR*BYCSU(J)

      END DO
      END DO

      DO I=2,NXLCYC
      DO J=J_0S,J_NYP
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))
      AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
      AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
      AA3=(ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*BYCSU(J)
      AA4=(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*BYCSU(J)
      AA5=((ZETA(I,J+1)-ETA(I,J+1))+(ZETA(I+1,J+1)-ETA(I+1,J+1))
     &-(ZETA(I,J)-ETA(I,J))-(ZETA(I+1,J)-ETA(I+1,J)))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      AV(I,J)=(-AA2*DELY2-(ZETAMEAN-ETAMEAN)*TNG(J-1)*DELYR
     &-ETAMEAN*2.0*TNG(J)*DELYR)*UVM(I,J)
      BV(I,J)=((AA1+AA2)*DELY2+AA5*DELYR+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CV(I,J)=(-AA1*DELY2+(ZETAMEAN-ETAMEAN)*TNG(J+1)*DELYR
     &+ETAMEAN*2.0*TNG(J)*DELYR)*UVM(I,J)
      END DO
      END DO

      if (hasSouthPole(grid_ICDYN)) then
        DO I=2,NXLCYC
          AV(I,2)=0.0
c         CV(2,I)=CV(2,I)/BV(2,I)  ! absorbed into TRIDIAG
        END DO
      end if

      if (hasNorthPole(grid_ICDYN)) then
        DO I=2,NXLCYC
          CV(I,NYPOLE)=0.0
        END DO
      end if

      DO I=2,NXLCYC
      DO J=J_0S,J_NYP
        DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
        DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
        DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)

        ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
        ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

        AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
        AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
        AA3=(ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*BYCSU(J)
        AA4=(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*BYCSU(J)
        AA5=((ZETA(I,J+1)-ETA(I,J+1))+(ZETA(I+1,J+1)-ETA(I+1,J+1))
     &  -(ZETA(I,J)-ETA(I,J))-(ZETA(I+1,J)-ETA(I+1,J)))*TNG(J)
        AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

        IF(J.EQ.NYPOLE) THEN
          AA9=(AA1*DELY2-(ZETAMEAN-ETAMEAN)*TNG(J+1)*DELYR
     &    -ETAMEAN*2.0*TNG(J)*DELYR)*VICEC(I,J+1)*UVM(I,J)*FLOAT(NPOL-0)
        ELSE
          AA9=0.0
        END IF

        VRT(I,J)=AA9+FXYa(I,J)-(AA3+AA4)*DELX2*VICEC(I,J)
     6 +((ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*VICEC(I+1,J)*DELX2
     7    +(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*VICEC(I-1,J)*DELX2)
     *       *BYCSU(J)
        VRT(I,J)=(VRT(I,J)+AMASS(I,J)*BYDTS*VICE(I,J,2)*2.0)*UVM(I,J)
      END DO
      END DO

       CALL TRIDIAG_new(AV, BV, CV, VRT, U_tmp, grid_ICDYN, 2, NYPOLE)

      DO I=2,NXLCYC
      DO J=J_0S,J_NYP
        VICE(I,J,1)=U_tmp(I,J)   ! VRT(J)   !
      END DO
      END DO

      DO I=2,NXLCYC
       DO J=J_0S,J_NYP
        VICE(I,J,3)=VICE(I,J,1)
       END DO
      END DO

C NOW THE SECOND HALF

      DO J=J_0S,J_NYP
       DO I=2,NXLCYC
        DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
        DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
        DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
        ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
        ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

        AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
        AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
        AA3=(ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*BYCSU(J)
        AA4=(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*BYCSU(J)
        AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

        AU(I,J)=-AA4*DELX2*UVM(I,J)
        BU(I,J)=((AA3+AA4)*DELX2+AA6*BYRAD2
     &  +AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
        CU(I,J)=-AA3*DELX2*UVM(I,J)
       END DO
      END DO

      CALL HALO_UPDATE(grid_ICDYN, VICE, FROM=SOUTH+NORTH)

      DO J=J_0S,J_NYP
      DO I=2,NXLCYC

        DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
        DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
        DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
        ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
        ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

        AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
        AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
        AA3=(ETA(I+1,J)*BYCSU(J)+ETA(I+1,J+1)*BYCSU(J))*BYCSU(J)
        AA4=(ETA(I,J)*BYCSU(J)+ETA(I,J+1)*BYCSU(J))*BYCSU(J)
        AA5=((ZETA(I,J+1)-ETA(I,J+1))+(ZETA(I+1,J+1)-ETA(I+1,J+1))
     &  -(ZETA(I,J)-ETA(I,J))-(ZETA(I+1,J)-ETA(I+1,J)))*TNG(J)
        AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

        FXY1(I,J)=AMASS(I,J)*BYDTS*VICE(I,J,1)*2.0
     1  -AA5*DELYR*VICE(I,J,1)
     1  -(AA1+AA2)*DELY2*VICE(I,J,1)
     1  +AA1*DELY2*VICE(I,J+1,1)-((ZETAMEAN-ETAMEAN)*TNG(J+1)*DELYR
     1  +ETAMEAN*2.0*TNG(J)*DELYR)*VICE(I,J+1,1)
     2  +AA2*DELY2*VICE(I,J-1,1)+((ZETAMEAN-ETAMEAN)*TNG(J-1)*DELYR
     2  +ETAMEAN*2.0*TNG(J)*DELYR)*VICE(I,J-1,1)

      END DO
      END DO

      DO 1201 J=J_0S,J_NYP
        DO I=2,NXLCYC
          URT(I)=FXYa(I,J)+FXY1(I,J)
          URT(I)=URT(I)*UVM(I,J)
        END DO

        CALL TRIDIAG_cyclic(AU(2,J),BU(2,J),CU(2,J),URT(2),VICE(2,J,1),
     &               NXLCYC-1)

 1201 CONTINUE

      DO J=J_0S,J_NYP
        DO I=2,NXLCYC
          UICE(I,J,1)=UICE(I,J,1)*UVM(I,J)
          VICE(I,J,1)=VICE(I,J,1)*UVM(I,J)
        END DO
      END DO

      RETURN
      END SUBROUTINE RELAX


      Subroutine GEOMICDYN
      use DOMAIN_DECOMP_1D, only : getDomainBounds,
     &     hasSouthPole, hasNorthPole
      USE CONSTANT, only : OMEGA,RADIUS,PI,TWOPI,radian
      IMPLICIT NONE
      INTEGER :: J_0,J_1,J_0S,J_1S
      INTEGER :: I,J!,K
      REAL*8 :: LAT1,COSP1,DXP1
      REAL*8 :: phit,phiu,fjeq,acor,acoru
      REAL*8 :: SINV,SINVm1
      REAL*8 :: FIM,BYIM

C****
C**** Extract useful local domain parameters from ice dynamics grid
C****
      call getDomainBounds(grid_ICDYN, J_STRT =J_0,    J_STOP =J_1,
     &     J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S)

C****
C**** calculate grid and initialise arrays
C****

c**** Geometry inherited from lat-lon B-grid (geomb)
C**** latitudinal spacing depends on whether you have even spacing or
C**** a partial box at the pole

      FIM=IMICDYN
      BYIM=1./FIM
      DLON = TWOPI*BYIM

      DLON_DG=360./REAL(IMICDYN)
      DLAT_DG=180./REAL(JMICDYN-1)   ! 1/2 box at pole (default)
      IF (JMICDYN.eq.90) DLAT_DG=2.  ! full box at pole for 2 deg res
      IF (JMICDYN.eq.180) DLAT_DG=1. ! full box at pole for 1 deg res
      IF (JMICDYN.eq.24) DLAT_DG=180./REAL(JMICDYN-1.5) ! 1/4 box at pole, 'real' 8x10
      DLATM=60.*DLAT_DG
      DLAT=DLAT_DG*radian

      acor=1.  ; acoru=1.
! full box at pole for 2 and 1 deg res
      if (ny1.eq.90 .or. ny1.eq.180) then   ! HARD-CODED DIMENSION
        acor=1.5d0 ; acoru=2.
      end if
! 1/4 box at pole, 'real' 8x10
      if (ny1.eq.24) then   ! HARD-CODED DIMENSION
        acor=0.75d0 ; acoru=0.5d0
      end if

c****
      do j = 1,ny1
        dyt(j) = dlat*radius
        dyu(j) = dlat*radius
      enddo

C**** polar box corrections
      IF (hasNorthPole(grid_ICDYN)) THEN
        dyt(ny1)=dyt(ny1)*acor
        dyu(ny1)=dyu(ny1)*acoru
      END IF
      IF (hasSouthPole(grid_ICDYN)) THEN
        dyt(1)=dyt(1)*acor
        dyu(1)=dyu(1)*acoru
      END IF

      do i=1,nx1
        dxt(i) = dlon*radius
        dxu(i) = dlon*radius
      enddo
      dxt(1) = dxt(nx1-1)
      dxt(nx1) = dxt(2)
      dxu(1) = dxu(nx1-1)
      dxu(nx1) = dxu(2)

      do i=1,nx1
        bydx2(i)=0.5/(dxu(i)*dxu(i))
        bydxr(i)=0.5/(dxu(i)*radius)
      end do
      do j=1,ny1
        bydy2(j)=0.5/(dyu(j)*dyu(j))
        bydyr(j)=0.5/(dyu(j)*radius)
      end do

      do j=j_0,j_1
        do i=1,nx1
          bydxdy(i,j) = 0.5/(dxu(i)*dyu(j))
        end do
      end do

      fjeq=0.5*(ny1+1)          ! equatorial index
      phit = 90.*radian
      cst(ny1) = cos(phit)
c      tngt(ny1)= sin(phit)/cst(ny1)
      phit = -90.*radian
      cst(1) = cos(phit)
c      tngt(1)= sin(phit)/cst(1)

      do j = 2,ny1-1
       phit = (j-fjeq)*dlat
       cst(j) = cos(phit)
       tngt(j)= sin(phit)/cst(j)
      enddo

      do j = 1,ny1
        phiu = (j-fjeq+0.5)*dlat
        csu(j) = cos(phiu)
        bycsu(j) = 1./csu(j)
        tng(j) = sin(phiu)/csu(j)
      end do

      do j = j_0,j_1
        phiu = (j-fjeq+0.5)*dlat
        do i=1,nx1
          sinen(i,j)=sin(phiu)
        enddo
      enddo

C**** fix some polar fields
C**** (can't use 'have_north_pole', needs adjacent boxes too!)
      TNGT(NY1)=TNGT(NY1-1)
      tng(ny1)=tng(ny1-1)
      csu(ny1)=csu(ny1-1)
      bycsu(ny1) = 1./csu(ny1)

      TNGT(1)=TNGT(2)
      tng(1) =tng(2)
      csu(1) =csu(2)
      bycsu(1) = 1./csu(1)

      LAT(1)  = -.25*TWOPI
      LAT(JMICDYN) = -LAT(1)
      LAT_DG(1,1:2)  = -90
      LAT_DG(JMICDYN,1) = +90
      !SINP(1)  = -1.
      !SINP(JMICDYN) = 1.
      !COSP(1)  = 0.
      !COSP(JMICDYN) = 0.
      DXP(1)  = 0.
      DXP(JMICDYN) = 0.
      DO J=2,JMICDYN-1
        LAT(J)  = DLAT*(J-FJEQ)
        LAT_DG(J,1)  = DLAT_DG*(J-FJEQ)
        !SINP(J) = SIN(LAT(J))
        !COSP(J) = COS(LAT(J))
        DXP(J)  = RADIUS*DLON*COS(LAT(J)) !COSP(J)
      END DO
      !BYDXP(2:JMICDYN-1) = 1.D0/DXP(2:JMICDYN-1)
      LAT1    = DLAT*(1.-FJEQ)
      COSP1   = COS(LAT1)
      DXP1    = RADIUS*DLON*COSP1
      DO J=2,JMICDYN
        LAT_DG(J,2)  = DLAT_DG*(J-FJEQ-.5)
        !COSV(J) = .5*(COSP(J-1)+COSP(J))
        DXV(J)  = .5*(DXP(J-1)+DXP(J))
        DYV(J)  = RADIUS*(LAT(J)-LAT(J-1))
C**** The following corrections have no effect for half polar boxes
C**** but are important for full and quarter polar box cases.
        IF (J.eq.2) THEN
          DXV(J)  = .5*(DXP1+DXP(J))
        END IF
        IF (J.eq.JMICDYN) THEN
          DXV(J)  = .5*(DXP(J-1)+DXP1)
        END IF
C****
      END DO
      DYP(1)  = RADIUS*(LAT(2)-LAT(1)-0.5*DLAT)
      DYP(JMICDYN) = RADIUS*(LAT(JMICDYN)-LAT(JMICDYN-1)-0.5*DLAT)

      SINV    = Sin (DLAT*(1+.5-FJEQ))
      DXYP(1) = RADIUS*RADIUS*DLON*(SINV+1)
      BYDXYP(1) = 1./DXYP(1)

      !aDXYP(:,1) = DXYP(1)

      SINVm1  = Sin (DLAT*(JMICDYN-.5-FJEQ))
      DXYP(JMICDYN)= RADIUS*RADIUS*DLON*(1-SINVm1)
      BYDXYP(JMICDYN) = 1./DXYP(JMICDYN)

      !aDXYP(:,JMICDYN) = DXYP(JMICDYN)

      DXYS(1)  = 0.
      DXYS(JMICDYN) = DXYP(JMICDYN)
      DXYN(1)  = DXYP(1)
      DXYN(JMICDYN) = 0.

      DO J=2,JMICDYN-1
        DYP(J)  =  radius*dlat
        SINVm1  = Sin (DLAT*(J-.5-FJEQ))
        SINV    = Sin (DLAT*(J+.5-FJEQ))
        DXYP(J) = RADIUS*RADIUS*DLON*(SINV-SINVm1)

        !aDXYP(:,J) = DXYP(J)

        BYDXYP(J) = 1./DXYP(J)
        DXYS(J) = .5*DXYP(J)
        DXYN(J) = .5*DXYP(J)
      END DO
      DO J=2,JMICDYN
        DXYV(J) = DXYN(J-1)+DXYS(J)
        !BYDXYV(J) = 1./DXYV(J)
      END DO
C**** imaxj
      !IMAXJ(1)=1
      !IMAXJ(JMICDYN)=1
      !IMAXJ(2:JMICDYN-1)=IMICDYN

C**** LONGITUDES (degrees) for diagnostics
      LON_DG(1,1) = -180.+0.5*dlon_dg
      LON_DG(1,2) = -180.+dlon_dg
      DO I=2,IMICDYN
        LON_DG(I,1) = LON_DG(I-1,1)+dlon_dg
        LON_DG(I,2) = LON_DG(I-1,2)+dlon_dg
      END DO
      lon(:) = lon_dg(:,1)*radian

c store b-grid lons/lats for interpolations
      lonb(:) = lon_dg(:,2)*radian + pi
      latb(1:jmicdyn-1) = lat_dg(2:jmicdyn,2)*radian
      latb(jmicdyn) = pi/2. ! not used

      DO I=1,IMICDYN
        SINIU(I) =SIN (I*TWOPI/IMICDYN)
        COSIU(I) =COS (I*TWOPI/IMICDYN)
      END DO
      SINIU(IMICDYN) = 0
      COSIU(IMICDYN) = 1

      end subroutine GEOMICDYN

      Subroutine ICDYN_MASKS
      use DOMAIN_DECOMP_1D, only : getDomainBounds
      use DOMAIN_DECOMP_1D, only : NORTH,SOUTH,HALO_UPDATE
      use pario, only : par_open,par_close,read_dist_data
      IMPLICIT NONE
      INTEGER :: J_0,J_1,J_0S,J_1S
      INTEGER :: I,J,K,fid
      real*8, allocatable :: msk(:,:)
      logical :: ex
C****
C**** Extract useful local domain parameters from ice dynamics grid
C****
      call getDomainBounds(grid_ICDYN, J_STRT =J_0,    J_STOP =J_1,
     &     J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S )

      allocate(
     &     msk(imicdyn,grid_icdyn%j_strt_halo:grid_icdyn%j_stop_halo))
      inquire(file='ICEDYN_MASKFAC',exist=ex)
      if(ex) then
        ! If this file exists, we multiply the velocity mask uvm by
        ! its contents.  This is a way of manually restricting flow in
        ! coastal regions where heffm=ceiling(focean) allows flow
        ! through thin strips of land.
        fid = par_open(grid_icdyn,'ICEDYN_MASKFAC','read')
        call read_dist_data(grid_icdyn,fid,'mask',msk)
        call par_close(grid_icdyn,fid)
      endif

      icedyn_processors_only: if(grid_ICDYN%have_domain) then

C**** sin/cos ice-ocean turning angle
      SINWAT=SIN(OIPHI)
      COSWAT=COS(OIPHI)

C**** Set land masks for tracer and velocity points
       do j=j_0,j_1
        do i=2,nx1-1

#if (defined CUBED_SPHERE) || (defined HYCOM1degRefine) ||\
    (defined HYCOM1degUnrefined) || (defined HYCOM2deg)
          heffm(i,j)=nint(focean(i-1,j))
#else
          heffm(i,j)=ceiling(focean(i-1,j))
#endif
        enddo
        heffm(1,j)=heffm(nx1-1,j)
        heffm(nx1,j)=heffm(2,j)
      enddo
C**** define velocity points (including exterior corners)
      CALL HALO_UPDATE(grid_ICDYN, HEFFM, FROM=NORTH)
      do j=j_0,j_1s
        do i=1,nx1-1
c          sumk=heffm(i,j)+heffm(i+1,j)+heffm(i,j+1)+heffm(i+1,j+1)
c          if (sumk.ge.3) uvm(i,j)=1  ! includes exterior corners
          uvm(i,j) = nint(min(heffm(i,j), heffm(i+1,j), heffm(i,j+1),
     *         heffm(i+1,j+1)))
        end do
      end do
C**** reset tracer points to surround velocity points (except for single
      CALL HALO_UPDATE(grid_ICDYN, UVM, FROM=SOUTH)
      do j=j_0s,j_1s
        do i=2,nx1-1
          k = nint(max (uvm(i,j), uvm(i-1,j), uvm(i,j-1), uvm(i-1,j-1)))
c         sumk = nint(uvm(i,j)+uvm(i+1,j)+uvm(i,j+1)+uvm(i+1,j+1))
c set to k except if an island
c         if (.not. (sumk.eq.4.and.focean(i-1,j).eq.0) ) then
            heffm(i,j) = k
c         end if
        enddo
      enddo
C**** final sweep to reinstate islands
c      do j=j_0s,j_1s
c        do i=2,nx1-1
c          sumk = nint(uvm(i,j)+uvm(i+1,j)+uvm(i,j+1)+uvm(i+1,j+1))
c          if (sumk.eq.4.and.heffm(i,j).eq.0.) then
c            uvm(i,j)=0 ; uvm(i+1,j)=0 ; uvm(i,j+1)=0 ; uvm(i+1,j+1)=0
c          end if
c        enddo
c      enddo
c set lateral boundary conditions
       do j=j_0,j_1
        heffm(1,j)   = heffm(nx1-1,j)
        heffm(nx1,j) = heffm(2,j)
      enddo

C**** Update halo of PHI for distributed memory implementation
      CALL HALO_UPDATE(grid_ICDYN, HEFFM, FROM=NORTH)
      do j=j_0,j_1s
        do i=1,nx1-1
          uvm(i,j) = nint(min(heffm(i,j), heffm(i+1,j), heffm(i,j+1),
     *         heffm(i+1,j+1)))
        end do
      end do

      if(ex) then
        do j=j_0,j_1s
          ! note: array from disk indexed w/o extra halo cells at IDL
          uvm(2:imicdyn+1,j) = uvm(2:imicdyn+1,j)*msk(1:imicdyn,j)
        enddo
      endif

c set cyclic conditions on eastern and western boundary
       do j=j_0,j_1S
        uvm(1,j) = uvm(nx1-1,j)
        uvm(nx1,j) = uvm(2,j)
      enddo

      endif icedyn_processors_only

      deallocate(msk)

      end subroutine ICDYN_MASKS

      END MODULE ICEDYN



      SUBROUTINE VPICEDYN
!@sum  vpicedyn is the entry point into the viscous-plastic ice
!@+    dynamics code
!@auth Gavin Schmidt (based on code from J. Zhang)
      USE DOMAIN_DECOMP_1D, only : getDomainBounds,NORTH,SOUTH,GLOBALSUM
      USE DOMAIN_DECOMP_1D, ONLY : HALO_UPDATE,am_i_root
      USE MODEL_COM, only : itime,qcheck
      USE ICEDYN, only : form,relax,uice,vice,uicec,vicec
     *        ,uvm,dxu,dyu,amass
      USE ICEDYN, only : NX1, NY1
      USE ICEDYN, only : grid=>grid_ICDYN
      IMPLICIT NONE
      REAL*8, DIMENSION(NX1,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &        USAVE,VSAVE
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &        rms_part,area_part
      REAL*8 rms0,rms,area
      INTEGER kki,i,j
      INTEGER :: J_0,J_1,J_0S,J_1S
      logical :: debug

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S )

      rms=0. ; rms0=0. ; debug=.false.
      area=0.
C KKI LOOP IS FOR PSEUDO-TIMESTEPPING
      KKI=0.
 10   KKI=KKI+1

C FIRST DO PREDICTOR
      DO J=J_0,J_1
      DO I=1,NX1
       UICE(I,J,3)=UICE(I,J,1)
       VICE(I,J,3)=VICE(I,J,1)
       UICEC(I,J)=UICE(I,J,1)
       VICEC(I,J)=VICE(I,J,1)
      END DO
      END DO

      CALL FORM
      CALL RELAX

       DO J=J_0,J_1
       UICE(1,J,1)=UICE(NX1-1,J,1)
       VICE(1,J,1)=VICE(NX1-1,J,1)
       UICE(NX1,J,1)=UICE(2,J,1)
       VICE(NX1,J,1)=VICE(2,J,1)
      END DO

C NOW DO REGULAR TIME STEP
C NOW DO MODIFIED EULER STEP
c
      DO J=J_0,J_1
      DO I=1,NX1
       UICE(I,J,1)=0.5*(UICE(I,J,1)+UICE(I,J,2))
       VICE(I,J,1)=0.5*(VICE(I,J,1)+VICE(I,J,2))
      END DO
      END DO

      CALL FORM

C NOW SET U(1)=U(2) AND SAME FOR V
      DO J=J_0,J_1
      DO I=1,NX1
       UICE(I,J,3)=UICE(I,J,1)
       VICE(I,J,3)=VICE(I,J,1)
       UICEC(I,J)=UICE(I,J,1)
       VICEC(I,J)=VICE(I,J,1)
       UICE(I,J,1)=UICE(I,J,2)
       VICE(I,J,1)=VICE(I,J,2)
      END DO
      END DO

      CALL RELAX

      DO J=J_0,J_1
       UICE(1,J,1)=UICE(NX1-1,J,1)
       VICE(1,J,1)=VICE(NX1-1,J,1)
       UICE(NX1,J,1)=UICE(2,J,1)
       VICE(NX1,J,1)=VICE(2,J,1)
      END DO

      if (kki.gt.1) then ! test convergence
        rms_part=0. ; area_part=0.
        do i=1,nx1
           do j=j_0,j_1
            if(amass(i,j)*UVM(i,j)>0) then
              rms_part(j)=rms_part(j)+  dxu(i)*dyu(j)*
     *         ((USAVE(i,j)-UICE(i,j,1))**2+(VSAVE(i,j)-VICE(i,j,1))**2)
              area_part(j)=area_part(j)+dxu(i)*dyu(j)
            end if
          end do
        end do
        CALL GLOBALSUM(grid, rms_part, rms, all=.true.)
        CALL GLOBALSUM(grid, area_part, area, all=.true.)
        if(debug.and.am_i_root()) write(6,*) itime,kki,4d4*rms/area
      end if

      if (kki > 2 .and. rms > rms0 .and. qcheck) then
        debug=.true.
        if(am_i_root())
     *    write(6,*) 'VPICEDYN rms rose, kki:',itime,kki,
     *      4d4*rms0/area,4d4*rms/area
      end if
      if (kki.eq.20) then
        rms=sqrt(rms/area)
        if(am_i_root())
     *    write(6,*) "20 iterations in VPICEDYN: itime,rms",itime,rms
      elseif (kki == 1 .or. rms > 25d-6*area) then
        USAVE=UICE(:,:,1)
        VSAVE=VICE(:,:,1)
        rms0=rms
        goto 10
      end if

      RETURN
      END SUBROUTINE VPICEDYN

      SUBROUTINE ALLOC_ICEDYN(IM_in,JM_in)
!@sum ALLOC_ICEDYN allocates arrays defined in the ICEDYN module.
!@auth Rosalinda de Fainchtein

C**** arrays allocated in this routine were originally dimensioned
C**** (..,ny1,..). Since ny1=jm (see above), the grid structure as defined
C**** in DOMAIN_DECOMP can be used in the calling routine.
C**** In the case that ny1 is NOT equal to JM, a structure appropriately
C**** modified to reflect the differences should be created in DOMAIN_DECOMP
C**** and used in the calling routine. No modification should be necesary
C**** to ALLOC_ICEDYN.
      !USE RESOLUTION, only : im,jm
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID,getDomainBounds
      USE DOMAIN_DECOMP_1D, ONLY : INIT_GRID
      USE ICEDYN, ONLY : FOCEAN
      USE ICEDYN, ONLY : PRESS,HEFFM,UVM,DWATN,COR,ZMAX,ZMIN,ETA,
     &                   ZETA,DRAGS,DRAGA,GAIRX,GAIRY,GWATX,GWATY,
     &                   PGFUB,PGFVB,FORCEX,FORCEY,AMASS,UICEC,
     &                   VICEC,DMU,DMV,HEFF,AREA,UICE,
     &                   VICE,SINEN,BYDXDY,DYT,DYU,BYDY2,BYDYR,
     &                   CST,CSU,TNGT,TNG,BYCSU,USI,VSI,
     &                   DXT,DXU,BYDX2,BYDXR
      USE ICEDYN, ONLY :
     &     LAT,LATB,LON,LONB,DXYP,BYDXYP,DXYV,DXP,DYP,DXV,DYV,DXYS,DXYN,
     &     LON_DG,LAT_DG,SINIU,COSIU
      USE ICEDYN, only : IMICDYN,JMICDYN,NX1,NY1,NYPOLE,NXLCYC
      USE ICEDYN, ONLY :  grid_NXY, grid_ICDYN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IM_in,JM_in
      LOGICAL, SAVE :: init = .false.

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

#ifdef CUBED_SPHERE
      IMICDYN = 2*IM_in
      JMICDYN = 2*IM_in
#else
      IMICDYN = IM_in
      JMICDYN = JM_in
#endif

      NX1=IMICDYN+2
      NY1=JMICDYN
      NYPOLE=NY1-1
      NXLCYC=NX1-1

      ALLOCATE(DXT(NX1),DXU(NX1),BYDX2(NX1),BYDXR(NX1))
      ALLOCATE(DYT(NY1),DYU(NY1),BYDY2(NY1),BYDYR(NY1))
      ALLOCATE(CST(NY1),CSU(NY1),TNGT(NY1),TNG(NY1),BYCSU(NY1))

      ALLOCATE(LAT(JMICDYN),LATB(JMICDYN))
      ALLOCATE(LON(IMICDYN),LONB(IMICDYN))
      ALLOCATE(DXYP(JMICDYN),BYDXYP(JMICDYN))
      ALLOCATE(DXYV(JMICDYN))
      ALLOCATE(DXP(JMICDYN),DYP(JMICDYN))
      ALLOCATE(DXV(JMICDYN),DYV(JMICDYN))
      ALLOCATE(DXYS(JMICDYN),DXYN(JMICDYN))
      ALLOCATE(LAT_DG(JMICDYN,2))
      ALLOCATE(LON_DG(IMICDYN,2))
      ALLOCATE(SINIU(IMICDYN),COSIU(IMICDYN))

c***   - grid_NXY is the ice dynamics grid (for the resolution of the momentum equation)
c***     it is a latlon grid with dimensions IMICDYN+2 and JMICDYN,
c***     both defined in the ICEDYN module
c***   - grid_ICDYN is the same grid as grid_NXY, with dimensions IMICDYN and JMICDYN,
c***     the two boundary ghost cells in the longitudinal direction have been removed

      CALL INIT_GRID(grid_NXY,NX1,NY1,1,npes_max=JMICDYN/3)
      CALL INIT_GRID(grid_ICDYN,IMICDYN,JMICDYN,1,npes_max=JMICDYN/3)

      call getDomainBounds(grid_NXY, I_STRT_HALO=I_0H, I_STOP_HALO=I_1H,
     &                J_STRT_HALO=J_0H, J_STOP_HALO=J_1H  )

      ALLOCATE( FOCEAN(NX1-2,J_0H:J_1H),
     $   STAT = IER)

      ALLOCATE(  PRESS(NX1,J_0H:J_1H),
     &           HEFFM(NX1,J_0H:J_1H),
     &           UVM(NX1,J_0H:J_1H),
     &           DWATN(NX1,J_0H:J_1H),
     &           COR(NX1,J_0H:J_1H),
     *           ZMAX(NX1,J_0H:J_1H),
     &           ZMIN(NX1,J_0H:J_1H),
     &           ETA(NX1,J_0H:J_1H),
     &           ZETA(NX1,J_0H:J_1H),
     &           DRAGS(NX1,J_0H:J_1H),
     &           DRAGA(NX1,J_0H:J_1H),
     &           GAIRX(NX1,J_0H:J_1H),
     &           GAIRY(NX1,J_0H:J_1H),
     *           GWATX(NX1,J_0H:J_1H),
     &           GWATY(NX1,J_0H:J_1H),
     &           PGFUB(NX1,J_0H:J_1H),
     &           PGFVB(NX1,J_0H:J_1H),
     &           FORCEX(NX1,J_0H:J_1H),
     &           FORCEY(NX1,J_0H:J_1H),
     &           AMASS(NX1,J_0H:J_1H),
     &           UICEC(NX1,J_0H:J_1H),
     &           VICEC(NX1,J_0H:J_1H),
     &           DMU(NX1,J_0H:J_1H),
     &           DMV(NX1,J_0H:J_1H),
     &           HEFF(NX1,J_0H:J_1H),
     &           AREA(NX1,J_0H:J_1H),
     &           USI(1:imicdyn,J_0H:J_1H),
     &           VSI(1:imicdyn,J_0H:J_1H),
     $   STAT = IER)
      ALLOCATE( UICE(NX1,J_0H:J_1H,3),
     &          VICE(NX1,J_0H:J_1H,3),
     $   STAT = IER)
C**** Geometry
      ALLOCATE( SINEN(NX1,J_0H:J_1H),
     &          BYDXDY(NX1,J_0H:J_1H),
     $   STAT = IER)
      RETURN
      END SUBROUTINE ALLOC_ICEDYN
