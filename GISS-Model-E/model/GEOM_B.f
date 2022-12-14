#include "rundeck_opts.h"
      MODULE GEOM
!@sum  GEOM contains spherical geometric variables and arrays
!@auth Original development team
!@cont GEOM_B
      USE CONSTANT, only : OMEGA,RADIUS,TWOPI,radian,AREAG
      USE RESOLUTION, only : IM,JM,LM
      IMPLICIT NONE
      PRIVATE

      public :: GEOM_ATM, areag
      public :: lonlat_to_ij
      public :: lonlat_to_tile
      public :: lon_to_I
      public :: lat_to_J

C**** The primary grid is the A grid (including both poles)
C**** The secondary grid is for the B grid velocities, located on the
C**** vertices of the A grid (Note: first velocity point is J=2)
C**** Polar boxes can have different latitudinal size and are treated
C**** as though they were 1/IM of their actual area
      SAVE

!@param IMH half the number of longitudinal boxes
      INTEGER, PARAMETER, public :: IMH=IM/2
!@param FIM,BYIM real values related to number of long. grid boxes
      REAL*8, PARAMETER, public :: FIM=IM, BYIM=1./FIM

!@param  DLON grid spacing in longitude (deg)
      REAL*8, PUBLIC, PARAMETER :: DLON = TWOPI*BYIM
C**** For the wonderland model set DLON=DLON/3
c      REAL*8, PUBLIC, PARAMETER :: DLON=TWOPI/(IM*3)
!@var DLAT,DLAT_DG,DLATM grid spacing in latitude (rad,deg,minutes)
      REAL*8, PUBLIC  :: DLAT,DLAT_DG, DLATM
!@param FJEQ equatorial value of J
      REAL*8, PUBLIC, PARAMETER :: FJEQ=.5*(1+JM)
!@var  J1U index of southernmost latitude (currently 2, later 1)
      INTEGER, parameter :: J1U = 2
!@var  JRANGE_HEMI lowest,highest lat index for SH,NH for A,B grid
      INTEGER, PUBLIC, parameter, dimension(2,2,2) :: JRANGE_HEMI = 
     *     reshape(
     *  (/1,JM/2,  1+JM/2,JM,  J1U,J1U-1+JM/2, J1U-1+JM/2,JM+J1U-2/),
     *  (/2,2,2/))
!@var  LAT latitude of mid point of primary grid box (radians)
!@var  LATV latitude of southern edge of primary grid box (radians)
      REAL*8, PUBLIC, DIMENSION(JM) :: LAT,LATV
!@var  LAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, PUBLIC, DIMENSION(JM,2) :: LAT_DG
!@var  LON longitude of mid points of primary grid box (radians)
!@var  LONV longitude of east edge of primary grid box (radians)
      REAL*8, PUBLIC, DIMENSION(IM) :: LON,LONV
!@var  LON_DG longitude of mid points of prim. and sec. grid boxes (deg)
      REAL*8, PUBLIC, DIMENSION(IM,2) :: LON_DG
!@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
C**** Note that this is not the exact area, but is what is required for
C**** some B-grid conservation quantities
      REAL*8, PUBLIC, DIMENSION(JM) :: DXYP,BYDXYP
      REAL*8, PUBLIC, DIMENSION(IM,JM) :: aDXYP

!@var  i-index (longitude) of each grid cell
      INTEGER, PUBLIC, ALLOCATABLE :: indx(:,:)
!@var  j-index (latitude) of each grid cell
      INTEGER, PUBLIC, ALLOCATABLE :: jndx(:,:)

      REAL*8, PUBLIC, DIMENSION(:,:), ALLOCATABLE ::
     &     AXYP,BYAXYP,LAT2D,LON2D,LAT2D_DG,LON2D_DG,SINLAT2D,COSLAT2D
     &    ,ddx_ci,ddx_cj,ddy_ci,ddy_cj
!@var WTJ area weighting used in JLMAP, JKMAP (for hemispheric means)
      REAL*8, PUBLIC, DIMENSION(JM,2,2) :: WTJ

!@var DXYV,BYDXYV area of grid box around velocity point (recip.)(m^2)
      REAL*8, PUBLIC, DIMENSION(JM) :: DXYV,BYDXYV

!@var  DXP,DYP,BYDXP,BYDYP distance between points on primary grid
!@+     (+inverse)
      REAL*8, PUBLIC, DIMENSION(JM) :: DXP,DYP,BYDXP,BYDYP
!@var  DXV,DYV distance between velocity points (secondary grid)
      REAL*8, PUBLIC, DIMENSION(JM) :: DXV,DYV
!@var  DXYN,DXYS half box areas to the North,South of primary grid point
      REAL*8, PUBLIC, DIMENSION(JM) :: DXYS,DXYN
!@var  SINP sin of latitude at primary grid points
!@var  SINLATV,COSLATV sin(latv), cos(latv)
!@var  COSP, COSV cos of latitude at primary, secondary latitudes
      REAL*8, PUBLIC, DIMENSION(JM) :: SINP,SINLATV
      REAL*8, PUBLIC, DIMENSION(JM) :: COSP,COSV,COSLATV,DXLATV
!@var  RAPVS,RAPVN,RAVPS,RAVPN area scalings for primary and sec. grid
      REAL*8, PUBLIC, DIMENSION(JM) :: RAPVS,RAPVN,RAVPS,RAVPN

!@var SINIV,COSIV,SINIP,COSIP longitud. sin,cos for wind,pressure grid
      REAL*8, PUBLIC, DIMENSION(IM) :: SINIV,COSIV,SINIP,COSIP,SINU,COSU
!@var  RAVJ scaling for A grid U/V to B grid points (func. of lat. j)
!@var  RAPJ scaling for B grid -> A grid conversion (1/4,1/im at poles)
      REAL*8, PUBLIC, DIMENSION(IM,JM) :: RAPJ,RAVJ
!@var  IDJJ J index of adjacent U/V points for A grid (func. of lat. j)
      INTEGER, PUBLIC, DIMENSION(IM,JM) :: IDJJ
!@var  IDIJ I index of adjacent U/V points for A grid (func. of lat/lon)
      INTEGER, PUBLIC, DIMENSION(:,:,:), allocatable :: IDIJ
!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, PUBLIC, DIMENSION(JM) :: KMAXJ
!@var  IMAXJ varying number of used longitudes
      INTEGER, PUBLIC, DIMENSION(JM) :: IMAXJ
!@var  FCOR latitudinally varying coriolis parameter
      REAL*8, PUBLIC, DIMENSION(JM) :: FCOR

!@var JG_U, JG_KE lat. grids on which U-wind and KE are defined
!@+   (1 for primary latitudes, 2 for secondary latitudes)
!@+   Information for diagnostics.
      integer, parameter :: jg_u=2, jg_ke=2

      real*8, public :: acor,acor2,polwt

      CONTAINS

      SUBROUTINE GEOM_ATM
!@sum  GEOM_ATM Calculate spherical geometry for B grid
!@auth Original development team (modifications by G. Schmidt)
      use domain_decomp_atm, only : grid, hasSouthPole, hasNorthPole,
     *                              Am_I_Root
      IMPLICIT NONE

      INTEGER :: I,J,K,IM1  !@var I,J,K,IM1  loop variables
      INTEGER :: JVPO,JMHALF
      REAL*8  :: RAVPO,LAT1,COSP1,DXP1, SINV,SINVm1

      integer :: i_0h,i_1h,j_0h,j_1h,i_0,i_1,j_0,j_1,j_0s,j_1s

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo
      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      j_0s = grid%j_strt_skp
      j_1s = grid%j_stop_skp

      allocate(idij(im,im,grid%j_strt_halo:grid%j_stop_halo))

      allocate(indx(i_0h:i_1h, j_0h:j_1h))
      allocate(jndx(i_0h:i_1h, j_0h:j_1h))

      If (Am_I_Root()) Write (6,900)
      Write (6,901) GRID%RANK,
     *              i_0h,i_0,i_1,i_1h, j_0h,j_0,j_0s,j_1s,j_1,j_1h,
     *              HasSouthPole(Grid), HasNorthPole(Grid)
  900 Format ('GEOM_B:   Rank  I_0H   I_0   I_1  I_1H',
     *               '   J_0H   J_0  J_0S  J_1S   J_1  J_1H   QSP  QNP') 
  901 Format ('GEOM_B: ',5I6,1X,6I6,2L5)

      allocate(
     &       axyp(i_0h:i_1h,j_0h:j_1h)
     &    ,byaxyp(i_0h:i_1h,j_0h:j_1h)
     &    ,lat2d(i_0h:i_1h,j_0h:j_1h)
     &    ,lat2d_dg(i_0h:i_1h,j_0h:j_1h)
     &    ,lon2d(i_0h:i_1h,j_0h:j_1h)
     &    ,lon2d_dg(i_0h:i_1h,j_0h:j_1h)
     &    ,sinlat2d(i_0h:i_1h,j_0h:j_1h)
     &    ,coslat2d(i_0h:i_1h,j_0h:j_1h)
     &    ,ddx_ci(i_0h:i_1h,j_0h:j_1h)
     &    ,ddx_cj(i_0h:i_1h,j_0h:j_1h)
     &    ,ddy_ci(i_0h:i_1h,j_0h:j_1h)
     &    ,ddy_cj(i_0h:i_1h,j_0h:j_1h)
     &     )

C**** latitudinal spacing depends on whether you have even spacing or
C**** a partial box at the pole

      DLAT_DG=180./REAL(JM)                   ! even spacing (default)
      IF (JM.eq.46) DLAT_DG=180./REAL(JM-1)   ! 1/2 box at pole for 4x5
cc    IF (JM.eq.24) DLAT_DG=180./REAL(JM-1)   ! 1/2 box at pole, orig 8x10
      IF (JM.eq.24) DLAT_DG=180./REAL(JM-1.5) ! 1/4 box at pole, 'real' 8x10
      DLATM=60.*DLAT_DG
      DLAT=DLAT_DG*radian
      LAT(1)  = -.25*TWOPI
      LAT(JM) = -LAT(1)
      SINP(1)  = -1.
      SINP(JM) = 1.
      COSP(1)  = 0.
      COSP(JM) = 0.
      DXP(1)  = 0.
      DXP(JM) = 0.
      DXV = 0.
      COSV = 0.
      DO J=2,JM-1
        LAT(J)  = DLAT*(J-FJEQ)
        SINP(J) = SIN(LAT(J))
        COSP(J) = COS(LAT(J))
        DXP(J)  = RADIUS*DLON*COSP(J)
      END DO
      BYDXP(2:JM-1) = 1.D0/DXP(2:JM-1)
      LAT1    = DLAT*(1.-FJEQ)
      COSP1   = COS(LAT1)
      DXP1    = RADIUS*DLON*COSP1
      DO J=2,JM
        COSV(J) = .5*(COSP(J-1)+COSP(J))
        DXV(J)  = .5*(DXP(J-1)+DXP(J))
        DYV(J)  = RADIUS*(LAT(J)-LAT(J-1))
C**** The following corrections have no effect for half polar boxes
C**** but are important for full and quarter polar box cases.
        IF (J.eq.2) THEN
          polwt = cosv(j)
          COSV(J) = .5*(COSP1+COSP(J))
          DXV(J)  = .5*(DXP1+DXP(J))
        END IF
        IF (J.eq.JM) THEN
          COSV(J) = .5*(COSP(J-1)+COSP1)
          DXV(J)  = .5*(DXP(J-1)+DXP1)
        END IF
C****
      END DO
      DYP(1)  = RADIUS*(LAT(2)-LAT(1)-0.5*DLAT)
      DYP(JM) = RADIUS*(LAT(JM)-LAT(JM-1)-0.5*DLAT)

      SINV    = Sin (DLAT*(1+.5-FJEQ))
      DXYP(1) = RADIUS*RADIUS*DLON*(SINV+1)
      BYDXYP(1) = 1./DXYP(1)

      aDXYP(:,1) = DXYP(1) 

      SINVm1  = Sin (DLAT*(JM-.5-FJEQ))
      DXYP(JM)= RADIUS*RADIUS*DLON*(1-SINVm1)
      BYDXYP(JM) = 1./DXYP(JM)

      aDXYP(:,JM) = DXYP(JM) 

      DXYS(1)  = 0.
      DXYS(JM) = DXYP(JM)
      DXYN(1)  = DXYP(1)
      DXYN(JM) = 0.
      polwt = (cosv(3)-cosv(2))/(cosv(3)-polwt)
      DO J=2,JM-1
        DYP(J)  =  radius*dlat !.5*(DYV(J)+DYV(J+1))
        SINVm1  = Sin (DLAT*(J-.5-FJEQ))
        SINV    = Sin (DLAT*(J+.5-FJEQ))
        DXYP(J) = RADIUS*RADIUS*DLON*(SINV-SINVm1)

        aDXYP(:,J) = DXYP(J) 

        BYDXYP(J) = 1./DXYP(J)
        DXYS(J) = .5*DXYP(J)
        DXYN(J) = .5*DXYP(J)
      END DO
      do j=2,jm
        latv(j) = dlat*(j-.5-fjeq)
        coslatv(j) = cos(latv(j))
        sinlatv(j) = sin(latv(j))
        dxlatv(j) = radius*dlon*coslatv(j)
      enddo
      BYDYP(:) = 1.D0/DYP(:)
      RAVPS(1)  = 0.
      RAPVS(1)  = 0.
      RAVPN(JM) = 0.
      RAPVN(JM) = 0.
      DO J=2,JM
        DXYV(J) = DXYN(J-1)+DXYS(J)
        BYDXYV(J) = 1./DXYV(J)
        RAPVS(J)   = .5*DXYS(J)/DXYV(J)
        RAPVN(J-1) = .5*DXYN(J-1)/DXYV(J)
        RAVPS(J)   = .5*DXYS(J)/DXYP(J)
        RAVPN(J-1) = .5*DXYN(J-1)/DXYP(J-1)
      END DO
      acor = dxyv(2)/(.5*dxp(2)*dyv(2)) ! gridbox area correction factor
      acor2 = dxyv(2)/(dxv(2)*dyv(2))
C**** LONGITUDES (degrees); used in ILMAP
      LON_DG(1,1) = -180.+360./(2.*FLOAT(IM))
      LON_DG(1,2) = -180.+360./    FLOAT(IM)
      DO I=2,IM
        LON_DG(I,1) = LON_DG(I-1,1)+360./FLOAT(IM)
        LON_DG(I,2) = LON_DG(I-1,2)+360./FLOAT(IM)
      END DO
C**** LATITUDES (degrees); used extensively in the diagn. print routines
      LAT_DG(1,1:2)=-90.
      LAT_DG(JM,1)=90.
      DO J=2,JM-1
        LAT_DG(J,1)=DLAT_DG*(J-FJEQ)    ! primary (tracer) latitudes
      END DO
      DO J=2,JM
        LAT_DG(J,2)=DLAT_DG*(J-JM/2-1)  ! secondary (velocity) latitudes
      END DO
C**** WTJ: area weighting for JKMAP, JLMAP hemispheres
      JMHALF= JM/2
      DO J=1,JM
        WTJ(J,1,1)=1.
        WTJ(J,2,1)=2.*FIM*DXYP(J)/AREAG
      END DO
      DO J=2,JM
        WTJ(J,1,2)=1.
        WTJ(J,2,2)=2.*FIM*DXYV(J)/AREAG
      END DO
cgsfc      WTJ(JMHALF+1,1,2)=.5
cgsfc      WTJ(JMHALF+1,2,2)=WTJ(JMHALF+1,2,2)/2.
      WTJ(1,1,2)=0.
      WTJ(1,2,2)=0.
C**** CALCULATE CORIOLIS PARAMETER
c      OMEGA = TWOPI*(EDPERD+EDPERY)/(EDPERD*EDPERY*SDAY)
      FCOR(1)  = -OMEGA*DXV(2)*DXV(2)/DLON
      FCOR(JM) = OMEGA*DXV(JM)*DXV(JM)/DLON
      DO J=2,JM-1
        FCOR(J) = OMEGA*(DXV(J)*DXV(J)-DXV(J+1)*DXV(J+1))/DLON
      END DO

C**** Set indexes and scalings for the influence of A grid points on
C**** adjacent velocity points

C**** Calculate relative directions of polar box to nearby U,V points
      DO I=1,IM
        SINIV(I)=SIN((I-1)*DLON)
        COSIV(I)=COS((I-1)*TWOPI*BYIM) ! DLON)
        LON(I)=DLON*(I-.5)
        LONV(I)=DLON*I
        SINIP(I)=SIN(LON(I))
        COSIP(I)=COS(LON(I))
        SINU(I) =SIN (I*TWOPI/IM)
        COSU(I) =COS (I*TWOPI/IM)
      END DO
      SINU(IM) = 0
      COSU(IM) = 1

C**** Conditions at the poles
      DO J=1,JM,JM-1
        IF(J.EQ.1) THEN
          JVPO=2
          RAVPO=2.*RAPVN(1)
        ELSE
          JVPO=JM
          RAVPO=2.*RAPVS(JM)
        END IF
        KMAXJ(J)=IM
        IMAXJ(J)=1
        RAVJ(1:KMAXJ(J),J)=RAVPO
        RAPJ(1:KMAXJ(J),J)=BYIM
        IDJJ(1:KMAXJ(J),J)=JVPO
        if(j.ge.grid%j_strt_halo .and. j.le.grid%j_stop_halo) then
          DO K=1,KMAXJ(J)
            IDIJ(K,1:IM,J)=K
          END DO
        endif
      END DO
C**** Conditions at non-polar points
      DO J=2,JM-1
        KMAXJ(J)=4
        IMAXJ(J)=IM
        DO K=1,2
          RAVJ(K,J)=RAPVS(J)
          RAPJ(K,J)=RAVPS(J)    ! = .25
          IDJJ(K,J)=J
          RAVJ(K+2,J)=RAPVN(J)
          RAPJ(K+2,J)=RAVPN(J)  ! = .25
#ifdef SCM
          IDJJ(K+2,J)=J
#else
          IDJJ(K+2,J)=J+1
#endif
        END DO
        if(j.ge.grid%j_strt_halo .and. j.le.grid%j_stop_halo) then
          IM1=IM
          DO I=1,IM
#ifdef SCM
            IDIJ(1,I,J)=I
            IDIJ(2,I,J)=I
            IDIJ(3,I,J)=I
            IDIJ(4,I,J)=I
#else 
            IDIJ(1,I,J)=IM1
            IDIJ(2,I,J)=I
            IDIJ(3,I,J)=IM1
            IDIJ(4,I,J)=I
            IM1=I
#endif
          END DO
        endif
      END DO

c
c temporary
c
      ! Fill in indices for grid cells, so we can
      ! halo-update them with our neighbors
      do j=j_0h,j_1h
      do i=i_0h,i_1h
        indx(i,j) = i
        jndx(i,j) = j
      end do    ! i
      end do    ! j


      do j=max(1,j_0h),min(jm,j_1h)
      do i=i_0h,i_1h
        axyp(i,j) = dxyp(j)
        byaxyp(i,j) = bydxyp(j)
        lat2d(i,j) = lat(j)
        lat2d_dg(i,j) = lat_dg(j,1)
        sinlat2d(i,j) = sin(lat(j))
        coslat2d(i,j) = cos(lat(j))
      enddo
      do i=i_0,i_1
        lon2d(i,j) = lon(i) ! i=1 is dlon/2
        lon2d_dg(i,j) = lon_dg(i,1) ! i=1 is -180 + dlon/2
      enddo
      enddo

      do j=j_0s,j_1s
      do i=i_0,i_1
        ddx_ci(i,j) =  .5/(radius*dlon*cosp(j))
        ddx_cj(i,j) = 0.
        ddy_ci(i,j) = 0.
        ddy_cj(i,j) =  .5/(radius*dlat)
c        dloni = lon2d(i+1,j)-lon2d(i-1,j)
c        dlati = lat2d(i+1,j)-lat2d(i-1,j)
c        dlonj = lon2d(i,j+1)-lon2d(i,j-1)
c        dlatj = lat2d(i,j+1)-lat2d(i,j-1)
c        if(dloni.lt.0.) dloni=dloni+twopi
c        if(dlonj.lt.0.) dlonj=dlonj+twopi
c        bydet = 1d0/(dloni*dlatj-dlonj*dlati)
c        ddx_ci(i,j) =  dlatj*bydet/(radius*coslat2d(i,j))
c        ddx_cj(i,j) = -dlonj*bydet/(radius*coslat2d(i,j))
c        ddy_ci(i,j) = -dlati*bydet/radius
c        ddy_cj(i,j) =  dloni*bydet/radius
      enddo
      enddo

C**** set up mapping arrays for budget/conserv diags
      call set_j_budg   !called after lat2d_dg is initialized
      call set_wtbudg() !sets area weights
 
      RETURN
      END SUBROUTINE GEOM_ATM

      subroutine lonlat_to_ij(ll,ij)
c converts lon,lat=ll(1:2) into model i,j=ij(1:2) (primary grid)
c this version is for the latlon grid.  ll are in degrees east.
      implicit none
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: ij(2)
      real*8 :: dlon_dg
      dlon_dg = 360./dble(im)

      ij(1) = nint( .5*(im + 1) + ll(1)/dlon_dg )
      ij(2) = nint( .5*(jm + 1) + ll(2)/dlat_dg )
      if(ij(2)>jm) ij(2)=jm
      if(ij(2)<1 ) ij(2)=1
      return
      end subroutine lonlat_to_ij

      subroutine lonlat_to_tile(ll,tile)
c returns tile = 1, latlon grid has only one face 
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: tile
      tile=1
      return
      end subroutine lonlat_to_tile


      function lon_to_I(deg)
c converts longitude into model i (primary grid)
c this version is for the latlon grid.  DEG is in degrees east.
      implicit none
      integer lon_to_I
      real*8, intent(in) :: DEG
      real*8 :: dlon_dg
      dlon_dg = 360./dble(im)
      lon_to_I = nint( .5*(im + 1) + DEG/dlon_dg )
      return
      end function lon_to_I

      function lat_to_J(deg)
c converts latitude into model j (primary grid)
c this version is for the latlon grid.
      implicit none
      integer lat_to_J
      real*8, intent(in) :: DEG
      lat_to_J = nint( .5*(jm + 1) + DEG/dlat_dg )
      if( lat_to_J>jm ) lat_to_J = jm
      if( lat_to_J<1 )  lat_to_J = 1
      return
      end function lat_to_J

      END MODULE GEOM

      module RAD_COSZ0
      use resolution, only : im,jm
      implicit none
      save
      private
      public :: cosz_init,coszs,coszt,daily_cosz

!@var SIND,COSD sin,cos of solar declination angle
!@+   (these are local copies that are set by daily_cosz)
      real*8 :: sind,cosd
!@var SINJ,COSJ sines and cosines for zenith angle calculation
      REAL*8, ALLOCATABLE, DIMENSION(:) :: SINJ,COSJ

      contains

      subroutine daily_cosz(SIND_in,COSD_in, COSZ_day,DUSK)
c Resets parameters needed for zenith angle calculations, and calculates
c daily-average cosz and the hour of dusk
      USE DOMAIN_DECOMP_ATM, ONLY: grid
      USE MODEL_COM
      USE CONSTANT, only : PI
      implicit none
      REAL*8 SIND_in,COSD_in
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     COSZ_day, ! average of cos(zenith angle) for the current day
     &     DUSK      ! hour of dusk, radians from local noon
      INTEGER :: J
      REAL*8 :: SJSD,CJCD,CJCD_SDUSK
      SIND = SIND_in
      COSD = COSD_in
      do j=grid%j_strt,grid%j_stop
        SJSD=SINJ(J)*SIND
        CJCD=COSJ(J)*COSD
        IF(CJCD-SJSD.LT.0.) THEN     ! constant daylight
          DUSK(:,J) = PI
          COSZ_day(:,J) = SJSD
        ELSEIF(CJCD+SJSD.LT.0.) THEN ! constant darkness
          DUSK(:,J) = 0D0
          COSZ_day(:,J) = 0D0
        ELSE
          DUSK(:,J) = ACOS(-SJSD/CJCD)
          CJCD_SDUSK = SQRT((CJCD-SJSD)*(CJCD+SJSD))
          COSZ_day(:,J) = max(0d0,(SJSD*DUSK(1,J)+CJCD_SDUSK)/PI)
        ENDIF
      enddo
      return
      end subroutine daily_cosz

      subroutine cosz_init(cosz_const)
      use domain_decomp_atm, only : grid, hasSouthPole, hasNorthPole
      use constant, only : twopi
      use resolution, only : jm
      use geom, only : dlat
      real*8, optional :: cosz_const ! for interface compatibility
      real*8 :: phis,cphis,sphis,phin,cphin,sphin,phim
      integer :: j_0,j_1,j_0s,j_1s,j_0h,j_1h,j

      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      j_0s = grid%j_strt_skp
      j_1s = grid%j_stop_skp
    
      allocate(sinj(j_0h:j_1h),cosj(j_0h:j_1h))

C**** COMPUTE THE AREA WEIGHTED LATITUDES AND THEIR SINES AND COSINES
      if (hasSouthPole(GRID)) then
        PHIS=-.25*TWOPI
        SPHIS=-1.
        CPHIS=0.
      else
        PHIS=DLAT*(J_0-1-.5*JM)
        SPHIS=SIN(PHIS)
        CPHIS=COS(PHIS)
      end if
      DO J=J_0,J_1S
        PHIN=DLAT*(J-.5*JM)
        SPHIN=SIN(PHIN)
        CPHIN=COS(PHIN)
        PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
        SINJ(J)=SIN(PHIM)
        COSJ(J)=COS(PHIM)
        PHIS=PHIN
        SPHIS=SPHIN
        CPHIS=CPHIN
      END DO
      IF (hasNorthPole(GRID)) THEN
        PHIN=.25*TWOPI
        SPHIN=1.
        CPHIN=0.
        PHIM=( PHIN*SPHIN + CPHIN
     *        -PHIS*SPHIS - CPHIS)
     *        /(SPHIN - SPHIS)
        SINJ(JM)=SIN(PHIM)
        COSJ(JM)=COS(PHIM)
      END IF
      return
      end subroutine cosz_init

      subroutine COSZS (ROT1,ROT2,COSZ,COSZA)
      REAL*8 ROT1,ROT2
      REAL*8, DIMENSION(:,:) :: COSZ,COSZA
      call COSZT (ROT1,ROT2,COSZ,COSZA)
      end subroutine COSZS

      subroutine COSZT (ROT1,ROT2,COSZ,COSZA)
!@sum  COSZ0 calculates Earth's zenith angle, weighted by time/sunlight
!@auth Original Development Team
      USE CONSTANT, only : twopi
      USE MODEL_COM
      USE GEOM, only : lon,sinip,cosip
c      USE RAD_COM, only : cosd,sind
      USE DOMAIN_DECOMP_ATM, ONLY: grid, hasSouthPole, hasNorthPole
      IMPLICIT NONE
      SAVE
      REAL*8 ROT1,ROT2
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     COSZ
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
     &     , optional :: COSZA
      REAL*8, DIMENSION(IM) :: LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
C**** ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
C**** COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      REAL*8, PARAMETER :: ZERO1=1.D-2
      INTEGER I,J
      REAL*8 S2DAWN,S2DUSK,ECOSZ,ECOSQZ,CLT1,CLT2,ZERO2,CDUSK,DUSK,DAWN
     *     ,SDUSK,SDAWN,CJCD,SJSD,SR1,CR1,SR2,CR2,DROT
      INTEGER :: J_0, J_1, J_0S, J_1S

C****
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP

      if ( present(COSZA) ) goto 777 ! COSZS
C****
!      ENTRY COSZT (ROT1,ROT2,COSZ)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE WEIGHTED BY DAYTIME
C**** HOURS FROM ROT1 TO ROT2, GREENWICH MEAN TIME IN RADIANS.  ROT1
C**** MUST BE BETWEEN 0 AND 2*PI.  ROT2 MUST BE BETWEEN ROT1 AND
C**** ROT1+2*PI.  I=1 MUST LIE ON THE INTERNATIONAL DATE LINE.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
      SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO I=1,IM
        LT1(I)=ROT1+LON(I)
        SLT1(I)=SR1*COSIP(I)+CR1*SINIP(I)
        LT2(I)=ROT2+LON(I)
        SLT2(I)=SR2*COSIP(I)+CR2*SINIP(I)
      END DO
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO J=1,JM,JM-1
        IF(((J .EQ. 1) .AND. (hasSouthPole(grid))) .OR.
     *     ((J .EQ. JM) .AND. (hasNorthPole(grid)))) THEN
           SJSD=SINJ(J)*SIND
           CJCD=COSJ(J)*COSD
           IF (SJSD+CJCD.GT.ZERO1) THEN
              IF (SJSD-CJCD.LT.0.) THEN
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
                 DUSK=ACOS(-SJSD/CJCD)
                 SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
                 DAWN=-DUSK
                 SDAWN=-SDUSK
                 COSZ(1,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/TWOPI
              ELSE
C**** CONSTANT DAYLIGHT NEAR THE POLES
                 COSZ(1,J)=SJSD
              END IF
           ELSE
C**** CONSTANT NIGHTIME NEAR THE POLES
              COSZ(1,J)=0.
           END IF
        END IF
      END DO
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 500 J=J_0S,J_1S
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.GT.ZERO1) THEN
      IF (SJSD-CJCD.LT.0.) THEN
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      DUSK=ACOS(-SJSD/CJCD)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      DAWN=-DUSK
      SDAWN=-SDUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 400 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 220
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  220 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 240
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      GO TO 400
  240 IF (DAWN.GE.LT1(I)) GO TO 300
      IF (DUSK.LT.LT2(I)) GO TO 260
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
      GO TO 400
  260 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 280
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I)))/DROT
      GO TO 400
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  280 COSZ(I,J)=(SJSD*(LT2(I)-DAWN-TWOPI+DUSK-LT1(I))+CJCD*
     *  (SLT2(I)-SDAWN+SDUSK-SLT1(I)))/DROT
      GO TO 400
  300 IF (DUSK.LT.LT2(I)) GO TO 320
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN))/DROT
      GO TO 400
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
  320 COSZ(I,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/DROT
  400 CONTINUE
      ELSE
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
        DO I=1,IM
          COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
        END DO
      END IF
      ELSE
C**** CONSTANT NIGHTIME AT THIS LATITUDE
        COSZ(1:IM,J)=0.
      END IF
  500 CONTINUE
      RETURN
C****
C****
!      ENTRY COSZS (ROT1,ROT2,COSZ,COSZA)
 777  continue
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE TWICE, FIRST WEIGHTED BY THE
C**** DAYTIME HOURS FROM ROT1 TO ROT2 AND SECONDLY WEIGHTED BY THE
C**** INCIDENT SUN LIGHT FROM ROT1 TO ROT2.  COSZT MUST HAVE BEEN
C**** CALLED JUST PREVIOUSLY.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
      SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO I=1,IM
        LT1(I)=ROT1+LON(I)
        SLT1(I)=SR1*COSIP(I)+CR1*SINIP(I)
        CLT1=CR1*COSIP(I)-SR1*SINIP(I)
        S2LT1(I)=2.*SLT1(I)*CLT1
        LT2(I)=ROT2+LON(I)
        SLT2(I)=SR2*COSIP(I)+CR2*SINIP(I)
        CLT2=CR2*COSIP(I)-SR2*SINIP(I)
        S2LT2(I)=2.*SLT2(I)*CLT2
      END DO
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      COSZ  = 0.0
      COSZA = 0.0
      DO J=1,JM,JM-1
        IF(((J .EQ. 1) .AND. (hasSouthPole(grid))) .OR.
     *     ((J .EQ. JM) .AND. (hasNorthPole(grid)))) THEN
           SJSD=SINJ(J)*SIND
           CJCD=COSJ(J)*COSD
           IF (SJSD+CJCD.GT.ZERO1) THEN
              IF (SJSD-CJCD.LT.0.) THEN
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
                 CDUSK=-SJSD/CJCD
                 DUSK=ACOS(CDUSK)
                 SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
                 S2DUSK=2.*SDUSK*CDUSK
                 DAWN=-DUSK
                 SDAWN=-SDUSK
                 S2DAWN=-S2DUSK
                 ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
                 ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *                .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
                 COSZ(1,J)=ECOSZ/TWOPI
                 COSZA(1,J)=ECOSQZ/ECOSZ
              ELSE
C**** CONSTANT DAYLIGHT NEAR THE POLES
                 ECOSZ=SJSD*TWOPI
                 ECOSQZ=SJSD*ECOSZ+.5*CJCD*CJCD*TWOPI
                 COSZ(1,J)=ECOSZ/TWOPI
                 COSZA(1,J)=ECOSQZ/ECOSZ
              END IF
           ELSE
C**** CONSTANT NIGHTIME NEAR THE POLES
              COSZ(1,J)=0.
              COSZA(1,J)=0.
           END IF
        END IF
      END DO
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 900 J=J_0S,J_1S
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.GT.ZERO1) THEN
      IF (SJSD-CJCD.LT.0.) THEN
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      CDUSK=-SJSD/CJCD
      DUSK=ACOS(CDUSK)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      S2DUSK=2.*SDUSK*CDUSK
      DAWN=-DUSK
      SDAWN=-SDUSK
      S2DAWN=-S2DUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 800 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 620
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  620 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 640
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      COSZA(I,J)=0.
      GO TO 800
  640 IF (DAWN.GE.LT1(I)) GO TO 700
      IF (DUSK.LT.LT2(I)) GO TO 660
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *  .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  660 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 680
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      ECOSZ=SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I))+
     *  .5*CJCD*(DUSK-LT1(I)+.5*(S2DUSK-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  680 ECOSZ=SJSD*(DROT-DAWN-TWOPI+DUSK)+
     *  CJCD*(SLT2(I)-SDAWN+SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I)+SLT2(I)-SDAWN)+
     *  .5*CJCD*(DUSK+DROT-DAWN-TWOPI+
     *  .5*(S2DUSK-S2LT1(I)+S2LT2(I)-S2DAWN)))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  700 IF (DUSK.GE.LT2(I)) THEN
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
        ECOSZ=SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SDAWN)+
     *       .5*CJCD*(LT2(I)-DAWN+.5*(S2LT2(I)-S2DAWN)))
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      ELSE
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
        ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *       .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      END IF
  800 CONTINUE
      ELSE
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
        DO I=1,IM
          ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
          ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *         .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
          COSZ(I,J)=ECOSZ/DROT
          COSZA(I,J)=ECOSQZ/ECOSZ
        END DO
      END IF
C**** CONSTANT NIGHTIME AT THIS LATITUDE
      ELSE
        COSZ(1:IM,J)=0.
        COSZA(1:IM,J)=0.
      END IF
  900 CONTINUE
      RETURN
      END subroutine COSZT

      end module RAD_COSZ0
