#include "rundeck_opts.h"

c
c Definitions of non-dimensional cube coordinates x,y:
c
c Except where indicated otherwise, all routines in this file assume
c that input/output x,y are defined to vary linearly with grid indices
c i,j (im=jm):
c x_edge(i) = -1+2*(i-1)/im    x_midpoint(i) = -1+2*(i-1/2)/im
c y_edge(j) = -1+2*(j-1)/jm    y_midpoint(j) = -1+2*(j-1/2)/jm
c
c Grid spacing is taken to give equal lengths along cube edges, i.e.
c y = lat/g when x=1; g is the latitude of cube corners.
c g = acos(1/3)/2 = asin(1/sqrt(3)) = atan(1/sqrt(2))
c
c For convenience, many routines internally define an alternative
c space denoted here as "capital" X,Y:
c X = sqrt(2)*tan(g*x)  Y = sqrt(2)*tan(g*y)
c X = tan(rotated_lon)  Y = tan(rotated_lat)/cos(rotated_lon)
c In addition to allowing many formulas to be written more compactly,
c this transformation has the following properties:
c - great circles are straight lines in X-Y space
c - latitude circles on polar tiles are X*X + Y*Y = const
c - latitude circles on equatorial tiles are Y*Y = const*(1+X*X)
c

      MODULE GEOM
!@sum  GEOM contains geometric variables and arrays
!@auth M. Kelley
!@cont GEOM_ATM

      USE RESOLUTION, only : IM,JM
      USE CONSTANT, only : radius,twopi,areag
      IMPLICIT NONE
      SAVE

!@var  lat2d latitude of mid point of primary grid box (radians)
      REAL*8, ALLOCATABLE :: lat2d(:,:)
!@var  lat2d_dg latitude of mid point of primary grid box (degrees)
      REAL*8, ALLOCATABLE :: lat2d_dg(:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sinlat2d, coslat2d
!@var  lon2d longitude of mid point of primary grid box (radians)
      REAL*8, ALLOCATABLE :: lon2d(:,:)
!@var  lon2d_dg longitude of mid point of primary grid box (degrees)
      REAL*8, ALLOCATABLE :: lon2d_dg(:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sinlon2d, coslon2d
!@var  lat2d_corner latitude  of lower left corner point (radians)
!@var  lon2d_corner longitude of lower left corner point (radians)
      REAL*8, ALLOCATABLE :: lat2d_corner(:,:)
      REAL*8, ALLOCATABLE :: lon2d_corner(:,:)

!@var lonbds,latbds (degrees)
!@+   lonbds,      latbds      (   1 |     2 |       3 |     4 ,i,j) =
!@+   lon2d_corner,lat2d_corner( i,j | i+1,j | i+1,j+1 | i,j+1 )/radian
      REAL*8, ALLOCATABLE :: latbds(:,:,:),lonbds(:,:,:)

!@var ddx_ci, ddx_cj,ddy_ci ddy_cj coeffs for obtaining N-S and E-W
!@+   gradients from centered differences in the I and J directions
      REAL*8, ALLOCATABLE :: ddx_ci(:,:),ddx_cj(:,:)
      REAL*8, ALLOCATABLE :: ddy_ci(:,:),ddy_cj(:,:)

!@var dlxsina, dlysina edge lengths multiplied by sin(alpha).
!@+   alpha is the local angle between gridlines of constant x and y.
!@+   dlxsina(i,j) is at the boundary between cells i,j-1 and i,j
!@+   dlysina(i,j) is at the boundary between cells i-1,j and i,j
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: dlxsina,dlysina
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: lonuc,latuc,lonvc,latvc

!@var ull2ucs,vll2ucs, ull2vcs,vll2vcs coeffs for projecting
!@+   a latlon-oriented vector with components ull,vll onto local
!@+   cubed-sphere basis vectors e1,e2 that are parallel to
!@+   gridlines of constant y,x:
!@+   ucs = ull*ull2ucs + vll*vll2ucs
!@+   vcs = ull*ull2vcs + vll*vll2vcs
!@+   These instances are defined at C-grid locations.
!@+   ull2ucs,vll2ucs (i,j) lie between cells i-1,j and i,j
!@+   ull2vcs,vll2vcs (i,j) lie between cells i,j-1 and i,j
      REAL*8, DIMENSION(:,:), ALLOCATABLE ::
     &     ull2ucs, vll2ucs, ull2vcs, vll2vcs

c     shift grid 10 degrees West to avoid corner over Japan
      real*8, parameter :: shiftwest = twopi/36.

!@var  axyp,byaxyp area of grid box (+inverse) (m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: axyp, byaxyp

      integer, allocatable, dimension(:) :: imaxj

!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, DIMENSION(JM) :: KMAXJ

      CONTAINS

      SUBROUTINE GEOM_ATM
!@sum  GEOM_ATM Calculate geometry for CS grid
!@auth M. Kelley
      USE CONSTANT, only : RADIUS,PI,TWOPI,radian
      use DOMAIN_DECOMP_ATM, only: grid, getDomainBounds, halo_update
      implicit none
      real*8 :: x,y,x1,x2,y1,y2,e1(2),e2(2)
      integer :: i0h, i1h, i0, i1
      integer :: j0h, j1h, j0, j1
      integer :: i,j,imin,imax,jmin,jmax
      real*8 :: dloni,dlonj,dlati,dlatj,bydet,distx,disty

      real*8, dimension(grid%i_strt:grid%i_stop+1,
     &                  grid%j_strt:grid%j_stop+1) :: axyp_int
      call getDomainBounds(grid, J_STRT_HALO=j0h, J_STOP_HALO=j1h,
     &     J_STRT=j0, J_STOP=j1,
     &     I_STRT_HALO=i0h, I_STOP_HALO=i1h, 
     &     I_STRT=i0,I_STOP=i1)

c      write(*,*) "geom  i0h,i1h,i1,j0h,j1h,j0,j1=",
c     &     i0h,i1h,i1,j0h,j1h,j0,j1

      allocate(lat2d(i0h:i1h, j0h:j1h))
      allocate(lon2d(i0h:i1h, j0h:j1h))

      allocate(lat2d_dg(i0h:i1h, j0h:j1h))
      allocate(lon2d_dg(i0h:i1h, j0h:j1h))

      allocate(sinlat2d(i0h:i1h, j0h:j1h))
      allocate(coslat2d(i0h:i1h, j0h:j1h))

      allocate(sinlon2d(i0h:i1h, j0h:j1h))
      allocate(coslon2d(i0h:i1h, j0h:j1h))

      allocate(lat2d_corner(i0h:i1h+1, j0h:j1h+1))
      allocate(lon2d_corner(i0h:i1h+1, j0h:j1h+1))

      allocate(latbds(4, i0h:i1h, j0h:j1h))
      allocate(lonbds(4, i0h:i1h, j0h:j1h))

      allocate(
     &     ddx_ci(i0h:i1h,j0h:j1h)
     &    ,ddx_cj(i0h:i1h,j0h:j1h)
     &    ,ddy_ci(i0h:i1h,j0h:j1h)
     &    ,ddy_cj(i0h:i1h,j0h:j1h))

      allocate(axyp(i0h:i1h, j0h:j1h))
      allocate(byaxyp(i0h:i1h, j0h:j1h))

      allocate(imaxj(j0:j1))

c
c calculate corner lons/lats, and areas integrated from the center
c of this cube face
c
      do j=j0,j1+1
      do i=i0,i1+1
        x = -1d0 + 2d0*(i-1)/im
        y = -1d0 + 2d0*(j-1)/im
        call csxy2ll(x,y,grid%tile,lon2d_corner(i,j),lat2d_corner(i,j))
        lon2d_corner(i,j)=lon2d_corner(i,j)-shiftwest
        if ( lon2d_corner(i,j) .lt. -pi) lon2d_corner(i,j)=
     &       lon2d_corner(i,j) + twopi
        axyp_int(i,j) = aint(x,y)
      enddo
      enddo
c
c cell areas, lons/lats at cell centers
c
      do j=j0,j1
      do i=i0,i1
        axyp(i,j) = axyp_int(i+1,j+1)-axyp_int(i,j+1)
     &             -axyp_int(i+1,j  )+axyp_int(i,j  )
        axyp(i,j) = axyp(i,j)*radius*radius
        byaxyp(i,j) = 1d0/axyp(i,j)
        x = -1d0 + 2d0*(dble(i)-.5d0)/im
        y = -1d0 + 2d0*(dble(j)-.5d0)/im
        call csxy2ll(x,y,grid%tile,lon2d(i,j),lat2d(i,j))
        lon2d(i,j) = lon2d(i,j)-shiftwest
        lat2d_dg(i,j) = lat2d(i,j)/radian
        lon2d_dg(i,j) = lon2d(i,j)/radian
        if(lon2d_dg(i,j) .lt. -180.) lon2d_dg(i,j)=lon2d_dg(i,j)+360.
        sinlat2d(i,j) = sin(lat2d(i,j))
        coslat2d(i,j) = cos(lat2d(i,j))
        lon2d(i,j) = lon2d(i,j) + pi ! IDL has a value of zero
        if(lon2d(i,j) .lt. 0.) lon2d(i,j)= lon2d(i,j) + twopi
        lonbds(1:2,i,j) = lon2d_corner(i:i+1:+1,j)/radian
        lonbds(3:4,i,j) = lon2d_corner(i+1:i:-1,j+1)/radian
        latbds(1:2,i,j) = lat2d_corner(i:i+1:+1,j)/radian
        latbds(3:4,i,j) = lat2d_corner(i+1:i:-1,j+1)/radian
      enddo
      enddo

c halo everything just in case
      call halo_update(grid,axyp)
      call halo_update(grid,byaxyp)
      call halo_update(grid,lon2d)
      call halo_update(grid,lat2d)
      call halo_update(grid,lon2d_dg)
      call halo_update(grid,lat2d_dg)
      call halo_update(grid,sinlat2d)
      call halo_update(grid,coslat2d)

      call set_j_budg   !called after lat2d_dg is initialized
      call set_wtbudg !sets area weights

      imaxj(:)=i1

      do j=j0,j1
        y1 = -1d0 + 2d0*(j-1-.5)/im
        y2 = -1d0 + 2d0*(j+1-.5)/im
        y  = -1d0 + 2d0*(j  -.5)/im
        do i=i0,i1
          if(j>1 .and. j<jm .and. i>1 .and. i<im) then
            x1 = -1d0 + 2d0*(i-1-.5)/im
            x2 = -1d0 + 2d0*(i+1-.5)/im
            x  = -1d0 + 2d0*(i  -.5)/im
            call e1e2(x,y,grid%tile,e1,e2)
            bydet = 1d0/(e1(1)*e2(2)-e1(2)*e2(1))
            distx = gcdist(x1,x2,y,y)*radius
            disty = gcdist(x,x,y1,y2)*radius
            ddx_ci(i,j) =  e2(2)*bydet/distx
            ddx_cj(i,j) = -e1(2)*bydet/disty
            ddy_ci(i,j) = -e2(1)*bydet/distx
            ddy_cj(i,j) =  e1(1)*bydet/disty
          else
            dloni = lon2d(i+1,j)-lon2d(i-1,j)
            dlati = lat2d(i+1,j)-lat2d(i-1,j)
            dlonj = lon2d(i,j+1)-lon2d(i,j-1)
            dlatj = lat2d(i,j+1)-lat2d(i,j-1)
c account for discontinuity of lon2d at the international date line
            if(abs(dloni).gt.pi) dloni=dloni-twopi*sign(1d0,dloni)
            if(abs(dlonj).gt.pi) dlonj=dlonj-twopi*sign(1d0,dlonj)
            bydet = 1d0/(dloni*dlatj-dlonj*dlati)
            ddx_ci(i,j) =  dlatj*bydet/(radius*coslat2d(i,j))
            ddx_cj(i,j) = -dlati*bydet/(radius*coslat2d(i,j))
            ddy_ci(i,j) = -dlonj*bydet/radius
            ddy_cj(i,j) =  dloni*bydet/radius
          endif
        enddo
      enddo

      kmaxj(:)= 2

c
c Factors for vector transformations and transport calculations.
c For now, ignore discontinuities and changes of orientation
c at the edges of cube faces since these factors are currently used
c only for sea ice transport on polar faces.  Will adapt as needed.
c
c at C-grid u locations
      jmin = max(1,j0h)
      jmax = min(jm,j1h)
      imin = max(1,i0h)
      imax = min(im+1,i1h+1)
      allocate(
     &     dlysina(imin:imax,jmin:jmax)
     &    ,lonuc(imin:imax,jmin:jmax)
     &    ,latuc(imin:imax,jmin:jmax)
     &    ,ull2ucs(imin:imax,jmin:jmax)
     &    ,vll2ucs(imin:imax,jmin:jmax))
      do j=jmin,jmax
        y1 = -1d0 + 2d0*(j-1)/im
        y2 = -1d0 + 2d0*(j  )/im
        y = .5*(y1+y2)
        do i=imin,imax
          x = -1d0 + 2d0*(i-1)/im
          call e1e2(x,y,grid%tile,e1,e2)
          bydet = 1d0/(e1(1)*e2(2)-e1(2)*e2(1))
          ull2ucs(i,j) = +e2(2)*bydet
          vll2ucs(i,j) = -e2(1)*bydet
          dlysina(i,j) = radius*
     &         gcdist(x,x,y1,y2)*sqrt(1d0-sum(e1*e2)**2)
          call csxy2ll(x,y,grid%tile,lonuc(i,j),latuc(i,j))
          lonuc(i,j) = lonuc(i,j)-shiftwest
          lonuc(i,j) = lonuc(i,j) + pi ! IDL has a value of zero
          if(lonuc(i,j) .lt. 0.) lonuc(i,j)= lonuc(i,j) + twopi
        enddo
      enddo
c at C-grid v locations
      jmin = max(1,j0h)
      jmax = min(jm+1,j1h+1)
      imin = max(1,i0h)
      imax = min(im,i1h)
      allocate(
     &     dlxsina(imin:imax,jmin:jmax)
     &    ,lonvc(imin:imax,jmin:jmax)
     &    ,latvc(imin:imax,jmin:jmax)
     &    ,ull2vcs(imin:imax,jmin:jmax)
     &    ,vll2vcs(imin:imax,jmin:jmax))
      do j=jmin,jmax
        y = -1d0 + 2d0*(j-1)/im
        do i=imin,imax
          x1 = -1d0 + 2d0*(i-1)/im
          x2 = -1d0 + 2d0*(i  )/im
          x = .5*(x1+x2)
          call e1e2(x,y,grid%tile,e1,e2)
          bydet = 1d0/(e1(1)*e2(2)-e1(2)*e2(1))
          ull2vcs(i,j) = -e1(2)*bydet
          vll2vcs(i,j) = +e1(1)*bydet
          dlxsina(i,j) = radius*
     &         gcdist(x1,x2,y,y)*sqrt(1d0-sum(e1*e2)**2)
          call csxy2ll(x,y,grid%tile,lonvc(i,j),latvc(i,j))
          lonvc(i,j) = lonvc(i,j)-shiftwest
          lonvc(i,j) = lonvc(i,j) + pi ! IDL has a value of zero
          if(lonvc(i,j) .lt. 0.) lonvc(i,j)= lonvc(i,j) + twopi
        enddo
      enddo

      return
      END SUBROUTINE GEOM_ATM

      subroutine csxy2ll(x,y,tile,lon,lat)
c converts x,y,tile to lon,lat (radians)
c This routine places the center of tile 1 at the IDL.
      USE CONSTANT, only : PI
      implicit none
      real*8 :: x,y ! input
      integer :: tile ! input
      real*8 :: lon,lat ! output
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: gx,gy,tmpx,tmpy,cosgx,tangx,tangy,coslon
      gx = g*x
      gy = g*y
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        tmpx = gx
        tmpy = gy
        gx = +tmpy
        gy = -tmpx
      elseif(tile.eq.6) then ! tile 6 = tile 3 flipped around the axis x+y=0
        tmpx = gx
        tmpy = gy
        gx = -tmpy
        gy = -tmpx
      endif
      if(tile.eq.3 .or. tile.eq.6) then
        tangx = tan(gx)
        tangy = tan(gy)
        lat = atan(1d0/sqrt(2d0*(tangx**2 +tangy**2)+1d-40))
        lon = atan2(tangy,tangx)
        if(tile.eq.6) lat = -lat
      else
        cosgx = cos(gx)
        coslon = cosgx/sqrt(2d0-cosgx*cosgx)
        lat = atan(coslon*sqrt(2d0)*tan(gy))
        lon = sign(acos(coslon),gx)
! add longitude offset to tiles 1,2,4,5. integer arithmetic.
        lon = lon + (mod(tile,3)-1)*pi/2. -pi*(1-(tile/3))
        if(lon.lt.-pi) lon=lon+2.*pi
      endif
      return
      end subroutine csxy2ll

      function aint(x,y)
c calculates the area integral from the center of a cube face
c to the point x,y.
      implicit none
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: x,y
      real*8 :: aint
      real*8 :: tangx,tangy
      tangx = tan(g*x)
      tangy = tan(g*y)
      aint = atan(2.*tangx*tangy/sqrt(1.+2.*(tangx*tangx+tangy*tangy)))
      return
      end function aint

      subroutine lonlat_to_ij(ll,ij)
c Converts lon,lat=ll(1:2) into model i,j=ij(1:2).
c This version is for the gnomonic grid.  ll are in degrees.
c If lon,lat lie on an adjacent tile, i,j correspond to the
c continuation of the local i,j index space to that tile.  If
c lon,lat lie on the opposite side of the cube, ij is set to -99
      use resolution, only : im,jm
      use constant, only : radian,pi,twopi
      use domain_decomp_atm, only : grid
      implicit none
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: ij(2)
      real*8 :: x,y,lonshift
      integer :: tile,n,absn,iorj
      lonshift = ll(1)*radian+shiftwest
      if(lonshift.gt.pi) lonshift=lonshift-twopi
      call ll2csxy(lonshift,ll(2)*radian,x,y,tile)
      ij(1) = min(1+int(.5*(x+1d0)*im),im)
      ij(2) = min(1+int(.5*(y+1d0)*jm),jm)
      absn = abs(tile-grid%tile)
      if(absn.eq.3) then  ! opposite side of the cube
        ij = -99
      elseif(absn.gt.0) then ! rotate/shift to ij space of the local tile
        n = (tile-grid%tile)*(1-2*mod(grid%tile,2))
        iorj = 1                         ! tiles 1-2/3-4/5-6: shift i
        if    ((n-2)*(n+4).eq.0) then    ! rotate clockwise, then shift i
          ij = (/ ij(2), 1+im-ij(1) /)
        elseif((n+2)*(n-4).eq.0) then    ! rotate counterclockwise,
          ij = (/ 1+im-ij(2), ij(1) /)   ! then shift j
          iorj = 2
        elseif((n-1)*(absn-5).eq.0) then ! tiles 1-6/2-3/4-5: shift j
          iorj = 2
        endif
        ij(iorj) = ij(iorj) + im*(1-2*(absn/4))*(tile-grid%tile)/absn
      endif
      return
      end subroutine lonlat_to_ij

      subroutine lonlat_to_tile(ll,tile)
c returns id of cube face containing (lon,lat) point where (lon,lat)=(ll(1),ll(2)) 
      use resolution, only : im,jm
      use constant, only : radian,pi,twopi
      use domain_decomp_atm, only : grid
      implicit none
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: tile 
      real*8 :: x,y,lonshift 
      lonshift = ll(1)*radian+shiftwest
      if(lonshift.gt.pi) lonshift=lonshift-twopi
      call ll2csxy(lonshift,ll(2)*radian,x,y,tile)
      return
      end subroutine lonlat_to_tile

      subroutine ll2csxy(lon,lat,x,y,tile)
c converts lon,lat (radians) to x,y,tile.
c This routine places the center of face 1 at the IDL.
      USE CONSTANT, only : PI
      implicit none
      real*8, parameter :: byg=1.62474893308877d0 ! byg=2/acos(1/3)
      real*8, parameter :: latpol=pi/2d0-1d-10
      real*8 :: lon,lat ! input
      real*8 :: x,y ! output
      integer :: tile ! output
      real*8 :: modlon,coslon,tanlat,tmpx,tmpy,bytanfac
      real*8, parameter :: ctrlon(5)=(/-1.,-.5,0.,.5,1./)*pi
      if(abs(lat).gt.latpol) then ! avoid calculating tan(pi/2)
        x = 0.
        y = 0.
        if(lat.gt.0.) then
          tile = 3
        else
          tile = 6
        endif
        return
      endif
      modlon = lon
      do while(modlon.lt.-pi/4d0)
        modlon = modlon + pi/2d0
      enddo
      do while(modlon.gt.+pi/4d0)
        modlon = modlon - pi/2d0
      enddo
      coslon = cos(modlon)
      tanlat = tan(lat)
      y = byg*atan(tanlat/(coslon*sqrt(2d0)))
      if(abs(y).le.1d0) then
c equatorial face
        x = byg*atan(tan(modlon)/sqrt(2d0))
c determine which face (1,2,4,5) we are on.  integer arithmetic
        tile = 1 + int((lon+1.25*pi)*2./pi)
c adjust x for roundoff errors
        if(abs(1d0-abs(x)).lt.1d-6) then
          if( (x.gt.0 .and. lon.lt.ctrlon(tile)) .or.
     &        (x.lt.0 .and. lon.gt.ctrlon(tile)) ) x = -x
        endif
        tile = tile +(tile/3) -5*(tile/5)
        if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
          tmpx = x
          tmpy = y
          x = -tmpy
          y = +tmpx
        endif
      else
c polar face
        bytanfac = sqrt(.5d0)/abs(tanlat)
        x = byg*atan(bytanfac*cos(lon))
        y = byg*atan(bytanfac*sin(lon))
        if(lat.gt.0.) then
          tile = 3
        else ! tile 6 x,y = tile 3 x,y flipped around the axis x+y=0
          tile = 6
          tmpx = x
          tmpy = y
          x = -tmpy
          y = -tmpx
        endif
      endif
      return
      end subroutine ll2csxy

      subroutine ll2csxy_vec(x_in,y_in,tile,ull,vll,uxy,vxy)
c
c Given latlon-oriented vector components ull,vll located at a
c position x_in,y_in on tile, return "vector" components uxy,vxy
c defining the directed great circle parallel to ull,vll.
c This great circle is a straight line in the X-Y space of this tile:
c         vxy*dX = uxy*dY
c See the beginning of this file for "capital" X,Y definitions.
c This routine could be made independent of the grid variant if
c position were specified as lon,lat or "capital" X,Y.
c
      implicit none
      real*8 :: x_in,y_in,ull,vll ! input
      integer :: tile ! input
      real*8 :: uxy,vxy ! output
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: x,y,rtxy,mag,bymag,utmp,vtmp
      x = x_in
      y = y_in
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        x = -y_in
        y = +x_in
      endif
      x = tan(g*x)*sqrt(2d0)
      y = tan(g*y)*sqrt(2d0)
      rtxy = sqrt(1d0+x*x+y*y)
      if(tile.eq.3 .or. tile.eq.6) then
        uxy = -y*ull - x*rtxy*vll
        vxy = +x*ull - y*rtxy*vll
        if(tile.eq.6) then ! tile 6 = tile 3 flipped around the axis x+y=0
          uxy = -uxy
          vxy = -vxy
        endif
      else
        uxy = (1.+x*x)*ull
        vxy = x*y*ull + rtxy*vll
        if(tile.eq.4 .or. tile.eq.5) then
          utmp = uxy; vtmp = vxy
          uxy = -vtmp
          vxy = +utmp
        endif
      endif
      mag = sqrt(uxy**2 + vxy**2)
      if (mag .gt. 1.d-8) then
         bymag=1./mag
         uxy = uxy*bymag
         vxy = vxy*bymag
      endif
      return
      end subroutine ll2csxy_vec

      subroutine e1e2(x,y,tile,e1,e2)
c
c this routine computes latlon-oriented basis vectors e1,e2 that
c are parallel to cubed-sphere gridlines at a position x,y on tile.
c e1(1:2) is parallel to the constant-y great circles and is
c  proportional to (/ cos(lat)*dlon/dx, dlat/dx /)
c e2(1:2) is parallel to the constant-x great circles and is
c  proportional to (/ cos(lat)*dlon/dy, dlat/dy /)
c
      implicit none
      real*8 :: x,y ! input
      integer :: tile ! input
      real*8, dimension(2) :: e1,e2
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: gx,gy,tangx,tangy,r2,tmpx,tmpy
      gx = g*x
      gy = g*y
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        tmpx = gx
        tmpy = gy
        gx = +tmpy
        gy = -tmpx
      endif
      tangx = tan(gx)
      tangy = tan(gy)
      r2 = 1d0/sqrt(1d0+2d0*(tangx**2+tangy**2))
      if(tile.eq.3 .or. tile.eq.6) then
        e1(1) = -tangy
        e1(2) = -tangx*r2
        e2(1) = tangx
        e2(2) = -tangy*r2
        if(tile.eq.6) then
          e1 = -e1
          e2 = -e2
        endif
      else
        e1(1)=1d0
        e1(2)=-2d0*tangx*tangy*r2
        e2(1)=0d0
        e2(2)=1d0
        if(tile.eq.4 .or. tile.eq.5) then
          e2(1)=1d0
          e2(2)=e1(2)
          e1(1)=0d0
          e1(2)=-1d0
        endif
      endif
      e1=e1/sqrt(sum(e1*e1))
      e2=e2/sqrt(sum(e2*e2))
      return
      end subroutine e1e2

      subroutine trigint(tile,x_in,y_in,a,cltcln,cltsln,slt)
c computes the integrals from 0 to x_in, 0 to y_in of
c area:                   a
c cos(lat)*cos(lon)*area: cltcln
c cos(lat)*sin(lon)*area: cltsln
c sin(lat)*area:          slt
c 
      implicit none
      integer :: tile
      real*8 :: x_in,y_in
      real*8 :: a,cltcln,cltsln,slt
      real*8 :: x,y,tmpx,tmpy,rtx,rty,atanx,atany
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      x = tan(g*x_in)*sqrt(2d0)
      y = tan(g*y_in)*sqrt(2d0)
      a = atan(x*y/sqrt(1.+x*x+y*y))
      tmpx = x; tmpy = y
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        x = +tmpy
        y = -tmpx
      elseif(tile.eq.6) then ! tile 6 = tile 3 flipped around the axis x+y=0
        x = -tmpy
        y = -tmpx
      endif
      rtx = sqrt(1.+x*x)
      rty = sqrt(1.+y*y)
      atanx = -.5*atan(x/rty)/rty
      atany = -.5*atan(y/rtx)/rtx
      if(tile.eq.3 .or. tile.eq.6) then ! polar tile
        cltsln = atanx
        cltcln = atany
        slt    = -(x*atany + y*atanx)
      else
        slt    = atanx
c longitude offsets; rotation on tiles 4/5
        if(tile.eq.1 .or. tile.eq.4) then
          cltsln = -atany
          cltcln = x*atany + y*atanx
        else
          cltsln = x*atany + y*atanx
          cltcln = atany
        endif
      endif
      if(tile.ge.4) slt = -slt
      return
      end subroutine trigint

      function gcdist(x1_in,x2_in,y1_in,y2_in)
c compute the great-circle distance between two x,y points.
      implicit none
      real*8 :: gcdist,x1_in,x2_in,y1_in,y2_in
      real*8 :: X1,X2,Y1,Y2,dot,vsqr1,vsqr2
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      X1 = tan(g*x1_in)
      Y1 = tan(g*y1_in)
      X2 = tan(g*x2_in)
      Y2 = tan(g*y2_in)
      dot   = 1d0+2d0*(X1*X2+Y1*Y2)
      vsqr1 = 1d0+2d0*(X1*X1+Y1*Y1)
      vsqr2 = 1d0+2d0*(X2*X2+Y2*Y2)
      gcdist = acos(dot/sqrt(vsqr1*vsqr2))
c note: in "capital" X,Y space,
c if Y1==Y2 gcdist = |atan(X2/sqrt(1+Y*Y))-atan(X1/sqrt(1+Y*Y))|
c if X1==X2 gcdist = |atan(Y2/sqrt(1+X*X))-atan(Y1/sqrt(1+X*X))|
      return
      end function gcdist

      END MODULE GEOM
