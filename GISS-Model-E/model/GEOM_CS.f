#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort_core(__LINE__,rc)
#define RETURN_(status) If (Present(rc)) rc=status; return

#include "rundeck_opts.h"
      MODULE GEOM
!@sum  GEOM contains spherical geometric variables and arrays
!@auth T. Clune
!@cont GEOM_CS

      USE MODEL_COM, only : IM,JM
      USE CONSTANT, only : radius, areag
      use FV_UTILS, only : abort_core
      use ESMF_Mod

      IMPLICIT NONE
      SAVE

!@var  i-index (longitude) of each grid cell
      INTEGER, ALLOCATABLE :: indx(:,:)
!@var  j-index (latitude) of each grid cell
      INTEGER, ALLOCATABLE :: jndx(:,:)

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
!@var  lat2d_corner latitude of corner point (radians)
      REAL*8, ALLOCATABLE :: lat2d_corner(:,:)
!@var  lon2d_corner longitude of mid point of primary grid box (radians)
      REAL*8, ALLOCATABLE :: lon2d_corner(:,:)
!@var ddx_ci, ddx_cj,ddy_ci ddy_cj 
      REAL*8, ALLOCATABLE :: ddx_ci(:,:),ddx_cj(:,:)
      REAL*8, ALLOCATABLE :: ddy_ci(:,:),ddy_cj(:,:)
      real*8 :: dloni,dlonj,dlati,dlatj,bydet

!@var  axyp,byaxyp area of grid box (+inverse) (m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: axyp, byaxyp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dxyp, bydxyp

      integer, allocatable, dimension(:) :: imaxj

!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, DIMENSION(JM) :: KMAXJ

!@var J_BUDG a mapping array that takes every grid point to the
!@+   zonal mean budget array
      integer, allocatable, dimension(:,:) :: J_BUDG
!@var j_0b, j_1b are the min/max zonal budget latitudes for this processor
      integer :: j_0b, j_1b


      CONTAINS

      SUBROUTINE GEOM_CS
!@sum  GEOM_CS Calculate spherical geometry for CS grid
!@auth T. Clune
      USE CONSTANT, only : RADIUS,PI,TWOPI,radian
      use fv_grid_tools_mod, only: agrid, corner_grid => grid, area
      use DOMAIN_DECOMP_ATM, only: grid, getDomainBounds, halo_update

      implicit none

      integer :: i0h, i1h, i0, i1
      integer :: j0h, j1h, j0, j1, j0s, j1s
      integer :: i,j

      call getDomainBounds(grid, J_STRT_HALO=j0h, J_STOP_HALO=j1h,
     &     J_STRT=j0, J_STOP=j1,
     &     I_STRT_HALO=i0h, I_STOP_HALO=i1h, 
     &     I_STRT=i0,I_STOP=i1, 
     &     J_STRT_SKP=j0s, J_STOP_SKP=j1s)

      write(*,*) "geom  i0h,i1h,i0,i1,j0h,j1h,j0,j1=",
     &     i0h,i1h,i0,i1,j0h,j1h,j0,j1

      allocate(indx(i0h:i1h, j0h:j1h))
      allocate(jndx(i0h:i1h, j0h:j1h))

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

      allocate(
     &    ddx_ci(i0h:i1h,j0h:j1h)
     &    ,ddx_cj(i0h:i1h,j0h:j1h)
     &    ,ddy_ci(i0h:i1h,j0h:j1h)
     &    ,ddy_cj(i0h:i1h,j0h:j1h))

      allocate(axyp(i0h:i1h, j0h:j1h))
      allocate(byaxyp(i0h:i1h, j0h:j1h))
      allocate(dxyp(j0h:j1h))
      allocate(bydxyp(j0h:j1h))

      allocate(imaxj(j0:j1))

      ! From FV CS grid
c      lon2d = agrid(:,:,1)
c      lat2d = agrid(:,:,2)

c      lon2d_corner = corner_grid(:,:,1)
c      lat2d_corner = corner_grid(:,:,2)

      call getLatLon_from_ESMFgrid(grid%esmf_grid,
     &         i0, i1, j0, j1, i0h, i1h, j0h, j1h,
     &         lat2d, lon2d, lon2d_corner, lat2d_corner)


      ! Fill in indices for grid cells, so we can
      ! halo-update them with our neighbors
      do j=j0,j1
      do i=i0,i1
        indx(i,j) = i
        jndx(i,j) = j
      end do    ! i
      end do    ! j

!    Update halo
      call halo_update(grid, lat2d)
      call halo_update(grid, lon2d)
      call halo_update(grid, lon2d_corner)
      call halo_update(grid, lat2d_corner)
      call halo_update(grid, indx)
      call halo_update(grid, jndx)


      lat2d_dg = lat2d/radian
      lon2d_dg = lon2d/radian
      sinlat2d = sin(lat2d)
      coslat2d = cos(lat2d)

c      AXYP = area(:,:)
c      BYAXYP = 1/AXYP
      axyp(:,:) =1.0  ! until problem with Initialize(fv) is fixed
      byaxyp(:,:) =1.0

      call set_j_budg   !called after lat2d_dg is initialized
      call set_wtbudg !sets area weights

      imaxj(:)=i1

      do j=j0,j1
      do i=i0,i1
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
      enddo
      enddo

      END SUBROUTINE GEOM_CS
!------------------------------------------------------------------------------
!
      subroutine getLatLon_from_ESMFgrid(ESMFgrid, 
     &         i0, i1, j0, j1,  i0h, i1h, j0h, j1h,
     &         lat2d, lon2d, lon2d_corner, lat2d_corner)

!@sum  Extracts the lat/lon information (cell center and corners) from the
!      ESMF grid subdomains.
!@fun  getLatLon_from_ESMFgrid
!
      use ESMF_Mod
!
      implicit none
!
!@var  i0, i1, j0, j1.      ESMF subdomain horizontal dimension limits.  
      integer         , intent(in) :: i0, i1, j0, j1

!@var  i0h, i1h, j0h, j1h.  ESMF subdomain horizontal dimension limits including the halo.  
      integer         , intent(in) :: i0h, i1h, j0h, j1h

!@var  ESMFgrid.            ESMF grid 
      type (ESMF_Grid), intent(in):: ESMFgrid
!
! INPUT/OUTPUT PARAMETERS:
!@var  lat2d                Contains latitudes for ESMF subdomain grid
      REAL*8 :: lat2d(i0h:i1h,j0h:j1h)

!@var  lat2d                Contains longitudes for ESMF subdomain grid
      REAL*8 :: lon2d(i0h:i1h,j0h:j1h)

!@var  lon2d_corner        Contains the longitudes at the corners of the 
!                           ESMF subdomain grid 
      REAL*8 :: lon2d_corner(i0h:i1h+1,j0h:j1h+1)

!@var  lat2d_corner         Contains the latitudes at the corners of the 
!                           ESMF subdomain grid 
      REAL*8 :: lat2d_corner(i0h:i1h+1,j0h:j1h+1)

! LOCAL VARIABLES:
      type (ESMF_Array),  pointer :: coords(:)
      real(ESMF_KIND_R8), pointer :: lons(:,:)
      real(ESMF_KIND_R8), pointer :: lats(:,:)
      type (ESMF_Array),  target  :: tarray(2)
      integer                     :: DIMS(2), STATUS
      integer			  :: counts(2), myPE
      type (ESMF_Array),  pointer     :: coordArray(:)
      real(ESMF_KIND_R8), pointer     :: lons3(:,:,:)
      real(ESMF_KIND_R8), pointer     :: lats3(:,:,:)

      call ESMF_GridGet(ESMFgrid, 
     &         horzRelLoc            = ESMF_CELL_CENTER, 
     &         globalCellCountPerDim = DIMS, RC=STATUS) 
      VERIFY_(STATUS)


!    Get ESMF number of subdomain grid points in each direction.  (variable counts)  
!    ESMF 2D subdomain variables are dimensioned as 1 to count(1) in the I-direction
!    and 1 to count(2) in the J-direction.   

      call ESMF_GridGetDELocalInfo(ESMFgrid, 
     &         horzRelLoc           = ESMF_CELL_CENTER, 
     &         myDE                 = myPE,
     &         localCellCountPerDim = counts,
     &         RC=STATUS)
      VERIFY_(STATUS)


!    Get lat/lons at cell center locations with Global I,J indices
!    Recall ESMF grid is indexed as IM x 6*JM and modelE is IM x JM

      coords=>tarray
      call ESMF_GridGetCoord(ESMFgrid, 
     &         horzRelloc  = ESMF_CELL_CENTER, 
     &         centerCoord = coords, rc=status)
      VERIFY_(STATUS)

!    Get Cell Center Locations
!    Longitudes

      call ESMF_ArrayGetData(coords(1), lons, RC=status)
      VERIFY_(STATUS)

!    Latitudes
      call ESMF_ArrayGetData(coords(2), lats, RC=status)
      VERIFY_(STATUS)

      lon2d(i0:i1,j0:j1) = lons(:,:)
      lat2d(i0:i1,j0:j1) = lats(:,:)
      
!    Get lat/lons at cell corners

      coordArray=>tarray

      call ESMF_GridGetCoord(ESMFgrid, 
     &          horzRelloc  = ESMF_CELL_CENTER, 
     &          cornerCoord = coordArray, rc=status)
      VERIFY_(STATUS)

!    We are counting on ESMF returning corners numbered
!    counter clockwise from the SW, so 1=SW, 2=SE, 3=NE, 4=NW

      call ESMF_ArrayGetData(coordArray(1), lons3, RC=status)
      VERIFY_(STATUS)

      call ESMF_ArrayGetData(coordArray(2), lats3, RC=status)
      VERIFY_(STATUS)

!    Extract the SW corner of each cell (covers the entire subdomain)
      lon2d_corner(i0:i1, j0:j1) = lons3(1,:,:)
      lat2d_corner(i0:i1, j0:j1) = lats3(1,:,:)

!    Extract only the NE corner of the subdomain.  ( two grid lines)

!    Right boundary of subdomain 
      lon2d_corner(i1+1,j0+1:j1+1) = lons3(3,counts(1),1:counts(2))
      lat2d_corner(i1+1,j0+1:j1+1) = lats3(3,counts(1),1:counts(2))

!    Top boundary of subdomain
      lon2d_corner(i0+1:i1,j1+1) = lons3(3,1:counts(1)-1,counts(2))
      lat2d_corner(i0+1:i1,j1+1) = lats3(3,1:counts(1)-1,counts(2))

!    SE: Extract only the lat/lons at the lower right corner of the subdomain

      lon2d_corner(i1+1,j0) = lons3(2,counts(1),1)
      lat2d_corner(i1+1,j0) = lats3(2,counts(1),1)

!    NW: Extract only the lat/lons at the upper left corner of the subdomain

      lon2d_corner(i0,j1+1) = lons3(4,1,counts(2))
      lat2d_corner(i0,j1+1) = lats3(4,1,counts(2))

      return

      end subroutine getLatLon_from_ESMFgrid
!------------------------------------------------------------------------------


      subroutine lonlat_to_ij(ll,ij)
c converts lon,lat=ll(1:2) into model i,j=ij(1:2)
c this is a dummy version to allow compilation of new river/lakes routine
c Will be fixed soon.
      implicit none
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: ij(2)
      ij(1) = 1
      ij(2) = 1
      return
      end subroutine lonlat_to_ij

      END MODULE GEOM


