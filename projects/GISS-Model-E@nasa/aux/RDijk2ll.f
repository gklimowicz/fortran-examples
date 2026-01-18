      program RDijk2ll_CS
!@sum On CS grid, convert i,j,k (with k=cube face index from 1..6)
!@+   to absolute lat-lon coordinates
!@auth Denis Gueyffier
      implicit none
      integer, parameter :: imt=90,jmt=90
      integer, parameter :: nrvr = 41 ! # of river mouths
      character*80 :: name,nameout,title,
     &     title1,title2,title3,title4,title5,title6,PROD
      real*8, parameter :: pi=3.141592653589793d0
      real*8, parameter :: shiftwest = 2.*3.141592653589793d0/36.
      real*8, parameter :: undef=-1.d30
      real*8, parameter :: radian = pi/180d0
      real*4, dimension(imt,jmt,6) :: down_lat,down_lon
      real*4, dimension(imt,jmt,6) :: down_lat_911,down_lon_911
      real*8, dimension(imt,jmt,6) :: lon2d,lat2d,lon2d_dg,lat2d_dg
      integer, dimension(imt,jmt,6) :: idown,jdown,kdown
      integer, dimension(imt,jmt,6) :: id911,jd911,kd911
      integer :: iu_RD,iu_TOPO,iu_MNAME,iu_MIJ,i,j,k,tile
      LOGICAL, dimension(imt,jmt) :: NODIR
      real*4 :: FOCEAN(imt,jmt,6)
      real*8:: x,y
      character*8,dimension(nrvr) :: namervr
      character*2,dimension(nrvr) :: mouthI,mouthJ
      character*1,dimension(nrvr) :: mouthK
      integer,dimension(nrvr) :: imouthI,imouthJ,imouthK
      real*4,dimension(nrvr) :: lat_rvr,lon_rvr
      iu_TOPO=30

c*    Read ocean fraction
         if (imt .eq. 32) then
            name="Z_CS32_4X5"
         else if (imt .eq. 90) then
            name="Z_C90fromZ1QX1N"
         endif
      call GetEnv ('PROD',PROD) !path to input file directory
      name = trim(PROD) // '/' // name 
      open(iu_TOPO,FILE=name,FORM='unformatted', STATUS='old')
      read(iu_TOPO) title,FOCEAN
      close(iu_TOPO)

c*    Read names of river mouths
      iu_MNAME=20
      if (imt .eq. 32) then
         name="mouthnames_DtoO_CS32"
      elseif (imt .eq. 90) then
         name="mouthnames_DtoO_CS90"
      endif
      name = trim(PROD) // '/' // name
      open(iu_MNAME,FILE=name,FORM='formatted', STATUS='old')
      READ (iu_MNAME,'(A8)') (namervr(I),I=1,nrvr) !Read mouths names
      write(*,*) namervr
      close(iu_MNAME)

c*    Read i,j,k coordinates of river mouths (k = index of cube face)
      if (imt .eq. 32) then
         name="mouthij_DtoO_CS32"   
      elseif (imt .eq. 90) then
         name="mouthij_DtoO_CS90"
      endif
      name = trim(PROD) // '/' // name
      open(iu_MIJ,FILE=name,FORM='formatted', STATUS='old')
      READ (iu_MIJ,'(A2,1X,A2,1X,A1)') (      
     &       mouthI(I),mouthJ(I),mouthK(I), I=1,nrvr)
      close(iu_MIJ)

c*     conversion char to int
      do i=1,nrvr
          read(mouthI(i),'(I2)') imouthI(i)
          read(mouthJ(i),'(I2)') imouthJ(i)
          read(mouthK(i),'(I1)') imouthK(i)
          write(*,*) imouthI(i),imouthJ(i),imouthK(i)
      enddo   
    
      iu_RD=20
      if (imt .eq. 32) then
         name="RDtoO.CS32.INT"
         nameout="RDdistocean_CS32.bin"
      elseif (imt .eq. 90) then
         name="RDtoO.CS90.INT"
         nameout="RDdistocean_CS90.bin"
      endif
      name = trim(PROD) // '/' // name
      
      write(*,*) name,imt,jmt
      
c* read i,j,k coordinates of downstream cells
      
      open( iu_RD, FILE=name,FORM='unformatted',
     &     STATUS='unknown')
      
      read(iu_RD) title,idown,jdown,kdown
      write(*,*) title
      close(iu_RD)
      write(*,*) "read RDijk2ll_CS"

      do tile=1,6
         do j=1,imt
            do i=1,imt
               x = -1d0 + 2d0*(dble(i)-.5d0)/imt
               y = -1d0 + 2d0*(dble(j)-.5d0)/imt
               call csxy2ll(x,y,tile,lon2d(i,j,tile),lat2d(i,j,tile))
               lon2d(i,j,tile) = lon2d(i,j,tile)-shiftwest
               lat2d_dg(i,j,tile) = lat2d(i,j,tile)/radian
               lon2d_dg(i,j,tile) = lon2d(i,j,tile)/radian
               if (lon2d_dg(i,j,tile) .le. -180.d0) then
                  lon2d_dg(i,j,tile)=lon2d_dg(i,j,tile)+360.d0
               endif
            enddo
         enddo
      enddo

      write(201,*) lon2d_dg 
      write(201,*) lat2d_dg
c*    convert i,j,k coordinates of river mouths to absolute lat-lon coordinates    
      do k=1,nrvr
         lat_rvr(k)=lat2d_dg(imouthI(k),imouthJ(k)
     &        ,imouthK(k))
         lon_rvr(k)=lon2d_dg(imouthI(k),imouthJ(k)
     &        ,imouthK(k))
      enddo

c*    convert i,j,k coordinates of downstream cells to absolute lat-lon coordinates    
      do k=1,6
        do j=1,jmt
          do i=1,imt
            if (FOCEAN(i,j,k) .lt. 1.e-6) then
              down_lat(i,j,k)=lat2d_dg(idown(i,j,k),jdown(i,j,k)
     &            ,kdown(i,j,k))
              down_lon(i,j,k)=lon2d_dg(idown(i,j,k),jdown(i,j,k)
     &            ,kdown(i,j,k))

c*    dummy emergency directions
              down_lat_911(i,j,k)=down_lat(i,j,k)
              down_lon_911(i,j,k)=down_lon(i,j,k)
            else
              down_lat(i,j,k)=undef
              down_lon(i,j,k)=undef
              down_lat_911(i,j,k)=undef
              down_lon_911(i,j,k)=undef
            end if

          enddo
        enddo
      enddo

c*    output everything 
      open(iu_RD,FILE=nameout,FORM='unformatted',
     &        STATUS='unknown')

      if (imt .eq. 32) then
         title1="river directions from dist. to ocean, CS32, April 09"
      elseif (imt .eq. 90) then
         title1="river dir. from dist. to ocean, CS90, August 2010"
      endif

      title2="Named River Mouths:"
      title3="Latitude of downstream river direction box"
      title4="Longitude of downstream river direction box"
      title5="Latitude of emergency downstream river direction box"
      title6="Longitude of emergency downstream river direction box"

      write(iu_RD) title1   
      write(iu_RD) title2,nrvr,namervr(1:nrvr),lat_rvr(1:nrvr)
     *     ,lon_rvr(1:nrvr)
      write(iu_RD) title3,down_lat
      write(iu_RD) title4,down_lon
      write(iu_RD) title5,down_lat_911
      write(iu_RD) title6,down_lon_911
      close(iu_RD)

      write(*,*) "wrote RD file"

 200  format(4(1X,f8.3))

      end program RDijk2ll_CS
c*

      subroutine csxy2ll(x,y,tile,lon,lat)
c converts x,y,tile to lon,lat (radians)
c This routine places the center of tile 1 at the IDL.
      implicit none
      real*8 :: x,y ! input
      integer :: tile ! input
      real*8 :: lon,lat ! output
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8,parameter :: pi=3.141592653589793d0
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


