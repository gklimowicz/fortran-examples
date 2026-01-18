#ifdef CUBED_SPHERE
      subroutine bilin_ll2cs_vec(grid,uin_glob,vin_glob,
     &     uout_loc,vout_loc,ims,jms,periodic)
!@sum Bilinearly interpolate a vector field from latlon grid to cubed-sphere
!@+   The vector field is defined on the sphere (Earth) and takes values in the
!@+   plane locally tangent to sphere at each point 
!@+   input field on latlon grid is global 
!@+   interpolation is performed in the local cs domain
!@+   output field is local to each cs domain 
!@+
!@+  (uin12,vin12)      (uin22,vin22)
!@+   x-------------------x  <--lat2
!@+   |                   |
!@+   |       +           |
!@+   | (lon2d, lat2d)    |
!@+   |                   |
!@+   x___________________x  <--lat1
!@+ (uin11,vin11)     (uin21,vin21)
!@+   ^                   ^
!@+   |                   |
!@+  lon1               lon2
!@+
!@var uin_glob input vector x-component (lat=cste) -- global
!@var vin_glob input vector y-component (lon=cste) -- global
!@var uout_loc output vector x-component (lat=cste) --local
!@var vout_loc output vector y-component (lon=cste) --local
!@auth Denis Gueyffier
      
      use regrid_com
      use geom, only : lat2d, lon2d, lat2d_dg, lon2d_dg
      implicit none
      type (dist_grid), intent(in) :: grid
      integer, intent(in) :: ims,jms
      real*8, intent(inout) :: 
c     &       uout_loc(:,:),vout_loc(:,:)
     &     uout_loc(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO),
     &     vout_loc(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)
      real*8, intent(in) :: 
     &     uin_glob(ims,jms),
     &     vin_glob(ims,jms)
      real*8 :: u11,u12,u21,u22
      real*8 :: v11,v12,v21,v22
      real*8 :: lon1,lon2,lat1,lat2, lon_curr, lat_curr
      real*8 :: dnm, pi, dlat_dg, dlon_dg
      integer :: i,j, i_lon,j_lat, im_lon, jm_lat, i1,i2, j1,j2, jeq,
     &      i0h,i1h,j0h,j1h
      logical, optional, intent(in) :: periodic

      i0h=grid%I_STRT_HALO
      i1h=grid%I_STOP_HALO
      j0h=grid%J_STRT_HALO
      j1h=grid%J_STOP_HALO

      pi = 4.0d0*atan(1.d0)

      dlat_dg=180./REAL(JMS)                   ! even spacing (default)
      IF (JMS.eq.46) dlat_dg=180./REAL(JMS-1)   ! 1/2 box at pole for 4x5
      dlon_dg = 360./dble(ims)
      jeq=0.5*(1+jms)
c      
c***  loop on CS points in local domain
c
      do j=grid%j_strt,grid%j_stop
         do i=grid%i_strt,grid%i_stop

c     find (ilon,jlat) indices of latlon cell which contains the center of the (i,j) CS cell  
            if (present(periodic)) then
               I_lon = ims/2  + (lon2d_dg(i,j)+.01)/dlon_dg
            else
               I_lon = ims/2 + 1  + (lon2d_dg(i,j)+.01)/dlon_dg
            endif
            J_lat = jms/2 + 1  + (lat2d_dg(i,j)+.01)/dlat_dg
            
c     find the indices of the 4 latlon cells surrounding the current CS cell
c     first reconstruct latlon coordinates of current latlon cell
            if (i_lon .ne. 1) then 
               lon_curr = -180.+(i_lon-0.5)*dlon_dg
            else
               lon_curr = -180.+0.5*dlon_dg
            endif
            if (j_lat .ne. 1 .and. j_lat .ne. jms) then
               lat_curr = dlat_dg*(j_lat - jeq)
            elseif (j_lat .eq. 1) then
               lat_curr = -90.
            elseif (j_lat .eq. jms) then
               lat_curr = 90
            endif

            if (lon2d_dg(i,j) .gt. lon_curr) then
               lon1 = lon_curr
               I1=i_lon
            else
               lon1 = lon_curr - dlon_dg
               I1 = i_lon - 1
            endif
            lon2 = lon1 + dlon_dg
            I2 = I1+1

            if (lat2d_dg(i,j) .gt. lat_curr) then
               lat1 = lat_curr
               J1 = j_lat
            else
               lat1 = lat_curr - dlat_dg
               J1 = j_lat -1
            endif
            lat2 = lat1 + dlat_dg
            J2 = J1 + 1


c            write(*,100) "lon1 lon lon2",lon1,lon2d_dg(i,j),lon2
c            write(*,100) "lat1 lat lat2",lat1,lat2d_dg(i,j),lat2
c 100        format(A,3(1X,F8.2))

            if ( I1 .ne. 0 .and. I2 .ne. ims+1 
     &           .and. J1 .ne. 0 .and. J2. ne. jms+1) then
            
               u11 = uin_glob(I1,J1)
               u21 = uin_glob(I2,J1)
               u12 = uin_glob(I1,J2)
               u22 = uin_glob(I2,J2)
               
               v11 = vin_glob(I1,J1)
               v21 = vin_glob(I2,J1)
               v12 = vin_glob(I1,J2)
               v22 = vin_glob(I2,J2)
               
c     Boundary conditions. As we do not have halo cells for global latlon arrays, 
c     we do boundary conditions by hand
            else

               if ( I1 .eq. 0) then
                  u11 = uin_glob(ims,J1)
                  u12 = uin_glob(ims,J2)
                  
                  v11 = vin_glob(ims,J1)
                  v12 = vin_glob(ims,J2)
               endif
               if ( I2 .eq. ims+1) then
                  u21 = uin_glob(1,J1)
                  u22 = uin_glob(1,J2)
                  
                  v21 = vin_glob(1,J1)
                  v22 = vin_glob(1,J2)
               endif
               if ( J1 .eq. 0) then
                  u11 = uin_glob(I1,jms)
                  u21 = uin_glob(I2,jms)
                  
                  v11 = vin_glob(I1,jms)
                  v21 = vin_glob(I2,jms)
               endif
               if ( J2 .eq. jms+1) then
                  u12 = uin_glob(I1,1)
                  u22 = uin_glob(I2,1)
                  
                  v12 = vin_glob(I1,1)
                  v22 = vin_glob(I2,1)
               endif
            endif

            dnm=1.d0/(dlon_dg*dlat_dg)
            
            uout_loc(i,j)=dnm*( 
     &           u11*( lon2-lon2d_dg(i,j) )*( lat2-lat2d_dg(i,j) ) 
     &           + u21*( lon2d_dg(i,j)-lon1 )*( lat2-lat2d_dg(i,j) ) 
     &           + u12*( lon2-lon2d_dg(i,j) )*( lat2d_dg(i,j)-lat1 )
     &           + u22*( lon2d_dg(i,j)-lon1 )*( lat2d_dg(i,j)-lat1 ) 
     &           )
            
            vout_loc(i,j)=dnm*( 
     &           v11*( lon2-lon2d_dg(i,j) )*( lat2-lat2d_dg(i,j) ) 
     &           + v21*( lon2d_dg(i,j)-lon1 )*( lat2-lat2d_dg(i,j) ) 
     &           + v12*( lon2-lon2d_dg(i,j) )*( lat2d_dg(i,j)-lat1 )
     &           + v22*( lon2d_dg(i,j)-lon1 )*( lat2d_dg(i,j)-lat1 ) 
     &           )

c            write(690+grid%gid+1,200) lon2d_dg(i,j),lat2d_dg(i,j),
c     &           uout_loc(i,j),vout_loc(i,j)
c 200        format(4(1X,f8.3))

c            write(630+grid%gid+1,200) lon2d_dg(i,j),lat2d_dg(i,j),
c     &           vout_loc(i,j)

         enddo
      enddo
      
      end subroutine bilin_ll2cs_vec
c*

c      subroutine parallel_bilin_latlon_B_2_CS_A(gridin,gridout,
c     &     ain_loc,aout_glob,IM_LL,JM_LL)
c!@sum linearly interpolate from Latlon B-grid to the CS A-grid
c!@auth Denis Gueyffier
c      use regrid_com
c      use constant, only : radian
c      use resolution, only : im_cs=>im,jm_cs=>jm
c      use geom, only : lat2d, lon2d, lat2d_dg, lon2d_dg,lonlat_to_ij
c      use domain_decomp_atm, only : atm_pack=>pack_data
c      use domain_decomp_1d, only : halo_update,pack_data
c      implicit none
c      include 'mpif.h'
c      type (dist_grid), intent(in) :: gridin,gridout
c      real*8 :: ain_loc(IM_LL,
c     &     gridin%J_STRT_HALO:gridin%J_STOP_HALO),
c     &     aout_glob(IM_CS,JM_CS,6)
c      real*8, allocatable :: lat_glob(:,:,:),lon_glob(:,:,:)
c      real*8, allocatable :: ain_glob(:,:)
c      real*4, allocatable :: a4glob(:,:)
c      real*8 :: a11,a12,a21,a22
c      real*8 :: lon1,lon2,lat1,lat2,lon_cs,lat_cs,lat_ll,lon_ll,fjeq
c      real*8 :: dnm, pi, dlat_dg, dlon_dg,lat_curr,lon_curr
c      real*8 :: ll(2)
c      integer :: im_ll,jm_ll,ierr,nij
c      integer :: i,j,k,i_lon,j_lat,im_lon,jm_lat,i1,i2,j1,j2, 
c     &     i0_ll,i1_ll,j0_ll,j1_ll,i0_cs,i1_cs,j0_cs,j1_cs,
c     &     j0h_ll,j1h_ll
c      integer :: ij(2),ij1(2),ij2(2),ij3(2),ij4(2),i_cs,j_cs
c      logical :: in1,in2,in3,in4
c      character*80 :: title,name
c      character*1 :: ci
c
c      i0_cs=gridout%I_STRT
c      i1_cs=gridout%I_STOP
c      j0_cs=gridout%J_STRT
c      j1_cs=gridout%J_STOP
c
c      i0_ll=gridin%I_STRT
c      i1_ll=gridin%I_STOP
c      j0h_ll=gridin%J_STRT_HALO
c      j1h_ll=gridin%J_STOP_HALO
c
c      aout_glob(:,:,:) = 0.
c
c      dlat_dg=180./REAL(JM_LL)                   ! even spacing (default)
c      IF (JM_LL.eq.46) dlat_dg=180./REAL(JM_LL-1)   ! 1/2 box at pole for 4x5
c      dlon_dg = 360./dble(IM_LL)
c      fjeq=0.5*(1+JM_LL)
c
c      call halo_update(gridin,ain_loc)
c
c      allocate(lat_glob(IM_CS,JM_CS,6),
c     &     lon_glob(IM_CS,JM_CS,6))
c
cc***  TODO: lon_glob/lat_glob should be precomputed and moved to init_icedyn
c      call atm_pack(gridout,lon2d_dg,lon_glob)
c      call atm_pack(gridout,lat2d_dg,lat_glob)
c
c      nij=im_cs*jm_cs*6
c
c      call MPI_BCAST( lat_glob, nij, MPI_DOUBLE_PRECISION,
c     *     0, MPI_COMM_WORLD, ierr )
c      call MPI_BCAST( lon_glob, nij, MPI_DOUBLE_PRECISION,
c     *     0, MPI_COMM_WORLD, ierr )
cc*
c
cc***  Loop on points in global target (CS) domain
cc***
cc***  Loop below could be optimized by pre-computing the subset of CS
cc***  indices corresponding to the latitude band of the 
cc***  current latlon processor. Instead of looping on the whole set of indices
cc***  we would loop on the subset
c      do k=1,6
c      do j=1,JM_CS
c         do i=1,IM_CS
c
cc     find (ilon,jlat) indices of latlon cell containing the current CS cell center
c            i_lon = nint(.5*(im_ll + 1)  + lon_glob(i,j,k)/dlon_dg )
c            j_lat = nint(.5*(jm_ll + 1)  + lat_glob(i,j,k)/dlat_dg )
c
cc     find the 4 corner velocities on latlon B-grid 
cc     first reconstruct latlon coordinates of current latlon cell center
c            if (i_lon .ne. 1) then 
c               lon_curr = -180.+(i_lon-0.5)*dlon_dg
c            else
c               lon_curr = -180.+0.5*dlon_dg
c            endif
c            if (j_lat .ne. 1 .and. j_lat .ne. jm_ll) then
c               lat_curr = dlat_dg*(j_lat - fjeq)
c            elseif (j_lat .eq. 1) then
c               lat_curr = -90.
c            elseif (j_lat .eq. jm_ll) then
c               lat_curr = 90.
c            endif
c
c            lon1 = lon_curr - 0.5*dlon_dg
c            I1 = i_lon - 1
c            lon2 = lon1 + dlon_dg
c            I2 = I1+1
c            lat1 = lat_curr - 0.5*dlat_dg
c            J1 = j_lat -1
c            lat2 = lat1 + dlat_dg
c            J2 = J1 + 1
c
cc            write(*,100) "lon1 lon lon2",lon1,lon_glob(i,j,k),lon2
cc            write(*,100) "lat1 lat lat2",lat1,lat_glob(i,j,k),lat2
cc 100        format(A,3(1X,F8.2))
c
c            if (J1 .ge. j0h_ll .and. J2 .le. j1h_ll) then
c
c               if (I1 .gt. 0 .and. I2 .le. im_ll)  then
c                  a11 = ain_loc(I1,J1)
c                  a21 = ain_loc(I2,J1)
c                  a12 = ain_loc(I1,J2)
c                  a22 = ain_loc(I2,J2)
c
c
cc     periodic boundary conditions in x direction
c               else
c                  if ( I1 .le. 0) then
c                     a11 = ain_loc(im_ll,J1) !periodicity
c                     a12 = ain_loc(im_ll,J2) !periodicity
c                     a21 = ain_loc(1,J1) !periodicity
c                     a22 = ain_loc(1,J2) !periodicity
c                  endif
c                  if ( I2 .gt. im_ll) then
c                     a21 = ain_loc(1,J1) !periodicity
c                     a22 = ain_loc(1,J2) !periodicity
c                     a11 = ain_loc(I1,J1)
c                     a12 = ain_loc(I1,J2)
c                  endif
c               endif
c                           
c               dnm=1.d0/(dlon_dg*dlat_dg)
c
c               if (J1 .le. j0h_ll .or. J2 .ge. j1h_ll .and. J1 .ne. 0) 
c     &              dnm=0.5*dnm
c               if (J1 .eq. 0) then 
c                  a11=a12;a21=a22;dnm=2*dnm
c               endif
c               
c               aout_glob(i,j,k)=dnm*( 
c     &              a11*(lon2-lon_glob(i,j,k))*(lat2-lat_glob(i,j,k))
c     &              +a21*(lon_glob(i,j,k)-lon1)*(lat2-lat_glob(i,j,k))
c     &              +a12*(lon2-lon_glob(i,j,k))*(lat_glob(i,j,k)-lat1)
c     &              +a22*(lon_glob(i,j,k)-lon1)*(lat_glob(i,j,k)-lat1) 
c     &              )
c            endif
c
c         enddo
c      enddo
c      enddo
c      
c      call sumxpe(aout_glob)
c
c
c      if (am_i_root()) then
c         do k=1,6
c            do j=1,JM_CS
c               do i=1,IM_CS
c                  write(1101,*) 
c     &                 lon_glob(i,j,k),lat_glob(i,j,k),
c     &                 aout_glob(i,j,k)
c               enddo
c            enddo
c         enddo
c      endif
c
c      deallocate(lat_glob,lon_glob)
c
c      end subroutine parallel_bilin_latlon_B_2_CS_A
c*
#endif

      subroutine init_xgrid_zonal(x2grids,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget)
!@sum Initialize exchange grid for constant latitude band zonal means (no zigzag)
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'

      type (x_2grids), intent(inout) :: x2grids
      integer, intent(in) :: imsource,jmsource,imtarget,jmtarget
      integer, intent(in) :: ntilessource,ntilestarget
      real*8,  allocatable, dimension(:) :: xgrid_area  !local variable
      integer, allocatable, dimension(:) :: tile        !local variable
      integer, allocatable, dimension(:,:) :: ijcub     !local variable
      integer, allocatable, dimension(:,:) :: ijlatlon  !local variable
      integer :: ncells                                 !local variable
      integer :: maxkey                                 !local variable
      integer :: mytile
      integer :: status,fid,n,vid,ikey,jlat,nc,nc2,
     *     ierr,icc,jcc,jlatm,ic,jc
      integer :: itile,j,idomain,iic,jjc,index,indexc
      character*200 :: exchfile
      real*8 :: checkarea
      character(len=10) :: imch,jmch,icch,jcch 


      x2grids%imsource=imsource
      x2grids%jmsource=jmsource
      x2grids%ntilessource=ntilessource
      x2grids%imtarget=imtarget
      x2grids%jmtarget=jmtarget
      x2grids%ntilestarget=ntilestarget

      write(imch,'(i10)') imsource
      write(jmch,'(i10)') jmsource
      write(icch,'(i10)') imtarget
      write(jcch,'(i10)') jmtarget
      imch=trim(adjustl(imch))
      jmch=trim(adjustl(jmch))
      icch=trim(adjustl(icch))
      jcch=trim(adjustl(jcch))
      exchfile="remap"//trim(imch)//"-"//trim(jmch)
     *     //"C"//trim(icch)//"-"//trim(jcch)//".nc"
      write(*,*) "filename=",exchfile
      imch=trim(adjustl(imch))
      jmch=trim(adjustl(jmch))
      icch=trim(adjustl(icch))
      jcch=trim(adjustl(jcch))
      
      exchfile="remap"//trim(imch)//"-"//trim(jmch)
     *     //"C"//trim(icch)//"-"//trim(jcch)//".nc"
      write(*,*) "-->filename=",exchfile

c      exchfile="exexch.nc"

c      
c Read weights
c
      status = nf_open(trim(exchfile),nf_nowrite,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN REMAP FILE"
 
      status = nf_inq_dimid(fid,'ncells',vid)
      status = nf_inq_dimlen(fid,vid,ncells)

c
c Allocate arrays with size depending on ncells
c
      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))
      
      if (AM_I_ROOT()) then   
         status = nf_inq_varid(fid,'xgrid_area',vid)
         status = nf_get_var_double(fid,vid,xgrid_area)

         
         status = nf_inq_varid(fid,'tile1',vid)
         status = nf_get_var_int(fid,vid,tile)
         
         status = nf_inq_varid(fid,'tile1_cell',vid)
         status = nf_get_var_int(fid,vid,ijcub)
         
         status = nf_inq_varid(fid,'tile2_cell',vid)
         status = nf_get_var_int(fid,vid,ijlatlon)
         
         status = nf_close(fid)
      endif

c
c     Broadcast x_grid area and indices to all procs
c
      nc2=2*ncells
     
c      write(6,*) "BEFMPI ROLLED"

      call MPI_BCAST( xgrid_area, ncells, MPI_DOUBLE_PRECISION,
     *     0, MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijcub, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijlatlon, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( tile, ncells, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

c      write(6,*) "AFT BCAST ROLLED"
c****
c*     ikey=ikey+1
c*     index=(iic-1)*jcmax+jjc
c*     array index_key : ikey->index
c*     array jlat_key : ikey->jlat
c*     array xarea_key : ikey->xgrid_area
c****

      
c     calculate maxkey = size of azonal 
 
      ikey=1
      jlatm=x2grids%jmtarget
      do jlat=1,jlatm
         do n=1,ncells
            itile=tile(n)
            if (itile .eq. mytile) then
               j=ijlatlon(2,n)
               if (j .eq. jlat) then
                  ikey=ikey+1
               endif
            endif
         enddo
      enddo
      
      maxkey=ikey-1
      x2grids%xgrid%maxkey=maxkey
      x2grids%xgrid%ncells=ncells
    
      allocate(x2grids%xgrid%index_key(maxkey))
      allocate(x2grids%xgrid%jlat_key(maxkey))
      allocate(x2grids%xgrid%xarea_key(maxkey))

c      write(6,*) "AFT ALOC ROLLED, maxkey",maxkey

      x2grids%xgrid%index_key(:)=0
      x2grids%xgrid%jlat_key(:)=0
      x2grids%xgrid%xarea_key(:)=0.d0

      checkarea=0.d0
      ikey=1
      ic=x2grids%imsource
      jc=x2grids%jmsource
      do iic=1,ic
         do jjc=1,jc
            index=(iic-1)*jc+jjc
            do n=1,ncells
               itile=tile(n)
               if (itile .eq. mytile) then
                  icc=ijcub(1,n)
                  jcc=ijcub(2,n)
                  indexc=(icc-1)*jc+jcc
c               write(6,*) "index indexc",index,indexc
                  if ( index .eq. indexc ) then
                     jlat=ijlatlon(2,n)
                     x2grids%xgrid%index_key(ikey)=index
                     x2grids%xgrid%jlat_key(ikey)=jlat
                     x2grids%xgrid%xarea_key(ikey)=xgrid_area(n)    
c                     write(6,*) "ikey=",ikey
                     ikey=ikey+1
                     checkarea=checkarea+xgrid_area(n)
                  endif
               endif
            enddo
         enddo
      enddo
c      write(*,*) "ikey=",ikey-1
c      write(*,*) "->checkarea=",checkarea

      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)
      
      end subroutine init_xgrid_zonal
c****


      subroutine zonalmean_band(x2grids,value,iinput,jinput,
     *     zonal_mean,area_latband)

!@sum  Calculate zonal mean using latitute band algorithm (not zigzags on budget grid)
!@auth Denis Gueyffier

      use regrid_com
      implicit none
      include 'netcdf.inc'
      type (x_2grids), intent(in) :: x2grids
      integer, intent(in) :: iinput,jinput
      real*8 :: value,zonal_mean(x2grids%jmtarget),
     &     area_latband(x2grids%jmtarget),
     &     zstore,astore
      integer :: ikey,jlat,current_index,index

      if ( ( (iinput .le. x2grids%ie) .and. (iinput .ge. x2grids%is)) 
     &     .and.
     &     ( (jinput .le. x2grids%je) .and. (jinput .ge. x2grids%js)) ) 
     &     then

      current_index=x2grids%jmsource*(iinput-1)+jinput

      do ikey=1,x2grids%xgrid%maxkey
         index=x2grids%xgrid%index_key(ikey)            
         if (index .eq. current_index) then
            jlat=x2grids%xgrid%jlat_key(ikey) 
            zstore=zonal_mean(jlat)
            zonal_mean(jlat)=zstore
     &           +x2grids%xgrid%xarea_key(ikey)*value
            astore=area_latband(jlat)
            area_latband(jlat)=astore
     &           +x2grids%xgrid%xarea_key(ikey)
         endif
      enddo

      endif

      end subroutine zonalmean_band
c*

      subroutine parallel_regrid(x2grids,tsource,   
     &     ttarget,atarget)

!@sum  Parallel regriding from source to target grid
!@+    x_2grids is used to determine type of source and target 
!@+    grids (cubed-sphere or lat-lon)
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tsource(x2grids%isd:x2grids%ied,
     &     x2grids%jsd:x2grids%jed)
      real*8 :: ttarget(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
     &     ,atarget(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
      integer :: n,icub,jcub,i,j,itile,ilon,jlat,ikey
      character*120:: ofi
      integer ::  status, fid, vid
       
      if ( (x2grids%ntilessource .eq. 6) .and.   !cs2ll
     &     (x2grids%ntilestarget .eq. 1) ) then   
         
         atarget(:,:,1) = 0.d0
         ttarget(:,:,1) = 0.d0

         do ikey=1,x2grids%xgrid%maxkey
            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)
            
            atarget(ilon,jlat,1) = atarget(ilon,jlat,1) + 
     &           x2grids%xgrid%xarea_key(ikey)
            ttarget(ilon,jlat,1) = ttarget(ilon,jlat,1) + 
     &           x2grids%xgrid%xarea_key(ikey)
     &           * tsource(icub,jcub) 
         enddo
         
c
c     sum all contributions
c     
         call SUMXPE(ttarget)
         call SUMXPE(atarget)
         
c     
c     root proc section
c     
         if (AM_I_ROOT()) then   
            ttarget(:,:,1) = ttarget(:,:,1)/atarget(:,:,1)
            
c         ofi='tstout.nc'
c         status = nf_open(trim(ofi),nf_write,fid)
c         if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
c         status = nf_inq_varid(fid,'lwup_sfc',vid)
c         write(*,*) NF_STRERROR(status)
c         status = nf_put_var_double(fid,vid,ttarget)
c         write(*,*) "STATUS",NF_STRERROR(status),"<<"
c         status = nf_close(fid)
         endif

      else if ( (x2grids%ntilessource .eq. 1) .and. !ll2cs
     &        (x2grids%ntilestarget .eq. 6) ) then
         atarget(:,:,:) = 0.d0
         ttarget(:,:,:) = 0.d0

         do ikey=1,x2grids%xgrid%maxkey
            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)
            
            atarget(icub,jcub,itile) = atarget(icub,jcub,itile) + 
     &           x2grids%xgrid%xarea_key(ikey)
            ttarget(icub,jcub,itile) = ttarget(icub,jcub,itile) + 
     &           x2grids%xgrid%xarea_key(ikey)
     &           * tsource(ilon,jlat) 
         enddo
         
c
c     sum all contributions
c     
         call SUMXPE(ttarget)
         call SUMXPE(atarget)  

c     rescale by local area
         if (AM_I_ROOT()) then
            ttarget(:,:,:) = ttarget(:,:,:)/atarget(:,:,:)       
         endif
         
      endif                     !ll2cs

      end subroutine parallel_regrid
c*


      subroutine parallel_regrid_wt(x2grids,wsource,missing,tsource,   
     &     ttarget,atarget)

!@sum  same as parallel regrid but we multiply by weight function wsource
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tsource(x2grids%isd:x2grids%ied,
     &     x2grids%jsd:x2grids%jed)
      real*8 :: wsource(x2grids%isd:x2grids%ied,
     &     x2grids%jsd:x2grids%jed)
      real*8 :: ttarget(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
     &     ,atarget(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
      real*8, intent(in) :: missing
      integer :: n,icub,jcub,i,j,k,itile,ilon,jlat,ikey
      character*120:: ofi
      integer ::  status, fid, vid

       
      if ( (x2grids%ntilessource .eq. 6) .and.   !cs2ll
     &     (x2grids%ntilestarget .eq. 1) ) then   
         
         atarget(:,:,1) = 0.d0
         ttarget(:,:,1) = 0.d0


         do ikey=1,x2grids%xgrid%maxkey
            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)
            
            atarget(ilon,jlat,1) = atarget(ilon,jlat,1) + 
     &           wsource(icub,jcub)*x2grids%xgrid%xarea_key(ikey)
            ttarget(ilon,jlat,1) = ttarget(ilon,jlat,1) + 
     &           x2grids%xgrid%xarea_key(ikey)
     &           *wsource(icub,jcub)*tsource(icub,jcub) 
         enddo

c     sum all contributions
c     
c         write(*,*) "bef sumxpe"
         call SUMXPE(ttarget)
         call SUMXPE(atarget)
c         write(*,*) "aft sumxpe"
c     
c     root proc section
c     
         if (AM_I_ROOT()) then   
          
            do j=1,x2grids%jmtarget
               do i=1,x2grids%imtarget
                  if (atarget(i,j,1) .le. 1.e-10) then
                     ttarget(i,j,1)=missing
                  else
                     ttarget(i,j,1) = ttarget(i,j,1)/atarget(i,j,1)
                  endif
               enddo
            enddo
          
           
c            ofi='tstout.nc'
c            write(*,*) ofi
c            status = nf_open(trim(ofi),nf_write,fid)
c            if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
c            status = nf_inq_varid(fid,'lwup_sfc',vid)
c            write(*,*) NF_STRERROR(status)
c            status = nf_put_var_double(fid,vid,ttarget)
c            write(*,*) "STATUS",NF_STRERROR(status),"<<"
c            status = nf_close(fid)
         endif

      else if ( (x2grids%ntilessource .eq. 1) .and. !ll2cs
     &        (x2grids%ntilestarget .eq. 6) ) then
         atarget(:,:,:) = 0.d0
         ttarget(:,:,:) = 0.d0
c         write(*,*) "IN regrid wt"

c         write(*,*) "maxkey regx=",x2grids%xgrid%maxkey
         do ikey=1,x2grids%xgrid%maxkey
            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)

c            write(*,*) "icub jcub itile xarea=",icub,jcub,itile,
c     &           x2grids%xgrid%xarea_key(ikey)
            if ( tsource(ilon,jlat) .gt. missing ) then
               atarget(icub,jcub,itile) = atarget(icub,jcub,itile)+
     &              wsource(ilon,jlat)*x2grids%xgrid%xarea_key(ikey)
               
               ttarget(icub,jcub,itile) = ttarget(icub,jcub,itile)+
     &              x2grids%xgrid%xarea_key(ikey)
     &              *wsource(ilon,jlat)*tsource(ilon,jlat) 
            endif
         enddo
     

c
c     sum all contributions
c    
   
         call SUMXPE(ttarget)
         call SUMXPE(atarget)

c     
c     root proc section
c     
         if (AM_I_ROOT()) then
               do k=1,6
               do j=1,x2grids%jmtarget
                  do i=1,x2grids%imtarget
                     if (atarget(i,j,k) .le. 1.e-10 ) then
                        ttarget(i,j,k)=missing
                     else
                        ttarget(i,j,k) = ttarget(i,j,k)
     &                       /atarget(i,j,k)
                     endif
                  enddo
               enddo
            enddo
         endif                  ! am_i_root
      endif                     !ll2cs
         
      end subroutine parallel_regrid_wt
c*


      subroutine repr_regrid(x2grids,tsource,ttarget,atarget)
 
!@sum  Parallel regriding using double-double arithmetics for reproducible results
!@+    independent of #procs
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tsource(x2grids%isd:x2grids%ied,
     &     x2grids%jsd:x2grids%jed)
      real*8 :: ttarget(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
     &     ,atarget(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
      complex*16, allocatable :: ttarget_dd(:,:,:),atarget_dd(:,:,:)

      complex*16 :: ttarget_tmp(1),tmp_dd,atarget_tmp(1),xarea_dd(1),
     &     tmp_dd_vec(1)
      integer :: n,icub,jcub,i,j,k,itile,ilon,jlat,ikey,imlon,jmlat,
     &     imcub,jmcub,ntcub
      character*120:: ofi
      integer ::  status, fid, vid
     

      allocate(
     &     ttarget_dd(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget),
     &     atarget_dd(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
     &     )

      if ( (x2grids%ntilessource .eq. 6) .and. !cs2ll
     &     (x2grids%ntilestarget .eq. 1) ) then   
         
         imlon=x2grids%imtarget 
         jmlat=x2grids%jmtarget
         
         do j=1,jmlat
            do i=1,imlon
               atarget_dd(i,j,1) = dcmplx(0.d0,0.d0)
               ttarget_dd(i,j,1) = dcmplx(0.d0,0.d0)
            enddo
         enddo

c         write(*,*) "maxkey==",x2grids%xgrid%maxkey

         do ikey=1,x2grids%xgrid%maxkey
            
            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)
          
            atarget_tmp(1)=atarget_dd(ilon,jlat,1)
            xarea_dd(1)=dcmplx(x2grids%xgrid%xarea_key(ikey)
     &           , 0.d0)
            call add_dd(xarea_dd,atarget_tmp,1) 
            atarget_dd(ilon,jlat,1)=atarget_tmp(1)
            
            ttarget_tmp(1)=ttarget_dd(ilon,jlat,1)
            call dxd(x2grids%xgrid%xarea_key(ikey),
     &           tsource(icub,jcub),tmp_dd)
            tmp_dd_vec(1)=tmp_dd
            call add_dd(tmp_dd_vec,ttarget_tmp,1)   
            ttarget_dd(ilon,jlat,1)=ttarget_tmp(1)
            
         enddo
         
         call sumxpe2d_exact(ttarget_dd,imlon,jmlat,1)
         call sumxpe2d_exact(atarget_dd,imlon,jmlat,1)
         
c     
c     root proc section
c     
         
         if (AM_I_ROOT()) then   
            do j=1,jmlat
               do i=1,imlon
                  ttarget(i,j,1)=real(ttarget_dd(i,j,1))
                  atarget(i,j,1)=real(atarget_dd(i,j,1))
               enddo
            enddo
            
            ttarget(:,:,1) = ttarget(:,:,1)/atarget(:,:,1)
            
c            ofi='tstout.nc'
c            status = nf_open(trim(ofi),nf_write,fid)
c            if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
c            status = nf_inq_varid(fid,'lwup_sfc',vid)
c            write(*,*) NF_STRERROR(status)
c            status = nf_put_var_double(fid,vid,ttarget)
c            write(*,*) "STATUS",NF_STRERROR(status),"<<"
c            status = nf_close(fid)
         endif


         
      else if ( (x2grids%ntilessource .eq. 1) .and. !ll2cs
     &        (x2grids%ntilestarget .eq. 6) ) then
         
         imcub=x2grids%imtarget
         jmcub=x2grids%jmtarget
         ntcub=x2grids%ntilestarget

         do k=1,ntcub
            do j=1,jmcub
               do i=1,imcub
                  atarget_dd(i,j,k) =  dcmplx(0.d0,0.d0)
                  ttarget_dd(i,j,k) =  dcmplx(0.d0,0.d0)
               enddo
            enddo
         enddo
         
c         write(*,*) "maxkey regx=",x2grids%xgrid%maxkey
         
         do ikey=1,x2grids%xgrid%maxkey
            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)
            
            atarget_tmp(1)=atarget_dd(icub,jcub,itile)
            xarea_dd(1)=dcmplx(x2grids%xgrid%xarea_key(ikey)
     &           , 0.d0)
            call add_dd(xarea_dd,atarget_tmp,1) 
            atarget_dd(icub,jcub,itile)=atarget_tmp(1) 
            
            ttarget_tmp(1)=ttarget_dd(icub,jcub,itile)
            call dxd(x2grids%xgrid%xarea_key(ikey),
     &           tsource(ilon,jlat),tmp_dd)
            tmp_dd_vec(1)=tmp_dd
            call add_dd(tmp_dd_vec,ttarget_tmp,1)   
            ttarget_dd(icub,jcub,itile)=ttarget_tmp(1)
            
         enddo
         
         call sumxpe2d_exact(ttarget_dd,imcub,jmcub,ntcub)
         call sumxpe2d_exact(atarget_dd,imcub,jmcub,ntcub)
         
c     
c     root proc section
c     
         
         if (AM_I_ROOT()) then   
            do k=1,ntcub
               do j=1,jmcub
                  do i=1,imcub
                     ttarget(i,j,k) = real(ttarget_dd(i,j,k))
                     atarget(i,j,k) = real(atarget_dd(i,j,k))
                     ttarget(i,j,k) = ttarget(i,j,k)/atarget(i,j,k)                    
                  enddo
               enddo
            enddo         
         endif
         
      endif

      deallocate(ttarget_dd,atarget_dd)

      end subroutine repr_regrid
c*


      subroutine repr_regrid_wt(x2grids,wsource,missing,tsource,
     &     ttarget,atarget)
!@sum  same as repr_regrid but we multiply by weight function wsource
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tsource(x2grids%isd:x2grids%ied,
     &     x2grids%jsd:x2grids%jed)
      real*8 :: wsource(x2grids%isd:x2grids%ied,
     &     x2grids%jsd:x2grids%jed)
      real*8 :: ttarget(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
     &     ,atarget(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
      real*8, intent(in) :: missing
      real *8 :: tmp_r8
      complex*16, allocatable :: ttarget_dd(:,:,:),atarget_dd(:,:,:)
      complex*16 :: ttarget_tmp(1),tmp_dd,atarget_tmp(1),xarea_dd(1),
     &     tmp_dd_vec(1),xaws_dd(1)
      integer :: n,icub,jcub,i,j,k,itile,ilon,jlat,ikey,imlon,jmlat,
     &     imcub,jmcub,ntcub
      character*120:: ofi
      integer ::  status, fid, vid

      allocate(
     &     ttarget_dd(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget),
     &     atarget_dd(x2grids%imtarget,x2grids%jmtarget
     &     ,x2grids%ntilestarget)
     &     )     

c      write(*,*) "bef test direction",x2grids%ntilessource

      if ( (x2grids%ntilessource .eq. 6) .and. !cs2ll
     &     (x2grids%ntilestarget .eq. 1) ) then   
         
c         write(*,*) "inside cs2ll"

         imlon=x2grids%imtarget 
         jmlat=x2grids%jmtarget
c         write(*,*) "imlon, jmlat=",imlon,jmlat
         do j=1,jmlat
            do i=1,imlon
               atarget_dd(i,j,1) = dcmplx(0.d0,0.d0)
               ttarget_dd(i,j,1) = dcmplx(0.d0,0.d0)
            enddo
         enddo

c         write(*,*) "maxkey==",x2grids%xgrid%maxkey

         do ikey=1,x2grids%xgrid%maxkey
            
            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)
          
            atarget_tmp(1)=atarget_dd(ilon,jlat,1)
            call dxd(x2grids%xgrid%xarea_key(ikey),
     &           wsource(icub,jcub),tmp_dd)
            xaws_dd(1)=tmp_dd
            call add_dd(xaws_dd,atarget_tmp,1) 
            atarget_dd(ilon,jlat,1)=atarget_tmp(1)
            
            ttarget_tmp(1)=ttarget_dd(ilon,jlat,1)
            call dxd(x2grids%xgrid%xarea_key(ikey),
     &           tsource(icub,jcub),tmp_dd)
            tmp_r8=real(tmp_dd)
            call dxd(tmp_r8,wsource(icub,jcub),tmp_dd)
            tmp_dd_vec(1)=tmp_dd
            call add_dd(tmp_dd_vec,ttarget_tmp,1)   
            ttarget_dd(ilon,jlat,1)=ttarget_tmp(1)
            
         enddo
         
         call sumxpe2d_exact(ttarget_dd,imlon,jmlat,1)
         call sumxpe2d_exact(atarget_dd,imlon,jmlat,1)
         
c     
c     root proc section
c     
         
         if (AM_I_ROOT()) then   
            do j=1,jmlat
               do i=1,imlon
                  ttarget(i,j,1)=real(ttarget_dd(i,j,1))
                  atarget(i,j,1)=real(atarget_dd(i,j,1))
                  if (atarget(i,j,1) .le. 1.e-10) then
                     ttarget(i,j,1)=missing
c                     write(*,*) "small=",i,j,atarget(i,j,1)
                  else
                     ttarget(i,j,1) = ttarget(i,j,1)/atarget(i,j,1)
                  endif
               enddo
            enddo
            
         endif


         
      else if ( (x2grids%ntilessource .eq. 1) .and. !ll2cs
     &        (x2grids%ntilestarget .eq. 6) ) then
         
         imcub=x2grids%imtarget
         jmcub=x2grids%jmtarget
         ntcub=x2grids%ntilestarget

         do k=1,ntcub
            do j=1,jmcub
               do i=1,imcub
                  atarget_dd(i,j,k) =  dcmplx(0.d0,0.d0)
                  ttarget_dd(i,j,k) =  dcmplx(0.d0,0.d0)
               enddo
            enddo
         enddo
         
c         write(*,*) "maxkey regx=",x2grids%xgrid%maxkey
         
         do ikey=1,x2grids%xgrid%maxkey
            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)
            

            atarget_tmp(1)=atarget_dd(icub,jcub,itile)
            call dxd(x2grids%xgrid%xarea_key(ikey),
     &           wsource(ilon,jlat),tmp_dd)
            xaws_dd(1)=tmp_dd
            call add_dd(xaws_dd,atarget_tmp,1) 
            atarget_dd(icub,jcub,itile)=atarget_tmp(1)
            
            ttarget_tmp(1)=ttarget_dd(icub,jcub,itile)
            call dxd(x2grids%xgrid%xarea_key(ikey),
     &           tsource(ilon,jlat),tmp_dd)
            tmp_r8=real(tmp_dd)
            call dxd(tmp_r8,wsource(ilon,jlat),tmp_dd)
            tmp_dd_vec(1)=tmp_dd
            call add_dd(tmp_dd_vec,ttarget_tmp,1)   
            ttarget_dd(icub,jcub,itile)=ttarget_tmp(1)
            
         enddo
         
         call sumxpe2d_exact(ttarget_dd,imcub,jmcub,ntcub)
         call sumxpe2d_exact(atarget_dd,imcub,jmcub,ntcub)
         
c     
c     root proc section
c     
         
         if (AM_I_ROOT()) then   
            do k=1,ntcub
               do j=1,jmcub
                  do i=1,imcub
                     ttarget(i,j,k) = real(ttarget_dd(i,j,k))
                     atarget(i,j,k) = real(atarget_dd(i,j,k))
                     if (atarget(i,j,k) .le. 1.e-10 ) then
                        ttarget(i,j,k)=missing
                     else
                        ttarget(i,j,k) = ttarget(i,j,k)
     &                       /atarget(i,j,k)
                     endif
                  enddo
               enddo
            enddo         
         endif
         
      endif

      deallocate(ttarget_dd,atarget_dd)  

      end subroutine repr_regrid_wt
c*

      
      
      subroutine root_regrid(x2gridsroot,tsource,ttarget)

!@sum  Root processor regrids data from source grid to target grid
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2gridsroot
      real*8 :: tsource(x2gridsroot%imsource,x2gridsroot%jmsource,
     &     x2gridsroot%ntilessource)
      real*8 :: ttarget(x2gridsroot%imtarget,x2gridsroot%jmtarget,
     &     x2gridsroot%ntilestarget)
     &     ,atarget(x2gridsroot%imtarget,x2gridsroot%jmtarget,
     &     x2gridsroot%ntilestarget)
      integer :: n,icub,jcub,i,j,itile,icc,jcc,il,jl


      if ((x2gridsroot%ntilessource .eq. 1) .and. 
     &     (x2gridsroot%ntilestarget .eq. 6)) then
         
c     ll2cs
         
         if (AM_I_ROOT()) then   
c            write(*,*) "ROOT REGRID LL2CS"
 
            atarget(:,:,:) = 0.d0
            ttarget(:,:,:) = 0.d0
            

            do n=1,x2gridsroot%xgridroot%ncells
               itile=x2gridsroot%xgridroot%tile(n)
               icc=x2gridsroot%xgridroot%ijcub(1,n)
               jcc=x2gridsroot%xgridroot%ijcub(2,n)
               il=x2gridsroot%xgridroot%ijlatlon(1,n)
               jl=x2gridsroot%xgridroot%ijlatlon(2,n)

               atarget(icc,jcc,itile) = atarget(icc,jcc,itile) 
     &              + x2gridsroot%xgridroot%xgrid_area(n)
               ttarget(icc,jcc,itile) = ttarget(icc,jcc,itile) 
     &              + x2gridsroot%xgridroot%xgrid_area(n)
     &              *tsource(il,jl,1)
            enddo
            
            do itile=1,x2gridsroot%ntilestarget
               do j=1,x2gridsroot%jmtarget
                  do i=1,x2gridsroot%imtarget
                     ttarget(i,j,itile) = ttarget(i,j,itile)
     &                    /atarget(i,j,itile)
                  enddo
               enddo
            enddo
c            write(6,*) "END ROOT REGRID LL2CS"
         endif
         
      endif
      
      end subroutine root_regrid
c*

      subroutine root_regrid_wt(x2gridsroot,wsource,missing,
     &     tsource,ttarget)

!@sum  Root processor regrids data from source grid to target grid
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2gridsroot
      real*8 :: tsource(x2gridsroot%imsource,x2gridsroot%jmsource,
     &     x2gridsroot%ntilessource)
      real*8 :: ttarget(x2gridsroot%imtarget,x2gridsroot%jmtarget,
     &     x2gridsroot%ntilestarget)
     &     ,atarget(x2gridsroot%imtarget,x2gridsroot%jmtarget,
     &     x2gridsroot%ntilestarget)
      real*8 :: wsource(x2gridsroot%imsource,x2gridsroot%jmsource,
     &     x2gridsroot%ntilessource)
      real *8 :: missing   ! missing value 
      integer :: n,icub,jcub,i,j,itile,icc,jcc,il,jl


      if ((x2gridsroot%ntilessource .eq. 1) .and. 
     &     (x2gridsroot%ntilestarget .eq. 6)) then
         
c     ll2cs
         
         if (AM_I_ROOT()) then   
 
            atarget(:,:,:) = 0.d0
            ttarget(:,:,:) = 0.d0
            

            do n=1,x2gridsroot%xgridroot%ncells
               itile=x2gridsroot%xgridroot%tile(n)
               icc=x2gridsroot%xgridroot%ijcub(1,n)
               jcc=x2gridsroot%xgridroot%ijcub(2,n)
               il=x2gridsroot%xgridroot%ijlatlon(1,n)
               jl=x2gridsroot%xgridroot%ijlatlon(2,n)

               if ( tsource(il,jl,1) .gt. missing ) then
               atarget(icc,jcc,itile) = atarget(icc,jcc,itile) 
     &              + wsource(il,jl,1)
     &              *x2gridsroot%xgridroot%xgrid_area(n)

               ttarget(icc,jcc,itile) = ttarget(icc,jcc,itile) 
     &              + x2gridsroot%xgridroot%xgrid_area(n)
     &              *wsource(il,jl,1)
     &              *tsource(il,jl,1)
               endif
            enddo
            
            do itile=1,x2gridsroot%ntilestarget
               do j=1,x2gridsroot%jmtarget
                  do i=1,x2gridsroot%imtarget
                     if (atarget(i,j,itile) .le. 0.) then
                        ttarget(i,j,itile) = missing
                     else
                        ttarget(i,j,itile) = ttarget(i,j,itile)
     &                    /atarget(i,j,itile)
                     endif
                  enddo
               enddo
            enddo
c            write(6,*) "END ROOT REGRID LL2CS"
         endif
         
      endif
      
      end subroutine root_regrid_wt
c*

         
      subroutine init_regrid(x2grids,grid,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget)

!@sum  Reads regriding file on root proc, broadcasts the 
!@+    x_grid data to all processes then instanciates locally
!@+    the x_2grids derived type (x_2grids type=x_grid plus info about
!@+    source and target grids). It also initializes domain decomposition
!@+    variables through dist_grid derived type. 
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'

      type (x_2grids), intent(inout) :: x2grids
      type (x_grid) :: xgrid
      type (dist_grid), intent(in) :: grid

      real*8,  allocatable, dimension(:) :: xgrid_area  !local variable
      integer, allocatable, dimension(:) :: tile        !local variable
      integer, allocatable, dimension(:,:) :: ijcub     !local variable
      integer, allocatable, dimension(:,:) :: ijlatlon  !local variable
      integer :: ncells                                 !local variable
      integer, intent(in) :: imsource,jmsource,imtarget,jmtarget
      integer, intent(in) :: ntilessource,ntilestarget
      integer :: status,fid,n,vid,ikey,jlat
      integer :: itile,j,idomain,iic,jjc,index,indexc,nc2
      integer :: ierr,i,icub,jcub,maxkey,gid,mytile
      character(len=30) :: fname

c     domain decomp. info on source grid is stored inside x2grid object
      x2grids%is=grid%is
      x2grids%ie=grid%ie
      x2grids%isd=grid%isd
      x2grids%ied=grid%ied

      x2grids%js=grid%js
      x2grids%je=grid%je
      x2grids%jsd=grid%jsd
      x2grids%jed=grid%jed

      x2grids%imsource=imsource
      x2grids%jmsource=jmsource
      x2grids%ntilessource=ntilessource
      x2grids%imtarget=imtarget
      x2grids%jmtarget=jmtarget
      x2grids%ntilestarget=ntilestarget

c     local variables
      gid=grid%gid
      mytile=grid%tile
c
 
      if (AM_I_ROOT()) then   

c     
c     Read weights
c     
         fname='REMAP'
         write(*,*) "REMAP file=",fname
         status = nf_open(trim(fname),nf_nowrite,fid)
         if (status .ne. NF_NOERR) write(*,*) 
     *        "UNABLE TO OPEN REMAP FILE"
         
         status = nf_inq_dimid(fid,'ncells',vid)
         status = nf_inq_dimlen(fid,vid,ncells)
      endif
            
c     Broadcast value of ncells & Allocate arrays with size 
c     depending on ncells on each processor

      call MPI_BCAST( ncells, 1, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))


      if (AM_I_ROOT()) then   
         status = nf_inq_varid(fid,'xgrid_area',vid)
         status = nf_get_var_double(fid,vid,xgrid_area)
         status = nf_inq_varid(fid,'tile1',vid)
         status = nf_get_var_int(fid,vid,tile)
         status = nf_inq_varid(fid,'tile1_cell',vid)
         status = nf_get_var_int(fid,vid,ijcub)
         status = nf_inq_varid(fid,'tile2_cell',vid)
         status = nf_get_var_int(fid,vid,ijlatlon)
         status = nf_close(fid)
      endif
      
c
c     Broadcast x_grid area and indices to all procs
c
      nc2=2*ncells
     
      call MPI_BCAST( xgrid_area, ncells, MPI_DOUBLE_PRECISION,
     *     0, MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( ijcub, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( ijlatlon, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( tile, ncells, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 
      
      
      if ((ntilessource .eq. 6) .and. (ntilestarget .eq. 1)) then   !cs2ll
         
c     fill x2grids with values local to the domain
c     first calculate maxkey
         
         ikey=1
         do n=1,ncells
            icub=ijcub(1,n)
            jcub=ijcub(2,n)
            itile=tile(n)
            if (itile .eq. mytile) then
               if ( ( (icub .le. x2grids%ie) .and. 
     &              (icub .ge. x2grids%is)) .and.
     &              ( (jcub .le. x2grids%je) .and. 
     &              (jcub .ge. x2grids%js)) ) then
                  ikey=ikey+1
               endif
            endif
         enddo
         
         maxkey=ikey-1
         
c         write(*,*) "gid maxkey=",gid,maxkey

         allocate(x2grids%xgrid%icub_key(maxkey),
     &        x2grids%xgrid%jcub_key(maxkey),
     &        x2grids%xgrid%ilon_key(maxkey),
     &        x2grids%xgrid%jlat_key(maxkey),
     &        x2grids%xgrid%xarea_key(maxkey),
     &        x2grids%xgrid%itile_key(maxkey)
     &        ) 
         
         ikey=1
         
         do n=1,ncells
            icub=ijcub(1,n)
            jcub=ijcub(2,n)
            itile=tile(n)
            
            if (itile .eq. mytile) then
               if ( ( (icub .le. x2grids%ie) .and. 
     &              (icub .ge. x2grids%is)) .and.
     &              ( (jcub .le. x2grids%je) .and. 
     &              (jcub .ge. x2grids%js)) ) then
                  x2grids%xgrid%icub_key(ikey)=icub
                  x2grids%xgrid%jcub_key(ikey)=jcub
                  x2grids%xgrid%ilon_key(ikey)=ijlatlon(1,n)
                  x2grids%xgrid%jlat_key(ikey)=ijlatlon(2,n)
                  x2grids%xgrid%xarea_key(ikey)=xgrid_area(n)
                  x2grids%xgrid%itile_key(ikey)=tile(n)
                  ikey=ikey+1
               endif
            endif
         enddo

         x2grids%xgrid%maxkey=maxkey
         
      else if ((ntilessource .eq. 1) .and. (ntilestarget .eq. 6)) then !ll2cs
         
         ikey=1
         do n=1,ncells
            jlat=ijlatlon(2,n)
            if ( (jlat .le. x2grids%je) .and. 
     &           (jlat .ge. x2grids%js) ) then !js = J_0, je=J_1
               ikey=ikey+1
            endif
         enddo
         
         maxkey=ikey-1
         
c         write(*,*) "maxkey=",maxkey
         
         allocate(x2grids%xgrid%icub_key(maxkey),
     &        x2grids%xgrid%jcub_key(maxkey),
     &        x2grids%xgrid%ilon_key(maxkey),
     &        x2grids%xgrid%jlat_key(maxkey),
     &        x2grids%xgrid%xarea_key(maxkey),
     &        x2grids%xgrid%itile_key(maxkey)
     &        ) 
         
         ikey=1
         do n=1,ncells
            icub=ijcub(1,n)
            jcub=ijcub(2,n)
            jlat=ijlatlon(2,n)
c            write(*,*) "jlat=",jlat
            if ( (jlat .le. x2grids%je) .and. 
     &           (jlat .ge. x2grids%js) ) then !js = J_0, je=J_1
               x2grids%xgrid%icub_key(ikey)=icub
               x2grids%xgrid%jcub_key(ikey)=jcub
               x2grids%xgrid%ilon_key(ikey)=ijlatlon(1,n)
               x2grids%xgrid%jlat_key(ikey)=ijlatlon(2,n)
               x2grids%xgrid%xarea_key(ikey)=xgrid_area(n)
               x2grids%xgrid%itile_key(ikey)=tile(n)
               ikey=ikey+1
            endif
         enddo

         x2grids%xgrid%maxkey=maxkey
         
c         write(*,*)  x2grids%xgrid%icub_key(1:ikey)
      endif    
      
      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)
      
      end subroutine init_regrid
c*


      subroutine dxd (da, db, ddc)
c  This subroutine computes ddc = da x db.
      use regrid_com
      implicit none
      real*8 da, db
      real*8 a1, a2, b1, b2, con, split, t1, t2
      complex*16 ddc
      parameter (split = 134217729.d0)

c   Split da and db into two parts with at most 26 bits each, using the 
c   Dekker-Veltkamp method.
      con = da * split
      a1 = con - (con - da)
      a2 = da - a1
      con = db * split
      b1 = con - (con - db)
      b2 = db - b1

c   Compute da * db using Dekker's method.
      t1 = da * db
      t2 = (((a1*b1 - t1) + a1*b2) + a2*b1) + a2*b2

      ddc = dcmplx(t1, t2)

      return
      end subroutine dxd
c*

      subroutine add_DD (dda, ddb, len)
!@sum  Compute dda + ddb using double-double arithmetics
!@auth Denis Gueyffier, adapted from  Yun He and Chris Ding, 
!@auth Journal of Supercomputing, 2001
      use regrid_com
      implicit none
      real*8 e, t1, t2
      integer i, len
      complex*16 :: dda(len), ddb(len)

      do i = 1, len
c   Compute dda + ddb using Knuth's trick.
      t1 = real(dda(i)) + real(ddb(i))
      e = t1 - real(dda(i))
      t2 = ((real(ddb(i)) - e) + (real(dda(i)) - (t1 - e)))
     &     +imag(dda(i)) + imag(ddb(i))

c   The result is t1 + t2, after normalization.
      ddb(i) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1) )
      enddo

      return
      end subroutine add_DD
c*

      subroutine test_globalsum_exact()
!@sum Testing reproducible globalsum 
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'mpif.h'
      integer i, imt, jmt
      integer :: myPE, totPEs, stat(MPI_STATUS_SIZE), ierr
      integer :: start, end, len, MPI_SUMDD, itype
      external add_DD 
      real*8, dimension(:), allocatable:: array, local_array
      complex*16 :: global_sum

      imt = 120
      jmt = 64

      call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, totPEs, ierr )

      len = imt*jmt/totPEs + 1
c      write(*,*) 'myPE=', myPE, ' len=', len
      allocate (local_array(len))

C  read data for root proc 
      if (myPE .eq. 0) then
         allocate(array(len*totPEs))
         open(10, file ='etaana.dat', status='unknown')
         do i= 1, imt*jmt
            read(10,*) array(i)
         enddo
         close(10)
         do i=imt*jmt+1,len*totPEs
            array(i)=0.0d0
         enddo
      endif 

C  scatter array
      call MPI_SCATTER(array, len, MPI_REAL8,
     &     local_array, len, MPI_REAL8, 0,
     &     MPI_COMM_WORLD, ierr)

      if (myPE .eq. 0) deallocate (array) 

      call globalsum_exact(global_sum,local_array,len)

      deallocate (local_array)

      end subroutine test_globalsum_exact
c*     



      subroutine globalsum_exact(global_sum,local_array,len)
!@sum This code calculates a reproducible globalsum 
!@+   independent on number of processors using double-double precision method
!@+   complex->complex*16, cmplx->dcmplx, 
!@+   MPI_REAL->MPI_REAL8, MPI_COMPLEX->MPI_DOUBLE_COMPLEX
!@auth Denis Gueyffier, based on Yun He and Chris Ding, 
!@+   Journal of Supercomputing, 2001

      use regrid_com
      implicit none
      include 'mpif.h'
      integer i, imt, jmt
      integer :: myPE, stat(MPI_STATUS_SIZE), ierr
      integer ::  len, MPI_SUMDD, itype
      external add_DD
      real*8, dimension(len):: local_array
      complex*16 :: local_sum, global_sum

c      write(*,*) "BEGIN gsum"

      call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )

C  operator MPI_SUMDD is created based on an external function add_dd
      call MPI_OP_CREATE(add_DD, .TRUE., MPI_SUMDD, ierr)

      
C  Each processor calculates the local_sum of its own section first. 
C  Complex number is defined to represent local_sum and local_sum err.
      local_sum = 0.0
      do i = 1, len
         call add_DD(dcmplx(local_array(i), 0.0d0), local_sum,1)
      enddo

c      write(*,*) "MyPE bef reduce=",myPE

C  add all local_sums on each PE to PE0 with MPI_SUMDD.
C  global_sum is a complex number, represents final (sum, error).
      call MPI_REDUCE (local_sum, global_sum, 1, MPI_DOUBLE_COMPLEX, 
     &     MPI_SUMDD, 0, MPI_COMM_WORLD, ierr)

c      write(*,*) "MyPE aft reduce=",myPE

      if (myPE.eq.0) then
         write(*,*)'quad precision sum, error= ', global_sum
      endif

      end subroutine globalsum_exact
c*

      subroutine sumxpe2d_exact(arr,iarr,jarr,ntile)
!@sum This code calculates a reproducible non rank reducing sum 
!@+   independent on number of processors using double-double arithmetics
!@auth Denis Gueyffier
      implicit none
      include 'mpif.h'
      integer :: ierr,arr_size,myPE
      integer :: MPI_SUMDD
      external add_DD
      integer, intent(in) :: iarr,jarr,ntile
      complex*16, dimension(iarr,jarr,ntile) :: arr
      complex*16, dimension(:),allocatable ::arr_rsh,arr_tmp

      call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
c     operator MPI_SUMDD is created based on an external function add_dd
      call MPI_OP_CREATE(add_DD, .TRUE., MPI_SUMDD, ierr)
      
c     reduction 
      arr_size = size(arr)

      allocate(arr_tmp(arr_size),arr_rsh(arr_size))
      arr_rsh=reshape(arr,(/arr_size/))

      call MPI_REDUCE (arr_rsh, arr_tmp, arr_size, MPI_DOUBLE_COMPLEX, 
     &     MPI_SUMDD, 0, MPI_COMM_WORLD, ierr)
      arr=reshape(arr_tmp,shape(arr))
      deallocate(arr_tmp,arr_rsh) 
    
      end subroutine sumxpe2d_exact
c*

c      subroutine init_csgrid_debug(dd2d,imc,jmc)
c
c!@sum  temporary instanciation of CS grid and domain decomp. using direct
c!@+    access to MPP. Used only for debugging purpose.
c!@auth Denis Gueyffier
c
c      use dd2d_utils
c      use mpp_mod
c      use mpp_domains_mod
c      use fv_mp_mod, only : mp_start,mp_stop
c     &     ,fv_domain_decomp=>domain_decomp
c      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
c      use fv_mp_mod, only : gid, domain, tile, npes_x, npes_y
c      use fv_grid_utils_mod, only : cosa_s,sina_s
c     &     ,grid_utils_init
cc     &     ,sw_corner,se_corner,ne_corner,nw_corner
c      use fv_arrays_mod, only: fv_atmos_type
c      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
c     &     dx, dy, dxa, dya, dxc, dyc, grid_type, dx_const, dy_const
c      use fv_control_mod, only : uniform_ppm,c2l_ord
c      
c      implicit none
c      integer :: i,j,npz
c      integer :: npes,ndims
c      integer :: commID
c      integer, dimension(:), allocatable :: pelist
c      integer :: ng, imc,jmc
c      real*8, dimension(:,:), allocatable :: P
c
c      type (dist_grid), intent(inout) :: dd2d
c      type(fv_atmos_type) :: atm
c      character*80 :: grid_name = 'Gnomonic'
c      character*120:: grid_file = 'Inline'
c      logical :: non_ortho
c
cc***  Temporarily use instanciation of grid through fv_grid_tools_mod's init_grid()
c
c      call mpp_init(MPP_VERBOSE)
cc code copied from fv_init:
c      npes = mpp_npes()
c      allocate(pelist(npes))
c      call mpp_get_current_pelist( pelist, commID=commID )
c      call mp_start(commID)
cc      write(6,*) 'commID ',commID
c      
c      ng = 3 ! number of ghost zones required
c
c      call fv_domain_decomp(imc+1,jmc+1,6,ng,grid_type)
cc      call mp_stop()
c
c      ndims = 2
c      npz = 5
c      call init_grid(atm,grid_name,grid_file,
c     &     imc+1, jmc+1, npz, ndims, 6, ng)
c
c      non_ortho=.true.
c      call grid_utils_init(Atm, imc+1, jmc+1, npz, Atm%grid, Atm%agrid,
c     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
c     &     uniform_ppm, grid_type, c2l_ord)
c
cc**   transfering domain decomposition info to Max's dd2d derived type
c      call init_dist_grid(
c     &     imc,jmc,6, 
c     &     is,ie,js,je,isd,ied,jsd,jed,dd2d)
c      
c      write(*,*) "is,ie,js,je,isd,ied,jsd,jed==",
c     &     is,ie,js,je,isd,ied,jsd,jed
c
c      end subroutine init_csgrid_debug
cc****



