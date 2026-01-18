      module REGRID_COM
!@sum Contains conservative regridding, interpolation and zonal mean routines
!@+   to exchange data between grids (latlon and cubed sphere) 
!@+ 
!@+    We define the x_grid container which stores information and data about the 
!@+    exchange grid. This data is common to the source and target grids
!@+    see http://www.gfdl.noaa.gov/~vb/gridstd/gridstdse2.html#x4-160002.7
!@     for precise definition of the concept of exchange grid
!@+    We also define the x_2grids type which contains an x_grid (exchange grid)
!@+    plus info about the direction of the remapping and about the resolution 
!@+    of the source and target grids
!@+ 
!@auth Denis Gueyffier


      use DOMAIN_DECOMP_1D, only: AM_I_ROOT,SUMXPE
      use dd2d_utils, only : dist_grid
!?      use cs2ll_utils, only : xgridremap_type

      integer, parameter :: nrecmax=200 ! max #records in input files

      private :: xgrid
      public :: x_2grids,x_2gridsroot,init_regrid_root

ccc   derived types
      type x_grid   ! stores x-grid information common to source and target grids
                    ! x_grid is distributed and contains data local to each domain
      integer, allocatable, dimension(:) :: index_key
      integer, allocatable, dimension(:) :: icub_key
      integer, allocatable, dimension(:) :: jcub_key
      integer, allocatable, dimension(:) :: ilon_key
      integer, allocatable, dimension(:) :: jlat_key
      real*8, allocatable, dimension(:) :: xarea_key
      integer, allocatable, dimension(:) :: itile_key
      integer :: maxkey      
      integer :: ncells
      end type x_grid

      type x_gridroot   ! stores x-grid information common to source and target grids
                        ! x_grid data is global and stays on root proc
      integer  :: ncells
      integer, allocatable, dimension(:,:) :: ijcub
      integer, allocatable, dimension(:,:) :: ijlatlon
      real*8, allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:) :: tile
      integer :: maxkey      
      end type x_gridroot

      type x_2grids   ! stores x-grid info plus source and target grid info 
      type (x_grid) :: xgrid
      integer :: ntilessource  ! #tiles of source grid (1 for latlon, 6 for cubed sphere)
      integer :: ntilestarget  ! #tiles of target grid (1 for latlon, 6 for cubed sphere)
      integer :: imsource      ! im for source grid
      integer :: jmsource      ! jm for source grid
      integer :: imtarget      ! im for target grid
      integer :: jmtarget      ! jm for target grid
      integer :: is,ie,js,je,isd,ied,jsd,jed ! domain decomp. for source grid
      end type x_2grids

      type x_2gridsroot   ! stores x-grid info plus source and target grid info 
      type (x_gridroot) :: xgridroot
      integer :: ntilessource  ! #tiles of source grid (1 for latlon, 6 for cubed sphere)
      integer :: ntilestarget  ! #tiles of target grid (1 for latlon, 6 for cubed sphere)
      integer :: imsource      ! im for source grid
      integer :: jmsource      ! jm for source grid
      integer :: imtarget      ! im for target grid
      integer :: jmtarget      ! jm for target grid
      end type x_2gridsroot

c***  create specializations of above types, may be moved to a different layer later
      type (x_2grids) :: xO2A    ! ocean to atmosphere
      type (x_2grids) :: xA2O    ! atmosphere to ocean
      type (x_2gridsroot) :: xA2O_root    ! atmosphere to ocean
      type (x_2gridsroot) :: xO2A_root    ! atmosphere to ocean

      contains

      subroutine init_regrid_root(x2gridsroot,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget,lremap)

!@sum  Reads regriding file on root proc, broadcasts the 
!@+    x_grid data to all processes then instanciates locally
!@+    the x_2grids derived type (x_2grids type=x_grid plus info about
!@+    source and target grids). It also initializes domain decomposition
!@+    variables through dist_grid derived type.  
!@auth Denis Gueyffier
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'
      type (x_2gridsroot), intent(inout) :: x2gridsroot
      type (x_gridroot) :: xgridroot
      real*8,  allocatable, dimension(:) :: xgrid_area  !local variable
      integer, allocatable, dimension(:) :: tile        !local variable
      integer, allocatable, dimension(:,:) :: ijcub     !local variable
      integer, allocatable, dimension(:,:) :: ijlatlon  !local variable
      integer :: ncells                                 !local variable
      integer, intent(in) :: imsource,jmsource,imtarget,jmtarget
      integer, intent(in) :: ntilessource,ntilestarget
      integer :: status,fid,n,vid,ikey,jlat
      integer :: itile,j,idomain,iic,jjc,index,indexc,nc2
      integer :: ierr,i
      logical, optional, intent(in) :: lremap
      logical :: lrem
      character(len=10) :: imch,jmch,icch,jcch	
      character(len=30) :: fname

      x2gridsroot%imsource=imsource
      x2gridsroot%jmsource=jmsource
      x2gridsroot%ntilessource=ntilessource
      x2gridsroot%imtarget=imtarget
      x2gridsroot%jmtarget=jmtarget
      x2gridsroot%ntilestarget=ntilestarget

      lrem = .false.
      if (present(lremap)) then
         if (lremap) lrem= .true.
      endif
c
      if (AM_I_ROOT()) then   

c     
c     Read weights
c     
         if (lrem) then
            write(imch,'(i10)') imsource	 
            write(jmch,'(i10)') jmsource	 
            write(icch,'(i10)') imtarget	 
            write(jcch,'(i10)') jmtarget	 
            imch=trim(adjustl(imch))	 
            jmch=trim(adjustl(jmch))	 
            icch=trim(adjustl(icch))	 
            jcch=trim(adjustl(jcch))	 
            fname="remap"//trim(imch)//"-"//trim(jmch)	 
     *           //"C"//trim(icch)//"-"//trim(jcch)//".nc"
         else
            fname='REMAP'
         endif
         write(*,*) "remap file=",fname,lrem
         status = nf_open(trim(fname),nf_nowrite,fid)

         if (status .ne. NF_NOERR) then
           write(*,*) "UNABLE TO OPEN REMAP FILE",trim(fname)
           call stop_model("UNABLE TO OPEN REMAP FILE",255)
         endif
         
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

      allocate(x2gridsroot%xgridroot%xgrid_area(ncells))
      allocate(x2gridsroot%xgridroot%ijcub(2,ncells))
      allocate(x2gridsroot%xgridroot%ijlatlon(2,ncells))
      allocate(x2gridsroot%xgridroot%tile(ncells))

      x2gridsroot%xgridroot%xgrid_area(:)=xgrid_area(:)
      x2gridsroot%xgridroot%ijcub(:,:)=ijcub(:,:)
      x2gridsroot%xgridroot%ijlatlon(:,:)=ijlatlon(:,:)
      x2gridsroot%xgridroot%tile(:)=tile(:)
      x2gridsroot%xgridroot%ncells=ncells

      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      end subroutine init_regrid_root

      end module REGRID_COM
