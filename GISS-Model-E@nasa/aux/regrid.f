      module regrid
!@auth D,Gueyffier
!@ver 1.0
      
      private
      public :: x_2grids, init_regrid, do_regrid, do_regrid_wt, 
     &     nrecmax

      integer, parameter :: nrecmax=200 ! max #records in input files

      type x_grid      ! stores exchange grid information 
      integer  :: ncells
      integer, allocatable, dimension(:,:) :: ijcub
      integer, allocatable, dimension(:,:) :: ijlatlon
      real*8, allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:) :: tile
      integer :: maxkey      
      end type x_grid

      type x_2grids   ! stores x-grid info plus source and target grid info 
      type (x_grid) :: xgrid
      integer :: ntilessource  ! #tiles of source grid (1 for latlon, 6 for cubed sphere)
      integer :: ntilestarget  ! #tiles of target grid (1 for latlon, 6 for cubed sphere)
      integer :: imsource      ! im for source grid
      integer :: jmsource      ! jm for source grid
      integer :: imtarget      ! im for target grid
      integer :: jmtarget      ! jm for target grid
      end type x_2grids

      contains 

      subroutine init_regrid(x2grids,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget,fname)

!@sum  Reads regriding file then instanciates the x_2grids derived type (x_2grids type=x_grid plus info about
!@+    source and target grids). 
!@auth Denis Gueyffier
      implicit none
      include 'netcdf.inc'

      type (x_2grids) :: x2grids
      type (x_grid) :: xgrid
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
      character(len=10) :: imch,jmch,icch,jcch
      character(len=30), optional :: fname
      character(len=30) :: fnamefinal

      x2grids%imsource=imsource
      x2grids%jmsource=jmsource
      x2grids%ntilessource=ntilessource
      x2grids%imtarget=imtarget
      x2grids%jmtarget=jmtarget
      x2grids%ntilestarget=ntilestarget



      if (.not.present(fname)) then
         write(imch,'(i10)') imsource
         write(jmch,'(i10)') jmsource
         write(icch,'(i10)') imtarget
         write(jcch,'(i10)') jmtarget
         imch=trim(adjustl(imch))
         jmch=trim(adjustl(jmch))
         icch=trim(adjustl(icch))
         jcch=trim(adjustl(jcch))
         fnamefinal="remap"//trim(imch)//"-"//trim(jmch)
     *        //"C"//trim(icch)//"-"//trim(jcch)//".nc"
         write(*,*) "remap file=",fnamefinal
      else
         fnamefinal=fname
      endif
c     
c     Read weights
c     
      status = nf_open(trim(fnamefinal),nf_nowrite,fid)
      
      if (status .ne. NF_NOERR) write(*,*) 
     *     "UNABLE TO OPEN REMAP FILE",trim(fnamefinal)
      
      status = nf_inq_dimid(fid,'ncells',vid)
      status = nf_inq_dimlen(fid,vid,ncells)
     

      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))

      status = nf_inq_varid(fid,'xgrid_area',vid)
      status = nf_get_var_double(fid,vid,xgrid_area)
      status = nf_inq_varid(fid,'tile1',vid)
      status = nf_get_var_int(fid,vid,tile)
      status = nf_inq_varid(fid,'tile1_cell',vid)
      status = nf_get_var_int(fid,vid,ijcub)
      status = nf_inq_varid(fid,'tile2_cell',vid)
      status = nf_get_var_int(fid,vid,ijlatlon)
      status = nf_close(fid)
    
      allocate(x2grids%xgrid%xgrid_area(ncells))
      allocate(x2grids%xgrid%ijcub(2,ncells))
      allocate(x2grids%xgrid%ijlatlon(2,ncells))
      allocate(x2grids%xgrid%tile(ncells))
      
      x2grids%xgrid%xgrid_area(:)=xgrid_area(:)
      x2grids%xgrid%ijcub(:,:)=ijcub(:,:)
      x2grids%xgrid%ijlatlon(:,:)=ijlatlon(:,:)
      x2grids%xgrid%tile(:)=tile(:)
      x2grids%xgrid%ncells=ncells

      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      write(*,*) "initialized regridding weights"

      end subroutine init_regrid

      subroutine do_regrid(x2grids,tsource,ttarget)
!@sum  regrid data from source grid to target grid
!@auth Denis Gueyffier
      implicit none
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tsource(x2grids%imsource,x2grids%jmsource,
     &     x2grids%ntilessource)
      real*8 :: ttarget(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget)
     &     ,atarget(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget)
      integer :: n,icub,jcub,i,j,itile,icc,jcc,il,jl

      if ((x2grids%ntilessource .eq. 1) .and.
     &     (x2grids%ntilestarget .eq. 6)) then

c     ll2cs

         atarget(:,:,:) = 0.d0
         ttarget(:,:,:) = 0.d0


         do n=1,x2grids%xgrid%ncells
            itile=x2grids%xgrid%tile(n)
            icc=x2grids%xgrid%ijcub(1,n)
            jcc=x2grids%xgrid%ijcub(2,n)
            il=x2grids%xgrid%ijlatlon(1,n)
            jl=x2grids%xgrid%ijlatlon(2,n)

            atarget(icc,jcc,itile) = atarget(icc,jcc,itile)
     &           + x2grids%xgrid%xgrid_area(n)
            ttarget(icc,jcc,itile) = ttarget(icc,jcc,itile)
     &           + x2grids%xgrid%xgrid_area(n)
     &           *tsource(il,jl,1)
         enddo

         do itile=1,x2grids%ntilestarget
            do j=1,x2grids%jmtarget
               do i=1,x2grids%imtarget
                  ttarget(i,j,itile) = ttarget(i,j,itile)
     &                 /atarget(i,j,itile)
               enddo
            enddo
         enddo

      endif

      end subroutine do_regrid

      subroutine do_regrid_wt(x2grids,wsource,missing,
     &     tsource,ttarget)
!@sum  regrid data from source grid to target grid
!@auth Denis Gueyffier
      implicit none
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tsource(x2grids%imsource,x2grids%jmsource,
     &     x2grids%ntilessource)
      real*8 :: wsource(x2grids%imsource,x2grids%jmsource,
     &     x2grids%ntilessource)
      real*8 :: ttarget(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget)
     &     ,atarget(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget)
      real*8, intent(in) :: missing
      integer :: n,icub,jcub,i,j,itile,icc,jcc,il,jl



      if ((x2grids%ntilessource .eq. 1) .and.
     &     (x2grids%ntilestarget .eq. 6)) then

c     ll2cs

         atarget(:,:,:) = 0.d0
         ttarget(:,:,:) = 0.d0


         do n=1,x2grids%xgrid%ncells
            itile=x2grids%xgrid%tile(n)
            icc=x2grids%xgrid%ijcub(1,n)
            jcc=x2grids%xgrid%ijcub(2,n)
            il=x2grids%xgrid%ijlatlon(1,n)
            jl=x2grids%xgrid%ijlatlon(2,n)

            if ( tsource(il,jl,1) .gt. missing ) then
               atarget(icc,jcc,itile) = atarget(icc,jcc,itile)
     &              + wsource(il,jl,1)
     &              *x2grids%xgrid%xgrid_area(n)
               
               ttarget(icc,jcc,itile) = ttarget(icc,jcc,itile)
     &              + x2grids%xgrid%xgrid_area(n)
     &              *wsource(il,jl,1)
     &              *tsource(il,jl,1)
            endif
         enddo

         do itile=1,x2grids%ntilestarget
            do j=1,x2grids%jmtarget
               do i=1,x2grids%imtarget
                  ttarget(i,j,itile) = ttarget(i,j,itile)
     &                 /atarget(i,j,itile)
               enddo
            enddo
         enddo
      endif

      end subroutine do_regrid_wt

      end module regrid
