      module ncio
!@sum  from Max Kelley's pario_nc library
!@+    containing serial versions of defvar and write_data
!@auth Max Kelley, minor modifications by Denis Gueyffier 

      implicit none
      save
      private

      include 'netcdf.inc'

c
c i/o interfaces
c
      public :: write_data
      interface write_data
        module procedure write_nc_0D
        module procedure write_nc_1D
        module procedure write_nc_2D
        module procedure write_nc_3D
        module procedure write_nc_4D
        module procedure write_nc_5D
        module procedure write_nc_0D_r8
        module procedure write_nc_1D_r8
        module procedure write_nc_2D_r8
        module procedure write_nc_3D_r8
        module procedure write_nc_4D_r8
        module procedure write_nc_5D_r8
        module procedure write_nc_0D_int
        module procedure write_nc_1D_int
        module procedure write_nc_2D_int
        module procedure write_nc_3D_int
        module procedure write_nc_2D_logical
        module procedure write_nc_1D_array_of_strings
      end interface write_data

      public :: defvar
      interface defvar
        module procedure defvar_0D
        module procedure defvar_1D
        module procedure defvar_2D
        module procedure defvar_3D
        module procedure defvar_4D
        module procedure defvar_5D
        module procedure defvar_0D_r8
        module procedure defvar_1D_r8
        module procedure defvar_2D_r8
        module procedure defvar_3D_r8
        module procedure defvar_4D_r8
        module procedure defvar_5D_r8
        module procedure defvar_0D_int
        module procedure defvar_1D_int
        module procedure defvar_2D_int
        module procedure defvar_3D_int
        module procedure defvar_4D_int
        module procedure defvar_5D_int
        module procedure defvar_2D_logical
        module procedure defvar_1D_array_of_strings
      end interface

      contains

      subroutine write_nc_0D(fid,varname,arr)
      real*4 :: arr
#include "do_write_nc.inc"
      end subroutine write_nc_0D
      subroutine write_nc_1D(fid,varname,arr)
      real*4 :: arr(:)
#include "do_write_nc.inc"
      end subroutine write_nc_1D
      subroutine write_nc_2D(fid,varname,arr)
      real*4 :: arr(:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_2D
      subroutine write_nc_3D(fid,varname,arr)
      real*4 :: arr(:,:,:)
#include "do_write_nc.inc" 
      end subroutine write_nc_3D
      subroutine write_nc_4D(fid,varname,arr)
      real*4 :: arr(:,:,:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_4D
      subroutine write_nc_5D(fid,varname,arr)
      real*4 :: arr(:,:,:,:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_5D

      subroutine write_nc_0D_r8(fid,varname,arr)
      real*8 :: arr
#include "do_write_nc_r8.inc"
      end subroutine write_nc_0D_r8
      subroutine write_nc_1D_r8(fid,varname,arr)
      real*8 :: arr(:)
#include "do_write_nc_r8.inc"
      end subroutine write_nc_1D_r8
      subroutine write_nc_2D_r8(fid,varname,arr)
      real*8 :: arr(:,:)
#include "do_write_nc_r8.inc"
      end subroutine write_nc_2D_r8
      subroutine write_nc_3D_r8(fid,varname,arr)
      real*8 :: arr(:,:,:)
#include "do_write_nc_r8.inc" 
      end subroutine write_nc_3D_r8
      subroutine write_nc_4D_r8(fid,varname,arr)
      real*8 :: arr(:,:,:,:)
#include "do_write_nc_r8.inc"
      end subroutine write_nc_4D_r8
      subroutine write_nc_5D_r8(fid,varname,arr)
      real*8 :: arr(:,:,:,:,:)
#include "do_write_nc_r8.inc"
      end subroutine write_nc_5D_r8

      subroutine write_nc_0D_int(fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      integer :: iarr
      real*4 :: arr
      arr = iarr
      call write_data(fid,varname,arr)
      end subroutine write_nc_0D_int
      subroutine write_nc_1D_int(fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      integer :: iarr(:)
      real*4 :: arr(size(iarr))
      arr = iarr
      call write_data(fid,varname,arr)
      end subroutine write_nc_1D_int
      subroutine write_nc_2D_int(fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      integer :: iarr(:,:)
      real*4 :: arr(size(iarr,1),size(iarr,2))
      arr = iarr
      call write_data(fid,varname,arr)
      end subroutine write_nc_2D_int
      subroutine write_nc_3D_int(fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      integer :: iarr(:,:,:)
      real*4 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      arr = iarr
      call write_data(fid,varname,arr)
      end subroutine write_nc_3D_int
      subroutine write_nc_2D_logical(fid,varname,larr)
      integer :: fid
      character(len=*) :: varname
      logical :: larr(:,:)
      real*4 :: arr(size(larr,1),size(larr,2))
      where(larr)
         arr = 1d0
      else where
         arr = 0d0
      end where
      call write_data(fid,varname,arr)
      end subroutine write_nc_2D_logical
      subroutine write_nc_1D_array_of_strings(fid,varname,arr)
      integer :: fid
      character(len=*) :: varname
      character(len=*) :: arr(:)
      integer :: rc,vid
      rc = nf_inq_varid(fid,trim(varname),vid)
      if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &     trim(varname),' not found in output file - stopping'
      rc = nf_put_var_text(fid,vid,arr)
      end subroutine write_nc_1D_array_of_strings

      subroutine defvar_0D(fid,npx,npy,ntiles,
     &     arr,varinfo)
      real*8 :: arr
      integer, parameter :: dtype = nf_float
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_0D
      subroutine defvar_1D(fid,npx,npy,ntiles,
     &     arr,varinfo)
      real*8 :: arr(:)
      integer, parameter :: dtype = nf_float
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_1D
      subroutine defvar_2D(fid,npx,npy,ntiles,
     &     arr,varinfo)
      real*8 :: arr(:,:)
      integer, parameter :: dtype = nf_float
      include 'do_defvar_nc.inc' 
      return
      end subroutine defvar_2D
      subroutine defvar_3D(fid,npx,npy,ntiles,
     &     arr,varinfo)
      real*8 :: arr(:,:,:)
      integer, parameter :: dtype = nf_float
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_3D
      subroutine defvar_4D(fid,npx,npy,ntiles,
     &     arr,varinfo)
      real*8 :: arr(:,:,:,:)
      integer, parameter :: dtype = nf_float
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_4D
      subroutine defvar_5D(fid,npx,npy,ntiles,
     &     arr,varinfo)
      real*8 :: arr(:,:,:,:,:)
      integer, parameter :: dtype = nf_float
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_5D

      subroutine defvar_0D_r8(fid,npx,npy,ntiles,
     &     arr,varinfo,dbl)
      real*8 :: arr
      logical :: dbl
      integer, parameter :: dtype = nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_0D_r8
      subroutine defvar_1D_r8(fid,npx,npy,ntiles,
     &     arr,varinfo,dbl)
      real*8 :: arr(:)
      logical :: dbl
      integer, parameter :: dtype = nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_1D_r8
      subroutine defvar_2D_r8(fid,npx,npy,ntiles,
     &     arr,varinfo,dbl)
      real*8 :: arr(:,:)
      logical :: dbl
      integer, parameter :: dtype = nf_double
      include 'do_defvar_nc.inc' 
      return
      end subroutine defvar_2D_r8
      subroutine defvar_3D_r8(fid,npx,npy,ntiles,
     &     arr,varinfo,dbl)
      real*8 :: arr(:,:,:)
      logical :: dbl
      integer, parameter :: dtype = nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_3D_r8
      subroutine defvar_4D_r8(fid,npx,npy,ntiles,
     &     arr,varinfo,dbl)
      real*8 :: arr(:,:,:,:)
      logical :: dbl
      integer, parameter :: dtype = nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_4D_r8
      subroutine defvar_5D_r8(fid,npx,npy,ntiles,
     &     arr,varinfo,dbl)
      real*8 :: arr(:,:,:,:,:)
      logical :: dbl
      integer, parameter :: dtype = nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_5D_r8

      subroutine defvar_0D_int(fid,npx,npy,ntiles,
     &     arr,varinfo)
      integer :: arr
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_0D_int
      subroutine defvar_1D_int(fid,npx,npy,ntiles,
     &     arr,varinfo)
      integer :: arr(:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_1D_int
      subroutine defvar_2D_int(fid,npx,npy,ntiles,
     &     arr,varinfo)
      integer :: arr(:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_int
      subroutine defvar_3D_int(fid,npx,npy,ntiles,
     &     arr,varinfo)
      integer :: arr(:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_3D_int
      subroutine defvar_4D_int(fid,npx,npy,ntiles,
     &     arr,varinfo)
      integer :: arr(:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_4D_int
      subroutine defvar_5D_int(fid,npx,npy,ntiles,
     &     arr,varinfo)
      integer :: arr(:,:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_5D_int

      subroutine defvar_2D_logical(fid,npx,npy,ntiles,
     &     arr,varinfo)
      logical :: arr(:,:)
c netcdf file will represent logical as 0/1 int
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_logical

      subroutine defvar_2D_char(fid,npx,npy,ntiles,
     &     arr,varinfo)
      character :: arr(:,:)
      integer, parameter :: dtype=nf_char
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_char

      subroutine defvar_1D_array_of_strings(fid,npx,npy,ntiles,
     &     arr,varinfo)
      integer :: fid,npx,npy,ntiles
      character(len=*) :: arr(:)
      character(len=*) :: varinfo
      character, allocatable :: arr2d(:,:)
      allocate(arr2d(len(arr(1)),size(arr,1)))
      call defvar_2D_char(fid,npx,npy,ntiles,
     &     arr2d,varinfo)
      deallocate(arr2d)
      return
      end subroutine defvar_1D_array_of_strings

      subroutine define_var(fid,dtype,
     &     varinfo_in,ndims_in,shp_in,im_in,jm_in,
     &     ntiles,rc,vid)

      integer :: fid,dtype,rc,vid
      character(len=*) :: varinfo_in
      integer :: ndims_in,shp_in(ndims_in),im_in,jm_in,ntiles
      character(len=40) :: vname,dname
      character(len=80) :: varinfo
      character*1, dimension(:), allocatable :: char_arr
      integer :: i,ndims,l,l1,l2,xdim,lv,dsize,status
     &     ,dsizx
      integer :: dids(7)
      logical :: is_dist
      rc = 0
      varinfo=trim(varinfo_in)
c remove spaces
      varinfo=''
      l=0
      do i=1,len_trim(varinfo_in)
        if(varinfo_in(i:i).eq.' ') cycle
        l=l+1
        varinfo(l:l)=varinfo_in(i:i)
      enddo
      lv=l
      allocate(char_arr(len_trim(varinfo)))
      do i=1,len_trim(varinfo)
        char_arr(i)=varinfo(i:i)
      enddo
      l=count(char_arr.eq.'('.or.char_arr.eq.')')
      if(l.eq.0) then
        vname=trim(varinfo)
        ndims=0
      elseif(l.eq.2) then
        vname=varinfo(1:index(varinfo,'(')-1)
        ndims=1+count(char_arr.eq.',')
      endif
      deallocate(char_arr)
      if(l.ne.0 .and. l.ne.2) then
        write(6,*) 'parsing error: ',trim(varinfo)
        rc = 1; return
      endif
      if(ndims.ne.ndims_in) then
        write(6,*) 'parsed ndims does not match actual ndims for ',
     &       trim(vname)
        rc = 1; return
      endif
c check whether variable is already defined
      if(nf_inq_varid(fid,trim(vname),vid).eq.nf_noerr) then
        write(6,*) 'error: variable ',trim(vname),
     &       ' is already defined'
        rc = 1; return
      endif
c
c loop through dimensions and define them if necessary
c
      is_dist = .false.
      if(ndims.gt.0) then
        l1=index(varinfo,'(')+1
        do xdim=1,ndims
          if(xdim.lt.ndims) then
            l2=l1+index(varinfo(l1:lv),',')-2
          else
            l2=l1+index(varinfo(l1:lv),')')-2
          endif
          dname=varinfo(l1:l2)
          dsize = shp_in(xdim)
          if(dname(1:6).eq.'dist_i') then
            dname=dname(6:len_trim(dname))
            dsize = im_in
            is_dist = .true.
          elseif(dname(1:6).eq.'dist_j') then
            dname=dname(6:len_trim(dname))
            dsize = jm_in
            is_dist = .true.
          endif
          if(nf_inq_dimid(fid,trim(dname),dids(xdim)).eq.nf_noerr) then
            status = nf_inq_dimlen(fid,dids(xdim),dsizx)
            if(dsizx.ne.dsize) then
              write(6,*) 'illegal operation: changing dimension size ',
     &             trim(dname),' for variable ',trim(vname)
              rc = 1; return
            endif
          else
            status = nf_def_dim(fid,trim(dname),dsize,dids(xdim))
          endif
          l1=l2+2
        enddo
      endif
c
c if more than one tile, add an extra dimension for distributed vars
c
      if(ntiles.gt.1 .and. is_dist) then
        ndims = ndims + 1
        xdim = ndims
        dname='tile'
        dsize = ntiles
        if(nf_inq_dimid(fid,trim(dname),dids(xdim)).eq.nf_noerr) then
          status = nf_inq_dimlen(fid,dids(xdim),dsizx)
          if(dsizx.ne.dsize) then
            write(6,*) 'illegal operation: changing dimension size ',
     &           trim(dname),' for variable ',trim(vname)
            rc = 1; return
          endif
        else
          status = nf_def_dim(fid,trim(dname),dsize,dids(xdim))
        endif
      endif
c
c define the variable
c
      status = nf_def_var(fid,trim(vname),dtype,ndims,dids,vid)
      return
      end subroutine define_var

      end module ncio
