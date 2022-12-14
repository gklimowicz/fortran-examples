      module pario
!@sum A version of the pario module built on the Parallel-NetCDF
!@+ library provided by ANL.
!@+   http://www.mcs.anl.gov/parallel-netcdf
!@+
!@+ This module has the same name and interface as pario_nc.f
!@+ See pario_nc.f for API documentation.


#ifdef OFFLINE
#else
c see whether model E is running in serial mode
#ifndef USE_MPI
#define SERIAL_MODE
#endif
#endif

      use dd2d_utils, only : dist_grid
#ifndef SERIAL_MODE
c these routines are only needed when running on multiple CPUs
      use dd2d_utils, only : get_nlnk
#endif
      implicit none
      save
      private

#ifndef SERIAL_MODE
      include 'mpif.h'
#endif
      include 'pnetcdf.inc'

c
c i/o interfaces
c
      public :: par_open,par_close,par_enddef,variable_exists
     &     ,get_record_dimlen,get_dimlen,get_dimlens,get_record_dimname

      public :: write_dist_data,read_dist_data
      interface write_dist_data
        module procedure par_write_nc_2D
        module procedure par_write_nc_3D
        module procedure par_write_nc_4D
        module procedure par_write_nc_5D
        module procedure par_write_nc_2D_int
        module procedure par_write_nc_3D_int
        module procedure par_write_nc_4D_int
        module procedure par_write_nc_2D_logical
        module procedure par_write_nc_3D_bundle
        module procedure par_write_nc_4D_bundle
        module procedure par_write_nc_5D_bundle
      end interface write_dist_data
      interface read_dist_data
        module procedure par_read_nc_2D
        module procedure par_read_nc_3D
        module procedure par_read_nc_4D
        module procedure par_read_nc_5D
        module procedure par_read_nc_2D_int
        module procedure par_read_nc_3D_int
        module procedure par_read_nc_4D_int
        module procedure par_read_nc_2D_logical
      end interface read_dist_data

      public :: write_data,read_data
      interface write_data
        module procedure write_nc_0D
        module procedure write_nc_1D
        module procedure write_nc_2D
        module procedure write_nc_3D
        module procedure write_nc_4D
        module procedure write_nc_5D
        module procedure write_nc_0D_int
        module procedure write_nc_1D_int
        module procedure write_nc_2D_int
        module procedure write_nc_3D_int
        module procedure write_nc_2D_logical
        module procedure write_nc_1D_array_of_strings
        module procedure write_nc_string
      end interface write_data
      interface read_data
        module procedure read_nc_0D
        module procedure read_nc_1D
        module procedure read_nc_2D
        module procedure read_nc_3D
        module procedure read_nc_4D
        module procedure read_nc_5D
        module procedure read_nc_0D_int
        module procedure read_nc_1D_int
        module procedure read_nc_2D_int
        module procedure read_nc_3D_int
        module procedure read_nc_2D_logical
      end interface read_data

      public :: defvar
      interface defvar
        module procedure defvar_0D
        module procedure defvar_1D
        module procedure defvar_2D
        module procedure defvar_3D
        module procedure defvar_4D
        module procedure defvar_5D
        module procedure defvar_0D_int
        module procedure defvar_1D_int
        module procedure defvar_2D_int
        module procedure defvar_3D_int
        module procedure defvar_4D_int
        module procedure defvar_5D_int
        module procedure defvar_2D_logical
        module procedure defvar_1D_array_of_strings
        module procedure defvar_string
      end interface

      public :: write_attr
      interface write_attr
        module procedure write_attr_text
        module procedure write_attr_0D_r8
        module procedure write_attr_1D_r8
        module procedure write_attr_0D_int
        module procedure write_attr_1D_int
      end interface

      public :: read_attr
      interface read_attr
        module procedure read_attr_text
        module procedure read_attr_0D_r8
        module procedure read_attr_1D_r8
        module procedure read_attr_0D_int
        module procedure read_attr_1D_int
      end interface

      public :: get_natts

      interface len_of_obj
        module procedure len_of_text
        module procedure len_of_int0D
        module procedure len_of_int1D
        module procedure len_of_r80D
        module procedure len_of_r81D
      end interface

      interface full_len_of_obj
        module procedure full_len_of_text
        module procedure len_of_int0D
        module procedure len_of_int1D
        module procedure len_of_r80D
        module procedure len_of_r81D
      end interface

      interface broadcast
        module procedure broadcast_0D_int
        module procedure broadcast_1D_int
        module procedure broadcast_0D_r8
        module procedure broadcast_1D_r8
      end interface broadcast

      integer, parameter :: success = 0, fail = -1

      real*8, parameter :: impossible_int=huge(1d0)

      public :: set_record_dimname

      contains

      function par_open(grid,fname,mode)
      type(dist_grid), intent(in) :: grid
      character(len=*) :: fname
      character(len=*) :: mode
      integer :: par_open
      integer :: rc,rc2,fid,vid,wc,idum
      INTEGER :: comm, info, mpierror
      logical :: am_root
      am_root = grid%am_i_globalroot
      comm = MPI_COMM_WORLD ! world for now
      info = MPI_INFO_NULL  ! for now

c create mpi hints
c      info=99
c      call MPI_INFO_CREATE(info,mpierror)
c      call MPI_INFO_SET(info,'nc_header_align_size', '1048567',mpierror)
c      call MPI_INFO_SET(info,'nc_var_align_size', '4194304',mpierror)
c      call MPI_INFO_SET(info,'striping_unit', '4194304',mpierror)

      if(trim(mode).eq.'create') then
        rc = nfmpi_create(comm,trim(fname),nf_64bit_offset,info,fid)
        if(am_root .and. rc.ne.nf_noerr) write(6,*)
     &       'error creating ',trim(fname)
      elseif(trim(mode).eq.'write') then
        rc = nfmpi_open(comm,trim(fname),nf_write,info,fid)
        if(am_root .and. rc.ne.nf_noerr) write(6,*)
     &       'error opening ',trim(fname)
      elseif(trim(mode).eq.'read') then
        rc = nfmpi_open(comm,trim(fname),nf_nowrite,info,fid)
        if(rc.ne.nf_noerr) then
          if(am_root) write(6,*) 'error opening ',trim(fname)
        else
          rc2 = nfmpi_inq_varid(fid,'write_status',vid)
          if(rc2.eq.nf_noerr) then
            rc2 = nfmpi_get_var_int_all(fid,vid,wc)
            if(am_root .and. wc.ne.success) then
              write(6,*) 'input file ',trim(fname),
     &             ' does not appear to have been written successfully:'
              write(6,*) 'write_status = ',wc
            endif
          else
            wc = success
          endif
        endif
      else
        if(am_root) then
          write(6,*) 'par_open: invalid mode ',trim(mode)
          write(6,*) 'mode must be one of [create write read]'
        endif
        rc = nf_noerr + 1
      endif
      call stoprc(rc,nf_noerr)
      if(trim(mode).eq.'read') call stoprc(wc,success)
c define/overwrite the success flag for error checking
      if(trim(mode).eq.'create') then
        rc = nfmpi_def_var(fid,'write_status',nf_int,0,idum,vid)
        rc = nfmpi_enddef(fid)
        call write_data(grid,fid,'write_status',fail)
        rc = nfmpi_redef(fid)
      elseif(trim(mode).eq.'write') then
        call write_data(grid,fid,'write_status',fail)
c        rc = nfmpi_sync(fid)
      endif
      par_open = fid
      return
      end function par_open

      subroutine par_close(grid,fid)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      integer :: rc,vid
      rc = nfmpi_inq_varid(fid,'write_status',vid)
      if(rc.eq.nf_noerr) then
        call write_data(grid,fid,'write_status',success)
      endif
      rc = nfmpi_close(fid)
      return
      end subroutine par_close

      subroutine par_enddef(grid,fid)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      integer :: rc
      rc = nfmpi_enddef(fid)
      return
      end subroutine par_enddef

      function variable_exists(grid,fid,varname)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: varname
      logical :: variable_exists
      integer :: vid
      variable_exists =
     &     nfmpi_inq_varid(fid,trim(varname),vid) == nf_noerr
      return
      end function variable_exists

      subroutine get_record_dimname(grid,fid,dname)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: dname
c
      integer :: rc,unlim_did
      character(len=64) :: dname_
c
      rc = nfmpi_inq_unlimdim(fid,unlim_did)
      rc = nfmpi_inq_dimname(fid,unlim_did,dname_)
      dname = trim(dname_)
      return
      end subroutine get_record_dimname

      function get_record_dimlen(grid,fid,checkvar)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*), optional :: checkvar
      integer :: get_record_dimlen
      integer*8 :: i8
c
      integer :: rc,vid,unlim_did,ndims,dids(7)
      rc = nfmpi_inq_unlimdim(fid,unlim_did)
      rc = nfmpi_inq_dimlen(fid,unlim_did,i8)
      if(rc.ne.nf_noerr) then
        if(grid%am_i_globalroot)
     &       write(6,*) 'get_record_dimlen: input file has no record '//
     &       ' dimension - stopping'
        call stoprc(0,1)
      endif
      if(present(checkvar)) then
        rc = nfmpi_inq_varid(fid,trim(checkvar),vid)
        if(rc.ne.nf_noerr) then
          if(grid%am_i_globalroot)
     &         write(6,*) 'get_record_dimlen: variable ',
     &         trim(checkvar),' not found in input file - stopping'
          call stoprc(0,1)
        endif
        rc = nfmpi_inq_varndims(fid,vid,ndims)
        rc = nfmpi_inq_vardimid(fid,vid,dids)
        if(dids(ndims).ne.unlim_did) then
          if(grid%am_i_globalroot)
     &         write(6,*) 'get_record_dimlen: variable ',
     &         trim(checkvar),' has no record dim - stopping'
          call stoprc(0,1)
        endif
      endif
      get_record_dimlen = i8
      return
      end function get_record_dimlen

      function get_dimlen(grid,fid,dname)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: dname
      integer :: get_dimlen
c
      integer*8 :: i8
      integer :: rc,vid,ndims,dids(7)
c
      rc = nfmpi_inq_dimid(fid,trim(dname),dids(1))
      if(rc.ne.nf_noerr) then
        if(grid%am_i_globalroot)
     &       write(6,*) 'get_dimlen: dimension ',
     &       trim(dname),' not found in input file - stopping'
        call stoprc(0,1)
      endif
      rc = nfmpi_inq_dimlen(fid,dids(1),i8)
      get_dimlen = i8
      return
      end function get_dimlen

      subroutine get_dimlens(grid,fid,vname,ndims,dlens)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: vname
      integer :: ndims,dlens(:)
c
      integer*8 :: i8
      integer :: idim,rc,vid,dids(7),unlim_did

      rc = nfmpi_inq_varid(fid,trim(vname),vid)
      if(rc.ne.nf_noerr) then
        if(grid%am_i_globalroot)
     &       write(6,*) 'get_dimlen: variable ',
     &       trim(vname),' not found in input file - stopping'
        call stoprc(0,1)
      endif
      rc = nfmpi_inq_varndims(fid,vid,ndims)
      rc = nfmpi_inq_vardimid(fid,vid,dids)
      do idim=1,ndims
        rc = nfmpi_inq_dimlen(fid,dids(idim),i8)
        dlens(idim) = i8
      enddo
      if(ndims.ge.3 .and. grid%ntiles.eq.6) then
        ! tile dimension is either the last or next to last
        rc = nfmpi_inq_unlimdim(fid,unlim_did)
        if(dids(ndims).eq.unlim_did) dlens(ndims-1) = dlens(ndims)
        ndims = ndims - 1
      endif
      return
      end subroutine get_dimlens

      subroutine par_write_nc_2D(grid,fid,varname,arr,jdim,record,
     &     no_xdim)
      real*8 :: arr(:,:)
#include "do_par_write_pnc.inc"
      end subroutine par_write_nc_2D
      subroutine par_write_nc_3D(grid,fid,varname,arr,jdim,record,
     &     no_xdim)
      real*8 :: arr(:,:,:)
#include "do_par_write_pnc.inc"
      end subroutine par_write_nc_3D
      subroutine par_write_nc_4D(grid,fid,varname,arr,jdim,record,
     &     no_xdim)
      real*8 :: arr(:,:,:,:)
#include "do_par_write_pnc.inc"
      end subroutine par_write_nc_4D
      subroutine par_write_nc_5D(grid,fid,varname,arr,jdim,record,
     &     no_xdim)
      real*8 :: arr(:,:,:,:,:)
#include "do_par_write_pnc.inc"
      end subroutine par_write_nc_5D

      subroutine par_write_nc_3D_bundle(grid,fid,varnames,arr,jdim,
     &     record)
      real*8 :: arr(:,:,:)
#include "do_par_write_bundle.inc"
      end subroutine par_write_nc_3D_bundle
      subroutine par_write_nc_4D_bundle(grid,fid,varnames,arr,jdim,
     &     record)
      real*8 :: arr(:,:,:,:)
#include "do_par_write_bundle.inc"
      end subroutine par_write_nc_4D_bundle
      subroutine par_write_nc_5D_bundle(grid,fid,varnames,arr,jdim,
     &     record)
      real*8 :: arr(:,:,:,:,:)
#include "do_par_write_bundle.inc"
      end subroutine par_write_nc_5D_bundle

      subroutine par_read_nc_2D(grid,fid,varname,arr,jdim,
     &     record,record1,no_xdim)
      real*8 :: arr(:,:)
#include "do_par_read_pnc.inc"
      end subroutine par_read_nc_2D
      subroutine par_read_nc_3D(grid,fid,varname,arr,jdim,
     &     record,record1,no_xdim)
      real*8 :: arr(:,:,:)
#include "do_par_read_pnc.inc"
      end subroutine par_read_nc_3D
      subroutine par_read_nc_4D(grid,fid,varname,arr,jdim,
     &     record,record1,no_xdim)
      real*8 :: arr(:,:,:,:)
#include "do_par_read_pnc.inc"
      end subroutine par_read_nc_4D
      subroutine par_read_nc_5D(grid,fid,varname,arr,jdim,
     &     record,record1,no_xdim)
      real*8 :: arr(:,:,:,:,:)
#include "do_par_read_pnc.inc"
      end subroutine par_read_nc_5D

      subroutine par_read_nc_2D_int(grid,fid,varname,iarr,
     &  no_xdim, jdim, record, record1)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2))
      integer, intent(in), optional :: jdim
      logical, intent(in), optional :: no_xdim
      integer, intent(in), optional :: record,record1
      arr = impossible_int
      call read_dist_data(grid,fid,varname,arr)
      where(arr.ne.impossible_int) iarr = arr
      end subroutine par_read_nc_2D_int

      subroutine par_read_nc_3D_int(grid,fid,varname,iarr,jdim,
     &     no_xdim, record, record1)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      integer, intent(in), optional :: jdim
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      logical, intent(in), optional :: no_xdim
      integer, intent(in), optional :: record,record1
      arr = impossible_int
      if(present(jdim)) then
        call read_dist_data(grid,fid,varname,arr,jdim=jdim)
      else
        call read_dist_data(grid,fid,varname,arr)
      endif
      where(arr.ne.impossible_int) iarr = arr
      end subroutine par_read_nc_3D_int
      subroutine par_read_nc_4D_int(grid,fid,varname,iarr,jdim,
     &     no_xdim, record, record1)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:,:)
      integer, intent(in), optional :: jdim
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3),size(iarr,4))
      logical, intent(in), optional :: no_xdim
      integer, intent(in), optional :: record,record1
      arr = impossible_int
      if(present(jdim)) then
        call read_dist_data(grid,fid,varname,arr,jdim=jdim)
      else
        call read_dist_data(grid,fid,varname,arr)
      endif
      where(arr.ne.impossible_int) iarr = arr
      end subroutine par_read_nc_4D_int

      subroutine par_write_nc_2D_int(grid,fid,varname,iarr,record,
     &     jdim, no_xdim)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:)
      integer, intent(in), optional :: record
      real*8 :: arr(size(iarr,1),size(iarr,2))
      integer, optional, intent(in) :: jdim
      logical, optional, intent(in) :: no_xdim
      arr = iarr
      if(present(record)) then
        call write_dist_data(grid,fid,varname,arr,record=record)
      else
        call write_dist_data(grid,fid,varname,arr)
      endif
      end subroutine par_write_nc_2D_int
      subroutine par_write_nc_3D_int(grid,fid,varname,iarr,jdim,record,
     &     no_xdim)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      integer, intent(in), optional :: jdim
      integer, intent(in), optional :: record
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      logical, optional, intent(in) :: no_xdim
      integer :: jdim_
      arr = iarr
      jdim_ = 2
      if(present(jdim)) jdim_ = jdim
      if(present(record)) then
        call write_dist_data(grid,fid,varname,arr,jdim=jdim_,
     &       record=record)
      else
        call write_dist_data(grid,fid,varname,arr,jdim=jdim_)
      endif
      end subroutine par_write_nc_3D_int
      subroutine par_write_nc_4D_int(grid,fid,varname,iarr,jdim,record,
     &     no_xdim)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:,:)
      integer, intent(in), optional :: jdim
      integer, intent(in), optional :: record
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3),size(iarr,4))
      logical, optional, intent(in) :: no_xdim
      integer :: jdim_
      arr = iarr
      jdim_ = 2
      if(present(jdim)) jdim_ = jdim
      if(present(record)) then
        call write_dist_data(grid,fid,varname,arr,jdim=jdim_,
     &       record=record)
      else
        call write_dist_data(grid,fid,varname,arr,jdim=jdim_)
      endif
      end subroutine par_write_nc_4D_int

      subroutine par_read_nc_2D_logical(grid,fid,varname,larr,
     &     jdim, no_xdim, record, record1)
      USE CONSTANT, only : UNDEF_VAL
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical :: larr(:,:)
      real*8 :: arr(size(larr,1),size(larr,2))
      integer, optional, intent(in) :: jdim
      logical, intent(in), optional :: no_xdim
      integer, optional, intent(in) :: record, record1
      arr = UNDEF_VAL
      call read_dist_data(grid,fid,varname,arr)
      larr = arr.eq.1d0
      end subroutine par_read_nc_2D_logical

      subroutine par_write_nc_2D_logical(grid,fid,varname,larr, jdim,
     &     no_xdim)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical :: larr(:,:)
      integer, optional, intent(in) :: jdim
      logical, intent(in), optional :: no_xdim
      real*8 :: arr(size(larr,1),size(larr,2))
      where(larr)
        arr = 1d0
      else where
        arr = 0d0
      end where
      call write_dist_data(grid,fid,varname,arr)
      end subroutine par_write_nc_2D_logical

      subroutine stoprc(rc,rc_ok)
      integer :: rc,rc_ok
      integer :: mpi_err
#ifndef SERIAL_MODE
c      call mpi_bcast(rc,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      if(rc.ne.rc_ok) then
        call mpi_finalize(mpi_err)
        call mpi_abort(MPI_COMM_WORLD,1,mpi_err)
      endif
#else
      if(rc.ne.rc_ok) stop
#endif
      return
      end subroutine stoprc

      subroutine broadcast_0D_int(i)
      integer :: i
#ifndef SERIAL_MODE
      integer :: ierr
      call mpi_bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      end subroutine broadcast_0D_int
      subroutine broadcast_1D_int(i)
      integer :: i(:)
#ifndef SERIAL_MODE
      integer :: ierr
      call mpi_bcast(i,size(i),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      end subroutine broadcast_1D_int
      subroutine broadcast_0D_r8(r8)
      real*8 :: r8
#ifndef SERIAL_MODE
      integer :: ierr
      call mpi_bcast(r8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
      end subroutine broadcast_0D_r8
      subroutine broadcast_1D_r8(r8)
      real*8 :: r8(:)
#ifndef SERIAL_MODE
      integer :: ierr
      call mpi_bcast(r8,size(r8),MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_WORLD,ierr)
#endif
      end subroutine broadcast_1D_r8

      subroutine write_nc_0D(grid,fid,varname,arr,record)
      real*8 :: arr
#include "do_write_pnc.inc"
      end subroutine write_nc_0D
      subroutine write_nc_1D(grid,fid,varname,arr,record)
      real*8 :: arr(:)
#include "do_write_pnc.inc"
      end subroutine write_nc_1D
      subroutine write_nc_2D(grid,fid,varname,arr,record)
      real*8 :: arr(:,:)
#include "do_write_pnc.inc"
      end subroutine write_nc_2D
      subroutine write_nc_3D(grid,fid,varname,arr,record)
      real*8 :: arr(:,:,:)
#include "do_write_pnc.inc"
      end subroutine write_nc_3D
      subroutine write_nc_4D(grid,fid,varname,arr,record)
      real*8 :: arr(:,:,:,:)
#include "do_write_pnc.inc"
      end subroutine write_nc_4D
      subroutine write_nc_5D(grid,fid,varname,arr,record)
      real*8 :: arr(:,:,:,:,:)
#include "do_write_pnc.inc"
      end subroutine write_nc_5D

      subroutine read_nc_0D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr
#include "do_read_pnc.inc"
      end subroutine read_nc_0D
      subroutine read_nc_1D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:)
#include "do_read_pnc.inc"
      end subroutine read_nc_1D
      subroutine read_nc_2D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:)
#include "do_read_pnc.inc"
      end subroutine read_nc_2D
      subroutine read_nc_3D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:,:)
#include "do_read_pnc.inc"
      end subroutine read_nc_3D
      subroutine read_nc_4D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:,:,:)
#include "do_read_pnc.inc"
      end subroutine read_nc_4D
      subroutine read_nc_5D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:,:,:,:)
#include "do_read_pnc.inc"
      end subroutine read_nc_5D

      subroutine write_nc_0D_int(grid,fid,varname,iarr,record)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr
      integer, intent(in), optional :: record
      real*8 :: arr
      if(grid%am_i_globalroot) arr = iarr
      if(present(record)) then
        call write_data(grid,fid,varname,arr,record=record)
      else
        call write_data(grid,fid,varname,arr)
      endif
      end subroutine write_nc_0D_int
      subroutine write_nc_1D_int(grid,fid,varname,iarr,record)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:)
      integer, intent(in), optional :: record
      real*8 :: arr(size(iarr))
      if(grid%am_i_globalroot) arr = iarr
      if(present(record)) then
        call write_data(grid,fid,varname,arr,record=record)
      else
        call write_data(grid,fid,varname,arr)
      endif
      end subroutine write_nc_1D_int
      subroutine write_nc_2D_int(grid,fid,varname,iarr,record)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:)
      integer, intent(in), optional :: record
      real*8 :: arr(size(iarr,1),size(iarr,2))
      if(grid%am_i_globalroot) arr = iarr
      if(present(record)) then
        call write_data(grid,fid,varname,arr,record=record)
      else
        call write_data(grid,fid,varname,arr)
      endif
      end subroutine write_nc_2D_int
      subroutine write_nc_3D_int(grid,fid,varname,iarr,record)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      integer, intent(in), optional :: record
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      if(grid%am_i_globalroot) arr = iarr
      if(present(record)) then
        call write_data(grid,fid,varname,arr,record=record)
      else
        call write_data(grid,fid,varname,arr)
      endif
      end subroutine write_nc_3D_int
      subroutine write_nc_2D_logical(grid,fid,varname,larr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical :: larr(:,:)
      real*8 :: arr(size(larr,1),size(larr,2))
      if(grid%am_i_globalroot) then
        where(larr)
          arr = 1d0
        else where
          arr = 0d0
        end where
      endif
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_2D_logical
      subroutine write_nc_1D_array_of_strings(grid,fid,varname,arr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      character(len=*) :: arr(:)
      integer*8 :: srt(2),cnt(2)
      integer :: rc,vid
      character :: char_dum
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(grid%am_i_globalroot .and. rc.ne.nf_noerr)
     &     write(6,*) 'variable ',
     &     trim(varname),' not found in output file - stopping'
      call stoprc(rc,nf_noerr)
      srt(:) = 1
      if(grid%am_i_globalroot) then
        cnt(:) = (/ len(arr(1)), size(arr,1) /)
        rc = nfmpi_put_vara_text_all(fid,vid,srt,cnt,arr)
      else
        cnt(:) = 0
        char_dum = 'x'
        rc = nfmpi_put_vara_text_all(fid,vid,srt,cnt,char_dum)
      endif
      end subroutine write_nc_1D_array_of_strings
      subroutine write_nc_string(grid,fid,varname,arr,record)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      character(len=*) :: arr
      integer, intent(in), optional :: record
      integer :: rc,vid,did
      integer*8 :: srt(2),cnt(2),nrecs8
      character :: char_dum
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(grid%am_i_globalroot .and. rc.ne.nf_noerr)
     &     write(6,*) 'variable ',
     &     trim(varname),' not found in output file - stopping'
      call stoprc(rc,nf_noerr)
      srt(:) = 1
      cnt(:) = 0
      if(present(record)) then
        nrecs8 = 0
        rc = nfmpi_inq_unlimdim(fid,did)
        rc = nfmpi_inq_dimlen(fid,did,nrecs8)
        if(record.le.0 .or. nrecs8+1.lt.record) then
          if(grid%am_i_globalroot) write(6,*)
     &         'error in record dim spec. for variable ',trim(varname)
          call stoprc(0,1)
        endif
        srt(2) = record
        cnt(2) = 1
      endif
      if(grid%am_i_globalroot) then
        cnt(1) = len(arr)
        rc = nfmpi_put_vara_text_all(fid,vid,srt,cnt,arr)
      else
        char_dum = 'x'
        rc = nfmpi_put_vara_text_all(fid,vid,srt,cnt,char_dum)
      endif
      end subroutine write_nc_string

      subroutine read_nc_0D_int(grid,fid,varname,iarr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr
      logical, intent(in), optional :: bcast_all
      real*8 :: arr
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      arr = impossible_int
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      if(arr.ne.impossible_int) iarr = arr
      end subroutine read_nc_0D_int
      subroutine read_nc_1D_int(grid,fid,varname,iarr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:)
      logical, intent(in), optional :: bcast_all
      real*8 :: arr(size(iarr))
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      arr = impossible_int
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      where(arr.ne.impossible_int) iarr = arr
      end subroutine read_nc_1D_int
      subroutine read_nc_2D_int(grid,fid,varname,iarr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:)
      logical, intent(in), optional :: bcast_all
      real*8 :: arr(size(iarr,1),size(iarr,2))
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      arr = impossible_int
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      where(arr.ne.impossible_int) iarr = arr
      end subroutine read_nc_2D_int
      subroutine read_nc_3D_int(grid,fid,varname,iarr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      logical, intent(in), optional :: bcast_all
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      arr = impossible_int
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      where(arr.ne.impossible_int) iarr = arr
      end subroutine read_nc_3D_int
      subroutine read_nc_2D_logical(grid,fid,varname,larr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical :: larr(:,:)
      logical, intent(in), optional :: bcast_all
      real*8 :: arr(size(larr,1),size(larr,2))
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      if(grid%am_i_globalroot .or. bc_all) larr = arr.eq.1d0
      end subroutine read_nc_2D_logical

      subroutine defvar_0D(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      real*8 :: arr
      integer, parameter :: dtype=nf_double
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_0D
      subroutine defvar_1D(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      real*8 :: arr(:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_1D
      subroutine defvar_2D(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      real*8 :: arr(:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_2D
      subroutine defvar_3D(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      real*8 :: arr(:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_3D
      subroutine defvar_4D(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      real*8 :: arr(:,:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_4D
      subroutine defvar_5D(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      real*8 :: arr(:,:,:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_5D

      subroutine defvar_0D_int(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      integer :: arr
      integer, parameter :: dtype=nf_int
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_0D_int
      subroutine defvar_1D_int(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      integer :: arr(:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_1D_int
      subroutine defvar_2D_int(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      integer :: arr(:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_2D_int
      subroutine defvar_3D_int(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      integer :: arr(:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_3D_int
      subroutine defvar_4D_int(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      integer :: arr(:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_4D_int
      subroutine defvar_5D_int(grid,fid,arr,varinfo,r4_on_disk,defby,
     &     with_record_dim)
      integer :: arr(:,:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_5D_int

      subroutine defvar_2D_logical(grid,fid,arr,varinfo,r4_on_disk,
     &     defby,with_record_dim)
      logical :: arr(:,:)
c netcdf file will represent logical as 0/1 int
      integer, parameter :: dtype=nf_int
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_2D_logical

      subroutine defvar_2D_char(grid,fid,arr,varinfo,r4_on_disk,
     &     defby, with_record_dim)
      character :: arr(:,:)
      integer, parameter :: dtype=nf_char
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_2D_char

      subroutine defvar_1D_array_of_strings(grid,fid,arr,varinfo)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: arr(:)
      character(len=*) :: varinfo
      character, allocatable :: arr2d(:,:)
      allocate(arr2d(len(arr(1)),size(arr,1)))
      call defvar_2D_char(grid,fid,arr2d,varinfo)
      deallocate(arr2d)
      return
      end subroutine defvar_1D_array_of_strings

      subroutine defvar_1D_char(grid,fid,arr,varinfo,r4_on_disk,
     &     defby, with_record_dim)
      character :: arr(:)
      integer, parameter :: dtype=nf_char
      include 'do_defvar_pnc.inc'
      return
      end subroutine defvar_1D_char

      subroutine defvar_string(grid,fid,arr,varinfo,with_record_dim)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: arr
      character(len=*) :: varinfo
      logical, intent(in), optional :: with_record_dim
      character, allocatable :: arr1d(:)
      allocate(arr1d(len(arr)))
      if(present(with_record_dim)) then
        call defvar_1D_char(grid,fid,arr1d,varinfo,
     &       with_record_dim=with_record_dim)
      else
        call defvar_1D_char(grid,fid,arr1d,varinfo)
      endif
      deallocate(arr1d)
      return
      end subroutine defvar_string

      subroutine write_attr_text(grid,fid,varname,attname,attval)
      character(len=*) :: attval
#include "setup_attput_pnc.inc"
      rc = nfmpi_put_att_text(fid,vid,trim(attname),attlen8,attval)
c        if(do_enddef) rc2 = nfmpi_enddef(fid)
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_text
      subroutine write_attr_0D_int(grid,fid,varname,attname,attval)
      integer :: attval
#include "setup_attput_pnc.inc"
      rc = nfmpi_put_att_int(fid,vid,trim(attname),nf_int,attlen8,
     &     attval)
c        if(do_enddef) rc2 = nfmpi_enddef(fid)
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_0D_int
      subroutine write_attr_1D_int(grid,fid,varname,attname,attval)
      integer :: attval(:)
#include "setup_attput_pnc.inc"
      rc = nfmpi_put_att_int(fid,vid,trim(attname),nf_int,attlen8,
     &     attval)
c        if(do_enddef) rc2 = nfmpi_enddef(fid)
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_1D_int
      subroutine write_attr_0D_r8(grid,fid,varname,attname,attval)
      real*8 :: attval
#include "setup_attput_pnc.inc"
      rc = nfmpi_put_att_double(fid,vid,trim(attname),nf_double,
     &     attlen8,attval)
c        if(do_enddef) rc2 = nfmpi_enddef(fid)
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_0D_r8
      subroutine write_attr_1D_r8(grid,fid,varname,attname,attval)
      real*8 :: attval(:)
#include "setup_attput_pnc.inc"
      rc = nfmpi_put_att_double(fid,vid,trim(attname),nf_double,
     &     attlen8,attval)
c        if(do_enddef) rc2 = nfmpi_enddef(fid)
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_1D_r8

      subroutine read_attr_text(grid,fid,varname,attname,attlen,
     &     attval,attnum)
      character(len=*) :: attval
      integer :: kpos
#include "setup_attget_pnc.inc"
      attval=''
      rc = nfmpi_get_att_text(fid,vid,trim(attname),attval)
      call stoprc(rc,nf_noerr)
      do kpos=1,len_trim(attval) ! remove extra NULL characters
        if(iachar(attval(kpos:kpos)).eq.0) attval(kpos:kpos)=' '
      enddo
      return
      end subroutine read_attr_text
      subroutine read_attr_0D_int(grid,fid,varname,attname,attlen,
     &     attval,attnum)
      integer :: attval
#include "setup_attget_pnc.inc"
      rc = nfmpi_get_att_int(fid,vid,trim(attname),attval)
      call stoprc(rc,nf_noerr)
      return
      end subroutine read_attr_0D_int
      subroutine read_attr_1D_int(grid,fid,varname,attname,attlen,
     &     attval,attnum)
      integer :: attval(:)
#include "setup_attget_pnc.inc"
      rc = nfmpi_get_att_int(fid,vid,trim(attname),attval)
      call stoprc(rc,nf_noerr)
      return
      end subroutine read_attr_1D_int
      subroutine read_attr_0D_r8(grid,fid,varname,attname,attlen,
     &     attval,attnum)
      real*8 :: attval
#include "setup_attget_pnc.inc"
      rc = nfmpi_get_att_double(fid,vid,trim(attname),attval)
      call stoprc(rc,nf_noerr)
      return
      end subroutine read_attr_0D_r8
      subroutine read_attr_1D_r8(grid,fid,varname,attname,attlen,
     &     attval,attnum)
      real*8 :: attval(:)
#include "setup_attget_pnc.inc"
      rc = nfmpi_get_att_double(fid,vid,trim(attname),attval)
      call stoprc(rc,nf_noerr)
      return
      end subroutine read_attr_1D_r8

      subroutine get_natts(grid,fid,varname,natts)
      type(dist_grid) :: grid
      integer :: fid
      character(len=*) :: varname
      integer :: natts
      integer :: rc,vid
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      call stoprc(rc,nf_noerr)
      rc = nfmpi_inq_varnatts(fid,vid,natts)
      return
      end subroutine get_natts

      function len_of_text(cstr)
      integer :: len_of_text
      character(len=*) :: cstr
      len_of_text = len_trim(cstr)
      return
      end function len_of_text
      function full_len_of_text(cstr)
      integer :: full_len_of_text
      character(len=*) :: cstr
      full_len_of_text = len(cstr)
      return
      end function full_len_of_text
      function len_of_int0D(i)
      integer :: len_of_int0D
      integer :: i
      len_of_int0D = 1
      return
      end function len_of_int0D
      function len_of_int1D(i)
      integer :: len_of_int1D
      integer :: i(:)
      len_of_int1D = size(i)
      return
      end function len_of_int1D
      function len_of_r80D(r8)
      integer :: len_of_r80D
      real*8 :: r8
      len_of_r80D = 1
      return
      end function len_of_r80D
      function len_of_r81D(r8)
      integer :: len_of_r81D
      real*8 :: r8(:)
      len_of_r81D = size(r8)
      return
      end function len_of_r81D

      subroutine define_var(fid,dtype,
     &     varinfo_in,ndims_in,shp_in,im_in,jm_in,
     &     ntiles,rc,vid,am_root,with_record_dim)
      integer :: fid,dtype,rc,vid
      character(len=*) :: varinfo_in
      integer :: ndims_in,shp_in(ndims_in),im_in,jm_in,ntiles
      logical :: am_root
      character(len=40) :: vname,dname
      character(len=80) :: varinfo
      logical :: with_record_dim
      character*1, dimension(:), allocatable :: char_arr
      integer :: i,ndims,l,l1,l2,xdim,lv,dsize,status
     &     ,dsizx,num_dist_dims,dist_dim_count
      integer :: dids(7)
      integer*8 :: dsize8
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
        if(am_root) write(6,*) 'parsing error: ',trim(varinfo)
        rc = 1; return
      endif
      if(ndims.ne.ndims_in) then
        if(am_root) write(6,*)
     &       'parsed ndims does not match actual ndims for ',
     &       trim(vname)
        rc = 1; return
      endif
c check whether variable is already defined
      if(nfmpi_inq_varid(fid,trim(vname),vid).eq.nf_noerr) then
        if(am_root) write(6,*) 'error: variable ',trim(vname),
     &       ' is already defined'
        rc = 1; return
      endif
c
c count the number of distributed dimensions
c
      num_dist_dims = 0
      if(index(varinfo,'dist_').gt.0) then
        do i=1,len_trim(varinfo)-4
          if(varinfo(i:i+4).eq.'dist_') num_dist_dims = num_dist_dims+1
        enddo
        if(num_dist_dims.gt.2) then
          if(am_root) write(6,*)
     &         'error: number of distributed dimensions for variable ',
     &         trim(vname),' exceeds 2'
          rc = 1; return
        endif
      endif
c
c loop through dimensions and define them if necessary
c
      dist_dim_count = 0
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
          if(dname(1:5).eq.'dist_') then
            dname=dname(6:len_trim(dname))
            dist_dim_count = dist_dim_count+1
            if(dist_dim_count.lt.num_dist_dims) then
              dsize = im_in
            else
              dsize = jm_in
            endif
          endif
          if(nfmpi_inq_dimid(fid,trim(dname),dids(xdim)).eq.nf_noerr)
     &         then
            status = nfmpi_inq_dimlen(fid,dids(xdim),dsize8)
            dsizx = dsize8
            if(dsizx.ne.dsize) then
              if(am_root) write(6,*)
     &             'illegal operation: changing dimension size ',
     &             trim(dname),' for variable ',trim(vname)
              rc = 1; return
            endif
          else
            dsize8 = dsize
            status = nfmpi_def_dim(fid,trim(dname),dsize8,dids(xdim))
          endif
          l1=l2+2
        enddo
      endif
c
c if more than one tile, add an extra dimension for distributed vars
c
      if(ntiles.gt.1 .and. num_dist_dims.gt.0) then
        ndims = ndims + 1
        xdim = ndims
        dname='tile'
        dsize = ntiles
        if(nfmpi_inq_dimid(fid,trim(dname),dids(xdim)).eq.nf_noerr) then
          status = nfmpi_inq_dimlen(fid,dids(xdim),dsize8)
          dsizx = dsize8
          if(dsizx.ne.dsize) then
            if(am_root) write(6,*)
     &           'illegal operation: changing dimension size ',
     &           trim(dname),' for variable ',trim(vname)
            rc = 1; return
          endif
        else
          dsize8 = dsize
          status = nfmpi_def_dim(fid,trim(dname),dsize8,dids(xdim))
        endif
      endif
c
c if record dimension requested, add it
c
      if(with_record_dim) then
        ndims = ndims + 1
        status = nfmpi_inq_unlimdim(fid,dids(ndims))
        if(dids(ndims).le.0) then ! record dim needs defining
          status = nfmpi_def_dim(fid,'record',nfmpi_unlimited,
     &         dids(ndims))
        endif
      endif
c
c define the variable
c
      status = nfmpi_def_var(fid,trim(vname),dtype,ndims,dids,vid)
      return
      end subroutine define_var

      subroutine set_record_dimname(grid,fid,record_dimname)
      type(dist_grid), intent(in) :: grid
      integer, intent(in) :: fid
      character(len=*), intent(in) :: record_dimname
      integer :: rc,did
      rc = nfmpi_inq_unlimdim(fid,did)
      if(did.gt.0) then
        if(grid%am_i_globalroot) write(6,*)
     &       'error in set_record_dimname: the record dimension ',
     &       'has already been named'
        call stoprc(0,1)
      endif
      rc = nfmpi_def_dim(fid,trim(record_dimname),nfmpi_unlimited,did)
      return
      end subroutine set_record_dimname

      end module pario

      subroutine copy_to_1D(arr,arr1d,nl,nk,
     &     isd,ied,jsd,jed, is,ie,js,je)
      implicit none
      real*8 arr(nl,isd:ied,jsd:jed,nk)
      real*8 arr1d(*)
      integer :: nl,nk, isd,ied,jsd,jed, is,ie,js,je
      integer :: i,j,k,l,n
      n = 0
      do k=1,nk
        do j=js,je
          do i=is,ie
            do l=1,nl
              n = n + 1
              arr1d(n) = arr(l,i,j,k)
            enddo
          enddo
        enddo
      enddo
      return
      end subroutine copy_to_1D

      subroutine copy_from_1D(arr1d,arr,nl,nk,
     &     isd,ied,jsd,jed, is,ie,js,je)
      implicit none
      real*8 arr(nl,isd:ied,jsd:jed,nk)
      real*8 arr1d(*)
      integer :: nl,nk, isd,ied,jsd,jed, is,ie,js,je
      integer :: i,j,k,l,n
      n = 0
      do k=1,nk
        do j=js,je
          do i=is,ie
            do l=1,nl
              n = n + 1
              arr(l,i,j,k) = arr1d(n)
            enddo
          enddo
        enddo
      enddo
      return
      end subroutine copy_from_1D

      subroutine copy_0D(in,out)
      implicit none
      real*8 in,out
      out = in
      return
      end subroutine copy_0D
