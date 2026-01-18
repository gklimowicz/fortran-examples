      module pario
!@sum These subroutines are used to read/write arrays (distributed and
!@+ local) to NetCDF files.
!@+
!@+ Simple usage example (to write); the same calls should be made from
!@+ all MPI nodes:
!@+
!@+   fid = par_open(grid,fname,'create')
!@+   call defvar(grid,fid,mdwnimp,'mdwnimp(dist_im,dist_jm)')
!@+   call write_dist_data(grid,fid,'mdwnimp',mdwnimp)
!@+   call par_close(grid, fid)
!@+
!@+ See also Parallelio.F90 for "combined" define/write/read subroutines.
!@+ This allows, in many cases, a SINGLE subroutine to be used to read,
!@+ write and define the contents of a NetCDF file.

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
      use dd2d_utils, only : pack_row,unpack_row,get_nlnk,pack_data
#endif
      implicit none
#ifndef SERIAL_MODE
      include 'mpif.h'
#endif
      save
      private

      include 'netcdf.inc'

! =====================================================================
! I/O Interfaces

      public :: par_open,par_close,par_enddef,variable_exists
     &     ,get_record_dimlen,get_dimlen,get_dimlens,get_record_dimname

      public :: write_dist_data,read_dist_data

      interface write_dist_data
      !@sum subroutine write_dist_data(grid,fid,varname,arr,jdim,no_xdim)
      !@+ ---------------------------------------------------
      !@+ Write a distributed array to a correesponding
      !@+ NetCDF variable on disk (after defvar() has been called).
      !@+
      !@var type(dist_grid), intent(in) :: grid
      !@+     The grid on which the array exists.
      !@+
      !@var integer :: fid
      !@+     Open file handle to write to (obtained via par_open())
      !@+
      !@var character(len=*) :: varname
      !@+     Name of NetCDF variable to write
      !@+
      !@var <type> :: arr(:,:,...)
      !@+     Array to write.
      !@+     <type> may be real*8, integer or logical
      !@+     Number of dimensions must be at least 2
      !@+
      !@+     NOTE: If the desired type/dimension implementation of this
      !@+           interface does not yet exist, it should be added.
      !@+
      !@var integer, intent(in), optional :: jdim = 2
      !@+     Specifies the index (starting from 1) of the LAST horizontal
      !@+     dimension.  If not specified, jdim=2; correct for model arrays
      !@+     like T(i,j,l).  To write an array dimensioned T(l,i,j) set jdim=3.
      !@+
      !@var logical, intent(in), optional :: no_xdim = .false.
      !@+     (WARNING: Negative logic; let has_xdim = .not. no_xdim)
      !@+     ?????
        module procedure par_write_nc_2D
        module procedure par_write_nc_3D
        module procedure par_write_nc_4D
        module procedure par_write_nc_5D
        module procedure par_write_nc_2D_int
        module procedure par_write_nc_3D_int
        module procedure par_write_nc_4D_int
        module procedure par_write_nc_2D_logical
      end interface write_dist_data

      interface read_dist_data
      !@sum subroutine read_dist_data(grid,fid,varname,arr,jdim,no_xdim,record,record1)
      !@+ ---------------------------------------------------
      !@+
      !@+ Reads a distributed array from a NetCDF variable on disk.
      !@+
      !@var type(dist_grid), intent(in) :: grid
      !@+     The grid on which the array exists.
      !@+
      !@var integer :: fid
      !@+     Open file handle to write to (obtained via par_open())
      !@+
      !@var character(len=*) :: varname
      !@+     Name of NetCDF variable to read
      !@+
      !@var <type> :: arr(:,:,...)
      !@+     Array to read.
      !@+     <type> may be real*8, integer or logical
      !@+     Number of dimensions must be at least 2
      !@+
      !@+     NOTE: If the desired type/dimension implementation of this
      !@+           interface does not yet exist, it should be added.
      !@+
      !@var integer, intent(in), optional :: jdim = 2
      !@+     Specifies the index (starting from 1) of the LAST horizontal
      !@+     dimension.  If not specified, jdim=2; correct for model arrays
      !@+     like T(i,j,l).  To write an array dimensioned T(l,i,j) set jdim=3.
      !@+
      !@var logical, intent(in), optional :: no_xdim = .false.
      !@+     (WARNING: Negative logic; let has_xdim = .not. no_xdim)
      !@+     ?????
      !@+
      !@var integer, intent(in), optional :: record,record1
      !@+     ????????

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
      !@sum subroutine write_data(grid,fid,varname,arr)
      !@+ --------------------------------
      !@+
      !@+ Write a non-dstributed array to a correesponding NetCDF variable
      !@+ on disk.  (after defvar() has been called).  The array is
      !@+ written from the root MPI node; the value of arr on other MPI
      !@+ nodes will have no effect on what is written.
      !@+
      !@var type(dist_grid), intent(in) :: grid
      !@+     The grid on which the array exists.
      !@+
      !@var integer :: fid
      !@+     Open file handle to write to (obtained via par_open())
      !@+
      !@var character(len=*) :: varname
      !@+     Name of NetCDF variable to write
      !@+
      !@var <type> :: arr(:,:,...)
      !@+     Array or scalar to write.
      !@+     <type> may be real*8, integer or logical
      !@+
      !@+     NOTE: If the desired type/dimension implementation of this
      !@+           interface does not yet exist, it should be added.
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
      end interface write_data


      interface read_data
      !@sum subroutine read_data(grid,fid,varname,arr,bcast_all)
      !@+ --------------------------------
      !@+
      !@+ Write a non-dstributed array to a correesponding NetCDF variable
      !@+ on disk.  (after defvar() has been called).  The array is
      !@+ written from the root MPI node; the value of arr on other MPI
      !@+ nodes will have no effect on what is written.
      !@+
      !@var type(dist_grid), intent(in) :: grid
      !@+     The grid on which the array exists.
      !@+
      !@var integer :: fid
      !@+     Open file handle to write to (obtained via par_open())
      !@+
      !@var character(len=*) :: varname
      !@+     Name of NetCDF variable to write
      !@+
      !@var <type> :: arr(:,:,...)
      !@+     Array or scalar to write.
      !@+     <type> may be real*8, integer or logical
      !@+
      !@+     NOTE: If the desired type/dimension implementation of this
      !@+           interface does not yet exist, it should be added.
      !@+
      !@var logical, intent(in), optional :: bcast_all = .false.
      !@+     If set, the value read into arr will be broadcast from the
      !@+     root MPI node to MPI nodes.  This is to be used if one
      !@+     wishes to read the same value into arr on all MPI nodes.

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
      !@sum subroutine defvar(grid,fid,arr,varinfo,r4_on_disk,defby)
      !@+ --------------------------------------------------------
      !@+
      !@+ Define a variable in a NetCDF files.  Variables must be defined
      !@+ before they can be written...  This subroutine defines both the
      !@+ variables AND any required dimensions used by the variable.  The
      !@+ sizes of dimensions are inferred from arr, whereas the NetCDF
      !@+ names of variables and dimensions are parsed from varinfo.
      !@+
      !@+ EXAMPLE:
      !@+     integer :: im,jm,lm
      !@+     real(real64), dimension(:,:,:) :: t
      !@+     allocate(t(im,jm,lm))
      !@+     call defvar(grid,fid,t,'t(im,jm,lm)')
      !@+
      !@+ Model variables share dimensions; it is not necessary to declare
      !@+ separate dimension names for each variable.  If a dimension is
      !@+ ever redeclared with a different size than previously, defvar()
      !@+ will abort.
      !@+
      !@+    NOTE: arr and varinfo are reversed, as compared to
      !@+          write_data() and write_dist_data()
      !@+
      !@var type(dist_grid), intent(in) :: grid
      !@+     The grid on which the array exists.
      !@+
      !@var integer :: fid
      !@+     Open file handle to write to (obtained via par_open())
      !@+
      !@var <type> :: arr(:,:,...)
      !@+     Array or scalar to write.
      !@+     <type> may be real*8, integer or logical
      !@+
      !@+     NOTE: If the desired type/dimension implementation of this
      !@+           interface does not yet exist, it should be added.
      !@+
      !@var character(*) :: varinfo
      !@+     String defining name of variable and its dimensions to
      !@+     define in NetCDF.
      !@+     Example: 't(im,jm,lm)'
      !@+
      !@var logical, intent(in), optional :: r4_on_disk
      !@+     Indicates the defined real variable should be a 4-byte float
      !@+     even though the passed array is 8-byte (which is what happens
      !@+     writing out diagnostic acc files).
      !@+
      !@var character(len=*), intent(in), optional :: defby
      !@+     Allows a given variable to have the attribute defby
      !@+     (“defined by”) set to whatever string is passed; this is not
      !@+     frequently used but allows the component “owner” of a
      !@+     variable to be defined if desired.

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

#ifndef SERIAL_MODE
c these routines are only needed when running on multiple CPUs
      interface pack_row_no_xdim
        module procedure pack_row_no_xdim_2d
        module procedure pack_row_no_xdim_3d
        module procedure pack_row_no_xdim_4d
        module procedure pack_row_no_xdim_5d
      end interface
      interface unpack_row_no_xdim
        module procedure unpack_row_no_xdim_2d
        module procedure unpack_row_no_xdim_3d
        module procedure unpack_row_no_xdim_4d
        module procedure unpack_row_no_xdim_5d
      end interface

      interface par_write_jdecomp_optimized
        module procedure par_write_ij
        module procedure par_write_ijx
        module procedure par_write_ijxx
        module procedure par_write_ijxxx
      end interface
#endif /* not SERIAL_MODE */
! End of I/O Interfaces
! =====================================================================

      integer, parameter :: success = 0, fail = -1

      real*8, parameter :: impossible_int=huge(1d0)

      contains

      function par_open(grid,fname,mode)
      type(dist_grid), intent(in) :: grid
      character(len=*) :: fname
      character(len=*) :: mode
      integer :: par_open
      integer :: rc,rc2,fid,vid,wc,idum
      integer :: chunksize
      if(grid%am_i_globalroot) then
        if(trim(mode).eq.'create') then
c          rc = nf_create(trim(fname),nf_clobber,fid)
          rc = nf_create(trim(fname),nf_64bit_offset,fid) ! when files get big
          if(rc.ne.nf_noerr) write(6,*)
     &         'error creating ',trim(fname)
        elseif(trim(mode).eq.'write') then
c          rc = nf_open(trim(fname),nf_write,fid)
          chunksize = 1024*1024*128
          rc = nf__open(trim(fname),nf_write,chunksize,fid)
          if(rc.ne.nf_noerr) write(6,*)
     &         'error opening ',trim(fname)
        elseif(trim(mode).eq.'read') then
c          rc = nf_open(trim(fname),nf_nowrite,fid)
          chunksize = 1024*1024*128
          rc = nf__open(trim(fname),nf_nowrite,chunksize,fid)
          if(rc.ne.nf_noerr) then
            write(6,*) 'error opening ',trim(fname)
          else
            rc2 = nf_inq_varid(fid,'write_status',vid)
            if(rc2.eq.nf_noerr) then
              rc2 = nf_get_var_int(fid,vid,wc)
              if(wc.ne.success) then
                write(6,*) 'input file ',trim(fname),
     &            ' does not appear to have been written successfully:'
                write(6,*) 'write_status = ',wc
              endif
            else
              wc = success
            endif
          endif
        else
          write(6,*) 'par_open: invalid mode ',trim(mode)
          write(6,*) 'mode must be one of [create write read]'
          rc = nf_noerr + 1
        endif
      endif
      call stoprc(rc,nf_noerr)
      if(trim(mode).eq.'read') call stoprc(wc,success)
c define/overwrite the success flag for error checking
      if(grid%am_i_globalroot) then
        if(trim(mode).eq.'create') then
          rc = nf_def_var(fid,'write_status',nf_int,0,idum,vid)
          rc = nf_enddef(fid)
          rc = nf_put_var_int(fid,vid,fail)
          rc = nf_redef(fid)
        elseif(trim(mode).eq.'write') then
          rc = nf_inq_varid(fid,'write_status',vid)
          rc = nf_put_var_int(fid,vid,fail)
          rc = nf_sync(fid)
        endif
        par_open = fid
      else
        par_open = -1
      endif
      return
      end function par_open

      subroutine par_close(grid,fid)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      integer :: rc,vid
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,'write_status',vid)
        if(rc.eq.nf_noerr) then
          rc = nf_put_var_int(fid,vid,success)
        endif
        rc = nf_close(fid)
        if(rc.ne.nf_noerr)
     &       write(6,*) 'error closing file'
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine par_close

      subroutine par_enddef(grid,fid)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      integer :: rc,omode
      if(grid%am_i_globalroot) then
        rc = nf_set_fill(fid,nf_nofill,omode)
        rc = nf_enddef(fid)
      endif
      return
      end subroutine par_enddef

      function variable_exists(grid,fid,varname)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: varname
      logical :: variable_exists
      integer :: vid,rc
      if(grid%am_i_globalroot) rc = nf_inq_varid(fid,trim(varname),vid)
      call broadcast(rc)
      variable_exists = rc == nf_noerr
      return
      end function variable_exists

      subroutine get_record_dimname(grid,fid,dname)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: dname
c
      integer :: rc,unlim_did,ierr,l,ll
      integer, parameter :: slen=64
      character :: dname_(slen)
      character(len=slen) :: dname__
c
      if(grid%am_i_globalroot) then
        rc = nf_inq_unlimdim(fid,unlim_did)
        dname__ = ''
        rc = nf_inq_dimname(fid,unlim_did,dname__)
        ll = len_trim(dname__)
        do l=1,ll
          dname_(l) = dname__(l:l)
        enddo
        if(rc.ne.nf_noerr) then
          if(grid%am_i_globalroot)
     &     write(6,*) 'get_record_dimname: input file has no record '//
     &         ' dimension - stopping'
        endif
      endif
      call stoprc(rc,nf_noerr)
#ifndef SERIAL_MODE
      call mpi_bcast(dname_,slen,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif
      call broadcast(ll)
      dname=''
      do l=1,ll
        dname(l:l) = dname_(l)
      enddo
      return
      end subroutine get_record_dimname

      function get_record_dimlen(grid,fid,checkvar)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*), optional :: checkvar
      integer :: get_record_dimlen
      integer :: i
c
      integer :: rc,vid,unlim_did,ndims,dids(7)
      if(grid%am_i_globalroot) then
        rc = nf_inq_unlimdim(fid,unlim_did)
        rc = nf_inq_dimlen(fid,unlim_did,i)
        if(rc.ne.nf_noerr) then
          if(grid%am_i_globalroot)
     &     write(6,*) 'get_record_dimlen: input file has no record '//
     &         ' dimension - stopping'
        endif
      endif
      call stoprc(rc,nf_noerr)
      if(present(checkvar)) then
        if(grid%am_i_globalroot) then
          rc = nf_inq_varid(fid,trim(checkvar),vid)
          if(rc.ne.nf_noerr) then
            if(grid%am_i_globalroot)
     &           write(6,*) 'get_record_dimlen: variable ',
     &           trim(checkvar),' not found in input file - stopping'
          endif
        endif
        call stoprc(rc,nf_noerr)
        if(grid%am_i_globalroot) then
          rc = nf_inq_varndims(fid,vid,ndims)
          rc = nf_inq_vardimid(fid,vid,dids)
          if(dids(ndims).ne.unlim_did) then
            if(grid%am_i_globalroot)
     &           write(6,*) 'get_record_dimlen: variable ',
     &           trim(checkvar),' has no record dim - stopping'
            rc = 0
          else
            rc = 1
          endif
        endif
        call stoprc(rc,1)
      endif
      call broadcast(i)
      get_record_dimlen = i
      return
      end function get_record_dimlen

      function get_dimlen(grid,fid,dname)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: dname
      integer :: get_dimlen
c
      integer :: rc,rc2,vid,ndims,dids(7),dlen

      if(grid%am_i_globalroot) then
        rc2 = 1
        rc = nf_inq_dimid(fid,trim(dname),dids(1))
        if(rc.ne.nf_noerr) then
          write(6,*) 'get_dimlen: dimension ',
     &         trim(dname),' not found in input file - stopping'
          rc2 = 0
        endif
        rc = nf_inq_dimlen(fid,dids(1),dlen)
      endif
      call stoprc(rc2,1)
      call broadcast(dlen)
      get_dimlen = dlen
      return
      end function get_dimlen

      subroutine get_dimlens(grid,fid,vname,ndims,dlens)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: vname
      integer :: ndims,dlens(:)
c
      integer :: rc,rc2,vid,dids(7),idim,unlim_did

      if(grid%am_i_globalroot) then
        rc2 = 1
        rc = nf_inq_varid(fid,trim(vname),vid)
        if(rc.ne.nf_noerr) then
          write(6,*) 'get_dimlen: variable ',
     &         trim(vname),' not found in input file - stopping'
          rc2 = 0
        else
          rc = nf_inq_varndims(fid,vid,ndims)
          rc = nf_inq_vardimid(fid,vid,dids)
          do idim=1,ndims
            rc = nf_inq_dimlen(fid,dids(idim),dlens(idim))
          enddo
          if(ndims.ge.3 .and. grid%ntiles.eq.6) then
            ! tile dimension is either the last or next to last
            rc = nf_inq_unlimdim(fid,unlim_did)
            if(dids(ndims).eq.unlim_did) dlens(ndims-1) = dlens(ndims)
            ndims = ndims - 1
          endif
        endif
      endif
      call stoprc(rc2,1)
      call broadcast(ndims)
      call broadcast(dlens)
      return
      end subroutine get_dimlens

      subroutine par_write_nc_2D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_2D
      subroutine par_write_nc_3D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_3D
      subroutine par_write_nc_4D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:,:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_4D
      subroutine par_write_nc_5D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:,:,:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_5D

      subroutine par_read_nc_2D(grid,fid,varname,arr,jdim,no_xdim,
     &     record,record1)
      real*8 :: arr(:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_2D
      subroutine par_read_nc_3D(grid,fid,varname,arr,jdim,no_xdim,
     &     record,record1)
      real*8 :: arr(:,:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_3D
      subroutine par_read_nc_4D(grid,fid,varname,arr,jdim,no_xdim,
     &     record,record1)
      real*8 :: arr(:,:,:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_4D
      subroutine par_read_nc_5D(grid,fid,varname,arr,jdim,no_xdim,
     &     record,record1)
      real*8 :: arr(:,:,:,:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_5D

      subroutine par_read_nc_2D_int(grid,fid,varname,iarr,
     &     jdim,no_xdim,record,record1)
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
     &     no_xdim,record,record1)
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
     &   no_xdim, record, record1)
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

      subroutine par_write_nc_2D_int(grid,fid,varname,iarr, jdim,
     &     no_xdim)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2))
      integer, optional, intent(in) :: jdim
      logical, intent(in), optional :: no_xdim
      arr = iarr
      call write_dist_data(grid,fid,varname,arr)
      end subroutine par_write_nc_2D_int

      subroutine par_write_nc_3D_int(grid,fid,varname,iarr,jdim,no_xdim)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      integer, intent(in), optional :: jdim
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      logical, intent(in), optional :: no_xdim
      arr = iarr
      if(present(jdim)) then
        call write_dist_data(grid,fid,varname,arr,jdim=jdim)
      else
        call write_dist_data(grid,fid,varname,arr)
      endif
      end subroutine par_write_nc_3D_int

      subroutine par_write_nc_4D_int(grid,fid,varname,iarr,jdim,no_xdim)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:,:)
      integer, intent(in), optional :: jdim
      logical, intent(in), optional :: no_xdim
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3),size(iarr,4))
      arr = iarr
      if(present(jdim)) then
        call write_dist_data(grid,fid,varname,arr,jdim=jdim)
      else
        call write_dist_data(grid,fid,varname,arr)
      endif
      end subroutine par_write_nc_4D_int

      subroutine par_read_nc_2D_logical(grid,fid,varname,larr,jdim,
     &     no_xdim, record, record1)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical :: larr(:,:)
      real*8 :: arr(size(larr,1),size(larr,2))
      integer, optional, intent(in) :: jdim
      logical, intent(in), optional :: no_xdim
      integer, intent(in), optional :: record,record1
      arr = 0.d0
      call read_dist_data(grid,fid,varname,arr)
      larr = arr.eq.1d0
      end subroutine par_read_nc_2D_logical
      subroutine par_write_nc_2D_logical(grid,fid,varname,larr,jdim,
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
      call mpi_bcast(rc,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      if(rc.ne.rc_ok) call mpi_abort(MPI_COMM_WORLD,1,mpi_err)
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

      subroutine write_nc_0D(grid,fid,varname,arr)
      real*8 :: arr
#include "do_write_nc.inc"
      end subroutine write_nc_0D
      subroutine write_nc_1D(grid,fid,varname,arr)
      real*8 :: arr(:)
#include "do_write_nc.inc"
      end subroutine write_nc_1D
      subroutine write_nc_2D(grid,fid,varname,arr)
      real*8 :: arr(:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_2D
      subroutine write_nc_3D(grid,fid,varname,arr)
      real*8 :: arr(:,:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_3D
      subroutine write_nc_4D(grid,fid,varname,arr)
      real*8 :: arr(:,:,:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_4D
      subroutine write_nc_5D(grid,fid,varname,arr)
      real*8 :: arr(:,:,:,:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_5D

      subroutine read_nc_0D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr
#include "do_read_nc.inc"
      end subroutine read_nc_0D
      subroutine read_nc_1D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:)
#include "do_read_nc.inc"
      end subroutine read_nc_1D
      subroutine read_nc_2D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:)
#include "do_read_nc.inc"
      end subroutine read_nc_2D
      subroutine read_nc_3D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:,:)
#include "do_read_nc.inc"
      end subroutine read_nc_3D
      subroutine read_nc_4D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:,:,:)
#include "do_read_nc.inc"
      end subroutine read_nc_4D
      subroutine read_nc_5D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:,:,:,:)
#include "do_read_nc.inc"
      end subroutine read_nc_5D

      subroutine write_nc_0D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr
      real*8 :: arr
      if(grid%am_i_globalroot) arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_0D_int
      subroutine write_nc_1D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:)
      real*8 :: arr(size(iarr))
      if(grid%am_i_globalroot) arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_1D_int
      subroutine write_nc_2D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2))
      if(grid%am_i_globalroot) arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_2D_int
      subroutine write_nc_3D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      if(grid%am_i_globalroot) arr = iarr
      call write_data(grid,fid,varname,arr)
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
      integer :: rc,vid
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,trim(varname),vid)
        if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in output file - stopping'
      endif
      call stoprc(rc,nf_noerr)
      if(grid%am_i_globalroot) then
        rc = nf_put_var_text(fid,vid,arr)
      endif
      end subroutine write_nc_1D_array_of_strings

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

      subroutine defvar_0D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_0D
      subroutine defvar_1D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_1D
      subroutine defvar_2D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D
      subroutine defvar_3D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_3D
      subroutine defvar_4D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:,:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_4D
      subroutine defvar_5D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:,:,:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_5D

      subroutine defvar_0D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_0D_int
      subroutine defvar_1D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_1D_int
      subroutine defvar_2D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_int
      subroutine defvar_3D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_3D_int
      subroutine defvar_4D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_4D_int
      subroutine defvar_5D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:,:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_5D_int

      subroutine defvar_2D_logical(grid,fid,arr,varinfo,r4_on_disk,
     &     defby)
      logical :: arr(:,:)
c netcdf file will represent logical as 0/1 int
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_logical

      subroutine defvar_2D_char(grid,fid,arr,varinfo,r4_on_disk,
     &     defby)
      character :: arr(:,:)
      integer, parameter :: dtype=nf_char
      include 'do_defvar_nc.inc'
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

      subroutine write_attr_text(grid,fid,varname,attname,attval)
      character(len=*) :: attval
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_text(fid,vid,trim(attname),attlen,attval)
c        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_text
      subroutine write_attr_0D_int(grid,fid,varname,attname,attval)
      integer :: attval
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_int(fid,vid,trim(attname),nf_int,attlen,attval)
c        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_0D_int
      subroutine write_attr_1D_int(grid,fid,varname,attname,attval)
      integer :: attval(:)
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_int(fid,vid,trim(attname),nf_int,attlen,attval)
c        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_1D_int
      subroutine write_attr_0D_r8(grid,fid,varname,attname,attval)
      real*8 :: attval
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_double(fid,vid,trim(attname),nf_double,attlen,
     &       attval)
c        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_0D_r8
      subroutine write_attr_1D_r8(grid,fid,varname,attname,attval)
      real*8 :: attval(:)
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_double(fid,vid,trim(attname),nf_double,attlen,
     &       attval)
c        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_1D_r8

      subroutine read_attr_text(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      character(len=*) :: attval
      integer :: kpos
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_text(fid,vid,trim(attname),tmpstr)
        do kpos=1,attlen ! remove extra NULL characters
          if(iachar(tmpstr(kpos)).eq.0) tmpstr(kpos)=' '
        enddo
      endif
      call stoprc(rc,nf_noerr)
#ifndef SERIAL_MODE
      call mpi_bcast(tmpstr,attlen,MPI_CHARACTER,0,
     &     MPI_COMM_WORLD,ierr)
#endif
      attval=''
      do l=1,attlen
        attval(l:l) = tmpstr(l)
      enddo
      return
      end subroutine read_attr_text
      subroutine read_attr_0D_int(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      integer :: attval
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_int(fid,vid,trim(attname),attval)
      endif
      call stoprc(rc,nf_noerr)
      call broadcast(attval)
      return
      end subroutine read_attr_0D_int
      subroutine read_attr_1D_int(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      integer :: attval(:)
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_int(fid,vid,trim(attname),attval)
      endif
      call stoprc(rc,nf_noerr)
      call broadcast(attval)
      return
      end subroutine read_attr_1D_int
      subroutine read_attr_0D_r8(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      real*8 :: attval
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_double(fid,vid,trim(attname),attval)
      endif
      call stoprc(rc,nf_noerr)
      call broadcast(attval)
      return
      end subroutine read_attr_0D_r8
      subroutine read_attr_1D_r8(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      real*8 :: attval(:)
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_double(fid,vid,trim(attname),attval)
      endif
      call stoprc(rc,nf_noerr)
      call broadcast(attval)
      return
      end subroutine read_attr_1D_r8

      subroutine get_natts(grid,fid,varname,natts)
      type(dist_grid) :: grid
      integer :: fid
      character(len=*) :: varname
      integer :: natts
      integer :: rc,vid
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,trim(varname),vid)
        if(rc.ne.nf_noerr)
     &       write(6,*) 'error: nonexistent variable ',trim(varname)
      endif
      call stoprc(rc,nf_noerr)
      if(grid%am_i_globalroot)
     &     rc = nf_inq_varnatts(fid,vid,natts)
      call broadcast(natts)
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

#ifndef SERIAL_MODE
      subroutine pack_row_no_xdim_2d(grid,arr,arr1d,jdim)
      real*8 arr(:,:)
      include 'row_setup_no_xdim.inc'
      call copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      return
      end subroutine pack_row_no_xdim_2d
      subroutine pack_row_no_xdim_3d(grid,arr,arr1d,jdim)
      real*8 arr(:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      return
      end subroutine pack_row_no_xdim_3d
      subroutine pack_row_no_xdim_4d(grid,arr,arr1d,jdim)
      real*8 arr(:,:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      return
      end subroutine pack_row_no_xdim_4d
      subroutine pack_row_no_xdim_5d(grid,arr,arr1d,jdim)
      real*8 arr(:,:,:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      return
      end subroutine pack_row_no_xdim_5d

      subroutine unpack_row_no_xdim_2d(grid,arr1d,arr,jdim)
      real*8 arr(:,:)
      include 'row_setup_no_xdim.inc'
      call copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      return
      end subroutine unpack_row_no_xdim_2d
      subroutine unpack_row_no_xdim_3d(grid,arr1d,arr,jdim)
      real*8 arr(:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      return
      end subroutine unpack_row_no_xdim_3d
      subroutine unpack_row_no_xdim_4d(grid,arr1d,arr,jdim)
      real*8 arr(:,:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      return
      end subroutine unpack_row_no_xdim_4d
      subroutine unpack_row_no_xdim_5d(grid,arr1d,arr,jdim)
      real*8 arr(:,:,:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      return
      end subroutine unpack_row_no_xdim_5d

      subroutine par_write_ij(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgij(:,:)
      integer :: rc
      if(grid%am_i_globalroot) then
        allocate(arrgij(grid%npx,grid%npy))
      else
        allocate(arrgij(1,1))
      endif
      call pack_data(grid,arr,arrgij)
      if(grid%am_i_globalroot) then
        rc = nf_put_var_double(fid,vid,arrgij)
      endif
      deallocate(arrgij)
      return
      end subroutine par_write_ij
      subroutine par_write_ijx(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgij(:,:)
      integer :: rc,k,srt(3),cnt(3)
      if(grid%am_i_globalroot) then
        allocate(arrgij(grid%npx,grid%npy))
      else
        allocate(arrgij(1,1))
      endif
      srt(1:2) = 1; cnt(1:3) = (/ grid%npx, grid%npy, 1 /)
      do k=1,size(arr,3)
        call pack_data(grid,arr(:,:,k),arrgij)
        srt(3) = k
        if(grid%am_i_globalroot)
     &       rc = nf_put_vara_double(fid,vid,srt,cnt,arrgij)
      enddo
      deallocate(arrgij)
      return
      end subroutine par_write_ijx
      subroutine par_write_ijxx(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:,:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgij(:,:)
      integer :: rc,k,l,srt(4),cnt(4)
      if(jdim.eq.3) then
        call par_write_xijx(grid,fid,vid,arr,jdim)
        return
      endif
      if(grid%am_i_globalroot) then
        allocate(arrgij(grid%npx,grid%npy))
      else
        allocate(arrgij(1,1))
      endif
      srt(1:2) = 1; cnt(1:4) = (/ grid%npx, grid%npy, 1, 1 /)
      do l=1,size(arr,4)
      srt(4) = l
      do k=1,size(arr,3)
        call pack_data(grid,arr(:,:,k,l),arrgij)
        srt(3) = k
        if(grid%am_i_globalroot)
     &       rc = nf_put_vara_double(fid,vid,srt,cnt,arrgij)
      enddo
      enddo
      deallocate(arrgij)
      return
      end subroutine par_write_ijxx
      subroutine par_write_ijxxx(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:,:,:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgij(:,:)
      integer :: rc,k,l,m,srt(5),cnt(5)
      if(grid%am_i_globalroot) then
        allocate(arrgij(grid%npx,grid%npy))
      else
        allocate(arrgij(1,1))
      endif
      srt(1:2) = 1; cnt(1:5) = (/ grid%npx, grid%npy, 1, 1, 1 /)
      do m=1,size(arr,5)
      srt(5) = m
      do l=1,size(arr,4)
      srt(4) = l
      do k=1,size(arr,3)
        call pack_data(grid,arr(:,:,k,l,m),arrgij)
        srt(3) = k
        if(grid%am_i_globalroot)
     &       rc = nf_put_vara_double(fid,vid,srt,cnt,arrgij)
      enddo
      enddo
      enddo
      deallocate(arrgij)
      return
      end subroutine par_write_ijxxx
      subroutine par_write_xijx(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:,:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgxij(:,:,:)
      integer :: rc,k,srt(4),cnt(4)
      if(grid%am_i_globalroot) then
        allocate(arrgxij(size(arr,1),grid%npx,grid%npy))
      else
        allocate(arrgxij(1,1,1))
      endif
      srt(1:3) = 1; cnt(1:4) = (/ size(arr,1), grid%npx, grid%npy, 1 /)
      do k=1,size(arr,4)
        call pack_data(grid,arr(:,:,:,k),arrgxij,jdim=3)
        srt(4) = k
        if(grid%am_i_globalroot)
     &       rc = nf_put_vara_double(fid,vid,srt,cnt,arrgxij)
      enddo
      deallocate(arrgxij)
      return
      end subroutine par_write_xijx
#endif /* not SERIAL_MODE */

      end module pario

      subroutine copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      implicit none
      integer :: nl,nj,nk,j1,j2
      real*8 arr1d(*)
      real*8 arr(nl,nj,nk)

      integer :: j,k,l,n

      n = 0
      do k=1,nk
      do j=j1,j2
      do l=1,nl
        n = n + 1
        arr1d(n) = arr(l,j,k)
      enddo
      enddo
      enddo
      return
      end subroutine copy_to_1D_no_xdim
      subroutine copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      implicit none
      integer :: nl,nj,nk,j1,j2
      real*8 arr1d(*)
      real*8 arr(nl,nj,nk)

      integer :: j,k,l,n

      n = 0
      do k=1,nk
      do j=j1,j2
      do l=1,nl
        n = n + 1
        arr(l,j,k) = arr1d(n)
      enddo
      enddo
      enddo
      return
      end subroutine copy_from_1D_no_xdim
