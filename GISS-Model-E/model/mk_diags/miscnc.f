!@sum Miscellaneous utilities for netcdf files.
!@auth M. Kelley
      subroutine copy_shared_vars(ifid,ofid)
      implicit none
      integer :: ifid,ofid
      include 'netcdf.inc'
      integer :: status,ivarid,ovarid,nvars,vtype,v4size,
     &     vsize_input
      character(len=80) :: vname
      real*4, dimension(:), allocatable :: v4
      status = nf_inq_nvars(ofid,nvars)
      do ovarid=1,nvars
         status = nf_inq_varname(ofid,ovarid,vname)
         status = nf_inq_varid(ifid,vname,ivarid)
         if(status.ne.nf_noerr) cycle
         status = nf_inq_vartype(ofid,ovarid,vtype)
         call get_varsize(ofid,vname,v4size)
         call get_varsize(ifid,vname,vsize_input)
         if(v4size.ne.vsize_input) cycle
         if(vtype.eq.nf_double) then ! real*8 data
           allocate(v4(2*v4size))
           status = nf_get_var_double(ifid,ivarid,v4)
           status = nf_put_var_double(ofid,ovarid,v4)
         elseif(vtype.eq.nf_char) then ! character data
           allocate(v4(v4size)) ! larger than necessary
           status = nf_get_var_text(ifid,ivarid,v4)
           status = nf_put_var_text(ofid,ovarid,v4)
         else                   ! floats and integers
           allocate(v4(v4size))
           status = nf_get_var_real(ifid,ivarid,v4)
           status = nf_put_var_real(ofid,ovarid,v4)
         endif
         deallocate(v4)
      enddo
      return
      end subroutine copy_shared_vars

      subroutine get_dimsize(fid,dim_name,dimsize)
      implicit none
      include 'netcdf.inc'
      integer :: fid,dimsize
      character(len=*) :: dim_name
      integer :: status,dimid
      status = nf_inq_dimid(fid,dim_name,dimid)
      if(status.ne.nf_noerr) then
        write(6,*) 'nonexistent dimension ',trim(dim_name)
        stop
      endif
      status = nf_inq_dimlen(fid,dimid,dimsize)
      return
      end subroutine get_dimsize

      subroutine get_varsize(fid,var_name,var_size)
      implicit none
      integer :: fid
      character(len=*) :: var_name
      integer :: var_size
      integer :: ndims,dimsizes(7)
      call get_vdimsizes(fid,var_name,ndims,dimsizes)
      if(ndims.eq.0) then
        var_size = 1
      else
        var_size = product(dimsizes(1:ndims))
      endif
      return
      end subroutine get_varsize

      subroutine get_varsize8(fid,var_name,var_size)
      implicit none
      integer :: fid
      character(len=*) :: var_name
      integer*8 :: var_size
      integer :: ndims,dimsizes(7)
      call get_vdimsizes(fid,var_name,ndims,dimsizes)
      if(ndims.eq.0) then
        var_size = 1
      else
        var_size = product(dimsizes(1:ndims))
      endif
      return
      end subroutine get_varsize8

      subroutine get_vdimsizes(fid,var_name,ndims,dimsizes)
      implicit none
      include 'netcdf.inc'
      integer :: fid
      character(len=*) :: var_name
      integer :: ndims,dimsizes(7)
      integer :: status,varid,dimids(7),dimlen,n
      status = nf_inq_varid(fid,var_name,varid)
      status = nf_inq_varndims(fid,varid,ndims)
      status = nf_inq_vardimid(fid,varid,dimids)
      do n=1,ndims
        status = nf_inq_dimlen(fid,dimids(n),dimsizes(n))
      enddo
      return
      end subroutine get_vdimsizes

      subroutine get_var_real(fid,var_name,var)
      implicit none
      include 'netcdf.inc'
      integer :: fid
      character(len=*) :: var_name
      real*4 :: var(1)
      integer :: status,varid
      status = nf_inq_varid(fid,var_name,varid)
      if(status.ne.nf_noerr) then
        write(6,*) 'nonexistent variable ',trim(var_name)
        stop
      endif
      status = nf_get_var_real(fid,varid,var)
      return
      end subroutine get_var_real

      subroutine put_var_real(fid,var_name,var)
      implicit none
      include 'netcdf.inc'
      integer :: fid
      character(len=*) :: var_name
      real*4 :: var(1)
      integer :: status,varid
      status = nf_inq_varid(fid,var_name,varid)
      if(status.ne.nf_noerr) then
        write(6,*) 'nonexistent variable ',trim(var_name)
        stop
      else
        status = nf_put_var_real(fid,varid,var)
      endif
      return
      end subroutine put_var_real

      subroutine get_var_int(fid,var_name,var)
      implicit none
      include 'netcdf.inc'
      integer :: fid
      character(len=*) :: var_name
      integer :: var(1)
      integer :: status,varid
      status = nf_inq_varid(fid,var_name,varid)
      if(status.ne.nf_noerr) then
        write(6,*) 'nonexistent variable ',trim(var_name)
        stop
      endif
      status = nf_get_var_int(fid,varid,var)
      return
      end subroutine get_var_int

      subroutine put_var_int(fid,var_name,var)
      implicit none
      include 'netcdf.inc'
      integer :: fid
      character(len=*) :: var_name
      integer :: var(1)
      integer :: status,varid
      status = nf_inq_varid(fid,var_name,varid)
      if(status.ne.nf_noerr) then
        write(6,*) 'nonexistent variable ',trim(var_name)
        stop
      else
        status = nf_put_var_int(fid,varid,var)
      endif
      return
      end subroutine put_var_int

      subroutine get_var_text(fid,var_name,var)
      implicit none
      include 'netcdf.inc'
      integer :: fid
      character(len=*) :: var_name
      character :: var(1)
      integer :: status,varid
      status = nf_inq_varid(fid,var_name,varid)
      if(status.ne.nf_noerr) then
        write(6,*) 'nonexistent variable ',trim(var_name)
        stop
      endif
      status = nf_get_var_text(fid,varid,var)
      return
      end subroutine get_var_text

      subroutine copy_file_structure(ifid,ofid)
      implicit none
      include 'netcdf.inc'
      integer :: ifid,ofid
      integer :: status,ndims,idimid,dimsiz
      character(len=40) :: dimname,vname,att_name
      integer :: ivarid,nvars,vtype,dimids(7),natts,ovid,n
c copy dimensions
      status = nf_inq_ndims(ifid,ndims)
      do idimid=1,ndims
        status = nf_inq_dimlen(ifid,idimid,dimsiz)
        status = nf_inq_dimname(ifid,idimid,dimname)
        status = nf_def_dim(ofid,dimname,dimsiz,idimid)
      enddo
c copy variables
      status = nf_inq_nvars(ifid,nvars)
      do ivarid=1,nvars
        status = nf_inq_var(ifid,ivarid,vname,vtype,ndims,dimids,natts)
        status = nf_def_var(ofid,vname,vtype,ndims,dimids,ovid)
        do n=1,natts
          status = nf_inq_attname(ifid,ivarid,n,att_name)
          status = nf_copy_att(ifid,ivarid,att_name,ofid,ovid)
        enddo
      enddo
c copy global attributes
      status = nf_inq_varnatts(ifid,nf_global,natts)
      do n=1,natts
        status = nf_inq_attname(ifid,nf_global,n,att_name)
        status = nf_copy_att(ifid,nf_global,att_name,ofid,nf_global)
      enddo
      status = nf_enddef(ofid)
      return
      end subroutine copy_file_structure

      subroutine handle_err(status,errmsg)
      implicit none
      integer :: status
      character(len=*) :: errmsg
      include 'netcdf.inc'
      if(status.ne.nf_noerr) then
        write(6,*) 'error '//trim(errmsg)
        stop
      endif
      end subroutine handle_err
