      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer, intent(in), optional :: jdim
      integer, intent(in), optional :: record
      logical, intent(in), optional :: no_xdim
#ifndef SERIAL_MODE
      real*8, dimension(:), allocatable :: arr_nohalo
      integer :: nl,nk,jd,nd,ierr,vid,did,rc,arrsize,td,rd
      integer*8, dimension(7) :: srt,cnt
      integer*8 :: nrecs8
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(grid%am_i_globalroot .and. rc.ne.nf_noerr)
     &     write(6,*) 'variable ',
     &     trim(varname),' not found in output file - stopping'
      call stoprc(rc,nf_noerr)
      if(present(jdim)) then
        jd=jdim
      else
        jd=2
      endif
      nd = size(shape(arr))
      td = nd+1
      cnt(1:nd) = shape(arr)
      srt(1:nd) = 1
      srt(jd-1:jd) = (/ grid%is, grid%js /)
      cnt(jd-1:jd) = 1+(/ grid%ie-grid%is, grid%je-grid%js /)
      srt(td) = grid%tile
      cnt(td) = 1
      if(present(record)) then
        nrecs8 = 0
        rc = nfmpi_inq_unlimdim(fid,did)
        rc = nfmpi_inq_dimlen(fid,did,nrecs8)
        if(record.le.0 .or. nrecs8+1.lt.record) then
          if(grid%am_i_globalroot) write(6,*)
     &         'error in record dim spec. for variable ',trim(varname)
          call stoprc(0,1)
        endif
        rd = td
        if(grid%ntiles.gt.1) rd = rd + 1
        srt(rd) = record
        cnt(rd) = 1
      endif
      if(grid%have_domain) then
        call get_nlnk(arr,jd,nl,nk)
      else
        srt=1; cnt=0   ! empty domain
      endif
      if(size(arr,jd).ne.cnt(jd)) then ! halo present
        arrsize = product(cnt(1:nd))
        allocate(arr_nohalo(arrsize))
        if(arrsize.gt.0) then
          call copy_to_1D( arr,arr_nohalo,nl,nk,
     &         grid%isd,grid%ied,grid%jsd,grid%jed,
     &         grid%is,grid%ie,grid%js,grid%je )
        endif
        rc = nfmpi_put_vara_double_all(fid,vid,srt,cnt,arr_nohalo)
        deallocate(arr_nohalo)
      else
        rc = nfmpi_put_vara_double_all(fid,vid,srt,cnt,arr)
      endif
#else
c call write_data when in serial mode
      call write_data(grid,fid,varname,arr)
#endif
