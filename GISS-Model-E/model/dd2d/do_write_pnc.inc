      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: varname
      integer, intent(in), optional :: record
      integer :: vid,did,rc,nd,ierr
      integer*8, dimension(7) :: srt,cnt
      integer*8 :: nrecs8
      real*8 :: arr_dum
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(grid%am_i_globalroot .and. rc.ne.nf_noerr) then
        write(6,*) 'variable ',
     &       trim(varname),' not found in output file - stopping'
      endif
      call stoprc(rc,nf_noerr)
      nd = size(shape(arr))
      srt(:) = 1
      cnt(:) = 0
      if(nd.eq.0) then ! to use collective i/o on zero-dimensional data,
        call copy_0D(arr,arr_dum) ! it must be identical on all processes
        call broadcast(arr_dum)
      else
        arr_dum = 0.
      endif
      if(present(record)) then
        nrecs8 = 0
        rc = nfmpi_inq_unlimdim(fid,did)
        rc = nfmpi_inq_dimlen(fid,did,nrecs8)
        if(record.le.0 .or. nrecs8+1.lt.record) then
          if(grid%am_i_globalroot) write(6,*)
     &         'error in record dim spec. for variable ',trim(varname)
          call stoprc(0,1)
        endif
        cnt(nd+1) = 1
        srt(nd+1) = record
      endif
      if(grid%am_i_globalroot) then
        if(nd.gt.0) cnt(1:nd) = shape(arr)
        rc = nfmpi_put_vara_double_all(fid,vid,srt,cnt,arr)
      else
        rc = nfmpi_put_vara_double_all(fid,vid,srt,cnt,arr_dum)
      endif
