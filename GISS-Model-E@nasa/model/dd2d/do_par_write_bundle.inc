      integer :: fid
      character(len=*) :: varnames(:)
      type(dist_grid), intent(in) :: grid
      integer, intent(in), optional :: jdim
      integer, intent(in), optional :: record
#ifndef SERIAL_MODE
      real*8, dimension(:), allocatable :: arr_nohalo
      integer :: nl,nk,jd,nd,n,nb,ierr,vid,did,rc,td,rd,vid_last
      integer, dimension(:), allocatable :: write_reqs,statuses
      integer*8, dimension(7) :: srt,cnt
      integer*8 :: nrecs8,m,arrsize
      integer :: isd,ied,jsd,jed
      nd = size(shape(arr))
      cnt(1:nd) = shape(arr)
      nb = cnt(nd)
      if(nb.ne.size(varnames)) then
        if(grid%am_i_globalroot) write(6,*)
     &       'do_par_write_bundle: size(varnames)!=size(arr,ndims)'
        call stoprc(0,1)
      endif
      if(present(jdim)) then
        jd=jdim
      else
        jd=2
      endif
      nd = nd-1
      td = nd+1
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
     &         'error in record dim spec. for bundle ',varnames
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
      arrsize = nb*product(cnt(1:nd))
      allocate(arr_nohalo(arrsize))
      if(arrsize.gt.0) then
        if(size(arr,jd).ne.cnt(jd)) then ! halo present
          isd=grid%isd; ied=grid%ied; jsd=grid%jsd; jed=grid%jed
        else
          isd=grid%is; ied=grid%ie; jsd=grid%js; jed=grid%je
        endif
        call copy_to_1D( arr,arr_nohalo,nl,nk,
     &       isd,ied,jsd,jed,
     &       grid%is,grid%ie,grid%js,grid%je )
      endif
      allocate(write_reqs(nb),statuses(nb))
      m = 1
      do n=1,nb
        rc = nfmpi_inq_varid(fid,trim(varnames(n)),vid)
        if(grid%am_i_globalroot .and. rc.ne.nf_noerr)
     &       write(6,*) 'variable ',
     &       trim(varnames(n)),' not found in output file - stopping'
        call stoprc(rc,nf_noerr)
        if(n.gt.1 .and. vid-vid_last.ne.1) then
          ! Some MPI-IO libraries currently bungle bundled
          ! writes of noncontiguous subsets of variables
          if(grid%am_i_globalroot) then
            write(6,*) 'do_par_write_bundle: a nonsequential collection'
            write(6,*) 'of IDs is possible but currently unsafe.'
            write(6,*) 'Define this bundle w/ consecutive variable IDs:'
            write(6,*) varnames
          endif
          call stoprc(0,1)
        endif
        vid_last = vid
c        rc = nfmpi_put_vara_double_all(fid,vid,srt,cnt,arr_nohalo(m))
        rc = nfmpi_iput_vara_double(fid,vid,srt,cnt,arr_nohalo(m),
     &       write_reqs(n))
        m = m + arrsize/nb
      enddo
      rc = nfmpi_wait_all(fid,nb,write_reqs,statuses)
      deallocate(arr_nohalo)
      deallocate(write_reqs,statuses)
#else
c call write_data when in serial mode
      call write_data(grid,fid,varname,arr)
#endif
