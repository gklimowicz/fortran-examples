      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer, intent(in), optional :: jdim
      integer, intent(in), optional :: record,record1
      logical, intent(in), optional :: no_xdim
#ifndef SERIAL_MODE
      real*8, dimension(:), allocatable :: arr_nohalo
      integer :: nl,nk,jd,nd,vid,unlim_did,rc,arrsize,td,rd,ndf
      integer :: itmp,idim
      integer*8, dimension(7) :: srt,cnt
      integer, dimension(7) :: dids
      integer*8 :: nrecs8,dlen
      logical :: var_has_unlimdim,partial_read
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(rc.ne.nf_noerr) then
        if(grid%am_i_globalroot) write(6,*) 'variable ',
     &       trim(varname),' not found in input file - skipping'
	return
      endif
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
      arrsize = product(cnt(1:nd))
      rc = nfmpi_inq_varndims(fid,vid,ndf)
      rc = nfmpi_inq_vardimid(fid,vid,dids)
      rc = nfmpi_inq_unlimdim(fid,unlim_did)
      var_has_unlimdim = dids(ndf).eq.unlim_did
      if(var_has_unlimdim) then
        rd = nd ! default
      else
        rd = -1
      endif
      if(present(record)) then
        nrecs8 = 0
        rc = nfmpi_inq_dimlen(fid,unlim_did,nrecs8)
        if(record.le.0 .or. nrecs8.lt.record) then
          if(grid%am_i_globalroot) write(6,*)
     &         'error in record dim spec. for variable ',trim(varname)
          call stoprc(0,1)
        endif
        rd = td
        if(grid%ntiles.gt.1) rd = rd + 1
        srt(rd) = record
        cnt(rd) = 1
        ndf = ndf - 1
      elseif(present(record1)) then
        nrecs8 = 0
        rc = nfmpi_inq_dimlen(fid,unlim_did,nrecs8)
        if(record1.le.0 .or. nrecs8.lt.record1) then
          if(grid%am_i_globalroot) write(6,*)
     &         'error in record dim spec. for variable ',trim(varname)
          call stoprc(0,1)
        endif
        if(grid%ntiles.gt.1) then
          ! swap the srt/cnt info for last two dimensions
          itmp = srt(rd); srt(rd) = srt(td); srt(td) = itmp
          itmp = cnt(rd); cnt(rd) = cnt(td); cnt(td) = itmp
          td = nd; rd = td + 1
        endif
        srt(rd) = record1
      elseif(grid%ntiles.gt.1 .and. var_has_unlimdim) then
      ! if the variable has the unlimited dimension, swap the
      ! srt/cnt info for last two dimensions
        itmp = srt(rd); srt(rd) = srt(td); srt(td) = itmp
        itmp = cnt(rd); cnt(rd) = cnt(td); cnt(td) = itmp
        td = nd; rd = td + 1
      endif
      ! shape checking
      partial_read = .false.
      do idim=1,ndf
        if(idim.eq.jd-1 .or. idim.eq.jd .or. idim.eq.td) cycle
        rc = nfmpi_inq_dimlen(fid,dids(idim),dlen)
        if(dlen.ne.cnt(idim)) then
          if(idim.eq.rd) then
            if(dlen.lt.cnt(idim)+srt(idim)-1) then
              if(grid%am_i_globalroot)
     &          write(6,*)
     &          'warning: only performing partial read of variable '//
     &          trim(varname)//' along the unlimited dimension'
              cnt(idim) = dlen-srt(idim)+1
              partial_read = .true.
            endif
          else
            if(grid%am_i_globalroot) then
              write(6,*) 'shape mismatch for variable ',trim(varname)
            endif
            call stoprc(0,1)
          endif
        endif
      enddo
      if(grid%have_domain) then
        call get_nlnk(arr,jd,nl,nk)
      else
        srt=1; cnt=0   ! empty domain
      endif
      if(size(arr,jd).ne.cnt(jd)) then ! halo present
        allocate(arr_nohalo(arrsize))
        if(partial_read) arr_nohalo = 0.
        rc = nfmpi_get_vara_double_all(fid,vid,srt,cnt,arr_nohalo)
        if(arrsize.gt.0) then
          call copy_from_1D( arr_nohalo,arr, nl,nk,
     &         grid%isd,grid%ied,grid%jsd,grid%jed,
     &         grid%is,grid%ie,grid%js,grid%je )
        endif
        deallocate(arr_nohalo)
      else
        rc = nfmpi_get_vara_double_all(fid,vid,srt,cnt,arr)
      endif
#else
c call read_data when in serial mode
      call read_data(grid,fid,varname,arr)
#endif
