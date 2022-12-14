      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical, intent(in), optional :: bcast_all
      integer :: vid,rc,ierr,nd,arrsize,shp(7)
      integer*8, dimension(7) :: srt,cnt
      real*8 :: arr_dum
      logical :: bcast_all_
      bcast_all_ = .false.
      if(present(bcast_all)) bcast_all_ = bcast_all
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(rc.ne.nf_noerr) then
        if(grid%am_i_globalroot) write(6,*) 'variable ',
     &       trim(varname),' not found in input file - skipping'
        return
      endif
      if(bcast_all_) then
        rc = nfmpi_get_var_double_all(fid,vid,arr)
      else
c        if(grid%am_i_globalroot) then
c          rc = nfmpi_get_var_double(fid,vid,arr)
c        else
c        endif
        nd = size(shape(arr))
        srt(:) = 1
        cnt(:) = 0
        if(grid%am_i_globalroot) then
          if(nd.gt.0) cnt(1:nd) = shape(arr)
          rc = nfmpi_get_vara_double_all(fid,vid,srt,cnt,arr)
        else
          rc = nfmpi_get_vara_double_all(fid,vid,srt,cnt,arr_dum)
        endif
      endif
