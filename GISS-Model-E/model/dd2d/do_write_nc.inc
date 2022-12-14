      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: varname
      integer :: vid,rc
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,trim(varname),vid)
        if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in output file - stopping'
      endif
      call stoprc(rc,nf_noerr)
      if(grid%am_i_globalroot) then
        rc = nf_put_var_double(fid,vid,arr)
      endif
