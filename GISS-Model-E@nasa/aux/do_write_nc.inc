      integer :: fid
      character(len=*) :: varname
      integer :: vid,rc
      rc = nf_inq_varid(fid,trim(varname),vid)
      if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &     trim(varname),' not found in output file - stopping'
      rc = nf_put_var_real(fid,vid,arr)
