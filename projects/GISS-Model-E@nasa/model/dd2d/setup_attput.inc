      type(dist_grid) :: grid
      integer :: fid
      character(len=*) :: varname,attname
      integer :: rc,rc2,vid,ierr,attlen,attlenx
      logical :: do_enddef
      if(grid%am_i_globalroot) then
        if(trim(varname).eq.'global') then
          rc = nf_noerr
          vid = nf_global
        else
          rc = nf_inq_varid(fid,trim(varname),vid)
          if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &         trim(varname),' not found in output file - stopping'
        endif
      endif
      call stoprc(rc,nf_noerr)

      if(grid%am_i_globalroot) then
c
c see whether the attribute needs to be defined.
c if the attribute exists but with a different size, delete it.
c if the file is not in define mode, put it into define mode
c
        attlen = len_of_obj(attval)
        do_enddef = .false.
        attlenx = 0
        rc = nf_inq_attlen(fid,vid,trim(attname),attlenx)
        if(rc.eq.nf_enotatt) then
          do_enddef = nf_redef(fid).eq.nf_noerr
        elseif(attlen.ne.attlenx) then
          do_enddef = nf_redef(fid).eq.nf_noerr
          rc = nf_del_att(fid,vid,trim(attname))
        endif
      endif
