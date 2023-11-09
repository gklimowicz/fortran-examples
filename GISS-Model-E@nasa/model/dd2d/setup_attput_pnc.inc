      type(dist_grid) :: grid
      integer :: fid
      character(len=*) :: varname,attname
      integer :: rc,rc2,vid,ierr,attlen,attlenx
      integer*8 :: attlen8
      logical :: do_enddef
      logical :: am_root
      am_root = grid%am_i_globalroot
      if(trim(varname).eq.'global') then
        vid = nf_global
      else
        rc = nfmpi_inq_varid(fid,trim(varname),vid)
        if(am_root .and. rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in output file - stopping'
        call stoprc(rc,nf_noerr)
      endif

c
c see whether the attribute needs to be defined.
c if the attribute exists but with a different size, delete it.
c if the file is not in define mode, put it into define mode
c
      attlen = len_of_obj(attval)
      do_enddef = .false.
      attlen8 = 0
      rc = nfmpi_inq_attlen(fid,vid,trim(attname),attlen8)
      attlenx = attlen8
      if(rc.eq.nf_enotatt) then
        do_enddef = nfmpi_redef(fid).eq.nf_noerr
      elseif(attlen.ne.attlenx) then
        do_enddef = nfmpi_redef(fid).eq.nf_noerr
        rc = nfmpi_del_att(fid,vid,trim(attname))
      endif
      attlen8 = attlen
