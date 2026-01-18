      type(dist_grid) :: grid
      integer :: fid
      character(len=*) :: varname,attname
      integer :: attlen
      integer, intent(in), optional :: attnum
      integer :: rc,vid,ierr,l,attlenx
      integer*8 :: attlen8
      logical :: am_root
      am_root = grid%am_i_globalroot
      if(trim(varname).eq.'global') then
        vid = nf_global
      else
        rc = nfmpi_inq_varid(fid,trim(varname),vid)
        if(am_root .and. rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in input file - skipping'
        if(rc.ne.nf_noerr) return
      endif
c
c Get the length of this attribute, and its name if it was
c requested by number
c
      if(present(attnum)) then
        rc = nfmpi_inq_attname(fid,vid,attnum,attname)
        rc = nfmpi_inq_attlen(fid,vid,trim(attname),attlen8)
        attlen = attlen8
      else
        rc = nfmpi_inq_attlen(fid,vid,trim(attname),attlen8)
        attlen = attlen8
        if(am_root .and. rc.ne.nf_noerr) then
          write(6,*) 'In read_attr:'
          write(6,*) 'attribute ',trim(attname),
     &         ' not found for variable ',trim(varname)
        endif
        if(rc.ne.nf_noerr) return
        !call stoprc(rc,nf_noerr)
      endif
      attlenx = full_len_of_obj(attval)
      if(attlen.gt.attlenx) then
        if(am_root) then
          write(6,*) 'In read_attr:'
          write(6,*) 'For ',trim(varname)//':'//trim(attname)
          write(6,*) 'Attribute length on disk is: ',attlen
          write(6,*) 'Data length in the program is: ',attlenx
        endif
        rc = 1
      else
        rc = 0
      endif
      call stoprc(rc,0)
