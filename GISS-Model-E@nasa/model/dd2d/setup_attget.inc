      type(dist_grid) :: grid
      integer :: fid
      character(len=*) :: varname,attname
      integer :: attlen
      integer, intent(in), optional :: attnum
      integer :: rc,vid,ierr,lname,l,attlenx
      character :: tmpstr(1500)
      if(grid%am_i_globalroot) then
        if(trim(varname).eq.'global') then
          rc = nf_noerr
          vid = nf_global
        else
          rc = nf_inq_varid(fid,trim(varname),vid)
          if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &         trim(varname),' not found in input file - skipping'
        endif
      endif
      call broadcast(rc)
      if(rc.ne.nf_noerr) return
      if(.not.present(attnum)) then
        if(grid%am_i_globalroot) then
          rc = nf_inq_attlen(fid,vid,trim(attname),attlen)
          if(rc.ne.nf_noerr) then
            write(6,*) 'In read_attr:'
            write(6,*) 'attribute ',trim(attname),
     &           ' not found for variable ',trim(varname)
          endif
        endif
        call broadcast(rc)
        if(rc.ne.nf_noerr) return
        !call stoprc(rc,nf_noerr)
      endif

c
c Get the name and length of this attribute and broadcast them.
c For some reason, mpi_bcast does not work on strings declared
c as character(len=xxx), so we have to copy attname into/out of
c a conventional array
c
      if(grid%am_i_globalroot) then
        if(present(attnum)) then
          rc = nf_inq_attname(fid,vid,attnum,attname)
          lname = len_trim(attname)
          do l=1,lname
            tmpstr(l) = attname(l:l)
          enddo
          rc = nf_inq_attlen(fid,vid,trim(attname),attlen)
        endif
        attlenx = full_len_of_obj(attval)
        if(attlen.gt.attlenx) then
          write(6,*) 'In read_attr:'
          write(6,*) 'For ',trim(varname)//':'//trim(attname)
          write(6,*) 'Attribute length on disk is: ',attlen
          write(6,*) 'Data length in the program is: ',attlenx
          rc = 1
        else
          rc = 0
        endif
      endif
      call stoprc(rc,0)
#ifndef SERIAL_MODE
      call mpi_bcast(attlen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(present(attnum)) then
        call mpi_bcast(lname,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(tmpstr,lname,MPI_CHARACTER,0,
     &       MPI_COMM_WORLD,ierr)
        attname=''
        do l=1,lname
          attname(l:l) = tmpstr(l)
        enddo
      endif
#endif
