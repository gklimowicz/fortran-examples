      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical, intent(in), optional :: bcast_all
      integer :: vid,rc,ierr,nd,arrsize,shp(7)
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,trim(varname),vid)
        if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in input file - skipping'
      endif
#ifndef SERIAL_MODE
      call mpi_bcast(rc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      if(rc.ne.nf_noerr) return
      if(grid%am_i_globalroot) rc = nf_get_var_double(fid,vid,arr)
#ifndef SERIAL_MODE
      if(present(bcast_all)) then
        if(bcast_all) then
          nd = size(shape(arr))
          if(nd.eq.0) then
            arrsize = 1
          else
            shp(1:nd) = shape(arr)
            arrsize = product(shp(1:nd))
          endif
          call mpi_bcast(arr,arrsize,
     &         MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif
      endif
#endif
