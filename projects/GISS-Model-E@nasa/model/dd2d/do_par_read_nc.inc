      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer, intent(in), optional :: jdim
      logical, intent(in), optional :: no_xdim
      integer, intent(in), optional :: record,record1
      integer :: jd,nd,jm,vid,unlim_did,rc,rc2,td,rd
      integer, dimension(7) :: srt,cnt,dids
      integer :: nrecs,dlen,itmp,idim,ndf
      logical :: var_has_unlimdim,partial_read
#ifndef SERIAL_MODE
      logical :: snglread
      integer :: i,j,k,l,m,tile
      real*8, dimension(:), allocatable :: arr1d
      real*8, dimension(:,:,:,:,:), allocatable :: arr5d
      integer :: nl,nk,nsend,ierr,dst_pe,tag,modpe
      integer :: mpi_status(mpi_status_size)
      logical :: has_xdim
#endif
      if(present(jdim)) then
        jd=jdim
      else
        jd=2
      endif
      nd = size(shape(arr))
      td = nd+1
      srt(:) = 1
      cnt(:) = 1
      cnt(1:nd) = shape(arr)
      jm = grid%npy
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,trim(varname),vid)
        if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in input file - skipping'
      endif
#ifndef SERIAL_MODE
      call mpi_bcast(rc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      if(rc.ne.nf_noerr) return
      if(grid%am_i_globalroot) then
        rc = nf_inq_varndims(fid,vid,ndf)
        rc = nf_inq_vardimid(fid,vid,dids)
        rc = nf_inq_unlimdim(fid,unlim_did)
        var_has_unlimdim = dids(ndf).eq.unlim_did
        if(var_has_unlimdim) then
          rd = nd ! default
        else
          rd = -1
        endif
      endif
      if(present(record)) then
        if(grid%am_i_globalroot) then
          if(grid%ntiles.gt.1) then
            rd = td + 1
          else
            rd = td
            td = rd + 1 ! to prevent srt(td)= from overwriting srt(rd)
          endif
          nrecs = 0
          rc = nf_inq_dimlen(fid,unlim_did,nrecs)
          if(record.le.0 .or. nrecs.lt.record) then
            write(6,*)
     &           'error in record dim spec. for variable ',trim(varname)
            rc = 0
          else
            rc = 1
          endif
          srt(rd) = record
          cnt(rd) = 1
          ndf = ndf - 1
        endif
        call stoprc(rc,1)
      elseif(present(record1)) then
        if(grid%am_i_globalroot) then
          nrecs = 0
          rc = nf_inq_dimlen(fid,unlim_did,nrecs)
          if(record1.le.0 .or. nrecs.lt.record1) then
            write(6,*)
     &           'error in record dim spec. for variable ',trim(varname)
            rc = 0
          else
            rc = 1
          endif
          srt(rd) = record1
        endif
        call stoprc(rc,1)
      endif
      if(grid%am_i_globalroot) then
        if(grid%ntiles.gt.1 .and. var_has_unlimdim
     &       .and. .not.present(record)
     &       ) then
         ! if the variable has the unlimited dimension, swap the
         ! srt/cnt info for last two dimensions
          itmp = srt(rd); srt(rd) = srt(td); srt(td) = itmp
          itmp = cnt(rd); cnt(rd) = cnt(td); cnt(td) = itmp
          td = nd; rd = td + 1
        endif
        ! shape checking
        partial_read = .false.
        rc2 = 1
        do idim=1,ndf
          if(idim.eq.jd-1 .or. idim.eq.jd .or. idim.eq.td) cycle
          rc = nf_inq_dimlen(fid,dids(idim),dlen)
          if(dlen.ne.cnt(idim)) then
            if(idim.eq.rd) then
              if(dlen.lt.cnt(idim)+srt(idim)-1) then
                write(6,*)
     &          'warning: only performing partial read of variable '//
     &          trim(varname)//' along the unlimited dimension'
                cnt(idim) = dlen-srt(idim)+1
                partial_read = .true.
              endif
            else
              write(6,*) 'shape mismatch for variable ',trim(varname)
              rc2 = 0
              exit
            endif
          endif
        enddo
      endif
      call stoprc(rc2,1)
#ifndef SERIAL_MODE
      if(grid%je.lt.grid%js) return ! empty domain
      tag=1024
      if(present(no_xdim)) then
        has_xdim = .not.no_xdim
      elseif(jd.eq.1) then
        has_xdim = .false.
      else
        has_xdim = .true.
      endif
c if latlon model, call routines better suited for very thin
c zonal strips if they are available
c      if(grid%ntiles.eq.1 .and. has_xdim .and. jd.ne.4 .and.
c     &     (jd.eq.2 .or. jd.lt.nd)) then
c        call par_read_jdecomp_optimized(grid,fid,vid,arr,jd)
c        return
c      endif
      if(grid%am_i_rowroot) then
        call get_nlnk(arr,jd,nl,nk)
        nsend = nl*nk*grid%maxnj
        if(has_xdim) then
          nsend = nsend*grid%npx
        elseif(jd.gt.1) then
          nsend = nsend*size(arr,jd-1)
        endif
        allocate(arr1d(nsend))
        arr1d = 0.
        if(has_xdim) cnt(jd-1) = grid%npx
        cnt(td) = 1
      endif
      if(grid%am_i_globalroot) then
        ! change snglread back to false if memory on root node is an issue
        ! todo: call get_varm to swap tile/record dims if necessary
        snglread = grid%ntiles.eq.1 .or. .not.var_has_unlimdim
        if(snglread) then ! read tne entire variable in 1 shot
          allocate(
     &       arr5d(nl,nsend/(nl*nk*grid%maxnj),grid%npy,nk,grid%ntiles))
          srt(jd) = 1
          cnt(jd) = grid%npy
          srt(td) = 1
          cnt(td) = grid%ntiles
          rc = nf_get_vara_double(fid,vid,srt,cnt,arr5d)
        endif
c loop in reverse order so that arr1d filled with roots row at end of loop
        do dst_pe=grid%nproc-grid%nprocx,0,-grid%nprocx
c read from netcdf file
          modpe=mod(dst_pe,grid%nproc_tile)+1
          srt(jd) = grid%jsr(modpe)
          cnt(jd) = 1+grid%jer(modpe)-grid%jsr(modpe)
          if(cnt(jd).le.0) cycle ! no data to send to this PE
          srt(td) = 1 + dst_pe/grid%nproc_tile
          if(snglread) then
            tile = srt(td)
            m = 0
            do k=1,size(arr5d,4)
            do j=grid%jsr(modpe),grid%jer(modpe)
            do i=1,size(arr5d,2)
            do l=1,size(arr5d,1)
              m = m + 1
              arr1d(m) = arr5d(l,i,j,k,tile)
            enddo
            enddo
            enddo
            enddo
          else
            rc = nf_get_vara_double(fid,vid,srt,cnt,arr1d)
          endif
          if(dst_pe.gt.0) then
            call mpi_send(arr1d,nsend,MPI_DOUBLE_PRECISION,dst_pe,
     &           tag,MPI_COMM_WORLD,ierr)
          endif
        enddo
      elseif(grid%am_i_rowroot) then
        call mpi_recv(arr1d,nsend,MPI_DOUBLE_PRECISION,0,
     &       tag,MPI_COMM_WORLD,mpi_status,ierr)
      endif
      if(has_xdim) then
        call unpack_row(grid,arr1d,arr,jdim=jd)
      else
        call unpack_row_no_xdim(grid,arr1d,arr,jd)
      endif
      if(grid%am_i_rowroot) deallocate(arr1d)
      if(grid%am_i_globalroot .and. snglread) deallocate(arr5d)
#else
      if(size(arr,jd).eq.jm .or. size(arr,jd).eq.1) then
        srt(jd-1) = grid%i_strt
        srt(jd)   = grid%j_strt
        rc = nf_get_vara_double(fid,vid,srt,cnt,arr)
      else
        call stop_model('unrecognized jdim size in do_par_read_nc',255)
      endif
#endif
