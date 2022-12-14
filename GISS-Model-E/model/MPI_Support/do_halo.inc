      type(dist_grid) :: grid
      integer, intent(in), optional :: jdim
      integer :: nl,nk
      integer :: xdim,jd,np,pep1,pem1,isize,jsize
     &     ,hwidth,pe_send,pe_recv
      integer :: i1,i2,j1,j2
      integer ::
     &     i1np,i2np,j1np,j2np,i1sr,i2sr,j1sr,j2sr
     &    ,i1sp,i2sp,j1sp,j2sp,i1nr,i2nr,j1nr,j2nr
     &    ,ip,ir,jp,jr
      if(present(jdim)) then
        jd = jdim
      else
        jd = 2
      endif

      call get_nlnk(arr,jd,nl,nk)
      jsize = size(arr,jd)

      hwidth = (jsize-(1+grid%je-grid%js))/2

      i1 = grid%is
      i2 = grid%ie
      j1 = grid%js-hwidth
      j2 = grid%je+hwidth

      i1np=grid%i1np
      i2np=grid%i2np
      j1np=grid%j1np
      j2np=grid%j2np
      i1sr=grid%i1sr
      i2sr=grid%i2sr
      j1sr=grid%j1sr
      j2sr=grid%j2sr
      i1sp=grid%i1sp
      i2sp=grid%i2sp
      j1sp=grid%j1sp
      j2sp=grid%j2sp
      i1nr=grid%i1nr
      i2nr=grid%i2nr
      j1nr=grid%j1nr
      j2nr=grid%j2nr

      if(grid%ntiles.eq.6) then

      i1 = i1-hwidth
      i2 = i2+hwidth
c
c east-west
c
      np = grid%nproc_comm_ew
      pep1 = mod(grid%rank_ew+1   ,np)
      pem1 = mod(grid%rank_ew-1+np,np)
c send to pep1, receive from pem1
      call sendrecv4D(arr,nl,nk,
     &     i1,i2,j1,j2,
     &     grid%i1ep,grid%i2ep,grid%j1ep,grid%j2ep,
     &     1,1,
     &     grid%i1wr,grid%i2wr,grid%j1wr,grid%j2wr,
     &     grid%comm_ew,pep1,pem1,
     &     bufsend,bufrecv
     &     )
c send to pem1, receive from pep1
      call sendrecv4D(arr,nl,nk,
     &     i1,i2,j1,j2,
     &     grid%i1wp,grid%i2wp,grid%j1wp,grid%j2wp,
     &     1,1,
     &     grid%i1er,grid%i2er,grid%j1er,grid%j2er,
     &     grid%comm_ew,pem1,pep1,
     &     bufsend,bufrecv
     &     )

c
c modify north-south bounds to get most corner halos "for free"
c
      if(grid%nprocx.gt.1) then
      if(j1np.eq.j2np) then
        i1np=i1np-grid%iincnp   !*hwidth
        i2np=i2np+grid%iincnp   !*hwidth
      else
        j1np=j1np-grid%jincnp   !*hwidth
        j2np=j2np+grid%jincnp   !*hwidth
      endif
      if(j1sr.eq.j2sr) then
        i1sr=i1sr-1             !*hwidth
        i2sr=i2sr+1             !*hwidth
      else
        j1sr=j1sr-1             !*hwidth
        j2sr=j2sr+1             !*hwidth
      endif
      if(j1sp.eq.j2sp) then
        i1sp=i1sp-grid%iincsp   !*hwidth
        i2sp=i2sp+grid%iincsp   !*hwidth
      else
        j1sp=j1sp-grid%jincsp   !*hwidth
        j2sp=j2sp+grid%jincsp   !*hwidth
      endif
      if(j1nr.eq.j2nr) then
        i1nr=i1nr-1             !*hwidth
        i2nr=i2nr+1             !*hwidth
      else
        j1nr=j1nr-1             !*hwidth
        j2nr=j2nr+1             !*hwidth
      endif
      endif

      endif ! if ntiles.eq.6

c
c north-south
c

      pep1 = grid%rank_ns+1
      pem1 = grid%rank_ns-1
      np = grid%nproc_comm_ns
      if(grid%ntiles.eq.1) then ! lat-lon
        if(grid%gid.eq.np-1) pep1 = MPI_PROC_NULL
        if(grid%gid.eq.0)    pem1 = MPI_PROC_NULL
      else ! periodic
        pep1 = mod(pep1   ,np)
        pem1 = mod(pem1+np,np)
      endif
c send to pep1, receive from pem1
      call sendrecv4D(arr,nl,nk,
     &     i1,i2,j1,j2,
     &     i1np,i2np,j1np,j2np,
     &     grid%iincnp,grid%jincnp,
     &     i1sr,i2sr,j1sr,j2sr,
     &     grid%comm_ns,pep1,pem1,
     &     bufsend,bufrecv
     &     )
c send to pem1, receive from pep1
      call sendrecv4D(arr,nl,nk,
     &     i1,i2,j1,j2,
     &     i1sp,i2sp,j1sp,j2sp,
     &     grid%iincsp,grid%jincsp,
     &     i1nr,i2nr,j1nr,j2nr,
     &     grid%comm_ns,pem1,pep1,
     &     bufsend,bufrecv
     &     )

      if(grid%nprocx.gt.1) then
c
c for corner halo:
c extra sendrecv in the sw-ne direction if next to the diagnonal
c
      pe_send = MPI_PROC_NULL
      if(grid%pe_send_sw.ne.MPI_PROC_NULL) then
        pe_send = grid%pe_send_sw
        ip=grid%is; jp=grid%js; ir=grid%is-1; jr=grid%js-1 ! fix later for hwidth>1
      elseif(grid%pe_send_ne.ne.MPI_PROC_NULL) then
        pe_send = grid%pe_send_ne
        ip=grid%ie; jp=grid%je; ir=grid%ie+1; jr=grid%je+1 ! fix later for hwidth>1
      endif
      if(pe_send.ne.MPI_PROC_NULL) then
        call sendrecv4D(arr,nl,nk,
     &       i1,i2,j1,j2,
     &       ip,ip,jp,jp,
     &       1,1,
     &       ir,ir,jr,jr,
     &       MPI_COMM_WORLD,pe_send,pe_send, ! send pe = recv pe
     &       bufsend,bufrecv
     &       )
      endif

c
c for corner halo:
c extra unidirectional send in the se-nw direction if on the diagnonal
c
      ip=grid%is; jp=grid%js    ! defaults
      pe_send = grid%pe_diag_s
      if(pe_send.ne.MPI_PROC_NULL) then
        if(pe_send.gt.grid%gid) then ! send to the NW
          ip=grid%is; jp=grid%je     ! fix later for hwidth>1
        else                         ! send to the SE
          ip=grid%ie; jp=grid%js     ! fix later for hwidth>1
        endif
      endif
      pe_recv = grid%pe_diag_r
      if(pe_recv.ne.MPI_PROC_NULL) then
        if(pe_recv.gt.grid%gid) then ! recv from the NW
          ir=grid%is-1; jr=grid%je+1 ! fix later for hwidth>1
        else                         ! recv from the SE
          ir=grid%ie+1; jr=grid%js-1 ! fix later for hwidth>1
        endif
      endif
      if(pe_send.ne.MPI_PROC_NULL.or.pe_recv.ne.MPI_PROC_NULL) then
        call sendrecv4D(arr,nl,nk,
     &       i1,i2,j1,j2,
     &       ip,ip,jp,jp,
     &       1,1,
     &       ir,ir,jr,jr,
     &       MPI_COMM_WORLD,pe_send,pe_recv,
     &       bufsend,bufrecv
     &       )
      endif

      endif
