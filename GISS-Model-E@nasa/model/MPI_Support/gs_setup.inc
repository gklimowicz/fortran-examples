     &     _
     &     _OPER_
     &     _
     &     _DOMAIN_
#ifdef _GATHER_
     &     (grid,local_arr,global_arr,jdim)
#endif
#ifdef _SCATTER_
     &     (grid,global_arr,local_arr,jdim)
#endif
      type(dist_grid) :: grid
      real*8 :: local_arr(_ARRDIMSL_)

#ifdef MULTITILE
      real*8 :: global_arr(_ARRDIMSL_,:)
#endif
#ifdef TILE
      real*8 :: global_arr(_ARRDIMSL_)
#endif
#ifdef ROW
      real*8 :: global_arr(:)
#endif
      integer, intent(in), optional :: jdim
      integer :: nt,nl,nk,i1,i2,j1,j2,i1g,i2g,j1g,j2g,n,nstart
      integer :: xdim,jd,isize,jsize,irem,jrem,itile
      logical :: has_halo
      integer, dimension(6) :: ! 6 is the max possible ntiles
     &      cntsijg,displsijg,cntslijg,displslijg
      integer, dimension(:), allocatable ::
     &     cntslij,displslij,cntsij,displsij
      logical :: am_i_gsroot
      integer :: comm_gs,nproc_comm
c
c
c

#ifdef MULTITILE
      nt=grid%ntiles
#else
      nt=1
#endif
      i1g=1
      i2g=grid%npx
#ifdef ROW
      j1g=grid%js
      j2g=grid%je
      am_i_gsroot = grid%am_i_rowroot
      comm_gs = grid%comm_row
      nproc_comm = grid%nprocx
#else
      j1g=1
      j2g=grid%npy
      am_i_gsroot = grid%am_i_tileroot
      comm_gs = grid%comm_tile
      nproc_comm = grid%nproc_tile
#endif

      if(present(jdim)) then
        jd=jdim
      else
        jd=2
      endif

c
c see whether the halo is present
c
      isize = size(local_arr,jd-1)
      irem = isize - (grid%ie-grid%is+1)
      jsize = size(local_arr,jd)
      jrem = jsize - (grid%je-grid%js+1)
      if(mod(irem,2).eq.0) then ! defined at cell centers
        i1 = grid%is - irem/2; i2 = grid%ie + irem/2
      else                      ! defined at cell edges
        i1 = grid%is - (irem-1)/2; i2 = grid%ie + 1 + (irem-1)/2
      endif
      if(mod(jrem,2).eq.0) then ! defined at cell centers
        j1 = grid%js - jrem/2; j2 = grid%je + jrem/2
      else                      ! defined at cell edges
        j1 = grid%js - (jrem-1)/2; j2 = grid%je + 1 + (jrem-1)/2
      endif
      if(irem.eq.0 .and. jrem.eq.0) then
        has_halo=.false.
      else
        has_halo=.true.
      endif

      call get_nlnk(local_arr,jd,nl,nk)
#ifndef SERIAL_MODE
      if(nproc_comm.gt.1 .or. nt.gt.1) then
        call alloc_gs_wksp(grid,nl,nk,1+j2g-j1g,nt,am_i_gsroot)
      else
        call alloc_gs_wksp_default()
      end if

      if(nproc_comm.gt.1 .or. nt.gt.1) then
c calculate gatherv/scatterv info
        if(nt.gt.1) then
          cntsijg(:) = grid%npx*grid%npy
          do itile=1,grid%ntiles
            displsijg(itile) = grid%npx*grid%npy*(itile-1)
          enddo
        endif
        if(nl.gt.1) then
          if(am_i_gsroot) then
            allocate(cntslij(grid%nproc_tile)
     &           ,displslij(grid%nproc_tile))
            do n=1,nproc_comm
              cntslij(n) = grid%cntsij(n+grid%rank_tile)*nl
            enddo
            displslij(1) = 0
            do n=2,nproc_comm
              displslij(n) = displslij(n-1) + cntslij(n-1)
            enddo
          endif
          if(nt.gt.1) then
            cntslijg(:) = cntsijg(:)*nl
            displslijg(:) = displsijg(:)*nl*nk
          endif
        else
          if(am_i_gsroot) then
            allocate(cntsij(grid%nproc_tile)
     &           ,displsij(grid%nproc_tile))
            do n=1,nproc_comm
              cntsij(n) = grid%cntsij(n+grid%rank_tile)
            enddo
            displsij(1) = 0
            do n=2,nproc_comm
              displsij(n) = displsij(n-1) + cntsij(n-1)
            enddo
          endif
        endif
      endif
#else
      allocate(cntsij(1), displsij(1)) ! unused but passed to routine with fixed size arrays
      call alloc_gs_wksp_default()
#endif
      if (.not. allocated(cntsij)) allocate(cntsij(1))
      if (.not. allocated(displsij)) allocate(displsij(1))
        
#ifdef _GATHER_
c
c gather
c
      if(nl.eq.1) then ! for efficiency
        call gather3D(grid,
     &        local_arr,global_arr
     &       ,i1,i2,j1,j2,nl,nk,nt,has_halo,nkmax
     &       ,i1g,i2g,j1g,j2g,am_i_gsroot,comm_gs,nproc_comm
     &       ,cntsij,displsij,cntsijg,displsijg
     &       ,buf1d_local_size,buf1d_tile_size,bufij_tile_size
     &       ,buf1d_local,buf1d_tile,bufij_tile
     &       )
      else
        call gather4D(grid,
     &        local_arr,global_arr
     &       ,i1,i2,j1,j2,nl,nk,nt,has_halo
     &       ,i1g,i2g,j1g,j2g,am_i_gsroot,comm_gs,nproc_comm
     &       ,cntslij,cntslijg,displslij,displslijg
     &       ,buf1d_local_size,buf1d_tile_size,bufij_tile_size
     &       ,buf1d_local,buf1d_tile,bufij_tile
     &       )
      endif

#else
c
c scatter
c
      if(nl.eq.1) then ! for efficiency
        call scatter3D(grid,
     &        local_arr,global_arr
     &       ,i1,i2,j1,j2,nl,nk,nt,has_halo,nkmax
     &       ,i1g,i2g,j1g,j2g,am_i_gsroot,comm_gs,nproc_comm
     &       ,cntsij,displsij,cntsijg,displsijg
     &       ,buf1d_local_size,buf1d_tile_size,bufij_tile_size
     &       ,buf1d_local,buf1d_tile,bufij_tile
     &       )
      else
        call scatter4D(grid,
     &        local_arr,global_arr
     &       ,i1,i2,j1,j2,nl,nk,nt,has_halo
     &       ,i1g,i2g,j1g,j2g,am_i_gsroot,comm_gs,nproc_comm
     &       ,cntslij,cntslijg,displslij,displslijg
     &       ,buf1d_local_size,buf1d_tile_size,bufij_tile_size
     &       ,buf1d_local,buf1d_tile,bufij_tile
     &       )
      endif
#endif

#ifndef SERIAL_MODE
      if(nproc_comm.gt.1) then
        if(am_i_gsroot) then
          if(nl.gt.1) then
            deallocate(cntslij,displslij)
          else
            deallocate(cntsij,displsij)
          endif
        else
          deallocate(cntsij,displsij)
        endif
      endif
#endif
