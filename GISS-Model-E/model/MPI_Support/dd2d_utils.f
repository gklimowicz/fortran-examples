      module dd2d_utils
!@auth M. Kelley
!@ver 1.0
!@sum  dd2d_utils provides communication procedures for domain
!@+    decomposition using MPI.  Currently, 2D cubed-sphere
!@+    and 1D lat-lon layouts are enabled.  These routines
!@+    were developed to meet the requirements of modelE
!@+    not met by other packages.
!@+
!@usage Access is through the following types/interfaces.
!@+
!@+     Below, grid is an instance of the dist_grid derived type,
!@+     and the optional argument jdim indicates which array index
!@+     corresponds to the last distributed dimension (usually, j).
!@+     If jdim is not specified, it is assumed to be 2.
!@+
!@+     dist_grid: a derived type containing domain decomp info,
!@+     initialized by calling subroutine
!@+     init_dist_grid(
!@+     &     npx,npy,ntiles, ! x,y size of tile and # of tiles
!@+     &     is,ie,js,je,    ! start, end of local i,j domain
!@+     &     isd,ied,jsd,jed, ! start, end of local data domain (not used)
!@+     &     grid)
!@+     ntiles = 1 for a latlon layout, and 6 for a cubed sphere layout
!@+     NOTE: all arguments save for grid are inputs - this routine does
!@+           NOT decide the domain decomposition is,ie,js,je!
!@+
!@+     halo_update(grid,array,jdim)
!@+        fills the halo cells of array with the appropriate values
!@+        from neighboring PEs.  The depth of the halo is inferred
!@+        by comparing the dimension sizes of array to the start,end
!@+        info in grid.  For a latlon grid, no halos are performed in
!@+        the I direction.
!@+        For the cubed-sphere layout, the current implementation
!@+        requires npx=npy and an N x N layout of PEs on each face.
!@+
!@+     pack_data(grid,local_array,global_array,jdim)
!@+        gathers local arrays into a global array on root
!@+     unpack_data(grid,global_array,local_array,jdim)
!@+        scatters a global array on root into local arrays
!@+     If global_array has the same rank as local_array, the
!@+        gather/scatter operations are over the local tile.
!@+        Each tile has a root PE with a valid copy of global_array.
!@+     If the rank of global_array = 1 + rank of local_array, the
!@+        gather/scatter operations are over all tiles.  The last
!@+        index of global_array is the tile index in this case.
!@+
!@+     globalsum(grid,local_arr,arrsum,all)
!@+        calculates the global sum of distributed array local_arr
!@+        and stores the result in arrsum; the result is
!@+        independent of the number of PEs.  If all=.true., all
!@+        PEs receive the result. Currently only implemented for
!@+        real*8 local_arr(i,j).
!@+
!@+     pack_row(grid,local_array,row_array_1d,jdim)
!@+     unpack_row(grid,row_array_1d,local_array,jdim)
!@+        These routines are like pack_data/unpack_data, but only
!@+        gather/scatter to/from one PE for each "row" of PEs.
!@+        For now, the row array with the contents of all
!@+        PEs in the row is 1D since it is only used for I/O.
!@+        Dimensioned versions of row_array will be added soon.
!@+
!@+

c
c Implementation notes for halo_update
c
c For simplicity, halo_update is programmed using mpi_sendrecv
c and predefined MPI communicators for "east-west" and
c "north-south" directions. For the latlon grid, "north-south"
c corresponds to geographic N-S.
c For the cubed sphere 2D decomposition, these communicators
c trace periodic paths.  Given the following example 2x2 PE
c layout on each face,
c
c              18 19 22 23
c              16 17 20 21
c        10 11 14 15
c         8  9 12 13
c   2  3  6  7
c   0  1  4  5
c
c  the PE groupings for the first periodic direction are
c
c  0,2,3,6,8,10,11,14,16,18,19,22
c  1,4,5,7,9,12,13,15,17,20,21,23
c
c  and for the second periodic direction they are
c
c  2,10,18
c  0,1,3,8,9,11,16,17,19
c  4,6,7,12,14,15,20,22,23
c  5,13,21
c
c Most corner halos can be obtained by adjusting halo bounds
c during the sendrecv along the second periodic direction.
c However, for PEs on or next to the UL-LR diagonal of each cubed
c sphere tile, one extra call to mpi_sendrecv is needed.
c

#ifdef OFFLINE
#else
c see whether model E is running in serial mode
#ifndef USE_MPI
#define SERIAL_MODE
#endif
#endif

c#ifdef USE_ESMF
c      use ESMF
c#endif
      use Hidden_mod, only: Hidden_type

      implicit none
      private
#ifndef SERIAL_MODE
      include 'mpif.h'
#endif
      save

      public :: init_dist_grid
      public :: dist_grid
      type dist_grid
        type (Hidden_type) :: private
        integer :: npx ! number of i cells
        integer :: npy ! number of j cells
        integer :: is  ! first i index of computational domain
        integer :: ie  ! last  i index of computational domain
        integer :: js  ! first j index of computational domain
        integer :: je  ! last  j index of computational domain
        integer :: isd ! first i index of data domain
        integer :: ied ! last  i index of data domain
        integer :: jsd ! first j index of data domain
        integer :: jed ! last  j index of data domain
        logical :: am_i_globalroot ! am I the root for all processes
        integer :: ntiles ! number of tiles globally

        integer :: gid ! rank of my processor
        integer :: nproc ! number of processors
        integer :: nproc_tile ! number of processors on my tile
        integer :: nprocx ! number of processors in the x direction
        integer :: nprocy ! number of processors in the y direction
        integer :: tile ! my tile ID
        integer :: rank_tile ! MPI rank on my tile
        integer :: comm_tile ! MPI communicator for my tile
        integer :: comm_intertile ! MPI communicator between tiles
        logical :: am_i_tileroot ! am I the root for my tile
        integer :: root_mytile ! who is root for my tile
        integer, dimension(:), allocatable :: isr ! first i of each PE on my tile
        integer, dimension(:), allocatable :: ier ! last i of each PE on my tile
        integer, dimension(:), allocatable :: jsr ! first j of each PE on my tile
        integer, dimension(:), allocatable :: jer ! last j of each PE on my tile
        integer, dimension(:), allocatable :: cntsij ! # cells for each PE on my tile
        integer :: maxnj  ! maximum value of 1+je-js on my tile

        integer :: rank_row ! MPI rank on my row
        integer :: comm_row ! MPI communicator for my row
        logical :: am_i_rowroot ! am I the root for my row

        integer :: comm_ew ! MPI communicator for east-west halos
        integer :: comm_ns ! MPI communicator for north-south halos
        integer :: rank_ew ! MPI rank for east-west halos
        integer :: rank_ns ! MPI rank for north-south halos
        integer :: nproc_comm_ew ! number of processors in comm_ew
        integer :: nproc_comm_ns ! number of processors in comm_ns

c indices for packing halo messages when halo width is 1
        integer :: i1wp ! first i of the data to be sent west
        integer :: i2wp ! last i of the data to be sent west
        integer :: j1wp ! first j of the data to be sent west
        integer :: j2wp ! last j of the data to be sent west
        integer :: i1ep ! first i of the data to be sent east
        integer :: i2ep ! last i of the data to be sent east
        integer :: j1ep ! first j of the data to be sent east
        integer :: j2ep ! last j of the data to be sent east
        integer :: i1sp ! first i of the data to be sent south
        integer :: i2sp ! last i of the data to be sent south
        integer :: j1sp ! first j of the data to be sent south
        integer :: j2sp ! last j of the data to be sent south
        integer :: i1np ! first i of the data to be sent north
        integer :: i2np ! last i of the data to be sent north
        integer :: j1np ! first j of the data to be sent north
        integer :: j2np ! last j of the data to be sent north
        integer :: iincsp ! i increment for south pack
        integer :: iincnp ! i increment for north pack
        integer :: jincsp ! j increment for south pack
        integer :: jincnp ! j increment for north pack

c indices for receiving halo messages when halo width is 1
        integer :: i1wr ! first i of the data to be received from west
        integer :: i2wr ! last i of the data to be received from west
        integer :: j1wr ! first j of the data to be received from west
        integer :: j2wr ! last j of the data to be received from west
        integer :: i1er ! first i of the data to be received from east
        integer :: i2er ! last i of the data to be received from east
        integer :: j1er ! first j of the data to be received from east
        integer :: j2er ! last j of the data to be received from east
        integer :: i1sr ! first i of the data to be received from south
        integer :: i2sr ! last i of the data to be received from south
        integer :: j1sr ! first j of the data to be received from south
        integer :: j2sr ! last j of the data to be received from south
        integer :: i1nr ! first i of the data to be received from north
        integer :: i2nr ! last i of the data to be received from north
        integer :: j1nr ! first j of the data to be received from north
        integer :: j2nr ! last j of the data to be received from north

c special processors for a cubed sphere layout
        integer :: pe_diag_s ! PE for extra diagonal send
        integer :: pe_diag_r ! PE for extra diagonal receive
        integer :: pe_send_ne ! PE for extra send to the NE
        integer :: pe_send_sw ! PE for extra send to the SW

#ifndef SERIAL_MODE
#endif /* not SERIAL_MODE */

c
c the following is for use in modelE only:
c
        logical :: have_domain ! whether this PE has any of the domain
        integer :: mpi_comm    ! MPI communicator for PEs having this domain
        integer :: npes_comm   ! number of PEs having this domain
        integer :: mpi_tag     ! for MPI book-keeping

c#ifdef USE_ESMF
c        TYPE (ESMF_Grid) :: ESMF_GRID
c#endif
        INTEGER :: NPES_WORLD
        INTEGER :: RANK
        INTEGER :: NPES_USED
         ! Parameters for Global domain
        INTEGER :: IM_WORLD     ! Number of Longitudes
        INTEGER :: JM_WORLD     ! Number of latitudes
         ! Parameters for local domain
        INTEGER :: I_STRT       ! Begin local domain longitude index
        INTEGER :: I_STOP       ! End   local domain longitude index
        INTEGER :: J_STRT       ! Begin local domain latitude  index
        INTEGER :: J_STOP       ! End   local domain latitude  index
        INTEGER :: J_STRT_SKP   ! Begin local domain exclusive of S pole
        INTEGER :: J_STOP_SKP   ! End   local domain exclusive of N pole
        INTEGER :: ni_loc       ! for transpose
         ! Parameters for halo of local domain
        INTEGER :: I_STRT_HALO  ! Begin halo longitude index
        INTEGER :: I_STOP_HALO  ! End   halo longitude index
        INTEGER :: J_STRT_HALO  ! Begin halo latitude  index
        INTEGER :: J_STOP_HALO  ! End   halo latitude  index
         ! Parameters for staggered "B" grid
         ! Note that global staggered grid begins at "2".
        INTEGER :: J_STRT_STGR  ! Begin local staggered domain
        INTEGER :: J_STOP_STGR  ! End   local staggered domain
         ! Controls for special cases

        INTEGER, DIMENSION(:), POINTER :: DJ_MAP
        INTEGER :: DJ
        INTEGER :: log_unit     ! for debugging
         !@var lookup_pet index of PET for a given J
        INTEGER, DIMENSION(:), POINTER :: lookup_pet
        LOGICAL :: BC_PERIODIC

      end type dist_grid

#ifndef SERIAL_MODE

c public interfaces
      public :: pack_row,unpack_row
      public :: halo_update
      public :: globalsum
      public :: isInLocalSubdomain

#endif /* not SERIAL_MODE */
      public :: pack_data,unpack_data

c
c pack/unpack interfaces
c
      interface pack_data
        module procedure do_2D_pack_1tile
        module procedure do_3D_pack_1tile
        module procedure do_4D_pack_1tile
        module procedure do_5D_pack_1tile
        module procedure do_2D_pack_multitile
        module procedure do_3D_pack_multitile
        module procedure do_4D_pack_multitile
        module procedure do_5D_pack_multitile
      end interface pack_data
      interface unpack_data
        module procedure do_2D_unpack_1tile
        module procedure do_3D_unpack_1tile
        module procedure do_4D_unpack_1tile
        module procedure do_5D_unpack_1tile
        module procedure do_2D_unpack_multitile
        module procedure do_3D_unpack_multitile
        module procedure do_4D_unpack_multitile
        module procedure do_5D_unpack_multitile
      end interface unpack_data

#ifndef SERIAL_MODE

      interface pack_row
        module procedure do_2D_pack_row
        module procedure do_3D_pack_row
        module procedure do_4D_pack_row
        module procedure do_5D_pack_row
      end interface pack_row
      interface unpack_row
        module procedure do_2D_unpack_row
        module procedure do_3D_unpack_row
        module procedure do_4D_unpack_row
        module procedure do_5D_unpack_row
      end interface unpack_row

c
c halo update interface
c
      interface halo_update
        module procedure halo_update_2D
        module procedure halo_update_3D
        module procedure halo_update_4D
c        module procedure halo_update_2D_int
      end interface halo_update

c
c globalsum interfaces
c
      interface globalsum
        module procedure globalsum_2D_r8
      end interface globalsum

#endif /* not SERIAL_MODE */

c
c get non-distributed dimensions
c
      public :: get_nlnk
      interface get_nlnk
        module procedure get_nlnk_2D
        module procedure get_nlnk_3D
        module procedure get_nlnk_4D
        module procedure get_nlnk_5D
      end interface


c
c buffers
c
      integer, parameter :: nkmax=100
      real*8, dimension(:), allocatable ::
     &     buf1d_local,buf1d_tile,bufij_tile
     &    ,bufsend,bufrecv
      logical :: reallocateBufsend
      integer ::
     &     buf1d_local_size=-1
     &    ,buf1d_tile_size=-1
     &    ,bufij_tile_size=-1
      integer :: communicator_


      contains

      subroutine init_dist_grid(
     &     npx,npy,ntiles,
     &     is,ie,js,je,
     &     isd,ied,jsd,jed,
     &     communicator, grid)
      integer, intent(in) :: npx,npy,ntiles,is,ie,js,je,isd,ied,jsd,jed
      integer, intent(in) :: communicator
      type(dist_grid) :: grid
      integer :: i,itile,ierr,iproc,n,ihem,ihem_sv,modrank,midp1
     &     ,group_world,group_mytile,group_intertile
     &     ,group_halo,halo_comm,group_row,row_comm,xpos,ypos
      integer, parameter :: ntiles_max=6
      integer :: nproc_comm
      integer, dimension(ntiles_max) ::
     &      tile_comms,tile_root_procs
      integer, dimension(:), allocatable :: proclist
      logical :: swap_ne,swap_sw

      grid%am_i_globalroot = .true. ! initialize to serial-mode default
      grid%npx = npx
      grid%npy = npy
      grid%ntiles = ntiles
      grid%is = is
      grid%ie = ie
      grid%js = js
      grid%je = je
      grid%isd = isd
      grid%ied = ied
      grid%jsd = jsd
      grid%jed = jed
      communicator_ = communicator

c      write(*,*) "init_dist_grid, is, ie, js, je"
c      write(*,*) "/ isd, ied, jsd, jed=",is,ie,js,je,isd,ied,jsd,jed
#ifdef SERIAL_MODE
      grid%gid =0
      grid%nproc = 1
#else
      call mpi_comm_rank(COMMUNICATOR,grid%gid,ierr)
      call mpi_comm_size(COMMUNICATOR,grid%nproc,ierr)
#endif
      grid%nproc_tile = grid%nproc/grid%ntiles
      grid%tile = 1+grid%gid/grid%nproc_tile

      allocate(grid%isr(grid%nproc_tile)
     &        ,grid%ier(grid%nproc_tile)
     &        ,grid%jsr(grid%nproc_tile)
     &        ,grid%jer(grid%nproc_tile)
     &        ,grid%cntsij(grid%nproc_tile)
     &        ,proclist(grid%nproc)
     &     )


c
c create communicators and define other necessary info
c
#ifndef SERIAL_MODE
      call mpi_comm_group(COMMUNICATOR,group_world,ierr)
      grid%comm_tile = COMMUNICATOR
      grid%am_i_rowroot = .true.
      grid%comm_row = MPI_COMM_NULL
      grid%comm_ew = MPI_COMM_NULL
      grid%comm_ns = COMMUNICATOR
      grid%pe_send_ne=MPI_PROC_NULL
      grid%pe_send_sw=MPI_PROC_NULL
      grid%pe_diag_s=MPI_PROC_NULL
      grid%pe_diag_r=MPI_PROC_NULL
#endif

      grid%am_i_globalroot = grid%gid.eq.0
      grid%rank_tile = grid%gid
      grid%am_i_tileroot = grid%am_i_globalroot
      grid%root_mytile = 0
      grid%rank_row = 0

      grid%rank_ew = 0
      grid%rank_ns = grid%gid
      grid%nprocx = 1
      grid%nprocy = grid%nproc
      grid%nproc_comm_ew = 1
      grid%nproc_comm_ns = grid%nprocy

      grid%i1wp=grid%is
      grid%i2wp=grid%is
      grid%j1wp=grid%js
      grid%j2wp=grid%je
      grid%i1ep=grid%ie
      grid%i2ep=grid%ie
      grid%j1ep=grid%js
      grid%j2ep=grid%je

      grid%i1sp=grid%is
      grid%i2sp=grid%ie
      grid%j1sp=grid%js
      grid%j2sp=grid%js
      grid%i1np=grid%is
      grid%i2np=grid%ie
      grid%j1np=grid%je
      grid%j2np=grid%je
      grid%iincsp=1
      grid%iincnp=1
      grid%jincsp=1
      grid%jincnp=1

      grid%i1wr=grid%is-1
      grid%i2wr=grid%is-1
      grid%j1wr=grid%js
      grid%j2wr=grid%je
      grid%i1er=grid%ie+1
      grid%i2er=grid%ie+1
      grid%j1er=grid%js
      grid%j2er=grid%je

      grid%i1sr=grid%is
      grid%i2sr=grid%ie
      grid%j1sr=grid%js-1
      grid%j2sr=grid%js-1
      grid%i1nr=grid%is
      grid%i2nr=grid%ie
      grid%j1nr=grid%je+1
      grid%j2nr=grid%je+1

#ifndef SERIAL_MODE
      if(grid%ntiles.eq.6) then ! assume cubed sphere
        grid%nprocx = sqrt(dble(grid%nproc_tile))
        grid%nprocy = grid%nprocx

c create intra-tile communicators
        do itile=1,grid%ntiles
          do i=1,grid%nproc_tile
            proclist(i) = grid%nproc_tile*(itile-1)+i-1
          enddo
          call mpi_group_incl(group_world,grid%nproc_tile,proclist,
     &         group_mytile,ierr)
          call mpi_comm_create(COMMUNICATOR,group_mytile,
     &         tile_comms(itile),ierr)
          if(itile.eq.grid%tile) then
            grid%root_mytile = proclist(1)
            grid%comm_tile = tile_comms(itile)
          endif
        enddo
        call mpi_comm_rank(grid%comm_tile,grid%rank_tile,ierr)
        grid%am_i_tileroot = grid%rank_tile.eq.0

c
c create row communicators
c
        do iproc=0,grid%nproc-1,grid%nprocx
          do i=1,grid%nprocx
            proclist(i) = iproc+i-1
          enddo
          call mpi_group_incl(group_world,grid%nprocx,proclist,
     &         group_row,ierr)
          call mpi_comm_create(COMMUNICATOR,group_row,
     &         row_comm,ierr)
          if(any(proclist(1:grid%nprocx).eq.grid%gid)) then
            grid%comm_row = row_comm
          endif
        enddo
        call mpi_comm_rank(grid%comm_row,grid%rank_row,ierr)
        grid%am_i_rowroot = grid%rank_row.eq.0

c
c halo update info
c
c
c first periodic direction, which we call "east-west"
c
        do n=1,grid%nprocx
          proclist=-1
          call genlist1(grid%nprocx,n,proclist)
          nproc_comm = count(proclist.ge.0)
          call mpi_group_incl(group_world,nproc_comm,proclist,
     &         group_halo,ierr)
          call mpi_comm_create(COMMUNICATOR,group_halo,
     &         halo_comm,ierr)
          if(any(proclist(1:nproc_comm).eq.grid%gid)) then
            grid%comm_ew = halo_comm
            grid%nproc_comm_ew = nproc_comm
          endif
        enddo
        call mpi_comm_rank(grid%comm_ew,grid%rank_ew,ierr)
        modrank=mod(grid%rank_ew,2*grid%nprocx)
        if(modrank.ge.grid%nprocx) then
c "east" is at j=je
          grid%i1ep=grid%is
          grid%i2ep=grid%ie
          grid%j1ep=grid%je
          grid%j2ep=grid%je 
          grid%i1er=grid%is
          grid%i2er=grid%ie
          grid%j1er=grid%je+1
          grid%j2er=grid%je+1
        endif
        if(modrank.eq.0 .or. modrank.gt.grid%nprocx) then
c "west" is at j=js
          grid%i1wp=grid%is
          grid%i2wp=grid%ie
          grid%j1wp=grid%js
          grid%j2wp=grid%js 
          grid%i1wr=grid%is
          grid%i2wr=grid%ie
          grid%j1wr=grid%js-1
          grid%j2wr=grid%js-1
        endif
c
c second periodic direction, which we call "north-south"
c
        do n=1,grid%nprocx
        do ihem=1,2
          proclist=-1
          if(ihem.eq.1) then
            call genlist2_odd(grid%nprocx,n,proclist)
          else
            call genlist2_even(grid%nprocx,n,proclist)
          endif
          nproc_comm = count(proclist.ge.0)
          call mpi_group_incl(group_world,nproc_comm,proclist,
     &         group_halo,ierr)
          call mpi_comm_create(COMMUNICATOR,group_halo,
     &         halo_comm,ierr)
          if(any(proclist(1:nproc_comm).eq.grid%gid)) then
            grid%comm_ns = halo_comm
            grid%nproc_comm_ns = nproc_comm
            ihem_sv=ihem
          endif
        enddo
        enddo
        call mpi_comm_rank(grid%comm_ns,grid%rank_ns,ierr)
        modrank=mod(grid%rank_ns,grid%nproc_comm_ns/3)
        midp1=(1+grid%nproc_comm_ns/3)/2
        swap_ne=.false.
        swap_sw=.false.
        if(ihem_sv.eq.1) then
          if(modrank.ge.midp1) swap_ne=.true.
          if(modrank.eq.0 .or. modrank.ge.midp1) swap_sw=.true.
        else
          if(modrank.lt.midp1) swap_ne=.true.
          if(modrank.gt.0 .and. modrank.lt.midp1) swap_sw=.true.
        endif
        if(swap_ne) then
c "north" is at i=ie
          grid%i1np=grid%ie
          grid%i2np=grid%ie
          grid%j1np=grid%js
          grid%j2np=grid%je 
          grid%i1nr=grid%ie+1
          grid%i2nr=grid%ie+1
          grid%j1nr=grid%js
          grid%j2nr=grid%je 
        endif
        if(swap_sw) then
c "south" is at i=is
          grid%i1sp=grid%is
          grid%i2sp=grid%is
          grid%j1sp=grid%js
          grid%j2sp=grid%je 
          grid%i1sr=grid%is-1
          grid%i2sr=grid%is-1
          grid%j1sr=grid%js
          grid%j2sr=grid%je
        endif
        if(modrank.eq.midp1-1) then
c change of orientation across "north" edge
          if(ihem_sv.eq.1) then
            grid%iincnp=-1
            call swap_int(grid%i1np,grid%i2np)
          else
            grid%jincnp=-1
            call swap_int(grid%j1np,grid%j2np)
          endif
        endif
        if(modrank.eq.midp1 .or. midp1.eq.1) then
c change of orientation across "south" edge
          if(ihem_sv.eq.1) then
            grid%jincsp=-1
            call swap_int(grid%j1sp,grid%j2sp)
          else
            grid%iincsp=-1
            call swap_int(grid%i1sp,grid%i2sp)
          endif
        endif

c
c extra communications for corner halos near the UL-LR diagonal
c
        if(grid%nprocx.gt.1) then
          ypos = 1+grid%rank_tile/grid%nprocx
          xpos = 1+grid%rank_tile-grid%nprocx*(ypos-1)
          if(on_ullr(xpos,ypos,grid%nprocx)) then
c extra sends needed in the nw-se direction along the diagonal
            if(grid%rank_tile.lt.grid%nproc_tile-grid%nprocx) then
              grid%pe_diag_s = grid%gid+grid%nprocx-1
            endif
            if(grid%rank_tile.gt.grid%nprocx-1) then
              grid%pe_diag_r = grid%gid-grid%nprocx+1
            endif
            if(mod(grid%tile,2).eq.0)
     &           call swap_int(grid%pe_diag_s,grid%pe_diag_r)
          else
c extra sendrecvs needed along sw-ne direction across the diagonal
            if(on_ullr(xpos+1,ypos,grid%nprocx)) then
              grid%pe_send_ne=grid%gid+grid%nprocx+1
            elseif(on_ullr(xpos-1,ypos,grid%nprocx)) then
              grid%pe_send_sw=grid%gid-grid%nprocx-1
            endif
          endif
        endif

      endif

      deallocate(proclist)

c create communicator among roots of each tile
      if(grid%ntiles.gt.1) then
        do itile=1,grid%ntiles
          tile_root_procs(itile) = (itile-1)*grid%nproc_tile
        enddo
        call mpi_group_incl(group_world,grid%ntiles,tile_root_procs,
     &       group_intertile,ierr)
        call mpi_comm_create(COMMUNICATOR,group_intertile,
     &       grid%comm_intertile,ierr)
      endif

c tabulate bounds info for each processor
      call mpi_allgather(grid%is,1,MPI_INTEGER,grid%isr,1,MPI_INTEGER,
     &     grid%comm_tile,ierr)
      call mpi_allgather(grid%ie,1,MPI_INTEGER,grid%ier,1,MPI_INTEGER,
     &     grid%comm_tile,ierr)
      call mpi_allgather(grid%js,1,MPI_INTEGER,grid%jsr,1,MPI_INTEGER,
     &     grid%comm_tile,ierr)
      call mpi_allgather(grid%je,1,MPI_INTEGER,grid%jer,1,MPI_INTEGER,
     &     grid%comm_tile,ierr)
      grid%maxnj = 1+maxval(grid%jer-grid%jsr)

#endif /* not SERIAL_MODE */

c tabulate counts and displacements for gather/scatter
      do iproc=1,grid%nproc_tile
        grid%cntsij(iproc) =
     &       (1+grid%ier(iproc)-grid%isr(iproc))*
     &       (1+grid%jer(iproc)-grid%jsr(iproc))
      enddo

c
c allocate send/receive buffers for halo updates
c assume a maximum halo width of 3, maximum nl*nk=1000
c
      n = 3*(7+grid%ie-grid%is)*1000
      reallocateBufsend = .false.
      if (.not. allocated(bufsend)) then
         reallocateBufsend = .true.
      else if (n.gt.size(bufsend)) then
         reallocateBufsend = .true.
         deallocate(bufsend)
         deallocate(bufrecv)
      end if
      if (reallocateBufsend) then
        allocate(bufsend(n),bufrecv(n))
      endif

      return
      end subroutine init_dist_grid

#undef _GATHER_
#undef _SCATTER_

c
c R8 pack routines
c
#define _GATHER_
#undef _OPER_
#define _OPER_ pack
#define MULTITILE
#include "gs_variants.inc"
#undef MULTITILE
#define TILE
#include "gs_variants.inc"
#undef TILE
#define ROW
#include "gs_variants.inc"
#undef ROW
#undef _GATHER_

c
c R8 unpack routines
c
#define _SCATTER_
#undef _OPER_
#define _OPER_ unpack
#define MULTITILE
#include "gs_variants.inc"
#undef MULTITILE
#define TILE
#include "gs_variants.inc"
#undef TILE
#define ROW
#include "gs_variants.inc"
#undef ROW
#undef _SCATTER_


#ifndef SERIAL_MODE
      subroutine halo_update_2D(grid,arr,jdim)
      real*8, dimension(:,:) :: arr
      include 'do_halo.inc'
      return
      end subroutine halo_update_2D
      subroutine halo_update_3D(grid,arr,jdim)
      real*8, dimension(:,:,:) :: arr
      include 'do_halo.inc'
      return
      end subroutine halo_update_3D
      subroutine halo_update_4D(grid,arr,jdim)
      real*8, dimension(:,:,:,:) :: arr
      include 'do_halo.inc'
      return
      end subroutine halo_update_4D

c      subroutine halo_update_2D_int(grid,iarr)
c      type(dist_grid) :: grid
c      integer, dimension(:,:) :: iarr
c      real*8 :: arr(size(iarr,1),size(iarr,2))
c      arr = iarr
c      call halo_update(grid,arr)
c      iarr = arr
c      return
c      end subroutine halo_update_2D_int

      subroutine globalsum_2D_r8(grid,arr,arrsum,all)
      type(dist_grid), intent(in) :: grid
      real*8, intent(in) :: arr(:,:)
      real*8 :: arrsum
      logical, intent(in), optional :: all
      real*8, dimension(:,:), allocatable :: arrtile
      real*8 :: arrsum_tile,arrsums(6)
      integer :: ierr
      if(grid%am_i_tileroot) then
        allocate(arrtile(grid%npx,grid%npy))
      endif
      call pack_data(grid,arr,arrtile)
      if(grid%am_i_tileroot) then
        arrsum_tile = sum(arrtile)
        deallocate(arrtile)
        if(grid%ntiles.gt.1) then
          call mpi_gather(arrsum_tile,1,MPI_DOUBLE_PRECISION,arrsums,
     &         1,MPI_DOUBLE_PRECISION,0,grid%comm_intertile,ierr)
          if(grid%am_i_globalroot) arrsum=sum(arrsums(1:grid%ntiles))
        else
          arrsum = arrsum_tile
        endif
      endif
      if(present(all)) then
        if(all) then
          call mpi_bcast(arrsum,1,MPI_DOUBLE_PRECISION,0,
     &         COMMUNICATOR_,ierr)
        endif
      endif
      return
      end subroutine globalsum_2D_r8

c      subroutine globalsum_2D_r8(grid,arr,arrsum,all)
c      type(dist_grid), intent(in) :: grid
c      real*8, intent(in) :: arr(:,:)
c      real*8 :: arrsum
c      logical, intent(in), optional :: all
c      real*16 :: arrsum_local, arrsum_global
c      real*16, dimension(:), allocatable :: arrsums
c      integer :: i,j,hi,hj,ierr
c      hi = (size(arr,1)-(1+grid%ie-grid%is))/2
c      hj = (size(arr,2)-(1+grid%je-grid%js))/2
c      arrsum_local = 0.
c      do j=1+hj,size(arr,2)-hj
c      do i=1+hi,size(arr,1)-hi
c        arrsum_local = arrsum_local + real(arr(i,j),kind=16)
c      enddo
c      enddo
c      if(grid%am_i_globalroot) allocate(arrsums(grid%nproc))
c      call mpi_gather(arrsum_local,2,MPI_DOUBLE_PRECISION,
c     &     arrsums, 2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
c      if(grid%am_i_globalroot) then
c        arrsum_global = sum(arrsums)
c        deallocate(arrsums)
c        arrsum =  real(arrsum_global,kind=8)
c      endif
c      if(present(all)) then
c        if(all) then
c          call mpi_bcast(arrsum,1,MPI_DOUBLE_PRECISION,0,
c     &         MPI_COMM_WORLD,ierr)
c        endif
c      endif
c      return
c      end subroutine globalsum_2D_r8

      subroutine genlist2_even(nproc1d,istart,pelist)
      implicit none
      integer :: nproc1d,istart
      integer, dimension(6*nproc1d) :: pelist
      integer :: n,nn,pe0
      call traverse_square_cw(nproc1d,istart,pelist)
      n=2*(nproc1d-istart)+1
      pelist(1:n) = pelist(1:n) + nproc1d**2
c copy lists for tile 2 to tiles 4,6
      pelist(n+1:2*n) = pelist(1:n) + 2*nproc1d**2
      pelist(2*n+1:3*n) = pelist(1:n) + 4*nproc1d**2
c circular shift the list so that the first element is
c on the diagonal
      n=3*n
      do nn=1,nproc1d-istart
        pe0=pelist(1)
        pelist(1:n-1)=pelist(2:n)
        pelist(n)=pe0
      enddo
      return
      end subroutine genlist2_even

      subroutine genlist2_odd(nproc1d,jstart,pelist)
      implicit none
      integer :: nproc1d,jstart
      integer, dimension(6*nproc1d) :: pelist
      integer :: n,nn,pe0
      call traverse_square_ccw(nproc1d,jstart,pelist)
      n=2*(nproc1d-jstart)+1
c copy lists for tile 1 to tiles 3,5
      pelist(n+1:2*n) = pelist(1:n) + 2*nproc1d**2
      pelist(2*n+1:3*n) = pelist(1:n) + 4*nproc1d**2
c circular shift the list so that the first element is
c on the diagonal
      n=3*n
      do nn=1,nproc1d-jstart
        pe0=pelist(1)
        pelist(1:n-1)=pelist(2:n)
        pelist(n)=pe0
      enddo
      return
      end subroutine genlist2_odd

      subroutine genlist1(nproc1d,istart,pelist)
      implicit none
      integer :: nproc1d,istart
      integer, dimension(6*nproc1d) :: pelist
      integer :: jstart,n,n2,nn,pe0
      call traverse_square_cw(nproc1d,istart,pelist(1))
      jstart=1+nproc1d-istart
      n=2*(nproc1d-istart+1)
      call traverse_square_ccw(nproc1d,jstart,pelist(n))
      n2=2*nproc1d
      pelist(n:n2) = pelist(n:n2) + nproc1d**2
c copy lists for tiles 1-2 to tiles 3-4,5-6
      n=n2
      pelist(n+1:2*n) = pelist(1:n) + 2*nproc1d**2
      pelist(2*n+1:3*n) = pelist(1:n) + 4*nproc1d**2
c circular shift the list so that the first element is
c on the diagonal
      n=6*nproc1d
      do nn=1,nproc1d-istart
        pe0=pelist(1)
        pelist(1:n-1)=pelist(2:n)
        pelist(n)=pe0
      enddo
      return
      end subroutine genlist1

      subroutine traverse_square_cw(np,istart,pelist)
      implicit none
      integer :: np,istart
      integer, dimension(2*(np-istart)+1) :: pelist
      integer :: n
      pelist(1)=istart-1
      do n=2,np-istart+1
        pelist(n)=pelist(n-1)+np
      enddo
      do n=np-istart+2,2*(np-istart)+1
        pelist(n)=pelist(n-1)+1
      enddo
      return
      end subroutine traverse_square_cw
      subroutine traverse_square_ccw(np,jstart,pelist)
      implicit none
      integer :: np,jstart
      integer, dimension(2*(np-jstart)+1) :: pelist
      integer :: n
      pelist(1)=(jstart-1)*np
      do n=2,np-jstart+1
        pelist(n)=pelist(n-1)+1
      enddo
      do n=np-jstart+2,2*(np-jstart)+1
        pelist(n)=pelist(n-1)+np
      enddo
      return
      end subroutine traverse_square_ccw

#endif /* not SERIAL_MODE */

      subroutine alloc_gs_wksp_default()
! Although these 3 arrays are not really used in a serial build, they
! are still passed to external libraries and/or interfaces in which they
! the corresponding dummy argument is assumed size.
      if (.not. allocated(buf1d_local)) then
        allocate(buf1d_local(1))
        buf1d_local_size=1
      endif
      if (.not. allocated(buf1d_tile)) then
        allocate(buf1d_tile(1))
        buf1d_tile_size=1
      endif
      if (.not. allocated(bufij_tile)) then
        allocate(bufij_tile(1))
        bufij_tile_size=1
      endif
      end subroutine alloc_gs_wksp_default

      subroutine alloc_gs_wksp(grid,nl,nk,nj,nt,am_i_gsroot)
c allocates gather/scatter workspace
      type(dist_grid) :: grid
      integer :: nl,nk,nj,nt
      integer :: lsize,tsize,nlk
      logical :: am_i_gsroot

      if(nl.eq.1) then
        nlk = min(nk,nkmax)     ! gs3D does nkmax at a time
      else
        nlk = nl                ! gs4D does one k at a time
      endif
      lsize = nlk*(1+grid%ie-grid%is)*(1+grid%je-grid%js)
      tsize = nlk*grid%npx*nj
      if(lsize.gt.buf1d_local_size) then
        if(allocated(buf1d_local)) deallocate(buf1d_local)
        allocate(buf1d_local(lsize))
        buf1d_local_size=lsize
      endif
      if(am_i_gsroot) then
        if(tsize.gt.buf1d_tile_size) then
          if(allocated(buf1d_tile)) deallocate(buf1d_tile)
          allocate(buf1d_tile(tsize))
          buf1d_tile_size=tsize
        endif
        if(nt.gt.1 .and. tsize.gt.bufij_tile_size) then
          if(allocated(bufij_tile)) deallocate(bufij_tile)
          allocate(bufij_tile(tsize))
          bufij_tile_size=tsize
        endif
      endif
      return
      end subroutine alloc_gs_wksp

      subroutine get_nlnk_2D(arr,jdim,nl,nk)
      real*8 :: arr(:,:)
      integer :: jdim
      integer :: nl,nk
      integer :: xdim
      nl = 1
      nk = 1
c      do xdim=1,jdim-2
c        nl = nl*size(arr,xdim)
c      enddo
      do xdim=jdim+1,2
        nk = nk*size(arr,xdim)
      enddo
      return
      end subroutine get_nlnk_2D
      subroutine get_nlnk_3D(arr,jdim,nl,nk)
      real*8 :: arr(:,:,:)
      integer :: jdim
      integer :: nl,nk
      integer :: xdim
      nl = 1
      nk = 1
      do xdim=1,jdim-2
        nl = nl*size(arr,xdim)
      enddo
      do xdim=jdim+1,3
        nk = nk*size(arr,xdim)
      enddo
      return
      end subroutine get_nlnk_3D
      subroutine get_nlnk_4D(arr,jdim,nl,nk)
      real*8 :: arr(:,:,:,:)
      integer :: jdim
      integer :: nl,nk
      integer :: xdim
      nl = 1
      nk = 1
      do xdim=1,jdim-2
        nl = nl*size(arr,xdim)
      enddo
      do xdim=jdim+1,4
        nk = nk*size(arr,xdim)
      enddo
      return
      end subroutine get_nlnk_4D
      subroutine get_nlnk_5D(arr,jdim,nl,nk)
      real*8 :: arr(:,:,:,:,:)
      integer :: jdim
      integer :: nl,nk
      integer :: xdim
      nl = 1
      nk = 1
      do xdim=1,jdim-2
        nl = nl*size(arr,xdim)
      enddo
      do xdim=jdim+1,5
        nk = nk*size(arr,xdim)
      enddo
      return
      end subroutine get_nlnk_5D

      subroutine swap_int(i1,i2)
      integer :: i1,i2
      integer :: itmp
      itmp = i1
      i1 = i2
      i2 = itmp
      return
      end subroutine swap_int

      function on_ullr(i,j,n)
      implicit none
      logical :: on_ullr
      integer :: i,j,n
      on_ullr = i+j.eq.n+1
      return
      end function on_ullr


! ----------------------------------------------------------------------
      logical function isInLocalSubdomain(distGrid, i, j)
! ----------------------------------------------------------------------
        type (DIST_GRID), intent(in) :: distGrid
        integer, intent(in) :: i
        integer, intent(in) :: j

        isInLocalSubdomain = 
     &       (i >= distGrid%i_strt .and. i <= distGrid%i_stop) .and.
     &       (j >= distGrid%j_strt .and. j <= distGrid%j_stop)

      end function isInLocalSubdomain

      end module dd2d_utils

      subroutine gather4D(grid,local_arr,global_arr
     &     ,i1,i2,j1,j2,nl,nk,nt,has_halo
     &     ,i1g,i2g,j1g,j2g,am_i_gsroot,comm_gs,nproc_comm
     &     ,cnts,cntsg,displs,displsg
     &     ,buf1d_local_size,buf1d_tile_size,bufij_tile_size
     &     ,buf1d_local,buf1d_tile,bufij_tile
     &     )
      use dd2d_utils, only : dist_grid
      implicit none
#ifndef SERIAL_MODE
      include 'mpif.h'
#endif
      type(dist_grid) :: grid
      integer :: i1,i2,j1,j2,nl,nk,nt
      integer :: i1g,i2g,j1g,j2g
      logical :: am_i_gsroot
      integer :: comm_gs,nproc_comm
      logical :: has_halo
      integer, dimension(nproc_comm) :: cnts,displs
      integer, dimension(nt) :: cntsg,displsg
      integer :: buf1d_local_size,buf1d_tile_size,bufij_tile_size
      real*8 :: buf1d_local(buf1d_local_size)
      real*8 :: buf1d_tile(buf1d_tile_size)
      real*8 :: bufij_tile(bufij_tile_size)
      real*8 local_arr(nl,i1:i2,j1:j2,nk)
      real*8 global_arr(nl,i1g:i2g,j1g:j2g,nk,nt)
c
      integer :: i,j,l,k,m,n,n0,nsend,iproc
      integer :: ierr
      real*8 :: r8dum

      if(nproc_comm.eq.1 .and. nt.eq.1) then
        global_arr(:,i1g:i2g,j1g:j2g,:,1)=local_arr(:,i1g:i2g,j1g:j2g,:)
        return
      endif

#ifndef SERIAL_MODE
      do k=1,nk
c
c first, collect into a buffer for this gather region with
c memory-contiguous segments from each PE.
c
        nsend = (grid%ie-grid%is+1)*(grid%je-grid%js+1)*nl
        if(has_halo) then ! get rid of halo first
          m = 0
          do j=grid%js,grid%je
          do i=grid%is,grid%ie
          do l=1,nl
            m = m + 1
            buf1d_local(m) = local_arr(l,i,j,k)
          enddo
          enddo
          enddo
          call mpi_gatherv(buf1d_local,nsend,MPI_DOUBLE_PRECISION,
     &         buf1d_tile,cnts,displs,MPI_DOUBLE_PRECISION,
     &         0,comm_gs,ierr)
        else ! use local halo-less array directly
          call mpi_gatherv(local_arr(1,i1,j1,k)
     &                                ,nsend,MPI_DOUBLE_PRECISION,
     &         buf1d_tile,cnts,displs,MPI_DOUBLE_PRECISION,
     &         0,comm_gs,ierr)
        endif

c
c Copy each segment of buf1d_tile into the corresponding i,j rectangle
c of the gather region
c
        if(am_i_gsroot) then
          if(nt.eq.1) then
            m = 0
            do iproc=1+grid%rank_tile,nproc_comm+grid%rank_tile
              do j=grid%jsr(iproc),grid%jer(iproc)
              do i=grid%isr(iproc),grid%ier(iproc)
              do l=1,nl
                m = m + 1
                global_arr(l,i,j,k,1) = buf1d_tile(m)
              enddo
              enddo
              enddo
            enddo
          else
            m = 0
            do iproc=1,grid%nproc_tile
              do j=grid%jsr(iproc),grid%jer(iproc)
              n0 = nl*grid%npx*(j-1)
              do i=grid%isr(iproc),grid%ier(iproc)
              n = n0 + nl*(i-1)
              do l=1,nl
                m = m + 1
                n = n + 1
                bufij_tile(n) = buf1d_tile(m)
              enddo
              enddo
              enddo
            enddo
c
c if a multi-tile gather, collect the tiles on global root
c
            nsend = nl*grid%npx*grid%npy
            if(grid%am_i_globalroot) then
            call mpi_gatherv(bufij_tile,nsend,MPI_DOUBLE_PRECISION,
     &         global_arr(1,1,1,k,1),cntsg,displsg,
     &           MPI_DOUBLE_PRECISION,
     &           0,grid%comm_intertile,ierr)
            else
            call mpi_gatherv(bufij_tile,nsend,MPI_DOUBLE_PRECISION,
     &         r8dum,                cntsg,displsg,
     &           MPI_DOUBLE_PRECISION,
     &           0,grid%comm_intertile,ierr)
            endif
          endif
        endif

      enddo                     ! k loop
      return
#endif /* not SERIAL_MODE */

      end subroutine gather4D

      subroutine scatter4D(grid,local_arr,global_arr
     &     ,i1,i2,j1,j2,nl,nk,nt,has_halo
     &     ,i1g,i2g,j1g,j2g,am_i_gsroot,comm_gs,nproc_comm
     &     ,cnts,cntsg,displs,displsg
     &     ,buf1d_local_size,buf1d_tile_size,bufij_tile_size
     &     ,buf1d_local,buf1d_tile,bufij_tile
     &     )
      use dd2d_utils, only : dist_grid
      implicit none
#ifndef SERIAL_MODE
      include 'mpif.h'
#endif
      type(dist_grid) :: grid
      integer :: i1,i2,j1,j2,nl,nk,nt
      integer :: i1g,i2g,j1g,j2g
      logical :: am_i_gsroot
      integer :: comm_gs,nproc_comm
      logical :: has_halo
      integer, dimension(nproc_comm) :: cnts,displs
      integer, dimension(nt) :: cntsg,displsg
      integer :: buf1d_local_size,buf1d_tile_size,bufij_tile_size
      real*8 :: buf1d_local(buf1d_local_size)
      real*8 :: buf1d_tile(buf1d_tile_size)
      real*8 :: bufij_tile(bufij_tile_size)
      real*8 local_arr(nl,i1:i2,j1:j2,nk)
      real*8 global_arr(nl,i1g:i2g,j1g:j2g,nk,nt)
c
      integer :: i,j,l,k,m,n,n0,nrecv,iproc
      integer :: ierr
      real*8 :: r8dum

      if(nproc_comm.eq.1 .and. nt.eq.1) then
        local_arr(:,i1g:i2g,j1g:j2g,:)=global_arr(:,i1g:i2g,j1g:j2g,:,1)
        return
      endif

#ifndef SERIAL_MODE
      do k=1,nk


c
c first, rearrange the data for this scatter region into
c memory-contiguous segments destined for each PE.
c

        if(am_i_gsroot) then

          if(nt.eq.1) then
            m = 0
            do iproc=1+grid%rank_tile,nproc_comm+grid%rank_tile
              do j=grid%jsr(iproc),grid%jer(iproc)
              do i=grid%isr(iproc),grid%ier(iproc)
              do l=1,nl
                m = m + 1
                buf1d_tile(m) = global_arr(l,i,j,k,1)
              enddo
              enddo
              enddo
            enddo

          else
c
c if a multi-tile scatter, global root first sends to tile roots.
c
            nrecv = nl*grid%npx*grid%npy
            if(grid%am_i_globalroot) then
            call mpi_scatterv(global_arr(1,1,1,k,1),cntsg,displsg,
     &           MPI_DOUBLE_PRECISION,
     &           bufij_tile,nrecv,MPI_DOUBLE_PRECISION,
     &           0,grid%comm_intertile,ierr)
            else
            call mpi_scatterv(r8dum,                cntsg,displsg,
     &           MPI_DOUBLE_PRECISION,
     &           bufij_tile,nrecv,MPI_DOUBLE_PRECISION,
     &           0,grid%comm_intertile,ierr)
            endif
            m = 0
            do iproc=1,grid%nproc_tile
              do j=grid%jsr(iproc),grid%jer(iproc)
              n0 = nl*grid%npx*(j-1)
              do i=grid%isr(iproc),grid%ier(iproc)
              n = n0 + nl*(i-1)
              do l=1,nl
                m = m + 1
                n = n + 1
                buf1d_tile(m) = bufij_tile(n)
              enddo
              enddo
              enddo
            enddo
          endif
        endif

        nrecv = nl*(grid%ie-grid%is+1)*(grid%je-grid%js+1)
        if(has_halo) then
          call mpi_scatterv(buf1d_tile,cnts,displs,
     &         MPI_DOUBLE_PRECISION,buf1d_local,nrecv,
     &         MPI_DOUBLE_PRECISION,
     &         0,comm_gs,ierr)
c copy the the local receive buffer into the local array
          m = 0
          do j=grid%js,grid%je
          do i=grid%is,grid%ie
          do l=1,nl
            m = m + 1
            local_arr(l,i,j,k) = buf1d_local(m)
          enddo
          enddo
          enddo
        else
          call mpi_scatterv(buf1d_tile,cnts,displs,
     &         MPI_DOUBLE_PRECISION,local_arr(1,i1,j1,k),nrecv,
     &         MPI_DOUBLE_PRECISION,
     &         0,comm_gs,ierr)
        endif

      enddo                     ! k loop

      return
#endif /* not SERIAL_MODE */

      end subroutine scatter4D

      subroutine gather3D(grid,local_arr,global_arr
     &     ,i1,i2,j1,j2,nl,nk,nt,has_halo,nkmax
     &     ,i1g,i2g,j1g,j2g,am_i_gsroot,comm_gs,nproc_comm
     &     ,cntsij,displsij,cntsijg,displsijg
     &     ,buf1d_local_size,buf1d_tile_size,bufij_tile_size
     &     ,buf1d_local,buf1d_tile,bufij_tile
     &     )
      use dd2d_utils, only : dist_grid
      implicit none
#ifndef SERIAL_MODE
      include 'mpif.h'
#endif
      type(dist_grid) :: grid
      integer :: i1,i2,j1,j2,nl,nk,nt
      integer :: i1g,i2g,j1g,j2g
      logical :: am_i_gsroot
      integer :: comm_gs,nproc_comm
      logical :: has_halo
      integer :: nkmax
      integer, dimension(nproc_comm) :: cntsij,displsij
      integer, dimension(nt) :: cntsijg,displsijg
      integer :: buf1d_local_size,buf1d_tile_size,bufij_tile_size
      real*8 :: buf1d_local(buf1d_local_size)
      real*8 :: buf1d_tile(buf1d_tile_size)
      real*8 :: bufij_tile(bufij_tile_size)
      real*8 local_arr(i1:i2,j1:j2,nk)
      real*8 global_arr(i1g:i2g,j1g:j2g,nk,nt)
c
      integer :: i,j,k,k1,k2,nk12,m,n,n0,nsend,iproc
      integer :: ierr
      integer, dimension(6) :: cntsijkg,displsijkg
      integer, dimension(:), allocatable :: cntsijk,displsijk
      real*8 :: r8dum

      if(nproc_comm.eq.1 .and. nt.eq.1) then
        global_arr(i1g:i2g,j1g:j2g,:,1)=local_arr(i1g:i2g,j1g:j2g,:)
        return
      endif

#ifndef SERIAL_MODE
      allocate(cntsijk(nproc_comm),displsijk(nproc_comm))

      k1 = 1
      k2 = min(nk,nkmax)

      do while(k1.le.nk)

        nk12 = 1 + k2 - k1

c calculate gatherv/scatterv info
        if(am_i_gsroot) then
          cntsijk(:) = cntsij(:)*nk12
          displsijk(:) = displsij(:)*nk12
        endif

c
c first, collect into a buffer for this gather region with
c memory-contiguous segments from each PE.
c
        m = 0
        do k=k1,k2
        do j=grid%js,grid%je
        do i=grid%is,grid%ie
          m = m + 1
          buf1d_local(m) = local_arr(i,j,k)
        enddo
        enddo
        enddo

        nsend = (grid%ie-grid%is+1)*(grid%je-grid%js+1)*nk12
        call mpi_gatherv(buf1d_local,nsend,MPI_DOUBLE_PRECISION,
     &       buf1d_tile,cntsijk,displsijk,MPI_DOUBLE_PRECISION,
     &       0,comm_gs,ierr)

c
c Copy each segment of buf1d_tile into the corresponding i,j rectangle
c of the gather region
c
        if(am_i_gsroot) then
          if(nt.eq.1) then
            m = 0
            do iproc=1+grid%rank_tile,nproc_comm+grid%rank_tile
              do k=k1,k2
              do j=grid%jsr(iproc),grid%jer(iproc)
              do i=grid%isr(iproc),grid%ier(iproc)
                m = m + 1
                global_arr(i,j,k,1) = buf1d_tile(m)
              enddo
              enddo
              enddo
            enddo
          else
            m = 0
            do iproc=1,grid%nproc_tile
              do k=k1,k2
              n0 = grid%npx*grid%npy*(k-k1)
              do j=grid%jsr(iproc),grid%jer(iproc)
              n = n0 + grid%npx*(j-1) + grid%isr(iproc)-1
              do i=grid%isr(iproc),grid%ier(iproc)
                m = m + 1
                n = n + 1
                bufij_tile(n) = buf1d_tile(m)
              enddo
              enddo
              enddo
            enddo
c
c if a multi-tile gather, collect the tiles on global root
c
            cntsijkg(1:nt) = cntsijg(1:nt)*nk12
            displsijkg(1:nt) = displsijg(1:nt)*nk
            nsend = nk12*grid%npx*grid%npy
            if(grid%am_i_globalroot) then
            call mpi_gatherv(bufij_tile,nsend,MPI_DOUBLE_PRECISION,
     &           global_arr(1,1,k1,1),cntsijkg,displsijkg,
     &           MPI_DOUBLE_PRECISION,
     &           0,grid%comm_intertile,ierr)
            else
            call mpi_gatherv(bufij_tile,nsend,MPI_DOUBLE_PRECISION,
     &           r8dum,               cntsijkg,displsijkg,
     &           MPI_DOUBLE_PRECISION,
     &           0,grid%comm_intertile,ierr)
            endif
          endif
        endif

        k1 = k2 + 1
        k2 = min(nk,k2+nkmax)

      enddo                     ! while k1.le.nk

      deallocate(cntsijk,displsijk)

      return
#endif /* not SERIAL_MODE */
      end subroutine gather3D

      subroutine scatter3D(grid,local_arr,global_arr
     &     ,i1,i2,j1,j2,nl,nk,nt,has_halo,nkmax
     &     ,i1g,i2g,j1g,j2g,am_i_gsroot,comm_gs,nproc_comm
     &     ,cntsij,displsij,cntsijg,displsijg
     &     ,buf1d_local_size,buf1d_tile_size,bufij_tile_size
     &     ,buf1d_local,buf1d_tile,bufij_tile
     &     )
      use dd2d_utils, only : dist_grid
      implicit none
#ifndef SERIAL_MODE
      include 'mpif.h'
#endif
      type(dist_grid) :: grid
      integer :: i1,i2,j1,j2,nl,nk,nt
      integer :: i1g,i2g,j1g,j2g
      logical :: am_i_gsroot
      integer :: comm_gs,nproc_comm
      logical :: has_halo
      integer :: nkmax
      integer, dimension(nproc_comm) :: cntsij,displsij
      integer, dimension(nt) :: cntsijg,displsijg
      integer :: buf1d_local_size,buf1d_tile_size,bufij_tile_size
      real*8 :: buf1d_local(buf1d_local_size)
      real*8 :: buf1d_tile(buf1d_tile_size)
      real*8 :: bufij_tile(bufij_tile_size)
      real*8 local_arr(i1:i2,j1:j2,nk)
      real*8 global_arr(i1g:i2g,j1g:j2g,nk,nt)
c
      integer :: i,j,l,k,k1,k2,nk12,m,n,n0,nrecv,iproc
      integer :: ierr
      integer, dimension(6) :: cntsijkg,displsijkg
      integer, dimension(:), allocatable :: cntsijk,displsijk
      real*8 :: r8dum

      if(nproc_comm.eq.1 .and. nt.eq.1) then
        local_arr(i1g:i2g,j1g:j2g,:)=global_arr(i1g:i2g,j1g:j2g,:,1)
        return
      endif

#ifndef SERIAL_MODE
      allocate(cntsijk(nproc_comm),displsijk(nproc_comm))

      k1 = 1
      k2 = min(nk,nkmax)

      do while(k1.le.nk)

        nk12 = 1 + k2 - k1

c
c first, rearrange the data of the scatter region into
c memory-contiguous segments destined for each PE.
c

        if(am_i_gsroot) then

          if(nt.eq.1) then
            m = 0
            do iproc=1+grid%rank_tile,nproc_comm+grid%rank_tile
              do k=k1,k2
              do j=grid%jsr(iproc),grid%jer(iproc)
              do i=grid%isr(iproc),grid%ier(iproc)
                m = m + 1
                buf1d_tile(m) = global_arr(i,j,k,1)
              enddo
              enddo
              enddo
            enddo

          else
c
c if a multi-tile scatter, global root first sends to tile roots.
c
            cntsijkg(1:nt) = cntsijg(1:nt)*nk12
            displsijkg(1:nt) = displsijg(1:nt)*nk
            nrecv = nk12*grid%npx*grid%npy
            if(grid%am_i_globalroot) then
            call mpi_scatterv(global_arr(1,1,k1,1),cntsijkg,displsijkg,
     &           MPI_DOUBLE_PRECISION,
     &           bufij_tile,nrecv,MPI_DOUBLE_PRECISION,
     &           0,grid%comm_intertile,ierr)
            else
            call mpi_scatterv(r8dum,               cntsijkg,displsijkg,
     &           MPI_DOUBLE_PRECISION,
     &           bufij_tile,nrecv,MPI_DOUBLE_PRECISION,
     &           0,grid%comm_intertile,ierr)
            endif
            m = 0
            do iproc=1,grid%nproc_tile
              do k=k1,k2
              n0 = grid%npx*grid%npy*(k-k1)
              do j=grid%jsr(iproc),grid%jer(iproc)
              n = n0 + grid%npx*(j-1) + grid%isr(iproc)-1
              do i=grid%isr(iproc),grid%ier(iproc)
                m = m + 1
                n = n + 1
                buf1d_tile(m) = bufij_tile(n)
              enddo
              enddo
              enddo
            enddo
          endif
        endif

c calculate gatherv/scatterv info
        if(am_i_gsroot) then
          cntsijk(:) = cntsij(:)*nk12
          displsijk(:) = displsij(:)*nk12
        endif

        nrecv = nk12*(grid%ie-grid%is+1)*(grid%je-grid%js+1)
        call mpi_scatterv(buf1d_tile,cntsijk,displsijk,
     &       MPI_DOUBLE_PRECISION,buf1d_local,nrecv,
     &       MPI_DOUBLE_PRECISION,
     &       0,comm_gs,ierr)

c copy the the receive buffer into the local array
        m = 0
        do k=k1,k2
        do j=grid%js,grid%je
        do i=grid%is,grid%ie
          m = m + 1
          local_arr(i,j,k) = buf1d_local(m)
        enddo
        enddo
        enddo

        k1 = k2 + 1
        k2 = min(nk,k2+nkmax)

      enddo                     ! while k1.le.nk

      deallocate(cntsijk,displsijk)

      return
#endif /* not SERIAL_MODE */
      end subroutine scatter3D

#ifndef SERIAL_MODE
      subroutine sendrecv4D(arr,nl,nk,i1,i2,j1,j2,
     &     i1p,i2p,j1p,j2p,iincp,jincp,i1r,i2r,j1r,j2r,
     &     mpi_comm,pe_send,pe_recv,
     &     bufsend,bufrecv
     &     )
      implicit none
#ifndef SERIAL_MODE
      include 'mpif.h'
#endif
      integer :: i1,i2,j1,j2,nl,nk,
     &     i1p,i2p,j1p,j2p,iincp,jincp,i1r,i2r,j1r,j2r,
     &     mpi_comm,pe_send,pe_recv
      real*8 :: bufsend(1),bufrecv(1)

      integer :: i,j,l,k,m,nsend
      integer :: mpi_status(MPI_STATUS_SIZE), ierr, send_tag, recv_tag
      real*8 arr(nl,i1:i2,j1:j2,nk)

c
c pack the send buffer
c
      m = 0
      if(nl.gt.1) then
        do k=1,nk
        do j=j1p,j2p,jincp
        do i=i1p,i2p,iincp
        do l=1,nl
          m = m + 1
          bufsend(m) = arr(l,i,j,k)
        enddo ! l
        enddo ! i
        enddo ! j
        enddo ! k
      else
        do k=1,nk
        do j=j1p,j2p,jincp
        do i=i1p,i2p,iincp
          m = m + 1
          bufsend(m) = arr(1,i,j,k)
        enddo ! i
        enddo ! j
        enddo ! k
      endif
      nsend = m

c
c call sendrecv
c
      send_tag = 1024
      recv_tag = 1024
      call mpi_sendrecv(
     &     bufsend,nsend,MPI_DOUBLE_PRECISION,pe_send,send_tag,
     &     bufrecv,nsend,MPI_DOUBLE_PRECISION,pe_recv,recv_tag,
     &     mpi_comm,mpi_status,ierr)

c
c unpack the receive buffer
c
      if(pe_recv.ne.MPI_PROC_NULL) then
      m = 0
      if(nl.gt.1) then
        do k=1,nk
        do j=j1r,j2r
        do i=i1r,i2r
        do l=1,nl
          m = m + 1
          arr(l,i,j,k) = bufrecv(m)
        enddo ! l
        enddo ! i
        enddo ! j
        enddo ! k
      else
        do k=1,nk
        do j=j1r,j2r
        do i=i1r,i2r
          m = m + 1
          arr(1,i,j,k) = bufrecv(m)
        enddo ! i
        enddo ! j
        enddo ! k
      endif
      endif
      return
      end subroutine sendrecv4D

#endif /* not SERIAL_MODE */
