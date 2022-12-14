#include "rundeck_opts.h"

      module ent_com
!@sum  ENT_COM contains the data needed for Dynamic Vegetation Model (ENT)
!@auth I. Aleinov
      use resolution, only : im,jm
      use ghy_com, only : ngm,imt,nlsn
      use ent_mod
      implicit none
      save

!@var entcells structures which keep the internal state of each Ent cell
      type(entcelltype_public), allocatable :: entcells(:,:)

!---  boundary conditions (read from file)
!@var vdata(:,:,k)  fraction of gridbox of veg.type k=1-12
!      real*8, ALLOCATABLE, dimension(:,:,:) :: vdata

!---  prognostic variables (saved to restart file)
!@var Cint Internal foliage CO2 concentration (mol/m3)
      real*8, ALLOCATABLE, dimension(:,:) :: Cint
!@var Qfol Foliage surface mixing ratio (kg/kg)
      real*8, ALLOCATABLE, dimension(:,:) :: Qfol
!@var cnc_ij canopy conductance
      real*8, ALLOCATABLE, dimension(:,:) :: cnc_ij
!@var excess_C extra land carbon accumulated due to structural
!@+   changes in LSM (kg/m^2). It should be redistributed in the 
!@+   atmosphere
      real*8, ALLOCATABLE, dimension(:,:) :: excess_C

!---  for I/O
      CHARACTER*80, parameter :: ENT_HEADER = "ENT01"
      integer, parameter :: ENT_IO_MAXBUF = 1023 !1023 after NK !924 ! 575


      contains

!/* the pario module needs a regular data layout */
#define ENT_IO_PLAIN_ARRAY

!#define ENT_IO_PLAIN_ARRAY
!!!#ifdef ENT_IO_PLAIN_ARRAY

      subroutine copy_array_to_ent_state( array )
!@sum get ent state from a simple k-IJ array
      use domain_decomp_atm, only : grid, getDomainBounds
      real*8, dimension(ENT_IO_MAXBUF,
     &     grid%i_strt_halo:grid%i_stop_halo,
     &     grid%j_strt_halo:grid%j_stop_halo) ::  array
      !---
      integer i, j, J_0, J_1, I_0, I_1

      call getDomainBounds(grid, 
     &     J_STRT=J_0, J_STOP=J_1, I_STRT=I_0, I_STOP=I_1)

      do j=J_0,J_1
        do i=I_0,I_1
          if ( array(1,i,j) > 0.d0 ) then ! the cell is present
            call ent_cell_construct( entcells(i,j) )
            call ent_cell_unpack(array(:,i,j), entcells(i,j))
          endif
        enddo
      enddo
 
      return
      end subroutine copy_array_to_ent_state

      subroutine copy_ent_state_to_array( array )
!@sum store ent state in a simple k-IJ array
      use domain_decomp_atm, only : grid, getDomainBounds
      real*8, intent(out), dimension(ENT_IO_MAXBUF,
     &     grid%i_strt_halo:grid%i_stop_halo,
     &     grid%j_strt_halo:grid%j_stop_halo) ::  array
      !---
      real*8, pointer :: cell_buf(:)
      integer i, j, J_0, J_1, I_0, I_1

      call getDomainBounds(grid, 
     &     J_STRT=J_0, J_STOP=J_1, I_STRT=I_0, I_STOP=I_1)

      nullify( cell_buf )

      do j=J_0,J_1
        do i=I_0,I_1
          array(:,i,j) = 0.d0
          call ent_cell_pack(cell_buf, entcells(i,j))
          if( size(cell_buf) > ENT_IO_MAXBUF) then
            print *,"ENT_IO_MAXBUF too small, set to",size(cell_buf)
            call stop_model("ent_write...: ENT_IO_MAXBUF too small",255)
          endif
          array(1:size(cell_buf),i,j) = cell_buf
          deallocate( cell_buf )
        enddo
      enddo

      return
      end subroutine copy_ent_state_to_array

      subroutine ent_read_state_plain( kunit, retcode )
!@sum read ent state from the file
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : am_i_root, getDomainBounds
      use domain_decomp_1d, only : UNPACK_COLUMN, broadcast
      !type(entcelltype_public), intent(out) :: entcells(:,:)
      integer, intent(in) :: kunit
      integer, intent(out) :: retcode
      !---
      CHARACTER*80 :: HEADER
      real*8, allocatable ::  buf(:,:,:), buf_glob(:,:,:)
      integer i, j, J_0H, J_1H, I_0H, I_1H 

      retcode = 0
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &     I_STRT_HALO=I_0H, I_STOP_HALO=I_1H)

      allocate( buf     (ENT_IO_MAXBUF, I_0H:I_1H, J_0H:J_1H) )

      if (AM_I_ROOT()) then
        allocate( buf_glob(ENT_IO_MAXBUF,im,jm) )
        HEADER = "        "
        READ(kunit,err=10,end=5) HEADER(1:4)
 5      continue
        if (HEADER(1:4) .ne. ENT_HEADER(1:4) ) retcode = 1
!     &       call stop_model("ent_read_state: incompatimle header",255)
        BACKSPACE kunit
        if ( retcode == 0 ) READ (kunit,err=10) HEADER, buf_glob
      endif
      call broadcast(grid, retcode)
      if ( retcode == 0 ) then
        CALL UNPACK_COLUMN(grid, buf_glob, buf)
        call copy_array_to_ent_state( buf )
      endif
      if (AM_I_ROOT()) deallocate( buf_glob )
      deallocate( buf )

      return
 10   continue
      call stop_model("ent_read_state: error reading",255)

      end subroutine ent_read_state_plain

      subroutine ent_write_state_plain( kunit )
!@sum write ent state to the file
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : am_i_root, getDomainBounds
      use domain_decomp_1d, only : PACK_COLUMN
      !use ent_com, only : entcells
      integer, intent(in) :: kunit
      !---
      real*8, allocatable ::  buf(:,:,:), buf_glob(:,:,:)
      integer i, j, J_0H, J_1H, I_0H, I_1H 

      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &     I_STRT_HALO=I_0H, I_STOP_HALO=I_1H)

      allocate( buf     (ENT_IO_MAXBUF, I_0H:I_1H, J_0H:J_1H) )
      if(AM_I_ROOT()) then
        allocate( buf_glob(ENT_IO_MAXBUF,im,jm) )
      else
        allocate( buf_glob(1,1,1) )
      end if

      call copy_ent_state_to_array( buf )
      CALL PACK_COLUMN(grid, buf, buf_glob)
      deallocate( buf )

      if (AM_I_ROOT()) then
        WRITE (kunit,err=10) ENT_HEADER, buf_glob
      endif
      deallocate( buf_glob )

      return
 10   continue
      call stop_model("ent_read_state: error writing",255)

      end subroutine ent_write_state_plain

!!!#else

      subroutine ent_read_state( kunit )
!@sum read ent state from the file
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : am_i_root, getDomainBounds
      use domain_decomp_1d, only : send_to_j, recv_from_j
      !type(entcelltype_public), intent(out) :: entcells(:,:)
      integer, intent(in) :: kunit
      !---
      !integer, parameter :: MAX_BUFFER=100000 ! need realistic estimate
      real*8, pointer ::  buffer(:)
      integer i, j, J_0, J_1
      integer tag, itag, bufsize

      nullify(buffer)

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

      !call openunit('ent_state',iu_entstate,.true.,.true.)
      do j=1,jm
        do i=1,im
          tag = (i-1) + (j-1)*IM
          itag = tag + IM*JM

          if (am_i_root()) then
            read(kunit) bufsize
            allocate( buffer(bufsize) )
            read(kunit) buffer
            !print *, "ent_read_state: i, j ", i, j
            !print *, buffer
            if ( j<J_0 .or. j>J_1 ) then  ! j is not on root
              call send_to_j(grid,size(buffer),j,itag)
              call send_to_j(grid,buffer,j,tag)
            endif
          endif

          if ( j>=J_0 .and. j<=J_1 ) then
            if ( .not.associated(buffer) ) then
              call recv_from_j(grid,bufsize,1,itag)
              allocate(buffer(bufsize))
              call recv_from_j(grid,buffer,1,tag)
            endif
            ! check length of buffer : if( buffer(1) > MAX_BUFFER ) ??
            if ( buffer(1) > 0.d0 ) then ! the cell is present
              call ent_cell_construct( entcells(i,j) )
              call ent_cell_unpack(buffer, entcells(i,j))
            endif
            deallocate( buffer )
            !this print will not work on mpi!!
            !call ent_cell_print(999,entcells(i,j))
          endif

        enddo
      enddo

      end subroutine ent_read_state


      subroutine ent_write_state( kunit )
!@sum write ent state to the file
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : am_i_root, getDomainBounds
      use domain_decomp_1d, only : send_to_j, recv_from_j
      !use ent_com, only : entcells
      integer, intent(in) :: kunit
      !---
      real*8, pointer :: buffer(:)
      integer i, j, J_0, J_1
      integer, save :: counter=0
      integer tag, itag, bufsize

      nullify(buffer)

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

      !ic = size(entcells,1)
      !jc = size(entcells,2)

      counter = mod(counter+1,2)

      !call openunit('ent_state_new',iu_entstate,.true.,.false.)
      do j=1,jm
        do i=1,im
          tag = (i-1) + (j-1)*IM
          itag = tag + IM*JM
          if ( j>=J_0 .and. j<=J_1 ) then
            call ent_cell_pack(buffer, entcells(i,j))
            if (.not.am_i_root()) then
              call send_to_j(grid,size(buffer),1,itag)
              call send_to_j(grid,buffer,1,tag)
              deallocate(buffer)
            endif
          endif
          if (am_i_root()) then
            if (.not. associated(buffer) ) then
              call recv_from_j(grid,bufsize,j,itag)
              allocate(buffer(bufsize))
              call recv_from_j(grid,buffer,j,tag)
            endif
            write(kunit) size(buffer)
            write(kunit) buffer
            deallocate(buffer)
            !this print will not work on mpi!!
            !call ent_cell_print(990+counter,entcells(i,j))
          endif
        enddo
      enddo

      end subroutine ent_write_state

!!!#endif

      end module ent_com

      SUBROUTINE ALLOC_ENT_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE ENT_MOD
      USE ENT_COM
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, getDomainBounds
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H, J_1, J_0, I_1H, I_0H, I_1, I_0
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &     J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE(    entcells(I_0:I_1,J_0:J_1),
     *         STAT=IER)
      ! initialize ent cells to something meaningful
      !call ent_cell_construct( entcells ) ! moved to init_module_ent
      call ent_cell_nullify( entcells )


      ALLOCATE(   ! vdata(I_0H:I_1H,J_0H:J_1H,12),
     *              Cint(I_0H:I_1H,J_0H:J_1H),
     *              Qfol(I_0H:I_1H,J_0H:J_1H),
     *            cnc_ij(I_0H:I_1H,J_0H:J_1H),
     *            excess_C(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)


      END SUBROUTINE ALLOC_ENT_COM

      subroutine def_rsf_vegetation(fid)
!@sum  def_rsf_vegetation defines vegetation array structure in restart files
!@auth M. Kelley
!@ver  beta
      use ent_com, only : Cint, Qfol, cnc_ij, excess_C, ent_io_maxbuf
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      real*8, dimension(ENT_IO_MAXBUF,
     &     grid%i_strt_halo:grid%i_stop_halo,
     &     grid%j_strt_halo:grid%j_stop_halo) ::  ent_array
      call defvar(grid,fid,cint,'cint(dist_im,dist_jm)')
      call defvar(grid,fid,qfol,'qfol(dist_im,dist_jm)')
      call defvar(grid,fid,cnc_ij,'cnc_ij(dist_im,dist_jm)')
      call defvar(grid,fid,excess_C,'excess_C(dist_im,dist_jm)')
      call defvar(grid,fid,ent_array,
     &     'ent_state(ent_io_maxbuf,dist_im,dist_jm)')
      return
      end subroutine def_rsf_vegetation

      subroutine new_io_vegetation(fid,iaction)
!@sum  new_io_vegetation read/write vegetation arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use ent_com, only : Cint, Qfol, cnc_ij, excess_C, ent_io_maxbuf,
     &     copy_ent_state_to_array,copy_array_to_ent_state
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      real*8, dimension(ENT_IO_MAXBUF,
     &     grid%i_strt_halo:grid%i_stop_halo,
     &     grid%j_strt_halo:grid%j_stop_halo) ::  ent_array
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid, fid, 'cint', cint)
        call write_dist_data(grid, fid, 'qfol', qfol)
        call write_dist_data(grid, fid, 'cnc_ij', cnc_ij)
        call write_dist_data(grid, fid, 'excess_C', excess_C)
        call copy_ent_state_to_array(ent_array)
        call write_dist_data(grid, fid, 'ent_state', ent_array, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'cint', cint)
        call read_dist_data(grid, fid, 'qfol', qfol)
        call read_dist_data(grid, fid, 'cnc_ij', cnc_ij)
        call read_dist_data(grid, fid, 'excess_C', excess_C)
        call read_dist_data(grid, fid, 'ent_state', ent_array, jdim=3)
        call copy_array_to_ent_state(ent_array)
      end select
      return
      end subroutine new_io_vegetation
