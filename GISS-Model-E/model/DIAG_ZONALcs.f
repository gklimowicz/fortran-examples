!@sum  DIAG_ZONAL defines the resolution and array bounds for zonal
!@+    diagnostics, provides gather/scatter procedures for
!@+    accumulation arrays, and will soon include other code
!@+    and better documentation.
!@auth Cubed Sphere Development Team
      module diag_zonal
      USE CONSTANT, only : twopi
      USE RESOLUTION, only : im,jm
      USE DOMAIN_DECOMP_ATM, only : dist_grid,sumxpe,am_i_root
      use precision_mod, only : reduce_precision
      implicit none
      private

!@param IMLON,JMLAT latlon grid sizes giving approx. equivalent res.
      INTEGER, PARAMETER, public :: JMLAT=2*JM,IMLON=4*IM,IMLONH=2*IM
!@param JM_BUDG number of zigzag bands for budget page diagnostics.
!@+     JM_BUDG is small compared to JMLAT to ensure that the longitudinal
!@+     variation of the latitude of the boundary between any two adjacent
!@+     bands is small compared to the latitudinal extent of those bands.
!@+     Since some lat-lon ocean conservation diagnostics are currently
!@+     dimensioned by JM_BUDG, the allowed values of JM_BUDG are
!@+     determined by the resolutions of the ocean model at this time.
! quantize JM_BUDG at 24, 46, or 90 using integer arithmetic:
      INTEGER, PARAMETER, public :: JM_BUDG = 24
     &     +(46-24)*(min(jm, 72)/ 72) ! 46 bands at C72  cube res.
     &     +(90-46)*(min(jm,144)/144) ! 90          C144

!@var XWON scale factor for diag. printout needed for Wonderland model
      REAL*8, public :: XWON = 1d0

      public :: pack_lc,unpack_lc,get_alloc_bounds

      interface pack_lc
        module procedure pack_lc_2d
        module procedure pack_lc_3d
        module procedure pack_lc_4d
      end interface pack_lc

      interface unpack_lc
        module procedure unpack_lc_2d
        module procedure unpack_lc_3d
        module procedure unpack_lc_4d
      end interface unpack_lc

      contains

c      subroutine get_bounds(grid,
c     &     j_strt_budg,j_stop_budg,
c     &     j_strt_jk,j_stop_jk
c     &     )
c      type(dist_grid), intent(in) :: grid
c      integer, optional, intent(out) :: j_strt_budg,j_stop_budg
c      integer, optional, intent(out) :: j_strt_jk,j_stop_jk
c      if(present(j_strt_budg)) j_strt_budg = grid%j_strt
c      if(present(j_stop_budg)) j_stop_budg = grid%j_stop
c      if(present(j_strt_jk)) j_strt_jk = grid%j_strt
c      if(present(j_stop_jk)) j_stop_jk = grid%j_stop
c      return
c      end subroutine get_bounds

      subroutine get_alloc_bounds(grid,
     &     j_strt_budg,j_stop_budg,
     &     j_strt_jk,j_stop_jk
     &     )
      type(dist_grid), intent(in) :: grid
      integer, optional, intent(out) :: j_strt_budg,j_stop_budg
      integer, optional, intent(out) :: j_strt_jk,j_stop_jk
      if(present(j_strt_budg)) j_strt_budg = 1
      if(present(j_stop_budg)) j_stop_budg = jm_budg
      if(present(j_strt_jk)) j_strt_jk = 1
      if(present(j_stop_jk)) j_stop_jk = jmlat
      return
      end subroutine get_alloc_bounds

      subroutine pack_lc_2d(grid,arr_loc,arr_glob)
      type(dist_grid), intent(in) :: grid
      real*8, intent(inout) :: arr_loc(:,:)
      real*8, intent(inout) :: arr_glob(:,:)
      call sumxpe(arr_loc, arr_glob, increment=.true.)
      if(am_i_root()) call reduce_precision(arr_glob,1d-8)
      arr_loc=0
      return
      end subroutine pack_lc_2d
      subroutine pack_lc_3d(grid,arr_loc,arr_glob)
      type(dist_grid), intent(in) :: grid
      real*8, intent(inout) :: arr_loc(:,:,:)
      real*8, intent(inout) :: arr_glob(:,:,:)
      call sumxpe(arr_loc, arr_glob, increment=.true.)
      if(am_i_root()) call reduce_precision(arr_glob,1d-8)
      arr_loc=0
      return
      end subroutine pack_lc_3d
      subroutine pack_lc_4d(grid,arr_loc,arr_glob)
      type(dist_grid), intent(in) :: grid
      real*8, intent(inout) :: arr_loc(:,:,:,:)
      real*8, intent(inout) :: arr_glob(:,:,:,:)
      call sumxpe(arr_loc, arr_glob, increment=.true.)
      if(am_i_root()) call reduce_precision(arr_glob,1d-8)
      arr_loc=0
      return
      end subroutine pack_lc_4d

      subroutine unpack_lc_2d(grid,arr_glob,arr_loc)
      type(dist_grid), intent(in) :: grid
      real*8, intent(in)   :: arr_glob(:,:)
      real*8, intent(out)  :: arr_loc(:,:)
      arr_loc=0
      return
      end subroutine unpack_lc_2d
      subroutine unpack_lc_3d(grid,arr_glob,arr_loc)
      type(dist_grid), intent(in) :: grid
      real*8, intent(in)   :: arr_glob(:,:,:)
      real*8, intent(out)  :: arr_loc(:,:,:)
      arr_loc=0
      return
      end subroutine unpack_lc_3d
      subroutine unpack_lc_4d(grid,arr_glob,arr_loc)
      type(dist_grid), intent(in) :: grid
      real*8, intent(in)   :: arr_glob(:,:,:,:)
      real*8, intent(out)  :: arr_loc(:,:,:,:)
      arr_loc=0
      return
      end subroutine unpack_lc_4d

      end module diag_zonal
