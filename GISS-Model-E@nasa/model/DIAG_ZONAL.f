!@sum  DIAG_ZONAL defines the resolution and array bounds for zonal
!@+    diagnostics, provides gather/scatter procedures for
!@+    accumulation arrays, and will soon include other code
!@+    and better documentation.
      module diag_zonal
      USE CONSTANT, only : twopi
      USE RESOLUTION, only : im,jm
      USE GEOM, only : dlon
      use domain_decomp_1d, only :
     &     pack_lc=>pack_dataj,unpack_lc=>unpack_dataj
      implicit none
      private

!@param JM_BUDG grid size for budget page diags
!@param IMLON,JMLAT latlon grid sizes
      INTEGER, PARAMETER, public :: JM_BUDG=JM,JMLAT=JM
      INTEGER, PARAMETER, public :: IMLON=IM,IMLONH=IM/2
!@var XWON scale factor for diag. printout needed for Wonderland model
      REAL*8, public :: XWON = TWOPI/(DLON*REAL(IM,KIND=8))

      public :: pack_lc,unpack_lc,get_bounds,get_alloc_bounds

      contains

      subroutine get_bounds(grid,
     &     j_strt_budg,j_stop_budg,
     &     j_strt_jk,j_stop_jk
     &     )
      use domain_decomp_1d, only : dist_grid
      type(dist_grid), intent(in) :: grid
      integer, optional, intent(out) :: j_strt_budg,j_stop_budg
      integer, optional, intent(out) :: j_strt_jk,j_stop_jk
      if(present(j_strt_budg)) j_strt_budg = grid%j_strt
      if(present(j_stop_budg)) j_stop_budg = grid%j_stop
      if(present(j_strt_jk)) j_strt_jk = grid%j_strt
      if(present(j_stop_jk)) j_stop_jk = grid%j_stop
      return
      end subroutine get_bounds

      subroutine get_alloc_bounds(grid,
     &     j_strt_budg,j_stop_budg,
     &     j_strt_jk,j_stop_jk
     &     )
      use domain_decomp_1d, only : dist_grid
      type(dist_grid), intent(in) :: grid
      integer, optional, intent(out) :: j_strt_budg,j_stop_budg
      integer, optional, intent(out) :: j_strt_jk,j_stop_jk
      if(present(j_strt_budg)) j_strt_budg = grid%j_strt_halo
      if(present(j_stop_budg)) j_stop_budg = grid%j_stop_halo
      if(present(j_strt_jk)) j_strt_jk = grid%j_strt_halo
      if(present(j_stop_jk)) j_stop_jk = grid%j_stop_halo
      return
      end subroutine get_alloc_bounds

      end module diag_zonal
