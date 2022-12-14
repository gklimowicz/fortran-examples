#include "rundeck_opts.h"
!!
!@sum  OCEANR_DIM  for ocean grid
!@auth Larissa Nazarenko
      MODULE OCEANR_DIM
      USE DOMAIN_DECOMP_1D, only : DIST_GRID
      USE OCEANRES,  only : imo,jmo, lmo

      implicit none

      private
      save
      public init_oceanr_grid  

      public I_0,  I_1,  J_0,  J_1, I_0H, I_1H, J_0H, J_1H
      public ogrid

c
      TYPE(DIST_GRID), target :: ogrid   ! ocean (Russell) grid
      ! domain bounds
      integer :: I_0,  I_1,  J_0,  J_1
      ! domain bounds with halos
      integer :: I_0H, I_1H, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      
      contains

      subroutine init_oceanr_grid
      USE DOMAIN_DECOMP_1D, only : init_grid, getDomainBounds
 
!      TYPE(DIST_GRID) :: ogrid   ! ocean (Russell) grid

      call init_grid( ogrid, imo, jmo, lmo )
      call getDomainBounds(ogrid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H ,
     &               I_STRT     =I_0,    I_STOP     =I_1,
     &               I_STRT_HALO=I_0H,   I_STOP_HALO=I_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
 
!      write (555,*) ' init_oceanr_grid: i0,i1,j0,j1', I_0,I_1,J_0,J_1
!      write (555,*) ' init_oceanr_grid: halos', I_0H,I_1H,J_0H,J_1H

      end subroutine init_oceanr_grid

      end module oceanr_dim
