#include "rundeck_opts.h"

      subroutine exp_loss_trop(ns,n)
!@sum Decay_Losstracers calculates tropospheric chemistry for st8025 tracer
!@+     by applying a pre-determined chemical loss rate below the tropopause
!@auth Clara Orbe
      USE RESOLUTION, only: im,jm,lm
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds, AM_I_ROOT
      USE ATM_COM, only: LTROPO
      USE MODEL_COM, only: dtsrc
      USE TRACER_COM
      use OldTracer_mod, only: trname
      USE FLUXES, only: tr3Dsource
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      implicit none
      integer n,ns,i,j,l
      REAL :: dcyrtst8025
      REAL :: expdecst8025
      INTEGER :: J_1, J_0, I_0, I_1
      INTEGER :: J_1H, J_0H 

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

       dcyrtst8025 = 4.6296d-7
       expdecst8025 = exp(-dcyrtst8025*dtsrc)

       do l=1,lm
        if (trname(n).eq.'st8025') then
         do j=J_0,J_1
          do i=I_0,I_1
             if (LTROPO(i,j).ge.l) then
               trm(i,j,l,n) = expdecst8025*trm(i,j,l,n)
               trmom(:,i,j,l,n) = expdecst8025*trmom(:,i,j,l,n)
               tr3Dsource(i,j,l,ns,n)=trm(i,j,l,n)
             else
              trm(i,j,l,n) = trm(i,j,l,n)
              trmom(:,i,j,l,n) = trmom(:,i,j,l,n)
              tr3Dsource(i,j,l,ns,n) = 0.d0
             endif
          enddo
          enddo

        endif
       enddo
      
      return
      end subroutine exp_loss_trop

