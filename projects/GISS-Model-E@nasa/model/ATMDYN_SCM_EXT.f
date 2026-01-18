#include "rundeck_opts.h"
#define JJ(J) (J)-J_0H+1
!
!dummy routines for running SCM   --   A. Wolf
!
      function getTotalEnergy() result(totalEnergy)
!@sum  getTotalEnergy Dummy
!@auth Tom Clune (SIVO)
!
      REAL*8 :: totalEnergy

      totalEnergy = 0.

      return
      end function getTotalEnergy


      subroutine addEnergyAsDiffuseHeat(deltaEnergy)
!@sum  addEnergyAsDiffuseHeat Dummy
!@auth Tom Clune (SIVO)
      real*8, intent(in) :: deltaEnergy

      return
      end subroutine addEnergyAsDiffuseHeat

      SUBROUTINE DISSIP
!@sum DISSIP adds in dissipated KE (m^2/s^2) as heat locally
!@auth Gavin Schmidt
      USE ATM_COM, only : dke,kea,t,pk
      IMPLICIT NONE

      return

      END SUBROUTINE DISSIP

C***** Add in dissipiated KE as heat locally
      subroutine addEnergyAsLocalHeat(deltaKE, T, PK, diagIndex)
!@sum  addEnergyAsLocalHeat Dummy
!@auth Tom Clune (SIVO)
      use DOMAIN_DECOMP_ATM, only: grid
      use DOMAIN_DECOMP_1D, only: getDomainBounds
      implicit none
      real*8 :: deltaKE(:,grid%j_strt_halo:,:)
      real*8 :: T(:,grid%j_strt_halo:,:)
      real*8 :: PK(:,:,grid%j_strt_halo:)
      integer, optional, intent(in) :: diagIndex

      return
      end subroutine addEnergyAsLocalHeat

      subroutine COMPUTE_WSAVE ! vertical velocity for diagnostics
!@sum COMPUTE_WSAVE Dummy

      return

      end subroutine COMPUTE_WSAVE
