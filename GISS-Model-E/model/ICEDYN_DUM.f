#include "rundeck_opts.h"

      MODULE ICEDYN_COM
      USE EXCHANGE_TYPES, only : iceocn_xchng_vars
      integer :: imic=0,kticij=0
      real*8,dimension(:),allocatable :: icij
      type(iceocn_xchng_vars) :: igice
      END MODULE ICEDYN_COM

      MODULE ICEDYN
      integer :: imicdyn=0
      END MODULE ICEDYN

      SUBROUTINE ICEDYN_DUM
!@sum ICEDYN_DUM dummy routines to replace ice dynamics
      !ENTRY alloc_icedyn
      ENTRY alloc_icedyn_com
      ENTRY gather_icdiags
      ENTRY io_icedyn
      ENTRY io_icdiag
      ENTRY reset_icdiag
      ENTRY diag_ICEDYN
      entry def_rsf_icedyn
      entry new_io_icedyn
      entry def_rsf_icdiag
      entry new_io_icdiag
      entry set_ioptrs_iceacc_default
      entry def_meta_icdiag
      entry write_meta_icdiag
      RETURN
      END SUBROUTINE ICEDYN_DUM

      subroutine alloc_icedyn(im,jm)
      use ICEDYN_COM
      integer :: im,jm
      allocate( icij(2) )
      icij=0
      return
      end subroutine alloc_icedyn

      SUBROUTINE DYNSI(atmice,iceocn,si_ocn)
!@sum DYNSI simple coding to estimate ice-ocean friction velocity
!@auth Gavin Schmidt
      USE CONSTANT, only : rhows
      USE RESOLUTION, only : im,jm
      USE FLUXES, only : focean
      USE MODEL_COM, only : kocean,dtsrc
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds
      USE GEOM, only : imaxj
      USE SEAICE, only : oi_ustar0
      USE EXCHANGE_TYPES, only : atmice_xchng_vars,iceocn_xchng_vars
      USE SEAICE_COM, only : icestate
      IMPLICIT NONE
      type(atmice_xchng_vars) :: atmice
      type(iceocn_xchng_vars) :: iceocn
      type(icestate) :: si_ocn
c
      INTEGER I,J
      REAL*8 ustar1
      INTEGER :: I_0,I_1, J_0,J_1
      real*8, dimension(:,:), pointer :: rsi,ui2rho,dmua,dmva

      rsi => si_ocn%rsi
      ui2rho => iceocn%ui2rho
      dmua => atmice%dmua
      dmva => atmice%dmva

      IF (KOCEAN.eq.1) THEN

        call getDomainBounds(grid, J_STRT=J_0,   J_STOP=J_1 )
        I_0 = GRID%I_STRT
        I_1 = GRID%I_STOP
        DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
c          UI2rho(I,J) = rhows*(oi_ustar0)**2  ! default
C**** with wind stress dependence
          if (rsi(i,j)*focean(i,j).gt.0) then
            ustar1= SQRT(SQRT(DMUA(I,J)**2+DMVA(I,J)**2)/
     &           (DTSRC*RHOWS))
            UI2rho(I,J)=rhows*(oi_ustar0*max(1d0,1d3*ustar1))**2
          end if
        END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE DYNSI

      SUBROUTINE ADVSI(atmice)
      USE EXCHANGE_TYPES, only : atmice_xchng_vars
      IMPLICIT NONE
      type(atmice_xchng_vars) :: atmice
      RETURN
      END SUBROUTINE ADVSI

      SUBROUTINE init_ICEDYN(iniOCEAN,atmice)
      USE EXCHANGE_TYPES, only : atmice_xchng_vars
      IMPLICIT NONE
      integer :: iniOCEAN
      type(atmice_xchng_vars) :: atmice
      RETURN
      END SUBROUTINE init_ICEDYN
