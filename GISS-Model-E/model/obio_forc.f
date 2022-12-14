#include "rundeck_opts.h"

      MODULE obio_forc

      USE timestream_mod, only : timestream
      implicit none


      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tirrq3d         !total mean irradiance in quanta
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: avgq            !mean daily irradiance in quanta
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: atmFe           !atm Fe in nM
      real, ALLOCATABLE, DIMENSION(:,:)  :: surfN
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: alk             !alkalinity in 'umol/kg'
#ifdef prescribe_o2sf
      real, ALLOCATABLE, DIMENSION(:,:,:) :: o2fc            !@PL array for prescribed surface ocean O2
#endif

      real solz               !mean cosine solar zenith angle
      real sunz               !solar zenith angle
      real, allocatable, dimension(:) ::  Ed, Es
      real wind               !surface wind from atmos
      real, allocatable, dimension(:) :: tirrq !total mean irradiance in quanta
      real, parameter ::  tirrq_critical=10. !in quanta threshold at compensation depth
      real rmud               !downwelling irradiance average cosine
      real rhosrf             !surface air density which comes from PBL.f
#ifdef OBIO_RUNOFF
      real river_runoff       
#endif

#ifdef STANDALONE_OCEAN
      real, ALLOCATABLE, DIMENSION(:,:,:,:,:):: Eda_glob,Esa_glob       !direct,diffuse downwelling irradiance
      real, ALLOCATABLE, DIMENSION(:,:,:,:,:):: Eda,Esa
      real, ALLOCATABLE, DIMENSION(:,:):: Eda2,Esa2
#endif
      type(timestream) :: stream_atmFe

      END MODULE obio_forc

!------------------------------------------------------------------------------
!NOT FOR HYCOM:idm and jdm were passed to the subroutine       
      subroutine alloc_obio_forc(kdm,ogrid,idm,jdm)

      use ocalbedo_mod, only: nlt
      USE obio_forc
      USE DOMAIN_DECOMP_1D, only :DIST_GRID



      implicit none

      integer, intent(in) :: kdm,idm,jdm
      type(DIST_grid), intent(in) :: ogrid

      INTEGER :: j_0,j_1,i_0,i_1

      I_0 = ogrid%I_STRT
      I_1 = ogrid%I_STOP
      J_0 = ogrid%J_STRT
      J_1 = ogrid%J_STOP

      ALLOCATE(tirrq3d(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(   avgq(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(    alk(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(atmFe(i_0:i_1,j_0:j_1,12))
      allocate(tirrq(kdm))
      allocate(Ed(nlt),Es(nlt))
#ifdef STANDALONE_OCEAN
      ALLOCATE(Eda_glob(idm,jdm,nlt,12,12),Esa_glob(idm,jdm,nlt,12,12))
      ALLOCATE(Eda2(nlt,12),Esa2(nlt,12))
      ALLOCATE (Eda(i_0:i_1,j_0:j_1,nlt,12,12))
      ALLOCATE (Esa(i_0:i_1,j_0:j_1,nlt,12,12))
#endif

      end subroutine alloc_obio_forc
