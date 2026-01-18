#include "rundeck_opts.h"
!------------------------------------------------------------------------------
module KochTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  KochTracersMetadata_mod encapsulates the KOCH tracers metadata
!@auth NCCS ASTG
  use sharedTracersMetadata_mod, only: DMS_setspec, &
    SO2_setspec, H2O2_s_setspec
  use sharedTracersMetadata_mod, only: convert_HSTAR
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: nPart, nGAS
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: set_HSTAR
  use OldTracer_mod, only: set_tr_RKD
  use OldTracer_mod, only: set_pm2p5fact
  use OldTracer_mod, only: set_pm10fact
  use OldTracer_mod, only: set_has_chemistry
  use OldTracer_mod, only: tr_RKD 
  use TRACER_COM, only:  n_MSA, n_SO4, n_DMS, &
    n_BCII,  n_BCIA,  n_BCB, n_OCII,  n_OCIA,  n_OCB, n_H2O2_s
  use TRACER_COM, only: n_vbsGm2, n_vbsGm1, n_vbsGz,  n_vbsGp1, n_vbsGp2, &
                        n_vbsGp3, n_vbsGp4, n_vbsGp5, n_vbsGp6, &
                        n_vbsAm2, n_vbsAm1, n_vbsAz,  n_vbsAp1, n_vbsAp2, &
                        n_vbsAp3, n_vbsAp4, n_vbsAp5, n_vbsAp6
  use Dictionary_mod, only: sync_param
  use RunTimeControls_mod, only: tracers_drydep
  use RunTimeControls_mod, only: sulf_only_aerosols
  use RunTimeControls_mod, only: tracers_special_shindell
  use RunTimeControls_mod, only: dynamic_biomass_burning
  use Tracer_mod, only: Tracer
#ifdef TRACERS_AEROSOLS_VBS
  use TRACERS_VBS, only: ivbs_m2,ivbs_m1,ivbs_m0,ivbs_p1,ivbs_p2,ivbs_p3,&
                         ivbs_p4,ivbs_p5,ivbs_p6
#endif /* TRACERS_AEROSOLS_VBS */

  implicit none
  private

  public Koch_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine KOCH_InitMetadata(pTracer)
    use TRACER_COM, only: coupled_chem
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  DMS_setSpec('DMS')
    call  MSA_setSpec('MSA')
    call  SO2_setSpec('SO2')
    call  SO4_setSpec('SO4')
    if (.not. tracers_special_shindell .or. coupled_chem.eq.0) then
      call  H2O2_s_setSpec('H2O2_s')
    end if
    if (.not. sulf_only_aerosols) then
      call  BCII_setSpec('BCII')
      call  BCIA_setSpec('BCIA')
      call  BCB_setSpec('BCB')
#ifdef TRACERS_AEROSOLS_VBS
      call  VBS_setSpec('vbsGm2', ivbs_m2,'igas')
      call  VBS_setSpec('vbsGm1', ivbs_m1,'igas')
      call  VBS_setSpec('vbsGz', ivbs_m0,'igas')
      call  VBS_setSpec('vbsGp1', ivbs_p1,'igas')
      call  VBS_setSpec('vbsGp2', ivbs_p2,'igas')
      call  VBS_setSpec('vbsGp3', ivbs_p3,'igas')
      call  VBS_setSpec('vbsGp4', ivbs_p4,'igas')
      call  VBS_setSpec('vbsGp5', ivbs_p5,'igas')
      call  VBS_setSpec('vbsGp6', ivbs_p6,'igas')

      call  VBS_setSpec('vbsAm2', ivbs_m2,'iaer')
      call  VBS_setSpec('vbsAm1', ivbs_m1,'iaer')
      call  VBS_setSpec('vbsAz', ivbs_m0,'iaer')
      call  VBS_setSpec('vbsAp1', ivbs_p1,'iaer')
      call  VBS_setSpec('vbsAp2', ivbs_p2,'iaer')
      call  VBS_setSpec('vbsAp3', ivbs_p3,'iaer')
      call  VBS_setSpec('vbsAp4', ivbs_p4,'iaer')
      call  VBS_setSpec('vbsAp5', ivbs_p5,'iaer')
      call  VBS_setSpec('vbsAp6', ivbs_p6,'iaer')
#else
      call  OCII_setSpec('OCII')   !Insoluble industrial organic mass
      call  OCIA_setSpec('OCIA')   !Aged industrial organic mass
      call  OCB_setSpec('OCB')     !Biomass organic mass
#endif /* TRACERS_AEROSOLS_VBS */
    end if
     
!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine MSA_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_MSA = n
      call set_ntm_power(n, -13)
      call set_tr_mm(n, 96.d+0) !(H2O2 34;SO2 64)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 5.d-7 ) !m (SO4 3;BC 1;OC 3)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine MSA_setSpec

    subroutine SO4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SO4 = n 
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d+0)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine SO4_setSpec

    subroutine BCII_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BCII = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine BCII_setSpec

    subroutine BCIA_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BCIA = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine BCIA_setSpec

    subroutine BCB_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BCB = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 0.8d0 ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
    end subroutine BCB_setSpec

#ifdef TRACERS_AEROSOLS_VBS
    subroutine VBS_setSpec(name, index, type)
      use OldTracer_mod, only: om2oc, set_om2oc, set_is_VBS_tracer
      use tracers_vbs, only: vbs_tr
      implicit none
      character(len=*), intent(in) :: name
      real*8 :: tmp
      integer, intent(in) :: index
      character(len=4), intent(in) :: type

      n = oldAddTracer(name)
      call set_is_VBS_tracer(n, .true.)

      select case(name)
        case("vbsGm2"); n_vbsGm2 = n
        case("vbsGm1"); n_vbsGm1 = n
        case("vbsGz"); n_vbsGz = n
        case("vbsGp1"); n_vbsGp1 = n
        case("vbsGp2"); n_vbsGp2 = n
        case("vbsGp3"); n_vbsGp3 = n
        case("vbsGp4"); n_vbsGp4 = n
        case("vbsGp5"); n_vbsGp5 = n
        case("vbsGp6"); n_vbsGp6 = n
        case("vbsAm2"); n_vbsAm2 = n
        case("vbsAm1"); n_vbsAm1 = n
        case("vbsAz"); n_vbsAz = n
        case("vbsAp1"); n_vbsAp1 = n
        case("vbsAp2"); n_vbsAp2 = n
        case("vbsAp3"); n_vbsAp3 = n
        case("vbsAp4"); n_vbsAp4 = n
        case("vbsAp5"); n_vbsAp5 = n
        case("vbsAp6"); n_vbsAp6 = n
      end select
      select case (type)
      case ('igas')
        vbs_tr%igas(index) = n
      case ('iaer')
        vbs_tr%iaer(index) = n
      end select

      call set_om2oc(n, 1.4d0)
      tmp = om2oc(n)
      call sync_param(trim(name)//"_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      select case(name)
      case ('vbsGm2', 'vbsGm1', 'vbsGz',  'vbsGp1', 'vbsGp2', &! VBS gas-phase
        'vbsGp3', 'vbsGp4', 'vbsGp5', 'vbsGp6')
        call set_tr_wd_type(n, ngas)
        call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
      case ('vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2', &! VBS aerosol-phase
            'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
        call set_tr_wd_type(n, npart)
        call set_trpdens(n, 1.5d3) !kg/m3
        call set_trradius(n, 3.d-7 ) !m
        call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
        call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
        call set_pm10fact(n, 1.d0) ! fraction that's PM10
        call set_has_chemistry(n, .true.)
      end select
    end subroutine VBS_setSpec
#endif /* TRACERS_AEROSOLS_VBS */

    subroutine OCII_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_OCII = n
      call set_om2oc(n, 1.4d0)
      tmp = om2oc(n)
      call sync_param("OCII_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine OCII_setSpec

    subroutine OCIA_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_OCIA = n
      call set_om2oc(n, 1.4d0)
      tmp = om2oc(n)
      call sync_param("OCIA_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine OCIA_setSpec

    subroutine OCB_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_OCB = n
      call set_om2oc(n, 1.4d0)
      tmp = om2oc(n)
      call sync_param("OCB_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 0.8d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
    end subroutine OCB_setSpec

  end subroutine KOCH_InitMetadata

end module KochTracersMetadata_mod



