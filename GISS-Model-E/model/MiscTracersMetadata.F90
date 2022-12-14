#include "rundeck_opts.h"
!------------------------------------------------------------------------------
module MiscTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  MiscTracersMetadata_mod encapsulates the tracers metadata for tracers
!@+    that have NOT been classified as belonging to any particular group.
!@auth NCCS ASTG
  USE CONSTANT, only: mair,mwat,pi,gasc
  use RunTimeControls_mod, only: tracers_ocean
  use RunTimeControls_mod, only: tracers_drydep
  use RunTimeControls_mod, only: tracers_water
  use OldTracer_mod, only: nWater
  use OldTracer_mod, only: set_needtrs
  use OldTracer_mod, only: nPart
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_trw0
  use OldTracer_mod, only: set_ntrocn
  use OldTracer_mod, only: set_isdust
  use OldTracer_mod, only: set_trglac
  use OldTracer_mod, only: set_rc_washt
  use OldTracer_mod, only: set_t_qlimit
  use OldTracer_mod, only: set_trli0
  use OldTracer_mod, only: set_trsi0
  use OldTracer_mod, only: set_tr_h2obych4
  use OldTracer_mod, only: set_trdecay
  use OldTracer_mod, only: set_pm2p5fact
  use OldTracer_mod, only: set_pm10fact
  use OldTracer_mod, only: set_has_chemistry
  use TRACER_COM, only: &
    n_Air, n_water, n_H2O18, n_H2O17, n_HDO, n_HTO, n_Pb210,n_Be7, &
    n_Be10, n_CFCn, n_CO2n, n_Age, n_SO4_d1, n_SO4_d2, n_SO4_d3, &
    n_N_d1, n_N_d2, n_N_d3, n_NH3,   n_NH4,   n_NO3p, &
    n_OCocean, n_clay, n_silt1, n_silt2, n_silt3, n_silt4, n_silt5
  implicit none
  integer :: n
 
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

#ifdef TRACERS_SPECIAL_O18
  subroutine H2O18_setSpec(name)
    use tracerconstants_mod, only: h2o18
    use oldtracer_mod, only : set_iso_index, trw0
    implicit none
    character(len=*), intent(in) :: name
    real*8 :: fracls

    n = oldAddTracer(name)
    n_H2O18 = n
    call set_ntm_power(n, -7)
    call set_tr_mm(n, H2O18%molMass)
    call set_needtrs(n,  .true.)
    call set_tr_wd_type(n, nwater)
    call set_iso_index(n, 2)          ! indexing for isotopic fractionation calcs
    call set_trw0(n, 2.228d-3   ) ! SMOW mass ratio of water molecules
    call set_trli0(n, 0.980d0*trw0(n)  ) ! d=-20
    call set_trsi0(n, fracls(n)*trw0(n))
    call set_tr_H2ObyCH4(n, trw0(n)*1.023d0 ) ! d=+23 (ie. no frac from O2)
    call set_ntrocn(n, -3)
#ifdef TRACERS_OCEAN
    if (tracers_ocean) call set_trglac(n, trw0(n)*0.98d0 ) ! d=-20
#endif

  end subroutine H2O18_setSpec

  subroutine HDO_setSpec(name)
    use oldtracer_mod, only : set_iso_index, trw0
    implicit none
    character(len=*), intent(in) :: name
    real*8 :: fracls

    n = oldAddTracer(name)
    n_HDO = n
    call set_ntm_power(n, -8)
    call set_tr_mm(n, 19d0)
    call set_needtrs(n,  .true.)
    call set_tr_wd_type(n, nwater)
    call set_iso_index(n, 3)          ! indexing for isotopic fractionation calcs
    call set_trw0(n, 3.29d-4    ) ! SMOW mass ratio of water molecules
    call set_trli0(n, 0.830d0*trw0(n)  ) ! d=-170
    call set_trsi0(n, fracls(n)*trw0(n))
    call set_tr_H2ObyCH4(n, trw0(n)*0.93d0  ) ! d=-70
    call set_ntrocn(n, -4)
#ifdef TRACERS_OCEAN
    call set_trglac(n, trw0(n)*0.84d0   ) ! d=-160
#endif

  end subroutine HDO_setSpec

  subroutine HTO_setSpec(name)
    use oldtracer_mod, only : set_iso_index
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_HTO = n
    call set_ntm_power(n, -18)
    call set_tr_mm(n, 20d0)
    call set_needtrs(n,  .true.)
    call set_tr_wd_type(n, nwater)
    call set_iso_index(n, 4)          ! indexing for isotopic fractionation calcs
    call set_trw0(n, 0d0)     !2.22d-18   ) ! SMOW mass ratio of water molecules
    call set_trli0(n, 0d0)
    call set_trsi0(n, 0d0)
    call set_tr_H2ObyCH4(n, 0d0)
    call set_trdecay(n,  1.77d-9) ! =5.59d-2 /yr
    call set_ntrocn(n, -18)
#ifdef TRACERS_OCEAN
    if (tracers_ocean) call set_trglac(n, 0d0)
#endif
  end subroutine HTO_setSpec

  subroutine H2O17_setSpec(name)
    use tracer_com, only: n_h2o17
    use oldtracer_mod, only : set_iso_index, trw0
    implicit none
    character(len=*), intent(in) :: name
    real*8 :: fracls

    n = oldAddTracer(name)
    n_H2O17 = n
    call set_ntm_power(n, -7)
    call set_tr_mm(n, 19d0)
    call set_needtrs(n,  .true.)
    call set_tr_wd_type(n, nwater)
    call set_iso_index(n, 5)          ! indexing for isotopic fractionation calcs
    call set_trw0(n, 4.020d-5   ) ! SMOW mass ratio of water molecules
    call set_trli0(n, 0.98937d0*trw0(n)  ) ! d=-10.63 D17O=0
    call set_trsi0(n, fracls(n)*trw0(n))
    call set_tr_H2ObyCH4(n, trw0(n)*1.011596d0 ) ! d=+11.596 (some frac from O2)
    call set_ntrocn(n, -3)
#ifdef TRACERS_OCEAN
    if (tracers_ocean) call set_trglac(n, trw0(n)*0.98937d0   ) ! d=-10.63 D17O=
#endif
  end subroutine H2O17_setSpec

#endif

  subroutine CFCn_setSpec(name)
    use tracer_com, only: gasex_index
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_CFCn = n
    call gasex_index%push_back(n_CFCn)
    call set_ntm_power(n, -12)
    call set_tr_mm(n, 137.37d0) !note units are in gr
    call set_needtrs(n, .true.)

  end subroutine CFCn_setSpec

  subroutine CO2n_setSpec(name)
    use tracer_com, only: gasex_index
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_CO2n = n
    call gasex_index%push_back(n_CO2n)
    call set_ntm_power(n, -6)
    call set_tr_mm(n, 44.d0)  !grams
    call set_t_qlimit(n,  .false.)
    call set_needtrs(n, .true.)

  end subroutine CO2n_setSpec

  subroutine OCocean_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_OCocean = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 15.6d0)
    call set_trpdens(n, 1.5d3) !kg/m3
    call set_trradius(n, 4.4d-7) !m, same as seasalt1
    call set_fq_aer(n, 1.0d0) ! same as seasalt
    call set_tr_wd_type(n, nPART)
    call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10

  end subroutine OCocean_setSpec

  subroutine Clay_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_clay=n
    call set_ntm_power(n, -9)
    call set_trpdens(n, 2.5d3)
    if (tracers_drydep) call set_trradius(n, 0.46D-06)
    call set_fq_aer(n, 5.D-1)
    call set_rc_washt(n, 5.D-1)
    call set_tr_wd_type(n, nPART)
    call set_tr_mm(n, 1.d+0)
    call set_isdust(n, 1)
    call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10

  end subroutine Clay_setSpec

  subroutine Silt1_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name
    n = oldAddTracer(name)
    n_Silt1=n
    call set_ntm_power(n, -9)
    call set_trpdens(n, 2.65d3)
    if (tracers_drydep) call set_trradius(n, 1.47D-06)
    call set_fq_aer(n, 5.D-1)
    call set_rc_washt(n, 5.D-1)
    call set_tr_wd_type(n, nPART)
    call set_tr_mm(n, 1.d+0)
    call set_isdust(n, 1)
    call set_pm2p5fact(n, 0.322d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10
  end subroutine Silt1_setSpec

  subroutine Silt2_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_Silt2=n
    call set_ntm_power(n, -9)
    call set_trpdens(n, 2.65d3)
    if (tracers_drydep) call set_trradius(n, 2.94D-06)
    call set_fq_aer(n, 5.D-1)
    call set_rc_washt(n, 5.D-1)
    call set_tr_wd_type(n, nPART)
    call set_tr_mm(n, 1.d+0)
    call set_isdust(n, 1)
    call set_pm2p5fact(n, 0.d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10

  end subroutine Silt2_setSpec

  subroutine Silt3_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name
    n = oldAddTracer(name)
    n_Silt3 = n
    call set_ntm_power(n, -9)
    call set_trpdens(n, 2.65d3)
    if (tracers_drydep) call set_trradius(n, 5.88D-06)
    call set_fq_aer(n, 5.D-1)
    call set_rc_washt(n, 5.D-1)
    call set_tr_wd_type(n, nPART)
    call set_tr_mm(n, 1.d+0)
    call set_isdust(n, 1)
    call set_pm2p5fact(n, 0.d0) ! fraction that's PM2.5
    call set_pm10fact(n, 0.322d0) ! fraction that's PM10
  end subroutine Silt3_setSpec

  subroutine Silt4_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name
    n = oldAddTracer(name)
    n_Silt4 = n
    call set_ntm_power(n, -9)
    call set_trpdens(n, 2.65d3)
    if (tracers_drydep) call set_trradius(n, 11.77D-06)
    call set_fq_aer(n, 5.D-1)
    call set_rc_washt(n, 5.D-1)
    call set_tr_wd_type(n, nPART)
    call set_tr_mm(n, 1.d+0)
    call set_isdust(n, 1)
  end subroutine Silt4_setSpec

  subroutine Silt5_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name
    n = oldAddTracer(name)
    n_Silt5 = n
    call set_ntm_power(n, -9)
    call set_trpdens(n, 2.65d3)
    if (tracers_drydep) call set_trradius(n, 23.53D-06)
    call set_fq_aer(n, 5.D-1)
    call set_rc_washt(n, 5.D-1)
    call set_tr_wd_type(n, nPART)
    call set_tr_mm(n, 1.d+0)
    call set_isdust(n, 1)
  end subroutine Silt5_setSpec

  subroutine NO3p_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name
    n = oldAddTracer(name)
    n_NO3p = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 62.d0)
    call set_trpdens(n, 1.7d3)
    call set_trradius(n, 3.d-7)
    call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart)
    call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10
    call set_has_chemistry(n, .true.)
  end subroutine NO3p_setSpec

  subroutine SO4_d1_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_SO4_d1 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 96.d0)  !!!! Sulfat
    call set_trpdens(n, 2.5d3) !kg/m3 this is clay density
    call set_trradius(n, 0.75D-06 ) !m
    call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart)
    call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10
    call set_has_chemistry(n, .true.)
  end subroutine SO4_d1_setSpec

  subroutine SO4_d2_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_SO4_d2 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 96.d0)
    call set_trpdens(n, 2.65d3) !kg/m3 this is Silt1 value
    call set_trradius(n, 2.2D-06 ) !m
    call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart)
    call set_pm2p5fact(n, 0.322d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10
    call set_has_chemistry(n, .true.)
  end subroutine SO4_d2_setSpec

  subroutine SO4_d3_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_SO4_d3 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 96.d0)
    call set_trpdens(n, 2.65d3) !this is Silt2 value
    call set_trradius(n, 4.4D-06 ) !m this is Silt2 value
    call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart)
    call set_pm2p5fact(n, 0.d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10
    call set_has_chemistry(n, .true.)
  end subroutine SO4_d3_setSpec

  subroutine N_d1_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name
    n = oldAddTracer(name)
    n_N_d1 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 62.d+0) ! NO3
    call set_trpdens(n, 2.5d3) !kg/m3 this is clay density
    call set_trradius(n, 0.75D-06 ) !m
    call set_fq_aer(n, 1.d0  ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart)
    call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10
    call set_has_chemistry(n, .true.)
  end subroutine N_d1_setSpec

  subroutine N_d2_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_N_d2 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 62.d+0)
    call set_trpdens(n, 2.65d3) !kg/m3 this is Silt1 value
    call set_trradius(n, 2.2D-06 ) !m
    call set_fq_aer(n, 1.d0  ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart)
    call set_pm2p5fact(n, 0.322d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10
    call set_has_chemistry(n, .true.)
  end subroutine N_d2_setSpec

  subroutine N_d3_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_N_d3 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 62.d0)
    call set_trpdens(n, 2.65d3) !this is Silt2 value
    call set_trradius(n, 4.4D-06 ) !m this is Silt2 value
    call set_fq_aer(n, 1.d0  ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart)
    call set_pm2p5fact(n, 0.d0) ! fraction that's PM2.5
    call set_pm10fact(n, 1.d0) ! fraction that's PM10
    call set_has_chemistry(n, .true.)
  end subroutine N_d3_setSpec

  subroutine Pb210_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_Pb210 = n
    call set_ntm_power(n, -23)
    call set_tr_mm(n, 210.d0)
    call set_trdecay(n,  9.85d-10)
    call set_trpdens(n, 1.7d3) !kg/m3 this is SO4 value
    call set_trradius(n, 3.d-7  ) !again S04 value
    call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart) ! same as SO4
    call set_has_chemistry(n, .true.)
  end subroutine Pb210_setSpec

  subroutine Be7_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_Be7 = n
    call set_ntm_power(n, -23) ! power of ten for tracer
    call set_tr_mm(n, 7.d0)
    call set_trdecay(n,  1.51d-7)
    call set_trpdens(n, 1.7d3) !kg/m3 this is SO4 value
    call set_trradius(n, 1.d-7  ) !appropriate for stratosphere
    call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart) ! same as SO4
    call set_has_chemistry(n, .true.)
  end subroutine Be7_setSpec

  subroutine Be10_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_Be10 = n
    call set_ntm_power(n, -23)
    call set_tr_mm(n, 10.d0)
    call set_trpdens(n, 1.7d3) !kg/m3 this is SO4 value
    call set_trradius(n, 1.d-7  ) !appropriate for stratosphere
    call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
    call set_tr_wd_type(n, npart) ! same as SO4
    call set_has_chemistry(n, .true.)
  end subroutine Be10_setSpec


  subroutine Air_setSpec(name)
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_Air = n
    call set_ntm_power(n, -2)
    call set_tr_mm(n, mair)

  end subroutine Air_setSpec

  subroutine Water_setSpec(name)
    use oldtracer_mod, only : set_iso_index
    use RunTimeControls_mod, only: tracers_special_o18
    implicit none
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_Water = n
    if (tracers_water .or. tracers_ocean) then
      call set_trw0(n, 1.d+0)
      call set_ntrocn(n, 0)
    end if
    if (tracers_water) then
      call set_ntm_power(n, -4)
      call set_tr_mm(n, mwat)
      call set_needtrs(n,  .true.)
      call set_tr_wd_type(n, nWater)
      call set_trli0(n, 1.d+0)
      call set_trsi0(n, 1.d+0)
      call set_tr_H2ObyCH4(n, 1.d+0)
    end if
    if (tracers_ocean) call set_trglac(n, 1.d+0)
#ifdef TRACERS_SPECIAL_O18
    if (tracers_special_o18) call set_iso_index(n, 1) ! indexing for isotopic fractionation calcs
#endif
  end subroutine Water_setSpec

end module MiscTracersMetadata_mod
