#include "rundeck_opts.h"
!------------------------------------------------------------------------------
module sharedTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  sharedTracersMetadata_mod encapsulates the metadata shared among various
!@+    tracers.
!@auth NCCS ASTG
  use OldTracer_mod, only: nPart
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_HSTAR
  use OldTracer_mod, only: set_F0
  use OldTracer_mod, only: set_tr_RKD
  use OldTracer_mod, only: set_tr_DHD
  use OldTracer_mod, only: set_trdecay
  use OldTracer_mod, only: tr_RKD 
  use OldTracer_mod, only: set_needtrs
#ifdef TRACERS_SPECIAL_Lerner
  use TRACERS_MPchem_COM, only: nMPtable
  use OldTracer_mod, only: set_iMPtable
  use OldTracer_mod, only: set_tcscale
#endif  /* TRACERS_SPECIAL_Lerner */
  use OldTracer_mod, only: dodrydep
  use OldTracer_mod, only: F0
  use OldTracer_mod, only: HSTAR
  use OldTracer_mod, only: ngas, nPART
  use OldTracer_mod, only: set_pm2p5fact
  use OldTracer_mod, only: set_pm10fact
  use OldTracer_mod, only: set_has_chemistry
  use TRACER_COM, only : seasonalNH3src
  use TRACER_COM, only: n_H2O2, n_NH3,  n_NH4, n_DMS, n_SO2, n_H2O2_s, &
    n_CH4, n_N2O, n_Rn222
  use Dictionary_mod, only: sync_param
  use RunTimeControls_mod, only: tracers_drydep
  use RunTimeControls_mod, only: tracers_special_lerner
  use RunTimeControls_mod, only: dynamic_biomass_burning  
  implicit none
  private 

  public DMS_setSpec
  public SO2_setSpec
  public H2O2_setSpec
  public NH3_setSpec
  public NH4_setSpec
  public H2O2_s_setSpec
  public CH4_setSpec
  public N2O_setSpec
  public Rn222_setSpec

!@param convert_HSTAR converts from mole/Joule to mole/(L*atm)
  real(8), parameter :: convert_HSTAR = 1.01325d2
  public convert_HSTAR

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  subroutine DMS_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_DMS = n
    call set_ntm_power(n, -12)
    call set_tr_mm(n, 62.d+0)
    call set_needtrs(n, .true.)
    call set_has_chemistry(n, .true.)
  end subroutine DMS_setSpec

  subroutine SO2_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_SO2 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 64.d+0)
    call set_tr_RKD(n, 0.0118d0 ) !mole/J or  1.2  M/atm
    call set_tr_DHD(n, -2.62d4) ! in J/mole= -6.27 kcal/mol
    call set_tr_wd_type(n, ngas)
    if (tracers_drydep) CALL SET_HSTAR(N, 1.D5)
    call set_has_chemistry(n, .true.)
  end subroutine SO2_setSpec

  subroutine H2O2_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_H2O2 = n
    call set_ntm_power(n, -11)
    call set_tr_mm(n, 34.016d0)
    call set_tr_RKD(n, 9.869d2    ) ! in mole/J = 1.d5 mole/(L atm)
    call set_tr_DHD(n, -5.52288d4 ) ! in J/mole = -13.2 kcal/mole.
    if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    call set_F0(n,  1.d0)
    call set_has_chemistry(n, .true.)
  end subroutine H2O2_setSpec

  subroutine NH3_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_NH3 = n
    call set_ntm_power(n, -10)
    call set_tr_mm(n, 17.d0)
    call set_tr_RKD(n, 100.d0) ! higher than nominal; effective Henry
!    call set_tr_RKD(n, 0.7303d0   ) !tr_RKD=74 M/atm
    call set_tr_DHD(n, -2.84d4  ) !tr_DHD=-6.80 kcal/mole
    call set_tr_wd_type(n, ngas)
    if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    call sync_param("seasonalNH3src", seasonalNH3src)
    call set_has_chemistry(n, .true.)
  end subroutine NH3_setSpec

  subroutine H2O2_s_setSpec(name)
    character(len=*), intent(in) :: name

    n = oldAddTracer(name)
    n_H2O2_s = n
    call set_ntm_power(n, -10)
    call set_tr_mm(n, 34.016d0)
    call set_tr_RKD(n, 986.9d0)
    call set_tr_DHD(n, -5.52288d4 ) ! in J/mole = -13.2 kcal/mole.
    call set_tr_wd_type(n, ngas)
    if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    call set_F0(n,  1.d0)
    call set_has_chemistry(n, .true.)
  end subroutine H2O2_s_setSpec

    subroutine CH4_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_CH4 = n
      call set_tr_mm(n, 16.d0)
#ifdef TRACERS_SPECIAL_Lerner
      if (tracers_special_lerner) then
        call set_ntm_power(n, -9)
        nMPtable=nMPtable+1
        call set_iMPtable(n, nMPtable)
        call set_tcscale(n, 1.d0)
      end if
#endif
      call set_ntm_power(n, -8)

    call set_has_chemistry(n, .true.)
    end subroutine CH4_setSpec

    subroutine N2O_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_N2O = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, 44.d0)
#ifdef TRACERS_SPECIAL_Lerner
      if (tracers_special_lerner) then
        nMPtable=nMPtable+1
        call set_iMPtable(n, nMPtable)
        call set_tcscale(n, 1.d0)
      end if
#endif
      call set_has_chemistry(n, .true.)
    end subroutine N2O_setSpec

    subroutine Rn222_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Rn222 = n 
      call set_ntm_power(n, -21)
      call set_tr_mm(n, 222.d0)
      call set_trdecay(n,  2.1d-6)
      call set_has_chemistry(n, .true.)
    end subroutine Rn222_setSpec

    subroutine NH4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_NH4 = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 18.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine NH4_setSpec


  ! TOMAS duplicates with nitrate: NH3_setSpec, NH4_setSpec

end module sharedTracersMetadata_mod
