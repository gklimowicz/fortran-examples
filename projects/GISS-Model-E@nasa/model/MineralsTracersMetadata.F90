!------------------------------------------------------------------------------
module MineralsTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  MineralsTracersMetadata_mod encapsulates the Minerals tracers metadata
!@auth NCCS ASTG, modified by Jan Perlwitz
  use sharedTracersMetadata_mod
  use TRACER_COM, only: &
    n_clayilli, n_claykaol, n_claysmec, n_claycalc, &
    n_clayquar, n_clayfeld, n_clayhema, n_claygyps, &
    n_clayilhe, n_claykahe, n_claysmhe, n_claycahe, &
    n_clayquhe, n_clayfehe, n_claygyhe,             &
    n_sil1illi, n_sil1kaol, n_sil1smec, n_sil1calc, &
    n_sil1quar, n_sil1feld, n_sil1hema, n_sil1gyps, &
    n_sil1ilhe, n_sil1kahe, n_sil1smhe, n_sil1cahe, &
    n_sil1quhe, n_sil1fehe, n_sil1gyhe,             &
    n_sil2illi, n_sil2kaol, n_sil2smec, n_sil2calc, &
    n_sil2quar, n_sil2feld, n_sil2hema, n_sil2gyps, &
    n_sil2ilhe, n_sil2kahe, n_sil2smhe, n_sil2cahe, &
    n_sil2quhe, n_sil2fehe, n_sil2gyhe,             &
    n_sil3illi, n_sil3kaol, n_sil3smec, n_sil3calc, &
    n_sil3quar, n_sil3feld, n_sil3hema, n_sil3gyps, &
    n_sil3ilhe, n_sil3kahe, n_sil3smhe, n_sil3cahe, &
    n_sil3quhe, n_sil3fehe, n_sil3gyhe,             &
    n_sil4illi, n_sil4kaol, n_sil4smec, n_sil4calc, &
    n_sil4quar, n_sil4feld, n_sil4hema, n_sil4gyps, &
    n_sil4ilhe, n_sil4kahe, n_sil4smhe, n_sil4cahe, &
    n_sil4quhe, n_sil4fehe, n_sil4gyhe,             &
    n_sil5illi, n_sil5kaol, n_sil5smec, n_sil5calc, &
    n_sil5quar, n_sil5feld, n_sil5hema, n_sil5gyps, &
    n_sil5ilhe, n_sil5kahe, n_sil5smhe, n_sil5cahe, &
    n_sil5quhe, n_sil5fehe, n_sil5gyhe
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_tr_mm, set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_tr_wd_TYPE
  use OldTracer_mod, only: set_isDust
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_rc_washt
  use OldTracer_mod, only: nPart
  use RunTimeControls_mod, only: tracers_dust_silt4
  use RunTimeControls_mod, only: tracers_dust_silt5
  use RunTimeControls_mod, only: tracers_drydep
  use trdust_mod, only : densityIllite, densityKaolinite, densitySmectite, &
       densityCalcite, densityQuartz, densityFeldspar, densityHematite, &
       densityGypsum, frIronOxideInAggregate
  use Tracer_mod, only: Tracer

  implicit none
  private

  public Minerals_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine Minerals_InitMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  ClayIlli_setSpec('ClayIlli')
    call  ClayKaol_setSpec('ClayKaol')
    call  ClaySmec_setSpec('ClaySmec')
    call  ClayCalc_setSpec('ClayCalc')
    call  ClayQuar_setSpec('ClayQuar')
    call  ClayFeld_setSpec('ClayFeld')
    call  ClayHema_setSpec('ClayHema')
    call  ClayGyps_setSpec('ClayGyps')
    call  ClayIlHe_setSpec('ClayIlHe')
    call  ClayKaHe_setSpec('ClayKaHe')
    call  ClaySmHe_setSpec('ClaySmHe')
    call  ClayCaHe_setSpec('ClayCaHe')
    call  ClayQuHe_setSpec('ClayQuHe')
    call  ClayFeHe_setSpec('ClayFeHe')
    call  ClayGyHe_setSpec('ClayGyHe')
    call  Sil1Illi_setSpec('Sil1Illi')
    call  Sil1Kaol_setSpec('Sil1Kaol')
    call  Sil1Smec_setSpec('Sil1Smec')
    call  Sil1Calc_setSpec('Sil1Calc')
    call  Sil1Quar_setSpec('Sil1Quar')
    call  Sil1Feld_setSpec('Sil1Feld')
    call  Sil1Hema_setSpec('Sil1Hema')
    call  Sil1Gyps_setSpec('Sil1Gyps')
    call  Sil1IlHe_setSpec('Sil1IlHe')
    call  Sil1KaHe_setSpec('Sil1KaHe')
    call  Sil1SmHe_setSpec('Sil1SmHe')
    call  Sil1CaHe_setSpec('Sil1CaHe')
    call  Sil1QuHe_setSpec('Sil1QuHe')
    call  Sil1FeHe_setSpec('Sil1FeHe')
    call  Sil1GyHe_setSpec('Sil1GyHe')
    call  Sil2Illi_setSpec('Sil2Illi')
    call  Sil2Kaol_setSpec('Sil2Kaol')
    call  Sil2Smec_setSpec('Sil2Smec')
    call  Sil2Calc_setSpec('Sil2Calc')
    call  Sil2Quar_setSpec('Sil2Quar')
    call  Sil2Feld_setSpec('Sil2Feld')
    call  Sil2Hema_setSpec('Sil2Hema')
    call  Sil2Gyps_setSpec('Sil2Gyps')
    call  Sil2IlHe_setSpec('Sil2IlHe')
    call  Sil2KaHe_setSpec('Sil2KaHe')
    call  Sil2SmHe_setSpec('Sil2SmHe')
    call  Sil2CaHe_setSpec('Sil2CaHe')
    call  Sil2QuHe_setSpec('Sil2QuHe')
    call  Sil2FeHe_setSpec('Sil2FeHe')
    call  Sil2GyHe_setSpec('Sil2GyHe')
    call  Sil3Illi_setSpec('Sil3Illi')
    call  Sil3Kaol_setSpec('Sil3Kaol')
    call  Sil3Smec_setSpec('Sil3Smec')
    call  Sil3Calc_setSpec('Sil3Calc')
    call  Sil3Quar_setSpec('Sil3Quar')
    call  Sil3Feld_setSpec('Sil3Feld')
    call  Sil3Hema_setSpec('Sil3Hema')
    call  Sil3Gyps_setSpec('Sil3Gyps')
    call  Sil3IlHe_setSpec('Sil3IlHe')
    call  Sil3KaHe_setSpec('Sil3KaHe')
    call  Sil3SmHe_setSpec('Sil3SmHe')
    call  Sil3CaHe_setSpec('Sil3CaHe')
    call  Sil3QuHe_setSpec('Sil3QuHe')
    call  Sil3FeHe_setSpec('Sil3FeHe')
    call  Sil3GyHe_setSpec('Sil3GyHe')
    if ( tracers_dust_silt4 ) then
      call  Sil4Illi_setSpec('Sil4Illi')
      call  Sil4Kaol_setSpec('Sil4Kaol')
      call  Sil4Smec_setSpec('Sil4Smec')
      call  Sil4Calc_setSpec('Sil4Calc')
      call  Sil4Quar_setSpec('Sil4Quar')
      call  Sil4Feld_setSpec('Sil4Feld')
      call  Sil4Hema_setSpec('Sil4Hema')
      call  Sil4Gyps_setSpec('Sil4Gyps')
      call  Sil4IlHe_setSpec('Sil4IlHe')
      call  Sil4KaHe_setSpec('Sil4KaHe')
      call  Sil4SmHe_setSpec('Sil4SmHe')
      call  Sil4CaHe_setSpec('Sil4CaHe')
      call  Sil4QuHe_setSpec('Sil4QuHe')
      call  Sil4FeHe_setSpec('Sil4FeHe')
      call  Sil4GyHe_setSpec('Sil4GyHe')
    end if
    if ( tracers_dust_silt5 ) then
      call  Sil5Illi_setSpec('Sil5Illi')
      call  Sil5Kaol_setSpec('Sil5Kaol')
      call  Sil5Smec_setSpec('Sil5Smec')
      call  Sil5Calc_setSpec('Sil5Calc')
      call  Sil5Quar_setSpec('Sil5Quar')
      call  Sil5Feld_setSpec('Sil5Feld')
      call  Sil5Hema_setSpec('Sil5Hema')
      call  Sil5Gyps_setSpec('Sil5Gyps')
      call  Sil5IlHe_setSpec('Sil5IlHe')
      call  Sil5KaHe_setSpec('Sil5KaHe')
      call  Sil5SmHe_setSpec('Sil5SmHe')
      call  Sil5CaHe_setSpec('Sil5CaHe')
      call  Sil5QuHe_setSpec('Sil5QuHe')
      call  Sil5FeHe_setSpec('Sil5FeHe')
      call  Sil5GyHe_setSpec('Sil5GyHe')
    end if

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine ClayIlli_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_clayilli = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityIllite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

    end subroutine ClayIlli_setSpec

    subroutine ClayKaol_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_claykaol = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityKaolinite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

    end subroutine ClayKaol_setSpec

    subroutine ClaySmec_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Claysmec = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densitySmectite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

    end subroutine ClaySmec_setSpec

    subroutine ClayCalc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_claycalc = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityCalcite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayCalc_setSpec

    subroutine ClayQuar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_clayquar = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityQuartz)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayQuar_setSpec

    subroutine ClayFeld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_clayfeld = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityFeldspar)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayFeld_setSpec

    subroutine ClayHema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_clayhema = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityHematite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayHema_setSpec

    subroutine ClayGyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_claygyps = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityGypsum)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayGyps_setSpec

    subroutine ClayIlHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_clayilhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) * densityIllite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayIlHe_setSpec

    subroutine ClayKaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_claykahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) *  densityKaolinite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayKaHe_setSpec

    subroutine ClaySmHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_claysmhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densitySmectite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClaySmHe_setSpec

    subroutine ClayCaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_claycahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityCalcite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayCaHe_setSpec

    subroutine ClayQuHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_clayquhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * DensityQuartz + &
           frIronOxideInAggregate * DensityHematite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayQuHe_setSpec

    subroutine ClayFeHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_clayfehe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityFeldspar + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayFeHe_setSpec

    subroutine ClayGyHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_claygyhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityGypsum + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine ClayGyHe_setSpec

    subroutine Sil1Illi_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1illi = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityIllite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1Illi_setSpec

    subroutine Sil1Kaol_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1kaol = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityKaolinite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1Kaol_setSpec

    subroutine Sil1Smec_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1smec = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densitySmectite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1Smec_setSpec

    subroutine Sil1Calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1calc = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityCalcite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1Calc_setSpec

    subroutine Sil1Quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1quar = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityQuartz)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1Quar_setSpec

    subroutine Sil1Feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1feld = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityFeldspar)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1Feld_setSpec

    subroutine Sil1Hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1hema = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1Hema_setSpec

    subroutine Sil1Gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1gyps = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityGypsum)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1Gyps_setSpec

    subroutine Sil1IlHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1ilhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) * densityIllite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1IlHe_setSpec

    subroutine Sil1KaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1kahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) *  densityKaolinite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1KaHe_setSpec

    subroutine Sil1SmHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1smhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densitySmectite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1SmHe_setSpec

    subroutine Sil1CaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1cahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityCalcite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1CaHe_setSpec

    subroutine Sil1QuHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1quhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * DensityQuartz + &
           frIronOxideInAggregate * DensityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1QuHe_setSpec

    subroutine Sil1FeHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1fehe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityFeldspar + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1FeHe_setSpec

    subroutine Sil1GyHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil1gyhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityGypsum + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1GyHe_setSpec

    subroutine Sil2Illi_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2illi = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityIllite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2Illi_setSpec

    subroutine Sil2Kaol_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2kaol = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityKaolinite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2Kaol_setSpec

    subroutine Sil2Smec_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2smec = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densitySmectite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2Smec_setSpec

    subroutine Sil2Calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2calc = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityCalcite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2Calc_setSpec

    subroutine Sil2Quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2quar = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityQuartz)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2Quar_setSpec

    subroutine Sil2Feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2feld = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityFeldspar)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2Feld_setSpec

    subroutine Sil2Hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2hema = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2Hema_setSpec

    subroutine Sil2Gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2gyps = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityGypsum)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2Gyps_setSpec

    subroutine Sil2IlHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2ilhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) * densityIllite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2IlHe_setSpec

    subroutine Sil2KaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2kahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) *  densityKaolinite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2KaHe_setSpec

    subroutine Sil2SmHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2smhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densitySmectite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2SmHe_setSpec

    subroutine Sil2CaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2cahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityCalcite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2CaHe_setSpec

    subroutine Sil2QuHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2quhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * DensityQuartz + &
           frIronOxideInAggregate * DensityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2QuHe_setSpec

    subroutine Sil2FeHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2fehe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityFeldspar + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2FeHe_setSpec

    subroutine Sil2GyHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil2gyhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityGypsum + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2GyHe_setSpec

    subroutine Sil3Illi_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3illi = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityIllite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3Illi_setSpec

    subroutine Sil3Kaol_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3kaol = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityKaolinite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3Kaol_setSpec

    subroutine Sil3Smec_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3smec = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densitySmectite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3Smec_setSpec

    subroutine Sil3Calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3calc = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityCalcite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3Calc_setSpec

    subroutine Sil3Quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3quar = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityQuartz)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3Quar_setSpec

    subroutine Sil3Feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3feld = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityFeldspar)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3Feld_setSpec

    subroutine Sil3Hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3hema = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3Hema_setSpec

    subroutine Sil3Gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3gyps = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityGypsum)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3Gyps_setSpec

    subroutine Sil3IlHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3ilhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) * densityIllite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3IlHe_setSpec

    subroutine Sil3KaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3kahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) *  densityKaolinite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3KaHe_setSpec

    subroutine Sil3SmHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3smhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densitySmectite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3SmHe_setSpec

    subroutine Sil3CaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3cahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityCalcite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3CaHe_setSpec

    subroutine Sil3QuHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3quhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * DensityQuartz + &
           frIronOxideInAggregate * DensityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3QuHe_setSpec

    subroutine Sil3FeHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3fehe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityFeldspar + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3FeHe_setSpec

    subroutine Sil3GyHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil3gyhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityGypsum + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3GyHe_setSpec

    subroutine Sil4Illi_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4illi = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityIllite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4Illi_setSpec

    subroutine Sil4Kaol_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4kaol = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityKaolinite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4Kaol_setSpec

    subroutine Sil4Smec_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4smec = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densitySmectite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4Smec_setSpec

    subroutine Sil4Calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4calc = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityCalcite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4Calc_setSpec

    subroutine Sil4Quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4quar = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityQuartz)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4Quar_setSpec

    subroutine Sil4Feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4feld = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityFeldspar)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4Feld_setSpec

    subroutine Sil4Hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4hema = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityHematite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4Hema_setSpec

    subroutine Sil4Gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4gyps = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityGypsum)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4Gyps_setSpec

    subroutine Sil4IlHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4ilhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) * densityIllite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4IlHe_setSpec

    subroutine Sil4KaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4kahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) *  densityKaolinite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4KaHe_setSpec

    subroutine Sil4SmHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4smhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densitySmectite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4SmHe_setSpec

    subroutine Sil4CaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4cahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityCalcite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4CaHe_setSpec

    subroutine Sil4QuHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4quhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * DensityQuartz + &
           frIronOxideInAggregate * DensityHematite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4QuHe_setSpec

    subroutine Sil4FeHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4fehe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityFeldspar + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4FeHe_setSpec

    subroutine Sil4GyHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil4gyhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityGypsum + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil4GyHe_setSpec

    subroutine Sil5Illi_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5illi = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityIllite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5Illi_setSpec

    subroutine Sil5Kaol_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5kaol = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityKaolinite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5Kaol_setSpec

    subroutine Sil5Smec_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5smec = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densitySmectite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5Smec_setSpec

    subroutine Sil5Calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5calc = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityCalcite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5Calc_setSpec

    subroutine Sil5Quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5quar = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityQuartz)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5Quar_setSpec

    subroutine Sil5Feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5feld = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityFeldspar)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5Feld_setSpec

    subroutine Sil5Hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5hema = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityHematite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5Hema_setSpec

    subroutine Sil5Gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5gyps = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, densityGypsum)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5Gyps_setSpec

    subroutine Sil5IlHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5ilhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) * densityIllite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5IlHe_setSpec

    subroutine Sil5KaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5kahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1d0 - frIronOxideInAggregate) *  densityKaolinite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5KaHe_setSpec

    subroutine Sil5SmHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5smhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densitySmectite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5SmHe_setSpec

    subroutine Sil5CaHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5cahe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityCalcite + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5CaHe_setSpec

    subroutine Sil5QuHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5quhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * DensityQuartz + &
           frIronOxideInAggregate * DensityHematite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5QuHe_setSpec

    subroutine Sil5FeHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5fehe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityFeldspar + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5FeHe_setSpec

    subroutine Sil5GyHe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_sil5gyhe = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frIronOxideInAggregate) * densityGypsum + &
           frIronOxideInAggregate * densityHematite)
      if (tracers_drydep) call set_trradius(n, 23.53D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil5GyHe_setSpec

  end subroutine Minerals_InitMetadata

end module MineralsTracersMetadata_mod

