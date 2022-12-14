!------------------------------------------------------------------------------
module SeasaltTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  SeasaltTracersMetadata_mod encapsulates the sea salt tracers metadata
!@auth NCCS ASTG, Kostas Tsigaridis
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: nPart, nGAS
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: set_pm2p5fact
  use OldTracer_mod, only: set_pm10fact
  use TRACER_COM, only:  n_seasalt1,  n_seasalt2
  use Tracer_mod, only: Tracer

  implicit none
  private

  public Seasalt_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine Seasalt_InitMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  seasalt1_setSpec('seasalt1')
    call  seasalt2_setSpec('seasalt2')
     
!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine seasalt1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_seasalt1 = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 75.d0)  !Na x 3.256
      call set_trpdens(n, 2.2d3) !kg/m3 This is for non-hydrated
      call set_trradius(n, 4.4d-7 ) ! This is non-hydrated
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
    end subroutine seasalt1_setSpec

    subroutine seasalt2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_seasalt2 = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, 75.d0)  !Na x 3.256
      call set_trpdens(n, 2.2d3) !kg/m3 This is for non-hydrated
      call set_trradius(n, 5.0d-6) ! This is non-hydrated
      call set_trradius(n, 1.7d-6 ) ! This is non-hydrated
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 0.5d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
    end subroutine seasalt2_setSpec

  end subroutine Seasalt_InitMetadata

end module SeasaltTracersMetadata_mod



