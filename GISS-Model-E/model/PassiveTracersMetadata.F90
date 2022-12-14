!------------------------------------------------------------------------------
module PassiveTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  PassiveTracersMetadata_mod encapsulates the Passive tracers metadata
!@auth NCCS ASTG
!  use sharedTracersMetadata_mod
  use TRACER_COM, only: n_SF6, n_SF6_c, n_nh5, n_nh50, n_e90
  use TRACER_COM, only: n_nh15, n_st8025, n_aoa, n_aoanh
  use TRACER_COM, only: n_tape_rec
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_tr_mm, set_ntm_power
  use OldTracer_mod, only: set_tr_mm, set_trdecay
  use Tracer_mod, only: Tracer
  use OldTracer_mod, only: set_has_chemistry
  use CONSTANT, only: mair
  implicit none

  private
  public Passive_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine Passive_InitMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  SF6_setSpec('SF6')
    call  SF6_c_setSpec('SF6_c')
    call  nh5_setSpec('nh5')
    call  nh50_setSpec('nh50')
    call  e90_setSpec('e90')
    call  st8025_setSpec('st8025')
    call  aoa_setSpec('aoa')
    call  aoanh_setSpec('aoanh')
    call  tape_rec_setSpec('tape_rec')
    call  nh15_setSpec('nh15')

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine SF6_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SF6 = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_has_chemistry(n, .true.)
    end subroutine SF6_setSpec

    subroutine SF6_c_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SF6_c = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_has_chemistry(n, .true.)
    end subroutine SF6_c_setSpec


    subroutine nh5_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_nh5 = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, mair)
      call set_trdecay(n, 2.3148d-6) ! 1/(5 days)
      call set_has_chemistry(n, .true.)
    end subroutine nh5_setSpec

    subroutine nh50_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_nh50 = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, mair)
      call set_trdecay(n, 2.3148d-7) ! 1/(50 days)
      call set_has_chemistry(n, .true.)
    end subroutine nh50_setSpec

    subroutine e90_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_e90 = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, mair)
      call set_trdecay(n, 1.2860d-7) ! 1/(90 days)
      call set_has_chemistry(n, .true.)
    end subroutine e90_setSpec

    subroutine st8025_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_st8025 = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, mair)
      call set_has_chemistry(n, .true.)
    end subroutine st8025_setSpec

    subroutine aoa_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_aoa = n
      call set_ntm_power(n, -1)
      call set_tr_mm(n, mair)
      call set_has_chemistry(n, .true.)
    end subroutine aoa_setSpec

    subroutine aoanh_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_aoanh = n
      call set_ntm_power(n, -1)
      call set_tr_mm(n, mair)
      call set_has_chemistry(n, .true.)
    end subroutine aoanh_setSpec

    subroutine tape_rec_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_tape_rec = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, mair)
      call set_has_chemistry(n, .true.)
    end subroutine tape_rec_setSpec

    subroutine nh15_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_nh15 = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, 146.01d0)
      call set_trdecay(n, 7.7160d-7) ! 1/(15 days)
      call set_has_chemistry(n, .true.)
    end subroutine nh15_setSpec

  end subroutine Passive_InitMetadata

end module PassiveTracersMetadata_mod
