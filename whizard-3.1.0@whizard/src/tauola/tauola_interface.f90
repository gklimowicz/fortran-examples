!
!  WHIZARD Tauola interface
!  Adapted from ilc_tauola_mod.f90
!  for Whizard1 developed by Timothy Barklow (SLAC)
!
!    Akiya Miyamoto
!    Bug fixes: Mikael Berggren, Juergen Reuter
!

module tauola_interface

  use kinds
  use iso_varying_string, string_t => varying_string
  use variables
  use model_data

  implicit none

  private

  public :: wo_tauola_pytaud
  public :: tauspin_pyjets
  public :: wo_tauola_get_helicity_mod
  public :: wo_tauola_init_call
  public :: wo_tauola_get_helicity
  public :: taudec_settings_t
  public :: pyjets_spin_t

  !!! THIS COMMON BLOCK IS USED FOR COMMUNICATION WITH TAUOLA
  common /taupos/ np1, np2
  integer ::  np1
  integer ::  np2

  !!! THIS COMMON BLOCK IS USED FOR COMMUNICATION WITH TAUOLA
  COMMON / MOMDEC / Q1,Q2,P1,P2,P3,P4
  double precision Q1(4),Q2(4),P1(4),P2(4),P3(4),P4(4)

  logical, save, public :: trans_spin
  logical, save, public :: tau_pol_vec

  integer, parameter :: n_pyjets_max = 4000
  double precision, dimension(n_pyjets_max) :: tauspin_pyjets

  double precision, external :: wthiggs

  type :: taudec_settings_t
    logical :: photos
    logical :: transverse
    logical :: dec_rad_cor
    integer :: dec_mode1
    integer :: dec_mode2
    real(default) :: mh
    real(default) :: mix_angle
    real(default) :: mtau
    logical :: use_pol_vec
  contains
    procedure :: init => taudec_settings_init
    procedure :: write => taudec_settings_write
  end type taudec_settings_t

  type :: pyjets_spin_t
    integer :: index_to_hepeup  ! =-1, if no matching entry in hepeup
    double precision :: helicity   ! copy of SPINUP
    integer :: pid        ! particle ID
    integer :: id_orig    ! pid of parent
    integer :: index_orig ! index of parent
    integer :: n_daughter ! number of daughter
    integer, dimension(10) :: index_daughter  ! index of daughter particles
  end type pyjets_spin_t

  interface
     module function wo_tauola_get_helicity_mod (ip) result (the_helicity)
       integer, intent(in) :: ip
       integer :: the_helicity
     end function wo_tauola_get_helicity_mod
     module subroutine wo_tauola_get_helicity (ip, the_helicity)
       integer, intent(in)  :: ip
       integer, intent(out) :: the_helicity
     end subroutine wo_tauola_get_helicity
     module subroutine wo_tauola_pytaud (itau, iorig, kforig, ndecay)
       integer, intent(in)  :: itau
       integer, intent(in)  :: iorig
       integer, intent(in)  :: kforig
       integer, intent(out) :: ndecay
     end subroutine wo_tauola_pytaud
     module subroutine wo_tauola_init_call (taudec_settings)
       type(taudec_settings_t), intent(in) :: taudec_settings
     end subroutine wo_tauola_init_call
     module subroutine taudec_settings_init (taudec_settings, var_list, model)
       class(taudec_settings_t), intent(out) :: taudec_settings
       type(var_list_t), intent(in) :: var_list
       class(model_data_t), intent(in) :: model
     end subroutine taudec_settings_init
     module subroutine taudec_settings_write (taudec_settings, unit)
       class(taudec_settings_t), intent(in) :: taudec_settings
       integer, intent(in), optional :: unit
     end subroutine taudec_settings_write
  end interface

  interface
     function pyr(i_dum)
       implicit none
       double precision :: pyr
       integer, intent(in) :: i_dum
     end function pyr
  end interface

end module tauola_interface


!*********************************************************************
!...PYTAUD
!...Routine to handle the decay of a polarized tau lepton.
!...Input:
!...ITAU is the position where the decaying tau is stored in /PYJETS/.
!...IORIG is the position where the mother of the tau is stored;
!...     is 0 when the mother is not stored.
!...KFORIG is the flavour of the mother of the tau;
!...     is 0 when the mother is not known.
!...Note that IORIG=0 does not necessarily imply KFORIG=0;
!...     e.g. in B hadron semileptonic decays the W  propagator
!...     is not explicitly stored but the W code is still unambiguous.
!...Output:
!...NDECAY is the number of decay products in the current tau decay.
!...These decay products should be added to the /PYJETS/ common block,
!...in positions N+1 through N+NDECAY. For each product I you must
!...give the flavour codes K(I,2) and the five-momenta P(I,1), P(I,2),
!...P(I,3), P(I,4) and P(I,5). The rest will be stored automatically.
subroutine pytaud (itau, iorig, kforig, ndecay)
  use tauola_interface !NODEP!
  implicit none
  integer itau,iorig,kforig
  integer ndecay
  !print *,"###############################################"
  !print *,"###### tauola pytaud was called ###############"
  !print *," itau,iorig,kforig=",itau,iorig,kforig
  !print *,"###############################################"

  call wo_tauola_pytaud (itau, iorig, kforig, ndecay)

end subroutine pytaud
