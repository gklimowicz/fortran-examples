! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

module phs_fks

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use constants
  use lorentz
  use phs_points
  use models, only: model_t
  use phs_base
  use resonances, only: resonance_contributors_t, resonance_history_t
  use phs_wood

  implicit none
  private

  public :: isr_kinematics_t
  public :: phs_point_set_t
  public :: real_jacobian_t
  public :: real_kinematics_t
  public :: get_boost_for_threshold_projection
  public :: threshold_projection_born
  public :: compute_dalitz_bounds
  public :: phs_fks_config_t
  public :: dalitz_plot_t
  public :: check_scalar_products
  public :: phs_fks_generator_t
  public :: phs_identifier_t
  public :: check_for_phs_identifier
  public :: phs_fks_t
  public :: get_xi_max_isr
  public :: compute_y_from_emitter
  public :: beta_emitter
  public :: get_filtered_resonance_histories

  integer, parameter, public :: UBF_FSR_SIMPLE = 1
  integer, parameter, public :: UBF_FSR_MASSIVE = 2
  integer, parameter, public :: UBF_FSR_MASSLESS_RECOIL = 3
  integer, parameter, public :: UBF_ISR = 4
  integer, parameter, public :: I_XI = 1
  integer, parameter, public :: I_Y = 2
  integer, parameter, public :: I_PHI = 3

  integer, parameter, public :: PHS_MODE_UNDEFINED = 0
  integer, parameter, public :: PHS_MODE_ADDITIONAL_PARTICLE = 1
  integer, parameter, public :: PHS_MODE_COLLINEAR_REMNANT = 2

  integer, parameter, public :: GEN_REAL_PHASE_SPACE = 1
  integer, parameter, public :: GEN_SOFT_MISMATCH = 2
  integer, parameter, public :: GEN_SOFT_LIMIT_TEST = 3
  integer, parameter, public :: GEN_COLL_LIMIT_TEST = 4
  integer, parameter, public :: GEN_ANTI_COLL_LIMIT_TEST = 5
  integer, parameter, public :: GEN_SOFT_COLL_LIMIT_TEST = 6
  integer, parameter, public :: GEN_SOFT_ANTI_COLL_LIMIT_TEST = 7

  integer, parameter, public :: SQRTS_FIXED = 1
  integer, parameter, public :: SQRTS_VAR = 2

  real(default), parameter :: xi_tilde_test_soft = 0.00001_default
  real(default), parameter :: xi_tilde_test_coll = 0.5_default
  real(default), parameter :: y_test_soft = 0.5_default
  real(default), parameter :: y_test_coll = 0.9999999_default
  !!! for testing EW singularities: y_test_coll = 0.99999999_default

  integer, parameter, public :: I_PLUS = 1
  integer, parameter, public :: I_MINUS = 2


  type :: isr_kinematics_t
    integer :: n_in
    real(default), dimension(2) :: x = one
    real(default), dimension(2) :: z = zero
    real(default) :: sqrts_born = zero
    real(default), dimension(:), allocatable :: beam_energy
    real(default) :: fac_scale = zero
    real(default), dimension(2) :: jacobian = one
    integer :: isr_mode = SQRTS_FIXED
  contains
    procedure :: write => isr_kinematics_write
  end type isr_kinematics_t

  type :: phs_point_set_t
     type(phs_point_t), dimension(:), allocatable :: phs_point
     logical :: initialized = .false.
  contains
    procedure :: init => phs_point_set_init
    procedure :: write => phs_point_set_write
    procedure :: get_n_momenta => phs_point_set_get_n_momenta
    procedure :: get_momenta => phs_point_set_get_momenta
    procedure :: get_momentum => phs_point_set_get_momentum
    procedure :: get_energy => phs_point_set_get_energy
    procedure :: get_sqrts => phs_point_set_get_sqrts
    generic :: set_momenta => set_momenta_p, set_momenta_phs_point
    procedure :: set_momenta_p => phs_point_set_set_momenta_p
    procedure :: set_momenta_phs_point => phs_point_set_set_momenta_phs_point
    procedure :: get_n_particles => phs_point_set_get_n_particles
    procedure :: get_n_phs => phs_point_set_get_n_phs
    procedure :: get_invariant_mass => phs_point_set_get_invariant_mass
    procedure :: write_phs_point => phs_point_set_write_phs_point
    procedure :: final => phs_point_set_final
  end type phs_point_set_t

  type :: real_jacobian_t
    real(default), dimension(4) :: jac = 1._default
  end type real_jacobian_t

  type :: real_kinematics_t
    logical :: supply_xi_max = .true.
    real(default) :: xi_tilde
    real(default) :: phi
    real(default), dimension(:), allocatable :: xi_max, y
    real(default) :: xi_mismatch, y_mismatch
    type(real_jacobian_t), dimension(:), allocatable :: jac
    real(default) :: jac_mismatch
    type(phs_point_set_t) :: p_born_cms
    type(phs_point_set_t) :: p_born_lab
    type(phs_point_set_t) :: p_real_cms
    type(phs_point_set_t) :: p_real_lab
    type(phs_point_set_t) :: p_born_onshell
    type(phs_point_set_t), dimension(2) :: p_real_onshell
    integer, dimension(:), allocatable :: alr_to_i_phs
    real(default), dimension(3) :: x_rad
    real(default), dimension(:), allocatable :: jac_rand
    real(default), dimension(:), allocatable :: y_soft
    real(default) :: cms_energy2
    type(vector4_t), dimension(:), allocatable :: xi_ref_momenta
  contains
    procedure :: init => real_kinematics_init
    procedure :: init_onshell => real_kinematics_init_onshell
    procedure :: write => real_kinematics_write
    procedure :: apply_threshold_projection_real => &
         real_kinematics_apply_threshold_projection_real
    procedure :: kt2 => real_kinematics_kt2
    procedure :: final => real_kinematics_final
  end type real_kinematics_t

  type, extends (phs_wood_config_t) :: phs_fks_config_t
    integer :: mode = PHS_MODE_UNDEFINED
    character(32) :: md5sum_born_config
    logical :: born_2_to_1 = .false.
    logical :: make_dalitz_plot = .false.
  contains
    procedure :: clear_phase_space => fks_config_clear_phase_space
    procedure :: write => phs_fks_config_write
    procedure :: set_mode => phs_fks_config_set_mode
    procedure :: configure => phs_fks_config_configure
    procedure :: startup_message => phs_fks_config_startup_message
    procedure, nopass :: allocate_instance => phs_fks_config_allocate_instance
    procedure :: generate_phase_space_extra => &
         phs_fks_config_generate_phase_space_extra
    procedure :: set_born_config => phs_fks_config_set_born_config
    procedure :: get_resonance_histories => &
         phs_fks_config_get_resonance_histories
  end type phs_fks_config_t

  type :: dalitz_plot_t
     integer :: unit = -1
     type(string_t) :: filename
     logical :: active = .false.
     logical :: inverse = .false.
  contains
    procedure :: init => dalitz_plot_init
    procedure :: write_header => dalitz_plot_write_header
    procedure :: register => dalitz_plot_register
    procedure :: final => dalitz_plot_final
  end type dalitz_plot_t

  type :: phs_fks_generator_t
    integer, dimension(:), allocatable :: emitters
    type(real_kinematics_t), pointer :: real_kinematics => null()
    type(isr_kinematics_t), pointer :: isr_kinematics => null()
    integer :: n_in
    real(default) :: xi_min
    real(default) :: y_max
    real(default) :: sqrts
    real(default) :: E_gluon
    real(default) :: mrec2
    real(default), dimension(:), allocatable :: m2
    logical :: massive_phsp = .false.
    logical, dimension(:), allocatable :: is_massive
    logical :: singular_jacobian = .false.
    integer :: i_fsr_first = -1
    type(resonance_contributors_t), dimension(:), allocatable :: resonance_contributors !!! Put somewhere else?
    integer :: mode = GEN_REAL_PHASE_SPACE
  contains
    procedure :: connect_kinematics => phs_fks_generator_connect_kinematics
    procedure :: compute_isr_kinematics => &
         phs_fks_generator_compute_isr_kinematics
    procedure :: final => phs_fks_generator_final
    generic :: generate_fsr => generate_fsr_default, generate_fsr_resonances
    procedure :: generate_fsr_default => phs_fks_generator_generate_fsr_default
    procedure :: generate_fsr_resonances => &
         phs_fks_generator_generate_fsr_resonances
    procedure :: generate_fsr_threshold => &
         phs_fks_generator_generate_fsr_threshold
    procedure :: generate_fsr_in => phs_fks_generator_generate_fsr_in
    procedure :: generate_fsr_out => phs_fks_generator_generate_fsr_out
    generic :: compute_emitter_kinematics => &
       compute_emitter_kinematics_massless, &
       compute_emitter_kinematics_massive
    procedure :: compute_emitter_kinematics_massless => &
       phs_fks_generator_compute_emitter_kinematics_massless
    procedure :: compute_emitter_kinematics_massive => &
       phs_fks_generator_compute_emitter_kinematics_massive
    procedure :: generate_isr_fixed_beam_energy => &
         phs_fks_generator_generate_isr_fixed_beam_energy
    procedure :: generate_isr_factorized => &
         phs_fks_generator_generate_isr_factorized
    procedure :: generate_isr => phs_fks_generator_generate_isr
    procedure :: set_sqrts_hat => phs_fks_generator_set_sqrts_hat
    procedure :: set_emitters => phs_fks_generator_set_emitters
    procedure :: setup_masses => phs_fks_generator_setup_masses
    procedure :: set_xi_and_y_bounds => phs_fks_generator_set_xi_and_y_bounds
    procedure :: set_isr_kinematics => phs_fks_generator_set_isr_kinematics
    procedure :: generate_radiation_variables => &
       phs_fks_generator_generate_radiation_variables
    procedure :: compute_xi_ref_momenta => &
         phs_fks_generator_compute_xi_ref_momenta
    procedure :: compute_xi_ref_momenta_threshold &
         => phs_fks_generator_compute_xi_ref_momenta_threshold
    procedure :: compute_cms_energy => phs_fks_generator_compute_cms_energy
    procedure :: compute_xi_max => phs_fks_generator_compute_xi_max
    procedure :: compute_xi_max_isr_factorized &
       => phs_fks_generator_compute_xi_max_isr_factorized
    procedure :: set_masses => phs_fks_generator_set_masses
    procedure :: compute_y_real_phs => phs_fks_generator_compute_y_real_phs
    procedure :: compute_y_mismatch => phs_fks_generator_compute_y_mismatch
    procedure :: compute_y_test => phs_fks_generator_compute_y_test
    procedure :: compute_xi_tilde => phs_fks_generator_compute_xi_tilde
    procedure :: prepare_generation => phs_fks_generator_prepare_generation
    procedure :: generate_fsr_from_xi_and_y => &
       phs_fks_generator_generate_fsr_from_xi_and_y
    procedure :: generate_isr_from_xi_and_y => &
       phs_fks_generator_generate_isr_from_xi_and_y
    procedure :: get_radiation_variables => &
       phs_fks_generator_get_radiation_variables
    procedure :: write => phs_fks_generator_write
  end type phs_fks_generator_t

  type :: phs_identifier_t
     integer, dimension(:), allocatable :: contributors
     integer :: emitter = -1
     logical :: evaluated = .false.
  contains
    generic :: init => init_from_emitter, init_from_emitter_and_contributors
    procedure :: init_from_emitter => phs_identifier_init_from_emitter
    procedure :: init_from_emitter_and_contributors &
       => phs_identifier_init_from_emitter_and_contributors
    procedure :: check => phs_identifier_check
    procedure :: write => phs_identifier_write
  end type phs_identifier_t

  type, extends (phs_wood_t) :: phs_fks_t
    integer :: mode = PHS_MODE_UNDEFINED
    type(vector4_t), dimension(:), allocatable :: p_born
    type(vector4_t), dimension(:), allocatable :: q_born
    type(vector4_t), dimension(:), allocatable :: p_real
    type(vector4_t), dimension(:), allocatable :: q_real
    type(vector4_t), dimension(:), allocatable :: p_born_tot
    type(phs_fks_generator_t) :: generator
    real(default) :: r_isr
    type(phs_identifier_t), dimension(:), allocatable :: phs_identifiers
  contains
    procedure :: write => phs_fks_write
    procedure :: init => phs_fks_init
    procedure :: compute_flux => phs_fks_compute_flux
    procedure :: allocate_momenta => phs_fks_allocate_momenta
    procedure :: evaluate_selected_channel => phs_fks_evaluate_selected_channel
    procedure :: evaluate_other_channels => phs_fks_evaluate_other_channels
    procedure :: get_mcpar => phs_fks_get_mcpar
    procedure :: set_beam_energy => phs_fks_set_beam_energy
    procedure :: set_emitters => phs_fks_set_emitters
    procedure :: set_momenta => phs_fks_set_momenta
    procedure :: setup_masses => phs_fks_setup_masses
    procedure :: get_born_momenta => phs_fks_get_born_momenta
    procedure :: get_outgoing_momenta => phs_fks_get_outgoing_momenta
    procedure :: get_incoming_momenta => phs_fks_get_incoming_momenta
    procedure :: set_isr_kinematics => phs_fks_set_isr_kinematics
    procedure :: generate_radiation_variables => &
       phs_fks_generate_radiation_variables
    procedure :: compute_xi_ref_momenta => phs_fks_compute_xi_ref_momenta
    procedure :: compute_xi_ref_momenta_threshold => &
         phs_fks_compute_xi_ref_momenta_threshold
    procedure :: compute_cms_energy => phs_fks_compute_cms_energy
    procedure :: set_reference_frames => phs_fks_set_reference_frames
    procedure :: i_phs_is_isr => phs_fks_i_phs_is_isr
    procedure :: generate_fsr_in => phs_fks_generate_fsr_in
    procedure :: generate_fsr => phs_fks_generate_fsr
    procedure :: get_onshell_projected_momenta => &
         phs_fks_get_onshell_projected_momenta
    procedure :: generate_fsr_threshold => phs_fks_generate_fsr_threshold
    generic :: compute_xi_max => &
         compute_xi_max_internal, compute_xi_max_with_output
    procedure :: compute_xi_max_internal => phs_fks_compute_xi_max_internal
    procedure :: compute_xi_max_with_output => phs_fks_compute_xi_max_with_output
    procedure :: generate_isr => phs_fks_generate_isr
    procedure :: compute_isr_kinematics => phs_fks_compute_isr_kinematics
    procedure :: final => phs_fks_final
  end type phs_fks_t


  interface
    module subroutine isr_kinematics_write (isr, unit)
      class(isr_kinematics_t), intent(in) :: isr
      integer, intent(in), optional :: unit
    end subroutine isr_kinematics_write
    module subroutine phs_point_set_init (phs_point_set, n_particles, n_phs)
      class(phs_point_set_t), intent(out) :: phs_point_set
      integer, intent(in) :: n_particles, n_phs
    end subroutine phs_point_set_init
    module subroutine phs_point_set_write (phs_point_set, i_phs, contributors, &
         unit, show_mass, testflag, check_conservation, ultra, n_in)
      class(phs_point_set_t), intent(in) :: phs_point_set
      integer, intent(in), optional :: i_phs
      integer, intent(in), dimension(:), optional :: contributors
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_mass
      logical, intent(in), optional :: testflag, ultra
      logical, intent(in), optional :: check_conservation
      integer, intent(in), optional :: n_in
    end subroutine phs_point_set_write
    elemental module function phs_point_set_get_n_momenta &
         (phs_point_set, i_res) result (n)
      integer :: n
      class(phs_point_set_t), intent(in) :: phs_point_set
      integer, intent(in) :: i_res
    end function phs_point_set_get_n_momenta
    pure module function phs_point_set_get_momenta &
         (phs_point_set, i_phs, n_in) result (p)
      type(vector4_t), dimension(:), allocatable :: p
      class(phs_point_set_t), intent(in) :: phs_point_set
      integer, intent(in) :: i_phs
      integer, intent(in), optional :: n_in
    end function phs_point_set_get_momenta
    pure module function phs_point_set_get_momentum &
         (phs_point_set, i_phs, i_mom) result (p)
      type(vector4_t) :: p
      class(phs_point_set_t), intent(in) :: phs_point_set
      integer, intent(in) :: i_phs, i_mom
    end function phs_point_set_get_momentum
    pure module function phs_point_set_get_energy &
         (phs_point_set, i_phs, i_mom) result (E)
      real(default) :: E
      class(phs_point_set_t), intent(in) :: phs_point_set
      integer, intent(in) :: i_phs, i_mom
    end function phs_point_set_get_energy
    module function phs_point_set_get_sqrts &
         (phs_point_set, i_phs) result (sqrts)
      real(default) :: sqrts
      class(phs_point_set_t), intent(in) :: phs_point_set
      integer, intent(in) :: i_phs
    end function phs_point_set_get_sqrts
    module subroutine phs_point_set_set_momenta_p (phs_point_set, i_phs, p)
      class(phs_point_set_t), intent(inout) :: phs_point_set
      integer, intent(in) :: i_phs
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine phs_point_set_set_momenta_p
    module subroutine phs_point_set_set_momenta_phs_point &
         (phs_point_set, i_phs, p)
      class(phs_point_set_t), intent(inout) :: phs_point_set
      integer, intent(in) :: i_phs
      type(phs_point_t), intent(in) :: p
    end subroutine phs_point_set_set_momenta_phs_point
    module function phs_point_set_get_n_particles &
         (phs_point_set, i) result (n_particles)
      integer :: n_particles
      class(phs_point_set_t), intent(in) :: phs_point_set
      integer, intent(in), optional :: i
    end function phs_point_set_get_n_particles
    module function phs_point_set_get_n_phs (phs_point_set) result (n_phs)
      integer :: n_phs
      class(phs_point_set_t), intent(in) :: phs_point_set
    end function phs_point_set_get_n_phs
    module function phs_point_set_get_invariant_mass &
         (phs_point_set, i_phs, i_part) result (m2)
      real(default) :: m2
      class(phs_point_set_t), intent(in) :: phs_point_set
      integer, intent(in) :: i_phs
      integer, intent(in), dimension(:) :: i_part
    end function phs_point_set_get_invariant_mass
    module subroutine phs_point_set_write_phs_point (phs_point_set, i_phs, &
         unit, show_mass, testflag, check_conservation, ultra, n_in)
      class(phs_point_set_t), intent(in) :: phs_point_set
      integer, intent(in) :: i_phs
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_mass
      logical, intent(in), optional :: testflag, ultra
      logical, intent(in), optional :: check_conservation
      integer, intent(in), optional :: n_in
    end subroutine phs_point_set_write_phs_point
    module subroutine phs_point_set_final (phs_point_set)
      class(phs_point_set_t), intent(inout) :: phs_point_set
    end subroutine phs_point_set_final
    module subroutine real_kinematics_init (r, n_tot, n_phs, n_alr, n_contr)
      class(real_kinematics_t), intent(inout) :: r
      integer, intent(in) :: n_tot, n_phs, n_alr, n_contr
    end subroutine real_kinematics_init
    module subroutine real_kinematics_init_onshell (r, n_tot, n_phs)
      class(real_kinematics_t), intent(inout) :: r
      integer, intent(in) :: n_tot, n_phs
    end subroutine real_kinematics_init_onshell
    module subroutine real_kinematics_write (r, unit)
      class(real_kinematics_t), intent(in) :: r
      integer, intent(in), optional :: unit
    end subroutine real_kinematics_write
    module function get_boost_for_threshold_projection &
         (p, sqrts, mtop) result (L)
      type(lorentz_transformation_t) :: L
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in) :: sqrts, mtop
    end function get_boost_for_threshold_projection
    module subroutine real_kinematics_apply_threshold_projection_real &
         (r, i_phs, mtop, L_to_cms, invert)
      class(real_kinematics_t), intent(inout) :: r
      integer, intent(in) :: i_phs
      real(default), intent(in) :: mtop
      type(lorentz_transformation_t), intent(in), dimension(:) :: L_to_cms
      logical, intent(in) :: invert
    end subroutine real_kinematics_apply_threshold_projection_real
    module subroutine threshold_projection_born &
         (mtop, L_to_cms, p_in, p_onshell)
      real(default), intent(in) :: mtop
      type(lorentz_transformation_t), intent(in) :: L_to_cms
      type(vector4_t), intent(in), dimension(:) :: p_in
      type(vector4_t), intent(out), dimension(:) :: p_onshell
    end subroutine threshold_projection_born
    pure module subroutine compute_dalitz_bounds &
         (q0, m2, mrec2, z1, z2, k0_rec_max)
      real(default), intent(in) :: q0, m2, mrec2
      real(default), intent(out) :: z1, z2, k0_rec_max
    end subroutine compute_dalitz_bounds
    module function real_kinematics_kt2 &
         (real_kinematics, i_phs, emitter, kt2_type, xi, y) result (kt2)
      real(default) :: kt2
      class(real_kinematics_t), intent(in) :: real_kinematics
      integer, intent(in) :: emitter, i_phs, kt2_type
      real(default), intent(in), optional :: xi, y
    end function real_kinematics_kt2
    module subroutine real_kinematics_final (real_kin)
      class(real_kinematics_t), intent(inout) :: real_kin
    end subroutine real_kinematics_final
    module subroutine fks_config_clear_phase_space (phs_config)
      class(phs_fks_config_t), intent(inout) :: phs_config
    end subroutine fks_config_clear_phase_space
    module subroutine phs_fks_config_write (object, unit, include_id)
      class(phs_fks_config_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: include_id
    end subroutine phs_fks_config_write
    module subroutine phs_fks_config_set_mode (phs_config, mode)
      class(phs_fks_config_t), intent(inout) :: phs_config
      integer, intent(in) :: mode
    end subroutine phs_fks_config_set_mode
    module subroutine phs_fks_config_configure (phs_config, sqrts, &
         sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
         ignore_mismatch, nlo_type, subdir)
      class(phs_fks_config_t), intent(inout) :: phs_config
      real(default), intent(in) :: sqrts
      logical, intent(in), optional :: sqrts_fixed
      logical, intent(in), optional :: lab_is_cm
      logical, intent(in), optional :: azimuthal_dependence
      logical, intent(in), optional :: rebuild
      logical, intent(in), optional :: ignore_mismatch
      integer, intent(in), optional :: nlo_type
      type(string_t), intent(in), optional :: subdir
    end subroutine phs_fks_config_configure
    module subroutine phs_fks_config_startup_message (phs_config, unit)
      class(phs_fks_config_t), intent(in) :: phs_config
      integer, intent(in), optional :: unit
    end subroutine phs_fks_config_startup_message
    module subroutine phs_fks_config_generate_phase_space_extra (phs_config)
      class(phs_fks_config_t), intent(inout) :: phs_config
    end subroutine phs_fks_config_generate_phase_space_extra
    module subroutine phs_fks_config_set_born_config (phs_config, phs_cfg_born)
      class(phs_fks_config_t), intent(inout) :: phs_config
      type(phs_wood_config_t), intent(in), target :: phs_cfg_born
    end subroutine phs_fks_config_set_born_config
    module function phs_fks_config_get_resonance_histories &
         (phs_config) result (resonance_histories)
      type(resonance_history_t), dimension(:), allocatable :: &
           resonance_histories
      class(phs_fks_config_t), intent(inout) :: phs_config
    end function phs_fks_config_get_resonance_histories
    module subroutine dalitz_plot_init (plot, unit, filename, inverse)
      class(dalitz_plot_t), intent(inout) :: plot
      integer, intent(in) :: unit
      type(string_t), intent(in) :: filename
      logical, intent(in) :: inverse
    end subroutine dalitz_plot_init
    module subroutine dalitz_plot_write_header (plot)
      class(dalitz_plot_t), intent(in) :: plot
    end subroutine dalitz_plot_write_header
    module subroutine dalitz_plot_register (plot, k0_n, k0_np1)
      class(dalitz_plot_t), intent(in) :: plot
      real(default), intent(in) :: k0_n, k0_np1
    end subroutine dalitz_plot_register
    module subroutine dalitz_plot_final (plot)
      class(dalitz_plot_t), intent(inout) :: plot
    end subroutine dalitz_plot_final
    module function check_scalar_products (p) result (valid)
      logical :: valid
      type(vector4_t), intent(in), dimension(:) :: p
    end function check_scalar_products
    module subroutine phs_fks_generator_connect_kinematics &
         (generator, isr_kinematics, real_kinematics, massive_phsp)
      class(phs_fks_generator_t), intent(inout) :: generator
      type(isr_kinematics_t), intent(in), pointer :: isr_kinematics
      type(real_kinematics_t), intent(in), pointer :: real_kinematics
      logical, intent(in) :: massive_phsp
    end subroutine phs_fks_generator_connect_kinematics
    module subroutine phs_fks_generator_compute_isr_kinematics &
         (generator, r, p_in)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in) :: r
      type(vector4_t), dimension(2), intent(in), optional :: p_in
    end subroutine phs_fks_generator_compute_isr_kinematics
    module subroutine phs_fks_generator_final (generator)
      class(phs_fks_generator_t), intent(inout) :: generator
    end subroutine phs_fks_generator_final
    module subroutine phs_identifier_init_from_emitter (phs_id, emitter)
      class(phs_identifier_t), intent(out) :: phs_id
      integer, intent(in) :: emitter
    end subroutine phs_identifier_init_from_emitter
    module subroutine phs_identifier_init_from_emitter_and_contributors &
       (phs_id, emitter, contributors)
       class(phs_identifier_t), intent(out) :: phs_id
       integer, intent(in) :: emitter
       integer, intent(in), dimension(:) :: contributors
    end subroutine phs_identifier_init_from_emitter_and_contributors
    module function phs_identifier_check &
         (phs_id, emitter, contributors) result (check)
      logical :: check
      class(phs_identifier_t), intent(in) :: phs_id
      integer, intent(in) :: emitter
      integer, intent(in), dimension(:), optional :: contributors
    end function phs_identifier_check
    module subroutine phs_identifier_write (phs_id, unit)
      class(phs_identifier_t), intent(in) :: phs_id
      integer, intent(in), optional :: unit
    end subroutine phs_identifier_write
    module subroutine check_for_phs_identifier &
         (phs_id, n_in, emitter, contributors, phs_exist, i_phs)
      type(phs_identifier_t), intent(in), dimension(:) :: phs_id
      integer, intent(in) :: n_in, emitter
      integer, intent(in), dimension(:), optional :: contributors
      logical, intent(out) :: phs_exist
      integer, intent(out) :: i_phs
    end subroutine check_for_phs_identifier
    module subroutine phs_fks_write (object, unit, verbose)
      class(phs_fks_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine phs_fks_write
    module subroutine phs_fks_init (phs, phs_config)
      class(phs_fks_t), intent(out) :: phs
      class(phs_config_t), intent(in), target :: phs_config
    end subroutine phs_fks_init
    module subroutine phs_fks_compute_flux (phs)
      class(phs_fks_t), intent(inout) :: phs
    end subroutine phs_fks_compute_flux
    module subroutine phs_fks_allocate_momenta (phs, phs_config, data_is_born)
      class(phs_fks_t), intent(inout) :: phs
      class(phs_config_t), intent(in) :: phs_config
      logical, intent(in) :: data_is_born
    end subroutine phs_fks_allocate_momenta
    module subroutine phs_fks_evaluate_selected_channel (phs, c_in, r_in)
      class(phs_fks_t), intent(inout) :: phs
      integer, intent(in) :: c_in
      real(default), intent(in), dimension(:) :: r_in
    end subroutine phs_fks_evaluate_selected_channel
    module subroutine phs_fks_evaluate_other_channels (phs, c_in)
      class(phs_fks_t), intent(inout) :: phs
      integer, intent(in) :: c_in
    end subroutine phs_fks_evaluate_other_channels
    module subroutine phs_fks_get_mcpar (phs, c, r)
      class(phs_fks_t), intent(in) :: phs
      integer, intent(in) :: c
      real(default), dimension(:), intent(out) :: r
    end subroutine phs_fks_get_mcpar
    module subroutine phs_fks_set_beam_energy (phs)
      class(phs_fks_t), intent(inout) :: phs
    end subroutine phs_fks_set_beam_energy
    module subroutine phs_fks_set_emitters (phs, emitters)
      class(phs_fks_t), intent(inout) :: phs
      integer, intent(in), dimension(:), allocatable :: emitters
    end subroutine phs_fks_set_emitters
    module subroutine phs_fks_set_momenta (phs, p)
      class(phs_fks_t), intent(inout) :: phs
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine phs_fks_set_momenta
    module subroutine phs_fks_setup_masses (phs, n_tot)
      class(phs_fks_t), intent(inout) :: phs
      integer, intent(in) :: n_tot
    end subroutine phs_fks_setup_masses
    module subroutine phs_fks_get_born_momenta (phs, p)
      class(phs_fks_t), intent(inout) :: phs
      type(vector4_t), intent(out), dimension(:) :: p
    end subroutine phs_fks_get_born_momenta
    module subroutine phs_fks_get_outgoing_momenta (phs, q)
      class(phs_fks_t), intent(in) :: phs
      type(vector4_t), intent(out), dimension(:) :: q
    end subroutine phs_fks_get_outgoing_momenta
    module subroutine phs_fks_get_incoming_momenta (phs, p)
      class(phs_fks_t), intent(in) :: phs
      type(vector4_t), intent(inout), dimension(:), allocatable :: p
    end subroutine phs_fks_get_incoming_momenta
    module subroutine phs_fks_set_isr_kinematics (phs, requires_boost)
      class(phs_fks_t), intent(inout) :: phs
      logical, intent(in) :: requires_boost
    end subroutine phs_fks_set_isr_kinematics
    module subroutine phs_fks_generate_radiation_variables &
         (phs, r_in, threshold)
      class(phs_fks_t), intent(inout) :: phs
      real(default), intent(in), dimension(:) :: r_in
      logical, intent(in) :: threshold
    end subroutine phs_fks_generate_radiation_variables
    module subroutine phs_fks_compute_xi_ref_momenta (phs, p_in, contributors)
      class(phs_fks_t), intent(inout) :: phs
      type(vector4_t), intent(in), dimension(:), optional :: p_in
      type(resonance_contributors_t), intent(in), dimension(:), optional :: &
           contributors
    end subroutine phs_fks_compute_xi_ref_momenta
    module subroutine phs_fks_compute_xi_ref_momenta_threshold (phs)
      class(phs_fks_t), intent(inout) :: phs
    end subroutine phs_fks_compute_xi_ref_momenta_threshold
    module subroutine phs_fks_compute_cms_energy (phs)
      class(phs_fks_t), intent(inout) :: phs
    end subroutine phs_fks_compute_cms_energy
    module subroutine phs_fks_set_reference_frames (phs, is_cms)
      class(phs_fks_t), intent(inout) :: phs
      logical, intent(in) :: is_cms
    end subroutine phs_fks_set_reference_frames
    module function phs_fks_i_phs_is_isr (phs, i_phs) result (is_isr)
      logical :: is_isr
      class(phs_fks_t), intent(in) :: phs
      integer, intent(in) :: i_phs
    end function phs_fks_i_phs_is_isr
    module subroutine phs_fks_generator_generate_fsr_default &
         (generator, emitter, i_phs, &
         p_born, p_real, xi_y_phi, no_jacobians)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: emitter, i_phs
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(out), dimension(:) :: p_real
      real(default), intent(in), dimension(3), optional :: xi_y_phi
      logical, intent(in), optional :: no_jacobians
    end subroutine phs_fks_generator_generate_fsr_default
    module subroutine phs_fks_generator_generate_fsr_resonances (generator, &
         emitter, i_phs, i_con, p_born, p_real, xi_y_phi, no_jacobians)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: emitter, i_phs
      integer, intent(in) :: i_con
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(inout), dimension(:) :: p_real
      real(default), intent(in), dimension(3), optional :: xi_y_phi
      logical, intent(in), optional :: no_jacobians
    end subroutine phs_fks_generator_generate_fsr_resonances
    module subroutine phs_fks_generator_generate_fsr_threshold (generator, &
         emitter, i_phs, p_born, p_real, xi_y_phi)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: emitter, i_phs
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(inout), dimension(:) :: p_real
      real(default), intent(in), dimension(3), optional :: xi_y_phi
    end subroutine phs_fks_generator_generate_fsr_threshold
    module subroutine phs_fks_generator_generate_fsr_in &
         (generator, p_born, p_real)
      class(phs_fks_generator_t), intent(inout) :: generator
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(out), dimension(:) :: p_real
    end subroutine phs_fks_generator_generate_fsr_in
    module subroutine phs_fks_generator_generate_fsr_out (generator, &
        emitter, i_phs, p_born, p_real, q0, p_emitter_index, &
        xi_y_phi, no_jacobians)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: emitter, i_phs
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(inout), dimension(:) :: p_real
      real(default), intent(in) :: q0
      integer, intent(in), optional :: p_emitter_index
      real(default), intent(in), dimension(3), optional :: xi_y_phi
      logical, intent(in), optional :: no_jacobians
    end subroutine phs_fks_generator_generate_fsr_out
    module subroutine phs_fks_generate_fsr_in (phs)
      class(phs_fks_t), intent(inout) :: phs
    end subroutine phs_fks_generate_fsr_in
    module subroutine phs_fks_generate_fsr (phs, emitter, i_phs, p_real, &
         i_con, xi_y_phi, no_jacobians)
      class(phs_fks_t), intent(inout) :: phs
      integer, intent(in) :: emitter, i_phs
      type(vector4_t), intent(out), dimension(:) :: p_real
      integer, intent(in), optional :: i_con
      real(default), intent(in), dimension(3), optional :: xi_y_phi
      logical, intent(in), optional :: no_jacobians
    end subroutine phs_fks_generate_fsr
    pure module function phs_fks_get_onshell_projected_momenta (phs) result (p)
      type(vector4_t), dimension(:), allocatable :: p
      class(phs_fks_t), intent(in) :: phs
    end function phs_fks_get_onshell_projected_momenta
    module subroutine phs_fks_generate_fsr_threshold &
         (phs, emitter, i_phs, p_real)
      class(phs_fks_t), intent(inout) :: phs
      integer, intent(in) :: emitter, i_phs
      type(vector4_t), intent(inout), dimension(:), optional :: p_real
    end subroutine phs_fks_generate_fsr_threshold
    module subroutine phs_fks_compute_xi_max_internal (phs, p, threshold)
      class(phs_fks_t), intent(inout) :: phs
      type(vector4_t), intent(in), dimension(:) :: p
      logical, intent(in) :: threshold
    end subroutine phs_fks_compute_xi_max_internal
    module subroutine phs_fks_compute_xi_max_with_output &
         (phs, emitter, i_phs, y, p, xi_max)
      class(phs_fks_t), intent(inout) :: phs
      integer, intent(in) :: i_phs, emitter
      real(default), intent(in) :: y
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(out) :: xi_max
    end subroutine phs_fks_compute_xi_max_with_output
    module subroutine phs_fks_generator_compute_emitter_kinematics_massless &
         (generator, y, q0, uk_em, uk)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in) :: y, q0
      real(default), intent(out) :: uk_em, uk
    end subroutine phs_fks_generator_compute_emitter_kinematics_massless
    module subroutine phs_fks_generator_compute_emitter_kinematics_massive &
         (generator, y, em, i_phs, q0, k0_em, uk_em, uk, compute_jac)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in) :: y
      integer, intent(in) :: em, i_phs
      real(default), intent(in) :: q0
      real(default), intent(inout) :: k0_em, uk_em, uk
      logical, intent(in) :: compute_jac
    end subroutine phs_fks_generator_compute_emitter_kinematics_massive
    module function get_xi_max_isr (xb, y) result (xi_max)
      real(default) :: xi_max
      real(default), dimension(2), intent(in) :: xb
      real(default), intent(in) :: y
    end function get_xi_max_isr
    module subroutine phs_fks_generate_isr (phs, i_phs, p_real)
      class(phs_fks_t), intent(inout) :: phs
      integer, intent(in) :: i_phs
      type(vector4_t), intent(out), dimension(:) :: p_real
    end subroutine phs_fks_generate_isr
    module subroutine phs_fks_generator_generate_isr_fixed_beam_energy &
         (generator, i_phs, p_born, p_real)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: i_phs
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(out), dimension(:) :: p_real
    end subroutine phs_fks_generator_generate_isr_fixed_beam_energy
    module subroutine phs_fks_generator_generate_isr_factorized &
         (generator, i_phs, emitter, p_born, p_real)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: i_phs, emitter
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(out), dimension(:) :: p_real
    end subroutine phs_fks_generator_generate_isr_factorized
    module subroutine phs_fks_generator_generate_isr &
         (generator, i_phs, p_born, p_real)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: i_phs
      type(vector4_t), intent(in) , dimension(:) :: p_born
      type(vector4_t), intent(out), dimension(:) :: p_real
    end subroutine phs_fks_generator_generate_isr
    module subroutine phs_fks_generator_set_sqrts_hat (generator, sqrts)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in) :: sqrts
    end subroutine phs_fks_generator_set_sqrts_hat
    module subroutine phs_fks_generator_set_emitters (generator, emitters)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in), dimension(:), allocatable ::  emitters
    end subroutine phs_fks_generator_set_emitters
    module subroutine phs_fks_generator_setup_masses (generator, n_tot)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: n_tot
    end subroutine phs_fks_generator_setup_masses
    module subroutine phs_fks_generator_set_xi_and_y_bounds &
         (generator, fks_xi_min, fks_y_max)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in), optional :: fks_xi_min, fks_y_max
    end subroutine phs_fks_generator_set_xi_and_y_bounds
    module subroutine phs_fks_generator_set_isr_kinematics (generator, p)
      class(phs_fks_generator_t), intent(inout) :: generator
      type(vector4_t), dimension(2), intent(in) :: p
    end subroutine phs_fks_generator_set_isr_kinematics
    module subroutine phs_fks_generator_generate_radiation_variables &
         (generator, r_in, p_born, phs_identifiers, threshold)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in), dimension(:) :: r_in
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
      logical, intent(in), optional :: threshold
    end subroutine phs_fks_generator_generate_radiation_variables
    module subroutine phs_fks_generator_compute_xi_ref_momenta &
         (generator, p_born, resonance_contributors)
      class(phs_fks_generator_t), intent(inout) :: generator
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(resonance_contributors_t), intent(in), dimension(:), optional &
         :: resonance_contributors
    end subroutine phs_fks_generator_compute_xi_ref_momenta
    module subroutine phs_fks_generator_compute_xi_ref_momenta_threshold &
         (generator, p_born)
      class(phs_fks_generator_t), intent(inout) :: generator
      type(vector4_t), intent(in), dimension(:) :: p_born
    end subroutine phs_fks_generator_compute_xi_ref_momenta_threshold
    module subroutine phs_fks_generator_compute_cms_energy (generator, p_born)
      class(phs_fks_generator_t), intent(inout) :: generator
      type(vector4_t), intent(in), dimension(:) :: p_born
    end subroutine phs_fks_generator_compute_cms_energy
    module subroutine phs_fks_generator_compute_xi_max (generator, emitter, &
         i_phs, p, xi_max, i_con, y_in)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: i_phs, emitter
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(out) :: xi_max
      integer, intent(in), optional :: i_con
      real(default), intent(in), optional :: y_in
    end subroutine phs_fks_generator_compute_xi_max
    module subroutine phs_fks_generator_compute_xi_max_isr_factorized &
       (generator, i_phs, p)
      class(phs_fks_generator_t), intent(inout) :: generator
      integer, intent(in) :: i_phs
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine phs_fks_generator_compute_xi_max_isr_factorized
    module subroutine phs_fks_generator_set_masses &
         (generator, p, phs_identifiers)
      class(phs_fks_generator_t), intent(inout) :: generator
      type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine phs_fks_generator_set_masses
    module subroutine compute_y_from_emitter (r_y, p, n_in, emitter, &
         massive, y_max, jac_rand, y, contributors, threshold)
      real(default), intent(in) :: r_y
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: n_in
      integer, intent(in) :: emitter
      logical, intent(in) :: massive
      real(default), intent(in) :: y_max
      real(default), intent(inout) :: jac_rand
      real(default), intent(out) :: y
      integer, intent(in), dimension(:), allocatable, optional :: contributors
      logical, intent(in), optional :: threshold
    end subroutine compute_y_from_emitter
    module subroutine phs_fks_generator_compute_y_real_phs &
         (generator, r_y, p, phs_identifiers, &
         jac_rand, y, threshold)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in) :: r_y
      type(vector4_t), intent(in), dimension(:) :: p
      type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
      real(default), intent(inout), dimension(:) :: jac_rand
      real(default), intent(out), dimension(:) :: y
      logical, intent(in), optional :: threshold
    end subroutine phs_fks_generator_compute_y_real_phs
    module subroutine phs_fks_generator_compute_y_mismatch &
         (generator, r_y, jac_rand, y, y_soft)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in) :: r_y
      real(default), intent(inout) :: jac_rand
      real(default), intent(out) :: y
      real(default), intent(out), dimension(:) :: y_soft
    end subroutine phs_fks_generator_compute_y_mismatch
    module subroutine phs_fks_generator_compute_y_test (generator, y)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(out), dimension(:):: y
    end subroutine phs_fks_generator_compute_y_test
    pure module function beta_emitter (q0, p) result (beta)
      real(default), intent(in) :: q0
      type(vector4_t), intent(in) :: p
      real(default) :: beta
    end function beta_emitter
    pure module subroutine phs_fks_generator_compute_xi_tilde (generator, r)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in) :: r
    end subroutine phs_fks_generator_compute_xi_tilde
    module subroutine phs_fks_generator_prepare_generation (generator, &
         r_in, i_phs, emitter, p_born, phs_identifiers, contributors, i_con)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), dimension(3), intent(in) :: r_in
      integer, intent(in) :: i_phs, emitter
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
      type(resonance_contributors_t), intent(in), dimension(:), optional :: &
           contributors
      integer, intent(in), optional :: i_con
    end subroutine phs_fks_generator_prepare_generation
    module subroutine phs_fks_generator_generate_fsr_from_xi_and_y &
         (generator, xi, y, &
       phi, emitter, i_phs, p_born, p_real)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in) :: xi, y, phi
      integer, intent(in) :: emitter, i_phs
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(out), dimension(:) :: p_real
    end subroutine phs_fks_generator_generate_fsr_from_xi_and_y
    module subroutine phs_fks_generator_generate_isr_from_xi_and_y &
         (generator, xi, xi_max, y, phi, i_phs, p_born, p_real)
      class(phs_fks_generator_t), intent(inout) :: generator
      real(default), intent(in) :: xi, xi_max, y, phi
      integer, intent(in) :: i_phs
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(out), dimension(:) :: p_real
    end subroutine phs_fks_generator_generate_isr_from_xi_and_y
    pure module subroutine phs_fks_generator_get_radiation_variables &
         (generator, i_phs, xi, y, phi)
      class(phs_fks_generator_t), intent(in) :: generator
      integer, intent(in) :: i_phs
      real(default), intent(out) :: xi, y
      real(default), intent(out), optional :: phi
    end subroutine phs_fks_generator_get_radiation_variables
    module subroutine phs_fks_generator_write (generator, unit)
      class(phs_fks_generator_t), intent(in) :: generator
      integer, intent(in), optional :: unit
    end subroutine phs_fks_generator_write
    module subroutine phs_fks_compute_isr_kinematics (phs, r)
      class(phs_fks_t), intent(inout) :: phs
      real(default), intent(in) :: r
    end subroutine phs_fks_compute_isr_kinematics
    module subroutine phs_fks_final (object)
      class(phs_fks_t), intent(inout) :: object
    end subroutine phs_fks_final
    module subroutine get_filtered_resonance_histories &
         (phs_config, n_in, flv_state, model, excluded_resonances, &
          resonance_histories_filtered, success)
      type(phs_fks_config_t), intent(inout) :: phs_config
      integer, intent(in) :: n_in
      integer, intent(in), dimension(:,:), allocatable :: flv_state
      type(model_t), intent(in) :: model
      type(string_t), intent(in), dimension(:), allocatable :: &
           excluded_resonances
      type(resonance_history_t), intent(out), dimension(:), &
         allocatable :: resonance_histories_filtered
      logical, intent(out) :: success
    end subroutine get_filtered_resonance_histories
  end interface

contains

  subroutine phs_fks_config_allocate_instance (phs)
    class(phs_t), intent(inout), pointer :: phs
    allocate (phs_fks_t :: phs)
  end subroutine phs_fks_config_allocate_instance


end module phs_fks

