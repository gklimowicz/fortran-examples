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

module powheg_matching

  use, intrinsic :: iso_fortran_env

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use diagnostics
  use constants, only: ZERO, ONE, TWO, THREE, FOUR, FIVE, TINY_07, PI, TWOPI
  use pdf, only: pdf_data_t
  use lorentz
  use phs_points, only: assignment(=), operator(*)
  use sm_qcd, only: qcd_t, alpha_qcd_from_scale_t, alpha_qcd_from_lambda_t
  use particles
  use grids
  use solver
  use rng_base
  use variables

  use phs_fks, only: phs_fks_generator_t, compute_dalitz_bounds, beta_emitter
  use phs_fks, only: phs_point_set_t, phs_identifier_t, phs_fks_t
  use phs_fks, only: get_xi_max_isr
  use phs_fks, only: I_XI, I_Y, I_PLUS, I_MINUS, UBF_FSR_SIMPLE, UBF_FSR_MASSIVE, UBF_FSR_MASSLESS_RECOIL, UBF_ISR
  use matching_base
  use instances, only: process_instance_t, process_instance_hook_t
  use pcm, only: pcm_nlo_t, pcm_nlo_workspace_t

  implicit none
  private

  public :: powheg_settings_t
  public :: radiation_t
  public :: process_deps_t
  public :: event_deps_t
  public :: sudakov_t
  public :: sudakov_wrapper_t
  public :: sudakov_simple_fsr_t
  public :: sudakov_eeqq_fsr_t
  public :: sudakov_massive_fsr_t
  public :: sudakov_isr_t
  public :: powheg_matching_t
  public :: powheg_matching_hook_t

  type :: powheg_settings_t
     real(default) :: pt2_min = zero
     real(default) :: lambda = zero
     integer :: size_grid_xi = 0
     integer :: size_grid_y = 0
     integer :: upper_bound_func_type = UBF_FSR_SIMPLE
     logical :: test_sudakov = .false.
     logical :: disable_sudakov = .false.
     logical :: singular_jacobian = .false.
  contains
     procedure :: init => powheg_settings_init
     procedure :: write => powheg_settings_write
  end type powheg_settings_t

  type :: radiation_t
    real(default) :: xi, y, phi, pt2
    integer :: alr
    logical :: valid = .false.
  contains
    procedure :: write => radiation_write
  end type radiation_t

  type :: process_deps_t
     real(default) :: lambda2_gen, lambda5MSB, sqrts
     integer :: n_alr
     integer :: alpha_power, alphas_power
     logical :: lab_is_cm = .true.
     type(pdf_data_t) :: pdf_data
     type(phs_identifier_t), dimension(:), allocatable :: phs_identifiers
     integer, dimension(:), allocatable :: alr_to_i_phs
     integer :: i_term_born
     integer, dimension(:), allocatable :: i_term_real
  contains
     procedure :: write => process_deps_write
  end type process_deps_t

  type :: event_deps_t
     real(default) :: s_hat
     real(default), dimension(2) :: x_born
     real(default) :: s_had
     type(phs_point_set_t) :: p_born_cms
     type(phs_point_set_t) :: p_born_lab
     type(phs_point_set_t) :: p_real_cms
     type(phs_point_set_t) :: p_real_lab
     real(default), dimension(:), allocatable :: sqme_born
  contains
     procedure :: write => event_deps_write
     procedure :: update => event_deps_update
     procedure :: set_cms => event_deps_set_cms
  end type event_deps_t

  type :: veto_counter_t
    integer :: n_ubf = 0
    integer :: n_xi_max = 0
    integer :: n_norm = 0
    integer :: n_sqme = 0
    integer :: n_veto_ubf = 0
    integer :: n_veto_xi_max = 0
    integer :: n_veto_norm = 0
    integer :: n_veto_sqme = 0
    integer :: n_fail_ubf = 0
  contains
    procedure :: record_ubf => veto_counter_record_ubf
    procedure :: record_xi_max => veto_counter_record_xi_max
    procedure :: record_norm => veto_counter_record_norm
    procedure :: record_sqme => veto_counter_record_sqme
    procedure :: record_ubf_fail => veto_counter_record_ubf_fail
    procedure :: write => veto_counter_write
  end type veto_counter_t

  type, abstract, extends (solver_function_t) :: sudakov_t
     type(process_deps_t), pointer :: process_deps => null()
     type(event_deps_t), pointer :: event_deps => null()
     type(powheg_settings_t), pointer :: powheg_settings => null()
     type(phs_fks_generator_t), pointer :: phs_fks_generator => null()
     type(qcd_t) :: qcd
     class(rng_t), pointer :: rng => null()
     real(default) :: xi2_max = zero
     real(default) :: norm_max = zero
     real(default) :: current_pt2_max = zero
     real(default) :: sum_log_rands = zero
     real(default) :: random = zero
     type(veto_counter_t) :: veto_counter
     integer :: i_phs = 0
  contains
     procedure :: write => sudakov_write
     procedure :: init => sudakov_init
     procedure :: set_normalization => sudakov_set_normalization
     procedure :: update_xi2_max => sudakov_update_xi2_max
     procedure (sudakov_upper_bound_func), deferred :: upper_bound_func
     procedure (sudakov_log_integrated_ubf), deferred :: log_integrated_ubf
     procedure (sudakov_generate_xi_and_y_and_phi), deferred :: generate_xi_and_y_and_phi
     procedure :: generate_phi => sudakov_generate_phi
     procedure (sudakov_kt2), deferred :: kt2
     procedure (sudakov_kt2_max), deferred :: kt2_max
     procedure (sudakov_reweight_ubf), deferred :: reweight_ubf
     procedure (sudakov_reweight_xi_max), deferred :: reweight_xi_max
     procedure :: alpha_s => sudakov_alpha_s
     procedure :: generate_pt2 => sudakov_generate_pt2
     procedure :: check_solution_interval => sudakov_check_solution_interval
     procedure :: generate_emission => sudakov_generate_emission
     procedure :: evaluate => sudakov_evaluate
     procedure :: associated_emitter => sudakov_associated_emitter
     procedure :: set_i_phs => sudakov_set_i_phs
     procedure :: alpha_s_rad => sudakov_alpha_s_rad
  end type sudakov_t

  type :: sudakov_wrapper_t
     class(sudakov_t), allocatable :: s
  end type sudakov_wrapper_t

  type, extends (sudakov_t) :: sudakov_simple_fsr_t
  contains
     procedure :: upper_bound_func => sudakov_simple_fsr_upper_bound_func
     procedure :: kt2 => sudakov_simple_fsr_kt2
     procedure :: kt2_max => sudakov_simple_fsr_kt2_max
     procedure :: log_integrated_ubf => sudakov_simple_fsr_log_integrated_ubf
     procedure :: reweight_ubf => sudakov_simple_fsr_reweight_ubf
     procedure :: reweight_xi_max => sudakov_simple_fsr_reweight_xi_max
     procedure :: generate_xi_and_y_and_phi => &
          sudakov_simple_fsr_generate_xi_and_y_and_phi
     procedure :: generate_xi => sudakov_simple_fsr_generate_xi
  end type sudakov_simple_fsr_t

  type, extends (sudakov_t) :: sudakov_eeqq_fsr_t
  contains
     procedure :: kt2 => sudakov_eeqq_fsr_kt2
     procedure :: kt2_max => sudakov_eeqq_fsr_kt2_max
     procedure :: upper_bound_func => sudakov_eeqq_fsr_upper_bound_func
     procedure :: log_integrated_ubf => sudakov_eeqq_fsr_log_integrated_ubf
     procedure :: reweight_ubf => sudakov_eeqq_fsr_reweight_ubf
     procedure :: reweight_xi_max => sudakov_eeqq_fsr_reweight_xi_max
     procedure :: generate_xi_and_y_and_phi => &
          sudakov_eeqq_fsr_generate_xi_and_y_and_phi
     procedure :: generate_y => sudakov_eeqq_fsr_generate_y
     procedure :: generate_xi => sudakov_eeqq_fsr_generate_xi
  end type sudakov_eeqq_fsr_t

  type, extends (sudakov_t) :: sudakov_massive_fsr_t
    real(default) :: z, z1, z2 = 0._default
    real(default) :: xi_1, xi_min, xi_m = 0._default
    real(default) :: xi_max_extended = 1._default
  contains
    procedure :: compute_xi_max_extended &
       => sudakov_massive_fsr_compute_xi_max_extended
    procedure :: generate_xi => sudakov_massive_fsr_generate_xi
    procedure :: generate_xi_and_y_and_phi => &
         sudakov_massive_fsr_generate_xi_and_y_and_phi
    procedure :: kt2 => sudakov_massive_fsr_kt2
    procedure :: kt2_max => sudakov_massive_fsr_kt2_max
    procedure :: upper_bound_func => sudakov_massive_fsr_upper_bound_func
    procedure :: log_integrated_ubf => sudakov_massive_fsr_log_integrated_ubf
    procedure :: reweight_ubf => sudakov_massive_fsr_reweight_ubf
    procedure :: reweight_xi_max => sudakov_massive_fsr_reweight_xi_max
  end type sudakov_massive_fsr_t

  type, extends (sudakov_t) :: sudakov_isr_t
    real(default) :: xi_max_extended = 1._default
    integer :: ubf_variant = 2
  contains
     procedure :: kt2 => sudakov_isr_kt2
     procedure :: kt2_max => sudakov_isr_kt2_max
     procedure :: upper_bound_func => sudakov_isr_upper_bound_func
     procedure :: log_integrated_ubf => sudakov_isr_log_integrated_ubf
     procedure :: reweight_ubf => sudakov_isr_reweight_ubf
     procedure :: reweight_xi_max => sudakov_isr_reweight_xi_max
     procedure :: compute_xi_max_extended &
        => sudakov_isr_compute_xi_max_extended
     procedure :: generate_xi_and_y_and_phi => &
          sudakov_isr_generate_xi_and_y_and_phi
     procedure :: generate_y => sudakov_isr_generate_y
     procedure :: generate_xi => sudakov_isr_generate_xi
     procedure :: compute_xi2_max => sudakov_isr_compute_xi2_max
  end type sudakov_isr_t

  type, extends(matching_t) :: powheg_matching_t
     type(grid_t) :: grid
     type(phs_fks_generator_t) :: phs_fks_generator
     type(powheg_settings_t) :: settings
     type(event_deps_t) :: event_deps
     type(process_deps_t) :: process_deps
     type(sudakov_wrapper_t), dimension(:), allocatable :: sudakov
     integer :: n_emissions = 0
     logical :: active = .true.
   contains
     procedure :: get_method => powheg_matching_get_method
     procedure :: before_shower => powheg_matching_before_shower
     procedure :: first_event => powheg_matching_first_event
     procedure :: after_shower => powheg_matching_after_shower
     procedure :: write => powheg_write
     procedure :: final => powheg_matching_final
     procedure :: setup_grids => powheg_matching_setup_grids
     procedure :: setup_sudakovs => powheg_matching_setup_sudakovs
     procedure :: init => powheg_matching_init
     generic :: update => update_momenta, &
                          update_particle_set
     procedure :: update_momenta => powheg_matching_update_momenta
     procedure :: update_particle_set => powheg_matching_update_particle_set
     procedure :: update_event_deps => powheg_matching_update_event_deps
     procedure :: boost_preal_to_lab_frame => &
          powheg_matching_boost_preal_to_lab_frame
     procedure :: boost_preal_to_cms => powheg_matching_boost_preal_to_cms
     procedure :: reweight_matrix_elements => &
          powheg_matching_reweight_matrix_elements
     procedure :: compute_sqme_real => powheg_matching_compute_sqme_real
     procedure :: set_scale => powheg_matching_set_scale
     procedure :: update_sudakovs => powheg_matching_update_sudakovs
     procedure :: import_norms_from_grid => powheg_matching_import_norms_from_grid
     procedure :: save_grids => powheg_matching_save_grids
     procedure :: load_grids => powheg_matching_load_grids
     procedure :: check_grids => powheg_matching_check_grids
     procedure :: generate_emission => powheg_matching_generate_emission
     procedure :: build_particle_set => powheg_matching_build_particle_set
     procedure :: reweight_norm => powheg_matching_reweight_norm
     procedure :: norm_from_xi_and_y => powheg_matching_norm_from_xi_and_y
     procedure :: prepare_for_events => powheg_matching_prepare_for_events
     procedure :: compute_lambda2_gen => powheg_matching_compute_lambda2_gen
     procedure :: compute_lambda5MSB => powheg_matching_compute_lambda5MSB
     procedure :: setup_nlo_environment => powheg_matching_setup_nlo_environment
     procedure :: copy_momenta => powheg_matching_copy_momenta
     procedure :: test_sudakov => powheg_test_sudakov
  end type powheg_matching_t

  type, extends(process_instance_hook_t) :: powheg_matching_hook_t
     type(string_t) :: process_name
     type(powheg_matching_t) :: powheg
  contains
       procedure :: init => powheg_matching_hook_init
       procedure :: final => powheg_matching_hook_final
       procedure :: evaluate => powheg_matching_hook_evaluate
  end type powheg_matching_hook_t


  abstract interface
     pure function sudakov_upper_bound_func (sudakov, xi, y, alpha_s) result (u)
       import
       real(default) :: u
       class(sudakov_t), intent(in) :: sudakov
       real(default), intent(in) :: xi, y, alpha_s
     end function sudakov_upper_bound_func
  end interface

  abstract interface
     pure function sudakov_log_integrated_ubf (sudakov, pt2) result (y)
       import
       real(default) :: y
       class(sudakov_t), intent(in) :: sudakov
       real(default), intent(in) :: pt2
     end function sudakov_log_integrated_ubf
  end interface

  abstract interface
     subroutine sudakov_generate_xi_and_y_and_phi (sudakov, r)
       import
       class(sudakov_t), intent(inout) :: sudakov
       type(radiation_t), intent(inout) :: r
     end subroutine sudakov_generate_xi_and_y_and_phi
  end interface

  abstract interface
     function sudakov_kt2 (sudakov, xi, y) result (kt2)
       import
       real(default) :: kt2
       class(sudakov_t), intent(in) :: sudakov
       real(default), intent(in) :: xi, y
     end function sudakov_kt2
  end interface

  abstract interface
     pure function sudakov_kt2_max (sudakov) result (kt2_max)
        import
        real(default) :: kt2_max
        class(sudakov_t), intent(in) :: sudakov
     end function sudakov_kt2_max
  end interface

  abstract interface
     function sudakov_reweight_ubf (sudakov, pt2) result (accepted)
       import
       logical :: accepted
       class(sudakov_t), intent(inout) :: sudakov
       real(default), intent(in) :: pt2
    end function sudakov_reweight_ubf
  end interface

  abstract interface
     function sudakov_reweight_xi_max (sudakov, xi) result (accepted)
       import
       logical :: accepted
       class(sudakov_t), intent(in) :: sudakov
       real(default), intent(in) :: xi
     end function sudakov_reweight_xi_max
  end interface


  interface
    module subroutine powheg_settings_init (settings, var_list)
      class(powheg_settings_t), intent(out) :: settings
      type(var_list_t), intent(in) :: var_list
    end subroutine powheg_settings_init
    module subroutine powheg_settings_write (powheg_settings, unit)
      class(powheg_settings_t), intent(in) :: powheg_settings
      integer, intent(in), optional :: unit
    end subroutine powheg_settings_write
    module subroutine radiation_write (radiation, unit)
      class(radiation_t), intent(in) :: radiation
      integer, intent(in), optional :: unit
    end subroutine radiation_write
    module subroutine process_deps_write (process_deps, unit)
      class(process_deps_t), intent(in) :: process_deps
      integer, intent(in), optional :: unit
    end subroutine process_deps_write
    module subroutine event_deps_write (event_deps, unit)
      class(event_deps_t), intent(in) :: event_deps
      integer, intent(in), optional :: unit
    end subroutine event_deps_write
    module subroutine event_deps_update (event_deps, sqme_born, &
         p_born, x_born, lt_lab_to_cms)
      class(event_deps_t), intent(inout) :: event_deps
      real(default), dimension(:), intent(in) :: sqme_born
      type(vector4_t), dimension(:), intent(in) :: p_born
      real(default), dimension(2), intent(in) :: x_born
      type(lorentz_transformation_t), intent(in), optional :: lt_lab_to_cms
    end subroutine event_deps_update
    module subroutine event_deps_set_cms (event_deps, lt_lab_to_cms)
      class(event_deps_t), intent(inout) :: event_deps
      type(lorentz_transformation_t), intent(in), optional :: lt_lab_to_cms
    end subroutine event_deps_set_cms
    pure module subroutine veto_counter_record_ubf (counter, vetoed)
      class(veto_counter_t), intent(inout) :: counter
      logical, intent(in) :: vetoed
    end subroutine veto_counter_record_ubf
    module subroutine veto_counter_record_xi_max (counter, vetoed)
      class(veto_counter_t), intent(inout) :: counter
      logical, intent(in) :: vetoed
    end subroutine veto_counter_record_xi_max
    module subroutine veto_counter_record_norm (counter, vetoed)
      class(veto_counter_t), intent(inout) :: counter
      logical, intent(in) :: vetoed
    end subroutine veto_counter_record_norm
    module subroutine veto_counter_record_sqme (counter, vetoed)
      class(veto_counter_t), intent(inout) :: counter
      logical, intent(in) :: vetoed
    end subroutine veto_counter_record_sqme
    module subroutine veto_counter_record_ubf_fail (counter)
      class(veto_counter_t), intent(inout) :: counter
    end subroutine veto_counter_record_ubf_fail
    module subroutine veto_counter_write (counter, unit)
      class(veto_counter_t), intent(in) :: counter
      integer, intent(in), optional :: unit
    end subroutine veto_counter_write
    module subroutine sudakov_write (sudakov, unit)
      class(sudakov_t), intent(in) :: sudakov
      integer, intent(in), optional :: unit
    end subroutine sudakov_write
    module subroutine sudakov_init (sudakov, process_deps, event_deps, &
           powheg_settings, qcd, phs_fks_generator, rng)
      class(sudakov_t), intent(out) :: sudakov
      type(process_deps_t), target, intent(in) :: process_deps
      type(event_deps_t), target, intent(in) :: event_deps
      type(powheg_settings_t), target, intent(in) :: powheg_settings
      type(qcd_t), intent(in) :: qcd
      type(phs_fks_generator_t), target, intent(in) :: phs_fks_generator
      class(rng_t), target, intent(in), optional :: rng
    end subroutine sudakov_init
    pure module subroutine sudakov_set_normalization (sudakov, norm_max)
      class(sudakov_t), intent(inout) :: sudakov
      real(default), intent(in) :: norm_max
    end subroutine sudakov_set_normalization
    pure module subroutine sudakov_update_xi2_max (sudakov, xi2_max)
      class(sudakov_t), intent(inout) :: sudakov
      real(default), intent(in) :: xi2_max
    end subroutine sudakov_update_xi2_max
    module subroutine sudakov_generate_phi (sudakov, r)
      class(sudakov_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_generate_phi
    module function sudakov_alpha_s (sudakov, kT2, use_correct, improve_nll) result (a)
      real(default) :: a
      class(sudakov_t), intent(in) :: sudakov
      real(default), intent(in) :: kT2
      logical, intent(in), optional :: use_correct, improve_nll
    end function sudakov_alpha_s
    module function sudakov_generate_pt2 (sudakov) result (pt2)
      real(default) :: pt2
      class(sudakov_t), intent(inout) :: sudakov
    end function sudakov_generate_pt2
    module subroutine sudakov_check_solution_interval (sudakov)
      class(sudakov_t), intent(inout) :: sudakov
    end subroutine sudakov_check_solution_interval
    module subroutine sudakov_generate_emission (sudakov, r, r_max)
      class(sudakov_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
      type(radiation_t), intent(in) :: r_max
    end subroutine sudakov_generate_emission
    module function sudakov_evaluate (solver_f, x) result (f)
      complex(default) :: f
      class(sudakov_t), intent(in) :: solver_f
      real(default), intent(in) :: x
    end function sudakov_evaluate
    elemental module function sudakov_associated_emitter &
         (sudakov) result (emitter)
      integer :: emitter
      class(sudakov_t), intent(in) :: sudakov
    end function sudakov_associated_emitter
    module subroutine sudakov_set_i_phs (sudakov, alr)
      class(sudakov_t), intent(inout) :: sudakov
      integer, intent(in) :: alr
    end subroutine sudakov_set_i_phs
    pure module function sudakov_simple_fsr_upper_bound_func &
         (sudakov, xi, y, alpha_s) result (u)
      real(default) :: u
      class(sudakov_simple_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi, y, alpha_s
    end function sudakov_simple_fsr_upper_bound_func
    module function sudakov_simple_fsr_kt2 (sudakov, xi, y) result (kt2)
      real(default) :: kt2
      class(sudakov_simple_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi, y
    end function sudakov_simple_fsr_kt2
    pure module function sudakov_simple_fsr_kt2_max (sudakov) result (kt2_max)
      real(default) :: kt2_max
      class(sudakov_simple_fsr_t), intent(in) :: sudakov
    end function sudakov_simple_fsr_kt2_max
    pure module function sudakov_simple_fsr_log_integrated_ubf &
         (sudakov, pt2) result (log_sudakov)
      real(default) :: log_sudakov
      class(sudakov_simple_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: pt2
    end function sudakov_simple_fsr_log_integrated_ubf
    module function sudakov_simple_fsr_reweight_ubf &
         (sudakov, pt2) result (accepted)
      logical :: accepted
      class(sudakov_simple_fsr_t), intent(inout) :: sudakov
      real(default), intent(in) :: pt2
      real(default) :: alpha_s_true, alpha_s_rad
      logical :: alpha_s_equal
    end function sudakov_simple_fsr_reweight_ubf
    module function sudakov_simple_fsr_reweight_xi_max &
         (sudakov, xi) result (accepted)
      logical :: accepted
      class(sudakov_simple_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi
    end function sudakov_simple_fsr_reweight_xi_max
    module subroutine sudakov_simple_fsr_generate_xi_and_y_and_phi (sudakov, r)
      class(sudakov_simple_fsr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_simple_fsr_generate_xi_and_y_and_phi
    module subroutine sudakov_simple_fsr_generate_xi (sudakov, r)
      class(sudakov_simple_fsr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_simple_fsr_generate_xi
    module function sudakov_eeqq_fsr_kt2 (sudakov, xi, y) result (kt2)
      real(default) :: kt2
      class(sudakov_eeqq_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi, y
    end function sudakov_eeqq_fsr_kt2
    pure module function sudakov_eeqq_fsr_kt2_max (sudakov) result (kt2_max)
      real(default) :: kt2_max
      class(sudakov_eeqq_fsr_t), intent(in) :: sudakov
    end function sudakov_eeqq_fsr_kt2_max
    pure module function sudakov_eeqq_fsr_upper_bound_func &
         (sudakov, xi, y, alpha_s) result (u)
      real(default) :: u
      class(sudakov_eeqq_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi, y, alpha_s
    end function sudakov_eeqq_fsr_upper_bound_func
    pure module function sudakov_eeqq_fsr_log_integrated_ubf &
         (sudakov, pt2) result (log_sudakov)
      real(default) :: log_sudakov
      class(sudakov_eeqq_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: pt2
    end function sudakov_eeqq_fsr_log_integrated_ubf
    module function sudakov_eeqq_fsr_reweight_ubf &
         (sudakov, pt2) result (accepted)
      logical :: accepted
      class(sudakov_eeqq_fsr_t), intent(inout) :: sudakov
      real(default), intent(in) :: pt2
    end function sudakov_eeqq_fsr_reweight_ubf
    module function sudakov_eeqq_fsr_reweight_xi_max &
         (sudakov, xi) result (accepted)
      logical :: accepted
      class(sudakov_eeqq_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi
    end function sudakov_eeqq_fsr_reweight_xi_max
    module subroutine sudakov_eeqq_fsr_generate_xi_and_y_and_phi (sudakov, r)
      class(sudakov_eeqq_fsr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_eeqq_fsr_generate_xi_and_y_and_phi
    module subroutine sudakov_eeqq_fsr_generate_y (sudakov, r)
      class(sudakov_eeqq_fsr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_eeqq_fsr_generate_y
    module subroutine sudakov_eeqq_fsr_generate_xi (sudakov, r)
      class(sudakov_eeqq_fsr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_eeqq_fsr_generate_xi
    module subroutine sudakov_massive_fsr_compute_xi_max_extended &
         (sudakov)
      class(sudakov_massive_fsr_t), intent(inout) :: sudakov
    end subroutine sudakov_massive_fsr_compute_xi_max_extended
    module subroutine sudakov_massive_fsr_generate_xi (sudakov, r)
      class(sudakov_massive_fsr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_massive_fsr_generate_xi
    module subroutine sudakov_massive_fsr_generate_xi_and_y_and_phi (sudakov, r)
      class(sudakov_massive_fsr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_massive_fsr_generate_xi_and_y_and_phi
    module function sudakov_massive_fsr_kt2 (sudakov, xi, y) result (kt2)
      real(default) :: kt2
      class(sudakov_massive_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi, y
    end function sudakov_massive_fsr_kt2
    pure module function sudakov_massive_fsr_kt2_max (sudakov) result (kt2_max)
      real(default) :: kt2_max
      class(sudakov_massive_fsr_t), intent(in) :: sudakov
    end function sudakov_massive_fsr_kt2_max
    pure module function sudakov_massive_fsr_upper_bound_func &
         (sudakov, xi, y, alpha_s) result (u)
      real(default) :: u
      class(sudakov_massive_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi, y, alpha_s
    end function sudakov_massive_fsr_upper_bound_func
    pure module function sudakov_massive_fsr_log_integrated_ubf &
         (sudakov, pt2) result (log_sudakov)
      real(default) :: log_sudakov
      class(sudakov_massive_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: pt2
    end function sudakov_massive_fsr_log_integrated_ubf
    module function sudakov_massive_fsr_reweight_ubf &
         (sudakov, pt2) result (accepted)
      logical :: accepted
      class(sudakov_massive_fsr_t), intent(inout) :: sudakov
      real(default), intent(in) :: pt2
    end function sudakov_massive_fsr_reweight_ubf
    module function sudakov_massive_fsr_reweight_xi_max &
         (sudakov, xi) result (accepted)
      logical :: accepted
      class(sudakov_massive_fsr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi
    end function sudakov_massive_fsr_reweight_xi_max
    module function sudakov_isr_kt2 (sudakov, xi, y) result (kt2)
      real(default) :: kt2
      class(sudakov_isr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi, y
    end function sudakov_isr_kt2
    pure module function sudakov_isr_kt2_max (sudakov) result (kt2_max)
      real(default) :: kt2_max
      class(sudakov_isr_t), intent(in) :: sudakov
    end function sudakov_isr_kt2_max
    pure module function sudakov_isr_upper_bound_func &
         (sudakov, xi, y, alpha_s) result (u)
      real(default) :: u
      class(sudakov_isr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi, y, alpha_s
    end function sudakov_isr_upper_bound_func
    pure module function sudakov_isr_log_integrated_ubf &
         (sudakov, pt2) result (log_sudakov)
      real(default) :: log_sudakov
      class(sudakov_isr_t), intent(in) :: sudakov
      real(default), intent(in) :: pt2
    end function sudakov_isr_log_integrated_ubf
    module function sudakov_isr_reweight_ubf (sudakov, pt2) result (accepted)
      logical :: accepted
      class(sudakov_isr_t), intent(inout) :: sudakov
      real(default), intent(in) :: pt2
    end function sudakov_isr_reweight_ubf
    module function sudakov_isr_reweight_xi_max (sudakov, xi) result (accepted)
      logical :: accepted
      class(sudakov_isr_t), intent(in) :: sudakov
      real(default), intent(in) :: xi
    end function sudakov_isr_reweight_xi_max
    module subroutine sudakov_isr_compute_xi_max_extended (sudakov)
      class(sudakov_isr_t), intent(inout) :: sudakov
      real(default) :: s_hat, s_had
    end subroutine sudakov_isr_compute_xi_max_extended
    module subroutine sudakov_isr_generate_xi_and_y_and_phi (sudakov, r)
      class(sudakov_isr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_isr_generate_xi_and_y_and_phi
    module subroutine sudakov_isr_generate_y (sudakov, r)
      class(sudakov_isr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_isr_generate_y
    module subroutine sudakov_isr_generate_xi (sudakov, r)
      class(sudakov_isr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_isr_generate_xi
    module subroutine sudakov_isr_compute_xi2_max (sudakov, r)
      class(sudakov_isr_t), intent(inout) :: sudakov
      type(radiation_t), intent(inout) :: r
    end subroutine sudakov_isr_compute_xi2_max
    module function powheg_matching_get_method (matching) result (method)
      type(string_t) :: method
      class(powheg_matching_t), intent(in) :: matching
    end function powheg_matching_get_method
    module subroutine powheg_matching_before_shower &
           (matching, particle_set, vetoed)
      class(powheg_matching_t), intent(inout) :: matching
      type(particle_set_t), intent(inout) :: particle_set
      logical, intent(out) :: vetoed
    end subroutine powheg_matching_before_shower
    module subroutine powheg_matching_first_event (matching)
      class(powheg_matching_t),  intent(inout), target :: matching
    end subroutine powheg_matching_first_event
    module subroutine powheg_matching_after_shower &
         (matching, particle_set, vetoed)
      class(powheg_matching_t), intent(inout) :: matching
      type(particle_set_t), intent(inout) :: particle_set
      logical, intent(out) :: vetoed
    end subroutine powheg_matching_after_shower
    module subroutine powheg_write (matching, unit)
      class(powheg_matching_t), intent(in) :: matching
      integer, intent(in), optional :: unit
    end subroutine powheg_write
    module subroutine powheg_matching_final (matching)
      class(powheg_matching_t), intent(in) :: matching
    end subroutine powheg_matching_final
    module subroutine powheg_matching_setup_grids (matching)
      class(powheg_matching_t), intent(inout), target :: matching
    end subroutine powheg_matching_setup_grids
    module subroutine powheg_matching_init (matching, var_list, process_name)
      class(powheg_matching_t), intent(out) :: matching
      type(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: process_name
    end subroutine powheg_matching_init
    module subroutine powheg_matching_update_momenta (powheg, p_born)
      class(powheg_matching_t), intent(inout) :: powheg
      type(vector4_t), dimension(:), intent(in) :: p_born
    end subroutine powheg_matching_update_momenta
    module subroutine powheg_matching_update_particle_set (powheg, particle_set)
      class(powheg_matching_t), intent(inout) :: powheg
      type(particle_set_t), intent(in) :: particle_set
    end subroutine powheg_matching_update_particle_set
    module subroutine powheg_matching_update_event_deps &
         (powheg, sqme_born, p_born, x_born, lt_lab_to_cms)
      class(powheg_matching_t), intent(inout) :: powheg
      real(default), dimension(:), intent(in) :: sqme_born
      type(vector4_t), dimension(:), intent(in) :: p_born
      real(default), dimension(2), intent(in) :: x_born
      type(lorentz_transformation_t), intent(in), optional :: lt_lab_to_cms
    end subroutine powheg_matching_update_event_deps
    module subroutine powheg_matching_boost_preal_to_lab_frame (powheg, i_phs)
      class(powheg_matching_t), intent(inout) :: powheg
      integer, intent(in) :: i_phs
    end subroutine powheg_matching_boost_preal_to_lab_frame
    module subroutine powheg_matching_boost_preal_to_cms (powheg, i_phs)
      class(powheg_matching_t), intent(inout) :: powheg
      integer, intent(in) :: i_phs
    end subroutine powheg_matching_boost_preal_to_cms
    module function powheg_matching_reweight_matrix_elements &
         (powheg, r) result (accepted)
      logical :: accepted
      class(powheg_matching_t), intent(inout) :: powheg
      type(radiation_t), intent(in) :: r
    end function powheg_matching_reweight_matrix_elements
    module function powheg_matching_compute_sqme_real &
         (powheg, alr, scale) result (sqme)
      real(default) :: sqme
      class(powheg_matching_t), intent(inout) :: powheg
      integer, intent(in) :: alr
      real(default), intent(in), optional :: scale
    end function powheg_matching_compute_sqme_real
    module subroutine powheg_matching_set_scale (powheg, pT2)
      class(powheg_matching_t), intent(inout) :: powheg
      real(default), intent(in) :: pT2
    end subroutine powheg_matching_set_scale
    module subroutine powheg_matching_update_sudakovs (powheg, alr, i_phs, y)
      class(powheg_matching_t), intent(inout) :: powheg
      integer, intent(in) :: alr, i_phs
      real(default), intent(in) :: y
    end subroutine powheg_matching_update_sudakovs
    module subroutine powheg_matching_import_norms_from_grid (powheg)
      class(powheg_matching_t), intent(inout) :: powheg
    end subroutine powheg_matching_import_norms_from_grid
    module subroutine powheg_matching_save_grids (powheg)
      class(powheg_matching_t), intent(inout) :: powheg
    end subroutine powheg_matching_save_grids
    module subroutine powheg_matching_load_grids (powheg)
      class(powheg_matching_t), intent(inout) :: powheg
    end subroutine powheg_matching_load_grids
    module function powheg_matching_check_grids (powheg) result (ok)
      logical :: ok
      class(powheg_matching_t), intent(in) :: powheg
    end function powheg_matching_check_grids
    module subroutine powheg_matching_generate_emission &
         (powheg, particle_set, pt2_generated)
      class(powheg_matching_t), intent(inout) :: powheg
      type(particle_set_t), intent(inout), optional :: particle_set
      real(default), intent(out), optional :: pt2_generated
    end subroutine powheg_matching_generate_emission
    module subroutine powheg_matching_build_particle_set &
         (powheg, particle_set, p_real, alr, y)
      class(powheg_matching_t), intent(inout) :: powheg
      type(particle_set_t), intent(inout) :: particle_set
      type(vector4_t), dimension(:), intent(in) :: p_real
      integer, intent(in) :: alr
      real(default), intent(in) :: y
    end subroutine powheg_matching_build_particle_set
    module function powheg_matching_reweight_norm (powheg, r) result (accepted)
      logical :: accepted
      class(powheg_matching_t), intent(inout) :: powheg
      type(radiation_t), intent(in) :: r
    end function powheg_matching_reweight_norm
    module function powheg_matching_norm_from_xi_and_y &
         (powheg, r) result (norm_true)
      real(default) :: norm_true
      class(powheg_matching_t), intent(inout) :: powheg
      type(radiation_t), intent(in) :: r
    end function powheg_matching_norm_from_xi_and_y
    module subroutine powheg_matching_prepare_for_events (matching)
      class(powheg_matching_t), intent(inout), target :: matching
    end subroutine powheg_matching_prepare_for_events
    module subroutine powheg_matching_compute_lambda2_gen (matching)
      class(powheg_matching_t), intent(inout) :: matching
    end subroutine powheg_matching_compute_lambda2_gen
    module subroutine powheg_matching_compute_lambda5MSB (matching)
      class(powheg_matching_t), intent(inout) :: matching
    end subroutine powheg_matching_compute_lambda5MSB
    module subroutine powheg_matching_setup_nlo_environment (matching)
      class(powheg_matching_t), intent(inout) :: matching
    end subroutine powheg_matching_setup_nlo_environment
    module subroutine powheg_matching_copy_momenta (matching, i_phs)
      class(powheg_matching_t), intent(inout) :: matching
      integer, intent(in) :: i_phs
    end subroutine powheg_matching_copy_momenta
    elemental module function sudakov_alpha_s_rad (sudakov, scale2) result (alpha_s_rad)
      real(default) :: alpha_s_rad
      class(sudakov_t), intent(in) :: sudakov
      real(default), intent(in) :: scale2
    end function sudakov_alpha_s_rad
    module subroutine powheg_matching_hook_init (hook, var_list, &
         instance, pdf_data)
      class(powheg_matching_hook_t), intent(inout), target :: hook
      type(var_list_t), intent(in) :: var_list
      class(process_instance_t), intent(in), target :: instance
      type(pdf_data_t), intent(in), optional :: pdf_data
    end subroutine powheg_matching_hook_init
    module subroutine powheg_matching_hook_final (hook)
      class(powheg_matching_hook_t), intent(inout) :: hook
    end subroutine powheg_matching_hook_final
    module subroutine powheg_matching_hook_evaluate (hook, instance)
      class(powheg_matching_hook_t), intent(inout) :: hook
      class(process_instance_t), intent(in), target :: instance
    end subroutine powheg_matching_hook_evaluate
    module subroutine powheg_test_sudakov (powheg)
      class(powheg_matching_t), intent(inout) :: powheg
    end subroutine powheg_test_sudakov
  end interface

contains

  subroutine powheg_matching_setup_sudakovs (powheg)
    class(powheg_matching_t), intent(inout), target :: powheg
    integer :: alr, ubf_type, n_in, emitter
    logical :: is_fsr, is_massive, is_eeqq
    allocate (powheg%sudakov (powheg%process_deps%n_alr))
    do alr = 1, powheg%process_deps%n_alr
       select type (pcm => powheg%process_instance%pcm)
       type is (pcm_nlo_t)
          associate(reg_data => pcm%region_data, &
               phs => powheg%phs_fks_generator)
               n_in = reg_data%get_n_in()
               emitter = reg_data%get_emitter (alr)
               is_fsr = emitter > n_in
               if (is_fsr) then
                  is_massive = phs%is_massive (emitter)
               else
                  if (emitter == 0) then
                     is_massive = phs%is_massive (1) .or. phs%is_massive (2)
                  else
                     is_massive = phs%is_massive (emitter)
                  end if
                  if (is_massive) call msg_bug ("setup_sudakovs: ISR " // &
                       "off massive emitters not implemented.")
               end if
               is_eeqq = n_in == 2 .and. &
                    reg_data%n_legs_born == 4 .and. &
                    .not. phs%is_massive(3) .and. &
                    .not. phs%is_massive(4)
               !  (PS 2021-05-28) This includes FSR regions of pp -> jj.
          end associate
       end select
       if (is_fsr) then
          if (is_eeqq) then
             ubf_type = UBF_FSR_MASSLESS_RECOIL
          else if (is_massive) then
             ubf_type = UBF_FSR_MASSIVE
          else
             ubf_type = powheg%settings%upper_bound_func_type
          end if
          select case (ubf_type)
          case (UBF_FSR_SIMPLE)
             allocate (sudakov_simple_fsr_t :: powheg%sudakov(alr)%s)
          case (UBF_FSR_MASSLESS_RECOIL)
             allocate (sudakov_eeqq_fsr_t :: powheg%sudakov(alr)%s)
          case (UBF_FSR_MASSIVE)
             allocate (sudakov_massive_fsr_t :: powheg%sudakov(alr)%s)
          case default
             call msg_fatal ("powheg_setup_sudakovs: Please choose " // &
                  "upper bounding function!")
          end select
       else
          allocate (sudakov_isr_t :: powheg%sudakov(alr)%s)
       end if
       if (allocated (powheg%rng)) then
          !!! generator mode
          call powheg%sudakov(alr)%s%init (powheg%process_deps, &
               powheg%event_deps, powheg%settings, &
               powheg%qcd, powheg%phs_fks_generator, powheg%rng)
       else
          !!! lookup mode
          call powheg%sudakov(alr)%s%init (powheg%process_deps, &
               powheg%event_deps, powheg%settings, &
               powheg%qcd, powheg%phs_fks_generator)
       end if
       call powheg%sudakov(alr)%s%set_i_phs (alr)
    end do
  end subroutine powheg_matching_setup_sudakovs


end module powheg_matching
