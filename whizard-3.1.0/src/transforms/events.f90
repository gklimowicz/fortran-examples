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

module events

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use diagnostics
  use variables
  use expr_base
  use model_data
  use state_matrices, only: FM_IGNORE_HELICITY, &
       FM_SELECT_HELICITY, FM_FACTOR_HELICITY, FM_CORRELATED_HELICITY
  use particles
  use subevt_expr
  use rng_base
  use process, only: process_t
  use instances, only: process_instance_t
  use process_stacks
  use event_base
  use event_transforms

  implicit none
  private

  public :: event_t
  public :: pacify

  type :: event_config_t
     logical :: unweighted = .false.
     integer :: norm_mode = NORM_UNDEFINED
     integer :: factorization_mode = FM_IGNORE_HELICITY
     logical :: keep_correlations = .false.
     logical :: colorize_subevt = .false.
     real(default) :: sigma = 1
     integer :: n = 1
     real(default) :: safety_factor = 1
     class(expr_factory_t), allocatable :: ef_selection
     class(expr_factory_t), allocatable :: ef_reweight
     class(expr_factory_t), allocatable :: ef_analysis
   contains
     procedure :: write => event_config_write
  end type event_config_t

  type, extends (generic_event_t) :: event_t
     type(event_config_t) :: config
     type(process_t), pointer :: process => null ()
     type(process_instance_t), pointer :: instance => null ()
     class(rng_t), allocatable :: rng
     integer :: selected_i_mci = 0
     integer :: selected_i_term = 0
     integer :: selected_channel = 0
     logical :: is_complete = .false.
     class(evt_t), pointer :: transform_first => null ()
     class(evt_t), pointer :: transform_last => null ()
     type(event_expr_t) :: expr
     logical :: selection_evaluated = .false.
     logical :: passed = .false.
     real(default), allocatable :: alpha_qcd_forced
     real(default), allocatable :: scale_forced
     real(default) :: reweight = 1
     logical :: analysis_flag = .false.
     integer :: i_event = 0
   contains
     procedure :: clone => event_clone
     procedure :: final => event_final
     procedure :: write => event_write
     procedure :: basic_init => event_init
     procedure :: set_sigma => event_set_sigma
     procedure :: set_n => event_set_n
     procedure :: import_transform => event_import_transform
     procedure :: connect => event_connect
     procedure :: set_selection => event_set_selection
     procedure :: set_reweight => event_set_reweight
     procedure :: set_analysis => event_set_analysis
     procedure :: setup_expressions => event_setup_expressions
     procedure :: evaluate_transforms => event_evaluate_transforms
     procedure :: set_index => event_set_index
     procedure :: increment_index => event_increment_index
     procedure :: evaluate_expressions => event_evaluate_expressions
     procedure :: passed_selection => event_passed_selection
     procedure :: store_alt_values => event_store_alt_values
     procedure :: is_nlo => event_is_nlo
     procedure :: reset_contents => event_reset_contents
     procedure :: reset_index => event_reset_index
     procedure :: import_instance_results => event_import_instance_results
     procedure :: accept_sqme_ref => event_accept_sqme_ref
     procedure :: accept_sqme_prc => event_accept_sqme_prc
     procedure :: accept_weight_ref => event_accept_weight_ref
     procedure :: accept_weight_prc => event_accept_weight_prc
     procedure :: update_normalization => event_update_normalization
     procedure :: check => event_check
     procedure :: generate => event_generate
     procedure :: get_hard_particle_set => event_get_hard_particle_set
     procedure :: select => event_select
     procedure :: set_hard_particle_set => event_set_hard_particle_set
     procedure :: set_alpha_qcd_forced => event_set_alpha_qcd_forced
     procedure :: set_scale_forced => event_set_scale_forced
     procedure :: recalculate => event_recalculate
     procedure :: get_process_ptr => event_get_process_ptr
     procedure :: get_process_instance_ptr => event_get_process_instance_ptr
     procedure :: get_model_ptr => event_get_model_ptr
     procedure :: get_i_mci => event_get_i_mci
     procedure :: get_i_term => event_get_i_term
     procedure :: get_channel => event_get_channel
     procedure :: has_transform => event_has_transform
     procedure :: get_norm_mode => event_get_norm_mode
     procedure :: get_kinematical_weight => event_get_kinematical_weight
     procedure :: has_index => event_has_index
     procedure :: get_index => event_get_index
     procedure :: get_fac_scale => event_get_fac_scale
     procedure :: get_alpha_s => event_get_alpha_s
     procedure :: get_sqrts => event_get_sqrts
     procedure :: get_polarization => event_get_polarization
     procedure :: get_beam_file => event_get_beam_file
     procedure :: get_process_name => event_get_process_name
     procedure :: get_actual_calls_total => event_get_actual_calls_total
  end type event_t


  interface pacify
     module procedure pacify_event
  end interface pacify

  interface
    module subroutine event_config_write (object, unit, show_expressions)
      class(event_config_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_expressions
    end subroutine event_config_write
    module subroutine event_clone (event, event_new)
      class(event_t), intent(in), target :: event
      class(event_t), intent(out), target:: event_new
    end subroutine event_clone
    module subroutine event_final (object)
      class(event_t), intent(inout) :: object
    end subroutine event_final
    module subroutine event_write (object, unit, show_process, &
         show_transforms, show_decay, verbose, testflag)
      class(event_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_process, show_transforms, show_decay
      logical, intent(in), optional :: verbose
      logical, intent(in), optional :: testflag
    end subroutine event_write
    elemental module subroutine event_set_sigma (event, sigma)
      class(event_t), intent(inout) :: event
      real(default), intent(in) :: sigma
    end subroutine event_set_sigma
    elemental module subroutine event_set_n (event, n)
      class(event_t), intent(inout) :: event
      integer, intent(in) :: n
    end subroutine event_set_n
    module subroutine event_import_transform (event, evt)
      class(event_t), intent(inout) :: event
      class(evt_t), intent(inout), pointer :: evt
    end subroutine event_import_transform
    module subroutine event_connect &
         (event, process_instance, model, process_stack)
      class(event_t), intent(inout), target :: event
      type(process_instance_t), intent(in), target :: process_instance
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
    end subroutine event_connect
    module subroutine event_set_selection (event, ef_selection)
      class(event_t), intent(inout) :: event
      class(expr_factory_t), intent(in) :: ef_selection
    end subroutine event_set_selection
    module subroutine event_set_reweight (event, ef_reweight)
      class(event_t), intent(inout) :: event
      class(expr_factory_t), intent(in) :: ef_reweight
    end subroutine event_set_reweight
    module subroutine event_set_analysis (event, ef_analysis)
      class(event_t), intent(inout) :: event
      class(expr_factory_t), intent(in) :: ef_analysis
    end subroutine event_set_analysis
    module subroutine event_setup_expressions (event)
      class(event_t), intent(inout), target :: event
    end subroutine event_setup_expressions
    module subroutine event_evaluate_transforms (event, r)
      class(event_t), intent(inout) :: event
      real(default), dimension(:), intent(in), optional :: r
    end subroutine event_evaluate_transforms
    module subroutine event_set_index (event, index)
      class(event_t), intent(inout) :: event
      integer, intent(in) :: index
    end subroutine event_set_index
    module subroutine event_increment_index (event, offset)
      class(event_t), intent(inout) :: event
      integer, intent(in), optional :: offset
    end subroutine event_increment_index
    module subroutine event_evaluate_expressions (event)
      class(event_t), intent(inout) :: event
    end subroutine event_evaluate_expressions
    module function event_passed_selection (event) result (flag)
      class(event_t), intent(in) :: event
      logical :: flag
    end function event_passed_selection
    module subroutine event_store_alt_values (event)
      class(event_t), intent(inout) :: event
    end subroutine event_store_alt_values
    module function event_is_nlo (event) result (is_nlo)
      logical :: is_nlo
      class(event_t), intent(in) :: event
    end function event_is_nlo
    module subroutine event_reset_contents (event)
      class(event_t), intent(inout) :: event
    end subroutine event_reset_contents
    module subroutine event_reset_index (event)
      class(event_t), intent(inout) :: event
    end subroutine event_reset_index
    module subroutine event_import_instance_results (event)
      class(event_t), intent(inout) :: event
    end subroutine event_import_instance_results
    module subroutine event_accept_sqme_ref (event)
      class(event_t), intent(inout) :: event
    end subroutine event_accept_sqme_ref
    module subroutine event_accept_sqme_prc (event)
      class(event_t), intent(inout) :: event
    end subroutine event_accept_sqme_prc
    module subroutine event_accept_weight_ref (event)
      class(event_t), intent(inout) :: event
    end subroutine event_accept_weight_ref
    module subroutine event_accept_weight_prc (event)
      class(event_t), intent(inout) :: event
    end subroutine event_accept_weight_prc
    module subroutine event_update_normalization (event, mode_ref)
      class(event_t), intent(inout) :: event
      integer, intent(in), optional :: mode_ref
    end subroutine event_update_normalization
    module subroutine event_check (event)
      class(event_t), intent(inout) :: event
    end subroutine event_check
    module subroutine event_generate (event, i_mci, r, i_nlo)
      class(event_t), intent(inout) :: event
      integer, intent(in) :: i_mci
      real(default), dimension(:), intent(in), optional :: r
      integer, intent(in), optional :: i_nlo
    end subroutine event_generate
    module subroutine event_get_hard_particle_set (event, pset)
      class(event_t), intent(in) :: event
      type(particle_set_t), intent(out) :: pset
    end subroutine event_get_hard_particle_set
    module subroutine event_select (event, i_mci, i_term, channel)
      class(event_t), intent(inout) :: event
      integer, intent(in) :: i_mci, i_term, channel
    end subroutine event_select
    module subroutine event_set_hard_particle_set (event, particle_set)
      class(event_t), intent(inout) :: event
      type(particle_set_t), intent(in) :: particle_set
    end subroutine event_set_hard_particle_set
    module subroutine event_set_alpha_qcd_forced (event, alpha_qcd)
      class(event_t), intent(inout) :: event
      real(default), intent(in) :: alpha_qcd
    end subroutine event_set_alpha_qcd_forced
    module subroutine event_set_scale_forced (event, scale)
      class(event_t), intent(inout) :: event
      real(default), intent(in) :: scale
    end subroutine event_set_scale_forced
    module subroutine event_recalculate (event, update_sqme, weight_factor, &
         recover_beams, recover_phs, check_match, success)
      class(event_t), intent(inout) :: event
      logical, intent(in) :: update_sqme
      real(default), intent(in), optional :: weight_factor
      logical, intent(in), optional :: recover_beams
      logical, intent(in), optional :: recover_phs
      logical, intent(in), optional :: check_match
      logical, intent(out), optional :: success
    end subroutine event_recalculate
    module function event_get_process_ptr (event) result (ptr)
      class(event_t), intent(in) :: event
      type(process_t), pointer :: ptr
    end function event_get_process_ptr
    module function event_get_process_instance_ptr (event) result (ptr)
      class(event_t), intent(in) :: event
      type(process_instance_t), pointer :: ptr
    end function event_get_process_instance_ptr
    module function event_get_model_ptr (event) result (model)
      class(event_t), intent(in) :: event
      class(model_data_t), pointer :: model
    end function event_get_model_ptr
    module function event_get_i_mci (event) result (i_mci)
      class(event_t), intent(in) :: event
      integer :: i_mci
    end function event_get_i_mci
    module function event_get_i_term (event) result (i_term)
      class(event_t), intent(in) :: event
      integer :: i_term
    end function event_get_i_term
    module function event_get_channel (event) result (channel)
      class(event_t), intent(in) :: event
      integer :: channel
    end function event_get_channel
    module function event_has_transform (event) result (flag)
      class(event_t), intent(in) :: event
      logical :: flag
    end function event_has_transform
    elemental module function event_get_norm_mode (event) result (norm_mode)
      class(event_t), intent(in) :: event
      integer :: norm_mode
    end function event_get_norm_mode
    module function event_get_kinematical_weight (event) result (f)
      class(event_t), intent(in) :: event
      real(default) :: f
    end function event_get_kinematical_weight
    module function event_has_index (event) result (flag)
      class(event_t), intent(in) :: event
      logical :: flag
    end function event_has_index
    module function event_get_index (event) result (index)
      class(event_t), intent(in) :: event
      integer :: index
    end function event_get_index
    module function event_get_fac_scale (event) result (fac_scale)
      class(event_t), intent(in) :: event
      real(default) :: fac_scale
    end function event_get_fac_scale
    module function event_get_alpha_s (event) result (alpha_s)
      class(event_t), intent(in) :: event
      real(default) :: alpha_s
    end function event_get_alpha_s
    module function event_get_sqrts (event) result (sqrts)
      class(event_t), intent(in) :: event
      real(default) :: sqrts
    end function event_get_sqrts
    module function event_get_polarization (event) result (pol)
      class(event_t), intent(in) :: event
      real(default), dimension(:), allocatable :: pol
    end function event_get_polarization
    module function event_get_beam_file (event) result (file)
      class(event_t), intent(in) :: event
      type(string_t) :: file
    end function event_get_beam_file
    module function event_get_process_name (event) result (name)
      class(event_t), intent(in) :: event
      type(string_t) :: name
    end function event_get_process_name
    elemental module function event_get_actual_calls_total (event) result (n)
      class(event_t), intent(in) :: event
      integer :: n
    end function event_get_actual_calls_total
    module subroutine pacify_event (event)
      class(event_t), intent(inout) :: event
      class(evt_t), pointer :: evt
    end subroutine pacify_event
  end interface

contains

  subroutine event_init (event, var_list, n_alt)
    class(event_t), intent(out) :: event
    type(var_list_t), intent(in), optional :: var_list
    integer, intent(in), optional :: n_alt
    type(string_t) :: norm_string, mode_string
    logical :: polarized_events
    if (present (n_alt)) then
       call event%base_init (n_alt)
       call event%expr%init (n_alt)
    else
       call event%base_init (0)
    end if
    if (present (var_list)) then
       event%config%unweighted = var_list%get_lval (&
            var_str ("?unweighted"))
       norm_string = var_list%get_sval (&
            var_str ("$sample_normalization"))
       event%config%norm_mode = &
            event_normalization_mode (norm_string, event%config%unweighted)
       polarized_events = &
            var_list%get_lval (var_str ("?polarized_events"))
       if (polarized_events) then
          mode_string = &
               var_list%get_sval (var_str ("$polarization_mode"))
          select case (char (mode_string))
          case ("ignore")
             event%config%factorization_mode = FM_IGNORE_HELICITY
          case ("helicity")
             event%config%factorization_mode = FM_SELECT_HELICITY
          case ("factorized")
             event%config%factorization_mode = FM_FACTOR_HELICITY
          case ("correlated")
             event%config%factorization_mode = FM_CORRELATED_HELICITY
          case default
             call msg_fatal ("Polarization mode " &
                  // char (mode_string) // " is undefined")
          end select
       else
          event%config%factorization_mode = FM_IGNORE_HELICITY
       end if
       event%config%colorize_subevt = &
            var_list%get_lval (var_str ("?colorize_subevt"))
       if (event%config%unweighted) then
          event%config%safety_factor = var_list%get_rval (&
               var_str ("safety_factor"))
       end if
    else
       event%config%norm_mode = NORM_SIGMA
    end if
    allocate (evt_trivial_t :: event%transform_first)
    event%transform_last => event%transform_first
  end subroutine event_init


end module events
