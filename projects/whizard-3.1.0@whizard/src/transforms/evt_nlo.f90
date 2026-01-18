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

module evt_nlo

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use constants
  use phs_points, only: phs_point_t
  use phs_points, only: assignment(=), operator(*), size
  use sm_qcd
  use model_data
  use particles
  use instances, only: process_instance_t
  use process_stacks
  use event_transforms
  use quantum_numbers, only: quantum_numbers_t

  use phs_fks, only: phs_fks_t, phs_fks_generator_t
  use phs_fks, only: phs_identifier_t, phs_point_set_t
  use resonances, only: resonance_contributors_t
  use fks_regions, only: region_data_t

  implicit none
  private

  public :: evt_nlo_t

  integer, parameter, public :: EVT_NLO_UNDEFINED = 0
  integer, parameter, public :: EVT_NLO_SEPARATE_BORNLIKE = 1
  integer, parameter, public :: EVT_NLO_SEPARATE_REAL = 2
  integer, parameter, public :: EVT_NLO_COMBINED = 3

  type :: nlo_event_deps_t
     logical :: lab_is_cm = .true.
     type(phs_point_set_t) :: p_born_cms
     type(phs_point_set_t) :: p_born_lab
     type(phs_point_set_t) :: p_real_cms
     type(phs_point_set_t) :: p_real_lab
     type(resonance_contributors_t), dimension(:), allocatable :: contributors
     type(phs_identifier_t), dimension(:), allocatable :: phs_identifiers
     integer, dimension(:), allocatable :: alr_to_i_con
     integer :: n_phs = 0
  end type nlo_event_deps_t

  type, extends (evt_t) :: evt_nlo_t
    type(phs_fks_generator_t) :: phs_fks_generator
    real(default) :: sqme_rad = zero
    integer :: i_evaluation = 0
    type(particle_set_t), dimension(:), allocatable :: particle_set_nlo
    type(qcd_t) :: qcd
    type(nlo_event_deps_t) :: event_deps
    integer :: mode = EVT_NLO_UNDEFINED
    integer, dimension(:), allocatable :: &
       i_evaluation_to_i_phs, i_evaluation_to_emitter, &
       i_evaluation_to_i_term
    logical :: keep_failed_events = .false.
    integer :: selected_i_flv = 0
  contains
    procedure :: write_name => evt_nlo_write_name
    procedure :: write => evt_nlo_write
    procedure :: connect => evt_nlo_connect
    procedure :: set_i_evaluation_mappings => evt_nlo_set_i_evaluation_mappings
    procedure :: get_i_phs => evt_nlo_get_i_phs
    procedure :: get_emitter => evt_nlo_get_emitter
    procedure :: get_i_term => evt_nlo_get_i_term
    procedure :: generate_weighted => evt_nlo_generate_weighted
    procedure :: reset_phs_identifiers => evt_nlo_reset_phs_identifiers
    procedure :: connected_set_real_IS_momenta => &
         evt_nlo_connected_set_real_IS_momenta
    procedure :: make_particle_set => evt_nlo_make_particle_set
    procedure :: evaluate_real_kinematics => evt_nlo_evaluate_real_kinematics
    procedure :: compute_subtraction_sqmes => evt_nlo_compute_subtraction_sqmes
    procedure :: compute_real => evt_nlo_compute_real
    procedure :: compute_all_sqme_rad => evt_nlo_compute_all_sqme_rad
    procedure :: boost_to_cms => evt_nlo_boost_to_cms
    procedure :: boost_to_lab => evt_nlo_boost_to_lab
    procedure :: setup_general_event_kinematics => &
         evt_nlo_setup_general_event_kinematics
    procedure :: setup_real_event_kinematics => &
         evt_nlo_setup_real_event_kinematics
    procedure :: set_mode => evt_nlo_set_mode
    procedure :: is_valid_event => evt_nlo_is_valid_event
    procedure :: get_selected_quantum_numbers => &
         evt_nlo_get_selected_quantum_numbers
    procedure :: prepare_new_event => evt_nlo_prepare_new_event
  end type evt_nlo_t

  type :: registered_triple_t
    integer, dimension(2) :: phs_em
    type(registered_triple_t), pointer :: next => null ()
  end type registered_triple_t

  interface
    module subroutine evt_nlo_write_name (evt, unit)
      class(evt_nlo_t), intent(in) :: evt
      integer, intent(in), optional :: unit
    end subroutine evt_nlo_write_name
    module subroutine evt_nlo_write (evt, unit, verbose, more_verbose, testflag)
      class(evt_nlo_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, more_verbose, testflag
    end subroutine evt_nlo_write
    module subroutine evt_nlo_connect &
         (evt, process_instance, model, process_stack)
      class(evt_nlo_t), intent(inout), target :: evt
      type(process_instance_t), intent(in), target :: process_instance
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
    end subroutine evt_nlo_connect
    module subroutine evt_nlo_set_i_evaluation_mappings &
         (evt, reg_data, alr_to_i_phs)
      class(evt_nlo_t), intent(inout) :: evt
      type(region_data_t), intent(in) :: reg_data
      integer, intent(in), dimension(:) :: alr_to_i_phs
    end subroutine evt_nlo_set_i_evaluation_mappings
    module function evt_nlo_get_i_phs (evt) result (i_phs)
      integer :: i_phs
      class(evt_nlo_t), intent(in) :: evt
    end function evt_nlo_get_i_phs  
    module function evt_nlo_get_emitter (evt) result (emitter)
      integer :: emitter
      class(evt_nlo_t), intent(in) :: evt
    end function evt_nlo_get_emitter
    module function evt_nlo_get_i_term (evt) result (i_term)
      integer :: i_term
      class(evt_nlo_t), intent(in) :: evt
    end function evt_nlo_get_i_term
    module subroutine evt_nlo_generate_weighted (evt, probability)
      class(evt_nlo_t), intent(inout) :: evt
      real(default), intent(inout) :: probability
    end subroutine evt_nlo_generate_weighted
    module subroutine evt_nlo_reset_phs_identifiers (evt)
       class(evt_nlo_t), intent(inout) :: evt
    end subroutine evt_nlo_reset_phs_identifiers
    module subroutine evt_nlo_connected_set_real_IS_momenta (evt)
      class(evt_nlo_t), intent(inout) :: evt
    end subroutine evt_nlo_connected_set_real_IS_momenta
    module subroutine evt_nlo_make_particle_set &
         (evt, factorization_mode, keep_correlations, r)
      class(evt_nlo_t), intent(inout) :: evt
      integer, intent(in) :: factorization_mode
      logical, intent(in) :: keep_correlations
      real(default), dimension(:), intent(in), optional :: r
    end subroutine evt_nlo_make_particle_set
    module subroutine evt_nlo_evaluate_real_kinematics (evt)
      class(evt_nlo_t), intent(inout) :: evt
    end subroutine evt_nlo_evaluate_real_kinematics
    module function evt_nlo_compute_subtraction_sqmes (evt) result (sqme)
      class(evt_nlo_t), intent(inout) :: evt
      real(default) :: sqme
    end function evt_nlo_compute_subtraction_sqmes
    module subroutine evt_nlo_compute_real (evt)
      class(evt_nlo_t), intent(inout) :: evt
    end subroutine evt_nlo_compute_real
    module function evt_nlo_compute_all_sqme_rad (evt) result (sqme)
      class(evt_nlo_t), intent(inout) :: evt
      real(default) :: sqme
    end function evt_nlo_compute_all_sqme_rad
    module function evt_nlo_boost_to_cms (evt, p_lab) result (p_cms)
      type(phs_point_t), intent(in) :: p_lab
      class(evt_nlo_t), intent(in) :: evt
      type(phs_point_t) :: p_cms
    end function evt_nlo_boost_to_cms
    module function evt_nlo_boost_to_lab (evt, p_cms) result (p_lab)
      type(phs_point_t) :: p_lab
      class(evt_nlo_t), intent(in) :: evt
      type(phs_point_t), intent(in) :: p_cms
    end function evt_nlo_boost_to_lab
    module subroutine evt_nlo_setup_general_event_kinematics &
         (evt, process_instance)
      class(evt_nlo_t), intent(inout) :: evt
      type(process_instance_t), intent(in) :: process_instance
    end subroutine evt_nlo_setup_general_event_kinematics
    module subroutine evt_nlo_setup_real_event_kinematics &
         (evt, process_instance)
      class(evt_nlo_t), intent(inout) :: evt
      type(process_instance_t), intent(in) :: process_instance
    end subroutine evt_nlo_setup_real_event_kinematics
    module subroutine evt_nlo_set_mode (evt, process_instance)
      class(evt_nlo_t), intent(inout) :: evt
      type(process_instance_t), intent(in) :: process_instance
    end subroutine evt_nlo_set_mode
    module function evt_nlo_is_valid_event (evt, i_term) result (valid)
      logical :: valid
      class(evt_nlo_t), intent(in) :: evt
      integer, intent(in) :: i_term
    end function evt_nlo_is_valid_event
    module function evt_nlo_get_selected_quantum_numbers &
         (evt, i_flv) result (qn_select)
      class(evt_nlo_t), intent(in) :: evt
      integer, intent(in) :: i_flv
      type(quantum_numbers_t), dimension(:), allocatable :: qn_select
    end function evt_nlo_get_selected_quantum_numbers
    module subroutine evt_nlo_prepare_new_event (evt, i_mci, i_term)
      class(evt_nlo_t), intent(inout) :: evt
      integer, intent(in) :: i_mci, i_term
    end subroutine evt_nlo_prepare_new_event
  end interface

end module evt_nlo

