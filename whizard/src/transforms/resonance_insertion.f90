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

module resonance_insertion

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use rng_base, only: rng_t
  use selectors, only: selector_t
  use particles, only: particle_t, particle_set_t
  use resonances, only: resonance_history_set_t
  use resonances, only: resonance_tree_t
  use instances, only: process_instance_ptr_t
  use event_transforms

  implicit none
  private

  public :: evt_resonance_t

  type, extends (evt_t) :: evt_resonance_t
     type(resonance_history_set_t), dimension(:), allocatable :: res_history_set
     integer, dimension(:), allocatable :: index_offset
     integer :: selected_component = 0
     type(string_t) :: libname
     type(string_t), dimension(:), allocatable :: proc_id
     real(default) :: on_shell_limit = 0
     real(default) :: on_shell_turnoff = 0
     real(default) :: background_factor = 1
     logical :: selector_active = .false.
     type(selector_t) :: selector
     integer :: selected_history = 0
     type(process_instance_ptr_t), dimension(:), allocatable :: instance
   contains
     procedure :: write_name => evt_resonance_write_name
     procedure :: write => evt_resonance_write
     procedure :: set_resonance_data => evt_resonance_set_resonance_data
     procedure :: set_library => evt_resonance_set_library
     procedure :: set_subprocess_instances &
          => evt_resonance_set_subprocess_instances
     procedure :: set_on_shell_limit => evt_resonance_set_on_shell_limit
     procedure :: set_on_shell_turnoff => evt_resonance_set_on_shell_turnoff
     procedure :: set_background_factor => evt_resonance_set_background_factor
     procedure :: import_rng => evt_resonance_import_rng
     procedure :: write_selector => evt_resonance_write_selector
     procedure :: init_selector => evt_resonance_init_selector
     procedure :: get_selector_weights => evt_resonance_get_selector_weights
     procedure :: fill_momenta => evt_resonance_fill_momenta
     procedure :: determine_on_shell_histories &
          => evt_resonance_determine_on_shell_histories
     procedure :: evaluate_subprocess => evt_resonance_evaluate_subprocess
     procedure :: get_master_sqme => evt_resonance_get_master_sqme
     procedure :: get_subprocess_sqme => evt_resonance_get_subprocess_sqme
     procedure :: apply_turnoff_factor => evt_resonance_apply_turnoff_factor
     procedure :: compute_probabilities => evt_resonance_compute_probabilities
     procedure :: select_component => evt_resonance_select_component
     procedure :: find_prt_invalid_color => evt_resonance_find_prt_invalid_color
     procedure :: prepare_new_event => evt_resonance_prepare_new_event
     procedure :: generate_weighted => evt_resonance_generate_weighted
     procedure :: make_particle_set => evt_resonance_make_particle_set
  end type evt_resonance_t


  interface
    module subroutine evt_resonance_write_name (evt, unit)
      class(evt_resonance_t), intent(in) :: evt
      integer, intent(in), optional :: unit
    end subroutine evt_resonance_write_name
    module subroutine evt_resonance_write &
         (evt, unit, verbose, more_verbose, testflag)
      class(evt_resonance_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, more_verbose, testflag
    end subroutine evt_resonance_write
    module subroutine evt_resonance_set_resonance_data (evt, res_history_set)
      class(evt_resonance_t), intent(inout) :: evt
      type(resonance_history_set_t), dimension(:), intent(in) :: res_history_set
    end subroutine evt_resonance_set_resonance_data
    module subroutine evt_resonance_set_library (evt, libname)
      class(evt_resonance_t), intent(inout) :: evt
      type(string_t), intent(in) :: libname
    end subroutine evt_resonance_set_library
    module subroutine evt_resonance_set_subprocess_instances (evt, instance)
      class(evt_resonance_t), intent(inout) :: evt
      type(process_instance_ptr_t), dimension(:), intent(in) :: instance
    end subroutine evt_resonance_set_subprocess_instances
    module subroutine evt_resonance_set_on_shell_limit (evt, on_shell_limit)
      class(evt_resonance_t), intent(inout) :: evt
      real(default), intent(in) :: on_shell_limit
    end subroutine evt_resonance_set_on_shell_limit
    module subroutine evt_resonance_set_on_shell_turnoff (evt, on_shell_turnoff)
      class(evt_resonance_t), intent(inout) :: evt
      real(default), intent(in) :: on_shell_turnoff
    end subroutine evt_resonance_set_on_shell_turnoff
    module subroutine evt_resonance_set_background_factor &
         (evt, background_factor)
      class(evt_resonance_t), intent(inout) :: evt
      real(default), intent(in) :: background_factor
    end subroutine evt_resonance_set_background_factor
    module subroutine evt_resonance_import_rng (evt, rng)
      class(evt_resonance_t), intent(inout) :: evt
      class(rng_t), allocatable, intent(inout) :: rng
    end subroutine evt_resonance_import_rng
    module subroutine evt_resonance_write_selector (evt, unit, testflag)
      class(evt_resonance_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine evt_resonance_write_selector
    module subroutine evt_resonance_init_selector (evt, weight, offset)
      class(evt_resonance_t), intent(inout) :: evt
      real(default), dimension(:), intent(in) :: weight
      integer, intent(in), optional :: offset
    end subroutine evt_resonance_init_selector
    module subroutine evt_resonance_get_selector_weights (evt, weight)
      class(evt_resonance_t), intent(in) :: evt
      real(default), dimension(0:), intent(out) :: weight
    end subroutine evt_resonance_get_selector_weights
    module subroutine evt_resonance_fill_momenta (evt)
      class(evt_resonance_t), intent(inout) :: evt
    end subroutine evt_resonance_fill_momenta
    module subroutine evt_resonance_determine_on_shell_histories &
         (evt, index_array)
      class(evt_resonance_t), intent(in) :: evt
      integer, dimension(:), allocatable, intent(out) :: index_array
    end subroutine evt_resonance_determine_on_shell_histories
    module subroutine evt_resonance_evaluate_subprocess (evt, index_array)
      class(evt_resonance_t), intent(inout) :: evt
      integer, dimension(:), intent(in) :: index_array
    end subroutine evt_resonance_evaluate_subprocess
    module function evt_resonance_get_master_sqme (evt) result (sqme)
      class(evt_resonance_t), intent(in) :: evt
      real(default) :: sqme
    end function evt_resonance_get_master_sqme
    module subroutine evt_resonance_get_subprocess_sqme (evt, sqme, index_array)
      class(evt_resonance_t), intent(in) :: evt
      real(default), dimension(:), intent(out) :: sqme
      integer, dimension(:), intent(in), optional :: index_array
    end subroutine evt_resonance_get_subprocess_sqme
    module subroutine evt_resonance_apply_turnoff_factor &
         (evt, sqme, index_array)
      class(evt_resonance_t), intent(in) :: evt
      real(default), dimension(:), intent(inout) :: sqme
      integer, dimension(:), intent(in) :: index_array
    end subroutine evt_resonance_apply_turnoff_factor
    module subroutine evt_resonance_compute_probabilities (evt)
      class(evt_resonance_t), intent(inout) :: evt
    end subroutine evt_resonance_compute_probabilities
    module subroutine evt_resonance_select_component (evt, i_component)
      class(evt_resonance_t), intent(inout) :: evt
      integer, intent(in) :: i_component
    end subroutine evt_resonance_select_component
    module subroutine evt_resonance_find_prt_invalid_color (evt, index, prt)
      class(evt_resonance_t), intent(in) :: evt
      integer, dimension(:), allocatable, intent(out) :: index
      type(particle_t), dimension(:), allocatable, intent(out), optional :: prt
    end subroutine evt_resonance_find_prt_invalid_color
    module subroutine evt_resonance_prepare_new_event (evt, i_mci, i_term)
      class(evt_resonance_t), intent(inout) :: evt
      integer, intent(in) :: i_mci, i_term
    end subroutine evt_resonance_prepare_new_event
    module subroutine evt_resonance_generate_weighted (evt, probability)
      class(evt_resonance_t), intent(inout) :: evt
      real(default), intent(inout) :: probability
    end subroutine evt_resonance_generate_weighted
    module subroutine evt_resonance_make_particle_set &
         (evt, factorization_mode, keep_correlations, r)
      class(evt_resonance_t), intent(inout) :: evt
      integer, intent(in) :: factorization_mode
      logical, intent(in) :: keep_correlations
      real(default), dimension(:), intent(in), optional :: r
    end subroutine evt_resonance_make_particle_set
  end interface

end module resonance_insertion
