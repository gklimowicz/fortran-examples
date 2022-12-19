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

module mci_vamp

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use diagnostics
  use phs_base
  use rng_base
  use rng_tao
  use vamp !NODEP!
  use exceptions !NODEP!

  use mci_base

  implicit none
  private

  public :: grid_parameters_t
  public :: history_parameters_t
  public :: mci_vamp_t
  public :: mci_vamp_instance_t

  type :: grid_parameters_t
     integer :: threshold_calls = 0
     integer :: min_calls_per_channel = 10
     integer :: min_calls_per_bin = 10
     integer :: min_bins = 3
     integer :: max_bins = 20
     logical :: stratified = .true.
     logical :: use_vamp_equivalences = .true.
     real(default) :: channel_weights_power = 0.25_default
     real(default) :: accuracy_goal = 0
     real(default) :: error_goal = 0
     real(default) :: rel_error_goal = 0
   contains
     procedure :: write => grid_parameters_write
  end type grid_parameters_t

  type :: history_parameters_t
     logical :: global = .true.
     logical :: global_verbose = .false.
     logical :: channel = .false.
     logical :: channel_verbose = .false.
   contains
     procedure :: write => history_parameters_write
  end type history_parameters_t

  type :: pass_t
     integer :: i_pass = 0
     integer :: i_first_it = 0
     integer :: n_it = 0
     integer :: n_calls = 0
     integer :: n_bins = 0
     logical :: adapt_grids = .false.
     logical :: adapt_weights = .false.
     logical :: is_final_pass = .false.
     logical :: integral_defined = .false.
     integer, dimension(:), allocatable :: calls
     integer, dimension(:), allocatable :: calls_valid
     real(default), dimension(:), allocatable :: integral
     real(default), dimension(:), allocatable :: error
     real(default), dimension(:), allocatable :: efficiency
     type(vamp_history), dimension(:), allocatable :: v_history
     type(vamp_history), dimension(:,:), allocatable :: v_histories
     type(pass_t), pointer :: next => null ()
   contains
     procedure :: final => pass_final
     procedure :: write => pass_write
     procedure :: read => pass_read
     procedure :: write_history => pass_write_history
     procedure :: configure => pass_configure
     procedure :: configure_history => pass_configure_history
     procedure :: update => pass_update
     procedure :: get_integration_index => pass_get_integration_index
     procedure :: get_calls => pass_get_calls
     procedure :: get_calls_valid => pass_get_calls_valid
     procedure :: get_integral => pass_get_integral
     procedure :: get_error => pass_get_error
     procedure :: get_efficiency => pass_get_efficiency
  end type pass_t

  type, extends (mci_t) :: mci_vamp_t
     logical, dimension(:), allocatable :: dim_is_flat
     type(grid_parameters_t) :: grid_par
     type(history_parameters_t) :: history_par
     integer :: min_calls = 0
     type(pass_t), pointer :: first_pass => null ()
     type(pass_t), pointer :: current_pass => null ()
     type(vamp_equivalences_t) :: equivalences
     logical :: rebuild = .true.
     logical :: check_grid_file = .true.
     logical :: grid_filename_set = .false.
     logical :: negative_weights = .false.
     logical :: verbose = .false.
     type(string_t) :: grid_filename
     character(32) :: md5sum_adapted = ""
   contains
     procedure :: reset => mci_vamp_reset
     procedure :: final => mci_vamp_final
     procedure :: write => mci_vamp_write
     procedure :: write_history_parameters => mci_vamp_write_history_parameters
     procedure :: write_history => mci_vamp_write_history
     procedure :: compute_md5sum => mci_vamp_compute_md5sum
     procedure :: get_md5sum => mci_vamp_get_md5sum
     procedure :: startup_message => mci_vamp_startup_message
     procedure :: write_log_entry => mci_vamp_write_log_entry
     procedure :: record_index => mci_vamp_record_index
     procedure :: set_grid_parameters => mci_vamp_set_grid_parameters
     procedure :: set_history_parameters => mci_vamp_set_history_parameters
     procedure :: set_rebuild_flag => mci_vamp_set_rebuild_flag
     procedure :: set_grid_filename => mci_vamp_set_grid_filename
     procedure :: prepend_grid_path => mci_vamp_prepend_grid_path
     procedure :: declare_flat_dimensions => mci_vamp_declare_flat_dimensions
     procedure :: declare_equivalences => mci_vamp_declare_equivalences
     procedure :: allocate_instance => mci_vamp_allocate_instance
     procedure :: add_pass => mci_vamp_add_pass
     procedure :: update_from_ref => mci_vamp_update_from_ref
     procedure :: update => mci_vamp_update
     procedure :: write_grids => mci_vamp_write_grids
     procedure :: read_grids_header => mci_vamp_read_grids_header
     procedure :: read_grids_data => mci_vamp_read_grids_data
     procedure :: read_grids => mci_vamp_read_grids
     procedure :: integrate => mci_vamp_integrate
     procedure :: check_goals => mci_vamp_check_goals
     procedure :: error_reached => mci_vamp_error_reached
     procedure :: rel_error_reached => mci_vamp_rel_error_reached
     procedure :: accuracy_reached => mci_vamp_accuracy_reached
     procedure :: prepare_simulation => mci_vamp_prepare_simulation
     procedure :: generate_weighted_event => mci_vamp_generate_weighted_event
     procedure :: generate_unweighted_event => &
          mci_vamp_generate_unweighted_event
     procedure :: rebuild_event => mci_vamp_rebuild_event
     procedure :: pacify => mci_vamp_pacify
  end type mci_vamp_t

  type, extends (vamp_data_t) :: mci_workspace_t
     class(mci_sampler_t), pointer :: sampler => null ()
     class(mci_vamp_instance_t), pointer :: instance => null ()
  end type mci_workspace_t

  type, extends (mci_instance_t) :: mci_vamp_instance_t
     type(mci_vamp_t), pointer :: mci => null ()
     logical :: grids_defined = .false.
     logical :: grids_from_file = .false.
     integer :: n_it = 0
     integer :: it = 0
     logical :: pass_complete = .false.
     integer :: n_calls = 0
     integer :: calls = 0
     integer :: calls_valid = 0
     logical :: it_complete = .false.
     logical :: enable_adapt_grids = .false.
     logical :: enable_adapt_weights = .false.
     logical :: allow_adapt_grids = .false.
     logical :: allow_adapt_weights = .false.
     integer :: n_adapt_grids = 0
     integer :: n_adapt_weights = 0
     logical :: generating_events = .false.
     real(default) :: safety_factor = 1
     type(vamp_grids) :: grids
     real(default) :: g = 0
     real(default), dimension(:), allocatable :: gi
     real(default) :: integral = 0
     real(default) :: error = 0
     real(default) :: efficiency = 0
     real(default), dimension(:), allocatable :: vamp_x
     logical :: vamp_weight_set = .false.
     real(default) :: vamp_weight = 0
     real(default) :: vamp_excess = 0
     logical :: allocate_global_history = .false.
     type(vamp_history), dimension(:), pointer :: v_history => null ()
     logical :: allocate_channel_history = .false.
     type(vamp_history), dimension(:,:), pointer :: v_histories => null ()
   contains
     procedure :: write => mci_vamp_instance_write
     procedure :: write_grids => mci_vamp_instance_write_grids
     procedure :: final => mci_vamp_instance_final
     procedure :: init => mci_vamp_instance_init
     procedure :: new_pass => mci_vamp_instance_new_pass
     procedure :: create_grids => mci_vamp_instance_create_grids
     procedure :: discard_integrals => mci_vamp_instance_discard_integrals
     procedure :: allow_adaptation => mci_vamp_instance_allow_adaptation
     procedure :: adapt_grids => mci_vamp_instance_adapt_grids
     procedure :: adapt_weights => mci_vamp_instance_adapt_weights
     procedure :: sample_grids => mci_vamp_instance_sample_grids
     procedure :: get_efficiency_array => mci_vamp_instance_get_efficiency_array
     procedure :: get_efficiency => mci_vamp_instance_get_efficiency
     procedure :: init_simulation => mci_vamp_instance_init_simulation
     procedure :: final_simulation => mci_vamp_instance_final_simulation
     procedure :: compute_weight => mci_vamp_instance_compute_weight
     procedure :: record_integrand => mci_vamp_instance_record_integrand
     procedure :: get_event_weight => mci_vamp_instance_get_event_weight
     procedure :: get_event_excess => mci_vamp_instance_get_event_excess
  end type mci_vamp_instance_t


  interface operator (.matches.)
     module procedure pass_matches
  end interface operator (.matches.)
  interface operator (.matches.)
     module procedure real_matches
  end interface operator (.matches.)

  interface
    module subroutine grid_parameters_write (object, unit)
      class(grid_parameters_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine grid_parameters_write
    module subroutine history_parameters_write (object, unit)
      class(history_parameters_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine history_parameters_write
    module subroutine pass_final (object)
      class(pass_t), intent(inout) :: object
    end subroutine pass_final
    module subroutine pass_write (object, unit, pacify)
      class(pass_t), intent(in) :: object
      integer, intent(in) :: unit
      logical, intent(in), optional :: pacify
    end subroutine pass_write
    module subroutine pass_read (object, u, n_pass, n_it)
      class(pass_t), intent(out) :: object
      integer, intent(in) :: u, n_pass, n_it
    end subroutine pass_read
    module subroutine pass_write_history (pass, unit)
      class(pass_t), intent(in) :: pass
      integer, intent(in), optional :: unit
    end subroutine pass_write_history
    module subroutine pass_configure (pass, n_it, n_calls, min_calls, &
         min_bins, max_bins, min_channel_calls)
      class(pass_t), intent(inout) :: pass
      integer, intent(in) :: n_it, n_calls, min_channel_calls
      integer, intent(in) :: min_calls, min_bins, max_bins
    end subroutine pass_configure
    module subroutine pass_configure_history (pass, n_channels, par)
      class(pass_t), intent(inout) :: pass
      integer, intent(in) :: n_channels
      type(history_parameters_t), intent(in) :: par
    end subroutine pass_configure_history
    module function pass_matches (pass, ref) result (ok)
      type(pass_t), intent(in) :: pass, ref
      integer :: n
      logical :: ok
    end function pass_matches
    module subroutine pass_update (pass, ref, ok)
      class(pass_t), intent(inout) :: pass
      type(pass_t), intent(in) :: ref
      logical, intent(out) :: ok
    end subroutine pass_update
    elemental module function real_matches (x, y) result (ok)
      real(default), intent(in) :: x, y
      logical :: ok
    end function real_matches
    module function pass_get_integration_index (pass) result (n)
      class (pass_t), intent(in) :: pass
      integer :: n
    end function pass_get_integration_index
    module function pass_get_calls (pass) result (calls)
      class(pass_t), intent(in) :: pass
      integer :: calls
    end function pass_get_calls
    module function pass_get_calls_valid (pass) result (calls_valid)
      class(pass_t), intent(in) :: pass
      integer :: calls_valid
    end function pass_get_calls_valid
    module function pass_get_integral (pass) result (integral)
      class(pass_t), intent(in) :: pass
      real(default) :: integral
    end function pass_get_integral
    module function pass_get_error (pass) result (error)
      class(pass_t), intent(in) :: pass
      real(default) :: error
    end function pass_get_error
    module function pass_get_efficiency (pass) result (efficiency)
      class(pass_t), intent(in) :: pass
      real(default) :: efficiency
    end function pass_get_efficiency
    module subroutine mci_vamp_reset (object)
      class(mci_vamp_t), intent(inout) :: object
    end subroutine mci_vamp_reset
    module subroutine mci_vamp_final (object)
      class(mci_vamp_t), intent(inout) :: object
    end subroutine mci_vamp_final
    module subroutine mci_vamp_write (object, unit, pacify, md5sum_version)
      class(mci_vamp_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacify
      logical, intent(in), optional :: md5sum_version
    end subroutine mci_vamp_write
    module subroutine mci_vamp_write_history_parameters (mci, unit)
      class(mci_vamp_t), intent(in) :: mci
      integer, intent(in), optional :: unit
    end subroutine mci_vamp_write_history_parameters
    module subroutine mci_vamp_write_history (mci, unit)
      class(mci_vamp_t), intent(in) :: mci
      integer, intent(in), optional :: unit
    end subroutine mci_vamp_write_history
    module subroutine mci_vamp_compute_md5sum (mci, pacify)
      class(mci_vamp_t), intent(inout) :: mci
      logical, intent(in), optional :: pacify
    end subroutine mci_vamp_compute_md5sum
    pure module function mci_vamp_get_md5sum (mci) result (md5sum)
      class(mci_vamp_t), intent(in) :: mci
      character(32) :: md5sum
    end function mci_vamp_get_md5sum
    module subroutine mci_vamp_startup_message (mci, unit, n_calls)
      class(mci_vamp_t), intent(in) :: mci
      integer, intent(in), optional :: unit, n_calls
    end subroutine mci_vamp_startup_message
    module subroutine mci_vamp_write_log_entry (mci, u)
      class(mci_vamp_t), intent(in) :: mci
      integer, intent(in) :: u
    end subroutine mci_vamp_write_log_entry
    module subroutine mci_vamp_record_index (mci, i_mci)
      class(mci_vamp_t), intent(inout) :: mci
      integer, intent(in) :: i_mci
    end subroutine mci_vamp_record_index
    module subroutine mci_vamp_set_grid_parameters (mci, grid_par)
      class(mci_vamp_t), intent(inout) :: mci
      type(grid_parameters_t), intent(in) :: grid_par
    end subroutine mci_vamp_set_grid_parameters
    module subroutine mci_vamp_set_history_parameters (mci, history_par)
      class(mci_vamp_t), intent(inout) :: mci
      type(history_parameters_t), intent(in) :: history_par
    end subroutine mci_vamp_set_history_parameters
    module subroutine mci_vamp_set_rebuild_flag (mci, rebuild, check_grid_file)
      class(mci_vamp_t), intent(inout) :: mci
      logical, intent(in) :: rebuild
      logical, intent(in) :: check_grid_file
    end subroutine mci_vamp_set_rebuild_flag
    module subroutine mci_vamp_set_grid_filename (mci, name, run_id)
      class(mci_vamp_t), intent(inout) :: mci
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: run_id
    end subroutine mci_vamp_set_grid_filename
    module subroutine mci_vamp_prepend_grid_path (mci, prefix)
      class(mci_vamp_t), intent(inout) :: mci
      type(string_t), intent(in) :: prefix
    end subroutine mci_vamp_prepend_grid_path
    module subroutine mci_vamp_declare_flat_dimensions (mci, dim_flat)
      class(mci_vamp_t), intent(inout) :: mci
      integer, dimension(:), intent(in) :: dim_flat
    end subroutine mci_vamp_declare_flat_dimensions
    module subroutine mci_vamp_declare_equivalences (mci, channel, dim_offset)
      class(mci_vamp_t), intent(inout) :: mci
      type(phs_channel_t), dimension(:), intent(in) :: channel
      integer, intent(in) :: dim_offset
    end subroutine mci_vamp_declare_equivalences
    module subroutine mci_vamp_add_pass &
         (mci, adapt_grids, adapt_weights, final_pass)
      class(mci_vamp_t), intent(inout) :: mci
      logical, intent(in), optional :: adapt_grids, adapt_weights, final_pass
    end subroutine mci_vamp_add_pass
    module subroutine mci_vamp_update_from_ref (mci, mci_ref, success)
      class(mci_vamp_t), intent(inout) :: mci
      class(mci_t), intent(in) :: mci_ref
      logical, intent(out) :: success
    end subroutine mci_vamp_update_from_ref
    module subroutine mci_vamp_update (mci, u, success)
      class(mci_vamp_t), intent(inout) :: mci
      integer, intent(in) :: u
      logical, intent(out) :: success
    end subroutine mci_vamp_update
    module subroutine mci_vamp_write_grids (mci, instance)
      class(mci_vamp_t), intent(in) :: mci
      class(mci_instance_t), intent(inout) :: instance
    end subroutine mci_vamp_write_grids
    module subroutine mci_vamp_read_grids_header (mci, success)
      class(mci_vamp_t), intent(inout) :: mci
      logical, intent(out) :: success
    end subroutine mci_vamp_read_grids_header
    module subroutine mci_vamp_read_grids_data (mci, instance, read_integrals)
      class(mci_vamp_t), intent(in) :: mci
      class(mci_instance_t), intent(inout) :: instance
      logical, intent(in), optional :: read_integrals
    end subroutine mci_vamp_read_grids_data
    module subroutine mci_vamp_read_grids (mci, instance, success)
      class(mci_vamp_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout) :: instance
      logical, intent(out) :: success
    end subroutine mci_vamp_read_grids
    module subroutine mci_vamp_integrate (mci, instance, sampler, &
         n_it, n_calls, results, pacify)
      class(mci_vamp_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout), target :: instance
      class(mci_sampler_t), intent(inout), target :: sampler
      integer, intent(in) :: n_it
      integer, intent(in) :: n_calls
      class(mci_results_t), intent(inout), optional :: results
      logical, intent(in), optional :: pacify
    end subroutine mci_vamp_integrate
    module subroutine mci_vamp_check_goals (mci, it, success)
      class(mci_vamp_t), intent(inout) :: mci
      integer, intent(in) :: it
      logical, intent(out) :: success
    end subroutine mci_vamp_check_goals
    module function mci_vamp_error_reached (mci, it) result (flag)
      class(mci_vamp_t), intent(in) :: mci
      integer, intent(in) :: it
      logical :: flag
    end function mci_vamp_error_reached
    module function mci_vamp_rel_error_reached (mci, it) result (flag)
      class(mci_vamp_t), intent(in) :: mci
      integer, intent(in) :: it
      logical :: flag
    end function mci_vamp_rel_error_reached
    module function mci_vamp_accuracy_reached (mci, it) result (flag)
      class(mci_vamp_t), intent(in) :: mci
      integer, intent(in) :: it
      logical :: flag
    end function mci_vamp_accuracy_reached
    module subroutine mci_vamp_prepare_simulation (mci)
      class(mci_vamp_t), intent(inout) :: mci
    end subroutine mci_vamp_prepare_simulation
    module subroutine mci_vamp_rebuild_event (mci, instance, sampler, state)
      class(mci_vamp_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout) :: instance
      class(mci_sampler_t), intent(inout) :: sampler
      class(mci_state_t), intent(in) :: state
    end subroutine mci_vamp_rebuild_event
    module subroutine mci_vamp_pacify (object, efficiency_reset, error_reset)
      class(mci_vamp_t), intent(inout) :: object
      logical, intent(in), optional :: efficiency_reset, error_reset
    end subroutine mci_vamp_pacify
    module subroutine mci_vamp_instance_write (object, unit, pacify)
      class(mci_vamp_instance_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacify
    end subroutine mci_vamp_instance_write
    module subroutine mci_vamp_instance_write_grids (object, unit)
      class(mci_vamp_instance_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine mci_vamp_instance_write_grids
    module subroutine mci_vamp_instance_final (object)
      class(mci_vamp_instance_t), intent(inout) :: object
    end subroutine mci_vamp_instance_final
    module subroutine mci_vamp_instance_init (mci_instance, mci)
      class(mci_vamp_instance_t), intent(out) :: mci_instance
      class(mci_t), intent(in), target :: mci
    end subroutine mci_vamp_instance_init
    module subroutine mci_vamp_instance_new_pass (instance, reshape)
      class(mci_vamp_instance_t), intent(inout) :: instance
      logical, intent(out) :: reshape
    end subroutine mci_vamp_instance_new_pass
    module subroutine mci_vamp_instance_create_grids (instance)
      class(mci_vamp_instance_t), intent(inout) :: instance
    end subroutine mci_vamp_instance_create_grids
    module subroutine mci_vamp_instance_discard_integrals (instance, reshape)
      class(mci_vamp_instance_t), intent(inout) :: instance
      logical, intent(in) :: reshape
    end subroutine mci_vamp_instance_discard_integrals
    module subroutine mci_vamp_instance_allow_adaptation (instance)
      class(mci_vamp_instance_t), intent(inout) :: instance
    end subroutine mci_vamp_instance_allow_adaptation
    module subroutine mci_vamp_instance_adapt_grids (instance)
      class(mci_vamp_instance_t), intent(inout) :: instance
    end subroutine mci_vamp_instance_adapt_grids
    module function mci_vamp_instance_get_efficiency_array &
         (mci) result (efficiency)
      class(mci_vamp_instance_t), intent(in) :: mci
      real(default), dimension(:), allocatable :: efficiency
    end function mci_vamp_instance_get_efficiency_array
    module function mci_vamp_instance_get_efficiency (mci) result (efficiency)
      class(mci_vamp_instance_t), intent(in) :: mci
      real(default) :: efficiency
    end function mci_vamp_instance_get_efficiency
    module subroutine mci_vamp_instance_init_simulation &
         (instance, safety_factor)
      class(mci_vamp_instance_t), intent(inout) :: instance
      real(default), intent(in), optional :: safety_factor
    end subroutine mci_vamp_instance_init_simulation
    module subroutine mci_vamp_instance_final_simulation (instance)
      class(mci_vamp_instance_t), intent(inout) :: instance
    end subroutine mci_vamp_instance_final_simulation
    module subroutine mci_vamp_instance_compute_weight (mci, c)
      class(mci_vamp_instance_t), intent(inout) :: mci
      integer, intent(in) :: c
    end subroutine mci_vamp_instance_compute_weight
    module subroutine mci_vamp_instance_record_integrand (mci, integrand)
      class(mci_vamp_instance_t), intent(inout) :: mci
      real(default), intent(in) :: integrand
    end subroutine mci_vamp_instance_record_integrand
    module function mci_vamp_instance_get_event_weight (mci) result (value)
      class(mci_vamp_instance_t), intent(in) :: mci
      real(default) :: value
    end function mci_vamp_instance_get_event_weight
    module function mci_vamp_instance_get_event_excess (mci) result (value)
      class(mci_vamp_instance_t), intent(in) :: mci
      real(default) :: value
    end function mci_vamp_instance_get_event_excess
  end interface

contains

  subroutine mci_vamp_allocate_instance (mci, mci_instance)
    class(mci_vamp_t), intent(in) :: mci
    class(mci_instance_t), intent(out), pointer :: mci_instance
    allocate (mci_vamp_instance_t :: mci_instance)
  end subroutine mci_vamp_allocate_instance

  function vamp_sampling_function &
       (xi, data, weights, channel, grids) result (f)
    real(default) :: f
    real(default), dimension(:), intent(in) :: xi
    class(vamp_data_t), intent(in) :: data
    real(default), dimension(:), intent(in), optional :: weights
    integer, intent(in), optional :: channel
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    type(exception) :: exc
    logical :: verbose
    character(*), parameter :: FN = "WHIZARD sampling function"
    class(mci_instance_t), pointer :: instance
    select type (data)
    type is (mci_workspace_t)
       instance => data%instance
       select type (instance)
       class is (mci_vamp_instance_t)
          verbose = instance%mci%verbose
          call instance%evaluate (data%sampler, channel, xi)
          if (signal_is_pending ()) then
             call raise_exception (exc, EXC_FATAL, FN, "signal received")
             call handle_vamp_exception (exc, verbose)
             call terminate_now_if_signal ()
          end if
          instance%calls = instance%calls + 1
          if (data%sampler%is_valid ()) &
               & instance%calls_valid = instance%calls_valid + 1
          f = instance%get_value ()
          call terminate_now_if_single_event ()
       class default
          call msg_bug("VAMP: " // FN // ": unknown MCI instance type")
       end select
    end select
  end function vamp_sampling_function

  pure function phi_trivial (xi, channel_dummy) result (x)
    real(default), dimension(:), intent(in) :: xi
    integer, intent(in) :: channel_dummy
    real(default), dimension(size(xi)) :: x
    x = xi
  end function phi_trivial

  subroutine mci_vamp_generate_weighted_event (mci, instance, sampler)
    class(mci_vamp_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout), target :: instance
    class(mci_sampler_t), intent(inout), target :: sampler
    class(vamp_data_t), allocatable :: data
    type(exception) :: vamp_exception
    select type (instance)
    type is (mci_vamp_instance_t)
       instance%vamp_weight_set = .false.
       allocate (mci_workspace_t :: data)
       select type (data)
       type is (mci_workspace_t)
          data%sampler => sampler
          data%instance => instance
       end select
       select type (rng => mci%rng)
       type is (rng_tao_t)
          if (instance%grids_defined) then
             call vamp_next_event ( &
                  instance%vamp_x, &
                  rng%state, &
                  instance%grids, &
                  vamp_sampling_function, &
                  data, &
                  phi = phi_trivial, &
                  weight = instance%vamp_weight, &
                  exc = vamp_exception)
             call handle_vamp_exception (vamp_exception, mci%verbose)
             instance%vamp_excess = 0
             instance%vamp_weight_set = .true.
          else
             call msg_bug ("VAMP: generate event: grids undefined")
          end if
       class default
         call msg_fatal ("VAMP event generation: &
               &random-number generator must be TAO")
       end select
    end select
  end subroutine mci_vamp_generate_weighted_event

  subroutine mci_vamp_generate_unweighted_event &
       (mci, instance, sampler)
    class(mci_vamp_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout), target :: instance
    class(mci_sampler_t), intent(inout), target :: sampler
    class(vamp_data_t), allocatable :: data
    logical :: positive
    type(exception) :: vamp_exception
    select type (instance)
    type is (mci_vamp_instance_t)
       instance%vamp_weight_set = .false.
       allocate (mci_workspace_t :: data)
       select type (data)
       type is (mci_workspace_t)
          data%sampler => sampler
          data%instance => instance
       end select
       select type (rng => mci%rng)
       type is (rng_tao_t)
          if (instance%grids_defined) then
             REJECTION: do
                call vamp_next_event ( &
                     instance%vamp_x, &
                     rng%state, &
                     instance%grids, &
                     vamp_sampling_function, &
                     data, &
                     phi = phi_trivial, &
                     excess = instance%vamp_excess, &
                     positive = positive, &
                     exc = vamp_exception)
                if (signal_is_pending ())  return
                if (sampler%is_valid ())  exit REJECTION
             end do REJECTION
             call handle_vamp_exception (vamp_exception, mci%verbose)
             if (positive) then
                instance%vamp_weight = 1
             else if (instance%negative_weights) then
                instance%vamp_weight = -1
             else
                call msg_fatal ("VAMP: event with negative weight generated")
                instance%vamp_weight = 0
             end if
             instance%vamp_weight_set = .true.
          else
             call msg_bug ("VAMP: generate event: grids undefined")
          end if
       class default
         call msg_fatal ("VAMP event generation: &
               &random-number generator must be TAO")
       end select
    end select
  end subroutine mci_vamp_generate_unweighted_event

  subroutine mci_vamp_instance_adapt_weights (instance)
    class(mci_vamp_instance_t), intent(inout) :: instance
    real(default) :: w_sum, w_avg_ch, sum_w_underflow, w_min
    real(default), dimension(:), allocatable :: weights
    integer :: n_ch, ch, n_underflow
    logical, dimension(:), allocatable :: mask, underflow
    type(exception) :: vamp_exception
    logical :: wsum_non_zero
    if (instance%enable_adapt_weights .and. instance%allow_adapt_weights) then
       associate (mci => instance%mci)
         if (instance%grids_defined) then
            allocate (weights (size (instance%grids%weights)))
            weights = instance%grids%weights &
                 * vamp_get_variance (instance%grids%grids) &
                 ** mci%grid_par%channel_weights_power
            w_sum = sum (weights)
            if (w_sum /= 0) then
               weights = weights / w_sum
               if (mci%n_chain /= 0) then
                  allocate (mask (mci%n_channel))
                  do ch = 1, mci%n_chain
                     mask = mci%chain == ch
                     n_ch = count (mask)
                     if (n_ch /= 0) then
                        w_avg_ch = sum (weights, mask) / n_ch
                        where (mask)  weights = w_avg_ch
                     end if
                  end do
               end if
               if (mci%grid_par%threshold_calls /= 0) then
                  w_min = &
                       real (mci%grid_par%threshold_calls, default) &
                       / instance%n_calls
                  allocate (underflow (mci%n_channel))
                  underflow = weights /= 0 .and. abs (weights) < w_min
                  n_underflow = count (underflow)
                  sum_w_underflow = sum (weights, mask=underflow)
                  if (sum_w_underflow /= 1) then
                     where (underflow)
                        weights = w_min
                     elsewhere
                        weights = weights &
                             * (1 - n_underflow * w_min) / (1 - sum_w_underflow)
                     end where
                  end if
               end if
            end if
            call instance%set_channel_weights (weights, wsum_non_zero)
            if (wsum_non_zero) call vamp_update_weights &
               (instance%grids, weights, exc = vamp_exception)
            call handle_vamp_exception (vamp_exception, mci%verbose)
         else
            call msg_bug ("VAMP: adapt weights: grids undefined")
         end if
       end associate
       instance%n_adapt_weights = instance%n_adapt_weights + 1
    end if
  end subroutine mci_vamp_instance_adapt_weights

  subroutine mci_vamp_instance_sample_grids &
       (instance, rng, sampler, eq)
    class(mci_vamp_instance_t), intent(inout), target :: instance
    class(rng_t), intent(inout) :: rng
    class(mci_sampler_t), intent(inout), target :: sampler
    type(vamp_equivalences_t), intent(in), optional :: eq
    class(vamp_data_t), allocatable :: data
    type(exception) :: vamp_exception
    allocate (mci_workspace_t :: data)
    select type (data)
    type is (mci_workspace_t)
       data%sampler => sampler
       data%instance => instance
    end select
    select type (rng)
    type is (rng_tao_t)
       instance%it = instance%it + 1
       instance%calls = 0
       if (instance%grids_defined) then
          call vamp_sample_grids ( &
               rng%state, &
               instance%grids, &
               vamp_sampling_function, &
               data, &
               1, &
               eq = eq, &
               history = instance%v_history(instance%it:), &
               histories = instance%v_histories(instance%it:,:), &
               integral = instance%integral, &
               std_dev = instance%error, &
               exc = vamp_exception, &
               negative_weights = instance%negative_weights)
          call handle_vamp_exception (vamp_exception, instance%mci%verbose)
          instance%efficiency = instance%get_efficiency ()
       else
          call msg_bug ("VAMP: sample grids: grids undefined")
       end if
    class default
       call msg_fatal ("VAMP integration: random-number generator must be TAO")
    end select
  end subroutine mci_vamp_instance_sample_grids

  subroutine handle_vamp_exception (exc, verbose)
    type(exception), intent(in) :: exc
    logical, intent(in) :: verbose
    integer :: exc_level
    if (verbose) then
       exc_level = EXC_INFO
    else
       exc_level = EXC_ERROR
    end if
    if (exc%level >= exc_level) then
       write (msg_buffer, "(A,':',1x,A)")  trim (exc%origin), trim (exc%message)
       select case (exc%level)
       case (EXC_INFO);  call msg_message ()
       case (EXC_WARN);  call msg_warning ()
       case (EXC_ERROR); call msg_error ()
       case (EXC_FATAL)
          if (signal_is_pending ()) then
             call msg_message ()
          else
             call msg_fatal ()
          end if
       end select
    end if
  end subroutine handle_vamp_exception


end module mci_vamp
