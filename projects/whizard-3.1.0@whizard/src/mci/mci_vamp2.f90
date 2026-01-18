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

module mci_vamp2

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use phs_base
  use rng_base
  use mci_base
  use vamp2

  implicit none
  private

  public :: mci_vamp2_config_t
  public :: mci_vamp2_t
  public :: mci_vamp2_instance_t

  type, extends (vamp2_func_t) :: mci_vamp2_func_t
     private
     real(default) :: integrand = 0._default
     class(mci_sampler_t), pointer :: sampler => null ()
     class(mci_vamp2_instance_t), pointer :: instance => null ()
   contains
     procedure, public :: set_workspace => mci_vamp2_func_set_workspace
     procedure, public :: get_probabilities => mci_vamp2_func_get_probabilities
     procedure, public :: get_weight => mci_vamp2_func_get_weight
     procedure, public :: set_integrand => mci_vamp2_func_set_integrand
     procedure, public :: evaluate_maps => mci_vamp2_func_evaluate_maps
     procedure, public :: evaluate_func => mci_vamp2_func_evaluate_func
  end type mci_vamp2_func_t

  type, extends (vamp2_config_t) :: mci_vamp2_config_t
     !
  end type mci_vamp2_config_t

  type :: list_pass_t
     type(pass_t), pointer :: first => null ()
     type(pass_t), pointer :: current => null ()
   contains
     procedure :: final => list_pass_final
     procedure :: add => list_pass_add
     procedure :: update_from_ref => list_pass_update_from_ref
     procedure :: has_last_integral => list_pass_has_last_integral
     procedure :: get_last_integral => list_pass_get_last_integral
     procedure :: write => list_pass_write
  end type list_pass_t

  type :: pass_t
     integer :: i_pass = 0
     integer :: i_first_it = 0
     integer :: n_it = 0
     integer :: n_calls = 0
     logical :: adapt_grids = .false.
     logical :: adapt_weights = .false.
     logical :: is_final_pass = .false.
     logical :: integral_defined = .false.
     integer, dimension(:), allocatable :: calls
     integer, dimension(:), allocatable :: calls_valid
     real(default), dimension(:), allocatable :: integral
     real(default), dimension(:), allocatable :: error
     real(default), dimension(:), allocatable :: efficiency
     type(pass_t), pointer :: next => null ()
   contains
     procedure :: write => pass_write
     procedure :: read => pass_read
     procedure :: configure => pass_configure
     procedure :: update => pass_update
     procedure :: get_integration_index => pass_get_integration_index
     procedure :: get_calls => pass_get_calls
     procedure :: get_calls_valid => pass_get_calls_valid
     procedure :: get_integral => pass_get_integral
     procedure :: get_error => pass_get_error
     procedure :: get_efficiency => pass_get_efficiency
  end type pass_t

  type, extends(mci_t) :: mci_vamp2_t
     type(mci_vamp2_config_t) :: config
     type(vamp2_t) :: integrator
     type(vamp2_equivalences_t) :: equivalences
     logical :: integrator_defined = .false.
     logical :: integrator_from_file = .false.
     logical :: adapt_grids = .false.
     logical :: adapt_weights = .false.
     integer :: n_adapt_grids = 0
     integer :: n_adapt_weights = 0
     integer :: n_calls = 0
     type(list_pass_t) :: list_pass
     logical :: rebuild = .true.
     logical :: check_grid_file = .true.
     logical :: grid_filename_set = .false.
     logical :: negative_weights = .false.
     logical :: verbose = .false.
     logical :: pass_complete = .false.
     logical :: it_complete = .false.
     type(string_t) :: grid_filename
     integer :: grid_checkpoint = 1
     logical :: binary_grid_format = .false.
     type(string_t) :: parallel_method
     character(32) :: md5sum_adapted = ""
   contains
     procedure, public :: final => mci_vamp2_final
     procedure, public :: write => mci_vamp2_write
     procedure, public :: compute_md5sum => mci_vamp2_compute_md5sum
     procedure, public :: get_md5sum => mci_vamp2_get_md5sum
     procedure, public :: startup_message => mci_vamp2_startup_message
     procedure, public :: write_log_entry => mci_vamp2_write_log_entry
     procedure, public :: record_index => mci_vamp2_record_index
     procedure, public :: set_config => mci_vamp2_set_config
     procedure, public :: set_rebuild_flag => mci_vamp2_set_rebuild_flag
     procedure, public :: set_grid_filename => mci_vamp2_set_grid_filename
     procedure, public :: get_grid_filename => mci_vamp2_get_grid_filename
     procedure :: prepend_grid_path => mci_vamp2_prepend_grid_path
     procedure, public :: declare_flat_dimensions => &
          mci_vamp2_declare_flat_dimensions
     procedure, public :: declare_equivalences => mci_vamp2_declare_equivalences
     procedure, public :: allocate_instance => mci_vamp2_allocate_instance
     procedure, public :: add_pass => mci_vamp2_add_pass
     procedure, public :: update_from_ref => mci_vamp2_update_from_ref
     procedure, public :: update => mci_vamp2_update
     procedure :: write_grids => mci_vamp2_write_grids
     procedure :: read_header => mci_vamp2_read_header
     procedure :: read_data => mci_vamp2_read_data
     procedure, private :: advance_to_data => mci_vamp2_advance_to_data
     procedure, public :: init_integrator => mci_vamp2_init_integrator
     procedure, public :: reset_result => mci_vamp2_reset_result
     procedure, public :: set_calls => mci_vamp2_set_calls
     procedure, private :: init_integration => mci_vamp2_init_integration
     procedure, public :: integrate => mci_vamp2_integrate
     procedure, public :: prepare_simulation => mci_vamp2_prepare_simulation
     procedure, public :: generate_weighted_event => &
          mci_vamp2_generate_weighted_event
     procedure, public :: generate_unweighted_event => &
          mci_vamp2_generate_unweighted_event
     procedure, public :: rebuild_event => mci_vamp2_rebuild_event
  end type mci_vamp2_t

  type, extends (mci_instance_t) :: mci_vamp2_instance_t
     class(mci_vamp2_func_t), allocatable :: func
     real(default), dimension(:), allocatable :: gi
     integer :: n_events = 0
     logical :: event_generated = .false.
     real(default) :: event_weight = 0._default
     real(default) :: event_excess = 0._default
     real(default) :: event_rescale_f_max = 1._default
     real(default), dimension(:), allocatable :: event_x
   contains
     procedure, public :: write => mci_vamp2_instance_write
     procedure, public :: final => mci_vamp2_instance_final
     procedure, public :: init => mci_vamp2_instance_init
     procedure, public :: set_workspace => mci_vamp2_instance_set_workspace
     procedure, public :: compute_weight => mci_vamp2_instance_compute_weight
     procedure, public :: record_integrand => mci_vamp2_instance_record_integrand
     procedure, public :: init_simulation => mci_vamp2_instance_init_simulation
     procedure, public :: final_simulation => mci_vamp2_instance_final_simulation
     procedure, public :: get_event_weight => mci_vamp2_instance_get_event_weight
     procedure, public :: get_event_excess => mci_vamp2_instance_get_event_excess
  end type mci_vamp2_instance_t


  interface operator (.matches.)
     module procedure pass_matches
  end interface operator (.matches.)
  interface operator (.matches.)
     module procedure real_matches
  end interface operator (.matches.)

  interface
    module subroutine mci_vamp2_func_set_workspace (self, instance, sampler)
      class(mci_vamp2_func_t), intent(inout) :: self
      class(mci_vamp2_instance_t), intent(inout), target :: instance
      class(mci_sampler_t), intent(inout), target :: sampler
    end subroutine mci_vamp2_func_set_workspace
    module function mci_vamp2_func_get_probabilities (self) result (gi)
      class(mci_vamp2_func_t), intent(inout) :: self
      real(default), dimension(self%n_channel) :: gi
    end function mci_vamp2_func_get_probabilities
    module function mci_vamp2_func_get_weight (self) result (g)
      class(mci_vamp2_func_t), intent(in) :: self
      real(default) :: g
    end function mci_vamp2_func_get_weight
    module subroutine mci_vamp2_func_set_integrand (self, integrand)
      class(mci_vamp2_func_t), intent(inout) :: self
      real(default), intent(in) :: integrand
    end subroutine mci_vamp2_func_set_integrand
    module subroutine mci_vamp2_func_evaluate_maps (self, x)
      class(mci_vamp2_func_t), intent(inout) :: self
      real(default), dimension(:), intent(in) :: x
    end subroutine mci_vamp2_func_evaluate_maps
    module function mci_vamp2_func_evaluate_func (self, x) result (f)
      class(mci_vamp2_func_t), intent(in) :: self
      real(default), dimension(:), intent(in) :: x
      real(default) :: f
    end function mci_vamp2_func_evaluate_func
    module subroutine list_pass_final (self)
      class(list_pass_t), intent(inout) :: self
    end subroutine list_pass_final
    module subroutine list_pass_add &
         (self, adapt_grids, adapt_weights, final_pass)
      class(list_pass_t), intent(inout) :: self
      logical, intent(in), optional :: adapt_grids, adapt_weights, final_pass
    end subroutine list_pass_add
    module subroutine list_pass_update_from_ref (self, ref, success)
      class(list_pass_t), intent(inout) :: self
      type(list_pass_t), intent(in) :: ref
      logical, intent(out) :: success
    end subroutine list_pass_update_from_ref
    module function list_pass_has_last_integral(self) result (flag)
      class(list_pass_t), intent(in) :: self
      logical :: flag
    end function list_pass_has_last_integral
    module subroutine list_pass_get_last_integral &
         (self, integral, error, efficiency)
      class(list_pass_t), intent(in) :: self
      real(default), intent(out) :: integral
      real(default), intent(out) :: error
      real(default), intent(out) :: efficiency
    end subroutine list_pass_get_last_integral
    module subroutine list_pass_write (self, unit, pacify)
      class(list_pass_t), intent(in) :: self
      integer, intent(in) :: unit
      logical, intent(in), optional :: pacify
    end subroutine list_pass_write
    module subroutine pass_write (self, unit, pacify)
      class(pass_t), intent(in) :: self
      integer, intent(in) :: unit
      logical, intent(in), optional :: pacify
    end subroutine pass_write
    module subroutine pass_read (self, u, n_pass, n_it)
      class(pass_t), intent(out) :: self
      integer, intent(in) :: u, n_pass, n_it
    end subroutine pass_read
    module subroutine pass_configure (pass, n_it, n_calls, n_calls_min)
      class(pass_t), intent(inout) :: pass
      integer, intent(in) :: n_it
      integer, intent(in) :: n_calls
      integer, intent(in) :: n_calls_min
    end subroutine pass_configure
    module function pass_matches (pass, ref) result (ok)
      type(pass_t), intent(in) :: pass, ref
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
    module function pass_get_calls_valid (pass) result (valid)
      class(pass_t), intent(in) :: pass
      integer :: valid
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
    module subroutine mci_vamp2_final (object)
      class(mci_vamp2_t), intent(inout) :: object
    end subroutine mci_vamp2_final
    module subroutine mci_vamp2_write (object, unit, pacify, md5sum_version)
      class(mci_vamp2_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacify
      logical, intent(in), optional :: md5sum_version
    end subroutine mci_vamp2_write
    module subroutine mci_vamp2_compute_md5sum (mci, pacify)
      class(mci_vamp2_t), intent(inout) :: mci
      logical, intent(in), optional :: pacify
    end subroutine mci_vamp2_compute_md5sum
    pure module function mci_vamp2_get_md5sum (mci) result (md5sum)
      class(mci_vamp2_t), intent(in) :: mci
      character(32) :: md5sum
    end function mci_vamp2_get_md5sum
    module subroutine mci_vamp2_startup_message (mci, unit, n_calls)
      class(mci_vamp2_t), intent(in) :: mci
      integer, intent(in), optional :: unit, n_calls
    end subroutine mci_vamp2_startup_message
    module subroutine mci_vamp2_write_log_entry (mci, u)
      class(mci_vamp2_t), intent(in) :: mci
      integer, intent(in) :: u
    end subroutine mci_vamp2_write_log_entry
    module subroutine mci_vamp2_record_index (mci, i_mci)
      class(mci_vamp2_t), intent(inout) :: mci
      integer, intent(in) :: i_mci
    end subroutine mci_vamp2_record_index
    module subroutine mci_vamp2_set_config (mci, config)
      class(mci_vamp2_t), intent(inout) :: mci
      type(mci_vamp2_config_t), intent(in) :: config
    end subroutine mci_vamp2_set_config
    module subroutine mci_vamp2_set_rebuild_flag (mci, rebuild, check_grid_file)
      class(mci_vamp2_t), intent(inout) :: mci
      logical, intent(in) :: rebuild
      logical, intent(in) :: check_grid_file
    end subroutine mci_vamp2_set_rebuild_flag
    module subroutine mci_vamp2_set_grid_filename (mci, name, run_id)
      class(mci_vamp2_t), intent(inout) :: mci
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: run_id
    end subroutine mci_vamp2_set_grid_filename
    module function mci_vamp2_get_grid_filename (mci, binary_grid_format) &
         result (filename)
      class(mci_vamp2_t), intent(in) :: mci
      logical, intent(in), optional :: binary_grid_format
      type(string_t) :: filename
    end function mci_vamp2_get_grid_filename
    module subroutine mci_vamp2_prepend_grid_path (mci, prefix)
      class(mci_vamp2_t), intent(inout) :: mci
      type(string_t), intent(in) :: prefix
    end subroutine mci_vamp2_prepend_grid_path
    module subroutine mci_vamp2_declare_flat_dimensions (mci, dim_flat)
      class(mci_vamp2_t), intent(inout) :: mci
      integer, dimension(:), intent(in) :: dim_flat
    end subroutine mci_vamp2_declare_flat_dimensions
    module subroutine mci_vamp2_declare_equivalences (mci, channel, dim_offset)
      class(mci_vamp2_t), intent(inout) :: mci
      type(phs_channel_t), dimension(:), intent(in) :: channel
      integer, intent(in) :: dim_offset
    end subroutine mci_vamp2_declare_equivalences
    module subroutine mci_vamp2_add_pass &
         (mci, adapt_grids, adapt_weights, final_pass)
      class(mci_vamp2_t), intent(inout) :: mci
      logical, intent(in), optional :: adapt_grids, adapt_weights, final_pass
    end subroutine mci_vamp2_add_pass
    module subroutine mci_vamp2_update_from_ref (mci, mci_ref, success)
      class(mci_vamp2_t), intent(inout) :: mci
      class(mci_t), intent(in) :: mci_ref
      logical, intent(out) :: success
    end subroutine mci_vamp2_update_from_ref
    module subroutine mci_vamp2_update (mci, u, success)
      class(mci_vamp2_t), intent(inout) :: mci
      integer, intent(in) :: u
      logical, intent(out) :: success
    end subroutine mci_vamp2_update
    module subroutine mci_vamp2_write_grids (mci)
      class(mci_vamp2_t), intent(in) :: mci
    end subroutine mci_vamp2_write_grids
    module subroutine mci_vamp2_read_header (mci, success)
      class(mci_vamp2_t), intent(inout) :: mci
      logical, intent(out) :: success
    end subroutine mci_vamp2_read_header
    module subroutine mci_vamp2_read_data (mci)
      class(mci_vamp2_t), intent(inout) :: mci
    end subroutine mci_vamp2_read_data
    module subroutine mci_vamp2_advance_to_data (mci, u, binary_grid_format)
      class(mci_vamp2_t), intent(in) :: mci
      integer, intent(in) :: u
      logical, intent(out) :: binary_grid_format
    end subroutine mci_vamp2_advance_to_data
    module subroutine mci_vamp2_init_integrator (mci)
      class(mci_vamp2_t), intent(inout) :: mci
    end subroutine mci_vamp2_init_integrator
    module subroutine mci_vamp2_reset_result (mci)
      class(mci_vamp2_t), intent(inout) :: mci
    end subroutine mci_vamp2_reset_result
    module subroutine mci_vamp2_set_calls (mci, n_calls)
      class(mci_vamp2_t), intent(inout) :: mci
      integer :: n_calls
    end subroutine mci_vamp2_set_calls
    module subroutine mci_vamp2_init_integration (mci, n_it, n_calls, instance)
      class(mci_vamp2_t), intent(inout) :: mci
      integer, intent(in) :: n_it
      integer, intent(in) :: n_calls
      class(mci_instance_t), intent(inout) :: instance
    end subroutine mci_vamp2_init_integration
    module subroutine mci_vamp2_integrate (mci, instance, sampler, &
         n_it, n_calls, results, pacify)
      class(mci_vamp2_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout), target :: instance
      class(mci_sampler_t), intent(inout), target :: sampler
      integer, intent(in) :: n_it
      integer, intent(in) :: n_calls
      class(mci_results_t), intent(inout), optional :: results
      logical, intent(in), optional :: pacify
    end subroutine mci_vamp2_integrate
    module subroutine mci_vamp2_prepare_simulation (mci)
      class(mci_vamp2_t), intent(inout) :: mci
    end subroutine mci_vamp2_prepare_simulation
    module subroutine mci_vamp2_generate_weighted_event (mci, instance, sampler)
      class(mci_vamp2_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout), target :: instance
      class(mci_sampler_t), intent(inout), target :: sampler
    end subroutine mci_vamp2_generate_weighted_event
    module subroutine mci_vamp2_generate_unweighted_event &
         (mci, instance, sampler)
      class(mci_vamp2_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout), target :: instance
      class(mci_sampler_t), intent(inout), target :: sampler
    end subroutine mci_vamp2_generate_unweighted_event
    module subroutine mci_vamp2_rebuild_event (mci, instance, sampler, state)
      class(mci_vamp2_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout) :: instance
      class(mci_sampler_t), intent(inout) :: sampler
      class(mci_state_t), intent(in) :: state
    end subroutine mci_vamp2_rebuild_event
    module subroutine mci_vamp2_instance_write (object, unit, pacify)
      class(mci_vamp2_instance_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacify
    end subroutine mci_vamp2_instance_write
    module subroutine mci_vamp2_instance_final (object)
      class(mci_vamp2_instance_t), intent(inout) :: object
    end subroutine mci_vamp2_instance_final
    module subroutine mci_vamp2_instance_set_workspace (instance, sampler)
      class(mci_vamp2_instance_t), intent(inout), target :: instance
      class(mci_sampler_t), intent(inout), target :: sampler
    end subroutine mci_vamp2_instance_set_workspace
    module subroutine mci_vamp2_instance_compute_weight (mci, c)
      class(mci_vamp2_instance_t), intent(inout) :: mci
      integer, intent(in) :: c
    end subroutine mci_vamp2_instance_compute_weight
    module subroutine mci_vamp2_instance_record_integrand (mci, integrand)
      class(mci_vamp2_instance_t), intent(inout) :: mci
      real(default), intent(in) :: integrand
    end subroutine mci_vamp2_instance_record_integrand
    module subroutine mci_vamp2_instance_init_simulation &
         (instance, safety_factor)
      class(mci_vamp2_instance_t), intent(inout) :: instance
      real(default), intent(in), optional :: safety_factor
    end subroutine mci_vamp2_instance_init_simulation
    module subroutine mci_vamp2_instance_final_simulation (instance)
      class(mci_vamp2_instance_t), intent(inout) :: instance
    end subroutine mci_vamp2_instance_final_simulation
    module function mci_vamp2_instance_get_event_weight (mci) result (weight)
      class(mci_vamp2_instance_t), intent(in) :: mci
      real(default) :: weight
    end function mci_vamp2_instance_get_event_weight
    module function mci_vamp2_instance_get_event_excess (mci) result (excess)
      class(mci_vamp2_instance_t), intent(in) :: mci
      real(default) :: excess
    end function mci_vamp2_instance_get_event_excess
  end interface

contains

  subroutine mci_vamp2_allocate_instance (mci, mci_instance)
    class(mci_vamp2_t), intent(in) :: mci
    class(mci_instance_t), intent(out), pointer :: mci_instance
    allocate (mci_vamp2_instance_t :: mci_instance)
  end subroutine mci_vamp2_allocate_instance

  subroutine mci_vamp2_instance_init (mci_instance, mci)
    class(mci_vamp2_instance_t), intent(out) :: mci_instance
    class(mci_t), intent(in), target :: mci
    call mci_instance%base_init (mci)
    allocate (mci_instance%gi(mci%n_channel), source=0._default)
    allocate (mci_instance%event_x(mci%n_dim), source=0._default)
    allocate (mci_vamp2_func_t :: mci_instance%func)
    call mci_instance%func%init (n_dim = mci%n_dim, n_channel = mci%n_channel)
  end subroutine mci_vamp2_instance_init


end module mci_vamp2
