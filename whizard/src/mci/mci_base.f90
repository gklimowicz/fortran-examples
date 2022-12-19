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

module mci_base

  use kinds
  use cputime
  use phs_base
  use rng_base

  implicit none
  private

  public :: mci_t
  public :: mci_instance_t
  public :: mci_state_t
  public :: mci_sampler_t
  public :: mci_results_t

  type, abstract :: mci_t
     integer :: n_dim = 0
     integer :: n_channel = 0
     integer :: n_chain = 0
     integer, dimension(:), allocatable :: chain
     real(default), dimension(:), allocatable :: chain_weights
     character(32) :: md5sum = ""
     logical :: integral_known = .false.
     logical :: error_known = .false.
     logical :: efficiency_known = .false.
     real(default) :: integral = 0
     real(default) :: error = 0
     real(default) :: efficiency = 0
     logical :: use_timer = .false.
     type(timer_t) :: timer
     class(rng_t), allocatable :: rng
   contains
     procedure :: base_final => mci_final
     procedure (mci_final), deferred :: final
     procedure :: base_write => mci_write
     procedure (mci_write), deferred :: write
     procedure (mci_startup_message), deferred :: startup_message
     procedure :: base_startup_message => mci_startup_message
     procedure(mci_write_log_entry), deferred :: write_log_entry
     procedure(mci_compute_md5sum), deferred :: compute_md5sum
     procedure :: record_index => mci_record_index
     procedure :: set_dimensions => mci_set_dimensions
     procedure :: base_set_dimensions => mci_set_dimensions
     procedure (mci_declare_flat_dimensions), deferred :: declare_flat_dimensions
     procedure (mci_declare_equivalences), deferred :: declare_equivalences
     procedure :: declare_chains => mci_declare_chains
     procedure :: collect_chain_weights => mci_collect_chain_weights
     procedure :: has_chains => mci_has_chains
     procedure :: write_chain_weights => mci_write_chain_weights
     procedure :: set_md5sum => mci_set_md5sum
     procedure :: add_pass => mci_add_pass
     procedure (mci_allocate_instance), deferred :: allocate_instance
     procedure :: import_rng => mci_import_rng
     procedure :: set_timer => mci_set_timer
     procedure :: start_timer => mci_start_timer
     procedure :: stop_timer => mci_stop_timer
     procedure :: sampler_test => mci_sampler_test
     procedure (mci_integrate), deferred :: integrate
     procedure (mci_prepare_simulation), deferred :: prepare_simulation
     procedure (mci_generate), deferred :: generate_weighted_event
     procedure (mci_generate), deferred :: generate_unweighted_event
     procedure (mci_rebuild), deferred :: rebuild_event
     procedure :: pacify => mci_pacify
     procedure :: get_integral => mci_get_integral
     procedure :: get_error => mci_get_error
     procedure :: get_efficiency => mci_get_efficiency
     procedure :: get_time => mci_get_time
     procedure :: get_md5sum => mci_get_md5sum
  end type mci_t

  type, abstract :: mci_instance_t
     logical :: valid = .false.
     real(default), dimension(:), allocatable :: w
     real(default), dimension(:), allocatable :: f
     real(default), dimension(:,:), allocatable :: x
     integer :: selected_channel = 0
     real(default) :: mci_weight = 0
     real(default) :: integrand  = 0
     logical :: negative_weights = .false.
     integer :: n_dropped = 0
   contains
     procedure (mci_instance_write), deferred :: write
     procedure (mci_instance_final), deferred :: final
     procedure (mci_instance_base_init), deferred :: init
     procedure :: base_init => mci_instance_base_init
     procedure :: set_channel_weights => mci_instance_set_channel_weights
     procedure (mci_instance_compute_weight), deferred :: compute_weight
     procedure (mci_instance_record_integrand), deferred :: record_integrand
     procedure :: evaluate => mci_instance_evaluate
     procedure (mci_instance_init_simulation), deferred :: init_simulation
     procedure (mci_instance_final_simulation), deferred :: final_simulation
     procedure :: fetch => mci_instance_fetch
     procedure :: get_value => mci_instance_get_value
     procedure :: get_event_weight => mci_instance_get_value
     procedure (mci_instance_get_event_excess), deferred :: get_event_excess
     procedure :: get_n_event_dropped => mci_instance_get_n_event_dropped
     procedure :: reset_n_event_dropped => mci_instance_reset_n_event_dropped
     procedure :: record_event_dropped => mci_instance_record_event_dropped
     procedure :: store => mci_instance_store
     procedure :: recall => mci_instance_recall
  end type mci_instance_t

  type :: mci_state_t
     integer :: selected_channel = 0
     real(default), dimension(:), allocatable :: x_in
     real(default) :: val
   contains
     procedure :: write => mci_state_write
  end type mci_state_t

  type, abstract :: mci_sampler_t
   contains
     procedure (mci_sampler_write), deferred :: write
     procedure (mci_sampler_evaluate), deferred :: evaluate
     procedure (mci_sampler_is_valid), deferred :: is_valid
     procedure (mci_sampler_rebuild), deferred :: rebuild
     procedure (mci_sampler_fetch), deferred :: fetch
  end type mci_sampler_t

  type, abstract :: mci_results_t
   contains
     procedure (mci_results_write), deferred :: write
     procedure (mci_results_write_verbose), deferred :: write_verbose
     generic :: record => record_simple, record_extended
     procedure (mci_results_record_simple), deferred :: record_simple
     procedure (mci_results_record_extended), deferred :: record_extended
  end type mci_results_t


  abstract interface
     subroutine mci_write_log_entry (mci, u)
       import
       class(mci_t), intent(in) :: mci
       integer, intent(in) :: u
     end subroutine mci_write_log_entry
  end interface

  abstract interface
     subroutine mci_compute_md5sum (mci, pacify)
       import
       class(mci_t), intent(inout) :: mci
       logical, intent(in), optional :: pacify
     end subroutine mci_compute_md5sum
  end interface

  abstract interface
     subroutine mci_declare_flat_dimensions (mci, dim_flat)
       import
       class(mci_t), intent(inout) :: mci
       integer, dimension(:), intent(in) :: dim_flat
     end subroutine mci_declare_flat_dimensions
  end interface

  abstract interface
     subroutine mci_declare_equivalences (mci, channel, dim_offset)
       import
       class(mci_t), intent(inout) :: mci
       type(phs_channel_t), dimension(:), intent(in) :: channel
       integer, intent(in) :: dim_offset
     end subroutine mci_declare_equivalences
  end interface

  abstract interface
     subroutine mci_allocate_instance (mci, mci_instance)
       import
       class(mci_t), intent(in) :: mci
       class(mci_instance_t), intent(out), pointer :: mci_instance
     end subroutine mci_allocate_instance
  end interface

  abstract interface
     subroutine mci_integrate (mci, instance, sampler, &
          n_it, n_calls, results, pacify)
       import
       class(mci_t), intent(inout) :: mci
       class(mci_instance_t), intent(inout), target :: instance
       class(mci_sampler_t), intent(inout), target :: sampler
       integer, intent(in) :: n_it
       integer, intent(in) :: n_calls
       logical, intent(in), optional :: pacify
       class(mci_results_t), intent(inout), optional :: results
     end subroutine mci_integrate
  end interface

  abstract interface
     subroutine mci_prepare_simulation (mci)
       import
       class(mci_t), intent(inout) :: mci
     end subroutine mci_prepare_simulation
  end interface

  abstract interface
     subroutine mci_generate (mci, instance, sampler)
       import
       class(mci_t), intent(inout) :: mci
       class(mci_instance_t), intent(inout), target :: instance
       class(mci_sampler_t), intent(inout), target :: sampler
     end subroutine mci_generate
  end interface

  abstract interface
     subroutine mci_rebuild (mci, instance, sampler, state)
       import
       class(mci_t), intent(inout) :: mci
       class(mci_instance_t), intent(inout) :: instance
       class(mci_sampler_t), intent(inout) :: sampler
       class(mci_state_t), intent(in) :: state
     end subroutine mci_rebuild
  end interface

  abstract interface
     subroutine mci_instance_write (object, unit, pacify)
       import
       class(mci_instance_t), intent(in) :: object
       integer, intent(in), optional :: unit
       logical, intent(in), optional :: pacify
     end subroutine mci_instance_write
  end interface

  abstract interface
     subroutine mci_instance_final (object)
       import
       class(mci_instance_t), intent(inout) :: object
     end subroutine mci_instance_final
  end interface

  abstract interface
     subroutine mci_instance_compute_weight (mci, c)
       import
       class(mci_instance_t), intent(inout) :: mci
       integer, intent(in) :: c
     end subroutine mci_instance_compute_weight
  end interface

  abstract interface
     subroutine mci_instance_record_integrand (mci, integrand)
       import
       class(mci_instance_t), intent(inout) :: mci
       real(default), intent(in) :: integrand
     end subroutine mci_instance_record_integrand
  end interface

  abstract interface
     subroutine mci_instance_init_simulation (instance, safety_factor)
       import
       class(mci_instance_t), intent(inout) :: instance
       real(default), intent(in), optional :: safety_factor
     end subroutine mci_instance_init_simulation
  end interface

  abstract interface
     subroutine mci_instance_final_simulation (instance)
       import
       class(mci_instance_t), intent(inout) :: instance
     end subroutine mci_instance_final_simulation
  end interface

  abstract interface
     function mci_instance_get_event_excess (mci) result (excess)
       import
       class(mci_instance_t), intent(in) :: mci
       real(default) :: excess
     end function mci_instance_get_event_excess
  end interface

  abstract interface
     subroutine mci_sampler_write (object, unit, testflag)
       import
       class(mci_sampler_t), intent(in) :: object
       integer, intent(in), optional :: unit
       logical, intent(in), optional :: testflag
     end subroutine mci_sampler_write
  end interface

  abstract interface
     subroutine mci_sampler_evaluate (sampler, c, x_in, val, x, f)
       import
       class(mci_sampler_t), intent(inout) :: sampler
       integer, intent(in) :: c
       real(default), dimension(:), intent(in) :: x_in
       real(default), intent(out) :: val
       real(default), dimension(:,:), intent(out) :: x
       real(default), dimension(:), intent(out) :: f
     end subroutine mci_sampler_evaluate
  end interface

  abstract interface
     function mci_sampler_is_valid (sampler) result (valid)
       import
       class(mci_sampler_t), intent(in) :: sampler
       logical :: valid
     end function mci_sampler_is_valid
  end interface

  abstract interface
     subroutine mci_sampler_rebuild (sampler, c, x_in, val, x, f)
       import
       class(mci_sampler_t), intent(inout) :: sampler
       integer, intent(in) :: c
       real(default), dimension(:), intent(in) :: x_in
       real(default), intent(in) :: val
       real(default), dimension(:,:), intent(out) :: x
       real(default), dimension(:), intent(out) :: f
     end subroutine mci_sampler_rebuild
  end interface

  abstract interface
     subroutine mci_sampler_fetch (sampler, val, x, f)
       import
       class(mci_sampler_t), intent(in) :: sampler
       real(default), intent(out) :: val
       real(default), dimension(:,:), intent(out) :: x
       real(default), dimension(:), intent(out) :: f
     end subroutine mci_sampler_fetch
  end interface

  abstract interface
     subroutine mci_results_write (object, unit, suppress)
       import
       class(mci_results_t), intent(in) :: object
       integer, intent(in), optional :: unit
       logical, intent(in), optional :: suppress
     end subroutine mci_results_write

     subroutine mci_results_write_verbose (object, unit)
       import
       class(mci_results_t), intent(in) :: object
       integer, intent(in), optional :: unit
     end subroutine mci_results_write_verbose
  end interface

  abstract interface
     subroutine mci_results_record_simple (object, n_it, &
          n_calls, integral, error, efficiency, chain_weights, suppress)
       import
       class(mci_results_t), intent(inout) :: object
       integer, intent(in) :: n_it
       integer, intent(in) :: n_calls
       real(default), intent(in) :: integral
       real(default), intent(in) :: error
       real(default), intent(in) :: efficiency
       real(default), dimension(:), intent(in), optional :: chain_weights
       logical, intent(in), optional :: suppress
     end subroutine mci_results_record_simple

     subroutine mci_results_record_extended (object, n_it, n_calls,&
          & n_calls_valid, integral, error, efficiency, efficiency_pos,&
          & efficiency_neg, chain_weights, suppress)
       import
       class(mci_results_t), intent(inout) :: object
       integer, intent(in) :: n_it
       integer, intent(in) :: n_calls
       integer, intent(in) :: n_calls_valid
       real(default), intent(in) :: integral
       real(default), intent(in) :: error
       real(default), intent(in) :: efficiency
       real(default), intent(in) :: efficiency_pos
       real(default), intent(in) :: efficiency_neg
       real(default), dimension(:), intent(in), optional :: chain_weights
       logical, intent(in), optional :: suppress
     end subroutine mci_results_record_extended
  end interface


  interface
    module subroutine mci_final (object)
      class(mci_t), intent(inout) :: object
    end subroutine mci_final
    module subroutine mci_write (object, unit, pacify, md5sum_version)
      class(mci_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacify
      logical, intent(in), optional :: md5sum_version
    end subroutine mci_write
    module subroutine mci_startup_message (mci, unit, n_calls)
      class(mci_t), intent(in) :: mci
      integer, intent(in), optional :: unit, n_calls
    end subroutine mci_startup_message
    module subroutine mci_record_index (mci, i_mci)
      class(mci_t), intent(inout) :: mci
      integer, intent(in) :: i_mci
    end subroutine mci_record_index
    module subroutine mci_set_dimensions (mci, n_dim, n_channel)
      class(mci_t), intent(inout) :: mci
      integer, intent(in) :: n_dim
      integer, intent(in) :: n_channel
    end subroutine mci_set_dimensions
    module subroutine mci_declare_chains (mci, chain)
      class(mci_t), intent(inout) :: mci
      integer, dimension(:), intent(in) :: chain
    end subroutine mci_declare_chains
    module subroutine mci_collect_chain_weights (mci, weight)
      class(mci_t), intent(inout) :: mci
      real(default), dimension(:), intent(in) :: weight
    end subroutine mci_collect_chain_weights
    module function mci_has_chains (mci) result (flag)
      class(mci_t), intent(in) :: mci
      logical :: flag
    end function mci_has_chains
    module subroutine mci_write_chain_weights (mci, unit)
      class(mci_t), intent(in) :: mci
      integer, intent(in), optional :: unit
    end subroutine mci_write_chain_weights
    module subroutine mci_set_md5sum (mci, md5sum)
      class(mci_t), intent(inout) :: mci
      character(32), intent(in) :: md5sum
    end subroutine mci_set_md5sum
    module subroutine mci_add_pass (mci, adapt_grids, adapt_weights, final_pass)
      class(mci_t), intent(inout) :: mci
      logical, intent(in), optional :: adapt_grids
      logical, intent(in), optional :: adapt_weights
      logical, intent(in), optional :: final_pass
    end subroutine mci_add_pass
    module subroutine mci_import_rng (mci, rng)
      class(mci_t), intent(inout) :: mci
      class(rng_t), intent(inout), allocatable :: rng
    end subroutine mci_import_rng
    module subroutine mci_set_timer (mci, active)
      class(mci_t), intent(inout) :: mci
      logical, intent(in) :: active
    end subroutine mci_set_timer
    module subroutine mci_start_timer (mci)
      class(mci_t), intent(inout) :: mci
    end subroutine mci_start_timer
    module subroutine mci_stop_timer (mci)
      class(mci_t), intent(inout) :: mci
    end subroutine mci_stop_timer
    module subroutine mci_sampler_test (mci, sampler, n_calls)
      class(mci_t), intent(inout) :: mci
      class(mci_sampler_t), intent(inout), target :: sampler
      integer, intent(in) :: n_calls
    end subroutine mci_sampler_test
    module subroutine mci_pacify (object, efficiency_reset, error_reset)
      class(mci_t), intent(inout) :: object
      logical, intent(in), optional :: efficiency_reset, error_reset
    end subroutine mci_pacify
    module function mci_get_integral (mci) result (integral)
      class(mci_t), intent(in) :: mci
      real(default) :: integral
    end function mci_get_integral
    module function mci_get_error (mci) result (error)
      class(mci_t), intent(in) :: mci
      real(default) :: error
    end function mci_get_error
    module function mci_get_efficiency (mci) result (efficiency)
      class(mci_t), intent(in) :: mci
      real(default) :: efficiency
    end function mci_get_efficiency
    module function mci_get_time (mci) result (time)
      class(mci_t), intent(in) :: mci
      real(default) :: time
    end function mci_get_time
    pure module function mci_get_md5sum (mci) result (md5sum)
      class(mci_t), intent(in) :: mci
      character(32) :: md5sum
    end function mci_get_md5sum
    module subroutine mci_instance_base_init (mci_instance, mci)
      class(mci_instance_t), intent(out) :: mci_instance
      class(mci_t), intent(in), target :: mci
    end subroutine mci_instance_base_init
    module subroutine mci_instance_set_channel_weights &
         (mci_instance, weights, sum_non_zero)
      class(mci_instance_t), intent(inout) :: mci_instance
      real(default), dimension(:), intent(in) :: weights
      logical, intent(out), optional :: sum_non_zero
    end subroutine mci_instance_set_channel_weights
    module subroutine mci_instance_evaluate (mci, sampler, c, x)
      class(mci_instance_t), intent(inout) :: mci
      class(mci_sampler_t), intent(inout) :: sampler
      integer, intent(in) :: c
      real(default), dimension(:), intent(in) :: x
    end subroutine mci_instance_evaluate
    module subroutine mci_instance_fetch (mci, sampler, c)
      class(mci_instance_t), intent(inout) :: mci
      class(mci_sampler_t), intent(in) :: sampler
      integer, intent(in) :: c
    end subroutine mci_instance_fetch
    module function mci_instance_get_value (mci) result (value)
      class(mci_instance_t), intent(in) :: mci
      real(default) :: value
    end function mci_instance_get_value
    module function mci_instance_get_n_event_dropped (mci) result (n_dropped)
      class(mci_instance_t), intent(in) :: mci
      integer :: n_dropped
    end function mci_instance_get_n_event_dropped
    module subroutine mci_instance_reset_n_event_dropped (mci)
      class(mci_instance_t), intent(inout) :: mci
    end subroutine mci_instance_reset_n_event_dropped
    module subroutine mci_instance_record_event_dropped (mci)
      class(mci_instance_t), intent(inout) :: mci
    end subroutine mci_instance_record_event_dropped
    module subroutine mci_state_write (object, unit)
      class(mci_state_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine mci_state_write
    module subroutine mci_instance_store (mci, state)
      class(mci_instance_t), intent(in) :: mci
      class(mci_state_t), intent(out) :: state
    end subroutine mci_instance_store
    module subroutine mci_instance_recall (mci, sampler, state)
      class(mci_instance_t), intent(inout) :: mci
      class(mci_sampler_t), intent(inout) :: sampler
      class(mci_state_t), intent(in) :: state
    end subroutine mci_instance_recall
  end interface

end module mci_base
