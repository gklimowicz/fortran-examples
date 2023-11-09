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

module process_mci

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use cputime
  use rng_base
  use mci_base
  use variables
  use integration_results
  use process_libraries
  use phs_base
  use process_counter
  use process_config

  implicit none
  private

  public :: process_mci_entry_t
  public :: mci_work_t

  integer, parameter, public :: REAL_FULL = 0
  integer, parameter, public :: REAL_SINGULAR = 1
  integer, parameter, public :: REAL_FINITE = 2

  type :: process_mci_entry_t
     integer :: i_mci = 0
     integer, dimension(:), allocatable :: i_component
     integer :: process_type = PRC_UNKNOWN
     integer :: n_par = 0
     integer :: n_par_sf = 0
     integer :: n_par_phs = 0
     character(32) :: md5sum = ""
     integer :: pass = 0
     integer :: n_it = 0
     integer :: n_calls = 0
     logical :: activate_timer = .false.
     real(default) :: error_threshold = 0
     class(mci_t), allocatable :: mci
     type(process_counter_t) :: counter
     type(integration_results_t) :: results
     logical :: negative_weights = .false.
     logical :: combined_integration = .false.
     integer :: real_partition_type = REAL_FULL
   contains
     procedure :: final => process_mci_entry_final
     procedure :: write => process_mci_entry_write
     procedure :: configure => process_mci_entry_configure
     procedure :: create_component_list => &
        process_mci_entry_create_component_list
     procedure :: set_parameters => process_mci_entry_set_parameters
     procedure :: compute_md5sum => process_mci_entry_compute_md5sum
     procedure :: sampler_test => process_mci_entry_sampler_test
     procedure :: integrate => process_mci_entry_integrate
     procedure :: final_integration => process_mci_entry_final_integration
     procedure :: get_time => process_mci_entry_get_time
     procedure :: time_message => process_mci_entry_time_message
     procedure :: prepare_simulation => process_mci_entry_prepare_simulation
     procedure :: generate_weighted_event => &
          process_mci_entry_generate_weighted_event
     procedure :: generate_unweighted_event => &
          process_mci_entry_generate_unweighted_event
     procedure :: has_integral => process_mci_entry_has_integral
     procedure :: get_integral => process_mci_entry_get_integral
     procedure :: get_error => process_mci_entry_get_error
     procedure :: get_accuracy => process_mci_entry_get_accuracy
     procedure :: get_chi2 => process_mci_entry_get_chi2
     procedure :: get_efficiency => process_mci_entry_get_efficiency
     procedure :: get_md5sum => process_mci_entry_get_md5sum
  end type process_mci_entry_t

  type :: mci_work_t
     type(process_mci_entry_t), pointer :: config => null ()
     real(default), dimension(:), allocatable :: x
     class(mci_instance_t), pointer :: mci => null ()
     type(process_counter_t) :: counter
     logical :: keep_failed_events = .false.
     integer :: n_event_dropped = 0
   contains
     procedure :: write => mci_work_write
     procedure :: final => mci_work_final
     procedure :: init => mci_work_init
     procedure :: set => mci_work_set
     procedure :: set_x_strfun => mci_work_set_x_strfun
     procedure :: set_x_process => mci_work_set_x_process
     procedure :: get_active_components => mci_work_get_active_components
     procedure :: get_x_strfun => mci_work_get_x_strfun
     procedure :: get_x_process => mci_work_get_x_process
     procedure :: init_simulation => mci_work_init_simulation
     procedure :: final_simulation => mci_work_final_simulation
     procedure :: reset_counter => mci_work_reset_counter
     procedure :: record_call => mci_work_record_call
     procedure :: get_counter => mci_work_get_counter
  end type mci_work_t


  interface
    module subroutine process_mci_entry_final (object)
      class(process_mci_entry_t), intent(inout) :: object
    end subroutine process_mci_entry_final
    module subroutine process_mci_entry_write (object, unit, pacify)
      class(process_mci_entry_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacify
    end subroutine process_mci_entry_write
    module subroutine process_mci_entry_configure (mci_entry, mci_template, &
         process_type, i_mci, i_component, component, &
         n_sfpar, rng_factory)
      class(process_mci_entry_t), intent(inout) :: mci_entry
      class(mci_t), intent(in), allocatable :: mci_template
      integer, intent(in) :: process_type
      integer, intent(in) :: i_mci
      integer, intent(in) :: i_component
      type(process_component_t), intent(in), target :: component
      integer, intent(in) :: n_sfpar
      class(rng_factory_t), intent(inout) :: rng_factory
    end subroutine process_mci_entry_configure
    module subroutine process_mci_entry_create_component_list (mci_entry, &
       i_component, component_config)
      class (process_mci_entry_t), intent(inout) :: mci_entry
      integer, intent(in) :: i_component
      type(process_component_def_t), intent(in) :: component_config
    end subroutine process_mci_entry_create_component_list
    module subroutine process_mci_entry_set_parameters (mci_entry, var_list)
      class(process_mci_entry_t), intent(inout) :: mci_entry
      type(var_list_t), intent(in) :: var_list
    end subroutine process_mci_entry_set_parameters
    module subroutine process_mci_entry_compute_md5sum (mci_entry, &
         config, component, beam_config)
      class(process_mci_entry_t), intent(inout) :: mci_entry
      type(process_config_data_t), intent(in) :: config
      type(process_component_t), dimension(:), intent(in) :: component
      type(process_beam_config_t), intent(in) :: beam_config
    end subroutine process_mci_entry_compute_md5sum
    module subroutine process_mci_entry_sampler_test &
         (mci_entry, mci_sampler, n_calls)
      class(process_mci_entry_t), intent(inout) :: mci_entry
      class(mci_sampler_t), intent(inout), target :: mci_sampler
      integer, intent(in) :: n_calls
    end subroutine process_mci_entry_sampler_test
    module subroutine process_mci_entry_integrate (mci_entry, mci_instance, &
           mci_sampler, n_it, n_calls, &
         adapt_grids, adapt_weights, final, pacify, &
         nlo_type)
      class(process_mci_entry_t), intent(inout) :: mci_entry
      class(mci_instance_t), intent(inout) :: mci_instance
      class(mci_sampler_t), intent(inout) :: mci_sampler
      integer, intent(in) :: n_it
      integer, intent(in) :: n_calls
      logical, intent(in), optional :: adapt_grids
      logical, intent(in), optional :: adapt_weights
      logical, intent(in), optional :: final, pacify
      integer, intent(in), optional :: nlo_type
    end subroutine process_mci_entry_integrate
    module subroutine process_mci_entry_final_integration (mci_entry)
      class(process_mci_entry_t), intent(inout) :: mci_entry
    end subroutine process_mci_entry_final_integration
    module subroutine process_mci_entry_get_time (mci_entry, time, sample)
      class(process_mci_entry_t), intent(in) :: mci_entry
      type(time_t), intent(out) :: time
      integer, intent(in) :: sample
    end subroutine process_mci_entry_get_time
    module subroutine process_mci_entry_time_message (mci_entry)
      class(process_mci_entry_t), intent(in) :: mci_entry
    end subroutine process_mci_entry_time_message
    module subroutine process_mci_entry_prepare_simulation (mci_entry)
      class(process_mci_entry_t), intent(inout) :: mci_entry
    end subroutine process_mci_entry_prepare_simulation
    module subroutine process_mci_entry_generate_weighted_event (mci_entry, &
        mci_instance, mci_sampler, keep_failed)
      class(process_mci_entry_t), intent(inout) :: mci_entry
      class(mci_instance_t), intent(inout) :: mci_instance
      class(mci_sampler_t), intent(inout) :: mci_sampler
      logical, intent(in) :: keep_failed
    end subroutine process_mci_entry_generate_weighted_event
    module subroutine process_mci_entry_generate_unweighted_event &
         (mci_entry, mci_instance, mci_sampler)
      class(process_mci_entry_t), intent(inout) :: mci_entry
      class(mci_instance_t), intent(inout) :: mci_instance
      class(mci_sampler_t), intent(inout) :: mci_sampler
    end subroutine process_mci_entry_generate_unweighted_event
    module function process_mci_entry_has_integral (mci_entry) result (flag)
      class(process_mci_entry_t), intent(in) :: mci_entry
      logical :: flag
    end function process_mci_entry_has_integral
    module function process_mci_entry_get_integral (mci_entry) result (integral)
      class(process_mci_entry_t), intent(in) :: mci_entry
      real(default) :: integral
    end function process_mci_entry_get_integral
    module function process_mci_entry_get_error (mci_entry) result (error)
      class(process_mci_entry_t), intent(in) :: mci_entry
      real(default) :: error
    end function process_mci_entry_get_error
    module function process_mci_entry_get_accuracy (mci_entry) result (accuracy)
      class(process_mci_entry_t), intent(in) :: mci_entry
      real(default) :: accuracy
    end function process_mci_entry_get_accuracy
    module function process_mci_entry_get_chi2 (mci_entry) result (chi2)
      class(process_mci_entry_t), intent(in) :: mci_entry
      real(default) :: chi2
    end function process_mci_entry_get_chi2
    module function process_mci_entry_get_efficiency &
         (mci_entry) result (efficiency)
      class(process_mci_entry_t), intent(in) :: mci_entry
      real(default) :: efficiency
    end function process_mci_entry_get_efficiency
    pure module function process_mci_entry_get_md5sum (entry) result (md5sum)
      class(process_mci_entry_t), intent(in) :: entry
      character(32) :: md5sum
    end function process_mci_entry_get_md5sum
    module subroutine mci_work_write (mci_work, unit, testflag)
      class(mci_work_t), intent(in) :: mci_work
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine mci_work_write
    module subroutine mci_work_final (mci_work)
      class(mci_work_t), intent(inout) :: mci_work
    end subroutine mci_work_final
    module subroutine mci_work_init (mci_work, mci_entry)
      class(mci_work_t), intent(out) :: mci_work
      type(process_mci_entry_t), intent(in), target :: mci_entry
    end subroutine mci_work_init
    module subroutine mci_work_set (mci_work, x)
      class(mci_work_t), intent(inout) :: mci_work
      real(default), dimension(:), intent(in) :: x
    end subroutine mci_work_set
    module subroutine mci_work_set_x_strfun (mci_work, x)
      class(mci_work_t), intent(inout) :: mci_work
      real(default), dimension(:), intent(in) :: x
    end subroutine mci_work_set_x_strfun
    module subroutine mci_work_set_x_process (mci_work, x)
      class(mci_work_t), intent(inout) :: mci_work
      real(default), dimension(:), intent(in) :: x
    end subroutine mci_work_set_x_process
    module function mci_work_get_active_components &
         (mci_work) result (i_component)
      class(mci_work_t), intent(in) :: mci_work
      integer, dimension(:), allocatable :: i_component
    end function mci_work_get_active_components
    pure module function mci_work_get_x_strfun (mci_work) result (x)
      class(mci_work_t), intent(in) :: mci_work
      real(default), dimension(mci_work%config%n_par_sf) :: x
    end function mci_work_get_x_strfun
    pure module function mci_work_get_x_process (mci_work) result (x)
      class(mci_work_t), intent(in) :: mci_work
      real(default), dimension(mci_work%config%n_par_phs) :: x
    end function mci_work_get_x_process
    module subroutine mci_work_final_simulation (mci_work)
      class(mci_work_t), intent(inout) :: mci_work
    end subroutine mci_work_final_simulation
    module subroutine mci_work_init_simulation &
         (mci_work, safety_factor, keep_failed_events)
      class(mci_work_t), intent(inout) :: mci_work
      real(default), intent(in), optional :: safety_factor
      logical, intent(in), optional :: keep_failed_events
    end subroutine mci_work_init_simulation
    module subroutine mci_work_reset_counter (mci_work)
      class(mci_work_t), intent(inout) :: mci_work
    end subroutine mci_work_reset_counter
    module subroutine mci_work_record_call (mci_work, status)
      class(mci_work_t), intent(inout) :: mci_work
      integer, intent(in) :: status
    end subroutine mci_work_record_call
    pure module function mci_work_get_counter (mci_work) result (counter)
      class(mci_work_t), intent(in) :: mci_work
      type(process_counter_t) :: counter
    end function mci_work_get_counter
  end interface

end module process_mci
