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

module mci_midpoint

  use kinds, only: default
  use phs_base

  use mci_base

  implicit none
  private

  public :: mci_midpoint_t
  public :: mci_midpoint_instance_t

  type, extends (mci_t) :: mci_midpoint_t
     integer :: n_dim_binned = 0
     logical, dimension(:), allocatable :: dim_is_binned
     logical :: calls_known = .false.
     integer :: n_calls = 0
     integer :: n_calls_pos = 0
     integer :: n_calls_nul = 0
     integer :: n_calls_neg = 0
     real(default) :: integral_pos = 0
     real(default) :: integral_neg = 0
     integer, dimension(:), allocatable :: n_bin
     logical :: max_known = .false.
     real(default) :: max = 0
     real(default) :: min = 0
     real(default) :: max_abs = 0
     real(default) :: min_abs = 0
   contains
     procedure :: final => mci_midpoint_final
     procedure :: write => mci_midpoint_write
     procedure :: startup_message => mci_midpoint_startup_message
     procedure :: write_log_entry => mci_midpoint_write_log_entry
     procedure :: compute_md5sum => mci_midpoint_compute_md5sum
     procedure :: set_dimensions => mci_midpoint_set_dimensions
     procedure :: declare_flat_dimensions => mci_midpoint_declare_flat_dimensions
     procedure :: declare_equivalences => mci_midpoint_ignore_equivalences
     procedure :: allocate_instance => mci_midpoint_allocate_instance
     procedure :: integrate => mci_midpoint_integrate
     procedure :: prepare_simulation => mci_midpoint_ignore_prepare_simulation
     procedure :: generate_weighted_event => mci_midpoint_generate_weighted_event
     procedure :: generate_unweighted_event => &
          mci_midpoint_generate_unweighted_event
     procedure :: rebuild_event => mci_midpoint_rebuild_event
  end type mci_midpoint_t

  type, extends (mci_instance_t) :: mci_midpoint_instance_t
     type(mci_midpoint_t), pointer :: mci => null ()
     logical :: max_known = .false.
     real(default) :: max = 0
     real(default) :: min = 0
     real(default) :: max_abs = 0
     real(default) :: min_abs = 0
     real(default) :: safety_factor = 1
     real(default) :: excess_weight = 0
   contains
     procedure :: write => mci_midpoint_instance_write
     procedure :: final => mci_midpoint_instance_final
     procedure :: init => mci_midpoint_instance_init
     procedure :: get_max => mci_midpoint_instance_get_max
     procedure :: set_max => mci_midpoint_instance_set_max
     procedure :: compute_weight => mci_midpoint_instance_compute_weight
     procedure :: record_integrand => mci_midpoint_instance_record_integrand
     procedure :: init_simulation => mci_midpoint_instance_init_simulation
     procedure :: final_simulation => mci_midpoint_instance_final_simulation
     procedure :: get_event_excess => mci_midpoint_instance_get_event_excess
  end type mci_midpoint_instance_t


  interface
    module subroutine mci_midpoint_final (object)
      class(mci_midpoint_t), intent(inout) :: object
    end subroutine mci_midpoint_final
    module subroutine mci_midpoint_write (object, unit, pacify, md5sum_version)
      class(mci_midpoint_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacify
      logical, intent(in), optional :: md5sum_version
    end subroutine mci_midpoint_write
    module subroutine mci_midpoint_startup_message (mci, unit, n_calls)
      class(mci_midpoint_t), intent(in) :: mci
      integer, intent(in), optional :: unit, n_calls
    end subroutine mci_midpoint_startup_message
    module subroutine mci_midpoint_write_log_entry (mci, u)
      class(mci_midpoint_t), intent(in) :: mci
      integer, intent(in) :: u
    end subroutine mci_midpoint_write_log_entry
    module subroutine mci_midpoint_compute_md5sum (mci, pacify)
      class(mci_midpoint_t), intent(inout) :: mci
      logical, intent(in), optional :: pacify
    end subroutine mci_midpoint_compute_md5sum
    module subroutine mci_midpoint_set_dimensions (mci, n_dim, n_channel)
      class(mci_midpoint_t), intent(inout) :: mci
      integer, intent(in) :: n_dim
      integer, intent(in) :: n_channel
    end subroutine mci_midpoint_set_dimensions
    module subroutine mci_midpoint_declare_flat_dimensions (mci, dim_flat)
      class(mci_midpoint_t), intent(inout) :: mci
      integer, dimension(:), intent(in) :: dim_flat
    end subroutine mci_midpoint_declare_flat_dimensions
    module subroutine mci_midpoint_ignore_equivalences &
         (mci, channel, dim_offset)
      class(mci_midpoint_t), intent(inout) :: mci
      type(phs_channel_t), dimension(:), intent(in) :: channel
      integer, intent(in) :: dim_offset
    end subroutine mci_midpoint_ignore_equivalences
    module subroutine mci_midpoint_integrate (mci, instance, sampler, n_it, &
         n_calls, results, pacify)
      class(mci_midpoint_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout), target :: instance
      class(mci_sampler_t), intent(inout), target :: sampler
      integer, intent(in) :: n_it
      integer, intent(in) :: n_calls
      logical, intent(in), optional :: pacify
      class(mci_results_t), intent(inout), optional :: results
    end subroutine mci_midpoint_integrate
    module subroutine mci_midpoint_ignore_prepare_simulation (mci)
      class(mci_midpoint_t), intent(inout) :: mci
    end subroutine mci_midpoint_ignore_prepare_simulation
    module subroutine mci_midpoint_generate_weighted_event &
         (mci, instance, sampler)
      class(mci_midpoint_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout), target :: instance
      class(mci_sampler_t), intent(inout), target :: sampler
    end subroutine mci_midpoint_generate_weighted_event
    module subroutine mci_midpoint_generate_unweighted_event &
         (mci, instance, sampler)
      class(mci_midpoint_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout), target :: instance
      class(mci_sampler_t), intent(inout), target :: sampler
    end subroutine mci_midpoint_generate_unweighted_event
    module subroutine mci_midpoint_rebuild_event (mci, instance, sampler, state)
      class(mci_midpoint_t), intent(inout) :: mci
      class(mci_instance_t), intent(inout) :: instance
      class(mci_sampler_t), intent(inout) :: sampler
      class(mci_state_t), intent(in) :: state
    end subroutine mci_midpoint_rebuild_event
    module subroutine mci_midpoint_instance_write (object, unit, pacify)
      class(mci_midpoint_instance_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacify
    end subroutine mci_midpoint_instance_write
    module subroutine mci_midpoint_instance_final (object)
      class(mci_midpoint_instance_t), intent(inout) :: object
    end subroutine mci_midpoint_instance_final
    module subroutine mci_midpoint_instance_init (mci_instance, mci)
      class(mci_midpoint_instance_t), intent(out) :: mci_instance
      class(mci_t), intent(in), target :: mci
    end subroutine mci_midpoint_instance_init
    module subroutine mci_midpoint_instance_get_max (instance)
      class(mci_midpoint_instance_t), intent(inout) :: instance
    end subroutine mci_midpoint_instance_get_max
    module subroutine mci_midpoint_instance_set_max (instance)
      class(mci_midpoint_instance_t), intent(inout) :: instance
    end subroutine mci_midpoint_instance_set_max
    module subroutine mci_midpoint_instance_compute_weight (mci, c)
      class(mci_midpoint_instance_t), intent(inout) :: mci
      integer, intent(in) :: c
    end subroutine mci_midpoint_instance_compute_weight
    module subroutine mci_midpoint_instance_record_integrand (mci, integrand)
      class(mci_midpoint_instance_t), intent(inout) :: mci
      real(default), intent(in) :: integrand
    end subroutine mci_midpoint_instance_record_integrand
    module subroutine mci_midpoint_instance_init_simulation &
         (instance, safety_factor)
      class(mci_midpoint_instance_t), intent(inout) :: instance
      real(default), intent(in), optional :: safety_factor
    end subroutine mci_midpoint_instance_init_simulation
    module subroutine mci_midpoint_instance_final_simulation (instance)
      class(mci_midpoint_instance_t), intent(inout) :: instance
    end subroutine mci_midpoint_instance_final_simulation
    module function mci_midpoint_instance_get_event_excess (mci) result (excess)
      class(mci_midpoint_instance_t), intent(in) :: mci
      real(default) :: excess
    end function mci_midpoint_instance_get_event_excess
  end interface

contains

  subroutine mci_midpoint_allocate_instance (mci, mci_instance)
    class(mci_midpoint_t), intent(in) :: mci
    class(mci_instance_t), intent(out), pointer :: mci_instance
    allocate (mci_midpoint_instance_t :: mci_instance)
  end subroutine mci_midpoint_allocate_instance


end module mci_midpoint
