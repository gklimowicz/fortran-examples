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

submodule (mci_none) mci_none_s

  use diagnostics, only: msg_message, msg_fatal

  implicit none

contains

  module subroutine mci_none_final (object)
    class(mci_none_t), intent(inout) :: object
  end subroutine mci_none_final

  module subroutine mci_none_write (object, unit, pacify, md5sum_version)
    class(mci_none_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: pacify
    logical, intent(in), optional :: md5sum_version
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)") "Integrator: non-functional dummy"
  end subroutine mci_none_write

  module subroutine mci_none_startup_message (mci, unit, n_calls)
    class(mci_none_t), intent(in) :: mci
    integer, intent(in), optional :: unit, n_calls
    call msg_message ("Integrator: none")
  end subroutine mci_none_startup_message

  module subroutine mci_none_write_log_entry (mci, u)
    class(mci_none_t), intent(in) :: mci
    integer, intent(in) :: u
    write (u, "(1x,A)")  "MC Integrator is none (no-op)"
  end subroutine mci_none_write_log_entry

  module subroutine mci_none_compute_md5sum (mci, pacify)
    class(mci_none_t), intent(inout) :: mci
    logical, intent(in), optional :: pacify
  end subroutine mci_none_compute_md5sum

  module subroutine mci_none_ignore_flat_dimensions (mci, dim_flat)
    class(mci_none_t), intent(inout) :: mci
    integer, dimension(:), intent(in) :: dim_flat
  end subroutine mci_none_ignore_flat_dimensions

  module subroutine mci_none_ignore_equivalences (mci, channel, dim_offset)
    class(mci_none_t), intent(inout) :: mci
    type(phs_channel_t), dimension(:), intent(in) :: channel
    integer, intent(in) :: dim_offset
  end subroutine mci_none_ignore_equivalences

  module subroutine mci_none_integrate (mci, instance, sampler, n_it, &
       n_calls, results, pacify)
    class(mci_none_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout), target :: instance
    class(mci_sampler_t), intent(inout), target :: sampler
    integer, intent(in) :: n_it
    integer, intent(in) :: n_calls
    logical, intent(in), optional :: pacify
    class(mci_results_t), intent(inout), optional :: results
    call msg_fatal &
         ("Integration: attempt to integrate with the 'mci_none' method")
  end subroutine mci_none_integrate

  module subroutine mci_none_ignore_prepare_simulation (mci)
    class(mci_none_t), intent(inout) :: mci
  end subroutine mci_none_ignore_prepare_simulation

  module subroutine mci_none_generate_no_event (mci, instance, sampler)
    class(mci_none_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout), target :: instance
    class(mci_sampler_t), intent(inout), target :: sampler
    call msg_fatal ("Integration: attempt to generate event " // &
         "with the 'mci_none' method")
  end subroutine mci_none_generate_no_event

  module subroutine mci_none_rebuild_event (mci, instance, sampler, state)
    class(mci_none_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout) :: instance
    class(mci_sampler_t), intent(inout) :: sampler
    class(mci_state_t), intent(in) :: state
  end subroutine mci_none_rebuild_event

  module subroutine mci_none_instance_write (object, unit, pacify)
    class(mci_none_instance_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: pacify
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)") "Integrator instance: non-functional dummy"
  end subroutine mci_none_instance_write

  module subroutine mci_none_instance_final (object)
    class(mci_none_instance_t), intent(inout) :: object
  end subroutine mci_none_instance_final

  module subroutine mci_none_instance_init (mci_instance, mci)
    class(mci_none_instance_t), intent(out) :: mci_instance
    class(mci_t), intent(in), target :: mci
  end subroutine mci_none_instance_init

  module subroutine mci_none_instance_compute_weight (mci, c)
    class(mci_none_instance_t), intent(inout) :: mci
    integer, intent(in) :: c
    call msg_fatal ("Integration: attempt to compute weight " // &
         "with the 'mci_none' method")
  end subroutine mci_none_instance_compute_weight

  module subroutine mci_none_instance_record_integrand (mci, integrand)
    class(mci_none_instance_t), intent(inout) :: mci
    real(default), intent(in) :: integrand
  end subroutine mci_none_instance_record_integrand

  module subroutine mci_none_instance_init_simulation (instance, safety_factor)
    class(mci_none_instance_t), intent(inout) :: instance
    real(default), intent(in), optional :: safety_factor
  end subroutine mci_none_instance_init_simulation

  module subroutine mci_none_instance_final_simulation (instance)
    class(mci_none_instance_t), intent(inout) :: instance
  end subroutine mci_none_instance_final_simulation

  module function mci_none_instance_get_event_excess (mci) result (excess)
    class(mci_none_instance_t), intent(in) :: mci
    real(default) :: excess
    excess = 0
  end function mci_none_instance_get_event_excess


end submodule mci_none_s

