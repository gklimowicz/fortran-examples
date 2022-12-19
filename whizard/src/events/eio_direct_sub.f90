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

submodule (eio_direct) eio_direct_s

  use io_units
  use diagnostics
  use cputime

  implicit none

contains

  module subroutine eio_direct_write (object, unit)
    class(eio_direct_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Event direct access:"
    if (object%i_evt_set) then
       write (u, "(3x,A,1x,I0)")  "i_evt =", object%i_evt
    else
       write (u, "(3x,A)")  "i_evt = [undefined]"
    end if
    write (u, "(3x,A,1x,I0)")  "i_prc =", object%i_prc
    write (u, "(3x,A,1x,I0)")  "i_mci =", object%i_prc
    write (u, "(3x,A,1x,I0)")  "i_term =", object%i_prc
    write (u, "(3x,A,1x,I0)")  "channel =", object%i_prc
    if (object%passed_set) then
       write (u, "(3x,A,1x,L1)")  "passed =", object%passed
    else
       write (u, "(3x,A)")  "passed = [N/A]"
    end if
    call object%pset%write (u)
  end subroutine eio_direct_write

  module subroutine eio_direct_final (object)
    class(eio_direct_t), intent(inout) :: object
    call object%pset%final ()
  end subroutine eio_direct_final

  module subroutine eio_direct_init_out &
       (eio, sample, data, success, extension)
    class(eio_direct_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    if (present (success))  success = .true.
  end subroutine eio_direct_init_out

  module subroutine eio_direct_init_in &
       (eio, sample, data, success, extension)
    class(eio_direct_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    if (present (success))  success = .true.
  end subroutine eio_direct_init_in

  module subroutine eio_direct_switch_inout (eio, success)
    class(eio_direct_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    if (present (success))  success = .true.
  end subroutine eio_direct_switch_inout

  module subroutine eio_direct_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_direct_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    integer, intent(in) :: i_prc
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    type(particle_set_t), pointer :: pset_ptr
    call eio%pset%final ()
    if (event%has_index ()) then
       call eio%set_event_index (event%get_index ())
    else
       call eio%reset_event_index ()
    end if
    if (present (passed)) then
       eio%passed = passed
       eio%passed_set = .true.
    else
       eio%passed_set = .false.
    end if
    pset_ptr => event%get_particle_set_ptr ()
    if (associated (pset_ptr)) then
       eio%i_prc = i_prc
       eio%pset = pset_ptr
    end if
  end subroutine eio_direct_output

  module subroutine eio_direct_input_i_prc (eio, i_prc, iostat)
    class(eio_direct_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    i_prc = eio%i_prc
    iostat = 0
  end subroutine eio_direct_input_i_prc

  module subroutine eio_direct_input_event (eio, event, iostat, event_handle)
    class(eio_direct_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    call event%select (eio%i_mci, eio%i_term, eio%channel)
    if (eio%has_event_index ()) then
       call event%set_index (eio%get_event_index ())
    else
       call event%reset_index ()
    end if
    call event%set_hard_particle_set (eio%pset)
  end subroutine eio_direct_input_event

  module subroutine eio_direct_skip (eio, iostat)
    class(eio_direct_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    iostat = 0
  end subroutine eio_direct_skip

  module function eio_direct_has_event_index (eio) result (flag)
    class(eio_direct_t), intent(in) :: eio
    logical :: flag
    flag = eio%i_evt_set
  end function eio_direct_has_event_index

  module function eio_direct_get_event_index (eio) result (index)
    class(eio_direct_t), intent(in) :: eio
    integer :: index
    if (eio%has_event_index ()) then
       index = eio%i_evt
    else
       index = 0
    end if
  end function eio_direct_get_event_index

  module function eio_direct_passed_known (eio) result (flag)
    class(eio_direct_t), intent(in) :: eio
    logical :: flag
    flag = eio%passed_set
  end function eio_direct_passed_known

  module function eio_direct_has_passed (eio) result (flag)
    class(eio_direct_t), intent(in) :: eio
    logical :: flag
    if (eio%passed_known ()) then
       flag = eio%passed
    else
       flag = .true.
    end if
  end function eio_direct_has_passed

  module function eio_direct_get_n_in (eio) result (n_in)
    class(eio_direct_t), intent(in) :: eio
    integer :: n_in
    n_in = eio%pset%get_n_in ()
  end function eio_direct_get_n_in

  module function eio_direct_get_n_out (eio) result (n_out)
    class(eio_direct_t), intent(in) :: eio
    integer :: n_out
    n_out = eio%pset%get_n_out ()
  end function eio_direct_get_n_out

  module function eio_direct_get_n_tot (eio) result (n_tot)
    class(eio_direct_t), intent(in) :: eio
    integer :: n_tot
    n_tot = eio%pset%get_n_tot ()
  end function eio_direct_get_n_tot

  module subroutine eio_direct_get_momentum_array (eio, p)
    class(eio_direct_t), intent(in) :: eio
    type(vector4_t), dimension(:), allocatable, intent(out) :: p
    integer :: n
    n = eio%get_n_tot ()
    allocate (p (n))
    p(:) = eio%pset%get_momenta ()
  end subroutine eio_direct_get_momentum_array

  module subroutine eio_direct_init_direct &
       (eio, n_beam, n_in, n_rem, n_vir, n_out, pdg, model)
    class(eio_direct_t), intent(out) :: eio
    integer, intent(in) :: n_beam
    integer, intent(in) :: n_in
    integer, intent(in) :: n_rem
    integer, intent(in) :: n_vir
    integer, intent(in) :: n_out
    integer, dimension(:), intent(in) :: pdg
    class(model_data_t), intent(in), target :: model
    call eio%pset%init_direct (n_beam, n_in, n_rem, n_vir, n_out, pdg, model)
  end subroutine eio_direct_init_direct

  module subroutine eio_direct_set_event_index (eio, index)
    class(eio_direct_t), intent(inout) :: eio
    integer, intent(in) :: index
    eio%i_evt = index
    eio%i_evt_set = .true.
  end subroutine eio_direct_set_event_index

  module subroutine eio_direct_reset_event_index (eio)
    class(eio_direct_t), intent(inout) :: eio
    eio%i_evt_set = .false.
  end subroutine eio_direct_reset_event_index

  module subroutine eio_direct_set_selection_indices &
       (eio, i_prc, i_mci, i_term, channel)
    class(eio_direct_t), intent(inout) :: eio
    integer, intent(in) :: i_prc
    integer, intent(in) :: i_mci
    integer, intent(in) :: i_term
    integer, intent(in) :: channel
    eio%i_prc = i_prc
    eio%i_mci = i_mci
    eio%i_term = i_term
    eio%channel = channel
  end subroutine eio_direct_set_selection_indices

  module subroutine eio_direct_set_momentum_single (eio, i, p, p2, on_shell)
    class(eio_direct_t), intent(inout) :: eio
    integer, intent(in) :: i
    type(vector4_t), intent(in) :: p
    real(default), intent(in), optional :: p2
    logical, intent(in), optional :: on_shell
    call eio%pset%set_momentum (i, p, p2, on_shell)
  end subroutine eio_direct_set_momentum_single

  module subroutine eio_direct_set_momentum_all (eio, p, p2, on_shell)
    class(eio_direct_t), intent(inout) :: eio
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), dimension(:), intent(in), optional :: p2
    logical, intent(in), optional :: on_shell
    call eio%pset%set_momentum (p, p2, on_shell)
  end subroutine eio_direct_set_momentum_all


end submodule eio_direct_s

