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

submodule (eio_callback) eio_callback_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine eio_callback_set_parameters &
       (eio, callback, count_interval)
    class(eio_callback_t), intent(inout) :: eio
    class(event_callback_t), intent(in) :: callback
    integer, intent(in) :: count_interval
    allocate (eio%callback, source = callback)
    eio%n_interval = count_interval
  end subroutine eio_callback_set_parameters

  module subroutine eio_callback_write (object, unit)
    class(eio_callback_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Event-sample callback:"
    write (u, "(3x,A,I0)")  "interval  = ", object%n_interval
    write (u, "(3x,A,I0)")  "evt count = ", object%i_evt
!        call object%timer%write (u)
  end subroutine eio_callback_write

  module subroutine eio_callback_final (object)
    class(eio_callback_t), intent(inout) :: object
  end subroutine eio_callback_final

  module subroutine eio_callback_init_out &
       (eio, sample, data, success, extension)
    class(eio_callback_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    eio%i_evt = 0
    eio%i_interval = 0
    if (present (success))  success = .true.
  end subroutine eio_callback_init_out

  module subroutine eio_callback_init_in &
       (eio, sample, data, success, extension)
    class(eio_callback_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    call msg_bug ("Event callback: event input not supported")
    if (present (success))  success = .false.
  end subroutine eio_callback_init_in

  module subroutine eio_callback_switch_inout (eio, success)
    class(eio_callback_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    call msg_bug ("Event callback: in-out switch not supported")
    if (present (success))  success = .false.
  end subroutine eio_callback_switch_inout

  module subroutine eio_callback_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_callback_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    integer, intent(in) :: i_prc
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    eio%i_evt = eio%i_evt + 1
    if (eio%n_interval > 0) then
       eio%i_interval = eio%i_interval + 1
       if (eio%i_interval >= eio%n_interval) then
          call eio%callback%proc (eio%i_evt, event)
          eio%i_interval = 0
       end if
    end if
  end subroutine eio_callback_output

  module subroutine eio_callback_input_i_prc (eio, i_prc, iostat)
    class(eio_callback_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    call msg_bug ("Event callback: event input not supported")
    i_prc = 0
    iostat = 1
  end subroutine eio_callback_input_i_prc

  module subroutine eio_callback_input_event &
       (eio, event, iostat, event_handle)
    class(eio_callback_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    call msg_bug ("Event callback: event input not supported")
    iostat = 1
  end subroutine eio_callback_input_event

  module subroutine eio_callback_skip (eio, iostat)
    class(eio_callback_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    iostat = 0
  end subroutine eio_callback_skip


end submodule eio_callback_s

