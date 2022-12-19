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

submodule (eio_checkpoints) eio_checkpoints_s

  use io_units
  use diagnostics

  implicit none

  character(*), parameter :: &
     checkpoint_head = "| % complete | events generated | events remaining &
     &| time remaining"
  character(*), parameter :: &
     checkpoint_bar  = "|==================================================&
     &=================|"
  character(*), parameter :: &
     checkpoint_fmt  = "('   ',F5.1,T16,I9,T35,I9,T58,A)"

contains

  module subroutine eio_checkpoints_set_parameters (eio, checkpoint, blank)
    class(eio_checkpoints_t), intent(inout) :: eio
    integer, intent(in) :: checkpoint
    logical, intent(in), optional :: blank
    eio%val = checkpoint
    if (present (blank))  eio%blank = blank
  end subroutine eio_checkpoints_set_parameters

  module subroutine eio_checkpoints_write (object, unit)
    class(eio_checkpoints_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (object%active) then
       write (u, "(1x,A)")  "Event-sample checkpoints:  active"
       write (u, "(3x,A,I0)")  "interval  = ", object%val
       write (u, "(3x,A,I0)")  "n_events  = ", object%n_events
       write (u, "(3x,A,I0)")  "n_read    = ", object%n_read
       write (u, "(3x,A,I0)")  "n_current = ", object%i_evt
       write (u, "(3x,A,L1)")  "blanking  = ", object%blank
       call object%timer%write (u)
    else
       write (u, "(1x,A)")  "Event-sample checkpoints:  off"
    end if
  end subroutine eio_checkpoints_write

  module subroutine eio_checkpoints_final (object)
    class(eio_checkpoints_t), intent(inout) :: object
    object%active = .false.
  end subroutine eio_checkpoints_final

  module subroutine eio_checkpoints_init_out &
       (eio, sample, data, success, extension)
    class(eio_checkpoints_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    if (present (data)) then
       if (eio%val > 0) then
          eio%active = .true.
          eio%i_evt = 0
          eio%n_read = 0
          eio%n_events = data%n_evt * data%nlo_multiplier
       end if
    end if
    if (present (success))  success = .true.
  end subroutine eio_checkpoints_init_out

  module subroutine eio_checkpoints_init_in &
       (eio, sample, data, success, extension)
    class(eio_checkpoints_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    call msg_bug ("Event checkpoints: event input not supported")
    if (present (success))  success = .false.
  end subroutine eio_checkpoints_init_in

  module subroutine eio_checkpoints_switch_inout (eio, success)
    class(eio_checkpoints_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    call msg_bug ("Event checkpoints: in-out switch not supported")
    if (present (success))  success = .false.
  end subroutine eio_checkpoints_switch_inout

  module subroutine eio_checkpoints_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_checkpoints_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    integer, intent(in) :: i_prc
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    logical :: rd
    rd = .false.;  if (present (reading))  rd = reading
    if (eio%active) then
       if (.not. eio%running)  call eio%startup ()
       if (eio%running) then
          eio%i_evt = eio%i_evt + 1
          if (rd) then
             eio%n_read = eio%n_read + 1
          else if (mod (eio%i_evt, eio%val) == 0) then
             call eio%message (eio%blank)
          end if
          if (eio%i_evt == eio%n_events)  call eio%shutdown ()
       end if
    end if
  end subroutine eio_checkpoints_output

  module subroutine eio_checkpoints_startup (eio)
    class(eio_checkpoints_t), intent(inout) :: eio
    if (eio%active .and. eio%i_evt < eio%n_events) then
       call msg_message ("")
       call msg_message (checkpoint_bar)
       call msg_message (checkpoint_head)
       call msg_message (checkpoint_bar)
       write (msg_buffer, checkpoint_fmt) 0., 0, eio%n_events - eio%i_evt, "???"
       call msg_message ()
       eio%running = .true.
       call eio%timer%start ()
    end if
  end subroutine eio_checkpoints_startup

  module subroutine eio_checkpoints_message (eio, testflag)
    class(eio_checkpoints_t), intent(inout) :: eio
    logical, intent(in), optional :: testflag
    real :: t
    type(time_t) :: time_remaining
    type(string_t) :: time_string
    call eio%timer%stop ()
    t = eio%timer
    call eio%timer%restart ()
    time_remaining = &
         nint (t / (eio%i_evt - eio%n_read) * (eio%n_events - eio%i_evt))
    time_string = time_remaining%to_string_ms (blank = testflag)
    write (msg_buffer, checkpoint_fmt) &
         100 * ((eio%i_evt - eio%n_read) / real (eio%n_events - eio%n_read)), &
         eio%i_evt - eio%n_read, &
         eio%n_events - eio%i_evt, &
         char (time_string)
    call msg_message ()
  end subroutine eio_checkpoints_message

  module subroutine eio_checkpoints_shutdown (eio)
    class(eio_checkpoints_t), intent(inout) :: eio
    if (mod (eio%i_evt, eio%val) /= 0) then
       write (msg_buffer, checkpoint_fmt) &
            100., eio%i_evt - eio%n_read, 0, "0m:00s"
       call msg_message ()
    end if
    call msg_message (checkpoint_bar)
    call msg_message ("")
    eio%running = .false.
  end subroutine eio_checkpoints_shutdown

  module subroutine eio_checkpoints_input_i_prc (eio, i_prc, iostat)
    class(eio_checkpoints_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    call msg_bug ("Event checkpoints: event input not supported")
    i_prc = 0
    iostat = 1
  end subroutine eio_checkpoints_input_i_prc

  module subroutine eio_checkpoints_input_event &
       (eio, event, iostat, event_handle)
    class(eio_checkpoints_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    call msg_bug ("Event checkpoints: event input not supported")
    iostat = 1
  end subroutine eio_checkpoints_input_event

  module subroutine eio_checkpoints_skip (eio, iostat)
    class(eio_checkpoints_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    iostat = 0
  end subroutine eio_checkpoints_skip


end submodule eio_checkpoints_s

