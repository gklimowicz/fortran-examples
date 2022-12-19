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

submodule (eio_weights) eio_weights_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine eio_weights_set_parameters (eio, pacify)
    class(eio_weights_t), intent(inout) :: eio
    logical, intent(in), optional :: pacify
    if (present (pacify))  eio%pacify = pacify
    eio%extension = "weights.dat"
  end subroutine eio_weights_set_parameters

  module subroutine eio_weights_write (object, unit)
    class(eio_weights_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Weight stream:"
    if (object%writing) then
       write (u, "(3x,A,A)")  "Writing to file   = ", char (object%filename)
       write (u, "(3x,A,L1)") "Reduced I/O prec. = ", object%pacify
    else
       write (u, "(3x,A)")  "[closed]"
    end if
  end subroutine eio_weights_write

  module subroutine eio_weights_final (object)
    class(eio_weights_t), intent(inout) :: object
    if (object%writing) then
       write (msg_buffer, "(A,A,A)")  "Events: closing weight stream file '", &
            char (object%filename), "'"
       call msg_message ()
       close (object%unit)
       object%writing = .false.
    end if
  end subroutine eio_weights_final

  module subroutine eio_weights_init_out &
       (eio, sample, data, success, extension)
    class(eio_weights_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    if (present(extension)) then
       eio%extension = extension
    else
       eio%extension = "weights.dat"
    end if
    eio%filename = sample // "." // eio%extension
    eio%unit = free_unit ()
    write (msg_buffer, "(A,A,A)")  "Events: writing to weight stream file '", &
         char (eio%filename), "'"
    call msg_message ()
    eio%writing = .true.
    open (eio%unit, file = char (eio%filename), &
         action = "write", status = "replace")
    if (present (success))  success = .true.
  end subroutine eio_weights_init_out

  module subroutine eio_weights_init_in &
       (eio, sample, data, success, extension)
    class(eio_weights_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    call msg_bug ("Weight stream: event input not supported")
    if (present (success))  success = .false.
  end subroutine eio_weights_init_in

  module subroutine eio_weights_switch_inout (eio, success)
    class(eio_weights_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    call msg_bug ("Weight stream: in-out switch not supported")
    if (present (success))  success = .false.
  end subroutine eio_weights_switch_inout

  module subroutine eio_weights_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_weights_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    integer, intent(in) :: i_prc
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    integer :: n_alt, i
    real(default) :: weight, sqme_ref, sqme_prc
    logical :: evt_pacify, evt_passed
    evt_pacify = eio%pacify;  if (present (pacify))  evt_pacify = pacify
    evt_passed = .true.;  if (present (passed))  evt_passed = passed
    if (eio%writing) then
       if (evt_passed) then
          weight = event%get_weight_prc ()
          sqme_ref = event%get_sqme_ref ()
          sqme_prc = event%get_sqme_prc ()
          n_alt = event%get_n_alt ()
1         format (I0,3(1x,ES17.10),3(1x,I0))
2         format (I0,3(1x,ES15.8),3(1x,I0))
          if (evt_pacify) then
             write (eio%unit, 2)  0, weight, sqme_prc, sqme_ref, &
                  i_prc
          else
             write (eio%unit, 1)  0, weight, sqme_prc, sqme_ref, &
                  i_prc
          end if
          do i = 1, n_alt
             weight = event%get_weight_alt(i)
             sqme_prc = event%get_sqme_alt(i)
             if (evt_pacify) then
                write (eio%unit, 2)  i, weight, sqme_prc
             else
                write (eio%unit, 1)  i, weight, sqme_prc
             end if
          end do
       end if
    else
       call eio%write ()
       call msg_fatal ("Weight stream file is not open for writing")
    end if
  end subroutine eio_weights_output

  module subroutine eio_weights_input_i_prc (eio, i_prc, iostat)
    class(eio_weights_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    call msg_bug ("Weight stream: event input not supported")
    i_prc = 0
    iostat = 1
  end subroutine eio_weights_input_i_prc

  module subroutine eio_weights_input_event &
       (eio, event, iostat, event_handle)
    class(eio_weights_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    call msg_bug ("Weight stream: event input not supported")
    iostat = 1
  end subroutine eio_weights_input_event

  module subroutine eio_weights_skip (eio, iostat)
    class(eio_weights_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    iostat = 0
  end subroutine eio_weights_skip


end submodule eio_weights_s

