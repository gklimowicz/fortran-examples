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

submodule (eio_dump) eio_dump_s

  use io_units
  use diagnostics
  use format_utils, only: write_separator
  use format_utils, only: pac_fmt
  use format_defs, only: FMT_16, FMT_19

  implicit none

contains

  module subroutine eio_dump_set_parameters (eio, extension, &
       pacify, weights, compressed, summary, screen, unit)
    class(eio_dump_t), intent(inout) :: eio
    type(string_t), intent(in), optional :: extension
    logical, intent(in), optional :: pacify
    logical, intent(in), optional :: weights
    logical, intent(in), optional :: compressed
    logical, intent(in), optional :: summary
    logical, intent(in), optional :: screen
    integer, intent(in), optional :: unit
    if (present (pacify))  eio%pacify = pacify
    if (present (weights))  eio%weights = weights
    if (present (compressed))  eio%compressed = compressed
    if (present (summary))  eio%summary = summary
    if (present (screen))  eio%screen = screen
    if (present (unit))  eio%unit = unit
    eio%extension = "pset.dat"
    if (present (extension))  eio%extension = extension
  end subroutine eio_dump_set_parameters

  module subroutine eio_dump_write (object, unit)
    class(eio_dump_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Dump event stream:"
    if (object%writing) then
       write (u, "(3x,A,L1)") "Screen output     = ", object%screen
       write (u, "(3x,A,A,A)")  "Writing to file   = '", char (object%filename), "'"
       write (u, "(3x,A,L1)") "Reduced I/O prec. = ", object%pacify
       write (u, "(3x,A,L1)") "Show weights/sqme = ", object%weights
       write (u, "(3x,A,L1)") "Compressed        = ", object%compressed
       write (u, "(3x,A,L1)") "Summary           = ", object%summary
    else
       write (u, "(3x,A)")  "[closed]"
    end if
  end subroutine eio_dump_write

  module subroutine eio_dump_final (object)
    class(eio_dump_t), intent(inout) :: object
    if (object%screen) then
       write (msg_buffer, "(A,A,A)")  "Events: display complete"
       call msg_message ()
       object%screen = .false.
    end if
    if (object%writing) then
       if (object%filename /= "") then
          write (msg_buffer, "(A,A,A)")  "Events: closing event dump file '", &
               char (object%filename), "'"
          call msg_message ()
          close (object%unit)
       end if
       object%writing = .false.
    end if
  end subroutine eio_dump_final

  module subroutine eio_dump_init_out &
       (eio, sample, data, success, extension)
    class(eio_dump_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    if (present(extension)) then
       eio%extension = extension
    else
       eio%extension = "pset.dat"
    end if
    if (sample == "" .and. eio%unit /= 0) then
       eio%filename = ""
       eio%writing = .true.
    else if (sample /= "") then
       eio%filename = sample // "." // eio%extension
       eio%unit = free_unit ()
       write (msg_buffer, "(A,A,A)")  "Events: writing to event dump file '", &
            char (eio%filename), "'"
       call msg_message ()
       eio%writing = .true.
       open (eio%unit, file = char (eio%filename), &
            action = "write", status = "replace")
    end if
    if (eio%screen) then
       write (msg_buffer, "(A,A,A)")  "Events: display on standard output"
       call msg_message ()
    end if
    eio%count = 0
    if (present (success))  success = .true.
  end subroutine eio_dump_init_out

  module subroutine eio_dump_init_in &
       (eio, sample, data, success, extension)
    class(eio_dump_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    call msg_bug ("Event dump: event input not supported")
    if (present (success))  success = .false.
  end subroutine eio_dump_init_in

  module subroutine eio_dump_switch_inout (eio, success)
    class(eio_dump_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    call msg_bug ("Event dump: in-out switch not supported")
    if (present (success))  success = .false.
  end subroutine eio_dump_switch_inout

  module subroutine eio_dump_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_dump_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    integer, intent(in) :: i_prc
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    character(len=7) :: fmt
    eio%count = eio%count + 1
    if (present (pacify)) then
       call pac_fmt (fmt, FMT_19, FMT_16, pacify)
    else
       call pac_fmt (fmt, FMT_19, FMT_16, eio%pacify)
    end if
    if (eio%writing)  call dump (eio%unit)
    if (eio%screen) then
       call dump (output_unit)
       if (logfile_unit () > 0)  call dump (logfile_unit ())
    end if
  contains
    subroutine dump (u)
      integer, intent(in) :: u
      integer :: i
      call write_separator (u, 2)
      write (u, "(1x,A,I0)", advance="no")  "Event"
      if (event%has_index ()) then
         write (u, "(1x,'#',I0)")  event%get_index ()
      else
         write (u, *)
      end if
      call write_separator (u, 2)
      write (u, "(1x,A,1x,I0)")  "count  =", eio%count
      if (present (passed)) then
         write (u, "(1x,A,1x,L1)")  "passed =", passed
      else
         write (u, "(1x,A)")  "passed = [N/A]"
      end if
      write (u, "(1x,A,1x,I0)")  "prc id =", i_prc
      if (eio%weights) then
         call write_separator (u)
         if (event%sqme_ref_known) then
            write (u, "(1x,A," // fmt // ")")  "sqme (ref)   = ", &
                 event%sqme_ref
         else
            write (u, "(1x,A)")  "sqme (ref)    = [undefined]"
         end if
         if (event%sqme_prc_known) then
            write (u, "(1x,A," // fmt // ")")  "sqme (prc)   = ", &
                 event%sqme_prc
         else
            write (u, "(1x,A)")  "sqme (prc)    = [undefined]"
         end if
         if (event%weight_ref_known) then
            write (u, "(1x,A," // fmt // ")")  "weight (ref) = ", &
                 event%weight_ref
         else
            write (u, "(1x,A)")  "weight (ref) = [undefined]"
         end if
         if (event%weight_prc_known) then
            write (u, "(1x,A," // fmt // ")")  "weight (prc) = ", &
                 event%weight_prc
         else
            write (u, "(1x,A)")  "weight (prc) = [undefined]"
         end if
         if (event%excess_prc_known) then
            write (u, "(1x,A," // fmt // ")")  "excess (prc) = ", &
                 event%excess_prc
         else
            write (u, "(1x,A)")  "excess (prc) = [undefined]"
         end if
         do i = 1, event%n_alt
            if (event%sqme_ref_known) then
               write (u, "(1x,A,I0,A," // fmt // ")")  "sqme (", i, ")     = ",&
                    event%sqme_prc
            else
               write (u, "(1x,A,I0,A)")  "sqme (", i, ")    = [undefined]"
            end if
            if (event%weight_prc_known) then
               write (u, "(1x,A,I0,A," // fmt // ")")  "weight (", i, ")   = ",&
                    event%weight_prc
            else
               write (u, "(1x,A,I0,A)")  "weight (", i, ")   = [undefined]"
            end if
         end do
      end if
      call write_separator (u)
      if (event%particle_set_is_valid) then
         call event%particle_set%write (unit = u, &
              summary = eio%summary, compressed = eio%compressed, &
              testflag = eio%pacify)
      else
         write (u, "(1x,A)")  "Particle set: [invalid]"
      end if
    end subroutine dump
  end subroutine eio_dump_output

  module subroutine eio_dump_input_i_prc (eio, i_prc, iostat)
    class(eio_dump_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    call msg_bug ("Dump stream: event input not supported")
    i_prc = 0
    iostat = 1
  end subroutine eio_dump_input_i_prc

  module subroutine eio_dump_input_event &
       (eio, event, iostat, event_handle)
    class(eio_dump_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    call msg_bug ("Dump stream: event input not supported")
    iostat = 1
  end subroutine eio_dump_input_event

  module subroutine eio_dump_skip (eio, iostat)
    class(eio_dump_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    iostat = 0
  end subroutine eio_dump_skip


end submodule eio_dump_s

