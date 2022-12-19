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

submodule (eio_ascii) eio_ascii_s

  use io_units
  use diagnostics
  use hep_common
  use hep_events

  implicit none

contains

  module subroutine eio_ascii_set_parameters (eio, &
       keep_beams, keep_remnants, ensure_order, extension, &
       show_process, show_transforms, show_decay, verbose)
    class(eio_ascii_t), intent(inout) :: eio
    logical, intent(in), optional :: keep_beams
    logical, intent(in), optional :: keep_remnants
    logical, intent(in), optional :: ensure_order
    type(string_t), intent(in), optional :: extension
    logical, intent(in), optional :: show_process, show_transforms, show_decay
    logical, intent(in), optional :: verbose
    if (present (keep_beams))  eio%keep_beams = keep_beams
    if (present (keep_remnants))  eio%keep_remnants = keep_remnants
    if (present (ensure_order))  eio%ensure_order = ensure_order
    if (present (extension)) then
       eio%extension = extension
    else
       select type (eio)
       type is (eio_ascii_ascii_t)
          eio%extension = "evt"
       type is (eio_ascii_athena_t)
          eio%extension = "athena.evt"
       type is (eio_ascii_debug_t)
          eio%extension = "debug"
       type is (eio_ascii_hepevt_t)
          eio%extension = "hepevt"
       type is (eio_ascii_hepevt_verb_t)
          eio%extension = "hepevt.verb"
       type is (eio_ascii_lha_t)
          eio%extension = "lha"
       type is (eio_ascii_lha_verb_t)
          eio%extension = "lha.verb"
       type is (eio_ascii_long_t)
          eio%extension = "long.evt"
       type is (eio_ascii_mokka_t)
          eio%extension = "mokka.evt"
       type is (eio_ascii_short_t)
          eio%extension = "short.evt"
       end select
    end if
    select type (eio)
    type is (eio_ascii_debug_t)
       if (present (show_process))  eio%show_process = show_process
       if (present (show_transforms))  eio%show_transforms = show_transforms
       if (present (show_decay))  eio%show_decay = show_decay
       if (present (verbose))  eio%verbose = verbose
    end select
  end subroutine eio_ascii_set_parameters

  module subroutine eio_ascii_write (object, unit)
    class(eio_ascii_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    select type (object)
    type is (eio_ascii_ascii_t)
       write (u, "(1x,A)")  "ASCII event stream (default format):"
    type is (eio_ascii_athena_t)
       write (u, "(1x,A)")  "ASCII event stream (ATHENA format):"
    type is (eio_ascii_debug_t)
       write (u, "(1x,A)")  "ASCII event stream (Debugging format):"
    type is (eio_ascii_hepevt_t)
       write (u, "(1x,A)")  "ASCII event stream (HEPEVT format):"
    type is (eio_ascii_hepevt_verb_t)
       write (u, "(1x,A)")  "ASCII event stream (verbose HEPEVT format):"
    type is (eio_ascii_lha_t)
       write (u, "(1x,A)")  "ASCII event stream (LHA format):"
    type is (eio_ascii_lha_verb_t)
       write (u, "(1x,A)")  "ASCII event stream (verbose LHA format):"
    type is (eio_ascii_long_t)
       write (u, "(1x,A)")  "ASCII event stream (long format):"
    type is (eio_ascii_mokka_t)
       write (u, "(1x,A)")  "ASCII event stream (MOKKA format):"
    type is (eio_ascii_short_t)
       write (u, "(1x,A)")  "ASCII event stream (short format):"
    end select
    if (object%writing) then
       write (u, "(3x,A,A)")  "Writing to file   = ", char (object%filename)
    else
       write (u, "(3x,A)")  "[closed]"
    end if
    write (u, "(3x,A,L1)")    "Keep beams        = ", object%keep_beams
    write (u, "(3x,A,L1)")    "Keep remnants     = ", object%keep_remnants
    select type (object)
    type is (eio_ascii_debug_t)
       write (u, "(3x,A,L1)")    "Show process      = ", object%show_process
       write (u, "(3x,A,L1)")    "Show transforms   = ", object%show_transforms
       write (u, "(3x,A,L1)")    "Show decay tree   = ", object%show_decay
       write (u, "(3x,A,L1)")    "Verbose output    = ", object%verbose
    end select
  end subroutine eio_ascii_write

  module subroutine eio_ascii_final (object)
    class(eio_ascii_t), intent(inout) :: object
    if (object%writing) then
       write (msg_buffer, "(A,A,A)")  "Events: closing ASCII file '", &
            char (object%filename), "'"
       call msg_message ()
       close (object%unit)
       object%writing = .false.
    end if
  end subroutine eio_ascii_final

  module subroutine eio_ascii_init_out &
       (eio, sample, data, success, extension)
    class(eio_ascii_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    logical :: is_width
    integer :: i
    if (.not. present (data)) &
         call msg_bug ("ASCII initialization: missing data")
    is_width = data%n_beam == 1
    eio%sample = sample
    call eio%check_normalization (data)
    call eio%set_splitting (data)
    call eio%set_filename ()
    eio%unit = free_unit ()
    write (msg_buffer, "(A,A,A)")  "Events: writing to ASCII file '", &
         char (eio%filename), "'"
    call msg_message ()
    eio%writing = .true.
    open (eio%unit, file = char (eio%filename), &
         action = "write", status = "replace")
    select type (eio)
    type is (eio_ascii_lha_t)
       call heprup_init( &
            data%pdg_beam, &
            data%energy_beam, &
            n_processes = data%n_proc, &
            unweighted = data%unweighted, &
            negative_weights = data%negative_weights)
       do i = 1, data%n_proc
          call heprup_set_process_parameters (i = i, &
               process_id = data%proc_num_id(i), &
               cross_section = data%cross_section(i), &
               error = data%error(i), &
               is_width = is_width)
       end do
       call heprup_write_ascii (eio%unit)
    type is (eio_ascii_lha_verb_t)
       call heprup_init( &
            data%pdg_beam, &
            data%energy_beam, &
            n_processes = data%n_proc, &
            unweighted = data%unweighted, &
            negative_weights = data%negative_weights)
       do i = 1, data%n_proc
          call heprup_set_process_parameters (i = i, &
               process_id = data%proc_num_id(i), &
               cross_section = data%cross_section(i), &
               error = data%error(i), &
               is_width = is_width)
       end do
       call heprup_write_verbose (eio%unit)
    end select
    if (present (success))  success = .true.
  end subroutine eio_ascii_init_out

  module subroutine eio_ascii_check_normalization (eio, data)
    class(eio_ascii_t), intent(in) :: eio
    type(event_sample_data_t), intent(in) :: data
    if (data%unweighted) then
    else
       select type (eio)
       type is (eio_ascii_athena_t);  call msg_fatal &
            ("Event output (Athena format): events must be unweighted.")
       type is (eio_ascii_hepevt_t);  call msg_fatal &
            ("Event output (HEPEVT format): events must be unweighted.")
       type is (eio_ascii_hepevt_verb_t);  call msg_fatal &
            ("Event output (HEPEVT format): events must be unweighted.")
       end select
       select case (data%norm_mode)
       case (NORM_SIGMA)
       case default
          select type (eio)
          type is (eio_ascii_lha_t)
             call msg_fatal &
                  ("Event output (LHA): normalization for weighted events &
                  &must be 'sigma'")
          type is (eio_ascii_lha_verb_t)
             call msg_fatal &
                  ("Event output (LHA): normalization for weighted events &
                  &must be 'sigma'")
          end select
       end select
    end if
  end subroutine eio_ascii_check_normalization

  module subroutine eio_ascii_init_in &
       (eio, sample, data, success, extension)
    class(eio_ascii_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    call msg_bug ("ASCII: event input not supported")
    if (present (success))  success = .false.
  end subroutine eio_ascii_init_in

  module subroutine eio_ascii_switch_inout (eio, success)
    class(eio_ascii_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    call msg_bug ("ASCII: in-out switch not supported")
    if (present (success))  success = .false.
  end subroutine eio_ascii_switch_inout

  module subroutine eio_ascii_split_out (eio)
    class(eio_ascii_t), intent(inout) :: eio
    if (eio%split) then
       eio%split_index = eio%split_index + 1
       call eio%set_filename ()
       write (msg_buffer, "(A,A,A)")  "Events: writing to ASCII file '", &
            char (eio%filename), "'"
       call msg_message ()
       close (eio%unit)
       open (eio%unit, file = char (eio%filename), &
            action = "write", status = "replace")
       select type (eio)
       type is (eio_ascii_lha_t)
          call heprup_write_ascii (eio%unit)
       type is (eio_ascii_lha_verb_t)
          call heprup_write_verbose (eio%unit)
       end select
    end if
  end subroutine eio_ascii_split_out

  module subroutine eio_ascii_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_ascii_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    integer, intent(in) :: i_prc
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    if (present (passed)) then
       if (.not. passed) then
          select type (eio)
          type is (eio_ascii_debug_t)
          type is (eio_ascii_ascii_t)
          class default
             return
          end select
       end if
    end if
    if (eio%writing) then
       select type (eio)
       type is (eio_ascii_lha_t)
          call hepeup_from_event (event, &
               process_index = i_prc, &
               keep_beams = eio%keep_beams, &
               keep_remnants = eio%keep_remnants)
          call hepeup_write_lha (eio%unit)
       type is (eio_ascii_lha_verb_t)
          call hepeup_from_event (event, &
               process_index = i_prc, &
               keep_beams = eio%keep_beams, &
               keep_remnants = eio%keep_remnants)
          call hepeup_write_verbose (eio%unit)
       type is (eio_ascii_ascii_t)
          call event%write (eio%unit, &
               show_process = .false., &
               show_transforms = .false., &
               show_decay = .false., &
               verbose = .false., testflag = pacify)
       type is (eio_ascii_athena_t)
          call hepevt_from_event (event, &
               keep_beams = eio%keep_beams, &
               keep_remnants = eio%keep_remnants, &
               ensure_order = eio%ensure_order)
          call hepevt_write_athena (eio%unit)
       type is (eio_ascii_debug_t)
          call event%write (eio%unit, &
               show_process = eio%show_process, &
               show_transforms = eio%show_transforms, &
               show_decay = eio%show_decay, &
               verbose = eio%verbose, &
               testflag = pacify)
       type is (eio_ascii_hepevt_t)
          call hepevt_from_event (event, &
               keep_beams = eio%keep_beams, &
               keep_remnants = eio%keep_remnants, &
               ensure_order = eio%ensure_order)
          call hepevt_write_hepevt (eio%unit)
       type is (eio_ascii_hepevt_verb_t)
          call hepevt_from_event (event, &
               keep_beams = eio%keep_beams, &
               keep_remnants = eio%keep_remnants, &
               ensure_order = eio%ensure_order)
          call hepevt_write_verbose (eio%unit)
       type is (eio_ascii_long_t)
          call hepevt_from_event (event, &
               keep_beams = eio%keep_beams, &
               keep_remnants = eio%keep_remnants, &
               ensure_order = eio%ensure_order)
          call hepevt_write_ascii (eio%unit, .true.)
       type is (eio_ascii_mokka_t)
          call hepevt_from_event (event, &
               keep_beams = eio%keep_beams, &
               keep_remnants = eio%keep_remnants, &
               ensure_order = eio%ensure_order)
          call hepevt_write_mokka (eio%unit)
       type is (eio_ascii_short_t)
          call hepevt_from_event (event, &
               keep_beams = eio%keep_beams, &
               keep_remnants = eio%keep_remnants, &
               ensure_order = eio%ensure_order)
          call hepevt_write_ascii (eio%unit, .false.)
       end select
    else
       call eio%write ()
       call msg_fatal ("ASCII file is not open for writing")
    end if
  end subroutine eio_ascii_output

  module subroutine eio_ascii_input_i_prc (eio, i_prc, iostat)
    class(eio_ascii_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    call msg_bug ("ASCII: event input not supported")
    i_prc = 0
    iostat = 1
  end subroutine eio_ascii_input_i_prc

  module subroutine eio_ascii_input_event &
       (eio, event, iostat, event_handle)
    class(eio_ascii_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    call msg_bug ("ASCII: event input not supported")
    iostat = 1
  end subroutine eio_ascii_input_event

  module subroutine eio_ascii_skip (eio, iostat)
     class(eio_ascii_t), intent(inout) :: eio
     integer, intent(out) :: iostat
     iostat = 0
  end subroutine eio_ascii_skip


end submodule eio_ascii_s

