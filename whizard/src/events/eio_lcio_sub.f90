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

submodule (eio_lcio) eio_lcio_s

  use io_units
  use string_utils
  use diagnostics
  use particles
  use hep_events

  implicit none

contains

  module subroutine eio_lcio_set_parameters &
       (eio, recover_beams, use_alphas_from_file, use_scale_from_file, &
       extension, proc_as_run_id, lcio_run_id)
    class(eio_lcio_t), intent(inout) :: eio
    logical, intent(in), optional :: recover_beams
    logical, intent(in), optional :: use_alphas_from_file
    logical, intent(in), optional :: use_scale_from_file
    logical, intent(in), optional :: proc_as_run_id
    integer, intent(in), optional :: lcio_run_id
    type(string_t), intent(in), optional :: extension
    if (present (recover_beams))  eio%recover_beams = recover_beams
    if (present (use_alphas_from_file)) &
         eio%use_alphas_from_file = use_alphas_from_file
    if (present (use_scale_from_file)) &
         eio%use_scale_from_file = use_scale_from_file
    if (present (proc_as_run_id)) &
         eio%proc_as_run_id = proc_as_run_id
    if (present (lcio_run_id)) &
         eio%lcio_run_id = lcio_run_id
    if (present (extension)) then
       eio%extension = extension
    else
       eio%extension = "slcio"
    end if
  end subroutine eio_lcio_set_parameters

  module subroutine eio_lcio_write (object, unit)
    class(eio_lcio_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "LCIO event stream:"
    if (object%writing) then
       write (u, "(3x,A,A)")  "Writing to file   = ", char (object%filename)
    else if (object%reading) then
       write (u, "(3x,A,A)")  "Reading from file = ", char (object%filename)
    else
       write (u, "(3x,A)")  "[closed]"
    end if
    write (u, "(3x,A,L1)")    "Recover beams     = ", object%recover_beams
    write (u, "(3x,A,L1)")    "Alpha_s from file = ", &
         object%use_alphas_from_file
    write (u, "(3x,A,L1)")    "Scale from file   = ", &
         object%use_scale_from_file
    write (u, "(3x,A,L1)")    "Process as run ID = ", &
         object%proc_as_run_id
    write (u, "(3x,A,I0)")    "LCIO run ID       = ", &
         object%lcio_run_id
    write (u, "(3x,A,A,A)")     "File extension    = '", &
         char (object%extension), "'"
    if (allocated (object%proc_num_id)) then
       write (u, "(3x,A)")  "Numerical process IDs:"
       do i = 1, size (object%proc_num_id)
          write (u, "(5x,I0,': ',I0)")  i, object%proc_num_id(i)
       end do
    end if
  end subroutine eio_lcio_write

  module subroutine eio_lcio_final (object)
    class(eio_lcio_t), intent(inout) :: object
    if (allocated (object%proc_num_id))  deallocate (object%proc_num_id)
    if (object%writing) then
       write (msg_buffer, "(A,A,A)")  "Events: closing LCIO file '", &
            char (object%filename), "'"
       call msg_message ()
       call lcio_writer_close (object%lcio_writer)
       object%writing = .false.
    else if (object%reading) then
       write (msg_buffer, "(A,A,A)")  "Events: closing LCIO file '", &
            char (object%filename), "'"
       call msg_message ()
       call lcio_reader_close (object%lcio_reader)
       object%reading = .false.
    end if
  end subroutine eio_lcio_final

  module subroutine eio_lcio_split_out (eio)
    class(eio_lcio_t), intent(inout) :: eio
    if (eio%split) then
       eio%split_index = eio%split_index + 1
       call eio%set_filename ()
       write (msg_buffer, "(A,A,A)")  "Events: writing to LCIO file '", &
            char (eio%filename), "'"
       call msg_message ()
       call lcio_writer_close (eio%lcio_writer)
       call lcio_writer_open_out (eio%lcio_writer, eio%filename)
    end if
  end subroutine eio_lcio_split_out

  module subroutine eio_lcio_common_init (eio, sample, data, extension)
    class(eio_lcio_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    if (.not. present (data)) &
         call msg_bug ("LCIO initialization: missing data")
    eio%data = data
    if (data%unweighted) then
       select case (data%norm_mode)
       case (NORM_UNIT)
       case default; call msg_fatal &
            ("LCIO: normalization for unweighted events must be '1'")
       end select
    else
       call msg_fatal ("LCIO: events must be unweighted")
    end if
    eio%n_alt = data%n_alt
    eio%sample = sample
    if (present (extension)) then
       eio%extension = extension
    end if
    call eio%set_filename ()
    allocate (eio%proc_num_id (data%n_proc), source = data%proc_num_id)
  end subroutine eio_lcio_common_init

  module subroutine eio_lcio_init_out &
       (eio, sample, data, success, extension)
    class(eio_lcio_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    call eio%set_splitting (data)
    call eio%common_init (sample, data, extension)
    write (msg_buffer, "(A,A,A)")  "Events: writing to LCIO file '", &
         char (eio%filename), "'"
    call msg_message ()
    eio%writing = .true.
    call lcio_writer_open_out (eio%lcio_writer, eio%filename)
    call lcio_run_header_init (eio%lcio_run_hdr)
    call lcio_run_header_write (eio%lcio_writer, eio%lcio_run_hdr)
    if (present (success))  success = .true.
  end subroutine eio_lcio_init_out

  module subroutine eio_lcio_init_in (eio, sample, data, success, extension)
    class(eio_lcio_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    logical :: exist
    eio%split = .false.
    call eio%common_init (sample, data, extension)
    write (msg_buffer, "(A,A,A)")  "Events: reading from LCIO file '", &
         char (eio%filename), "'"
    call msg_message ()
    inquire (file = char (eio%filename), exist = exist)
    if (.not. exist)  call msg_fatal ("Events: LCIO file not found.")
    eio%reading = .true.
    call lcio_open_file (eio%lcio_reader, eio%filename)
    if (present (success))  success = .true.
  end subroutine eio_lcio_init_in

  module subroutine eio_lcio_switch_inout (eio, success)
    class(eio_lcio_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    call msg_bug ("LCIO: in-out switch not supported")
    if (present (success))  success = .false.
  end subroutine eio_lcio_switch_inout

  module subroutine eio_lcio_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_lcio_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    integer, intent(in) :: i_prc
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    type(particle_set_t), pointer :: pset_ptr
    real(default) :: sqme_prc, weight
    real(default), dimension(:), allocatable :: pol
    integer :: i
    if (present (passed)) then
       if (.not. passed)  return
    end if
    if (eio%writing) then
       pset_ptr => event%get_particle_set_ptr ()
       if (eio%proc_as_run_id) then
          call lcio_event_init (eio%lcio_event, &
               proc_id = eio%proc_num_id (i_prc), &
               event_id = event%get_index (), &
               run_id = eio%proc_num_id (i_prc))
       else
          call lcio_event_init (eio%lcio_event, &
               proc_id = eio%proc_num_id (i_prc), &
               event_id = event%get_index (), &
               run_id = eio%lcio_run_id)
       end if
       call lcio_event_from_particle_set (eio%lcio_event, pset_ptr)
       call lcio_event_set_weight (eio%lcio_event, event%weight_prc)
       call lcio_event_set_sqrts (eio%lcio_event, event%get_sqrts ())
       call lcio_event_set_sqme (eio%lcio_event, event%get_sqme_prc ())
       call lcio_event_set_scale (eio%lcio_event, event%get_fac_scale ())
       call lcio_event_set_alpha_qcd (eio%lcio_event, event%get_alpha_s ())
       if (eio%data%n_beam == 2) then
       call lcio_event_set_xsec (eio%lcio_event, eio%data%cross_section(i_prc), &
            eio%data%error(i_prc))
       end if
       pol = event%get_polarization ()
       do i = 1, eio%data%n_beam
          call lcio_event_set_polarization (eio%lcio_event, pol(i), i)
       end do
       call lcio_event_set_beam_file (eio%lcio_event, &
            event%get_beam_file ())
       call lcio_event_set_process_name (eio%lcio_event, &
            event%get_process_name ())
       do i = 1, eio%n_alt
           sqme_prc = event%get_sqme_alt(i)
           weight = event%get_weight_alt(i)
           call lcio_event_set_alt_sqme (eio%lcio_event, sqme_prc, i)
           call lcio_event_set_alt_weight (eio%lcio_event, weight, i)
       end do
       call lcio_event_write (eio%lcio_writer, eio%lcio_event)
       call maybe_transfer_event_to_handle (eio%lcio_event, &
            event_handle, .true.)
       call lcio_event_final (eio%lcio_event, .true.)
    else
       call eio%write ()
       call msg_fatal ("LCIO file is not open for writing")
    end if
  end subroutine eio_lcio_output

  module subroutine eio_lcio_input_i_prc (eio, i_prc, iostat)
    class(eio_lcio_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    logical :: ok
    integer :: i, proc_num_id
    iostat = 0
    call lcio_read_event (eio%lcio_reader, eio%lcio_event, ok)
    if (.not. ok) then
       iostat = -1
       return
    end if
    proc_num_id = lcio_event_get_process_id (eio%lcio_event)
    i_prc = 0
    FIND_I_PRC: do i = 1, size (eio%proc_num_id)
       if (eio%proc_num_id(i) == proc_num_id) then
          i_prc = i
          exit FIND_I_PRC
       end if
    end do FIND_I_PRC
    if (i_prc == 0)  call err_index
  contains
    subroutine err_index
      call msg_error ("LCIO: reading events: undefined process ID " &
           // char (str (proc_num_id)) // ", aborting read")
      iostat = 1
    end subroutine err_index
  end subroutine eio_lcio_input_i_prc

  module subroutine eio_lcio_input_event (eio, event, iostat, event_handle)
    class(eio_lcio_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    iostat = 0
    call event%reset_contents ()
    call event%select (0, 0, 1)
    call event%set_index (lcio_event_get_event_index (eio%lcio_event))
    call lcio_to_event (event, eio%lcio_event, eio%fallback_model, &
         recover_beams = eio%recover_beams, &
         use_alpha_s = eio%use_alphas_from_file, &
         use_scale = eio%use_scale_from_file)
    call maybe_transfer_event_to_handle (eio%lcio_event, &
         event_handle, .false.)
    call lcio_event_final (eio%lcio_event, .false.)
  end subroutine eio_lcio_input_event

  subroutine maybe_transfer_event_to_handle (lcio_event, event_handle, delete)
    type(lcio_event_t), intent(inout) :: lcio_event
    class(event_handle_t), intent(inout), optional :: event_handle
    logical, intent(in) :: delete
    if (present (event_handle)) then
       select type (event_handle)
       type is (lcio_event_t)
          call lcio_event_final (event_handle, delete)  ! just in case
          event_handle = lcio_event
          call lcio_event_nullify (lcio_event) ! avoid destructor call
       end select
    end if
  end subroutine maybe_transfer_event_to_handle

  module subroutine eio_lcio_skip (eio, iostat)
    class(eio_lcio_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    iostat = 0
  end subroutine eio_lcio_skip


end submodule eio_lcio_s

