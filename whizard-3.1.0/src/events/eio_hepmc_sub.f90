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

submodule (eio_hepmc) eio_hepmc_s

  use io_units
  use string_utils
  use diagnostics
  use particles
  use model_data
  use hep_events

  implicit none

contains

  module subroutine eio_hepmc_set_parameters &
       (eio, recover_beams, use_alphas_from_file, &
       use_scale_from_file, extension, output_cross_section, &
       hepmc3_mode, hepmc3_write_flows)
    class(eio_hepmc_t), intent(inout) :: eio
    logical, intent(in), optional :: recover_beams
    logical, intent(in), optional :: use_alphas_from_file
    logical, intent(in), optional :: use_scale_from_file
    logical, intent(in), optional :: output_cross_section
    type(string_t), intent(in), optional :: extension
    integer, intent(in), optional :: hepmc3_mode
    logical ,intent(in), optional :: hepmc3_write_flows
    if (present (recover_beams)) &
         eio%recover_beams = recover_beams
    if (present (use_alphas_from_file)) &
         eio%use_alphas_from_file = use_alphas_from_file
    if (present (use_scale_from_file)) &
         eio%use_scale_from_file = use_scale_from_file
    if (present (extension)) then
       eio%extension = extension
    else
       eio%extension = "hepmc"
    end if
    if (present (output_cross_section)) &
         eio%output_cross_section = output_cross_section
    if (present (hepmc3_mode)) &
         eio%hepmc3_mode = hepmc3_mode
    if (present (hepmc3_write_flows)) &
         eio%hepmc3_flows = hepmc3_write_flows
  end subroutine eio_hepmc_set_parameters

  module subroutine eio_hepmc_write (object, unit)
    class(eio_hepmc_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "HepMC event stream:"
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
    write (u, "(3x,A,A,A)")     "File extension    = '", &
         char (object%extension), "'"
    write (u, "(3x,A,I0)")    "HepMC3 mode       = ", object%hepmc3_mode
    write (u, "(3x,A,L1)")    "HepMC3 flows      = ", object%hepmc3_flows
    if (allocated (object%proc_num_id)) then
       write (u, "(3x,A)")  "Numerical process IDs:"
       do i = 1, size (object%proc_num_id)
          write (u, "(5x,I0,': ',I0)")  i, object%proc_num_id(i)
       end do
    end if
  end subroutine eio_hepmc_write

  module subroutine eio_hepmc_final (object)
    class(eio_hepmc_t), intent(inout) :: object
    if (allocated (object%proc_num_id))  deallocate (object%proc_num_id)
    if (object%writing) then
       write (msg_buffer, "(A,A,A)")  "Events: closing HepMC file '", &
            char (object%filename), "'"
       call msg_message ()
       call hepmc_iostream_close (object%iostream)
       object%writing = .false.
    else if (object%reading) then
       write (msg_buffer, "(A,A,A)")  "Events: closing HepMC file '", &
            char (object%filename), "'"
       call msg_message ()
       call hepmc_iostream_close (object%iostream)
       object%reading = .false.
    end if
  end subroutine eio_hepmc_final

  module subroutine eio_hepmc_split_out (eio)
    class(eio_hepmc_t), intent(inout) :: eio
    if (eio%split) then
       eio%split_index = eio%split_index + 1
       call eio%set_filename ()
       write (msg_buffer, "(A,A,A)")  "Events: writing to HepMC file '", &
            char (eio%filename), "'"
       call msg_message ()
       call hepmc_iostream_close (eio%iostream)
       call hepmc_iostream_open_out (eio%iostream, &
            eio%filename, eio%hepmc3_mode)
    end if
  end subroutine eio_hepmc_split_out

  module subroutine eio_hepmc_common_init (eio, sample, data, extension)
    class(eio_hepmc_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    if (.not. present (data)) &
         call msg_bug ("HepMC initialization: missing data")
    eio%data = data
    ! We could relax this condition now with weighted hepmc events
    if (data%unweighted) then
       select case (data%norm_mode)
       case (NORM_UNIT)
       case default; call msg_fatal &
            ("HepMC: normalization for unweighted events must be '1'")
       end select
    end if
    eio%sample = sample
    if (present (extension)) then
       eio%extension = extension
    end if
    call eio%set_filename ()
    allocate (eio%proc_num_id (data%n_proc), source = data%proc_num_id)
  end subroutine eio_hepmc_common_init

  module subroutine eio_hepmc_init_out (eio, sample, data, success, extension)
    class(eio_hepmc_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    call eio%set_splitting (data)
    call eio%common_init (sample, data, extension)
    write (msg_buffer, "(A,A,A)")  "Events: writing to HepMC file '", &
         char (eio%filename), "'"
    call msg_message ()
    eio%writing = .true.
    call hepmc_iostream_open_out (eio%iostream, &
         eio%filename, eio%hepmc3_mode)
    if (present (success))  success = .true.
  end subroutine eio_hepmc_init_out

  module subroutine eio_hepmc_init_in (eio, sample, data, success, extension)
    class(eio_hepmc_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    logical :: exist
    eio%split = .false.
    call eio%common_init (sample, data, extension)
    write (msg_buffer, "(A,A,A)")  "Events: reading from HepMC file '", &
         char (eio%filename), "'"
    call msg_message ()
    inquire (file = char (eio%filename), exist = exist)
    if (.not. exist)  call msg_fatal ("Events: HepMC file not found.")
    eio%reading = .true.
    call hepmc_iostream_open_in (eio%iostream, &
         eio%filename, eio%hepmc3_mode)
    if (present (success))  success = .true.
  end subroutine eio_hepmc_init_in

  module subroutine eio_hepmc_switch_inout (eio, success)
    class(eio_hepmc_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    call msg_bug ("HepMC: in-out switch not supported")
    if (present (success))  success = .false.
  end subroutine eio_hepmc_switch_inout

  module subroutine eio_hepmc_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_hepmc_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    integer, intent(in) :: i_prc
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    type(particle_set_t), pointer :: pset_ptr
    if (present (passed)) then
       if (.not. passed)  return
    end if
    if (eio%writing) then
       pset_ptr => event%get_particle_set_ptr ()
       call hepmc_event_init (eio%hepmc_event, &
            proc_id = eio%proc_num_id(i_prc), &
            event_id = event%get_index ())
       if (eio%output_cross_section .and. eio%data%n_beam == 2) then
          call hepmc_event_from_particle_set (eio%hepmc_event, pset_ptr, &
               eio%data%cross_section(i_prc), eio%data%error(i_prc), &
               color = eio%hepmc3_flows)
       else
          call hepmc_event_from_particle_set (eio%hepmc_event, pset_ptr, &
               color = eio%hepmc3_flows)
       end if
       call hepmc_event_set_scale (eio%hepmc_event, event%get_fac_scale ())
       call hepmc_event_set_alpha_qcd (eio%hepmc_event, event%get_alpha_s ())
       call hepmc_event_set_alpha_qed (eio%hepmc_event, -1._default)
       if (.not. eio%data%unweighted .or. eio%data%negative_weights) then
          select case (eio%data%norm_mode)
          case (NORM_UNIT,NORM_N_EVT)
             call hepmc_event_add_weight &
                  (eio%hepmc_event, event%weight_prc, .false.)
          case default
             call hepmc_event_add_weight &
                  (eio%hepmc_event, event%weight_prc, .true.)
          end select
       end if
       call hepmc_iostream_write_event (eio%iostream, &
            eio%hepmc_event, eio%hepmc3_mode)
       call maybe_transfer_event_to_handle (eio%hepmc_event, event_handle)
       call hepmc_event_final (eio%hepmc_event)
    else
       call eio%write ()
       call msg_fatal ("HepMC file is not open for writing")
    end if
  end subroutine eio_hepmc_output

  module subroutine eio_hepmc_input_i_prc (eio, i_prc, iostat)
    class(eio_hepmc_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    logical :: ok
    integer :: i, proc_num_id
    iostat = 0
    call hepmc_event_init (eio%hepmc_event)
    call hepmc_iostream_read_event (eio%iostream, &
         eio%hepmc_event, ok=ok)
    proc_num_id = hepmc_event_get_process_id (eio%hepmc_event)
    if (.not. ok) then
       iostat = -1
       return
    end if
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
      call msg_error ("HepMC: reading events: undefined process ID " &
           // char (str (proc_num_id)) // ", aborting read")
      iostat = 1
    end subroutine err_index
  end subroutine eio_hepmc_input_i_prc

  module subroutine eio_hepmc_input_event (eio, event, iostat, event_handle)
    class(eio_hepmc_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    iostat = 0
    call event%reset_contents ()
    call event%select (1, 1, 1)
    call hepmc_to_event (event, eio%hepmc_event, &
         eio%fallback_model, &
         recover_beams = eio%recover_beams, &
         use_alpha_s = eio%use_alphas_from_file, &
         use_scale = eio%use_scale_from_file)
    call maybe_transfer_event_to_handle (eio%hepmc_event, event_handle)
    call hepmc_event_final (eio%hepmc_event)
  end subroutine eio_hepmc_input_event

  subroutine maybe_transfer_event_to_handle (hepmc_event, event_handle)
    type(hepmc_event_t), intent(inout) :: hepmc_event
    class(event_handle_t), intent(inout), optional :: event_handle
    if (present (event_handle)) then
       select type (event_handle)
       type is (hepmc_event_t)
          call hepmc_event_final (event_handle)  ! just in case
          event_handle = hepmc_event
          call hepmc_event_nullify (hepmc_event) ! avoid destructor call
       end select
    end if
  end subroutine maybe_transfer_event_to_handle

  module subroutine eio_hepmc_skip (eio, iostat)
    class(eio_hepmc_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    iostat = 0
  end subroutine eio_hepmc_skip


end submodule eio_hepmc_s

