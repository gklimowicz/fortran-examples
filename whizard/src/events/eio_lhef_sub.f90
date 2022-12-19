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

submodule (eio_lhef) eio_lhef_s

  use io_units
  use string_utils
  use numeric_utils
  use diagnostics
  use os_interface
  use hep_common
  use hep_events

  implicit none

contains

  module subroutine eio_lhef_set_parameters (eio, &
       keep_beams, keep_remnants, recover_beams, &
       use_alphas_from_file, use_scale_from_file, &
       version, extension, write_sqme_ref, write_sqme_prc, write_sqme_alt)
    class(eio_lhef_t), intent(inout) :: eio
    logical, intent(in), optional :: keep_beams
    logical, intent(in), optional :: keep_remnants
    logical, intent(in), optional :: recover_beams
    logical, intent(in), optional :: use_alphas_from_file
    logical, intent(in), optional :: use_scale_from_file
    character(*), intent(in), optional :: version
    type(string_t), intent(in), optional :: extension
    logical, intent(in), optional :: write_sqme_ref
    logical, intent(in), optional :: write_sqme_prc
    logical, intent(in), optional :: write_sqme_alt
    if (present (keep_beams))  eio%keep_beams = keep_beams
    if (present (keep_remnants))  eio%keep_remnants = keep_remnants
    if (present (recover_beams))  eio%recover_beams = recover_beams
    if (present (use_alphas_from_file)) &
         eio%use_alphas_from_file = use_alphas_from_file
    if (present (use_scale_from_file)) &
         eio%use_scale_from_file = use_scale_from_file
    if (present (version)) then
       select case (version)
       case ("1.0", "2.0", "3.0")
          eio%version = version
       case default
          call msg_error ("LHEF version " // version &
               // " is not supported.  Inserting 2.0")
          eio%version = "2.0"
       end select
    end if
    if (present (extension)) then
       eio%extension = extension
    else
       eio%extension = "lhe"
    end if
    if (present (write_sqme_ref))  eio%write_sqme_ref = write_sqme_ref
    if (present (write_sqme_prc))  eio%write_sqme_prc = write_sqme_prc
    if (present (write_sqme_alt))  eio%write_sqme_alt = write_sqme_alt
  end subroutine eio_lhef_set_parameters

  module subroutine eio_lhef_write (object, unit)
    class(eio_lhef_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "LHEF event stream:"
    if (object%writing) then
       write (u, "(3x,A,A)")  "Writing to file   = ", char (object%filename)
    else if (object%reading) then
       write (u, "(3x,A,A)")  "Reading from file = ", char (object%filename)
    else
       write (u, "(3x,A)")  "[closed]"
    end if
    write (u, "(3x,A,L1)")    "Keep beams        = ", object%keep_beams
    write (u, "(3x,A,L1)")    "Keep remnants     = ", object%keep_remnants
    write (u, "(3x,A,L1)")    "Recover beams     = ", object%recover_beams
    write (u, "(3x,A,L1)")    "Alpha_s from file = ", &
         object%use_alphas_from_file
    write (u, "(3x,A,L1)")    "Scale from file   = ", &
         object%use_scale_from_file
    write (u, "(3x,A,A)")     "Version           = ", object%version
    write (u, "(3x,A,A,A)")     "File extension    = '", &
         char (object%extension), "'"
    if (allocated (object%proc_num_id)) then
       write (u, "(3x,A)")  "Numerical process IDs:"
       do i = 1, size (object%proc_num_id)
          write (u, "(5x,I0,': ',I0)")  i, object%proc_num_id(i)
       end do
    end if
  end subroutine eio_lhef_write

  module subroutine eio_lhef_final (object)
    class(eio_lhef_t), intent(inout) :: object
    if (allocated (object%proc_num_id))  deallocate (object%proc_num_id)
    if (object%writing) then
       write (msg_buffer, "(A,A,A)")  "Events: closing LHEF file '", &
            char (object%filename), "'"
       call msg_message ()
       call object%write_footer ()
       close (object%unit)
       object%writing = .false.
    else if (object%reading) then
       write (msg_buffer, "(A,A,A)")  "Events: closing LHEF file '", &
            char (object%filename), "'"
       call msg_message ()
       call object%cstream%final ()
       close (object%unit)
       object%reading = .false.
    end if
  end subroutine eio_lhef_final

  module subroutine eio_lhef_common_init (eio, sample, data, extension)
    class(eio_lhef_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    if (.not. present (data)) &
         call msg_bug ("LHEF initialization: missing data")
    eio%data = data
    eio%unweighted = data%unweighted
    if (eio%unweighted) then
       select case (data%norm_mode)
       case (NORM_UNIT)
       case default;  call msg_fatal &
            ("LHEF: normalization for unweighted events must be '1'")
       end select
    else
       select case (data%norm_mode)
       case (NORM_SIGMA)
       case default;  call msg_fatal &
            ("LHEF: normalization for weighted events must be 'sigma'")
       end select
    end if
    eio%n_alt = data%n_alt
    eio%sample = sample
    if (present (extension)) then
       eio%extension = extension
    end if
    call eio%set_filename ()
    eio%unit = free_unit ()
    call eio%init_tags (data)
    allocate (eio%proc_num_id (data%n_proc), source = data%proc_num_id)
  end subroutine eio_lhef_common_init

  module subroutine eio_lhef_init_tags (eio, data)
    class(eio_lhef_t), intent(inout) :: eio
    type(event_sample_data_t), intent(in) :: data
    real(default), parameter :: pb_per_fb = 1.e-3_default
    real(default) :: xsec_width
    integer :: i
    call eio%tag_lhef%init ( &
         var_str ("LesHouchesEvents"), &
         [xml_attribute (var_str ("version"), var_str (eio%version))], &
         .true.)
    call eio%tag_head%init ( &
         var_str ("header"), &
         .true.)
    call eio%tag_init%init ( &
         var_str ("init"), &
         .true.)
    call eio%tag_event%init (var_str ("event"), &
         .true.)
    allocate (eio%tag_whiz_info)
    call eio%tag_whiz_info%init (var_str ("WhizardInfo"), .true.)
    select case (eio%version)
    case ("1.0")
       allocate (eio%tag_gen_n)
       call eio%tag_gen_n%init ( &
            var_str ("generator_name"), &
            .true.)
       allocate (eio%tag_gen_v)
       call eio%tag_gen_v%init ( &
            var_str ("generator_version"), &
            .true.)
    end select
    select case (eio%version)
    case ("2.0", "3.0")
       allocate (eio%tag_generator)
       call eio%tag_generator%init ( &
            var_str ("generator"), &
            [xml_attribute (var_str ("version"), var_str ("3.1.0"))], &
            .true.)
       allocate (eio%tag_xsecinfo)
       if (data%n_beam == 2) then
          xsec_width = data%total_cross_section * pb_per_fb
       else
          xsec_width = data%total_cross_section
       end if
       call eio%tag_xsecinfo%init ( &
            var_str ("xsecinfo"), &
            [xml_attribute (var_str ("neve"), str (data%n_evt)), &
             xml_attribute (var_str ("totxsec"), &
                            str (xsec_width))])
    end select
    select case (eio%version)
    case ("2.0")
       allocate (eio%tag_weight)
       call eio%tag_weight%init (var_str ("weight"), &
            [xml_attribute (var_str ("name"))])
       if (eio%write_sqme_ref) then
          allocate (eio%tag_sqme_ref)
          call eio%tag_sqme_ref%init (var_str ("weight"), &
               [xml_attribute (var_str ("name"), var_str ("sqme_ref"))], &
               .true.)
       end if
       if (eio%write_sqme_prc) then
          allocate (eio%tag_sqme_prc)
          call eio%tag_sqme_prc%init (var_str ("weight"), &
               [xml_attribute (var_str ("name"), var_str ("sqme_prc"))], &
               .true.)
       end if
       if (eio%n_alt > 0) then
          if (eio%write_sqme_alt) then
             allocate (eio%tag_sqme_alt (1))
             call eio%tag_sqme_alt(1)%init (var_str ("weight"), &
                  [xml_attribute (var_str ("name"), var_str ("sqme_alt"))], &
                  .true.)
          end if
          allocate (eio%tag_wgts_alt (1))
          call eio%tag_wgts_alt(1)%init (var_str ("weight"), &
               [xml_attribute (var_str ("name"), var_str ("wgts_alt"))], &
               .true.)
       end if
    case ("3.0")
       if (eio%write_sqme_ref) then
          allocate (eio%tag_sqme_ref)
          call eio%tag_sqme_ref%init (var_str ("weightinfo"), &
               [xml_attribute (var_str ("name"), var_str ("sqme_ref"))])
       end if
       if (eio%write_sqme_prc) then
          allocate (eio%tag_sqme_prc)
          call eio%tag_sqme_prc%init (var_str ("weightinfo"), &
               [xml_attribute (var_str ("name"), var_str ("sqme_prc"))])
       end if
       if (eio%n_alt > 0) then
          if (eio%write_sqme_alt) then
             allocate (eio%tag_sqme_alt (eio%n_alt))
             do i = 1, eio%n_alt
                call eio%tag_sqme_alt(i)%init (var_str ("weightinfo"), &
                     [xml_attribute (var_str ("name"), &
                                     var_str ("sqme_alt") // str (i))])
             end do
          end if
          allocate (eio%tag_wgts_alt (eio%n_alt))
          do i = 1, eio%n_alt
             call eio%tag_wgts_alt(i)%init (var_str ("weightinfo"), &
                  [xml_attribute (var_str ("name"), &
                                  var_str ("wgts_alt") // str (i))])
          end do
       end if
       allocate (eio%tag_weightinfo)
       call eio%tag_weightinfo%init (var_str ("weightinfo"), &
            [xml_attribute (var_str ("name"))])
       allocate (eio%tag_weights)
       call eio%tag_weights%init (var_str ("weights"), .true.)
    end select
  end subroutine eio_lhef_init_tags

  module subroutine eio_lhef_init_out (eio, sample, data, success, extension)
    class(eio_lhef_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    integer :: u, i
    logical :: is_width
    is_width = data%n_beam == 1
    call eio%set_splitting (data)
    call eio%common_init (sample, data, extension)
    write (msg_buffer, "(A,A,A)")  "Events: writing to LHEF file '", &
         char (eio%filename), "'"
    call msg_message ()
    eio%writing = .true.
    u = eio%unit
    open (u, file = char (eio%filename), &
         action = "write", status = "replace")
    call eio%write_header (is_width)
    call heprup_init &
         (data%pdg_beam, &
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
    call eio%tag_init%write (u);  write (u, *)
    call heprup_write_lhef (u)
    select case (eio%version)
    case ("2.0");  call eio%write_init_20 (data)
    case ("3.0");  call eio%write_init_30 (data)
    end select
    call eio%tag_init%close (u);  write (u, *)
    if (present (success))  success = .true.
  end subroutine eio_lhef_init_out

  module subroutine eio_lhef_init_in (eio, sample, data, success, extension)
    class(eio_lhef_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(string_t), intent(in), optional :: extension
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    logical :: exist, ok, closing
    type(event_sample_data_t) :: data_file
    type(string_t) :: string
    integer :: u
    eio%split = .false.
    call eio%common_init (sample, data, extension)
    write (msg_buffer, "(A,A,A)")  "Events: reading from LHEF file '", &
         char (eio%filename), "'"
    call msg_message ()
    inquire (file = char (eio%filename), exist = exist)
    if (.not. exist)  call msg_fatal ("Events: LHEF file not found.")
    eio%reading = .true.
    u = eio%unit
    open (u, file = char (eio%filename), &
         action = "read", status = "old")
    call eio%cstream%init (u)
    call eio%read_header ()
    call eio%tag_init%read (eio%cstream, ok)
    if (.not. ok)  call err_init
    select case (eio%version)
    case ("1.0");  call eio%read_init_10 (data_file)
       call eio%tag_init%read_content (eio%cstream, string, closing)
       if (string /= "" .or. .not. closing)  call err_init
    case ("2.0");  call eio%read_init_20 (data_file)
    case ("3.0");  call eio%read_init_30 (data_file)
    end select
    call eio%merge_data (data, data_file)
    if (present (success))  success = .true.

  contains

    subroutine err_init
      call msg_fatal ("LHEF: syntax error in init tag")
    end subroutine err_init

  end subroutine eio_lhef_init_in

  module subroutine eio_merge_data (eio, data, data_file)
    class(eio_lhef_t), intent(inout) :: eio
    type(event_sample_data_t), intent(inout) :: data
    type(event_sample_data_t), intent(in) :: data_file
    real, parameter :: tolerance = 1000 * epsilon (1._default)
    if (data%unweighted .neqv. data_file%unweighted)  call err_weights
    if (data%negative_weights .neqv. data_file%negative_weights) &
         call err_weights
    if (data%norm_mode /= data_file%norm_mode)  call err_norm
    if (data%n_beam /= data_file%n_beam)  call err_beams
    if (any (data%pdg_beam /= data_file%pdg_beam))  call err_beams
    if (any (abs ((data%energy_beam - data_file%energy_beam)) &
         > (data%energy_beam + data_file%energy_beam) * tolerance)) &
         call err_beams
    if (data%n_proc /= data_file%n_proc)  call err_proc
    if (any (data%proc_num_id /= data_file%proc_num_id))  call err_proc
    where (data%cross_section == 0)
       data%cross_section = data_file%cross_section
       data%error = data_file%error
    end where
    data%total_cross_section = sum (data%cross_section)
    if (data_file%n_evt > 0) then
       if (data%n_evt > 0 .and. data_file%n_evt /= data%n_evt)  call err_n_evt
       data%n_evt = data_file%n_evt
    end if
  contains
    subroutine err_weights
      call msg_fatal ("LHEF: mismatch in event weight properties")
    end subroutine err_weights
    subroutine err_norm
      call msg_fatal ("LHEF: mismatch in event normalization")
    end subroutine err_norm
    subroutine err_beams
      call msg_fatal ("LHEF: mismatch in beam properties")
    end subroutine err_beams
    subroutine err_proc
      call msg_fatal ("LHEF: mismatch in process definitions")
    end subroutine err_proc
    subroutine err_n_evt
      call msg_error ("LHEF: mismatch in specified number of events (ignored)")
    end subroutine err_n_evt
  end subroutine eio_merge_data

  module subroutine eio_lhef_switch_inout (eio, success)
    class(eio_lhef_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    call msg_bug ("LHEF: in-out switch not supported")
    if (present (success))  success = .false.
  end subroutine eio_lhef_switch_inout

  module subroutine eio_lhef_split_out (eio)
    class(eio_lhef_t), intent(inout) :: eio
    integer :: u
    if (eio%split) then
       eio%split_index = eio%split_index + 1
       call eio%set_filename ()
       write (msg_buffer, "(A,A,A)")  "Events: writing to LHEF file '", &
            char (eio%filename), "'"
       call msg_message ()
       call eio%write_footer ()
       u = eio%unit
       close (u)
       open (u, file = char (eio%filename), &
            action = "write", status = "replace")
       call eio%write_header ()
       call eio%tag_init%write (u);  write (u, *)
       call heprup_write_lhef (u)
       select case (eio%version)
       case ("2.0");  call eio%write_init_20 (eio%data)
       case ("3.0");  call eio%write_init_30 (eio%data)
       end select
       call eio%tag_init%close (u);  write (u, *)
    end if
  end subroutine eio_lhef_split_out

  module subroutine eio_lhef_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_lhef_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    integer, intent(in) :: i_prc
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    integer :: u
    u = given_output_unit (eio%unit);  if (u < 0)  return
    if (present (passed)) then
       if (.not. passed)  return
    end if
    if (eio%writing) then
       call hepeup_from_event (event, &
            process_index = eio%proc_num_id (i_prc), &
            keep_beams = eio%keep_beams, &
            keep_remnants = eio%keep_remnants)
       write (u, '(A)') "<event>"
       call hepeup_write_lhef (eio%unit)
       select case (eio%version)
         case ("2.0");  call eio%write_event_20 (event)
         case ("3.0");  call eio%write_event_30 (event)
       end select
       write (u, '(A)') "</event>"
    else
       call eio%write ()
       call msg_fatal ("LHEF file is not open for writing")
    end if
  end subroutine eio_lhef_output

  module subroutine eio_lhef_input_i_prc (eio, i_prc, iostat)
    class(eio_lhef_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    integer :: i, proc_num_id
    type(string_t) :: s
    logical :: ok
    iostat = 0
    call eio%tag_lhef%read_content (eio%cstream, s, ok)
    if (ok) then
       if (s == "") then
          iostat = -1
       else
          call err_close
       end if
       return
    else
       call eio%cstream%revert_record (s)
    end if
    call eio%tag_event%read (eio%cstream, ok)
    if (.not. ok) then
       call err_evt1
       return
    end if
    call hepeup_read_lhef (eio%unit)
    call hepeup_get_event_parameters (proc_id = proc_num_id)
    i_prc = 0
    FIND_I_PRC: do i = 1, size (eio%proc_num_id)
       if (eio%proc_num_id(i) == proc_num_id) then
          i_prc = i
          exit FIND_I_PRC
       end if
    end do FIND_I_PRC
    if (i_prc == 0)  call err_index
  contains
    subroutine err_close
      call msg_error ("LHEF: reading events: syntax error in closing tag")
      iostat = 1
    end subroutine
    subroutine err_evt1
      call msg_error ("LHEF: reading events: invalid event tag, &
           &aborting read")
      iostat = 2
    end subroutine err_evt1
    subroutine err_index
      call msg_error ("LHEF: reading events: undefined process ID " &
           // char (str (proc_num_id)) // ", aborting read")
      iostat = 3
    end subroutine err_index
  end subroutine eio_lhef_input_i_prc

  module subroutine eio_lhef_input_event (eio, event, iostat, event_handle)
    class(eio_lhef_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    type(string_t) :: s
    logical :: closing
    iostat = 0
    call event%reset_contents ()
    call event%select (1, 1, 1)
    call hepeup_to_event (event, eio%fallback_model, &
         recover_beams = eio%recover_beams, &
         use_alpha_s = eio%use_alphas_from_file, &
         use_scale = eio%use_scale_from_file)
    select case (eio%version)
    case ("1.0")
       call eio%tag_event%read_content (eio%cstream, s, closing = closing)
       if (s /= "" .or. .not. closing)  call err_evt2
    case ("2.0");  call eio%read_event_20 (event)
    case ("3.0");  call eio%read_event_30 (event)
    end select
    call event%increment_index ()
  contains
    subroutine err_evt2
      call msg_error ("LHEF: reading events: syntax error in event record, &
           &aborting read")
      iostat = 2
    end subroutine err_evt2

  end subroutine eio_lhef_input_event

  module subroutine eio_lhef_skip (eio, iostat)
    class(eio_lhef_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    if (eio%reading) then
       read (eio%unit, iostat = iostat)
    else
       call eio%write ()
       call msg_fatal ("Raw event file is not open for reading")
    end if
  end subroutine eio_lhef_skip

  module subroutine eio_lhef_write_header (eio, is_width)
    class(eio_lhef_t), intent(in) :: eio
    logical, intent(in), optional :: is_width
    logical :: is_w
    integer :: u
    u = given_output_unit (eio%unit);  if (u < 0)  return
    is_w = .false.
    if (present (is_width))  is_w = is_width
    call eio%tag_lhef%write (u);  write (u, *)
    call eio%tag_head%write (u);  write (u, *)
    select case (eio%version)
    case ("1.0")
       write (u, "(2x)", advance = "no")
       call eio%tag_gen_n%write (var_str ("WHIZARD"), u)
       write (u, *)
       write (u, "(2x)", advance = "no")
       call eio%tag_gen_v%write (var_str ("3.1.0"), u)
       write (u, *)
    end select
    if (is_w) then
       call eio%tag_whiz_info%write (u); write (u, *)
       write (u, "(A)") &
            "#  Special LHE event setup for decays, units in GeV"
       call eio%tag_whiz_info%close (u); write (u, *)
    end if
    call eio%tag_head%close (u);  write (u, *)
  end subroutine eio_lhef_write_header

  module subroutine eio_lhef_write_footer (eio)
    class(eio_lhef_t), intent(in) :: eio
    integer :: u
    u = given_output_unit (eio%unit);  if (u < 0)  return
    call eio%tag_lhef%close (u)
  end subroutine eio_lhef_write_footer

  module subroutine eio_lhef_read_header (eio)
    class(eio_lhef_t), intent(inout) :: eio
    logical :: success, closing
    type(string_t) :: content
    call eio%tag_lhef%read (eio%cstream, success)
    if (.not. success .or. .not. eio%tag_lhef%has_content)  call err_lhef
    if (eio%tag_lhef%get_attribute (1) /= eio%version)  call err_version
    call eio%tag_head%read (eio%cstream, success)
    if (.not. success)  call err_header
    if (eio%tag_head%has_content) then
       SKIP_HEADER_CONTENT: do
          call eio%tag_head%read_content (eio%cstream, content, closing)
          if (closing)  exit SKIP_HEADER_CONTENT
       end do SKIP_HEADER_CONTENT
    end if
  contains
    subroutine err_lhef
      call msg_fatal ("LHEF: LesHouchesEvents tag absent or corrupted")
    end subroutine err_lhef
    subroutine err_header
      call msg_fatal ("LHEF: header tag absent or corrupted")
    end subroutine err_header
    subroutine err_version
       call msg_error ("LHEF: version mismatch: expected " &
            // eio%version // ", found " &
            // char (eio%tag_lhef%get_attribute (1)))
    end subroutine err_version
  end subroutine eio_lhef_read_header

  module subroutine eio_lhef_read_init_10 (eio, data)
    class(eio_lhef_t), intent(in) :: eio
    type(event_sample_data_t), intent(out) :: data
    integer :: n_proc, i
    logical :: is_width
    call heprup_read_lhef (eio%unit)
    call heprup_get_run_parameters (n_processes = n_proc)
    call data%init (n_proc)
    if (IDBMUP(2) == 0) then
       data%n_beam = 1
       is_width = .true.
    else
       data%n_beam = 2
       is_width = .false.
    end if
    call heprup_get_run_parameters ( &
         unweighted = data%unweighted, &
         negative_weights = data%negative_weights, &
         beam_pdg = data%pdg_beam, &
         beam_energy = data%energy_beam)
    if (data%unweighted) then
       data%norm_mode = NORM_UNIT
    else
       data%norm_mode = NORM_SIGMA
    end if
    do i = 1, n_proc
       call heprup_get_process_parameters (i, &
            process_id = data%proc_num_id(i), &
            cross_section = data%cross_section(i), &
            error = data%error(i), &
            is_width = is_width)
    end do
  end subroutine eio_lhef_read_init_10

  module subroutine eio_lhef_write_init_20 (eio, data)
    class(eio_lhef_t), intent(in) :: eio
    type(event_sample_data_t), intent(in) :: data
    integer :: u
    u = eio%unit
    call eio%tag_generator%write (u)
    write (u, "(A)", advance="no")  "WHIZARD"
    call eio%tag_generator%close (u);  write (u, *)
    call eio%tag_xsecinfo%write (u);  write (u, *)
  end subroutine eio_lhef_write_init_20

  module subroutine eio_lhef_read_init_20 (eio, data)
    class(eio_lhef_t), intent(inout) :: eio
    type(event_sample_data_t), intent(out) :: data
    real(default), parameter :: pb_per_fb = 1.e-3_default
    type(string_t) :: content
    logical :: found, closing
    call eio_lhef_read_init_10 (eio, data)
    SCAN_INIT_TAGS: do
       call eio%tag_generator%read (eio%cstream, found)
       if (found) then
          if (.not. eio%tag_generator%has_content)  call err_generator
          call eio%tag_generator%read_content (eio%cstream, content, closing)
          call msg_message ("LHEF: Event file has been generated by " &
               // char (content) // " " &
               // char (eio%tag_generator%get_attribute (1)))
          cycle SCAN_INIT_TAGS
       end if
       call eio%tag_xsecinfo%read (eio%cstream, found)
       if (found) then
          if (eio%tag_xsecinfo%has_content)  call err_xsecinfo
          cycle SCAN_INIT_TAGS
       end if
       call eio%tag_init%read_content (eio%cstream, content, closing)
       if (closing) then
          if (content /= "")  call err_init
          exit SCAN_INIT_TAGS
       end if
    end do SCAN_INIT_TAGS
    data%n_evt = &
         read_ival (eio%tag_xsecinfo%get_attribute (1))
    if (data%n_beam == 1) then
       data%total_cross_section = &
            read_rval (eio%tag_xsecinfo%get_attribute (2))
    else
       data%total_cross_section = &
            read_rval (eio%tag_xsecinfo%get_attribute (2)) / pb_per_fb
    end if
  contains
    subroutine err_generator
      call msg_fatal ("LHEF: invalid generator tag")
    end subroutine err_generator
    subroutine err_xsecinfo
      call msg_fatal ("LHEF: invalid xsecinfo tag")
    end subroutine err_xsecinfo
    subroutine err_init
      call msg_fatal ("LHEF: syntax error after init tag")
    end subroutine err_init
  end subroutine eio_lhef_read_init_20

  module subroutine eio_lhef_write_event_20 (eio, event)
    class(eio_lhef_t), intent(in) :: eio
    class(generic_event_t), intent(in) :: event
    type(string_t) :: s
    integer :: i, u
    u = eio%unit
    if (eio%write_sqme_ref) then
       s = str (event%get_sqme_ref ())
       call eio%tag_sqme_ref%write (s, u);  write (u, *)
    end if
    if (eio%write_sqme_prc) then
       s = str (event%get_sqme_prc ())
       call eio%tag_sqme_prc%write (s, u);  write (u, *)
    end if
    if (eio%n_alt > 0) then
       if (eio%write_sqme_alt) then
          s = str (event%get_sqme_alt(1))
          do i = 2, eio%n_alt
             s = s // " " // str (event%get_sqme_alt(i));  write (u, *)
          end do
          call eio%tag_sqme_alt(1)%write (s, u)
       end if
       s = str (event%get_weight_alt(1))
       do i = 2, eio%n_alt
          s = s // " " // str (event%get_weight_alt(i));  write (u, *)
       end do
       call eio%tag_wgts_alt(1)%write (s, u)
    end if
  end subroutine eio_lhef_write_event_20

  module subroutine eio_lhef_read_event_20 (eio, event)
    class(eio_lhef_t), intent(inout) :: eio
    class(generic_event_t), intent(inout) :: event
    type(string_t) :: content
    logical :: found, closing
    SCAN_EVENT_TAGS: do
       call eio%tag_weight%read (eio%cstream, found)
       if (found) then
          if (.not. eio%tag_weight%has_content)  call err_weight
          call eio%tag_weight%read_content (eio%cstream, content, closing)
          if (.not. closing)  call err_weight
          if (eio%tag_weight%get_attribute (1) == "sqme_prc") then
             call event%set_sqme_ref (read_rval (content))
          end if
          cycle SCAN_EVENT_TAGS
       end if
       call eio%tag_event%read_content (eio%cstream, content, closing)
       if (closing) then
          if (content /= "")  call err_event
          exit SCAN_EVENT_TAGS
       end if
    end do SCAN_EVENT_TAGS
  contains
    subroutine err_weight
      call msg_fatal ("LHEF: invalid weight tag in event record")
    end subroutine err_weight
    subroutine err_event
      call msg_fatal ("LHEF: syntax error after event tag")
    end subroutine err_event
  end subroutine eio_lhef_read_event_20

  module subroutine eio_lhef_write_init_30 (eio, data)
    class(eio_lhef_t), intent(in) :: eio
    type(event_sample_data_t), intent(in) :: data
    integer :: u, i
    u = given_output_unit (eio%unit)
    call eio%tag_generator%write (u)
    write (u, "(A)", advance="no")  "WHIZARD"
    call eio%tag_generator%close (u);  write (u, *)
    call eio%tag_xsecinfo%write (u);  write (u, *)
    if (eio%write_sqme_ref) then
       call eio%tag_sqme_ref%write (u);  write (u, *)
    end if
    if (eio%write_sqme_prc) then
       call eio%tag_sqme_prc%write (u);  write (u, *)
    end if
    if (eio%write_sqme_alt) then
       do i = 1, eio%n_alt
          call eio%tag_sqme_alt(i)%write (u);  write (u, *)
       end do
    end if
    do i = 1, eio%n_alt
       call eio%tag_wgts_alt(i)%write (u);  write (u, *)
    end do
  end subroutine eio_lhef_write_init_30

  module subroutine eio_lhef_read_init_30 (eio, data)
    class(eio_lhef_t), intent(inout) :: eio
    type(event_sample_data_t), intent(out) :: data
    real(default), parameter :: pb_per_fb = 1.e-3_default
    type(string_t) :: content
    logical :: found, closing
    integer :: n_weightinfo
    call eio_lhef_read_init_10 (eio, data)
    n_weightinfo = 0
    eio%i_weight_sqme = 0
    SCAN_INIT_TAGS: do
       call eio%tag_generator%read (eio%cstream, found)
       if (found) then
          if (.not. eio%tag_generator%has_content)  call err_generator
          call eio%tag_generator%read_content (eio%cstream, content, closing)
          call msg_message ("LHEF: Event file has been generated by " &
               // char (content) // " " &
               // char (eio%tag_generator%get_attribute (1)))
          cycle SCAN_INIT_TAGS
       end if
       call eio%tag_xsecinfo%read (eio%cstream, found)
       if (found) then
          if (eio%tag_xsecinfo%has_content)  call err_xsecinfo
          cycle SCAN_INIT_TAGS
       end if
       call eio%tag_weightinfo%read (eio%cstream, found)
       if (found) then
          if (eio%tag_weightinfo%has_content)  call err_xsecinfo
          n_weightinfo = n_weightinfo + 1
          if (eio%tag_weightinfo%get_attribute (1) == "sqme_prc") then
             eio%i_weight_sqme = n_weightinfo
          end if
          cycle SCAN_INIT_TAGS
       end if
       call eio%tag_init%read_content (eio%cstream, content, closing)
       if (closing) then
          if (content /= "")  call err_init
          exit SCAN_INIT_TAGS
       end if
    end do SCAN_INIT_TAGS
    data%n_evt = &
         read_ival (eio%tag_xsecinfo%get_attribute (1))
    if (data%n_beam == 1) then
       data%total_cross_section = &
            read_rval (eio%tag_xsecinfo%get_attribute (2))
    else
       data%total_cross_section = &
            read_rval (eio%tag_xsecinfo%get_attribute (2)) / pb_per_fb
    end if
  contains
    subroutine err_generator
      call msg_fatal ("LHEF: invalid generator tag")
    end subroutine err_generator
    subroutine err_xsecinfo
      call msg_fatal ("LHEF: invalid xsecinfo tag")
    end subroutine err_xsecinfo
    subroutine err_init
      call msg_fatal ("LHEF: syntax error after init tag")
    end subroutine err_init
  end subroutine eio_lhef_read_init_30

  module subroutine eio_lhef_write_event_30 (eio, event)
    class(eio_lhef_t), intent(in) :: eio
    class(generic_event_t), intent(in) :: event
    type(string_t) :: s
    integer :: u, i
    u = eio%unit
    s = ""
    if (eio%write_sqme_ref) then
       s = s // str (event%get_sqme_ref ()) // " "
    end if
    if (eio%write_sqme_prc) then
       s = s // str (event%get_sqme_prc ()) // " "
    end if
    if (eio%n_alt > 0) then
       if (eio%write_sqme_alt) then
          s = s // str (event%get_sqme_alt(1)) // " "
          do i = 2, eio%n_alt
             s = s // str (event%get_sqme_alt(i)) // " "
          end do
       end if
       s = s // str (event%get_weight_alt(1)) // " "
       do i = 2, eio%n_alt
          s = s // str (event%get_weight_alt(i)) // " "
       end do
    end if
    if (len_trim (s) > 0) then
       call eio%tag_weights%write (trim (s), u);  write (u, *)
    end if
  end subroutine eio_lhef_write_event_30

  module subroutine eio_lhef_read_event_30 (eio, event)
    class(eio_lhef_t), intent(inout) :: eio
    class(generic_event_t), intent(inout) :: event
    type(string_t) :: content, string
    logical :: found, closing
    integer :: i
    SCAN_EVENT_TAGS: do
       call eio%tag_weights%read (eio%cstream, found)
       if (found) then
          if (.not. eio%tag_weights%has_content)  call err_weights
          call eio%tag_weights%read_content (eio%cstream, content, closing)
          if (.not. closing)  call err_weights
          if (eio%i_weight_sqme > 0) then
             SCAN_WEIGHTS: do i = 1, eio%i_weight_sqme
                call split (content, string, " ")
                content = adjustl (content)
                if (i == eio%i_weight_sqme) then
                   call event%set_sqme_ref (read_rval (string))
                   exit SCAN_WEIGHTS
                end if
             end do SCAN_WEIGHTS
          end if
          cycle SCAN_EVENT_TAGS
       end if
       call eio%tag_event%read_content (eio%cstream, content, closing)
       if (closing) then
          if (content /= "")  call err_event
          exit SCAN_EVENT_TAGS
       end if
    end do SCAN_EVENT_TAGS
  contains
    subroutine err_weights
      call msg_fatal ("LHEF: invalid weights tag in event record")
    end subroutine err_weights
    subroutine err_event
      call msg_fatal ("LHEF: syntax error after event tag")
    end subroutine err_event
  end subroutine eio_lhef_read_event_30


end submodule eio_lhef_s

