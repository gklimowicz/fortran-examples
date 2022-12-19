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

submodule (eio_raw) eio_raw_s

  use io_units
  use diagnostics
  use model_data
  use particles

  implicit none

contains

  module subroutine eio_raw_write (object, unit)
    class(eio_raw_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Raw event stream:"
    write (u, "(3x,A,L1)")  "Check MD5 sum     = ", object%check
    if (object%n_alt > 0) then
       write (u, "(3x,A,I0)")  "Alternate weights = ", object%n_alt
    end if
    write (u, "(3x,A,L1)")    "Alpha_s from file = ", &
         object%use_alphas_from_file
    write (u, "(3x,A,L1)")    "Scale from file   = ", &
         object%use_scale_from_file
    write (u, "(3x,A,L1)")    "Events for fNLO   = ", &
         object%fixed_order_nlo
    if (object%reading) then
       write (u, "(3x,A,A)")  "Reading from file = ", char (object%filename)
    else if (object%writing) then
       write (u, "(3x,A,A)")  "Writing to file   = ", char (object%filename)
    else
       write (u, "(3x,A)")  "[closed]"
    end if
  end subroutine eio_raw_write

  module subroutine eio_raw_final (object)
    class(eio_raw_t), intent(inout) :: object
    if (object%reading .or. object%writing) then
       write (msg_buffer, "(A,A,A)")  "Events: closing raw file '", &
            char (object%filename), "'"
       call msg_message ()
       close (object%unit)
       object%reading = .false.
       object%writing = .false.
    end if
  end subroutine eio_raw_final

  module subroutine eio_raw_set_parameters (eio, check, use_alphas_from_file, &
       use_scale_from_file, fixed_order_nlo, version_string, extension)
    class(eio_raw_t), intent(inout) :: eio
    logical, intent(in), optional :: check, use_alphas_from_file, &
         use_scale_from_file, fixed_order_nlo
    type(string_t), intent(in), optional :: version_string
    type(string_t), intent(in), optional :: extension
    if (present (check))  eio%check = check
    if (present (use_alphas_from_file))  eio%use_alphas_from_file = &
         use_alphas_from_file
    if (present (use_scale_from_file))  eio%use_scale_from_file = &
         use_scale_from_file
    if (present (fixed_order_nlo))  eio%fixed_order_nlo = &
         fixed_order_nlo
    if (present (version_string)) then
       select case (char (version_string))
       case ("", "2.2.4")
          eio%file_version = CURRENT_FILE_VERSION
       case ("2.2")
          eio%file_version = 1
       case default
          call msg_fatal ("Raw event I/O: unsupported version '" &
               // char (version_string) // "'")
          eio%file_version = 0
       end select
    end if
    if (present (extension)) then
       eio%extension = extension
    else
       eio%extension = "evx"
    end if
  end subroutine eio_raw_set_parameters

  module subroutine eio_raw_init_out (eio, sample, data, success, extension)
    class(eio_raw_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    type(string_t), intent(in), optional :: extension
    character(32) :: md5sum_prc, md5sum_cfg
    character(32), dimension(:), allocatable :: md5sum_alt
    integer :: i
    if (present (extension)) then
       eio%extension  = extension
    else
       eio%extension = "evx"
    end if
    eio%filename = sample // "." // eio%extension
    eio%unit = free_unit ()
    write (msg_buffer, "(A,A,A)")  "Events: writing to raw file '", &
         char (eio%filename), "'"
    call msg_message ()
    eio%writing = .true.
    if (present (data)) then
       md5sum_prc = data%md5sum_prc
       md5sum_cfg = data%md5sum_cfg
       eio%norm_mode = data%norm_mode
       eio%sigma = data%total_cross_section
       eio%n = data%n_evt
       eio%n_alt = data%n_alt
       if (eio%n_alt > 0) then
          !!! !!! !!! Workaround for gfortran 5.0 ICE
          allocate (md5sum_alt (data%n_alt))
          md5sum_alt = data%md5sum_alt
          !!! allocate (md5sum_alt (data%n_alt), source = data%md5sum_alt)
       end if
    else
       md5sum_prc = ""
       md5sum_cfg = ""
    end if
    open (eio%unit, file = char (eio%filename), form = "unformatted", &
         action = "write", status = "replace")
    select case (eio%file_version)
    case (2:);  write (eio%unit)  eio%file_version
    end select
    write (eio%unit)  md5sum_prc
    write (eio%unit)  md5sum_cfg
    write (eio%unit)  eio%norm_mode
    write (eio%unit)  eio%n_alt
    if (allocated (md5sum_alt)) then
       do i = 1, eio%n_alt
          write (eio%unit)  md5sum_alt(i)
       end do
    end if
    if (present (success))  success = .true.
  end subroutine eio_raw_init_out

  module subroutine eio_raw_init_in (eio, sample, data, success, extension)
    class(eio_raw_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    type(string_t), intent(in), optional :: extension
    character(32) :: md5sum_prc, md5sum_cfg
    character(32), dimension(:), allocatable :: md5sum_alt
    integer :: i, file_version
    if (present (success))  success = .true.
    if (present (extension)) then
       eio%extension = extension
    else
       eio%extension = "evx"
    end if
    eio%filename = sample // "." // eio%extension
    eio%unit = free_unit ()
    if (present (data)) then
       eio%sigma = data%total_cross_section
       eio%n = data%n_evt
    end if
    write (msg_buffer, "(A,A,A)")  "Events: reading from raw file '", &
         char (eio%filename), "'"
    call msg_message ()
    eio%reading = .true.
    open (eio%unit, file = char (eio%filename), form = "unformatted", &
         action = "read", status = "old")
    select case (eio%file_version)
    case (2:);  read (eio%unit)  file_version
    case default;  file_version = 1
    end select
    if (file_version /= eio%file_version) then
       call msg_error ("Reading event file: raw-file version mismatch.")
       if (present (success))  success = .false.
       return
    else if (file_version /= CURRENT_FILE_VERSION) then
       call msg_warning ("Reading event file: compatibility mode.")
    end if
    read (eio%unit)  md5sum_prc
    read (eio%unit)  md5sum_cfg
    read (eio%unit)  eio%norm_mode
    read (eio%unit)  eio%n_alt
    if (present (data)) then
       if (eio%n_alt /= data%n_alt) then
          if (present (success))  success = .false.
          return
       end if
    end if
    allocate (md5sum_alt (eio%n_alt))
    do i = 1, eio%n_alt
       read (eio%unit)  md5sum_alt(i)
    end do
    if (present (success)) then
       if (present (data)) then
          if (eio%check) then
             if (data%md5sum_prc /= "") then
                success = success .and. md5sum_prc == data%md5sum_prc
             end if
             if (data%md5sum_cfg /= "") then
                success = success .and. md5sum_cfg == data%md5sum_cfg
             end if
             do i = 1, eio%n_alt
                if (data%md5sum_alt(i) /= "") then
                   success = success .and. md5sum_alt(i) == data%md5sum_alt(i)
                end if
             end do
          else
             call msg_warning ("Reading event file: MD5 sum check disabled")
          end if
       end if
    end if
  end subroutine eio_raw_init_in

  module subroutine eio_raw_switch_inout (eio, success)
    class(eio_raw_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    write (msg_buffer, "(A,A,A)")  "Events: appending to raw file '", &
         char (eio%filename), "'"
    call msg_message ()
    close (eio%unit, status = "keep")
    eio%reading = .false.
    open (eio%unit, file = char (eio%filename), form = "unformatted", &
         action = "write", position = "append", status = "old")
    eio%writing = .true.
    if (present (success))  success = .true.
  end subroutine eio_raw_switch_inout

  module subroutine eio_raw_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_raw_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    integer, intent(in) :: i_prc
    type(particle_set_t), pointer :: pset
    integer :: i
    if (eio%writing) then
       if (event%has_valid_particle_set ()) then
          select type (event)
          type is (event_t)
             write (eio%unit)  i_prc
             write (eio%unit)  event%get_index ()
             write (eio%unit)  event%get_i_mci ()
             write (eio%unit)  event%get_i_term ()
             write (eio%unit)  event%get_channel ()
             write (eio%unit)  event%expr%weight_prc
             write (eio%unit)  event%expr%excess_prc
             write (eio%unit)  event%get_n_dropped ()
             write (eio%unit)  event%expr%sqme_prc
             do i = 1, eio%n_alt
                write (eio%unit)  event%expr%weight_alt(i)
                write (eio%unit)  event%expr%sqme_alt(i)
             end do
             allocate (pset)
             call event%get_hard_particle_set (pset)
             call pset%write_raw (eio%unit)
             call pset%final ()
             deallocate (pset)
             select case (eio%file_version)
             case (2:)
                if (event%has_transform ()) then
                   write (eio%unit)  .true.
                   pset => event%get_particle_set_ptr ()
                   call pset%write_raw (eio%unit)
                else
                   write (eio%unit)  .false.
                end if
             end select
          class default
             call msg_bug ("Event: write raw: defined only for full event_t")
          end select
       else
          call msg_bug ("Event: write raw: particle set is undefined")
       end if
    else
       call eio%write ()
       call msg_fatal ("Raw event file is not open for writing")
    end if
  end subroutine eio_raw_output

  module subroutine eio_raw_input_i_prc (eio, i_prc, iostat)
    class(eio_raw_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    if (eio%reading) then
       read (eio%unit, iostat = iostat)  i_prc
    else
       call eio%write ()
       call msg_fatal ("Raw event file is not open for reading")
    end if
  end subroutine eio_raw_input_i_prc

  module subroutine eio_raw_input_event (eio, event, iostat, event_handle)
    class(eio_raw_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    integer :: event_index, i_mci, i_term, channel, i
    real(default) :: weight, excess, sqme
    integer :: n_dropped
    real(default), dimension(:), allocatable :: weight_alt, sqme_alt
    logical :: has_transform
    type(particle_set_t), pointer :: pset
    class(model_data_t), pointer :: model
    if (eio%reading) then
       select type (event)
       type is (event_t)
          read (eio%unit, iostat = iostat)  event_index
          if (iostat /= 0)  return
          read (eio%unit, iostat = iostat)  i_mci
          if (iostat /= 0)  return
          read (eio%unit, iostat = iostat)  i_term
          if (iostat /= 0)  return
          read (eio%unit, iostat = iostat)  channel
          if (iostat /= 0)  return
          read (eio%unit, iostat = iostat)  weight
          if (iostat /= 0)  return
          read (eio%unit, iostat = iostat)  excess
          if (iostat /= 0)  return
          read (eio%unit, iostat = iostat)  n_dropped
          if (iostat /= 0)  return
          read (eio%unit, iostat = iostat)  sqme
          if (iostat /= 0)  return
          call event%reset_contents ()
          call event%set_index (event_index)
          call event%select (i_mci, i_term, channel)
          if (eio%norm_mode /= NORM_UNDEFINED) then
             call event_normalization_update (weight, &
                  eio%sigma, eio%n, event%get_norm_mode (), eio%norm_mode)
             call event_normalization_update (excess, &
                  eio%sigma, eio%n, event%get_norm_mode (), eio%norm_mode)
          end if
          call event%set (sqme_ref = sqme, weight_ref = weight, &
               excess_prc = excess, &
               n_dropped = n_dropped)
          if (eio%n_alt /= 0) then
             allocate (sqme_alt (eio%n_alt), weight_alt (eio%n_alt))
             do i = 1, eio%n_alt
                read (eio%unit, iostat = iostat)  weight_alt(i)
                if (iostat /= 0)  return
                read (eio%unit, iostat = iostat)  sqme_alt(i)
                if (iostat /= 0)  return
             end do
             call event%set (sqme_alt = sqme_alt, weight_alt = weight_alt)
          end if
          model => null ()
          if (associated (event%process)) then
             model => event%process%get_model_ptr ()
          end if
          allocate (pset)
          call pset%read_raw (eio%unit, iostat)
          if (iostat /= 0)  return
          if (associated (model))  call pset%set_model (model)
          call event%set_hard_particle_set (pset)
          if (eio%use_alphas_from_file .or. eio%use_scale_from_file) then
             call event%recalculate (update_sqme = .true.)
             if (eio%fixed_order_nlo) then
                if (event%weight_prc /= event%weight_ref .and. &
                     event%weight_prc == 0) then
                   event%weight_prc = event%weight_ref
                end if
             end if
          end if
          call pset%final ()
          deallocate (pset)
          select case (eio%file_version)
          case (2:)
             read (eio%unit, iostat = iostat)  has_transform
             if (iostat /= 0)  return
             if (has_transform) then
                allocate (pset)
                call pset%read_raw (eio%unit, iostat)
                if (iostat /= 0)  return
                if (associated (model)) &
                     call pset%set_model (model)
                call event%link_particle_set (pset)
             end if
          end select
       class default
          call msg_bug ("Event: read raw: defined only for full event_t")
       end select
    else
       call eio%write ()
       call msg_fatal ("Raw event file is not open for reading")
    end if
  end subroutine eio_raw_input_event

  module subroutine eio_raw_skip (eio, iostat)
    class(eio_raw_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    if (eio%reading) then
       read (eio%unit, iostat = iostat)
    else
       call eio%write ()
       call msg_fatal ("Raw event file is not open for reading")
    end if
  end subroutine eio_raw_skip


end submodule eio_raw_s

