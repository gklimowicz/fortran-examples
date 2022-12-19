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

module event_streams_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use model_data
  use eio_data
  use process, only: process_t
  use instances, only: process_instance_t
  use models
  use rt_data
  use events

  use event_streams

  implicit none
  private

  public :: event_streams_1
  public :: event_streams_2
  public :: event_streams_3
  public :: event_streams_4

contains

  subroutine event_streams_1 (u)
    integer, intent(in) :: u
    type(event_stream_array_t) :: es_array
    type(rt_data_t) :: global
    type(event_t) :: event
    type(string_t) :: sample
    type(string_t), dimension(0) :: empty_string_array

    write (u, "(A)")  "* Test output: event_streams_1"
    write (u, "(A)")  "*   Purpose: handle empty event stream array"
    write (u, "(A)")

    sample = "event_streams_1"

    call es_array%init (sample, empty_string_array, global)
    call es_array%output (event, 42, 1)
    call es_array%write (u)
    call es_array%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: event_streams_1"

  end subroutine event_streams_1

  subroutine event_streams_2 (u)
    use processes_ut, only: prepare_test_process
    integer, intent(in) :: u
    type(event_stream_array_t) :: es_array
    type(rt_data_t) :: global
    type(model_data_t), target :: model
    type(event_t), allocatable, target :: event
    type(process_t), allocatable, target :: process
    type(process_instance_t), allocatable, target :: process_instance
    type(string_t) :: sample
    type(string_t), dimension(0) :: empty_string_array
    integer :: i_prc, iostat

    write (u, "(A)")  "* Test output: event_streams_2"
    write (u, "(A)")  "*   Purpose: handle empty event stream array"
    write (u, "(A)")

    call syntax_model_file_init ()
    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call model%init_test ()

    write (u, "(A)")  "* Generate test process event"
    write (u, "(A)")

    allocate (process)
    allocate (process_instance)
    call prepare_test_process (process, process_instance, model, &
         run_id = var_str ("run_test"))
    call process_instance%setup_event_data ()

    allocate (event)
    call event%basic_init ()
    call event%connect (process_instance, process%get_model_ptr ())
    call event%generate (1, [0.4_default, 0.4_default])
    call event%set_index (42)
    call event%evaluate_expressions ()
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Allocate raw eio stream and write event to file"
    write (u, "(A)")

    sample = "event_streams_2"

    call es_array%init (sample, [var_str ("raw")], global)
    call es_array%output (event, 1, 1)
    call es_array%write (u)
    call es_array%final ()

    write (u, "(A)")
    write (u, "(A)") "* Reallocate raw eio stream for reading"
    write (u, "(A)")

    sample = "foo"
    call es_array%init (sample, empty_string_array, global, &
         input = var_str ("raw"), input_sample = var_str ("event_streams_2"))
    call es_array%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Reread event"
    write (u, "(A)")

    call es_array%input_i_prc (i_prc, iostat)

    write (u, "(1x,A,I0)")  "i_prc = ", i_prc
    write (u, "(A)")
    call es_array%input_event (event, iostat)
    call es_array%final ()

    call event%write (u)

    call global%final ()

    call model%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: event_streams_2"

  end subroutine event_streams_2

  subroutine event_streams_3 (u)
    use processes_ut, only: prepare_test_process
    integer, intent(in) :: u
    type(event_stream_array_t) :: es_array
    type(rt_data_t) :: global
    type(model_data_t), target :: model
    type(event_t), allocatable, target :: event
    type(process_t), allocatable, target :: process
    type(process_instance_t), allocatable, target :: process_instance
    type(string_t) :: sample
    type(string_t), dimension(0) :: empty_string_array
    integer :: i_prc, iostat

    write (u, "(A)")  "* Test output: event_streams_3"
    write (u, "(A)")  "*   Purpose: handle in/out switching"
    write (u, "(A)")

    call syntax_model_file_init ()
    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call model%init_test ()

    write (u, "(A)")  "* Generate test process event"
    write (u, "(A)")

    allocate (process)
    allocate (process_instance)
    call prepare_test_process (process, process_instance, model, &
         run_id = var_str ("run_test"))
    call process_instance%setup_event_data ()

    allocate (event)
    call event%basic_init ()
    call event%connect (process_instance, process%get_model_ptr ())
    call event%generate (1, [0.4_default, 0.4_default])
    call event%increment_index ()
    call event%evaluate_expressions ()

    write (u, "(A)") "* Allocate raw eio stream and write event to file"
    write (u, "(A)")

    sample = "event_streams_3"

    call es_array%init (sample, [var_str ("raw")], global)
    call es_array%output (event, 1, 1)
    call es_array%write (u)
    call es_array%final ()

    write (u, "(A)")
    write (u, "(A)") "* Reallocate raw eio stream for reading"
    write (u, "(A)")

    call es_array%init (sample, empty_string_array, global, &
         input = var_str ("raw"))
    call es_array%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Reread event"
    write (u, "(A)")

    call es_array%input_i_prc (i_prc, iostat)
    call es_array%input_event (event, iostat)

    write (u, "(A)") "* Attempt to read another event (fail), then generate"
    write (u, "(A)")

    call es_array%input_i_prc (i_prc, iostat)
    if (iostat < 0) then
       call es_array%switch_inout ()
       call event%generate (1, [0.3_default, 0.3_default])
       call event%increment_index ()
       call event%evaluate_expressions ()
       call es_array%output (event, 1, 2)
    end if
    call es_array%write (u)
    call es_array%final ()

    write (u, "(A)")
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Reallocate raw eio stream for reading"
    write (u, "(A)")

    call es_array%init (sample, empty_string_array, global, &
         input = var_str ("raw"))
    call es_array%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Reread two events and display 2nd event"
    write (u, "(A)")

    call es_array%input_i_prc (i_prc, iostat)
    call es_array%input_event (event, iostat)
    call es_array%input_i_prc (i_prc, iostat)

    call es_array%input_event (event, iostat)
    call es_array%final ()

    call event%write (u)

    call global%final ()

    call model%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: event_streams_3"

  end subroutine event_streams_3

  subroutine event_streams_4 (u)
    integer, intent(in) :: u
    type(event_stream_array_t) :: es_array
    type(rt_data_t) :: global
    type(process_t), allocatable, target :: process
    type(string_t) :: sample
    type(string_t), dimension(0) :: empty_string_array
    type(event_sample_data_t) :: data

    write (u, "(A)")  "* Test output: event_streams_4"
    write (u, "(A)")  "*   Purpose: handle in/out switching"
    write (u, "(A)")

    write (u, "(A)")  "* Generate test process event"
    write (u, "(A)")

    call syntax_model_file_init ()
    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call global%set_log (var_str ("?check_event_file"), &
         .true., is_known = .true.)

    allocate (process)

    write (u, "(A)") "* Allocate raw eio stream for writing"
    write (u, "(A)")

    sample = "event_streams_4"
    data%md5sum_cfg = "1234567890abcdef1234567890abcdef"

    call es_array%init (sample, [var_str ("raw")], global, data)
    call es_array%write (u)
    call es_array%final ()

    write (u, "(A)")
    write (u, "(A)") "* Reallocate raw eio stream for reading"
    write (u, "(A)")

    call es_array%init (sample, empty_string_array, global, &
         data, input = var_str ("raw"))
    call es_array%write (u)
    call es_array%final ()

    write (u, "(A)")
    write (u, "(A)") "* Reallocate modified raw eio stream for reading (fail)"
    write (u, "(A)")

    data%md5sum_cfg = "1234567890______1234567890______"
    call es_array%init (sample, empty_string_array, global, &
         data, input = var_str ("raw"))
    call es_array%write (u)
    call es_array%final ()

    write (u, "(A)")
    write (u, "(A)") "* Repeat ignoring checksum"
    write (u, "(A)")

    call global%set_log (var_str ("?check_event_file"), &
         .false., is_known = .true.)
    call es_array%init (sample, empty_string_array, global, &
         data, input = var_str ("raw"))
    call es_array%write (u)
    call es_array%final ()

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: event_streams_4"

  end subroutine event_streams_4


end module event_streams_uti

