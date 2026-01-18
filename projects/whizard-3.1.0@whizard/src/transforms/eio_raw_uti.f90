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

module eio_raw_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use model_data
  use variables
  use events
  use eio_data
  use eio_base

  use eio_raw

  use process, only: process_t
  use instances, only: process_instance_t

  implicit none
  private

  public :: eio_raw_1
  public :: eio_raw_2

contains

  subroutine eio_raw_1 (u)
    use processes_ut, only: prepare_test_process, cleanup_test_process
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(event_t), allocatable, target :: event
    type(process_t), allocatable, target :: process
    type(process_instance_t), allocatable, target :: process_instance
    class(eio_t), allocatable :: eio
    integer :: i_prc, iostat
    type(string_t) :: sample

    write (u, "(A)")  "* Test output: eio_raw_1"
    write (u, "(A)")  "*   Purpose: generate and read/write an event"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call model%init_test ()

    allocate (process)
    allocate (process_instance)
    call prepare_test_process (process, process_instance, model, &
         run_id = var_str ("run_test"))
    call process_instance%setup_event_data ()

    allocate (event)
    call event%basic_init ()
    call event%connect (process_instance, process%get_model_ptr ())

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_raw_1"

    allocate (eio_raw_t :: eio)

    call eio%init_out (sample)
    call event%generate (1, [0._default, 0._default])
    call event%increment_index ()
    call event%evaluate_expressions ()
    call event%write (u)
    write (u, "(A)")

    call eio%output (event, i_prc = 42)
    call eio%write (u)
    call eio%final ()

    call event%final ()
    deallocate (event)
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Re-read the event"
    write (u, "(A)")

    call eio%init_in (sample)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data ()
    allocate (event)
    call event%basic_init ()
    call event%connect (process_instance, process%get_model_ptr ())

    call eio%input_i_prc (i_prc, iostat)
    if (iostat /= 0)  write (u, "(A,I0)")  "I/O error (i_prc):", iostat
    call eio%input_event (event, iostat)
    if (iostat /= 0)  write (u, "(A,I0)")  "I/O error (event):", iostat
    call eio%write (u)

    write (u, "(A)")
    write (u, "(1x,A,I0)")  "i_prc = ", i_prc
    write (u, "(A)")
    call event%write (u)
    write (u, "(A)")
    write (u, "(A)")  "* Generate and append another event"
    write (u, "(A)")

    call eio%switch_inout ()
    call event%generate (1, [0._default, 0._default])
    call event%increment_index ()
    call event%evaluate_expressions ()
    call event%write (u)
    write (u, "(A)")

    call eio%output (event, i_prc = 5)
    call eio%write (u)
    call eio%final ()

    call event%final ()
    deallocate (event)
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Re-read both events"
    write (u, "(A)")

    call eio%init_in (sample)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data ()
    allocate (event)
    call event%basic_init ()
    call event%connect (process_instance, process%get_model_ptr ())

    call eio%input_i_prc (i_prc, iostat)
    if (iostat /= 0)  write (u, "(A,I0)")  "I/O error (i_prc/1):", iostat
    call eio%input_event (event, iostat)
    if (iostat /= 0)  write (u, "(A,I0)")  "I/O error (event/1):", iostat
    call eio%input_i_prc (i_prc, iostat)
    if (iostat /= 0)  write (u, "(A,I0)")  "I/O error (i_prc/2):", iostat
    call eio%input_event (event, iostat)
    if (iostat /= 0)  write (u, "(A,I0)")  "I/O error (event/2):", iostat
    call eio%write (u)

    write (u, "(A)")
    write (u, "(1x,A,I0)")  "i_prc = ", i_prc
    write (u, "(A)")
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio%final ()
    deallocate (eio)

    call event%final ()
    deallocate (event)

    call cleanup_test_process (process, process_instance)
    deallocate (process_instance)
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_raw_1"

  end subroutine eio_raw_1

  subroutine eio_raw_2 (u)
    use processes_ut, only: prepare_test_process, cleanup_test_process
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(var_list_t) :: var_list
    type(event_t), allocatable, target :: event
    type(process_t), allocatable, target :: process
    type(process_instance_t), allocatable, target :: process_instance
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    integer :: i_prc, iostat
    type(string_t) :: sample

    write (u, "(A)")  "* Test output: eio_raw_2"
    write (u, "(A)")  "*   Purpose: generate and read/write an event"
    write (u, "(A)")  "*            with multiple weights"
    write (u, "(A)")

    call model%init_test ()

    write (u, "(A)")  "* Initialize test process"

    allocate (process)
    allocate (process_instance)
    call prepare_test_process (process, process_instance, model, &
         run_id = var_str ("run_test"))
    call process_instance%setup_event_data ()

    call data%init (n_proc = 1, n_alt = 2)

    call var_list%append_log (var_str ("?unweighted"), .false., &
         intrinsic = .true.)
    call var_list%append_string (var_str ("$sample_normalization"), &
         var_str ("auto"), intrinsic = .true.)
    call var_list%append_real (var_str ("safety_factor"), &
         1._default, intrinsic = .true.)

    allocate (event)
    call event%basic_init (var_list, n_alt = 2)
    call event%connect (process_instance, process%get_model_ptr ())

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_raw_2"

    allocate (eio_raw_t :: eio)

    call eio%init_out (sample, data)
    call event%generate (1, [0._default, 0._default])
    call event%increment_index ()
    call event%evaluate_expressions ()
    call event%set (sqme_alt = [2._default, 3._default])
    call event%set (weight_alt = &
         [2 * event%get_weight_ref (), 3 * event%get_weight_ref ()])
    call event%store_alt_values ()
    call event%check ()

    call event%write (u)
    write (u, "(A)")

    call eio%output (event, i_prc = 42)
    call eio%write (u)
    call eio%final ()

    call event%final ()
    deallocate (event)
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Re-read the event"
    write (u, "(A)")

    call eio%init_in (sample, data)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data ()
    allocate (event)
    call event%basic_init (var_list, n_alt = 2)
    call event%connect (process_instance, process%get_model_ptr ())

    call eio%input_i_prc (i_prc, iostat)
    if (iostat /= 0)  write (u, "(A,I0)")  "I/O error (i_prc):", iostat
    call eio%input_event (event, iostat)
    if (iostat /= 0)  write (u, "(A,I0)")  "I/O error (event):", iostat
    call eio%write (u)

    write (u, "(A)")
    write (u, "(1x,A,I0)")  "i_prc = ", i_prc
    write (u, "(A)")
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio%final ()
    deallocate (eio)

    call event%final ()
    deallocate (event)

    call cleanup_test_process (process, process_instance)
    deallocate (process_instance)
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_raw_2"

  end subroutine eio_raw_2


end module eio_raw_uti

