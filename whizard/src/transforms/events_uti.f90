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

module events_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use model_data
  use particles
  use process_libraries
  use process_stacks
  use event_transforms
  use decays
  use decays_ut, only: prepare_testbed

  use process, only: process_t
  use instances, only: process_instance_t

  use events

  implicit none
  private

  public :: events_1
  public :: events_2
  public :: events_4
  public :: events_5
  public :: events_6
  public :: events_7

contains

  subroutine events_1 (u)
    integer, intent(in) :: u
    type(event_t), target :: event

    write (u, "(A)")  "* Test output: events_1"
    write (u, "(A)")  "*   Purpose: display an empty event object"
    write (u, "(A)")

    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: events_1"

  end subroutine events_1

  subroutine events_2 (u)
    use processes_ut, only: prepare_test_process, cleanup_test_process
    integer, intent(in) :: u
    type(event_t), allocatable, target :: event
    type(process_t), allocatable, target :: process
    type(process_instance_t), allocatable, target :: process_instance
    type(model_data_t), target :: model

    write (u, "(A)")  "* Test output: events_2"
    write (u, "(A)")  "*   Purpose: generate and display an event"
    write (u, "(A)")

    call model%init_test ()

    write (u, "(A)")  "* Generate test process event"

    allocate (process)
    allocate (process_instance)
    call prepare_test_process (process, process_instance, model)
    call process_instance%setup_event_data ()

    write (u, "(A)")
    write (u, "(A)")  "* Initialize event object"

    allocate (event)
    call event%basic_init ()
    call event%connect (process_instance, process%get_model_ptr ())

    write (u, "(A)")
    write (u, "(A)")  "* Generate test process event"

    call process_instance%generate_weighted_event (1)

    write (u, "(A)")
    write (u, "(A)")  "* Fill event object"
    write (u, "(A)")

    call event%generate (1, [0.4_default, 0.4_default])
    call event%increment_index ()
    call event%evaluate_expressions ()
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call event%final ()
    deallocate (event)

    call cleanup_test_process (process, process_instance)
    deallocate (process_instance)
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: events_2"

  end subroutine events_2

  subroutine events_4 (u)
    use processes_ut, only: prepare_test_process, cleanup_test_process
    integer, intent(in) :: u
    type(event_t), allocatable, target :: event
    type(process_t), allocatable, target :: process
    type(process_instance_t), allocatable, target :: process_instance
    type(process_t), allocatable, target :: process2
    type(process_instance_t), allocatable, target :: process2_instance
    type(particle_set_t) :: particle_set
    type(model_data_t), target :: model

    write (u, "(A)")  "* Test output: events_4"
    write (u, "(A)")  "*   Purpose: generate and recover an event"
    write (u, "(A)")

    call model%init_test ()

    write (u, "(A)")  "* Generate test process event and save particle set"
    write (u, "(A)")

    allocate (process)
    allocate (process_instance)
    call prepare_test_process (process, process_instance, model)
    call process_instance%setup_event_data ()

    allocate (event)
    call event%basic_init ()
    call event%connect (process_instance, process%get_model_ptr ())

    call event%generate (1, [0.4_default, 0.4_default])
    call event%increment_index ()
    call event%evaluate_expressions ()
    call event%write (u)

    particle_set = event%get_particle_set_ptr ()
    ! NB: 'particle_set' contains pointers to the model within 'process'

    call event%final ()
    deallocate (event)

    write (u, "(A)")
    write (u, "(A)")  "* Recover event from particle set"
    write (u, "(A)")

    allocate (process2)
    allocate (process2_instance)
    call prepare_test_process (process2, process2_instance, model)
    call process2_instance%setup_event_data ()

    allocate (event)
    call event%basic_init ()
    call event%connect (process2_instance, process2%get_model_ptr ())

    call event%select (1, 1, 1)
    call event%set_hard_particle_set (particle_set)
    call event%recalculate (update_sqme = .true.)
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Transfer sqme and evaluate expressions"
    write (u, "(A)")

    call event%accept_sqme_prc ()
    call event%accept_weight_prc ()
    call event%check ()
    call event%set_index (1)
    call event%evaluate_expressions ()
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Reset contents"
    write (u, "(A)")

    call event%reset_contents ()
    call event%reset_index ()
    event%transform_first%particle_set_exists = .false.
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call particle_set%final ()

    call event%final ()
    deallocate (event)

    call cleanup_test_process (process2, process2_instance)
    deallocate (process2_instance)
    deallocate (process2)

    call cleanup_test_process (process, process_instance)
    deallocate (process_instance)
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: events_4"

  end subroutine events_4

  subroutine events_5 (u)
    use processes_ut, only: prepare_test_process, cleanup_test_process
    integer, intent(in) :: u
    type(event_t), allocatable, target :: event
    type(process_t), allocatable, target :: process
    type(process_instance_t), allocatable, target :: process_instance
    type(process_t), allocatable, target :: process2
    type(process_instance_t), allocatable, target :: process2_instance
    type(particle_set_t) :: particle_set
    real(default) :: sqme, weight
    type(model_data_t), target :: model

    write (u, "(A)")  "* Test output: events_5"
    write (u, "(A)")  "*   Purpose: generate and recover an event"
    write (u, "(A)")

    call model%init_test ()

    write (u, "(A)")  "* Generate test process event and save particle set"
    write (u, "(A)")

    allocate (process)
    allocate (process_instance)
    call prepare_test_process (process, process_instance, model)
    call process_instance%setup_event_data ()

    allocate (event)
    call event%basic_init ()
    call event%connect (process_instance, process%get_model_ptr ())

    call event%generate (1, [0.4_default, 0.4_default])
    call event%increment_index ()
    call event%evaluate_expressions ()
    call event%write (u)

    particle_set = event%get_particle_set_ptr ()
    sqme = event%get_sqme_ref ()
    weight = event%get_weight_ref ()

    call event%final ()
    deallocate (event)

    write (u, "(A)")
    write (u, "(A)")  "* Recover event from particle set"
    write (u, "(A)")

    allocate (process2)
    allocate (process2_instance)
    call prepare_test_process (process2, process2_instance, model)
    call process2_instance%setup_event_data ()

    allocate (event)
    call event%basic_init ()
    call event%connect (process2_instance, process2%get_model_ptr ())

    call event%select (1, 1, 1)
    call event%set_hard_particle_set (particle_set)
    call event%recalculate (update_sqme = .false.)
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually set sqme and evaluate expressions"
    write (u, "(A)")

    call event%set (sqme_ref = sqme, weight_ref = weight)
    call event%accept_sqme_ref ()
    call event%accept_weight_ref ()
    call event%set_index (1)
    call event%evaluate_expressions ()
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call particle_set%final ()

    call event%final ()
    deallocate (event)

    call cleanup_test_process (process2, process2_instance)
    deallocate (process2_instance)
    deallocate (process2)

    call cleanup_test_process (process, process_instance)
    deallocate (process_instance)
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: events_5"

  end subroutine events_5

  subroutine events_6 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(string_t) :: prefix, procname1, procname2
    type(process_library_t), target :: lib
    type(process_stack_t) :: process_stack
    class(evt_t), pointer :: evt_decay
    type(event_t), allocatable, target :: event
    type(process_t), pointer :: process
    type(process_instance_t), allocatable, target :: process_instance

    write (u, "(A)")  "* Test output: events_6"
    write (u, "(A)")  "*   Purpose: generate an event with subsequent decays"
    write (u, "(A)")

    write (u, "(A)")  "* Generate test process and decay"
    write (u, "(A)")

    call os_data%init ()

    prefix = "events_6"
    procname1 = prefix // "_p"
    procname2 = prefix // "_d"
    call prepare_testbed &
         (lib, process_stack, prefix, os_data, &
         scattering=.true., decay=.true.)

    write (u, "(A)")  "* Initialize decay process"

    process => process_stack%get_process_ptr (procname1)
    model => process%get_model_ptr ()
    call model%set_unstable (25, [procname2])

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data ()
    call process_instance%init_simulation (1)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize event transform: decay"

    allocate (evt_decay_t :: evt_decay)
    call evt_decay%connect (process_instance, model, process_stack)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize event object"
    write (u, "(A)")

    allocate (event)
    call event%basic_init ()
    call event%connect (process_instance, model)
    call event%import_transform (evt_decay)

    call event%write (u, show_decay = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Generate event"
    write (u, "(A)")

    call event%generate (1, [0.4_default, 0.4_default])
    call event%increment_index ()
    call event%evaluate_expressions ()
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call event%final ()
    deallocate (event)

    call process_instance%final ()
    deallocate (process_instance)

    call process_stack%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: events_6"

  end subroutine events_6

  subroutine events_7 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(string_t) :: prefix, procname2
    type(process_library_t), target :: lib
    type(process_stack_t) :: process_stack
    type(process_t), pointer :: process
    type(process_instance_t), allocatable, target :: process_instance

    write (u, "(A)")  "* Test output: events_7"
    write (u, "(A)")  "*   Purpose: check decay options"
    write (u, "(A)")

    write (u, "(A)")  "* Prepare test process"
    write (u, "(A)")

    call os_data%init ()

    prefix = "events_7"
    procname2 = prefix // "_d"
    call prepare_testbed &
         (lib, process_stack, prefix, os_data, &
         scattering=.false., decay=.true.)

    write (u, "(A)")  "* Generate decay event, default options"
    write (u, "(A)")

    process => process_stack%get_process_ptr (procname2)
    model => process%get_model_ptr ()
    call model%set_unstable (25, [procname2])

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data (model)
    call process_instance%init_simulation (1)

    call process_instance%generate_weighted_event (1)
    call process_instance%write (u)

    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Generate decay event, helicity-diagonal decay"
    write (u, "(A)")

    process => process_stack%get_process_ptr (procname2)
    model => process%get_model_ptr ()
    call model%set_unstable (25, [procname2], diagonal = .true.)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data (model)
    call process_instance%init_simulation (1)

    call process_instance%generate_weighted_event (1)
    call process_instance%write (u)

    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Generate decay event, isotropic decay, &
         &polarized final state"
    write (u, "(A)")

    process => process_stack%get_process_ptr (procname2)
    model => process%get_model_ptr ()
    call model%set_unstable (25, [procname2], isotropic = .true.)
    call model%set_polarized (6)
    call model%set_polarized (-6)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data (model)
    call process_instance%init_simulation (1)

    call process_instance%generate_weighted_event (1)
    call process_instance%write (u)

    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process_stack%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: events_7"

  end subroutine events_7


end module events_uti

