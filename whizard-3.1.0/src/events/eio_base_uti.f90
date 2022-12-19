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

module eio_base_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use lorentz
  use model_data
  use particles
  use event_base
  use event_handles, only: event_handle_t
  use eio_data

  use eio_base

  implicit none
  private

  public :: eio_prepare_test, eio_cleanup_test
  public :: eio_prepare_fallback_model, eio_cleanup_fallback_model

  public :: eio_base_1

  abstract interface
     subroutine eio_prepare_event (event, unweighted, n_alt, sample_norm)
       import
       class(generic_event_t), intent(inout), pointer :: event
       logical, intent(in), optional :: unweighted
       integer, intent(in), optional :: n_alt
       type(string_t), intent(in), optional :: sample_norm
     end subroutine eio_prepare_event
  end interface

  abstract interface
     subroutine eio_cleanup_event (event)
       import
       class(generic_event_t), intent(inout), pointer :: event
     end subroutine eio_cleanup_event
  end interface

  procedure(eio_prepare_event), pointer :: eio_prepare_test => null ()
  procedure(eio_cleanup_event), pointer :: eio_cleanup_test => null ()

  abstract interface
     subroutine eio_prepare_model (model)
       import
       class(model_data_t), intent(inout), pointer :: model
     end subroutine eio_prepare_model
  end interface

  abstract interface
     subroutine eio_cleanup_model (model)
       import
       class(model_data_t), intent(inout), target :: model
     end subroutine eio_cleanup_model
  end interface

  type, extends (eio_t) :: eio_test_t
     integer :: event_n = 0
     integer :: event_i = 0
     integer :: i_prc = 0
     type(vector4_t), dimension(:,:), allocatable :: event_p
   contains
     procedure :: write => eio_test_write
     procedure :: final => eio_test_final
     procedure :: init_out => eio_test_init_out
     procedure :: init_in => eio_test_init_in
     procedure :: switch_inout => eio_test_switch_inout
     procedure :: output => eio_test_output
     procedure :: input_i_prc => eio_test_input_i_prc
     procedure :: input_event => eio_test_input_event
     procedure :: skip => eio_test_skip
  end type eio_test_t


  procedure(eio_prepare_model), pointer :: eio_prepare_fallback_model => null ()
  procedure(eio_cleanup_model), pointer :: eio_cleanup_fallback_model => null ()


contains

  subroutine eio_base_1 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    class(eio_t), allocatable :: eio
    integer :: i_prc,  iostat
    type(string_t) :: sample

    write (u, "(A)")  "* Test output: eio_base_1"
    write (u, "(A)")  "*   Purpose: generate and read/write an event"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event, unweighted = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_test1"

    allocate (eio_test_t :: eio)

    call eio%init_out (sample)
    call event%generate (1, [0._default, 0._default])
    call eio%output (event, 42)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Re-read the event"
    write (u, "(A)")

    call eio%init_in (sample)
    call eio%input_i_prc (i_prc, iostat)
    call eio%input_event (event, iostat)
    call eio%write (u)
    write (u, "(A)")
    write (u, "(1x,A,I0)")  "i = ", i_prc

    write (u, "(A)")
    write (u, "(A)")  "* Generate and append another event"
    write (u, "(A)")

    call eio%switch_inout ()
    call event%generate (1, [0._default, 0._default])
    call eio%output (event, 5)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Re-read both events"
    write (u, "(A)")

    call eio%init_in (sample)
    call eio%input_i_prc (i_prc, iostat)
    call eio%input_event (event, iostat)
    call eio%input_i_prc (i_prc, iostat)
    call eio%input_event (event, iostat)
    call eio%write (u)
    write (u, "(A)")
    write (u, "(1x,A,I0)")  "i = ", i_prc

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio%final ()
    deallocate (eio)

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_base_1"

  end subroutine eio_base_1


  subroutine eio_test_write (object, unit)
    class(eio_test_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Test event stream"
    if (object%event_i /= 0) then
       write (u, "(1x,A,I0,A)")  "Event #", object%event_i, ":"
       do i = 1, size (object%event_p, 1)
          call vector4_write (object%event_p(i, object%event_i), u)
       end do
    end if
  end subroutine eio_test_write

  subroutine eio_test_final (object)
    class(eio_test_t), intent(inout) :: object
    object%event_i = 0
  end subroutine eio_test_final

  subroutine eio_test_init_out (eio, sample, data, success, extension)
    class(eio_test_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(event_sample_data_t), intent(in), optional :: data
    logical, intent(out), optional :: success
    type(string_t), intent(in), optional :: extension
    eio%sample = sample
    eio%event_n = 0
    eio%event_i = 0
    allocate (eio%event_p (2, 10))
    if (present (success))  success = .true.
  end subroutine eio_test_init_out

  subroutine eio_test_init_in (eio, sample, data, success, extension)
    class(eio_test_t), intent(inout) :: eio
    type(string_t), intent(in) :: sample
    type(event_sample_data_t), intent(inout), optional :: data
    logical, intent(out), optional :: success
    type(string_t), intent(in), optional :: extension
    if (present (success))  success = .true.
  end subroutine eio_test_init_in

  subroutine eio_test_switch_inout (eio, success)
    class(eio_test_t), intent(inout) :: eio
    logical, intent(out), optional :: success
    if (present (success))  success = .true.
  end subroutine eio_test_switch_inout

  subroutine eio_test_output &
       (eio, event, i_prc, reading, passed, pacify, event_handle)
    class(eio_test_t), intent(inout) :: eio
    class(generic_event_t), intent(in), target :: event
    logical, intent(in), optional :: reading, passed, pacify
    class(event_handle_t), intent(inout), optional :: event_handle
    integer, intent(in) :: i_prc
    type(particle_set_t), pointer :: pset
    type(particle_t) :: prt
    eio%event_n = eio%event_n + 1
    eio%event_i = eio%event_n
    eio%i_prc = i_prc
    pset => event%get_particle_set_ptr ()
    prt = pset%get_particle (3)
    eio%event_p(1, eio%event_i) = prt%get_momentum ()
    prt = pset%get_particle (4)
    eio%event_p(2, eio%event_i) = prt%get_momentum ()
  end subroutine eio_test_output

  subroutine eio_test_input_i_prc (eio, i_prc, iostat)
    class(eio_test_t), intent(inout) :: eio
    integer, intent(out) :: i_prc
    integer, intent(out) :: iostat
    i_prc = eio%i_prc
    iostat = 0
  end subroutine eio_test_input_i_prc

  subroutine eio_test_input_event (eio, event, iostat, event_handle)
    class(eio_test_t), intent(inout) :: eio
    class(generic_event_t), intent(inout), target :: event
    integer, intent(out) :: iostat
    class(event_handle_t), intent(inout), optional :: event_handle
    eio%event_i = eio%event_i + 1
    iostat = 0
  end subroutine eio_test_input_event

  subroutine eio_test_skip (eio, iostat)
    class(eio_test_t), intent(inout) :: eio
    integer, intent(out) :: iostat
    iostat = 0
  end subroutine eio_test_skip


end module eio_base_uti
