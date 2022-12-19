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

module eio_weights_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use event_base
  use eio_data
  use eio_base

  use eio_weights

  use eio_base_ut, only: eio_prepare_test, eio_cleanup_test

  implicit none
  private

  public :: eio_weights_1
  public :: eio_weights_2
  public :: eio_weights_3

contains

  subroutine eio_weights_1 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file
    character(80) :: buffer

    write (u, "(A)")  "* Test output: eio_weights_1"
    write (u, "(A)")  "*   Purpose: generate an event and write weight to file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event, unweighted = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_weights_1"

    allocate (eio_weights_t :: eio)

    call eio%init_out (sample)
    call event%generate (1, [0._default, 0._default])

    call eio%output (event, i_prc = 42)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* File contents: &
         &(weight, sqme(evt), sqme(prc), i_prc)"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = "eio_weights_1.weights.dat", &
         action = "read", status = "old")
    read (u_file, "(A)")  buffer
    write (u, "(A)") trim (buffer)
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_weights_1"
  end subroutine eio_weights_1

  subroutine eio_weights_2 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, i
    character(80) :: buffer

    write (u, "(A)")  "* Test output: eio_weights_2"
    write (u, "(A)")  "*   Purpose: generate an event and write weight to file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event, unweighted = .false., n_alt = 2)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_weights_2"

    allocate (eio_weights_t :: eio)

    call eio%init_out (sample)
    select type (eio)
    type is (eio_weights_t)
       call eio%set_parameters (pacify = .true.)
    end select
    call event%generate (1, [0._default, 0._default])
    call event%set (sqme_alt = [2._default, 3._default])
    call event%set (weight_alt = &
         [2 * event%get_weight_prc (), 3 * event%get_weight_prc ()])

    call eio%output (event, i_prc = 42)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* File contents: &
         &(weight, sqme(evt), sqme(prc), i_prc)"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = "eio_weights_2.weights.dat", &
         action = "read", status = "old")
    do i = 1, 3
       read (u_file, "(A)")  buffer
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_weights_2"

  end subroutine eio_weights_2

  subroutine eio_weights_3 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, iostat
    character(80) :: buffer

    write (u, "(A)")  "* Test output: eio_weights_3"
    write (u, "(A)")  "*   Purpose: generate three events and write to file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event, unweighted = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write events"
    write (u, "(A)")

    sample = "eio_weights_3"

    allocate (eio_weights_t :: eio)
    select type (eio)
    type is (eio_weights_t)
       call eio%set_parameters (pacify = .true.)
    end select

    call eio%init_out (sample)

    call event%generate (1, [0._default, 0._default])
    call eio%output (event, i_prc = 1)

    call event%generate (1, [0.1_default, 0._default])
    call eio%output (event, i_prc = 1, passed = .false.)

    call event%generate (1, [0.2_default, 0._default])
    call eio%output (event, i_prc = 1, passed = .true.)

    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* File contents: &
         &(weight, sqme(evt), sqme(prc), i_prc), should be just two entries"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = "eio_weights_3.weights.dat", &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat=iostat)  buffer
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_weights_3"
  end subroutine eio_weights_3


end module eio_weights_uti
