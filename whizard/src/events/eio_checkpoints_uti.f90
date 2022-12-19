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

module eio_checkpoints_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use event_base
  use eio_data
  use eio_base

  use eio_checkpoints

  use eio_base_ut, only: eio_prepare_test, eio_cleanup_test

  implicit none
  private

  public :: eio_checkpoints_1

contains

  subroutine eio_checkpoints_1 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    class(eio_t), allocatable :: eio
    type(event_sample_data_t) :: data
    type(string_t) :: sample
    integer :: i, n_events

    write (u, "(A)")  "* Test output: eio_checkpoints_1"
    write (u, "(A)")  "*   Purpose: generate a number of events &
         &with screen output"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Generate events"
    write (u, "(A)")

    sample = "eio_checkpoints_1"

    allocate (eio_checkpoints_t :: eio)

    n_events = 10
    call data%init (1, 0)
    data%n_evt = n_events

    select type (eio)
    type is (eio_checkpoints_t)
       call eio%set_parameters (checkpoint = 4)
    end select

    call eio%init_out (sample, data)

    do i = 1, n_events
       call event%generate (1, [0._default, 0._default])
       call eio%output (event, i_prc = 0)
    end do

    write (u, "(A)")  "* Checkpointing status"
    write (u, "(A)")

    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_checkpoints_1"

  end subroutine eio_checkpoints_1


end module eio_checkpoints_uti
