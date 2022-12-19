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

module eio_dump_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use event_base
  use eio_data
  use eio_base

  use eio_dump

  use eio_base_ut, only: eio_prepare_test, eio_cleanup_test

  implicit none
  private

  public :: eio_dump_1

contains

  subroutine eio_dump_1 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    class(eio_t), allocatable :: eio
    integer :: i_prc
    integer :: u_file

    write (u, "(A)")  "* Test output: eio_dump_1"
    write (u, "(A)")  "*   Purpose: generate events and write essentials to output"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event, unweighted = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write three events (two passed)"
    write (u, "(A)")

    allocate (eio_dump_t :: eio)
    select type (eio)
    type is (eio_dump_t)
       call eio%set_parameters (unit = u, weights = .true., pacify = .true.)
    end select

    i_prc = 42

    call eio%init_out (var_str (""))

    call event%generate (1, [0._default, 0._default])
    call eio%output (event, i_prc = i_prc)

    call event%generate (1, [0.1_default, 0._default])
    call event%set_index (99)
    call eio%output (event, i_prc = i_prc, passed = .false.)

    call event%generate (1, [0.2_default, 0._default])
    call event%increment_index ()
    call eio%output (event, i_prc = i_prc, passed = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of eio_dump object"
    write (u, "(A)")

    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    select type (eio)
    type is (eio_dump_t)
       eio%writing = .false.
    end select
    call eio%final ()

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_dump_1"
  end subroutine eio_dump_1


end module eio_dump_uti
