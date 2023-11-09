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

module eio_direct_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz, only: vector4_t
  use model_data, only: model_data_t
  use event_base
  use eio_data
  use eio_base

  use eio_direct

  use eio_base_ut, only: eio_prepare_test, eio_cleanup_test

  implicit none
  private

  public :: eio_direct_1

contains

  subroutine eio_direct_1 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    class(eio_t), allocatable :: eio
    type(event_sample_data_t) :: data
    type(string_t) :: sample
    type(vector4_t), dimension(:), allocatable :: p
    class(model_data_t), pointer :: model
    integer :: i, n_events, iostat, i_prc

    write (u, "(A)")  "* Test output: eio_direct_1"
    write (u, "(A)")  "*   Purpose: generate and read/write an event"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event, unweighted = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Initial state"
    write (u, "(A)")

    allocate (eio_direct_t :: eio)
    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Extract an empty event"
    write (u, "(A)")

    call eio%output (event, 1)
    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Retrieve contents"
    write (u, "(A)")

    select type (eio)
    class is (eio_direct_t)
       if (eio%has_event_index ())  write (u, "(A,1x,I0)")  "index =", eio%get_event_index ()
       if (eio%passed_known ())  write (u, "(A,1x,L1)")  "passed =", eio%has_passed ()
       write (u, "(A,1x,I0)")  "n_in =", eio%get_n_in ()
       write (u, "(A,1x,I0)")  "n_out =", eio%get_n_out ()
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Generate and extract an event"
    write (u, "(A)")

    call event%generate (1, [0._default, 0._default])
    call event%set_index (42)
    model => event%get_model_ptr ()

    sample = ""
    call eio%init_out (sample)
    call eio%output (event, 1, passed = .true.)
    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Retrieve contents"
    write (u, "(A)")

    select type (eio)
    class is (eio_direct_t)
       if (eio%has_event_index ())  write (u, "(A,1x,I0)")  "index =", eio%get_event_index ()
       if (eio%passed_known ())  write (u, "(A,1x,L1)")  "passed =", eio%has_passed ()
       write (u, "(A,1x,I0)")  "n_in =", eio%get_n_in ()
       write (u, "(A,1x,I0)")  "n_out =", eio%get_n_out ()
    end select

    select type (eio)
    class is (eio_direct_t)
       call eio%get_momentum_array (p)
       if (allocated (p)) then
          write (u, "(A)")  "p[3] ="
          call p(3)%write (u)
       end if
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Re-create an eio event record: initialization"
    write (u, "(A)")

    call eio%final ()

    select type (eio)
    class is (eio_direct_t)
       call eio%init_direct ( &
            n_beam = 0, n_in = 2, n_rem = 0, n_vir = 0, n_out = 2, &
            pdg = [25, 25, 25, 25], model = model)
       call eio%set_event_index (42)
       call eio%set_selection_indices (1, 1, 1, 1)
       call eio%write (u)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Re-create an eio event record: &
         &set momenta, interchanged"
    write (u, "(A)")

    select type (eio)
    class is (eio_direct_t)
       call eio%set_momentum (p([1,2,4,3]), on_shell=.true.)
       call eio%write (u)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* 'read' i_prc"
    write (u, "(A)")

    call eio%input_i_prc (i_prc, iostat)
    write (u, "(1x,A,1x,I0)")  "i_prc =", i_prc
    write (u, "(1x,A,1x,I0)")  "iostat =", iostat

    write (u, "(A)")
    write (u, "(A)")  "* 'read' (fill) event"
    write (u, "(A)")

    call eio%input_event (event, iostat)
    write (u, "(1x,A,1x,I0)")  "iostat =", iostat
    write (u, "(A)")

    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio%final ()
    deallocate (eio)

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_direct_1"

  end subroutine eio_direct_1


end module eio_direct_uti
