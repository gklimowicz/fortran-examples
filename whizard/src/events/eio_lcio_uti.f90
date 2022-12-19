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

module eio_lcio_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use model_data
  use particles
  use event_base
  use eio_data
  use eio_base
  use hep_events
  use lcio_interface

  use eio_lcio

  use eio_base_ut, only: eio_prepare_test, eio_cleanup_test
  use eio_base_ut, only: eio_prepare_fallback_model, eio_cleanup_fallback_model

  implicit none
  private

  public :: eio_lcio_1
  public :: eio_lcio_2

contains

  subroutine eio_lcio_1 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(particle_set_t), pointer :: pset_ptr
    type(string_t) :: sample
    integer :: u_file, iostat
    character(215) :: buffer

    write (u, "(A)")  "* Test output: eio_lcio_1"
    write (u, "(A)")  "*   Purpose: write a LCIO file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event)

    call data%init (1)
    data%n_beam = 2
    data%unweighted = .true.
    data%norm_mode = NORM_UNIT
    data%pdg_beam = 25
    data%energy_beam = 500
    data%proc_num_id = [42]
    data%cross_section(1) = 100
    data%error(1) = 1
    data%total_cross_section = sum (data%cross_section)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_lcio_1"

    allocate (eio_lcio_t :: eio)
    select type (eio)
    type is (eio_lcio_t)
       call eio%set_parameters ()
    end select

    call eio%init_out (sample, data)

    call event%generate (1, [0._default, 0._default])
    call event%set_index (77)
    call event%pacify_particle_set ()

    call eio%output (event, i_prc = 1)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Reset data"
    write (u, "(A)")

    deallocate (eio)
    allocate (eio_lcio_t :: eio)

    select type (eio)
    type is (eio_lcio_t)
       call eio%set_parameters ()
    end select
    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write LCIO file contents to ASCII file"
    write (u, "(A)")

    select type (eio)
    type is (eio_lcio_t)
       call lcio_event_init (eio%lcio_event, &
            proc_id = 42, &
            event_id = event%get_index ())
       pset_ptr => event%get_particle_set_ptr ()
       call lcio_event_from_particle_set &
            (eio%lcio_event,  pset_ptr)
       call write_lcio_event (eio%lcio_event, var_str ("test_file.slcio"))
       call lcio_event_final (eio%lcio_event, .true.)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Read in ASCII contents of LCIO file"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = "test_file.slcio", &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (iostat /= 0)  exit
       if (trim (buffer) == "")  cycle
       if (buffer(1:12) == " - timestamp")  cycle
       if (buffer(1:6) == " date:")  cycle
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_lcio_1"

  end subroutine eio_lcio_1

  subroutine eio_lcio_2 (u)
    integer, intent(in) :: u
    class(model_data_t), pointer :: fallback_model
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: iostat, i_prc

    write (u, "(A)")  "* Test output: eio_lcio_2"
    write (u, "(A)")  "*   Purpose: read a LCIO event"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    allocate (fallback_model)
    call eio_prepare_fallback_model (fallback_model)
    call eio_prepare_test (event)

    call data%init (1)
    data%n_beam = 2
    data%unweighted = .true.
    data%norm_mode = NORM_UNIT
    data%pdg_beam = 25
    data%energy_beam = 500
    data%proc_num_id = [42]
    data%cross_section(1) = 100
    data%error(1) = 1
    data%total_cross_section = sum (data%cross_section)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_lcio_2"

    allocate (eio_lcio_t :: eio)
    select type (eio)
    type is (eio_lcio_t)
       call eio%set_parameters (recover_beams = .false.)
    end select
    call eio%set_fallback_model (fallback_model)

    call eio%init_out (sample, data)
    call event%generate (1, [0._default, 0._default])
    call event%set_index (88)
    call event%evaluate_expressions ()
    call event%pacify_particle_set ()

    call eio%output (event, i_prc = 1)
    call eio%write (u)
    call eio%final ()
    deallocate (eio)

    call event%reset_contents ()
    call event%reset_index ()

    write (u, "(A)")
    write (u, "(A)")  "* Initialize"
    write (u, "(A)")

    allocate (eio_lcio_t :: eio)
    select type (eio)
    type is (eio_lcio_t)
       call eio%set_parameters (recover_beams = .false.)
    end select
    call eio%set_fallback_model (fallback_model)

    call data%init (1)
    data%n_beam = 2
    data%unweighted = .true.
    data%norm_mode = NORM_UNIT
    data%pdg_beam = 25
    data%energy_beam = 500
    data%proc_num_id = [42]
    call data%write (u)
    write (u, *)

    write (u, "(A)")  "* Initialize"
    write (u, "(A)")

    call eio%init_in (sample, data)
    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Read event"
    write (u, "(A)")

    call eio%input_i_prc (i_prc, iostat)

    select type (eio)
    type is (eio_lcio_t)
       write (u, "(A,I0,A,I0)")  "Found process #", i_prc, &
            " with ID = ", eio%proc_num_id(i_prc)
    end select

    call eio%input_event (event, iostat)
    call event%select (1, 1, 1)
    call event%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Read closing"
    write (u, "(A)")

    call eio%input_i_prc (i_prc, iostat)
    write (u, "(A,I0)")  "iostat = ", iostat

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio%final ()

    call eio_cleanup_test (event)
    call eio_cleanup_fallback_model (fallback_model)
    deallocate (fallback_model)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_lcio_2"

  end subroutine eio_lcio_2


end module eio_lcio_uti
