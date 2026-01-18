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

module eio_stdhep_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use model_data
  use event_base
  use eio_data
  use eio_base
  use xdr_wo_stdhep

  use eio_stdhep

  use eio_base_ut, only: eio_prepare_test, eio_cleanup_test
  use eio_base_ut, only: eio_prepare_fallback_model, eio_cleanup_fallback_model

  implicit none
  private

  public :: eio_stdhep_1
  public :: eio_stdhep_2
  public :: eio_stdhep_3
  public :: eio_stdhep_4

contains

  subroutine eio_stdhep_1 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, iostat
    character(215) :: buffer

    write (u, "(A)")  "* Test output: eio_stdhep_1"
    write (u, "(A)")  "*   Purpose: generate an event in STDHEP HEPEVT format"
    write (u, "(A)")  "*      and write weight to file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event)

    call data%init (1)
    data%n_evt = 1
    data%n_beam = 2
    data%pdg_beam = 25
    data%energy_beam = 500
    data%proc_num_id = [42]
    data%cross_section(1) = 100
    data%error(1) = 1
    data%total_cross_section = sum (data%cross_section)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_stdhep_1"

    allocate (eio_stdhep_hepevt_t :: eio)
    select type (eio)
    type is (eio_stdhep_hepevt_t)
       call eio%set_parameters ()
    end select

    call eio%init_out (sample, data)
    call event%generate (1, [0._default, 0._default])
    call event%set_index (61)  ! not supported by reader, actually
    call event%evaluate_expressions ()
    call event%pacify_particle_set ()

    call eio%output (event, i_prc = 1)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Write STDHEP file contents to ASCII file"
    write (u, "(A)")

    call write_stdhep_event &
         (sample // ".hep", var_str ("eio_stdhep_1.hep.out"), 1)

    write (u, "(A)")
    write (u, "(A)")  "* Read in ASCII contents of STDHEP file"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = "eio_stdhep_1.hep.out", &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (iostat /= 0)  exit
       if (trim (buffer) == "")  cycle
       if (buffer(1:18) == "    total blocks: ")  &
            buffer = "    total blocks: [...]"
       if (buffer(1:25) == "           title: WHIZARD")  &
            buffer = "           title: WHIZARD [version]"
       if (buffer(1:17) == "            date:")  &
            buffer = "            date: [...]"
       if (buffer(1:17) == "    closing date:")  &
            buffer = "    closing date: [...]"
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Reset data"
    write (u, "(A)")

    deallocate (eio)
    allocate (eio_stdhep_hepevt_t :: eio)

    select type (eio)
    type is (eio_stdhep_hepevt_t)
       call eio%set_parameters (keep_beams = .true.)
    end select
    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_stdhep_1"

  end subroutine eio_stdhep_1

  subroutine eio_stdhep_2 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(model_data_t), pointer :: fallback_model
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, iostat
    character(215) :: buffer

    write (u, "(A)")  "* Test output: eio_stdhep_2"
    write (u, "(A)")  "*   Purpose: generate an event in STDHEP HEPEUP format"
    write (u, "(A)")  "*      and write weight to file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    allocate (fallback_model)
    call eio_prepare_fallback_model (fallback_model)
    call eio_prepare_test (event, unweighted = .false.)

    call data%init (1)
    data%n_evt = 1
    data%n_beam = 2
    data%pdg_beam = 25
    data%energy_beam = 500
    data%proc_num_id = [42]
    data%cross_section(1) = 100
    data%error(1) = 1
    data%total_cross_section = sum (data%cross_section)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_stdhep_2"

    allocate (eio_stdhep_hepeup_t :: eio)
    select type (eio)
    type is (eio_stdhep_hepeup_t)
       call eio%set_parameters ()
    end select
    call eio%set_fallback_model (fallback_model)

    call eio%init_out (sample, data)
    call event%generate (1, [0._default, 0._default])
    call event%set_index (62)  ! not supported by reader, actually
    call event%evaluate_expressions ()

    call eio%output (event, i_prc = 1)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Write STDHEP file contents to ASCII file"
    write (u, "(A)")

    call write_stdhep_event &
         (sample // ".up.hep", var_str ("eio_stdhep_2.hep.out"), 2)

    write (u, "(A)")
    write (u, "(A)")  "* Read in ASCII contents of STDHEP file"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = "eio_stdhep_2.hep.out", &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (iostat /= 0)  exit
       if (trim (buffer) == "")  cycle
       if (buffer(1:18) == "    total blocks: ")  &
            buffer = "    total blocks: [...]"
       if (buffer(1:25) == "           title: WHIZARD")  &
            buffer = "           title: WHIZARD [version]"
       if (buffer(1:17) == "            date:")  &
            buffer = "            date: [...]"
       if (buffer(1:17) == "    closing date:")  &
            buffer = "    closing date: [...]"
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Reset data"
    write (u, "(A)")

    deallocate (eio)
    allocate (eio_stdhep_hepeup_t :: eio)

    select type (eio)
    type is (eio_stdhep_hepeup_t)
       call eio%set_parameters (keep_beams = .true.)
    end select
    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)
    call eio_cleanup_fallback_model (fallback_model)
    deallocate (fallback_model)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_stdhep_2"

  end subroutine eio_stdhep_2

  subroutine eio_stdhep_3 (u)
    integer, intent(in) :: u
    class(model_data_t), pointer :: fallback_model
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: iostat, i_prc

    write (u, "(A)")  "* Test output: eio_stdhep_3"
    write (u, "(A)")  "*   Purpose: read a StdHep file, HEPEVT block"
    write (u, "(A)")

    write (u, "(A)")  "* Write a StdHep data file, HEPEVT block"
    write (u, "(A)")

    allocate (fallback_model)
    call eio_prepare_fallback_model (fallback_model)
    call eio_prepare_test (event)

    call data%init (1)
    data%n_evt = 1
    data%n_beam = 2
    data%pdg_beam = 25
    data%energy_beam = 500
    data%proc_num_id = [42]
    data%cross_section(1) = 100
    data%error(1) = 1
    data%total_cross_section = sum (data%cross_section)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event"
    write (u, "(A)")

    sample = "eio_stdhep_3"

    allocate (eio_stdhep_hepevt_t :: eio)
    select type (eio)
    type is (eio_stdhep_hepevt_t)
       call eio%set_parameters ()
    end select
    call eio%set_fallback_model (fallback_model)

    call eio%init_out (sample, data)
    call event%generate (1, [0._default, 0._default])
    call event%set_index (63)  ! not supported by reader, actually
    call event%evaluate_expressions ()

    call eio%output (event, i_prc = 1)
    call eio%write (u)
    call eio%final ()

    call eio_cleanup_test (event)
    call eio_cleanup_fallback_model (fallback_model)
    deallocate (eio)
    deallocate (fallback_model)

    write (u, "(A)")  "* Initialize test process"
    write (u, "(A)")

    allocate (fallback_model)
    call eio_prepare_fallback_model (fallback_model)
    call eio_prepare_test (event, unweighted = .false.)

    allocate (eio_stdhep_hepevt_t :: eio)
    select type (eio)
    type is (eio_stdhep_hepevt_t)
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
    type is (eio_stdhep_hepevt_t)
       write (u, "(A,I0,A,I0)")  "Found process #", i_prc, &
            " with ID = ", eio%proc_num_id(i_prc)
    end select

    call eio%input_event (event, iostat)

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
    write (u, "(A)")  "* Test output end: eio_stdhep_3"

  end subroutine eio_stdhep_3

  subroutine eio_stdhep_4 (u)
    integer, intent(in) :: u
    class(model_data_t), pointer :: fallback_model
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: iostat, i_prc

    write (u, "(A)")  "* Test output: eio_stdhep_3"
    write (u, "(A)")  "*   Purpose: read a StdHep file, HEPRUP/HEPEUP block"
    write (u, "(A)")

    write (u, "(A)")  "* Write a StdHep data file, HEPRUP/HEPEUP block"
    write (u, "(A)")

    allocate (fallback_model)
    call eio_prepare_fallback_model (fallback_model)
    call eio_prepare_test (event)

    call data%init (1)
    data%n_evt = 1
    data%n_beam = 2
    data%pdg_beam = 25
    data%energy_beam = 500
    data%proc_num_id = [42]
    data%cross_section(1) = 100
    data%error(1) = 1
    data%total_cross_section = sum (data%cross_section)

    write (u, "(A)")
    write (u, "(A)")  "* Generate and write an event, HEPEUP/HEPRUP"
    write (u, "(A)")

    sample = "eio_stdhep_4"

    allocate (eio_stdhep_hepeup_t :: eio)
    select type (eio)
    type is (eio_stdhep_hepeup_t)
       call eio%set_parameters ()
    end select
    call eio%set_fallback_model (fallback_model)

    call eio%init_out (sample, data)
    call event%generate (1, [0._default, 0._default])
    call event%set_index (64)   ! not supported by reader, actually
    call event%evaluate_expressions ()
    call event%pacify_particle_set ()

    call eio%output (event, i_prc = 1)
    call eio%write (u)
    call eio%final ()

    call eio_cleanup_test (event)
    call eio_cleanup_fallback_model (fallback_model)
    deallocate (eio)
    deallocate (fallback_model)

    write (u, "(A)")  "* Initialize test process"
    write (u, "(A)")

    allocate (fallback_model)
    call eio_prepare_fallback_model (fallback_model)
    call eio_prepare_test (event, unweighted = .false.)

    allocate (eio_stdhep_hepeup_t :: eio)
    select type (eio)
    type is (eio_stdhep_hepeup_t)
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
    type is (eio_stdhep_hepeup_t)
       write (u, "(A,I0,A,I0)")  "Found process #", i_prc, &
            " with ID = ", eio%proc_num_id(i_prc)
    end select

    call eio%input_event (event, iostat)

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
    write (u, "(A)")  "* Test output end: eio_stdhep_4"

  end subroutine eio_stdhep_4


end module eio_stdhep_uti
