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

module eio_lhef_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use model_data
  use event_base
  use eio_data
  use eio_base

  use eio_lhef

  use eio_base_ut, only: eio_prepare_test, eio_cleanup_test
  use eio_base_ut, only: eio_prepare_fallback_model, eio_cleanup_fallback_model

  implicit none
  private

  public :: eio_lhef_1
  public :: eio_lhef_2
  public :: eio_lhef_3
  public :: eio_lhef_4
  public :: eio_lhef_5
  public :: eio_lhef_6

contains

  subroutine eio_lhef_1 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, iostat
    character(80) :: buffer

    write (u, "(A)")  "* Test output: eio_lhef_1"
    write (u, "(A)")  "*   Purpose: generate an event and write weight to file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event, unweighted = .false.)

    call data%init (1)
    data%n_evt = 1
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

    sample = "eio_lhef_1"

    allocate (eio_lhef_t :: eio)
    select type (eio)
    type is (eio_lhef_t)
       call eio%set_parameters ()
    end select

    call eio%init_out (sample, data)
    call event%generate (1, [0._default, 0._default])

    call eio%output (event, i_prc = 1)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* File contents:"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = char (sample // "." // eio%extension), &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (buffer(1:21) == "  <generator_version>")  buffer = "[...]"
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Reset data"
    write (u, "(A)")

    deallocate (eio)
    allocate (eio_lhef_t :: eio)
    select type (eio)
    type is (eio_lhef_t)
       call eio%set_parameters ()
    end select

    select type (eio)
    type is (eio_lhef_t)
       call eio%set_parameters (keep_beams = .true.)
    end select
    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_lhef_1"

  end subroutine eio_lhef_1

  subroutine eio_lhef_2 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, iostat
    character(80) :: buffer

    write (u, "(A)")  "* Test output: eio_lhef_2"
    write (u, "(A)")  "*   Purpose: generate an event and write weight to file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event, unweighted = .false.)

    call data%init (1)
    data%unweighted = .false.
    data%norm_mode = NORM_SIGMA
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

    sample = "eio_lhef_2"

    allocate (eio_lhef_t :: eio)
    select type (eio)
    type is (eio_lhef_t)
       call eio%set_parameters (version = "2.0", write_sqme_prc = .true.)
    end select

    call eio%init_out (sample, data)
    call event%generate (1, [0._default, 0._default])


    call eio%output (event, i_prc = 1)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* File contents:"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = char (sample // "." // eio%extension), &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (buffer(1:10) == "<generator")  buffer = "[...]"
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_lhef_2"

  end subroutine eio_lhef_2

  subroutine eio_lhef_3 (u)
    integer, intent(in) :: u
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, iostat
    character(80) :: buffer

    write (u, "(A)")  "* Test output: eio_lhef_3"
    write (u, "(A)")  "*   Purpose: generate an event and write weight to file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize test process"

    call eio_prepare_test (event, unweighted = .false.)

    call data%init (1)
    data%unweighted = .false.
    data%norm_mode = NORM_SIGMA
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

    sample = "eio_lhef_3"

    allocate (eio_lhef_t :: eio)
    select type (eio)
    type is (eio_lhef_t)
       call eio%set_parameters (version = "3.0", write_sqme_prc = .true.)
    end select

    call eio%init_out (sample, data)
    call event%generate (1, [0._default, 0._default])


    call eio%output (event, i_prc = 1)
    call eio%write (u)
    call eio%final ()

    write (u, "(A)")
    write (u, "(A)")  "* File contents:"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = char (sample // ".lhe"), &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (buffer(1:10) == "<generator")  buffer = "[...]"
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio_cleanup_test (event)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_lhef_3"

  end subroutine eio_lhef_3

  subroutine eio_lhef_4 (u)
    integer, intent(in) :: u
    class(model_data_t), pointer :: fallback_model
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, iostat, i_prc

    write (u, "(A)")  "* Test output: eio_lhef_4"
    write (u, "(A)")  "*   Purpose: read a LHEF 1.0 file"
    write (u, "(A)")

    write (u, "(A)")  "* Write a LHEF data file"
    write (u, "(A)")

    u_file = free_unit ()
    sample = "eio_lhef_4"
    open (u_file, file = char (sample // ".lhe"), &
         status = "replace", action = "readwrite")

    write (u_file, "(A)")  '<LesHouchesEvents version="1.0">'
    write (u_file, "(A)")  '<header>'
    write (u_file, "(A)")  '  <arbitrary_tag opt="foo">content</arbitrary_tag>'
    write (u_file, "(A)")  '  Text'
    write (u_file, "(A)")  '  <another_tag />'
    write (u_file, "(A)")  '</header>'
    write (u_file, "(A)")  '<init>'
    write (u_file, "(A)")  ' 25 25  5.0000000000E+02  5.0000000000E+02 &
         & -1 -1 -1 -1 3 1'
    write (u_file, "(A)")  '  1.0000000000E-01  1.0000000000E-03 &
         & 1.0000000000E+00 42'
    write (u_file, "(A)")  '</init>'
    write (u_file, "(A)")  '<event>'
    write (u_file, "(A)")  ' 4 42  3.0574068604E+08  1.0000000000E+03 &
         & -1.0000000000E+00 -1.0000000000E+00'
    write (u_file, "(A)")  ' 25 -1 0 0 0 0  0.0000000000E+00  0.0000000000E+00 &
         & 4.8412291828E+02  5.0000000000E+02  1.2500000000E+02 &
         & 0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  ' 25 -1 0 0 0 0  0.0000000000E+00  0.0000000000E+00 &
         &-4.8412291828E+02  5.0000000000E+02  1.2500000000E+02 &
         & 0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  ' 25 1 1 2 0 0 -1.4960220911E+02 -4.6042825611E+02 &
         & 0.0000000000E+00  5.0000000000E+02  1.2500000000E+02 &
         & 0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  ' 25 1 1 2 0 0  1.4960220911E+02  4.6042825611E+02 &
         & 0.0000000000E+00  5.0000000000E+02  1.2500000000E+02 &
         & 0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  '</event>'
    write (u_file, "(A)")  '</LesHouchesEvents>'
    close (u_file)


    write (u, "(A)")  "* Initialize test process"
    write (u, "(A)")

    allocate (fallback_model)
    call eio_prepare_fallback_model (fallback_model)
    call eio_prepare_test (event, unweighted = .false.)

    allocate (eio_lhef_t :: eio)
    select type (eio)
    type is (eio_lhef_t)
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

    write (u, "(A)")  "* Initialize and read header"
    write (u, "(A)")

    call eio%init_in (sample, data)
    call eio%write (u)

    write (u, *)

    select type (eio)
    type is (eio_lhef_t)
       call eio%tag_lhef%write (u);  write (u, *)
    end select

    write (u, *)
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Read event"
    write (u, "(A)")

    call eio%input_i_prc (i_prc, iostat)

    select type (eio)
    type is (eio_lhef_t)
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
    write (u, "(A)")  "* Test output end: eio_lhef_4"

  end subroutine eio_lhef_4

  subroutine eio_lhef_5 (u)
    integer, intent(in) :: u
    class(model_data_t), pointer :: fallback_model
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, iostat, i_prc

    write (u, "(A)")  "* Test output: eio_lhef_5"
    write (u, "(A)")  "*   Purpose: read a LHEF 2.0 file"
    write (u, "(A)")

    write (u, "(A)")  "* Write a LHEF data file"
    write (u, "(A)")

    u_file = free_unit ()
    sample = "eio_lhef_5"
    open (u_file, file = char (sample // ".lhe"), &
         status = "replace", action = "readwrite")

    write (u_file, "(A)")  '<LesHouchesEvents version="2.0">'
    write (u_file, "(A)")  '<header>'
    write (u_file, "(A)")  '</header>'
    write (u_file, "(A)")  '<init>'
    write (u_file, "(A)")  ' 25 25  5.0000000000E+02  5.0000000000E+02 &
         &-1 -1 -1 -1 4 1'
    write (u_file, "(A)")  '  1.0000000000E-01  1.0000000000E-03 &
         & 0.0000000000E+00 42'
    write (u_file, "(A)")  '<generator version="2.2.3">WHIZARD&
         &</generator>'
    write (u_file, "(A)")  '<xsecinfo neve="1" totxsec="1.0000000000E-01" />'
    write (u_file, "(A)")  '</init>'
    write (u_file, "(A)")  '<event>'
    write (u_file, "(A)")  ' 4 42  3.0574068604E+08  1.0000000000E+03 &
         &-1.0000000000E+00 -1.0000000000E+00'
    write (u_file, "(A)")  ' 25 -1 0 0 0 0  0.0000000000E+00 &
         & 0.0000000000E+00  4.8412291828E+02  5.0000000000E+02 &
         & 1.2500000000E+02  0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  ' 25 -1 0 0 0 0  0.0000000000E+00 &
         & 0.0000000000E+00 -4.8412291828E+02  5.0000000000E+02 &
         & 1.2500000000E+02  0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  ' 25 1 1 2 0 0 -1.4960220911E+02 &
         &-4.6042825611E+02  0.0000000000E+00  5.0000000000E+02 &
         & 1.2500000000E+02  0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  ' 25 1 1 2 0 0  1.4960220911E+02 &
         & 4.6042825611E+02  0.0000000000E+00  5.0000000000E+02 &
         & 1.2500000000E+02  0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  '<weight name="sqme_prc">1.0000000000E+00</weight>'
    write (u_file, "(A)")  '</event>'
    write (u_file, "(A)")  '</LesHouchesEvents>'
    close (u_file)

    write (u, "(A)")  "* Initialize test process"
    write (u, "(A)")

    allocate (fallback_model)
    call eio_prepare_fallback_model (fallback_model)
    call eio_prepare_test (event, unweighted = .false.)

    allocate (eio_lhef_t :: eio)
    select type (eio)
    type is (eio_lhef_t)
       call eio%set_parameters (version = "2.0", recover_beams = .false.)
    end select
    call eio%set_fallback_model (fallback_model)

    call data%init (1)
    data%unweighted = .false.
    data%norm_mode = NORM_SIGMA
    data%n_beam = 2
    data%pdg_beam = 25
    data%energy_beam = 500
    data%proc_num_id = [42]
    call data%write (u)
    write (u, *)

    write (u, "(A)")  "* Initialize and read header"
    write (u, "(A)")

    call eio%init_in (sample, data)
    call eio%write (u)

    write (u, *)

    select type (eio)
    type is (eio_lhef_t)
       call eio%tag_lhef%write (u);  write (u, *)
    end select

    write (u, *)
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Read event"
    write (u, "(A)")

    call eio%input_i_prc (i_prc, iostat)

    select type (eio)
    type is (eio_lhef_t)
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
    write (u, "(A)")  "* Test output end: eio_lhef_5"

  end subroutine eio_lhef_5

  subroutine eio_lhef_6 (u)
    integer, intent(in) :: u
    class(model_data_t), pointer :: fallback_model
    class(generic_event_t), pointer :: event
    type(event_sample_data_t) :: data
    class(eio_t), allocatable :: eio
    type(string_t) :: sample
    integer :: u_file, iostat, i_prc

    write (u, "(A)")  "* Test output: eio_lhef_6"
    write (u, "(A)")  "*   Purpose: read a LHEF 3.0 file"
    write (u, "(A)")

    write (u, "(A)")  "* Write a LHEF data file"
    write (u, "(A)")

    u_file = free_unit ()
    sample = "eio_lhef_6"
    open (u_file, file = char (sample // ".lhe"), &
         status = "replace", action = "readwrite")

    write (u_file, "(A)")  '<LesHouchesEvents version="3.0">'
    write (u_file, "(A)")  '<header>'
    write (u_file, "(A)")  '</header>'
    write (u_file, "(A)")  '<init>'
    write (u_file, "(A)")  ' 25 25  5.0000000000E+02  5.0000000000E+02 &
         &-1 -1 -1 -1 4 1'
    write (u_file, "(A)")  '  1.0000000000E-01  1.0000000000E-03 &
         & 0.0000000000E+00 42'
    write (u_file, "(A)")  '<generator version="2.2.3">WHIZARD&
         &</generator>'
    write (u_file, "(A)")  '<xsecinfo neve="1" totxsec="1.0000000000E-01" />'
    write (u_file, "(A)")  '<weightinfo name="sqme_prc" />'
    write (u_file, "(A)")  '</init>'
    write (u_file, "(A)")  '<event>'
    write (u_file, "(A)")  ' 4 42  3.0574068604E+08  1.0000000000E+03 &
         &-1.0000000000E+00 -1.0000000000E+00'
    write (u_file, "(A)")  ' 25 -1 0 0 0 0  0.0000000000E+00 &
         & 0.0000000000E+00  4.8412291828E+02  5.0000000000E+02 &
         & 1.2500000000E+02  0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  ' 25 -1 0 0 0 0  0.0000000000E+00 &
         & 0.0000000000E+00 -4.8412291828E+02  5.0000000000E+02 &
         & 1.2500000000E+02  0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  ' 25 1 1 2 0 0 -1.4960220911E+02 &
         &-4.6042825611E+02  0.0000000000E+00  5.0000000000E+02 &
         & 1.2500000000E+02  0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  ' 25 1 1 2 0 0  1.4960220911E+02 &
         & 4.6042825611E+02  0.0000000000E+00  5.0000000000E+02 &
         & 1.2500000000E+02  0.0000000000E+00  9.0000000000E+00'
    write (u_file, "(A)")  '<weights>1.0000000000E+00</weights>'
    write (u_file, "(A)")  '</event>'
    write (u_file, "(A)")  '</LesHouchesEvents>'
    close (u_file)

    write (u, "(A)")  "* Initialize test process"
    write (u, "(A)")

    allocate (fallback_model)
    call eio_prepare_fallback_model (fallback_model)
    call eio_prepare_test (event, unweighted = .false.)

    allocate (eio_lhef_t :: eio)
    select type (eio)
    type is (eio_lhef_t)
       call eio%set_parameters (version = "3.0", recover_beams = .false.)
    end select
    call eio%set_fallback_model (fallback_model)

    call data%init (1)
    data%unweighted = .false.
    data%norm_mode = NORM_SIGMA
    data%n_beam = 2
    data%pdg_beam = 25
    data%energy_beam = 500
    data%proc_num_id = [42]
    call data%write (u)
    write (u, *)

    write (u, "(A)")  "* Initialize and read header"
    write (u, "(A)")

    call eio%init_in (sample, data)
    call eio%write (u)

    write (u, *)

    select type (eio)
    type is (eio_lhef_t)
       call eio%tag_lhef%write (u);  write (u, *)
    end select

    write (u, *)
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Read event"
    write (u, "(A)")

    call eio%input_i_prc (i_prc, iostat)

    select type (eio)
    type is (eio_lhef_t)
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
    write (u, "(A)")  "* Test output end: eio_lhef_6"

  end subroutine eio_lhef_6


end module eio_lhef_uti
