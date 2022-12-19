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

module api_uti

  use iso_fortran_env, only: int32, real64 !NODEP!
  use diagnostics, only: msg_message

  use api

  implicit none
  private

  public :: api_1
  public :: api_2
  public :: api_3
  public :: api_4
  public :: api_5
  public :: api_6
  public :: api_7
  public :: api_8

contains

  subroutine api_1 (u)
    integer, intent(in) :: u
    type(whizard_api_t) :: whizard

    character(:), allocatable :: logfile
    integer :: u_log
    integer :: iostat
    character(80) :: buffer

    write (u, "(A)")  "* Test output: api_1"
    write (u, "(A)")  "*   Purpose:  call init/final"
    write (u, "(A)")

    logfile = "api_1_log.out"
    call whizard%option ("logfile", logfile)
    call whizard%init ()
    call msg_message ("Intentional error: double init")
    call whizard%init ()
    call whizard%final ()
    call msg_message ("Intentional error: double final")
    call whizard%final ()

    open (newunit = u_log, file = logfile, action = "read", status = "old")
    do
       read (u_log, "(A)", iostat=iostat)  buffer
       if (iostat /= 0)  exit
       if (buffer(1:10) == "| WHIZARD:")  write (u, "(A)")  trim (buffer)
       if (buffer(1:10) == "*** ERROR:")  write (u, "(A)")  trim (buffer)
    end do
    close (u_log)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: api_1"

  end subroutine api_1

  subroutine api_2 (u)
    integer, intent(in) :: u
    type(whizard_api_t) :: whizard

    character(:), allocatable :: job_id
    character(:), allocatable :: sample
    logical :: unweighted
    integer :: n_events
    real(real64) :: sqrts
    logical :: known

    write (u, "(A)")  "* Test output: api_2"
    write (u, "(A)")  "*   Purpose:  access Sindarin variables"
    write (u, "(A)")

    call whizard%option ("logfile", "api_2_log.out")
    call whizard%option ("job_id", "api_2_ID")
    call whizard%init ()

    call whizard%get_var ("sqrts", sqrts, known)
    write (u, "(A,1x,L1)")  "sqrts is known =", known
    write (u, "(A,1x,F5.1)")  "sqrts =", sqrts

    call whizard%set_var ("sqrts", 100._real64)
    call whizard%set_var ("n_events", 3_int32)
    call whizard%set_var ("?unweighted", .false.)
    call whizard%set_var ("$sample", "foo")

    call whizard%get_var ("sqrts", sqrts, known)
    call whizard%get_var ("$job_id", job_id)
    call whizard%get_var ("n_events", n_events)
    call whizard%get_var ("?unweighted", unweighted)
    call whizard%get_var ("$sample", sample)

    write (u, "(A,1x,L1)")  "sqrts is known =", known
    write (u, "(A,1x,F5.1)")  "sqrts =", sqrts
    write (u, "(A,1x,A)")  "$job_id =", job_id
    write (u, "(A,1x,I0)")  "n_events =", n_events
    write (u, "(A,1x,L1)")  "?unweighted =", unweighted
    write (u, "(A,1x,A)")  "$sample =", sample

    write (u, *)

    call whizard%set_var ("?unweighted", .true.)
    call whizard%get_var ("?unweighted", unweighted)
    write (u, "(A,1x,L1)")  "?unweighted =", unweighted

    call whizard%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: api_2"

  end subroutine api_2

  subroutine api_3 (u)
    integer, intent(in) :: u
    type(whizard_api_t) :: whizard

    character(:), allocatable :: model_name

    write (u, "(A)")  "* Test output: api_3"
    write (u, "(A)")  "*   Purpose:  set model in advance and via Sindarin string"
    write (u, "(A)")

    call whizard%option ("model", "QCD")
    call whizard%option ("logfile", "api_3_log.out")
    call whizard%init ()

    call whizard%get_var ("$model_name", model_name)
    write (u, "(A,1x,A)")  "model =", model_name

    call whizard%command ("model = QED")

    call whizard%get_var ("$model_name", model_name)
    write (u, "(A,1x,A)")  "model =", model_name

    call whizard%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: api_3"

  end subroutine api_3

  subroutine api_4 (u)
    integer, intent(in) :: u
    type(whizard_api_t) :: whizard

    character(:), allocatable :: flv_string

    write (u, "(A)")  "* Test output: api_4"
    write (u, "(A)")  "*   Purpose:  translate PDG code(s) to flavor string"
    write (u, "(A)")

    call whizard%option ("model", "QED")
    call whizard%option ("logfile", "api_4_log.out")
    call whizard%init ()

    write (u, "(A,1x,A)")  "electron =", whizard%flv_string (11)
    write (u, "(A,1x,A)")  "leptons  =", whizard%flv_string ([11,13,15])

    call whizard%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: api_4"

  end subroutine api_4

  subroutine api_5 (u)
    integer, intent(in) :: u
    type(whizard_api_t) :: whizard

    real(real64) :: sqrts, integral, error
    logical :: known

    write (u, "(A)")  "* Test output: api_5"
    write (u, "(A)")  "*   Purpose:  integrate and retrieve results"
    write (u, "(A)")

    call whizard%option ("model", "QED")
    call whizard%option ("library", "api_5_lib")
    call whizard%option ("logfile", "api_5_log.out")
    call whizard%option ("rebuild", "T")
    call whizard%init ()

    write (u, "(A)")  "* Process setup"
    write (u, "(A)")

    call whizard%command ("process api_5_p = e1, E1 => e2, E2")
    call whizard%command ("sqrts = 10")
    call whizard%command ("iterations = 1:100")
    call whizard%set_var ("seed", 0_int32)
    call whizard%get_integration_result ("api_5_p", integral, error, known)
    write (u, 2)  "integral is known =", known

    call whizard%command ("integrate (api_5_p)")

    write (u, "(A)")
    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call whizard%get_integration_result ("api_5_p", integral, error, known)
    write (u, 2)  "integral is known =", known

    call whizard%get_var ("sqrts", sqrts)
    call whizard%get_integration_result ("api_5_p", integral, error)
    write (u, 1)  "sqrt(s)       =", sqrts, "GeV"
    write (u, 1)  "cross section =", integral / 1000, "pb"
    write (u, 1)  "error         =", error / 1000, "pb"
1   format (2x,A,1x,F5.1,1x,A)
2   format (2x,A,1x,L1)

    call whizard%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: api_5"

  end subroutine api_5

  subroutine api_6 (u)
    integer, intent(in) :: u

    type(whizard_api_t) :: whizard
    type(simulation_api_t) :: sample

    integer :: it_begin, it_end
    integer :: i
    integer(int32) :: idx
    real(real64) :: sqme
    real(real64) :: weight

    write (u, "(A)")  "* Test output: api_6"
    write (u, "(A)")  "*   Purpose:  generate events"
    write (u, "(A)")

    call whizard%option ("model", "QED")
    call whizard%option ("library", "api_6_lib")
    call whizard%option ("logfile", "api_6_log.out")
    call whizard%option ("rebuild", "T")
    call whizard%init ()

    call whizard%command ("process api_6_p = e1, E1 => e2, E2")

    call whizard%set_var ("sqrts", 10._real64)
    call whizard%command ("iterations = 1:100")
    call whizard%set_var ("seed", 0_int32)
    call whizard%command ("integrate (api_6_p)")

    call whizard%set_var ("?unweighted", .false.)
    call whizard%set_var ("$sample", "api_6_evt")
    call whizard%command ("sample_format = dump")
    call whizard%set_var ("n_events", 2_int32)
    call whizard%set_var ("event_index_offset", 4_int32)

    call whizard%new_sample ("api_6_p", sample)
    call sample%open (it_begin, it_end)
    do i = it_begin, it_end
       call sample%next_event ()
       call sample%get_event_index (idx)
       call sample%get_weight (weight)
       call sample%get_sqme (sqme)
       write (u, "(A,I0)")  "Event #", idx
       write (u, 3)  "sqme    =", sqme
       write (u, 3)  "weight  =", weight
3      format (2x,A,1x,ES10.3)
    end do
    call sample%close ()

    call whizard%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: api_6"

  end subroutine api_6

  subroutine api_7 (u)
    integer, intent(in) :: u

    type(whizard_api_t) :: whizard
    type(simulation_api_t) :: sample

    integer :: it_begin, it_end
    integer :: i
    integer(int32) :: idx
    integer(int32) :: i_proc
    character(:), allocatable :: proc_id
    real(real64) :: sqrts
    real(real64) :: scale
    real(real64) :: alpha_s
    real(real64) :: sqme
    real(real64) :: weight

    write (u, "(A)")  "* Test output: api_7"
    write (u, "(A)")  "*   Purpose:  generate events"
    write (u, "(A)")

    call whizard%option ("model", "QCD")
    call whizard%option ("library", "api_7_lib")
    call whizard%option ("logfile", "api_7_log.out")
    call whizard%option ("rebuild", "T")
    call whizard%init ()

    call whizard%command ("process api_7_p1 = u, U => t, T")
    call whizard%command ("process api_7_p2 = d, D => t, T")
    call whizard%command ("process api_7_p3 = s, S => t, T")

    call whizard%set_var ("sqrts", 1000._real64)
    call whizard%command ("beams = p, p => pdf_builtin")
    call whizard%set_var ("?alphas_is_fixed", .false.)
    call whizard%set_var ("?alphas_from_pdf_builtin", .true.)
    call whizard%command ("iterations = 1:100")
    call whizard%set_var ("seed", 0_int32)
    call whizard%command ("integrate (api_7_p1)")
    call whizard%command ("integrate (api_7_p2)")
    call whizard%command ("integrate (api_7_p3)")

    call whizard%set_var ("?unweighted", .false.)
    call whizard%set_var ("$sample", "api_7_evt")
    call whizard%command ("sample_format = dump")
    call whizard%set_var ("n_events", 10_int32)

    call whizard%new_sample ("api_7_p1, api_7_p2 ,api_7_p3", sample)
    call sample%open (it_begin, it_end)
    do i = it_begin, it_end
       call sample%next_event ()
       call sample%get_event_index (idx)
       call sample%get_process_index (i_proc)
       call sample%get_process_id (proc_id)
       call sample%get_sqrts (sqrts)
       call sample%get_fac_scale (scale)
       call sample%get_alpha_s (alpha_s)
       call sample%get_weight (weight)
       call sample%get_sqme (sqme)
       write (u, "(A,I0)")  "Event #", idx
       write (u, 1)  "process #", i_proc
       write (u, 2)  "proc_id =", proc_id
       write (u, 3)  "sqrts   =", sqrts
       write (u, 3)  "f_scale =", scale
       write (u, 3)  "alpha_s =", alpha_s
       write (u, 3)  "sqme    =", sqme
       write (u, 3)  "weight  =", weight
1      format (2x,A,I0)
2      format (2x,A,1x,A)
3      format (2x,A,1x,ES10.3)
    end do
    call sample%close ()

    call whizard%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: api_7"

  end subroutine api_7

  subroutine api_8 (u)
    integer, intent(in) :: u
    type(whizard_api_t) :: whizard

    character(:), allocatable :: logfile
    integer :: u_log
    integer :: iostat
    character(80) :: buffer

    character(:), allocatable :: pack_args
    character(:), allocatable :: unpack_args

    write (u, "(A)")  "* Test output: api_8"
    write (u, "(A)")  "*   Purpose:  check pack/unpack options"
    write (u, "(A)")

    logfile = "api_8_log.out"
    call whizard%option ("logfile", logfile)
    call whizard%option ("pack", "api_8_foo")
    call whizard%option ("unpack", "api_8_bar.tgz, api_8_gee.tgz")
    call whizard%init ()
    call msg_message ("WHIZARD: Intentional errors: pack/unpack files do not exist")

    call whizard%final ()

    open (newunit = u_log, file = logfile, action = "read", status = "old")
    do
       read (u_log, "(A)", iostat=iostat)  buffer
       if (iostat /= 0)  exit
       if (buffer(1:10) == "| WHIZARD:")  write (u, "(A)")  trim (buffer)
       if (buffer(1:10) == "*** ERROR:")  write (u, "(A)")  trim (buffer)
    end do
    close (u_log)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: api_8"

  end subroutine api_8


end module api_uti
