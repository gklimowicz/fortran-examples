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

module api_hepmc_uti

  use iso_fortran_env, only: int32, real64 !NODEP!
  use hepmc_interface, only: hepmc_event_t

  use api

  implicit none
  private

  public :: api_hepmc_1

contains

  subroutine api_hepmc_1 (u)
    use hepmc_interface, only: hepmc_event_get_event_index
    use hepmc_interface, only: hepmc_event_get_n_particles
    use hepmc_interface, only: hepmc_event_print
    use hepmc_interface, only: hepmc_event_final
    integer, intent(in) :: u

    type(whizard_api_t) :: whizard
    type(simulation_api_t) :: sample

    integer :: it_begin, it_end
    integer :: i
    integer(int32) :: idx, npt

    type(hepmc_event_t) :: hepmc_event

    write (u, "(A)")  "* Test output: api_hepmc_1"
    write (u, "(A)")  "*   Purpose:  generate events"
    write (u, "(A)")

    call whizard%option ("model", "QED")
    call whizard%option ("library", "api_hepmc_1_lib")
    call whizard%option ("logfile", "api_hepmc_1_log.out")
    call whizard%option ("rebuild", "T")
    call whizard%init ()

    call whizard%command ("process api_hepmc_1_p = e1, E1 => e2, E2")

    call whizard%set_var ("sqrts", 10._real64)
    call whizard%command ("iterations = 1:100")
    call whizard%set_var ("seed", 0_int32)
    call whizard%command ("integrate (api_hepmc_1_p)")

    call whizard%set_var ("?unweighted", .false.)
    call whizard%set_var ("$sample", "api_hepmc_1_evt")
    call whizard%command ("sample_format = hepmc")
    call whizard%set_var ("n_events", 2_int32)

    call whizard%new_sample ("api_hepmc_1_p", sample)
    call sample%open (it_begin, it_end)
    do i = it_begin, it_end
       call sample%next_event (hepmc_event)
       idx = hepmc_event_get_event_index (hepmc_event)
       npt = hepmc_event_get_n_particles (hepmc_event)
       call hepmc_event_print (hepmc_event)
       write (u, "(A,I0)")  "Event #", idx
       write (u, "(2x,A,1x,I0)")  "n_particles =", npt
       call hepmc_event_final (hepmc_event)
    end do
    call sample%close ()

    call whizard%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: api_hepmc_1"

  end subroutine api_hepmc_1


end module api_hepmc_uti
