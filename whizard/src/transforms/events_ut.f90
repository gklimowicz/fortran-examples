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

module events_ut
  use unit_tests
  use events_uti

  implicit none
  private

  public :: events_test

contains

  subroutine events_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (events_1, "events_1", &
         "empty event record", &
         u, results)
    call test (events_2, "events_2", &
         "generate event", &
         u, results)
    call test (events_4, "events_4", &
         "recover event", &
         u, results)
    call test (events_5, "events_5", &
         "partially recover event", &
         u, results)
    call test (events_6, "events_6", &
         "decays", &
         u, results)
    call test (events_7, "events_7", &
         "decay options", &
         u, results)
  end subroutine events_test


end module events_ut
