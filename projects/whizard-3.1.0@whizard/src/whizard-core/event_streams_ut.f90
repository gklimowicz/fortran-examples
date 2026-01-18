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

module event_streams_ut
  use unit_tests
  use event_streams_uti

  implicit none
  private

  public :: event_streams_test

contains

  subroutine event_streams_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (event_streams_1, "event_streams_1", &
         "empty event stream array", &
         u, results)
    call test (event_streams_2, "event_streams_2", &
         "nontrivial event stream array", &
         u, results)
    call test (event_streams_3, "event_streams_3", &
         "switch input/output", &
         u, results)
    call test (event_streams_4, "event_streams_4", &
         "check MD5 sum", &
         u, results)
  end subroutine event_streams_test


end module event_streams_ut
