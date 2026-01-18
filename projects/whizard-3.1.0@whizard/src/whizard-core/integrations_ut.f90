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

module integrations_ut
  use unit_tests
  use integrations_uti

  implicit none
  private

  public :: integrations_test
  public :: integrations_history_test

contains

  subroutine integrations_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (integrations_1, "integrations_1", &
         "intrinsic test process", &
         u, results)
    call test (integrations_2, "integrations_2", &
         "intrinsic test process with cut", &
         u, results)
    call test (integrations_3, "integrations_3", &
         "standard phase space", &
         u, results)
    call test (integrations_4, "integrations_4", &
         "VAMP integration (one iteration)", &
         u, results)
    call test (integrations_5, "integrations_5", &
         "VAMP integration (three iterations)", &
         u, results)
    call test (integrations_6, "integrations_6", &
         "VAMP integration (three passes)", &
         u, results)
    call test (integrations_7, "integrations_7", &
         "VAMP integration with wood phase space", &
         u, results)
    call test (integrations_8, "integrations_8", &
         "integration with structure function", &
         u, results)
    call test (integrations_9, "integrations_9", &
         "handle sign change", &
         u, results)
  end subroutine integrations_test

  subroutine integrations_history_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (integrations_history_1, "integrations_history_1", &
         "Test integration history files", &
         u, results)
  end subroutine integrations_history_test


end module integrations_ut
