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

module integration_results_ut
  use unit_tests
  use integration_results_uti

  implicit none
  private

  public :: integration_results_test

contains

  subroutine integration_results_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
      call test (integration_results_1, "integration_results_1", &
           "record single line and write to log", &
           u, results)
      call test (integration_results_2, "integration_results_2", &
           "record single result and write to log", &
           u, results)
      call test (integration_results_3, "integration_results_3", &
           "initialize display and add/display each entry", &
           u, results)
      call test (integration_results_4, "integration_results_4", &
           "record extended results and display", &
           u, results)
      call test (integration_results_5, "integration_results_5", &
           "record extended results and display", &
           u, results)
  end subroutine integration_results_test


end module integration_results_ut
