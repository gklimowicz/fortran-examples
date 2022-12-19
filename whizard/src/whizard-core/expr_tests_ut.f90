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
module expr_tests_ut

  use unit_tests
  use expr_tests_uti

  implicit none
  private

  public :: subevt_expr_test

contains

  subroutine subevt_expr_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (subevt_expr_1, "subevt_expr_1", &
         "parton-event expressions", &
         u, results)
    call test (subevt_expr_2, "subevt_expr_2", &
         "parton-event expressions", &
         u, results)
    call test (processes_5, "processes_5", &
         "handle cuts (partonic event)", &
         u, results)
    call test (processes_6, "processes_6", &
         "handle scales and weight (partonic event)", &
         u, results)
    call test (events_3, "events_3", &
         "expression evaluation", &
         u, results)
end subroutine subevt_expr_test


end module expr_tests_ut
