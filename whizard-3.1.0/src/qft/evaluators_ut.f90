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

module evaluators_ut
  use unit_tests
  use evaluators_uti

  implicit none
  private

  public :: evaluator_test

contains

  subroutine evaluator_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (evaluator_1, "evaluator_1", &
         "check evaluators (1)", &
         u, results)
    call test (evaluator_2, "evaluator_2", &
         "check evaluators (2)", &
         u, results)
    call test (evaluator_3, "evaluator_3", &
         "check evaluators (3)", &
         u, results)
    call test (evaluator_4, "evaluator_4", &
         "check evaluator product with filter", &
         u, results)
  end subroutine evaluator_test


end module evaluators_ut
