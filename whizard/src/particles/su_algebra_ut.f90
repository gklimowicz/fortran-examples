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

module su_algebra_ut
  use unit_tests
  use su_algebra_uti

  implicit none
  private

  public :: su_algebra_test

contains

  subroutine su_algebra_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (su_algebra_1, "su_algebra_1", &
         "generator ordering", &
         u, results)
    call test (su_algebra_2, "su_algebra_2", &
         "Cartan generator representation", &
         u, results)
    call test (su_algebra_3, "su_algebra_3", &
         "Cartan generator mapping", &
         u, results)
    call test (su_algebra_4, "su_algebra_4", &
         "Root-helicity mapping", &
         u, results)
  end subroutine su_algebra_test


end module su_algebra_ut
