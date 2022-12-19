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

module compilations_ut
  use unit_tests
  use compilations_uti

  implicit none
  private

  public :: compilations_test
  public :: compilations_static_test

contains

  subroutine compilations_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (compilations_1, "compilations_1", &
         "intrinsic test processes", &
         u, results)
    call test (compilations_2, "compilations_2", &
         "external process (omega)", &
         u, results)
    call test (compilations_3, "compilations_3", &
         "static executable: driver", &
         u, results)
end subroutine compilations_test

  subroutine compilations_static_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (compilations_static_1, "compilations_static_1", &
         "static executable: compilation", &
         u, results)
    call test (compilations_static_2, "compilations_static_2", &
         "static executable: shortcut", &
         u, results)
end subroutine compilations_static_test


end module compilations_ut
