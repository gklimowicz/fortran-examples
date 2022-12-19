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

module dispatch_ut
  use unit_tests
  use dispatch_uti

  implicit none
  private

  public :: dispatch_test

  public :: dispatch_sf_data_test

contains

  subroutine dispatch_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (dispatch_1, "dispatch_1", &
         "process configuration method", &
         u, results)
    call test (dispatch_2, "dispatch_2", &
         "process core", &
         u, results)
    call test (dispatch_7, "dispatch_7", &
         "structure-function data", &
         u, results)
    call test (dispatch_8, "dispatch_8", &
         "beam structure", &
         u, results)
    call test (dispatch_10, "dispatch_10", &
         "process core update", &
         u, results)
    call test (dispatch_11, "dispatch_11", &
         "QCD coupling", &
         u, results)
  end subroutine dispatch_test


end module dispatch_ut
