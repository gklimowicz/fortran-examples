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

module api_ut
  use unit_tests
  use api_uti

  implicit none
  private

  public :: api_test

contains

  subroutine api_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (api_1, "api_1", &
         "init/final", &
         u, results)
    call test (api_2, "api_2", &
         "set/get Sindarin values", &
         u, results)
    call test (api_3, "api_3", &
         "preload model and execute command (model)", &
         u, results)
    call test (api_4, "api_4", &
         "flavor string translation", &
         u, results)
    call test (api_5, "api_5", &
         "integration", &
         u, results)
    call test (api_6, "api_6", &
         "event generation", &
         u, results)
    call test (api_7, "api_7", &
         "more event generation", &
         u, results)
    call test (api_8, "api_8", &
         "check pack/unpack options", &
         u, results)
  end subroutine api_test


end module api_ut
