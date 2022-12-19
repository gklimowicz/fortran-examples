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

module ttv_formfactors_ut
  use unit_tests
  use ttv_formfactors_uti

  implicit none
  private

  public ::ttv_formfactors_test

contains

  subroutine ttv_formfactors_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test(ttv_formfactors_1, "ttv_formfactors_1", &
              "Basic setup", u, results)
    call test(ttv_formfactors_2, "ttv_formfactors_2", &
              "Test flags", u, results)
  end subroutine ttv_formfactors_test


end module ttv_formfactors_ut
