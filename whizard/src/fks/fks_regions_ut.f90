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

module fks_regions_ut
  use unit_tests
  use fks_regions_uti

  implicit none
  private

  public :: fks_regions_test

contains

  subroutine fks_regions_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test(fks_regions_1, "fks_regions_1", &
         "Test flavor structure utilities", u, results)
    call test(fks_regions_2, "fks_regions_2", &
         "Test singular regions for final-state radiation for n = 2", &
         u, results)
    call test(fks_regions_3, "fks_regions_3", &
         "Test singular regions for final-state radiation for n = 3", &
         u, results)
    call test(fks_regions_4, "fks_regions_4", &
         "Test singular regions for final-state radiation for n = 4", &
         u, results)
    call test(fks_regions_5, "fks_regions_5", &
         "Test singular regions for final-state radiation for n = 5", &
         u, results)
    call test(fks_regions_6, "fks_regions_6", &
         "Test singular regions for initial-state radiation", &
         u, results)
    call test(fks_regions_7, "fks_regions_7", &
         "Check Latex output", u, results)
    call test(fks_regions_8, "fks_regions_8", &
         "Test singular regions for initial-state photon contributions", &
         u, results)
  end subroutine fks_regions_test


end module fks_regions_ut
