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

module mci_vamp2_ut
  use unit_tests
  use mci_vamp2_uti

  implicit none
  private

  public :: mci_vamp2_test

contains

  subroutine mci_vamp2_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (mci_vamp2_1, "mci_vamp2_1", "one-dimensional integral", u, results)
    call test (mci_vamp2_2, "mci_vamp2_2", &
         "multiple iterations", &
         u, results)
    call test (mci_vamp2_3, "mci_vamp2_3", &
         "grid adaptation", &
         u, results)
  end subroutine mci_vamp2_test


end module mci_vamp2_ut
