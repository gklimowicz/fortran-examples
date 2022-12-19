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

module mci_vamp_ut
  use unit_tests
  use mci_vamp_uti

  implicit none
  private

  public :: mci_vamp_test

contains

  subroutine mci_vamp_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (mci_vamp_1, "mci_vamp_1", &
         "one-dimensional integral", &
         u, results)
    call test (mci_vamp_2, "mci_vamp_2", &
         "multiple iterations", &
         u, results)
    call test (mci_vamp_3, "mci_vamp_3", &
         "grid adaptation", &
         u, results)
    call test (mci_vamp_4, "mci_vamp_4", &
         "two-dimensional integration", &
         u, results)
    call test (mci_vamp_5, "mci_vamp_5", &
         "two-dimensional integration", &
         u, results)
    call test (mci_vamp_6, "mci_vamp_6", &
         "weight adaptation", &
         u, results)
    call test (mci_vamp_7, "mci_vamp_7", &
         "use channel equivalences", &
         u, results)
    call test (mci_vamp_8, "mci_vamp_8", &
         "integration passes", &
         u, results)
    call test (mci_vamp_9, "mci_vamp_9", &
         "weighted event", &
         u, results)
    call test (mci_vamp_10, "mci_vamp_10", &
         "grids I/O", &
         u, results)
   call test (mci_vamp_11, "mci_vamp_11", &
        "weighted events with grid I/O", &
        u, results)
   call test (mci_vamp_12, "mci_vamp_12", &
        "unweighted events with grid I/O", &
        u, results)
    call test (mci_vamp_13, "mci_vamp_13", &
         "updating integration results", &
         u, results)
    call test (mci_vamp_14, "mci_vamp_14", &
         "accuracy goal", &
         u, results)
    call test (mci_vamp_15, "mci_vamp_15", &
         "VAMP history", &
         u, results)
    call test (mci_vamp_16, "mci_vamp_16", &
         "1-D integral with sign change", &
         u, results)
  end subroutine mci_vamp_test


end module mci_vamp_ut
