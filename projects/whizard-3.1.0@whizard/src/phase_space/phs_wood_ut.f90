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

module phs_wood_ut
  use unit_tests
  use phs_wood_uti

  implicit none
  private

  public :: phs_wood_test
  public :: phs_wood_vis_test

  public :: write_test_phs_file

contains

  subroutine phs_wood_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (phs_wood_1, "phs_wood_1", &
         "phase-space configuration", &
         u, results)
    call test (phs_wood_2, "phs_wood_2", &
         "phase-space evaluation", &
         u, results)
    call test (phs_wood_3, "phs_wood_3", &
         "phase-space generation", &
         u, results)
    call test (phs_wood_4, "phs_wood_4", &
         "nontrivial process", &
         u, results)
    call test (phs_wood_5, "phs_wood_5", &
         "equivalences", &
         u, results)
    call test (phs_wood_6, "phs_wood_6", &
         "phase-space generation", &
         u, results)
  end subroutine phs_wood_test

  subroutine phs_wood_vis_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (phs_wood_vis_1, "phs_wood_vis_1", &
         "visualizing phase space channels", &
         u, results)
  end subroutine phs_wood_vis_test


end module phs_wood_ut
