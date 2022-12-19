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

module resonance_insertion_ut
  use unit_tests
  use resonance_insertion_uti

  implicit none
  private

  public :: resonance_insertion_test

contains

  subroutine resonance_insertion_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (resonance_insertion_1, "resonance_insertion_1", &
         "simple resonance insertion", &
         u, results)
    call test (resonance_insertion_2, "resonance_insertion_2", &
         "resonance color mismatch", &
         u, results)
    call test (resonance_insertion_3, "resonance_insertion_3", &
         "complex resonance history", &
         u, results)
    call test (resonance_insertion_4, "resonance_insertion_4", &
         "resonance history selection", &
         u, results)
    call test (resonance_insertion_5, "resonance_insertion_5", &
         "resonance history on/off", &
         u, results)
    call test (resonance_insertion_6, "resonance_insertion_6", &
         "resonance insertion with beam structure", &
         u, results)
  end subroutine resonance_insertion_test


end module resonance_insertion_ut
