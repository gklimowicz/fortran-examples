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

module beam_structures_ut
  use unit_tests
  use beam_structures_uti

  implicit none
  private

  public :: beam_structures_test

contains

  subroutine beam_structures_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (beam_structures_1, "beam_structures_1", &
         "empty beam structure record", &
         u, results)
    call test (beam_structures_2, "beam_structures_2", &
         "beam structure records", &
         u, results)
    call test (beam_structures_3, "beam_structures_3", &
         "beam structure expansion", &
         u, results)
    call test (beam_structures_4, "beam_structures_4", &
         "beam structure contents", &
         u, results)
    call test (beam_structures_5, "beam_structures_5", &
         "polarization", &
         u, results)
    call test (beam_structures_6, "beam_structures_6", &
         "momenta", &
         u, results)
end subroutine beam_structures_test


end module beam_structures_ut
