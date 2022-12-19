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

module resonances_ut
  use unit_tests
  use resonances_uti

  implicit none
  private

  public :: resonances_test

contains

  subroutine resonances_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (resonances_1, "resonances_1", &
         "check resonance history setup", &
         u, results)
    call test (resonances_2, "resonances_2", &
         "check O'Mega restriction strings", &
         u, results)
    call test (resonances_3, "resonances_3", &
         "check resonance history set", &
         u, results)
    call test (resonances_4, "resonances_4", &
         "resonance history: distance evaluation", &
         u, results)
    call test (resonances_5, "resonances_5", &
         "resonance history: on-shell test", &
         u, results)
    call test (resonances_6, "resonances_6", &
         "check resonance history setup", &
         u, results)
    call test (resonances_7, "resonances_7", &
         "display tree format of history set elements", &
         u, results)
  end subroutine resonances_test


end module resonances_ut
