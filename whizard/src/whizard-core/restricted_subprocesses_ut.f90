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

module restricted_subprocesses_ut
  use unit_tests
  use restricted_subprocesses_uti

  implicit none
  private

  public :: restricted_subprocesses_test

contains

  subroutine restricted_subprocesses_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (restricted_subprocesses_1, "restricted_subprocesses_1", &
         "single subprocess", &
         u, results)
    call test (restricted_subprocesses_2, "restricted_subprocesses_2", &
         "subprocess library", &
         u, results)
    call test (restricted_subprocesses_3, "restricted_subprocesses_3", &
         "resonance kinematics and probability", &
         u, results)
    call test (restricted_subprocesses_4, "restricted_subprocesses_4", &
         "event transform", &
         u, results)
    call test (restricted_subprocesses_5, "restricted_subprocesses_5", &
         "event transform with gaussian turnoff", &
         u, results)
    call test (restricted_subprocesses_6, "restricted_subprocesses_6", &
         "event transform with background switched off", &
         u, results)
  end subroutine restricted_subprocesses_test


end module restricted_subprocesses_ut
