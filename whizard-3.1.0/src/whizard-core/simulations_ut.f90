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

module simulations_ut
  use unit_tests
  use simulations_uti

  implicit none
  private

  public :: simulations_test

contains

  subroutine simulations_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (simulations_1, "simulations_1", &
         "initialization", &
         u, results)
    call test (simulations_2, "simulations_2", &
         "weighted events", &
         u, results)
    call test (simulations_3, "simulations_3", &
         "unweighted events", &
         u, results)
    call test (simulations_4, "simulations_4", &
         "process with structure functions", &
         u, results)
    call test (simulations_5, "simulations_5", &
         "raw event I/O", &
         u, results)
    call test (simulations_6, "simulations_6", &
         "raw event I/O with structure functions", &
         u, results)
    call test (simulations_7, "simulations_7", &
         "automatic raw event I/O", &
         u, results)
    call test (simulations_8, "simulations_8", &
         "rescan raw event file", &
         u, results)
    call test (simulations_9, "simulations_9", &
         "rescan mismatch", &
         u, results)
    call test (simulations_10, "simulations_10", &
         "alternative weight", &
         u, results)
    call test (simulations_11, "simulations_11", &
         "decay", &
         u, results)
    call test (simulations_12, "simulations_12", &
         "split event files", &
         u, results)
    call test (simulations_13, "simulations_13", &
         "callback", &
         u, results)
    call test (simulations_14, "simulations_14", &
         "resonant subprocesses evaluation", &
         u, results)
    call test (simulations_15, "simulations_15", &
         "resonant subprocesses in simulation", &
         u, results)
  end subroutine simulations_test


end module simulations_ut
