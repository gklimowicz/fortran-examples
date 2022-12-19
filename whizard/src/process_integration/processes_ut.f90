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

module processes_ut
  use unit_tests
  use processes_uti

  implicit none
  private

  public :: processes_test

  public :: prepare_test_process
  public :: cleanup_test_process

contains

  subroutine processes_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (processes_1, "processes_1", &
         "write an empty process object", &
         u, results)
    call test (processes_2, "processes_2", &
         "initialize a simple process object", &
         u, results)
    call test (processes_3, "processes_3", &
         "retrieve a trivial matrix element", &
         u, results)
    call test (processes_4, "processes_4", &
         "create and fill a process instance (partonic event)", &
         u, results)
    call test (processes_7, "processes_7", &
         "process configuration with structure functions", &
         u, results)
    call test (processes_8, "processes_8", &
         "process evaluation with structure functions", &
         u, results)
    call test (processes_9, "processes_9", &
         "multichannel kinematics and structure functions", &
         u, results)
    call test (processes_10, "processes_10", &
         "event generation", &
         u, results)
    call test (processes_11, "processes_11", &
         "integration", &
         u, results)
    call test (processes_12, "processes_12", &
         "event post-processing", &
         u, results)
    call test (processes_13, "processes_13", &
         "colored interaction", &
         u, results)
    call test (processes_14, "processes_14", &
         "process configuration and MD5 sum", &
         u, results)
    call test (processes_15, "processes_15", &
         "decay process", &
         u, results)
    call test (processes_16, "processes_16", &
         "decay integration", &
         u, results)
    call test (processes_17, "processes_17", &
         "decay of moving particle", &
         u, results)
    call test (processes_18, "processes_18", &
         "extract resonance history set", &
         u, results)
    call test (processes_19, "processes_19", &
         "add trivial hooks to a process instance ", &
         u, results)
  end subroutine processes_test


end module processes_ut
