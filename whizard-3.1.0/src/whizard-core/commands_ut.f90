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

module commands_ut
  use unit_tests
  use system_dependencies, only: MPOST_AVAILABLE
  use commands_uti

  implicit none
  private

  public :: commands_test

contains

  subroutine commands_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (commands_1, "commands_1", &
         "empty command list", &
         u, results)
    call test (commands_2, "commands_2", &
         "model", &
         u, results)
    call test (commands_3, "commands_3", &
         "process declaration", &
         u, results)
    call test (commands_4, "commands_4", &
         "compilation", &
         u, results)
    call test (commands_5, "commands_5", &
         "integration", &
         u, results)
    call test (commands_6, "commands_6", &
         "variables", &
         u, results)
    call test (commands_7, "commands_7", &
         "process library", &
         u, results)
    call test (commands_8, "commands_8", &
         "event generation", &
         u, results)
    call test (commands_9, "commands_9", &
         "cuts", &
         u, results)
    call test (commands_10, "commands_10", &
         "beams", &
         u, results)
    call test (commands_11, "commands_11", &
         "structure functions", &
         u, results)
    call test (commands_12, "commands_12", &
         "event rescanning", &
         u, results)
    call test (commands_13, "commands_13", &
         "event output formats", &
         u, results)
    call test (commands_14, "commands_14", &
         "empty libraries", &
         u, results)
    call test (commands_15, "commands_15", &
         "compilation", &
         u, results)
    call test (commands_16, "commands_16", &
         "observables", &
         u, results)
    call test (commands_17, "commands_17", &
         "histograms", &
         u, results)
    call test (commands_18, "commands_18", &
         "plots", &
         u, results)
    call test (commands_19, "commands_19", &
         "graphs", &
         u, results)
    call test (commands_20, "commands_20", &
         "record data", &
         u, results)
    call test (commands_21, "commands_21", &
         "analysis expression", &
         u, results)
    call test (commands_22, "commands_22", &
         "write analysis", &
         u, results)
    if (MPOST_AVAILABLE) then
       call test (commands_23, "commands_23", &
            "compile analysis", &
            u, results)
    end if
    call test (commands_24, "commands_24", &
         "drawing options", &
         u, results)
    call test (commands_25, "commands_25", &
         "local process environment", &
         u, results)
    call test (commands_26, "commands_26", &
         "alternative setups", &
         u, results)
    call test (commands_27, "commands_27", &
         "unstable and polarized particles", &
         u, results)
    call test (commands_28, "commands_28", &
         "quit", &
         u, results)
    call test (commands_29, "commands_29", &
         "SLHA interface", &
         u, results)
    call test (commands_30, "commands_30", &
         "scales", &
         u, results)
    call test (commands_31, "commands_31", &
         "event weights/reweighting", &
         u, results)
    call test (commands_32, "commands_32", &
         "event selection", &
         u, results)
    call test (commands_33, "commands_33", &
         "execute shell command", &
         u, results)
    call test (commands_34, "commands_34", &
         "analysis via callback", &
         u, results)
  end subroutine commands_test


end module commands_ut
