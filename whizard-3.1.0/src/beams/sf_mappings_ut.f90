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

module sf_mappings_ut
  use unit_tests
  use sf_mappings_uti

  implicit none
  private

  public :: sf_mappings_test

contains

  subroutine sf_mappings_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (sf_mappings_1, "sf_mappings_1", &
         "standard pair mapping", &
         u, results)
    call test (sf_mappings_2, "sf_mappings_2", &
         "structure-function mapping channels", &
         u, results)
    call test (sf_mappings_3, "sf_mappings_3", &
         "resonant pair mapping", &
         u, results)
    call test (sf_mappings_4, "sf_mappings_4", &
         "on-shell pair mapping", &
         u, results)
    call test (sf_mappings_5, "sf_mappings_5", &
         "endpoint pair mapping", &
         u, results)
    call test (sf_mappings_6, "sf_mappings_6", &
         "endpoint resonant mapping", &
         u, results)
    call test (sf_mappings_7, "sf_mappings_7", &
         "endpoint on-shell mapping", &
         u, results)
    call test (sf_mappings_8, "sf_mappings_8", &
         "power pair mapping", &
         u, results)
    call test (sf_mappings_9, "sf_mappings_9", &
         "power resonance mapping", &
         u, results)
    call test (sf_mappings_10, "sf_mappings_10", &
         "power on-shell mapping", &
         u, results)
    call test (sf_mappings_11, "sf_mappings_11", &
         "endpoint/power combined mapping", &
         u, results)
    call test (sf_mappings_12, "sf_mappings_12", &
         "endpoint/power resonant combined mapping", &
         u, results)
    call test (sf_mappings_13, "sf_mappings_13", &
         "endpoint/power on-shell combined mapping", &
         u, results)
    call test (sf_mappings_14, "sf_mappings_14", &
         "rescaled on-shell mapping", &
         u, results)
    call test (sf_mappings_15, "sf_mappings_15", &
         "resonant single mapping", &
         u, results)
    call test (sf_mappings_16, "sf_mappings_16", &
         "on-shell single mapping", &
         u, results)
  end subroutine sf_mappings_test


end module sf_mappings_ut