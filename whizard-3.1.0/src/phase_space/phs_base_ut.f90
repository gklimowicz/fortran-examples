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

module phs_base_ut
  use unit_tests
  use phs_base_uti

  implicit none
  private

  public :: phs_base_test

  public :: init_test_process_data
  public :: init_test_decay_data
  public :: phs_test_config_t
  public :: phs_test_t

contains

  subroutine phs_base_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (phs_base_1, "phs_base_1", &
         "phase-space configuration", &
         u, results)
    call test (phs_base_2, "phs_base_2", &
         "phase-space evaluation", &
         u, results)
    call test (phs_base_3, "phs_base_3", &
         "channel equivalences", &
         u, results)
    call test (phs_base_4, "phs_base_4", &
         "MD5 sum", &
         u, results)
    call test (phs_base_5, "phs_base_5", &
         "channel collection", &
         u, results)
  end subroutine phs_base_test


end module phs_base_ut
