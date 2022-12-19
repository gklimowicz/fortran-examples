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

module isr_epa_handler_ut
  use unit_tests
  use isr_epa_handler_uti

  implicit none
  private

  public :: isr_handler_test
  public :: epa_handler_test

contains

  subroutine isr_handler_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (isr_handler_1, "isr_handler_1", &
         "collinear case, no modification", &
         u, results)
    call test (isr_handler_2, "isr_handler_2", &
         "two-photon recoil", &
         u, results)
    call test (isr_handler_3, "isr_handler_3", &
         "two-photon recoil with boost", &
         u, results)
  end subroutine isr_handler_test

  subroutine epa_handler_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (epa_handler_1, "epa_handler_1", &
         "collinear case, no modification", &
         u, results)
    call test (epa_handler_2, "epa_handler_2", &
         "two-beam recoil", &
         u, results)
    call test (epa_handler_3, "epa_handler_3", &
         "two-beam recoil with boost", &
         u, results)
  end subroutine epa_handler_test


end module isr_epa_handler_ut
