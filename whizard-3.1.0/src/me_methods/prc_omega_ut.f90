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

module prc_omega_ut
  use unit_tests
  use prc_omega_uti

  implicit none
  private

  public :: prc_omega_test
  public :: prc_omega_diags_test

contains

  subroutine prc_omega_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (prc_omega_1, "prc_omega_1", &
         "build and load simple OMega process", &
         u, results)
    call test (prc_omega_2, "prc_omega_2", &
         "OMega option passing", &
         u, results)
    call test (prc_omega_3, "prc_omega_3", &
         "helicity selection", &
         u, results)
    call test (prc_omega_4, "prc_omega_4", &
         "update QCD alpha", &
         u, results)
    call test (prc_omega_5, "prc_omega_5", &
         "running QCD alpha", &
         u, results)
    call test (prc_omega_6, "prc_omega_6", &
         "OMega UFO support", &
         u, results)
end subroutine prc_omega_test

  subroutine prc_omega_diags_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (prc_omega_diags_1, "prc_omega_diags_1", &
         "generate Feynman diagrams", &
         u, results)
end subroutine prc_omega_diags_test


end module prc_omega_ut
