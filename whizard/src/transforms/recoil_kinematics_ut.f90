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

module recoil_kinematics_ut
  use unit_tests
  use recoil_kinematics_uti

  implicit none
  private

  public :: recoil_kinematics_test

contains

  subroutine recoil_kinematics_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (recoil_kinematics_1, "recoil_kinematics_1", &
         "iterative solution of non-collinear kinematics", &
         u, results)
    call test (recoil_kinematics_2, "recoil_kinematics_2", &
         "Q distribution", &
         u, results)
    call test (recoil_kinematics_3, "recoil_kinematics_3", &
         "generate recoil event", &
         u, results)
    call test (recoil_kinematics_4, "recoil_kinematics_4", &
         "reference frame", &
         u, results)
    call test (recoil_kinematics_5, "recoil_kinematics_5", &
         "initial reference frame", &
         u, results)
    call test (recoil_kinematics_6, "recoil_kinematics_6", &
         "massless/massive on-shell projection", &
         u, results)
  end subroutine recoil_kinematics_test


end module recoil_kinematics_ut
