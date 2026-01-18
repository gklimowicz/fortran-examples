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

module hepmc_interface_ut
  use unit_tests
  use system_dependencies, only: HEPMC2_AVAILABLE
  use system_dependencies, only: HEPMC3_AVAILABLE
  use hepmc_interface_uti

  implicit none
  private

  public :: hepmc_interface_test

contains

  subroutine hepmc_interface_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    if (HEPMC2_AVAILABLE) then
       call test (hepmc_interface_1, "hepmc2_interface_1", &
            "check HepMC2 interface", &
            u, results)
    else if (HEPMC3_AVAILABLE) then
       call test (hepmc_interface_1, "hepmc3_interface_1", &
            "check HepMC3 interface", &
            u, results)
    end if
  end subroutine hepmc_interface_test


end module hepmc_interface_ut
