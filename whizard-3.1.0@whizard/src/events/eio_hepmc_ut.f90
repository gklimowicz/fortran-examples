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

module eio_hepmc_ut
  use unit_tests
  use system_dependencies, only: HEPMC2_AVAILABLE
  use system_dependencies, only: HEPMC3_AVAILABLE
  use eio_hepmc_uti

  implicit none
  private

  public :: eio_hepmc_test

contains

  subroutine eio_hepmc_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    if (HEPMC2_AVAILABLE) then
       call test (eio_hepmc_1, "eio_hepmc2_1", &
            "write event contents", &
            u, results)
    else if (HEPMC3_AVAILABLE) then
       call test (eio_hepmc_1, "eio_hepmc3_1", &
            "write event contents", &
            u, results)
    end if
    if (HEPMC2_AVAILABLE) then
       call test (eio_hepmc_2, "eio_hepmc2_2", &
            "read event contents", &
            u, results)
    else if (HEPMC3_AVAILABLE) then
       call test (eio_hepmc_2, "eio_hepmc3_2", &
            "read event contents", &
            u, results)
    end if
    if (HEPMC2_AVAILABLE) then
       call test (eio_hepmc_3, "eio_hepmc2_3", &
            "write event contents", &
            u, results)
    else if (HEPMC3_AVAILABLE) then
       call test (eio_hepmc_3, "eio_hepmc3_3", &
            "event contents weighted, '1' normalization", &
            u, results)
    end if
  end subroutine eio_hepmc_test


end module eio_hepmc_ut
