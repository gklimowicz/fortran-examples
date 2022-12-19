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

module radiation_generator_ut
  use unit_tests
  use radiation_generator_uti

  implicit none
  private

  public :: radiation_generator_test

contains

  subroutine radiation_generator_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test(radiation_generator_1, "radiation_generator_1", &
         "Test the generator of N+1-particle flavor structures in QCD", &
         u, results)
    call test(radiation_generator_2, "radiation_generator_2", &
         "Test multiple splittings in QCD", &
         u, results)
    call test(radiation_generator_3, "radiation_generator_3", &
         "Test the generator of N+1-particle flavor structures in QED", &
         u, results)
    call test(radiation_generator_4, "radiation_generator_4", &
         "Test multiple splittings in QED", &
         u, results)
  end subroutine radiation_generator_test


end module radiation_generator_ut
