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

module vamp2_ut
  use unit_tests
  use vamp2_uti

  implicit none
  private

  public :: vamp2_test

contains
  subroutine vamp2_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
      call test (vamp2_1, "vamp2_1", "VAMP2 initialisation and&
           & grid preparation", u, results)
      call test (vamp2_2, "vamp2_2", "VAMP2 intgeration of two-dimensional &
           & function with two channels", u, results)
      call test (vamp2_3, "vamp2_3", "VAMP2 intgeration of two-dimensional &
           & function with two channels", u, results)
      call test (vamp2_4, "vamp2_4", "VAMP2 intgeration of two-dimensional &
           & function with two channels with chains", u, results)
      call test (vamp2_5, "vamp2_5", "VAMP2 intgeration of two-dimensional &
           & function with two channels with equivalences", u, results)
  end subroutine vamp2_test

end module vamp2_ut
