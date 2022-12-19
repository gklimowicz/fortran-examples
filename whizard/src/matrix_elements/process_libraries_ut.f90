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

module process_libraries_ut
  use unit_tests
  use process_libraries_uti

  implicit none
  private

  public :: process_libraries_test

contains

  subroutine process_libraries_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (process_libraries_1, "process_libraries_1", &
         "empty process list", &
         u, results)
    call test (process_libraries_2, "process_libraries_2", &
         "process definition list", &
         u, results)
    call test (process_libraries_3, "process_libraries_3", &
         "recover process definition list from file", &
         u, results)
    call test (process_libraries_4, "process_libraries_4", &
         "build and load internal process library", &
         u, results)
    call test (process_libraries_5, "process_libraries_5", &
         "build external process library", &
         u, results)
    call test (process_libraries_6, "process_libraries_6", &
         "build and load external process library", &
         u, results)
    call test (process_libraries_7, "process_libraries_7", &
         "process definition list", &
         u, results)
    call test (process_libraries_8, "process_libraries_8", &
         "library status checks", &
         u, results)
  end subroutine process_libraries_test


end module process_libraries_ut
