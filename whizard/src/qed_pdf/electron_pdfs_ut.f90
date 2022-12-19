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

module electron_pdfs_ut
  use unit_tests
  use electron_pdfs_uti

  implicit none
  private

  public :: electron_pdfs_test

contains

  subroutine electron_pdfs_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test (electron_pdfs_1, "electron_pdfs_1", &
         "Electron PDFs: auxiliary functions", &
         u, results)
    call test (electron_pdfs_2, "electron_pdfs_2", &
         "Electron PDFs: auxiliary functions (2)", &
         u, results)
    call test (electron_pdfs_3, "electron_pdfs_3", &
         "Electron PDFs: auxiliary functions (3)", &
         u, results)
    call test (electron_pdfs_4, "electron_pdfs_4", &
         "Electron PDFs: auxiliary functions (4)", &
         u, results)
    call test (electron_pdfs_5, "electron_pdfs_5", &
         "Electron PDFs: auxiliary functions (5)", &
         u, results)
    call test (electron_pdfs_6, "electron_pdfs_6", &
         "Electron PDFs: full electron PDFs", &
         u, results)
  end subroutine electron_pdfs_test


end module electron_pdfs_ut
