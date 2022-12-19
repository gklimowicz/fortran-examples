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

module pdf

  use kinds, only: default, double
  use beam_structures
  use lhapdf !NODEP!
  use pdf_builtin !NODEP!

  implicit none
  private

  public :: pdf_data_t

  integer, parameter, public :: STRF_NONE = 0
  integer, parameter, public :: STRF_LHAPDF6 = 1
  integer, parameter, public :: STRF_LHAPDF5 = 2
  integer, parameter, public :: STRF_PDF_BUILTIN = 3


  type :: pdf_data_t
    type(lhapdf_pdf_t) :: pdf
    real(default) :: xmin, xmax, qmin, qmax
    integer :: type = STRF_NONE
    integer :: set = 0
   contains
     procedure :: init => pdf_data_init
     procedure :: write => pdf_data_write
     procedure :: setup => pdf_data_setup
     procedure :: evolve => pdf_data_evolve
  end type pdf_data_t


  interface
    module subroutine pdf_data_init (pdf_data, pdf_data_in)
      class(pdf_data_t), intent(out) :: pdf_data
      type(pdf_data_t), target, intent(in) :: pdf_data_in
    end subroutine pdf_data_init
    module subroutine pdf_data_write (pdf_data, unit)
      class(pdf_data_t), intent(in) :: pdf_data
      integer, intent(in), optional :: unit
    end subroutine pdf_data_write
    module subroutine pdf_data_setup &
         (pdf_data, caller, beam_structure, lhapdf_member, set)
      class(pdf_data_t), intent(inout) :: pdf_data
      character(len=*), intent(in) :: caller
      type(beam_structure_t), intent(in) :: beam_structure
      integer, intent(in) :: lhapdf_member, set
    end subroutine pdf_data_setup
    module subroutine pdf_data_evolve (pdf_data, x, q_in, f)
      class(pdf_data_t), intent(inout) :: pdf_data
      real(double), intent(in) :: x, q_in
      real(double), dimension(-6:6), intent(out) :: f
    end subroutine pdf_data_evolve
  end interface

end module pdf
