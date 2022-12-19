!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

! Wrap a common interface around the different PDF sets.

module pdf_builtin
  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string

  implicit none
  save

  private

  public :: pdf_init
  public :: pdf_get_name
  public :: pdf_evolve
  public :: pdf_evolve_LHAPDF
  public :: pdf_provides_photon
  public :: pdf_get_id
  public :: pdf_alphas
  public :: pdf_alphas_LHAPDF
  public :: pdf_getmass

  interface
    module subroutine pdf_init (pdftype, prefix, verbose)
      integer, intent(in) :: pdftype
      type(string_t), intent(in), optional :: prefix
      logical, intent(in), optional :: verbose
    end subroutine pdf_init
    module subroutine pdf_evolve (pdftype, x, q, f, fphoton)
      integer, intent(in) :: pdftype
      real(kind=default), intent(in) :: x, q
      real(kind=default), intent(out), optional :: f(-6:6), fphoton
    end subroutine pdf_evolve
    module subroutine pdf_evolve_LHAPDF (set, x, q, ff)
      integer, intent(in) :: set
      double precision, intent(in) :: x, q
      double precision, dimension(-6:6), intent(out) :: ff
    end subroutine pdf_evolve_LHAPDF
    module function pdf_get_name (pdftype) result (name)
      integer, intent(in) :: pdftype
      type(string_t) :: name
    end function pdf_get_name
    module function pdf_get_id (name) result (id)
      type(string_t), intent(in) :: name
      integer :: id
    end function pdf_get_id
    module function pdf_provides_photon (pdftype) result (flag)
      integer, intent(in) :: pdftype
      logical :: flag
    end function pdf_provides_photon
    module function pdf_getmass (nf) result (mass)
      integer, intent(in) :: nf
      real(kind=double) :: mass
    end function pdf_getmass
    module function pdf_alphas (pdftype, q) result (alphas)
      real(kind=default), intent(in) :: q
      real(kind=default) :: alphas
      integer, intent(in) :: pdftype
    end function pdf_alphas
    module function pdf_alphas_LHAPDF (pdftype, q) result (alphas)
      integer, intent(in) :: pdftype
      real(kind=double), intent(in) :: q
      real(kind=double) :: alphas
    end function pdf_alphas_LHAPDF
  end interface

end module pdf_builtin
