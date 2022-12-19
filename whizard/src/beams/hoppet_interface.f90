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

module hoppet_interface
  use lhapdf !NODEP!

  implicit none
  private

  public :: hoppet_init, hoppet_eval

contains

  subroutine hoppet_init (pdf_builtin, pdf, pdf_id)
    logical, intent(in) :: pdf_builtin
    type(lhapdf_pdf_t), intent(inout), optional :: pdf
    integer, intent(in), optional :: pdf_id
    external InitForWhizard
    call InitForWhizard (pdf_builtin, pdf, pdf_id)
  end subroutine hoppet_init

  subroutine hoppet_eval (x, q, f)
    double precision, intent(in)  :: x, q
    double precision, intent(out) :: f(-6:6)
    external EvalForWhizard
    call EvalForWhizard (x, q, f)
  end subroutine hoppet_eval

end module hoppet_interface
