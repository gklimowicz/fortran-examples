!$Id: lhapdf.f90 6133 2014-09-17 14:42:33Z kilian $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     with contributions from
!     Bijan Chokoufe <bijan.chokoufe@desy.de>
!     Fabian Bach <fabian.bach@t-online.de>
!     Christian Speckner <cnspeckn@googlemail.com>
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

module lhapdf

  use, intrinsic :: iso_c_binding
  use kinds

  implicit none
  private

  ! Public types
  public :: lhapdf_pdf_t

  public :: lhapdf_copy_pointer

  type :: lhapdf_pdf_t
     private
     type(c_ptr) :: cptr
   contains
     procedure :: init => lhapdf_pdf_init
     procedure :: is_associated => lhapdf_is_associated
     procedure :: getxmin => lhapdf_getxmin
     procedure :: getxmax => lhapdf_getxmax
     procedure :: getq2min => lhapdf_getq2min
     procedure :: getq2max => lhapdf_getq2max
     procedure :: has_photon => lhapdf_hasphoton
     procedure :: evolve_pdfm => lhapdf_evolve_pdfm
     procedure :: evolve_pdfphotonm => lhapdf_evolve_pdfphotonm
     procedure :: evolve_pdfpm => lhapdf_evolve_pdfpm
     procedure :: get_qmass => lhapdf_get_qmass
     procedure :: get_order => lhapdf_get_order
     procedure :: num_pdfm => lhapdf_num_pdfm
     procedure :: alphas_pdf => lhapdf_alphas_pdf
     procedure :: final => lhapdf_final
  end type lhapdf_pdf_t

  ! Interface for generic operators

  interface
     function lhapdf_init_pdf (setname, imem) bind (C) result (pdf)
       import
       integer(c_int), intent(in), value :: imem
       character(len=1, kind=c_char), dimension(*), intent(in) :: setname
       type(c_ptr) :: pdf
     end function lhapdf_init_pdf
  end interface

  interface
     subroutine lhapdf_pdf_delete (pdf) bind (C)
       import
       type(c_ptr), value :: pdf
     end subroutine lhapdf_pdf_delete
  end interface

  interface
     function lhapdf_pdf_getxmin (pdf) bind (C) result (xmin)
       import
       type(c_ptr), intent(in), value :: pdf
       real(c_double) :: xmin
     end function lhapdf_pdf_getxmin
  end interface

  interface
     function lhapdf_pdf_getxmax (pdf) bind (C) result (xmax)
       import
       type(c_ptr), intent(in), value :: pdf
       real(c_double) :: xmax
     end function lhapdf_pdf_getxmax
  end interface

  interface
     function lhapdf_pdf_getq2min (pdf) bind (C) result (q2min)
       import
       type(c_ptr), intent(in), value :: pdf
       real(c_double) :: q2min
     end function lhapdf_pdf_getq2min
  end interface

  interface
     function lhapdf_pdf_getq2max (pdf) bind (C) result (q2max)
       import
       type(c_ptr), intent(in), value :: pdf
       real(c_double) :: q2max
     end function lhapdf_pdf_getq2max
  end interface

  interface
     function lhapdf_has_photon (pdf) bind (C) result (flag)
       import
       type(c_ptr), intent(in), value :: pdf
       logical(c_bool) :: flag
     end function lhapdf_has_photon
  end interface

  interface
     subroutine lhapdf_evolvepdfm (pdf, x, q, ff) bind (C)
       import
       type(c_ptr), intent(in), value :: pdf
       real(c_double), intent(in), value :: x, q
       real(c_double), dimension(-6:6), intent(out) :: ff
     end subroutine lhapdf_evolvepdfm
  end interface

  interface
     subroutine lhapdf_evolvepdfphotonm (pdf, x, q, ff, fphot) bind (C)
       import
       type(c_ptr), intent(in), value :: pdf
       real(c_double), intent(in), value :: x, q
       real(c_double), dimension(-6:6), intent(out) :: ff
       real(c_double), intent(out) :: fphot
     end subroutine lhapdf_evolvepdfphotonm
  end interface

  interface
     subroutine lhapdf_evolvepdfpm (pdf, x, q, s, scheme, ff) bind (C)
       import
       type(c_ptr), intent(in), value :: pdf
       integer(c_int), intent(in), value :: scheme
       real(c_double), intent(in), value :: x, q, s
       real(c_double), dimension(-6:6), intent(out) :: ff
     end subroutine lhapdf_evolvepdfpm
  end interface

  interface
     function lhapdf_getqmass (pdf, nf) bind (C) result (mass)
       import
       type(c_ptr), intent(in), value :: pdf
       integer(c_int), intent(in), value :: nf
       real(c_double) :: mass
     end function lhapdf_getqmass
  end interface

  interface
     function lhapdf_getorder (pdf) bind (C) result (order)
       import
       type(c_ptr), intent(in), value :: pdf
       integer(c_int) :: order
     end function lhapdf_getorder
  end interface

  interface
     function lhapdf_numpdfm (pdf) bind (C) result (numpdf)
       import
       type(c_ptr), intent(in), value :: pdf
       integer(c_int) :: numpdf
     end function lhapdf_numpdfm
  end interface

  interface
     function lhapdf_alphaspdf (pdf, q) bind (C) result (as)
       import
       type(c_ptr), intent(in), value :: pdf
       real(c_double), intent(in), value :: q
       real(c_double) :: as
     end function lhapdf_alphaspdf
  end interface

contains

  subroutine lhapdf_pdf_init (pdf, setname, imem)
    class(lhapdf_pdf_t), intent(out) :: pdf
    character(*), intent(in) :: setname
    integer, intent(in) :: imem
    character(len=1, kind=c_char), dimension(len(setname)+1) :: pdf_setname
    integer(c_int) :: pdf_imem
    integer :: i, strlen
    strlen = len(setname)
    forall (i=1:strlen)
       pdf_setname(i) = setname(i:i)
    end forall
    pdf_setname(strlen+1) = c_null_char
    pdf_imem = imem
    pdf%cptr = lhapdf_init_pdf (pdf_setname, pdf_imem)
  end subroutine lhapdf_pdf_init

  function lhapdf_is_associated (pdf) result (flag)
    class(lhapdf_pdf_t), intent(in) :: pdf
    logical :: flag
    flag = c_associated (pdf%cptr)
  end function lhapdf_is_associated

  function lhapdf_getxmin (pdf) result (xmin)
    class(lhapdf_pdf_t), intent(in) :: pdf
    real(double) :: xmin
    xmin = lhapdf_pdf_getxmin (pdf%cptr)
  end function lhapdf_getxmin

  function lhapdf_getxmax (pdf) result (xmax)
    class(lhapdf_pdf_t), intent(in) :: pdf
    real(double) :: xmax
    xmax = lhapdf_pdf_getxmax (pdf%cptr)
  end function lhapdf_getxmax

  function lhapdf_getq2min (pdf) result (q2min)
    class(lhapdf_pdf_t), intent(in) :: pdf
    real(double) :: q2min
    q2min = lhapdf_pdf_getq2min (pdf%cptr)
  end function lhapdf_getq2min

  function lhapdf_getq2max (pdf) result (q2max)
    class(lhapdf_pdf_t), intent(in) :: pdf
    real(double) :: q2max
    q2max = lhapdf_pdf_getq2max (pdf%cptr)
  end function lhapdf_getq2max

  function lhapdf_hasphoton (pdf) result (flag)
    class(lhapdf_pdf_t), intent(in) :: pdf
    logical :: flag
    flag = lhapdf_has_photon (pdf%cptr)
  end function lhapdf_hasphoton

  subroutine lhapdf_evolve_pdfm (pdf, x, q, ff)
    class(lhapdf_pdf_t), intent(inout) :: pdf
    real(double), intent(in) :: x, q
    real(double), dimension(-6:6), intent(out) :: ff
    real(c_double) :: c_x, c_q
    c_x = x
    c_q = q
    call lhapdf_evolvepdfm (pdf%cptr, c_x, c_q, ff)
  end subroutine lhapdf_evolve_pdfm

  subroutine lhapdf_evolve_pdfphotonm (pdf, x, q, ff, fphot)
    class(lhapdf_pdf_t), intent(inout) :: pdf
    real(double), intent(in) :: x, q
    real(double), dimension(-6:6), intent(out) :: ff
    real(double), intent(out) :: fphot
    real(c_double) :: c_x, c_q
    c_x = x
    c_q = q
    call lhapdf_evolvepdfphotonm &
         (pdf%cptr, c_x, c_q, ff, fphot)
  end subroutine lhapdf_evolve_pdfphotonm

  subroutine lhapdf_evolve_pdfpm (pdf, x, q, s, scheme, ff)
    class(lhapdf_pdf_t), intent(inout) :: pdf
    real(double), intent(in) :: x, q, s
    integer, intent(in) :: scheme
    real(double), dimension(-6:6), intent(out) :: ff
    real(c_double) :: c_x, c_q, c_s
    integer(c_int) :: c_scheme
    c_x = x
    c_q = q
    c_s = s
    c_scheme = scheme
    call lhapdf_evolvepdfpm (pdf%cptr, &
         c_x, c_q, c_s, c_scheme, ff)
  end subroutine lhapdf_evolve_pdfpm

  function lhapdf_get_qmass (pdf, nf) result (mass)
    class(lhapdf_pdf_t), intent(in) :: pdf
    integer, intent(in) :: nf
    real(double) :: mass
    integer(c_int) :: c_nf
    c_nf = nf
    mass = lhapdf_getqmass (pdf%cptr, c_nf)
  end function lhapdf_get_qmass

  function lhapdf_get_order (pdf) result (order)
    class(lhapdf_pdf_t), intent(in) :: pdf
    integer(c_int) :: order
    order = lhapdf_getorder (pdf%cptr)
  end function lhapdf_get_order

  function lhapdf_num_pdfm (pdf) result (numpdf)
    class(lhapdf_pdf_t), intent(inout) :: pdf
    integer :: numpdf
    numpdf = lhapdf_numpdfm (pdf%cptr)
  end function lhapdf_num_pdfm

  function lhapdf_alphas_pdf (pdf, q) result (as)
    class(lhapdf_pdf_t), intent(in), target :: pdf
    real(double), intent(in) :: q
    real(double) :: as
    real(c_double) :: c_q = 0
    c_q = q
    as = lhapdf_alphaspdf (pdf%cptr, c_q)
  end function lhapdf_alphas_pdf

  subroutine lhapdf_final (pdf)
    class(lhapdf_pdf_t), intent(inout) :: pdf
    call lhapdf_pdf_delete (pdf%cptr)
  end subroutine lhapdf_final

  subroutine lhapdf_copy_pointer (pdf_in, pdf_out)
    type(lhapdf_pdf_t), intent(in), target :: pdf_in
    type(lhapdf_pdf_t), intent(out), target :: pdf_out
    pdf_out%cptr = transfer (pdf_in%cptr, pdf_out%cptr)
  end subroutine lhapdf_copy_pointer

end module lhapdf
