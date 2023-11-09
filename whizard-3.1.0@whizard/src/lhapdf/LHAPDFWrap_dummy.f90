function lhapdf_init_pdf (setname, imem) bind (C) result (pdf)
  use iso_c_binding
  character(len=1, kind=c_char), dimension(*), intent(in) :: setname
  integer(c_int), intent(in), value :: imem
  type(c_ptr) :: pdf
end function lhapdf_init_pdf

subroutine lhapdf_pdf_delete (pdf) bind (C)
  use iso_c_binding
  type(c_ptr), value :: pdf
end subroutine lhapdf_pdf_delete

function lhapdf_pdf_getxmin (pdf) bind (C) result (xmin)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  real(c_double) :: xmin
end function lhapdf_pdf_getxmin

function lhapdf_pdf_getxmax (pdf) bind (C) result (xmax)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  real(c_double) :: xmax
end function lhapdf_pdf_getxmax

function lhapdf_pdf_getq2min (pdf) bind (C) result (q2min)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  real(c_double) :: q2min
end function lhapdf_pdf_getq2min

function lhapdf_pdf_getq2max (pdf) bind (C) result (q2max)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  real(c_double) :: q2max
end function lhapdf_pdf_getq2max

function lhapdf_has_photon (pdf) bind (C) result (flag)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  logical(c_bool) :: flag
end function lhapdf_has_photon

subroutine lhapdf_evolvepdfm (pdf, x, q, ff) bind (C)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  real(c_double), intent(in), value :: x, q
  real(c_double), dimension(-6:6), intent(out) :: ff
end subroutine lhapdf_evolvepdfm

subroutine lhapdf_evolvepdfphotonm (pdf, x, q, ff, fphot) bind (C)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  real(c_double), intent(in), value :: x, q
  real(c_double), dimension(-6:6), intent(out) :: ff
  real(c_double), intent(out) :: fphot
end subroutine lhapdf_evolvepdfphotonm

subroutine lhapdf_evolvepdfpm (pdf, x, q, s, scheme, ff) bind (C)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  integer(c_int), intent(in), value :: scheme
  real(c_double), intent(in), value :: x, q, s
  real(c_double), dimension(-6:6), intent(out) :: ff
end subroutine lhapdf_evolvepdfpm

function lhapdf_getqmass (pdf, nf) bind (C) result (mass)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  integer(c_int), intent(in), value :: nf
  real(c_double) :: mass
end function lhapdf_getqmass

function lhapdf_numpdfm (pdf) bind (C) result (numpdf)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  integer(c_int) :: numpdf
end function lhapdf_numpdfm

function lhapdf_alphaspdf (pdf) bind (C) result (as)
  use iso_c_binding
  type(c_ptr), intent(in), value :: pdf
  real(c_double) :: as
end function lhapdf_alphaspdf

function lhapdf_getorder (pdf) bind (C) result (order)
  use iso_c_binding
  type(c_ptr), intent(in) :: pdf
  integer(c_int) :: order
end function lhapdf_getorder
