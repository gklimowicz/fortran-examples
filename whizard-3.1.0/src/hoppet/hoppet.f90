!$Id: hoppet.f90 6133 2014-09-17 14:42:33Z kilian $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     with contributions from
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

! HOPPET wrapper

!----------------------------------------------------------------------
! The module containing the HOPPET pdf
module whizard_hoppet_extras
  use system_dependencies, only: LHAPDF5_AVAILABLE !NODEP!
  use system_dependencies, only: LHAPDF6_AVAILABLE !NODEP!
  use streamlined_interface
  use dglap_holders
  implicit none

  type(pdf_table), save :: whizard_table
end module whizard_hoppet_extras


!----------------------------------------------------------------------
! Get a copy of the LHAPDF structure function, and subtract the LO gluon
! splitting contribution from the b quark pdf
subroutine InitForWhizard (pdf_builtin, pdf, pdf_id)
  use whizard_hoppet_extras
  use pdf_builtin
  use lhapdf

  integer :: iQ, lcl_nf, nloop, order, factscheme, id
  double precision :: mb, Q, dy, ymax, Qmin, Qmax, dlnlnQ
  logical, intent(in) :: pdf_builtin
  type(lhapdf_pdf_t), intent(inout), optional :: pdf
  integer, intent(in), optional :: pdf_id
  external :: evolvePDF

  ! set hoppet start parameters
  dy    = 0.1d0      ! the internal grid spacing (smaller->higher accuarcy)
                     ! 0.1 should provide at least 10^{-3} accuracy 
  nloop = 2          ! the number of loops to initialise (max=3!)
  !--------------------------------------
  ! set extended start parameters
  ymax = 15.0d0      ! highest value of ln1/x user wants to access (default: 12.0d0)
  Qmin = 1.0d0       ! (default)
  Qmax = 28000d0     ! twice LHC c.o.m. (default)
  dlnlnQ = dy/4.0_dp ! min(0.5_dp*dy,0.07_dp) (default)
  order = -6         ! order of numerical interpolation (default)
  factscheme = 1     ! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar (default)

  !!! Default builtin PDF is CTEQ6L
  if (present (pdf_id)) then
     id = pdf_id
  else
     id = 3
  end if

  ! start the hoppet evolution/convolution package 
!  call hoppetStart(dy, nloop)
  call hoppetStartExtended &
       (ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme)

  ! initialise our PDF using the LHAPDF subroutine for PDF-access
  ! (any other subroutine with same interface can be used in its place)
  if (pdf_builtin) then
     call hoppetAssign (pdf_builtin_evolvePDF)
  else if (LHAPDF5_AVAILABLE) then
     call hoppetAssign (evolvePDF)
  else if (LHAPDF6_AVAILABLE) then
     call hoppetAssign (lhapdf6_evolvePDF)
  end if

  ! allocate a PDF table from from the streamlined interface (table)
  call AllocPdfTable (whizard_table, tables(0))
  whizard_table%tab = tables(0)%tab

  lcl_nf=5
  call SetNfDglapHolder (dh, lcl_nf)

  ! use lhapdf's b mass
  if (pdf_builtin) then
     mb = pdf_getmass (5)
  else if (LHAPDF5_AVAILABLE) then
     call GetQmass (5, mb)
  else if (LHAPDF6_AVAILABLE) then
     mb = pdf%get_qmass (5)
  end if

  do iQ = lbound(whizard_table%tab,dim=3), ubound(whizard_table%tab,dim=3)

    Q = whizard_table%Q_vals (iQ)
    if (Q > mb) then

      ! NB 1/(2*lcl_nf) factor is because Pqg includes a 2nf factor
      whizard_table%tab(:,5,iQ) = whizard_table%tab(:,5,iQ) &
           &  - (dh%P_LO%qg * whizard_table%tab(:,0,iQ))  &
           &    * ( log(Q**2/mb**2) * hoppet_alphasQCD(Q)/ (2*pi*2*lcl_nf) )
      whizard_table%tab(:,-5,iQ) = whizard_table%tab(:,-5,iQ) &
           &  - (dh%P_LO%qg * whizard_table%tab(:,0,iQ))  &
           &    * ( log(Q**2/mb**2) * hoppet_alphasQCD(Q)/ (2*pi*2*lcl_nf) )
    end if
  end do
contains
  subroutine lhapdf6_evolvePDF (x, q, ff)
    double precision, intent(in) :: x, q
    double precision, dimension(-6:6), intent(out) :: ff
    call pdf%evolve_pdfm (x, q, ff)
  end subroutine lhapdf6_evolvePDF
  
  subroutine pdf_builtin_evolvePDF (x, q, ff)
    double precision, intent(in) :: x, q
    double precision, dimension(-6:6), intent(out) :: ff
    call pdf_evolve_LHAPDF (id, x, q, ff)
  end subroutine pdf_builtin_evolvePDF

  function hoppet_alphasQCD (q) result (as)
    double precision, intent(in) :: q
    double precision :: as
    double precision :: alphasPDF
    if (pdf_builtin) then
       as = pdf_alphas_LHAPDF (id, q)
    else if (LHAPDF5_AVAILABLE) then
       as = alphasPDF (q)
    else if (LHAPDF6_AVAILABLE) then
       as = pdf%alphas_pdf (q)
    end if
  end function hoppet_alphasQCD

end subroutine InitForWhizard


!----------------------------------------------------------------------
! Return in f(-6:6) the value of the internally stored pdf at the
! given x,Q, with the usual LHApdf meanings for the indices -6:6.
subroutine EvalForWhizard (x,Q,f)
  use whizard_hoppet_extras
  implicit none
  double precision, intent(in)  :: x, Q
  double precision, intent(out) :: f(-6:6)
  
  call EvalPdfTable_xQ (whizard_table,x,Q,f)
end subroutine EvalForWhizard
