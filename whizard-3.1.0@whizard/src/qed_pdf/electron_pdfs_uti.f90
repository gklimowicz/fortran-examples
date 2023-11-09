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

module electron_pdfs_uti

  use kinds, only: default
  use numeric_utils
  use format_defs, only: FMT_15
  use constants
  use physics_defs, only: ME_REF, ALPHA_QED_ME_REF
  use electron_pdfs

  implicit none
  private

  public :: electron_pdfs_1
  public :: electron_pdfs_2
  public :: electron_pdfs_3
  public :: electron_pdfs_4
  public :: electron_pdfs_5
  public :: electron_pdfs_6

contains

  subroutine electron_pdfs_1 (u)
    integer, intent(in) :: u
    type(qed_pdf_t) :: pdf
    real(default) :: Q, alpha
    real(default) :: x1, x2, x3, x4
    integer :: n_lep
    Q = 10._default
    x1 = 0.1_default
    x2 = 0.5_default
    x3 = 0.9_default
    x4 = 0.999_default
    alpha = ALPHA_QED_ME_REF
    n_lep = 1

    write (u, "(A)")  "* Test output: electron_pdfs_1"
    write (u, "(A)")  "*   Purpose: check analytic properties"
    write (u, "(A)")

    write (u, "(A)")  "* Auxiliary functions I:"
    write (u, "(A)")

    write (u, "(A)")  "* Q = 10 GeV, elec_asym, LL+NLL, alpha fixed:"
    write (u, "(A)")

    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, 1)
    write (u, "(1x,A,F9.6)")  " elec_asym (LL,x=0.1)   = ", &
         elec_asym (pdf, x1, Q, alpha, .false.)
    write (u, "(1x,A,F9.6)")  " elec_asym (LL,x=0.5)   = ", &
         elec_asym (pdf, x2, Q, alpha, .false.)
    write (u, "(1x,A,F9.6)")  " elec_asym (LL,x=0.9)   = ", &
         elec_asym (pdf, x3, Q, alpha, .false.)
    write (u, "(1x,A,F9.6)")  " elec_asym (LL,x=0.999) = ", &
         elec_asym (pdf, x4, Q, alpha, .false.)
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, 1)
    write (u, "(A)")
    write (u, "(1x,A,F9.6)")  " elec_asym (NLL,x=0.1)   = ", &
         elec_asym (pdf, x1, Q, alpha, .false.)
    write (u, "(1x,A,F9.6)")  " elec_asym (NLL,x=0.5)   = ", &
         elec_asym (pdf, x2, Q, alpha, .false.)
    write (u, "(1x,A,F9.6)")  " elec_asym (NLL,x=0.9)   = ", &
         elec_asym (pdf, x3, Q, alpha, .false.)
    write (u, "(1x,A,F9.6)")  " elec_asym (NLL,x=0.999) = ", &
         elec_asym (pdf, x4, Q, alpha, .false.)

    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, 1)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(A)")
    call pdf%write (u, with_qed = .true.)
    write (u, "(A)")

    write (u, "(A)")  "* Q = 10 GeV, elec_asym, LL+NLL, alpha running:"
    write (u, "(A)")

    write (u, "(1x,A,F9.6)")  " Integrator t (10 GeV)   = ", &
         t_alpha (pdf, Q)
    write (u, "(A)")

    write (u, "(1x,A,F9.6)")  " elec_asym (LL,x=0.1)   = ", &
         elec_asym (pdf, x1, Q, alpha, .true.)
    write (u, "(1x,A,F9.6)")  " elec_asym (LL,x=0.5)   = ", &
         elec_asym (pdf, x2, Q, alpha, .true.)
    write (u, "(1x,A,F9.6)")  " elec_asym (LL,x=0.9)   = ", &
         elec_asym (pdf, x3, Q, alpha, .true.)
    write (u, "(1x,A,F9.6)")  " elec_asym (LL,x=0.999) = ", &
         elec_asym (pdf, x4, Q, alpha, .true.)
    write (u, "(A)")

    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, 1)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F9.6)")  " elec_asym (NLL,x=0.1)   = ", &
         elec_asym (pdf, x1, Q, alpha, .true.)
    write (u, "(1x,A,F9.6)")  " elec_asym (NLL,x=0.5)   = ", &
         elec_asym (pdf, x2, Q, alpha, .true.)
    write (u, "(1x,A,F9.6)")  " elec_asym (NLL,x=0.9)   = ", &
         elec_asym (pdf, x3, Q, alpha, .true.)
    write (u, "(1x,A,F9.6)")  " elec_asym (NLL,x=0.999) = ", &
         elec_asym (pdf, x4, Q, alpha, .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, phot_asym, LL+NLL, alpha fixed:"
    write (u, "(A)")

    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, 1)
    write (u, "(1x,A,F9.6)")  " phot_asym (LL,x=0.1)   = ", &
         phot_asym (pdf, x1, Q, alpha, n_lep, .false.)
    write (u, "(1x,A,F9.6)")  " phot_asym (LL,x=0.5)   = ", &
         phot_asym (pdf, x2, Q, alpha, n_lep, .false.)
    write (u, "(1x,A,F9.6)")  " phot_asym (LL,x=0.9)   = ", &
         phot_asym (pdf, x3, Q, alpha, n_lep, .false.)
    write (u, "(1x,A,F9.6)")  " phot_asym (LL,x=0.999) = ", &
         phot_asym (pdf, x4, Q, alpha, n_lep, .false.)
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, 1)
    write (u, "(A)")
    write (u, "(1x,A,F9.6)")  " phot_asym (NLL,x=0.1)   = ", &
         phot_asym (pdf, x1, Q, alpha, n_lep, .false.)
    write (u, "(1x,A,F9.6)")  " phot_asym (NLL,x=0.5)   = ", &
         phot_asym (pdf, x2, Q, alpha, n_lep, .false.)
    write (u, "(1x,A,F9.6)")  " phot_asym (NLL,x=0.9)   = ", &
         phot_asym (pdf, x3, Q, alpha, n_lep, .false.)
    write (u, "(1x,A,F9.6)")  " phot_asym (NLL,x=0.999) = ", &
         phot_asym (pdf, x4, Q, alpha, n_lep, .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, phot_asym, LL+NLL, alpha running:"
    write (u, "(A)")

    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, 1)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F9.6)")  " phot_asym (LL,x=0.1)   = ", &
         phot_asym (pdf, x1, Q, alpha, n_lep, .true.)
    write (u, "(1x,A,F9.6)")  " phot_asym (LL,x=0.5)   = ", &
         phot_asym (pdf, x2, Q, alpha, n_lep, .true.)
    write (u, "(1x,A,F9.6)")  " phot_asym (LL,x=0.9)   = ", &
         phot_asym (pdf, x3, Q, alpha, n_lep, .true.)
    write (u, "(1x,A,F9.6)")  " phot_asym (LL,x=0.999) = ", &
         phot_asym (pdf, x4, Q, alpha, n_lep, .true.)
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, 1)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(A)")
    write (u, "(1x,A,F9.6)")  " phot_asym (NLL,x=0.1)   = ", &
         phot_asym (pdf, x1, Q, alpha, n_lep, .true.)
    write (u, "(1x,A,F9.6)")  " phot_asym (NLL,x=0.5)   = ", &
         phot_asym (pdf, x2, Q, alpha, n_lep, .true.)
    write (u, "(1x,A,F9.6)")  " phot_asym (NLL,x=0.9)   = ", &
         phot_asym (pdf, x3, Q, alpha, n_lep, .true.)
    write (u, "(1x,A,F9.6)")  " phot_asym (NLL,x=0.999) = ", &
         phot_asym (pdf, x4, Q, alpha, n_lep, .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: electron_pdfs_1"

  end subroutine electron_pdfs_1

  subroutine electron_pdfs_2 (u)
    integer, intent(in) :: u
    real(default) :: Q, alpha
    real(default) :: x1, x2, x3, x4, ln0
    real(default), dimension(6):: jll_nll
    logical, dimension(6) :: order
    integer :: n_lep
    Q = 10._default
    x1 = 0.1_default
    x2 = 0.5_default
    x3 = 0.9_default
    x4 = 0.999_default
    ln0 = 0._default
    n_lep = 1
    order = .true.

    write (u, "(A)")  "* Test output: electron_pdfs_2"
    write (u, "(A)")  "*   Purpose: check analytic properties"
    write (u, "(A)")

    write (u, "(A)")  "* Auxiliary functions II:"
    write (u, "(A)")

    write (u, "(A)")  "* Q = 10 GeV, elecbar_asym_p, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call elecbar_asym_p (x1, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " elecbar_asym_p (x=0.100)   = ", jll_nll
    call elecbar_asym_p (x2, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " elecbar_asym_p (x=0.500)   = ", jll_nll
    call elecbar_asym_p (x3, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " elecbar_asym_p (x=0.900)   = ", jll_nll
    call elecbar_asym_p (x4, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " elecbar_asym_p (x=0.999)   = ", jll_nll
    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, elecbar_asym_p, LL+NLL, alpha running:"
    write (u, "(A)")
    call elecbar_asym_p (x1, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " elecbar_asym_p (x=0.100)   = ", jll_nll
    call elecbar_asym_p (x2, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " elecbar_asym_p (x=0.500)   = ", jll_nll
    call elecbar_asym_p (x3, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " elecbar_asym_p (x=0.900)   = ", jll_nll
    call elecbar_asym_p (x4, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " elecbar_asym_p (x=0.999)   = ", jll_nll

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, photbar_asym_p, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call photbar_asym_p (x1, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " photbar_asym_p (x=0.100)   = ", jll_nll
    call photbar_asym_p (x2, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " photbar_asym_p (x=0.500)   = ", jll_nll
    call photbar_asym_p (x3, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " photbar_asym_p (x=0.900)   = ", jll_nll
    call photbar_asym_p (x4, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " photbar_asym_p (x=0.999)   = ", jll_nll
    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, photbar_asym_p, LL+NLL, alpha running:"
    write (u, "(A)")
    call photbar_asym_p (x1, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " photbar_asym_p (x=0.100)   = ", jll_nll
    call photbar_asym_p (x2, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " photbar_asym_p (x=0.500)   = ", jll_nll
    call photbar_asym_p (x3, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " photbar_asym_p (x=0.900)   = ", jll_nll
    call photbar_asym_p (x4, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " photbar_asym_p (x=0.999)   = ", jll_nll

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat_s, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call rechat_singlet (x1, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_singlet (x=0.100)   = ", jll_nll
    call rechat_singlet (x2, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_singlet (x=0.500)   = ", jll_nll
    call rechat_singlet (x3, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_singlet (x=0.900)   = ", jll_nll
    call rechat_singlet (x4, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_singlet (x=0.999)   = ", jll_nll
    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat_s, LL+NLL, alpha running:"
    write (u, "(A)")
    call rechat_singlet (x1, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_singlet (x=0.100)   = ", jll_nll
    call rechat_singlet (x2, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_singlet (x=0.500)   = ", jll_nll
    call rechat_singlet (x3, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_singlet (x=0.900)   = ", jll_nll
    call rechat_singlet (x4, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_singlet (x=0.999)   = ", jll_nll

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat_ns, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call rechat_nonsinglet (x1, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_nonsinglet (x=0.100)   = ", jll_nll
    call rechat_nonsinglet (x2, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_nonsinglet (x=0.500)   = ", jll_nll
    call rechat_nonsinglet (x3, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_nonsinglet (x=0.900)   = ", jll_nll
    call rechat_nonsinglet (x4, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_nonsinglet (x=0.999)   = ", jll_nll
    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat_ns, LL+NLL, alpha running:"
    write (u, "(A)")
    call rechat_nonsinglet (x1, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_nonsinglet (x=0.100)   = ", jll_nll
    call rechat_nonsinglet (x2, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_nonsinglet (x=0.500)   = ", jll_nll
    call rechat_nonsinglet (x3, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_nonsinglet (x=0.900)   = ", jll_nll
    call rechat_nonsinglet (x4, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_nonsinglet (x=0.999)   = ", jll_nll

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat_photon, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call rechat_photon (x1, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_photon (x=0.100)   = ", jll_nll
    call rechat_photon (x2, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_photon (x=0.500)   = ", jll_nll
    call rechat_photon (x3, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_photon (x=0.900)   = ", jll_nll
    call rechat_photon (x4, jll_nll, n_lep, ln0, order, running=.false.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_photon (x=0.999)   = ", jll_nll
    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat_photon, LL+NLL, alpha running:"
    write (u, "(A)")
    call rechat_photon (x1, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_photon (x=0.100)   = ", jll_nll
    call rechat_photon (x2, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_photon (x=0.500)   = ", jll_nll
    call rechat_photon (x3, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_photon (x=0.900)   = ", jll_nll
    call rechat_photon (x4, jll_nll, n_lep, ln0, order, running=.true.)
    write (u, "(1x,A,6(1x,ES11.4))")  " rechat_photon (x=0.999)   = ", jll_nll

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: electron_pdfs_2"

  end subroutine electron_pdfs_2

  subroutine electron_pdfs_3 (u)
    integer, intent(in) :: u
    type(qed_pdf_t) :: pdf
    real(default) :: Q, alpha
    real(default) :: x1, x2, x3, x4
    integer :: n_lep
    Q = 10._default
    x1 = 0.1_default
    x2 = 0.5_default
    x3 = 0.9_default
    x4 = 0.999_default
    alpha = ALPHA_QED_ME_REF
    n_lep = 1

    write (u, "(A)")  "* Test output: electron_pdfs_3"
    write (u, "(A)")  "*   Purpose: check analytic properties"
    write (u, "(A)")

    write (u, "(A)")  "* Auxiliary functions III:"
    write (u, "(A)")

    write (u, "(A)")  "* Q = 10 GeV, bar_asym, e+-, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,LL,x=0.100)   = ", &
         bar_asym (pdf, EPDF_ELE, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,LL,x=0.500)   = ", &
         bar_asym (pdf, EPDF_ELE, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,LL,x=0.900)   = ", &
         bar_asym (pdf, EPDF_ELE, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,LL,x=0.999)   = ", &
         bar_asym (pdf, EPDF_ELE, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,NLL,x=0.100)  = ", &
         bar_asym (pdf, EPDF_ELE, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,NLL,x=0.500)  = ", &
         bar_asym (pdf, EPDF_ELE, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,NLL,x=0.900)  = ", &
         bar_asym (pdf, EPDF_ELE, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,NLL,x=0.999)  = ", &
         bar_asym (pdf, EPDF_ELE, x4, Q, alpha, running=.false.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, bar_asym, e+-, LL+NLL, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,LL,x=0.100)   = ", &
         bar_asym (pdf, EPDF_ELE, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,LL,x=0.500)   = ", &
         bar_asym (pdf, EPDF_ELE, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,LL,x=0.900)   = ", &
         bar_asym (pdf, EPDF_ELE, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,LL,x=0.999)   = ", &
         bar_asym (pdf, EPDF_ELE, x4, Q, alpha, running=.true.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,NLL,x=0.100)  = ", &
         bar_asym (pdf, EPDF_ELE, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,NLL,x=0.500)  = ", &
         bar_asym (pdf, EPDF_ELE, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,NLL,x=0.900)  = ", &
         bar_asym (pdf, EPDF_ELE, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (ELE,NLL,x=0.999)  = ", &
         bar_asym (pdf, EPDF_ELE, x4, Q, alpha, running=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, bar_asym, gam, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,LL,x=0.100)   = ", &
         bar_asym (pdf, EPDF_G, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,LL,x=0.500)   = ", &
         bar_asym (pdf, EPDF_G, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,LL,x=0.900)   = ", &
         bar_asym (pdf, EPDF_G, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,LL,x=0.999)   = ", &
         bar_asym (pdf, EPDF_G, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,NLL,x=0.100)  = ", &
         bar_asym (pdf, EPDF_G, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,NLL,x=0.500)  = ", &
         bar_asym (pdf, EPDF_G, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,NLL,x=0.900)  = ", &
         bar_asym (pdf, EPDF_G, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,NLL,x=0.999)  = ", &
         bar_asym (pdf, EPDF_G, x4, Q, alpha, running=.false.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, bar_asym, gam, LL+NLL, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,LL,x=0.100)   = ", &
         bar_asym (pdf, EPDF_G, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,LL,x=0.500)   = ", &
         bar_asym (pdf, EPDF_G, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,LL,x=0.900)   = ", &
         bar_asym (pdf, EPDF_G, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,LL,x=0.999)   = ", &
         bar_asym (pdf, EPDF_G, x4, Q, alpha, running=.true.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,NLL,x=0.100)  = ", &
         bar_asym (pdf, EPDF_G, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,NLL,x=0.500)  = ", &
         bar_asym (pdf, EPDF_G, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,NLL,x=0.900)  = ", &
         bar_asym (pdf, EPDF_G, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " bar_asym (GAM,NLL,x=0.999)  = ", &
         bar_asym (pdf, EPDF_G, x4, Q, alpha, running=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, recbar, S, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    write (u, "(1x,A,F11.6)")  " recbar (S,LL,x=0.100)   = ", &
         recbar (pdf, EPDF_S, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (S,LL,x=0.500)   = ", &
         recbar (pdf, EPDF_S, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (S,LL,x=0.900)   = ", &
         recbar (pdf, EPDF_S, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (S,LL,x=0.999)   = ", &
         recbar (pdf, EPDF_S, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    write (u, "(1x,A,F11.6)")  " recbar (S,NLL,x=0.100)  = ", &
         recbar (pdf, EPDF_S, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (S,NLL,x=0.500)  = ", &
         recbar (pdf, EPDF_S, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (S,NLL,x=0.900)  = ", &
         recbar (pdf, EPDF_S, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (S,NLL,x=0.999)  = ", &
         recbar (pdf, EPDF_S, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, recbar, S, LL+NLL, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " recbar (S,LL,x=0.100)   = ", &
         recbar (pdf, EPDF_S, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (S,LL,x=0.500)   = ", &
         recbar (pdf, EPDF_S, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (S,LL,x=0.900)   = ", &
         recbar (pdf, EPDF_S, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (S,LL,x=0.999)   = ", &
         recbar (pdf, EPDF_S, x4, Q, alpha, running=.true.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " recbar (S,NLL,x=0.100)  = ", &
         recbar (pdf, EPDF_S, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (S,NLL,x=0.500)  = ", &
         recbar (pdf, EPDF_S, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (S,NLL,x=0.900)  = ", &
         recbar (pdf, EPDF_S, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (S,NLL,x=0.999)  = ", &
         recbar (pdf, EPDF_S, x4, Q, alpha, running=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, recbar, NS, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    write (u, "(1x,A,F11.6)")  " recbar (NS,LL,x=0.100)   = ", &
         recbar (pdf, EPDF_NS, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,LL,x=0.500)   = ", &
         recbar (pdf, EPDF_NS, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,LL,x=0.900)   = ", &
         recbar (pdf, EPDF_NS, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,LL,x=0.999)   = ", &
         recbar (pdf, EPDF_NS, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    write (u, "(1x,A,F11.6)")  " recbar (NS,NLL,x=0.100)  = ", &
         recbar (pdf, EPDF_NS, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,NLL,x=0.500)  = ", &
         recbar (pdf, EPDF_NS, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,NLL,x=0.900)  = ", &
         recbar (pdf, EPDF_NS, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,NLL,x=0.999)  = ", &
         recbar (pdf, EPDF_NS, x4, Q, alpha, running=.false.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, recbar, NS, LL+NLL, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,LL,x=0.100)   = ", &
         recbar (pdf, EPDF_NS, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,LL,x=0.500)   = ", &
         recbar (pdf, EPDF_NS, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,LL,x=0.900)   = ", &
         recbar (pdf, EPDF_NS, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,LL,x=0.999)   = ", &
         recbar (pdf, EPDF_NS, x4, Q, alpha, running=.true.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,NLL,x=0.100)  = ", &
         recbar (pdf, EPDF_NS, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,NLL,x=0.500)  = ", &
         recbar (pdf, EPDF_NS, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,NLL,x=0.900)  = ", &
         recbar (pdf, EPDF_NS, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (NS,NLL,x=0.999)  = ", &
         recbar (pdf, EPDF_NS, x4, Q, alpha, running=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, recbar, GAM, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,LL,x=0.100)   = ", &
         recbar (pdf, EPDF_G, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,LL,x=0.500)   = ", &
         recbar (pdf, EPDF_G, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,LL,x=0.900)   = ", &
         recbar (pdf, EPDF_G, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,LL,x=0.999)   = ", &
         recbar (pdf, EPDF_G, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,NLL,x=0.100)  = ", &
         recbar (pdf, EPDF_G, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,NLL,x=0.500)  = ", &
         recbar (pdf, EPDF_G, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,NLL,x=0.900)  = ", &
         recbar (pdf, EPDF_G, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,NLL,x=0.999)  = ", &
         recbar (pdf, EPDF_G, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, recbar, GAM, LL+NLL, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,LL,x=0.100)   = ", &
         recbar (pdf, EPDF_G, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,LL,x=0.500)   = ", &
         recbar (pdf, EPDF_G, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,LL,x=0.900)   = ", &
         recbar (pdf, EPDF_G, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,LL,x=0.999)   = ", &
         recbar (pdf, EPDF_G, x4, Q, alpha, running=.true.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,NLL,x=0.100)  = ", &
         recbar (pdf, EPDF_G, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,NLL,x=0.500)  = ", &
         recbar (pdf, EPDF_G, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,NLL,x=0.900)  = ", &
         recbar (pdf, EPDF_G, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " recbar (GAM,NLL,x=0.999)  = ", &
         recbar (pdf, EPDF_G, x4, Q, alpha, running=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: electron_pdfs_3"

  end subroutine electron_pdfs_3

  subroutine electron_pdfs_4 (u)
    integer, intent(in) :: u
    type(qed_pdf_t) :: pdf
    real(default) :: Q, alpha
    real(default) :: x1, x2, x3, x4
    integer :: n_lep
    Q = 10._default
    x1 = 0.1_default
    x2 = 0.5_default
    x3 = 0.9_default
    x4 = 0.999_default
    alpha = ALPHA_QED_ME_REF
    n_lep = 1

    write (u, "(A)")  "* Test output: electron_pdfs_4"
    write (u, "(A)")  "*   Purpose: check analytic properties"
    write (u, "(A)")

    write (u, "(A)")  "* Auxiliary functions IV:"
    write (u, "(A)")

    write (u, "(A)")  "* Q = 10 GeV, rechat, S, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    write (u, "(1x,A,F11.6)")  " rechat (S,LL,x=0.100)   = ", &
         rechat (pdf, EPDF_S, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (S,LL,x=0.500)   = ", &
         rechat (pdf, EPDF_S, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (S,LL,x=0.900)   = ", &
         rechat (pdf, EPDF_S, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (S,LL,x=0.999)   = ", &
         rechat (pdf, EPDF_S, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    write (u, "(1x,A,F11.6)")  " rechat (S,NLL,x=0.100)  = ", &
         rechat (pdf, EPDF_S, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (S,NLL,x=0.500)  = ", &
         rechat (pdf, EPDF_S, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (S,NLL,x=0.900)  = ", &
         rechat (pdf, EPDF_S, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (S,NLL,x=0.999)  = ", &
         rechat (pdf, EPDF_S, x4, Q, alpha, running=.false.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat, S, LL+NLL, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " rechat (S,LL,x=0.100)   = ", &
         rechat (pdf, EPDF_S, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (S,LL,x=0.500)   = ", &
         rechat (pdf, EPDF_S, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (S,LL,x=0.900)   = ", &
         rechat (pdf, EPDF_S, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (S,LL,x=0.999)   = ", &
         rechat (pdf, EPDF_S, x4, Q, alpha, running=.true.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " rechat (S,NLL,x=0.100)  = ", &
         rechat (pdf, EPDF_S, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (S,NLL,x=0.500)  = ", &
         rechat (pdf, EPDF_S, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (S,NLL,x=0.900)  = ", &
         rechat (pdf, EPDF_S, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (S,NLL,x=0.999)  = ", &
         rechat (pdf, EPDF_S, x4, Q, alpha, running=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat, NS, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    write (u, "(1x,A,F11.6)")  " rechat (NS,LL,x=0.100)   = ", &
         rechat (pdf, EPDF_NS, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,LL,x=0.500)   = ", &
         rechat (pdf, EPDF_NS, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,LL,x=0.900)   = ", &
         rechat (pdf, EPDF_NS, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,LL,x=0.999)   = ", &
         rechat (pdf, EPDF_NS, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    write (u, "(1x,A,F11.6)")  " rechat (NS,NLL,x=0.100)  = ", &
         rechat (pdf, EPDF_NS, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,NLL,x=0.500)  = ", &
         rechat (pdf, EPDF_NS, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,NLL,x=0.900)  = ", &
         rechat (pdf, EPDF_NS, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,NLL,x=0.999)  = ", &
         rechat (pdf, EPDF_NS, x4, Q, alpha, running=.false.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat, NS, LL+NLL, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,LL,x=0.100)   = ", &
         rechat (pdf, EPDF_NS, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,LL,x=0.500)   = ", &
         rechat (pdf, EPDF_NS, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,LL,x=0.900)   = ", &
         rechat (pdf, EPDF_NS, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,LL,x=0.999)   = ", &
         rechat (pdf, EPDF_NS, x4, Q, alpha, running=.true.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,NLL,x=0.100)  = ", &
         rechat (pdf, EPDF_NS, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,NLL,x=0.500)  = ", &
         rechat (pdf, EPDF_NS, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,NLL,x=0.900)  = ", &
         rechat (pdf, EPDF_NS, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (NS,NLL,x=0.999)  = ", &
         rechat (pdf, EPDF_NS, x4, Q, alpha, running=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat, GAM, LL+NLL, alpha fixed:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,LL,x=0.100)   = ", &
         rechat (pdf, EPDF_G, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,LL,x=0.500)   = ", &
         rechat (pdf, EPDF_G, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,LL,x=0.900)   = ", &
         rechat (pdf, EPDF_G, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,LL,x=0.999)   = ", &
         rechat (pdf, EPDF_G, x4, Q, alpha, running=.false.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,NLL,x=0.100)  = ", &
         rechat (pdf, EPDF_G, x1, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,NLL,x=0.500)  = ", &
         rechat (pdf, EPDF_G, x2, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,NLL,x=0.900)  = ", &
         rechat (pdf, EPDF_G, x3, Q, alpha, running=.false.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,NLL,x=0.999)  = ", &
         rechat (pdf, EPDF_G, x4, Q, alpha, running=.false.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, rechat, GAM, LL+NLL, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 0, n_lep)
    call pdf%allocate_aqed (order = 0, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,LL,x=0.100)   = ", &
         rechat (pdf, EPDF_G, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,LL,x=0.500)   = ", &
         rechat (pdf, EPDF_G, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,LL,x=0.900)   = ", &
         rechat (pdf, EPDF_G, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,LL,x=0.999)   = ", &
         rechat (pdf, EPDF_G, x4, Q, alpha, running=.true.)
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,NLL,x=0.100)  = ", &
         rechat (pdf, EPDF_G, x1, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,NLL,x=0.500)  = ", &
         rechat (pdf, EPDF_G, x2, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,NLL,x=0.900)  = ", &
         rechat (pdf, EPDF_G, x3, Q, alpha, running=.true.)
    write (u, "(1x,A,F11.6)")  " rechat (GAM,NLL,x=0.999)  = ", &
         rechat (pdf, EPDF_G, x4, Q, alpha, running=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: electron_pdfs_4"

  end subroutine electron_pdfs_4

  subroutine electron_pdfs_5 (u)
    integer, intent(in) :: u
    type(qed_pdf_t) :: pdf
    real(default) :: Q, alpha
    real(default) :: x1, x2, x3, x4
    integer :: n_lep
    Q = 10._default
    x1 = 0.1_default
    x2 = 0.5_default
    x3 = 0.9_default
    x4 = 0.999_default
    alpha = ALPHA_QED_ME_REF
    n_lep = 1

    write (u, "(A)")  "* Test output: electron_pdfs_5"
    write (u, "(A)")  "*   Purpose: check analytic properties"
    write (u, "(A)")

    write (u, "(A)")  "* Auxiliary functions V:"
    write (u, "(A)")

    write (u, "(A)")  "* Integrals over endpoint_func_NS, interval [0,1]:"
    write (u, "(A)")
    write (u, "(1x,A,F11.6)")  " endpoint_func_NS (0.100) = ", &
         endpoint_func_NS (x1)
    write (u, "(1x,A,F11.6)")  " endpoint_func_NS (0.500) = ", &
         endpoint_func_NS (x2)
    write (u, "(1x,A,F11.6)")  " endpoint_func_NS (0.900) = ", &
         endpoint_func_NS (x3)
    write (u, "(1x,A,F11.6)")  " endpoint_func_NS (0.999) = ", &
         endpoint_func_NS (x4)

    write (u, "(A)")
    write (u, "(A)")  "* Integrals over endpoint_func_S, interval [0,1]:"
    write (u, "(A)")
    write (u, "(1x,A,F11.6)")  " endpoint_func_S (0.100) = ", &
         endpoint_func_S (x1, n_lep)
    write (u, "(1x,A,F11.6)")  " endpoint_func_S (0.500) = ", &
         endpoint_func_S (x2, n_lep)
    write (u, "(1x,A,F11.6)")  " endpoint_func_S (0.900) = ", &
         endpoint_func_S (x3, n_lep)
    write (u, "(1x,A,F11.6)")  " endpoint_func_S (0.999) = ", &
         endpoint_func_S (x4, n_lep)

    write (u, "(A)")
    write (u, "(A)")  "* Integrals over endpoint_func_GAM, interval [0,1]:"
    write (u, "(A)")
    write (u, "(1x,A,F11.6)")  " endpoint_func_GAM (0.100) = ", &
         endpoint_func_GAM (x1)
    write (u, "(1x,A,F11.6)")  " endpoint_func_GAM (0.500) = ", &
         endpoint_func_GAM (x2)
    write (u, "(1x,A,F11.6)")  " endpoint_func_GAM (0.900) = ", &
         endpoint_func_GAM (x3)
    write (u, "(1x,A,F11.6)")  " endpoint_func_GAM (0.999) = ", &
         endpoint_func_GAM (x4)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, recnum, alpha fixed:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    write (u, "(1x,A,ES11.4)")  " recnum (S,   0.100) = ", &
         rec_num (pdf, EPDF_S, x1, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (S,   0.500) = ", &
         rec_num (pdf, EPDF_S, x2, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (S,   0.900) = ", &
         rec_num (pdf, EPDF_S, x3, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (S,   0.999) = ", &
         rec_num (pdf, EPDF_S, x4, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (NS,  0.100) = ", &
         rec_num (pdf, EPDF_NS, x1, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (NS,  0.500) = ", &
         rec_num (pdf, EPDF_NS, x2, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (NS,  0.900) = ", &
         rec_num (pdf, EPDF_NS, x3, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (NS,  0.999) = ", &
         rec_num (pdf, EPDF_NS, x4, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (GAM, 0.100) = ", &
         rec_num (pdf, EPDF_G, x1, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (GAM, 0.500) = ", &
         rec_num (pdf, EPDF_G, x2, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (GAM, 0.900) = ", &
         rec_num (pdf, EPDF_G, x3, Q, alpha, .false.)
    write (u, "(1x,A,ES11.4)")  " recnum (GAM, 0.999) = ", &
         rec_num (pdf, EPDF_G, x4, Q, alpha, .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, recnum, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (S,   0.100) = ", &
         rec_num (pdf, EPDF_S, x1, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (S,   0.500) = ", &
         rec_num (pdf, EPDF_S, x2, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (S,   0.900) = ", &
         rec_num (pdf, EPDF_S, x3, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (S,   0.999) = ", &
         rec_num (pdf, EPDF_S, x4, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (NS,  0.100) = ", &
         rec_num (pdf, EPDF_NS, x1, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (NS,  0.500) = ", &
         rec_num (pdf, EPDF_NS, x2, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (NS,  0.900) = ", &
         rec_num (pdf, EPDF_NS, x3, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (NS,  0.999) = ", &
         rec_num (pdf, EPDF_NS, x4, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (GAM, 0.100) = ", &
         rec_num (pdf, EPDF_G, x1, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (GAM, 0.500) = ", &
         rec_num (pdf, EPDF_G, x2, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (GAM, 0.900) = ", &
         rec_num (pdf, EPDF_G, x3, Q, alpha, .true.)
    write (u, "(1x,A,ES11.4)")  " recnum (GAM, 0.999) = ", &
         rec_num (pdf, EPDF_G, x4, Q, alpha, .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: electron_pdfs_5"

  end subroutine electron_pdfs_5

  subroutine electron_pdfs_6 (u)
    integer, intent(in) :: u
    type(qed_pdf_t) :: pdf
    real(default) :: Q, alpha
    real(default), dimension(12) :: x
    integer :: n_lep
    Q = 10._default
    x( 1) = 0.1_default
    x( 2) = 0.2_default
    x( 3) = 0.3_default
    x( 4) = 0.4_default    
    x( 5) = 0.5_default
    x( 6) = 0.6_default
    x( 7) = 0.7_default
    x( 8) = 0.8_default
    x( 9) = 0.9_default
    x(10) = 0.95_default
    x(11) = 0.99_default
    x(12) = 0.999_default
    alpha = ALPHA_QED_ME_REF
    n_lep = 1

    write (u, "(A)")  "* Test output: electron_pdfs_6"
    write (u, "(A)")  "*   Purpose: full electron PDFs"
    write (u, "(A)")

    write (u, "(A)")  "* Full NLL electron PDFs:"
    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, NLL, alpha fixed:"
    write (u, "(A)")
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)    
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.100, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(1), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(1), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(1), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.200, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(2), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(2), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(2), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.300, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(3), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(3), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(3), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.400, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(4), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(4), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(4), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.500, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(5), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(5), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(5), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.600, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(6), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(6), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(6), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.700, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(7), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(7), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(7), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.800, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(8), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(8), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(8), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.900, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(9), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(9), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(9), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.950, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(10), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(10), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(10), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.990, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(11), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(11), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(11), Q, alpha, .false., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.999, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(12), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_G, x(12), Q, alpha, .false., .true.), &
         elec_pdf (pdf, EPDF_NS, x(12), Q, alpha, .false., .true.)
    write (u, "(A)")
    write (u, "(A)")  "* Q = 10 GeV, NLL, alpha running:"
    write (u, "(A)")
    call pdf%init (ME_REF, ALPHA_QED_ME_REF, -1._default, &
         3000._default, 3, 1, n_lep)
    call pdf%allocate_aqed (order = 1, n_f = 0, n_lep = 1, running = .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.100, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(1), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(1), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(1), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.200, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(2), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(2), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(2), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.300, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(3), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(3), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(3), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.400, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(4), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(4), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(4), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.500, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(5), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(5), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(5), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.600, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(6), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(6), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(6), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.700, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(7), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(7), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(7), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.800, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(8), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(8), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(8), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.900, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(9), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(9), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(9), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.950, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(10), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(10), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(10), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.990, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(11), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(11), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(11), Q, alpha, .true., .true.)
    write (u, "(1x,A,3(1x,F11.6))")  " ePDF (x = 0.999, S/GAM/NS) = ", &
         elec_pdf (pdf, EPDF_S, x(12), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_G, x(12), Q, alpha, .true., .true.), &
         elec_pdf (pdf, EPDF_NS, x(12), Q, alpha, .true., .true.)
    write (u, "(A)")
    write (u, "(A)")  "* Check singlet-nonsinglet linear combination"
    write (u, "(A)")
    write (u, "(1x,A,F11.6)")  " ePDF (x = 0.950, e- - [S + NS]/2) = ", &
         abs(elec_pdf (pdf, EPDF_ELE, x(10), Q, alpha, .true., .true.) - &
         (elec_pdf (pdf, EPDF_S, x(10), Q, alpha, .true., .true.) + &
          elec_pdf (pdf, EPDF_NS, x(10), Q, alpha, .true., .true.))/two)
    write (u, "(1x,A,F11.6)")  " ePDF (x = 0.950, e+ - [S - NS]/2) = ", &
         abs(elec_pdf (pdf, EPDF_POS, x(10), Q, alpha, .true., .true.) - &
         (elec_pdf (pdf, EPDF_S, x(10), Q, alpha, .true., .true.) - &
          elec_pdf (pdf, EPDF_NS, x(10), Q, alpha, .true., .true.))/two)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: electron_pdfs_6"

  end subroutine electron_pdfs_6

end module electron_pdfs_uti
