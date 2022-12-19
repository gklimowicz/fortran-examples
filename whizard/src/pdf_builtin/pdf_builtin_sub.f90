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

submodule (pdf_builtin) pdf_builtin_s

  use string_utils, only: upper_case
  use constants, only: PI
  use diagnostics
  use mrst2004qed
  use cteq6pdf
  use mstwpdf
  use ct10pdf
  use CJ_pdf
  use ct14pdf
  use ct18pdf

  implicit none
  save

! The available sets
  integer, parameter :: nsets = 24
  integer, parameter :: &
       CTEQ6M = 1, CTEQ6D = 2, CTEQ6L = 3, CTEQ6L1 = 4, &
       MRST2004QEDp = 5, MRST2004QEDn = 6, MSTW2008LO = 7, MSTW2008NLO = 8, &
       MSTW2008NNLO = 9, CT10 = 10, CJ12_max = 11, CJ12_mid = 12, &
       CJ12_min = 13, MMHT2014LO = 14, MMHT2014NLO = 15, MMHT2014NNLO = 16, &
       CT14LL = 17, CT14L = 18, CT14N = 19, CT14NN = 20, CJ15_LO = 21, &
       CJ15_NLO = 22, CT18N = 23, CT18NN = 24

! Limits
  real(kind=default), parameter :: &
       cteq6_q_min = 1.3, cteq6_q_max = 10.E3, &
       cteq6_x_min = 1.E-6, cteq6_x_max = 1., &
       mrst2004qed_q_min = sqrt (1.26_default), mrst2004qed_q_max = sqrt (0.99E7_default), &
       mrst2004qed_x_min = 1.01E-5, mrst2004qed_x_max = 1., &
       mstw2008_q_min = 1.01, mstw2008_q_max = sqrt (0.99E9_default), &
       mstw2008_x_min = 1.01E-6, mstw2008_x_max = 1., &
       ct10_q_min = 1.3, ct10_q_max = 1.E5, &
       ct10_x_min = 1.E-8, ct10_x_max = 1., &
       cj12_q_min = 1.3, cj12_q_max = 1.E5, &
       cj12_x_min = 1.E-6, cj12_x_max = 1., &
       mmht2014_q_min = mstw2008_q_min, mmht2014_q_max = mstw2008_q_max, &
       mmht2014_x_min = mstw2008_x_min, mmht2014_x_max = mstw2008_x_max, &
       ct14_q_min = 1.3, ct14_q_max = 1.E5, ct14_x_min = 1.E-9, &
       ct14_x_max = 1., ct18_q_min = 1.3, ct18_q_max = 1.E5, &
       ct18_x_min = 1.E-9, ct18_x_max = 1.

! Lambda_QCD and quark masses

  real(kind=double), dimension(6), parameter :: alam_ct10 = &
       [  0.37219423568859566_double, 0.37219423568859566_double, &
          0.37219423568859566_double, 0.37219423568859566_double, &
          0.32560730624033124_double, 0.22600000000000001_double ]

  real(kind=double), dimension(6), parameter :: alam_cteq6m = &
       [  0.37248648555333957_double, 0.37248648555333957_double, &
          0.37248648555333957_double, 0.37248648555333957_double, &
          0.32590307925177625_double, 0.22623045868581182_double ]

  real(kind=double), dimension(6), parameter :: alam_cteq6l = &
       [  0.37218603304249098_double, 0.37218603304249098_double, &
          0.37218603304249098_double, 0.37218603304249098_double, &
          0.32559900534407232_double, 0.22599353258372407_double ]

  real(kind=double), dimension(6), parameter :: alam_cteq6ll = &
       [  0.24560434095679567_double, 0.24560434095679567_double, &
          0.24560434095679567_double, 0.24560434095679567_double, &
          0.21495099184302788_double, 0.16499885346541945_double ]

  real(kind=double), dimension(0:5), parameter :: amhat = &
       [  0.0000000000000000_double,  0.0000000000000000_double, &
          0.0000000000000000_double,  0.0000000000000000_double, &
          1.3000000000000000_double,  4.500000000000000_double ]

  !!! Global parameter for minimal value of evolution
  real(kind=double), parameter :: amn = 0.375_double, &
    tol = .0005_double

  integer, parameter :: nhq = 2

  real(kind=double), parameter :: &
       ca = 3.D0, cf = 4./3.D0, tr = 0.5D0
  real(kind=double), parameter :: &
       b00 =  11./3.d0 * ca, b01 =  -4./3.d0 * tr, b10 =  34./3.d0 * ca**2, &
       b11 = -20./3.d0 * ca*tr - 4.* cf*tr

  ! Init flags
  integer :: cteq6_initialized = -1
  integer :: mstw2008_initialized = -1
  integer :: mmht2014_initialized = -1
  logical :: ct10_initialized = .false.
  integer :: cj12_initialized = -1
  integer :: ct14_initialized = -1
  integer :: cj15_initialized = -1
  integer :: ct18_initialized = -1
  logical :: &
       mrst2004qedp_initialized =  .false., &
       mrst2004qedn_initialized =  .false.
  type(string_t) :: mrst2004qedp_prefix, mrst2004qedn_prefix, &
       mstw2008_prefix, mmht2014_prefix

contains

! Get PDF name
  module function pdf_get_name (pdftype) result (name)
    integer, intent(in) :: pdftype
    type(string_t) :: name
    select case (pdftype)
    case (CTEQ6M)
       name = var_str ("CTEQ6M")
    case (CTEQ6D)
       name = var_str ("CTEQ6D")
    case (CTEQ6L)
       name = var_str ("CTEQ6L")
    case (CTEQ6L1)
       name = var_str ("CTEQ6L1")
    case (MRST2004QEDp)
       name = var_str ("MRST2004QEDp")
    case (MRST2004QEDn)
       name = var_str ("MRST2004QEDn")
    case (MSTW2008LO)
       name = var_str ("MSTW2008LO")
    case (MSTW2008NLO)
       name = var_str ("MSTW2008NLO")
    case (MSTW2008NNLO)
       name = var_str ("MSTW2008NNLO")
    case (CT10)
       name = var_str ("CT10")
    case (CJ12_max)
       name = var_str ("CJ12_max")
    case (CJ12_mid)
       name = var_str ("CJ12_mid")
    case (CJ12_min)
       name = var_str ("CJ12_min")
    case (CJ15_LO)
       name = var_str ("CJ15LO")
    case (CJ15_NLO)
       name = var_str ("CJ15NLO")
    case (MMHT2014LO)
       name = var_str ("MMHT2014LO")
    case (MMHT2014NLO)
       name = var_str ("MMHT2014NLO")
    case (MMHT2014NNLO)
       name = var_str ("MMHT2014NNLO")
    case (CT14LL)
       name = var_str ("CT14LL")
    case (CT14L)
       name = var_str ("CT14L")
    case (CT14N)
       name = var_str ("CT14N")
    case (CT14NN)
       name = var_str ("CT14NN")
    case (CT18N)
       name = var_str ("CT18N")
    case (CT18NN)
       name = var_str ("CT18NN")
    case default
       call msg_fatal ("pdf_builtin: internal: invalid PDF set!")
    end select
  end function pdf_get_name

! Get the ID of a PDF set
  module function pdf_get_id (name) result (id)
    type(string_t), intent(in) :: name
    integer :: id
    do id = 1, nsets
       if (upper_case (pdf_get_name (id)) == upper_case (name)) return
    end do
    id = -1
  end function pdf_get_id

! Query whether a PDF supplies a photon distribution
  module function pdf_provides_photon (pdftype) result (flag)
    integer, intent(in) :: pdftype
    logical :: flag
    select case (pdftype)
    case (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1)
       flag = .false.
    case (MRST2004QEDp, MRST2004QEDn)
       flag = .true.
    case (MSTW2008LO, MSTW2008NLO, MSTW2008NNLO)
       flag = .false.
    case (CT10)
       flag = .false.
    case (CJ12_max, CJ12_mid, CJ12_min)
       flag = .false.
    case (CJ15_LO, CJ15_NLO)
       flag = .false.
    case (MMHT2014LO, MMHT2014NLO, MMHT2014NNLO)
       flag = .false.
    case (CT14LL, CT14L, CT14N, CT14NN)
       flag = .false.
    case (CT18N, CT18NN)
       flag = .false.
    case default
       call msg_fatal ("pdf_builtin: internal: invalid PDF set!")
    end select
  end function pdf_provides_photon

! Initialize a PDF
  module subroutine pdf_init (pdftype, prefix, verbose)
    integer, intent(in) :: pdftype
    type(string_t), intent(in), optional :: prefix
    type(string_t) :: mprefix
    logical, intent(in), optional :: verbose
    logical :: mverbose
    if (present (prefix)) then
       mprefix = prefix
    else
       mprefix = ""
    end if
    if (present (verbose)) then
       mverbose = verbose
    else
       mverbose = .true.
    end if
    select case (pdftype)
    case (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1)
       if (cteq6_initialized == pdftype) return
       call setctq6 (pdftype, char (mprefix))
       cteq6_initialized = pdftype
    case (MRST2004QEDp)
       if (mrst2004qedp_initialized) return
       mrst2004qedp_initialized = .true.
       mrst2004qedp_prefix = mprefix
    case (MRST2004QEDn)
       if (mrst2004qedn_initialized) return
       mrst2004qedn_initialized = .true.
       mrst2004qedn_prefix = mprefix
    case (MSTW2008LO, MSTW2008NLO, MSTW2008NNLO)
       if (mstw2008_initialized == pdftype) return
       mstw2008_initialized = pdftype
       mstw2008_prefix = mprefix
    case (CT10)
       if (ct10_initialized) return
       call setct10 (char (mprefix), 100)
       ct10_initialized = .true.
    case (CJ12_max)
       if (cj12_initialized == pdftype) return
       call setCJ (char (mprefix), 300)
       cj12_initialized = pdftype
    case (CJ12_mid)
       if (cj12_initialized == pdftype) return
       call setCJ (char (mprefix), 200)
       cj12_initialized = pdftype
    case (CJ12_min)
       if (cj12_initialized == pdftype) return
       call setCJ (char (mprefix), 100)
       cj12_initialized = pdftype
    case (CJ15_LO)
       if (cj15_initialized == pdftype) return
       call setCJ (char (mprefix), 400)
       cj15_initialized = pdftype
    case (CJ15_NLO)
       if (cj15_initialized == pdftype) return
       call setCJ (char (mprefix), 500)
       cj15_initialized = pdftype
    case (MMHT2014LO, MMHT2014NLO, MMHT2014NNLO)
       if (mmht2014_initialized == pdftype) return
       mmht2014_initialized = pdftype
       mmht2014_prefix = mprefix
    case (CT14LL)
       if (ct14_initialized == pdftype) return
       ct14_initialized = pdftype
       call setct14 (char (mprefix), 1)
    case (CT14L)
       if (ct14_initialized == pdftype) return
       ct14_initialized = pdftype
       call setct14 (char (mprefix), 2)
    case (CT14N)
       if (ct14_initialized == pdftype) return
       ct14_initialized = pdftype
       call setct14 (char (mprefix), 3)
    case (CT14NN)
       if (ct14_initialized == pdftype) return
       ct14_initialized = pdftype
       call setct14 (char (mprefix), 4)
    case (CT18N)
       if (ct18_initialized == pdftype) return
       ct18_initialized = pdftype
       call setct18 (char (prefix), 1)
    case (CT18NN)
       if (ct18_initialized == pdftype) return
       ct18_initialized = pdftype
       call setct18 (char (prefix), 2)
    case default
       call msg_fatal ("pdf_builtin: internal: invalid PDF set!")
    end select
    if (mverbose) call msg_message ("Initialized builtin PDF " // &
         char (pdf_get_name (pdftype)))
    !!! Up to now only proton
  end subroutine pdf_init

! Evolve PDF
  module subroutine pdf_evolve (pdftype, x, q, f, fphoton)
    integer, intent(in) :: pdftype
    real(kind=default), intent(in) :: x, q
    real(kind=default), intent(out), optional :: f(-6:6), fphoton
    real(kind=double) :: mx, mq
    real(kind=double) :: upv, dnv, ups, dns, str, chm, bot, glu, phot, &
         sbar, bbar, cbar
    type(string_t) :: setname
    select case (pdftype)
    case (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1)
       if (cteq6_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (cteq6_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (cteq6_initialized)) // &
            " cannot be used simultaneously")
       mx = max (min (x, cteq6_x_max), cteq6_x_min)
       mq = max (min (q, cteq6_q_max), cteq6_q_min)
       if (present (f)) f = [ 0._double, &
            ctq6pdf (-5, mx, mq), ctq6pdf (-4, mx, mq), &
            ctq6pdf (-3, mx, mq), ctq6pdf (-1, mx, mq), &
            ctq6pdf (-2, mx, mq), ctq6pdf ( 0, mx, mq), &
            ctq6pdf ( 2, mx, mq), ctq6pdf ( 1, mx, mq), &
            ctq6pdf ( 3, mx, mq), ctq6pdf ( 4, mx, mq), &
            ctq6pdf ( 5, mx, mq), 0._double ]
       if (present (fphoton)) call msg_fatal ("photon pdf requested for " // &
            char (pdf_get_name (pdftype)) // " which does not provide it!")
    case (MRST2004QEDp)
       if (.not. mrst2004qedp_initialized) call msg_fatal ( &
            "pdf_builtin: internal: PDF set MRST2004QEDp requested without " // &
            "initialization!")
       mx = max (min (x, mrst2004qed_x_max), mrst2004qed_x_min)
       mq = max (min (q, mrst2004qed_q_max), mrst2004qed_q_min)
       call mrstqed (mx, mq, 1, upv, dnv, ups, dns, str, chm, bot, glu, phot, &
            char (mrst2004qedp_prefix))
       if (present (f)) f = &
            [ 0._double, bot, chm, str, ups, dns, glu, dns + dnv, ups + upv, &
            str, chm, bot, 0._double ] / mx
       if (present (fphoton)) fphoton = phot / mx
    case (MRST2004QEDn)
       if (.not. mrst2004qedn_initialized) call msg_fatal ( &
            "pdf_builtin: internal: PDF set MRST2004QEDn requested without " // &
            "initialization!")
       mx = max (min (x, mrst2004qed_x_max), mrst2004qed_x_min)
       mq = max (min (q, mrst2004qed_q_max), mrst2004qed_q_min)
       call mrstqed (mx, mq, 2, upv, dnv, ups, dns, str, chm, bot, glu, phot, &
            char (mrst2004qedn_prefix))
       if (present (f)) f = &
            [ 0._double, bot, chm, str, ups, dns, glu, dns + dnv, ups + upv, &
            str, chm, bot, 0._double ] / mx
       if (present (fphoton)) fphoton = phot / mx
    case (MSTW2008LO, MSTW2008NLO, MSTW2008NNLO)
       if (mstw2008_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (mstw2008_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (mstw2008_initialized)) // &
            " cannot be used simultaneously")
       select case (pdftype)
       case (MSTW2008LO)
          setname = var_str ("mstw2008lo")
       case (MSTW2008NLO)
          setname = var_str ("mstw2008nlo")
       case (MSTW2008NNLO)
          setname = var_str ("mstw2008nnlo")
       end select
       mx = max (min (x, mstw2008_x_max), mstw2008_x_min)
       mq = max (min (q, mstw2008_q_max), mstw2008_q_min)
       call getmstwpdf (char (mstw2008_prefix), char (setname), 0, mx, mq, &
            upv, dnv, ups, dns, str, sbar, chm, cbar, bot, bbar, glu, phot)
       if (present (f)) f = &
            [ 0._double, bbar, cbar, sbar, ups, dns, glu, dns + dnv, ups + upv, &
            str, chm, bot, 0._double ] / mx
       if (present (fphoton)) call msg_fatal ("photon pdf requested for " // &
            char (pdf_get_name (pdftype)) // " which does not provide it!")
    case (CT10)
       if (.not. ct10_initialized) call msg_fatal ( &
            "pdf_builtin: internal: PDF set CT10 requested without " // &
            "initialization!")
       mx = max (min (x, ct10_x_max), ct10_x_min)
       mq = max (min (q, ct10_q_max), ct10_q_min)
       if (present (f)) f = [ 0._double, &
            getct10pdf (-5, mx, mq), getct10pdf (-4, mx, mq), &
            getct10pdf (-3, mx, mq), getct10pdf (-1, mx, mq), &
            getct10pdf (-2, mx, mq), getct10pdf ( 0, mx, mq), &
            getct10pdf ( 2, mx, mq), getct10pdf ( 1, mx, mq), &
            getct10pdf ( 3, mx, mq), getct10pdf ( 4, mx, mq), &
            getct10pdf ( 5, mx, mq), 0._double ]
       if (present (fphoton)) call msg_fatal ("photon pdf requested for " // &
            "CT10 which does not provide it!")
    case (CJ12_max, CJ12_mid, CJ12_min)
       if (cj12_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (cj12_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (cj12_initialized)) // &
            " cannot be used simultaneously")
       mx = max (min (x, cj12_x_max), cj12_x_min)
       mq = max (min (q, cj12_q_max), cj12_q_min)
       if (present (f)) f = [ 0._double, &
            CJpdf (-5, mx, mq), CJpdf (-4, mx, mq), CJpdf (-3, mx, mq), &
            CJpdf (-2, mx, mq), CJpdf (-1, mx, mq), CJpdf ( 0, mx, mq), &
            CJpdf ( 1, mx, mq), CJpdf ( 2, mx, mq), CJpdf ( 3, mx, mq), &
            CJpdf ( 4, mx, mq), CJpdf ( 5, mx, mq), 0._double ]
       if (present (fphoton)) call msg_fatal ("photon pdf requested for " // &
            "CJ12 which does not provide it!")
    case (CJ15_LO, CJ15_NLO)
       if (cj15_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (cj15_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (cj15_initialized)) // &
            " cannot be used simultaneously")
       mx = max (min (x, cj12_x_max), cj12_x_min)
       mq = max (min (q, cj12_q_max), cj12_q_min)
       if (present (f)) f = [ 0._double, &
            CJpdf (-5, mx, mq), CJpdf (-4, mx, mq), CJpdf (-3, mx, mq), &
            CJpdf (-2, mx, mq), CJpdf (-1, mx, mq), CJpdf ( 0, mx, mq), &
            CJpdf ( 1, mx, mq), CJpdf ( 2, mx, mq), CJpdf ( 3, mx, mq), &
            CJpdf ( 4, mx, mq), CJpdf ( 5, mx, mq), 0._double ]
       if (present (fphoton)) call msg_fatal ("photon pdf requested for " // &
            "CJ15 which does not provide it!")
    case (MMHT2014LO, MMHT2014NLO, MMHT2014NNLO)
       if (mmht2014_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (mmht2014_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (mmht2014_initialized)) // &
            " cannot be used simultaneously")
       select case (pdftype)
       case (MMHT2014LO)
          setname = var_str ("mmht2014lo")
       case (MMHT2014NLO)
          setname = var_str ("mmht2014nlo")
       case (MMHT2014NNLO)
          setname = var_str ("mmht2014nnlo")
       end select
       mx = max (min (x, mmht2014_x_max), mmht2014_x_min)
       mq = max (min (q, mmht2014_q_max), mmht2014_q_min)
       call getmstwpdf (char (mmht2014_prefix), char (setname), 0, mx, mq, &
            upv, dnv, ups, dns, str, sbar, chm, cbar, bot, bbar, glu, phot)
       if (present (f)) f = &
            [ 0._double, bbar, cbar, sbar, ups, dns, glu, dns + dnv, ups + upv, &
            str, chm, bot, 0._double ] / mx
       if (present (fphoton)) call msg_fatal ("photon pdf requested for " // &
            char (pdf_get_name (pdftype)) // " which does not provide it!")
    case (CT14LL, CT14L, CT14N, CT14NN)
       if (ct14_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (ct14_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (ct14_initialized)) // &
            " cannot be used simultaneously")
       mx = max (min (x, ct14_x_max), ct14_x_min)
       mq = max (min (q, ct14_q_max), ct14_q_min)
       if (present (f)) f = [ 0._double, &
            getct14pdf (-5, mx, mq), getct14pdf (-4, mx, mq), &
            getct14pdf (-3, mx, mq), getct14pdf (-1, mx, mq), &
            getct14pdf (-2, mx, mq), getct14pdf ( 0, mx, mq), &
            getct14pdf ( 2, mx, mq), getct14pdf ( 1, mx, mq), &
            getct14pdf ( 3, mx, mq), getct14pdf ( 4, mx, mq), &
            getct14pdf ( 5, mx, mq), 0._double ]
    case (CT18N, CT18NN)
       if (ct18_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (ct18_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (ct18_initialized)) // &
            " cannot be used simultaneously")
       mx = max (min (x, ct18_x_max), ct18_x_min)
       mq = max (min (q, ct18_q_max), ct18_q_min)
       if (present (f)) f = [ 0._double, &
            getct18pdf (-5, mx, mq), getct18pdf (-4, mx, mq), &
            getct18pdf (-3, mx, mq), getct18pdf (-1, mx, mq), &
            getct18pdf (-2, mx, mq), getct18pdf ( 0, mx, mq), &
            getct18pdf ( 2, mx, mq), getct18pdf ( 1, mx, mq), &
            getct18pdf ( 3, mx, mq), getct18pdf ( 4, mx, mq), &
            getct18pdf ( 5, mx, mq), 0._double ]
       if (present (fphoton)) call msg_fatal ("photon pdf requested for " // &
            char (pdf_get_name (pdftype)) // " which does not provide it!")
    case default
       call msg_fatal ("pdf_builtin: internal: invalid PDF set!")
    end select
  end subroutine pdf_evolve

! included for compatibility with LHAPDF
! use a double precision array for the pdfs instead of a
! default precision
  module subroutine pdf_evolve_LHAPDF (set, x, q, ff)
    integer, intent(in) :: set
    double precision, intent(in) :: x, q
    real(kind=default) :: dx, dq
    double precision, dimension(-6:6), intent(out) :: ff

    real(kind=default) :: f(-6:6)
    dx = x
    dq = q

    call pdf_evolve (set, dx, dq, f)
    ff = f
  end subroutine pdf_evolve_LHAPDF

! PDF-specific running alphas
  module function pdf_alphas (pdftype, q) result (alphas)
    real(kind=double) :: as
    real(kind=default), intent(in) :: q
    real(kind=double) :: qdummy
    real(kind=default) :: alphas
    integer, intent(in) :: pdftype
    integer, parameter :: nset_dummy = 1
    qdummy = q
    select case (pdftype)
    case (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1)
       if (cteq6_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (cteq6_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (cteq6_initialized)) // &
            " cannot be used simultaneously")
       select case (pdftype)
          case (CTEQ6M, CTEQ6D)
             as = PI*pdf_builtin_alpi (qdummy,2,5,alam_cteq6m)
          case (CTEQ6L)
             as = PI*pdf_builtin_alpi (qdummy,2,5,alam_cteq6l)
          case (CTEQ6L1)
             as = PI*pdf_builtin_alpi (qdummy,1,5,alam_cteq6ll)
          case default
             call msg_fatal ("Unrecognized PDF setup in evolution of alphas.")
       end select
    case (CT10)
       if (.not. ct10_initialized) call msg_fatal ( &
            "pdf_builtin: internal: PDF set CT10 requested without " // &
            "initialization!")
        as = PI*pdf_builtin_alpi (qdummy,2,5,alam_ct10)
    case (MRST2004QEDp)
       if (.not. mrst2004qedp_initialized) call msg_fatal ( &
            "pdf_builtin: internal: PDF set MRST2004QEDp requested without " // &
            "initialization!")
       call mrst_alphas (as,qdummy)
    case (MRST2004QEDn)
       if (.not. mrst2004qedn_initialized) call msg_fatal ( &
            "pdf_builtin: internal: PDF set MRST2004QEDn requested without " // &
            "initialization!")
       call mrst_alphas (as,qdummy)
    case (MSTW2008LO, MSTW2008NLO, MSTW2008NNLO)
       if (mstw2008_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (mstw2008_initialized /= pdftype) &
            call msg_fatal ( &
              "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
              char (pdf_get_name (mstw2008_initialized)) // &
              " cannot be used simultaneously")
       select case (pdftype)
          case (MSTW2008lo)
             call pdf_builtin_alfa (as,qdummy,0)
          case (MSTW2008nlo)
             call pdf_builtin_alfa (as,qdummy,1)
          case (MSTW2008nnlo)
             call pdf_builtin_alfa (as,qdummy,2)
       end select
    case (CJ12_max, CJ12_mid, CJ12_min)
       if (cj12_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (cj12_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (cj12_initialized)) // &
            " cannot be used simultaneously")
       as = PI*pdf_builtin_alpi (qdummy,2,5,alam_cteq6m)
    case (CJ15_LO, CJ15_NLO)
       if (cj15_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (cj15_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (cj15_initialized)) // &
            " cannot be used simultaneously")
       select case (pdftype)
       case (CJ15_NLO)
          as = PI*pdf_builtin_alpi (qdummy,2,5,alam_cteq6m)
       case (CJ15_LO)
          as = PI*pdf_builtin_alpi (qdummy,2,5,alam_cteq6l)
       end select
    case (MMHT2014LO, MMHT2014NLO, MMHT2014NNLO)
       if (mmht2014_initialized < 0) &
            call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (mmht2014_initialized /= pdftype) &
            call msg_fatal ( &
              "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
              char (pdf_get_name (mmht2014_initialized)) // &
              " cannot be used simultaneously")
       select case (pdftype)
          case (MMHT2014lo)
             call pdf_builtin_alfa (as,qdummy,0)
          case (MMHT2014nlo)
             call pdf_builtin_alfa (as,qdummy,1)
          case (MMHT2014nnlo)
             call pdf_builtin_alfa (as,qdummy,2)
       end select
    case (CT14LL, CT14L, CT14N, CT14NN)
       call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (ct14_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (ct14_initialized)) // &
            " cannot be used simultaneously")
       as = CT14Alphas(qdummy)
    case (CT18N, CT18NN)
       call msg_fatal ("pdf_builtin: internal: PDF set " // &
            char (pdf_get_name (pdftype)) // " requested without initialization!")
       if (ct18_initialized /= pdftype) &
            call msg_fatal ( &
            "PDF sets " // char (pdf_get_name (pdftype)) // " and " // &
            char (pdf_get_name (ct18_initialized)) // &
            " cannot be used simultaneously")
       as = CT18Alphas(qdummy)
    case default
         call msg_fatal ("pdf_builtin: internal: invalid PDF set!")
    end select
    alphas = as
  end function pdf_alphas

  module function pdf_alphas_LHAPDF (pdftype, q) result (alphas)
    integer, intent(in) :: pdftype
    real(kind=double), intent(in) :: q
    real(kind=double) :: alphas
    real(kind=default) :: q_def
    q_def = q
    alphas = pdf_alphas (pdftype, q_def)
  end function pdf_alphas_LHAPDF

  module function pdf_getmass (nf) result (mass)
    integer, intent(in) :: nf
    real(kind=double) :: mass
    select case (abs (nf))
       case (1:3)
          mass = 0._double
       case (4)
          mass = 1.3_double
       case (5)
          mass = 4.5_double
       case default
          call msg_fatal ("PDF builtin: invalid PDG code for quark mass.")
       end select
  end function pdf_getmass

  function pdf_builtin_alpi (AMU,NORDER,NF,alam)
    integer, intent(in) :: nf, norder
    real(kind=double), intent(in) :: amu
    real(kind=double), dimension(6), intent(in) :: alam
    real(kind=double) :: alm
    real(kind=double) :: pdf_builtin_alpi
    integer :: irt, n_eff
    n_eff = num_flavor (AMU,NF)
    alm  = alam (n_eff+1)
    pdf_builtin_alpi = pdf_builtin_alpqcd (norder, n_eff, amu/alm, irt)
    select case (irt)
       case (1)
          call msg_warning ("AMU < ALAM in pdf_builtin_alpi")
       case (2)
          call msg_warning ("pdf_builtin_alpi > 3; Be aware!")
    end select
  end function pdf_builtin_alpi

  double precision function pdf_builtin_ALPQCD (IRDR, NF, RML, IRT)
    real(kind=double) :: b0, b1
    integer, intent(in) :: irdr, nf
    integer :: irt
    real(kind=double), intent(in) :: rml
    real(kind=double) :: al, aln
    real(kind=double) :: rm2
    IRT = 0
    if (IRDR .LT. 1 .OR. IRDR .GT. 2) THEN
       print *,                                                        &
            &  'Order out of range in pdf_builtin_ALPQCD: IRDR = ', IRDR
       STOP
    end if
    b0 = (11.d0*CA  - 2.* NF) / 3.d0
    b1 = (34.d0*CA**2 - 10.d0*CA*NF - 6.d0*CF*NF) / 3.d0
    rm2 = rml**2
    if (RM2 .LE. 1.) THEN
       IRT = 1
       pdf_builtin_ALPQCD = 99.
       RETURN
    end if
    aln = log (RM2)
    al = 4.d0/ B0 / aln
    IF (IRDR .GE. 2) AL = AL * (1.d0-B1*log(ALN) / ALN / B0**2)
    if (AL .GE. 3.) THEN
       IRT = 2
    end if
    pdf_builtin_ALPQCD = AL
  end function pdf_builtin_ALPQCD

  function num_flavor (amu,nf)
    real(kind=double), intent(in) :: amu
    integer, intent(in) :: nf
    integer :: num_flavor, i
    num_flavor = nf - nhq
    IF ((num_flavor .EQ. nf) .OR. (amu .LE. amn)) GOTO 20
    do 10 I = nf - NHQ + 1, nf
       if (amu .GE. amhat(I)) THEN
          num_flavor = I
       else
          goto 20
       end if
10     continue
20     return
  end function num_flavor

  subroutine pdf_builtin_alfa (alfas,Qalfa,norder)
    integer, intent(in) :: norder
    real(kind=double), intent(in) :: Qalfa
    real(kind=double), intent(out) :: alfas
    real(kind=double) :: alphasq0
    real(kind=double) :: mCharmSave, mBottomSave, mTopSave
    mTopSave = 10000000000.0_double
    mBottomSave = 4.75_double
    mCharmSave = 1.4_double
    select case (norder)
       case (0)
          alphasq0 = 0.68183_double
       case (1)
          alphasq0 = 0.491278_double
       case (2)
          alphasq0 = 0.45077_double
       case default
          call msg_fatal ("Wrong MSTW order!")
    end select
    if (Qalfa.GT.mTopSave) call msg_fatal &
         ("No MSTW alphas for such high scales.")
    alfas = pdf_builtin_mstw_alphas (Qalfa,norder,mTopSave**2)
  end subroutine pdf_builtin_alfa

  function pdf_builtin_mstw_alphas (mur, norder, m2t) result (alfas)
    real(kind=double), intent(in) :: mur
    real(kind=double) :: alfas
    integer, intent(in) :: norder
    integer :: nf
    real(kind=double), intent(in) :: m2t
    real(kind=double) :: as0, asc, asb
    real(kind=double) :: r2, m2, asi, asf, r20, r2b, r2c, r2t
    real(kind=double), parameter :: m20 = 1.0_double, m2c = 1.96_double, &
         m2b = 22.5625_double
    real(kind=double) :: logfr
    integer :: nff
    real(kind=double) :: ast
    nff = 4
    select case (norder)
    case (0)
       as0 = 5.4258307424173556E-002_double
       asc = 4.0838232991033577E-002_double
       asb = 2.2297506503539639E-002_double
    case (1)
       as0 = 3.9094820221093209E-002_double
       asc = 3.0202751103489092E-002_double
       asb = 1.7751676069301441E-002_double
    case (2)
       as0 = 3.5871136848766867E-002_double
       asc = 2.8092349518054758E-002_double
       asb = 1.6936586710596203E-002_double
    case default
       call msg_fatal ("Wrong evolution order in MSTW alphas evolution.")
    end select
    r2 = mur**2
    m2 = r2 * exp(+logfr)
    if (m2 .gt. m2t) then
       nf = 6
       r2t = m2t * r2/m2
       asi = ast
       asf = mstw_alphas (r2,r2t,ast,nf,norder)
    else
       if (m2 .gt. m2b) then
          nf = 5
          r2b = m2b * r2/m2
          asi = asb
          asf = mstw_alphas (r2,r2b,asb,nf,norder)
       else
          if (m2 .gt. m2c) then
             nf = 4
             r2c = m2c * r2/m2
             asi = asc
             asf = mstw_alphas (r2,r2c,asc,nf,norder)
          else
             nf = 3
             r20 = m20 * r2/m2
             asi = as0
             asf = mstw_alphas (r2,r20,as0,nf,norder)
          end if
       end if
    end if
    alfas = 4.D0*PI*ASF
  end function pdf_builtin_mstw_alphas

  function mstw_alphas (r2, r20, as0, nf, naord) result (as)
    real(kind=double) :: as
    real(kind=double), intent(in) :: r2, r20, as0
    integer, intent(in) :: nf, naord
    integer, parameter :: nastps = 20
    integer :: k1
    real(kind=double), parameter :: sxth = 1./6.
    real(kind=double) :: lrrat, dlr, xk0, xk1, xk2, xk3
    as = as0
    lrrat = log (r2/r20)
    dlr = lrrat / nastps
    ! ..Solution of the evolution equation depending on  NAORD
    !   (fourth-order Runge-Kutta beyond the leading order)
    select case (naord)
    case (0)
       as = as0 / (1.+ beta0(nf) * as0 * lrrat)
    case (1)
       do k1 = 1, nastps
          xk0 = dlr * fbeta1 (as, nf)
          xk1 = dlr * fbeta1 (as + 0.5 * xk0, nf)
          xk2 = dlr * fbeta1 (as + 0.5 * xk1, nf)
          xk3 = dlr * fbeta1 (as + xk2, nf)
          as = as + sxth * (xk0 + 2.* xk1 + 2.* xk2 + xk3)
       end do
    case (2)
       do k1 = 1, nastps
          xk0 = dlr * fbeta2 (as, nf)
          xk1 = dlr * fbeta2 (as + 0.5 * xk0, nf)
          xk2 = dlr * fbeta2 (as + 0.5 * xk1, nf)
          xk3 = dlr * fbeta2 (as + xk2, nf)
          as = as + sxth * (xk0 + 2.* xk1 + 2.* xk2 + xk3)
       end do
    case (3)
       do k1 = 1, nastps
          xk0 = dlr * fbeta3 (as, nf)
          xk1 = dlr * fbeta3 (as + 0.5 * xk0, nf)
          xk2 = dlr * fbeta3 (as + 0.5 * xk1, nf)
          xk3 = dlr * fbeta3 (as + xk2, nf)
          as = as + sxth * (xk0 + 2.* xk1 + 2.* xk2 + xk3)
       end do
    end select
  end function mstw_alphas

  double precision function beta0 (nf)
    integer, intent(in) :: nf
    beta0 = b00 + b01 * nf
  end function beta0

  double precision function beta1 (nf)
    integer, intent(in) :: nf
    beta1 = b10 + b11 * nf
  end function beta1

  double precision function beta2 (nf)
    integer, intent(in) :: nf
    beta2 = 1428.50 - 279.611 * nf + 6.01852 * nf**2
  end function beta2

  double precision function beta3 (nf)
    integer, intent(in) :: nf
    beta3 = 29243.0 - 6946.30 * nf + 405.089 * nf**2 + &
         1.49931 * nf**3
  end function beta3

  ! ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
  double precision function fbeta1 (a, nf)
    real(kind=double), intent(in) :: a
    integer, intent(in) :: nf
    fbeta1 = - a**2 * ( beta0(nf) + a * beta1(nf) )
  end function fbeta1

  double precision function fbeta2 (a, nf)
    real(kind=double), intent(in) :: a
    integer, intent(in) :: nf
    fbeta2 = - a**2 * ( beta0(nf) + a * ( beta1(nf) &
         + a * beta2(nf) ) )
  end function fbeta2

  double precision function fbeta3 (a, nf)
    real(kind=double), intent(in) :: a
    integer, intent(in) :: nf
    fbeta3 = - a**2 * ( beta0(nf) + a * ( beta1(nf) &
         + a * ( beta2(nf) + a * beta3(nf)) ) )
  end function fbeta3

  subroutine mrst_alphas (alpha, q)
    real(kind=double), intent(in) :: q
    real(kind=double), intent(out) :: alpha
    integer :: idir
    integer, parameter :: nflav = 4
    real(kind=double) :: alphas, qz2, qz, alambda, astep, &
         tol2
    alphas = 0.1205_double
    qz2=8315._double
    qz = sqrt(qz2)
    alambda = 0.3000_double
    astep = 0.010_double
    tol2 = 0.0000001_double
    idir = +1
100  continue
    alambda = alambda + idir*astep
    call mrst_lambda (nflav,alpha,qz,alambda)
    if (idir*(alphas-alpha).gt.0.0) then
       goto 200
    else
       astep = 0.5*astep
       idir = -1*idir
       goto 200
    end if
200  continue
    if(abs(alpha-alphas).gt.tol2) goto 100
    ! alambda found  -- save it !!!!
    ! next call mrst_lambda to get alphas at q with the correct alambda
    call mrst_lambda (nflav,alpha,q,alambda)
    RETURN
  end subroutine mrst_alphas

  subroutine mrst_lambda (nflav,alpha,Q,alambda)
    integer, intent(in) :: nflav
    integer :: ith
    real(kind=double), intent(in) :: q, alambda
    real(kind=double), intent(out) :: alpha
    real(kind=double) :: al, al2, alfinv, alfqc3, alfqc4, alfqc5, &
         alfqs3, alfqs5
    real(kind=double) :: x1, x2, t, tt
    real(kind=double) :: qsdt, qsdtt, qsct, qsctt, qs, q2
    real(kind=double) :: b0, b1, del, f, fp, as, as2
    real(kind=double) :: flav
    !  The value of Lambda required corresponds to nflav=4
    qsdt=8.18    !!  This is the value of 4m_c^2
    qsct=74.0    !!  This is the value of 4m_b^2
    al2=alambda*alambda
    q2=q*q
    t = log (q2/al2)
    ! CHECK: explicitly initialising ALFQC{3,4,5} (by AB)
    alfqc3 = 0
    alfqc4 = 0
    alfqc5 = 0
    ith=0
    tt=t
    qsdtt=qsdt/4.
    qsctt=qsct/4.
    al = alambda
    al2=al*al
    flav=4.
    qs=al2*exp(T)
    if (qs.lt.0.5d0) then   !!  running stops below 0.5
       qs=0.5d0
       t = log (qs/al2)
       tt= t
    end if
    if (qs.lt.qsdtt .and. nflav.gt.3) GOTO 312
11  continue
    b0 = 11._double - 2./3._double * flav
    x1 = 4.*PI/b0
    b1=102.-38.*FLAV/3.
    X2=B1/B0**2
    AS2=X1/T*(1.-X2* log(T)/T)
5   AS=AS2
    F=-T+X1/AS-X2*log(X1/AS+X2)
    FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
    AS2=AS-F/FP
    DEL=ABS(F/FP/AS)
    IF((DEL-tol).GT.0.) goto 5
    ALPHA=AS2
51  continue
    IF(ITH.EQ.0) RETURN
    GOTO (13,14,15) ITH
    ! GOTO 5
12  ITH=1
    T=log(QSCTT/AL2)
    GOTO 11
13  ALFQC4=ALPHA
    FLAV=5.
    ITH=2
    GOTO 11
14  ALFQC5=ALPHA
    ITH=3
    T=TT
    GOTO 11
15  ALFQS5=ALPHA
    ALFINV=1./ALFQS5+1./ALFQC4-1./ALFQC5
    ALPHA=1./ALFINV
    RETURN
311 CONTINUE
    B0=11-2.*FLAV/3.
    X1=4.*PI/B0
    B1=102.-38.*FLAV/3.
    X2=B1/B0**2
    AS2=X1/T*(1.-X2*log(T)/T)
35  AS=AS2
    F=-T+X1/AS-X2*log(X1/AS+X2)
    FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
    AS2=AS-F/FP
    DEL=ABS(F/FP/AS)
    IF((DEL-tol).GT.0.) goto 35
    ALPHA=AS2
351 continue
    IF(ITH.EQ.0) RETURN
    GOTO (313,314,315) ITH
312 ITH=1
    T=log(QSDTT/AL2)
    GOTO 311
313 ALFQC4=ALPHA
    FLAV=3.
    ITH=2
    GOTO 311
314 ALFQC3=ALPHA
    ITH=3
    T=TT
    GOTO 311
315 ALFQS3=ALPHA
    ALFINV=1./ALFQS3+1./ALFQC4-1./ALFQC3
    ALPHA=1./ALFINV
  end subroutine mrst_lambda

end submodule pdf_builtin_s
