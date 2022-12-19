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

submodule (electron_pdfs) electron_pdfs_s

  use iso_varying_string, string_t => varying_string
  use io_units
  use sm_qed
  use diagnostics
  use constants
  use format_defs, only: FMT_12
  use numeric_utils
  use physics_defs
  use sm_physics

  implicit none

contains

  module subroutine qed_pdf_init &
       (qed_pdf, mass, alpha, charge, q_max, order, log_order, n_lep)
    class(qed_pdf_t), intent(out) :: qed_pdf
    real(default), intent(in) :: mass, alpha, q_max, charge
    integer, intent(in) :: order, log_order, n_lep
    qed_pdf%mass = mass
    qed_pdf%q_max = q_max
    qed_pdf%alpha = alpha
    qed_pdf%order = order
    qed_pdf%log_order = log_order
    qed_pdf%n_lep = n_lep
    qed_pdf%eps = alpha/pi * charge**2 &
         * (2 * log (q_max / mass) - 1)
  end subroutine qed_pdf_init

  module subroutine qed_pdf_write (qed_pdf, unit, with_qed)
    class(qed_pdf_t), intent(in) :: qed_pdf
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: with_qed
    integer :: u
    logical :: show_qed
    u = given_output_unit (unit)
    show_qed = .false.
    if (present (with_qed))  show_qed = with_qed
    write (u, "(3x,A)")  "QED structure function (PDF):"
    write (u, "(5x,A,I0)") "Flavor    = ", qed_pdf%flv
    write (u, "(5x,A," // FMT_12 // ")") "Mass      = ", qed_pdf%mass
    write (u, "(5x,A," // FMT_12 // ")") "q_max     = ", qed_pdf%q_max
    write (u, "(5x,A," // FMT_12 // ")") "alpha     = ", qed_pdf%alpha
    write (u, "(5x,A,I0)") "Order     = ", qed_pdf%order
    write (u, "(5x,A,I0)") "Log. ord. = ", qed_pdf%log_order
    write (u, "(5x,A,I0)") "# leptons = ", qed_pdf%n_lep
    write (u, "(5x,A," // FMT_12 // ")") "epsilon   = ", qed_pdf%eps
    if (show_qed) then
       call qed_pdf%aqed%write (u)
    end if
  end subroutine qed_pdf_write

  module subroutine qed_pdf_set_order (qed_pdf, order)
    class(qed_pdf_t), intent(inout) :: qed_pdf
    integer, intent(in) :: order
    qed_pdf%order = order
  end subroutine qed_pdf_set_order

  module subroutine qed_pdf_evolve_qed_pdf (qed_pdf, x, xb, rb, ff)
    class(qed_pdf_t), intent(inout) :: qed_pdf
    real(default), intent(in) :: x, xb, rb
    real(default), intent(inout) :: ff
    real(default), parameter :: &
         & xmin = 0.00714053329734592839549879772019_default   
    real(default), parameter :: &
         g1 = 3._default / 4._default, &
         g2 = (27 - 8 * pi**2) / 96._default, &
         g3 = (27 - 24 * pi**2 + 128 * zeta3) / 384._default
    real(default) :: x_2, log_x, log_xb
    if (ff > 0 .and. qed_pdf%order > 0) then
       ff = ff * (1 + g1 * qed_pdf%eps)
       x_2 = x * x
       if (rb > 0)  ff = ff * (1 - (1-x_2) / (2 * rb))
       if (qed_pdf%order > 1) then
          ff = ff * (1 + g2 * qed_pdf%eps**2)
          if (rb > 0 .and. xb > 0 .and. x > xmin) then
             log_x  = log_prec (x, xb)
             log_xb = log_prec (xb, x)
             ff = ff * (1 - ((1 + 3 * x_2) * log_x + xb * (4 * (1 + x) * &
                  log_xb + 5 + x)) / (8 * rb) * qed_pdf%eps)
          end if
          if (qed_pdf%order > 2) then
             ff = ff * (1 + g3 * qed_pdf%eps**3)
             if (rb > 0 .and. xb > 0 .and. x > xmin) then
                ff = ff * (1 - ((1 + x) * xb &
                     * (6 * Li2(x) + 12 * log_xb**2 - 3 * pi**2) &
                     + 1.5_default * (1 + 8 * x + 3 * x_2) * log_x &
                     + 6 * (x + 5) * xb * log_xb &
                     + 12 * (1 + x_2) * log_x * log_xb &
                     - (1 + 7 * x_2) * log_x**2 / 2 &
                     + (39 - 24 * x - 15 * x_2) / 4) &
                     / (48 * rb) * qed_pdf%eps**2)
             end if
          end if
       end if
    end if
  end subroutine qed_pdf_evolve_qed_pdf

  module function elec_asym (epdf, x, scale, alpha, running) result (elec_as)
    type(qed_pdf_t), intent(in) :: epdf
    real(default) :: elec_as
    real(default), intent(in) :: x
    real(default), intent(in) :: scale
    real(default), intent(in) :: alpha
    logical, intent(in) :: running
    real(default) :: lambda0, lambda1, xi1, xihat1, fac
    real(default) :: al0_2pi, al2pi, ca, cb, psixi1
    real(default) :: ln0, eta0, t, xi0, b0, b1
    integer :: nf, nlep
    lambda0 = 3._default/4._default
    al0_2pi = alpha / two / Pi
    if (allocated (epdf%q_in)) then
       ln0 = log(epdf%q_in**2/epdf%mass**2)
       eta0 = alpha/Pi * log(scale**2/epdf%q_in**2)
    else
       ln0 = zero
       eta0 = alpha/Pi * log(scale**2/epdf%mass**2)
    end if
    if (running) then
       t = t_alpha (epdf, scale)
       xi0 = two*t
    else
       xi0 = eta0
    end if
    select case (epdf%log_order)
    case (EPDF_LL)
       elec_as = exp((lambda0 - eulerc) * xi0) / gamma (one + xi0) * &
            xi0 * (1 - x)**(-one + xi0)
    case (EPDF_NLL)
       lambda1 = 3._default/8._default - Pi2/two + 6._default*zeta3 - &
            epdf%n_lep/18._default * (three + four*Pi2)
       if (running) then
          select type (aqed => epdf%aqed)
          type is (alpha_qed_from_scale_t)
             nf = aqed%nf
             nlep = aqed%nlep
             al2pi = aqed%get (scale) / two / Pi
          type is (alpha_qed_fixed_t)
             call msg_fatal &
                  ("elec_asym: has to be called with running alpha.")
          end select
          b0 = coeffqed_b0(nf,nlep)
          b1 = coeffqed_b1(nf,nlep)
          xi1 = xi0 + al2pi / two / Pi / b0 * (one - exp (Pi * b0 * xi0)) * &
               (20._default/9._default * epdf%n_lep + four * Pi * b1 / b0)
          xihat1 = xi0 * lambda0 - al2pi / two / Pi / b0 * &
               (one - exp (Pi * b0 * xi0)) * (lambda1 - three * Pi * b1 / b0)
       else
          al2pi = al0_2pi
          xi1 = xi0 * ( one - 10._default/9._default * al2pi * epdf%n_lep)
          xihat1 = xi0 * ( lambda0 + al2pi/two * lambda1 )
       end if
       fac = exp( xihat1 - eulerc * xi1 ) / gamma (1 + xi1) * &
            xi1 * (1 - x)**(-one + xi1)
       psixi1 = psir(xi1)
       ca = - eulerc - psixi1
       cb = 0.5_default * eulerc**2 + Pi2/12. + eulerc*psixi1 + &
            0.5_default * psixi1**2 - 0.5_default * psimr(xi1,1)
       elec_as = fac * (one + two*al0_2pi * &
            ( (ln0 - one)*( ca + 3._default/4._default ) - &
            2*cb + 7._default/4._default + &
            (ln0 - one - two*ca)*log(1-x) - log(1-x)**2 ))
    case default
       elec_as = 0
    end select
  end function elec_asym

  module function phot_asym &
       (epdf, x, scale, alpha, nlep, running) result (phot_as)
    type(qed_pdf_t), intent(in) :: epdf
    real(default) :: phot_as
    real(default), intent(in) :: x
    real(default), intent(in) :: scale
    real(default), intent(in) :: alpha
    integer, intent(in) :: nlep
    logical, intent(in) :: running
    real(default) :: ln0, eta0, t, xi0, xihat0, den
    real(default) :: al0_2pi, c_b0, c_b1ob0
    real(default) :: lambda1, xi10, xihat10
    real(default) :: chi10, mf1k, mf10
    real(default) :: gam1, gam2, gam3, gam4, gam5
    real(default) :: d11, d21, d12, d22, d13, d23, d14, d24
    real(default) :: c11, c21, c31, c12, c22, c32, c13, c23, c33, &
         c14, c24, c34, k1, k2, k3, k4
    al0_2pi = alpha / two / Pi
    if (allocated (epdf%q_in)) then
       ln0 = log(epdf%q_in**2/epdf%mass**2)
       eta0 = alpha/Pi * log(scale**2/epdf%q_in**2)
    else
       ln0 = zero
       eta0 = alpha/Pi * log(scale**2/epdf%mass**2)
    end if
    if (running) then
       c_b0 = - coeffqed_b0 (0, nlep)
       c_b1ob0 = coeffqed_b1 (0, nlep) / coeffqed_b0 (0, nlep)
       t = t_alpha (epdf, scale)
       chi10 = zero
    else
       c_b0 = zero
       c_b1ob0 = zero
       t = eta0 / two
       chi10 = al0_2pi * nlep
    end if
    select case (epdf%log_order)
    case (EPDF_LL)
       xi0 = two * t
       xihat0 = three * t / two
       den = -two/three * nlep - three/two - two * log(1-x)
       mf1k = exp(-eulerc * xi0) * (1-x)**xi0 / gamma(1+xi0) * &
            (one / den - (Pi2*xi0 - 6._default * zeta3 * xi0**2) / &
            three / den**2 - (30._default * Pi2 - 360._default * &
            zeta3 * xi0 + Pi**4 * xi0**2) / 45._default / den**3)
       mf10 = one / den - 2._default * Pi2 / three / den**3
       phot_as = - exp(xihat0) * mf1k + &
            exp(-two/three * nlep * t) * mf10
    case (EPDF_NLL)
       lambda1 = three/8._default - Pi2/two + 6._default*zeta3 - &
            nlep/18._default * (three + four*Pi2)
       xi10 = two - four*al0_2pi * (5._default/9._default * nlep + Pi*c_b1ob0)
       xihat10 = three/two * ( one + al0_2pi * two * &
            (lambda1/three - Pi * c_b1ob0))
       d11 = xi10
       d21 = - (two/three * nlep + two * Pi * c_b0 + xihat10 + chi10)
       c11 = al0_2pi * exp(-D21*t)
       c21 = - al0_2pi * (five + four/three * nlep) * exp(-d21*t)
       c31 = al0_2pi * (6._default + Pi2/6._default + 32._default/9._default * &
            nlep + two * Pi * c_b1ob0) * exp(-d21*t)
       d12 = d11
       d22 = d21
       c12 = - al0_2pi
       c22 = al0_2pi * (5._default + four/three * nlep)
       c32 = - al0_2pi * (6._default + Pi2/6._default + &
            32._default/9._default * nlep + two *Pi *c_b1ob0)
       d13 = d11
       d23 = - (two/three * nlep + xihat10 + chi10)
       c13 = zero
       c23 = zero
       c33 = - exp(-d23*t)
       d14 = d13
       d24 = d23
       c14 = zero
       c24 = zero
       c34 = one
       k1 = xi10 * t
       k2 = zero
       k3 = xi10 * t
       k4 = zero
       gam1 = sum_rm (x, al0_2pi, ln0, &
            c11, c21, c31, d21/d11, d11, k1, d11, d21)
       gam2 = sum_rm (x, al0_2pi, ln0, &
            c12, c22, c32, d22/d12, d12, k2, d12, d22)
       gam3 = sum_rm (x, al0_2pi, ln0, &
            c13, c23, c33, d23/d13, d13, k3, d13, d23)
       gam4 = sum_rm (x, al0_2pi, ln0, &
            c14, c24, c34, d24/d14, d14, k4, d14, d24)
       gam5 = al0_2pi * (one + (1-x)**2)/x * (ln0 - two*log(x) - one)
       phot_as = exp( - ( two/three * nlep + chi10 ) * t ) * &
            (gam1 + gam2 + gam3 + gam4 + gam5)
    end select
  end function phot_asym

  module function bar_asym &
       (epdf, flv, x, scale, alpha, running) result (bar_as)
    type(qed_pdf_t), intent(in) :: epdf
    integer, intent(in) :: flv
    real(default) :: bar_as
    real(default), intent(in) :: x
    real(default), intent(in) :: scale
    real(default), intent(in) :: alpha
    logical, intent(in) :: running
    logical, dimension(6) :: order
    integer :: nlep
    real(default) :: ln0, eta0, al_2pi, p
    real(default), dimension(6) :: jll_nll
    if (allocated (epdf%q_in)) then
       ln0 = log(epdf%q_in**2/epdf%mass**2)
       eta0 = alpha/Pi * log(scale**2/epdf%q_in**2) 
    else
       ln0 = zero
       eta0 = alpha/Pi * log(scale**2/epdf%mass**2)
    end if
    if (running) then
       p = t_alpha (epdf, scale)
    else
       p = eta0 / two
    end if
    order = .false.
    order(1:3) = .true.
    select case (epdf%log_order)
    case (EPDF_LL)
       nlep = epdf%n_lep
       al_2pi = zero
    case (EPDF_NLL)
       if (running) then
          select type (aqed => epdf%aqed)
          type is (alpha_qed_from_scale_t)
             nlep = aqed%nlep
             al_2pi = aqed%get (scale) / two / Pi
          type is (alpha_qed_fixed_t)
             call msg_fatal &
                  ("bar_asym: has to be called with running alpha.")
          end select
       else
          nlep = epdf%n_lep
          al_2pi = alpha / two / Pi
       end if
       order(4:6) = .true.
    end select
    select case (flv)
    case (EPDF_ELE,EPDF_POS,EPDF_S,EPDF_NS)
       call elecbar_asym_p (x, jll_nll, nlep, ln0, order, running)
    case (EPDF_G)
       call photbar_asym_p (x, jll_nll, nlep, ln0, order, running)
    case default
       call msg_fatal &
            ("bar_asym: wrong lepton flavor.")
    end select
    bar_as = rec_series (p, al_2pi, jll_nll)
  end function bar_asym

  module function rec_num &
       (epdf, flv, x, scale, alpha, running) result (recnum)
    type(qed_pdf_t), intent(in) :: epdf
    integer, intent(in) :: flv
    real(default) :: recnum
    real(default), intent(in) :: x
    real(default), intent(in) :: scale
    real(default), intent(in) :: alpha
    logical, intent(in) :: running
    integer :: nlep
    real(default) :: ln0, eta0, al_2pi, p, prefac
    if (allocated (epdf%q_in)) then
       ln0 = log(epdf%q_in**2/epdf%mass**2)
       eta0 = alpha/Pi * log(scale**2/epdf%q_in**2) 
    else
       ln0 = zero
       eta0 = alpha/Pi * log(scale**2/epdf%mass**2)
    end if
    if (running) then
       p = t_alpha (epdf, scale)
    else
       p = eta0 / two
    end if
    select case (epdf%log_order)
    case (EPDF_LL)
       recnum = 0
       return
    case (EPDF_NLL)
       if (running) then
          select type (aqed => epdf%aqed)
          type is (alpha_qed_from_scale_t)
             nlep = aqed%nlep
             al_2pi = aqed%get (scale) / two / Pi
          type is (alpha_qed_fixed_t)
             call msg_fatal &
                  ("bar_asym: has to be called with running alpha.")
          end select
       else
          nlep = epdf%n_lep
          al_2pi = alpha / two / Pi
       end if
       prefac = al_2pi * p**2 / two
    end select
    select case (flv)
    case (EPDF_S)
       recnum = prefac * endpoint_func_S (x, nlep)
    case (EPDF_NS)
       recnum = prefac * endpoint_func_NS (x)
    case (EPDF_G)
       recnum = prefac * endpoint_func_GAM (x)
    case default
       call msg_fatal &
            ("rec_num: wrong lepton flavor.")
    end select
  end function rec_num

  module function recbar &
       (epdf, flv, x, scale, alpha, running) result (bar)
    type(qed_pdf_t), intent(in) :: epdf
    integer, intent(in) :: flv
    real(default) :: bar
    real(default), intent(in) :: x
    real(default), intent(in) :: scale
    real(default), intent(in) :: alpha
    logical, intent(in) :: running
    logical, dimension(6) :: order
    integer :: nlep
    real(default) :: ln0, eta0, al_2pi, p
    real(default), dimension(6) :: jll_nll
    if (allocated (epdf%q_in)) then
       ln0 = log(epdf%q_in**2/epdf%mass**2)
       eta0 = alpha/Pi * log(scale**2/epdf%q_in**2) 
    else
       ln0 = zero
       eta0 = alpha/Pi * log(scale**2/epdf%mass**2)
    end if
    if (running) then
       p = t_alpha (epdf, scale)
    else
       p = eta0 / two
    end if
    order = .false.
    order(1:3) = .true.
    select case (epdf%log_order)
    case (EPDF_LL)
       nlep = epdf%n_lep
       al_2pi = zero
    case (EPDF_NLL)
       if (running) then
          select type (aqed => epdf%aqed)
          type is (alpha_qed_from_scale_t)
             nlep = aqed%nlep
             al_2pi = aqed%get (scale) / two / Pi
          type is (alpha_qed_fixed_t)
             call msg_fatal &
                  ("recbar: has to be called with running alpha.")
          end select
       else
          nlep = epdf%n_lep
          al_2pi = alpha / two / Pi
       end if
       order(4:6) = .true.
    end select
    select case (flv)
    case (EPDF_S)
       call recbar_singlet (x, jll_nll, nlep, ln0, order, running)
    case (EPDF_G)
       call recbar_photon (x, jll_nll, nlep, ln0, order, running)
    case (EPDF_NS)
       call recbar_nonsinglet (x, jll_nll, nlep, ln0, order, running)
    case default
       call msg_fatal &
            ("recbar: wrong lepton flavor.")
    end select
    bar = rec_series (p, al_2pi, jll_nll)
  end function recbar

  module function rechat &
       (epdf, flv, x, scale, alpha, running) result (hat)
    type(qed_pdf_t), intent(in) :: epdf
    integer, intent(in) :: flv
    real(default) :: hat
    real(default), intent(in) :: x
    real(default), intent(in) :: scale
    real(default), intent(in) :: alpha
    logical, intent(in) :: running
    logical, dimension(6) :: order
    integer :: nlep
    real(default) :: ln0, eta0, al_2pi, p
    real(default), dimension(6) :: jll_nll
    if (allocated (epdf%q_in)) then
       ln0 = log(epdf%q_in**2/epdf%mass**2)
       eta0 = alpha/Pi * log(scale**2/epdf%q_in**2) 
    else
       ln0 = zero
       eta0 = alpha/Pi * log(scale**2/epdf%mass**2)
    end if
    if (running) then
       p = t_alpha (epdf, scale)
    else
       p = eta0 / two
    end if
    order = .false.
    order(1:3) = .true.
    select case (epdf%log_order)
    case (EPDF_LL)
       nlep = epdf%n_lep
       al_2pi = zero
    case (EPDF_NLL)
       if (running) then
          select type (aqed => epdf%aqed)
          type is (alpha_qed_from_scale_t)
             nlep = aqed%nlep
             al_2pi = aqed%get (scale) / two / Pi
          type is (alpha_qed_fixed_t)
             call msg_fatal &
                  ("recbar: has to be called with running alpha.")
          end select
       else
          nlep = epdf%n_lep
          al_2pi = alpha / two / Pi
       end if
       order(4:6) = .true.
    end select
    select case (flv)
    case (EPDF_S)
       call rechat_singlet (x, jll_nll, nlep, ln0, order, running)
    case (EPDF_G)
       call rechat_photon (x, jll_nll, nlep, ln0, order, running)
    case (EPDF_NS)
       call rechat_nonsinglet (x, jll_nll, nlep, ln0, order, running)
    case default
       call msg_fatal &
            ("rechat: wrong lepton flavor.")
    end select
    hat = rec_series (p, al_2pi, jll_nll)
  end function rechat

  module subroutine elecbar_asym_p (x, jll_nll, nlep, ln0, order, running)
    real(default), intent(in) :: x
    real(default), dimension(6), intent(out) :: jll_nll
    real(default), intent(in) :: ln0
    integer, intent(in) :: nlep
    logical, dimension(6), intent(in) :: order
    logical, intent(in) :: running
    real(default) :: c_b0, c_b1
    jll_nll = 0._default
    c_b0 = - coeffqed_b0 (0, nlep)
    c_b1 = - coeffqed_b1 (0, nlep)
    if (order(1)) then
       jll_nll(1) = two/(one - x)
    end if
    if (order(2)) then
       jll_nll(2) = (two*(-three - four*log(1-x))) / (-one + x)
    end if
    if (order(3)) then
       jll_nll(3) = (-27._default + 8._default*Pi2 - 72._default*log(1-x) &
            - 48*log(1-x)**2) / (two*(-one + x))
    end if
    if (order(4)) then
       jll_nll(4) = (two*(-one + ln0)) / (one-x) - &
            (four*log(1-x))/(one-x)
    end if
    if (order(5)) then
       if (running) then
          jll_nll(5) = (one - 20._default * nlep/9._default + four*c_b0*Pi &
               - (four*c_b1*Pi)/c_b0 + (four*Pi2)/three + &
               ln0*(6._default-four*c_b0*Pi))/(one-x) - (two*(7._default - &
               four*ln0 - four*c_b0*Pi)*log(1-x))/(one-x) - &
               (12._default*log(1-x)**2)/(one-x)
       else
          jll_nll(5) = (one + 6._default*ln0 - (20._default * nlep)/9._default &
               + (four*Pi2)/three)/(one-x) - (two*(7._default-four*ln0) * &
               log(1-x))/(one-x) - (12._default*log(1-x)**2)/(one-x);
       end if
    end if
    if (order(6)) then
       if (running) then
          jll_nll(6) = two * (-(((8.5_default + (80._default*nlep)/9._default - &
               28._default*c_b0*Pi + (16._default*c_b1*Pi)/c_b0 &
               - (20._default*Pi2)/three + 8._default*c_b0**2*Pi2 - &
               two*ln0*(9._default - 8._default*c_b0*Pi))*log(1-x))/(one-x)) - &
               (6._default*(five - two*ln0 - four*c_b0*Pi) * log(1-x)**2) / &
               (one-x) - (16._default*log(1-x)**3) / (one-x) + (4.5_default - &
               (12._default*c_b1*Pi)/c_b0 - four*c_b0**2*Pi2 + (three + &
               four*c_b1)*Pi2 - nlep*(11._default/3._default - &
               (20._default*c_b0*Pi)/9._default + (four*Pi2)/9._default) - &
               (two*c_b0*Pi*(three + four*Pi2))/three + ln0*(6.75_default - &
               12._default*c_b0*Pi - two*Pi2 + four*c_b0**2*Pi2) - &
               20._default*zeta3) / (one-x))
       else
          jll_nll(6) = two * (-(((8.5_default - 18._default*ln0 + (80._default &
               * nlep)/9._default - (20._default*Pi2)/three) * log(1-x))/(one-x)) - &
               (6._default*(five-two*ln0) * log(1-x)**2)/(one-x) - &
               (16._default*log(1-x)**3)/(one-x) + (4.5_default + three*Pi2 + &
               ln0*(6.75_default - two*Pi2) - nlep*(11._default/3._default + &
               (four*Pi2)/9._default) - 20._default*zeta3) / (one-x))
       end if
    end if
  end subroutine elecbar_asym_p

  module subroutine photbar_asym_p (x, jll_nll, nlep, ln0, order, running)
    real(default), intent(in) :: x
    real(default), dimension(6), intent(out) :: jll_nll
    real(default), intent(in) :: ln0
    integer, intent(in) :: nlep
    logical, dimension(6), intent(in) :: order
    logical, intent(in) :: running
    real(default) :: c_b0, c_b1ob0
    jll_nll = 0._default
    if (running) then
       c_b0 = - coeffqed_b0 (0, nlep)
       c_b1ob0 = coeffqed_b1 (0, nlep) / coeffqed_b0 (0, nlep)
    else
       c_b0 = zero
       c_b1ob0 = zero
    end if
    if (order(1)) then
       jll_nll(1) = one
    end if
    if (order(2)) then
       jll_nll(2) = 1.5_default - (two * nlep)/three +two*log(1-x)
    end if
    if (order(3)) then
       jll_nll(3) = 2.25_default - nlep + four * nlep**2 / 9._default - &
            (two*Pi2)/three + (6._default - (four*nlep)/three)*log(1-x) + &
            four*log(1-x)**2
    end if
    if (order(4)) then
       jll_nll(4) = -one + ln0
    end if
    if (order(5)) then
       jll_nll(5) = -four + (three*ln0)/two - (two*(13._default + &
            three*ln0)*nlep)/9._default - two*c_b1ob0*Pi - 2*c_b0*(-one + &
            ln0)*Pi + (-7._default + two*ln0 - (four*nlep)/three)*log(1-x) &
            - three*log(1-x)**2
    end if
    if (order(6)) then
       jll_nll(6) = -5.625_default - (23._default * nlep)/6._default + &
            (52._default * nlep**2)/27._default + four*c_b0*Pi  - &
            6._default*c_b1ob0*Pi + (40._default*c_b0*nlep*Pi)/9._default + &
            (8._default*c_b1ob0*nlep*Pi)/three + (11._default*Pi2)/6._default &
            - four*c_b0**2*Pi2 + four*c_b0*c_b1ob0*Pi2 + &
            (two*nlep*Pi2)/9._default + ln0*(2.25_default + &
            (four*nlep**2)/9._default - 6._default*c_b0*Pi - (two*Pi2)/three + &
            four*c_b0**2*Pi2 + nlep*(-one + (8._default*c_b0*Pi)/three)) + &
            (-18.5_default + (8._default*nlep**2)/9._default + &
            18._default*c_b0*Pi - 8._default*c_b1ob0*Pi + two*Pi2 + &
            ln0*(6._default - (four*nlep)/three - 8._default*c_b0*Pi) + &
            (four*nlep*(-5._default + two*c_b0*Pi))/three)*log(1-x) + &
            (-18.5_default + four*ln0 - (two*nlep)/three + &
            10._default*c_b0*Pi)*log(1-x)**2 - 6._default*log(1-x)**3 - &
            6._default*zeta3
    end if
  end subroutine photbar_asym_p

  module subroutine recbar_singlet (x, jll_nll, nlep, ln0, order, running)
    real(default), intent(in) :: x
    real(default), dimension(6), intent(out) :: jll_nll
    real(default), intent(in) :: ln0
    integer, intent(in) :: nlep
    logical, dimension(6), intent(in) :: order
    logical, intent(in) :: running
    real(default) :: c_b0, c_b1ob0
    jll_nll = 0._default
    if (running) then
       c_b0 = - coeffqed_b0 (0, nlep)
       c_b1ob0 = coeffqed_b1 (0, nlep) / coeffqed_b0 (0, nlep)
    else
       c_b0 = zero
       c_b1ob0 = zero
    end if
    if (order(1)) then
       jll_nll(1) = -two + two/(one-x)
    end if
    if (order(2)) then
       jll_nll(2) = -two + 6._default/(one-x) - &
            8._default*log(1-x) + (8._default*log(1-x))/(one-x)
    end if
    if (order(3)) then
       jll_nll(3) = 4.5_default + four*Pi2 + (13.5_default - &
            four*Pi2)/(one-x) - 12._default*log(1-x) + &
            (36._default*log(1-x))/(one-x) - 24._default*log(1-x)**2 + &
            (24._default*log(1-x)**2)/(one-x)
    end if
    if (order(4)) then
       jll_nll(4) = two - two*ln0 + (-two+two*ln0)/(one-x) + &
            four*log(1-x) - (four*log(1-x))/(one-x)
    end if
    if (order(5)) then
       jll_nll(5) = (9._default - 12._default*nlep - 18._default*x &
            + 32._default*nlep*x - 36._default*c_b0*Pi*x + &
            36._default*c_b1ob0*Pi*x - 12._default*Pi2*x + &
            18._default*ln0*(-two-x+two*c_b0*Pi*x) - &
            18._default*(-two+(-five + four*ln0 + four*c_b0*Pi)*x)*log(1-x) &
            + 108._default*x*log(1-x)**2)/(9._default*(-one+x))
    end if
    if (order(6)) then
       jll_nll(6) = (-9._default*ln0*(36._default + 16._default*c_b0**2*Pi2*x - &
            (9._default + 8._default*Pi2)*x - 16._default*c_b0*Pi*(two+x)) + &
            two*nlep*(44._default + (22._default + 8._default*Pi2)*x - &
            8._default*c_b0*Pi*(-three + 8._default*x)) + &
            two*(16._default*nlep*(-three + 13._default*x) + 36._default*ln0* &
            (-6._default + (-three + 8._default*c_b0*Pi)*x) + three*(72._default + &
            48._default*c_b0**2*Pi2*x + (-21._default + 96._default*c_b1ob0*Pi - &
            40._default*Pi2)*x - 24._default*c_b0*Pi*(three + four*x)))*log(1-x) - &
            216._default*(-two + (-three + two*ln0 + four*c_b0*Pi)*x)*log(1-x)**2 + &
            576._default*x*log(1-x)**3 + &
            6._default*(-15._default - 8._default*Pi2 - 12._default*x - &
            10._default*Pi2*x + 24._default*c_b0**2*Pi2*x + &
            24._default*c_b1ob0*Pi*(two+x) + two*c_b0*Pi*(-15._default + &
            (21._default - 12._default*c_b1ob0*Pi + 8._default*Pi2)*x) + &
            120._default*x*zeta3)) / (18.*(-one+x))
    end if
  end subroutine recbar_singlet

  module subroutine recbar_nonsinglet (x, jll_nll, nlep, ln0, order, running)
    real(default), intent(in) :: x
    real(default), dimension(6), intent(out) :: jll_nll
    real(default), intent(in) :: ln0
    integer, intent(in) :: nlep
    logical, dimension(6), intent(in) :: order
    logical, intent(in) :: running
    real(default) :: c_b0, c_b1ob0
    jll_nll = 0._default
    if (running) then
       c_b0 = - coeffqed_b0 (0, nlep)
       c_b1ob0 = coeffqed_b1 (0, nlep) / coeffqed_b0 (0, nlep)
    else
       c_b0 = zero
       c_b1ob0 = zero
    end if
    if (order(1)) then
       jll_nll(1) = -two + two/(one - x)
    end if
    if (order(2)) then
       jll_nll(2) = -two + 6._default/(one - x) &
            - 8._default*log(1-x) + (8._default*log(1-x))/(one - x)
    end if
    if (order(3)) then
       jll_nll(3) = 4.5_default + 4._default*Pi2 + (13.5_default &
            - four*Pi2)/(one - x) - 12._default*log(1-x) + &
            (36._default*log(1-x))/(one - x) - 24._default*log(1-x)**2 + &
            (24._default*log(1-x)**2)/(one - x)
    end if
    if (order(4)) then
       jll_nll(4) = two - two*ln0 + (-two + two*ln0)/(one-x) + &
            four*log(1-x) - (four*log(1-x))/(one-x)
    end if
    if (order(5)) then
       jll_nll(5) = -two - two*ln0+(32._default*nlep)/9._default + &
            four*c_b1ob0*Pi - (four*Pi2)/three + one/(one-x) + &
            (6._default*ln0)/(one-x) - (20._default*nlep)/(9._default*(one-x)) - &
            (four*c_b1ob0*Pi)/(one-x) + (four*Pi2)/(three*(one-x)) + &
            10._default*log(1-x) - 8._default*ln0*log(1-x) - &
            (14._default*log(1-x))/(one-x) + (8._default*ln0*log(1-x))/(one-x) + &
            12._default*log(1-x)**2 - (12._default*log(1-x)**2)/(one-x) + &
            c_b0*(-four*Pi + four*ln0*Pi + (four*Pi)/(one-x) - &
            (four*ln0*Pi)/(one-x) - 8._default*Pi*log(1-x) + &
            (8._default*Pi*log(1-x))/(one-x))
    end if
    if (order(6)) then
       jll_nll(6) = -four + (9._default * ln0)/two + (22._default*nlep)/9._default + &
            8._default*c_b1ob0*Pi - (10._default*Pi2)/three + four*ln0*Pi2 + &
            (8._default*nlep*Pi2)/9._default + 9._default/(one-x) + &
            (27._default*ln0)/(two*(one-x)) - (22._default*nlep)/(three*(one-x)) - &
            (24._default*c_b1ob0*Pi)/(one-x) + (6._default*Pi2)/(one-x) - &
            (four*ln0*Pi2)/(one-x) - (8._default*nlep*Pi2)/(9._default*(one-x)) - &
            7._default*log(1-x) - 12._default*ln0*log(1-x) + (208._default* &
            nlep*log(1-x))/9._default + 32._default*c_b1ob0*Pi*log(1-x) - &
            (40._default*Pi2*log(1-x))/three - (17._default*log(1-x))/(one-x) + &
            (36._default*ln0*log(1-x))/(one-x) - (160._default*nlep*log(1-x))/ &
            (9._default*(one-x)) - (32._default*c_b1ob0*Pi*log(1-x))/(one-x) + &
            (40._default*Pi2*log(1-x))/(three*(one-x)) + 36._default*log(1-x)**2 - &
            24._default*ln0*log(1-x)**2 - (60._default*log(1-x)**2)/(one-x) + &
            (24._default*ln0*log(1-x)**2)/(one-x) + 32._default*log(1-x)**3 - &
            (32*log(1-x)**3)/(one-x) + c_b0**2*(8._default*Pi2 - &
            8._default*ln0*Pi2 - (8._default*Pi2)/(one-x) + &
            (8._default*ln0*Pi2)/(one-x) + 16._default*Pi2*log(1-x) - &
            (16._default*Pi2*log(1-x))/(one-x)) + c_b0*(14._default*Pi + &
            8._default*ln0*Pi - (64._default*nlep*Pi)/9._default - &
            8._default*c_b1ob0*Pi2 + (16._default*Pi**3)/three - &
            (24._default*ln0*Pi)/(one-x) + (40._default*nlep*Pi)/(9._default*(one-x)) + &
            (8._default*c_b1ob0*Pi2)/(one-x) + (-four*Pi - &
            (16._default*Pi**3)/three)/(one-x) - 32._default*Pi*log(1-x) + &
            32._default*ln0*Pi*log(1-x) + (56._default*Pi*log(1-x))/(one-x) - &
            (32._default*ln0*Pi*log(1-x))/(one-x) - 48._default*Pi*log(1-x)**2 + &
            (48*Pi*log(1-x)**2)/(1-x)) + 40._default*zeta3 - (40._default*zeta3)/(one-x) 
    end if
  end subroutine recbar_nonsinglet

  module subroutine recbar_photon (x, jll_nll, nlep, ln0, order, running)
    real(default), intent(in) :: x
    real(default), dimension(6), intent(out) :: jll_nll
    real(default), intent(in) :: ln0
    integer, intent(in) :: nlep
    logical, dimension(6), intent(in) :: order
    logical, intent(in) :: running
    real(default) :: c_b0, c_b1ob0
    jll_nll = 0._default
    if (running) then
       c_b0 = - coeffqed_b0 (0, nlep)
       c_b1ob0 = coeffqed_b1 (0, nlep) / coeffqed_b0 (0, nlep)
    else
       c_b0 = zero
       c_b1ob0 = zero
    end if
    if (order(1)) then
       jll_nll(1) = one
    end if
    if (order(2)) then
       jll_nll(2) = 1.5_default - (two*nlep)/three + two*log(1-x)
    end if
    if (order(3)) then
       jll_nll(3) = 2.25_default - nlep + (four*nlep**2)/9._default - &
            (two*Pi2)/three + (6._default - (four*nlep)/three)*log(1-x) + &
            four*log(1-x)**2
    end if
    if (order(4)) then
       jll_nll(4) = -one + ln0
    end if
    if (order(5)) then
       jll_nll(5) = -four + (three*ln0)/two - (two*(13._default + &
            three*ln0)*nlep)/9._default - two*c_b1ob0*Pi - &
            two*c_b0*(-one+ln0)*Pi + (-7._default + two*ln0 - &
            (four*nlep)/three)*log(1-x) - three*log(1-x)**2
    end if
    if (order(6)) then
       jll_nll(6) = -5.625_default - (23._default*nlep)/6._default + &
            (52._default*nlep**2)/27._default + four*c_b0*Pi - &
            6._default*c_b1ob0*Pi + (40._default*c_b0*nlep*Pi)/9._default + &
            (8._default*c_b1ob0*nlep*Pi)/three + (11._default*Pi2)/6._default - &
            four*c_b0**2*Pi2 + four*c_b0*c_b1ob0*Pi2 + (two*nlep*Pi2)/9._default + &
            ln0*(2.25_default + (four*nlep**2)/9._default - 6._default*c_b0*PI - &
            (two*Pi2)/three + four*c_b0**2*Pi2 + nlep*(-one + &
            (8._default*c_b0*Pi)/three)) + (-18.5_default + &
            (8._default*nlep**2)/9._default + 18._default*c_b0*Pi - &
            8._default*c_b1ob0*Pi + two*Pi2 + ln0*(6._default - (four*nlep)/three - &
            8._default*c_b0*Pi) + (four*nlep*(-five + two*c_b0*Pi))/three)*log(1-x) + &
            (-18.5_default + four*ln0 - (two*nlep)/three + 10._default*c_b0*Pi) * &
            log(1-x)**2 - 6._default*log(1-x)**3 - 6._default*zeta3
    end if
  end subroutine recbar_photon

  module subroutine rechat_singlet (x, jll_nll, nlep, ln0, order, running)
    real(default), intent(in) :: x
    real(default), dimension(6), intent(out) :: jll_nll
    real(default), intent(in) :: ln0
    integer, intent(in) :: nlep
    logical, dimension(6), intent(in) :: order
    logical, intent(in) :: running
    real(default) :: c_b0, c_b1ob0
    real(default) :: j6_1, j6_2, j6_3
    jll_nll = 0._default
    if (running) then
       c_b0 = - coeffqed_b0 (0, nlep)
       c_b1ob0 = coeffqed_b1 (0, nlep) / coeffqed_b0 (0, nlep)
    else
       c_b0 = zero
       c_b1ob0 = zero
    end if
     if (order(1)) then
       jll_nll(1) = 1-x
    end if
    if (order(2)) then
       jll_nll(2) = -three-x - (two*nlep*(-one+x)*(four+x*(7+4*x)))/(three*x) - &
            four*(-one+x)*log(1-x) + ((one-four*nlep+(three+four*nlep)*x**2)*log(x))/(-one+x)
    end if
    if (order(3)) then
       jll_nll(3) = -(-64*nlep**2 - 513*x -552*nlep*x + 16*nlep**2*x - 36*Pi2*x + &
            96*nlep*Pi2*x + 378*x**2 + 1104*nlep*x**2 + 96*nlep**2*x**2 + 144*Pi2*x**2 + &
            135*x**3 - 552*nlep*x**3 + 16*nlep**2*x**3 - 108*Pi2*x**3 -96*nlep*Pi2*x**3 - &
            64*nlep**2*x**4 + 432*(-1+x)**2*x*log(1-x)**2 -54*x*log(x) - &
            96*nlep**2*x*log(x) - 432*x**2*log(x) + 288*nlep*x**2*log(x) - &
            162*x**3*log(x) +96*nlep*x**3*log(x) + 96*nlep**2*x**3*log(x) - &
            384*nlep*x**4*log(x) + 18*x*log(x)**2 - 144*nlep*x*log(x)**2 + &
            126*x**3*log(x)**2 + 144*nlep*x**3*log(x)**2 + 24*log(1-x)*((-1+x) * &
            (9*x*(3+x) + 4*nlep*(-4-3*x+3*x**2+4*x**3)) - 18*(x+x**3)*log(x)) + &
            72*(3+8*nlep)*x*(-1+x**2)*polylog(2,x)) / (36._default*(-1+x)*x)
    end if
    if (order(4)) then
       jll_nll(4) = -one + ln0*(1-x) + x - 2*(1-x)*log(1-x)
    end if
    if (order(5)) then
       jll_nll(5) = -(-72*nlep + 48*ln0*nlep+9*x - 54*ln0*x - 176*nlep*x + 36*ln0*nlep*x + &
            36*c_b0*Pi*x - 36*c_b1ob0*Pi*x - 36*c_b0*ln0*Pi*x + 42*Pi2*x + 9*x**2 - &
            18*ln0*x**2 + 368*nlep*x**2 -84*ln0*nlep*x**2 - 36*c_b0*Pi*x**2 + &
            36*c_b1ob0*Pi*x**2 + 36*c_b0*ln0*Pi*x**2 - 42*Pi2*x**2 - 9*x**3 + 54*ln0*x**3 + &
            104*nlep*x**3 - 84*ln0*nlep*x**3 - 36*c_b0*Pi*x**3 + 36*c_b1ob0*Pi*x**3 + &
            36*c_b0*ln0*Pi*x**3 - 18*Pi2*x**3 - 9*x**4 + 18*ln0*x**4 - 296*nlep*x**4 + &
            36*ln0*nlep*x**4 + 36*c_b0*Pi*x**4 - 36*c_b1ob0*Pi*x**4 - 36*c_b0*ln0*Pi*x**4 + &
            18*Pi2*x**4 + 72*nlep*x**5 + 48*ln0*nlep*x**5 + 18*x*log(1-x) + &
            72*ln0*x*log(1-x) + 72*c_b0*Pi*x*log(1-x) + 54*x**2*log(1-x) - &
            72*ln0*x**2*log(1-x) - 72*c_b0*Pi*x**2*log(1-x) - 18*x**3*log(1-x) - &
            72*ln0*x**3*log(1-x) - 72*c_b0*Pi*x**3*log(1-x) - 54*x**4*log(1-x) + &
            72*ln0*x**4*log(1-x) + 72*c_b0*Pi*x**4*log(1-x) - 108*x*log(1-x)**2 + &
            108*x**2*log(1-x)**2 + 108*x**3*log(1-x)**2 - 108*x**4*log(1-x)**2 - &
            96*nlep*log(x) - 18*ln0*x*log(x) - 192*nlep*x*log(x) + 72*ln0*nlep*x*log(x) - &
            18*ln0*x**2*log(x) + 120*nlep*x**2*log(x) + 72*ln0*nlep*x**2*log(x) + &
            18*x**3*log(x) - 54*ln0*x**3*log(x) + 264*nlep*x**3*log(x) - &
            72*ln0*nlep*x**3*log(x) + 18*x**4*log(x) - 54*ln0*x**4*log(x) - &
            48*nlep*x**4*log(x) - 72*ln0*nlep*x**4*log(x) - 96*nlep*x**5*log(x) + &
            72*x**3*log(1-x)*log(x) + 72*x**4*log(1-x)*log(x) + 9*x*log(x)**2 - &
            108*nlep*x*log(x)**2 - 27*x**2*log(x)**2 - 108*nlep*x**2*log(x)**2 + &
            27*x**3*log(x)**2 + 108*nlep*x**3*log(x)**2 - 9*x**4*log(x)**2 + &
            108*nlep*x**4*log(x)**2 + 144*x*log(x)*log(1+x) - 144*x**2*log(x)*log(1+x) - &
            72*x**3*log(x)*log(1+x) + 72*x**4*log(x)*log(1+x) - 108*x*log(1+x)**2 + &
            108*x**2*log(1+x)**2 + 36*(-one+x)*x*(1+x)**2*polylog(2,1-x) + &
            72*x*(2-2*x-x**2+x**3)*polylog(2,-x) - 216*x*polylog(2,1/(1+x)) + &
            216*x**2*polylog(2,1/(1+x))) / (18.*x*(-one+x**2))
    end if
    if (order(6)) then
       j6_1 = -1152*nlep + 32*nlep**2 + 192*ln0*nlep**2 + 1152*c_b1ob0*nlep*Pi + 96*nlep*Pi2 - &
            1350*x +1539*ln0*x - 3660*nlep*x + 1656*ln0*nlep*x + 720*nlep**2*x + 144*ln0*nlep**2*x - &
            1296*c_b1ob0*Pi*x + 864*c_b1ob0*nlep*Pi*x + 72*Pi2*x + 216*ln0*Pi2*x + 120*nlep*Pi2*x + &
            486*x**2 + 405*ln0*x**2 + 5004*nlep*x**2 - 1656*ln0*nlep*x**2 - 176*nlep**2*x**2 - &
            336*ln0*nlep**2*x**2 - 432*c_b1ob0*Pi*x**2 - 2016*c_b1ob0*nlep*Pi*x**2 + 108*Pi2*x**2 - &
            216*ln0*Pi2*x**2 - 216*nlep*Pi2*x**2 + 1350*x**3 - 1539*ln0*x**3 + 4092*nlep*x**3 - &
            1656*ln0*nlep*x**3 - 1328*nlep**2*x**3 - 336*ln0*nlep**2*x**3 + 1296*c_b1ob0*Pi*x**3 - &
            2016*c_b1ob0*nlep*Pi*x**3 - 504*Pi2*x**3 - 216*ln0*Pi2*x**3 - 216*nlep*Pi2*x**3 - &
            486*x**4 - 405*ln0*x**4 - 3852*nlep*x**4 + 1656*ln0*nlep*x**4 + 144*nlep**2*x**4 + &
            144*ln0*nlep**2*x**4 + 432*c_b1ob0*Pi*x**4 + 864*c_b1ob0*nlep*Pi*x**4 + &
            324*Pi2*x**4 + 216*ln0*Pi2*x**4 + 120*nlep*Pi2*x**4 - 432*nlep*x**5 + &
            608*nlep**2*x**5 + 192*ln0*nlep**2*x**5 + 1152*c_b1ob0*nlep*Pi*x**5 + &
            96*nlep*Pi2*x**5 + 288*Pi2*x*(-1+2*x+x**3)*log(2._default) - &
            216*(-1+x)*(1+x)*(3+x*(7+9*x))*log(2._default)**2 - &
            144*x*(-11+x*(16+x+4*x**2))*log(2._default)**3 - &
            432*c_b0**2*Pi2*(-1+x)**2*x*(1+x)*(-1+ln0-2*log(1-x)) + &
            1728*c_b1ob0*Pi*(-1+x)**2*x*(1+x)*log(1-x) + 6*(16*nlep**2*(-3+x)*x*(1+x)*(-1+3*x) + &
            4*nlep*(8+12*ln0*(-4+x*(-3+7*x*(1+x))) + x*(181+x*(-237+x*(-117+229*x)))) - &
            3*(-1+x)*x*(36*ln0*(1+x)*(3+x)+Pi2*(-56-4*x+44*x**2) + 3*(-37+(-36+x)*x - &
            40*log(2._default)**2)))*log(1-x)
       j6_2 = 1512*(-1+x)*x*(1+x)*log(1-x)**2 + 1296*ln0*(-1+x)*x*(1+x)*log(1-x)**2 + &
            1080*(-1+x)*x**2*(1+x)*log(1-x)**2 - 1296*ln0*(-1+x)*x**2*(1+x)*log(1-x)**2 + &
            144*nlep*(-1+x)**2*(1+x)*(4+7*x+4*x**2)*log(1-x)**2 + &
            1728*(-1+x)**2*x*(1+x)*log(1-x)**3 + &
            96*nlep*(4*nlep-9*ln0*x**4+4*(-4-3*ln0+nlep)*x**5)*(log(1-x)-log(x)) - &
            432*c_b1ob0*Pi*x*(1+x)*(1-4*nlep+(3+4*nlep)*x**2)*log(x) - &
            6*(3*ln0*x*(16*nlep**2*(-1+x)*(1+x)**2 + 16*nlep*x*(3+4*x) - 9*(1+x)*(1+x*(8+3*x))) - &
            216*log(2._default) + 2*x*(8*nlep**2*(1+x)*(-7+x*(-13+16*x)) + &
            3*Pi2*(-13+x*(9+(5-9*x)*x)) + 2*nlep*(238 + 18*Pi2*(-1+x)*(1+x)*(2+x) + &
            x*(-83+x*(-152+105*x))) + 144*x*log(2._default) + &
            9*(1+x)*(3-x**2+8*(-1+x)*log(2._default) + 4*(-1+x)**2*log(2._default)**2)))*log(x) + &
            216*x*log(1-x)*log(x) - 576*nlep**2*(-1+x)*x*(1+x)**2*log(1-x)*log(x) - &
            216*x**2*(7+10*x)*log(1-x)*log(x) + 648*ln0*x*(1+x)*(1+3*x**2)*log(1-x)*log(x) + &
            144*nlep*(1+x)*(8+x*(15+8*x**3-9*x*(2+x)+12*ln0*(-1+x**2)))*log(1-x)*log(x) - &
            216*(1+x)*(x-4*nlep*x+(11+4*nlep)*x**3)*log(1-x)**2*log(x) + 648*log(x)**2 + &
            432*x*log(x)**2 - 504*nlep*x*(1+x)*log(x)**2 - 1440*nlep*x**4*(1+x)*log(x)**2 + &
            288*nlep**2*(-1+x)*x*(1+x)**2*log(x)**2 + 72*nlep*x**2*(1+x)*(15+8*x)*log(x)**2 - &
            54*ln0*x*(1+x)*(1-8*nlep+(7+8*nlep)*x**2)*log(x)**2 + &
            108*x**2*(-13+x*(-8+7*x))*log(x)**2 - 216*(-1+x)*x*(-5-x*(5+2*x) + &
            2*nlep*(-1+x+2*x**2))*log(1-x)*log(x)**2 + &
            36*x*(-1+x**2)*(-1+5*x+12*nlep*(2+x))*log(x)**3 + 72*x*(Pi2*(-1+x)*(1+x**2) + &
            (-x**3+2*(-1+x)*(4+x**2)*log(2._default))*log(64._default))*log(1+x) + &
            1728*(-1+x)*x*log(1-x)*log(x/four)*log(1+x) - 1296*(-1+x)*x*log(x)*log(1+x) - &
            1728*(-1+x)*x**2*log(x)*log(1+x) + 432*x*(1+x)*(2+x*(-3+2*x))*log(x)**2*log(1+x) + &
            432*(-3+2*x*(-1+x*(2+x))+ 2*(-1+x)*x**3*log(2._default) - &
            2*(-1+x)*x*log(2-2*x))*log(2*x)*log(1+x) - 1296*(-1+x)*log(1+x)**2 - &
            2484*(-1+x)*x*log(1+x)**2 + 108*(-1+x)*x**2*(29+14*x)*log(1+x)**2 - &
            1728*x**4*log(2._default)*log(1+x)**2 + 432*x*(-13+15*x)*log(x)*log(1+x)**2 + &
            864*x**3*log(4*x)*log(1+x)**2 - 432*x**4*log(x)*(log((1-x)/2.)+log(1+x)) + &
            288*(-1+x)*x*log(1+x)**2*(3*(log(4-4*x)+log(x)) + (-16+x**2)*log(1+x)) + &
            864*(-1+x)*x*(-((-1+x**2)*log(1-x)) + x**2*log(1+x) -log(32*(1+x)))*polylog(2,(1-x)/2.) 
       j6_3 = 648*x*(1+x)*polylog(2,1-x) + 648*ln0*x*(1+x)*(-1+x**2)*polylog(2,1-x) - &
            576*nlep**2*x*(1+x)*(-1+x**2)*polylog(2,1-x) + 288*nlep*(-1+x)*(1+x)*(-2 + &
            6*ln0*x*(1+x)+(-3+x)*x*(3+2*x))*polylog(2,1-x) - 432*(-1+x)*x*(1+x)*(4*(2 + &
            nlep+x+nlep*x)*log(1-x)+(-1+2*nlep*(-2+x)-x)*log(x)+2*(-1+x)*log(1+x))*polylog(2,1-x) - &
            432*(-1+x**2)*(3+x*(2+x*(-1+log(4._default))-log(4._default)) - &
            2*(-1+x)*x*log(1-x))*polylog(2,(-1+x)/(2.*x)) -1296*(-1+x)*polylog(2,-x) - &
            4320*(-1+x)*x*polylog(2,-x) + 216*(-1+x)*x**2*(-1+15*x)*polylog(2,-x) + &
            864*x*((-1+x)*log(1-x) + (1+x**3)*log(x) + (-5+6*x+x**2)*log(1+x))*polylog(2,-x) + &
            1296*(-1+x)*x*polylog(2,1/(1+x)) + 864*(-1+x)*x**2*polylog(2,1/(1+x)) - &
            1728*(-1+x)*x**3*polylog(2,1/(1+x)) - 2592*(-1+x)*x*log(1-x)*polylog(2,1/(1+x)) - &
            864*x*(-4+5*x+x**3)*log(1+x)*polylog(2,1/(1+x)) + 12*c_b0*Pi*(-224*nlep+96*ln0*nlep + &
            90*x-108*ln0*x - 212*nlep*x+72*ln0*nlep*x - 36*c_b1ob0*Pi*x+54*Pi2*x - 36*ln0*x**2 + &
            556*nlep*x**2-168*ln0*nlep*x**2 + 36*c_b1ob0*Pi*x**2 - 54*Pi2*x**2 - 90*x**3 + &
            108*ln0*x**3+292*nlep*x**3 - 168*ln0*nlep*x**3 + 36*c_b1ob0*Pi*x**3 - 30*Pi2*x**3 + &
            36*ln0*x**4 - 332*nlep*x**4 + 72*ln0*nlep*x**4 - 36*c_b1ob0*Pi*x**4 + &
            30*Pi2*x**4 - 80*nlep*x**5 + 96*ln0*nlep*x**5 + 36*x*(-1+x**2)*(4*ln0*(-1+x) - &
            3*(1+x))*log(1-x) + 216*(-1+x)*x*(1+x)*log(1-x)**2 - 216*(-1+x)*x**2*(1+x)*log(1-x)**2 - &
            6*(1+x)*(3*x*(-1-2*x*(1+x)+ln0*(2+6*x**2)) + 4*nlep*(8+x*(9+4*(-3+x)*x*(1+x) + &
            6*ln0*(-1+x**2))))*log(x) + 36*x*(1+x)*(1+5*x**2)*log(1-x)*log(x) - &
            9*(-1+x)**3*x*log(x)**2 + 180*nlep*x*(1+x)*(-1+x**2)*log(x)**2 + &
            72*(-1+x)*x*(-2+x**2)*log(x)*log(1+x) + 108*(-1+x)*x*log(1+x)**2 + &
            72*x*(1+x)*(-1+x**2)*polylog(2,1-x) + 72*(-1+x)*x*(-2+x**2)*polylog(2,-x) + &
            216*(-1+x)*x*polylog(2,1/(1+x))) + &
            216*(-1+x**2)*(-6+x*(-8+5*x)+10*(-1+x)*x*log(1+x))*polylog(2,x/(1+x)) + &
            432*(-1+x)*(1+x)*(-3+x*(-2+x+log(4._default)-x*log(4._default))+2*(-1+x)*x*log(1+x)) * & 
            polylog(2,-1+2/(1+x)) + 864*(-1+x)**2*x*(1+x)*polylog(3,(1-x)/2.) + &
            432*x*(-1+x**2)*(8+3*x+4*nlep*(1+x))*polylog(3,1-x) - &
            864*(-1+x)**2*x*(1+x)*polylog(3,(-1+x)/(2.*x)) - &
            432*x*(6*nlep+x)*(-1+x**2)*polylog(3,(-1+x)/x) + &
            432*x*(1+x)*(7+x*(-16+7*x))*polylog(3,-x) + &
            216*x*(9+20*nlep+3*x+8*nlep*x)*(-1+x**2)*polylog(3,x) + &
            864*x*(1+x)*(3+x*(-4+3*x))*polylog(3,1/(1+x)) + &
            3024*(-1+x)**2*x*(1+x)*polylog(3,x/(1+x)) - &
            864*(-1+x)**2*x*(1+x)*polylog(3,(2*x)/(1+x)) + &
            864*(-1+x)*x*(4+x**2)*polylog(3,(1+x)/2.) - &
            54*x*(-107+16*nlep*(-1+x)*(1+x)*(5+2*x)+x*(99+x*(67+21*x)))*zeta3
       jll_nll(6) = (j6_1 + j6_2 + j6_3) / (108._default*x*(-one+x**2))
    end if
  end subroutine rechat_singlet

  module subroutine rechat_nonsinglet (x, jll_nll, nlep, ln0, order, running)
    real(default), intent(in) :: x
    real(default), dimension(6), intent(out) :: jll_nll
    real(default), intent(in) :: ln0
    integer, intent(in) :: nlep
    logical, dimension(6), intent(in) :: order
    logical, intent(in) :: running
    real(default) :: c_b0, c_b1ob0
    real(default) :: j6_1, j6_2, j6_3
    jll_nll = 0._default
    if (running) then
       c_b0 = - coeffqed_b0 (0, nlep)
       c_b1ob0 = coeffqed_b1 (0, nlep) / coeffqed_b0 (0, nlep)
    else
       c_b0 = zero
       c_b1ob0 = zero
    end if
     if (order(1)) then
       jll_nll(1) = 1-x
    end if
    if (order(2)) then
       jll_nll(2) = ((1 - x)*(3 + x) - 4*(-1 + x)**2*log(1-x) + log(x) + 3*x**2*log(x)) / (-1 + x)
    end if
    if (order(3)) then
       jll_nll(3) = ((-1 + x)*(4*Pi2*(-1 + 3*x) - 3*(19 + 5*x)) - 48*(-1 + x)**2*log(1-x)**2 - &
            2*log(x)*(-3 - 9*x**2 + log(x) + 7*x**2*log(x)) + 24*log(1-x)*(3 - x**2 + &
            2*(1 + x**2)*log(x)) - 48*x*(log(2 - 2*x) - log(2*x)) - 24*(-1 + x**2)*polylog(2,x)) / &
            (four*(-1 + x))
    end if
    if (order(4)) then
       jll_nll(4) = (1-x)*(-1+ln0-2*log(1-x))
    end if
    if (order(5)) then
       jll_nll(5) = (135 + 54*ln0 + 68*nlep - 36*c_b0*Pi + 36*c_b1ob0*Pi + 36*c_b0*ln0*Pi + &
            18*Pi2 - 153*x + 18*ln0*x - 44*nlep*x + 36*c_b0*Pi*x - 36*c_b1ob0*Pi*x - &
            36*c_b0*ln0*Pi*x - 18*Pi2*x - 135*x**2 - 54*ln0*x**2 - 68*nlep*x**2 + 36*c_b0*Pi*x**2 - &
            36*c_b1ob0*Pi*x**2 - 36*c_b0*ln0*Pi*x**2 + 6*Pi2*x**2 + 153*x**3 - 18*ln0*x**3 + &
            44*nlep*x**3 - 36*c_b0*Pi*x**3 + 36*c_b1ob0*Pi*x**3 + 36*c_b0*ln0*Pi*x**3 - &
            6*Pi2*x**3 - 18*log(1-x) - 72*ln0*log(1-x) - 72*c_b0*Pi*log(1-x) - 54*x*log(1-x) + &
            72*ln0*x*log(1-x) + 72*c_b0*Pi*x*log(1-x) + 18*x**2*log(1-x) + 72*ln0*x**2*log(1-x) + &
            72*c_b0*Pi*x**2*log(1-x) + 54*x**3*log(1-x) - 72*ln0*x**3*log(1-x) - &
            72*c_b0*Pi*x**3*log(1-x) + 108*log(1-x)**2 - 108*x*log(1-x)**2 - 108*x**2*log(1-x)**2 + &
            108*x**3*log(1-x)**2 + 72*log(x) + 18*ln0*log(x) + 12*nlep*log(x) + 72*x*log(x) + &
            18*ln0*x*log(x) + 12*nlep*x*log(x) - 90*x**2*log(x) + 54*ln0*x**2*log(x) + &
            12*nlep*x**2*log(x) - 90*x**3*log(x) + 54*ln0*x**3*log(x) + 12*nlep*x**3*log(x) - &
            72*x**2*log(1-x)*log(x) - 72*x**3*log(1-x)*log(x) + 27*log(x)**2 - 9*x*log(x)**2 + &
            9*x**2*log(x)**2 - 27*x**3*log(x)**2 + 144*log(x)*log(1+x) - 144*x*log(x)*log(1+x) - &
            72*x**2*log(x)*log(1+x) + 72*x**3*log(x)*log(1+x) - 108*log(1+x)**2 + 108*x*log(1+x)**2 - &
            36*(-1+x)*(1+x)**2*polylog(2,1-x) + 72*(2-2*x-x**2+x**3)*polylog(2,-x) - &
            216*polylog(2,1/(1+x)) + 216*x*polylog(2,1/(1+x))) / (18._default*(-1+x**2)) 
    end if
    if (order(6)) then
       j6_1 = 24*Pi2*(-1+x)*x + 4*nlep*(63+4*Pi2*(-1+x)-19*x)*(-1+x)*x*(1+x) + &
            144*c_b1ob0*Pi*(-1+x)*x*(1+x)*(3+x) - 12*Pi2*(-1+x)*x**2*(3+17*x) - &
            18*(-1+x)*x*(1+x)*(-73+57*x) + 9*ln0*(-1+x)*x*(1+x)*(8*Pi2*(-1+x)-3*(19+5*x)) + &
            24*log(2._default)*(-4*Pi2*x*(-1+2*x+x**3)+2*x*(-11+x*(16+x+4*x**2))*log(2._default)**2 + &
            (-1+x)*(1+x)*(3+x*(7+9*x))*log(8._default)) - &
            144*c_b0**2*Pi2*(-1+x)**2*x*(1+x)*(-1+ln0-2*log(1-x)) + &
            576*c_b1ob0*Pi*(-1+x)**2*x*(1+x)*log(1-x) - &
            2*(-1+x)*x*(108*ln0*(1+x)*(3+x)+32*nlep*(11+(3-8*x)*x)+3*(81-27*x*(4+7*x) + &
            4*Pi2*(-6+x+9*x**2)+120*log(2._default)**2))*log(1-x) + &
            504*(-1+x)*x*(1+x)*log(1-x)**2 - 432*ln0*(-1+x)**2*x*(1+x)*log(1-x)**2 + &
            360*(-1+x)*x**2*(1+x)*log(1-x)**2 + 576*(-1+x)**2*x*(1+x)*log(1-x)**3 - &
            540*x*log(x) + 54*ln0*x*log(x) - 152*nlep*x*log(x) - 144*c_b1ob0*Pi*x*log(x) - &
            84*Pi2*x*log(x) + 612*x**2*log(x) + 486*ln0*x**2*log(x) + 136*nlep*x**2*log(x) - &
            144*c_b1ob0*Pi*x**2*log(x) + 132*Pi2*x**2*log(x) + 468*x**3*log(x) + &
            594*ln0*x**3*log(x)-24*nlep*x**3*log(x) - 432*c_b1ob0*Pi*x**3*log(x) + &
            84*Pi2*x**3*log(x) - 684*x**4*log(x)+162*ln0*x**4*log(x) - &
            312*nlep*x**4*log(x) - 432*c_b1ob0*Pi*x**4*log(x) - 36*Pi2*x**4*log(x) - &
            432*log(2._default)*log(x) - 288*x*log(2._default)*log(x) + &
            576*x**2*log(2._default)*log(x)+288*x**3*log(2._default)*log(x) - &
            144*x**4*log(2._default)*log(x)
       j6_2 = 144*x*log(2._default)**2*log(x) - 144*x**2*log(2._default)**2*log(x) - &
            144*x**3*log(2._default)**2*log(x) + 144*x**4*log(2._default)**2*log(x) + &
            648*x*(1+x)*log(1-x)*log(x) + 216*ln0*x*(1+x)*log(1-x)*log(x) + &
            96*nlep*x*(1+x)*log(1-x)*log(x) - 576*x**2*(1+x)*log(1-x)*log(x) - &
            720*x**3*(1+x)*log(1-x)*log(x) + 648*ln0*x**3*(1+x)*log(1-x)*log(x) + &
            96*nlep*x**3*(1+x)*log(1-x)*log(x) - 72*x*(1+x)*(1+11*x**2)*log(1-x)**2*log(x) - &
            216*log(x)**2-288*x*log(x)**2 - 18*ln0*x*log(x)**2-24*nlep*x*log(x)**2 + &
            108*x**2*log(x)**2 - 18*ln0*x**2*log(x)**2-24*nlep*x**2*log(x)**2 + &
            288*x**3*log(x)**2 - 126*ln0*x**3*log(x)**2 - 72*nlep*x**3*log(x)**2 - &
            36*x**4*log(x)**2 - 126*ln0*x**4*log(x)**2 - 72*nlep*x**4*log(x)**2 + &
            72*(-1+x)*x*(1+x*(5+2*x))*log(1-x)*log(x)**2 - 36*x*log(x)**3 - &
            12*x**2*log(x)**3+36*x**3*log(x)**3 + 12*x**4*log(x)**3 - &
            576*(-1+x)*x*log(1-x)*log(x/four)*log(1+x) + 432*(-1+x)*x*log(x)*log(1+x) + &
            576*(-1+x)*x**2*log(x)*log(1+x) - 144*x*(1+x)*(2+x*(-3+2*x))*log(x)**2*log(1+x) - &
            24*(-1+x)*(Pi2*(x+x**3)+12*x*(4+x**2)*log(2._default)**2 + &
            6*(3+x*(5+x+x**2*(-1+log(4._default)))-2*x*log(2-2*x))*log(2*x))*log(1+x) + &
            432*(-1+x)*log(1+x)**2+828*(-1+x)*x*log(1+x)**2 - 36*(-1+x)*x**2*(29+14*x)*log(1+x)**2
       j6_3 = 576*x**4*log(2._default)*log(1+x)**2 - 144*x*(-13+15*x)*log(x)*log(1+x)**2 - &
            288*(-1+x)*x*(log(4-4*x)+log(x))*log(1+x)**2 - 288*x**3*log(4*x)*log(1+x)**2 - &
            96*(-4+x)*(-1+x)*x*(4+x)*log(1+x)**3 + &
            288*(-1+x)*x*(log(32._default)+(-1+x**2)*log(1-x)-(-1+x**2)*log(1+x))*polylog(2,(1-x)/two) + &
            792*x*(1+x)*polylog(2,1-x) - 216*ln0*x*(1+x)*polylog(2,1-x) - &
            576*x**3*(1+x)*polylog(2,1-x) + 216*ln0*x**3*(1+x)*polylog(2,1-x) + &
            72*x*(1+x)*log((-1+x)**8)*polylog(2,1-x) - 144*x*(1+x)*log(x)*polylog(2,1-x) + &
            288*x*(1+x)*log(1+x)*polylog(2,1-x) + &
            144*x**2*(1+x)*((4-8*x)*log(1-x)+x*log(x)+2*(-2+x)*log(1+x))*polylog(2,1-x) + &
            144*(-1+x**2)*(3+x*(2+x*(-1+log(4._default))-log(4._default)) - &
            2*(-1+x)*x*log(1-x))*polylog(2,(-1+x)/(two*x)) + 432*(-1+x)*polylog(2,-x) + &
            1440*(-1+x)*x*polylog(2,-x) - 72*(-1+x)*x**2*(-1+15*x)*polylog(2,-x) - &
            288*x*((-1+x)*log(1-x)+(1+x**3)*log(x)+(-5+6*x+x**2)*log(1+x))*polylog(2,-x) - &
            432*(-1+x)*x*polylog(2,1/(1+x)) - 288*(-1+x)*x**2*polylog(2,1/(1+x)) + &
            576*(-1+x)*x**3*polylog(2,1/(1+x)) + 864*(-1+x)*x*log(1-x)*polylog(2,1/(1+x)) + &
            288*x*(-4+5*x+x**3)*log(1+x)*polylog(2,1/(1+x)) + &
            4*c_b0*Pi*x*(-36*c_b1ob0*Pi*(-1+x)*(-1+x**2)-2*(1-x)*(18*ln0*(1+x)*(3+x) - &
            2*nlep*(1+x)*(-17+11*x)+3*(9+Pi2-15*x+3*(-8+Pi2)*x**2))+36*(3-3*(-1+x)*x + &
            4*ln0*(-1+x)**2*(1+x))*log(1-x)-216*(-1+x)**2*(1+x)*log(1-x)**2 - &
            108*x**3*(log(1-x)-log(x))-54*log(x)-36*ln0*log(x) - &
            12*nlep*log(x)-18*x*log(x)-36*ln0*x*log(x)-12*nlep*x*log(x)+144*x**2*log(x) - &
            108*ln0*x**2*log(x)-12*nlep*x**2*log(x)-108*ln0*x**3*log(x)-12*nlep*x**3*log(x) + &
            36*(1+x)*log(1-x)*log(x) + &
            180*x**2*(1+x)*log(1-x)*log(x)-27*log(x)**2+9*x*log(x)**2-9*x**2*log(x)**2 + &
            27*x**3*log(x)**2 - 72*(-1+x)*(-2+x**2)*log(x)*log(1+x)-108*(-1+x)*log(1+x)**2 - &
            72*(1+x)*polylog(2,1-x) + 72*x**2*(1+x)*polylog(2,1-x)-72*(-1+x)*(-2+x**2)*polylog(2,-x) - &
            216*(-1+x)*polylog(2,one/(1+x))) - &
            72*(-1+x**2)*(-6+x*(-8+5*x)+10*(-1+x)*x*log(1+x))*polylog(2,x/(1+x)) + &
            144*(-1+x**2)*(3+x*(2+x*(-1+log(4._default))-log(4._default)) - &
            2*(-1+x)*x*log(1+x))*polylog(2,-one+two/(1+x)) - 288*(-1+x)**2*x*(1+x)*polylog(3,(1-x)/two) + &
            144*(-1+x)*x*(1+x)*(4+7*x)*polylog(3,1-x) + 288*(-1+x)**2*x*(1+x)*polylog(3,(-1+x)/(two*x)) - &
            144*x**2*(-1+x**2)*polylog(3,(-1+x)/x) - 144*x*(1+x)*(7+x*(-16+7*x))*polylog(3,-x) + &
            72*(-1+x)*x*(1+x)*(5+7*x)*polylog(3,x) + 288*x*(-3+x+x**2-3*x**3)*polylog(3,one/(1+x)) - &
            1008*(-1+x)**2*x*(1+x)*polylog(3,x/(1+x)) + 288*(-1+x)**2*x*(1+x)*polylog(3,(2*x)/(1+x)) - &
            288*(-1+x)*x*(4+x**2)*polylog(3,(1+x)/two) + 18*x*(29+x*(59+x*(-69+61*x)))*zeta3
       jll_nll(6) = (j6_1 + j6_2 + j6_3) / (36._default*x*(-1+x**2))
    end if
  end subroutine rechat_nonsinglet

  module subroutine rechat_photon (x, jll_nll, nlep, ln0, order, running)
    real(default), intent(in) :: x
    real(default), dimension(6), intent(out) :: jll_nll
    real(default), intent(in) :: ln0
    integer, intent(in) :: nlep
    logical, dimension(6), intent(in) :: order
    logical, intent(in) :: running
    real(default) :: c_b0, c_b1ob0
    real(default) :: j6_1, j6_2, j6_3
    jll_nll = 0._default
    if (running) then
       c_b0 = - coeffqed_b0 (0, nlep)
       c_b1ob0 = coeffqed_b1 (0, nlep) / coeffqed_b0 (0, nlep)
    else
       c_b0 = zero
       c_b1ob0 = zero
    end if
     if (order(1)) then
       jll_nll(1) = -3 + 2/x + x
    end if
    if (order(2)) then
       jll_nll(2) = -((-1+x)*(4*nlep*(-2+x)+3*x)-12*(2-3*x+x**2)*log(1-x) + 6*(-2+x)*x*log(x)) / &
            (6._default*x)
    end if
    if (order(3)) then
       jll_nll(3) = (-99 + 1068*nlep - 48*nlep**2 + 72*Pi2 - (496*nlep)/x + (32*nlep**2)/x + &
            99*x - 636*nlep*x + 16*nlep**2*x - 24*Pi2*x + 64*nlep*x**2 + &
            (144*(2-3*x+x**2)*log(1-x)**2)/x - 180*log(x) - 48*nlep*log(x) - (192*nlep*log(x))/x + &
            72*x*log(x) + 240*nlep*x*log(x) - 36*log(x)**2 + 144*nlep*log(x)**2 + 18*x*log(x)**2 - &
            72*nlep*x*log(x)**2 - (24*log(1-x)*((-1+x)*(2*nlep*(-2+x)+3*x) + &
            6*(2-2*x+x**2)*log(x))) / x - (288*polylog(2,x))/x) / 36._default
    end if
    if (order(4)) then
       jll_nll(4) = ((-1+ln0)*(2-3*x+x**2) - 2*(2-2*x+x**2)*log(x)) / x
    end if
    if (order(5)) then
       jll_nll(5) = (45 + 9*ln0 + 108*nlep + 36*ln0*nlep + 108*c_b1ob0*Pi - (56*nlep)/x - &
            (24*ln0*nlep)/x - (72*c_b1ob0*Pi)/x - 45*x - 9*ln0*x - 52*nlep*x - 12*ln0*nlep*x - &
            36*c_b1ob0*Pi*x - (54*(2-3*x+x**2)*log(1-x)**2)/x + 36*ln0*log(x) - 48*nlep*log(x) + &
            (48*nlep*log(x))/x + 45*x*log(x) - 18*ln0*x*log(x) + 24*nlep*x*log(x) - 18*log(x)**2 + &
            9*x*log(x)**2 + (6*log(1-x)*((-1+x)*(12+8*nlep+6*ln0*(-2+x)-9*x-4*nlep*x) + &
            6*(-2+x)*x*log(x)))/x - (36*c_b0*Pi*((-1+ln0)*(2-3*x+x**2) - 2*(2-2*x+x**2)*log(x))) / x + &
            36*(-2+x)*polylog(2,1-x))/18._default
    end if
    if (order(6)) then
       j6_1 = -2808+2176*nlep + 2976*ln0*nlep - 448*nlep**2 - 192*ln0*nlep**2 - 768*c_b0*nlep*Pi - &
            1152*c_b1ob0*nlep*Pi - 1152*c_b0*ln0*nlep*Pi + 144*Pi2 + 1728*c_b0**2*Pi2 - &
            1728*c_b0*c_b1ob0*Pi2 + 288*ln0*Pi2 - 1728*c_b0**2*ln0*Pi2 - 96*nlep*Pi2 + &
            3645*x+594*ln0*x + 12420*nlep*x - 6408*ln0*nlep*x + 864*nlep**2*x + 288*ln0*nlep**2*x + &
            648*c_b0*Pi*x + 432*c_b1ob0*Pi*x + 432*c_b0*ln0*Pi*x + 1728*c_b0*nlep*Pi*x + &
            1728*c_b1ob0*nlep*Pi*x + 1728*c_b0*ln0*nlep*Pi*x + 324*Pi2*x - 2592*c_b0**2*Pi2*x + &
            2592*c_b0*c_b1ob0*Pi2*x - 432*ln0*Pi2*x + 2592*c_b0**2*ln0*Pi2*x + 144*nlep*Pi2*x + &
            1971*x**2 - 594*ln0*x**2 - 16324*nlep*x**2 + 840*ln0*nlep*x**2 + 32*nlep**2*x**2 + &
            96*ln0*nlep**2*x**2 - 648*c_b0*Pi*x**2 - 432*c_b1ob0*Pi*x**2 - 432*c_b0*ln0*Pi*x**2 - &
            192*c_b0*nlep*Pi*x**2 + 576*c_b1ob0*nlep*Pi*x**2 + 576*c_b0*ln0*nlep*Pi*x**2 - &
            180*Pi2*x**2-864*c_b0**2*Pi2*x**2 + 864*c_b0*c_b1ob0*Pi2*x**2 - 144*ln0*Pi2*x**2 + &
            864*c_b0**2*ln0*Pi2*x**2 + 48*nlep*Pi2*x**2 - 3645*x**3-594*ln0*x**3 - &
            12868*nlep*x**3+6024*ln0*nlep*x**3 - 864*nlep**2*x**3 - 288*ln0*nlep**2*x**3 - &
            648*c_b0*Pi*x**3 - 432*c_b1ob0*Pi*x**3 - 432*c_b0*ln0*Pi*x**3 - 1728*c_b0*nlep*Pi*x**3 - &
            1728*c_b1ob0*nlep*Pi*x**3 - 1728*c_b0*ln0*nlep*Pi*x**3 - 468*Pi2*x**3 + &
            2592*c_b0**2*Pi2*x**3 - 2592*c_b0*c_b1ob0*Pi2*x**3 + 432*ln0*Pi2*x**3 - &
            2592*c_b0**2*ln0*Pi2*x**3 - 144*nlep*Pi2*x**3 + 837*x**4 + 594*ln0*x**4 + 14148*nlep*x**4 - &
            3816*ln0*nlep*x**4 + 416*nlep**2*x**4 + 96*ln0*nlep**2*x**4 + 648*c_b0*Pi*x**4 + &
            432*c_b1ob0*Pi*x**4 + 432*c_b0*ln0*Pi*x**4 + 960*c_b0*nlep*Pi*x**4 + &
            576*c_b1ob0*nlep*Pi*x**4 + 576*c_b0*ln0*nlep*Pi*x**4 + 180*Pi2*x**4 - 864*c_b0**2*Pi2*x**4 + &
            864*c_b0*c_b1ob0*Pi2*x**4 - 144*ln0*Pi2*x**4 + 864*c_b0**2*ln0*Pi2*x**4 + 48*nlep*Pi2*x**4 + &
            448*nlep*x**5 + 384*ln0*nlep*x**5 - 864*Pi2*x*log(2._default) + &
            864*Pi2*x**3*log(2._default) + 2016*x*log(2._default)**3 - 2016*x**3*log(2._default)**3 - &
            384*nlep**2*log(2-2*x) + 432*ln0*x**2*log(2-2*x) + 192*nlep**2*x**2*log(2-2*x) - &
            432*ln0*x**4*log(2-2*x) + 192*nlep**2*x**4*log(2-2*x) + 1080*log(1-x) + 2304*nlep*log(1-x) + &
            576*ln0*nlep*log(1-x) - 864*c_b0*Pi*log(1-x) + 3456*c_b1ob0*Pi*log(1-x) + &
            3456*c_b0*ln0*Pi*log(1-x) - 1152*c_b0*nlep*Pi*log(1-x) - 576*Pi2*log(1-x) - 3564*x*log(1-x) - &
            432*ln0*x*log(1-x) - 3744*nlep*x*log(1-x) - 864*ln0*nlep*x*log(1-x) + &
            576*nlep**2*x*log(1-x) + 1296*c_b0*Pi*x*log(1-x) - 5184*c_b1ob0*Pi*x*log(1-x) - &
            5184*c_b0*ln0*Pi*x*log(1-x) + 1728*c_b0*nlep*Pi*x*log(1-x) + 1008*Pi2*x*log(1-x)
       j6_2 = 1404*x**2*log(1-x) - 864*nlep*x**2*log(1-x) - 288*ln0*nlep*x**2*log(1-x) + &
            432*c_b0*Pi*x**2*log(1-x) - 1728*c_b1ob0*Pi*x**2*log(1-x) - 1728*c_b0*ln0*Pi*x**2*log(1-x) + &
            576*c_b0*nlep*Pi*x**2*log(1-x) + 144*Pi2*x**2*log(1-x) + 3564*x**3*log(1-x) + &
            432*ln0*x**3*log(1-x) + 3744*nlep*x**3*log(1-x) + 864*ln0*nlep*x**3*log(1-x) - &
            576*nlep**2*x**3*log(1-x) - 1296*c_b0*Pi*x**3*log(1-x) + 5184*c_b1ob0*Pi*x**3*log(1-x) + &
            5184*c_b0*ln0*Pi*x**3*log(1-x) - 1728*c_b0*nlep*Pi*x**3*log(1-x) - 1008*Pi2*x**3*log(1-x) - &
            2484*x**4*log(1-x) - 1440*nlep*x**4*log(1-x) - 288*ln0*nlep*x**4*log(1-x) + &
            432*c_b0*Pi*x**4*log(1-x) - 1728*c_b1ob0*Pi*x**4*log(1-x) - 1728*c_b0*ln0*Pi*x**4*log(1-x) + &
            576*c_b0*nlep*Pi*x**4*log(1-x) + 432*Pi2*x**4*log(1-x) + 1728*log(1-x)**2 - &
            1728*ln0*log(1-x)**2 + 288*nlep*log(1-x)**2 - 4320*c_b0*Pi*log(1-x)**2 - &
            2700*x*log(1-x)**2 + 2592*ln0*x*log(1-x)**2 - 432*nlep*x*log(1-x)**2 + &
            6480*c_b0*Pi*x*log(1-x)**2 - 756*x**2*log(1-x)**2 + 864*ln0*x**2*log(1-x)**2 - &
            144*nlep*x**2*log(1-x)**2 + 2160*c_b0*Pi*x**2*log(1-x)**2 + 2700*x**3*log(1-x)**2 - &
            2592*ln0*x**3*log(1-x)**2 + 432*nlep*x**3*log(1-x)**2 - 6480*c_b0*Pi*x**3*log(1-x)**2 - &
            972*x**4*log(1-x)**2 + 864*ln0*x**4*log(1-x)**2 - 144*nlep*x**4*log(1-x)**2 + &
            2160*c_b0*Pi*x**4*log(1-x)**2 + 1728*log(1-x)**3 - 3024*x*log(1-x)**3 - &
            432*x**2*log(1-x)**3 + 3024*x**3*log(1-x)**3 - 1296*x**4*log(1-x)**3 - 2880*nlep*log(x) + &
            1152*ln0*nlep*log(x) + 2304*c_b0*nlep*Pi*log(x) + 3456*c_b0**2*Pi2*log(x) - &
            2916*x*log(x) + 1080*ln0*x*log(x) + 18288*nlep*x*log(x) + 288*ln0*nlep*x*log(x) - &
            384*nlep**2*x*log(x) - 864*c_b0*Pi*x*log(x) + 1728*c_b1ob0*Pi*x*log(x) + &
            1728*c_b0*ln0*Pi*x*log(x) - 2304*c_b0*nlep*Pi*x*log(x)+288*Pi2*x*log(x) - &
            3456*c_b0**2*Pi2*x*log(x) - 702*x**2*log(x) + 3960*nlep*x**2*log(x) - &
            2592*ln0*nlep*x**2*log(x) + 648*c_b0*Pi*x**2*log(x) - 864*c_b1ob0*Pi*x**2*log(x) - &
            864*c_b0*ln0*Pi*x**2*log(x) - 1152*c_b0*nlep*Pi*x**2*log(x) - 864*Pi2*x**2*log(x) - &
            1728*c_b0**2*Pi2*x**2*log(x) + 2916*x**3*log(x) - 1080*ln0*x**3*log(x) - &
            17520*nlep*x**3*log(x) - 288*ln0*nlep*x**3*log(x) + 384*nlep**2*x**3*log(x) + &
            864*c_b0*Pi*x**3*log(x) - 1728*c_b1ob0*Pi*x**3*log(x) - 1728*c_b0*ln0*Pi*x**3*log(x) + &
            2304*c_b0*nlep*Pi*x**3*log(x) - 864*Pi2*x**3*log(x) + 3456*c_b0**2*Pi2*x**3*log(x) + &
            702*x**4*log(x) - 1080*nlep*x**4*log(x) + 1440*ln0*nlep*x**4*log(1-x) - &
            648*c_b0*Pi*x**4*log(x) + 864*c_b1ob0*Pi*x**4*log(x) + 864*c_b0*ln0*Pi*x**4*log(x) - &
            1152*c_b0*nlep*Pi*x**4*log(x) + 288*Pi2*x**4*log(x) - 1728*c_b0**2*Pi2*x**4*log(x) - &
            768*nlep*x**5*log(x) - 1296*x*log(1-x)*log(x) - 1728*ln0*x*log(1-x)*log(x) - &
            3456*c_b0*Pi*x*log(1-x)*log(x) - 648*x**2*log(1-x)*log(x) + &
            864*ln0*x**2*log(1-x)*log(x) + 1728*c_b0*Pi*x**2*log(1-x)*log(x) + &
            1296*x**3*log(1-x)*log(x) + 1728*ln0*x**3*log(1-x)*log(x)
       j6_3 =  3456*c_b0*Pi*x**3*log(1-x)*log(x) + 648*x**4*log(1-x)*log(x) - &
            864*ln0*x**4*log(1-x)*log(x) - 1728*c_b0*Pi*x**4*log(1-x)*log(x) + &
            3024*x*log(1-x)**2*log(x) - 1512*x**2*log(1-x)**2*log(x) - &
            3024*x**3*log(1-x)**2*log(x) + 1512*x**4*log(1-x)**2*log(x) - 1152*nlep*log(x)**2 - &
            432*x*log(x)**2 + 216*ln0*x*log(x)**2 + 1728*nlep*x*log(x)**2 - &
            864*ln0*nlep*x*log(x)**2 - 432*c_b0*Pi*x*log(x)**2 - 432*x**2*log(x)**2 - &
            108*ln0*x**2*log(x)**2 + 2016*nlep*x**2*log(x)**2 + 432*ln0*nlep*x**2*log(x)**2 + &
            216*c_b0*Pi*x**2*log(x)**2 + 432*x**3*log(x)**2 - 216*ln0*x**3*log(x)**2 - &
            1728*nlep*x**3*log(x)**2 + 864*ln0*nlep*x**3*log(x)**2 + 432*c_b0*Pi*x**3*log(x)**2 + &
            432*x**4*log(x)**2 + 108*ln0*x**4*log(x)**2 - 864*nlep*x**4*log(x)**2 - &
            432*ln0*nlep*x**4*log(x)**2 - 216*c_b0*Pi*x**4*log(x)**2 - 2592*x*log(1-x)*log(x)**2 + &
            3024*x**2*log(1-x)*log(x)**2 + 4320*x**3*log(1-x)*log(x)**2 - &
            1296*x**4*log(1-x)*log(x)**2 - 144*x*log(x)**3 + 1152*nlep*x*log(x)**3 - &
            936*x**2*log(x)**3 - 576*nlep*x**2*log(x)**3 - 432*x**3*log(x)**3 - &
            1152*nlep*x**3*log(x)**3 + 360*x**4*log(x)**3 + 576*nlep*x**4*log(x)**3 + &
            384*nlep**2*log(2*x) - 432*ln0*x**2*log(2*x) - 192*nlep**2*x**2*log(2*x) + &
            432*ln0*x**4*log(2*x) - 192*nlep**2*x**4*log(2*x) + 576*Pi2*log(1+x) + &
            2016*Pi2*x*log(1+x) - 2016*Pi2*x**2*log(1+x) - 576*Pi2*x**3*log(1+x) + &
            2592*log(x)*log(1+x) + 2592*x*log(x)*log(1+x) - 864*x**2*log(x)*log(1+x) - &
            3456*x**3*log(x)*log(1+x) - 864*x**4*log(x)*log(1+x) + 864*log(x)**2*log(1+x) + &
            864*x*log(x)**2*log(1+x) - 432*x**2*log(x)**2*log(1+x) - 864*x**3*log(x)**2*log(1+x) - &
            432*x**4*log(x)**2*log(1+x)+1296*log(1+x)**2+3456*x*log(1+x)**2 - &
            648*x**2*log(1+x)**2-3456*x**3*log(1+x)**2 - 648*x**4*log(1+x)**2+864*log(x)*log(1+x)**2 + &
            10368*x*log(x)*log(1+x)**2 - 9504*x**2*log(x)*log(1+x)**2 - &
            1728*x**3*log(x)*log(1+x)**2 - 1440*log(1+x)**3 - 4608*x*log(1+x)**3 + &
            3888*x**2*log(1+x)**3 + 288*x**3*log(1+x)**3 + 1872*x**4*log(1+x)**3 - &
            144*(-1+x**2)*(12-12*ln0+4*nlep-15*x-4*nlep*x-24*c_b0*Pi*x+6*x**2+2*nlep*x**2 + &
            12*c_b0*Pi*x**2 + (12+18*x-9*x**2)*log(1-x)+3*(-6+x)*x*log(x))*polylog(2,1-x) + &
            432*(-1+x)*(-12-30*x-25*x**2-3*x**3-8*(1+x)**2*log(x) + &
            2*x*(-12+x+3*x**2)*log(1+x))*polylog(2,-x) + 1728*x*polylog(2,1/(1+x)) - &
            864*x**3*polylog(2,1/(1+x)) - 864*x**4*polylog(2,1/(1+x)) - &
            10368*x*log(1+x)*polylog(2,1/(1+x)) + 7776*x**2*log(1+x)*polylog(2,1/(1+x)) + &
            1728*x**3*log(1+x)*polylog(2,1/(1+x)) + 864*x**4*log(1+x)*polylog(2,1/(1+x)) + &
            2592*polylog(2,x/(1+x)) + 5184*x*polylog(2,x/(1+x)) - 1296*x**2*polylog(2,x/(1+x)) - &
            6048*x**3*polylog(2,x/(1+x)) - 432*x**4*polylog(2,x/(1+x)) + &
            1728*log(x)*polylog(2,x/(1+x)) - 1728*x**2*log(x)*polylog(2,x/(1+x)) + &
            3456*x*log(1+x)*polylog(2,x/(1+x)) - 3456*x**2*log(1+x)*polylog(2,x/(1+x)) - &
            3456*x**3*log(1+x)*polylog(2,x/(1+x)) + 3456*x**4*log(1+x)*polylog(2,x/(1+x)) - &
            4320*x*polylog(3,1-x) + 6480*x**2*polylog(3,1-x) + 7776*x**3*polylog(3,1-x) - &
            3024*x**4*polylog(3,1-x) + 5184*x**2*polylog(3,(-1+x)/x) + 3456*x**3*polylog(3,(-1+x)/x) - &
            1728*x**4*polylog(3,(-1+x)/x) + 1728*polylog(3,x) - 4320*x*polylog(3,x) + &
            5616*x**2*polylog(3,x) + 7776*x**3*polylog(3,x) - 3888*x**4*polylog(3,x) - &
            432*polylog(3,x**2) - 432*x*polylog(3,x**2) + 216*x**2*polylog(3,x**2) + &
            432*x**3*polylog(3,x**2) + 216*x**4*polylog(3,x**2) + 3456*x*polylog(3,x/(1+x)) - &
            1728*x**2*polylog(3,x/(1+x)) - 3456*x**3*polylog(3,x/(1+x)) + 1728*x**4*polylog(3,x/(1+x)) - &
            3456*zeta3 + 2592*x*zeta3 - 864*x**2*zeta3 - 6048*x**3*zeta3 + 864*x**4*zeta3
       jll_nll(6) = (j6_1 + j6_2 + j6_3) / (216._default*x*(-1+x**2))
    end if
  end subroutine rechat_photon

  function sum_rm (x, al0_2pi, ln0, cc1, cc2, cc3, cc4, cc5, &
       k, m1, m2) result (s_rm)
    real(default) :: s_rm
    real(default), intent(in) :: x, k, al0_2pi, ln0, m1, m2
    real(default), intent(in) :: cc1, cc2, cc3, cc4, cc5
    real(default) :: f0, f1, f2, fac, den
    real(default) :: mf1, mf2, mf3, mf4, mf5
    real(default) :: rr1, rr2, rr3, rr4, rr5
    f0 = two - Pi2/three + three/two * ln0
    f1 = two*(one - ln0)
    f2 = - two
    den = m2 - m1 * log(1-x)
    fac = exp(-eulerc*k) * (1-x)**k / gamma(1+k)
    mf1 = one/den - (Pi2*k - 6._default*zeta3 * k**2)*m1/6._default / den**2 &
         - (30._default*Pi2 - 360._default*zeta3*k + &
         Pi**4 * k**2)*m1**2/180._default / den**3
    mf2 = one
    mf3 = - log(1-x) + Pi2*k/6._default - zeta3 * k**2
    mf4 = log(1-x)**2 - Pi2/6._default + k*(-Pi2/three*log(1-x) + two*zeta3) + &
         k**2 * (two*zeta3*log(1-x) - Pi**4/180._default)
    mf5 = - log(1-x)**3 + Pi2/two*log(1-x) - two*zeta3 + &
         k*(Pi2/two*log(1-x)**2 - 6._default*zeta3*log(1-x) - &
         Pi**4/60._default) + k**2 * (-three*zeta3*log(1-x)**2 + &
         Pi**4/60._default*log(1-x) + three/two*Pi2*zeta3 - 12._default*zeta5)
    rr1 = (cc3 - cc4*cc2 + cc4**2 * cc1) * &
         (one + al0_2pi * (f0 - cc4*f1 + cc4**2 * f2))
    rr2 = (cc2 - cc4*cc1) / cc5 + al0_2pi/cc5 * (cc2*f0 + cc3*f1 - &
         cc4*(cc1*f0 + cc2*f1 + cc3*f2) + cc4**2 * (cc1*f1 + cc2*f2) - &
         cc4**3 * cc1 * f2)
    rr3 = cc1/cc5 + al0_2pi/cc5 * (cc1*f0 + cc2*f1 + cc3*f2 - &
         cc4*(cc1*f1 + cc2*f2) + cc4**2 * cc1 * f2)
    rr4 = al0_2pi * 1/cc5 * (cc1*f1 + cc2*f2 - cc4*cc1*f2)
    rr5 = al0_2pi * cc1/cc5 * f2
    s_rm = fac * (rr1*mf1 + rr2*mf2 + rr3*mf3 + rr4*mf4 + rr5*mf5)
  end function sum_rm

  module function t_alpha (epdf, scale) result (t)
    real(default) :: t
    type(qed_pdf_t), intent(in) :: epdf
    real(default), intent(in) :: scale
    real(default) :: alphamu, alpharef 
    select type (alpha => epdf%aqed)
    type is (alpha_qed_from_scale_t)
       alpharef = alpha%ref
       alphamu = alpha%get (scale)
    type is (alpha_qed_fixed_t)
       call msg_fatal &
            ("t integrator: has to be called with running alpha.")
    end select
    t = - log (alphamu/alpharef) / two / Pi / coeffqed_b0 (0, epdf%n_lep)
  end function t_alpha

  function rec_series (p, al_2pi, expansion) result (rec)
    real(default) :: rec
    real(default), intent(in) :: p, al_2pi
    real(default), dimension(6) :: expansion
    rec = expansion(1) * p + expansion(2) * p**2 / two + &
         expansion(3) * p**3/6._default + al_2pi * (expansion(4) + &
         expansion(5) * p + expansion(6) * p**2 / two)
  end function rec_series

  function f_lim_1 (x) result (f)
    real(default), intent(in) :: x
    real(default) :: xb, f
    xb = one - x
    f = two * log(xb)**2 + (two/three)*Pi2 * log(xb) + &
         (two/three)*Pi2 - four * zeta3
  end function f_lim_1

  function f_lim_2 (x) result (f)
    real(default), intent(in) :: x
    real(default) :: xb, f
    xb = one - x
    f = -two*log(xb)**2 + (two - two/three*Pi2) * log(xb) + &
         (Pi2/three) + four*zeta3
  end function f_lim_2

  function func_1 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = ((2 + (-1 + x)*y*(2 + (-1 + x)*y)) * &
         ((1 + (-1 + x)*y)*log(1 - x)*  &
         ((-1 + x)*log(1 - x) + 2*x*log(x)) - &
         2*x*log(x/(1 + (-1 + x)*y))* &
         log(1 - x/(1 + (-1 + x)*y)) + &
         (-1 + x)*(-1 + y)* &
         log(1 - x/(1 + (-1 + x)*y))**2))/ &
         ((-1 + x)*x*y*(1 + (-1 + x)*y)) + &
         2*(1 + (-1 + x)*y)* &
         (log(1 - x)**2*log(x) - &
         log(x/(1 + (-1 + x)*y))* &
         log(one - x/(one + (-one + x)*y))**2) - &
         f_lim_1 (x)
  end function func_1

  function func_2 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = -(((2 + 2*y*(-1 + x) + y**2*(-1 + x)**2) * &
         ((-1 + x)*((1 + y*(-1 + x))*log(1 - x)**2 + &
         (-1 + y)*log(((-1 + y)*(-1 + x))/(1 + y*(-1 + x)))**2) &
         - (x + y*(-1 + x)*x)*polylog(2,x) + &
         x*polylog(2,x/(1 + y*(-1 + x)))))/(y*(1 + &
         y*(-1 + x))*(-1 + x)*x)) + 2*(1 + y*(-1 + x)) * &
         (log(1 - x)*polylog(2,x) - log(((-1 + y)*(-1 + x)) / &
         (1 + y*(-1 + x)))*polylog(2,x/(1 + y*(-1 + x)))) - &
         f_lim_2(x)
  end function func_2

  function func_3 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = ((2 + 2*y*(-1 + x) + y**2*(-1 + x)**2) * &
         ((1 + y*(-1 + x))*polylog(2,1 - x) - &
         polylog(2,((-1 + y)*(-1 + x))/(1 + y*(-1 + x))))) / &
         (y*(1 + y*(-1 + x))*(-1 + x)) + 2*(1 + y*(-1 + x)) * &
         (polylog(3,1 - x) - polylog(3,((-1 + y)*(-1 + x)) / &
         (1 + y*(-1 + x)))) + two
  end function func_3

  function func_4 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = two*(1 + y*(-1 + x))*(log(x)**2*log(1 + x) - &
         log(x/(1 + y*(-1 + x)))**2*log(1 + x/(1 + y*(-1 + x)))) - &
         ((1 + (1 + y*(-1 + x))**2) * (-(log(x)**2/(1 + x)) - &
         (2*log(x)*log(1 + x))/x - ((-1 + y)*log(x/(1 + y*(-1 + x))) * &
         (x*log(x/(1 + y*(-1 + x))) + 2*(1 + y*(-1 + x) + x) * &
         log(1 + x/(1 + y*(-1 + x))))) / ((1 + y*(-1 + x)) * x * &
         (1 + y*(-1 + x) + x))))/y
  end function func_4

  function func_5 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f =  -(((1 + (1 + y*(-1 + x))**2) * ((-2*log(x)*log(1 + x)) / &
         (1 + x) - log(1 + x)**2/x - (2*(-1 + y)*log(x/(1 + &
         y*(-1 + x)))*log(1 + x/(1 + y*(-1 + x)))) / ((1 + &
         y*(-1 + x))*(1 + y*(-1 + x) + x)) - ((-1 + y) * &
         log(1 + x/(1 + y*(-1 + x)))**2)/(x + y*(-1 + x)*x)))/y) + &
         2*(1 + y*(-1 + x))*(log(x)*log(1 + x)**2 - log(x/(1 + &
         y*(-1 + x))) * log(1 + x/(1 + y*(-1 + x)))**2) - &
         two*log(two)**2
  end function func_5

  function func_8 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = -(((2 + 2*y*(-1 + x) + y**2*(-1 + x)**2)*(log(x)*log(1 + x) - &
         y*log(x)*log(1 + x) + y*x*log(x)*log(1 + x) - log(x/(1 + &
         y*(-1 + x)))*log(1 + x/(1 + y*(-1 + x))) + y*log(x/(1 + &
         y*(-1 + x)))*log(1 + x/(1 + y*(-1 + x))) + (-1 + y - y*x) * &
         polylog(2,-x) - (-1 + y)*polylog(2,-(x/(1 + y*(-1 + x)))))) / &
         (y*(1 + y*(-1 + x))*x)) + 2*(1 + y*(-1 + x))*(log(x) * &
         polylog(2,-x) - log(x/(1 + y*(-1 + x))) * &
         polylog(2,-(x/(1 + y*(-1 + x))))) + Pi2/6._default
  end function func_8

  function func_9 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = -(((1 + (1 + y*(-1 + x))**2) * (log(1 + x)**2 / x + &
         ((-1 + y)*log(1 + x/(1 + y*(-1 + x)))**2) / (x + &
         y*(-1 + x)*x) - polylog(2,-x)/(1 + x) - ((-1 + y) * &
         polylog(2,-(x/(1 + y*(-1 + x))))) / ((1 + y*(-1 + x)) * &
         (1 + y*(-1 + x) + x))))/y) + 2*(1 + y*(-1 + x)) * &
         (log(1 + x)*polylog(2,-x) - log(1 + x/(1 + y*(-1 + x))) * &
         polylog(2,-(x/(1 + y*(-1 + x))))) + Pi2/12._default + &
         two*log(two)**2
  end function func_9

  function func_10 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = -(((1 + (1 + y*(-1 + x))**2) * (-((log(x/(1 + x)) * &
         log(1 + x))/(1 + x)) - ((-1 + y) * &
         log(x/(1 + y*(-1 + x) + x))*log(1 + x/(1 + y*(-1 + x)))) / &
         ((1 + y*(-1 + x))*(1 + y*(-1 + x) + x)) - &
         polylog(2,1/(1 + x))/(1 + x) - ((-1 + y) * &
         polylog(2,1/(1 + x/(1 + y*(-1 + x))))) / &
         ((1 + y*(-1 + x))*(1 + y*(-1 + x) + x))))/y) + &
         2*(1 + y*(-1 + x))*(log(1 + x)*polylog(2,one/(1 + x)) - &
         log(1 + x/(1 + y*(-1 + x))) * polylog(2,one/(1 + x/(1 + y*(-1 + x))))) - &
         Pi2/12._default + three/two*log(two)**2
  end function func_10

  function func_11 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = ((2 + 2*y*(-1 + x) + y**2*(-1 + x)**2)*((1 + &
         y*(-1 + x))*polylog(2,-x) + (-1 + y) * &
         polylog(2,-(x/(1 + y*(-1 + x)))))) / (y*(1 + &
         y*(-1 + x))*x) + 2*(1 + y*(-1 + x))*(polylog(3,-x) - &
         polylog(3,-(x/(1 + y*(-1 + x))))) + Pi2/6._default
  end function func_11

  function func_12 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = -(((1 + (1 + y*(-1 + x))**2)*(polylog(2,1/(1 + x))/(1 + x) + &
         ((-1 + y)*polylog(2,1/(1 + x/(1 + y*(-1 + x))))) / &
         ((1 + y*(-1 + x))*(1 + y*(-1 + x) + x))))/y) + 2*(1 + &
         y*(-1 + x))*(polylog(3,1/(1 + x)) - &
         polylog(3,1/(1 + x/(1 + y*(-1 + x))))) + Pi2/12._default - &
         0.5_default*log(two)**2
  end function func_12

  function func_13 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f =  -(((2 + 2*y*(-1 + x) + y**2*(-1 + x)**2)*log((1 - y) / &
         (1 + y*(-1 + x)))*log(1 + y*(-1 + x)))/(y*(-1 + x))) - &
         ((2 + 2*y*(-1 + x) + y**2*(-1 + x)**2)*log((1 - y) / &
         (1 + y*(-1 + x)))*log(y - y*x))/(1 + y*(-1 + x)) + &
         ((2 + 2*y*(-1 + x) + y**2*(-1 + x)**2)*log(1 + &
         y*(-1 + x))*log(y - y*x))/(1 + y*(-1 + x)) - 2*(1 + &
         y*(-1 + x))*log((1 - y)/(1 + y*(-1 + x)))*log(1 + &
         y*(-1 + x))*log(y - y*x) - two*log(1-x) + two - Pi2/three
  end function func_13

  function func_14 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f =  (log(1 + y*(-1 + x)) * (((2 + y*(-1 + x))*(2 + &
         2*y*(-1 + x) + y**2*(-1 + x)**2)*log(1 + &
         y*(-1 + x)))/(1 + y*(-1 + x)) + (2*y*(2 + &
         y*(-1 + x))*(2 + 2*y*(-1 + x) + y**2*(-1 + x)**2) * &
         (-1 + x) * log(((-1 + y)*(-1 + x))/(1 + &
         y*(-1 + x)))) / (1 + y*(-1 + x)) + (1 + (1 + &
         y*(-1 + x))**2)*(2 + y*(-1 + x)) * log(1 + y*(-1 + x)) * &
         log(((-1 + y)*(-1 + x))/(1 + y*(-1 + x))) + &
         y*(1 + (1 + y*(-1 + x))**2)*(1 - x)* log(1 + y*(-1 + x)) * &
         log(((-1 + y)*(-1 + x))/(1 + y*(-1 + x))) + &
         2*y*(1 + y*(-1 + x))*(2 + y*(-1 + x))*(-1 + x) * &
         log(1 + y*(-1 + x)) * log(((-1 + y)*(-1 + x)) / &
         (1 + y*(-1 + x))))) / (2 + y*(-1 + x))**2
  end function func_14

  function func_15 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = (((2 + y*(-1 + x))*(2 + 2*y*(-1 + x) + &
         y**2*(-1 + x)**2)*log(1 + y*(-1 + x))*log(2 + &
         y*(-1 + x)))/(1 + y*(-1 + x)) - y*(1 + (1 + &
         y*(-1 + x))**2)*(1 - x)*log(1 + y*(-1 + x)) * &
         log(((-1 + y)*(-1 + x))/(1 + y*(-1 + x))) + &
         (y*(2 + y*(-1 + x))*(2 + 2*y*(-1 + x) + &
         y**2*(-1 + x)**2)*(-1 + x)*log(2 + y*(-1 + x)) * &
         log(((-1 + y)*(-1 + x))/(1 + y*(-1 + x)))) / &
         (1 + y*(-1 + x)) + (1 + (1 + y*(-1 + x))**2) * &
         (2 + y*(-1 + x))*log(1 + y*(-1 + x))*log(2 + &
         y*(-1 + x))*log(((-1 + y)*(-1 + x))/(1 + y*(-1 + x))) + &
         y*(1 + (1 + y*(-1 + x))**2)*(1 - x)*log(1 + &
         y*(-1 + x))*log(2 + y*(-1 + x))*log(((-1 + y) * &
         (-1 + x))/(1 + y*(-1 + x))) + 2*y*(1 + y*(-1 + &
         x))*(2 + y*(-1 + x))*(-1 + x)*log(1 + y*(-1 + x)) * &
         log(2 + y*(-1 + x))*log(((-1 + y)*(-1 + x)) / &
         (1 + y*(-1 + x)))) / (2 + y*(-1 + x))**2
  end function func_15

  function func_16 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f =  (-((y*(2 + y*(-1 + x))*(2 + 2*y*(-1 + x) + y**2*(-1 + x)**2) * &
         (-1 + x)*log(2 + y*(-1 + x))*log(((-1 + y)*(-1 + x))/(1 + &
         y*(-1 + x))))/(1 + y*(-1 + x))) + ((2 + y*(-1 + x))*(2 + &
         2*y*(-1 + x) + y**2*(-1 + x)**2)*polylog(2,-1 + y - y*x)) / &
         (1 + y*(-1 + x)) + (1 + (1 + y*(-1 + x))**2)*(2 + y*(-1 + &
         x))*log(((-1 + y)*(-1 + x))/(1 + y*(-1 + x))) * &
         polylog(2,-1 + y - y*x) + y*(1 + (1 + y*(-1 + x))**2) * &
         (1 - x)*log(((-1 + y)*(-1 + x))/(1 + y*(-1 + x))) * &
         polylog(2,-1 + y - y*x) + 2*y*(1 + y*(-1 + x))*(2 + &
         y*(-1 + x))*(-1 + x)*log(((-1 + y)*(-1 + x))/(1 + &
         y*(-1 + x)))*polylog(2,-1 + y - y*x))/(2 + y*(-1 + x))**2 + &
         Pi2/12._default*log(xb)
  end function func_16

  function func_20 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = polylog(2,((-1 + y)*(-1 + x))/(1 + y*(-1 + x))) / (1 + y*(-1 + x)) + &
         polylog(3,((-1 + y)*(-1 + x))/(1 + y*(-1 + x)))
  end function func_20

  function func_24 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = (log(1 + x/(1 + y*(-1 + x)))*(-((-1 + y)*(-1 + x)*(1 + y*(-1 + x) + x) * &
         log(1 + x/(1 + y*(-1 + x)))) + x*log(x/(1 + y*(-1 + x)))*(2*(-1 + y + &
         x - y*x) + (1 + y*(-1 + x) + x)*log(1 + x/(1 + y*(-1 + x)))))) / &
         ((1 + y*(-1 + x))**2*x*(1 + y*(-1 + x) + x))
  end function func_24

  function func_25 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = (((-1 + y)*(-1 + x)*log(1 + x/(1 + y*(-1 + x)))**2)/x + ((-1 + y + x - &
         y*x + (1 + y*(-1 + x) + x)*log(1 + x/(1 + y*(-1 + x)))) * &
         polylog(2,-(x/(1 + y*(-1 + x))))) / (1 + y*(-1 + x) + x)) / &
         (1 + y*(-1 + x))**2 + Pi2/12._default * log(two)
  end function func_25

  function func_30 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = ((-1 + y + x - y*x)*log(x/(1 + y*(-1 + x) + x)) * log(1 + x/(1 + &
         y*(-1 + x))) + (-1 + y + x - y*x + (1 + y*(-1 + x) + x)*log(1 + &
         x/(1 + y*(-1 + x))))*polylog(2,1/(1 + x/(1 + y*(-1 + x))))) / &
         ((1 + y*(-1 + x))**2*(1 + y*(-1 + x) + x)) - Pi2/12._default*log(two) + &
         0.5*log(two)**3
  end function func_30

  function func_32 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = (polylog(2,((-1 + y)*(-1 + x))/(1 + y*(-1 + x))) + &
         polylog(3,((-1 + y)*(-1 + x))/(1 + y*(-1 + x)))) / &
         (1 + y*(-1 + x))**2
  end function func_32

  function func_34 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = ((-1 + y)*(-1 + x)*polylog(2,1/(1 + x/(1 + y*(-1 + x))))) / &
         ((1 + y*(-1 + x))*(1 + y*(-1 + x) + x)) + polylog(3,1/(1 + &
         x/(1 + y*(-1 + x)))) - one/6._default*log(two)**3 + &
         Pi2/12._default*log(two) - 7._default/8._default*zeta3
  end function func_34

  function func_35 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = ((-1 + y)*(-1 + x)*polylog(2,1/(1 + x/(1 + y*(-1 + x)))) + &
         (1 + y*(-1 + x) + x)*polylog(3,1/(1 + x/(1 + y*(-1 + x))))) / &
         ((1 + y*(-1 + x))**2*(1 + y*(-1 + x) + x)) - &
         one/6._default*log(two)**3 + Pi2/12._default*log(two) - &
         7._default/8._default*zeta3
  end function func_35

  function func_37 (y, x) result (f)
    real(default), intent(in) :: x, y
    real(default) :: xb, f
    xb = one - x
    f = (log(y - y*x)*(log(y - y*x) + log(((-1 + y)*(-1 + x)) / &
         (1 + y*(-1 + x)))*(2 + 2*y*(-1 + x) + log(y - y*x)))) / &
         (1 + y*(-1 + x))**2 - log(xb)**3 + Pi2/three*log(xb) - two*zeta3
  end function func_37

  module function endpoint_func_S (x, nlep) result (f)
    real(default), intent(in) :: x
    integer, intent(in) :: nlep
    real(default) :: f
    real(default) :: result, abserr
    real(default), parameter :: epsabs = 0.001_default, &
         epsrel = 0.001_default
    real(default), parameter :: a = 0._default, b = 1._default
    integer, parameter :: limit = 10000
    call gauss_kronrod (GAUSS_KRONROD_41, int_fun, a, b, &
         limit, result, abserr, epsabs, epsrel)
    f = result
  contains
    function int_fun (y) result (d_f)
      real(default) :: d_f
      real(default), intent(in) :: y
      d_f = four*func_1(y,x) + four*func_2(y,x) + four*func_3(y,x) - &
           two*func_4(y,x) - four*func_5(y,x) - four*func_8(y,x) - &
           four*func_9(y,x) + four*func_10(y,x) + 4.*func_11(y,x) - &
           8._default*func_12(y,x) - four*func_13(y,x) + two*func_14(y,x) - &
           8._default*func_15(y,x) - 8._default*func_16(y,x) - &
           24._default*real(nlep,kind=default)*func_20(y,x);
    end function int_fun
  end function endpoint_func_S

  module function endpoint_func_NS (x) result (f)
    real(default), intent(in) :: x
    real(default) :: f
    real(default) :: result, abserr
    real(default), parameter :: epsabs = 0.001_default, &
         epsrel = 0.001_default
    real(default), parameter :: a = 0._default, b = 1._default
    integer, parameter :: limit = 10000
    call gauss_kronrod (GAUSS_KRONROD_41, int_fun, a, b, &
         limit, result, abserr, epsabs, epsrel)
    f = result
  contains
    function int_fun (y) result (d_f)
      real(default) :: d_f
      real(default), intent(in) :: y
      d_f = four*func_1(y,x) + four*func_2(y,x) + four*func_3(y,x) + &
           two*func_4(y,x) + four*func_5(y,x) + four*func_8(y,x) + &
           four*func_9(y,x) - four*func_10(y,x) - four*func_11(y,x) + &
           8._default*func_12(y,x) - four*func_13(y,x) - two*func_14(y,x) + &
           8._default*func_15(y,x) + 8._default*func_16(y,x)
    end function int_fun
  end function endpoint_func_NS

  module function endpoint_func_GAM (x) result (f)
    real(default), intent(in) :: x
    real(default) :: f
    real(default) :: result, abserr
    real(default), parameter :: epsabs = 0.001_default, &
         epsrel = 0.001_default
    real(default), parameter :: a = 0._default, b = 1._default
    integer, parameter :: limit = 10000
    call gauss_kronrod (GAUSS_KRONROD_41, int_fun, a, b, &
         limit, result, abserr, epsabs, epsrel)
    f = result
  contains
    function int_fun (y) result (d_f)
      real(default) :: d_f
      real(default), intent(in) :: y
      d_f = -8._default*func_24(y,x) - 8._default*func_25(y,x) + &
           8._default*func_30(y,x) + 8._default*func_32(y,x) + 16._default* &
           func_34(y,x) - 16._default*func_35(y,x) - four*func_37(y,x)
    end function int_fun
  end function endpoint_func_GAM

  function photon_matching (x, x0, x1, p) result (p_match)
    real(default) :: p_match
    real(default), intent(in) :: x, x0, x1, p
    real(default) :: xm, logx
    logx = - log10(1-x)
    if (logx < x0) then
       p_match = zero
    else if (logx > x1) then
       p_match = one
    else
       xm = (logx - x0) / (x1 - x0)
       p_match = xm**p / (xm**p + (1-xm)**p)
    end if
  end function photon_matching
  module function elec_pdf (epdf, flv, x, scale, alpha, &
       running, w_num) result (e_pdf)
    type(qed_pdf_t), intent(in) :: epdf
    integer, intent(in) :: flv
    real(default) :: e_pdf
    real(default), intent(in) :: x
    real(default), intent(in) :: scale
    real(default), intent(in) :: alpha    
    logical, intent(in) :: running, w_num
    integer :: nlep
    real(default) :: ln0, eta0, al_2pi, p
    real(default), parameter :: x0gam = 2.0_default, x1gam = &
         6._default, pgam = 2._default
    if (allocated (epdf%q_in)) then
       ln0 = log(epdf%q_in**2/epdf%mass**2)
       eta0 = alpha/Pi * log(scale**2/epdf%q_in**2) 
    else
       ln0 = zero
       eta0 = alpha/Pi * log(scale**2/epdf%mass**2)
    end if
    if (running) then
       p = t_alpha (epdf, scale)
    else
       p = eta0 / two
    end if
    select case (epdf%log_order)
    case (EPDF_LL)
       nlep = epdf%n_lep
       al_2pi = zero
    case (EPDF_NLL)
       if (running) then
          select type (aqed => epdf%aqed)
          type is (alpha_qed_from_scale_t)
             nlep = aqed%nlep
             al_2pi = aqed%get (scale) / two / Pi
          type is (alpha_qed_fixed_t)
             call msg_fatal &
                  ("elec_pdf: has to be called with running alpha.")
          end select
       else
          nlep = epdf%n_lep
          al_2pi = alpha / two / Pi
       end if
    end select
    select case (flv)
    case (EPDF_S)
       e_pdf = elec_asym (epdf, x, scale, alpha, running) + &
            recbar (epdf, flv, x, scale, alpha, running) + &
            rechat (epdf, flv, x, scale, alpha, running) - &
            bar_asym (epdf, flv, x, scale, alpha, running)
       if (w_num) then
          e_pdf = e_pdf + &
               rec_num (epdf, flv, x, scale, alpha, running)
       end if
    case (EPDF_NS)
       e_pdf = elec_asym (epdf, x, scale, alpha, running) + &
            recbar (epdf, flv, x, scale, alpha, running) + &
            rechat (epdf, flv, x, scale, alpha, running) - &
            bar_asym (epdf, flv, x, scale, alpha, running)
       if (w_num) then
          e_pdf = e_pdf + &
               rec_num (epdf, flv, x, scale, alpha, running)
       end if
    case (EPDF_G)
       e_pdf = recbar (epdf, flv, x, scale, alpha, running) + &
            rechat (epdf, flv, x, scale, alpha, running) + &
            photon_matching (x, x0gam, x1gam, pgam) * &
            (phot_asym (epdf, x, scale, alpha, nlep, running) - &
            recbar (epdf, flv, x, scale, alpha, running))
       if (w_num) then
          e_pdf = e_pdf + &
               rec_num (epdf, flv, x, scale, alpha, running)
       end if
    case (EPDF_ELE)
       e_pdf = elec_asym (epdf, x, scale, alpha, running) - &
            bar_asym (epdf, EPDF_S, x, scale, alpha, running) + &
            0.5_default * (recbar (epdf, EPDF_S, x, scale, alpha, running) + &
            recbar (epdf, EPDF_NS, x, scale, alpha, running) + &
            rechat (epdf, EPDF_S, x, scale, alpha, running) + &
            rechat (epdf, EPDF_NS, x, scale, alpha, running))
       if (w_num) then
          e_pdf = e_pdf + 0.5_default * ( &
               rec_num (epdf, EPDF_S, x, scale, alpha, running) + &
               rec_num (epdf, EPDF_NS, x, scale, alpha, running))
       end if
    case (EPDF_POS)
       e_pdf = 0.5_default * (recbar (epdf, EPDF_S, x, scale, alpha, running) - &
            recbar (epdf, EPDF_NS, x, scale, alpha, running) + &
            rechat (epdf, EPDF_S, x, scale, alpha, running) - &
            rechat (epdf, EPDF_NS, x, scale, alpha, running))
       if (w_num) then
          e_pdf = e_pdf + 0.5_default * ( &
            rec_num (epdf, EPDF_S, x, scale, alpha, running) - &
            rec_num (epdf, EPDF_NS, x, scale, alpha, running))
       end if
    case default
       call msg_fatal &
            ("elec_pdf: wrong lepton flavor.")
    end select
  end function elec_pdf


end submodule electron_pdfs_s
