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

submodule (sm_physics) sm_physics_s

  use io_units
  use numeric_utils
  use diagnostics
  use permutations, only: factorial

  implicit none

contains

  pure module function beta0 (nf)
    real(default), intent(in) :: nf
    real(default) :: beta0
    beta0 = 11.0_default - two/three * nf
  end function beta0

  pure module function beta1 (nf)
    real(default), intent(in) :: nf
    real(default) :: beta1
    beta1 = 51.0_default - 19.0_default/three * nf
  end function beta1

  pure module function beta2 (nf)
    real(default), intent(in) :: nf
    real(default) :: beta2
    beta2 = 2857.0_default - 5033.0_default / 9.0_default * &
                    nf + 325.0_default/27.0_default * nf**2
  end function beta2

  pure module function coeff_b0 (nf)
    real(default), intent(in) :: nf
    real(default) :: coeff_b0
    coeff_b0 = (11.0_default * CA - two * nf) / (12.0_default * pi)
  end function coeff_b0

  pure module function coeff_b1 (nf)
    real(default), intent(in) :: nf
    real(default) :: coeff_b1
    coeff_b1 = (17.0_default * CA**2 - five * CA * nf - three * CF * nf) / &
               (24.0_default * pi**2)
  end function coeff_b1

  pure module function coeff_b2 (nf)
    real(default), intent(in) :: nf
    real(default) :: coeff_b2
    coeff_b2 = (2857.0_default/54.0_default * CA**3 - &
                    1415.0_default/54.0_default * &
                    CA**2 * nf - 205.0_default/18.0_default * CA*CF*nf &
                    + 79.0_default/54.0_default * CA*nf**2 + &
                    11.0_default/9.0_default * CF * nf**2) / (four*pi)**3
  end function coeff_b2

  pure module function coeffqed_b0 (nf, nlep)
    integer, intent(in) :: nf, nlep
    real(default) :: n_lep, coeffqed_b0
    n_lep = real(nlep, kind=default)
    coeffqed_b0 = - (three * sumQ2q (nf) + n_lep) / (three*pi)
  end function coeffqed_b0

  pure module function coeffqed_b1 (nf, nlep)
    integer, intent(in) :: nf, nlep
    real(default) :: n_lep, coeffqed_b1
    n_lep = real(nlep, kind=default)
    coeffqed_b1 = - (three * sumQ4q (nf) + n_lep) / (four*pi**2)
  end function coeffqed_b1

  pure function sumQ2q (nf)
    integer, intent(in) :: nf
    real(default) :: sumQ2q
    select case (nf)
    case (0)
       sumQ2q = zero
    case (1)
       sumQ2q = 1.0_default/9.0_default
    case (2)
       sumQ2q = 5.0_default/9.0_default
    case (3)
       sumQ2q = 2.0_default/3.0_default
    case (4)
       sumQ2q = 10.0_default/9.0_default
    case (5)
       sumQ2q = 11.0_default/9.0_default
    case (6:)
       sumQ2q = 5.0_default/3.0_default
    end select
  end function sumQ2q

  pure function sumQ4q (nf)
    integer, intent(in) :: nf
    real(default) :: sumQ4q
    select case (nf)
    case (0)
       sumQ4q = zero
    case (1)
       sumQ4q = 1.0_default/81.0_default
    case (2)
       sumQ4q = 17.0_default/81.0_default
    case (3)
       sumQ4q = 2.0_default/9.0_default
    case (4)
       sumQ4q = 34.0_default/81.0_default
    case (5)
       sumQ4q = 35.0_default/81.0_default
    case (6:)
       sumQ4q = 17.0_default/27.0_default
    end select
  end function sumQ4q

  pure module function running_as (scale, al_mz, mz, order, nf) result (ascale)
    real(default), intent(in) :: scale
    real(default), intent(in), optional :: al_mz, nf, mz
    integer, intent(in), optional :: order
    integer :: ord
    real(default) :: az, m_z, as_log, n_f, b0, b1, b2, ascale
    real(default) :: as0, as1
    if (present (mz)) then
       m_z = mz
    else
       m_z = MZ_REF
    end if
    if (present (order)) then
       ord = order
    else
       ord = 0
    end if
    if (present (al_mz)) then
       az = al_mz
    else
       az = ALPHA_QCD_MZ_REF
    end if
    if (present (nf)) then
       n_f = nf
    else
       n_f = 5
    end if
    b0 = coeff_b0 (n_f)
    b1 = coeff_b1 (n_f)
    b2 = coeff_b2 (n_f)
    as_log = one + b0 * az * log(scale**2/m_z**2)
    as0 = az / as_log
    as1 = as0 - as0**2 * b1/b0 * log(as_log)
    select case (ord)
    case (0)
       ascale = as0
    case (1)
       ascale = as1
    case (2)
       ascale = as1 + as0**3 * (b1**2/b0**2 * ((log(as_log))**2 - &
            log(as_log) + as_log - one) - b2/b0 * (as_log - one))
    case default
       ascale = as0
    end select
  end function running_as

  pure module function running_as_lam (nf, scale, lambda, order) result (ascale)
    real(default), intent(in) :: nf, scale
    real(default), intent(in), optional :: lambda
    integer, intent(in), optional :: order
    real(default) :: lambda_qcd
    real(default) :: as0, as1, logmul, b0, b1, b2, ascale
    integer :: ord
    if (present (lambda)) then
       lambda_qcd = lambda
    else
       lambda_qcd = LAMBDA_QCD_REF
    end if
    if (present (order)) then
       ord = order
    else
       ord = 0
    end if
    b0 = beta0(nf)
    logmul = log(scale**2/lambda_qcd**2)
    as0 = four*pi / b0 / logmul
    if (ord > 0) then
       b1 = beta1(nf)
       as1 = as0 * (one - two* b1 / b0**2 * log(logmul) / logmul)
    end if
    select case (ord)
    case (0)
       ascale = as0
    case (1)
       ascale = as1
    case (2)
       b2 = beta2(nf)
       ascale = as1 + as0 * four * b1**2/b0**4/logmul**2 * &
            ((log(logmul) - 0.5_default)**2 + &
             b2*b0/8.0_default/b1**2 - five/four)
    case default
       ascale = as0
    end select
  end function running_as_lam

  pure module function running_alpha &
       (scale, al_me, me, order, nf, nlep) result (ascale)
    real(default), intent(in) :: scale
    real(default), intent(in), optional :: al_me, me
    integer, intent(in), optional :: order, nf, nlep
    integer :: ord, n_f, n_lep
    real(default) :: ae, m_e, a_log, b0, b1, ascale
    real(default) :: a0, a1
    if (present (me)) then
       m_e = me
    else
       m_e = ME_REF
    end if
    if (present (order)) then
       ord = order
    else
       ord = 0
    end if
    if (present (al_me)) then
       ae = al_me
    else
       ae = ALPHA_QED_ME_REF
    end if
    if (present (nf)) then
       n_f = nf
    else
       n_f = 5
    end if
    if (present (nlep)) then
       n_lep = nlep
    else
       n_lep = 1
    end if
    b0 = coeffqed_b0 (n_f, n_lep)
    b1 = coeffqed_b1 (n_f, n_lep)
    a_log = one + b0 * ae * log(scale**2/m_e**2)
    a0 = ae / a_log
    a1 = ae / (a_log + ae * b1/b0 * &
         log((a_log + ae * b1/b0)/(one + ae * b1/b0)))
    select case (ord)
    case (0)
       ascale = a0
    case (1)
       ascale = a1
    case default
       ascale = a0
    end select
  end function running_alpha

  pure module function running_alpha_num &
       (scale, al_me, me, order, nf, nlep) result (ascale)
    real(default), intent(in) :: scale
    real(default), intent(in), optional :: al_me, me
    integer, intent(in), optional :: order, nf, nlep
    integer, parameter :: n_steps = 20
    integer :: ord, n_f, n_lep, k1
    real(default), parameter :: sxth = 1._default/6._default
    real(default) :: ae, ascale, m_e, log_q, dlr, &
         b0, b1, xk0, xk1, xk2, xk3
    if (present (order)) then
       ord = order
    else
       ord = 0
    end if
    if (present (al_me)) then
       ae = al_me
    else
       ae = ALPHA_QED_ME_REF
    end if
    if (present (me)) then
       m_e = me
    else
       m_e = ME_REF
    end if
    if (present (nf)) then
       n_f = nf
    else
       n_f = 5
    end if
    if (present (nlep)) then
       n_lep = nlep
    else
       n_lep = 1
    end if
    ascale = ae
    log_q = log (scale**2/m_e**2)
    dlr = log_q / n_steps
    b0 = coeffqed_b0 (n_f, n_lep)
    b1 = coeffqed_b1 (n_f, n_lep)
    ! ..Solution of the evolution equation depending on ORD
    !   (fourth-order Runge-Kutta beyond the leading order)
    select case (ord)
    case (0)
       ascale = ae / (one + b0 * ae * log_q)
    case (1:)
       do k1 = 1, n_steps
          xk0 = dlr * beta_qed (ascale)
          xk1 = dlr * beta_qed (ascale + 0.5 * xk0)
          xk2 = dlr * beta_qed (ascale + 0.5 * xk1)
          xk3 = dlr * beta_qed (ascale + xk2)
          ascale = ascale + sxth * (xk0 + 2._default * xk1 + &
               2._default * xk2 + xk3)
       end do
    end select
  contains
    pure function beta_qed (alpha)
      real(default), intent(in) :: alpha
      real(default) :: beta_qed
      beta_qed = - alpha**2 * (b0 + alpha * b1)
    end function beta_qed
  end function running_alpha_num

  module function lambda_qcd (as_q, q, nf, order) result (lambda)
    real(default), intent(in) :: as_q, q
    integer, intent(in) :: nf, order
    real(default) :: lambda
    real(default), parameter :: acc = 1e-8_default
    if (order == 0) then
       lambda = lambda_qcd_lo (as_q, q, nf)
    else if (order == 1) then
       lambda = lambda_qcd_nlo (as_q, q, nf)
    else if (order == 2) then
       lambda = lambda_qcd_nnlo (as_q, q, nf)
    else
       call msg_error ("lambda_qcd: order unknown")
    end if
  contains
    function lambda_qcd_lo (as_q, q, nf) result (lambda)
      real(default) :: lambda
      real(default), intent(in) :: as_q, q
      integer, intent(in) :: nf
      real(default) :: b0, t0, t1, as0, as1
      b0 = coeff_b0(real(nf, default))
      t1 = one/b0/as_q
      FIND_ROOT: do
         if (signal_is_pending ())  return
         t0 = t1
         as0 = one/b0/t1
         as1 = - one/b0/t1**2
         t1 = (as_q-as0)/as1 + t1
         if (abs(t0-t1)/t0 < acc) exit FIND_ROOT
      end do FIND_ROOT
      lambda = q * exp(-t1/two)
    end function lambda_qcd_lo
    function lambda_qcd_nlo (as_q, q, nf) result (lambda)
      real(default) :: lambda
      real(default), intent(in) :: as_q, q
      integer, intent(in) :: nf
      real(default) :: b0, b1, t0, t1, as0, as1, logt
      b0 = coeff_b0(real(nf, default))
      b1 = coeff_b1(real(nf, default))
      t1 = one/b0/as_q
      FIND_ROOT: do
         if (signal_is_pending ())  return
         logt = log(t1)
         t0 = t1
         as0 = one/b0/t1 - b1/b0 * logt/(b0 * t1)**2
         as1 = - one/b0/t1**2 - b1/b0**3 * (one - two*logt)/t1**3
         t1 = (as_q-as0)/as1 + t1
         if (abs(t0-t1)/t0 < acc) exit FIND_ROOT
      end do FIND_ROOT
      lambda = q * exp(-t1/two)
    end function lambda_qcd_nlo
    function lambda_qcd_nnlo (as_q, q, nf) result (lambda)
      real(default) :: lambda
      real(default), intent(in) :: as_q, q
      integer, intent(in) :: nf
      real(default) :: b0, b1, b2, t0, t1, as0, as1, logt
      b0 = coeff_b0(real(nf, default))
      b1 = coeff_b1(real(nf, default))
      b2 = coeff_b2(real(nf, default))
      t1 = one/b0/as_q
      FIND_ROOT: do
         if (signal_is_pending ())  return
         logt = log(t1)
         t0 = t1
         as0 = one/b0/t1 * (one - b1/b0**2 * logt/t1 + (b1/b0**2 * logt/t1)**2 &
              - (b1**2 * (logt + one) - b0*b2)/b0**4/t1**2)
         as1 = one/b0/t1 * (-two*b1**2 * logt**2/(b0**4 * t1**3) &
              + two*(b1**2 * (logt + one) - b0*b2)/(b0**4 * t1**3) &
              + b1 * logt/(b0**2 * t1**2) + two*b1**2 * logt/(b0**4 * t1**3) &
              - b1/(b0**2 * t1**2) - b1**2/(b0**4 * t1**3)) &
              - (b1**2 * logt**2/(b0**4 * t1**2) - (b1**2 * (logt + one) &
              - b0*b2)/(b0**4 * t1**2) - b1 * logt/(b0**2 * t1) + one)/(b0 * t1**2)
         t1 = (as_q-as0)/as1 + t1
         if (abs(t0-t1)/t0 < acc) exit FIND_ROOT
      end do FIND_ROOT
      lambda = q * exp(-t1/two)
    end function lambda_qcd_nnlo
  end function lambda_qcd

  elemental module function gamma_g (nf) result (gg)
    real(default), intent(in) :: nf
    real(default) :: gg
    gg = 11.0_default/6.0_default * CA - two/three * TR * nf
  end function gamma_g

  elemental module function k_g (nf) result (kg)
    real(default), intent(in) :: nf
    real(default) :: kg
    kg = (67.0_default/18.0_default - pi**2/6.0_default) * CA - &
         10.0_default/9.0_default * TR * nf
  end function k_g

  elemental module function Li2 (x)
    real(default), intent(in) :: x
    real(default) :: Li2
    Li2 = real( Li2_double (real(x, kind=double)), kind=default)
  end function Li2

  elemental function Li2_double (x)  result (Li2)
    real(double), intent(in) :: x
    real(double) :: Li2
    real(double), parameter :: pi2_6 = pi**2/6
    if (abs(1-x) < tiny_07) then
       Li2 = pi2_6
    else if (abs(1-x) <  0.5_double) then
       Li2 = pi2_6 - log(1-x) * log(x) - Li2_restricted (1-x)
    else if (abs(x) > 1.d0) then
       ! Li2 = 0
       ! call msg_bug (" Dilogarithm called outside of defined range.")
       !!! Reactivate Dilogarithm identity
        Li2 = -pi2_6 - 0.5_default * log(-x) * log(-x) - Li2_restricted (1/x)
    else
       Li2 = Li2_restricted (x)
    end if
  contains
    elemental function Li2_restricted (x) result (Li2)
      real(double), intent(in) :: x
      real(double) :: Li2
      real(double) :: tmp, z, z2
      z = - log (1-x)
      z2 = z**2
! Horner's rule for the powers z^3 through z^19
      tmp = 43867._double/798._double
      tmp = tmp * z2 /342._double - 3617._double/510._double
      tmp = tmp * z2 /272._double + 7._double/6._double
      tmp = tmp * z2 /210._double - 691._double/2730._double
      tmp = tmp * z2 /156._double + 5._double/66._double
      tmp = tmp * z2 /110._double - 1._double/30._double
      tmp = tmp * z2 / 72._double + 1._double/42._double
      tmp = tmp * z2 / 42._double - 1._double/30._double
      tmp = tmp * z2 / 20._double + 1._double/6._double
! The first three terms of the power series
      Li2 = z2 * z * tmp / 6._double - 0.25_double * z2 + z
    end function Li2_restricted
  end function Li2_double

  elemental module function psic (z) result (psi)
    complex(default), intent(in) :: z
    complex(default) :: psi
    complex(default) :: shift, zz, zi, zi2
    shift = 0
    zz = z
    if (abs (aimag(zz)) < 10._default) then
       do while (abs (zz) < 10._default)
          shift = shift - 1 / zz
          zz = zz + 1
       end do
    end if
    zi  = 1/zz
    zi2 = zi*zi
    psi = shift + log(zz) - zi/2 - zi2 / 5040._default * ( 420._default + &
         zi2 * ( -42._default + zi2 * (20._default - 21._default * zi2)))
  end function psic

  elemental module function psir (x) result (psi)
    real(default), intent(in) :: x
    real(default) :: psi
    psi = real (psic (cmplx (x,0,kind=default)), kind=default)
  end function psir

  elemental module function psim (z, m) result (psi)
    complex(default), intent(in) :: z
    integer, intent(in) :: m
    complex(default) :: psi
    complex(default) :: shift, rec, zz, zi, zi2
    real(default) :: c1, c2, c3, c4, c5, c6, c7
    integer :: i
    if (m < 1) then
       psi = psic(z)
    else
       shift = 0
       zz = z
       if (abs (aimag (zz)) < 10._default) then
          CHECK_ABS: do
             rec = (-1)**m * factorial (m) / zz**(m+1)
             shift = shift - rec
             zz = zz + 1
             if (abs (zz) > 10._default) exit CHECK_ABS
          end do CHECK_ABS
       end if
       c1 =   1._default
       c2 =   1._default / 2._default
       c3 =   1._default / 6._default
       c4 = - 1._default / 30._default
       c5 =   1._default / 42._default
       c6 = - 1._default / 30._default
       c7 =   5._default / 66._default
       do i = 2, m
          c1 = c1 * (i-1)
          c2 = c2 * i
          c3 = c3 * (i+1)
          c4 = c4 * (i+3)
          c5 = c5 * (i+5)
          c6 = c6 * (i+7)
          c7 = c7 * (i+9)
       end do
       zi  = 1/zz
       zi2 = zi*zi
       psi = shift + (-1)**(m-1) * zi**m * ( c1 + zi * ( c2 + zi * ( &
            c3 + zi2 * ( c4 + zi2 * ( c5 + zi2 * ( c6 + ( c7 * zi2)))))))
    end if
  end function psim

  elemental module function psimr (x, m) result (psi)
    real(default), intent(in) :: x
    integer, intent(in) :: m
    real(default) :: psi
    psi = real (psim (cmplx (x,0,kind=default), m), kind=default)
  end function psimr

  module function cnielsen (n, m, x) result (nplog)
    integer, intent(in) :: n, m
    real(default), intent(in) :: x
    complex(default) :: nplog
    real(default), parameter :: c1 = 4._default/3._default, &
         c2 = 1._default/3._default
    real(default), dimension(0:4), parameter :: &
         fct = [1.0_default,1.0_default,2.0_default,6.0_default,24.0_default]
    real(default), dimension(4,4) :: s1, cc
    real(default), dimension(0:30,10) :: aa
    complex(default), dimension(0:5) :: vv
    real(default), dimension(0:5) :: uu
    real(default) :: x1, h, alfa, b0, b1, b2, qq, rr
    complex(default) :: sj, sk
    integer, dimension(10), parameter :: &
         nc = [24,26,28,30,22,24,26,19,22,17]
    integer, dimension(31), parameter :: &
         index = [1,2,3,4,0,0,0,0,0,0,5,6,7,0,0,0,0,0,0,0, &
                  8,9,0,0,0,0,0,0,0,0,10]
    real(default), dimension(0:4), parameter :: &
         sgn = [1._default, -1._default, 1._default, -1._default, 1._default]
    integer :: it, j, k, l, m1, n1
    if ((n<1) .or. (n>4) .or. (m<1) .or. (m>4) .or. (n+m > 5)) then
       call msg_fatal &
            ("The Nielsen dilogarithms cannot be applied for these values.")
    end if
    s1 = 0._default
    s1(1,1) = 1.6449340668482_default
    s1(1,2) = 1.2020569031596_default
    s1(1,3) = 1.0823232337111_default
    s1(1,4) = 1.0369277551434_default
    s1(2,1) = 1.2020569031596_default
    s1(2,2) = 2.7058080842778e-1_default
    s1(2,3) = 9.6551159989444e-2_default
    s1(3,1) = 1.0823232337111_default
    s1(3,2) = 9.6551159989444e-2_default
    s1(4,1) = 1.0369277551434_default
    cc = 0._default
    cc(1,1) = 1.6449340668482_default
    cc(1,2) = 1.2020569031596_default
    cc(1,3) = 1.0823232337111_default
    cc(1,4) = 1.0369277551434_default
    cc(2,2) =-1.8940656589945_default
    cc(2,3) =-3.0142321054407_default
    cc(3,1) = 1.8940656589945_default
    cc(3,2) = 3.0142321054407_default
    aa = 0._default
    aa( 0,1) = 0.96753215043498_default
    aa( 1,1) = 0.16607303292785_default
    aa( 2,1) = 0.02487932292423_default
    aa( 3,1) = 0.00468636195945_default
    aa( 4,1) = 0.00100162749616_default
    aa( 5,1) = 0.00023200219609_default
    aa( 6,1) = 0.00005681782272_default
    aa( 7,1) = 0.00001449630056_default
    aa( 8,1) = 0.00000381632946_default
    aa( 9,1) = 0.00000102990426_default
    aa(10,1) = 0.00000028357538_default
    aa(11,1) = 0.00000007938705_default
    aa(12,1) = 0.00000002253670_default
    aa(13,1) = 0.00000000647434_default
    aa(14,1) = 0.00000000187912_default
    aa(15,1) = 0.00000000055029_default
    aa(16,1) = 0.00000000016242_default
    aa(17,1) = 0.00000000004827_default
    aa(18,1) = 0.00000000001444_default
    aa(19,1) = 0.00000000000434_default
    aa(20,1) = 0.00000000000131_default
    aa(21,1) = 0.00000000000040_default
    aa(22,1) = 0.00000000000012_default
    aa(23,1) = 0.00000000000004_default
    aa(24,1) = 0.00000000000001_default

    aa( 0,2) = 0.95180889127832_default
    aa( 1,2) = 0.43131131846532_default
    aa( 2,2) = 0.10002250714905_default
    aa( 3,2) = 0.02442415595220_default
    aa( 4,2) = 0.00622512463724_default
    aa( 5,2) = 0.00164078831235_default
    aa( 6,2) = 0.00044407920265_default
    aa( 7,2) = 0.00012277494168_default
    aa( 8,2) = 0.00003453981284_default
    aa( 9,2) = 0.00000985869565_default
    aa(10,2) = 0.00000284856995_default
    aa(11,2) = 0.00000083170847_default
    aa(12,2) = 0.00000024503950_default
    aa(13,2) = 0.00000007276496_default
    aa(14,2) = 0.00000002175802_default
    aa(15,2) = 0.00000000654616_default
    aa(16,2) = 0.00000000198033_default
    aa(17,2) = 0.00000000060204_default
    aa(18,2) = 0.00000000018385_default
    aa(19,2) = 0.00000000005637_default
    aa(20,2) = 0.00000000001735_default
    aa(21,2) = 0.00000000000536_default
    aa(22,2) = 0.00000000000166_default
    aa(23,2) = 0.00000000000052_default
    aa(24,2) = 0.00000000000016_default
    aa(25,2) = 0.00000000000005_default
    aa(26,2) = 0.00000000000002_default

    aa( 0,3) = 0.98161027991365_default
    aa( 1,3) = 0.72926806320726_default
    aa( 2,3) = 0.22774714909321_default
    aa( 3,3) = 0.06809083296197_default
    aa( 4,3) = 0.02013701183064_default
    aa( 5,3) = 0.00595478480197_default
    aa( 6,3) = 0.00176769013959_default
    aa( 7,3) = 0.00052748218502_default
    aa( 8,3) = 0.00015827461460_default
    aa( 9,3) = 0.00004774922076_default
    aa(10,3) = 0.00001447920408_default
    aa(11,3) = 0.00000441154886_default
    aa(12,3) = 0.00000135003870_default
    aa(13,3) = 0.00000041481779_default
    aa(14,3) = 0.00000012793307_default
    aa(15,3) = 0.00000003959070_default
    aa(16,3) = 0.00000001229055_default
    aa(17,3) = 0.00000000382658_default
    aa(18,3) = 0.00000000119459_default
    aa(19,3) = 0.00000000037386_default
    aa(20,3) = 0.00000000011727_default
    aa(21,3) = 0.00000000003687_default
    aa(22,3) = 0.00000000001161_default
    aa(23,3) = 0.00000000000366_default
    aa(24,3) = 0.00000000000116_default
    aa(25,3) = 0.00000000000037_default
    aa(26,3) = 0.00000000000012_default
    aa(27,3) = 0.00000000000004_default
    aa(28,3) = 0.00000000000001_default

    aa( 0,4) = 1.0640521184614_default
    aa( 1,4) = 1.0691720744981_default
    aa( 2,4) = 0.41527193251768_default
    aa( 3,4) = 0.14610332936222_default
    aa( 4,4) = 0.04904732648784_default
    aa( 5,4) = 0.01606340860396_default
    aa( 6,4) = 0.00518889350790_default
    aa( 7,4) = 0.00166298717324_default
    aa( 8,4) = 0.00053058279969_default
    aa( 9,4) = 0.00016887029251_default
    aa(10,4) = 0.00005368328059_default
    aa(11,4) = 0.00001705923313_default
    aa(12,4) = 0.00000542174374_default
    aa(13,4) = 0.00000172394082_default
    aa(14,4) = 0.00000054853275_default
    aa(15,4) = 0.00000017467795_default
    aa(16,4) = 0.00000005567550_default
    aa(17,4) = 0.00000001776234_default
    aa(18,4) = 0.00000000567224_default
    aa(19,4) = 0.00000000181313_default
    aa(20,4) = 0.00000000058012_default
    aa(21,4) = 0.00000000018579_default
    aa(22,4) = 0.00000000005955_default
    aa(23,4) = 0.00000000001911_default
    aa(24,4) = 0.00000000000614_default
    aa(25,4) = 0.00000000000197_default
    aa(26,4) = 0.00000000000063_default
    aa(27,4) = 0.00000000000020_default
    aa(28,4) = 0.00000000000007_default
    aa(29,4) = 0.00000000000002_default
    aa(30,4) = 0.00000000000001_default

    aa( 0,5) = 0.97920860669175_default
    aa( 1,5) = 0.08518813148683_default
    aa( 2,5) = 0.00855985222013_default
    aa( 3,5) = 0.00121177214413_default
    aa( 4,5) = 0.00020722768531_default
    aa( 5,5) = 0.00003996958691_default
    aa( 6,5) = 0.00000838064065_default
    aa( 7,5) = 0.00000186848945_default
    aa( 8,5) = 0.00000043666087_default
    aa( 9,5) = 0.00000010591733_default
    aa(10,5) = 0.00000002647892_default
    aa(11,5) = 0.00000000678700_default
    aa(12,5) = 0.00000000177654_default
    aa(13,5) = 0.00000000047342_default
    aa(14,5) = 0.00000000012812_default
    aa(15,5) = 0.00000000003514_default
    aa(16,5) = 0.00000000000975_default
    aa(17,5) = 0.00000000000274_default
    aa(18,5) = 0.00000000000077_default
    aa(19,5) = 0.00000000000022_default
    aa(20,5) = 0.00000000000006_default
    aa(21,5) = 0.00000000000002_default
    aa(22,5) = 0.00000000000001_default

    aa( 0,6) = 0.95021851963952_default
    aa( 1,6) = 0.29052529161433_default
    aa( 2,6) = 0.05081774061716_default
    aa( 3,6) = 0.00995543767280_default
    aa( 4,6) = 0.00211733895031_default
    aa( 5,6) = 0.00047859470550_default
    aa( 6,6) = 0.00011334321308_default
    aa( 7,6) = 0.00002784733104_default
    aa( 8,6) = 0.00000704788108_default
    aa( 9,6) = 0.00000182788740_default
    aa(10,6) = 0.00000048387492_default
    aa(11,6) = 0.00000013033842_default
    aa(12,6) = 0.00000003563769_default
    aa(13,6) = 0.00000000987174_default
    aa(14,6) = 0.00000000276586_default
    aa(15,6) = 0.00000000078279_default
    aa(16,6) = 0.00000000022354_default
    aa(17,6) = 0.00000000006435_default
    aa(18,6) = 0.00000000001866_default
    aa(19,6) = 0.00000000000545_default
    aa(20,6) = 0.00000000000160_default
    aa(21,6) = 0.00000000000047_default
    aa(22,6) = 0.00000000000014_default
    aa(23,6) = 0.00000000000004_default
    aa(24,6) = 0.00000000000001_default

    aa( 0,7) = 0.95064032186777_default
    aa( 1,7) = 0.54138285465171_default
    aa( 2,7) = 0.13649979590321_default
    aa( 3,7) = 0.03417942328207_default
    aa( 4,7) = 0.00869027883583_default
    aa( 5,7) = 0.00225284084155_default
    aa( 6,7) = 0.00059516089806_default
    aa( 7,7) = 0.00015995617766_default
    aa( 8,7) = 0.00004365213096_default
    aa( 9,7) = 0.00001207474688_default
    aa(10,7) = 0.00000338018176_default
    aa(11,7) = 0.00000095632476_default
    aa(12,7) = 0.00000027313129_default
    aa(13,7) = 0.00000007866968_default
    aa(14,7) = 0.00000002283195_default
    aa(15,7) = 0.00000000667205_default
    aa(16,7) = 0.00000000196191_default
    aa(17,7) = 0.00000000058018_default
    aa(18,7) = 0.00000000017246_default
    aa(19,7) = 0.00000000005151_default
    aa(20,7) = 0.00000000001545_default
    aa(21,7) = 0.00000000000465_default
    aa(22,7) = 0.00000000000141_default
    aa(23,7) = 0.00000000000043_default
    aa(24,7) = 0.00000000000013_default
    aa(25,7) = 0.00000000000004_default
    aa(26,7) = 0.00000000000001_default

    aa( 0,8) = 0.98800011672229_default
    aa( 1,8) = 0.04364067609601_default
    aa( 2,8) = 0.00295091178278_default
    aa( 3,8) = 0.00031477809720_default
    aa( 4,8) = 0.00004314846029_default
    aa( 5,8) = 0.00000693818230_default
    aa( 6,8) = 0.00000124640350_default
    aa( 7,8) = 0.00000024293628_default
    aa( 8,8) = 0.00000005040827_default
    aa( 9,8) = 0.00000001099075_default
    aa(10,8) = 0.00000000249467_default
    aa(11,8) = 0.00000000058540_default
    aa(12,8) = 0.00000000014127_default
    aa(13,8) = 0.00000000003492_default
    aa(14,8) = 0.00000000000881_default
    aa(15,8) = 0.00000000000226_default
    aa(16,8) = 0.00000000000059_default
    aa(17,8) = 0.00000000000016_default
    aa(18,8) = 0.00000000000004_default
    aa(19,8) = 0.00000000000001_default

    aa( 0,9) = 0.95768506546350_default
    aa( 1,9) = 0.19725249679534_default
    aa( 2,9) = 0.02603370313918_default
    aa( 3,9) = 0.00409382168261_default
    aa( 4,9) = 0.00072681707110_default
    aa( 5,9) = 0.00014091879261_default
    aa( 6,9) = 0.00002920458914_default
    aa( 7,9) = 0.00000637631144_default
    aa( 8,9) = 0.00000145167850_default
    aa( 9,9) = 0.00000034205281_default
    aa(10,9) = 0.00000008294302_default
    aa(11,9) = 0.00000002060784_default
    aa(12,9) = 0.00000000522823_default
    aa(13,9) = 0.00000000135066_default
    aa(14,9) = 0.00000000035451_default
    aa(15,9) = 0.00000000009436_default
    aa(16,9) = 0.00000000002543_default
    aa(17,9) = 0.00000000000693_default
    aa(18,9) = 0.00000000000191_default
    aa(19,9) = 0.00000000000053_default
    aa(20,9) = 0.00000000000015_default
    aa(21,9) = 0.00000000000004_default
    aa(22,9) = 0.00000000000001_default

    aa( 0,10) = 0.99343651671347_default
    aa( 1,10) = 0.02225770126826_default
    aa( 2,10) = 0.00101475574703_default
    aa( 3,10) = 0.00008175156250_default
    aa( 4,10) = 0.00000899973547_default
    aa( 5,10) = 0.00000120823987_default
    aa( 6,10) = 0.00000018616913_default
    aa( 7,10) = 0.00000003174723_default
    aa( 8,10) = 0.00000000585215_default
    aa( 9,10) = 0.00000000114739_default
    aa(10,10) = 0.00000000023652_default
    aa(11,10) = 0.00000000005082_default
    aa(12,10) = 0.00000000001131_default
    aa(13,10) = 0.00000000000259_default
    aa(14,10) = 0.00000000000061_default
    aa(15,10) = 0.00000000000015_default
    aa(16,10) = 0.00000000000004_default
    aa(17,10) = 0.00000000000001_default

    if (x == 1._default) then
       nplog = s1(n,m)
    else if (x > 2._default .or. x < -1.0_default) then
       x1 = 1._default / x
       h = c1 * x1 + c2
       alfa = h + h
       vv(0) = 1._default
       if (x < -1.0_default) then
          vv(1) = log(-x)
       else if (x > 2._default) then
          vv(1) = log(cmplx(-x,0._default,kind=default))
       end if
       do l = 2, n+m
          vv(l) = vv(1) * vv(l-1)/l
       end do
       sk = 0._default
       do k = 0, m-1
          m1 = m-k
          rr = x1**m1 / (fct(m1) * fct(n-1))
          sj = 0._default
          do j = 0, k
             n1 = n+k-j
             l = index(10*n1+m1-10)
             b1 = 0._default
             b2 = 0._default
             do it = nc(l), 0, -1
                b0 = aa(it,l) + alfa*b1 - b2
                b2 = b1
                b1 = b0
             end do
             qq = (fct(n1-1) / fct(k-j)) * (b0 - h*b2) * rr / m1**n1
             sj = sj + vv(j) * qq
          end do
          sk = sk + sgn(k) * sj
       end do
       sj = 0._default
       do j = 0, n-1
          sj = sj + vv(j) * cc(n-j,m)
       end do
       nplog = sgn(n) * sk + sgn(m) * (sj + vv(n+m))
    else if (x > 0.5_default) then
       x1 = 1._default - x
       h = c1 * x1 + c2
       alfa = h + h
       vv(0) = 1._default
       uu(0) = 1._default
       vv(1) = log(cmplx(x1,0._default,kind=default))
       uu(1) = log(x)
       do l = 2, m
          vv(l) = vv(1) * vv(l-1) / l
       end do
       do l = 2, n
          uu(l) = uu(1) * uu(l-1) / l
       end do
       sk = 0._default
       do k = 0, n-1
          m1 = n-k
          rr = x1**m1 / fct(m1)
          sj = 0._default
          do j = 0, m-1
             n1 = m-j
             l = index(10*n1 + m1 - 10)
             b1 = 0._default
             b2 = 0._default
             do it = nc(l), 0, -1
                b0 = aa(it,l) + alfa*b1 - b2
                b2 = b1
                b1 = b0
             end do
             qq = sgn(j) * (b0 - h*b2) * rr / m1**n1
             sj = sj + vv(j) * qq
          end do
          sk = sk + uu(k) * (s1(m1,m) - sj)
       end do
       nplog = sk + sgn(m) * uu(n) * vv(m)
    else
       l = index(10*n + m - 10)
       h = c1 * x + c2
       alfa = h + h
       b1 = 0._default
       b2 = 0._default
       do it = nc(l), 0, -1
          b0 = aa(it,l) + alfa*b1 - b2
          b2 = b1
          b1 = b0
       end do
       nplog = (b0 - h*b2) * x**m / (fct(m) * m**n)
    end if
  end function cnielsen

  module function nielsen (n, m, x) result (nplog)
    integer, intent(in) :: n, m
    real(default), intent(in) :: x
    real(default) :: nplog
    nplog = real (cnielsen (n, m, x))
  end function nielsen

  module function polylog (n, x) result (plog)
    integer, intent(in) :: n
    real(default), intent(in) :: x
    real(default) :: plog
    plog = nielsen (n-1,1,x)
  end function polylog

  module function dilog (x) result (dlog)
    real(default), intent(in) :: x
    real(default) :: dlog
    dlog = polylog (2,x)
  end function dilog

  module function trilog (x) result (tlog)
    real(default), intent(in) :: x
    real(default) :: tlog
    tlog = polylog (3,x)
  end function trilog

  elemental module function faux (x) result (y)
    real(default), intent(in) :: x
    complex(default) :: y
    if (1 <= x) then
       y = asin(sqrt(1/x))**2
    else
       y = - 1/4.0_default * (log((1 + sqrt(1 - x))/ &
            (1 - sqrt(1 - x))) - cmplx (0.0_default, pi, kind=default))**2
    end if
  end function faux

  elemental module function fonehalf (x) result (y)
    real(default), intent(in) :: x
    complex(default) :: y
    if (abs(x) < eps0) then
       y = 0
    else
       y = - 2.0_default * x * (1 + (1 - x) * faux(x))
    end if
  end function fonehalf

  module function fonehalf_pseudo (x) result (y)
    real(default), intent(in) :: x
    complex(default) :: y
    if (abs(x) < eps0) then
       y = 0
    else
       y = - 2.0_default * x * faux(x)
    end if
  end function fonehalf_pseudo

  elemental module function fone (x) result  (y)
    real(default), intent(in) :: x
    complex(default) :: y
    if (abs(x) < eps0) then
       y = 2.0_default
    else
       y = 2.0_default + 3.0_default * x + &
            3.0_default * x * (2.0_default - x) * &
            faux(x)
    end if
  end function fone

  elemental module function gaux (x) result (y)
    real(default), intent(in) :: x
    complex(default) :: y
    if (1 <= x) then
       y = sqrt(x - 1) * asin(sqrt(1/x))
    else
       y = sqrt(1 - x) * (log((1 + sqrt(1 - x)) / &
            (1 - sqrt(1 - x))) - &
            cmplx (0.0_default, pi, kind=default)) / 2.0_default
    end if
  end function gaux

  elemental module function tri_i1 (a,b) result (y)
    real(default), intent(in) :: a,b
    complex(default) :: y
    if (a < eps0 .or. b < eps0) then
       y = 0
    else
       y = a*b/2.0_default/(a-b) + a**2 * b**2/2.0_default/(a-b)**2 * &
            (faux(a) - faux(b)) + &
            a**2 * b/(a-b)**2 * (gaux(a) - gaux(b))
    end if
  end function tri_i1

  elemental module function tri_i2 (a,b) result (y)
    real(default), intent(in) :: a,b
    complex(default) :: y
    if (a < eps0 .or. b < eps0) then
       y = 0
    else
       y = - a * b / 2.0_default / (a-b) * (faux(a) - faux(b))
    end if
  end function tri_i2

  elemental module function run_b0 (nf) result (bnull)
    integer, intent(in) :: nf
    real(default) :: bnull
    bnull = 33.0_default - 2.0_default * nf
  end function run_b0

  elemental module function run_b1 (nf) result (bone)
    integer, intent(in) :: nf
    real(default) :: bone
    bone = 6.0_default * (153.0_default - 19.0_default * nf)/run_b0(nf)**2
  end function run_b1

  elemental module function run_aa (nf) result (aaa)
    integer, intent(in) :: nf
    real(default) :: aaa
    aaa = 12.0_default * PI / run_b0(nf)
  end function run_aa

  elemental function run_bb (nf) result (bbb)
    integer, intent(in) :: nf
    real(default) :: bbb
    bbb = run_b1(nf) / run_aa(nf)
  end function run_bb

  pure module subroutine ff_dipole (v_ijk, y_ijk, p_ij, pp_k, p_i, p_j, p_k)
    type(vector4_t), intent(in) :: p_i, p_j, p_k
    type(vector4_t), intent(out) :: p_ij, pp_k
    real(default), intent(out) :: y_ijk
    real(default) :: z_i
    real(default), intent(out) :: v_ijk
    z_i   = (p_i*p_k) / ((p_k*p_j) + (p_k*p_i))
    y_ijk = (p_i*p_j) / ((p_i*p_j) + (p_i*p_k) + (p_j*p_k))
    p_ij  = p_i + p_j - y_ijk/(1.0_default - y_ijk) * p_k
    pp_k  = (1.0/(1.0_default - y_ijk)) * p_k
    !!! We don't multiply by alpha_s right here:
    v_ijk = 8.0_default * PI * CF * &
         (2.0 / (1.0 - z_i*(1.0 - y_ijk)) - (1.0 + z_i))
  end subroutine ff_dipole

  pure module subroutine fi_dipole (v_ija, x_ija, p_ij, pp_a, p_i, p_j, p_a)
    type(vector4_t), intent(in) :: p_i, p_j, p_a
    type(vector4_t), intent(out) :: p_ij, pp_a
    real(default), intent(out) :: x_ija
    real(default) :: z_i
    real(default), intent(out) :: v_ija
    z_i   = (p_i*p_a) / ((p_a*p_j) + (p_a*p_i))
    x_ija = ((p_i*p_a) + (p_j*p_a) - (p_i*p_j)) &
         / ((p_i*p_a) + (p_j*p_a))
    p_ij  = p_i + p_j - (1.0_default - x_ija) * p_a
    pp_a  = x_ija * p_a
    !!! We don't not multiply by alpha_s right here:
    v_ija = 8.0_default * PI * CF * &
         (2.0 / (1.0 - z_i + (1.0 - x_ija)) - (1.0 + z_i)) / x_ija
  end subroutine fi_dipole

  pure module subroutine if_dipole (v_kja, u_j, p_aj, pp_k, p_k, p_j, p_a)
    type(vector4_t), intent(in) :: p_k, p_j, p_a
    type(vector4_t), intent(out) :: p_aj, pp_k
    real(default), intent(out) :: u_j
    real(default) :: x_kja
    real(default), intent(out) :: v_kja
    u_j   = (p_a*p_j) / ((p_a*p_j) + (p_a*p_k))
    x_kja = ((p_a*p_k) + (p_a*p_j) - (p_j*p_k)) &
         / ((p_a*p_j) + (p_a*p_k))
    p_aj  = x_kja * p_a
    pp_k  = p_k + p_j - (1.0_default - x_kja) * p_a
    v_kja = 8.0_default * PI * CF * &
         (2.0 / (1.0 - x_kja + u_j) - (1.0 + x_kja)) / x_kja
  end subroutine if_dipole

  pure module subroutine ii_dipole (v_jab, v_j, p_in, p_out, flag_1or2)
    type(vector4_t), dimension(:), intent(in) :: p_in
    type(vector4_t), dimension(size(p_in)-1), intent(out) :: p_out
    logical, intent(in) :: flag_1or2
    real(default), intent(out) :: v_j
    real(default), intent(out) :: v_jab
    type(vector4_t) :: p_a, p_b, p_j
    type(vector4_t) :: k, kk
    type(vector4_t) :: p_aj
    real(default) :: x_jab
    integer :: i
    !!! flag_1or2 decides whether this a 12 or 21 dipole
    if (flag_1or2) then
       p_a = p_in(1)
       p_b = p_in(2)
    else
       p_b = p_in(1)
       p_a = p_in(2)
    end if
    !!! We assume that the unresolved particle has always the last
    !!! momentum
    p_j = p_in(size(p_in))
    x_jab = ((p_a*p_b) - (p_a*p_j) - (p_b*p_j)) / (p_a*p_b)
    v_j = (p_a*p_j) / (p_a * p_b)
    p_aj  = x_jab * p_a
    k     = p_a + p_b - p_j
    kk    = p_aj + p_b
    do i = 3, size(p_in)-1
       p_out(i) = p_in(i) - 2.0*((k+kk)*p_in(i))/((k+kk)*(k+kk)) * (k+kk) + &
            (2.0 * (k*p_in(i)) / (k*k)) * kk
    end do
    if (flag_1or2) then
       p_out(1) = p_aj
       p_out(2) = p_b
    else
       p_out(1) = p_b
       p_out(2) = p_aj
    end if
    v_jab = 8.0_default * PI * CF * &
         (2.0 / (1.0 - x_jab) - (1.0 + x_jab)) / x_jab
  end subroutine ii_dipole
  elemental module function delta (x,eps) result (z)
     real(default), intent(in) :: x, eps
     real(default) :: z
     if (x > one - eps) then
        z = one / eps
     else
        z = 0
     end if
  end function delta

  elemental module function plus_distr (x,eps) result (plusd)
    real(default), intent(in) :: x, eps
    real(default) :: plusd
    if (x > one - eps) then
       plusd = log(eps) / eps
    else
       plusd = one / (one - x)
    end if
  end function plus_distr

  elemental module function pqq (x,eps) result (pqqx)
    real(default), intent(in) :: x, eps
    real(default) :: pqqx
    if (x > (1.0_default - eps)) then
       pqqx = (eps - one) / two + two * log(eps) / eps - &
            three * (eps - one) / eps / two
    else
       pqqx = (one + x**2) / (one - x)
    end if
    pqqx = CF * pqqx
  end function pqq

  elemental module function pgq (x) result (pgqx)
    real(default), intent(in) :: x
    real(default) :: pgqx
    pgqx = TR * (x**2 + (one - x)**2)
  end function pgq

  elemental module function pqg (x) result (pqgx)
    real(default), intent(in) :: x
    real(default) :: pqgx
    pqgx = CF * (one + (one - x)**2) / x
  end function pqg

  elemental module function pgg (x, nf, eps) result (pggx)
    real(default), intent(in) :: x, nf, eps
    real(default) :: pggx
    pggx = two * CA * ( plus_distr (x, eps) + (one-x)/x - one + &
                   x*(one-x)) + delta (x, eps)  * gamma_g(nf)
  end function pgg

  elemental module function pqq_reg (x) result (pqqregx)
     real(default), intent(in) :: x
     real(default) :: pqqregx
     pqqregx = - CF * (one + x)
  end function pqq_reg

  elemental module function pgg_reg (x) result (pggregx)
     real(default), intent(in) :: x
     real(default) :: pggregx
     pggregx = two * CA * ((one - x)/x - one + x*(one - x))
  end function pgg_reg

  module function kbarqg (x) result (kbarqgx)
    real(default), intent(in) :: x
    real(default) :: kbarqgx
    kbarqgx = pqg(x) * log((one-x)/x) + CF * x
  end function kbarqg

  module function kbargq (x) result (kbargqx)
    real(default), intent(in) :: x
    real(default) :: kbargqx
    kbargqx = pgq(x) * log((one-x)/x) + two * TR * x * (one - x)
  end function kbargq

  module function kbarqq (x,eps) result (kbarqqx)
    real(default), intent(in) :: x, eps
    real(default) :: kbarqqx
    kbarqqx = CF*(log_plus_distr(x,eps) - (one+x) * log((one-x)/x) + (one - &
         x) - (five - pi**2) * delta(x,eps))
  end function kbarqq

  module function kbargg (x,eps,nf) result (kbarggx)
    real(default), intent(in) :: x, eps, nf
    real(default) :: kbarggx
    kbarggx = CA * (log_plus_distr(x,eps) + two * ((one-x)/x - one + &
                         x*(one-x) * log((1-x)/x))) - delta(x,eps) * &
                         ((50.0_default/9.0_default - pi**2) * CA - &
                         16.0_default/9.0_default * TR * nf)
  end function kbargg

  module function ktildeqq (x,eps) result (ktildeqqx)
    real(default), intent(in) :: x, eps
    real(default) :: ktildeqqx
    ktildeqqx = pqq_reg (x) * log(one-x) + CF * ( - log2_plus_distr (x,eps) &
                          - pi**2/three * delta(x,eps))
  end function ktildeqq

  module function ktildeqg (x,eps) result (ktildeqgx)
    real(default), intent(in) :: x, eps
    real(default) :: ktildeqgx
    ktildeqgx = pqg (x) * log(one-x)
  end function ktildeqg

  module function ktildegq (x,eps) result (ktildegqx)
    real(default), intent(in) :: x, eps
    real(default) :: ktildegqx
    ktildegqx = pgq (x) * log(one-x)
  end function ktildegq

  module function ktildegg (x,eps) result (ktildeggx)
    real(default), intent(in) :: x, eps
    real(default) :: ktildeggx
    ktildeggx = pgg_reg (x) * log(one-x) + CA * ( - &
       log2_plus_distr (x,eps) - pi**2/three * delta(x,eps))
  end function ktildegg

  pure module function insert_q () result (i_q)
    real(default), dimension(0:2) :: i_q
    i_q(0) = gamma_q + k_q - pi**2/three * CF
    i_q(1) = gamma_q
    i_q(2) = CF
  end function insert_q

  pure module function insert_g (nf) result (i_g)
    real(default), intent(in) :: nf
    real(default), dimension(0:2) :: i_g
    i_g(0) = gamma_g (nf) + k_g (nf) - pi**2/three * CA
    i_g(1) = gamma_g (nf)
    i_g(2) = CA
  end function insert_g

  pure module function k_q_al (alpha)
    real(default), intent(in) :: alpha
    real(default) :: k_q_al
    k_q_al = k_q - CF * (log(alpha))**2 + gamma_q * &
                      (alpha - one - log(alpha))
  end function k_q_al

  pure module function k_g_al (alpha, nf)
    real(default), intent(in) :: alpha, nf
    real(default) :: k_g_al
    k_g_al = k_g (nf) - CA * (log(alpha))**2 + gamma_g (nf) * &
                     (alpha - one - log(alpha))
  end function k_g_al

  module function plus_distr_al (x,alpha,eps) result (plusd_al)
    real(default), intent(in) :: x,  eps, alpha
    real(default) :: plusd_al
    if ((one - alpha) >= (one - eps)) then
       plusd_al = zero
       call msg_fatal ('sm_physics, plus_distr_al: alpha and epsilon chosen wrongly')
    elseif (x < (1.0_default - alpha)) then
       plusd_al = 0
    else if (x > (1.0_default - eps)) then
       plusd_al = log(eps/alpha)/eps
    else
       plusd_al = one/(one-x)
    end if
  end function plus_distr_al

  module function kbarqg_al (x,alpha,eps) result (kbarqgx)
    real(default), intent(in) :: x, alpha, eps
    real(default) :: kbarqgx
    kbarqgx = pqg (x) * log(alpha*(one-x)/x) + CF * x
  end function kbarqg_al
  module function kbargq_al (x,alpha,eps) result (kbargqx)
    real(default), intent(in) :: x, alpha, eps
    real(default) :: kbargqx
    kbargqx = pgq (x) * log(alpha*(one-x)/x) + two * TR * x * (one-x)
  end function kbargq_al
  module function kbarqq_al (x,alpha,eps) result (kbarqqx)
    real(default), intent(in) :: x, alpha, eps
    real(default) :: kbarqqx
    kbarqqx = CF * (one - x) + pqq_reg(x) * log(alpha*(one-x)/x) &
         + CF * log_plus_distr(x,eps) &
         - (gamma_q + k_q_al(alpha) - CF * &
         five/6.0_default  * pi**2 - CF * (log(alpha))**2) * &
         delta(x,eps) + &
         CF * two/(one -x)*log(alpha*(two-x)/(one+alpha-x))
    if (x < (one-alpha)) then
       kbarqqx = kbarqqx - CF * two/(one-x) * log((two-x)/(one-x))
    end if
  end function kbarqq_al

  module function kbargg_al (x,alpha,eps,nf) result (kbarggx)
    real(default), intent(in) :: x, alpha, eps, nf
    real(default) :: kbarggx
    kbarggx = pgg_reg(x) * log(alpha*(one-x)/x) &
         + CA * log_plus_distr(x,eps) &
         - (gamma_g(nf) + k_g_al(alpha,nf) - CA * &
         five/6.0_default  * pi**2 - CA * (log(alpha))**2) * &
         delta(x,eps) + &
         CA * two/(one -x)*log(alpha*(two-x)/(one+alpha-x))
    if (x < (one-alpha)) then
       kbarggx = kbarggx - CA * two/(one-x) * log((two-x)/(one-x))
    end if
  end function kbargg_al

  module function ktildeqq_al (x,alpha,eps) result (ktildeqqx)
    real(default), intent(in) :: x, eps, alpha
    real(default) :: ktildeqqx
    ktildeqqx = pqq_reg(x) * log((one-x)/alpha) + CF*( &
         - log2_plus_distr_al(x,alpha,eps) - Pi**2/three * delta(x,eps) &
         + (one+x**2)/(one-x) * log(min(one,(alpha/(one-x)))) &
         + two/(one-x) * log((one+alpha-x)/alpha))
    if (x > (one-alpha)) then
       ktildeqqx = ktildeqqx - CF*two/(one-x)*log(two-x)
    end if
  end function ktildeqq_al

  module function log_plus_distr (x,eps) result (lpd)
     real(default), intent(in) :: x, eps
     real(default) :: lpd, eps2
     eps2 = min (eps, 0.1816_default)
     if (x > (1.0_default - eps2)) then
        lpd = ((log(eps2))**2 + two*Li2(eps2) - pi**2/three)/eps2
     else
        lpd = two*log((one-x)/x)/(one-x)
     end if
  end function log_plus_distr

  module function log2_plus_distr (x,eps) result (lpd)
    real(default), intent(in) :: x, eps
    real(default) :: lpd
    if (x > (1.0_default - eps)) then
       lpd = - (log(eps))**2/eps
    else
       lpd = two*log(one/(one-x))/(one-x)
    end if
  end function log2_plus_distr

  module function log2_plus_distr_al (x,alpha,eps) result (lpd_al)
    real(default), intent(in) :: x, eps, alpha
    real(default) :: lpd_al
    if ((one - alpha) >= (one - eps)) then
       lpd_al = zero
       call msg_fatal ('alpha and epsilon chosen wrongly')
    elseif (x < (one - alpha)) then
       lpd_al = 0
    elseif (x > (1.0_default - eps)) then
       lpd_al = - ((log(eps))**2 - (log(alpha))**2)/eps
    else
       lpd_al = two*log(one/(one-x))/(one-x)
    end if
  end function log2_plus_distr_al

  elemental module function p_qqg (z) result (P)
    real(default), intent(in) :: z
    real(default) :: P
    P = CF * (one + z**2) / (one - z)
  end function p_qqg
  elemental module function p_gqq (z) result (P)
    real(default), intent(in) :: z
    real(default) :: P
    P = TR * (z**2 + (one - z)**2)
  end function p_gqq
  elemental module function p_ggg (z) result (P)
    real(default), intent(in) :: z
    real(default) :: P
    P = NC * ((one - z) / z + z / (one - z) + z * (one - z))
  end function p_ggg

  pure module function integral_over_p_qqg (zmin, zmax) result (integral)
    real(default), intent(in) :: zmin, zmax
    real(default) :: integral
    integral = (two / three) * (- zmax**2 + zmin**2 - &
         two * (zmax - zmin) + four * log((one - zmin) / (one - zmax)))
  end function integral_over_p_qqg

  pure module function integral_over_p_gqq (zmin, zmax) result (integral)
    real(default), intent(in) :: zmin, zmax
    real(default) :: integral
    integral = 0.5_default * ((two / three) * &
         (zmax**3 - zmin**3) - (zmax**2 - zmin**2) + (zmax - zmin))
  end function integral_over_p_gqq

  pure module function integral_over_p_ggg (zmin, zmax) result (integral)
    real(default), intent(in) :: zmin, zmax
    real(default) :: integral
    integral = three * ((log(zmax) - two * zmax - &
         log(one - zmax) + zmax**2 / two - zmax**3 / three) - &
         (log(zmin) - zmin - zmin - log(one - zmin) + zmin**2 &
         / two - zmin**3 / three) )
  end function integral_over_p_ggg

  elemental module function p_qqg_pol (z, l_a, l_b, l_c) result (P)
    real(default), intent(in) :: z
    integer, intent(in) :: l_a, l_b, l_c
    real(default) :: P
    if (l_a /= l_b) then
       P = zero
       return
    end if
    if (l_c == -1) then
       P = one - z
    else
       P = (one + z)**2 / (one - z)
    end if
    P = P * CF
  end function p_qqg_pol

  module function pqqm (n, c_f) result (pqq_m)
    integer, intent(in) :: n
    real(default), intent(in) :: c_f
    complex(default) :: pqq_m
    pqq_m = three - four * (eulerc + &
         psic(cmplx(N+1,zero,kind=default))) + two/N/(N+1)
  end function pqqm

  elemental module function top_width_sm_lo (alpha, sinthw, vtb, mtop, mw, mb) &
         result (gamma)
    real(default) :: gamma
    real(default), intent(in) :: alpha, sinthw, vtb, mtop, mw, mb
    real(default) :: kappa
    kappa = sqrt ((mtop**2 - (mw + mb)**2) * (mtop**2 - (mw - mb)**2))
    gamma = alpha / four * mtop / (two * sinthw**2) * &
         vtb**2 * kappa / mtop**2 * &
         ((mtop**2 + mb**2) / (two * mtop**2) + &
          (mtop**2 - mb**2)**2 / (two * mtop**2 * mw**2) - &
           mw**2 / mtop**2)
  end function top_width_sm_lo

  elemental module function g_mu_from_alpha (alpha, mw, sinthw) result (g_mu)
    real(default) :: g_mu
    real(default), intent(in) :: alpha, mw, sinthw
    g_mu = pi * alpha / sqrt(two) / mw**2 / sinthw**2
  end function g_mu_from_alpha

  elemental module function alpha_from_g_mu (g_mu, mw, sinthw) result (alpha)
    real(default) :: alpha
    real(default), intent(in) :: g_mu, mw, sinthw
    alpha = g_mu * sqrt(two) / pi * mw**2 * sinthw**2
  end function alpha_from_g_mu

  elemental module function top_width_sm_qcd_nlo_massless_b &
         (alpha, sinthw, vtb, mtop, mw, alphas) result (gamma)
    real(default) :: gamma
    real(default), intent(in) :: alpha, sinthw, vtb, mtop, mw, alphas
    real(default) :: prefac, g_mu, w2
    g_mu = g_mu_from_alpha (alpha, mw, sinthw)
    prefac = g_mu * mtop**3 * vtb**2 / (16 * sqrt(two) * pi)
    w2 = mw**2 / mtop**2
    gamma = prefac * (f0 (w2) - (two * alphas) / (3 * Pi) * f1 (w2))
  end function top_width_sm_qcd_nlo_massless_b

  elemental module function f0 (w2) result (f)
    real(default) :: f
    real(default), intent(in) :: w2
    f = two * (one - w2)**2 * (1 + 2 * w2)
  end function f0

  elemental module function f1 (w2) result (f)
    real(default) :: f
    real(default), intent(in) :: w2
    f = f0 (w2) * (pi**2 + two * Li2 (w2) - two * Li2 (one - w2)) &
         + four * w2 * (one - w2 - two * w2**2) * log (w2) &
         + two * (one - w2)**2 * (five + four * w2) * log (one - w2) &
         - (one - w2) * (five + 9 * w2 - 6 * w2**2)
  end function f1

  elemental module function top_width_sm_qcd_nlo_jk &
         (alpha, sinthw, vtb, mtop, mw, mb, alphas) result (gamma)
    real(default) :: gamma
    real(default), intent(in) :: alpha, sinthw, vtb, mtop, mw, mb, alphas
    real(default) :: prefac, g_mu, eps2, i_xi
    g_mu = g_mu_from_alpha (alpha, mw, sinthw)
    prefac = g_mu * mtop**3 * vtb**2 / (16 * sqrt(two) * pi)
    eps2 = (mb / mtop)**2
    i_xi = (mw / mtop)**2
    gamma = prefac * (ff0 (eps2, i_xi) - &
         (two * alphas) / (3 * Pi) * ff1 (eps2, i_xi))
  end function top_width_sm_qcd_nlo_jk

  elemental module function top_width_sm_qcd_nlo_ce &
       (alpha, sinthw, vtb, mtop, mw, mb, alpha_s) result (gamma)
    real(default) :: gamma
    real(default), intent(in) :: alpha, sinthw, vtb, mtop, mw, mb, alpha_s
    real(default) :: pm, pp, p0, p3
    real(default) :: yw, yp
    real(default) :: W0, Wp, Wm, w2
    real(default) :: beta2
    real(default) :: f
    real(default) :: g_mu, gamma0
    beta2 = (mb / mtop)**2
    w2 = (mw / mtop)**2
    p0 = (one - w2 + beta2) / two
    p3 = sqrt (lambda (one, w2, beta2)) / two
    pp = p0 + p3
    pm = p0 - p3
    W0 = (one + w2 - beta2) / two
    Wp = W0 + p3
    Wm = W0 - p3
    yp = log (pp / pm) / two
    yw = log (Wp / Wm) / two
    f = (one - beta2)**2 + w2 * (one + beta2) - two * w2**2
    g_mu = g_mu_from_alpha (alpha, mw, sinthw)
    gamma0 = g_mu * mtop**3 * vtb**2 / (8 * pi * sqrt(two))
    gamma = gamma0 * alpha_s / twopi * CF * &
         (8 * f * p0 * (Li2(one - pm) - Li2(one - pp) - two * Li2(one - pm / pp) &
         + yp * log((four * p3**2) / (pp**2 * Wp)) + yw * log (pp)) &
         + four * (one - beta2) * ((one - beta2)**2 + w2 * (one + beta2) - four * w2**2) * yw &
         + (3 - beta2 + 11 * beta2**2 - beta2**3 + w2 * (6 - 12 * beta2 + two * beta2**2) &
         - w2**2 * (21 + 5 * beta2) + 12 * w2**3) * yp &
         + 8 * f * p3 * log (sqrt(w2) / (four * p3**2)) &
         + 6 * (one - four * beta2 + 3 * beta2**2 + w2 * (3 + beta2) - four * w2**2) * p3 * log(sqrt(beta2)) &
         + (5 - 22 * beta2 + 5 * beta2**2 + 9 * w2 * (one + beta2) - 6 * w2**2) * p3)
  end function top_width_sm_qcd_nlo_ce

  elemental module function ff0 (eps2, w2) result (f)
    real(default) :: f
    real(default), intent(in) :: eps2, w2
    f = one / two * sqrt(ff_lambda (eps2, w2)) * ff_f0 (eps2, w2)
  end function ff0

  elemental module function ff_f0 (eps2, w2) result (f)
    real(default) :: f
    real(default), intent(in) :: eps2, w2
    f = four * ((1 - eps2)**2 + w2 * (1 + eps2) - 2 * w2**2)
  end function ff_f0

  elemental module function ff_lambda (eps2, w2) result (l)
    real(default) :: l
    real(default), intent(in) :: eps2, w2
    l = one + w2**2 + eps2**2 - two * (w2 + eps2 + w2 * eps2)
  end function ff_lambda

  elemental module function ff1 (eps2, w2) result (f)
    real(default) :: f
    real(default), intent(in) :: eps2, w2
    real(default) :: uq, uw, sq_lam, fff
    sq_lam = sqrt (ff_lambda (eps2, w2))
    fff = ff_f0 (eps2, w2)
    uw = (one - eps2 + w2 - sq_lam) / &
         (one - eps2 + w2 + sq_lam)
    uq = (one + eps2 - w2 - sq_lam) / &
         (one + eps2 - w2 + sq_lam)
    f = one / two * fff * (one + eps2 - w2) * &
         (pi**2 + two * Li2 (uw) - two * Li2 (one - uw) - four * Li2 (uq) &
          - four * Li2 (uq * uw) + log ((one - uq) / w2) * log (one - uq) &
          - log (one - uq * uw)**2 + one / four * log (w2 / uw)**2 &
          - log (uw) * log ((one - uq * uw)**2 / (one - uq)) &
          - two * log (uq) * log ((one - uq) * (one - uq * uw))) &
          - sq_lam * fff * (two * log (sqrt (w2)) &
          + three * log (sqrt (eps2)) - two * log (sq_lam**2)) &
          + four * (one - eps2) * ((one - eps2)**2 + w2 * (one + eps2) &
          - four * w2**2) * log (uw) &
          + (three - eps2 + 11 * eps2**2 - eps2**3 + w2 * &
            (6 - 12 * eps2 + 2 * eps2**2) - w2**2 * (21 + five * eps2) &
          + 12 * w2**3) * log (uq) &
          + 6 * sq_lam * (one - eps2) * &
            (one + eps2 - w2) * log (sqrt (eps2)) &
          + sq_lam * (- five + 22 * eps2 - five * eps2**2 - 9 * w2 * &
            (one + eps2) + 6 * w2**2)
  end function ff1


end submodule sm_physics_s

