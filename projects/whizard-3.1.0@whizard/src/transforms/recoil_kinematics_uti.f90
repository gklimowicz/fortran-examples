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

module recoil_kinematics_uti

  use kinds, only: default
  use constants, only: twopi
  use constants, only: degree
  use lorentz, only: vector4_t
  use lorentz, only: vector4_moving
  use lorentz, only: lorentz_transformation_t
  use lorentz, only: inverse
  use lorentz, only: operator(+)
  use lorentz, only: operator(*)
  use lorentz, only: operator(**)
  use lorentz, only: pacify

  use recoil_kinematics, only: solve_recoil
  use recoil_kinematics, only: recoil_momenta
  use recoil_kinematics, only: recoil_transformation
  use recoil_kinematics, only: initial_transformation
  use recoil_kinematics, only: generate_q2_recoil
  use recoil_kinematics, only: generate_recoil

  implicit none
  private

  public :: recoil_kinematics_1
  public :: recoil_kinematics_2
  public :: recoil_kinematics_3
  public :: recoil_kinematics_4
  public :: recoil_kinematics_5
  public :: recoil_kinematics_6

contains

  subroutine recoil_kinematics_1 (u)
    integer, intent(in) :: u

    real(default) :: sqrts
    real(default), dimension(2) :: xc, xcb
    real(default), dimension(2) :: q
    real(default), dimension(2) :: phi
    real(default), dimension(2) :: mo
    real(default), dimension(2) :: cos_th, sin_th
    real(default), dimension(2) :: x
    real(default), dimension(2) :: xb
    type(vector4_t), dimension(2) :: km
    type(vector4_t), dimension(2) :: qm
    type(vector4_t), dimension(2) :: qo
    integer :: i
    logical :: ok

    character(*), parameter :: FMT1 = "(1x,A,9(1x,F15.10))"
    character(*), parameter :: FMT2 = "(1x,A,9(1x,F10.5))"
    character(*), parameter :: FMT4 = "(3x,ES8.1,9(1x,ES19.12))"

    write (u, "(A)")  "* Test output: recoil_kinematics_1"
    write (u, "(A)")  "*   Purpose: compute kinematics for various input data"
    write (u, "(A)")

    sqrts = 100
    write (u, FMT1)  "sqrts =", sqrts

    write (u, "(A)")
    write (u, "(A)") "*** collinear data set"
    write (u, "(A)")

    xc = [0.6_default, 0.9_default]
    xcb = 1 - xc
    phi = [0.1_default, 0.2_default] * twopi
    q = 0
    mo = 0

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    call show_results

    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    call show_momenta

    write (u, "(A)")
    write (u, "(A)") "*** moderate data set"
    write (u, "(A)")

    xc = [0.6_default, 0.9_default]
    xcb = 1 - xc
    phi = [0.1_default, 0.2_default] * twopi
    q = [0.2_default, 0.05_default] * sqrts

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    call show_results
    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    call show_momenta

    write (u, "(A)")
    write (u, "(A)") "*** semi-soft data set"
    write (u, "(A)")

    xcb= [0.1_default, 0.0001_default]
    xc = 1 - xcb
    phi = [0.1_default, 0.2_default] * twopi
    q = [0.2_default, 0.00001_default] * sqrts

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    call show_results

    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    call show_momenta

    write (u, "(A)")
    write (u, "(A)") "*** hard-soft data set"
    write (u, "(A)")

    xcb= [0.1_default, 1.e-30_default]
    xc = 1 - xcb
    phi = [0.1_default, 0.2_default] * twopi
    q = [0.2_default, 1.e-35_default] * sqrts

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    call show_results

    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    call show_momenta

    write (u, "(A)")
    write (u, "(A)") "*** hard data set"
    write (u, "(A)")

    xc = [0.2_default, 0.4_default]
    xcb = 1 - xc
    phi = [0.1_default, 0.8_default] * twopi
    q = [0.74_default, 0.3_default] * sqrts

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    call show_results

    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    call show_momenta

    write (u, "(A)")
    write (u, "(A)") "*** failing data set"
    write (u, "(A)")

    xc = [0.2_default, 0.4_default]
    xcb = 1 - xc
    phi = [0.1_default, 0.8_default] * twopi
    q = [0.9_default, 0.3_default] * sqrts

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    if (.not. ok) then
       write (u, "(A)")
       write (u, "(A)") "Failed as expected."
    end if

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: recoil_kinematics_1"

  contains

    subroutine show_data
      write (u, FMT1) "sqs_h =", sqrt (xc(1) * xc(2)) * sqrts
      write (u, FMT1) "xc    =", xc
      write (u, FMT1) "xcb   =", xcb
      write (u, FMT1) "Q     =", Q
      write (u, FMT1) "phi/D =", phi / degree
    end subroutine show_data

    subroutine show_results
      write (u, "(A)")
      write (u, "(A)") "Result:"
      write (u, FMT1) "th/D  =", atan2 (sin_th, cos_th) / degree
      write (u, FMT1) "x     =", x
      write (u, "(A)")
    end subroutine show_results

    subroutine show_momenta
      type(vector4_t) :: qm0, qo0
      real(default), parameter :: tol = 1.e-7_default
      call pacify (km, tol)
      call pacify (qm, tol)
      call pacify (qo, tol)
      write (u, "(A)") "Momenta: k"
      call km(1)%write (u, testflag=.true.)
      call km(2)%write (u, testflag=.true.)
      write (u, FMT1) "k^2 =", abs (km(1)**2), abs (km(2)**2)
      write (u, "(A)")
      write (u, "(A)") "Momenta: q"
      call qm(1)%write (u, testflag=.true.)
      call qm(2)%write (u, testflag=.true.)
      write (u, "(A)")
      write (u, "(A)") "Momenta: q(os)"
      call qo(1)%write (u, testflag=.true.)
      call qo(2)%write (u, testflag=.true.)
      write (u, "(A)")
      write (u, "(A)") "Check: parton momentum sum: q vs q(os)"
      qm0 = qm(1) + qm(2)
      call qm0%write (u, testflag=.true.)
      qo0 = qo(1) + qo(2)
      call qo0%write (u, testflag=.true.)
      write (u, "(A)")
      write (u, "(A)") "* Check: momentum transfer (off-shell/on-shell)"
      write (u, FMT2) "|q| =", abs (qm(1)**1), abs (qm(2)**1)
      write (u, FMT2) "Q   =", q
      write (u, FMT2) "|qo|=", abs (qo(1)**1), abs (qo(2)**1)
      write (u, "(A)")
      write (u, "(A)") "* Check: sqrts, sqrts_hat"
      write (u, FMT1) "|p| =", (km(1)+km(2)+qm(1)+qm(2))**1, (qm(1)+qm(2))**1
      write (u, FMT1) "sqs =", sqrts, sqrt (product (xc)) * sqrts
      write (u, FMT1) "|po|=", abs ((km(1)+km(2)+qo(1)+qo(2))**1), abs ((qo(1)+qo(2))**1)
    end subroutine show_momenta

  end subroutine recoil_kinematics_1

  subroutine recoil_kinematics_2 (u)
    integer, intent(in) :: u

    real(default) :: sqrts
    real(default) :: q_max
    real(default) :: m
    real(default) :: x_bar
    real(default) :: r
    real(default) :: q2, q2_old
    integer :: i
    integer :: n_bin

    character(*), parameter :: FMT1 = "(1x,A,9(1x,F15.10))"
    character(*), parameter :: FMT3 = "(2x,9(1x,F10.5))"

    write (u, "(A)")  "* Test output: recoil_kinematics_2"
    write (u, "(A)")  "*   Purpose: compute Q distribution"
    write (u, "(A)")

    n_bin = 20

    write (u, "(A)") "* No Q cutoff, xbar = 1"
    write (u, "(A)")

    sqrts = 100
    q_max = sqrts
    m = 0.511e-3_default
    x_bar = 1._default
    call show_table

    write (u, "(A)")
    write (u, "(A)") "* With Q cutoff, xbar = 1"
    write (u, "(A)")

    q_max = 10
    call show_table

    write (u, "(A)")
    write (u, "(A)") "* No Q cutoff, xbar = 0.01"
    write (u, "(A)")

    q_max = sqrts
    x_bar = 0.01_default
    call show_table

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: recoil_kinematics_2"

  contains

    subroutine show_table
      write (u, FMT1)  "sqrts =", sqrts
      write (u, FMT1)  "q_max =", q_max
      write (u, FMT1)  "m     =", m
      write (u, FMT1)  "x_bar =", x_bar
      write (u, "(A)")
      write (u, "(1x,A)")  "Table:  r         |Q|    |Q_i/Q_(i-1)|"
      q2_old = 0
      do i = 0, n_bin
         r = real (i, default) / n_bin
         call generate_q2_recoil (sqrts**2, x_bar, q_max**2, m**2, r, q2)
         if (q2_old > 0) then
            write (u, FMT3)  r, sqrt (q2), sqrt (q2 / q2_old)
         else
            write (u, FMT3)  r, sqrt (q2)
         end if
         q2_old = q2
      end do
    end subroutine show_table

  end subroutine recoil_kinematics_2

  subroutine recoil_kinematics_3 (u)
    integer, intent(in) :: u

    real(default) :: sqrts
    real(default), dimension(2) :: q_max
    real(default), dimension(2) :: m, mo
    real(default), dimension(2) :: xc, xcb
    real(default), dimension(4) :: r
    type(vector4_t), dimension(2) :: km
    type(vector4_t), dimension(2) :: qm
    type(vector4_t), dimension(2) :: qo
    logical :: ok

    character(*), parameter :: FMT1 = "(1x,A,9(1x,F15.10))"
    character(*), parameter :: FMT2 = "(1x,A,9(1x,F10.5))"

    write (u, "(A)")  "* Test output: recoil_kinematics_3"
    write (u, "(A)")  "*   Purpose: generate momenta from RNG parameters"
    write (u, "(A)")

    write (u, "(A)") "*** collinear data set"
    write (u, "(A)")

    sqrts = 100
    q_max = sqrts
    m     = 0.511e-3_default
    mo    = 0

    xc = [0.6_default, 0.9_default]
    xcb = 1 - xc

    r = [0._default, 0._default, 0._default, 0._default]

    call show_data
    call generate_recoil (sqrts, q_max, m, mo, xc, xcb, r, km, qm, qo, ok)
    call show_momenta

    write (u, "(A)")
    write (u, "(A)") "*** moderate data set"
    write (u, "(A)")

    xc = [0.6_default, 0.9_default]
    xcb = 1 - xc

    r = [0.8_default, 0.2_default, 0.1_default, 0.2_default]

    call show_data
    call generate_recoil (sqrts, q_max, m, mo, xc, xcb, r, km, qm, qo, ok)
    call show_momenta

    write (u, "(A)")
    write (u, "(A)") "*** failing data set"
    write (u, "(A)")

    xc = [0.2_default, 0.4_default]
    xcb = 1 - xc

    r = [0.9999_default, 0.3_default, 0.1_default, 0.8_default]

    call show_data
    call generate_recoil (sqrts, q_max, m, mo, xc, xcb, r, km, qm, qo, ok)
    if (.not. ok) then
       write (u, "(A)")
       write (u, "(A)") "Failed as expected."
    else
       call show_momenta
    end if

  contains

    subroutine show_data
      write (u, FMT1)  "sqrts =", sqrts
      write (u, FMT1)  "q_max =", q_max
      write (u, FMT1)  "m     =", m
      write (u, FMT1)  "xc    =", xc
      write (u, FMT1)  "xcb   =", xcb
      write (u, FMT1)  "r     =", r
    end subroutine show_data

    subroutine show_momenta
      real(default), parameter :: tol = 1.e-7_default
      call pacify (km, tol)
      call pacify (qo, tol)
      write (u, "(A)")
      write (u, "(A)") "* Momenta: k"
      call km(1)%write (u, testflag=.true.)
      call km(2)%write (u, testflag=.true.)
      write (u, FMT1) "k^2 =", abs (km(1)**2), abs (km(2)**2)
      write (u, "(A)")
      write (u, "(A)") "* Momenta: q(os)"
      call qo(1)%write (u, testflag=.true.)
      call qo(2)%write (u, testflag=.true.)
      write (u, FMT1) "q^2 =", abs (qo(1)**2), abs (qo(2)**2)
      write (u, "(A)")
      write (u, "(A)") "* Check: momentum transfer (off-shell/on-shell)"
      write (u, FMT2) "Q   =", q_check (1), q_check (2)
      write (u, FMT2) "|q| =", abs (qm(1)**1), abs (qm(2)**1)
      write (u, "(A)")
      write (u, "(A)") "* Check: sqrts, sqrts_hat"
      write (u, FMT1) "sqs =", sqrts, sqrt (product (xc)) * sqrts
      write (u, FMT1) "|po|=", abs ((km(1)+km(2)+qo(1)+qo(2))**1), abs ((qo(1)+qo(2))**1)
    end subroutine show_momenta

    function q_check (i) result (q)
      integer, intent(in) :: i
      real(default) :: q
      real(default) :: q2
      call generate_q2_recoil (sqrts**2, xcb(i), q_max(i)**2, m(i)**2, r(i), q2)
      q = sqrt (q2)
    end function q_check

  end subroutine recoil_kinematics_3

  subroutine recoil_kinematics_4 (u)
    integer, intent(in) :: u

    real(default) :: sqrts
    real(default), dimension(2) :: xc, xcb
    real(default), dimension(2) :: q
    real(default), dimension(2) :: phi
    real(default), dimension(2) :: cos_th, sin_th
    real(default), dimension(2) :: mo
    real(default), dimension(2) :: x
    real(default), dimension(2) :: xb
    type(vector4_t), dimension(2) :: km
    type(vector4_t), dimension(2) :: qm
    type(vector4_t), dimension(2) :: qo
    type(lorentz_transformation_t) :: lt
    logical :: ok

    character(*), parameter :: FMT1 = "(1x,A,9(1x,F15.10))"
    character(*), parameter :: FMT2 = "(1x,A,9(1x,F10.5))"

    write (u, "(A)")  "* Test output: recoil_kinematics_4"
    write (u, "(A)")  "*   Purpose: check Lorentz transformation for recoil"
    write (u, "(A)")

    sqrts = 100
    write (u, FMT1)  "sqrts =", sqrts

    write (u, "(A)")
    write (u, "(A)") "*** collinear data set"
    write (u, "(A)")

    xc = [0.6_default, 0.9_default]
    xcb = 1 - xc
    phi = [0.1_default, 0.2_default] * twopi
    q = 0
    mo = 0

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    call recoil_transformation (sqrts, xc, qo, lt)
    call show_transformation

    write (u, "(A)")
    write (u, "(A)") "*** moderate data set"
    write (u, "(A)")

    xc = [0.6_default, 0.9_default]
    xcb = 1 - xc
    phi = [0.1_default, 0.2_default] * twopi
    q = [0.2_default, 0.05_default] * sqrts

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    call recoil_transformation (sqrts, xc, qo, lt)
    call show_transformation

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: recoil_kinematics_4"

  contains

    subroutine show_data
      write (u, FMT1) "sqs_h =", sqrt (xc(1) * xc(2)) * sqrts
      write (u, FMT1) "xc    =", xc
      write (u, FMT1) "xcb   =", xcb
      write (u, FMT1) "Q     =", Q
      write (u, FMT1) "phi/D =", phi / degree
    end subroutine show_data

    subroutine show_transformation
      type(vector4_t), dimension(2) :: qc
      type(vector4_t), dimension(2) :: qct
      real(default), parameter :: tol = 1.e-7_default
      qc(1) = xc(1) * vector4_moving (sqrts/2, sqrts/2, 3)
      qc(2) = xc(2) * vector4_moving (sqrts/2,-sqrts/2, 3)
      qct = lt * qc
      call pacify (qct, tol)
      write (u, "(A)")
      write (u, "(A)") "Momenta: q(os)"
      call qo(1)%write (u, testflag=.true.)
      call qo(2)%write (u, testflag=.true.)
      write (u, "(A)")
      write (u, "(A)") "Momenta: LT * qc"
      call qct(1)%write (u, testflag=.true.)
      call qct(2)%write (u, testflag=.true.)
    end subroutine show_transformation

  end subroutine recoil_kinematics_4

  subroutine recoil_kinematics_5 (u)
    integer, intent(in) :: u

    real(default) :: sqrts
    real(default) :: sqrtsi
    real(default), dimension(2) :: x
    type(vector4_t), dimension(2) :: p
    type(vector4_t), dimension(2) :: pi
    type(vector4_t), dimension(2) :: p0
    type(lorentz_transformation_t) :: lt
    logical :: ok

    character(*), parameter :: FMT1 = "(1x,A,9(1x,F15.10))"
    character(*), parameter :: FMT2 = "(1x,A,9(1x,F10.5))"

    write (u, "(A)")  "* Test output: recoil_kinematics_5"
    write (u, "(A)")  "*   Purpose: determine initial Lorentz transformation"
    write (u, "(A)")

    sqrts = 100
    write (u, FMT1)  "sqrts =", sqrts

    x = [0.6_default, 0.9_default]

    p(1) = x(1) * vector4_moving (sqrts/2, sqrts/2, 3)
    p(2) = x(2) * vector4_moving (sqrts/2,-sqrts/2, 3)

    call show_data
    call initial_transformation (p, sqrtsi, lt, ok)

    pi(1) = vector4_moving (sqrtsi/2, sqrtsi/2, 3)
    pi(2) = vector4_moving (sqrtsi/2,-sqrtsi/2, 3)

    p0 = inverse (lt) * p

    call show_momenta

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: recoil_kinematics_5"

  contains

    subroutine show_data
      write (u, FMT1) "sqrts =", sqrts
      write (u, FMT1) "x     =", x
    end subroutine show_data


    subroutine show_momenta
      real(default), parameter :: tol = 1.e-7_default
      write (u, "(A)")
      write (u, "(A)") "* Momenta: p_in(c.m.)"
      call pi(1)%write (u, testflag=.true.)
      call pi(2)%write (u, testflag=.true.)
      write (u, "(A)")
      write (u, "(A)") "* Momenta: inv(LT) * p_in(lab)"
      call p0(1)%write (u, testflag=.true.)
      call p0(2)%write (u, testflag=.true.)
    end subroutine show_momenta

  end subroutine recoil_kinematics_5

  subroutine recoil_kinematics_6 (u)
    integer, intent(in) :: u

    real(default) :: sqrts
    real(default), dimension(2) :: xc, xcb
    real(default), dimension(2) :: q
    real(default), dimension(2) :: phi
    real(default), dimension(2) :: cos_th, sin_th
    real(default), dimension(2) :: x
    real(default), dimension(2) :: xb
    real(default), dimension(2) :: mo, z
    type(vector4_t), dimension(2) :: km
    type(vector4_t), dimension(2) :: qm
    type(vector4_t), dimension(2) :: qo
    type(lorentz_transformation_t) :: lt
    logical :: ok

    character(*), parameter :: FMT1 = "(1x,A,9(1x,F15.10))"
    character(*), parameter :: FMT2 = "(1x,A,9(1x,F11.6))"

    write (u, "(A)")  "* Test output: recoil_kinematics_6"
    write (u, "(A)")  "*   Purpose: check effect of mass in on-shell projection"
    write (u, "(A)")

    sqrts = 10
    write (u, FMT1)  "sqrts =", sqrts
    z = 0
    mo = 0.511e-3
    write (u, FMT1)  "mass  =", mo

    write (u, "(A)")
    write (u, "(A)") "*** collinear data set"
    write (u, "(A)")

    xc = [0.6_default, 0.9_default]
    xcb = 1 - xc
    phi = [0.1_default, 0.2_default] * twopi
    q = 0

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, z, km, qm, qo, ok)
    call recoil_transformation (sqrts, xc, qo, lt)
    write (u, "(A)")
    write (u, "(A)") "Massless projection:"
    call show_momenta

    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    call recoil_transformation (sqrts, xc, qo, lt)
    write (u, "(A)")
    write (u, "(A)") "Massive projection:"
    call show_momenta

    write (u, "(A)")
    write (u, "(A)") "*** moderate data set"
    write (u, "(A)")

    xc = [0.6_default, 0.9_default]
    xcb = 1 - xc
    phi = [0.1_default, 0.2_default] * twopi
    q = [0.2_default, 0.05_default] * sqrts

    call show_data
    call solve_recoil (sqrts, xc, xcb, phi, q**2, x, xb, cos_th, sin_th, ok)
    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, z, km, qm, qo, ok)
    call recoil_transformation (sqrts, xc, qo, lt)
    write (u, "(A)")
    write (u, "(A)") "Massless projection:"
    call show_momenta

    call recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    call recoil_transformation (sqrts, xc, qo, lt)
    write (u, "(A)")
    write (u, "(A)") "Massive projection:"
    call show_momenta

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: recoil_kinematics_6"

  contains

    subroutine show_data
      write (u, FMT1) "sqs_h =", sqrt (xc(1) * xc(2)) * sqrts
      write (u, FMT1) "xc    =", xc
      write (u, FMT1) "xcb   =", xcb
      write (u, FMT1) "Q     =", Q
      write (u, FMT1) "phi/D =", phi / degree
    end subroutine show_data

    subroutine show_momenta
      write (u, "(A)") "Momenta: q(os)"
      call qo(1)%write (u, testflag=.true.)
      write (u, FMT2)  "m = ", abs (qo(1)**1)
      call qo(2)%write (u, testflag=.true.)
      write (u, FMT2)  "m = ", abs (qo(2)**1)
    end subroutine show_momenta

  end subroutine recoil_kinematics_6


end module recoil_kinematics_uti

