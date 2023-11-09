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

submodule (recoil_kinematics) recoil_kinematics_s

  use constants, only: twopi

  implicit none

contains

  elemental subroutine generate_phi_recoil (r, phi)
    real(default), intent(in) :: r
    real(default), intent(out) :: phi

    phi = r * twopi

  end subroutine generate_phi_recoil

  elemental module subroutine generate_q2_recoil (s, x_bar, q2_max, m2, r, q2)
    real(default), intent(in) :: s
    real(default), intent(in) :: q2_max
    real(default), intent(in) :: x_bar
    real(default), intent(in) :: m2
    real(default), intent(in) :: r
    real(default), intent(out) :: q2

    real(default) :: q2_max_evt

    q2_max_evt = q2_max_event (s, x_bar, q2_max)

    q2 = m2 * (exp (r * log (1 + (q2_max_evt / m2))) - 1)

  end subroutine generate_q2_recoil

  elemental function q2_max_event (s, x_bar, q2_max) result (q2)
    real(default), intent(in) :: s
    real(default), intent(in) :: x_bar
    real(default), intent(in) :: q2_max
    real(default) :: q2

    q2 = min (x_bar * s, q2_max)

  end function q2_max_event

  subroutine polar_angles (s, xb, rho, ee, q2, sin_th, cos_th, ok)
    real(default), intent(in) :: s
    real(default), intent(in) :: xb
    real(default), intent(in) :: rho
    real(default), dimension(2), intent(in) :: ee
    real(default), dimension(2), intent(in) :: q2
    real(default), dimension(2), intent(out) :: sin_th
    real(default), dimension(2), intent(out) :: cos_th
    logical, intent(out) :: ok

    real(default), dimension(2) :: sin2_th_2

    sin2_th_2 = q2 / (ee * rho * xb * s)

    if (all (sin2_th_2 <= 1)) then
       sin_th = 2 * sqrt (sin2_th_2 * (1 - sin2_th_2))
       cos_th = 1 - 2 * sin2_th_2
       ok = .true.
    else
       sin_th = 0
       cos_th = 1
       ok = .false.
    end if

  end subroutine polar_angles

  function lambda_factor (sin_th, cos_th, cphi) result (lambda)
    real(default), dimension(2), intent(in) :: sin_th
    real(default), dimension(2), intent(in) :: cos_th
    real(default), intent(in) :: cphi
    real(default) :: lambda

    lambda = (1 - cos_th(1) * cos_th(2) - cphi * sin_th(1) * sin_th(2)) / 2

  end function lambda_factor

  function scale_factor (che, lambda, xb0, approximate) result (rho)
    real(default), intent(in) :: che
    real(default), intent(in) :: lambda
    real(default), intent(in) :: xb0
    logical, intent(in), optional :: approximate
    real(default) :: rho

    real(default), parameter :: &
         e0 = (100 * epsilon (1._default)) ** (0.3_default)
    logical :: approx

    if (present (approximate)) then
       approx = approximate
    else
       approx = (xb0/che) < e0
    end if

    if (approx) then
       rho = 1 - lambda * (xb0/(2*che)) * (1 + (1-lambda) * (xb0/che))
    else
       rho =  (che / ((1-lambda)*xb0)) &
         * (1 - sqrt (1 - 2 * (1-lambda) * (xb0/che) &
         &    + (1-lambda) * (xb0 / che)**2))
    end if

  end function scale_factor

  subroutine scaled_x (rho, ee, xb0, x, xb)
    real(default), intent(in) :: rho
    real(default), dimension(2), intent(in) :: ee
    real(default), intent(in) :: xb0
    real(default), dimension(2), intent(out) :: x
    real(default), dimension(2), intent(out) :: xb

    xb = rho * ee * xb0
    x = 1 - xb

  end subroutine scaled_x

  module subroutine solve_recoil &
       (sqrts, xc, xcb, phi, q2, x, xb, cos_th, sin_th, ok)
    real(default), intent(in) :: sqrts
    real(default), dimension(2), intent(in) :: xc
    real(default), dimension(2), intent(in) :: xcb
    real(default), dimension(2), intent(in) :: phi
    real(default), dimension(2), intent(in) :: q2
    real(default), dimension(2), intent(out) :: x
    real(default), dimension(2), intent(out) :: xb
    real(default), dimension(2), intent(out) :: cos_th
    real(default), dimension(2), intent(out) :: sin_th
    logical, intent(out) :: ok

    real(default) :: s
    real(default), dimension(2) :: ee
    real(default), dimension(2) :: th
    real(default) :: xb0, cphi
    real(default) :: che, lambda
    real(default) :: rho_new, rho, rho_old
    real(default) :: dr_old, dr_new
    real(default), parameter :: dr_limit = 100 * epsilon (1._default)
    integer, parameter :: n_it_max = 20
    integer :: i

    ok = .true.

    s = sqrts**2
    ee = sqrt ([xcb(1)/xcb(2), xcb(2)/xcb(1)])
    che = sum (ee) / 2
    xb0 = sqrt (xcb(1) * xcb(2))
    cphi = cos (phi(1) - phi(2))

    rho_old = 10
    rho = 1
    th = 0
    sin_th = sin (th)
    cos_th = cos (th)
    lambda = lambda_factor (sin_th, cos_th, cphi)
    call scaled_x (rho, ee, xb0, x, xb)

    iterate_loop: do i = 1, n_it_max

       call polar_angles (s, xb0, rho, ee, q2, sin_th, cos_th, ok)
       if (.not. ok)  return
       th = atan2 (sin_th, cos_th)

       lambda = lambda_factor (sin_th, cos_th, cphi)
       rho_new = scale_factor (che, lambda, xb0)
       call scaled_x (rho_new, ee, xb0, x, xb)

       dr_old = abs (rho - rho_old)
       dr_new = abs (rho_new - rho)

       rho_old = rho
       rho = rho_new

       if (dr_new < dr_limit .or. dr_new >= dr_old)  exit iterate_loop

    end do iterate_loop

  end subroutine solve_recoil

  module subroutine recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, &
       km, qm, qo, ok)
    real(default), intent(in) :: sqrts
    real(default), dimension(2), intent(in) :: xc
    real(default), dimension(2), intent(in) :: xb
    real(default), dimension(2), intent(in) :: cos_th
    real(default), dimension(2), intent(in) :: sin_th
    real(default), dimension(2), intent(in) :: phi
    real(default), dimension(2), intent(in) :: mo
    type(vector4_t), dimension(2), intent(out) :: km
    type(vector4_t), dimension(2), intent(out) :: qm
    type(vector4_t), dimension(2), intent(out) :: qo
    logical, intent(out) :: ok

    type(vector4_t), dimension(2) :: pm
    type(lorentz_transformation_t) :: lt
    real(default) :: sqsh
    real(default) :: po4, po2
    real(default), dimension(2) :: p0, p3

    pm(1) = &
         vector4_moving (sqrts/2, &
         vector3_moving ([0._default, 0._default, sqrts/2]))
    pm(2) = &
         vector4_moving (sqrts/2, &
         vector3_moving ([0._default, 0._default,-sqrts/2]))

    km(1) = xb(1) * (sqrts/2) * vector4_moving ( &
         1._default, &
         vector3_moving ([ &
         &     sin_th(1) * cos (phi(1)), &
         &     sin_th(1) * sin (phi(1)), &
         &     cos_th(1)]) &
         )
    km(2) = xb(2) * (sqrts/2) * vector4_moving ( &
         1._default, &
         vector3_moving ([ &
         &    -sin_th(2) * cos (phi(2)), &
         &    -sin_th(2) * sin (phi(2)), &
         &    -cos_th(2)]) &
         )

    qm(1) = pm(1) - km(1)
    qm(2) = pm(2) - km(2)

    sqsh = sqrt (xc(1)*xc(2)) * sqrts
    lt = transformation (3, qm(1), qm(2), sqsh)

    po4 = lambda (sqsh**2, mo(1)**2, mo(2)**2)
    ok = po4 > 0
    if (ok) then
       po2 = sqrt (po4)/4
       p0 = sqrt (po2 + mo**2)
       p3 = [sqrt (po2), -sqrt (po2)]
       qo = lt * vector4_moving (p0, p3, 3)
    else
       qo = vector4_null
    end if

  end subroutine recoil_momenta

  module subroutine recoil_transformation (sqrts, xc, qo, lt)
    real(default), intent(in) :: sqrts
    real(default), dimension(2), intent(in) :: xc
    type(vector4_t), dimension(2), intent(in) :: qo
    type(lorentz_transformation_t), intent(out) :: lt

    real(default) :: sqsh
    type(vector4_t), dimension(2) :: qc
    type(lorentz_transformation_t) :: ltc, lto

    qc(1) = xc(1) * vector4_moving (sqrts/2, sqrts/2, 3)
    qc(2) = xc(2) * vector4_moving (sqrts/2,-sqrts/2, 3)

    sqsh = sqrt (xc(1) * xc(2)) * sqrts
    ltc = transformation (3, qc(1), qc(2), sqsh)
    lto = transformation (3, qo(1), qo(2), sqsh)
    lt = lto * inverse (ltc)

  end subroutine recoil_transformation

  module subroutine initial_transformation (p, sqrts, lt, ok)
    type(vector4_t), dimension(2), intent(in) :: p
    real(default), intent(out) :: sqrts
    type(lorentz_transformation_t), intent(out) :: lt
    logical, intent(out) :: ok

    ok = all (transverse_part (p) == 0)

    sqrts = (p(1) + p(2)) ** 1
    lt = boost (p(1) + p(2), sqrts)

  end subroutine initial_transformation

  module subroutine generate_recoil &
       (sqrts, q_max, m, mo, xc, xcb, r, km, qm, qo, ok)
    real(default), intent(in) :: sqrts
    real(default), intent(in), dimension(2) :: q_max
    real(default), intent(in), dimension(2) :: m
    real(default), intent(in), dimension(2) :: mo
    real(default), intent(in), dimension(2) :: xc
    real(default), intent(in), dimension(2) :: xcb
    real(default), intent(in), dimension(4) :: r
    type(vector4_t), dimension(2), intent(out) :: km
    type(vector4_t), dimension(2), intent(out) :: qm
    type(vector4_t), dimension(2), intent(out) :: qo
    logical, intent(out) :: ok

    real(default), dimension(2) :: q2
    real(default), dimension(2) :: phi
    real(default), dimension(2) :: x
    real(default), dimension(2) :: xb
    real(default), dimension(2) :: cos_th
    real(default), dimension(2) :: sin_th

    call generate_q2_recoil (sqrts**2, xcb, q_max**2, m**2, r(1:2), q2)
    call generate_phi_recoil (r(3:4), phi)

    call solve_recoil (sqrts, xc, xcb, phi, q2, x, xb, cos_th, sin_th, ok)
    if (ok) then
       call recoil_momenta &
            (sqrts, xc, xb, cos_th, sin_th, phi, mo, km, qm, qo, ok)
    end if

  end subroutine generate_recoil


end submodule recoil_kinematics_s

