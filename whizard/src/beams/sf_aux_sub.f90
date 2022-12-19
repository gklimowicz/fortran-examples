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

submodule (sf_aux) sf_aux_s

  use io_units
  use numeric_utils

  implicit none

contains

  module subroutine splitting_data_write (d, unit)
    class(splitting_data_t), intent(in) :: d
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A)") "Splitting data:"
    write (u, "(2x,A,L1)")  "collinear = ", d%collinear
1   format (2x,A,1x,ES15.8)
    write (u, 1) "x0   =", d%x0
    write (u, 1) "x    =", d%x
    write (u, 1) "xb   =", d%xb
    write (u, 1) "x1   =", d%x1
    write (u, 1) "t0   =", d%t0
    write (u, 1) "t    =", d%t
    write (u, 1) "t1   =", d%t1
    write (u, 1) "phi0 =", d%phi0
    write (u, 1) "phi  =", d%phi
    write (u, 1) "phi1 =", d%phi1
    write (u, 1) "E    =", d%E
    write (u, 1) "p    =", d%p
    write (u, 1) "pb   =", d%pb
    write (u, 1) "s    =", d%s
    write (u, 1) "u    =", d%u
    write (u, 1) "m2   =", d%m2
  end subroutine splitting_data_write

  module subroutine splitting_data_init (d, k, mk2, mr2, mo2, collinear)
    class(splitting_data_t), intent(out) :: d
    type(vector4_t), intent(in) :: k
    real(default), intent(in) :: mk2, mr2, mo2
    logical, intent(in), optional :: collinear
    if (present (collinear))  d%collinear = collinear
    d%E = energy (k)
    d%x1 = 1 - sqrt (max (mr2, 0._default)) / d%E
    d%p = sqrt (d%E**2 - mk2)
    d%s = mk2
    d%u = mr2
    d%m2 = mo2
  end subroutine splitting_data_init

  module function splitting_get_x_bounds (d) result (x)
    class(splitting_data_t), intent(in) :: d
    real(default), dimension(2) :: x
    x = [ d%x0, d%x1 ]
  end function splitting_get_x_bounds

  elemental module subroutine splitting_set_t_bounds (d, x, xb)
    class(splitting_data_t), intent(inout) :: d
    real(default), intent(in), optional :: x, xb
    real(default) :: tp, tm
    if (present (x))  d%x = x
    if (present (xb)) d%xb = xb
    if (vanishes (d%u)) then
       d%pb = d%E
    else
       if (.not. vanishes (d%xb)) then
          d%pb = sqrt (max (d%E**2 - d%u / d%xb**2, 0._default))
       else
          d%pb = 0
       end if
    end if
    tp = -2 * d%xb * d%E**2 + d%s + d%u
    tm = -2 * d%xb * d%p * d%pb
    d%t0 = tp + tm
    d%t1 = tp - tm
    d%t = d%t1
  end subroutine splitting_set_t_bounds

  module subroutine splitting_sample_t (d, r, t0, t1)
    class(splitting_data_t), intent(inout) :: d
    real(default), intent(in) :: r
    real(default), intent(in), optional :: t0, t1
    real(default) :: tt0, tt1, tt0m, tt1m
    if (d%collinear) then
       d%t = d%t1
    else
       tt0 = d%t0;  if (present (t0))  tt0 = max (t0, tt0)
       tt1 = d%t1;  if (present (t1))  tt1 = min (t1, tt1)
       tt0m = tt0 - d%m2
       tt1m = tt1 - d%m2
       if (tt0m < 0 .and. tt1m < 0 .and. abs(tt0m) > &
            epsilon(tt0m) .and. abs(tt1m) > epsilon(tt0m)) then
          d%t = d%m2 + tt0m * exp (r * log (tt1m / tt0m))
       else
          d%t = tt1
       end if
    end if
  end subroutine splitting_sample_t

  module subroutine splitting_inverse_t (d, r, t0, t1)
    class(splitting_data_t), intent(in) :: d
    real(default), intent(out) :: r
    real(default), intent(in), optional :: t0, t1
    real(default) :: tt0, tt1, tt0m, tt1m
    if (d%collinear) then
       r = 0
    else
       tt0 = d%t0;  if (present (t0))  tt0 = max (t0, tt0)
       tt1 = d%t1;  if (present (t1))  tt1 = min (t1, tt1)
       tt0m = tt0 - d%m2
       tt1m = tt1 - d%m2
       if (tt0m < 0 .and. tt1m < 0) then
          r = log ((d%t - d%m2) / tt0m) / log (tt1m / tt0m)
       else
          r = 0
       end if
    end if
  end subroutine splitting_inverse_t

  module subroutine splitting_sample_phi (d, r)
    class(splitting_data_t), intent(inout) :: d
    real(default), intent(in) :: r
    if (d%collinear) then
       d%phi = 0
    else
       d%phi = (1-r) * d%phi0 + r * d%phi1
    end if
  end subroutine splitting_sample_phi

  module subroutine splitting_inverse_phi (d, r)
    class(splitting_data_t), intent(in) :: d
    real(default), intent(out) :: r
    if (d%collinear) then
       r = 0
    else
       r = (d%phi - d%phi0) / (d%phi1 - d%phi0)
    end if
  end subroutine splitting_inverse_phi

  module function splitting_split_momentum (d, k) result (q)
    class(splitting_data_t), intent(in) :: d
    type(vector4_t), dimension(2) :: q
    type(vector4_t), intent(in) :: k
    real(default) :: st2, ct2, st, ct, cp, sp
    type(lorentz_transformation_t) :: rot
    real(default) :: tt0, tt1, den
    type(vector3_t) :: kk, q1, q2
    if (d%collinear) then
       if (vanishes (d%s) .and. vanishes(d%u)) then
          q(1) = d%xb * k
          q(2) = d%x * k
       else
          kk = space_part (k)
          q1 = d%xb * (d%pb / d%p) * kk
          q2 = kk - q1
          q(1) = vector4_moving (d%xb * d%E, q1)
          q(2) = vector4_moving (d%x * d%E, q2)
       end if
    else
       den = 2 * d%xb * d%p * d%pb
       tt0 = max (d%t - d%t0, 0._default)
       tt1 = min (d%t - d%t1, 0._default)
       if (den**2 <= epsilon(den)) then
          st2 = 0
       else
          st2 = - (tt0 * tt1) / den ** 2
       end if
       if (st2 > 1) then
          st2 = 1
       end if
       ct2 = 1 - st2
       st = sqrt (max (st2, 0._default))
       ct = sqrt (max (ct2, 0._default))
       if ((d%t - d%t0 + d%t - d%t1) < 0) then
          ct = - ct
       end if
       sp = sin (d%phi)
       cp = cos (d%phi)
       rot = rotation_to_2nd (3, space_part (k))
       q1 = vector3_moving (d%xb * d%pb * [st * cp, st * sp, ct])
       q2 = vector3_moving (d%p, 3) - q1
       q(1) = rot * vector4_moving (d%xb * d%E, q1)
       q(2) = rot * vector4_moving (d%x * d%E, q2)
    end if
  end function splitting_split_momentum

  elemental module subroutine on_shell (p, m2, keep)
    type(vector4_t), intent(inout) :: p
    real(default), intent(in) :: m2
    integer, intent(in) :: keep
    real(default) :: E, E2, pn
    select case (keep)
    case (KEEP_ENERGY)
       E = energy (p)
       E2 = E ** 2
       if (E2 >= m2) then
          pn = sqrt (E2 - m2)
          p = vector4_moving (E, pn * direction (space_part (p)))
       else
          p = vector4_null
       end if
    case (KEEP_MOMENTUM)
       E = sqrt (space_part (p) ** 2 + m2)
       p = vector4_moving (E, space_part (p))
    end select
  end subroutine on_shell

  module subroutine splitting_recover (d, k, q, keep)
    class(splitting_data_t), intent(inout) :: d
    type(vector4_t), intent(in) :: k
    type(vector4_t), dimension(2), intent(in) :: q
    integer, intent(in) :: keep
    type(lorentz_transformation_t) :: rot
    type(vector4_t) :: k0
    type(vector4_t), dimension(2) :: q0
    real(default) :: p1, p2, p3, pt2, pp2, pl
    real(default) :: aux, den, norm
    real(default) :: st2, ct2, ct
    rot = inverse (rotation_to_2nd (3, space_part (k)))
    q0 = rot * q
    p1 = vector4_get_component (q0(2), 1)
    p2 = vector4_get_component (q0(2), 2)
    p3 = vector4_get_component (q0(2), 3)
    pt2 = p1 ** 2 + p2 ** 2
    pp2 = p1 ** 2 + p2 ** 2 + p3 ** 2
    pl = abs (p3)
    k0 = vector4_moving (d%E, d%p, 3)
    select case (keep)
    case (KEEP_ENERGY)
       d%x = energy (q0(2)) / d%E
       d%xb = energy (q0(1)) / d%E
       call d%set_t_bounds ()
       if (.not. d%collinear) then
          aux = (d%xb * d%pb) ** 2 * pp2 - d%p ** 2 * pt2
          den = d%p ** 2 - (d%xb * d%pb) ** 2
          if (aux >= 0 .and. den > 0) then
             norm = (d%p * pl + sqrt (aux)) / den
          else
             norm = 1
          end if
       end if
    case (KEEP_MOMENTUM)
       d%xb = sqrt (space_part (q0(1)) ** 2 + d%u) / d%E
       d%x = 1 - d%xb
       call d%set_t_bounds ()
       norm = 1
    end select
    if (d%collinear) then
       d%t = d%t1
       d%phi = 0
    else
       if ((d%xb * d%pb * norm)**2 < epsilon(d%xb)) then
          st2 = 1
       else
          st2 = pt2 / (d%xb * d%pb * norm ) ** 2
       end if
       if (st2 > 1) then
          st2 = 1
       end if
       ct2 = 1 - st2
       ct = sqrt (max (ct2, 0._default))
       if (.not. vanishes (1 + ct)) then
          d%t = d%t1 - 2 * d%xb * d%p * d%pb * st2 / (1 + ct)
       else
          d%t = d%t0
       end if
       if (.not. vanishes (p1) .or. .not. vanishes (p2)) then
          d%phi = atan2 (-p2, -p1)
       else
          d%phi = 0
       end if
    end if
  end subroutine splitting_recover

  module function splitting_get_x (sd) result (x)
    class(splitting_data_t), intent(in) :: sd
    real(default) :: x
    x = sd%x
  end function splitting_get_x

  module function splitting_get_xb (sd) result (xb)
    class(splitting_data_t), intent(in) :: sd
    real(default) :: xb
    xb = sd%xb
  end function splitting_get_xb


end submodule sf_aux_s

