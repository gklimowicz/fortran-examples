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

submodule (lorentz) lorentz_s

  use constants, only: pi, twopi, degree, two, tiny_07, eps0
  use numeric_utils
  use io_units
  use format_defs, only: FMT_11, FMT_13, FMT_15, FMT_19
  use format_utils, only: pac_fmt
  use diagnostics

  implicit none

contains

  module subroutine vector3_write (p, unit, testflag)
    type(vector3_t), intent(in) :: p
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    character(len=7) :: fmt
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    call pac_fmt (fmt, FMT_19, FMT_15, testflag)
    write(u, "(1x,A,3(1x," // fmt // "))") 'P = ', p%p
  end subroutine vector3_write

  elemental module function vector3_canonical (k) result (p)
    type(vector3_t) :: p
    integer, intent(in) :: k
    p = vector3_null
    p%p(k) = 1
  end function vector3_canonical

  elemental module function vector3_moving_canonical (p, k) result(q)
    type(vector3_t) :: q
    real(default), intent(in) :: p
    integer, intent(in) :: k
    q = vector3_null
    q%p(k) = p
  end function vector3_moving_canonical
  pure module function vector3_moving_generic (p) result(q)
    real(default), dimension(3), intent(in) :: p
    type(vector3_t) :: q
    q%p = p
  end function vector3_moving_generic

  elemental module function vector3_eq (p, q) result (r)
    logical :: r
    type(vector3_t), intent(in) :: p,q
    r = all (abs (p%p - q%p) < eps0)
  end function vector3_eq
  elemental module function vector3_neq (p, q) result (r)
    logical :: r
    type(vector3_t), intent(in) :: p,q
    r = any (abs(p%p - q%p) > eps0)
  end function vector3_neq

  elemental module function add_vector3 (p, q) result (r)
    type(vector3_t) :: r
    type(vector3_t), intent(in) :: p,q
    r%p = p%p + q%p
  end function add_vector3
  elemental module function sub_vector3 (p, q) result (r)
    type(vector3_t) :: r
    type(vector3_t), intent(in) :: p,q
    r%p = p%p - q%p
  end function sub_vector3

  elemental module function prod_real_vector3 (s, p) result (q)
    type(vector3_t) :: q
    real(default), intent(in) :: s
    type(vector3_t), intent(in) :: p
    q%p = s * p%p
  end function prod_real_vector3
  elemental module function prod_vector3_real (p, s) result (q)
    type(vector3_t) :: q
    real(default), intent(in) :: s
    type(vector3_t), intent(in) :: p
    q%p = s * p%p
  end function prod_vector3_real
  elemental module function div_vector3_real (p, s) result (q)
    type(vector3_t) :: q
    real(default), intent(in) :: s
    type(vector3_t), intent(in) :: p
    q%p = p%p/s
  end function div_vector3_real
  elemental module function prod_integer_vector3 (s, p) result (q)
    type(vector3_t) :: q
    integer, intent(in) :: s
    type(vector3_t), intent(in) :: p
    q%p = s * p%p
  end function prod_integer_vector3
  elemental module function prod_vector3_integer (p, s) result (q)
    type(vector3_t) :: q
    integer, intent(in) :: s
    type(vector3_t), intent(in) :: p
    q%p = s * p%p
  end function prod_vector3_integer
  elemental module function div_vector3_integer (p, s) result (q)
    type(vector3_t) :: q
    integer, intent(in) :: s
    type(vector3_t), intent(in) :: p
    q%p = p%p/s
  end function div_vector3_integer

  elemental module function prod_vector3 (p, q) result (s)
    real(default) :: s
    type(vector3_t), intent(in) :: p,q
    s = dot_product (p%p, q%p)
  end function prod_vector3

  elemental module function vector3_cross_product (p, q) result (r)
    type(vector3_t) :: r
    type(vector3_t), intent(in) :: p,q
    integer :: i
    do i=1,3
       r%p(i) = dot_product (p%p, matmul(epsilon_three(i,:,:), q%p))
    end do
  end function vector3_cross_product

  elemental module function power_vector3 (p, e) result (s)
    real(default) :: s
    type(vector3_t), intent(in) :: p
    integer, intent(in) :: e
    s = dot_product (p%p, p%p)
    if (e/=2) then
       if (mod(e,2)==0) then
          s = s**(e/2)
       else
          s = sqrt(s)**e
       end if
    end if
  end function power_vector3

  elemental module function negate_vector3 (p) result (q)
    type(vector3_t) :: q
    type(vector3_t), intent(in) :: p
    integer :: i
    do i = 1, 3
       if (abs (p%p(i)) < eps0) then
          q%p(i) = 0
       else
          q%p(i) = -p%p(i)
       end if
    end do
  end function negate_vector3

  module subroutine vector3_set_component (p, i, value)
    type(vector3_t), intent(inout) :: p
    integer, intent(in) :: i
    real(default), intent(in) :: value
    p%p(i) = value
  end subroutine vector3_set_component

  pure module function sum_vector3 (p) result (q)
    type(vector3_t) :: q
    type(vector3_t), dimension(:), intent(in) :: p
    integer :: i
    do i=1, 3
       q%p(i) = sum (p%p(i))
    end do
  end function sum_vector3

  elemental module function vector3_get_component (p, k) result (c)
    type(vector3_t), intent(in) :: p
    integer, intent(in) :: k
    real(default) :: c
    c = p%p(k)
  end function vector3_get_component

  pure module function vector3_get_components (p) result (a)
    type(vector3_t), intent(in) :: p
    real(default), dimension(3) :: a
    a = p%p
  end function vector3_get_components

  elemental module function vector3_get_direction (p) result (q)
    type(vector3_t) :: q
    type(vector3_t), intent(in) :: p
    real(default) :: pp
    pp = p**1
    if (pp > eps0) then
       q%p = p%p / pp
    else
       q%p = 0
    end if
  end function vector3_get_direction

  module subroutine vector4_write &
         (p, unit, show_mass, testflag, compressed, ultra)
    class(vector4_t), intent(in) :: p
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_mass, testflag, compressed, ultra
    logical :: comp, sm, tf, extreme
    integer :: u
    character(len=7) :: fmt
    real(default) :: m
    comp = .false.; if (present (compressed))  comp = compressed
    sm = .false.;  if (present (show_mass))  sm = show_mass
    tf = .false.;  if (present (testflag))  tf = testflag
    extreme = .false.; if (present (ultra))  extreme = ultra
    if (extreme) then
       call pac_fmt (fmt, FMT_19, FMT_11, testflag)
    else
       call pac_fmt (fmt, FMT_19, FMT_13, testflag)
    end if
    u = given_output_unit (unit);  if (u < 0)  return
    if (comp) then
       write (u, "(4(F12.3,1X))", advance="no")  p%p(0:3)
    else
       write (u, "(1x,A,1x," // fmt // ")") 'E = ', p%p(0)
       write (u, "(1x,A,3(1x," // fmt // "))") 'P = ', p%p(1:)
       if (sm) then
          m = p**1
          if (tf)  call pacify (m, tolerance = 1E-6_default)
          write (u, "(1x,A,1x," // fmt // ")") 'M = ', m
       end if
    end if
  end subroutine vector4_write

  module subroutine vector4_write_raw (p, u)
    type(vector4_t), intent(in) :: p
    integer, intent(in) :: u
    write (u) p%p
  end subroutine vector4_write_raw

  module subroutine vector4_read_raw (p, u, iostat)
    type(vector4_t), intent(out) :: p
    integer, intent(in) :: u
    integer, intent(out), optional :: iostat
    read (u, iostat=iostat) p%p
  end subroutine vector4_read_raw

  elemental module function vector4_canonical (k) result (p)
    type(vector4_t) :: p
    integer, intent(in) :: k
    p = vector4_null
    p%p(k) = 1
  end function vector4_canonical

  elemental module function vector4_at_rest (m) result (p)
    type(vector4_t) :: p
    real(default), intent(in) :: m
    p = vector4_t ([ m, zero, zero, zero ])
  end function vector4_at_rest

  elemental module function vector4_moving_canonical (E, p, k) result (q)
    type(vector4_t) :: q
    real(default), intent(in) :: E, p
    integer, intent(in) :: k
    q = vector4_at_rest(E)
    q%p(k) = p
  end function vector4_moving_canonical
  elemental module function vector4_moving_generic (E, p) result (q)
    type(vector4_t) :: q
    real(default), intent(in) :: E
    type(vector3_t), intent(in) :: p
    q%p(0) = E
    q%p(1:) = p%p
  end function vector4_moving_generic

  elemental module function vector4_eq (p, q) result (r)
    logical :: r
    type(vector4_t), intent(in) :: p,q
    r = all (abs (p%p - q%p) < eps0)
  end function vector4_eq
  elemental module function vector4_neq (p, q) result (r)
    logical :: r
    type(vector4_t), intent(in) :: p,q
    r = any (abs (p%p - q%p) > eps0)
  end function vector4_neq

  elemental module function add_vector4 (p,q) result (r)
    type(vector4_t) :: r
    type(vector4_t), intent(in) :: p,q
    r%p = p%p + q%p
  end function add_vector4
  elemental module function sub_vector4 (p,q) result (r)
    type(vector4_t) :: r
    type(vector4_t), intent(in) :: p,q
    r%p = p%p - q%p
  end function sub_vector4

  elemental module function prod_real_vector4 (s, p) result (q)
    type(vector4_t) :: q
    real(default), intent(in) :: s
    type(vector4_t), intent(in) :: p
    q%p = s * p%p
  end function prod_real_vector4
  elemental module function prod_vector4_real (p, s) result (q)
    type(vector4_t) :: q
    real(default), intent(in) :: s
    type(vector4_t), intent(in) :: p
    q%p = s * p%p
  end function prod_vector4_real
  elemental module function div_vector4_real (p, s) result (q)
    type(vector4_t) :: q
    real(default), intent(in) :: s
    type(vector4_t), intent(in) :: p
    q%p = p%p/s
  end function div_vector4_real
  elemental module function prod_integer_vector4 (s, p) result (q)
    type(vector4_t) :: q
    integer, intent(in) :: s
    type(vector4_t), intent(in) :: p
    q%p = s * p%p
  end function prod_integer_vector4
  elemental module function prod_vector4_integer (p, s) result (q)
    type(vector4_t) :: q
    integer, intent(in) :: s
    type(vector4_t), intent(in) :: p
    q%p = s * p%p
  end function prod_vector4_integer
  elemental module function div_vector4_integer (p, s) result (q)
    type(vector4_t) :: q
    integer, intent(in) :: s
    type(vector4_t), intent(in) :: p
    q%p = p%p/s
  end function div_vector4_integer

  elemental module function prod_vector4 (p, q) result (s)
    real(default) :: s
    type(vector4_t), intent(in) :: p,q
    s = p%p(0)*q%p(0) - dot_product(p%p(1:), q%p(1:))
  end function prod_vector4

  elemental module function power_vector4 (p, e) result (s)
    real(default) :: s
    type(vector4_t), intent(in) :: p
    integer, intent(in) :: e
    s = p * p
    if (e /= 2) then
       if (mod(e, 2) == 0) then
          s = s**(e / 2)
       else if (s >= 0) then
          s = sqrt(s)**e
       else
          s = -(sqrt(abs(s))**e)
       end if
    end if
  end function power_vector4

  elemental module function negate_vector4 (p) result (q)
    type(vector4_t) :: q
    type(vector4_t), intent(in) :: p
    integer :: i
    do i = 0, 3
       if (abs (p%p(i)) < eps0) then
          q%p(i) = 0
       else
          q%p(i) = -p%p(i)
       end if
    end do
  end function negate_vector4

  pure module function sum_vector4 (p) result (q)
    type(vector4_t) :: q
    type(vector4_t), dimension(:), intent(in) :: p
    integer :: i
    do i = 0, 3
       q%p(i) = sum (p%p(i))
    end do
  end function sum_vector4

  pure module function sum_vector4_mask (p, mask) result (q)
    type(vector4_t) :: q
    type(vector4_t), dimension(:), intent(in) :: p
    logical, dimension(:), intent(in) :: mask
    integer :: i
    do i = 0, 3
       q%p(i) = sum (p%p(i), mask=mask)
    end do
  end function sum_vector4_mask

  module subroutine vector4_set_component (p, k, c)
    type(vector4_t), intent(inout) :: p
    integer, intent(in) :: k
    real(default), intent(in) :: c
    p%p(k) = c
  end subroutine vector4_set_component

  elemental module function vector4_get_component (p, k) result (c)
    real(default) :: c
    type(vector4_t), intent(in) :: p
    integer, intent(in) :: k
    c = p%p(k)
  end function vector4_get_component

  pure module function vector4_get_components (p) result (a)
    real(default), dimension(0:3) :: a
    type(vector4_t), intent(in) :: p
    a = p%p
  end function vector4_get_components

  elemental module function vector4_get_space_part (p) result (q)
    type(vector3_t) :: q
    type(vector4_t), intent(in) :: p
    q%p = p%p(1:)
  end function vector4_get_space_part

  elemental module function vector4_get_direction (p) result (q)
    type(vector3_t) :: q
    type(vector4_t), intent(in) :: p
    real(default) :: qq
    q%p = p%p(1:)
    qq = q**1
    if (abs(qq) > eps0) then
       q%p = q%p / qq
    else
       q%p = 0
    end if
  end function vector4_get_direction

  elemental module subroutine vector4_invert_direction (p)
    type(vector4_t), intent(inout) :: p
    p%p(1:3) = -p%p(1:3)
  end subroutine vector4_invert_direction

  pure module subroutine array_from_vector4_1 (a, p)
    real(default), dimension(:), intent(out) :: a
    type(vector4_t), intent(in) :: p
    a = p%p
  end subroutine array_from_vector4_1

  pure module subroutine array_from_vector4_2 (a, p)
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), dimension(:,:), intent(out) :: a
    integer :: i
    forall (i=1:size(p))
       a(:,i) = p(i)%p
    end forall
  end subroutine array_from_vector4_2

  pure module subroutine array_from_vector3_1 (a, p)
    real(default), dimension(:), intent(out) :: a
    type(vector3_t), intent(in) :: p
    a = p%p
  end subroutine array_from_vector3_1

  pure module subroutine array_from_vector3_2 (a, p)
    type(vector3_t), dimension(:), intent(in) :: p
    real(default), dimension(:,:), intent(out) :: a
    integer :: i
    forall (i=1:size(p))
       a(:,i) = p(i)%p
    end forall
  end subroutine array_from_vector3_2

  pure module subroutine vector4_from_array (p, a)
    type(vector4_t), intent(out) :: p
    real(default), dimension(:), intent(in) :: a
    p%p(0:3) = a
  end subroutine vector4_from_array

  pure module subroutine vector3_from_array (p, a)
    type(vector3_t), intent(out) :: p
    real(default), dimension(:), intent(in) :: a
    p%p(1:3) = a
  end subroutine vector3_from_array

  pure module function vector4 (a) result (p)
    type(vector4_t) :: p
    real(default), intent(in), dimension(4) :: a
    p%p = a
  end function vector4

  pure module function vector4_to_pythia6 (vector4, m) result (p)
    real(double), dimension(1:5) :: p
    class(vector4_t), intent(in) :: vector4
    real(default), intent(in), optional :: m
    p(1:3) = vector4%p(1:3)
    p(4) = vector4%p(0)
    if (present (m)) then
       p(5) = m
    else
       p(5) = vector4 ** 1
    end if
  end function vector4_to_pythia6

  pure module subroutine vector4_from_c_prt (p, c_prt)
    type(vector4_t), intent(out) :: p
    type(c_prt_t), intent(in) :: c_prt
    p%p(0) = c_prt%pe
    p%p(1) = c_prt%px
    p%p(2) = c_prt%py
    p%p(3) = c_prt%pz
  end subroutine vector4_from_c_prt

  pure module subroutine c_prt_from_vector4 (c_prt, p)
    type(c_prt_t), intent(out) :: c_prt
    type(vector4_t), intent(in) :: p
    c_prt%pe = p%p(0)
    c_prt%px = p%p(1)
    c_prt%py = p%p(2)
    c_prt%pz = p%p(3)
    c_prt%p2 = p ** 2
  end subroutine c_prt_from_vector4

  elemental module function vector4_to_c_prt (p, p2) result (c_prt)
    type(c_prt_t) :: c_prt
    type(vector4_t), intent(in) :: p
    real(default), intent(in), optional :: p2
    c_prt%pe = p%p(0)
    c_prt%px = p%p(1)
    c_prt%py = p%p(2)
    c_prt%pz = p%p(3)
    if (present (p2)) then
       c_prt%p2 = p2
    else
       c_prt%p2 = p ** 2
    end if
  end function vector4_to_c_prt

  elemental module function vector3_azimuthal_angle (p) result (phi)
    real(default) :: phi
    type(vector3_t), intent(in) :: p
    if (any (abs (p%p(1:2)) > 0)) then
       phi = atan2(p%p(2), p%p(1))
       if (phi < 0) phi = phi + twopi
    else
       phi = 0
    end if
  end function vector3_azimuthal_angle
  elemental module function vector4_azimuthal_angle (p) result (phi)
    real(default) :: phi
    type(vector4_t), intent(in) :: p
    phi = vector3_azimuthal_angle (space_part (p))
  end function vector4_azimuthal_angle

  elemental module function vector3_azimuthal_angle_deg (p) result (phi)
    real(default) :: phi
    type(vector3_t), intent(in) :: p
    phi = vector3_azimuthal_angle (p) / degree
  end function vector3_azimuthal_angle_deg
  elemental module function vector4_azimuthal_angle_deg (p) result (phi)
    real(default) :: phi
    type(vector4_t), intent(in) :: p
    phi = vector4_azimuthal_angle (p) / degree
  end function vector4_azimuthal_angle_deg

  elemental module function vector3_azimuthal_distance (p, q) result (dphi)
    real(default) :: dphi
    type(vector3_t), intent(in) :: p,q
    dphi = vector3_azimuthal_angle (q) - vector3_azimuthal_angle (p)
    if (dphi <= -pi) then
       dphi = dphi + twopi
    else if (dphi > pi) then
       dphi = dphi - twopi
    end if
  end function vector3_azimuthal_distance
  elemental module function vector4_azimuthal_distance (p, q) result (dphi)
    real(default) :: dphi
    type(vector4_t), intent(in) :: p,q
    dphi = vector3_azimuthal_distance &
         (space_part (p), space_part (q))
  end function vector4_azimuthal_distance

  elemental module function vector3_azimuthal_distance_deg (p, q) result (dphi)
    real(default) :: dphi
    type(vector3_t), intent(in) :: p,q
    dphi = vector3_azimuthal_distance (p, q) / degree
  end function vector3_azimuthal_distance_deg
  elemental module function vector4_azimuthal_distance_deg (p, q) result (dphi)
    real(default) :: dphi
    type(vector4_t), intent(in) :: p,q
    dphi = vector4_azimuthal_distance (p, q) / degree
  end function vector4_azimuthal_distance_deg

  elemental module function polar_angle_vector3 (p) result (theta)
    real(default) :: theta
    type(vector3_t), intent(in) :: p
    if (any (abs (p%p) > 0)) then
       theta = atan2 (sqrt(p%p(1)**2 + p%p(2)**2), p%p(3))
    else
       theta = 0
    end if
  end function polar_angle_vector3
  elemental module function polar_angle_vector4 (p) result (theta)
    real(default) :: theta
    type(vector4_t), intent(in) :: p
    theta = polar_angle (space_part (p))
  end function polar_angle_vector4

  elemental module function polar_angle_ct_vector3 (p) result (ct)
    real(default) :: ct
    type(vector3_t), intent(in) :: p
    if (any (abs (p%p) > 0)) then
       ct = p%p(3) / p**1
    else
       ct = 1
    end if
  end function polar_angle_ct_vector3
  elemental module function polar_angle_ct_vector4 (p) result (ct)
    real(default) :: ct
    type(vector4_t), intent(in) :: p
    ct = polar_angle_ct (space_part (p))
  end function polar_angle_ct_vector4

  elemental module function polar_angle_deg_vector3 (p) result (theta)
    real(default) :: theta
    type(vector3_t), intent(in) :: p
    theta = polar_angle (p) / degree
  end function polar_angle_deg_vector3
  elemental module function polar_angle_deg_vector4 (p) result (theta)
    real(default) :: theta
    type(vector4_t), intent(in) :: p
    theta = polar_angle (p) / degree
  end function polar_angle_deg_vector4

  elemental module function enclosed_angle_vector3 (p, q) result (theta)
    real(default) :: theta
    type(vector3_t), intent(in) :: p, q
    theta = acos (enclosed_angle_ct (p, q))
  end function enclosed_angle_vector3
  elemental module function enclosed_angle_vector4 (p, q) result (theta)
    real(default) :: theta
    type(vector4_t), intent(in) :: p, q
    theta = enclosed_angle (space_part (p), space_part (q))
  end function enclosed_angle_vector4

  elemental module function enclosed_angle_ct_vector3 (p, q) result (ct)
    real(default) :: ct
    type(vector3_t), intent(in) :: p, q
    if (any (abs (p%p) > 0) .and. any (abs (q%p) > 0)) then
       ct = p*q / (p**1 * q**1)
       if (ct>1) then
          ct = 1
       else if (ct<-1) then
          ct = -1
       end if
    else
       ct = 1
    end if
  end function enclosed_angle_ct_vector3
  elemental module function enclosed_angle_ct_vector4 (p, q) result (ct)
    real(default) :: ct
    type(vector4_t), intent(in) :: p, q
    ct = enclosed_angle_ct (space_part (p), space_part (q))
  end function enclosed_angle_ct_vector4

  elemental module function enclosed_angle_deg_vector3 (p, q) result (theta)
    real(default) :: theta
    type(vector3_t), intent(in) :: p, q
    theta = enclosed_angle (p, q) / degree
  end function enclosed_angle_deg_vector3
  elemental module function enclosed_angle_deg_vector4 (p, q) result (theta)
    real(default) :: theta
    type(vector4_t), intent(in) :: p, q
    theta = enclosed_angle (p, q) / degree
  end function enclosed_angle_deg_vector4

  elemental module function enclosed_angle_rest_frame_vector4 (p, q) result (theta)
    type(vector4_t), intent(in) :: p, q
    real(default) :: theta
    theta = acos (enclosed_angle_ct_rest_frame (p, q))
  end function enclosed_angle_rest_frame_vector4
  elemental module function enclosed_angle_ct_rest_frame_vector4 (p, q) result (ct)
    type(vector4_t), intent(in) :: p, q
    real(default) :: ct
    if (invariant_mass(q) > 0) then
       ct = enclosed_angle_ct ( &
            space_part (boost(-q, invariant_mass (q)) * p), &
            space_part (q))
    else
       ct = 1
    end if
  end function enclosed_angle_ct_rest_frame_vector4
  elemental module function enclosed_angle_deg_rest_frame_vector4 (p, q) &
       result (theta)
    type(vector4_t), intent(in) :: p, q
    real(default) :: theta
    theta = enclosed_angle_rest_frame (p, q) / degree
  end function enclosed_angle_deg_rest_frame_vector4

  elemental module function transverse_part_vector4_beam_axis (p) result (pT)
    real(default) :: pT
    type(vector4_t), intent(in) :: p
    pT = sqrt(p%p(1)**2 + p%p(2)**2)
  end function transverse_part_vector4_beam_axis

  elemental module function transverse_part_vector4_vector4 (p1, p2) result (pT)
    real(default) :: pT
    type(vector4_t), intent(in) :: p1, p2
    real(default) :: p1_norm, p2_norm, p1p2, pT2
    p1_norm = space_part_norm(p1)**2
    p2_norm = space_part_norm(p2)**2
!    p1p2 = p1%p(1:3)*p2%p(1:3)
    p1p2 = vector4_get_space_part(p1) * vector4_get_space_part(p2)
    pT2 = (p1_norm*p2_norm - p1p2)/p1_norm
    pT = sqrt (pT2)
  end function transverse_part_vector4_vector4

  elemental module function longitudinal_part_vector4 (p) result (pL)
    real(default) :: pL
    type(vector4_t), intent(in) :: p
    pL = p%p(3)
  end function longitudinal_part_vector4

  elemental module function space_part_norm_vector4 (p) result (p3)
    real(default) :: p3
    type(vector4_t), intent(in) :: p
    p3 = sqrt (p%p(1)**2 + p%p(2)**2 + p%p(3)**2)
  end function space_part_norm_vector4

  elemental module function energy_vector4 (p) result (E)
    real(default) :: E
    type(vector4_t), intent(in) :: p
    E = p%p(0)
  end function energy_vector4

  elemental module function energy_vector3 (p, mass) result (E)
    real(default) :: E
    type(vector3_t), intent(in) :: p
    real(default), intent(in), optional :: mass
    if (present (mass)) then
       E = sqrt (p**2 + mass**2)
    else
       E = p**1
    end if
  end function energy_vector3

  elemental module function energy_real (p, mass) result (E)
    real(default) :: E
    real(default), intent(in) :: p
    real(default), intent(in), optional :: mass
    if (present (mass)) then
       E = sqrt (p**2 + mass**2)
    else
       E = abs (p)
    end if
  end function energy_real

  elemental module function invariant_mass_vector4 (p) result (m)
    real(default) :: m
    type(vector4_t), intent(in) :: p
    real(default) :: msq
    msq = p*p
    if (msq >= 0) then
       m = sqrt (msq)
    else
       m = - sqrt (abs (msq))
    end if
  end function invariant_mass_vector4
  elemental module function invariant_mass_squared_vector4 (p) result (msq)
    real(default) :: msq
    type(vector4_t), intent(in) :: p
    msq = p*p
  end function invariant_mass_squared_vector4

  elemental module function transverse_mass_vector4 (p) result (m)
    real(default) :: m
    type(vector4_t), intent(in) :: p
    real(default) :: msq
    msq = p%p(0)**2 - p%p(1)**2 - p%p(2)**2
    if (msq >= 0) then
       m = sqrt (msq)
    else
       m = - sqrt (abs (msq))
    end if
  end function transverse_mass_vector4

  elemental module function rapidity_vector4 (p) result (y)
    real(default) :: y
    type(vector4_t), intent(in) :: p
    y = .5 * log( (energy (p) + longitudinal_part (p)) &
         &       /(energy (p) - longitudinal_part (p)))
  end function rapidity_vector4

  elemental module function pseudorapidity_vector4 (p) result (eta)
    real(default) :: eta
    type(vector4_t), intent(in) :: p
    eta = -log( tan (.5 * polar_angle (p)))
  end function pseudorapidity_vector4

  elemental module function rapidity_distance_vector4 (p, q) result (dy)
    type(vector4_t), intent(in) :: p, q
    real(default) :: dy
    dy = rapidity (q) - rapidity (p)
  end function rapidity_distance_vector4

  elemental module function pseudorapidity_distance_vector4 (p, q) result (deta)
    real(default) :: deta
    type(vector4_t), intent(in) :: p, q
    deta = pseudorapidity (q) - pseudorapidity (p)
  end function pseudorapidity_distance_vector4

  elemental module function eta_phi_distance_vector4 (p, q) result (dr)
    type(vector4_t), intent(in) :: p, q
    real(default) :: dr
    dr = sqrt ( &
         pseudorapidity_distance (p, q)**2 &
         + azimuthal_distance (p, q)**2)
  end function eta_phi_distance_vector4

  module subroutine lorentz_transformation_write (L, unit, testflag, ultra)
    class(lorentz_transformation_t), intent(in) :: L
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag, ultra
    integer :: u, i
    logical :: ult
    character(len=7) :: fmt
    ult = .false.; if (present (ultra)) ult = ultra
    if (ult) then
       call pac_fmt (fmt, FMT_19, FMT_11, ultra)
    else
       call pac_fmt (fmt, FMT_19, FMT_13, testflag)
    end if
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A,3(1x," // fmt // "))")  "L00 = ", L%L(0,0)
    write (u, "(1x,A,3(1x," // fmt // "))")  "L0j = ", L%L(0,1:3)
    do i = 1, 3
       write (u, "(1x,A,I0,A,3(1x," // fmt // "))")  &
            "L", i, "0 = ", L%L(i,0)
       write (u, "(1x,A,I0,A,3(1x," // fmt // "))")  &
            "L", i, "j = ", L%L(i,1:3)
    end do
  end subroutine lorentz_transformation_write

  pure module function lorentz_transformation_get_components (L) result (a)
    type(lorentz_transformation_t), intent(in) :: L
    real(default), dimension(0:3,0:3) :: a
    a = L%L
  end function lorentz_transformation_get_components

  elemental module function lorentz_transformation_inverse (L) result (IL)
    type(lorentz_transformation_t) :: IL
    type(lorentz_transformation_t), intent(in) :: L
    IL%L(0,0) = L%L(0,0)
    IL%L(0,1:) = -L%L(1:,0)
    IL%L(1:,0) = -L%L(0,1:)
    IL%L(1:,1:) = transpose(L%L(1:,1:))
  end function lorentz_transformation_inverse

  module function create_orthogonal (p_in) result (p_out)
    type(vector3_t), intent(in) :: p_in
    type(vector3_t) :: p_out
    real(default) :: ab
    ab = sqrt (p_in%p(1)**2 + p_in%p(2)**2)
    if (abs (ab) < eps0) then
      p_out%p(1) = 1
      p_out%p(2) = 0
      p_out%p(3) = 0
    else
      p_out%p(1) = p_in%p(2)
      p_out%p(2) = -p_in%p(1)
      p_out%p(3) = 0
      p_out = p_out / ab
    end if
  end function create_orthogonal

  module function create_unit_vector (p_in) result (p_out)
    type(vector4_t), intent(in) :: p_in
    type(vector3_t) :: p_out
    p_out%p = p_in%p(1:3) / space_part_norm (p_in)
  end function create_unit_vector

  module function normalize(p) result (p_norm)
    type(vector3_t) :: p_norm
    type(vector3_t), intent(in) :: p
    real(default) :: abs
    abs = sqrt (p%p(1)**2 + p%p(2)**2 + p%p(3)**2)
    p_norm = p / abs
  end function normalize

  pure module function compute_resonance_mass (p, i_res_born, i_gluon) result (m)
    real(default) :: m
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in), dimension(:) :: i_res_born
    integer, intent(in), optional :: i_gluon
    type(vector4_t) :: p_res
    p_res = get_resonance_momentum (p, i_res_born, i_gluon)
    m = p_res**1
  end function compute_resonance_mass

  pure module function get_resonance_momentum &
       (p, i_res_born, i_gluon) result (p_res)
    type(vector4_t) :: p_res
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in), dimension(:) :: i_res_born
    integer, intent(in), optional :: i_gluon
    integer :: i
    p_res = vector4_null
    do i = 1, size (i_res_born)
       p_res = p_res + p (i_res_born(i))
    end do
    if (present (i_gluon)) p_res = p_res + p (i_gluon)
  end function get_resonance_momentum

  module function create_two_particle_decay (s, p1, p2) result (p_rest)
    type(vector4_t), dimension(3) :: p_rest
    real(default), intent(in) :: s
    type(vector4_t), intent(in) :: p1, p2
    real(default) :: m1_sq, m2_sq
    real(default) :: E1, E2, p
    m1_sq = p1**2; m2_sq = p2**2
    p = sqrt (lambda (s, m1_sq, m2_sq)) / (two * sqrt (s))
    E1 = sqrt (m1_sq + p**2); E2 = sqrt (m2_sq + p**2)
    p_rest(1)%p = [sqrt (s), zero, zero, zero]
    p_rest(2)%p(0) = E1
    p_rest(2)%p(1:3) = p * p1%p(1:3) / space_part_norm (p1)
    p_rest(3)%p(0) = E2; p_rest(3)%p(1:3) = -p_rest(2)%p(1:3)
  end function create_two_particle_decay

  module function create_three_particle_decay (p1, p2, p3) result (p_rest)
    type(vector4_t), dimension(4) :: p_rest
    type(vector4_t), intent(in) :: p1, p2, p3
    real(default) :: E1, E2, E3
    real(default) :: pr1, pr2, pr3
    real(default) :: s, s1, s2, s3
    real(default) :: m1_sq, m2_sq, m3_sq
    real(default) :: cos_theta_12
    type(vector3_t) :: v3_unit
    type(lorentz_transformation_t) :: rot
    m1_sq = p1**2
    m2_sq = p2**2
    m3_sq = p3**2
    s1 = (p1 + p2)**2
    s2 = (p2 + p3)**2
    s3 = (p3 + p1)**2
    s = s1 + s2 + s3 - m1_sq - m2_sq - m3_sq
    E1 = (s - s2 - m1_sq) / (two * sqrt (s2))
    E2 = (s2 + m2_sq - m3_sq) / (two * sqrt (s2))
    E3 = (s2 + m3_sq - m2_sq) / (two * sqrt (s2))
    pr1 = sqrt (lambda (s, s2, m1_sq)) / (two * sqrt (s2))
    pr2 = sqrt (lambda (s2, m2_sq, m3_sq)) / (two * sqrt(s2))
    pr3 = pr2
    cos_theta_12 = ((s - s2 - m1_sq) * (s2 + m2_sq - m3_sq) + two * s2 * (m1_sq + m2_sq - s1)) / &
         sqrt (lambda (s, s2, m1_sq) * lambda (s2, m2_sq, m3_sq))
    v3_unit%p = [zero, zero, one]
    p_rest(1)%p(0) = E1
    p_rest(1)%p(1:3) = v3_unit%p * pr1
    p_rest(2)%p(0) = E2
    p_rest(2)%p(1:3) = v3_unit%p * pr2
    p_rest(3)%p(0) = E3
    p_rest(3)%p(1:3) = v3_unit%p * pr3
    p_rest(4)%p(0) = (s + s2 - m1_sq) / (2 * sqrt (s2))
    p_rest(4)%p(1:3) = - p_rest(1)%p(1:3)
    rot = rotation (cos_theta_12, sqrt (one - cos_theta_12**2), 2)
    p_rest(2) = rot * p_rest(2)
    p_rest(3)%p(1:3) = - p_rest(2)%p(1:3)
  end function create_three_particle_decay

  recursive module subroutine generate_on_shell_decay (p_dec, &
      p_in, p_out, i_real, msq_in, jac, evaluate_special)
    type(vector4_t), intent(in) :: p_dec
    type(vector4_t), intent(in), dimension(:) :: p_in
    type(vector4_t), intent(inout), dimension(:) :: p_out
    integer, intent(in) :: i_real
    real(default), intent(in), optional :: msq_in
    real(default), intent(inout), optional :: jac
    procedure(evaluate_one_to_two_splitting_special), intent(in), &
          pointer, optional :: evaluate_special
    type(vector4_t) :: p_dec_new
    integer :: n_recoil
    n_recoil = size (p_in) - 1
    if (n_recoil > 1) then
       if (present (evaluate_special)) then
          call evaluate_special (p_dec, p_in(1), sum (p_in (2 : n_recoil + 1)), &
               p_out(i_real), p_dec_new)
          call generate_on_shell_decay (p_dec_new, p_in (2 : ), p_out, &
               i_real + 1, msq_in, jac, evaluate_special)
       else
          call evaluate_one_to_two_splitting (p_dec, p_in(1), &
               sum (p_in (2 : n_recoil + 1)), p_out(i_real), p_dec_new, msq_in, jac)
          call generate_on_shell_decay (p_dec_new, p_in (2 : ), p_out, &
               i_real + 1, msq_in, jac)
       end if
    else
       call evaluate_one_to_two_splitting (p_dec, p_in(1), p_in(2), &
            p_out(i_real), p_out(i_real + 1), msq_in, jac)
    end if

  end subroutine generate_on_shell_decay

  subroutine evaluate_one_to_two_splitting (p_origin, &
      p1_in, p2_in, p1_out, p2_out, msq_in, jac)
    type(vector4_t), intent(in) :: p_origin
    type(vector4_t), intent(in) :: p1_in, p2_in
    type(vector4_t), intent(inout) :: p1_out, p2_out
    real(default), intent(in), optional :: msq_in
    real(default), intent(inout), optional :: jac
    type(lorentz_transformation_t) :: L
    type(vector4_t) :: p1_rest, p2_rest
    real(default) :: m, msq, msq1, msq2
    real(default) :: E1, E2, p
    real(default) :: lda, rlda_soft

    call get_rest_frame (p1_in, p2_in, p1_rest, p2_rest)

    msq = p_origin**2; m = sqrt(msq)
    msq1 = p1_in**2; msq2 = p2_in**2

    lda = lambda (msq, msq1, msq2)
    if (lda < zero) then
       print *, 'Encountered lambda < 0 in 1 -> 2 splitting! '
       print *, 'lda: ', lda
       print *, 'm: ', m, 'msq: ', msq
       print *, 'm1: ', sqrt (msq1), 'msq1: ', msq1
       print *, 'm2: ', sqrt (msq2), 'msq2: ', msq2
       stop
    end if
    p = sqrt (lda) / (two * m)

    E1 = sqrt (msq1 + p**2)
    E2 = sqrt (msq2 + p**2)

    p1_out = shift_momentum (p1_rest, E1, p)
    p2_out = shift_momentum (p2_rest, E2, p)

    L = boost (p_origin, p_origin**1)
    p1_out = L  * p1_out
    p2_out = L  * p2_out

    if (present (jac) .and. present (msq_in)) then
       jac = jac * sqrt(lda) / msq
       rlda_soft = sqrt (lambda (msq_in, msq1, msq2))
       !!! We have to undo the Jacobian which has already been
       !!! supplied by the Born phase space.
       jac = jac * msq_in / rlda_soft
    end if

  contains

   subroutine get_rest_frame (p1_in, p2_in, p1_out, p2_out)
     type(vector4_t), intent(in) :: p1_in, p2_in
     type(vector4_t), intent(out) :: p1_out, p2_out
     type(lorentz_transformation_t) :: L
     L = inverse (boost (p1_in + p2_in, (p1_in + p2_in)**1))
     p1_out = L * p1_in; p2_out = L * p2_in
   end subroutine get_rest_frame

   function shift_momentum (p_in, E, p) result (p_out)
     type(vector4_t) :: p_out
     type(vector4_t), intent(in) :: p_in
     real(default), intent(in) :: E, p
     type(vector3_t) :: vec
     vec%p(1:3) = p_in%p(1:3) / space_part_norm (p_in)
     p_out = vector4_moving (E, p * vec)
   end function shift_momentum

  end subroutine evaluate_one_to_two_splitting

  elemental module function boost_from_rest_frame (p, m) result (L)
    type(lorentz_transformation_t) :: L
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: m
    L = boost_from_rest_frame_vector3 (space_part (p), m)
  end function boost_from_rest_frame
  elemental module function boost_from_rest_frame_vector3 (p, m) result (L)
    type(lorentz_transformation_t) :: L
    type(vector3_t), intent(in) :: p
    real(default), intent(in) :: m
    type(vector3_t) :: beta_gamma
    real(default) :: bg2, g, c
    integer :: i,j
    if (m > eps0) then
       beta_gamma = p / m
       bg2 = beta_gamma**2
    else
       bg2 = 0
       L = identity
       return
    end if
    if (bg2 > eps0) then
       g = sqrt(1 + bg2);  c = (g-1)/bg2
    else
       g = one + bg2 / two
       c = one / two
    end if
    L%L(0,0)  = g
    L%L(0,1:) = beta_gamma%p
    L%L(1:,0) = L%L(0,1:)
    do i=1,3
       do j=1,3
          L%L(i,j) = delta_three(i,j) + c*beta_gamma%p(i)*beta_gamma%p(j)
       end do
    end do
  end function boost_from_rest_frame_vector3
  elemental module function boost_canonical (beta_gamma, k) result (L)
    type(lorentz_transformation_t) :: L
    real(default), intent(in) :: beta_gamma
    integer, intent(in) :: k
    real(default) :: g
    g = sqrt(1 + beta_gamma**2)
    L = identity
    L%L(0,0) = g
    L%L(0,k) = beta_gamma
    L%L(k,0) = L%L(0,k)
    L%L(k,k) = L%L(0,0)
  end function boost_canonical
  elemental module function boost_generic (beta_gamma, axis) result (L)
    type(lorentz_transformation_t) :: L
    real(default), intent(in) :: beta_gamma
    type(vector3_t), intent(in) :: axis
    if (any (abs (axis%p) > 0)) then
       L = boost_from_rest_frame_vector3 (beta_gamma * axis, axis**1)
    else
       L = identity
    end if
  end function boost_generic

  elemental module function rotation_generic_cs (cp, sp, axis) result (R)
    type(lorentz_transformation_t) :: R
    real(default), intent(in) :: cp, sp
    type(vector3_t), intent(in) :: axis
    integer :: i,j
    R = identity
    do i=1,3
       do j=1,3
          R%L(i,j) = cp*delta_three(i,j) + (1-cp)*axis%p(i)*axis%p(j)  &
               &   - sp*dot_product(epsilon_three(i,j,:), axis%p)
       end do
    end do
  end function rotation_generic_cs
  elemental module function rotation_generic (axis) result (R)
    type(lorentz_transformation_t) :: R
    type(vector3_t), intent(in) :: axis
    real(default) :: phi
    if (any (abs(axis%p) > 0)) then
       phi = abs(axis**1)
       R = rotation_generic_cs (cos(phi), sin(phi), axis/phi)
    else
       R = identity
    end if
  end function rotation_generic
  elemental module function rotation_canonical_cs (cp, sp, k) result (R)
    type(lorentz_transformation_t) :: R
    real(default), intent(in) :: cp, sp
    integer, intent(in) :: k
    integer :: i,j
    R = identity
    do i=1,3
       do j=1,3
          R%L(i,j) = -sp*epsilon_three(i,j,k)
       end do
       R%L(i,i) = cp
    end do
    R%L(k,k) = 1
  end function rotation_canonical_cs
  elemental module function rotation_canonical (phi, k) result (R)
    type(lorentz_transformation_t) :: R
    real(default), intent(in) :: phi
    integer, intent(in) :: k
    R = rotation_canonical_cs(cos(phi), sin(phi), k)
  end function rotation_canonical
  elemental module function rotation_to_2nd_generic (p, q) result (R)
    type(lorentz_transformation_t) :: R
    type(vector3_t), intent(in) :: p, q
    type(vector3_t) :: a, b, ab
    real(default) :: ct, st
    if (any (abs (p%p) > 0) .and. any (abs (q%p) > 0)) then
       a = direction (p)
       b = direction (q)
       ab = cross_product(a,b)
       ct = a * b;  st = ab**1
       if (abs(st) > eps0) then
          R = rotation_generic_cs (ct, st, ab / st)
       else if (ct < 0) then
          R = space_reflection
       else
          R = identity
       end if
    else
       R = identity
    end if
  end function rotation_to_2nd_generic
  elemental module function rotation_to_2nd_canonical (k, p) result (R)
    type(lorentz_transformation_t) :: R
    integer, intent(in) :: k
    type(vector3_t), intent(in) :: p
    type(vector3_t) :: b, ab
    real(default) :: ct, st
    integer :: i, j
    if (any (abs (p%p) > 0)) then
       b = direction (p)
       ab%p = 0
       do i = 1, 3
          do j = 1, 3
             ab%p(j) = ab%p(j) + b%p(i) * epsilon_three(i,j,k)
          end do
       end do
       ct = b%p(k);  st = ab**1
       if (abs(st) > eps0) then
          R = rotation_generic_cs (ct, st, ab / st)
       else if (ct < 0) then
          R = space_reflection
       else
          R = identity
       end if
    else
       R = identity
    end if
  end function rotation_to_2nd_canonical

  elemental module function transformation_rec_generic (axis, p1, p2, m) result (L)
    type(vector3_t), intent(in) :: axis
    type(vector4_t), intent(in) :: p1, p2
    real(default), intent(in) :: m
    type(lorentz_transformation_t) :: L
    L = boost (p1 + p2, m)
    L = L * rotation_to_2nd (axis, space_part (inverse (L) * p1))
  end function transformation_rec_generic
  elemental module function transformation_rec_canonical (k, p1, p2, m) result (L)
    integer, intent(in) :: k
    type(vector4_t), intent(in) :: p1, p2
    real(default), intent(in) :: m
    type(lorentz_transformation_t) :: L
    L = boost (p1 + p2, m)
    L = L * rotation_to_2nd (k, space_part (inverse (L) * p1))
  end function transformation_rec_canonical
  elemental module function prod_LT_vector4 (L, p) result (np)
    type(vector4_t) :: np
    type(lorentz_transformation_t), intent(in) :: L
    type(vector4_t), intent(in) :: p
    np%p = matmul (L%L, p%p)
  end function prod_LT_vector4
  elemental module function prod_LT_LT (L1, L2) result (NL)
    type(lorentz_transformation_t) :: NL
    type(lorentz_transformation_t), intent(in) :: L1,L2
    NL%L = matmul (L1%L, L2%L)
  end function prod_LT_LT
  elemental module function prod_vector4_LT (p, L) result (np)
    type(vector4_t) :: np
    type(vector4_t), intent(in) :: p
    type(lorentz_transformation_t), intent(in) :: L
    np%p = matmul (p%p, L%L)
  end function prod_vector4_LT

  elemental module function LT_compose_r3_r2_b3 &
       (cp, sp, ct, st, beta_gamma) result (L)
    type(lorentz_transformation_t) :: L
    real(default), intent(in) :: cp, sp, ct, st, beta_gamma
    real(default) :: gamma
    if (abs(beta_gamma) < eps0) then
       L%L(0,0)  = 1
       L%L(1:,0) = 0
       L%L(0,1:) = 0
       L%L(1,1:) = [  ct*cp, -ct*sp, st ]
       L%L(2,1:) = [     sp,     cp,  zero ]
       L%L(3,1:) = [ -st*cp,  st*sp, ct ]
    else
       gamma = sqrt(1 + beta_gamma**2)
       L%L(0,0)  = gamma
       L%L(1,0)  = 0
       L%L(2,0)  = 0
       L%L(3,0)  = beta_gamma
       L%L(0,1:) = beta_gamma * [ -st*cp,  st*sp, ct ]
       L%L(1,1:) =              [  ct*cp, -ct*sp, st ]
       L%L(2,1:) =              [     sp,     cp, zero ]
       L%L(3,1:) = gamma      * [ -st*cp,  st*sp, ct ]
    end if
  end function LT_compose_r3_r2_b3

  elemental module function LT_compose_r2_r3_b3 &
       (ct, st, cp, sp, beta_gamma) result (L)
    type(lorentz_transformation_t) :: L
    real(default), intent(in) :: ct, st, cp, sp, beta_gamma
    real(default) :: gamma
    if (abs(beta_gamma) < eps0) then
       L%L(0,0)  = 1
       L%L(1:,0) = 0
       L%L(0,1:) = 0
       L%L(1,1:) = [  ct*cp,    -sp,     st*cp ]
       L%L(2,1:) = [  ct*sp,     cp,     st*sp ]
       L%L(3,1:) = [ -st   ,   zero,     ct    ]
    else
       gamma = sqrt(1 + beta_gamma**2)
       L%L(0,0)  = gamma
       L%L(1,0)  = 0
       L%L(2,0)  = 0
       L%L(3,0)  = beta_gamma
       L%L(0,1:) = beta_gamma * [ -st   ,   zero,     ct    ]
       L%L(1,1:) =              [  ct*cp,    -sp,     st*cp ]
       L%L(2,1:) =              [  ct*sp,     cp,     st*sp ]
       L%L(3,1:) = gamma      * [ -st   ,   zero,     ct    ]
    end if
  end function LT_compose_r2_r3_b3

  elemental module function axis_from_p_r3_r2_b3 &
       (p, cp, sp, ct, st, beta_gamma) result (n)
    type(vector3_t) :: n
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: cp, sp, ct, st, beta_gamma
    real(default) :: gamma, px, py
    px = cp * p%p(1) - sp * p%p(2)
    py = sp * p%p(1) + cp * p%p(2)
    n%p(1) =  ct * px + st * p%p(3)
    n%p(2) = py
    n%p(3) = -st * px + ct * p%p(3)
    if (abs(beta_gamma) > eps0) then
       gamma = sqrt(1 + beta_gamma**2)
       n%p(3) = n%p(3) * gamma + p%p(0) * beta_gamma
    end if
  end function axis_from_p_r3_r2_b3

  elemental module function axis_from_p_b3 (p, beta_gamma) result (n)
    type(vector3_t) :: n
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: beta_gamma
    real(default) :: gamma
    n%p = p%p(1:3)
    if (abs(beta_gamma) > eps0) then
       gamma = sqrt(1 + beta_gamma**2)
       n%p(3) = n%p(3) * gamma + p%p(0) * beta_gamma
    end if
  end function axis_from_p_b3

  elemental module function lambda (m1sq, m2sq, m3sq)
    real(default) :: lambda
    real(default), intent(in) :: m1sq, m2sq, m3sq
    lambda = (m1sq - m2sq - m3sq)**2 - 4*m2sq*m3sq
  end function lambda

  module function colliding_momenta (sqrts, m, p_cm) result (p)
    type(vector4_t), dimension(2) :: p
    real(default), intent(in) :: sqrts
    real(default), dimension(2), intent(in), optional :: m
    real(default), intent(in), optional :: p_cm
    real(default), dimension(2) :: dmsq
    real(default) :: ch, sh
    real(default), dimension(2) :: E0, p0
    integer, dimension(2), parameter :: sgn = [1, -1]
    if (abs(sqrts) < eps0) then
       call msg_fatal (" Colliding beams: sqrts is zero (please set sqrts)")
       p = vector4_null;  return
    else if (sqrts <= 0) then
       call msg_fatal (" Colliding beams: sqrts is negative")
       p = vector4_null;  return
    end if
    if (present (m)) then
       dmsq = sgn * (m(1)**2-m(2)**2)
       E0 = (sqrts + dmsq/sqrts) / 2
       if (any (E0 < m)) then
          call msg_fatal &
               (" Colliding beams: beam energy is less than particle mass")
          p = vector4_null;  return
       end if
       p0 = sgn * sqrt (E0**2 - m**2)
    else
       E0 = sqrts / 2
       p0 = sgn * E0
    end if
    if (present (p_cm)) then
       sh = p_cm / sqrts
       ch = sqrt (1 + sh**2)
       p = vector4_moving (E0 * ch + p0 * sh, E0 * sh + p0 * ch, 3)
    else
       p = vector4_moving (E0, p0, 3)
    end if
  end function colliding_momenta

  elemental module subroutine pacify_vector3 (p, tolerance)
    type(vector3_t), intent(inout) :: p
    real(default), intent(in) :: tolerance
    where (abs (p%p) < tolerance)  p%p = zero
  end subroutine pacify_vector3

  elemental module subroutine pacify_vector4 (p, tolerance)
    type(vector4_t), intent(inout) :: p
    real(default), intent(in) :: tolerance
    where (abs (p%p) < tolerance)  p%p = zero
  end subroutine pacify_vector4

  elemental module subroutine pacify_LT (LT, tolerance)
    type(lorentz_transformation_t), intent(inout) :: LT
    real(default), intent(in) :: tolerance
    where (abs (LT%L) < tolerance)  LT%L = zero
  end subroutine pacify_LT

  module subroutine vector_set_reshuffle (p1, list, p2)
    type(vector4_t), intent(in), dimension(:), allocatable :: p1
    integer, intent(in), dimension(:), allocatable :: list
    type(vector4_t), intent(out), dimension(:), allocatable :: p2
    integer :: n, n_p
    n_p = size (p1)
    if (size (list) /= n_p) return
    allocate (p2 (n_p))
    do n = 1, n_p
      p2(n) = p1(list(n))
    end do
  end subroutine vector_set_reshuffle

  module function vector_set_is_cms (p, n_in) result (is_cms)
    logical :: is_cms
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in) :: n_in
    integer :: i
    type(vector4_t) :: p_sum
    p_sum%p = 0._default
    do i = 1, n_in
       p_sum = p_sum + p(i)
    end do
    is_cms = all (abs (p_sum%p(1:3)) < tiny_07)
  end function vector_set_is_cms

  module subroutine vector4_write_set (p, unit, show_mass, testflag, &
        check_conservation, ultra, n_in)
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_mass
    logical, intent(in), optional :: testflag, ultra
    logical, intent(in), optional :: check_conservation
    integer, intent(in), optional :: n_in
    logical :: extreme
    integer :: i, j
    real(default), dimension(0:3) :: p_tot
    character(len=7) :: fmt
    integer :: u
    logical :: yorn, is_test
    integer :: n
    extreme = .false.; if (present (ultra))  extreme = ultra
    is_test = .false.; if (present (testflag)) is_test = testflag
    u = given_output_unit (unit);  if (u < 0)  return
    n = 2; if (present (n_in)) n = n_in
    p_tot = 0
    yorn = .false.; if (present (check_conservation)) yorn = check_conservation
    do i = 1, size (p)
      if (yorn .and. i > n) then
         forall (j=0:3) p_tot(j) = p_tot(j) - p(i)%p(j)
      else
         forall (j=0:3) p_tot(j) = p_tot(j) + p(i)%p(j)
      end if
      call vector4_write (p(i), u, show_mass=show_mass, &
           testflag=testflag, ultra=ultra)
    end do
    if (extreme) then
       call pac_fmt (fmt, FMT_19, FMT_11, testflag)
    else
       call pac_fmt (fmt, FMT_19, FMT_15, testflag)
    end if
    if (is_test)  call pacify (p_tot, 1.E-9_default)
    if (.not. is_test) then
       write (u, "(A5)") 'Total: '
       write (u, "(1x,A,1x," // fmt // ")")    "E = ", p_tot(0)
       write (u, "(1x,A,3(1x," // fmt // "))") "P = ", p_tot(1:)
    end if
  end subroutine vector4_write_set

  module subroutine vector4_check_momentum_conservation (p, n_in, unit, &
     abs_smallness, rel_smallness, verbose)
    type(vector4_t), dimension(:), intent(in) :: p
    integer, intent(in) :: n_in
    integer, intent(in), optional :: unit
    real(default), intent(in), optional :: abs_smallness, rel_smallness
    logical, intent(in), optional :: verbose
    integer :: u, i
    type(vector4_t) :: psum_in, psum_out
    logical, dimension(0:3) :: p_diff
    logical :: verb
    u = given_output_unit (unit);  if (u < 0)  return
    verb = .false.; if (present (verbose)) verb = verbose
    psum_in = vector4_null
    do i = 1, n_in
       psum_in = psum_in + p(i)
    end do
    psum_out = vector4_null
    do i = n_in + 1, size (p)
       psum_out = psum_out + p(i)
    end do
    p_diff = vanishes (psum_in%p - psum_out%p, &
         abs_smallness = abs_smallness, rel_smallness = rel_smallness)
    if (.not. all (p_diff)) then
       call msg_warning ("Momentum conservation: FAIL", unit = u)
       if (verb) then
          write (u, "(A)") "Incoming:"
          call vector4_write (psum_in, u)
          write (u, "(A)") "Outgoing:"
          call vector4_write (psum_out, u)
       end if
    else
       if (verb) then
          write (u, "(A)") "Momentum conservation: CHECK"
       end if
    end if
  end subroutine vector4_check_momentum_conservation

  module subroutine spinor_product (p1, p2, prod1, prod2)
    type(vector4_t), intent(in) :: p1, p2
    complex(default), intent(out) :: prod1, prod2
    real(default) :: sij
    complex(default) :: phase
    real(default) :: pp_1, pp_2
    pp_1 = p1%p(0) + p1%p(3)
    pp_2 = p2%p(0) + p2%p(3)
    sij = (p1+p2)**2
    phase = cmplx ((p1%p(1)*pp_2 - p2%p(1)*pp_1)/sqrt (sij*pp_1*pp_2), &
                   (p1%p(2)*pp_2 - p2%p(2)*pp_1)/sqrt (sij*pp_1*pp_2), &
                    default)
    !!! <ij>
    prod1 = sqrt (sij) * phase
    !!! [ij]
    if (abs(prod1) > 0) then
       prod2 = - sij / prod1
    else
       prod2 = 0
    end if
  end subroutine spinor_product


end submodule lorentz_s

