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

submodule (sf_mappings) sf_mappings_s

  use io_units
  use constants, only: pi, zero, one
  use numeric_utils
  use diagnostics

  implicit none

contains

  module subroutine sf_mapping_base_init (mapping, n_par)
    class(sf_mapping_t), intent(out) :: mapping
    integer, intent(in) :: n_par
    allocate (mapping%i (n_par))
    mapping%i = 0
  end subroutine sf_mapping_base_init

  module subroutine sf_mapping_set_index (mapping, j, i)
    class(sf_mapping_t), intent(inout) :: mapping
    integer, intent(in) :: j, i
    mapping%i(j) = i
  end subroutine sf_mapping_set_index

  module function sf_mapping_get_index (mapping, j) result (i)
    class(sf_mapping_t), intent(inout) :: mapping
    integer, intent(in) :: j
    integer :: i
    i = mapping%i(j)
  end function sf_mapping_get_index

  module function sf_mapping_get_n_dim (mapping) result (n)
    class(sf_mapping_t), intent(in) :: mapping
    integer :: n
    n = size (mapping%i)
  end function sf_mapping_get_n_dim

  module subroutine sf_mapping_check (mapping, u, p_in, pb_in, fmt_p, fmt_f)
    class(sf_mapping_t), intent(inout) :: mapping
    integer, intent(in) :: u
    real(default), dimension(:), intent(in) :: p_in, pb_in
    character(*), intent(in) :: fmt_p
    character(*), intent(in), optional :: fmt_f
    real(default), dimension(size(p_in)) :: p, pb, r, rb
    real(default) :: f, tolerance
    tolerance = 1.5E-17_default
    p = p_in
    pb= pb_in
    call mapping%compute (r, rb, f, p, pb)
    call pacify (p, tolerance)
    call pacify (pb, tolerance)
    call pacify (r, tolerance)
    call pacify (rb, tolerance)
    write (u, "(3x,A,9(1x," // fmt_p // "))")  "p =", p
    write (u, "(3x,A,9(1x," // fmt_p // "))")  "pb=", pb
    write (u, "(3x,A,9(1x," // fmt_p // "))")  "r =", r
    write (u, "(3x,A,9(1x," // fmt_p // "))")  "rb=", rb
    if (present (fmt_f)) then
       write (u, "(3x,A,9(1x," // fmt_f // "))")  "f =", f
    else
       write (u, "(3x,A,9(1x," // fmt_p // "))")  "f =", f
    end if
    write (u, *)
    call mapping%inverse (r, rb, f, p, pb)
    call pacify (p, tolerance)
    call pacify (pb, tolerance)
    call pacify (r, tolerance)
    call pacify (rb, tolerance)
    write (u, "(3x,A,9(1x," // fmt_p // "))")  "p =", p
    write (u, "(3x,A,9(1x," // fmt_p // "))")  "pb=", pb
    write (u, "(3x,A,9(1x," // fmt_p // "))")  "r =", r
    write (u, "(3x,A,9(1x," // fmt_p // "))")  "rb=", rb
    if (present (fmt_f)) then
       write (u, "(3x,A,9(1x," // fmt_f // "))")  "f =", f
    else
       write (u, "(3x,A,9(1x," // fmt_p // "))")  "f =", f
    end if
    write (u, *)
    write (u, "(3x,A,9(1x," // fmt_p // "))")  "*r=", product (r)
  end subroutine sf_mapping_check

  module function sf_mapping_integral (mapping, n_calls) result (integral)
    class(sf_mapping_t), intent(inout) :: mapping
    integer, intent(in) :: n_calls
    real(default) :: integral
    integer :: n_dim, n_bin, k
    real(default), dimension(:), allocatable :: p, pb, r, rb
    integer, dimension(:), allocatable :: ii
    real(default) :: dx, f, s

    n_dim = mapping%get_n_dim ()
    allocate (p (n_dim))
    allocate (pb(n_dim))
    allocate (r (n_dim))
    allocate (rb(n_dim))
    allocate (ii(n_dim))
    n_bin = nint (real (n_calls, default) ** (1._default / n_dim))
    dx = 1._default / n_bin
    s = 0
    ii = 1

    SAMPLE: do
       do k = 1, n_dim
          p(k)  = ii(k) * dx - dx/2
          pb(k) = (n_bin - ii(k)) * dx + dx/2
       end do
       call mapping%compute (r, rb, f, p, pb)
       s = s + f
       INCR: do k = 1, n_dim
          ii(k) = ii(k) + 1
          if (ii(k) <= n_bin) then
             exit INCR
          else if (k < n_dim) then
             ii(k) = 1
          else
             exit SAMPLE
          end if
       end do INCR
    end do SAMPLE

    integral = s / real (n_bin, default) ** n_dim

  end function sf_mapping_integral

  module subroutine sf_s_mapping_write (object, unit)
    class(sf_s_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,',',I0,')')", advance="no")  object%i
    end if
    write (u, "(A,F7.5,A)")  ": standard (", object%power, ")"
  end subroutine sf_s_mapping_write

  module subroutine sf_s_mapping_init (mapping, power)
    class(sf_s_mapping_t), intent(out) :: mapping
    real(default), intent(in), optional :: power
    call mapping%base_init (2)
    if (present (power)) then
       mapping%power_set = .true.
       mapping%power = power
    end if
  end subroutine sf_s_mapping_init

  module subroutine sf_s_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_s_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: r2
    integer :: j
    if (mapping%power_set) then
       call map_unit_square (r2, f, p(mapping%i), mapping%power)
    else
       call map_unit_square (r2, f, p(mapping%i))
    end if
    r = p
    rb= pb
    do j = 1, 2
       r (mapping%i(j)) = r2(j)
       rb(mapping%i(j)) = 1 - r2(j)
    end do
  end subroutine sf_s_mapping_compute

  module subroutine sf_s_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_s_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: p2
    integer :: j
    if (mapping%power_set) then
       call map_unit_square_inverse (r(mapping%i), f, p2, mapping%power)
    else
       call map_unit_square_inverse (r(mapping%i), f, p2)
    end if
    p = r
    pb= rb
    do j = 1, 2
       p (mapping%i(j)) = p2(j)
       pb(mapping%i(j)) = 1 - p2(j)
    end do
  end subroutine sf_s_mapping_inverse

  module subroutine sf_res_mapping_write (object, unit)
    class(sf_res_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,',',I0,')')", advance="no")  object%i
    end if
    write (u, "(A,F7.5,', ',F7.5,A)")  ": resonance (", object%m, object%w, ")"
  end subroutine sf_res_mapping_write

  module subroutine sf_res_mapping_init (mapping, m, w)
    class(sf_res_mapping_t), intent(out) :: mapping
    real(default), intent(in) :: m, w
    call mapping%base_init (2)
    mapping%m = m
    mapping%w = w
  end subroutine sf_res_mapping_init

  module subroutine sf_res_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_res_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: r2, p2
    real(default) :: fbw, f2, p1m
    integer :: j
    p2 = p(mapping%i)
    call map_breit_wigner &
         (p1m, fbw, p2(1), mapping%m, mapping%w, x_free)
    call map_unit_square (r2, f2, [p1m, p2(2)])
    f = fbw * f2
    r = p
    rb= pb
    do j = 1, 2
       r (mapping%i(j)) = r2(j)
       rb(mapping%i(j)) = 1 - r2(j)
    end do
  end subroutine sf_res_mapping_compute

  module subroutine sf_res_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_res_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: p2
    real(default) :: fbw, f2, p1m
    call map_unit_square_inverse (r(mapping%i), f2, p2)
    call map_breit_wigner_inverse &
         (p2(1), fbw, p1m, mapping%m, mapping%w, x_free)
    p = r
    pb= rb
    p (mapping%i(1)) = p1m
    pb(mapping%i(1)) = 1 - p1m
    p (mapping%i(2)) = p2(2)
    pb(mapping%i(2)) = 1 - p2(2)
    f = fbw * f2
  end subroutine sf_res_mapping_inverse

  module subroutine sf_res_mapping_single_write (object, unit)
    class(sf_res_mapping_single_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,')')", advance="no")  object%i
    end if
    write (u, "(A,F7.5,', ',F7.5,A)")  ": resonance (", object%m, object%w, ")"
  end subroutine sf_res_mapping_single_write

  module subroutine sf_res_mapping_single_init (mapping, m, w)
    class(sf_res_mapping_single_t), intent(out) :: mapping
    real(default), intent(in) :: m, w
    call mapping%base_init (1)
    mapping%m = m
    mapping%w = w
  end subroutine sf_res_mapping_single_init

  module subroutine sf_res_mapping_single_compute &
       (mapping, r, rb, f, p, pb, x_free)
    class(sf_res_mapping_single_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(1) :: r2, p2
    real(default) :: fbw
    integer :: j
    p2 = p(mapping%i)
    call map_breit_wigner &
         (r2(1), fbw, p2(1), mapping%m, mapping%w, x_free)
    f = fbw
    r = p
    rb= pb
    r (mapping%i(1)) = r2(1)
    rb(mapping%i(1)) = 1 - r2(1)
  end subroutine sf_res_mapping_single_compute

  module subroutine sf_res_mapping_single_inverse &
       (mapping, r, rb, f, p, pb, x_free)
    class(sf_res_mapping_single_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(1) :: p2
    real(default) :: fbw
    call map_breit_wigner_inverse &
         (r(mapping%i(1)), fbw, p2(1), mapping%m, mapping%w, x_free)
    p = r
    pb= rb
    p (mapping%i(1)) = p2(1)
    pb(mapping%i(1)) = 1 - p2(1)
    f = fbw
  end subroutine sf_res_mapping_single_inverse

  module subroutine sf_os_mapping_write (object, unit)
    class(sf_os_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,',',I0,')')", advance="no")  object%i
    end if
    write (u, "(A,F7.5,A)")  ": on-shell (", object%m, ")"
  end subroutine sf_os_mapping_write

  module subroutine sf_os_mapping_init (mapping, m)
    class(sf_os_mapping_t), intent(out) :: mapping
    real(default), intent(in) :: m
    call mapping%base_init (2)
    mapping%m = m
    mapping%lm2 = abs (2 * log (mapping%m))
  end subroutine sf_os_mapping_init

  module subroutine sf_os_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_os_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: r2, p2
    integer :: j
    p2 = p(mapping%i)
    call map_on_shell (r2, f, p2, mapping%lm2, x_free)
    r = p
    rb= pb
    do j = 1, 2
       r (mapping%i(j)) = r2(j)
       rb(mapping%i(j)) = 1 - r2(j)
    end do
  end subroutine sf_os_mapping_compute

  module subroutine sf_os_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_os_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: p2, r2
    r2 = r(mapping%i)
    call map_on_shell_inverse (r2, f, p2, mapping%lm2, x_free)
    p = r
    pb= rb
    p (mapping%i(1)) = p2(1)
    pb(mapping%i(1)) = 1 - p2(1)
    p (mapping%i(2)) = p2(2)
    pb(mapping%i(2)) = 1 - p2(2)
  end subroutine sf_os_mapping_inverse

  module subroutine sf_os_mapping_single_write (object, unit)
    class(sf_os_mapping_single_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,')')", advance="no")  object%i
    end if
    write (u, "(A,F7.5,A)")  ": on-shell (", object%m, ")"
  end subroutine sf_os_mapping_single_write

  module subroutine sf_os_mapping_single_init (mapping, m)
    class(sf_os_mapping_single_t), intent(out) :: mapping
    real(default), intent(in) :: m
    call mapping%base_init (1)
    mapping%m = m
    mapping%lm2 = abs (2 * log (mapping%m))
  end subroutine sf_os_mapping_single_init

  module subroutine sf_os_mapping_single_compute &
       (mapping, r, rb, f, p, pb, x_free)
    class(sf_os_mapping_single_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(1) :: r2, p2
    integer :: j
    p2 = p(mapping%i)
    call map_on_shell_single (r2, f, p2, mapping%lm2, x_free)
    r = p
    rb= pb
    r (mapping%i(1)) = r2(1)
    rb(mapping%i(1)) = 1 - r2(1)
  end subroutine sf_os_mapping_single_compute

  module subroutine sf_os_mapping_single_inverse &
       (mapping, r, rb, f, p, pb, x_free)
    class(sf_os_mapping_single_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(1) :: p2, r2
    r2 = r(mapping%i)
    call map_on_shell_single_inverse (r2, f, p2, mapping%lm2, x_free)
    p = r
    pb= rb
    p (mapping%i(1)) = p2(1)
    pb(mapping%i(1)) = 1 - p2(1)
  end subroutine sf_os_mapping_single_inverse

  module subroutine sf_ep_mapping_write (object, unit)
    class(sf_ep_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,',',I0,')')", advance="no")  object%i
    end if
    write (u, "(A,ES12.5,A)")  ": endpoint (a =", object%a, ")"
  end subroutine sf_ep_mapping_write

  module subroutine sf_ep_mapping_init (mapping, a)
    class(sf_ep_mapping_t), intent(out) :: mapping
    real(default), intent(in), optional :: a
    call mapping%base_init (2)
    if (present (a))  mapping%a = a
  end subroutine sf_ep_mapping_init

  module subroutine sf_ep_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_ep_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: px, r2
    real(default) :: f1, f2
    integer :: j
    call map_endpoint_1 (px(1), f1, p(mapping%i(1)), mapping%a)
    call map_endpoint_01 (px(2), f2, p(mapping%i(2)), mapping%a)
    call map_unit_square (r2, f, px)
    f = f * f1 * f2
    r = p
    rb= pb
    do j = 1, 2
       r (mapping%i(j)) = r2(j)
       rb(mapping%i(j)) = 1 - r2(j)
    end do
  end subroutine sf_ep_mapping_compute

  module subroutine sf_ep_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_ep_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: r2, px, p2
    real(default) :: f1, f2
    integer :: j
    do j = 1, 2
       r2(j) = r(mapping%i(j))
    end do
    call map_unit_square_inverse (r2, f, px)
    call map_endpoint_inverse_1 (px(1), f1, p2(1), mapping%a)
    call map_endpoint_inverse_01 (px(2), f2, p2(2), mapping%a)
    f = f * f1 * f2
    p = r
    pb= rb
    do j = 1, 2
       p (mapping%i(j)) = p2(j)
       pb(mapping%i(j)) = 1 - p2(j)
    end do
  end subroutine sf_ep_mapping_inverse

  module subroutine sf_epr_mapping_write (object, unit)
    class(sf_epr_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,',',I0,')')", advance="no")  object%i
    end if
    if (object%resonance) then
       write (u, "(A,F7.5,A,F7.5,', ',F7.5,A)")  ": ep/res (a = ", object%a, &
            " | ", object%m, object%w, ")"
    else
       write (u, "(A,F7.5,A)")  ": ep/nores (a = ", object%a, ")"
    end if
  end subroutine sf_epr_mapping_write

  module subroutine sf_epr_mapping_init (mapping, a, m, w)
    class(sf_epr_mapping_t), intent(out) :: mapping
    real(default), intent(in) :: a
    real(default), intent(in), optional :: m, w
    call mapping%base_init (2)
    mapping%a = a
    if (present (m) .and. present (w)) then
       mapping%m = m
       mapping%w = w
    else
       mapping%resonance = .false.
    end if
  end subroutine sf_epr_mapping_init

  module subroutine sf_epr_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_epr_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: px, r2
    real(default) :: f1, f2
    integer :: j
    if (mapping%resonance) then
       call map_breit_wigner &
            (px(1), f1, p(mapping%i(1)), mapping%m, mapping%w, x_free)
    else
       px(1) = p(mapping%i(1))
       f1 = 1
    end if
    call map_endpoint_01 (px(2), f2, p(mapping%i(2)), mapping%a)
    call map_unit_square (r2, f, px)
    f = f * f1 * f2
    r = p
    rb= pb
    do j = 1, 2
       r (mapping%i(j)) = r2(j)
       rb(mapping%i(j)) = 1 - r2(j)
    end do
  end subroutine sf_epr_mapping_compute

  module subroutine sf_epr_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_epr_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: px, p2
    real(default) :: f1, f2
    integer :: j
    call map_unit_square_inverse (r(mapping%i), f, px)
    if (mapping%resonance) then
       call map_breit_wigner_inverse &
            (px(1), f1, p2(1), mapping%m, mapping%w, x_free)
    else
       p2(1) = px(1)
       f1 = 1
    end if
    call map_endpoint_inverse_01 (px(2), f2, p2(2), mapping%a)
    f = f * f1 * f2
    p = r
    pb= rb
    do j = 1, 2
       p (mapping%i(j)) = p2(j)
       pb(mapping%i(j)) = 1 - p2(j)
    end do
  end subroutine sf_epr_mapping_inverse

  module subroutine sf_epo_mapping_write (object, unit)
    class(sf_epo_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,',',I0,')')", advance="no")  object%i
    end if
    write (u, "(A,F7.5,A,F7.5,A)")  ": ep/on-shell (a = ", object%a, &
         " | ", object%m, ")"
  end subroutine sf_epo_mapping_write

  module subroutine sf_epo_mapping_init (mapping, a, m)
    class(sf_epo_mapping_t), intent(out) :: mapping
    real(default), intent(in) :: a, m
    call mapping%base_init (2)
    mapping%a = a
    mapping%m = m
    mapping%lm2 = abs (2 * log (mapping%m))
  end subroutine sf_epo_mapping_init

  module subroutine sf_epo_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_epo_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: px, r2
    real(default) :: f2
    integer :: j
    px(1) = 0
    call map_endpoint_01 (px(2), f2, p(mapping%i(2)), mapping%a)
    call map_on_shell (r2, f, px, mapping%lm2)
    f = f * f2
    r = p
    rb= pb
    do j = 1, 2
       r (mapping%i(j)) = r2(j)
       rb(mapping%i(j)) = 1 - r2(j)
    end do
  end subroutine sf_epo_mapping_compute

  module subroutine sf_epo_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_epo_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: px, p2
    real(default) :: f2
    integer :: j
    call map_on_shell_inverse (r(mapping%i), f, px, mapping%lm2)
    p2(1) = 0
    call map_endpoint_inverse_01 (px(2), f2, p2(2), mapping%a)
    f = f * f2
    p = r
    pb= rb
    do j = 1, 2
       p (mapping%i(j)) = p2(j)
       pb(mapping%i(j)) = 1 - p2(j)
    end do
  end subroutine sf_epo_mapping_inverse

  module subroutine sf_ip_mapping_write (object, unit)
    class(sf_ip_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,',',I0,')')", advance="no")  object%i
    end if
    write (u, "(A,ES12.5,A)")  ": isr (eps =", object%eps, ")"
  end subroutine sf_ip_mapping_write

  module subroutine sf_ip_mapping_init (mapping, eps)
    class(sf_ip_mapping_t), intent(out) :: mapping
    real(default), intent(in), optional :: eps
    call mapping%base_init (2)
    if (present (eps))  mapping%eps = eps
    if (mapping%eps <= 0) &
         call msg_fatal ("ISR mapping: regulator epsilon must not be zero")
  end subroutine sf_ip_mapping_init

  module subroutine sf_ip_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_ip_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: px, pxb, r2, r2b
    real(default) :: f1, f2, xb, y, yb
    integer :: j
    call map_power_1 (xb, f1, pb(mapping%i(1)), 2 * mapping%eps)
    call map_power_01 (y, yb, f2, pb(mapping%i(2)), mapping%eps)
    px(1)  = 1 - xb
    pxb(1) = xb
    px(2)  = y
    pxb(2) = yb
    call map_unit_square_prec (r2, r2b, f, px, pxb)
    f = f * f1 * f2
    r = p
    rb= pb
    do j = 1, 2
       r (mapping%i(j)) = r2 (j)
       rb(mapping%i(j)) = r2b(j)
    end do
  end subroutine sf_ip_mapping_compute

  module subroutine sf_ip_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_ip_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: r2, r2b, px, pxb, p2, p2b
    real(default) :: f1, f2, xb, y, yb
    integer :: j
    do j = 1, 2
       r2 (j) = r (mapping%i(j))
       r2b(j) = rb(mapping%i(j))
    end do
    call map_unit_square_inverse_prec (r2, r2b, f, px, pxb)
    xb = pxb(1)
    if (px(1) > 0) then
       y  = px(2)
       yb = pxb(2)
    else
       y  = 0.5_default
       yb = 0.5_default
    end if
    call map_power_inverse_1 (xb, f1, p2b(1), 2 * mapping%eps)
    call map_power_inverse_01 (y, yb, f2, p2b(2), mapping%eps)
    p2 = 1 - p2b
    f = f * f1 * f2
    p = r
    pb= rb
    do j = 1, 2
       p (mapping%i(j)) = p2(j)
       pb(mapping%i(j)) = p2b(j)
    end do
  end subroutine sf_ip_mapping_inverse

  module subroutine sf_ipr_mapping_write (object, unit)
    class(sf_ipr_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,',',I0,')')", advance="no")  object%i
    end if
    if (object%resonance) then
       write (u, "(A,F7.5,A,F7.5,', ',F7.5,A)")  ": isr/res (eps = ", &
            object%eps, " | ", object%m, object%w, ")"
    else
       write (u, "(A,F7.5,A)")  ": isr/res (eps = ", object%eps, ")"
    end if
  end subroutine sf_ipr_mapping_write

  module subroutine sf_ipr_mapping_init (mapping, eps, m, w)
    class(sf_ipr_mapping_t), intent(out) :: mapping
    real(default), intent(in), optional :: eps, m, w
    call mapping%base_init (2)
    if (present (eps))  mapping%eps = eps
    if (mapping%eps <= 0) &
         call msg_fatal ("ISR mapping: regulator epsilon must not be zero")
    if (present (m) .and. present (w)) then
       mapping%m = m
       mapping%w = w
    else
       mapping%resonance = .false.
    end if
  end subroutine sf_ipr_mapping_init

  module subroutine sf_ipr_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_ipr_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: px, pxb, r2, r2b
    real(default) :: f1, f2, y, yb
    integer :: j
    if (mapping%resonance) then
       call map_breit_wigner &
            (px(1), f1, p(mapping%i(1)), mapping%m, mapping%w, x_free)
    else
       px(1) = p(mapping%i(1))
       f1 = 1
    end if
    call map_power_01 (y, yb, f2, pb(mapping%i(2)), mapping%eps)
    pxb(1) = 1 - px(1)
    px(2)  = y
    pxb(2) = yb
    call map_unit_square_prec (r2, r2b, f, px, pxb)
    f = f * f1 * f2
    r = p
    rb= pb
    do j = 1, 2
       r (mapping%i(j)) = r2 (j)
       rb(mapping%i(j)) = r2b(j)
    end do
  end subroutine sf_ipr_mapping_compute

  module subroutine sf_ipr_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_ipr_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: r2, r2b, px, pxb, p2, p2b
    real(default) :: f1, f2, y, yb
    integer :: j
    do j = 1, 2
       r2 (j) = r (mapping%i(j))
       r2b(j) = rb(mapping%i(j))
    end do
    call map_unit_square_inverse_prec (r2, r2b, f, px, pxb)
    if (px(1) > 0) then
       y  = px(2)
       yb = pxb(2)
    else
       y  = 0.5_default
       yb = 0.5_default
    end if
    if (mapping%resonance) then
       call map_breit_wigner_inverse &
            (px(1), f1, p2(1), mapping%m, mapping%w, x_free)
    else
       p2(1) = px(1)
       f1 = 1
    end if
    call map_power_inverse_01 (y, yb, f2, p2b(2), mapping%eps)
    p2b(1) = 1 - p2(1)
    p2 (2) = 1 - p2b(2)
    f = f * f1 * f2
    p = r
    pb= rb
    do j = 1, 2
       p (mapping%i(j)) = p2(j)
       pb(mapping%i(j)) = p2b(j)
    end do
  end subroutine sf_ipr_mapping_inverse

  module subroutine sf_ipo_mapping_write (object, unit)
    class(sf_ipo_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,',',I0,')')", advance="no")  object%i
    end if
    write (u, "(A,F7.5,A,F7.5,A)")  ": isr/os (eps = ", object%eps, &
         " | ", object%m, ")"
  end subroutine sf_ipo_mapping_write

  module subroutine sf_ipo_mapping_init (mapping, eps, m)
    class(sf_ipo_mapping_t), intent(out) :: mapping
    real(default), intent(in), optional :: eps, m
    call mapping%base_init (2)
    if (present (eps))  mapping%eps = eps
    if (mapping%eps <= 0) &
         call msg_fatal ("ISR mapping: regulator epsilon must not be zero")
    mapping%m = m
  end subroutine sf_ipo_mapping_init

  module subroutine sf_ipo_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_ipo_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: px, pxb, r2, r2b
    real(default) :: f1, f2, y, yb
    integer :: j
    call map_power_01 (y, yb, f2, pb(mapping%i(2)), mapping%eps)
    px(1)  = mapping%m ** 2
    if (present (x_free))  px(1) = px(1) / x_free
    pxb(1) = 1 - px(1)
    px(2)  = y
    pxb(2) = yb
    call map_unit_square_prec (r2, r2b, f1, px, pxb)
    f = f1 * f2
    r = p
    rb= pb
    do j = 1, 2
       r (mapping%i(j)) = r2 (j)
       rb(mapping%i(j)) = r2b(j)
    end do
  end subroutine sf_ipo_mapping_compute

  module subroutine sf_ipo_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_ipo_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(2) :: r2, r2b, px, pxb, p2, p2b
    real(default) :: f1, f2, y, yb
    integer :: j
    do j = 1, 2
       r2 (j) = r (mapping%i(j))
       r2b(j) = rb(mapping%i(j))
    end do
    call map_unit_square_inverse_prec (r2, r2b, f1, px, pxb)
    y  = px(2)
    yb = pxb(2)
    call map_power_inverse_01 (y, yb, f2, p2b(2), mapping%eps)
    p2(1) = 0
    p2b(1)= 1
    p2(2) = 1 - p2b(2)
    f = f1 * f2
    p = r
    pb= rb
    do j = 1, 2
       p (mapping%i(j)) = p2(j)
       pb(mapping%i(j)) = p2b(j)
    end do
  end subroutine sf_ipo_mapping_inverse

  module subroutine sf_ei_mapping_write (object, unit)
    class(sf_ei_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,3(',',I0),')')", advance="no")  object%i
    end if
    write (u, "(A,ES12.5,A,ES12.5,A)")  ": ep/isr (a =", object%ep%a, &
         ", eps =", object%ip%eps, ")"
  end subroutine sf_ei_mapping_write

  module subroutine sf_ei_mapping_init (mapping, a, eps)
    class(sf_ei_mapping_t), intent(out) :: mapping
    real(default), intent(in), optional :: a, eps
    call mapping%base_init (4)
    call mapping%ep%init (a)
    call mapping%ip%init (eps)
  end subroutine sf_ei_mapping_init

  module subroutine sf_ei_mapping_set_index (mapping, j, i)
    class(sf_ei_mapping_t), intent(inout) :: mapping
    integer, intent(in) :: j, i
    mapping%i(j) = i
    select case (j)
    case (1:2);  call mapping%ep%set_index (j, i)
    case (3:4);  call mapping%ip%set_index (j-2, i)
    end select
  end subroutine sf_ei_mapping_set_index

  module subroutine sf_ei_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_ei_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(size(p)) :: q, qb
    real(default) :: f1, f2
    call mapping%ep%compute (q, qb, f1, p, pb, x_free)
    call mapping%ip%compute (r, rb, f2, q, qb, x_free)
    f = f1 * f2
  end subroutine sf_ei_mapping_compute

  module subroutine sf_ei_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_ei_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(size(p)) :: q, qb
    real(default) :: f1, f2
    call mapping%ip%inverse (r, rb, f2, q, qb, x_free)
    call mapping%ep%inverse (q, qb, f1, p, pb, x_free)
    f = f1 * f2
  end subroutine sf_ei_mapping_inverse

  module subroutine sf_eir_mapping_write (object, unit)
    class(sf_eir_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,3(',',I0),')')", advance="no")  object%i
    end if
    write (u, "(A,F7.5,A,F7.5,A,F7.5,', ',F7.5,A)")  &
         ": ep/isr/res (a =", object%ep%a, &
         ", eps =", object%ip%eps, " | ", object%res%m, object%res%w, ")"
  end subroutine sf_eir_mapping_write

  module subroutine sf_eir_mapping_init (mapping, a, eps, m, w)
    class(sf_eir_mapping_t), intent(out) :: mapping
    real(default), intent(in) :: a, eps, m, w
    call mapping%base_init (4)
    call mapping%res%init (m, w)
    call mapping%ep%init (a)
    call mapping%ip%init (eps)
  end subroutine sf_eir_mapping_init

  module subroutine sf_eir_mapping_set_index (mapping, j, i)
    class(sf_eir_mapping_t), intent(inout) :: mapping
    integer, intent(in) :: j, i
    mapping%i(j) = i
    select case (j)
    case (1);  call mapping%res%set_index (1, i)
    case (3);  call mapping%res%set_index (2, i)
    end select
    select case (j)
    case (1:2);  call mapping%ep%set_index (j, i)
    case (3:4);  call mapping%ip%set_index (j-2, i)
    end select
  end subroutine sf_eir_mapping_set_index

  module subroutine sf_eir_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_eir_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(size(p)) :: px, pxb, q, qb
    real(default) :: f0, f1, f2
    call mapping%res%compute (px, pxb, f0, p, pb, x_free)
    call mapping%ep%compute (q, qb, f1, px, pxb, x_free)
    call mapping%ip%compute (r, rb, f2, q, qb, x_free)
    f = f0 * f1 * f2
  end subroutine sf_eir_mapping_compute

  module subroutine sf_eir_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_eir_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(size(p)) :: px, pxb, q, qb
    real(default) :: f0, f1, f2
    call mapping%ip%inverse (r, rb, f2, q, qb, x_free)
    call mapping%ep%inverse (q, qb, f1, px, pxb, x_free)
    call mapping%res%inverse (px, pxb, f0, p, pb, x_free)
    f = f0 * f1 * f2
  end subroutine sf_eir_mapping_inverse

  module subroutine sf_eio_mapping_write (object, unit)
    class(sf_eio_mapping_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "map"
    if (any (object%i /= 0)) then
       write (u, "('(',I0,3(',',I0),')')", advance="no")  object%i
    end if
    write (u, "(A,F7.5,A,F7.5,A,F7.5,A)")  ": ep/isr/os (a =", object%ep%a, &
         ", eps =", object%ip%eps, " | ", object%os%m, ")"
  end subroutine sf_eio_mapping_write

  module subroutine sf_eio_mapping_init (mapping, a, eps, m)
    class(sf_eio_mapping_t), intent(out) :: mapping
    real(default), intent(in), optional :: a, eps, m
    call mapping%base_init (4)
    call mapping%os%init (m)
    call mapping%ep%init (a)
    call mapping%ip%init (eps)
  end subroutine sf_eio_mapping_init

  module subroutine sf_eio_mapping_set_index (mapping, j, i)
    class(sf_eio_mapping_t), intent(inout) :: mapping
    integer, intent(in) :: j, i
    mapping%i(j) = i
    select case (j)
    case (1);  call mapping%os%set_index (1, i)
    case (3);  call mapping%os%set_index (2, i)
    end select
    select case (j)
    case (1:2);  call mapping%ep%set_index (j, i)
    case (3:4);  call mapping%ip%set_index (j-2, i)
    end select
  end subroutine sf_eio_mapping_set_index

  module subroutine sf_eio_mapping_compute (mapping, r, rb, f, p, pb, x_free)
    class(sf_eio_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(size(p)) :: px, pxb, q, qb
    real(default) :: f0, f1, f2
    call mapping%os%compute (px, pxb, f0, p, pb, x_free)
    call mapping%ep%compute (q, qb, f1, px, pxb, x_free)
    call mapping%ip%compute (r, rb, f2, q, qb, x_free)
    f = f0 * f1 * f2
  end subroutine sf_eio_mapping_compute

  module subroutine sf_eio_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
    class(sf_eio_mapping_t), intent(inout) :: mapping
    real(default), dimension(:), intent(in) :: r, rb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: p, pb
    real(default), intent(inout), optional :: x_free
    real(default), dimension(size(p)) :: px, pxb, q, qb
    real(default) :: f0, f1, f2
    call mapping%ip%inverse (r, rb, f2, q, qb, x_free)
    call mapping%ep%inverse (q, qb, f1, px, pxb, x_free)
    call mapping%os%inverse (px, pxb, f0, p, pb, x_free)
    f = f0 * f1 * f2
  end subroutine sf_eio_mapping_inverse

  subroutine map_unit_square (r, factor, p, power)
    real(default), dimension(2), intent(out) :: r
    real(default), intent(out) :: factor
    real(default), dimension(2), intent(in) :: p
    real(default), intent(in), optional :: power
    real(default) :: xx, yy
    factor = 1
    xx = p(1)
    yy = p(2)
    if (present(power)) then
       if (p(1) > 0 .and. power > 1) then
          xx = p(1)**power
          factor = factor * power * xx / p(1)
       end if
    end if
    if (.not. vanishes (xx)) then
       r(1) = xx ** yy
       r(2) = xx / r(1)
       factor = factor * abs (log (xx))
    else
       r = 0
    end if
  end subroutine map_unit_square

  subroutine map_unit_square_inverse (r, factor, p, power)
    real(kind=default), dimension(2), intent(in) :: r
    real(kind=default), intent(out) :: factor
    real(kind=default), dimension(2), intent(out) :: p
    real(kind=default), intent(in), optional :: power
    real(kind=default) :: lg, xx, yy
    factor = 1
    xx = r(1) * r(2)
    if (.not. vanishes (xx)) then
       lg = log (xx)
       if (.not. vanishes (lg)) then
          yy = log (r(1)) / lg
       else
          yy = 0
       end if
       p(2) = yy
       factor = factor * abs (lg)
       if (present(power)) then
          p(1) = xx**(1._default/power)
          factor = factor * power * xx / p(1)
       else
          p(1) = xx
       end if
    else
       p = 0
    end if
  end subroutine map_unit_square_inverse

  subroutine map_unit_square_prec (r, rb, factor, p, pb)
    real(default), dimension(2), intent(out) :: r
    real(default), dimension(2), intent(out) :: rb
    real(default), intent(out) :: factor
    real(default), dimension(2), intent(in) :: p
    real(default), dimension(2), intent(in) :: pb
    if (p(1) > 0.5_default) then
       call compute_prec_xy_1 (r(1), rb(1), p(1), pb(1), p (2))
       call compute_prec_xy_1 (r(2), rb(2), p(1), pb(1), pb(2))
       factor = - log_prec (p(1), pb(1))
    else if (.not. vanishes (p(1))) then
       call compute_prec_xy_0 (r(1), rb(1), p(1), pb(1), p (2))
       call compute_prec_xy_0 (r(2), rb(2), p(1), pb(1), pb(2))
       factor = - log_prec (p(1), pb(1))
    else
       r  = 0
       rb = 1
       factor = 0
    end if
  end subroutine map_unit_square_prec

  subroutine map_unit_square_inverse_prec (r, rb, factor, p, pb)
    real(default), dimension(2), intent(in) :: r
    real(default), dimension(2), intent(in) :: rb
    real(default), intent(out) :: factor
    real(default), dimension(2), intent(out) :: p
    real(default), dimension(2), intent(out) :: pb
    call inverse_prec_x (r, rb, p(1), pb(1))
    if (all (r > 0)) then
       if (rb(1) < rb(2)) then
          call inverse_prec_y (r, rb, p(2), pb(2))
       else
          call inverse_prec_y ([r(2),r(1)], [rb(2),rb(1)], pb(2), p(2))
       end if
       factor = - log_prec (p(1), pb(1))
    else
       p(1)  = 0
       pb(1) = 1
       p(2)  = 0.5_default
       pb(2) = 0.5_default
       factor = 0
    end if
  end subroutine map_unit_square_inverse_prec

  subroutine compute_prec_xy_1 (z, zb, x, xb, y)
    real(default), intent(out) :: z, zb
    real(default), intent(in) :: x, xb, y
    real(default) :: a1, a2, a3
    a1 = y * xb
    a2 = a1 * (1 - y) * xb / 2
    a3 = a2 * (2 - y) * xb / 3
    if (abs (a3) < epsilon (a3)) then
       zb = a1 + a2 + a3
       z = 1 - zb
    else
       z = x ** y
       zb = 1 - z
    end if
  end subroutine compute_prec_xy_1

  subroutine compute_prec_xy_0 (z, zb, x, xb, y)
    real(default), intent(out) :: z, zb
    real(default), intent(in) :: x, xb, y
    real(default) :: a1, a2, a3, lx
    lx = -log (x)
    a1 = y * lx
    a2 = a1 * y * lx / 2
    a3 = a2 * y * lx / 3
    if (abs (a3) < epsilon (a3)) then
       zb = a1 + a2 + a3
       z = 1 - zb
    else
       z = x ** y
       zb = 1 - z
    end if
  end subroutine compute_prec_xy_0

  subroutine inverse_prec_x (r, rb, x, xb)
    real(default), dimension(2), intent(in) :: r, rb
    real(default), intent(out) :: x, xb
    real(default) :: a0, a1
    a0 = rb(1) + rb(2)
    a1 = rb(1) * rb(2)
    if (a0 > 0.5_default) then
       xb = a0 - a1
       x = 1 - xb
    else
       x = r(1) * r(2)
       xb = 1 - x
    end if
  end subroutine inverse_prec_x

  subroutine inverse_prec_y (r, rb, y, yb)
    real(default), dimension(2), intent(in) :: r, rb
    real(default), intent(out) :: y, yb
    real(default) :: log1, log2, a1, a2, a3
    log1 = log_prec (r(1), rb(1))
    log2 = log_prec (r(2), rb(2))
    if (abs (log2**3) < epsilon (one)) then
       if (abs(log1) < epsilon (one)) then
          y = zero
       else
          y = one / (one + log2 / log1)
       end if
       if (abs(log2) < epsilon (one)) then
          yb = zero
       else
          yb = one / (one + log1 / log2)
       end if
       return
    end if
    a1 = - rb(1) / log2
    a2 = - rb(1) ** 2 * (one / log2**2 + one / (2 * log2))
    a3 = - rb(1) ** 3 * (one / log2**3 + one / log2**2 + one/(3 * log2))
    if (abs (a3) < epsilon (a3)) then
       y  = a1 + a2 + a3
       yb = one - y
    else
       y  = one / (one + log2 / log1)
       yb = one / (one + log1 / log2)
    end if
  end subroutine inverse_prec_y

  module subroutine map_on_shell (r, factor, p, lm2, x_free)
    real(default), dimension(2), intent(out) :: r
    real(default), intent(out) :: factor
    real(default), dimension(2), intent(in) :: p
    real(default), intent(in) :: lm2
    real(default), intent(in), optional :: x_free
    real(default) :: lx
    lx = lm2;  if (present (x_free))  lx = lx + log (x_free)
    r(1) = exp (- p(2) * lx)
    r(2) = exp (- (1 - p(2)) * lx)
    factor = lx
  end subroutine map_on_shell

  module subroutine map_on_shell_inverse (r, factor, p, lm2, x_free)
    real(default), dimension(2), intent(in) :: r
    real(default), intent(out) :: factor
    real(default), dimension(2), intent(out) :: p
    real(default), intent(in) :: lm2
    real(default), intent(in), optional :: x_free
    real(default) :: lx
    lx = lm2;  if (present (x_free))  lx = lx + log (x_free)
    p(1) = 0
    p(2) = abs (log (r(1))) / lx
    factor = lx
  end subroutine map_on_shell_inverse

  module subroutine map_on_shell_single (r, factor, p, lm2, x_free)
    real(default), dimension(1), intent(out) :: r
    real(default), intent(out) :: factor
    real(default), dimension(1), intent(in) :: p
    real(default), intent(in) :: lm2
    real(default), intent(in), optional :: x_free
    real(default) :: lx
    lx = lm2;  if (present (x_free))  lx = lx + log (x_free)
    r(1) = exp (- lx)
    factor = 1
  end subroutine map_on_shell_single

  module subroutine map_on_shell_single_inverse (r, factor, p, lm2, x_free)
    real(default), dimension(1), intent(in) :: r
    real(default), intent(out) :: factor
    real(default), dimension(1), intent(out) :: p
    real(default), intent(in) :: lm2
    real(default), intent(in), optional :: x_free
    real(default) :: lx
    lx = lm2;  if (present (x_free))  lx = lx + log (x_free)
    p(1) = 0
    factor = 1
  end subroutine map_on_shell_single_inverse

  subroutine map_breit_wigner (r, factor, p, m, w, x_free)
    real(default), intent(out) :: r
    real(default), intent(out) :: factor
    real(default), intent(in) :: p
    real(default), intent(in) :: m
    real(default), intent(in) :: w
    real(default), intent(in), optional :: x_free
    real(default) :: m2, mw, a1, a2, a3, z, tmp
    m2 = m ** 2
    mw = m * w
    if (present (x_free)) then
       m2 = m2 / x_free
       mw = mw / x_free
    end if
    a1 = atan (- m2 / mw)
    a2 = atan ((1 - m2) / mw)
    a3 = (a2 - a1) * mw
    z = (1-p) * a1 + p * a2
    if (-pi/2 < z .and. z < pi/2) then
       tmp = tan (z)
       r = max (m2 + mw * tmp, 0._default)
       factor = a3 * (1 + tmp ** 2)
    else
       r = 0
       factor = 0
    end if
  end subroutine map_breit_wigner

  subroutine map_breit_wigner_inverse (r, factor, p, m, w, x_free)
    real(default), intent(in) :: r
    real(default), intent(out) :: factor
    real(default), intent(out) :: p
    real(default), intent(in) :: m
    real(default), intent(in) :: w
    real(default) :: m2, mw, a1, a2, a3, tmp
    real(default), intent(in), optional :: x_free
    m2 = m ** 2
    mw = m * w
    if (present (x_free)) then
       m2 = m2 / x_free
       mw = mw / x_free
    end if
    a1 = atan (- m2 / mw)
    a2 = atan ((1 - m2) / mw)
    a3 = (a2 - a1) * mw
    tmp = (r - m2) / mw
    p = (atan (tmp) - a1) / (a2 - a1)
    factor = a3 * (1 + tmp ** 2)
  end subroutine map_breit_wigner_inverse

  subroutine map_endpoint_1 (x3, factor, x1, a)
    real(default), intent(out) :: x3, factor
    real(default), intent(in) :: x1
    real(default), intent(in) :: a
    real(default) :: x2
    if (abs (x1) < 1) then
       x2 = tan (x1 * pi / 2)
       x3 = tanh (a * x2)
       factor = a * pi/2 * (1 + x2 ** 2) * (1 - x3 ** 2)
    else
       x3 = x1
       factor = 0
    end if
  end subroutine map_endpoint_1

  subroutine map_endpoint_inverse_1 (x3, factor, x1, a)
    real(default), intent(in) :: x3
    real(default), intent(out) :: x1, factor
    real(default), intent(in) :: a
    real(default) :: x2
    if (abs (x3) < 1) then
       x2 = atanh (x3) / a
       x1 = 2 / pi * atan (x2)
       factor = a * pi/2 * (1 + x2 ** 2) * (1 - x3 ** 2)
    else
       x1 = x3
       factor = 0
    end if
  end subroutine map_endpoint_inverse_1

  subroutine map_endpoint_01 (x4, factor, x0, a)
    real(default), intent(out) :: x4, factor
    real(default), intent(in) :: x0
    real(default), intent(in) :: a
    real(default) :: x1, x3
    x1 = 2 * x0 - 1
    call map_endpoint_1 (x3, factor, x1, a)
    x4 = (x3 + 1) / 2
  end subroutine map_endpoint_01

  subroutine map_endpoint_inverse_01 (x4, factor, x0, a)
    real(default), intent(in) :: x4
    real(default), intent(out) :: x0, factor
    real(default), intent(in) :: a
    real(default) :: x1, x3
    x3 = 2 * x4 - 1
    call map_endpoint_inverse_1 (x3, factor, x1, a)
    x0 = (x1 + 1) / 2
  end subroutine map_endpoint_inverse_01

  module subroutine map_power_1 (xb, factor, rb, eps)
    real(default), intent(out) :: xb, factor
    real(default), intent(in) :: rb
    real(double) :: rb_db, factor_db, eps_db, xb_db
    real(default), intent(in) :: eps
    rb_db = real (rb, kind=double)
    eps_db = real (eps, kind=double)
    xb_db = rb_db ** (1 / eps_db)
    if (rb_db > 0) then
       factor_db = xb_db / rb_db / eps_db
       factor = real (factor_db, kind=default)
    else
       factor = 0
    end if
    xb = real (xb_db, kind=default)
  end subroutine map_power_1

  module subroutine map_power_inverse_1 (xb, factor, rb, eps)
    real(default), intent(in) :: xb
    real(default), intent(out) :: rb, factor
    real(double) :: xb_db, factor_db, eps_db, rb_db
    real(default), intent(in) :: eps
    xb_db = real (xb, kind=double)
    eps_db = real (eps, kind=double)
    rb_db = xb_db ** eps_db
    if (xb_db > 0) then
       factor_db = xb_db / rb_db / eps_db
       factor = real (factor_db, kind=default)
    else
       factor = 0
    end if
    rb = real (rb_db, kind=default)
  end subroutine map_power_inverse_1

  subroutine map_power_01 (y, yb, factor, r, eps)
    real(default), intent(out) :: y, yb, factor
    real(default), intent(in) :: r
    real(default), intent(in) :: eps
    real(default) :: u, ub, zp, zm
    u = 2 * r - 1
    if (u > 0) then
       ub = 2 * (1 - r)
       call map_power_1 (zm, factor, ub, eps)
       zp = 2 - zm
    else if (u < 0) then
       ub = 2 * r
       call map_power_1 (zp, factor, ub, eps)
       zm = 2 - zp
    else
       factor = 1 / eps
       zp = 1
       zm = 1
    end if
    y  = zp / 2
    yb = zm / 2
  end subroutine map_power_01

  subroutine map_power_inverse_01 (y, yb, factor, r, eps)
    real(default), intent(in) :: y, yb
    real(default), intent(out) :: r, factor
    real(default), intent(in) :: eps
    real(default) :: ub, zp, zm
    zp = 2 * y
    zm = 2 * yb
    if (zm < zp) then
       call map_power_inverse_1 (zm, factor, ub, eps)
       r = 1 - ub / 2
    else if (zp < zm) then
       call map_power_inverse_1 (zp, factor, ub, eps)
       r = ub / 2
    else
       factor = 1 / eps
       ub = 1
       r = ub / 2
    end if
  end subroutine map_power_inverse_01

  module subroutine sf_channel_write (object, unit)
    class(sf_channel_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    if (allocated (object%map_code)) then
       do i = 1, size (object%map_code)
          select case (object%map_code (i))
          case (SFMAP_NONE)
             write (u, "(1x,A)", advance="no") "-"
          case (SFMAP_SINGLE)
             write (u, "(1x,A)", advance="no") "+"
          case (SFMAP_MULTI_S)
             write (u, "(1x,A)", advance="no") "s"
          case (SFMAP_MULTI_RES, SFMAP_MULTI_SRS)
             write (u, "(1x,A)", advance="no") "r"
          case (SFMAP_MULTI_ONS, SFMAP_MULTI_SON)
             write (u, "(1x,A)", advance="no") "o"
          case (SFMAP_MULTI_EP)
             write (u, "(1x,A)", advance="no") "e"
          case (SFMAP_MULTI_EPR)
             write (u, "(1x,A)", advance="no") "p"
          case (SFMAP_MULTI_EPO)
             write (u, "(1x,A)", advance="no") "q"
          case (SFMAP_MULTI_IP)
             write (u, "(1x,A)", advance="no") "i"
          case (SFMAP_MULTI_IPR)
             write (u, "(1x,A)", advance="no") "i"
          case (SFMAP_MULTI_IPO)
             write (u, "(1x,A)", advance="no") "i"
          case (SFMAP_MULTI_EI)
             write (u, "(1x,A)", advance="no") "i"
          case default
             write (u, "(1x,A)", advance="no") "?"
          end select
       end do
    else
       write (u, "(1x,A)", advance="no") "-"
    end if
    if (allocated (object%multi_mapping)) then
       write (u, "(1x,'/')", advance="no")
       call object%multi_mapping%write (u)
    else
       write (u, *)
    end if
  end subroutine sf_channel_write

  module subroutine sf_channel_init (channel, n_strfun)
    class(sf_channel_t), intent(out) :: channel
    integer, intent(in) :: n_strfun
    allocate (channel%map_code (n_strfun))
    channel%map_code = SFMAP_NONE
  end subroutine sf_channel_init

  module subroutine sf_channel_assign (copy, original)
    class(sf_channel_t), intent(out) :: copy
    type(sf_channel_t), intent(in) :: original
    allocate (copy%map_code (size (original%map_code)))
    copy%map_code = original%map_code
    if (allocated (original%multi_mapping)) then
         allocate (copy%multi_mapping, source = original%multi_mapping)
    end if
  end subroutine sf_channel_assign

  module subroutine allocate_sf_channels (channel, n_channel, n_strfun)
    type(sf_channel_t), dimension(:), intent(out), allocatable :: channel
    integer, intent(in) :: n_channel
    integer, intent(in) :: n_strfun
    integer :: c
    allocate (channel (n_channel))
    do c = 1, n_channel
       call channel(c)%init (n_strfun)
    end do
  end subroutine allocate_sf_channels

  module subroutine sf_channel_activate_mapping (channel, i_sf)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    channel%map_code(i_sf) = SFMAP_SINGLE
  end subroutine sf_channel_activate_mapping

  module function sf_channel_is_single_mapping (channel, i_sf) result (flag)
    class(sf_channel_t), intent(in) :: channel
    integer, intent(in) :: i_sf
    logical :: flag
    flag = channel%map_code(i_sf) == SFMAP_SINGLE
  end function sf_channel_is_single_mapping

  module function sf_channel_is_multi_mapping (channel, i_sf) result (flag)
    class(sf_channel_t), intent(in) :: channel
    integer, intent(in) :: i_sf
    logical :: flag
    select case (channel%map_code(i_sf))
    case (SFMAP_NONE, SFMAP_SINGLE)
       flag = .false.
    case default
       flag = .true.
    end select
  end function sf_channel_is_multi_mapping

  module function sf_channel_get_multi_mapping_n_par (channel) result (n_par)
    class(sf_channel_t), intent(in) :: channel
    integer :: n_par
    if (allocated (channel%multi_mapping)) then
       n_par = channel%multi_mapping%get_n_dim ()
    else
       n_par = 0
    end if
  end function sf_channel_get_multi_mapping_n_par

  module function any_sf_channel_has_mapping (channel) result (flag)
    type(sf_channel_t), dimension(:), intent(in) :: channel
    logical :: flag
    integer :: c
    flag = .false.
    do c = 1, size (channel)
       flag = flag .or. any (channel(c)%map_code /= SFMAP_NONE)
    end do
  end function any_sf_channel_has_mapping

  module subroutine sf_channel_set_par_index (channel, j, i_par)
    class(sf_channel_t), intent(inout) :: channel
    integer, intent(in) :: j
    integer, intent(in) :: i_par
    associate (mapping => channel%multi_mapping)
      if (j >= 1 .and. j <= mapping%get_n_dim ()) then
         if (mapping%get_index (j) == 0) then
            call channel%multi_mapping%set_index (j, i_par)
         else
            call msg_bug ("Structure-function setup: mapping index set twice")
         end if
      else
         call msg_bug ("Structure-function setup: mapping index out of range")
      end if
    end associate
  end subroutine sf_channel_set_par_index


end submodule sf_mappings_s

