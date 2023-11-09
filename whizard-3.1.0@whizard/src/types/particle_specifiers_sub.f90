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

submodule (particle_specifiers) particle_specifiers_s

  use io_units
  use diagnostics

  implicit none

contains

  recursive module function prt_expr_to_string (object) result (string)
    class(prt_expr_t), intent(in) :: object
    type(string_t) :: string
    if (allocated (object%x)) then
       string = object%x%to_string ()
    else
       string = ""
    end if
  end function prt_expr_to_string

  module function prt_expr_get_n_terms (object) result (n)
    class(prt_expr_t), intent(in) :: object
    integer :: n
    if (allocated (object%x)) then
       select type (x => object%x)
       type is (prt_spec_sum_t)
          n = size (x%expr)
       class default
          n = 1
       end select
    else
       n = 0
    end if
  end function prt_expr_get_n_terms

  recursive module subroutine prt_expr_term_to_array (object, array, i)
    class(prt_expr_t), intent(in) :: object
    type(prt_spec_t), dimension(:), intent(inout), allocatable :: array
    integer, intent(in) :: i
    integer :: j
    if (allocated (array))  deallocate (array)
    select type (x => object%x)
    type is (prt_spec_t)
       allocate (array (1))
       array(1) = x
    type is (prt_spec_list_t)
       allocate (array (size (x%expr)))
       do j = 1, size (array)
          select type (y => x%expr(j)%x)
          type is (prt_spec_t)
             array(j) = y
          end select
       end do
    type is (prt_spec_sum_t)
       call x%expr(i)%term_to_array (array, 1)
    end select
  end subroutine prt_expr_term_to_array

  module subroutine prt_spec_write1 (object, unit, advance)
    type(prt_spec_t), intent(in) :: object
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: advance
    character(3) :: adv
    integer :: u
    u = given_output_unit (unit)
    adv = "yes";  if (present (advance))  adv = advance
    write (u, "(A)", advance = adv)  char (object%to_string ())
  end subroutine prt_spec_write1

  module subroutine prt_spec_write2 (prt_spec, unit, advance)
    type(prt_spec_t), dimension(:), intent(in) :: prt_spec
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: advance
    character(3) :: adv
    integer :: u, i
    u = given_output_unit (unit)
    adv = "yes";  if (present (advance))  adv = advance
    do i = 1, size (prt_spec)
       if (i > 1)  write (u, "(A)", advance="no")  ", "
       call prt_spec_write (prt_spec(i), u, advance="no")
    end do
    write (u, "(A)", advance = adv)
  end subroutine prt_spec_write2

  pure module subroutine prt_spec_read1 (prt_spec, string)
    type(prt_spec_t), intent(out) :: prt_spec
    type(string_t), intent(in) :: string
    type(string_t) :: arg, buffer
    integer :: b1, b2, c, n, i
    b1 = scan (string, "(")
    b2 = scan (string, ")")
    if (b1 == 0) then
       prt_spec%name = trim (adjustl (string))
    else
       prt_spec%name = trim (adjustl (extract (string, 1, b1-1)))
       arg = trim (adjustl (extract (string, b1+1, b2-1)))
       if (arg == "*") then
          prt_spec%polarized = .true.
       else
          n = 0
          buffer = arg
          do
             if (verify (buffer, " ") == 0)  exit
             n = n + 1
             c = scan (buffer, "+")
             if (c == 0)  exit
             buffer = extract (buffer, c+1)
          end do
          allocate (prt_spec%decay (n))
          buffer = arg
          do i = 1, n
             c = scan (buffer, "+")
             if (c == 0)  c = len (buffer) + 1
             prt_spec%decay(i) = trim (adjustl (extract (buffer, 1, c-1)))
             buffer = extract (buffer, c+1)
          end do
       end if
    end if
  end subroutine prt_spec_read1

  pure module subroutine prt_spec_read2 (prt_spec, string)
    type(prt_spec_t), dimension(:), intent(out), allocatable :: prt_spec
    type(string_t), intent(in) :: string
    type(string_t) :: buffer
    integer :: c, i, n
    n = 0
    buffer = string
    do
       n = n + 1
       c = scan (buffer, ",")
       if (c == 0)  exit
       buffer = extract (buffer, c+1)
    end do
    allocate (prt_spec (n))
    buffer = string
    do i = 1, size (prt_spec)
       c = scan (buffer, ",")
       if (c == 0)  c = len (buffer) + 1
       call prt_spec_read (prt_spec(i), &
            trim (adjustl (extract (buffer, 1, c-1))))
       buffer = extract (buffer, c+1)
    end do
  end subroutine prt_spec_read2

  elemental module function new_prt_spec_ (name) result (prt_spec)
    type(string_t), intent(in) :: name
    type(prt_spec_t) :: prt_spec
    prt_spec%name = name
  end function new_prt_spec_

  elemental module function new_prt_spec_polarized (name, polarized) result (prt_spec)
    type(string_t), intent(in) :: name
    logical, intent(in) :: polarized
    type(prt_spec_t) :: prt_spec
    prt_spec%name = name
    prt_spec%polarized = polarized
  end function new_prt_spec_polarized

  pure module function new_prt_spec_unstable (name, decay) result (prt_spec)
    type(string_t), intent(in) :: name
    type(string_t), dimension(:), intent(in) :: decay
    type(prt_spec_t) :: prt_spec
    prt_spec%name = name
    allocate (prt_spec%decay (size (decay)))
    prt_spec%decay = decay
  end function new_prt_spec_unstable

  elemental module function prt_spec_get_name (prt_spec) result (name)
    class(prt_spec_t), intent(in) :: prt_spec
    type(string_t) :: name
    name = prt_spec%name
  end function prt_spec_get_name

  module function prt_spec_to_string (object) result (string)
    class(prt_spec_t), intent(in) :: object
    type(string_t) :: string
    integer :: i
    string = object%name
    if (allocated (object%decay)) then
       string = string // "("
       do i = 1, size (object%decay)
          if (i > 1)  string = string // " + "
          string = string // object%decay(i)
       end do
       string = string // ")"
    else if (object%polarized) then
       string = string // "(*)"
    end if
  end function prt_spec_to_string

  elemental module function prt_spec_is_polarized (prt_spec) result (flag)
    class(prt_spec_t), intent(in) :: prt_spec
    logical :: flag
    flag = prt_spec%polarized
  end function prt_spec_is_polarized

  elemental module function prt_spec_is_unstable (prt_spec) result (flag)
    class(prt_spec_t), intent(in) :: prt_spec
    logical :: flag
    flag = allocated (prt_spec%decay)
  end function prt_spec_is_unstable

  elemental module function prt_spec_get_n_decays (prt_spec) result (n)
    class(prt_spec_t), intent(in) :: prt_spec
    integer :: n
    if (allocated (prt_spec%decay)) then
       n = size (prt_spec%decay)
    else
       n = 0
    end if
  end function prt_spec_get_n_decays

  module subroutine prt_spec_get_decays (prt_spec, decay)
    class(prt_spec_t), intent(in) :: prt_spec
    type(string_t), dimension(:), allocatable, intent(out) :: decay
    if (allocated (prt_spec%decay)) then
       allocate (decay (size (prt_spec%decay)))
       decay = prt_spec%decay
    else
       allocate (decay (0))
    end if
  end subroutine prt_spec_get_decays

  module subroutine prt_spec_expand_sub (object)
    class(prt_spec_t), intent(inout) :: object
  end subroutine prt_spec_expand_sub

  recursive module function prt_spec_list_to_string (object) result (string)
    class(prt_spec_list_t), intent(in) :: object
    type(string_t) :: string
    integer :: i
    string = ""
    if (allocated (object%expr)) then
       do i = 1, size (object%expr)
          if (i > 1)  string = string // ", "
          select type (x => object%expr(i)%x)
          type is (prt_spec_list_t)
             string = string // "(" // x%to_string () // ")"
          class default
             string = string // x%to_string ()
          end select
       end do
    end if
  end function prt_spec_list_to_string

  module subroutine prt_spec_list_flatten (object)
    class(prt_spec_list_t), intent(inout) :: object
    type(prt_expr_t), dimension(:), allocatable :: tmp_expr
    integer :: i, n_flat, i_flat
    n_flat = 0
    do i = 1, size (object%expr)
       select type (y => object%expr(i)%x)
       type is (prt_spec_list_t)
          n_flat = n_flat + size (y%expr)
       class default
          n_flat = n_flat + 1
       end select
    end do
    if (n_flat > size (object%expr)) then
       allocate (tmp_expr (n_flat))
       i_flat = 0
       do i = 1, size (object%expr)
          select type (y => object%expr(i)%x)
          type is (prt_spec_list_t)
             tmp_expr (i_flat + 1 : i_flat + size (y%expr)) = y%expr
             i_flat = i_flat + size (y%expr)
          class default
             tmp_expr (i_flat + 1) = object%expr(i)
             i_flat = i_flat + 1
          end select
       end do
    end if
    if (allocated (tmp_expr)) &
         call move_alloc (from = tmp_expr, to = object%expr)
  end subroutine prt_spec_list_flatten

  recursive module subroutine prt_spec_list_expand_sub (object)
    class(prt_spec_list_t), intent(inout) :: object
    integer :: i
    if (allocated (object%expr)) then
       do i = 1, size (object%expr)
          call object%expr(i)%expand ()
       end do
    end if
  end subroutine prt_spec_list_expand_sub

  recursive module function prt_spec_sum_to_string (object) result (string)
    class(prt_spec_sum_t), intent(in) :: object
    type(string_t) :: string
    integer :: i
    string = ""
    if (allocated (object%expr)) then
       do i = 1, size (object%expr)
          if (i > 1)  string = string // " + "
          select type (x => object%expr(i)%x)
          type is (prt_spec_list_t)
             string = string // "(" // x%to_string () // ")"
          type is (prt_spec_sum_t)
             string = string // "(" // x%to_string () // ")"
          class default
             string = string // x%to_string ()
          end select
       end do
    end if
  end function prt_spec_sum_to_string

  module subroutine prt_spec_sum_flatten (object)
    class(prt_spec_sum_t), intent(inout) :: object
    type(prt_expr_t), dimension(:), allocatable :: tmp_expr
    integer :: i, n_flat, i_flat
    n_flat = 0
    do i = 1, size (object%expr)
       select type (y => object%expr(i)%x)
       type is (prt_spec_sum_t)
          n_flat = n_flat + size (y%expr)
          class default
          n_flat = n_flat + 1
       end select
    end do
    if (n_flat > size (object%expr)) then
       allocate (tmp_expr (n_flat))
       i_flat = 0
       do i = 1, size (object%expr)
          select type (y => object%expr(i)%x)
          type is (prt_spec_sum_t)
             tmp_expr (i_flat + 1 : i_flat + size (y%expr)) = y%expr
             i_flat = i_flat + size (y%expr)
          class default
             tmp_expr (i_flat + 1) = object%expr(i)
             i_flat = i_flat + 1
          end select
       end do
    end if
    if (allocated (tmp_expr)) &
         call move_alloc (from = tmp_expr, to = object%expr)
  end subroutine prt_spec_sum_flatten

  recursive module subroutine prt_spec_sum_expand_sub (object)
    class(prt_spec_sum_t), intent(inout) :: object
    integer :: i
    if (allocated (object%expr)) then
       do i = 1, size (object%expr)
          call object%expr(i)%expand ()
       end do
    end if
  end subroutine prt_spec_sum_expand_sub


end submodule particle_specifiers_s

