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

submodule (selectors) selectors_s

  use io_units
  use diagnostics
  use format_defs, only: FMT_14, FMT_19

  implicit none

contains

  module subroutine selector_write (object, unit, testflag)
    class(selector_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u, i
    logical :: truncate
    u = given_output_unit (unit)
    truncate = .false.;  if (present (testflag))  truncate = testflag
    write (u, "(1x,A)")  "Selector: i, weight, acc. weight"
    if (allocated (object%weight)) then
       do i = 1, size (object%weight)
          if (truncate) then
             write (u, "(3x,I0,2(1x," // FMT_14 // "))") &
                  object%map(i), object%weight(i), object%acc(i)
          else
             write (u, "(3x,I0,2(1x," // FMT_19 // "))") &
                  object%map(i), object%weight(i), object%acc(i)
          end if
       end do
    else
       write (u, "(3x,A)")  "[undefined]"
    end if
  end subroutine selector_write

  module subroutine selector_init (selector, weight, negative_weights, offset)
    class(selector_t), intent(out) :: selector
    real(default), dimension(:), intent(in) :: weight
    logical, intent(in), optional :: negative_weights
    integer, intent(in), optional :: offset
    real(default) :: s
    integer :: n, i
    logical :: neg_wgt
    logical, dimension(:), allocatable :: mask
    if (present (offset))  selector%offset = offset
    if (size (weight) == 0) &
         call msg_bug ("Selector init: zero-size weight array")
    neg_wgt = .false.
    if (present (negative_weights))  neg_wgt = negative_weights
    if (.not. neg_wgt .and. any (weight < 0)) &
         call msg_fatal ("Selector init: negative weight encountered")
    s = sum (weight)
    allocate (mask (size (weight)), &
         source = weight /= 0)
    n = count (mask)
    if (n > 0) then
       allocate (selector%map (n), &
            source = pack ([(i + selector%offset, i = 1, size (weight))], mask))
       allocate (selector%weight (n), &
            source = pack (abs (weight) / s, mask))
       allocate (selector%acc (n))
       selector%acc(1) = selector%weight(1)
       do i = 2, n - 1
          selector%acc(i) = selector%acc(i-1) + selector%weight(i)
       end do
       selector%acc(n) = 1
    else
       allocate (selector%map (1), source = 1)
       allocate (selector%weight (1), source = 0._default)
       allocate (selector%acc (1), source = 1._default)
    end if
  end subroutine selector_init

  module function selector_select (selector, x) result (n)
    class(selector_t), intent(in) :: selector
    real(default), intent(in) :: x
    integer :: n
    integer :: i
    if (x < 0 .or. x > 1) &
         call msg_bug ("Selector: random number out of range")
    do i = 1, size (selector%acc)
       if (x <= selector%acc(i))  exit
    end do
    n = selector%map(i)
  end function selector_select

  module subroutine selector_generate (selector, rng, n)
    class(selector_t), intent(in) :: selector
    class(rng_t), intent(inout) :: rng
    integer, intent(out) :: n
    real(default) :: x
    select case (size (selector%acc))
    case (1)
       n = selector%map(1)
    case default
       call rng%generate (x)
       n = selector%select (x)
    end select
  end subroutine selector_generate

  module function selector_get_weight (selector, n) result (weight)
    class(selector_t), intent(in) :: selector
    integer, intent(in) :: n
    real(default) :: weight
    integer :: i
    do i = 1, size (selector%weight)
       if (selector%map(i) == n) then
          weight = selector%weight(i)
          return
       end if
    end do
    weight = 0
  end function selector_get_weight


end submodule selectors_s

