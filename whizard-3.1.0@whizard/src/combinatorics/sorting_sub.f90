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

submodule (sorting) sorting_s

  use diagnostics

  implicit none

contains

  module function sort_int (val_in) result (val)
    integer, dimension(:), intent(in) :: val_in
    integer, dimension(size(val_in)) :: val
    val = val_in( order (val_in) )
  end function sort_int

  module function sort_real (val_in) result (val)
    real(default), dimension(:), intent(in) :: val_in
    real(default), dimension(size(val_in)) :: val
    val = val_in( order (val_in) )
  end function sort_real

  module function sort_int_abs (val_in) result (val)
    integer, dimension(:), intent(in) :: val_in
    integer, dimension(size(val_in)) :: val
    val = val_in( order_abs (val_in) )
  end function sort_int_abs

  module function order_int (val) result (idx)
    integer, dimension(:), intent(in) :: val
    integer, dimension(size(val)) :: idx
    integer :: n, i, s, b1, b2, e1, e2
    n = size (idx)
    do i = 1, n
       idx(i) = i
    end do
    s = 1
    do while (s < n)
       do b1 = 1, n-s, 2*s
          b2 = b1 + s
          e1 = b2 - 1
          e2 = min (e1 + s, n)
          call merge (idx(b1:e2), idx(b1:e1), idx(b2:e2), val)
       end do
       s = 2 * s
    end do
  end function order_int

  module function order_real (val) result (idx)
    real(default), dimension(:), intent(in) :: val
    integer, dimension(size(val)) :: idx
    integer :: n, i, s, b1, b2, e1, e2
    n = size (idx)
    do i = 1, n
       idx(i) = i
    end do
    s = 1
    do while (s < n)
       do b1 = 1, n-s, 2*s
          b2 = b1 + s
          e1 = b2 - 1
          e2 = min (e1 + s, n)
          call merge (idx(b1:e2), idx(b1:e1), idx(b2:e2), val)
       end do
       s = 2 * s
    end do
  end function order_real

  module function order_int_abs (val) result (idx)
    integer, dimension(:), intent(in) :: val
    integer, dimension(size(val)) :: idx
    integer :: n, i, s, b1, b2, e1, e2
    n = size (idx)
    do i = 1, n
       idx(i) = i
    end do
    s = 1
    do while (s < n)
       do b1 = 1, n-s, 2*s
          b2 = b1 + s
          e1 = b2 - 1
          e2 = min (e1 + s, n)
          call merge_abs (idx(b1:e2), idx(b1:e1), idx(b2:e2), val)
       end do
       s = 2 * s
    end do
  end function order_int_abs

  module subroutine merge_int (res, src1, src2, val)
    integer, dimension(:), intent(out) :: res
    integer, dimension(:), intent(in) :: src1, src2
    integer, dimension(:), intent(in) :: val
    integer, dimension(size(res)) :: tmp
    integer :: i1, i2, i
    i1 = 1
    i2 = 1
    do i = 1, size (tmp)
       if (val(src1(i1)) <= val(src2(i2))) then
          tmp(i) = src1(i1);  i1 = i1 + 1
          if (i1 > size (src1)) then
             tmp(i+1:) = src2(i2:)
             exit
          end if
       else
          tmp(i) = src2(i2);  i2 = i2 + 1
          if (i2 > size (src2)) then
             tmp(i+1:) = src1(i1:)
             exit
          end if
       end if
    end do
    res = tmp
  end subroutine merge_int

  module subroutine merge_real (res, src1, src2, val)
    integer, dimension(:), intent(out) :: res
    integer, dimension(:), intent(in) :: src1, src2
    real(default), dimension(:), intent(in) :: val
    integer, dimension(size(res)) :: tmp
    integer :: i1, i2, i
    i1 = 1
    i2 = 1
    do i = 1, size (tmp)
       if (val(src1(i1)) <= val(src2(i2))) then
          tmp(i) = src1(i1);  i1 = i1 + 1
          if (i1 > size (src1)) then
             tmp(i+1:) = src2(i2:)
             exit
          end if
       else
          tmp(i) = src2(i2);  i2 = i2 + 1
          if (i2 > size (src2)) then
             tmp(i+1:) = src1(i1:)
             exit
          end if
       end if
    end do
    res = tmp
  end subroutine merge_real

  module subroutine merge_int_abs (res, src1, src2, val)
    integer, dimension(:), intent(out) :: res
    integer, dimension(:), intent(in) :: src1, src2
    integer, dimension(:), intent(in) :: val
    integer, dimension(size(res)) :: tmp
    integer :: i1, i2, i
    i1 = 1
    i2 = 1
    do i = 1, size (tmp)
       if (abs (val(src1(i1))) < abs (val(src2(i2))) .or. &
          (abs (val(src1(i1))) == abs (val(src2(i2))) .and. &
          val(src1(i1)) >= val(src2(i2)))) then
          tmp(i) = src1(i1);  i1 = i1 + 1
          if (i1 > size (src1)) then
             tmp(i+1:) = src2(i2:)
             exit
          end if
       else
          tmp(i) = src2(i2);  i2 = i2 + 1
          if (i2 > size (src2)) then
             tmp(i+1:) = src1(i1:)
             exit
          end if
       end if
    end do
    res = tmp
  end subroutine merge_int_abs


end submodule sorting_s

