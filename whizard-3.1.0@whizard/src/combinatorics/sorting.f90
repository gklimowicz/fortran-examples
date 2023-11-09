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

module sorting

  use kinds, only: default

  implicit none
  private

  public :: sort
  public :: sort_abs
  public :: order
  public :: order_abs

  interface sort
     module procedure sort_int
     module procedure sort_real
  end interface

  interface sort_abs
     module procedure sort_int_abs
  end interface

  interface order
     module procedure order_int
     module procedure order_real
  end interface

  interface order_abs
     module procedure order_int_abs
  end interface

  interface merge
     module procedure merge_int
     module procedure merge_real
  end interface

  interface merge_abs
     module procedure merge_int_abs
  end interface


  interface
    module function sort_int (val_in) result (val)
      integer, dimension(:), intent(in) :: val_in
      integer, dimension(size(val_in)) :: val
    end function sort_int
    module function sort_real (val_in) result (val)
      real(default), dimension(:), intent(in) :: val_in
      real(default), dimension(size(val_in)) :: val
    end function sort_real
    module function sort_int_abs (val_in) result (val)
      integer, dimension(:), intent(in) :: val_in
      integer, dimension(size(val_in)) :: val
    end function sort_int_abs
    module function order_int (val) result (idx)
      integer, dimension(:), intent(in) :: val
      integer, dimension(size(val)) :: idx
    end function order_int
    module function order_real (val) result (idx)
      real(default), dimension(:), intent(in) :: val
      integer, dimension(size(val)) :: idx
    end function order_real
    module function order_int_abs (val) result (idx)
      integer, dimension(:), intent(in) :: val
      integer, dimension(size(val)) :: idx
    end function order_int_abs
    module subroutine merge_int (res, src1, src2, val)
      integer, dimension(:), intent(out) :: res
      integer, dimension(:), intent(in) :: src1, src2
      integer, dimension(:), intent(in) :: val
    end subroutine merge_int
    module subroutine merge_real (res, src1, src2, val)
      integer, dimension(:), intent(out) :: res
      integer, dimension(:), intent(in) :: src1, src2
      real(default), dimension(:), intent(in) :: val
    end subroutine merge_real
    module subroutine merge_int_abs (res, src1, src2, val)
      integer, dimension(:), intent(out) :: res
      integer, dimension(:), intent(in) :: src1, src2
      integer, dimension(:), intent(in) :: val
    end subroutine merge_int_abs
  end interface

end module sorting
