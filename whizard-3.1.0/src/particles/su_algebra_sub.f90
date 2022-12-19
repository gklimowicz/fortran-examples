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

submodule (su_algebra) su_algebra_s

  use physics_defs, only: SCALAR, SPINOR, VECTOR, VECTORSPINOR, TENSOR

  implicit none

contains

  module function algebra_dimension (s) result (n)
    integer :: n
    integer, intent(in) :: s
    n = fundamental_dimension (s) ** 2 - 1
  end function algebra_dimension

  module function fundamental_dimension (s) result (d)
    integer :: d
    integer, intent(in) :: s
    d = s
  end function fundamental_dimension

  module function helicity_value (s, i) result (h)
    integer :: h
    integer, intent(in) :: s, i
    integer, dimension(1), parameter :: hh1 = [0]
    integer, dimension(2), parameter :: hh2 = [1, -1]
    integer, dimension(3), parameter :: hh3 = [1,  0, -1]
    integer, dimension(4), parameter :: hh4 = [2,  1, -1, -2]
    integer, dimension(5), parameter :: hh5 = [2,  1,  0, -1, -2]
    h = 0
    select case (s)
    case (SCALAR)
       select case (i)
       case (1:1);  h = hh1(i)
       end select
    case (SPINOR)
       select case (i)
       case (1:2);  h = hh2(i)
       end select
    case (VECTOR)
       select case (i)
       case (1:3);  h = hh3(i)
       end select
    case (VECTORSPINOR)
       select case (i)
       case (1:4);  h = hh4(i)
       end select
    case (TENSOR)
       select case (i)
       case (1:5);  h = hh5(i)
       end select
    end select
  end function helicity_value

  module function helicity_index (s, h) result (i)
    integer, intent(in) :: s, h
    integer :: i
    integer, dimension(0:0), parameter :: hi1 = [1]
    integer, dimension(-1:1), parameter :: hi2 = [2, 0, 1]
    integer, dimension(-1:1), parameter :: hi3 = [3, 2, 1]
    integer, dimension(-2:2), parameter :: hi4 = [4, 3, 0, 2, 1]
    integer, dimension(-2:2), parameter :: hi5 = [5, 4, 3, 2, 1]
    select case (s)
    case (SCALAR)
       i = hi1(h)
    case (SPINOR)
       i = hi2(h)
    case (VECTOR)
       i = hi3(h)
    case (VECTORSPINOR)
       i = hi4(h)
    case (TENSOR)
       i = hi5(h)
    end select
  end function helicity_index

  elemental module function is_cartan_generator (s, i) result (cartan)
    logical :: cartan
    integer, intent(in) :: s, i
    select case (s)
    case (SCALAR)
    case (SPINOR)
       select case (i)
       case (3);  cartan = .true.
       case default
          cartan = .false.
       end select
    case (VECTOR)
       select case (i)
       case (3,8);  cartan = .true.
       case default
          cartan = .false.
       end select
    case (VECTORSPINOR)
       select case (i)
       case (3,6,15);  cartan = .true.
       case default
          cartan = .false.
       end select
    case (TENSOR)
       select case (i)
       case (3,6,15,24);  cartan = .true.
       case default
          cartan = .false.
       end select
    case default
       cartan = .false.
    end select
  end function is_cartan_generator

  elemental module function cartan_index (s, k) result (ci)
    integer :: ci
    integer, intent(in) :: s, k
    integer, dimension(1), parameter :: ci2 = [3]
    integer, dimension(2), parameter :: ci3 = [3,8]
    integer, dimension(3), parameter :: ci4 = [3,6,15]
    integer, dimension(4), parameter :: ci5 = [3,6,15,24]
    select case (s)
    case (SPINOR)
       ci = ci2(k)
    case (VECTOR)
       ci = ci3(k)
    case (VECTORSPINOR)
       ci = ci4(k)
    case (TENSOR)
       ci = ci5(k)
    case default
       ci = 0
    end select
  end function cartan_index

  module function cartan_element (s, h) result (a)
    real(default), dimension(:), allocatable :: a
    integer, intent(in) :: s, h
    real(default), parameter :: sqrt2 = sqrt (2._default)
    real(default), parameter :: sqrt3 = sqrt (3._default)
    real(default), parameter :: sqrt10 = sqrt (10._default)
    allocate (a (algebra_dimension (s)), source = 0._default)
    select case (s)
    case (SCALAR)
    case (SPINOR)
       select case (h)
       case (1)
          a(3) =  1._default / 2
       case (-1)
          a(3) = -1._default / 2
       end select
    case (VECTOR)
       select case (h)
       case (1)
          a(3) =  1._default / 2
          a(8) =  1._default / (2 * sqrt3)
       case (-1)
          a(3) = -1._default / 2
          a(8) =  1._default / (2 * sqrt3)
       case (0)
          a(8) = -1._default / sqrt3
       end select
    case (VECTORSPINOR)
       select case (h)
       case (2)
          a(3)  =  1._default / 2
          a(15) =  1._default / (2 * sqrt2)
       case (-2)
          a(3)  = -1._default / 2
          a(15) =  1._default / (2 * sqrt2)
       case (1)
          a(6)  =  1._default / 2
          a(15) = -1._default / (2 * sqrt2)
       case (-1)
          a(6)  = -1._default / 2
          a(15) = -1._default / (2 * sqrt2)
       end select
    case (TENSOR)
       select case (h)
       case (2)
          a(3)  =  1._default / 2
          a(15) =  1._default / (2 * sqrt2)
          a(24) =  1._default / (2 * sqrt10)
       case (-2)
          a(3)  = -1._default / 2
          a(15) =  1._default / (2 * sqrt2)
          a(24) =  1._default / (2 * sqrt10)
       case (1)
          a(6)  =  1._default / 2
          a(15) = -1._default / (2 * sqrt2)
          a(24) =  1._default / (2 * sqrt10)
       case (-1)
          a(6)  = -1._default / 2
          a(15) = -1._default / (2 * sqrt2)
          a(24) =  1._default / (2 * sqrt10)
       case (0)
          a(24) = -4._default / (2 * sqrt10)
       end select
    end select
  end function cartan_element

  module function cartan_coeff (s, rd) result (a)
    real(default), dimension(:), allocatable :: a
    integer, intent(in) :: s
    real(default), dimension(:), intent(in) :: rd
    real(default), parameter :: sqrt2 = sqrt (2._default)
    real(default), parameter :: sqrt3 = sqrt (3._default)
    real(default), parameter :: sqrt10 = sqrt (10._default)
    integer :: n
    n = algebra_dimension (s)
    allocate (a (n), source = 0._default)
    select case (s)
    case (SPINOR)
       a(3) = rd(1) - rd(2)
    case (VECTOR)
       a(3) = rd(1) - rd(3)
       a(8) = (rd(1) - 2 * rd(2) + rd(3)) / sqrt3
    case (VECTORSPINOR)
       a(3) = rd(1) - rd(4)
       a(6) = rd(2) - rd(3)
       a(15) = (rd(1) - rd(2) - rd(3) + rd(4)) / sqrt2
    case (TENSOR)
       a(3) = rd(1) - rd(5)
       a(6) = rd(2) - rd(4)
       a(15) = (rd(1) - rd(2) - rd(4) + rd(5)) / sqrt2
       a(24) = (rd(1) + rd(2) - 4 * rd(3) + rd(4) + rd(5)) / sqrt10
    end select
  end function cartan_coeff

  module function root_index (s, h1, h2, r) result (ai)
    integer :: ai
    integer, intent(in) :: s, h1, h2
    logical :: r
    ai = 0
    select case (s)
    case (SCALAR)
    case (SPINOR)
       select case (h1)
       case (1)
          select case (h2)
          case (-1);  ai = 1
          end select
       end select
    case (VECTOR)
       select case (h1)
       case (1)
          select case (h2)
          case (-1);  ai = 1
          case (0);   ai = 4
          end select
       case (0)
          select case (h2)
          case (-1);  ai = 6
          end select
       end select
    case (VECTORSPINOR)
       select case (h1)
       case (2)
          select case (h2)
          case (-2);  ai = 1
          case (1);   ai = 7
          case (-1);  ai = 11
          end select
       case (1)
          select case (h2)
          case (-1);  ai = 4
          case (-2);  ai = 13
          end select
       case (-1)
          select case (h2)
          case (-2);  ai = 9
          end select
       end select
    case (TENSOR)
       select case (h1)
       case (2)
          select case (h2)
          case (-2);  ai = 1
          case (1);   ai = 7
          case (-1);  ai = 11
          case (0);   ai = 16
          end select
       case (1)
          select case (h2)
          case (-1);  ai = 4
          case (-2);  ai = 13
          case (0);   ai = 20
          end select
       case (-1)
          select case (h2)
          case (-2);  ai = 9
          end select
       case (0)
          select case (h2)
          case (-2);  ai = 18
          case (-1);  ai = 22
          end select
       end select
    end select
    if (ai /= 0 .and. .not. r)  ai = ai + 1
  end function root_index

  module subroutine root_helicity (s, i, h1, h2, r)
    integer, intent(in) :: s, i
    integer, intent(out) :: h1, h2
    logical, intent(out) :: r
    h1 = 0
    h2 = 0
    r  = .false.
    select case (s)
    case (SCALAR)
    case (SPINOR)
       select case (i)
       case ( 1, 2);  h1 =  1;  h2 = -1;  r = i == 1
       end select
    case (VECTOR)
       select case (i)
       case ( 1, 2);  h1 =  1;  h2 = -1;  r = i == 1
       case ( 4, 5);  h1 =  1;  h2 =  0;  r = i == 4
       case ( 6, 7);  h1 =  0;  h2 = -1;  r = i == 6
       end select
    case (VECTORSPINOR)
       select case (i)
       case ( 1, 2);  h1 =  2;  h2 = -2;  r = i == 1
       case ( 4, 5);  h1 =  1;  h2 = -1;  r = i == 4
       case ( 7, 8);  h1 =  2;  h2 =  1;  r = i == 7
       case ( 9,10);  h1 = -1;  h2 = -2;  r = i == 9
       case (11,12);  h1 =  2;  h2 = -1;  r = i ==11
       case (13,14);  h1 =  1;  h2 = -2;  r = i ==13
       end select
    case (TENSOR)
       select case (i)
       case ( 1, 2);  h1 =  2;  h2 = -2;  r = i == 1
       case ( 4, 5);  h1 =  1;  h2 = -1;  r = i == 4
       case ( 7, 8);  h1 =  2;  h2 =  1;  r = i == 7
       case ( 9,10);  h1 = -1;  h2 = -2;  r = i == 9
       case (11,12);  h1 =  2;  h2 = -1;  r = i ==11
       case (13,14);  h1 =  1;  h2 = -2;  r = i ==13
       case (16,17);  h1 =  2;  h2 =  0;  r = i ==16
       case (18,19);  h1 =  0;  h2 = -2;  r = i ==18
       case (20,21);  h1 =  1;  h2 =  0;  r = i ==20
       case (22,23);  h1 =  0;  h2 = -1;  r = i ==22
       end select
    end select
  end subroutine root_helicity


end submodule su_algebra_s

