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

module su_algebra_uti

  use kinds, only: default
  use physics_defs, only: SCALAR, SPINOR, VECTOR, VECTORSPINOR, TENSOR
  use su_algebra

  implicit none
  private

  public :: su_algebra_1
  public :: su_algebra_2
  public :: su_algebra_3
  public :: su_algebra_4

contains

  subroutine su_algebra_1 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: su_algebra_1"
    write (u, "(A)")  "*   Purpose: test su(N) algebra implementation"
    write (u, "(A)")

    write (u, "(A)")  "* su(N) generators: &
         &list and mark Cartan subalgebra"

    write (u, "(A)")
    write (u, "(A)")  "* s = 0"
    call cartan_check (SCALAR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call cartan_check (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call cartan_check (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call cartan_check (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call cartan_check (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: su_algebra_1"

  contains

    subroutine cartan_check (s)
      integer, intent(in) :: s
      integer :: i
      write (u, *)
      do i = 1, algebra_dimension (s)
         write (u, "(1x,L1)", advance="no")  is_cartan_generator (s, i)
      end do
      write (u, *)

    end subroutine cartan_check

  end subroutine su_algebra_1

  subroutine su_algebra_2 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: su_algebra_2"
    write (u, "(A)")  "*   Purpose: test su(N) algebra implementation"
    write (u, "(A)")

    write (u, "(A)")  "* diagonal su(N) generators: &
         &show explicit representation"
    write (u, "(A)")  "* and check trace and Killing form"

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call cartan_show (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call cartan_show (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call cartan_show (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call cartan_show (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: su_algebra_2"

  contains

    subroutine cartan_show (s)
      integer, intent(in) :: s
      real(default), dimension(:,:), allocatable :: rd
      integer, dimension(:), allocatable :: ci
      integer :: n, d, h, i, j, k, l

      n = algebra_dimension (s)
      d = fundamental_dimension (s)

      write (u, *)
      write (u, "(A2,5X)", advance="no")  "h:"
      do i = 1, d
         j = helicity_index (s, helicity_value (s, i))
         write (u, "(1x,I2,5X)", advance="no")  helicity_value (s, j)
      end do
      write (u, "(8X)", advance="no")
      write (u, "(1X,A)")  "tr"

      allocate (rd (n,d), source = 0._default)
      do i = 1, d
         h = helicity_value (s, i)
         rd(:,i) = cartan_element (s, h)
      end do

      allocate (ci (d-1), source = 0)
      do k = 1, d-1
         ci(k) = cartan_index (s, k)
      end do

      write (u, *)
      do k = 1, d-1
         write (u, "('T',I2,':',1X)", advance="no")  ci(k)
         do i = 1, d
            write (u, 1, advance="no")  rd(ci(k),i)
         end do
         write (u, "(8X)", advance="no")
         write (u, 1)  sum (rd(ci(k),:))
      end do

      write (u, *)
      write (u, "(6X)", advance="no")
      do k = 1, d-1
         write (u, "(2X,'T',I2,3X)", advance="no")  ci(k)
      end do
      write (u, *)

      do k = 1, d-1
         write (u, "('T',I2,2X)", advance="no")  ci(k)
         do l = 1, d-1
            write (u, 1, advance="no")  dot_product (rd(ci(k),:), rd(ci(l),:))
         end do
         write (u, *)
      end do

1     format (1x,F7.4)

    end subroutine cartan_show

  end subroutine su_algebra_2

  subroutine su_algebra_3 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: su_algebra_3"
    write (u, "(A)")  "*   Purpose: test su(N) algebra implementation"
    write (u, "(A)")

    write (u, "(A)")  "* diagonal su(N) generators: &
         &transform to matrix and back"

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call cartan_expand (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call cartan_expand (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call cartan_expand (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call cartan_expand (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: su_algebra_3"

  contains

    subroutine cartan_expand (s)
      integer, intent(in) :: s
      real(default), dimension(:,:), allocatable :: rd
      integer, dimension(:), allocatable :: ci
      real(default), dimension(:), allocatable :: a
      logical, dimension(:), allocatable :: mask
      integer :: n, d, h, i, k, l

      n = algebra_dimension (s)
      d = fundamental_dimension (s)

      allocate (rd (n,d), source = 0._default)
      do i = 1, d
         h = helicity_value (s, i)
         rd(:,i) = cartan_element (s, h)
      end do

      allocate (ci (d-1), source = 0)
      do k = 1, d-1
         ci(k) = cartan_index (s, k)
      end do

      allocate (a (n))

      write (u, *)
      do k = 1, d-1
         a(:) = cartan_coeff (s, rd(ci(k),:))
         write (u, "('T',I2,':',1X)", advance="no")  ci(k)
         do i = 1, n
            if (is_cartan_generator (s, i)) then
               write (u, 1, advance="no")  a(i)
            else if (a(i) /= 0) then
               ! this should not happen (nonzero non-Cartan entry)
               write (u, "(1X,':',I2,':',3X)", advance="no")  i
            end if
         end do
         write (u, *)
      end do

1     format (1X,F7.4)

    end subroutine cartan_expand

  end subroutine su_algebra_3

  subroutine su_algebra_4 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: su_algebra_4"
    write (u, "(A)")  "*   Purpose: test su(N) algebra implementation"
    write (u, "(A)")

    write (u, "(A)")  "* off-diagonal su(N) generators: &
         &mapping from/to helicity pair"

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call root_expand (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call root_expand (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call root_expand (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call root_expand (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: su_algebra_4"

  contains

    subroutine root_expand (s)
      integer, intent(in) :: s
      integer :: n, d, i, j, h1, h2
      logical :: r

      n = algebra_dimension (s)

      write (u, *)
      do i = 1, n
         if (is_cartan_generator (s, i))  cycle
         call root_helicity (s, i, h1, h2, r)
         j = root_index (s, h1, h2, r)
         write (u, "('T',I2,':')", advance="no")  j
         write (u, "(2(1x,I2))", advance="no")  h1, h2
         if (r) then
            write (u, *)
         else
            write (u, "('*')")
         end if
      end do

    end subroutine root_expand

  end subroutine su_algebra_4


end module su_algebra_uti
