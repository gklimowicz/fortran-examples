! utils.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module utils
  use kinds
  implicit none
  private
  public :: create_array_pointer
  private :: create_integer_array_pointer
  private :: create_real_array_pointer
  private :: create_integer_array2_pointer
  private :: create_real_array2_pointer
  public :: copy_array_pointer
  private :: copy_integer_array_pointer
  private :: copy_real_array_pointer
  private :: copy_integer_array2_pointer
  private :: copy_real_array2_pointer
  public :: swap
  private :: swap_integer, swap_real
  public :: sort
  private :: sort_real, sort_real_and_real_array, sort_real_and_integer
  public :: outer_product
  public :: factorize, gcd, lcm
  private :: gcd_internal
  public :: find_free_unit
  integer, dimension(13), parameter, private :: &
       PRIMES = (/ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41 /)
  integer, parameter, private :: MIN_UNIT = 11, MAX_UNIT = 99
  interface create_array_pointer
     module procedure &
          create_integer_array_pointer, &
          create_real_array_pointer, &
          create_integer_array2_pointer, &
          create_real_array2_pointer
  end interface
  interface copy_array_pointer
     module procedure &
          copy_integer_array_pointer, &
          copy_real_array_pointer, &
          copy_integer_array2_pointer, &
          copy_real_array2_pointer
  end interface
  interface swap
     module procedure swap_integer, swap_real
  end interface
  interface sort
     module procedure &
          sort_real, sort_real_and_real_array, &
          sort_real_and_integer
  end interface
contains
  pure subroutine create_integer_array_pointer (lhs, n, lb)
    integer, dimension(:), pointer :: lhs
    integer, intent(in) :: n
    integer, intent(in), optional :: lb
    if (associated (lhs)) then
       if (size (lhs) /= n) then
          deallocate (lhs)
          if (present (lb)) then
             allocate (lhs(lb:n+lb-1))
          else
             allocate (lhs(n))
          end if
       end if
    else
       if (present (lb)) then
          allocate (lhs(lb:n+lb-1))
       else
          allocate (lhs(n))
       end if
    end if
    lhs = 0
  end subroutine create_integer_array_pointer
  pure subroutine create_real_array_pointer (lhs, n, lb)
    real(kind=default), dimension(:), pointer :: lhs
    integer, intent(in) :: n
    integer, intent(in), optional :: lb
    if (associated (lhs)) then
       if (size (lhs) /= n) then
          deallocate (lhs)
          if (present (lb)) then
             allocate (lhs(lb:n+lb-1))
          else
             allocate (lhs(n))
          end if
       end if
    else
       if (present (lb)) then
          allocate (lhs(lb:n+lb-1))
       else
          allocate (lhs(n))
       end if
    end if
    lhs = 0
  end subroutine create_real_array_pointer
  pure subroutine create_integer_array2_pointer (lhs, n, lb)
    integer, dimension(:,:), pointer :: lhs
    integer, dimension(:), intent(in) :: n
    integer, dimension(:), intent(in), optional :: lb
    if (associated (lhs)) then
       if (any (ubound (lhs) /= n)) then
          deallocate (lhs)
          if (present (lb)) then
             allocate (lhs(lb(1):n(1)+lb(1)-1,lb(2):n(2)+lb(2)-1))
          else
             allocate (lhs(n(1),n(2)))
          end if
       end if
    else
       if (present (lb)) then
          allocate (lhs(lb(1):n(1)+lb(1)-1,lb(2):n(2)+lb(2)-1))
       else
          allocate (lhs(n(1),n(2)))
       end if
    end if
    lhs = 0
  end subroutine create_integer_array2_pointer
  pure subroutine create_real_array2_pointer (lhs, n, lb)
    real(kind=default), dimension(:,:), pointer :: lhs
    integer, dimension(:), intent(in) :: n
    integer, dimension(:), intent(in), optional :: lb
    if (associated (lhs)) then
       if (any (ubound (lhs) /= n)) then
          deallocate (lhs)
          if (present (lb)) then
             allocate (lhs(lb(1):n(1)+lb(1)-1,lb(2):n(2)+lb(2)-1))
          else
             allocate (lhs(n(1),n(2)))
          end if
       end if
    else
       if (present (lb)) then
          allocate (lhs(lb(1):n(1)+lb(1)-1,lb(2):n(2)+lb(2)-1))
       else
          allocate (lhs(n(1),n(2)))
       end if
    end if
    lhs = 0
  end subroutine create_real_array2_pointer
  pure subroutine copy_integer_array_pointer (lhs, rhs, lb)
    integer, dimension(:), pointer :: lhs
    integer, dimension(:), intent(in) :: rhs
    integer, intent(in), optional :: lb
    call create_integer_array_pointer (lhs, size (rhs), lb)
    lhs = rhs
  end subroutine copy_integer_array_pointer
  pure subroutine copy_real_array_pointer (lhs, rhs, lb)
    real(kind=default), dimension(:), pointer :: lhs
    real(kind=default), dimension(:), intent(in) :: rhs
    integer, intent(in), optional :: lb
    call create_real_array_pointer (lhs, size (rhs), lb)
    lhs = rhs
  end subroutine copy_real_array_pointer
  pure subroutine copy_integer_array2_pointer (lhs, rhs, lb)
    integer, dimension(:,:), pointer :: lhs
    integer, dimension(:,:), intent(in) :: rhs
    integer, dimension(:), intent(in), optional :: lb
    call create_integer_array2_pointer &
         (lhs, (/ size (rhs, dim=1), size (rhs, dim=2) /), lb)
    lhs = rhs
  end subroutine copy_integer_array2_pointer
  pure subroutine copy_real_array2_pointer (lhs, rhs, lb)
    real(kind=default), dimension(:,:), pointer :: lhs
    real(kind=default), dimension(:,:), intent(in) :: rhs
    integer, dimension(:), intent(in), optional :: lb
    call create_real_array2_pointer &
         (lhs, (/ size (rhs, dim=1), size (rhs, dim=2) /), lb)
    lhs = rhs
  end subroutine copy_real_array2_pointer
  elemental subroutine swap_integer (a, b)
    integer, intent(inout) :: a, b
    integer :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swap_integer
  elemental subroutine swap_real (a, b)
    real(kind=default), intent(inout) :: a, b
    real(kind=default) :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swap_real
  pure subroutine sort_real (key, reverse)
    real(kind=default), dimension(:), intent(inout) :: key
    logical, intent(in), optional :: reverse
    logical :: rev
    integer :: i, j
    if (present (reverse)) then
       rev = reverse
    else
       rev = .false.
    end if
    do i = 1, size (key) - 1
       if (rev) then
          j = sum (maxloc (key(i:))) + i - 1
       else
          j = sum (minloc (key(i:))) + i - 1
       end if
       if (j /= i) then
          call swap (key(i), key(j))
       end if
    end do
  end subroutine sort_real
  pure subroutine sort_real_and_real_array (key, table, reverse)
    real(kind=default), dimension(:), intent(inout) :: key
    real(kind=default), dimension(:,:), intent(inout) :: table
    logical, intent(in), optional :: reverse
    logical :: rev
    integer :: i, j
    if (present (reverse)) then
       rev = reverse
    else
       rev = .false.
    end if
    do i = 1, size (key) - 1
       if (rev) then
          j = sum (maxloc (key(i:))) + i - 1
       else
          j = sum (minloc (key(i:))) + i - 1
       end if
       if (j /= i) then
          call swap (key(i), key(j))
          call swap (table(:,i), table(:,j))
       end if
    end do
  end subroutine sort_real_and_real_array
  pure subroutine sort_real_and_integer (key, table, reverse)
    real(kind=default), dimension(:), intent(inout) :: key
    integer, dimension(:), intent(inout) :: table
    logical, intent(in), optional :: reverse
    logical :: rev
    integer :: i, j
    if (present (reverse)) then
       rev = reverse
    else
       rev = .false.
    end if
    do i = 1, size (key) - 1
       if (rev) then
          j = sum (maxloc (key(i:))) + i - 1
       else
          j = sum (minloc (key(i:))) + i - 1
       end if
       if (j /= i) then
          call swap (key(i), key(j))
          call swap (table(i), table(j))
       end if
    end do
  end subroutine sort_real_and_integer
  pure function outer_product (x, y) result (xy)
    real(kind=default), dimension(:), intent(in) :: x, y
    real(kind=default), dimension(size(x),size(y)) :: xy
    xy = spread (x, dim=2, ncopies=size(y)) &
           * spread (y, dim=1, ncopies=size(x))
  end function outer_product
  pure recursive function gcd_internal (m, n) result (gcd_m_n)
    integer, intent(in) :: m, n
    integer :: gcd_m_n
    if (n <= 0) then
       gcd_m_n = m
    else
       gcd_m_n = gcd_internal (n, modulo (m, n))
    end if
  end function gcd_internal
  elemental function gcd (m, n) result (gcd_m_n)
    integer, intent(in) :: m, n
    integer :: gcd_m_n
    gcd_m_n = gcd_internal (m, n)
  end function gcd
  elemental function lcm (m, n) result (lcm_m_n)
    integer, intent(in) :: m, n
    integer :: lcm_m_n
    lcm_m_n = (m * n) / gcd (m, n)
  end function lcm
  pure subroutine factorize (n, factors, i)
    integer, intent(in) :: n
    integer, dimension(:), intent(out) :: factors
    integer, intent(out) :: i
    integer :: nn, p
    nn = n
    i = 0
    do p = 1, size (PRIMES)
       try: do
          if (modulo (nn, PRIMES(p)) == 0) then
             i = i + 1
             factors(i) = PRIMES(p)
             nn = nn / PRIMES(p)
             if (i >= size (factors)) then
                factors(i) = nn
                return
             end if
          else
             exit try
          end if
       end do try
       if (nn == 1) then
          return
       end if
    end do
  end subroutine factorize
  subroutine find_free_unit (u, iostat)
    integer, intent(out) :: u
    integer, intent(out), optional :: iostat
    logical :: exists, is_open
    integer :: i, status
    do i = MIN_UNIT, MAX_UNIT
       inquire (unit = i, exist = exists, opened = is_open, &
                iostat = status)
       if (status == 0) then
          if (exists .and. .not. is_open) then
             u = i
             if (present (iostat)) then
                iostat = 0
             end if
             return
          end if
       end if
    end do
    if (present (iostat)) then
       iostat = -1
    end if
    u = -1
  end subroutine find_free_unit
end module utils
