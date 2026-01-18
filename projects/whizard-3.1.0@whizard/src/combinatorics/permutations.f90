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

module permutations

  use kinds, only: TC

  implicit none
  private

  public :: permutation_t
  public :: permutation_init
  public :: permutation_final
  public :: permutation_write
  public :: permutation_size
  public :: permute
  public :: permutation_ok
  public :: permutation_find
  public :: permutation_array_make
  public :: factorial
  public :: tc_permute
  public :: tc_decay_level

  type :: permutation_t
     private
     integer, dimension(:), allocatable :: p
  end type permutation_t


  interface tc_decay_level
     module procedure decay_level_simple
     module procedure decay_level_complement
  end interface

  interface
    elemental module subroutine permutation_init (p, size)
      type(permutation_t), intent(inout) :: p
      integer, intent(in) :: size
    end subroutine permutation_init
    elemental module subroutine permutation_final (p)
      type(permutation_t), intent(inout) :: p
    end subroutine permutation_final
    module subroutine permutation_write (p, u)
      type(permutation_t), intent (in) :: p
      integer, intent(in) :: u
    end subroutine permutation_write
    elemental module function permutation_size (perm) result (s)
      type(permutation_t), intent(in) :: perm
      integer :: s
    end function permutation_size
    elemental module function permute (i, p) result (j)
      integer, intent(in) :: i
      type(permutation_t), intent(in) :: p
      integer :: j
    end function permute
    elemental module function permutation_ok (perm) result (ok)
      type(permutation_t), intent(in) :: perm
      logical :: ok
    end function permutation_ok
    module subroutine permutation_find (perm, a1, a2)
      type(permutation_t), intent(inout) :: perm
      integer, dimension(:), intent(in) :: a1, a2
    end subroutine permutation_find
    module subroutine permutation_array_make (pa, code)
      type(permutation_t), dimension(:), allocatable, intent(out) :: pa
      integer, dimension(:), intent(in) :: code
    end subroutine permutation_array_make
    elemental module function factorial (n) result (f)
      integer, intent(in) :: n
      integer :: f
    end function factorial
    module function tc_permute (k, perm, mask_in) result (pk)
      integer(TC), intent(in) :: k, mask_in
      type(permutation_t), intent(in) :: perm
      integer(TC) :: pk
    end function tc_permute
    module function decay_level_complement (k, mask) result (l)
      integer(TC), intent(in) :: k, mask
      integer :: l
    end function decay_level_complement
    module function decay_level_simple (k) result(l)
      integer(TC), intent(in) :: k
      integer :: l
    end function decay_level_simple
  end interface

end module permutations
