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

submodule (permutations) permutations_s

  implicit none



contains

  elemental module subroutine permutation_init (p, size)
    type(permutation_t), intent(inout) :: p
    integer, intent(in) :: size
    integer :: i
    allocate (p%p (size))
    forall (i = 1:size)
       p%p(i) = i
    end forall
  end subroutine permutation_init

  elemental module subroutine permutation_final (p)
    type(permutation_t), intent(inout) :: p
    deallocate (p%p)
  end subroutine permutation_final

  module subroutine permutation_write (p, u)
    type(permutation_t), intent (in) :: p
    integer, intent(in) :: u
    integer :: i
    do i = 1, size (p%p)
       if (size (p%p) < 10) then
          write (u,"(1x,I1)", advance="no") p%p(i)
       else
          write (u,"(1x,I3)", advance="no") p%p(i)
       end if
    end do
    write (u, *)
  end subroutine permutation_write

  elemental module function permutation_size (perm) result (s)
    type(permutation_t), intent(in) :: perm
    integer :: s
    s = size (perm%p)
  end function permutation_size

  elemental module function permute (i, p) result (j)
    integer, intent(in) :: i
    type(permutation_t), intent(in) :: p
    integer :: j
    if (i > 0 .and. i <= size (p%p)) then
       j = p%p(i)
    else
       j = 0
    end if
  end function permute

  elemental module function permutation_ok (perm) result (ok)
    type(permutation_t), intent(in) :: perm
    logical :: ok
    integer :: i
    logical, dimension(:), allocatable :: set
    ok = .true.
    allocate (set (size (perm%p)))
    set = .false.
    do i = 1, size (perm%p)
       ok = (perm%p(i) > 0 .and. perm%p(i) <= size (perm%p))
       if (.not.ok) return
       set(perm%p(i)) = .true.
    end do
    ok = all (set)
  end function permutation_ok

  module subroutine permutation_find (perm, a1, a2)
    type(permutation_t), intent(inout) :: perm
    integer, dimension(:), intent(in) :: a1, a2
    integer :: i, j
    if (allocated (perm%p))  deallocate (perm%p)
    allocate (perm%p (size (a1)))
    do i = 1, size (a1)
       do j = 1, size (a2)
          if (a1(i) == a2(j)) then
             perm%p(i) = j
             exit
          end if
          perm%p(i) = 0
       end do
    end do
  end subroutine permutation_find

  module subroutine permutation_array_make (pa, code)
    type(permutation_t), dimension(:), allocatable, intent(out) :: pa
    integer, dimension(:), intent(in) :: code
    logical, dimension(size(code)) :: mask
    logical, dimension(:,:), allocatable :: imask
    integer, dimension(:), allocatable :: n_i
    type(permutation_t) :: p_init
    type(permutation_t), dimension(:), allocatable :: p_tmp
    integer :: psize, i, j, k, n_different, n, nn_k
    psize = size (code)
    mask = .true.
    n_different = 0
    do i=1, psize
       if (mask(i)) then
          n_different = n_different + 1
          mask = mask .and. (code /= code(i))
       end if
    end do
    allocate (imask(psize, n_different), n_i(n_different))
    mask = .true.
    k = 0
    do i=1, psize
       if (mask(i)) then
          k = k + 1
          imask(:,k) = (code == code(i))
          n_i(k) = factorial (count(imask(:,k)))
          mask = mask .and. (code /= code(i))
       end if
    end do
    n = product (n_i)
    allocate (pa (n))
    call permutation_init (p_init, psize)
    pa(1) = p_init
    nn_k = 1
    do k = 1, n_different
       allocate (p_tmp (n_i(k)))
       do i = nn_k, 1, -1
          call permutation_array_with_mask (p_tmp, imask(:,k), pa(i))
          do j = n_i(k), 1, -1
             pa((i-1)*n_i(k) + j) = p_tmp(j)
          end do
       end do
       deallocate (p_tmp)
       nn_k = nn_k * n_i(k)
    end do
    call permutation_final (p_init)
    deallocate (imask, n_i)
  end subroutine permutation_array_make

  subroutine permutation_array_with_mask (pa, mask, p_init)
    type(permutation_t), dimension(:), intent(inout) :: pa
    logical, dimension(:), intent(in) :: mask
    type(permutation_t), intent(in) :: p_init
    integer :: plen
    integer :: i, ii, j, fac_i, k, x
    integer, dimension(:), allocatable :: index
    plen = size (pa)
    allocate (index(count(mask)))
    ii = 0
    do i = 1, size (mask)
       if (mask(i)) then
          ii = ii + 1
          index(ii) = i
       end if
    end do
    pa = p_init
    ii = 0
    fac_i = 1
    do i = 1, size (mask)
       if (mask(i)) then
          ii = ii + 1
          fac_i = fac_i * ii
          x = permute (i, p_init)
          do j = 1, plen
             k = ii - mod (((j-1)*fac_i)/plen, ii)
             call insert (pa(j), x, k, ii, index)
          end do
       end if
    end do
    deallocate (index)
  contains
    subroutine insert (p, x, k, n, index)
      type(permutation_t), intent(inout) :: p
      integer, intent(in) :: x, k, n
      integer, dimension(:), intent(in) :: index
      integer :: i
      do i = n, k+1, -1
         p%p(index(i)) = p%p(index(i-1))
      end do
      p%p(index(k)) = x
    end subroutine insert
  end subroutine permutation_array_with_mask

  elemental module function factorial (n) result (f)
    integer, intent(in) :: n
    integer :: f
    integer :: i
    f = 1
    do i=2, abs(n)
       f = f*i
    end do
  end function factorial

  module function tc_permute (k, perm, mask_in) result (pk)
    integer(TC), intent(in) :: k, mask_in
    type(permutation_t), intent(in) :: perm
    integer(TC) :: pk
    integer :: i
    pk = iand (k, mask_in)
    do i = 1, size (perm%p)
       if (btest(k,i-1))  pk = ibset (pk, perm%p(i)-1)
    end do
  end function tc_permute

  module function decay_level_complement (k, mask) result (l)
    integer(TC), intent(in) :: k, mask
    integer :: l
    l = min (decay_level_simple (k), &
         &   decay_level_simple (ieor (k, mask)) + 1)
  end function decay_level_complement

  module function decay_level_simple (k) result(l)
    integer(TC), intent(in) :: k
    integer :: l
    integer :: i
    l = 0
    do i=0, bit_size(k)-1
       if (btest(k,i)) l = l+1
    end do
  end function decay_level_simple


end submodule permutations_s

