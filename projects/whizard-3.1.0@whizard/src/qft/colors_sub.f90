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

submodule (colors) colors_s

  use io_units
  use diagnostics

  implicit none

contains

  pure module subroutine color_init_trivial (col)
    class(color_t), intent(inout) :: col
    col%defined = .true.
    col%c1 = 0
    col%c2 = 0
    col%ghost = .false.
  end subroutine color_init_trivial

  pure module subroutine color_init_trivial_ghost (col, ghost)
    class(color_t), intent(inout) :: col
    logical, intent(in) :: ghost
    col%defined = .true.
    col%c1 = 0
    col%c2 = 0
    col%ghost = ghost
  end subroutine color_init_trivial_ghost

  pure module subroutine color_init_array (col, c1)
    class(color_t), intent(inout) :: col
    integer, dimension(:), intent(in) :: c1
    col%defined = .true.
    col%c1 = pack (c1, c1 /= 0, [0,0])
    col%c2 = col%c1
    col%ghost = .false.
  end subroutine color_init_array

  pure module subroutine color_init_array_ghost (col, c1, ghost)
    class(color_t), intent(inout) :: col
    integer, dimension(:), intent(in) :: c1
    logical, intent(in) :: ghost
    call color_init_array (col, c1)
    col%ghost = ghost
  end subroutine color_init_array_ghost

  pure module subroutine color_init_arrays (col, c1, c2)
    class(color_t), intent(inout) :: col
    integer, dimension(:), intent(in) :: c1, c2
    col%defined = .true.
    if (size (c1) == size (c2)) then
       col%c1 = pack (c1, c1 /= 0, [0,0])
       col%c2 = pack (c2, c2 /= 0, [0,0])
    else if (size (c1) /= 0) then
       col%c1 = pack (c1, c1 /= 0, [0,0])
       col%c2 = col%c1
    else if (size (c2) /= 0) then
       col%c1 = pack (c2, c2 /= 0, [0,0])
       col%c2 = col%c1
    end if
    col%ghost = .false.
  end subroutine color_init_arrays

  pure module subroutine color_init_arrays_ghost (col, c1, c2, ghost)
    class(color_t), intent(inout) :: col
    integer, dimension(:), intent(in) :: c1, c2
    logical, intent(in) :: ghost
    call color_init_arrays (col, c1, c2)
    col%ghost = ghost
  end subroutine color_init_arrays_ghost

  elemental module subroutine color_init_col_acl (col, col_in, acl_in)
    class(color_t), intent(inout) :: col
    integer, intent(in) :: col_in, acl_in
    integer, dimension(0) :: null_array
    select case (col_in)
    case (0)
       select case (acl_in)
       case (0)
          call color_init_array (col, null_array)
       case default
          call color_init_array (col, [-acl_in])
       end select
    case default
       select case (acl_in)
       case (0)
          call color_init_array (col, [col_in])
       case default
          call color_init_array (col, [col_in, -acl_in])
       end select
    end select
  end subroutine color_init_col_acl

  pure module subroutine color_init_from_array1 (col, c1)
    type(color_t), intent(inout) :: col
    integer, dimension(:), intent(in) :: c1
    logical, dimension(size(c1)) :: mask
    mask = c1 /= 0
    col%defined = .true.
    col%c1 = pack (c1, mask, col%c1)
    col%c2 = col%c1
    col%ghost = .false.
  end subroutine color_init_from_array1

  pure module subroutine color_init_from_array1g (col, c1, ghost)
    type(color_t), intent(inout) :: col
    integer, dimension(:), intent(in) :: c1
    logical, intent(in) :: ghost
    call color_init_from_array1 (col, c1)
    col%ghost = ghost
  end subroutine color_init_from_array1g

  pure module subroutine color_init_from_array2 (col, c1)
    integer, dimension(:,:), intent(in) :: c1
    type(color_t), dimension(:), intent(inout) :: col
    integer :: i
    do i = 1, size (c1,2)
       call color_init_from_array1 (col(i), c1(:,i))
    end do
  end subroutine color_init_from_array2

  pure module subroutine color_init_from_array2g (col, c1, ghost)
    integer, dimension(:,:), intent(in) :: c1
    type(color_t), dimension(:), intent(out) :: col
    logical, intent(in), dimension(:) :: ghost
    call color_init_from_array2 (col, c1)
    col%ghost = ghost
  end subroutine color_init_from_array2g

  elemental module subroutine color_set_ghost (col, ghost)
    class(color_t), intent(inout) :: col
    logical, intent(in) :: ghost
    col%ghost = ghost
  end subroutine color_set_ghost

  elemental module subroutine color_undefine (col, undefine_ghost)
    class(color_t), intent(inout) :: col
    logical, intent(in), optional :: undefine_ghost
    col%defined = .false.
    if (present (undefine_ghost)) then
       if (undefine_ghost)  col%ghost = .false.
    else
       col%ghost = .false.
    end if
  end subroutine color_undefine

  module subroutine color_write_single (col, unit)
    class(color_t), intent(in) :: col
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    if (col%ghost) then
       write (u, "(A)", advance="no")  "c*"
    else if (col%defined) then
       write (u, "(A)", advance="no")  "c("
       if (col%c1(1) /= 0)  write (u, "(I0)", advance="no")  col%c1(1)
       if (any (col%c1 /= 0))  write (u, "(1x)", advance="no")
       if (col%c1(2) /= 0)  write (u, "(I0)", advance="no")  col%c1(2)
       if (.not. col%is_diagonal ()) then
          write (u, "(A)", advance="no")  "|"
          if (col%c2(1) /= 0)  write (u, "(I0)", advance="no")  col%c2(1)
          if (any (col%c2 /= 0))  write (u, "(1x)", advance="no")
          if (col%c2(2) /= 0)  write (u, "(I0)", advance="no")  col%c2(2)
       end if
       write (u, "(A)", advance="no") ")"
    end if
  end subroutine color_write_single

  module subroutine color_write_array (col, unit)
    type(color_t), dimension(:), intent(in) :: col
    integer, intent(in), optional :: unit
    integer :: u
    integer :: i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A)", advance="no") "["
    do i = 1, size (col)
       if (i > 1)  write (u, "(1x)", advance="no")
       call color_write_single (col(i), u)
    end do
    write (u, "(A)", advance="no") "]"
  end subroutine color_write_array

  module subroutine color_write_raw (col, u)
    class(color_t), intent(in) :: col
    integer, intent(in) :: u
    logical :: defined
    defined = col%is_defined () .or. col%is_ghost ()
    write (u) defined
    if (defined) then
       write (u) col%c1, col%c2
       write (u) col%ghost
    end if
  end subroutine color_write_raw

  module subroutine color_read_raw (col, u, iostat)
    class(color_t), intent(inout) :: col
    integer, intent(in) :: u
    integer, intent(out), optional :: iostat
    logical :: defined
    read (u, iostat=iostat) col%defined
    if (col%defined) then
       read (u, iostat=iostat) col%c1, col%c2
       read (u, iostat=iostat) col%ghost
    end if
  end subroutine color_read_raw

  elemental module function color_is_defined (col) result (defined)
    logical :: defined
    class(color_t), intent(in) :: col
    defined = col%defined
  end function color_is_defined

  elemental module function color_is_nonzero (col) result (flag)
    logical :: flag
    class(color_t), intent(in) :: col
    flag = col%defined &
         .and. .not. col%ghost &
         .and. any (col%c1 /= 0 .or. col%c2 /= 0)
  end function color_is_nonzero

  elemental module function color_is_diagonal (col) result (diagonal)
    logical :: diagonal
    class(color_t), intent(in) :: col
    if (col%defined) then
       diagonal = all (col%c1 == col%c2)
    else
       diagonal = .true.
    end if
  end function color_is_diagonal

  elemental module function color_is_ghost (col) result (ghost)
    logical :: ghost
    class(color_t), intent(in) :: col
    ghost = col%ghost
  end function color_is_ghost

  pure function color_ghost_parity (col) result (parity)
    type(color_t), dimension(:), intent(in) :: col
    logical :: parity
    parity = mod (count (col%ghost), 2) == 1
  end function color_ghost_parity

  elemental module function color_get_type (col) result (ctype)
    class(color_t), intent(in) :: col
    integer :: ctype
    if (col%defined) then
       ctype = -1
       if (col%ghost) then
          if (all (col%c1 == 0 .and. col%c2 == 0)) then
             ctype = 8
          end if
       else
          if (all ((col%c1 == 0 .and. col%c2 == 0) &
               & .or. (col%c1 > 0 .and. col%c2 > 0) &
               & .or. (col%c1 < 0 .and. col%c2 < 0))) then
             if (all (col%c1 == 0)) then
                ctype = 1
             else if ((col%c1(1) > 0 .and. col%c1(2) == 0)) then
                ctype = 3
             else if ((col%c1(1) < 0 .and. col%c1(2) == 0)) then
                ctype = -3
             else if ((col%c1(1) > 0 .and. col%c1(2) < 0) &
                  .or.(col%c1(1) < 0 .and. col%c1(2) > 0)) then
                ctype = 8
             end if
          end if
       end if
    else
       ctype = 0
    end if
  end function color_get_type

  elemental module function color_get_number_of_indices (col) result (n)
    integer :: n
    class(color_t), intent(in) :: col
    if (col%defined .and. .not. col%ghost) then
       n = count (col%c1 /= 0)
    else
       n = 0
    end if
  end function color_get_number_of_indices

  elemental module function color_get_col (col) result (c)
    integer :: c
    class(color_t), intent(in) :: col
    integer :: i
    if (col%defined .and. .not. col%ghost) then
       do i = 1, size (col%c1)
          if (col%c1(i) > 0) then
             c = col%c1(i)
             return
          end if
       end do
    end if
    c = 0
  end function color_get_col

  elemental module function color_get_acl (col) result (c)
    integer :: c
    class(color_t), intent(in) :: col
    integer :: i
    if (col%defined .and. .not. col%ghost) then
       do i = 1, size (col%c1)
          if (col%c1(i) < 0) then
             c = - col%c1(i)
             return
          end if
       end do
    end if
    c = 0
  end function color_get_acl

  elemental module function color_get_max_value0 (col) result (cmax)
    integer :: cmax
    type(color_t), intent(in) :: col
    if (col%defined .and. .not. col%ghost) then
       cmax = maxval (abs (col%c1))
    else
       cmax = 0
    end if
  end function color_get_max_value0

  pure module function color_get_max_value1 (col) result (cmax)
    integer :: cmax
    type(color_t), dimension(:), intent(in) :: col
    cmax = maxval (color_get_max_value0 (col))
  end function color_get_max_value1

  pure module function color_get_max_value2 (col) result (cmax)
    integer :: cmax
    type(color_t), dimension(:,:), intent(in) :: col
    integer, dimension(size(col, 2)) :: cm
    integer :: i
    forall (i = 1:size(col, 2))
       cm(i) = color_get_max_value1 (col(:,i))
    end forall
    cmax = maxval (cm)
  end function color_get_max_value2

  elemental module function color_match (col1, col2) result (eq)
    logical :: eq
    class(color_t), intent(in) :: col1, col2
    if (col1%defined .and. col2%defined) then
       if (col1%ghost .and. col2%ghost) then
          eq = .true.
       else if (.not. col1%ghost .and. .not. col2%ghost) then
          eq = all (col1%c1 == col2%c1) .and. all (col1%c2 == col2%c2)
       else
          eq = .false.
       end if
    else
       eq = .true.
    end if
  end function color_match

  elemental module function color_eq (col1, col2) result (eq)
    logical :: eq
    class(color_t), intent(in) :: col1, col2
    if (col1%defined .and. col2%defined) then
       if (col1%ghost .and. col2%ghost) then
          eq = .true.
       else if (.not. col1%ghost .and. .not. col2%ghost) then
          eq = all (col1%c1 == col2%c1) .and. all (col1%c2 == col2%c2)
       else
          eq = .false.
       end if
    else if (.not. col1%defined &
       .and. .not. col2%defined) then
       eq = col1%ghost .eqv. col2%ghost
    else
       eq = .false.
    end if
  end function color_eq

  elemental module function color_neq (col1, col2) result (neq)
    logical :: neq
    class(color_t), intent(in) :: col1, col2
    if (col1%defined .and. col2%defined) then
       if (col1%ghost .and. col2%ghost) then
          neq = .false.
       else if (.not. col1%ghost .and. .not. col2%ghost) then
          neq = any (col1%c1 /= col2%c1) .or. any (col1%c2 /= col2%c2)
       else
          neq = .true.
       end if
    else if (.not. col1%defined &
         .and. .not. col2%defined) then
       neq = col1%ghost .neqv. col2%ghost
    else
       neq = .true.
    end if
  end function color_neq

  elemental module subroutine color_add_offset (col, offset)
    class(color_t), intent(inout) :: col
    integer, intent(in) :: offset
    if (col%defined .and. .not. col%ghost) then
       where (col%c1 /= 0)  col%c1 = col%c1 + sign (offset, col%c1)
       where (col%c2 /= 0)  col%c2 = col%c2 + sign (offset, col%c2)
    end if
  end subroutine color_add_offset

  module subroutine color_canonicalize (col)
    type(color_t), dimension(:), intent(inout) :: col
    integer, dimension(2*size(col)) :: map
    integer :: n_col, i, j, k
    n_col = 0
    do i = 1, size (col)
       if (col(i)%defined .and. .not. col(i)%ghost) then
          do j = 1, size (col(i)%c1)
             if (col(i)%c1(j) /= 0) then
                k = find (abs (col(i)%c1(j)), map(:n_col))
                if (k == 0) then
                   n_col = n_col + 1
                   map(n_col) = abs (col(i)%c1(j))
                   k = n_col
                end if
                col(i)%c1(j) = sign (k, col(i)%c1(j))
             end if
             if (col(i)%c2(j) /= 0) then
                k = find (abs (col(i)%c2(j)), map(:n_col))
                if (k == 0) then
                   n_col = n_col + 1
                   map(n_col) = abs (col(i)%c2(j))
                   k = n_col
                end if
                col(i)%c2(j) = sign (k, col(i)%c2(j))
             end if
          end do
       end if
    end do
  contains
    function find (c, array) result (k)
      integer :: k
      integer, intent(in) :: c
      integer, dimension(:), intent(in) :: array
      integer :: i
      k = 0
      do i = 1, size (array)
         if (c == array (i)) then
            k = i
            return
         end if
      end do
    end function find
  end subroutine color_canonicalize

  subroutine extract_color_line_indices (col, c_index, col_pos)
    type(color_t), dimension(:), intent(in) :: col
    integer, dimension(:), intent(out), allocatable :: c_index
    type(color_t), dimension(size(col)), intent(out) :: col_pos
    integer, dimension(:), allocatable :: c_tmp
    integer :: i, j, k, n, c
    allocate (c_tmp (sum (col%get_number_of_indices ())), source=0)
    n = 0
    SCAN1: do i = 1, size (col)
       if (col(i)%defined .and. .not. col(i)%ghost) then
          SCAN2: do j = 1, 2
             c = abs (col(i)%c1(j))
             if (c /= 0) then
                do k = 1, n
                   if (c_tmp(k) == c) then
                      col_pos(i)%c1(j) = k
                      cycle SCAN2
                   end if
                end do
                n = n + 1
                c_tmp(n) = c
                col_pos(i)%c1(j) = n
             end if
          end do SCAN2
       end if
    end do SCAN1
    allocate (c_index (n))
    c_index = c_tmp(1:n)
  end subroutine extract_color_line_indices

  module subroutine color_array_make_contractions (col_in, col_out)
    type(color_t), dimension(:), intent(in) :: col_in
    type(color_t), dimension(:,:), intent(out), allocatable :: col_out
    type(list_t) :: list
    type(entry_t), pointer :: entry
    integer, dimension(:), allocatable :: c_index
    type(color_t), dimension(size(col_in)) :: col_pos
    integer :: n_prt, n_c_index
    integer, dimension(:), allocatable :: map
    integer :: i, j, c
    n_prt = size (col_in)
    call extract_color_line_indices (col_in, c_index, col_pos)
    n_c_index = size (c_index)
    allocate (map (n_c_index))
    map = 0
    call list_append_if_valid (list, map)
    entry => list%first
    do while (associated (entry))
       do i = 1, n_c_index
          if (entry%map(i) == 0) then
             c = c_index(i)
             do j = i + 1, n_c_index
                if (entry%map(j) == 0) then
                   map = entry%map
                   map(i) = c
                   map(j) = c
                   call list_append_if_valid (list, map)
                end if
             end do
          end if
       end do
       entry => entry%next
    end do
    call list_to_array (list, col_out)
  contains
    subroutine list_append_if_valid (list, map)
      type(list_t), intent(inout) :: list
      integer, dimension(:), intent(in) :: map
      type(entry_t), pointer :: entry
      integer :: i, j, c, p
      entry => list%first
      do while (associated (entry))
         if (all (map == entry%map))  return
         entry => entry%next
      end do
      allocate (entry)
      allocate (entry%map (n_c_index))
      entry%map = map
      allocate (entry%col (n_prt))
      do i = 1, n_prt
         do j = 1, 2
            c = col_in(i)%c1(j)
            if (c /= 0) then
               p = col_pos(i)%c1(j)
               entry%col(i)%defined = .true.
               if (map(p) /= 0) then
                  entry%col(i)%c1(j) = sign (map(p), c)
               else
                  entry%col(i)%c1(j) = c
               endif
               entry%col(i)%c2(j) = entry%col(i)%c1(j)
            end if
         end do
         if (any (entry%col(i)%c1 /= 0) .and. &
              entry%col(i)%c1(1) == - entry%col(i)%c1(2))  return
      end do
      if (associated (list%last)) then
         list%last%next => entry
      else
         list%first => entry
      end if
      list%last => entry
      list%n = list%n + 1
    end subroutine list_append_if_valid
    subroutine list_to_array (list, col)
      type(list_t), intent(inout) :: list
      type(color_t), dimension(:,:), intent(out), allocatable :: col
      type(entry_t), pointer :: entry
      integer :: i
      allocate (col (n_prt, list%n - 1))
      do i = 0, list%n - 1
         entry => list%first
         list%first => list%first%next
         if (i /= 0)  col(:,i) = entry%col
         deallocate (entry)
      end do
      list%last => null ()
    end subroutine list_to_array
  end subroutine color_array_make_contractions

  elemental module subroutine color_invert (col)
    class(color_t), intent(inout) :: col
    if (col%defined .and. .not. col%ghost) then
       col%c1 = - col%c1
       col%c2 = - col%c2
       if (col%c1(1) < 0 .and. col%c1(2) > 0) then
          col%c1 = col%c1(2:1:-1)
          col%c2 = col%c2(2:1:-1)
       end if
    end if
  end subroutine color_invert

  module subroutine color_make_color_map (map, col1, col2)
    integer, dimension(:,:), intent(out), allocatable :: map
    type(color_t), dimension(:), intent(in) :: col1, col2
    integer, dimension(:,:), allocatable :: map1
    integer :: i, j, k
    allocate (map1 (2, 2 * sum (col1%get_number_of_indices ())))
    k = 0
    do i = 1, size (col1)
       if (col1(i)%defined .and. .not. col1(i)%ghost) then
          do j = 1, size (col1(i)%c1)
             if (col1(i)%c1(j) /= 0 &
                  .and. all (map1(1,:k) /= abs (col1(i)%c1(j)))) then
                k = k + 1
                map1(1,k) = abs (col1(i)%c1(j))
                map1(2,k) = abs (col2(i)%c1(j))
             end if
             if (col1(i)%c2(j) /= 0 &
                  .and. all (map1(1,:k) /= abs (col1(i)%c2(j)))) then
                k = k + 1
                map1(1,k) = abs (col1(i)%c2(j))
                map1(2,k) = abs (col2(i)%c2(j))
             end if
          end do
       end if
    end do
    allocate (map (2, k))
    map(:,:) = map1(:,:k)
  end subroutine color_make_color_map

  module subroutine color_translate0 (col, map)
    type(color_t), intent(inout) :: col
    integer, dimension(:,:), intent(in) :: map
    type(color_t) :: col_tmp
    integer :: i
    if (col%defined .and. .not. col%ghost) then
       col_tmp = col
       do i = 1, size (map,2)
          where (abs (col%c1) == map(1,i))
             col_tmp%c1 = sign (map(2,i), col%c1)
          end where
          where (abs (col%c2) == map(1,i))
             col_tmp%c2 = sign (map(2,i), col%c2)
          end where
       end do
       col = col_tmp
    end if
  end subroutine color_translate0

  module subroutine color_translate0_offset (col, map, offset)
    type(color_t), intent(inout) :: col
    integer, dimension(:,:), intent(in) :: map
    integer, intent(in) :: offset
    logical, dimension(size(col%c1)) :: mask1, mask2
    type(color_t) :: col_tmp
    integer :: i
    if (col%defined .and. .not. col%ghost) then
       col_tmp = col
       mask1 = col%c1 /= 0
       mask2 = col%c2 /= 0
       do i = 1, size (map,2)
          where (abs (col%c1) == map(1,i))
             col_tmp%c1 = sign (map(2,i), col%c1)
             mask1 = .false.
          end where
          where (abs (col%c2) == map(1,i))
             col_tmp%c2 = sign (map(2,i), col%c2)
             mask2 = .false.
          end where
       end do
       col = col_tmp
       where (mask1)  col%c1 = sign (abs (col%c1) + offset, col%c1)
       where (mask2)  col%c2 = sign (abs (col%c2) + offset, col%c2)
    end if
  end subroutine color_translate0_offset

  module subroutine color_translate1 (col, map, offset)
    type(color_t), dimension(:), intent(inout) :: col
    integer, dimension(:,:), intent(in) :: map
    integer, intent(in), optional :: offset
    integer :: i
    if (present (offset)) then
       do i = 1, size (col)
          call color_translate0_offset (col(i), map, offset)
       end do
    else
       do i = 1, size (col)
          call color_translate0 (col(i), map)
       end do
    end if
  end subroutine color_translate1

  elemental module function merge_colors (col1, col2) result (col)
    type(color_t) :: col
    class(color_t), intent(in) :: col1, col2
    if (color_is_defined (col1) .and. color_is_defined (col2)) then
       if (color_is_ghost (col1) .and. color_is_ghost (col2)) then
          call color_init_trivial_ghost (col, .true.)
       else
          call color_init_arrays (col, col1%c1, col2%c1)
       end if
    else if (color_is_defined (col1)) then
       call color_init_array (col, col1%c1)
    else if (color_is_defined (col2)) then
       call color_init_array (col, col2%c1)
    end if
  end function merge_colors

  module function color_fusion (col1, col2) result (col)
    class(color_t), intent(in) :: col1, col2
    type(color_t) :: col
    integer, dimension(2) :: ctype
    if (col1%is_defined () .and. col2%is_defined ()) then
       if (col1%is_diagonal () .and. col2%is_diagonal ()) then
          ctype = [col1%get_type (), col2%get_type ()]
          select case (ctype(1))
          case (1)
             select case (ctype(2))
             case (1,3,-3,8)
                col = col2
             end select
          case (3)
             select case (ctype(2))
             case (1)
                col = col1
             case (-3)
                call t_a (col1%get_col (), col2%get_acl ())
             case (8)
                call t_o (col1%get_col (), col2%get_acl (), &
                     &    col2%get_col ())
             end select
          case (-3)
             select case (ctype(2))
             case (1)
                col = col1
             case (3)
                call t_a (col2%get_col (), col1%get_acl ())
             case (8)
                call a_o (col1%get_acl (), col2%get_col (), &
                     &    col2%get_acl ())
             end select
          case (8)
             select case (ctype(2))
             case (1)
                col = col1
             case (3)
                call t_o (col2%get_col (), col1%get_acl (), &
                     &    col1%get_col ())
             case (-3)
                call a_o (col2%get_acl (), col1%get_col (), &
                     &    col1%get_acl ())
             case (8)
                call o_o (col1%get_col (), col1%get_acl (), &
                     &    col2%get_col (), col2%get_acl ())
             end select
          end select
       end if
    end if
  contains
    subroutine t_a (c1, c2)
      integer, intent(in) :: c1, c2
      if (c1 == c2) then
         call col%init_col_acl (0, 0)
      else
         call col%init_col_acl (c1, c2)
      end if
    end subroutine t_a
    subroutine t_o (c1, c2, c3)
      integer, intent(in) :: c1, c2, c3
      if (c1 == c2) then
         call col%init_col_acl (c3, 0)
      else if (c2 == 0 .and. c3 == 0) then
         call col%init_col_acl (c1, 0)
      end if
    end subroutine t_o
    subroutine a_o (c1, c2, c3)
      integer, intent(in) :: c1, c2, c3
      if (c1 == c2) then
         call col%init_col_acl (0, c3)
      else if (c2 == 0 .and. c3 == 0) then
         call col%init_col_acl (0, c1)
      end if
    end subroutine a_o
    subroutine o_o (c1, c2, c3, c4)
      integer, intent(in) :: c1, c2, c3, c4
      if (all ([c1,c2,c3,c4] /= 0)) then
         if (c2 == c3 .and. c4 == c1) then
            call col%init_col_acl (0, 0)
         else if (c2 == c3) then
            call col%init_col_acl (c1, c4)
         else if (c4 == c1) then
            call col%init_col_acl (c3, c2)
         end if
      end if
    end subroutine o_o
  end function color_fusion

  module function compute_color_factor (col1, col2, nc) result (factor)
    real(default) :: factor
    type(color_t), dimension(:), intent(in) :: col1, col2
    integer, intent(in), optional :: nc
    type(color_t), dimension(size(col1)) :: col
    integer :: ncol, nloops, nghost
    ncol = 3;  if (present (nc))  ncol = nc
    col = col1 .merge. col2
    nloops = count_color_loops (col)
    nghost = count (col%is_ghost ())
    factor = real (ncol, default) ** (nloops - nghost)
    if (color_ghost_parity (col))  factor = - factor
  end function compute_color_factor

  module function count_color_loops (col) result (count)
    integer :: count
    type(color_t), dimension(:), intent(in) :: col
    type(color_t), dimension(size(col)) :: cc
    integer :: i, n, offset
    cc = col
    n = size (cc)
    offset = n
    call color_add_offset (cc, offset)
    count = 0
    SCAN_LOOPS: do
       do i = 1, n
          if (color_is_nonzero (cc(i))) then
             if (any (cc(i)%c1 > offset)) then
                count = count + 1
                call follow_line1 (pick_new_line (cc(i)%c1, count, 1))
                cycle SCAN_LOOPS
             end if
          end if
       end do
       exit SCAN_LOOPS
    end do SCAN_LOOPS
  contains
    function pick_new_line (c, reset_val, sgn) result (line)
      integer :: line
      integer, dimension(:), intent(inout) :: c
      integer, intent(in) :: reset_val
      integer, intent(in) :: sgn
      integer :: i
      if (any (c == count)) then
         line = count
      else
         do i = 1, size (c)
            if (sign (1, c(i)) == sgn .and. abs (c(i)) > offset) then
               line = c(i)
               c(i) = reset_val
               return
            end if
         end do
         call color_mismatch
      end if
    end function pick_new_line

    subroutine reset_line (c, line)
      integer, dimension(:), intent(inout) :: c
      integer, intent(in) :: line
      integer :: i
      do i = 1, size (c)
         if (c(i) == line) then
            c(i) = 0
            return
         end if
      end do
    end subroutine reset_line

    recursive subroutine follow_line1 (line)
      integer, intent(in) :: line
      integer :: i
      if (line == count) return
      do i = 1, n
         if (any (cc(i)%c1 == -line)) then
            call reset_line (cc(i)%c1, -line)
            call follow_line2 (pick_new_line (cc(i)%c2, 0, sign (1, -line)))
            return
         end if
      end do
      call color_mismatch ()
    end subroutine follow_line1

    recursive subroutine follow_line2 (line)
      integer, intent(in) :: line
      integer :: i
      do i = 1, n
         if (any (cc(i)%c2 == -line)) then
            call reset_line (cc(i)%c2, -line)
            call follow_line1 (pick_new_line (cc(i)%c1, 0, sign (1, -line)))
            return
         end if
      end do
      call color_mismatch ()
    end subroutine follow_line2

    subroutine color_mismatch ()
      call color_write (col)
      print *
      call msg_fatal ("Color flow mismatch: Non-closed color lines appear during ", &
         [var_str ("the evaluation of color correlations. This can happen if there "), &
          var_str ("are different color structures in the initial or final state of "), &
          var_str ("the process definition. If so, please use separate processes for "), &
          var_str ("the different initial / final states. In a future WHIZARD version "), &
          var_str ("this will be fixed.")])
    end subroutine color_mismatch
  end function count_color_loops


end submodule colors_s

