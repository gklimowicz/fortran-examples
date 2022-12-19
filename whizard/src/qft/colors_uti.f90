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

module colors_uti

  use colors

  implicit none
  private

  public :: color_1
  public :: color_2

contains

  subroutine color_1 (u)
    integer, intent(in) :: u
    type(color_t), dimension(4) :: col1, col2, col
    type(color_t), dimension(:), allocatable :: col3
    type(color_t), dimension(:,:), allocatable :: col_array
    integer :: count, i
    call col1%init_col_acl ([1, 0, 2, 3], [0, 1, 3, 2])
    col2 = col1
    call color_write (col1, u)
    write (u, "(A)")
    call color_write (col2, u)
    write (u, "(A)")
    col = col1 .merge. col2
    call color_write (col, u)
    write (u, "(A)")
    count = count_color_loops (col)
    write (u, "(A,I1)") "Number of color loops (3): ", count
    call col2%init_col_acl ([1, 0, 2, 3], [0, 2, 3, 1])
    call color_write (col1, u)
    write (u, "(A)")
    call color_write (col2, u)
    write (u, "(A)")
    col = col1 .merge. col2
    call color_write (col, u)
    write (u, "(A)")
    count = count_color_loops (col)
    write (u, "(A,I1)")  "Number of color loops (2): ", count
    write (u, "(A)")
    allocate (col3 (4))
    call color_init_from_array (col3, &
         reshape ([1, 0,   0, -1,  2, -3,  3, -2], &
                  [2, 4]))
    call color_write (col3, u)
    write (u, "(A)")
    call color_array_make_contractions (col3, col_array)
    write (u, "(A)")  "Contractions:"
    do i = 1, size (col_array, 2)
       call color_write (col_array(:,i), u)
       write (u, "(A)")
    end do
    deallocate (col3)
    write (u, "(A)")
    allocate (col3 (6))
    call color_init_from_array (col3, &
         reshape ([1, -2,   3, 0,  0, -1,  2, -4,  -3, 0,  4, 0], &
                  [2, 6]))
    call color_write (col3, u)
    write (u, "(A)")
    call color_array_make_contractions (col3, col_array)
    write (u, "(A)")  "Contractions:"
    do i = 1, size (col_array, 2)
       call color_write (col_array(:,i), u)
       write (u, "(A)")
    end do
  end subroutine color_1

  subroutine color_2 (u)
    integer, intent(in) :: u
    type(color_t) :: s1, t1, t2, a1, a2, o1, o2, o3, o4, g1

    write (u, "(A)")  "* Test output: color_2"
    write (u, "(A)")  "*   Purpose: test all combinations for color-object fusion"
    write (u, "(A)")

    call s1%init_col_acl (0,0)
    call t1%init_col_acl (1,0)
    call t2%init_col_acl (2,0)
    call a1%init_col_acl (0,1)
    call a2%init_col_acl (0,2)
    call o1%init_col_acl (1,2)
    call o2%init_col_acl (1,3)
    call o3%init_col_acl (2,3)
    call o4%init_col_acl (2,1)
    call g1%init (ghost=.true.)

    call wrt ("s1", s1)
    call wrt ("t1", t1)
    call wrt ("t2", t2)
    call wrt ("a1", a1)
    call wrt ("a2", a2)
    call wrt ("o1", o1)
    call wrt ("o2", o2)
    call wrt ("o3", o3)
    call wrt ("o4", o4)
    call wrt ("g1", g1)

    write (u, *)

    call wrt ("s1 * s1", s1 .fuse. s1)

    write (u, *)

    call wrt ("s1 * t1", s1 .fuse. t1)
    call wrt ("s1 * a1", s1 .fuse. a1)
    call wrt ("s1 * o1", s1 .fuse. o1)

    write (u, *)

    call wrt ("t1 * s1", t1 .fuse. s1)
    call wrt ("a1 * s1", a1 .fuse. s1)
    call wrt ("o1 * s1", o1 .fuse. s1)

    write (u, *)

    call wrt ("t1 * t1", t1 .fuse. t1)

    write (u, *)

    call wrt ("t1 * t2", t1 .fuse. t2)
    call wrt ("t1 * a1", t1 .fuse. a1)
    call wrt ("t1 * a2", t1 .fuse. a2)
    call wrt ("t1 * o1", t1 .fuse. o1)
    call wrt ("t2 * o1", t2 .fuse. o1)

    write (u, *)

    call wrt ("t2 * t1", t2 .fuse. t1)
    call wrt ("a1 * t1", a1 .fuse. t1)
    call wrt ("a2 * t1", a2 .fuse. t1)
    call wrt ("o1 * t1", o1 .fuse. t1)
    call wrt ("o1 * t2", o1 .fuse. t2)

    write (u, *)

    call wrt ("a1 * a1", a1 .fuse. a1)

    write (u, *)

    call wrt ("a1 * a2", a1 .fuse. a2)
    call wrt ("a1 * o1", a1 .fuse. o1)
    call wrt ("a2 * o2", a2 .fuse. o2)

    write (u, *)

    call wrt ("a2 * a1", a2 .fuse. a1)
    call wrt ("o1 * a1", o1 .fuse. a1)
    call wrt ("o2 * a2", o2 .fuse. a2)

    write (u, *)

    call wrt ("o1 * o1", o1 .fuse. o1)

    write (u, *)

    call wrt ("o1 * o2", o1 .fuse. o2)
    call wrt ("o1 * o3", o1 .fuse. o3)
    call wrt ("o1 * o4", o1 .fuse. o4)

    write (u, *)

    call wrt ("o2 * o1", o2 .fuse. o1)
    call wrt ("o3 * o1", o3 .fuse. o1)
    call wrt ("o4 * o1", o4 .fuse. o1)

    write (u, *)

    call wrt ("g1 * g1", g1 .fuse. g1)

    write (u, *)

    call wrt ("g1 * s1", g1 .fuse. s1)
    call wrt ("g1 * t1", g1 .fuse. t1)
    call wrt ("g1 * a1", g1 .fuse. a1)
    call wrt ("g1 * o1", g1 .fuse. o1)

    write (u, *)

    call wrt ("s1 * g1", s1 .fuse. g1)
    call wrt ("t1 * g1", t1 .fuse. g1)
    call wrt ("a1 * g1", a1 .fuse. g1)
    call wrt ("o1 * g1", o1 .fuse. g1)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: color_2"

  contains

    subroutine wrt (s, col)
      character(*), intent(in) :: s
      class(color_t), intent(in) :: col
      write (u, "(A,1x,'=',1x)", advance="no")  s
      call col%write (u)
      write (u, *)
    end subroutine wrt

  end subroutine color_2


end module colors_uti
