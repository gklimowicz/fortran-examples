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

module grids_uti

  use kinds, only: default
  use constants, only: zero, one, two, three, four, tiny_07
  use file_utils, only: delete_file
  use numeric_utils

  use grids

  implicit none
  private

  public :: grids_1
  public :: grids_2
  public :: grids_3
  public :: grids_4
  public :: grids_5

contains

  subroutine grids_1 (u)
    integer, intent(in) :: u
    type(grid_t) :: grid
    write (u, "(A)")  "* Test output: grids_1"
    write (u, "(A)")  "*   Purpose: Test Index Function"
    write (u, "(A)")

    call grid%init ([3])
    call grid%write(u)
    call assert (u, grid%get_index([1]) == 1, "grid%get_index(1) == 1")
    call assert (u, grid%get_index([2]) == 2, "grid%get_index(2) == 2")
    call assert (u, grid%get_index([3]) == 3, "grid%get_index(3) == 3")
    call grid%final ()

    call grid%init ([3,3])
    call grid%write(u)
    call assert (u, grid%get_index([1,1]) == 1, "grid%get_index(1,1) == 1")
    call assert (u, grid%get_index([2,1]) == 2, "grid%get_index(2,1) == 2")
    call assert (u, grid%get_index([3,1]) == 3, "grid%get_index(3,1) == 3")
    call assert (u, grid%get_index([1,2]) == 4, "grid%get_index(1,2) == 4")
    call assert (u, grid%get_index([2,2]) == 5, "grid%get_index(2,2) == 5")
    call assert (u, grid%get_index([3,2]) == 6, "grid%get_index(3,2) == 6")
    call assert (u, grid%get_index([1,3]) == 7, "grid%get_index(1,3) == 7")
    call assert (u, grid%get_index([2,3]) == 8, "grid%get_index(2,3) == 8")
    call assert (u, grid%get_index([3,3]) == 9, "grid%get_index(3,3) == 9")
    call grid%final ()

    call grid%init ([3,3,2])
    call grid%write(u)
    call assert (u, grid%get_index([1,1,1]) == 1,   "grid%get_index(1,1,1) == 1")
    call assert (u, grid%get_index([2,1,2]) == 2+9, "grid%get_index(2,1,2) == 2+9")
    call assert (u, grid%get_index([3,3,1]) == 9,   "grid%get_index(3,3,1) == 3")
    call assert (u, grid%get_index([3,1,2]) == 3+9, "grid%get_index(3,1,2) == 4+9")
    call assert (u, grid%get_index([2,2,1]) == 5,   "grid%get_index(2,2,1) == 5")
    call assert (u, grid%get_index([3,2,2]) == 6+9, "grid%get_index(3,2,2) == 6+9")
    call assert (u, grid%get_index([1,3,1]) == 7,   "grid%get_index(1,3,1) == 7")
    call assert (u, grid%get_index([2,3,2]) == 8+9, "grid%get_index(2,3,2) == 8+9")
    call assert (u, grid%get_index([3,3,2]) == 9+9, "grid%get_index(3,3,2) == 9+9")
    call grid%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: grids_1"
  end subroutine grids_1

  subroutine grids_2 (u)
    integer, intent(in) :: u
    type(grid_t) :: grid
    write (u, "(A)")  "* Test output: grids_2"
    write (u, "(A)")  "*   Purpose: Saving and Loading"
    write (u, "(A)")

    call grid%init ([3])
    call grid%set_values ([one, two, three])
    call grid%save_to_file ('grids_2_test')
    call grid%final ()

    call assert (u, verify_points_for_grid('grids_2_test', [3]), &
         "verify_points_for_grid")
    call grid%load_from_file ('grids_2_test')
    call grid%write (u)
    call assert (u, nearly_equal (grid%get_value([1]), one),   "grid%get_value(1) == 1")
    call assert (u, nearly_equal (grid%get_value([2]), two),   "grid%get_value(2) == 2")
    call assert (u, nearly_equal (grid%get_value([3]), three), "grid%get_value(3) == 3")
    call grid%final ()

    call grid%init ([3,3])
    call grid%set_values ([one, two, three, four, zero, zero, zero, zero, zero])
    call grid%save_to_file ('grids_2_test')
    call grid%final ()

    call assert (u, verify_points_for_grid('grids_2_test', [3,3]), &
         "verify_points_for_grid")
    call grid%load_from_file ('grids_2_test')
    call grid%write (u)
    call assert (u, nearly_equal (grid%get_value([1,1]), one),   "grid%get_value(1,1) == 1")
    call assert (u, nearly_equal (grid%get_value([2,1]), two),   "grid%get_value(2,1) == 2")
    call assert (u, nearly_equal (grid%get_value([3,1]), three), "grid%get_value(3,1) == 3")
    call assert (u, nearly_equal (grid%get_value([1,2]), four),  "grid%get_value(1,2) == 4")
    call delete_file ('grids_2_test')

    call grid%load_from_file ('grids_2_test')
    call assert (u, .not. verify_points_for_grid('grids_2_test', [3,3]), &
         "verify_points_for_grid")
    call grid%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: grids_2"
  end subroutine grids_2

  subroutine grids_3 (u)
    integer, intent(in) :: u
    type(grid_t) :: grid
    integer, dimension(2) :: fail
    write (u, "(A)")  "* Test output: grids_3"
    write (u, "(A)")  "*   Purpose: Get Segments"
    write (u, "(A)")

    call grid%init ([3])
    call assert (u, all(grid%get_segment([0.00_default]) == [1]), &
                   "all(grid%get_segment([0.00_default]) == [1])")
    call assert (u, all(grid%get_segment([0.32_default]) == [1]), &
                   "all(grid%get_segment([0.32_default]) == [1])")
    call assert (u, all(grid%get_segment([0.52_default]) == [2]), &
                   "all(grid%get_segment([0.52_default]) == [2])")
    call assert (u, all(grid%get_segment([1.00_default]) == [3]), &
                   "all(grid%get_segment([1.00_default]) == [3])")
    call grid%final ()

    call grid%init ([3,3])
    call assert (u, all(grid%get_segment([0.00_default,0.00_default]) == [1,1]), &
                   "all(grid%get_segment([0.00_default,0.00_default]) == [1,1])")
    call assert (u, all(grid%get_segment([0.32_default,0.32_default]) == [1,1]), &
                   "all(grid%get_segment([0.32_default,0.32_default]) == [1,1])")
    call assert (u, all(grid%get_segment([0.52_default,0.52_default]) == [2,2]), &
                   "all(grid%get_segment([0.52_default,0.52_default]) == [2,2])")
    call assert (u, all(grid%get_segment([1.00_default,1.00_default]) == [3,3]), &
                   "all(grid%get_segment([1.00_default,1.00_default]) == [3,3])")
    write (u, "(A)")  "* A double error is expected"
    fail = grid%get_segment([1.10_default,1.10_default], u)
    call grid%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: grids_3"
  end subroutine grids_3

  subroutine grids_4 (u)
    integer, intent(in) :: u
    type(grid_t) :: grid
    write (u, "(A)")  "* Test output: grids_4"
    write (u, "(A)")  "*   Purpose: Update Maxima"
    write (u, "(A)")

    call grid%init ([4,4])
    call grid%update_maxima ([0.1_default, 0.0_default], 0.3_default)
    call grid%update_maxima ([0.9_default, 0.95_default], 1.7_default)
    call grid%write (u)
    call assert_equal (u, grid%get_value([1,1]), 0.3_default, &
               "grid%get_value([1,1]")
    call assert_equal (u, grid%get_value([2,2]), 0.0_default, &
               "grid%get_value([2,2]")
    call assert_equal (u, grid%get_value([4,4]), 1.7_default, &
               "grid%get_value([4,4]")

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: grids_4"
  end subroutine grids_4

  subroutine grids_5 (u)
    integer, intent(in) :: u
    type(grid_t) :: grid
    real(default) :: first, second
    write (u, "(A)")  "* Test output: grids_5"
    write (u, "(A)")  "*   Purpose: Finding and checking"
    write (u, "(A)")

    call grid%init ([2,2,2])
    first = one / two - tiny_07
    second = two / two - tiny_07
    call grid%update_maxima ([0.1_default, 0.0_default, first], 0.3_default)
    call grid%update_maxima ([0.9_default, 0.95_default, second], 1.7_default)
    call grid%write (u)
    call assert (u, .not. grid%is_non_zero_everywhere (), &
               ".not. grid%is_non_zero_everywhere (")
    call assert_equal (u, grid%get_maximum_in_3d (1), 0.3_default, &
         "grid%get_maximum_in_3d (1)")
    call assert_equal (u, grid%get_maximum_in_3d (2), 1.7_default, &
         "grid%get_maximum_in_3d (2)")

    call grid%update_maxima ([0.9_default, 0.95_default, first], 1.8_default)
    call grid%update_maxima ([0.1_default, 0.95_default, first], 1.5_default)
    call grid%update_maxima ([0.9_default, 0.15_default, first], 1.5_default)
    call grid%update_maxima ([0.1_default, 0.0_default, second], 0.2_default)
    call grid%update_maxima ([0.1_default, 0.9_default, second], 0.2_default)
    call grid%update_maxima ([0.9_default, 0.0_default, second], 0.2_default)
    call grid%write (u)
    call assert (u, grid%is_non_zero_everywhere (), &
               "grid%is_non_zero_everywhere (")
    call assert_equal (u, grid%get_maximum_in_3d (1), 1.8_default, &
         "grid%get_maximum_in_3d (1)")
    call assert_equal (u, grid%get_maximum_in_3d (2), 1.7_default, &
         "grid%get_maximum_in_3d (2)")

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: grids_5"
  end subroutine grids_5


end module grids_uti
