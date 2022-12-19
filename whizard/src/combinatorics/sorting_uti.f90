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

module sorting_uti

  use kinds, only: default

  use sorting

  implicit none
  private

  public :: sorting_1

contains

  subroutine sorting_1 (u)
    integer, intent(in) :: u
    integer, parameter :: NMAX = 10
    real(default), dimension(NMAX) :: rval
    integer, dimension(NMAX) :: ival
    real, dimension(NMAX,NMAX) :: harvest_r
    integer, dimension(NMAX,NMAX) :: harvest_i
    integer, dimension(NMAX,NMAX) :: harvest_a
    integer :: i, j
    harvest_r(:, 1) = [0.9976, 0., 0., 0., 0., 0., 0., 0., 0., 0.]
    harvest_r(:, 2) = [0.5668, 0.9659, 0., 0., 0., 0., 0., 0., 0., 0.]
    harvest_r(:, 3) = [0.7479, 0.3674, 0.4806, 0., 0., 0., 0., 0., 0., &
         0.]
    harvest_r(:, 4) = [0.0738, 0.0054, 0.3471, 0.3422, 0., 0., 0., 0., &
         0., 0.]
    harvest_r(:, 5) = [0.2180, 0.1332, 0.9005, 0.3868, 0.4455, 0., 0., &
         0., 0., 0.]
    harvest_r(:, 6) = [0.6619, 0.0161, 0.6509, 0.6464, 0.3230, &
         0.8557, 0., 0., 0., 0.]
    harvest_r(:, 7) = [0.4013, 0.2069, 0.9685, 0.5984, 0.6730, &
         0.4569, 0.3300, 0., 0., 0.]
    harvest_r(:, 8) = [0.1004, 0.7555, 0.6057, 0.7190, 0.8973, &
         0.6582, 0.1507, 0.6123, 0., 0.]
    harvest_r(:, 9) = [0.9787, 0.9991, 0.2568, 0.5509, 0.6590, &
         0.5540, 0.9778, 0.9019, 0.6579, 0.]
    harvest_r(:,10) = [0.7289, 0.4025, 0.9286, 0.1478, 0.6745, &
         0.7696, 0.3393, 0.1158, 0.6144, 0.8206]

    harvest_i(:, 1) = [18, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    harvest_i(:, 2) = [14, 9, 0, 0, 0, 0, 0, 0, 0, 0]
    harvest_i(:, 3) = [ 7, 8,11, 0, 0, 0, 0, 0, 0, 0]
    harvest_i(:, 4) = [19,19,14,19, 0, 0, 0, 0, 0, 0]
    harvest_i(:, 5) = [ 1,14,15,18,14, 0, 0, 0, 0, 0]
    harvest_i(:, 6) = [16,11, 1, 9,11, 2, 0, 0, 0, 0]
    harvest_i(:, 7) = [11,10,17, 6,13,13,10, 0, 0, 0]
    harvest_i(:, 8) = [ 5, 1, 2,10, 7, 0,15,12, 0, 0]
    harvest_i(:, 9) = [15,19, 2, 6,11, 0, 2, 4, 2, 0]
    harvest_i(:,10) = [ 1, 4, 8, 4,11, 0, 8, 7,19,13]

    harvest_a(:, 1) = [-6,  0,  0,  0,  0,  0,  0,  0,  0,  0]
    harvest_a(:, 2) = [-8, -9,  0,  0,  0,  0,  0,  0,  0,  0]
    harvest_a(:, 3) = [ 4, -3,  3,  0,  0,  0,  0,  0,  0,  0]
    harvest_a(:, 4) = [-6,  6,  2, -2,  0,  0,  0,  0,  0,  0]
    harvest_a(:, 5) = [ 1, -2,  0, -6,  8,  0,  0,  0,  0,  0]
    harvest_a(:, 6) = [-2, -1, -8, -5,  8, -5,  0,  0,  0,  0]
    harvest_a(:, 7) = [-9,  0, -6,  2,  5,  3,  2,  0,  0,  0]
    harvest_a(:, 8) = [-5, -7,  6,  7, -3,  0, -7,  4,  0,  0]
    harvest_a(:, 9) = [ 5,  0, -1, -7,  5,  2,  7, -3,  3,  0]
    harvest_a(:,10) = [-9,  2, -6,  3, -9,  5,  5,  7,  5, -9]


    write (u, "(A)")  "* Test output: Sorting"
    write (u, "(A)")  "*   Purpose: test sorting routines"
    write (u, "(A)")

    write (u, "(A)")  "* Sorting real values:"

    do i = 1, NMAX
       write (u, "(A)")
       rval(:i) = harvest_r(:i,i)
       write (u, "(10(1x,F7.4))") rval(:i)
       rval(:i) = sort (rval(:i))
       write (u, "(10(1x,F7.4))") rval(:i)
       do j = i, 2, -1
          if (rval(j)-rval(j-1) < 0) &
             write (u, "(A)") "*** Sorting failure. ***"
       end do
    end do

    write (u, "(A)")
    write (u, "(A)") "* Sorting integer values:"

    do i = 1, NMAX
       write (u, "(A)")
       ival(:i) = harvest_i(:i,i)
       write (u, "(10(1x,I2))") ival(:i)
       ival(:i) = sort (ival(:i))
       write (u, "(10(1x,I2))") ival(:i)
       do j = i, 2, -1
          if (ival(j)-ival(j-1) < 0) &
             write (u, "(A)")  "*** Sorting failure. ***"
       end do
    end do

    write (u, "(A)")
    write (u, "(A)") "* Sorting integer values by absolute value:"

    do i = 1, NMAX
       write (u, "(A)")
       ival(:i) = harvest_a(:i,i)
       write (u, "(10(1x,I2))") ival(:i)
       ival(:i) = sort_abs (ival(:i))
       write (u, "(10(1x,I2))") ival(:i)
       do j = i, 2, -1
          if (abs(ival(j))-abs(ival(j-1)) < 0 .or. &
               (abs(ival(j))==abs(ival(j-1))) .and. ival(j)>ival(j-1)) &
             write (u, "(A)")  "*** Sorting failure. ***"
       end do
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sorting_1"

  end subroutine sorting_1


end module sorting_uti
