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

module selectors_uti

  use kinds, only: default
  use rng_base

  use selectors

  use rng_base_ut, only: rng_test_t

  implicit none
  private

  public :: selectors_1
  public :: selectors_2

contains

  subroutine selectors_1 (u)
    integer, intent(in) :: u
    type(selector_t) :: selector
    class(rng_t), allocatable, target :: rng
    integer :: i, n

    write (u, "(A)")  "* Test output: selectors_1"
    write (u, "(A)")  "*   Purpose: initialize a selector and test it"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize selector"
    write (u, "(A)")

    call selector%init &
         ([2._default, 3.5_default, 0._default, &
         2._default, 0.5_default, 2._default])
    call selector%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Select numbers using predictable test generator"
    write (u, "(A)")

    allocate (rng_test_t :: rng)
    call rng%init (1)

    do i = 1, 5
       call selector%generate (rng, n)
       write (u, "(1x,I0)")  n
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Select numbers using real input number"
    write (u, "(A)")

    write (u, "(1x,A,I0)")  "select(0.00) = ", selector%select (0._default)
    write (u, "(1x,A,I0)")  "select(0.77) = ", selector%select (0.77_default)
    write (u, "(1x,A,I0)")  "select(1.00) = ", selector%select (1._default)

    write (u, "(A)")
    write (u, "(A)")  "* Get weight"
    write (u, "(A)")

    write (u, "(1x,A,ES19.12)")  "weight(2) =", selector%get_weight(2)
    write (u, "(1x,A,ES19.12)")  "weight(3) =", selector%get_weight(3)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call rng%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: selectors_1"

  end subroutine selectors_1

  subroutine selectors_2 (u)
    integer, intent(in) :: u
    type(selector_t) :: selector
    class(rng_t), allocatable, target :: rng
    integer :: i, n

    write (u, "(A)")  "* Test output: selectors_2"
    write (u, "(A)")  "*   Purpose: initialize and use a selector &
         &with offset index"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize selector"
    write (u, "(A)")

    call selector%init &
         ([2._default, 3.5_default, 0._default, &
         2._default, 0.5_default, 2._default], &
         offset = -1)
    call selector%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Select numbers using predictable test generator"
    write (u, "(A)")

    allocate (rng_test_t :: rng)
    call rng%init (1)

    do i = 1, 5
       call selector%generate (rng, n)
       write (u, "(1x,I0)")  n
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Select numbers using real input number"
    write (u, "(A)")

    write (u, "(1x,A,I0)")  "select(0.00) = ", selector%select (0._default)
    write (u, "(1x,A,I0)")  "select(0.77) = ", selector%select (0.77_default)
    write (u, "(1x,A,I0)")  "select(1.00) = ", selector%select (1._default)

    write (u, "(A)")
    write (u, "(A)")  "* Get weight"
    write (u, "(A)")

    write (u, "(1x,A,ES19.12)")  "weight(1) =", selector%get_weight(1)
    write (u, "(1x,A,ES19.12)")  "weight(2) =", selector%get_weight(2)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call rng%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: selectors_2"

  end subroutine selectors_2


end module selectors_uti
