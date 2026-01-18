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

module rng_stream_uti

  use kinds, only: default
  use kinds, only: i16
  use rng_base

  use rng_stream

  implicit none
  private

  public :: rng_stream_1
  public :: rng_stream_2
  public :: rng_stream_3

contains

  subroutine rng_stream_1 (u)
    integer, intent(in) :: u
    class(rng_t), allocatable, target :: rng

    real(default) :: x
    real(default), dimension(2) :: x2

    write (u, "(A)")  "* Test output: rng_stream_1"
    write (u, "(A)")  "*   Purpose: initialize and call the RNGstream random-number &
         &generator"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize generator (default seed)"
    write (u, "(A)")

    allocate (rng_stream_t :: rng)
    call rng%init ()

    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Get random number"
    write (u, "(A)")

    call rng%generate (x)
    write (u, "(A,2(1x,F9.7))")  "x =", x

    write (u, "(A)")
    write (u, "(A)")  "* Get random number pair"
    write (u, "(A)")

    call rng%generate (x2)
    write (u, "(A,2(1x,F9.7))")  "x =", x2

    write (u, "(A)")
    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Jump to next substream and get random number"
    write (u, "(A)")

    select type (rng)
    type is (rng_stream_t)
       call rng%next_substream ()
    end select
    call rng%generate (x)
    write (u, "(A,2(1x,F9.7))")  "x =", x

    write (u, "(A)")
    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call rng%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: rng_stream_1"

  end subroutine rng_stream_1

  subroutine rng_stream_2 (u)
    integer, intent(in) :: u
    type(rng_stream_factory_t) :: rng_factory
    class(rng_t), allocatable :: rng
    real(default) :: x

    write (u, "(A)")  "* Test output: rng_stream_2"
    write (u, "(A)")  "*   Purpose: initialize and use a rng factory"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize factory"
    write (u, "(A)")

    call rng_factory%init ()
    call rng_factory%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Make a generator"
    write (u, "(A)")

    call rng_factory%make (rng)
    call rng%write (u)
    call rng%generate (x)
    write (u, *)
    write (u, "(1x,A,F7.5)")  "x = ", x
    call rng%final ()
    deallocate (rng)

    write (u, "(A)")
    write (u, "(A)")  "* Repeat"
    write (u, "(A)")

    call rng_factory%make (rng)
    call rng%write (u)
    call rng%generate (x)
    write (u, *)
    write (u, "(1x,A,F7.5)")  "x = ", x
    call rng%final ()
    deallocate (rng)

    write (u, *)
    call rng_factory%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize factory with different seed"
    write (u, "(A)")

    call rng_factory%init (1_i16)
    call rng_factory%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Make a generator"
    write (u, "(A)")

    call rng_factory%make (rng)
    call rng%write (u)
    call rng%generate (x)
    write (u, *)
    write (u, "(1x,A,F7.5)")  "x = ", x
    call rng%final ()
    deallocate (rng)

    write (u, *)
    call rng_factory%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: rng_stream_2"

  end subroutine rng_stream_2

  subroutine rng_stream_3 (u)
    integer, intent(in) :: u
    class(rng_t), allocatable, target :: rng

    integer :: i
    real(default) :: x
    real(default), dimension(2) :: x2

    write (u, "(A)")  "* Test output: rng_stream_3"
    write (u, "(A)")  "*   Purpose: initialize and advance the state of the random-number &
         &generator"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize generator (default seed)"
    write (u, "(A)")

    allocate (rng_stream_t :: rng)
    call rng%init ()

    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Advance state by 15389 by iteration"
    write (u, "(A)")

    do i = 1, 15389
       call rng%generate (x)
    end do

    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Get random number"
    write (u, "(A)")

    call rng%generate (x)
    write (u, "(A,2(1x,F9.7))")  "x =", x

    write (u, "(A)")
    write (u, "(A)")  "* Advance state by 15389 by procedure"
    write (u, "(A)")

    select type (rng)
    type is (rng_stream_t)
       call rng%reset_substream ()
       call rng%advance_state (15389)
    end select

    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Get random number"
    write (u, "(A)")

    call rng%generate (x)
    write (u, "(A,2(1x,F9.7))")  "x =", x


    write (u, "(A)")
    write (u, "(A)")  "* Get random number pair"
    write (u, "(A)")

    call rng%generate (x2)
    write (u, "(A,2(1x,F9.7))")  "x =", x2

    write (u, "(A)")
    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Jump to next substream and get random number"
    write (u, "(A)")

    select type (rng)
    type is (rng_stream_t)
       call rng%next_substream ()
    end select
    call rng%generate (x)
    write (u, "(A,2(1x,F9.7))")  "x =", x

    write (u, "(A)")
    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call rng%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: rng_stream_3"

  end subroutine rng_stream_3


end module rng_stream_uti
