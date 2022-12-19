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

module rng_tao_uti

  use kinds, only: default
  use kinds, only: i16
  use rng_base

  use rng_tao

  implicit none
  private

  public :: rng_tao_1
  public :: rng_tao_2

contains

  subroutine rng_tao_1 (u)
    integer, intent(in) :: u
    class(rng_t), allocatable, target :: rng

    real(default) :: x
    real(default), dimension(2) :: x2

    write (u, "(A)")  "* Test output: rng_tao_1"
    write (u, "(A)")  "*   Purpose: initialize and call the TAO random-number &
         &generator"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize generator (default seed)"
    write (u, "(A)")

    allocate (rng_tao_t :: rng)
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
    write (u, "(A)")  "* Cleanup"

    call rng%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: rng_tao_1"

  end subroutine rng_tao_1

  subroutine rng_tao_2 (u)
    integer, intent(in) :: u
    type(rng_tao_factory_t) :: rng_factory
    class(rng_t), allocatable :: rng
    real(default) :: x

    write (u, "(A)")  "* Test output: rng_tao_2"
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
    write (u, "(A)")  "* Test output end: rng_tao_2"

  end subroutine rng_tao_2


end module rng_tao_uti
