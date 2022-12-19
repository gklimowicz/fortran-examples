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

submodule (hashes) hashes_s

  use bytes

  implicit none

contains

  module function hash (key) result (hashval)
    integer(i32) :: hashval
    integer(i8), dimension(:), intent(in) :: key
    type(word32_t) :: w
    integer :: i
    w = 0_i32
    do i = 1, size (key)
       w = w + key(i)
       w = w + ishft (w, 10)
       w = ieor (w, ishft (w, -6))
    end do
    w = w + ishft (w, 3)
    w = ieor (w, ishft (w, -11))
    w = w + ishft (w, 15)
    hashval = w
  end function hash


end submodule hashes_s

