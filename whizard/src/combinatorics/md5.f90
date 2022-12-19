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

module md5

  use bytes

  implicit none
  private

  public :: md5sum

  type :: block_t
     private
     type(word32_t), dimension(0:15) :: w
     type(block_t), pointer :: next => null ()
     integer :: fill = 0
  end type block_t

  type :: message_t
     private
     type(block_t), pointer :: first => null ()
     type(block_t), pointer :: last => null ()
     integer :: n_blocks = 0
  end type message_t


  interface block_write
     module procedure block_write_unit
  end interface
  interface message_write
     module procedure message_write_unit
  end interface
  interface md5sum
     module procedure md5sum_from_string
     module procedure md5sum_from_unit
  end interface

  interface
    module subroutine block_write_unit (b, unit, bytes, decimal)
      type(block_t), intent(in) :: b
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: bytes, decimal
    end subroutine block_write_unit
    module subroutine message_write_unit (m, unit, bytes, decimal)
      type(message_t), intent(in) :: m
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: bytes, decimal
    end subroutine message_write_unit
    module function md5sum_from_string (s) result (digest)
      character(len=*), intent(in) :: s
      character(len=32) :: digest
    end function md5sum_from_string
    module function md5sum_from_unit (u) result (digest)
      integer, intent(in) :: u
      character(len=32) :: digest
    end function md5sum_from_unit
  end interface

end module md5
