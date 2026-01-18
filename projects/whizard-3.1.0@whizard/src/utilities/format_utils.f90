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

module format_utils

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: write_separator
  public :: write_indent
  public :: write_integer_array
  public :: quote_underscore
  public :: tex_format
  public :: mp_format
  public :: pac_fmt
  public :: refmt_tiny
  public :: write_compressed_integer_array

  interface
    module subroutine write_separator (u, mode)
      integer, intent(in) :: u
      integer, intent(in), optional :: mode
    end subroutine write_separator
    module subroutine write_indent (unit, indent)
      integer, intent(in) :: unit
      integer, intent(in), optional :: indent
    end subroutine write_indent
    module subroutine write_integer_array (array, unit, n_max, no_skip)
      integer, intent(in), dimension(:) :: array
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: n_max
      logical, intent(in), optional :: no_skip
    end subroutine write_integer_array
    module function quote_underscore (string) result (quoted)
      type(string_t) :: quoted
      type(string_t), intent(in) :: string
    end function quote_underscore
    module function tex_format (rval, n_digits) result (string)
      type(string_t) :: string
      real(default), intent(in) :: rval
      integer, intent(in) :: n_digits
    end function tex_format
    module function mp_format (rval) result (string)
      type(string_t) :: string
      real(default), intent(in) :: rval
    end function mp_format
    module subroutine pac_fmt (fmt, fmt_orig, fmt_pac, pacify)
      character(*), intent(in) :: fmt_orig, fmt_pac
      character(*), intent(out) :: fmt
      logical, intent(in), optional :: pacify
    end subroutine pac_fmt
    elemental module function refmt_tiny (val) result (trunc_val)
      real(default), intent(in) :: val
      real(default) :: trunc_val
    end function refmt_tiny
    module subroutine write_compressed_integer_array (chars, array)
      character(len=*), intent(out) :: chars
      integer, intent(in), allocatable, dimension(:) :: array
    end subroutine write_compressed_integer_array
  end interface

end module format_utils
