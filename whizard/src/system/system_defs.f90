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

module system_defs

  use, intrinsic :: iso_fortran_env, only: iostat_end, iostat_eor !NODEP!

  implicit none
  private

  integer, parameter, public :: VERSION_STRLEN = 255
  character(len=VERSION_STRLEN), parameter, public :: &
       & VERSION_STRING = "WHIZARD version 3.1.0 (Dec 14 2022)"

  integer, parameter, public :: BUFFER_SIZE = 1000

  integer, parameter, public :: EOF = iostat_end,  EOR = iostat_eor

  character, parameter, public :: BLANK = ' '
  character, parameter, public :: TAB = achar(9)
  character, parameter, public :: CR = achar(13)
  character, parameter, public :: LF = achar(10)
  character, parameter, public :: BACKSLASH = achar(92)

  character(*), parameter, public :: WHITESPACE_CHARS = BLANK// TAB // CR // LF
  character(*), parameter, public :: LCLETTERS = "abcdefghijklmnopqrstuvwxyz"
  character(*), parameter, public :: UCLETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  character(*), parameter, public :: DIGIT_CHARS = "0123456789"

  integer, parameter, public :: MAX_ERRORS = 10
  integer, parameter, public :: ENVVAR_LEN = 1000
  integer, parameter, public :: DLERROR_LEN = 160

end module system_defs
