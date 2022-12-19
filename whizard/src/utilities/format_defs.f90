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

module format_defs

  implicit none
  private

  character(*), parameter, public :: FMT_19 = "ES19.12"
  character(*), parameter, public :: FMT_18 = "ES18.11"
  character(*), parameter, public :: FMT_17 = "ES17.10"
  character(*), parameter, public :: FMT_16 = "ES16.9"
  character(*), parameter, public :: FMT_15 = "ES15.8"
  character(*), parameter, public :: FMT_14 = "ES14.7"
  character(*), parameter, public :: FMT_13 = "ES13.6"
  character(*), parameter, public :: FMT_12 = "ES12.5"
  character(*), parameter, public :: FMT_11 = "ES11.4"
  character(*), parameter, public :: FMT_10 = "ES10.3"

  character(*), parameter, public :: FMF_12 = "F12.9"

end module format_defs
