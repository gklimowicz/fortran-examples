! WHIZARD <<Version>> <<Date>>
! 
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     with contributions from
!     Christian Speckner <cnspeckn@googlemail.com>
!     and  Fabian Bach, Felix Braam, Sebastian Schmidt, Daniel Wiesler
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
module io_units

  use, intrinsic :: iso_fortran_env, only: output_unit

  implicit none
  private

  ! ----------------------------------------------------------------------
  ! Public signature

  ! This function returns a free I/O unit.
  ! (Obsolete in F2008 where the NEWUNIT specifier has been introduced)
  ! Usage: u = free_unit ()

  public :: free_unit

  ! This function returns the standard output unit if no argument is given,
  ! otherwise the argument.
  ! Usage: u = given_output_unit (unit)

  public :: given_output_unit

  ! End of public signature
  ! ----------------------------------------------------------------------

  integer, parameter, private :: MIN_UNIT = 11, MAX_UNIT = 99

contains

  function free_unit () result (unit)
    integer :: unit
    logical :: exists, is_open
    integer :: i, status
    do i = MIN_UNIT, MAX_UNIT
       inquire (unit=i, exist=exists, opened=is_open, iostat=status)
       if (status == 0) then
          if (exists .and. .not. is_open) then
            unit = i; return
          end if
       end if
    end do
    unit = -1
  end function free_unit

  function given_output_unit (unit) result (u)
    integer, intent(in), optional :: unit
    integer :: u
    if (present (unit)) then
       u = unit
    else
       u = output_unit
    end if
  end function given_output_unit

end module io_units

