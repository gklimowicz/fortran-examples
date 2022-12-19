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

module process_counter

  implicit none
  private

  public :: process_counter_t

  integer, parameter, public :: STAT_UNDEFINED = 0
  integer, parameter, public :: STAT_INITIAL = 1
  integer, parameter, public :: STAT_ACTIVATED = 2
  integer, parameter, public :: STAT_BEAM_MOMENTA = 3
  integer, parameter, public :: STAT_FAILED_KINEMATICS = 4
  integer, parameter, public :: STAT_SEED_KINEMATICS = 5
  integer, parameter, public :: STAT_HARD_KINEMATICS = 6
  integer, parameter, public :: STAT_EFF_KINEMATICS = 7
  integer, parameter, public :: STAT_FAILED_CUTS = 8
  integer, parameter, public :: STAT_PASSED_CUTS = 9
  integer, parameter, public :: STAT_EVALUATED_TRACE = 10
  integer, parameter, public :: STAT_EVENT_COMPLETE = 11


  type :: process_counter_t
     integer :: total = 0
     integer :: failed_kinematics = 0
     integer :: failed_cuts = 0
     integer :: has_passed = 0
     integer :: evaluated = 0
     integer :: complete = 0
   contains
     procedure :: write => process_counter_write
     procedure :: reset => process_counter_reset
     procedure :: record => process_counter_record
  end type process_counter_t


  interface
    module subroutine process_counter_write (object, unit)
      class(process_counter_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine process_counter_write
    module subroutine process_counter_reset (counter)
      class(process_counter_t), intent(out) :: counter
    end subroutine process_counter_reset
    module subroutine process_counter_record (counter, status)
      class(process_counter_t), intent(inout) :: counter
      integer, intent(in) :: status
    end subroutine process_counter_record
  end interface

end module process_counter
