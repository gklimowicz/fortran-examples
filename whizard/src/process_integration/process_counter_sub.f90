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

submodule (process_counter) process_counter_s

  use io_units

  implicit none

contains

  module subroutine process_counter_write (object, unit)
    class(process_counter_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (object%total > 0) then
       write (u, "(1x,A)")  "Call statistics (current run):"
       write (u, "(3x,A,I0)")  "total       = ", object%total
       write (u, "(3x,A,I0)")  "failed kin. = ", object%failed_kinematics
       write (u, "(3x,A,I0)")  "failed cuts = ", object%failed_cuts
       write (u, "(3x,A,I0)")  "passed cuts = ", object%has_passed
       write (u, "(3x,A,I0)")  "evaluated   = ", object%evaluated
    else
       write (u, "(1x,A)")  "Call statistics (current run): [no calls]"
    end if
  end subroutine process_counter_write

  module subroutine process_counter_reset (counter)
    class(process_counter_t), intent(out) :: counter
    counter%total = 0
    counter%failed_kinematics = 0
    counter%failed_cuts = 0
    counter%has_passed = 0
    counter%evaluated = 0
    counter%complete = 0
  end subroutine process_counter_reset

  module subroutine process_counter_record (counter, status)
    class(process_counter_t), intent(inout) :: counter
    integer, intent(in) :: status
    if (status <= STAT_FAILED_KINEMATICS) then
       counter%failed_kinematics = counter%failed_kinematics + 1
    else if (status <= STAT_FAILED_CUTS) then
       counter%failed_cuts = counter%failed_cuts + 1
    else if (status <= STAT_PASSED_CUTS) then
       counter%has_passed = counter%has_passed + 1
    else
       counter%evaluated = counter%evaluated + 1
    end if
    counter%total = counter%total + 1
  end subroutine process_counter_record


end submodule process_counter_s

