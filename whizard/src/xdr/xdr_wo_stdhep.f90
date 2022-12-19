!$Id: xdr_wo_stdhep.f90 6133 2014-09-17 14:42:33Z kilian $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     with contributions from
!     Fabian Bach <fabian.bach@t-online.de>
!     Christian Speckner <cnspeckn@googlemail.com>
!     Christian Weiss <christian.weiss@desy.de>
!     and Hans-Werner Boschmann, Felix Braam,
!     Sebastian Schmidt, Daniel Wiesler 
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

module xdr_wo_stdhep

  use, intrinsic :: iso_c_binding
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: write_stdhep_event

  interface
     subroutine print_stdhep_to_file (infile, outfile, maxevt) bind(C)
       import
       character(c_char), dimension(*), intent(in) :: infile, outfile       
       integer(c_int), value :: maxevt
     end subroutine print_stdhep_to_file
  end interface

contains

  subroutine write_stdhep_event (infile, outfile, maxevt)
    type(string_t), intent(in) :: infile, outfile
    integer, intent(in) :: maxevt    
    call print_stdhep_to_file (char (infile) // c_null_char, &
          char (outfile) // c_null_char, int (maxevt, c_int))
  end subroutine write_stdhep_event
  
end module xdr_wo_stdhep

