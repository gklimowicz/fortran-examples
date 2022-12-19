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

submodule (c_particles) c_particles_s

  use io_units
  use format_defs, only: FMT_14, FMT_19

  implicit none

contains

  module subroutine c_prt_write (prt, unit)
    type(c_prt_t), intent(in) :: prt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)", advance="no")  "prt("
    write (u, "(I0,':')", advance="no")  prt%type
    if (prt%polarized /= 0) then
       write (u, "(I0,'/',I0,'|')", advance="no")  prt%pdg, prt%h
    else
       write (u, "(I0,'|')", advance="no") prt%pdg
    end if
    write (u, "(" // FMT_14 // ",';'," // FMT_14 // ",','," // &
         FMT_14 // ",','," // FMT_14 // ")", advance="no") &
         prt%pe, prt%px, prt%py, prt%pz
    write (u, "('|'," // FMT_19 // ")", advance="no")  prt%p2
    write (u, "(A)")  ")"
  end subroutine c_prt_write


end submodule c_particles_s

