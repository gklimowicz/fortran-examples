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

submodule (physics_defs) physics_defs_s

  implicit none

contains

  elemental module function component_status_of_string (string) result (i)
    integer :: i
    type(string_t), intent(in) :: string
    select case (char(string))
    case ("born")
       i = BORN
    case ("real")
       i = NLO_REAL
    case ("virtual")
       i = NLO_VIRTUAL
    case ("mismatch")
       i = NLO_MISMATCH
    case ("dglap")
       i = NLO_DGLAP
    case ("subtraction")
       i = NLO_SUBTRACTION
    case ("full")
       i = NLO_FULL
    case ("GKS")
       i = GKS
    case default
       i = COMPONENT_UNDEFINED
    end select
  end function component_status_of_string

  elemental module function component_status_to_string (i) result (string)
    type(string_t) :: string
    integer, intent(in) :: i
    select case (i)
    case (BORN)
       string = "born"
    case (NLO_REAL)
       string = "real"
    case (NLO_VIRTUAL)
       string = "virtual"
    case (NLO_MISMATCH)
       string = "mismatch"
    case (NLO_DGLAP)
       string = "dglap"
    case (NLO_SUBTRACTION)
       string = "subtraction"
    case (NLO_FULL)
       string = "full"
    case (GKS)
       string = "GKS"
    case default
       string = "undefined"
    end select
  end function component_status_to_string

  elemental module function is_nlo_component (comp) result (is_nlo)
    logical :: is_nlo
    integer, intent(in) :: comp
    select case (comp)
    case (BORN : GKS)
       is_nlo = .true.
    case default
       is_nlo = .false.
    end select
  end function is_nlo_component

  module function is_subtraction_component (emitter, nlo_type) result (is_subtraction)
    logical :: is_subtraction
    integer, intent(in) :: emitter, nlo_type
    is_subtraction = nlo_type == NLO_REAL .and. emitter < 0
  end function is_subtraction_component

  module function thr_leg (emitter) result (leg)
    integer :: leg
    integer, intent(in) :: emitter
    leg = emitter - THR_EMITTER_OFFSET
  end function thr_leg


end submodule physics_defs_s

