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

submodule (matching_base) matching_base_s

  use debug_master, only: debug_on
  use diagnostics

  implicit none

contains

  pure module subroutine matching_import_rng (matching, rng)
    class(matching_t), intent(inout) :: matching
    class(rng_t), allocatable, intent(inout) :: rng
    call move_alloc (from = rng, to = matching%rng)
  end subroutine matching_import_rng

  module subroutine matching_connect (matching, process_instance, model, shower)
    class(matching_t), intent(inout) :: matching
    type(process_instance_t), intent(in), target :: process_instance
    class(model_data_t), intent(in), target, optional :: model
    class(shower_base_t), intent(in), target, optional :: shower
    if (debug_on) call msg_debug (D_MATCHING, "matching_connect")
    matching%process_instance => process_instance
    if (present (model))  matching%model => model
    if (present (shower))  matching%shower => shower
  end subroutine matching_connect

  module subroutine matching_prepare_for_events (matching)
    class(matching_t), intent(inout), target :: matching
  end subroutine matching_prepare_for_events

  module subroutine matching_first_event (matching)
    class(matching_t), intent(inout), target :: matching
  end subroutine matching_first_event

  module subroutine matching_final (matching)
    class(matching_t), intent(in) :: matching
  end subroutine matching_final

  elemental module function matching_method_of_string (string) result (i)
    integer :: i
    type(string_t), intent(in) :: string
    select case (char (string))
    case ("MLM")
       i = MATCH_MLM
    case ("CKKW")
       i = MATCH_CKKW
    case ("POWHEG")
       i = MATCH_POWHEG
    case default
       i = MATCH_UNDEFINED
    end select
  end function matching_method_of_string

  elemental module function matching_method_to_string (i) result (string)
    type(string_t) :: string
    integer, intent(in) :: i
    select case (i)
    case (MATCH_MLM)
       string = "MLM"
    case (MATCH_CKKW)
       string = "CKKW"
    case (MATCH_POWHEG)
       string = "POWHEG"
    case default
       string = "UNDEFINED"
    end select
  end function matching_method_to_string


end submodule matching_base_s

