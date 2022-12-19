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

module matching_base

  use iso_varying_string, string_t => varying_string
  use sm_qcd
  use model_data
  use particles
  use variables
  use shower_base
  use instances, only: process_instance_t
  use rng_base

  implicit none
  private

  public :: matching_t
  public :: MATCH_MLM, MATCH_CKKW, MATCH_POWHEG
  public :: matching_method

  integer, parameter :: MATCH_MLM = 1
  integer, parameter :: MATCH_CKKW = 2
  integer, parameter :: MATCH_POWHEG = 3
  integer, parameter :: MATCH_UNDEFINED = 17

  type, abstract :: matching_t
    logical :: is_hadron_collision = .false.
    type(qcd_t) :: qcd
    class(shower_base_t), pointer :: shower => null ()
    type(process_instance_t), pointer :: process_instance => null ()
    class(model_data_t), pointer :: model => null ()
    class(rng_t), allocatable :: rng
    type(string_t) :: process_name
   contains
     procedure (matching_init), deferred :: init
     procedure (matching_write), deferred :: write
     procedure :: import_rng => matching_import_rng
     procedure :: connect => matching_connect
     procedure :: base_connect => matching_connect
     procedure (matching_before_shower), deferred :: before_shower
     procedure (matching_after_shower), deferred :: after_shower
     procedure :: prepare_for_events => matching_prepare_for_events
     procedure :: first_event => matching_first_event
     procedure (matching_get_method), deferred :: get_method
     procedure :: final => matching_final
  end type matching_t


  abstract interface
     subroutine matching_init (matching, var_list, process_name)
       import
       class(matching_t), intent(out) :: matching
       type(var_list_t), intent(in) :: var_list
       type(string_t), intent(in) :: process_name
     end subroutine matching_init
  end interface

  abstract interface
     subroutine matching_write (matching, unit)
       import
       class(matching_t), intent(in) :: matching
       integer, intent(in), optional :: unit
     end subroutine matching_write
  end interface

  abstract interface
     subroutine matching_before_shower (matching, particle_set, vetoed)
       import
       class(matching_t), intent(inout) :: matching
       type(particle_set_t), intent(inout) :: particle_set
       logical, intent(out) :: vetoed
     end subroutine matching_before_shower
  end interface

  abstract interface
     subroutine matching_after_shower (matching, particle_set, vetoed)
       import
       class(matching_t), intent(inout) :: matching
       type(particle_set_t), intent(inout) :: particle_set
       logical, intent(out) :: vetoed
     end subroutine matching_after_shower
  end interface

  abstract interface
     function matching_get_method (matching) result (method)
       import
       type(string_t) :: method
       class(matching_t), intent(in) :: matching
     end function matching_get_method
  end interface

  interface matching_method
     module procedure matching_method_of_string
     module procedure matching_method_to_string
  end interface

  interface
    pure module subroutine matching_import_rng (matching, rng)
      class(matching_t), intent(inout) :: matching
      class(rng_t), allocatable, intent(inout) :: rng
    end subroutine matching_import_rng
    module subroutine matching_connect &
         (matching, process_instance, model, shower)
      class(matching_t), intent(inout) :: matching
      type(process_instance_t), intent(in), target :: process_instance
      class(model_data_t), intent(in), target, optional :: model
      class(shower_base_t), intent(in), target, optional :: shower
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
    end function matching_method_of_string
    elemental module function matching_method_to_string (i) result (string)
      type(string_t) :: string
      integer, intent(in) :: i
    end function matching_method_to_string
  end interface

end module matching_base
