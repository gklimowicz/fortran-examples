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

submodule (event_base) event_base_s

  use io_units
  use string_utils, only: lower_case
  use diagnostics

  implicit none

contains

  module subroutine generic_event_init (event, n_alt)
    class(generic_event_t), intent(out) :: event
    integer, intent(in) :: n_alt
    event%n_alt = n_alt
    allocate (event%sqme_alt (n_alt))
    allocate (event%weight_alt (n_alt))
  end subroutine generic_event_init

  module function generic_event_has_valid_particle_set (event) result (flag)
    class(generic_event_t), intent(in) :: event
    logical :: flag
    flag = event%particle_set_is_valid
  end function generic_event_has_valid_particle_set

  module subroutine generic_event_accept_particle_set (event)
    class(generic_event_t), intent(inout) :: event
    event%particle_set_is_valid = .true.
  end subroutine generic_event_accept_particle_set

  module subroutine generic_event_discard_particle_set (event)
    class(generic_event_t), intent(inout) :: event
    event%particle_set_is_valid = .false.
  end subroutine generic_event_discard_particle_set

  module function generic_event_get_particle_set_ptr (event) result (ptr)
    class(generic_event_t), intent(in) :: event
    type(particle_set_t), pointer :: ptr
    ptr => event%particle_set
  end function generic_event_get_particle_set_ptr

  module subroutine generic_event_link_particle_set (event, particle_set)
    class(generic_event_t), intent(inout) :: event
    type(particle_set_t), intent(in), target :: particle_set
    event%particle_set => particle_set
    call event%accept_particle_set ()
  end subroutine generic_event_link_particle_set

  module function generic_event_sqme_prc_is_known (event) result (flag)
    class(generic_event_t), intent(in) :: event
    logical :: flag
    flag = event%sqme_prc_known
  end function generic_event_sqme_prc_is_known

  module function generic_event_sqme_ref_is_known (event) result (flag)
    class(generic_event_t), intent(in) :: event
    logical :: flag
    flag = event%sqme_ref_known
  end function generic_event_sqme_ref_is_known

  module function generic_event_sqme_alt_is_known (event) result (flag)
    class(generic_event_t), intent(in) :: event
    logical :: flag
    flag = event%sqme_alt_known
  end function generic_event_sqme_alt_is_known

  module function generic_event_weight_prc_is_known (event) result (flag)
    class(generic_event_t), intent(in) :: event
    logical :: flag
    flag = event%weight_prc_known
  end function generic_event_weight_prc_is_known

  module function generic_event_weight_ref_is_known (event) result (flag)
    class(generic_event_t), intent(in) :: event
    logical :: flag
    flag = event%weight_ref_known
  end function generic_event_weight_ref_is_known

  module function generic_event_weight_alt_is_known (event) result (flag)
    class(generic_event_t), intent(in) :: event
    logical :: flag
    flag = event%weight_alt_known
  end function generic_event_weight_alt_is_known

  module function generic_event_excess_prc_is_known (event) result (flag)
    class(generic_event_t), intent(in) :: event
    logical :: flag
    flag = event%excess_prc_known
  end function generic_event_excess_prc_is_known

  module function generic_event_get_n_alt (event) result (n)
    class(generic_event_t), intent(in) :: event
    integer :: n
    n = event%n_alt
  end function generic_event_get_n_alt

  module function generic_event_get_sqme_prc (event) result (sqme)
    class(generic_event_t), intent(in) :: event
    real(default) :: sqme
    if (event%sqme_prc_known) then
       sqme = event%sqme_prc
    else
       sqme = 0
    end if
  end function generic_event_get_sqme_prc

  module function generic_event_get_sqme_ref (event) result (sqme)
    class(generic_event_t), intent(in) :: event
    real(default) :: sqme
    if (event%sqme_ref_known) then
       sqme = event%sqme_ref
    else
       sqme = 0
    end if
  end function generic_event_get_sqme_ref

  module function generic_event_get_sqme_alt_0 (event, i) result (sqme)
    class(generic_event_t), intent(in) :: event
    integer, intent(in) :: i
    real(default) :: sqme
    if (event%sqme_alt_known) then
       sqme = event%sqme_alt(i)
    else
       sqme = 0
    end if
  end function generic_event_get_sqme_alt_0

  module function generic_event_get_sqme_alt_1 (event) result (sqme)
    class(generic_event_t), intent(in) :: event
    real(default), dimension(event%n_alt) :: sqme
    sqme = event%sqme_alt
  end function generic_event_get_sqme_alt_1

  module function generic_event_get_weight_prc (event) result (weight)
    class(generic_event_t), intent(in) :: event
    real(default) :: weight
    if (event%weight_prc_known) then
       weight = event%weight_prc
    else
       weight = 0
    end if
  end function generic_event_get_weight_prc

  module function generic_event_get_weight_ref (event) result (weight)
    class(generic_event_t), intent(in) :: event
    real(default) :: weight
    if (event%weight_ref_known) then
       weight = event%weight_ref
    else
       weight = 0
    end if
  end function generic_event_get_weight_ref

  module function generic_event_get_weight_alt_0 (event, i) result (weight)
    class(generic_event_t), intent(in) :: event
    integer, intent(in) :: i
    real(default) :: weight
    if (event%weight_alt_known) then
       weight = event%weight_alt(i)
    else
       weight = 0
    end if
  end function generic_event_get_weight_alt_0

  module function generic_event_get_weight_alt_1 (event) result (weight)
    class(generic_event_t), intent(in) :: event
    real(default), dimension(event%n_alt) :: weight
    weight = event%weight_alt
  end function generic_event_get_weight_alt_1

  module function generic_event_get_excess_prc (event) result (excess)
    class(generic_event_t), intent(in) :: event
    real(default) :: excess
    if (event%excess_prc_known) then
       excess = event%excess_prc
    else
       excess = 0
    end if
  end function generic_event_get_excess_prc

  module function generic_event_get_n_dropped (event) result (n_dropped)
    class(generic_event_t), intent(in) :: event
    integer :: n_dropped
    if (event%n_dropped_known) then
       n_dropped = event%n_dropped
    else
       n_dropped = 0
    end if
  end function generic_event_get_n_dropped

  module subroutine generic_event_set_sqme_prc (event, sqme)
    class(generic_event_t), intent(inout) :: event
    real(default), intent(in) :: sqme
    event%sqme_prc = sqme
    event%sqme_prc_known = .true.
  end subroutine generic_event_set_sqme_prc

  module subroutine generic_event_set_sqme_ref (event, sqme)
    class(generic_event_t), intent(inout) :: event
    real(default), intent(in) :: sqme
    event%sqme_ref = sqme
    event%sqme_ref_known = .true.
  end subroutine generic_event_set_sqme_ref

  module subroutine generic_event_set_sqme_alt (event, sqme)
    class(generic_event_t), intent(inout) :: event
    real(default), dimension(:), intent(in) :: sqme
    event%sqme_alt = sqme
    event%sqme_alt_known = .true.
  end subroutine generic_event_set_sqme_alt

  module subroutine generic_event_set_weight_prc (event, weight)
    class(generic_event_t), intent(inout) :: event
    real(default), intent(in) :: weight
    event%weight_prc = weight
    event%weight_prc_known = .true.
  end subroutine generic_event_set_weight_prc

  module subroutine generic_event_set_weight_ref (event, weight)
    class(generic_event_t), intent(inout) :: event
    real(default), intent(in) :: weight
    event%weight_ref = weight
    event%weight_ref_known = .true.
  end subroutine generic_event_set_weight_ref

  module subroutine generic_event_set_weight_alt (event, weight)
    class(generic_event_t), intent(inout) :: event
    real(default), dimension(:), intent(in) :: weight
    event%weight_alt = weight
    event%weight_alt_known = .true.
  end subroutine generic_event_set_weight_alt

  module subroutine generic_event_set_excess_prc (event, excess)
    class(generic_event_t), intent(inout) :: event
    real(default), intent(in) :: excess
    event%excess_prc = excess
    event%excess_prc_known = .true.
  end subroutine generic_event_set_excess_prc

  module subroutine generic_event_set_n_dropped (event, n_dropped)
    class(generic_event_t), intent(inout) :: event
    integer, intent(in) :: n_dropped
    event%n_dropped = n_dropped
    event%n_dropped_known = .true.
  end subroutine generic_event_set_n_dropped

  module subroutine generic_event_set (event, &
       weight_ref, weight_prc, weight_alt, &
       excess_prc, n_dropped, &
       sqme_ref, sqme_prc, sqme_alt)
    class(generic_event_t), intent(inout) :: event
    real(default), intent(in), optional :: weight_ref, weight_prc
    real(default), intent(in), optional :: sqme_ref, sqme_prc
    real(default), dimension(:), intent(in), optional :: sqme_alt, weight_alt
    real(default), intent(in), optional :: excess_prc
    integer, intent(in), optional :: n_dropped
    if (present (sqme_prc)) then
       call event%set_sqme_prc (sqme_prc)
    end if
    if (present (sqme_ref)) then
       call event%set_sqme_ref (sqme_ref)
    end if
    if (present (sqme_alt)) then
       call event%set_sqme_alt (sqme_alt)
    end if
    if (present (weight_prc)) then
       call event%set_weight_prc (weight_prc)
    end if
    if (present (weight_ref)) then
       call event%set_weight_ref (weight_ref)
    end if
    if (present (weight_alt)) then
       call event%set_weight_alt (weight_alt)
    end if
    if (present (excess_prc)) then
       call event%set_excess_prc (excess_prc)
    end if
    if (present (n_dropped)) then
       call event%set_n_dropped (n_dropped)
    end if
  end subroutine generic_event_set

  module subroutine generic_event_reset_contents (event)
    class(generic_event_t), intent(inout) :: event
    call event%discard_particle_set ()
    event%sqme_ref_known = .false.
    event%sqme_prc_known = .false.
    event%sqme_alt_known = .false.
    event%weight_ref_known = .false.
    event%weight_prc_known = .false.
    event%weight_alt_known = .false.
    event%excess_prc_known = .false.
  end subroutine generic_event_reset_contents

  module subroutine generic_event_pacify_particle_set (event)
    class(generic_event_t), intent(inout) :: event
    if (event%has_valid_particle_set ())  call pacify (event%particle_set)
  end subroutine generic_event_pacify_particle_set

  module function event_normalization_mode (string, unweighted) result (mode)
    integer :: mode
    type(string_t), intent(in) :: string
    logical, intent(in) :: unweighted
    select case (lower_case (char (string)))
    case ("auto")
       if (unweighted) then
          mode = NORM_UNIT
       else
          mode = NORM_SIGMA
       end if
    case ("1")
       mode = NORM_UNIT
    case ("1/n")
       mode = NORM_N_EVT
    case ("sigma")
       mode = NORM_SIGMA
    case ("sigma/n")
       mode = NORM_S_N
    case default
       call msg_fatal ("Event normalization: unknown value '" &
            // char (string) // "'")
    end select
  end function event_normalization_mode

  module function event_normalization_string (norm_mode) result (string)
    integer, intent(in) :: norm_mode
    type(string_t) :: string
    select case (norm_mode)
    case (NORM_UNDEFINED); string = "[undefined]"
    case (NORM_UNIT);      string = "'1'"
    case (NORM_N_EVT);     string = "'1/n'"
    case (NORM_SIGMA);     string = "'sigma'"
    case (NORM_S_N);       string = "'sigma/n'"
    case default;          string = "???"
    end select
  end function event_normalization_string

  module subroutine event_normalization_update &
       (weight, sigma, n, mode_new, mode_old)
    real(default), intent(inout) :: weight
    real(default), intent(in) :: sigma
    integer, intent(in) :: n
    integer, intent(in) :: mode_new, mode_old
    if (mode_new /= mode_old) then
       if (sigma > 0 .and. n > 0) then
          weight = weight / factor (mode_old) * factor (mode_new)
       else
          call msg_fatal ("Event normalization update: null sample")
       end if
    end if
  contains
    function factor (mode)
      real(default) :: factor
      integer, intent(in) :: mode
      select case (mode)
      case (NORM_UNIT);   factor = 1._default
      case (NORM_N_EVT);  factor = 1._default / n
      case (NORM_SIGMA);  factor = sigma
      case (NORM_S_N);    factor = sigma / n
      case default
         call msg_fatal ("Event normalization update: undefined mode")
      end select
    end function factor
  end subroutine event_normalization_update

  module subroutine event_callback_nop_write (event_callback, unit)
    class(event_callback_nop_t), intent(in) :: event_callback
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "NOP"
  end subroutine event_callback_nop_write

  module subroutine event_callback_nop (event_callback, i, event)
    class(event_callback_nop_t), intent(in) :: event_callback
    integer(i64), intent(in) :: i
    class(generic_event_t), intent(in) :: event
  end subroutine event_callback_nop


end submodule event_base_s

