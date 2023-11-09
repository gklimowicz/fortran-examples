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

submodule (iterations) iterations_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine iterations_list_init &
       (it_list, n_it, n_calls, adapt, adapt_code, adapt_grids, adapt_weights)
    class(iterations_list_t), intent(inout) :: it_list
    integer, dimension(:), intent(in) :: n_it, n_calls
    logical, dimension(:), intent(in), optional :: adapt
    type(string_t), dimension(:), intent(in), optional :: adapt_code
    logical, dimension(:), intent(in), optional :: adapt_grids, adapt_weights
    integer :: i
    it_list%n_pass = size (n_it)
    if (allocated (it_list%pass)) deallocate (it_list%pass)
    allocate (it_list%pass (it_list%n_pass))
    it_list%pass%n_it = n_it
    it_list%pass%n_calls = n_calls
    if (present (adapt)) then
       it_list%pass%custom_adaptation = adapt
       do i = 1, it_list%n_pass
          if (adapt(i)) then
             if (verify (adapt_code(i), "wg") /= 0) then
                call msg_error ("iteration specification: " &
                     // "adaptation code letters must be 'w' or 'g'")
             end if
             it_list%pass(i)%adapt_grids = scan (adapt_code(i), "g") /= 0
             it_list%pass(i)%adapt_weights = scan (adapt_code(i), "w") /= 0
          end if
       end do
    else if (present (adapt_grids) .and. present (adapt_weights)) then
       it_list%pass%custom_adaptation = .true.
       it_list%pass%adapt_grids = adapt_grids
       it_list%pass%adapt_weights = adapt_weights
    end if
    do i = 1, it_list%n_pass - 1
       if (.not. it_list%pass(i)%custom_adaptation) then
          it_list%pass(i)%adapt_grids = .true.
          it_list%pass(i)%adapt_weights = .true.
       end if
    end do
  end subroutine iterations_list_init

  module subroutine iterations_list_clear (it_list)
    class(iterations_list_t), intent(inout) :: it_list
    it_list%n_pass = 0
    deallocate (it_list%pass)
  end subroutine iterations_list_clear

  module subroutine iterations_list_write (it_list, unit)
    class(iterations_list_t), intent(in) :: it_list
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(A)")  char (it_list%to_string ())
  end subroutine iterations_list_write

  module function iterations_list_to_string (it_list) result (buffer)
    class(iterations_list_t), intent(in) :: it_list
    type(string_t) :: buffer
    character(30) :: ibuf
    integer :: i
    buffer = "iterations = "
    if (it_list%n_pass > 0) then
       do i = 1, it_list%n_pass
          if (i > 1)  buffer = buffer // ", "
          write (ibuf, "(I0,':',I0)") &
               it_list%pass(i)%n_it, it_list%pass(i)%n_calls
          buffer = buffer // trim (ibuf)
          if (it_list%pass(i)%custom_adaptation &
               .or. it_list%pass(i)%adapt_grids &
               .or. it_list%pass(i)%adapt_weights) then
             buffer = buffer // ':"'
             if (it_list%pass(i)%adapt_grids)  buffer = buffer // "g"
             if (it_list%pass(i)%adapt_weights)  buffer = buffer // "w"
             buffer = buffer // '"'
          end if
       end do
    else
       buffer = buffer // "[undefined]"
    end if
  end function iterations_list_to_string

  module function iterations_list_get_n_pass (it_list) result (n_pass)
    class(iterations_list_t), intent(in) :: it_list
    integer :: n_pass
    n_pass = it_list%n_pass
  end function iterations_list_get_n_pass

  module function iterations_list_get_n_calls (it_list, pass) result (n_calls)
    class(iterations_list_t), intent(in) :: it_list
    integer :: n_calls
    integer, intent(in) :: pass
    if (pass <= it_list%n_pass) then
       n_calls = it_list%pass(pass)%n_calls
    else
       n_calls = 0
    end if
  end function iterations_list_get_n_calls

  module subroutine iterations_list_set_n_calls (it_list, pass, n_calls)
    class(iterations_list_t), intent(inout) :: it_list
    integer, intent(in) :: pass, n_calls
    it_list%pass(pass)%n_calls = n_calls
  end subroutine iterations_list_set_n_calls

  module function iterations_list_adapt_grids (it_list, pass) result (flag)
    logical :: flag
    class(iterations_list_t), intent(in) :: it_list
    integer, intent(in) :: pass
    if (pass <= it_list%n_pass) then
       flag = it_list%pass(pass)%adapt_grids
    else
       flag = .false.
    end if
  end function iterations_list_adapt_grids

  module function iterations_list_adapt_weights (it_list, pass) result (flag)
    logical :: flag
    class(iterations_list_t), intent(in) :: it_list
    integer, intent(in) :: pass
    if (pass <= it_list%n_pass) then
       flag = it_list%pass(pass)%adapt_weights
    else
       flag = .false.
    end if
  end function iterations_list_adapt_weights

  module function iterations_list_get_n_it (it_list, pass) result (n_it)
    class(iterations_list_t), intent(in) :: it_list
    integer :: n_it
    integer, intent(in) :: pass
    if (pass <= it_list%n_pass) then
       n_it = it_list%pass(pass)%n_it
    else
       n_it = 0
    end if
  end function iterations_list_get_n_it


end submodule iterations_s

