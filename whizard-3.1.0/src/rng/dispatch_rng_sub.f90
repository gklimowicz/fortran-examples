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

submodule (dispatch_rng) dispatch_rng_s

  use diagnostics
  use rng_tao
  use rng_stream

  implicit none



contains

  module subroutine dispatch_rng_factory &
       (rng_factory, var_list, next_rng_seed)
    class(rng_factory_t), allocatable, intent(inout) :: rng_factory
    type(var_list_t), intent(in) :: var_list
    integer, intent(out) :: next_rng_seed
    type(var_list_t) :: local
    type(string_t) :: rng_method
    integer :: seed
    character(30) :: buffer
    integer(i16) :: s
    rng_method = var_list%get_sval (var_str ("$rng_method"))
    seed = var_list%get_ival (var_str ("seed"))
    s = int (mod (seed, 32768), i16)
    select case (char (rng_method))
    case ("tao")
       allocate (rng_tao_factory_t :: rng_factory)
       call msg_message ("RNG: Initializing TAO random-number generator")
       next_rng_seed = seed + 1
    case ("rng_stream")
       allocate (rng_stream_factory_t :: rng_factory)
       call msg_message ("RNG: Initializing RNG Stream random-number generator")
       next_rng_seed = seed + 1
    case default
       if (associated (dispatch_rng_factory_fallback)) then
          call dispatch_rng_factory_fallback &
               (rng_factory, var_list, next_rng_seed)
       end if
       if (.not. allocated (rng_factory)) then
          call msg_fatal ("Random-number generator '" &
               // char (rng_method) // "' not implemented")
       end if
    end select
    write (buffer, "(I0)")  s
    call msg_message ("RNG: Setting seed for random-number generator to " &
            // trim (buffer))
    call rng_factory%init (s)
  end subroutine dispatch_rng_factory

  module subroutine update_rng_seed_in_var_list (var_list, next_rng_seed)
    type(var_list_t), intent(inout), optional :: var_list
    integer, intent(in) :: next_rng_seed
    if (present (var_list)) then 
       call var_list%set_int (var_str ("seed"), next_rng_seed, is_known=.true.)
    end if
  end subroutine update_rng_seed_in_var_list


end submodule dispatch_rng_s

