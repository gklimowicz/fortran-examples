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

module dispatch_rng

  use kinds, only: i16
  use iso_varying_string, string_t => varying_string
  use variables
  use rng_base

  implicit none
  private

  public :: dispatch_rng_factory
  public :: dispatch_rng_factory_fallback
  public :: update_rng_seed_in_var_list

  procedure (dispatch_rng_factory), pointer :: &
       dispatch_rng_factory_fallback => null ()

  interface
    module subroutine dispatch_rng_factory &
         (rng_factory, var_list, next_rng_seed)
      class(rng_factory_t), allocatable, intent(inout) :: rng_factory
      type(var_list_t), intent(in) :: var_list
      integer, intent(out) :: next_rng_seed
    end subroutine dispatch_rng_factory
    module subroutine update_rng_seed_in_var_list (var_list, next_rng_seed)
      type(var_list_t), intent(inout), optional :: var_list
      integer, intent(in) :: next_rng_seed
    end subroutine update_rng_seed_in_var_list
  end interface

end module dispatch_rng
