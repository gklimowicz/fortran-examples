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

module rng_tao

  use kinds, only: default
  use tao_random_numbers !NODEP!

  use rng_base

  implicit none
  private

  public :: rng_tao_t
  public :: rng_tao_factory_t

  type, extends (rng_t) :: rng_tao_t
     integer :: seed = 0
     integer :: n_calls = 0
     type(tao_random_state) :: state
   contains
     procedure :: write => rng_tao_write
     procedure :: init => rng_tao_init
     procedure :: final => rng_tao_final
     procedure :: generate_single => rng_tao_generate_single
     procedure :: generate_array => rng_tao_generate_array
  end type rng_tao_t

  type, extends (rng_factory_t) :: rng_tao_factory_t
     integer(i16) :: s = 0
     integer(i16) :: i = 0
   contains
     procedure :: write => rng_tao_factory_write
     procedure :: init => rng_tao_factory_init
     procedure :: make => rng_tao_factory_make
  end type rng_tao_factory_t


  interface
    module subroutine rng_tao_write (rng, unit, indent)
      class(rng_tao_t), intent(in) :: rng
      integer, intent(in), optional :: unit, indent
    end subroutine rng_tao_write
    module subroutine rng_tao_init (rng, seed)
      class(rng_tao_t), intent(out) :: rng
      integer, intent(in), optional :: seed
    end subroutine rng_tao_init
    module subroutine rng_tao_final (rng)
      class(rng_tao_t), intent(inout) :: rng
    end subroutine rng_tao_final
    module subroutine rng_tao_generate_single (rng, x)
      class(rng_tao_t), intent(inout) :: rng
      real(default), intent(out) :: x
    end subroutine rng_tao_generate_single
    module subroutine rng_tao_generate_array (rng, x)
      class(rng_tao_t), intent(inout) :: rng
      real(default), dimension(:), intent(out) :: x
    end subroutine rng_tao_generate_array
    module subroutine rng_tao_factory_write (object, unit)
      class(rng_tao_factory_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine rng_tao_factory_write
    module subroutine rng_tao_factory_init (factory, seed)
      class(rng_tao_factory_t), intent(out) :: factory
      integer(i16), intent(in), optional :: seed
    end subroutine rng_tao_factory_init
  end interface

contains

  subroutine rng_tao_factory_make (factory, rng)
    class(rng_tao_factory_t), intent(inout) :: factory
    class(rng_t), intent(out), allocatable :: rng
    allocate (rng_tao_t :: rng)
    select type (rng)
    type is (rng_tao_t)
       call rng%init (factory%s * 65536 + factory%i)
       factory%i = int (factory%i + 1, kind = i16)
    end select
  end subroutine rng_tao_factory_make


end module rng_tao
