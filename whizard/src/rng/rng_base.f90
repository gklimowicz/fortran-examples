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

module rng_base

  use kinds, only: default
  use kinds, only: i16

  implicit none
  private

  public :: rng_t
  public :: rng_factory_t

  type, abstract :: rng_t
   contains
     procedure (rng_init), deferred :: init
     procedure (rng_final), deferred :: final
     procedure (rng_write), deferred :: write
     generic :: generate => generate_single, generate_array
     procedure (rng_generate_single), deferred :: generate_single
     procedure (rng_generate_array), deferred :: generate_array
     generic :: generate_gaussian => &
          rng_generate_gaussian_single, rng_generate_gaussian_array
     procedure, private :: rng_generate_gaussian_single
     procedure, private :: rng_generate_gaussian_array
  end type rng_t

  type, abstract :: rng_factory_t
   contains
     procedure (rng_factory_write), deferred :: write
     procedure (rng_factory_init), deferred :: init
     procedure (rng_factory_make), deferred :: make
  end type rng_factory_t


  abstract interface
     subroutine rng_init (rng, seed)
       import
       class(rng_t), intent(out) :: rng
       integer, intent(in), optional :: seed
     end subroutine rng_init
  end interface

  abstract interface
     subroutine rng_final (rng)
       import
       class(rng_t), intent(inout) :: rng
     end subroutine rng_final
  end interface

  abstract interface
     subroutine rng_write (rng, unit, indent)
       import
       class(rng_t), intent(in) :: rng
       integer, intent(in), optional :: unit, indent
     end subroutine rng_write
  end interface

  abstract interface
     subroutine rng_generate_single (rng, x)
       import
       class(rng_t), intent(inout) :: rng
       real(default), intent(out) :: x
     end subroutine rng_generate_single
  end interface

  abstract interface
     subroutine rng_generate_array (rng, x)
       import
       class(rng_t), intent(inout) :: rng
       real(default), dimension(:), intent(out) :: x
     end subroutine rng_generate_array
  end interface

  abstract interface
     subroutine rng_factory_write (object, unit)
       import
       class(rng_factory_t), intent(in) :: object
       integer, intent(in), optional :: unit
     end subroutine rng_factory_write
  end interface

  abstract interface
     subroutine rng_factory_init (factory, seed)
       import
       class(rng_factory_t), intent(out) :: factory
       integer(i16), intent(in), optional :: seed
     end subroutine rng_factory_init
  end interface

  abstract interface
     subroutine rng_factory_make (factory, rng)
       import
       class(rng_factory_t), intent(inout) :: factory
       class(rng_t), intent(out), allocatable :: rng
     end subroutine rng_factory_make
  end interface


  interface
    module subroutine rng_generate_gaussian_single (rng, x)
      class(rng_t), intent(inout) :: rng
      real(default), intent(out) :: x
    end subroutine rng_generate_gaussian_single
    module subroutine rng_generate_gaussian_array (rng, x)
      class(rng_t), intent(inout) :: rng
      real(default), dimension(:), intent(out) :: x
    end subroutine rng_generate_gaussian_array
  end interface

end module rng_base
