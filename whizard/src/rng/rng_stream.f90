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

module rng_stream

  use kinds, only: default
  use kinds, only: i16
  use rng_base

  implicit none
  private

  public :: rng_stream_t
  public :: rng_stream_factory_t

  type, extends(rng_t) :: rng_stream_t
     private
     logical :: antithetic = .false.
     logical :: increased_precision = .false.
     real(default), dimension(3, 2) :: current_position = 1234._default
     real(default), dimension(3, 2) :: beginning_substream = 1234._default
     real(default), dimension(3, 2) :: initial_stream = 1234._default
   contains
     procedure, public :: write => rng_stream_write
     procedure, public :: init => rng_stream_init
     procedure, public :: final => rng_stream_final
     procedure, public :: set_initial_stream => rng_stream_set_initial_stream
     procedure, public :: generate_single => rng_stream_generate_single
     procedure, public :: generate_array => rng_stream_generate_array
     procedure, public :: reset_stream => rng_stream_reset_stream
     procedure, public :: reset_substream => rng_stream_reset_substream
     procedure, public :: next_substream => rng_stream_next_substream
     procedure, public :: advance_state => rng_stream_advance_state
  end type rng_stream_t

  type, extends(rng_factory_t) :: rng_stream_factory_t
     real(default), dimension(3, 2) :: seed = 1234._default
   contains
     procedure, public :: write => rng_stream_factory_write
     procedure, public :: init => rng_stream_factory_init
     procedure, public :: make => rng_stream_factory_make
     procedure, private :: advance_seed => rng_stream_factory_advance_seed
  end type rng_stream_factory_t


  interface
    module subroutine rng_stream_write (rng, unit, indent)
      class(rng_stream_t), intent(in) :: rng
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: indent
    end subroutine rng_stream_write
    module subroutine rng_stream_init (rng, seed)
      class(rng_stream_t), intent(out) :: rng
      integer, intent(in), optional :: seed
    end subroutine rng_stream_init
    module subroutine rng_stream_final (rng)
      class(rng_stream_t), intent(inout) :: rng
    end subroutine rng_stream_final
    module subroutine rng_stream_set_initial_stream (rng, seed)
      class(rng_stream_t), intent(inout) :: rng
      real(default), dimension(3, 2), intent(in) :: seed
    end subroutine rng_stream_set_initial_stream
    module subroutine rng_stream_generate_single (rng, x)
      class(rng_stream_t), intent(inout) :: rng
      real(default), intent(out) :: x
    end subroutine rng_stream_generate_single
    module subroutine rng_stream_generate_array (rng, x)
      class(rng_stream_t), intent(inout) :: rng
      real(default), dimension(:), intent(out) :: x
    end subroutine rng_stream_generate_array
    module subroutine rng_stream_reset_stream (rng)
      class(rng_stream_t), intent(inout) :: rng
    end subroutine rng_stream_reset_stream
    module subroutine rng_stream_reset_substream (rng)
      class(rng_stream_t), intent(inout) :: rng
    end subroutine rng_stream_reset_substream
    module subroutine rng_stream_next_substream (rng)
      class(rng_stream_t), intent(inout) :: rng
    end subroutine rng_stream_next_substream
    module subroutine rng_stream_advance_state (rng, n)
      class(rng_stream_t), intent(inout) :: rng
      integer, intent(in) :: n
    end subroutine rng_stream_advance_state
    module subroutine rng_stream_factory_write (object, unit)
      class(rng_stream_factory_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine rng_stream_factory_write
    module subroutine rng_stream_factory_init (factory, seed)
      class(rng_stream_factory_t), intent(out) :: factory
      integer(i16), intent(in), optional :: seed
    end subroutine rng_stream_factory_init
    module subroutine rng_stream_factory_advance_seed (factory)
      class(rng_stream_factory_t), intent(inout) :: factory
    end subroutine rng_stream_factory_advance_seed
  end interface

contains

  subroutine rng_stream_factory_make (factory, rng)
    class(rng_stream_factory_t), intent(inout) :: factory
    class(rng_t), intent(out), allocatable :: rng
    allocate (rng_stream_t :: rng)
    select type (rng)
    type is (rng_stream_t)
       call rng%init ()
       call rng%set_initial_stream (factory%seed)
    end select
    call factory%advance_seed ()
  end subroutine rng_stream_factory_make


end module rng_stream
