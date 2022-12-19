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

module event_streams

  use iso_varying_string, string_t => varying_string
  use events
  use event_handles, only: event_handle_t
  use eio_data
  use eio_base
  use rt_data

  use dispatch_transforms, only: dispatch_eio

  implicit none
  private

  public :: event_stream_array_t

  type :: event_stream_entry_t
     class(eio_t), allocatable :: eio
  end type event_stream_entry_t

  type :: event_stream_array_t
     type(event_stream_entry_t), dimension(:), allocatable :: entry
     integer :: i_in = 0
   contains
     procedure :: write => event_stream_array_write
     procedure :: is_valid => event_stream_array_is_valid
     procedure :: final => event_stream_array_final
     procedure :: init => event_stream_array_init
     procedure :: switch_inout => event_stream_array_switch_inout
     procedure :: output => event_stream_array_output
     procedure :: input_i_prc => event_stream_array_input_i_prc
     procedure :: input_event => event_stream_array_input_event
     procedure :: skip_eio_entry => event_stream_array_skip_eio_entry
     procedure :: has_input => event_stream_array_has_input
  end type event_stream_array_t


  interface
    module subroutine event_stream_array_write (object, unit)
      class(event_stream_array_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine event_stream_array_write
    module function event_stream_array_is_valid (es_array) result (flag)
      class(event_stream_array_t), intent(in) :: es_array
      logical :: flag
    end function event_stream_array_is_valid
    module subroutine event_stream_array_final (es_array)
      class(event_stream_array_t), intent(inout) :: es_array
    end subroutine event_stream_array_final
    module subroutine event_stream_array_init &
         (es_array, sample, stream_fmt, global, &
         data, input, input_sample, input_data, allow_switch, &
         checkpoint, callback, &
         error)
      class(event_stream_array_t), intent(out) :: es_array
      type(string_t), intent(in) :: sample
      type(string_t), dimension(:), intent(in) :: stream_fmt
      type(rt_data_t), intent(in) :: global
      type(event_sample_data_t), intent(inout), optional :: data
      type(string_t), intent(in), optional :: input
      type(string_t), intent(in), optional :: input_sample
      type(event_sample_data_t), intent(inout), optional :: input_data
      logical, intent(in), optional :: allow_switch
      integer, intent(in), optional :: checkpoint
      integer, intent(in), optional :: callback
      logical, intent(out), optional :: error
    end subroutine event_stream_array_init
    module subroutine event_stream_array_switch_inout (es_array)
      class(event_stream_array_t), intent(inout) :: es_array
    end subroutine event_stream_array_switch_inout
    module subroutine event_stream_array_output &
         (es_array, event, i_prc, event_index, passed, pacify, event_handle)
      class(event_stream_array_t), intent(inout) :: es_array
      type(event_t), intent(in), target :: event
      integer, intent(in) :: i_prc, event_index
      logical, intent(in), optional :: passed, pacify
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine event_stream_array_output
    module subroutine event_stream_array_input_i_prc (es_array, i_prc, iostat)
      class(event_stream_array_t), intent(inout) :: es_array
      integer, intent(out) :: i_prc
      integer, intent(out) :: iostat
    end subroutine event_stream_array_input_i_prc
    module subroutine event_stream_array_input_event &
         (es_array, event, iostat, event_handle)
      class(event_stream_array_t), intent(inout) :: es_array
      type(event_t), intent(inout), target :: event
      integer, intent(out) :: iostat
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine event_stream_array_input_event
    module subroutine event_stream_array_skip_eio_entry (es_array, iostat)
      class(event_stream_array_t), intent(inout) :: es_array
      integer, intent(out) :: iostat
    end subroutine event_stream_array_skip_eio_entry
    module function event_stream_array_has_input (es_array) result (flag)
      class(event_stream_array_t), intent(in) :: es_array
      logical :: flag
    end function event_stream_array_has_input
  end interface

end module event_streams
