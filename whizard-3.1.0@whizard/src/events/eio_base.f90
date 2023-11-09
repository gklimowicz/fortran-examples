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

module eio_base

  use kinds, only: i64
  use iso_varying_string, string_t => varying_string
  use model_data
  use event_base
  use event_handles, only: event_handle_t
  use eio_data

  implicit none
  private

  public :: eio_t

  type, abstract :: eio_t
     type(string_t) :: sample
     type(string_t) :: extension
     type(string_t) :: filename
     logical :: has_file = .false.
     logical :: split = .false.
     integer :: split_n_evt = 0
     integer :: split_n_kbytes = 0
     integer :: split_index = 0
     integer :: split_count = 0
     class(model_data_t), pointer :: fallback_model => null ()
   contains
     procedure (eio_write), deferred :: write
     procedure (eio_final), deferred :: final
     procedure :: set_splitting => eio_set_splitting
     procedure :: update_split_count => eio_update_split_count
     procedure :: set_filename => eio_set_filename
     procedure :: set_fallback_model => eio_set_fallback_model
     procedure (eio_init_out), deferred :: init_out
     procedure (eio_init_in), deferred :: init_in
     procedure (eio_switch_inout), deferred :: switch_inout
     procedure :: split_out => eio_split_out
     procedure :: file_size_kbytes => eio_file_size_kbytes
     procedure (eio_output), deferred :: output
     procedure (eio_input_i_prc), deferred :: input_i_prc
     procedure (eio_input_event), deferred :: input_event
     procedure (eio_skip), deferred :: skip
  end type eio_t


  abstract interface
     subroutine eio_write (object, unit)
       import
       class(eio_t), intent(in) :: object
       integer, intent(in), optional :: unit
     end subroutine eio_write
  end interface

  abstract interface
     subroutine eio_final (object)
       import
       class(eio_t), intent(inout) :: object
     end subroutine eio_final
  end interface

  abstract interface
     subroutine eio_init_out (eio, sample, data, success, extension)
       import
       class(eio_t), intent(inout) :: eio
       type(string_t), intent(in) :: sample
       type(event_sample_data_t), intent(in), optional :: data
       logical, intent(out), optional :: success
       type(string_t), intent(in), optional :: extension
     end subroutine eio_init_out
  end interface

  abstract interface
     subroutine eio_init_in (eio, sample, data, success, extension)
       import
       class(eio_t), intent(inout) :: eio
       type(string_t), intent(in) :: sample
       type(event_sample_data_t), intent(inout), optional :: data
       logical, intent(out), optional :: success
       type(string_t), intent(in), optional :: extension
     end subroutine eio_init_in
  end interface

  abstract interface
     subroutine eio_switch_inout (eio, success)
       import
       class(eio_t), intent(inout) :: eio
       logical, intent(out), optional :: success
     end subroutine eio_switch_inout
  end interface

  abstract interface
     subroutine eio_output &
          (eio, event, i_prc, reading, passed, pacify, event_handle)
       import
       class(eio_t), intent(inout) :: eio
       class(generic_event_t), intent(in), target :: event
       integer, intent(in) :: i_prc
       logical, intent(in), optional :: reading, passed, pacify
       class(event_handle_t), intent(inout), optional :: event_handle
     end subroutine eio_output
  end interface

  abstract interface
     subroutine eio_input_i_prc (eio, i_prc, iostat)
       import
       class(eio_t), intent(inout) :: eio
       integer, intent(out) :: i_prc
       integer, intent(out) :: iostat
     end subroutine eio_input_i_prc
  end interface

  abstract interface
     subroutine eio_input_event (eio, event, iostat, event_handle)
       import
       class(eio_t), intent(inout) :: eio
       class(generic_event_t), intent(inout), target :: event
       integer, intent(out) :: iostat
       class(event_handle_t), intent(inout), optional :: event_handle
     end subroutine eio_input_event
  end interface

  abstract interface
     subroutine eio_skip (eio, iostat)
       import
       class(eio_t), intent(inout) :: eio
       integer, intent(out) :: iostat
     end subroutine eio_skip
  end interface


  interface
    module subroutine eio_set_splitting (eio, data)
      class(eio_t), intent(inout) :: eio
      type(event_sample_data_t), intent(in) :: data
    end subroutine eio_set_splitting
    module subroutine eio_update_split_count (eio, increased)
      class(eio_t), intent(inout) :: eio
      logical, intent(out) :: increased
    end subroutine eio_update_split_count
    module subroutine eio_set_filename (eio)
      class(eio_t), intent(inout) :: eio
    end subroutine eio_set_filename
    module subroutine eio_set_fallback_model (eio, model)
      class(eio_t), intent(inout) :: eio
      class(model_data_t), intent(in), target :: model
    end subroutine eio_set_fallback_model
    module subroutine eio_split_out (eio)
      class(eio_t), intent(inout) :: eio
    end subroutine eio_split_out
    module function eio_file_size_kbytes (eio) result (kbytes)
      class(eio_t), intent(in) :: eio
      integer :: kbytes
    end function eio_file_size_kbytes
  end interface

end module eio_base
