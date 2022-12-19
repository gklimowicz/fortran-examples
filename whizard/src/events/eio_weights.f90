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

module eio_weights

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use event_base
  use event_handles, only: event_handle_t
  use eio_data
  use eio_base

  implicit none
  private

  public :: eio_weights_t

  type, extends (eio_t) :: eio_weights_t
     logical :: writing = .false.
     integer :: unit = 0
     logical :: pacify = .false.
   contains
     procedure :: set_parameters => eio_weights_set_parameters
     procedure :: write => eio_weights_write
     procedure :: final => eio_weights_final
     procedure :: init_out => eio_weights_init_out
     procedure :: init_in => eio_weights_init_in
     procedure :: switch_inout => eio_weights_switch_inout
     procedure :: output => eio_weights_output
     procedure :: input_i_prc => eio_weights_input_i_prc
     procedure :: input_event => eio_weights_input_event
     procedure :: skip => eio_weights_skip
  end type eio_weights_t


  interface
    module subroutine eio_weights_set_parameters (eio, pacify)
      class(eio_weights_t), intent(inout) :: eio
      logical, intent(in), optional :: pacify
    end subroutine eio_weights_set_parameters
    module subroutine eio_weights_write (object, unit)
      class(eio_weights_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine eio_weights_write
    module subroutine eio_weights_final (object)
      class(eio_weights_t), intent(inout) :: object
    end subroutine eio_weights_final
    module subroutine eio_weights_init_out &
         (eio, sample, data, success, extension)
      class(eio_weights_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_weights_init_out
    module subroutine eio_weights_init_in &
         (eio, sample, data, success, extension)
      class(eio_weights_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(inout), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_weights_init_in
    module subroutine eio_weights_switch_inout (eio, success)
      class(eio_weights_t), intent(inout) :: eio
      logical, intent(out), optional :: success
    end subroutine eio_weights_switch_inout
    module subroutine eio_weights_output &
         (eio, event, i_prc, reading, passed, pacify, event_handle)
      class(eio_weights_t), intent(inout) :: eio
      class(generic_event_t), intent(in), target :: event
      integer, intent(in) :: i_prc
      logical, intent(in), optional :: reading, passed, pacify
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_weights_output
    module subroutine eio_weights_input_i_prc (eio, i_prc, iostat)
      class(eio_weights_t), intent(inout) :: eio
      integer, intent(out) :: i_prc
      integer, intent(out) :: iostat
    end subroutine eio_weights_input_i_prc
    module subroutine eio_weights_input_event &
         (eio, event, iostat, event_handle)
      class(eio_weights_t), intent(inout) :: eio
      class(generic_event_t), intent(inout), target :: event
      integer, intent(out) :: iostat
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_weights_input_event
    module subroutine eio_weights_skip (eio, iostat)
      class(eio_weights_t), intent(inout) :: eio
      integer, intent(out) :: iostat
    end subroutine eio_weights_skip
  end interface

end module eio_weights
