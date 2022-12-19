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

module eio_ascii

  use iso_varying_string, string_t => varying_string
  use event_base
  use event_handles, only: event_handle_t
  use eio_data
  use eio_base

  implicit none
  private

  public :: eio_ascii_t
  public :: eio_ascii_ascii_t
  public :: eio_ascii_athena_t
  public :: eio_ascii_debug_t
  public :: eio_ascii_hepevt_t
  public :: eio_ascii_hepevt_verb_t
  public :: eio_ascii_lha_t
  public :: eio_ascii_lha_verb_t
  public :: eio_ascii_long_t
  public :: eio_ascii_mokka_t
  public :: eio_ascii_short_t

  type, abstract, extends (eio_t) :: eio_ascii_t
     logical :: writing = .false.
     integer :: unit = 0
     logical :: keep_beams = .false.
     logical :: keep_remnants = .true.
     logical :: ensure_order = .false.
   contains
     procedure :: set_parameters => eio_ascii_set_parameters
     procedure :: write => eio_ascii_write
     procedure :: final => eio_ascii_final
     procedure :: init_out => eio_ascii_init_out
     procedure :: check_normalization => eio_ascii_check_normalization
     procedure :: init_in => eio_ascii_init_in
     procedure :: switch_inout => eio_ascii_switch_inout
     procedure :: split_out => eio_ascii_split_out
     procedure :: output => eio_ascii_output
     procedure :: input_i_prc => eio_ascii_input_i_prc
     procedure :: input_event => eio_ascii_input_event
     procedure :: skip => eio_ascii_skip
  end type eio_ascii_t

  type, extends (eio_ascii_t) :: eio_ascii_ascii_t
  end type eio_ascii_ascii_t

  type, extends (eio_ascii_t) :: eio_ascii_athena_t
  end type eio_ascii_athena_t

  type, extends (eio_ascii_t) :: eio_ascii_debug_t
     logical :: show_process = .true.
     logical :: show_transforms = .true.
     logical :: show_decay = .true.
     logical :: verbose = .true.
  end type eio_ascii_debug_t

  type, extends (eio_ascii_t) :: eio_ascii_hepevt_t
  end type eio_ascii_hepevt_t

  type, extends (eio_ascii_t) :: eio_ascii_hepevt_verb_t
  end type eio_ascii_hepevt_verb_t

  type, extends (eio_ascii_t) :: eio_ascii_lha_t
  end type eio_ascii_lha_t

  type, extends (eio_ascii_t) :: eio_ascii_lha_verb_t
  end type eio_ascii_lha_verb_t

  type, extends (eio_ascii_t) :: eio_ascii_long_t
  end type eio_ascii_long_t

  type, extends (eio_ascii_t) :: eio_ascii_mokka_t
  end type eio_ascii_mokka_t

  type, extends (eio_ascii_t) :: eio_ascii_short_t
  end type eio_ascii_short_t


  interface
    module subroutine eio_ascii_set_parameters (eio, &
         keep_beams, keep_remnants, ensure_order, extension, &
         show_process, show_transforms, show_decay, verbose)
      class(eio_ascii_t), intent(inout) :: eio
      logical, intent(in), optional :: keep_beams
      logical, intent(in), optional :: keep_remnants
      logical, intent(in), optional :: ensure_order
      type(string_t), intent(in), optional :: extension
      logical, intent(in), optional :: show_process, show_transforms, show_decay
      logical, intent(in), optional :: verbose
    end subroutine eio_ascii_set_parameters
    module subroutine eio_ascii_write (object, unit)
      class(eio_ascii_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine eio_ascii_write
    module subroutine eio_ascii_final (object)
      class(eio_ascii_t), intent(inout) :: object
    end subroutine eio_ascii_final
    module subroutine eio_ascii_init_out &
         (eio, sample, data, success, extension)
      class(eio_ascii_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_ascii_init_out
    module subroutine eio_ascii_check_normalization (eio, data)
      class(eio_ascii_t), intent(in) :: eio
      type(event_sample_data_t), intent(in) :: data
    end subroutine eio_ascii_check_normalization
    module subroutine eio_ascii_init_in &
         (eio, sample, data, success, extension)
      class(eio_ascii_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(inout), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_ascii_init_in
    module subroutine eio_ascii_switch_inout (eio, success)
      class(eio_ascii_t), intent(inout) :: eio
      logical, intent(out), optional :: success
    end subroutine eio_ascii_switch_inout
    module subroutine eio_ascii_split_out (eio)
      class(eio_ascii_t), intent(inout) :: eio
    end subroutine eio_ascii_split_out
    module subroutine eio_ascii_output &
         (eio, event, i_prc, reading, passed, pacify, event_handle)
      class(eio_ascii_t), intent(inout) :: eio
      class(generic_event_t), intent(in), target :: event
      integer, intent(in) :: i_prc
      logical, intent(in), optional :: reading, passed, pacify
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_ascii_output
    module subroutine eio_ascii_input_i_prc (eio, i_prc, iostat)
      class(eio_ascii_t), intent(inout) :: eio
      integer, intent(out) :: i_prc
      integer, intent(out) :: iostat
    end subroutine eio_ascii_input_i_prc
    module subroutine eio_ascii_input_event &
         (eio, event, iostat, event_handle)
      class(eio_ascii_t), intent(inout) :: eio
      class(generic_event_t), intent(inout), target :: event
      integer, intent(out) :: iostat
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_ascii_input_event
    module subroutine eio_ascii_skip (eio, iostat)
       class(eio_ascii_t), intent(inout) :: eio
       integer, intent(out) :: iostat
    end subroutine eio_ascii_skip
  end interface

end module eio_ascii
