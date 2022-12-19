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

module eio_stdhep

  use kinds, only: i32, i64
  use iso_varying_string, string_t => varying_string
  use event_base
  use event_handles, only: event_handle_t
  use eio_data
  use eio_base

  implicit none
  private

  public :: eio_stdhep_t
  public :: eio_stdhep_hepevt_t
  public :: eio_stdhep_hepeup_t
  public :: eio_stdhep_hepev4_t
  public :: stdhep_init_out
  public :: stdhep_init_in
  public :: stdhep_write
  public :: stdhep_end

  type, abstract, extends (eio_t) :: eio_stdhep_t
     logical :: writing = .false.
     logical :: reading = .false.
     integer :: unit = 0
     logical :: keep_beams = .false.
     logical :: keep_remnants = .true.
     logical :: ensure_order = .false.
     logical :: recover_beams = .false.
     logical :: use_alphas_from_file = .false.
     logical :: use_scale_from_file = .false.
     integer, dimension(:), allocatable :: proc_num_id
     integer(i64) :: n_events_expected = 0
   contains
     procedure :: set_parameters => eio_stdhep_set_parameters
     procedure :: write => eio_stdhep_write
     procedure :: final => eio_stdhep_final
     procedure :: common_init => eio_stdhep_common_init
     procedure :: split_out => eio_stdhep_split_out
     procedure :: init_out => eio_stdhep_init_out
     procedure :: init_in => eio_stdhep_init_in
     procedure :: switch_inout => eio_stdhep_switch_inout
     procedure :: output => eio_stdhep_output
     procedure :: input_i_prc => eio_stdhep_input_i_prc
     procedure :: input_event => eio_stdhep_input_event
     procedure :: skip => eio_stdhep_skip
  end type eio_stdhep_t

  type, extends (eio_stdhep_t) :: eio_stdhep_hepevt_t
  end type eio_stdhep_hepevt_t

  type, extends (eio_stdhep_t) :: eio_stdhep_hepeup_t
  end type eio_stdhep_hepeup_t

  type, extends (eio_stdhep_t) :: eio_stdhep_hepev4_t
  end type eio_stdhep_hepev4_t


  interface
    module subroutine eio_stdhep_set_parameters (eio, &
         keep_beams, keep_remnants, ensure_order, recover_beams, &
         use_alphas_from_file, use_scale_from_file, extension)
      class(eio_stdhep_t), intent(inout) :: eio
      logical, intent(in), optional :: keep_beams
      logical, intent(in), optional :: keep_remnants
      logical, intent(in), optional :: ensure_order
      logical, intent(in), optional :: recover_beams
      logical, intent(in), optional :: use_alphas_from_file
      logical, intent(in), optional :: use_scale_from_file
      type(string_t), intent(in), optional :: extension
    end subroutine eio_stdhep_set_parameters
    module subroutine eio_stdhep_write (object, unit)
      class(eio_stdhep_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine eio_stdhep_write
    module subroutine eio_stdhep_final (object)
      class(eio_stdhep_t), intent(inout) :: object
    end subroutine eio_stdhep_final
    module subroutine eio_stdhep_common_init (eio, sample, data, extension)
      class(eio_stdhep_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
    end subroutine eio_stdhep_common_init
    module subroutine eio_stdhep_split_out (eio)
      class(eio_stdhep_t), intent(inout) :: eio
    end subroutine eio_stdhep_split_out
    module subroutine eio_stdhep_init_out &
         (eio, sample, data, success, extension)
      class(eio_stdhep_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_stdhep_init_out
    module subroutine eio_stdhep_init_in &
          (eio, sample, data, success, extension)
      class(eio_stdhep_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(inout), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_stdhep_init_in
    module subroutine eio_stdhep_switch_inout (eio, success)
      class(eio_stdhep_t), intent(inout) :: eio
      logical, intent(out), optional :: success
    end subroutine eio_stdhep_switch_inout
    module subroutine eio_stdhep_output &
         (eio, event, i_prc, reading, passed, pacify, event_handle)
      class(eio_stdhep_t), intent(inout) :: eio
      class(generic_event_t), intent(in), target :: event
      integer, intent(in) :: i_prc
      logical, intent(in), optional :: reading, passed, pacify
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_stdhep_output
    module subroutine eio_stdhep_input_i_prc (eio, i_prc, iostat)
      class(eio_stdhep_t), intent(inout) :: eio
      integer, intent(out) :: i_prc
      integer, intent(out) :: iostat
    end subroutine eio_stdhep_input_i_prc
    module subroutine eio_stdhep_input_event (eio, event, iostat, event_handle)
      class(eio_stdhep_t), intent(inout) :: eio
      class(generic_event_t), intent(inout), target :: event
      integer, intent(out) :: iostat
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_stdhep_input_event
    module subroutine eio_stdhep_skip (eio, iostat)
      class(eio_stdhep_t), intent(inout) :: eio
      integer, intent(out) :: iostat
    end subroutine eio_stdhep_skip
    module subroutine stdhep_init_out (file, title, nevt)
      character(len=*), intent(in) :: file, title
      integer(i64), intent(in) :: nevt
    end subroutine stdhep_init_out
    module subroutine stdhep_init_in (file, nevt)
      character(len=*), intent(in) :: file
      integer(i64), intent(out) :: nevt
    end subroutine stdhep_init_in
    module subroutine stdhep_write (ilbl)
      integer, intent(in) :: ilbl
    end subroutine stdhep_write
    module subroutine stdhep_read (ilbl, lok)
      integer, intent(out) :: ilbl, lok
    end subroutine stdhep_read
    module subroutine stdhep_end ()
    end subroutine stdhep_end
  end interface

end module eio_stdhep
