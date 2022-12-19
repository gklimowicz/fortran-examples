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

module eio_lcio

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use event_base
  use event_handles, only: event_handle_t
  use eio_data
  use eio_base
  use lcio_interface

  implicit none
  private

  public :: eio_lcio_t

  type, extends (eio_t) :: eio_lcio_t
     logical :: writing = .false.
     logical :: reading = .false.
     type(event_sample_data_t) :: data
     logical :: recover_beams = .false.
     logical :: use_alphas_from_file = .false.
     logical :: use_scale_from_file = .false.
     logical :: proc_as_run_id = .true.
     integer :: n_alt = 0
     integer :: lcio_run_id = 0
     type(lcio_writer_t) :: lcio_writer
     type(lcio_reader_t) :: lcio_reader
     type(lcio_run_header_t) :: lcio_run_hdr
     type(lcio_event_t) :: lcio_event
     integer, dimension(:), allocatable :: proc_num_id
   contains
     procedure :: set_parameters => eio_lcio_set_parameters
     procedure :: write => eio_lcio_write
     procedure :: final => eio_lcio_final
     procedure :: split_out => eio_lcio_split_out
     procedure :: common_init => eio_lcio_common_init
     procedure :: init_out => eio_lcio_init_out
     procedure :: init_in => eio_lcio_init_in
     procedure :: switch_inout => eio_lcio_switch_inout
     procedure :: output => eio_lcio_output
     procedure :: input_i_prc => eio_lcio_input_i_prc
     procedure :: input_event => eio_lcio_input_event
     procedure :: skip => eio_lcio_skip
  end type eio_lcio_t


  interface
    module subroutine eio_lcio_set_parameters &
         (eio, recover_beams, use_alphas_from_file, use_scale_from_file, &
         extension, proc_as_run_id, lcio_run_id)
      class(eio_lcio_t), intent(inout) :: eio
      logical, intent(in), optional :: recover_beams
      logical, intent(in), optional :: use_alphas_from_file
      logical, intent(in), optional :: use_scale_from_file
      logical, intent(in), optional :: proc_as_run_id
      integer, intent(in), optional :: lcio_run_id
      type(string_t), intent(in), optional :: extension
    end subroutine eio_lcio_set_parameters
    module subroutine eio_lcio_write (object, unit)
      class(eio_lcio_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine eio_lcio_write
    module subroutine eio_lcio_final (object)
      class(eio_lcio_t), intent(inout) :: object
    end subroutine eio_lcio_final
    module subroutine eio_lcio_split_out (eio)
      class(eio_lcio_t), intent(inout) :: eio
    end subroutine eio_lcio_split_out
    module subroutine eio_lcio_common_init (eio, sample, data, extension)
      class(eio_lcio_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
    end subroutine eio_lcio_common_init
    module subroutine eio_lcio_init_out &
         (eio, sample, data, success, extension)
      class(eio_lcio_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_lcio_init_out
    module subroutine eio_lcio_init_in (eio, sample, data, success, extension)
      class(eio_lcio_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(inout), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_lcio_init_in
    module subroutine eio_lcio_switch_inout (eio, success)
      class(eio_lcio_t), intent(inout) :: eio
      logical, intent(out), optional :: success
    end subroutine eio_lcio_switch_inout
    module subroutine eio_lcio_output &
         (eio, event, i_prc, reading, passed, pacify, event_handle)
      class(eio_lcio_t), intent(inout) :: eio
      class(generic_event_t), intent(in), target :: event
      integer, intent(in) :: i_prc
      logical, intent(in), optional :: reading, passed, pacify
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_lcio_output
    module subroutine eio_lcio_input_i_prc (eio, i_prc, iostat)
      class(eio_lcio_t), intent(inout) :: eio
      integer, intent(out) :: i_prc
      integer, intent(out) :: iostat
    end subroutine eio_lcio_input_i_prc
    module subroutine eio_lcio_input_event (eio, event, iostat, event_handle)
      class(eio_lcio_t), intent(inout) :: eio
      class(generic_event_t), intent(inout), target :: event
      integer, intent(out) :: iostat
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_lcio_input_event
    module subroutine eio_lcio_skip (eio, iostat)
      class(eio_lcio_t), intent(inout) :: eio
      integer, intent(out) :: iostat
    end subroutine eio_lcio_skip
  end interface

end module eio_lcio
