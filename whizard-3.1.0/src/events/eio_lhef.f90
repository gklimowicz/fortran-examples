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

module eio_lhef

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use xml
  use event_base
  use event_handles, only: event_handle_t
  use eio_data
  use eio_base

  implicit none
  private

  public :: eio_lhef_t

  type, extends (eio_t) :: eio_lhef_t
     logical :: writing = .false.
     logical :: reading = .false.
     integer :: unit = 0
     type(event_sample_data_t) :: data
     type(cstream_t) :: cstream
     character(3) :: version = "1.0"
     logical :: keep_beams = .false.
     logical :: keep_remnants = .true.
     logical :: keep_virtuals = .false.
     logical :: recover_beams = .true.
     logical :: unweighted = .true.
     logical :: write_sqme_ref = .false.
     logical :: write_sqme_prc = .false.
     logical :: write_sqme_alt = .false.
     logical :: use_alphas_from_file = .false.
     logical :: use_scale_from_file = .false.
     integer :: n_alt = 0
     integer, dimension(:), allocatable :: proc_num_id
     integer :: i_weight_sqme = 0
     type(xml_tag_t) :: tag_lhef, tag_head, tag_init, tag_event
     type(xml_tag_t), allocatable :: tag_whiz_info
     type(xml_tag_t), allocatable :: tag_gen_n, tag_gen_v
     type(xml_tag_t), allocatable :: tag_generator, tag_xsecinfo
     type(xml_tag_t), allocatable :: tag_sqme_ref, tag_sqme_prc
     type(xml_tag_t), dimension(:), allocatable :: tag_sqme_alt, tag_wgts_alt
     type(xml_tag_t), allocatable :: tag_weight, tag_weightinfo, tag_weights
   contains
     procedure :: set_parameters => eio_lhef_set_parameters
     procedure :: write => eio_lhef_write
     procedure :: final => eio_lhef_final
     procedure :: common_init => eio_lhef_common_init
     procedure :: init_tags => eio_lhef_init_tags
     procedure :: init_out => eio_lhef_init_out
     procedure :: init_in => eio_lhef_init_in
     procedure :: merge_data => eio_merge_data
     procedure :: switch_inout => eio_lhef_switch_inout
     procedure :: split_out => eio_lhef_split_out
     procedure :: output => eio_lhef_output
     procedure :: input_i_prc => eio_lhef_input_i_prc
     procedure :: input_event => eio_lhef_input_event
     procedure :: skip => eio_lhef_skip
     procedure :: write_header => eio_lhef_write_header
     procedure :: write_footer => eio_lhef_write_footer
     procedure :: read_header => eio_lhef_read_header
     procedure :: read_init_10 => eio_lhef_read_init_10
     procedure :: write_init_20 => eio_lhef_write_init_20
     procedure :: read_init_20 => eio_lhef_read_init_20
     procedure :: write_event_20 => eio_lhef_write_event_20
     procedure :: read_event_20 => eio_lhef_read_event_20
     procedure :: write_init_30 => eio_lhef_write_init_30
     procedure :: read_init_30 => eio_lhef_read_init_30
     procedure :: write_event_30 => eio_lhef_write_event_30
     procedure :: read_event_30 => eio_lhef_read_event_30
  end type eio_lhef_t


  interface
    module subroutine eio_lhef_set_parameters (eio, &
         keep_beams, keep_remnants, recover_beams, &
         use_alphas_from_file, use_scale_from_file, &
         version, extension, write_sqme_ref, write_sqme_prc, write_sqme_alt)
      class(eio_lhef_t), intent(inout) :: eio
      logical, intent(in), optional :: keep_beams
      logical, intent(in), optional :: keep_remnants
      logical, intent(in), optional :: recover_beams
      logical, intent(in), optional :: use_alphas_from_file
      logical, intent(in), optional :: use_scale_from_file
      character(*), intent(in), optional :: version
      type(string_t), intent(in), optional :: extension
      logical, intent(in), optional :: write_sqme_ref
      logical, intent(in), optional :: write_sqme_prc
      logical, intent(in), optional :: write_sqme_alt
    end subroutine eio_lhef_set_parameters
    module subroutine eio_lhef_write (object, unit)
      class(eio_lhef_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine eio_lhef_write
    module subroutine eio_lhef_final (object)
      class(eio_lhef_t), intent(inout) :: object
    end subroutine eio_lhef_final
    module subroutine eio_lhef_common_init (eio, sample, data, extension)
      class(eio_lhef_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
    end subroutine eio_lhef_common_init
    module subroutine eio_lhef_init_tags (eio, data)
      class(eio_lhef_t), intent(inout) :: eio
      type(event_sample_data_t), intent(in) :: data
    end subroutine eio_lhef_init_tags
    module subroutine eio_lhef_init_out (eio, sample, data, success, extension)
      class(eio_lhef_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_lhef_init_out
    module subroutine eio_lhef_init_in (eio, sample, data, success, extension)
      class(eio_lhef_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(inout), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_lhef_init_in
    module subroutine eio_merge_data (eio, data, data_file)
      class(eio_lhef_t), intent(inout) :: eio
      type(event_sample_data_t), intent(inout) :: data
      type(event_sample_data_t), intent(in) :: data_file
    end subroutine eio_merge_data
    module subroutine eio_lhef_switch_inout (eio, success)
      class(eio_lhef_t), intent(inout) :: eio
      logical, intent(out), optional :: success
    end subroutine eio_lhef_switch_inout
    module subroutine eio_lhef_split_out (eio)
      class(eio_lhef_t), intent(inout) :: eio
    end subroutine eio_lhef_split_out
    module subroutine eio_lhef_output &
         (eio, event, i_prc, reading, passed, pacify, event_handle)
      class(eio_lhef_t), intent(inout) :: eio
      class(generic_event_t), intent(in), target :: event
      integer, intent(in) :: i_prc
      logical, intent(in), optional :: reading, passed, pacify
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_lhef_output
    module subroutine eio_lhef_input_i_prc (eio, i_prc, iostat)
      class(eio_lhef_t), intent(inout) :: eio
      integer, intent(out) :: i_prc
      integer, intent(out) :: iostat
    end subroutine eio_lhef_input_i_prc
    module subroutine eio_lhef_input_event (eio, event, iostat, event_handle)
      class(eio_lhef_t), intent(inout) :: eio
      class(generic_event_t), intent(inout), target :: event
      integer, intent(out) :: iostat
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_lhef_input_event
    module subroutine eio_lhef_skip (eio, iostat)
      class(eio_lhef_t), intent(inout) :: eio
      integer, intent(out) :: iostat
    end subroutine eio_lhef_skip
    module subroutine eio_lhef_write_header (eio, is_width)
      class(eio_lhef_t), intent(in) :: eio
      logical, intent(in), optional :: is_width
    end subroutine eio_lhef_write_header
    module subroutine eio_lhef_write_footer (eio)
      class(eio_lhef_t), intent(in) :: eio
    end subroutine eio_lhef_write_footer
    module subroutine eio_lhef_read_header (eio)
      class(eio_lhef_t), intent(inout) :: eio
    end subroutine eio_lhef_read_header
    module subroutine eio_lhef_read_init_10 (eio, data)
      class(eio_lhef_t), intent(in) :: eio
      type(event_sample_data_t), intent(out) :: data
    end subroutine eio_lhef_read_init_10
    module subroutine eio_lhef_write_init_20 (eio, data)
      class(eio_lhef_t), intent(in) :: eio
      type(event_sample_data_t), intent(in) :: data
    end subroutine eio_lhef_write_init_20
    module subroutine eio_lhef_read_init_20 (eio, data)
      class(eio_lhef_t), intent(inout) :: eio
      type(event_sample_data_t), intent(out) :: data
    end subroutine eio_lhef_read_init_20
    module subroutine eio_lhef_write_event_20 (eio, event)
      class(eio_lhef_t), intent(in) :: eio
      class(generic_event_t), intent(in) :: event
    end subroutine eio_lhef_write_event_20
    module subroutine eio_lhef_read_event_20 (eio, event)
      class(eio_lhef_t), intent(inout) :: eio
      class(generic_event_t), intent(inout) :: event
    end subroutine eio_lhef_read_event_20
    module subroutine eio_lhef_write_init_30 (eio, data)
      class(eio_lhef_t), intent(in) :: eio
      type(event_sample_data_t), intent(in) :: data
    end subroutine eio_lhef_write_init_30
    module subroutine eio_lhef_read_init_30 (eio, data)
      class(eio_lhef_t), intent(inout) :: eio
      type(event_sample_data_t), intent(out) :: data
    end subroutine eio_lhef_read_init_30
    module subroutine eio_lhef_write_event_30 (eio, event)
      class(eio_lhef_t), intent(in) :: eio
      class(generic_event_t), intent(in) :: event
    end subroutine eio_lhef_write_event_30
    module subroutine eio_lhef_read_event_30 (eio, event)
      class(eio_lhef_t), intent(inout) :: eio
      class(generic_event_t), intent(inout) :: event
    end subroutine eio_lhef_read_event_30
  end interface

end module eio_lhef
