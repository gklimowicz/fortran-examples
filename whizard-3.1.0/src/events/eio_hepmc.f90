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

module eio_hepmc

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use event_base
  use event_handles, only: event_handle_t
  use eio_data
  use eio_base
  use hepmc_interface

  implicit none
  private

  public :: eio_hepmc_t

  type, extends (eio_t) :: eio_hepmc_t
     logical :: writing = .false.
     logical :: reading = .false.
     type(event_sample_data_t) :: data
     ! logical :: keep_beams = .false.
     logical :: recover_beams = .false.
     logical :: use_alphas_from_file = .false.
     logical :: use_scale_from_file = .false.
     logical :: output_cross_section = .false.
     integer :: hepmc3_mode = HEPMC3_MODE_HEPMC3
     logical :: hepmc3_flows = .false.
     type(hepmc_iostream_t) :: iostream
     type(hepmc_event_t) :: hepmc_event
     integer, dimension(:), allocatable :: proc_num_id
   contains
     procedure :: set_parameters => eio_hepmc_set_parameters
     procedure :: write => eio_hepmc_write
     procedure :: final => eio_hepmc_final
     procedure :: split_out => eio_hepmc_split_out
     procedure :: common_init => eio_hepmc_common_init
     procedure :: init_out => eio_hepmc_init_out
     procedure :: init_in => eio_hepmc_init_in
     procedure :: switch_inout => eio_hepmc_switch_inout
     procedure :: output => eio_hepmc_output
     procedure :: input_i_prc => eio_hepmc_input_i_prc
     procedure :: input_event => eio_hepmc_input_event
     procedure :: skip => eio_hepmc_skip
  end type eio_hepmc_t


  interface
    module subroutine eio_hepmc_set_parameters &
         (eio, recover_beams, use_alphas_from_file, &
         use_scale_from_file, extension, output_cross_section, &
         hepmc3_mode, hepmc3_write_flows)
      class(eio_hepmc_t), intent(inout) :: eio
      logical, intent(in), optional :: recover_beams
      logical, intent(in), optional :: use_alphas_from_file
      logical, intent(in), optional :: use_scale_from_file
      logical, intent(in), optional :: output_cross_section
      type(string_t), intent(in), optional :: extension
      integer, intent(in), optional :: hepmc3_mode
      logical ,intent(in), optional :: hepmc3_write_flows
    end subroutine eio_hepmc_set_parameters
    module subroutine eio_hepmc_write (object, unit)
      class(eio_hepmc_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine eio_hepmc_write
    module subroutine eio_hepmc_final (object)
      class(eio_hepmc_t), intent(inout) :: object
    end subroutine eio_hepmc_final
    module subroutine eio_hepmc_split_out (eio)
      class(eio_hepmc_t), intent(inout) :: eio
    end subroutine eio_hepmc_split_out
    module subroutine eio_hepmc_common_init (eio, sample, data, extension)
      class(eio_hepmc_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
    end subroutine eio_hepmc_common_init
    module subroutine eio_hepmc_init_out (eio, sample, data, success, extension)
      class(eio_hepmc_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(in), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_hepmc_init_out
    module subroutine eio_hepmc_init_in (eio, sample, data, success, extension)
      class(eio_hepmc_t), intent(inout) :: eio
      type(string_t), intent(in) :: sample
      type(string_t), intent(in), optional :: extension
      type(event_sample_data_t), intent(inout), optional :: data
      logical, intent(out), optional :: success
    end subroutine eio_hepmc_init_in
    module subroutine eio_hepmc_switch_inout (eio, success)
      class(eio_hepmc_t), intent(inout) :: eio
      logical, intent(out), optional :: success
    end subroutine eio_hepmc_switch_inout
    module subroutine eio_hepmc_output &
         (eio, event, i_prc, reading, passed, pacify, event_handle)
      class(eio_hepmc_t), intent(inout) :: eio
      class(generic_event_t), intent(in), target :: event
      integer, intent(in) :: i_prc
      logical, intent(in), optional :: reading, passed, pacify
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_hepmc_output
    module subroutine eio_hepmc_input_i_prc (eio, i_prc, iostat)
      class(eio_hepmc_t), intent(inout) :: eio
      integer, intent(out) :: i_prc
      integer, intent(out) :: iostat
    end subroutine eio_hepmc_input_i_prc
    module subroutine eio_hepmc_input_event (eio, event, iostat, event_handle)
      class(eio_hepmc_t), intent(inout) :: eio
      class(generic_event_t), intent(inout), target :: event
      integer, intent(out) :: iostat
      class(event_handle_t), intent(inout), optional :: event_handle
    end subroutine eio_hepmc_input_event
    module subroutine eio_hepmc_skip (eio, iostat)
      class(eio_hepmc_t), intent(inout) :: eio
      integer, intent(out) :: iostat
    end subroutine eio_hepmc_skip
  end interface

end module eio_hepmc
