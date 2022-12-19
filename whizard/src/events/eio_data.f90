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

module eio_data

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  use event_base

  implicit none
  private

  public :: event_sample_data_t

  type :: event_sample_data_t
     character(32) :: md5sum_prc = ""
     character(32) :: md5sum_cfg = ""
     logical :: unweighted = .true.
     logical :: negative_weights = .false.
     integer :: norm_mode = NORM_UNDEFINED
     integer :: n_beam = 0
     integer, dimension(2) :: pdg_beam = 0
     real(default), dimension(2) :: energy_beam = 0
     integer :: n_proc = 0
     integer :: n_evt = 0
     integer :: nlo_multiplier = 1
     integer :: split_n_evt = 0
     integer :: split_n_kbytes = 0
     integer :: split_index = 0
     real(default) :: total_cross_section = 0
     integer, dimension(:), allocatable :: proc_num_id
     integer :: n_alt = 0
     character(32), dimension(:), allocatable :: md5sum_alt
     real(default), dimension(:), allocatable :: cross_section
     real(default), dimension(:), allocatable :: error
   contains
     procedure :: init => event_sample_data_init
     procedure :: write => event_sample_data_write
  end type event_sample_data_t


  interface
    module subroutine event_sample_data_init (data, n_proc, n_alt)
      class(event_sample_data_t), intent(out) :: data
      integer, intent(in) :: n_proc
      integer, intent(in), optional :: n_alt
    end subroutine event_sample_data_init
    module subroutine event_sample_data_write (data, unit)
      class(event_sample_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
    end subroutine event_sample_data_write
  end interface

end module eio_data
