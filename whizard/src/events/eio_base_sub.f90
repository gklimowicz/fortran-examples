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

submodule (eio_base) eio_base_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine eio_set_splitting (eio, data)
    class(eio_t), intent(inout) :: eio
    type(event_sample_data_t), intent(in) :: data
    eio%split = data%split_n_evt > 0 .or. data%split_n_kbytes > 0
    if (eio%split) then
       eio%split_n_evt = data%split_n_evt
       eio%split_n_kbytes = data%split_n_kbytes
       eio%split_index = data%split_index
       eio%split_count = 0
    end if
  end subroutine eio_set_splitting

  module subroutine eio_update_split_count (eio, increased)
    class(eio_t), intent(inout) :: eio
    logical, intent(out) :: increased
    integer :: split_count_old
    if (eio%split_n_kbytes > 0) then
       split_count_old = eio%split_count
       eio%split_count = eio%file_size_kbytes () / eio%split_n_kbytes
       increased = eio%split_count > split_count_old
    end if
  end subroutine eio_update_split_count

  module subroutine eio_set_filename (eio)
    class(eio_t), intent(inout) :: eio
    character(32) :: buffer
    if (eio%split) then
       write (buffer, "(I0,'.')")  eio%split_index
       eio%filename = eio%sample // "." // trim (buffer) // eio%extension
       eio%has_file = .true.
    else
       eio%filename = eio%sample // "." // eio%extension
       eio%has_file = .true.
    end if
  end subroutine eio_set_filename

  module subroutine eio_set_fallback_model (eio, model)
    class(eio_t), intent(inout) :: eio
    class(model_data_t), intent(in), target :: model
    eio%fallback_model => model
  end subroutine eio_set_fallback_model

  module subroutine eio_split_out (eio)
    class(eio_t), intent(inout) :: eio
  end subroutine eio_split_out

  module function eio_file_size_kbytes (eio) result (kbytes)
    class(eio_t), intent(in) :: eio
    integer :: kbytes
    integer(i64) :: bytes
    if (eio%has_file) then
       inquire (file = char (eio%filename), size = bytes)
       if (bytes > 0) then
          kbytes = bytes / 1024
       else
          kbytes = 0
       end if
    else
       kbytes = 0
    end if
  end function eio_file_size_kbytes


end submodule eio_base_s

