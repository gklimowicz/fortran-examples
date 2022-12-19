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

submodule (prc_test) prc_test_s

  implicit none

contains

  module function prc_test_def_type_string () result (string)
    type(string_t) :: string
    string = "test_me"
  end function prc_test_def_type_string

  module subroutine prc_test_def_get_features (features)
    type(string_t), dimension(:), allocatable, intent(out) :: features
    allocate (features (0))
  end subroutine prc_test_def_get_features

  module subroutine prc_test_def_init (object, model_name, prt_in, prt_out)
    class(prc_test_def_t), intent(out) :: object
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: prt_in
    type(string_t), dimension(:), intent(in) :: prt_out
    object%model_name = model_name
    allocate (object%prt_in (size (prt_in)))
    object%prt_in = prt_in
    allocate (object%prt_out (size (prt_out)))
    object%prt_out = prt_out
  end subroutine prc_test_def_init

  module subroutine prc_test_def_write (object, unit)
    class(prc_test_def_t), intent(in) :: object
    integer, intent(in) :: unit
  end subroutine prc_test_def_write

  module subroutine prc_test_def_read (object, unit)
    class(prc_test_def_t), intent(out) :: object
    integer, intent(in) :: unit
  end subroutine prc_test_def_read

  module subroutine prc_test_def_connect (def, lib_driver, i, proc_driver)
    class(prc_test_def_t), intent(in) :: def
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    class(prc_core_driver_t), intent(inout) :: proc_driver
  end subroutine prc_test_def_connect

  module function prc_test_get_amplitude (p) result (amp)
    complex(default) :: amp
    real(default), dimension(:,:), intent(in) :: p
    amp = 1
  end function prc_test_get_amplitude

  module function prc_test_type_name () result (string)
    type(string_t) :: string
    string = "test_me"
  end function prc_test_type_name

  module subroutine prc_test_fill_constants (driver, data)
    class(prc_test_t), intent(in) :: driver
    type(process_constants_t), intent(out) :: data
    data%id = driver%id
    data%model_name = driver%model_name
    if (driver%scattering) then
       data%n_in  = 2
       data%n_out = 2
       data%n_flv = 1
       data%n_hel = 1
       data%n_col = 1
       data%n_cin = 2
       data%n_cf  = 1
       allocate (data%flv_state (4, 1))
       data%flv_state = 25
       allocate (data%hel_state (4, 1))
       data%hel_state = 0
       allocate (data%col_state (2, 4, 1))
       data%col_state = 0
       allocate (data%ghost_flag (4, 1))
       data%ghost_flag = .false.
       allocate (data%color_factors (1))
       data%color_factors = 1
       allocate (data%cf_index (2, 1))
       data%cf_index = 1
    else
       data%n_in  = 1
       data%n_out = 2
       data%n_flv = 1
       data%n_hel = 2
       data%n_col = 1
       data%n_cin = 2
       data%n_cf  = 1
       allocate (data%flv_state (3, 1))
       data%flv_state(:,1) = [25, 6, -6]
       allocate (data%hel_state (3, 2))
       data%hel_state(:,1) = [0, 1,-1]
       data%hel_state(:,2) = [0,-1, 1]
       allocate (data%col_state (2, 3, 1))
       data%col_state = reshape ([0,0, 1,0, 0,-1], [2,3,1])
       allocate (data%ghost_flag (3, 1))
       data%ghost_flag = .false.
       allocate (data%color_factors (1))
       data%color_factors = 3
       allocate (data%cf_index (2, 1))
       data%cf_index = 1
    end if
  end subroutine prc_test_fill_constants


end submodule prc_test_s

