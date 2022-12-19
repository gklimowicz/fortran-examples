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

submodule (phs_none) phs_none_s

  use io_units, only: given_output_unit
  use diagnostics, only: msg_message, msg_fatal

  implicit none

contains

  module subroutine phs_none_config_final (object)
    class(phs_none_config_t), intent(inout) :: object
  end subroutine phs_none_config_final

  module subroutine phs_none_config_write (object, unit, include_id)
    class(phs_none_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: include_id
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  &
         "Partonic phase-space configuration: non-functional dummy"
  end subroutine phs_none_config_write

  module subroutine phs_none_config_configure (phs_config, sqrts, &
       sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
       ignore_mismatch, nlo_type, subdir)
    class(phs_none_config_t), intent(inout) :: phs_config
    real(default), intent(in) :: sqrts
    logical, intent(in), optional :: sqrts_fixed
    logical, intent(in), optional :: lab_is_cm
    logical, intent(in), optional :: azimuthal_dependence
    logical, intent(in), optional :: rebuild
    logical, intent(in), optional :: ignore_mismatch
    integer, intent(in), optional :: nlo_type
    type(string_t), intent(in), optional :: subdir
  end subroutine phs_none_config_configure

  module subroutine phs_none_config_startup_message (phs_config, unit)
    class(phs_none_config_t), intent(in) :: phs_config
    integer, intent(in), optional :: unit
    call msg_message ("Phase space: none")
  end subroutine phs_none_config_startup_message

  module subroutine phs_none_write (object, unit, verbose)
    class(phs_none_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    u = given_output_unit (unit)
    write (u, "(A)")  "Partonic phase space: none"
  end subroutine phs_none_write

  module subroutine phs_none_final (object)
    class(phs_none_t), intent(inout) :: object
  end subroutine phs_none_final

  module subroutine phs_none_init (phs, phs_config)
    class(phs_none_t), intent(out) :: phs
    class(phs_config_t), intent(in), target :: phs_config
    call phs%base_init (phs_config)
  end subroutine phs_none_init

  module subroutine phs_none_evaluate_selected_channel (phs, c_in, r_in)
    class(phs_none_t), intent(inout) :: phs
    integer, intent(in) :: c_in
    real(default), intent(in), dimension(:) :: r_in
    call msg_fatal &
         ("Phase space: attempt to evaluate with the 'phs_none' method")
  end subroutine phs_none_evaluate_selected_channel

  module subroutine phs_none_evaluate_other_channels (phs, c_in)
    class(phs_none_t), intent(inout) :: phs
    integer, intent(in) :: c_in
  end subroutine phs_none_evaluate_other_channels

  module subroutine phs_none_inverse (phs)
    class(phs_none_t), intent(inout) :: phs
    call msg_fatal ("Phase space: attempt to evaluate inverse " // &
         "with the 'phs_none' method")
  end subroutine phs_none_inverse


end submodule phs_none_s

