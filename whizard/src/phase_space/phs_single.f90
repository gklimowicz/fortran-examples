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

module phs_single

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use phs_base

  implicit none
  private

  public :: phs_single_config_t
  public :: phs_single_t

  type, extends (phs_config_t) :: phs_single_config_t
  contains
    procedure :: final => phs_single_config_final
    procedure :: write => phs_single_config_write
    procedure :: configure => phs_single_config_configure
    procedure :: startup_message => phs_single_config_startup_message
    procedure, nopass :: allocate_instance => phs_single_config_allocate_instance
  end type phs_single_config_t

  type, extends (phs_t) :: phs_single_t
   contains
     procedure :: write => phs_single_write
     procedure :: final => phs_single_final
     procedure :: init => phs_single_init
     procedure :: compute_factor => phs_single_compute_factor
     procedure :: evaluate_selected_channel => phs_single_evaluate_selected_channel
     procedure :: evaluate_other_channels => phs_single_evaluate_other_channels
     procedure :: decay_p => phs_single_decay_p
     procedure :: inverse => phs_single_inverse
  end type phs_single_t


  interface
    module subroutine phs_single_config_final (object)
      class(phs_single_config_t), intent(inout) :: object
    end subroutine phs_single_config_final
    module subroutine phs_single_config_write (object, unit, include_id)
      class(phs_single_config_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: include_id
    end subroutine phs_single_config_write
    module subroutine phs_single_config_configure (phs_config, sqrts, &
         sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
         ignore_mismatch, nlo_type, subdir)
      class(phs_single_config_t), intent(inout) :: phs_config
      real(default), intent(in) :: sqrts
      logical, intent(in), optional :: sqrts_fixed
      logical, intent(in), optional :: lab_is_cm
      logical, intent(in), optional :: azimuthal_dependence
      logical, intent(in), optional :: rebuild
      logical, intent(in), optional :: ignore_mismatch
      integer, intent(in), optional :: nlo_type
      type(string_t), intent(in), optional :: subdir
    end subroutine phs_single_config_configure
    module subroutine phs_single_config_startup_message (phs_config, unit)
      class(phs_single_config_t), intent(in) :: phs_config
      integer, intent(in), optional :: unit
    end subroutine phs_single_config_startup_message
    module subroutine phs_single_write (object, unit, verbose)
      class(phs_single_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine phs_single_write
    module subroutine phs_single_final (object)
      class(phs_single_t), intent(inout) :: object
    end subroutine phs_single_final
    module subroutine phs_single_init (phs, phs_config)
      class(phs_single_t), intent(out) :: phs
      class(phs_config_t), intent(in), target :: phs_config
    end subroutine phs_single_init
    module subroutine phs_single_compute_factor (phs)
      class(phs_single_t), intent(inout) :: phs
    end subroutine phs_single_compute_factor
    module subroutine phs_single_evaluate_selected_channel (phs, c_in, r_in)
      class(phs_single_t), intent(inout) :: phs
      integer, intent(in) :: c_in
      real(default), intent(in), dimension(:) :: r_in
    end subroutine phs_single_evaluate_selected_channel
    module subroutine phs_single_evaluate_other_channels (phs, c_in)
      class(phs_single_t), intent(inout) :: phs
      integer, intent(in) :: c_in
    end subroutine phs_single_evaluate_other_channels
    module function phs_single_decay_p (phs) result (p)
      class(phs_single_t), intent(in) :: phs
      type(vector4_t), dimension(2) :: p
    end function phs_single_decay_p
    module subroutine phs_single_inverse (phs)
      class(phs_single_t), intent(inout) :: phs
    end subroutine phs_single_inverse
  end interface

contains

  subroutine phs_single_config_allocate_instance (phs)
    class(phs_t), intent(inout), pointer :: phs
    allocate (phs_single_t :: phs)
  end subroutine phs_single_config_allocate_instance


end module phs_single
