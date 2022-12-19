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

module phs_rambo

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use phs_base

  implicit none
  private

  type, extends (phs_config_t) :: phs_rambo_config_t
  contains
    procedure :: final => phs_rambo_config_final
    procedure :: write => phs_rambo_config_write
    procedure :: configure => phs_rambo_config_configure
    procedure :: startup_message => phs_rambo_config_startup_message
    procedure, nopass :: allocate_instance => phs_rambo_config_allocate_instance
  end type phs_rambo_config_t

  type, extends (phs_t) :: phs_rambo_t
     real(default), dimension(:), allocatable :: k
     real(default), dimension(:), allocatable :: m
   contains
    procedure :: write => phs_rambo_write
    procedure :: final => phs_rambo_final
    procedure :: init => phs_rambo_init
    procedure :: evaluate_selected_channel => phs_rambo_evaluate_selected_channel
    procedure :: evaluate_other_channels => phs_rambo_evaluate_other_channels
    procedure, private :: decay_intermediate => phs_rambo_decay_intermediate
    procedure, private :: generate_intermediates => &
         phs_rambo_generate_intermediates
    procedure, private :: invert_intermediates => phs_rambo_invert_intermediates
    procedure :: inverse => phs_rambo_inverse
  end type phs_rambo_t


  public :: phs_rambo_config_t
  public :: phs_rambo_t

  interface
    module subroutine phs_rambo_config_final (object)
      class(phs_rambo_config_t), intent(inout) :: object
    end subroutine phs_rambo_config_final
    module subroutine phs_rambo_config_write (object, unit, include_id)
      class(phs_rambo_config_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: include_id
    end subroutine phs_rambo_config_write
    module subroutine phs_rambo_config_configure (phs_config, sqrts, &
         sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
         ignore_mismatch, nlo_type, subdir)
      class(phs_rambo_config_t), intent(inout) :: phs_config
      real(default), intent(in) :: sqrts
      logical, intent(in), optional :: sqrts_fixed
      logical, intent(in), optional :: lab_is_cm
      logical, intent(in), optional :: azimuthal_dependence
      logical, intent(in), optional :: rebuild
      logical, intent(in), optional :: ignore_mismatch
      integer, intent(in), optional :: nlo_type
      type(string_t), intent(in), optional :: subdir
    end subroutine phs_rambo_config_configure
    module subroutine phs_rambo_config_startup_message (phs_config, unit)
      class(phs_rambo_config_t), intent(in) :: phs_config
      integer, intent(in), optional :: unit
    end subroutine phs_rambo_config_startup_message
    module subroutine phs_rambo_write (object, unit, verbose)
      class(phs_rambo_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine phs_rambo_write
    module subroutine phs_rambo_final (object)
      class(phs_rambo_t), intent(inout) :: object
    end subroutine phs_rambo_final
    module subroutine phs_rambo_init (phs, phs_config)
      class(phs_rambo_t), intent(out) :: phs
      class(phs_config_t), intent(in), target :: phs_config
    end subroutine phs_rambo_init
    module subroutine phs_rambo_evaluate_selected_channel (phs, c_in, r_in)
      class(phs_rambo_t), intent(inout) :: phs
      integer, intent(in) :: c_in
      real(default), intent(in), dimension(:) :: r_in
    end subroutine phs_rambo_evaluate_selected_channel
    module subroutine phs_rambo_evaluate_other_channels (phs, c_in)
      class(phs_rambo_t), intent(inout) :: phs
      integer, intent(in) :: c_in
    end subroutine phs_rambo_evaluate_other_channels
    module subroutine phs_rambo_decay_intermediate (phs, i, r_angle, p)
      class(phs_rambo_t), intent(in) :: phs
      integer, intent(in) :: i
      real(default), dimension(2), intent(in) :: r_angle
      type(vector4_t), dimension(2), intent(out) :: p
    end subroutine phs_rambo_decay_intermediate
    module subroutine phs_rambo_generate_intermediates (phs, r)
      class(phs_rambo_t), intent(inout) :: phs
      real(default), dimension(:), intent(in) :: r
    end subroutine phs_rambo_generate_intermediates
    module subroutine phs_rambo_invert_intermediates (phs)
      class(phs_rambo_t), intent(inout) :: phs
    end subroutine phs_rambo_invert_intermediates
    module subroutine phs_rambo_inverse (phs)
      class(phs_rambo_t), intent(inout) :: phs
    end subroutine phs_rambo_inverse
  end interface

contains

  subroutine phs_rambo_config_allocate_instance (phs)
    class(phs_t), intent(inout), pointer :: phs
    allocate (phs_rambo_t :: phs)
  end subroutine phs_rambo_config_allocate_instance


end module phs_rambo
