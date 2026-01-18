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

module mappings

  use kinds, only: default
  use kinds, only: TC
  use iso_varying_string, string_t => varying_string
  use model_data
  use flavors

  implicit none
  private

  public :: mapping_defaults_t
  public :: mapping_defaults_md5sum
  public :: mapping_t
  public :: operator(==)

  integer, parameter :: &
       & EXTERNAL_PRT = -1, &
       & NO_MAPPING = 0, S_CHANNEL = 1, T_CHANNEL =  2, U_CHANNEL = 3, &
       & RADIATION = 4, COLLINEAR = 5, INFRARED = 6, &
       & STEP_MAPPING_E = 11, STEP_MAPPING_H = 12, &
       & ON_SHELL = 99

  type :: mapping_defaults_t
     real(default) :: energy_scale = 10
     real(default) :: invariant_mass_scale = 10
     real(default) :: momentum_transfer_scale = 10
     logical :: step_mapping = .true.
     logical :: step_mapping_exp = .true.
     logical :: enable_s_mapping = .false.
   contains
     procedure :: write => mapping_defaults_write
  end type mapping_defaults_t

  type :: mapping_t
     private
     integer :: type = NO_MAPPING
     integer(TC) :: bincode
     type(flavor_t) :: flv
     real(default) :: mass = 0
     real(default) :: width = 0
     logical :: a_unknown = .true.
     real(default) :: a1 = 0
     real(default) :: a2 = 0
     real(default) :: a3 = 0
     logical :: b_unknown = .true.
     real(default) :: b1 = 0
     real(default) :: b2 = 0
     real(default) :: b3 = 0
     logical :: variable_limits = .true.
   contains
     procedure :: write => mapping_write
     procedure :: init => mapping_init
     procedure :: set_parameters => mapping_set_parameters
     procedure :: set_step_mapping_parameters => &
          mapping_set_step_mapping_parameters
     procedure :: is_set => mapping_is_set
     procedure :: is_s_channel => mapping_is_s_channel
     procedure :: is_on_shell => mapping_is_on_shell
     procedure :: get_bincode => mapping_get_bincode
     procedure :: get_flv => mapping_get_flv
     procedure :: get_mass => mapping_get_mass
     procedure :: get_width => mapping_get_width
     procedure :: compute_msq_from_x => mapping_compute_msq_from_x
     procedure :: compute_x_from_msq => mapping_compute_x_from_msq
     procedure :: compute_ct_from_x => mapping_compute_ct_from_x
     procedure :: compute_x_from_ct => mapping_compute_x_from_ct
  end type mapping_t


  interface operator(==)
     module procedure mapping_equal
  end interface

  interface
    module subroutine mapping_defaults_write (object, unit)
      class(mapping_defaults_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine mapping_defaults_write
    module function mapping_defaults_md5sum &
         (mapping_defaults) result (md5sum_map)
      character(32) :: md5sum_map
      type(mapping_defaults_t), intent(in) :: mapping_defaults
    end function mapping_defaults_md5sum
    module subroutine mapping_write (map, unit, verbose)
      class(mapping_t), intent(in) :: map
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine mapping_write
    module subroutine mapping_init (mapping, bincode, type, f, model)
      class(mapping_t), intent(inout) :: mapping
      integer(TC), intent(in) :: bincode
      type(string_t), intent(in) :: type
      integer, intent(in), optional :: f
      class(model_data_t), intent(in), optional, target :: model
    end subroutine mapping_init
    module subroutine mapping_set_parameters &
         (map, mapping_defaults, variable_limits)
      class(mapping_t), intent(inout) :: map
      type(mapping_defaults_t), intent(in) :: mapping_defaults
      logical, intent(in) :: variable_limits
    end subroutine mapping_set_parameters
    module subroutine mapping_set_step_mapping_parameters (map, &
         mass, width, variable_limits)
      class(mapping_t), intent(inout) :: map
      real(default), intent(in) :: mass, width
      logical, intent(in) :: variable_limits
    end subroutine mapping_set_step_mapping_parameters
  module function mapping_is_set (mapping) result (flag)
    class(mapping_t), intent(in) :: mapping
    logical :: flag
  end function mapping_is_set
  module function mapping_is_s_channel (mapping) result (flag)
    class(mapping_t), intent(in) :: mapping
    logical :: flag
  end function mapping_is_s_channel
  module function mapping_is_on_shell (mapping) result (flag)
    class(mapping_t), intent(in) :: mapping
    logical :: flag
  end function mapping_is_on_shell
    module function mapping_get_bincode (mapping) result (bincode)
      class(mapping_t), intent(in) :: mapping
      integer(TC) :: bincode
    end function mapping_get_bincode
    module function mapping_get_flv (mapping) result (flv)
      class(mapping_t), intent(in) :: mapping
      type(flavor_t) :: flv
    end function mapping_get_flv
    module function mapping_get_mass (mapping) result (mass)
      class(mapping_t), intent(in) :: mapping
      real(default) :: mass
    end function mapping_get_mass
    module function mapping_get_width (mapping) result (width)
      class(mapping_t), intent(in) :: mapping
      real(default) :: width
    end function mapping_get_width
    module function mapping_equal (m1, m2) result (equal)
      type(mapping_t), intent(in) :: m1, m2
      logical :: equal
    end function mapping_equal
    module subroutine mapping_compute_msq_from_x &
         (map, s, msq_min, msq_max, msq, f, x)
      class(mapping_t), intent(inout) :: map
      real(default), intent(in) :: s, msq_min, msq_max
      real(default), intent(out) :: msq, f
      real(default), intent(in) :: x
    end subroutine mapping_compute_msq_from_x
    module subroutine mapping_compute_x_from_msq &
         (map, s, msq_min, msq_max, msq, f, x)
      class(mapping_t), intent(inout) :: map
      real(default), intent(in) :: s, msq_min, msq_max
      real(default), intent(in) :: msq
      real(default), intent(out) :: f, x
    end subroutine mapping_compute_x_from_msq
    module subroutine mapping_compute_ct_from_x (map, s, ct, st, f, x)
      class(mapping_t), intent(inout) :: map
      real(default), intent(in) :: s
      real(default), intent(out) :: ct, st, f
      real(default), intent(in) :: x
    end subroutine mapping_compute_ct_from_x
    module subroutine mapping_compute_x_from_ct (map, s, ct, f, x)
      class(mapping_t), intent(inout) :: map
      real(default), intent(in) :: s
      real(default), intent(in) :: ct
      real(default), intent(out) :: f, x
    end subroutine mapping_compute_x_from_ct
  end interface

end module mappings
