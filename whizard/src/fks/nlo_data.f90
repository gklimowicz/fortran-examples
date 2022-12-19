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

module nlo_data

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use constants, only: zero
  use lorentz
  use variables, only: var_list_t
  use physics_defs, only: NO_FACTORIZATION, FACTORIZATION_THRESHOLD

  implicit none
  private

  public :: fks_template_t
  public :: get_threshold_momenta
  public :: nlo_settings_t

  integer, parameter, public :: FKS_DEFAULT = 1
  integer, parameter, public :: FKS_RESONANCES = 2

  integer, dimension(2), parameter, public :: ASSOCIATED_LEG_PAIR = [1, 3]


  type :: fks_template_t
    logical :: subtraction_disabled = .false.
    integer :: mapping_type = FKS_DEFAULT
    logical :: count_kinematics = .false.
    real(default) :: fks_dij_exp1
    real(default) :: fks_dij_exp2
    real(default) :: xi_min
    real(default) :: y_max
    real(default) :: xi_cut, delta_o, delta_i
    type(string_t), dimension(:), allocatable :: excluded_resonances
    integer :: n_f
  contains
    procedure :: write => fks_template_write
    procedure :: set_parameters => fks_template_set_parameters
    procedure :: set_mapping_type => fks_template_set_mapping_type
    procedure :: set_counter => fks_template_set_counter
  end type fks_template_t

  type :: nlo_settings_t
     logical :: use_internal_color_correlations = .true.
     logical :: use_internal_spin_correlations = .false.
     logical :: use_resonance_mappings = .false.
     logical :: combined_integration = .false.
     logical :: fixed_order_nlo = .false.
     logical :: test_soft_limit = .false.
     logical :: test_coll_limit = .false.
     logical :: test_anti_coll_limit = .false.
     integer, dimension(:), allocatable :: selected_alr
     integer :: factorization_mode = NO_FACTORIZATION
     !!! Probably not the right place for this. Revisit after refactoring
     real(default) :: powheg_damping_scale = zero
     type(fks_template_t) :: fks_template
     type(string_t) :: virtual_selection
     logical :: virtual_resonance_aware_collinear = .true.
     logical :: use_born_scale = .false.
     logical :: cut_all_real_sqmes = .false.
     type(string_t) :: nlo_correction_type
     logical :: reuse_amplitudes_fks = .false.
  contains
    procedure :: init => nlo_settings_init
    procedure :: write => nlo_settings_write
  end type nlo_settings_t


  interface
    module subroutine fks_template_write (template, unit)
      class(fks_template_t), intent(in) :: template
      integer, intent(in), optional :: unit
    end subroutine fks_template_write
    module subroutine fks_template_set_parameters (template, &
         exp1, exp2, xi_min, y_max, xi_cut, delta_o, delta_i)
      class(fks_template_t), intent(inout) :: template
      real(default), intent(in) :: exp1, exp2
      real(default), intent(in) :: xi_min, y_max, &
           xi_cut, delta_o, delta_i
    end subroutine fks_template_set_parameters
    module subroutine fks_template_set_mapping_type (template, val)
      class(fks_template_t), intent(inout) :: template
      integer, intent(in) :: val
    end subroutine fks_template_set_mapping_type
    module subroutine fks_template_set_counter (template)
      class(fks_template_t), intent(inout) :: template
    end subroutine fks_template_set_counter
    module function get_threshold_momenta (p) result (p_thr)
      type(vector4_t), dimension(4) :: p_thr
      type(vector4_t), intent(in), dimension(:) :: p
    end function get_threshold_momenta
    module subroutine nlo_settings_init (nlo_settings, var_list, fks_template)
      class(nlo_settings_t), intent(inout) :: nlo_settings
      type(var_list_t), intent(in) :: var_list
      type(fks_template_t), intent(in), optional :: fks_template
    end subroutine nlo_settings_init
    module subroutine nlo_settings_write (nlo_settings, unit)
      class(nlo_settings_t), intent(in) :: nlo_settings
      integer, intent(in), optional :: unit
    end subroutine nlo_settings_write
  end interface

end module nlo_data

