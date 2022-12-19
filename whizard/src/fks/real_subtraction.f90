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

module real_subtraction

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use physics_defs
  use lorentz
  use flavors
  use phs_fks, only: real_kinematics_t, isr_kinematics_t
  use phs_fks, only: I_PLUS, I_MINUS
  use phs_fks, only: SQRTS_VAR, SQRTS_FIXED
  use phs_fks, only: phs_point_set_t

  use fks_regions
  use nlo_data

  implicit none
  private

  public :: soft_mismatch_t
  public :: coll_subtraction_t
  public :: real_subtraction_t
  public :: real_subtraction_compute_k_perp_fsr
  public :: real_subtraction_compute_k_perp_isr
  public :: real_partition_t
  public :: powheg_damping_simple_t
  public :: real_partition_fixed_order_t

  integer, parameter, public :: INTEGRATION = 0
  integer, parameter, public :: FIXED_ORDER_EVENTS = 1


  type :: soft_subtraction_t
    type(region_data_t), pointer :: reg_data => null ()
    real(default), dimension(:,:), allocatable :: momentum_matrix
    logical :: use_resonance_mappings = .false.
    type(vector4_t) :: p_soft = vector4_null
    logical :: use_internal_color_correlations = .true.
    logical :: use_internal_spin_correlations = .false.
    logical :: xi2_expanded = .true.
    integer :: factorization_mode = NO_FACTORIZATION
  contains
    procedure :: init => soft_subtraction_init
    procedure :: requires_boost => soft_subtraction_requires_boost
    procedure :: create_softvec_fsr => soft_subtraction_create_softvec_fsr
    procedure :: create_softvec_isr => soft_subtraction_create_softvec_isr
    procedure :: create_softvec_mismatch => &
       soft_subtraction_create_softvec_mismatch
    procedure :: compute => soft_subtraction_compute
    procedure :: evaluate_factorization_default => &
       soft_subtraction_evaluate_factorization_default
    procedure :: compute_momentum_matrix => &
         soft_subtraction_compute_momentum_matrix
    procedure :: evaluate_factorization_threshold => &
       soft_subtraction_evaluate_factorization_threshold
    procedure :: i_xi_ref => soft_subtraction_i_xi_ref
    procedure :: final => soft_subtraction_final
  end type soft_subtraction_t

  type :: soft_mismatch_t
    type(region_data_t), pointer :: reg_data => null ()
    real(default), dimension(:), allocatable :: sqme_born
    real(default), dimension(:,:,:), allocatable :: sqme_born_color_c
    real(default), dimension(:,:,:), allocatable :: sqme_born_charge_c
    type(real_kinematics_t), pointer :: real_kinematics => null ()
    type(soft_subtraction_t) :: sub_soft
  contains
    procedure :: init => soft_mismatch_init
    procedure :: evaluate => soft_mismatch_evaluate
    procedure :: compute => soft_mismatch_compute
    procedure :: final => soft_mismatch_final
  end type soft_mismatch_t

  type :: coll_subtraction_t
    integer :: n_in, n_alr
    logical :: use_resonance_mappings = .false.
    real(default) :: CA = 0, CF = 0, TR = 0
  contains
    procedure :: init => coll_subtraction_init
    procedure :: set_parameters => coll_subtraction_set_parameters
    procedure :: compute_fsr => coll_subtraction_compute_fsr
    procedure :: compute_isr => coll_subtraction_compute_isr
    procedure :: final => coll_subtraction_final
  end type coll_subtraction_t

  type :: real_subtraction_t
     type(nlo_settings_t), pointer :: settings => null ()
     type(region_data_t), pointer :: reg_data => null ()
     type(real_kinematics_t), pointer :: real_kinematics => null ()
     type(isr_kinematics_t), pointer :: isr_kinematics => null ()
     real(default), dimension(:,:), allocatable :: sqme_real_non_sub
     real(default), dimension(:), allocatable :: sqme_born
     real(default), dimension(:,:), allocatable :: sf_factors
     real(default), dimension(:), allocatable :: sqme_real_arr
     real(default), dimension(:,:,:), allocatable :: sqme_born_color_c
     real(default), dimension(:,:,:), allocatable :: sqme_born_charge_c
     real(default), dimension(:,:,:,:), allocatable :: sqme_born_spin_c
     type(soft_subtraction_t) :: sub_soft
     type(coll_subtraction_t) :: sub_coll
     logical, dimension(:), allocatable :: sc_required
     logical :: subtraction_deactivated = .false.
     integer :: purpose = INTEGRATION
     logical :: radiation_event = .true.
     logical :: subtraction_event = .false.
     integer, dimension(:), allocatable :: selected_alr
  contains
    procedure :: init => real_subtraction_init
    procedure :: set_real_kinematics => real_subtraction_set_real_kinematics
    procedure :: set_isr_kinematics => real_subtraction_set_isr_kinematics
    procedure :: get_i_res => real_subtraction_get_i_res
    procedure :: compute => real_subtraction_compute
    procedure :: evaluate_emitter_region => &
         real_subtraction_evaluate_emitter_region
    procedure :: sum_up_s_alpha => real_subtraction_sum_up_s_alpha
    procedure :: sum_up_s_alpha_soft => real_subtraction_sum_up_s_alpha_soft
    procedure :: evaluate_region_fsr => real_subtraction_evaluate_region_fsr
    procedure :: evaluate_region_isr => real_subtraction_evaluate_region_isr
    procedure :: evaluate_subtraction_terms_fsr => &
         real_subtraction_evaluate_subtraction_terms_fsr
    procedure :: evaluate_subtraction_terms_isr => &
        real_subtraction_evaluate_subtraction_terms_isr
    procedure :: get_phs_factor => real_subtraction_get_phs_factor
    procedure :: get_i_contributor => real_subtraction_get_i_contributor
    procedure :: compute_sub_soft => real_subtraction_compute_sub_soft
    procedure :: get_spin_correlation_term => &
         real_subtraction_get_spin_correlation_term
    procedure :: compute_sub_coll => real_subtraction_compute_sub_coll
    procedure :: compute_sub_coll_soft => real_subtraction_compute_sub_coll_soft
    procedure :: requires_spin_correlations => &
         real_subtraction_requires_spin_correlations
    procedure :: final => real_subtraction_final
  end type real_subtraction_t

  type, abstract :: real_partition_t
  contains
    procedure (real_partition_init), deferred :: init
    procedure (real_partition_write), deferred :: write
    procedure (real_partition_get_f), deferred :: get_f
  end type real_partition_t

  type, extends (real_partition_t) :: powheg_damping_simple_t
     real(default) :: h2 = 5._default
     integer :: emitter
  contains
    procedure :: get_f => powheg_damping_simple_get_f
    procedure :: init => powheg_damping_simple_init
    procedure :: write => powheg_damping_simple_write
  end type powheg_damping_simple_t

  type, extends (real_partition_t) :: real_partition_fixed_order_t
     real(default) :: scale
     type(ftuple_t), dimension(:), allocatable :: fks_pairs
  contains
    procedure :: init => real_partition_fixed_order_init
    procedure :: write => real_partition_fixed_order_write
    procedure :: get_f => real_partition_fixed_order_get_f
  end type real_partition_fixed_order_t


  abstract interface
     subroutine real_partition_init (partition, scale, reg_data)
       import
       class(real_partition_t), intent(out) :: partition
       real(default), intent(in) :: scale
       type(region_data_t), intent(in) :: reg_data
     end subroutine real_partition_init
  end interface

  abstract interface
     subroutine real_partition_write (partition, unit)
       import
       class(real_partition_t), intent(in) :: partition
       integer, intent(in), optional :: unit
     end subroutine real_partition_write
  end interface

  abstract interface
    function real_partition_get_f (partition, p) result (f)
       import
       real(default) :: f
       class(real_partition_t), intent(in) :: partition
       type(vector4_t), intent(in), dimension(:) :: p
    end function real_partition_get_f
  end interface


  interface
    module subroutine soft_subtraction_init (sub_soft, reg_data)
      class(soft_subtraction_t), intent(inout) :: sub_soft
      type(region_data_t), intent(in), target :: reg_data
    end subroutine soft_subtraction_init
    module function soft_subtraction_requires_boost &
         (sub_soft, sqrts) result (requires_boost)
      logical :: requires_boost
      class(soft_subtraction_t), intent(in) :: sub_soft
      real(default), intent(in) :: sqrts
    end function soft_subtraction_requires_boost
    module subroutine soft_subtraction_create_softvec_fsr &
         (sub_soft, p_born, y, phi, emitter, xi_ref_momentum)
      class(soft_subtraction_t), intent(inout) :: sub_soft
      type(vector4_t), intent(in), dimension(:) :: p_born
      real(default), intent(in) :: y, phi
      integer, intent(in) :: emitter
      type(vector4_t), intent(in) :: xi_ref_momentum
    end subroutine soft_subtraction_create_softvec_fsr
    module subroutine soft_subtraction_create_softvec_isr (sub_soft, y, phi)
      class(soft_subtraction_t), intent(inout) :: sub_soft
      real(default), intent(in) :: y, phi
    end subroutine soft_subtraction_create_softvec_isr
    module subroutine soft_subtraction_create_softvec_mismatch &
         (sub_soft, E, y, phi, p_em)
      class(soft_subtraction_t), intent(inout) :: sub_soft
      real(default), intent(in) :: E, phi, y
      type(vector4_t), intent(in) :: p_em
    end subroutine soft_subtraction_create_softvec_mismatch
    module function soft_subtraction_compute (sub_soft, p_born, &
         born_ij, y, q2, alpha_coupling, alr, emitter, i_res) result (sqme)
      real(default) :: sqme
      class(soft_subtraction_t), intent(inout) :: sub_soft
      type(vector4_t), intent(in), dimension(:) :: p_born
      real(default), intent(in), dimension(:,:) :: born_ij
      real(default), intent(in) :: y
      real(default), intent(in) :: q2, alpha_coupling
      integer, intent(in) :: alr, emitter, i_res
    end function soft_subtraction_compute
    module function soft_subtraction_evaluate_factorization_default &
         (sub_soft, p, born_ij) result (kb)
      real(default) :: kb
      class(soft_subtraction_t), intent(inout) :: sub_soft
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in), dimension(:,:) :: born_ij
    end function soft_subtraction_evaluate_factorization_default
    module subroutine soft_subtraction_compute_momentum_matrix &
         (sub_soft, p_born)
      class(soft_subtraction_t), intent(inout) :: sub_soft
      type(vector4_t), intent(in), dimension(:) :: p_born
    end subroutine soft_subtraction_compute_momentum_matrix
    module function soft_subtraction_evaluate_factorization_threshold &
         (sub_soft, leg, p_born, born_ij) result (kb)
      real(default) :: kb
      class(soft_subtraction_t), intent(inout) :: sub_soft
      integer, intent(in) :: leg
      type(vector4_t), intent(in), dimension(:) :: p_born
      real(default), intent(in), dimension(:,:) :: born_ij
    end function soft_subtraction_evaluate_factorization_threshold
    module function soft_subtraction_i_xi_ref &
         (sub_soft, alr, i_phs) result (i_xi_ref)
      integer :: i_xi_ref
      class(soft_subtraction_t), intent(in) :: sub_soft
      integer, intent(in) :: alr, i_phs
    end function soft_subtraction_i_xi_ref
    module subroutine soft_subtraction_final (sub_soft)
      class(soft_subtraction_t), intent(inout) :: sub_soft
    end subroutine soft_subtraction_final
    module subroutine soft_mismatch_init (soft_mismatch, reg_data, &
         real_kinematics, factorization_mode)
      class(soft_mismatch_t), intent(inout) :: soft_mismatch
      type(region_data_t), intent(in), target :: reg_data
      type(real_kinematics_t), intent(in), target :: real_kinematics
      integer, intent(in) :: factorization_mode
    end subroutine soft_mismatch_init
    module function soft_mismatch_evaluate &
         (soft_mismatch, alpha_s) result (sqme_mismatch)
      real(default) :: sqme_mismatch
      class(soft_mismatch_t), intent(inout) :: soft_mismatch
      real(default), intent(in) :: alpha_s
    end function soft_mismatch_evaluate
    module function soft_mismatch_compute &
         (soft_mismatch, alr, xi, y, p_em, p_res, p_soft, &
       sqme_born, sqme_soft, alpha_s, s) result (sqme_mismatch)
      real(default) :: sqme_mismatch
      class(soft_mismatch_t), intent(in) :: soft_mismatch
      integer, intent(in) :: alr
      real(default), intent(in) :: xi, y
      type(vector4_t), intent(in) :: p_em, p_res, p_soft
      real(default), intent(in) :: sqme_born, sqme_soft
      real(default), intent(in) :: alpha_s, s
    end function soft_mismatch_compute
    module subroutine soft_mismatch_final (soft_mismatch)
      class(soft_mismatch_t), intent(inout) :: soft_mismatch
    end subroutine soft_mismatch_final
    module subroutine coll_subtraction_init (coll_sub, n_alr, n_in)
      class(coll_subtraction_t), intent(inout) :: coll_sub
      integer, intent(in) :: n_alr, n_in
    end subroutine coll_subtraction_init
    module subroutine coll_subtraction_set_parameters (coll_sub, CA, CF, TR)
      class(coll_subtraction_t), intent(inout) :: coll_sub
      real(default), intent(in) :: CA, CF, TR
    end subroutine coll_subtraction_set_parameters
    module function coll_subtraction_compute_fsr (coll_sub, emitter, &
         flst, p_res, p_born, sqme_born, mom_times_sqme_spin_c, &
         xi, alpha_coupling, double_fsr) result (sqme)
      real(default) :: sqme
      class(coll_subtraction_t), intent(in) :: coll_sub
      integer, intent(in) :: emitter
      integer, dimension(:), intent(in) :: flst
      type(vector4_t), intent(in) :: p_res
      type(vector4_t), intent(in), dimension(:) :: p_born
      real(default), intent(in) :: sqme_born, mom_times_sqme_spin_c
      real(default), intent(in) :: xi, alpha_coupling
      logical, intent(in) :: double_fsr
    end function coll_subtraction_compute_fsr
    module function coll_subtraction_compute_isr &
       (coll_sub, emitter, flst, p_born, sqme_born, mom_times_sqme_spin_c, &
       xi, alpha_coupling, isr_mode) result (sqme)
      real(default) :: sqme
      class(coll_subtraction_t), intent(in) :: coll_sub
      integer, intent(in) :: emitter
      integer, dimension(:), intent(in) :: flst
      type(vector4_t), intent(in), dimension(:) :: p_born
      real(default), intent(in) :: sqme_born
      real(default), intent(in) :: mom_times_sqme_spin_c
      real(default), intent(in) :: xi, alpha_coupling
      integer, intent(in) :: isr_mode
    end function coll_subtraction_compute_isr
    module subroutine coll_subtraction_final (sub_coll)
      class(coll_subtraction_t), intent(inout) :: sub_coll
    end subroutine coll_subtraction_final
    module subroutine real_subtraction_init (rsub, reg_data, settings)
      class(real_subtraction_t), intent(inout), target :: rsub
      type(region_data_t), intent(in), target :: reg_data
      type(nlo_settings_t), intent(in), target :: settings
    end subroutine real_subtraction_init
    module subroutine real_subtraction_set_real_kinematics &
         (rsub, real_kinematics)
      class(real_subtraction_t), intent(inout) :: rsub
      type(real_kinematics_t), intent(in), target :: real_kinematics
    end subroutine real_subtraction_set_real_kinematics
    module subroutine real_subtraction_set_isr_kinematics (rsub, fractions)
      class(real_subtraction_t), intent(inout) :: rsub
      type(isr_kinematics_t), intent(in), target :: fractions
    end subroutine real_subtraction_set_isr_kinematics
    module function real_subtraction_get_i_res (rsub, alr) result (i_res)
      integer :: i_res
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr
    end function real_subtraction_get_i_res
    module subroutine real_subtraction_compute (rsub, emitter, &
         i_phs, alpha_s, alpha_qed, separate_alrs, sqme)
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: emitter, i_phs
      logical, intent(in) :: separate_alrs
      real(default), intent(inout), dimension(:) :: sqme
      real(default), intent(in) :: alpha_s, alpha_qed
    end subroutine real_subtraction_compute
    module function real_subtraction_evaluate_emitter_region (rsub, alr, &
         emitter, i_phs, i_res, alpha_coupling) result (sqme)
      real(default) :: sqme
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, emitter, i_phs, i_res
      real(default), intent(in) :: alpha_coupling
    end function real_subtraction_evaluate_emitter_region
    module function real_subtraction_sum_up_s_alpha &
         (rsub, alr, i_phs) result (sum_s_alpha)
      real(default) :: sum_s_alpha
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, i_phs
    end function real_subtraction_sum_up_s_alpha
    module function real_subtraction_sum_up_s_alpha_soft &
         (rsub, alr, i_phs) result (sum_s_alpha_soft)
      real(default) :: sum_s_alpha_soft
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, i_phs
    end function real_subtraction_sum_up_s_alpha_soft
    module function real_subtraction_evaluate_region_fsr (rsub, alr, &
         emitter, i_phs, i_res, alpha_coupling) result (sqme_tot)
      real(default) :: sqme_tot
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, emitter, i_phs, i_res
      real(default), intent(in) :: alpha_coupling
    end function real_subtraction_evaluate_region_fsr
    module function real_subtraction_evaluate_region_isr (rsub, alr, &
         emitter, i_phs, i_res, alpha_coupling) result (sqme_tot)
      real(default) :: sqme_tot
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, emitter, i_phs, i_res
      real(default), intent(in) :: alpha_coupling
    end function real_subtraction_evaluate_region_isr
    module subroutine real_subtraction_evaluate_subtraction_terms_fsr &
         (rsub, alr, emitter, i_phs, i_res, alpha_coupling, sqme_soft, &
          sqme_coll, sqme_cs)
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, emitter, i_phs, i_res
      real(default), intent(in) :: alpha_coupling
      real(default), intent(out) :: sqme_soft, sqme_coll, sqme_cs
    end subroutine real_subtraction_evaluate_subtraction_terms_fsr
    module subroutine real_subtraction_evaluate_subtraction_terms_isr (rsub, &
        alr, emitter, i_phs, i_res, alpha_coupling, sqme_soft, sqme_coll_plus, &
        sqme_coll_minus, sqme_cs_plus, sqme_cs_minus)
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, emitter, i_phs, i_res
      real(default), intent(in) :: alpha_coupling
      real(default), intent(out) :: sqme_soft
      real(default), intent(out) :: sqme_coll_plus, sqme_coll_minus
      real(default), intent(out) :: sqme_cs_plus, sqme_cs_minus
    end subroutine real_subtraction_evaluate_subtraction_terms_isr
    module function real_subtraction_get_phs_factor &
         (rsub, i_con) result (factor)
      real(default) :: factor
      class(real_subtraction_t), intent(in) :: rsub
      integer, intent(in) :: i_con
    end function real_subtraction_get_phs_factor
    module function real_subtraction_get_i_contributor &
         (rsub, alr) result (i_con)
      integer :: i_con
      class(real_subtraction_t), intent(in) :: rsub
      integer, intent(in) :: alr
    end function real_subtraction_get_i_contributor
    module function real_subtraction_compute_sub_soft (rsub, alr, emitter, &
         i_phs, i_res, alpha_coupling) result (sqme_soft)
      real(default) :: sqme_soft
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, emitter, i_phs, i_res
      real(default), intent(in) :: alpha_coupling
    end function real_subtraction_compute_sub_soft
    module function real_subtraction_get_spin_correlation_term &
         (rsub, alr, i_born, emitter) result (mom_times_sqme)
      real(default) :: mom_times_sqme
      class(real_subtraction_t), intent(in) :: rsub
      integer, intent(in) :: alr, i_born, emitter
    end function real_subtraction_get_spin_correlation_term
    module function real_subtraction_compute_k_perp_fsr &
         (p, phi) result (k_perp_fsr)
      real(default), dimension(0:3) :: k_perp_fsr
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: phi
    end function real_subtraction_compute_k_perp_fsr
    module function real_subtraction_compute_k_perp_isr &
         (p, phi) result (k_perp_isr)
      real(default), dimension(0:3) :: k_perp_isr
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: phi
    end function real_subtraction_compute_k_perp_isr
    module function real_subtraction_compute_sub_coll &
         (rsub, alr, em, i_phs, alpha_coupling) result (sqme_coll)
      real(default) :: sqme_coll
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, em, i_phs
      real(default), intent(in) :: alpha_coupling
    end function real_subtraction_compute_sub_coll
    module function real_subtraction_compute_sub_coll_soft &
         (rsub, alr, em, i_phs, alpha_coupling) result (sqme_cs)
      real(default) :: sqme_cs
      class(real_subtraction_t), intent(inout) :: rsub
      integer, intent(in) :: alr, em, i_phs
      real(default), intent(in) :: alpha_coupling
    end function real_subtraction_compute_sub_coll_soft
    module function real_subtraction_requires_spin_correlations &
         (rsub) result (val)
      logical :: val
      class(real_subtraction_t), intent(in) :: rsub
    end function real_subtraction_requires_spin_correlations
    module subroutine real_subtraction_final (rsub)
      class(real_subtraction_t), intent(inout) :: rsub
    end subroutine real_subtraction_final
    module function powheg_damping_simple_get_f (partition, p) result (f)
      real(default) :: f
      class(powheg_damping_simple_t), intent(in) :: partition
      type(vector4_t), intent(in), dimension(:) :: p
    end function powheg_damping_simple_get_f
    module subroutine powheg_damping_simple_init (partition, scale, reg_data)
      class(powheg_damping_simple_t), intent(out) :: partition
      real(default), intent(in) :: scale
      type(region_data_t), intent(in) :: reg_data
    end subroutine powheg_damping_simple_init
    module subroutine powheg_damping_simple_write (partition, unit)
      class(powheg_damping_simple_t), intent(in) :: partition
      integer, intent(in), optional :: unit
    end subroutine powheg_damping_simple_write
    module subroutine real_partition_fixed_order_init &
         (partition, scale, reg_data)
      class(real_partition_fixed_order_t), intent(out) :: partition
      real(default), intent(in) :: scale
      type(region_data_t), intent(in) :: reg_data
    end subroutine real_partition_fixed_order_init
    module subroutine real_partition_fixed_order_write (partition, unit)
      class(real_partition_fixed_order_t), intent(in) :: partition
      integer, intent(in), optional :: unit
    end subroutine real_partition_fixed_order_write
    module function real_partition_fixed_order_get_f (partition, p) result (f)
      real(default) :: f
      class(real_partition_fixed_order_t), intent(in) :: partition
      type(vector4_t), intent(in), dimension(:) :: p
    end function real_partition_fixed_order_get_f
  end interface

end module real_subtraction
