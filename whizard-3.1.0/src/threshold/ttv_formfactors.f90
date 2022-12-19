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

module ttv_formfactors

  use, intrinsic :: iso_fortran_env !NODEP!

  use kinds
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use constants
  use sm_physics
  use lorentz
  use interpolation
  use nr_tools
  use diagnostics

  implicit none
  private
  save

  public :: onshell_projection_t
  public :: settings_t
  public :: formfactor_t
  public :: width_t
  public :: threshold
  public :: threshold_t
  public :: GAM, GAM_M1S
  public :: AS_SOFT
  public :: AS_LL_SOFT
  public :: AS_USOFT
  public :: AS_HARD
  public :: SWITCHOFF_RESUMMED
  public :: TOPPIK_RESUMMED
  public :: m1s_to_mpole
  public :: phase_space_point_t
  public :: init_parameters
  public :: init_threshold_grids
  public :: v_matching
  public :: f_switch_off
  public :: alphas_notsohard
  public :: generate_on_shell_decay_threshold

  integer, parameter, public :: &
       MATCHED_EXPANDED_NOTSOHARD = -5, &
       MATCHED_NOTSOHARD = -4, &
       MATCHED_EXPANDED = - 3, &
       RESUMMED_SWITCHOFF = - 2, &
       MATCHED = -1, &
       RESUMMED = 1, &
       EXPANDED_HARD = 4, &
       EXPANDED_SOFT = 5, &
       EXPANDED_SOFT_SWITCHOFF = 6, &
       RESUMMED_ANALYTIC_LL = 7, &
       EXPANDED_NOTSOHARD = 8, &
       TREE = 9

  type :: onshell_projection_t
     logical :: production
     logical :: decay
     logical :: width
     logical :: boost_decay
  contains
     procedure :: debug_write => onshell_projection_debug_write
     procedure :: set_all => onshell_projection_set_all
     procedure :: active => onshell_projection_active
  end type onshell_projection_t

  type :: helicity_approximation_t
    logical :: simple = .false.
    logical :: extra = .false.
    logical :: ultra = .false.
  contains
  
  end type helicity_approximation_t

  type :: settings_t
     ! look what is set by initialized_parameters, bundle them in a class and rename to initialized
     logical :: initialized_parameters
     ! this belongs to init_threshold_phase_space_grid in phase_space_grid_t
     logical :: initialized_ps
     ! this belongs to the ff_grid_t, its usefulness is doubtful
     logical :: initialized_ff
     logical :: mpole_dynamic
     integer :: offshell_strategy
     logical :: factorized_computation
     logical :: interference
     logical :: only_interference_term
     logical :: nlo
     logical :: no_nlo_width_in_signal_propagators
     logical :: force_minus_one
     logical :: flip_relative_sign
     integer :: sel_hel_top = 0
     integer :: sel_hel_topbar = 0
     logical :: Z_disabled
     type(onshell_projection_t) :: onshell_projection
     type(helicity_approximation_t) :: helicity_approximation
  contains
     procedure :: setup_flags => settings_setup_flags
     procedure :: write => settings_write
     procedure :: use_nlo_width => settings_use_nlo_width
  end type settings_t

  type :: formfactor_t
     logical :: active
  contains
     procedure :: activate => formfactor_activate
     procedure :: disable => formfactor_disable
     procedure :: compute => formfactor_compute
  end type formfactor_t

  type :: width_t
    real(default) :: aem
    real(default) :: sw
    real(default) :: mw
    real(default) :: mb
    real(default) :: vtb
    real(default) :: gam_inv
  contains
     procedure :: init => width_init
     procedure :: compute => width_compute
  end type width_t

  type :: threshold_t
     type(settings_t) :: settings
     type(formfactor_t) :: formfactor
     type(width_t) :: width
  contains
   
  end type threshold_t

  type :: phase_space_point_t
    real(default) :: p2 = 0, k2 = 0, q2 = 0
    real(default) :: sqrts = 0, p = 0, p0 = 0
    real(default) :: mpole = 0, en = 0
    logical :: inside_grid = .false., onshell = .false.
  contains
    procedure :: init => phase_space_point_init_rel
    procedure :: init_nonrel => phase_space_point_init_nonrel
    procedure :: is_onshell => phase_space_point_is_onshell
    procedure :: write => phase_space_point_write
  end type phase_space_point_t


  type(threshold_t) :: threshold
  real(default) :: M1S, GAM, GAM_M1S
  integer :: NRQCD_ORDER
  real(default) :: MTPOLE = - one
  real(default) :: mtpole_init
  real(default) :: RESCALE_H, MU_HARD, AS_HARD
  real(default) :: AS_MZ, MASS_Z
  real(default) :: MU_USOFT, AS_USOFT

  real(default) :: RESCALE_F, MU_SOFT, AS_SOFT, AS_LL_SOFT, NUSTAR_FIXED
  logical :: NUSTAR_DYNAMIC, SWITCHOFF_RESUMMED, TOPPIK_RESUMMED
  real(default) :: B0
  real(default) :: B1

  real(default), dimension(2) :: aa2, aa3, aa4, aa5, aa8, aa0
  character(len=200) :: parameters_ref
  type(nr_spline_t) :: ff_p_spline
  real(default) :: v1, v2

  integer :: POINTS_SQ, POINTS_P, POINTS_P0, n_q
  real(default), dimension(:), allocatable :: sq_grid, p_grid, p0_grid, q_grid
  complex(default), dimension(:,:,:,:), allocatable :: ff_grid
  complex(single), dimension(:,:,:,:,:), allocatable :: Vmatrix

  real(default) :: sqrts_min, sqrts_max, sqrts_it


  interface
    module subroutine onshell_projection_debug_write (onshell_projection)
      class(onshell_projection_t), intent(in) :: onshell_projection
    end subroutine onshell_projection_debug_write
    pure module subroutine onshell_projection_set_all (onshell_projection, flag)
      class(onshell_projection_t), intent(inout) :: onshell_projection
      logical, intent(in) :: flag
    end subroutine onshell_projection_set_all
    pure module function onshell_projection_active (onshell_projection) result (active)
      logical :: active
      class(onshell_projection_t), intent(in) :: onshell_projection
    end function onshell_projection_active
    module subroutine settings_setup_flags (settings, ff_in, offshell_strategy_in, &
             top_helicity_selection)
      class(settings_t), intent(inout) :: settings
      integer, intent(in) :: ff_in, offshell_strategy_in, top_helicity_selection
    end subroutine settings_setup_flags
    module subroutine settings_write (settings, unit)
      class(settings_t), intent(in) :: settings
      integer, intent(in), optional :: unit
    end subroutine settings_write
    pure module function settings_use_nlo_width (settings, ff) result (nlo)
      logical :: nlo
      class(settings_t), intent(in) :: settings
      integer, intent(in) :: ff
    end function settings_use_nlo_width
    pure module subroutine formfactor_activate (formfactor)
      class(formfactor_t), intent(inout) :: formfactor
    end subroutine formfactor_activate
    pure module subroutine formfactor_disable (formfactor)
      class(formfactor_t), intent(inout) :: formfactor
    end subroutine formfactor_disable
    module function formfactor_compute (formfactor, ps, vec_type, FF_mode) result (FF)
      complex(default) :: FF
      class(formfactor_t), intent(in) :: formfactor
      type(phase_space_point_t), intent(in) :: ps
      integer, intent(in) :: vec_type, FF_mode
    end function formfactor_compute
    pure module subroutine width_init (width, aemi, sw, mw, mb, vtb, gam_inv)
      class(width_t), intent(inout) :: width
      real(default), intent(in) :: aemi, sw, mw, mb, vtb, gam_inv
    end subroutine width_init
    pure module function width_compute (width, top_mass, sqrts, initial) result (gamma)
      real(default) :: gamma
      class(width_t), intent(in) :: width
      real(default), intent(in) :: top_mass, sqrts
      logical, intent(in), optional :: initial
    end function width_compute
    pure module subroutine phase_space_point_init_rel (ps_point, p2, k2, q2, m)
      class(phase_space_point_t), intent(inout) :: ps_point
      real(default), intent(in) :: p2
      real(default), intent(in) :: k2
      real(default), intent(in) :: q2
      real(default), intent(in), optional :: m
    end subroutine phase_space_point_init_rel
    pure module subroutine phase_space_point_init_nonrel (ps_point, sqrts, p, p0, m)
      class(phase_space_point_t), intent(inout) :: ps_point
      real(default), intent(in) :: sqrts
      real(default), intent(in) :: p
      real(default), intent(in) :: p0
      real(default), intent(in), optional :: m
    end subroutine phase_space_point_init_nonrel
    pure module function phase_space_point_is_onshell (ps_point, m) result (flag)
      logical :: flag
      class(phase_space_point_t), intent(in) :: ps_point
      real(default), intent(in) :: m
    end function phase_space_point_is_onshell
    module subroutine phase_space_point_write (psp, unit)
      class(phase_space_point_t), intent(in) :: psp
      integer, intent(in), optional :: unit
    end subroutine phase_space_point_write
    module subroutine init_parameters (mpole_out, gam_out, m1s_in, Vtb, gam_inv, &
           aemi, sw, az, mz, mw, mb, h_in, f_in, nrqcd_order_in, ff_in, &
           offshell_strategy_in, v1_in, v2_in, scan_sqrts_min, &
           scan_sqrts_max, scan_sqrts_stepsize, mpole_fixed, top_helicity_selection)
      real(default), intent(out) :: mpole_out
      real(default), intent(out) :: gam_out
      real(default), intent(in) :: m1s_in
      real(default), intent(in) :: Vtb
      real(default), intent(in) :: gam_inv
      real(default), intent(in) :: aemi
      real(default), intent(in) :: sw
      real(default), intent(in) :: az
      real(default), intent(in) :: mz
      real(default), intent(in) :: mw
      real(default), intent(in) :: mb
      real(default), intent(in) :: h_in
      real(default), intent(in) :: f_in
      real(default), intent(in) :: nrqcd_order_in
      real(default), intent(in) :: ff_in
      real(default), intent(in) :: offshell_strategy_in
      real(default), intent(in) :: v1_in
      real(default), intent(in) :: v2_in
      real(default), intent(in) :: scan_sqrts_min
      real(default), intent(in) :: scan_sqrts_max
      real(default), intent(in) :: scan_sqrts_stepsize
      logical, intent(in) :: mpole_fixed
      real(default), intent(in) :: top_helicity_selection
    end subroutine init_parameters
    module subroutine init_threshold_grids (test)
      real(default), intent(in) :: test
    end subroutine init_threshold_grids
    pure module function v_matching (sqrts, gamma) result (v)
      real(default) :: v
      real(default), intent(in) :: sqrts, gamma
    end function v_matching
    pure module function f_switch_off (v) result (fval)
      real(default), intent(in) :: v
      real(default) :: fval
    end function f_switch_off
    pure module function alphas_notsohard (sqrts) result (a_soft)
      real(default) :: a_soft
      real(default), intent(in) :: sqrts
    end function alphas_notsohard
    pure module function m1s_to_mpole (sqrts) result (mpole)
      real(default), intent(in) :: sqrts
      real(default) :: mpole
    end function m1s_to_mpole
    pure module function int_to_char (i) result (c)
      integer, intent(in) :: i
      character(len=len(trim(int2fixed(i)))) :: c
    end function int_to_char
    pure module function real_to_char (r) result (c)
      real(default), intent(in) :: r
      character(len=len(trim(real2fixed(r)))) :: c
    end function real_to_char
    pure module function complex_to_char (z) result (c)
      complex(default), intent(in) :: z
      character(len=len(trim(real2fixed(real(z))))+len(trim(real2fixed(aimag(z))))+5) :: c
    end function complex_to_char
    pure module function logical_to_char (l) result (c)
      logical, intent(in) :: l
      character(len=1) :: c
    end function logical_to_char
    module subroutine generate_on_shell_decay_threshold (p_decay, p_top, p_decay_onshell)
      !!! Gluon must be on first position in this array
      type(vector4_t), intent(in), dimension(:) :: p_decay
      type(vector4_t), intent(inout) :: p_top
      type(vector4_t), intent(inout), dimension(:) :: p_decay_onshell
    end subroutine generate_on_shell_decay_threshold
  end interface

end module ttv_formfactors
