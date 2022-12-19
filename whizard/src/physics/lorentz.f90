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

module lorentz

  use kinds, only: default, double
  use constants, only: zero, one
  use c_particles

  implicit none
  private

  public :: vector3_t
  public :: vector3_write
  public :: vector3_null
  public :: vector3_canonical
  public :: vector3_moving
  public :: vector3_set_component
  public :: vector3_get_component
  public :: vector3_get_components
  public :: vector4_t
  public :: vector4_write
  public :: vector4_write_raw
  public :: vector4_read_raw
  public :: vector4_null
  public :: vector4_canonical
  public :: vector4_at_rest
  public :: vector4_moving
  public :: vector4_set_component
  public :: vector4_get_component
  public :: vector4_get_components
  public :: vector4_invert_direction
  public :: assignment (=)
  public :: vector4
  public :: vector4_to_c_prt
  public :: lorentz_transformation_t
  public :: lorentz_transformation_write
  public :: lorentz_transformation_get_components
  public :: identity
  public :: space_reflection
  public :: compute_resonance_mass
  public :: get_resonance_momentum
  public :: create_two_particle_decay
  public :: create_three_particle_decay
  public :: evaluate_one_to_two_splitting_special
  public :: generate_on_shell_decay
  public :: vector_set_reshuffle
  public :: vector_set_is_cms
  public :: vector4_write_set
  public :: vector4_check_momentum_conservation
  public :: spinor_product

  public :: operator(==), operator(/=)
  public :: operator(+), operator(-)
  public :: operator(*), operator(/)
  public :: operator(**)

  public :: cross_product
  public :: sum
  public :: direction
  public :: space_part
  public :: azimuthal_angle
  public :: azimuthal_angle_deg
  public :: azimuthal_distance
  public :: azimuthal_distance_deg
  public :: polar_angle
  public :: polar_angle_ct
  public :: polar_angle_deg
  public :: enclosed_angle
  public :: enclosed_angle_ct
  public :: enclosed_angle_deg
  public :: enclosed_angle_rest_frame
  public :: enclosed_angle_ct_rest_frame
  public :: enclosed_angle_deg_rest_frame
  public :: transverse_part
  public :: longitudinal_part
  public :: space_part_norm
  public :: energy
  public :: invariant_mass
  public :: invariant_mass_squared
  public :: transverse_mass
  public :: rapidity
  public :: pseudorapidity
  public :: rapidity_distance
  public :: pseudorapidity_distance
  public :: eta_phi_distance
  public :: inverse
  public :: create_orthogonal
  public :: create_unit_vector
  public :: normalize
  public :: boost
  public :: rotation
  public :: rotation_to_2nd
  public :: transformation
  public :: LT_compose_r3_r2_b3
  public :: LT_compose_r2_r3_b3
  public :: axis_from_p_r3_r2_b3, axis_from_p_b3
  public :: lambda
  public :: colliding_momenta
  public :: pacify

  type :: vector3_t
     real(default), dimension(3) :: p
  end type vector3_t

  type :: vector4_t
     real(default), dimension(0:3) :: p = &
        [zero, zero, zero, zero]
  contains
    procedure :: write => vector4_write
    procedure :: to_pythia6 => vector4_to_pythia6
  end type vector4_t
  type :: lorentz_transformation_t
     private
     real(default), dimension(0:3, 0:3) :: L
   contains
     procedure :: write => lorentz_transformation_write
  end type lorentz_transformation_t


  type(vector3_t), parameter :: vector3_null = &
       vector3_t ([ zero, zero, zero ])

  type(vector4_t), parameter :: vector4_null = &
       vector4_t ([ zero, zero, zero, zero ])

  integer, dimension(3,3), parameter :: delta_three = &
       & reshape( source = [ 1,0,0, 0,1,0, 0,0,1 ], &
       &          shape  = [3,3] )
  integer, dimension(3,3,3), parameter :: epsilon_three = &
       & reshape( source = [ 0, 0,0,  0,0,-1,   0,1,0, &
       &                     0, 0,1,  0,0, 0,  -1,0,0, &
       &                     0,-1,0,  1,0, 0,   0,0,0 ],&
       &          shape = [3,3,3] )
  type(lorentz_transformation_t), parameter :: &
       & identity = &
       & lorentz_transformation_t ( &
       & reshape( source = [ one, zero, zero, zero, &
       &                     zero, one, zero, zero, &
       &                     zero, zero, one, zero, &
       &                     zero, zero, zero, one ],&
       &          shape = [4,4] ) )
  type(lorentz_transformation_t), parameter :: &
       & space_reflection = &
       & lorentz_transformation_t ( &
       & reshape( source = [ one, zero, zero, zero, &
       &                     zero,-one, zero, zero, &
       &                     zero, zero,-one, zero, &
       &                     zero, zero, zero,-one ],&
       &          shape = [4,4] ) )

  interface vector3_moving
     module procedure vector3_moving_canonical
     module procedure vector3_moving_generic
  end interface
  interface operator(==)
     module procedure vector3_eq
  end interface
  interface operator(/=)
     module procedure vector3_neq
  end interface
  interface operator(+)
     module procedure add_vector3
  end interface
  interface operator(-)
     module procedure sub_vector3
  end interface
  interface operator(*)
     module procedure prod_integer_vector3, prod_vector3_integer
     module procedure prod_real_vector3, prod_vector3_real
  end interface
  interface operator(/)
     module procedure div_vector3_real, div_vector3_integer
  end interface
  interface operator(*)
     module procedure prod_vector3
  end interface
  interface cross_product
     module procedure vector3_cross_product
  end interface
  interface operator(**)
     module procedure power_vector3
  end interface
  interface operator(-)
     module procedure negate_vector3
  end interface
  interface sum
     module procedure sum_vector3
  end interface
  interface direction
     module procedure vector3_get_direction
  end interface
  interface vector4_moving
     module procedure vector4_moving_canonical
     module procedure vector4_moving_generic
  end interface
  interface operator(==)
     module procedure vector4_eq
  end interface
  interface operator(/=)
     module procedure vector4_neq
  end interface
  interface operator(+)
     module procedure add_vector4
  end interface
  interface operator(-)
     module procedure sub_vector4
  end interface
  interface operator(*)
     module procedure prod_real_vector4, prod_vector4_real
     module procedure prod_integer_vector4, prod_vector4_integer
  end interface
  interface operator(/)
     module procedure div_vector4_real
     module procedure div_vector4_integer
  end interface
  interface operator(*)
     module procedure prod_vector4
  end interface
  interface operator(**)
     module procedure power_vector4
  end interface
  interface operator(-)
     module procedure negate_vector4
  end interface
  interface sum
     module procedure sum_vector4, sum_vector4_mask
  end interface
  interface space_part
     module procedure vector4_get_space_part
  end interface
  interface direction
     module procedure vector4_get_direction
  end interface
  interface assignment (=)
     module procedure array_from_vector4_1, array_from_vector4_2, &
            array_from_vector3_1, array_from_vector3_2, &
            vector4_from_array, vector3_from_array
  end interface
  interface assignment (=)
     module procedure vector4_from_c_prt, c_prt_from_vector4
  end interface
  interface azimuthal_angle
     module procedure vector3_azimuthal_angle
     module procedure vector4_azimuthal_angle
  end interface
  interface azimuthal_angle_deg
     module procedure vector3_azimuthal_angle_deg
     module procedure vector4_azimuthal_angle_deg
  end interface
  interface azimuthal_distance
     module procedure vector3_azimuthal_distance
     module procedure vector4_azimuthal_distance
  end interface
  interface azimuthal_distance_deg
     module procedure vector3_azimuthal_distance_deg
     module procedure vector4_azimuthal_distance_deg
  end interface
  interface polar_angle
     module procedure polar_angle_vector3
     module procedure polar_angle_vector4
  end interface
  interface polar_angle_ct
     module procedure polar_angle_ct_vector3
     module procedure polar_angle_ct_vector4
  end interface
  interface polar_angle_deg
     module procedure polar_angle_deg_vector3
     module procedure polar_angle_deg_vector4
  end interface
  interface enclosed_angle
     module procedure enclosed_angle_vector3
     module procedure enclosed_angle_vector4
  end interface
  interface enclosed_angle_ct
     module procedure enclosed_angle_ct_vector3
     module procedure enclosed_angle_ct_vector4
  end interface
  interface enclosed_angle_deg
     module procedure enclosed_angle_deg_vector3
     module procedure enclosed_angle_deg_vector4
  end interface
  interface enclosed_angle_rest_frame
     module procedure enclosed_angle_rest_frame_vector4
  end interface
  interface enclosed_angle_ct_rest_frame
     module procedure enclosed_angle_ct_rest_frame_vector4
  end interface
  interface enclosed_angle_deg_rest_frame
     module procedure enclosed_angle_deg_rest_frame_vector4
  end interface
  interface transverse_part
     module procedure transverse_part_vector4_beam_axis
     module procedure transverse_part_vector4_vector4
  end interface
  interface longitudinal_part
     module procedure longitudinal_part_vector4
  end interface
  interface space_part_norm
     module procedure space_part_norm_vector4
  end interface
  interface energy
     module procedure energy_vector4
     module procedure energy_vector3
     module procedure energy_real
  end interface
  interface invariant_mass
     module procedure invariant_mass_vector4
  end interface
  interface invariant_mass_squared
     module procedure invariant_mass_squared_vector4
  end interface
  interface transverse_mass
     module procedure transverse_mass_vector4
  end interface
  interface rapidity
     module procedure rapidity_vector4
  end interface
  interface pseudorapidity
     module procedure pseudorapidity_vector4
  end interface
  interface rapidity_distance
     module procedure rapidity_distance_vector4
  end interface
  interface pseudorapidity_distance
     module procedure pseudorapidity_distance_vector4
  end interface
  interface eta_phi_distance
     module procedure eta_phi_distance_vector4
  end interface
  interface inverse
     module procedure lorentz_transformation_inverse
  end interface
  abstract interface
     subroutine evaluate_one_to_two_splitting_special (p_origin, &
          p1_in, p2_in, p1_out, p2_out, msq_in, jac)
       import
       type(vector4_t), intent(in) :: p_origin
       type(vector4_t), intent(in) :: p1_in, p2_in
       type(vector4_t), intent(inout) :: p1_out, p2_out
       real(default), intent(in), optional :: msq_in
       real(default), intent(inout), optional :: jac
     end subroutine evaluate_one_to_two_splitting_special
  end interface

  interface boost
     module procedure boost_from_rest_frame
     module procedure boost_from_rest_frame_vector3
     module procedure boost_generic
     module procedure boost_canonical
  end interface
  interface rotation
     module procedure rotation_generic
     module procedure rotation_canonical
     module procedure rotation_generic_cs
     module procedure rotation_canonical_cs
  end interface
  interface rotation_to_2nd
     module procedure rotation_to_2nd_generic
     module procedure rotation_to_2nd_canonical
  end interface
  interface transformation
     module procedure transformation_rec_generic
     module procedure transformation_rec_canonical
  end interface
  interface operator(*)
     module procedure prod_LT_vector4
     module procedure prod_LT_LT
     module procedure prod_vector4_LT
  end interface
  interface pacify
     module procedure pacify_vector3
     module procedure pacify_vector4
     module procedure pacify_LT
  end interface pacify


  interface
    module subroutine vector3_write (p, unit, testflag)
      type(vector3_t), intent(in) :: p
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine vector3_write
    elemental module function vector3_canonical (k) result (p)
      type(vector3_t) :: p
      integer, intent(in) :: k
    end function vector3_canonical
    elemental module function vector3_moving_canonical (p, k) result(q)
      type(vector3_t) :: q
      real(default), intent(in) :: p
      integer, intent(in) :: k
    end function vector3_moving_canonical
    pure module function vector3_moving_generic (p) result(q)
      real(default), dimension(3), intent(in) :: p
      type(vector3_t) :: q
    end function vector3_moving_generic
    elemental module function vector3_eq (p, q) result (r)
      logical :: r
      type(vector3_t), intent(in) :: p,q
    end function vector3_eq
    elemental module function vector3_neq (p, q) result (r)
      logical :: r
      type(vector3_t), intent(in) :: p,q
    end function vector3_neq
    elemental module function add_vector3 (p, q) result (r)
      type(vector3_t) :: r
      type(vector3_t), intent(in) :: p,q
    end function add_vector3
    elemental module function sub_vector3 (p, q) result (r)
      type(vector3_t) :: r
      type(vector3_t), intent(in) :: p,q
    end function sub_vector3
    elemental module function prod_real_vector3 (s, p) result (q)
      type(vector3_t) :: q
      real(default), intent(in) :: s
      type(vector3_t), intent(in) :: p
    end function prod_real_vector3
    elemental module function prod_vector3_real (p, s) result (q)
      type(vector3_t) :: q
      real(default), intent(in) :: s
      type(vector3_t), intent(in) :: p
    end function prod_vector3_real
    elemental module function div_vector3_real (p, s) result (q)
      type(vector3_t) :: q
      real(default), intent(in) :: s
      type(vector3_t), intent(in) :: p
    end function div_vector3_real
    elemental module function prod_integer_vector3 (s, p) result (q)
      type(vector3_t) :: q
      integer, intent(in) :: s
      type(vector3_t), intent(in) :: p
    end function prod_integer_vector3
    elemental module function prod_vector3_integer (p, s) result (q)
      type(vector3_t) :: q
      integer, intent(in) :: s
      type(vector3_t), intent(in) :: p
    end function prod_vector3_integer
    elemental module function div_vector3_integer (p, s) result (q)
      type(vector3_t) :: q
      integer, intent(in) :: s
      type(vector3_t), intent(in) :: p
    end function div_vector3_integer
    elemental module function prod_vector3 (p, q) result (s)
      real(default) :: s
      type(vector3_t), intent(in) :: p,q
    end function prod_vector3
    elemental module function vector3_cross_product (p, q) result (r)
      type(vector3_t) :: r
      type(vector3_t), intent(in) :: p,q
    end function vector3_cross_product
    elemental module function power_vector3 (p, e) result (s)
      real(default) :: s
      type(vector3_t), intent(in) :: p
      integer, intent(in) :: e
    end function power_vector3
    elemental module function negate_vector3 (p) result (q)
      type(vector3_t) :: q
      type(vector3_t), intent(in) :: p
    end function negate_vector3
    module subroutine vector3_set_component (p, i, value)
      type(vector3_t), intent(inout) :: p
      integer, intent(in) :: i
      real(default), intent(in) :: value
    end subroutine vector3_set_component
    pure module function sum_vector3 (p) result (q)
      type(vector3_t) :: q
      type(vector3_t), dimension(:), intent(in) :: p
    end function sum_vector3
    elemental module function vector3_get_component (p, k) result (c)
      type(vector3_t), intent(in) :: p
      integer, intent(in) :: k
      real(default) :: c
    end function vector3_get_component
    pure module function vector3_get_components (p) result (a)
      type(vector3_t), intent(in) :: p
      real(default), dimension(3) :: a
    end function vector3_get_components
    elemental module function vector3_get_direction (p) result (q)
      type(vector3_t) :: q
      type(vector3_t), intent(in) :: p
    end function vector3_get_direction
    module subroutine vector4_write &
           (p, unit, show_mass, testflag, compressed, ultra)
      class(vector4_t), intent(in) :: p
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_mass, testflag, compressed, ultra
    end subroutine vector4_write
    module subroutine vector4_write_raw (p, u)
      type(vector4_t), intent(in) :: p
      integer, intent(in) :: u
    end subroutine vector4_write_raw
    module subroutine vector4_read_raw (p, u, iostat)
      type(vector4_t), intent(out) :: p
      integer, intent(in) :: u
      integer, intent(out), optional :: iostat
    end subroutine vector4_read_raw
    elemental module function vector4_canonical (k) result (p)
      type(vector4_t) :: p
      integer, intent(in) :: k
    end function vector4_canonical
    elemental module function vector4_at_rest (m) result (p)
      type(vector4_t) :: p
      real(default), intent(in) :: m
    end function vector4_at_rest
    elemental module function vector4_moving_canonical (E, p, k) result (q)
      type(vector4_t) :: q
      real(default), intent(in) :: E, p
      integer, intent(in) :: k
    end function vector4_moving_canonical
    elemental module function vector4_moving_generic (E, p) result (q)
      type(vector4_t) :: q
      real(default), intent(in) :: E
      type(vector3_t), intent(in) :: p
    end function vector4_moving_generic
    elemental module function vector4_eq (p, q) result (r)
      logical :: r
      type(vector4_t), intent(in) :: p,q
    end function vector4_eq
    elemental module function vector4_neq (p, q) result (r)
      logical :: r
      type(vector4_t), intent(in) :: p,q
    end function vector4_neq
    elemental module function add_vector4 (p,q) result (r)
      type(vector4_t) :: r
      type(vector4_t), intent(in) :: p,q
    end function add_vector4
    elemental module function sub_vector4 (p,q) result (r)
      type(vector4_t) :: r
      type(vector4_t), intent(in) :: p,q
    end function sub_vector4
    elemental module function prod_real_vector4 (s, p) result (q)
      type(vector4_t) :: q
      real(default), intent(in) :: s
      type(vector4_t), intent(in) :: p
    end function prod_real_vector4
    elemental module function prod_vector4_real (p, s) result (q)
      type(vector4_t) :: q
      real(default), intent(in) :: s
      type(vector4_t), intent(in) :: p
    end function prod_vector4_real
    elemental module function div_vector4_real (p, s) result (q)
      type(vector4_t) :: q
      real(default), intent(in) :: s
      type(vector4_t), intent(in) :: p
    end function div_vector4_real
    elemental module function prod_integer_vector4 (s, p) result (q)
      type(vector4_t) :: q
      integer, intent(in) :: s
      type(vector4_t), intent(in) :: p
    end function prod_integer_vector4
    elemental module function prod_vector4_integer (p, s) result (q)
      type(vector4_t) :: q
      integer, intent(in) :: s
      type(vector4_t), intent(in) :: p
    end function prod_vector4_integer
    elemental module function div_vector4_integer (p, s) result (q)
      type(vector4_t) :: q
      integer, intent(in) :: s
      type(vector4_t), intent(in) :: p
    end function div_vector4_integer
    elemental module function prod_vector4 (p, q) result (s)
      real(default) :: s
      type(vector4_t), intent(in) :: p,q
    end function prod_vector4
    elemental module function power_vector4 (p, e) result (s)
      real(default) :: s
      type(vector4_t), intent(in) :: p
      integer, intent(in) :: e
    end function power_vector4
    elemental module function negate_vector4 (p) result (q)
      type(vector4_t) :: q
      type(vector4_t), intent(in) :: p
    end function negate_vector4
    pure module function sum_vector4 (p) result (q)
      type(vector4_t) :: q
      type(vector4_t), dimension(:), intent(in) :: p
    end function sum_vector4
    pure module function sum_vector4_mask (p, mask) result (q)
      type(vector4_t) :: q
      type(vector4_t), dimension(:), intent(in) :: p
      logical, dimension(:), intent(in) :: mask
    end function sum_vector4_mask
    module subroutine vector4_set_component (p, k, c)
      type(vector4_t), intent(inout) :: p
      integer, intent(in) :: k
      real(default), intent(in) :: c
    end subroutine vector4_set_component
    elemental module function vector4_get_component (p, k) result (c)
      real(default) :: c
      type(vector4_t), intent(in) :: p
      integer, intent(in) :: k
    end function vector4_get_component
    pure module function vector4_get_components (p) result (a)
      real(default), dimension(0:3) :: a
      type(vector4_t), intent(in) :: p
    end function vector4_get_components
    elemental module function vector4_get_space_part (p) result (q)
      type(vector3_t) :: q
      type(vector4_t), intent(in) :: p
    end function vector4_get_space_part
    elemental module function vector4_get_direction (p) result (q)
      type(vector3_t) :: q
      type(vector4_t), intent(in) :: p
    end function vector4_get_direction
    elemental module subroutine vector4_invert_direction (p)
      type(vector4_t), intent(inout) :: p
    end subroutine vector4_invert_direction
    pure module subroutine array_from_vector4_1 (a, p)
      real(default), dimension(:), intent(out) :: a
      type(vector4_t), intent(in) :: p
    end subroutine array_from_vector4_1
    pure module subroutine array_from_vector4_2 (a, p)
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), dimension(:,:), intent(out) :: a
    end subroutine array_from_vector4_2
    pure module subroutine array_from_vector3_1 (a, p)
      real(default), dimension(:), intent(out) :: a
      type(vector3_t), intent(in) :: p
    end subroutine array_from_vector3_1
    pure module subroutine array_from_vector3_2 (a, p)
      type(vector3_t), dimension(:), intent(in) :: p
      real(default), dimension(:,:), intent(out) :: a
    end subroutine array_from_vector3_2
    pure module subroutine vector4_from_array (p, a)
      type(vector4_t), intent(out) :: p
      real(default), dimension(:), intent(in) :: a
    end subroutine vector4_from_array
    pure module subroutine vector3_from_array (p, a)
      type(vector3_t), intent(out) :: p
      real(default), dimension(:), intent(in) :: a
    end subroutine vector3_from_array
    pure module function vector4 (a) result (p)
      type(vector4_t) :: p
      real(default), intent(in), dimension(4) :: a
    end function vector4
    pure module function vector4_to_pythia6 (vector4, m) result (p)
      real(double), dimension(1:5) :: p
      class(vector4_t), intent(in) :: vector4
      real(default), intent(in), optional :: m
    end function vector4_to_pythia6
    pure module subroutine vector4_from_c_prt (p, c_prt)
      type(vector4_t), intent(out) :: p
      type(c_prt_t), intent(in) :: c_prt
    end subroutine vector4_from_c_prt
    pure module subroutine c_prt_from_vector4 (c_prt, p)
      type(c_prt_t), intent(out) :: c_prt
      type(vector4_t), intent(in) :: p
    end subroutine c_prt_from_vector4
    elemental module function vector4_to_c_prt (p, p2) result (c_prt)
      type(c_prt_t) :: c_prt
      type(vector4_t), intent(in) :: p
      real(default), intent(in), optional :: p2
    end function vector4_to_c_prt
    elemental module function vector3_azimuthal_angle (p) result (phi)
      real(default) :: phi
      type(vector3_t), intent(in) :: p
    end function vector3_azimuthal_angle
    elemental module function vector4_azimuthal_angle (p) result (phi)
      real(default) :: phi
      type(vector4_t), intent(in) :: p
    end function vector4_azimuthal_angle
    elemental module function vector3_azimuthal_angle_deg (p) result (phi)
      real(default) :: phi
      type(vector3_t), intent(in) :: p
    end function vector3_azimuthal_angle_deg
    elemental module function vector4_azimuthal_angle_deg (p) result (phi)
      real(default) :: phi
      type(vector4_t), intent(in) :: p
    end function vector4_azimuthal_angle_deg
    elemental module function vector3_azimuthal_distance (p, q) result (dphi)
      real(default) :: dphi
      type(vector3_t), intent(in) :: p,q
    end function vector3_azimuthal_distance
    elemental module function vector4_azimuthal_distance (p, q) result (dphi)
      real(default) :: dphi
      type(vector4_t), intent(in) :: p,q
    end function vector4_azimuthal_distance
    elemental module function vector3_azimuthal_distance_deg (p, q) result (dphi)
      real(default) :: dphi
      type(vector3_t), intent(in) :: p,q
    end function vector3_azimuthal_distance_deg
    elemental module function vector4_azimuthal_distance_deg (p, q) result (dphi)
      real(default) :: dphi
      type(vector4_t), intent(in) :: p,q
    end function vector4_azimuthal_distance_deg
    elemental module function polar_angle_vector3 (p) result (theta)
      real(default) :: theta
      type(vector3_t), intent(in) :: p
    end function polar_angle_vector3
    elemental module function polar_angle_vector4 (p) result (theta)
      real(default) :: theta
      type(vector4_t), intent(in) :: p
    end function polar_angle_vector4
    elemental module function polar_angle_ct_vector3 (p) result (ct)
      real(default) :: ct
      type(vector3_t), intent(in) :: p
    end function polar_angle_ct_vector3
    elemental module function polar_angle_ct_vector4 (p) result (ct)
      real(default) :: ct
      type(vector4_t), intent(in) :: p
    end function polar_angle_ct_vector4
    elemental module function polar_angle_deg_vector3 (p) result (theta)
      real(default) :: theta
      type(vector3_t), intent(in) :: p
    end function polar_angle_deg_vector3
    elemental module function polar_angle_deg_vector4 (p) result (theta)
      real(default) :: theta
      type(vector4_t), intent(in) :: p
    end function polar_angle_deg_vector4
    elemental module function enclosed_angle_vector3 (p, q) result (theta)
      real(default) :: theta
      type(vector3_t), intent(in) :: p, q
    end function enclosed_angle_vector3
    elemental module function enclosed_angle_vector4 (p, q) result (theta)
      real(default) :: theta
      type(vector4_t), intent(in) :: p, q
    end function enclosed_angle_vector4
    elemental module function enclosed_angle_ct_vector3 (p, q) result (ct)
      real(default) :: ct
      type(vector3_t), intent(in) :: p, q
    end function enclosed_angle_ct_vector3
    elemental module function enclosed_angle_ct_vector4 (p, q) result (ct)
      real(default) :: ct
      type(vector4_t), intent(in) :: p, q
    end function enclosed_angle_ct_vector4
    elemental module function enclosed_angle_deg_vector3 (p, q) result (theta)
      real(default) :: theta
      type(vector3_t), intent(in) :: p, q
    end function enclosed_angle_deg_vector3
    elemental module function enclosed_angle_deg_vector4 (p, q) result (theta)
      real(default) :: theta
      type(vector4_t), intent(in) :: p, q
    end function enclosed_angle_deg_vector4
    elemental module function enclosed_angle_rest_frame_vector4 (p, q) result (theta)
      type(vector4_t), intent(in) :: p, q
      real(default) :: theta
    end function enclosed_angle_rest_frame_vector4
    elemental module function enclosed_angle_ct_rest_frame_vector4 (p, q) result (ct)
      type(vector4_t), intent(in) :: p, q
      real(default) :: ct
    end function enclosed_angle_ct_rest_frame_vector4
    elemental module function enclosed_angle_deg_rest_frame_vector4 (p, q) &
         result (theta)
      type(vector4_t), intent(in) :: p, q
      real(default) :: theta
    end function enclosed_angle_deg_rest_frame_vector4
    elemental module function transverse_part_vector4_beam_axis (p) result (pT)
      real(default) :: pT
      type(vector4_t), intent(in) :: p
    end function transverse_part_vector4_beam_axis
    elemental module function transverse_part_vector4_vector4 (p1, p2) result (pT)
      real(default) :: pT
      type(vector4_t), intent(in) :: p1, p2
    end function transverse_part_vector4_vector4
    elemental module function longitudinal_part_vector4 (p) result (pL)
      real(default) :: pL
      type(vector4_t), intent(in) :: p
    end function longitudinal_part_vector4
    elemental module function space_part_norm_vector4 (p) result (p3)
      real(default) :: p3
      type(vector4_t), intent(in) :: p
    end function space_part_norm_vector4
    elemental module function energy_vector4 (p) result (E)
      real(default) :: E
      type(vector4_t), intent(in) :: p
    end function energy_vector4
    elemental module function energy_vector3 (p, mass) result (E)
      real(default) :: E
      type(vector3_t), intent(in) :: p
      real(default), intent(in), optional :: mass
    end function energy_vector3
    elemental module function energy_real (p, mass) result (E)
      real(default) :: E
      real(default), intent(in) :: p
      real(default), intent(in), optional :: mass
    end function energy_real
    elemental module function invariant_mass_vector4 (p) result (m)
      real(default) :: m
      type(vector4_t), intent(in) :: p
    end function invariant_mass_vector4
    elemental module function invariant_mass_squared_vector4 (p) result (msq)
      real(default) :: msq
      type(vector4_t), intent(in) :: p
    end function invariant_mass_squared_vector4
    elemental module function transverse_mass_vector4 (p) result (m)
      real(default) :: m
      type(vector4_t), intent(in) :: p
    end function transverse_mass_vector4
    elemental module function rapidity_vector4 (p) result (y)
      real(default) :: y
      type(vector4_t), intent(in) :: p
    end function rapidity_vector4
    elemental module function pseudorapidity_vector4 (p) result (eta)
      real(default) :: eta
      type(vector4_t), intent(in) :: p
    end function pseudorapidity_vector4
    elemental module function rapidity_distance_vector4 (p, q) result (dy)
      type(vector4_t), intent(in) :: p, q
      real(default) :: dy
    end function rapidity_distance_vector4
    elemental module function pseudorapidity_distance_vector4 (p, q) result (deta)
      real(default) :: deta
      type(vector4_t), intent(in) :: p, q
    end function pseudorapidity_distance_vector4
    elemental module function eta_phi_distance_vector4 (p, q) result (dr)
      type(vector4_t), intent(in) :: p, q
      real(default) :: dr
    end function eta_phi_distance_vector4
    module subroutine lorentz_transformation_write (L, unit, testflag, ultra)
      class(lorentz_transformation_t), intent(in) :: L
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag, ultra
    end subroutine lorentz_transformation_write
    pure module function lorentz_transformation_get_components (L) result (a)
      type(lorentz_transformation_t), intent(in) :: L
      real(default), dimension(0:3,0:3) :: a
    end function lorentz_transformation_get_components
    elemental module function lorentz_transformation_inverse (L) result (IL)
      type(lorentz_transformation_t) :: IL
      type(lorentz_transformation_t), intent(in) :: L
    end function lorentz_transformation_inverse
    module function create_orthogonal (p_in) result (p_out)
      type(vector3_t), intent(in) :: p_in
      type(vector3_t) :: p_out
    end function create_orthogonal
    module function create_unit_vector (p_in) result (p_out)
      type(vector4_t), intent(in) :: p_in
      type(vector3_t) :: p_out
    end function create_unit_vector
    module function normalize(p) result (p_norm)
      type(vector3_t) :: p_norm
      type(vector3_t), intent(in) :: p
    end function normalize
    pure module function compute_resonance_mass (p, i_res_born, i_gluon) result (m)
      real(default) :: m
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in), dimension(:) :: i_res_born
      integer, intent(in), optional :: i_gluon
    end function compute_resonance_mass
    pure module function get_resonance_momentum &
         (p, i_res_born, i_gluon) result (p_res)
      type(vector4_t) :: p_res
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in), dimension(:) :: i_res_born
      integer, intent(in), optional :: i_gluon
    end function get_resonance_momentum
    module function create_two_particle_decay (s, p1, p2) result (p_rest)
      type(vector4_t), dimension(3) :: p_rest
      real(default), intent(in) :: s
      type(vector4_t), intent(in) :: p1, p2
    end function create_two_particle_decay
    module function create_three_particle_decay (p1, p2, p3) result (p_rest)
      type(vector4_t), dimension(4) :: p_rest
      type(vector4_t), intent(in) :: p1, p2, p3
    end function create_three_particle_decay
    recursive module subroutine generate_on_shell_decay (p_dec, &
        p_in, p_out, i_real, msq_in, jac, evaluate_special)
      type(vector4_t), intent(in) :: p_dec
      type(vector4_t), intent(in), dimension(:) :: p_in
      type(vector4_t), intent(inout), dimension(:) :: p_out
      integer, intent(in) :: i_real
      real(default), intent(in), optional :: msq_in
      real(default), intent(inout), optional :: jac
      procedure(evaluate_one_to_two_splitting_special), intent(in), &
            pointer, optional :: evaluate_special
    end subroutine generate_on_shell_decay
    elemental module function boost_from_rest_frame (p, m) result (L)
      type(lorentz_transformation_t) :: L
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: m
    end function boost_from_rest_frame
    elemental module function boost_from_rest_frame_vector3 (p, m) result (L)
      type(lorentz_transformation_t) :: L
      type(vector3_t), intent(in) :: p
      real(default), intent(in) :: m
    end function boost_from_rest_frame_vector3
    elemental module function boost_canonical (beta_gamma, k) result (L)
      type(lorentz_transformation_t) :: L
      real(default), intent(in) :: beta_gamma
      integer, intent(in) :: k
    end function boost_canonical
    elemental module function boost_generic (beta_gamma, axis) result (L)
      type(lorentz_transformation_t) :: L
      real(default), intent(in) :: beta_gamma
      type(vector3_t), intent(in) :: axis
    end function boost_generic
    elemental module function rotation_generic_cs (cp, sp, axis) result (R)
      type(lorentz_transformation_t) :: R
      real(default), intent(in) :: cp, sp
      type(vector3_t), intent(in) :: axis
    end function rotation_generic_cs
    elemental module function rotation_generic (axis) result (R)
      type(lorentz_transformation_t) :: R
      type(vector3_t), intent(in) :: axis
    end function rotation_generic
    elemental module function rotation_canonical_cs (cp, sp, k) result (R)
      type(lorentz_transformation_t) :: R
      real(default), intent(in) :: cp, sp
      integer, intent(in) :: k
    end function rotation_canonical_cs
    elemental module function rotation_canonical (phi, k) result (R)
      type(lorentz_transformation_t) :: R
      real(default), intent(in) :: phi
      integer, intent(in) :: k
    end function rotation_canonical
    elemental module function rotation_to_2nd_generic (p, q) result (R)
      type(lorentz_transformation_t) :: R
      type(vector3_t), intent(in) :: p, q
    end function rotation_to_2nd_generic
    elemental module function rotation_to_2nd_canonical (k, p) result (R)
      type(lorentz_transformation_t) :: R
      integer, intent(in) :: k
      type(vector3_t), intent(in) :: p
    end function rotation_to_2nd_canonical
    elemental module function transformation_rec_generic (axis, p1, p2, m) result (L)
      type(vector3_t), intent(in) :: axis
      type(vector4_t), intent(in) :: p1, p2
      real(default), intent(in) :: m
      type(lorentz_transformation_t) :: L
    end function transformation_rec_generic
    elemental module function transformation_rec_canonical (k, p1, p2, m) result (L)
      integer, intent(in) :: k
      type(vector4_t), intent(in) :: p1, p2
      real(default), intent(in) :: m
      type(lorentz_transformation_t) :: L
    end function transformation_rec_canonical
    elemental module function prod_LT_vector4 (L, p) result (np)
      type(vector4_t) :: np
      type(lorentz_transformation_t), intent(in) :: L
      type(vector4_t), intent(in) :: p
    end function prod_LT_vector4
    elemental module function prod_LT_LT (L1, L2) result (NL)
      type(lorentz_transformation_t) :: NL
      type(lorentz_transformation_t), intent(in) :: L1,L2
    end function prod_LT_LT
    elemental module function prod_vector4_LT (p, L) result (np)
      type(vector4_t) :: np
      type(vector4_t), intent(in) :: p
      type(lorentz_transformation_t), intent(in) :: L
    end function prod_vector4_LT
    elemental module function LT_compose_r3_r2_b3 &
         (cp, sp, ct, st, beta_gamma) result (L)
      type(lorentz_transformation_t) :: L
      real(default), intent(in) :: cp, sp, ct, st, beta_gamma
    end function LT_compose_r3_r2_b3
    elemental module function LT_compose_r2_r3_b3 &
         (ct, st, cp, sp, beta_gamma) result (L)
      type(lorentz_transformation_t) :: L
      real(default), intent(in) :: ct, st, cp, sp, beta_gamma
    end function LT_compose_r2_r3_b3
    elemental module function axis_from_p_r3_r2_b3 &
         (p, cp, sp, ct, st, beta_gamma) result (n)
      type(vector3_t) :: n
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: cp, sp, ct, st, beta_gamma
    end function axis_from_p_r3_r2_b3
    elemental module function axis_from_p_b3 (p, beta_gamma) result (n)
      type(vector3_t) :: n
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: beta_gamma
    end function axis_from_p_b3
    elemental module function lambda (m1sq, m2sq, m3sq)
      real(default) :: lambda
      real(default), intent(in) :: m1sq, m2sq, m3sq
    end function lambda
    module function colliding_momenta (sqrts, m, p_cm) result (p)
      type(vector4_t), dimension(2) :: p
      real(default), intent(in) :: sqrts
      real(default), dimension(2), intent(in), optional :: m
      real(default), intent(in), optional :: p_cm
    end function colliding_momenta
    elemental module subroutine pacify_vector3 (p, tolerance)
      type(vector3_t), intent(inout) :: p
      real(default), intent(in) :: tolerance
    end subroutine pacify_vector3
    elemental module subroutine pacify_vector4 (p, tolerance)
      type(vector4_t), intent(inout) :: p
      real(default), intent(in) :: tolerance
    end subroutine pacify_vector4
    elemental module subroutine pacify_LT (LT, tolerance)
      type(lorentz_transformation_t), intent(inout) :: LT
      real(default), intent(in) :: tolerance
    end subroutine pacify_LT
    module subroutine vector_set_reshuffle (p1, list, p2)
      type(vector4_t), intent(in), dimension(:), allocatable :: p1
      integer, intent(in), dimension(:), allocatable :: list
      type(vector4_t), intent(out), dimension(:), allocatable :: p2
    end subroutine vector_set_reshuffle
    module function vector_set_is_cms (p, n_in) result (is_cms)
      logical :: is_cms
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: n_in
    end function vector_set_is_cms
    module subroutine vector4_write_set (p, unit, show_mass, testflag, &
          check_conservation, ultra, n_in)
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_mass
      logical, intent(in), optional :: testflag, ultra
      logical, intent(in), optional :: check_conservation
      integer, intent(in), optional :: n_in
    end subroutine vector4_write_set
    module subroutine vector4_check_momentum_conservation (p, n_in, unit, &
       abs_smallness, rel_smallness, verbose)
      type(vector4_t), dimension(:), intent(in) :: p
      integer, intent(in) :: n_in
      integer, intent(in), optional :: unit
      real(default), intent(in), optional :: abs_smallness, rel_smallness
      logical, intent(in), optional :: verbose
    end subroutine vector4_check_momentum_conservation
    module subroutine spinor_product (p1, p2, prod1, prod2)
      type(vector4_t), intent(in) :: p1, p2
      complex(default), intent(out) :: prod1, prod2
    end subroutine spinor_product
  end interface

end module lorentz
