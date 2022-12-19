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

module lorentz_uti

  use kinds, only: default
  use constants, only: zero, Pi
  use format_defs, only: FMT_12  
  use lorentz

  implicit none
  private

  public :: lorentz_1
  public :: lorentz_2
  public :: lorentz_3
  public :: lorentz_4
  public :: lorentz_5

contains

  subroutine lorentz_1 (u)
    integer, intent(in) :: u
    type(vector3_t) :: v3_1, v3_2

    write (u, "(A)")  "* Test output: lorentz_1"
    write (u, "(A)")  "*   Purpose: testing vector3_t"
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)")  "* Null 3-vector"
    write (u, "(A)")
    call vector3_write (vector3_null, u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Canonical 3-vector"
    write (u, "(A)")
    call vector3_write (vector3_canonical (1), u, testflag = .true.)
    call vector3_write (vector3_canonical (2), u, testflag = .true.)
    call vector3_write (vector3_canonical (3), u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Canonical moving 3-vector"
    write (u, "(A)")
    call vector3_write (vector3_moving (42._default, 1), u, testflag = .true.)
    call vector3_write (vector3_moving (42._default, 2), u, testflag = .true.)
    call vector3_write (vector3_moving (42._default, 3), u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Generic moving 3-vector"
    write (u, "(A)")
    call vector3_write (vector3_moving ([3._default, 4._default, 5._default]), &
         u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Simple algebra with 3-vectors"
    write (u, "(A)")
    v3_1 = vector3_moving ([3._default, 4._default, 5._default])
    v3_2 = vector3_moving ([-2._default, 5._default, -1._default])
    write (u, "(1x,A)")  "v3_1:"
    call vector3_write (v3_1, u, testflag=.true.)
    write (u, "(1x,A)")  "v3_2:"
    call vector3_write (v3_2, u, testflag=.true.)
    write (u, "(1x,A)")  "-v3_1:"
    call vector3_write (-v3_1, u, testflag=.true.)
    write (u, "(1x,A)")  "v3_1 / |v3_1|:"
    call vector3_write (direction (v3_1), u, testflag=.true.)
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1(x): ", vector3_get_component (v3_1, 1)
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1(y): ", vector3_get_component (v3_1, 2)
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1(z): ", vector3_get_component (v3_1, 3)
    write (u, "(1x,A)")  "v3_1 + v3_2:"
    call vector3_write (v3_1 + v3_2, u, testflag=.true.)
    write (u, "(1x,A)")  "v3_1 - v3_2:"
    call vector3_write (v3_1 - v3_2, u, testflag=.true.)
    write (u, "(1x,A,L1)")  "v3_1 == v3_2: ", v3_1 == v3_2
    write (u, "(1x,A,L1)")  "v3_1 /= v3_2: ", v3_1 /= v3_2
    write (u, "(1x,A)")  "2 * v3_1:"
    call vector3_write (2._default * v3_1, u, testflag=.true.)
    write (u, "(1x,A)")  "v3_2 / 4:"
    call vector3_write (v3_2 / 4, u, testflag=.true.)
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1, azimuth (radians):", azimuthal_angle (v3_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1, azimuth (degrees):", &
         azimuthal_angle_deg (v3_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1, polar (radians)  :", polar_angle (v3_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1, polar (degrees)  :", &
         polar_angle_deg (v3_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1, cosine polar     :", &
         polar_angle_ct (v3_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1, energy w. mass=1 :", &
         energy (v3_1, 1._default)
    write (u, "(1x,A)")  "3-vector orthogonal to v3_1:"
    call vector3_write (create_orthogonal (v3_1), u, testflag=.true.)
    write (u, "(1x,A)")  "unit 3-vector from v3_1:"

    write (u, "(A)")
    write (u, "(A)")  "* Dot and cross product"
    write (u, "(A)")
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1 * v3_2: ", v3_1 * v3_2
    write (u, "(1x,A," // FMT_12 // ")")  "v3_1**3    : ", v3_1**3
    write (u, "(1x,A)")  "v3_1 x v3_2:"
    call vector3_write (cross_product (v3_1, v3_2), u, testflag=.true.)
    write (u, "(1x,A," // FMT_12 // ")")  "enclosed angle (radians):", &
         enclosed_angle (v3_1, v3_2)
    write (u, "(1x,A," // FMT_12 // ")")  "enclosed angle (degrees):", &
         enclosed_angle_deg (v3_1, v3_2)
    write (u, "(1x,A," // FMT_12 // ")")  "cosine (enclosed angle) :", &
         enclosed_angle_ct (v3_1, v3_2)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: lorentz_1"

  end subroutine lorentz_1

  subroutine lorentz_2 (u)
    integer, intent(in) :: u
    type(vector3_t) :: v3_1, v3_2
    type(vector4_t) :: v4_1, v4_2, v4_1_inv

    write (u, "(A)")  "* Test output: lorentz_2"
    write (u, "(A)")  "*   Purpose: testing vector4_t"
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)")  "* Null 4-vector"
    write (u, "(A)")
    call vector4_write (vector4_null, u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Canonical 4-vector"
    write (u, "(A)")
    call vector4_write (vector4_canonical (0), u, testflag = .true.)
    call vector4_write (vector4_canonical (1), u, testflag = .true.)
    call vector4_write (vector4_canonical (2), u, testflag = .true.)
    call vector4_write (vector4_canonical (3), u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* 4-vector at rest with mass m = 17"
    write (u, "(A)")
    call vector4_write (vector4_at_rest (17._default), u, testflag = .true.)
    
    write (u, "(A)")
    write (u, "(A)")  "* Canonical moving 4-vector"
    write (u, "(A)")
    call vector4_write (vector4_moving (17._default, 42._default, 1), u, testflag = .true.)
    call vector4_write (vector4_moving (17._default, 42._default, 2), u, testflag = .true.)
    call vector4_write (vector4_moving (17._default, 42._default, 3), u, testflag = .true.)
     
    write (u, "(A)")
    write (u, "(A)")  "* Generic moving 4-vector"
    write (u, "(A)")
    v3_1 = [3._default, 4._default, 5._default]
    call vector4_write (vector4_moving (17._default, v3_1), u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Simple algebra with 4-vectors"
    write (u, "(A)")
    v3_2 = [-2._default, 5._default, -1._default]
    v4_1 = vector4_moving (8._default, v3_1)    
    v4_2 = vector4_moving (zero, v3_2)
    write (u, "(1x,A)")  "v4_1:"
    call vector4_write (v4_1, u, testflag=.true.)
    write (u, "(1x,A)")  "v4_2:"
    call vector4_write (v4_2, u, testflag=.true.)
    write (u, "(1x,A)")  "-v4_1:"
    call vector4_write (-v4_1, u, testflag=.true.)
    v4_1_inv = v4_1
    call vector4_invert_direction (v4_1_inv)
    write (u, "(1x,A)")  "v4_1, inverted direction:"
    call vector4_write (v4_1_inv, u, testflag=.true.)
    write (u, "(1x,A)")  "(v4_1)_spatial / |(v4_1)_spatial|:"
    call vector3_write (direction (v4_1), u, testflag=.true.)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1(E): ", energy (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1(x): ", vector4_get_component (v4_1, 1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1(y): ", vector4_get_component (v4_1, 2)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1(z): ", vector4_get_component (v4_1, 3)
    write (u, "(1x,A)")  "space_part (v4_1):"
    call vector3_write (space_part (v4_1), u, testflag=.true.)
    write (u, "(1x,A," // FMT_12 // ")")  "norm space_part (v4_1): ", &
         space_part_norm (v4_1)
    write (u, "(1x,A)") "unit vector from v4_1:"
    call vector3_write (create_unit_vector (v4_1), u, testflag = .true.)
    write (u, "(1x,A)")  "v4_1 + v4_2:"
    call vector4_write (v4_1 + v4_2, u, testflag=.true.)
    write (u, "(1x,A)")  "v4_1 - v4_2:"
    call vector4_write (v4_1 - v4_2, u, testflag=.true.)
    write (u, "(1x,A,L1)")  "v4_1 == v4_2: ", v4_1 == v4_2
    write (u, "(1x,A,L1)")  "v4_1 /= v4_2: ", v4_1 /= v4_2
    write (u, "(1x,A)")  "2 * v4_1:"
    call vector4_write (2._default * v4_1, u, testflag=.true.)
    write (u, "(1x,A)")  "v4_2 / 4:"
    call vector4_write (v4_2 / 4, u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Angles and kinematic properties of 4-vectors"
    write (u, "(A)")
   
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, azimuth (radians):", azimuthal_angle (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, azimuth (degrees):", &
         azimuthal_angle_deg (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, polar (radians)  :", polar_angle (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, polar (degrees)  :", &
         polar_angle_deg (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, cosine polar     :", &
         polar_angle_ct (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, invariant mass   :", &
         invariant_mass (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, invariant mass sq:", &
         invariant_mass_squared (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_2, invariant mass   :", &
         invariant_mass (v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_2, invariant mass sq:", &
         invariant_mass_squared (v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, transverse mass  :", &
         transverse_mass (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, rapidity         :", &
         rapidity (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, pseudorapidity   :", &
         pseudorapidity (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, pT               :", &
         transverse_part (v4_1)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1, pL               :", &
         longitudinal_part (v4_1)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: lorentz_2"

  end subroutine lorentz_2

  subroutine lorentz_3 (u)
    integer, intent(in) :: u
    type(vector3_t) :: v3_1, v3_2
    type(vector4_t) :: v4_1, v4_2

    write (u, "(A)")  "* Test output: lorentz_3"
    write (u, "(A)")  "*   Purpose: testing bilinear functions of 4-vectors"
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)")  "* Products and distances of 4-vectors"
    write (u, "(A)")
    v3_1 = [3._default, 4._default, 5._default]
    v3_2 = [-2._default, 5._default, -1._default]
    v4_1 = vector4_moving (8._default, v3_1)    
    v4_2 = vector4_moving (6._default, v3_2)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1 * v4_2: ", v4_1 * v4_2
    write (u, "(1x,A," // FMT_12 // ")")  "rapidity distance       :", &
         rapidity_distance (v4_1, v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "pseudorapidity distance :", &
         pseudorapidity_distance (v4_1, v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "eta phi distance        :", &
         eta_phi_distance (v4_1, v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "enclosed angle (radians):", &
         enclosed_angle (v4_1, v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "enclosed angle (degrees):", &
         enclosed_angle_deg (v4_1, v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "cosine (enclosed angle) :", &
         enclosed_angle_ct (v4_1, v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "rest frame theta (rad)  :", &
         enclosed_angle_rest_frame (v4_1, v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "rest frame theta (deg)  :", &
         enclosed_angle_deg_rest_frame (v4_1, v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "rest frame cosine(theta):", &
         enclosed_angle_ct_rest_frame (v4_1, v4_2)
    write (u, "(1x,A," // FMT_12 // ")")  "v4_1_T w.r.t. v4_2      :", &
         transverse_part (v4_1, v4_2)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: lorentz_3"

  end subroutine lorentz_3

  subroutine lorentz_4 (u)
    integer, intent(in) :: u
    type(vector3_t) :: v3_1, v3_2
    type(vector4_t) :: v4
    type(lorentz_transformation_t) :: LT
    real(default) :: tol

    write (u, "(A)")  "* Test output: lorentz_4"
    write (u, "(A)")  "*   Purpose: testing Lorentz transformations"
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)")  "* Basic Lorentz transformatios"
    write (u, "(A)")
    write (u, "(1x,A)")  "LT = 1:"
    call lorentz_transformation_write (identity, u, testflag=.true.)
    write (u, "(A)")
    write (u, "(1x,A)")  "LT = space reflection:"
    call lorentz_transformation_write (space_reflection, u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Lorentz transformations: rotations"
    write (u, "(A)")
    v3_1 = [1._default, 2._default, 3._default]
    v3_2 = [-2._default, 1._default, -5._default]
    tol = 1.e-12_default
    write (u, "(1x,A)")  "Rotation of Pi/4 around 1-axis, def. by cos and sin:"
    LT = rotation (0.707107_default, 0.707107_default, 1)
    call pacify (LT, tol)
    call lorentz_transformation_write (LT, u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation of Pi/4 around 2-axis, def. by cos and sin:"
    LT = rotation (0.707107_default, 0.707107_default, 2)
    call pacify (LT, tol)
    call lorentz_transformation_write (LT, u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation of Pi/4 around 3-axis, def. by cos and sin:"
    LT = rotation (0.707107_default, 0.707107_default, 3)
    call pacify (LT, tol)
    call lorentz_transformation_write (LT, u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation of Pi/4 around 1-axis, def. by angle:"
    LT = rotation (Pi/4._default, 1)
    call pacify (LT, tol)
    call lorentz_transformation_write (LT, u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation of Pi/4 around 2-axis, def. by angle:"
    LT = rotation (Pi/4._default, 2)
    call pacify (LT, tol)
    call lorentz_transformation_write (LT, u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation of Pi/4 around 3-axis, def. by angle:"
    LT = rotation (Pi/4._default, 3)
    call pacify (LT, tol)
    call lorentz_transformation_write (LT, u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation of Pi/4 around axis = (1,2,3):"
    call lorentz_transformation_write (rotation (0.707107_default, 0.707107_default, &
         normalize (v3_1)), u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation in plane to axis = (1,2,3), angle given by length of axis:"
    call lorentz_transformation_write (rotation (v3_1), u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation from v3_1=(1,2,3) to v3_2=(-2,1,-5):"
    call lorentz_transformation_write (rotation_to_2nd (v3_1,v3_2), u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation from 1-axis to v3_2=(-2,1,-5):"
    call lorentz_transformation_write (rotation_to_2nd (1,v3_2), u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation from 2-axis to v3_2=(-2,1,-5):"
    call lorentz_transformation_write (rotation_to_2nd (2,v3_2), u, testflag=.true.)
    write (u, "(1x,A)")  "Rotation from 3-axis to v3_2=(-2,1,-5):"
    call lorentz_transformation_write (rotation_to_2nd (3,v3_2), u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Lorentz transformations: boosts"
    write (u, "(A)")
    write (u, "(1x,A)")  "Boost from rest frame to 3-vector, mass m=10:"
    call lorentz_transformation_write (boost (v3_1, 10._default), u, testflag=.true.)
    write (u, "(1x,A)")  "Boost from rest frame to 4-vector, mass m=10:"
    v4 = vector4_moving (42._default, v3_1)
    call lorentz_transformation_write (boost (v4, 10._default), u, testflag=.true.)
    write (u, "(1x,A)")  "Boost along 1-axis, beta*gamma = 12"
    call lorentz_transformation_write (boost (12._default, 1), u, testflag=.true.)
    write (u, "(1x,A)")  "Boost along 2-axis, beta*gamma = 12"
    call lorentz_transformation_write (boost (12._default, 2), u, testflag=.true.)
    write (u, "(1x,A)")  "Boost along 3-axis, beta*gamma = 12"
    call lorentz_transformation_write (boost (12._default, 3), u, testflag=.true.)
    write (u, "(1x,A)")  "Boost along axis=(1,2,3), beta*gamma = 12"
    call lorentz_transformation_write (boost (12._default, v3_1), u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: lorentz_4"

  end subroutine lorentz_4

  subroutine lorentz_5 (u)
    integer, intent(in) :: u
    type(vector4_t), dimension(2) :: p
    real(default), dimension(2) :: m
    real(default) :: sqrts
    type(vector4_t), dimension(8) :: tt_mom
    type(vector4_t), dimension(:), allocatable :: tin, tout
    integer, dimension(:), allocatable :: shuffle
    
    write (u, "(A)")  "* Test output: lorentz_5"
    write (u, "(A)")  "*   Purpose: testing additional kinematics and sets of 4-vectors"
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)")  "* Colliding momenta, 13 TeV, massless"
    write (u, "(A)")
    sqrts = 13000._default
    p = colliding_momenta (sqrts)
    call vector4_write (p(1), u, testflag=.true.)
    call vector4_write (p(2), u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Colliding momenta, 10 GeV, massive muons"
    write (u, "(A)")
    sqrts = 10._default
    m = [0.1057_default, 0.1057_default]
    p = colliding_momenta (sqrts, m)
    call vector4_write (p(1), u, testflag=.true.)
    call vector4_write (p(2), u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Kinematical function lambda"
    write (u, "(A)")

    write (u, "(1x,A," // FMT_12 // ")")  "s = 172.3**2, m1 = 4.2, m2 = 80.418:", &
         lambda (172.3_default**2, 4.2_default**2, 80.418_default**2)

    write (u, "(A)")
    write (u, "(A)")  "* Test vector_set"
    write (u, "(A)")

    tt_mom(1) = [2.5000000000000000e+02_default, zero, zero, 2.4999999999947777e+02_default]
    tt_mom(2) = [2.5000000000000000e+02_default, zero, zero, -2.4999999999947777e+02_default]
    tt_mom(3) = [1.1557492413664579e+02_default, 3.9011599241011098e+01_default, &
         -6.4278142734963140e+01_default, 8.7671766153043137e+01_default]
    tt_mom(4) = [1.4617918132729235e+02_default, -1.0947970597860679e+02_default, &
         1.5484441802571380e+01_default, -9.5525593923398418e+01_default]
    tt_mom(5) = [5.2637589215119526e+01_default, -4.7413198564695762e+01_default, &
         1.0087885417286579e+01_default, 2.0516525153079229e+01_default]
    tt_mom(6) = [5.4760292922264796e+01_default, 1.5197406985690520e+01_default, &
         5.1527071739328015e+01_default, -1.0615525413924287e+01_default]
    tt_mom(7) = [3.2415057664609684e+01_default, 7.5539389341684711e+00_default, &
         -1.5935831743946720e+01_default, -2.7139737100881156e+01_default]
    tt_mom(8) = [9.8432954734067863E+01_default, 9.5129959382432399e+01_default, &
         3.1145755197238966e+00_default, 2.5092565132081496e+01_default]
    write (u, "(1x,A)")  "Write routine for vector sets, maximal compression:"
    call vector4_write_set (tt_mom, u, show_mass=.true., testflag=.true., &
         check_conservation=.true., ultra=.true.)
    write (u, "(1x,A,L1)")  "Vector set is CMS frame: ", vector_set_is_cms (tt_mom, 2)
    write (u, "(1x,A)")  "Reshuffle vector set, final state inverted:"
    allocate (tin (8))
    tin = tt_mom
    allocate (shuffle (8), source = [1,2,8,7,6,5,4,3])
    call vector_set_reshuffle (tin, shuffle, tout)
    call vector4_write_set (tout, u, show_mass=.true., testflag=.true., &
         check_conservation=.true., ultra=.true.)
    write (u, "(1x,A)")  "Vector set, check momentum conservation:"
    call vector4_check_momentum_conservation (tt_mom, 2, u, &
         abs_smallness = 1.e-12_default, verbose=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: lorentz_5"

  end subroutine lorentz_5


end module lorentz_uti
