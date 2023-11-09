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

module recoil_kinematics

  use kinds, only: default
  use lorentz, only: vector4_t
  use lorentz, only: vector4_null
  use lorentz, only: vector4_moving
  use lorentz, only: vector3_moving
  use lorentz, only: transverse_part
  use lorentz, only: lorentz_transformation_t
  use lorentz, only: inverse
  use lorentz, only: boost
  use lorentz, only: transformation
  use lorentz, only: operator(+)
  use lorentz, only: operator(-)
  use lorentz, only: operator(*)
  use lorentz, only: operator(**)
  use lorentz, only: lambda

  implicit none
  private

  public :: generate_q2_recoil
  public :: solve_recoil
  public :: recoil_momenta
  public :: recoil_transformation
  public :: initial_transformation
  public :: generate_recoil





  interface
    elemental module subroutine generate_q2_recoil (s, x_bar, q2_max, m2, r, q2)
      real(default), intent(in) :: s
      real(default), intent(in) :: q2_max
      real(default), intent(in) :: x_bar
      real(default), intent(in) :: m2
      real(default), intent(in) :: r
      real(default), intent(out) :: q2
    end subroutine generate_q2_recoil
    module subroutine solve_recoil &
         (sqrts, xc, xcb, phi, q2, x, xb, cos_th, sin_th, ok)
      real(default), intent(in) :: sqrts
      real(default), dimension(2), intent(in) :: xc
      real(default), dimension(2), intent(in) :: xcb
      real(default), dimension(2), intent(in) :: phi
      real(default), dimension(2), intent(in) :: q2
      real(default), dimension(2), intent(out) :: x
      real(default), dimension(2), intent(out) :: xb
      real(default), dimension(2), intent(out) :: cos_th
      real(default), dimension(2), intent(out) :: sin_th
      logical, intent(out) :: ok
    end subroutine solve_recoil
    module subroutine recoil_momenta (sqrts, xc, xb, cos_th, sin_th, phi, mo, &
         km, qm, qo, ok)
      real(default), intent(in) :: sqrts
      real(default), dimension(2), intent(in) :: xc
      real(default), dimension(2), intent(in) :: xb
      real(default), dimension(2), intent(in) :: cos_th
      real(default), dimension(2), intent(in) :: sin_th
      real(default), dimension(2), intent(in) :: phi
      real(default), dimension(2), intent(in) :: mo
      type(vector4_t), dimension(2), intent(out) :: km
      type(vector4_t), dimension(2), intent(out) :: qm
      type(vector4_t), dimension(2), intent(out) :: qo
      logical, intent(out) :: ok
    end subroutine recoil_momenta
    module subroutine recoil_transformation (sqrts, xc, qo, lt)
      real(default), intent(in) :: sqrts
      real(default), dimension(2), intent(in) :: xc
      type(vector4_t), dimension(2), intent(in) :: qo
      type(lorentz_transformation_t), intent(out) :: lt
    end subroutine recoil_transformation
    module subroutine initial_transformation (p, sqrts, lt, ok)
      type(vector4_t), dimension(2), intent(in) :: p
      real(default), intent(out) :: sqrts
      type(lorentz_transformation_t), intent(out) :: lt
      logical, intent(out) :: ok
    end subroutine initial_transformation
    module subroutine generate_recoil &
         (sqrts, q_max, m, mo, xc, xcb, r, km, qm, qo, ok)
      real(default), intent(in) :: sqrts
      real(default), intent(in), dimension(2) :: q_max
      real(default), intent(in), dimension(2) :: m
      real(default), intent(in), dimension(2) :: mo
      real(default), intent(in), dimension(2) :: xc
      real(default), intent(in), dimension(2) :: xcb
      real(default), intent(in), dimension(4) :: r
      type(vector4_t), dimension(2), intent(out) :: km
      type(vector4_t), dimension(2), intent(out) :: qm
      type(vector4_t), dimension(2), intent(out) :: qo
      logical, intent(out) :: ok
    end subroutine generate_recoil
  end interface

end module recoil_kinematics
