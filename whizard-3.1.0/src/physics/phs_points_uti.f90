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

module phs_points_uti

  use kinds, only: default
  use constants, only: zero
  use format_defs, only: FMT_12
  use lorentz
  use phs_points

  implicit none
  private

  public :: phs_points_1

contains

  subroutine phs_points_1 (u)
    integer, intent(in) :: u
    type(vector4_t), dimension(8) :: tt_mom
    type(phs_point_t) :: phs_p
    type(vector4_t) :: p_sum
    type(vector4_t), dimension(:), allocatable :: p_tau, p_out

    write (u, "(A)")  "* Test output: phs_points_1"
    write (u, "(A)")  "*   Purpose: handling a 2->6 PSP"
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)")  "*   Setting up a 2->6 off-shell top PSP"
    write (u, "(A)")

    tt_mom(1) = [2.5000000000000000e+02_default, zero, zero, 2.4999999999947775e+02_default]
    tt_mom(2) = [2.5000000000000000e+02_default, zero, zero, -2.4999999999947775e+02_default]
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
    tt_mom(8) = [9.8432954734067863e+01_default, 9.5129959382432389e+01_default, &
         3.1145755197238953e+00_default, 2.5092565132081493e+01_default]
    phs_p = tt_mom

    write (u, "(A)")
    write (u, "(A)")  "*   Retrieving the size of PSP"
    write (u, "(A)")
    write (u, "(3x,A,I0)")  "Size PSP  = ", size (phs_p)

    write (u, "(A)")
    write (u, "(A)")  "*   Returning the set of 4-momenta from PSP"
    write (u, "(A)")
    p_out = phs_p%get ()
    write (u, "(3x,A)")  "set 4-mom.  = "
    call vector4_write_set (p_out, u, testflag = .true., ultra = .true.)

    write (u, "(A)")
    write (u, "(A)")  "*   Sum of momenta of PSP"
    write (u, "(A)")
    p_sum = sum (phs_p)
    call pacify (p_sum, tolerance = 1.e-12_default)
    write (u, "(3x,A)")  "Sum:"
    call p_sum%write (u)

    write (u, "(A)")
    write (u, "(A)")  "*   Reconstructing top/antitop candidate invariant masses from PSP"
    write (u, "(A)")
    write (u, "(3x,A," // FMT_12 // ")")  "m2(top)   = ", sqrt (phs_p%get_msq ([3,6,8]))
    write (u, "(3x,A," // FMT_12 // ")")  "m2(a-top) = ", sqrt (phs_p%get_msq ([4,5,7]))

    write (u, "(A)")
    write (u, "(A)")  "*   Select a specific 4-vector from PSP, here for a tau"
    write (u, "(A)")
    p_tau = phs_p%select ([7])
    write (u, "(3x,A)")  "p(tau):"
    call p_tau(1)%write (u, show_mass = .true., testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_points_1"

  end subroutine phs_points_1


end module phs_points_uti
