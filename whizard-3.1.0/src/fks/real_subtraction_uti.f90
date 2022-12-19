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

module real_subtraction_uti
  use kinds, only: default

  use physics_defs
  use lorentz
  use numeric_utils
  use real_subtraction

  implicit none
  private

  public :: real_subtraction_1

contains

  subroutine real_subtraction_1 (u)
    integer, intent(in) :: u

    type(coll_subtraction_t) :: coll_sub
    real(default) :: sqme_coll
    type(vector4_t) :: p_res
    type(vector4_t), dimension(5) :: p_born
    real(default), dimension(4) :: k_perp
    real(default), dimension(4,4) :: b_munu
    integer :: mu, nu
    real(default) :: born, born_c
    integer, dimension(6) :: flst

    p_born(1)%p = [500, 0, 0, 500]
    p_born(2)%p = [500, 0, 0, -500]
    p_born(3)%p = [3.7755E+02, 2.2716E+02, -95.4172, 2.8608E+02]
    p_born(4)%p = [4.9529E+02, -2.739E+02, 84.8535, -4.0385E+02]
    p_born(5)%p = [1.2715E+02, 46.7375, 10.5637, 1.1778E+02]
    p_res = p_born(1) + p_born(2)
    flst = [11, -11 , -2, 2, -2, 2]

    b_munu(1, :) = [0., 0., 0., 0.]
    b_munu(2, :) = [0., 1., 1., 1.]
    b_munu(3, :) = [0., 1., 1., 1.]
    b_munu(4, :) = [0., 1., 1., 1.]

    k_perp = real_subtraction_compute_k_perp_fsr (p = p_born(5), phi = 0.5_default)
    born = - b_munu(1, 1) + b_munu(2, 2) + b_munu(3, 3) + b_munu(4, 4)
    born_c = 0.
    do mu = 1, 4
       do nu = 1, 4
          born_c = born_c + k_perp(mu) * k_perp(nu) * b_munu(mu, nu)
       end do
    end do

    write (u, "(A)") "* Test output: real_subtraction_1"
    write (u, "(A)") "*   Purpose: final-state collinear subtraction"
    write (u, "(A)")

    write (u, "(A, L1)") "* vanishing scalar-product of 3-momenta k_perp and p_born(emitter): ", &
         nearly_equal (dot_product (p_born(5)%p(1:3), k_perp(2:4)), 0._default)

    call coll_sub%init (n_alr = 1, n_in = 2)
    call coll_sub%set_parameters (CA, CF, TR)

    write (u, "(A)")
    write (u, "(A)") "* g -> qq splitting"
    write (u, "(A)")

    sqme_coll = coll_sub%compute_fsr(5, flst, p_res, p_born, &
         born, born_c, 0.5_default, 0.25_default, .false.)

    write (u, "(A,F15.12)") "ME: ", sqme_coll

    write (u, "(A)")
    write (u, "(A)") "* g -> gg splitting"
    write (u, "(A)")

    b_munu(1, :) = [0., 0., 0., 0.]
    b_munu(2, :) = [0., 0., 0., 1.]
    b_munu(3, :) = [0., 0., 1., 1.]
    b_munu(4, :) = [0., 0., 1., 1.]

    k_perp = real_subtraction_compute_k_perp_fsr (p = p_born(5), phi = 0.5_default)
    born = - b_munu(1, 1) + b_munu(2, 2) + b_munu(3, 3) + b_munu(4, 4)
    born_c = 0.
    do mu = 1, 4
       do nu = 1, 4
          born_c = born_c + k_perp(mu) * k_perp(nu) * b_munu(mu, nu)
       end do
    end do

    flst = [11, -11, 2, -2, 21, 21]
    sqme_coll = coll_sub%compute_fsr(5, flst, p_res, p_born, &
         born, born_c, 0.5_default, 0.25_default, .true.)

    write (u, "(A,F15.12)") "ME: ", sqme_coll

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: real_subtraction_1"
    write (u, "(A)")
  end subroutine real_subtraction_1


end module real_subtraction_uti
