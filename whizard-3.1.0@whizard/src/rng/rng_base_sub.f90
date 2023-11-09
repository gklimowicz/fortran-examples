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

submodule (rng_base) rng_base_s

  use constants, only: TWOPI

  implicit none

contains

  module subroutine rng_generate_gaussian_single (rng, x)
    class(rng_t), intent(inout) :: rng
    real(default), intent(out) :: x
    real(default), dimension(2) :: u
    call rng%generate (u)
    x = sin (twopi * u(1)) * sqrt (- 2 * log (u(2)))
  end subroutine rng_generate_gaussian_single

  module subroutine rng_generate_gaussian_array (rng, x)
    class(rng_t), intent(inout) :: rng
    real(default), dimension(:), intent(out) :: x
    integer :: i
    do i = 1, size (x)
       call rng%generate_gaussian (x(i))
    end do
  end subroutine rng_generate_gaussian_array


end submodule rng_base_s

