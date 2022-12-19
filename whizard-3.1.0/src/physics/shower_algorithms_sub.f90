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

submodule (shower_algorithms) shower_algorithms_s

  use diagnostics
  use constants

  implicit none

contains

  module subroutine generate_vetoed (x, overestimator, true_function, &
         sudakov, inverse_sudakov, scale_min)
    real(default), dimension(:), intent(out) :: x
    !class(rng_t), intent(inout) :: rng
    procedure(XXX_function), pointer, intent(in) :: overestimator, true_function
    procedure(sudakov_p), pointer, intent(in) :: sudakov, inverse_sudakov
    real(default), intent(in) :: scale_min
    real(default) :: random, scale_max, scale
    scale_max = inverse_sudakov (one)
    do while (scale_max > scale_min)
       !call rng%generate (random)
       scale = inverse_sudakov (random * sudakov (scale_max))
       call generate_on_hypersphere (x, overestimator, scale)
       !call rng%generate (random)
       if (random < true_function (x) / overestimator (x)) then
          return !!! accept x
       end if
       scale_max = scale
    end do
  end subroutine generate_vetoed

  subroutine generate_on_hypersphere (x, overestimator, scale)
    real(default), dimension(:), intent(out) :: x
    procedure(XXX_function), pointer, intent(in) :: overestimator
    real(default), intent(in) :: scale
    call msg_bug ("generate_on_hypersphere: not implemented")
  end subroutine generate_on_hypersphere




end submodule shower_algorithms_s

