! parameters.Test.f90
!
! Copyright (C) 1999-2022 by 
!
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
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
module parameters_test
  use kinds

  implicit none
  private
  public :: import_from_whizard, model_update_alpha_s

  real(default), dimension(2), public :: mass, width
  real(default), public :: gy

contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(3), intent(in) :: par_array
    integer, intent(in) :: scheme
    type :: parameter_set
       real(default) :: gy
       real(default) :: ms
       real(default) :: mf
    end type parameter_set
    type(parameter_set) :: par
    par%gy = par_array(1)
    par%ms = par_array(2)
    par%mf = par_array(3)
    mass = 0
    width = 0
    mass(1) = par%ms
    mass(2) = par%mf
    gy = par%gy
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
  end subroutine model_update_alpha_s
end module parameters_test
