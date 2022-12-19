! parameters.QED.omega.f90 --
!
! Copyright (C) 1999-2022 by 
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

module parameters_qcd
  use kinds
  use constants

  implicit none
  private
  public :: init_parameters

  real(default), dimension(22), public :: mass, width
  complex(default), public :: gs

  complex(default), parameter :: ALPHAS = 0.12

contains

  subroutine init_parameters
    mass(1:22) = 0
    width(1:22) = 0
    gs = sqrt(4.0_default*PI*ALPHAS) / sqrt(2.0_default)
  end subroutine init_parameters

end module parameters_qcd
