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

module parameters_sym
  use kinds
  use constants

  implicit none
  private
  public :: init_parameters

  integer, parameter, public :: NG = 9
  complex(default), dimension(NG,NG), public :: g_saa
  complex(default), dimension(NG,NG,NG), public :: g3, ig3, g_saaa
  complex(default), dimension(NG,NG,NG,NG), public :: g4

  complex(default), parameter :: ALPHAS = 0.12

contains

  subroutine init_parameters
    integer :: i
    complex(default) :: g, ghgg
    g = sqrt(4.0_default*PI*ALPHAS) / sqrt(2.0_default)
    ghgg = ALPHAS
    g_saa = 0
    g_saaa = 0
    g3 = 0
    ig3 = 0
    g4 = 0
    forall (i = 1:NG)
       g_saa(i,i) = ghgg
       g_saaa(i,i,i) = ghgg*g
       g3(i,i,i) = g
       ig3(i,i,i) = (0, 1) * g
       g4(i,i,i,i) = g**2
    end forall
  end subroutine init_parameters

end module parameters_sym
