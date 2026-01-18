! parameters.QCD.f90
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
  real(default), dimension(21), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs
  public :: import_from_whizard, model_update_alpha_s
contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(8), intent(in) :: par_array
    integer, intent(in) :: scheme
    integer, parameter :: scheme_massive = 2
    type :: parameter_set
       real(default) :: alphas
       real(default) :: md = 0 !!! par_array(1) locked in default scheme
       real(default) :: mu = 0 !!! par_array(2) locked in default scheme
       real(default) :: ms
       real(default) :: mc
       real(default) :: mb
       real(default) :: mtop
       real(default) :: wtop
    end type parameter_set
    type(parameter_set) :: par
    par%alphas = par_array(1)
    select case (scheme)
    case (scheme_massive)
       par%md = par_array(2)
       par%mu = par_array(3)
    end select
    par%ms = par_array(4)
    par%mc = par_array(5)
    par%mb = par_array(6)
    par%mtop = par_array(7)
    par%wtop = par_array(8)
    mass(1:21) = 0
    width(1:21) = 0
    mass(1) = par%md
    mass(2) = par%mu
    mass(3) = par%ms
    mass(4) = par%mc
    mass(5) = par%mb
    mass(6) = par%mtop
    width(6) = par%wtop    
    ! color-flow basis: gs is divided by sqrt2
    gs = sqrt(2.0_default*PI*par%alphas)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs    
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs        
  end subroutine model_update_alpha_s
end module parameters_qcd
