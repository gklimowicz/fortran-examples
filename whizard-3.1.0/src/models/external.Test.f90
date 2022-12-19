! $Id: external_Test.f90 2364 2010-04-20 12:47:06Z cnspeckn $
!
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     Christian Speckner <cnspeckn@googlemail.com>
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
! External parameters for 'test' model.
! We compute the mass of the fermion in terms of the mass of the scalar.
! This could be done by a derived parameter within WHIZARD, but this
! code is just for testing the feature
subroutine init_external_parameters (par) bind (C)
  use iso_c_binding
  use kinds
  real(c_default_float), dimension(*), intent(inout) :: par
  real(default) :: ms, mf, ff
  ! Take the parameter(s) that are needed
  ms = par(2)
  ff = par(3)
  ! Compute the external parameter(s)
  mf = ff * ms
  ! Put the results back into the parameter array
  par(4) = mf
end subroutine init_external_parameters
