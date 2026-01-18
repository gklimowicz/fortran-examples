! parameters.QED.f90
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
module parameters_qed
  use kinds
  use constants

  implicit none
  private
  public :: import_from_whizard, model_update_alpha_s

  real(default), dimension(22), public :: mass, width
  complex(default), public :: qlep

contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(4), intent(in) :: par_array
    integer, intent(in) :: scheme
    real(default) :: e, qelep
    type :: parameter_set
       real(default) :: ee
       real(default) :: me
       real(default) :: mmu
       real(default) :: mtau
    end type parameter_set
    type(parameter_set) :: par
    par%ee = par_array(1)
    par%me = par_array(2)
    par%mmu = par_array(3)
    par%mtau = par_array(4)
    mass(1:22) = 0
    width(1:22) = 0
    mass(11) = par%me
    mass(13) = par%mmu
    mass(15) = par%mtau
    e = par%ee
    qelep = - 1
    qlep = - e * qelep
  end subroutine import_from_whizard
  
  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s    
  end subroutine model_update_alpha_s
end module parameters_qed
