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

module su_algebra

  use kinds, only: default

  implicit none
  private

  public :: algebra_dimension
  public :: fundamental_dimension
  public :: helicity_value
  public :: helicity_index
  public :: is_cartan_generator
  public :: cartan_index
  public :: cartan_element
  public :: cartan_coeff
  public :: root_index
  public :: root_helicity

  interface
    module function algebra_dimension (s) result (n)
      integer :: n
      integer, intent(in) :: s
    end function algebra_dimension
    module function fundamental_dimension (s) result (d)
      integer :: d
      integer, intent(in) :: s
    end function fundamental_dimension
    module function helicity_value (s, i) result (h)
      integer :: h
      integer, intent(in) :: s, i
    end function helicity_value
    module function helicity_index (s, h) result (i)
      integer, intent(in) :: s, h
      integer :: i
    end function helicity_index
    elemental module function is_cartan_generator (s, i) result (cartan)
      logical :: cartan
      integer, intent(in) :: s, i
    end function is_cartan_generator
    elemental module function cartan_index (s, k) result (ci)
      integer :: ci
      integer, intent(in) :: s, k
    end function cartan_index
    module function cartan_element (s, h) result (a)
      real(default), dimension(:), allocatable :: a
      integer, intent(in) :: s, h
    end function cartan_element
    module function cartan_coeff (s, rd) result (a)
      real(default), dimension(:), allocatable :: a
      integer, intent(in) :: s
      real(default), dimension(:), intent(in) :: rd
    end function cartan_coeff
    module function root_index (s, h1, h2, r) result (ai)
      integer :: ai
      integer, intent(in) :: s, h1, h2
      logical :: r
    end function root_index
    module subroutine root_helicity (s, i, h1, h2, r)
      integer, intent(in) :: s, i
      integer, intent(out) :: h1, h2
      logical, intent(out) :: r
    end subroutine root_helicity
  end interface

end module su_algebra
