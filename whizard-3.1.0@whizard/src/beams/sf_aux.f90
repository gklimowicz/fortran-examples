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

module sf_aux

  use kinds, only: default
  use constants, only: twopi
  use lorentz

  implicit none
  private

  public :: splitting_data_t
  public :: on_shell

  integer, parameter, public :: KEEP_ENERGY = 0, KEEP_MOMENTUM = 1

  type :: splitting_data_t
!     private
     logical :: collinear = .false.
     real(default) :: x0 = 0
     real(default) :: x1
     real(default) :: t0
     real(default) :: t1
     real(default) :: phi0 = 0
     real(default) :: phi1 = twopi
     real(default) :: E, p, s, u, m2
     real(default) :: x, xb, pb
     real(default) :: t = 0
     real(default) :: phi = 0
   contains
     procedure :: write => splitting_data_write
     procedure :: init => splitting_data_init
     procedure :: get_x_bounds => splitting_get_x_bounds
     procedure :: set_t_bounds => splitting_set_t_bounds
     procedure :: sample_t => splitting_sample_t
     procedure :: inverse_t => splitting_inverse_t
     procedure :: sample_phi => splitting_sample_phi
     procedure :: inverse_phi => splitting_inverse_phi
     procedure :: split_momentum => splitting_split_momentum
     procedure :: recover => splitting_recover
     procedure :: get_x => splitting_get_x
     procedure :: get_xb => splitting_get_xb
  end type splitting_data_t


  interface
    module subroutine splitting_data_write (d, unit)
      class(splitting_data_t), intent(in) :: d
      integer, intent(in), optional :: unit
    end subroutine splitting_data_write
    module subroutine splitting_data_init (d, k, mk2, mr2, mo2, collinear)
      class(splitting_data_t), intent(out) :: d
      type(vector4_t), intent(in) :: k
      real(default), intent(in) :: mk2, mr2, mo2
      logical, intent(in), optional :: collinear
    end subroutine splitting_data_init
    module function splitting_get_x_bounds (d) result (x)
      class(splitting_data_t), intent(in) :: d
      real(default), dimension(2) :: x
    end function splitting_get_x_bounds
    elemental module subroutine splitting_set_t_bounds (d, x, xb)
      class(splitting_data_t), intent(inout) :: d
      real(default), intent(in), optional :: x, xb
    end subroutine splitting_set_t_bounds
    module subroutine splitting_sample_t (d, r, t0, t1)
      class(splitting_data_t), intent(inout) :: d
      real(default), intent(in) :: r
      real(default), intent(in), optional :: t0, t1
    end subroutine splitting_sample_t
    module subroutine splitting_inverse_t (d, r, t0, t1)
      class(splitting_data_t), intent(in) :: d
      real(default), intent(out) :: r
      real(default), intent(in), optional :: t0, t1
    end subroutine splitting_inverse_t
    module subroutine splitting_sample_phi (d, r)
      class(splitting_data_t), intent(inout) :: d
      real(default), intent(in) :: r
    end subroutine splitting_sample_phi
    module subroutine splitting_inverse_phi (d, r)
      class(splitting_data_t), intent(in) :: d
      real(default), intent(out) :: r
    end subroutine splitting_inverse_phi
    module function splitting_split_momentum (d, k) result (q)
      class(splitting_data_t), intent(in) :: d
      type(vector4_t), dimension(2) :: q
      type(vector4_t), intent(in) :: k
    end function splitting_split_momentum
    elemental module subroutine on_shell (p, m2, keep)
      type(vector4_t), intent(inout) :: p
      real(default), intent(in) :: m2
      integer, intent(in) :: keep
    end subroutine on_shell
    module subroutine splitting_recover (d, k, q, keep)
      class(splitting_data_t), intent(inout) :: d
      type(vector4_t), intent(in) :: k
      type(vector4_t), dimension(2), intent(in) :: q
      integer, intent(in) :: keep
    end subroutine splitting_recover
    module function splitting_get_x (sd) result (x)
      class(splitting_data_t), intent(in) :: sd
      real(default) :: x
    end function splitting_get_x
    module function splitting_get_xb (sd) result (xb)
      class(splitting_data_t), intent(in) :: sd
      real(default) :: xb
    end function splitting_get_xb
  end interface

end module sf_aux

