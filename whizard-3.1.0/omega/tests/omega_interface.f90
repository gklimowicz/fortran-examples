! omega_interface.f90 --
! omega_interface.f90 -- package the O'Mega interface functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

module omega_interface

  implicit none
  private

  type omega_procedures
     procedure(number_particles_in), nopass, pointer :: number_particles_in => NULL()
     procedure(number_particles_out), nopass, pointer :: number_particles_out => NULL()
     procedure(number_spin_states), nopass, pointer :: number_spin_states => NULL()
     procedure(spin_states), nopass, pointer :: spin_states => NULL()
     procedure(number_flavor_states), nopass, pointer :: number_flavor_states => NULL()
     procedure(flavor_states), nopass, pointer :: flavor_states => NULL()
     procedure(external_masses), nopass, pointer :: external_masses => NULL()
     procedure(number_color_indices), nopass, pointer :: number_color_indices => NULL()
     procedure(number_color_flows), nopass, pointer :: number_color_flows => NULL()
     procedure(color_flows), nopass, pointer :: color_flows => NULL()
     procedure(number_color_factors), nopass, pointer :: number_color_factors => NULL()
     procedure(color_factors), nopass, pointer :: color_factors => NULL()
     procedure(color_sum), nopass, pointer :: color_sum => NULL()
     procedure(new_event), nopass, pointer :: new_event => NULL()
     procedure(reset_helicity_selection), nopass, pointer :: reset_helicity_selection => NULL()
     procedure(is_allowed), nopass, pointer :: is_allowed => NULL()
     procedure(get_amplitude), nopass, pointer :: get_amplitude => NULL()
  end type omega_procedures

  public :: omega_procedures

  abstract interface

     pure function number_particles_in () result (n)
       integer :: n
     end function number_particles_in

     pure function number_particles_out () result (n)
       integer :: n
     end function number_particles_out

     pure function number_spin_states () result (n)
       integer :: n
     end function number_spin_states

     pure subroutine spin_states (a)
       integer, dimension(:,:), intent(out) :: a
     end subroutine spin_states

     pure function number_flavor_states () result (n)
       integer :: n
     end function number_flavor_states

     pure subroutine flavor_states (a)
       integer, dimension(:,:), intent(out) :: a
     end subroutine flavor_states

     pure subroutine external_masses (m, flv)
       use kinds
       real(kind=default), dimension(:), intent(out) :: m
       integer, intent(in) :: flv
     end subroutine external_masses

     pure function number_color_indices () result (n)
       integer :: n
     end function number_color_indices

     pure function number_color_flows () result (n)
       integer :: n
     end function number_color_flows

     pure subroutine color_flows (a, g)
       integer, dimension(:,:,:), intent(out) :: a
       logical, dimension(:,:), intent(out) :: g
     end subroutine color_flows

     pure function number_color_factors () result (n)
       integer :: n
     end function number_color_factors

     pure subroutine color_factors (cf)
       use omega_color
       type(omega_color_factor), dimension(:), intent(out) :: cf
     end subroutine color_factors

     !pure unless OpenMP
     !pure function color_sum (flv, hel) result (amp2)
     function color_sum (flv, hel) result (amp2)
       use kinds
       integer, intent(in) :: flv, hel
       real(kind=default) :: amp2
     end function color_sum

     subroutine new_event (p)
       use kinds
       real(kind=default), dimension(0:3,*), intent(in) :: p
     end subroutine new_event

     subroutine reset_helicity_selection (threshold, cutoff)
       use kinds
       real(kind=default), intent(in) :: threshold
       integer, intent(in) :: cutoff
     end subroutine reset_helicity_selection

     pure function is_allowed (flv, hel, col) result (yorn)
       logical :: yorn
       integer, intent(in) :: flv, hel, col
     end function is_allowed

     pure function get_amplitude (flv, hel, col) result (amp_result)
       use kinds
       complex(kind=default) :: amp_result
       integer, intent(in) :: flv, hel, col
     end function get_amplitude

  end interface

end module omega_interface

