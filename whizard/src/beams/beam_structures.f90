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

module beam_structures

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use polarizations

  implicit none
  private

  public :: beam_structure_t

  type :: beam_structure_entry_t
     logical :: is_valid = .false.
     type(string_t) :: name
   contains
     procedure :: to_string => beam_structure_entry_to_string
  end type beam_structure_entry_t

  type :: beam_structure_record_t
     type(beam_structure_entry_t), dimension(:), allocatable :: entry
  end type beam_structure_record_t

  type :: beam_structure_t
     private
     integer :: n_beam = 0
     type(string_t), dimension(:), allocatable :: prt
     type(beam_structure_record_t), dimension(:), allocatable :: record
     type(smatrix_t), dimension(:), allocatable :: smatrix
     real(default), dimension(:), allocatable :: pol_f
     real(default), dimension(:), allocatable :: p
     real(default), dimension(:), allocatable :: theta
     real(default), dimension(:), allocatable :: phi
   contains
     procedure :: final_sf => beam_structure_final_sf
     procedure :: write => beam_structure_write
     procedure :: to_string => beam_structure_to_string
     procedure :: init_sf => beam_structure_init_sf
     procedure :: set_sf => beam_structure_set_sf
     procedure :: expand => beam_structure_expand
     procedure :: final_pol => beam_structure_final_pol
     procedure :: init_pol => beam_structure_init_pol
     procedure :: has_polarized_beams => beam_structure_has_polarized_beams
     procedure :: set_smatrix => beam_structure_set_smatrix
     procedure :: init_smatrix => beam_structure_init_smatrix
     procedure :: set_sentry => beam_structure_set_sentry
     procedure :: set_pol_f => beam_structure_set_pol_f
     procedure :: final_mom => beam_structure_final_mom
     procedure :: set_momentum => beam_structure_set_momentum
     procedure :: set_theta => beam_structure_set_theta
     procedure :: set_phi => beam_structure_set_phi
     procedure :: is_set => beam_structure_is_set
     procedure :: get_n_beam => beam_structure_get_n_beam
     procedure :: get_prt => beam_structure_get_prt
     procedure :: get_n_record => beam_structure_get_n_record
     procedure :: get_i_entry => beam_structure_get_i_entry
     procedure :: get_name => beam_structure_get_name
     procedure :: has_pdf => beam_structure_has_pdf
     procedure :: contains => beam_structure_contains
     procedure :: polarized => beam_structure_polarized
     procedure :: get_smatrix => beam_structure_get_smatrix
     procedure :: get_pol_f => beam_structure_get_pol_f
     procedure :: asymmetric => beam_structure_asymmetric
     procedure :: get_momenta => beam_structure_get_momenta
     procedure :: check_against_n_in => beam_structure_check_against_n_in
  end type beam_structure_t


  abstract interface
     function strfun_mode_fun (name) result (n)
       import
       type(string_t), intent(in) :: name
       integer :: n
     end function strfun_mode_fun
  end interface


  interface
    module function beam_structure_entry_to_string (object) result (string)
      class(beam_structure_entry_t), intent(in) :: object
      type(string_t) :: string
    end function beam_structure_entry_to_string
    module subroutine beam_structure_final_sf (object)
      class(beam_structure_t), intent(inout) :: object
    end subroutine beam_structure_final_sf
    module subroutine beam_structure_write (object, unit)
      class(beam_structure_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine beam_structure_write
    module function beam_structure_to_string (object, sf_only) result (string)
      class(beam_structure_t), intent(in) :: object
      logical, intent(in), optional :: sf_only
      type(string_t) :: string
    end function beam_structure_to_string
    module subroutine beam_structure_init_sf (beam_structure, prt, dim_array)
      class(beam_structure_t), intent(inout) :: beam_structure
      type(string_t), dimension(:), intent(in) :: prt
      integer, dimension(:), intent(in), optional :: dim_array
    end subroutine beam_structure_init_sf
    module subroutine beam_structure_set_sf (beam_structure, i, j, name)
      class(beam_structure_t), intent(inout) :: beam_structure
      integer, intent(in) :: i, j
      type(string_t), intent(in) :: name
    end subroutine beam_structure_set_sf
    module subroutine beam_structure_expand (beam_structure, strfun_mode)
      class(beam_structure_t), intent(inout) :: beam_structure
      procedure(strfun_mode_fun) :: strfun_mode
    end subroutine beam_structure_expand
    module subroutine beam_structure_final_pol (beam_structure)
      class(beam_structure_t), intent(inout) :: beam_structure
    end subroutine beam_structure_final_pol
    module subroutine beam_structure_init_pol (beam_structure, n)
      class(beam_structure_t), intent(inout) :: beam_structure
      integer, intent(in) :: n
    end subroutine beam_structure_init_pol
    elemental module function beam_structure_has_polarized_beams &
         (beam_structure) result (pol)
      logical :: pol
      class(beam_structure_t), intent(in) :: beam_structure
    end function beam_structure_has_polarized_beams
    module subroutine beam_structure_set_smatrix (beam_structure, i, smatrix)
      class(beam_structure_t), intent(inout) :: beam_structure
      integer, intent(in) :: i
      type(smatrix_t), intent(in) :: smatrix
    end subroutine beam_structure_set_smatrix
    module subroutine beam_structure_init_smatrix (beam_structure, i, n_entry)
      class(beam_structure_t), intent(inout) :: beam_structure
      integer, intent(in) :: i
      integer, intent(in) :: n_entry
    end subroutine beam_structure_init_smatrix
    module subroutine beam_structure_set_sentry &
         (beam_structure, i, i_entry, index, value)
      class(beam_structure_t), intent(inout) :: beam_structure
      integer, intent(in) :: i
      integer, intent(in) :: i_entry
      integer, dimension(:), intent(in) :: index
      complex(default), intent(in) :: value
    end subroutine beam_structure_set_sentry
    module subroutine beam_structure_set_pol_f (beam_structure, f)
      class(beam_structure_t), intent(inout) :: beam_structure
      real(default), dimension(:), intent(in) :: f
    end subroutine beam_structure_set_pol_f
    module subroutine beam_structure_final_mom (beam_structure)
      class(beam_structure_t), intent(inout) :: beam_structure
    end subroutine beam_structure_final_mom
    module subroutine beam_structure_set_momentum (beam_structure, p)
      class(beam_structure_t), intent(inout) :: beam_structure
      real(default), dimension(:), intent(in) :: p
    end subroutine beam_structure_set_momentum
    module subroutine beam_structure_set_theta (beam_structure, theta)
      class(beam_structure_t), intent(inout) :: beam_structure
      real(default), dimension(:), intent(in) :: theta
    end subroutine beam_structure_set_theta
    module subroutine beam_structure_set_phi (beam_structure, phi)
      class(beam_structure_t), intent(inout) :: beam_structure
      real(default), dimension(:), intent(in) :: phi
    end subroutine beam_structure_set_phi
    module function beam_structure_is_set (beam_structure) result (flag)
      class(beam_structure_t), intent(in) :: beam_structure
      logical :: flag
    end function beam_structure_is_set
    module function beam_structure_get_n_beam (beam_structure) result (n)
      class(beam_structure_t), intent(in) :: beam_structure
      integer :: n
    end function beam_structure_get_n_beam
    module function beam_structure_get_prt (beam_structure) result (prt)
      class(beam_structure_t), intent(in) :: beam_structure
      type(string_t), dimension(:), allocatable :: prt
    end function beam_structure_get_prt
    module function beam_structure_get_n_record (beam_structure) result (n)
      class(beam_structure_t), intent(in) :: beam_structure
      integer :: n
    end function beam_structure_get_n_record
    module function beam_structure_get_i_entry &
         (beam_structure, i) result (i_entry)
      class(beam_structure_t), intent(in) :: beam_structure
      integer, intent(in) :: i
      integer, dimension(:), allocatable :: i_entry
    end function beam_structure_get_i_entry
    module function beam_structure_get_name (beam_structure, i) result (name)
      type(string_t) :: name
      class(beam_structure_t), intent(in) :: beam_structure
      integer, intent(in) :: i
    end function beam_structure_get_name
    module function beam_structure_has_pdf (beam_structure) result (has_pdf)
      logical :: has_pdf
      class(beam_structure_t), intent(in) :: beam_structure
    end function beam_structure_has_pdf
    module function beam_structure_contains (beam_structure, name) result (flag)
      class(beam_structure_t), intent(in) :: beam_structure
      character(*), intent(in) :: name
      logical :: flag
    end function beam_structure_contains
    module function beam_structure_polarized (beam_structure) result (flag)
      class(beam_structure_t), intent(in) :: beam_structure
      logical :: flag
    end function beam_structure_polarized
    module function beam_structure_get_smatrix (beam_structure) result (smatrix)
      class(beam_structure_t), intent(in) :: beam_structure
      type(smatrix_t), dimension(:), allocatable :: smatrix
    end function beam_structure_get_smatrix
    module function beam_structure_get_pol_f (beam_structure) result (pol_f)
      class(beam_structure_t), intent(in) :: beam_structure
      real(default), dimension(:), allocatable :: pol_f
    end function beam_structure_get_pol_f
    module function beam_structure_asymmetric (beam_structure) result (flag)
      class(beam_structure_t), intent(in) :: beam_structure
      logical :: flag
    end function beam_structure_asymmetric
    module function beam_structure_get_momenta (beam_structure) result (p)
      class(beam_structure_t), intent(in) :: beam_structure
      type(vector3_t), dimension(:), allocatable :: p
    end function beam_structure_get_momenta
    module subroutine beam_structure_check_against_n_in &
         (beam_structure, n_in, applies)
      class(beam_structure_t), intent(in) :: beam_structure
      integer, intent(in) :: n_in
      logical, intent(out) :: applies
    end subroutine beam_structure_check_against_n_in
  end interface

end module beam_structures
