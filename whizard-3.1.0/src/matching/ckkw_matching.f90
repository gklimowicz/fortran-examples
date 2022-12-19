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

module ckkw_matching

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use constants
  use lorentz
  use particles
  use rng_base
  use shower_base
  use shower_partons
  use variables
  use matching_base

  implicit none
  private

  public :: ckkw_matching_settings_t
  public :: ckkw_pseudo_shower_weights_t
  public :: ckkw_matching_t
  public :: ckkw_matching_apply

  type :: ckkw_matching_settings_t
     real(default) :: alphaS = 0.118_default
     real(default) :: Qmin = one
     integer :: n_max_jets = 0
   contains
     procedure :: init => ckkw_matching_settings_init
     procedure :: write => ckkw_matching_settings_write
  end type ckkw_matching_settings_t

  type :: ckkw_pseudo_shower_weights_t
     real(default) :: alphaS
     real(default), dimension(:), allocatable :: weights
     real(default), dimension(:,:), allocatable :: weights_by_type
   contains
    procedure :: init => ckkw_pseudo_shower_weights_init
    procedure :: write => ckkw_pseudo_shower_weights_write
    procedure :: fake => ckkw_pseudo_shower_weights_fake
  end type ckkw_pseudo_shower_weights_t

  type, extends (matching_t) :: ckkw_matching_t
     type(ckkw_matching_settings_t) :: settings
     type(ckkw_pseudo_shower_weights_t) :: weights
   contains
       procedure :: init => ckkw_matching_init
       procedure :: write => ckkw_matching_write
       procedure :: get_method => ckkw_matching_get_method
       procedure :: before_shower => ckkw_matching_before_shower
       procedure :: after_shower => ckkw_matching_after_shower
  end type ckkw_matching_t


  interface
    module subroutine ckkw_matching_settings_init (settings, var_list)
      class(ckkw_matching_settings_t), intent(out) :: settings
      type(var_list_t), intent(in) :: var_list
    end subroutine ckkw_matching_settings_init
    module subroutine ckkw_matching_settings_write (settings, unit)
      class(ckkw_matching_settings_t), intent(in) :: settings
      integer, intent(in), optional :: unit
    end subroutine ckkw_matching_settings_write
    module subroutine ckkw_pseudo_shower_weights_init (weights)
      class(ckkw_pseudo_shower_weights_t), intent(out) :: weights
    end subroutine ckkw_pseudo_shower_weights_init
    module subroutine ckkw_pseudo_shower_weights_write (weights, unit)
      class(ckkw_pseudo_shower_weights_t), intent(in) :: weights
      integer, intent(in), optional :: unit
    end subroutine ckkw_pseudo_shower_weights_write
    pure module subroutine ckkw_pseudo_shower_weights_fake &
         (weights, particle_set)
      class(ckkw_pseudo_shower_weights_t), intent(inout) :: weights
      type(particle_set_t), intent(in) :: particle_set
    end subroutine ckkw_pseudo_shower_weights_fake
    module subroutine ckkw_matching_init (matching, var_list, process_name)
      class(ckkw_matching_t), intent(out) :: matching
      type(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: process_name
    end subroutine ckkw_matching_init
    module subroutine ckkw_matching_write (matching, unit)
      class(ckkw_matching_t), intent(in) :: matching
      integer, intent(in), optional :: unit
    end subroutine ckkw_matching_write
    module function ckkw_matching_get_method (matching) result (method)
      type(string_t) :: method
      class(ckkw_matching_t), intent(in) :: matching
    end function ckkw_matching_get_method
    module subroutine ckkw_matching_before_shower &
         (matching, particle_set, vetoed)
      class(ckkw_matching_t), intent(inout) :: matching
      type(particle_set_t), intent(inout) :: particle_set
      logical, intent(out) :: vetoed
    end subroutine ckkw_matching_before_shower
    module subroutine ckkw_matching_apply &
         (partons, settings, weights, rng, vetoed)
      type(parton_pointer_t), dimension(:), intent(inout), allocatable :: &
           partons
      type(ckkw_matching_settings_t), intent(in) :: settings
      type(ckkw_pseudo_shower_weights_t), intent(in) :: weights
      class(rng_t), intent(inout), allocatable :: rng
      logical, intent(out) :: vetoed
    end subroutine ckkw_matching_apply
    module subroutine ckkw_matching_after_shower &
         (matching, particle_set, vetoed)
      class(ckkw_matching_t), intent(inout) :: matching
      type(particle_set_t), intent(inout) :: particle_set
      logical, intent(out) :: vetoed
    end subroutine ckkw_matching_after_shower
  end interface

end module ckkw_matching
