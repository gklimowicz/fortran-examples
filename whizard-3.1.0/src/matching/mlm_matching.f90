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

module mlm_matching

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use constants
  use lorentz
  use particles
  use variables
  use matching_base

  implicit none
  private

  public :: mlm_matching_settings_t
  public :: mlm_matching_t

  type :: mlm_matching_settings_t
     real(default) :: mlm_Qcut_ME = one
     real(default) :: mlm_Qcut_PS = one
     real(default) :: mlm_ptmin, mlm_etamax, mlm_Rmin, mlm_Emin
     real(default) :: mlm_ETclusfactor = 0.2_default
     real(default) :: mlm_ETclusminE = five
     real(default) :: mlm_etaclusfactor = one
     real(default) :: mlm_Rclusfactor = one
     real(default) :: mlm_Eclusfactor = one
     integer :: kt_imode_hadronic = 4313
     integer :: kt_imode_leptonic = 1111
     integer :: mlm_nmaxMEjets = 0
   contains
     procedure :: init => mlm_matching_settings_init
     procedure :: write => mlm_matching_settings_write
  end type mlm_matching_settings_t

  type, extends (matching_t) :: mlm_matching_t
     type(vector4_t), dimension(:), allocatable, public :: P_ME
     type(vector4_t), dimension(:), allocatable, public :: P_PS
     type(vector4_t), dimension(:), allocatable, private :: JETS_ME
     type(vector4_t), dimension(:), allocatable, private :: JETS_PS
     type(mlm_matching_settings_t) :: settings
   contains
     procedure :: init => mlm_matching_init
     procedure :: write => mlm_matching_write
     procedure :: get_method => mlm_matching_get_method
     procedure :: before_shower => mlm_matching_before_shower
     procedure :: after_shower => mlm_matching_after_shower
     procedure :: fill_P_PS => mlm_matching_fill_P_PS
     procedure :: apply => mlm_matching_apply
  end type mlm_matching_t


  interface
    module subroutine mlm_matching_settings_init (settings, var_list)
      class(mlm_matching_settings_t), intent(out) :: settings
      type(var_list_t), intent(in) :: var_list
    end subroutine mlm_matching_settings_init
    module subroutine mlm_matching_settings_write (settings, unit)
      class(mlm_matching_settings_t), intent(in) :: settings
      integer, intent(in), optional :: unit
    end subroutine mlm_matching_settings_write
  module subroutine mlm_matching_init (matching, var_list, process_name)
    class(mlm_matching_t), intent(out) :: matching
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_name
  end subroutine mlm_matching_init
    module subroutine mlm_matching_write (matching, unit)
      class(mlm_matching_t), intent(in) :: matching
      integer, intent(in), optional :: unit
    end subroutine mlm_matching_write
    module function mlm_matching_get_method (matching) result (method)
       type(string_t) :: method
       class(mlm_matching_t), intent(in) :: matching
    end function mlm_matching_get_method
    module subroutine mlm_matching_before_shower &
         (matching, particle_set, vetoed)
      class(mlm_matching_t), intent(inout) :: matching
      type(particle_set_t), intent(inout) :: particle_set
      logical, intent(out) :: vetoed
    end subroutine mlm_matching_before_shower
    module subroutine mlm_matching_after_shower (matching, particle_set, vetoed)
      class(mlm_matching_t), intent(inout) :: matching
      type(particle_set_t), intent(inout) :: particle_set
      logical, intent(out) :: vetoed
    end subroutine mlm_matching_after_shower
    module subroutine mlm_matching_fill_P_PS (matching, particle_set)
      class(mlm_matching_t), intent(inout) :: matching
      type(particle_set_t), intent(in) :: particle_set
    end subroutine mlm_matching_fill_P_PS
    module subroutine mlm_matching_apply (matching, vetoed)
      class(mlm_matching_t), intent(inout) :: matching
      logical, intent(out) :: vetoed
    end subroutine mlm_matching_apply
  end interface

end module mlm_matching
