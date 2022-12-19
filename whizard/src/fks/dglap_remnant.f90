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

module dglap_remnant

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use phs_fks, only: isr_kinematics_t
  use fks_regions, only: region_data_t

  use nlo_data

  implicit none
  private

  public :: dglap_remnant_t

  type :: dglap_remnant_t
     type(nlo_settings_t), pointer :: settings => null ()
     type(region_data_t), pointer :: reg_data => null ()
     type(isr_kinematics_t), pointer :: isr_kinematics => null ()
     real(default) :: CA = 0, CF = 0, TR = 0
     real(default), dimension(:), allocatable :: sqme_born
     real(default), dimension(:,:), allocatable :: sf_factors
     real(default), dimension(:,:,:), allocatable :: sqme_color_c_extra
   contains
     procedure :: init => dglap_remnant_init
     procedure :: set_parameters => dglap_remnant_set_parameters
     procedure :: evaluate => dglap_remnant_evaluate
     procedure :: final => dglap_remnant_final
  end type dglap_remnant_t


  interface
    module subroutine dglap_remnant_init &
         (dglap, settings, reg_data, isr_kinematics)
      class(dglap_remnant_t), intent(inout) :: dglap
      type(nlo_settings_t), intent(in), target :: settings
      type(region_data_t), intent(in), target :: reg_data
      type(isr_kinematics_t), intent(in), target :: isr_kinematics
    end subroutine dglap_remnant_init
    module subroutine dglap_remnant_set_parameters (dglap, CA, CF, TR)
      class(dglap_remnant_t), intent(inout) :: dglap
      real(default), intent(in) :: CA, CF, TR
    end subroutine dglap_remnant_set_parameters
    module subroutine dglap_remnant_evaluate &
         (dglap, alpha_coupling, separate_uborns, sqme_dglap)
      class(dglap_remnant_t), intent(inout) :: dglap
      real(default), dimension(2), intent(in) :: alpha_coupling
      logical, intent(in) :: separate_uborns
      real(default), intent(inout), dimension(:) :: sqme_dglap
    end subroutine dglap_remnant_evaluate
    module subroutine dglap_remnant_final (dglap)
      class(dglap_remnant_t), intent(inout) :: dglap
    end subroutine dglap_remnant_final
  end interface

end module dglap_remnant

