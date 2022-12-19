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

module shower_pythia8

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use os_interface
  use lorentz
  use shower_base
  use particles
  use model_data
  use pdf
  use whizard_lha
  use pythia8
  use tauola_interface, only: taudec_settings_t

  implicit none
  private

  public :: shower_pythia8_t



  type, extends(shower_base_t) :: shower_pythia8_t
     private
     type(pythia8_t) :: pythia
     type(whizard_lha_t) :: lhaup
     logical :: user_process_set = .false.
     logical :: pythia_initialized = .false., &
          lhaup_initialized = .false.
   contains
     procedure :: init => shower_pythia8_init
     procedure, private :: set_user_process => shower_pythia8_set_user_process
     procedure :: import_particle_set => shower_pythia8_import_particle_set
     procedure :: generate_emissions => shower_pythia8_generate_emissions
     procedure :: make_particle_set => shower_pythia8_make_particle_set
     procedure, private :: transfer_settings => shower_pythia8_transfer_settings
     procedure :: get_final_colored_ME_momenta => &
          shower_pythia8_get_final_colored_ME_momenta
  end type shower_pythia8_t


  interface
    module subroutine shower_pythia8_init &
         (shower, settings, taudec_settings, pdf_data, os_data)
      class(shower_pythia8_t), intent(out) :: shower
      type(shower_settings_t), intent(in) :: settings
      type(taudec_settings_t), intent(in) :: taudec_settings
      type(pdf_data_t), intent(in) :: pdf_data
      type(os_data_t), intent(in) :: os_data
    end subroutine shower_pythia8_init
    module subroutine shower_pythia8_set_user_process (shower, pset)
      class(shower_pythia8_t), intent(inout) :: shower
      type(particle_set_t), intent(in) :: pset
    end subroutine shower_pythia8_set_user_process
    module subroutine shower_pythia8_import_particle_set &
         (shower, particle_set)
      class(shower_pythia8_t), target, intent(inout) :: shower
      type(particle_set_t), intent(in) :: particle_set
    end subroutine shower_pythia8_import_particle_set
    module subroutine shower_pythia8_generate_emissions &
          (shower, valid, number_of_emissions)
      class(shower_pythia8_t), intent(inout), target :: shower
      logical, intent(out) :: valid
      integer, optional, intent(in) :: number_of_emissions
    end subroutine shower_pythia8_generate_emissions
    module subroutine shower_pythia8_make_particle_set &
         (shower, particle_set, model, model_hadrons)
      class(shower_pythia8_t), intent(in) :: shower
      type(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model
      class(model_data_t), intent(in), target :: model_hadrons
    end subroutine shower_pythia8_make_particle_set
    module subroutine shower_pythia8_transfer_settings (shower)
      class(shower_pythia8_t), intent(inout), target :: shower
    end subroutine shower_pythia8_transfer_settings
    module subroutine shower_pythia8_get_final_colored_ME_momenta &
           (shower, momenta)
      class(shower_pythia8_t), intent(in) :: shower
      type(vector4_t), dimension(:), allocatable, intent(out) :: momenta
    end subroutine shower_pythia8_get_final_colored_ME_momenta
  end interface

end module shower_pythia8

