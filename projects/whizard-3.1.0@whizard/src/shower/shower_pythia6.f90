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

module shower_pythia6

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use os_interface
  use lorentz
  use shower_base
  use particles
  use model_data
  use pdf
  use tauola_interface

  implicit none
  private

  public :: shower_pythia6_t
  public :: pythia6_combine_with_particle_set
  public :: pythia6_setup_lhe_io_units
  public :: pythia6_set_config
  public :: pythia6_set_error
  public :: pythia6_get_error
  public :: pythia6_tauola_active
  public :: pythia6_handle_errors
  public :: pythia6_set_verbose
  public :: pythia6_set_last_treated_line

  type, extends (shower_base_t) :: shower_pythia6_t
     integer :: initialized_for_NPRUP = 0
     logical :: warning_given = .false.
   contains
       procedure :: init => shower_pythia6_init
       procedure :: import_particle_set => shower_pythia6_import_particle_set
       procedure :: generate_emissions => shower_pythia6_generate_emissions
       procedure :: make_particle_set => shower_pythia6_make_particle_set
       procedure :: transfer_settings => shower_pythia6_transfer_settings
       procedure :: combine_with_particle_set => &
            shower_pythia6_combine_with_particle_set
       procedure :: get_final_colored_ME_momenta => &
            shower_pythia6_get_final_colored_ME_momenta
  end type shower_pythia6_t


  interface
    module subroutine shower_pythia6_init &
         (shower, settings, taudec_settings, pdf_data, os_data)
      class(shower_pythia6_t), intent(out) :: shower
      type(shower_settings_t), intent(in) :: settings
      type(taudec_settings_t), intent(in) :: taudec_settings
      type(pdf_data_t), intent(in) :: pdf_data
      type(os_data_t), intent(in) :: os_data
    end subroutine shower_pythia6_init
    module subroutine shower_pythia6_import_particle_set &
         (shower, particle_set)
      class(shower_pythia6_t), target, intent(inout) :: shower
      type(particle_set_t), intent(in) :: particle_set
    end subroutine shower_pythia6_import_particle_set
    module subroutine shower_pythia6_generate_emissions &
         (shower, valid, number_of_emissions)
      class(shower_pythia6_t), intent(inout), target :: shower
      logical, intent(out) :: valid
      integer, optional, intent(in) :: number_of_emissions
    end subroutine shower_pythia6_generate_emissions
    module subroutine shower_pythia6_make_particle_set &
         (shower, particle_set, model, model_hadrons)
      class(shower_pythia6_t), intent(in) :: shower
      type(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model
      class(model_data_t), intent(in), target :: model_hadrons
    end subroutine shower_pythia6_make_particle_set
    module subroutine shower_pythia6_transfer_settings (shower)
      class(shower_pythia6_t), intent(inout) :: shower
    end subroutine shower_pythia6_transfer_settings
    module subroutine shower_pythia6_combine_with_particle_set &
           (shower, particle_set, model_in, model_hadrons)
      class(shower_pythia6_t), intent(in) :: shower
      type(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model_in
      class(model_data_t), intent(in), target :: model_hadrons
    end subroutine shower_pythia6_combine_with_particle_set
    module subroutine pythia6_combine_with_particle_set (particle_set, model_in, &
         model_hadrons, settings)
      type(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model_in
      class(model_data_t), intent(in), target :: model_hadrons
      type(shower_settings_t), intent(in) :: settings
    end subroutine pythia6_combine_with_particle_set
    module subroutine shower_pythia6_get_final_colored_ME_momenta &
           (shower, momenta)
      class(shower_pythia6_t), intent(in) :: shower
      type(vector4_t), dimension(:), allocatable, intent(out) :: momenta
    end subroutine shower_pythia6_get_final_colored_ME_momenta
    module subroutine pythia6_setup_lhe_io_units (u_W2P, u_P2W)
      integer, intent(out) :: u_W2P
      integer, intent(out), optional :: u_P2W
    end subroutine pythia6_setup_lhe_io_units
    module subroutine pythia6_set_config (pygive_all)
      type(string_t), intent(in) :: pygive_all
    end subroutine pythia6_set_config
    module subroutine pythia6_set_error (mstu23)
      integer, intent(in) :: mstu23
    end subroutine pythia6_set_error
    module function pythia6_get_error () result (mstu23)
      integer :: mstu23
    end function pythia6_get_error
    module function pythia6_tauola_active () result (active)
      logical :: active
    end function pythia6_tauola_active
    module function pythia6_handle_errors () result (valid)
      logical :: valid
    end function pythia6_handle_errors
    module subroutine pythia6_set_verbose (verbose)
      logical, intent(in) :: verbose
    end subroutine pythia6_set_verbose
    module subroutine pythia6_set_last_treated_line (last_line)
      integer,intent(in) :: last_line
    end subroutine pythia6_set_last_treated_line
  end interface

end module shower_pythia6

