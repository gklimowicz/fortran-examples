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

module shower_base

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use particles
  use os_interface
  use rng_base
  use variables
  use model_data
  use pdf
  use tauola_interface

  implicit none
  private

  public :: PS_WHIZARD, PS_PYTHIA6, PS_PYTHIA8, PS_UNDEFINED
  public :: shower_method_of_string
  public :: shower_method_to_string
  public :: shower_settings_t
  public :: shower_base_t
  public :: D_alpha_s_isr
  public :: D_alpha_s_fsr
  public :: mass_type
  public :: mass_squared_type
  public :: number_of_flavors

  integer, parameter :: PS_UNDEFINED = 0
  integer, parameter :: PS_WHIZARD = 1
  integer, parameter :: PS_PYTHIA6 = 2
  integer, parameter :: PS_PYTHIA8 = 3
  real(default), public :: D_min_scale = 0.5_default
  logical, public :: treat_light_quarks_massless = .true.
  logical, public :: treat_duscb_quarks_massless = .false.
  real(default), public :: scalefactor1 = 0.02_default
  real(default), public :: scalefactor2 = 0.02_default

  type :: shower_settings_t
     logical :: active = .false.
     logical :: isr_active = .false.
     logical :: fsr_active = .false.
     logical :: muli_active = .false.
     logical :: hadronization_active = .false.
     logical :: tau_dec = .false.
     logical :: verbose = .false.
     integer :: method = PS_UNDEFINED
     logical :: hadron_collision = .false.
     logical :: mlm_matching = .false.
     logical :: ckkw_matching = .false.
     logical :: powheg_matching = .false.
     type(string_t) :: pythia6_pygive
     type(string_t) :: pythia8_config
     type(string_t) :: pythia8_config_file
     real(default) :: min_virtuality = 1._default   ! PARJ(82)^2
     real(default) :: fsr_lambda = 0.29_default     ! PARP(72)
     real(default) :: isr_lambda = 0.29_default     ! PARP(61)
     integer :: max_n_flavors = 5                   ! MSTJ(45)
     logical :: isr_alphas_running = .true.        ! MSTP(64)
     logical :: fsr_alphas_running = .true.        ! MSTJ(44)
     real(default) :: fixed_alpha_s = 0.2_default   ! PARU(111)
     logical :: alpha_s_fudged = .true.
     logical :: isr_pt_ordered = .false.
     logical :: isr_angular_ordered = .true.        ! MSTP(62)
     real(default) :: isr_primordial_kt_width = 1.5_default  ! PARP(91)
     real(default) :: isr_primordial_kt_cutoff = 5._default  ! PARP(93)
     real(default) :: isr_z_cutoff = 0.999_default  ! 1-PARP(66)
     real(default) :: isr_minenergy = 2._default    ! PARP(65)
     real(default) :: isr_tscalefactor = 1._default
     logical :: isr_only_onshell_emitted_partons = .true.   ! MSTP(63)
   contains
     procedure :: init => shower_settings_init
     procedure :: write => shower_settings_write
  end type shower_settings_t

  type, abstract :: shower_base_t
     class(rng_t), allocatable :: rng
     type(string_t) :: name
     type(pdf_data_t) :: pdf_data
     type(shower_settings_t) :: settings
     type(taudec_settings_t) :: taudec_settings
     type(os_data_t) :: os_data
     real(default) :: fac_scale
     real(default) :: alpha_s
   contains
     procedure :: write_msg => shower_base_write_msg
     procedure :: import_rng => shower_base_import_rng
     procedure (shower_base_init), deferred :: init
     procedure :: prepare_new_event => shower_base_prepare_new_event
     procedure (shower_base_import_particle_set), deferred :: import_particle_set
     procedure (shower_base_generate_emissions), deferred :: generate_emissions
     procedure (shower_base_make_particle_set), deferred :: make_particle_set
     procedure (shower_base_get_final_colored_ME_momenta), deferred :: &
          get_final_colored_ME_momenta
  end type shower_base_t


  abstract interface
    subroutine shower_base_init (shower, settings, taudec_settings, pdf_data, os_data)
      import
      class(shower_base_t), intent(out) :: shower
      type(shower_settings_t), intent(in) :: settings
      type(taudec_settings_t), intent(in) :: taudec_settings
      type(pdf_data_t), intent(in) :: pdf_data
      type(os_data_t), intent(in) :: os_data
    end subroutine shower_base_init
   end interface

  abstract interface
     subroutine shower_base_import_particle_set &
            (shower, particle_set)
       import
       class(shower_base_t), target, intent(inout) :: shower
       type(particle_set_t), intent(in) :: particle_set
     end subroutine shower_base_import_particle_set
  end interface

  abstract interface
     subroutine shower_base_generate_emissions &
            (shower, valid, number_of_emissions)
      import
      class(shower_base_t), intent(inout), target :: shower
       logical, intent(out) :: valid
      integer, optional, intent(in) :: number_of_emissions
    end subroutine shower_base_generate_emissions
   end interface

  abstract interface
     subroutine shower_base_make_particle_set &
         (shower, particle_set, model, model_hadrons)
       import
       class(shower_base_t), intent(in) :: shower
       type(particle_set_t), intent(inout) :: particle_set
       class(model_data_t), intent(in), target :: model
       class(model_data_t), intent(in), target :: model_hadrons
     end subroutine shower_base_make_particle_set
  end interface

  abstract interface
     subroutine shower_base_get_final_colored_ME_momenta &
            (shower, momenta)
       import
       class(shower_base_t), intent(in) :: shower
       type(vector4_t), dimension(:), allocatable, intent(out) :: momenta
     end subroutine shower_base_get_final_colored_ME_momenta
  end interface


  interface
    elemental module function shower_method_of_string (string) result (i)
      integer :: i
      type(string_t), intent(in) :: string
    end function shower_method_of_string
    elemental module function shower_method_to_string (i) result (string)
      type(string_t) :: string
      integer, intent(in) :: i
    end function shower_method_to_string
    module subroutine shower_settings_init (settings, var_list)
      class(shower_settings_t), intent(out) :: settings
      type(var_list_t), intent(in) :: var_list
    end subroutine shower_settings_init
    module subroutine shower_settings_write (settings, unit)
      class(shower_settings_t), intent(in) :: settings
      integer, intent(in), optional :: unit
    end subroutine shower_settings_write
    module subroutine shower_base_write_msg (shower)
      class(shower_base_t), intent(inout) :: shower
    end subroutine shower_base_write_msg
    pure module subroutine shower_base_import_rng (shower, rng)
      class(shower_base_t), intent(inout) :: shower
      class(rng_t), intent(inout), allocatable :: rng
    end subroutine shower_base_import_rng
    module subroutine shower_base_prepare_new_event (shower, fac_scale, alpha_s)
      class(shower_base_t), intent(inout) :: shower
      real(default), intent(in) :: fac_scale, alpha_s
    end subroutine shower_base_prepare_new_event
    module function D_alpha_s_isr (tin, settings) result (alpha_s)
      real(default), intent(in) :: tin
      type(shower_settings_t), intent(in) :: settings
      real(default) :: alpha_s
    end function D_alpha_s_isr
    module function D_alpha_s_fsr (tin, settings) result (alpha_s)
      real(default), intent(in) :: tin
      type(shower_settings_t), intent(in) :: settings
      real(default) :: alpha_s
    end function D_alpha_s_fsr
    elemental module function mass_type (type, m2_default) result (mass)
      integer, intent(in) :: type
      real(default), intent(in) :: m2_default
      real(default) :: mass
    end function mass_type
    elemental module function mass_squared_type &
         (type, m2_default) result (mass2)
      integer, intent(in) :: type
      real(default), intent(in) :: m2_default
      real(default) :: mass2
    end function mass_squared_type
    elemental module function number_of_flavors &
         (t, d_nf, min_virtuality) result (nr)
      real(default), intent(in) :: t, min_virtuality
      integer, intent(in) :: d_nf
      real(default) :: nr
    end function number_of_flavors
  end interface

end module shower_base
