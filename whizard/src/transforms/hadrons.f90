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

module hadrons

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use event_transforms
  use lorentz
  use model_data
  use models
  use particles
  use physics_defs
  use process, only: process_t
  use instances, only: process_instance_t
  use process_stacks
  use whizard_lha
  use pythia8
  use rng_base, only: rng_t
  use shower_base
  use shower_pythia6
  use sm_qcd
  use variables

  implicit none
  private

  public :: HADRONS_UNDEFINED, HADRONS_WHIZARD, HADRONS_PYTHIA6, HADRONS_PYTHIA8
  public :: hadrons_method
  public :: hadron_settings_t
  public :: hadrons_hadrons_t
  public :: had_flav_t
  public :: lund_end
  public :: lund_pt_t
  public :: hadrons_pythia6_t
  public :: hadrons_pythia8_t
  public :: evt_hadrons_t

  integer, parameter :: HADRONS_UNDEFINED = 0
  integer, parameter :: HADRONS_WHIZARD = 1
  integer, parameter :: HADRONS_PYTHIA6 = 2
  integer, parameter :: HADRONS_PYTHIA8 = 3

  type :: hadron_settings_t
     logical :: active = .false.
     integer :: method = HADRONS_UNDEFINED
     real(default) :: enhanced_fraction = 0
     real(default) :: enhanced_width = 0
   contains
     procedure :: init => hadron_settings_init
     procedure :: write => hadron_settings_write
  end type hadron_settings_t

  type, abstract :: hadrons_t
     class(rng_t), allocatable :: rng
     type(shower_settings_t) :: shower_settings
     type(hadron_settings_t) :: hadron_settings
     type(model_t), pointer :: model => null()
   contains
     procedure (hadrons_init), deferred :: init
     procedure (hadrons_hadronize), deferred :: hadronize
     procedure (hadrons_make_particle_set), deferred :: make_particle_set
     procedure :: import_rng => hadrons_import_rng
  end type hadrons_t

  type, extends (hadrons_t) :: hadrons_hadrons_t
     contains
         procedure :: init => hadrons_hadrons_init
         procedure :: hadronize => hadrons_hadrons_hadronize
         procedure :: make_particle_set => hadrons_hadrons_make_particle_set
    end type hadrons_hadrons_t

  type had_flav_t
  end type had_flav_t

  type lund_end
     logical :: from_pos
     integer :: i_end
     integer :: i_max
     integer :: id_had
     integer :: i_pos_old
     integer :: i_neg_old
     integer :: i_pos_new
     integer :: i_neg_new
     real(default) :: px_old
     real(default) :: py_old
     real(default) :: px_new
     real(default) :: py_new
     real(default) :: px_had
     real(default) :: py_had
     real(default) :: m_had
     real(default) :: mT2_had
     real(default) :: z_had
     real(default) :: gamma_old
     real(default) :: gamma_new
     real(default) :: x_pos_old
     real(default) :: x_pos_new
     real(default) :: x_pos_had
     real(default) :: x_neg_old
     real(default) :: x_neg_new
     real(default) :: x_neg_had
     type(had_flav_t) :: old_flav
     type(had_flav_t) :: new_flav
     type(vector4_t) :: p_had
     type(vector4_t) :: p_pre
  end type lund_end

  type lund_pt_t
     real(default) :: sigma_min
     real(default) :: sigma_q
     real(default) :: enhanced_frac
     real(default) :: enhanced_width
     real(default) :: sigma_to_had
     class(rng_t), allocatable :: rng
   contains
     procedure :: init => lund_pt_init
  end type lund_pt_t

  type, extends (hadrons_t) :: hadrons_pythia6_t
   contains
     procedure :: init => hadrons_pythia6_init
     procedure :: hadronize => hadrons_pythia6_hadronize
     procedure :: make_particle_set => hadrons_pythia6_make_particle_set
  end type hadrons_pythia6_t

  type, extends (hadrons_t) :: hadrons_pythia8_t
     type(pythia8_t) :: pythia
     type(whizard_lha_t) :: lhaup
     logical :: user_process_set = .false.
     logical :: pythia_initialized = .false., &
          lhaup_initialized = .false.
  contains
     procedure :: init => hadrons_pythia8_init
     procedure, private :: transfer_settings => hadrons_pythia8_transfer_settings
     procedure, private :: set_user_process => hadrons_pythia8_set_user_process
     procedure, private :: import_particle_set => &
          hadrons_pythia8_import_particle_set
     procedure :: hadronize => hadrons_pythia8_hadronize
     procedure :: make_particle_set => hadrons_pythia8_make_particle_set
  end type hadrons_pythia8_t

  type, extends (evt_t) :: evt_hadrons_t
     class(hadrons_t), allocatable :: hadrons
     type(model_t), pointer :: model_hadrons => null()
     type(qcd_t) :: qcd
     logical :: is_first_event
   contains
     procedure :: init => evt_hadrons_init
     procedure :: write_name => evt_hadrons_write_name
     procedure :: write => evt_hadrons_write
     procedure :: first_event => evt_hadrons_first_event
     procedure :: generate_weighted => evt_hadrons_generate_weighted
     procedure :: make_particle_set => evt_hadrons_make_particle_set
     procedure :: connect => evt_hadrons_connect
     procedure :: make_rng => evt_hadrons_make_rng
     procedure :: prepare_new_event => evt_hadrons_prepare_new_event
  end type evt_hadrons_t


  interface hadrons_method
     module procedure hadrons_method_of_string
     module procedure hadrons_method_to_string
  end interface
  abstract interface
    subroutine hadrons_init &
         (hadrons, shower_settings, hadron_settings, model_hadrons)
      import
      class(hadrons_t), intent(out) :: hadrons
      type(shower_settings_t), intent(in) :: shower_settings
      type(hadron_settings_t), intent(in) :: hadron_settings
      type(model_t), target, intent(in) :: model_hadrons
    end subroutine hadrons_init
   end interface

  abstract interface
     subroutine hadrons_hadronize (hadrons, particle_set, valid)
       import
       class(hadrons_t), intent(inout) :: hadrons
       type(particle_set_t), intent(in) :: particle_set
       logical, intent(out) :: valid
     end subroutine hadrons_hadronize
  end interface
  abstract interface
     subroutine hadrons_make_particle_set (hadrons, particle_set, &
          model, valid)
       import
       class(hadrons_t), intent(in) :: hadrons
       type(particle_set_t), intent(inout) :: particle_set
       class(model_data_t), intent(in), target :: model
       logical, intent(out) :: valid
     end subroutine hadrons_make_particle_set
  end interface


  interface
    elemental module function hadrons_method_of_string (string) result (i)
      integer :: i
      type(string_t), intent(in) :: string
    end function hadrons_method_of_string
    elemental module function hadrons_method_to_string (i) result (string)
      type(string_t) :: string
      integer, intent(in) :: i
    end function hadrons_method_to_string
    module subroutine hadron_settings_init (hadron_settings, var_list)
      class(hadron_settings_t), intent(out) :: hadron_settings
      type(var_list_t), intent(in) :: var_list
    end subroutine hadron_settings_init
    module subroutine hadron_settings_write (settings, unit)
      class(hadron_settings_t), intent(in) :: settings
      integer, intent(in), optional :: unit
    end subroutine hadron_settings_write
    pure module subroutine hadrons_import_rng (hadrons, rng)
      class(hadrons_t), intent(inout) :: hadrons
      class(rng_t), intent(inout), allocatable :: rng
    end subroutine hadrons_import_rng
    module subroutine hadrons_hadrons_init &
         (hadrons, shower_settings, hadron_settings, model_hadrons)
      class(hadrons_hadrons_t), intent(out) :: hadrons
      type(shower_settings_t), intent(in) :: shower_settings
      type(hadron_settings_t), intent(in) :: hadron_settings
      type(model_t), intent(in), target :: model_hadrons
    end subroutine hadrons_hadrons_init
    module subroutine hadrons_hadrons_hadronize (hadrons, particle_set, valid)
      class(hadrons_hadrons_t), intent(inout) :: hadrons
      type(particle_set_t), intent(in) :: particle_set
      logical, intent(out) :: valid
    end subroutine hadrons_hadrons_hadronize
    module subroutine lund_pt_init (lund_pt, settings)
      class (lund_pt_t), intent(out) :: lund_pt
      type(hadron_settings_t), intent(in) :: settings
    end subroutine lund_pt_init
    module subroutine hadrons_hadrons_make_particle_set &
           (hadrons, particle_set, model, valid)
      class(hadrons_hadrons_t), intent(in) :: hadrons
      type(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model
      logical, intent(out) :: valid
    end subroutine hadrons_hadrons_make_particle_set
    module subroutine hadrons_pythia6_init &
         (hadrons, shower_settings, hadron_settings, model_hadrons)
      class(hadrons_pythia6_t), intent(out) :: hadrons
      type(shower_settings_t), intent(in) :: shower_settings
      type(hadron_settings_t), intent(in) :: hadron_settings
      type(model_t), intent(in), target :: model_hadrons
    end subroutine hadrons_pythia6_init
    module subroutine hadrons_pythia6_hadronize (hadrons, particle_set, valid)
      class(hadrons_pythia6_t), intent(inout) :: hadrons
      type(particle_set_t), intent(in) :: particle_set
      logical, intent(out) :: valid
    end subroutine hadrons_pythia6_hadronize
    module subroutine hadrons_pythia6_make_particle_set &
         (hadrons, particle_set, model, valid)
      class(hadrons_pythia6_t), intent(in) :: hadrons
      type(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model
      logical, intent(out) :: valid
    end subroutine hadrons_pythia6_make_particle_set
    module subroutine hadrons_pythia8_init &
         (hadrons, shower_settings, hadron_settings, model_hadrons)
      class(hadrons_pythia8_t), intent(out) :: hadrons
      type(shower_settings_t), intent(in) :: shower_settings
      type(hadron_settings_t), intent(in) :: hadron_settings
      type(model_t), intent(in), target :: model_hadrons
    end subroutine hadrons_pythia8_init
    module subroutine hadrons_pythia8_transfer_settings (hadrons)
      class(hadrons_pythia8_t), intent(inout), target :: hadrons
    end subroutine hadrons_pythia8_transfer_settings
    module subroutine hadrons_pythia8_set_user_process (hadrons, pset)
      class(hadrons_pythia8_t), intent(inout) :: hadrons
      type(particle_set_t), intent(in) :: pset
    end subroutine hadrons_pythia8_set_user_process
    module subroutine hadrons_pythia8_import_particle_set &
         (hadrons, particle_set)
      class(hadrons_pythia8_t), target, intent(inout) :: hadrons
      type(particle_set_t), intent(in) :: particle_set
    end subroutine hadrons_pythia8_import_particle_set
    module subroutine hadrons_pythia8_hadronize (hadrons, particle_set, valid)
      class(hadrons_pythia8_t), intent(inout) :: hadrons
      type(particle_set_t), intent(in) :: particle_set
      logical, intent(out) :: valid
    end subroutine hadrons_pythia8_hadronize
    module subroutine hadrons_pythia8_make_particle_set &
           (hadrons, particle_set, model, valid)
      class(hadrons_pythia8_t), intent(in) :: hadrons
      type(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model
      logical, intent(out) :: valid
    end subroutine hadrons_pythia8_make_particle_set
    module subroutine evt_hadrons_init (evt, model_hadrons)
      class(evt_hadrons_t), intent(out) :: evt
      type(model_t), intent(in), target :: model_hadrons
    end subroutine evt_hadrons_init
    module subroutine evt_hadrons_write_name (evt, unit)
      class(evt_hadrons_t), intent(in) :: evt
      integer, intent(in), optional :: unit
    end subroutine evt_hadrons_write_name
    module subroutine evt_hadrons_write &
         (evt, unit, verbose, more_verbose, testflag)
      class(evt_hadrons_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, more_verbose, testflag
    end subroutine evt_hadrons_write
    module subroutine evt_hadrons_first_event (evt)
      class(evt_hadrons_t), intent(inout) :: evt
    end subroutine evt_hadrons_first_event
    module subroutine evt_hadrons_generate_weighted (evt, probability)
      class(evt_hadrons_t), intent(inout) :: evt
      real(default), intent(inout) :: probability
    end subroutine evt_hadrons_generate_weighted
    module subroutine evt_hadrons_make_particle_set &
         (evt, factorization_mode, keep_correlations, r)
      class(evt_hadrons_t), intent(inout) :: evt
      integer, intent(in) :: factorization_mode
      logical, intent(in) :: keep_correlations
      real(default), dimension(:), intent(in), optional :: r
    end subroutine evt_hadrons_make_particle_set
    module subroutine evt_hadrons_connect &
         (evt, process_instance, model, process_stack)
      class(evt_hadrons_t), intent(inout), target :: evt
      type(process_instance_t), intent(in), target :: process_instance
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
    end subroutine evt_hadrons_connect  
    module subroutine evt_hadrons_make_rng (evt, process)
      class(evt_hadrons_t), intent(inout) :: evt
      type(process_t), intent(inout) :: process
      class(rng_t), allocatable :: rng
    end subroutine evt_hadrons_make_rng
    module subroutine evt_hadrons_prepare_new_event (evt, i_mci, i_term)
      class(evt_hadrons_t), intent(inout) :: evt
      integer, intent(in) :: i_mci, i_term
    end subroutine evt_hadrons_prepare_new_event
  end interface

end module hadrons
