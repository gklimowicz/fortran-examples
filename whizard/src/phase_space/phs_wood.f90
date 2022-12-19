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

module phs_wood

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use lorentz
  use model_data
  use flavors
  use phs_base
  use mappings
  use resonances, only: resonance_history_set_t
  use phs_forests
  use cascades
  use cascades2

  implicit none
  private

  public :: phs_wood_config_t
  public :: phs_wood_t

  type, extends (phs_config_t) :: phs_wood_config_t
     character(32) :: md5sum_forest = ""
     type(string_t) :: phs_path
     integer :: io_unit = 0
     logical :: io_unit_keep_open = .false.
     logical :: use_equivalences = .false.
     logical :: fatal_beam_decay = .true.
     type(mapping_defaults_t) :: mapping_defaults
     type(phs_parameters_t) :: par
     type(string_t) :: run_id
     type(cascade_set_t), allocatable :: cascade_set
     logical :: use_cascades2 = .false.
     type(feyngraph_set_t), allocatable :: feyngraph_set
     type(phs_forest_t) :: forest
     type(os_data_t) :: os_data
     logical :: is_combined_integration = .false.
   contains
     procedure :: final => phs_wood_config_final
     procedure :: increase_n_par => phs_wood_config_increase_n_par
     procedure :: write => phs_wood_config_write
     procedure :: write_forest => phs_wood_config_write_forest
     procedure :: set_parameters => phs_wood_config_set_parameters
     procedure :: enable_equivalences => phs_wood_config_enable_equivalences
     procedure :: set_mapping_defaults => phs_wood_config_set_mapping_defaults
     procedure :: set_input => phs_wood_config_set_input
     procedure :: generate_phase_space => phs_wood_config_generate_phase_space
     procedure :: write_phase_space => phs_wood_config_write_phase_space
     procedure :: clear_phase_space => phs_wood_config_clear_phase_space
     procedure :: extract_resonance_history_set &
          => phs_wood_config_extract_resonance_history_set
     procedure :: configure => phs_wood_config_configure
     procedure :: compute_md5sum_forest => phs_wood_config_compute_md5sum_forest
     procedure :: make_phs_filename => phs_wood_make_phs_filename
     procedure :: reshuffle_flavors => phs_wood_config_reshuffle_flavors
     procedure :: set_momentum_links => phs_wood_config_set_momentum_links
     procedure :: record_s_mappings => phs_wood_config_record_s_mappings
     procedure :: record_on_shell => phs_wood_config_record_on_shell
     procedure :: get_md5sum => phs_wood_config_get_md5sum
     procedure :: read_phs_file => phs_wood_read_phs_file
     procedure :: startup_message => phs_wood_config_startup_message
     procedure, nopass :: allocate_instance => phs_wood_config_allocate_instance
  end type phs_wood_config_t

  type, extends (phs_t) :: phs_wood_t
     real(default) :: sqrts = 0
     type(phs_forest_t) :: forest
     real(default), dimension(3) :: r_real
     integer :: n_r_born = 0
   contains
     procedure :: write => phs_wood_write
     procedure :: write_forest => phs_wood_write_forest
     procedure :: final => phs_wood_final
     procedure :: init => phs_wood_init
     procedure :: evaluate_selected_channel => phs_wood_evaluate_selected_channel
     procedure :: evaluate_other_channels => phs_wood_evaluate_other_channels
     procedure :: inverse => phs_wood_inverse
  end type phs_wood_t


  interface
    module subroutine phs_wood_config_final (object)
      class(phs_wood_config_t), intent(inout) :: object
    end subroutine phs_wood_config_final
    module subroutine phs_wood_config_increase_n_par (phs_config)
      class(phs_wood_config_t), intent(inout) :: phs_config
    end subroutine phs_wood_config_increase_n_par
    module subroutine phs_wood_config_write (object, unit, include_id)
      class(phs_wood_config_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: include_id
    end subroutine phs_wood_config_write
    module subroutine phs_wood_config_write_forest (object, unit)
      class(phs_wood_config_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine phs_wood_config_write_forest
    module subroutine phs_wood_config_set_parameters (phs_config, par)
      class(phs_wood_config_t), intent(inout) :: phs_config
      type(phs_parameters_t), intent(in) :: par
    end subroutine phs_wood_config_set_parameters
    module subroutine phs_wood_config_enable_equivalences (phs_config)
      class(phs_wood_config_t), intent(inout) :: phs_config
    end subroutine phs_wood_config_enable_equivalences
    module subroutine phs_wood_config_set_mapping_defaults &
         (phs_config, mapping_defaults)
      class(phs_wood_config_t), intent(inout) :: phs_config
      type(mapping_defaults_t), intent(in) :: mapping_defaults
    end subroutine phs_wood_config_set_mapping_defaults
    module subroutine phs_wood_config_set_input (phs_config, unit)
      class(phs_wood_config_t), intent(inout) :: phs_config
      integer, intent(in) :: unit
    end subroutine phs_wood_config_set_input
    module subroutine phs_wood_config_generate_phase_space (phs_config)
      class(phs_wood_config_t), intent(inout) :: phs_config
    end subroutine phs_wood_config_generate_phase_space
    module subroutine phs_wood_config_write_phase_space (phs_config, &
         filename_vis, unit)
      class(phs_wood_config_t), intent(in) :: phs_config
      integer, intent(in), optional :: unit
      type(string_t), intent(in), optional :: filename_vis
    end subroutine phs_wood_config_write_phase_space
    module subroutine phs_wood_config_clear_phase_space (phs_config)
      class(phs_wood_config_t), intent(inout) :: phs_config
    end subroutine phs_wood_config_clear_phase_space
    module subroutine phs_wood_config_extract_resonance_history_set &
         (phs_config, res_set, include_trivial)
      class(phs_wood_config_t), intent(in) :: phs_config
      type(resonance_history_set_t), intent(out) :: res_set
      logical, intent(in), optional :: include_trivial
    end subroutine phs_wood_config_extract_resonance_history_set
    module subroutine phs_wood_config_configure (phs_config, sqrts, &
         sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
         ignore_mismatch, nlo_type, subdir)
      class(phs_wood_config_t), intent(inout) :: phs_config
      real(default), intent(in) :: sqrts
      logical, intent(in), optional :: sqrts_fixed
      logical, intent(in), optional :: lab_is_cm
      logical, intent(in), optional :: azimuthal_dependence
      logical, intent(in), optional :: rebuild
      logical, intent(in), optional :: ignore_mismatch
      integer, intent(in), optional :: nlo_type
      type(string_t), intent(in), optional :: subdir
    end subroutine phs_wood_config_configure
    module subroutine phs_wood_config_compute_md5sum_forest (phs_config)
      class(phs_wood_config_t), intent(inout) :: phs_config
    end subroutine phs_wood_config_compute_md5sum_forest
    module function phs_wood_make_phs_filename &
         (phs_config, subdir) result (filename)
      class(phs_wood_config_t), intent(in) :: phs_config
      type(string_t), intent(in), optional :: subdir
      type(string_t) :: filename
    end function phs_wood_make_phs_filename
    module subroutine phs_wood_config_reshuffle_flavors &
         (phs_config, reshuffle, flv_extra)
      class(phs_wood_config_t), intent(inout) :: phs_config
      integer, intent(in), dimension(:), allocatable :: reshuffle
      type(flavor_t), intent(in) :: flv_extra
    end subroutine phs_wood_config_reshuffle_flavors
    module subroutine phs_wood_config_set_momentum_links (phs_config, reshuffle)
      class(phs_wood_config_t), intent(inout) :: phs_config
      integer, intent(in), dimension(:), allocatable :: reshuffle
    end subroutine phs_wood_config_set_momentum_links
    module subroutine phs_wood_config_record_s_mappings (phs_config)
      class(phs_wood_config_t), intent(inout) :: phs_config
    end subroutine phs_wood_config_record_s_mappings
    module subroutine phs_wood_config_record_on_shell (phs_config)
      class(phs_wood_config_t), intent(inout) :: phs_config
    end subroutine phs_wood_config_record_on_shell
    module function phs_wood_config_get_md5sum (phs_config) result (md5sum)
      class(phs_wood_config_t), intent(in) :: phs_config
      character(32) :: md5sum
    end function phs_wood_config_get_md5sum
    module subroutine phs_wood_read_phs_file &
         (phs_config, exist, found, match, subdir)
      class(phs_wood_config_t), intent(inout) :: phs_config
      logical, intent(out) :: exist
      logical, intent(out) :: found
      logical, intent(out), optional :: match
      type(string_t), intent(in), optional :: subdir
    end subroutine phs_wood_read_phs_file
    module subroutine phs_wood_config_startup_message (phs_config, unit)
      class(phs_wood_config_t), intent(in) :: phs_config
      integer, intent(in), optional :: unit
    end subroutine phs_wood_config_startup_message
    module subroutine phs_wood_write (object, unit, verbose)
      class(phs_wood_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine phs_wood_write
    module subroutine phs_wood_write_forest (object, unit)
      class(phs_wood_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine phs_wood_write_forest
    module subroutine phs_wood_final (object)
      class(phs_wood_t), intent(inout) :: object
    end subroutine phs_wood_final
    module subroutine phs_wood_init (phs, phs_config)
      class(phs_wood_t), intent(out) :: phs
      class(phs_config_t), intent(in), target :: phs_config
    end subroutine phs_wood_init
    module subroutine phs_wood_evaluate_selected_channel (phs, c_in, r_in)
      class(phs_wood_t), intent(inout) :: phs
      real(default), intent(in), dimension(:) :: r_in
      integer, intent(in) :: c_in
    end subroutine phs_wood_evaluate_selected_channel
    module subroutine phs_wood_evaluate_other_channels (phs, c_in)
      class(phs_wood_t), intent(inout) :: phs
      integer, intent(in) :: c_in
    end subroutine phs_wood_evaluate_other_channels
    module subroutine phs_wood_inverse (phs)
      class(phs_wood_t), intent(inout) :: phs
    end subroutine phs_wood_inverse
  end interface

contains

  subroutine phs_wood_config_allocate_instance (phs)
    class(phs_t), intent(inout), pointer :: phs
    allocate (phs_wood_t :: phs)
  end subroutine phs_wood_config_allocate_instance


end module phs_wood
