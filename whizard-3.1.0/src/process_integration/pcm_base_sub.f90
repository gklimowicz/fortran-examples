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

submodule (pcm_base) pcm_base_s

  use io_units
  use diagnostics
  use format_utils, only: write_integer_array
  use format_utils, only: write_separator
  use physics_defs, only: BORN, NLO_REAL

  implicit none

contains

  module function core_entry_get_core_ptr (core_entry) result (core)
    class(core_entry_t), intent(in), target :: core_entry
    class(prc_core_t), pointer :: core
    if (allocated (core_entry%core)) then
       core => core_entry%core
    else
       core => null ()
    end if
  end function core_entry_get_core_ptr

  module subroutine core_entry_configure (core_entry, lib, id)
    class(core_entry_t), intent(inout) :: core_entry
    type(process_library_t), intent(in), target :: lib
    type(string_t), intent(in) :: id
    call core_entry%core%init &
         (core_entry%core_def, lib, id, core_entry%i_component)
  end subroutine core_entry_configure

  module subroutine pcm_set_blha_defaults (pcm, polarized_beams, var_list)
    class(pcm_t), intent(inout) :: pcm
    type(var_list_t), intent(in) :: var_list
    logical, intent(in) :: polarized_beams
    logical :: muon_yukawa_off
    real(default) :: top_yukawa
    type(string_t) :: ew_scheme
    muon_yukawa_off = &
         var_list%get_lval (var_str ("?openloops_switch_off_muon_yukawa"))
    top_yukawa = &
         var_list%get_rval (var_str ("blha_top_yukawa"))
    ew_scheme = &
         var_list%get_sval (var_str ("$blha_ew_scheme"))
    if (ew_scheme == "")  ew_scheme = "Gmu"
    call pcm%blha_defaults%init &
         (polarized_beams, muon_yukawa_off, top_yukawa, ew_scheme)
  end subroutine pcm_set_blha_defaults

  module subroutine pcm_allocate_components (pcm, comp, meta)
    class(pcm_t), intent(inout) :: pcm
    type(process_component_t), dimension(:), allocatable, intent(out) :: comp
    type(process_metadata_t), intent(in) :: meta
    pcm%n_components = meta%n_components
    allocate (comp (pcm%n_components))
    allocate (pcm%component_selected (pcm%n_components), source = .false.)
    allocate (pcm%component_active (pcm%n_components), source = .true.)
  end subroutine pcm_allocate_components

  module function pcm_get_i_core (pcm, i_component) result (i_core)
    class(pcm_t), intent(in) :: pcm
    integer, intent(in) :: i_component
    integer :: i_core
    if (allocated (pcm%i_core)) then
       i_core = pcm%i_core(i_component)
    else
       i_core = 0
    end if
  end function pcm_get_i_core

  module subroutine pcm_record_inactive_components (pcm, component, meta)
    class(pcm_t), intent(inout) :: pcm
    type(process_component_t), dimension(:), intent(in) :: component
    type(process_metadata_t), intent(inout) :: meta
    integer :: i
    pcm%component_active = component%active
    do i = 1, pcm%n_components
       if (.not. component(i)%active)  call meta%deactivate_component (i)
    end do
  end subroutine pcm_record_inactive_components

  module function pcm_work_is_valid (pcm_work) result (valid)
    logical :: valid
    class(pcm_workspace_t), intent(in) :: pcm_work
    valid = .not. pcm_work%bad_point
  end function pcm_work_is_valid

  pure module subroutine pcm_work_set_bad_point (pcm_work, bad_point)
    class(pcm_workspace_t), intent(inout) :: pcm_work
    logical, intent(in) :: bad_point
    pcm_work%bad_point = pcm_work%bad_point .or. bad_point
  end subroutine pcm_work_set_bad_point


end submodule pcm_base_s

