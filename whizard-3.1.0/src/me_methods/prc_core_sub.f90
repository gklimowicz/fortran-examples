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

submodule (prc_core) prc_core_s

  implicit none

contains

  module subroutine prc_core_init (object, def, lib, id, i_component)
    class(prc_core_t), intent(inout) :: object
    class(prc_core_def_t), intent(in), target :: def
    type(process_library_t), intent(in), target :: lib
    type(string_t), intent(in) :: id
    integer, intent(in) :: i_component
    object%def => def
    call lib%connect_process (id, i_component, object%data, object%driver)
    object%data_known = .true.
  end subroutine prc_core_init

  module function prc_core_has_matrix_element (object) result (flag)
    class(prc_core_t), intent(in) :: object
    logical :: flag
    flag = object%data%n_flv /= 0
  end function prc_core_has_matrix_element

  module function prc_core_needs_external_code () result (flag)
    logical :: flag
    flag = .false.
  end function prc_core_needs_external_code

  module subroutine prc_core_prepare_external_code &
       (core, flv_states, var_list, os_data, libname, model, i_core, is_nlo)
    class(prc_core_t), intent(inout) :: core
    integer, intent(in), dimension(:,:), allocatable :: flv_states
    type(var_list_t), intent(in) :: var_list
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in) :: libname
    type(model_data_t), intent(in), target :: model
    integer, intent(in) :: i_core
    logical, intent(in) :: is_nlo
    call core%write ()
    call msg_bug ("prc_core_prepare_external_code called &
         &but not overridden")
  end subroutine prc_core_prepare_external_code

  module function prc_core_uses_blha () result (flag)
    logical :: flag
    flag = .false.
  end function prc_core_uses_blha

  module subroutine prc_core_get_constants (object, data, i_term)
    class(prc_core_t), intent(in) :: object
    type(process_constants_t), intent(out) :: data
    integer, intent(in) :: i_term
    data = object%data
  end subroutine prc_core_get_constants

  module function prc_core_get_alpha_s (object, core_state) result (alpha_qcd)
    class(prc_core_t), intent(in) :: object
    class(prc_core_state_t), intent(in), allocatable :: core_state
    real(default) :: alpha_qcd
    alpha_qcd = -1
  end function prc_core_get_alpha_s

  module function prc_core_get_alpha_qed &
       (object, core_state) result (alpha_qed)
    class(prc_core_t), intent(in) :: object
    class(prc_core_state_t), intent(in), allocatable :: core_state
    real(default) :: alpha_qed
    alpha_qed = -1
  end function prc_core_get_alpha_qed

  module subroutine prc_core_set_equivalent_flv_hel_indices (object)
    class(prc_core_t), intent(inout) :: object
    integer :: i, n_flv, n_hel
    n_flv = object%data%n_flv
    n_hel = object%data%n_hel
    if (.not. allocated (object%data%eqv_flv_index)) &
         allocate (object%data%eqv_flv_index(n_flv))
    if (.not. allocated (object%data%eqv_hel_index)) &
         allocate (object%data%eqv_hel_index(n_hel))
    if (size (object%data%eqv_flv_index) /= n_flv) &
         call msg_bug ("BLHA Core: Size mismatch between eqv_flv_index and number of flavors.")
    if (size (object%data%eqv_hel_index) /= n_hel) &
         call msg_bug ("BLHA Core: Size mismatch between eqv_hel_index and number of helicities.")
    object%data%eqv_flv_index = [(i, i = 1, n_flv)]
    object%data%eqv_hel_index = [(i, i = 1, n_hel)]
  end subroutine prc_core_set_equivalent_flv_hel_indices

  module function prc_core_get_equivalent_flv_index &
       (object) result (eqv_flv_index)
    class(prc_core_t), intent(in) :: object
    integer, dimension(:), allocatable :: eqv_flv_index
    eqv_flv_index = object%data%eqv_flv_index
  end function prc_core_get_equivalent_flv_index

  module function prc_core_get_equivalent_hel_index &
       (object) result (eqv_hel_index)
    class(prc_core_t), intent(in) :: object
    integer, dimension(:), allocatable :: eqv_hel_index
    eqv_hel_index = object%data%eqv_hel_index
  end function prc_core_get_equivalent_hel_index

  module subroutine prc_core_ignore_workspace (object, core_state)
    class(prc_core_t), intent(in) :: object
    class(prc_core_state_t), intent(inout), allocatable :: core_state
  end subroutine prc_core_ignore_workspace

  module subroutine helicity_selection_write (object, unit)
    class(helicity_selection_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (object%active) then
       write (u, "(3x,A)")  "Helicity selection data:"
       write (u, "(5x,A,ES17.10)") &
            "threshold =", object%threshold
       write (u, "(5x,A,I0)") &
            "cutoff    = ", object%cutoff
    end if
  end subroutine helicity_selection_write


end submodule prc_core_s

