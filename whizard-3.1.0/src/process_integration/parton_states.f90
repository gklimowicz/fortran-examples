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
module parton_states

  use kinds, only: default
  use variables
  use expr_base
  use model_data
  use flavors
  use quantum_numbers
  use state_matrices
  use interactions
  use evaluators

  use beams
  use sf_base
  use prc_core
  use subevt_expr

  implicit none
  private

  public :: isolated_state_t
  public :: connected_state_t
  public :: refill_evaluator

  type, abstract :: parton_state_t
     logical :: has_trace = .false.
     logical :: has_matrix = .false.
     logical :: has_flows = .false.
     type(evaluator_t) :: trace
     type(evaluator_t) :: matrix
     type(evaluator_t) :: flows
   contains
     procedure :: write => parton_state_write
     procedure :: final => parton_state_final
     procedure :: receive_kinematics => parton_state_receive_kinematics
     procedure :: send_kinematics => parton_state_send_kinematics
     procedure :: evaluate_trace => parton_state_evaluate_trace
     procedure :: evaluate_matrix => parton_state_evaluate_matrix
     procedure :: evaluate_event_data => parton_state_evaluate_event_data
     procedure :: normalize_matrix_by_trace => &
          parton_state_normalize_matrix_by_trace
     procedure :: get_trace_int_ptr => parton_state_get_trace_int_ptr
     procedure :: get_matrix_int_ptr => parton_state_get_matrix_int_ptr
     procedure :: get_flows_int_ptr => parton_state_get_flows_int_ptr
     procedure :: get_n_out => parton_state_get_n_out
  end type parton_state_t

  type, extends (parton_state_t) :: isolated_state_t
     logical :: sf_chain_is_allocated = .false.
     type(sf_chain_instance_t), pointer :: sf_chain_eff => null ()
     logical :: int_is_allocated = .false.
     type(interaction_t), pointer :: int_eff => null ()
   contains
     procedure :: init => isolated_state_init
     procedure :: setup_square_trace => isolated_state_setup_square_trace
     procedure :: setup_identity_trace => isolated_state_setup_identity_trace
     procedure :: setup_square_matrix => isolated_state_setup_square_matrix
     procedure :: setup_square_flows => isolated_state_setup_square_flows
     procedure :: evaluate_sf_chain => isolated_state_evaluate_sf_chain
  end type isolated_state_t

  type, extends (parton_state_t) :: connected_state_t
     type(state_flv_content_t) :: state_flv
     logical :: has_flows_sf = .false.
     type(evaluator_t) :: flows_sf
     logical :: has_expr = .false.
     type(parton_expr_t) :: expr
   contains
     procedure :: setup_connected_trace => connected_state_setup_connected_trace
     procedure :: setup_connected_matrix => connected_state_setup_connected_matrix
     procedure :: setup_connected_flows => connected_state_setup_connected_flows
     procedure :: setup_state_flv => connected_state_setup_state_flv
     procedure :: get_state_flv => connected_state_get_state_flv
     procedure :: setup_subevt => connected_state_setup_subevt
     procedure :: renew_flv_content_subevt => &
          connected_state_renew_flv_content_subevt
     procedure :: setup_var_list => connected_state_setup_var_list
     procedure :: setup_cuts => connected_state_setup_cuts
     procedure :: setup_scale => connected_state_setup_scale
     procedure :: setup_fac_scale => connected_state_setup_fac_scale
     procedure :: setup_ren_scale => connected_state_setup_ren_scale
     procedure :: setup_weight => connected_state_setup_weight
     procedure :: reset_expressions => connected_state_reset_expressions
     procedure :: evaluate_expressions => connected_state_evaluate_expressions
     procedure :: get_beam_index => connected_state_get_beam_index
     procedure :: get_in_index => connected_state_get_in_index
  end type connected_state_t


  interface
    module subroutine parton_state_write (state, unit, testflag)
      class(parton_state_t), intent(in) :: state
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine parton_state_write
    module subroutine parton_state_final (state)
      class(parton_state_t), intent(inout) :: state
    end subroutine parton_state_final
    module subroutine isolated_state_init (state, sf_chain, int)
      class(isolated_state_t), intent(out) :: state
      type(sf_chain_instance_t), intent(in), target :: sf_chain
      type(interaction_t), intent(in), target :: int
    end subroutine isolated_state_init
    module subroutine isolated_state_setup_square_trace (state, core, &
         qn_mask_in, col, keep_fs_flavor)
      class(isolated_state_t), intent(inout), target :: state
      class(prc_core_t), intent(in) :: core
      type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask_in
      integer, intent(in), dimension(:), allocatable :: col
      logical, intent(in) :: keep_fs_flavor
    end subroutine isolated_state_setup_square_trace
    module subroutine isolated_state_setup_identity_trace (state, core, &
         qn_mask_in, keep_fs_flavors, keep_colors)
      class(isolated_state_t), intent(inout), target :: state
      class(prc_core_t), intent(in) :: core
      type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask_in
      logical, intent(in), optional :: keep_fs_flavors, keep_colors
    end subroutine isolated_state_setup_identity_trace
    module subroutine isolated_state_setup_square_matrix &
         (state, core, model, qn_mask_in, col)
      class(isolated_state_t), intent(inout), target :: state
      class(prc_core_t), intent(in) :: core
      class(model_data_t), intent(in), target :: model
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask_in
      integer, dimension(:), intent(in) :: col
    end subroutine isolated_state_setup_square_matrix
    module subroutine isolated_state_setup_square_flows &
         (state, core, model, qn_mask_in)
      class(isolated_state_t), intent(inout), target :: state
      class(prc_core_t), intent(in) :: core
      class(model_data_t), intent(in), target :: model
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask_in
    end subroutine isolated_state_setup_square_flows
    module subroutine connected_state_setup_connected_trace &
         (state, isolated, int, resonant, undo_helicities, &
          keep_fs_flavors, requires_extended_sf)
      class(connected_state_t), intent(inout), target :: state
      type(isolated_state_t), intent(in), target :: isolated
      type(interaction_t), intent(in), optional, target :: int
      logical, intent(in), optional :: resonant
      logical, intent(in), optional :: undo_helicities
      logical, intent(in), optional :: keep_fs_flavors
      logical, intent(in), optional :: requires_extended_sf
    end subroutine connected_state_setup_connected_trace
    module subroutine connected_state_setup_connected_matrix &
         (state, isolated, int, resonant, qn_filter_conn)
      class(connected_state_t), intent(inout), target :: state
      type(isolated_state_t), intent(in), target :: isolated
      type(interaction_t), intent(in), optional, target :: int
      logical, intent(in), optional :: resonant
      type(quantum_numbers_t), intent(in), optional :: qn_filter_conn
    end subroutine connected_state_setup_connected_matrix
    module subroutine connected_state_setup_connected_flows &
         (state, isolated, int, resonant, qn_filter_conn, mask_color)
      class(connected_state_t), intent(inout), target :: state
      type(isolated_state_t), intent(in), target :: isolated
      type(interaction_t), intent(in), optional, target :: int
      logical, intent(in), optional :: resonant, mask_color
      type(quantum_numbers_t), intent(in), optional :: qn_filter_conn
    end subroutine connected_state_setup_connected_flows
    module subroutine connected_state_setup_state_flv (state, n_out_hard)
      class(connected_state_t), intent(inout), target :: state
      integer, intent(in) :: n_out_hard
    end subroutine connected_state_setup_state_flv
    module function connected_state_get_state_flv (state) result (state_flv)
      class(connected_state_t), intent(in) :: state
      type(state_flv_content_t) :: state_flv
    end function connected_state_get_state_flv
    module subroutine connected_state_setup_subevt &
         (state, sf_chain, f_beam, f_in, f_out)
      class(connected_state_t), intent(inout), target :: state
      type(sf_chain_instance_t), intent(in), target :: sf_chain
      type(flavor_t), dimension(:), intent(in) :: f_beam, f_in, f_out
    end subroutine connected_state_setup_subevt
    module subroutine connected_state_renew_flv_content_subevt &
         (state, sf_chain, f_beam, f_in, f_out)
      class(connected_state_t), intent(inout), target :: state
      type(sf_chain_instance_t), intent(in), target :: sf_chain
      type(flavor_t), dimension(:), intent(in) :: f_beam, f_in, f_out
    end subroutine connected_state_renew_flv_content_subevt
    module subroutine connected_state_setup_var_list &
         (state, process_var_list, beam_data)
      class(connected_state_t), intent(inout), target :: state
      type(var_list_t), intent(in), target :: process_var_list
      type(beam_data_t), intent(in) :: beam_data
    end subroutine connected_state_setup_var_list
    module subroutine connected_state_setup_cuts (state, ef_cuts)
      class(connected_state_t), intent(inout), target :: state
      class(expr_factory_t), intent(in) :: ef_cuts
    end subroutine connected_state_setup_cuts
    module subroutine connected_state_setup_scale (state, ef_scale)
      class(connected_state_t), intent(inout), target :: state
      class(expr_factory_t), intent(in) :: ef_scale
    end subroutine connected_state_setup_scale
    module subroutine connected_state_setup_fac_scale (state, ef_fac_scale)
      class(connected_state_t), intent(inout), target :: state
      class(expr_factory_t), intent(in) :: ef_fac_scale
    end subroutine connected_state_setup_fac_scale
    module subroutine connected_state_setup_ren_scale (state, ef_ren_scale)
      class(connected_state_t), intent(inout), target :: state
      class(expr_factory_t), intent(in) :: ef_ren_scale
    end subroutine connected_state_setup_ren_scale
    module subroutine connected_state_setup_weight (state, ef_weight)
      class(connected_state_t), intent(inout), target :: state
      class(expr_factory_t), intent(in) :: ef_weight
    end subroutine connected_state_setup_weight
    module subroutine connected_state_reset_expressions (state)
      class(connected_state_t), intent(inout) :: state
    end subroutine connected_state_reset_expressions
    module subroutine parton_state_receive_kinematics (state)
      class(parton_state_t), intent(inout), target :: state
    end subroutine parton_state_receive_kinematics
    module subroutine parton_state_send_kinematics (state)
      class(parton_state_t), intent(inout), target :: state
    end subroutine parton_state_send_kinematics
    module subroutine connected_state_evaluate_expressions (state, passed, &
         scale, fac_scale, ren_scale, weight, scale_forced, force_evaluation)
      class(connected_state_t), intent(inout) :: state
      logical, intent(out) :: passed
      real(default), intent(out) :: scale, weight
      real(default), intent(out), allocatable :: fac_scale, ren_scale
      real(default), intent(in), allocatable, optional :: scale_forced
      logical, intent(in), optional :: force_evaluation
    end subroutine connected_state_evaluate_expressions
    module subroutine isolated_state_evaluate_sf_chain (state, fac_scale)
      class(isolated_state_t), intent(inout) :: state
      real(default), intent(in) :: fac_scale
    end subroutine isolated_state_evaluate_sf_chain
    module subroutine parton_state_evaluate_trace (state)
      class(parton_state_t), intent(inout) :: state
    end subroutine parton_state_evaluate_trace
    module subroutine parton_state_evaluate_matrix (state)
      class(parton_state_t), intent(inout) :: state
    end subroutine parton_state_evaluate_matrix
    module subroutine parton_state_evaluate_event_data (state, only_momenta)
      class(parton_state_t), intent(inout) :: state
      logical, intent(in), optional :: only_momenta
    end subroutine parton_state_evaluate_event_data
    module subroutine parton_state_normalize_matrix_by_trace (state)
      class(parton_state_t), intent(inout) :: state
    end subroutine parton_state_normalize_matrix_by_trace
    module function parton_state_get_trace_int_ptr (state) result (ptr)
      class(parton_state_t), intent(in), target :: state
      type(interaction_t), pointer :: ptr
    end function parton_state_get_trace_int_ptr
    module function parton_state_get_matrix_int_ptr (state) result (ptr)
      class(parton_state_t), intent(in), target :: state
      type(interaction_t), pointer :: ptr
    end function parton_state_get_matrix_int_ptr
    module function parton_state_get_flows_int_ptr (state) result (ptr)
      class(parton_state_t), intent(in), target :: state
      type(interaction_t), pointer :: ptr
    end function parton_state_get_flows_int_ptr
    module subroutine connected_state_get_beam_index (state, i_beam)
      class(connected_state_t), intent(in) :: state
      integer, dimension(:), intent(out) :: i_beam
    end subroutine connected_state_get_beam_index
    module subroutine connected_state_get_in_index (state, i_in)
      class(connected_state_t), intent(in) :: state
      integer, dimension(:), intent(out) :: i_in
    end subroutine connected_state_get_in_index
    module subroutine refill_evaluator (sqme, qn, flv_index, evaluator)
      complex(default), intent(in), dimension(:) :: sqme
      type(quantum_numbers_t), intent(in), dimension(:,:) :: qn
      integer, intent(in), dimension(:), optional :: flv_index
      type(evaluator_t), intent(inout) :: evaluator
    end subroutine refill_evaluator
    module function parton_state_get_n_out (state) result (n)
      class(parton_state_t), intent(in), target :: state
      integer :: n
    end function parton_state_get_n_out
  end interface

end module parton_states
