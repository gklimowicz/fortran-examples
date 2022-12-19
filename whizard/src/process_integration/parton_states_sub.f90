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

submodule (parton_states) parton_states_s

  use debug_master, only: debug_on
  use io_units
  use format_utils, only: write_separator
  use diagnostics
  use lorentz
  use subevents
  use helicities
  use colors
  use polarizations
  use process_constants

  implicit none

contains

  module subroutine parton_state_write (state, unit, testflag)
    class(parton_state_t), intent(in) :: state
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    select type (state)
    class is (isolated_state_t)
       if (state%sf_chain_is_allocated) then
          call write_separator (u)
          call state%sf_chain_eff%write (u)
       end if
       if (state%int_is_allocated) then
          call write_separator (u)
          write (u, "(1x,A)") &
               "Effective interaction:"
          call write_separator (u)
          call state%int_eff%basic_write (u, testflag = testflag)
       end if
    class is (connected_state_t)
       if (state%has_flows_sf) then
          call write_separator (u)
          write (u, "(1x,A)") &
               "Evaluator (extension of the beam evaluator &
               &with color contractions):"
          call write_separator (u)
          call state%flows_sf%write (u, testflag = testflag)
       end if
    end select
    if (state%has_trace) then
       call write_separator (u)
       write (u, "(1x,A)") &
            "Evaluator (trace of the squared transition matrix):"
       call write_separator (u)
       call state%trace%write (u, testflag = testflag)
    end if
    if (state%has_matrix) then
       call write_separator (u)
       write (u, "(1x,A)") &
            "Evaluator (squared transition matrix):"
       call write_separator (u)
       call state%matrix%write (u, testflag = testflag)
    end if
    if (state%has_flows) then
       call write_separator (u)
       write (u, "(1x,A)") &
            "Evaluator (squared color-flow matrix):"
       call write_separator (u)
       call state%flows%write (u, testflag = testflag)
    end if
    select type (state)
    class is (connected_state_t)
       if (state%has_expr) then
          call write_separator (u)
          call state%expr%write (u)
       end if
    end select
  end subroutine parton_state_write

  module subroutine parton_state_final (state)
    class(parton_state_t), intent(inout) :: state
    if (state%has_flows) then
       call state%flows%final ()
       state%has_flows = .false.
    end if
    if (state%has_matrix) then
       call state%matrix%final ()
       state%has_matrix = .false.
    end if
    if (state%has_trace) then
       call state%trace%final ()
       state%has_trace = .false.
    end if
    select type (state)
    class is (connected_state_t)
       if (state%has_flows_sf) then
          call state%flows_sf%final ()
          state%has_flows_sf = .false.
       end if
       call state%expr%final ()
    class is (isolated_state_t)
       if (state%int_is_allocated) then
          call state%int_eff%final ()
          deallocate (state%int_eff)
          state%int_is_allocated = .false.
       end if
       if (state%sf_chain_is_allocated) then
          call state%sf_chain_eff%final ()
       end if
    end select
  end subroutine parton_state_final

  module subroutine isolated_state_init (state, sf_chain, int)
    class(isolated_state_t), intent(out) :: state
    type(sf_chain_instance_t), intent(in), target :: sf_chain
    type(interaction_t), intent(in), target :: int
    state%sf_chain_eff => sf_chain
    state%int_eff => int
  end subroutine isolated_state_init

  module subroutine isolated_state_setup_square_trace (state, core, &
       qn_mask_in, col, keep_fs_flavor)
    class(isolated_state_t), intent(inout), target :: state
    class(prc_core_t), intent(in) :: core
    type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask_in
    !!! Actually need allocatable attribute here for once because col might
    !!! enter the subroutine non-allocated.
    integer, intent(in), dimension(:), allocatable :: col
    logical, intent(in) :: keep_fs_flavor
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask
    associate (data => core%data)
      allocate (qn_mask (data%n_in + data%n_out))
      qn_mask( : data%n_in) = &
              quantum_numbers_mask (.false., .true., .false.) &
              .or. qn_mask_in
      qn_mask(data%n_in + 1 : ) = &
              quantum_numbers_mask (.not. keep_fs_flavor, .true., .true.)
      if (core%use_color_factors) then
         call state%trace%init_square (state%int_eff, qn_mask, &
              col_flow_index = data%cf_index, &
              col_factor = data%color_factors, &
              col_index_hi = col, &
              nc = core%nc)
      else
         call state%trace%init_square (state%int_eff, qn_mask, nc = core%nc)
      end if
    end associate
    state%has_trace = .true.
  end subroutine isolated_state_setup_square_trace

  module subroutine isolated_state_setup_identity_trace (state, core, &
       qn_mask_in, keep_fs_flavors, keep_colors)
    class(isolated_state_t), intent(inout), target :: state
    class(prc_core_t), intent(in) :: core
    type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask_in
    logical, intent(in), optional :: keep_fs_flavors, keep_colors
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask
    logical :: fs_flv_flag, col_flag
    fs_flv_flag = .true.; col_flag = .true.
    if (present(keep_fs_flavors)) fs_flv_flag = .not. keep_fs_flavors
    if (present(keep_colors)) col_flag = .not. keep_colors
    associate (data => core%data)
      allocate (qn_mask (data%n_in + data%n_out))
      qn_mask( : data%n_in) = &
           quantum_numbers_mask (.false., col_flag, .false.) .or. qn_mask_in
      qn_mask(data%n_in + 1 : ) = &
           quantum_numbers_mask (fs_flv_flag, col_flag, .true.)
    end associate
    call state%int_eff%set_mask (qn_mask)
    call state%trace%init_identity (state%int_eff)
    state%has_trace = .true.
  end subroutine isolated_state_setup_identity_trace

  module subroutine isolated_state_setup_square_matrix &
       (state, core, model, qn_mask_in, col)
    class(isolated_state_t), intent(inout), target :: state
    class(prc_core_t), intent(in) :: core
    class(model_data_t), intent(in), target :: model
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask_in
    integer, dimension(:), intent(in) :: col
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask
    type(flavor_t), dimension(:), allocatable :: flv
    integer :: i
    logical :: helmask, helmask_hd
    associate (data => core%data)
      allocate (qn_mask (data%n_in + data%n_out))
      allocate (flv (data%n_flv))
      do i = 1, data%n_in + data%n_out
         call flv%init (data%flv_state(i,:), model)
         if ((data%n_in == 1 .or. i > data%n_in) &
              .and. any (.not. flv%is_stable ())) then
            helmask = all (flv%decays_isotropically ())
            helmask_hd = all (flv%decays_diagonal ())
            qn_mask(i) = quantum_numbers_mask (.false., .true., helmask, &
                 mask_hd = helmask_hd)
         else if (i > data%n_in) then
            helmask = all (.not. flv%is_polarized ())
            qn_mask(i) = quantum_numbers_mask (.false., .true., helmask)
         else
            qn_mask(i) = quantum_numbers_mask (.false., .true., .false.) &
              .or. qn_mask_in(i)
         end if
      end do
      if (core%use_color_factors) then
         call state%matrix%init_square (state%int_eff, qn_mask, &
              col_flow_index = data%cf_index, &
              col_factor = data%color_factors, &
              col_index_hi = col, &
              nc = core%nc)
      else
         call state%matrix%init_square (state%int_eff, &
              qn_mask, &
              nc = core%nc)
      end if
    end associate
    state%has_matrix = .true.
  end subroutine isolated_state_setup_square_matrix

  module subroutine isolated_state_setup_square_flows &
       (state, core, model, qn_mask_in)
    class(isolated_state_t), intent(inout), target :: state
    class(prc_core_t), intent(in) :: core
    class(model_data_t), intent(in), target :: model
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask_in
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask
    type(flavor_t), dimension(:), allocatable :: flv
    integer :: i
    logical :: helmask, helmask_hd
    associate (data => core%data)
      allocate (qn_mask (data%n_in + data%n_out))
      allocate (flv (data%n_flv))
      do i = 1, data%n_in + data%n_out
         call flv%init (data%flv_state(i,:), model)
         if ((data%n_in == 1 .or. i > data%n_in) &
              .and. any (.not. flv%is_stable ())) then
            helmask = all (flv%decays_isotropically ())
            helmask_hd = all (flv%decays_diagonal ())
            qn_mask(i) = quantum_numbers_mask (.false., .false., helmask, &
                 mask_hd = helmask_hd)
         else if (i > data%n_in) then
            helmask = all (.not. flv%is_polarized ())
            qn_mask(i) = quantum_numbers_mask (.false., .false., helmask)
         else
            qn_mask(i) = quantum_numbers_mask (.false., .false., .false.) &
              .or. qn_mask_in(i)
         end if
      end do
      call state%flows%init_square (state%int_eff, qn_mask, &
           expand_color_flows = .true.)
    end associate
    state%has_flows = .true.
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
    type(quantum_numbers_mask_t) :: mask
    type(interaction_t), pointer :: src_int, beam_int
    logical :: reduce, fs_flv_flag
    if (debug_on) call msg_debug (D_PROCESS_INTEGRATION, &
         "connected_state_setup_connected_trace")
    reduce = .false.; fs_flv_flag = .true.
    if (present (undo_helicities)) reduce = undo_helicities
    if (present (keep_fs_flavors)) fs_flv_flag = .not. keep_fs_flavors
    mask = quantum_numbers_mask (fs_flv_flag, .true., .true.)
    if (present (int)) then
       src_int => int
    else
       src_int => isolated%sf_chain_eff%get_out_int_ptr ()
    end if

    if (debug2_active (D_PROCESS_INTEGRATION)) then
       call src_int%basic_write ()
    end if

    call state%trace%init_product (src_int, isolated%trace, &
         qn_mask_conn = mask, &
         qn_mask_rest = mask, &
         connections_are_resonant = resonant, &
         ignore_sub_for_qn = requires_extended_sf)

    if (reduce) then
       beam_int => isolated%sf_chain_eff%get_beam_int_ptr ()
       call undo_qn_hel (beam_int, mask, beam_int%get_n_tot ())
       call undo_qn_hel (src_int, mask, src_int%get_n_tot ())
       call beam_int%set_matrix_element (cmplx (1, 0, default))
       call src_int%set_matrix_element (cmplx (1, 0, default))
    end if

    state%has_trace = .true.
  contains
    subroutine undo_qn_hel (int_in, mask, n_tot)
      type(interaction_t), intent(inout) :: int_in
      type(quantum_numbers_mask_t), intent(in) :: mask
      integer, intent(in) :: n_tot
      type(quantum_numbers_mask_t), dimension(n_tot) :: mask_in
      mask_in = mask
      call int_in%set_mask (mask_in)
    end subroutine undo_qn_hel
  end subroutine connected_state_setup_connected_trace

  module subroutine connected_state_setup_connected_matrix &
       (state, isolated, int, resonant, qn_filter_conn)
    class(connected_state_t), intent(inout), target :: state
    type(isolated_state_t), intent(in), target :: isolated
    type(interaction_t), intent(in), optional, target :: int
    logical, intent(in), optional :: resonant
    type(quantum_numbers_t), intent(in), optional :: qn_filter_conn
    type(quantum_numbers_mask_t) :: mask
    type(interaction_t), pointer :: src_int
    mask = quantum_numbers_mask (.false., .true., .true.)
    if (present (int)) then
       src_int => int
    else
       src_int => isolated%sf_chain_eff%get_out_int_ptr ()
    end if
    call state%matrix%init_product &
         (src_int, isolated%matrix, mask, &
          qn_filter_conn = qn_filter_conn, &
          connections_are_resonant = resonant)
    state%has_matrix = .true.
  end subroutine connected_state_setup_connected_matrix

  module subroutine connected_state_setup_connected_flows &
       (state, isolated, int, resonant, qn_filter_conn, mask_color)
    class(connected_state_t), intent(inout), target :: state
    type(isolated_state_t), intent(in), target :: isolated
    type(interaction_t), intent(in), optional, target :: int
    logical, intent(in), optional :: resonant, mask_color
    type(quantum_numbers_t), intent(in), optional :: qn_filter_conn
    type(quantum_numbers_mask_t) :: mask
    type(quantum_numbers_mask_t), dimension(:), allocatable :: mask_sf
    type(interaction_t), pointer :: src_int
    logical :: mask_c
    mask_c = .false.
    if (present (mask_color))  mask_c = mask_color
    mask = quantum_numbers_mask (.false., .false., .true.)
    if (present (int)) then
       src_int => int
    else
       src_int => isolated%sf_chain_eff%get_out_int_ptr ()
       call state%flows_sf%init_color_contractions (src_int)
       state%has_flows_sf = .true.
       src_int => state%flows_sf%interaction_t
       if (mask_c) then
          allocate (mask_sf (src_int%get_n_tot ()))
          mask_sf = quantum_numbers_mask (.false., .true., .false.)
          call src_int%reduce_state_matrix (mask_sf, keep_order = .true.)
       end if
    end if
    call state%flows%init_product (src_int, isolated%flows, mask, &
         qn_filter_conn = qn_filter_conn, &
         connections_are_resonant = resonant)
    state%has_flows = .true.
  end subroutine connected_state_setup_connected_flows

  module subroutine connected_state_setup_state_flv (state, n_out_hard)
    class(connected_state_t), intent(inout), target :: state
    integer, intent(in) :: n_out_hard
    call state%matrix%get_flv_content (state%state_flv, n_out_hard)
  end subroutine connected_state_setup_state_flv

  module function connected_state_get_state_flv (state) result (state_flv)
    class(connected_state_t), intent(in) :: state
    type(state_flv_content_t) :: state_flv
    state_flv = state%state_flv
  end function connected_state_get_state_flv

  module subroutine connected_state_setup_subevt &
       (state, sf_chain, f_beam, f_in, f_out)
    class(connected_state_t), intent(inout), target :: state
    type(sf_chain_instance_t), intent(in), target :: sf_chain
    type(flavor_t), dimension(:), intent(in) :: f_beam, f_in, f_out
    integer :: n_beam, n_in, n_out, n_vir, n_tot, i, j
    integer, dimension(:), allocatable :: i_beam, i_in, i_out
    integer :: sf_out_i
    type(interaction_t), pointer :: sf_int
    sf_int => sf_chain%get_out_int_ptr ()
    n_beam = size (f_beam)
    n_in = size (f_in)
    n_out = size (f_out)
    n_vir = state%trace%get_n_vir ()
    n_tot = state%trace%get_n_tot ()
    allocate (i_beam (n_beam), i_in (n_in), i_out (n_out))
    i_beam = [(i, i = 1, n_beam)]
    do j = 1, n_in
       sf_out_i = sf_chain%get_out_i (j)
       i_in(j) = interaction_find_link &
            (state%trace%interaction_t, sf_int, sf_out_i)
    end do
    i_out = [(i, i = n_vir + 1, n_tot)]
    call state%expr%setup_subevt (state%trace%interaction_t, &
         i_beam, i_in, i_out, f_beam, f_in, f_out)
    state%has_expr = .true.
  end subroutine connected_state_setup_subevt

  module subroutine connected_state_renew_flv_content_subevt &
       (state, sf_chain, f_beam, f_in, f_out)
    class(connected_state_t), intent(inout), target :: state
    type(sf_chain_instance_t), intent(in), target :: sf_chain
    type(flavor_t), dimension(:), intent(in) :: f_beam, f_in, f_out
    integer :: n_beam, n_in, n_out, n_vir, n_tot, i, j
    integer, dimension(:), allocatable :: i_beam, i_in, i_out
    integer :: sf_out_i
    type(interaction_t), pointer :: sf_int
    sf_int => sf_chain%get_out_int_ptr ()
    n_beam = size (f_beam)
    n_in = size (f_in)
    n_out = size (f_out)
    n_vir = state%trace%get_n_vir ()
    n_tot = state%trace%get_n_tot ()
    allocate (i_beam (n_beam), i_in (n_in), i_out (n_out))
    i_beam = [(i, i = 1, n_beam)]
    do j = 1, n_in
       sf_out_i = sf_chain%get_out_i (j)
       i_in(j) = interaction_find_link &
            (state%trace%interaction_t, sf_int, sf_out_i)
    end do
    i_out = [(i, i = n_vir + 1, n_tot)]
    call state%expr%renew_flv_content_subevt (state%trace%interaction_t, &
         i_beam, i_in, i_out, f_beam, f_in, f_out)
    state%has_expr = .true.
  end subroutine connected_state_renew_flv_content_subevt

  module subroutine connected_state_setup_var_list &
       (state, process_var_list, beam_data)
    class(connected_state_t), intent(inout), target :: state
    type(var_list_t), intent(in), target :: process_var_list
    type(beam_data_t), intent(in) :: beam_data
    call state%expr%setup_vars (beam_data%get_sqrts ())
    call state%expr%link_var_list (process_var_list)
  end subroutine connected_state_setup_var_list

  module subroutine connected_state_setup_cuts (state, ef_cuts)
    class(connected_state_t), intent(inout), target :: state
    class(expr_factory_t), intent(in) :: ef_cuts
    call state%expr%setup_selection (ef_cuts)
  end subroutine connected_state_setup_cuts

  module subroutine connected_state_setup_scale (state, ef_scale)
    class(connected_state_t), intent(inout), target :: state
    class(expr_factory_t), intent(in) :: ef_scale
    call state%expr%setup_scale (ef_scale)
  end subroutine connected_state_setup_scale

  module subroutine connected_state_setup_fac_scale (state, ef_fac_scale)
    class(connected_state_t), intent(inout), target :: state
    class(expr_factory_t), intent(in) :: ef_fac_scale
    call state%expr%setup_fac_scale (ef_fac_scale)
  end subroutine connected_state_setup_fac_scale

  module subroutine connected_state_setup_ren_scale (state, ef_ren_scale)
    class(connected_state_t), intent(inout), target :: state
    class(expr_factory_t), intent(in) :: ef_ren_scale
    call state%expr%setup_ren_scale (ef_ren_scale)
  end subroutine connected_state_setup_ren_scale

  module subroutine connected_state_setup_weight (state, ef_weight)
    class(connected_state_t), intent(inout), target :: state
    class(expr_factory_t), intent(in) :: ef_weight
    call state%expr%setup_weight (ef_weight)
  end subroutine connected_state_setup_weight

  module subroutine connected_state_reset_expressions (state)
    class(connected_state_t), intent(inout) :: state
    if (state%has_expr)  call state%expr%reset_contents ()
  end subroutine connected_state_reset_expressions

  module subroutine parton_state_receive_kinematics (state)
    class(parton_state_t), intent(inout), target :: state
    if (state%has_trace) then
       call state%trace%receive_momenta ()
       select type (state)
       class is (connected_state_t)
          if (state%has_expr) then
             call state%expr%fill_subevt (state%trace%interaction_t)
          end if
       end select
    end if
  end subroutine parton_state_receive_kinematics

  module subroutine parton_state_send_kinematics (state)
    class(parton_state_t), intent(inout), target :: state
    if (state%has_trace) then
       call state%trace%send_momenta ()
       select type (state)
       class is (connected_state_t)
          call state%expr%fill_subevt (state%trace%interaction_t)
       end select
    end if
  end subroutine parton_state_send_kinematics

  module subroutine connected_state_evaluate_expressions (state, passed, &
       scale, fac_scale, ren_scale, weight, scale_forced, force_evaluation)
    class(connected_state_t), intent(inout) :: state
    logical, intent(out) :: passed
    real(default), intent(out) :: scale, weight
    real(default), intent(out), allocatable :: fac_scale, ren_scale
    real(default), intent(in), allocatable, optional :: scale_forced
    logical, intent(in), optional :: force_evaluation
    if (state%has_expr) then
       call state%expr%evaluate (passed, scale, fac_scale, ren_scale, weight, &
            scale_forced, force_evaluation)
    end if
  end subroutine connected_state_evaluate_expressions

  module subroutine isolated_state_evaluate_sf_chain (state, fac_scale)
    class(isolated_state_t), intent(inout) :: state
    real(default), intent(in) :: fac_scale
    if (state%sf_chain_is_allocated)  &
         call state%sf_chain_eff%evaluate (fac_scale)
  end subroutine isolated_state_evaluate_sf_chain

  module subroutine parton_state_evaluate_trace (state)
    class(parton_state_t), intent(inout) :: state
    if (state%has_trace) call state%trace%evaluate ()
  end subroutine parton_state_evaluate_trace

  module subroutine parton_state_evaluate_matrix (state)
    class(parton_state_t), intent(inout) :: state
    if (state%has_matrix) call state%matrix%evaluate ()
  end subroutine parton_state_evaluate_matrix

  module subroutine parton_state_evaluate_event_data (state, only_momenta)
    class(parton_state_t), intent(inout) :: state
    logical, intent(in), optional :: only_momenta
    logical :: only_mom
    only_mom = .false.; if (present (only_momenta)) only_mom = only_momenta
    select type (state)
    type is (connected_state_t)
       if (state%has_flows_sf) then
          call state%flows_sf%receive_momenta ()
          if (.not. only_mom) call state%flows_sf%evaluate ()
       end if
    end select
    if (state%has_matrix) then
       call state%matrix%receive_momenta ()
       if (.not. only_mom) call state%matrix%evaluate ()
    end if
    if (state%has_flows) then
       call state%flows%receive_momenta ()
       if (.not. only_mom) call state%flows%evaluate ()
    end if
  end subroutine parton_state_evaluate_event_data

  module subroutine parton_state_normalize_matrix_by_trace (state)
    class(parton_state_t), intent(inout) :: state
    if (state%has_matrix) call state%matrix%normalize_by_trace ()
  end subroutine parton_state_normalize_matrix_by_trace

  module function parton_state_get_trace_int_ptr (state) result (ptr)
    class(parton_state_t), intent(in), target :: state
    type(interaction_t), pointer :: ptr
    if (state%has_trace) then
       ptr => state%trace%interaction_t
    else
       ptr => null ()
    end if
  end function parton_state_get_trace_int_ptr

  module function parton_state_get_matrix_int_ptr (state) result (ptr)
    class(parton_state_t), intent(in), target :: state
    type(interaction_t), pointer :: ptr
    if (state%has_matrix) then
       ptr => state%matrix%interaction_t
    else
       ptr => null ()
    end if
  end function parton_state_get_matrix_int_ptr

  module function parton_state_get_flows_int_ptr (state) result (ptr)
    class(parton_state_t), intent(in), target :: state
    type(interaction_t), pointer :: ptr
    if (state%has_flows) then
       ptr => state%flows%interaction_t
    else
       ptr => null ()
    end if
  end function parton_state_get_flows_int_ptr

  module subroutine connected_state_get_beam_index (state, i_beam)
    class(connected_state_t), intent(in) :: state
    integer, dimension(:), intent(out) :: i_beam
    call state%expr%get_beam_index (i_beam)
  end subroutine connected_state_get_beam_index

  module subroutine connected_state_get_in_index (state, i_in)
    class(connected_state_t), intent(in) :: state
    integer, dimension(:), intent(out) :: i_in
    call state%expr%get_in_index (i_in)
  end subroutine connected_state_get_in_index

  module subroutine refill_evaluator (sqme, qn, flv_index, evaluator)
    complex(default), intent(in), dimension(:) :: sqme
    type(quantum_numbers_t), intent(in), dimension(:,:) :: qn
    integer, intent(in), dimension(:), optional :: flv_index
    type(evaluator_t), intent(inout) :: evaluator
    integer :: i, i_flv
    do i = 1, size (sqme)
       if (present (flv_index)) then
          i_flv = flv_index(i)
       else
          i_flv = i
       end if
       call evaluator%add_to_matrix_element (qn(:,i_flv), sqme(i), &
            match_only_flavor = .true.)
    end do
  end subroutine refill_evaluator

  module function parton_state_get_n_out (state) result (n)
    class(parton_state_t), intent(in), target :: state
    integer :: n
    n = state%trace%get_n_out ()
  end function parton_state_get_n_out


end submodule parton_states_s

