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

submodule (resonance_insertion) resonance_insertion_s

  use io_units
  use format_utils, only: write_separator
  use format_defs, only: FMT_12
  use interactions, only: interaction_t
  use subevents, only: PRT_RESONANT

  implicit none

contains

  module subroutine evt_resonance_write_name (evt, unit)
    class(evt_resonance_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Event transform: resonance insertion"
  end subroutine evt_resonance_write_name

  module subroutine evt_resonance_write &
       (evt, unit, verbose, more_verbose, testflag)
    class(evt_resonance_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, more_verbose, testflag
    integer :: u, i
    u = given_output_unit (unit)
    call write_separator (u, 2)
    call evt%write_name (u)
    call write_separator (u, 2)
    write (u, "(1x,A,A,A)")  "Process library = '", char (evt%libname), "'"
    if (allocated (evt%res_history_set)) then
       do i = 1, size (evt%res_history_set)
          if (i == evt%selected_component) then
             write (u, "(1x,A,I0,A)")  "Component #", i, ": *"
          else
             write (u, "(1x,A,I0,A)")  "Component #", i, ":"
          end if
          call evt%res_history_set(i)%write (u, indent=1)
       end do
    end if
    call write_separator (u)
    if (allocated (evt%instance)) then
       write (u, "(1x,A)")  "Subprocess instances: allocated"
    else
       write (u, "(1x,A)")  "Subprocess instances: not allocated"
    end if
    if (evt%particle_set_exists) then
       if (evt%selected_history > 0) then
          write (u, "(1x,A,I0)")  "Selected: resonance history #", &
               evt%selected_history
       else
          write (u, "(1x,A)")  "Selected: no resonance history"
       end if
    else
       write (u, "(1x,A)")  "Selected: [none]"
    end if
    write (u, "(1x,A,1x," // FMT_12 // ")")  &
         "On-shell limit    =", evt%on_shell_limit
    write (u, "(1x,A,1x," // FMT_12 // ")")  &
         "On-shell turnoff  =", evt%on_shell_turnoff
    write (u, "(1x,A,1x," // FMT_12 // ")")  &
         "Background factor =", evt%background_factor
    call write_separator (u)
    if (evt%selector_active) then
       write (u, "(2x)", advance="no")
       call evt%selector%write (u, testflag=testflag)
       call write_separator (u)
    end if
    call evt%base_write (u, testflag = testflag, show_set = .false.)
    call write_separator (u)
    if (evt%particle_set_exists) then
       call evt%particle_set%write &
            (u, summary = .true., compressed = .true., testflag = testflag)
       call write_separator (u)
    end if
  end subroutine evt_resonance_write

  module subroutine evt_resonance_set_resonance_data (evt, res_history_set)
    class(evt_resonance_t), intent(inout) :: evt
    type(resonance_history_set_t), dimension(:), intent(in) :: res_history_set
    integer :: i
    evt%res_history_set = res_history_set
    allocate (evt%index_offset (size (evt%res_history_set)), source = 0)
    do i = 2, size (evt%res_history_set)
       evt%index_offset(i) = &
            evt%index_offset(i-1) + evt%res_history_set(i-1)%get_n_history ()
    end do
  end subroutine evt_resonance_set_resonance_data

  module subroutine evt_resonance_set_library (evt, libname)
    class(evt_resonance_t), intent(inout) :: evt
    type(string_t), intent(in) :: libname
    evt%libname = libname
  end subroutine evt_resonance_set_library

  module subroutine evt_resonance_set_subprocess_instances (evt, instance)
    class(evt_resonance_t), intent(inout) :: evt
    type(process_instance_ptr_t), dimension(:), intent(in) :: instance
    evt%instance = instance
  end subroutine evt_resonance_set_subprocess_instances

  module subroutine evt_resonance_set_on_shell_limit (evt, on_shell_limit)
    class(evt_resonance_t), intent(inout) :: evt
    real(default), intent(in) :: on_shell_limit
    evt%on_shell_limit = on_shell_limit
  end subroutine evt_resonance_set_on_shell_limit

  module subroutine evt_resonance_set_on_shell_turnoff (evt, on_shell_turnoff)
    class(evt_resonance_t), intent(inout) :: evt
    real(default), intent(in) :: on_shell_turnoff
    evt%on_shell_turnoff = on_shell_turnoff
  end subroutine evt_resonance_set_on_shell_turnoff

  module subroutine evt_resonance_set_background_factor (evt, background_factor)
    class(evt_resonance_t), intent(inout) :: evt
    real(default), intent(in) :: background_factor
    evt%background_factor = background_factor
  end subroutine evt_resonance_set_background_factor

  module subroutine evt_resonance_import_rng (evt, rng)
    class(evt_resonance_t), intent(inout) :: evt
    class(rng_t), allocatable, intent(inout) :: rng
    call move_alloc (from = rng, to = evt%rng)
  end subroutine evt_resonance_import_rng

  module subroutine evt_resonance_write_selector (evt, unit, testflag)
    class(evt_resonance_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (evt%selector_active) then
       call evt%selector%write (u, testflag)
    else
       write (u, "(1x,A)")  "Selector: [inactive]"
    end if
  end subroutine evt_resonance_write_selector

  module subroutine evt_resonance_init_selector (evt, weight, offset)
    class(evt_resonance_t), intent(inout) :: evt
    real(default), dimension(:), intent(in) :: weight
    integer, intent(in), optional :: offset
    if (any (weight > 0)) then
       call evt%selector%init (weight, offset = offset)
       evt%selector_active = .true.
    else
       evt%selector_active = .false.
    end if
  end subroutine evt_resonance_init_selector

  module subroutine evt_resonance_get_selector_weights (evt, weight)
    class(evt_resonance_t), intent(in) :: evt
    real(default), dimension(0:), intent(out) :: weight
    integer :: i
    do i = 0, ubound (weight,1)
       weight(i) = evt%selector%get_weight (i)
    end do
  end subroutine evt_resonance_get_selector_weights

  module subroutine evt_resonance_fill_momenta (evt)
    class(evt_resonance_t), intent(inout) :: evt
    integer :: i, n
    if (associated (evt%previous)) then
       evt%particle_set = evt%previous%particle_set
    else if (associated (evt%process_instance)) then
       ! this branch only for unit test
       call evt%process_instance%get_trace &
            (evt%particle_set, i_term=1, n_incoming=evt%process%get_n_in ())
    end if
  end subroutine evt_resonance_fill_momenta

  module subroutine evt_resonance_determine_on_shell_histories &
       (evt, index_array)
    class(evt_resonance_t), intent(in) :: evt
    integer, dimension(:), allocatable, intent(out) :: index_array
    integer :: i
    i = evt%selected_component
    call evt%res_history_set(i)%determine_on_shell_histories &
         (evt%particle_set%get_outgoing_momenta (), &
         evt%on_shell_limit, &
         index_array)
  end subroutine evt_resonance_determine_on_shell_histories

  module subroutine evt_resonance_evaluate_subprocess (evt, index_array)
    class(evt_resonance_t), intent(inout) :: evt
    integer, dimension(:), intent(in) :: index_array
    integer :: k, i
    if (allocated (evt%instance)) then
       do k = 1, size (index_array)
          i = index_array(k)
          associate (instance => evt%instance(i)%p)
            call instance%choose_mci (1)
            call instance%set_trace (evt%particle_set, 1, check_match=.false.)
            call instance%recover (channel = 1, i_term = 1, &
                 update_sqme = .true., recover_phs = .false.)
          end associate
       end do
    end if
  end subroutine evt_resonance_evaluate_subprocess

  module function evt_resonance_get_master_sqme (evt) result (sqme)
    class(evt_resonance_t), intent(in) :: evt
    real(default) :: sqme
    sqme = evt%process_instance%get_sqme ()
  end function evt_resonance_get_master_sqme

  module subroutine evt_resonance_get_subprocess_sqme (evt, sqme, index_array)
    class(evt_resonance_t), intent(in) :: evt
    real(default), dimension(:), intent(out) :: sqme
    integer, dimension(:), intent(in), optional :: index_array
    integer :: k, i
    if (present (index_array)) then
       sqme = 0
       do k = 1, size (index_array)
          call get_sqme (index_array(k))
       end do
    else
       do i = 1, size (evt%instance)
          call get_sqme (i)
       end do
    end if
  contains
    subroutine get_sqme (i)
      integer, intent(in) :: i
      associate (instance => evt%instance(i)%p)
        sqme(i) = instance%get_sqme ()
      end associate
    end subroutine get_sqme
  end subroutine evt_resonance_get_subprocess_sqme

  module subroutine evt_resonance_apply_turnoff_factor (evt, sqme, index_array)
    class(evt_resonance_t), intent(in) :: evt
    real(default), dimension(:), intent(inout) :: sqme
    integer, dimension(:), intent(in) :: index_array
    integer :: k, i_res, i_prc
    do k = 1, size (index_array)
       i_res = evt%selected_component
       i_prc = index_array(k) + evt%index_offset(i_res)
       sqme(i_prc) = sqme(i_prc) &
            * evt%res_history_set(i_res)%evaluate_gaussian &
            &   (evt%particle_set%get_outgoing_momenta (), &
            &    evt%on_shell_turnoff, index_array(k))
    end do
  end subroutine evt_resonance_apply_turnoff_factor

  module subroutine evt_resonance_compute_probabilities (evt)
    class(evt_resonance_t), intent(inout) :: evt
    integer, dimension(:), allocatable :: index_array
    real(default) :: sqme_master, sqme_sum, sqme_bg
    real(default), dimension(:), allocatable :: sqme_res
    integer :: n, ic
    if (.not. associated (evt%process_instance))  return
    n = size (evt%instance)
    call evt%select_component (0)
    FIND_ACTIVE_COMPONENT: do ic = 1, evt%process%get_n_components ()
       if (evt%process%component_is_selected (ic)) then
          call evt%select_component (ic)
          exit FIND_ACTIVE_COMPONENT
       end if
    end do FIND_ACTIVE_COMPONENT
    if (evt%selected_component > 0) then
       call evt%determine_on_shell_histories (index_array)
    else
       allocate (index_array (0))
    end if
    call evt%evaluate_subprocess &
         (index_array + evt%index_offset(evt%selected_component))
    allocate (sqme_res (n), source = 0._default)
    call evt%get_subprocess_sqme &
         (sqme_res, index_array + evt%index_offset(evt%selected_component))
    sqme_master = evt%get_master_sqme ()
    sqme_sum = sum (sqme_res)
    sqme_bg = abs (sqme_master - sqme_sum)
    if (evt%on_shell_turnoff > 0) then
       call evt%apply_turnoff_factor (sqme_res, index_array)
    end if
    if (any (sqme_res > 0)) then
       sqme_bg = sqme_bg * evt%background_factor
    end if
    call evt%init_selector ([sqme_bg, sqme_res], offset = -1)
  end subroutine evt_resonance_compute_probabilities

  module subroutine evt_resonance_select_component (evt, i_component)
    class(evt_resonance_t), intent(inout) :: evt
    integer, intent(in) :: i_component
    evt%selected_component = i_component
  end subroutine evt_resonance_select_component

  module subroutine evt_resonance_find_prt_invalid_color (evt, index, prt)
    class(evt_resonance_t), intent(in) :: evt
    integer, dimension(:), allocatable, intent(out) :: index
    type(particle_t), dimension(:), allocatable, intent(out), optional :: prt
    if (evt%particle_set_exists) then
       call evt%particle_set%find_prt_invalid_color (index, prt)
    else
       allocate (prt (0))
    end if
  end subroutine evt_resonance_find_prt_invalid_color

  module subroutine evt_resonance_prepare_new_event (evt, i_mci, i_term)
    class(evt_resonance_t), intent(inout) :: evt
    integer, intent(in) :: i_mci, i_term
    call evt%reset ()
  end subroutine evt_resonance_prepare_new_event

  module subroutine evt_resonance_generate_weighted (evt, probability)
    class(evt_resonance_t), intent(inout) :: evt
    real(default), intent(inout) :: probability
    integer :: i_term
    call evt%fill_momenta ()
    call evt%compute_probabilities ()
    if (associated (evt%process_instance)) then
       i_term = evt%process_instance%get_first_active_i_term ()
       if (evt%process_instance%term(i_term)%passed) then
          call evt%selector%generate (evt%rng, evt%selected_history)
       end if
    else
       call evt%selector%generate (evt%rng, evt%selected_history)
    end if
    probability = 1
  end subroutine evt_resonance_generate_weighted

  module subroutine evt_resonance_make_particle_set &
       (evt, factorization_mode, keep_correlations, r)
    class(evt_resonance_t), intent(inout) :: evt
    integer, intent(in) :: factorization_mode
    logical, intent(in) :: keep_correlations
    real(default), dimension(:), intent(in), optional :: r
    type(particle_set_t), target :: prt_set
    type(particle_t), dimension(:), allocatable :: prt
    integer :: n_beam, n_in, n_vir, n_res, n_out, i, i_res, i_term, i_tree
    type(interaction_t), pointer :: int_matrix, int_flows
    integer, dimension(:), allocatable :: map
    type(resonance_tree_t) :: res_tree
    if (associated (evt%previous)) then
       if (evt%previous%particle_set_exists) then
          if (evt%selected_history > 0) then
             if (allocated (evt%instance)) then
                associate (instance => evt%instance(evt%selected_history)%p)
                  call instance%evaluate_event_data (weight = 1._default)
                  i_term = 1
                  int_matrix => instance%get_matrix_int_ptr (i_term)
                  int_flows  => instance%get_flows_int_ptr (i_term)
                  call evt%factorize_interactions (int_matrix, int_flows, &
                       factorization_mode, keep_correlations, r)
                  call evt%tag_incoming ()
                end associate
             else  ! this branch only for unit test
                evt%particle_set = evt%previous%particle_set
             end if
             i_tree = evt%selected_history &
                  - evt%index_offset(evt%selected_component)
             call evt%res_history_set(evt%selected_component)%get_tree &
                  (i_tree, res_tree)
             n_beam = evt%particle_set%get_n_beam ()
             n_in = evt%particle_set%get_n_in ()
             n_vir = evt%particle_set%get_n_vir ()
             n_out = evt%particle_set%get_n_out ()
             n_res = res_tree%get_n_resonances ()
             allocate (map (n_beam + n_in + n_vir + n_out))
             map(1:n_beam+n_in+n_vir) &
                  = [(i, i = 1, n_beam+n_in+n_vir)]
             map(n_beam+n_in+n_vir+1:n_beam+n_in+n_vir+n_out) &
                  = [(i + n_res, &
                  &   i = n_beam+n_in+n_vir+1, &
                  &       n_beam+n_in+n_vir+n_out)]
             call prt_set%transfer (evt%particle_set, n_res, map)
             do i = 1, n_res
                i_res = n_beam + n_in + n_vir + i
                call prt_set%insert (i_res, &
                     PRT_RESONANT, &
                     res_tree%get_flv (i), &
                     res_tree%get_children (i, &
                     &   n_beam+n_in+n_vir, n_beam+n_in+n_vir+n_res))
             end do
             do i = n_res, 1, -1
                i_res = n_beam + n_in + n_vir + i
                call prt_set%recover_color (i_res)
             end do
             call prt_set%set_momentum &
                  (map(:), evt%particle_set%get_momenta (), on_shell = .true.)
             do i = n_res, 1, -1
                i_res = n_beam + n_in + n_vir + i
                call prt_set%recover_momentum (i_res)
             end do
             call evt%particle_set%final ()
             evt%particle_set = prt_set
             call prt_set%final ()
             evt%particle_set_exists = .true.
          else  ! retain particle set, as copied from previous evt
             evt%particle_set_exists = .true.
          end if
       else
          evt%particle_set_exists = .false.
       end if
    else
       evt%particle_set_exists = .false.
    end if
  end subroutine evt_resonance_make_particle_set


end submodule resonance_insertion_s

