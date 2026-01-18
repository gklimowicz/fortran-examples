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

submodule (decays) decays_s

  use io_units
  use format_utils, only: write_indent, write_separator
  use format_defs, only: FMT_15
  use numeric_utils
  use helicities
  use quantum_numbers

  implicit none

contains

  recursive module subroutine decay_term_config_final (object)
    class(decay_term_config_t), intent(inout) :: object
    integer :: i
    if (allocated (object%prt)) then
       do i = 1, size (object%prt)
          if (allocated (object%prt(i)%c))  call object%prt(i)%c%final ()
       end do
    end if
  end subroutine decay_term_config_final

  recursive module subroutine decay_term_config_write &
       (object, unit, indent, verbose)
    class(decay_term_config_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    logical, intent(in), optional :: verbose
    integer :: i, j, u, ind
    logical :: verb
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    verb = .true.;  if (present (verbose))  verb = verbose
    call write_indent (u, ind)
    write (u, "(1x,A)", advance="no")  "Final state:"
    do i = 1, size (object%prt)
       select type (prt_config => object%prt(i)%c)
       type is (stable_config_t)
          write (u, "(1x,A)", advance="no") &
               char (prt_config%flv(1)%get_name ())
          do j = 2, size (prt_config%flv)
             write (u, "(':',A)", advance="no") &
                  char (prt_config%flv(j)%get_name ())
          end do
       type is (unstable_config_t)
          write (u, "(1x,A)", advance="no") &
               char (prt_config%flv%get_name ())
       end select
    end do
    write (u, *)
    if (verb) then
       do i = 1, size (object%prt)
          call object%prt(i)%c%write (u, ind)
       end do
    end if
  end subroutine decay_term_config_write

  recursive module subroutine decay_term_config_compute (term)
    class(decay_term_config_t), intent(inout) :: term
    integer :: i
    do i = 1, size (term%prt)
       select type (unstable_config => term%prt(i)%c)
       type is (unstable_config_t)
          call unstable_config%compute ()
       end select
    end do
  end subroutine decay_term_config_compute

  recursive module subroutine decay_term_final (object)
    class(decay_term_t), intent(inout) :: object
    integer :: i
    if (allocated (object%particle_out)) then
       do i = 1, size (object%particle_out)
          call object%particle_out(i)%c%final ()
       end do
    end if
  end subroutine decay_term_final

  recursive module subroutine decay_term_write (object, unit, indent)
    class(decay_term_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    integer :: i, u, ind
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    call object%config%write (u, ind, verbose = .false.)
    do i = 1, size (object%particle_out)
       call object%particle_out(i)%c%write (u, ind)
    end do
  end subroutine decay_term_write

  recursive module subroutine decay_term_write_process_instances &
       (term, unit, verbose)
    class(decay_term_t), intent(in) :: term
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: i
    do i = 1, size (term%particle_out)
       select type (unstable => term%particle_out(i)%c)
       type is (unstable_t)
          call unstable%write_process_instances (unit, verbose)
       end select
    end do
  end subroutine decay_term_write_process_instances

  module subroutine decay_term_make_rng (term, process)
    class(decay_term_t), intent(inout) :: term
    type(process_t), intent(inout) :: process
    class(rng_t), allocatable :: rng
    integer :: i
    do i = 1, size (term%particle_out)
       select type (unstable => term%particle_out(i)%c)
       type is (unstable_t)
          call process%make_rng (rng)
          call unstable%import_rng (rng)
       end select
    end do
  end subroutine decay_term_make_rng

  recursive module subroutine decay_term_link_interactions (term, trace)
    class(decay_term_t), intent(inout) :: term
    type(interaction_t), intent(in), target :: trace
    integer :: i
    do i = 1, size (term%particle_out)
       select type (unstable => term%particle_out(i)%c)
       type is (unstable_t)
          call unstable%link_interactions (i, trace)
       end select
    end do
  end subroutine decay_term_link_interactions

  recursive module subroutine decay_term_select_chain (term)
    class(decay_term_t), intent(inout) :: term
    integer :: i
    do i = 1, size (term%particle_out)
       select type (unstable => term%particle_out(i)%c)
       type is (unstable_t)
          call unstable%select_chain ()
       end select
    end do
  end subroutine decay_term_select_chain

  recursive module subroutine decay_term_generate (term)
    class(decay_term_t), intent(inout) :: term
    integer :: i
    do i = 1, size (term%particle_out)
       select type (unstable => term%particle_out(i)%c)
       type is (unstable_t)
          call unstable%generate ()
       end select
    end do
  end subroutine decay_term_generate

  recursive module subroutine decay_root_config_final (object)
    class(decay_root_config_t), intent(inout) :: object
    integer :: i
    if (allocated (object%term_config)) then
       do i = 1, size (object%term_config)
          call object%term_config(i)%final ()
       end do
    end if
  end subroutine decay_root_config_final

  recursive module subroutine decay_root_config_write &
       (object, unit, indent, verbose)
    class(decay_root_config_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    logical, intent(in), optional :: verbose
    integer :: u, ind
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    call write_indent (u, ind)
    write (u, "(1x,A)")  "Final-state decay tree:"
    call object%write_header (unit, indent)
    call object%write_terms (unit, indent, verbose)
  end subroutine decay_root_config_write

  module subroutine decay_root_config_write_header (object, unit, indent)
    class(decay_root_config_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    integer :: u, ind
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    call write_indent (u, ind)
    if (associated (object%process)) then
       write (u, 3)  "process ID      =", char (object%process_id), "*"
    else
       write (u, 3)  "process ID      =", char (object%process_id)
    end if
3   format (3x,A,2(1x,A))
  end subroutine decay_root_config_write_header

  module recursive subroutine decay_root_config_write_terms &
       (object, unit, indent, verbose)
    class(decay_root_config_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    logical, intent(in), optional :: verbose
    integer :: i, u, ind
    logical :: verb
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    verb = .true.;  if (present (verbose))  verb = verbose
    if (verb .and. allocated (object%term_config)) then
       do i = 1, size (object%term_config)
          call object%term_config(i)%write (u, ind + 1)
       end do
    end if
  end subroutine decay_root_config_write_terms

  module subroutine decay_root_config_init (decay, model, process_id, n_terms)
    class(decay_root_config_t), intent(out) :: decay
    class(model_data_t), intent(in), target :: model
    type(string_t), intent(in) :: process_id
    integer, intent(in), optional :: n_terms
    decay%model => model
    decay%process_id = process_id
    if (present (n_terms)) then
       allocate (decay%term_config (n_terms))
    end if
  end subroutine decay_root_config_init

  recursive module subroutine decay_root_config_init_term &
       (decay, i, flv, stable, model, process_stack, var_list)
    class(decay_root_config_t), intent(inout) :: decay
    integer, intent(in) :: i
    type(flavor_t), dimension(:,:), intent(in) :: flv
    logical, dimension(:), intent(in) :: stable
    class(model_data_t), intent(in), target :: model
    type(process_stack_t), intent(in), optional :: process_stack
    type(var_list_t), intent(in), optional, target :: var_list
    call decay%term_config(i)%init (flv, stable, model, process_stack, var_list)
  end subroutine decay_root_config_init_term

  recursive module subroutine decay_root_config_connect &
       (decay, process, model, process_stack, process_instance, var_list)
    class(decay_root_config_t), intent(out) :: decay
    type(process_t), intent(in), target :: process
    class(model_data_t), intent(in), target :: model
    type(process_stack_t), intent(in), optional :: process_stack
    type(process_instance_t), intent(in), optional, target :: process_instance
    type(var_list_t), intent(in), optional, target :: var_list
    type(connected_state_t), pointer :: connected_state
    type(interaction_t), pointer :: int
    type(flavor_t), dimension(:,:), allocatable :: flv
    logical, dimension(:), allocatable :: stable
    real(default), dimension(:), allocatable :: m_prod, m_dec
    integer :: i
    call decay%init (model, process%get_id (), process%get_n_terms ())
    do i = 1, size (decay%term_config)
       if (present (process_instance)) then
          connected_state => process_instance%get_connected_state_ptr (i)
          int => connected_state%get_matrix_int_ptr ()
          call int%get_flv_out (flv)
       else
          call process%get_term_flv_out (i, flv)
       end if
       allocate (m_prod (size (flv(:,1)%get_mass ())))
       m_prod = flv(:,1)%get_mass ()
       call flv%set_model (model)
       allocate (m_dec (size (flv(:,1)%get_mass ())))
       m_dec = flv(:,1)%get_mass ()
       allocate (stable (size (flv, 1)))
       stable = flv(:,1)%is_stable ()
       call check_masses ()
       call decay%init_term (i, flv, stable, model, process_stack, var_list)
       deallocate (flv, stable, m_prod, m_dec)
    end do
    decay%process => process
  contains
    subroutine check_masses ()
      integer :: i
      logical :: ok
      ok = .true.
      do i = 1, size (m_prod)
         if (.not. stable(i)) then
            if (.not. nearly_equal (m_prod(i), m_dec(i))) then
               write (msg_buffer, "(A,A,A)")  "particle '", &
                    char (flv(i,1)%get_name ()), "':"
               call msg_message
               write (msg_buffer, &
                    "(2x,A,1x," // FMT_15 // ",3x,A,1x," // FMT_15 // ")") &
                    "m_prod =", m_prod(i), "m_dec =", m_dec(i)
               call msg_message
               ok = .false.
            end if
         end if
      end do
      if (.not. ok)  call msg_fatal &
           ("Particle mass mismatch between production and decay")
    end subroutine check_masses
  end subroutine decay_root_config_connect

  recursive module subroutine decay_root_config_compute (decay)
    class(decay_root_config_t), intent(inout) :: decay
    integer :: i
    do i = 1, size (decay%term_config)
       call decay%term_config(i)%compute ()
    end do
  end subroutine decay_root_config_compute

  recursive module subroutine decay_gen_final (object)
    class(decay_gen_t), intent(inout) :: object
    integer :: i
    if (allocated (object%term)) then
       do i = 1, size (object%term)
          call object%term(i)%final ()
       end do
    end if
  end subroutine decay_gen_final

  module subroutine decay_root_final (object)
    class(decay_root_t), intent(inout) :: object
    call object%base_final ()
  end subroutine decay_root_final

  module subroutine decay_root_write (object, unit)
    class(decay_root_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%config)) then
       call object%config%write (unit, verbose = .false.)
    else
       write (u, "(1x,A)")  "Final-state decay tree: [not configured]"
    end if
    if (object%selected_mci > 0) then
       write (u, "(3x,A,I0)")  "Selected MCI    = ", object%selected_mci
    else
       write (u, "(3x,A)")  "Selected MCI    = [undefined]"
    end if
    if (object%selected_term > 0) then
       write (u, "(3x,A,I0)")  "Selected term   = ", object%selected_term
       call object%term(object%selected_term)%write (u, 1)
    else
       write (u, "(3x,A)")  "Selected term   = [undefined]"
    end if
  end subroutine decay_root_write

  recursive module subroutine decay_gen_write_process_instances &
       (decay, unit, verbose)
    class(decay_gen_t), intent(in) :: decay
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    logical :: verb
    verb = .true.;  if (present (verbose))  verb = verbose
    if (associated (decay%process_instance)) then
       if (verb) then
          call decay%process_instance%write (unit)
       else
          call decay%process_instance%write_header (unit)
       end if
    end if
    if (decay%selected_term > 0) then
       call decay%term(decay%selected_term)%write_process_instances (unit, verb)
    end if
  end subroutine decay_gen_write_process_instances

  recursive module subroutine decay_gen_init (decay, term_config)
    class(decay_gen_t), intent(out) :: decay
    type(decay_term_config_t), dimension(:), intent(in), target :: term_config
    integer :: i
    allocate (decay%term (size (term_config)))
    do i = 1, size (decay%term)
       call decay%term(i)%init (term_config(i))
    end do
  end subroutine decay_gen_init

  module subroutine decay_root_init (decay_root, config, process_instance)
    class(decay_root_t), intent(out) :: decay_root
    type(decay_root_config_t), intent(in), target :: config
    type(process_instance_t), intent(in), target :: process_instance
    call decay_root%base_init (config%term_config)
    decay_root%config => config
    decay_root%process_instance => process_instance
    call decay_root%make_term_rng (config%process)
    call decay_root%link_term_interactions ()
  end subroutine decay_root_init

  module subroutine decay_gen_set_mci (decay, i)
    class(decay_gen_t), intent(inout) :: decay
    integer, intent(in) :: i
    decay%selected_mci = i
  end subroutine decay_gen_set_mci

  module subroutine decay_gen_set_term (decay, i)
    class(decay_gen_t), intent(inout) :: decay
    integer, intent(in) :: i
    decay%selected_term = i
  end subroutine decay_gen_set_term

  module function decay_gen_get_mci (decay) result (i)
    class(decay_gen_t), intent(inout) :: decay
    integer :: i
    i = decay%selected_mci
  end function decay_gen_get_mci

  module function decay_gen_get_term (decay) result (i)
    class(decay_gen_t), intent(inout) :: decay
    integer :: i
    i = decay%selected_term
  end function decay_gen_get_term

  module subroutine decay_gen_make_term_rng (decay, process)
    class(decay_gen_t), intent(inout) :: decay
    type(process_t), intent(in), pointer :: process
    integer :: i
    do i = 1, size (decay%term)
       call decay%term(i)%make_rng (process)
    end do
  end subroutine decay_gen_make_term_rng

  recursive module subroutine decay_gen_link_term_interactions (decay)
    class(decay_gen_t), intent(inout) :: decay
    integer :: i
    type(interaction_t), pointer :: trace
    associate (instance => decay%process_instance)
      do i = 1, size (decay%term)
         trace => instance%get_trace_int_ptr (i)
         call decay%term(i)%link_interactions (trace)
      end do
    end associate
  end subroutine decay_gen_link_term_interactions

  module subroutine decay_root_select_chain (decay_root)
    class(decay_root_t), intent(inout) :: decay_root
    if (decay_root%selected_term > 0) then
       call decay_root%term(decay_root%selected_term)%select_chain ()
    else
       call msg_bug ("Decays: no term selected for parent process")
    end if
  end subroutine decay_root_select_chain

  module subroutine decay_root_generate (decay_root)
    class(decay_root_t), intent(inout) :: decay_root
    type(connected_state_t), pointer :: connected_state
    if (decay_root%selected_term > 0) then
       connected_state => decay_root%process_instance%get_connected_state_ptr &
            (decay_root%selected_term)
       call connected_state%normalize_matrix_by_trace ()
       call decay_root%term(decay_root%selected_term)%generate ()
    else
       call msg_bug ("Decays: no term selected for parent process")
    end if
  end subroutine decay_root_generate

  recursive module subroutine decay_config_write (object, unit, indent, verbose)
    class(decay_config_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    logical, intent(in), optional :: verbose
    integer :: u, ind
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    call write_indent (u, ind)
    write (u, "(1x,A)")  "Decay:"
    call object%write_header (unit, indent)
    call write_indent (u, ind)
    write (u, 2)  "branching ratio =", object%weight * 100
    call write_indent (u, ind)
    write (u, 1)  "partial width   =", object%integral
    call write_indent (u, ind)
    write (u, 1)  "error (abs)     =", object%abs_error
    call write_indent (u, ind)
    write (u, 1)  "error (rel)     =", object%rel_error
1   format (3x,A,ES19.12)
2   format (3x,A,F11.6,1x,'%')
    call object%write_terms (unit, indent, verbose)
  end subroutine decay_config_write

  recursive module subroutine decay_config_connect &
       (decay, process, model, process_stack, process_instance, var_list)
    class(decay_config_t), intent(out) :: decay
    type(process_t), intent(in), target :: process
    class(model_data_t), intent(in), target :: model
    type(process_stack_t), intent(in), optional :: process_stack
    type(process_instance_t), intent(in), optional, target :: process_instance
    type(var_list_t), intent(in), optional, target :: var_list
    real(default), dimension(:), allocatable :: integral_mci
    type(string_t) :: process_id
    integer :: i, n_mci
    call decay%decay_root_config_t%connect &
         (process, model, process_stack, var_list=var_list)
    process_id = process%get_id ()
    if (process%lab_is_cm ()) then
       call msg_fatal ("Decay process " // char (process_id) &
            // ": unusable because rest frame is fixed.")
    end if
    decay%integral = process%get_integral ()
    decay%abs_error = process%get_error ()
    if (present (var_list)) then
       call update (decay%integral, "integral(" // process_id // ")")
       call update (decay%abs_error, "error(" // process_id // ")")
    end if
    n_mci = process%get_n_mci ()
    allocate (integral_mci (n_mci))
    do i = 1, n_mci
       integral_mci(i) = process%get_integral_mci (i)
    end do
    call decay%mci_selector%init (integral_mci)
  contains
    subroutine update (var, var_name)
      real(default), intent(inout) :: var
      type(string_t), intent(in) :: var_name
      if (var_list%contains (var_name)) then
         var = var_list%get_rval (var_name)
      end if
    end subroutine update
  end subroutine decay_config_connect

  module subroutine decay_config_set_flv (decay, flv)
    class(decay_config_t), intent(inout) :: decay
    type(flavor_t), intent(in) :: flv
    decay%flv = flv
  end subroutine decay_config_set_flv

  recursive module subroutine decay_config_compute (decay)
    class(decay_config_t), intent(inout) :: decay
    call decay%decay_root_config_t%compute ()
    if (.not. vanishes (decay%integral)) then
       decay%rel_error = decay%abs_error / decay%integral
    else
       decay%rel_error = 0
    end if
  end subroutine decay_config_compute

  recursive module subroutine decay_final (object)
    class(decay_t), intent(inout) :: object
    integer :: i
    call object%base_final ()
    do i = 1, object%config%process%get_n_mci ()
       call object%process_instance%final_simulation (i)
    end do
    call object%process_instance%final ()
    deallocate (object%process_instance)
  end subroutine decay_final

  recursive module subroutine decay_write (object, unit, indent, recursive)
    class(decay_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent, recursive
    integer :: u, ind
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    call object%config%write (unit, indent, verbose = .false.)
    if (allocated (object%rng)) then
       call object%rng%write (u, ind + 1)
    end if
    call write_indent (u, ind)
    if (object%selected_mci > 0) then
       write (u, "(3x,A,I0)")  "Selected MCI    = ", object%selected_mci
    else
       write (u, "(3x,A)")  "Selected MCI    = [undefined]"
    end if
    call write_indent (u, ind)
    if (object%selected_term > 0) then
       write (u, "(3x,A,I0)")  "Selected term   = ", object%selected_term
       call object%term(object%selected_term)%write (u, ind + 1)
    else
       write (u, "(3x,A)")  "Selected term   = [undefined]"
    end if
  end subroutine decay_write

  recursive module subroutine decay_init (decay, config)
    class(decay_t), intent(out) :: decay
    type(decay_config_t), intent(in), target :: config
    integer :: i
    call decay%base_init (config%term_config)
    decay%config => config
    allocate (decay%process_instance)
    call decay%process_instance%init (decay%config%process)
    call decay%process_instance%setup_event_data (decay%config%model)
    do i = 1, decay%config%process%get_n_mci ()
       call decay%process_instance%init_simulation (i)
    end do
    call decay%config%process%make_rng (decay%rng)
    call decay%make_term_rng (decay%config%process)
  end subroutine decay_init

  recursive module subroutine decay_link_interactions (decay, i_prt, trace)
    class(decay_t), intent(inout) :: decay
    integer, intent(in) :: i_prt
    type(interaction_t), intent(in), target :: trace
    type(interaction_t), pointer :: beam_int
    integer :: n_in, n_vir
    beam_int => decay%process_instance%get_beam_int_ptr ()
    n_in = trace%get_n_in ()
    n_vir = trace%get_n_vir ()
    call beam_int%set_source_link (1, trace, &
         n_in + n_vir + i_prt)
    call decay%link_term_interactions ()
  end subroutine decay_link_interactions

  recursive module subroutine decay_select_chain (decay)
    class(decay_t), intent(inout) :: decay
    real(default) :: x
    integer :: i
    call decay%rng%generate (x)
    decay%selected_mci = decay%config%mci_selector%select (x)
    call decay%process_instance%choose_mci (decay%selected_mci)
    decay%selected_term = decay%process_instance%select_i_term ()
    do i = 1, size (decay%term)
       call decay%term(i)%select_chain ()
    end do
  end subroutine decay_select_chain

  recursive module subroutine decay_generate (decay)
    class(decay_t), intent(inout) :: decay
    type(isolated_state_t), pointer :: isolated_state
    integer :: i
    call decay%process_instance%receive_beam_momenta ()
    call decay%process_instance%generate_unweighted_event (decay%selected_mci)
    if (signal_is_pending ())  return
    call decay%process_instance%evaluate_event_data ()
    isolated_state => &
         decay%process_instance%get_isolated_state_ptr (decay%selected_term)
    call isolated_state%normalize_matrix_by_trace ()
    do i = 1, size (decay%term)
       call decay%term(i)%generate ()
       if (signal_is_pending ())  return
    end do
  end subroutine decay_generate

  module subroutine stable_config_final (object)
    class(stable_config_t), intent(inout) :: object
  end subroutine stable_config_final

  recursive module subroutine stable_config_write &
       (object, unit, indent, verbose)
    class(stable_config_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    logical, intent(in), optional :: verbose
    integer :: u, i, ind
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    call write_indent (u, ind)
    write (u, "(1x,'+',1x,A)", advance = "no")  "Stable:"
    write (u, "(1x,A)", advance = "no")  char (object%flv(1)%get_name ())
    do i = 2, size (object%flv)
       write (u, "(':',A)", advance = "no") &
            char (object%flv(i)%get_name ())
    end do
    write (u, *)
  end subroutine stable_config_write

  module subroutine stable_config_init (config, flv)
    class(stable_config_t), intent(out) :: config
    type(flavor_t), dimension(:), intent(in) :: flv
    integer, dimension (size (flv)) :: pdg
    logical, dimension (size (flv)) :: mask
    integer :: i
    pdg = flv%get_pdg ()
    mask(1) = .true.
    forall (i = 2 : size (pdg))
       mask(i) = all (pdg(i) /= pdg(1:i-1))
    end forall
    allocate (config%flv (count (mask)))
    config%flv = pack (flv, mask)
  end subroutine stable_config_init

  module subroutine stable_final (object)
    class(stable_t), intent(inout) :: object
  end subroutine stable_final

  module subroutine stable_write (object, unit, indent)
    class(stable_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    call object%config%write (unit, indent)
  end subroutine stable_write

  module subroutine stable_init (stable, config)
    class(stable_t), intent(out) :: stable
    type(stable_config_t), intent(in), target :: config
    stable%config => config
  end subroutine stable_init

  recursive module subroutine unstable_config_final (object)
    class(unstable_config_t), intent(inout) :: object
    integer :: i
    if (allocated (object%decay_config)) then
       do i = 1, size (object%decay_config)
          call object%decay_config(i)%final ()
       end do
    end if
  end subroutine unstable_config_final

  recursive module subroutine unstable_config_write &
       (object, unit, indent, verbose)
    class(unstable_config_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    logical, intent(in), optional :: verbose
    integer :: u, i, ind
    logical :: verb
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    verb = .true.;  if (present (verbose))  verb = verbose
    call write_indent (u, ind)
    write (u, "(1x,'+',1x,A,1x,A)")  "Unstable:", &
         char (object%flv%get_name ())
    call write_indent (u, ind)
    write (u, 1)  "total width =", object%integral
    call write_indent (u, ind)
    write (u, 1)  "error (abs) =", object%abs_error
    call write_indent (u, ind)
    write (u, 1)  "error (rel) =", object%rel_error
1   format (5x,A,ES19.12)
    if (verb .and. allocated (object%decay_config)) then
       do i = 1, size (object%decay_config)
          call object%decay_config(i)%write (u, ind + 1)
       end do
    end if
  end subroutine unstable_config_write

  module subroutine unstable_config_init (unstable, flv, set_decays, model)
    class(unstable_config_t), intent(out) :: unstable
    type(flavor_t), intent(in) :: flv
    logical, intent(in), optional :: set_decays
    class(model_data_t), intent(in), optional, target :: model
    type(string_t), dimension(:), allocatable :: decay
    unstable%flv = flv
    if (present (set_decays)) then
       call unstable%flv%get_decays (decay)
       call unstable%init_decays (decay, model)
    end if
  end subroutine unstable_config_init

  recursive module subroutine unstable_config_init_decays &
       (unstable, decay_id, model, process_stack, var_list)
    class(unstable_config_t), intent(inout) :: unstable
    type(string_t), dimension(:), intent(in) :: decay_id
    class(model_data_t), intent(in), target :: model
    type(process_stack_t), intent(in), optional :: process_stack
    type(var_list_t), intent(in), optional :: var_list
    integer :: i
    allocate (unstable%decay_config (size (decay_id)))
    do i = 1, size (decay_id)
       associate (decay => unstable%decay_config(i))
         if (present (process_stack)) then
            call decay%connect (process_stack%get_process_ptr (decay_id(i)), &
                 model, process_stack, var_list=var_list)
         else
            call decay%init (model, decay_id(i))
         end if
         call decay%set_flv (unstable%flv)
       end associate
    end do
  end subroutine unstable_config_init_decays

  module subroutine unstable_config_connect_decay (unstable, i, process, model)
    class(unstable_config_t), intent(inout) :: unstable
    integer, intent(in) :: i
    type(process_t), intent(in), target :: process
    class(model_data_t), intent(in), target :: model
    associate (decay => unstable%decay_config(i))
      call decay%connect (process, model)
    end associate
  end subroutine unstable_config_connect_decay

  recursive module subroutine unstable_config_compute (unstable)
    class(unstable_config_t), intent(inout) :: unstable
    integer :: i
    do i = 1, size (unstable%decay_config)
       call unstable%decay_config(i)%compute ()
    end do
    unstable%integral = sum (unstable%decay_config%integral)
    if (unstable%integral <= 0) then
       call unstable%write ()
       call msg_fatal ("Decay configuration: computed total width is zero")
    end if
    unstable%abs_error = sqrt (sum (unstable%decay_config%abs_error ** 2))
    unstable%rel_error = unstable%abs_error / unstable%integral
    call unstable%selector%init (unstable%decay_config%integral)
    do i = 1, size (unstable%decay_config)
       unstable%decay_config(i)%weight &
            = unstable%selector%get_weight (i)
    end do
  end subroutine unstable_config_compute

  recursive module subroutine unstable_final (object)
    class(unstable_t), intent(inout) :: object
    integer :: i
    if (allocated (object%decay)) then
       do i = 1, size (object%decay)
          call object%decay(i)%final ()
       end do
    end if
  end subroutine unstable_final

  recursive module subroutine unstable_write (object, unit, indent)
    class(unstable_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    integer :: u, ind
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    call object%config%write (u, ind, verbose=.false.)
    if (allocated (object%rng)) then
       call object%rng%write (u, ind + 2)
    end if
    call write_indent (u, ind)
    if (object%selected_decay > 0) then
       write (u, "(5x,A,I0)") "Sel. decay  = ", object%selected_decay
       call object%decay(object%selected_decay)%write (u, ind + 1)
    else
       write (u, "(5x,A)")  "Sel. decay  = [undefined]"
    end if
  end subroutine unstable_write

  recursive module subroutine unstable_write_process_instances &
       (unstable, unit, verbose)
    class(unstable_t), intent(in) :: unstable
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    if (unstable%selected_decay > 0) then
       call unstable%decay(unstable%selected_decay)% &
            write_process_instances (unit, verbose)
    end if
  end subroutine unstable_write_process_instances

  recursive module subroutine unstable_init (unstable, config)
    class(unstable_t), intent(out) :: unstable
    type(unstable_config_t), intent(in), target :: config
    integer :: i
    unstable%config => config
    allocate (unstable%decay (size (config%decay_config)))
    do i = 1, size (config%decay_config)
       call unstable%decay(i)%init (config%decay_config(i))
    end do
  end subroutine unstable_init

  recursive module subroutine unstable_link_interactions &
       (unstable, i_prt, trace)
    class(unstable_t), intent(inout) :: unstable
    integer, intent(in) :: i_prt
    type(interaction_t), intent(in), target :: trace
    integer :: i
    do i = 1, size (unstable%decay)
       call unstable%decay(i)%link_interactions (i_prt, trace)
    end do
  end subroutine unstable_link_interactions

  module subroutine unstable_import_rng (unstable, rng)
    class(unstable_t), intent(inout) :: unstable
    class(rng_t), intent(inout), allocatable :: rng
    call move_alloc (from = rng, to = unstable%rng)
  end subroutine unstable_import_rng

  recursive module subroutine unstable_select_chain (unstable)
    class(unstable_t), intent(inout) :: unstable
    real(default) :: x
    call unstable%rng%generate (x)
    unstable%selected_decay = unstable%config%selector%select (x)
    call unstable%decay(unstable%selected_decay)%select_chain ()
  end subroutine unstable_select_chain

  recursive module subroutine unstable_generate (unstable)
    class(unstable_t), intent(inout) :: unstable
    call unstable%decay(unstable%selected_decay)%generate ()
  end subroutine unstable_generate

  module subroutine decay_chain_final (object)
    class(decay_chain_t), intent(inout) :: object
    type(decay_chain_entry_t), pointer :: entry
    do while (associated (object%last))
       entry => object%last
       object%last => entry%previous
       call entry%final ()
       deallocate (entry)
    end do
    call object%correlated_trace%final ()
  end subroutine decay_chain_final

  module subroutine decay_chain_write (object, unit)
    class(decay_chain_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    call write_separator (u, 2)
    write (u, "(1x,A)")  "Decay chain:"
    call write_entries (object%last)
    call write_separator (u, 2)
    write (u, "(1x,A)")  "Evaluator (correlated trace of the decay chain):"
    call write_separator (u)
    call object%correlated_trace%write (u)
    call write_separator (u, 2)
  contains
    recursive subroutine write_entries (entry)
      type(decay_chain_entry_t), intent(in), pointer :: entry
      if (associated (entry)) then
         call write_entries (entry%previous)
         call write_separator (u, 2)
         write (u, "(1x,A,I0)")  "Decay #", entry%index
         call entry%config%write_header (u)
         write (u, "(3x,A,I0)")  "Selected MCI    = ", entry%selected_mci
         write (u, "(3x,A,I0)")  "Selected term   = ", entry%selected_term
         call entry%config%term_config(entry%selected_term)%write (u, indent=1)
         call entry%write (u)
      end if
    end subroutine write_entries
  end subroutine decay_chain_write

  module subroutine decay_chain_build (chain, decay_root)
    class(decay_chain_t), intent(inout), target :: chain
    type(decay_root_t), intent(in) :: decay_root
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask
    type(interaction_t), pointer :: int_last_decay
    call chain%final ()
    if (decay_root%selected_term > 0) then
       chain%process_instance => decay_root%process_instance
       chain%selected_term = decay_root%selected_term
       call chain%build_term_entries (decay_root%term(decay_root%selected_term))
    end if
    int_last_decay => chain%last%get_matrix_int_ptr ()
    allocate (qn_mask (int_last_decay%get_n_tot ()))
    call qn_mask%init (mask_f = .true., mask_c = .true., mask_h = .true.)
    call chain%correlated_trace%init_qn_sum (int_last_decay, qn_mask)
  end subroutine decay_chain_build

  recursive module subroutine decay_chain_build_term_entries (chain, term)
    class(decay_chain_t), intent(inout) :: chain
    type(decay_term_t), intent(in) :: term
    integer :: i
    do i = 1, size (term%particle_out)
       select type (unstable => term%particle_out(i)%c)
       type is (unstable_t)
          if (unstable%selected_decay > 0) then
             call chain%build_decay_entries &
                  (unstable%decay(unstable%selected_decay))
          end if
       end select
    end do
  end subroutine decay_chain_build_term_entries

  recursive module subroutine decay_chain_build_decay_entries (chain, decay)
    class(decay_chain_t), intent(inout) :: chain
    type(decay_t), intent(in) :: decay
    type(decay_chain_entry_t), pointer :: entry
    type(connected_state_t), pointer :: previous_state
    type(isolated_state_t), pointer :: current_decay
    type(helicity_t) :: hel
    type(quantum_numbers_t) :: qn_filter_conn
    allocate (entry)
    if (associated (chain%last)) then
       entry%previous => chain%last
       entry%index = entry%previous%index + 1
       previous_state => entry%previous%connected_state_t
    else
       entry%index = 1
       previous_state => &
            chain%process_instance%get_connected_state_ptr (chain%selected_term)
    end if
    entry%config => decay%config
    entry%selected_mci = decay%selected_mci
    entry%selected_term = decay%selected_term
    current_decay => decay%process_instance%get_isolated_state_ptr &
         (decay%selected_term)
    call entry%setup_connected_trace &
         (current_decay, previous_state%get_trace_int_ptr (), resonant=.true.)
    if (entry%config%flv%has_decay_helicity ()) then
       call hel%init (entry%config%flv%get_decay_helicity ())
       call qn_filter_conn%init (hel)
       call entry%setup_connected_matrix &
            (current_decay, previous_state%get_matrix_int_ptr (), &
            resonant=.true., qn_filter_conn = qn_filter_conn)
       call entry%setup_connected_flows &
            (current_decay, previous_state%get_flows_int_ptr (), &
            resonant=.true., qn_filter_conn = qn_filter_conn)
    else
       call entry%setup_connected_matrix &
            (current_decay, previous_state%get_matrix_int_ptr (), &
            resonant=.true.)
       call entry%setup_connected_flows &
            (current_decay, previous_state%get_flows_int_ptr (), &
            resonant=.true.)
    end if
    chain%last => entry
    call chain%build_term_entries (decay%term(decay%selected_term))
  end subroutine decay_chain_build_decay_entries

  module subroutine decay_chain_evaluate (chain)
    class(decay_chain_t), intent(inout) :: chain
    call evaluate (chain%last)
    call chain%correlated_trace%receive_momenta ()
    call chain%correlated_trace%evaluate ()
  contains
    recursive subroutine evaluate (entry)
      type(decay_chain_entry_t), intent(inout), pointer :: entry
      if (associated (entry)) then
         call evaluate (entry%previous)
         call entry%receive_kinematics ()
         call entry%evaluate_trace ()
         call entry%evaluate_event_data ()
      end if
    end subroutine evaluate
  end subroutine decay_chain_evaluate

  module function decay_chain_get_probability (chain) result (x)
    class(decay_chain_t), intent(in) :: chain
    real(default) :: x
    x = real (chain%correlated_trace%get_matrix_element (1))
  end function decay_chain_get_probability

  module subroutine evt_decay_write_name (evt, unit)
    class(evt_decay_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Event transform: partonic decays"
  end subroutine evt_decay_write_name

  module subroutine evt_decay_write (evt, unit, verbose, more_verbose, testflag)
    class(evt_decay_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, more_verbose, testflag
    logical :: verb, verb2
    integer :: u
    u = given_output_unit (unit)
    verb = .true.;  if (present (verbose))  verb = verbose
    verb2 = .false.;  if (present (more_verbose))  verb2 = more_verbose
    call write_separator (u, 2)
    call evt%write_name (u)
    call write_separator (u, 2)
    call evt%base_write (u, testflag = testflag)
    if (associated (evt%var_list)) then
       call write_separator (u)
       write (u, "(1x,A)")  "Variable list for simulation: &
            &[associated, not shown]"
    end if
    if (verb) then
       call write_separator (u)
       call evt%decay_root%write (u)
       if (verb2) then
          call evt%decay_chain%write (u)
          call evt%decay_root%write_process_instances (u, verb)
       end if
    else
       call write_separator (u, 2)
    end if
  end subroutine evt_decay_write

  module subroutine evt_decay_set_var_list (evt, var_list)
    class(evt_decay_t), intent(inout) :: evt
    type(var_list_t), intent(in), target :: var_list
    evt%var_list => var_list
  end subroutine evt_decay_set_var_list

  module subroutine evt_decay_connect &
       (evt, process_instance, model, process_stack)
    class(evt_decay_t), intent(inout), target :: evt
    type(process_instance_t), intent(in), target :: process_instance
    class(model_data_t), intent(in), target :: model
    type(process_stack_t), intent(in), optional :: process_stack
    call evt%base_connect (process_instance, model)
    if (associated (evt%var_list)) then
       call evt%decay_root_config%connect (process_instance%process, &
            model, process_stack, process_instance, evt%var_list)
    else
       call evt%decay_root_config%connect (process_instance%process, &
            model, process_stack, process_instance)
    end if
    call evt%decay_root_config%compute ()
    call evt%decay_root%init (evt%decay_root_config, evt%process_instance)
  end subroutine evt_decay_connect

  module subroutine evt_decay_prepare_new_event (evt, i_mci, i_term)
    class(evt_decay_t), intent(inout) :: evt
    integer, intent(in) :: i_mci, i_term
    call evt%reset ()
    evt%decay_root%selected_mci = i_mci
    evt%decay_root%selected_term = i_term
    call evt%decay_root%select_chain ()
    call evt%decay_chain%build (evt%decay_root)
  end subroutine evt_decay_prepare_new_event

  module subroutine evt_decay_generate_weighted (evt, probability)
    class(evt_decay_t), intent(inout) :: evt
    real(default), intent(inout) :: probability
    call evt%decay_root%generate ()
    if (signal_is_pending ())  return
    call evt%decay_chain%evaluate ()
    probability = evt%decay_chain%get_probability ()
  end subroutine evt_decay_generate_weighted

  module subroutine evt_decay_make_particle_set &
       (evt, factorization_mode, keep_correlations, r)
    class(evt_decay_t), intent(inout) :: evt
    integer, intent(in) :: factorization_mode
    logical, intent(in) :: keep_correlations
    real(default), dimension(:), intent(in), optional :: r
    type(interaction_t), pointer :: int_matrix, int_flows
    type(decay_chain_entry_t), pointer :: last_entry
    last_entry => evt%decay_chain%last
    int_matrix => last_entry%get_matrix_int_ptr ()
    int_flows  => last_entry%get_flows_int_ptr ()
    call evt%factorize_interactions (int_matrix, int_flows, &
         factorization_mode, keep_correlations, r)
    call evt%tag_incoming ()
  end subroutine evt_decay_make_particle_set

  module subroutine pacify_decay (evt)
    class(evt_decay_t), intent(inout) :: evt
    call pacify_decay_gen (evt%decay_root)
  end subroutine pacify_decay

  recursive module subroutine pacify_decay_gen (decay)
    class(decay_gen_t), intent(inout) :: decay
    if (associated (decay%process_instance)) then
       call pacify (decay%process_instance)
    end if
    if (decay%selected_term > 0) then
       call pacify_term (decay%term(decay%selected_term))
    end if
  end subroutine pacify_decay_gen

  recursive module subroutine pacify_term (term)
    class(decay_term_t), intent(inout) :: term
    integer :: i
    do i = 1, size (term%particle_out)
       select type (unstable => term%particle_out(i)%c)
       type is (unstable_t);  call pacify_unstable (unstable)
       end select
    end do
  end subroutine pacify_term

  recursive module subroutine pacify_unstable (unstable)
    class(unstable_t), intent(inout) :: unstable
    if (unstable%selected_decay > 0) then
       call pacify_decay_gen (unstable%decay(unstable%selected_decay))
    end if
  end subroutine pacify_unstable

  module subroutine init_test_case1 (unstable, i, flv, integral, relerr, model)
    class(unstable_config_t), intent(inout) :: unstable
    integer, intent(in) :: i
    type(flavor_t), dimension(:,:), intent(in) :: flv
    real(default), intent(in) :: integral
    real(default), intent(in) :: relerr
    class(model_data_t), intent(in), target :: model
    associate (decay => unstable%decay_config(i))
      allocate (decay%term_config (1))
      call decay%init_term (1, flv, stable = [.true., .true.], model=model)
      decay%integral = integral
      decay%abs_error = integral * relerr
    end associate
  end subroutine init_test_case1

  module subroutine init_test_case2 (unstable, flv1, flv21, flv22, model)
    class(unstable_config_t), intent(inout) :: unstable
    type(flavor_t), dimension(:,:), intent(in) :: flv1, flv21, flv22
    class(model_data_t), intent(in), target :: model
    associate (decay => unstable%decay_config(1))
      decay%integral = 1.e-3_default
      decay%abs_error = decay%integral * .01_default

      allocate (decay%term_config (1))
      call decay%init_term (1, flv1, stable = [.false., .true.], model=model)

      select type (w => decay%term_config(1)%prt(1)%c)
      type is (unstable_config_t)

         associate (w_decay => w%decay_config(1))
           w_decay%integral = 2._default
           allocate (w_decay%term_config (1))
           call w_decay%init_term (1, flv21, stable = [.true., .true.], &
                model=model)
         end associate
         associate (w_decay => w%decay_config(2))
           w_decay%integral = 1._default
           allocate (w_decay%term_config (1))
           call w_decay%init_term (1, flv22, stable = [.true., .true.], &
                model=model)
         end associate
         call w%compute ()

      end select
    end associate
  end subroutine init_test_case2


end submodule decays_s

