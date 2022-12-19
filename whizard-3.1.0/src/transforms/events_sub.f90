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

submodule (events) events_s

  use debug_master, only: debug_on
  use constants, only: one
  use io_units
  use format_utils, only: pac_fmt, write_separator
  use format_defs, only: FMT_12, FMT_19
  use numeric_utils
  use pcm, only: pcm_nlo_workspace_t
  use decays
  use evt_nlo
  use phs_fks, only: SQRTS_FIXED, SQRTS_VAR

  implicit none

contains

  module subroutine event_config_write (object, unit, show_expressions)
    class(event_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_expressions
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A,L1)")  "Unweighted         = ", object%unweighted
    write (u, "(3x,A,A)")   "Normalization      = ", &
         char (event_normalization_string (object%norm_mode))
    write (u, "(3x,A)", advance="no")  "Helicity handling  = "
    select case (object%factorization_mode)
    case (FM_IGNORE_HELICITY)
       write (u, "(A)")  "drop"
    case (FM_SELECT_HELICITY)
       write (u, "(A)")  "select"
    case (FM_FACTOR_HELICITY)
       write (u, "(A)")  "factorize"
    end select
    write (u, "(3x,A,L1)")  "Keep correlations  = ", object%keep_correlations
    if (object%colorize_subevt) then
       write (u, "(3x,A,L1)")  "Colorize subevent  = ", object%colorize_subevt
    end if
    if (.not. nearly_equal (object%safety_factor, one)) then
       write (u, "(3x,A," // FMT_12 // ")")  &
            "Safety factor      = ", object%safety_factor
    end if
    if (present (show_expressions)) then
       if (show_expressions) then
          if (allocated (object%ef_selection)) then
             call write_separator (u)
             write (u, "(3x,A)") "Event selection expression:"
             call object%ef_selection%write (u)
          end if
          if (allocated (object%ef_reweight)) then
             call write_separator (u)
             write (u, "(3x,A)") "Event reweighting expression:"
             call object%ef_reweight%write (u)
          end if
          if (allocated (object%ef_analysis)) then
             call write_separator (u)
             write (u, "(3x,A)") "Analysis expression:"
             call object%ef_analysis%write (u)
          end if
       end if
    end if
  end subroutine event_config_write

  module subroutine event_clone (event, event_new)
    class(event_t), intent(in), target :: event
    class(event_t), intent(out), target:: event_new
    type(string_t) :: id
    integer :: num_id
    event_new%config = event%config
    event_new%process => event%process
    event_new%instance => event%instance
    if (allocated (event%rng)) &
         allocate(event_new%rng, source=event%rng)
    event_new%selected_i_mci = event%selected_i_mci
    event_new%selected_i_term = event%selected_i_term
    event_new%selected_channel = event%selected_channel
    event_new%is_complete = event%is_complete
    event_new%transform_first => event%transform_first
    event_new%transform_last => event%transform_last
    event_new%selection_evaluated = event%selection_evaluated
    event_new%passed = event%passed
    if (allocated (event%alpha_qcd_forced)) &
         allocate(event_new%alpha_qcd_forced, source=event%alpha_qcd_forced)
    if (allocated (event%scale_forced)) &
         allocate(event_new%scale_forced, source=event%scale_forced)
    event_new%reweight = event%reweight
    event_new%analysis_flag = event%analysis_flag
    event_new%i_event = event%i_event
    id = event_new%process%get_id ()
    if (id /= "")  call event_new%expr%set_process_id (id)
    num_id = event_new%process%get_num_id ()
    if (num_id /= 0)  call event_new%expr%set_process_num_id (num_id)
    call event_new%expr%setup_vars (event_new%process%get_sqrts ())
    call event_new%expr%link_var_list (event_new%process%get_var_list_ptr ())
  end subroutine event_clone

  module subroutine event_final (object)
    class(event_t), intent(inout) :: object
    class(evt_t), pointer :: evt
    if (allocated (object%rng))  call object%rng%final ()
    call object%expr%final ()
    do while (associated (object%transform_first))
       evt => object%transform_first
       object%transform_first => evt%next
       call evt%final ()
       deallocate (evt)
    end do
  end subroutine event_final

  module subroutine event_write (object, unit, show_process, &
       show_transforms, show_decay, verbose, testflag)
    class(event_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_process, show_transforms, show_decay
    logical, intent(in), optional :: verbose
    logical, intent(in), optional :: testflag
    logical :: prc, trans, dec, verb
    class(evt_t), pointer :: evt
    character(len=7) :: fmt
    integer :: u, i
    call pac_fmt (fmt, FMT_19, FMT_12, testflag)
    u = given_output_unit (unit)
    prc = .true.;  if (present (show_process))  prc = show_process
    trans = .true.;  if (present (show_transforms))  trans = show_transforms
    dec = .true.;  if (present (show_decay))  dec = show_decay
    verb = .false.;  if (present (verbose))  verb = verbose
    call write_separator (u, 2)
    write (u, "(1x,A)", advance="no")  "Event"
    if (object%has_index ()) then
       write (u, "(1x,'#',I0)", advance="no")  object%get_index ()
    end if
    if (object%is_complete) then
       write (u, *)
    else
       write (u, "(1x,A)")  "[incomplete]"
    end if
    call write_separator (u)
    call object%config%write (u)
    if (object%sqme_ref_is_known () .or. object%weight_ref_is_known ()) then
       call write_separator (u)
    end if
    if (object%sqme_ref_is_known ()) then
       write (u, "(3x,A," // fmt // ")") &
            "Squared matrix el. (ref) = ", object%get_sqme_ref ()
       if (object%sqme_alt_is_known ()) then
          do i = 1, object%get_n_alt ()
             write (u, "(5x,A," // fmt // ",1x,I0)")  &
                  "alternate sqme   = ", object%get_sqme_alt(i), i
          end do
       end if
    end if
    if (object%sqme_prc_is_known ()) &
       write (u, "(3x,A," // fmt // ")") &
            "Squared matrix el. (prc) = ", object%get_sqme_prc ()

    if (object%weight_ref_is_known ()) then
       write (u, "(3x,A," // fmt // ")") &
            "Event weight (ref)       = ", object%get_weight_ref ()
       if (object%weight_alt_is_known ()) then
          do i = 1, object%get_n_alt ()
             write (u, "(5x,A," // fmt // ",1x,I0)")  &
                  "alternate weight       = ", object%get_weight_alt(i), i
          end do
       end if
    end if
    if (object%weight_prc_is_known ()) &
       write (u, "(3x,A," // fmt // ")") &
            "Event weight (prc)       = ", object%get_weight_prc ()

    if (object%selected_i_mci /= 0) then
       call write_separator (u)
       write (u, "(3x,A,I0)")  "Selected MCI group = ", object%selected_i_mci
       write (u, "(3x,A,I0)")  "Selected term      = ", object%selected_i_term
       write (u, "(3x,A,I0)")  "Selected channel   = ", object%selected_channel
    end if
    if (object%selection_evaluated) then
       call write_separator (u)
       write (u, "(3x,A,L1)")  "Passed selection   = ", object%passed
       if (object%passed) then
          write (u, "(3x,A," // fmt // ")") &
               "Reweighting factor = ", object%reweight
          write (u, "(3x,A,L1)") &
               "Analysis flag      = ", object%analysis_flag
       end if
    end if
    if (associated (object%instance)) then
       if (prc) then
          if (verb) then
             call object%instance%write (u, testflag)
          else
             call object%instance%write_header (u)
          end if
       end if
       if (trans) then
          evt => object%transform_first
          do while (associated (evt))
             select type (evt)
             type is (evt_decay_t)
                call evt%write (u, verbose = dec, more_verbose = verb, &
                     testflag = testflag)
             class default
                call evt%write (u, verbose = verb, testflag = testflag)
             end select
             call write_separator (u, 2)
             evt => evt%next
          end do
       else
          call write_separator (u, 2)
       end if
       if (object%expr%subevt_filled) then
          call object%expr%write (u, pacified = testflag)
          call write_separator (u, 2)
       end if
    else
       call write_separator (u, 2)
       write (u, "(1x,A)")  "Process instance: [undefined]"
       call write_separator (u, 2)
    end if
  end subroutine event_write

  elemental module subroutine event_set_sigma (event, sigma)
    class(event_t), intent(inout) :: event
    real(default), intent(in) :: sigma
    event%config%sigma = sigma
  end subroutine event_set_sigma

  elemental module subroutine event_set_n (event, n)
    class(event_t), intent(inout) :: event
    integer, intent(in) :: n
    event%config%n = n
  end subroutine event_set_n

  module subroutine event_import_transform (event, evt)
    class(event_t), intent(inout) :: event
    class(evt_t), intent(inout), pointer :: evt
    event%transform_last%next => evt
    evt%previous => event%transform_last
    event%transform_last => evt
    evt => null ()
  end subroutine event_import_transform

  module subroutine event_connect &
       (event, process_instance, model, process_stack)
    class(event_t), intent(inout), target :: event
    type(process_instance_t), intent(in), target :: process_instance
    class(model_data_t), intent(in), target :: model
    type(process_stack_t), intent(in), optional :: process_stack
    type(string_t) :: id
    integer :: num_id
    class(evt_t), pointer :: evt
    event%process => process_instance%process
    event%instance => process_instance
    id = event%process%get_id ()
    if (id /= "")  call event%expr%set_process_id (id)
    num_id = event%process%get_num_id ()
    if (num_id /= 0)  call event%expr%set_process_num_id (num_id)
    call event%expr%setup_vars (event%process%get_sqrts ())
    call event%expr%link_var_list (event%process%get_var_list_ptr ())
    call event%process%make_rng (event%rng)
    evt => event%transform_first
    do while (associated (evt))
       call evt%connect (process_instance, model, process_stack)
       evt => evt%next
    end do
  end subroutine event_connect

  module subroutine event_set_selection (event, ef_selection)
    class(event_t), intent(inout) :: event
    class(expr_factory_t), intent(in) :: ef_selection
    allocate (event%config%ef_selection, source = ef_selection)
  end subroutine event_set_selection

  module subroutine event_set_reweight (event, ef_reweight)
    class(event_t), intent(inout) :: event
    class(expr_factory_t), intent(in) :: ef_reweight
    allocate (event%config%ef_reweight, source = ef_reweight)
  end subroutine event_set_reweight

  module subroutine event_set_analysis (event, ef_analysis)
    class(event_t), intent(inout) :: event
    class(expr_factory_t), intent(in) :: ef_analysis
    allocate (event%config%ef_analysis, source = ef_analysis)
  end subroutine event_set_analysis

  module subroutine event_setup_expressions (event)
    class(event_t), intent(inout), target :: event
    call event%expr%setup_selection (event%config%ef_selection)
    call event%expr%setup_analysis (event%config%ef_analysis)
    call event%expr%setup_reweight (event%config%ef_reweight)
    call event%expr%colorize (event%config%colorize_subevt)
  end subroutine event_setup_expressions

  module subroutine event_evaluate_transforms (event, r)
    class(event_t), intent(inout) :: event
    real(default), dimension(:), intent(in), optional :: r
    class(evt_t), pointer :: evt
    real(default) :: weight_over_sqme
    integer :: i_term, emitter, n_in
    logical :: failed_but_keep
    failed_but_keep = .false.
    if (debug_on) call msg_debug (D_TRANSFORMS, "event_evaluate_transforms")
    call event%discard_particle_set ()
    call event%check ()
    if (event%instance%is_complete_event ()) then
       i_term = event%instance%select_i_term ()
       event%selected_i_term = i_term
       evt => event%transform_first
       do while (associated (evt))
          call evt%prepare_new_event &
               (event%selected_i_mci, event%selected_i_term)
          evt => evt%next
       end do
       if (debug_on) call msg_debug &
            (D_TRANSFORMS, "Before event transformations")
       if (debug_on) call msg_debug &
            (D_TRANSFORMS, "event%weight_prc", event%weight_prc)
       if (debug_on) call msg_debug &
            (D_TRANSFORMS, "event%sqme_prc", event%sqme_prc)
       evt => event%transform_first
       do while (associated (evt))
          call print_transform_name_if_debug ()
          if (evt%only_weighted_events) then
             select type (evt)
             type is (evt_nlo_t)
                i_term = evt%get_i_term ()
                failed_but_keep = .not. evt%is_valid_event (i_term) &
                     .and. evt%keep_failed_events
                if (.not. any(evt%process_instance%term%passed .and. evt%process_instance%term%active) &
                    .and. .not. evt%keep_failed_events)  return
             end select
             if (abs (event%weight_prc) > 0._default) then
                weight_over_sqme = event%weight_prc / event%sqme_prc
                call evt%generate_weighted (event%sqme_prc)
                event%weight_prc = weight_over_sqme * event%sqme_prc
                select type (evt)
                type is (evt_nlo_t)
                   if (.not. evt%is_valid_event (i_term)) event%weight_prc = 0
                end select
             else
                if (.not. failed_but_keep) exit
             end if
          else
             call evt%generate_unweighted ()
          end if
          if (signal_is_pending ())  return
          select type (evt)
          type is (evt_nlo_t)
             if (evt%i_evaluation > 0) then
                emitter = evt%process_instance%kin(i_term)%emitter
                n_in = evt%process_instance%kin(i_term)%n_in
                select type (pcm_work => evt%process_instance%term(i_term)%pcm_work)
                type is (pcm_nlo_workspace_t)
                   if (emitter <= n_in .and. pcm_work%isr_kinematics%isr_mode == SQRTS_VAR) then
                      call evt%connected_set_real_IS_momenta ()
                   end if
                end select
             end if
          end select
          call evt%make_particle_set (event%config%factorization_mode, &
               event%config%keep_correlations)
          if (signal_is_pending ())  return
          if (.not. evt%particle_set_exists) exit
          evt => evt%next
       end do
       evt => event%transform_last
       if ((associated (evt) .and. evt%particle_set_exists) .or. failed_but_keep) then
          if (event%is_nlo ()) then
             select type (evt)
             type is (evt_nlo_t)
                evt%particle_set_nlo (event%i_event + 1) = evt%particle_set
                evt%i_evaluation = evt%i_evaluation + 1
                call event%link_particle_set &
                     (evt%particle_set_nlo(event%i_event + 1))
             end select
          else
             call event%link_particle_set (evt%particle_set)
          end if
       end if
       if (debug_on) call msg_debug &
            (D_TRANSFORMS, "After event transformations")
       if (debug_on) call msg_debug &
            (D_TRANSFORMS, "event%weight_prc", event%weight_prc)
       if (debug_on) call msg_debug &
            (D_TRANSFORMS, "event%sqme_prc", event%sqme_prc)
       if (debug_on) call msg_debug &
            (D_TRANSFORMS, "evt%particle_set_exists", evt%particle_set_exists)
    end if
  contains
    subroutine print_transform_name_if_debug ()
       if (debug_active (D_TRANSFORMS)) then
          print *, 'Current event transform: '
          call evt%write_name ()
       end if
    end subroutine print_transform_name_if_debug
  end subroutine event_evaluate_transforms

  module subroutine event_set_index (event, index)
    class(event_t), intent(inout) :: event
    integer, intent(in) :: index
    call event%expr%set_event_index (index)
  end subroutine event_set_index

  module subroutine event_increment_index (event, offset)
    class(event_t), intent(inout) :: event
    integer, intent(in), optional :: offset
    call event%expr%increment_event_index (offset)
  end subroutine event_increment_index

  module subroutine event_evaluate_expressions (event)
    class(event_t), intent(inout) :: event
    if (event%has_valid_particle_set ()) then
       call event%expr%fill_subevt (event%get_particle_set_ptr ())
    end if
    if (event%weight_ref_is_known ()) then
       call event%expr%set (weight_ref = event%get_weight_ref ())
    end if
    if (event%weight_prc_is_known ()) then
       call event%expr%set (weight_prc = event%get_weight_prc ())
    end if
    if (event%excess_prc_is_known ()) then
       call event%expr%set (excess_prc = event%get_excess_prc ())
    end if
    if (event%sqme_ref_is_known ()) then
       call event%expr%set (sqme_ref = event%get_sqme_ref ())
    end if
    if (event%sqme_prc_is_known ()) then
       call event%expr%set (sqme_prc = event%get_sqme_prc ())
    end if
    if (event%has_valid_particle_set ()) then
       call event%expr%evaluate &
            (event%passed, event%reweight, event%analysis_flag)
       event%selection_evaluated = .true.
    end if
  end subroutine event_evaluate_expressions

  module function event_passed_selection (event) result (flag)
    class(event_t), intent(in) :: event
    logical :: flag
    flag = event%passed
  end function event_passed_selection

  module subroutine event_store_alt_values (event)
    class(event_t), intent(inout) :: event
    if (event%weight_alt_is_known ()) then
       call event%expr%set (weight_alt = event%get_weight_alt ())
    end if
    if (event%sqme_alt_is_known ()) then
       call event%expr%set (sqme_alt = event%get_sqme_alt ())
    end if
  end subroutine event_store_alt_values

  module function event_is_nlo (event) result (is_nlo)
    logical :: is_nlo
    class(event_t), intent(in) :: event
    if (associated (event%instance)) then
       select type (pcm_work => event%instance%pcm_work)
       type is (pcm_nlo_workspace_t)
          is_nlo = pcm_work%is_fixed_order_nlo_events ()
       class default
          is_nlo = .false.
       end select
    else
       is_nlo = .false.
    end if
  end function event_is_nlo

  module subroutine event_reset_contents (event)
    class(event_t), intent(inout) :: event
    class(evt_t), pointer :: evt
    call event%base_reset_contents ()
    event%selected_i_mci = 0
    event%selected_i_term = 0
    event%selected_channel = 0
    event%is_complete = .false.
    call event%expr%reset_contents ()
    event%selection_evaluated = .false.
    event%passed = .false.
    event%analysis_flag = .false.
    if (associated (event%instance)) then
       call event%instance%reset (reset_mci = .true.)
    end if
    if (allocated (event%alpha_qcd_forced))  deallocate (event%alpha_qcd_forced)
    if (allocated (event%scale_forced))  deallocate (event%scale_forced)
    evt => event%transform_first
    do while (associated (evt))
       call evt%reset ()
       evt => evt%next
    end do
  end subroutine event_reset_contents

  module subroutine event_reset_index (event)
    class(event_t), intent(inout) :: event
    call event%expr%reset_event_index ()
  end subroutine event_reset_index

  module subroutine event_import_instance_results (event)
    class(event_t), intent(inout) :: event
    if (associated (event%instance)) then
       if (event%instance%has_evaluated_trace ()) then
          call event%set ( &
               sqme_prc = event%instance%get_sqme (), &
               weight_prc = event%instance%get_weight (), &
               excess_prc = event%instance%get_excess (), &
               n_dropped = event%instance%get_n_dropped () &
               )
       end if
    end if
  end subroutine event_import_instance_results

  module subroutine event_accept_sqme_ref (event)
    class(event_t), intent(inout) :: event
    if (event%sqme_ref_is_known ()) then
       call event%set (sqme_prc = event%get_sqme_ref ())
    end if
  end subroutine event_accept_sqme_ref

  module subroutine event_accept_sqme_prc (event)
    class(event_t), intent(inout) :: event
    if (event%sqme_prc_is_known ()) then
       call event%set (sqme_ref = event%get_sqme_prc ())
    end if
  end subroutine event_accept_sqme_prc

  module subroutine event_accept_weight_ref (event)
    class(event_t), intent(inout) :: event
    if (event%weight_ref_is_known ()) then
       call event%set (weight_prc = event%get_weight_ref ())
    end if
  end subroutine event_accept_weight_ref

  module subroutine event_accept_weight_prc (event)
    class(event_t), intent(inout) :: event
    if (event%weight_prc_is_known ()) then
       call event%set (weight_ref = event%get_weight_prc ())
    end if
  end subroutine event_accept_weight_prc

  module subroutine event_update_normalization (event, mode_ref)
    class(event_t), intent(inout) :: event
    integer, intent(in), optional :: mode_ref
    integer :: mode_old
    real(default) :: weight, excess
    if (present (mode_ref)) then
       mode_old = mode_ref
    else if (event%config%unweighted) then
       mode_old = NORM_UNIT
    else
       mode_old = NORM_SIGMA
    end if
    weight = event%get_weight_prc ()
    call event_normalization_update (weight, &
         event%config%sigma, event%config%n, &
         mode_new = event%config%norm_mode, &
         mode_old = mode_old)
    call event%set_weight_prc (weight)
    excess = event%get_excess_prc ()
    call event_normalization_update (excess, &
         event%config%sigma, event%config%n, &
         mode_new = event%config%norm_mode, &
         mode_old = mode_old)
    call event%set_excess_prc (excess)
  end subroutine event_update_normalization

  module subroutine event_check (event)
    class(event_t), intent(inout) :: event
    event%is_complete = event%has_valid_particle_set () &
         .and. event%sqme_ref_is_known () &
         .and. event%sqme_prc_is_known () &
         .and. event%weight_ref_is_known () &
         .and. event%weight_prc_is_known ()
    if (event%get_n_alt () /= 0) then
       event%is_complete = event%is_complete &
            .and. event%sqme_alt_is_known () &
            .and. event%weight_alt_is_known ()
    end if
  end subroutine event_check

  module subroutine event_generate (event, i_mci, r, i_nlo)
    class(event_t), intent(inout) :: event
    integer, intent(in) :: i_mci
    real(default), dimension(:), intent(in), optional :: r
    integer, intent(in), optional :: i_nlo
    logical :: generate_new
    generate_new = .true.
    if (present (i_nlo)) generate_new = (i_nlo == 1)
    if (generate_new) call event%reset_contents ()
    event%selected_i_mci = i_mci
    if (event%config%unweighted) then
       call event%instance%generate_unweighted_event (i_mci)
       if (signal_is_pending ()) return
       call event%instance%evaluate_event_data ()
       call event%instance%normalize_weight ()
    else
       if (generate_new) &
          call event%instance%generate_weighted_event (i_mci)
       if (signal_is_pending ()) return
       call event%instance%evaluate_event_data ()
    end if
    event%selected_channel = event%instance%get_channel ()
    call event%import_instance_results ()
    call event%accept_sqme_prc ()
    call event%update_normalization ()
    call event%accept_weight_prc ()
    call event%evaluate_transforms (r)
    if (signal_is_pending ())  return
    call event%check ()
  end subroutine event_generate

  module subroutine event_get_hard_particle_set (event, pset)
    class(event_t), intent(in) :: event
    type(particle_set_t), intent(out) :: pset
    class(evt_t), pointer :: evt
    evt => event%transform_first
    pset = evt%particle_set
  end subroutine event_get_hard_particle_set

  module subroutine event_select (event, i_mci, i_term, channel)
    class(event_t), intent(inout) :: event
    integer, intent(in) :: i_mci, i_term, channel
    if (associated (event%instance)) then
       event%selected_i_mci = i_mci
       event%selected_i_term = i_term
       event%selected_channel = channel
    else
       event%selected_i_mci = 0
       event%selected_i_term = 0
       event%selected_channel = 0
    end if
  end subroutine event_select

  module subroutine event_set_hard_particle_set (event, particle_set)
    class(event_t), intent(inout) :: event
    type(particle_set_t), intent(in) :: particle_set
    class(evt_t), pointer :: evt
    evt => event%transform_first
    call evt%set_particle_set (particle_set, &
         event%selected_i_mci, event%selected_i_term)
    call event%link_particle_set (evt%particle_set)
    evt => evt%next
    do while (associated (evt))
       call evt%reset ()
       evt => evt%next
    end do
  end subroutine event_set_hard_particle_set

  module subroutine event_set_alpha_qcd_forced (event, alpha_qcd)
    class(event_t), intent(inout) :: event
    real(default), intent(in) :: alpha_qcd
    if (allocated (event%alpha_qcd_forced)) then
       event%alpha_qcd_forced = alpha_qcd
    else
       allocate (event%alpha_qcd_forced, source = alpha_qcd)
    end if
  end subroutine event_set_alpha_qcd_forced

  module subroutine event_set_scale_forced (event, scale)
    class(event_t), intent(inout) :: event
    real(default), intent(in) :: scale
    if (allocated (event%scale_forced)) then
       event%scale_forced = scale
    else
       allocate (event%scale_forced, source = scale)
    end if
  end subroutine event_set_scale_forced

  module subroutine event_recalculate (event, update_sqme, weight_factor, &
       recover_beams, recover_phs, check_match, success)
    class(event_t), intent(inout) :: event
    logical, intent(in) :: update_sqme
    real(default), intent(in), optional :: weight_factor
    logical, intent(in), optional :: recover_beams
    logical, intent(in), optional :: recover_phs
    logical, intent(in), optional :: check_match
    logical, intent(out), optional :: success
    type(particle_set_t), pointer :: particle_set
    integer :: i_mci, i_term, channel
    logical :: rec_phs_mci
    rec_phs_mci = .true.;  if (present (recover_phs))  rec_phs_mci = recover_phs
    if (present (success))  success = .false.
    if (event%has_valid_particle_set ()) then
       particle_set => event%get_particle_set_ptr ()
       i_mci = event%selected_i_mci
       i_term = event%selected_i_term
       channel = event%selected_channel
       if (i_mci == 0 .or. i_term == 0 .or. channel == 0) then
          call msg_bug ("Event: recalculate: undefined selection parameters")
       end if
       call event%instance%choose_mci (i_mci)
       call event%instance%set_trace &
            (particle_set, i_term, recover_beams, check_match, success)
       if (present (success)) then
          if (.not. success)  return
       end if
       if (allocated (event%alpha_qcd_forced)) then
          call event%instance%set_alpha_qcd_forced &
               (i_term, event%alpha_qcd_forced)
       end if
       call event%instance%recover (channel, i_term, &
            update_sqme, rec_phs_mci, event%scale_forced)
       if (signal_is_pending ())  return
       if (update_sqme .and. present (weight_factor)) then
          call event%instance%evaluate_event_data &
               (weight = event%instance%get_sqme () * weight_factor)
       else if (event%weight_ref_is_known ()) then
          call event%instance%evaluate_event_data &
               (weight = event%get_weight_ref ())
       else if (rec_phs_mci) then
          call event%instance%recover_event ()
          if (signal_is_pending ())  return
          call event%instance%evaluate_event_data ()
          if (event%config%unweighted) then
             call event%instance%normalize_weight ()
          end if
       end if
       if (signal_is_pending ())  return
       if (update_sqme) then
          call event%import_instance_results ()
       else
          call event%accept_sqme_ref ()
          call event%accept_weight_ref ()
       end if
    else
       call msg_bug ("Event: can't recalculate, particle set is undefined")
    end if
  end subroutine event_recalculate

  module function event_get_process_ptr (event) result (ptr)
    class(event_t), intent(in) :: event
    type(process_t), pointer :: ptr
    ptr => event%process
  end function event_get_process_ptr

  module function event_get_process_instance_ptr (event) result (ptr)
    class(event_t), intent(in) :: event
    type(process_instance_t), pointer :: ptr
    ptr => event%instance
  end function event_get_process_instance_ptr

  module function event_get_model_ptr (event) result (model)
    class(event_t), intent(in) :: event
    class(model_data_t), pointer :: model
    if (associated (event%process)) then
       model => event%process%get_model_ptr ()
    else
       model => null ()
    end if
  end function event_get_model_ptr

  module function event_get_i_mci (event) result (i_mci)
    class(event_t), intent(in) :: event
    integer :: i_mci
    i_mci = event%selected_i_mci
  end function event_get_i_mci

  module function event_get_i_term (event) result (i_term)
    class(event_t), intent(in) :: event
    integer :: i_term
    i_term = event%selected_i_term
  end function event_get_i_term

  module function event_get_channel (event) result (channel)
    class(event_t), intent(in) :: event
    integer :: channel
    channel = event%selected_channel
  end function event_get_channel

  module function event_has_transform (event) result (flag)
    class(event_t), intent(in) :: event
    logical :: flag
    if (associated (event%transform_first)) then
       flag = associated (event%transform_first%next)
    else
       flag = .false.
    end if
  end function event_has_transform

  elemental module function event_get_norm_mode (event) result (norm_mode)
    class(event_t), intent(in) :: event
    integer :: norm_mode
    norm_mode = event%config%norm_mode
  end function event_get_norm_mode

  module function event_get_kinematical_weight (event) result (f)
    class(event_t), intent(in) :: event
    real(default) :: f
    if (event%sqme_ref_is_known () .and. event%weight_ref_is_known () &
         .and. abs (event%get_sqme_ref ()) > 0) then
       f = event%get_weight_ref () / event%get_sqme_ref ()
    else
       f = 0
    end if
  end function event_get_kinematical_weight

  module function event_has_index (event) result (flag)
    class(event_t), intent(in) :: event
    logical :: flag
    flag = event%expr%has_event_index ()
  end function event_has_index

  module function event_get_index (event) result (index)
    class(event_t), intent(in) :: event
    integer :: index
    index = event%expr%get_event_index ()
  end function event_get_index

  module function event_get_fac_scale (event) result (fac_scale)
    class(event_t), intent(in) :: event
    real(default) :: fac_scale
    fac_scale = event%instance%get_fac_scale (event%selected_i_term)
  end function event_get_fac_scale

  module function event_get_alpha_s (event) result (alpha_s)
    class(event_t), intent(in) :: event
    real(default) :: alpha_s
    alpha_s = event%instance%get_alpha_s (event%selected_i_term)
  end function event_get_alpha_s

  module function event_get_sqrts (event) result (sqrts)
    class(event_t), intent(in) :: event
    real(default) :: sqrts
    sqrts = event%instance%get_sqrts ()
  end function event_get_sqrts

  module function event_get_polarization (event) result (pol)
    class(event_t), intent(in) :: event
    real(default), dimension(:), allocatable :: pol
    pol = event%instance%get_polarization ()
  end function event_get_polarization

  module function event_get_beam_file (event) result (file)
    class(event_t), intent(in) :: event
    type(string_t) :: file
    file = event%instance%get_beam_file ()
  end function event_get_beam_file

  module function event_get_process_name (event) result (name)
    class(event_t), intent(in) :: event
    type(string_t) :: name
    name = event%instance%get_process_name ()
  end function event_get_process_name

  elemental module function event_get_actual_calls_total (event) result (n)
    class(event_t), intent(in) :: event
    integer :: n
    if (associated (event%instance)) then
       n = event%instance%get_actual_calls_total ()
    else
       n = 0
    end if
  end function event_get_actual_calls_total

  module subroutine pacify_event (event)
    class(event_t), intent(inout) :: event
    class(evt_t), pointer :: evt
    call event%pacify_particle_set ()
    if (event%expr%subevt_filled)  call pacify (event%expr)
    evt => event%transform_first
    do while (associated (evt))
       select type (evt)
       type is (evt_decay_t);  call pacify (evt)
       end select
       evt => evt%next
    end do
  end subroutine pacify_event


end submodule events_s

