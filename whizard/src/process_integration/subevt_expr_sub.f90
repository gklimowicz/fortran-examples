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

submodule (subevt_expr) subevt_expr_s

  use constants, only: zero, one
  use io_units
  use format_utils, only: write_separator
  use diagnostics

  implicit none

contains

  module subroutine subevt_expr_write (object, unit, pacified)
    class(subevt_expr_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: pacified
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Local variables:"
    call write_separator (u)
    call object%var_list%write (u, follow_link=.false., &
         pacified = pacified)
    call write_separator (u)
    if (object%subevt_filled) then
       call object%subevt_t%write (u, pacified = pacified)
       if (object%has_selection) then
          call write_separator (u)
          write (u, "(1x,A)")  "Selection expression:"
          call write_separator (u)
          call object%selection%write (u)
       end if
    else
       write (u, "(1x,A)")  "subevt: [undefined]"
    end if
  end subroutine subevt_expr_write

  module subroutine subevt_expr_final (object)
    class(subevt_expr_t), intent(inout) :: object
    call object%var_list%final ()
    if (object%has_selection) then
       call object%selection%final ()
    end if
  end subroutine subevt_expr_final

  module subroutine subevt_expr_setup_vars (expr, sqrts)
    class(subevt_expr_t), intent(inout), target :: expr
    real(default), intent(in) :: sqrts
    call expr%var_list%final ()
    call expr%var_list%append_real (var_str ("sqrts"), sqrts, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_real_ptr (var_str ("sqrts_hat"), &
         expr%sqrts_hat, is_known = expr%subevt_filled, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_int_ptr (var_str ("n_in"), expr%n_in, &
         is_known = expr%subevt_filled, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_int_ptr (var_str ("n_out"), expr%n_out, &
         is_known = expr%subevt_filled, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_int_ptr (var_str ("n_tot"), expr%n_tot, &
         is_known = expr%subevt_filled, &
         locked = .true., verbose = .false., intrinsic = .true.)
  end subroutine subevt_expr_setup_vars

  module subroutine subevt_expr_setup_var_self (expr)
    class(subevt_expr_t), intent(inout), target :: expr
    if (.not. expr%var_list%contains (var_str ("@evt"))) then
       call expr%var_list%append_subevt_ptr &
            (var_str ("@evt"), expr%subevt_t, &
            is_known = expr%subevt_filled, &
            locked = .true., verbose = .false., intrinsic=.true.)
    end if
  end subroutine subevt_expr_setup_var_self

  module subroutine subevt_expr_link_var_list (expr, var_list)
    class(subevt_expr_t), intent(inout) :: expr
    type(var_list_t), intent(in), target :: var_list
    call expr%var_list%link (var_list)
  end subroutine subevt_expr_link_var_list

  module subroutine subevt_expr_setup_selection (expr, ef_cuts)
    class(subevt_expr_t), intent(inout), target :: expr
    class(expr_factory_t), intent(in) :: ef_cuts
    call ef_cuts%build (expr%selection)
    if (allocated (expr%selection)) then
       call expr%setup_var_self ()
       call expr%selection%setup_lexpr (expr%var_list)
       expr%has_selection = .true.
    end if
  end subroutine subevt_expr_setup_selection

  module subroutine subevt_expr_colorize (expr, colorize_subevt)
    class(subevt_expr_t), intent(inout), target :: expr
    logical, intent(in) :: colorize_subevt
    expr%colorize_subevt = colorize_subevt
  end subroutine subevt_expr_colorize

  module subroutine subevt_expr_reset_contents (expr)
    class(subevt_expr_t), intent(inout) :: expr
    expr%subevt_filled = .false.
  end subroutine subevt_expr_reset_contents

  module subroutine subevt_expr_evaluate (expr, passed)
    class(subevt_expr_t), intent(inout) :: expr
    logical, intent(out) :: passed
    if (expr%has_selection) then
       call expr%selection%evaluate ()
       if (expr%selection%is_known ()) then
          passed = expr%selection%get_log ()
       else
          call msg_error ("Evaluate selection expression: result undefined")
          passed = .false.
       end if
    else
       passed = .true.
    end if
  end subroutine subevt_expr_evaluate

  module subroutine parton_expr_final (object)
    class(parton_expr_t), intent(inout) :: object
    call object%base_final ()
    if (object%has_scale) then
       call object%scale%final ()
    end if
    if (object%has_fac_scale) then
       call object%fac_scale%final ()
    end if
    if (object%has_ren_scale) then
       call object%ren_scale%final ()
    end if
    if (object%has_weight) then
       call object%weight%final ()
    end if
  end subroutine parton_expr_final

  module subroutine parton_expr_write (object, unit, prefix, pacified)
    class(parton_expr_t), intent(in) :: object
    integer, intent(in), optional :: unit
    character(*), intent(in), optional :: prefix
    logical, intent(in), optional :: pacified
    integer :: u
    u = given_output_unit (unit)
    call object%base_write (u, pacified = pacified)
    if (object%subevt_filled) then
       if (object%has_scale) then
          call write_separator (u)
          write (u, "(1x,A)")  "Scale expression:"
          call write_separator (u)
          call object%scale%write (u)
       end if
       if (object%has_fac_scale) then
          call write_separator (u)
          write (u, "(1x,A)")  "Factorization scale expression:"
          call write_separator (u)
          call object%fac_scale%write (u)
       end if
       if (object%has_ren_scale) then
          call write_separator (u)
          write (u, "(1x,A)")  "Renormalization scale expression:"
          call write_separator (u)
          call object%ren_scale%write (u)
       end if
       if (object%has_weight) then
          call write_separator (u)
          write (u, "(1x,A)")  "Weight expression:"
          call write_separator (u)
          call object%weight%write (u)
       end if
    end if
  end subroutine parton_expr_write

  module subroutine parton_expr_setup_vars (expr, sqrts)
    class(parton_expr_t), intent(inout), target :: expr
    real(default), intent(in) :: sqrts
    call expr%base_setup_vars (sqrts)
  end subroutine parton_expr_setup_vars

  module subroutine parton_expr_setup_scale (expr, ef_scale)
    class(parton_expr_t), intent(inout), target :: expr
    class(expr_factory_t), intent(in) :: ef_scale
    call ef_scale%build (expr%scale)
    if (allocated (expr%scale)) then
       call expr%setup_var_self ()
       call expr%scale%setup_expr (expr%var_list)
       expr%has_scale = .true.
    end if
  end subroutine parton_expr_setup_scale

  module subroutine parton_expr_setup_fac_scale (expr, ef_fac_scale)
    class(parton_expr_t), intent(inout), target :: expr
    class(expr_factory_t), intent(in) :: ef_fac_scale
    call ef_fac_scale%build (expr%fac_scale)
    if (allocated (expr%fac_scale)) then
       call expr%setup_var_self ()
       call expr%fac_scale%setup_expr (expr%var_list)
       expr%has_fac_scale = .true.
    end if
  end subroutine parton_expr_setup_fac_scale

  module subroutine parton_expr_setup_ren_scale (expr, ef_ren_scale)
    class(parton_expr_t), intent(inout), target :: expr
    class(expr_factory_t), intent(in) :: ef_ren_scale
    call ef_ren_scale%build (expr%ren_scale)
    if (allocated (expr%ren_scale)) then
       call expr%setup_var_self ()
       call expr%ren_scale%setup_expr (expr%var_list)
       expr%has_ren_scale = .true.
    end if
  end subroutine parton_expr_setup_ren_scale

  module subroutine parton_expr_setup_weight (expr, ef_weight)
    class(parton_expr_t), intent(inout), target :: expr
    class(expr_factory_t), intent(in) :: ef_weight
    call ef_weight%build (expr%weight)
    if (allocated (expr%weight)) then
       call expr%setup_var_self ()
       call expr%weight%setup_expr (expr%var_list)
       expr%has_weight = .true.
    end if
  end subroutine parton_expr_setup_weight

  module subroutine parton_expr_setup_subevt (expr, int, &
       i_beam, i_in, i_out, f_beam, f_in, f_out)
    class(parton_expr_t), intent(inout) :: expr
    type(interaction_t), intent(in), target :: int
    integer, dimension(:), intent(in) :: i_beam, i_in, i_out
    type(flavor_t), dimension(:), intent(in) :: f_beam, f_in, f_out
    allocate (expr%i_beam (size (i_beam)))
    allocate (expr%i_in (size (i_in)))
    allocate (expr%i_out (size (i_out)))
    expr%i_beam = i_beam
    expr%i_in = i_in
    expr%i_out = i_out
    call interaction_to_subevt (int, &
         expr%i_beam, expr%i_in, expr%i_out, expr%subevt_t)
    call expr%set_pdg_beam (f_beam%get_pdg ())
    call expr%set_pdg_incoming (f_in%get_pdg ())
    call expr%set_pdg_outgoing (f_out%get_pdg ())
    call expr%set_p2_beam     (f_beam%get_mass () ** 2)
    call expr%set_p2_incoming (f_in%get_mass ()   ** 2)
    call expr%set_p2_outgoing (f_out%get_mass ()  ** 2)
    expr%n_in  = size (i_in)
    expr%n_out = size (i_out)
    expr%n_tot = expr%n_in + expr%n_out
  end subroutine parton_expr_setup_subevt

  module subroutine parton_expr_renew_flv_content_subevt (expr, int, &
       i_beam, i_in, i_out, f_beam, f_in, f_out)
    class(parton_expr_t), intent(inout) :: expr
    type(interaction_t), intent(in), target :: int
    integer, dimension(:), intent(in) :: i_beam, i_in, i_out
    type(flavor_t), dimension(:), intent(in) :: f_beam, f_in, f_out
    expr%i_beam = i_beam
    expr%i_in = i_in
    expr%i_out = i_out
    call expr%set_pdg_beam (f_beam%get_pdg ())
    call expr%set_pdg_incoming (f_in%get_pdg ())
    call expr%set_pdg_outgoing (f_out%get_pdg ())
    expr%n_in  = size (i_in)
    expr%n_out = size (i_out)
    expr%n_tot = expr%n_in + expr%n_out
  end subroutine parton_expr_renew_flv_content_subevt

  subroutine interaction_to_subevt (int, j_beam, j_in, j_out, subevt)
    type(interaction_t), intent(in), target :: int
    integer, dimension(:), intent(in) :: j_beam, j_in, j_out
    type(subevt_t), intent(out) :: subevt
    type(flavor_t), dimension(:), allocatable :: flv
    integer :: n_beam, n_in, n_out, i, j
    allocate (flv (int%get_n_tot ()))
    flv = quantum_numbers_get_flavor (int%get_quantum_numbers (1))
    n_beam = size (j_beam)
    n_in = size (j_in)
    n_out = size (j_out)
    call subevt_init (subevt, n_beam + n_in + n_out)
    do i = 1, n_beam
       j = j_beam(i)
       call subevt%set_beam (i, flv(j)%get_pdg (), &
            vector4_null, flv(j)%get_mass () ** 2)
    end do
    do i = 1, n_in
       j = j_in(i)
       call subevt%set_incoming (n_beam + i, flv(j)%get_pdg (), &
            vector4_null, flv(j)%get_mass () ** 2)
    end do
    do i = 1, n_out
       j = j_out(i)
       call subevt%set_outgoing (n_beam + n_in + i, &
            flv(j)%get_pdg (), vector4_null, &
            flv(j)%get_mass () ** 2)
    end do
  end subroutine interaction_to_subevt

  module subroutine interaction_momenta_to_subevt_id &
       (int, j_beam, j_in, j_out, subevt)
    type(interaction_t), intent(in) :: int
    integer, dimension(:), intent(in) :: j_beam, j_in, j_out
    type(subevt_t), intent(inout) :: subevt
    call subevt%set_p_beam (- int%get_momenta (j_beam))
    call subevt%set_p_incoming (- int%get_momenta (j_in))
    call subevt%set_p_outgoing (int%get_momenta (j_out))
  end subroutine interaction_momenta_to_subevt_id

  module subroutine interaction_momenta_to_subevt_tr &
       (int, j_beam, j_in, j_out, lt, subevt)
    type(interaction_t), intent(in) :: int
    integer, dimension(:), intent(in) :: j_beam, j_in, j_out
    type(subevt_t), intent(inout) :: subevt
    type(lorentz_transformation_t), intent(in) :: lt
    call subevt%set_p_beam (- lt * int%get_momenta (j_beam))
    call subevt%set_p_incoming (- lt * int%get_momenta (j_in))
    call subevt%set_p_outgoing (lt * int%get_momenta (j_out))
  end subroutine interaction_momenta_to_subevt_tr

  module subroutine parton_expr_fill_subevt (expr, int)
    class(parton_expr_t), intent(inout) :: expr
    type(interaction_t), intent(in), target :: int
    call interaction_momenta_to_subevt (int, &
         expr%i_beam, expr%i_in, expr%i_out, expr%subevt_t)
    expr%sqrts_hat = expr%get_sqrts_hat ()
    expr%subevt_filled = .true.
  end subroutine parton_expr_fill_subevt

  module subroutine parton_expr_evaluate (expr, passed, scale, fac_scale, &
       ren_scale, weight, scale_forced, force_evaluation)
    class(parton_expr_t), intent(inout) :: expr
    logical, intent(out) :: passed
    real(default), intent(out) :: scale
    real(default), allocatable, intent(out) :: fac_scale
    real(default), allocatable, intent(out) :: ren_scale
    real(default), intent(out) :: weight
    real(default), intent(in), allocatable, optional :: scale_forced
    logical, intent(in), optional :: force_evaluation
    logical :: force_scale, force_eval
    force_scale = .false.; force_eval = .false.
    if (present (scale_forced))  force_scale = allocated (scale_forced)
    if (present (force_evaluation)) force_eval = force_evaluation
    call expr%base_evaluate (passed)
    if (passed .or. force_eval) then
       if (force_scale) then
          scale = scale_forced
       else if (expr%has_scale) then
          call expr%scale%evaluate ()
          if (expr%scale%is_known ()) then
             scale = expr%scale%get_real ()
          else
             call msg_error ("Evaluate scale expression: result undefined")
             scale = zero
          end if
       else
          scale = expr%sqrts_hat
       end if
       if (expr%has_fac_scale) then
          call expr%fac_scale%evaluate ()
          if (expr%fac_scale%is_known ()) then
             if (.not. allocated (fac_scale)) then
                allocate (fac_scale, source = expr%fac_scale%get_real ())
             else
                fac_scale = expr%fac_scale%get_real ()
             end if
          else
             call msg_error ("Evaluate factorization scale expression: &
                  &result undefined")
          end if
       end if
       if (expr%has_ren_scale) then
          call expr%ren_scale%evaluate ()
          if (expr%ren_scale%is_known ()) then
             if (.not. allocated (ren_scale)) then
                allocate (ren_scale, source = expr%ren_scale%get_real ())
             else
                ren_scale = expr%ren_scale%get_real ()
             end if
          else
             call msg_error ("Evaluate renormalization scale expression: &
                  &result undefined")
          end if
       end if
       if (expr%has_weight) then
          call expr%weight%evaluate ()
          if (expr%weight%is_known ()) then
             weight = expr%weight%get_real ()
          else
             call msg_error ("Evaluate weight expression: result undefined")
             weight = zero
          end if
       else
          weight = one
       end if
    else
       weight = zero
    end if
  end subroutine parton_expr_evaluate

  module subroutine parton_expr_get_beam_index (expr, i_beam)
    class(parton_expr_t), intent(in) :: expr
    integer, dimension(:), intent(out) :: i_beam
    i_beam = expr%i_beam
  end subroutine parton_expr_get_beam_index

  module subroutine parton_expr_get_in_index (expr, i_in)
    class(parton_expr_t), intent(in) :: expr
    integer, dimension(:), intent(out) :: i_in
    i_in = expr%i_in
  end subroutine parton_expr_get_in_index

  module subroutine event_expr_final (object)
    class(event_expr_t), intent(inout) :: object
    call object%base_final ()
    if (object%has_reweight) then
       call object%reweight%final ()
    end if
    if (object%has_analysis) then
       call object%analysis%final ()
    end if
  end subroutine event_expr_final

  module subroutine event_expr_write (object, unit, prefix, pacified)
    class(event_expr_t), intent(in) :: object
    integer, intent(in), optional :: unit
    character(*), intent(in), optional :: prefix
    logical, intent(in), optional :: pacified
    integer :: u
    u = given_output_unit (unit)
    call object%base_write (u, pacified = pacified)
    if (object%subevt_filled) then
       if (object%has_reweight) then
          call write_separator (u)
          write (u, "(1x,A)")  "Reweighting expression:"
          call write_separator (u)
          call object%reweight%write (u)
       end if
       if (object%has_analysis) then
          call write_separator (u)
          write (u, "(1x,A)")  "Analysis expression:"
          call write_separator (u)
          call object%analysis%write (u)
       end if
    end if
  end subroutine event_expr_write

  module subroutine event_expr_init (expr, n_alt)
    class(event_expr_t), intent(out) :: expr
    integer, intent(in), optional :: n_alt
    if (present (n_alt)) then
       expr%n_alt = n_alt
       allocate (expr%sqme_alt (n_alt), source = 0._default)
       allocate (expr%weight_alt (n_alt), source = 0._default)
    end if
  end subroutine event_expr_init

  module subroutine event_expr_setup_vars (expr, sqrts)
    class(event_expr_t), intent(inout), target :: expr
    real(default), intent(in) :: sqrts
    call expr%base_setup_vars (sqrts)
    call expr%var_list%append_string_ptr (var_str ("$process_id"), &
         expr%id, is_known = expr%has_id, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_int_ptr (var_str ("process_num_id"), &
         expr%num_id, is_known = expr%has_num_id, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_real_ptr (var_str ("sqme"), &
         expr%sqme_prc, is_known = expr%has_sqme_prc, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_real_ptr (var_str ("sqme_ref"), &
         expr%sqme_ref, is_known = expr%has_sqme_ref, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_int_ptr (var_str ("event_index"), &
         expr%index, is_known = expr%has_index, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_real_ptr (var_str ("event_weight"), &
         expr%weight_prc, is_known = expr%has_weight_prc, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_real_ptr (var_str ("event_weight_ref"), &
         expr%weight_ref, is_known = expr%has_weight_ref, &
         locked = .true., verbose = .false., intrinsic = .true.)
    call expr%var_list%append_real_ptr (var_str ("event_excess"), &
         expr%excess_prc, is_known = expr%has_excess_prc, &
         locked = .true., verbose = .false., intrinsic = .true.)
  end subroutine event_expr_setup_vars

  module subroutine event_expr_setup_analysis (expr, ef_analysis)
    class(event_expr_t), intent(inout), target :: expr
    class(expr_factory_t), intent(in) :: ef_analysis
    call ef_analysis%build (expr%analysis)
    if (allocated (expr%analysis)) then
       call expr%setup_var_self ()
       call expr%analysis%setup_lexpr (expr%var_list)
       expr%has_analysis = .true.
    end if
  end subroutine event_expr_setup_analysis

  module subroutine event_expr_setup_reweight (expr, ef_reweight)
    class(event_expr_t), intent(inout), target :: expr
    class(expr_factory_t), intent(in) :: ef_reweight
    call ef_reweight%build (expr%reweight)
    if (allocated (expr%reweight)) then
       call expr%setup_var_self ()
       call expr%reweight%setup_expr (expr%var_list)
       expr%has_reweight = .true.
    end if
  end subroutine event_expr_setup_reweight

  module subroutine event_expr_set_process_id (expr, id)
    class(event_expr_t), intent(inout) :: expr
    type(string_t), intent(in) :: id
    expr%id = id
    expr%has_id = .true.
  end subroutine event_expr_set_process_id

  module subroutine event_expr_set_process_num_id (expr, num_id)
    class(event_expr_t), intent(inout) :: expr
    integer, intent(in) :: num_id
    expr%num_id = num_id
    expr%has_num_id = .true.
  end subroutine event_expr_set_process_num_id

  module subroutine event_expr_reset_contents (expr)
    class(event_expr_t), intent(inout) :: expr
    call expr%base_reset_contents ()
    expr%has_sqme_ref = .false.
    expr%has_sqme_prc = .false.
    expr%has_sqme_alt = .false.
    expr%has_weight_ref = .false.
    expr%has_weight_prc = .false.
    expr%has_weight_alt = .false.
    expr%has_excess_prc = .false.
  end subroutine event_expr_reset_contents

  module subroutine event_expr_set (expr, &
       weight_ref, weight_prc, weight_alt, &
       excess_prc, &
       sqme_ref, sqme_prc, sqme_alt)
    class(event_expr_t), intent(inout) :: expr
    real(default), intent(in), optional :: weight_ref, weight_prc
    real(default), intent(in), optional :: excess_prc
    real(default), intent(in), optional :: sqme_ref, sqme_prc
    real(default), dimension(:), intent(in), optional :: sqme_alt, weight_alt
    if (present (sqme_ref)) then
       expr%has_sqme_ref = .true.
       expr%sqme_ref = sqme_ref
    end if
    if (present (sqme_prc)) then
       expr%has_sqme_prc = .true.
       expr%sqme_prc = sqme_prc
    end if
    if (present (sqme_alt)) then
       expr%has_sqme_alt = .true.
       expr%sqme_alt = sqme_alt
    end if
    if (present (weight_ref)) then
       expr%has_weight_ref = .true.
       expr%weight_ref = weight_ref
    end if
    if (present (weight_prc)) then
       expr%has_weight_prc = .true.
       expr%weight_prc = weight_prc
    end if
    if (present (weight_alt)) then
       expr%has_weight_alt = .true.
       expr%weight_alt = weight_alt
    end if
    if (present (excess_prc)) then
       expr%has_excess_prc = .true.
       expr%excess_prc = excess_prc
    end if
  end subroutine event_expr_set

  module function event_expr_has_event_index (expr) result (flag)
    class(event_expr_t), intent(in) :: expr
    logical :: flag
    flag = expr%has_index
  end function event_expr_has_event_index

  module function event_expr_get_event_index (expr) result (index)
    class(event_expr_t), intent(in) :: expr
    integer :: index
    if (expr%has_index) then
       index = expr%index
    else
       index = 0
    end if
  end function event_expr_get_event_index

  module subroutine event_expr_set_event_index (expr, index)
    class(event_expr_t), intent(inout) :: expr
    integer, intent(in) :: index
    expr%index = index
    expr%has_index = .true.
  end subroutine event_expr_set_event_index

  module subroutine event_expr_reset_event_index (expr)
    class(event_expr_t), intent(inout) :: expr
    expr%has_index = .false.
  end subroutine event_expr_reset_event_index

  module subroutine event_expr_increment_event_index (expr, offset)
    class(event_expr_t), intent(inout) :: expr
    integer, intent(in), optional :: offset
    if (expr%has_index) then
       expr%index = expr%index + 1
    else if (present (offset)) then
       call expr%set_event_index (offset + 1)
    else
       call expr%set_event_index (1)
    end if
  end subroutine event_expr_increment_event_index

  module subroutine event_expr_fill_subevt (expr, particle_set)
    class(event_expr_t), intent(inout) :: expr
    type(particle_set_t), intent(in) :: particle_set
    call particle_set%to_subevt (expr%subevt_t, expr%colorize_subevt)
    expr%sqrts_hat = expr%get_sqrts_hat ()
    expr%n_in  = expr%get_n_in  ()
    expr%n_out = expr%get_n_out ()
    expr%n_tot = expr%n_in + expr%n_out
    expr%subevt_filled = .true.
  end subroutine event_expr_fill_subevt

  module subroutine event_expr_evaluate (expr, passed, reweight, analysis_flag)
    class(event_expr_t), intent(inout) :: expr
    logical, intent(out) :: passed
    real(default), intent(out) :: reweight
    logical, intent(out) :: analysis_flag
    call expr%base_evaluate (passed)
    if (passed) then
       if (expr%has_reweight) then
          call expr%reweight%evaluate ()
          if (expr%reweight%is_known ()) then
             reweight = expr%reweight%get_real ()
          else
             call msg_error ("Evaluate reweight expression: &
                  &result undefined")
             reweight = 0
          end if
       else
          reweight = 1
       end if
       if (expr%has_analysis) then
          call expr%analysis%evaluate ()
          if (expr%analysis%is_known ()) then
             analysis_flag = expr%analysis%get_log ()
          else
             call msg_error ("Evaluate analysis expression: &
                  &result undefined")
             analysis_flag = .false.
          end if
       else
          analysis_flag = .true.
       end if
    end if
  end subroutine event_expr_evaluate


end submodule subevt_expr_s

