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

module api

  use iso_fortran_env, only: int32, real64 !NODEP!
  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use diagnostics, only: logfile_init
  use diagnostics, only: logging
  use diagnostics, only: logfile_final
  use diagnostics, only: mask_term_signals
  use diagnostics, only: release_term_signals
  use diagnostics, only: msg_message
  use diagnostics, only: msg_error
  use diagnostics, only: msg_fatal
  use diagnostics, only: msg_summary
  use string_utils, only: split_string
  use lexers, only: stream_t
  use lexers, only: lexer_t
  use parser, only: parse_tree_t
  use parser, only: parse_node_t
  use flavors, only: flavor_t
  use event_handles, only: event_handle_t
  use events, only: event_t
  use event_streams, only: event_stream_array_t
  use simulations, only: simulation_t
  use commands, only: lexer_init_cmd_list
  use commands, only: syntax_cmd_list
  use commands, only: command_list_t
  use whizard, only: whizard_t
  use whizard, only: whizard_options_t

  implicit none
  private

  public :: whizard_api_t
  public :: simulation_api_t

  type :: whizard_api_t
     private
     type(whizard_t), pointer :: master => null ()
     character(:), allocatable :: logfile
     type(whizard_options_t) :: options
     type(lexer_t) :: sindarin_lexer
   contains
     procedure :: option
     procedure :: init
     procedure :: final
   !   final :: finalize
     generic :: set_var => set_var_real64
     generic :: set_var => set_var_int32
     generic :: set_var => set_var_logical
     generic :: set_var => set_var_character
     procedure :: set_var_real64
     procedure :: set_var_int32
     procedure :: set_var_logical
     procedure :: set_var_character
     generic :: get_var => get_var_real64
     generic :: get_var => get_var_int32
     generic :: get_var => get_var_logical
     generic :: get_var => get_var_character
     procedure :: get_var_real64
     procedure :: get_var_int32
     procedure :: get_var_logical
     procedure :: get_var_character
     procedure :: get_var_character_length
     procedure :: get_integration_result
     procedure :: command
     generic :: flv_string => flv_string_single, flv_string_array
     procedure :: flv_string_single
     procedure :: flv_string_array
     procedure :: new_sample
  end type whizard_api_t

  type :: simulation_api_t
     private
     type(simulation_t), pointer :: sim => null ()
     type(event_stream_array_t) :: esa
     integer :: it_begin = 0
     integer :: it_end = 0
   contains
     procedure :: final => sample_final
   !   final :: sample_finalize
     procedure :: open => sample_open
     procedure :: next_event => sample_next_event
     procedure :: close => sample_close
     procedure :: get_event_index => sample_get_event_index
     procedure :: get_process_index => sample_get_process_index
     procedure :: get_process_id => sample_get_process_id
     procedure :: get_sqrts => sample_get_sqrts
     procedure :: get_fac_scale => sample_get_fac_scale
     procedure :: get_alpha_s => sample_get_alpha_s
     procedure :: get_sqme => sample_get_sqme
     procedure :: get_weight => sample_get_weight
  end type simulation_api_t




  logical :: whizard_banner_shown = .false.

contains

  subroutine option (whizard, key, value)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: key
    character(*), intent(in) :: value

    logical :: rebuild

    if (associated (whizard%master)) then
       call msg_error ("WHIZARD: options must be set before initialization; &
            &extra option '" // key // "' ignored.")
       return
    end if

    associate (options => whizard%options)
      select case (key)
      case ("model")
         options%preload_model = value
      case ("library")
         options%preload_libraries = value
      case ("logfile")
         whizard%logfile = value
      case ("job_id")
         options%job_id = value
      case ("pack")
         call split_string (var_str (value), var_str (","), &
              options%pack_args)
         options%pack_args = trim (adjustl (options%pack_args))
      case ("unpack")
         call split_string (var_str (value), var_str (","), &
              options%unpack_args)
         options%unpack_args = trim (adjustl (options%unpack_args))
      case ("rebuild")
         read (value, "(L1)")  rebuild
         options%rebuild_library = rebuild
         options%recompile_library = rebuild
         options%rebuild_phs = rebuild
         options%rebuild_grids = rebuild
         options%rebuild_events = rebuild
      case ("rebuild_library")
         read (value, "(L1)")  options%rebuild_library
      case ("recompile")
         read (value, "(L1)")  options%recompile_library
      case ("rebuild_phase_space")
         read (value, "(L1)")  options%rebuild_phs
      case ("rebuild_grids")
         read (value, "(L1)")  options%rebuild_grids
      case ("rebuild_events")
         read (value, "(L1)")  options%rebuild_events
      case default
         call msg_error ("WHIZARD: option '" // key // "' not recognized.")
      end select
    end associate

  end subroutine option

  subroutine init (whizard)
    class(whizard_api_t), intent(inout) :: whizard

    whizard%options%default_lib = "whizard_processes"
    if (.not. associated (whizard%master)) then
       allocate (whizard%master)
       if (allocated (whizard%logfile)) then
          logging = .true.
          call logfile_init (var_str (whizard%logfile))
          call whizard%master%init (whizard%options, &
               logfile = var_str (whizard%logfile))
       else
          call whizard%master%init (whizard%options)
       end if
       call lexer_init_cmd_list (whizard%sindarin_lexer)
       call msg_message ("WHIZARD: master object initialized.")
    else
       call msg_error ("WHIZARD: extra call to initializer is ignored")
    end if

  end subroutine init

  subroutine final (whizard)
    class(whizard_api_t), intent(inout) :: whizard

    if (associated (whizard%master)) then
       call whizard%sindarin_lexer%final ()
       call whizard%master%final ()
       call msg_message ("WHIZARD: master object finalized.")
       deallocate (whizard%master)
       call msg_summary ()
       call logfile_final ()
    else
       call msg_error ("WHIZARD: extra call to finalizer is ignored")
    end if

  end subroutine final

  subroutine finalize (whizard)
    type(whizard_api_t), intent(inout), target :: whizard

    if (associated (whizard%master))  call whizard%final ()

  end subroutine finalize

  subroutine sanity_check (whizard)
    class(whizard_api_t), intent(in) :: whizard

    if (.not. associated (whizard%master)) then
       call msg_fatal ("WHIZARD: method call without initialization")
    end if

  end subroutine sanity_check

  subroutine set_var_real64 (whizard, name, value)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: name
    real(real64), intent(in) :: value

    call sanity_check (whizard)
    call whizard%master%global%set_real &
         (var_str (name), real (value, default), &
         verbose=.true., is_known=.true.)

  end subroutine set_var_real64

  subroutine set_var_int32 (whizard, name, value)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: name
    integer(int32), intent(in) :: value

    call sanity_check (whizard)
    call whizard%master%global%set_int &
         (var_str (name), int (value), &
         verbose=.true., is_known=.true.)

  end subroutine set_var_int32

  subroutine set_var_logical (whizard, name, value)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: name
    logical, intent(in) :: value

    call sanity_check (whizard)
    call whizard%master%global%set_log &
         (var_str (name), value, &
         verbose=.true., is_known=.true.)

  end subroutine set_var_logical

  subroutine set_var_character (whizard, name, value)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: name
    character(*), intent(in) :: value

    call sanity_check (whizard)
    call whizard%master%global%set_string &
         (var_str (name), var_str (value), &
         verbose=.true., is_known=.true.)

  end subroutine set_var_character

  subroutine get_var_real64 (whizard, name, value, known)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: name
    real(real64), intent(out) :: value
    logical, intent(out), optional :: known

    call sanity_check (whizard)
    value = whizard%master%global%get_rval (var_str (name))
    if (present (known)) &
         known = whizard%master%global%is_known (var_str (name))

  end subroutine get_var_real64

  subroutine get_var_int32 (whizard, name, value, known)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: name
    integer(int32), intent(out) :: value
    logical, intent(out), optional :: known

    call sanity_check (whizard)
    value = whizard%master%global%get_ival (var_str (name))
    if (present (known)) &
         known = whizard%master%global%is_known (var_str (name))

  end subroutine get_var_int32

  subroutine get_var_logical (whizard, name, value, known)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: name
    logical, intent(out) :: value
    logical, intent(out), optional :: known

    call sanity_check (whizard)
    value = whizard%master%global%get_lval (var_str (name))
    if (present (known)) &
         known = whizard%master%global%is_known (var_str (name))

  end subroutine get_var_logical

  subroutine get_var_character (whizard, name, value, known, strlen)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: name
    character(:), allocatable, intent(out) :: value
    logical, intent(out), optional :: known
    integer, intent(in), optional :: strlen

    call sanity_check (whizard)
    value = char (whizard%master%global%get_sval (var_str (name)))
    if (present (known)) &
         known = whizard%master%global%is_known (var_str (name))
    if (present (strlen)) then
       if (len (value) > strlen) value = value(1:strlen)
    end if

  end subroutine get_var_character

  function get_var_character_length (whizard, name) result (strlen)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: name
    integer :: strlen

    call sanity_check (whizard)
    strlen = len (whizard%master%global%get_sval (var_str (name)))

  end function get_var_character_length

  subroutine get_integration_result (whizard, proc_id, integral, error, known)
    class(whizard_api_t), intent(in) :: whizard
    character(*), intent(in) :: proc_id
    real(real64), intent(out) :: integral
    real(real64), intent(out) :: error
    logical, intent(out), optional :: known

    character(:), allocatable :: integral_var
    character(:), allocatable :: error_var

    integral_var = "integral(" // proc_id // ")"
    error_var = "error(" // proc_id // ")"

    integral = whizard%master%global%get_rval (var_str (integral_var))
    error = whizard%master%global%get_rval (var_str (error_var))
    if (present (known)) &
         known = whizard%master%global%is_known (var_str (integral_var))

  end subroutine get_integration_result

  subroutine command (whizard, code)
    class(whizard_api_t), intent(inout) :: whizard
    character(*), intent(in) :: code

    type(stream_t), target :: stream
    type(parse_tree_t) :: parse_tree
    type(parse_node_t), pointer :: pn_root
    type(command_list_t), target :: cmd_list

    call sanity_check (whizard)

    call stream%init (var_str (code))
    call whizard%sindarin_lexer%assign_stream (stream)
    call parse_tree%parse (syntax_cmd_list, whizard%sindarin_lexer)
    pn_root => parse_tree%get_root_ptr ()
    if (associated (pn_root)) then
       call cmd_list%compile (pn_root, whizard%master%global)
       call mask_term_signals ()
       call cmd_list%execute (whizard%master%global)
       call release_term_signals ()
    end if
    call stream%final ()
    call whizard%sindarin_lexer%clear ()

  end subroutine command

  function flv_string_single (whizard, pdg) result (string)
    class(whizard_api_t), intent(in) :: whizard
    integer(int32), intent(in) :: pdg
    character(:), allocatable :: string

    type(flavor_t) :: flv

    if (associated (whizard%master%global%model)) then
       call flv%init (int (pdg), whizard%master%global%model)
    end if
    string = '"' // char (flv%get_name ()) // '"'

  end function flv_string_single

  function flv_string_array (whizard, pdg) result (string)
    class(whizard_api_t), intent(in) :: whizard
    integer(int32), dimension(:), intent(in) :: pdg
    character(:), allocatable :: string

    integer :: i

    if (size (pdg) > 0) then
       string = whizard%flv_string (pdg(1))
       do i = 2, size (pdg)
          string = string // ":" // whizard%flv_string (pdg(i))
       end do
    else
       string = ""
    end if

  end function flv_string_array

  subroutine new_sample (whizard, proc_id_string, sample)
    class(whizard_api_t), intent(in) :: whizard
    character(*), intent(in) :: proc_id_string
    type(simulation_api_t), intent(out) :: sample

    integer :: i, n_events
    type(string_t), dimension(:), allocatable :: proc_id

    call sanity_check (whizard)

    call split_string (var_str (proc_id_string), var_str (","), proc_id)
    proc_id = trim (adjustl (proc_id))

    allocate (sample%sim)
    call sample%sim%init (proc_id, .true., .true., whizard%master%global)
    if (sample%sim%is_valid ()) then

       call sample%sim%init_process_selector ()
       call sample%sim%setup_openmp ()
       call sample%sim%compute_n_events (n_events)
       call sample%sim%set_n_events_requested (n_events)
       call sample%sim%activate_extra_logging ()
       call sample%sim%prepare_event_streams (sample%esa)

    end if

  end subroutine new_sample

  subroutine sample_final (sample)
    class(simulation_api_t), intent(inout) :: sample

    call sample%esa%final ()

    if (associated (sample%sim)) then
       call sample%sim%final ()
       deallocate (sample%sim)
    end if

  end subroutine sample_final

  subroutine sample_finalize (sample)
    type(simulation_api_t), intent(inout) :: sample

    call sample%final ()

  end subroutine sample_finalize

  subroutine sample_open (sample, it_begin, it_end)
    class(simulation_api_t), intent(inout) :: sample
    integer, intent(out) :: it_begin
    integer, intent(out) :: it_end

    if (sample%sim%is_valid ()) then
       if (sample%esa%is_valid ()) then
          call sample%sim%before_first_event (it_begin, it_end, sample%esa)
       else
          call sample%sim%before_first_event (it_begin, it_end)
       end if
    end if
    sample%it_begin = it_begin
    sample%it_end = it_end

  end subroutine sample_open

  subroutine sample_next_event (sample, event_handle_out, event_handle_in)
    class(simulation_api_t), intent(inout) :: sample
    class(event_handle_t), intent(inout), optional :: event_handle_out
    class(event_handle_t), intent(inout), optional :: event_handle_in

    if (sample%sim%is_valid ()) then
       call mask_term_signals ()
       if (sample%esa%is_valid ()) then
          call sample%sim%next_event &
               (sample%esa, event_handle_out, event_handle_in)
       else
          call sample%sim%next_event ()
       end if
       call release_term_signals ()
    end if

  end subroutine sample_next_event

  subroutine sample_close (sample)
    class(simulation_api_t), intent(inout) :: sample

    if (sample%sim%is_valid ()) then
       call sample%sim%after_last_event (sample%it_begin, sample%it_end)
    end if
    call sample%final ()

  end subroutine sample_close

  subroutine sample_get_event_index (sample, i)
    class(simulation_api_t), intent(in) :: sample
    integer(int32), intent(out) :: i

    i = sample%sim%get_event_index ()

  end subroutine sample_get_event_index

  subroutine sample_get_process_index (sample, i)
    class(simulation_api_t), intent(in) :: sample
    integer(int32), intent(out) :: i

    i = sample%sim%get_process_index ()

  end subroutine sample_get_process_index

  subroutine sample_get_process_id (sample, proc_id)
    class(simulation_api_t), intent(in) :: sample
    character(:), allocatable, intent(out) :: proc_id

    class(event_t), pointer :: event

    event => sample%sim%get_event_ptr ()
    proc_id = char (event%get_process_name ())

  end subroutine sample_get_process_id

  subroutine sample_get_sqrts (sample, sqrts)
    class(simulation_api_t), intent(in) :: sample
    real(real64), intent(out) :: sqrts

    class(event_t), pointer :: event

    event => sample%sim%get_event_ptr ()
    sqrts = event%get_sqrts ()

  end subroutine sample_get_sqrts

  subroutine sample_get_fac_scale (sample, fac_scale)
    class(simulation_api_t), intent(in) :: sample
    real(real64), intent(out) :: fac_scale

    class(event_t), pointer :: event

    event => sample%sim%get_event_ptr ()
    fac_scale = event%get_fac_scale ()

  end subroutine sample_get_fac_scale

  subroutine sample_get_alpha_s (sample, alpha_s)
    class(simulation_api_t), intent(in) :: sample
    real(real64), intent(out) :: alpha_s

    class(event_t), pointer :: event

    event => sample%sim%get_event_ptr ()
    alpha_s = event%get_alpha_s ()

  end subroutine sample_get_alpha_s

  subroutine sample_get_sqme (sample, sqme)
    class(simulation_api_t), intent(in) :: sample
    real(real64), intent(out) :: sqme

    class(event_t), pointer :: event

    event => sample%sim%get_event_ptr ()
    sqme = event%get_sqme_prc ()

  end subroutine sample_get_sqme

  subroutine sample_get_weight (sample, weight)
    class(simulation_api_t), intent(in) :: sample
    real(real64), intent(out) :: weight

    class(event_t), pointer :: event

    event => sample%sim%get_event_ptr ()
    weight = event%get_weight_prc ()

  end subroutine sample_get_weight


end module api
