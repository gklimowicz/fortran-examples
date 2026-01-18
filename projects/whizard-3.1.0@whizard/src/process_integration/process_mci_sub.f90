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

submodule (process_mci) process_mci_s

  use debug_master, only: debug_on
  use io_units
  use diagnostics
  use physics_defs
  use md5

  implicit none

contains

  module subroutine process_mci_entry_final (object)
    class(process_mci_entry_t), intent(inout) :: object
    if (allocated (object%mci))  call object%mci%final ()
  end subroutine process_mci_entry_final

  module subroutine process_mci_entry_write (object, unit, pacify)
    class(process_mci_entry_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: pacify
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A,I0)")  "Associated components = ", object%i_component
    write (u, "(3x,A,I0)")  "MC input parameters   = ", object%n_par
    write (u, "(3x,A,I0)")  "MC parameters (SF)    = ", object%n_par_sf
    write (u, "(3x,A,I0)")  "MC parameters (PHS)   = ", object%n_par_phs
    if (object%pass > 0) then
       write (u, "(3x,A,I0)")  "Current pass          = ", object%pass
       write (u, "(3x,A,I0)")  "Number of iterations  = ", object%n_it
       write (u, "(3x,A,I0)")  "Number of calls       = ", object%n_calls
    end if
    if (object%md5sum /= "") then
       write (u, "(3x,A,A,A)") "MD5 sum (components)  = '", object%md5sum, "'"
    end if
    if (allocated (object%mci)) then
       call object%mci%write (u)
    end if
    call object%counter%write (u)
    if (object%results%exist ()) then
       call object%results%write (u, suppress = pacify)
       call object%results%write_chain_weights (u)
    end if
  end subroutine process_mci_entry_write

  module subroutine process_mci_entry_configure (mci_entry, mci_template, &
       process_type, i_mci, i_component, component, &
       n_sfpar, rng_factory)
    class(process_mci_entry_t), intent(inout) :: mci_entry
    class(mci_t), intent(in), allocatable :: mci_template
    integer, intent(in) :: process_type
    integer, intent(in) :: i_mci
    integer, intent(in) :: i_component
    type(process_component_t), intent(in), target :: component
    integer, intent(in) :: n_sfpar
    class(rng_factory_t), intent(inout) :: rng_factory
    class(rng_t), allocatable :: rng
    associate (phs_config => component%phs_config)
      mci_entry%i_mci = i_mci
      call mci_entry%create_component_list (i_component, component%get_config ())
      mci_entry%n_par_sf = n_sfpar
      mci_entry%n_par_phs = phs_config%get_n_par ()
      mci_entry%n_par = mci_entry%n_par_sf + mci_entry%n_par_phs
      mci_entry%process_type = process_type
      if (allocated (mci_template)) then
         allocate (mci_entry%mci, source = mci_template)
         call mci_entry%mci%record_index (mci_entry%i_mci)
         call mci_entry%mci%set_dimensions &
              (mci_entry%n_par, phs_config%get_n_channel ())
         call mci_entry%mci%declare_flat_dimensions &
              (phs_config%get_flat_dimensions ())
         if (phs_config%provides_equivalences) then
            call mci_entry%mci%declare_equivalences &
                 (phs_config%channel, mci_entry%n_par_sf)
         end if
         if (phs_config%provides_chains) then
            call mci_entry%mci%declare_chains (phs_config%chain)
         end if
         call rng_factory%make (rng)
         call mci_entry%mci%import_rng (rng)
      end if
      call mci_entry%results%init (process_type)
    end associate
  end subroutine process_mci_entry_configure

  module subroutine process_mci_entry_create_component_list (mci_entry, &
     i_component, component_config)
    class (process_mci_entry_t), intent(inout) :: mci_entry
    integer, intent(in) :: i_component
    type(process_component_def_t), intent(in) :: component_config
    integer, dimension(:), allocatable :: i_list
    integer :: n
    integer, save :: i_rfin_offset = 0
    if (debug_on)  call msg_debug &
         (D_PROCESS_INTEGRATION, "process_mci_entry_create_component_list")
    if (mci_entry%combined_integration) then
       if (debug_on) call msg_debug (D_PROCESS_INTEGRATION, &
            "mci_entry%real_partition_type", mci_entry%real_partition_type)
       n = get_n_components (mci_entry%real_partition_type)
       allocate (i_list (n))
       select case (mci_entry%real_partition_type)
       case (REAL_FULL)
          i_list = component_config%get_association_list ()
          allocate (mci_entry%i_component (size (i_list)))
          mci_entry%i_component = i_list
       case (REAL_SINGULAR)
          i_list = component_config%get_association_list (ASSOCIATED_REAL_FIN)
          allocate (mci_entry%i_component (size(i_list)))
          mci_entry%i_component = i_list
       case (REAL_FINITE)
          allocate (mci_entry%i_component (1))
          mci_entry%i_component(1) = &
               component_config%get_associated_real_fin () + i_rfin_offset
          i_rfin_offset = i_rfin_offset + 1
       end select
    else
       allocate (mci_entry%i_component (1))
       mci_entry%i_component(1) = i_component
    end if
  contains
    function get_n_components (real_partition_type) result (n_components)
      integer :: n_components
      integer, intent(in) :: real_partition_type
      select case (real_partition_type)
      case (REAL_FULL)
         n_components = size (component_config%get_association_list ())
      case (REAL_SINGULAR)
         n_components = size (component_config%get_association_list &
            (ASSOCIATED_REAL_FIN))
      end select
      if (debug_on)  call msg_debug &
           (D_PROCESS_INTEGRATION, "n_components", n_components)
    end function get_n_components
  end subroutine process_mci_entry_create_component_list

  module subroutine process_mci_entry_set_parameters (mci_entry, var_list)
    class(process_mci_entry_t), intent(inout) :: mci_entry
    type(var_list_t), intent(in) :: var_list
    integer :: integration_results_verbosity
    real(default) :: error_threshold
    integration_results_verbosity = &
         var_list%get_ival (var_str ("integration_results_verbosity"))
    error_threshold = &
         var_list%get_rval (var_str ("error_threshold"))
    mci_entry%activate_timer = &
         var_list%get_lval (var_str ("?integration_timer"))
    call mci_entry%results%set_verbosity (integration_results_verbosity)
    call mci_entry%results%set_error_threshold (error_threshold)
  end subroutine process_mci_entry_set_parameters

  module subroutine process_mci_entry_compute_md5sum (mci_entry, &
       config, component, beam_config)
    class(process_mci_entry_t), intent(inout) :: mci_entry
    type(process_config_data_t), intent(in) :: config
    type(process_component_t), dimension(:), intent(in) :: component
    type(process_beam_config_t), intent(in) :: beam_config
    type(string_t) :: buffer
    integer :: i
    if (mci_entry%md5sum == "") then
       buffer = config%get_md5sum () // beam_config%get_md5sum ()
       do i = 1, size (component)
          if (component(i)%is_active ()) then
             buffer = buffer // component(i)%get_md5sum ()
          end if
       end do
       mci_entry%md5sum = md5sum (char (buffer))
    end if
    if (allocated (mci_entry%mci)) then
       call mci_entry%mci%set_md5sum (mci_entry%md5sum)
    end if
  end subroutine process_mci_entry_compute_md5sum

  module subroutine process_mci_entry_sampler_test &
       (mci_entry, mci_sampler, n_calls)
    class(process_mci_entry_t), intent(inout) :: mci_entry
    class(mci_sampler_t), intent(inout), target :: mci_sampler
    integer, intent(in) :: n_calls
    call mci_entry%mci%sampler_test (mci_sampler, n_calls)
  end subroutine process_mci_entry_sampler_test

  module subroutine process_mci_entry_integrate (mci_entry, mci_instance, &
         mci_sampler, n_it, n_calls, &
       adapt_grids, adapt_weights, final, pacify, &
       nlo_type)
    class(process_mci_entry_t), intent(inout) :: mci_entry
    class(mci_instance_t), intent(inout) :: mci_instance
    class(mci_sampler_t), intent(inout) :: mci_sampler
    integer, intent(in) :: n_it
    integer, intent(in) :: n_calls
    logical, intent(in), optional :: adapt_grids
    logical, intent(in), optional :: adapt_weights
    logical, intent(in), optional :: final, pacify
    integer, intent(in), optional :: nlo_type
    integer :: u_log
    u_log = logfile_unit ()
    mci_entry%pass = mci_entry%pass + 1
    mci_entry%n_it = n_it
    mci_entry%n_calls = n_calls
    if (mci_entry%pass == 1)  &
         call mci_entry%mci%startup_message (n_calls = n_calls)
    call mci_entry%mci%set_timer (active = mci_entry%activate_timer)
    call mci_entry%results%display_init (screen = .true., unit = u_log)
    call mci_entry%results%new_pass ()
    if (present (nlo_type)) then
       select case (nlo_type)
       case (NLO_VIRTUAL, NLO_REAL, NLO_MISMATCH, NLO_DGLAP)
          mci_instance%negative_weights = .true.
       end select
    end if
    call mci_entry%mci%add_pass (adapt_grids, adapt_weights, final)
    call mci_entry%mci%start_timer ()
    call mci_entry%mci%integrate (mci_instance, mci_sampler, n_it, &
         n_calls, mci_entry%results, pacify = pacify)
    call mci_entry%mci%stop_timer ()
    if (signal_is_pending ())  return
  end subroutine process_mci_entry_integrate

  module subroutine process_mci_entry_final_integration (mci_entry)
    class(process_mci_entry_t), intent(inout) :: mci_entry
    call mci_entry%results%display_final ()
    call mci_entry%time_message ()
  end subroutine process_mci_entry_final_integration

  module subroutine process_mci_entry_get_time (mci_entry, time, sample)
    class(process_mci_entry_t), intent(in) :: mci_entry
    type(time_t), intent(out) :: time
    integer, intent(in) :: sample
    real(default) :: time_last_pass, efficiency, calls
    time_last_pass = mci_entry%mci%get_time ()
    calls = mci_entry%results%get_n_calls ()
    efficiency = mci_entry%mci%get_efficiency ()
    if (time_last_pass > 0 .and. calls > 0 .and. efficiency > 0) then
       time = nint (time_last_pass / calls / efficiency * sample)
    end if
  end subroutine process_mci_entry_get_time

  module subroutine process_mci_entry_time_message (mci_entry)
    class(process_mci_entry_t), intent(in) :: mci_entry
    type(time_t) :: time
    integer :: sample
    sample = 10000
    call mci_entry%get_time (time, sample)
    if (time%is_known ()) then
       call msg_message ("Time estimate for generating 10000 events: " &
            // char (time%to_string_dhms ()))
    end if
  end subroutine process_mci_entry_time_message

  module subroutine process_mci_entry_prepare_simulation (mci_entry)
    class(process_mci_entry_t), intent(inout) :: mci_entry
    call mci_entry%mci%prepare_simulation ()
  end subroutine process_mci_entry_prepare_simulation

  module subroutine process_mci_entry_generate_weighted_event (mci_entry, &
      mci_instance, mci_sampler, keep_failed)
    class(process_mci_entry_t), intent(inout) :: mci_entry
    class(mci_instance_t), intent(inout) :: mci_instance
    class(mci_sampler_t), intent(inout) :: mci_sampler
    logical, intent(in) :: keep_failed
    logical :: generate_new
    generate_new = .true.
    call mci_instance%reset_n_event_dropped ()
    REJECTION: do while (generate_new)
       call mci_entry%mci%generate_weighted_event (mci_instance, mci_sampler)
       if (signal_is_pending ())  return
       if (.not. mci_sampler%is_valid()) then
          if (keep_failed) then
             generate_new = .false.
          else
             call mci_instance%record_event_dropped ()
             generate_new = .true.
          end if
       else
          generate_new = .false.
       end if
    end do REJECTION
  end subroutine process_mci_entry_generate_weighted_event

  module subroutine process_mci_entry_generate_unweighted_event &
       (mci_entry, mci_instance, mci_sampler)
    class(process_mci_entry_t), intent(inout) :: mci_entry
    class(mci_instance_t), intent(inout) :: mci_instance
    class(mci_sampler_t), intent(inout) :: mci_sampler
    call mci_entry%mci%generate_unweighted_event (mci_instance, mci_sampler)
  end subroutine process_mci_entry_generate_unweighted_event

  module function process_mci_entry_has_integral (mci_entry) result (flag)
    class(process_mci_entry_t), intent(in) :: mci_entry
    logical :: flag
    flag = mci_entry%results%exist ()
  end function process_mci_entry_has_integral

  module function process_mci_entry_get_integral (mci_entry) result (integral)
    class(process_mci_entry_t), intent(in) :: mci_entry
    real(default) :: integral
    integral = mci_entry%results%get_integral ()
  end function process_mci_entry_get_integral

  module function process_mci_entry_get_error (mci_entry) result (error)
    class(process_mci_entry_t), intent(in) :: mci_entry
    real(default) :: error
    error = mci_entry%results%get_error ()
  end function process_mci_entry_get_error

  module function process_mci_entry_get_accuracy (mci_entry) result (accuracy)
    class(process_mci_entry_t), intent(in) :: mci_entry
    real(default) :: accuracy
    accuracy = mci_entry%results%get_accuracy ()
  end function process_mci_entry_get_accuracy

  module function process_mci_entry_get_chi2 (mci_entry) result (chi2)
    class(process_mci_entry_t), intent(in) :: mci_entry
    real(default) :: chi2
    chi2 = mci_entry%results%get_chi2 ()
  end function process_mci_entry_get_chi2

  module function process_mci_entry_get_efficiency &
       (mci_entry) result (efficiency)
    class(process_mci_entry_t), intent(in) :: mci_entry
    real(default) :: efficiency
    efficiency = mci_entry%results%get_efficiency ()
  end function process_mci_entry_get_efficiency

  pure module function process_mci_entry_get_md5sum (entry) result (md5sum)
    class(process_mci_entry_t), intent(in) :: entry
    character(32) :: md5sum
    md5sum = entry%mci%get_md5sum ()
  end function process_mci_entry_get_md5sum

  module subroutine mci_work_write (mci_work, unit, testflag)
    class(mci_work_t), intent(in) :: mci_work
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A,I0,A)")  "Active MCI instance #", &
         mci_work%config%i_mci, " ="
    write (u, "(2x)", advance="no")
    do i = 1, mci_work%config%n_par
       write (u, "(1x,F7.5)", advance="no")  mci_work%x(i)
       if (i == mci_work%config%n_par_sf) &
            write (u, "(1x,'|')", advance="no")
    end do
    write (u, *)
    if (associated (mci_work%mci)) then
       call mci_work%mci%write (u, pacify = testflag)
       call mci_work%counter%write (u)
    end if
  end subroutine mci_work_write

  module subroutine mci_work_final (mci_work)
    class(mci_work_t), intent(inout) :: mci_work
    if (associated (mci_work%mci)) then
       call mci_work%mci%final ()
       deallocate (mci_work%mci)
    end if
  end subroutine mci_work_final

  module subroutine mci_work_init (mci_work, mci_entry)
    class(mci_work_t), intent(out) :: mci_work
    type(process_mci_entry_t), intent(in), target :: mci_entry
    mci_work%config => mci_entry
    allocate (mci_work%x (mci_entry%n_par))
    if (allocated (mci_entry%mci)) then
       call mci_entry%mci%allocate_instance (mci_work%mci)
       call mci_work%mci%init (mci_entry%mci)
    end if
  end subroutine mci_work_init

  module subroutine mci_work_set (mci_work, x)
    class(mci_work_t), intent(inout) :: mci_work
    real(default), dimension(:), intent(in) :: x
    mci_work%x = x
  end subroutine mci_work_set

  module subroutine mci_work_set_x_strfun (mci_work, x)
    class(mci_work_t), intent(inout) :: mci_work
    real(default), dimension(:), intent(in) :: x
    mci_work%x(1 : mci_work%config%n_par_sf) = x
  end subroutine mci_work_set_x_strfun

  module subroutine mci_work_set_x_process (mci_work, x)
    class(mci_work_t), intent(inout) :: mci_work
    real(default), dimension(:), intent(in) :: x
    mci_work%x(mci_work%config%n_par_sf + 1 : mci_work%config%n_par) = x
  end subroutine mci_work_set_x_process

  module function mci_work_get_active_components (mci_work) result (i_component)
    class(mci_work_t), intent(in) :: mci_work
    integer, dimension(:), allocatable :: i_component
    allocate (i_component (size (mci_work%config%i_component)))
    i_component = mci_work%config%i_component
  end function mci_work_get_active_components

  pure module function mci_work_get_x_strfun (mci_work) result (x)
    class(mci_work_t), intent(in) :: mci_work
    real(default), dimension(mci_work%config%n_par_sf) :: x
    x = mci_work%x(1 : mci_work%config%n_par_sf)
  end function mci_work_get_x_strfun

  pure module function mci_work_get_x_process (mci_work) result (x)
    class(mci_work_t), intent(in) :: mci_work
    real(default), dimension(mci_work%config%n_par_phs) :: x
    x = mci_work%x(mci_work%config%n_par_sf + 1 : mci_work%config%n_par)
  end function mci_work_get_x_process

  module subroutine mci_work_init_simulation &
       (mci_work, safety_factor, keep_failed_events)
    class(mci_work_t), intent(inout) :: mci_work
    real(default), intent(in), optional :: safety_factor
    logical, intent(in), optional :: keep_failed_events
    call mci_work%mci%init_simulation (safety_factor)
    call mci_work%counter%reset ()
    if (present (keep_failed_events)) &
       mci_work%keep_failed_events = keep_failed_events
  end subroutine mci_work_init_simulation

  module subroutine mci_work_final_simulation (mci_work)
    class(mci_work_t), intent(inout) :: mci_work
    call mci_work%mci%final_simulation ()
  end subroutine mci_work_final_simulation

  module subroutine mci_work_reset_counter (mci_work)
    class(mci_work_t), intent(inout) :: mci_work
    call mci_work%counter%reset ()
  end subroutine mci_work_reset_counter

  module subroutine mci_work_record_call (mci_work, status)
    class(mci_work_t), intent(inout) :: mci_work
    integer, intent(in) :: status
    call mci_work%counter%record (status)
  end subroutine mci_work_record_call

  pure module function mci_work_get_counter (mci_work) result (counter)
    class(mci_work_t), intent(in) :: mci_work
    type(process_counter_t) :: counter
    counter = mci_work%counter
  end function mci_work_get_counter


end submodule process_mci_s

