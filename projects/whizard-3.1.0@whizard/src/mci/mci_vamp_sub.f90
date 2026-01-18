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

submodule (mci_vamp) mci_vamp_s

  use io_units
  use constants, only: zero
  use format_utils, only: pac_fmt
  use format_utils, only: write_separator
  use format_defs, only: FMT_12, FMT_14, FMT_17, FMT_19
  use md5

  implicit none

contains

  module subroutine grid_parameters_write (object, unit)
    class(grid_parameters_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A,I0)") "threshold_calls       = ", &
         object%threshold_calls
    write (u, "(3x,A,I0)") "min_calls_per_channel = ", &
         object%min_calls_per_channel
    write (u, "(3x,A,I0)") "min_calls_per_bin     = ", &
         object%min_calls_per_bin
    write (u, "(3x,A,I0)") "min_bins              = ", &
         object%min_bins
    write (u, "(3x,A,I0)") "max_bins              = ", &
         object%max_bins
    write (u, "(3x,A,L1)") "stratified            = ", &
         object%stratified
    write (u, "(3x,A,L1)") "use_vamp_equivalences = ", &
         object%use_vamp_equivalences
    write (u, "(3x,A,F10.7)") "channel_weights_power = ", &
         object%channel_weights_power
    if (object%accuracy_goal > 0) then
       write (u, "(3x,A,F10.7)") "accuracy_goal         = ", &
            object%accuracy_goal
    end if
    if (object%error_goal > 0) then
       write (u, "(3x,A,F10.7)") "error_goal            = ", &
            object%error_goal
    end if
    if (object%rel_error_goal > 0) then
       write (u, "(3x,A,F10.7)") "rel_error_goal        = ", &
            object%rel_error_goal
    end if
  end subroutine grid_parameters_write

  module subroutine history_parameters_write (object, unit)
    class(history_parameters_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A,L1)") "history(global)       = ", object%global
    write (u, "(3x,A,L1)") "history(global) verb. = ", object%global_verbose
    write (u, "(3x,A,L1)") "history(channels)     = ", object%channel
    write (u, "(3x,A,L1)") "history(chann.) verb. = ", object%channel_verbose
  end subroutine history_parameters_write

  module subroutine pass_final (object)
    class(pass_t), intent(inout) :: object
    if (allocated (object%v_history)) then
       call vamp_delete_history (object%v_history)
    end if
    if (allocated (object%v_histories)) then
       call vamp_delete_history (object%v_histories)
    end if
  end subroutine pass_final

  module subroutine pass_write (object, unit, pacify)
    class(pass_t), intent(in) :: object
    integer, intent(in) :: unit
    logical, intent(in), optional :: pacify
    integer :: u, i
    character(len=7) :: fmt
    call pac_fmt (fmt, FMT_17, FMT_14, pacify)
    u = given_output_unit (unit)
    write (u, "(3x,A,I0)")  "n_it          = ", object%n_it
    write (u, "(3x,A,I0)")  "n_calls       = ", object%n_calls
    write (u, "(3x,A,I0)")  "n_bins        = ", object%n_bins
    write (u, "(3x,A,L1)")  "adapt grids   = ", object%adapt_grids
    write (u, "(3x,A,L1)")  "adapt weights = ", object%adapt_weights
    if (object%integral_defined) then
       write (u, "(3x,A)")  "Results:  [it, calls, valid, integral, error, efficiency]"
       do i = 1, object%n_it
          write (u, "(5x,I0,2(1x,I0),3(1x," // fmt // "))") &
               i, object%calls(i), object%calls_valid(i), object%integral(i), object%error(i), &
               object%efficiency(i)
       end do
    else
       write (u, "(3x,A)")  "Results: [undefined]"
    end if
  end subroutine pass_write

  module subroutine pass_read (object, u, n_pass, n_it)
    class(pass_t), intent(out) :: object
    integer, intent(in) :: u, n_pass, n_it
    integer :: i, j
    character(80) :: buffer
    object%i_pass = n_pass + 1
    object%i_first_it = n_it + 1
    call read_ival (u, object%n_it)
    call read_ival (u, object%n_calls)
    call read_ival (u, object%n_bins)
    call read_lval (u, object%adapt_grids)
    call read_lval (u, object%adapt_weights)
    allocate (object%calls (object%n_it), source = 0)
    allocate (object%calls_valid (object%n_it), source = 0)
    allocate (object%integral (object%n_it), source = 0._default)
    allocate (object%error (object%n_it), source = 0._default)
    allocate (object%efficiency (object%n_it), source = 0._default)
    read (u, "(A)")  buffer
    select case (trim (adjustl (buffer)))
    case ("Results:  [it, calls, valid, integral, error, efficiency]")
       do i = 1, object%n_it
          read (u, *) &
               j, object%calls(i), object%calls_valid(i), object%integral(i), object%error(i), &
               object%efficiency(i)
       end do
       object%integral_defined = .true.
    case ("Results: [undefined]")
       object%integral_defined = .false.
    case default
       call msg_fatal ("Reading integration pass: corrupted file")
    end select
  end subroutine pass_read

  module subroutine pass_write_history (pass, unit)
    class(pass_t), intent(in) :: pass
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (allocated (pass%v_history)) then
       call vamp_write_history (u, pass%v_history)
    else
       write (u, "(1x,A)")  "Global history: [undefined]"
    end if
    if (allocated (pass%v_histories)) then
       write (u, "(1x,A)")  "Channel histories:"
       call vamp_write_history (u, pass%v_histories)
    else
       write (u, "(1x,A)")  "Channel histories: [undefined]"
    end if
  end subroutine pass_write_history

  module subroutine pass_configure (pass, n_it, n_calls, min_calls, &
       min_bins, max_bins, min_channel_calls)
    class(pass_t), intent(inout) :: pass
    integer, intent(in) :: n_it, n_calls, min_channel_calls
    integer, intent(in) :: min_calls, min_bins, max_bins
    pass%n_it = n_it
    if (min_calls /= 0) then
       pass%n_bins =  max (min_bins, &
            min (n_calls / min_calls, max_bins))
    else
       pass%n_bins = max_bins
    end if
    pass%n_calls = max (n_calls, max (min_calls, min_channel_calls))
    if (pass%n_calls /= n_calls) then
       write (msg_buffer, "(A,I0)")  "VAMP: too few calls, resetting " &
            // "n_calls to ", pass%n_calls
       call msg_warning ()
    end if
    allocate (pass%calls (n_it), source = 0)
    allocate (pass%calls_valid (n_it), source = 0)
    allocate (pass%integral (n_it), source = 0._default)
    allocate (pass%error (n_it), source = 0._default)
    allocate (pass%efficiency (n_it), source = 0._default)
  end subroutine pass_configure

  module subroutine pass_configure_history (pass, n_channels, par)
    class(pass_t), intent(inout) :: pass
    integer, intent(in) :: n_channels
    type(history_parameters_t), intent(in) :: par
    if (par%global) then
       allocate (pass%v_history (pass%n_it))
       call vamp_create_history (pass%v_history, &
            verbose = par%global_verbose)
    end if
    if (par%channel) then
       allocate (pass%v_histories (pass%n_it, n_channels))
       call vamp_create_history (pass%v_histories, &
            verbose = par%channel_verbose)
    end if
  end subroutine pass_configure_history

  module function pass_matches (pass, ref) result (ok)
    type(pass_t), intent(in) :: pass, ref
    integer :: n
    logical :: ok
    ok = .true.
    if (ok)  ok = pass%i_pass == ref%i_pass
    if (ok)  ok = pass%i_first_it == ref%i_first_it
    if (ok)  ok = pass%n_it == ref%n_it
    if (ok)  ok = pass%n_calls == ref%n_calls
    if (ok)  ok = pass%n_bins == ref%n_bins
    if (ok)  ok = pass%adapt_grids .eqv. ref%adapt_grids
    if (ok)  ok = pass%adapt_weights .eqv. ref%adapt_weights
    if (ok)  ok = pass%integral_defined .eqv. ref%integral_defined
    if (pass%integral_defined) then
       n = pass%n_it
       if (ok)  ok = all (pass%calls(:n) == ref%calls(:n))
       if (ok)  ok = all (pass%calls_valid(:n) == ref%calls_valid (:n))
       if (ok)  ok = all (pass%integral(:n) .matches. ref%integral(:n))
       if (ok)  ok = all (pass%error(:n) .matches. ref%error(:n))
       if (ok)  ok = all (pass%efficiency(:n) .matches. ref%efficiency(:n))
    end if
  end function pass_matches

  module subroutine pass_update (pass, ref, ok)
    class(pass_t), intent(inout) :: pass
    type(pass_t), intent(in) :: ref
    logical, intent(out) :: ok
    integer :: n, n_ref
    ok = .true.
    if (ok)  ok = pass%i_pass == ref%i_pass
    if (ok)  ok = pass%i_first_it == ref%i_first_it
    if (ok)  ok = pass%n_calls == ref%n_calls
    if (ok)  ok = pass%n_bins == ref%n_bins
    if (ok)  ok = pass%adapt_grids .eqv. ref%adapt_grids
    if (ok)  ok = pass%adapt_weights .eqv. ref%adapt_weights
    if (ok) then
       if (ref%integral_defined) then
          if (.not. allocated (pass%calls)) then
             allocate (pass%calls (pass%n_it), source = 0)
             allocate (pass%calls_valid (pass%n_it), source = 0)
             allocate (pass%integral (pass%n_it), source = 0._default)
             allocate (pass%error (pass%n_it), source = 0._default)
             allocate (pass%efficiency (pass%n_it), source = 0._default)
          end if
          n = count (pass%calls /= 0)
          n_ref = count (ref%calls /= 0)
          ok = n <= n_ref .and. n_ref <= pass%n_it
          if (ok)  ok = all (pass%calls(:n) == ref%calls(:n))
          if (ok)  ok = all (pass%calls_valid(:n) == ref%calls_valid(:n))
          if (ok)  ok = all (pass%integral(:n) .matches. ref%integral(:n))
          if (ok)  ok = all (pass%error(:n) .matches. ref%error(:n))
          if (ok)  ok = all (pass%efficiency(:n) .matches. ref%efficiency(:n))
          if (ok) then
             pass%calls(n+1:n_ref) = ref%calls(n+1:n_ref)
             pass%calls_valid(n+1:n_ref) = ref%calls_valid(n+1:n_ref)
             pass%integral(n+1:n_ref) = ref%integral(n+1:n_ref)
             pass%error(n+1:n_ref) = ref%error(n+1:n_ref)
             pass%efficiency(n+1:n_ref) = ref%efficiency(n+1:n_ref)
             pass%integral_defined = any (pass%calls /= 0)
          end if
       end if
    end if
  end subroutine pass_update

  elemental module function real_matches (x, y) result (ok)
    real(default), intent(in) :: x, y
    logical :: ok
    real(default), parameter :: tolerance = 1.e-8_default
    ok = abs (x - y) <= tolerance * max (abs (x), abs (y))
  end function real_matches

  module function pass_get_integration_index (pass) result (n)
    class (pass_t), intent(in) :: pass
    integer :: n
    integer :: i
    n = 0
    if (allocated (pass%calls)) then
       do i = 1, pass%n_it
          if (pass%calls(i) == 0)  exit
          n = i
       end do
    end if
  end function pass_get_integration_index

  module function pass_get_calls (pass) result (calls)
    class(pass_t), intent(in) :: pass
    integer :: calls
    integer :: n
    n = pass%get_integration_index ()
    if (n /= 0) then
       calls = pass%calls(n)
    else
       calls = 0
    end if
  end function pass_get_calls

  module function pass_get_calls_valid (pass) result (calls_valid)
    class(pass_t), intent(in) :: pass
    integer :: calls_valid
    integer :: n
    n = pass%get_integration_index ()
    if (n /= 0) then
       calls_valid = pass%calls_valid(n)
    else
       calls_valid = 0
    end if
  end function pass_get_calls_valid

  module function pass_get_integral (pass) result (integral)
    class(pass_t), intent(in) :: pass
    real(default) :: integral
    integer :: n
    n = pass%get_integration_index ()
    if (n /= 0) then
       integral = pass%integral(n)
    else
       integral = 0
    end if
  end function pass_get_integral

  module function pass_get_error (pass) result (error)
    class(pass_t), intent(in) :: pass
    real(default) :: error
    integer :: n
    n = pass%get_integration_index ()
    if (n /= 0) then
       error = pass%error(n)
    else
       error = 0
    end if
  end function pass_get_error

  module function pass_get_efficiency (pass) result (efficiency)
    class(pass_t), intent(in) :: pass
    real(default) :: efficiency
    integer :: n
    n = pass%get_integration_index ()
    if (n /= 0) then
       efficiency = pass%efficiency(n)
    else
       efficiency = 0
    end if
  end function pass_get_efficiency

  module subroutine mci_vamp_reset (object)
    class(mci_vamp_t), intent(inout) :: object
    type(pass_t), pointer :: current_pass
    do while (associated (object%first_pass))
       current_pass => object%first_pass
       object%first_pass => current_pass%next
       call current_pass%final ()
       deallocate (current_pass)
    end do
    object%current_pass => null ()
  end subroutine mci_vamp_reset

  module subroutine mci_vamp_final (object)
    class(mci_vamp_t), intent(inout) :: object
    call object%reset ()
    call vamp_equivalences_final (object%equivalences)
    call object%base_final ()
  end subroutine mci_vamp_final

  module subroutine mci_vamp_write (object, unit, pacify, md5sum_version)
    class(mci_vamp_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: pacify
    logical, intent(in), optional :: md5sum_version
    type(pass_t), pointer :: current_pass
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "VAMP integrator:"
    call object%base_write (u, pacify, md5sum_version)
    if (allocated (object%dim_is_flat)) then
       write (u, "(3x,A,999(1x,I0))")  "Flat dimensions    =", &
            pack ([(i, i = 1, object%n_dim)], object%dim_is_flat)
    end if
    write (u, "(1x,A)")  "Grid parameters:"
    call object%grid_par%write (u)
    write (u, "(3x,A,I0)") "min_calls             = ", object%min_calls
    write (u, "(3x,A,L1)") "negative weights      = ", &
         object%negative_weights
    write (u, "(3x,A,L1)") "verbose               = ", &
         object%verbose
    if (object%grid_par%use_vamp_equivalences) then
       call vamp_equivalences_write (object%equivalences, u)
    end if
    current_pass => object%first_pass
    do while (associated (current_pass))
       write (u, "(1x,A,I0,A)")  "Integration pass:"
       call current_pass%write (u, pacify)
       current_pass => current_pass%next
    end do
    if (object%md5sum_adapted /= "") then
       write (u, "(1x,A,A,A)")  "MD5 sum (including results) = '", &
            object%md5sum_adapted, "'"
    end if
  end subroutine mci_vamp_write

  module subroutine mci_vamp_write_history_parameters (mci, unit)
    class(mci_vamp_t), intent(in) :: mci
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "VAMP history parameters:"
    call mci%history_par%write (unit)
  end subroutine mci_vamp_write_history_parameters

  module subroutine mci_vamp_write_history (mci, unit)
    class(mci_vamp_t), intent(in) :: mci
    integer, intent(in), optional :: unit
    type(pass_t), pointer :: current_pass
    integer :: i_pass
    integer :: u
    u = given_output_unit (unit)
    if (associated (mci%first_pass)) then
       write (u, "(1x,A)")  "VAMP history (global):"
       i_pass = 0
       current_pass => mci%first_pass
       do while (associated (current_pass))
          i_pass = i_pass + 1
          write (u, "(1x,A,I0,':')")  "Pass #", i_pass
          call current_pass%write_history (u)
          current_pass => current_pass%next
       end do
    end if
  end subroutine mci_vamp_write_history

  module subroutine mci_vamp_compute_md5sum (mci, pacify)
    class(mci_vamp_t), intent(inout) :: mci
    logical, intent(in), optional :: pacify
    integer :: u
    mci%md5sum_adapted = ""
    u = free_unit ()
    open (u, status = "scratch", action = "readwrite")
    write (u, "(A)")  mci%md5sum
    call mci%write (u, pacify, md5sum_version = .true.)
    rewind (u)
    mci%md5sum_adapted = md5sum (u)
    close (u)
  end subroutine mci_vamp_compute_md5sum

  pure module function mci_vamp_get_md5sum (mci) result (md5sum)
    class(mci_vamp_t), intent(in) :: mci
    character(32) :: md5sum
    if (mci%md5sum_adapted /= "") then
       md5sum = mci%md5sum_adapted
    else
       md5sum = mci%md5sum
    end if
  end function mci_vamp_get_md5sum

  module subroutine mci_vamp_startup_message (mci, unit, n_calls)
    class(mci_vamp_t), intent(in) :: mci
    integer, intent(in), optional :: unit, n_calls
    integer :: num_calls, n_bins
    if (present (n_calls)) then
       num_calls = n_calls
    else
       num_calls = 0
    end if
    if (mci%min_calls /= 0) then
       n_bins =  max (mci%grid_par%min_bins, &
            min (num_calls / mci%min_calls, &
            mci%grid_par%max_bins))
    else
       n_bins = mci%grid_par%max_bins
    end if
    call mci%base_startup_message (unit = unit, n_calls = n_calls)
    if (mci%grid_par%use_vamp_equivalences) then
       write (msg_buffer, "(A,2(1x,I0,1x,A))") &
            "Integrator: Using VAMP channel equivalences"
       call msg_message (unit = unit)
    end if
    write (msg_buffer, "(A,2(1x,I0,1x,A),L1)") &
         "Integrator:", num_calls, &
         "initial calls,", n_bins, &
         "bins, stratified = ", &
         mci%grid_par%stratified
    call msg_message (unit = unit)
    write (msg_buffer, "(A,2(1x,I0,1x,A))") &
         "Integrator: VAMP"
    call msg_message (unit = unit)
  end subroutine mci_vamp_startup_message

  module subroutine mci_vamp_write_log_entry (mci, u)
    class(mci_vamp_t), intent(in) :: mci
    integer, intent(in) :: u
    write (u, "(1x,A)")  "MC Integrator is VAMP"
    call write_separator (u)
    call mci%write_history (u)
    call write_separator (u)
    if (mci%grid_par%use_vamp_equivalences) then
       call vamp_equivalences_write (mci%equivalences, u)
    else
       write (u, "(3x,A)") "No VAMP equivalences have been used"
    end if
    call write_separator (u)
    call mci%write_chain_weights (u)
  end subroutine mci_vamp_write_log_entry

  module subroutine mci_vamp_record_index (mci, i_mci)
    class(mci_vamp_t), intent(inout) :: mci
    integer, intent(in) :: i_mci
    type(string_t) :: basename, suffix
    character(32) :: buffer
    if (mci%grid_filename_set) then
       basename = mci%grid_filename
       call split (basename, suffix, ".", back=.true.)
       write (buffer, "(I0)")  i_mci
       if (basename /= "") then
          mci%grid_filename = basename // ".m" // trim (buffer) // "." // suffix
       else
          mci%grid_filename = suffix // ".m" // trim (buffer) // ".vg"
       end if
    end if
  end subroutine mci_vamp_record_index

  module subroutine mci_vamp_set_grid_parameters (mci, grid_par)
    class(mci_vamp_t), intent(inout) :: mci
    type(grid_parameters_t), intent(in) :: grid_par
    mci%grid_par = grid_par
    mci%min_calls = grid_par%min_calls_per_bin * mci%n_channel
  end subroutine mci_vamp_set_grid_parameters

  module subroutine mci_vamp_set_history_parameters (mci, history_par)
    class(mci_vamp_t), intent(inout) :: mci
    type(history_parameters_t), intent(in) :: history_par
    mci%history_par = history_par
  end subroutine mci_vamp_set_history_parameters

  module subroutine mci_vamp_set_rebuild_flag (mci, rebuild, check_grid_file)
    class(mci_vamp_t), intent(inout) :: mci
    logical, intent(in) :: rebuild
    logical, intent(in) :: check_grid_file
    mci%rebuild = rebuild
    mci%check_grid_file = check_grid_file
  end subroutine mci_vamp_set_rebuild_flag

  module subroutine mci_vamp_set_grid_filename (mci, name, run_id)
    class(mci_vamp_t), intent(inout) :: mci
    type(string_t), intent(in) :: name
    type(string_t), intent(in), optional :: run_id
    if (present (run_id)) then
       mci%grid_filename = name // "." // run_id // ".vg"
    else
       mci%grid_filename = name // ".vg"
    end if
    mci%grid_filename_set = .true.
  end subroutine mci_vamp_set_grid_filename

  module subroutine mci_vamp_prepend_grid_path (mci, prefix)
    class(mci_vamp_t), intent(inout) :: mci
    type(string_t), intent(in) :: prefix
    if (mci%grid_filename_set) then
       mci%grid_filename = prefix // "/" // mci%grid_filename
    else
       call msg_warning ("Cannot add prefix to invalid grid filename!")
    end if
  end subroutine mci_vamp_prepend_grid_path

  module subroutine mci_vamp_declare_flat_dimensions (mci, dim_flat)
    class(mci_vamp_t), intent(inout) :: mci
    integer, dimension(:), intent(in) :: dim_flat
    integer :: d
    allocate (mci%dim_is_flat (mci%n_dim), source = .false.)
    do d = 1, size (dim_flat)
       mci%dim_is_flat(dim_flat(d)) = .true.
    end do
  end subroutine mci_vamp_declare_flat_dimensions

  module subroutine mci_vamp_declare_equivalences (mci, channel, dim_offset)
    class(mci_vamp_t), intent(inout) :: mci
    type(phs_channel_t), dimension(:), intent(in) :: channel
    integer, intent(in) :: dim_offset
    integer, dimension(:), allocatable :: perm, mode
    integer :: n_channels, n_dim, n_equivalences
    integer :: c, i, j, left, right
    integer :: n_dim_perm
    n_channels = mci%n_channel
    n_dim = mci%n_dim
    n_equivalences = 0
    do c = 1, n_channels
       n_equivalences = n_equivalences + size (channel(c)%eq)
    end do
    call vamp_equivalences_init (mci%equivalences, &
         n_equivalences, n_channels, n_dim)
    allocate (perm (n_dim))
    allocate (mode (n_dim))
    perm = [(i, i = 1, n_dim)]
    mode = VEQ_IDENTITY
    c = 1
    j = 0
    do i = 1, n_equivalences
       if (j < size (channel(c)%eq)) then
          j = j + 1
       else
          c = c + 1
          j = 1
       end if
       associate (eq => channel(c)%eq(j))
         left = c
         right = eq%c
         n_dim_perm = size (eq%perm)
         perm(dim_offset + 1:dim_offset + n_dim_perm) = eq%perm + dim_offset
         mode(dim_offset + 1:dim_offset + n_dim_perm) = eq%mode
         call vamp_equivalence_set (mci%equivalences, &
              i, left, right, perm, mode)
       end associate
    end do
    call vamp_equivalences_complete (mci%equivalences)
  end subroutine mci_vamp_declare_equivalences

  module subroutine mci_vamp_add_pass &
       (mci, adapt_grids, adapt_weights, final_pass)
    class(mci_vamp_t), intent(inout) :: mci
    logical, intent(in), optional :: adapt_grids, adapt_weights, final_pass
    integer :: i_pass, i_it
    type(pass_t), pointer :: new
    allocate (new)
    if (associated (mci%current_pass)) then
       i_pass = mci%current_pass%i_pass + 1
       i_it   = mci%current_pass%i_first_it + mci%current_pass%n_it
       mci%current_pass%next => new
    else
       i_pass = 1
       i_it = 1
       mci%first_pass => new
    end if
    mci%current_pass => new
    new%i_pass = i_pass
    new%i_first_it = i_it
    if (present (adapt_grids)) then
       new%adapt_grids = adapt_grids
    else
       new%adapt_grids = .false.
    end if
    if (present (adapt_weights)) then
       new%adapt_weights = adapt_weights
    else
       new%adapt_weights = .false.
    end if
    if (present (final_pass)) then
       new%is_final_pass = final_pass
    else
       new%is_final_pass = .false.
    end if
  end subroutine mci_vamp_add_pass

  module subroutine mci_vamp_update_from_ref (mci, mci_ref, success)
    class(mci_vamp_t), intent(inout) :: mci
    class(mci_t), intent(in) :: mci_ref
    logical, intent(out) :: success
    type(pass_t), pointer :: current_pass, ref_pass
    select type (mci_ref)
    type is (mci_vamp_t)
       current_pass => mci%first_pass
       ref_pass => mci_ref%first_pass
       success = .true.
       do while (success .and. associated (current_pass))
          if (associated (ref_pass)) then
             if (associated (current_pass%next)) then
                success = current_pass .matches. ref_pass
             else
                call current_pass%update (ref_pass, success)
                if (current_pass%integral_defined) then
                   mci%integral = current_pass%get_integral ()
                   mci%error = current_pass%get_error ()
                   mci%efficiency = current_pass%get_efficiency ()
                   mci%integral_known = .true.
                   mci%error_known = .true.
                   mci%efficiency_known = .true.
                end if
             end if
             current_pass => current_pass%next
             ref_pass => ref_pass%next
          else
             success = .false.
          end if
       end do
    end select
  end subroutine mci_vamp_update_from_ref

  module subroutine mci_vamp_update (mci, u, success)
    class(mci_vamp_t), intent(inout) :: mci
    integer, intent(in) :: u
    logical, intent(out) :: success
    character(80) :: buffer
    character(32) :: md5sum_file
    type(mci_vamp_t) :: mci_file
    integer :: n_pass, n_it
    call read_sval (u, md5sum_file)
    if (mci%check_grid_file) then
       success = md5sum_file == mci%md5sum
    else
       success = .true.
    end if
    if (success) then
       read (u, *)
       read (u, "(A)")  buffer
       if (trim (adjustl (buffer)) == "VAMP integrator:") then
          n_pass = 0
          n_it = 0
          do
             read (u, "(A)")  buffer
             select case (trim (adjustl (buffer)))
             case ("")
                exit
             case ("Integration pass:")
                call mci_file%add_pass ()
                call mci_file%current_pass%read (u, n_pass, n_it)
                n_pass = n_pass + 1
                n_it = n_it + mci_file%current_pass%n_it
             end select
          end do
          call mci%update_from_ref (mci_file, success)
          call mci_file%final ()
       else
          call msg_fatal ("VAMP: reading grid file: corrupted data")
       end if
    end if
  end subroutine mci_vamp_update

  module subroutine mci_vamp_write_grids (mci, instance)
    class(mci_vamp_t), intent(in) :: mci
    class(mci_instance_t), intent(inout) :: instance
    integer :: u
    select type (instance)
    type is (mci_vamp_instance_t)
       if (mci%grid_filename_set) then
          if (instance%grids_defined) then
             u = free_unit ()
             open (u, file = char (mci%grid_filename), &
                  action = "write", status = "replace")
             write (u, "(1x,A,A,A)")  "MD5sum = '", mci%md5sum, "'"
             write (u, *)
             call mci%write (u)
             write (u, *)
             write (u, "(1x,A)")  "VAMP grids:"
             call vamp_write_grids (instance%grids, u, &
                  write_integrals = .true.)
             close (u)
          else
             call msg_bug ("VAMP: write grids: grids undefined")
          end if
       else
          call msg_bug ("VAMP: write grids: filename undefined")
       end if
    end select
  end subroutine mci_vamp_write_grids

  module subroutine mci_vamp_read_grids_header (mci, success)
    class(mci_vamp_t), intent(inout) :: mci
    logical, intent(out) :: success
    logical :: exist
    integer :: u
    success = .false.
    if (mci%grid_filename_set) then
       inquire (file = char (mci%grid_filename), exist = exist)
       if (exist) then
          u = free_unit ()
          open (u, file = char (mci%grid_filename), &
               action = "read", status = "old")
          call mci%update (u, success)
          close (u)
          if (.not. success) then
             write (msg_buffer, "(A,A,A)") &
                  "VAMP: parameter mismatch, discarding grid file '", &
                  char (mci%grid_filename), "'"
             call msg_message ()
          end if
       end if
    else
       call msg_bug ("VAMP: read grids: filename undefined")
    end if
  end subroutine mci_vamp_read_grids_header

  module subroutine mci_vamp_read_grids_data (mci, instance, read_integrals)
    class(mci_vamp_t), intent(in) :: mci
    class(mci_instance_t), intent(inout) :: instance
    logical, intent(in), optional :: read_integrals
    integer :: u
    character(80) :: buffer
    select type (instance)
    type is (mci_vamp_instance_t)
       if (.not. instance%grids_defined) then
          u = free_unit ()
          open (u, file = char (mci%grid_filename), &
               action = "read", status = "old")
          do
             read (u, "(A)")  buffer
             if (trim (adjustl (buffer)) == "VAMP grids:")  exit
          end do
          call vamp_read_grids (instance%grids, u, read_integrals)
          close (u)
          call instance%set_channel_weights (instance%grids%weights)
          instance%grids_defined = .true.
       else
          call msg_bug ("VAMP: read grids: grids already defined")
       end if
    end select
  end subroutine mci_vamp_read_grids_data

  module subroutine mci_vamp_read_grids (mci, instance, success)
    class(mci_vamp_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout) :: instance
    logical, intent(out) :: success
    logical :: exist
    integer :: u
    character(80) :: buffer
    select type (instance)
    type is (mci_vamp_instance_t)
       success = .false.
       if (mci%grid_filename_set) then
          if (.not. instance%grids_defined) then
             inquire (file = char (mci%grid_filename), exist = exist)
             if (exist) then
                u = free_unit ()
                open (u, file = char (mci%grid_filename), &
                     action = "read", status = "old")
                call mci%update (u, success)
                if (success) then
                   read (u, "(A)")  buffer
                   if (trim (adjustl (buffer)) == "VAMP grids:") then
                      call vamp_read_grids (instance%grids, u)
                   else
                      call msg_fatal ("VAMP: reading grid file: &
                           &corrupted grid data")
                   end if
                else
                   write (msg_buffer, "(A,A,A)") &
                        "VAMP: parameter mismatch, discarding grid file '", &
                        char (mci%grid_filename), "'"
                   call msg_message ()
                end if
                close (u)
                instance%grids_defined = success
             end if
          else
             call msg_bug ("VAMP: read grids: grids already defined")
          end if
       else
          call msg_bug ("VAMP: read grids: filename undefined")
       end if
    end select
  end subroutine mci_vamp_read_grids

  subroutine read_rval (u, rval)
    integer, intent(in) :: u
    real(default), intent(out) :: rval
    character(80) :: buffer
    read (u, "(A)")  buffer
    buffer = adjustl (buffer(scan (buffer, "=") + 1:))
    read (buffer, *)  rval
  end subroutine read_rval

  subroutine read_ival (u, ival)
    integer, intent(in) :: u
    integer, intent(out) :: ival
    character(80) :: buffer
    read (u, "(A)")  buffer
    buffer = adjustl (buffer(scan (buffer, "=") + 1:))
    read (buffer, *)  ival
  end subroutine read_ival

  subroutine read_sval (u, sval)
    integer, intent(in) :: u
    character(*), intent(out) :: sval
    character(80) :: buffer
    read (u, "(A)")  buffer
    buffer = adjustl (buffer(scan (buffer, "=") + 1:))
    read (buffer, *)  sval
  end subroutine read_sval

  subroutine read_lval (u, lval)
    integer, intent(in) :: u
    logical, intent(out) :: lval
    character(80) :: buffer
    read (u, "(A)")  buffer
    buffer = adjustl (buffer(scan (buffer, "=") + 1:))
    read (buffer, *)  lval
  end subroutine read_lval

  module subroutine mci_vamp_integrate (mci, instance, sampler, &
       n_it, n_calls, results, pacify)
    class(mci_vamp_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout), target :: instance
    class(mci_sampler_t), intent(inout), target :: sampler
    integer, intent(in) :: n_it
    integer, intent(in) :: n_calls
    class(mci_results_t), intent(inout), optional :: results
    logical, intent(in), optional :: pacify
    integer :: it
    logical :: reshape, from_file, success
    select type (instance)
    type is (mci_vamp_instance_t)
       if (associated (mci%current_pass)) then
          mci%current_pass%integral_defined = .false.
          call mci%current_pass%configure (n_it, n_calls, &
               mci%min_calls, mci%grid_par%min_bins, &
               mci%grid_par%max_bins, &
               mci%grid_par%min_calls_per_channel * mci%n_channel)
          call mci%current_pass%configure_history &
               (mci%n_channel, mci%history_par)
          instance%pass_complete = .false.
          instance%it_complete = .false.
          call instance%new_pass (reshape)
          if (.not. instance%grids_defined .or. instance%grids_from_file) then
             if (mci%grid_filename_set .and. .not. mci%rebuild) then
                call mci%read_grids_header (success)
                from_file = success
                if (.not. instance%grids_defined .and. success) then
                   call mci%read_grids_data (instance)
                end if
             else
                from_file = .false.
             end if
          else
             from_file = .false.
          end if
          if (from_file) then
             if (.not. mci%check_grid_file) &
                  call msg_warning ("Reading grid file: MD5 sum check disabled")
             call msg_message ("VAMP: " &
                  // "using grids and results from file '" &
                  // char (mci%grid_filename) // "'")
          else if (.not. instance%grids_defined) then
             call instance%create_grids ()
          end if
          do it = 1, instance%n_it
             if (signal_is_pending ())  return
             reshape = reshape .or. &
                  (instance%grids_from_file .and. n_it > mci%current_pass%get_integration_index ())
             instance%grids_from_file = from_file .and. &
                  it <= mci%current_pass%get_integration_index ()
             if (.not. instance%grids_from_file) then
                instance%it_complete = .false.
                call instance%adapt_grids ()
                if (signal_is_pending ())  return
                call instance%adapt_weights ()
                if (signal_is_pending ())  return
                call instance%discard_integrals (reshape)
                if (mci%grid_par%use_vamp_equivalences) then
                   call instance%sample_grids (mci%rng, sampler, &
                        mci%equivalences)
                else
                   call instance%sample_grids (mci%rng, sampler)
                end if
                if (signal_is_pending ())  return
                instance%it_complete = .true.
                if (instance%integral /= 0) then
                   mci%current_pass%calls(it) = instance%calls
                   mci%current_pass%calls_valid(it) = instance%calls_valid
                   mci%current_pass%integral(it) = instance%integral
                   if (abs (instance%error / instance%integral) &
                        > epsilon (1._default)) then
                      mci%current_pass%error(it) = instance%error
                   end if
                   mci%current_pass%efficiency(it) = instance%efficiency
                end if
                mci%current_pass%integral_defined = .true.
             end if
             if (present (results)) then
                if (mci%has_chains ()) then
                   call mci%collect_chain_weights (instance%w)
                   call results%record (1, &
                        n_calls        = mci%current_pass%calls(it), &
                        n_calls_valid  = mci%current_pass%calls_valid(it), &
                        integral       = mci%current_pass%integral(it), &
                        error          = mci%current_pass%error(it), &
                        efficiency     = mci%current_pass%efficiency(it), &
                        ! TODO Insert pos. and neg. Efficiency from VAMP.
                        efficiency_pos = 0._default, &
                        efficiency_neg = 0._default, &
                        chain_weights  = mci%chain_weights, &
                        suppress = pacify)
                else
                   call results%record (1, &
                        n_calls        = mci%current_pass%calls(it), &
                        n_calls_valid  = mci%current_pass%calls_valid(it), &
                        integral       = mci%current_pass%integral(it), &
                        error          = mci%current_pass%error(it), &
                        efficiency     = mci%current_pass%efficiency(it), &
                        ! TODO Insert pos. and neg. Efficiency from VAMP.
                        efficiency_pos = 0._default, &
                        efficiency_neg = 0._default, &
                        suppress = pacify)
                end if
             end if
             if (.not. instance%grids_from_file &
                  .and. mci%grid_filename_set) then
                call mci%write_grids (instance)
             end if
             call instance%allow_adaptation ()
             reshape = .false.
             if (.not. mci%current_pass%is_final_pass) then
                call mci%check_goals (it, success)
                if (success)  exit
             end if
          end do
          if (signal_is_pending ())  return
          instance%pass_complete = .true.
          mci%integral = mci%current_pass%get_integral()
          mci%error = mci%current_pass%get_error()
          mci%efficiency = mci%current_pass%get_efficiency()
          mci%integral_known = .true.
          mci%error_known = .true.
          mci%efficiency_known = .true.
          call mci%compute_md5sum (pacify)
       else
          call msg_bug ("MCI integrate: current_pass object not allocated")
       end if
    end select
  end subroutine mci_vamp_integrate

  module subroutine mci_vamp_check_goals (mci, it, success)
    class(mci_vamp_t), intent(inout) :: mci
    integer, intent(in) :: it
    logical, intent(out) :: success
    success = .false.
    if (mci%error_reached (it)) then
       mci%current_pass%n_it = it
       call msg_message ("VAMP: error goal reached; &
            &skipping iterations")
       success = .true.
       return
    end if
    if (mci%rel_error_reached (it)) then
       mci%current_pass%n_it = it
       call msg_message ("VAMP: relative error goal reached; &
            &skipping iterations")
       success = .true.
       return
    end if
    if (mci%accuracy_reached (it)) then
       mci%current_pass%n_it = it
       call msg_message ("VAMP: accuracy goal reached; &
            &skipping iterations")
       success = .true.
       return
    end if
  end subroutine mci_vamp_check_goals

  module function mci_vamp_error_reached (mci, it) result (flag)
    class(mci_vamp_t), intent(in) :: mci
    integer, intent(in) :: it
    logical :: flag
    real(default) :: error_goal, error
    error_goal = mci%grid_par%error_goal
    if (error_goal > 0) then
       associate (pass => mci%current_pass)
         if (pass%integral_defined) then
            error = abs (pass%error(it))
            flag = error < error_goal
         else
            flag = .false.
         end if
       end associate
    else
       flag = .false.
    end if
  end function mci_vamp_error_reached

  module function mci_vamp_rel_error_reached (mci, it) result (flag)
    class(mci_vamp_t), intent(in) :: mci
    integer, intent(in) :: it
    logical :: flag
    real(default) :: rel_error_goal, rel_error
    rel_error_goal = mci%grid_par%rel_error_goal
    if (rel_error_goal > 0) then
       associate (pass => mci%current_pass)
         if (pass%integral_defined) then
            if (pass%integral(it) /= 0) then
               rel_error = abs (pass%error(it) / pass%integral(it))
               flag = rel_error < rel_error_goal
            else
               flag = .true.
            end if
         else
            flag = .false.
         end if
       end associate
    else
       flag = .false.
    end if
  end function mci_vamp_rel_error_reached

  module function mci_vamp_accuracy_reached (mci, it) result (flag)
    class(mci_vamp_t), intent(in) :: mci
    integer, intent(in) :: it
    logical :: flag
    real(default) :: accuracy_goal, accuracy
    accuracy_goal = mci%grid_par%accuracy_goal
    if (accuracy_goal > 0) then
       associate (pass => mci%current_pass)
         if (pass%integral_defined) then
            if (pass%integral(it) /= 0) then
               accuracy = abs (pass%error(it) / pass%integral(it)) &
                    * sqrt (real (pass%calls(it), default))
               flag = accuracy < accuracy_goal
            else
               flag = .true.
            end if
         else
            flag = .false.
         end if
       end associate
    else
       flag = .false.
    end if
  end function mci_vamp_accuracy_reached

  module subroutine mci_vamp_prepare_simulation (mci)
    class(mci_vamp_t), intent(inout) :: mci
    logical :: success
    if (mci%grid_filename_set) then
       call mci%read_grids_header (success)
       call mci%compute_md5sum ()
       if (.not. success) then
          call msg_fatal ("Simulate: " &
               // "reading integration grids from file '" &
               // char (mci%grid_filename) // "' failed")
       end if
    else
       call msg_bug ("VAMP: simulation: no grids, no grid filename")
    end if
  end subroutine mci_vamp_prepare_simulation

  module subroutine mci_vamp_rebuild_event (mci, instance, sampler, state)
    class(mci_vamp_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout) :: instance
    class(mci_sampler_t), intent(inout) :: sampler
    class(mci_state_t), intent(in) :: state
    call msg_bug ("MCI vamp rebuild event not implemented yet")
  end subroutine mci_vamp_rebuild_event

  module subroutine mci_vamp_pacify (object, efficiency_reset, error_reset)
    class(mci_vamp_t), intent(inout) :: object
    logical, intent(in), optional :: efficiency_reset, error_reset
    logical :: err_reset
    type(pass_t), pointer :: current_pass
    err_reset = .false.
    if (present (error_reset))  err_reset = error_reset
    current_pass => object%first_pass
    do while (associated (current_pass))
       if (allocated (current_pass%error) .and. err_reset) then
          current_pass%error = 0
       end if
       if (allocated (current_pass%efficiency) .and. err_reset) then
          current_pass%efficiency = 1
       end if
       current_pass => current_pass%next
    end do
  end subroutine mci_vamp_pacify

  module subroutine mci_vamp_instance_write (object, unit, pacify)
    class(mci_vamp_instance_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: pacify
    integer :: u, i
    character(len=7) :: fmt
    call pac_fmt (fmt, FMT_17, FMT_14, pacify)
    u = given_output_unit (unit)
    write (u, "(3x,A," // FMT_19 // ")") "Integrand = ", object%integrand
    write (u, "(3x,A," // FMT_19 // ")") "Weight    = ", object%mci_weight
    if (object%vamp_weight_set) then
       write (u, "(3x,A," // FMT_19 // ")") "VAMP wgt  = ", object%vamp_weight
       if (object%vamp_excess /= 0) then
          write (u, "(3x,A," // FMT_19 // ")") "VAMP exc  = ", &
               object%vamp_excess
       end if
    end if
    write (u, "(3x,A,L1)")  "adapt grids   = ", object%enable_adapt_grids
    write (u, "(3x,A,L1)")  "adapt weights = ", object%enable_adapt_weights
    if (object%grids_defined) then
       if (object%grids_from_file) then
          write (u, "(3x,A)")  "VAMP grids: read from file"
       else
          write (u, "(3x,A)")  "VAMP grids: defined"
       end if
    else
       write (u, "(3x,A)")  "VAMP grids: [undefined]"
    end if
    write (u, "(3x,A,I0)")  "n_it          = ", object%n_it
    write (u, "(3x,A,I0)")  "it            = ", object%it
    write (u, "(3x,A,L1)")  "pass complete = ", object%it_complete
    write (u, "(3x,A,I0)")  "n_calls       = ", object%n_calls
    write (u, "(3x,A,I0)")  "calls         = ", object%calls
    write (u, "(3x,A,I0)")  "calls_valid   = ", object%calls_valid
    write (u, "(3x,A,L1)")  "it complete   = ", object%it_complete
    write (u, "(3x,A,I0)")  "n adapt.(g)   = ", object%n_adapt_grids
    write (u, "(3x,A,I0)")  "n adapt.(w)   = ", object%n_adapt_weights
    write (u, "(3x,A,L1)")  "gen. events   = ", object%generating_events
    write (u, "(3x,A,L1)")  "neg. weights  = ", object%negative_weights
    if (object%safety_factor /= 1)  write &
          (u, "(3x,A," // fmt // ")")  "safety f = ", object%safety_factor
    write (u, "(3x,A," // fmt // ")")  "integral = ", object%integral
    write (u, "(3x,A," // fmt // ")")  "error    = ", object%error
    write (u, "(3x,A," // fmt // ")")  "eff.     = ", object%efficiency
    write (u, "(3x,A)")  "weights:"
    do i = 1, size (object%w)
       write (u, "(5x,I0,1x," // FMT_12 // ")")  i, object%w(i)
    end do
  end subroutine mci_vamp_instance_write

  module subroutine mci_vamp_instance_write_grids (object, unit)
    class(mci_vamp_instance_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (object%grids_defined) then
       call vamp_write_grids (object%grids, u, write_integrals = .true.)
    end if
  end subroutine mci_vamp_instance_write_grids

  module subroutine mci_vamp_instance_final (object)
    class(mci_vamp_instance_t), intent(inout) :: object
    if (object%allocate_global_history) then
       if (associated (object%v_history)) then
          call vamp_delete_history (object%v_history)
          deallocate (object%v_history)
       end if
    end if
    if (object%allocate_channel_history) then
       if (associated (object%v_histories)) then
          call vamp_delete_history (object%v_histories)
          deallocate (object%v_histories)
       end if
    end if
    if (object%grids_defined) then
       call vamp_delete_grids (object%grids)
       object%grids_defined = .false.
    end if
  end subroutine mci_vamp_instance_final

  module subroutine mci_vamp_instance_init (mci_instance, mci)
    class(mci_vamp_instance_t), intent(out) :: mci_instance
    class(mci_t), intent(in), target :: mci
    call mci_instance%base_init (mci)
    select type (mci)
    type is (mci_vamp_t)
       mci_instance%mci => mci
       allocate (mci_instance%gi (mci%n_channel))
       mci_instance%allocate_global_history = .not. mci%history_par%global
       mci_instance%allocate_channel_history = .not. mci%history_par%channel
       mci_instance%negative_weights = mci%negative_weights
    end select
  end subroutine mci_vamp_instance_init

  module subroutine mci_vamp_instance_new_pass (instance, reshape)
    class(mci_vamp_instance_t), intent(inout) :: instance
    logical, intent(out) :: reshape
    type(pass_t), pointer :: current
    associate (mci => instance%mci)
      current => mci%current_pass
      instance%n_it = current%n_it
      if (instance%n_calls == 0) then
         reshape = .false.
         instance%n_calls = current%n_calls
      else if (instance%n_calls == current%n_calls) then
         reshape = .false.
      else
         reshape = .true.
         instance%n_calls = current%n_calls
      end if
      instance%it = 0
      instance%calls = 0
      instance%calls_valid = 0
      instance%enable_adapt_grids = current%adapt_grids
      instance%enable_adapt_weights = current%adapt_weights
      instance%generating_events = .false.
      if (instance%allocate_global_history) then
         if (associated (instance%v_history)) then
            call vamp_delete_history (instance%v_history)
            deallocate (instance%v_history)
         end if
         allocate (instance%v_history (instance%n_it))
         call vamp_create_history (instance%v_history, verbose = .false.)
      else
         instance%v_history => current%v_history
      end if
      if (instance%allocate_channel_history) then
         if (associated (instance%v_histories)) then
            call vamp_delete_history (instance%v_histories)
            deallocate (instance%v_histories)
         end if
         allocate (instance%v_histories (instance%n_it, mci%n_channel))
         call vamp_create_history (instance%v_histories, verbose = .false.)
      else
         instance%v_histories => current%v_histories
      end if
    end associate
  end subroutine mci_vamp_instance_new_pass

  module subroutine mci_vamp_instance_create_grids (instance)
    class(mci_vamp_instance_t), intent(inout) :: instance
    type (pass_t), pointer :: current
    integer, dimension(:), allocatable :: num_div
    real(default), dimension(:,:), allocatable :: region
    associate (mci => instance%mci)
      current => mci%current_pass
      allocate (num_div (mci%n_dim))
      allocate (region (2, mci%n_dim))
      region(1,:) = 0
      region(2,:) = 1
      num_div = current%n_bins
      instance%n_adapt_grids = 0
      instance%n_adapt_weights = 0
      if (.not. instance%grids_defined) then
         call vamp_create_grids (instance%grids, &
              region, &
              current%n_calls, &
              weights = instance%w, &
              num_div = num_div, &
              stratified = mci%grid_par%stratified)
         instance%grids_defined = .true.
      else
         call msg_bug ("VAMP: create grids: grids already defined")
      end if
    end associate
  end subroutine mci_vamp_instance_create_grids

  module subroutine mci_vamp_instance_discard_integrals (instance, reshape)
    class(mci_vamp_instance_t), intent(inout) :: instance
    logical, intent(in) :: reshape
    instance%calls = 0
    instance%calls_valid = 0
    instance%integral = 0
    instance%error = 0
    instance%efficiency = 0
    associate (mci => instance%mci)
      if (instance%grids_defined) then
         if (mci%grid_par%use_vamp_equivalences) then
            if (reshape) then
               call vamp_discard_integrals (instance%grids, &
                    num_calls = instance%n_calls, &
                    stratified = mci%grid_par%stratified, &
                    eq = mci%equivalences)
            else
               call vamp_discard_integrals (instance%grids, &
                    stratified = mci%grid_par%stratified, &
                    eq = mci%equivalences)
            end if
         else
            if (reshape) then
               call vamp_discard_integrals (instance%grids, &
                    num_calls = instance%n_calls, &
                    stratified = mci%grid_par%stratified)
            else
               call vamp_discard_integrals (instance%grids, &
                    stratified = mci%grid_par%stratified)
            end if
         end if
      else
         call msg_bug ("VAMP: discard integrals: grids undefined")
      end if
    end associate
  end subroutine mci_vamp_instance_discard_integrals

  module subroutine mci_vamp_instance_allow_adaptation (instance)
    class(mci_vamp_instance_t), intent(inout) :: instance
    instance%allow_adapt_grids = .true.
    instance%allow_adapt_weights = .true.
  end subroutine mci_vamp_instance_allow_adaptation

  module subroutine mci_vamp_instance_adapt_grids (instance)
    class(mci_vamp_instance_t), intent(inout) :: instance
    if (instance%enable_adapt_grids .and. instance%allow_adapt_grids) then
       if (instance%grids_defined) then
          call vamp_refine_grids (instance%grids)
          instance%n_adapt_grids = instance%n_adapt_grids + 1
      else
         call msg_bug ("VAMP: adapt grids: grids undefined")
      end if
    end if
  end subroutine mci_vamp_instance_adapt_grids

  module function mci_vamp_instance_get_efficiency_array &
       (mci) result (efficiency)
    class(mci_vamp_instance_t), intent(in) :: mci
    real(default), dimension(:), allocatable :: efficiency
    allocate (efficiency (mci%mci%n_channel))
    if (.not. mci%negative_weights) then
       where (mci%grids%grids%f_max /= 0)
          efficiency = mci%grids%grids%mu(1) / abs (mci%grids%grids%f_max)
       elsewhere
          efficiency = 0
       end where
    else
       where (mci%grids%grids%f_max /= 0)
          efficiency = &
               (mci%grids%grids%mu_plus(1) - mci%grids%grids%mu_minus(1)) &
               / abs (mci%grids%grids%f_max)
       elsewhere
          efficiency = 0
       end where
    end if
  end function mci_vamp_instance_get_efficiency_array

  module function mci_vamp_instance_get_efficiency (mci) result (efficiency)
    class(mci_vamp_instance_t), intent(in) :: mci
    real(default) :: efficiency
    real(default), dimension(:), allocatable :: weight
    real(default) :: norm
    allocate (weight (mci%mci%n_channel))
    weight = mci%grids%weights * abs (mci%grids%grids%f_max)
    norm = sum (weight)
    if (norm /= 0) then
       efficiency = dot_product (mci%get_efficiency_array (), weight) / norm
    else
       efficiency = 1
    end if
  end function mci_vamp_instance_get_efficiency

  module subroutine mci_vamp_instance_init_simulation &
       (instance, safety_factor)
    class(mci_vamp_instance_t), intent(inout) :: instance
    real(default), intent(in), optional :: safety_factor
    associate (mci => instance%mci)
      allocate (instance%vamp_x (mci%n_dim))
      instance%it = 0
      instance%calls = 0
      instance%generating_events = .true.
      if (present (safety_factor))  instance%safety_factor = safety_factor
      if (.not. instance%grids_defined) then
         if (mci%grid_filename_set) then
            if (.not. mci%check_grid_file) &
                 call msg_warning ("Reading grid file: MD5 sum check disabled")
            call msg_message ("Simulate: " &
                 // "using integration grids from file '" &
                 // char (mci%grid_filename) // "'")
            call mci%read_grids_data (instance)
            if (instance%safety_factor /= 1) then
               write (msg_buffer, "(A,ES10.3,A)")  "Simulate: &
                    &applying safety factor", instance%safety_factor, &
                    " to event rejection"
               call msg_message ()
               instance%grids%grids%f_max = &
                    instance%grids%grids%f_max * instance%safety_factor
            end if
         else
            call msg_bug ("VAMP: simulation: no grids, no grid filename")
         end if
      end if
    end associate
  end subroutine mci_vamp_instance_init_simulation

  module subroutine mci_vamp_instance_final_simulation (instance)
    class(mci_vamp_instance_t), intent(inout) :: instance
    if (allocated (instance%vamp_x))  deallocate (instance%vamp_x)
  end subroutine mci_vamp_instance_final_simulation

  module subroutine mci_vamp_instance_compute_weight (mci, c)
    class(mci_vamp_instance_t), intent(inout) :: mci
    integer, intent(in) :: c
    integer :: i
    mci%selected_channel = c
    !$OMP PARALLEL PRIVATE(i) SHARED(mci)
    !$OMP DO
    do i = 1, mci%mci%n_channel
       if (mci%w(i) /= 0) then
          mci%gi(i) = vamp_probability (mci%grids%grids(i), mci%x(:,i))
       else
          mci%gi(i) = 0
       end if
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    mci%g = 0
    if (mci%gi(c) /= 0) then
       do i = 1, mci%mci%n_channel
          if (mci%w(i) /= 0 .and. mci%f(i) /= 0) then
             mci%g = mci%g + mci%w(i) * mci%gi(i) / mci%f(i)
          end if
       end do
    end if
    if (mci%g /= 0) then
       mci%mci_weight = mci%gi(c) / mci%g
    else
       mci%mci_weight = 0
    end if
  end subroutine mci_vamp_instance_compute_weight

  module subroutine mci_vamp_instance_record_integrand (mci, integrand)
    class(mci_vamp_instance_t), intent(inout) :: mci
    real(default), intent(in) :: integrand
    mci%integrand = integrand
  end subroutine mci_vamp_instance_record_integrand

  module function mci_vamp_instance_get_event_weight (mci) result (value)
    class(mci_vamp_instance_t), intent(in) :: mci
    real(default) :: value
    if (mci%vamp_weight_set) then
       value = mci%vamp_weight
    else
       call msg_bug ("VAMP: attempt to read undefined event weight")
    end if
  end function mci_vamp_instance_get_event_weight

  module function mci_vamp_instance_get_event_excess (mci) result (value)
    class(mci_vamp_instance_t), intent(in) :: mci
    real(default) :: value
    if (mci%vamp_weight_set) then
       value = mci%vamp_excess
    else
       call msg_bug ("VAMP: attempt to read undefined event excess weight")
    end if
  end function mci_vamp_instance_get_event_excess


end submodule mci_vamp_s

