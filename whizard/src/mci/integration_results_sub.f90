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

submodule (integration_results) integration_results_s

  use io_units
  use format_utils, only: mp_format, pac_fmt
  use format_defs, only: FMT_10, FMT_14
  use numeric_utils, only: pacify
  use diagnostics
  use md5

  implicit none

contains

  module function integration_entry_init (process_type, pass,&
       & it, n_it, n_calls, n_calls_valid, improved, integral, error,&
       & efficiency, efficiency_pos, efficiency_neg, chi2, chain_weights)&
       & result (entry)
    type(integration_entry_t) :: entry
    integer, intent(in) :: process_type, pass, it, n_it, n_calls, n_calls_valid
    logical, intent(in) :: improved
    real(default), intent(in) :: integral, error, efficiency, &
         efficiency_pos, efficiency_neg
    real(default), intent(in), optional :: chi2
    real(default), dimension(:), intent(in), optional :: chain_weights
    entry%process_type = process_type
    entry%pass = pass
    entry%it = it
    entry%n_it = n_it
    entry%n_calls = n_calls
    entry%n_calls_valid = n_calls_valid
    entry%improved = improved
    entry%integral = integral
    entry%error = error
    entry%efficiency = efficiency
    entry%efficiency_pos = efficiency_pos
    entry%efficiency_neg = efficiency_neg
    if (present (chi2)) entry%chi2 = chi2
    if (present (chain_weights)) then
       allocate (entry%chain_weights (size (chain_weights)))
       entry%chain_weights = chain_weights
    end if
  end function integration_entry_init

  elemental module function integration_entry_get_pass (entry) result (n)
    integer :: n
    class(integration_entry_t), intent(in) :: entry
    n = entry%pass
  end function integration_entry_get_pass

  elemental module function integration_entry_get_n_calls (entry) result (n)
    integer :: n
    class(integration_entry_t), intent(in) :: entry
    n = entry%n_calls
  end function integration_entry_get_n_calls

  elemental module function integration_entry_get_n_calls_valid &
       (entry) result (n)
    integer :: n
    class(integration_entry_t), intent(in) :: entry
    n = entry%n_calls_valid
  end function integration_entry_get_n_calls_valid

  elemental module function integration_entry_get_integral (entry) result (int)
    real(default) :: int
    class(integration_entry_t), intent(in) :: entry
    int = entry%integral
  end function integration_entry_get_integral

  elemental module function integration_entry_get_error (entry) result (err)
    real(default) :: err
    class(integration_entry_t), intent(in) :: entry
    err = entry%error
  end function integration_entry_get_error

  elemental module function integration_entry_get_relative_error &
       (entry) result (err)
    real(default) :: err
    class(integration_entry_t), intent(in) :: entry
    err = 0
    if (entry%integral /= 0) then
       err = entry%error / entry%integral
    end if
  end function integration_entry_get_relative_error

  elemental module function integration_entry_get_accuracy &
       (entry) result (acc)
    real(default) :: acc
    class(integration_entry_t), intent(in) :: entry
    acc = accuracy (entry%integral, entry%error, entry%n_calls)
  end function integration_entry_get_accuracy

  elemental module function accuracy (integral, error, n_calls) result (acc)
    real(default) :: acc
    real(default), intent(in) :: integral, error
    integer, intent(in) :: n_calls
    acc = 0
    if (integral /= 0) then
       acc = error / integral * sqrt (real (n_calls, default))
    end if
  end function accuracy

  elemental module function integration_entry_get_efficiency &
       (entry) result (eff)
    real(default) :: eff
    class(integration_entry_t), intent(in) :: entry
    eff = entry%efficiency
  end function integration_entry_get_efficiency

  elemental module function integration_entry_get_efficiency_pos &
       (entry) result (eff)
    real(default) :: eff
    class(integration_entry_t), intent(in) :: entry
    eff = entry%efficiency_pos
  end function integration_entry_get_efficiency_pos

  elemental module function integration_entry_get_efficiency_neg &
       (entry) result (eff)
    real(default) :: eff
    class(integration_entry_t), intent(in) :: entry
    eff = entry%efficiency_neg
  end function integration_entry_get_efficiency_neg

  elemental module function integration_entry_get_chi2 (entry) result (chi2)
    real(default) :: chi2
    class(integration_entry_t), intent(in) :: entry
    chi2 = entry%chi2
  end function integration_entry_get_chi2

  elemental module function integration_entry_has_improved &
       (entry) result (flag)
    logical :: flag
    class(integration_entry_t), intent(in) :: entry
    flag = entry%improved
  end function integration_entry_has_improved

  elemental module function integration_entry_get_n_groves &
       (entry) result (n_groves)
    integer :: n_groves
    class(integration_entry_t), intent(in) :: entry
    n_groves = 0
    if (allocated (entry%chain_weights)) then
       n_groves = size (entry%chain_weights, 1)
    end if
  end function integration_entry_get_n_groves

  module subroutine integration_entry_write (entry, unit, verbosity, suppress)
    class(integration_entry_t), intent(in) :: entry
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: verbosity
    logical, intent(in), optional :: suppress
    integer :: u
    character(1) :: star
    character(12) :: fmt
    character(7) :: fmt2
    character(120) :: buffer
    integer :: verb
    logical :: supp
    u = given_output_unit (unit);  if (u < 0)  return
    verb = 0; if (present (verbosity)) verb = verbosity
    supp = .false.;  if (present (suppress)) supp = suppress
    if (entry%process_type /= PRC_UNKNOWN) then
       if (entry%improved .and. .not. supp) then
          star = "*"
       else
          star = " "
       end if
       call pac_fmt (fmt, FMT_14, "3x," // FMT_10 // ",1x", suppress)
       call pac_fmt (fmt2, "1x,F6.2", "2x,F5.1", suppress)
       write (buffer, "(1x,I3,1x,I10)") entry%it, entry%n_calls
       if (verb > 1) then
          write (buffer, "(A,1x,I10)") trim (buffer), entry%n_calls_valid
       end if
       write (buffer, "(A,1x," // fmt // ",1x,ES9.2,1x,F7.2," // &
            "1x,F7.2,A1," // fmt2 // ")") &
            trim (buffer), &
            entry%integral, &
            abs(entry%error), &
            abs(integration_entry_get_relative_error (entry)) * 100, &
            abs(integration_entry_get_accuracy (entry)), &
            star, &
            entry%efficiency * 100
       if (verb > 2) then
          write (buffer, "(A,1X," // fmt2 // ",1X," // fmt2 // ")") &
               trim (buffer), &
               entry%efficiency_pos * 100, &
               entry%efficiency_neg * 100
       end if
       if (entry%n_it /= 1) then
          write (buffer, "(A,1x,F7.2,1x,I3)") &
               trim (buffer), &
               entry%chi2, &
               entry%n_it
       end if
       write (u, "(A)") trim (buffer)
    end if
    flush (u)
  end subroutine integration_entry_write

  module subroutine integration_entry_write_verbose (entry, unit)
    class(integration_entry_t), intent(in) :: entry
    integer, intent(in) :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, *)  "  process_type = ", entry%process_type
    write (u, *)  "          pass = ", entry%pass
    write (u, *)  "            it = ", entry%it
    write (u, *)  "          n_it = ", entry%n_it
    write (u, *)  "       n_calls = ", entry%n_calls
    write (u, *)  " n_calls_valid = ", entry%n_calls_valid
    write (u, *)  "      improved = ", entry%improved
    write (u, *)  "      integral = ", entry%integral
    write (u, *)  "         error = ", entry%error
    write (u, *)  "    efficiency = ", entry%efficiency
    write (u, *)  "efficiency_pos = ", entry%efficiency_pos
    write (u, *)  "efficiency_neg = ", entry%efficiency_neg
    write (u, *)  "          chi2 = ", entry%chi2
    if (allocated (entry%chain_weights)) then
       write (u, *)  "    n_groves = ", size (entry%chain_weights)
       write (u, *)  "chain_weights = ", entry%chain_weights
    else
       write (u, *)  "    n_groves = 0"
    end if
    flush (u)
  end subroutine integration_entry_write_verbose

  module subroutine integration_entry_read (entry, unit)
    class(integration_entry_t), intent(out) :: entry
    integer, intent(in) :: unit
    character(30) :: dummy
    character :: equals
    integer :: n_groves
    read (unit, *)  dummy, equals, entry%process_type
    read (unit, *)  dummy, equals, entry%pass
    read (unit, *)  dummy, equals, entry%it
    read (unit, *)  dummy, equals, entry%n_it
    read (unit, *)  dummy, equals, entry%n_calls
    read (unit, *)  dummy, equals, entry%n_calls_valid
    read (unit, *)  dummy, equals, entry%improved
    read (unit, *)  dummy, equals, entry%integral
    read (unit, *)  dummy, equals, entry%error
    read (unit, *)  dummy, equals, entry%efficiency
    read (unit, *)  dummy, equals, entry%efficiency_pos
    read (unit, *)  dummy, equals, entry%efficiency_neg
    read (unit, *)  dummy, equals, entry%chi2
    read (unit, *)  dummy, equals, n_groves
    if (n_groves /= 0) then
       allocate (entry%chain_weights (n_groves))
       read (unit, *)  dummy, equals, entry%chain_weights
    end if
  end subroutine integration_entry_read

  module subroutine integration_entry_write_chain_weights (entry, unit)
    class(integration_entry_t), intent(in) :: entry
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    if (allocated (entry%chain_weights)) then
       do i = 1, size (entry%chain_weights)
          write (u, "(1x,I3)", advance="no")  &
               nint (entry%chain_weights(i) * 100)
       end do
       write (u, *)
    end if
  end subroutine integration_entry_write_chain_weights

  module subroutine integration_results_init (results, process_type)
    class(integration_results_t), intent(out) :: results
    integer, intent(in) :: process_type
    results%process_type = process_type
    results%n_pass = 0
    results%n_it = 0
    allocate (results%entry (RESULTS_CHUNK_SIZE))
    allocate (results%average (RESULTS_CHUNK_SIZE))
  end subroutine integration_results_init

  module subroutine integration_results_set_verbosity (results, verbosity)
    class(integration_results_t), intent(inout) :: results
    integer, intent(in) :: verbosity
    results%verbosity = verbosity
  end subroutine integration_results_set_verbosity

  module subroutine integration_results_set_error_threshold &
       (results, error_threshold)
    class(integration_results_t), intent(inout) :: results
    real(default), intent(in) :: error_threshold
    results%error_threshold = error_threshold
  end subroutine integration_results_set_error_threshold

  module subroutine integration_results_write (object, unit, suppress)
    class(integration_results_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: suppress
    logical :: verb
    integer :: u, n
    u = given_output_unit (unit);  if (u < 0)  return
    call object%write_dline (unit)
    if (object%n_it /= 0) then
       call object%write_header (unit, logfile = .false.)
       call object%write_dline (unit)
       do n = 1, object%n_it
          if (n > 1) then
             if (object%entry(n)%pass /= object%entry(n-1)%pass) then
                call object%write_hline (unit)
                call object%average(object%entry(n-1)%pass)%write ( &
                     & unit, suppress = suppress)
                call object%write_hline (unit)
             end if
          end if
          call object%entry(n)%write (unit, &
               suppress = suppress)
       end do
       call object%write_hline(unit)
       call object%average(object%n_pass)%write (unit, suppress = suppress)
    else
       call msg_message ("[WHIZARD integration results: empty]", unit)
    end if
    call object%write_dline (unit)
    flush (u)
  end subroutine integration_results_write

  module subroutine integration_results_write_verbose (object, unit)
    class(integration_results_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, n
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, *)  "begin(integration_results)"
    write (u, *)  "  n_pass = ", object%n_pass
    write (u, *)  "    n_it = ", object%n_it
    if (object%n_it > 0) then
       write (u, *)  "begin(integration_pass)"
       do n = 1, object%n_it
          if (n > 1) then
             if (object%entry(n)%pass /= object%entry(n-1)%pass) then
                write (u, *)  "end(integration_pass)"
                write (u, *)  "begin(integration_pass)"
             end if
          end if
          write (u, *)  "begin(iteration)"
          call object%entry(n)%write_verbose (unit)
          write (u, *)  "end(iteration)"
       end do
       write (u, *)  "end(integration_pass)"
    end if
    write (u, *)  "end(integration_results)"
    flush (u)
  end subroutine integration_results_write_verbose

  module subroutine integration_results_write_chain_weights (results, unit)
    class(integration_results_t), intent(in) :: results
    integer, intent(in), optional :: unit
    integer :: u, i, n
    u = given_output_unit (unit);  if (u < 0)  return
    if (allocated (results%entry(1)%chain_weights) .and. results%n_it /= 0) then
       call msg_message ("Phase-space chain (grove) weight history: " &
            // "(numbers in %)", unit)
       write (u, "(A9)", advance="no")  "| chain |"
    do i = 1, integration_entry_get_n_groves (results%entry(1))
          write (u, "(1x,I3)", advance="no")  i
       end do
       write (u, *)
       call results%write_dline (unit)
       do n = 1, results%n_it
          if (n > 1) then
             if (results%entry(n)%pass /= results%entry(n-1)%pass) then
                call results%write_hline (unit)
             end if
          end if
          write (u, "(1x,I6,1x,A1)", advance="no")  n, "|"
          call results%entry(n)%write_chain_weights (unit)
       end do
       flush (u)
       call results%write_dline(unit)
    end if
  end subroutine integration_results_write_chain_weights

  module subroutine integration_results_read (results, unit)
    class(integration_results_t), intent(out) :: results
    integer, intent(in) :: unit
    character(80) :: buffer
    character :: equals
    integer :: pass, it
    read (unit, *)  buffer
    if (trim (adjustl (buffer)) /= "begin(integration_results)") then
       call read_err ();  return
    end if
    read (unit, *)  buffer, equals, results%n_pass
    read (unit, *)  buffer, equals, results%n_it
    allocate (results%entry (results%n_it + RESULTS_CHUNK_SIZE))
    allocate (results%average (results%n_it + RESULTS_CHUNK_SIZE))
    it = 0
    do pass = 1, results%n_pass
       read (unit, *)  buffer
       if (trim (adjustl (buffer)) /= "begin(integration_pass)") then
          call read_err ();  return
       end if
       READ_ENTRIES: do
          read (unit, *)  buffer
          if (trim (adjustl (buffer)) /= "begin(iteration)") then
             exit READ_ENTRIES
          end if
          it = it + 1
          call results%entry(it)%read (unit)
          read (unit, *)  buffer
          if (trim (adjustl (buffer)) /= "end(iteration)") then
             call read_err (); return
          end if
       end do READ_ENTRIES
       if (trim (adjustl (buffer)) /= "end(integration_pass)") then
          call read_err (); return
       end if
       results%average(pass) = compute_average (results%entry, pass)
    end do
    read (unit, *)  buffer
    if (trim (adjustl (buffer)) /= "end(integration_results)") then
       call read_err (); return
    end if
  contains
    subroutine read_err ()
      call msg_fatal ("Reading integration results from file: syntax error")
    end subroutine read_err
  end subroutine integration_results_read

  module subroutine write_header (results, unit, logfile)
    class(integration_results_t), intent(in) :: results
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: logfile
    character(5) :: phys_unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    select case (results%process_type)
    case (PRC_DECAY);      phys_unit = "[GeV]"
    case (PRC_SCATTERING); phys_unit = "[fb] "
    case default
       phys_unit = "    "
    end select
    write (msg_buffer, "(A, A)") &
         "It      Calls"
    if (results%verbosity > 1) then
       write (msg_buffer, "(A, A)") trim (msg_buffer), &
            "      Valid"
    end if
    write (msg_buffer, "(A, A)") trim (msg_buffer), &
         "  Integral" // phys_unit // &
         " Error" // phys_unit // &
         "  Err[%]    Acc  Eff[%]"
    if (results%verbosity > 2) then
       write (msg_buffer, "(A, A)") trim (msg_buffer), &
            "  (+)[%]  (-)[%]"
    end if
    write (msg_buffer, "(A, A)") trim (msg_buffer), &
         "   Chi2 N[It] |"
    call msg_message (unit=u, logfile=logfile)
  end subroutine write_header

  module subroutine write_hline (results, unit)
    class(integration_results_t), intent(in) :: results
    integer, intent(in), optional :: unit
    integer :: u, len
    u = given_output_unit (unit);  if (u < 0)  return
    len = 77
    if (results%verbosity > 1) len = len + 11
    if (results%verbosity > 2) len = len + 16
    write (u, "(A)")  "|" // (repeat ("-", len)) // "|"
    flush (u)
  end subroutine write_hline

  module subroutine write_dline (results, unit)
    class(integration_results_t), intent(in) :: results
    integer, intent(in), optional :: unit
    integer :: u, len
    u = given_output_unit (unit);  if (u < 0)  return
    len = 77
    if (results%verbosity > 1) len = len + 11
    if (results%verbosity > 2) len = len + 16
    write (u, "(A)")  "|" // (repeat ("=", len)) // "|"
    flush (u)
  end subroutine write_dline

  module subroutine integration_results_display_init &
       (results, screen, unit)
    class(integration_results_t), intent(inout) :: results
    logical, intent(in) :: screen
    integer, intent(in), optional :: unit
    integer :: u
    if (present (unit))  results%unit = unit
    u = given_output_unit ()
    results%screen = screen
    if (results%n_it == 0) then
       if (results%screen) then
          call results%write_dline (u)
          call results%write_header (u, &
               logfile=.false.)
          call results%write_dline (u)
       end if
       if (results%unit /= 0) then
          call results%write_dline (results%unit)
          call results%write_header (results%unit, &
               logfile=.false.)
          call results%write_dline (results%unit)
       end if
    else
       if (results%screen) then
          call results%write_hline (u)
       end if
       if (results%unit /= 0) then
          call results%write_hline (results%unit)
       end if
    end if
  end subroutine integration_results_display_init

  module subroutine integration_results_display_current (results, pacify)
    class(integration_results_t), intent(in) :: results
    integer :: u
    logical, intent(in), optional :: pacify
    u = given_output_unit ()
    if (results%screen) then
       call results%entry(results%n_it)%write (u, &
            verbosity = results%verbosity, suppress = pacify)
    end if
    if (results%unit /= 0) then
       call results%entry(results%n_it)%write ( &
            results%unit, verbosity = results%verbosity, suppress = pacify)
    end if
  end subroutine integration_results_display_current

  module subroutine integration_results_display_pass (results, pacify)
    class(integration_results_t), intent(in) :: results
    logical, intent(in), optional :: pacify
    integer :: u
    u = given_output_unit ()
    if (results%screen) then
       call results%write_hline (u)
       call results%average(results%entry(results%n_it)%pass)%write ( &
                u, verbosity = results%verbosity, suppress = pacify)
    end if
    if (results%unit /= 0) then
       call results%write_hline (results%unit)
       call results%average(results%entry(results%n_it)%pass)%write ( &
                results%unit, verbosity = results%verbosity, suppress = pacify)
    end if
  end subroutine integration_results_display_pass

  module subroutine integration_results_display_final (results)
    class(integration_results_t), intent(inout) :: results
    integer :: u
    u = given_output_unit ()
    if (results%screen) then
       call results%write_dline (u)
    end if
    if (results%unit /= 0) then
       call results%write_dline (results%unit)
    end if
    results%screen = .false.
    results%unit = 0
  end subroutine integration_results_display_final

  module subroutine integration_results_expand (results)
    class(integration_results_t), intent(inout) :: results
    type(integration_entry_t), dimension(:), allocatable :: entry_tmp
    if (results%n_it == size (results%entry)) then
       allocate (entry_tmp (results%n_it))
       entry_tmp = results%entry
       deallocate (results%entry)
       allocate (results%entry (results%n_it + RESULTS_CHUNK_SIZE))
       results%entry(:results%n_it) = entry_tmp
       deallocate (entry_tmp)
    end if
    if (results%n_pass == size (results%average)) then
       allocate (entry_tmp (results%n_pass))
       entry_tmp = results%average
       deallocate (results%average)
       allocate (results%average (results%n_it + RESULTS_CHUNK_SIZE))
       results%average(:results%n_pass) = entry_tmp
       deallocate (entry_tmp)
    end if
  end subroutine integration_results_expand

  module subroutine integration_results_new_pass (results)
    class(integration_results_t), intent(inout) :: results
    results%current_pass = results%current_pass + 1
  end subroutine integration_results_new_pass

  module subroutine integration_results_append (results, &
       n_it, n_calls, n_calls_valid, &
       integral, error, efficiency, efficiency_pos, efficiency_neg, &
       chain_weights)
    class(integration_results_t), intent(inout) :: results
    integer, intent(in) :: n_it, n_calls, n_calls_valid
    real(default), intent(in) :: integral, error, efficiency, &
         efficiency_pos, efficiency_neg
    real(default), dimension(:), intent(in), optional :: chain_weights
    logical :: improved
    type(integration_entry_t) :: entry
    real(default) :: err_checked
    improved = .true.
    if (results%n_it /= 0) improved = abs(accuracy (integral, error, n_calls)) &
            < abs(results%entry(results%n_it)%get_accuracy ())
    err_checked = 0
    if (abs (error) >= results%error_threshold) err_checked = error
    entry = integration_entry_t ( &
         results%process_type, results%current_pass, &
         results%n_it+1, n_it, n_calls, n_calls_valid, improved, &
         integral, err_checked, efficiency, efficiency_pos, efficiency_neg, &
         chain_weights=chain_weights)
    if (results%n_it == 0) then
       results%n_it = 1
       results%n_pass = 1
    else
       call results%expand ()
       if (entry%pass /= results%entry(results%n_it)%pass) &
            results%n_pass = results%n_pass + 1
       results%n_it = results%n_it + 1
    end if
    results%entry(results%n_it) = entry
    results%average(results%n_pass) = &
         compute_average (results%entry, entry%pass)
  end subroutine integration_results_append

  module subroutine integration_results_record_simple &
       (object, n_it, n_calls, integral, error, efficiency, &
        chain_weights, suppress)
    class(integration_results_t), intent(inout) :: object
    integer, intent(in) :: n_it, n_calls
    real(default), intent(in) :: integral, error, efficiency
    real(default), dimension(:), intent(in), optional :: chain_weights
    logical, intent(in), optional :: suppress
    real(default) :: err
    err = 0._default
    if (abs (error) >= abs (integral) * INTEGRATION_ERROR_TOLERANCE) then
       err = error
    end if
    call object%append (n_it, n_calls, 0, integral, err, efficiency, &
         0._default, 0._default, chain_weights)
    call object%display_current (suppress)
  end subroutine integration_results_record_simple

  module subroutine integration_results_record_extended (object, n_it, &
       n_calls, n_calls_valid, integral, error, efficiency, efficiency_pos, &
       efficiency_neg, chain_weights, suppress)
    class(integration_results_t), intent(inout) :: object
    integer, intent(in) :: n_it, n_calls, n_calls_valid
    real(default), intent(in) :: integral, error, efficiency, &
         efficiency_pos, efficiency_neg
    real(default), dimension(:), intent(in), optional :: chain_weights
    logical, intent(in), optional :: suppress
    real(default) :: err
    err = 0._default
    if (abs (error) >= abs (integral) * INTEGRATION_ERROR_TOLERANCE) then
       err = error
    end if
    call object%append (n_it, n_calls, n_calls_valid, integral, err, &
         efficiency, efficiency_pos, efficiency_neg, chain_weights)
    call object%display_current (suppress)
  end subroutine integration_results_record_extended

  type(integration_entry_t) function compute_average (entry, pass) &
       & result (result)
    type(integration_entry_t), dimension(:), intent(in) :: entry
    integer, intent(in) :: pass
    integer :: i
    logical, dimension(size(entry)) :: mask
    real(default), dimension(size(entry)) :: ivar
    real(default) :: sum_ivar, variance
    result%process_type = entry(1)%process_type
    result%pass = pass
    mask = entry%pass == pass .and. entry%process_type /= PRC_UNKNOWN
    result%it = maxval (entry%it, mask)
    result%n_it = count (mask)
    result%n_calls = sum (entry%n_calls, mask)
    result%n_calls_valid = sum (entry%n_calls_valid, mask)
    if (.not. any (mask .and. entry%error == 0)) then
       where (mask)
          ivar = 1 / entry%error ** 2
       elsewhere
          ivar = 0
       end where
       sum_ivar = sum (ivar, mask)
       variance = 0
       if (sum_ivar /= 0) then
          variance = 1 / sum_ivar
       end if
       result%integral = sum (entry%integral * ivar, mask) * variance
       if (result%n_it > 1) then
          result%chi2 = &
               sum ((entry%integral - result%integral)**2 * ivar, mask) &
               / (result%n_it - 1)
       end if
    else if (result%n_it /= 0) then
       result%integral = sum (entry%integral, mask) / result%n_it
       variance = 0
       if (result%n_it > 1) then
          variance = &
               sum ((entry%integral - result%integral)**2, mask) &
               / (result%n_it - 1)
          if (result%integral /= 0) then
             if (abs (variance / result%integral) &
                  < 100 * epsilon (1._default)) then
                variance = 0
             end if
          end if
       end if
       result%chi2 = variance / result%n_it
    end if
    result%error = sqrt (variance)
    result%efficiency = entry(last_index (mask))%efficiency
    result%efficiency_pos = entry(last_index (mask))%efficiency_pos
    result%efficiency_neg = entry(last_index (mask))%efficiency_neg
  contains
    integer function last_index (mask) result (index)
      logical, dimension(:), intent(in) :: mask
      integer :: i
      do i = size (mask), 1, -1
         if (mask(i)) exit
      end do
      index = i
    end function last_index
  end function compute_average

  module function integration_results_exist (results) result (flag)
    logical :: flag
    class(integration_results_t), intent(in) :: results
    flag = results%n_pass > 0
  end function integration_results_exist

  module function results_get_entry (results, last, it, pass) result (entry)
    class(integration_results_t), intent(in) :: results
    type(integration_entry_t) :: entry
    logical, intent(in), optional :: last
    integer, intent(in), optional :: it, pass
    if (present (last)) then
       if (allocated (results%entry) .and. results%n_it > 0) then
          entry = results%entry(results%n_it)
       else
          call error ()
       end if
    else if (present (it)) then
       if (allocated (results%entry) .and. it > 0 .and. it <= results%n_it) then
          entry = results%entry(it)
       else
          call error ()
       end if
    else if (present (pass)) then
       if (allocated (results%average) &
            .and. pass > 0 .and. pass <= results%n_pass) then
          entry = results%average (pass)
       else
          call error ()
       end if
    else
       if (allocated (results%average) .and. results%n_pass > 0) then
          entry = results%average (results%n_pass)
       else
          call error ()
       end if
    end if
  contains
    subroutine error ()
      call msg_fatal ("Requested integration result is not available")
    end subroutine error
  end function results_get_entry

  module function integration_results_get_n_calls (results, last, it, pass) &
       result (n_calls)
    class(integration_results_t), intent(in), target :: results
    integer :: n_calls
    logical, intent(in), optional :: last
    integer, intent(in), optional :: it, pass
    type(integration_entry_t) :: entry
    entry = results%get_entry (last, it, pass)
    n_calls = entry%get_n_calls ()
  end function integration_results_get_n_calls

  module function integration_results_get_integral (results, last, it, pass) &
       result (integral)
    class(integration_results_t), intent(in), target :: results
    real(default) :: integral
    logical, intent(in), optional :: last
    integer, intent(in), optional :: it, pass
    type(integration_entry_t) :: entry
    entry = results%get_entry (last, it, pass)
    integral = entry%get_integral ()
  end function integration_results_get_integral

  module function integration_results_get_error (results, last, it, pass) &
       result (error)
    class(integration_results_t), intent(in), target :: results
    real(default) :: error
    logical, intent(in), optional :: last
    integer, intent(in), optional :: it, pass
    type(integration_entry_t) :: entry
    entry = results%get_entry (last, it, pass)
    error = entry%get_error ()
  end function integration_results_get_error

  module function integration_results_get_accuracy (results, last, it, pass) &
       result (accuracy)
    class(integration_results_t), intent(in), target :: results
    real(default) :: accuracy
    logical, intent(in), optional :: last
    integer, intent(in), optional :: it, pass
    type(integration_entry_t) :: entry
    entry = results%get_entry (last, it, pass)
    accuracy = entry%get_accuracy ()
  end function integration_results_get_accuracy

  module function integration_results_get_chi2 (results, last, it, pass) &
       result (chi2)
    class(integration_results_t), intent(in), target :: results
    real(default) :: chi2
    logical, intent(in), optional :: last
    integer, intent(in), optional :: it, pass
    type(integration_entry_t) :: entry
    entry = results%get_entry (last, it, pass)
    chi2 = entry%get_chi2 ()
  end function integration_results_get_chi2

  module function integration_results_get_efficiency (results, last, it, pass) &
       result (efficiency)
    class(integration_results_t), intent(in), target :: results
    real(default) :: efficiency
    logical, intent(in), optional :: last
    integer, intent(in), optional :: it, pass
    type(integration_entry_t) :: entry
    entry = results%get_entry (last, it, pass)
    efficiency = entry%get_efficiency ()
  end function integration_results_get_efficiency

  function integration_results_get_current_pass (results) result (pass)
    integer :: pass
    type(integration_results_t), intent(in) :: results
    pass = results%n_pass
  end function integration_results_get_current_pass

  function integration_results_get_current_it (results) result (it)
    integer :: it
    type(integration_results_t), intent(in) :: results
    it = 0
    if (allocated (results%entry)) then
       it = count (results%entry(1:results%n_it)%pass == results%n_pass)
    end if
  end function integration_results_get_current_it

  function integration_results_get_last_it (results) result (it)
    integer :: it
    type(integration_results_t), intent(in) :: results
    it = results%n_it
  end function integration_results_get_last_it

  function integration_results_get_best_it (results) result (it)
    integer :: it
    type(integration_results_t), intent(in) :: results
    integer :: i
    real(default) :: acc, acc_best
    acc_best = -1
    it = 0
    do i = 1, results%n_it
       if (results%entry(i)%pass == results%n_pass) then
          acc = integration_entry_get_accuracy (results%entry(i))
          if (acc_best < 0 .or. acc <= acc_best) then
             acc_best = acc
             it = i
          end if
       end if
    end do
  end function integration_results_get_best_it

  function integration_results_get_md5sum (results) result (md5sum_results)
    character(32) :: md5sum_results
    type(integration_results_t), intent(in) :: results
    integer :: u
    u = free_unit ()
    open (unit = u, status = "scratch", action = "readwrite")
    call results%write_verbose (u)
    rewind (u)
    md5sum_results = md5sum (u)
    close (u)
  end function integration_results_get_md5sum

  module subroutine integration_results_pacify (results, efficiency_reset)
    class(integration_results_t), intent(inout) :: results
    logical, intent(in), optional :: efficiency_reset
    integer :: i
    logical :: reset
    reset = .false.
    if (present (efficiency_reset))  reset = efficiency_reset
    if (allocated (results%entry)) then
       do i = 1, size (results%entry)
          call pacify (results%entry(i)%error, &
               results%entry(i)%integral * 1.E-9_default)
          if (reset)  results%entry(i)%efficiency = 1
       end do
    end if
    if (allocated (results%average)) then
       do i = 1, size (results%average)
          call pacify (results%average(i)%error, &
               results%average(i)%integral * 1.E-9_default)
          if (reset)  results%average(i)%efficiency = 1
       end do
    end if
  end subroutine integration_results_pacify

  module subroutine integration_results_record_correction (object, corr, err)
    class(integration_results_t), intent(inout) :: object
    real(default), intent(in) :: corr, err
    integer :: u
    u = given_output_unit ()
    if (object%screen) then
      call object%write_hline (u)
      call msg_message ("NLO Correction: [O(alpha_s+1)/O(alpha_s)]")
      write(msg_buffer,'(1X,A1,F7.2,A4,F6.2,1X,A3)') '(', corr, ' +- ', err, ') %'
      call msg_message ()
    end if
  end subroutine integration_results_record_correction

  module subroutine integration_results_write_driver &
       (results, filename, eff_reset)
    type(integration_results_t), intent(inout) :: results
    type(string_t), intent(in) :: filename
    logical, intent(in), optional :: eff_reset
    type(string_t) :: file_tex
    integer :: unit
    integer :: n, i, n_pass, pass
    integer, dimension(:), allocatable :: ipass
    real(default) :: ymin, ymax, yavg, ydif, y0, y1
    real(default), dimension(results%n_it) :: ymin_arr, ymax_arr
    logical :: reset
    file_tex = filename // ".tex"
    unit = free_unit ()
    open (unit=unit, file=char(file_tex), action="write", status="replace")
    reset = .false.; if (present (eff_reset))  reset = eff_reset
    n = results%n_it
    n_pass = results%n_pass
    allocate (ipass (results%n_pass))
    ipass(1) = 0
    pass = 2
    do i = 1, n-1
       if (integration_entry_get_pass (results%entry(i)) &
           /= integration_entry_get_pass (results%entry(i+1))) then
          ipass(pass) = i
          pass = pass + 1
       end if
    end do
    ymin_arr = integration_entry_get_integral (results%entry(:n)) &
                   - integration_entry_get_error (results%entry(:n))
    ymin = minval (ymin_arr)
    ymax_arr = integration_entry_get_integral (results%entry(:n)) &
                   + integration_entry_get_error (results%entry(:n))
    ymax = maxval (ymax_arr)
    yavg = (ymax + ymin) / 2
    ydif = (ymax - ymin)
    if (ydif * 1.5 > GML_MIN_RANGE_RATIO * yavg) then
       y0 = yavg - ydif * 0.75
       y1 = yavg + ydif * 0.75
    else
       y0 = yavg * (1 - GML_MIN_RANGE_RATIO / 2)
       y1 = yavg * (1 + GML_MIN_RANGE_RATIO / 2)
    end if
    write (unit, "(A)") "\documentclass{article}"
    write (unit, "(A)") "\usepackage{a4wide}"
    write (unit, "(A)") "\usepackage{gamelan}"
    write (unit, "(A)") "\usepackage{amsmath}"
    write (unit, "(A)") ""
    write (unit, "(A)") "\begin{document}"
    write (unit, "(A)") "\begin{gmlfile}"
    write (unit, "(A)") "\section*{Integration Results Display}"
    write (unit, "(A)") ""
    write (unit, "(A)") "Process: \verb|" // char (filename) // "|"
    write (unit, "(A)") ""
    write (unit, "(A)") "\vspace*{2\baselineskip}"
    write (unit, "(A)") "\unitlength 1mm"
    write (unit, "(A)") "\begin{gmlcode}"
    write (unit, "(A)") "  picture sym;  sym = fshape (circle scaled 1mm)();"
    write (unit, "(A)") "  color col.band;  col.band = 0.9white;"
    write (unit, "(A)") "  color col.eband;  col.eband = 0.98white;"
    write (unit, "(A)") "\end{gmlcode}"
    write (unit, "(A)") "\begin{gmlgraph*}(130,180)[history]"
    write (unit, "(A)") "  setup (linear, linear);"
    write (unit, "(A,I0,A)") "  history.n_pass = ", n_pass, ";"
    write (unit, "(A,I0,A)") "  history.n_it   = ", n, ";"
    write (unit, "(A,A,A)")  "  history.y0 = #""", char (mp_format (y0)), """;"
    write (unit, "(A,A,A)")  "  history.y1 = #""", char (mp_format (y1)), """;"
    write (unit, "(A)") &
         "  graphrange (#0.5, history.y0), (#(n+0.5), history.y1);"
    do pass = 1, n_pass
       write (unit, "(A,I0,A,I0,A)") &
            "  history.pass[", pass, "] = ", ipass(pass), ";"
       write (unit, "(A,I0,A,A,A)") &
            "  history.avg[", pass, "] = #""", &
            char (mp_format &
               (integration_entry_get_integral (results%average(pass)))), &
            """;"
       write (unit, "(A,I0,A,A,A)") &
            "  history.err[", pass, "] = #""", &
            char (mp_format &
               (integration_entry_get_error (results%average(pass)))), &
            """;"
       write (unit, "(A,I0,A,A,A)") &
            "  history.chi[", pass, "] = #""", &
            char (mp_format &
               (integration_entry_get_chi2 (results%average(pass)))), &
            """;"
    end do
    write (unit, "(A,I0,A,I0,A)") &
         "  history.pass[", n_pass + 1, "] = ", n, ";"
    write (unit, "(A)")  "  for i = 1 upto history.n_pass:"
    write (unit, "(A)")  "    if history.chi[i] greater one:"
    write (unit, "(A)")  "    fill plot ("
    write (unit, "(A)")  &
         "      (#(history.pass[i]  +.5), " &
         // "history.avg[i] minus history.err[i] times history.chi[i]),"
    write (unit, "(A)")  &
         "      (#(history.pass[i+1]+.5), " &
         // "history.avg[i] minus history.err[i] times history.chi[i]),"
    write (unit, "(A)")  &
         "      (#(history.pass[i+1]+.5), " &
         // "history.avg[i] plus history.err[i] times history.chi[i]),"
    write (unit, "(A)")  &
         "      (#(history.pass[i]  +.5), " &
         // "history.avg[i] plus history.err[i] times history.chi[i])"
    write (unit, "(A)")  "    ) withcolor col.eband fi;"
    write (unit, "(A)")  "    fill plot ("
    write (unit, "(A)")  &
         "      (#(history.pass[i]  +.5), history.avg[i] minus history.err[i]),"
    write (unit, "(A)")  &
         "      (#(history.pass[i+1]+.5), history.avg[i] minus history.err[i]),"
    write (unit, "(A)")  &
         "      (#(history.pass[i+1]+.5), history.avg[i] plus history.err[i]),"
    write (unit, "(A)")  &
         "      (#(history.pass[i]  +.5), history.avg[i] plus history.err[i])"
    write (unit, "(A)")  "    ) withcolor col.band;"
    write (unit, "(A)")  "    draw plot ("
    write (unit, "(A)")  &
         "      (#(history.pass[i]  +.5), history.avg[i]),"
    write (unit, "(A)")  &
         "      (#(history.pass[i+1]+.5), history.avg[i])"
    write (unit, "(A)")  "      ) dashed evenly;"
    write (unit, "(A)")  "  endfor"
    write (unit, "(A)")  "  for i = 1 upto history.n_pass + 1:"
    write (unit, "(A)")  "    draw plot ("
    write (unit, "(A)")  &
         "      (#(history.pass[i]+.5), history.y0),"
    write (unit, "(A)")  &
         "      (#(history.pass[i]+.5), history.y1)"
    write (unit, "(A)")  "      ) dashed withdots;"
    write (unit, "(A)")  "  endfor"
    do i = 1, n
       write (unit, "(A,I0,A,A,A,A,A)") "  plot (history) (#", &
          i, ", #""", &
    char (mp_format (integration_entry_get_integral (results%entry(i)))),&
          """) vbar #""", &
          char (mp_format (integration_entry_get_error (results%entry(i)))), &
          """;"
    end do
    write (unit, "(A)") "  draw piecewise from (history) " &
      // "withsymbol sym;"
    write (unit, "(A)") "  fullgrid.lr (5,20);"
    write (unit, "(A)") "  standardgrid.bt (n);"
    write (unit, "(A)")  "  begingmleps ""Whizard-Logo.eps"";"
    write (unit, "(A)")  "    base := (120*unitlength,170*unitlength);"
    write (unit, "(A)")  "    height := 9.6*unitlength;"
    write (unit, "(A)")  "    width := 11.2*unitlength;"
    write (unit, "(A)")  "  endgmleps;"
    write (unit, "(A)") "\end{gmlgraph*}"
    write (unit, "(A)") "\end{gmlfile}"
    write (unit, "(A)") "\clearpage"
    write (unit, "(A)") "\begin{verbatim}"
    if (reset) then
      call results%pacify (reset)
    end if
    call integration_results_write (results, unit)
    write (unit, "(A)") "\end{verbatim}"
    write (unit, "(A)") "\end{document}"
    close (unit)
  end subroutine integration_results_write_driver

  module subroutine integration_results_compile_driver &
       (results, filename, os_data)
    type(integration_results_t), intent(in) :: results
    type(string_t), intent(in) :: filename
    type(os_data_t), intent(in) :: os_data
    integer :: unit_dev, status
    type(string_t) :: file_tex, file_dvi, file_ps, file_pdf, file_mp
    type(string_t) :: setenv_tex, setenv_mp, pipe, pipe_dvi
    if (.not. os_data%event_analysis) then
       call msg_warning ("Skipping integration history display " &
           // "because latex or mpost is not available")
       return
    end if
    file_tex = filename // ".tex"
    file_dvi = filename // ".dvi"
    file_ps = filename // ".ps"
    file_pdf = filename // ".pdf"
    file_mp = filename // ".mp"
    call msg_message ("Creating integration history display "&
         // char (file_ps) // " and " // char (file_pdf))
    BLOCK: do
       unit_dev = free_unit ()
       open (file = "/dev/null", unit = unit_dev, &
              action = "write", iostat = status)
       if (status /= 0) then
          pipe = ""
          pipe_dvi = ""
       else
          pipe = " > /dev/null"
          pipe_dvi = " 2>/dev/null 1>/dev/null"
       end if
       close (unit_dev)
       if (os_data%whizard_texpath /= "") then
          setenv_tex = &
               "TEXINPUTS=" // os_data%whizard_texpath // ":$TEXINPUTS "
          setenv_mp = &
               "MPINPUTS=" // os_data%whizard_texpath // ":$MPINPUTS "
       else
          setenv_tex = ""
          setenv_mp = ""
       end if
       call os_system_call (setenv_tex // os_data%latex // " " // &
            file_tex // pipe, status)
       if (status /= 0)  exit BLOCK
       if (os_data%gml /= "") then
          call os_system_call (setenv_mp // os_data%gml // " " // &
               file_mp // pipe, status)
       else
          call msg_error ("Could not use GAMELAN/MetaPOST.")
          exit BLOCK
       end if
       if (status /= 0)  exit BLOCK
       call os_system_call (setenv_tex // os_data%latex // " " // &
             file_tex // pipe, status)
       if (status /= 0)  exit BLOCK
       if (os_data%event_analysis_ps) then
          call os_system_call (os_data%dvips // " " // &
             file_dvi // pipe_dvi, status)
          if (status /= 0)  exit BLOCK
       else
          call msg_warning ("Skipping PostScript generation because dvips " &
               // "is not available")
          exit BLOCK
       end if
       if (os_data%event_analysis_pdf) then
          call os_system_call (os_data%ps2pdf // " " // &
                  file_ps, status)
          if (status /= 0)  exit BLOCK
       else
          call msg_warning ("Skipping PDF generation because ps2pdf " &
               // "is not available")
          exit BLOCK
       end if
       exit BLOCK
    end do BLOCK
    if (status /= 0) then
       call msg_error ("Unable to compile integration history display")
    end if
  end subroutine integration_results_compile_driver


end submodule integration_results_s

