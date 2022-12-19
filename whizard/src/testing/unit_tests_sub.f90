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

submodule (unit_tests) unit_tests_s

  use io_units

  implicit none

  character(*), parameter :: ref_prefix = "ref-output/"
  character(*), parameter :: ref = ".ref"

  character(*), parameter :: err_prefix = "err-output/"
  character(*), parameter :: err = ".out"


contains

  module subroutine test_results_add (list, name, description, success)
    class(test_results_t), intent(inout) :: list
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: description
    logical, intent(in) :: success
    type(test_result_t), pointer :: result
    allocate (result)
    result%success = success
    result%name = name
    result%description = description
    if (associated (list%first)) then
       list%last%next => result
    else
       list%first => result
    end if
    list%last => result
    if (success) then
       list%n_success = list%n_success + 1
    else
       list%n_failure = list%n_failure + 1
    end if
  end subroutine test_results_add

  module subroutine test_results_write (list, u)
    class(test_results_t), intent(in) :: list
    integer, intent(in) :: u
    type(test_result_t), pointer :: result
    write (u, "(A)")  "*** Test Summary ***"
    if (list%n_success > 0) then
       write (u, "(2x,A)")  "Success:"
       result => list%first
       do while (associated (result))
          if (result%success)  write (u, "(4x,A,': ',A)")  &
               char (result%name), char (result%description)
          result => result%next
       end do
    end if
    if (list%n_failure > 0) then
       write (u, "(2x,A)")  "Failure:"
       result => list%first
       do while (associated (result))
          if (.not. result%success)  write (u, "(4x,A,': ',A)")  &
               char (result%name), char (result%description)
          result => result%next
       end do
    end if
    write (u, "(A,I0)")  "Total   = ", list%n_success + list%n_failure
    write (u, "(A,I0)")  "Success = ", list%n_success
    write (u, "(A,I0)")  "Failure = ", list%n_failure
    write (u, "(A)")  "*** End of test Summary ***"
  end subroutine test_results_write

  module subroutine test_results_report (list, success)
    class(test_results_t), intent(in) :: list
    logical, intent(out) :: success
    success = list%n_failure == 0
  end subroutine test_results_report

  module function test_results_get_n_pass (list) result (n)
    class(test_results_t), intent(in) :: list
    integer :: n
    n = list%n_success
  end function test_results_get_n_pass

  module function test_results_get_n_fail (list) result (n)
    class(test_results_t), intent(in) :: list
    integer :: n
    n = list%n_failure
  end function test_results_get_n_fail

  module function test_results_get_n_total (list) result (n)
    class(test_results_t), intent(in) :: list
    integer :: n
    n = list%n_success + list%n_failure
  end function test_results_get_n_total

  module subroutine test_results_final (list)
    class(test_results_t), intent(inout) :: list
    type(test_result_t), pointer :: result
    do while (associated (list%first))
       result => list%first
       list%first => result%next
       deallocate (result)
    end do
    list%last => null ()
    list%n_success = 0
    list%n_failure = 0
  end subroutine test_results_final

  module subroutine test_results_wrapup (list, u, success)
    class(test_results_t), intent(inout) :: list
    integer, intent(in) :: u
    logical, intent(out), optional :: success
    call list%write (u)
    if (present (success))  call list%report (success)
    call list%final ()
  end subroutine test_results_wrapup

  module subroutine test (test_proc, name, description, u_log, results)
    procedure(unit_test) :: test_proc
    character(*), intent(in) :: name
    character(*), intent(in) :: description
    integer, intent(in) :: u_log
    type(test_results_t), intent(inout) :: results
    integer :: u_test
    logical :: success

    call start_test (u_log, name)

    u_test = free_unit ()
    open (unit = u_test, status="scratch", action="readwrite")
    call fatal_force_crash ()
    call test_proc (u_test)
    rewind (u_test)

    call compare_test_results (u_test, u_log, name, success)
    call results%add (name, description, success)
    close (u_test)
    
  end subroutine test

  module subroutine start_test (u_log, name)
    integer, intent(in) :: u_log
    character(*), intent(in) :: name

    write (*, "(A)", advance="no")  "Running test: " // name
    write (u_log, "(A)")  "Test: " // name
    
  end subroutine start_test
  
  module subroutine compare_test_results (u_test, u_log, name, success)
    integer, intent(in) :: u_test
    integer, intent(in) :: u_log
    character(*), intent(in) :: name
    logical, intent(out) :: success
    
    integer :: u_ref, u_err
    logical :: exist
    character(256) :: buffer1, buffer2
    integer :: iostat1, iostat2

    inquire (file=ref_prefix//name//ref, exist=exist)
    if (exist) then
       open (newunit = u_ref, file=ref_prefix//name//ref, &
            status="old", action="read")
       COMPARE_FILES: do
          read (u_test, "(A)", iostat=iostat1)  buffer1
          read (u_ref, "(A)", iostat=iostat2)  buffer2
          if (iostat1 /= iostat2) then
             success = .false.
             exit COMPARE_FILES
          else if (iostat1 < 0) then
             success = .true.
             exit COMPARE_FILES
          else if (buffer1 /= buffer2) then
             success = .false.
             exit COMPARE_FILES
          end if
       end do COMPARE_FILES
       close (u_ref)
    else
       write (*, "(A)", advance="no") " ... no reference output available"
       write (u_log, "(A)") "  No reference output available."
       success = .false.
    end if
    if (success) then
       write (*, "(A)") " ... success."
       write (u_log, "(A)")  "  Success."
    else
       write (*, "(A)") " ... failure.  See: " // err_prefix//name//err
       write (u_log, "(A)")  "  Failure."
       rewind (u_test)
       open (newunit = u_err, file=err_prefix//name//err, &
            action="write", status="replace")
       WRITE_OUTPUT: do
          read (u_test, "(A)", end=1) buffer1
          write (u_err, "(A)")  trim (buffer1)
       end do WRITE_OUTPUT
1      close (u_err)
    end if

  end subroutine compare_test_results


end submodule unit_tests_s

