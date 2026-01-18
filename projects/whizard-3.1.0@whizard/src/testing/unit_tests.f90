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

module unit_tests

  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: test_results_t
  public :: unit_test
  public :: test
  public :: start_test
  public :: compare_test_results

  type :: test_result_t
     logical :: success = .false.
     type(string_t) :: name
     type(string_t) :: description
     type(test_result_t), pointer :: next => null ()
  end type test_result_t

  type :: test_results_t
     private
     type(test_result_t), pointer :: first => null ()
     type(test_result_t), pointer :: last => null ()
     integer :: n_success = 0
     integer :: n_failure = 0
   contains
     procedure :: add => test_results_add
     procedure, private :: write => test_results_write
     procedure, private :: report => test_results_report
     procedure :: get_n_pass => test_results_get_n_pass
     procedure :: get_n_fail => test_results_get_n_fail
     procedure :: get_n_total => test_results_get_n_total
     procedure, private :: final => test_results_final
     procedure :: wrapup => test_results_wrapup
  end type test_results_t


  abstract interface
     subroutine unit_test (u)
       integer, intent(in) :: u
     end subroutine unit_test
  end interface


  interface
    module subroutine test_results_add (list, name, description, success)
      class(test_results_t), intent(inout) :: list
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: description
      logical, intent(in) :: success
    end subroutine test_results_add
    module subroutine test_results_write (list, u)
      class(test_results_t), intent(in) :: list
      integer, intent(in) :: u
    end subroutine test_results_write
    module subroutine test_results_report (list, success)
      class(test_results_t), intent(in) :: list
      logical, intent(out) :: success
    end subroutine test_results_report
    module function test_results_get_n_pass (list) result (n)
      class(test_results_t), intent(in) :: list
      integer :: n
    end function test_results_get_n_pass
    module function test_results_get_n_fail (list) result (n)
      class(test_results_t), intent(in) :: list
      integer :: n
    end function test_results_get_n_fail
    module function test_results_get_n_total (list) result (n)
      class(test_results_t), intent(in) :: list
      integer :: n
    end function test_results_get_n_total
    module subroutine test_results_final (list)
      class(test_results_t), intent(inout) :: list
    end subroutine test_results_final
    module subroutine test_results_wrapup (list, u, success)
      class(test_results_t), intent(inout) :: list
      integer, intent(in) :: u
      logical, intent(out), optional :: success
    end subroutine test_results_wrapup
    module subroutine test (test_proc, name, description, u_log, results)
      procedure(unit_test) :: test_proc
      character(*), intent(in) :: name
      character(*), intent(in) :: description
      integer, intent(in) :: u_log
      type(test_results_t), intent(inout) :: results
    end subroutine test
    module subroutine start_test (u_log, name)
      integer, intent(in) :: u_log
      character(*), intent(in) :: name
    end subroutine start_test
    module subroutine compare_test_results (u_test, u_log, name, success)
      integer, intent(in) :: u_test
      integer, intent(in) :: u_log
      character(*), intent(in) :: name
      logical, intent(out) :: success
    end subroutine compare_test_results
  end interface

end module unit_tests
