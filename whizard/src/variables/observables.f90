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

module observables

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use subevents
  use variables

  implicit none
  private

  public :: var_list_init_num_id
  public :: var_list_init_process_results
  public :: var_list_set_observables_unary
  public :: var_list_set_observables_binary
  public :: var_list_set_observables_sev
  public :: var_list_check_observable
  public :: var_list_check_result_var

  interface
    module subroutine var_list_init_num_id (var_list, proc_id, num_id)
      type(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: proc_id
      integer, intent(in), optional :: num_id
    end subroutine var_list_init_num_id
    module subroutine var_list_init_process_results (var_list, proc_id, &
         n_calls, integral, error, accuracy, chi2, efficiency)
      type(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: proc_id
      integer, intent(in), optional :: n_calls
      real(default), intent(in), optional :: integral, error, accuracy
      real(default), intent(in), optional :: chi2, efficiency
    end subroutine var_list_init_process_results
    module subroutine var_list_set_observables_unary (var_list, prt1)
      type(var_list_t), intent(inout) :: var_list
      type(prt_t), intent(in), target :: prt1
    end subroutine var_list_set_observables_unary
    module subroutine var_list_set_observables_binary (var_list, prt1, prt2)
      type(var_list_t), intent(inout) :: var_list
      type(prt_t), intent(in), target :: prt1
      type(prt_t), intent(in), optional, target :: prt2
    end subroutine var_list_set_observables_binary
    module subroutine var_list_set_observables_sev (var_list, pval)
      type(var_list_t), intent(inout) :: var_list
      type(subevt_t), intent(in), target:: pval
    end subroutine var_list_set_observables_sev
    module subroutine var_list_check_observable (var_list, name, type)
      class(var_list_t), intent(in), target :: var_list
      type(string_t), intent(in) :: name
      integer, intent(inout) :: type
    end subroutine var_list_check_observable
    module subroutine var_list_check_result_var (var_list, name, type)
      class(var_list_t), intent(in), target :: var_list
      type(string_t), intent(in) :: name
      integer, intent(inout) :: type
    end subroutine var_list_check_result_var
  end interface

end module observables
