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

module dispatch_mci_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use variables
  use mci_base
  use mci_none
  use mci_midpoint
  use mci_vamp
  use dispatch_mci

  implicit none
  private

  public :: dispatch_mci_1

contains

  subroutine dispatch_mci_1 (u)
    integer, intent(in) :: u
    type(var_list_t) :: var_list
    class(mci_t), allocatable :: mci
    type(string_t) :: process_id

    write (u, "(A)")  "* Test output: dispatch_mci_1"
    write (u, "(A)")  "*   Purpose: select integration method"
    write (u, "(A)")

    call var_list%init_defaults (0)

    process_id = "dispatch_mci_1"

    write (u, "(A)")  "* Allocate MCI as none_t"
    write (u, "(A)")

    call var_list%set_string (&
         var_str ("$integration_method"), &
         var_str ("none"), is_known = .true.)
    call dispatch_mci_setup (mci, var_list, process_id)
    select type (mci)
    type is (mci_none_t)
       call mci%write (u)
    end select

    call mci%final ()
    deallocate (mci)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate MCI as midpoint_t"
    write (u, "(A)")

    call var_list%set_string (&
         var_str ("$integration_method"), &
         var_str ("midpoint"), is_known = .true.)
    call dispatch_mci_setup (mci, var_list, process_id)
    select type (mci)
    type is (mci_midpoint_t)
       call mci%write (u)
    end select

    call mci%final ()
    deallocate (mci)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate MCI as vamp_t"
    write (u, "(A)")

    call var_list%set_string (&
         var_str ("$integration_method"), &
         var_str ("vamp"), is_known = .true.)
    call var_list%set_int (var_str ("threshold_calls"), &
         1, is_known = .true.)
    call var_list%set_int (var_str ("min_calls_per_channel"), &
         2, is_known = .true.)
    call var_list%set_int (var_str ("min_calls_per_bin"), &
         3, is_known = .true.)
    call var_list%set_int (var_str ("min_bins"), &
         4, is_known = .true.)
    call var_list%set_int (var_str ("max_bins"), &
         5, is_known = .true.)
    call var_list%set_log (var_str ("?stratified"), &
         .false., is_known = .true.)
    call var_list%set_log (var_str ("?use_vamp_equivalences"),&
         .false., is_known = .true.)
    call var_list%set_real (var_str ("channel_weights_power"),&
         4._default, is_known = .true.)
    call var_list%set_log (&
         var_str ("?vamp_history_global_verbose"), &
         .true., is_known = .true.)
    call var_list%set_log (&
         var_str ("?vamp_history_channels"), &
         .true., is_known = .true.)
    call var_list%set_log (&
         var_str ("?vamp_history_channels_verbose"), &
         .true., is_known = .true.)
    call var_list%set_log (var_str ("?stratified"), &
         .false., is_known = .true.)

    call dispatch_mci_setup (mci, var_list, process_id)
    select type (mci)
    type is (mci_vamp_t)
       call mci%write (u)
       call mci%write_history_parameters (u)
    end select

    call mci%final ()
    deallocate (mci)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate MCI as vamp_t, allow for negative weights"
    write (u, "(A)")

    call var_list%set_string (&
         var_str ("$integration_method"), &
         var_str ("vamp"), is_known = .true.)
    call var_list%set_log (var_str ("?negative_weights"), &
         .true., is_known = .true.)

    call dispatch_mci_setup (mci, var_list, process_id)
    select type (mci)
    type is (mci_vamp_t)
       call mci%write (u)
       call mci%write_history_parameters (u)
    end select

    call mci%final ()
    deallocate (mci)

    call var_list%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_mci_1"

  end subroutine dispatch_mci_1


end module dispatch_mci_uti
