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

module integration_results_uti

  use kinds, only: default

  use integration_results

  implicit none
  private

  public :: integration_results_1
  public :: integration_results_2
  public :: integration_results_3
  public :: integration_results_4
  public :: integration_results_5

contains

  subroutine integration_results_1 (u)
    integer, intent(in) :: u
    type(integration_entry_t) :: entry

    write (u, "(A)")  "* Test output: integration_results_1"
    write (u, "(A)")  "*   Purpose: record single entry and write to log"
    write (u, "(A)")

    write (u, "(A)")  "* Write single line output"
    write (u, "(A)")

    entry = integration_entry_t ( &
         & process_type = 1, &
         & pass = 1, &
         & it = 1, &
         & n_it = 10, &
         & n_calls = 1000, &
         & n_calls_valid = 500, &
         & improved = .true., &
         & integral = 1.0_default, &
         & error = 0.5_default, &
         & efficiency = 0.25_default, &
         & efficiency_pos = 0.22_default, &
         & efficiency_neg = 0.03_default)
    call entry%write (u, 3)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integration_results_1"

  end subroutine integration_results_1

  subroutine integration_results_2 (u)
    integer, intent(in) :: u
    type(integration_results_t) :: results

    write (u, "(A)")  "* Test output: integration_results_2"
    write (u, "(A)")  "*   Purpose: record single result and write to log"
    write (u, "(A)")

    write (u, "(A)")  "* Write single line output"
    write (u, "(A)")

    call results%init (PRC_DECAY)
    call results%append (1, 250, 0, 1.0_default, 0.5_default, 0.25_default,&
         & 0._default, 0._default)

    call results%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integration_results_2"

  end subroutine integration_results_2

  subroutine integration_results_3 (u)
    integer, intent(in) :: u
    type(integration_results_t) :: results

    write (u, "(A)")  "* Test output: integration_results_2"
    write (u, "(A)")  "*   Purpose: intialize display, record three entries,&
         & display pass average and finalize display"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize display and add entry"
    write (u, "(A)")

    call results%init (PRC_DECAY)
    call results%set_verbosity (1)
    call results%display_init (screen = .false., unit = u)
    call results%new_pass ()

    call results%record (1, 250, 1.0_default, 0.5_default, 0.25_default)
    call results%record (1, 250, 1.1_default, 0.5_default, 0.25_default)
    call results%record (1, 250, 0.9_default, 0.5_default, 0.25_default)

    write (u, "(A)")
    write (u, "(A)")  "* Display pass"
    write (u, "(A)")

    call results%display_pass ()

    write (u, "(A)")
    write (u, "(A)")  "* Finalize displays"
    write (u, "(A)")

    call results%display_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integration_results_3"

  end subroutine integration_results_3

  subroutine integration_results_4 (u)
    integer, intent(in) :: u
    type(integration_results_t) :: results

    write (u, "(A)")  "* Test output: integration_results_4"
    write (u, "(A)")  "*   Purpose: record extended results and display with verbosity = 2"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize display and record extended result"
    write (u, "(A)")

    call results%init (PRC_DECAY)
    call results%set_verbosity (2)
    call results%display_init (screen = .false., unit = u)
    call results%new_pass ()

    call results%record (1, 250, 150, 1.0_default, 0.5_default, 0.25_default,&
         & 0.22_default, 0.03_default)
    call results%record (1, 250, 180, 1.1_default, 0.5_default, 0.25_default,&
         & 0.23_default, 0.02_default)
    call results%record (1, 250, 130, 0.9_default, 0.5_default, 0.25_default,&
         & 0.25_default, 0.00_default)

    write (u, "(A)")
    write (u, "(A)")  "* Display pass"
    write (u, "(A)")

    call results%display_pass ()

    write (u, "(A)")
    write (u, "(A)")  "* Finalize displays"
    write (u, "(A)")

    call results%display_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integration_results_4"

  end subroutine integration_results_4

  subroutine integration_results_5 (u)
    integer, intent(in) :: u
    type(integration_results_t) :: results

    write (u, "(A)")  "* Test output: integration_results_5"
    write (u, "(A)")  "*   Purpose: record extended results and display  with verbosity = 3"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize display and record extended result"
    write (u, "(A)")

    call results%init (PRC_DECAY)
    call results%set_verbosity (3)
    call results%display_init (screen = .false., unit = u)
    call results%new_pass ()

    call results%record (1, 250, 150, 1.0_default, 0.5_default, 0.25_default,&
         & 0.22_default, 0.03_default)
    call results%record (1, 250, 180, 1.1_default, 0.5_default, 0.25_default,&
         & 0.23_default, 0.02_default)
    call results%record (1, 250, 130, 0.9_default, 0.5_default, 0.25_default,&
         & 0.25_default, 0.00_default)
    call results%display_pass ()
    call results%display_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: integration_results_5"

  end subroutine integration_results_5


end module integration_results_uti
