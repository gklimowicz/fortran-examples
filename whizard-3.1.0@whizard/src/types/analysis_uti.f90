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

module analysis_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use format_defs, only: FMT_19

  use analysis

  implicit none
  private

  public :: analysis_1

contains

  subroutine analysis_1 (u)
    integer, intent(in) :: u
    type(string_t) :: id1, id2, id3, id4
    integer :: i
    id1 = "foo"
    id2 = "bar"
    id3 = "hist"
    id4 = "plot"

    write (u, "(A)")  "* Test output: Analysis"
    write (u, "(A)")  "*   Purpose: test the analysis routines"
    write (u, "(A)")

    call analysis_init_observable (id1)
    call analysis_init_observable (id2)
    call analysis_init_histogram &
         (id3, 0.5_default, 5.5_default, 1._default, normalize_bins=.false.)
    call analysis_init_plot (id4)
    do i = 1, 3
       write (u, "(A,1x," // FMT_19 // ")")  "data = ", real(i,default)
       call analysis_record_data (id1, real(i,default))
       call analysis_record_data (id2, real(i,default), &
                                        weight=real(i,default))
       call analysis_record_data (id3, real(i,default))
       call analysis_record_data (id4, real(i,default), real(i,default)**2)
    end do
    write (u, "(A,10(1x,I5))") "n_entries = ", &
         analysis_get_n_entries (id1), &
         analysis_get_n_entries (id2), &
         analysis_get_n_entries (id3), &
         analysis_get_n_entries (id3, within_bounds = .true.), &
         analysis_get_n_entries (id4), &
         analysis_get_n_entries (id4, within_bounds = .true.)
    write (u, "(A,10(1x," // FMT_19 // "))")  "average   = ", &
         analysis_get_average (id1), &
         analysis_get_average (id2), &
         analysis_get_average (id3), &
         analysis_get_average (id3, within_bounds = .true.)
    write (u, "(A,10(1x," // FMT_19 // "))")  "error     = ", &
         analysis_get_error (id1), &
         analysis_get_error (id2), &
         analysis_get_error (id3), &
         analysis_get_error (id3, within_bounds = .true.)

    write (u, "(A)")
    write (u, "(A)") "* Clear analysis #2"
    write (u, "(A)")

    call analysis_clear (id2)
    do i = 4, 6
       print *, "data = ", real(i,default)
       call analysis_record_data (id1, real(i,default))
       call analysis_record_data (id2, real(i,default), &
                                        weight=real(i,default))
       call analysis_record_data (id3, real(i,default))
       call analysis_record_data (id4, real(i,default), real(i,default)**2)
    end do
    write (u, "(A,10(1x,I5))")  "n_entries = ", &
         analysis_get_n_entries (id1), &
         analysis_get_n_entries (id2), &
         analysis_get_n_entries (id3), &
         analysis_get_n_entries (id3, within_bounds = .true.), &
         analysis_get_n_entries (id4), &
         analysis_get_n_entries (id4, within_bounds = .true.)
    write (u, "(A,10(1x," // FMT_19 // "))")  "average   = ", &
         analysis_get_average (id1), &
         analysis_get_average (id2), &
         analysis_get_average (id3), &
         analysis_get_average (id3, within_bounds = .true.)
    write (u, "(A,10(1x," // FMT_19 // "))")  "error     = ", &
         analysis_get_error (id1), &
         analysis_get_error (id2), &
         analysis_get_error (id3), &
         analysis_get_error (id3, within_bounds = .true.)
    write (u, "(A)")
    call analysis_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call analysis_clear ()
    call analysis_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: analysis_1"
  end subroutine analysis_1


end module analysis_uti
