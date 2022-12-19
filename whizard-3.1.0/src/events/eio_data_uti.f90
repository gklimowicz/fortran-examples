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

module eio_data_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use event_base

  use eio_data

  implicit none
  private

  public :: eio_data_1
  public :: eio_data_2

contains

  subroutine eio_data_1 (u)
    integer, intent(in) :: u
    type(event_sample_data_t) :: data

    write (u, "(A)")  "* Test output: eio_data_1"
    write (u, "(A)")  "*   Purpose:  display event sample data"
    write (u, "(A)")

    write (u, "(A)")  "* Decay process, one component"
    write (u, "(A)")

    call data%init (1, 1)
    data%n_beam = 1
    data%pdg_beam(1) = 25
    data%energy_beam(1) = 125

    data%norm_mode = NORM_UNIT

    data%proc_num_id = [42]
    data%cross_section = [1.23e-4_default]
    data%error = 5e-6_default

    data%md5sum_prc = "abcdefghijklmnopabcdefghijklmnop"
    data%md5sum_cfg = "12345678901234561234567890123456"
    data%md5sum_alt(1) = "uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu"

    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Scattering process, two components"
    write (u, "(A)")

    call data%init (2)
    data%n_beam = 2
    data%pdg_beam = [2212, -2212]
    data%energy_beam = [8._default, 10._default]

    data%norm_mode = NORM_SIGMA

    data%proc_num_id = [12, 34]
    data%cross_section = [100._default, 88._default]
    data%error = [1._default, 0.1_default]

    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_data_1"

  end subroutine eio_data_1

  subroutine eio_data_2 (u)
    integer, intent(in) :: u
    type(string_t) :: s
    logical :: unweighted
    real(default) :: w, w0, sigma
    integer :: n

    write (u, "(A)")  "* Test output: eio_data_2"
    write (u, "(A)")  "*   Purpose:  handle event normalization"
    write (u, "(A)")

    write (u, "(A)")  "* Normalization strings"
    write (u, "(A)")

    s = "auto"
    unweighted = .true.
    write (u, "(1x,A,1x,L1,1x,A)")  char (s), unweighted, &
         char (event_normalization_string &
         (event_normalization_mode (s, unweighted)))
    s = "AUTO"
    unweighted = .false.
    write (u, "(1x,A,1x,L1,1x,A)")  char (s), unweighted, &
         char (event_normalization_string &
         (event_normalization_mode (s, unweighted)))

    unweighted = .true.

    s = "1"
    write (u, "(2(1x,A))") char (s), char (event_normalization_string &
         (event_normalization_mode (s, unweighted)))
    s = "1/n"
    write (u, "(2(1x,A))") char (s), char (event_normalization_string &
         (event_normalization_mode (s, unweighted)))
    s = "Sigma"
    write (u, "(2(1x,A))") char (s), char (event_normalization_string &
         (event_normalization_mode (s, unweighted)))
    s = "sigma/N"
    write (u, "(2(1x,A))") char (s), char (event_normalization_string &
         (event_normalization_mode (s, unweighted)))

    write (u, "(A)")
    write (u, "(A)")  "* Normalization update"
    write (u, "(A)")

    sigma = 5
    n = 2

    w0 = 1

    w = w0
    call event_normalization_update (w, sigma, n, NORM_UNIT, NORM_UNIT)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_N_EVT, NORM_UNIT)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_SIGMA, NORM_UNIT)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_S_N, NORM_UNIT)
    write (u, "(2(F6.3))")  w0, w

    write (u, *)

    w0 = 0.5

    w = w0
    call event_normalization_update (w, sigma, n, NORM_UNIT, NORM_N_EVT)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_N_EVT, NORM_N_EVT)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_SIGMA, NORM_N_EVT)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_S_N, NORM_N_EVT)
    write (u, "(2(F6.3))")  w0, w

    write (u, *)

    w0 = 5.0

    w = w0
    call event_normalization_update (w, sigma, n, NORM_UNIT, NORM_SIGMA)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_N_EVT, NORM_SIGMA)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_SIGMA, NORM_SIGMA)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_S_N, NORM_SIGMA)
    write (u, "(2(F6.3))")  w0, w

    write (u, *)

    w0 = 2.5

    w = w0
    call event_normalization_update (w, sigma, n, NORM_UNIT, NORM_S_N)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_N_EVT, NORM_S_N)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_SIGMA, NORM_S_N)
    write (u, "(2(F6.3))")  w0, w
    w = w0
    call event_normalization_update (w, sigma, n, NORM_S_N, NORM_S_N)
    write (u, "(2(F6.3))")  w0, w

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: eio_data_2"

  end subroutine eio_data_2


end module eio_data_uti
