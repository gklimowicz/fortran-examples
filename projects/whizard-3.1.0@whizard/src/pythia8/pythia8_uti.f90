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

module pythia8_uti
  use kinds, only: default
  use io_units
  use iso_varying_string, string_t => varying_string
  use model_data, only: model_data_t
  use particles, only: particle_t, PRT_DEFINITE_HELICITY, PRT_GENERIC_POLARIZATION
  use rng_base, only: rng_t
  use rng_stream, only: rng_stream_t
  use whizard_lha
  use pythia8

  implicit none
  private

  public :: pythia8_1
  public :: pythia8_2

contains
  subroutine pythia8_1 (u)
    integer, intent(in) :: u
    type(pythia8_t) :: pythia
    type(whizard_lha_t) :: lha

    write (u, "(A)") "* Test output: pythia8_1"
    write (u, "(A)") "*   Purpose: Construct and destruct a Pythia8 object."
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)") "* Construct Pythia8 object."
    write (u, "(A)")

    call pythia%init ()

    write (u, "(A)")
    write (u, "(A)") "* Destruct Pythia8 object."
    write (u, "(A)")

    call pythia%final ()
  end subroutine pythia8_1

  subroutine pythia8_2 (u)
    integer, intent(in) :: u
    type(pythia8_t) :: pythia
    type(whizard_lha_t) :: lha
    integer :: i
    integer, parameter :: N_PROC = 5
    real(default), dimension(N_PROC) :: xsec, xerror, max_weight

    write (u, "(A)") "* Test output: pythia8_2"
    write (u, "(A)") &
         "*   Purpose: Initialize Pythia8 with a LHA User Process object.."
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)") "* Construct Pythia8 object."
    write (u, "(A)")
    call pythia%init ()

    write (u, "(A)")
    write (u, "(A)") "* Read string 'Beam:frameType = 5' into " // &
         "Pythia8 allowing for LHA user processes."
    write (u, "(A)")

    call pythia%read_string (var_str ("Beams:frameType = 5"))
    call pythia%read_string (var_str ("Random:setSeed = on"))
    call pythia%read_string (var_str ("Random:seed = 1234"))

    write (u, "(A)")
    write (u, "(A)") &
         "* Setup LHA User Process object and let Pythia8 point to it."
    write (u, "(A)")

    call lha%init ()
    call lha%set_init &
         ([2212, 2212], [6500._default, 6500._default], 1, .false., .false.)

    xsec = [1.0, 1.2, 1.4, 1.6, 1.8] * 1e3_default ! fb
    xerror = 0.05_default * xsec
    max_weight = 1e-3_default * xsec
    do i = 1, N_PROC
       call lha%set_process_parameters (process_id = i, &
            cross_section = xsec(i), error = xerror(i), &
            max_weight = max_weight(i))
    end do
    call pythia%set_lhaup_ptr (lha)

    write (u, "(A)")
    write (u, "(A)") "* Initialize Pythia8."
    write (u, "(A)")

    call pythia%init_pythia ()

    write (u, "(A)")
    write (u, "(A)") "* Destruct Pythia8 object."
    write (u, "(A)")

    call pythia%final ()
  end subroutine pythia8_2

end module pythia8_uti
