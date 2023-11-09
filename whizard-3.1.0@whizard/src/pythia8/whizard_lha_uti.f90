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

module whizard_lha_uti
  use kinds, only: default
  use io_units
  use whizard_lha
  use flavors, only: flavor_t
  use lorentz, only: vector4_at_rest
  use subevents, only: PRT_BEAM, PRT_INCOMING, PRT_OUTGOING, &
       PRT_UNDEFINED, PRT_VIRTUAL, PRT_RESONANT, PRT_BEAM_REMNANT
  use particles, only: particle_set_t

  implicit none
  private

  public :: whizard_lha_1
  public :: whizard_lha_2

contains
  subroutine whizard_lha_1 (u)
    integer, intent(in) :: u
    type(whizard_lha_t) :: lha
    integer :: i
    integer, parameter :: N_PROC = 5
    real(default), dimension(N_PROC) :: xsec, xerror, max_weight
    write (u, "(A)") "* Test output: whizard_lha_1"
    write (u, "(A)") "*   Purpose: Construct LHAupWhizard object and initialize the beams."
    write (u, *)

    xsec = [1.0, 1.2, 1.4, 1.6, 1.8] * 1e3 ! fb
    xerror = 0.05 * xsec
    max_weight = 1e-3 * xsec

    call lha%init ()

    write (u, "(A)")
    write (u, "(A)") "* Set initialisation (Beams) and weighting strategy."
    write (u, "(A)")

    call lha%set_init &
         ([2212, 2212], [6500._default, 6500._default], 1, .true., .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Set process parameters for 5 different processes."
    write (u, "(A)")
    do i = 1, N_PROC
       call lha%set_process_parameters (process_id = i, &
            cross_section = xsec(i), error = xerror(i), &
            max_weight = max_weight(i))
    end do

    call lha%list_init ()

    write (u, "(A)")
    write (u, "(A)") "* Cleanup"

    call lha%final ()
  end subroutine whizard_lha_1

  subroutine whizard_lha_2 (u)
    integer, intent(in) :: u
    type(whizard_lha_t) :: lha
    integer :: i
    integer, parameter :: N_PROC = 1
    real(default), dimension(N_PROC) :: xsec, xerror, max_weight
    type(particle_set_t) :: pset
    type(flavor_t), dimension(:), allocatable :: flv
    pset%n_beam = 2
    pset%n_in   = 2
    pset%n_vir  = 2
    pset%n_out  = 3
    pset%n_tot  = 9

    write (u, "(A)") "* Test output: whizard_lha_2"
    write (u, "(A)") "*   Purpose: Setup LHAupWhizard and set event record."
    write (u, "(A)")

    xsec = [1.0] * 1e3
    xerror = 0.05 * xsec
    max_weight = 1e-3 * xsec

    call lha%init ()

    write (u, "(A)")
    write (u, "(A)") "* Set initialisation (Beams) and weighting strategy."
    write (u, "(A)")

    call lha%set_init &
         ([2212, 2212], [6500._default, 6500._default], 1, .true., .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Set process parameters for 5 different processes."
    write (u, "(A)")
    do i = 1, N_PROC
       call lha%set_process_parameters (process_id = i, &
            cross_section = xsec(i), error = xerror(i), &
            max_weight = max_weight(i))
    end do

    write (u, "(A)")
    write (u, "(A)") "* Set event record."
    write (u, "(A)")

    allocate (pset%prt (pset%n_tot))
    call pset%prt(1)%reset_status (PRT_BEAM)
    call pset%prt(2)%reset_status (PRT_BEAM)
    call pset%prt(3)%reset_status (PRT_INCOMING)
    call pset%prt(4)%reset_status (PRT_INCOMING)
    call pset%prt(5)%reset_status (PRT_BEAM_REMNANT)
    call pset%prt(6)%reset_status (PRT_BEAM_REMNANT)
    call pset%prt(7)%reset_status (PRT_OUTGOING)
    call pset%prt(8)%reset_status (PRT_OUTGOING)
    call pset%prt(9)%reset_status (PRT_OUTGOING)

    call pset%prt(1)%set_children ([3,5])
    call pset%prt(2)%set_children ([4,6])
    call pset%prt(3)%set_children ([7,8,9])
    call pset%prt(4)%set_children ([7,8,9])

    call pset%prt(3)%set_parents ([1])
    call pset%prt(4)%set_parents ([2])
    call pset%prt(5)%set_parents ([1])
    call pset%prt(6)%set_parents ([2])
    call pset%prt(7)%set_parents ([3,4])
    call pset%prt(8)%set_parents ([3,4])
    call pset%prt(9)%set_parents ([3,4])

    call pset%prt(1)%set_momentum (vector4_at_rest (1._default))
    call pset%prt(2)%set_momentum (vector4_at_rest (2._default))
    call pset%prt(3)%set_momentum (vector4_at_rest (4._default))
    call pset%prt(4)%set_momentum (vector4_at_rest (6._default))
    call pset%prt(5)%set_momentum (vector4_at_rest (3._default))
    call pset%prt(6)%set_momentum (vector4_at_rest (5._default))
    call pset%prt(7)%set_momentum (vector4_at_rest (7._default))
    call pset%prt(8)%set_momentum (vector4_at_rest (8._default))
    call pset%prt(9)%set_momentum (vector4_at_rest (9._default))

    allocate (flv (9))
    call flv%init ([2011, 2012, 11, 12, 91, 92, 3, 4, 5])
    do i = 1, 9
       call pset%prt(i)%set_flavor (flv(i))
    end do

    call lha%set_event_process (1, 1000._default, 0.1_default, 1._default / 127., 1._default)
    call lha%set_event (1, pset)
    call lha%list_init ()

    write (u, "(A)")
    write (u, "(A)") "* Cleanup"

    call lha%final ()
  end subroutine whizard_lha_2

end module whizard_lha_uti
