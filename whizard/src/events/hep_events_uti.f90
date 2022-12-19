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

module hep_events_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use flavors
  use colors
  use helicities
  use quantum_numbers
  use state_matrices, only: FM_SELECT_HELICITY, FM_FACTOR_HELICITY
  use interactions
  use evaluators
  use model_data
  use particles
  use subevents
  use hepmc_interface

  use hep_events

  implicit none
  private

  public :: hep_events_1

contains

  subroutine hep_events_1 (u)
    use os_interface
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t), dimension(3) :: flv
    type(color_t), dimension(3) :: col
    type(helicity_t), dimension(3) :: hel
    type(quantum_numbers_t), dimension(3) :: qn
    type(vector4_t), dimension(3) :: p
    type(interaction_t), target :: int1, int2
    type(quantum_numbers_mask_t) :: qn_mask_conn
    type(evaluator_t), target :: eval
    type(interaction_t), pointer :: int
    type(particle_set_t) :: particle_set1, particle_set2
    type(hepmc_event_t) :: hepmc_event
    type(hepmc_iostream_t) :: iostream
    real(default) :: cross_section, error, weight
    logical :: ok

    write (u, "(A)")  "* Test output: HEP events"
    write (u, "(A)")  "*   Purpose: test HepMC event routines"
    write (u, "(A)")

    write (u, "(A)")  "* Reading model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Initializing production process"

    call int1%basic_init (2, 0, 1, set_relations=.true.)
    call flv%init ([1, -1, 23], model)
    call col%init_col_acl ([0, 0, 0], [0, 0, 0])
    call hel(3)%init ( 1, 1)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0.25_default, 0._default))
    call hel(3)%init ( 1,-1)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0._default, 0.25_default))
    call hel(3)%init (-1, 1)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0._default,-0.25_default))
    call hel(3)%init (-1,-1)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0.25_default, 0._default))
    call hel(3)%init ( 0, 0)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0.5_default, 0._default))
    call int1%freeze ()
    p(1) = vector4_moving (45._default, 45._default, 3)
    p(2) = vector4_moving (45._default,-45._default, 3)
    p(3) = p(1) + p(2)
    call int1%set_momenta (p)

    write (u, "(A)")
    write (u, "(A)")  "* Setup decay process"

    call int2%basic_init (1, 0, 2, set_relations=.true.)
    call flv%init ([23, 1, -1], model)
    call col%init_col_acl ([0, 501, 0], [0, 0, 501])
    call hel%init ([1, 1, 1], [1, 1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(1._default, 0._default))
    call hel%init ([1, 1, 1], [-1,-1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(0._default, 0.1_default))
    call hel%init ([-1,-1,-1], [1, 1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(0._default,-0.1_default))
    call hel%init ([-1,-1,-1], [-1,-1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(1._default, 0._default))
    call hel%init ([0, 1,-1], [0, 1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(4._default, 0._default))
    call hel%init ([0,-1, 1], [0, 1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(2._default, 0._default))
    call hel%init ([0, 1,-1], [0,-1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(2._default, 0._default))
    call hel%init ([0,-1, 1], [0,-1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(4._default, 0._default))
    call flv%init ([23, 2, -2], model)
    call hel%init ([0, 1,-1], [0, 1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(0.5_default, 0._default))
    call hel%init ([0,-1, 1], [0,-1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(0.5_default, 0._default))
    call int2%freeze ()
    p(2) = vector4_moving (45._default, 45._default, 2)
    p(3) = vector4_moving (45._default,-45._default, 2)
    call int2%set_momenta (p)
    call int2%set_source_link (1, int1, 3)
    call int1%basic_write (u)
    call int2%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Concatenate production and decay"

    call eval%init_product (int1, int2, qn_mask_conn, &
         connections_are_resonant=.true.)
    call eval%receive_momenta ()
    call eval%evaluate ()
    call eval%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Factorize as subevent (complete, polarized)"
    write (u, "(A)")

    int => eval%interaction_t
    call particle_set1%init &
         (ok, int, int, FM_FACTOR_HELICITY, &
          [0.2_default, 0.2_default], .false., .true.)
    call particle_set1%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Factorize as subevent (in/out only, selected helicity)"
    write (u, "(A)")

    int => eval%interaction_t
    call particle_set2%init &
         (ok, int, int, FM_SELECT_HELICITY, &
          [0.9_default, 0.9_default], .false., .false.)
    call particle_set2%write (u)
    call particle_set2%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Factorize as subevent (complete, selected helicity)"
    write (u, "(A)")

    int => eval%interaction_t
    call particle_set2%init &
         (ok, int, int, FM_SELECT_HELICITY, &
          [0.7_default, 0.7_default], .false., .true.)
    call particle_set2%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Transfer particle_set to HepMC, print, and output to"
    write (u, "(A)")  "        hep_events.hepmc.dat"
    write (u, "(A)")

    cross_section = 42.0_default
    error = 17.0_default
    weight = 1.0_default
    call hepmc_event_init (hepmc_event, 11, 127)
    call hepmc_event_from_particle_set (hepmc_event, particle_set2, &
         cross_section, error, .true.)
    call hepmc_event_add_weight (hepmc_event, weight, .true.)
    call hepmc_event_print (hepmc_event)
    call hepmc_iostream_open_out &
         (iostream , var_str ("hep_events.hepmc.dat"), 2)
    call hepmc_iostream_write_event (iostream, hepmc_event)
    call hepmc_iostream_close (iostream)

    write (u, "(A)")
    write (u, "(A)")  "* Recover from HepMC file"
    write (u, "(A)")

    call particle_set2%final ()
    call hepmc_event_final (hepmc_event)
    call hepmc_event_init (hepmc_event)
    call hepmc_iostream_open_in &
         (iostream , var_str ("hep_events.hepmc.dat"), HEPMC3_MODE_HEPMC3)
    call hepmc_iostream_read_event (iostream, hepmc_event, ok=ok)
    call hepmc_iostream_close (iostream)
    call hepmc_event_to_particle_set (particle_set2, &
         hepmc_event, model, model, PRT_DEFINITE_HELICITY)
    call particle_set2%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call particle_set1%final ()
    call particle_set2%final ()
    call eval%final ()
    call int1%final ()
    call int2%final ()
    call hepmc_event_final (hepmc_event)
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: hep_events_1"

  end subroutine hep_events_1


end module hep_events_uti
