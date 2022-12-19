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

module hepmc_interface_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use lorentz
  use flavors
  use colors
  use polarizations

  use hepmc_interface

  implicit none
  private

  public :: hepmc_interface_1

contains

  subroutine hepmc_interface_1 (u)
    use physics_defs, only: VECTOR
    use model_data, only: field_data_t
    integer, intent(in) :: u
    integer :: u_file, iostat
    type(hepmc_event_t) :: evt
    type(hepmc_vertex_t) :: v1, v2, v3, v4
    type(hepmc_particle_t) :: prt1, prt2, prt3, prt4, prt5, prt6, prt7, prt8
    type(hepmc_iostream_t) :: iostream
    type(flavor_t) :: flv
    type(color_t) :: col
    type(polarization_t) :: pol
    type(field_data_t), target :: photon_data
    character(80) :: buffer

    write (u, "(A)")  "* Test output: HepMC interface"
    write (u, "(A)")  "*   Purpose: test HepMC interface"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization"
    write (u, "(A)")

    ! Initialize a photon flavor object and some polarization
    call photon_data%init (var_str ("PHOTON"), 22)
    call photon_data%set (spin_type=VECTOR)
    call photon_data%freeze ()
    call flv%init (photon_data)
    call pol%init_angles &
         (flv, 0.6_default, 1._default, 0.5_default)

    ! Event initialization
    call hepmc_event_init (evt, 20, 1)

    write (u, "(A)")  "* p -> q splitting"
    write (u, "(A)")

    ! $p\to q$ splittings
    call hepmc_vertex_init (v1)
    call hepmc_event_add_vertex (evt, v1)
    call hepmc_vertex_init (v2)
    call hepmc_event_add_vertex (evt, v2)
    call particle_init (prt1, &
         0._default, 0._default, 7000._default, 7000._default, &
         2212, 3)
    call hepmc_vertex_add_particle_in (v1, prt1)
    call particle_init (prt2, &
         0._default, 0._default,-7000._default, 7000._default, &
         2212, 3)
    call hepmc_vertex_add_particle_in (v2, prt2)
    call particle_init (prt3, &
         .750_default, -1.569_default, 32.191_default, 32.238_default, &
         1, 3)
    call color_init_from_array (col, [501])
    call hepmc_particle_set_color (prt3, col)
    call hepmc_vertex_add_particle_out (v1, prt3)
    call particle_init (prt4, &
         -3.047_default, -19._default, -54.629_default, 57.920_default, &
         -2, 3)
    call color_init_from_array (col, [-501])
    call hepmc_particle_set_color (prt4, col)
    call hepmc_vertex_add_particle_out (v2, prt4)

    write (u, "(A)")  "* Hard interaction"
    write (u, "(A)")

    ! Hard interaction
    call hepmc_vertex_init (v3)
    call hepmc_event_add_vertex (evt, v3)
    call hepmc_vertex_add_particle_in (v3, prt3)
    call hepmc_vertex_add_particle_in (v3, prt4)
    call particle_init (prt6, &
         -3.813_default, 0.113_default, -1.833_default, 4.233_default, &
         22, 1)
    call hepmc_particle_set_polarization (prt6, pol)
    call hepmc_vertex_add_particle_out (v3, prt6)
    call particle_init (prt5, &
         1.517_default, -20.68_default, -20.605_default, 85.925_default, &
         -24, 3)
    call hepmc_vertex_add_particle_out (v3, prt5)
    call hepmc_event_set_signal_process_vertex (evt, v3)

    ! $W^-$ decay
    call vertex_init_pos (v4, &
         0.12_default, -0.3_default, 0.05_default, 0.004_default)
    call hepmc_event_add_vertex (evt, v4)
    call hepmc_vertex_add_particle_in (v4, prt5)
    call particle_init (prt7, &
         -2.445_default, 28.816_default, 6.082_default, 29.552_default, &
         1, 1)
    call hepmc_vertex_add_particle_out (v4, prt7)
    call particle_init (prt8, &
         3.962_default, -49.498_default, -26.687_default, 56.373_default, &
         -2, 1)
    call hepmc_vertex_add_particle_out (v4, prt8)

    ! Event output
    call hepmc_event_print (evt)
    write (u, "(A)")  "Writing to file 'hepmc_test.hepmc'"
    write (u, "(A)")

    call hepmc_iostream_open_out (iostream , var_str ("hepmc_test.hepmc"), 2)
    call hepmc_iostream_write_event (iostream, evt)
    call hepmc_iostream_close (iostream)

    write (u, "(A)")  "Writing completed"

    write (u, "(A)")
    write (u, "(A)")  "* File contents:"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = "hepmc_test.hepmc", &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (buffer(1:14) == "HepMC::Version")  buffer = "[...]"
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    ! Wrapup
    ! call pol%final ()
    call hepmc_event_final (evt)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: hepmc_interface_1"

  contains

    subroutine vertex_init_pos (v, x, y, z, t)
      type(hepmc_vertex_t), intent(out) :: v
      real(default), intent(in) :: x, y, z, t
      type(vector4_t) :: xx
      xx = vector4_moving (t, vector3_moving ([x, y, z]))
      call hepmc_vertex_init (v, xx)
    end subroutine vertex_init_pos

    subroutine particle_init (prt, px, py, pz, E, pdg, status)
      type(hepmc_particle_t), intent(out) :: prt
      real(default), intent(in) :: px, py, pz, E
      integer, intent(in) :: pdg, status
      type(vector4_t) :: p
      p = vector4_moving (E, vector3_moving ([px, py, pz]))
      call hepmc_particle_init (prt, p, pdg, status)
    end subroutine particle_init

  end subroutine hepmc_interface_1


end module hepmc_interface_uti
