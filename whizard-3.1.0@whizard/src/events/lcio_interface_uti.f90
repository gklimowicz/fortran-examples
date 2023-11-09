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

module lcio_interface_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use lorentz
  use flavors
  use colors
  use polarizations

  use lcio_interface

  implicit none
  private

  public :: lcio_interface_1

contains

  subroutine lcio_interface_1 (u)
    use physics_defs, only: VECTOR
    use model_data, only: field_data_t
    integer, intent(in) :: u
    integer :: u_file, iostat
    type(lcio_event_t) :: evt
    type(lcio_particle_t) :: prt1, prt2, prt3, prt4, prt5, prt6, prt7, prt8
    type(flavor_t) :: flv
    type(color_t) :: col
    type(polarization_t) :: pol
    type(field_data_t), target :: photon_data
    character(220) :: buffer

    write (u, "(A)")  "* Test output: LCIO interface"
    write (u, "(A)")  "*   Purpose: test LCIO interface"
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
    call lcio_event_init (evt, 20, 1, 42)

    write (u, "(A)")  "* p -> q splitting"
    write (u, "(A)")

    ! $p\to q$ splittings
    call particle_init (prt1, &
         0._default, 0._default, 7000._default, 7000._default, &
         2212, 1._default, 3)
    call particle_init (prt2, &
         0._default, 0._default,-7000._default, 7000._default, &
         2212, 1._default, 3)
    call particle_init (prt3, &
          .750_default, -1.569_default, 32.191_default, 32.238_default, &
          1, -1._default/3._default, 3)
    call color_init_from_array (col, [501])
    call lcio_particle_set_color (prt3, col)
    call lcio_particle_set_parent (prt3, prt1)
    call lcio_particle_set_parent (prt3, prt2)
    call particle_init (prt4, &
         -3.047_default, -19._default, -54.629_default, 57.920_default, &
         -2, -2._default/3._default, 3)
    call color_init_from_array (col, [-501])
    call lcio_particle_set_color (prt4, col)
    call lcio_particle_set_parent (prt4, prt1)
    call lcio_particle_set_parent (prt4, prt2)

    write (u, "(A)")  "* Hard interaction"
    write (u, "(A)")

    ! Hard interaction
    call particle_init (prt6, &
         -3.813_default, 0.113_default, -1.833_default, 4.233_default, &
         22, 0._default, 1)
    call lcio_polarization_init (prt6, pol)
    call particle_init (prt5, &
         1.517_default, -20.68_default, -20.605_default, 85.925_default, &
         -24, -1._default, 3)
    call lcio_particle_set_parent (prt5, prt3)
    call lcio_particle_set_parent (prt5, prt4)
    call lcio_particle_set_parent (prt6, prt3)
    call lcio_particle_set_parent (prt6, prt4)

    ! $W^-$ decay
    call particle_init (prt7, &
         -2.445_default, 28.816_default, 6.082_default, 29.552_default, &
         1, -1._default/3._default, 1)
    call particle_init (prt8, &
         3.962_default, -49.498_default, -26.687_default, 56.373_default, &
         -2, -2._default/3._default, 1)
    call lcio_particle_set_t (prt7, 0.12_default)
    call lcio_particle_set_t (prt8, 0.12_default)
    call lcio_particle_set_vtx &
         (prt7, vector3_moving ([-0.3_default, 0.05_default, 0.004_default]))
    call lcio_particle_set_vtx &
         (prt8, vector3_moving ([-0.3_default, 0.05_default, 0.004_default]))
    call lcio_particle_set_parent (prt7, prt5)
    call lcio_particle_set_parent (prt8, prt5)
    call lcio_particle_add_to_evt_coll (prt1, evt)
    call lcio_particle_add_to_evt_coll (prt2, evt)
    call lcio_particle_add_to_evt_coll (prt3, evt)
    call lcio_particle_add_to_evt_coll (prt4, evt)
    call lcio_particle_add_to_evt_coll (prt5, evt)
    call lcio_particle_add_to_evt_coll (prt6, evt)
    call lcio_particle_add_to_evt_coll (prt7, evt)
    call lcio_particle_add_to_evt_coll (prt8, evt)
    call lcio_event_add_coll (evt)

    ! Event output
    write (u, "(A)")  "Writing in ASCII form to file 'lcio_test.slcio'"
    write (u, "(A)")

    call write_lcio_event (evt, var_str ("lcio_test.slcio"))

    write (u, "(A)")  "Writing completed"

    write (u, "(A)")
    write (u, "(A)")  "* File contents:"
    write (u, "(A)")

    u_file = free_unit ()
    open (u_file, file = "lcio_test.slcio", &
         action = "read", status = "old")
    do
       read (u_file, "(A)", iostat = iostat)  buffer
       if (trim (buffer) == "")  cycle
       if (buffer(1:12) == " - timestamp")  buffer = "[...]"
       if (buffer(1:6) == " date:")  buffer = "[...]"
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do
    close (u_file)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    ! Wrapup
    ! call pol%final ()
    call lcio_event_final (evt, .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: lcio_interface_1"

  contains

    subroutine particle_init &
         (prt, px, py, pz, E, pdg, charge, status)
      type(lcio_particle_t), intent(out) :: prt
      real(default), intent(in) :: px, py, pz, E, charge
      integer, intent(in) :: pdg, status
      type(vector4_t) :: p
      p = vector4_moving (E, vector3_moving ([px, py, pz]))
      call lcio_particle_init (prt, p, pdg, charge, status)
    end subroutine particle_init

  end subroutine lcio_interface_1


end module lcio_interface_uti
