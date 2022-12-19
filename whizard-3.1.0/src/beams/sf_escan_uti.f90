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

module sf_escan_uti

  use kinds, only: default
  use physics_defs, only: ELECTRON
  use lorentz
  use pdg_arrays
  use flavors
  use interactions, only: reset_interaction_counter
  use model_data
  use sf_aux
  use sf_base

  use sf_escan

  implicit none
  private

  public :: sf_escan_1
  public :: sf_escan_2

contains

  subroutine sf_escan_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t), dimension(2) :: pdg_in
    type(pdg_array_t), dimension(2) :: pdg_out
    integer, dimension(:), allocatable :: pdg1, pdg2
    class(sf_data_t), allocatable :: data

    write (u, "(A)")  "* Test output: sf_escan_1"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &energy-scan structure function data"
    write (u, "(A)")

    call model%init_qed_test ()
    pdg_in(1) = ELECTRON
    pdg_in(2) = -ELECTRON

    allocate (escan_data_t :: data)
    select type (data)
    type is (escan_data_t)
       call data%init (model, pdg_in, norm = 2._default)
    end select

    call data%write (u)

    write (u, "(A)")

    write (u, "(1x,A)")  "Outgoing particle codes:"
    call data%get_pdg_out (pdg_out)
    pdg1 = pdg_out(1)
    pdg2 = pdg_out(2)
    write (u, "(2x,99(1x,I0))")  pdg1, pdg2

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_escan_1"

  end subroutine sf_escan_1

  subroutine sf_escan_2 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t), dimension(2) :: flv
    type(pdg_array_t), dimension(2) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k1, k2
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: x_free, f

    write (u, "(A)")  "* Test output: sf_escan_2"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &beam-events structure function data"
    write (u, "(A)")

    call model%init_qed_test ()
    call flv(1)%init (ELECTRON, model)
    call flv(2)%init (-ELECTRON, model)
    pdg_in(1) = ELECTRON
    pdg_in(2) = -ELECTRON

    call reset_interaction_counter ()

    allocate (escan_data_t :: data)
    select type (data)
    type is (escan_data_t)
       call data%init (model, pdg_in)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1,2])

    write (u, "(A)")  "* Initialize incoming momentum with E=500"
    write (u, "(A)")
    E = 250
    k1 = vector4_moving (E, sqrt (E**2 - flv(1)%get_mass ()**2), 3)
    k2 = vector4_moving (E,-sqrt (E**2 - flv(2)%get_mass ()**2), 3)
    call vector4_write (k1, u)
    call vector4_write (k2, u)
    call sf_int%seed_kinematics ([k1, k2])

    write (u, "(A)")
    write (u, "(A)")  "* Set dummy parameters and generate x"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r  = 0.8
    rb = 1 - r
    x_free = 1

    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)

    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f
    write (u, "(A,9(1x,F10.7))")  "xf=", x_free

    write (u, "(A)")
    write (u, "(A)")  "* Inverse kinematics"
    write (u, "(A)")

    call sf_int%recover_x (x, xb, x_free)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)

    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f
    write (u, "(A,9(1x,F10.7))")  "xf=", x_free

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate"
    write (u, "(A)")

    call sf_int%apply (scale = 0._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_escan_2"

  end subroutine sf_escan_2


end module sf_escan_uti
