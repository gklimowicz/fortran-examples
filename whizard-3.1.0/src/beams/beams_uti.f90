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

module beams_uti

  use kinds, only: default
  use lorentz
  use flavors
  use interactions, only: reset_interaction_counter
  use polarizations, only: smatrix_t
  use model_data
  use beam_structures

  use beams

  implicit none
  private

  public :: beam_1
  public :: beam_2
  public :: beam_3

contains

  subroutine beam_1 (u)
    integer, intent(in) :: u
    type(beam_data_t), target :: beam_data
    type(beam_t) :: beam
    real(default) :: sqrts
    type(flavor_t), dimension(2) :: flv
    type(smatrix_t), dimension(2) :: smatrix
    real(default), dimension(2) :: pol_f
    type(model_data_t), target :: model

    write (u, "(A)")  "* Test output: beam_1"
    write (u, "(A)")  "*   Purpose: test basic beam setup"
    write (u, "(A)")

    write (u, "(A)")  "* Reading model file"
    write (u, "(A)")

    call reset_interaction_counter ()

    call model%init_sm_test ()

    write (u, "(A)")  "* Unpolarized scattering, massless fermions"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([1,-1], model)

    call beam_data%init_sqrts (sqrts, flv)
    call beam_data%write (u)
    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Unpolarized scattering, massless bosons"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([22,22], model)

    call beam_data%init_sqrts (sqrts, flv)
    call beam_data%write (u)
    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Unpolarized scattering, massive bosons"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([24,-24], model)

    call beam_data%init_sqrts (sqrts, flv)
    call beam_data%write (u)
    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Polarized scattering, massless fermions"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([1,-1], model)

    call smatrix(1)%init (2, 1)
    call smatrix(1)%set_entry (1, [1,1], (1._default, 0._default))
    pol_f(1) = 0.5_default

    call smatrix(2)%init (2, 3)
    call smatrix(2)%set_entry (1, [1,1], (1._default, 0._default))
    call smatrix(2)%set_entry (2, [-1,-1], (1._default, 0._default))
    call smatrix(2)%set_entry (3, [-1,1], (1._default, 0._default))
    pol_f(2) = 1._default

    call beam_data%init_sqrts (sqrts, flv, smatrix, pol_f)
    call beam_data%write (u)
    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Semi-polarized scattering, massless bosons"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([22,22], model)

    call smatrix(1)%init (2, 0)
    pol_f(1) = 0._default

    call smatrix(2)%init (2, 1)
    call smatrix(2)%set_entry (1, [1,1], (1._default, 0._default))
    pol_f(2) = 1._default

    call beam_data%init_sqrts (sqrts, flv, smatrix, pol_f)
    call beam_data%write (u)
    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Semi-polarized scattering, massive bosons"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([24,-24], model)

    call smatrix(1)%init (2, 0)
    pol_f(1) = 0._default

    call smatrix(2)%init (2, 1)
    call smatrix(2)%set_entry (1, [0,0], (1._default, 0._default))
    pol_f(2) = 1._default

    call beam_data%init_sqrts (sqrts, flv, smatrix, pol_f)
    call beam_data%write (u)
    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Unpolarized decay, massive boson"
    write (u, "(A)")

    call reset_interaction_counter ()
    call flv(1)%init (23, model)

    call beam_data%init_decay (flv(1:1))
    call beam_data%write (u)
    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)

    write (u, "(A)")
    write (u, "(A)")  "* Polarized decay, massive boson"
    write (u, "(A)")

    call reset_interaction_counter ()
    call flv(1)%init (23, model)
    call smatrix(1)%init (2, 1)
    call smatrix(1)%set_entry (1, [0,0], (1._default, 0._default))
    pol_f(1) = 0.4_default

    call beam_data%init_decay (flv(1:1), smatrix(1:1), pol_f(1:1))
    call beam_data%write (u)
    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call beam_final (beam)
    call beam_data%final ()

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: beam_1"

  end subroutine beam_1

  subroutine beam_2 (u)
    integer, intent(in) :: u
    type(beam_data_t), target :: beam_data
    type(beam_t) :: beam
    real(default) :: sqrts
    type(flavor_t), dimension(2) :: flv
    integer, dimension(0) :: no_records
    type(beam_structure_t) :: beam_structure
    type(model_data_t), target :: model

    write (u, "(A)")  "* Test output: beam_2"
    write (u, "(A)")  "*   Purpose: transfer beam polarization using &
         &beam structure"
    write (u, "(A)")

    write (u, "(A)")  "* Reading model file"
    write (u, "(A)")

    call model%init_sm_test ()

    write (u, "(A)")  "* Unpolarized scattering, massless fermions"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([1,-1], model)
    call beam_structure%init_sf (flv%get_name (), no_records)
    call beam_structure%final_pol ()

    call beam_structure%write (u)
    write (u, *)

    call beam_data%init_structure (beam_structure, sqrts, model)
    call beam_data%write (u)

    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Unpolarized scattering, massless bosons"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([22,22], model)

    call beam_structure%init_sf (flv%get_name (), no_records)
    call beam_structure%final_pol ()

    call beam_structure%write (u)
    write (u, *)

    call beam_data%init_structure (beam_structure, sqrts, model)
    call beam_data%write (u)

    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Unpolarized scattering, massive bosons"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([24,-24], model)

    call beam_structure%init_sf (flv%get_name (), no_records)
    call beam_structure%final_pol ()

    call beam_structure%write (u)
    write (u, *)

    call beam_data%init_structure (beam_structure, sqrts, model)
    call beam_data%write (u)

    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Polarized scattering, massless fermions"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([1,-1], model)
    call beam_structure%init_sf (flv%get_name (), no_records)
    call beam_structure%init_pol (2)

    call beam_structure%init_smatrix (1, 1)
    call beam_structure%set_sentry (1, 1, [1,1], (1._default, 0._default))

    call beam_structure%init_smatrix (2, 3)
    call beam_structure%set_sentry (2, 1, [1,1], (1._default, 0._default))
    call beam_structure%set_sentry (2, 2, [-1,-1], (1._default, 0._default))
    call beam_structure%set_sentry (2, 3, [-1,1], (1._default, 0._default))

    call beam_structure%set_pol_f ([0.5_default, 1._default])
    call beam_structure%write (u)
    write (u, *)

    call beam_data%init_structure (beam_structure, sqrts, model)
    call beam_data%write (u)
    write (u, *)

    call beam_init (beam, beam_data)
    call beam_write (beam, u)

    call beam_final (beam)
    call beam_data%final ()
    call beam_structure%final_pol ()
    call beam_structure%final_sf ()

    write (u, "(A)")
    write (u, "(A)")  "* Semi-polarized scattering, massless bosons"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([22,22], model)

    call beam_structure%init_sf (flv%get_name (), no_records)
    call beam_structure%init_pol (2)

    call beam_structure%init_smatrix (1, 0)

    call beam_structure%init_smatrix (2, 1)
    call beam_structure%set_sentry (2, 1, [1,1], (1._default, 0._default))

    call beam_structure%set_pol_f ([0._default, 1._default])
    call beam_structure%write (u)
    write (u, *)

    call beam_data%init_structure (beam_structure, sqrts, model)
    call beam_data%write (u)

    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Semi-polarized scattering, massive bosons"
    write (u, "(A)")

    call reset_interaction_counter ()
    sqrts = 500
    call flv%init ([24,-24], model)

    call beam_structure%init_sf (flv%get_name (), no_records)
    call beam_structure%init_pol (2)

    call beam_structure%init_smatrix (1, 0)

    call beam_structure%init_smatrix (2, 1)
    call beam_structure%set_sentry (2, 1, [0,0], (1._default, 0._default))
    call beam_structure%write (u)

    write (u, "(A)")
    call beam_data%init_structure (beam_structure, sqrts, model)
    call beam_data%write (u)

    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)
    call beam_final (beam)
    call beam_data%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Unpolarized decay, massive boson"
    write (u, "(A)")

    call reset_interaction_counter ()
    call flv(1)%init (23, model)

    call beam_structure%init_sf ([flv(1)%get_name ()], no_records)
    call beam_structure%final_pol ()
    call beam_structure%write (u)

    write (u, "(A)")
    call beam_data%init_structure (beam_structure, sqrts, model)
    call beam_data%write (u)

    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)

    write (u, "(A)")
    write (u, "(A)")  "* Polarized decay, massive boson"
    write (u, "(A)")

    call reset_interaction_counter ()
    call flv(1)%init (23, model)
    call beam_structure%init_sf ([flv(1)%get_name ()], no_records)

    call beam_structure%init_pol (1)

    call beam_structure%init_smatrix (1, 1)
    call beam_structure%set_sentry (1, 1, [0,0], (1._default, 0._default))
    call beam_structure%set_pol_f ([0.4_default])
    call beam_structure%write (u)
    write (u, *)

    call beam_data%init_structure (beam_structure, sqrts, model)
    call beam_data%write (u)
    write (u, "(A)")
    call beam_init (beam, beam_data)
    call beam_write (beam, u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call beam_final (beam)
    call beam_data%final ()

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: beam_2"

  end subroutine beam_2

  subroutine beam_3 (u)
    integer, intent(in) :: u
    type(beam_data_t), target :: beam_data
    type(beam_t) :: beam
    type(flavor_t), dimension(2) :: flv
    integer, dimension(0) :: no_records
    type(model_data_t), target :: model
    type(beam_structure_t) :: beam_structure
    type(vector3_t), dimension(2) :: p3
    type(vector4_t), dimension(2) :: p

    write (u, "(A)")  "* Test output: beam_3"
    write (u, "(A)")  "*   Purpose: set up beams with generic momenta"
    write (u, "(A)")

    write (u, "(A)")  "* Reading model file"
    write (u, "(A)")

    call reset_interaction_counter ()

    call model%init_sm_test ()

    write (u, "(A)")  "* 1: Scattering process"
    write (u, "(A)")

    call flv%init ([2212,2212], model)

    p3(1) = vector3_moving ([5._default, 0._default, 10._default])
    p3(2) = -vector3_moving ([1._default, 1._default, -10._default])

    call beam_structure%init_sf (flv%get_name (), no_records)
    call beam_structure%set_momentum (p3 ** 1)
    call beam_structure%set_theta (polar_angle (p3))
    call beam_structure%set_phi (azimuthal_angle (p3))
    call beam_structure%write (u)
    write (u, *)

    call beam_data%init_structure (beam_structure, 0._default, model)
    call pacify (beam_data%l_cm_to_lab, 1e-20_default)
    call beam_data%compute_md5sum ()
    call beam_data%write (u, verbose = .true.)
    write (u, *)

    write (u, "(1x,A)")  "Beam momenta reconstructed from LT:"
    p = beam_data%L_cm_to_lab * beam_data%p_cm
    call pacify (p, 1e-12_default)
    call vector4_write (p(1), u)
    call vector4_write (p(2), u)
    write (u, "(A)")

    call beam_init (beam, beam_data)
    call beam_write (beam, u)

    call beam_final (beam)
    call beam_data%final ()
    call beam_structure%final_sf ()
    call beam_structure%final_mom ()

    write (u, "(A)")
    write (u, "(A)")  "* 2: Decay"
    write (u, "(A)")

    call flv(1)%init (23, model)
    p3(1) = vector3_moving ([10._default, 5._default, 50._default])

    call beam_structure%init_sf ([flv(1)%get_name ()], no_records)
    call beam_structure%set_momentum ([p3(1) ** 1])
    call beam_structure%set_theta ([polar_angle (p3(1))])
    call beam_structure%set_phi ([azimuthal_angle (p3(1))])
    call beam_structure%write (u)
    write (u, *)

    call beam_data%init_structure (beam_structure, 0._default, model)
    call beam_data%write (u, verbose = .true.)
    write (u, "(A)")

    write (u, "(1x,A)")  "Beam momentum reconstructed from LT:"
    p(1) = beam_data%L_cm_to_lab * beam_data%p_cm(1)
    call pacify (p(1), 1e-12_default)
    call vector4_write (p(1), u)
    write (u, "(A)")

    call beam_init (beam, beam_data)
    call beam_write (beam, u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call beam_final (beam)
    call beam_data%final ()
    call beam_structure%final_sf ()
    call beam_structure%final_mom ()

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: beam_3"

  end subroutine beam_3


end module beams_uti
