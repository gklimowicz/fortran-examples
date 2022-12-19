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

module sf_ewa_uti

  use kinds, only: default
  use lorentz
  use pdg_arrays
  use flavors
  use interactions, only: reset_interaction_counter
  use model_data
  use sf_aux
  use sf_base

  use sf_ewa

  implicit none
  private

  public :: sf_ewa_1
  public :: sf_ewa_2
  public :: sf_ewa_3
  public :: sf_ewa_4
  public :: sf_ewa_5

contains

  subroutine sf_ewa_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(pdg_array_t), dimension(1) :: pdg_out
    integer, dimension(:), allocatable :: pdg1
    class(sf_data_t), allocatable :: data

    write (u, "(A)")  "* Test output: sf_ewa_1"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &test structure function data"
    write (u, "(A)")

    write (u, "(A)")  "* Create empty data object"
    write (u, "(A)")

    call model%init_sm_test ()
    pdg_in = 2

    allocate (ewa_data_t :: data)
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize for Z boson"
    write (u, "(A)")

    select type (data)
    type is (ewa_data_t)
       call data%init (model, pdg_in, 0.01_default, &
            500._default, 5000._default, .false., .false.)
       call data%set_id (23)
    end select

    call data%write (u)

    write (u, "(A)")

    write (u, "(1x,A)")  "Outgoing particle codes:"
    call data%get_pdg_out (pdg_out)
    pdg1 = pdg_out(1)
    write (u, "(2x,99(1x,I0))")  pdg1

    write (u, "(A)")
    write (u, "(A)")  "* Initialize for W boson"
    write (u, "(A)")

    deallocate (data)
    allocate (ewa_data_t :: data)
    select type (data)
    type is (ewa_data_t)
       call data%init (model, pdg_in, 0.01_default, &
            500._default, 5000._default, .false., .false.)
       call data%set_id (24)
    end select

    call data%write (u)

    write (u, "(A)")

    write (u, "(1x,A)")  "Outgoing particle codes:"
    call data%get_pdg_out (pdg_out)
    pdg1 = pdg_out(1)
    write (u, "(2x,99(1x,I0))")  pdg1

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_ewa_1"

  end subroutine sf_ewa_1

  subroutine sf_ewa_2 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k
    type(vector4_t), dimension(2) :: q
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f

    write (u, "(A)")  "* Test output: sf_ewa_2"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_sm_test ()
    call flv%init (2, model)
    pdg_in = 2

    call reset_interaction_counter ()

    allocate (ewa_data_t :: data)
    select type (data)
    type is (ewa_data_t)
       call data%init (model, pdg_in, 0.01_default, &
            500._default, 3000._default, .false., .true.)
       call data%set_id (24)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])
    call sf_int%setup_constants ()

    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize incoming momentum with E=1500"
    write (u, "(A)")
    E = 1500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call pacify (k, 1e-10_default)
    call vector4_write (k, u)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for r=0.4, no EWA mapping, collinear"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = 0.4_default
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)

    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Recover x from momenta"
    write (u, "(A)")

    q = sf_int%get_momenta (outgoing=.true.)
    call sf_int%final ()
    deallocate (sf_int)

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])
    call sf_int%setup_constants ()

    call sf_int%seed_kinematics ([k])
    call sf_int%set_momenta (q, outgoing=.true.)
    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false., &
         set_momenta=.true.)

    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate EWA structure function"
    write (u, "(A)")

    call sf_int%apply (scale = 100._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_ewa_2"

  end subroutine sf_ewa_2

  subroutine sf_ewa_3 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k
    type(vector4_t), dimension(2) :: q
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f

    write (u, "(A)")  "* Test output: sf_ewa_3"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_sm_test ()
    call flv%init (2, model)
    pdg_in = 2

    call reset_interaction_counter ()

    allocate (ewa_data_t :: data)
    select type (data)
    type is (ewa_data_t)
       call data%init (model, pdg_in, 0.01_default, &
            500._default, 3000._default, .false., .true.)
       call data%set_id (24)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])
    call sf_int%setup_constants ()

    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize incoming momentum with E=1500"
    write (u, "(A)")
    E = 1500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call pacify (k, 1e-10_default)
    call vector4_write (k, u)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for r=0.4, with EWA mapping, collinear"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = 0.4_default
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.true.)

    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Recover x from momenta"
    write (u, "(A)")

    q = sf_int%get_momenta (outgoing=.true.)
    call sf_int%final ()
    deallocate (sf_int)

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])
    call sf_int%setup_constants ()

    call sf_int%seed_kinematics ([k])
    call sf_int%set_momenta (q, outgoing=.true.)
    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.true., &
         set_momenta=.true.)

    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate EWA structure function"
    write (u, "(A)")

    call sf_int%apply (scale = 100._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_ewa_3"

  end subroutine sf_ewa_3

  subroutine sf_ewa_4 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k
    type(vector4_t), dimension(2) :: q
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f

    write (u, "(A)")  "* Test output: sf_ewa_4"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call modeL%init_sm_test ()
    call flv%init (2, model)
    pdg_in = 2

    call reset_interaction_counter ()

    allocate (ewa_data_t :: data)
    select type (data)
    type is (ewa_data_t)
       call data%init (model, pdg_in, 0.01_default, &
            500._default, 3000.0_default, .true., .true.)
       call data%set_id (24)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])
    call sf_int%setup_constants ()

    write (u, "(A)")  "* Initialize incoming momentum with E=1500"
    write (u, "(A)")
    E = 1500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call pacify (k, 1e-10_default)
    call vector4_write (k, u)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for r=0.5/0.5/0.25, with EWA mapping, "
    write (u, "(A)")  "          non-coll., keeping energy"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = [0.5_default, 0.5_default, 0.25_default]
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_ENERGY
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.true.)
    call sf_int%pacify_momenta (1e-10_default)

    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r from momenta"
    write (u, "(A)")

    q = sf_int%get_momenta (outgoing=.true.)
    call sf_int%final ()
    deallocate (sf_int)

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])
    call sf_int%setup_constants ()

    call sf_int%seed_kinematics ([k])
    call sf_int%set_momenta (q, outgoing=.true.)
    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.true., &
         set_momenta=.true.)
    call sf_int%pacify_momenta (1e-10_default)

    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate EWA structure function"
    write (u, "(A)")

    call sf_int%apply (scale = 1500._default)
    call sf_int%write (u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_ewa_4"

  end subroutine sf_ewa_4

  subroutine sf_ewa_5 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f

    write (u, "(A)")  "* Test output: sf_ewa_5"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_sm_test ()
    call flv%init (2, model)
    pdg_in = [1, 2, -1, -2]

    call reset_interaction_counter ()

    allocate (ewa_data_t :: data)
    select type (data)
    type is (ewa_data_t)
       call data%init (model, pdg_in, 0.01_default, &
            500._default, 3000._default, .false., .true.)
       call data%set_id (24)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])
    call sf_int%setup_constants ()

    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize incoming momentum with E=1500"
    write (u, "(A)")
    E = 1500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call pacify (k, 1e-10_default)
    call vector4_write (k, u)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for r=0.4, no EWA mapping, collinear"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = 0.4_default
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)

    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate EWA structure function"
    write (u, "(A)")

    call sf_int%apply (scale = 100._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_ewa_5"

  end subroutine sf_ewa_5


end module sf_ewa_uti
