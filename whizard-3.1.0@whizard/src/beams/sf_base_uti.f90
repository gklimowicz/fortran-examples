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

module sf_base_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use format_defs, only: FMT_19
  use format_utils, only: write_separator
  use diagnostics
  use lorentz
  use pdg_arrays
  use flavors
  use colors
  use helicities
  use quantum_numbers
  use state_matrices, only: FM_IGNORE_HELICITY
  use interactions
  use particles
  use model_data
  use beams
  use sf_aux
  use sf_mappings

  use sf_base

  implicit none
  private

  public :: sf_base_1
  public :: sf_base_2
  public :: sf_base_3
  public :: sf_base_4
  public :: sf_base_5
  public :: sf_base_6
  public :: sf_base_7
  public :: sf_base_8
  public :: sf_base_9
  public :: sf_base_10
  public :: sf_base_11
  public :: sf_base_12
  public :: sf_base_13
  public :: sf_base_14

  public :: sf_test_data_t

  type, extends (sf_data_t) :: sf_test_data_t
     class(model_data_t), pointer :: model => null ()
     integer :: mode = 0
     type(flavor_t) :: flv_in
     type(flavor_t) :: flv_out
     type(flavor_t) :: flv_rad
     real(default) :: m = 0
     logical :: collinear = .true.
     real(default), dimension(:), allocatable :: qbounds
   contains
     procedure :: write => sf_test_data_write
     procedure :: init => sf_test_data_init
     procedure :: get_n_par => sf_test_data_get_n_par
     procedure :: get_pdg_out => sf_test_data_get_pdg_out
     procedure :: allocate_sf_int => sf_test_data_allocate_sf_int
  end type sf_test_data_t

  type, extends (sf_int_t) :: sf_test_t
     type(sf_test_data_t), pointer :: data => null ()
     real(default) :: x = 0
   contains
     procedure :: type_string => sf_test_type_string
     procedure :: write => sf_test_write
     procedure :: init => sf_test_init
     procedure :: complete_kinematics => sf_test_complete_kinematics
     procedure :: inverse_kinematics => sf_test_inverse_kinematics
     procedure :: apply => sf_test_apply
  end type sf_test_t

  type, extends (sf_data_t) :: sf_test_spectrum_data_t
     class(model_data_t), pointer :: model => null ()
     type(flavor_t) :: flv_in
     type(flavor_t) :: flv_out
     type(flavor_t) :: flv_rad
     logical :: with_radiation = .true.
     real(default) :: m = 0
   contains
     procedure :: write => sf_test_spectrum_data_write
     procedure :: init => sf_test_spectrum_data_init
     procedure :: get_n_par => sf_test_spectrum_data_get_n_par
     procedure :: get_pdg_out => sf_test_spectrum_data_get_pdg_out
     procedure :: allocate_sf_int => &
          sf_test_spectrum_data_allocate_sf_int
  end type sf_test_spectrum_data_t

  type, extends (sf_int_t) :: sf_test_spectrum_t
     type(sf_test_spectrum_data_t), pointer :: data => null ()
   contains
     procedure :: type_string => sf_test_spectrum_type_string
     procedure :: write => sf_test_spectrum_write
     procedure :: init => sf_test_spectrum_init
     procedure :: complete_kinematics => sf_test_spectrum_complete_kinematics
     procedure :: inverse_kinematics => sf_test_spectrum_inverse_kinematics
     procedure :: apply => sf_test_spectrum_apply
  end type sf_test_spectrum_t

  type, extends (sf_data_t) :: sf_test_generator_data_t
     class(model_data_t), pointer :: model => null ()
     type(flavor_t) :: flv_in
     type(flavor_t) :: flv_out
     type(flavor_t) :: flv_rad
     real(default) :: m = 0
   contains
     procedure :: write => sf_test_generator_data_write
     procedure :: init => sf_test_generator_data_init
     procedure :: is_generator => sf_test_generator_data_is_generator
     procedure :: get_n_par => sf_test_generator_data_get_n_par
     procedure :: get_pdg_out => sf_test_generator_data_get_pdg_out
     procedure :: allocate_sf_int => &
          sf_test_generator_data_allocate_sf_int
  end type sf_test_generator_data_t

  type, extends (sf_int_t) :: sf_test_generator_t
     type(sf_test_generator_data_t), pointer :: data => null ()
   contains
     procedure :: type_string => sf_test_generator_type_string
     procedure :: write => sf_test_generator_write
     procedure :: init => sf_test_generator_init
     procedure :: is_generator => sf_test_generator_is_generator
     procedure :: generate_free => sf_test_generator_generate_free
     procedure :: recover_x => sf_test_generator_recover_x
     procedure :: complete_kinematics => sf_test_generator_complete_kinematics
     procedure :: inverse_kinematics => sf_test_generator_inverse_kinematics
     procedure :: apply => sf_test_generator_apply
  end type sf_test_generator_t


contains

  subroutine sf_base_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(pdg_array_t), dimension(1) :: pdg_out
    integer, dimension(:), allocatable :: pdg1
    class(sf_data_t), allocatable :: data

    write (u, "(A)")  "* Test output: sf_base_1"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &test structure function data"
    write (u, "(A)")

    call model%init_test ()
    pdg_in = 25

    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (model, pdg_in)
    end select

    call data%write (u)

    write (u, "(A)")

    write (u, "(1x,A)")  "Outgoing particle code:"
    call data%get_pdg_out (pdg_out)
    pdg1 = pdg_out(1)
    write (u, "(2x,99(1x,I0))")  pdg1

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_1"

  end subroutine sf_base_1

  subroutine sf_base_2 (u)
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

    write (u, "(A)")  "* Test output: sf_base_2"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    pdg_in = 25
    call flv%init (25, model)

    call reset_interaction_counter ()

    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (model, pdg_in)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])

    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize incoming momentum with E=500"
    write (u, "(A)")
    E = 500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call vector4_write (k, u)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = 0
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=1"
    write (u, "(A)")

    r = 1
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5"
    write (u, "(A)")

    r = 0.5_default
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics with mapping for r=0.8"
    write (u, "(A)")

    r = 0.8_default
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.true.)
    call sf_int%write (u)

    write (u, "(A)")
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

    call sf_int%seed_kinematics ([k])
    call sf_int%set_momenta (q, outgoing=.true.)
    call sf_int%recover_x (x, xb)

    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb

    write (u, "(A)")
    write (u, "(A)")  "* Compute inverse kinematics for x=0.64 and evaluate"
    write (u, "(A)")

    x = 0.64_default
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.true.)
    call sf_int%apply (scale=0._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_2"

  end subroutine sf_base_2

  subroutine sf_base_3 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(flavor_t) :: flv
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f

    write (u, "(A)")  "* Test output: sf_base_3"
    write (u, "(A)")  "*   Purpose: check various kinematical setups"
    write (u, "(A)")  "*            for collinear structure-function splitting."
    write (u, "(A)")  "             (two masses equal, one zero)"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    pdg_in = 25
    call flv%init (25, model)

    call reset_interaction_counter ()

    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (model, pdg_in)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    call sf_int%write (u)

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    write (u, "(A)")
    write (u, "(A)")  "* Initialize incoming momentum with E=500"

    E = 500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set radiated mass to zero"

    sf_int%mr2 = 0
    sf_int%mo2 = sf_int%mi2

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5, keeping energy"
    write (u, "(A)")

    r = 0.5_default
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_ENERGY
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5, keeping momentum"
    write (u, "(A)")

    r = 0.5_default
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_MOMENTUM
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set outgoing mass to zero"

    sf_int%mr2 = sf_int%mi2
    sf_int%mo2 = 0

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5, keeping energy"
    write (u, "(A)")

    r = 0.5_default
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_ENERGY
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5, keeping momentum"
    write (u, "(A)")

    r = 0.5_default
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_MOMENTUM
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set incoming mass to zero"

    k = vector4_moving (E, E, 3)
    call sf_int%seed_kinematics ([k])

    sf_int%mr2 = sf_int%mi2
    sf_int%mo2 = sf_int%mi2
    sf_int%mi2 = 0

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5, keeping energy"
    write (u, "(A)")

    r = 0.5_default
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_ENERGY
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5, keeping momentum"
    write (u, "(A)")

    r = 0.5_default
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_MOMENTUM
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set all masses to zero"

    sf_int%mr2 = 0
    sf_int%mo2 = 0
    sf_int%mi2 = 0

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5, keeping energy"
    write (u, "(A)")

    r = 0.5_default
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_ENERGY
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5, keeping momentum"
    write (u, "(A)")

    r = 0.5_default
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_MOMENTUM
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_3"

  end subroutine sf_base_3

  subroutine sf_base_4 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(flavor_t) :: flv
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f

    write (u, "(A)")  "* Test output: sf_base_4"
    write (u, "(A)")  "*   Purpose: check various kinematical setups"
    write (u, "(A)")  "*            for free structure-function splitting."
    write (u, "(A)")  "             (two masses equal, one zero)"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    pdg_in = 25
    call flv%init (25, model)

    call reset_interaction_counter ()

    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (model, pdg_in, collinear=.false.)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    call sf_int%write (u)

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    write (u, "(A)")
    write (u, "(A)")  "* Initialize incoming momentum with E=500"

    E = 500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set radiated mass to zero"

    sf_int%mr2 = 0
    sf_int%mo2 = sf_int%mi2

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5/0.5/0.125, keeping energy"
    write (u, "(A)")

    r = [0.5_default, 0.5_default, 0.125_default]
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_ENERGY
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5/0.5/0.125, keeping momentum"
    write (u, "(A)")

    r = [0.5_default, 0.5_default, 0.125_default]
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_MOMENTUM
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set outgoing mass to zero"

    sf_int%mr2 = sf_int%mi2
    sf_int%mo2 = 0

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5/0.5/0.125, keeping energy"
    write (u, "(A)")

    r = [0.5_default, 0.5_default, 0.125_default]
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_ENERGY
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5/0.5/0.125, keeping momentum"
    write (u, "(A)")

    r = [0.5_default, 0.5_default, 0.125_default]
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_MOMENTUM
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set incoming mass to zero"

    k = vector4_moving (E, E, 3)
    call sf_int%seed_kinematics ([k])

    sf_int%mr2 = sf_int%mi2
    sf_int%mo2 = sf_int%mi2
    sf_int%mi2 = 0

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5/0.5/0.125, keeping energy"
    write (u, "(A)")

    r = [0.5_default, 0.5_default, 0.125_default]
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_ENERGY
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5/0.5/0.125, keeping momentum"
    write (u, "(A)")

    r = [0.5_default, 0.5_default, 0.125_default]
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_MOMENTUM
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set all masses to zero"

    sf_int%mr2 = 0
    sf_int%mo2 = 0
    sf_int%mi2 = 0

    write (u, "(A)")
    write (u, "(A)")  "* Re-Initialize structure-function object with Q bounds"

    call reset_interaction_counter ()

    select type (data)
    type is (sf_test_data_t)
       call data%init (model, pdg_in, collinear=.false., &
            qbounds = [1._default, 100._default])
    end select

    call sf_int%init (data)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5/0.5/0.125, keeping energy"
    write (u, "(A)")

    r = [0.5_default, 0.5_default, 0.125_default]
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_ENERGY
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5/0.5/0.125, keeping momentum"
    write (u, "(A)")

    r = [0.5_default, 0.5_default, 0.125_default]
    rb = 1 - r
    sf_int%on_shell_mode = KEEP_MOMENTUM
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover x and r"
    write (u, "(A)")

    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_4"

  end subroutine sf_base_4

  subroutine sf_base_5 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(pdg_array_t), dimension(2) :: pdg_out
    integer, dimension(:), allocatable :: pdg1, pdg2
    type(flavor_t) :: flv
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t), dimension(2) :: k
    type(vector4_t), dimension(4) :: q
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f

    write (u, "(A)")  "* Test output: sf_base_5"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &a pair spectrum object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    allocate (sf_test_spectrum_data_t :: data)
    select type (data)
    type is (sf_test_spectrum_data_t)
       call data%init (model, pdg_in, with_radiation=.true.)
    end select

    write (u, "(1x,A)")  "Outgoing particle codes:"
    call data%get_pdg_out (pdg_out)
    pdg1 = pdg_out(1)
    pdg2 = pdg_out(2)
    write (u, "(2x,99(1x,I0))")  pdg1, pdg2

    write (u, "(A)")
    write (u, "(A)")  "* Initialize spectrum object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize incoming momenta with sqrts=1000"

    E = 500
    k(1) = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    k(2) = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call sf_int%seed_kinematics (k)

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.4,0.8"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = [0.4_default, 0.8_default]
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics with mapping for r=0.6,0.8"
    write (u, "(A)")

    r = [0.6_default, 0.8_default]
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.true.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Recover x from momenta"
    write (u, "(A)")

    q = sf_int%get_momenta (outgoing=.true.)
    call sf_int%final ()
    deallocate (sf_int)

    call reset_interaction_counter ()
    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    call sf_int%seed_kinematics (k)
    call sf_int%set_momenta (q, outgoing=.true.)
    call sf_int%recover_x (x, xb)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb

    write (u, "(A)")
    write (u, "(A)")  "* Compute inverse kinematics for x=0.36,0.64 &
         &and evaluate"
    write (u, "(A)")

    x = [0.36_default, 0.64_default]
    xb = 1 - x
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.true.)
    call sf_int%apply (scale=0._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_5"

  end subroutine sf_base_5

  subroutine sf_base_6 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(flavor_t) :: flv
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t), dimension(2) :: k
    type(vector4_t), dimension(2) :: q
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f

    write (u, "(A)")  "* Test output: sf_base_6"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &a pair spectrum object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    allocate (sf_test_spectrum_data_t :: data)
    select type (data)
    type is (sf_test_spectrum_data_t)
       call data%init (model, pdg_in, with_radiation=.false.)
    end select

    write (u, "(A)")  "* Initialize spectrum object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    write (u, "(A)")  "* Initialize incoming momenta with sqrts=1000"

    E = 500
    k(1) = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    k(2) = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call sf_int%seed_kinematics (k)

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.4,0.8"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = [0.4_default, 0.8_default]
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Recover x from momenta"
    write (u, "(A)")

    q = sf_int%get_momenta (outgoing=.true.)
    call sf_int%final ()
    deallocate (sf_int)

    call reset_interaction_counter ()
    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    call sf_int%seed_kinematics (k)
    call sf_int%set_momenta (q, outgoing=.true.)
    call sf_int%recover_x (x, xb)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb

    write (u, "(A)")
    write (u, "(A)")  "* Compute inverse kinematics for x=0.4,0.8 &
         &and evaluate"
    write (u, "(A)")

    x = [0.4_default, 0.8_default]
    xb = 1 - x
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%apply (scale=0._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_6"

  end subroutine sf_base_6

  subroutine sf_base_7 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(flavor_t) :: flv
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    real(default), dimension(:), allocatable :: value

    write (u, "(A)")  "* Test output: sf_base_7"
    write (u, "(A)")  "*   Purpose: check direct access method"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (model, pdg_in)
    end select

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    write (u, "(A)")  "* Probe structure function: states"
    write (u, "(A)")

    write (u, "(A,I0)")  "n_states = ", sf_int%get_n_states ()
    write (u, "(A,I0)")  "n_in     = ", sf_int%get_n_in ()
    write (u, "(A,I0)")  "n_rad    = ", sf_int%get_n_rad ()
    write (u, "(A,I0)")  "n_out    = ", sf_int%get_n_out ()
    write (u, "(A)")
    write (u, "(A)", advance="no")  "state(1)  = "
    call quantum_numbers_write (sf_int%get_state (1), u)
    write (u, *)

    allocate (value (sf_int%get_n_states ()))
    call sf_int%compute_values (value, &
         E=[500._default], x=[0.5_default], xb=[0.5_default], scale=0._default)

    write (u, "(A)")
    write (u, "(A)", advance="no")  "value (E=500, x=0.5) ="
    write (u, "(9(1x," // FMT_19 // "))")  value

    call sf_int%compute_values (value, &
         x=[0.1_default], xb=[0.9_default], scale=0._default)

    write (u, "(A)")
    write (u, "(A)", advance="no")  "value (E=500, x=0.1) ="
    write (u, "(9(1x," // FMT_19 // "))")  value


    write (u, "(A)")
    write (u, "(A)")  "* Initialize spectrum object"
    write (u, "(A)")

    deallocate (value)
    call sf_int%final ()
    deallocate (sf_int)
    deallocate (data)

    allocate (sf_test_spectrum_data_t :: data)
    select type (data)
    type is (sf_test_spectrum_data_t)
       call data%init (model, pdg_in, with_radiation=.false.)
    end select

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    write (u, "(A)")  "* Probe spectrum: states"
    write (u, "(A)")

    write (u, "(A,I0)")  "n_states = ", sf_int%get_n_states ()
    write (u, "(A,I0)")  "n_in     = ", sf_int%get_n_in ()
    write (u, "(A,I0)")  "n_rad    = ", sf_int%get_n_rad ()
    write (u, "(A,I0)")  "n_out    = ", sf_int%get_n_out ()
    write (u, "(A)")
    write (u, "(A)", advance="no")  "state(1)  = "
    call quantum_numbers_write (sf_int%get_state (1), u)
    write (u, *)

    allocate (value (sf_int%get_n_states ()))
    call sf_int%compute_value (1, value(1), &
         E = [500._default, 500._default], &
         x = [0.5_default, 0.6_default], &
         xb= [0.5_default, 0.4_default], &
         scale = 0._default)

    write (u, "(A)")
    write (u, "(A)", advance="no")  "value (E=500,500, x=0.5,0.6) ="
    write (u, "(9(1x," // FMT_19 // "))")  value

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_7"

  end subroutine sf_base_7

  subroutine sf_base_8 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    type(beam_data_t), target :: beam_data
    class(sf_data_t), allocatable, target :: data_strfun
    class(sf_data_t), allocatable, target :: data_spectrum
    type(sf_config_t), dimension(:), allocatable :: sf_config
    type(sf_chain_t) :: sf_chain

    write (u, "(A)")  "* Test output: sf_base_8"
    write (u, "(A)")  "*   Purpose: set up a structure-function chain"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    call beam_data%init_sqrts (1000._default, [flv, flv])

    allocate (sf_test_data_t :: data_strfun)
    select type (data_strfun)
    type is (sf_test_data_t)
       call data_strfun%init (model, pdg_in)
    end select

    allocate (sf_test_spectrum_data_t :: data_spectrum)
    select type (data_spectrum)
    type is (sf_test_spectrum_data_t)
       call data_spectrum%init (model, pdg_in, with_radiation=.true.)
    end select

    write (u, "(A)")  "* Set up chain with beams only"
    write (u, "(A)")

    call sf_chain%init (beam_data)
    call write_separator (u, 2)
    call sf_chain%write (u)
    call write_separator (u, 2)
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Set up chain with structure function"
    write (u, "(A)")

    allocate (sf_config (1))
    call sf_config(1)%init ([1], data_strfun)
    call sf_chain%init (beam_data, sf_config)

    call write_separator (u, 2)
    call sf_chain%write (u)
    call write_separator (u, 2)
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Set up chain with spectrum and structure function"
    write (u, "(A)")

    deallocate (sf_config)
    allocate (sf_config (2))
    call sf_config(1)%init ([1,2], data_spectrum)
    call sf_config(2)%init ([2], data_strfun)
    call sf_chain%init (beam_data, sf_config)

    call write_separator (u, 2)
    call sf_chain%write (u)
    call write_separator (u, 2)
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_8"

  end subroutine sf_base_8

  subroutine sf_base_9 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    type(beam_data_t), target :: beam_data
    class(sf_data_t), allocatable, target :: data_strfun
    class(sf_data_t), allocatable, target :: data_spectrum
    type(sf_config_t), dimension(:), allocatable, target :: sf_config
    type(sf_chain_t), target :: sf_chain
    type(sf_chain_instance_t), target :: sf_chain_instance
    type(sf_channel_t), dimension(2) :: sf_channel
    type(vector4_t), dimension(2) :: p
    integer :: j

    write (u, "(A)")  "* Test output: sf_base_9"
    write (u, "(A)")  "*   Purpose: set up a structure-function chain &
         &and create an instance"
    write (u, "(A)")  "*            compute kinematics"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    call beam_data%init_sqrts (1000._default, [flv, flv])

    allocate (sf_test_data_t :: data_strfun)
    select type (data_strfun)
    type is (sf_test_data_t)
       call data_strfun%init (model, pdg_in)
    end select

    allocate (sf_test_spectrum_data_t :: data_spectrum)
    select type (data_spectrum)
    type is (sf_test_spectrum_data_t)
       call data_spectrum%init (model, pdg_in, with_radiation=.true.)
    end select

    write (u, "(A)")  "* Set up chain with beams only"
    write (u, "(A)")

    call sf_chain%init (beam_data)

    call sf_chain_instance%init (sf_chain, n_channel = 1)

    call sf_chain_instance%link_interactions ()
    sf_chain_instance%status = SF_DONE_CONNECTIONS
    call sf_chain_instance%compute_kinematics (1, [real(default) ::])

    call write_separator (u, 2)
    call sf_chain%write (u)
    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    call sf_chain_instance%get_out_momenta (p)

    write (u, "(A)")
    write (u, "(A)")  "* Outgoing momenta:"

    do j = 1, 2
       write (u, "(A)")
       call vector4_write (p(j), u)
    end do

    call sf_chain_instance%final ()
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Set up chain with structure function"
    write (u, "(A)")

    allocate (sf_config (1))
    call sf_config(1)%init ([1], data_strfun)
    call sf_chain%init (beam_data, sf_config)

    call sf_chain_instance%init (sf_chain, n_channel = 1)

    call sf_channel(1)%init (1)
    call sf_channel(1)%activate_mapping ([1])
    call sf_chain_instance%set_channel (1, sf_channel(1))

    call sf_chain_instance%link_interactions ()
    sf_chain_instance%status = SF_DONE_CONNECTIONS
    call sf_chain_instance%compute_kinematics (1, [0.8_default])

    call write_separator (u, 2)
    call sf_chain%write (u)
    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    call sf_chain_instance%get_out_momenta (p)

    write (u, "(A)")
    write (u, "(A)")  "* Outgoing momenta:"

    do j = 1, 2
       write (u, "(A)")
       call vector4_write (p(j), u)
    end do

    call sf_chain_instance%final ()
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Set up chain with spectrum and structure function"
    write (u, "(A)")

    deallocate (sf_config)
    allocate (sf_config (2))
    call sf_config(1)%init ([1,2], data_spectrum)
    call sf_config(2)%init ([2], data_strfun)
    call sf_chain%init (beam_data, sf_config)

    call sf_chain_instance%init (sf_chain, n_channel = 1)

    call sf_channel(2)%init (2)
    call sf_channel(2)%activate_mapping ([2])
    call sf_chain_instance%set_channel (1, sf_channel(2))

    call sf_chain_instance%link_interactions ()
    sf_chain_instance%status = SF_DONE_CONNECTIONS
    call sf_chain_instance%compute_kinematics &
         (1, [0.5_default, 0.6_default, 0.8_default])

    call write_separator (u, 2)
    call sf_chain%write (u)
    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    call sf_chain_instance%get_out_momenta (p)

    write (u, "(A)")
    write (u, "(A)")  "* Outgoing momenta:"

    do j = 1, 2
       write (u, "(A)")
       call vector4_write (p(j), u)
    end do

    call sf_chain_instance%final ()
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_9"

  end subroutine sf_base_9

  subroutine sf_base_10 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    type(beam_data_t), target :: beam_data
    class(sf_data_t), allocatable, target :: data_strfun
    type(sf_config_t), dimension(:), allocatable, target :: sf_config
    type(sf_chain_t), target :: sf_chain
    type(sf_chain_instance_t), target :: sf_chain_instance
    type(sf_channel_t), dimension(2) :: sf_channel
    real(default), dimension(2) :: x_saved

    write (u, "(A)")  "* Test output: sf_base_10"
    write (u, "(A)")  "*   Purpose: set up a structure-function chain"
    write (u, "(A)")  "*            and check mappings"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    call beam_data%init_sqrts (1000._default, [flv, flv])

    allocate (sf_test_data_t :: data_strfun)
    select type (data_strfun)
    type is (sf_test_data_t)
       call data_strfun%init (model, pdg_in)
    end select

    write (u, "(A)")  "* Set up chain with structure function pair &
         &and standard mapping"
    write (u, "(A)")

    allocate (sf_config (2))
    call sf_config(1)%init ([1], data_strfun)
    call sf_config(2)%init ([2], data_strfun)
    call sf_chain%init (beam_data, sf_config)

    call sf_chain_instance%init (sf_chain, n_channel = 1)

    call sf_channel(1)%init (2)
    call sf_channel(1)%set_s_mapping ([1,2])
    call sf_chain_instance%set_channel (1, sf_channel(1))

    call sf_chain_instance%link_interactions ()
    sf_chain_instance%status = SF_DONE_CONNECTIONS
    call sf_chain_instance%compute_kinematics (1, [0.8_default, 0.6_default])

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Invert the kinematics calculation"
    write (u, "(A)")

    x_saved = sf_chain_instance%x

    call sf_chain_instance%init (sf_chain, n_channel = 1)

    call sf_channel(2)%init (2)
    call sf_channel(2)%set_s_mapping ([1, 2])
    call sf_chain_instance%set_channel (1, sf_channel(2))

    call sf_chain_instance%link_interactions ()
    sf_chain_instance%status = SF_DONE_CONNECTIONS
    call sf_chain_instance%inverse_kinematics (x_saved, 1 - x_saved)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)


    call sf_chain_instance%final ()
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_10"

  end subroutine sf_base_10

  subroutine sf_base_11 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    type(beam_data_t), target :: beam_data
    class(sf_data_t), allocatable, target :: data_strfun
    class(sf_data_t), allocatable, target :: data_spectrum
    type(sf_config_t), dimension(:), allocatable, target :: sf_config
    type(sf_chain_t), target :: sf_chain
    type(sf_chain_instance_t), target :: sf_chain_instance
    type(sf_channel_t), dimension(2) :: sf_channel
    type(particle_set_t) :: pset
    type(interaction_t), pointer :: int
    logical :: ok

    write (u, "(A)")  "* Test output: sf_base_11"
    write (u, "(A)")  "*   Purpose: set up a structure-function chain"
    write (u, "(A)")  "*            create an instance and evaluate"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    call beam_data%init_sqrts (1000._default, [flv, flv])

    allocate (sf_test_data_t :: data_strfun)
    select type (data_strfun)
    type is (sf_test_data_t)
       call data_strfun%init (model, pdg_in)
    end select

    allocate (sf_test_spectrum_data_t :: data_spectrum)
    select type (data_spectrum)
    type is (sf_test_spectrum_data_t)
       call data_spectrum%init (model, pdg_in, with_radiation=.true.)
    end select

    write (u, "(A)")  "* Set up chain with beams only"
    write (u, "(A)")

    call sf_chain%init (beam_data)

    call sf_chain_instance%init (sf_chain, n_channel = 1)
    call sf_chain_instance%link_interactions ()
    call sf_chain_instance%exchange_mask ()
    call sf_chain_instance%init_evaluators ()

    call sf_chain_instance%compute_kinematics (1, [real(default) ::])
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    int => sf_chain_instance%get_out_int_ptr ()
    call pset%init (ok, int, int, FM_IGNORE_HELICITY, &
         [0._default, 0._default], .false., .true.)
    call sf_chain_instance%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Particle content:"
    write (u, "(A)")

    call write_separator (u)
    call pset%write (u)
    call write_separator (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover chain:"
    write (u, "(A)")

    call sf_chain_instance%init (sf_chain, n_channel = 1)
    call sf_chain_instance%link_interactions ()
    call sf_chain_instance%exchange_mask ()
    call sf_chain_instance%init_evaluators ()

    int => sf_chain_instance%get_out_int_ptr ()
    call pset%fill_interaction (int, 2, check_match=.false.)

    call sf_chain_instance%recover_kinematics (1)
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    call pset%final ()
    call sf_chain_instance%final ()
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")
    write (u, "(A)")
    write (u, "(A)")  "* Set up chain with structure function"
    write (u, "(A)")

    allocate (sf_config (1))
    call sf_config(1)%init ([1], data_strfun)
    call sf_chain%init (beam_data, sf_config)

    call sf_chain_instance%init (sf_chain, n_channel = 1)
    call sf_channel(1)%init (1)
    call sf_channel(1)%activate_mapping ([1])
    call sf_chain_instance%set_channel (1, sf_channel(1))
    call sf_chain_instance%link_interactions ()
    call sf_chain_instance%exchange_mask ()
    call sf_chain_instance%init_evaluators ()

    call sf_chain_instance%compute_kinematics (1, [0.8_default])
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    int => sf_chain_instance%get_out_int_ptr ()
    call pset%init (ok, int, int, FM_IGNORE_HELICITY, &
         [0._default, 0._default], .false., .true.)
    call sf_chain_instance%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Particle content:"
    write (u, "(A)")

    call write_separator (u)
    call pset%write (u)
    call write_separator (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover chain:"
    write (u, "(A)")

    call sf_chain_instance%init (sf_chain, n_channel = 1)
    call sf_channel(1)%init (1)
    call sf_channel(1)%activate_mapping ([1])
    call sf_chain_instance%set_channel (1, sf_channel(1))
    call sf_chain_instance%link_interactions ()
    call sf_chain_instance%exchange_mask ()
    call sf_chain_instance%init_evaluators ()

    int => sf_chain_instance%get_out_int_ptr ()
    call pset%fill_interaction (int, 2, check_match=.false.)

    call sf_chain_instance%recover_kinematics (1)
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    call pset%final ()
    call sf_chain_instance%final ()
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")
    write (u, "(A)")
    write (u, "(A)")  "* Set up chain with spectrum and structure function"
    write (u, "(A)")

    deallocate (sf_config)
    allocate (sf_config (2))
    call sf_config(1)%init ([1,2], data_spectrum)
    call sf_config(2)%init ([2], data_strfun)
    call sf_chain%init (beam_data, sf_config)

    call sf_chain_instance%init (sf_chain, n_channel = 1)
    call sf_channel(2)%init (2)
    call sf_channel(2)%activate_mapping ([2])
    call sf_chain_instance%set_channel (1, sf_channel(2))
    call sf_chain_instance%link_interactions ()
    call sf_chain_instance%exchange_mask ()
    call sf_chain_instance%init_evaluators ()

    call sf_chain_instance%compute_kinematics &
         (1, [0.5_default, 0.6_default, 0.8_default])
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    int => sf_chain_instance%get_out_int_ptr ()
    call pset%init (ok, int, int, FM_IGNORE_HELICITY, &
         [0._default, 0._default], .false., .true.)
    call sf_chain_instance%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Particle content:"
    write (u, "(A)")

    call write_separator (u)
    call pset%write (u)
    call write_separator (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover chain:"
    write (u, "(A)")

    call sf_chain_instance%init (sf_chain, n_channel = 1)
    call sf_channel(2)%init (2)
    call sf_channel(2)%activate_mapping ([2])
    call sf_chain_instance%set_channel (1, sf_channel(2))
    call sf_chain_instance%link_interactions ()
    call sf_chain_instance%exchange_mask ()
    call sf_chain_instance%init_evaluators ()

    int => sf_chain_instance%get_out_int_ptr ()
    call pset%fill_interaction (int, 2, check_match=.false.)

    call sf_chain_instance%recover_kinematics (1)
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    call pset%final ()
    call sf_chain_instance%final ()
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_11"

  end subroutine sf_base_11

  subroutine sf_base_12 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    type(beam_data_t), target :: beam_data
    class(sf_data_t), allocatable, target :: data
    type(sf_config_t), dimension(:), allocatable, target :: sf_config
    type(sf_chain_t), target :: sf_chain
    type(sf_chain_instance_t), target :: sf_chain_instance
    real(default), dimension(2) :: x_saved
    real(default), dimension(2,3) :: p_saved
    type(sf_channel_t), dimension(:), allocatable :: sf_channel

    write (u, "(A)")  "* Test output: sf_base_12"
    write (u, "(A)")  "*   Purpose: set up and evaluate a multi-channel &
         &structure-function chain"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    call beam_data%init_sqrts (1000._default, [flv, flv])

    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (model, pdg_in)
    end select

    write (u, "(A)")  "* Set up chain with structure function pair &
         &and three different mappings"
    write (u, "(A)")

    allocate (sf_config (2))
    call sf_config(1)%init ([1], data)
    call sf_config(2)%init ([2], data)
    call sf_chain%init (beam_data, sf_config)

    call sf_chain_instance%init (sf_chain, n_channel = 3)

    call allocate_sf_channels (sf_channel, n_channel = 3, n_strfun = 2)

    ! channel 1: no mapping
    call sf_chain_instance%set_channel (1, sf_channel(1))

    ! channel 2: single-particle mappings
    call sf_channel(2)%activate_mapping ([1,2])
    ! call sf_chain_instance%activate_mapping (2, [1,2])
    call sf_chain_instance%set_channel (2, sf_channel(2))

    ! channel 3: two-particle mapping
    call sf_channel(3)%set_s_mapping ([1,2])
    ! call sf_chain_instance%set_s_mapping (3, [1, 2])
    call sf_chain_instance%set_channel (3, sf_channel(3))

    call sf_chain_instance%link_interactions ()
    call sf_chain_instance%exchange_mask ()
    call sf_chain_instance%init_evaluators ()

    write (u, "(A)")  "* Compute kinematics in channel 1 and evaluate"
    write (u, "(A)")

    call sf_chain_instance%compute_kinematics (1, [0.8_default, 0.6_default])
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Invert the kinematics calculation"
    write (u, "(A)")

    x_saved = sf_chain_instance%x

    call sf_chain_instance%inverse_kinematics (x_saved, 1 - x_saved)
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Compute kinematics in channel 2 and evaluate"
    write (u, "(A)")

    p_saved = sf_chain_instance%p

    call sf_chain_instance%compute_kinematics (2, p_saved(:,2))
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Compute kinematics in channel 3 and evaluate"
    write (u, "(A)")

    call sf_chain_instance%compute_kinematics (3, p_saved(:,3))
    call sf_chain_instance%evaluate (scale=0._default)

    call write_separator (u, 2)
    call sf_chain_instance%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_chain_instance%final ()
    call sf_chain%final ()

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_12"

  end subroutine sf_base_12

  subroutine sf_base_13 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t), dimension(2) :: k
    type(vector4_t), dimension(2) :: q
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f, x_free

    write (u, "(A)")  "* Test output: sf_base_13"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &a pair generator object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    allocate (sf_test_generator_data_t :: data)
    select type (data)
    type is (sf_test_generator_data_t)
       call data%init (model, pdg_in)
    end select

    write (u, "(A)")  "* Initialize generator object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    write (u, "(A)")  "* Generate free r values"
    write (u, "(A)")

    x_free = 1
    call sf_int%generate_free (r, rb, x_free)

    write (u, "(A)")  "* Initialize incoming momenta with sqrts=1000"

    E = 500
    k(1) = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    k(2) = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call sf_int%seed_kinematics (k)

    write (u, "(A)")
    write (u, "(A)")  "* Complete kinematics"
    write (u, "(A)")

    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f
    write (u, "(A,9(1x,F10.7))")  "xf=", x_free

    write (u, "(A)")
    write (u, "(A)")  "* Recover x from momenta"
    write (u, "(A)")

    q = sf_int%get_momenta (outgoing=.true.)
    call sf_int%final ()
    deallocate (sf_int)

    call reset_interaction_counter ()
    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)

    call sf_int%seed_kinematics (k)
    call sf_int%set_momenta (q, outgoing=.true.)
    x_free = 1
    call sf_int%recover_x (x, xb, x_free)
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "xf=", x_free

    write (u, "(A)")
    write (u, "(A)")  "* Compute inverse kinematics &
         &and evaluate"
    write (u, "(A)")

    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%apply (scale=0._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "r =", r
    write (u, "(A,9(1x,F10.7))")  "rb=", rb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_13"

  end subroutine sf_base_13

  subroutine sf_base_14 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    type(beam_data_t), target :: beam_data
    class(sf_data_t), allocatable, target :: data_strfun
    class(sf_data_t), allocatable, target :: data_generator
    type(sf_config_t), dimension(:), allocatable, target :: sf_config
    real(default), dimension(:), allocatable :: p_in
    type(sf_chain_t), target :: sf_chain
    type(sf_chain_instance_t), target :: sf_chain_instance

    write (u, "(A)")  "* Test output: sf_base_14"
    write (u, "(A)")  "*   Purpose: set up a structure-function chain"
    write (u, "(A)")  "*            create an instance and evaluate"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)
    pdg_in = 25

    call reset_interaction_counter ()

    call beam_data%init_sqrts (1000._default, [flv, flv])

    allocate (sf_test_data_t :: data_strfun)
    select type (data_strfun)
    type is (sf_test_data_t)
       call data_strfun%init (model, pdg_in)
    end select

    allocate (sf_test_generator_data_t :: data_generator)
    select type (data_generator)
    type is (sf_test_generator_data_t)
       call data_generator%init (model, pdg_in)
    end select

    write (u, "(A)")  "* Set up chain with generator and structure function"
    write (u, "(A)")

    allocate (sf_config (2))
    call sf_config(1)%init ([1,2], data_generator)
    call sf_config(2)%init ([2], data_strfun)
    call sf_chain%init (beam_data, sf_config)

    call sf_chain_instance%init (sf_chain, n_channel = 1)
    call sf_chain_instance%link_interactions ()
    call sf_chain_instance%exchange_mask ()
    call sf_chain_instance%init_evaluators ()

    write (u, "(A)")  "* Inject integration parameter"
    write (u, "(A)")

    allocate (p_in (sf_chain%get_n_bound ()), source = 0.9_default)
    write (u, "(A,9(1x,F10.7))")  "p_in =", p_in

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate"
    write (u, "(A)")

    call sf_chain_instance%compute_kinematics (1, p_in)
    call sf_chain_instance%evaluate (scale=0._default)

    call sf_chain_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Extract integration parameter"
    write (u, "(A)")

    call sf_chain_instance%get_mcpar (1, p_in)
    write (u, "(A,9(1x,F10.7))")  "p_in =", p_in

    call sf_chain_instance%final ()
    call sf_chain%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_base_14"

  end subroutine sf_base_14


  subroutine sf_test_data_write (data, unit, verbose)
    class(sf_test_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "SF test data:"
    write (u, "(3x,A,A)") "model     = ", char (data%model%get_name ())
    write (u, "(3x,A)", advance="no") "incoming  = "
    call data%flv_in%write (u);  write (u, *)
    write (u, "(3x,A)", advance="no") "outgoing  = "
    call data%flv_out%write (u);  write (u, *)
    write (u, "(3x,A)", advance="no") "radiated  = "
    call data%flv_rad%write (u);  write (u, *)
    write (u, "(3x,A," // FMT_19 // ")")  "mass      = ", data%m
    write (u, "(3x,A,L1)")  "collinear = ", data%collinear
    if (.not. data%collinear .and. allocated (data%qbounds)) then
       write (u, "(3x,A," // FMT_19 // ")")  "qmin      = ", data%qbounds(1)
       write (u, "(3x,A," // FMT_19 // ")")  "qmax      = ", data%qbounds(2)
    end if
  end subroutine sf_test_data_write

  subroutine sf_test_data_init (data, model, pdg_in, collinear, qbounds, mode)
    class(sf_test_data_t), intent(out) :: data
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), intent(in) :: pdg_in
    logical, intent(in), optional :: collinear
    real(default), dimension(2), intent(in), optional :: qbounds
    integer, intent(in), optional :: mode
    data%model => model
    if (present (mode))  data%mode = mode
    if (pdg_in%get (1) /= 25) then
       call msg_fatal ("Test spectrum function: input flavor must be 's'")
    end if
    call data%flv_in%init (25, model)
    data%m = data%flv_in%get_mass ()
    if (present (collinear))  data%collinear = collinear
    call data%flv_out%init (25, model)
    call data%flv_rad%init (25, model)
    if (present (qbounds)) then
       allocate (data%qbounds (2))
       data%qbounds = qbounds
    end if
  end subroutine sf_test_data_init

  function sf_test_data_get_n_par (data) result (n)
    class(sf_test_data_t), intent(in) :: data
    integer :: n
    if (data%collinear) then
       n = 1
    else
       n = 3
    end if
  end function sf_test_data_get_n_par

  subroutine sf_test_data_get_pdg_out (data, pdg_out)
    class(sf_test_data_t), intent(in) :: data
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    pdg_out(1) = 25
  end subroutine sf_test_data_get_pdg_out

  subroutine sf_test_data_allocate_sf_int (data, sf_int)
    class(sf_test_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    if (allocated (sf_int)) deallocate (sf_int)
    allocate (sf_test_t :: sf_int)
  end subroutine sf_test_data_allocate_sf_int

  function sf_test_type_string (object) result (string)
    class(sf_test_t), intent(in) :: object
    type(string_t) :: string
    string = "Test"
  end function sf_test_type_string

  subroutine sf_test_write (object, unit, testflag)
    class(sf_test_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%data)) then
       call object%data%write (u)
       call object%base_write (u, testflag)
    else
       write (u, "(1x,A)")  "SF test data: [undefined]"
    end if
  end subroutine sf_test_write

  subroutine sf_test_init (sf_int, data)
    class(sf_test_t), intent(out) :: sf_int
    class(sf_data_t), intent(in), target :: data
    type(quantum_numbers_mask_t), dimension(3) :: mask
    type(helicity_t) :: hel0
    type(color_t) :: col0
    type(quantum_numbers_t), dimension(3) :: qn
    mask = quantum_numbers_mask (.false., .false., .false.)
    select type (data)
    type is (sf_test_data_t)
       if (allocated (data%qbounds)) then
          call sf_int%base_init (mask, &
               [data%m**2], [0._default], [data%m**2], &
               [data%qbounds(1)], [data%qbounds(2)])
       else
          call sf_int%base_init (mask, &
               [data%m**2], [0._default], [data%m**2])
       end if
       sf_int%data => data
       call hel0%init (0)
       call col0%init ()
       call qn(1)%init (data%flv_in,  col0, hel0)
       call qn(2)%init (data%flv_rad, col0, hel0)
       call qn(3)%init (data%flv_out, col0, hel0)
       call sf_int%add_state (qn)
       call sf_int%freeze ()
       call sf_int%set_incoming ([1])
       call sf_int%set_radiated ([2])
       call sf_int%set_outgoing ([3])
    end select
    sf_int%status = SF_INITIAL
  end subroutine sf_test_init

  subroutine sf_test_complete_kinematics (sf_int, x, xb, f, r, rb, map)
    class(sf_test_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: r
    real(default), dimension(:), intent(in) :: rb
    logical, intent(in) :: map
    if (map) then
       x(1) = r(1)**2
       f = 2 * r(1)
    else
       x(1) = r(1)
       f = 1
    end if
    xb(1) = 1 - x(1)
    if (size (x) == 3) then
       x(2:3)  = r(2:3)
       xb(2:3) = rb(2:3)
    end if
    call sf_int%split_momentum (x, xb)
    sf_int%x = x(1)
    select case (sf_int%status)
    case (SF_FAILED_KINEMATICS);  f = 0
    end select
  end subroutine sf_test_complete_kinematics

  subroutine sf_test_inverse_kinematics (sf_int, x, xb, f, r, rb, map, set_momenta)
    class(sf_test_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: r
    real(default), dimension(:), intent(out) :: rb
    logical, intent(in) :: map
    logical, intent(in), optional :: set_momenta
    logical :: set_mom
    set_mom = .false.;  if (present (set_momenta))  set_mom = set_momenta
    if (map) then
       r(1) = sqrt (x(1))
       f = 2 * r(1)
    else
       r(1) = x(1)
       f = 1
    end if
    if (size (x) == 3)  r(2:3) = x(2:3)
    rb = 1 - r
    sf_int%x = x(1)
    if (set_mom) then
       call sf_int%split_momentum (x, xb)
       select case (sf_int%status)
       case (SF_FAILED_KINEMATICS);  f = 0
       end select
    end if
  end subroutine sf_test_inverse_kinematics

  subroutine sf_test_apply (sf_int, scale, negative_sf, rescale, i_sub)
    class(sf_test_t), intent(inout) :: sf_int
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(in), optional :: rescale
    integer, intent(in), optional :: i_sub
    select case (sf_int%data%mode)
    case (0)
       call sf_int%set_matrix_element &
            (cmplx (1._default, kind=default))
    case (1)
       call sf_int%set_matrix_element &
            (cmplx (sf_int%x, kind=default))
    end select
    sf_int%status = SF_EVALUATED
  end subroutine sf_test_apply

  subroutine sf_test_spectrum_data_write (data, unit, verbose)
    class(sf_test_spectrum_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "SF test spectrum data:"
    write (u, "(3x,A,A)") "model     = ", char (data%model%get_name ())
    write (u, "(3x,A)", advance="no") "incoming  = "
    call data%flv_in%write (u);  write (u, *)
    write (u, "(3x,A)", advance="no") "outgoing  = "
    call data%flv_out%write (u);  write (u, *)
    write (u, "(3x,A)", advance="no") "radiated  = "
    call data%flv_rad%write (u);  write (u, *)
    write (u, "(3x,A," // FMT_19 // ")")  "mass      = ", data%m
  end subroutine sf_test_spectrum_data_write

  subroutine sf_test_spectrum_data_init (data, model, pdg_in, with_radiation)
    class(sf_test_spectrum_data_t), intent(out) :: data
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), intent(in) :: pdg_in
    logical, intent(in) :: with_radiation
    data%model => model
    data%with_radiation = with_radiation
    if (pdg_in%get (1) /= 25) then
       call msg_fatal ("Test structure function: input flavor must be 's'")
    end if
    call data%flv_in%init (25, model)
    data%m = data%flv_in%get_mass ()
    call data%flv_out%init (25, model)
    if (with_radiation) then
       call data%flv_rad%init (25, model)
    end if
  end subroutine sf_test_spectrum_data_init

  function sf_test_spectrum_data_get_n_par (data) result (n)
    class(sf_test_spectrum_data_t), intent(in) :: data
    integer :: n
    n = 2
  end function sf_test_spectrum_data_get_n_par

  subroutine sf_test_spectrum_data_get_pdg_out (data, pdg_out)
    class(sf_test_spectrum_data_t), intent(in) :: data
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    pdg_out(1) = 25
    pdg_out(2) = 25
  end subroutine sf_test_spectrum_data_get_pdg_out

  subroutine sf_test_spectrum_data_allocate_sf_int (data, sf_int)
    class(sf_test_spectrum_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (sf_test_spectrum_t :: sf_int)
  end subroutine sf_test_spectrum_data_allocate_sf_int

  function sf_test_spectrum_type_string (object) result (string)
    class(sf_test_spectrum_t), intent(in) :: object
    type(string_t) :: string
    string = "Test Spectrum"
  end function sf_test_spectrum_type_string

  subroutine sf_test_spectrum_write (object, unit, testflag)
    class(sf_test_spectrum_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%data)) then
       call object%data%write (u)
       call object%base_write (u, testflag)
    else
       write (u, "(1x,A)")  "SF test spectrum data: [undefined]"
    end if
  end subroutine sf_test_spectrum_write

  subroutine sf_test_spectrum_init (sf_int, data)
    class(sf_test_spectrum_t), intent(out) :: sf_int
    class(sf_data_t), intent(in), target :: data
    type(quantum_numbers_mask_t), dimension(6) :: mask
    type(helicity_t) :: hel0
    type(color_t) :: col0
    type(quantum_numbers_t), dimension(6) :: qn
    mask = quantum_numbers_mask (.false., .false., .false.)
    select type (data)
    type is (sf_test_spectrum_data_t)
       if (data%with_radiation) then
          call sf_int%base_init (mask(1:6), &
               [data%m**2, data%m**2], &
               [0._default, 0._default], &
               [data%m**2, data%m**2])
          sf_int%data => data
          call hel0%init (0)
          call col0%init ()
          call qn(1)%init (data%flv_in,  col0, hel0)
          call qn(2)%init (data%flv_in,  col0, hel0)
          call qn(3)%init (data%flv_rad, col0, hel0)
          call qn(4)%init (data%flv_rad, col0, hel0)
          call qn(5)%init (data%flv_out, col0, hel0)
          call qn(6)%init (data%flv_out, col0, hel0)
          call sf_int%add_state (qn(1:6))
          call sf_int%set_incoming ([1,2])
          call sf_int%set_radiated ([3,4])
          call sf_int%set_outgoing ([5,6])
       else
          call sf_int%base_init (mask(1:4), &
               [data%m**2, data%m**2], &
               [real(default) :: ], &
               [data%m**2, data%m**2])
          sf_int%data => data
          call hel0%init (0)
          call col0%init ()
          call qn(1)%init (data%flv_in,  col0, hel0)
          call qn(2)%init (data%flv_in,  col0, hel0)
          call qn(3)%init (data%flv_out, col0, hel0)
          call qn(4)%init (data%flv_out, col0, hel0)
          call sf_int%add_state (qn(1:4))
          call sf_int%set_incoming ([1,2])
          call sf_int%set_outgoing ([3,4])
       end if
       call sf_int%freeze ()
    end select
    sf_int%status = SF_INITIAL
  end subroutine sf_test_spectrum_init

  subroutine sf_test_spectrum_complete_kinematics (sf_int, x, xb, f, r, rb, map)
    class(sf_test_spectrum_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: r
    real(default), dimension(:), intent(in) :: rb
    logical, intent(in) :: map
    real(default), dimension(2) :: xb1
    if (map) then
       x = r**2
       f = 4 * r(1) * r(2)
    else
       x = r
       f = 1
    end if
    xb = 1 - x
    if (sf_int%data%with_radiation) then
       call sf_int%split_momenta (x, xb)
    else
       call sf_int%reduce_momenta (x)
    end if
    select case (sf_int%status)
    case (SF_FAILED_KINEMATICS);  f = 0
    end select
  end subroutine sf_test_spectrum_complete_kinematics

  subroutine sf_test_spectrum_inverse_kinematics &
       (sf_int, x, xb, f, r, rb, map, set_momenta)
    class(sf_test_spectrum_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
     real(default), dimension(:), intent(in) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: r
    real(default), dimension(:), intent(out) :: rb
    logical, intent(in) :: map
    logical, intent(in), optional :: set_momenta
    real(default), dimension(2) :: xb1
    logical :: set_mom
    set_mom = .false.;  if (present (set_momenta))  set_mom = set_momenta
    if (map) then
       r = sqrt (x)
       f = 4 * r(1) * r(2)
    else
       r = x
       f = 1
    end if
    rb = 1 - r
    if (set_mom)  then
       if (sf_int%data%with_radiation) then
          call sf_int%split_momenta (x, xb)
       else
          call sf_int%reduce_momenta (x)
       end if
       select case (sf_int%status)
       case (SF_FAILED_KINEMATICS);  f = 0
       end select
    end if
  end subroutine sf_test_spectrum_inverse_kinematics

  subroutine sf_test_spectrum_apply (sf_int, scale, negative_sf, rescale, i_sub)
    class(sf_test_spectrum_t), intent(inout) :: sf_int
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(in), optional :: rescale
    integer, intent(in), optional :: i_sub
    call sf_int%set_matrix_element &
         (cmplx (1._default, kind=default))
    sf_int%status = SF_EVALUATED
  end subroutine sf_test_spectrum_apply

  subroutine sf_test_generator_data_write (data, unit, verbose)
    class(sf_test_generator_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "SF test generator data:"
    write (u, "(3x,A,A)") "model     = ", char (data%model%get_name ())
    write (u, "(3x,A)", advance="no") "incoming  = "
    call data%flv_in%write (u);  write (u, *)
    write (u, "(3x,A)", advance="no") "outgoing  = "
    call data%flv_out%write (u);  write (u, *)
    write (u, "(3x,A," // FMT_19 // ")")  "mass      = ", data%m
  end subroutine sf_test_generator_data_write

  subroutine sf_test_generator_data_init (data, model, pdg_in)
    class(sf_test_generator_data_t), intent(out) :: data
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), intent(in) :: pdg_in
    data%model => model
    if (pdg_in%get (1) /= 25) then
       call msg_fatal ("Test generator: input flavor must be 's'")
    end if
    call data%flv_in%init (25, model)
    data%m = data%flv_in%get_mass ()
    call data%flv_out%init (25, model)
  end subroutine sf_test_generator_data_init

  function sf_test_generator_data_is_generator (data) result (flag)
    class(sf_test_generator_data_t), intent(in) :: data
    logical :: flag
    flag = .true.
  end function sf_test_generator_data_is_generator

  function sf_test_generator_data_get_n_par (data) result (n)
    class(sf_test_generator_data_t), intent(in) :: data
    integer :: n
    n = 2
  end function sf_test_generator_data_get_n_par

  subroutine sf_test_generator_data_get_pdg_out (data, pdg_out)
    class(sf_test_generator_data_t), intent(in) :: data
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    pdg_out(1) = 25
    pdg_out(2) = 25
  end subroutine sf_test_generator_data_get_pdg_out

  subroutine sf_test_generator_data_allocate_sf_int (data, sf_int)
    class(sf_test_generator_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (sf_test_generator_t :: sf_int)
  end subroutine sf_test_generator_data_allocate_sf_int

  function sf_test_generator_type_string (object) result (string)
    class(sf_test_generator_t), intent(in) :: object
    type(string_t) :: string
    string = "Test Generator"
  end function sf_test_generator_type_string

  subroutine sf_test_generator_write (object, unit, testflag)
    class(sf_test_generator_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%data)) then
       call object%data%write (u)
       call object%base_write (u, testflag)
    else
       write (u, "(1x,A)")  "SF test generator data: [undefined]"
    end if
  end subroutine sf_test_generator_write

  subroutine sf_test_generator_init (sf_int, data)
    class(sf_test_generator_t), intent(out) :: sf_int
    class(sf_data_t), intent(in), target :: data
    type(quantum_numbers_mask_t), dimension(4) :: mask
    type(helicity_t) :: hel0
    type(color_t) :: col0
    type(quantum_numbers_t), dimension(4) :: qn
    mask = quantum_numbers_mask (.false., .false., .false.)
    select type (data)
    type is (sf_test_generator_data_t)
       call sf_int%base_init (mask(1:4), &
            [data%m**2, data%m**2], &
            [real(default) :: ], &
            [data%m**2, data%m**2])
       sf_int%data => data
       call hel0%init (0)
       call col0%init ()
       call qn(1)%init (data%flv_in,  col0, hel0)
       call qn(2)%init (data%flv_in,  col0, hel0)
       call qn(3)%init (data%flv_out, col0, hel0)
       call qn(4)%init (data%flv_out, col0, hel0)
       call sf_int%add_state (qn(1:4))
       call sf_int%set_incoming ([1,2])
       call sf_int%set_outgoing ([3,4])
       call sf_int%freeze ()
    end select
    sf_int%status = SF_INITIAL
  end subroutine sf_test_generator_init

  function sf_test_generator_is_generator (sf_int) result (flag)
    class(sf_test_generator_t), intent(in) :: sf_int
    logical :: flag
    flag = sf_int%data%is_generator ()
  end function sf_test_generator_is_generator

  subroutine sf_test_generator_generate_free (sf_int, r, rb,  x_free)
    class(sf_test_generator_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(inout) :: x_free
    r = [0.8, 0.5]
    rb= 1 - r
    x_free = x_free * product (r)
  end subroutine sf_test_generator_generate_free

  subroutine sf_test_generator_recover_x (sf_int, x, xb, x_free)
    class(sf_test_generator_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(inout), optional :: x_free
    call sf_int%base_recover_x (x, xb)
    if (present (x_free))  x_free = x_free * product (x)
  end subroutine sf_test_generator_recover_x

  subroutine sf_test_generator_complete_kinematics (sf_int, x, xb, f, r, rb, map)
    class(sf_test_generator_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: r
    real(default), dimension(:), intent(in) :: rb
    logical, intent(in) :: map
    x = r
    xb= rb
    f = 1
    call sf_int%reduce_momenta (x)
  end subroutine sf_test_generator_complete_kinematics

  subroutine sf_test_generator_inverse_kinematics &
       (sf_int, x, xb, f, r, rb, map, set_momenta)
    class(sf_test_generator_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: r
    real(default), dimension(:), intent(out) :: rb
    logical, intent(in) :: map
    logical, intent(in), optional :: set_momenta
    logical :: set_mom
    set_mom = .false.;  if (present (set_momenta))  set_mom = set_momenta
    r = x
    rb= xb
    f = 1
    if (set_mom)  call sf_int%reduce_momenta (x)
  end subroutine sf_test_generator_inverse_kinematics

  subroutine sf_test_generator_apply (sf_int, scale, negative_sf, rescale, i_sub)
    class(sf_test_generator_t), intent(inout) :: sf_int
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(in), optional :: rescale
    integer, intent(in), optional :: i_sub
    call sf_int%set_matrix_element &
         (cmplx (1._default, kind=default))
    sf_int%status = SF_EVALUATED
  end subroutine sf_test_generator_apply


end module sf_base_uti
