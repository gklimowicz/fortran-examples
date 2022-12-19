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

module sf_isr_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use format_defs, only: FMT_12
  use physics_defs, only: ELECTRON
  use lorentz
  use pdg_arrays
  use flavors
  use interactions, only: reset_interaction_counter
  use interactions, only: interaction_t
  use model_data
  use sf_aux, only: KEEP_ENERGY
  use sf_mappings
  use sf_base

  use sf_isr

  implicit none
  private

  public :: sf_isr_1
  public :: sf_isr_2
  public :: sf_isr_3
  public :: sf_isr_4
  public :: sf_isr_5

contains

  subroutine sf_isr_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(pdg_array_t), dimension(1) :: pdg_out
    integer, dimension(:), allocatable :: pdg1
    class(sf_data_t), allocatable :: data

    write (u, "(A)")  "* Test output: sf_isr_1"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &test structure function data"
    write (u, "(A)")

    write (u, "(A)")  "* Create empty data object"
    write (u, "(A)")

    call model%init_qed_test ()
    pdg_in = ELECTRON

    allocate (isr_data_t :: data)
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize"
    write (u, "(A)")

    select type (data)
    type is (isr_data_t)
       call data%init (model, pdg_in, 1./137._default, 10._default, &
            0.000511_default, order = 3, recoil = .false.)
    end select

    call data%write (u)

    write (u, "(A)")

    write (u, "(1x,A)")  "Outgoing particle codes:"
    call data%get_pdg_out (pdg_out)
    pdg1 = pdg_out(1)
    write (u, "(2x,99(1x,I0))")  pdg1

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_isr_1"

  end subroutine sf_isr_1

  subroutine sf_isr_2 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(flavor_t) :: flv
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f, f_isr

    write (u, "(A)")  "* Test output: sf_isr_2"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_qed_test ()
    pdg_in = ELECTRON
    call flv%init (ELECTRON, model)

    call reset_interaction_counter ()

    allocate (isr_data_t :: data)
    select type (data)
    type is (isr_data_t)
       call data%init (model, pdg_in, 1./137._default, 500._default, &
            0.000511_default, order = 3, recoil = .false.)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])

    write (u, "(A)")  "* Initialize incoming momentum with E=500"
    write (u, "(A)")
    E = 500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call pacify (k, 1e-10_default)
    call vector4_write (k, u)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for r=0.9, no ISR mapping, &
         &collinear"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = 0.9_default
    rb = 1 - r
    write (u, "(A,9(1x," // FMT_12 // "))")  "r =", r
    write (u, "(A,9(1x," // FMT_12 // "))")  "rb=", rb

    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)

    write (u, "(A)")
    write (u, "(A,9(1x," // FMT_12 // "))")  "x =", x
    write (u, "(A,9(1x," // FMT_12 // "))")  "xb=", xb
    write (u, "(A,9(1x," // FMT_12 // "))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Invert kinematics"
    write (u, "(A)")

    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.false.)
    write (u, "(A,9(1x," // FMT_12 // "))")  "r =", r
    write (u, "(A,9(1x," // FMT_12 // "))")  "rb=", rb
    write (u, "(A,9(1x," // FMT_12 // "))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate ISR structure function"
    write (u, "(A)")

    call sf_int%apply (scale = 100._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Structure-function value, default order"
    write (u, "(A)")

    f_isr = sf_int%get_matrix_element (1)

    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr         =", f_isr
    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr * f_map =", f_isr * f

    write (u, "(A)")
    write (u, "(A)")  "* Re-evaluate structure function, leading order"
    write (u, "(A)")

    select type (sf_int)
    type is (isr_t)
       call sf_int%set_order (0)
    end select
    call sf_int%apply (scale = 100._default)
    f_isr = sf_int%get_matrix_element (1)

    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr         =", f_isr
    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr * f_map =", f_isr * f

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_isr_2"

  end subroutine sf_isr_2

  subroutine sf_isr_3 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f, f_isr

    write (u, "(A)")  "* Test output: sf_isr_3"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_qed_test ()
    call flv%init (ELECTRON, model)
    pdg_in = ELECTRON

    call reset_interaction_counter ()

    allocate (isr_data_t :: data)
    select type (data)
    type is (isr_data_t)
       call data%init (model, pdg_in, 1./137._default, 500._default, &
            0.000511_default, order = 3, recoil = .false.)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])

    write (u, "(A)")  "* Initialize incoming momentum with E=500"
    write (u, "(A)")
    E = 500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call pacify (k, 1e-10_default)
    call vector4_write (k, u)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for r=0.7, with ISR mapping, &
         &collinear"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = 0.7_default
    rb = 1 - r
    write (u, "(A,9(1x," // FMT_12 // "))")  "r =", r
    write (u, "(A,9(1x," // FMT_12 // "))")  "rb=", rb

    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.true.)

    write (u, "(A)")
    write (u, "(A,9(1x," // FMT_12 // "))")  "x =", x
    write (u, "(A,9(1x," // FMT_12 // "))")  "xb=", xb
    write (u, "(A,9(1x," // FMT_12 // "))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Invert kinematics"
    write (u, "(A)")

    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.true.)
    write (u, "(A,9(1x," // FMT_12 // "))")  "r =", r
    write (u, "(A,9(1x," // FMT_12 // "))")  "rb=", rb
    write (u, "(A,9(1x," // FMT_12 // "))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate ISR structure function"
    write (u, "(A)")

    call sf_int%apply (scale = 100._default)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Structure-function value, default order"
    write (u, "(A)")

    f_isr = sf_int%get_matrix_element (1)

    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr         =", f_isr
    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr * f_map =", f_isr * f

    write (u, "(A)")
    write (u, "(A)")  "* Re-evaluate structure function, leading order"
    write (u, "(A)")

    select type (sf_int)
    type is (isr_t)
       call sf_int%set_order (0)
    end select
    call sf_int%apply (scale = 100._default)
    f_isr = sf_int%get_matrix_element (1)

    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr         =", f_isr
    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr * f_map =", f_isr * f

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_isr_3"

  end subroutine sf_isr_3

  subroutine sf_isr_4 (u)
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
    real(default) :: f, f_isr
    character(len=80) :: buffer
    integer :: u_scratch, iostat

    write (u, "(A)")  "* Test output: sf_isr_4"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_qed_test ()
    call flv%init (ELECTRON, model)
    pdg_in = ELECTRON

    call reset_interaction_counter ()

    write (u, "(A)")
    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    allocate (isr_data_t :: data)
    select type (data)
    type is (isr_data_t)
       call data%init (model, pdg_in, 1./137._default, 500._default, &
            0.000511_default, order = 3, recoil = .true.)
    end select

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])

    write (u, "(A)")
    write (u, "(A)")  "* Initialize incoming momentum with E=500"
    write (u, "(A)")
    E = 500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call pacify (k, 1e-10_default)
    call vector4_write (k, u)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5/0.5/0.25, with ISR mapping, "
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

    call sf_int%seed_kinematics ([k])
    call sf_int%set_momenta (q, outgoing=.true.)
    call sf_int%recover_x (x, xb)
    call sf_int%inverse_kinematics (x, xb, f, r, rb, map=.true.)

    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "r =", r

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate ISR structure function"
    write (u, "(A)")

    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.true.)
    call sf_int%pacify_momenta (1e-10_default)
    call sf_int%apply (scale = 10._default)
    u_scratch = free_unit ()
    open (u_scratch, status="scratch", action = "readwrite")
    call sf_int%write (u_scratch, testflag = .true.)
    rewind (u_scratch)
    do
       read (u_scratch, "(A)", iostat=iostat) buffer
       if (iostat /= 0) exit
       if (buffer(1:25) == " P =   0.000000E+00  9.57") then
          buffer = replace (buffer, 26, "XXXX")
       end if
       if (buffer(1:25) == " P =   0.000000E+00 -9.57") then
          buffer = replace (buffer, 26, "XXXX")
       end if
       write (u, "(A)") buffer
    end do
    close (u_scratch)

    write (u, "(A)")
    write (u, "(A)")  "* Structure-function value"
    write (u, "(A)")

    f_isr = sf_int%get_matrix_element (1)

    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr         =", f_isr
    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr * f_map =", f_isr * f

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_isr_4"

  end subroutine sf_isr_4

  subroutine sf_isr_5 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    class(sf_mapping_t), allocatable :: mapping
    class(sf_int_t), dimension(:), allocatable :: sf_int
    type(vector4_t), dimension(2) :: k
    real(default) :: E, f_map
    real(default), dimension(:), allocatable :: p, pb, r, rb, x, xb
    real(default), dimension(2) :: f, f_isr
    integer :: i

    write (u, "(A)")  "* Test output: sf_isr_5"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_qed_test ()
    call flv%init (ELECTRON, model)
    pdg_in = ELECTRON

    call reset_interaction_counter ()

    allocate (isr_data_t :: data)
    select type (data)
    type is (isr_data_t)
       call data%init (model, pdg_in, 1./137._default, 500._default, &
            0.000511_default, order = 3, recoil = .false.)
    end select

    allocate (sf_ip_mapping_t :: mapping)
    select type (mapping)
    type is (sf_ip_mapping_t)
       select type (data)
       type is (isr_data_t)
          call mapping%init (eps = data%get_eps ())
       end select
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    allocate (isr_t :: sf_int (2))

    do i = 1, 2
       call sf_int(i)%init (data)
       call sf_int(i)%set_beam_index ([i])
    end do

    write (u, "(A)")  "* Initialize incoming momenta with E=500"
    write (u, "(A)")
    E = 500
    k(1) = vector4_moving (E,   sqrt (E**2 - flv%get_mass ()**2), 3)
    k(2) = vector4_moving (E, - sqrt (E**2 - flv%get_mass ()**2), 3)
    call pacify (k, 1e-10_default)
    do i = 1, 2
       call vector4_write (k(i), u)
       call sf_int(i)%seed_kinematics (k(i:i))
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for p=[0.7,0.4], collinear"
    write (u, "(A)")

    allocate (p (2 * data%get_n_par ()))
    allocate (pb(size (p)))
    allocate (r (size (p)))
    allocate (rb(size (p)))
    allocate (x (size (p)))
    allocate (xb(size (p)))

    p = [0.7_default, 0.4_default]
    pb= 1 - p
    call mapping%compute (r, rb, f_map, p, pb)

    write (u, "(A,9(1x," // FMT_12 // "))")  "p =", p
    write (u, "(A,9(1x," // FMT_12 // "))")  "pb=", pb
    write (u, "(A,9(1x," // FMT_12 // "))")  "r =", r
    write (u, "(A,9(1x," // FMT_12 // "))")  "rb=", rb
    write (u, "(A,9(1x," // FMT_12 // "))")  "fm=", f_map

    do i = 1, 2
       call sf_int(i)%complete_kinematics (x(i:i), xb(i:i), f(i), r(i:i), rb(i:i), &
            map=.false.)
    end do

    write (u, "(A)")
    write (u, "(A,9(1x," // FMT_12 // "))")  "x =", x
    write (u, "(A,9(1x," // FMT_12 // "))")  "xb=", xb
    write (u, "(A,9(1x," // FMT_12 // "))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Invert kinematics"
    write (u, "(A)")

    do i = 1, 2
       call sf_int(i)%inverse_kinematics (x(i:i), xb(i:i), f(i), r(i:i), rb(i:i), &
            map=.false.)
    end do
    call mapping%inverse (r, rb, f_map, p, pb)

    write (u, "(A,9(1x," // FMT_12 // "))")  "p =", p
    write (u, "(A,9(1x," // FMT_12 // "))")  "pb=", pb
    write (u, "(A,9(1x," // FMT_12 // "))")  "r =", r
    write (u, "(A,9(1x," // FMT_12 // "))")  "rb=", rb
    write (u, "(A,9(1x," // FMT_12 // "))")  "fm=", f_map

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate ISR structure function"

    call sf_int(1)%apply (scale = 100._default)
    call sf_int(2)%apply (scale = 100._default)

    write (u, "(A)")
    write (u, "(A)")  "* Structure function #1"
    write (u, "(A)")
    call sf_int(1)%write (u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Structure function #2"
    write (u, "(A)")
    call sf_int(2)%write (u, testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Structure-function value, default order"
    write (u, "(A)")

    do i = 1, 2
       f_isr(i) = sf_int(i)%get_matrix_element (1)
    end do

    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr         =", &
         product (f_isr)
    write (u, "(A,9(1x," // FMT_12 // "))")  "f_isr * f_map =", &
         product (f_isr * f) * f_map

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    do i = 1, 2
       call sf_int(i)%final ()
    end do
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_isr_5"

  end subroutine sf_isr_5


end module sf_isr_uti
