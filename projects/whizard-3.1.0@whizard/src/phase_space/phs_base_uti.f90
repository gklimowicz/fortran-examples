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

module phs_base_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use diagnostics
  use io_units
  use format_defs, only: FMT_19
  use physics_defs, only: BORN
  use lorentz
  use flavors
  use model_data
  use process_constants

  use phs_base

  implicit none
  private

  public :: init_test_process_data
  public :: init_test_decay_data
  public :: phs_test_config_t
  public :: phs_test_t

  public :: phs_base_1
  public :: phs_base_2
  public :: phs_base_3
  public :: phs_base_4
  public :: phs_base_5

  type, extends (phs_config_t) :: phs_test_config_t
     logical :: create_equivalences = .false.
   contains
     procedure :: final => phs_test_config_final
     procedure :: write => phs_test_config_write
     procedure :: configure => phs_test_config_configure
     procedure :: startup_message => phs_test_config_startup_message
     procedure, nopass :: allocate_instance => phs_test_config_allocate_instance
  end type phs_test_config_t

  type, extends (phs_t) :: phs_test_t
     real(default) :: m = 0
     real(default), dimension(:), allocatable :: x
   contains
     procedure :: write => phs_test_write
     procedure :: final => phs_test_final
     procedure :: init => phs_test_init
     procedure :: evaluate_selected_channel => phs_test_evaluate_selected_channel
     procedure :: evaluate_other_channels => phs_test_evaluate_other_channels
     procedure :: inverse => phs_test_inverse
  end type phs_test_t


contains

  subroutine phs_base_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    class(phs_config_t), allocatable :: phs_data

    write (u, "(A)")  "* Test output: phs_base_1"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &test phase-space configuration data"
    write (u, "(A)")

    call model%init_test ()

    write (u, "(A)")  "* Initialize a process and a matching &
         &phase-space configuration"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_base_1"), process_data)

    allocate (phs_test_config_t :: phs_data)
    call phs_data%init (process_data, model)

    call phs_data%write (u)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_base_1"

  end subroutine phs_base_1

  subroutine phs_base_2 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(process_constants_t) :: process_data
    real(default) :: sqrts, E
    class(phs_config_t), allocatable, target :: phs_data
    class(phs_t), pointer :: phs => null ()
    type(vector4_t), dimension(2) :: p, q

    write (u, "(A)")  "* Test output: phs_base_2"
    write (u, "(A)")  "*   Purpose: test simple two-channel phase space"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)

    write (u, "(A)")  "* Initialize a process and a matching &
         &phase-space configuration"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_base_2"), process_data)

    allocate (phs_test_config_t :: phs_data)
    call phs_data%init (process_data, model)

    sqrts = 1000._default
    call phs_data%configure (sqrts)

    call phs_data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize the phase-space instance"
    write (u, "(A)")

    call phs_data%allocate_instance (phs)
    select type (phs)
    type is (phs_test_t)
       call phs%init (phs_data)
    end select

    call phs%write (u, verbose=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Set incoming momenta"
    write (u, "(A)")

    E = sqrts / 2
    p(1) = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    p(2) = vector4_moving (E,-sqrt (E**2 - flv%get_mass ()**2), 3)

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute phase-space point in channel 1 &
         &for x = 0.5, 0.125"
    write (u, "(A)")

    call phs%evaluate_selected_channel (1, [0.5_default, 0.125_default])
    call phs%evaluate_other_channels (1)
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute phase-space point in channel 2 &
         &for x = 0.125, 0.125"
    write (u, "(A)")

    call phs%evaluate_selected_channel (2, [0.125_default, 0.125_default])
    call phs%evaluate_other_channels (2)
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Inverse kinematics"
    write (u, "(A)")

    call phs%get_outgoing_momenta (q)
    deallocate (phs)
    call phs_data%allocate_instance (phs)
    call phs%init (phs_data)

    sqrts = 1000._default
    select type (phs_data)
    type is (phs_test_config_t)
       call phs_data%configure (sqrts)
    end select

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()
    call phs%set_outgoing_momenta (q)

    call phs%inverse ()
    call phs%write (u)

    call phs%final ()
    deallocate (phs)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_base_2"

  end subroutine phs_base_2

  subroutine phs_base_3 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    class(phs_config_t), allocatable :: phs_data

    write (u, "(A)")  "* Test output: phs_base_3"
    write (u, "(A)")  "*   Purpose: construct phase-space configuration data &
         &with equivalences"
    write (u, "(A)")

    call model%init_test ()

    write (u, "(A)")  "* Initialize a process and a matching &
         &phase-space configuration"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_base_3"), process_data)

    allocate (phs_test_config_t :: phs_data)
    call phs_data%init (process_data, model)
    select type (phs_data)
    type is (phs_test_config_t)
       phs_data%create_equivalences = .true.
    end select

    call phs_data%configure (1000._default)
    call phs_data%write (u)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_base_3"

  end subroutine phs_base_3

  subroutine phs_base_4 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    class(phs_config_t), allocatable :: phs_data

    write (u, "(A)")  "* Test output: phs_base_4"
    write (u, "(A)")  "*   Purpose: compute and compare MD5 sums"
    write (u, "(A)")

    call model%init_test ()

    write (u, "(A)")  "* Model parameters"
    write (u, "(A)")

    call model%write (unit = u, &
         show_parameters = .true., &
         show_particles = .false., show_vertices = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize a process and a matching &
         &phase-space configuration"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_base_4"), process_data)
    process_data%md5sum = "test_process_data_m6sum_12345678"

    allocate (phs_test_config_t :: phs_data)
    call phs_data%init (process_data, model)

    call phs_data%compute_md5sum ()
    call phs_data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Modify model parameter"
    write (u, "(A)")

    call model%set_par (var_str ("ms"), 100._default)
    call model%write (show_parameters = .true., &
         show_particles = .false., show_vertices = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* PHS configuration"
    write (u, "(A)")

    call phs_data%compute_md5sum ()
    call phs_data%write (u)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_base_4"

  end subroutine phs_base_4

  subroutine phs_base_5 (u)
    integer, intent(in) :: u
    type(phs_channel_t), dimension(:), allocatable :: channel
    type(phs_channel_collection_t) :: coll
    integer :: i, n

    write (u, "(A)")  "* Test output: phs_base_5"
    write (u, "(A)")  "*   Purpose: collect channel properties"
    write (u, "(A)")

    write (u, "(A)")  "* Set up an array of channels"
    write (u, "(A)")

    n = 6

    allocate (channel (n))
    call channel(2)%set_resonant (75._default, 3._default)
    call channel(4)%set_resonant (130._default, 1._default)
    call channel(5)%set_resonant (75._default, 3._default)
    call channel(6)%set_on_shell (33._default)

    do i = 1, n
       write (u, "(1x,I0)", advance="no")  i
       call channel(i)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Collect distinct properties"
    write (u, "(A)")

    do i = 1, n
       call coll%push (channel(i))
    end do

    write (u, "(1x,A,I0)")  "n = ", coll%get_n ()
    write (u, "(A)")

    call coll%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Channel array with collection index assigned"
    write (u, "(A)")

    do i = 1, n
       write (u, "(1x,I0)", advance="no")  i
       call channel(i)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call coll%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_base_5"

  end subroutine phs_base_5


  subroutine init_test_process_data (id, data)
    type(process_constants_t), intent(out) :: data
    type(string_t), intent(in), optional :: id
    if (present (id)) then
       data%id = id
    else
       data%id = "testproc"
    end if
    data%model_name = "Test"
    data%n_in = 2
    data%n_out = 2
    data%n_flv = 1
    allocate (data%flv_state (data%n_in + data%n_out, data%n_flv))
    data%flv_state = 25
  end subroutine init_test_process_data

  subroutine init_test_decay_data (id, data)
    type(process_constants_t), intent(out) :: data
    type(string_t), intent(in), optional :: id
    if (present (id)) then
       data%id = id
    else
       data%id = "testproc"
    end if
    data%model_name = "Test"
    data%n_in = 1
    data%n_out = 2
    data%n_flv = 1
    allocate (data%flv_state (data%n_in + data%n_out, data%n_flv))
    data%flv_state(:,1) = [25, 6, -6]
  end subroutine init_test_decay_data

  subroutine phs_test_config_final (object)
    class(phs_test_config_t), intent(inout) :: object
  end subroutine phs_test_config_final

  subroutine phs_test_config_write (object, unit, include_id)
    class(phs_test_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: include_id
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Partonic phase-space configuration:"
    call object%base_write (unit)
  end subroutine phs_test_config_write

  subroutine phs_test_config_configure (phs_config, sqrts, &
       sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
       ignore_mismatch, nlo_type, subdir)
    class(phs_test_config_t), intent(inout) :: phs_config
    real(default), intent(in) :: sqrts
    logical, intent(in), optional :: sqrts_fixed
    logical, intent(in), optional :: lab_is_cm
    logical, intent(in), optional :: azimuthal_dependence
    logical, intent(in), optional :: rebuild
    logical, intent(in), optional :: ignore_mismatch
    integer, intent(in), optional :: nlo_type
    type(string_t), intent(in), optional :: subdir
    phs_config%n_channel = 2
    phs_config%n_par = 2
    phs_config%sqrts = sqrts
    if (.not. present (nlo_type)) &
      phs_config%nlo_type = BORN
    if (present (sqrts_fixed)) then
       phs_config%sqrts_fixed = sqrts_fixed
    end if
    if (present (lab_is_cm)) then
       phs_config%lab_is_cm = lab_is_cm
    end if
    if (present (azimuthal_dependence)) then
       phs_config%azimuthal_dependence = azimuthal_dependence
    end if
    if (allocated (phs_config%channel))  deallocate (phs_config%channel)
    allocate (phs_config%channel (phs_config%n_channel))
    if (phs_config%create_equivalences) then
       call setup_test_equivalences (phs_config)
       call setup_test_channel_props (phs_config)
    end if
    call phs_config%compute_md5sum ()
  end subroutine phs_test_config_configure

  subroutine setup_test_equivalences (phs_config)
    class(phs_test_config_t), intent(inout) :: phs_config
    integer :: i
    associate (channel => phs_config%channel(1))
      allocate (channel%eq (2))
      do i = 1, size (channel%eq)
         call channel%eq(i)%init (phs_config%n_par)
      end do
      associate (eq => channel%eq(1))
        eq%c = 1;  eq%perm = [1, 2];  eq%mode = [EQ_IDENTITY, EQ_SYMMETRIC]
      end associate
      associate (eq => channel%eq(2))
        eq%c = 2;  eq%perm = [2, 1];  eq%mode = [EQ_INVARIANT, EQ_IDENTITY]
      end associate
    end associate
  end subroutine setup_test_equivalences

  subroutine setup_test_channel_props (phs_config)
    class(phs_test_config_t), intent(inout) :: phs_config
    associate (channel => phs_config%channel(2))
      call channel%set_resonant (140._default, 3.1415_default)
    end associate
  end subroutine setup_test_channel_props

  subroutine phs_test_config_startup_message (phs_config, unit)
    class(phs_test_config_t), intent(in) :: phs_config
    integer, intent(in), optional :: unit
    call phs_config%base_startup_message (unit)
    write (msg_buffer, "(A)") "Phase space: Test"
    call msg_message (unit = unit)
  end subroutine phs_test_config_startup_message

  subroutine phs_test_config_allocate_instance (phs)
    class(phs_t), intent(inout), pointer :: phs
    allocate (phs_test_t :: phs)
  end subroutine phs_test_config_allocate_instance

  subroutine phs_test_write (object, unit, verbose)
    class(phs_test_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    logical :: verb
    u = given_output_unit (unit)
    verb = .false.;  if (present (verbose))  verb = verbose
    if (verb) then
       write (u, "(1x,A)")  "Partonic phase space: data"
       write (u, "(3x,A," // FMT_19 // ")")  "m = ", object%m
    end if
    call object%base_write (u)
  end subroutine phs_test_write

  subroutine phs_test_final (object)
    class(phs_test_t), intent(inout) :: object
  end subroutine phs_test_final

  subroutine phs_test_init (phs, phs_config)
    class(phs_test_t), intent(out) :: phs
    class(phs_config_t), intent(in), target :: phs_config
    call phs%base_init (phs_config)
    phs%m = phs%config%flv(1,1)%get_mass ()
    allocate (phs%x (phs_config%n_par), source = 0._default)
  end subroutine phs_test_init

  subroutine phs_test_evaluate_selected_channel (phs, c_in, r_in)
    class(phs_test_t), intent(inout) :: phs
    integer, intent(in) :: c_in
    real(default), intent(in), dimension(:) :: r_in
    if (phs%p_defined) then
       call phs%select_channel (c_in)
       phs%r(:,c_in) = r_in
       select case (c_in)
       case (1)
          phs%x = r_in
       case (2)
          phs%x(1) = r_in(1) ** (1 / 3._default)
          phs%x(2) = r_in(2)
       end select
       call compute_kinematics_solid_angle (phs%p, phs%q, phs%x)
       phs%volume = 1
       phs%q_defined = .true.
    end if
  end subroutine phs_test_evaluate_selected_channel

  subroutine phs_test_evaluate_other_channels (phs, c_in)
    class(phs_test_t), intent(inout) :: phs
    integer, intent(in) :: c_in
    integer :: c, n_channel
    if (phs%p_defined) then
       n_channel = phs%config%n_channel
       do c = 1, n_channel
          if (c /= c_in) then
             call inverse_kinematics_solid_angle (phs%p, phs%q, phs%x)
             select case (c)
             case (1)
                phs%r(:,c) = phs%x
             case (2)
                phs%r(1,c) = phs%x(1) ** 3
                phs%r(2,c) = phs%x(2)
             end select
          end if
       end do
       phs%f(1) = 1
       if (phs%r(1,2) /= 0) then
          phs%f(2) = 1 / (3 * phs%r(1,2) ** (2/3._default))
       else
          phs%f(2) = 0
       end if
       phs%r_defined = .true.
    end if
  end subroutine phs_test_evaluate_other_channels

  subroutine phs_test_inverse (phs)
    class(phs_test_t), intent(inout) :: phs
    integer :: c, n_channel
    real(default), dimension(:), allocatable :: x
    if (phs%p_defined .and. phs%q_defined) then
       call phs%select_channel ()
       n_channel = phs%config%n_channel
       allocate (x (phs%config%n_par))
       do c = 1, n_channel
          call inverse_kinematics_solid_angle (phs%p, phs%q, x)
          select case (c)
          case (1)
             phs%r(:,c) = x
          case (2)
             phs%r(1,c) = x(1) ** 3
             phs%r(2,c) = x(2)
          end select
       end do
       phs%f(1) = 1
       if (phs%r(1,2) /= 0) then
          phs%f(2) = 1 / (3 * phs%r(1,2) ** (2/3._default))
       else
          phs%f(2) = 0
       end if
       phs%volume = 1
       phs%r_defined = .true.
    end if
  end subroutine phs_test_inverse


end module phs_base_uti
