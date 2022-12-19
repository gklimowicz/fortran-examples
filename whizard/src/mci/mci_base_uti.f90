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

module mci_base_uti

  use kinds, only: default
  use io_units
  use diagnostics
  use phs_base
  use rng_base

  use mci_base

  use rng_base_ut, only: rng_test_t

  implicit none
  private

  public :: mci_test_t

  public :: mci_base_1
  public :: mci_base_2
  public :: mci_base_3
  public :: mci_base_4
  public :: mci_base_5
  public :: mci_base_6
  public :: mci_base_7
  public :: mci_base_8

  type, extends (mci_t) :: mci_test_t
     integer :: divisions = 0
     integer :: tries = 0
     real(default) :: max_factor = 1
   contains
     procedure :: final => mci_test_final
     procedure :: write => mci_test_write
     procedure :: startup_message => mci_test_startup_message
     procedure :: write_log_entry => mci_test_write_log_entry
     procedure :: compute_md5sum => mci_test_compute_md5sum
     procedure :: declare_flat_dimensions => mci_test_ignore_flat_dimensions
     procedure :: declare_equivalences => mci_test_ignore_equivalences
     procedure :: set_divisions => mci_test_set_divisions
     procedure :: set_max_factor => mci_test_set_max_factor
     procedure :: allocate_instance => mci_test_allocate_instance
     procedure :: integrate => mci_test_integrate
     procedure :: prepare_simulation => mci_test_ignore_prepare_simulation
     procedure :: generate_weighted_event => mci_test_generate_weighted_event
     procedure :: generate_unweighted_event => &
          mci_test_generate_unweighted_event
     procedure :: rebuild_event => mci_test_rebuild_event
  end type mci_test_t

  type, extends (mci_instance_t) :: mci_test_instance_t
     type(mci_test_t), pointer :: mci => null ()
     real(default) :: g = 0
     real(default), dimension(:), allocatable :: gi
     real(default) :: value = 0
     real(default) :: rel_value = 0
     real(default), dimension(:), allocatable :: max
   contains
     procedure :: write => mci_test_instance_write
     procedure :: final => mci_test_instance_final
     procedure :: init => mci_test_instance_init
     procedure :: compute_weight => mci_test_instance_compute_weight
     procedure :: record_integrand => mci_test_instance_record_integrand
     procedure :: init_simulation => mci_test_instance_init_simulation
     procedure :: final_simulation => mci_test_instance_final_simulation
     procedure :: get_event_excess => mci_test_instance_get_event_excess
  end type mci_test_instance_t

  type, extends (mci_sampler_t) :: test_sampler_t
     real(default) :: integrand = 0
     integer :: selected_channel = 0
     real(default), dimension(:,:), allocatable :: x
     real(default), dimension(:), allocatable :: f
   contains
     procedure :: init => test_sampler_init
     procedure :: write => test_sampler_write
     procedure :: compute => test_sampler_compute
     procedure :: is_valid => test_sampler_is_valid
     procedure :: evaluate => test_sampler_evaluate
     procedure :: rebuild => test_sampler_rebuild
     procedure :: fetch => test_sampler_fetch
  end type test_sampler_t

  type, extends (mci_results_t) :: mci_test_results_t
     integer :: n_it = 0
     integer :: n_calls = 0
     real(default) :: integral = 0
     real(default) :: error = 0
     real(default) :: efficiency = 0
   contains
     procedure :: write => mci_test_results_write
     procedure :: write_verbose => mci_test_results_write_verbose
     procedure :: record_simple => mci_test_results_record_simple
     procedure :: record_extended => mci_test_results_record_extended
  end type mci_test_results_t


contains

  subroutine mci_test_final (object)
    class(mci_test_t), intent(inout) :: object
    call object%base_final ()
  end subroutine mci_test_final

  subroutine mci_test_write (object, unit, pacify, md5sum_version)
    class(mci_test_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: pacify
    logical, intent(in), optional :: md5sum_version
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Test integrator:"
    call object%base_write (u, pacify, md5sum_version)
    if (object%divisions /= 0) then
       write (u, "(3x,A,I0)")  "Number of divisions  = ", object%divisions
    end if
    if (allocated (object%rng))  call object%rng%write (u)
  end subroutine mci_test_write

  subroutine mci_test_startup_message (mci, unit, n_calls)
    class(mci_test_t), intent(in) :: mci
    integer, intent(in), optional :: unit, n_calls
    call mci%base_startup_message (unit = unit, n_calls = n_calls)
    write (msg_buffer, "(A,1x,I0,1x,A)") &
         "Integrator: Test:", mci%divisions, "divisions"
    call msg_message (unit = unit)
  end subroutine mci_test_startup_message

  subroutine mci_test_write_log_entry (mci, u)
    class(mci_test_t), intent(in) :: mci
    integer, intent(in) :: u
  end subroutine mci_test_write_log_entry

  subroutine mci_test_compute_md5sum (mci, pacify)
    class(mci_test_t), intent(inout) :: mci
    logical, intent(in), optional :: pacify
  end subroutine mci_test_compute_md5sum

  subroutine mci_test_ignore_flat_dimensions (mci, dim_flat)
    class(mci_test_t), intent(inout) :: mci
    integer, dimension(:), intent(in) :: dim_flat
  end subroutine mci_test_ignore_flat_dimensions

  subroutine mci_test_ignore_equivalences (mci, channel, dim_offset)
    class(mci_test_t), intent(inout) :: mci
    type(phs_channel_t), dimension(:), intent(in) :: channel
    integer, intent(in) :: dim_offset
  end subroutine mci_test_ignore_equivalences

  subroutine mci_test_set_divisions (object, divisions)
    class(mci_test_t), intent(inout) :: object
    integer, intent(in) :: divisions
    object%divisions = divisions
  end subroutine mci_test_set_divisions

  subroutine mci_test_set_max_factor (object, max_factor)
    class(mci_test_t), intent(inout) :: object
    real(default), intent(in) :: max_factor
    object%max_factor = max_factor
  end subroutine mci_test_set_max_factor

  subroutine mci_test_allocate_instance (mci, mci_instance)
    class(mci_test_t), intent(in) :: mci
    class(mci_instance_t), intent(out), pointer :: mci_instance
    allocate (mci_test_instance_t :: mci_instance)
  end subroutine mci_test_allocate_instance

  subroutine mci_test_integrate (mci, instance, sampler, &
       n_it, n_calls, results, pacify)
    class(mci_test_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout), target :: instance
    class(mci_sampler_t), intent(inout), target :: sampler
    integer, intent(in) :: n_it
    integer, intent(in) :: n_calls
    logical, intent(in), optional :: pacify
    class(mci_results_t), intent(inout), optional :: results
    real(default), dimension(:), allocatable :: integral
    real(default), dimension(:), allocatable :: x
    integer :: i, j, c
    select type (instance)
    type is (mci_test_instance_t)
       allocate (integral (mci%n_channel))
       integral = 0
       allocate (x (mci%n_dim))
       select case (mci%n_dim)
       case (1)
          do c = 1, mci%n_channel
             do i = 1, mci%divisions
                x(1) = (i - 0.5_default) / mci%divisions
                call instance%evaluate (sampler, c, x)
                integral(c) = integral(c) + instance%get_value ()
             end do
          end do
          mci%integral = dot_product (instance%w, integral) &
               / mci%divisions
          mci%integral_known = .true.
       case (2)
          do c = 1, mci%n_channel
             do i = 1, mci%divisions
                x(1) = (i - 0.5_default) / mci%divisions
                do j = 1, mci%divisions
                   x(2) = (j - 0.5_default) / mci%divisions
                   call instance%evaluate (sampler, c, x)
                   integral(c) = integral(c) + instance%get_value ()
                end do
             end do
          end do
          mci%integral = dot_product (instance%w, integral) &
               / mci%divisions / mci%divisions
          mci%integral_known = .true.
       end select
       if (present (results)) then
          call results%record (n_it, n_calls, &
               mci%integral, mci%error, &
               efficiency = 0._default)
       end if
    end select
  end subroutine mci_test_integrate

  subroutine mci_test_ignore_prepare_simulation (mci)
    class(mci_test_t), intent(inout) :: mci
  end subroutine mci_test_ignore_prepare_simulation

  subroutine mci_test_generate_weighted_event (mci, instance, sampler)
    class(mci_test_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout), target :: instance
    class(mci_sampler_t), intent(inout), target :: sampler
    real(default) :: r
    real(default), dimension(:), allocatable :: x
    integer :: c
    select type (instance)
    type is (mci_test_instance_t)
       allocate (x (mci%n_dim))
       select case (mci%n_channel)
       case (1)
          c = 1
          call mci%rng%generate (x(1))
       case (2)
          call mci%rng%generate (r)
          if (r < instance%w(1)) then
             c = 1
          else
             c = 2
          end if
          call mci%rng%generate (x)
       end select
       call instance%evaluate (sampler, c, x)
    end select
  end subroutine mci_test_generate_weighted_event

  subroutine mci_test_generate_unweighted_event (mci, instance, sampler)
    class(mci_test_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout), target :: instance
    class(mci_sampler_t), intent(inout), target :: sampler
    real(default) :: r
    integer :: i
    select type (instance)
    type is (mci_test_instance_t)
       mci%tries = 0
       do i = 1, 10
          call mci%generate_weighted_event (instance, sampler)
          mci%tries = mci%tries + 1
          call mci%rng%generate (r)
          if (r < instance%rel_value)  exit
       end do
    end select
  end subroutine mci_test_generate_unweighted_event

  subroutine mci_test_rebuild_event (mci, instance, sampler, state)
    class(mci_test_t), intent(inout) :: mci
    class(mci_instance_t), intent(inout) :: instance
    class(mci_sampler_t), intent(inout) :: sampler
    class(mci_state_t), intent(in) :: state
    select type (instance)
    type is (mci_test_instance_t)
       call instance%recall (sampler, state)
    end select
  end subroutine mci_test_rebuild_event

  subroutine mci_test_instance_write (object, unit, pacify)
    class(mci_test_instance_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: pacify
    integer :: u, c
    u = given_output_unit (unit)
    write (u, "(1x,A,ES13.7)") "Result value = ", object%value
    write (u, "(1x,A,ES13.7)") "Rel. weight  = ", object%rel_value
    write (u, "(1x,A,ES13.7)") "Integrand    = ", object%integrand
    write (u, "(1x,A,ES13.7)") "MCI weight   = ", object%mci_weight
    write (u, "(3x,A,I0)")  "c = ", object%selected_channel
    write (u, "(3x,A,ES13.7)") "g = ", object%g
    write (u, "(1x,A)")  "Channel parameters:"
    do c = 1, object%mci%n_channel
       write (u, "(1x,I0,A,4(1x,ES13.7))")  c, ": w/f/g/m =", &
            object%w(c), object%f(c), object%gi(c), object%max(c)
       write (u, "(4x,A,9(1x,F9.7))")  "x =", object%x(:,c)
    end do
  end subroutine mci_test_instance_write

  subroutine mci_test_instance_final (object)
    class(mci_test_instance_t), intent(inout) :: object
  end subroutine mci_test_instance_final

  subroutine mci_test_instance_init (mci_instance, mci)
    class(mci_test_instance_t), intent(out) :: mci_instance
    class(mci_t), intent(in), target :: mci
    call mci_instance%base_init (mci)
    select type (mci)
    type is (mci_test_t)
       mci_instance%mci => mci
    end select
    allocate (mci_instance%gi (mci%n_channel))
    mci_instance%gi = 0
    allocate (mci_instance%max (mci%n_channel))
    select case (mci%n_channel)
    case (1)
       mci_instance%max = 1._default
    case (2)
       mci_instance%max = 2._default
    end select
  end subroutine mci_test_instance_init

  subroutine mci_test_instance_compute_weight (mci, c)
    class(mci_test_instance_t), intent(inout) :: mci
    integer, intent(in) :: c
    integer :: i
    mci%selected_channel = c
    select case (mci%mci%n_dim)
    case (1)
       mci%gi(1) = 1
    case (2)
       mci%gi(1) = 1
       mci%gi(2) = 2 * mci%x(2,2)
    end select
    mci%g = 0
    do i = 1, mci%mci%n_channel
       mci%g = mci%g + mci%w(i) * mci%gi(i) / mci%f(i)
    end do
    mci%mci_weight = mci%gi(c) / mci%g
  end subroutine mci_test_instance_compute_weight

  subroutine mci_test_instance_record_integrand (mci, integrand)
    class(mci_test_instance_t), intent(inout) :: mci
    real(default), intent(in) :: integrand
    mci%integrand = integrand
    mci%value = mci%integrand * mci%mci_weight
    mci%rel_value = mci%value / mci%max(mci%selected_channel) &
         / mci%mci%max_factor
  end subroutine mci_test_instance_record_integrand

  subroutine mci_test_instance_init_simulation (instance, safety_factor)
    class(mci_test_instance_t), intent(inout) :: instance
    real(default), intent(in), optional :: safety_factor
  end subroutine mci_test_instance_init_simulation

  subroutine mci_test_instance_final_simulation (instance)
    class(mci_test_instance_t), intent(inout) :: instance
  end subroutine mci_test_instance_final_simulation

  function mci_test_instance_get_event_excess (mci) result (excess)
    class(mci_test_instance_t), intent(in) :: mci
    real(default) :: excess
    excess = 0
  end function mci_test_instance_get_event_excess

  subroutine test_sampler_init (sampler, n)
    class(test_sampler_t), intent(out) :: sampler
    integer, intent(in) :: n
    allocate (sampler%x (n, n))
    allocate (sampler%f (n))
  end subroutine test_sampler_init

  subroutine test_sampler_write (object, unit, testflag)
    class(test_sampler_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u, c
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Test sampler:"
    write (u, "(3x,A,ES13.7)")  "Integrand = ", object%integrand
    write (u, "(3x,A,I0)")      "Channel   = ", object%selected_channel
    do c = 1, size (object%f)
       write (u, "(1x,I0,':',1x,A,ES13.7)") c, "f = ", object%f(c)
       write (u, "(4x,A,9(1x,F9.7))") "x =", object%x(:,c)
    end do
  end subroutine test_sampler_write

  subroutine test_sampler_compute (sampler, c, x_in)
    class(test_sampler_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    sampler%selected_channel = c
    select case (size (sampler%f))
    case (1)
       sampler%x(:,1) = x_in
       sampler%f = 1
    case (2)
       select case (c)
       case (1)
          sampler%x(:,1) = x_in
          sampler%x(1,2) = sqrt (x_in(1))
          sampler%x(2,2) = x_in(2)
       case (2)
          sampler%x(1,1) = x_in(1) ** 2
          sampler%x(2,1) = x_in(2)
          sampler%x(:,2) = x_in
       end select
       sampler%f(1) = 1
       sampler%f(2) = 2 * sampler%x(1,2)
    end select
  end subroutine test_sampler_compute

  function test_sampler_is_valid (sampler) result (valid)
    class(test_sampler_t), intent(in) :: sampler
    logical :: valid
    valid = .true.
  end function test_sampler_is_valid

  subroutine test_sampler_evaluate (sampler, c, x_in, val, x, f)
    class(test_sampler_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    call sampler%compute (c, x_in)
    sampler%integrand = 1
    val = sampler%integrand
    x = sampler%x
    f = sampler%f
  end subroutine test_sampler_evaluate

  subroutine test_sampler_rebuild (sampler, c, x_in, val, x, f)
    class(test_sampler_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(in) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    call sampler%compute (c, x_in)
    sampler%integrand = val
    x = sampler%x
    f = sampler%f
  end subroutine test_sampler_rebuild

  subroutine test_sampler_fetch (sampler, val, x, f)
    class(test_sampler_t), intent(in) :: sampler
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    val = sampler%integrand
    x = sampler%x
    f = sampler%f
  end subroutine test_sampler_fetch

  subroutine mci_test_results_write (object, unit, suppress)
    class(mci_test_results_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: suppress
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A,1x,I0)") "Iterations = ", object%n_it
    write (u, "(3x,A,1x,I0)") "Calls      = ", object%n_calls
    write (u, "(3x,A,1x,F12.10)")  "Integral   = ", object%integral
    write (u, "(3x,A,1x,F12.10)")  "Error      = ", object%error
    write (u, "(3x,A,1x,F12.10)")  "Efficiency = ", object%efficiency
  end subroutine mci_test_results_write

  subroutine mci_test_results_write_verbose (object, unit)
    class(mci_test_results_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A,1x,I0)") "Iterations = ", object%n_it
    write (u, "(3x,A,1x,I0)") "Calls      = ", object%n_calls
    write (u, "(3x,A,1x,F12.10)")  "Integral   = ", object%integral
    write (u, "(3x,A,1x,F12.10)")  "Error      = ", object%error
    write (u, "(3x,A,1x,F12.10)")  "Efficiency = ", object%efficiency
  end subroutine mci_test_results_write_verbose

  subroutine mci_test_results_record_simple (object, n_it, n_calls, &
       integral, error, efficiency, chain_weights, suppress)
    class(mci_test_results_t), intent(inout) :: object
    integer, intent(in) :: n_it
    integer, intent(in) :: n_calls
    real(default), intent(in) :: integral
    real(default), intent(in) :: error
    real(default), intent(in) :: efficiency
    real(default), dimension(:), intent(in), optional :: chain_weights
    logical, intent(in), optional :: suppress
    object%n_it = n_it
    object%n_calls = n_calls
    object%integral = integral
    object%error = error
    object%efficiency = efficiency
  end subroutine mci_test_results_record_simple

  subroutine mci_test_results_record_extended (object, n_it, n_calls, &
       & n_calls_valid, integral, error, efficiency, efficiency_pos, &
       & efficiency_neg, chain_weights, suppress)
    class(mci_test_results_t), intent(inout) :: object
    integer, intent(in) :: n_it
    integer, intent(in) :: n_calls
    integer, intent(in) :: n_calls_valid
    real(default), intent(in) :: integral
    real(default), intent(in) :: error
    real(default), intent(in) :: efficiency
    real(default), intent(in) :: efficiency_pos
    real(default), intent(in) :: efficiency_neg
    real(default), dimension(:), intent(in), optional :: chain_weights
    logical, intent(in), optional :: suppress
    object%n_it = n_it
    object%n_calls = n_calls
    object%integral = integral
    object%error = error
    object%efficiency = efficiency
  end subroutine mci_test_results_record_extended

  subroutine mci_base_1 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler

    real(default) :: integrand

    write (u, "(A)")  "* Test output: mci_base_1"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &test integrator"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_test_t :: mci)
    call mci%set_dimensions (2, 2)

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_t :: sampler)
    select type (sampler)
    type is (test_sampler_t)
       call sampler%init (2)
    end select

    write (u, "(A)")  "* Evaluate sampler for given point and channel"
    write (u, "(A)")

    call sampler%evaluate (1, [0.25_default, 0.8_default], &
         integrand, mci_instance%x, mci_instance%f)

    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute MCI weight"
    write (u, "(A)")

    call mci_instance%compute_weight (1)
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Get integrand and compute weight for another point"
    write (u, "(A)")

    call mci_instance%evaluate (sampler, 2, [0.5_default, 0.6_default])
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recall results, again"
    write (u, "(A)")

    call mci_instance%final ()
    deallocate (mci_instance)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    call mci_instance%fetch (sampler, 2)
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Retrieve value"
    write (u, "(A)")

    write (u, "(1x,A,ES13.7)")  "Weighted integrand = ", &
         mci_instance%get_value ()

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_base_1"

  end subroutine mci_base_1

  subroutine mci_base_2 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler

    write (u, "(A)")  "* Test output: mci_base_2"
    write (u, "(A)")  "*   Purpose: perform a test integral"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_test_t :: mci)
    call mci%set_dimensions (1, 1)
    select type (mci)
    type is (mci_test_t)
       call mci%set_divisions (10)
    end select

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_t :: sampler)
    select type (sampler)
    type is (test_sampler_t)
       call sampler%init (1)
    end select

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call mci%integrate (mci_instance, sampler, 0, 0)

    call mci%write (u)

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_base_2"

  end subroutine mci_base_2

  subroutine mci_base_3 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler

    write (u, "(A)")  "* Test output: mci_base_3"
    write (u, "(A)")  "*   Purpose: perform a nontrivial test integral"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_test_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_test_t)
       call mci%set_divisions (10)
    end select

    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_t :: sampler)
    select type (sampler)
    type is (test_sampler_t)
       call sampler%init (2)
    end select

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call mci%integrate (mci_instance, sampler, 0, 0)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with higher resolution"
    write (u, "(A)")

    select type (mci)
    type is (mci_test_t)
       call mci%set_divisions (100)
    end select

    call mci%integrate (mci_instance, sampler, 0, 0)
    call mci%write (u)

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_base_3"

  end subroutine mci_base_3

  subroutine mci_base_4 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_base_4"
    write (u, "(A)")  "*   Purpose: generate events"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator, instance, sampler"
    write (u, "(A)")

    allocate (mci_test_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_test_t)
       call mci%set_divisions (10)
    end select

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_t :: sampler)
    select type (sampler)
    type is (test_sampler_t)
       call sampler%init (2)
    end select

    allocate (rng_test_t :: rng)
    call mci%import_rng (rng)

    write (u, "(A)")  "* Generate weighted event"
    write (u, "(A)")

    call mci%generate_weighted_event (mci_instance, sampler)

    call sampler%write (u)
    write (u, *)
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate unweighted event"
    write (u, "(A)")

    call mci%generate_unweighted_event (mci_instance, sampler)

    select type (mci)
    type is (mci_test_t)
       write (u, "(A,I0)")  " Success in try ", mci%tries
       write (u, "(A)")
    end select

    call sampler%write (u)
    write (u, *)
    call mci_instance%write (u)

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_base_4"

  end subroutine mci_base_4

  subroutine mci_base_5 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng
    class(mci_state_t), allocatable :: state

    write (u, "(A)")  "* Test output: mci_base_5"
    write (u, "(A)")  "*   Purpose: store and recall an event"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator, instance, sampler"
    write (u, "(A)")

    allocate (mci_test_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_test_t)
       call mci%set_divisions (10)
    end select

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_t :: sampler)
    select type (sampler)
    type is (test_sampler_t)
       call sampler%init (2)
    end select

    allocate (rng_test_t :: rng)
    call mci%import_rng (rng)

    write (u, "(A)")  "* Generate weighted event"
    write (u, "(A)")

    call mci%generate_weighted_event (mci_instance, sampler)

    call sampler%write (u)
    write (u, *)
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Store data"
    write (u, "(A)")

    allocate (state)
    call mci_instance%store (state)
    call mci_instance%final ()
    deallocate (mci_instance)

    call state%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recall data and rebuild event"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)
    call mci%rebuild_event (mci_instance, sampler, state)

    call sampler%write (u)
    write (u, *)
    call mci_instance%write (u)

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_base_5"

  end subroutine mci_base_5

  subroutine mci_base_6 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci

    write (u, "(A)")  "* Test output: mci_base_6"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &test integrator with chains"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_test_t :: mci)
    call mci%set_dimensions (1, 5)

    write (u, "(A)")  "* Introduce chains"
    write (u, "(A)")

    call mci%declare_chains ([1, 2, 2, 1, 2])

    call mci%write (u)

    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_base_6"

  end subroutine mci_base_6

  subroutine mci_base_7 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(mci_results_t), allocatable :: results

    write (u, "(A)")  "* Test output: mci_base_7"
    write (u, "(A)")  "*   Purpose: perform a nontrivial test integral &
         &and record results"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_test_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_test_t)
       call mci%set_divisions (10)
    end select

    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_t :: sampler)
    select type (sampler)
    type is (test_sampler_t)
       call sampler%init (2)
    end select

    allocate (mci_test_results_t :: results)

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call mci%integrate (mci_instance, sampler, 1, 1000, results)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Display results"
    write (u, "(A)")

    call results%write (u)

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_base_7"

  end subroutine mci_base_7

  subroutine mci_base_8 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci

    real(default) :: dummy

    write (u, "(A)")  "* Test output: mci_base_8"
    write (u, "(A)")  "*   Purpose: check timer availability"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator with timer"
    write (u, "(A)")

    allocate (mci_test_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_test_t)
       call mci%set_divisions (10)
    end select

    call mci%set_timer (active = .true.)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Start timer"
    write (u, "(A)")

    call mci%start_timer ()
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Stop timer"
    write (u, "(A)")

    call mci%stop_timer ()
    write (u, "(A)")  " (ok)"

    write (u, "(A)")
    write (u, "(A)")  "* Readout"
    write (u, "(A)")

    dummy = mci%get_time ()
    write (u, "(A)")  " (ok)"

    write (u, "(A)")
    write (u, "(A)")  "* Deactivate timer"
    write (u, "(A)")

    call mci%set_timer (active = .false.)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_base_8"

  end subroutine mci_base_8


end module mci_base_uti
