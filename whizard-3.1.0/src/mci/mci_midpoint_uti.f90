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

module mci_midpoint_uti

  use kinds, only: default
  use io_units
  use rng_base
  use mci_base

  use mci_midpoint

  use rng_base_ut, only: rng_test_t

  implicit none
  private

  public :: mci_midpoint_1
  public :: mci_midpoint_2
  public :: mci_midpoint_3
  public :: mci_midpoint_4
  public :: mci_midpoint_5
  public :: mci_midpoint_6
  public :: mci_midpoint_7

  type, extends (mci_sampler_t) :: test_sampler_1_t
     real(default), dimension(:), allocatable :: x
     real(default) :: val
   contains
     procedure :: write => test_sampler_1_write
     procedure :: evaluate => test_sampler_1_evaluate
     procedure :: is_valid => test_sampler_1_is_valid
     procedure :: rebuild => test_sampler_1_rebuild
     procedure :: fetch => test_sampler_1_fetch
  end type test_sampler_1_t

  type, extends (mci_sampler_t) :: test_sampler_2_t
     real(default) :: val
     real(default), dimension(2) :: x
   contains
     procedure :: write => test_sampler_2_write
     procedure :: evaluate => test_sampler_2_evaluate
     procedure :: is_valid => test_sampler_2_is_valid
     procedure :: rebuild => test_sampler_2_rebuild
     procedure :: fetch => test_sampler_2_fetch
  end type test_sampler_2_t

  type, extends (mci_sampler_t) :: test_sampler_4_t
     real(default) :: val
     real(default), dimension(:), allocatable :: x
   contains
     procedure :: write => test_sampler_4_write
     procedure :: evaluate => test_sampler_4_evaluate
     procedure :: is_valid => test_sampler_4_is_valid
     procedure :: rebuild => test_sampler_4_rebuild
     procedure :: fetch => test_sampler_4_fetch
  end type test_sampler_4_t


contains

  subroutine test_sampler_1_write (object, unit, testflag)
    class(test_sampler_1_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)") "Test sampler: f(x) = 3 x^2"
  end subroutine test_sampler_1_write

  subroutine test_sampler_1_evaluate (sampler, c, x_in, val, x, f)
    class(test_sampler_1_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    if (allocated (sampler%x))  deallocate (sampler%x)
    allocate (sampler%x (size (x_in)))
    sampler%x = x_in
    sampler%val = 3 * x_in(1) ** 2
    call sampler%fetch (val, x, f)
  end subroutine test_sampler_1_evaluate

  function test_sampler_1_is_valid (sampler) result (valid)
    class(test_sampler_1_t), intent(in) :: sampler
    logical :: valid
    valid = .true.
  end function test_sampler_1_is_valid

  subroutine test_sampler_1_rebuild (sampler, c, x_in, val, x, f)
    class(test_sampler_1_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(in) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    if (allocated (sampler%x))  deallocate (sampler%x)
    allocate (sampler%x (size (x_in)))
    sampler%x = x_in
    sampler%val = val
    x(:,1) = sampler%x
    f = 1
  end subroutine test_sampler_1_rebuild

  subroutine test_sampler_1_fetch (sampler, val, x, f)
    class(test_sampler_1_t), intent(in) :: sampler
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    val = sampler%val
    x(:,1) = sampler%x
    f = 1
  end subroutine test_sampler_1_fetch

  subroutine test_sampler_2_write (object, unit, testflag)
    class(test_sampler_2_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)") "Test sampler: f(x) = 3 x^2 + 2 y"
  end subroutine test_sampler_2_write

  subroutine test_sampler_2_evaluate (sampler, c, x_in, val, x, f)
    class(test_sampler_2_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    sampler%x = x_in
    sampler%val = 3 * x_in(1) ** 2 + 2 * x_in(2)
    call sampler%fetch (val, x, f)
  end subroutine test_sampler_2_evaluate

  function test_sampler_2_is_valid (sampler) result (valid)
    class(test_sampler_2_t), intent(in) :: sampler
    logical :: valid
    valid = .true.
  end function test_sampler_2_is_valid

  subroutine test_sampler_2_rebuild (sampler, c, x_in, val, x, f)
    class(test_sampler_2_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(in) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    sampler%x = x_in
    sampler%val = val
    x(:,1) = sampler%x
    f = 1
  end subroutine test_sampler_2_rebuild

  subroutine test_sampler_2_fetch (sampler, val, x, f)
    class(test_sampler_2_t), intent(in) :: sampler
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    val = sampler%val
    x(:,1) = sampler%x
    f = 1
  end subroutine test_sampler_2_fetch

  subroutine test_sampler_4_write (object, unit, testflag)
    class(test_sampler_4_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)") "Test sampler: f(x) = 1 - 3 x^2"
  end subroutine test_sampler_4_write

  subroutine test_sampler_4_evaluate (sampler, c, x_in, val, x, f)
    class(test_sampler_4_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    if (x_in(1) >= .5_default) then
       sampler%val = 1 - 3 * x_in(1) ** 2
    else
       sampler%val = 0
    end if
    if (.not. allocated (sampler%x))  allocate (sampler%x (size (x_in)))
    sampler%x = x_in
    call sampler%fetch (val, x, f)
  end subroutine test_sampler_4_evaluate

  function test_sampler_4_is_valid (sampler) result (valid)
    class(test_sampler_4_t), intent(in) :: sampler
    logical :: valid
    valid = .true.
  end function test_sampler_4_is_valid

  subroutine test_sampler_4_rebuild (sampler, c, x_in, val, x, f)
    class(test_sampler_4_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(in) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    sampler%x = x_in
    sampler%val = val
    x(:,1) = sampler%x
    f = 1
  end subroutine test_sampler_4_rebuild

  subroutine test_sampler_4_fetch (sampler, val, x, f)
    class(test_sampler_4_t), intent(in) :: sampler
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    val = sampler%val
    x(:,1) = sampler%x
    f = 1
  end subroutine test_sampler_4_fetch

  subroutine mci_midpoint_1 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler

    write (u, "(A)")  "* Test output: mci_midpoint_1"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_midpoint_t :: mci)
    call mci%set_dimensions (1, 1)

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_1_t :: sampler)
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate for x = 0.8"
    write (u, "(A)")

    call mci_instance%evaluate (sampler, 1, [0.8_default])
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate for x = 0.7"
    write (u, "(A)")

    call mci_instance%evaluate (sampler, 1, [0.7_default])
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate for x = 0.9"
    write (u, "(A)")

    call mci_instance%evaluate (sampler, 1, [0.9_default])
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_calls = 1000"
    write (u, "(A)")

    call mci%integrate (mci_instance, sampler, 1, 1000)
    call mci%write (u)

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_midpoint_1"

  end subroutine mci_midpoint_1

  subroutine mci_midpoint_2 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler

    write (u, "(A)")  "* Test output: mci_midpoint_2"
    write (u, "(A)")  "*   Purpose: integrate function in two dimensions"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_midpoint_t :: mci)
    call mci%set_dimensions (2, 1)

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_2_t :: sampler)
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate for x = 0.8, y = 0.2"
    write (u, "(A)")

    call mci_instance%evaluate (sampler, 1, [0.8_default, 0.2_default])
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_calls = 1000"
    write (u, "(A)")

    call mci%integrate (mci_instance, sampler, 1, 1000)
    call mci%write (u)

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_midpoint_2"

  end subroutine mci_midpoint_2

  subroutine mci_midpoint_3 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler

    write (u, "(A)")  "* Test output: mci_midpoint_3"
    write (u, "(A)")  "*   Purpose: integrate function with one flat dimension"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_midpoint_t :: mci)
    select type (mci)
    type is (mci_midpoint_t)
       call mci%set_dimensions (2, 1)
       call mci%declare_flat_dimensions ([2])
    end select

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_1_t :: sampler)
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate for x = 0.8, y = 0.2"
    write (u, "(A)")

    call mci_instance%evaluate (sampler, 1, [0.8_default, 0.2_default])
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_calls = 1000"
    write (u, "(A)")

    call mci%integrate (mci_instance, sampler, 1, 1000)
    call mci%write (u)

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_midpoint_3"

  end subroutine mci_midpoint_3

  subroutine mci_midpoint_4 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler

    write (u, "(A)")  "* Test output: mci_midpoint_4"
    write (u, "(A)")  "*   Purpose: integrate function with sign flip &
         &in one dimension"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_midpoint_t :: mci)
    call mci%set_dimensions (1, 1)

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_4_t :: sampler)
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate for x = 0.8"
    write (u, "(A)")

    call mci_instance%evaluate (sampler, 1, [0.8_default])
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_calls = 1000"
    write (u, "(A)")

    call mci%integrate (mci_instance, sampler, 1, 1000)
    call mci%write (u)

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_midpoint_4"

  end subroutine mci_midpoint_4

  subroutine mci_midpoint_5 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng
    class(mci_state_t), allocatable :: state

    write (u, "(A)")  "* Test output: mci_midpoint_5"
    write (u, "(A)")  "*   Purpose: generate weighted events"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_midpoint_t :: mci)
    call mci%set_dimensions (2, 1)

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_2_t :: sampler)

    write (u, "(A)")  "* Initialize random-number generator"
    write (u, "(A)")

    allocate (rng_test_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    write (u, "(A)")  "* Generate weighted event"
    write (u, "(A)")

    call mci%generate_weighted_event (mci_instance, sampler)
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate weighted event"
    write (u, "(A)")

    call mci%generate_weighted_event (mci_instance, sampler)
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

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    deallocate (mci_instance)
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_midpoint_5"

  end subroutine mci_midpoint_5

  subroutine mci_midpoint_6 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_midpoint_6"
    write (u, "(A)")  "*   Purpose: generate unweighted events"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_midpoint_t :: mci)
    call mci%set_dimensions (1, 1)

    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_4_t :: sampler)

    write (u, "(A)")  "* Initialize random-number generator"
    write (u, "(A)")

    allocate (rng_test_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    write (u, "(A)")  "* Integrate (determine maximum of integrand"
    write (u, "(A)")
    call mci%integrate (mci_instance, sampler, 1, 1000)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate unweighted event"
    write (u, "(A)")

    call mci%generate_unweighted_event (mci_instance, sampler)
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    deallocate (mci_instance)
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_midpoint_6"

  end subroutine mci_midpoint_6

  subroutine mci_midpoint_7 (u)
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_midpoint_7"
    write (u, "(A)")  "*   Purpose: generate unweighted event &
         &with excess weight"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_midpoint_t :: mci)
    call mci%set_dimensions (1, 1)

    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_4_t :: sampler)

    write (u, "(A)")  "* Initialize random-number generator"
    write (u, "(A)")

    allocate (rng_test_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    write (u, "(A)")  "* Integrate (determine maximum of integrand"
    write (u, "(A)")
    call mci%integrate (mci_instance, sampler, 1, 2)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate unweighted event"
    write (u, "(A)")

    call mci_instance%init_simulation ()
    call mci%generate_unweighted_event (mci_instance, sampler)
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Use getter methods"
    write (u, "(A)")

    write (u, "(1x,A,1x,ES19.12)")  "weight =", mci_instance%get_event_weight ()
    write (u, "(1x,A,1x,ES19.12)")  "excess =", mci_instance%get_event_excess ()

    write (u, "(A)")
    write (u, "(A)")  "* Apply safety factor"
    write (u, "(A)")

    call mci_instance%init_simulation (safety_factor = 2.1_default)

    write (u, "(A)")  "* Generate unweighted event"
    write (u, "(A)")

    call mci%generate_unweighted_event (mci_instance, sampler)
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Use getter methods"
    write (u, "(A)")

    write (u, "(1x,A,1x,ES19.12)")  "weight =", mci_instance%get_event_weight ()
    write (u, "(1x,A,1x,ES19.12)")  "excess =", mci_instance%get_event_excess ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    deallocate (mci_instance)
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_midpoint_7"

  end subroutine mci_midpoint_7


end module mci_midpoint_uti
