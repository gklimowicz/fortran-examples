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

module mci_vamp_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use constants, only: PI, TWOPI
  use rng_base
  use rng_tao
  use phs_base
  use mci_base
  use vamp, only: vamp_write_grids !NODEP!

  use mci_vamp

  implicit none
  private

  public :: mci_vamp_1
  public :: mci_vamp_2
  public :: mci_vamp_3
  public :: mci_vamp_4
  public :: mci_vamp_5
  public :: mci_vamp_6
  public :: mci_vamp_7
  public :: mci_vamp_8
  public :: mci_vamp_9
  public :: mci_vamp_10
  public :: mci_vamp_11
  public :: mci_vamp_12
  public :: mci_vamp_13
  public :: mci_vamp_14
  public :: mci_vamp_15
  public :: mci_vamp_16

  type, extends (mci_sampler_t) :: test_sampler_1_t
     real(default), dimension(:), allocatable :: x
     real(default) :: val
     integer :: mode = 1
   contains
     procedure :: write => test_sampler_1_write
     procedure :: evaluate => test_sampler_1_evaluate
     procedure :: is_valid => test_sampler_1_is_valid
     procedure :: rebuild => test_sampler_1_rebuild
     procedure :: fetch => test_sampler_1_fetch
  end type test_sampler_1_t

  type, extends (mci_sampler_t) :: test_sampler_2_t
     real(default), dimension(:,:), allocatable :: x
     real(default), dimension(:), allocatable :: f
     real(default) :: val
   contains
     procedure :: write => test_sampler_2_write
     procedure :: compute => test_sampler_2_compute
     procedure :: evaluate => test_sampler_2_evaluate
     procedure :: is_valid => test_sampler_2_is_valid
     procedure :: rebuild => test_sampler_2_rebuild
     procedure :: fetch => test_sampler_2_fetch
  end type test_sampler_2_t

  type, extends (mci_sampler_t) :: test_sampler_3_t
     real(default), dimension(:,:), allocatable :: x
     real(default), dimension(:), allocatable :: f
     real(default) :: val
     real(default) :: a = 1
     real(default) :: b = 1
   contains
     procedure :: write => test_sampler_3_write
     procedure :: compute => test_sampler_3_compute
     procedure :: evaluate => test_sampler_3_evaluate
     procedure :: is_valid => test_sampler_3_is_valid
     procedure :: rebuild => test_sampler_3_rebuild
     procedure :: fetch => test_sampler_3_fetch
  end type test_sampler_3_t


contains

  subroutine test_sampler_1_write (object, unit, testflag)
    class(test_sampler_1_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    select case (object%mode)
    case (1)
       write (u, "(1x,A)") "Test sampler: f(x) = 3 x^2"
    case (2)
       write (u, "(1x,A)") "Test sampler: f(x) = 11 x^10"
    case (3)
       write (u, "(1x,A)") "Test sampler: f(x) = 11 x^10 * 2 * cos^2 (2 pi y)"
    case (4)
       write (u, "(1x,A)") "Test sampler: f(x) = (1 - 3 x^2) theta(x - 1/2)"
    end select
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
    select case (sampler%mode)
    case (1)
       sampler%val = 3 * x_in(1) ** 2
    case (2)
       sampler%val = 11 * x_in(1) ** 10
    case (3)
       sampler%val = 11 * x_in(1) ** 10 * 2 * cos (twopi * x_in(2)) ** 2
    case (4)
       if (x_in(1) >= .5_default) then
          sampler%val = 1 - 3 * x_in(1) ** 2
       else
          sampler%val = 0
       end if
    end select
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
    write (u, "(1x,A)") "Two-channel test sampler 2"
  end subroutine test_sampler_2_write

  subroutine test_sampler_2_compute (sampler, c, x_in)
    class(test_sampler_2_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default) :: xx, yy, uu, vv
    if (.not. allocated (sampler%x)) &
         allocate (sampler%x (size (x_in), 2))
    if (.not. allocated (sampler%f)) &
         allocate (sampler%f (2))
    select case (c)
    case (1)
       xx = x_in(1)
       yy = x_in(2)
       uu = xx * yy
       vv = (1 + log (xx/yy) / log (xx*yy)) / 2
    case (2)
       uu = x_in(1)
       vv = x_in(2)
       xx = uu ** vv
       yy = uu ** (1 - vv)
    end select
    sampler%val = (2 * sin (pi * xx) * sin (pi * yy)) ** 2 &
         + 2 * sin (pi * vv) ** 2
    sampler%f(1) = 1
    sampler%f(2) = abs (log (uu))
    sampler%x(:,1) = [xx, yy]
    sampler%x(:,2) = [uu, vv]
  end subroutine test_sampler_2_compute

  subroutine test_sampler_2_evaluate (sampler, c, x_in, val, x, f)
    class(test_sampler_2_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    call sampler%compute (c, x_in)
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
    call sampler%compute (c, x_in)
    x = sampler%x
    f = sampler%f
  end subroutine test_sampler_2_rebuild

  subroutine test_sampler_2_fetch (sampler, val, x, f)
    class(test_sampler_2_t), intent(in) :: sampler
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    val = sampler%val
    x = sampler%x
    f = sampler%f
  end subroutine test_sampler_2_fetch

  subroutine test_sampler_3_write (object, unit, testflag)
    class(test_sampler_3_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)") "Two-channel test sampler 3"
    write (u, "(3x,A,F5.2)")  "a = ", object%a
    write (u, "(3x,A,F5.2)")  "b = ", object%b
  end subroutine test_sampler_3_write

  subroutine test_sampler_3_compute (sampler, c, x_in)
    class(test_sampler_3_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default) :: u, v, xx
    if (.not. allocated (sampler%x)) &
         allocate (sampler%x (size (x_in), 2))
    if (.not. allocated (sampler%f)) &
         allocate (sampler%f (2))
    select case (c)
    case (1)
       u = x_in(1)
       xx = u ** 0.2_default
       v = (1 - xx) ** 5._default
    case (2)
       v = x_in(1)
       xx = 1 - v ** 0.2_default
       u = xx ** 5._default
    end select
    sampler%val = sampler%a * 5 * xx ** 4 + sampler%b * 5 * (1 - xx) ** 4
    sampler%f(1) = 0.2_default * u ** (-0.8_default)
    sampler%f(2) = 0.2_default * v ** (-0.8_default)
    sampler%x(:,1) = [u]
    sampler%x(:,2) = [v]
  end subroutine test_sampler_3_compute

  subroutine test_sampler_3_evaluate (sampler, c, x_in, val, x, f)
    class(test_sampler_3_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    call sampler%compute (c, x_in)
    call sampler%fetch (val, x, f)
  end subroutine test_sampler_3_evaluate

  function test_sampler_3_is_valid (sampler) result (valid)
    class(test_sampler_3_t), intent(in) :: sampler
    logical :: valid
    valid = .true.
  end function test_sampler_3_is_valid

  subroutine test_sampler_3_rebuild (sampler, c, x_in, val, x, f)
    class(test_sampler_3_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(in) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    call sampler%compute (c, x_in)
    x = sampler%x
    f = sampler%f
  end subroutine test_sampler_3_rebuild

  subroutine test_sampler_3_fetch (sampler, val, x, f)
    class(test_sampler_3_t), intent(in) :: sampler
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    val = sampler%val
    x = sampler%x
    f = sampler%f
  end subroutine test_sampler_3_fetch

  subroutine mci_vamp_1 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_1"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(single channel)"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 1)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

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
    write (u, "(A)")  "* Integrate with n_calls = 1000"
    write (u, "(A)")  "   (lower precision to avoid"
    write (u, "(A)")  "      numerical noise)"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass ()
    end select
    call mci%integrate (mci_instance, sampler, 1, 1000, pacify = .true.)
    call mci%write (u, .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u, .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_1"

  end subroutine mci_vamp_1

  subroutine mci_vamp_2 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_2"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(single channel)"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 1)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_1_t :: sampler)
    select type (sampler)
    type is (test_sampler_1_t)
       sampler%mode = 2
    end select
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 100"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .false.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 100)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_2"

  end subroutine mci_vamp_2

  subroutine mci_vamp_3 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_3"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(single channel)"
    write (u, "(A)")  "*            and adapt grid"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 1)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_1_t :: sampler)
    select type (sampler)
    type is (test_sampler_1_t)
       sampler%mode = 2
    end select
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 100"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 100)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_3"

  end subroutine mci_vamp_3

  subroutine mci_vamp_4 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_4"
    write (u, "(A)")  "*   Purpose: integrate function in two dimensions &
         &(single channel)"
    write (u, "(A)")  "*            and adapt grid"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (2, 1)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_1_t :: sampler)
    select type (sampler)
    type is (test_sampler_1_t)
       sampler%mode = 3
    end select
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 1000"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_4"

  end subroutine mci_vamp_4

  subroutine mci_vamp_5 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_5"
    write (u, "(A)")  "*   Purpose: integrate function in two dimensions &
         &(two channels)"
    write (u, "(A)")  "*            and adapt grid"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_2_t :: sampler)
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 1000"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_5"

  end subroutine mci_vamp_5

  subroutine mci_vamp_6 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_6"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(two channels)"
    write (u, "(A)")  "*            and adapt weights"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_3_t :: sampler)
    select type (sampler)
    type is (test_sampler_3_t)
       sampler%a = 0.9_default
       sampler%b = 0.1_default
    end select
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 1000"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_weights = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()
    deallocate (mci_instance)
    deallocate (mci)

    write (u, "(A)")
    write (u, "(A)")  "* Re-initialize with chained channels"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 2)
    call mci%declare_chains ([1,1])
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 1000"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_weights = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_6"

  end subroutine mci_vamp_6

  subroutine mci_vamp_7 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    type(phs_channel_t), dimension(:), allocatable :: channel
    class(rng_t), allocatable :: rng
    real(default), dimension(:,:), allocatable :: x
    integer :: u_grid, iostat, i, div, ch
    character(16) :: buffer

    write (u, "(A)")  "* Test output: mci_vamp_7"
    write (u, "(A)")  "*   Purpose: check effect of channel equivalences"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_3_t :: sampler)
    select type (sampler)
    type is (test_sampler_3_t)
       sampler%a = 0.7_default
       sampler%b = 0.3_default
    end select
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 2 and n_calls = 1000, &
         &adapt grids"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 2, 1000)

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Write grids and extract binning"
    write (u, "(A)")

    u_grid = free_unit ()
    open (u_grid, status = "scratch", action = "readwrite")
    select type (mci_instance)
    type is (mci_vamp_instance_t)
       call vamp_write_grids (mci_instance%grids, u_grid)
    end select
    rewind (u_grid)
    allocate (x (0:20, 2))
    do div = 1, 2
       FIND_BINS1: do
          read (u_grid, "(A)")  buffer
          if (trim (adjustl (buffer)) == "begin d%x") then
             do
                read (u_grid, *, iostat = iostat)  i, x(i,div)
                if (iostat /= 0)  exit FIND_BINS1
             end do
          end if
       end do FIND_BINS1
    end do
    close (u_grid)

    write (u, "(1x,A,L1)")  "Equal binning in both channels = ", &
         all (x(:,1) == x(:,2))
    deallocate (x)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()
    deallocate (mci_instance)
    deallocate (mci)

    write (u, "(A)")
    write (u, "(A)")  "* Re-initialize integrator, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .true.
       call mci%set_grid_parameters (grid_par)
    end select

    write (u, "(A)")  "* Define equivalences"
    write (u, "(A)")

    allocate (channel (2))
    do ch = 1, 2
       allocate (channel(ch)%eq (2))
       do i = 1, 2
          associate (eq => channel(ch)%eq(i))
            call eq%init (1)
            eq%c = i
            eq%perm = [1]
            eq%mode = [0]
          end associate
       end do
       write (u, "(1x,I0,':')", advance = "no")  ch
       call channel(ch)%write (u)
    end do
    call mci%declare_equivalences (channel, dim_offset = 0)

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 2 and n_calls = 1000, &
         &adapt grids"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 2, 1000)

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Write grids and extract binning"
    write (u, "(A)")

    u_grid = free_unit ()
    open (u_grid, status = "scratch", action = "readwrite")
    select type (mci_instance)
    type is (mci_vamp_instance_t)
       call vamp_write_grids (mci_instance%grids, u_grid)
    end select
    rewind (u_grid)
    allocate (x (0:20, 2))
    do div = 1, 2
       FIND_BINS2: do
          read (u_grid, "(A)")  buffer
          if (trim (adjustl (buffer)) == "begin d%x") then
             do
                read (u_grid, *, iostat = iostat)  i, x(i,div)
                if (iostat /= 0)  exit FIND_BINS2
             end do
          end if
       end do FIND_BINS2
    end do
    close (u_grid)

    write (u, "(1x,A,L1)")  "Equal binning in both channels = ", &
         all (x(:,1) == x(:,2))
    deallocate (x)


    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_7"

  end subroutine mci_vamp_7

  subroutine mci_vamp_8 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_8"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(two channels)"
    write (u, "(A)")  "*            in three passes"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_3_t :: sampler)
    select type (sampler)
    type is (test_sampler_3_t)
       sampler%a = 0.9_default
       sampler%b = 0.1_default
    end select
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with grid and weight adaptation"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true., adapt_weights = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with grid adaptation"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate without adaptation"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass ()
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_8"

  end subroutine mci_vamp_8

  subroutine mci_vamp_9 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_9"
    write (u, "(A)")  "*   Purpose: integrate function in two dimensions &
         &(two channels)"
    write (u, "(A)")  "*            and generate a weighted event"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_2_t :: sampler)
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 1000"
    write (u, "(A)")

    call mci%add_pass ()
    call mci%integrate (mci_instance, sampler, 1, 1000)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate a weighted event"
    write (u, "(A)")

    call mci_instance%init_simulation ()
    call mci%generate_weighted_event (mci_instance, sampler)

    write (u, "(1x,A)")  "MCI instance:"
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final_simulation ()
    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_9"

  end subroutine mci_vamp_9

  subroutine mci_vamp_10 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng
    type(string_t) :: file1, file2
    character(80) :: buffer1, buffer2
    integer :: u1, u2, iostat1, iostat2
    logical :: equal, success

    write (u, "(A)")  "* Test output: mci_vamp_10"
    write (u, "(A)")  "*   Purpose: write and read VAMP grids"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    mci%md5sum = "1234567890abcdef1234567890abcdef"

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_2_t :: sampler)
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 1000"
    write (u, "(A)")

    call mci%add_pass ()
    call mci%integrate (mci_instance, sampler, 1, 1000)

    write (u, "(A)")  "* Write grids to file"
    write (u, "(A)")

    file1 = "mci_vamp_10.1"
    select type (mci)
    type is (mci_vamp_t)
       call mci%set_grid_filename (file1)
       call mci%write_grids (mci_instance)
    end select

    call mci_instance%final ()
    call mci%final ()
    deallocate (mci)

    write (u, "(A)")  "* Read grids from file"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_vamp_t)
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    mci%md5sum = "1234567890abcdef1234567890abcdef"

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    select type (mci)
    type is (mci_vamp_t)
       call mci%set_grid_filename (file1)
       call mci%add_pass ()
       call mci%current_pass%configure (1, 1000, &
            mci%min_calls, &
            mci%grid_par%min_bins, mci%grid_par%max_bins, &
            mci%grid_par%min_calls_per_channel * mci%n_channel)
       call mci%read_grids_header (success)
       call mci%compute_md5sum ()
       call mci%read_grids_data (mci_instance, read_integrals = .true.)
    end select
    write (u, "(1x,A,L1)")  "success = ", success

    write (u, "(A)")
    write (u, "(A)")  "* Write grids again"
    write (u, "(A)")

    file2 = "mci_vamp_10.2"
    select type (mci)
    type is (mci_vamp_t)
       call mci%set_grid_filename (file2)
       call mci%write_grids (mci_instance)
    end select

    u1 = free_unit ()
    open (u1, file = char (file1) // ".vg", action = "read", status = "old")
    u2 = free_unit ()
    open (u2, file = char (file2) // ".vg", action = "read", status = "old")

    equal = .true.
    iostat1 = 0
    iostat2 = 0
    do while (equal .and. iostat1 == 0 .and. iostat2 == 0)
       read (u1, "(A)", iostat = iostat1)  buffer1
       read (u2, "(A)", iostat = iostat2)  buffer2
       equal = buffer1 == buffer2 .and. iostat1 == iostat2
    end do
    close (u1)
    close (u2)

    if (equal) then
       write (u, "(1x,A)")  "Success: grid files are identical"
    else
       write (u, "(1x,A)")  "Failure: grid files differ"
    end if

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_10"

  end subroutine mci_vamp_10

  subroutine mci_vamp_11 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_11"
    write (u, "(A)")  "*   Purpose: integrate function in two dimensions &
         &(two channels)"
    write (u, "(A)")  "*            and generate a weighted event"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
       call mci%set_grid_filename (var_str ("mci_vamp_11"))
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_2_t :: sampler)

    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 1000"
    write (u, "(A)")

    call mci%add_pass ()
    call mci%integrate (mci_instance, sampler, 1, 1000)

    write (u, "(A)")  "* Reset instance"
    write (u, "(A)")

    call mci_instance%final ()
    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Generate a weighted event"
    write (u, "(A)")

    call mci_instance%init_simulation ()
    call mci%generate_weighted_event (mci_instance, sampler)

    write (u, "(A)")  "* Cleanup"

    call mci_instance%final_simulation ()
    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_11"

  end subroutine mci_vamp_11

  subroutine mci_vamp_12 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_12"
    write (u, "(A)")  "*   Purpose: integrate function in two dimensions &
         &(two channels)"
    write (u, "(A)")  "*            and generate an unweighted event"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
       call mci%set_grid_filename (var_str ("mci_vamp_12"))
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_2_t :: sampler)

    write (u, "(A)")  "* Integrate with n_it = 3 and n_calls = 1000"
    write (u, "(A)")

    call mci%add_pass ()
    call mci%integrate (mci_instance, sampler, 1, 1000)

    write (u, "(A)")  "* Reset instance"
    write (u, "(A)")

    call mci_instance%final ()
    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Generate an unweighted event"
    write (u, "(A)")

    call mci_instance%init_simulation ()
    call mci%generate_unweighted_event (mci_instance, sampler)

    write (u, "(1x,A)")  "MCI instance:"
    call mci_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final_simulation ()
    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_12"

  end subroutine mci_vamp_12

  subroutine mci_vamp_13 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci, mci_ref
    logical :: success

    write (u, "(A)")  "* Test output: mci_vamp_13"
    write (u, "(A)")  "*   Purpose: match and update integrators"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator with no passes"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (2, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
    end select
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize reference"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci_ref)
    call mci_ref%set_dimensions (2, 2)
    select type (mci_ref)
    type is (mci_vamp_t)
       call mci_ref%set_grid_parameters (grid_par)
    end select

    select type (mci_ref)
    type is (mci_vamp_t)
       call mci_ref%add_pass (adapt_grids = .true.)
       call mci_ref%current_pass%configure (2, 1000, 0, 1, 5, 0)
       mci_ref%current_pass%calls = [77, 77]
       mci_ref%current_pass%integral = [1.23_default, 3.45_default]
       mci_ref%current_pass%error = [0.23_default, 0.45_default]
       mci_ref%current_pass%efficiency = [0.1_default, 0.6_default]
       mci_ref%current_pass%integral_defined = .true.

       call mci_ref%add_pass ()
       call mci_ref%current_pass%configure (2, 2000, 0, 1, 7, 0)
       mci_ref%current_pass%calls = [99, 0]
       mci_ref%current_pass%integral = [7.89_default, 0._default]
       mci_ref%current_pass%error = [0.89_default, 0._default]
       mci_ref%current_pass%efficiency = [0.86_default, 0._default]
       mci_ref%current_pass%integral_defined = .true.
    end select

    call mci_ref%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Update integrator (no-op, should succeed)"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%update_from_ref (mci_ref, success)
    end select

    write (u, "(1x,A,L1)")  "success = ", success
    write (u, "(A)")
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Add pass to integrator"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true.)
       call mci%current_pass%configure (2, 1000, 0, 1, 5, 0)
       mci%current_pass%calls = [77, 77]
       mci%current_pass%integral = [1.23_default, 3.45_default]
       mci%current_pass%error = [0.23_default, 0.45_default]
       mci%current_pass%efficiency = [0.1_default, 0.6_default]
       mci%current_pass%integral_defined = .true.
    end select

    write (u, "(A)")  "* Update integrator (no-op, should succeed)"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%update_from_ref (mci_ref, success)
    end select

    write (u, "(1x,A,L1)")  "success = ", success
    write (u, "(A)")
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Add pass to integrator, wrong parameters"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass ()
       call mci%current_pass%configure (2, 1000, 0, 1, 7, 0)
    end select

    write (u, "(A)")  "* Update integrator (should fail)"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%update_from_ref (mci_ref, success)
    end select

    write (u, "(1x,A,L1)")  "success = ", success
    write (u, "(A)")
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Reset and add passes to integrator"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%reset ()
       call mci%add_pass (adapt_grids = .true.)
       call mci%current_pass%configure (2, 1000, 0, 1, 5, 0)
       mci%current_pass%calls = [77, 77]
       mci%current_pass%integral = [1.23_default, 3.45_default]
       mci%current_pass%error = [0.23_default, 0.45_default]
       mci%current_pass%efficiency = [0.1_default, 0.6_default]
       mci%current_pass%integral_defined = .true.

       call mci%add_pass ()
       call mci%current_pass%configure (2, 2000, 0, 1, 7, 0)
    end select

    write (u, "(A)")  "* Update integrator (should succeed)"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%update_from_ref (mci_ref, success)
    end select

    write (u, "(1x,A,L1)")  "success = ", success
    write (u, "(A)")
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Update again (no-op, should succeed)"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%update_from_ref (mci_ref, success)
    end select

    write (u, "(1x,A,L1)")  "success = ", success
    write (u, "(A)")
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Add extra result to integrator"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       mci%current_pass%calls(2) = 1234
    end select

    write (u, "(A)")  "* Update integrator (should fail)"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%update_from_ref (mci_ref, success)
    end select

    write (u, "(1x,A,L1)")  "success = ", success
    write (u, "(A)")
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci%final ()
    call mci_ref%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_13"

  end subroutine mci_vamp_13

  subroutine mci_vamp_14 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_14"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(single channel)"
    write (u, "(A)")  "*            and check accuracy goal"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 1)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%use_vamp_equivalences = .false.
       grid_par%accuracy_goal = 5E-2_default
       call mci%set_grid_parameters (grid_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_1_t :: sampler)
    select type (sampler)
    type is (test_sampler_1_t)
       sampler%mode = 2
    end select
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_it = 5 and n_calls = 100"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 5, 100)
    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_14"

  end subroutine mci_vamp_14

  subroutine mci_vamp_15 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    type(history_parameters_t) :: history_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_15"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(two channels)"
    write (u, "(A)")  "*            in three passes, show history"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    history_par%channel = .true.

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 2)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%stratified = .false.
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
       call mci%set_history_parameters (history_par)
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    allocate (test_sampler_3_t :: sampler)
    select type (sampler)
    type is (test_sampler_3_t)
       sampler%a = 0.9_default
       sampler%b = 0.1_default
    end select
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Pass 1: grid and weight adaptation"

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true., adapt_weights = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)

    write (u, "(A)")
    write (u, "(A)")  "* Pass 2: grid adaptation"

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass (adapt_grids = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)

    write (u, "(A)")
    write (u, "(A)")  "* Pass 3: without adaptation"

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass ()
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of MCI record, with history"
    write (u, "(A)")

    call mci%write (u)
    select type (mci)
    type is (mci_vamp_t)
       call mci%write_history (u)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_15"

  end subroutine mci_vamp_15

  subroutine mci_vamp_16 (u)
    integer, intent(in) :: u
    type(grid_parameters_t) :: grid_par
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng

    write (u, "(A)")  "* Test output: mci_vamp_16"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(single channel)"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator"
    write (u, "(A)")

    allocate (mci_vamp_t :: mci)
    call mci%set_dimensions (1, 1)
    select type (mci)
    type is (mci_vamp_t)
       grid_par%use_vamp_equivalences = .false.
       call mci%set_grid_parameters (grid_par)
       mci%negative_weights = .true.
    end select

    allocate (rng_tao_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")  "* Initialize test sampler"
    write (u, "(A)")

    allocate (test_sampler_1_t :: sampler)
    select type (sampler)
    type is (test_sampler_1_t)
       sampler%mode = 4
    end select
    call sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_calls = 1000"
    write (u, "(A)")  "   (lower precision to avoid"
    write (u, "(A)")  "      numerical noise)"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp_t)
       call mci%add_pass ()
    end select
    call mci%integrate (mci_instance, sampler, 1, 1000, pacify = .true.)
    call mci%write (u, .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u, .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp_16"

  end subroutine mci_vamp_16


end module mci_vamp_uti
