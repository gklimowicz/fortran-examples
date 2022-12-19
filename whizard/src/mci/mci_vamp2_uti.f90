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

module mci_vamp2_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  use io_units
  use constants, only: PI, TWOPI
  use rng_base
  use rng_tao
  use rng_stream
  use mci_base

  use mci_vamp2

  implicit none
  private

  public :: mci_vamp2_1
  public :: mci_vamp2_2
  public :: mci_vamp2_3

  type, extends (mci_sampler_t) :: test_sampler_1_t
     real(default), dimension(:), allocatable :: x
     real(default) :: val
     integer :: mode = 1
   contains
     procedure, public :: write => test_sampler_1_write
     procedure, public :: evaluate => test_sampler_1_evaluate
     procedure, public :: is_valid => test_sampler_1_is_valid
     procedure, public :: rebuild => test_sampler_1_rebuild
     procedure, public :: fetch => test_sampler_1_fetch
  end type test_sampler_1_t

  type, extends (mci_sampler_t) :: test_sampler_2_t
     real(default), dimension(:,:), allocatable :: x
     real(default), dimension(:), allocatable :: f
     real(default) :: val
   contains
     procedure, public :: write => test_sampler_2_write
     procedure, public :: compute => test_sampler_2_compute
     procedure, public :: evaluate => test_sampler_2_evaluate
     procedure, public :: is_valid => test_sampler_2_is_valid
     procedure, public :: rebuild => test_sampler_2_rebuild
     procedure, public :: fetch => test_sampler_2_fetch
  end type test_sampler_2_t


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

  subroutine mci_vamp2_1 (u)
    integer, intent(in) :: u
    type(mci_vamp2_config_t) :: config
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable, target :: mci_sampler
    class(rng_t), allocatable :: rng
    type(string_t) :: filename

    write (u, "(A)") "* Test output: mci_vamp2_1"
    write (u, "(A)") "*   Purpose: integrate function in one dimension (single channel)"

    write (u, "(A)")
    write (u, "(A)") "* Initialise integrator"
    write (u, "(A)")

    allocate (mci_vamp2_t :: mci)
    call mci%set_dimensions (1, 1)

    filename = "mci_vamp2_1"
    select type (mci)
    type is (mci_vamp2_t)
       call mci%set_config (config)
       call mci%set_grid_filename (filename)
    end select

    allocate (rng_stream_t :: rng)
    call rng%init ()
    call mci%import_rng (rng)

    call mci%write (u, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)") "* Initialise instance"
    write (u, "(A)")

    call mci%allocate_instance (mci_instance)
    call mci_instance%init (mci)

    write (u, "(A)")
    write (u, "(A)") "* Initialise test sampler"
    write (u, "(A)")

    allocate (test_sampler_1_t :: mci_sampler)
    call mci_sampler%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Integrate with n_calls = 1000"
    write (u, "(A)")  "   (lower precision to avoid"
    write (u, "(A)")  "      numerical noise)"
    write (u, "(A)")

    select type (mci)
    type is (mci_vamp2_t)
       call mci%add_pass ()
    end select
    call mci%integrate (mci_instance, mci_sampler, 1, 1000, pacify = .true.)
    call mci%write (u, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)") "* Dump channel weights and grids to file"
    write (u, "(A)")

    mci%md5sum = "1234567890abcdef1234567890abcdef"
    select type (mci)
    type is (mci_vamp2_t)
       call mci%write_grids ()
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp2_1"

  end subroutine mci_vamp2_1

  subroutine mci_vamp2_2 (u)
    type(mci_vamp2_config_t) :: config
    integer, intent(in) :: u
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng
    type(string_t) :: filename

    write (u, "(A)")  "* Test output: mci_vamp2_2"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(single channel), but multiple iterations."

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp2_t :: mci)
    call mci%set_dimensions (1, 1)
    filename = "mci_vamp2_2"
    select type (mci)
    type is (mci_vamp2_t)
       call mci%set_config (config)
       call mci%set_grid_filename (filename)
    end select

    allocate (rng_stream_t :: rng)
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
    type is (mci_vamp2_t)
       call mci%add_pass (adapt_grids = .false.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000, pacify = .true.)
    call mci%write (u, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)") "* Dump channel weights and grids to file"
    write (u, "(A)")

    mci%md5sum = "1234567890abcdef1234567890abcdef"
    select type (mci)
    type is (mci_vamp2_t)
       call mci%write_grids ()
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp2_2"

  end subroutine mci_vamp2_2

  subroutine mci_vamp2_3 (u)
    integer, intent(in) :: u
    type(mci_vamp2_config_t) :: config
    class(mci_t), allocatable, target :: mci
    class(mci_instance_t), pointer :: mci_instance => null ()
    class(mci_sampler_t), allocatable :: sampler
    class(rng_t), allocatable :: rng
    type(string_t) :: filename

    write (u, "(A)")  "* Test output: mci_vamp2_3"
    write (u, "(A)")  "*   Purpose: integrate function in one dimension &
         &(single channel)"
    write (u, "(A)")  "*            and adapt grid"

    write (u, "(A)")
    write (u, "(A)")  "* Initialize integrator, sampler, instance"
    write (u, "(A)")

    allocate (mci_vamp2_t :: mci)
    call mci%set_dimensions (1, 1)
    filename = "mci_vamp2_3"
    select type (mci)
    type is (mci_vamp2_t)
       call mci%set_grid_filename (filename)
       call mci%set_config (config)
    end select

    allocate (rng_stream_t :: rng)
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
    type is (mci_vamp2_t)
       call mci%add_pass (adapt_grids = .true.)
    end select
    call mci%integrate (mci_instance, sampler, 3, 1000, pacify = .true.)
    call mci%write (u, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Contents of mci_instance:"
    write (u, "(A)")

    call mci_instance%write (u, pacify = .true.)

    write (u, "(A)")
    write (u, "(A)") "* Dump channel weights and grids to file"
    write (u, "(A)")

    mci%md5sum = "1234567890abcdef1234567890abcdef"
    select type (mci)
    type is (mci_vamp2_t)
       call mci%write_grids ()
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    call mci_instance%final ()
    call mci%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: mci_vamp2_3"

  end subroutine mci_vamp2_3


end module mci_vamp2_uti
