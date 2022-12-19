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

module vamp2_uti

  use kinds, only: default
  use io_units
  use constants, only: pi
  use numeric_utils, only: nearly_equal
  use format_defs, only: FMT_12
  use rng_base
  use rng_stream
  use vegas, only: vegas_func_t, vegas_grid_t, operator(==)
  use vamp2

  implicit none
  private

  public :: vamp2_1
  public :: vamp2_2
  public :: vamp2_3
  public :: vamp2_4
  public :: vamp2_5

   type, extends (vamp2_func_t) :: vamp2_test_func_t
     !
   contains
     procedure, public :: evaluate_maps => vamp2_test_func_evaluate_maps
     procedure, public :: evaluate_func => vamp2_test_func_evaluate
end type vamp2_test_func_t

  type, extends(vamp2_func_t) :: vamp2_test_func_2_t
     !
   contains
     procedure :: evaluate_maps => vamp2_test_func_2_evaluate_maps
     procedure :: evaluate_func => vamp2_test_func_2_evaluate_func
  end type vamp2_test_func_2_t

  type, extends(vamp2_func_t) :: vamp2_test_func_3_t
     !
   contains
     procedure :: evaluate_maps => vamp2_test_func_3_evaluate_maps
     procedure :: evaluate_func => vamp2_test_func_3_evaluate_func
  end type vamp2_test_func_3_t


contains
  subroutine vamp2_test_func_evaluate_maps (self, x)
    class(vamp2_test_func_t), intent(inout) :: self
    real(default), dimension(:), intent(in) :: x
    self%xi(:, 1) = x
    self%det(1) = 1
    self%valid_x = .true.
  end subroutine vamp2_test_func_evaluate_maps

  real(default) function vamp2_test_func_evaluate (self, x) result (f)
    class(vamp2_test_func_t), intent(in) :: self
    real(default), dimension(:), intent(in) :: x
    f = 1.0 / (pi**3)
    f = f / ( 1.0 - cos (x(1)) * cos (x(2)) * cos (x(3)))
  end function vamp2_test_func_evaluate

  subroutine vamp2_test_func_2_evaluate_maps (self, x)
    class(vamp2_test_func_2_t), intent(inout) :: self
    real(default), dimension(:), intent(in) :: x
    select case (self%current_channel)
    case (1)
       self%xi(:, 1) = x
       self%xi(1, 2) = x(1) * x(2)
       self%xi(2, 2) = 0.5 * ( 1. + log(x(1) / x(2)) / log(x(1) * x(2)))
    case (2)
       self%xi(1, 1) = x(1)**x(2)
       self%xi(2, 1) = x(1)**(1. - x(2))
       self%xi(:, 2) = x
    end select
    self%det(1) = 1.
    self%det(2) = abs (log(self%xi(1, 2)))
    self%valid_x = .true.
  end subroutine vamp2_test_func_2_evaluate_maps

  real(default) function vamp2_test_func_2_evaluate_func (self, x) result (f)
    class(vamp2_test_func_2_t), intent(in) :: self
    real(default), dimension(:), intent(in) :: x
    f = 4. * sin(pi * self%xi(1, 1))**2 * sin(pi * self%xi(2, 1))**2 &
         + 2. * sin(pi * self%xi(2, 2))**2
  end function vamp2_test_func_2_evaluate_func

  subroutine vamp2_test_func_3_evaluate_maps (self, x)
    class(vamp2_test_func_3_t), intent(inout) :: self
    real(default), dimension(:), intent(in) :: x
    real(default) :: u, v, xx
    select case (self%current_channel)
    case (1)
       u = x(1)
       xx = u**0.2_default
       v = (1 - xx)**5._default
    case (2)
       v = x(1)
       xx = 1 - v**0.2_default
       u = xx**5._default
    end select
    self%det(1) = 0.2_default * u**(-0.8_default)
    self%det(2) = 0.2_default * v**(-0.8_default)
    self%xi(:, 1) = [u]
    self%xi(:, 2) = [v]
    self%valid_x = .true.
  end subroutine vamp2_test_func_3_evaluate_maps

  real(default) function vamp2_test_func_3_evaluate_func (self, x) result (f)
    class(vamp2_test_func_3_t), intent(in) :: self
    real(default), dimension(:), intent(in) :: x
    real(default) :: xx
    select case (self%current_channel)
    case (1)
       xx = x(1)**0.2_default
    case (2)
       xx = 1 - x(1)**0.2_default
    end select
    f = 5 * xx**4 + 5 * (1 - xx)**4
  end function vamp2_test_func_3_evaluate_func

  subroutine vamp2_1 (u)
    integer, intent(in) :: u
    type(vamp2_t) :: mc_integrator
    class(rng_t), allocatable :: rng
    class(vamp2_func_t), allocatable :: func
    real(default), dimension(3), parameter :: x_lower = 0., &
         x_upper = pi
    real(default) :: result, abserr

    write (u, "(A)") "* Test output: vamp2_1"
    write (u, "(A)") "*   Purpose: initialise the VAMP2 MC integrator and the grid"
    write (u, "(A)")

    write (u, "(A)") "* Initialise random number generator (default seed)"
    write (u, "(A)")

    allocate (rng_stream_t :: rng)
    call rng%init ()

    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise MC integrator with n_channel = 1 and n_dim = 3"
    write (u, "(A)")

    allocate (vamp2_test_func_t :: func)
    call func%init (n_dim = 3, n_channel = 1)
    mc_integrator = vamp2_t (1, 3)
    call mc_integrator%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise grid with n_calls = 10000"
    write (u, "(A)")

    call mc_integrator%set_limits (x_lower, x_upper)
    call mc_integrator%set_calls (10000)

    write (u, "(A)")
    write (u, "(A)") "* Integrate with n_it = 3 and n_calls = 10000 (Adaptation)"
    write (u, "(A)")

    call mc_integrator%integrate (func, rng, 3, result=result, abserr=abserr)
    write (u, "(2x,A," // FMT_12 // ",A," // FMT_12 // ")") "Result: ", result, " +/- ", abserr

    write (u, "(A)")
    write (u, "(A)") "* Integrate with n_it = 3 and n_calls = 2000 (Precision)"
    write (u, "(A)")

    call mc_integrator%set_calls (2000)
    call mc_integrator%integrate (func, rng, 3, result=result, abserr=abserr)
    write (u,  "(2x,A," // FMT_12 // ",A," // FMT_12 // ")") "Result: ", result, " +/- ", abserr

    write (u, "(A)")
    write (u, "(A)") "* Cleanup"

    call mc_integrator%final ()
    call rng%final ()
    deallocate (rng)
  end subroutine vamp2_1

  subroutine vamp2_2 (u)
    integer, intent(in) :: u
    type(vamp2_t) :: mc_integrator
    class(rng_t), allocatable :: rng
    class(vamp2_func_t), allocatable :: func
    real(default), dimension(2), parameter :: x_lower = 0., &
         x_upper = 1.
    real(default) :: result, abserr

    write (u, "(A)") "* Test output: vamp2_2"
    write (u, "(A)") "*   Purpose:  intgeration of two-dimensional &
       & function with two channels"
    write (u, "(A)")

    write (u, "(A)") "* Initialise random number generator (default seed)"
    write (u, "(A)")

    allocate (rng_stream_t :: rng)
    call rng%init ()

    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise MC integrator with n_channel = 1 and n_dim = 3"
    write (u, "(A)")

    allocate (vamp2_test_func_2_t :: func)
    call func%init (n_dim = 2, n_channel = 2)
    mc_integrator = vamp2_t (2, 2)
    call mc_integrator%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise grid with n_calls = 10000"
    write (u, "(A)")

    call mc_integrator%set_limits (x_lower, x_upper)
    call mc_integrator%set_calls (1000)

    write (u, "(A)")
    write (u, "(A)") "* Integrate with n_it = 3 and n_calls = 10000 (Adaptation)"
    write (u, "(A)")

    call mc_integrator%integrate (func, rng, 3, verbose = .true., result=result, abserr=abserr)
    write (u,  "(2x,A," // FMT_12 // ",A," // FMT_12 // ")") "Result: ", result, " +/- ", abserr

    write (u, "(A)")
    write (u, "(A)") "* Integrate with n_it = 3 and n_calls = 2000 (Precision)"
    write (u, "(A)")

    call mc_integrator%set_calls (200)
    call mc_integrator%integrate (func, rng, 3, verbose = .true., result=result, abserr=abserr)
    write (u,  "(2x,A," // FMT_12 // ",A," // FMT_12 // ")") "Result: ", result, " +/- ", abserr

    write (u, "(A)")
    write (u, "(A)") "* Cleanup"

    call mc_integrator%final ()
    call rng%final ()
    deallocate (rng)
  end subroutine vamp2_2

  subroutine vamp2_3 (u)
    integer, intent(in) :: u
    type(vamp2_t) :: mc_integrator
    class(rng_t), allocatable :: rng
    class(vamp2_func_t), allocatable :: func
    real(default), dimension(2), parameter :: x_lower = 0., &
         x_upper = 1.
    real(default) :: result, abserr
    integer :: unit

    write (u, "(A)") "* Test output: vamp2_3"
    write (u, "(A)") "*   Purpose:  intgeration of two-dimensional &
       & function with two channels"
    write (u, "(A)")

    write (u, "(A)") "* Initialise random number generator (default seed)"
    write (u, "(A)")

    allocate (rng_stream_t :: rng)
    call rng%init ()

    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise MC integrator with n_channel = 1 and n_dim = 3"
    write (u, "(A)")

    allocate (vamp2_test_func_2_t :: func)
    call func%init (n_dim = 2, n_channel = 2)
    mc_integrator = vamp2_t (2, 2)
    call mc_integrator%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise grid with n_calls = 20000"
    write (u, "(A)")

    call mc_integrator%set_limits (x_lower, x_upper)
    call mc_integrator%set_calls (20000)

    write (u, "(A)")
    write (u, "(A)") "* Integrate with n_it = 3 and n_calls = 20000 (Adaptation)"
    write (u, "(A)")

    call mc_integrator%integrate (func, rng, 3, result=result, abserr=abserr)
    write (u,  "(2x,A," // FMT_12 // ",A," // FMT_12 // ")") "Result: ", result, " +/- ", abserr

    write (u, "(A)")
    write (u, "(A)") "* Write grid to file vamp2_3.grids"
    write (u, "(A)")

    unit = free_unit ()
    open (unit, file = "vamp2_3.grids", &
         action = "write", status = "replace")
    call mc_integrator%write_grids (unit)
    close (unit)

    write (u, "(A)")
    write (u, "(A)") "* Read grid from file vamp2_3.grids"
    write (u, "(A)")

    call mc_integrator%final ()

    unit = free_unit ()
    open (unit, file = "vamp2_3.grids", &
         action = "read", status = "old")
    call mc_integrator%read_grids (unit)
    close (unit)

    write (u, "(A)")
    write (u, "(A)") "* Integrate with n_it = 3 and n_calls = 5000 (Precision)"
    write (u, "(A)")

    call mc_integrator%set_calls (5000)
    call mc_integrator%integrate (func, rng, 3, result=result, abserr=abserr)
    write (u,  "(2x,A," // FMT_12 // ",A," // FMT_12 // ")") "Result: ", result, " +/- ", abserr

    write (u, "(A)")
    write (u, "(A)") "* Cleanup"

    call mc_integrator%final ()
    call rng%final ()
    deallocate (rng)
  end subroutine vamp2_3

  subroutine vamp2_4 (u)
    integer, intent(in) :: u
    type(vamp2_t) :: mc_integrator
    class(rng_t), allocatable :: rng
    class(vamp2_func_t), allocatable :: func
    real(default), dimension(2), parameter :: x_lower = 0., &
         x_upper = 1.
    real(default) :: result, abserr
    integer :: unit

    write (u, "(A)") "* Test output: vamp2_4"
    write (u, "(A)") "*   Purpose:  intgeration of two-dimensional &
       & function with two channels with chains"
    write (u, "(A)")

    write (u, "(A)") "* Initialise random number generator (default seed)"
    write (u, "(A)")

    allocate (rng_stream_t :: rng)
    call rng%init ()

    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise MC integrator with n_channel = 2 and n_dim = 2"
    write (u, "(A)")

    allocate (vamp2_test_func_2_t :: func)
    call func%init (n_dim = 2, n_channel = 2)
    mc_integrator = vamp2_t (2, 2)
    call mc_integrator%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise grid with n_calls = 20000 and set chains"
    write (u, "(A)")

    call mc_integrator%set_limits (x_lower, x_upper)
    call mc_integrator%set_calls (20000)
    call mc_integrator%set_chain (2, [1, 2])

    write (u, "(A)")
    write (u, "(A)") "* Integrate with n_it = 3 and n_calls = 10000 (Adaptation)"
    write (u, "(A)")

    call mc_integrator%integrate (func, rng, 3, result=result, abserr=abserr)
    write (u,  "(2x,A," // FMT_12 // ",A," // FMT_12 // ")") "Result: ", result, " +/- ", abserr

    write (u, "(A)")
    write (u, "(A)") "* Write grid to file vamp2_4.grids"
    write (u, "(A)")

    unit = free_unit ()
    open (unit, file = "vamp2_4.grids", &
         action = "write", status = "replace")
    call mc_integrator%write_grids (unit)
    close (unit)

    write (u, "(A)")
    write (u, "(A)") "* Read grid from file vamp2_4.grids"
    write (u, "(A)")

    call mc_integrator%final ()

    unit = free_unit ()
    open (unit, file = "vamp2_4.grids", &
         action = "read", status = "old")
    call mc_integrator%read_grids (unit)
    close (unit)

    write (u, "(A)")
    write (u, "(A)") "* Integrate with n_it = 3 and n_calls = 5000 (Precision)"
    write (u, "(A)")

    call mc_integrator%set_calls (5000)
    call mc_integrator%integrate (func, rng, 3, result=result, abserr=abserr)
    write (u,  "(2x,A," // FMT_12 // ",A," // FMT_12 // ")") "Result: ", result, " +/- ", abserr

    write (u, "(A)")
    write (u, "(A)") "* Cleanup"

    call mc_integrator%final ()
    call rng%final ()
    deallocate (rng)
  end subroutine vamp2_4

  subroutine vamp2_5 (u)
    integer, intent(in) :: u
    type(vamp2_t) :: mc_integrator
    class(rng_t), allocatable :: rng
    class(vamp2_func_t), allocatable :: func
    real(default), dimension(1), parameter :: x_lower = 0., &
         x_upper = 1.
    real(default) :: result, abserr
    integer :: unit
    type(vamp2_config_t) :: config
    type(vamp2_equivalences_t) :: eqv
    type(vegas_grid_t), dimension(2) :: grid

    write (u, "(A)") "* Test output: vamp2_5"
    write (u, "(A)") "*   Purpose:  intgeration of two-dimensional &
       & function with two channels with equivalences"
    write (u, "(A)")

    write (u, "(A)") "* Initialise random number generator (default seed)"
    write (u, "(A)")

    allocate (rng_stream_t :: rng)
    call rng%init ()

    call rng%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise MC integrator with n_channel = 2 and n_dim = 1"
    write (u, "(A)")

    allocate (vamp2_test_func_3_t :: func)
    call func%init (n_dim = 1, n_channel = 2)
    config%equivalences = .true.
    mc_integrator = vamp2_t (n_channel = 2, n_dim = 1)
    call mc_integrator%set_config (config)
    call mc_integrator%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Initialise grid with n_calls = 20000 and set chains"
    write (u, "(A)")

    call mc_integrator%set_limits (x_lower, x_upper)
    call mc_integrator%set_calls (20000)

    write (u, "(A)")
    write (u, "(A)") "* Initialise equivalences"
    write (u, "(A)")

    eqv = vamp2_equivalences_t (n_eqv = 4, n_channel = 2, n_dim = 1)
    call eqv%set_equivalence &
         (i_eqv = 1, dest = 2, src = 1, perm = [1], mode = [VEQ_IDENTITY])
    call eqv%set_equivalence &
         (i_eqv = 2, dest = 1, src = 2, perm = [1], mode = [VEQ_IDENTITY])
    call eqv%set_equivalence &
         (i_eqv = 3, dest = 1, src = 1, perm = [1], mode = [VEQ_IDENTITY])
    call eqv%set_equivalence &
         (i_eqv = 4, dest = 2, src = 2, perm = [1], mode = [VEQ_IDENTITY])
    call eqv%write (u)
    call mc_integrator%set_equivalences (eqv)

    write (u, "(A)")
    write (u, "(A)") &
         "* Integrate with n_it = 3 and n_calls = 10000 (Grid-only Adaptation)"
    write (u, "(A)")

    call mc_integrator%integrate (func, rng, 3, &
         adapt_weights = .false., result=result, abserr=abserr)
    if (nearly_equal &
         (result, 2.000_default, rel_smallness = 0.003_default)) then
       write (u,  "(2x,A)") "Result: 2.000 [ok]"
    else
       write (u,  "(2x,A," // FMT_12 // ",A," // FMT_12 // ",A)") &
            "Result: ", result, " +/- ", abserr, " [not ok]"
    end if

    write (u, "(A)")
    write (u, "(A)") "* Compare the grids of both channels"
    write (u, "(A)")

    grid(1) = mc_integrator%get_grid(channel = 1)
    grid(2) = mc_integrator%get_grid(channel = 2)

    write (u, "(2X,A,1X,L1)") "Equal grids =", (grid(1) == grid(2))

    write (u, "(A)")
    write (u, "(A)") "* Write grid to file vamp2_5.grids"
    write (u, "(A)")

    unit = free_unit ()
    open (unit, file = "vamp2_5.grids", &
         action = "write", status = "replace")
    call mc_integrator%write_grids (unit)
    close (unit)

    write (u, "(A)")
    write (u, "(A)") "* Integrate with n_it = 3 and n_calls = 5000 (Precision)"
    write (u, "(A)")

    call mc_integrator%set_calls (5000)
    call mc_integrator%integrate (func, rng, 3, adapt_weights = .false., &
         refine_grids = .false., result=result, abserr=abserr)
    if (nearly_equal &
         (result, 2.000_default, rel_smallness = 0.003_default)) then
       write (u,  "(2x,A)") "Result: 2.000 [ok]"
    else
       write (u,  "(2x,A," // FMT_12 // ",A," // FMT_12 // ",A)") &
            "Result: ", result, " +/- ", abserr, " [not ok]"
    end if

    write (u, "(A)")
    write (u, "(A)") "* Cleanup"

    call mc_integrator%final ()
    call rng%final ()
    deallocate (rng)
  end subroutine vamp2_5

end module vamp2_uti
