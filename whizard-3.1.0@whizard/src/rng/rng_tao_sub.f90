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

submodule (rng_tao) rng_tao_s

  use io_units
  use format_utils, only: write_indent

  implicit none

contains

  module subroutine rng_tao_write (rng, unit, indent)
    class(rng_tao_t), intent(in) :: rng
    integer, intent(in), optional :: unit, indent
    integer :: u, ind
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    call write_indent (u, ind)
    write (u, "(A)")  "TAO random-number generator:"
    call write_indent (u, ind)
    write (u, "(2x,A,I0)")  "seed  = ", rng%seed
    call write_indent (u, ind)
    write (u, "(2x,A,I0)")  "calls = ", rng%n_calls
  end subroutine rng_tao_write

  module subroutine rng_tao_init (rng, seed)
    class(rng_tao_t), intent(out) :: rng
    integer, intent(in), optional :: seed
    if (present (seed))  rng%seed = seed
    call tao_random_create (rng%state, rng%seed)
  end subroutine rng_tao_init

  module subroutine rng_tao_final (rng)
    class(rng_tao_t), intent(inout) :: rng
    call tao_random_destroy (rng%state)
  end subroutine rng_tao_final

  module subroutine rng_tao_generate_single (rng, x)
    class(rng_tao_t), intent(inout) :: rng
    real(default), intent(out) :: x
    real(default) :: r
    call tao_random_number (rng%state, r)
    x = r
    rng%n_calls = rng%n_calls + 1
  end subroutine rng_tao_generate_single

  module subroutine rng_tao_generate_array (rng, x)
    class(rng_tao_t), intent(inout) :: rng
    real(default), dimension(:), intent(out) :: x
    real(default) :: r
    integer :: i
    do i = 1, size (x)
       call tao_random_number (rng%state, r)
       x(i) = r
    end do
    rng%n_calls = rng%n_calls + size (x)
  end subroutine rng_tao_generate_array

  module subroutine rng_tao_factory_write (object, unit)
    class(rng_tao_factory_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A,2(I0,A))") &
         "RNG factory: tao (", object%s, ",", object%i, ")"
  end subroutine rng_tao_factory_write

  module subroutine rng_tao_factory_init (factory, seed)
    class(rng_tao_factory_t), intent(out) :: factory
    integer(i16), intent(in), optional :: seed
    if (present (seed))  factory%s = seed
  end subroutine rng_tao_factory_init


end submodule rng_tao_s

