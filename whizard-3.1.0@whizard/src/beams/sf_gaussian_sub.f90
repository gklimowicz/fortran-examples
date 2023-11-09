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

submodule (sf_gaussian) sf_gaussian_s

  use io_units
  use format_defs, only: FMT_12
  use file_registries
  use diagnostics
  use lorentz

  implicit none

contains

  module subroutine gaussian_data_init &
       (data, model, pdg_in, spread, rng_factory)
    class(gaussian_data_t), intent(out) :: data
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), dimension(2), intent(in) :: pdg_in
    real(default), dimension(2), intent(in) :: spread
    class(rng_factory_t), intent(inout), allocatable :: rng_factory
    if (any (spread < 0)) then
       call msg_fatal ("Gaussian beam spread: must not be negative")
    end if
    call data%flv_in(1)%init (pdg_in(1)%get (1), model)
    call data%flv_in(2)%init (pdg_in(2)%get (1), model)
    data%spread = spread
    call move_alloc (from = rng_factory, to = data%rng_factory)
  end subroutine gaussian_data_init

  module function gaussian_data_is_generator (data) result (flag)
    class(gaussian_data_t), intent(in) :: data
    logical :: flag
    flag = .true.
  end function gaussian_data_is_generator

  module function gaussian_data_get_n_par (data) result (n)
    class(gaussian_data_t), intent(in) :: data
    integer :: n
    n = 2
  end function gaussian_data_get_n_par

  module subroutine gaussian_data_get_pdg_out (data, pdg_out)
    class(gaussian_data_t), intent(in) :: data
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    integer :: i, n
    n = 2
    do i = 1, n
       pdg_out(i) = data%flv_in(i)%get_pdg ()
    end do
  end subroutine gaussian_data_get_pdg_out

  module subroutine gaussian_data_write (data, unit, verbose)
    class(gaussian_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)") "Gaussian beam spread data:"
    write (u, "(3x,A,A,A,A)") "prt_in = ", &
         char (data%flv_in(1)%get_name ()), &
         ", ", char (data%flv_in(2)%get_name ())
    write (u, "(3x,A,2(1x," // FMT_12 // "))") "spread =", data%spread
    call data%rng_factory%write (u)
  end subroutine gaussian_data_write

  module function gaussian_type_string (object) result (string)
    class(gaussian_t), intent(in) :: object
    type(string_t) :: string
    if (associated (object%data)) then
       string = "Gaussian: gaussian beam-energy spread"
    else
       string = "Gaussian: [undefined]"
    end if
  end function gaussian_type_string

  module subroutine gaussian_write (object, unit, testflag)
    class(gaussian_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%data)) then
       call object%data%write (u)
       call object%rng%write (u)
       call object%base_write (u, testflag)
    else
       write (u, "(1x,A)")  "gaussian data: [undefined]"
    end if
  end subroutine gaussian_write

  module subroutine gaussian_init (sf_int, data)
    class(gaussian_t), intent(out) :: sf_int
    class(sf_data_t), intent(in), target :: data
    real(default), dimension(2) :: m2
    real(default), dimension(0) :: mr2
    type(quantum_numbers_mask_t), dimension(4) :: mask
    integer, dimension(4) :: hel_lock
    type(quantum_numbers_t), dimension(4) :: qn_fc, qn_hel, qn
    type(polarization_t), target :: pol1, pol2
    type(polarization_iterator_t) :: it_hel1, it_hel2
    integer :: i
    select type (data)
    type is (gaussian_data_t)
       m2 = data%flv_in%get_mass () ** 2
       hel_lock = [3, 4, 1, 2]
       mask = quantum_numbers_mask (.false., .false., .false.)
       call sf_int%base_init (mask, m2, mr2, m2, hel_lock = hel_lock)
       sf_int%data => data
       do i = 1, 2
          call qn_fc(i)%init ( &
               flv = data%flv_in(i), &
               col = color_from_flavor (data%flv_in(i)))
          call qn_fc(i+2)%init ( &
               flv = data%flv_in(i), &
               col = color_from_flavor (data%flv_in(i)))
       end do
       call pol1%init_generic (data%flv_in(1))
       call it_hel1%init (pol1)
       do while (it_hel1%is_valid ())
          qn_hel(1) = it_hel1%get_quantum_numbers ()
          qn_hel(3) = it_hel1%get_quantum_numbers ()
          call pol2%init_generic (data%flv_in(2))
          call it_hel2%init (pol2)
          do while (it_hel2%is_valid ())
             qn_hel(2) = it_hel2%get_quantum_numbers ()
             qn_hel(4) = it_hel2%get_quantum_numbers ()
             qn = qn_hel .merge. qn_fc
             call sf_int%add_state (qn)
             call it_hel2%advance ()
          end do
          ! call pol2%final ()
          call it_hel1%advance ()
       end do
       ! call pol1%final ()
       call sf_int%freeze ()
       call sf_int%set_incoming ([1,2])
       call sf_int%set_outgoing ([3,4])
       sf_int%status = SF_INITIAL
    end select
    call sf_int%data%rng_factory%make (sf_int%rng)
  end subroutine gaussian_init

  module subroutine sf_gaussian_final (object)
    class(gaussian_t), intent(inout) :: object
    call object%interaction_t%final ()
  end subroutine sf_gaussian_final

  module function gaussian_is_generator (sf_int) result (flag)
    class(gaussian_t), intent(in) :: sf_int
    logical :: flag
    flag = sf_int%data%is_generator ()
  end function gaussian_is_generator

  module subroutine gaussian_generate_free (sf_int, r, rb, x_free)
    class(gaussian_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(inout) :: x_free
    real(default), dimension(size(r)) :: z
    associate (data => sf_int%data)
      do
         call sf_int%rng%generate_gaussian (z)
         rb = z * data%spread
         r = 1 - rb
         x_free = x_free * product (r)
         if (all (r > 0))  exit
      end do
    end associate
  end subroutine gaussian_generate_free

  module subroutine gaussian_complete_kinematics &
       (sf_int, x, xb, f, r, rb, map)
    class(gaussian_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: r
    real(default), dimension(:), intent(in) :: rb
    logical, intent(in) :: map
    if (map) then
       call msg_fatal ("gaussian: map flag not supported")
    else
       x = r
       xb= rb
       f = 1
    end if
    call sf_int%reduce_momenta (x)
  end subroutine gaussian_complete_kinematics

  module subroutine gaussian_inverse_kinematics &
       (sf_int, x, xb, f, r, rb, map, set_momenta)
    class(gaussian_t), intent(inout) :: sf_int
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
       call msg_fatal ("gaussian: map flag not supported")
    else
       r = x
       rb= xb
       f = 1
    end if
    if (set_mom) then
       call sf_int%reduce_momenta (x)
    end if
  end subroutine gaussian_inverse_kinematics

  module subroutine gaussian_apply &
       (sf_int, scale, negative_sf, rescale, i_sub)
    class(gaussian_t), intent(inout) :: sf_int
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(in), optional :: rescale
    integer, intent(in), optional :: i_sub
    real(default) :: f
    f = 1
    call sf_int%set_matrix_element (cmplx (f, kind=default))
    sf_int%status = SF_EVALUATED
  end subroutine gaussian_apply


end submodule sf_gaussian_s

