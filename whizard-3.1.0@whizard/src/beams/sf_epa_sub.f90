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

submodule (sf_epa) sf_epa_s

  use io_units
  use constants, only: pi
  use format_defs, only: FMT_17, FMT_19
  use numeric_utils
  use diagnostics
  use physics_defs, only: PHOTON
  use colors

  implicit none

contains

  module subroutine epa_data_init (data, model, mode, pdg_in, alpha, &
       x_min, q_min, q_max, mass, recoil, keep_energy)
    class(epa_data_t), intent(inout) :: data
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), intent(in) :: pdg_in
    integer, intent(in) :: mode
    real(default), intent(in) :: alpha, x_min, q_min, q_max
    real(default), intent(in), optional :: mass
    logical, intent(in), optional :: recoil
    logical, intent(in), optional :: keep_energy
    integer :: n_flv, i
    data%model => model
    data%mode = mode
    n_flv = pdg_in%get_length ()
    allocate (data%flv_in (n_flv))
    do i = 1, n_flv
       call data%flv_in(i)%init (pdg_in%get (i), model)
    end do
    data%alpha = alpha
    data%E_max = q_max / 2
    data%x_min = x_min
    data%x_max = 1
    if (vanishes (data%x_min)) then
       data%error = ZERO_XMIN;  return
    end if
    data%q_min = q_min
    data%q_max = q_max
    select case (char (data%model%get_name ()))
    case ("QCD","Test")
       data%error = NO_EPA;  return
    end select
    if (present (recoil)) then
       data%recoil = recoil
    end if
    if (present (keep_energy)) then
       data%keep_energy = keep_energy
    end if
    if (present (mass)) then
       data%mass = mass
    else
       data%mass = data%flv_in(1)%get_mass ()
       if (any (data%flv_in%get_mass () /= data%mass)) then
          data%error = MASS_MIX;  return
       end if
    end if
    if (max (data%mass, data%q_min) == 0) then
       data%error = ZERO_QMIN;  return
    else if (max (data%mass, data%q_min) >= data%E_max) then
       data%error = Q_MAX_TOO_SMALL;  return
    end if
    data%log = log ((data%q_max / max (data%mass, data%q_min)) ** 2 )
    data%a  = data%alpha / pi
    data%c0 = log (data%x_min) * (data%log - log (data%x_min))
    data%c1 = log (data%x_max) * (data%log - log (data%x_max))
    data%dc = data%c1 - data%c0
  end subroutine epa_data_init

  module subroutine epa_data_check (data)
    class(epa_data_t), intent(in) :: data
    select case (data%error)
    case (NO_EPA)
       call msg_fatal ("EPA structure function not available for model " &
            // char (data%model%get_name ()) // ".")
    case (ZERO_QMIN)
       call msg_fatal ("EPA: Particle mass is zero")
    case (Q_MAX_TOO_SMALL)
       call msg_fatal ("EPA: Particle mass exceeds Qmax")
    case (ZERO_XMIN)
       call msg_fatal ("EPA: x_min must be larger than zero")
    case (MASS_MIX)
       call msg_fatal ("EPA: incoming particle masses must be uniform")
    end select
  end subroutine epa_data_check

  module subroutine epa_data_write (data, unit, verbose)
    class(epa_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)") "EPA data:"
    if (allocated (data%flv_in)) then
       write (u, "(3x,A)", advance="no") "  flavor =  "
       do i = 1, size (data%flv_in)
          if (i > 1)  write (u, "(',',1x)", advance="no")
          call data%flv_in(i)%write (u)
       end do
       write (u, *)
       write (u, "(3x,A," // FMT_19 // ")") "  alpha    = ", data%alpha
       write (u, "(3x,A," // FMT_19 // ")") "  x_min    = ", data%x_min
       write (u, "(3x,A," // FMT_19 // ")") "  x_max    = ", data%x_max
       write (u, "(3x,A," // FMT_19 // ")") "  q_min    = ", data%q_min
       write (u, "(3x,A," // FMT_19 // ")") "  q_max    = ", data%q_max
       write (u, "(3x,A," // FMT_19 // ")") "  E_max    = ", data%e_max
       write (u, "(3x,A," // FMT_19 // ")") "  mass     = ", data%mass
       write (u, "(3x,A," // FMT_19 // ")") "  a        = ", data%a
       write (u, "(3x,A," // FMT_19 // ")") "  c0       = ", data%c0
       write (u, "(3x,A," // FMT_19 // ")") "  c1       = ", data%c1
       write (u, "(3x,A," // FMT_19 // ")") "  log      = ", data%log
       write (u, "(3x,A,L2)")      "  recoil   = ", data%recoil
       write (u, "(3x,A,L2)")      "  keep en. = ", data%keep_energy
    else
       write (u, "(3x,A)") "[undefined]"
    end if
  end subroutine epa_data_write

  module function epa_data_get_n_par (data) result (n)
    class(epa_data_t), intent(in) :: data
    integer :: n
    if (data%recoil) then
       n = 3
    else
       n = 1
    end if
  end function epa_data_get_n_par

  module subroutine epa_data_get_pdg_out (data, pdg_out)
    class(epa_data_t), intent(in) :: data
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    pdg_out(1) = PHOTON
  end subroutine epa_data_get_pdg_out

  module function epa_type_string (object) result (string)
    class(epa_t), intent(in) :: object
    type(string_t) :: string
    if (associated (object%data)) then
       string = "EPA: equivalent photon approx."
    else
       string = "EPA: [undefined]"
    end if
  end function epa_type_string

  module subroutine epa_write (object, unit, testflag)
    class(epa_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%data)) then
       call object%data%write (u)
       if (object%status >= SF_DONE_KINEMATICS) then
          write (u, "(1x,A)")  "SF parameters:"
          write (u, "(3x,A," // FMT_17 // ")")  "x =", object%x
          if (object%status >= SF_FAILED_EVALUATION) then
             write (u, "(3x,A," // FMT_17 // ")")  "E =", object%E
          end if
       end if
       call object%base_write (u, testflag)
    else
       write (u, "(1x,A)")  "EPA data: [undefined]"
    end if
  end subroutine epa_write

  module subroutine epa_init (sf_int, data)
    class(epa_t), intent(out) :: sf_int
    class(sf_data_t), intent(in), target :: data
    type(quantum_numbers_mask_t), dimension(3) :: mask
    integer, dimension(3) :: hel_lock
    type(polarization_t), target :: pol
    type(quantum_numbers_t), dimension(1) :: qn_fc
    type(flavor_t) :: flv_photon
    type(color_t) :: col_photon
    type(quantum_numbers_t) :: qn_hel, qn_photon, qn, qn_rad
    type(polarization_iterator_t) :: it_hel
    integer :: i
    mask = quantum_numbers_mask (.false., .false., &
         mask_h = [.false., .false., .true.])
    hel_lock = [2, 1, 0]
    select type (data)
    type is (epa_data_t)
       call sf_int%base_init (mask, [data%mass**2], &
            [data%mass**2], [0._default], hel_lock = hel_lock)
       sf_int%data => data
       call flv_photon%init (PHOTON, data%model)
       call col_photon%init ()
       call qn_photon%init (flv_photon, col_photon)
       do i = 1, size (data%flv_in)
          call pol%init_generic (data%flv_in(i))
          call qn_fc(1)%init ( &
               flv = data%flv_in(i), &
               col = color_from_flavor (data%flv_in(i), 1))
          call it_hel%init (pol)
          do while (it_hel%is_valid ())
             qn_hel = it_hel%get_quantum_numbers ()
             qn = qn_hel .merge. qn_fc(1)
             qn_rad = qn
             call qn_rad%tag_radiated ()
             call sf_int%add_state ([qn, qn_rad, qn_photon])
             call it_hel%advance ()
          end do
          !  call pol%final ()
       end do
       call sf_int%freeze ()
       if (data%keep_energy) then
          sf_int%on_shell_mode = KEEP_ENERGY
       else
          sf_int%on_shell_mode = KEEP_MOMENTUM
       end if
       call sf_int%set_incoming ([1])
       call sf_int%set_radiated ([2])
       call sf_int%set_outgoing ([3])
    end select
  end subroutine epa_init

  module subroutine epa_setup_constants (sf_int)
    class(epa_t), intent(inout), target :: sf_int
    type(state_iterator_t) :: it
    type(flavor_t) :: flv
    integer :: i, n_me
    n_me = sf_int%get_n_matrix_elements ()
    allocate (sf_int%charge2 (n_me))
    call it%init (sf_int%interaction_t%get_state_matrix_ptr ())
    do while (it%is_valid ())
       i = it%get_me_index ()
       flv = it%get_flavor (1)
       sf_int%charge2(i) = flv%get_charge () ** 2
       call it%advance ()
    end do
    sf_int%status = SF_INITIAL
  end subroutine epa_setup_constants

  module subroutine epa_complete_kinematics (sf_int, x, xb, f, r, rb, map)
    class(epa_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: r
    real(default), dimension(:), intent(in) :: rb
    logical, intent(in) :: map
    real(default) :: delta, sqrt_delta, lx
    if (map) then
       associate (data => sf_int%data)
         delta = data%log ** 2 -  4 * (r(1) * data%c1 + rb(1) * data%c0)
         if (delta > 0) then
            sqrt_delta = sqrt (delta)
            lx = (data%log - sqrt_delta) / 2
         else
            sf_int%status = SF_FAILED_KINEMATICS
            f = 0
            return
         end if
         x(1) = exp (lx)
         f = x(1) * data%dc / sqrt_delta
       end associate
    else
       x(1) = r(1)
       if (sf_int%data%x_min < x(1) .and. x(1) < sf_int%data%x_max) then
          f = 1
       else
          sf_int%status = SF_FAILED_KINEMATICS
          f = 0
          return
       end if
    end if
    xb(1) = 1 - x(1)
    if (size(x) == 3) then
       x(2:3) = r(2:3)
       xb(2:3) = rb(2:3)
    end if
    call sf_int%split_momentum (x, xb)
    select case (sf_int%status)
    case (SF_DONE_KINEMATICS)
       sf_int%x = x(1)
       sf_int%xb= xb(1)
       sf_int%E  = energy (sf_int%get_momentum (1))
    case (SF_FAILED_KINEMATICS)
       sf_int%x = 0
       sf_int%xb= 0
       f = 0
    end select
  end subroutine epa_complete_kinematics

  module subroutine sf_epa_recover_x (sf_int, x, xb, x_free)
    class(epa_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(inout), optional :: x_free
    call sf_int%base_recover_x (x, xb, x_free)
    sf_int%x  = x(1)
    sf_int%xb = xb(1)
  end subroutine sf_epa_recover_x

  module subroutine epa_inverse_kinematics &
       (sf_int, x, xb, f, r, rb, map, set_momenta)
    class(epa_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: r
    real(default), dimension(:), intent(out) :: rb
    logical, intent(in) :: map
    logical, intent(in), optional :: set_momenta
    real(default) :: lx, delta, sqrt_delta, c
    logical :: set_mom
    set_mom = .false.;  if (present (set_momenta))  set_mom = set_momenta
    if (map) then
       associate (data => sf_int%data)
         lx = log (x(1))
         sqrt_delta = data%log - 2 * lx
         delta = sqrt_delta ** 2
         c = (data%log ** 2 - delta) / 4
         r (1) = (c - data%c0) / data%dc
         rb(1) = (data%c1 - c) / data%dc
         f = x(1) * data%dc / sqrt_delta
       end associate
    else
       r (1) = x(1)
       rb(1) = xb(1)
       if (sf_int%data%x_min < x(1) .and. x(1) < sf_int%data%x_max) then
          f = 1
       else
          f = 0
       end if
    end if
    if (size(r) == 3) then
       r (2:3) = x(2:3)
       rb(2:3) = xb(2:3)
    end if
    if (set_mom) then
       call sf_int%split_momentum (x, xb)
       select case (sf_int%status)
       case (SF_FAILED_KINEMATICS);  f = 0
       end select
    end if
    sf_int%E  = energy (sf_int%get_momentum (1))
  end subroutine epa_inverse_kinematics

  module subroutine epa_apply (sf_int, scale, negative_sf, rescale, i_sub)
    class(epa_t), intent(inout) :: sf_int
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(in), optional :: rescale
    integer, intent(in), optional :: i_sub
    real(default) :: x, xb, qminsq, qmaxsq, f, E, m2
    associate (data => sf_int%data)
      x = sf_int%x
      xb= sf_int%xb
      E = sf_int%E
      m2 = data%mass ** 2
      qminsq = max (x ** 2 / xb * data%mass ** 2, data%q_min ** 2)
      select case (data%mode)
      case (0)
         qmaxsq = min (4 * xb * E ** 2, data%q_max ** 2)
         if (qminsq < qmaxsq) then
            f = data%a / x &
                 * ((xb + x ** 2 / 2) * log (qmaxsq / qminsq) &
                 - (1 - x / 2) ** 2 &
                 * log ((x**2 + qmaxsq / E ** 2) / (x**2 + qminsq / E ** 2)) &
                 - x ** 2 * data%mass ** 2 / qminsq * (1 - qminsq / qmaxsq))
         else
            f = 0
         end if
      case (1)
         qmaxsq = min (4 * xb * E ** 2, data%q_max ** 2)
         if (qminsq < qmaxsq) then
            f = data%a / x &
                 * ((xb + x ** 2 / 2) * log (qmaxsq / qminsq) &
                 - x ** 2 * data%mass ** 2 / qminsq * (1 - qminsq / qmaxsq))
         else
            f = 0
         end if
      case (2)
         qmaxsq = data%q_max ** 2
         if (data%mass ** 2 < qmaxsq) then
            f = data%a / x &
                 * ((xb + x ** 2 / 2) * log (qmaxsq / m2) &
                 - x ** 2 * data%mass ** 2 / qminsq * (1 - qminsq / qmaxsq))
         else
            f = 0
         end if
      case (3)
         qmaxsq = data%q_max ** 2
         if (data%mass ** 2 < qmaxsq) then
            f = data%a / x &
                 * ((xb + x ** 2 / 2) * log (qmaxsq / m2) &
                 - x ** 2 * (1 - m2 / qmaxsq))
         else
            f = 0
         end if
      case (4)
         qmaxsq = data%q_max ** 2
         if (data%mass ** 2 < qmaxsq) then
            f = data%a / x &
                 * ((xb + x ** 2 / 2) * log (qmaxsq / m2))
         else
            f = 0
         end if
      end select
      if (sf_int%get_n_matrix_elements () > 1) then
         f = f / 2
      end if
      call sf_int%set_matrix_element &
           (cmplx (f, kind=default) * sf_int%charge2)
    end associate
    sf_int%status = SF_EVALUATED
  end subroutine epa_apply


end submodule sf_epa_s

