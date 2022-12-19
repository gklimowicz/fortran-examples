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

submodule (sf_isr) sf_isr_s

  use io_units
  use constants, only: pi
  use format_defs, only: FMT_15, FMT_19
  use numeric_utils
  use diagnostics
  use physics_defs, only: PHOTON
  use sm_physics, only: Li2
  use lorentz
  use colors
  use quantum_numbers
  use polarizations

  implicit none

contains

  module subroutine isr_data_init (data, model, pdg_in, alpha, q_max, &
       mass, order, recoil, keep_energy)
    class(isr_data_t), intent(out) :: data
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), intent(in) :: pdg_in
    real(default), intent(in) :: alpha
    real(default), intent(in) :: q_max
    real(default), intent(in), optional :: mass
    integer, intent(in), optional :: order
    logical, intent(in), optional :: recoil
    logical, intent(in), optional :: keep_energy
    integer :: i, n_flv
    real(default) :: charge
    data%model => model
    n_flv = pdg_in%get_length ()
    allocate (data%flv_in (n_flv))
    do i = 1, n_flv
       call data%flv_in(i)%init (pdg_in%get (i), model)
    end do
    data%alpha = alpha
    data%q_max = q_max
    if (present (order)) then
       call data%set_order (order)
    end if
    if (present (recoil)) then
       data%recoil = recoil
    end if
    if (present (keep_energy)) then
       data%keep_energy = keep_energy
    end if
    data%real_mass = data%flv_in(1)%get_mass ()
    if (present (mass)) then
       if (mass > 0) then
          data%mass = mass
       else
          data%mass = data%real_mass
          if (any (data%flv_in%get_mass () /= data%mass)) then
             data%error = MASS_MIX;  return
          end if
       end if
    else
       data%mass = data%real_mass
       if (any (data%flv_in%get_mass () /= data%mass)) then
          data%error = MASS_MIX;  return
       end if
    end if
    if (vanishes (data%mass)) then
       data%error = ZERO_MASS;  return
    else if (data%mass >= data%q_max) then
       data%error = Q_MAX_TOO_SMALL;  return
    end if
    data%log = log (1 + (data%q_max / data%mass)**2)
    charge = data%flv_in(1)%get_charge ()
    if (any (abs (data%flv_in%get_charge ()) /= abs (charge))) then
       data%error = CHARGE_MIX;  return
    else if (charge == 0) then
       data%error = CHARGE_ZERO;  return
    end if
    data%eps = data%alpha / pi * charge ** 2 &
         * (2 * log (data%q_max / data%mass) - 1)
    if (data%eps > 1) then
       data%error = EPS_TOO_LARGE;  return
    end if
    call data%pdf%init (data%mass, data%alpha, charge, data%q_max, data%order, &
         0, 1)
  end subroutine isr_data_init

  elemental module subroutine isr_data_set_order (data, order)
    class(isr_data_t), intent(inout) :: data
    integer, intent(in) :: order
    if (order < 0 .or. order > 3) then
       data%error = INVALID_ORDER
    else
       data%order = order
    end if
  end subroutine isr_data_set_order

  module subroutine isr_data_check (data)
    class(isr_data_t), intent(in) :: data
    select case (data%error)
    case (ZERO_MASS)
       call msg_fatal ("ISR: Particle mass is zero")
    case (Q_MAX_TOO_SMALL)
       call msg_fatal ("ISR: Particle mass exceeds Qmax")
    case (EPS_TOO_LARGE)
       call msg_fatal ("ISR: Expansion parameter too large, " // &
            "perturbative expansion breaks down")
    case (INVALID_ORDER)
       call msg_error ("ISR: LLA order invalid (valid values are 0,1,2,3)")
    case (MASS_MIX)
       call msg_fatal ("ISR: Incoming particle masses must be uniform")
    case (CHARGE_MIX)
       call msg_fatal ("ISR: Incoming particle charges must be uniform")
    case (CHARGE_ZERO)
       call msg_fatal ("ISR: Incoming particle must be charged")
    end select
  end subroutine isr_data_check

  module subroutine isr_data_write (data, unit, verbose)
    class(isr_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)") "ISR data:"
    if (allocated (data%flv_in)) then
       write (u, "(3x,A)", advance="no") "  flavor =  "
       do i = 1, size (data%flv_in)
          if (i > 1)  write (u, "(',',1x)", advance="no")
          call data%flv_in(i)%write (u)
       end do
       write (u, *)
       write (u, "(3x,A," // FMT_19 // ")") "  alpha    = ", data%alpha
       write (u, "(3x,A," // FMT_19 // ")") "  q_max    = ", data%q_max
       write (u, "(3x,A," // FMT_19 // ")") "  mass     = ", data%mass
       write (u, "(3x,A," // FMT_19 // ")") "  eps      = ", data%eps
       write (u, "(3x,A," // FMT_19 // ")") "  log      = ", data%log
       write (u, "(3x,A,I2)")      "  order    = ", data%order
       write (u, "(3x,A,L2)")      "  recoil   = ", data%recoil
       write (u, "(3x,A,L2)")      "  keep en. = ", data%keep_energy
    else
       write (u, "(3x,A)") "[undefined]"
    end if
  end subroutine isr_data_write

  module function isr_data_get_n_par (data) result (n)
    class(isr_data_t), intent(in) :: data
    integer :: n
    if (data%recoil) then
       n = 3
    else
       n = 1
    end if
  end function isr_data_get_n_par

  module subroutine isr_data_get_pdg_out (data, pdg_out)
    class(isr_data_t), intent(in) :: data
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    pdg_out(1) = data%flv_in%get_pdg ()
  end subroutine isr_data_get_pdg_out

  module function isr_data_get_eps (data) result (eps)
    class(isr_data_t), intent(in) :: data
    real(default) :: eps
    eps = data%eps
  end function isr_data_get_eps

  module function isr_type_string (object) result (string)
    class(isr_t), intent(in) :: object
    type(string_t) :: string
    if (associated (object%data)) then
       string = "ISR: e+ e- ISR spectrum"
    else
       string = "ISR: [undefined]"
    end if
  end function isr_type_string

  module subroutine isr_write (object, unit, testflag)
    class(isr_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%data)) then
       call object%data%write (u)
       if (object%status >= SF_DONE_KINEMATICS) then
          write (u, "(1x,A)")  "SF parameters:"
          write (u, "(3x,A," // FMT_15 // ")")  "x =", object%x
          write (u, "(3x,A," // FMT_15 // ")")  "xb=", object%xb
       end if
       call object%base_write (u, testflag)
    else
       write (u, "(1x,A)")  "ISR data: [undefined]"
    end if
  end subroutine isr_write

  module subroutine isr_set_order (object, order)
    class(isr_t), intent(inout) :: object
    integer, intent(in) :: order
    call object%data%set_order (order)
    call object%data%pdf%set_order (order)
  end subroutine isr_set_order

  module subroutine isr_complete_kinematics (sf_int, x, xb, f, r, rb, map)
    class(isr_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: r
    real(default), dimension(:), intent(in) :: rb
    logical, intent(in) :: map
    real(default) :: eps
    eps = sf_int%data%eps
    if (map) then
       call map_power_1 (sf_int%xb, f, rb(1), eps)
    else
       sf_int%xb = rb(1)
       if (rb(1) > 0) then
          f = 1
       else
          f = 0
       end if
    end if
    sf_int%x = 1 - sf_int%xb
    x(1) = sf_int%x
    xb(1) = sf_int%xb
    if (size (x) == 3) then
       x(2:3) = r(2:3)
       xb(2:3) = rb(2:3)
    end if
    call sf_int%split_momentum (x, xb)
    select case (sf_int%status)
    case (SF_FAILED_KINEMATICS)
       sf_int%x = 0
       sf_int%xb= 0
       f = 0
    end select
  end subroutine isr_complete_kinematics

  module subroutine sf_isr_recover_x (sf_int, x, xb, x_free)
    class(isr_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(inout), optional :: x_free
    call sf_int%base_recover_x (x, xb, x_free)
    sf_int%x  = x(1)
    sf_int%xb = xb(1)
  end subroutine sf_isr_recover_x

  module subroutine isr_inverse_kinematics &
       (sf_int, x, xb, f, r, rb, map, set_momenta)
    class(isr_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: r
    real(default), dimension(:), intent(out) :: rb
    logical, intent(in) :: map
    logical, intent(in), optional :: set_momenta
    real(default) :: eps
    logical :: set_mom
    set_mom = .false.;  if (present (set_momenta))  set_mom = set_momenta
    eps = sf_int%data%eps
    if (map) then
       call map_power_inverse_1 (xb(1), f, rb(1), eps)
    else
       rb(1) = xb(1)
       if (rb(1) > 0) then
          f = 1
       else
          f = 0
       end if
    end if
    r(1) = 1 - rb(1)
    if (size(r) == 3) then
       r(2:3) = x(2:3)
       rb(2:3)= xb(2:3)
    end if
    if (set_mom) then
       call sf_int%split_momentum (x, xb)
       select case (sf_int%status)
       case (SF_FAILED_KINEMATICS)
          r = 0
          rb= 0
          f = 0
       end select
    end if
  end subroutine isr_inverse_kinematics

  module subroutine isr_init (sf_int, data)
    class(isr_t), intent(out) :: sf_int
    class(sf_data_t), intent(in), target :: data
    type(quantum_numbers_mask_t), dimension(3) :: mask
    integer, dimension(3) :: hel_lock
    type(polarization_t), target :: pol
    type(quantum_numbers_t), dimension(1) :: qn_fc
    type(flavor_t) :: flv_photon
    type(color_t) :: col_photon
    type(quantum_numbers_t) :: qn_hel, qn_photon, qn
    type(polarization_iterator_t) :: it_hel
    real(default) :: m2
    integer :: i
    mask = quantum_numbers_mask (.false., .false., &
         mask_h = [.false., .true., .false.])
    hel_lock = [3, 0, 1]
    select type (data)
    type is (isr_data_t)
       m2 = data%mass**2
       call sf_int%base_init (mask, [m2], [0._default], [m2], &
            hel_lock = hel_lock)
       sf_int%data => data
       call flv_photon%init (PHOTON, data%model)
       call col_photon%init ()
       call qn_photon%init (flv_photon, col_photon)
       call qn_photon%tag_radiated ()
       do i = 1, size (data%flv_in)
          call pol%init_generic (data%flv_in(i))
          call qn_fc(1)%init (&
               flv = data%flv_in(i), &
               col = color_from_flavor (data%flv_in(i), 1))
          call it_hel%init (pol)
          do while (it_hel%is_valid ())
             qn_hel = it_hel%get_quantum_numbers ()
             qn = qn_hel .merge. qn_fc(1)
             call sf_int%add_state ([qn, qn_photon, qn])
             call it_hel%advance ()
          end do
          ! call pol%final ()  !!! Obsolete
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
       sf_int%status = SF_INITIAL
    end select
  end subroutine isr_init

  module subroutine isr_apply (sf_int, scale, negative_sf, rescale, i_sub)
    class(isr_t), intent(inout) :: sf_int
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(in), optional :: rescale
    integer, intent(in), optional :: i_sub
    real(default) :: f, finv, x, xb, eps, rb
    real(default) :: log_x, log_xb, x_2
    associate (data => sf_int%data)
      eps = sf_int%data%eps
      x = sf_int%x
      xb = sf_int%xb
      call map_power_inverse_1 (xb, finv, rb, eps)
      if (finv > 0) then
         f = 1 / finv
      else
         f = 0
      end if
      call data%pdf%evolve_qed_pdf (x, xb, rb, f)
    end associate
    call sf_int%set_matrix_element (cmplx (f, kind=default))
    sf_int%status = SF_EVALUATED
  end subroutine isr_apply


end submodule sf_isr_s

