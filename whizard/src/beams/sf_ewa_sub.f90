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

submodule (sf_ewa) sf_ewa_s

  use io_units
  use constants, only: pi
  use format_defs, only: FMT_17, FMT_19
  use numeric_utils
  use diagnostics
  use physics_defs, only: W_BOSON, Z_BOSON
  use lorentz
  use colors

  implicit none

contains

  module subroutine ewa_data_init (data, model, pdg_in, x_min, pt_max, &
        sqrts, recoil, keep_energy, mass)
    class(ewa_data_t), intent(inout) :: data
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), intent(in) :: pdg_in
    real(default), intent(in) :: x_min, pt_max, sqrts
    logical, intent(in) :: recoil, keep_energy
    real(default), intent(in), optional :: mass
    real(default) :: g, ee
    integer :: n_flv, i
    data%model => model
    if (.not. any (pdg_in .match. &
         [1,2,3,4,5,6,11,13,15,-1,-2,-3,-4,-5,-6,-11,-13,-15])) then
       data%error = WRONG_PRT;  return
    end if
    n_flv = pdg_in%get_length ()
    allocate (data%flv_in (n_flv))
    allocate (data%flv_out(n_flv))
    do i = 1, n_flv
       call data%flv_in(i)%init (pdg_in%get (i), model)
    end do
    data%pt_max = pt_max
    data%sqrts = sqrts
    data%x_min = x_min
    data%x_max = 1
    if (vanishes (data%x_min)) then
       data%error = ZERO_XMIN;  return
    end if
    select case (char (data%model%get_name ()))
    case ("QCD","QED","Test")
       data%error = NO_EWA;  return
    end select
    ee = data%model%get_real (var_str ("ee"))
    data%sinthw = data%model%get_real (var_str ("sw"))
    data%costhw = data%model%get_real (var_str ("cw"))
    data%mZ = data%model%get_real (var_str ("mZ"))
    data%mW = data%model%get_real (var_str ("mW"))
    if (data%sinthw /= 0) then
       g = ee / data%sinthw
    else
       data%error = ZERO_SW;  return
    end if
    data%cv = g / 2._default
    data%ca = g / 2._default
    data%coeff = 1._default / (8._default * PI**2)
    data%recoil = recoil
    data%keep_energy = keep_energy
    if (present (mass)) then
       data%mass = mass
       data%m_out = mass
       data%mass_set = .true.
    else
       data%mass = data%flv_in(1)%get_mass ()
       if (any (data%flv_in%get_mass () /= data%mass)) then
          data%error = MASS_MIX;  return
       end if
    end if
  end subroutine ewa_data_init

  module subroutine ewa_set_id (data, id)
    class(ewa_data_t), intent(inout) :: data
    integer, intent(in) :: id
    integer :: i, isospin, pdg
    if (.not. allocated (data%flv_in)) &
         call msg_bug ("EWA: incoming particles not set")
    data%id = id
    select case (data%id)
    case (23)
       data%m_out = data%mass
       data%flv_out = data%flv_in
    case (24)
       do i = 1, size (data%flv_in)
          pdg = data%flv_in(i)%get_pdg ()
          isospin = data%flv_in(i)%get_isospin_type ()
          if (isospin > 0) then
             !!! up-type quark or neutrinos
             if (data%flv_in(i)%is_antiparticle ()) then
                call data%flv_out(i)%init (pdg + 1, data%model)
             else
                call data%flv_out(i)%init (pdg - 1, data%model)
             end if
          else
             !!! down-type quark or lepton
             if (data%flv_in(i)%is_antiparticle ()) then
                call data%flv_out(i)%init (pdg - 1, data%model)
             else
                call data%flv_out(i)%init (pdg + 1, data%model)
             end if
          end if
       end do
       if (.not. data%mass_set) then
          data%m_out = data%flv_out(1)%get_mass ()
          if (any (data%flv_out%get_mass () /= data%m_out)) then
             data%error = MASS_MIX_OUT;  return
          end if
       end if
    end select
  end subroutine ewa_set_id

  module subroutine ewa_data_check (data)
    class(ewa_data_t), intent(in) :: data
    select case (data%error)
    case (WRONG_PRT)
       call msg_fatal ("EWA structure function only accessible for " &
            // "SM quarks and leptons.")
    case (NO_EWA)
       call msg_fatal ("EWA structure function not available for model " &
            // char (data%model%get_name ()))
    case (ZERO_SW)
       call msg_fatal ("EWA: Vanishing value of sin(theta_w)")
    case (ZERO_QMIN)
       call msg_fatal ("EWA: Particle mass is zero")
    case (Q_MAX_TOO_SMALL)
       call msg_fatal ("EWA: Particle mass exceeds Qmax")
    case (ZERO_XMIN)
       call msg_fatal ("EWA: x_min must be larger than zero")
    case (MASS_MIX)
       call msg_fatal ("EWA: incoming particle masses must be uniform")
    case (MASS_MIX_OUT)
       call msg_fatal ("EWA: outgoing particle masses must be uniform")
    case (ISOSPIN_MIX)
       call msg_fatal ("EWA: incoming particle isospins must be uniform")
    end select
  end subroutine ewa_data_check

  module subroutine ewa_data_write (data, unit, verbose)
    class(ewa_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)") "EWA data:"
    if (allocated (data%flv_in) .and. allocated (data%flv_out)) then
       write (u, "(3x,A)", advance="no") "  flavor(in)  =  "
       do i = 1, size (data%flv_in)
          if (i > 1)  write (u, "(',',1x)", advance="no")
          call data%flv_in(i)%write (u)
       end do
       write (u, *)
       write (u, "(3x,A)", advance="no") "  flavor(out) =  "
       do i = 1, size (data%flv_out)
          if (i > 1)  write (u, "(',',1x)", advance="no")
          call data%flv_out(i)%write (u)
       end do
       write (u, *)
       write (u, "(3x,A," // FMT_19 // ")") "  x_min     = ", data%x_min
       write (u, "(3x,A," // FMT_19 // ")") "  x_max     = ", data%x_max
       write (u, "(3x,A," // FMT_19 // ")") "  pt_max    = ", data%pt_max
       write (u, "(3x,A," // FMT_19 // ")") "  sqrts     = ", data%sqrts
       write (u, "(3x,A," // FMT_19 // ")") "  mass      = ", data%mass
       write (u, "(3x,A," // FMT_19 // ")") "  cv        = ", data%cv
       write (u, "(3x,A," // FMT_19 // ")") "  ca        = ", data%ca
       write (u, "(3x,A," // FMT_19 // ")") "  coeff     = ", data%coeff
       write (u, "(3x,A," // FMT_19 // ")") "  costhw    = ", data%costhw
       write (u, "(3x,A," // FMT_19 // ")") "  sinthw    = ", data%sinthw
       write (u, "(3x,A," // FMT_19 // ")") "  mZ        = ", data%mZ
       write (u, "(3x,A," // FMT_19 // ")") "  mW        = ", data%mW
       write (u, "(3x,A,L2)")      "  recoil    = ", data%recoil
       write (u, "(3x,A,L2)")      "  keep en.  = ", data%keep_energy
       write (u, "(3x,A,I2)")      "  PDG (VB)  = ", data%id
    else
       write (u, "(3x,A)") "[undefined]"
    end if
  end subroutine ewa_data_write

  module function ewa_data_get_n_par (data) result (n)
    class(ewa_data_t), intent(in) :: data
    integer :: n
    if (data%recoil) then
       n = 3
    else
       n = 1
    end if
  end function ewa_data_get_n_par

  module subroutine ewa_data_get_pdg_out (data, pdg_out)
    class(ewa_data_t), intent(in) :: data
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    integer, dimension(:), allocatable :: pdg1
    integer :: i, n_flv
    if (allocated (data%flv_out)) then
       n_flv = size (data%flv_out)
    else
       n_flv = 0
    end if
    allocate (pdg1 (n_flv))
    do i = 1, n_flv
       pdg1(i) = data%flv_out(i)%get_pdg ()
    end do
    pdg_out(1) = pdg1
  end subroutine ewa_data_get_pdg_out

  module function ewa_type_string (object) result (string)
    class(ewa_t), intent(in) :: object
    type(string_t) :: string
    if (associated (object%data)) then
       string = "EWA: equivalent W/Z approx."
    else
       string = "EWA: [undefined]"
    end if
  end function ewa_type_string

  module subroutine ewa_write (object, unit, testflag)
    class(ewa_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%data)) then
       call object%data%write (u)
       if (object%status >= SF_DONE_KINEMATICS) then
          write (u, "(1x,A)")  "SF parameters:"
          write (u, "(3x,A," // FMT_17 // ")")  "x =", object%x
          write (u, "(3x,A," // FMT_17 // ")")  "xb=", object%xb
       end if
       call object%base_write (u, testflag)
    else
       write (u, "(1x,A)")  "EWA data: [undefined]"
    end if
  end subroutine ewa_write

  module subroutine ewa_init (sf_int, data)
    class(ewa_t), intent(out) :: sf_int
    class(sf_data_t), intent(in), target :: data
    type(quantum_numbers_mask_t), dimension(3) :: mask
    integer, dimension(3) :: hel_lock
    type(polarization_t), target :: pol
    type(quantum_numbers_t), dimension(1) :: qn_fc, qn_fc_fin
    type(flavor_t) :: flv_z, flv_wp, flv_wm
    type(color_t) :: col0
    type(quantum_numbers_t) :: qn_hel, qn_z, qn_wp, qn_wm, qn, qn_rad, qn_w
    type(polarization_iterator_t) :: it_hel
    integer :: i, isospin
    select type (data)
    type is (ewa_data_t)
       mask = quantum_numbers_mask (.false., .false., &
            mask_h = [.false., .false., .true.])
       hel_lock = [2, 1, 0]
       call col0%init ()
       select case (data%id)
       case (23)
          !!! Z boson, flavor is not changing
          call sf_int%base_init (mask, [data%mass**2], [data%mass**2], &
               [data%mZ**2], hel_lock = hel_lock)
          sf_int%data => data
          call flv_z%init (Z_BOSON, data%model)
          call qn_z%init (flv_z, col0)
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
                call sf_int%add_state ([qn, qn_rad, qn_z])
                call it_hel%advance ()
             end do
             !  call pol%final ()
          end do
       case (24)
          call sf_int%base_init (mask, [data%mass**2], [data%m_out**2], &
               [data%mW**2], hel_lock = hel_lock)
          sf_int%data => data
          call flv_wp%init (W_BOSON, data%model)
          call flv_wm%init (- W_BOSON, data%model)
          call qn_wp%init (flv_wp, col0)
          call qn_wm%init (flv_wm, col0)
          do i = 1, size (data%flv_in)
             isospin = data%flv_in(i)%get_isospin_type ()
             if (isospin > 0) then
                !!! up-type quark or neutrinos
                if (data%flv_in(i)%is_antiparticle ()) then
                   qn_w = qn_wm
                else
                   qn_w = qn_wp
                end if
             else
                !!! down-type quark or lepton
                if (data%flv_in(i)%is_antiparticle ()) then
                   qn_w = qn_wp
                else
                   qn_w = qn_wm
                end if
             end if
             call pol%init_generic (data%flv_in(i))
             call qn_fc(1)%init ( &
                  flv = data%flv_in(i), &
                  col = color_from_flavor (data%flv_in(i), 1))
             call qn_fc_fin(1)%init ( &
                  flv = data%flv_out(i), &
                  col = color_from_flavor (data%flv_out(i), 1))
             call it_hel%init (pol)
             do while (it_hel%is_valid ())
                qn_hel = it_hel%get_quantum_numbers ()
                qn = qn_hel .merge. qn_fc(1)
                qn_rad = qn_hel .merge. qn_fc_fin(1)
                call qn_rad%tag_radiated ()
                call sf_int%add_state ([qn, qn_rad, qn_w])
                call it_hel%advance ()
             end do
             ! call pol%final ()
          end do
       case default
          call msg_fatal ("EWA initialization failed: wrong particle type.")
       end select
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
  end subroutine ewa_init

  module subroutine ewa_setup_constants (sf_int)
    class(ewa_t), intent(inout), target :: sf_int
    type(state_iterator_t) :: it
    type(flavor_t) :: flv
    real(default) :: q, t3
    integer :: i
    sf_int%n_me = sf_int%get_n_matrix_elements ()
    allocate (sf_int%cv (sf_int%n_me))
    allocate (sf_int%ca (sf_int%n_me))
    associate (data => sf_int%data)
      select case (data%id)
      case (23)
         call it%init (sf_int%interaction_t%get_state_matrix_ptr ())
         do while (it%is_valid ())
            i = it%get_me_index ()
            flv = it%get_flavor (1)
            q = flv%get_charge ()
            t3 = flv%get_isospin ()
            if (flv%is_antiparticle ()) then
               sf_int%cv(i) = - data%cv &
                    * (t3 - 2._default * q * data%sinthw**2) / data%costhw
               sf_int%ca(i) = data%ca *  t3 / data%costhw
            else
               sf_int%cv(i) = data%cv &
                    * (t3 - 2._default * q * data%sinthw**2) / data%costhw
               sf_int%ca(i) = data%ca *  t3 / data%costhw
            end if
            call it%advance ()
         end do
      case (24)
         call it%init (sf_int%interaction_t%get_state_matrix_ptr ())
         do while (it%is_valid ())
            i = it%get_me_index ()
            flv = it%get_flavor (1)
            if (flv%is_antiparticle ()) then
               sf_int%cv(i) = data%cv / sqrt(2._default)
               sf_int%ca(i) = - data%ca / sqrt(2._default)
            else
               sf_int%cv(i) = data%cv / sqrt(2._default)
               sf_int%ca(i) = data%ca / sqrt(2._default)
            end if
            call it%advance ()
         end do
      end select
    end associate
    sf_int%status = SF_INITIAL
  end subroutine ewa_setup_constants

  module subroutine ewa_complete_kinematics (sf_int, x, xb, f, r, rb, map)
    class(ewa_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: r
    real(default), dimension(:), intent(in) :: rb
    logical, intent(in) :: map
    real(default) :: e_1
    real(default) :: x0, x1, lx0, lx1, lx
    e_1 = energy (sf_int%get_momentum (1))
    if (sf_int%data%recoil) then
       select case (sf_int%data%id)
       case (23)
          x0 = max (sf_int%data%x_min, sf_int%data%mz / e_1)
       case (24)
          x0 = max (sf_int%data%x_min, sf_int%data%mw / e_1)
       end select
    else
       x0 = sf_int%data%x_min
    end if
    x1 = sf_int%data%x_max
    if ( x0 >= x1) then
       f = 0
       sf_int%status = SF_FAILED_KINEMATICS
       return
    end if
    if (map) then
       lx0 = log (x0)
       lx1 = log (x1)
       lx = lx1 * r(1) + lx0 * rb(1)
       x(1) = exp(lx)
       f = x(1) * (lx1 - lx0)
    else
       x(1) = r(1)
       if (x0 < x(1) .and. x(1) < x1) then
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
       sf_int%x  = x(1)
       sf_int%xb = xb(1)
    case (SF_FAILED_KINEMATICS)
       sf_int%x  = 0
       sf_int%xb = 0
       f = 0
    end select
  end subroutine ewa_complete_kinematics

  module subroutine sf_ewa_recover_x (sf_int, x, xb, x_free)
    class(ewa_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(inout), optional :: x_free
    call sf_int%base_recover_x (x, xb, x_free)
    sf_int%x  = x(1)
    sf_int%xb = xb(1)
  end subroutine sf_ewa_recover_x

  module subroutine ewa_inverse_kinematics &
       (sf_int, x, xb, f, r, rb, map, set_momenta)
    class(ewa_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: r
    real(default), dimension(:), intent(out) :: rb
    logical, intent(in) :: map
    logical, intent(in), optional :: set_momenta
    real(default) :: x0, x1, lx0, lx1, lx, e_1
    logical :: set_mom
    set_mom = .false.;  if (present (set_momenta))  set_mom = set_momenta
    e_1 = energy (sf_int%get_momentum (1))
    if (sf_int%data%recoil) then
       select case (sf_int%data%id)
       case (23)
          x0 = max (sf_int%data%x_min, sf_int%data%mz / e_1)
       case (24)
          x0 = max (sf_int%data%x_min, sf_int%data%mw / e_1)
       end select
    else
       x0 = sf_int%data%x_min
    end if
    x1 = sf_int%data%x_max
    if (map) then
       lx0 = log (x0)
       lx1 = log (x1)
       lx = log (x(1))
       r(1)  = (lx - lx0) / (lx1 - lx0)
       rb(1) = (lx1 - lx) / (lx1 - lx0)
       f = x(1) * (lx1 - lx0)
    else
       r (1) = x(1)
       rb(1) = 1 - x(1)
       if (x0 < x(1) .and. x(1) < x1) then
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
  end subroutine ewa_inverse_kinematics

  module subroutine ewa_apply (sf_int, scale, negative_sf, rescale, i_sub)
    class(ewa_t), intent(inout) :: sf_int
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(in), optional :: rescale
    integer, intent(in), optional :: i_sub
    real(default) :: x, xb, pt2, c1, c2
    real(default) :: cv, ca
    real(default) :: f, fm, fp, fL
    integer :: i
    associate (data => sf_int%data)
      x  = sf_int%x
      xb = sf_int%xb
      pt2 = min ((data%pt_max)**2, (xb * data%sqrts / 2)**2)
      select case (data%id)
      case (23)
         !!! Z boson structure function
         c1 = log (1 + pt2 / (xb * (data%mZ)**2))
         c2 = 1 / (1 + (xb * (data%mZ)**2) / pt2)
      case (24)
         !!! W boson structure function
         c1 = log (1 + pt2 / (xb * (data%mW)**2))
         c2 = 1 / (1 + (xb * (data%mW)**2) / pt2)
      end select
      do i = 1, sf_int%n_me
         cv = sf_int%cv(i)
         ca = sf_int%ca(i)
         fm = data%coeff * &
              ((cv + ca)**2 + ((cv - ca) * xb)**2) * (c1 - c2) / (2 * x)
         fp = data%coeff * &
              ((cv - ca)**2 + ((cv + ca) * xb)**2) * (c1 - c2) / (2 * x)
         fL = data%coeff * &
              (cv**2 + ca**2) * (2 * xb / x) * c2
         f = fp + fm + fL
         if (.not. vanishes (f)) then
            fp = fp / f
            fm = fm / f
            fL = fL / f
         end if
         call sf_int%set_matrix_element (i, cmplx (f, kind=default))
      end do
    end associate
    sf_int%status = SF_EVALUATED
  end subroutine ewa_apply


end submodule sf_ewa_s

