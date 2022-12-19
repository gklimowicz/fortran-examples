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

submodule (sf_base) sf_base_s

  use io_units
  use format_utils, only: write_separator
  use format_defs, only: FMT_17, FMT_19
  use constants, only: one, two
  use diagnostics
  use physics_defs, only: n_beams_rescaled

  implicit none

contains

  module subroutine sf_rescale_set_i_beam (func, i_beam)
    class(sf_rescale_t), intent(inout) :: func
    integer, intent(in) :: i_beam
    func%i_beam = i_beam
  end subroutine sf_rescale_set_i_beam

  module subroutine sf_rescale_collinear_apply (func, x)
    class(sf_rescale_collinear_t), intent(in) :: func
    real(default), intent(inout) :: x
    real(default) :: xi
    if (debug2_active (D_BEAMS)) then
       print *, 'Rescaling function - Collinear: '
       print *, 'Input, unscaled x: ', x
       print *, 'xi_tilde: ', func%xi_tilde
    end if
    xi = func%xi_tilde * (one - x)
    x = x / (one - xi)
    if (debug2_active (D_BEAMS))  print *, 'rescaled x: ', x
  end subroutine sf_rescale_collinear_apply

  module subroutine sf_rescale_collinear_set (func, xi_tilde)
    class(sf_rescale_collinear_t), intent(inout) :: func
    real(default), intent(in) :: xi_tilde
    func%xi_tilde = xi_tilde
  end subroutine sf_rescale_collinear_set

  module subroutine sf_rescale_real_apply (func, x)
    class(sf_rescale_real_t), intent(in) :: func
    real(default), intent(inout) :: x
    real(default) :: onepy, onemy
    if (debug2_active (D_BEAMS)) then
       print *, 'Rescaling function - Real: '
       print *, 'Input, unscaled: ', x
       print *, 'Beam index: ', func%i_beam
       print *, 'xi: ', func%xi, 'y: ', func%y
    end if
    x = x / sqrt (one - func%xi)
    onepy = one + func%y; onemy = one - func%y
    if (func%i_beam == 1) then
       x = x * sqrt ((two - func%xi * onemy) / (two - func%xi * onepy))
    else if (func%i_beam == 2) then
       x = x * sqrt ((two - func%xi * onepy) / (two - func%xi * onemy))
    else
       call msg_fatal ("sf_rescale_real_apply - invalid beam index")
    end if
    if (debug2_active (D_BEAMS))  print *, 'rescaled x: ', x
  end subroutine sf_rescale_real_apply

  module subroutine sf_rescale_real_set (func, xi, y)
    class(sf_rescale_real_t), intent(inout) :: func
    real(default), intent(in) :: xi, y
    func%xi = xi; func%y = y
  end subroutine sf_rescale_real_set

  module subroutine sf_rescale_dglap_apply (func, x)
    class(sf_rescale_dglap_t), intent(in) :: func
    real(default), intent(inout) :: x
    if (debug2_active (D_BEAMS))  then
       print *, "Rescaling function - DGLAP:"
       print *, "Input: ", x
       print *, "Beam index: ", func%i_beam
       print *, "z: ", func%z
    end if
    x = x / func%z(func%i_beam)
    if (debug2_active (D_BEAMS)) print *, "scaled x: ", x
  end subroutine sf_rescale_dglap_apply

  module subroutine sf_rescale_dglap_set (func, z)
    class(sf_rescale_dglap_t), intent(inout) :: func
    real(default), dimension(:), intent(in) :: z
    func%z = z
  end subroutine sf_rescale_dglap_set

  module function sf_data_is_generator (data) result (flag)
    class(sf_data_t), intent(in) :: data
    logical :: flag
    flag = .false.
  end function sf_data_is_generator

  elemental module function sf_data_get_pdf_set (data) result (pdf_set)
    class(sf_data_t), intent(in) :: data
    integer :: pdf_set
    pdf_set = 0
  end function sf_data_get_pdf_set

  module function sf_data_get_beam_file (data) result (file)
    class(sf_data_t), intent(in) :: data
    type(string_t) :: file
    file = ""
  end function sf_data_get_beam_file

  module subroutine sf_config_write (object, unit, verbose)
    class(sf_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    u = given_output_unit (unit)
    if (allocated (object%i)) then
       write (u, "(1x,A,2(1x,I0))")  "Structure-function configuration: &
            &beam(s)", object%i
       if (allocated (object%data)) &
            call object%data%write (u, verbose = verbose)
    else
       write (u, "(1x,A)")  "Structure-function configuration: [undefined]"
    end if
  end subroutine sf_config_write

  module subroutine sf_config_init (sf_config, i_beam, sf_data)
    class(sf_config_t), intent(out) :: sf_config
    integer, dimension(:), intent(in) :: i_beam
    class(sf_data_t), intent(in) :: sf_data
    allocate (sf_config%i (size (i_beam)), source = i_beam)
    allocate (sf_config%data, source = sf_data)
  end subroutine sf_config_init

  elemental module function sf_config_get_pdf_set (sf_config) result (pdf_set)
    class(sf_config_t), intent(in) :: sf_config
    integer :: pdf_set
    pdf_set = sf_config%data%get_pdf_set ()
  end function sf_config_get_pdf_set

  module function sf_config_get_beam_file (sf_config) result (file)
    class(sf_config_t), intent(in) :: sf_config
    type(string_t) :: file
    file = sf_config%data%get_beam_file ()
  end function sf_config_get_beam_file

  subroutine write_sf_status (status, u)
    integer, intent(in) :: status
    integer, intent(in) :: u
    select case (status)
    case (SF_UNDEFINED)
       write (u, "(1x,'[',A,']')")  "undefined"
    case (SF_INITIAL)
       write (u, "(1x,'[',A,']')")  "initialized"
    case (SF_DONE_LINKS)
       write (u, "(1x,'[',A,']')")  "links set"
    case (SF_FAILED_MASK)
       write (u, "(1x,'[',A,']')")  "mask mismatch"
    case (SF_DONE_MASK)
       write (u, "(1x,'[',A,']')")  "mask set"
    case (SF_FAILED_CONNECTIONS)
       write (u, "(1x,'[',A,']')")  "connections failed"
    case (SF_DONE_CONNECTIONS)
       write (u, "(1x,'[',A,']')")  "connections set"
    case (SF_SEED_KINEMATICS)
       write (u, "(1x,'[',A,']')")  "incoming momenta set"
    case (SF_FAILED_KINEMATICS)
       write (u, "(1x,'[',A,']')")  "kinematics failed"
    case (SF_DONE_KINEMATICS)
       write (u, "(1x,'[',A,']')")  "kinematics set"
    case (SF_FAILED_EVALUATION)
       write (u, "(1x,'[',A,']')")  "evaluation failed"
    case (SF_EVALUATED)
       write (u, "(1x,'[',A,']')")  "evaluated"
    end select
  end subroutine write_sf_status

  module subroutine sf_int_base_write (object, unit, testflag)
    class(sf_int_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "SF instance:"
    call write_sf_status (object%status, u)
    if (allocated (object%beam_index)) &
         write (u, "(3x,A,2(1x,I0))")  "beam      =", object%beam_index
    if (allocated (object%incoming)) &
         write (u, "(3x,A,2(1x,I0))")  "incoming  =", object%incoming
    if (allocated (object%radiated)) &
         write (u, "(3x,A,2(1x,I0))")  "radiated  =", object%radiated
    if (allocated (object%outgoing)) &
         write (u, "(3x,A,2(1x,I0))")  "outgoing  =", object%outgoing
    if (allocated (object%par_index)) &
         write (u, "(3x,A,2(1x,I0))")  "parameter =", object%par_index
    if (object%qmin_defined) &
         write (u, "(3x,A,1x," // FMT_19 // ")")  "q_min     =", object%qmin
    if (object%qmax_defined) &
         write (u, "(3x,A,1x," // FMT_19 // ")")  "q_max     =", object%qmax
    call object%interaction_t%basic_write (u, testflag = testflag)
  end subroutine sf_int_base_write

  module subroutine sf_int_base_init &
       (sf_int, mask, mi2, mr2, mo2, qmin, qmax, hel_lock)
    class(sf_int_t), intent(out) :: sf_int
    type (quantum_numbers_mask_t), dimension(:), intent(in) :: mask
    real(default), dimension(:), intent(in) :: mi2, mr2, mo2
    real(default), dimension(:), intent(in), optional :: qmin, qmax
    integer, dimension(:), intent(in), optional :: hel_lock
    allocate (sf_int%mi2 (size (mi2)))
    sf_int%mi2 = mi2
    allocate (sf_int%mr2 (size (mr2)))
    sf_int%mr2 = mr2
    allocate (sf_int%mo2 (size (mo2)))
    sf_int%mo2 = mo2
    if (present (qmin)) then
       sf_int%qmin_defined = .true.
       allocate (sf_int%qmin (size (qmin)))
       sf_int%qmin = qmin
    end if
    if (present (qmax)) then
       sf_int%qmax_defined = .true.
       allocate (sf_int%qmax (size (qmax)))
       sf_int%qmax = qmax
    end if
    call sf_int%interaction_t%basic_init &
         (size (mi2), 0, size (mr2) + size (mo2), &
         mask = mask, hel_lock = hel_lock, set_relations = .true.)
  end subroutine sf_int_base_init

  module subroutine sf_int_set_incoming (sf_int, incoming)
    class(sf_int_t), intent(inout) :: sf_int
    integer, dimension(:), intent(in) :: incoming
    allocate (sf_int%incoming (size (incoming)))
    sf_int%incoming = incoming
  end subroutine sf_int_set_incoming

  module subroutine sf_int_set_radiated (sf_int, radiated)
    class(sf_int_t), intent(inout) :: sf_int
    integer, dimension(:), intent(in) :: radiated
    allocate (sf_int%radiated (size (radiated)))
    sf_int%radiated = radiated
  end subroutine sf_int_set_radiated

  module subroutine sf_int_set_outgoing (sf_int, outgoing)
    class(sf_int_t), intent(inout) :: sf_int
    integer, dimension(:), intent(in) :: outgoing
    allocate (sf_int%outgoing (size (outgoing)))
    sf_int%outgoing = outgoing
  end subroutine sf_int_set_outgoing

  module subroutine sf_int_setup_constants (sf_int)
    class(sf_int_t), intent(inout), target :: sf_int
  end subroutine sf_int_setup_constants

  module subroutine sf_int_set_beam_index (sf_int, beam_index)
    class(sf_int_t), intent(inout) :: sf_int
    integer, dimension(:), intent(in) :: beam_index
    allocate (sf_int%beam_index (size (beam_index)))
    sf_int%beam_index = beam_index
  end subroutine sf_int_set_beam_index

  module subroutine sf_int_set_par_index (sf_int, par_index)
    class(sf_int_t), intent(inout) :: sf_int
    integer, dimension(:), intent(in) :: par_index
    allocate (sf_int%par_index (size (par_index)))
    sf_int%par_index = par_index
  end subroutine sf_int_set_par_index

  module subroutine sf_int_receive_momenta (sf_int)
    class(sf_int_t), intent(inout) :: sf_int
    if (sf_int%status >= SF_INITIAL) then
       call sf_int%receive_momenta ()
       sf_int%status = SF_SEED_KINEMATICS
    end if
  end subroutine sf_int_receive_momenta

  module subroutine sf_int_seed_momenta (sf_int, k)
    class(sf_int_t), intent(inout) :: sf_int
    type(vector4_t), dimension(:), intent(in) :: k
    if (sf_int%status >= SF_INITIAL) then
       call sf_int%set_momenta (k, outgoing=.false.)
       sf_int%status = SF_SEED_KINEMATICS
    end if
  end subroutine sf_int_seed_momenta

  module subroutine sf_int_seed_energies (sf_int, E)
    class(sf_int_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: E
    type(vector4_t), dimension(:), allocatable :: k
    integer :: j
    if (sf_int%status >= SF_INITIAL) then
       allocate (k (size (E)))
       if (all (E**2 >= sf_int%mi2)) then
          do j = 1, size (E)
             k(j) = vector4_moving (E(j), &
                  (3-2*j) * sqrt (E(j)**2 - sf_int%mi2(j)), 3)
          end do
          call sf_int%seed_kinematics (k)
       end if
    end if
  end subroutine sf_int_seed_energies

  module function sf_int_is_generator (sf_int) result (flag)
    class(sf_int_t), intent(in) :: sf_int
    logical :: flag
    flag = .false.
  end function sf_int_is_generator

  module subroutine sf_int_generate_free (sf_int, r, rb,  x_free)
    class(sf_int_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(inout) :: x_free
    r = 0
    rb= 1
  end subroutine sf_int_generate_free

  module subroutine sf_int_split_momentum (sf_int, x, xb)
    class(sf_int_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    type(vector4_t) :: k
    type(vector4_t), dimension(2) :: q
    type(splitting_data_t) :: sd
    real(default) :: E1, E2
    logical :: fail
    if (sf_int%status >= SF_SEED_KINEMATICS) then
       k = sf_int%get_momentum (1)
       call sd%init (k, &
            sf_int%mi2(1), sf_int%mr2(1), sf_int%mo2(1), &
            collinear = size (x) == 1)
       call sd%set_t_bounds (x(1), xb(1))
       select case (size (x))
       case (1)
       case (3)
          if (sf_int%qmax_defined) then
             if (sf_int%qmin_defined) then
                call sd%sample_t (x(2), &
                     t0 = - sf_int%qmax(1) ** 2, t1 = - sf_int%qmin(1) ** 2)
             else
                call sd%sample_t (x(2), &
                     t0 = - sf_int%qmax(1) ** 2)
             end if
          else
             if (sf_int%qmin_defined) then
                call sd%sample_t (x(2), t1 = - sf_int%qmin(1) ** 2)
             else
                call sd%sample_t (x(2))
             end if
          end if
          call sd%sample_phi (x(3))
       case default
          call msg_bug ("Structure function: impossible number of parameters")
       end select
       q = sd%split_momentum (k)
       call on_shell (q, [sf_int%mr2, sf_int%mo2], &
            sf_int%on_shell_mode)
       call sf_int%set_momenta (q, outgoing=.true.)
       E1 = energy (q(1))
       E2 = energy (q(2))
       fail = E1 < 0 .or. E2 < 0 &
            .or. E1 ** 2 < sf_int%mr2(1) &
            .or. E2 ** 2 < sf_int%mo2(1)
       if (fail) then
          sf_int%status = SF_FAILED_KINEMATICS
       else
          sf_int%status = SF_DONE_KINEMATICS
       end if
    end if
  end subroutine sf_int_split_momentum

  module subroutine sf_int_split_momenta (sf_int, x, xb)
    class(sf_int_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    type(vector4_t), dimension(2) :: k
    type(vector4_t), dimension(4) :: q
    real(default), dimension(4) :: E
    logical :: fail
    if (sf_int%status >= SF_SEED_KINEMATICS) then
       select case (size (x))
       case (2)
       case default
          call msg_bug ("Pair structure function: recoil requested &
               &but not implemented yet")
       end select
       k(1) = sf_int%get_momentum (1)
       k(2) = sf_int%get_momentum (2)
       q(1:2) = xb * k
       q(3:4) = x * k
       select case (size (sf_int%mr2))
       case (2)
          call on_shell (q, &
               [sf_int%mr2(1), sf_int%mr2(2), &
               sf_int%mo2(1), sf_int%mo2(2)], &
               sf_int%on_shell_mode)
          call sf_int%set_momenta (q, outgoing=.true.)
          E = energy (q)
          fail = any (E < 0) &
               .or. any (E(1:2) ** 2 < sf_int%mr2) &
               .or. any (E(3:4) ** 2 < sf_int%mo2)
       case default;  call msg_bug ("split momenta: incorrect use")
       end select
       if (fail) then
          sf_int%status = SF_FAILED_KINEMATICS
       else
          sf_int%status = SF_DONE_KINEMATICS
       end if
    end if
  end subroutine sf_int_split_momenta

  module subroutine sf_int_reduce_momenta (sf_int, x)
    class(sf_int_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    type(vector4_t), dimension(2) :: k
    type(vector4_t), dimension(2) :: q
    real(default), dimension(2) :: E
    logical :: fail
    if (sf_int%status >= SF_SEED_KINEMATICS) then
       select case (size (x))
       case (2)
       case default
          call msg_bug ("Pair spectrum: recoil requested &
               &but not implemented yet")
       end select
       k(1) = sf_int%get_momentum (1)
       k(2) = sf_int%get_momentum (2)
       q = x * k
       call on_shell (q, &
            [sf_int%mo2(1), sf_int%mo2(2)], &
            sf_int%on_shell_mode)
       call sf_int%set_momenta (q, outgoing=.true.)
       E = energy (q)
       fail = any (E < 0) &
            .or. any (E ** 2 < sf_int%mo2)
       if (fail) then
          sf_int%status = SF_FAILED_KINEMATICS
       else
          sf_int%status = SF_DONE_KINEMATICS
       end if
    end if
  end subroutine sf_int_reduce_momenta

  module subroutine sf_int_recover_x (sf_int, x, xb, x_free)
    class(sf_int_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(inout), optional :: x_free
    type(vector4_t), dimension(:), allocatable :: k
    type(vector4_t), dimension(:), allocatable :: q
    type(splitting_data_t) :: sd
    if (sf_int%status >= SF_SEED_KINEMATICS) then
       allocate (k (sf_int%interaction_t%get_n_in ()))
       allocate (q (sf_int%interaction_t%get_n_out ()))
       k = sf_int%get_momenta (outgoing=.false.)
       q = sf_int%get_momenta (outgoing=.true.)
       select case (size (k))
       case (1)
          call sd%init (k(1), &
               sf_int%mi2(1), sf_int%mr2(1), sf_int%mo2(1), &
               collinear = size (x) == 1)
          call sd%recover (k(1), q, sf_int%on_shell_mode)
          x(1)  = sd%get_x ()
          xb(1) = sd%get_xb ()
          select case (size (x))
          case (1)
          case (3)
             if (sf_int%qmax_defined) then
                if (sf_int%qmin_defined) then
                   call sd%inverse_t (x(2), &
                        t0 = - sf_int%qmax(1) ** 2, t1 = - sf_int%qmin(1) ** 2)
                else
                   call sd%inverse_t (x(2), &
                        t0 = - sf_int%qmax(1) ** 2)
                end if
             else
                if (sf_int%qmin_defined) then
                   call sd%inverse_t (x(2), t1 = - sf_int%qmin(1) ** 2)
                else
                   call sd%inverse_t (x(2))
                end if
             end if
             call sd%inverse_phi (x(3))
             xb(2:3) = 1 - x(2:3)
          case default
             call msg_bug ("Structure function: impossible number &
                  &of parameters")
          end select
       case (2)
          select case (size (x))
          case (2)
          case default
             call msg_bug ("Pair structure function: recoil requested &
                  &but not implemented yet")
          end select
          select case (sf_int%on_shell_mode)
          case (KEEP_ENERGY)
             select case (size (q))
             case (4)
                x = energy (q(3:4)) / energy (k)
                xb= energy (q(1:2)) / energy (k)
             case (2)
                x = energy (q) / energy (k)
                xb= 1 - x
             end select
          case (KEEP_MOMENTUM)
             select case (size (q))
             case (4)
                x = longitudinal_part (q(3:4)) / longitudinal_part (k)
                xb= longitudinal_part (q(1:2)) / longitudinal_part (k)
             case (2)
                x = longitudinal_part (q) / longitudinal_part (k)
                xb= 1 - x
             end select
          end select
       end select
    end if
  end subroutine sf_int_recover_x

  pure module function sf_int_get_n_in (object) result (n_in)
    class(sf_int_t), intent(in) :: object
    integer :: n_in
    n_in = object%interaction_t%get_n_in ()
  end function sf_int_get_n_in

  pure module function sf_int_get_n_rad (object) result (n_rad)
    class(sf_int_t), intent(in) :: object
    integer :: n_rad
    n_rad = object%interaction_t%get_n_out () &
         - object%interaction_t%get_n_in ()
  end function sf_int_get_n_rad

  pure module function sf_int_get_n_out (object) result (n_out)
    class(sf_int_t), intent(in) :: object
    integer :: n_out
    n_out = object%interaction_t%get_n_in ()
  end function sf_int_get_n_out

  module function sf_int_get_n_states (sf_int) result (n_states)
    class(sf_int_t), intent(in) :: sf_int
    integer :: n_states
    n_states = sf_int%get_n_matrix_elements ()
  end function sf_int_get_n_states

  module function sf_int_get_state (sf_int, i) result (qn)
    class(sf_int_t), intent(in) :: sf_int
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer, intent(in) :: i
    allocate (qn (sf_int%get_n_tot ()))
    qn = sf_int%get_quantum_numbers (i)
  end function sf_int_get_state

  module subroutine sf_int_get_values (sf_int, value)
    class(sf_int_t), intent(in) :: sf_int
    real(default), dimension(:), intent(out) :: value
    integer :: i
    if (sf_int%status >= SF_EVALUATED) then
       do i = 1, size (value)
          value(i) = real (sf_int%get_matrix_element (i))
       end do
    else
       value = 0
    end if
  end subroutine sf_int_get_values

  module subroutine sf_int_compute_values (sf_int, value, x, xb, scale, E)
    class(sf_int_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: value
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    real(default), intent(in) :: scale
    real(default), dimension(:), intent(in), optional :: E
    real(default), dimension(size (x)) :: xx, xxb
    real(default) :: f
    if (present (E))  call sf_int%seed_kinematics (E)
    if (sf_int%status >= SF_SEED_KINEMATICS) then
       call sf_int%complete_kinematics (xx, xxb, f, x, xb, map=.false.)
       call sf_int%apply (scale)
       call sf_int%get_values (value)
       value = value * f
    else
       value = 0
    end if
  end subroutine sf_int_compute_values

  module subroutine sf_int_compute_value &
       (sf_int, i_state, value, x, xb, scale, E)
    class(sf_int_t), intent(inout) :: sf_int
    integer, intent(in) :: i_state
    real(default), intent(out) :: value
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    real(default), intent(in) :: scale
    real(default), dimension(:), intent(in), optional :: E
    real(default), dimension(:), allocatable :: value_array
    if (sf_int%status >= SF_INITIAL) then
       allocate (value_array (sf_int%get_n_states ()))
       call sf_int%compute_values (value_array, x, xb, scale, E)
       value = value_array(i_state)
    else
       value = 0
    end if
  end subroutine sf_int_compute_value

  module subroutine sf_chain_final (object)
    class(sf_chain_t), intent(inout) :: object
    integer :: i
    call object%final_tracing ()
    if (allocated (object%sf)) then
       do i = 1, size (object%sf, 1)
          associate (sf => object%sf(i))
            if (allocated (sf%int)) then
               call sf%int%final ()
            end if
          end associate
       end do
    end if
    call beam_final (object%beam_t)
  end subroutine sf_chain_final

  module subroutine sf_chain_write (object, unit)
    class(sf_chain_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Incoming particles / structure-function chain:"
    if (associated (object%beam_data)) then
       write (u, "(3x,A,I0)")  "n_in      = ", object%n_in
       write (u, "(3x,A,I0)")  "n_strfun  = ", object%n_strfun
       write (u, "(3x,A,I0)")  "n_par     = ", object%n_par
       if (object%n_par /= object%n_bound) then
          write (u, "(3x,A,I0)")  "n_bound   = ", object%n_bound
       end if
       call object%beam_data%write (u)
       call write_separator (u)
       call beam_write (object%beam_t, u)
       if (allocated (object%sf)) then
          do i = 1, object%n_strfun
             associate (sf => object%sf(i))
               call write_separator (u)
               if (allocated (sf%int)) then
                  call sf%int%write (u)
               else
                  write (u, "(1x,A)")  "SF instance: [undefined]"
               end if
             end associate
          end do
       end if
    else
       write (u, "(3x,A)")  "[undefined]"
    end if
  end subroutine sf_chain_write

  module subroutine sf_chain_init (sf_chain, beam_data, sf_config)
    class(sf_chain_t), intent(out) :: sf_chain
    type(beam_data_t), intent(in), target :: beam_data
    type(sf_config_t), dimension(:), intent(in), optional, target :: sf_config
    integer :: i
    sf_chain%beam_data => beam_data
    sf_chain%n_in = beam_data%get_n_in ()
    call beam_init (sf_chain%beam_t, beam_data)
    if (present (sf_config)) then
       sf_chain%n_strfun = size (sf_config)
       allocate (sf_chain%sf (sf_chain%n_strfun))
       do i = 1, sf_chain%n_strfun
          call sf_chain%set_strfun (i, sf_config(i)%i, sf_config(i)%data)
       end do
    end if
  end subroutine sf_chain_init

  module subroutine sf_chain_receive_beam_momenta (sf_chain)
    class(sf_chain_t), intent(inout), target :: sf_chain
    type(interaction_t), pointer :: beam_int
    beam_int => sf_chain%get_beam_int_ptr ()
    call beam_int%receive_momenta ()
  end subroutine sf_chain_receive_beam_momenta

  module subroutine sf_chain_set_beam_momenta (sf_chain, p)
    class(sf_chain_t), intent(inout) :: sf_chain
    type(vector4_t), dimension(:), intent(in) :: p
    call beam_set_momenta (sf_chain%beam_t, p)
  end subroutine sf_chain_set_beam_momenta

  module subroutine sf_chain_set_strfun (sf_chain, i, beam_index, data)
    class(sf_chain_t), intent(inout) :: sf_chain
    integer, intent(in) :: i
    integer, dimension(:), intent(in) :: beam_index
    class(sf_data_t), intent(in), target :: data
    integer :: n_par, j
    n_par = data%get_n_par ()
    call data%allocate_sf_int (sf_chain%sf(i)%int)
    associate (sf_int => sf_chain%sf(i)%int)
      call sf_int%init (data)
      call sf_int%set_beam_index (beam_index)
      call sf_int%set_par_index &
           ([(j, j = sf_chain%n_par + 1, sf_chain%n_par + n_par)])
      sf_chain%n_par = sf_chain%n_par + n_par
      if (.not. data%is_generator ()) then
         sf_chain%n_bound = sf_chain%n_bound + n_par
      end if
    end associate
  end subroutine sf_chain_set_strfun

  module function sf_chain_get_n_par (sf_chain) result (n)
    class(sf_chain_t), intent(in) :: sf_chain
    integer :: n
    n = sf_chain%n_par
  end function sf_chain_get_n_par

  module function sf_chain_get_n_bound (sf_chain) result (n)
    class(sf_chain_t), intent(in) :: sf_chain
    integer :: n
    n = sf_chain%n_bound
  end function sf_chain_get_n_bound

  module function sf_chain_get_beam_int_ptr (sf_chain) result (int)
    type(interaction_t), pointer :: int
    class(sf_chain_t), intent(in), target :: sf_chain
    int => beam_get_int_ptr (sf_chain%beam_t)
  end function sf_chain_get_beam_int_ptr

  module subroutine sf_chain_setup_tracing (sf_chain, file)
    class(sf_chain_t), intent(inout) :: sf_chain
    type(string_t), intent(in) :: file
    if (sf_chain%n_strfun > 0) then
       sf_chain%trace_enable = .true.
       sf_chain%trace_unit = free_unit ()
       open (sf_chain%trace_unit, file = char (file), action = "write", &
            status = "replace")
       call sf_chain%write_trace_header ()
    else
       call msg_error ("Beam structure: no structure functions, tracing &
            &disabled")
    end if
  end subroutine sf_chain_setup_tracing

  module subroutine sf_chain_final_tracing (sf_chain)
    class(sf_chain_t), intent(inout) :: sf_chain
    if (sf_chain%trace_enable) then
       close (sf_chain%trace_unit)
       sf_chain%trace_enable = .false.
    end if
  end subroutine sf_chain_final_tracing

  module subroutine sf_chain_write_trace_header (sf_chain)
    class(sf_chain_t), intent(in) :: sf_chain
    integer :: u
    if (sf_chain%trace_enable) then
       u = sf_chain%trace_unit
       write (u, "('# ',A)")  "WHIZARD output: &
            &structure-function sampling data"
       write (u, "('# ',A,1x,I0)")  "Number of sf records:", sf_chain%n_strfun
       write (u, "('# ',A,1x,I0)")  "Number of parameters:", sf_chain%n_par
       write (u, "('# ',A)")  "Columns: channel, p(n_par), x(n_par), f, Jac * f"
    end if
  end subroutine sf_chain_write_trace_header

  module subroutine sf_chain_trace (sf_chain, c_sel, p, x, f, sf_sum)
    class(sf_chain_t), intent(in) :: sf_chain
    integer, intent(in) :: c_sel
    real(default), dimension(:,:), intent(in) :: p
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: f
    real(default), intent(in) :: sf_sum
    real(default) :: sf_sum_pac, f_sf_sum_pac
    integer :: u, i
    if (sf_chain%trace_enable) then
       u = sf_chain%trace_unit
       write (u, "(1x,I0)", advance="no")  c_sel
       write (u, "(2x)", advance="no")
       do i = 1, sf_chain%n_par
          write (u, "(1x," // FMT_17 // ")", advance="no")  p(i,c_sel)
       end do
       write (u, "(2x)", advance="no")
       do i = 1, sf_chain%n_par
          write (u, "(1x," // FMT_17 // ")", advance="no")  x(i)
       end do
       write (u, "(2x)", advance="no")
       sf_sum_pac = sf_sum
       f_sf_sum_pac = f(c_sel) * sf_sum
       call pacify (sf_sum_pac, 1.E-28_default)
       call pacify (f_sf_sum_pac, 1.E-28_default)
       write (u, "(2(1x," // FMT_17 // "))")  sf_sum_pac, f_sf_sum_pac
    end if
  end subroutine sf_chain_trace

  module subroutine sf_chain_instance_final (object)
    class(sf_chain_instance_t), intent(inout) :: object
    integer :: i
    if (allocated (object%sf)) then
       do i = 1, size (object%sf, 1)
          associate (sf => object%sf(i))
            if (allocated (sf%int)) then
               call sf%eval%final ()
               call sf%int%final ()
            end if
          end associate
       end do
    end if
    call beam_final (object%beam_t)
  end subroutine sf_chain_instance_final

  module subroutine sf_chain_instance_write (object, unit, col_verbose)
    class(sf_chain_instance_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: col_verbose
    integer :: u, i, c
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "Structure-function chain instance:"
    call write_sf_status (object%status, u)
    if (allocated (object%out_sf)) then
       write (u, "(3x,A)", advance="no")  "outgoing (interactions) ="
       do i = 1, size (object%out_sf)
          write (u, "(1x,I0,':',I0)", advance="no") &
               object%out_sf(i), object%out_sf_i(i)
       end do
       write (u, *)
    end if
    if (object%out_eval /= 0) then
       write (u, "(3x,A)", advance="no")  "outgoing (evaluators)   ="
       do i = 1, size (object%out_sf)
          write (u, "(1x,I0,':',I0)", advance="no") &
               object%out_eval, object%out_eval_i(i)
       end do
       write (u, *)
    end if
    if (allocated (object%sf)) then
       if (size (object%sf) /= 0) then
          write (u, "(1x,A)")  "Structure-function parameters:"
          do c = 1, size (object%f)
             write (u, "(1x,A,I0,A)", advance="no")  "Channel #", c, ":"
             if (c == object%selected_channel) then
                write (u, "(1x,A)")  "[selected]"
             else
                write (u, *)
             end if
             write (u, "(3x,A,9(1x,F9.7))")  "p =", object%p(:,c)
             write (u, "(3x,A,9(1x,F9.7))")  "pb=", object%pb(:,c)
             write (u, "(3x,A,9(1x,F9.7))")  "r =", object%r(:,c)
             write (u, "(3x,A,9(1x,F9.7))")  "rb=", object%rb(:,c)
             write (u, "(3x,A,9(1x,ES13.7))")  "f =", object%f(c)
             write (u, "(3x,A)", advance="no") "m ="
             call object%channel(c)%write (u)
          end do
          write (u, "(3x,A,9(1x,F9.7))")  "x =", object%x
          write (u, "(3x,A,9(1x,F9.7))")  "xb=", object%xb
          if (.not. all (object%bound)) then
             write (u, "(3x,A,9(1x,L1))")  "bound =", object%bound
          end if
       end if
    end if
    call write_separator (u)
    call beam_write (object%beam_t, u, col_verbose = col_verbose)
    if (allocated (object%sf)) then
       do i = 1, size (object%sf)
          associate (sf => object%sf(i))
            call write_separator (u)
            if (allocated (sf%int)) then
               if (allocated (sf%r)) then
                  write (u, "(1x,A)")  "Structure-function parameters:"
                  do c = 1, size (sf%f)
                     write (u, "(1x,A,I0,A)", advance="no")  "Channel #", c, ":"
                     if (c == object%selected_channel) then
                        write (u, "(1x,A)")  "[selected]"
                     else
                        write (u, *)
                     end if
                     write (u, "(3x,A,9(1x,F9.7))")  "r =", sf%r(:,c)
                     write (u, "(3x,A,9(1x,F9.7))")  "rb=", sf%rb(:,c)
                     write (u, "(3x,A,9(1x,ES13.7))")  "f =", sf%f(c)
                     write (u, "(3x,A,9(1x,L1,7x))") "m =", sf%m(c)
                  end do
                  write (u, "(3x,A,9(1x,F9.7))")  "x =", sf%x
                  write (u, "(3x,A,9(1x,F9.7))")  "xb=", sf%xb
               end if
               call sf%int%write(u)
               if (.not. sf%eval%is_empty ()) then
                     call sf%eval%write (u, col_verbose = col_verbose)
               end if
            end if
          end associate
       end do
    end if
  end subroutine sf_chain_instance_write

  module subroutine sf_chain_instance_init (chain, config, n_channel)
    class(sf_chain_instance_t), intent(out), target :: chain
    type(sf_chain_t), intent(in), target :: config
    integer, intent(in) :: n_channel
    integer :: i, j
    integer :: n_par_tot, n_par, n_strfun
    chain%config => config
    n_strfun = config%n_strfun
    chain%beam_t = config%beam_t
    allocate (chain%out_sf (config%n_in), chain%out_sf_i (config%n_in))
    allocate (chain%out_eval_i (config%n_in))
    chain%out_sf = 0
    chain%out_sf_i = [(i, i = 1, config%n_in)]
    chain%out_eval_i = chain%out_sf_i
    n_par_tot = 0
    if (n_strfun /= 0) then
       allocate (chain%sf (n_strfun))
       do i = 1, n_strfun
          associate (sf => chain%sf(i))
            allocate (sf%int, source=config%sf(i)%int)
            sf%int%interaction_t = config%sf(i)%int%interaction_t
            n_par = size (sf%int%par_index)
            allocate (sf%r (n_par, n_channel));  sf%r = 0
            allocate (sf%rb(n_par, n_channel));  sf%rb= 0
            allocate (sf%f (n_channel));         sf%f = 0
            allocate (sf%m (n_channel));         sf%m = .false.
            allocate (sf%x (n_par));             sf%x = 0
            allocate (sf%xb(n_par));             sf%xb= 0
            n_par_tot = n_par_tot + n_par
          end associate
       end do
       allocate (chain%p (n_par_tot, n_channel));  chain%p = 0
       allocate (chain%pb(n_par_tot, n_channel));  chain%pb= 0
       allocate (chain%r (n_par_tot, n_channel));  chain%r = 0
       allocate (chain%rb(n_par_tot, n_channel));  chain%rb= 0
       allocate (chain%f (n_channel));             chain%f = 0
       allocate (chain%x (n_par_tot));             chain%x = 0
       allocate (chain%xb(n_par_tot));             chain%xb= 0
       call allocate_sf_channels &
            (chain%channel, n_channel=n_channel, n_strfun=n_strfun)
    end if
    allocate (chain%bound (n_par_tot), source = .true.)
    do i = 1, n_strfun
       associate (sf => chain%sf(i))
         if (sf%int%is_generator ()) then
            do j = 1, size (sf%int%par_index)
               chain%bound(sf%int%par_index(j)) = .false.
            end do
         end if
       end associate
    end do
    chain%status = SF_INITIAL
  end subroutine sf_chain_instance_init

  module subroutine sf_chain_instance_select_channel (chain, channel)
    class(sf_chain_instance_t), intent(inout) :: chain
    integer, intent(in), optional :: channel
    if (present (channel)) then
       chain%selected_channel = channel
    else
       chain%selected_channel = 0
    end if
  end subroutine sf_chain_instance_select_channel

  module subroutine sf_chain_instance_set_channel (chain, c, channel)
    class(sf_chain_instance_t), intent(inout) :: chain
    integer, intent(in) :: c
    type(sf_channel_t), intent(in) :: channel
    integer :: i, j, k
    if (chain%status >= SF_INITIAL) then
       chain%channel(c) = channel
       j = 0
       do i = 1, chain%config%n_strfun
          associate (sf => chain%sf(i))
            sf%m(c) = channel%is_single_mapping (i)
            if (channel%is_multi_mapping (i)) then
               do k = 1, size (sf%int%beam_index)
                  j = j + 1
                  call chain%channel(c)%set_par_index &
                       (j, sf%int%par_index(k))
               end do
            end if
          end associate
       end do
       if (j /= chain%channel(c)%get_multi_mapping_n_par ()) then
          print *, "index last filled    = ", j
          print *, "number of parameters = ", &
               chain%channel(c)%get_multi_mapping_n_par ()
          call msg_bug ("Structure-function setup: mapping index mismatch")
       end if
       chain%status = SF_INITIAL
    end if
  end subroutine sf_chain_instance_set_channel

  module subroutine sf_chain_instance_link_interactions (chain)
    class(sf_chain_instance_t), intent(inout), target :: chain
    type(interaction_t), pointer :: int
    integer :: i, j, b
    if (chain%status >= SF_INITIAL) then
       do b = 1, chain%config%n_in
          int => beam_get_int_ptr (chain%beam_t)
          call interaction_set_source_link_beam (int, b, &
               chain%config%beam_t, b)
       end do
       if (allocated (chain%sf)) then
          do i = 1, size (chain%sf)
             associate (sf_int => chain%sf(i)%int)
               do j = 1, size (sf_int%beam_index)
                  b = sf_int%beam_index(j)
                  call link (sf_int%interaction_t, b, sf_int%incoming(j))
                  chain%out_sf(b) = i
                  chain%out_sf_i(b) = sf_int%outgoing(j)
               end do
             end associate
          end do
       end if
       chain%status = SF_DONE_LINKS
    end if
  contains
    subroutine link (int, b, in_index)
      type(interaction_t), intent(inout) :: int
      integer, intent(in) :: b, in_index
      integer :: i
      i = chain%out_sf(b)
      select case (i)
      case (0)
         call interaction_set_source_link_beam (int, in_index, &
              chain%beam_t, chain%out_sf_i(b))
      case default
         call int%set_source_link (in_index, &
              chain%sf(i)%int, chain%out_sf_i(b))
      end select
    end subroutine link
  end subroutine sf_chain_instance_link_interactions

  module subroutine sf_chain_exchange_mask (chain)
    class(sf_chain_instance_t), intent(inout), target :: chain
    type(interaction_t), pointer :: int
    type(quantum_numbers_mask_t), dimension(:), allocatable :: mask
    integer :: i
    if (chain%status >= SF_DONE_LINKS) then
       if (allocated (chain%sf)) then
          int => beam_get_int_ptr (chain%beam_t)
          allocate (mask (int%get_n_out ()))
          mask = int%get_mask ()
          if (size (chain%sf) /= 0) then
             do i = 1, size (chain%sf) - 1
                call chain%sf(i)%int%exchange_mask ()
             end do
             do i = size (chain%sf), 1, -1
                call chain%sf(i)%int%exchange_mask ()
             end do
             if (any (mask .neqv. int%get_mask ())) then
                chain%status = SF_FAILED_MASK
                return
             end if
             do i = 1, size (chain%sf)
                call chain%sf(i)%int%setup_constants ()
             end do
          end if
       end if
       chain%status = SF_DONE_MASK
    end if
  end subroutine sf_chain_exchange_mask

  module subroutine sf_chain_instance_init_evaluators (chain, extended_sf)
    class(sf_chain_instance_t), intent(inout), target :: chain
    logical, intent(in), optional :: extended_sf
    type(interaction_t), pointer :: int
    type(quantum_numbers_mask_t) :: mask
    integer :: i
    logical :: yorn
    yorn = .false.; if (present (extended_sf)) yorn = extended_sf
    if (chain%status >= SF_DONE_MASK) then
       if (allocated (chain%sf)) then
          if (size (chain%sf) /= 0) then
             mask = quantum_numbers_mask (.false., .false., .true.)
             int => beam_get_int_ptr (chain%beam_t)
             do i = 1, size (chain%sf)
                associate (sf => chain%sf(i))
                   if (yorn) then
                      if (int%get_n_sub () == 0) then
                         call int%declare_subtraction (n_beams_rescaled)
                      end if
                      if (sf%int%interaction_t%get_n_sub () == 0) then
                         call sf%int%interaction_t%declare_subtraction (n_beams_rescaled)
                      end if
                   end if
                   call sf%eval%init_product (int, sf%int%interaction_t, mask,&
                        & ignore_sub_for_qn = .true.)
                   if (sf%eval%is_empty ()) then
                      chain%status = SF_FAILED_CONNECTIONS
                      return
                   end if
                   int => sf%eval%interaction_t
                end associate
             end do
             call find_outgoing_particles ()
          end if
       else if (chain%out_eval == 0) then
          int => beam_get_int_ptr (chain%beam_t)
          call int%tag_hard_process ()
       end if
       chain%status = SF_DONE_CONNECTIONS
    end if
  contains
    subroutine find_outgoing_particles ()
      type(interaction_t), pointer :: int, int_next
      integer :: i, j, out_sf, out_i
      chain%out_eval = size (chain%sf)
      do j = 1, size (chain%out_eval_i)
         out_sf = chain%out_sf(j)
         out_i = chain%out_sf_i(j)
         if (out_sf == 0) then
            int => beam_get_int_ptr (chain%beam_t)
            out_sf = 1
         else
            int => chain%sf(out_sf)%int%interaction_t
         end if
         do i = out_sf, chain%out_eval
            int_next => chain%sf(i)%eval%interaction_t
            out_i = interaction_find_link (int_next, int, out_i)
            int => int_next
         end do
         chain%out_eval_i(j) = out_i
      end do
      call int%tag_hard_process (chain%out_eval_i)
    end subroutine find_outgoing_particles

  end subroutine sf_chain_instance_init_evaluators

  module subroutine sf_chain_instance_write_interaction &
       (chain, i_sf, i_int, unit)
    class(sf_chain_instance_t), intent(in) :: chain
    integer, intent(in) :: i_sf, i_int
    integer, intent(in) :: unit
    class(interaction_t), pointer :: int_in1 => null ()
    class(interaction_t), pointer :: int_in2 => null ()
    integer :: u
    u = given_output_unit (unit); if (u < 0) return
    if (chain%status >= SF_DONE_MASK) then
       if (allocated (chain%sf)) then
          int_in1 => evaluator_get_int_in_ptr (chain%sf(i_sf)%eval, 1)
          int_in2 => evaluator_get_int_in_ptr (chain%sf(i_sf)%eval, 2)
          if (int_in1%get_tag () == i_int) then
             call int_in1%basic_write (u)
          else if (int_in2%get_tag () == i_int) then
             call int_in2%basic_write (u)
          else
             write (u, "(A,1x,I0,1x,A,1x,I0)") 'No tag of sf', i_sf, 'matches' , i_int
          end if
       else
          write (u, "(A)") 'No sf_chain allocated!'
       end if
    else
       write (u, "(A)") 'sf_chain not ready!'
    end if
  end subroutine sf_chain_instance_write_interaction

  module subroutine sf_chain_instance_compute_kinematics (chain, c_sel, p_in)
    class(sf_chain_instance_t), intent(inout), target :: chain
    integer, intent(in) :: c_sel
    real(default), dimension(:), intent(in) :: p_in
    type(interaction_t), pointer :: int
    real(default) :: f_mapping
    integer :: i, j, c
    if (chain%status >= SF_DONE_CONNECTIONS) then
       call chain%select_channel (c_sel)
       int => beam_get_int_ptr (chain%beam_t)
       call int%receive_momenta ()
       if (allocated (chain%sf)) then
          if (size (chain%sf) /= 0) then
             forall (i = 1:size (chain%sf))  chain%sf(i)%int%status = SF_INITIAL
             chain%p (:,c_sel) = unpack (p_in, chain%bound, 0._default)
             chain%pb(:,c_sel) = 1 - chain%p(:,c_sel)
             chain%f = 1
             chain%x_free = 1
             do i = 1, size (chain%sf)
                associate (sf => chain%sf(i))
                  call sf%int%generate_free (sf%r(:,c_sel), sf%rb(:,c_sel), &
                       chain%x_free)
                  do j = 1, size (sf%x)
                     if (.not. chain%bound(sf%int%par_index(j))) then
                        chain%p (sf%int%par_index(j),c_sel) = sf%r (j,c_sel)
                        chain%pb(sf%int%par_index(j),c_sel) = sf%rb(j,c_sel)
                     end if
                  end do
                end associate
             end do
             if (allocated (chain%channel(c_sel)%multi_mapping)) then
                call chain%channel(c_sel)%multi_mapping%compute &
                     (chain%r(:,c_sel), chain%rb(:,c_sel), &
                      f_mapping, &
                      chain%p(:,c_sel), chain%pb(:,c_sel), &
                      chain%x_free)
                chain%f(c_sel) = f_mapping
             else
                chain%r (:,c_sel) = chain%p (:,c_sel)
                chain%rb(:,c_sel) = chain%pb(:,c_sel)
                chain%f(c_sel) = 1
             end if
             do i = 1, size (chain%sf)
                associate (sf => chain%sf(i))
                  call sf%int%seed_kinematics ()
                  do j = 1, size (sf%x)
                     sf%r (j,c_sel) = chain%r (sf%int%par_index(j),c_sel)
                     sf%rb(j,c_sel) = chain%rb(sf%int%par_index(j),c_sel)
                  end do
                  call sf%int%complete_kinematics &
                       (sf%x, sf%xb, sf%f(c_sel), sf%r(:,c_sel), sf%rb(:,c_sel), &
                        sf%m(c_sel))
                  do j = 1, size (sf%x)
                     chain%x(sf%int%par_index(j)) = sf%x(j)
                     chain%xb(sf%int%par_index(j)) = sf%xb(j)
                  end do
                  if (sf%int%status <= SF_FAILED_KINEMATICS) then
                     chain%status = SF_FAILED_KINEMATICS
                     return
                  end if
                  do c = 1, size (sf%f)
                     if (c /= c_sel) then
                        call sf%int%inverse_kinematics &
                             (sf%x, sf%xb, sf%f(c), sf%r(:,c), sf%rb(:,c), sf%m(c))
                        do j = 1, size (sf%x)
                           chain%r (sf%int%par_index(j),c) = sf%r (j,c)
                           chain%rb(sf%int%par_index(j),c) = sf%rb(j,c)
                        end do
                     end if
                     chain%f(c) = chain%f(c) * sf%f(c)
                  end do
                  if (.not. sf%eval%is_empty ()) then
                     call sf%eval%receive_momenta ()
                  end if
                end associate
             end do
             do c = 1, size (chain%f)
                if (c /= c_sel) then
                   if (allocated (chain%channel(c)%multi_mapping)) then
                      call chain%channel(c)%multi_mapping%inverse &
                           (chain%r(:,c), chain%rb(:,c), &
                            f_mapping, &
                            chain%p(:,c), chain%pb(:,c), &
                            chain%x_free)
                      chain%f(c) = chain%f(c) * f_mapping
                   else
                      chain%p (:,c) = chain%r (:,c)
                      chain%pb(:,c) = chain%rb(:,c)
                   end if
                end if
             end do
          end if
       end if
       chain%status = SF_DONE_KINEMATICS
    end if
  end subroutine sf_chain_instance_compute_kinematics

  module subroutine sf_chain_instance_inverse_kinematics (chain, x, xb)
    class(sf_chain_instance_t), intent(inout), target :: chain
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    type(interaction_t), pointer :: int
    real(default) :: f_mapping
    integer :: i, j, c
    if (chain%status >= SF_DONE_CONNECTIONS) then
       call chain%select_channel ()
       int => beam_get_int_ptr (chain%beam_t)
       call int%receive_momenta ()
       if (allocated (chain%sf)) then
          chain%f = 1
          if (size (chain%sf) /= 0) then
             forall (i = 1:size (chain%sf))  chain%sf(i)%int%status = SF_INITIAL
             chain%x = x
             chain%xb= xb
             do i = 1, size (chain%sf)
                associate (sf => chain%sf(i))
                  call sf%int%seed_kinematics ()
                  do j = 1, size (sf%x)
                     sf%x(j)  = chain%x(sf%int%par_index(j))
                     sf%xb(j) = chain%xb(sf%int%par_index(j))
                  end do
                  do c = 1, size (sf%f)
                     call sf%int%inverse_kinematics &
                          (sf%x, sf%xb, sf%f(c), sf%r(:,c), sf%rb(:,c), sf%m(c), &
                          set_momenta = c==1)
                     chain%f(c) = chain%f(c) * sf%f(c)
                     do j = 1, size (sf%x)
                        chain%r (sf%int%par_index(j),c) = sf%r (j,c)
                        chain%rb(sf%int%par_index(j),c) = sf%rb(j,c)
                     end do
                  end do
                  if (.not. sf%eval%is_empty ()) then
                     call sf%eval%receive_momenta ()
                  end if
                end associate
             end do
             do c = 1, size (chain%f)
                if (allocated (chain%channel(c)%multi_mapping)) then
                   call chain%channel(c)%multi_mapping%inverse &
                        (chain%r(:,c), chain%rb(:,c), &
                        f_mapping, &
                        chain%p(:,c), chain%pb(:,c), &
                        chain%x_free)
                   chain%f(c) = chain%f(c) * f_mapping
                else
                   chain%p (:,c) = chain%r (:,c)
                   chain%pb(:,c) = chain%rb(:,c)
                end if
             end do
          end if
       end if
       chain%status = SF_DONE_KINEMATICS
    end if
  end subroutine sf_chain_instance_inverse_kinematics

  module subroutine sf_chain_instance_recover_kinematics (chain, c_sel)
    class(sf_chain_instance_t), intent(inout), target :: chain
    integer, intent(in) :: c_sel
    real(default) :: f_mapping
    integer :: i, j, c
    if (chain%status >= SF_DONE_CONNECTIONS) then
       call chain%select_channel (c_sel)
       if (allocated (chain%sf)) then
          do i = size (chain%sf), 1, -1
             associate (sf => chain%sf(i))
               if (.not. sf%eval%is_empty ()) then
                  call sf%eval%send_momenta ()
               end if
             end associate
          end do
          chain%f = 1
          if (size (chain%sf) /= 0) then
             forall (i = 1:size (chain%sf))  chain%sf(i)%int%status = SF_INITIAL
             chain%x_free = 1
             do i = 1, size (chain%sf)
                associate (sf => chain%sf(i))
                  call sf%int%seed_kinematics ()
                  call sf%int%recover_x (sf%x, sf%xb, chain%x_free)
                  do j = 1, size (sf%x)
                     chain%x(sf%int%par_index(j))  = sf%x(j)
                     chain%xb(sf%int%par_index(j)) = sf%xb(j)
                  end do
                  do c = 1, size (sf%f)
                     call sf%int%inverse_kinematics &
                          (sf%x, sf%xb, sf%f(c), sf%r(:,c), sf%rb(:,c), sf%m(c), &
                          set_momenta = .false.)
                     chain%f(c) = chain%f(c) * sf%f(c)
                     do j = 1, size (sf%x)
                        chain%r (sf%int%par_index(j),c) = sf%r (j,c)
                        chain%rb(sf%int%par_index(j),c) = sf%rb(j,c)
                     end do
                  end do
                end associate
             end do
             do c = 1, size (chain%f)
                if (allocated (chain%channel(c)%multi_mapping)) then
                   call chain%channel(c)%multi_mapping%inverse &
                        (chain%r(:,c), chain%rb(:,c), &
                        f_mapping, &
                        chain%p(:,c), chain%pb(:,c), &
                        chain%x_free)
                   chain%f(c) = chain%f(c) * f_mapping
                else
                   chain%p (:,c) = chain%r (:,c)
                   chain%pb(:,c) = chain%rb(:,c)
                end if
             end do
          end if
       end if
       chain%status = SF_DONE_KINEMATICS
    end if
  end subroutine sf_chain_instance_recover_kinematics

  module subroutine sf_chain_instance_return_beam_momenta (chain)
    class(sf_chain_instance_t), intent(in), target :: chain
    type(interaction_t), pointer :: int
    if (chain%status >= SF_DONE_KINEMATICS) then
       int => beam_get_int_ptr (chain%beam_t)
       call int%send_momenta ()
    end if
  end subroutine sf_chain_instance_return_beam_momenta

  module subroutine sf_chain_instance_evaluate &
       (chain, scale, negative_sf, sf_rescale)
    class(sf_chain_instance_t), intent(inout), target :: chain
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(inout), optional :: sf_rescale
    type(interaction_t), pointer :: out_int
    real(default) :: sf_sum
    integer :: i_beam, i_sub, n_sub
    logical :: rescale
    n_sub = 0
    rescale = .false.; if (present (sf_rescale)) rescale = .true.
    if (rescale) then
       n_sub = chain%get_n_sub ()
    end if
    if (chain%status >= SF_DONE_KINEMATICS) then
       if (allocated (chain%sf)) then
          if (size (chain%sf) /= 0) then
             do i_beam = 1, size (chain%sf)
                associate (sf => chain%sf(i_beam))
                   if (rescale) then
                      call sf_rescale%set_i_beam (i_beam)
                      do i_sub = 0, n_sub
                         select case (i_sub)
                         case (0)
                            if (n_sub == 0) then
                               call sf%int%apply (scale, negative_sf, sf_rescale, i_sub = i_sub)
                            else
                               call sf%int%apply (scale, negative_sf, i_sub = i_sub)
                            end if
                         case default
                            if (i_beam == i_sub) then
                               call sf%int%apply (scale, negative_sf, sf_rescale, i_sub = i_sub)
                            else
                               call sf%int%apply (scale, negative_sf, i_sub = i_sub)
                            end if
                         end select
                      end do
                   else
                      call sf%int%apply (scale, negative_sf, i_sub = n_sub)
                   end if
                   if (sf%int%status <= SF_FAILED_EVALUATION) then
                      chain%status = SF_FAILED_EVALUATION
                      return
                   end if
                   if (.not. sf%eval%is_empty ())  call sf%eval%evaluate ()
                end associate
             end do
             out_int => chain%get_out_int_ptr ()
             sf_sum = real (out_int%sum ())
             call chain%config%trace &
                  (chain%selected_channel, chain%p, chain%x, chain%f, sf_sum)
          end if
       end if
       chain%status = SF_EVALUATED
    end if
  end subroutine sf_chain_instance_evaluate

  module subroutine sf_chain_instance_get_out_momenta (chain, p)
    class(sf_chain_instance_t), intent(in), target :: chain
    type(vector4_t), dimension(:), intent(out) :: p
    type(interaction_t), pointer :: int
    integer :: i, j
    if (chain%status >= SF_DONE_KINEMATICS) then
       do j = 1, size (chain%out_sf)
          i = chain%out_sf(j)
          select case (i)
          case (0)
             int => beam_get_int_ptr (chain%beam_t)
          case default
             int => chain%sf(i)%int%interaction_t
          end select
          p(j) = int%get_momentum (chain%out_sf_i(j))
       end do
    end if
  end subroutine sf_chain_instance_get_out_momenta

  module function sf_chain_instance_get_out_int_ptr (chain) result (int)
    class(sf_chain_instance_t), intent(in), target :: chain
    type(interaction_t), pointer :: int
    if (chain%out_eval == 0) then
       int => beam_get_int_ptr (chain%beam_t)
    else
       int => chain%sf(chain%out_eval)%eval%interaction_t
    end if
  end function sf_chain_instance_get_out_int_ptr

  module function sf_chain_instance_get_out_i (chain, j) result (i)
    class(sf_chain_instance_t), intent(in) :: chain
    integer, intent(in) :: j
    integer :: i
    i = chain%out_eval_i(j)
  end function sf_chain_instance_get_out_i

  module function sf_chain_instance_get_out_mask (chain) result (mask)
    class(sf_chain_instance_t), intent(in), target :: chain
    type(quantum_numbers_mask_t), dimension(:), allocatable :: mask
    type(interaction_t), pointer :: int
    allocate (mask (chain%config%n_in))
    int => chain%get_out_int_ptr ()
    mask = int%get_mask (chain%out_eval_i)
  end function sf_chain_instance_get_out_mask

  module subroutine sf_chain_instance_get_mcpar (chain, c, r)
    class(sf_chain_instance_t), intent(in) :: chain
    integer, intent(in) :: c
    real(default), dimension(:), intent(out) :: r
    if (allocated (chain%p))  r = pack (chain%p(:,c), chain%bound)
  end subroutine sf_chain_instance_get_mcpar

  module function sf_chain_instance_get_f (chain, c) result (f)
    class(sf_chain_instance_t), intent(in) :: chain
    integer, intent(in) :: c
    real(default) :: f
    if (allocated (chain%f)) then
       f = chain%f(c)
    else
       f = 1
    end if
  end function sf_chain_instance_get_f

  module function sf_chain_instance_get_status (chain) result (status)
    class(sf_chain_instance_t), intent(in) :: chain
    integer :: status
    status = chain%status
  end function sf_chain_instance_get_status

  module subroutine sf_chain_instance_get_matrix_elements (chain, i, ff)
     class(sf_chain_instance_t), intent(in) :: chain
     integer, intent(in) :: i
     real(default), intent(out), dimension(:), allocatable :: ff

     associate (sf => chain%sf(i))
        ff = real (sf%int%get_matrix_element ())
     end associate
  end subroutine sf_chain_instance_get_matrix_elements

  module function sf_chain_instance_get_beam_int_ptr (chain) result (int)
    type(interaction_t), pointer :: int
    class(sf_chain_instance_t), intent(in), target :: chain
    int => beam_get_int_ptr (chain%beam_t)
  end function sf_chain_instance_get_beam_int_ptr

  module function sf_chain_instance_get_n_sub (chain) result (n_sub)
    type(interaction_t), pointer :: int
    class(sf_chain_instance_t), intent(in), target :: chain
    integer :: n_sub
    int => beam_get_int_ptr (chain%beam_t)
    n_sub = int%get_n_sub ()
  end function sf_chain_instance_get_n_sub


end submodule sf_base_s

