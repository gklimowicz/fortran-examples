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

submodule (phs_single) phs_single_s

  use io_units
  use constants
  use numeric_utils
  use diagnostics
  use physics_defs

  implicit none

contains

  module subroutine phs_single_config_final (object)
    class(phs_single_config_t), intent(inout) :: object
  end subroutine phs_single_config_final

  module subroutine phs_single_config_write (object, unit, include_id)
    class(phs_single_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: include_id
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Partonic phase-space configuration (single-particle):"
    call object%base_write (unit)
  end subroutine phs_single_config_write

  module subroutine phs_single_config_configure (phs_config, sqrts, &
       sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
       ignore_mismatch, nlo_type, subdir)
    class(phs_single_config_t), intent(inout) :: phs_config
    real(default), intent(in) :: sqrts
    logical, intent(in), optional :: sqrts_fixed
    logical, intent(in), optional :: lab_is_cm
    logical, intent(in), optional :: azimuthal_dependence
    logical, intent(in), optional :: rebuild
    logical, intent(in), optional :: ignore_mismatch
    integer, intent(in), optional :: nlo_type
    type(string_t), intent(in), optional :: subdir
    if (.not. present (nlo_type)) &
      phs_config%nlo_type = BORN
    if (phs_config%n_out == 2) then
       phs_config%n_channel = 1
       phs_config%n_par = 2
       phs_config%sqrts = sqrts
       if (present (sqrts_fixed))  phs_config%sqrts_fixed = sqrts_fixed
       if (present (lab_is_cm))  phs_config%lab_is_cm = lab_is_cm
       if (present (azimuthal_dependence)) then
          phs_config%azimuthal_dependence = azimuthal_dependence
          if (.not. azimuthal_dependence) then
             allocate (phs_config%dim_flat (1))
             phs_config%dim_flat(1) = 2
          end if
       end if
       if (allocated (phs_config%channel))  deallocate (phs_config%channel)
       allocate (phs_config%channel (1))
       call phs_config%compute_md5sum ()
    else
       call msg_fatal ("Single-particle phase space requires n_out = 2")
    end if
  end subroutine phs_single_config_configure

  module subroutine phs_single_config_startup_message (phs_config, unit)
    class(phs_single_config_t), intent(in) :: phs_config
    integer, intent(in), optional :: unit
    call phs_config%base_startup_message (unit)
    write (msg_buffer, "(A,2(1x,I0,1x,A))") &
         "Phase space: single-particle"
    call msg_message (unit = unit)
  end subroutine phs_single_config_startup_message

  module subroutine phs_single_write (object, unit, verbose)
    class(phs_single_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    u = given_output_unit (unit)
    call object%base_write (u)
  end subroutine phs_single_write

  module subroutine phs_single_final (object)
    class(phs_single_t), intent(inout) :: object
  end subroutine phs_single_final

  module subroutine phs_single_init (phs, phs_config)
    class(phs_single_t), intent(out) :: phs
    class(phs_config_t), intent(in), target :: phs_config
    call phs%base_init (phs_config)
    phs%volume = 1 / (4 * twopi5)
    call phs%compute_factor ()
  end subroutine phs_single_init

  module subroutine phs_single_compute_factor (phs)
    class(phs_single_t), intent(inout) :: phs
    real(default) :: s_hat
    select case (phs%config%n_in)
    case (1)
       if (.not. phs%p_defined) then
          if (sum (phs%m_out) < phs%m_in(1)) then
             s_hat = phs%m_in(1) ** 2
             phs%f(1) = 1 / s_hat &
                  * sqrt (lambda (s_hat, phs%m_out(1)**2, phs%m_out(2)**2))
          else
             print *, "m_in  = ", phs%m_in
             print *, "m_out = ", phs%m_out
             call msg_fatal ("Decay is kinematically forbidden")
          end if
       end if
    case (2)
       if (phs%config%sqrts_fixed) then
          if (phs%p_defined)  return
          s_hat = phs%config%sqrts ** 2
       else
          if (.not. phs%p_defined)  return
          s_hat = sum (phs%p) ** 2
       end if
       if (sum (phs%m_in)**2 < s_hat .and. sum (phs%m_out)**2 < s_hat) then
          phs%f(1) = 1 / s_hat * &
               ( lambda (s_hat, phs%m_in (1)**2, phs%m_in (2)**2)   &
               * lambda (s_hat, phs%m_out(1)**2, phs%m_out(2)**2) ) &
               ** 0.25_default
       else
          phs%f(1) = 0
       end if
    end select
  end subroutine phs_single_compute_factor

  module subroutine phs_single_evaluate_selected_channel (phs, c_in, r_in)
    class(phs_single_t), intent(inout) :: phs
    integer, intent(in) :: c_in
    real(default), intent(in), dimension(:) :: r_in
    if (phs%p_defined) then
       call phs%select_channel (c_in)
       phs%r(:,c_in) = r_in
       select case (phs%config%n_in)
       case (2)
          if (all (phs%m_in == phs%m_out)) then
             call compute_kinematics_solid_angle (phs%p, phs%q, r_in)
          else
             call msg_bug ("PHS single: inelastic scattering not implemented")
          end if
       case (1)
          call compute_kinematics_solid_angle (phs%decay_p (), phs%q, r_in)
       end select
       call phs%compute_factor ()
       phs%q_defined = .true.
       phs%r_defined = .true.
    end if
  end subroutine phs_single_evaluate_selected_channel

  module subroutine phs_single_evaluate_other_channels (phs, c_in)
    class(phs_single_t), intent(inout) :: phs
    integer, intent(in) :: c_in
  end subroutine phs_single_evaluate_other_channels

  module function phs_single_decay_p (phs) result (p)
    class(phs_single_t), intent(in) :: phs
    type(vector4_t), dimension(2) :: p
    real(default) :: k
    real(default), dimension(2) :: E
    k = sqrt (lambda (phs%m_in(1) ** 2, phs%m_out(1) ** 2, phs%m_out(2) ** 2)) &
         / (2 * phs%m_in(1))
    E = sqrt (phs%m_out ** 2 + k ** 2)
    p(1) = vector4_moving (E(1), k, 3)
    p(2) = vector4_moving (E(2),-k, 3)
  end function phs_single_decay_p

  module subroutine phs_single_inverse (phs)
    class(phs_single_t), intent(inout) :: phs
    real(default), dimension(:), allocatable :: x
    if (phs%p_defined .and. phs%q_defined) then
       call phs%select_channel ()
       allocate (x (phs%config%n_par))
       call inverse_kinematics_solid_angle (phs%p, phs%q, x)
       phs%r(:,1) = x
       call phs%compute_factor ()
       phs%r_defined = .true.
    end if
  end subroutine phs_single_inverse


end submodule phs_single_s

