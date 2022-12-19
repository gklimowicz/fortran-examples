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

submodule (observables) observables_s

  use io_units
  use diagnostics
  use lorentz

  implicit none

contains

  integer function obs_pdg1 (prt1) result (pdg)
    type(prt_t), intent(in) :: prt1
    pdg = prt_get_pdg (prt1)
  end function obs_pdg1

  integer function obs_helicity1 (prt1) result (h)
    type(prt_t), intent(in) :: prt1
    if (prt_is_polarized (prt1)) then
       h = prt_get_helicity (prt1)
    else
       h = -9
    end if
  end function obs_helicity1

  integer function obs_n_col1 (prt1) result (n)
    type(prt_t), intent(in) :: prt1
    if (prt_is_colorized (prt1)) then
       n = prt_get_n_col (prt1)
    else
       n = 0
    end if
  end function obs_n_col1

  integer function obs_n_acl1 (prt1) result (n)
    type(prt_t), intent(in) :: prt1
    if (prt_is_colorized (prt1)) then
       n = prt_get_n_acl (prt1)
    else
       n = 0
    end if
  end function obs_n_acl1

  real(default) function obs_mass_squared1 (prt1) result (p2)
    type(prt_t), intent(in) :: prt1
    p2 = prt_get_msq (prt1)
  end function obs_mass_squared1

  real(default) function obs_signed_mass1 (prt1) result (m)
    type(prt_t), intent(in) :: prt1
    real(default) :: msq
    msq = prt_get_msq (prt1)
    m = sign (sqrt (abs (msq)), msq)
  end function obs_signed_mass1

  real(default) function obs_energy1 (prt1) result (e)
    type(prt_t), intent(in) :: prt1
    e = energy (prt_get_momentum (prt1))
  end function obs_energy1

  real(default) function obs_px1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = vector4_get_component (prt_get_momentum (prt1), 1)
  end function obs_px1

  real(default) function obs_py1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = vector4_get_component (prt_get_momentum (prt1), 2)
  end function obs_py1

  real(default) function obs_pz1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = vector4_get_component (prt_get_momentum (prt1), 3)
  end function obs_pz1

  real(default) function obs_p1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = space_part_norm (prt_get_momentum (prt1))
  end function obs_p1

  real(default) function obs_pl1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = longitudinal_part (prt_get_momentum (prt1))
  end function obs_pl1

  real(default) function obs_pt1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = transverse_part (prt_get_momentum (prt1))
  end function obs_pt1

  real(default) function obs_theta1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = polar_angle (prt_get_momentum (prt1))
  end function obs_theta1

  real(default) function obs_phi1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = azimuthal_angle (prt_get_momentum (prt1))
  end function obs_phi1

  real(default) function obs_rap1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = rapidity (prt_get_momentum (prt1))
  end function obs_rap1

  real(default) function obs_eta1 (prt1) result (p)
    type(prt_t), intent(in) :: prt1
    p = pseudorapidity (prt_get_momentum (prt1))
  end function obs_eta1

  real(default) function obs_theta_star1 (prt1) result (dist)
    type(prt_t), intent(in) :: prt1
    call msg_fatal (" 'Theta_star' is undefined as unary observable")
    dist = 0
  end function obs_theta_star1

  real(default) function obs_dist1 (prt1) result (dist)
    type(prt_t), intent(in) :: prt1
    call msg_fatal (" 'Dist' is undefined as unary observable")
    dist = 0
  end function obs_dist1

  integer function obs_pdg2 (prt1, prt2) result (pdg)
    type(prt_t), intent(in) :: prt1, prt2
    call msg_fatal (" PDG_Code is undefined as binary observable")
    pdg = 0
  end function obs_pdg2

  integer function obs_helicity2 (prt1, prt2) result (h)
    type(prt_t), intent(in) :: prt1, prt2
    call msg_fatal (" Helicity is undefined as binary observable")
    h = 0
  end function obs_helicity2

  integer function obs_n_col2 (prt1, prt2) result (n)
    type(prt_t), intent(in) :: prt1, prt2
    call msg_fatal (" Ncol is undefined as binary observable")
    n = 0
  end function obs_n_col2

  integer function obs_n_acl2 (prt1, prt2) result (n)
    type(prt_t), intent(in) :: prt1, prt2
    call msg_fatal (" Nacl is undefined as binary observable")
    n = 0
  end function obs_n_acl2

  real(default) function obs_mass_squared2 (prt1, prt2) result (p2)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    p2 = prt_get_msq (prt)
  end function obs_mass_squared2

  real(default) function obs_signed_mass2 (prt1, prt2) result (m)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    real(default) :: msq
    call prt_init_combine (prt, prt1, prt2)
    msq = prt_get_msq (prt)
    m = sign (sqrt (abs (msq)), msq)
  end function obs_signed_mass2

  real(default) function obs_energy2 (prt1, prt2) result (e)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    e = energy (prt_get_momentum (prt))
  end function obs_energy2

  real(default) function obs_px2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    p = vector4_get_component (prt_get_momentum (prt), 1)
  end function obs_px2

  real(default) function obs_py2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    p = vector4_get_component (prt_get_momentum (prt), 2)
  end function obs_py2

  real(default) function obs_pz2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    p = vector4_get_component (prt_get_momentum (prt), 3)
  end function obs_pz2

  real(default) function obs_p2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    p = space_part_norm (prt_get_momentum (prt))
  end function obs_p2

  real(default) function obs_pl2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    p = longitudinal_part (prt_get_momentum (prt))
  end function obs_pl2

  real(default) function obs_pt2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    p = transverse_part (prt_get_momentum (prt))
  end function obs_pt2

  real(default) function obs_theta2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    p = enclosed_angle (prt_get_momentum (prt1), prt_get_momentum (prt2))
  end function obs_theta2

  real(default) function obs_phi2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    p = azimuthal_distance (prt_get_momentum (prt1), prt_get_momentum (prt2))
  end function obs_phi2

  real(default) function obs_rap2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    p = rapidity_distance &
         (prt_get_momentum (prt1), prt_get_momentum (prt2))
  end function obs_rap2

  real(default) function obs_eta2 (prt1, prt2) result (p)
    type(prt_t), intent(in) :: prt1, prt2
    type(prt_t) :: prt
    call prt_init_combine (prt, prt1, prt2)
    p = pseudorapidity_distance &
         (prt_get_momentum (prt1), prt_get_momentum (prt2))
  end function obs_eta2

  real(default) function obs_theta_star2 (prt1, prt2) result (theta)
    type(prt_t), intent(in) :: prt1, prt2
    theta = enclosed_angle_rest_frame &
         (prt_get_momentum (prt1), &
         prt_get_momentum (prt1) + prt_get_momentum (prt2))
  end function obs_theta_star2

  real(default) function obs_dist2 (prt1, prt2) result (dist)
    type(prt_t), intent(in) :: prt1, prt2
    dist = eta_phi_distance &
         (prt_get_momentum (prt1), prt_get_momentum (prt2))
  end function obs_dist2

  real(default) function obs_ktmeasure (prt1, prt2) result (kt)
    type(prt_t), intent(in) :: prt1, prt2
    real (default) :: q2, e1, e2
    ! Normalized scale to one for now! (#67)
    q2 = 1
    e1 = energy (prt_get_momentum (prt1))
    e2 = energy (prt_get_momentum (prt2))
    kt = (2/q2) * min(e1**2,e2**2) *  &
         (1 - enclosed_angle_ct(prt_get_momentum (prt1), &
         prt_get_momentum (prt2)))
  end function obs_ktmeasure

  real(default) function obs_ht (sev) result (ht)
    type(subevt_t), intent(in) :: sev
    integer :: i, n
    type(prt_t) :: prt
    n = sev%get_length ()
    ht = 0
    do i = 1, n
       prt = sev%get_prt (i)
       ht = ht + &
            sqrt (obs_pt1(prt)**2 + obs_mass_squared1(prt))
    end do
  end function obs_ht

  module subroutine var_list_init_num_id (var_list, proc_id, num_id)
    type(var_list_t), intent(inout) :: var_list
    type(string_t), intent(in) :: proc_id
    integer, intent(in), optional :: num_id
    call var_list%set_procvar_int (proc_id, var_str ("num_id"), num_id)
  end subroutine var_list_init_num_id

  module subroutine var_list_init_process_results (var_list, proc_id, &
       n_calls, integral, error, accuracy, chi2, efficiency)
    type(var_list_t), intent(inout) :: var_list
    type(string_t), intent(in) :: proc_id
    integer, intent(in), optional :: n_calls
    real(default), intent(in), optional :: integral, error, accuracy
    real(default), intent(in), optional :: chi2, efficiency
    call var_list%set_procvar_real (proc_id, var_str ("integral"), integral)
    call var_list%set_procvar_real (proc_id, var_str ("error"), error)
  end subroutine var_list_init_process_results

  module subroutine var_list_set_observables_unary (var_list, prt1)
    type(var_list_t), intent(inout) :: var_list
    type(prt_t), intent(in), target :: prt1
    call var_list%append_obs1_iptr (var_str ("PDG"), obs_pdg1, prt1)
    call var_list%append_obs1_iptr (var_str ("Hel"), obs_helicity1, prt1)
    call var_list%append_obs1_iptr (var_str ("Ncol"), obs_n_col1, prt1)
    call var_list%append_obs1_iptr (var_str ("Nacl"), obs_n_acl1, prt1)
    call var_list%append_obs1_rptr &
         (var_str ("M"), obs_signed_mass1, prt1)
    call var_list%append_obs1_rptr &
         (var_str ("M2"), obs_mass_squared1, prt1)
    call var_list%append_obs1_rptr (var_str ("E"), obs_energy1, prt1)
    call var_list%append_obs1_rptr (var_str ("Px"), obs_px1, prt1)
    call var_list%append_obs1_rptr (var_str ("Py"), obs_py1, prt1)
    call var_list%append_obs1_rptr (var_str ("Pz"), obs_pz1, prt1)
    call var_list%append_obs1_rptr (var_str ("P"), obs_p1, prt1)
    call var_list%append_obs1_rptr (var_str ("Pl"), obs_pl1, prt1)
    call var_list%append_obs1_rptr (var_str ("Pt"), obs_pt1, prt1)
    call var_list%append_obs1_rptr (var_str ("Theta"), obs_theta1, prt1)
    call var_list%append_obs1_rptr (var_str ("Phi"), obs_phi1, prt1)
    call var_list%append_obs1_rptr (var_str ("Rap"), obs_rap1, prt1)
    call var_list%append_obs1_rptr (var_str ("Eta"), obs_eta1, prt1)
    call var_list%append_obs1_rptr &
         (var_str ("Theta_star"), obs_theta_star1, prt1)
    call var_list%append_obs1_rptr (var_str ("Dist"), obs_dist1, prt1)
    call var_list%append_uobs_real (var_str ("_User_obs_real"), prt1)
    call var_list%append_uobs_int (var_str ("_User_obs_int"), prt1)
  end subroutine var_list_set_observables_unary

  module subroutine var_list_set_observables_binary (var_list, prt1, prt2)
    type(var_list_t), intent(inout) :: var_list
    type(prt_t), intent(in), target :: prt1
    type(prt_t), intent(in), optional, target :: prt2
    call var_list%append_obs2_iptr (var_str ("PDG"), obs_pdg2, prt1, prt2)
    call var_list%append_obs2_iptr (var_str ("Hel"), obs_helicity2, prt1, prt2)
    call var_list%append_obs2_iptr (var_str ("Ncol"), obs_n_col2, prt1, prt2)
    call var_list%append_obs2_iptr (var_str ("Nacl"), obs_n_acl2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("M"), obs_signed_mass2, prt1, prt2)
    call var_list%append_obs2_rptr &
         (var_str ("M2"), obs_mass_squared2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("E"), obs_energy2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Px"), obs_px2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Py"), obs_py2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Pz"), obs_pz2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("P"), obs_p2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Pl"), obs_pl2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Pt"), obs_pt2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Theta"), obs_theta2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Phi"), obs_phi2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Rap"), obs_rap2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Eta"), obs_eta2, prt1, prt2)
    call var_list%append_obs2_rptr &
         (var_str ("Theta_star"), obs_theta_star2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("Dist"), obs_dist2, prt1, prt2)
    call var_list%append_obs2_rptr (var_str ("kT"), obs_ktmeasure, prt1, prt2)
    call var_list%append_uobs_real (var_str ("_User_obs_real"), prt1, prt2)
    call var_list%append_uobs_int (var_str ("_User_obs_int"), prt1, prt2)
  end subroutine var_list_set_observables_binary

  module subroutine var_list_set_observables_sev (var_list, pval)
    type(var_list_t), intent(inout) :: var_list
    type(subevt_t), intent(in), target:: pval
    call var_list%append_obsev_rptr (var_str ("Ht"), obs_ht, pval)
  end subroutine var_list_set_observables_sev

  module subroutine var_list_check_observable (var_list, name, type)
    class(var_list_t), intent(in), target :: var_list
    type(string_t), intent(in) :: name
    integer, intent(inout) :: type
    if (string_is_observable_id (name)) then
       call msg_fatal ("Variable name '" // char (name) &
            // "' is reserved for an observable")
       type = V_NONE
       return
    end if
  end subroutine var_list_check_observable

  function string_is_observable_id (string) result (flag)
    logical :: flag
    type(string_t), intent(in) :: string
    select case (char (string))
    case ("PDG", "Hel", "Ncol", "Nacl", &
         "M", "M2", "E", "Px", "Py", "Pz", "P", "Pl", "Pt", &
         "Theta", "Phi", "Rap", "Eta", "Theta_star", "Dist", "kT", &
         "Ht")
       flag = .true.
    case default
       flag = .false.
    end select
  end function string_is_observable_id

  module subroutine var_list_check_result_var (var_list, name, type)
    class(var_list_t), intent(in), target :: var_list
    type(string_t), intent(in) :: name
    integer, intent(inout) :: type
    if (string_is_integer_result_var (name))  type = V_INT
    if (.not. var_list%contains (name)) then
       if (string_is_result_var (name)) then
          call msg_fatal ("Result variable '" // char (name) // "' " &
               // "set without prior integration")
          type = V_NONE
          return
       else if (string_is_num_id (name)) then
          call msg_fatal ("Numeric process ID '" // char (name) // "' " &
               // "set without process declaration")
          type = V_NONE
          return
       end if
    end if
  end subroutine var_list_check_result_var

  function string_is_integer_result_var (string) result (flag)
    logical :: flag
    type(string_t), intent(in) :: string
    type(string_t) :: buffer, name, separator
    buffer = string
    call split (buffer, name, "(", separator=separator)  ! ")"
    if (separator == "(") then
       select case (char (name))
       case ("num_id", "n_calls")
          flag = .true.
       case default
          flag = .false.
       end select
    else
       flag = .false.
    end if
  end function string_is_integer_result_var

  function string_is_result_var (string) result (flag)
    logical :: flag
    type(string_t), intent(in) :: string
    type(string_t) :: buffer, name, separator
    buffer = string
    call split (buffer, name, "(", separator=separator)  ! ")"
    if (separator == "(") then
       select case (char (name))
       case ("integral", "error")
          flag = .true.
       case default
          flag = .false.
       end select
    else
       flag = .false.
    end if
  end function string_is_result_var

  function string_is_num_id (string) result (flag)
    logical :: flag
    type(string_t), intent(in) :: string
    type(string_t) :: buffer, name, separator
    buffer = string
    call split (buffer, name, "(", separator=separator)  ! ")"
    if (separator == "(") then
       select case (char (name))
       case ("num_id")
          flag = .true.
       case default
          flag = .false.
       end select
    else
       flag = .false.
    end if
  end function string_is_num_id


end submodule observables_s

