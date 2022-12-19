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

module recola_wrapper

  use recola !NODEP!

  use kinds
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use constants, only: zero
  use diagnostics, only: msg_fatal, msg_message, msg_debug, msg_debug2, D_ME_METHODS
  use io_units, only: given_output_unit

  implicit none
  private

  public :: get_recola_particle_string
  public :: rclwrap_get_new_recola_id
  public :: rclwrap_get_current_recola_id
  public :: rclwrap_request_generate_processes
  public :: rclwrap_add_process
  public :: rclwrap_define_processes
  public :: rclwrap_generate_processes
  public :: rclwrap_compute_process
  public :: rclwrap_get_amplitude
  public :: rclwrap_get_squared_amplitude
  public :: rclwrap_set_pole_mass
  public :: rclwrap_set_onshell_mass
  public :: rclwrap_use_gfermi_scheme
  public :: rclwrap_set_light_fermions
  public :: rclwrap_set_light_fermion
  public :: rclwrap_unset_light_fermion
  public :: rclwrap_set_onshell_scheme
  public :: rclwrap_set_alpha_s
  public :: rclwrap_get_alpha_s
  public :: rclwrap_get_alpha
  public :: rclwrap_get_helicity_configurations
  public :: rclwrap_get_color_configurations
  public :: rclwrap_use_dim_reg_soft
  public :: rclwrap_use_mass_reg_soft
  public :: rclwrap_set_delta_uv
  public :: rclwrap_set_mu_uv
  public :: rclwrap_set_delta_ir
  public :: rclwrap_set_mu_ir
  public :: rclwrap_get_renormalization_scale
  public :: rclwrap_get_flavor_scheme
  public :: rclwrap_use_alpha0_scheme
  public :: rclwrap_use_alphaz_scheme
  public :: rclwrap_set_complex_mass_scheme
  public :: rclwrap_set_resonant_particle
  public :: rclwrap_switch_on_resonant_self_energies
  public :: rclwrap_switch_off_resonant_self_energies
  public :: rclwrap_set_draw_level_branches
  public :: rclwrap_set_print_level_amplitude
  public :: rclwrap_set_print_level_squared_amplitude
  public :: rclwrap_set_print_level_correlations
  public :: rclwrap_set_print_level_RAM
  public :: rclwrap_scale_coupling3
  public :: rclwrap_scale_coupling4
  public :: rclwrap_switch_off_coupling3
  public :: rclwrap_switch_off_coupling4
  public :: rclwrap_set_ifail
  public :: rclwrap_get_ifail
  public :: rclwrap_set_output_file
  public :: rclwrap_set_gs_power
  public :: rclwrap_select_gs_power_born_amp
  public :: rclwrap_unselect_gs_power_born_amp
  public :: rclwrap_select_gs_power_loop_amp
  public :: rclwrap_unselect_gs_power_loop_amp
  public :: rclwrap_select_all_gs_powers_born_amp
  public :: rclwrap_unselect_all_gs_powers_loop_amp
  public :: rclwrap_select_all_gs_powers_loop_amp
  public :: rclwrap_unselect_all_gs_powers_born_amp
  public :: rclwrap_set_resonant_squared_momentum
  public :: rclwrap_compute_running_alpha_s
  public :: rclwrap_set_dynamic_settings
  public :: rclwrap_rescale_process
  public :: rclwrap_get_polarized_squared_amplitude
  public :: rclwrap_compute_color_correlation
  public :: rclwrap_compute_all_color_correlations
  public :: rclwrap_rescale_color_correlation
  public :: rclwrap_rescale_all_color_correlations
  public :: rclwrap_get_color_correlation
  public :: rclwrap_compute_spin_correlation
  public :: rclwrap_rescale_spin_correlation
  public :: rclwrap_get_spin_correlation
  public :: rclwrap_compute_spin_color_correlation
  public :: rclwrap_rescale_spin_color_correlation
  public :: rclwrap_get_spin_color_correlation
  public :: rclwrap_get_momenta
  public :: rclwrap_reset_recola

  public :: rclwrap_is_active
  logical, parameter :: rclwrap_is_active = .true.


  type :: rcl_process_t
     private
     integer :: id
     type(string_t) :: process_string
     type(string_t) :: order
   contains
     procedure :: get_params => rcl_process_get_params
     procedure :: write => rcl_process_write
  end type rcl_process_t

  interface rcl_process_t
     module procedure new_rcl_process_t
  end interface

  type :: rcl_controller_t
     private
     logical :: active = .false.
     logical :: defined = .false.
     logical :: done = .false.
     integer :: recola_id = 0
     type(rcl_process_t), dimension (:), allocatable :: processes
     integer :: n_processes = 0
   contains
     procedure :: add_process => rcl_controller_add_process
     procedure :: define_processes => rcl_controller_define_processes
     procedure :: reset => rcl_controller_reset
     procedure :: write => rcl_controller_write
     procedure :: get_new_id => rcl_controller_get_new_id
     procedure :: get_current_id => rcl_controller_get_current_id
     procedure :: activate => rcl_controller_activate
     procedure :: generate_processes => rcl_controller_generate_processes
  end type rcl_controller_t


  type(rcl_controller_t), target, save :: rcl_controller


contains

  elemental function get_recola_particle_string (pdg) result (name)
    type(string_t) :: name
    integer, intent(in) :: pdg
    select case (pdg)
    case (1)
       name = var_str ("d")
    case (-1)
       name = var_str ("d~")
    case (2)
       name = var_str ("u")
    case (-2)
       name = var_str ("u~")
    case (3)
       name = var_str ("s")
    case (-3)
       name = var_str ("s~")
    case (4)
       name = var_str ("c")
    case (-4)
       name = var_str ("c~")
    case (5)
       name = var_str ("b")
    case (-5)
       name = var_str ("b~")
    case (6)
       name = var_str ("t")
    case (-6)
       name = var_str ("t~")
    case (11)
       name = var_str ("e-")
    case (-11)
       name = var_str ("e+")
    case (12)
       name = var_str ("nu_e")
    case (-12)
       name = var_str ("nu_e~")
    case (13)
       name = var_str ("mu-")
    case (-13)
       name = var_str ("mu+")
    case (14)
       name = var_str ("nu_mu")
    case (-14)
       name = var_str ("nu_mu~")
    case (15)
       name = var_str ("tau-")
    case (-15)
       name = var_str ("tau+")
    case (16)
       name = var_str ("nu_tau")
    case (-16)
       name = var_str ("nu_tau~")
    case (21)
       name = var_str ("g")
    case (22)
       name = var_str ("A")
    case (23)
       name = var_str ("Z")
    case (24)
       name = var_str ("W+")
    case (-24)
       name = var_str ("W-")
    case (25)
       name = var_str ("H")
    end select
  end function get_recola_particle_string

  subroutine rclwrap_define_process (id, process_string, order)
    integer, intent(in) :: id
    type(string_t), intent(in) :: process_string
    type(string_t), intent(in) :: order
    if (debug_on) call msg_debug2 (D_ME_METHODS, "define_process_rcl")
    call define_process_rcl (id, char (process_string), char (order))
  end subroutine rclwrap_define_process

  function new_rcl_process_t (id, process_string, order)
    integer, intent(in) :: id
    type(string_t), intent(in) :: process_string, order
    type(rcl_process_t) :: new_rcl_process_t
    new_rcl_process_t%id = id
    new_rcl_process_t%process_string = process_string
    new_rcl_process_t%order = order
  end function new_rcl_process_t

  subroutine rcl_process_get_params (prc, id, process_string, order)
    class(rcl_process_t), intent(in) :: prc
    integer, intent(out) :: id
    type(string_t), intent(out) :: process_string
    type(string_t), intent(out) :: order
    id = prc%id
    process_string = prc%process_string
    order = prc%order
  end subroutine rcl_process_get_params

  subroutine rcl_process_write (object, unit)
    class(rcl_process_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A,I0,2(1x,A,1x))") "RECOLA process:", &
         "id=", object%id, "process_string=", char(object%process_string), &
         "order=", char(object%order)
  end subroutine rcl_process_write

  subroutine rcl_controller_add_process (rcl, process)
    class(rcl_controller_t), intent(inout) :: rcl
    type(rcl_process_t), intent(in) :: process
    type(rcl_process_t), dimension (:), allocatable :: temp
    if (rcl%n_processes == size(rcl%processes)) then
      allocate( temp(2 * rcl%n_processes) )
      temp(:rcl%n_processes) = rcl%processes
      call move_alloc(temp, rcl%processes)
    end if
    rcl%processes(rcl%n_processes + 1) = process
    rcl%n_processes = rcl%n_processes + 1
  end subroutine rcl_controller_add_process

  subroutine rcl_controller_define_processes (rcl)
    class(rcl_controller_t), intent(inout) :: rcl
    integer :: id, i
    type(string_t) :: process_string
    type(string_t) :: order
    if (rcl%defined) then
      if (.not. rcl%done) call rclwrap_generate_processes ()
      if (debug_on) call msg_debug2 (D_ME_METHODS, "reset_recola_rcl")
      call reset_recola_rcl ()
    end if
    do i = 1, rcl%n_processes
      call rcl%processes(i)%get_params(id, process_string, order)
      call rclwrap_define_process (id, process_string, order)
    end do
    rcl%defined = .true.
    rcl%done = .false.
  end subroutine rcl_controller_define_processes

  subroutine rcl_controller_reset (rcl)
    class(rcl_controller_t), intent(inout) :: rcl
    if (rcl%active .or. rcl%done) then
       if (debug_on) call msg_debug2 (D_ME_METHODS, "reset_recola_rcl")
       if (allocated (rcl%processes)) deallocate (rcl%processes)
       call reset_recola_rcl ()
    end if
    rcl%active = .false.
    rcl%defined = .false.
    rcl%done = .false.
    rcl%recola_id = 0
    rcl%n_processes = 0
  end subroutine rcl_controller_reset

  subroutine rcl_controller_write (object, unit)
    class(rcl_controller_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A,2(1x,A,L1),2(1x,A,I0))")  "RECOLA controller:", &
         "active=", object%active, "done=", object%done, &
         "id=", object%recola_id, "n_processes=", object%n_processes
  end subroutine rcl_controller_write

  subroutine rcl_controller_get_new_id (object, id)
    class(rcl_controller_t), intent(inout) :: object
    integer, intent(out) :: id
    object%recola_id = object%recola_id + 1
    id = object%recola_id
  end subroutine rcl_controller_get_new_id

  subroutine rcl_controller_get_current_id (object, id)
    class(rcl_controller_t), intent(inout) :: object
    integer, intent(out) :: id
    id = object%recola_id
  end subroutine rcl_controller_get_current_id

  subroutine rcl_controller_activate (rcl)
    class(rcl_controller_t), intent(inout) :: rcl
    if ( .not. allocated(rcl%processes) ) allocate ( rcl%processes(10) )
    rcl_controller%active = .true.
  end subroutine rcl_controller_activate

  subroutine rcl_controller_generate_processes (rcl)
    class(rcl_controller_t), intent(inout) :: rcl
    if (rcl_controller%active) then
       if (.not. rcl_controller%done) then
          call msg_message ("Recola: preparing processes for integration")
          call generate_processes_rcl ()
          rcl_controller%done = .true.
       end if
    end if
  end subroutine rcl_controller_generate_processes

  subroutine rclwrap_get_new_recola_id (id)
    integer, intent(out) :: id
    call rcl_controller%get_new_id (id)
  end subroutine rclwrap_get_new_recola_id

  function rclwrap_get_current_recola_id () result (n)
    integer :: n
    call rcl_controller%get_current_id (n)
  end function rclwrap_get_current_recola_id

  subroutine rclwrap_request_generate_processes ()
    if (debug_on) call msg_debug2 (D_ME_METHODS, "request_generate_processes_rcl")
    call rcl_controller%activate ()
  end subroutine rclwrap_request_generate_processes

  subroutine rclwrap_add_process (id, process_string, order)
    integer, intent(in) :: id
    type(string_t), intent(in) :: process_string, order
    type(rcl_process_t) :: prc
    if (debug_on) call msg_debug2 (D_ME_METHODS, "add_process_rcl: id", id)
    prc = rcl_process_t (id, process_string, order)
    call rcl_controller%add_process (prc)
  end subroutine rclwrap_add_process

  subroutine rclwrap_define_processes ()
    if (debug_on) call msg_debug2 (D_ME_METHODS, "define_processes_rcl")
    call rcl_controller%define_processes ()
  end subroutine rclwrap_define_processes

  subroutine rclwrap_generate_processes ()
    if (debug_on) call msg_debug2 (D_ME_METHODS, "generate_processes_rcl")
    call rcl_controller%generate_processes ()
  end subroutine rclwrap_generate_processes

  subroutine rclwrap_compute_process (id, p, order, sqme)
    integer, intent(in) :: id
    real(double), intent(in), dimension(:,:) :: p
    character(len=*), intent(in) :: order
    real(double), intent(out), dimension(0:1), optional :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "compute_process_rcl")
    call compute_process_rcl (id, p, order, sqme)
  end subroutine rclwrap_compute_process

  subroutine rclwrap_get_amplitude (id, g_power, order, col, hel, amp)
    integer, intent(in) :: id, g_power
    character(len=*), intent(in) :: order
    integer, dimension(:), intent(in) :: col, hel
    complex(double), intent(out) :: amp
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_amplitude_rcl")
    call get_amplitude_rcl (id, g_power, order, col, hel, amp)
  end subroutine rclwrap_get_amplitude

  subroutine rclwrap_get_squared_amplitude (id, alphas_power, order, sqme)
    integer, intent(in) :: id, alphas_power
    character(len=*), intent(in) :: order
    real(double), intent(out) :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_squared_amplitude_rcl")
    call get_squared_amplitude_rcl (id, alphas_power, order, sqme)
  end subroutine rclwrap_get_squared_amplitude

  subroutine rclwrap_set_pole_mass (pdg_id, mass, width)
    integer, intent(in) :: pdg_id
    real(double), intent(in) :: mass, width
    if (debug_on) call msg_debug2 (D_ME_METHODS, "rclwrap_set_pole_mass of ", pdg_id)
    select case (abs(pdg_id))
    case (11)
       if (width > zero) &
          call msg_fatal ("Recola pole mass: Attempting to set non-zero electron width!")
       call set_pole_mass_electron_rcl (mass)
    case (13)
       call set_pole_mass_muon_rcl (mass, width)
    case (15)
       call set_pole_mass_tau_rcl (mass, width)
    case (1)
       if (width > zero) &
          call msg_fatal ("Recola pole mass: Attempting to set non-zero down-quark width!")
       call set_pole_mass_down_rcl (mass)
    case (2)
       if (width > zero) &
          call msg_fatal ("Recola pole mass: Attempting to set non-zero up-quark width!")
       call set_pole_mass_up_rcl (mass)
    case (3)
       if (width > zero) &
          call msg_fatal ("Recola pole mass: Attempting to set non-zero strange-quark width!")
       call set_pole_mass_strange_rcl (mass)
    case (4)
       call set_pole_mass_charm_rcl (mass, width)
    case (5)
       call set_pole_mass_bottom_rcl (mass, width)
    case (6)
       call set_pole_mass_top_rcl (mass, width)
    case (23)
       call set_pole_mass_z_rcl (mass, width)
    case (24)
       call set_pole_mass_w_rcl (mass, width)
    case (25)
       call set_pole_mass_h_rcl (mass, width)
    case default
       call msg_fatal ("Recola pole mass: Unsupported particle")
    end select
  end subroutine rclwrap_set_pole_mass

  subroutine rclwrap_set_onshell_mass (pdg_id, mass, width)
    integer, intent(in) :: pdg_id
    real(double), intent(in) :: mass, width
    if (debug_on) call msg_debug2 (D_ME_METHODS, "rclwrap_set_onshell_mass of ", pdg_id)
    select case (abs(pdg_id))
    case (23)
       call set_onshell_mass_z_rcl (mass, width)
    case (24)
       call set_onshell_mass_w_rcl (mass, width)
    case default
       call msg_fatal ("Recola onshell mass: Only for W and Z")
    end select
  end subroutine rclwrap_set_onshell_mass

  subroutine rclwrap_use_gfermi_scheme (gf)
    real(double), intent(in), optional :: gf
    if (debug_on) call msg_debug2 (D_ME_METHODS, "use_gfermi_scheme_rcl", &
         real(gf, kind=default))
    call use_gfermi_scheme_rcl (gf)
  end subroutine rclwrap_use_gfermi_scheme

  subroutine rclwrap_set_light_fermions (m)
    real(double), intent(in) :: m
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_light_fermions_rcl", &
         real(m, kind=default))
    call set_light_fermions_rcl (m)
  end subroutine rclwrap_set_light_fermions

  subroutine rclwrap_set_light_fermion (pdg_id)
    integer, intent(in) :: pdg_id
    if (debug_on) call msg_debug2 (D_ME_METHODS, "rclwrap_set_light_fermion", pdg_id)
    select case (abs(pdg_id))
    case (1)
       call set_light_down_rcl ()
    case (2)
       call set_light_up_rcl ()
    case (3)
       call set_light_strange_rcl ()
    case (4)
       call set_light_charm_rcl ()
    case (5)
       call set_light_bottom_rcl ()
    case (6)
       call set_light_top_rcl ()
    case (11)
       call set_light_electron_rcl ()
    case (13)
       call set_light_muon_rcl ()
    case (15)
       call set_light_tau_rcl ()
    end select
  end subroutine rclwrap_set_light_fermion

  subroutine rclwrap_unset_light_fermion (pdg_id)
    integer, intent(in) :: pdg_id
    if (debug_on) call msg_debug2 (D_ME_METHODS, "rclwrap_unset_light_fermion", pdg_id)
    select case (abs(pdg_id))
    case (1)
       call unset_light_down_rcl ()
    case (2)
       call unset_light_up_rcl ()
    case (3)
       call unset_light_strange_rcl ()
    case (4)
       call unset_light_charm_rcl ()
    case (5)
       call unset_light_bottom_rcl ()
    case (6)
       call unset_light_top_rcl ()
    case (11)
       call unset_light_electron_rcl ()
    case (13)
       call unset_light_muon_rcl ()
    case (15)
       call unset_light_tau_rcl ()
    end select
  end subroutine rclwrap_unset_light_fermion

  subroutine rclwrap_set_onshell_scheme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_on_shell_scheme_rcl")
    call set_on_shell_scheme_rcl ()
  end subroutine rclwrap_set_onshell_scheme

  subroutine rclwrap_set_alpha_s (alpha_s, mu, nf)
    real(double), intent(in) :: alpha_s, mu
    integer, intent(in) :: nf
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_alphas_rcl")
    call set_alphas_rcl (alpha_s, mu, nf)
  end subroutine rclwrap_set_alpha_s

  function rclwrap_get_alpha_s () result (alpha_s)
    real(double) :: alpha_s
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_alphas_rcl")
    call get_alphas_rcl (alpha_s)
  end function rclwrap_get_alpha_s

  function rclwrap_get_alpha () result (alpha)
    real(double) :: alpha
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_alpha_rcl")
    call get_alpha_rcl (alpha)
  end function rclwrap_get_alpha

  subroutine rclwrap_get_helicity_configurations (id, hel)
    integer, intent(in) :: id
    integer, intent(inout), dimension(:,:), allocatable :: hel
    call get_helicity_configurations_rcl (id, hel)
  end subroutine rclwrap_get_helicity_configurations

  subroutine rclwrap_get_color_configurations (id, col)
    integer, intent(in) :: id
    integer, intent(out), dimension(:,:), allocatable :: col
    call get_colour_configurations_rcl (id, col)
  end subroutine rclwrap_get_color_configurations

  subroutine rclwrap_use_dim_reg_soft ()
    if (debug_on) call msg_debug2 (D_ME_METHODS, "use_dim_reg_soft_rcl")
    call use_dim_reg_soft_rcl ()
  end subroutine rclwrap_use_dim_reg_soft

  subroutine rclwrap_use_mass_reg_soft (m)
    real(double), intent(in) :: m
    if (debug_on) call msg_debug2 (D_ME_METHODS, "use_mass_reg_soft_rcl")
    call use_mass_reg_soft_rcl (m)
  end subroutine rclwrap_use_mass_reg_soft

  subroutine rclwrap_set_delta_uv (d)
    real(double), intent(in) :: d
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_delta_uv_rcl")
    call set_delta_uv_rcl (d)
  end subroutine rclwrap_set_delta_uv

  subroutine rclwrap_set_mu_uv (mu)
    real(double), intent(in) :: mu
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_mu_uv_rcl")
    call set_mu_uv_rcl (mu)
  end subroutine rclwrap_set_mu_uv

  subroutine rclwrap_set_delta_ir (d, d2)
    real(double), intent(in) :: d, d2
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_delta_ir_rcl", &
         real(d, kind=default))
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_delta_ir_rcl", &
         real(d2, kind=default))
    call set_delta_ir_rcl (d, d2)
  end subroutine rclwrap_set_delta_ir

  subroutine rclwrap_set_mu_ir (mu)
    real(double), intent(in) :: mu
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_mu_ir_rcl")
    call set_mu_ir_rcl (mu)
  end subroutine rclwrap_set_mu_ir

  subroutine rclwrap_get_renormalization_scale (mu)
    real(double), intent(out) :: mu
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_renormalization_scale_rcl")
    call get_renormalization_scale_rcl (mu)
  end subroutine rclwrap_get_renormalization_scale

  subroutine rclwrap_get_flavor_scheme (nf)
    integer, intent(out) :: nf
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_flavour_scheme_rcl")
    call get_flavour_scheme_rcl (nf)
  end subroutine rclwrap_get_flavor_scheme

  subroutine rclwrap_use_alpha0_scheme (al0)
    real(double), intent(in), optional :: al0
    if (debug_on) call msg_debug2 (D_ME_METHODS, "use_alpha0_scheme_rcl")
    call use_alpha0_scheme_rcl (al0)
  end subroutine rclwrap_use_alpha0_scheme

  subroutine rclwrap_use_alphaz_scheme (alz)
    real(double), intent(in), optional :: alz
    if (debug_on) call msg_debug2 (D_ME_METHODS, "use_alphaz_scheme_rcl")
    call use_alphaz_scheme_rcl (alz)
  end subroutine rclwrap_use_alphaz_scheme

  subroutine rclwrap_set_complex_mass_scheme ()
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_complex_mass_scheme_rcl")
    call set_complex_mass_scheme_rcl ()
  end subroutine rclwrap_set_complex_mass_scheme

  subroutine rclwrap_set_resonant_particle (pdg_id)
    integer, intent(in) :: pdg_id
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_resonant_particle_rcl")
    call set_resonant_particle_rcl (char(get_recola_particle_string (pdg_id)))
  end subroutine rclwrap_set_resonant_particle

  subroutine rclwrap_switch_on_resonant_self_energies ()
    if (debug_on) call msg_debug2 (D_ME_METHODS, "switchon_resonant_selfenergies_rcl")
    call switchon_resonant_selfenergies_rcl ()
  end subroutine rclwrap_switch_on_resonant_self_energies

  subroutine rclwrap_switch_off_resonant_self_energies ()
    if (debug_on) call msg_debug2 (D_ME_METHODS, "switchoff_resonant_selfenergies_rcl")
    call switchoff_resonant_selfenergies_rcl ()
  end subroutine rclwrap_switch_off_resonant_self_energies

  subroutine rclwrap_set_draw_level_branches (n)
    integer, intent(in) :: n
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_draw_level_branches_rcl")
    call set_draw_level_branches_rcl (n)
  end subroutine rclwrap_set_draw_level_branches

  subroutine rclwrap_set_print_level_amplitude (n)
    integer, intent(in) :: n
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_print_level_amplitude_rcl")
    call set_print_level_amplitude_rcl (n)
  end subroutine rclwrap_set_print_level_amplitude

  subroutine rclwrap_set_print_level_squared_amplitude (n)
    integer, intent(in) :: n
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_print_level_squared_amplitude_rcl")
    call set_print_level_squared_amplitude_rcl (n)
  end subroutine rclwrap_set_print_level_squared_amplitude

  subroutine rclwrap_set_print_level_correlations (n)
    integer, intent(in) :: n
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_print_level_correlations_rcl")
    call set_print_level_correlations_rcl (n)
  end subroutine rclwrap_set_print_level_correlations

  subroutine rclwrap_set_print_level_RAM (n)
    integer, intent(in) :: n
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_print_level_RAM_rcl")
    call set_print_level_RAM_rcl (n)
  end subroutine rclwrap_set_print_level_RAM

  subroutine rclwrap_scale_coupling3 (pdg_id1, pdg_id2, pdg_id3, factor)
    integer, intent(in) :: pdg_id1, pdg_id2, pdg_id3
    complex(double), intent(in) :: factor
    if (debug_on) call msg_debug2 (D_ME_METHODS, "scale_coupling3_rcl")
    call scale_coupling3_rcl (factor, char(get_recola_particle_string (pdg_id1)), &
         char(get_recola_particle_string (pdg_id2)), char(get_recola_particle_string (pdg_id3)))
  end subroutine rclwrap_scale_coupling3

  subroutine rclwrap_scale_coupling4 (pdg_id1, pdg_id2, pdg_id3, pdg_id4, factor)
    integer, intent(in) :: pdg_id1, pdg_id2, pdg_id3, pdg_id4
    complex(double), intent(in) :: factor
    if (debug_on) call msg_debug2 (D_ME_METHODS, "scale_coupling4_rcl")
    call scale_coupling4_rcl (factor, char(get_recola_particle_string (pdg_id1)), &
         char(get_recola_particle_string (pdg_id2)), char(get_recola_particle_string (pdg_id3)), &
         char(get_recola_particle_string (pdg_id4)))
  end subroutine rclwrap_scale_coupling4

  subroutine rclwrap_switch_off_coupling3 (pdg_id1, pdg_id2, pdg_id3)
    integer, intent(in) :: pdg_id1, pdg_id2, pdg_id3
    if (debug_on) call msg_debug2 (D_ME_METHODS, "switchoff_coupling3_rcl")
    call switchoff_coupling3_rcl (char(get_recola_particle_string (pdg_id1)), &
         char(get_recola_particle_string (pdg_id2)), char(get_recola_particle_string (pdg_id3)))
  end subroutine rclwrap_switch_off_coupling3

  subroutine rclwrap_switch_off_coupling4 (pdg_id1, pdg_id2, pdg_id3, pdg_id4)
    integer, intent(in) :: pdg_id1, pdg_id2, pdg_id3, pdg_id4
    if (debug_on) call msg_debug2 (D_ME_METHODS, "switchoff_coupling4_rcl")
    call switchoff_coupling4_rcl &
         (char(get_recola_particle_string (pdg_id1)), &
          char(get_recola_particle_string (pdg_id2)), &
          char(get_recola_particle_string (pdg_id3)), &
          char(get_recola_particle_string (pdg_id4)))
  end subroutine rclwrap_switch_off_coupling4

  subroutine rclwrap_set_ifail (i)
    integer, intent(in) :: i
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_ifail_rcl")
    call set_ifail_rcl (i)
  end subroutine rclwrap_set_ifail

  subroutine rclwrap_get_ifail (i)
    integer, intent(out) :: i
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_ifail_rcl")
    call get_ifail_rcl (i)
  end subroutine rclwrap_get_ifail

  subroutine rclwrap_set_output_file (filename)
    character(len=*), intent(in) :: filename
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_output_file_rcl")
    call set_output_file_rcl (filename)
  end subroutine rclwrap_set_output_file

  subroutine rclwrap_set_gs_power (id, gs_array)
    integer, intent(in) :: id
    integer, dimension(:,:), intent(in) :: gs_array
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_gs_power_rcl")
    call set_gs_power_rcl (id, gs_array)
  end subroutine rclwrap_set_gs_power

  subroutine rclwrap_select_gs_power_born_amp (id, gs_power)
    integer, intent(in) :: id, gs_power
    if (debug_on) call msg_debug2 (D_ME_METHODS, "select_gs_power_BornAmpl_rcl")
    call select_gs_power_BornAmpl_rcl (id, gs_power)
 end subroutine rclwrap_select_gs_power_born_amp

  subroutine rclwrap_unselect_gs_power_born_amp (id, gs_power)
    integer, intent(in) :: id, gs_power
    if (debug_on) call msg_debug2 (D_ME_METHODS, "unselect_gs_power_BornAmpl_rcl")
    call unselect_gs_power_BornAmpl_rcl (id, gs_power)
 end subroutine rclwrap_unselect_gs_power_born_amp

  subroutine rclwrap_select_gs_power_loop_amp (id, gs_power)
    integer, intent(in) :: id, gs_power
    if (debug_on) call msg_debug2 (D_ME_METHODS, "select_gs_power_LoopAmpl_rcl")
    call select_gs_power_LoopAmpl_rcl (id, gs_power)
 end subroutine rclwrap_select_gs_power_loop_amp

  subroutine rclwrap_unselect_gs_power_loop_amp (id, gs_power)
    integer, intent(in) :: id, gs_power
    if (debug_on) call msg_debug2 (D_ME_METHODS, "unselect_gs_power_LoopAmpl_rcl")
    call unselect_gs_power_LoopAmpl_rcl (id, gs_power)
 end subroutine rclwrap_unselect_gs_power_loop_amp

  subroutine rclwrap_select_all_gs_powers_born_amp (id)
    integer, intent(in) :: id
    if (debug_on) call msg_debug2 (D_ME_METHODS, "select_all_gs_powers_BornAmpl_rcl")
    call select_all_gs_powers_BornAmpl_rcl (id)
  end subroutine rclwrap_select_all_gs_powers_born_amp

  subroutine rclwrap_unselect_all_gs_powers_loop_amp (id)
    integer, intent(in) :: id
    if (debug_on) call msg_debug2 (D_ME_METHODS, "unselect_all_gs_powers_BornAmpl_rcl")
    call unselect_all_gs_powers_BornAmpl_rcl (id)
  end subroutine rclwrap_unselect_all_gs_powers_loop_amp

  subroutine rclwrap_select_all_gs_powers_loop_amp (id)
    integer, intent(in) :: id
    if (debug_on) call msg_debug2 (D_ME_METHODS, "select_all_gs_powers_LoopAmpl_rcl")
    call select_all_gs_powers_LoopAmpl_rcl (id)
  end subroutine rclwrap_select_all_gs_powers_loop_amp

  subroutine rclwrap_unselect_all_gs_powers_born_amp (id)
    integer, intent(in) :: id
    if (debug_on) call msg_debug2 (D_ME_METHODS, "unselect_all_gs_powers_LoopAmpl_rcl")
    call unselect_all_gs_powers_LoopAmpl_rcl (id)
  end subroutine rclwrap_unselect_all_gs_powers_born_amp

  subroutine rclwrap_set_resonant_squared_momentum (id, i_res, p2)
    integer, intent(in) :: id, i_res
    real(double), intent(in) :: p2
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_resonant_squared_momentum_rcl")
    call set_resonant_squared_momentum_rcl (id, i_res, p2)
  end subroutine rclwrap_set_resonant_squared_momentum

  subroutine rclwrap_compute_running_alpha_s (Q, nf, n_loops)
    real(double), intent(in) :: Q
    integer, intent(in) :: nf, n_loops
    if (debug_on) call msg_debug2 (D_ME_METHODS, "compute_running_alphas_rcl")
    call compute_running_alphas_rcl (Q, nf, n_loops)
  end subroutine rclwrap_compute_running_alpha_s

  subroutine rclwrap_set_dynamic_settings ()
    if (debug_on) call msg_debug2 (D_ME_METHODS, "set_dynamic_settings_rcl")
    call set_dynamic_settings_rcl (1)
  end subroutine rclwrap_set_dynamic_settings

  subroutine rclwrap_rescale_process (id, order, sqme)
    integer, intent(in) :: id
    character(len=*), intent(in) :: order
    real(double), dimension(0:1), intent(out), optional :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "rescale_process_rcl")
    call rescale_process_rcl (id, order, sqme)
  end subroutine rclwrap_rescale_process

  subroutine rclwrap_get_polarized_squared_amplitude (id, &
     alphas_power, order, hel, sqme)
    integer, intent(in) :: id, alphas_power
    character(len=*), intent(in) :: order
    integer, dimension(:), intent(in) :: hel
    real(double), intent(out) :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_polarized_squared_amplitude_rcl")
    call get_polarized_squared_amplitude_rcl (id, alphas_power, &
         order, hel, sqme)
  end subroutine rclwrap_get_polarized_squared_amplitude

  subroutine rclwrap_compute_color_correlation (id, p, &
     i1, i2, sqme)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(in) :: p
    integer, intent(in) :: i1, i2
    real(double), intent(out), optional :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "compute_colour_correlation_rcl")
    call compute_colour_correlation_rcl (id, p, i1, i2, sqme)
  end subroutine rclwrap_compute_color_correlation

  subroutine rclwrap_compute_all_color_correlations (id, p)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(in) :: p
    if (debug_on) call msg_debug2 (D_ME_METHODS, "compute_all_colour_correlations_rcl")
    call compute_all_colour_correlations_rcl (id, p)
  end subroutine rclwrap_compute_all_color_correlations

  subroutine rclwrap_rescale_color_correlation (id, i1, i2, sqme)
    integer, intent(in) :: id, i1, i2
    real(double), intent(out), optional :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "rescale_colour_correlation_rcl")
    call rescale_colour_correlation_rcl (id, i1, i2, sqme)
  end subroutine rclwrap_rescale_color_correlation

  subroutine rclwrap_rescale_all_color_correlations (id)
    integer, intent(in) :: id
    if (debug_on) call msg_debug2 (D_ME_METHODS, "rescale_all_colour_correlations_rcl")
    call rescale_all_colour_correlations_rcl (id)
  end subroutine rclwrap_rescale_all_color_correlations

  subroutine rclwrap_get_color_correlation (id, alphas_power, i1, i2, sqme)
    integer, intent(in) :: id, alphas_power, i1, i2
    real(double), intent(out) :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_colour_correlation_rcl")
    call get_colour_correlation_rcl (id, alphas_power, i1, i2, sqme)
  end subroutine rclwrap_get_color_correlation

  subroutine rclwrap_compute_spin_correlation (id, p, i_photon, pol, sqme)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(in) :: p
    integer, intent(in) :: i_photon
    complex(double), dimension(:), intent(in) :: pol
    real(double), intent(out), optional :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "compute_spin_correlation_rcl")
    call compute_spin_correlation_rcl (id, p, i_photon, pol, sqme)
  end subroutine rclwrap_compute_spin_correlation

  subroutine rclwrap_rescale_spin_correlation (id, i_photon, pol, sqme)
    integer, intent(in) :: id, i_photon
    complex(double), dimension(:), intent(in) :: pol
    real(double), intent(out), optional :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "rescale_spin_correlation_rcl")
    call rescale_spin_correlation_rcl (id, i_photon, pol, sqme)
  end subroutine rclwrap_rescale_spin_correlation

  subroutine rclwrap_get_spin_correlation (id, alphas_power, sqme)
    integer, intent(in) :: id, alphas_power
    real(double), intent(out) :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_spin_correlation_rcl")
    call get_spin_correlation_rcl (id, alphas_power, sqme)
  end subroutine rclwrap_get_spin_correlation

  subroutine rclwrap_compute_spin_color_correlation (id, p, &
       i_gluon, i_spectator, pol, sqme)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(in) :: p
    integer, intent(in) :: i_gluon, i_spectator
    complex(double), dimension(:), intent(in) :: pol
    real(double), intent(out), optional :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "compute_spin_colour_correlation_rcl")
    call compute_spin_colour_correlation_rcl (id, p, &
         i_gluon, i_spectator, pol, sqme)
  end subroutine rclwrap_compute_spin_color_correlation

  subroutine rclwrap_rescale_spin_color_correlation (id, i_gluon, &
       i_spectator, pol, sqme)
    integer, intent(in) :: id, i_gluon, i_spectator
    complex(double), dimension(:), intent(in) :: pol
    real(double), intent(out), optional :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "rescale_spin_colour_correlation_rcl")
    call rescale_spin_colour_correlation_rcl (id, i_gluon, &
         i_spectator, pol, sqme)
  end subroutine rclwrap_rescale_spin_color_correlation

  subroutine rclwrap_get_spin_color_correlation (id, alphas_power, &
       i_gluon, i_spectator, sqme)
    integer, intent(in) :: id, alphas_power, i_gluon, i_spectator
    real(double), intent(out) :: sqme
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_spin_colour_correlation_rcl")
    call get_spin_colour_correlation_rcl (id, alphas_power, &
         i_gluon, i_spectator, sqme)
  end subroutine rclwrap_get_spin_color_correlation

  subroutine rclwrap_get_momenta (id, p)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(out) :: p
    if (debug_on) call msg_debug2 (D_ME_METHODS, "get_momenta_rcl")
    call get_momenta_rcl (id, p)
  end subroutine rclwrap_get_momenta

  subroutine rclwrap_reset_recola
    if (debug_on) call msg_debug (D_ME_METHODS, "rclwrap_reset_recola")
    call rcl_controller%reset ()
  end subroutine rclwrap_reset_recola


end module recola_wrapper
