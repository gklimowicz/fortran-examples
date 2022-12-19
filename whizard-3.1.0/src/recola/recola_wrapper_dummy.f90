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

  use kinds
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: rclwrap_is_active
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

  logical, parameter :: rclwrap_is_active = .false.


contains

  elemental function get_recola_particle_string (pdg) result (name)
    type(string_t) :: name
    integer, intent(in) :: pdg
    name = var_str ("?")
  end function get_recola_particle_string

  subroutine rclwrap_get_new_recola_id (id)
    integer, intent(out) :: id
    id = 0
  end subroutine rclwrap_get_new_recola_id

  function rclwrap_get_current_recola_id () result (n)
    integer :: n
    n = 0
  end function rclwrap_get_current_recola_id

  subroutine rclwrap_request_generate_processes ()
  end subroutine rclwrap_request_generate_processes

   subroutine rclwrap_add_process (id, process_string, order)
     integer, intent(in) :: id
     type(string_t), intent(in) :: process_string, order
   end subroutine rclwrap_add_process

   subroutine rclwrap_define_processes ()
   end subroutine rclwrap_define_processes

  subroutine rclwrap_generate_processes ()
  end subroutine rclwrap_generate_processes

  subroutine rclwrap_compute_process (id, p, order, sqme)
    integer, intent(in) :: id
    real(double), intent(in), dimension(:,:) :: p
    character(len=*), intent(in) :: order
    real(double), intent(out), dimension(0:1), optional :: sqme
  end subroutine rclwrap_compute_process

  subroutine rclwrap_get_amplitude (id, g_power, order, col, hel, amp)
    integer, intent(in) :: id, g_power
    character(len=*), intent(in) :: order
    integer, dimension(:), intent(in) :: col, hel
    complex(double), intent(out) :: amp
  end subroutine rclwrap_get_amplitude

  subroutine rclwrap_get_squared_amplitude (id, alphas_power, order, sqme)
    integer, intent(in) :: id, alphas_power
    character(len=*), intent(in) :: order
    real(double), intent(out) :: sqme
  end subroutine rclwrap_get_squared_amplitude

  subroutine rclwrap_set_pole_mass (pdg_id, mass, width)
    integer, intent(in) :: pdg_id
    real(double), intent(in) :: mass, width
  end subroutine rclwrap_set_pole_mass

  subroutine rclwrap_set_onshell_mass (pdg_id, mass, width)
    integer, intent(in) :: pdg_id
    real(double), intent(in) :: mass, width
  end subroutine rclwrap_set_onshell_mass

  subroutine rclwrap_use_gfermi_scheme (gf)
    real(double), intent(in), optional :: gf
  end subroutine rclwrap_use_gfermi_scheme

  subroutine rclwrap_set_light_fermions (m)
    real(double), intent(in) :: m
  end subroutine rclwrap_set_light_fermions

  subroutine rclwrap_set_light_fermion (pdg_id)
    integer, intent(in) :: pdg_id
  end subroutine rclwrap_set_light_fermion

  subroutine rclwrap_unset_light_fermion (pdg_id)
    integer, intent(in) :: pdg_id
  end subroutine rclwrap_unset_light_fermion

  subroutine rclwrap_set_onshell_scheme
  end subroutine rclwrap_set_onshell_scheme

  subroutine rclwrap_set_alpha_s (alpha_s, mu, nf)
    real(double), intent(in) :: alpha_s, mu
    integer, intent(in) :: nf
  end subroutine rclwrap_set_alpha_s

  function rclwrap_get_alpha_s () result (alpha_s)
    real(double) :: alpha_s
  end function rclwrap_get_alpha_s

  function rclwrap_get_alpha () result (alpha)
    real(double) :: alpha
  end function rclwrap_get_alpha

  subroutine rclwrap_get_helicity_configurations (id, hel)
    integer, intent(in) :: id
    integer, intent(inout), dimension(:,:), allocatable :: hel
  end subroutine rclwrap_get_helicity_configurations

  subroutine rclwrap_get_color_configurations (id, col)
    integer, intent(in) :: id
    integer, intent(out), dimension(:,:), allocatable :: col
  end subroutine rclwrap_get_color_configurations

  subroutine rclwrap_use_dim_reg_soft ()
  end subroutine rclwrap_use_dim_reg_soft

  subroutine rclwrap_use_mass_reg_soft (m)
    real(double), intent(in) :: m
  end subroutine rclwrap_use_mass_reg_soft

  subroutine rclwrap_set_delta_uv (d)
    real(double), intent(in) :: d
  end subroutine rclwrap_set_delta_uv

  subroutine rclwrap_set_mu_uv (mu)
    real(double), intent(in) :: mu
  end subroutine rclwrap_set_mu_uv

  subroutine rclwrap_set_delta_ir (d, d2)
    real(double), intent(in) :: d, d2
  end subroutine rclwrap_set_delta_ir

  subroutine rclwrap_set_mu_ir (mu)
    real(double), intent(in) :: mu
  end subroutine rclwrap_set_mu_ir

  subroutine rclwrap_get_renormalization_scale (mu)
    real(double), intent(out) :: mu
  end subroutine rclwrap_get_renormalization_scale

  subroutine rclwrap_get_flavor_scheme (nf)
    integer, intent(out) :: nf
  end subroutine rclwrap_get_flavor_scheme

  subroutine rclwrap_use_alpha0_scheme (al0)
    real(double), intent(in), optional :: al0
  end subroutine rclwrap_use_alpha0_scheme

  subroutine rclwrap_use_alphaz_scheme (alz)
    real(double), intent(in), optional :: alz
  end subroutine rclwrap_use_alphaz_scheme

  subroutine rclwrap_set_complex_mass_scheme ()
  end subroutine rclwrap_set_complex_mass_scheme

  subroutine rclwrap_set_resonant_particle (pdg_id)
    integer, intent(in) :: pdg_id
  end subroutine rclwrap_set_resonant_particle

  subroutine rclwrap_switch_on_resonant_self_energies ()
  end subroutine rclwrap_switch_on_resonant_self_energies

  subroutine rclwrap_switch_off_resonant_self_energies ()
  end subroutine rclwrap_switch_off_resonant_self_energies

  subroutine rclwrap_set_draw_level_branches (n)
    integer, intent(in) :: n
  end subroutine rclwrap_set_draw_level_branches

  subroutine rclwrap_set_print_level_amplitude (n)
    integer, intent(in) :: n
  end subroutine rclwrap_set_print_level_amplitude

  subroutine rclwrap_set_print_level_squared_amplitude (n)
    integer, intent(in) :: n
  end subroutine rclwrap_set_print_level_squared_amplitude

  subroutine rclwrap_set_print_level_correlations (n)
    integer, intent(in) :: n
  end subroutine rclwrap_set_print_level_correlations

  subroutine rclwrap_set_print_level_RAM (n)
    integer, intent(in) :: n
  end subroutine rclwrap_set_print_level_RAM

  subroutine rclwrap_scale_coupling3 (pdg_id1, pdg_id2, pdg_id3, factor)
    integer, intent(in) :: pdg_id1, pdg_id2, pdg_id3
    complex(double), intent(in) :: factor
  end subroutine rclwrap_scale_coupling3

  subroutine rclwrap_scale_coupling4 (pdg_id1, pdg_id2, pdg_id3, pdg_id4, factor)
    integer, intent(in) :: pdg_id1, pdg_id2, pdg_id3, pdg_id4
    complex(double), intent(in) :: factor
  end subroutine rclwrap_scale_coupling4

  subroutine rclwrap_switch_off_coupling3 (pdg_id1, pdg_id2, pdg_id3)
    integer, intent(in) :: pdg_id1, pdg_id2, pdg_id3
  end subroutine rclwrap_switch_off_coupling3

  subroutine rclwrap_switch_off_coupling4 (pdg_id1, pdg_id2, pdg_id3, pdg_id4)
    integer, intent(in) :: pdg_id1, pdg_id2, pdg_id3, pdg_id4
  end subroutine rclwrap_switch_off_coupling4

  subroutine rclwrap_set_ifail (i)
    integer, intent(in) :: i
  end subroutine rclwrap_set_ifail

  subroutine rclwrap_get_ifail (i)
    integer, intent(out) :: i
  end subroutine rclwrap_get_ifail

  subroutine rclwrap_set_output_file (filename)
    character(len=*), intent(in) :: filename
  end subroutine rclwrap_set_output_file

  subroutine rclwrap_set_gs_power (id, gs_array)
    integer, intent(in) :: id
    integer, dimension(:,:), intent(in) :: gs_array
  end subroutine rclwrap_set_gs_power

  subroutine rclwrap_select_gs_power_born_amp (id, gs_power)
    integer, intent(in) :: id, gs_power
 end subroutine rclwrap_select_gs_power_born_amp

  subroutine rclwrap_unselect_gs_power_born_amp (id, gs_power)
    integer, intent(in) :: id, gs_power
 end subroutine rclwrap_unselect_gs_power_born_amp

  subroutine rclwrap_select_gs_power_loop_amp (id, gs_power)
    integer, intent(in) :: id, gs_power
 end subroutine rclwrap_select_gs_power_loop_amp

  subroutine rclwrap_unselect_gs_power_loop_amp (id, gs_power)
    integer, intent(in) :: id, gs_power
 end subroutine rclwrap_unselect_gs_power_loop_amp

  subroutine rclwrap_select_all_gs_powers_born_amp (id)
    integer, intent(in) :: id
  end subroutine rclwrap_select_all_gs_powers_born_amp

  subroutine rclwrap_unselect_all_gs_powers_loop_amp (id)
    integer, intent(in) :: id
  end subroutine rclwrap_unselect_all_gs_powers_loop_amp

  subroutine rclwrap_select_all_gs_powers_loop_amp (id)
    integer, intent(in) :: id
  end subroutine rclwrap_select_all_gs_powers_loop_amp

  subroutine rclwrap_unselect_all_gs_powers_born_amp (id)
    integer, intent(in) :: id
  end subroutine rclwrap_unselect_all_gs_powers_born_amp

  subroutine rclwrap_set_resonant_squared_momentum (id, i_res, p2)
    integer, intent(in) :: id, i_res
    real(double), intent(in) :: p2
  end subroutine rclwrap_set_resonant_squared_momentum

  subroutine rclwrap_compute_running_alpha_s (Q, nf, n_loops)
    real(double), intent(in) :: Q
    integer, intent(in) :: nf, n_loops
  end subroutine rclwrap_compute_running_alpha_s

  subroutine rclwrap_set_dynamic_settings ()
  end subroutine rclwrap_set_dynamic_settings

  subroutine rclwrap_rescale_process (id, order, sqme)
    integer, intent(in) :: id
    character(len=*), intent(in) :: order
    real(double), dimension(0:1), intent(out), optional :: sqme
  end subroutine rclwrap_rescale_process

  subroutine rclwrap_get_polarized_squared_amplitude (id, &
     alphas_power, order, hel, sqme)
    integer, intent(in) :: id, alphas_power
    character(len=*), intent(in) :: order
    integer, dimension(:), intent(in) :: hel
    real(double), intent(out) :: sqme
  end subroutine rclwrap_get_polarized_squared_amplitude

  subroutine rclwrap_compute_color_correlation (id, p, &
     i1, i2, sqme)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(in) :: p
    integer, intent(in) :: i1, i2
    real(double), intent(out), optional :: sqme
  end subroutine rclwrap_compute_color_correlation

  subroutine rclwrap_compute_all_color_correlations (id, p)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(in) :: p
  end subroutine rclwrap_compute_all_color_correlations

  subroutine rclwrap_rescale_color_correlation (id, i1, i2, sqme)
    integer, intent(in) :: id, i1, i2
    real(double), intent(out), optional :: sqme
  end subroutine rclwrap_rescale_color_correlation

  subroutine rclwrap_rescale_all_color_correlations (id)
    integer, intent(in) :: id
  end subroutine rclwrap_rescale_all_color_correlations

  subroutine rclwrap_get_color_correlation (id, alphas_power, i1, i2, sqme)
    integer, intent(in) :: id, alphas_power, i1, i2
    real(double), intent(out) :: sqme
  end subroutine rclwrap_get_color_correlation

  subroutine rclwrap_compute_spin_correlation (id, p, i_photon, pol, sqme)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(in) :: p
    integer, intent(in) :: i_photon
    complex(double), dimension(:), intent(in) :: pol
    real(double), intent(out), optional :: sqme
  end subroutine rclwrap_compute_spin_correlation

  subroutine rclwrap_rescale_spin_correlation (id, i_photon, pol, sqme)
    integer, intent(in) :: id, i_photon
    complex(double), dimension(:), intent(in) :: pol
    real(double), intent(out), optional :: sqme
  end subroutine rclwrap_rescale_spin_correlation

  subroutine rclwrap_get_spin_correlation (id, alphas_power, sqme)
    integer, intent(in) :: id, alphas_power
    real(double), intent(out) :: sqme
  end subroutine rclwrap_get_spin_correlation

  subroutine rclwrap_compute_spin_color_correlation (id, p, &
     i_gluon, i_spectator, pol, sqme)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(in) :: p
    integer, intent(in) :: i_gluon, i_spectator
    complex(double), dimension(:), intent(in) :: pol
    real(double), intent(out), optional :: sqme
  end subroutine rclwrap_compute_spin_color_correlation

  subroutine rclwrap_rescale_spin_color_correlation (id, i_gluon, &
     i_spectator, pol, sqme)
    integer, intent(in) :: id, i_gluon, i_spectator
    complex(double), dimension(:), intent(in) :: pol
    real(double), intent(out), optional :: sqme
  end subroutine rclwrap_rescale_spin_color_correlation

  subroutine rclwrap_get_spin_color_correlation (id, alphas_power, &
     i_gluon, i_spectator, sqme)
    integer, intent(in) :: id, alphas_power, i_gluon, i_spectator
    real(double), intent(out) :: sqme
  end subroutine rclwrap_get_spin_color_correlation

  subroutine rclwrap_get_momenta (id, p)
    integer, intent(in) :: id
    real(double), dimension(:,:), intent(out) :: p
  end subroutine rclwrap_get_momenta

  subroutine rclwrap_reset_recola
  end subroutine rclwrap_reset_recola


end module recola_wrapper
