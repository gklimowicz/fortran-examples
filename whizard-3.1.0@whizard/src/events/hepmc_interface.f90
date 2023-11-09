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

module hepmc_interface

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use flavors
  use colors
  use helicities
  use polarizations
  use event_handles, only: event_handle_t

  implicit none
  private

  public :: hepmc_is_available
  public :: hepmc_four_vector_t
  public :: hepmc_four_vector_init
  public :: hepmc_four_vector_final
  public :: hepmc_four_vector_to_vector4
  public :: hepmc_polarization_t
  public :: hepmc_polarization_init
  public :: hepmc_polarization_final
  public :: hepmc_polarization_to_pol
  public :: hepmc_polarization_to_hel
  public :: hepmc_particle_t
  public :: hepmc_particle_init
  public :: hepmc_particle_set_color
  public :: hepmc_particle_set_polarization
  public :: hepmc_particle_get_barcode
  public :: hepmc_particle_get_momentum
  public :: hepmc_particle_get_mass_squared
  public :: hepmc_particle_get_pdg
  public :: hepmc_particle_get_status
  public :: hepmc_particle_is_beam
  public :: hepmc_particle_get_production_vertex
  public :: hepmc_particle_get_decay_vertex
  public :: hepmc_particle_get_parent_barcodes
  public :: hepmc_particle_get_child_barcodes
  public :: hepmc_particle_get_polarization
  public :: hepmc_particle_get_color
  public :: hepmc_vertex_to_vertex
  public :: hepmc_vertex_t
  public :: hepmc_vertex_init
  public :: hepmc_vertex_is_valid
  public :: hepmc_vertex_add_particle_in
  public :: hepmc_vertex_add_particle_out
  public :: hepmc_vertex_get_n_in
  public :: hepmc_vertex_get_n_out
  public :: hepmc_particle_get_parents
  public :: hepmc_particle_get_children
  public :: hepmc_particle_get_n_parents
  public :: hepmc_particle_get_n_children
  public :: hepmc_vertex_particle_in_iterator_t
  public :: hepmc_vertex_particle_in_iterator_init
  public :: hepmc_vertex_particle_in_iterator_final
  public :: hepmc_vertex_particle_in_iterator_advance
  public :: hepmc_vertex_particle_in_iterator_reset
  public :: hepmc_vertex_particle_in_iterator_is_valid
  public :: hepmc_vertex_particle_in_iterator_get
  public :: hepmc_vertex_get_nth_particle_in
  public :: hepmc_vertex_get_nth_particle_out
  public :: hepmc_vertex_particle_out_iterator_t
  public :: hepmc_vertex_particle_out_iterator_init
  public :: hepmc_vertex_particle_out_iterator_final
  public :: hepmc_vertex_particle_out_iterator_advance
  public :: hepmc_vertex_particle_out_iterator_reset
  public :: hepmc_vertex_particle_out_iterator_is_valid
  public :: hepmc_vertex_particle_out_iterator_get
  public :: hepmc_event_t
  public :: hepmc_event_init
  public :: hepmc_event_final
  public :: hepmc_event_nullify
  public :: hepmc_event_get_c_ptr
  public :: hepmc_event_print
  public :: hepmc_event_get_event_index
  public :: hepmc_event_get_n_particles
  public :: hepmc_event_get_n_beams
  public :: hepmc_event_set_process_id
  public :: hepmc_event_get_process_id
  public :: hepmc_event_set_scale
  public :: hepmc_event_get_scale
  public :: hepmc_event_set_alpha_qcd
  public :: hepmc_event_get_alpha_qcd
  public :: hepmc_event_set_alpha_qed
  public :: hepmc_event_get_alpha_qed
  public :: hepmc_event_clear_weights
  public :: hepmc_event_add_weight
  public :: hepmc_event_get_weights_size
  public :: hepmc_event_get_weight
  public :: hepmc_event_add_vertex
  public :: hepmc_event_set_signal_process_vertex
  public :: hepmc_event_get_signal_process_vertex
  public :: hepmc_event_set_beam_particles
  public :: hepmc_event_set_cross_section
  public :: hepmc_event_particle_iterator_t
  public :: hepmc_event_particle_iterator_init
  public :: hepmc_event_particle_iterator_final
  public :: hepmc_event_particle_iterator_advance
  public :: hepmc_event_particle_iterator_reset
  public :: hepmc_event_particle_iterator_is_valid
  public :: hepmc_event_particle_iterator_get
  public :: hepmc_event_get_nth_particle
  public :: hepmc_event_get_nth_beam
  public :: hepmc_iostream_t
  public :: hepmc_iostream_open_out
  public :: hepmc_iostream_open_in
  public :: hepmc_iostream_close
  public :: hepmc_iostream_write_event
  public :: hepmc_iostream_read_event

  type :: hepmc_four_vector_t
     private
     type(c_ptr) :: obj
  end type hepmc_four_vector_t

  type :: hepmc_polarization_t
     private
     logical :: polarized = .false.
     type(c_ptr) :: obj
  end type hepmc_polarization_t

  type :: hepmc_particle_t
     private
     type(c_ptr) :: obj
  end type hepmc_particle_t

  type :: hepmc_vertex_t
     private
     type(c_ptr) :: obj
  end type hepmc_vertex_t

  type :: hepmc_vertex_particle_in_iterator_t
     private
     type(c_ptr) :: obj
     type(c_ptr) :: v_obj
  end type hepmc_vertex_particle_in_iterator_t

  type :: hepmc_vertex_particle_out_iterator_t
     private
     type(c_ptr) :: obj
     type(c_ptr) :: v_obj
  end type hepmc_vertex_particle_out_iterator_t

  type, extends (event_handle_t) :: hepmc_event_t
     private
     type(c_ptr) :: obj = c_null_ptr
  end type hepmc_event_t

  type :: hepmc_event_particle_iterator_t
     private
     type(c_ptr) :: obj
     type(c_ptr) :: evt_obj
  end type hepmc_event_particle_iterator_t

  type :: hepmc_iostream_t
     private
     type(c_ptr) :: obj
  end type hepmc_iostream_t


  integer, parameter, public :: HEPMC3_MODE_HEPMC2 = 1
  integer, parameter, public :: HEPMC3_MODE_HEPMC3 = 2
  integer, parameter, public :: HEPMC3_MODE_ROOT = 3
  integer, parameter, public :: HEPMC3_MODE_ROOTTREE = 4
  integer, parameter, public :: HEPMC3_MODE_HEPEVT = 5

  interface
     logical(c_bool) function hepmc_available () bind(C)
       import
     end function hepmc_available
  end interface
  interface
     type(c_ptr) function new_four_vector_xyz (x, y, z) bind(C)
       import
       real(c_double), value :: x, y, z
     end function new_four_vector_xyz
  end interface
  interface
     type(c_ptr) function new_four_vector_xyzt (x, y, z, t) bind(C)
       import
       real(c_double), value :: x, y, z, t
     end function new_four_vector_xyzt
  end interface
  interface hepmc_four_vector_init
     module procedure hepmc_four_vector_init_v4
     module procedure hepmc_four_vector_init_v3
     module procedure hepmc_four_vector_init_hepmc_prt
  end interface
  interface
     subroutine four_vector_delete (p_obj) bind(C)
       import
       type(c_ptr), value :: p_obj
     end subroutine four_vector_delete
  end interface
  interface
     function four_vector_px (p_obj) result (px) bind(C)
       import
       real(c_double) :: px
       type(c_ptr), value :: p_obj
     end function four_vector_px
  end interface
  interface
     function four_vector_py (p_obj) result (py) bind(C)
       import
       real(c_double) :: py
       type(c_ptr), value :: p_obj
     end function four_vector_py
  end interface
  interface
     function four_vector_pz (p_obj) result (pz) bind(C)
       import
       real(c_double) :: pz
       type(c_ptr), value :: p_obj
     end function four_vector_pz
  end interface
  interface
     function four_vector_e (p_obj) result (e) bind(C)
       import
       real(c_double) :: e
       type(c_ptr), value :: p_obj
     end function four_vector_e
  end interface
  interface
     type(c_ptr) function new_polarization (theta, phi) bind(C)
       import
       real(c_double), value :: theta, phi
     end function new_polarization
  end interface
  interface hepmc_polarization_init
     module procedure hepmc_polarization_init_pol
     module procedure hepmc_polarization_init_hel
     module procedure hepmc_polarization_init_int
  end interface
  interface
     subroutine polarization_delete (pol_obj) bind(C)
       import
       type(c_ptr), value :: pol_obj
     end subroutine polarization_delete
  end interface
  interface
     function polarization_theta (pol_obj) result (theta) bind(C)
       import
       real(c_double) :: theta
       type(c_ptr), value :: pol_obj
     end function polarization_theta
  end interface
  interface
     function polarization_phi (pol_obj) result (phi) bind(C)
       import
       real(c_double) :: phi
       type(c_ptr), value :: pol_obj
     end function polarization_phi
  end interface
  interface
     type(c_ptr) function new_gen_particle (prt_obj, pdg_id, status) bind(C)
       import
       type(c_ptr), value :: prt_obj
       integer(c_int), value :: pdg_id, status
     end function new_gen_particle
  end interface
  interface
     subroutine gen_particle_set_flow (prt_obj, code_index, code) bind(C)
       import
       type(c_ptr), value :: prt_obj
       integer(c_int), value :: code_index, code
     end subroutine gen_particle_set_flow
  end interface
  interface hepmc_particle_set_color
     module procedure hepmc_particle_set_color_col
     module procedure hepmc_particle_set_color_int
  end interface hepmc_particle_set_color
  interface
     subroutine gen_particle_set_polarization (prt_obj, pol_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj, pol_obj
     end subroutine gen_particle_set_polarization
  end interface
  interface hepmc_particle_set_polarization
     module procedure hepmc_particle_set_polarization_pol
     module procedure hepmc_particle_set_polarization_hel
     module procedure hepmc_particle_set_polarization_int
  end interface
  interface
     function gen_particle_barcode (prt_obj) result (barcode) bind(C)
       import
       integer(c_int) :: barcode
       type(c_ptr), value :: prt_obj
     end function gen_particle_barcode
  end interface
  interface
     type(c_ptr) function gen_particle_momentum (prt_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj
     end function gen_particle_momentum
  end interface
  interface
     function gen_particle_generated_mass (prt_obj) result (mass) bind(C)
       import
       real(c_double) :: mass
       type(c_ptr), value :: prt_obj
     end function gen_particle_generated_mass
  end interface
  interface
     function gen_particle_pdg_id (prt_obj) result (pdg_id) bind(C)
       import
       integer(c_int) :: pdg_id
       type(c_ptr), value :: prt_obj
     end function gen_particle_pdg_id
  end interface
  interface
     function gen_particle_status (prt_obj) result (status) bind(C)
       import
       integer(c_int) :: status
       type(c_ptr), value :: prt_obj
     end function gen_particle_status
  end interface
  interface
     function gen_particle_is_beam (prt_obj) result (is_beam) bind(C)
       import
       logical(c_bool) :: is_beam
       type(c_ptr), value :: prt_obj
     end function gen_particle_is_beam
  end interface
  interface
     type(c_ptr) function gen_particle_production_vertex (prt_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj
     end function gen_particle_production_vertex
  end interface
  interface
     type(c_ptr) function gen_particle_end_vertex (prt_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj
     end function gen_particle_end_vertex
  end interface
  interface
     type(c_ptr) function gen_particle_polarization (prt_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj
     end function gen_particle_polarization
  end interface
  interface
     function gen_particle_flow (prt_obj, code_index) result (code) bind(C)
       import
       integer(c_int) :: code
       type(c_ptr), value :: prt_obj
       integer(c_int), value :: code_index
     end function gen_particle_flow
  end interface
  interface
     function gen_vertex_pos_x (v_obj) result (x) bind(C)
       import
       type(c_ptr), value :: v_obj
       real(c_double) :: x
     end function gen_vertex_pos_x
  end interface
  interface
     function gen_vertex_pos_y (v_obj) result (y) bind(C)
       import
       type(c_ptr), value :: v_obj
       real(c_double) :: y
     end function gen_vertex_pos_y
  end interface
  interface
     function gen_vertex_pos_z (v_obj) result (z) bind(C)
       import
       type(c_ptr), value :: v_obj
       real(c_double) :: z
     end function gen_vertex_pos_z
  end interface
  interface
     function gen_vertex_time (v_obj) result (t) bind(C)
       import
       type(c_ptr), value :: v_obj
       real(c_double) :: t
     end function gen_vertex_time
  end interface
  interface
     type(c_ptr) function new_gen_vertex () bind(C)
       import
     end function new_gen_vertex
  end interface
  interface
     type(c_ptr) function new_gen_vertex_pos (prt_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj
     end function new_gen_vertex_pos
  end interface
  interface
     function gen_vertex_is_valid (v_obj) result (flag) bind(C)
       import
       logical(c_bool) :: flag
       type(c_ptr), value :: v_obj
     end function gen_vertex_is_valid
  end interface
  interface
     subroutine gen_vertex_add_particle_in (v_obj, prt_obj) bind(C)
       import
       type(c_ptr), value :: v_obj, prt_obj
     end subroutine gen_vertex_add_particle_in
  end interface
  interface
     subroutine gen_vertex_add_particle_out (v_obj, prt_obj) bind(C)
       import
       type(c_ptr), value :: v_obj, prt_obj
     end subroutine gen_vertex_add_particle_out
  end interface
  interface
     function gen_vertex_particles_in_size (v_obj) result (size) bind(C)
       import
       integer(c_int) :: size
       type(c_ptr), value :: v_obj
     end function gen_vertex_particles_in_size
  end interface
  interface
     function gen_vertex_particles_out_size (v_obj) result (size) bind(C)
       import
       integer(c_int) :: size
       type(c_ptr), value :: v_obj
     end function gen_vertex_particles_out_size
  end interface
  interface
     function gen_particle_get_n_parents (p_obj) result (size) bind(C)
       import
       integer(c_int) :: size
       type(c_ptr), value :: p_obj
     end function gen_particle_get_n_parents
  end interface
  interface
     function gen_particle_get_n_children (p_obj) result (size) bind(C)
       import
       integer(c_int) :: size
       type(c_ptr), value :: p_obj
     end function gen_particle_get_n_children
  end interface
  interface
     type(c_ptr) function &
          new_vertex_particles_in_const_iterator (v_obj) bind(C)
       import
       type(c_ptr), value :: v_obj
     end function new_vertex_particles_in_const_iterator
  end interface
  interface
     subroutine vertex_particles_in_const_iterator_delete (it_obj) bind(C)
       import
       type(c_ptr), value :: it_obj
     end subroutine vertex_particles_in_const_iterator_delete
  end interface
  interface
     subroutine vertex_particles_in_const_iterator_advance (it_obj) bind(C)
       import
       type(c_ptr), value :: it_obj
     end subroutine vertex_particles_in_const_iterator_advance
  end interface
  interface
     subroutine vertex_particles_in_const_iterator_reset &
          (it_obj, v_obj) bind(C)
       import
       type(c_ptr), value :: it_obj, v_obj
     end subroutine vertex_particles_in_const_iterator_reset
  end interface
  interface
     function vertex_particles_in_const_iterator_is_valid &
          (it_obj, v_obj) result (flag) bind(C)
       import
       logical(c_bool) :: flag
       type(c_ptr), value :: it_obj, v_obj
     end function vertex_particles_in_const_iterator_is_valid
  end interface
  interface
     type(c_ptr) function &
          vertex_particles_in_const_iterator_get (it_obj) bind(C)
       import
       type(c_ptr), value :: it_obj
     end function vertex_particles_in_const_iterator_get
  end interface
  interface
     type(c_ptr) function vertex_get_nth_particle_in (vtx_obj, n) bind(C)
       import
       type(c_ptr), value :: vtx_obj
       integer(c_int), value :: n
     end function vertex_get_nth_particle_in
  end interface
  interface
     type(c_ptr) function vertex_get_nth_particle_out (vtx_obj, n) bind(C)
       import
       type(c_ptr), value :: vtx_obj
       integer(c_int), value :: n
     end function vertex_get_nth_particle_out
  end interface
  interface
     type(c_ptr) function &
          new_vertex_particles_out_const_iterator (v_obj) bind(C)
       import
       type(c_ptr), value :: v_obj
     end function new_vertex_particles_out_const_iterator
  end interface
  interface
     subroutine vertex_particles_out_const_iterator_delete (it_obj) bind(C)
       import
       type(c_ptr), value :: it_obj
     end subroutine vertex_particles_out_const_iterator_delete
  end interface
  interface
     subroutine vertex_particles_out_const_iterator_advance (it_obj) bind(C)
       import
       type(c_ptr), value :: it_obj
     end subroutine vertex_particles_out_const_iterator_advance
  end interface
  interface
     subroutine vertex_particles_out_const_iterator_reset &
          (it_obj, v_obj) bind(C)
       import
       type(c_ptr), value :: it_obj, v_obj
     end subroutine vertex_particles_out_const_iterator_reset
  end interface
  interface
     function vertex_particles_out_const_iterator_is_valid &
          (it_obj, v_obj) result (flag) bind(C)
       import
       logical(c_bool) :: flag
       type(c_ptr), value :: it_obj, v_obj
     end function vertex_particles_out_const_iterator_is_valid
  end interface
  interface
     type(c_ptr) function &
          vertex_particles_out_const_iterator_get (it_obj) bind(C)
       import
       type(c_ptr), value :: it_obj
     end function vertex_particles_out_const_iterator_get
  end interface
  interface
     type(c_ptr) function new_gen_event (proc_id, event_id) bind(C)
       import
       integer(c_int), value :: proc_id, event_id
     end function new_gen_event
  end interface
  interface
     subroutine gen_event_delete (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end subroutine gen_event_delete
  end interface
  interface
     subroutine gen_event_print (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end subroutine gen_event_print
  end interface
  interface
     integer(c_int) function gen_event_event_number (evt_obj) bind(C)
       use iso_c_binding !NODEP!
       type(c_ptr), value :: evt_obj
     end function gen_event_event_number
  end interface
  interface
     integer(c_int) function gen_event_get_n_particles &
          (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function gen_event_get_n_particles
  end interface
    interface
     integer(c_int) function gen_event_get_n_beams &
          (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function gen_event_get_n_beams
  end interface
  interface
     subroutine gen_event_set_signal_process_id (evt_obj, proc_id) bind(C)
       import
       type(c_ptr), value :: evt_obj
       integer(c_int), value :: proc_id
     end subroutine gen_event_set_signal_process_id
  end interface
  interface
     integer(c_int) function gen_event_signal_process_id (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function gen_event_signal_process_id
  end interface
  interface
     subroutine gen_event_set_event_scale (evt_obj, scale) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: scale
     end subroutine gen_event_set_event_scale
  end interface
  interface
     real(c_double) function gen_event_event_scale (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function gen_event_event_scale
  end interface
  interface
     subroutine gen_event_set_alpha_qcd (evt_obj, a) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: a
     end subroutine gen_event_set_alpha_qcd
  end interface
  interface
     real(c_double) function gen_event_alpha_qcd (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function gen_event_alpha_qcd
  end interface
  interface
     subroutine gen_event_set_alpha_qed (evt_obj, a) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: a
     end subroutine gen_event_set_alpha_qed
  end interface
  interface
     real(c_double) function gen_event_alpha_qed (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function gen_event_alpha_qed
  end interface
  interface
     subroutine gen_event_clear_weights (evt_obj) bind(C)
       use iso_c_binding !NODEP!
       type(c_ptr), value :: evt_obj
     end subroutine gen_event_clear_weights
  end interface
  interface
     subroutine gen_event_add_weight (evt_obj, w) bind(C)
       use iso_c_binding !NODEP!
       type(c_ptr), value :: evt_obj
       real(c_double), value :: w
     end subroutine gen_event_add_weight
  end interface
  interface
     integer(c_int) function gen_event_weights_size (evt_obj) bind(C)
       use iso_c_binding !NODEP!
       type(c_ptr), value :: evt_obj
     end function gen_event_weights_size
  end interface
  interface
     real(c_double) function gen_event_weight (evt_obj, i) bind(C)
       use iso_c_binding !NODEP!
       type(c_ptr), value :: evt_obj
       integer(c_int), value :: i
     end function gen_event_weight
  end interface
  interface
     subroutine gen_event_add_vertex (evt_obj, v_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
       type(c_ptr), value :: v_obj
     end subroutine gen_event_add_vertex
  end interface
  interface
     subroutine gen_event_set_signal_process_vertex (evt_obj, v_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
       type(c_ptr), value :: v_obj
     end subroutine gen_event_set_signal_process_vertex
  end interface
  interface
     function gen_event_get_signal_process_vertex (evt_obj) &
          result (v_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
       type(c_ptr) :: v_obj
     end function gen_event_get_signal_process_vertex
  end interface
  interface
     logical(c_bool) function gen_event_set_beam_particles &
          (evt_obj, prt1_obj, prt2_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj, prt1_obj, prt2_obj
     end function gen_event_set_beam_particles
  end interface

  interface
     subroutine gen_event_set_cross_section (evt_obj, xs, xs_err) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: xs, xs_err
     end subroutine gen_event_set_cross_section
  end interface

  interface
     type(c_ptr) function new_event_particle_const_iterator (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function new_event_particle_const_iterator
  end interface
  interface
     subroutine event_particle_const_iterator_delete (it_obj) bind(C)
       import
       type(c_ptr), value :: it_obj
     end subroutine event_particle_const_iterator_delete
  end interface
  interface
     subroutine event_particle_const_iterator_advance (it_obj) bind(C)
       import
       type(c_ptr), value :: it_obj
     end subroutine event_particle_const_iterator_advance
  end interface
  interface
     subroutine event_particle_const_iterator_reset (it_obj, evt_obj) bind(C)
       import
       type(c_ptr), value :: it_obj, evt_obj
     end subroutine event_particle_const_iterator_reset
  end interface
  interface
     function event_particle_const_iterator_is_valid &
          (it_obj, evt_obj) result (flag) bind(C)
       import
       logical(c_bool) :: flag
       type(c_ptr), value :: it_obj, evt_obj
     end function event_particle_const_iterator_is_valid
  end interface
  interface
     type(c_ptr) function event_particle_const_iterator_get (it_obj) bind(C)
       import
       type(c_ptr), value :: it_obj
     end function event_particle_const_iterator_get
  end interface
  interface
     type(c_ptr) function gen_event_get_nth_particle (evt_obj, n) bind(C)
       import
       type(c_ptr), value :: evt_obj
       integer(c_int), value :: n
     end function gen_event_get_nth_particle
  end interface
    interface
     integer(c_int) function gen_event_get_nth_beam (evt_obj, n) bind(C)
       import
       type(c_ptr), value :: evt_obj
       integer(c_int), value :: n
     end function gen_event_get_nth_beam
  end interface
  interface
     type(c_ptr) function new_io_gen_event_out (hepmc3_mode, filename) bind(C)
       import
       integer(c_int), intent(in) :: hepmc3_mode
       character(c_char), dimension(*), intent(in) :: filename
     end function new_io_gen_event_out
  end interface
  interface
     type(c_ptr) function new_io_gen_event_in (hepmc3_mode, filename) bind(C)
       import
       integer(c_int), intent(in) :: hepmc3_mode
       character(c_char), dimension(*), intent(in) :: filename
     end function new_io_gen_event_in
  end interface
  interface
     subroutine io_gen_event_delete (io_obj) bind(C)
       import
       type(c_ptr), value :: io_obj
     end subroutine io_gen_event_delete
  end interface
  interface
     subroutine io_gen_event_write_event (io_obj, evt_obj) bind(C)
       import
       type(c_ptr), value :: io_obj, evt_obj
     end subroutine io_gen_event_write_event
  end interface
  interface
     subroutine io_gen_event_write_event_hepmc2 (io_obj, evt_obj) bind(C)
       import
       type(c_ptr), value :: io_obj, evt_obj
     end subroutine io_gen_event_write_event_hepmc2
  end interface
  interface
     logical(c_bool) function io_gen_event_read_event (io_obj, evt_obj) bind(C)
       import
       type(c_ptr), value :: io_obj, evt_obj
     end function io_gen_event_read_event
  end interface

  interface
    module function hepmc_is_available () result (flag)
      logical :: flag
    end function hepmc_is_available
    module subroutine hepmc_four_vector_init_v4 (pp, p)
      type(hepmc_four_vector_t), intent(out) :: pp
      type(vector4_t), intent(in) :: p
    end subroutine hepmc_four_vector_init_v4
    module subroutine hepmc_four_vector_init_v3 (pp, p)
      type(hepmc_four_vector_t), intent(out) :: pp
      type(vector3_t), intent(in) :: p
    end subroutine hepmc_four_vector_init_v3
    module subroutine hepmc_four_vector_init_hepmc_prt (pp, prt)
      type(hepmc_four_vector_t), intent(out) :: pp
      type(hepmc_particle_t), intent(in) :: prt
    end subroutine hepmc_four_vector_init_hepmc_prt
    module subroutine hepmc_four_vector_final (p)
      type(hepmc_four_vector_t), intent(inout) :: p
    end subroutine hepmc_four_vector_final
    module subroutine hepmc_four_vector_to_vector4 (pp, p)
      type(hepmc_four_vector_t), intent(in) :: pp
      type(vector4_t), intent(out) :: p
    end subroutine hepmc_four_vector_to_vector4
    module subroutine hepmc_polarization_init_pol (hpol, pol)
      type(hepmc_polarization_t), intent(out) :: hpol
      type(polarization_t), intent(in) :: pol
    end subroutine hepmc_polarization_init_pol
    module subroutine hepmc_polarization_init_hel (hpol, hel)
      type(hepmc_polarization_t), intent(out) :: hpol
      type(helicity_t), intent(in) :: hel
    end subroutine hepmc_polarization_init_hel
    module subroutine hepmc_polarization_init_int (hpol, hel)
      type(hepmc_polarization_t), intent(out) :: hpol
      integer, intent(in) :: hel
    end subroutine hepmc_polarization_init_int
    module subroutine hepmc_polarization_final (hpol)
      type(hepmc_polarization_t), intent(inout) :: hpol
    end subroutine hepmc_polarization_final
    module subroutine hepmc_polarization_to_pol (hpol, flv, pol)
      type(hepmc_polarization_t), intent(in) :: hpol
      type(flavor_t), intent(in) :: flv
      type(polarization_t), intent(out) :: pol
    end subroutine hepmc_polarization_to_pol
    module subroutine hepmc_polarization_to_hel (hpol, flv, hel)
      type(hepmc_polarization_t), intent(in) :: hpol
      type(flavor_t), intent(in) :: flv
      type(helicity_t), intent(out) :: hel
    end subroutine hepmc_polarization_to_hel
    module subroutine hepmc_particle_init (prt, p, pdg, status)
      type(hepmc_particle_t), intent(out) :: prt
      type(vector4_t), intent(in) :: p
      integer, intent(in) :: pdg, status
    end subroutine hepmc_particle_init
    module subroutine hepmc_particle_set_color_col (prt, col)
      type(hepmc_particle_t), intent(inout) :: prt
      type(color_t), intent(in) :: col
    end subroutine hepmc_particle_set_color_col
    module subroutine hepmc_particle_set_color_int (prt, col)
      type(hepmc_particle_t), intent(inout) :: prt
      integer, dimension(2), intent(in) :: col
    end subroutine hepmc_particle_set_color_int
    module subroutine hepmc_particle_set_polarization_pol (prt, pol)
      type(hepmc_particle_t), intent(inout) :: prt
      type(polarization_t), intent(in) :: pol
    end subroutine hepmc_particle_set_polarization_pol
    module subroutine hepmc_particle_set_polarization_hel (prt, hel)
      type(hepmc_particle_t), intent(inout) :: prt
      type(helicity_t), intent(in) :: hel
    end subroutine hepmc_particle_set_polarization_hel
    module subroutine hepmc_particle_set_polarization_int (prt, hel)
      type(hepmc_particle_t), intent(inout) :: prt
      integer, intent(in) :: hel
    end subroutine hepmc_particle_set_polarization_int
    module function hepmc_particle_get_barcode (prt) result (barcode)
      integer :: barcode
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_barcode
    module function hepmc_particle_get_momentum (prt) result (p)
      type(vector4_t) :: p
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_momentum
    module function hepmc_particle_get_mass_squared (prt) result (m2)
      real(default) :: m2
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_mass_squared
    module function hepmc_particle_get_pdg (prt) result (pdg)
      integer :: pdg
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_pdg
    module function hepmc_particle_get_status (prt) result (status)
      integer :: status
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_status
    module function hepmc_particle_is_beam (prt) result (is_beam)
      logical :: is_beam
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_is_beam
    module function hepmc_particle_get_production_vertex (prt) result (v)
      type(hepmc_vertex_t) :: v
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_production_vertex
    module function hepmc_particle_get_decay_vertex (prt) result (v)
      type(hepmc_vertex_t) :: v
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_decay_vertex
    module function hepmc_particle_get_parent_barcodes &
         (prt) result (parent_barcode)
      type(hepmc_particle_t), intent(in) :: prt
      integer, dimension(:), allocatable :: parent_barcode
    end function hepmc_particle_get_parent_barcodes
    module function hepmc_particle_get_child_barcodes &
         (prt) result (child_barcode)
      type(hepmc_particle_t), intent(in) :: prt
      integer, dimension(:), allocatable :: child_barcode
    end function hepmc_particle_get_child_barcodes
    module function hepmc_particle_get_polarization (prt) result (pol)
      type(hepmc_polarization_t) :: pol
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_polarization
    module function hepmc_particle_get_color (prt) result (col)
      integer, dimension(2) :: col
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_color
    module function hepmc_vertex_to_vertex (vtx) result (v)
      type(hepmc_vertex_t), intent(in) :: vtx
      type(vector4_t) :: v
    end function hepmc_vertex_to_vertex
    module subroutine hepmc_vertex_init (v, x)
      type(hepmc_vertex_t), intent(out) :: v
      type(vector4_t), intent(in), optional :: x
    end subroutine hepmc_vertex_init
    module function hepmc_vertex_is_valid (v) result (flag)
      logical :: flag
      type(hepmc_vertex_t), intent(in) :: v
    end function hepmc_vertex_is_valid
    module subroutine hepmc_vertex_add_particle_in (v, prt)
      type(hepmc_vertex_t), intent(inout) :: v
      type(hepmc_particle_t), intent(in) :: prt
    end subroutine hepmc_vertex_add_particle_in
    module subroutine hepmc_vertex_add_particle_out (v, prt)
      type(hepmc_vertex_t), intent(inout) :: v
      type(hepmc_particle_t), intent(in) :: prt
    end subroutine hepmc_vertex_add_particle_out
    module function hepmc_vertex_get_n_in (v) result (n_in)
      integer :: n_in
      type(hepmc_vertex_t), intent(in) :: v
    end function hepmc_vertex_get_n_in
    module function hepmc_vertex_get_n_out (v) result (n_out)
      integer :: n_out
      type(hepmc_vertex_t), intent(in) :: v
    end function hepmc_vertex_get_n_out
    module function hepmc_particle_get_parents (p) result (n_p)
      integer :: n_p
      type(hepmc_particle_t), intent(in) :: p
    end function hepmc_particle_get_parents
    module function hepmc_particle_get_children (p) result (n_ch)
      integer :: n_ch
      type(hepmc_particle_t), intent(in) :: p
    end function hepmc_particle_get_children
    module function hepmc_particle_get_n_parents (prt) result (n_parents)
      integer :: n_parents
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_n_parents
    module function hepmc_particle_get_n_children (prt) result (n_children)
      integer :: n_children
      type(hepmc_particle_t), intent(in) :: prt
    end function hepmc_particle_get_n_children
    module subroutine hepmc_vertex_particle_in_iterator_init (it, v)
      type(hepmc_vertex_particle_in_iterator_t), intent(out) :: it
      type(hepmc_vertex_t), intent(in) :: v
    end subroutine hepmc_vertex_particle_in_iterator_init
    module subroutine hepmc_vertex_particle_in_iterator_final (it)
      type(hepmc_vertex_particle_in_iterator_t), intent(inout) :: it
    end subroutine hepmc_vertex_particle_in_iterator_final
    module subroutine hepmc_vertex_particle_in_iterator_advance (it)
      type(hepmc_vertex_particle_in_iterator_t), intent(inout) :: it
    end subroutine hepmc_vertex_particle_in_iterator_advance
    module subroutine hepmc_vertex_particle_in_iterator_reset (it)
      type(hepmc_vertex_particle_in_iterator_t), intent(inout) :: it
    end subroutine hepmc_vertex_particle_in_iterator_reset
    module function hepmc_vertex_particle_in_iterator_is_valid &
         (it) result (flag)
      logical :: flag
      type(hepmc_vertex_particle_in_iterator_t), intent(in) :: it
    end function hepmc_vertex_particle_in_iterator_is_valid
    module function hepmc_vertex_particle_in_iterator_get (it) result (prt)
      type(hepmc_particle_t) :: prt
      type(hepmc_vertex_particle_in_iterator_t), intent(in) :: it
    end function hepmc_vertex_particle_in_iterator_get
    module function hepmc_vertex_get_nth_particle_in (vtx, n) result (prt)
      type(hepmc_particle_t) :: prt
      type(hepmc_vertex_t), intent(in) :: vtx
      integer, intent(in) :: n
    end function hepmc_vertex_get_nth_particle_in
    module function hepmc_vertex_get_nth_particle_out (vtx, n) result (prt)
      type(hepmc_particle_t) :: prt
      type(hepmc_vertex_t), intent(in) :: vtx
      integer, intent(in) :: n
    end function hepmc_vertex_get_nth_particle_out
    module subroutine hepmc_vertex_particle_out_iterator_init (it, v)
      type(hepmc_vertex_particle_out_iterator_t), intent(out) :: it
      type(hepmc_vertex_t), intent(in) :: v
    end subroutine hepmc_vertex_particle_out_iterator_init
    module subroutine hepmc_vertex_particle_out_iterator_final (it)
      type(hepmc_vertex_particle_out_iterator_t), intent(inout) :: it
    end subroutine hepmc_vertex_particle_out_iterator_final
    module subroutine hepmc_vertex_particle_out_iterator_advance (it)
      type(hepmc_vertex_particle_out_iterator_t), intent(inout) :: it
    end subroutine hepmc_vertex_particle_out_iterator_advance
    module subroutine hepmc_vertex_particle_out_iterator_reset (it)
      type(hepmc_vertex_particle_out_iterator_t), intent(inout) :: it
    end subroutine hepmc_vertex_particle_out_iterator_reset
    module function hepmc_vertex_particle_out_iterator_is_valid &
         (it) result (flag)
      logical :: flag
      type(hepmc_vertex_particle_out_iterator_t), intent(in) :: it
    end function hepmc_vertex_particle_out_iterator_is_valid
    module function hepmc_vertex_particle_out_iterator_get (it) result (prt)
      type(hepmc_particle_t) :: prt
      type(hepmc_vertex_particle_out_iterator_t), intent(in) :: it
    end function hepmc_vertex_particle_out_iterator_get
    module subroutine hepmc_event_init (evt, proc_id, event_id)
      type(hepmc_event_t), intent(out) :: evt
      integer, intent(in), optional :: proc_id, event_id
    end subroutine hepmc_event_init
    module subroutine hepmc_event_final (evt)
      type(hepmc_event_t), intent(inout) :: evt
    end subroutine hepmc_event_final
    module subroutine hepmc_event_nullify (evt)
      type(hepmc_event_t), intent(inout) :: evt
    end subroutine hepmc_event_nullify
    module function hepmc_event_get_c_ptr (evt) result (p)
      type(hepmc_event_t), intent(in) :: evt
      type(c_ptr) :: p
    end function hepmc_event_get_c_ptr
    module subroutine hepmc_event_print (evt)
      type(hepmc_event_t), intent(in) :: evt
    end subroutine hepmc_event_print
    module function hepmc_event_get_event_index (evt) result (i_proc)
      integer :: i_proc
      type(hepmc_event_t), intent(in) :: evt
    end function hepmc_event_get_event_index
    module function hepmc_event_get_n_particles (evt) result (n_tot)
      integer :: n_tot
      type(hepmc_event_t), intent(in) :: evt
    end function hepmc_event_get_n_particles
    module function hepmc_event_get_n_beams (evt) result (n_tot)
      integer :: n_tot
      type(hepmc_event_t), intent(in) :: evt
    end function hepmc_event_get_n_beams
    module subroutine hepmc_event_set_process_id (evt, proc)
      type(hepmc_event_t), intent(in) :: evt
      integer, intent(in) :: proc
    end subroutine hepmc_event_set_process_id
    module function hepmc_event_get_process_id (evt) result (i_proc)
      integer :: i_proc
      type(hepmc_event_t), intent(in) :: evt
    end function hepmc_event_get_process_id
    module subroutine hepmc_event_set_scale (evt, scale)
      type(hepmc_event_t), intent(in) :: evt
      real(default), intent(in) :: scale
    end subroutine hepmc_event_set_scale
    module function hepmc_event_get_scale (evt) result (scale)
      real(default) :: scale
      type(hepmc_event_t), intent(in) :: evt
    end function hepmc_event_get_scale
    module subroutine hepmc_event_set_alpha_qcd (evt, alpha)
      type(hepmc_event_t), intent(in) :: evt
      real(default), intent(in) :: alpha
    end subroutine hepmc_event_set_alpha_qcd
    module function hepmc_event_get_alpha_qcd (evt) result (alpha)
      real(default) :: alpha
      type(hepmc_event_t), intent(in) :: evt
    end function hepmc_event_get_alpha_qcd
    module subroutine hepmc_event_set_alpha_qed (evt, alpha)
      type(hepmc_event_t), intent(in) :: evt
      real(default), intent(in) :: alpha
    end subroutine hepmc_event_set_alpha_qed
    module function hepmc_event_get_alpha_qed (evt) result (alpha)
      real(default) :: alpha
      type(hepmc_event_t), intent(in) :: evt
    end function hepmc_event_get_alpha_qed
    module subroutine hepmc_event_clear_weights (evt)
      type(hepmc_event_t), intent(in) :: evt
    end subroutine hepmc_event_clear_weights
    module subroutine hepmc_event_add_weight (evt, weight, rescale)
      type(hepmc_event_t), intent(in) :: evt
      real(default), intent(in) :: weight
      logical, intent(in) :: rescale
    end subroutine hepmc_event_add_weight
    module function hepmc_event_get_weights_size (evt) result (n)
      integer :: n
      type(hepmc_event_t), intent(in) :: evt
    end function hepmc_event_get_weights_size
    module function hepmc_event_get_weight (evt, index, rescale) result (weight)
      real(default) :: weight
      type(hepmc_event_t), intent(in) :: evt
      integer, intent(in) :: index
      logical, intent(in) :: rescale
    end function hepmc_event_get_weight
    module subroutine hepmc_event_add_vertex (evt, v)
      type(hepmc_event_t), intent(inout) :: evt
      type(hepmc_vertex_t), intent(in) :: v
    end subroutine hepmc_event_add_vertex
    module subroutine hepmc_event_set_signal_process_vertex (evt, v)
      type(hepmc_event_t), intent(inout) :: evt
      type(hepmc_vertex_t), intent(in) :: v
    end subroutine hepmc_event_set_signal_process_vertex
    module function hepmc_event_get_signal_process_vertex (evt) result (v)
      type(hepmc_event_t), intent(in) :: evt
      type(hepmc_vertex_t) :: v
    end function hepmc_event_get_signal_process_vertex
    module subroutine hepmc_event_set_beam_particles (evt, prt1, prt2)
      type(hepmc_event_t), intent(inout) :: evt
      type(hepmc_particle_t), intent(in) :: prt1, prt2
      logical(c_bool) :: flag
    end subroutine hepmc_event_set_beam_particles
    module subroutine hepmc_event_set_cross_section (evt, xsec, xsec_err)
      type(hepmc_event_t), intent(inout) :: evt
      real(default), intent(in) :: xsec, xsec_err
    end subroutine hepmc_event_set_cross_section
    module subroutine hepmc_event_particle_iterator_init (it, evt)
      type(hepmc_event_particle_iterator_t), intent(out) :: it
      type(hepmc_event_t), intent(in) :: evt
    end subroutine hepmc_event_particle_iterator_init
    module subroutine hepmc_event_particle_iterator_final (it)
      type(hepmc_event_particle_iterator_t), intent(inout) :: it
    end subroutine hepmc_event_particle_iterator_final
    module subroutine hepmc_event_particle_iterator_advance (it)
      type(hepmc_event_particle_iterator_t), intent(inout) :: it
    end subroutine hepmc_event_particle_iterator_advance
    module subroutine hepmc_event_particle_iterator_reset (it)
      type(hepmc_event_particle_iterator_t), intent(inout) :: it
    end subroutine hepmc_event_particle_iterator_reset
    module function hepmc_event_particle_iterator_is_valid (it) result (flag)
      logical :: flag
      type(hepmc_event_particle_iterator_t), intent(in) :: it
    end function hepmc_event_particle_iterator_is_valid
    module function hepmc_event_particle_iterator_get (it) result (prt)
      type(hepmc_particle_t) :: prt
      type(hepmc_event_particle_iterator_t), intent(in) :: it
    end function hepmc_event_particle_iterator_get
    module function hepmc_event_get_nth_particle (evt, n) result (prt)
      type(hepmc_particle_t) :: prt
      type(hepmc_event_t), intent(in) :: evt
      integer, intent(in) :: n
    end function hepmc_event_get_nth_particle
    module function hepmc_event_get_nth_beam (evt, n) result (beam_barcode)
      integer :: beam_barcode
      type(hepmc_event_t), intent(in) :: evt
      integer, intent(in) :: n
    end function hepmc_event_get_nth_beam
    module subroutine hepmc_iostream_open_out (iostream, filename, hepmc3_mode)
      type(hepmc_iostream_t), intent(out) :: iostream
      type(string_t), intent(in) :: filename
      integer, intent(in) :: hepmc3_mode
    end subroutine hepmc_iostream_open_out
    module subroutine hepmc_iostream_open_in (iostream, filename, hepmc3_mode)
      type(hepmc_iostream_t), intent(out) :: iostream
      type(string_t), intent(in) :: filename
      integer, intent(in) :: hepmc3_mode
    end subroutine hepmc_iostream_open_in
    module subroutine hepmc_iostream_close (iostream)
      type(hepmc_iostream_t), intent(inout) :: iostream
    end subroutine hepmc_iostream_close
    module subroutine hepmc_iostream_write_event (iostream, evt, hepmc3_mode)
      type(hepmc_iostream_t), intent(inout) :: iostream
      type(hepmc_event_t), intent(in) :: evt
      integer, intent(in), optional :: hepmc3_mode
    end subroutine hepmc_iostream_write_event
    module subroutine hepmc_iostream_read_event (iostream, evt, ok)
      type(hepmc_iostream_t), intent(inout) :: iostream
      type(hepmc_event_t), intent(inout) :: evt
      logical, intent(out) :: ok
    end subroutine hepmc_iostream_read_event
  end interface

end module hepmc_interface
