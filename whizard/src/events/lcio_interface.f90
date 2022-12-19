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

module lcio_interface

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

  public :: lcio_is_available
  public :: lcio_run_header_t
  public :: lcio_run_header_init
  public :: lcio_run_header_write
  public :: lccollection_t
  public :: lcio_event_t
  public :: lcio_event_init
  public :: show_lcio_event
  public :: write_lcio_event
  public :: lcio_event_final
  public :: lcio_event_nullify
  public :: lcio_event_get_c_ptr
  public :: lcio_event_set_weight
  public :: lcio_event_set_alpha_qcd
  public :: lcio_event_set_scale
  public :: lcio_event_set_sqrts
  public :: lcio_event_set_xsec
  public :: lcio_event_set_beam
  public :: lcio_event_set_polarization
  public :: lcio_event_set_beam_file
  public :: lcio_event_set_process_name
  public :: lcio_event_set_alt_sqme
  public :: lcio_event_set_sqme
  public :: lcio_event_set_alt_weight
  public :: lcio_event_add_coll
  public :: lcio_particle_t
  public :: lcio_particle_add_to_evt_coll
  public :: lcio_particle_init
  public :: lcio_particle_set_color
  public :: lcio_particle_get_flow
  public :: lcio_particle_get_momentum
  public :: lcio_particle_get_mass_squared
  public :: lcio_particle_get_vertex
  public :: lcio_particle_get_time
  public :: lcio_polarization_init
  public :: lcio_particle_to_pol
  public :: lcio_particle_to_hel
  public :: lcio_particle_set_vtx
  public :: lcio_particle_set_t
  public :: lcio_particle_set_parent
  public :: lcio_particle_get_status
  public :: lcio_particle_get_pdg
  public :: lcio_particle_get_n_parents
  public :: lcio_particle_get_n_children
  public :: lcio_get_n_parents
  public :: lcio_get_n_children
  public :: lcio_writer_t
  public :: lcio_writer_open_out
  public :: lcio_writer_close
  public :: lcio_event_write
  public :: lcio_reader_t
  public :: lcio_open_file
  public :: lcio_reader_close
  public :: lcio_read_event
  public :: lcio_event_get_event_index
  public :: lcio_event_get_process_id
  public :: lcio_event_get_n_tot
  public :: lcio_event_get_alphas
  public :: lcio_event_get_scaleval
  public :: lcio_event_get_particle

  type :: lcio_run_header_t
     private
     type(c_ptr) :: obj
  end type lcio_run_header_t

  type :: lccollection_t
     private
     type(c_ptr) :: obj = c_null_ptr
  end type lccollection_t

  type, extends (event_handle_t) :: lcio_event_t
     private
     type(c_ptr) :: obj = c_null_ptr
     type(lccollection_t) :: lccoll
  end type lcio_event_t

  type :: lcio_particle_t
     private
     type(c_ptr) :: obj
  end type lcio_particle_t

  type :: lcio_writer_t
     private
     type(c_ptr) :: obj
  end type lcio_writer_t

  type :: lcio_reader_t
     private
     type(c_ptr) :: obj
  end type lcio_reader_t


  interface
     logical(c_bool) function lcio_available () bind(C)
       import
     end function lcio_available
  end interface
  interface
     type(c_ptr) function new_lcio_run_header (proc_id) bind(C)
       import
       integer(c_int), value :: proc_id
     end function new_lcio_run_header
  end interface
  interface
     subroutine run_header_set_simstring &
          (runhdr_obj, simstring) bind(C)
       import
       type(c_ptr), value :: runhdr_obj
       character(c_char), dimension(*), intent(in) :: simstring
     end subroutine run_header_set_simstring
  end interface
  interface
     subroutine write_run_header (lcwrt_obj, runhdr_obj) bind(C)
       import
       type(c_ptr), value :: lcwrt_obj
       type(c_ptr), value :: runhdr_obj
     end subroutine write_run_header
  end interface
  interface
     type(c_ptr) function new_lccollection () bind(C)
       import
     end function new_lccollection
  end interface
  interface
     type(c_ptr) function new_lcio_event (proc_id, event_id, run_id) bind(C)
       import
       integer(c_int), value :: proc_id, event_id, run_id
     end function new_lcio_event
  end interface
  interface
     subroutine lcio_event_delete (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end subroutine lcio_event_delete
  end interface
  interface
     subroutine dump_lcio_event (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end subroutine dump_lcio_event
  end interface
  interface
     subroutine lcio_event_to_file (evt_obj, filename) bind(C)
       import
       type(c_ptr), value :: evt_obj
       character(c_char), dimension(*), intent(in) :: filename
     end subroutine lcio_event_to_file
  end interface
  interface
     subroutine lcio_set_weight (evt_obj, weight) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: weight
     end subroutine lcio_set_weight
  end interface
  interface
     subroutine lcio_set_alpha_qcd (evt_obj, alphas) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: alphas
     end subroutine lcio_set_alpha_qcd
  end interface
  interface
     subroutine lcio_set_scale (evt_obj, scale) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: scale
     end subroutine lcio_set_scale
  end interface
  interface
     subroutine lcio_set_sqrts (evt_obj, sqrts) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: sqrts
     end subroutine lcio_set_sqrts
  end interface
  interface
     subroutine lcio_set_xsec (evt_obj, xsec, xsec_err) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: xsec, xsec_err
     end subroutine lcio_set_xsec
  end interface
  interface
     subroutine lcio_set_beam (evt_obj, pdg, beam) bind(C)
       import
       type(c_ptr), value :: evt_obj
       integer(c_int), value :: pdg, beam
     end subroutine lcio_set_beam
  end interface
  interface
     subroutine lcio_set_pol (evt_obj, pol, beam) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: pol
       integer(c_int), value :: beam
     end subroutine lcio_set_pol
  end interface
  interface
     subroutine lcio_set_beam_file (evt_obj, file) bind(C)
       import
       type(c_ptr), value :: evt_obj
       character(len=1, kind=c_char), dimension(*), intent(in) :: file
     end subroutine lcio_set_beam_file
  end interface
  interface
     subroutine lcio_set_process_name (evt_obj, name) bind(C)
       import
       type(c_ptr), value :: evt_obj
       character(len=1, kind=c_char), dimension(*), intent(in) :: name
     end subroutine lcio_set_process_name
  end interface
  interface
     subroutine lcio_set_sqme (evt_obj, sqme) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: sqme
     end subroutine lcio_set_sqme
  end interface
  interface
     subroutine lcio_set_alt_sqme (evt_obj, sqme, index) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: sqme
       integer(c_int), value :: index
     end subroutine lcio_set_alt_sqme
  end interface
  interface
     subroutine lcio_set_alt_weight (evt_obj, weight, index) bind(C)
       import
       type(c_ptr), value :: evt_obj
       real(c_double), value :: weight
       integer(c_int), value :: index
     end subroutine lcio_set_alt_weight
  end interface
  interface
     subroutine lcio_event_add_collection &
          (evt_obj, lccoll_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj, lccoll_obj
     end subroutine lcio_event_add_collection
  end interface
  interface
     type(c_ptr) function new_lcio_particle &
          (px, py, pz, pdg_id, mass, charge, status) bind(C)
       import
       integer(c_int), value :: pdg_id, status
       real(c_double), value :: px, py, pz, mass, charge
     end function new_lcio_particle
  end interface
  interface
     subroutine add_particle_to_collection &
          (prt_obj, lccoll_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj, lccoll_obj
     end subroutine add_particle_to_collection
  end interface
  interface
     subroutine lcio_set_color_flow (prt_obj, col1, col2) bind(C)
       import
       type(c_ptr), value :: prt_obj
       integer(c_int), value :: col1, col2
     end subroutine lcio_set_color_flow
  end interface
  interface lcio_particle_set_color
     module procedure lcio_particle_set_color_col
     module procedure lcio_particle_set_color_int
  end interface lcio_particle_set_color
  interface
     integer(c_int) function lcio_particle_flow (prt_obj, col_index) bind(C)
       use iso_c_binding !NODEP!
       type(c_ptr), value :: prt_obj
       integer(c_int), value :: col_index
     end function lcio_particle_flow
  end interface
  interface
     real(c_double) function lcio_three_momentum (prt_obj, p_index) bind(C)
       use iso_c_binding !NODEP!
       type(c_ptr), value :: prt_obj
       integer(c_int), value :: p_index
     end function lcio_three_momentum
  end interface
  interface
     real(c_double) function lcio_energy (prt_obj) bind(C)
       import
       type(c_ptr), intent(in), value :: prt_obj
     end function lcio_energy
  end interface
  interface
     function lcio_mass (prt_obj) result (mass) bind(C)
       import
       real(c_double) :: mass
       type(c_ptr), value :: prt_obj
     end function lcio_mass
  end interface
  interface
     real(c_double) function lcio_vtx_x (prt) bind(C)
       import
       type(c_ptr), value :: prt
     end function lcio_vtx_x
  end interface
  interface
     real(c_double) function lcio_vtx_y (prt) bind(C)
       import
       type(c_ptr), value :: prt
     end function lcio_vtx_y
  end interface
  interface
     real(c_double) function lcio_vtx_z (prt) bind(C)
       import
       type(c_ptr), value :: prt
     end function lcio_vtx_z
  end interface
  interface
     real(c_float) function lcio_prt_time (prt) bind(C)
       import
       type(c_ptr), value :: prt
     end function lcio_prt_time
  end interface
  interface
     subroutine lcio_particle_set_spin (prt_obj, s1, s2, s3) bind(C)
       import
       type(c_ptr), value :: prt_obj
       real(c_double), value :: s1, s2, s3
     end subroutine lcio_particle_set_spin
  end interface
  interface lcio_polarization_init
     module procedure lcio_polarization_init_pol
     module procedure lcio_polarization_init_hel
     module procedure lcio_polarization_init_int
  end interface
  interface
     function lcio_polarization_degree (prt_obj) result (degree) bind(C)
       import
       real(c_double) :: degree
       type(c_ptr), value :: prt_obj
     end function lcio_polarization_degree
  end interface
  interface
     function lcio_polarization_theta (prt_obj) result (theta) bind(C)
       import
       real(c_double) :: theta
       type(c_ptr), value :: prt_obj
     end function lcio_polarization_theta
  end interface
  interface
     function lcio_polarization_phi (prt_obj) result (phi) bind(C)
       import
       real(c_double) :: phi
       type(c_ptr), value :: prt_obj
     end function lcio_polarization_phi
  end interface
  interface
     subroutine lcio_particle_set_vertex (prt_obj, vx, vy, vz) bind(C)
       import
       type(c_ptr), value :: prt_obj
       real(c_double), value :: vx, vy, vz
     end subroutine lcio_particle_set_vertex
  end interface
  interface
     subroutine lcio_particle_set_time (prt_obj, t) bind(C)
       import
       type(c_ptr), value :: prt_obj
       real(c_float), value :: t
     end subroutine lcio_particle_set_time
  end interface

  interface
     subroutine lcio_particle_add_parent (prt_obj1, prt_obj2) bind(C)
       import
       type(c_ptr), value :: prt_obj1, prt_obj2
     end subroutine lcio_particle_add_parent
  end interface
  interface
     integer(c_int) function lcio_particle_get_generator_status &
          (prt_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj
     end function lcio_particle_get_generator_status
  end interface
  interface
     integer(c_int) function lcio_particle_get_pdg_code (prt_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj
     end function lcio_particle_get_pdg_code
  end interface
  interface
     integer(c_int) function lcio_n_parents (prt_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj
     end function lcio_n_parents
  end interface
  interface
     integer(c_int) function lcio_n_daughters (prt_obj) bind(C)
       import
       type(c_ptr), value :: prt_obj
     end function lcio_n_daughters
  end interface
  interface
     integer(c_int) function lcio_event_parent_k &
          (evt_obj, num_part, k_parent) bind (C)
       use iso_c_binding !NODEP!
       type(c_ptr), value :: evt_obj
       integer(c_int), value :: num_part, k_parent
     end function lcio_event_parent_k
  end interface
  interface
     integer(c_int) function lcio_event_daughter_k &
          (evt_obj, num_part, k_daughter) bind (C)
       use iso_c_binding !NODEP!
       type(c_ptr), value :: evt_obj
       integer(c_int), value :: num_part, k_daughter
     end function lcio_event_daughter_k
  end interface
  interface
     type(c_ptr) function open_lcio_writer_new (filename, complevel) bind(C)
       import
       character(c_char), dimension(*), intent(in) :: filename
       integer(c_int), intent(in) :: complevel
     end function open_lcio_writer_new
  end interface
  interface
     subroutine lcio_writer_delete (io_obj) bind(C)
       import
       type(c_ptr), value :: io_obj
     end subroutine lcio_writer_delete
  end interface
  interface
     subroutine lcio_write_event (io_obj, evt_obj) bind(C)
       import
       type(c_ptr), value :: io_obj, evt_obj
     end subroutine lcio_write_event
  end interface
  interface
     type(c_ptr) function open_lcio_reader (filename) bind(C)
       import
       character(c_char), dimension(*), intent(in) :: filename
     end function open_lcio_reader
  end interface
  interface
     subroutine lcio_reader_delete (io_obj) bind(C)
       import
       type(c_ptr), value :: io_obj
     end subroutine lcio_reader_delete
  end interface
  interface
     type(c_ptr) function read_lcio_event (io_obj) bind(C)
       import
       type(c_ptr), value :: io_obj
     end function read_lcio_event
  end interface
  interface
     integer(c_int) function lcio_event_get_event_number (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function lcio_event_get_event_number
  end interface

  interface
     integer(c_int) function lcio_event_signal_process_id (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function lcio_event_signal_process_id
  end interface
  interface
     integer(c_int) function lcio_event_get_n_particles (evt_obj) bind(C)
       import
       type(c_ptr), value :: evt_obj
     end function lcio_event_get_n_particles
  end interface
  interface
     function lcio_event_get_alpha_qcd (evt_obj) result (as) bind(C)
       import
       real(c_double) :: as
       type(c_ptr), value :: evt_obj
     end function lcio_event_get_alpha_qcd
  end interface
  interface
     function lcio_event_get_scale (evt_obj) result (scale) bind(C)
       import
       real(c_double) :: scale
       type(c_ptr), value :: evt_obj
     end function lcio_event_get_scale
  end interface
  interface
     type(c_ptr) function lcio_event_particle_k (evt_obj, k) bind(C)
       import
       type(c_ptr), value :: evt_obj
       integer(c_int), value :: k
     end function lcio_event_particle_k
  end interface

  interface
    module function lcio_is_available () result (flag)
      logical :: flag
    end function lcio_is_available
    module subroutine lcio_run_header_init (runhdr, proc_id, run_id)
      type(lcio_run_header_t), intent(out) :: runhdr
      integer, intent(in), optional :: proc_id, run_id
    end subroutine lcio_run_header_init
    module subroutine lcio_run_header_write (wrt, hdr)
      type(lcio_writer_t), intent(inout) :: wrt
      type(lcio_run_header_t), intent(inout) :: hdr
    end subroutine lcio_run_header_write
    module subroutine lcio_event_init (evt, proc_id, event_id, run_id)
      type(lcio_event_t), intent(out) :: evt
      integer, intent(in), optional :: proc_id, event_id, run_id
    end subroutine lcio_event_init
    module subroutine show_lcio_event (evt)
      type(lcio_event_t), intent(in) :: evt
    end subroutine show_lcio_event
    module subroutine write_lcio_event (evt, filename)
      type(lcio_event_t), intent(in) :: evt
      type(string_t), intent(in) :: filename
    end subroutine write_lcio_event
    module subroutine lcio_event_final (evt, delete)
      type(lcio_event_t), intent(inout) :: evt
      logical, intent(in) :: delete
    end subroutine lcio_event_final
    module subroutine lcio_event_nullify (evt)
      type(lcio_event_t), intent(inout) :: evt
    end subroutine lcio_event_nullify
    module function lcio_event_get_c_ptr (evt) result (p)
      type(lcio_event_t), intent(in) :: evt
      type(c_ptr) :: p
    end function lcio_event_get_c_ptr
    module subroutine lcio_event_set_weight (evt, weight)
      type(lcio_event_t), intent(inout) :: evt
      real(default), intent(in) :: weight
    end subroutine lcio_event_set_weight
    module subroutine lcio_event_set_alpha_qcd (evt, alphas)
      type(lcio_event_t), intent(inout) :: evt
      real(default), intent(in) :: alphas
    end subroutine lcio_event_set_alpha_qcd
    module subroutine lcio_event_set_scale (evt, scale)
      type(lcio_event_t), intent(inout) :: evt
      real(default), intent(in) :: scale
    end subroutine lcio_event_set_scale
    module subroutine lcio_event_set_sqrts (evt, sqrts)
      type(lcio_event_t), intent(inout) :: evt
      real(default), intent(in) :: sqrts
    end subroutine lcio_event_set_sqrts
    module subroutine lcio_event_set_xsec (evt, xsec, xsec_err)
      type(lcio_event_t), intent(inout) :: evt
      real(default), intent(in) :: xsec, xsec_err
    end subroutine lcio_event_set_xsec
    module subroutine lcio_event_set_beam (evt, pdg, beam)
      type(lcio_event_t), intent(inout) :: evt
      integer, intent(in) :: pdg, beam
    end subroutine lcio_event_set_beam
    module subroutine lcio_event_set_polarization (evt, pol, beam)
      type(lcio_event_t), intent(inout) :: evt
      real(default), intent(in) :: pol
      integer, intent(in) :: beam
    end subroutine lcio_event_set_polarization
    module subroutine lcio_event_set_beam_file (evt, file)
      type(lcio_event_t), intent(inout) :: evt
      type(string_t), intent(in) :: file
    end subroutine lcio_event_set_beam_file
    module subroutine lcio_event_set_process_name (evt, name)
      type(lcio_event_t), intent(inout) :: evt
      type(string_t), intent(in) :: name
    end subroutine lcio_event_set_process_name
    module subroutine lcio_event_set_alt_sqme (evt, sqme, index)
      type(lcio_event_t), intent(inout) :: evt
      real(default), intent(in) :: sqme
      integer, intent(in) :: index
    end subroutine lcio_event_set_alt_sqme
    module subroutine lcio_event_set_sqme (evt, sqme)
      type(lcio_event_t), intent(inout) :: evt
      real(default), intent(in) :: sqme
    end subroutine lcio_event_set_sqme
    module subroutine lcio_event_set_alt_weight (evt, weight, index)
      type(lcio_event_t), intent(inout) :: evt
      real(default), intent(in) :: weight
      integer, intent(in) :: index
    end subroutine lcio_event_set_alt_weight
    module subroutine lcio_event_add_coll (evt)
      type(lcio_event_t), intent(inout) :: evt
    end subroutine lcio_event_add_coll
    module subroutine lcio_particle_add_to_evt_coll &
         (lprt, evt)
      type(lcio_particle_t), intent(in) :: lprt
      type(lcio_event_t), intent(inout) :: evt
    end subroutine lcio_particle_add_to_evt_coll
    module subroutine lcio_particle_init (prt, p, pdg, charge, status)
      type(lcio_particle_t), intent(out) :: prt
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: charge
      integer, intent(in) :: pdg, status
    end subroutine lcio_particle_init
    module subroutine lcio_particle_set_color_col (prt, col)
      type(lcio_particle_t), intent(inout) :: prt
      type(color_t), intent(in) :: col
    end subroutine lcio_particle_set_color_col
    module subroutine lcio_particle_set_color_int (prt, col)
      type(lcio_particle_t), intent(inout) :: prt
      integer, dimension(2), intent(in) :: col
    end subroutine lcio_particle_set_color_int
    module function lcio_particle_get_flow (prt) result (col)
      integer, dimension(2) :: col
      type(lcio_particle_t), intent(in) :: prt
    end function lcio_particle_get_flow
    module function lcio_particle_get_momentum (prt) result (p)
      type(vector4_t) :: p
      type(lcio_particle_t), intent(in) :: prt
    end function lcio_particle_get_momentum
    module function lcio_particle_get_mass_squared (prt) result (m2)
      real(default) :: m2
      type(lcio_particle_t), intent(in) :: prt
    end function lcio_particle_get_mass_squared
    module function lcio_particle_get_vertex (prt) result (vtx)
      type(vector3_t) :: vtx
      type(lcio_particle_t), intent(in) :: prt
    end function lcio_particle_get_vertex
    module function lcio_particle_get_time (prt) result (time)
      real(default) :: time
      type(lcio_particle_t), intent(in) :: prt
    end function lcio_particle_get_time
    module subroutine lcio_polarization_init_pol (prt, pol)
      type(lcio_particle_t), intent(inout) :: prt
      type(polarization_t), intent(in) :: pol
    end subroutine lcio_polarization_init_pol
    module subroutine lcio_polarization_init_hel (prt, hel)
      type(lcio_particle_t), intent(inout) :: prt
      type(helicity_t), intent(in) :: hel
    end subroutine lcio_polarization_init_hel
    module subroutine lcio_polarization_init_int (prt, hel)
      type(lcio_particle_t), intent(inout) :: prt
      integer, intent(in) :: hel
    end subroutine lcio_polarization_init_int
    module subroutine lcio_particle_to_pol (prt, flv, pol)
      type(lcio_particle_t), intent(in) :: prt
      type(flavor_t), intent(in) :: flv
      type(polarization_t), intent(out) :: pol
    end subroutine lcio_particle_to_pol
    module subroutine lcio_particle_to_hel (prt, flv, hel)
      type(lcio_particle_t), intent(in) :: prt
      type(flavor_t), intent(in) :: flv
      type(helicity_t), intent(out) :: hel
    end subroutine lcio_particle_to_hel
    module subroutine lcio_particle_set_vtx (prt, vtx)
      type(lcio_particle_t), intent(inout) :: prt
      type(vector3_t), intent(in) :: vtx
    end subroutine lcio_particle_set_vtx
    module subroutine lcio_particle_set_t (prt, t)
      type(lcio_particle_t), intent(inout) :: prt
      real(default), intent(in) :: t
    end subroutine lcio_particle_set_t
    module subroutine lcio_particle_set_parent (daughter, parent)
      type(lcio_particle_t), intent(inout) :: daughter, parent
    end subroutine lcio_particle_set_parent
    module function lcio_particle_get_status (lptr) result (status)
      integer :: status
      type(lcio_particle_t), intent(in) :: lptr
    end function lcio_particle_get_status
    module function lcio_particle_get_pdg (lptr) result (pdg)
      integer :: pdg
      type(lcio_particle_t), intent(in) :: lptr
    end function lcio_particle_get_pdg
    module function lcio_particle_get_n_parents (lptr) result (n_parents)
      integer :: n_parents
      type(lcio_particle_t), intent(in) :: lptr
    end function lcio_particle_get_n_parents
    module function lcio_particle_get_n_children (lptr) result (n_children)
      integer :: n_children
      type(lcio_particle_t), intent(in) :: lptr
    end function lcio_particle_get_n_children
    module function lcio_get_n_parents &
         (evt, num_part, k_parent) result (index_parent)
      type(lcio_event_t), intent(in) :: evt
      integer, intent(in) :: num_part, k_parent
      integer :: index_parent
    end function lcio_get_n_parents
    module function lcio_get_n_children &
         (evt, num_part, k_daughter) result (index_daughter)
      type(lcio_event_t), intent(in) :: evt
      integer, intent(in) :: num_part, k_daughter
      integer :: index_daughter
    end function lcio_get_n_children
    module subroutine lcio_writer_open_out (lcio_writer, filename)
      type(lcio_writer_t), intent(out) :: lcio_writer
      type(string_t), intent(in) :: filename
    end subroutine lcio_writer_open_out
    module subroutine lcio_writer_close (lciowriter)
      type(lcio_writer_t), intent(inout) :: lciowriter
    end subroutine lcio_writer_close
    module subroutine lcio_event_write (wrt, evt)
      type(lcio_writer_t), intent(inout) :: wrt
      type(lcio_event_t), intent(in) :: evt
    end subroutine lcio_event_write
    module subroutine lcio_open_file (lcio_reader, filename)
      type(lcio_reader_t), intent(out) :: lcio_reader
      type(string_t), intent(in) :: filename
    end subroutine lcio_open_file
    module subroutine lcio_reader_close (lcioreader)
      type(lcio_reader_t), intent(inout) :: lcioreader
    end subroutine lcio_reader_close
    module subroutine lcio_read_event (lcrdr, evt, ok)
      type(lcio_reader_t), intent(inout) :: lcrdr
      type(lcio_event_t), intent(out) :: evt
      logical, intent(out) :: ok
    end subroutine lcio_read_event
    module function lcio_event_get_event_index (evt) result (i_evt)
      integer :: i_evt
      type(lcio_event_t), intent(in) :: evt
    end function lcio_event_get_event_index
    module function lcio_event_get_process_id (evt) result (i_proc)
      integer :: i_proc
      type(lcio_event_t), intent(in) :: evt
    end function lcio_event_get_process_id
    module function lcio_event_get_n_tot (evt) result (n_tot)
      integer :: n_tot
      type(lcio_event_t), intent(in) :: evt
    end function lcio_event_get_n_tot
    module function lcio_event_get_alphas (evt) result (as)
      type(lcio_event_t), intent(in) :: evt
      real(default) :: as
    end function lcio_event_get_alphas
    module function lcio_event_get_scaleval (evt) result (scale)
      type(lcio_event_t), intent(in) :: evt
      real(default) :: scale
    end function lcio_event_get_scaleval
    module function lcio_event_get_particle (evt, n) result (prt)
      type(lcio_event_t), intent(in) :: evt
      integer, intent(in) :: n
      type(lcio_particle_t) :: prt
    end function lcio_event_get_particle
  end interface

end module lcio_interface
