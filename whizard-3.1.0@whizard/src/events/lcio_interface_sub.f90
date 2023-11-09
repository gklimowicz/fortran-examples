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

submodule (lcio_interface) lcio_interface_s

  use constants, only: PI
  use physics_defs, only: ns_per_mm
  use diagnostics

  implicit none

contains

  module function lcio_is_available () result (flag)
    logical :: flag
    flag = lcio_available ()
  end function lcio_is_available

  module subroutine lcio_run_header_init (runhdr, proc_id, run_id)
    type(lcio_run_header_t), intent(out) :: runhdr
    integer, intent(in), optional :: proc_id, run_id
    integer(c_int) :: rid
    rid = 0; if (present (run_id))  rid = run_id
    runhdr%obj = new_lcio_run_header (rid)
    call run_header_set_simstring (runhdr%obj, &
         "WHIZARD version:" // "3.1.0")
  end subroutine lcio_run_header_init

  module subroutine lcio_run_header_write (wrt, hdr)
    type(lcio_writer_t), intent(inout) :: wrt
    type(lcio_run_header_t), intent(inout) :: hdr
    call write_run_header (wrt%obj, hdr%obj)
  end subroutine lcio_run_header_write

  module subroutine lcio_event_init (evt, proc_id, event_id, run_id)
    type(lcio_event_t), intent(out) :: evt
    integer, intent(in), optional :: proc_id, event_id, run_id
    integer(c_int) :: pid, eid, rid
    pid = 0;  if (present (proc_id))  pid = proc_id
    eid = 0;  if (present (event_id)) eid = event_id
    rid = 0;  if (present (run_id))   rid = run_id
    evt%obj = new_lcio_event (pid, eid, rid)
    evt%lccoll%obj = new_lccollection ()
  end subroutine lcio_event_init

  module subroutine show_lcio_event (evt)
    type(lcio_event_t), intent(in) :: evt
    if (c_associated (evt%obj)) then
       call dump_lcio_event (evt%obj)
    else
       call msg_error ("LCIO event is not allocated.")
    end if
  end subroutine show_lcio_event

  module subroutine write_lcio_event (evt, filename)
    type(lcio_event_t), intent(in) :: evt
    type(string_t), intent(in) :: filename
    call lcio_event_to_file (evt%obj, char (filename) // c_null_char)
  end subroutine write_lcio_event

  module subroutine lcio_event_final (evt, delete)
    type(lcio_event_t), intent(inout) :: evt
    logical, intent(in) :: delete
    if (c_associated (evt%obj)) then
       if (delete)  call lcio_event_delete (evt%obj)
       evt%obj = c_null_ptr
       evt%lccoll%obj = c_null_ptr
    end if
  end subroutine lcio_event_final

  module subroutine lcio_event_nullify (evt)
    type(lcio_event_t), intent(inout) :: evt
    evt%obj = c_null_ptr
    evt%lccoll%obj = c_null_ptr
  end subroutine lcio_event_nullify

  module function lcio_event_get_c_ptr (evt) result (p)
    type(lcio_event_t), intent(in) :: evt
    type(c_ptr) :: p
    p = evt%obj
  end function lcio_event_get_c_ptr

  module subroutine lcio_event_set_weight (evt, weight)
    type(lcio_event_t), intent(inout) :: evt
    real(default), intent(in) :: weight
    call lcio_set_weight (evt%obj, real (weight, c_double))
  end subroutine lcio_event_set_weight

  module subroutine lcio_event_set_alpha_qcd (evt, alphas)
    type(lcio_event_t), intent(inout) :: evt
    real(default), intent(in) :: alphas
    call lcio_set_alpha_qcd (evt%obj, real (alphas, c_double))
  end subroutine lcio_event_set_alpha_qcd

  module subroutine lcio_event_set_scale (evt, scale)
    type(lcio_event_t), intent(inout) :: evt
    real(default), intent(in) :: scale
    call lcio_set_scale (evt%obj, real (scale, c_double))
  end subroutine lcio_event_set_scale

  module subroutine lcio_event_set_sqrts (evt, sqrts)
    type(lcio_event_t), intent(inout) :: evt
    real(default), intent(in) :: sqrts
    call lcio_set_sqrts (evt%obj, real (sqrts, c_double))
  end subroutine lcio_event_set_sqrts

  module subroutine lcio_event_set_xsec (evt, xsec, xsec_err)
    type(lcio_event_t), intent(inout) :: evt
    real(default), intent(in) :: xsec, xsec_err
    call lcio_set_xsec (evt%obj, &
         real (xsec, c_double), real (xsec_err, c_double))
  end subroutine lcio_event_set_xsec

  module subroutine lcio_event_set_beam (evt, pdg, beam)
    type(lcio_event_t), intent(inout) :: evt
    integer, intent(in) :: pdg, beam
    call lcio_set_beam (evt%obj, &
         int (pdg, c_int), int (beam, c_int))
  end subroutine lcio_event_set_beam

  module subroutine lcio_event_set_polarization (evt, pol, beam)
    type(lcio_event_t), intent(inout) :: evt
    real(default), intent(in) :: pol
    integer, intent(in) :: beam
    call lcio_set_pol (evt%obj, real (pol, c_double), &
         int (beam, c_int))
  end subroutine lcio_event_set_polarization

  module subroutine lcio_event_set_beam_file (evt, file)
    type(lcio_event_t), intent(inout) :: evt
    type(string_t), intent(in) :: file
    call lcio_set_beam_file (evt%obj, &
         char (file) // c_null_char)
  end subroutine lcio_event_set_beam_file

  module subroutine lcio_event_set_process_name (evt, name)
    type(lcio_event_t), intent(inout) :: evt
    type(string_t), intent(in) :: name
    call lcio_set_process_name (evt%obj, &
         char (name) // c_null_char)
  end subroutine lcio_event_set_process_name

  module subroutine lcio_event_set_alt_sqme (evt, sqme, index)
    type(lcio_event_t), intent(inout) :: evt
    real(default), intent(in) :: sqme
    integer, intent(in) :: index
    call lcio_set_alt_sqme (evt%obj, real (sqme, c_double), &
         int (index, c_int))
  end subroutine lcio_event_set_alt_sqme

  module subroutine lcio_event_set_sqme (evt, sqme)
    type(lcio_event_t), intent(inout) :: evt
    real(default), intent(in) :: sqme
    call lcio_set_sqme (evt%obj, real (sqme, c_double))
  end subroutine lcio_event_set_sqme

  module subroutine lcio_event_set_alt_weight (evt, weight, index)
    type(lcio_event_t), intent(inout) :: evt
    real(default), intent(in) :: weight
    integer, intent(in) :: index
    call lcio_set_alt_weight (evt%obj, real (weight, c_double), &
         int (index, c_int))
  end subroutine lcio_event_set_alt_weight

  module subroutine lcio_event_add_coll (evt)
    type(lcio_event_t), intent(inout) :: evt
    call lcio_event_add_collection (evt%obj, &
         evt%lccoll%obj)
  end subroutine lcio_event_add_coll

  module subroutine lcio_particle_add_to_evt_coll &
       (lprt, evt)
    type(lcio_particle_t), intent(in) :: lprt
    type(lcio_event_t), intent(inout) :: evt
    call add_particle_to_collection (lprt%obj, evt%lccoll%obj)
  end subroutine lcio_particle_add_to_evt_coll

  module subroutine lcio_particle_init (prt, p, pdg, charge, status)
    type(lcio_particle_t), intent(out) :: prt
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: charge
    real(default) :: mass
    real(default) :: px, py, pz
    integer, intent(in) :: pdg, status
    px = vector4_get_component (p, 1)
    py = vector4_get_component (p, 2)
    pz = vector4_get_component (p, 3)
    mass = p**1
    prt%obj = new_lcio_particle (real (px, c_double), real (py, c_double), &
         real (pz, c_double), int (pdg, c_int), &
         real (mass, c_double), real (charge, c_double), int (status, c_int))
  end subroutine lcio_particle_init

  module subroutine lcio_particle_set_color_col (prt, col)
    type(lcio_particle_t), intent(inout) :: prt
    type(color_t), intent(in) :: col
    integer(c_int), dimension(2) :: c
    c(1) = col%get_col ()
    c(2) = col%get_acl ()
    if (c(1) /= 0 .or. c(2) /= 0)  then
       call lcio_set_color_flow (prt%obj, c(1), c(2))
    end  if
  end subroutine lcio_particle_set_color_col

  module subroutine lcio_particle_set_color_int (prt, col)
    type(lcio_particle_t), intent(inout) :: prt
    integer, dimension(2), intent(in) :: col
    integer(c_int), dimension(2) :: c
    c = col
    if (c(1) /= 0 .or. c(2) /= 0) then
       call lcio_set_color_flow (prt%obj, c(1), c(2))
    end if
  end subroutine lcio_particle_set_color_int

  module function lcio_particle_get_flow (prt) result (col)
    integer, dimension(2) :: col
    type(lcio_particle_t), intent(in) :: prt
    col(1) = lcio_particle_flow (prt%obj, 0_c_int)
    col(2) = - lcio_particle_flow (prt%obj, 1_c_int)
  end function lcio_particle_get_flow

  module function lcio_particle_get_momentum (prt) result (p)
    type(vector4_t) :: p
    type(lcio_particle_t), intent(in) :: prt
    real(default) :: E, px, py, pz
    E = lcio_energy (prt%obj)
    px = lcio_three_momentum (prt%obj, 0_c_int)
    py = lcio_three_momentum (prt%obj, 1_c_int)
    pz = lcio_three_momentum (prt%obj, 2_c_int)
    p = vector4_moving ( E, vector3_moving ([ px, py, pz ]))
  end function lcio_particle_get_momentum

  module function lcio_particle_get_mass_squared (prt) result (m2)
    real(default) :: m2
    type(lcio_particle_t), intent(in) :: prt
    real(default) :: m
    m = lcio_mass (prt%obj)
    m2 = sign (m**2, m)
  end function lcio_particle_get_mass_squared

  module function lcio_particle_get_vertex (prt) result (vtx)
    type(vector3_t) :: vtx
    type(lcio_particle_t), intent(in) :: prt
    real(default) :: vx, vy, vz
    vx = lcio_vtx_x (prt%obj)
    vy = lcio_vtx_y (prt%obj)
    vz = lcio_vtx_z (prt%obj)
    vtx = vector3_moving ([vx, vy, vz])
  end function lcio_particle_get_vertex

  module function lcio_particle_get_time (prt) result (time)
    real(default) :: time
    type(lcio_particle_t), intent(in) :: prt
    time = lcio_prt_time (prt%obj)
    time = time / ns_per_mm
  end function lcio_particle_get_time

  module subroutine lcio_polarization_init_pol (prt, pol)
    type(lcio_particle_t), intent(inout) :: prt
    type(polarization_t), intent(in) :: pol
    real(default) :: r, theta, phi
    if (pol%is_polarized ()) then
       call pol%to_angles (r, theta, phi)
       call lcio_particle_set_spin (prt%obj, &
            real(r, c_double), real (theta, c_double), real (phi, c_double))
    end if
  end subroutine lcio_polarization_init_pol

  module subroutine lcio_polarization_init_hel (prt, hel)
    type(lcio_particle_t), intent(inout) :: prt
    type(helicity_t), intent(in) :: hel
    integer, dimension(2) :: h
    if (hel%is_defined ()) then
       h = hel%to_pair ()
       select case (h(1))
       case (1:)
          call lcio_particle_set_spin (prt%obj, 1._c_double, &
               0._c_double, 0._c_double)
       case (:-1)
          call lcio_particle_set_spin (prt%obj, 1._c_double, &
               real (pi, c_double), 0._c_double)
       case (0)
          call lcio_particle_set_spin (prt%obj, 1._c_double, &
               real (pi/2, c_double), 0._c_double)
       end select
    end if
  end subroutine lcio_polarization_init_hel

  module subroutine lcio_polarization_init_int (prt, hel)
    type(lcio_particle_t), intent(inout) :: prt
    integer, intent(in) :: hel
    call lcio_particle_set_spin (prt%obj, 0._c_double, &
         0._c_double, real (hel, c_double))
  end subroutine lcio_polarization_init_int

  module subroutine lcio_particle_to_pol (prt, flv, pol)
    type(lcio_particle_t), intent(in) :: prt
    type(flavor_t), intent(in) :: flv
    type(polarization_t), intent(out) :: pol
    real(default) :: degree, theta, phi
    degree = lcio_polarization_degree (prt%obj)
    theta = lcio_polarization_theta (prt%obj)
    phi = lcio_polarization_phi (prt%obj)
    call pol%init_angles (flv, degree, theta, phi)
  end subroutine lcio_particle_to_pol

  module subroutine lcio_particle_to_hel (prt, flv, hel)
    type(lcio_particle_t), intent(in) :: prt
    type(flavor_t), intent(in) :: flv
    type(helicity_t), intent(out) :: hel
    real(default) :: theta
    integer :: hmax
    theta = lcio_polarization_theta (prt%obj)
    hmax = flv%get_spin_type () / 2
    call hel%init (sign (hmax, nint (cos (theta))))
  end subroutine lcio_particle_to_hel

  module subroutine lcio_particle_set_vtx (prt, vtx)
    type(lcio_particle_t), intent(inout) :: prt
    type(vector3_t), intent(in) :: vtx
    call lcio_particle_set_vertex (prt%obj, real(vtx%p(1), c_double), &
         real(vtx%p(2), c_double), real(vtx%p(3), c_double))
  end subroutine lcio_particle_set_vtx

  module subroutine lcio_particle_set_t (prt, t)
    type(lcio_particle_t), intent(inout) :: prt
    real(default), intent(in) :: t
    real(default) :: ns_from_t_mm
    ns_from_t_mm = ns_per_mm * t
    call lcio_particle_set_time (prt%obj, real(ns_from_t_mm, c_float))
  end subroutine lcio_particle_set_t

  module subroutine lcio_particle_set_parent (daughter, parent)
    type(lcio_particle_t), intent(inout) :: daughter, parent
    call lcio_particle_add_parent (daughter%obj, parent%obj)
  end subroutine lcio_particle_set_parent

  module function lcio_particle_get_status (lptr) result (status)
    integer :: status
    type(lcio_particle_t), intent(in) :: lptr
    status = lcio_particle_get_generator_status (lptr%obj)
  end function lcio_particle_get_status

  module function lcio_particle_get_pdg (lptr) result (pdg)
    integer :: pdg
    type(lcio_particle_t), intent(in) :: lptr
    pdg = lcio_particle_get_pdg_code (lptr%obj)
  end function lcio_particle_get_pdg

  module function lcio_particle_get_n_parents (lptr) result (n_parents)
    integer :: n_parents
    type(lcio_particle_t), intent(in) :: lptr
    n_parents = lcio_n_parents (lptr%obj)
  end function lcio_particle_get_n_parents

  module function lcio_particle_get_n_children (lptr) result (n_children)
    integer :: n_children
    type(lcio_particle_t), intent(in) :: lptr
    n_children = lcio_n_daughters (lptr%obj)
  end function lcio_particle_get_n_children

  module function lcio_get_n_parents &
       (evt, num_part, k_parent) result (index_parent)
    type(lcio_event_t), intent(in) :: evt
    integer, intent(in) :: num_part, k_parent
    integer :: index_parent
    index_parent = lcio_event_parent_k (evt%obj, int (num_part, c_int), &
         int (k_parent, c_int))
  end function lcio_get_n_parents

  module function lcio_get_n_children &
       (evt, num_part, k_daughter) result (index_daughter)
    type(lcio_event_t), intent(in) :: evt
    integer, intent(in) :: num_part, k_daughter
    integer :: index_daughter
    index_daughter = lcio_event_daughter_k (evt%obj, int (num_part, c_int), &
         int (k_daughter, c_int))
  end function lcio_get_n_children

  module subroutine lcio_writer_open_out (lcio_writer, filename)
    type(lcio_writer_t), intent(out) :: lcio_writer
    type(string_t), intent(in) :: filename
    lcio_writer%obj = open_lcio_writer_new (char (filename) // &
         c_null_char, 9_c_int)
  end subroutine lcio_writer_open_out

  module subroutine lcio_writer_close (lciowriter)
    type(lcio_writer_t), intent(inout) :: lciowriter
    call lcio_writer_delete (lciowriter%obj)
  end subroutine lcio_writer_close

  module subroutine lcio_event_write (wrt, evt)
    type(lcio_writer_t), intent(inout) :: wrt
    type(lcio_event_t), intent(in) :: evt
    call lcio_write_event (wrt%obj, evt%obj)
  end subroutine lcio_event_write

  module subroutine lcio_open_file (lcio_reader, filename)
    type(lcio_reader_t), intent(out) :: lcio_reader
    type(string_t), intent(in) :: filename
    lcio_reader%obj = open_lcio_reader (char (filename) // c_null_char)
  end subroutine lcio_open_file

  module subroutine lcio_reader_close (lcioreader)
    type(lcio_reader_t), intent(inout) :: lcioreader
    call lcio_reader_delete (lcioreader%obj)
  end subroutine lcio_reader_close

  module subroutine lcio_read_event (lcrdr, evt, ok)
    type(lcio_reader_t), intent(inout) :: lcrdr
    type(lcio_event_t), intent(out) :: evt
    logical, intent(out) :: ok
    evt%obj = read_lcio_event (lcrdr%obj)
    ok = c_associated (evt%obj)
  end subroutine lcio_read_event

  module function lcio_event_get_event_index (evt) result (i_evt)
    integer :: i_evt
    type(lcio_event_t), intent(in) :: evt
    i_evt = lcio_event_get_event_number (evt%obj)
  end function lcio_event_get_event_index

  module function lcio_event_get_process_id (evt) result (i_proc)
    integer :: i_proc
    type(lcio_event_t), intent(in) :: evt
    i_proc = lcio_event_signal_process_id (evt%obj)
  end function lcio_event_get_process_id

  module function lcio_event_get_n_tot (evt) result (n_tot)
    integer :: n_tot
    type(lcio_event_t), intent(in) :: evt
    n_tot = lcio_event_get_n_particles (evt%obj)
  end function lcio_event_get_n_tot

  module function lcio_event_get_alphas (evt) result (as)
    type(lcio_event_t), intent(in) :: evt
    real(default) :: as
    as = lcio_event_get_alpha_qcd (evt%obj)
  end function lcio_event_get_alphas

  module function lcio_event_get_scaleval (evt) result (scale)
    type(lcio_event_t), intent(in) :: evt
    real(default) :: scale
    scale = lcio_event_get_scale (evt%obj)
  end function lcio_event_get_scaleval

  module function lcio_event_get_particle (evt, n) result (prt)
    type(lcio_event_t), intent(in) :: evt
    integer, intent(in) :: n
    type(lcio_particle_t) :: prt
    prt%obj = lcio_event_particle_k (evt%obj, int (n, c_int))
  end function lcio_event_get_particle


end submodule lcio_interface_s

