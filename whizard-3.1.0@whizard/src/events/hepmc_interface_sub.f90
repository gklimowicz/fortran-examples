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

submodule (hepmc_interface) hepmc_interface_s

  use constants, only: PI
  use physics_defs, only: pb_per_fb
  use system_dependencies, only: HEPMC2_AVAILABLE
  use system_dependencies, only: HEPMC3_AVAILABLE
  use diagnostics

  implicit none

contains

  module function hepmc_is_available () result (flag)
    logical :: flag
    flag = hepmc_available ()
  end function hepmc_is_available

  module subroutine hepmc_four_vector_init_v4 (pp, p)
    type(hepmc_four_vector_t), intent(out) :: pp
    type(vector4_t), intent(in) :: p
    real(default), dimension(0:3) :: pa
    pa = vector4_get_components (p)
    pp%obj = new_four_vector_xyzt &
         (real (pa(1), c_double), &
          real (pa(2), c_double), &
          real (pa(3), c_double), &
          real (pa(0), c_double))
  end subroutine hepmc_four_vector_init_v4

  module subroutine hepmc_four_vector_init_v3 (pp, p)
    type(hepmc_four_vector_t), intent(out) :: pp
    type(vector3_t), intent(in) :: p
    real(default), dimension(3) :: pa
    pa = vector3_get_components (p)
    pp%obj = new_four_vector_xyz &
         (real (pa(1), c_double), &
          real (pa(2), c_double), &
          real (pa(3), c_double))
  end subroutine hepmc_four_vector_init_v3

  module subroutine hepmc_four_vector_init_hepmc_prt (pp, prt)
    type(hepmc_four_vector_t), intent(out) :: pp
    type(hepmc_particle_t), intent(in) :: prt
    pp%obj = gen_particle_momentum (prt%obj)
  end subroutine hepmc_four_vector_init_hepmc_prt

  module subroutine hepmc_four_vector_final (p)
    type(hepmc_four_vector_t), intent(inout) :: p
    call four_vector_delete (p%obj)
  end subroutine hepmc_four_vector_final

  module subroutine hepmc_four_vector_to_vector4 (pp, p)
    type(hepmc_four_vector_t), intent(in) :: pp
    type(vector4_t), intent(out) :: p
    real(default) :: E
    real(default), dimension(3) :: p3
    E = four_vector_e (pp%obj)
    p3(1) = four_vector_px (pp%obj)
    p3(2) = four_vector_py (pp%obj)
    p3(3) = four_vector_pz (pp%obj)
    p = vector4_moving (E, vector3_moving (p3))
  end subroutine hepmc_four_vector_to_vector4

  module subroutine hepmc_polarization_init_pol (hpol, pol)
    type(hepmc_polarization_t), intent(out) :: hpol
    type(polarization_t), intent(in) :: pol
    real(default) :: r, theta, phi
    if (pol%is_polarized ()) then
       call pol%to_angles (r, theta, phi)
       if (r >= 0.5) then
          hpol%polarized = .true.
          hpol%obj = new_polarization &
               (real (theta, c_double), real (phi, c_double))
       end if
    end if
  end subroutine hepmc_polarization_init_pol

  module subroutine hepmc_polarization_init_hel (hpol, hel)
    type(hepmc_polarization_t), intent(out) :: hpol
    type(helicity_t), intent(in) :: hel
    integer, dimension(2) :: h
    if (hel%is_defined ()) then
       h = hel%to_pair ()
       select case (h(1))
       case (1:)
          hpol%polarized = .true.
          hpol%obj = new_polarization (0._c_double, 0._c_double)
       case (:-1)
          hpol%polarized = .true.
          hpol%obj = new_polarization (real (pi, c_double), 0._c_double)
       case (0)
          hpol%polarized = .true.
          hpol%obj = new_polarization (real (pi/2, c_double), 0._c_double)
       end select
    end if
  end subroutine hepmc_polarization_init_hel

  module subroutine hepmc_polarization_init_int (hpol, hel)
    type(hepmc_polarization_t), intent(out) :: hpol
    integer, intent(in) :: hel
    select case (hel)
    case (1:)
       hpol%polarized = .true.
       hpol%obj = new_polarization (0._c_double, 0._c_double)
    case (:-1)
       hpol%polarized = .true.
       hpol%obj = new_polarization (real (pi, c_double), 0._c_double)
    case (0)
       hpol%polarized = .true.
       hpol%obj = new_polarization (real (pi/2, c_double), 0._c_double)
    end select
  end subroutine hepmc_polarization_init_int

  module subroutine hepmc_polarization_final (hpol)
    type(hepmc_polarization_t), intent(inout) :: hpol
    if (hpol%polarized)  call polarization_delete (hpol%obj)
  end subroutine hepmc_polarization_final

  module subroutine hepmc_polarization_to_pol (hpol, flv, pol)
    type(hepmc_polarization_t), intent(in) :: hpol
    type(flavor_t), intent(in) :: flv
    type(polarization_t), intent(out) :: pol
    real(default) :: theta, phi
    theta = polarization_theta (hpol%obj)
    phi = polarization_phi (hpol%obj)
    call pol%init_angles (flv, 1._default, theta, phi)
  end subroutine hepmc_polarization_to_pol

  module subroutine hepmc_polarization_to_hel (hpol, flv, hel)
    type(hepmc_polarization_t), intent(in) :: hpol
    type(flavor_t), intent(in) :: flv
    type(helicity_t), intent(out) :: hel
    real(default) :: theta
    integer :: hmax
    theta = polarization_theta (hpol%obj)
    hmax = flv%get_spin_type () / 2
    call hel%init (sign (hmax, nint (cos (theta))))
  end subroutine hepmc_polarization_to_hel

  module subroutine hepmc_particle_init (prt, p, pdg, status)
    type(hepmc_particle_t), intent(out) :: prt
    type(vector4_t), intent(in) :: p
    integer, intent(in) :: pdg, status
    type(hepmc_four_vector_t) :: pp
    call hepmc_four_vector_init (pp, p)
    prt%obj = new_gen_particle (pp%obj, int (pdg, c_int), int (status, c_int))
    call hepmc_four_vector_final (pp)
  end subroutine hepmc_particle_init

  module subroutine hepmc_particle_set_color_col (prt, col)
    type(hepmc_particle_t), intent(inout) :: prt
    type(color_t), intent(in) :: col
    integer(c_int) :: c
    c = col%get_col ()
    if (c /= 0)  call gen_particle_set_flow (prt%obj, 1_c_int, c)
    c = col%get_acl ()
    if (c /= 0)  call gen_particle_set_flow (prt%obj, 2_c_int, c)
  end subroutine hepmc_particle_set_color_col

  module subroutine hepmc_particle_set_color_int (prt, col)
    type(hepmc_particle_t), intent(inout) :: prt
    integer, dimension(2), intent(in) :: col
    integer(c_int) :: c
    c = col(1)
    if (c /= 0)  call gen_particle_set_flow (prt%obj, 1_c_int, c)
    c = col(2)
    if (c /= 0)  call gen_particle_set_flow (prt%obj, 2_c_int, c)
  end subroutine hepmc_particle_set_color_int

  module subroutine hepmc_particle_set_polarization_pol (prt, pol)
    type(hepmc_particle_t), intent(inout) :: prt
    type(polarization_t), intent(in) :: pol
    type(hepmc_polarization_t) :: hpol
    call hepmc_polarization_init (hpol, pol)
    if (hpol%polarized)  call gen_particle_set_polarization (prt%obj, hpol%obj)
    call hepmc_polarization_final (hpol)
  end subroutine hepmc_particle_set_polarization_pol

  module subroutine hepmc_particle_set_polarization_hel (prt, hel)
    type(hepmc_particle_t), intent(inout) :: prt
    type(helicity_t), intent(in) :: hel
    type(hepmc_polarization_t) :: hpol
    call hepmc_polarization_init (hpol, hel)
    if (hpol%polarized)  call gen_particle_set_polarization (prt%obj, hpol%obj)
    call hepmc_polarization_final (hpol)
  end subroutine hepmc_particle_set_polarization_hel

  module subroutine hepmc_particle_set_polarization_int (prt, hel)
    type(hepmc_particle_t), intent(inout) :: prt
    integer, intent(in) :: hel
    type(hepmc_polarization_t) :: hpol
    call hepmc_polarization_init (hpol, hel)
    if (hpol%polarized)  call gen_particle_set_polarization (prt%obj, hpol%obj)
    call hepmc_polarization_final (hpol)
  end subroutine hepmc_particle_set_polarization_int

  module function hepmc_particle_get_barcode (prt) result (barcode)
    integer :: barcode
    type(hepmc_particle_t), intent(in) :: prt
    barcode = gen_particle_barcode (prt%obj)
  end function hepmc_particle_get_barcode

  module function hepmc_particle_get_momentum (prt) result (p)
    type(vector4_t) :: p
    type(hepmc_particle_t), intent(in) :: prt
    type(hepmc_four_vector_t) :: pp
    call hepmc_four_vector_init (pp, prt)
    call hepmc_four_vector_to_vector4 (pp, p)
    call hepmc_four_vector_final (pp)
  end function hepmc_particle_get_momentum

  module function hepmc_particle_get_mass_squared (prt) result (m2)
    real(default) :: m2
    type(hepmc_particle_t), intent(in) :: prt
    real(default) :: m
    m = gen_particle_generated_mass (prt%obj)
    m2 = sign (m**2, m)
  end function hepmc_particle_get_mass_squared

  module function hepmc_particle_get_pdg (prt) result (pdg)
    integer :: pdg
    type(hepmc_particle_t), intent(in) :: prt
    pdg = gen_particle_pdg_id (prt%obj)
  end function hepmc_particle_get_pdg

  module function hepmc_particle_get_status (prt) result (status)
    integer :: status
    type(hepmc_particle_t), intent(in) :: prt
    status = gen_particle_status (prt%obj)
  end function hepmc_particle_get_status

  module function hepmc_particle_is_beam (prt) result (is_beam)
    logical :: is_beam
    type(hepmc_particle_t), intent(in) :: prt
    is_beam = gen_particle_is_beam (prt%obj)
  end function hepmc_particle_is_beam

  module function hepmc_particle_get_production_vertex (prt) result (v)
    type(hepmc_vertex_t) :: v
    type(hepmc_particle_t), intent(in) :: prt
    v%obj = gen_particle_production_vertex (prt%obj)
  end function hepmc_particle_get_production_vertex

  module function hepmc_particle_get_decay_vertex (prt) result (v)
    type(hepmc_vertex_t) :: v
    type(hepmc_particle_t), intent(in) :: prt
    v%obj = gen_particle_end_vertex (prt%obj)
  end function hepmc_particle_get_decay_vertex

  module function hepmc_particle_get_parent_barcodes &
       (prt) result (parent_barcode)
    type(hepmc_particle_t), intent(in) :: prt
    integer, dimension(:), allocatable :: parent_barcode
    type(hepmc_vertex_t) :: v
    type(hepmc_vertex_particle_in_iterator_t) :: it
    integer :: i
    v = hepmc_particle_get_production_vertex (prt)
    if (hepmc_vertex_is_valid (v)) then
       allocate (parent_barcode (hepmc_vertex_get_n_in (v)))
       if (size (parent_barcode) /= 0) then
          if (HEPMC2_AVAILABLE) then
             call hepmc_vertex_particle_in_iterator_init (it, v)
             do i = 1, size (parent_barcode)
                parent_barcode(i) = hepmc_particle_get_barcode &
                     (hepmc_vertex_particle_in_iterator_get (it))
                call hepmc_vertex_particle_in_iterator_advance (it)
             end do
             call hepmc_vertex_particle_in_iterator_final (it)
          else if (HEPMC3_AVAILABLE) then
             do i = 1, size (parent_barcode)
                parent_barcode(i) = hepmc_particle_get_barcode &
                     (hepmc_vertex_get_nth_particle_in (v, i))
             end do
          end if
       end if
    else
       allocate (parent_barcode (0))
    end if
  end function hepmc_particle_get_parent_barcodes

  module function hepmc_particle_get_child_barcodes &
       (prt) result (child_barcode)
    type(hepmc_particle_t), intent(in) :: prt
    integer, dimension(:), allocatable :: child_barcode
    type(hepmc_vertex_t) :: v
    type(hepmc_vertex_particle_out_iterator_t) :: it
    integer :: i
    v = hepmc_particle_get_decay_vertex (prt)
    if (hepmc_vertex_is_valid (v)) then
       allocate (child_barcode (hepmc_vertex_get_n_out (v)))
       if (size (child_barcode) /= 0) then
          if (HEPMC2_AVAILABLE) then
             call hepmc_vertex_particle_out_iterator_init (it, v)
             do i = 1, size (child_barcode)
                child_barcode(i) = hepmc_particle_get_barcode &
                     (hepmc_vertex_particle_out_iterator_get (it))
                call hepmc_vertex_particle_out_iterator_advance (it)
             end do
             call hepmc_vertex_particle_out_iterator_final (it)
          else if (HEPMC3_AVAILABLE) then
             do i = 1, size (child_barcode)
                child_barcode(i) = hepmc_particle_get_barcode &
                     (hepmc_vertex_get_nth_particle_out (v, i))
             end do
          end if
       end if
    else
       allocate (child_barcode (0))
    end if
  end function hepmc_particle_get_child_barcodes

  module function hepmc_particle_get_polarization (prt) result (pol)
    type(hepmc_polarization_t) :: pol
    type(hepmc_particle_t), intent(in) :: prt
    pol%obj = gen_particle_polarization (prt%obj)
  end function hepmc_particle_get_polarization

  module function hepmc_particle_get_color (prt) result (col)
    integer, dimension(2) :: col
    type(hepmc_particle_t), intent(in) :: prt
    col(1) = gen_particle_flow (prt%obj, 1)
    col(2) = - gen_particle_flow (prt%obj, 2)
  end function hepmc_particle_get_color

  module function hepmc_vertex_to_vertex (vtx) result (v)
    type(hepmc_vertex_t), intent(in) :: vtx
    type(vector4_t) :: v
    real(default) :: t, vx, vy, vz
    if (hepmc_vertex_is_valid (vtx)) then
       t = gen_vertex_time (vtx%obj)
       vx = gen_vertex_pos_x (vtx%obj)
       vy = gen_vertex_pos_y (vtx%obj)
       vz = gen_vertex_pos_z (vtx%obj)
       v = vector4_moving (t, &
            vector3_moving ([vx, vy, vz]))
    end if
  end function hepmc_vertex_to_vertex

  module subroutine hepmc_vertex_init (v, x)
    type(hepmc_vertex_t), intent(out) :: v
    type(vector4_t), intent(in), optional :: x
    type(hepmc_four_vector_t) :: pos
    if (present (x)) then
       call hepmc_four_vector_init (pos, x)
       v%obj = new_gen_vertex_pos (pos%obj)
       call hepmc_four_vector_final (pos)
    else
       v%obj = new_gen_vertex ()
    end if
  end subroutine hepmc_vertex_init

  module function hepmc_vertex_is_valid (v) result (flag)
    logical :: flag
    type(hepmc_vertex_t), intent(in) :: v
    flag = gen_vertex_is_valid (v%obj)
  end function hepmc_vertex_is_valid

  module subroutine hepmc_vertex_add_particle_in (v, prt)
    type(hepmc_vertex_t), intent(inout) :: v
    type(hepmc_particle_t), intent(in) :: prt
    call gen_vertex_add_particle_in (v%obj, prt%obj)
  end subroutine hepmc_vertex_add_particle_in

  module subroutine hepmc_vertex_add_particle_out (v, prt)
    type(hepmc_vertex_t), intent(inout) :: v
    type(hepmc_particle_t), intent(in) :: prt
    call gen_vertex_add_particle_out (v%obj, prt%obj)
  end subroutine hepmc_vertex_add_particle_out

  module function hepmc_vertex_get_n_in (v) result (n_in)
    integer :: n_in
    type(hepmc_vertex_t), intent(in) :: v
    n_in = gen_vertex_particles_in_size (v%obj)
  end function hepmc_vertex_get_n_in

  module function hepmc_vertex_get_n_out (v) result (n_out)
    integer :: n_out
    type(hepmc_vertex_t), intent(in) :: v
    n_out = gen_vertex_particles_out_size (v%obj)
  end function hepmc_vertex_get_n_out

  module function hepmc_particle_get_parents (p) result (n_p)
    integer :: n_p
    type(hepmc_particle_t), intent(in) :: p
    n_p = gen_particle_get_n_parents (p%obj)
  end function hepmc_particle_get_parents

  module function hepmc_particle_get_children (p) result (n_ch)
    integer :: n_ch
    type(hepmc_particle_t), intent(in) :: p
    n_ch = gen_particle_get_n_children (p%obj)
  end function hepmc_particle_get_children

  module function hepmc_particle_get_n_parents (prt) result (n_parents)
    integer :: n_parents
    type(hepmc_particle_t), intent(in) :: prt
    type(hepmc_vertex_t) :: v
    if (HEPMC2_AVAILABLE) then
       v = hepmc_particle_get_production_vertex (prt)
       if (hepmc_vertex_is_valid (v)) then
          n_parents = hepmc_vertex_get_n_in (v)
       else
          n_parents = 0
       end if
    else if (HEPMC3_AVAILABLE) then
       n_parents = hepmc_particle_get_parents (prt)
    end if
  end function hepmc_particle_get_n_parents

  module function hepmc_particle_get_n_children (prt) result (n_children)
    integer :: n_children
    type(hepmc_particle_t), intent(in) :: prt
    type(hepmc_vertex_t) :: v
    if (HEPMC2_AVAILABLE) then
       v = hepmc_particle_get_decay_vertex (prt)
       if (hepmc_vertex_is_valid (v)) then
          n_children = hepmc_vertex_get_n_out (v)
       else
          n_children = 0
       end if
    else if (HEPMC3_AVAILABLE) then
       n_children = hepmc_particle_get_children (prt)
    end if
  end function hepmc_particle_get_n_children

  module subroutine hepmc_vertex_particle_in_iterator_init (it, v)
    type(hepmc_vertex_particle_in_iterator_t), intent(out) :: it
    type(hepmc_vertex_t), intent(in) :: v
    it%obj = new_vertex_particles_in_const_iterator (v%obj)
    it%v_obj = v%obj
  end subroutine hepmc_vertex_particle_in_iterator_init

  module subroutine hepmc_vertex_particle_in_iterator_final (it)
    type(hepmc_vertex_particle_in_iterator_t), intent(inout) :: it
    call vertex_particles_in_const_iterator_delete (it%obj)
  end subroutine hepmc_vertex_particle_in_iterator_final

  module subroutine hepmc_vertex_particle_in_iterator_advance (it)
    type(hepmc_vertex_particle_in_iterator_t), intent(inout) :: it
    call vertex_particles_in_const_iterator_advance (it%obj)
  end subroutine hepmc_vertex_particle_in_iterator_advance

  module subroutine hepmc_vertex_particle_in_iterator_reset (it)
    type(hepmc_vertex_particle_in_iterator_t), intent(inout) :: it
    call vertex_particles_in_const_iterator_reset (it%obj, it%v_obj)
  end subroutine hepmc_vertex_particle_in_iterator_reset

  module function hepmc_vertex_particle_in_iterator_is_valid &
       (it) result (flag)
    logical :: flag
    type(hepmc_vertex_particle_in_iterator_t), intent(in) :: it
    flag = vertex_particles_in_const_iterator_is_valid (it%obj, it%v_obj)
  end function hepmc_vertex_particle_in_iterator_is_valid

  module function hepmc_vertex_particle_in_iterator_get (it) result (prt)
    type(hepmc_particle_t) :: prt
    type(hepmc_vertex_particle_in_iterator_t), intent(in) :: it
    prt%obj = vertex_particles_in_const_iterator_get (it%obj)
  end function hepmc_vertex_particle_in_iterator_get

  module function hepmc_vertex_get_nth_particle_in (vtx, n) result (prt)
    type(hepmc_particle_t) :: prt
    type(hepmc_vertex_t), intent(in) :: vtx
    integer, intent(in) :: n
    integer(c_int) :: nth
    nth = n
    prt%obj = vertex_get_nth_particle_in (vtx%obj, nth)
  end function hepmc_vertex_get_nth_particle_in

  module function hepmc_vertex_get_nth_particle_out (vtx, n) result (prt)
    type(hepmc_particle_t) :: prt
    type(hepmc_vertex_t), intent(in) :: vtx
    integer, intent(in) :: n
    integer(c_int) :: nth
    nth = n
    prt%obj = vertex_get_nth_particle_out (vtx%obj, nth)
  end function hepmc_vertex_get_nth_particle_out

  module subroutine hepmc_vertex_particle_out_iterator_init (it, v)
    type(hepmc_vertex_particle_out_iterator_t), intent(out) :: it
    type(hepmc_vertex_t), intent(in) :: v
    it%obj = new_vertex_particles_out_const_iterator (v%obj)
    it%v_obj = v%obj
  end subroutine hepmc_vertex_particle_out_iterator_init

  module subroutine hepmc_vertex_particle_out_iterator_final (it)
    type(hepmc_vertex_particle_out_iterator_t), intent(inout) :: it
    call vertex_particles_out_const_iterator_delete (it%obj)
  end subroutine hepmc_vertex_particle_out_iterator_final

  module subroutine hepmc_vertex_particle_out_iterator_advance (it)
    type(hepmc_vertex_particle_out_iterator_t), intent(inout) :: it
    call vertex_particles_out_const_iterator_advance (it%obj)
  end subroutine hepmc_vertex_particle_out_iterator_advance

  module subroutine hepmc_vertex_particle_out_iterator_reset (it)
    type(hepmc_vertex_particle_out_iterator_t), intent(inout) :: it
    call vertex_particles_out_const_iterator_reset (it%obj, it%v_obj)
  end subroutine hepmc_vertex_particle_out_iterator_reset

  module function hepmc_vertex_particle_out_iterator_is_valid &
       (it) result (flag)
    logical :: flag
    type(hepmc_vertex_particle_out_iterator_t), intent(in) :: it
    flag = vertex_particles_out_const_iterator_is_valid (it%obj, it%v_obj)
  end function hepmc_vertex_particle_out_iterator_is_valid

  module function hepmc_vertex_particle_out_iterator_get (it) result (prt)
    type(hepmc_particle_t) :: prt
    type(hepmc_vertex_particle_out_iterator_t), intent(in) :: it
    prt%obj = vertex_particles_out_const_iterator_get (it%obj)
  end function hepmc_vertex_particle_out_iterator_get

  module subroutine hepmc_event_init (evt, proc_id, event_id)
    type(hepmc_event_t), intent(out) :: evt
    integer, intent(in), optional :: proc_id, event_id
    integer(c_int) :: pid, eid
    pid = 0;  if (present (proc_id))  pid = proc_id
    eid = 0;  if (present (event_id)) eid = event_id
    evt%obj = new_gen_event (pid, eid)
  end subroutine hepmc_event_init

  module subroutine hepmc_event_final (evt)
    type(hepmc_event_t), intent(inout) :: evt
    if (c_associated (evt%obj)) then
       call gen_event_delete (evt%obj)
       evt%obj = c_null_ptr
    end if
  end subroutine hepmc_event_final

  module subroutine hepmc_event_nullify (evt)
    type(hepmc_event_t), intent(inout) :: evt
    evt%obj = c_null_ptr
  end subroutine hepmc_event_nullify

  module function hepmc_event_get_c_ptr (evt) result (p)
    type(hepmc_event_t), intent(in) :: evt
    type(c_ptr) :: p
    p = evt%obj
  end function hepmc_event_get_c_ptr

  module subroutine hepmc_event_print (evt)
    type(hepmc_event_t), intent(in) :: evt
    call gen_event_print (evt%obj)
  end subroutine hepmc_event_print

  module function hepmc_event_get_event_index (evt) result (i_proc)
    integer :: i_proc
    type(hepmc_event_t), intent(in) :: evt
    i_proc = gen_event_event_number (evt%obj)
  end function hepmc_event_get_event_index

  module function hepmc_event_get_n_particles (evt) result (n_tot)
    integer :: n_tot
    type(hepmc_event_t), intent(in) :: evt
    n_tot = gen_event_get_n_particles (evt%obj)
  end function hepmc_event_get_n_particles

  module function hepmc_event_get_n_beams (evt) result (n_tot)
    integer :: n_tot
    type(hepmc_event_t), intent(in) :: evt
    n_tot = gen_event_get_n_beams (evt%obj)
  end function hepmc_event_get_n_beams

  module subroutine hepmc_event_set_process_id (evt, proc)
    type(hepmc_event_t), intent(in) :: evt
    integer, intent(in) :: proc
    integer(c_int) :: i_proc
    i_proc = proc
    call gen_event_set_signal_process_id (evt%obj, i_proc)
  end subroutine hepmc_event_set_process_id

  module function hepmc_event_get_process_id (evt) result (i_proc)
    integer :: i_proc
    type(hepmc_event_t), intent(in) :: evt
    i_proc = gen_event_signal_process_id (evt%obj)
  end function hepmc_event_get_process_id

  module subroutine hepmc_event_set_scale (evt, scale)
    type(hepmc_event_t), intent(in) :: evt
    real(default), intent(in) :: scale
    real(c_double) :: cscale
    cscale = scale
    call gen_event_set_event_scale (evt%obj, cscale)
  end subroutine hepmc_event_set_scale

  module function hepmc_event_get_scale (evt) result (scale)
    real(default) :: scale
    type(hepmc_event_t), intent(in) :: evt
    scale = gen_event_event_scale (evt%obj)
  end function hepmc_event_get_scale

  module subroutine hepmc_event_set_alpha_qcd (evt, alpha)
    type(hepmc_event_t), intent(in) :: evt
    real(default), intent(in) :: alpha
    real(c_double) :: a
    a = alpha
    call gen_event_set_alpha_qcd (evt%obj, a)
  end subroutine hepmc_event_set_alpha_qcd

  module function hepmc_event_get_alpha_qcd (evt) result (alpha)
    real(default) :: alpha
    type(hepmc_event_t), intent(in) :: evt
    alpha = gen_event_alpha_qcd (evt%obj)
  end function hepmc_event_get_alpha_qcd

  module subroutine hepmc_event_set_alpha_qed (evt, alpha)
    type(hepmc_event_t), intent(in) :: evt
    real(default), intent(in) :: alpha
    real(c_double) :: a
    a = alpha
    call gen_event_set_alpha_qed (evt%obj, a)
  end subroutine hepmc_event_set_alpha_qed

  module function hepmc_event_get_alpha_qed (evt) result (alpha)
    real(default) :: alpha
    type(hepmc_event_t), intent(in) :: evt
    alpha = gen_event_alpha_qed (evt%obj)
  end function hepmc_event_get_alpha_qed

  module subroutine hepmc_event_clear_weights (evt)
    type(hepmc_event_t), intent(in) :: evt
    call gen_event_clear_weights (evt%obj)
  end subroutine hepmc_event_clear_weights

  module subroutine hepmc_event_add_weight (evt, weight, rescale)
    type(hepmc_event_t), intent(in) :: evt
    real(default), intent(in) :: weight
    logical, intent(in) :: rescale
    real(c_double) :: w
    if (rescale) then
       w = weight * pb_per_fb
    else
       w = weight
    end if
    call gen_event_add_weight (evt%obj, w)
  end subroutine hepmc_event_add_weight

  module function hepmc_event_get_weights_size (evt) result (n)
    integer :: n
    type(hepmc_event_t), intent(in) :: evt
    n = gen_event_weights_size (evt%obj)
  end function hepmc_event_get_weights_size

  module function hepmc_event_get_weight (evt, index, rescale) result (weight)
    real(default) :: weight
    type(hepmc_event_t), intent(in) :: evt
    integer, intent(in) :: index
    logical, intent(in) :: rescale
    integer(c_int) :: i
    i = index - 1
    if (rescale) then
       weight = gen_event_weight (evt%obj, i) / pb_per_fb
    else
       weight = gen_event_weight (evt%obj, i)
    end if
  end function hepmc_event_get_weight

  module subroutine hepmc_event_add_vertex (evt, v)
    type(hepmc_event_t), intent(inout) :: evt
    type(hepmc_vertex_t), intent(in) :: v
    call gen_event_add_vertex (evt%obj, v%obj)
  end subroutine hepmc_event_add_vertex

  module subroutine hepmc_event_set_signal_process_vertex (evt, v)
    type(hepmc_event_t), intent(inout) :: evt
    type(hepmc_vertex_t), intent(in) :: v
    call gen_event_set_signal_process_vertex (evt%obj, v%obj)
  end subroutine hepmc_event_set_signal_process_vertex

  module function hepmc_event_get_signal_process_vertex (evt) result (v)
    type(hepmc_event_t), intent(in) :: evt
    type(hepmc_vertex_t) :: v
    v%obj = gen_event_get_signal_process_vertex (evt%obj)
  end function hepmc_event_get_signal_process_vertex

  module subroutine hepmc_event_set_beam_particles (evt, prt1, prt2)
    type(hepmc_event_t), intent(inout) :: evt
    type(hepmc_particle_t), intent(in) :: prt1, prt2
    logical(c_bool) :: flag
    flag = gen_event_set_beam_particles (evt%obj, prt1%obj, prt2%obj)
  end subroutine hepmc_event_set_beam_particles

  module subroutine hepmc_event_set_cross_section (evt, xsec, xsec_err)
    type(hepmc_event_t), intent(inout) :: evt
    real(default), intent(in) :: xsec, xsec_err
    call gen_event_set_cross_section &
         (evt%obj, &
         real (xsec * 1e-3_default, c_double), &
         real (xsec_err * 1e-3_default, c_double))
  end subroutine hepmc_event_set_cross_section

  module subroutine hepmc_event_particle_iterator_init (it, evt)
    type(hepmc_event_particle_iterator_t), intent(out) :: it
    type(hepmc_event_t), intent(in) :: evt
    it%obj = new_event_particle_const_iterator (evt%obj)
    it%evt_obj = evt%obj
  end subroutine hepmc_event_particle_iterator_init

  module subroutine hepmc_event_particle_iterator_final (it)
    type(hepmc_event_particle_iterator_t), intent(inout) :: it
    call event_particle_const_iterator_delete (it%obj)
  end subroutine hepmc_event_particle_iterator_final

  module subroutine hepmc_event_particle_iterator_advance (it)
    type(hepmc_event_particle_iterator_t), intent(inout) :: it
    call event_particle_const_iterator_advance (it%obj)
  end subroutine hepmc_event_particle_iterator_advance

  module subroutine hepmc_event_particle_iterator_reset (it)
    type(hepmc_event_particle_iterator_t), intent(inout) :: it
    call event_particle_const_iterator_reset (it%obj, it%evt_obj)
  end subroutine hepmc_event_particle_iterator_reset

  module function hepmc_event_particle_iterator_is_valid (it) result (flag)
    logical :: flag
    type(hepmc_event_particle_iterator_t), intent(in) :: it
    flag = event_particle_const_iterator_is_valid (it%obj, it%evt_obj)
  end function hepmc_event_particle_iterator_is_valid

  module function hepmc_event_particle_iterator_get (it) result (prt)
    type(hepmc_particle_t) :: prt
    type(hepmc_event_particle_iterator_t), intent(in) :: it
    prt%obj = event_particle_const_iterator_get (it%obj)
  end function hepmc_event_particle_iterator_get

  module function hepmc_event_get_nth_particle (evt, n) result (prt)
    type(hepmc_particle_t) :: prt
    type(hepmc_event_t), intent(in) :: evt
    integer, intent(in) :: n
    integer :: n_tot
    integer(c_int) :: nth
    nth = n
    n_tot = gen_event_get_n_particles (evt%obj)
    if (n > n_tot .or. n < 1) then
       prt%obj = c_null_ptr
       call msg_error ("HepMC interface called for wrong particle ID.")
    else
       prt%obj = gen_event_get_nth_particle (evt%obj, nth)
    end if
  end function hepmc_event_get_nth_particle

  module function hepmc_event_get_nth_beam (evt, n) result (beam_barcode)
    integer :: beam_barcode
    type(hepmc_event_t), intent(in) :: evt
    integer, intent(in) :: n
    integer(c_int) :: bc
    bc = gen_event_get_nth_beam (evt%obj, n)
    beam_barcode = bc
  end function hepmc_event_get_nth_beam

  module subroutine hepmc_iostream_open_out (iostream, filename, hepmc3_mode)
    type(hepmc_iostream_t), intent(out) :: iostream
    type(string_t), intent(in) :: filename
    integer, intent(in) :: hepmc3_mode
    integer(c_int) :: mode
    mode = hepmc3_mode
    iostream%obj = &
         new_io_gen_event_out (mode, char (filename) // c_null_char)
  end subroutine hepmc_iostream_open_out

  module subroutine hepmc_iostream_open_in (iostream, filename, hepmc3_mode)
    type(hepmc_iostream_t), intent(out) :: iostream
    type(string_t), intent(in) :: filename
    integer, intent(in) :: hepmc3_mode
    integer(c_int) :: mode
    mode = hepmc3_mode
    iostream%obj = &
            new_io_gen_event_in (mode, char (filename) // c_null_char)
  end subroutine hepmc_iostream_open_in

  module subroutine hepmc_iostream_close (iostream)
    type(hepmc_iostream_t), intent(inout) :: iostream
    call io_gen_event_delete (iostream%obj)
  end subroutine hepmc_iostream_close

  module subroutine hepmc_iostream_write_event (iostream, evt, hepmc3_mode)
    type(hepmc_iostream_t), intent(inout) :: iostream
    type(hepmc_event_t), intent(in) :: evt
    integer, intent(in), optional :: hepmc3_mode
    integer :: mode
    mode = HEPMC3_MODE_HEPMC3
    if (present (hepmc3_mode))  mode = hepmc3_mode
    call io_gen_event_write_event (iostream%obj, evt%obj)
  end subroutine hepmc_iostream_write_event

  module subroutine hepmc_iostream_read_event (iostream, evt, ok)
    type(hepmc_iostream_t), intent(inout) :: iostream
    type(hepmc_event_t), intent(inout) :: evt
    logical, intent(out) :: ok
    ok = io_gen_event_read_event (iostream%obj, evt%obj)
  end subroutine hepmc_iostream_read_event


end submodule hepmc_interface_s

