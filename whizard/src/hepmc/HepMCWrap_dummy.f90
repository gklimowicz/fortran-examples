! WHIZARD <<Version>> <<Date>>
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
!! Dummy interface for non-existent HEPMC library
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Tell the caller that this is not the true HepMC library
logical(c_bool) function hepmc_available () bind(C)
  use iso_c_binding
  hepmc_available = .false.
end function hepmc_available

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GenEvent functions

! extern "C" void* new_gen_event( int proc_id, int event_id ) {}
type(c_ptr) function new_gen_event (proc_id, event_id) bind(C)
  use iso_c_binding
  integer(c_int), value :: proc_id, event_id
  new_gen_event = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_gen_event

! extern "C" int gen_event_get_n_particles( GenEvent* evt) {}
integer(c_int) function gen_event_get_n_particles (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  gen_event_get_n_particles = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_get_n_particles

! extern "C" int gen_event_get_n_beams( GenEvent* evt) {}
integer(c_int) function gen_event_get_n_beams (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  gen_event_get_n_beams = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_get_n_beams

! extern "C" GenParticlePtr gen_event_get_nth_particle( GenEvent* evt, int n) {
type(c_ptr) function gen_event_get_nth_particle( evt_obj, n_part) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  integer(c_int), value :: n_part
  gen_event_get_nth_particle = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_get_nth_particle

! extern "C" int gen_event_get_nth_beam( GenEvent* evt, int n) {
integer(c_int) function gen_event_get_nth_beam( evt_obj, n_beam) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  integer(c_int), value :: n_beam
  gen_event_get_nth_beam = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_get_nth_beam

! extern "C" void gen_event_delete( void* evt) {}
subroutine gen_event_delete (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_delete

! extern "C" void gen_event_print( void* evt ) {}
subroutine gen_event_print (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_print

! extern "C" int gen_event_event_number( GenEvent* evt ) {}
integer(c_int) function gen_event_event_number (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  gen_event_event_number = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_event_number

! extern "C" void gen_event_set_signal_process_id( GenEvent* evt, int id ) {}
subroutine gen_event_set_signal_process_id (evt_obj, id) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  integer(c_int), value :: id
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_set_signal_process_id

! extern "C" int gen_event_signal_process_id( GenEvent* evt ) {}
integer(c_int) function gen_event_signal_process_id (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  gen_event_signal_process_id = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_signal_process_id

! extern "C" void gen_event_set_event_scale( GenEvent* evt, double scale ) {}
subroutine gen_event_set_event_scale (evt_obj, scale) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: scale
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_set_event_scale

! extern "C" double gen_event_event_scale( GenEvent* evt) {}
real(c_double) function gen_event_event_scale (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  gen_event_event_scale = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_event_scale

! extern "C" void gen_event_set_alpha_qcd( GenEvent* evt, double a ) {}
subroutine gen_event_set_alpha_qcd (evt_obj, a) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: a
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_set_alpha_qcd

! extern "C" double gen_event_alpha_qcd( GenEvent* evt) {}
real(c_double) function gen_event_alpha_qcd (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  gen_event_alpha_qcd = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_alpha_qcd

! extern "C" void gen_event_set_alpha_qed( GenEvent* evt, double a ) {}
subroutine gen_event_set_alpha_qed (evt_obj, a) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: a
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_set_alpha_qed

! extern "C" double gen_event_alpha_qed( GenEvent* evt) {}
real(c_double) function gen_event_alpha_qed (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  gen_event_alpha_qed = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_alpha_qed

! extern "C" void gen_event_clear_weights( GenEvent* evt ) {
subroutine gen_event_clear_weights (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_clear_weights

! extern "C" void gen_event_add_weight( GenEvent* evt, double w ) {}
subroutine gen_event_add_weight (evt_obj, w) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: w
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_add_weight

! extern "C" int gen_event_weights_size( GenEvent* evt ) {}
integer(c_int) function gen_event_weights_size (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  gen_event_weights_size = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_weights_size

! extern "C" double gen_event_weight( GenEvent* evt, int i ) {}
real(c_double) function gen_event_weight (evt_obj, i) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  integer(c_int), value :: i
  gen_event_weight = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_weight

! extern "C" void gen_event_add_vertex( void* evt, void* v ) {}
subroutine gen_event_add_vertex (evt_obj, v_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  type(c_ptr), value :: v_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_add_vertex

! extern "C" void gen_event_set_signal_process_vertex( void* evt, void* v ) {}
subroutine gen_event_set_signal_process_vertex (evt_obj, v_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  type(c_ptr), value :: v_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_set_signal_process_vertex

! extern "C" GenVertex* gen_event_get_signal_process_vertex( void* evt ) {}
type(c_ptr) function gen_event_get_signal_process_vertex &
     (evt_obj) bind (C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_get_signal_process_vertex

! extern "C" bool gen_event_set_beam_particles( void* evt, void* prt1, void* prt2) {}
logical(c_bool) function gen_event_set_beam_particles &
     (evt_obj, prt1_obj, prt2_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  type(c_ptr), value :: prt1_obj, prt2_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_event_set_beam_particles

! extern "C" void gen_event_set_cross_section( GenEvent* evt, double xs, double xs_err) {}
subroutine gen_event_set_cross_section (evt_obj, xs, xs_err) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: xs, xs_err
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_event_set_cross_section

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GenEvent particle iterator functions

! extern "C" void* new_event_particle_const_iterator( void* evt )
type(c_ptr) function new_event_particle_const_iterator (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  new_event_particle_const_iterator = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_event_particle_const_iterator

! extern "C" void event_particle_const_iterator_delete( void* it ) {}
subroutine event_particle_const_iterator_delete (it_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine event_particle_const_iterator_delete

! extern "C" void event_particle_const_iterator_advance( void* it ) {}
subroutine event_particle_const_iterator_advance (it_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine event_particle_const_iterator_advance

! extern "C" void event_particle_const_iterator_reset( void* it, void* evt ) {}
subroutine event_particle_const_iterator_reset (it_obj, evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj, evt_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine event_particle_const_iterator_reset

! extern "C" bool event_particle_const_iterator_is_valid( void* it, void* evt )
function event_particle_const_iterator_is_valid &
     (it_obj, evt_obj) result (flag) bind(C)
  use iso_c_binding
  logical(c_bool) :: flag
  type(c_ptr), value :: it_obj, evt_obj
  flag = .false.
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function event_particle_const_iterator_is_valid

! extern "C" void* event_particle_const_iterator_get( void* it )
type(c_ptr) function event_particle_const_iterator_get (it_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj
  event_particle_const_iterator_get = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function event_particle_const_iterator_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GenVertex functions

! extern "C" void* new_gen_vertex()
type(c_ptr) function new_gen_vertex () bind(C)
  use iso_c_binding
  new_gen_vertex = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_gen_vertex

! extern "C" void new_gen_vertex_pos( void* pos ) {}
type(c_ptr) function new_gen_vertex_pos (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  new_gen_vertex_pos = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_gen_vertex_pos

! extern "C" void gen_vertex_delete( void* v ) {}
! [not implemented]

! extern "C" void gen_vertex_add_particle_in( void* v, void* p ) {}
subroutine gen_vertex_add_particle_in (v_obj, prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: v_obj, prt_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_vertex_add_particle_in

! extern "C" void gen_vertex_add_particle_out( void* v, void* p ) {}
subroutine gen_vertex_add_particle_out (v_obj, prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: v_obj, prt_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_vertex_add_particle_out

! extern "C" bool gen_vertex_is_valid( void* v )
function gen_vertex_is_valid (v_obj) result (flag) bind(C)
  use iso_c_binding
  logical(c_bool) :: flag
  type(c_ptr), value :: v_obj
  flag = .false.
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_vertex_is_valid

! extern "C" int gen_vertex_particles_in_size( void* v )
function gen_vertex_particles_in_size (v_obj) result (size) bind(C)
  use iso_c_binding
  integer(c_int) :: size
  type(c_ptr), value :: v_obj
  size = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_vertex_particles_in_size

! extern "C" int gen_vertex_particles_out_size( void* v )
function gen_vertex_particles_out_size (v_obj) result (size) bind(C)
  use iso_c_binding
  integer(c_int) :: size
  type(c_ptr), value :: v_obj
  size = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_vertex_particles_out_size

! extern "C" double gen_vertex_pos_x( GenVertex* v ) 
function gen_vertex_pos_x (v_obj) result (x) bind(C)
  use iso_c_binding
  real(c_double) :: x
  type(c_ptr), value :: v_obj
  x = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_vertex_pos_x

! extern "C" double gen_vertex_pos_y( GenVertex* v ) 
function gen_vertex_pos_y (v_obj) result (y) bind(C)
  use iso_c_binding
  real(c_double) :: y
  type(c_ptr), value :: v_obj
  y = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_vertex_pos_y

! extern "C" double gen_vertex_pos_z( GenVertex* v ) 
function gen_vertex_pos_z (v_obj) result (z) bind(C)
  use iso_c_binding
  real(c_double) :: z
  type(c_ptr), value :: v_obj
  z = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_vertex_pos_z

! extern "C" double gen_vertex_time( GenVertex* v ) 
function gen_vertex_time (v_obj) result (t) bind(C)
  use iso_c_binding
  real(c_double) :: t
  type(c_ptr), value :: v_obj
  t = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_vertex_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GenVertex iterator over in-particles

! extern "C" void* new_vertex_particles_in_const_iterator( void* v )
type(c_ptr) function &
     new_vertex_particles_in_const_iterator (v_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: v_obj
  new_vertex_particles_in_const_iterator = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_vertex_particles_in_const_iterator

! extern "C" void vertex_particles_in_const_iterator_delete( void* it ) {}
subroutine vertex_particles_in_const_iterator_delete (it_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine vertex_particles_in_const_iterator_delete

! extern "C" void vertex_particles_in_const_iterator_advance( void* it ) {}
subroutine vertex_particles_in_const_iterator_advance (it_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine vertex_particles_in_const_iterator_advance

! extern "C" void vertex_particles_in_const_iterator_reset( void* it, void* v )
subroutine vertex_particles_in_const_iterator_reset &
     (it_obj, v_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj, v_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine vertex_particles_in_const_iterator_reset

! extern "C" bool vertex_particles_in_const_iterator_is_valid
! ( void* it, void* v )
function vertex_particles_in_const_iterator_is_valid &
     (it_obj, v_obj) result (flag) bind(C)
  use iso_c_binding
  logical(c_bool) :: flag
  type(c_ptr), value :: it_obj, v_obj
  flag = .false.
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function vertex_particles_in_const_iterator_is_valid

! extern "C" void* vertex_particles_in_const_iterator_get( void* it )
type(c_ptr) function &
     vertex_particles_in_const_iterator_get (it_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj
  vertex_particles_in_const_iterator_get = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function vertex_particles_in_const_iterator_get

! extern "C" GenParticle* vertex_get_nth_particle_in( GenVertex::particles_in_const_iterator* it, int n)
type(c_ptr) function &
     vertex_get_nth_particle_in (vtx_obj, n_part) bind (C)
  use iso_c_binding
  type(c_ptr), value :: vtx_obj
  integer(c_int), value :: n_part
  vertex_get_nth_particle_in = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function vertex_get_nth_particle_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GenVertex iterator over out-particles

! extern "C" void* new_vertex_particles_out_const_iterator( void* v )
type(c_ptr) function &
     new_vertex_particles_out_const_iterator (v_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: v_obj
  new_vertex_particles_out_const_iterator = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_vertex_particles_out_const_iterator

! extern "C" void vertex_particles_out_const_iterator_delete( void* it ) {}
subroutine vertex_particles_out_const_iterator_delete (it_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine vertex_particles_out_const_iterator_delete

! extern "C" void vertex_particles_out_const_iterator_advance( void* it ) {}
subroutine vertex_particles_out_const_iterator_advance (it_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine vertex_particles_out_const_iterator_advance

! extern "C" void vertex_particles_out_const_iterator_reset
! ( void* it, void* v )
subroutine vertex_particles_out_const_iterator_reset &
     (it_obj, v_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj, v_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine vertex_particles_out_const_iterator_reset

! extern "C" bool vertex_particles_out_const_iterator_is_valid( void* )
function vertex_particles_out_const_iterator_is_valid &
     (it_obj, v_obj) result (flag) bind(C)
  use iso_c_binding
  logical(c_bool) :: flag
  type(c_ptr), value :: it_obj, v_obj
  flag = .false.
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function vertex_particles_out_const_iterator_is_valid

! extern "C" void* vertex_particles_out_const_iterator_get( void* it )
type(c_ptr) function &
     vertex_particles_out_const_iterator_get (it_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: it_obj
  vertex_particles_out_const_iterator_get = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function vertex_particles_out_const_iterator_get

! extern "C" GenParticle* vertex_get_nth_particle_out( GenVertex::particles_out_const_iterator* it, int n)
type(c_ptr) function &
     vertex_get_nth_particle_out (vtx_obj, n_part) bind (C)
  use iso_c_binding
  type(c_ptr), value :: vtx_obj
  integer(c_int), value :: n_part
  vertex_get_nth_particle_out = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function vertex_get_nth_particle_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GenParticle functions

! extern "C" void* new_gen_particle(void* momentum, int pdg_id, int status)
type(c_ptr) function new_gen_particle (prt_obj, pdg_id, status) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  integer(c_int), value :: pdg_id, status
  new_gen_particle = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_gen_particle

! extern "C" void gen_particle_delete( void* prt ) {}
! [not implemented]

! extern "C" void gen_particle_set_flow( void* prt, int code_index, int code )
subroutine gen_particle_set_flow (prt_obj, code_index, code) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  integer(c_int), value :: code_index, code
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_particle_set_flow

! extern "C" void gen_particle_set_polarization( void* prt, void* pol) {}
subroutine gen_particle_set_polarization (prt_obj, pol_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj, pol_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine gen_particle_set_polarization

! extern "C" int gen_particle_barcode( void* prt )
function gen_particle_barcode (prt_obj) result (barcode) bind(C)
  use iso_c_binding
  integer(c_int) :: barcode
  type(c_ptr), value :: prt_obj
  barcode = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_barcode

! extern "C" void* gen_particle_momentum( void* prt )
type(c_ptr) function gen_particle_momentum (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  gen_particle_momentum = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_momentum

! extern "C" double gen_particle_generated_mass( void* prt )
function gen_particle_generated_mass (prt_obj) result (mass) bind(C)
  use iso_c_binding
  real(c_double) :: mass
  type(c_ptr), value :: prt_obj
  mass = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_generated_mass

! extern "C" int gen_particle_pdg_id( void* prt )
function gen_particle_pdg_id (prt_obj) result (pdg_id) bind(C)
  use iso_c_binding
  integer(c_int) :: pdg_id
  type(c_ptr), value :: prt_obj
  pdg_id = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_pdg_id

! extern "C" int gen_particle_get_n_children( GenParticle* prt ) {}
function gen_particle_get_n_children (prt_obj) result (n_ch) bind(C)
  use iso_c_binding
  integer(c_int) :: n_ch
  type(c_ptr), value :: prt_obj
  n_ch = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_get_n_children

! extern "C" int gen_particle_get_n_parents( GenParticle* prt ) {}
function gen_particle_get_n_parents (prt_obj) result (n_p) bind(C)
  use iso_c_binding
  integer(c_int) :: n_p
  type(c_ptr), value :: prt_obj
  n_p = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_get_n_parents

! extern "C" int gen_particle_status( void* prt )
function gen_particle_status (prt_obj) result (status) bind(C)
  use iso_c_binding
  integer(c_int) :: status
  type(c_ptr), value :: prt_obj
  status = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_status

! extern "C" int gen_particle_is_beam( void* prt )
function gen_particle_is_beam (prt_obj) result (is_beam) bind(C)
  use iso_c_binding
  logical(c_bool) :: is_beam
  type(c_ptr), value :: prt_obj
  is_beam = .false.
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_is_beam

! extern "C" void* gen_particle_production_vertex( void* prt )
type(c_ptr) function gen_particle_production_vertex (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  gen_particle_production_vertex = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_production_vertex

! extern "C" void* gen_particle_end_vertex( void* prt )
type(c_ptr) function gen_particle_end_vertex (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  gen_particle_end_vertex = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_end_vertex

! extern "C" void* gen_particle_polarization( void* prt )
type(c_ptr) function gen_particle_polarization (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  gen_particle_polarization = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_polarization

! extern "C" int gen_particle_flow( void* prt, int code_index )
function gen_particle_flow (prt_obj, code_index) result (code) bind(C)
  use iso_c_binding
  integer(c_int) :: code
  type(c_ptr), value :: prt_obj
  integer(c_int), value :: code_index
  code = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function gen_particle_flow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FourVector functions

! extern "C" void* new_four_vector_xyzt
! ( double x, double y, double z, double t)
type(c_ptr) function new_four_vector_xyz (x, y, z) bind(C)
  use iso_c_binding
  real(c_double), value :: x, y, z
  new_four_vector_xyz = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_four_vector_xyz

! extern "C" void* new_four_vector_xyz( double x, double y, double z)
type(c_ptr) function new_four_vector_xyzt (x, y, z, t) bind(C)
  use iso_c_binding
  real(c_double), value :: x, y, z, t
  new_four_vector_xyzt = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_four_vector_xyzt

! extern "C" void four_vector_delete( void* p ) {}
subroutine four_vector_delete (p_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: p_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine four_vector_delete

! extern "C" double four_vector_px( void* p )
function four_vector_px (p_obj) result (px) bind(C)
  use iso_c_binding
  real(c_double) :: px
  type(c_ptr), value :: p_obj
  px = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function four_vector_px

! extern "C" double four_vector_py( void* p )
function four_vector_py (p_obj) result (py) bind(C)
  use iso_c_binding
  real(c_double) :: py
  type(c_ptr), value :: p_obj
  py = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function four_vector_py

! extern "C" double four_vector_pz( void* p )
function four_vector_pz (p_obj) result (pz) bind(C)
  use iso_c_binding
  real(c_double) :: pz
  type(c_ptr), value :: p_obj
  pz = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function four_vector_pz

! extern "C" double four_vector_e( void* p )
function four_vector_e (p_obj) result (e) bind(C)
  use iso_c_binding
  real(c_double) :: e
  type(c_ptr), value :: p_obj
  e = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function four_vector_e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Polarization functions

! extern "C" void* new_polarization( double theta, double phi )
type(c_ptr) function new_polarization (theta, phi) bind(C)
  use iso_c_binding
  real(c_double), value :: theta, phi
  new_polarization = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_polarization

! extern "C" void polarization_delete( void* pol ) {}
subroutine polarization_delete (pol_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: pol_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine polarization_delete


! extern "C" double polarization_theta( void* pol )
function polarization_theta (pol_obj) result (theta) bind(C)
  use iso_c_binding
  real(c_double) :: theta
  type(c_ptr), value :: pol_obj
  theta = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function polarization_theta

! extern "C" double polarization_phi( void* pol )
function polarization_phi (pol_obj) result (phi) bind(C)
  use iso_c_binding
  real(c_double) :: phi
  type(c_ptr), value :: pol_obj
  phi = 0
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function polarization_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! IO_GenEvent functions

! extern "C" void* new_io_gen_event_in( int* io_format, char* filename )
type(c_ptr) function new_io_gen_event_in (io_format, filename) bind(C)
  use iso_c_binding
  character(c_char), dimension(*), intent(in) :: filename
  integer(c_int), intent(in) :: io_format
  new_io_gen_event_in = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_io_gen_event_in

! extern "C" void* new_io_gen_event_out( char* filename )
type(c_ptr) function new_io_gen_event_out (filename) bind(C)
  use iso_c_binding
  character(c_char), dimension(*), intent(in) :: filename
  new_io_gen_event_out = c_null_ptr
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function new_io_gen_event_out

! extern "C" void io_gen_event_delete( void* iostream ) {}
subroutine io_gen_event_delete (io_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine io_gen_event_delete

! extern "C" void io_gen_event_write_event
! ( void* iostream, const void* evt) {}
subroutine io_gen_event_write_event (io_obj, evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj, evt_obj
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end subroutine io_gen_event_write_event

! extern "C" bool io_gen_event_read_event
! ( void* iostream, void* evt) {}
logical(c_bool) function io_gen_event_read_event (io_obj, evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj, evt_obj
  io_gen_event_read_event = .false.
  write (0, "(A)")  "************************************************************"
  write (0, "(A)")  "*** HepMC: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "************************************************************"
  stop
end function io_gen_event_read_event
