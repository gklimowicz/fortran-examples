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
!! Dummy interface for non-existent LCIO library
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Tell the caller that this is not the true LCIO library
logical(c_bool) function lcio_available () bind(C)
  use iso_c_binding
  lcio_available = .false.
end function lcio_available

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LCEventImpl functions

! extern "C" void* new_lcio_event( int proc_id, int event_id ) {}
type(c_ptr) function new_lcio_event (proc_id, event_id, run_id) bind(C)
  use iso_c_binding
  integer(c_int), value :: proc_id, event_id, run_id
  new_lcio_event = c_null_ptr
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop              
end function new_lcio_event

! extern "C" void lcio_set_weight( LCEventImpl* evt, double wgt )
subroutine lcio_set_weight (evt_obj, weight) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: weight
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop
end subroutine lcio_set_weight

! extern "C" void lcio_set_sqme( LCEventImpl* evt, double sqme ) {
subroutine lcio_set_sqme (evt_obj, sqme) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: sqme
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop
end subroutine lcio_set_sqme

! extern "C" void lcio_set_alt_weight( LCEventImpl* evt, double wgt, int index )
subroutine lcio_set_alt_weight (evt_obj, weight, index) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: weight
  integer(c_int), value :: index
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop
end subroutine lcio_set_alt_weight

! extern "C" void lcio_set_alt_sqme( LCEventImpl* evt, double sqme, int index )
subroutine lcio_set_alt_sqme (evt_obj, sqme, index) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: sqme
  integer(c_int), value :: index
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop
end subroutine lcio_set_alt_sqme

! extern "C" void lcio_set_alpha_qcd ( LCEventImpl* evt, double alphas )
subroutine lcio_set_alpha_qcd (evt_obj, alphas) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: alphas
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_set_alpha_qcd

! extern "C" void lcio_set_scale ( LCEventImpl* evt, double scale )
subroutine lcio_set_scale (evt_obj, scale) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: scale
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_set_scale

! extern "C" void lcio_set_sqrts ( LCEventImpl* evt, double sqrts )
subroutine lcio_set_sqrts (evt_obj, sqrts) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: sqrts
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_set_sqrts

! extern "C" void lcio_set_xsec ( LCEventImpl* evt, double xsec, double xsec_err )
subroutine lcio_set_xsec (evt_obj, xsec, xsec_err) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: xsec, xsec_err
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_set_xsec

! extern "C" void lcio_set_beam ( LCEventImpl* evt, int pdg, int beam )
subroutine lcio_set_beam (evt_obj, pdg, beam) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  integer(c_int), value :: pdg, beam
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_set_beam

! extern "C" void lcio_set_pol ( LCEventImpl* evt, double pol1, double pol2 )
subroutine lcio_set_pol (evt_obj, pol1, pol2) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  real(c_double), value :: pol1, pol2
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_set_pol

! extern "C" void lcio_set_beam_file ( LCEventImpl* evt, char* file ) {
subroutine lcio_set_beam_file (evt_obj, file) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  character(len=1, kind=c_char), dimension(*), intent(in) :: file  
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_set_beam_file

! extern "C" void lcio_set_process_name ( LCEventImpl* evt, char* name ) {
subroutine lcio_set_process_name (evt_obj, name) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  character(len=1, kind=c_char), dimension(*), intent(in) :: name
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_set_process_name

! extern "C" void lcio_event_delete( void* evt) {}
subroutine lcio_event_delete (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_event_delete

! extern "C" LCEvent* read_lcio_event ( LCReader* lcRdr ) 
type (c_ptr) function read_lcio_event (io_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj
  read_lcio_event = c_null_ptr
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function read_lcio_event

! extern "C" void dump_lcio_event ( LCEventImpl* evt )
subroutine dump_lcio_event (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine dump_lcio_event

! extern "C" int lcio_get_event_number (LCEvent* evt)
integer(c_int) function lcio_event_get_event_number (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  lcio_event_get_event_number = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_event_get_event_number


! extern "C" int lcio_event_signal_process_id (LCEvent* evt)
integer(c_int) function lcio_event_signal_process_id (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  lcio_event_signal_process_id = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_event_signal_process_id
     
! extern "C" int lcio_event_get_n_particles (LCEvent* evt) 
integer(c_int) function lcio_event_get_n_particles (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  lcio_event_get_n_particles = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_event_get_n_particles

! extern "C" double lcio_event_get_alpha_qcd (LCEvent* evt) 
real(c_double) function lcio_event_get_alpha_qcd (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  lcio_event_get_alpha_qcd = 0._c_double
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_event_get_alpha_qcd
     
! extern "C" double lcio_event_get_scale (LCEvent* evt) 
real(c_double) function lcio_event_get_scale (evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  lcio_event_get_scale = 0._c_double
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_event_get_scale
     
! extern "C" void lcio_event_to_file ( LCEvent* evt, char* filename )
subroutine lcio_event_to_file (evt_obj, filename) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  character(c_char), dimension(*), intent(in) :: filename
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_event_to_file
     
! extern "C" void lcio_event_add_collection ( LCEventImpl* evt, LCCollectionVec* mcVec )     
subroutine lcio_event_add_collection (evt_obj, lccoll_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj, lccoll_obj
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_event_add_collection

! extern "C" MCParticleImpl* lcio_event_particle_k ( LCEventImpl* evt, int k ) 
type(c_ptr) function lcio_event_particle_k (evt_obj, k) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  integer(c_int), value :: k
  lcio_event_particle_k = c_null_ptr
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_event_particle_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MCParticleImpl and LCCollectionVec functions

! extern "C" LCCollectionVec* new_lccollection()
type(c_ptr) function new_lccollection () bind(C)
  use iso_c_binding
  new_lccollection = c_null_ptr
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function new_lccollection

! extern "C" void add_particle_to_collection (MCParticleImpl* mcp, LCCollectionVec* mcVec)
subroutine add_particle_to_collection (prt_obj, lccoll_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj, lccoll_obj
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine add_particle_to_collection
     
! extern "C" MCParticleImpl* new_lcio_particle(void* momentum, int pdg_id, int status)
type(c_ptr) function new_lcio_particle &
     (px, py, pz, pdg_id, mass, charge, status) bind(C)
  use iso_c_binding
  integer(c_int), value :: pdg_id, status
  real(c_double), value :: px, py, pz, mass, charge
  new_lcio_particle = c_null_ptr
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function new_lcio_particle

! extern "C" void lcio_particle_add_parent
subroutine lcio_particle_add_parent (io_obj1, io_obj2) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj1, io_obj2
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_particle_add_parent

! extern "C" MCParticleImpl* lcio_set_color_flow
subroutine lcio_set_color_flow (prt_obj, cflow1, cflow2) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  integer(c_int), value :: cflow1, cflow2
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_set_color_flow

! extern "C" MCParticleImpl* lcio_particle_set_spin
!      (MCParticleImpl* mcp, int s1, int s2, int s3)
subroutine lcio_particle_set_spin (prt_obj, s1, s2, s3) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  real(c_double), value :: s1, s2, s3
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_particle_set_spin

! extern "C" MCParticleImpl* lcio_particle_set_time 
subroutine lcio_particle_set_time (prt_obj, t) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  real(c_double), value :: t
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_particle_set_time

! extern "C" MCParticleImpl* lcio_particle_set_vertex
! (MCParticleImpl* mcp, const double vx, const double vy, const double vz) 
subroutine lcio_particle_set_vertex (prt_obj, vx, vy, vz) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  real(c_double), value :: vx, vy, vz
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_particle_set_vertex

! extern "C" double lcio_polarization_degree ( MCParticleImpl* mcp)
real(c_double) function lcio_polarization_degree (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_polarization_degree = 0._c_double
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_polarization_degree

! extern "C" double lcio_polarization_theta ( MCParticleImpl* mcp)
real(c_double) function lcio_polarization_theta (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_polarization_theta = 0._c_double
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_polarization_theta

! extern "C" double lcio_polarization_phi ( MCParticleImpl* mcp)
real(c_double) function lcio_polarization_phi (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_polarization_phi = 0._c_double
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_polarization_phi
     
! extern "C" const int* lcio_particle_flow
integer(c_int) function lcio_particle_flow (evt_obj, col_index) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  integer(c_int) :: col_index
  lcio_particle_flow = 0_c_int
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_particle_flow

! extern "C" int lcio_particle_get_generator_status ( MCParticleImpl* mcp)
integer(c_int) function lcio_particle_get_generator_status (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_particle_get_generator_status = 0_c_int
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_particle_get_generator_status

! extern "C" int lcio_particle_get_pdg_code ( MCParticleImpl* mcp)
integer(c_int) function lcio_particle_get_pdg_code (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_particle_get_pdg_code = 0_c_int
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_particle_get_pdg_code

! extern "C" double lcio_three_momentum ( MCParticleImpl* mcp, int p_index ) {
real(c_double) function lcio_three_momentum (prt_obj, p_index) bind (C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  integer(c_int), value :: p_index
  lcio_three_momentum = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_three_momentum

! extern "C" double lcio_energy ( MCParticleImpl* mcp) {
real(c_double) function lcio_energy (prt_obj) bind (C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_energy = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_energy

! extern "C" double lcio_mass ( MCParticleImpl* mcp) {
real(c_double) function lcio_mass (prt_obj) bind (C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_mass = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_mass

! extern "C" int lcio_n_parents ( MCParticleImpl* mcp)
integer(c_int) function lcio_n_parents (prt_obj) bind (C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_n_parents = 0_c_int
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_n_parents
     
! extern "C" int lcio_n_daughters ( MCParticleImpl* mcp)
integer(c_int) function lcio_n_daughters (prt_obj) bind (C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_n_daughters = 0_c_int
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_n_daughters

! extern "C" double lcio_vtx_x (MCParticleImpl* mcp) 
real(c_double) function lcio_vtx_x (prt_obj) bind (C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_vtx_x = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_vtx_x

! extern "C" double lcio_vtx_y (MCParticleImpl* mcp) 
real(c_double) function lcio_vtx_y (prt_obj) bind (C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_vtx_y = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_vtx_y

! extern "C" double lcio_vtx_z (MCParticleImpl* mcp) 
real(c_double) function lcio_vtx_z (prt_obj) bind (C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_vtx_z = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_vtx_z

! extern "C" double lcio_prt_time (MCParticleImpl* mcp) {
real(c_float) function lcio_prt_time (prt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: prt_obj
  lcio_prt_time = 0
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_prt_time

! extern "C" int lcio_event_daughter_k ( LCEventImpl* evt, int num_part, int k_daughter)
integer(c_int) function lcio_event_daughter_k &
     (evt_obj, num_part, k_daughter) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  integer(c_int), value :: num_part, k_daughter
  lcio_event_daughter_k = 0_c_int
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_event_daughter_k

! extern "C" int lcio_event_parent_k ( LCEventImpl* evt, int num_part, int k_parent)
integer(c_int) function lcio_event_parent_k &
     (evt_obj, num_part, k_parent) bind(C)
  use iso_c_binding
  type(c_ptr), value :: evt_obj
  integer(c_int), value :: num_part, k_parent
  lcio_event_parent_k = 0_c_int
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_event_parent_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LCWriter functions

! extern "C" LCWriter* open_lcio_writer_new 
type (c_ptr) function open_lcio_writer_new (filename, complevel) bind(C)
  use iso_c_binding
  character(c_char), dimension(*), intent(in) :: filename
  integer(c_int), value :: complevel
  open_lcio_writer_new = c_null_ptr
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function open_lcio_writer_new

! extern "C" LCWriter* lcio_writer_delete
subroutine lcio_writer_delete (io_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_writer_delete

! extern "C" LCWriter* lcio_write_event
subroutine lcio_write_event (io_obj, evt_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj, evt_obj
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_write_event

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LCReader functions

! extern "C" LCReader* open_lcio_reader () 
type(c_ptr) function open_lcio_reader (filename) bind(C)
  use iso_c_binding
  character(c_char), dimension(*), intent(in) :: filename
  open_lcio_reader = c_null_ptr
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function open_lcio_reader

! extern "C" LCReader* open_lcio_reader () 
type(c_ptr) function open_lcio_reader_direct_access (filename) bind(C)
  use iso_c_binding
  character(c_char), dimension(*), intent(in) :: filename
  open_lcio_reader_direct_access = c_null_ptr
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function open_lcio_reader_direct_access

! extern "C" void lcio_get_n_runs ( LCReader* lcRdr ) 
integer(c_int) function lcio_get_n_runs (io_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj
  lcio_get_n_runs = 0_c_int
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_get_n_runs

! extern "C" void lcio_get_n_events ( LCReader* lcRdr ) 
integer(c_int) function lcio_get_n_events (io_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj
  lcio_get_n_events = 0_c_int
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function lcio_get_n_events

! extern "C" void lcio_reader_delete ( LCReader* lcRdr ) 
subroutine lcio_reader_delete (io_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: io_obj
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine lcio_reader_delete
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LCRunHeader functions

! extern "C" LCRunHeaderImpl* new_lcio_run_header( int run_id ) {
type(c_ptr) function new_lcio_run_header (run_id) bind(C)
  use iso_c_binding
  integer(c_int), value :: run_id
  new_lcio_run_header = c_null_ptr
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end function new_lcio_run_header

! extern "C" void run_header_set_simstring (LCRunHeaderImpl* runHdr, char* simstring)      
subroutine run_header_set_simstring (runhdr_obj, simstring) bind(C)
  use iso_c_binding
  type(c_ptr), value :: runhdr_obj
  character(c_char), dimension(*), intent(in) :: simstring
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine run_header_set_simstring

! extern "C" void dump_run_header ( LCRunHeaderImpl* runHdr )
subroutine dump_run_header (runhdr_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: runhdr_obj
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine dump_run_header

! extern "C" void write_run_header (LCWriter* lcWrt, const LCRunHeaderImpl* runHdr)
subroutine write_run_header (lcwrt_obj, runhdr_obj) bind(C)
  use iso_c_binding
  type(c_ptr), value :: lcwrt_obj, runhdr_obj
  write (0, "(A)")  "***********************************************************"
  write (0, "(A)")  "*** LCIO: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "***********************************************************"
  stop                     
end subroutine write_run_header
