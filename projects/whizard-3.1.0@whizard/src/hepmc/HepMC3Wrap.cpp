//////////////////////////////////////////////////////////////////////////
// Interface for building HEPMC events
//////////////////////////////////////////////////////////////////////////

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/Print.h"
#include "HepMC3/Units.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/ReaderHEPEVT.h"
#include "HepMC3/WriterHEPEVT.h"
#include "HepMC3_WHIZARD_Polarization.h"

using namespace HepMC3;

// Tell the caller that this is the true HepMC library
extern "C" bool hepmc_available() {
  return true;
}

//////////////////////////////////////////////////////////////////////////
// Polarization functions (no longer existent (yet) in HepMC3

Polarization::Polarization( )
  : m_theta( 0. ),
    m_phi( 0. ),
    m_defined( false )
{ }

Polarization::Polarization( double theta, double phi )
  : m_theta( valid_theta(theta) ),
    m_phi  ( valid_phi(phi) ),
    m_defined( true )
{ }

double Polarization::valid_theta( double theta ) {
  // this is just absolute value.
  theta = ( theta>0 ? theta : -theta );
  // translate to 0 < theta < 2pi
  theta = ( theta/(2*HepMC_pi) - int(theta/(2*HepMC_pi)) ) 
    * 2*HepMC_pi;
  // now translate to 0 < theta < pi
  if ( theta > HepMC_pi ) theta = 2*HepMC_pi - theta;
  return theta;
}

double Polarization::valid_phi( double phi ) {
  //
  // translate to -2pi < phi < 2pi
  phi = ( phi/(2*HepMC_pi) - int(phi/(2*HepMC_pi)) ) * 2*HepMC_pi;
  // translates to 0 < phi < 2pi
  if ( phi < 0 ) phi = 2*HepMC_pi + phi;
  return phi;
}

//////////////////////////////////////////////////////////////////////////
// GenEvent functions

extern "C" GenEvent* new_gen_event( int proc_id, int event_id ) {
  GenEvent* evt = new GenEvent(Units::GEV,Units::MM);
  evt->set_event_number(event_id);
  evt->add_attribute("signal_process_id",
  		     make_shared<IntAttribute>( proc_id ));
  return evt;
}

extern "C" int gen_event_get_n_particles( GenEvent* evt) {
  return evt->particles().size();
}

extern "C" int gen_event_get_n_beams( GenEvent* evt) {
  return evt->beams().size();
}

// Deletion leads to a segmentation fault in the hepmc_interface unit test (?)
extern "C" void gen_event_delete( GenEvent* evt) {
  // delete evt;
}

// This seems also deprecated
extern "C" void gen_event_print( GenEvent* evt ) {
  Print::listing (*evt);
}

extern "C" int gen_event_event_number( GenEvent* evt ) {
  return evt->event_number();
}

// The following two are no standard attributes of the HepMC3
// event record any more
extern "C" void gen_event_set_signal_process_id( GenEvent* evt, int id ) {
  evt->add_attribute("signal_process_id",
		     make_shared<IntAttribute>( id ));
}

extern "C" int gen_event_signal_process_id( GenEvent* evt ) {
  shared_ptr<IntAttribute> A_signal_process_id =
    evt->attribute<IntAttribute>("signal_process_id");
  int signal_process_id=A_signal_process_id?(A_signal_process_id->value()):0;
  return signal_process_id;
}
 
extern "C" void gen_event_set_event_scale( GenEvent* evt, double scale ) {
  evt->add_attribute("event_scale",
		     make_shared<DoubleAttribute>( scale ));
}

extern "C" double gen_event_event_scale( GenEvent* evt) {
  shared_ptr<DoubleAttribute> A_event_scale =
    evt->attribute<DoubleAttribute>("event_scale");
  double event_scale=A_event_scale?(A_event_scale->value()):0.0;
  return event_scale;
}

extern "C" void gen_event_set_alpha_qcd( GenEvent* evt, double a ) {
  evt->add_attribute("alphaQCD",
		     make_shared<DoubleAttribute>( a ));
}

extern "C" double gen_event_alpha_qcd( GenEvent* evt) {
  shared_ptr<DoubleAttribute> A_alpha_qcd = 
    evt->attribute<DoubleAttribute>("alphaQCD");
  double alpha_qcd=A_alpha_qcd?(A_alpha_qcd->value()):0.0;
  return alpha_qcd;
}

extern "C" void gen_event_set_alpha_qed( GenEvent* evt, double a ) {
  evt->add_attribute("alphaQED",
		     make_shared<DoubleAttribute>( a ));
}

extern "C" double gen_event_alpha_qed( GenEvent* evt) {
  shared_ptr<DoubleAttribute> A_alpha_qed = 
    evt->attribute<DoubleAttribute>("alphaQED");
  double alpha_qed=A_alpha_qed?(A_alpha_qed->value()):0.0;
  return alpha_qed;
}

extern "C" void gen_event_clear_weights( GenEvent* evt ) {
  evt->weights().clear();
}

extern "C" void gen_event_add_weight( GenEvent* evt, double w ) {
  evt->weights().push_back( w );
}

extern "C" int gen_event_weights_size( GenEvent* evt ) {
  return evt->weights().size();
}

extern "C" double gen_event_weight( GenEvent* evt, int i ) {
  if (0 <= i && i <= evt->weights().size()) {
    return evt->weights()[i];
  } else {
   return 0;
  }
}

extern "C" void gen_event_add_vertex( GenEvent* evt, GenVertex* v ) {
  evt->add_vertex( v );
}

extern "C" void gen_event_set_signal_process_vertex
( GenEvent* evt, GenVertex* v ) {
  // No longer existent.
  // evt->set_signal_process_vertex( v );
}

extern "C" GenVertex* gen_event_get_signal_process_vertex
( GenEvent* evt ) {
  // No longer existent.
  // return evt->signal_process_vertex();
}

extern "C" void gen_event_set_beam_particles
( GenEvent* evt, GenParticlePtr prt1, GenParticlePtr prt2) {
  evt->set_beam_particles( prt1, prt2 );
}

extern "C" void gen_event_set_cross_section
( GenEvent* evt, double xs, double xs_err) {
  shared_ptr<GenCrossSection> xsec = make_shared<GenCrossSection>();
  xsec->set_cross_section (xs, xs_err);
  evt->set_cross_section( xsec );
}

//////////////////////////////////////////////////////////////////////////
// GenEvent particle iterator functions

extern "C" GenEvent* 
new_event_particle_const_iterator( GenEvent* evt ) {
  new GenEvent();
}

extern "C" void event_particle_const_iterator_delete
( GenEvent* it ) {
  delete it;
}

extern "C" void event_particle_const_iterator_advance
( GenEvent* it ) {
  delete it;
}

extern "C" void event_particle_const_iterator_reset
( GenEvent* it, GenEvent* evt ) {
  delete it;
}

extern "C" bool event_particle_const_iterator_is_valid
( GenEvent* it, GenEvent* evt ) {
  return it != 0;
}

extern "C" GenParticlePtr event_particle_const_iterator_get
( GenEvent* it ) {
  new GenParticlePtr;
}

//////////////////////////////////////////////////////////////////////////
// GenVertex functions

extern "C" GenVertex* new_gen_vertex() {
  return new GenVertex( FourVector::ZERO_VECTOR() );
}

extern "C" GenVertex* new_gen_vertex_pos( FourVector* pos ) {
  return new GenVertex( *pos );
}
 
extern "C" void gen_vertex_delete( GenVertex* v ) {
  delete v;
}

extern "C" void gen_vertex_add_particle_in( GenVertex* v, GenParticle* p ) {
  v->add_particle_in( p );
}

extern "C" void gen_vertex_add_particle_out( GenVertex* v, GenParticle* p ) {
  v->add_particle_out( p );
}

extern "C" bool gen_vertex_is_valid( GenVertex* v ) {
  return v != 0;
}

extern "C" int gen_vertex_particles_in_size( GenVertex* v ) {
  return v->particles_in().size();
}

extern "C" int gen_vertex_particles_out_size( GenVertex* v ) {
  return v->particles_out().size();
}

extern "C" double gen_vertex_pos_x( GenVertex* v ) {
  return v->position().x();
}

extern "C" double gen_vertex_pos_y( GenVertex* v ) {
  return v->position().y();
}

extern "C" double gen_vertex_pos_z( GenVertex* v ) {
  return v->position().z();
}

extern "C" double gen_vertex_time( GenVertex* v ) {
  return v->position().t();
}

//////////////////////////////////////////////////////////////////////////
// GenVertex iterator over in-particles
//  iterators do not exist anymore in HepMCv3
// the following are all dummy routines

extern "C" GenVertex* 
new_vertex_particles_in_const_iterator( GenVertex* v ) {
  new GenVertex();  
}

extern "C" void vertex_particles_in_const_iterator_delete
( GenVertex* it ) {
  delete it;
}

extern "C" void vertex_particles_in_const_iterator_advance
( GenVertex* it ) {
  delete it;
}

extern "C" void vertex_particles_in_const_iterator_reset
( GenVertex* it, GenVertex* v ) {
  delete it;
}

extern "C" bool vertex_particles_in_const_iterator_is_valid
( GenVertex* it, GenVertex* v ) {
  return it != 0;
}

extern "C" GenParticlePtr vertex_particles_in_const_iterator_get
( GenVertex* it ) {
  new GenParticlePtr();
}

extern "C" GenParticle* vertex_get_nth_particle_in( GenVertex* vtx, int n) {
  return vtx->particles_in()[n-1].get();
}

//////////////////////////////////////////////////////////////////////////
// GenVertex iterator over out-particles
//  iterators do not exist anymore in HepMCv3
// the following are all dummy routines

extern "C" GenVertex* 
new_vertex_particles_out_const_iterator( GenVertex* v ) {
  new GenVertex();
}

extern "C" void vertex_particles_out_const_iterator_delete
( GenVertex* it ) {
  delete it;
}

extern "C" void vertex_particles_out_const_iterator_advance
( GenVertex* it ) {
  delete it;
}

extern "C" void vertex_particles_out_const_iterator_reset
( GenVertex* it, GenVertex* v ) {
  delete it;
}

extern "C" bool vertex_particles_out_const_iterator_is_valid
( GenVertex* it, GenVertex* v ) {
  return it != 0;
}

extern "C" GenParticlePtr vertex_particles_out_const_iterator_get
( GenVertex* it ) {
  new GenParticlePtr();
}

extern "C" GenParticle* vertex_get_nth_particle_out( GenVertex* vtx, int n) {
  return vtx->particles_out()[n-1].get();
}

extern "C" GenParticle* gen_event_get_nth_particle( GenEvent* evt, int n) {
  return evt->particles()[n-1].get();
}

extern "C" int gen_event_get_nth_beam( GenEvent* evt, int n) {
  return evt->beams()[n-1].get()->id();
}

//////////////////////////////////////////////////////////////////////////
// GenParticle functions
 
extern "C" GenParticle* new_gen_particle
(FourVector* momentum, int pdg_id, int status) {
  return new GenParticle( *momentum, pdg_id, status );
}
   
extern "C" void gen_particle_delete( GenParticle* prt ) {
  delete prt;
}

extern "C" void gen_particle_set_flow
( GenParticle* prt, int code_index, int code ) {
  prt->add_attribute("flow"+to_string(code_index), make_shared<IntAttribute>( code ));
}

extern "C" void gen_particle_set_polarization
( GenParticle* prt, Polarization* pol) {
  double theta = pol->theta ();
  double phi = pol->phi ();
  prt->add_attribute("theta",
		     make_shared<DoubleAttribute>( theta ));
  prt->add_attribute("phi",
		     make_shared<DoubleAttribute>( phi ));
}

extern "C" int gen_particle_barcode( GenParticle* prt ) {
  return prt->id();
}

extern "C" FourVector* gen_particle_momentum( GenParticle* prt ) {
  return new FourVector( prt->momentum() );
}

extern "C" double gen_particle_generated_mass( GenParticle* prt ) {
  return prt->generated_mass();
}

extern "C" int gen_particle_pdg_id( GenParticle* prt ) {
  return prt->pdg_id();
}

extern "C" int gen_particle_get_n_children( GenParticle* prt ) {
  return prt->children().size();
}

extern "C" int gen_particle_get_n_parents( GenParticle* prt ) {
  return prt->parents().size();
}

extern "C" int gen_particle_status( GenParticle* prt ) {
  return prt->status();
}

// No longer exists. Return false for now. 

extern "C" bool gen_particle_is_beam( GenParticle* prt ) {
  return false;
}

extern "C" GenVertex* gen_particle_production_vertex( GenParticle* prt ) {
  return prt->production_vertex().get();
}

extern "C" GenVertex* gen_particle_end_vertex( GenParticle* prt ) {
  return prt->end_vertex().get();
}

extern "C" Polarization* gen_particle_polarization( GenParticle* prt ) {
  shared_ptr<DoubleAttribute> A_theta =
    prt->attribute<DoubleAttribute>("theta");
  shared_ptr<DoubleAttribute> A_phi =
    prt->attribute<DoubleAttribute>("phi");
  double theta=A_theta?(A_theta->value()):0.0;
  double phi=A_phi?(A_phi->value()):0.0;
  return new Polarization(theta, phi);
}

extern "C" int gen_particle_flow( GenParticle* prt, int code_index ) {
  shared_ptr<IntAttribute> A_flow =
    prt->attribute<IntAttribute>("flow"+std::to_string(code_index));
  int flow=A_flow?(A_flow->value()):0;
  return flow;
}

//////////////////////////////////////////////////////////////////////////
// FourVector functions

extern "C" FourVector* new_four_vector_xyzt
( double x, double y, double z, double t) {
  return new FourVector( x, y, z, t);
}
 
extern "C" FourVector* new_four_vector_xyz( double x, double y, double z) {
  return new FourVector( x, y, z, 0);
}

extern "C" void four_vector_delete( FourVector* p ) {
  delete p;
}

extern "C" double four_vector_px( FourVector* p ) {
  return p->px();
}

extern "C" double four_vector_py( FourVector* p ) {
  return p->py();
}

extern "C" double four_vector_pz( FourVector* p ) {
  return p->pz();
}

extern "C" double four_vector_e( FourVector* p ) {
  return p->e();
}

//////////////////////////////////////////////////////////////////////////
// Polarization functions: the Polarization type does not exist in HepMC3!
// We inserted a stub here

extern "C" Polarization* new_polarization( double theta, double phi ) {
  return new Polarization( theta, phi );
}

extern "C" void polarization_delete( Polarization* pol ) {
  delete pol;
}

extern "C" double polarization_theta( Polarization* pol ) {
  return pol->theta();
}

extern "C" double polarization_phi( Polarization* pol ) {
  return pol->phi();
}
//////////////////////////////////////////////////////////////////////////
/// // IO_GenEvent functions

extern "C" Writer* new_io_gen_event_out( int* io_format, char* filename ) {
  switch (*io_format) {
  case 1:
    return (Writer*)(new WriterAsciiHepMC2( filename ));
  case 2:
    return (Writer*)(new WriterAscii( filename));
  case 4:
    return (Writer*)(new WriterHEPEVT( filename));
  default:
    return (Writer*)(new WriterAscii( filename));
  }
}

extern "C" Reader* new_io_gen_event_in( int* io_format, char* filename ) {
  switch (*io_format) {
  case 1:
    return (Reader*)(new ReaderAsciiHepMC2( filename));
  case 2:
    return (Reader*)(new ReaderAscii( filename));
  case 4:
    return (Reader*)(new ReaderHEPEVT( filename));
  default:
    return (Reader*)(new ReaderAscii( filename));
  }
}

extern "C" void io_gen_event_delete( Writer* iostream ) {
  iostream->close();
  delete iostream;
}

extern "C" void io_gen_event_write_event
( Writer* writer, const GenEvent* evt) {
  writer->write_event( *evt);
}

extern "C" bool io_gen_event_read_event
( Reader* reader, GenEvent* evt) {
  bool ok;
  ok = reader->read_event( *evt);
  if (reader->failed()) {
    return false;
  }
  else {
    return ok;
  }
}


