//////////////////////////////////////////////////////////////////////////
// Interface for building LCIO events
//////////////////////////////////////////////////////////////////////////
#include<stdio.h>
#include<string>
#include<iostream>
#include<fstream>

#include "lcio.h"
#include "IO/LCWriter.h"
#include "IO/LCReader.h"
#include "EVENT/LCIO.h"
#include "EVENT/MCParticle.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCRunHeaderImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/MCParticleImpl.h"
#include "IMPL/LCTOOLS.h"
#include "UTIL/LCTime.h"

using namespace std;
using namespace lcio;
using namespace IMPL;
using namespace EVENT;

// Tell the caller that this is the true LCIO library
extern "C" bool lcio_available() {
  return true;
}

//////////////////////////////////////////////////////////////////////////
// LCEventImpl functions

// The run number at the moment is set to the process ID. We add the process
// ID in addition as an event variable.

extern "C" LCEventImpl* new_lcio_event ( int proc_id, int event_id, int run_id ) {
  LCEventImpl* evt = new LCEventImpl();
  evt->setRunNumber ( run_id );
  evt->parameters().setValue("Run ID", run_id );
  evt->parameters().setValue("ProcessID", proc_id );
  evt->setEventNumber ( event_id );
  evt->parameters().setValue("Event Number", event_id );
  LCTime now;
  evt->setTimeStamp ( now.timeStamp() );
  return evt;
}

extern "C" void lcio_event_delete( LCEventImpl* evt) {
  delete evt;
}

extern "C" void lcio_set_weight( LCEventImpl* evt, double wgt ) {
  evt->setWeight ( wgt );
}

extern "C" void lcio_set_sqme( LCEventImpl* evt, double sqme ) {
#if LCIO_VERSION_GE (2, 17)
  double sqme_dble = sqme;
#else
  float sqme_dble = sqme;
#endif
  evt->parameters().setValue ( "sqme" , sqme_dble );
}

extern "C" void lcio_set_alt_weight( LCEventImpl* evt, double wgt, int index ) {
#if LCIO_VERSION_GE (2, 17)
  double weight = wgt;
#else
  float weight = wgt;
#endif
  evt->parameters().setValue ( "weight_alt"+to_string(index), weight );
}

extern "C" void lcio_set_alt_sqme( LCEventImpl* evt, double sqme, int index ) {
#if LCIO_VERSION_GE (2, 17)
  double sqme_dble = sqme;
#else
  float sqme_dble = sqme;
#endif
  evt->parameters().setValue ( "sqme_alt"+to_string(index) , sqme_dble );
}

extern "C" void lcio_set_alpha_qcd ( LCEventImpl* evt, double alphas ) {
  float alpha_qcd = alphas;
  evt->parameters().setValue ( "alphaQCD", alpha_qcd );
}

extern "C" void lcio_set_scale ( LCEventImpl* evt, double scale ) {
  float scale_f = scale;
  evt->parameters().setValue ( "scale", scale_f );
}

extern "C" void lcio_set_sqrts ( LCEventImpl* evt, double sqrts ) {
  float sqrts_f = sqrts;
  evt->parameters().setValue ( "Energy", sqrts_f );
}

extern "C" void lcio_set_xsec ( LCEventImpl* evt, double xsec, double xsec_err ) {
  float xsec_f = xsec;
  float xsec_err_f = xsec_err;
  evt->parameters().setValue ( "crossSection", xsec_f );
  evt->parameters().setValue ( "crossSectionError", xsec_err_f );
}

extern "C" void lcio_set_beam ( LCEventImpl* evt, int pdg, int beam ) {
  if (beam == 1){    
    evt->parameters().setValue ( "beamPDG1", pdg );
  }
  else if (beam == 2){
    evt->parameters().setValue ( "beamPDG2", pdg );
  }
}

extern "C" void lcio_set_pol ( LCEventImpl* evt, double pol, int beam ) {
  float pol_f = pol;
  if (beam == 1){
    evt->parameters().setValue ( "beamPol1", pol_f );
  }
  else if (beam == 2){
    evt->parameters().setValue ( "beamPol2", pol_f );
  }
}

extern "C" void lcio_set_beam_file ( LCEventImpl* evt, char* file ) {
  evt->parameters().setValue ( "BeamSpectrum", file );
}

extern "C" void lcio_set_process_name ( LCEventImpl* evt, char* name ) {
  evt->parameters().setValue ( "processName", name );
}

extern "C" EVENT::LCEvent* read_lcio_event ( IO::LCReader* lcRdr) {
  return lcRdr->readNextEvent ();
}

// dump the event to the screen

extern "C" void dump_lcio_event ( LCEvent* evt) {
  LCTOOLS::dumpEventDetailed ( evt );
}

extern "C" int lcio_event_get_event_number (LCEvent* evt) {
  return evt->getEventNumber();
}

extern "C" int lcio_event_signal_process_id (LCEvent* evt) {
  return evt->getParameters().getIntVal("ProcessID");
}

extern "C" int lcio_event_get_n_particles (LCEvent* evt) {
  LCCollection* col = evt->getCollection( LCIO::MCPARTICLE );
  int n = col->getNumberOfElements();
  return n;
}

extern "C" double lcio_event_get_alpha_qcd (LCEvent* evt) {
  double alphas;
  return alphas = evt->parameters().getFloatVal( "alphaQCD" );
}

extern "C" double lcio_event_get_scale (LCEvent* evt) {
  double scale;
  return scale = evt->parameters().getFloatVal( "scale" );
}

// Write parameters in LCIO event in ASCII form to a stream

extern "C" std::ostream& printParameters
( const EVENT::LCParameters& params, std::ofstream &out){
  StringVec intKeys ;
  int nIntParameters = params.getIntKeys( intKeys ).size() ;
  for(int i=0; i< nIntParameters ; i++ ){
    IntVec intVec ;
    params.getIntVals(  intKeys[i], intVec ) ;
    int nInt  = intVec.size()  ;   
    out << " parameter " << intKeys[i] << " [int]: " ;     
    if( nInt == 0 ){ 
      out << " [empty] " << std::endl ;
    }
    for(int j=0; j< nInt ; j++ ){
      out << intVec[j] << ", " ;
    }
    out << endl ;
  }
  StringVec floatKeys ;
  int nFloatParameters = params.getFloatKeys( floatKeys ).size() ;
  for(int i=0; i< nFloatParameters ; i++ ){
    FloatVec floatVec ;
    params.getFloatVals(  floatKeys[i], floatVec ) ;
    int nFloat  = floatVec.size()  ;   
    out << " parameter " << floatKeys[i] << " [float]: " ;
    if( nFloat == 0 ){
      out << " [empty] " << std::endl ;
    }
    for(int j=0; j< nFloat ; j++ ){
      out << floatVec[j] << ", " ;
    }
    out << endl ;
  }
#if LCIO_VERSION_GE (2, 17)
  StringVec doubleKeys ;
  int nDoubleParameters = params.getDoubleKeys( doubleKeys ).size() ;
  for(int i=0; i< nDoubleParameters ; i++ ){
    DoubleVec doubleVec ;
    params.getDoubleVals(  doubleKeys[i], doubleVec ) ;
    int nDouble  = doubleVec.size()  ;
    out << " parameter " << doubleKeys[i] << " [double]: " ;
    if( nDouble == 0 ){
      out << " [empty] " << std::endl ;
    }
    for(int j=0; j< nDouble ; j++ ){
      out << doubleVec[j] << ", " ;
    }
    out << endl ;
  }
#endif
  StringVec stringKeys ;
  int nStringParameters = params.getStringKeys( stringKeys ).size() ;
  for(int i=0; i< nStringParameters ; i++ ){
    StringVec stringVec ;
    params.getStringVals(  stringKeys[i], stringVec ) ;
    int nString  = stringVec.size()  ;   
    out << " parameter " << stringKeys[i] << " [string]: " ; 
    if( nString == 0 ){ 
      out << " [empty] " << std::endl ;
            }
            for(int j=0; j< nString ; j++ ){
                out << stringVec[j] << ", " ;
            }
            out << endl ;
        }
  return out;
    }		

// Write MCParticles as ASCII to stream

extern "C" std::ostream& printMCParticles
(const EVENT::LCCollection* col,  std::ofstream &out) {
        out << endl 
            << "--------------- " << "print out of "  << LCIO::MCPARTICLE
	    << " collection " << "--------------- " << endl ;
        out << endl 
            << "  flag:  0x" << hex  << col->getFlag() << dec << endl ;
        printParameters( col->getParameters(), out ) ;
        int nParticles =  col->getNumberOfElements() ;
        out << "  " << LCTOOLS::getSimulatorStatusString() << std::endl ;
        // fill map with particle pointers and collection indices
        typedef std::map< MCParticle*, int > PointerToIndexMap ;
        PointerToIndexMap p2i_map ;
        for( int k=0; k<nParticles; k++){
	  MCParticle* part =  static_cast<MCParticle*>( col->getElementAt( k ) ) ;
	  p2i_map[ part ] = k ; 
        }	
        out << endl
             <<  "[   id   ]index|      PDG |    px,     py,        pz    | energy  |gen|[simstat ]| vertex x,     y   ,   z     | endpoint x,    y  ,   z     |    mass |  charge |            spin             | colorflow | [parents] - [daughters]"    
             << endl 
             << endl ;	 
        // loop over collection - preserve order
        for(  int index = 0 ; index < nParticles ; index++){
	  char buff[215];
	  MCParticle* part =  static_cast<MCParticle*>( col->getElementAt( index ) ) ;
#if LCIO_VERSION_GE (2, 13)
	  int part_id = part->id() - 1;
#else
	  int part_id = part->id();
#endif
	  sprintf(buff, "[%8.8d]%5d|%10d|% 1.2e,% 1.2e,% 1.2e|% 1.2e| %1d |%s|% 1.2e,% 1.2e,% 1.2e|% 1.2e,% 1.2e,% 1.2e|% 1.2e|% 1.2e|% 1.2e,% 1.2e,% 1.2e|  (%d, %d)   | [",
		  part_id, index, part->getPDG(),
	          part->getMomentum()[0], part->getMomentum()[1], 
	          part->getMomentum()[2], part->getEnergy(),
		  part->getGeneratorStatus(),
		  LCTOOLS::getSimulatorStatusString( part ).c_str(),
		  part->getVertex()[0], part->getVertex()[1], part->getVertex()[2],
		  part->getEndpoint()[0], part->getEndpoint()[1],  
		  part->getEndpoint()[2], part->getMass(), part->getCharge(),
		  part->getSpin()[0], part->getSpin()[1], part->getSpin()[2],
		  part->getColorFlow()[0], part->getColorFlow()[1] );
	  out  << buff;	  
	  for(unsigned int k=0;k<part->getParents().size();k++){
	    if(k>0) out << "," ;
	    out << p2i_map[ part->getParents()[k] ]  ;
	  }
	  out << "] - [" ;
	  for(unsigned int k=0;k<part->getDaughters().size();k++){
	    if(k>0) out << "," ;
	    out << p2i_map[ part->getDaughters()[k] ]  ;
	  }
	  out << "] " << endl ;	     
        }
        out << endl 
            << "-------------------------------------------------------------------------------- " 
            << endl ;
	return out;
    }
	    
// Write LCIO event to ASCII file

extern "C" void lcio_event_to_file ( LCEvent* evt, char* filename ) {
  ofstream myfile;
  myfile.open ( filename );
  myfile << endl
	 << "=========================================" << endl;
  myfile << " - Event  : " << evt->getEventNumber() << endl;
  myfile << " - run:  "         << evt->getRunNumber() << endl;
  myfile << " - timestamp "     << evt->getTimeStamp() << endl;
  myfile << " - weight "        << evt->getWeight() << endl;
  myfile << "=========================================" << endl;    
  LCTime evtTime( evt->getTimeStamp() ) ;
  myfile << " date:      "      << evtTime.getDateString() << endl ;     
  myfile << " detector : "      << evt->getDetectorName() << endl ;
  myfile << " event parameters: " << endl ;  
  printParameters (evt->getParameters(), myfile );
  const std::vector< std::string >* strVec = evt->getCollectionNames() ;
  // loop over all collections:
  std::vector< std::string >::const_iterator name ;
  for( name = strVec->begin() ; name != strVec->end() ; name++){
    LCCollection* col = evt->getCollection( *name ) ;
    myfile << endl
	   << " collection name : " << *name 
	   << endl 
	   << " parameters: " << endl ;
            // call the detailed print functions depending on type name
    if( evt->getCollection( *name )->getTypeName() == LCIO::MCPARTICLE ){            
      if( col->getTypeName() != LCIO::MCPARTICLE ){	  
	  myfile << " collection not of type " << LCIO::MCPARTICLE << endl ;
	  return ;
        }
      printMCParticles (col, myfile);
    }
    myfile.close();
  }
}

// add collection to LCIO event

extern "C" void lcio_event_add_collection
( LCEventImpl* evt, LCCollectionVec* mcVec ) {
  evt->addCollection( mcVec, LCIO::MCPARTICLE );
}

extern "C" MCParticle* lcio_event_particle_k ( LCEventImpl* evt, int k ) {
  LCCollection* col = evt->getCollection( LCIO::MCPARTICLE );  
  MCParticle* mcp =  static_cast<MCParticle*>(col->getElementAt ( k ));
  return mcp;
}

// returns the index of the parent / daughter with incoming index

extern "C" int lcio_event_parent_k
( LCEventImpl* evt, int num_part, int k_parent) {
  LCCollection* col = evt->getCollection( LCIO::MCPARTICLE );
  int nParticles = col->getNumberOfElements() ;
  std::vector<int> p_parents[nParticles];
  typedef std::map< MCParticle*, int > PointerToIndexMap ;
  PointerToIndexMap p2i_map ;
  for( int k=0; k<nParticles; k++){
    MCParticle* part = static_cast<MCParticle*>( col->getElementAt ( k ) );
    p2i_map[ part ] = k;
  }
  for( int index = 0; index < nParticles ; index++){
    MCParticle* part = static_cast<MCParticle*>( col->getElementAt ( index ) );    
    for(unsigned int k =0;k<part->getParents().size();k++){
      p_parents[index].push_back( p2i_map[ part -> getParents()[k] ]) ;
    }
  }
  return p_parents[num_part-1][k_parent-1] + 1;
}

extern "C" int lcio_event_daughter_k
( LCEventImpl* evt, int num_part, int k_daughter) {
  LCCollection* col = evt->getCollection( LCIO::MCPARTICLE );
  int nParticles = col->getNumberOfElements() ;
  std::vector<int> p_daughters[nParticles];
  typedef std::map< MCParticle*, int > PointerToIndexMap ;
  PointerToIndexMap p2i_map ;
  for( int k=0; k<nParticles; k++){
    MCParticle* part = static_cast<MCParticle*>( col->getElementAt ( k ) );
    p2i_map[ part ] = k;
  }
  for( int index = 0; index < nParticles ; index++){
    MCParticle* part = static_cast<MCParticle*>( col->getElementAt ( index ) );    
    for(unsigned int k =0;k<part->getDaughters().size();k++){
      p_daughters[index].push_back( p2i_map[ part -> getDaughters()[k] ]) ;
    }
  }
  return p_daughters[num_part-1][k_daughter-1] + 1;
}

//////////////////////////////////////////////////////////////////////////
// MCParticle and LCCollectionVec functions

extern "C" LCCollectionVec* new_lccollection() {
  LCCollectionVec* mcVec = new LCCollectionVec(LCIO::MCPARTICLE);
  return mcVec;
}

extern "C" void add_particle_to_collection 
(MCParticleImpl* mcp, LCCollectionVec* mcVec) {
  mcVec->push_back( mcp );
  
}

extern "C" MCParticleImpl* new_lcio_particle 
(double px, double py, double pz, int pdg, double mass, double charge, int status) {
  MCParticleImpl* mcp = new MCParticleImpl() ;
  double p[3] =  { px, py, pz };
  mcp->setPDG ( pdg );
  mcp->setMomentum ( p );
  mcp->setMass ( mass );
  mcp->setCharge ( charge );
  mcp->setGeneratorStatus ( status );
  mcp->setCreatedInSimulation (false);
  return mcp; 
}

extern "C" MCParticleImpl* lcio_set_color_flow
(MCParticleImpl* mcp, int cflow1, int cflow2) {
  int cflow[2] = { cflow1, cflow2 };
  mcp->setColorFlow ( cflow );
  return mcp;
}

extern "C" MCParticleImpl* lcio_particle_set_spin
(MCParticleImpl* mcp, const double spin1, const double spin2, const double spin3) {
  float spin1_fl = spin1;
  float spin2_fl = spin2;
  float spin3_fl = spin3;
  float spin[3] = { spin1_fl, spin2_fl, spin3_fl };
  mcp->setSpin( spin );
  return mcp;
}

extern "C" MCParticleImpl* lcio_particle_set_time
(MCParticleImpl* mcp, const double t) {
  mcp->setTime( t );
  return mcp;
}

extern "C" MCParticleImpl* lcio_particle_set_vertex
(MCParticleImpl* mcp, const double vx, const double vy, const double vz) {
  double vtx[3] = { vx, vy, vz };
  mcp->setVertex( vtx );
  return mcp;
}

extern "C" void lcio_particle_add_parent
( MCParticleImpl* daughter , MCParticleImpl* parent) {
  daughter->addParent( parent );
}

extern "C" int lcio_particle_get_generator_status ( MCParticleImpl* mcp) {
  return mcp->getGeneratorStatus();
}

extern "C" int lcio_particle_get_pdg_code ( MCParticleImpl* mcp) {
  return mcp->getPDG();
}

extern "C" int lcio_particle_flow ( MCParticleImpl* mcp, int col_index ) {
  return mcp->getColorFlow()[ col_index ];
}

extern "C" double lcio_polarization_degree ( MCParticleImpl* mcp) {
  return mcp->getSpin()[ 0 ];
}

extern "C" double lcio_polarization_theta ( MCParticleImpl* mcp) {
  return mcp->getSpin()[ 1 ];
}

extern "C" double lcio_polarization_phi ( MCParticleImpl* mcp) {
  return mcp->getSpin()[ 2 ];
}

extern "C" double lcio_three_momentum ( MCParticleImpl* mcp, int p_index ) {
  return mcp->getMomentum()[ p_index ];  
}

extern "C" double lcio_energy ( MCParticleImpl* mcp ) {
  return mcp->getEnergy();
}

extern "C" double lcio_mass ( MCParticleImpl* mcp ) {
  return mcp->getMass();
}

extern "C" int lcio_n_parents ( MCParticleImpl* mcp) {
  return mcp->getParents().size();
}

extern "C" int lcio_n_daughters ( MCParticleImpl* mcp) {
  return mcp->getDaughters().size();
}  

extern "C" double lcio_vtx_x (MCParticleImpl* mcp) {
  return mcp->getVertex()[0];
}

extern "C" double lcio_vtx_y (MCParticleImpl* mcp) {
  return mcp->getVertex()[1];
}

extern "C" double lcio_vtx_z (MCParticleImpl* mcp) {
  return mcp->getVertex()[2];
}

extern "C" float lcio_prt_time (MCParticleImpl* mcp) {
  return mcp->getTime();
}

//////////////////////////////////////////////////////////////////////////
// LCWriter functions

extern "C" LCWriter* open_lcio_writer_new 
( char* filename, int complevel ) {
  LCWriter* lcWrt = LCFactory::getInstance()->createLCWriter();
  lcWrt->setCompressionLevel (complevel);  
  lcWrt->open( filename, LCIO::WRITE_NEW );
  return lcWrt;
}

extern "C" LCWriter* open_lcio_writer_append
( char* filename ) {
  LCWriter* lcWrt = LCFactory::getInstance()->createLCWriter();    
  lcWrt->open( filename, LCIO::WRITE_APPEND );
  return lcWrt;
}

// write the event

extern "C" LCWriter* lcio_write_event
( LCWriter* lcWrt, LCEventImpl* evt) {
  lcWrt->writeEvent( evt );
  return lcWrt;
}

// destructor

extern "C" void lcio_writer_delete ( LCWriter* lcWrt ) {
  lcWrt->close();
  delete lcWrt;
}

//////////////////////////////////////////////////////////////////////////
// LCReader functions

extern "C" LCReader* open_lcio_reader ( char* filename) {
  LCReader* lcRdr = LCFactory::getInstance()->createLCReader();
  lcRdr->open ( filename );
  return lcRdr;
}

extern "C" LCReader* open_lcio_reader_direct_access ( char* filename) {
  LCReader* lcRdr = LCFactory::getInstance()->createLCReader(LCReader::directAccess);
  lcRdr->open ( filename );
  return lcRdr;
}

extern "C" int lcio_get_n_runs ( LCReader* lcRdr ) {
  return lcRdr->getNumberOfRuns();
}

extern "C" int lcio_get_n_events ( LCReader* lcRdr ) {
  return lcRdr->getNumberOfEvents();
}

extern "C" void lcio_reader_delete ( LCReader* lcRdr ) {
  lcRdr->close();
  delete lcRdr;
}

//////////////////////////////////////////////////////////////////////////
// LCRunHeader functions

// We set the process ID equal to the run number, and at it also as an
// explicit parameter.

extern "C" LCRunHeaderImpl* new_lcio_run_header( int rn ) {
  LCRunHeaderImpl* runHdr = new LCRunHeaderImpl;
  runHdr->setRunNumber (rn);
  return runHdr;
}

extern "C" void run_header_set_simstring
(LCRunHeaderImpl* runHdr, char* simstring) {
  runHdr->parameters().setValue ( "SimulationProgram", simstring );
}

extern "C" bool read_run_header ( LCReader* lcRdr , LCRunHeader* runHdr ) {
  return ((runHdr = lcRdr->readNextRunHeader ()) != 0);
}  

extern "C" void dump_run_header ( LCRunHeaderImpl* runHdr ) {
  LCTOOLS::dumpRunHeader( runHdr );
}
    
extern "C" void write_run_header 
(LCWriter* lcWrt, const LCRunHeaderImpl* runHdr) {
  lcWrt->writeRunHeader (runHdr);
}


