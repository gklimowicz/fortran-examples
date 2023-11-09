#include "WOStdHepRdr.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
//#include <stdlib.h>
#include <cmath>	// for pow()

#define IDRUP_NAME "_idrup" 
#define EVTWGT_NAME "_weight" 

extern "C" void print_stdhep_to_file (char* infile, char* outfile, int maxEvt) {    
  std::ofstream myfile;
  WOUTIL::WOStdHepRdr rdr( infile );
  myfile.open ( outfile );
  rdr.printHeader ( myfile );
  int count = 0; 
  
  while( maxEvt < 0 || count < maxEvt ) {
    
    rdr.readEvent ( myfile, count );
      
    ++count ;
  }

  myfile.close() ;
  
}

namespace WOUTIL{

  WOStdHepRdr::WOStdHepRdr(const char* evfile){
    //
    //   Use Willie's reader from LELAPS, and open the input file
    //
    _reader = new WOUTIL::WOStdHep(evfile,false);
    if(_reader->getError()) {
      std::stringstream description ; 
      description << "WOStdHepRdr: no stdhep file: " << evfile << std::ends ;
    }

  }
  WOStdHepRdr::~WOStdHepRdr(){

    delete _reader ;
  }
  
  void WOStdHepRdr::printHeader(std::ostream& os ) {

      _reader->printFileHeader( os ) ;
      // }
  }


  //
  // Read an event and return a LCCollectionVec of MCParticles
  //
  int WOStdHepRdr::readEvent(std::ostream& os, int count )
  {
    double c_light = 299.792;// mm/ns
    //  Read the event, check for errors
    
    int errorcode = _reader->readEvent() ;
    
    if( errorcode != LSH_SUCCESS ){
      
      if(  errorcode != LSH_ENDOFFILE ) {
	
	std::stringstream description ; 
	description << "WOStdHepRdr::readEvent: error when reading event: " << errorcode << std::ends ;
	
	}
      
      else {
	
	return 0 ;
      }
      }
    
    //
    //  Loop over particles
    //

    int blockid = _reader->blockId();

    // user defined process  id
    long idrup = _reader->idrup() ;  

    double evtWeight = _reader->eventweight() ;
    

    if (blockid == 204) {
      os << "Block HEPRUP [content not shown]" << std::endl;
    }
    else if ( blockid == 101 || blockid == 201 ) {

      // Print out metadata information if HEPEV4 block is filled
      if( _reader->isStdHepEv4() ){

	int pid = _reader->idrup();
	float weight = _reader->eventweight();
	float scale = _reader->scale(0,0);
	float aqed = _reader->alphaQED();
	float aqcd = _reader->alphaQCD();

	os << "Process ID = " <<  pid << "   Weight = " << weight << 
	  "     Scale = " << scale << "     aQED = " << aqed << "     aQCD = " << 
	  aqcd << std::endl;
      }

      os << " Event #" << count << std::endl;

      int NHEP = _reader->nTracks();

      if (count == 0) { 
	os << "[   id   ]index|      PDG |q(el)|    px,      py,        pz    | energy  |   mass  |sta| vertex x,     y   ,   z     |              spin           | colorflow | [daughters]" << 
	  std::endl; }

      for( int IHEP=0; IHEP<NHEP; IHEP++ )
	{
	  //  Fill particle from StdHep info
	  
	  char buff[215];
	  
	  //  PDGID
	  int pdgid = _reader->pid(IHEP);
	  
	  //  charge
	  float charge = threeCharge( pdgid ) / 3.  ;
	  
	  //  Momentum vector
	  float p0[3] = {static_cast<float>(_reader->Px(IHEP)),
			 static_cast<float>(_reader->Py(IHEP)),
			 static_cast<float>(_reader->Pz(IHEP))};
	  
	  //  Mass
	  float mass = _reader->M(IHEP);
	  
	  // Energy
	  float energy = _reader->E(IHEP);
	  
	  //  Generator status
	  int status = _reader->status(IHEP);
	  
	  //  Vertex
	  double v0[3] = {_reader->X(IHEP),_reader->Y(IHEP),_reader->Z(IHEP)};
	  
	  //  Creation time (note the units)
	  double time = _reader->T(IHEP)/c_light;
	  
	  float spin[3] = { 0, 0, 0 };
	  int colorFlow[2] = { 0, 0 };

	  // add spin and color flow information if available 
	  if( _reader->isStdHepEv4() ){

	    spin[0] = _reader->spinX( IHEP ); 
	    spin[1] = _reader->spinY( IHEP ); 
	    spin[2] = _reader->spinZ( IHEP );
	  
	    colorFlow[0] = _reader->colorflow( IHEP , 0 );
	    colorFlow[1] = _reader->colorflow( IHEP , 1 );

	  }
	  
	  int fd = _reader->daughter1(IHEP)%10000 - 1;
	  int ld = _reader->daughter2(IHEP)%10000 - 1;
	  if (fd < 0) { fd = 0; };
	  if (ld < 0) { ld = 0; };
	  
	  sprintf(buff, "[%8.8d]%5d|%10d|% 1.2f|% 1.2e,% 1.2e, % 1.2e|% 1.2e|% 1.2e| %1d |% 1.2e,% 1.2e,% 1.2e|% 1.2e,% 1.2e,% 1.2e|  (%d, %d)   | [ %d, %d ]", 
		  IHEP, IHEP, pdgid, charge, p0[0], p0[1], p0[2], energy, mass, status, v0[0], v0[1], v0[2], spin[0], spin[1], spin[2], colorFlow[0], colorFlow[1], fd, ld);
	  os << buff << std::endl;
	  
	}// End loop over particles
    }
    else if ( blockid == 203 ) {

      os << " Event #" << count - 1 << std::endl;

      int NHEP = _reader->eup_nhep();
      int pid = _reader->eup_pid();
      float weight = _reader->eup_weight();
      float scale = _reader->eup_scale();
      float aqed = _reader->eup_aQED();
      float aqcd = _reader->eup_aQCD();

      os << "Process ID = " <<  pid << "   Weight = " << weight << 
	"     Scale = " << scale << "     aQED = " << aqed << "     aQCD = " << 
	aqcd << std::endl;

      if (count == 1) { 
	os << "[   id   ]index|      PDG |q(el)|    px,      py,        pz    | energy  |   mass  |stat|   spin  |lifetime | colorflow | [mothers]" << 
	  std::endl; }

      for( int IHEP=0; IHEP<NHEP; IHEP++ )
	{
	  //  Fill particle from StdHep info
	  
	  char buff[215];

	  //  PDGID
	  int pdgid = _reader->eup_pid(IHEP);
	   
	  //  charge
	  float charge = threeCharge( pdgid ) / 3.  ;

	  //  Momentum vector
	  float p0[3] = {static_cast<float>(_reader->eup_Px(IHEP)),
			 static_cast<float>(_reader->eup_Py(IHEP)),
			 static_cast<float>(_reader->eup_Pz(IHEP))};
	  
	  //  Mass
	  float mass = _reader->eup_M(IHEP);
	   
	  // Energy
	  float energy = _reader->eup_E(IHEP);

	  // add spin and color flow information  
	  float spin = _reader->eup_spin( IHEP );
	  float vtimeup = _reader->eup_vtime( IHEP );
	  int colorFlow[2] = {
	    static_cast<int>(_reader->eup_colflow( IHEP , 0 )), 
	    static_cast<int>(_reader->eup_colflow( IHEP , 1 ))};
	  // 
	  int fm = _reader->eup_mothup1(IHEP) - 1;
	  int lm = _reader->eup_mothup2(IHEP) - 1;
	  if (fm < 0) { fm = 0; };
	  if (lm < 0) { lm = 0; };

	  //  Generator status
	  int status = _reader->eup_status(IHEP);
	  if (status > 0) {	  
	    sprintf(buff, "[%8.8d]%5d|%10d|% 1.2f|% 1.2e,% 1.2e, % 1.2e|% 1.2e|% 1.2e| +%1u |% 1.2e|% 1.2e|  (%d, %d)   | [ %d, %d ]",
		    IHEP, IHEP, pdgid, charge, p0[0], p0[1], p0[2], energy, mass, status, 
		    spin, vtimeup, colorFlow[0], colorFlow[1], fm, lm); 
	    os << buff << std::endl;      
	  } else {
	    status = - status;
	    sprintf(buff, "[%8.8d]%5d|%10d|% 1.2f|% 1.2e,% 1.2e, % 1.2e|% 1.2e|% 1.2e| -%1u |% 1.2e|% 1.2e|  (%d, %d)   | [ %d, %d ]",
		    IHEP, IHEP, pdgid, charge, p0[0], p0[1], p0[2], energy, mass, status, 
		    spin, vtimeup, colorFlow[0], colorFlow[1], fm, lm); 
	    os << buff << std::endl;  
	  }    
	}
    }
    return 0;
  }


  int WOStdHepRdr::threeCharge( int pdgID ) const {
    //
    // code copied from HepPDT package, author L.Garren
    // modified to take pdg
    
    ///  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj
    ///  The location enum provides a convenient index into the PID.
    enum location { nj=1, nq3, nq2, nq1, nl, nr, n, n8, n9, n10 };

    int charge=0;
    int ida, sid;
    unsigned short q1, q2, q3;
    static int ch100[100] = { -1, 2,-1, 2,-1, 2,-1, 2, 0, 0,
			      -3, 0,-3, 0,-3, 0,-3, 0, 0, 0,
			      0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
			      0, 0, 0, 3, 0, 0, 3, 0, 0, 0,
			      0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
			      0, 6, 3, 6, 0, 0, 0, 0, 0, 0,
			      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			      0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    ida = (pdgID < 0) ? -pdgID : pdgID ;
    
    //     q1 = digit(nq1);
    //     q2 = digit(nq2);
    //     q3 = digit(nq3);

    q1 =  ( ida / ( (int) std::pow( 10.0, (nq1 -1) ) )  ) % 10 ;
    q2 =  ( ida / ( (int) std::pow( 10.0, (nq2 -1) ) )  ) % 10 ;
    q3 =  ( ida / ( (int) std::pow( 10.0, (nq3 -1) ) )  ) % 10 ;
    
//     sid = fundamentalID();
    //---- ParticleID::fundamentalID -------
    short dig_n9 =  ( ida / ( (int) std::pow( 10.0, (n9 -1) ) )  ) % 10 ;
    short dig_n10 =  ( ida / ( (int) std::pow( 10.0, (n10 -1) ) )  ) % 10 ;
    
    if( ( dig_n10 == 1 ) && ( dig_n9 == 0 ) ) {
      
      sid = 0 ;
    } 
    else if( q2 == 0 && q1 == 0) {
      
      sid = ida % 10000;
    } 
    else if( ida <= 102 ) {
      
      sid = ida ; 
    } 
    else {

      sid = 0;
    }
    //----------------

    int extraBits = ida / 10000000 ;
    // everything beyond the 7th digit (e.g. outside the numbering scheme)

    short dig_nj =  ( ida / ( (int) std::pow( 10.0, (nj -1) ) )  ) % 10 ;

    if( ida == 0 || extraBits > 0 ) {      // ion or illegal
      return 0;
    } else if( sid > 0 && sid <= 100 ) {	// use table
      charge = ch100[sid-1];
      if(ida==1000017 || ida==1000018) { charge = 0; }
      if(ida==1000034 || ida==1000052) { charge = 0; }
      if(ida==1000053 || ida==1000054) { charge = 0; }
      if(ida==5100061 || ida==5100062) { charge = 6; }
    } else if( dig_nj == 0 ) { 		// KL, Ks, or undefined
      return 0;
    } else if( q1 == 0 ) {			// mesons
      if( q2 == 3 || q2 == 5 ) {
	charge = ch100[q3-1] - ch100[q2-1];
      } else {
	charge = ch100[q2-1] - ch100[q3-1];
      }
    } else if( q3 == 0 ) {			// diquarks
      charge = ch100[q2-1] + ch100[q1-1];
    } else { 					// baryons
      charge = ch100[q3-1] + ch100[q2-1] + ch100[q1-1];
    }
    if( charge == 0 ) {
      return 0;
    } else if( pdgID < 0 ) {
      charge = -charge; 
    }
    return charge;
  }

} // namespace UTIL

