#include "lcio.h"
#include <stdio.h>
#include <limits>

#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCRunHeader.h" 

#include "EVENT/SimCalorimeterHit.h" 
#include "EVENT/CalorimeterHit.h" 
#include "EVENT/RawCalorimeterHit.h" 

#include "UTIL/CellIDDecoder.h"
#include "UTIL/LCTime.h"

#include <cstdlib>

using namespace std ;
using namespace lcio ;

std::string getSimulatorStatusString(const MCParticle* mcp=0){

  if( mcp == 0  ) {

    std::stringstream str ;
    
    str << "simulator status bits: [sbvtcls] "
	<< " s: created in simulation"
	<< " b: backscatter"
	<< " v: vertex is not endpoint of parent" 
	<< " t: decayed in tracker"
	<< " c: decayed in calorimeter"
	<< " l: has left detector"
	<< " s: stopped"
	<< " o: overlay"
	<< std::endl ;
    
    return str.str() ;
  }
  
  std::string s("[    0   ]") ;
  
  if( mcp->getSimulatorStatus() == 0 ) 
    
    return s ;
  
  if( mcp->isCreatedInSimulation() )
    s[1]='s' ;
  else
    s[1]=' ' ;
  if( mcp->isBackscatter() )
    s[2]='b' ;
  else
    s[2]=' ' ;
  if( mcp->vertexIsNotEndpointOfParent() )
    s[3]='v' ;
  else
    s[3]=' ' ;
  if( mcp->isDecayedInTracker() )
    s[4]='t' ;
  else
    s[4]=' ' ;
  if( mcp->isDecayedInCalorimeter() )
    s[5]='c' ;
  else
    s[5]=' ' ;
  if( mcp->hasLeftDetector() )
    s[6]='l' ;
  else
    s[6]=' ' ;
  if( mcp->isStopped() )
    s[7]='s' ;
  else
    s[7]=' ' ;
  if( mcp->isOverlay() )
    s[8]='o' ;
  else
    s[8]=' ' ;
  
  return s ;
}

// For backwards compatibility of the reference output in the functional tests
// output of the double parameters are output as "float"

void printParameters( const EVENT::LCParameters& params ){
  
  StringVec intKeys ;
  int nIntParameters = params.getIntKeys( intKeys ).size() ;
  for(int i=0; i< nIntParameters ; i++ ){
    IntVec intVec ;
    params.getIntVals(  intKeys[i], intVec ) ;
    int nInt  = intVec.size()  ;   
    cout << " parameter " << intKeys[i] << " [int]: " ; 
    
    if( nInt == 0 ){ 
      cout << " [empty] " << std::endl ;
    }
    for(int j=0; j< nInt ; j++ ){
      cout << intVec[j] << ", " ;
    }
    cout << endl ;
  }
  StringVec floatKeys ;
  int nFloatParameters = params.getFloatKeys( floatKeys ).size() ;
  for(int i=0; i< nFloatParameters ; i++ ){
    FloatVec floatVec ;
    params.getFloatVals(  floatKeys[i], floatVec ) ;
    int nFloat  = floatVec.size()  ;   
    cout << " parameter " << floatKeys[i] << " [float]: " ; 
    if( nFloat == 0 ){ 
      cout << " [empty] " << std::endl ;
    }
    for(int j=0; j< nFloat ; j++ ){
      cout << floatVec[j] << ", " ;
    }
    cout << endl ;
  }
#if LCIO_VERSION_GE (2, 17)
  StringVec doubleKeys ;
  int nDoubleParameters = params.getDoubleKeys( doubleKeys ).size() ;
  for(int i=0; i< nDoubleParameters ; i++ ){
    DoubleVec doubleVec ;
    params.getDoubleVals(  doubleKeys[i], doubleVec ) ;
    int nDouble  = doubleVec.size()  ;
    cout << " parameter " << doubleKeys[i] << " [float]: " ;
    if( nDouble == 0 ){
      cout << " [empty] " << std::endl ;
    }
    for(int j=0; j< nDouble ; j++ ){
      cout << doubleVec[j] << ", " ;
    }
    cout << endl ;
  }
#endif
  StringVec stringKeys ;
  int nStringParameters = params.getStringKeys( stringKeys ).size() ;
  for(int i=0; i< nStringParameters ; i++ ){
    StringVec stringVec ;
    params.getStringVals(  stringKeys[i], stringVec ) ;
    int nString  = stringVec.size()  ;   
    cout << " parameter " << stringKeys[i] << " [string]: " ; 
    if( nString == 0 ){ 
      cout << " [empty] " << std::endl ;
    }
    for(int j=0; j< nString ; j++ ){
      cout << stringVec[j] << ", " ;
    }
    cout << endl ;
  }
  
}


void printMCParticles(const EVENT::LCCollection* col ) {
  
  if( col->getTypeName() != LCIO::MCPARTICLE ){
    
    cout << " collection not of type " << LCIO::MCPARTICLE << endl ;
    return ;
  }
  
  cout << endl 
       << "--------------- " << "print out of "  << LCIO::MCPARTICLE << " collection "
       << "--------------- " << endl ;
  
  cout << endl 
       << "  flag:  0x" << hex  << col->getFlag() << dec << endl ;
  
  printParameters( col->getParameters() ) ;
  
  int nParticles =  col->getNumberOfElements() ;
  
  
  cout << "  " << getSimulatorStatusString() << std::endl ;
  
  // fill map with particle pointers and collection indices
  typedef std::map< MCParticle*, int > PointerToIndexMap ;
  PointerToIndexMap p2i_map ;
  std::vector<MCParticle*> moms ;
  
  /*
    cout << endl;
    MCParticle* part=NULL;
    cout << header(part);
    cout << tail(part);
    
    for( int i=0 ; i< nParticles ; i++ ){
    part = dynamic_cast<MCParticle*>( col->getElementAt( i ));
    cout << lcio_short<EVENT::MCParticle>(part, col); //WARNING! 
    //not exact the same output as the code below:
    //<< lcio_short prints the id's of the parents/daughters particles, 
    //the "old" version prints the position in the collection
    }
  */
  
  for( int k=0; k<nParticles; k++){
    MCParticle* part =  dynamic_cast<MCParticle*>( col->getElementAt( k ) ) ;
    p2i_map[ part ] = k ; 
    
    moms.push_back( part ) ;
  }
  
  std::cout << endl
            <<  "[   id   ]index|      PDG |    px,     py,        pz    | energy  |gen|[simstat ]| vertex x,     y   ,   z     |    mass |  charge |            spin             | colorflow | [parents] - [daughters]"    
            << endl 
            << endl ;
  
  // loop over collection - preserve order
  for(  int index = 0 ; index < nParticles ; index++){
    
    MCParticle* part =  dynamic_cast<MCParticle*>( col->getElementAt( index ) ) ;

#if LCIO_VERSION_GE (2, 13)
    printf("[%8.8d]", part->id() - 1);
#else
    printf("[%8.8d]", part->id() );
#endif
    printf("%5d|"   , index );
    printf("%10d|" , part->getPDG() );
    printf("% 1.2e,% 1.2e,% 1.2e|" , 
	   part->getMomentum()[0] ,
	   part->getMomentum()[1] , 
	   part->getMomentum()[2] );
    printf("% 1.2e|" , part->getEnergy() ) ; 
    
    printf(" %1d |" , part->getGeneratorStatus()  );
    printf("%s|" , getSimulatorStatusString( part ).c_str()  ); 
    printf("% 1.2e,% 1.2e,% 1.2e|" , 
	   part->getVertex()[0] , 
	   part->getVertex()[1] , 
	   part->getVertex()[2] );
    float mass;
    if (abs (part->getMass ()) < part->getEnergy () * std::numeric_limits<float>::epsilon())
      mass = 0;
    else
      mass = part->getMass ();
    printf("% 1.2e|" , mass ) ; 
    printf("% 1.2e|" , part->getCharge() ) ; 
    
    printf("% 1.2e,% 1.2e,% 1.2e|" , 
	   part->getSpin()[0] ,
	   part->getSpin()[1] , 
	   part->getSpin()[2] );
    
    printf("  (%d, %d)   |" , 
	   part->getColorFlow()[0] ,
	   part->getColorFlow()[1] );
    
    cout << " [" ;
    
    for(unsigned int k=0;k<part->getParents().size();k++){
      if(k>0) cout << "," ;
      cout << p2i_map[ part->getParents()[k] ]  ;
    }
    cout << "] - [" ;
    for(unsigned int k=0;k<part->getDaughters().size();k++){
      if(k>0) cout << "," ;
      cout << p2i_map[ part->getDaughters()[k] ]  ;
    }
    cout << "] " << endl ;
  }
  
  cout << endl 
       << "-------------------------------------------------------------------------------- " 
       << endl ;
}

int main(int argc, char** argv ){


  char* FILEN ;
  int runNumber=0 ;
  int evtNumber=0 ;
  int nthEvent=1 ;

  // read file name from command line (only argument) 
  if( argc < 3 ) {

    cout << " usage: dumpevent filename runNum evtNum " << endl ;
    cout << "    or: dumpevent filename n      " << endl ;
    cout << "  where the first dumps the event with the specified run and event number" << endl ;
    cout << "  and the second simply dumps the n-th event in the file" << endl << endl ;
    cout << "  set the environment variable LCIO_READ_COL_NAMES to a space separated list" << endl ;
    cout << "  of collection names that you would like to dump (all are dumped if not set)" << endl ;

    exit(1) ;
  }
  
  FILEN = argv[1] ;

  bool dumpNthEvent( argc == 3 ) ;
 

  if( dumpNthEvent ) {

    nthEvent  = atoi( argv[2] ) ;

    if( nthEvent < 1 ) {

      cout << " usage: dumpevent filename n   -   whith  n > 0 !  " << endl ;
      
      exit(1) ;
    }

  }else{

    runNumber = atoi( argv[2] ) ;
    evtNumber = atoi( argv[3] ) ;
  }
  

//  // set the default encoding for cellid's according to the old Mokka convention
//  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;
//  CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;
//  CellIDDecoder<RawCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

  LCReader* lcReader ;
  if( dumpNthEvent ) 
    lcReader = LCFactory::getInstance()->createLCReader() ;
  else
    lcReader = LCFactory::getInstance()->createLCReader(LCReader::directAccess) ;
  

//  // ------ check if LCIO_READ_COL_NAMES is set -------------
//  
//  char* rColChar = getenv ("LCIO_READ_COL_NAMES");
//  
//  if ( rColChar != 0 ) {
//
//    std::vector< std::string > colSubset ;
//    std::stringstream sts( rColChar ) ;
//    std::string colName;
//
//    while( sts >> colName) {
//    
//      colSubset.push_back( colName ) ;
//    }
//
//    lcReader->setReadCollectionNames(  colSubset ) ;
//  }
//  //-----------------------------------------------------------



  LCEvent* evt(0) ;

  try{
    
     lcReader->open( FILEN ) ;
     
     if( dumpNthEvent ) {
       
     if( nthEvent > 0 )
	   lcReader->skipNEvents(  nthEvent - 1 ) ;

       evt = lcReader->readNextEvent() ; 
       
     }else{
       
       evt = lcReader->readEvent(runNumber,  evtNumber) ; 
     }
   
     if( !evt  ){

       if(dumpNthEvent){

	 cout << " less than " << nthEvent << "  events in  file " << FILEN << endl ;    
	 
       }else{

	 cout << " couldn't find event " << evtNumber << " - run " << runNumber 
	      << " in file " << FILEN << endl ;    
       } 
       
       exit(1) ;
     }

     // the event:
     cout << endl 
	  << "============================================================================" << endl ;
     cout << "        Event  : " << evt->getEventNumber() 
	  << " - run:  "         << evt->getRunNumber()
	  << " - timestamp "     << evt->getTimeStamp()   
	  << " - weight "        << evt->getWeight()   
	  << endl ;
     cout << "============================================================================" << endl ;    

     LCTime evtTime( evt->getTimeStamp() ) ;
     cout << " date:      "      << evtTime.getDateString() << endl ;     
     cout << " detector : "      << evt->getDetectorName() << endl ;
     
     cout << " event parameters: " << endl ; 
     
     printParameters( evt->getParameters() ) ;
     
     
     const std::vector< std::string >* strVec = evt->getCollectionNames() ;
     
     // loop over all collections:
     std::vector< std::string >::const_iterator name ;
     
     for( name = strVec->begin() ; name != strVec->end() ; name++){
       
       LCCollection* col = evt->getCollection( *name ) ;
       
       cout << endl 
	    << " collection name : " << *name 
	    << endl 
	    << " parameters: " << endl ;
       
       //      printParameters( col->getParameters() ) ;
       
       
       
       // call the detailed print functions depending on type name
       if( evt->getCollection( *name )->getTypeName() == LCIO::MCPARTICLE ){
	 
	 printMCParticles( col ) ;
	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::SIMTRACKERHIT ){
//	 
//	 printSimTrackerHits( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::TPCHIT ){
//	 
//	 printTPCHits( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::TRACKERHIT ){
//	 
//	 printTrackerHits( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::TRACKERHITPLANE ){
//	 
//	 printTrackerHitPlane( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::TRACKERHITZCYLINDER ){
//	 
//	 printTrackerHitZCylinder( col ) ;
//
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::SIMCALORIMETERHIT ){
//	 
//	 printSimCalorimeterHits( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::CALORIMETERHIT ){
//	 
//	 printCalorimeterHits( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::RAWCALORIMETERHIT ){
//	 
//	 printRawCalorimeterHits( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::LCFLOATVEC ){
//	 
//	 printLCFloatVecs( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::LCINTVEC ){
//	 
//	 printLCIntVecs( col ) ;                               
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::LCSTRVEC ){
//	 
//	 printLCStrVecs( col ) ;                                 
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::TRACK ){
//	 
//	 printTracks( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::CLUSTER ){
//	 
//	 printClusters( col ) ;
//
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::RECONSTRUCTEDPARTICLE ){
//	 
//	 printReconstructedParticles( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::VERTEX ){
//	 
//	 printVertices( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::LCGENERICOBJECT ){
//	 
//	 printLCGenericObjects( col ) ;
//	 
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::LCRELATION ){
//	 
//	 printRelation( col ) ;
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::TRACKERRAWDATA ){
//	 
//	 printTrackerRawData( col ) ;
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::TRACKERDATA ){
//	 
//	 printTrackerData( col ) ;
//      }
//      else if( evt->getCollection( *name )->getTypeName() == LCIO::TRACKERPULSE ){
//	 
//	 printTrackerPulse( col ) ;
//      }
       
     }
     
  }

  //     LCTOOLS::dumpEventDetailed( evt ) ;
     
     
     lcReader->close() ;
     
   }
   catch( IOException& e) {
     cout << e.what() << endl ;
     exit(1) ;
   }
   return 0 ;
}

