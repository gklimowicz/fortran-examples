#include "WOStdHepRdr.h"

#include <stdio.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>

using namespace std;
using namespace WOUTIL;

/** Simple test program to demonstrate reading of binary .stdhep generator files.
 */

int main(int argc, char** argv ){
  
  if(argc < 3 || argc > 5) {
    
     std::cout << " usage: stdhep_rd infile.stdhep outfile.slcio maxEvt " <<  std::endl 
	       << "   infile.hep    - input file name "  <<  std::endl 
	       << "   maxEvt        - max number of events to read [-1: all]"  <<  std::endl
	       << "   outfile       - ouput file name [optional]"  <<  std::endl; 
     
    return 1;
  }
  
  ofstream myfile;
  std::string outFile;
  
  std::string inFile  = argv[1] ;
  int maxEvt = std::atoi( argv[2] ) ; 
  if ( argc > 3 ) {
    std::string outFile = argv[3] ;
    myfile.open ( outFile.c_str() ) ;
  }

  std::cout << "==================================================== " << std::endl
	    << " WHIZARD StdHep Reader : " << std::endl  ;

  // Open an instance of the StdHep Reader with the given filename
  WOStdHepRdr rdr( inFile.c_str()  ) ;
  
  std::cout << " opened file : " << inFile << std::endl ;
  
  if ( argc > 3 ) {
    rdr.printHeader( myfile );
  }
  else {
    rdr.printHeader( cout );
  }

  std::stringstream description ; 

  description << " file generated with WHIZARD StdHep Reader from "  << inFile  ;

  int count = 0; 

  while( maxEvt < 0  || count < maxEvt ){
    
    // read the next stdhep event and add an MCParticle collection to the event
    if ( argc > 3 ) {
      rdr.readEvent ( myfile, count );
    } 
    else {
      rdr.readEvent ( std::cout, count );	
    }
    
    ++count ;
    
  } // evt loop

  if ( argc > 3 ) {
    std::cout << "  converted " << count << " events - written to file " << outFile  << std::endl ;
  } else {
    std::cout << "  converted " << count << " events - written to stdout"  << std::endl ;
  }
  
  std::cout << "==================================================== " 
	    << std::endl << std::endl ;
  
  if ( argc > 3 ) {
    myfile.close() ;
  }
  
return 0 ;
 
}

