//////////////////////////////////////////////////////////////////////////
// HepMC reader for test case
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace HepMC;

void read_event_file()
{
    std::cout << "Reading HepMC event file:" << std::endl;
    {
        // Open input stream
	std::ifstream istr( "hepmc_6_p.hepmc" );
	if( !istr ) {
	  std::cerr << "Cannot open HepMC event file" << std::endl;
	  exit(-1);
	}
	HepMC::IO_GenEvent ascii_in(istr);

	// Now read the file
	int icount=0;
	HepMC::GenEvent* evt = ascii_in.read_next_event();
	while ( evt ) {
	  icount++;
	  // evt->print();  // [would print all]

	  // Extract some info from this event
	  std::cout << "Read event #" << evt->event_number()
		    << std::endl;
	  std::cout << "  scale    = " << evt->event_scale()
		    << std::endl;
	  std::cout << "  alphaQCD = " << evt->alphaQCD()
		    << std::endl;
	  std::cout << "  alphaQED = " << evt->alphaQED()
		    << std::endl;

	  // Done event analysis
	  delete evt;
	  ascii_in >> evt;
	}

	// Summary printout
	std::cout << icount << " events found. Finished." << std::endl;
    } // ascii_out and istr destructors are called here
}


int main() {
  read_event_file();
  return 0;
}
