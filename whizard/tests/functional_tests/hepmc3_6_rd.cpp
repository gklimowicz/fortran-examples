//////////////////////////////////////////////////////////////////////////
// HepMC reader for test case
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Print.h"

using namespace HepMC3;

void read_event_file()
{
    std::cout << "Reading HepMC event file:" << std::endl;
     {
       // Open input stream
       std::ifstream stream( "hepmc_6_p.hepmc" );
       if( !stream ) {
	 std::cerr << "Cannot open HepMC event file" << std::endl;
	 exit(-1);
       }
       HepMC3::ReaderAscii ascii_in(stream);

       // Now read the file
       int icount=0;
       GenEvent evt = GenEvent(Units::GEV,Units::MM);
       bool ok;
       int num = 1;
       while ( num > 0) {
	 ok = ascii_in.read_event(evt);
	 num = evt.event_number();
	 if (num > 0) {
	   icount++;
	   // Print::listing (evt); // [would print all]
	   // Extract some info from this event
	   std::cout << "Read event #" << num
		     << std::endl;
	   std::cout << "  scale    = " <<
	     evt.attribute<DoubleAttribute>("event_scale")->value()
		     << std::endl;
	   std::cout << "  alphaQCD = " <<
	     evt.attribute<DoubleAttribute>("alphaQCD")->value()
		     << std::endl;
	   std::cout << "  alphaQED = " <<
	     evt.attribute<DoubleAttribute>("alphaQED")->value()
		     << std::endl;
	 }
	 // Done event analysis
	}
	// Summary printout
	std::cout << icount << " events found. Finished." << std::endl;
     }
}

int main() {
  read_event_file();
  return 0;
}
