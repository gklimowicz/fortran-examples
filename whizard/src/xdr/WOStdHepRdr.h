#ifndef WOUTIL_LCStdHepRdr_H
#define WOUTIL_LCStdHepRdr_H 1

#include "WOStdHep.hh"
#include <iostream>

namespace WOUTIL{
  
  /**Basic utility for reading a binary stdhep file and filling
   * a LCCollectionVec with MCParticles containing the stdhep
   * file information.
   * 
   * @author cassell
   * @version $Id: WOStdHepRdr.h,v 1.4 2007-11-12 16:39:04 gaede Exp $
   */
  class WOStdHepRdr{
    
  public:

	/** Open the stdhep input file in the constructer
	 */
    WOStdHepRdr(const char* evfile) ;

	/** noop
	 */
	~WOStdHepRdr() ;


	/** Get number of events in the stdhep file.
     *  This number is read from the file header (no guarantee that it is correct)
    */
    long getNumberOfEvents() const {
      return _reader->numEvents() ;
    }
    
    /** Get total number of expected events in the whole set of stdhep files
     *  from which this stdhep file belongs to.
     *  This number is read from the file header (no guarantee that it is correct)
     */
    long getNumberOfTotalEventsExpected() const {
      return _reader->numEventsExpected() ;
    }

    /** Read an event and return an LCCollectionVec of MCParticles.
     * @deprecated please use updateEvent()
     */
    int readEvent(std::ostream& os = std::cout, int count = 0 ) ;

    /** Print the file header to the given ostream.
     */
    void printHeader(std::ostream& os = std::cout ) ; 


    /** Return the charge of the particle times 3  - code copied from HepPDT package.
     */
    int threeCharge( int pdgID ) const ;

  private:
    
	WOStdHep* _reader;
    

  }; // class

} // namespace WOUTIL

#endif /* ifndef WOUTIL_WOStdHepRdr_H */
