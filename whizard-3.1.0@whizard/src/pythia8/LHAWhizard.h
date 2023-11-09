/* LHAWHIZARD.h

   Copyright Notice of WHIZARD

   Author: Simon Braß, May 2018
*/

#ifndef Pythia8_LHAWhizard_H
#define Pythia8_LHAWhizard_H

#include "Pythia8/Pythia.h"
#include "LHAWhizard_events.h"
#include <iostream>

namespace Pythia8 {
  /* A derived class from LHAupWhizard which retrieves an event from WHIZARD.

     We provide an interface to WHIZARD to Pythia.
  */
  class LHAupWhizard : public LHAup {

  public:

    ~LHAupWhizard () {}

    /* Define the basic setInit and overload it. */
    bool setInit () {return true;}
    bool setInit (const int beam_pdg[2], const double beam_energy[2], const int n_processes, const bool unweighted, bool negative_weights);
    bool setProcessParameters (const int idProcess, const double cross_section, const double error, const double max_weight);
    bool setEvent (int idProcess) {return true;}
    void setEventProcess (const int idProcess, const double scale, const double qlpha_qcd, const double alpha_qed, const double weight);
    bool setEvent (const int idProcess, const int n_particles, const lha_particle_t particle_set[]);

    /* We use events provided directly by a Fortran-to-C++ interfeace from WHIZARD and do not need to read events from file. */
    bool fileFound () {return true;}
    bool useExternal () {return true;}

    /* We do not need the member functions introduced for LHEF. Just ignoring them. */
    void newEventFile (const char*) {}
    bool skipEvent(int nSkip) {return true;}
    bool openLHEF (string filenameIn) {return false;}
    bool closeLHEF (bool updateInit = false) {return false;}
  };
}
#endif
