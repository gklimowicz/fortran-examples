/* Pythia8Wrap.h

   Copyright Notice of WHIZARD

   Author: Simon Braß, Octobre 2018

 */

#ifndef Pythia8_Wrap_H
#define Pythia8_Wrap_H
#include "Pythia8/Pythia.h"
#include "LHAWhizard_events.h"
#include "LHAWhizard.h"

extern "C" {
  double whizard_rndm_generate (void* rndm);
}

namespace Pythia8 {

  // --------------------------------------------------
  // Random Engine
  //
  // Import pointer.

  class WhizardRndm : public RndmEngine {

  public:

    WhizardRndm (void* rndmIn) { rndm = rndmIn; }

    ~WhizardRndm () {};

    double flat () { return whizard_rndm_generate (rndm); }

  private:
    // (void) pointer to a Fortran object
    void* rndm;
  };
}
#endif
