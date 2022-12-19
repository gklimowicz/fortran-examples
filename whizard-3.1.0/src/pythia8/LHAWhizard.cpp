// Function definitions (not found in the header) for the LHAupWhizard class.

#include "LHAWhizard.h"

namespace Pythia8 {
  bool LHAupWhizard::setInit (const int beam_pdg[2], const double beam_energy[2], const int n_processes, const bool unweighted, bool negative_weights) {
    /* WHIZARD provides the full hadronic matrix element. Hence, we reset the PDF convolution from Pythia8. */
    // We support only unweighted events with unit weight (correspond to strategy 3) and weighted events with weight in pb.
    // As we are only supporting strategy 3 or 4, Pythia8 handles all processes commonly, no need to provide a meaningfull n_processes.
    // We take care about the unit convention as WHIZARD uses fb.
    int strategy = (unweighted) ? 3 : 4;
    if (negative_weights) strategy *= -1;
    setBeamA(beam_pdg[0], beam_energy[0], -1, -1);
    setBeamB(beam_pdg[1], beam_energy[1], -1, -1);
    setStrategy(strategy);
    return true;
  }

  /* Fill the characteristics for a particular process.
     WHIZARD only supports strategy 3 and 4.
  */
  bool LHAupWhizard::setProcessParameters(const int idProcess, const double cross_section, const double error, const double max_weight) {
    // Althought we accept a process id, Pythia internally discards the number.
    // The only accpeted trategies 3 or 4 (unweighted or weighted events) don't care about the process specifics and handle all events commonly.
    int max_weight_ = 0;
    if (abs(strategy ()) == 3) {
      max_weight_ = 1;
    } else if (abs(strategy () == 4)) {
      max_weight_ = max_weight;
    }
    addProcess(idProcess, cross_section, error, max_weight_);
    return true;
  }

  void LHAupWhizard::setEventProcess(const int idProcess, const double scale, const double alpha_qcd, const double alpha_qed, const double weight) {
    setProcess (idProcess, weight, scale, alpha_qed, alpha_qcd);
  }

  bool LHAupWhizard::setEvent(const int idProcess, const int n_particles, const lha_particle_t particle_set[]) {
    for (int ip = 0; ip < n_particles; ++ip) {
      addParticle(particle_set[ip].idPart, particle_set[ip].statusPart,
                  particle_set[ip].motherPart[0], particle_set[ip].motherPart[1],
                  particle_set[ip].colorPart[0], particle_set[ip].colorPart[1],
                  particle_set[ip].pPart[1], particle_set[ip].pPart[2], particle_set[ip].pPart[3], particle_set[ip].pPart[0],
                  particle_set[ip].mPart, particle_set[ip].tauPart, particle_set[ip].spinPart, -1.);
    }

    int id1 = id(1); int id2 = id(2);
    double x1 = (eBeamA() > 0.) ? e(1) / eBeamA() : 0.;
    double x2 = (eBeamB() > 0.) ? e(2) / eBeamB() : 0.;
    setIdX (id1, id2, x1, x2);
    // TODO sbrass What happens if we have pp beams?
    // Set id and x values even when not supplied.
    setPdf (id1, id2, x1, x2, 0., 0., 0., 0.);
    return true;
  }
}
