#include "Pythia8/Pythia.h"
#include "LHAWhizard.h"
#include "LHAWhizard_events.h"

using namespace Pythia8;

extern "C" {

  // Tell the caller that this is the true LHAup interface and the true Pythia8 library
  bool lhaup_whizard_available() {
    return true;
  }

  LHAupWhizard* new_whizard_lha () {
    LHAupWhizard* whizard_lha = new LHAupWhizard ();
    return whizard_lha;
  }

  void lhaup_whizard_delete (LHAupWhizard* whizard_lha) {
    delete whizard_lha;
  }

  bool lhaup_whizard_set_init (LHAupWhizard* whizard_lha, const int beam_pdg[2], const double beam_energy[2], const int n_processes, const bool unweighted, const bool negative_weights ) {
    return whizard_lha->setInit (beam_pdg, beam_energy, n_processes, unweighted, negative_weights);
  }

  bool lhaup_whizard_set_process_parameters (LHAupWhizard* whizard_lha, const int process_id, const double cross_section, const double error, const double max_weight) {
    return whizard_lha->setProcessParameters (process_id, cross_section, error, max_weight);
  }

  void lhaup_whizard_list_init (LHAupWhizard* whizard_lha) {
    whizard_lha->listInit ();
  }

  void lhaup_whizard_list_event (LHAupWhizard* whizard_lha) {
    whizard_lha->listEvent ();
  }

  void lhaup_whizard_set_event_process (LHAupWhizard* whizard_lha, const int idProcess, const double scale, const double alpha_qcd, const double alpha_qed, const double weight) {
    whizard_lha->setEventProcess (idProcess, scale, alpha_qcd, alpha_qed, weight);
  }

  bool lhaup_whizard_set_event (LHAupWhizard* whizard_lha, const int idProcess, const int n_particles, const lha_particle_t particle_set[]) {
    return whizard_lha->setEvent (idProcess, n_particles, particle_set);
  }

}
