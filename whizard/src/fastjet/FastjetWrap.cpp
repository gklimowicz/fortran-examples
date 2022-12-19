#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include "fastjet/EECambridgePlugin.hh"

using namespace fastjet;
using namespace std;

extern "C" {

  // Tell the caller that this is the true Fastjet library
  bool fastjet_available() {
    return true;
  }

  void fastjet_print_banner () {
    ClusterSequence::print_banner ();
  }

  PseudoJet* new_pseudojet (const double px_in, const double py_in, const double pz_in, const double E_in) {
    PseudoJet* j = new PseudoJet(px_in, py_in, pz_in, E_in);
    return j;
  }

  void pseudojet_delete (PseudoJet *j) {
    delete j;
  }

  double pseudojet_get_e (const PseudoJet* j) {
    return j->e();
  }

  double pseudojet_get_px (const PseudoJet* j) {
    return j->px();
  }

  double pseudojet_get_py (const PseudoJet* j) {
    return j->py();
  }

  double pseudojet_get_pz (const PseudoJet* j) {
    return j->pz();
  }

  double pseudojet_get_perp (const PseudoJet* j) {
    return j->perp();
  }

  double pseudojet_get_rap (const PseudoJet* j) {
    return j->rap();
  }

  double pseudojet_get_phi (const PseudoJet* j) {
    return j->phi();
  }

  vector<PseudoJet>* pseudojet_get_constituents (const PseudoJet* j) {
    vector<PseudoJet>* cv = new vector<PseudoJet>;
    *cv = j->constituents();
    return cv;
  }

  bool pseudojet_contains (const PseudoJet* j, const PseudoJet* p) {
    return j->contains(*p);
  }

  vector<PseudoJet>* new_pseudojet_vector (const PseudoJet *j[], const int n) {
    vector<PseudoJet>* jv = new vector<PseudoJet>;
    for (int i = 0; i < n; i++) {
      jv->push_back( *j[i] );
    }
    return jv;
  }

  void pseudojet_vector_delete (vector<PseudoJet>* jv) {
    delete jv;
  }

  int pseudojet_vector_get_size (const vector<PseudoJet>* jv) {
    return jv->size();
  }

  PseudoJet* pseudojet_vector_get_jet (const vector<PseudoJet>* jv, const int i) {
    PseudoJet* j = new PseudoJet;
    *j = (*jv)[i];
    return j;
  }

  vector<PseudoJet>* pseudojet_vector_sorted_by_pt (vector<PseudoJet>* jv) {
    vector<PseudoJet>* sjv = new vector<PseudoJet>;
    *sjv = sorted_by_pt(*jv);
    return sjv;
  }

  JetDefinition* new_jet_definition (const JetAlgorithm jet_alg, const double R,
      const double p, const double jet_ycut) {
    JetDefinition *  jet_def;
    if (jet_alg==plugin_algorithm) {
      EECambridgePlugin *eec = new EECambridgePlugin (jet_ycut);
      jet_def = new JetDefinition (eec);
    }
    else if ((jet_alg == genkt_algorithm) || (jet_alg == ee_genkt_algorithm)) {
      jet_def = new JetDefinition (jet_alg, R, p);
    }
    else if (jet_alg == ee_kt_algorithm) {
      jet_def = new JetDefinition (jet_alg);
    }
    else{
      jet_def = new JetDefinition (jet_alg, R);
    }
    return jet_def;
  }

  void jet_definition_delete (JetDefinition* jet_def) {
    delete jet_def;
  }

  int jet_definition_description_strlen (const JetDefinition* jet_def) {
    return jet_def->description().length();
  }

  string* jet_definition_get_description (const JetDefinition* jet_def) {
    string* s = new string (jet_def->description());
    return s;
  }

  ClusterSequence* new_cluster_sequence (const vector<PseudoJet>* particles,
      JetDefinition* jet_def) {
    // run the clustering, extract the jets
    ClusterSequence* cs = new ClusterSequence (*particles, *jet_def);
    return cs;
  }

  void cluster_sequence_delete (ClusterSequence* cs) {
    delete cs;
  }

  vector<PseudoJet>* cluster_sequence_get_inclusive_jets (const ClusterSequence *cs) {
    vector<PseudoJet>* jets = new vector<PseudoJet>;
    *jets = cs->inclusive_jets();
    return jets;
  }

  vector<PseudoJet>* cluster_sequence_get_exclusive_jets (const ClusterSequence *cs,
      const double dcut) {
    vector<PseudoJet>* jets = new vector<PseudoJet>;
    *jets = cs->exclusive_jets(dcut);
    return jets;
  }

  vector<int>* cluster_sequence_get_jet_indices (const ClusterSequence* cs,
      const vector<PseudoJet>* jets) {
    vector<int>* idx = new vector<int>;
    *idx = cs->particle_jet_indices(*jets);
    return idx;
  }

  int int_vector_get (const vector<int>* iv, const int i) {
    return (*iv)[i];
  }

  void int_vector_delete (vector<int>* iv) {
    delete iv;
  }

}
