#include "Pythia8Wrap.h"

using namespace Pythia8;

// Define external C-function.

extern "C" {

  bool pythia8_available () {
    return true;
  }

  Pythia* new_pythia8 (bool printBanner = true) {
    Pythia* pythia = new Pythia ("", printBanner);
    return pythia;
  }

  void pythia8_delete (Pythia* pythia) {
    delete pythia;
  }

#if PYTHIA_VERSION_INTEGER > 8300
  bool pythia8_set_lhaup_ptr (Pythia* pythia, LHAupWhizard* whizard_lha) {
    return pythia->setLHAupPtr ((LHAupPtr)whizard_lha);
  }
#else
  bool pythia8_set_lhaup_ptr (Pythia* pythia, LHAupWhizard* whizard_lha) {
    return pythia->setLHAupPtr (whizard_lha);
  }
#endif

  bool pythia8_set_rndm_engine_ptr (Pythia* pythia, void* rndm) {
    WhizardRndm* whizard_rndm = new WhizardRndm (rndm);
    return pythia->setRndmEnginePtr (whizard_rndm);
  }

  bool pythia8_read_string (Pythia* pythia, const char* str) {
    return pythia->readString (string (str));
  }

  bool pythia8_read_file (Pythia* pythia, const char* filename, const int subrun) {
    return pythia->readFile (string (filename), subrun);
  }

  bool pythia8_init (Pythia* pythia) {
    return pythia->init ();
  }

  bool pythia8_next (Pythia* pythia) {
    return pythia->next ();
  }

  int pythia8_get_event_size (Pythia* pythia) {
    return pythia->event.size ();
  }

  lha_particle_t pythia8_get_single_event (Pythia* pythia, const int i) {
    Particle event = pythia->event[i];
    lha_particle_t lha_event;
    lha_event.idPart = event.id();
    lha_event.statusPart = event.status();
    lha_event.motherPart[0] = event.mother1();
    lha_event.motherPart[1] = event.mother2();
    lha_event.colorPart[0] = event.col();
    lha_event.colorPart[1] = event.acol();
    lha_event.pPart[0] = event.e();
    lha_event.pPart[1] = event.px();
    lha_event.pPart[2] = event.py();
    lha_event.pPart[3] = event.pz();
    lha_event.mPart = event.m();
    lha_event.tauPart = event.tau();
    lha_event.spinPart = event.pol();
    // TODO sbrass Implement scale of parton
    return lha_event;
  }

  int pythia8_get_particle_status (Pythia* pythia, const int i) {
    return pythia->event[i].status();
  }

  int pythia8_get_particle_id (Pythia* pythia, const int i) {
    return pythia->event[i].id();
  }

  void pythia8_get_particle_momentum (Pythia* pythia, const int i, double* p) {
    // TODO sbrass implement a check on size of p
    Particle event = pythia->event[i];
    p[0] = event.e();
    p[1] = event.px();
    p[2] = event.py();
    p[3] = event.pz();
  }

  int pythia8_get_n_mothers (Pythia* pythia, const int i) {
    vector<int> motherList = pythia->event[i].motherList ();
    return motherList.size ();
  }

  int pythia8_get_n_daughters (Pythia* pythia, const int i) {
    vector<int> daughterList = pythia->event[i].daughterList ();
    return daughterList.size ();
  }

  void pythia8_get_mother_array (Pythia* pythia, const int i, const int n_mothers, int mother[]) {
    vector<int> motherList = pythia->event[i].motherList ();
    if (std::size_t (n_mothers) != motherList.size ()) {
      std::cerr << "[pythia8_get_mother_array] mismatch in array size." << endl;
      exit(2);
    }
    std::copy (motherList.begin (), motherList.end (), mother);
  }

  void pythia8_get_daughter_array (Pythia* pythia, const int i, const int n_daughters, int* daughter) {
    vector<int> daughterList = pythia->event[i].daughterList ();
    if (std::size_t (n_daughters) != daughterList.size ()) {
      std::cerr << "[pythia8_get_mother_array] mismatch in array size." << endl;
      exit(2);
    }
    std::copy (daughterList.begin (), daughterList.end (), daughter);
  }

  int pythia8_get_status_hepmc (Pythia* pythia, const int i) {
    return pythia->event[i].statusHepMC();
  }

  void pythia8_get_decay_vertex (Pythia* pythia, const int i, double* time, double space_vertex[3]) {
    *time = pythia->event[i].tDec();
    space_vertex[0] = pythia->event[i].xDec();
    space_vertex[1] = pythia->event[i].yDec();
    space_vertex[2] = pythia->event[i].zDec();
  }

  void pythia8_get_production_vertex (Pythia* pythia, const int i, double* time, double space_vertex[3]) {
    *time = pythia->event[i].tProd();
    space_vertex[0] = pythia->event[i].xProd();
    space_vertex[1] = pythia->event[i].yProd();
    space_vertex[2] = pythia->event[i].zProd();
  }

  void pythia8_get_event_info (Pythia* pythia, double* alpha_s, double* alpha_em, double* scale) {
    *alpha_s = pythia->info.alphaS();
    *alpha_em = pythia->info.alphaEM();
    *scale = pythia->info.scalup();
    // See: http://home.thep.lu.se/~torbjorn/pythia82html/EventInformation.html
  }

  // Debug Output
  void pythia8_list_lha_event (Pythia* pythia) {
    pythia->LHAeventList ();
  }

  void pythia8_list_event (Pythia* pythia) {
    pythia->event.list ();
  }
 
}


