#ifndef LHAWhizard_Events_H
#define LHAWhizard_Events_H

extern "C" {
  typedef struct {
    int idPart, statusPart;
    int motherPart[2];
    int colorPart[2];
    double pPart[4];
    double mPart, tauPart, spinPart;
  } lha_particle_t;
}

#endif
