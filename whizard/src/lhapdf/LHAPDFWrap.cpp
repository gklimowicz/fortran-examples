#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/Info.h"
#include "LHAPDF/Exceptions.h"

using namespace LHAPDF;
using namespace std;

extern "C" {

  LHAPDF::PDF* lhapdf_init_pdf (char* setname, const int imem) {
    LHAPDF::PDF* pdf = LHAPDF::mkPDF(setname, imem);
    return pdf;
  }

  void lhapdf_pdf_delete (LHAPDF::PDF* pdf) {
    delete pdf;
  }

  double lhapdf_pdf_getxmin (LHAPDF::PDF* pdf) {
    return pdf->xMin();
  }

  double lhapdf_pdf_getxmax (LHAPDF::PDF* pdf) {
    return pdf->xMax();
  }

  double lhapdf_pdf_getq2min (LHAPDF::PDF* pdf) {
    return pdf->q2Min();
  }

  double lhapdf_pdf_getq2max (LHAPDF::PDF* pdf) {
    return pdf->q2Max();
  }

  bool lhapdf_has_photon (const LHAPDF::PDF* pdf) {
    return pdf->hasFlavor(22);
  }

  /// Get xf(x) values for common partons from PDF pdf
  //  Evaluate for the 13 LHAPDF5 standard partons
  void lhapdf_evolvepdfm (const LHAPDF::PDF* pdf, const double x, const double q, double* fxq) {
    for (size_t i = 0; i < 13; ++i) {
      fxq[i] = pdf->xfxQ(i-6, x, q);
    }
  }

  /// Get xfx values from current PDF, including an extra photon flavor
  void lhapdf_evolvepdfphotonm (const LHAPDF::PDF* pdf, const double x, const double q, double* fxq, double &photonfxq) {
    lhapdf_evolvepdfm (pdf, x, q, fxq);
    photonfxq = pdf->xfxQ(22, x, q);
  }

  void lhapdf_evolvepdfpm (const LHAPDF::PDF* pdf, const double x, const double q, const double s, const int scheme, double fxq) {
    throw LHAPDF::NotImplementedError("Photon structure function are not yet supported");
  }

  double lhapdf_getqmass (const LHAPDF::PDF* pdf, const int nf) {
    if (nf*nf == 1) return pdf->info().get_entry_as<double>("MDown");
    else if (nf*nf == 4) return pdf->info().get_entry_as<double>("MUp");
    else if (nf*nf == 9) return pdf->info().get_entry_as<double>("MStrange");
    else if (nf*nf == 16) return pdf->info().get_entry_as<double>("MCharm");
    else if (nf*nf == 25) return pdf->info().get_entry_as<double>("MBottom");
    else if (nf*nf == 36) return pdf->info().get_entry_as<double>("MTop");
    else throw LHAPDF::UserError("Trying to get quark mass for invalid quark ID #" + LHAPDF::to_str(nf));
  }

  int lhapdf_numpdfm (const LHAPDF::PDF* pdf, int numpdf) {
    numpdf = pdf->info().get_entry_as<int>("NumMembers");
    if (numpdf > 1)  numpdf-=1;
    return numpdf;
  }

  double lhapdf_alphaspdf (const LHAPDF::PDF* pdf, const double q) {
    return pdf->alphasQ(q);
  }

  int lhapdf_getorder (const LHAPDF::PDF* pdf, int order) {
    order = pdf->info().get_entry_as<int>("OrderQCD");
    return order;
  }

}
