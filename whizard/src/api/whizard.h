/* Public API */

/* **************************************************************** */
/* Plain C part, interfaces with Fortran bind(C) procedures */
#ifdef __cplusplus
extern "C" {
#endif

  typedef void* whizard_t;
  typedef void* sample_handle_t;

void whizard_create (void* wh);
void whizard_option (void* wh, const char* key, const char* value);
void whizard_init (void* wh);
void whizard_final (void* wh);

void whizard_set_double (void* wh, const char* var, const double value);
void whizard_set_int (void* wh, const char* var, const int value);
void whizard_set_bool (void* wh, const char* var, const int value);
void whizard_set_char (void* wh, const char* var, const char* value);
int whizard_get_double (void* wh, const char* var, double* value);
int whizard_get_int (void* wh, const char* var, int* value);
int whizard_get_bool (void* wh, const char* var, int* value);
int whizard_get_char (void* wh, const char* var, char* value, const int strlen);
int whizard_get_char_len (void* wh, const char* var);

int whizard_flv_string (void* wh, const int pdg, char* fstr, const int strlen);
int whizard_flv_string_len (void* wh, const int pdg);
int whizard_flv_array_string (void* wh, const int* pdg, const int nf, char* fstr, const int strlen);
int whizard_flv_array_string_len (void* wh, const int* pdg, const int nf);

void whizard_command (void* wh, const char* cmd);
int whizard_get_integration_result (void* wh, const char* proc_id, double* integral, double* error);

void whizard_new_sample (void* wh, const char* name, void* sample);
void whizard_sample_open (void* sample, int* it_begin, int* it_end);
void whizard_sample_next_event (void* sample);
void whizard_sample_close (void* sample);

void whizard_sample_get_event_index (void* sample, int* idx);
void whizard_sample_get_process_index (void* sample, int* i_proc);
void whizard_sample_get_process_id (void* sample, char* proc_id, const int strlen);
int whizard_sample_get_process_id_len (void* sample);
void whizard_sample_get_fac_scale (void* sample, double* f_scale);
void whizard_sample_get_alpha_s (void* sample, double* alpha_s);
void whizard_sample_get_weight (void* sample, double* weight);
void whizard_sample_get_sqme (void* sample, double* sqme);

#ifdef __cplusplus
}
#endif

/* **************************************************************** */
/* C++ part: wrapper classes for Whizard API objects                */

#ifdef __cplusplus

#ifdef WHIZARD_WITH_HEPMC3
#include "HepMC3/GenEvent.h"
extern "C" {
HepMC3::GenEvent* whizard_sample_next_event_hepmc( void* sample );
}
#endif
#ifdef WHIZARD_WITH_HEPMC2
#include "HepMC/GenEvent.h"
extern "C" {
HepMC::GenEvent* whizard_sample_next_event_hepmc( void* sample );
}
#endif

#ifdef WHIZARD_WITH_LCIO
#include "lcio.h"
extern "C" {
lcio::LCEvent* whizard_sample_next_event_lcio( void* sample );
}
#endif


#include <string>
#include <vector>

using namespace std;

class WhizardSample;  // pre-declaration for use inside Whizard

class Whizard {

 private:
  void* wh;

 public:
  Whizard();
  void option( const string key, const string value );
  void init();
  ~Whizard();
  void set_double( const string var, const double value );
  void set_int( const string var, const int value );
  void set_bool( const string var, const int value );
  void set_string( const string var, const string value );
  int get_double( const string var, double* value );
  int get_int( const string var, int* value );
  int get_bool( const string var, int* value );
  int get_string( const string var, string** value );
  string* flv_string( const int f );
  string* flv_array_string( const vector<int> fa );
  void command( const string cmd );
  int get_integration_result( const string proc_id, double* integral, double* error );
  WhizardSample* new_sample( const string name );

};

class WhizardSample {

 private:
  void* sample;

  friend WhizardSample* Whizard::new_sample( const string );
  WhizardSample( void* wh, const string name );

 public:
  void open( int* it_begin, int* it_end );
  void next_event();
  void close();
  int get_event_index();
  int get_process_index();
  string* get_process_id();
  double get_fac_scale();
  double get_alpha_s();
  double get_weight();
  double get_sqme();
#ifdef WHIZARD_WITH_HEPMC2
  void next_event( HepMC::GenEvent** evt );
#endif
#ifdef WHIZARD_WITH_HEPMC3
  void next_event( HepMC3::GenEvent** evt );
#endif
#ifdef WHIZARD_WITH_LCIO
  void next_event( lcio::LCEvent** evt );
#endif

};


#endif
