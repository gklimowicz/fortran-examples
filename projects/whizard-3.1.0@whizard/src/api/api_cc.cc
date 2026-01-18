#include <string>
#include <vector>
#include "whizard.h"

#ifdef WHIZARD_WITH_HEPMC3
#include "HepMC3/GenEvent.h"
#endif
#ifdef WHIZARD_WITH_HEPMC2
#include "HepMC/GenEvent.h"
#endif

using namespace std;

Whizard::Whizard() {
  whizard_create( &wh );
};

void Whizard::option( const string key, const string value ) {
  whizard_option( &wh, &key[0], &value[0] );
};

void Whizard::init() {
  whizard_init( &wh );
};

Whizard::~Whizard() {
  whizard_final( &wh );
};

void Whizard::set_double( const string var, const double value ) {
  whizard_set_double ( &wh, &var[0], value );
};

void Whizard::set_int( const string var, const int value ) {
  whizard_set_int ( &wh, &var[0], value );
};

void Whizard::set_bool( const string var, const int value ) {
  whizard_set_bool ( &wh, &var[0], value );
};

void Whizard::set_string( const string var, const string value ) {
  whizard_set_char ( &wh, &var[0], &value[0] );
};

int Whizard::get_double( const string var, double* value ) {
  return whizard_get_double ( &wh, &var[0], value );
};

int Whizard::get_int( const string var, int* value ) {
  return whizard_get_int ( &wh, &var[0], value );
};

int Whizard::get_bool( const string var, int* value ) {
  return whizard_get_bool ( &wh, &var[0], value );
};

int Whizard::get_string( const string var, string** value ) {
  int strlen = whizard_get_char_len( &wh, &var[0] );
  char* tmp = new char ( strlen );
  int known = whizard_get_char ( &wh, &var[0], tmp, strlen );
  *value = new string (tmp);
  delete tmp;
  return known;
};

string* Whizard::flv_string( const int f ) {
  int strlen = whizard_flv_string_len( &wh, f );
  char* tmp = new char ( strlen );
  whizard_flv_string ( &wh, f, tmp, strlen+1 );
  string *value = new string (tmp);
  delete tmp;
  return value;
};

string* Whizard::flv_array_string( const vector<int> fa ) {
  int strlen = whizard_flv_array_string_len( &wh, &fa[0], fa.size() );
  char* tmp = new char ( strlen );
  whizard_flv_array_string ( &wh, &fa[0], fa.size(), tmp, strlen+1 );
  string *value = new string (tmp);
  delete tmp;
  return value;
};

void Whizard::command( const string cmd ) {
  whizard_command( &wh, &cmd[0] );
};

int Whizard::get_integration_result( const string proc_id, double* integral, double* error ) {
  return whizard_get_integration_result( &wh, &proc_id[0], integral, error );
};

WhizardSample* Whizard::new_sample( const string name ) {
  return new WhizardSample( wh, name );
};


WhizardSample::WhizardSample( void* wh, const string name ) {
  whizard_new_sample( &wh, &name[0], &sample );
};

void WhizardSample::open( int* it_begin, int* it_end ) {
  whizard_sample_open( &sample, it_begin, it_end );
};

void WhizardSample::next_event() {
  whizard_sample_next_event( &sample );
};

void WhizardSample::close() {
  whizard_sample_close( &sample );
};

int WhizardSample::get_event_index() {
  int idx;
  whizard_sample_get_event_index( &sample, &idx );
  return idx;
};

int WhizardSample::get_process_index() {
  int i_proc;
  whizard_sample_get_process_index( &sample, &i_proc );
  return i_proc;
};

string* WhizardSample::get_process_id() {
  int strlen = whizard_sample_get_process_id_len( &sample );
  char* tmp = new char ( strlen );
  whizard_sample_get_process_id ( &sample, tmp, strlen+1 );
  string *proc_id = new string (tmp);
  delete tmp;
  return proc_id;
};

double WhizardSample::get_fac_scale() {
  double f_scale;
  whizard_sample_get_fac_scale( &sample, &f_scale );
  return f_scale;
};

double WhizardSample::get_alpha_s() {
  double alpha_s;
  whizard_sample_get_alpha_s( &sample, &alpha_s );
  return alpha_s;
};

double WhizardSample::get_weight() {
  double weight;
  whizard_sample_get_weight( &sample, &weight );
  return weight;
};

double WhizardSample::get_sqme() {
  double sqme;
  whizard_sample_get_sqme( &sample, &sqme );
  return sqme;
};

#ifdef WHIZARD_WITH_HEPMC2
void WhizardSample::next_event( HepMC::GenEvent** evt ) {
  *evt = whizard_sample_next_event_hepmc( &sample );
};
#endif
#ifdef WHIZARD_WITH_HEPMC3
void WhizardSample::next_event( HepMC3::GenEvent** evt ) {
  *evt = whizard_sample_next_event_hepmc( &sample );
};
#endif

#ifdef WHIZARD_WITH_LCIO
void WhizardSample::next_event( lcio::LCEvent** evt ) {
  *evt = whizard_sample_next_event_lcio( &sample );
};
#endif

