/* API for executing WHIZARD unit tests                             */

/* **************************************************************** */
/* Plain C part, interfaces with Fortran bind(C) procedures         */

#ifdef __cplusplus
extern "C" {
#endif

void whizard_ut_setup( const char* ut_name, int* u_log, void* results );
void whizard_ut_wrapup( const int u_log, void* results );
int whizard_ut_get_n_pass( void* results );
int whizard_ut_get_n_fail( void* results );
int whizard_ut_get_n_total( void* results );
void whizard_ut_start( const int u_log, const char* name );
void whizard_ut_end( const int u_log, const char* name, const char* description, void* results );

#ifdef __cplusplus
}
#endif

/* **************************************************************** */
/* C++ part: wrapper class for Fortran test utilities               */

#ifdef __cplusplus

#include <string>
#include <cstdio>

using namespace std;

class WhizardCCTest {

 private:
  string name;
  int u_log;
  void* results;

 public:
  WhizardCCTest( const char* test_collection_name );
  int get_n_fail();
  ~WhizardCCTest();
  void run_test( void (*f)(FILE*), const char* test_name, const char* test_description );

};

#endif
