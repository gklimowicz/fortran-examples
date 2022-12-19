#include <iostream>
#include <string>
#include <cstdio>
#include "whizard.h"
#include "whizard_ut.h"

using namespace std;

WhizardCCTest::WhizardCCTest( const char* test_collection_name ) {

  name = test_collection_name;

  cout << "| ============================================================================\n";
  cout << "| Running WHIZARD self-test: " << name << " with C++ driver program\n";
  cout << "| ----------------------------------------------------------------------------\n";

  whizard_ut_setup (test_collection_name, &u_log, &results);

}

int WhizardCCTest::get_n_fail() {

  return whizard_ut_get_n_fail (&results);

}

WhizardCCTest::~WhizardCCTest() {

  whizard_ut_wrapup (u_log, &results);

  cout << "| ----------------------------------------------------------------------------\n";
  cout << "| Finished WHIZARD self-test: " << name << " with C++ driver program\n";
  cout << "| ============================================================================\n";

}

void WhizardCCTest::run_test( void (*f)(FILE*), const char* test_name, const char* test_description ) {

  whizard_ut_start( u_log, test_name );

  string outfile_name;
  outfile_name = test_name;
  outfile_name += ".out";

  FILE *outfile;
  outfile = fopen( &outfile_name[0], "w" );

  f( outfile );

  fclose( outfile );
  whizard_ut_end( u_log, test_name, test_description, &results );

}


