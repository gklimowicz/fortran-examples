// Example C++ code for calling WHIZARD via the C++ API
// For detailed explanations, cf. the WHIZARD manual
// 
// ########################################################################
// #
// # Copyright (C) 1999-2020 by 
// #     Wolfgang Kilian <kilian@physik.uni-siegen.de>
// #     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
// #     Juergen Reuter <juergen.reuter@desy.de>
// #     with contributions from
// #     cf. main AUTHORS file
// #
// # WHIZARD is free software; you can redistribute it and/or modify it
// # under the terms of the GNU General Public License as published by 
// # the Free Software Foundation; either version 2, or (at your option)
// # any later version.
// #
// # WHIZARD is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the implied warranty of
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
// # GNU General Public License for more details.
// #
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// #
// ########################################################################

#include <cstdio>
#include <string>
#include "whizard.h"

int main( int argc, char* argv[] )
{
  // WHIZARD and event-sample objects
  Whizard* whizard;
  WhizardSample* sample;

  // Local variables
  double integral, error;
  double sqme, weight;
  int idx;
  int it, it_begin, it_end;

  // Initialize WHIZARD, setting some global option
  whizard = new Whizard();
  whizard->option( "model", "QED" );
  whizard->init();

  // Define a process, set some variables
  whizard->command( "process mupair = e1, E1 => e2, E2" );
  whizard->set_double( "sqrts", 10. );
  whizard->set_int( "seed", 0 );

  // Generate matrix-element code, integrate and retrieve result
  whizard->command( "integrate (mupair)" );

  // Print result
  whizard->get_integration_result( "mupair", &integral, &error );
  printf( " cross section = %5.1f pb\n", integral / 1000. );
  printf( " error = %5.1f pb\n", error / 1000. );

  // Settings for event generation
  whizard->set_string( "$sample", "mupair_events" );
  whizard->set_int( "n_events", 2 );

  // Create an event-sample object and generate events
  sample = whizard->new_sample( "mupair" );
  sample->open( &it_begin, &it_end );
  for (it=it_begin; it<=it_end; it++) {
    sample->next_event();
    idx = sample->get_event_index();
    weight = sample->get_weight();
    sqme = sample->get_sqme();
    printf( "Event #%d\n", idx );
    printf( " sqme = %10.3e\n", sqme );
    printf( " weight = %10.3e\n", weight );
  }

  // Finalize the event-sample object
  sample->close();
  delete sample;

  // Finalize the WHIZARD object
  delete whizard;
}
