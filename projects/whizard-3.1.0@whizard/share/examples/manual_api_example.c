/* Example C code for calling WHIZARD via the C API  */
/* For detailed explanations, cf. the WHIZARD manual */
/*   
########################################################################
#
# Copyright (C) 1999-2020 by 
#     Wolfgang Kilian <kilian@physik.uni-siegen.de>
#     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
#     Juergen Reuter <juergen.reuter@desy.de>
#     with contributions from
#     cf. main AUTHORS file
#
# WHIZARD is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by 
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# WHIZARD is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
########################################################################
*/

#include <stdio.h>
#include "whizard.h"

int main( int argc, char* argv[] )
{
  /* WHIZARD and event-sample objects */
  void* wh;
  void* sample;
  
  /* Local variables */
  double integral, error;
  double sqme, weight;
  int idx;
  int it, it_begin, it_end;

  /* Initialize WHIZARD, setting some global option */
  whizard_create( &wh );
  whizard_option( &wh, "model", "QED" );
  whizard_init( &wh );
  
  /* Define a process, set some variables */
  whizard_command( &wh, "process mupair = e1, E1 => e2, E2" );
  whizard_set_double( &wh, "sqrts", 10. );
  whizard_set_int( &wh, "seed", 0 );
  
  /* Generate matrix-element code, integrate and retrieve result */
  whizard_command( &wh, "integrate (mupair)" );

  /* Print result */
  whizard_get_integration_result( &wh, "mupair", &integral, &error);
  printf( " cross section = %5.1f pb\n", integral / 1000. );
  printf( " error = %5.1f pb\n", error / 1000. );

  /* Settings for event generation */
  whizard_set_char( &wh, "$sample", "mupair_events" );
  whizard_set_int( &wh, "n_events", 2 );

  /* Create an event-sample object and generate events */
  whizard_new_sample( &wh, "mupair", &sample );
  whizard_sample_open( &sample, &it_begin, &it_end );
  for (it=it_begin; it<=it_end; it++) {
    whizard_sample_next_event( &sample );
    whizard_sample_get_event_index( &sample, &idx );
    whizard_sample_get_weight( &sample, &weight );
    whizard_sample_get_sqme( &sample, &sqme );
    printf( "Event #%d\n", idx );
    printf( " sqme = %10.3e\n", sqme );
    printf( " weight = %10.3e\n", weight );
  }

  /* Finalize the event-sample object */
  whizard_sample_close( &sample );

  /* Finalize the WHIZARD object */
  whizard_final( &wh );
}
