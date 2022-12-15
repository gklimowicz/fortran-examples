#include <iostream>
#include <string>
#include <cstdio>

using namespace std;

/*   Create a simple time series with nonzero mean and auotcovariance function C(t) = C(0)*a^t with a=.6 .

    Usage:  acorTest l file 
    here acorTest   = name of the program
	     l          = the length of the sequence produced
		 file       = the name of the file where the numbers are written.  It is overwritten if it is there already.
		 
   The output file has the format needed for the acor programs: one floating point number per line.  
   		 
	Jonathan Goodman, goodman@cims.nyu.edu, March, 2009.
		 
*/

int main(int argc, char *argv[]) {

// Open the input file and complain if something goes wrong.

  if (argc != 3 ) {
    cout << "The correct syntax is: ms length 'filename'" << endl;
    return 1;}

  FILE *Xoutput = fopen( argv[2], "w" ) ;
  
  char* ls = argv[1];   // The first argument is l, which is a character rather than an integer ... 
  int l = atoi( ls );   //   ... must convert.

   
  double a = .9;
  double x =  0;
  srandom(72437);
  for ( int i = 0; i < l; i++) {
    x = a*x + random()/( (double) RAND_MAX );
    fprintf( Xoutput, " %10.4f\n", x); }
   

  if ( fclose( Xoutput ) )
    cout << "Problem closing file" << endl;
   
  return 0;
}