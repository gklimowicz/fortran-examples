/*  The code for the command line version of acor.  See the README file for more information */

#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <math.h>

#include "acor.h"     // The header file specific for this code.  Should be in the same directory with the source.

/* The input file has unknown length so we don't know what size array to allocate.  Instead, 
   allocate arrays of size VECSIZE one at a time as they are needed.  NVECS is the maximum
   number of such vectors to allocate.  The total storage used would be VECSIZE * NVECS.  */

#define NVECS   10000    /*  This works up to 50 Million samples */   
#define VECSIZE 5000

using namespace std;


int main(int argc, char *argv[]) {


/*--------------------   Read the data file  --------------------------------------*/



// Step one: Open the input file and read the time series into short vectors.  Complain if something goes wrong.

  if (argc != 2 ) {
    cout << "The correct syntax is: acor 'filename'" << endl;
    return 1;}

   FILE *Xinput = fopen( argv[1], "r" );
      
   //   First read the numbers in to a collection of not very big vectors
   
   float xt;                      // A dummy variable because I can't pass x[k] to fscanf.
   double* vecs[NVECS];           // An array of pointers to the vectors that hold the input, 
                                  // before its size is known and it is assembled into a single
								  // vector.
   int scanfReturnValue;          // self evident.
   int L;                         // The length of the sequence, unknown until you hit the 
                                  // end of the sequence.
   int done = 0;                  // Set to a nonzero value when you hit the last number in 
								  // input file
   int i, j, k;     // loop indices, declared here because I don't understand scoping rules.
   for ( i = 0; i < NVECS; i++ ) {
      vecs[i] = new double [VECSIZE];   // Allocate the individual vectors as needed.
	  for ( j = 0; j < VECSIZE; j++) {
            scanfReturnValue =               // The scanf return value, to check for doneness.
         fscanf( Xinput, "%f", &xt);         // Read the dummy variable 
	     if ( scanfReturnValue != 1 ) {      // If done . . . 
	        L = i*VECSIZE + j;
		    done = 1; break; }
		 if ( done )  break;
		 vecs[i][j] = (double) xt; }          // Copy it to the one of the vectors, for now.
      if (done) break;}
	  
   if ( fclose( Xinput ) )
      cout << "Problem closing file" << endl;
	  
//  Step two: assemble the short vectors into one vector, X

   double* X;
   X = new double[L];
   k    = 0;
   done = 0;
   for (i=0; i < NVECS; i++) {
      for ( j = 0; j < VECSIZE; j++) {
	     X[k++] = vecs[i][j];
		 if ( k >= L ) {
		    done = 1;
			break;}
		}
      if ( done ) break; }
	  
//  Step three, return to sender

   int ngb;   // The number of vectors to give back.
   if ( (L / VECSIZE ) * VECSIZE == L ) ngb = L / VECSIZE;   // The unlikely case where you used the last vector fully.
   else ngb = ( L / VECSIZE ) + 1;                           // Probably the last vector is not full.
   
   for ( i = 0; i < ngb; i++ ) delete[] vecs[i];
   
   
/*-------------------   do the math  ---------------------------------*/

   double Xbar, sigma, tau;   // The numbers we came for.
   acor( &Xbar, &sigma, &tau, X, L);
   
   cout << "sample mean = "          << Xbar  << 
		",  standard deviation = "   << sigma << 
		",  autocorrelation time = " << tau   << 
		", series length = " << L << endl;
   
   
   return 0;
  }
