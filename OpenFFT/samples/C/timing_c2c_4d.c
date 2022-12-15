/*******************************************************
    This file is part of OpenFFT. 

    OpenFFT is an open-source parallel package for 3-D FFTs
    built on a communication-optimal domain decomposition method.
 
    Copyright (C) 2013-2015  Truong Vinh Truong Duy and Taisuke Ozaki, 
                             The University of Tokyo.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

******************************************************/

/**********************************************************************
  timing_c2c_4d.c:
  
   This program is used for benchmarking the performance of OpenFFT 
   with timing and GLOPS results.
   It can be executed with an arbitrary number of processes.
   Time is measured by MPI_Wtime().
   A numeric input parameter can be provided for specifying the size of 
   the 4 dimensions. If no input parameter is provided, it will be executed 
   with a default size of 32^4 data points.     
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openfft.h"

int main(int argc, char* argv[])
{ 
  dcomplex *Rhor, *Rhok;
  double time2,time;
  int numprocs,myid;
  int N1,N2,N3,N4;
  int offt_measure,measure_time,print_memory;
  int My_Max_NumGrid,My_NumGrid_In,My_NumGrid_Out;
  unsigned long long int i,l;
  double n1,flops;
  int NTEST=10;
  int My_Index_In[8],My_Index_Out[8];

  /* MPI */
  MPI_Init(&argc, &argv); 
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  /* Set global input grid size based on the program's parameters */

  if(argc==2){
    N1 = N2 = N3 = N4 = atoi(argv[1]);
  }
  else if(argc==4){
    N1 = atoi(argv[1]);
    N2 = atoi(argv[2]);
    N3 = atoi(argv[3]);
    N4 = atoi(argv[4]);
  }
  else{
    N1 = N2 = N3 = N4 = 32;
  }

  /* Select auto-tuning of communication */

  offt_measure = 0;

  /* Set whether to use the timing and print memory functions of OpenFFT 
     or not. Default=0 (not use) */

  measure_time = 0;
  print_memory = 1;

  /* Initialize OpenFFT */ 

  openfft_init_c2c_4d(N1,N2,N3,N4,
		     &My_Max_NumGrid,&My_NumGrid_In,My_Index_In,
		     &My_NumGrid_Out,My_Index_Out,
		     offt_measure,measure_time,print_memory);

  /* Allocate local input and output arrays */

  Rhor = (dcomplex*)malloc(sizeof(dcomplex)*My_Max_NumGrid); 
  Rhok = (dcomplex*)malloc(sizeof(dcomplex)*My_Max_NumGrid); 

  /* Set local input */

  for (i=0; i<My_NumGrid_In; i++){
    Rhor[i].r = 0.1; 
    Rhor[i].i = 0.2;
  }

  /* Repeat FFT transform for NTEST times */

  MPI_Barrier(MPI_COMM_WORLD);
  time2 = 0.0;
  for(l=0;l<NTEST;l++){
    time = MPI_Wtime();
    openfft_exec_c2c_4d(Rhor, Rhok);
    time2 = time2 + MPI_Wtime() - time;

    /* Re-set local input */

    for (i=0; i<My_NumGrid_In; i++){
      Rhor[i].r = 0.1; 
      Rhor[i].i = 0.2;
    }

  }

  /* Average time of all processes */

  MPI_Allreduce(&time2,&time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  time = time/NTEST;

  /* Calculate flops */

  n1 = (double) (N1*N2*N3);
  n1 = n1*N4;
  n1 = pow(n1,1.0/4.0);
  flops = 5.0*n1*log(n1)/log(2);
  flops = flops*4*pow(n1,3);
  flops = flops/time;

  /* Finalize OpenFFT */

  openfft_finalize();

  /* Free arrays */
  
  free(Rhor);
  free(Rhok);

  /* Print results */

  if(myid==0){
    printf("=======OpenFFT complex2complex========\n");
    printf("N1=%5d,N2=%5d,N3=%5d,N4=%5d,numprocs=%5d\n",N1,N2,N3,N4,numprocs);
    printf("Time openfft_execute Total = %15.12f\n",time);
    printf("Gflops                     = %15.12f\n",flops/pow(1000,3));
  }
  MPI_Finalize();

}

