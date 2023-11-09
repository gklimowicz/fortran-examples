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
  breaktime_c2c_4d.c:
  
   This program is used for benchmarking the performance of OpenFFT 
   with timing and GLOPS results. Total elapsed time is broken down into 
   several parts.
   It can be executed with an arbitrary number of processes.
   Time is measured by OpenFFT.   
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
  double timers[10],time1,time3;
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

  /* Select communication pattern number 6 for a correct breakdown of time */

  offt_measure = 6;

  /* Set whether to use the timing and print memory functions of OpenFFT 
     or not. Default=0 (not use) */

  measure_time = 1;
  print_memory = 1;

  /* Initialize OpenFFT */ 

  time1 = openfft_init_c2c_4d(N1,N2,N3,N4,
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

  for(i=0;i<10;i++) timers[i] = 0.0;

  for(l=0;l<NTEST;l++){
    timers[0] += openfft_exec_c2c_4d(Rhor, Rhok);
    for(i=1;i<=6;i++) timers[i] += OFFT_timers[i];

    /* Re-set local input */

    for (i=0; i<My_NumGrid_In; i++){
      Rhor[i].r = 0.1; 
      Rhor[i].i = 0.2;
    }

  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  /* Average time of all processes */

  for(i=0;i<=6;i++) timers[i] /= NTEST;

  /* Calculate flops */

  n1 = (double) (N1*N2*N3);
  n1 = n1*N4;
  n1 = pow(n1,1.0/4.0);
  flops = 5.0*n1*log(n1)/log(2);
  flops = flops*4*pow(n1,3);
  flops = flops/timers[0];

  /* Finalize OpenFFT */

  time3 = openfft_finalize();

  /* Free arrays */

  free(Rhor);
  free(Rhok);

  /* Print results */

  if(myid==0){
    printf("=======OpenFFT complex2complex========\n");
    printf("N1=%5d,N2=%5d,N3=%5d,N4=%d,numprocs=%5d\n",N1,N2,N3,N4,numprocs);
    printf("Time openfft_execute Total  = %15.12f\n",timers[0]);
    printf("Time openfft_execute Comm   = %15.12f\n",timers[1]);
    printf("Time openfft_execute FFT    = %15.12f\n",timers[2]);
    printf("Time openfft_execute Copy   = %15.12f\n",timers[3]);
    printf("Time openfft_execute Others = %15.12f\n",timers[4]);
    printf("Time openfft_execute Copy F = %15.12f\n",timers[5]);
    printf("Time openfft_execute Copy C = %15.12f\n",timers[6]);
    printf("Gflops                      = %15.12f\n",flops/pow(1000,3));
  }

  MPI_Finalize();

}

