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
  check_r2c_3d.c:

   This program transforms input data values to output data values.
   It can be executed with an arbitrary number of processes.
   Its input and output should match the corresponding values in check_r2c_3d.dat. 
   This program does not require any input parameter.

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <openfft.h>

int main(int argc, char* argv[])
{ 
  double *Rhor;
  dcomplex *Rhok;
  int numprocs,myid;
  int const N1=2,N2=3,N3=4;
  int N3r=N3/2+1;
  double Input[N1][N2][N3];
  dcomplex Out[N1][N2][N3r],Output[N1][N2][N3r],Output_ref[N1][N2][N3r];
  int offt_measure,measure_time,print_memory;
  int My_Max_NumGrid,My_NumGrid_In,My_NumGrid_Out;
  int i,j,k,l;
  double factor;
  int My_Index_In[6],My_Index_Out[6];

  /* MPI */
  MPI_Init(&argc, &argv); 
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  /* Set global input */

  Input[0][0][0] = 1.000;
  Input[0][0][1] = 0.999;
  Input[0][0][2] = 0.987;
  Input[0][0][3] = 0.936;
  Input[0][1][0] = 0.994;
  Input[0][1][1] = 0.989;
  Input[0][1][2] = 0.963;
  Input[0][1][3] = 0.891;
  Input[0][2][0] = 0.903;
  Input[0][2][1] = 0.885;
  Input[0][2][2] = 0.823;
  Input[0][2][3] = 0.694;
  Input[1][0][0] = 0.500;
  Input[1][0][1] = 0.499;
  Input[1][0][2] = 0.487;
  Input[1][0][3] = 0.436;
  Input[1][1][0] = 0.494;
  Input[1][1][1] = 0.489;
  Input[1][1][2] = 0.463;
  Input[1][1][3] = 0.391;
  Input[1][2][0] = 0.403;
  Input[1][2][1] = 0.385;
  Input[1][2][2] = 0.323;
  Input[1][2][3] = 0.194;

  /* Select auto-tuning of communication */

  offt_measure = 0;

  /* Set whether to use the timing and print memory functions of OpenFFT 
     or not. Default=0 (not use) */

  measure_time = 0;
  print_memory = 0;

  /* Initialize OpenFFT */ 

  openfft_init_r2c_3d(N1,N2,N3,
		     &My_Max_NumGrid,&My_NumGrid_In,My_Index_In,
		     &My_NumGrid_Out,My_Index_Out,
		     offt_measure,measure_time,print_memory);

  /* Allocate local input and output arrays */

  Rhor = (double*)malloc(sizeof(double)*My_Max_NumGrid); 
  Rhok = (dcomplex*)malloc(sizeof(dcomplex)*My_Max_NumGrid); 

  /* Set local input */

  MPI_Barrier(MPI_COMM_WORLD);

  printf("myid=%4d: Input in the ABC(XYZ) order with %d grid points ",
	 myid,My_NumGrid_In);
  if(My_NumGrid_In > 0){
    printf("from (A=%d,B=%d,C=%d) to (A=%d,B=%d,C=%d)\n",
	 My_Index_In[0],My_Index_In[1],My_Index_In[2],
	 My_Index_In[3],My_Index_In[4],My_Index_In[5]);
  }
  else{
    printf("\n");
  }

  for(i=0;i<My_Max_NumGrid;i++){
    Rhor[i] = 0.0;
    Rhok[i].r = 0.0;
    Rhok[i].i = 0.0;
  }

  l=0;

  if(My_NumGrid_In > 0){
   if(My_Index_In[0]==My_Index_In[3]){
     i=My_Index_In[0];
     for(j=My_Index_In[1];j<=My_Index_In[4];j++){
       for(k=My_Index_In[2];k<=My_Index_In[5];k++){
	 Rhor[l] = Input[i][j][k];  
	 l++;
       }
     }
   }
   else if(My_Index_In[0]<My_Index_In[3]){
     for(i=My_Index_In[0];i<=My_Index_In[3];i++){
       if(i==My_Index_In[0]){
	 for(j=My_Index_In[1];j<N2;j++){
	   for(k=My_Index_In[2];k<=My_Index_In[5];k++){
	     Rhor[l] = Input[i][j][k];  
	     l++;
	   }
	 }
       }
       else if(My_Index_In[0]<i && i<My_Index_In[3]){
	 for(j=0;j<N2;j++){
	   for(k=My_Index_In[2];k<=My_Index_In[5];k++){
	     Rhor[l] = Input[i][j][k];  
	     l++;
	   }
	 }
       }
       else if(i==My_Index_In[3]){
	 for(j=0;j<=My_Index_In[4];j++){
	   for(k=My_Index_In[2];k<=My_Index_In[5];k++){
	     Rhor[l] = Input[i][j][k];  
	     l++;
	   }
	 }
       }
     }
   }
  }

  /* Print global input */
  
  if(!myid){
    printf("Input values\n\n");
    for(i=0;i<N1;i++){
      printf("Input(i,j,k) for i = %d\n\n",i); 
      for(j=0;j<N2;j++){
	printf("Real\t"); 
	for(k=0;k<N3;k++){
	  printf("% 3.3f\t",Input[i][j][k]);
	}
	printf("\n\n");
      }
    }
  }


  /* FFT transform */

  openfft_exec_r2c_3d(Rhor, Rhok);

  /* Get local output */

  MPI_Barrier(MPI_COMM_WORLD);

  printf("myid=%4d: Output in the CBA(ZYX) order with %d grid points ",
	 myid,My_NumGrid_Out);
  if(My_NumGrid_Out > 0){
    printf("from (C=%d,B=%d,A=%d) to (C=%d,B=%d,A=%d)\n",
	   My_Index_Out[0],My_Index_Out[1],My_Index_Out[2],
	 My_Index_Out[3],My_Index_Out[4],My_Index_Out[5]);
  }
  else{
    printf("\n");
  }

  factor = sqrt(N1*N2*N3);

  for(i=0;i<N1;i++){
    for(j=0;j<N2;j++){
      for(k=0;k<N3r;k++){
	Out[i][j][k].r = 0.0;
	Out[i][j][k].i = 0.0;
	Output[i][j][k].r = 0.0;
	Output[i][j][k].i = 0.0;
      }
    }
  }

  l=0;

  if(My_NumGrid_Out > 0){
   if(My_Index_Out[0]==My_Index_Out[3]){
     i=My_Index_Out[0];
     for(j=My_Index_Out[1];j<=My_Index_Out[4];j++){
       for(k=My_Index_Out[2];k<=My_Index_Out[5];k++){
	 Out[k][j][i].r = Rhok[l].r/factor;
	 Out[k][j][i].i = Rhok[l].i/factor;
	 l++;
       }
     }
   }
   else if(My_Index_Out[0]<My_Index_Out[3]){
     for(i=My_Index_Out[0];i<=My_Index_Out[3];i++){
       if(i==My_Index_Out[0]){
	 for(j=My_Index_Out[1];j<N2;j++){
	   for(k=My_Index_Out[2];k<=My_Index_Out[5];k++){
	     Out[k][j][i].r = Rhok[l].r/factor;
	     Out[k][j][i].i = Rhok[l].i/factor;
	     l++;
	   }
	 }
       }
       else if(My_Index_Out[0]<i && i<My_Index_Out[3]){
	 for(j=0;j<N2;j++){
	   for(k=My_Index_Out[2];k<=My_Index_Out[5];k++){
	     Out[k][j][i].r = Rhok[l].r/factor;
	     Out[k][j][i].i = Rhok[l].i/factor;
	     l++;
	   }
	 }
       }
       else if(i==My_Index_Out[3]){
	 for(j=0;j<=My_Index_Out[4];j++){
	   for(k=My_Index_Out[2];k<=My_Index_Out[5];k++){
	     Out[k][j][i].r = Rhok[l].r/factor;
	     Out[k][j][i].i = Rhok[l].i/factor;
	     l++;
	   }
	 }
       }
     }
   }
  }

  /* Gather results from all processes */
  
  MPI_Allreduce(Out,Output,N1*N2*N3r,
		MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);

  /* Print global output */

  if(!myid){
    printf("Output values\n\n");
    for(i=0;i<N1;i++){
      printf("Output(i,j,k) for i = %d\n\n",i); 
      for(j=0;j<N2;j++){
	printf("Real\t"); 
	for(k=0;k<N3r;k++){
	  printf("% 3.3f\t",Output[i][j][k].r);
	}
	printf("\nImag\t"); 
	for(k=0;k<N3r;k++){
	  printf("% 3.3f\t",Output[i][j][k].i);
	}
	printf("\n\n");
      }
    }
  }

  /* Error check */

  Output_ref[0][0][0].r = 3.292; Output_ref[0][0][0].i = 0.000;
  Output_ref[0][0][1].r = 0.051; Output_ref[0][0][1].i =-0.144;
  Output_ref[0][0][2].r = 0.113; Output_ref[0][0][2].i = 0.000;
  Output_ref[0][1][0].r = 0.143; Output_ref[0][1][0].i =-0.188;
  Output_ref[0][1][1].r = 0.016; Output_ref[0][1][1].i = 0.051;
  Output_ref[0][1][2].r =-0.024; Output_ref[0][1][2].i = 0.025;
  Output_ref[0][2][0].r = 0.143; Output_ref[0][2][0].i = 0.188;
  Output_ref[0][2][1].r =-0.050; Output_ref[0][2][1].i = 0.016;
  Output_ref[0][2][2].r =-0.024; Output_ref[0][2][2].i =-0.025;
  Output_ref[1][0][0].r = 1.225; Output_ref[1][0][0].i = 0.000;
  Output_ref[1][0][1].r = 0.000; Output_ref[1][0][1].i = 0.000;
  Output_ref[1][0][2].r = 0.000; Output_ref[1][0][2].i = 0.000;
  Output_ref[1][1][0].r = 0.000; Output_ref[1][1][0].i = 0.000;
  Output_ref[1][1][1].r = 0.000; Output_ref[1][1][1].i = 0.000;
  Output_ref[1][1][2].r = 0.000; Output_ref[1][1][2].i = 0.000;
  Output_ref[1][2][0].r = 0.000; Output_ref[1][2][0].i =-0.000;
  Output_ref[1][2][1].r =-0.000; Output_ref[1][2][1].i = 0.000;
  Output_ref[1][2][2].r = 0.000; Output_ref[1][2][2].i =-0.000;

  l = 0;
  if(!myid){
    for(i=0;i<N1;i++){
      for(j=0;j<N2;j++){
	for(k=0;k<N3r;k++){
	  if(fabs(Output[i][j][k].r-Output_ref[i][j][k].r)>0.001 || 
	     fabs(Output[i][j][k].i-Output_ref[i][j][k].i)>0.001){
	    l = 1;
	    printf("ERROR Output[%d,%d,%d] (% 3.3f,% 3.3f)\n",
		   i,j,k,Output[i][j][k].r,Output[i][j][k].i);
	  }
	}
      }
    }
    if (l==0)
      printf("Check done. All output elements are correct.\n");
    else
      printf("Check done. Some output elements are incorrect.\n");
  }

  /* Free arrays */

  free(Rhor);
  free(Rhok);

  /* Finalize OpenFFT */

  openfft_finalize(MPI_COMM_WORLD);

  MPI_Finalize();

}

