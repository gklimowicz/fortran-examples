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
  check_c2c_4d.c:

   This program transforms input data values to output data values.
   It can be executed with an arbitrary number of processes.
   Its input and output should match the corresponding values in check_c2c_4d.dat. 
   This program does not require any input parameter.

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <openfft.h>

#define RANGE 4.5
 
int main(int argc, char* argv[])
{ 
  dcomplex *Rhor, *Rhok;
  int numprocs,myid;
  int const N1=3,N2=4,N3=5,N4=6;
  dcomplex Input[N1][N2][N3][N4],Output[N1][N2][N3][N4];
  dcomplex Out[N1][N2][N3][N4],Output_ref[N1][N2][N3][N4];
  dcomplex In[N1][N2][N3][N4];
  int offt_measure,measure_time,print_memory;
  int My_Max_NumGrid,My_NumGrid_In,My_NumGrid_Out,Tot_In,Tot_Out;
  int i,j,k,l,m,x;
  double factor;
  int My_Index_In[8],My_Index_Out[8];
  FILE *file_in;
  char BUF[1000];

  /* MPI */
  MPI_Init(&argc, &argv); 
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if(!myid) printf("Executing with %d processes\n",numprocs);

  /* Set global input */

  file_in = fopen("check_c2c_4d.din","r");
  fgets(BUF,sizeof(BUF),file_in);
  fgets(BUF,sizeof(BUF),file_in);
  fgets(BUF,sizeof(BUF),file_in);
  fgets(BUF,sizeof(BUF),file_in);
  for(m=0;m<N1;m++){
    for(i=0;i<N2;i++){
      for(j=0;j<N3;j++){
	for(k=0;k<N4;k++){
	  fgets(BUF,sizeof(BUF),file_in);
	  sscanf(BUF,"%lf  %lf\n",&Input[m][i][j][k].r,&Input[m][i][j][k].i);
	}
      }
    }
  }
  fclose(file_in);

  /* Select auto-tuning of communication */

  if(argc==2){
    offt_measure = atoi(argv[1]);
  }
  else{
    offt_measure = 0;
  }

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

  MPI_Barrier(MPI_COMM_WORLD);

  for(i=0;i<My_Max_NumGrid;i++){
    Rhor[i].r = 0.0;
    Rhor[i].i = 0.0;
    Rhok[i].r = 0.0;
    Rhok[i].i = 0.0;
  }

  x=0;

  if(My_NumGrid_In > 0){
    if(My_Index_In[0]==My_Index_In[4]){
      m = My_Index_In[0];
      if(My_Index_In[1]==My_Index_In[5]){
	i=My_Index_In[1];
	for(j=My_Index_In[2];j<=My_Index_In[6];j++){
	  for(k=My_Index_In[3];k<=My_Index_In[7];k++){
	    Rhor[x].r = Input[m][i][j][k].r;  
	    Rhor[x].i = Input[m][i][j][k].i;  
	    x++;
	  }
	}
      }
      else if(My_Index_In[1]<My_Index_In[5]){
	for(i=My_Index_In[1];i<=My_Index_In[5];i++){
	  if(i==My_Index_In[1]){
	    for(j=My_Index_In[2];j<N3;j++){
	      for(k=My_Index_In[3];k<=My_Index_In[7];k++){
		Rhor[x].r = Input[m][i][j][k].r;  
		Rhor[x].i = Input[m][i][j][k].i;  
		x++;
	      }
	    }
	  }
	  else if(My_Index_In[1]<i && i<My_Index_In[5]){
	    for(j=0;j<N3;j++){
	      for(k=My_Index_In[3];k<=My_Index_In[7];k++){
		Rhor[x].r = Input[m][i][j][k].r;  
		Rhor[x].i = Input[m][i][j][k].i;  
		x++;
	      }
	    }
	  }
	  else if(i==My_Index_In[5]){
	    for(j=0;j<=My_Index_In[6];j++){
	      for(k=My_Index_In[3];k<=My_Index_In[7];k++){
		Rhor[x].r = Input[m][i][j][k].r;  
		Rhor[x].i = Input[m][i][j][k].i;  
		x++;
	      }
	    }
	  }
	}
      }
    } /* if */

    else if(My_Index_In[0]<My_Index_In[4]){
      for(m=My_Index_In[0];m<=My_Index_In[4];m++){
	if(m==My_Index_In[0]){
	  for(i=My_Index_In[1];i<N2;i++){
	    if(i==My_Index_In[1]){
	      for(j=My_Index_In[2];j<N3;j++){
		for(k=My_Index_In[3];k<=My_Index_In[7];k++){
		  Rhor[x].r = Input[m][i][j][k].r;  
		  Rhor[x].i = Input[m][i][j][k].i;  
		  x++;
		}
	      }
	    }
	    else{
	      for(j=0;j<N3;j++){
		for(k=My_Index_In[3];k<=My_Index_In[7];k++){
		  Rhor[x].r = Input[m][i][j][k].r;  
		  Rhor[x].i = Input[m][i][j][k].i;  
		  x++;
		}
	      }
	    }
	  }
	}
	else if(My_Index_In[0]<m && m<My_Index_In[4]){
	  for(i=0;i<N2;i++){
	    for(j=0;j<N3;j++){
	      for(k=My_Index_In[3];k<=My_Index_In[7];k++){
		Rhor[x].r = Input[m][i][j][k].r;  
		Rhor[x].i = Input[m][i][j][k].i;  
		x++;
	      }
	    }
	  }
	}
	else if(m==My_Index_In[4]){
	  for(i=0;i<=My_Index_In[5];i++){
	    if(i==My_Index_In[5]){
	      for(j=0;j<=My_Index_In[6];j++){
		for(k=My_Index_In[3];k<=My_Index_In[7];k++){
		  Rhor[x].r = Input[m][i][j][k].r;  
		  Rhor[x].i = Input[m][i][j][k].i;  
		  x++;
		}
	      }
	    }
	    else{
	      for(j=0;j<N3;j++){
		for(k=My_Index_In[3];k<=My_Index_In[7];k++){
		  Rhor[x].r = Input[m][i][j][k].r;  
		  Rhor[x].i = Input[m][i][j][k].i;  
		  x++;
		}
	      }
	    }
	  }
	}
      } 
    }/* else */


  } /* Gridin */


  printf("myid=%4d: Input in the ABCD(XYZU) order with %d grid points %d ",
	 myid,My_NumGrid_In,x);
  if(My_NumGrid_In > 0){
    printf("from (A=%d,B=%d,C=%d,D=%d) to (A=%d,B=%d,C=%d,D=%d)\n",
	   My_Index_In[0],My_Index_In[1],My_Index_In[2],My_Index_In[3],
	   My_Index_In[4],My_Index_In[5],My_Index_In[6],My_Index_In[7]);
  }
  else{
    printf("\n");
  }

  if(My_NumGrid_In != x){
    printf("ERROR\n");
  }

  MPI_Allreduce(&x,&Tot_In,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  if(Tot_In != N1*N2*N3*N4){
    printf("ERROR\n");
  } 

  /* FFT transform */

  openfft_exec_c2c_4d(Rhor, Rhok);

  /* Get local output */

  MPI_Barrier(MPI_COMM_WORLD);

  factor = sqrt(N1*N2*N3*N4);

  for(m=0;m<N1;m++){
    for(i=0;i<N2;i++){
      for(j=0;j<N3;j++){
	for(k=0;k<N4;k++){
	  Out[m][i][j][k].r = 0.0;
	  Out[m][i][j][k].i = 0.0;
	  Output[m][i][j][k].r = 0.0;
	  Output[m][i][j][k].i = 0.0;
	}
      }
    }
  }

  x=0;

  if(My_NumGrid_Out > 0){
    if(My_Index_Out[0]==My_Index_Out[4]){
      m = My_Index_Out[0];
      if(My_Index_Out[1]==My_Index_Out[5]){
	i=My_Index_Out[1];
	for(j=My_Index_Out[2];j<=My_Index_Out[6];j++){
	  for(k=My_Index_Out[3];k<=My_Index_Out[7];k++){
	    Out[k][j][i][m].r = Rhok[x].r/factor;
	    Out[k][j][i][m].i = Rhok[x].i/factor;
	    x++;
	  }
	}
      }
      else if(My_Index_Out[1]<My_Index_Out[5]){
	for(i=My_Index_Out[1];i<=My_Index_Out[5];i++){
	  if(i==My_Index_Out[1]){
	    for(j=My_Index_Out[2];j<N2;j++){
	      for(k=My_Index_Out[3];k<=My_Index_Out[7];k++){
		Out[k][j][i][m].r = Rhok[x].r/factor;
		Out[k][j][i][m].i = Rhok[x].i/factor;
		x++;
	      }
	    }
	  }
	  else if(My_Index_Out[1]<i && i<My_Index_Out[5]){
	    for(j=0;j<N2;j++){
	      for(k=My_Index_Out[3];k<=My_Index_Out[7];k++){
		Out[k][j][i][m].r = Rhok[x].r/factor;
		Out[k][j][i][m].i = Rhok[x].i/factor;
		x++;
	      }
	    }
	  }
	  else if(i==My_Index_Out[5]){
	    for(j=0;j<=My_Index_Out[6];j++){
	      for(k=My_Index_Out[3];k<=My_Index_Out[7];k++){
		Out[k][j][i][m].r = Rhok[x].r/factor;
		Out[k][j][i][m].i = Rhok[x].i/factor;
		x++;
	      }
	    }
	  }
	}
      }
    } /* if */

    else if(My_Index_Out[0]<My_Index_Out[4]){
      for(m=My_Index_Out[0];m<=My_Index_Out[4];m++){
	if(m==My_Index_Out[0]){
	  for(i=My_Index_Out[1];i<N3;i++){
	    if(i==My_Index_Out[1]){
	      for(j=My_Index_Out[2];j<N2;j++){
		for(k=My_Index_Out[3];k<=My_Index_Out[7];k++){
		  Out[k][j][i][m].r = Rhok[x].r/factor;
		  Out[k][j][i][m].i = Rhok[x].i/factor;
		  x++;
		}
	      }
	    }
	    else{
	      for(j=0;j<N2;j++){
		for(k=My_Index_Out[3];k<=My_Index_Out[7];k++){
		  Out[k][j][i][m].r = Rhok[x].r/factor;
		  Out[k][j][i][m].i = Rhok[x].i/factor;
		  x++;
		}
	      }
	    }
	  }
	}
	else if(My_Index_Out[0]<m && m<My_Index_Out[4]){
	  for(i=0;i<N3;i++){
	    for(j=0;j<N2;j++){
	      for(k=My_Index_Out[3];k<=My_Index_Out[7];k++){
		Out[k][j][i][m].r = Rhok[x].r/factor;
		Out[k][j][i][m].i = Rhok[x].i/factor;
		x++;
	      }
	    }
	  }
	}
	else if(m==My_Index_Out[4]){
	  for(i=0;i<=My_Index_Out[5];i++){
	    if(i==My_Index_Out[5]){
	      for(j=0;j<=My_Index_Out[6];j++){
		for(k=My_Index_Out[3];k<=My_Index_Out[7];k++){
		  Out[k][j][i][m].r = Rhok[x].r/factor;
		  Out[k][j][i][m].i = Rhok[x].i/factor;
		  x++;
		}
	      }
	    }
	    else{
	      for(j=0;j<N2;j++){
		for(k=My_Index_Out[3];k<=My_Index_Out[7];k++){
		  Out[k][j][i][m].r = Rhok[x].r/factor;
		  Out[k][j][i][m].i = Rhok[x].i/factor;
		  x++;
		}
	      }
	    }
	  }
	}
      } 
    }/* else */

  } /* Gridout */


  printf("myid=%4d: Output in the DCBA(UZYX) order with %d grid points %d ",
	 myid,My_NumGrid_Out,x);
  if(My_NumGrid_Out > 0){
    printf("from (D=%d,C=%d,B=%d,A=%d) to (D=%d,C=%d,B=%d,A=%d)\n",
	   My_Index_Out[0],My_Index_Out[1],My_Index_Out[2],My_Index_Out[3],
	   My_Index_Out[4],My_Index_Out[5],My_Index_Out[6],My_Index_Out[7]);
  }
  else{
    printf("\n");
  }

  if(My_NumGrid_Out != x){
    printf("ERROR\n");
  }

  MPI_Allreduce(&x,&Tot_Out,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  if(Tot_Out != N1*N2*N3*N4){
    printf("ERROR\n");
  } 

  /* Gather results from all processes */

  MPI_Allreduce(Out,Output,N1*N2*N3*N4,
		MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);

  /* Print global output */

  if(!myid){
    sprintf(BUF,"check_c2c_4dx%d.dout",numprocs);
    file_in = fopen(BUF,"w");
    for(m=0;m<N1;m++){
      for(i=0;i<N2;i++){
	for(j=0;j<N3;j++){
	  for(k=0;k<N4;k++){
	    fprintf(file_in,"%10.3f  %10.3f\n",
		    Output[m][i][j][k].r,Output[m][i][j][k].i);
	  }
	}
      }
    }
    fclose(file_in);
  }

  
  if(!myid){
    file_in = fopen("check_c2c_4d.dout","r");
    for(m=0;m<N1;m++){
      for(i=0;i<N2;i++){
	for(j=0;j<N3;j++){
	  for(k=0;k<N4;k++){
	    fgets(BUF,sizeof(BUF),file_in);
	    sscanf(BUF,"%lf  %lf\n",
		   &Output_ref[m][i][j][k].r,&Output_ref[m][i][j][k].i);
	  }
	}
      }
    }
    fclose(file_in);
  }

  /* Error check */

  l = 0;
  if(!myid){
    for(m=0;m<N1;m++){
      for(i=0;i<N2;i++){
	for(j=0;j<N3;j++){
	  for(k=0;k<N4;k++){
	    if(fabs(Output[m][i][j][k].r-Output_ref[m][i][j][k].r)>0.001 || 
	       fabs(Output[m][i][j][k].i-Output_ref[m][i][j][k].i)>0.001){
	      l = 1;
	      printf("ERROR Output[%d,%d,%d,%d] (% 3.3f,% 3.3f)\n",
		     m,i,j,k,Output[m][i][j][k].r,Output[m][i][j][k].i);
	    }
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

