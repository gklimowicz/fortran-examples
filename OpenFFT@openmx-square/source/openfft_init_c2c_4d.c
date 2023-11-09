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

/*******************************************************
  openfft_init_c2c_4d prepares data for the c2c transformation of 4D data.
ABCD -> ABDC -> DCAB -> DCBA 

1. Input: 

 - 4 dimensions of data: N1, N2, N3, N4.

 - offf_measure: for choosing the communication performing method.

   + 0: auto-tuning of communication, where OpenFFT automatically performs tests 
        with all of the following patterns and picks the best performer in run time 
        (recommended for high performance).

   + 1: MPI_Alltoallv.

   + 2: MPI_Isend and MPI_Irecv within sub-groups of processes.

   + 3: MPI_Isend and MPI_Irecv with communication-computation overlap. 

   + 4: MPI_Isend and MPI_Irecv within sub-groups of processes 
        with communication-computation overlap.

   + 5: MPI_Sendrecv.

   + 6: MPI_Isend and MPI_Irecv.

   + Others: default communication, which is 3.

 - measure_time and print_memory (0: disabled (default), 1: enabled).

2. Output: arrays allocated and variables initialized.

 - My_Max_NumGrid: the maximum number of grid points allocated to a process.

 - My_NumGrid_In : the number of grid points allocated to a process 
                   for the first transformation.

 - My_Index_In   : the indexes of grid points allocated to a process
                   for the first transformation.

 - My_NumGrid_Out: the number of grid points allocated to a process 
                   for the last transformation.

 - My_Index_Out  : the indexes of grid points allocated to a process 
                   for the last transformation.

******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openfft.h"
#include <values.h>

double openfft_init_c2c_4d(int N1, int N2, int N3, int N4,
			   int *My_Max_NumGrid, 
			   int *My_NumGrid_In, int *My_Index_In,
			   int *My_NumGrid_Out, int *My_Index_Out,
			   int offt_measure,
			   int measure_time, int print_memory)
{

  int numprocs,myid,tag=999,ID,IDS,IDR;
  long long int BN_ABC,BN_DCA,BN_ABD,BN_DCB;
  long long int GNs,n3D,N3D;
  long long int GN_ABC,GN_DCA,GN_ABD,GN_DCB;
  long long int n1,n2,n3,n4,blocksize,i,comm_patt_no,l;
  MPI_Status stat;
  MPI_Request request;
  double time0,TStime,TEtime;
  long long int myStart,myEnd;
  double MB,time,time1,time_min;
  long long int size_array0,size_array1,size_array;
  MPI_Comm mpi_comm_level1 = MPI_COMM_WORLD;
  dcomplex *Rhor, *Rhok;
  int numthreads,mytid,nt;
  int N_A[] = {N1}, N_B[] = {N2}, N_C[] = {N3}, N_D[] = {N4};
  int *inembed_A = N_A, *onembed_A = N_A;
  int *inembed_B = N_B, *onembed_B = N_B;
  int *inembed_C = N_C, *onembed_C = N_C;
  int *inembed_D = N_D, *onembed_D = N_D;
  
  time0 = 0.0;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if(!myid){
    printf("\n");
    printf("\n============================================================");
    printf("\n   OpenFFT: An open-source parallel package for 3-D FFTs    ");
    printf("\n   built on a communication-optimal decomposition method.   ");
    printf("\n============================================================");
    printf("\n");
    printf("\n=================OpenFFT1.2 c2c 4-D Interface===============");
    printf("\n");
  } 

  OFFT_MEASURE = offt_measure;
  OFFT_measure_time = measure_time;
  OFFT_print_memory = print_memory;

  if(OFFT_measure_time){
      MPI_Barrier(mpi_comm_level1);
      openfft_dtime(&TStime);
  }

  OFFT_alloc_first = 1;
  OFFT_alloc_first1 = 1;
  OFFT_alloc_first2 = 1;

  OFFT_Ngrid1 = N1;
  OFFT_Ngrid2 = N2;
  OFFT_Ngrid3 = N3;
  OFFT_Ngrid4 = N4;

  Max_OneD_Grids = OFFT_Ngrid1 > OFFT_Ngrid2 ? OFFT_Ngrid1 : OFFT_Ngrid2;
  Max_OneD_Grids =  Max_OneD_Grids > OFFT_Ngrid3 ?  Max_OneD_Grids : OFFT_Ngrid3;
  Max_OneD_Grids =  Max_OneD_Grids > OFFT_Ngrid4 ?  Max_OneD_Grids : OFFT_Ngrid4;

  numthreads = omp_get_max_threads();

  /****************************************************************
      ABC to ABD for MPI communication  
  ****************************************************************/

  /* set My_NumGrid_ABC */

  N3D = OFFT_Ngrid1*OFFT_Ngrid2*OFFT_Ngrid3;
  My_NumGrid_ABC =  (((myid+1)*N3D+numprocs-1)/numprocs)*OFFT_Ngrid4
                  - ((myid*N3D+numprocs-1)/numprocs)*OFFT_Ngrid4;
  OFFT_My_NumGrid_In = My_NumGrid_ABC; 

  /* set myStart,myEnd,OFFT_My_Index_In */

  /* 1D */
  N3D = OFFT_Ngrid1*OFFT_Ngrid2*OFFT_Ngrid3;
  myStart = ((myid*N3D+numprocs-1)/numprocs)*OFFT_Ngrid4;
  myEnd   = (((myid+1)*N3D+numprocs-1)/numprocs)*OFFT_Ngrid4 - 1;

  /* 3D */
  N3D = OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4;
  
  My_Index_In[0] = myStart/N3D;
  My_Index_In[1] = (myStart-My_Index_In[0]*N3D)/(OFFT_Ngrid3*OFFT_Ngrid4);
  My_Index_In[2] = (myStart-My_Index_In[0]*N3D-My_Index_In[1]*OFFT_Ngrid3*OFFT_Ngrid4)/OFFT_Ngrid4;
  My_Index_In[3] = myStart-My_Index_In[0]*N3D-My_Index_In[1]*OFFT_Ngrid3*OFFT_Ngrid4-My_Index_In[2]*OFFT_Ngrid4;

  My_Index_In[4] = myEnd/N3D;
  My_Index_In[5] = (myEnd-My_Index_In[4]*N3D)/(OFFT_Ngrid3*OFFT_Ngrid4);
  My_Index_In[6] = (myEnd-My_Index_In[4]*N3D-My_Index_In[5]*OFFT_Ngrid3*OFFT_Ngrid4)/OFFT_Ngrid4;
  My_Index_In[7] = myEnd-My_Index_In[4]*N3D-My_Index_In[5]*OFFT_Ngrid3*OFFT_Ngrid4-My_Index_In[6]*OFFT_Ngrid4;
  
  /* set GNs (starting point in 1D ABCD) */

  N3D = OFFT_Ngrid1*OFFT_Ngrid2*OFFT_Ngrid3;
  GNs = ((myid*N3D+numprocs-1)/numprocs)*OFFT_Ngrid4;

  /* find Num_Snd_Grid_ABC2ABD */

  Num_Snd_Grid_ABC2ABD = (int*)malloc(sizeof(int)*numprocs);
  Num_Rcv_Grid_ABC2ABD = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_ABC2ABD[ID] = 0;

  N3D = OFFT_Ngrid1*OFFT_Ngrid2*OFFT_Ngrid4; /* for ABDC */

  for (BN_ABC=0; BN_ABC<My_NumGrid_ABC; BN_ABC++){

    /* position in 1D ABCD */
    GN_ABC = BN_ABC + GNs;

    /* position in 4D ABCD n1n2n3n4*/
    n1 = GN_ABC/(OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4);
    n2 = (GN_ABC - n1*OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4)/(OFFT_Ngrid3*OFFT_Ngrid4);
    n3 = (GN_ABC - n1*OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4-n2*OFFT_Ngrid3*OFFT_Ngrid4)/OFFT_Ngrid4;
    n4 = GN_ABC - n1*OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4-n2*OFFT_Ngrid3*OFFT_Ngrid4-n3*OFFT_Ngrid4;

    /* n1n2n4n3 (ABDC) */
    n3D = n1*OFFT_Ngrid2*OFFT_Ngrid4 + n2*OFFT_Ngrid4 + n4;
    ID = n3D*numprocs/N3D;
    Num_Snd_Grid_ABC2ABD[ID]++;
  }

  /* MPI: Num_Snd_Grid_ABC2ABD */  
  
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Snd_Grid_ABC2ABD[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Rcv_Grid_ABC2ABD[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* allocation of arrays */  

  if (OFFT_alloc_first==0){

    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_ABC2ABD[ID]);
    }  
    free(Index_Snd_Grid_ABC2ABD);
  
    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_ABC2ABD[ID]);
    }  
    free(Index_Rcv_Grid_ABC2ABD);

    free(ID_NN_ABC2ABD_S);
    free(ID_NN_ABC2ABD_R);
    free(GP_ABC2ABD_S);
    free(GP_ABC2ABD_R);
    free(GP_ABC2ABD_S_A);
    free(GP_ABC2ABD_R_A);
  }

  size_Index_Snd_Grid_ABC2ABD = 0;
  Index_Snd_Grid_ABC2ABD = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_ABC2ABD[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid_ABC2ABD[ID]);
    size_Index_Snd_Grid_ABC2ABD += Num_Snd_Grid_ABC2ABD[ID];
  }  
  
  size_Index_Rcv_Grid_ABC2ABD = 0; 
  Index_Rcv_Grid_ABC2ABD = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_ABC2ABD[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid_ABC2ABD[ID]);
    size_Index_Rcv_Grid_ABC2ABD += Num_Rcv_Grid_ABC2ABD[ID]; 
  }  

  OFFT_alloc_first = 0;

  /* construct Index_Snd_Grid_ABC2ABD */  

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_ABC2ABD[ID] = 0;

  N3D = OFFT_Ngrid1*OFFT_Ngrid2*OFFT_Ngrid4; /* for ABDC */

  for (BN_ABC=0; BN_ABC<My_NumGrid_ABC; BN_ABC++){

    /* position in 1D ABCD */
    GN_ABC = BN_ABC + GNs;

    /* position in 3D ABCD n1n2n3n4*/
    n1 = GN_ABC/(OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4);
    n2 = (GN_ABC - n1*OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4)/(OFFT_Ngrid3*OFFT_Ngrid4);
    n3 = (GN_ABC - n1*OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4-n2*OFFT_Ngrid3*OFFT_Ngrid4)/OFFT_Ngrid4;
    n4 = GN_ABC - n1*OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4-n2*OFFT_Ngrid3*OFFT_Ngrid4-n3*OFFT_Ngrid4;

    /* n1n2n4n3 (ABDC) */
    n3D = n1*OFFT_Ngrid2*OFFT_Ngrid4 + n2*OFFT_Ngrid4 + n4;
    ID = n3D*numprocs/N3D;
    
    /* position in 1D ABDC */
    GN_ABD = n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3 + n2*OFFT_Ngrid4*OFFT_Ngrid3 + n4*OFFT_Ngrid3 + n3;
    /* offset from the starting point in 1D ABDC */
    BN_ABD = GN_ABD - ((ID*N3D+numprocs-1)/numprocs)*OFFT_Ngrid3;
    Index_Snd_Grid_ABC2ABD[ID][Num_Snd_Grid_ABC2ABD[ID]] = BN_ABD;
    Num_Snd_Grid_ABC2ABD[ID]++;
  }

  /* MPI: Index_Snd_Grid_ABC2ABD */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
   
    if (Num_Snd_Grid_ABC2ABD[IDS]!=0){
      MPI_Isend( &Index_Snd_Grid_ABC2ABD[IDS][0], Num_Snd_Grid_ABC2ABD[IDS], 
		 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Rcv_Grid_ABC2ABD[IDR]!=0){
      MPI_Recv( &Index_Rcv_Grid_ABC2ABD[IDR][0], Num_Rcv_Grid_ABC2ABD[IDR], 
		MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Snd_Grid_ABC2ABD[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* reset: Index_Snd_Grid_ABC2ABD

  Index_Snd_Grid_ABC2ABD:  BN_ABC
  Index_Rcv_Grid_ABC2ABD:  BN_ABD
  */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_ABC2ABD[ID] = 0;

  N3D = OFFT_Ngrid1*OFFT_Ngrid2*OFFT_Ngrid4; /* for ABDC */

  for (BN_ABC=0; BN_ABC<My_NumGrid_ABC; BN_ABC++){

    /* position in 1D ABCD */
    GN_ABC = BN_ABC + GNs;

    /* position in 4D ABCD n1n2n3n4*/
    n1 = GN_ABC/(OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4);
    n2 = (GN_ABC - n1*OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4)/(OFFT_Ngrid3*OFFT_Ngrid4);
    n3 = (GN_ABC - n1*OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4-n2*OFFT_Ngrid3*OFFT_Ngrid4)/OFFT_Ngrid4;
    n4 = GN_ABC - n1*OFFT_Ngrid2*OFFT_Ngrid3*OFFT_Ngrid4-n2*OFFT_Ngrid3*OFFT_Ngrid4-n3*OFFT_Ngrid4;

    /* n1n2n4n3 (ABDC) */
    n3D = n1*OFFT_Ngrid2*OFFT_Ngrid4 + n2*OFFT_Ngrid4 + n4;
    ID = n3D*numprocs/N3D;
    Index_Snd_Grid_ABC2ABD[ID][Num_Snd_Grid_ABC2ABD[ID]] = BN_ABC;
    Num_Snd_Grid_ABC2ABD[ID]++;
  }

  /* find the maximum Num_Snd_Grid_ABC2ABD and Num_Rcv_Grid_ABC2ABD */

  Max_Num_Snd_Grid_ABC2ABD = 0;
  Max_Num_Rcv_Grid_ABC2ABD = 0;

  for (ID=0; ID<numprocs; ID++){
    if (Max_Num_Snd_Grid_ABC2ABD<Num_Snd_Grid_ABC2ABD[ID]){ 
      Max_Num_Snd_Grid_ABC2ABD = Num_Snd_Grid_ABC2ABD[ID];
    }

    if (Max_Num_Rcv_Grid_ABC2ABD<Num_Rcv_Grid_ABC2ABD[ID]){ 
      Max_Num_Rcv_Grid_ABC2ABD = Num_Rcv_Grid_ABC2ABD[ID];
    }
  }

  /* find NN_ABC2ABD_S and NN_ABC2ABD_R 
     and set ID_NN_ABC2ABD_S, 
             ID_NN_ABC2ABD_R,
             GP_ABC2ABD_S,
             GP_ABC2ABD_R
  */

  NN_ABC2ABD_S = 0;
  NN_ABC2ABD_R = 0;
 
  for (ID=0; ID<numprocs; ID++){
    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;
    if (Num_Snd_Grid_ABC2ABD[IDS]!=0) NN_ABC2ABD_S++;
    if (Num_Rcv_Grid_ABC2ABD[IDR]!=0) NN_ABC2ABD_R++;
  }

  ID_NN_ABC2ABD_S = (int*)malloc(sizeof(int)*NN_ABC2ABD_S);
  ID_NN_ABC2ABD_R = (int*)malloc(sizeof(int)*NN_ABC2ABD_R);
  GP_ABC2ABD_S = (int*)malloc(sizeof(int)*(NN_ABC2ABD_S+1));
  GP_ABC2ABD_R = (int*)malloc(sizeof(int)*(NN_ABC2ABD_R+1));
  GP_ABC2ABD_S_A = (int*)malloc(sizeof(int)*numprocs);
  GP_ABC2ABD_R_A = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++){
    GP_ABC2ABD_S_A[ID] = -1;
    GP_ABC2ABD_R_A[ID] = -1;
  }

  NN_ABC2ABD_S = 0;
  NN_ABC2ABD_R = 0;
  GP_ABC2ABD_S[0] = 0;
  GP_ABC2ABD_R[0] = 0;

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_ABC2ABD[IDS]!=0){
      ID_NN_ABC2ABD_S[NN_ABC2ABD_S] = IDS;
      NN_ABC2ABD_S++;
      GP_ABC2ABD_S[NN_ABC2ABD_S] = GP_ABC2ABD_S[NN_ABC2ABD_S-1] + Num_Snd_Grid_ABC2ABD[IDS];
      GP_ABC2ABD_S_A[IDS] = GP_ABC2ABD_S[NN_ABC2ABD_S-1];
    }

    if (Num_Rcv_Grid_ABC2ABD[IDR]!=0){
      ID_NN_ABC2ABD_R[NN_ABC2ABD_R] = IDR;
      NN_ABC2ABD_R++;
      GP_ABC2ABD_R[NN_ABC2ABD_R] = GP_ABC2ABD_R[NN_ABC2ABD_R-1] + Num_Rcv_Grid_ABC2ABD[IDR];
      GP_ABC2ABD_R_A[IDR] = GP_ABC2ABD_R[NN_ABC2ABD_R-1];
    }
  }

  /* set My_NumGrid_ABD */

  N3D = OFFT_Ngrid1*OFFT_Ngrid2*OFFT_Ngrid4;
  My_NumGrid_ABD = (((myid+1)*N3D+numprocs-1)/numprocs)*OFFT_Ngrid3
                  - ((myid*N3D+numprocs-1)/numprocs)*OFFT_Ngrid3;

  /****************************************************************
      ABD to DCA for MPI communication  
  ****************************************************************/

  /* set GNs (starting point in 1D ABDC) */

  N3D = OFFT_Ngrid1*OFFT_Ngrid2*OFFT_Ngrid4;
  GNs = ((myid*N3D+numprocs-1)/numprocs)*OFFT_Ngrid3;

  /* find Num_Snd_Grid_ABD2DCA */

  Num_Snd_Grid_ABD2DCA = (int*)malloc(sizeof(int)*numprocs);
  Num_Rcv_Grid_ABD2DCA = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_ABD2DCA[ID] = 0;

  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid1; /* for DCAB */

  for (BN_ABD=0; BN_ABD<My_NumGrid_ABD; BN_ABD++){

    /* position in 1D ABDC */
    GN_ABD = GNs + BN_ABD;    

    /* position in 4D ABDC n1n2n4n3 */
    n1 = GN_ABD/(OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3);
    n2 = (GN_ABD - n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3)/(OFFT_Ngrid4*OFFT_Ngrid3);
    n4 = (GN_ABD - n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3-n2*OFFT_Ngrid4*OFFT_Ngrid3)/OFFT_Ngrid3;
    n3 = GN_ABD - n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3-n2*OFFT_Ngrid4*OFFT_Ngrid3-n4*OFFT_Ngrid3;

    /* n4n3n1n2 (DCAB) */
    n3D = n4*OFFT_Ngrid3*OFFT_Ngrid1 + n3*OFFT_Ngrid1 + n1;
    ID = n3D*numprocs/N3D;
    Num_Snd_Grid_ABD2DCA[ID]++;
  }

  /* MPI: Num_Snd_Grid_ABD2DCA */  

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Snd_Grid_ABD2DCA[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Rcv_Grid_ABD2DCA[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* allocation of arrays */  

  if (OFFT_alloc_first1==0){

    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_ABD2DCA[ID]);
    }  
    free(Index_Snd_Grid_ABD2DCA);
  
    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_ABD2DCA[ID]);
    }  
    free(Index_Rcv_Grid_ABD2DCA);

    free(ID_NN_ABD2DCA_S);
    free(ID_NN_ABD2DCA_R);
    free(GP_ABD2DCA_S);
    free(GP_ABD2DCA_R);
    free(GP_ABD2DCA_S_A);
    free(GP_ABD2DCA_R_A);
  }

  size_Index_Snd_Grid_ABD2DCA = 0;
  Index_Snd_Grid_ABD2DCA = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_ABD2DCA[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid_ABD2DCA[ID]);
    size_Index_Snd_Grid_ABD2DCA += Num_Snd_Grid_ABD2DCA[ID];
  }  
  
  size_Index_Rcv_Grid_ABD2DCA = 0; 
  Index_Rcv_Grid_ABD2DCA = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_ABD2DCA[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid_ABD2DCA[ID]);
    size_Index_Rcv_Grid_ABD2DCA += Num_Rcv_Grid_ABD2DCA[ID]; 
  }  

  OFFT_alloc_first1 = 0;

  /* construct Index_Snd_Grid_ABD2DCA */  

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_ABD2DCA[ID] = 0;

  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid1; /* for DCAB */

  for (BN_ABD=0; BN_ABD<My_NumGrid_ABD; BN_ABD++){

    /* position in 1D ABDC */
    GN_ABD = GNs + BN_ABD;    

    /* position in 4D ABDC n1n2n4n3*/
    n1 = GN_ABD/(OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3);
    n2 = (GN_ABD - n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3)/(OFFT_Ngrid4*OFFT_Ngrid3);
    n4 = (GN_ABD - n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3-n2*OFFT_Ngrid4*OFFT_Ngrid3)/OFFT_Ngrid3;
    n3 = GN_ABD - n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3-n2*OFFT_Ngrid4*OFFT_Ngrid3-n4*OFFT_Ngrid3;

    /* n4n3n1n2 (DCAB) */
    n3D = n4*OFFT_Ngrid3*OFFT_Ngrid1 + n3*OFFT_Ngrid1 + n1;
    ID = n3D*numprocs/N3D;

    /* position in 1D DCAB */
    GN_DCA = n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2 + n3*OFFT_Ngrid1*OFFT_Ngrid2 + n1*OFFT_Ngrid2 + n2;
    /* offset from the starting point in 1D DCAB */
    BN_DCA = GN_DCA - ((ID*N3D+numprocs-1)/numprocs)*OFFT_Ngrid2;
    Index_Snd_Grid_ABD2DCA[ID][Num_Snd_Grid_ABD2DCA[ID]] = BN_DCA;
    Num_Snd_Grid_ABD2DCA[ID]++;
  }

  /* MPI: Index_Snd_Grid_ABD2DCA */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
   
    if (Num_Snd_Grid_ABD2DCA[IDS]!=0){
      MPI_Isend( &Index_Snd_Grid_ABD2DCA[IDS][0], Num_Snd_Grid_ABD2DCA[IDS], 
                 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Rcv_Grid_ABD2DCA[IDR]!=0){
      MPI_Recv( &Index_Rcv_Grid_ABD2DCA[IDR][0], Num_Rcv_Grid_ABD2DCA[IDR], 
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Snd_Grid_ABD2DCA[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* reset: Index_Snd_Grid_ABD2DCA

     Index_Snd_Grid_DCA2ABD:  BN_ABD
     Index_Rcv_Grid_DCA2ABD:  BN_DCA
  */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_ABD2DCA[ID] = 0;

  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid1; /* for DCAB */

  for (BN_ABD=0; BN_ABD<My_NumGrid_ABD; BN_ABD++){

    /* position in 1D ABDC */
    GN_ABD = GNs + BN_ABD;    

    /* position in 4D ABDC n1n2n4n3 */
    n1 = GN_ABD/(OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3);
    n2 = (GN_ABD - n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3)/(OFFT_Ngrid4*OFFT_Ngrid3);
    n4 = (GN_ABD - n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3-n2*OFFT_Ngrid4*OFFT_Ngrid3)/OFFT_Ngrid3;
    n3 = GN_ABD - n1*OFFT_Ngrid2*OFFT_Ngrid4*OFFT_Ngrid3-n2*OFFT_Ngrid4*OFFT_Ngrid3-n4*OFFT_Ngrid3;

    /* n4n3n1n2 (DCAB) */
    n3D = n4*OFFT_Ngrid3*OFFT_Ngrid1 + n3*OFFT_Ngrid1 + n1;
    ID = n3D*numprocs/N3D;
    Index_Snd_Grid_ABD2DCA[ID][Num_Snd_Grid_ABD2DCA[ID]] = BN_ABD;
    Num_Snd_Grid_ABD2DCA[ID]++;
  }

  /* find the maximum Num_Snd_Grid_ABD2DCA and Num_Rcv_Grid_ABD2DCA */

  Max_Num_Snd_Grid_ABD2DCA = 0;
  Max_Num_Rcv_Grid_ABD2DCA = 0;

  for (ID=0; ID<numprocs; ID++){
    if (Max_Num_Snd_Grid_ABD2DCA<Num_Snd_Grid_ABD2DCA[ID]){ 
      Max_Num_Snd_Grid_ABD2DCA = Num_Snd_Grid_ABD2DCA[ID];
    }

    if (Max_Num_Rcv_Grid_ABD2DCA<Num_Rcv_Grid_ABD2DCA[ID]){ 
      Max_Num_Rcv_Grid_ABD2DCA = Num_Rcv_Grid_ABD2DCA[ID];
    }
  }

  /* find NN_ABD2DCA_S and NN_ABD2DCA_R 
     and set ID_NN_ABD2DCA_S, 
             ID_NN_ABD2DCA_R,
             GP_ABD2DCA_S,
             GP_ABD2DCA_R
  */

  NN_ABD2DCA_S = 0;
  NN_ABD2DCA_R = 0;

  for (ID=0; ID<numprocs; ID++){
    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;
    if (Num_Snd_Grid_ABD2DCA[IDS]!=0) NN_ABD2DCA_S++;
    if (Num_Rcv_Grid_ABD2DCA[IDR]!=0) NN_ABD2DCA_R++;
  }

  ID_NN_ABD2DCA_S = (int*)malloc(sizeof(int)*NN_ABD2DCA_S);
  ID_NN_ABD2DCA_R = (int*)malloc(sizeof(int)*NN_ABD2DCA_R);
  GP_ABD2DCA_S = (int*)malloc(sizeof(int)*(NN_ABD2DCA_S+1));
  GP_ABD2DCA_R = (int*)malloc(sizeof(int)*(NN_ABD2DCA_R+1));
  GP_ABD2DCA_S_A = (int*)malloc(sizeof(int)*numprocs);
  GP_ABD2DCA_R_A = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++){
    GP_ABD2DCA_S_A[ID] = -1;
    GP_ABD2DCA_R_A[ID] = -1;
  }

  NN_ABD2DCA_S = 0;
  NN_ABD2DCA_R = 0;
  GP_ABD2DCA_S[0] = 0;
  GP_ABD2DCA_R[0] = 0;

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_ABD2DCA[IDS]!=0){
      ID_NN_ABD2DCA_S[NN_ABD2DCA_S] = IDS;
      NN_ABD2DCA_S++;
      GP_ABD2DCA_S[NN_ABD2DCA_S] = GP_ABD2DCA_S[NN_ABD2DCA_S-1] + Num_Snd_Grid_ABD2DCA[IDS];
      GP_ABD2DCA_S_A[IDS] = GP_ABD2DCA_S[NN_ABD2DCA_S-1]; 
    }

    if (Num_Rcv_Grid_ABD2DCA[IDR]!=0){
      ID_NN_ABD2DCA_R[NN_ABD2DCA_R] = IDR;
      NN_ABD2DCA_R++;
      GP_ABD2DCA_R[NN_ABD2DCA_R] = GP_ABD2DCA_R[NN_ABD2DCA_R-1] + Num_Rcv_Grid_ABD2DCA[IDR];
      GP_ABD2DCA_R_A[IDR] = GP_ABD2DCA_R[NN_ABD2DCA_R-1]; 
    }
  }

  /* set My_NumGrid_DCA */

  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid1;
  My_NumGrid_DCA = (((myid+1)*N3D+numprocs-1)/numprocs)*OFFT_Ngrid2
                  - ((myid*N3D+numprocs-1)/numprocs)*OFFT_Ngrid2;

  /****************************************************************
      DCA to DCB for MPI communication  
  ****************************************************************/

  /* set GNs (starting point in 1D DCAB) */

  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid1;
  GNs = ((myid*N3D+numprocs-1)/numprocs)*OFFT_Ngrid2;

  /* find Num_Snd_Grid_DCA2DCB */

  Num_Snd_Grid_DCA2DCB = (int*)malloc(sizeof(int)*numprocs);
  Num_Rcv_Grid_DCA2DCB = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_DCA2DCB[ID] = 0;

  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid2; /* for DCBA */

  for (BN_DCA=0; BN_DCA<My_NumGrid_DCA; BN_DCA++){

    /* position in 1D DCAB */
    GN_DCA = GNs + BN_DCA;    

    /* position in 4D DCAB n4n3n1n2 */
    n4 = GN_DCA/(OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2);
    n3 = (GN_DCA - n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2)/(OFFT_Ngrid1*OFFT_Ngrid2);
    n1 = (GN_DCA - n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2-n3*OFFT_Ngrid1*OFFT_Ngrid2)/OFFT_Ngrid2;
    n2 = GN_DCA - n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2-n3*OFFT_Ngrid1*OFFT_Ngrid2-n1*OFFT_Ngrid2;

    /* n4n3n2n1 (DCBA) */
    n3D = n4*OFFT_Ngrid3*OFFT_Ngrid2 + n3*OFFT_Ngrid2 + n2;
    ID = n3D*numprocs/N3D;
    Num_Snd_Grid_DCA2DCB[ID]++;
  }

  /* MPI: Num_Snd_Grid_DCA2DCB */  
  
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Snd_Grid_DCA2DCB[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Rcv_Grid_DCA2DCB[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* allocation of arrays */  

  if (OFFT_alloc_first2==0){

    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_DCA2DCB[ID]);
    }  
    free(Index_Snd_Grid_DCA2DCB);
  
    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_DCA2DCB[ID]);
    }  
    free(Index_Rcv_Grid_DCA2DCB);

    free(ID_NN_DCA2DCB_S);
    free(ID_NN_DCA2DCB_R);
    free(GP_DCA2DCB_S);
    free(GP_DCA2DCB_R);
    free(GP_DCA2DCB_S_A);
    free(GP_DCA2DCB_R_A);
  }
  
  size_Index_Snd_Grid_DCA2DCB = 0;
  Index_Snd_Grid_DCA2DCB = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_DCA2DCB[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid_DCA2DCB[ID]);
    size_Index_Snd_Grid_DCA2DCB += Num_Snd_Grid_DCA2DCB[ID];
  }  
  
  size_Index_Rcv_Grid_DCA2DCB = 0; 
  Index_Rcv_Grid_DCA2DCB = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_DCA2DCB[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid_DCA2DCB[ID]);
    size_Index_Rcv_Grid_DCA2DCB += Num_Rcv_Grid_DCA2DCB[ID]; 
  }  

  OFFT_alloc_first2 = 0;

  /* construct Index_Snd_Grid_DCA2DCB */  

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_DCA2DCB[ID] = 0;

  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid2; /* for DCBA */

  for (BN_DCA=0; BN_DCA<My_NumGrid_DCA; BN_DCA++){

    /* position in 1D DCAB */
    GN_DCA = GNs + BN_DCA;    

    /* position in 4D DCAB n4n3n1n2 */
    n4 = GN_DCA/(OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2);
    n3 = (GN_DCA - n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2)/(OFFT_Ngrid1*OFFT_Ngrid2);
    n1 = (GN_DCA - n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2-n3*OFFT_Ngrid1*OFFT_Ngrid2)/OFFT_Ngrid2;
    n2 = GN_DCA - n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2-n3*OFFT_Ngrid1*OFFT_Ngrid2-n1*OFFT_Ngrid2;

    /* n4n3n2n1 (DCBA) */
    n3D = n4*OFFT_Ngrid3*OFFT_Ngrid2 + n3*OFFT_Ngrid2 + n2;
    ID = n3D*numprocs/N3D;

    /* position in 1D DCBA */
    GN_DCB = n4*OFFT_Ngrid3*OFFT_Ngrid2*OFFT_Ngrid1 + n3*OFFT_Ngrid2*OFFT_Ngrid1 + n2*OFFT_Ngrid1 + n1;
    /* offset from the starting point in 1D DCBA */
    BN_DCB = GN_DCB - ((ID*N3D+numprocs-1)/numprocs)*OFFT_Ngrid1;
    Index_Snd_Grid_DCA2DCB[ID][Num_Snd_Grid_DCA2DCB[ID]] = BN_DCB;
    Num_Snd_Grid_DCA2DCB[ID]++;
  }

  /* MPI: Index_Snd_Grid_DCA2DCB */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
   
    if (Num_Snd_Grid_DCA2DCB[IDS]!=0){
      MPI_Isend( &Index_Snd_Grid_DCA2DCB[IDS][0], Num_Snd_Grid_DCA2DCB[IDS], 
                 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Rcv_Grid_DCA2DCB[IDR]!=0){
      MPI_Recv( &Index_Rcv_Grid_DCA2DCB[IDR][0], Num_Rcv_Grid_DCA2DCB[IDR], 
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Snd_Grid_DCA2DCB[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* reset: Index_Snd_Grid_DCA2DCB

     Index_Snd_Grid_DCA2DCB:  BN_ABD
     Index_Rcv_Grid_DCA2DCB:  BN_DCA
  */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_DCA2DCB[ID] = 0;

  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid2; /* for DCBA */

  for (BN_DCA=0; BN_DCA<My_NumGrid_DCA; BN_DCA++){

    /* position in 1D DCAB */
    GN_DCA = GNs + BN_DCA;    

    /* position in 4D DCAB n4n3n1n2 */
    n4 = GN_DCA/(OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2);
    n3 = (GN_DCA - n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2)/(OFFT_Ngrid1*OFFT_Ngrid2);
    n1 = (GN_DCA - n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2-n3*OFFT_Ngrid1*OFFT_Ngrid2)/OFFT_Ngrid2;
    n2 = GN_DCA - n4*OFFT_Ngrid3*OFFT_Ngrid1*OFFT_Ngrid2-n3*OFFT_Ngrid1*OFFT_Ngrid2-n1*OFFT_Ngrid2;

    /* n4n3n2n1 (DCBA) */
    n3D = n4*OFFT_Ngrid3*OFFT_Ngrid2 + n3*OFFT_Ngrid2 + n2;
    ID = n3D*numprocs/N3D;
    Index_Snd_Grid_DCA2DCB[ID][Num_Snd_Grid_DCA2DCB[ID]] = BN_DCA;
    Num_Snd_Grid_DCA2DCB[ID]++;
  }

  /* find the maximum Num_Snd_Grid_DCA2DCB and Num_Rcv_Grid_DCA2DCB */

  Max_Num_Snd_Grid_DCA2DCB = 0;
  Max_Num_Rcv_Grid_DCA2DCB = 0;

  for (ID=0; ID<numprocs; ID++){
    if (Max_Num_Snd_Grid_DCA2DCB<Num_Snd_Grid_DCA2DCB[ID]){ 
      Max_Num_Snd_Grid_DCA2DCB = Num_Snd_Grid_DCA2DCB[ID];
    }

    if (Max_Num_Rcv_Grid_DCA2DCB<Num_Rcv_Grid_DCA2DCB[ID]){ 
      Max_Num_Rcv_Grid_DCA2DCB = Num_Rcv_Grid_DCA2DCB[ID];
    }
  }

  /* find NN_DCA2DCB_S and NN_DCA2DCB_R 
     and set ID_NN_DCA2DCB_S, 
             ID_NN_DCA2DCB_R,
             GP_DCA2DCB_S,
             GP_DCA2DCB_R
  */

  NN_DCA2DCB_S = 0;
  NN_DCA2DCB_R = 0;

  for (ID=0; ID<numprocs; ID++){
    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;
    if (Num_Snd_Grid_DCA2DCB[IDS]!=0) NN_DCA2DCB_S++;
    if (Num_Rcv_Grid_DCA2DCB[IDR]!=0) NN_DCA2DCB_R++;
  }

  ID_NN_DCA2DCB_S = (int*)malloc(sizeof(int)*NN_DCA2DCB_S);
  ID_NN_DCA2DCB_R = (int*)malloc(sizeof(int)*NN_DCA2DCB_R);
  GP_DCA2DCB_S = (int*)malloc(sizeof(int)*(NN_DCA2DCB_S+1));
  GP_DCA2DCB_R = (int*)malloc(sizeof(int)*(NN_DCA2DCB_R+1));
  GP_DCA2DCB_S_A = (int*)malloc(sizeof(int)*numprocs);
  GP_DCA2DCB_R_A = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++){
    GP_DCA2DCB_S_A[ID] = -1;
    GP_DCA2DCB_R_A[ID] = -1;
  }

  NN_DCA2DCB_S = 0;
  NN_DCA2DCB_R = 0;
  GP_DCA2DCB_S[0] = 0;
  GP_DCA2DCB_R[0] = 0;

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_DCA2DCB[IDS]!=0){
      ID_NN_DCA2DCB_S[NN_DCA2DCB_S] = IDS;
      NN_DCA2DCB_S++;
      GP_DCA2DCB_S[NN_DCA2DCB_S] = GP_DCA2DCB_S[NN_DCA2DCB_S-1] + Num_Snd_Grid_DCA2DCB[IDS];
      GP_DCA2DCB_S_A[IDS] = GP_DCA2DCB_S[NN_DCA2DCB_S-1]; 
    }

    if (Num_Rcv_Grid_DCA2DCB[IDR]!=0){
      ID_NN_DCA2DCB_R[NN_DCA2DCB_R] = IDR;
      NN_DCA2DCB_R++;
      GP_DCA2DCB_R[NN_DCA2DCB_R] = GP_DCA2DCB_R[NN_DCA2DCB_R-1] + Num_Rcv_Grid_DCA2DCB[IDR];
      GP_DCA2DCB_R_A[IDR] = GP_DCA2DCB_R[NN_DCA2DCB_R-1]; 
    }
  }


  /* set My_NumGrid_DCB */

  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid2;
  My_NumGrid_DCB = (((myid+1)*N3D+numprocs-1)/numprocs)*OFFT_Ngrid1
                  - ((myid*N3D+numprocs-1)/numprocs)*OFFT_Ngrid1;
  OFFT_My_NumGrid_Out = My_NumGrid_DCB;

  /* set myStart,myEnd,OFFT_My_Index_Out */

  /* 1D */
  N3D = OFFT_Ngrid4*OFFT_Ngrid3*OFFT_Ngrid2;
  myStart = ((myid*N3D+numprocs-1)/numprocs)*OFFT_Ngrid1;
  myEnd   = (((myid+1)*N3D+numprocs-1)/numprocs)*OFFT_Ngrid1 - 1;

  /* 3D */
  N3D = OFFT_Ngrid3*OFFT_Ngrid2*OFFT_Ngrid1;
  
  My_Index_Out[0] = myStart/N3D;
  My_Index_Out[1] = (myStart-My_Index_Out[0]*N3D)/(OFFT_Ngrid2*OFFT_Ngrid1);
  My_Index_Out[2] = (myStart-My_Index_Out[0]*N3D-My_Index_Out[1]*OFFT_Ngrid2*OFFT_Ngrid1)/OFFT_Ngrid1;
  My_Index_Out[3] = myStart-My_Index_Out[0]*N3D-My_Index_Out[1]*OFFT_Ngrid2*OFFT_Ngrid1-My_Index_Out[2]*OFFT_Ngrid1;

  My_Index_Out[4] = myEnd/N3D;
  My_Index_Out[5] = (myEnd-My_Index_Out[4]*N3D)/(OFFT_Ngrid2*OFFT_Ngrid1);
  My_Index_Out[6] = (myEnd-My_Index_Out[4]*N3D-My_Index_Out[5]*OFFT_Ngrid2*OFFT_Ngrid1)/OFFT_Ngrid1;
  My_Index_Out[7] = myEnd-My_Index_Out[4]*N3D-My_Index_Out[5]*OFFT_Ngrid2*OFFT_Ngrid1-My_Index_Out[6]*OFFT_Ngrid1;

  /* find My_Max_NumGrid */

  OFFT_My_Max_NumGrid = 0;
  if (OFFT_My_Max_NumGrid<My_NumGrid_ABC) OFFT_My_Max_NumGrid = My_NumGrid_ABC;
  if (OFFT_My_Max_NumGrid<My_NumGrid_ABD) OFFT_My_Max_NumGrid = My_NumGrid_ABD;
  if (OFFT_My_Max_NumGrid<My_NumGrid_DCA) OFFT_My_Max_NumGrid = My_NumGrid_DCA;
  if (OFFT_My_Max_NumGrid<My_NumGrid_DCB) OFFT_My_Max_NumGrid = My_NumGrid_DCB;

  *My_Max_NumGrid = OFFT_My_Max_NumGrid;
  *My_NumGrid_In = OFFT_My_NumGrid_In;
  *My_NumGrid_Out = OFFT_My_NumGrid_Out;

  /* opt 1 */
  count_snd = (int*)malloc(sizeof(int)*numprocs);
  count_rcv = (int*)malloc(sizeof(int)*numprocs);
  displ_snd = (int*)malloc(sizeof(int)*numprocs);
  displ_rcv = (int*)malloc(sizeof(int)*numprocs);

  /* opt 2 */
  blocksize = BLOCK_SIZE;
  if(blocksize > numprocs) blocksize = numprocs;
  stat_arr = (MPI_Status*)malloc(sizeof(MPI_Status)*blocksize*2);
  request_arr = (MPI_Request*)malloc(sizeof(MPI_Request)*blocksize*2);

  size_array = NN_ABC2ABD_S > NN_ABD2DCA_S ? NN_ABC2ABD_S : NN_ABD2DCA_S;
  size_array = size_array > NN_DCA2DCB_S ? size_array : NN_DCA2DCB_S;
  size_array = size_array > blocksize ? size_array : blocksize;
  recvID = (int*)malloc(sizeof(int)*size_array);

  /* opt 3 */
  RQ_ABC2ABD_S = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_ABC2ABD_S);
  RQ_ABC2ABD_R = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_ABC2ABD_R);
  ST_ABC2ABD_S = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_ABC2ABD_S);
  ST_ABC2ABD_R = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_ABC2ABD_R);
  RQ_ABD2DCA_S = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_ABD2DCA_S);
  RQ_ABD2DCA_R = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_ABD2DCA_R);
  ST_ABD2DCA_S = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_ABD2DCA_S);
  ST_ABD2DCA_R = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_ABD2DCA_R);
  RQ_DCA2DCB_S = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_DCA2DCB_S);
  RQ_DCA2DCB_R = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_DCA2DCB_R);
  ST_DCA2DCB_S = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_DCA2DCB_S);
  ST_DCA2DCB_R = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_DCA2DCB_R);

  /* opt 4 */
  request_send = (MPI_Request*)malloc(sizeof(MPI_Request)*blocksize);
  request_recv = (MPI_Request*)malloc(sizeof(MPI_Request)*blocksize);
  stat_send = (MPI_Status*)malloc(sizeof(MPI_Status)*blocksize);
  stat_recv = (MPI_Status*)malloc(sizeof(MPI_Status)*blocksize);

  request_a2a = (MPI_Request*)malloc(sizeof(MPI_Request)*1);
  stat_a2a = (MPI_Status*)malloc(sizeof(MPI_Status)*numprocs);
  
  size_array0 = GP_ABC2ABD_S[NN_ABC2ABD_S] > GP_ABD2DCA_S[NN_ABD2DCA_S] ? 
    GP_ABC2ABD_S[NN_ABC2ABD_S] : GP_ABD2DCA_S[NN_ABD2DCA_S];
  size_array0 = size_array0 > GP_DCA2DCB_S[NN_DCA2DCB_S] ? 
    size_array0 : GP_DCA2DCB_S[NN_DCA2DCB_S];
  size_array1 = GP_ABC2ABD_R[NN_ABC2ABD_R] > GP_ABD2DCA_R[NN_ABD2DCA_R] ? 
    GP_ABC2ABD_R[NN_ABC2ABD_R] : GP_ABD2DCA_R[NN_ABD2DCA_R];
  size_array1 = size_array1 > GP_DCA2DCB_R[NN_DCA2DCB_R] ? 
    size_array1 : GP_DCA2DCB_R[NN_DCA2DCB_R];

  nt = fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());

  OFFT_min  = fftw_malloc(sizeof(fftw_complex)*OFFT_My_Max_NumGrid); 
  OFFT_mout = fftw_malloc(sizeof(fftw_complex)*OFFT_My_Max_NumGrid); 

  OFFT_p_mD = fftw_plan_many_dft(1,N_D,My_NumGrid_ABC/OFFT_Ngrid4,
				 OFFT_min,inembed_D,
				 1,OFFT_Ngrid4,
				 OFFT_mout,onembed_D,
				 1,OFFT_Ngrid4,
				 -1,FFTW_ESTIMATE);

  OFFT_p_mC = fftw_plan_many_dft(1,N_C,My_NumGrid_ABD/OFFT_Ngrid3,
				 OFFT_min,inembed_C,
				 1,OFFT_Ngrid3,
				 OFFT_mout,onembed_C,
				 1,OFFT_Ngrid3,
				 -1,FFTW_ESTIMATE);

  OFFT_p_mB = fftw_plan_many_dft(1,N_B,My_NumGrid_DCA/OFFT_Ngrid2,
				 OFFT_min,inembed_B,
				 1,OFFT_Ngrid2,
				 OFFT_mout,onembed_B,
				 1,OFFT_Ngrid2,
				 -1,FFTW_ESTIMATE);

  OFFT_p_mA = fftw_plan_many_dft(1,N_A,My_NumGrid_DCB/OFFT_Ngrid1,
				 OFFT_min,inembed_A,
				 1,OFFT_Ngrid1,
				 OFFT_mout,onembed_A,
				 1,OFFT_Ngrid1,
				 -1,FFTW_ESTIMATE);

  array0 = (dcomplex*)malloc(sizeof(dcomplex)*size_array0); 
  array1 = (dcomplex*)malloc(sizeof(dcomplex)*size_array1); 

  /* OFFT_MEASURE */

#ifdef kcomp

  if(OFFT_MEASURE == -1){
    COMM_PATT = COMM_DEFAULT;
  }
  else if(OFFT_MEASURE == 0){
    if(!myid) 
      printf("\nIn auto-tuning of communication pattern...\n");

    Rhor = (dcomplex*)malloc(sizeof(dcomplex)*OFFT_My_Max_NumGrid); 
    Rhok = (dcomplex*)malloc(sizeof(dcomplex)*OFFT_My_Max_NumGrid); 
   
    for (i=0; i<OFFT_My_NumGrid_In; i++){
      Rhor[i].r = 0.1; 
      Rhor[i].i = 0.2;
    }

    time_min = DBL_MAX;
    for(i=1;i<COMM_NUM;i++){
      if(i!=3){
	COMM_PATT = i;
	MPI_Barrier(mpi_comm_level1);
	time1 = 0.0;
	for(l=0;l<NumTest;l++){
	  time = MPI_Wtime();
	  openfft_exec_c2c_4d(Rhor, Rhok);
	  time1 = time1 + MPI_Wtime() - time;
	}
	MPI_Allreduce(&time1,&time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	time = time/NumTest;
	if(time_min > time){
	  time_min = time;
	  comm_patt_no = COMM_PATT;
	}
	if(!myid){
	  printf("Time of Pattern number %d = %10.15lf.\n",COMM_PATT,time);
	  printf("     Current Min Time = %10.15lf of Pattern number %d\n",
		 time_min,comm_patt_no);
	}
      }
    }
    COMM_PATT = comm_patt_no;

    free(Rhor);
    free(Rhok);
  }
  else{
    switch(OFFT_MEASURE){
    case 1: COMM_PATT = 1; break;
    case 2: COMM_PATT = 2; break;
    case 3: COMM_PATT = 3; break;
    case 4: COMM_PATT = 4; break;
    case 5: COMM_PATT = 5; break;
    case 6: COMM_PATT = 6; break;
    default: COMM_PATT = COMM_DEFAULT; break;
    }
  }
  if(!myid) 
    printf("\nExecuting with %d process(es), %d thread(s), and communication pattern number %d.\n\n",numprocs,numthreads,COMM_PATT);

#else

  if(OFFT_MEASURE == -1){
    COMM_PATT = COMM_DEFAULT;
  }
  else if(OFFT_MEASURE == 0){
    if(!myid) 
      printf("\nIn auto-tuning of communication pattern...\n");

    Rhor = (dcomplex*)malloc(sizeof(dcomplex)*OFFT_My_Max_NumGrid); 
    Rhok = (dcomplex*)malloc(sizeof(dcomplex)*OFFT_My_Max_NumGrid); 
   
    for (i=0; i<OFFT_My_NumGrid_In; i++){
      Rhor[i].r = 0.1; 
      Rhor[i].i = 0.2;
    }

    time_min = DBL_MAX;
    for(i=0;i<COMM_NUM;i++){
      COMM_PATT = i;
      MPI_Barrier(mpi_comm_level1);
      time1 = 0.0;
      for(l=0;l<NumTest;l++){
	time = MPI_Wtime();
	openfft_exec_c2c_4d(Rhor, Rhok);
	time1 = time1 + MPI_Wtime() - time;
      }
      MPI_Allreduce(&time1,&time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      time = time/NumTest;
      if(time_min > time){
	time_min = time;
	comm_patt_no = COMM_PATT;
      }
      if(!myid){
	printf("Time of Pattern number %d = %10.15lf.\n",COMM_PATT,time);
	printf("     Current Min Time = %10.15lf of Pattern number %d\n",
	       time_min,comm_patt_no);
      }
    }
    COMM_PATT = comm_patt_no;

    free(Rhor);
    free(Rhok);
  }
  else{
    switch(OFFT_MEASURE){
    case 1: COMM_PATT = 1; break;
    case 2: COMM_PATT = 2; break;
    case 3: COMM_PATT = 3; break;
    case 4: COMM_PATT = 4; break;
    case 5: COMM_PATT = 5; break;
    case 6: COMM_PATT = 6; break;
    default: COMM_PATT = COMM_DEFAULT; break;
    }
  }
  if(!myid) 
    printf("\nExecuting with %d process(es), %d thread(s), and communication pattern number %d.\n\n",numprocs,numthreads,COMM_PATT);

#endif
    
  if(OFFT_measure_time){
      MPI_Barrier(mpi_comm_level1);
      openfft_dtime(&TEtime);
      time0 = TEtime - TStime;
  }

  /* Print memory usage */

  if(OFFT_print_memory==1 && myid==0){
    MB = 1024*1024;

    printf("=========OpenFFT: Memory Usage per Process============\n");
    printf("Index_Snd_Grid_ABC2ABD: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Snd_Grid_ABC2ABD/MB));
    printf("Index_Rcv_Grid_ABC2ABD: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Rcv_Grid_ABC2ABD/MB));
    printf("Index_Snd_Grid_ABD2DCA: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Snd_Grid_ABD2DCA/MB));
    printf("Index_Rcv_Grid_ABD2DCA: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Rcv_Grid_ABD2DCA/MB));
    printf("Index_Snd_Grid_DCA2DCB: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Snd_Grid_DCA2DCB/MB));
    printf("Index_Rcv_Grid_DCA2DCB: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Rcv_Grid_DCA2DCB/MB));
    printf("ID_NN_ABC2ABD_S       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_ABC2ABD_S/MB));
    printf("ID_NN_ABC2ABD_R       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_ABC2ABD_R/MB));
    printf("ID_NN_ABD2DCA_S       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_ABD2DCA_S/MB));
    printf("ID_NN_ABD2DCA_R       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_ABD2DCA_R/MB));
    printf("ID_NN_DCA2DCB_S       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_DCA2DCB_S/MB));
    printf("ID_NN_DCA2DCB_R       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_DCA2DCB_R/MB));
    printf("GP_ABC2ABD_S          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_ABC2ABD_S+1)/MB));
    printf("GP_ABC2ABD_R          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_ABC2ABD_R+1)/MB));
    printf("GP_ABC2ABD_S_A        : %10.2f MB\n",
	   (double)(sizeof(int)*numprocs/MB));
    printf("GP_ABC2ABD_R_A        : %10.2f MB\n",
	   (double)(sizeof(int)*numprocs/MB));
    printf("GP_ABD2DCA_S          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_ABD2DCA_S+1)/MB));
    printf("GP_ABD2DCA_R          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_ABD2DCA_R+1)/MB));
    printf("GP_ABD2DCA_S_A        : %10.2f MB\n",
	   (double)(sizeof(int)*numprocs/MB));
    printf("GP_ABD2DCA_R_A        : %10.2f MB\n",
	   (double)(sizeof(int)*numprocs/MB));
    printf("GP_DCA2DCB_S          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_DCA2DCB_S+1)/MB));
    printf("GP_DCA2DCB_R          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_DCA2DCB_R+1)/MB));
    printf("GP_DCA2DCB_S_A        : %10.2f MB\n",
	   (double)(sizeof(int)*numprocs/MB));
    printf("GP_DCA2DCB_R_A        : %10.2f MB\n",
	   (double)(sizeof(int)*numprocs/MB));
    printf("Temporary array0    : %10.2f MB\n",
	   (double)(sizeof(dcomplex)*size_array0/MB));
    printf("Temporary array1    : %10.2f MB\n",
	   (double)(sizeof(dcomplex)*size_array1/MB));
    printf("Temporary OFFT_in   : %10.2f MB\n",
	   (double)(sizeof(fftw_complex)*Max_OneD_Grids/MB));
    printf("Temporary OFFT_out  : %10.2f MB\n",
	   (double)(sizeof(fftw_complex)*Max_OneD_Grids/MB));
    printf("=====================================================\n");
  }

  return time0;
}  


