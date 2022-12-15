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
  openfft_init_c2c_3d prepares data for the c2c transformation of 3D data.

1. Input: 

 - 3 dimensions of data: N1, N2, N3.

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

double openfft_init_c2c_3d(int N1, int N2, int N3,
			   int *My_Max_NumGrid, 
			   int *My_NumGrid_In, int *My_Index_In,
			   int *My_NumGrid_Out, int *My_Index_Out,
			   int offt_measure,
			   int measure_time, int print_memory)
{

  int numprocs,myid,tag=999,ID,IDS,IDR;
  long long int BN_AB,BN_CB,BN_CA;
  long long int GNs,n2D,N2D;
  long long int GN_AB,GN_CB,GN_CA;
  long long int n1,n2,n3,blocksize,i,comm_patt_no,l;
  MPI_Status stat;
  MPI_Request request;
  double time0,TStime,TEtime;
  long long int myStart,myEnd;
  double MB,time,time1,time_min;
  long long int size_array0,size_array1,size_array;
  MPI_Comm mpi_comm_level1 = MPI_COMM_WORLD;
  dcomplex *Rhor, *Rhok;
  int numthreads;
  int N_A[] = {N1}, N_B[] = {N2}, N_C[] = {N3};
  int *inembed_A = N_A, *onembed_A = N_A;
  int *inembed_B = N_B, *onembed_B = N_B;
  int *inembed_C = N_C, *onembed_C = N_C;
  
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
    printf("\n=================OpenFFT1.2 c2c 3-D Interface===============");
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

  OFFT_Ngrid1 = N1;
  OFFT_Ngrid2 = N2;
  OFFT_Ngrid3 = N3;

  Max_OneD_Grids = OFFT_Ngrid1 > OFFT_Ngrid2 ? OFFT_Ngrid1 : OFFT_Ngrid2;
  Max_OneD_Grids =  Max_OneD_Grids > OFFT_Ngrid3 ?  Max_OneD_Grids : OFFT_Ngrid3;

  /****************************************************************
      AB to CA for MPI communication  
  ****************************************************************/

  /* set My_NumGrid_AB */

  N2D = OFFT_Ngrid1*OFFT_Ngrid2;
  My_NumGrid_AB =  (((myid+1)*N2D+numprocs-1)/numprocs)*OFFT_Ngrid3
                  - ((myid*N2D+numprocs-1)/numprocs)*OFFT_Ngrid3;
  OFFT_My_NumGrid_In = My_NumGrid_AB; 

  /* set myStart,myEnd,OFFT_My_Index_In */

  N2D = OFFT_Ngrid1*OFFT_Ngrid2;
  myStart = ((myid*N2D+numprocs-1)/numprocs)*OFFT_Ngrid3;
  myEnd   = (((myid+1)*N2D+numprocs-1)/numprocs)*OFFT_Ngrid3 - 1;

  N2D = OFFT_Ngrid2*OFFT_Ngrid3;
  My_Index_In[0] = myStart/N2D;
  My_Index_In[1] = (myStart-My_Index_In[0]*N2D)/OFFT_Ngrid3;
  My_Index_In[2] = myStart-My_Index_In[0]*N2D-My_Index_In[1]*OFFT_Ngrid3;
  
  My_Index_In[3] = myEnd/N2D;
  My_Index_In[4] = (myEnd-My_Index_In[3]*N2D)/OFFT_Ngrid3;
  My_Index_In[5] = myEnd-My_Index_In[3]*N2D-My_Index_In[4]*OFFT_Ngrid3;

  /* set GNs */

  N2D = OFFT_Ngrid1*OFFT_Ngrid2;
  GNs = ((myid*N2D+numprocs-1)/numprocs)*OFFT_Ngrid3;

  /* find Num_Snd_Grid_AB2CA */

  Num_Snd_Grid_AB2CA = (int*)malloc(sizeof(int)*numprocs);
  Num_Rcv_Grid_AB2CA = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_AB2CA[ID] = 0;

  N2D = OFFT_Ngrid3*OFFT_Ngrid1; /* for CA */

  for (BN_AB=0; BN_AB<My_NumGrid_AB; BN_AB++){

    GN_AB = BN_AB + GNs;

    n1 = GN_AB/(OFFT_Ngrid2*OFFT_Ngrid3);
    n2 = (GN_AB - n1*(OFFT_Ngrid2*OFFT_Ngrid3))/OFFT_Ngrid3;
    n3 = GN_AB - n1*(OFFT_Ngrid2*OFFT_Ngrid3) - n2*OFFT_Ngrid3;
    n2D = n3*OFFT_Ngrid1 + n1;
    ID = n2D*numprocs/N2D;
    Num_Snd_Grid_AB2CA[ID]++;
  }

  /* MPI: Num_Snd_Grid_AB2CA */  

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Snd_Grid_AB2CA[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Rcv_Grid_AB2CA[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* allocation of arrays */  

  if (OFFT_alloc_first==0){

    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_AB2CA[ID]);
    }  
    free(Index_Snd_Grid_AB2CA);
  
    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_AB2CA[ID]);
    }  
    free(Index_Rcv_Grid_AB2CA);

    free(ID_NN_AB2CA_S);
    free(ID_NN_AB2CA_R);
    free(GP_AB2CA_S);
    free(GP_AB2CA_R);
    free(GP_AB2CA_S_A);
    free(GP_AB2CA_R_A);
  }

  size_Index_Snd_Grid_AB2CA = 0;
  Index_Snd_Grid_AB2CA = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_AB2CA[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid_AB2CA[ID]);
    size_Index_Snd_Grid_AB2CA += Num_Snd_Grid_AB2CA[ID];
  }  
  
  size_Index_Rcv_Grid_AB2CA = 0; 
  Index_Rcv_Grid_AB2CA = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_AB2CA[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid_AB2CA[ID]);
    size_Index_Rcv_Grid_AB2CA += Num_Rcv_Grid_AB2CA[ID]; 
  }  

  OFFT_alloc_first = 0;

  /* construct Index_Snd_Grid_AB2CA */  

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_AB2CA[ID] = 0;

  N2D = OFFT_Ngrid3*OFFT_Ngrid1;  /* for CA */

  for (BN_AB=0; BN_AB<My_NumGrid_AB; BN_AB++){

    GN_AB = BN_AB + GNs;

    n1 = GN_AB/(OFFT_Ngrid2*OFFT_Ngrid3);
    n2 = (GN_AB - n1*(OFFT_Ngrid2*OFFT_Ngrid3))/OFFT_Ngrid3;
    n3 = GN_AB - n1*(OFFT_Ngrid2*OFFT_Ngrid3) - n2*OFFT_Ngrid3;
    n2D = n3*OFFT_Ngrid1 + n1;
    ID = n2D*numprocs/N2D;
    GN_CA = n3*OFFT_Ngrid1*OFFT_Ngrid2 + n1*OFFT_Ngrid2 + n2;
    BN_CA = GN_CA - ((ID*N2D+numprocs-1)/numprocs)*OFFT_Ngrid2;
    Index_Snd_Grid_AB2CA[ID][Num_Snd_Grid_AB2CA[ID]] = BN_CA;
    Num_Snd_Grid_AB2CA[ID]++;
  }

  /* MPI: Index_Snd_Grid_AB2CA */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
   
    if (Num_Snd_Grid_AB2CA[IDS]!=0){
      MPI_Isend( &Index_Snd_Grid_AB2CA[IDS][0], Num_Snd_Grid_AB2CA[IDS], 
		 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Rcv_Grid_AB2CA[IDR]!=0){
      MPI_Recv( &Index_Rcv_Grid_AB2CA[IDR][0], Num_Rcv_Grid_AB2CA[IDR], 
		MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Snd_Grid_AB2CA[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* reset: Index_Snd_Grid_AB2CA

  Index_Snd_Grid_AB2CA:  BN_AB
  Index_Rcv_Grid_AB2CA:  BN_CA
  */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_AB2CA[ID] = 0;

  N2D = OFFT_Ngrid3*OFFT_Ngrid1; /* for CA */

  for (BN_AB=0; BN_AB<My_NumGrid_AB; BN_AB++){

    GN_AB = BN_AB + GNs;
    n1 = GN_AB/(OFFT_Ngrid2*OFFT_Ngrid3);
    n2 = (GN_AB - n1*(OFFT_Ngrid2*OFFT_Ngrid3))/OFFT_Ngrid3;
    n3 = GN_AB - n1*(OFFT_Ngrid2*OFFT_Ngrid3) - n2*OFFT_Ngrid3;
    n2D = n3*OFFT_Ngrid1 + n1;
    ID = n2D*numprocs/N2D;
    Index_Snd_Grid_AB2CA[ID][Num_Snd_Grid_AB2CA[ID]] = BN_AB;
    Num_Snd_Grid_AB2CA[ID]++;
  }

  /* find the maximum Num_Snd_Grid_AB2CA and Num_Rcv_Grid_AB2CA */

  Max_Num_Snd_Grid_AB2CA = 0;
  Max_Num_Rcv_Grid_AB2CA = 0;

  for (ID=0; ID<numprocs; ID++){
    if (Max_Num_Snd_Grid_AB2CA<Num_Snd_Grid_AB2CA[ID]){ 
      Max_Num_Snd_Grid_AB2CA = Num_Snd_Grid_AB2CA[ID];
    }

    if (Max_Num_Rcv_Grid_AB2CA<Num_Rcv_Grid_AB2CA[ID]){ 
      Max_Num_Rcv_Grid_AB2CA = Num_Rcv_Grid_AB2CA[ID];
    }
  }

  /* find NN_AB2CA_S and NN_AB2CA_R 
     and set ID_NN_AB2CA_S, 
             ID_NN_AB2CA_R,
             GP_AB2CA_S,
             GP_AB2CA_R
  */

  NN_AB2CA_S = 0;
  NN_AB2CA_R = 0;
 
  for (ID=0; ID<numprocs; ID++){
    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;
    if (Num_Snd_Grid_AB2CA[IDS]!=0) NN_AB2CA_S++;
    if (Num_Rcv_Grid_AB2CA[IDR]!=0) NN_AB2CA_R++;
  }

  ID_NN_AB2CA_S = (int*)malloc(sizeof(int)*NN_AB2CA_S);
  ID_NN_AB2CA_R = (int*)malloc(sizeof(int)*NN_AB2CA_R);
  GP_AB2CA_S = (int*)malloc(sizeof(int)*(NN_AB2CA_S+1));
  GP_AB2CA_R = (int*)malloc(sizeof(int)*(NN_AB2CA_R+1));
  GP_AB2CA_S_A = (int*)malloc(sizeof(int)*numprocs);
  GP_AB2CA_R_A = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++){
    GP_AB2CA_S_A[ID] = -1;
    GP_AB2CA_R_A[ID] = -1;
  }

  NN_AB2CA_S = 0;
  NN_AB2CA_R = 0;
  GP_AB2CA_S[0] = 0;
  GP_AB2CA_R[0] = 0;

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_AB2CA[IDS]!=0){
      ID_NN_AB2CA_S[NN_AB2CA_S] = IDS;
      NN_AB2CA_S++;
      GP_AB2CA_S[NN_AB2CA_S] = GP_AB2CA_S[NN_AB2CA_S-1] + Num_Snd_Grid_AB2CA[IDS];
      GP_AB2CA_S_A[IDS] = GP_AB2CA_S[NN_AB2CA_S-1];
    }

    if (Num_Rcv_Grid_AB2CA[IDR]!=0){
      ID_NN_AB2CA_R[NN_AB2CA_R] = IDR;
      NN_AB2CA_R++;
      GP_AB2CA_R[NN_AB2CA_R] = GP_AB2CA_R[NN_AB2CA_R-1] + Num_Rcv_Grid_AB2CA[IDR];
      GP_AB2CA_R_A[IDR] = GP_AB2CA_R[NN_AB2CA_R-1];
    }
  }

  /* set My_NumGrid_CA */

  N2D = OFFT_Ngrid3*OFFT_Ngrid1;
  My_NumGrid_CA = (((myid+1)*N2D+numprocs-1)/numprocs)*OFFT_Ngrid2
                  - ((myid*N2D+numprocs-1)/numprocs)*OFFT_Ngrid2;

  /****************************************************************
      CA to CB for MPI communication  
  ****************************************************************/

  /* set GNs */

  N2D = OFFT_Ngrid3*OFFT_Ngrid1;
  GNs = ((myid*N2D+numprocs-1)/numprocs)*OFFT_Ngrid2;

  /* find Num_Snd_Grid_CA2CB */

  Num_Snd_Grid_CA2CB = (int*)malloc(sizeof(int)*numprocs);
  Num_Rcv_Grid_CA2CB = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_CA2CB[ID] = 0;

  N2D = OFFT_Ngrid3*OFFT_Ngrid2; /* for CB */

  for (BN_CA=0; BN_CA<My_NumGrid_CA; BN_CA++){

    GN_CA = GNs + BN_CA;    
    n3 = GN_CA/(OFFT_Ngrid1*OFFT_Ngrid2);
    n1 = (GN_CA - n3*(OFFT_Ngrid1*OFFT_Ngrid2))/OFFT_Ngrid2;
    n2 = GN_CA - n3*(OFFT_Ngrid1*OFFT_Ngrid2) - n1*OFFT_Ngrid2;
    n2D = n3*OFFT_Ngrid2 + n2;
    ID = n2D*numprocs/N2D;
    Num_Snd_Grid_CA2CB[ID]++;
  }

  /* MPI: Num_Snd_Grid_CA2CB */  

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Snd_Grid_CA2CB[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Rcv_Grid_CA2CB[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* allocation of arrays */  

  if (OFFT_alloc_first1==0){

    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_CA2CB[ID]);
    }  
    free(Index_Snd_Grid_CA2CB);
  
    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_CA2CB[ID]);
    }  
    free(Index_Rcv_Grid_CA2CB);

    free(ID_NN_CA2CB_S);
    free(ID_NN_CA2CB_R);
    free(GP_CA2CB_S);
    free(GP_CA2CB_R);
    free(GP_CA2CB_S_A);
    free(GP_CA2CB_R_A);
  }

  size_Index_Snd_Grid_CA2CB = 0;
  Index_Snd_Grid_CA2CB = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_CA2CB[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid_CA2CB[ID]);
    size_Index_Snd_Grid_CA2CB += Num_Snd_Grid_CA2CB[ID];
  }  
  
  size_Index_Rcv_Grid_CA2CB = 0; 
  Index_Rcv_Grid_CA2CB = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_CA2CB[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid_CA2CB[ID]);
    size_Index_Rcv_Grid_CA2CB += Num_Rcv_Grid_CA2CB[ID]; 
  }  

  OFFT_alloc_first1 = 0;

  /* construct Index_Snd_Grid_CA2CB */  

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_CA2CB[ID] = 0;

  N2D = OFFT_Ngrid3*OFFT_Ngrid2; /* for CB */

  for (BN_CA=0; BN_CA<My_NumGrid_CA; BN_CA++){

    GN_CA = GNs + BN_CA;    
    n3 = GN_CA/(OFFT_Ngrid1*OFFT_Ngrid2);
    n1 = (GN_CA - n3*(OFFT_Ngrid1*OFFT_Ngrid2))/OFFT_Ngrid2;
    n2 = GN_CA - n3*(OFFT_Ngrid1*OFFT_Ngrid2) - n1*OFFT_Ngrid2;
    n2D = n3*OFFT_Ngrid2 + n2;
    ID = n2D*numprocs/N2D;
    GN_CB = n3*OFFT_Ngrid2*OFFT_Ngrid1 + n2*OFFT_Ngrid1 + n1;
    BN_CB = GN_CB - ((ID*N2D+numprocs-1)/numprocs)*OFFT_Ngrid1;
    Index_Snd_Grid_CA2CB[ID][Num_Snd_Grid_CA2CB[ID]] = BN_CB;
    Num_Snd_Grid_CA2CB[ID]++;
  }

  /* MPI: Index_Snd_Grid_CA2CB */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
   
    if (Num_Snd_Grid_CA2CB[IDS]!=0){
      MPI_Isend( &Index_Snd_Grid_CA2CB[IDS][0], Num_Snd_Grid_CA2CB[IDS], 
                 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Rcv_Grid_CA2CB[IDR]!=0){
      MPI_Recv( &Index_Rcv_Grid_CA2CB[IDR][0], Num_Rcv_Grid_CA2CB[IDR], 
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Snd_Grid_CA2CB[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* reset: Index_Snd_Grid_CA2CB

     Index_Snd_Grid_CB2CA:  BN_CA
     Index_Rcv_Grid_CB2CA:  BN_CB
  */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_CA2CB[ID] = 0;

  N2D = OFFT_Ngrid3*OFFT_Ngrid2; /* for CB */

  for (BN_CA=0; BN_CA<My_NumGrid_CA; BN_CA++){

    GN_CA = GNs + BN_CA;    
    n3 = GN_CA/(OFFT_Ngrid1*OFFT_Ngrid2);
    n1 = (GN_CA - n3*(OFFT_Ngrid1*OFFT_Ngrid2))/OFFT_Ngrid2;
    n2 = GN_CA - n3*(OFFT_Ngrid1*OFFT_Ngrid2) - n1*OFFT_Ngrid2;
    n2D = n3*OFFT_Ngrid2 + n2;
    ID = n2D*numprocs/N2D;
    Index_Snd_Grid_CA2CB[ID][Num_Snd_Grid_CA2CB[ID]] = BN_CA;
    Num_Snd_Grid_CA2CB[ID]++;
  }

  /* find the maximum Num_Snd_Grid_CA2CB and Num_Rcv_Grid_CA2CB */

  Max_Num_Snd_Grid_CA2CB = 0;
  Max_Num_Rcv_Grid_CA2CB = 0;

  for (ID=0; ID<numprocs; ID++){
    if (Max_Num_Snd_Grid_CA2CB<Num_Snd_Grid_CA2CB[ID]){ 
      Max_Num_Snd_Grid_CA2CB = Num_Snd_Grid_CA2CB[ID];
    }

    if (Max_Num_Rcv_Grid_CA2CB<Num_Rcv_Grid_CA2CB[ID]){ 
      Max_Num_Rcv_Grid_CA2CB = Num_Rcv_Grid_CA2CB[ID];
    }
  }

  /* find NN_CA2CB_S and NN_CA2CB_R 
     and set ID_NN_CA2CB_S, 
             ID_NN_CA2CB_R,
             GP_CA2CB_S,
             GP_CA2CB_R
  */

  NN_CA2CB_S = 0;
  NN_CA2CB_R = 0;

  for (ID=0; ID<numprocs; ID++){
    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;
    if (Num_Snd_Grid_CA2CB[IDS]!=0) NN_CA2CB_S++;
    if (Num_Rcv_Grid_CA2CB[IDR]!=0) NN_CA2CB_R++;
  }

  ID_NN_CA2CB_S = (int*)malloc(sizeof(int)*NN_CA2CB_S);
  ID_NN_CA2CB_R = (int*)malloc(sizeof(int)*NN_CA2CB_R);
  GP_CA2CB_S = (int*)malloc(sizeof(int)*(NN_CA2CB_S+1));
  GP_CA2CB_R = (int*)malloc(sizeof(int)*(NN_CA2CB_R+1));
  GP_CA2CB_S_A = (int*)malloc(sizeof(int)*numprocs);
  GP_CA2CB_R_A = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++){
    GP_CA2CB_S_A[ID] = -1;
    GP_CA2CB_R_A[ID] = -1;
  }

  NN_CA2CB_S = 0;
  NN_CA2CB_R = 0;
  GP_CA2CB_S[0] = 0;
  GP_CA2CB_R[0] = 0;

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_CA2CB[IDS]!=0){
      ID_NN_CA2CB_S[NN_CA2CB_S] = IDS;
      NN_CA2CB_S++;
      GP_CA2CB_S[NN_CA2CB_S] = GP_CA2CB_S[NN_CA2CB_S-1] + Num_Snd_Grid_CA2CB[IDS];
      GP_CA2CB_S_A[IDS] = GP_CA2CB_S[NN_CA2CB_S-1]; 
    }

    if (Num_Rcv_Grid_CA2CB[IDR]!=0){
      ID_NN_CA2CB_R[NN_CA2CB_R] = IDR;
      NN_CA2CB_R++;
      GP_CA2CB_R[NN_CA2CB_R] = GP_CA2CB_R[NN_CA2CB_R-1] + Num_Rcv_Grid_CA2CB[IDR];
      GP_CA2CB_R_A[IDR] = GP_CA2CB_R[NN_CA2CB_R-1]; 
    }
  }

  /* set My_NumGrid_CB */

  N2D = OFFT_Ngrid3*OFFT_Ngrid2;
  My_NumGrid_CB = (((myid+1)*N2D+numprocs-1)/numprocs)*OFFT_Ngrid1
                  - ((myid*N2D+numprocs-1)/numprocs)*OFFT_Ngrid1;
  OFFT_My_NumGrid_Out = My_NumGrid_CB;

  /* set myStart,myEnd,OFFT_My_Index_Out */

  N2D = OFFT_Ngrid3*OFFT_Ngrid2;
  myStart = ((myid*N2D+numprocs-1)/numprocs)*OFFT_Ngrid1;
  myEnd   = (((myid+1)*N2D+numprocs-1)/numprocs)*OFFT_Ngrid1 - 1;

  N2D = OFFT_Ngrid2*OFFT_Ngrid1;
  My_Index_Out[0] = myStart/N2D;
  My_Index_Out[1] = (myStart-My_Index_Out[0]*N2D)/OFFT_Ngrid1;
  My_Index_Out[2] = myStart-My_Index_Out[0]*N2D-My_Index_Out[1]*OFFT_Ngrid1;
  
  My_Index_Out[3] = myEnd/N2D;
  My_Index_Out[4] = (myEnd-My_Index_Out[3]*N2D)/OFFT_Ngrid1;
  My_Index_Out[5] = myEnd-My_Index_Out[3]*N2D-My_Index_Out[4]*OFFT_Ngrid1;

  /* find My_Max_NumGrid */

  OFFT_My_Max_NumGrid = 0;
  if (OFFT_My_Max_NumGrid<My_NumGrid_AB) OFFT_My_Max_NumGrid = My_NumGrid_AB;
  if (OFFT_My_Max_NumGrid<My_NumGrid_CA) OFFT_My_Max_NumGrid = My_NumGrid_CA;
  if (OFFT_My_Max_NumGrid<My_NumGrid_CB) OFFT_My_Max_NumGrid = My_NumGrid_CB;

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

  size_array = NN_AB2CA_S > NN_CA2CB_S ? NN_AB2CA_S : NN_CA2CB_S;
  size_array = size_array > blocksize ? size_array : blocksize;
  recvID = (int*)malloc(sizeof(int)*size_array);

  /* opt 3 */
  RQ_AB2CA_S = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_AB2CA_S);
  RQ_AB2CA_R = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_AB2CA_R);
  ST_AB2CA_S = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_AB2CA_S);
  ST_AB2CA_R = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_AB2CA_R);
  RQ_CA2CB_S = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_CA2CB_S);
  RQ_CA2CB_R = (MPI_Request*)malloc(sizeof(MPI_Request)*NN_CA2CB_R);
  ST_CA2CB_S = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_CA2CB_S);
  ST_CA2CB_R = (MPI_Status*)malloc(sizeof(MPI_Status)*NN_CA2CB_R);

  /* opt 4 */
  request_send = (MPI_Request*)malloc(sizeof(MPI_Request)*blocksize);
  request_recv = (MPI_Request*)malloc(sizeof(MPI_Request)*blocksize);
  stat_send = (MPI_Status*)malloc(sizeof(MPI_Status)*blocksize);
  stat_recv = (MPI_Status*)malloc(sizeof(MPI_Status)*blocksize);

  size_array0 = GP_AB2CA_S[NN_AB2CA_S] > GP_CA2CB_S[NN_CA2CB_S] ? 
    GP_AB2CA_S[NN_AB2CA_S] : GP_CA2CB_S[NN_CA2CB_S];
  size_array1 = GP_AB2CA_R[NN_AB2CA_R] > GP_CA2CB_R[NN_CA2CB_R] ? 
    GP_AB2CA_R[NN_AB2CA_R] : GP_CA2CB_R[NN_CA2CB_R];
  array0 = (dcomplex*)malloc(sizeof(dcomplex)*size_array0); 
  array1 = (dcomplex*)malloc(sizeof(dcomplex)*size_array1); 

  numthreads = omp_get_max_threads();
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());

  OFFT_min  = fftw_malloc(sizeof(fftw_complex)*OFFT_My_Max_NumGrid); 
  OFFT_mout = fftw_malloc(sizeof(fftw_complex)*OFFT_My_Max_NumGrid); 

  OFFT_p_mC = fftw_plan_many_dft(1,N_C,My_NumGrid_AB/OFFT_Ngrid3,
				 OFFT_min,inembed_C,
				 1,OFFT_Ngrid3,
				 OFFT_mout,onembed_C,
				 1,OFFT_Ngrid3,
				 -1,FFTW_ESTIMATE);

  OFFT_p_mB = fftw_plan_many_dft(1,N_B,My_NumGrid_CA/OFFT_Ngrid2,
				 OFFT_min,inembed_B,
				 1,OFFT_Ngrid2,
				 OFFT_mout,onembed_B,
				 1,OFFT_Ngrid2,
				 -1,FFTW_ESTIMATE);

  OFFT_p_mA = fftw_plan_many_dft(1,N_A,My_NumGrid_CB/OFFT_Ngrid1,
				 OFFT_min,inembed_A,
				 1,OFFT_Ngrid1,
				 OFFT_mout,onembed_A,
				 1,OFFT_Ngrid1,
				 -1,FFTW_ESTIMATE);

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
	  openfft_exec_c2c_3d(Rhor, Rhok);
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
	openfft_exec_c2c_3d(Rhor, Rhok);
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
    printf("Index_Snd_Grid_AB2CA: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Snd_Grid_AB2CA/MB));
    printf("Index_Rcv_Grid_AB2CA: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Rcv_Grid_AB2CA/MB));
    printf("Index_Snd_Grid_CA2CB: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Snd_Grid_CA2CB/MB));
    printf("Index_Rcv_Grid_CA2CB: %10.2f MB\n",
	   (double)(sizeof(int)*size_Index_Rcv_Grid_CA2CB/MB));
    printf("ID_NN_AB2CA_S       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_AB2CA_S/MB));
    printf("ID_NN_AB2CA_R       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_AB2CA_R/MB));
    printf("ID_NN_CA2CB_S       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_CA2CB_S/MB));
    printf("ID_NN_CA2CB_R       : %10.2f MB\n",
	   (double)(sizeof(int)*NN_CA2CB_R/MB));
    printf("GP_AB2CA_S          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_AB2CA_S+1)/MB));
    printf("GP_AB2CA_R          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_AB2CA_R+1)/MB));
    printf("GP_AB2CA_S_A        : %10.2f MB\n",
	   (double)(sizeof(int)*numprocs/MB));
    printf("GP_AB2CA_R_A        : %10.2f MB\n",
	   (double)(sizeof(int)*numprocs/MB));
    printf("GP_CA2CB_S          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_CA2CB_S+1)/MB));
    printf("GP_CA2CB_R          : %10.2f MB\n",
	   (double)(sizeof(int)*(NN_CA2CB_R+1)/MB));
    printf("GP_CA2CB_S_A        : %10.2f MB\n",
	   (double)(sizeof(int)*numprocs/MB));
    printf("GP_CA2CB_R_A        : %10.2f MB\n",
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


