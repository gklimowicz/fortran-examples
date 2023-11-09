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
    MERCHANTABCILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

******************************************************/

/**********************************************************************
  openfft_exec_c2c_4d performs the c2c transformation of the 4D data. 
  Must be called after openfft_init_c2c_4d.

     Input : Rhor (double complex).

     Output: Rhok (double complex).

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "openfft.h"

double openfft_exec_c2c_4d(dcomplex *Rhor, dcomplex *Rhok) 
{
  long long int i,j,BN_ABC,BN_DCA,BN_ABD,BN_DCB,gp,NN_S,NN_R,NN_RS;
  dcomplex *sendbuf,*recvbuf;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  int blocksize,n,bsize,b,recvCount,recvNum,sendcnt,recvcnt;
  double Stime_proc,Etime_proc;
  double time0,time1,time00,time11,timeC,timeF,timeDf,timeDc,timeD,timeT;
  MPI_Status stat;
  MPI_Comm mpi_comm_level1 = MPI_COMM_WORLD;
  int numthreads,mytid;

  time0 = 0.0;
  timeC = 0.0;
  timeF = 0.0;
  timeD = 0.0;
  timeDf = 0.0;
  timeDc = 0.0;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  

  /*------------------ FFT along the D-axis in the ABC partition ------------------*/

  if (OFFT_measure_time==1){
    openfft_dtime(&Stime_proc);
    openfft_dtime(&time0);
  }

  fftw_execute_dft(OFFT_p_mD,(fftw_complex*)Rhor,(fftw_complex*)Rhok);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeF += time1 - time0;
  }

  /*
  if(!myid){
    printf("\nABC\n");
    for(i=0;i<My_NumGrid_ABC;i++){
      printf("%d %f\n",i,Rhor[i]);
    }
  }
  /**/

  /*------------------ MPI: ABC to ABD partitions ------------------*/

switch (COMM_PATT){

case 1: {

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }

  for(i=0;i<numprocs;i+=1){
    count_snd[i] = 0;
    count_rcv[i] = 0;
    displ_snd[i] = 0;
    displ_rcv[i] = 0;
  }

  for (ID=0; ID<NN_ABC2ABD_S; ID+=1){

    IDS = ID_NN_ABC2ABD_S[ID];
    gp = GP_ABC2ABD_S[ID];

    count_snd[IDS] = Num_Snd_Grid_ABC2ABD[IDS]*2;
    displ_snd[IDS] = gp*2; 

    if(IDS!=myid){
      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDS]; i++){
	BN_ABC = Index_Snd_Grid_ABC2ABD[IDS][i];
	array0[gp+i].r = Rhok[BN_ABC].r;
	array0[gp+i].i = Rhok[BN_ABC].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      } 
    }
  }

  for (ID=0; ID<NN_ABC2ABD_R; ID+=1){

    IDR = ID_NN_ABC2ABD_R[ID];
    gp = GP_ABC2ABD_R[ID];

    count_rcv[IDR] = Num_Rcv_Grid_ABC2ABD[IDR]*2;
    displ_rcv[IDR] = gp*2; 
  }

  count_snd[myid] = 0;
  count_rcv[myid] = 0;

  if (OFFT_measure_time==1) openfft_dtime(&time00);

  MPI_Alltoallv(array0, count_snd, displ_snd, MPI_DOUBLE,
		array1, count_rcv, displ_rcv, MPI_DOUBLE, mpi_comm_level1);

  if (OFFT_measure_time==1){
    openfft_dtime(&time11);
    timeT += time11 - time00;
    timeC += time11 - time00;
  }

  /* copy them to Rhor */

  for (ID=0; ID<NN_ABC2ABD_R; ID+=1){

    IDR = ID_NN_ABC2ABD_R[ID];
    gp = GP_ABC2ABD_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABC2ABD[IDR]; i++){
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = array1[gp+i].r;
	Rhor[BN_ABD].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_ABD],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    else{
      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDR]; i++){
	BN_ABC = Index_Snd_Grid_ABC2ABD[IDR][i];
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = Rhok[BN_ABC].r;
	Rhor[BN_ABD].i = Rhok[BN_ABC].i;
	/*
	memcpy(&Rhor[BN_ABD],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      }
    }
  }

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += time1 - time0 - timeT;
  }

break;

} /* opt == 1 */

case 2: {

  blocksize = BLOCK_SIZE;
  if(blocksize > numprocs) blocksize = numprocs;
  
  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for(n=0; n<numprocs; n+=blocksize){
   
  NN_RS = 0;
  bsize = numprocs-n < blocksize ? numprocs-n : blocksize;

  /* MPI_Irecv */

  for (b=0; b<bsize; b++){

    IDR = (myid+b+n) % numprocs;
    gp = GP_ABC2ABD_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_ABC2ABD[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_arr[NN_RS]);
      NN_RS++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_ABC2ABD_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDS]; i++){
        BN_ABC = Index_Snd_Grid_ABC2ABD[IDS][i];
        array0[gp+i].r = Rhok[BN_ABC].r;
        array0[gp+i].i = Rhok[BN_ABC].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_ABC2ABD[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_arr[NN_RS]);
      NN_RS++;  
    }
  }
  
  /* MPI_Waitall */

  MPI_Waitall(NN_RS,request_arr,stat_arr);

  } /* n */

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

  /* copy them to Rhor */

  for (ID=0; ID<NN_ABC2ABD_R; ID++){

    IDR = ID_NN_ABC2ABD_R[ID];
    gp = GP_ABC2ABD_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABC2ABD[IDR]; i++){
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = array1[gp+i].r;
	Rhor[BN_ABD].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_ABD],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDR]; i++){
	BN_ABC = Index_Snd_Grid_ABC2ABD[IDR][i];
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = Rhok[BN_ABC].r;
	Rhor[BN_ABD].i = Rhok[BN_ABC].i;
	/*
	memcpy(&Rhor[BN_ABD],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt = 2*/

case 3: {

  NN_S = 0;
  NN_R = 0;
  recvCount = 0;

  if (OFFT_measure_time==1) openfft_dtime(&time0);
   
  /* MPI_Irecv */

  for (ID=0; ID<NN_ABC2ABD_R; ID++){

    IDR = ID_NN_ABC2ABD_R[ID];
    gp = GP_ABC2ABD_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_ABC2ABD[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_ABC2ABD_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_ABC2ABD_S; ID++){

    IDS = ID_NN_ABC2ABD_S[ID];
    gp = GP_ABC2ABD_S[ID];

    if (IDS!=myid){
      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDS]; i++){
        BN_ABC = Index_Snd_Grid_ABC2ABD[IDS][i];
        array0[gp+i].r = Rhok[BN_ABC].r;
        array0[gp+i].i = Rhok[BN_ABC].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      }
      MPI_Isend( &array0[gp], Num_Snd_Grid_ABC2ABD[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_ABC2ABD_S[NN_S]);
      NN_S++;  
    }
  }

  /* MPI_Waitsome */ 

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,RQ_ABC2ABD_R,&recvNum,recvID,ST_ABC2ABD_R);
    for(j=0;j<recvNum;j++){
      IDR = ST_ABC2ABD_R[j].MPI_SOURCE;
      gp = GP_ABC2ABD_R_A[IDR];

      for (i=0; i<Num_Rcv_Grid_ABC2ABD[IDR]; i++){
        BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = array1[gp+i].r;
	Rhor[BN_ABD].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_ABD],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  /* Local data */
 
  IDR = myid;
  for (i=0; i<Num_Snd_Grid_ABC2ABD[IDR]; i++){
    BN_ABC = Index_Snd_Grid_ABC2ABD[IDR][i];
    BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
    Rhor[BN_ABD].r = Rhok[BN_ABC].r;
    Rhor[BN_ABD].i = Rhok[BN_ABC].i;
    /*
    memcpy(&Rhor[BN_ABD],&Rhok[BN_ABC],SIZE_DCOMPLEX);
    */
  }

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_ABC2ABD_S,ST_ABC2ABD_S);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

break;

} /* opt == 3 */

case 4: {

  blocksize = BLOCK_SIZE;
  if(blocksize > numprocs) blocksize = numprocs;
  
  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for(n=0; n<numprocs; n+=blocksize){
   
  NN_S = 0;
  NN_R = 0;
  recvCount = 0;
  bsize = numprocs-n < blocksize ? numprocs-n : blocksize;

  /* MPI_Irecv */

  for (b=0; b<bsize; b++){

    IDR = (myid+b+n) % numprocs;
    gp = GP_ABC2ABD_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_ABC2ABD[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_ABC2ABD_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDS]; i++){
        BN_ABC = Index_Snd_Grid_ABC2ABD[IDS][i];
	array0[gp+i].r = Rhok[BN_ABC].r;
	array0[gp+i].i = Rhok[BN_ABC].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_ABC2ABD[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;  
    }
  }
  
  /* MPI_Waitsome */

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,request_recv,&recvNum,recvID,stat_recv);
    for(j=0;j<recvNum;j++){
      IDR = stat_recv[j].MPI_SOURCE;
      gp = GP_ABC2ABD_R_A[IDR];

      for (i=0; i<Num_Rcv_Grid_ABC2ABD[IDR]; i++){
        BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = array1[gp+i].r;
	Rhor[BN_ABD].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_ABD],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);

  } /* n */

  /* Local data */
 
  IDR = myid;
  for (i=0; i<Num_Snd_Grid_ABC2ABD[IDR]; i++){
    BN_ABC = Index_Snd_Grid_ABC2ABD[IDR][i];
    BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
    Rhor[BN_ABD].r = Rhok[BN_ABC].r;
    Rhor[BN_ABD].i = Rhok[BN_ABC].i;
    /*
    memcpy(&Rhor[BN_ABD],&Rhok[BN_ABC],SIZE_DCOMPLEX);
    */
  }

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

break;

} /* opt = 4*/


case 5: {

  if (OFFT_measure_time==1) openfft_dtime(&time0);
  
  for(n=0; n<numprocs; n++){
    IDR = (myid - n + numprocs) % numprocs;
    IDS = (myid + n) % numprocs;

    if(GP_ABC2ABD_R_A[IDR]>-1 && IDR != myid){
      recvbuf = &array1[GP_ABC2ABD_R_A[IDR]];
      recvcnt = Num_Rcv_Grid_ABC2ABD[IDR]*2;
    }
    else{
      IDR = MPI_PROC_NULL;
      recvcnt = 0;
      recvbuf = NULL;
    }

    if(GP_ABC2ABD_S_A[IDS]>-1 && IDS != myid){
      sendbuf = &array0[GP_ABC2ABD_S_A[IDS]];
      sendcnt = Num_Snd_Grid_ABC2ABD[IDS]*2;

      gp = GP_ABC2ABD_S_A[IDS];
      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDS]; i++){
        BN_ABC = Index_Snd_Grid_ABC2ABD[IDS][i];
	array0[gp+i].r = Rhok[BN_ABC].r;
	array0[gp+i].i = Rhok[BN_ABC].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      } 
    }
    else{
      IDS = MPI_PROC_NULL;
      sendbuf = NULL;
      sendcnt = 0;
    }
    MPI_Sendrecv(sendbuf, sendcnt, MPI_DOUBLE, IDS, tag,
		 recvbuf, recvcnt, MPI_DOUBLE, IDR, tag, 
		 MPI_COMM_WORLD, &stat);
  }
   
  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

  /* copy them to Rhor */

  for (ID=0; ID<NN_ABC2ABD_R; ID++){

    IDR = ID_NN_ABC2ABD_R[ID];
    gp = GP_ABC2ABD_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABC2ABD[IDR]; i++){
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = array1[gp+i].r;
	Rhor[BN_ABD].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_ABD],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDR]; i++){
	BN_ABC = Index_Snd_Grid_ABC2ABD[IDR][i];
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = Rhok[BN_ABC].r;
	Rhor[BN_ABD].i = Rhok[BN_ABC].i;
	/*
	memcpy(&Rhor[BN_ABD],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt == 5 */

case 60: {

  for (ID=0; ID<NN_ABC2ABD_S; ID++){

    IDS = ID_NN_ABC2ABD_S[ID];
    gp = IDS*Max_Num_Snd_Grid_ABC2ABD_A2A;

    if(IDS!=myid){
      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDS]; i++){
	BN_ABC = Index_Snd_Grid_ABC2ABD[IDS][i];
	array0[gp+i].r = Rhok[BN_ABC].r;
	array0[gp+i].i = Rhok[BN_ABC].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      } 
    }
  }

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  MPI_Alltoall(array0, Max_Num_Snd_Grid_ABC2ABD_A2A, MPI_DOUBLE_COMPLEX,
	       array1, Max_Num_Rcv_Grid_ABC2ABD_A2A, MPI_DOUBLE_COMPLEX, 
	       mpi_comm_level1);
  /*
  MPI_Alltoall(array0, Max_Num_Snd_Grid_ABC2ABD_A2A*2, MPI_DOUBLE,
	       array1, Max_Num_Rcv_Grid_ABC2ABD_A2A*2, MPI_DOUBLE, 
	       mpi_comm_level1);
  */
  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

  /* copy them to Rhor */

  for (ID=0; ID<NN_ABC2ABD_R; ID++){

    IDR = ID_NN_ABC2ABD_R[ID];
    gp = IDR*Max_Num_Rcv_Grid_ABC2ABD_A2A;

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABC2ABD[IDR]; i++){
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = array1[gp+i].r;
	Rhor[BN_ABD].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_ABD],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    else{
      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDR]; i++){
	BN_ABC = Index_Snd_Grid_ABC2ABD[IDR][i];
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = Rhok[BN_ABC].r;
	Rhor[BN_ABD].i = Rhok[BN_ABC].i;
	/*
	memcpy(&Rhor[BN_ABD],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt == 60 */

default: {

  NN_S = 0;
  NN_R = 0;

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }
   
  /* MPI_Irecv */

  for (ID=0; ID<NN_ABC2ABD_R; ID++){

    IDR = ID_NN_ABC2ABD_R[ID];
    gp = GP_ABC2ABD_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_ABC2ABD[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_ABC2ABD_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */
  

  for (ID=0; ID<NN_ABC2ABD_S; ID++){

    IDS = ID_NN_ABC2ABD_S[ID];
    gp = GP_ABC2ABD_S[ID];

    if (IDS!=myid){

      if (OFFT_measure_time==1) openfft_dtime(&time00);

      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDS]; i++){
        BN_ABC = Index_Snd_Grid_ABC2ABD[IDS][i];
	array0[gp+i].r = Rhok[BN_ABC].r;
	array0[gp+i].i = Rhok[BN_ABC].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      }     

      if (OFFT_measure_time==1){
	openfft_dtime(&time11);
	timeT += time11 - time00;
      }

      MPI_Isend( &array0[gp], Num_Snd_Grid_ABC2ABD[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_ABC2ABD_S[NN_S]);
      NN_S++;  
    }
  }
  
  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_ABC2ABD_S,ST_ABC2ABD_S);
  if (NN_R!=0) MPI_Waitall(NN_R,RQ_ABC2ABD_R,ST_ABC2ABD_R);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += timeT;
    timeC += time1 - time0 - timeT;
  }

  /* copy them to Rhor */

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for (ID=0; ID<NN_ABC2ABD_R; ID++){

    IDR = ID_NN_ABC2ABD_R[ID];
    gp = GP_ABC2ABD_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABC2ABD[IDR]; i++){
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = array1[gp+i].r;
	Rhor[BN_ABD].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_ABD],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_ABC2ABD[IDR]; i++){
	BN_ABC = Index_Snd_Grid_ABC2ABD[IDR][i];
	BN_ABD = Index_Rcv_Grid_ABC2ABD[IDR][i];
	Rhor[BN_ABD].r = Rhok[BN_ABC].r;
	Rhor[BN_ABD].i = Rhok[BN_ABC].i;
	/*
	memcpy(&Rhor[BN_ABD],&Rhok[BN_ABC],SIZE_DCOMPLEX);
	*/
      }
    }
  }

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += time1 - time0;
  }


break;

} /* opt else */

} /* switch */

  /*------------------ FFT along the C-axis in the ABD partition ------------------*/

 if (OFFT_measure_time==1){
   openfft_dtime(&time0);    
 }

 fftw_execute_dft(OFFT_p_mC,(fftw_complex*)Rhor,(fftw_complex*)Rhok);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeF += time1 - time0;
  }

  /*
  if(!myid){
    printf("\nABD\n");
    for(i=0;i<My_NumGrid_ABD;i++){
      printf("%d %f\n",i,Rhok[i]);
    }
    printf("\nAAAAAAAAAA\n");
  }


  /**/

  /*------------------ MPI: ABD to DCA partitions ------------------*/

switch(COMM_PATT){

case 1:{

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }

  for(i=0;i<numprocs;i+=1){
    count_snd[i] = 0;
    count_rcv[i] = 0;
    displ_snd[i] = 0;
    displ_rcv[i] = 0;
  }

  for (ID=0; ID<NN_ABD2DCA_S; ID+=1){

    IDS = ID_NN_ABD2DCA_S[ID];
    gp = GP_ABD2DCA_S[ID];

    count_snd[IDS] = Num_Snd_Grid_ABD2DCA[IDS]*2;
    displ_snd[IDS] = gp*2; 

    if(IDS!=myid){
      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDS]; i++){
	BN_ABD = Index_Snd_Grid_ABD2DCA[IDS][i];
	array0[gp+i].r = Rhok[BN_ABD].r;
	array0[gp+i].i = Rhok[BN_ABD].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      } 
    }
  }
 
  for (ID=0; ID<NN_ABD2DCA_R; ID+=1){
    
    IDR = ID_NN_ABD2DCA_R[ID];
    gp = GP_ABD2DCA_R[ID];

    count_rcv[IDR] = Num_Rcv_Grid_ABD2DCA[IDR]*2;
    displ_rcv[IDR] = gp*2; 

  }

  count_snd[myid] = 0;
  count_rcv[myid] = 0;

  if (OFFT_measure_time==1) openfft_dtime(&time00);

  MPI_Alltoallv(array0, count_snd, displ_snd, MPI_DOUBLE,
                       array1, count_rcv, displ_rcv, MPI_DOUBLE, mpi_comm_level1);

  if (OFFT_measure_time==1){
    openfft_dtime(&time11);
    timeT += time11 - time00;
    timeC += time11 - time00;
  }

  /* copy them to Rhor */

  for (ID=0; ID<NN_ABD2DCA_R; ID+=1){

    IDR = ID_NN_ABD2DCA_R[ID];
    gp = GP_ABD2DCA_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABD2DCA[IDR]; i++){
	BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = array1[gp+i].r;
	Rhor[BN_DCA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    else{
      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDR]; i++){
	BN_ABD = Index_Snd_Grid_ABD2DCA[IDR][i];
	BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = Rhok[BN_ABD].r;
	Rhor[BN_DCA].i = Rhok[BN_ABD].i;
	/*
	memcpy(&Rhor[BN_DCA],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }
    }
  }

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += time1 - time0 - timeT;
  }

break;

} /* opt == 1 */

case 2: {

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for(n=0; n<numprocs; n+=blocksize){
   
  NN_RS = 0;
  bsize = numprocs-n < blocksize ? numprocs-n : blocksize;

  /* MPI_Irecv */

  for (b=0; b<bsize; b++){

    IDR = (myid+b+n) % numprocs;
    gp = GP_ABD2DCA_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_ABD2DCA[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_arr[NN_RS]);
      NN_RS++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_ABD2DCA_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDS]; i++){
        BN_ABD = Index_Snd_Grid_ABD2DCA[IDS][i];
	array0[gp+i].r = Rhok[BN_ABD].r;
	array0[gp+i].i = Rhok[BN_ABD].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_ABD2DCA[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_arr[NN_RS]);
      NN_RS++;  
    }
  }
  
  /* MPI_Waitall */

  MPI_Waitall(NN_RS,request_arr,stat_arr);

  } /* n */

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

  for (ID=0; ID<NN_ABD2DCA_R; ID++){

    IDR = ID_NN_ABD2DCA_R[ID];
    gp = GP_ABD2DCA_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABD2DCA[IDR]; i++){
        BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = array1[gp+i].r;
	Rhor[BN_DCA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDR]; i++){
	BN_ABD = Index_Snd_Grid_ABD2DCA[IDR][i];
	BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = Rhok[BN_ABD].r;
	Rhor[BN_DCA].i = Rhok[BN_ABD].i;
	/*
	memcpy(&Rhor[BN_DCA],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt = 2 */

case 3: {

  NN_S = 0;
  NN_R = 0;
  recvCount = 0;

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  /* MPI_Irecv */

  for (ID=0; ID<NN_ABD2DCA_R; ID++){

    IDR = ID_NN_ABD2DCA_R[ID];
    gp = GP_ABD2DCA_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_ABD2DCA[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_ABD2DCA_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_ABD2DCA_S; ID++){

    IDS = ID_NN_ABD2DCA_S[ID];
    gp = GP_ABD2DCA_S[ID];

    if (IDS!=myid){
      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDS]; i++){
        BN_ABD = Index_Snd_Grid_ABD2DCA[IDS][i];
        array0[gp+i].r = Rhok[BN_ABD].r;
        array0[gp+i].i = Rhok[BN_ABD].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }     
      MPI_Isend( &array0[gp], Num_Snd_Grid_ABD2DCA[IDS]*2, 
  	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_ABD2DCA_S[NN_S]);
      NN_S++;  
    }
  }

  /* MPI_Waitsome */

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,RQ_ABD2DCA_R,&recvNum,recvID,ST_ABD2DCA_R);
    for(j=0;j<recvNum;j++){
      IDR = ST_ABD2DCA_R[j].MPI_SOURCE;
      gp = GP_ABD2DCA_R_A[IDR];
      for (i=0; i<Num_Rcv_Grid_ABD2DCA[IDR]; i++){
        BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
        Rhor[BN_DCA].r = array1[gp+i].r;
        Rhor[BN_DCA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  IDR = myid;
  for (i=0; i<Num_Snd_Grid_ABD2DCA[IDR]; i++){
    BN_ABD = Index_Snd_Grid_ABD2DCA[IDR][i];
    BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
    Rhor[BN_DCA].r = Rhok[BN_ABD].r;
    Rhor[BN_DCA].i = Rhok[BN_ABD].i;
    /*
    memcpy(&Rhor[BN_DCA],&Rhok[BN_ABD],SIZE_DCOMPLEX);
    */
  }

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_ABD2DCA_S,ST_ABD2DCA_S);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

break;

} /* opt == 3 */

case 4: {

  blocksize = BLOCK_SIZE;
  if(blocksize > numprocs) blocksize = numprocs;
  
  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for(n=0; n<numprocs; n+=blocksize){
   
  NN_S = 0;
  NN_R = 0;
  recvCount = 0;
  bsize = numprocs-n < blocksize ? numprocs-n : blocksize;

  /* MPI_Irecv */

  for (b=0; b<bsize; b++){

    IDR = (myid+b+n) % numprocs;
    gp = GP_ABD2DCA_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_ABD2DCA[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_ABD2DCA_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDS]; i++){
        BN_ABD = Index_Snd_Grid_ABD2DCA[IDS][i];
	array0[gp+i].r = Rhok[BN_ABD].r;
	array0[gp+i].i = Rhok[BN_ABD].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_ABD2DCA[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;  
    }
  }
  
  /* MPI_Waitsome */

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,request_recv,&recvNum,recvID,stat_recv);
    for(j=0;j<recvNum;j++){
      IDR = stat_recv[j].MPI_SOURCE;
      gp = GP_ABD2DCA_R_A[IDR];

      for (i=0; i<Num_Rcv_Grid_ABD2DCA[IDR]; i++){
        BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = array1[gp+i].r;
	Rhor[BN_DCA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);

  } /* n */

  /* Local data */

  IDR = myid;
  for (i=0; i<Num_Snd_Grid_ABD2DCA[IDR]; i++){
    BN_ABD = Index_Snd_Grid_ABD2DCA[IDR][i];
    BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
    Rhor[BN_DCA].r = Rhok[BN_ABD].r;
    Rhor[BN_DCA].i = Rhok[BN_ABD].i;
    /*
    memcpy(&Rhor[BN_DCA],&Rhok[BN_ABD],SIZE_DCOMPLEX);
    */
  }

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

break;

} /* opt = 4 */

case 5: {

  if (OFFT_measure_time==1) openfft_dtime(&time0);
  
  for(n=0; n<numprocs; n++){
    IDR = (myid - n + numprocs) % numprocs;
    IDS = (myid + n) % numprocs;

    if(GP_ABD2DCA_R_A[IDR]>-1 && IDR != myid){
      recvbuf = &array1[GP_ABD2DCA_R_A[IDR]];
      recvcnt = Num_Rcv_Grid_ABD2DCA[IDR]*2;
    }
    else{
      IDR = MPI_PROC_NULL;
      recvbuf = NULL;
      recvcnt = 0;
    }

    if(GP_ABD2DCA_S_A[IDS]>-1 && IDS != myid){
      sendbuf = &array0[GP_ABD2DCA_S_A[IDS]];
      sendcnt = Num_Snd_Grid_ABD2DCA[IDS]*2;
      
      gp = GP_ABD2DCA_S_A[IDS];
      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDS]; i++){
        BN_ABD = Index_Snd_Grid_ABD2DCA[IDS][i];
	array0[gp+i].r = Rhok[BN_ABD].r;
	array0[gp+i].i = Rhok[BN_ABD].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }     
    }
    else{
      IDS = MPI_PROC_NULL;
      sendbuf = NULL;
      sendcnt = 0;
    }
    MPI_Sendrecv(sendbuf, sendcnt, MPI_DOUBLE, IDS, tag,
		 recvbuf, recvcnt, MPI_DOUBLE, IDR, tag, 
		 MPI_COMM_WORLD, &stat);
  }
   
  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

  /* copy them to Rhok */

  for (ID=0; ID<NN_ABD2DCA_R; ID++){

    IDR = ID_NN_ABD2DCA_R[ID];
    gp = GP_ABD2DCA_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABD2DCA[IDR]; i++){
        BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = array1[gp+i].r;
	Rhor[BN_DCA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDR]; i++){
	BN_ABD = Index_Snd_Grid_ABD2DCA[IDR][i];
	BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = Rhok[BN_ABD].r;
	Rhor[BN_DCA].i = Rhok[BN_ABD].i;
	/*
	memcpy(&Rhor[BN_DCA],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt == 5 */

case 60: {

  for (ID=0; ID<NN_ABD2DCA_S; ID++){

    IDS = ID_NN_ABD2DCA_S[ID];
    gp = IDS*Max_Num_Snd_Grid_ABD2DCA_A2A;

    if(IDS!=myid){
      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDS]; i++){
	BN_ABD = Index_Snd_Grid_ABD2DCA[IDS][i];
	array0[gp+i].r = Rhok[BN_ABD].r;
	array0[gp+i].i = Rhok[BN_ABD].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      } 
    }
  }

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  MPI_Alltoall(array0, Max_Num_Snd_Grid_ABD2DCA_A2A, MPI_DOUBLE_COMPLEX,
	       array1, Max_Num_Rcv_Grid_ABD2DCA_A2A, MPI_DOUBLE_COMPLEX, 
	       mpi_comm_level1);
  /*
  MPI_Alltoall(array0, Max_Num_Snd_Grid_ABD2DCA_A2A*2, MPI_DOUBLE,
	       array1, Max_Num_Rcv_Grid_ABD2DCA_A2A*2, MPI_DOUBLE, 
	       mpi_comm_level1);
  */

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

  /* copy them to Rhok */

  for (ID=0; ID<NN_ABD2DCA_R; ID++){

    IDR = ID_NN_ABD2DCA_R[ID];
    gp = IDR*Max_Num_Rcv_Grid_ABD2DCA_A2A;

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABD2DCA[IDR]; i++){
	BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = array1[gp+i].r;
	Rhor[BN_DCA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    else{
      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDR]; i++){
	BN_ABD = Index_Snd_Grid_ABD2DCA[IDR][i];
	BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = Rhok[BN_ABD].r;
	Rhor[BN_DCA].i = Rhok[BN_ABD].i;
	/*
	memcpy(&Rhor[BN_DCA],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt == 60 */

default: {

  NN_S = 0;
  NN_R = 0;

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }

  /* MPI_Irecv */

  for (ID=0; ID<NN_ABD2DCA_R; ID++){

    IDR = ID_NN_ABD2DCA_R[ID];
    gp = GP_ABD2DCA_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_ABD2DCA[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_ABD2DCA_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_ABD2DCA_S; ID++){

    IDS = ID_NN_ABD2DCA_S[ID];
    gp = GP_ABD2DCA_S[ID];

    if (IDS!=myid){

      if (OFFT_measure_time==1) openfft_dtime(&time00);

      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDS]; i++){
        BN_ABD = Index_Snd_Grid_ABD2DCA[IDS][i];
	array0[gp+i].r = Rhok[BN_ABD].r;
	array0[gp+i].i = Rhok[BN_ABD].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }     

      if (OFFT_measure_time==1){
	openfft_dtime(&time11);
	timeT += time11 - time00;
      }

      MPI_Isend( &array0[gp], Num_Snd_Grid_ABD2DCA[IDS]*2, 
  	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_ABD2DCA_S[NN_S]);
      NN_S++;  
    }
  }

  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_ABD2DCA_S,ST_ABD2DCA_S);
  if (NN_R!=0) MPI_Waitall(NN_R,RQ_ABD2DCA_R,ST_ABD2DCA_R);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += timeT;
    timeC += time1 - time0 - timeT;
  }

  /* copy them to Rhor */

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for (ID=0; ID<NN_ABD2DCA_R; ID++){

    IDR = ID_NN_ABD2DCA_R[ID];
    gp = GP_ABD2DCA_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_ABD2DCA[IDR]; i++){
        BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = array1[gp+i].r;
	Rhor[BN_DCA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_ABD2DCA[IDR]; i++){
	BN_ABD = Index_Snd_Grid_ABD2DCA[IDR][i];
	BN_DCA = Index_Rcv_Grid_ABD2DCA[IDR][i];
	Rhor[BN_DCA].r = Rhok[BN_ABD].r;
	Rhor[BN_DCA].i = Rhok[BN_ABD].i;
	/*
	memcpy(&Rhor[BN_DCA],&Rhok[BN_ABD],SIZE_DCOMPLEX);
	*/
      }
    }
  }

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += time1 - time0;
  }

break;

} /* opt else */

} /* switch */

  /*------------------ FFT along the B-axis in the DCA partition ------------------*/

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);    
  }
   
  fftw_execute_dft(OFFT_p_mB,(fftw_complex*)Rhor,(fftw_complex*)Rhok);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeF += time1 - time0;
  }

  /*
  if(!myid){
    printf("\nDCA\n");
    for(i=0;i<My_NumGrid_DCA;i++){
      printf("%d %f\n",i,Rhor[i]);
    }
  }
  /**/

  /*------------------ MPI: DCA to DCB partitions ------------------*/

switch (COMM_PATT){

case 1: {

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }

  for(i=0;i<numprocs;i+=1){
    count_snd[i] = 0;
    count_rcv[i] = 0;
    displ_snd[i] = 0;
    displ_rcv[i] = 0;
  }

  for (ID=0; ID<NN_DCA2DCB_S; ID+=1){

    IDS = ID_NN_DCA2DCB_S[ID];
    gp = GP_DCA2DCB_S[ID];

    count_snd[IDS] = Num_Snd_Grid_DCA2DCB[IDS]*2;
    displ_snd[IDS] = gp*2; 

    if(IDS!=myid){
      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDS]; i++){
	BN_DCA = Index_Snd_Grid_DCA2DCB[IDS][i];
	array0[gp+i].r = Rhok[BN_DCA].r;
	array0[gp+i].i = Rhok[BN_DCA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      } 
    }
  }

  for (ID=0; ID<NN_DCA2DCB_R; ID+=1){

    IDR = ID_NN_DCA2DCB_R[ID];
    gp = GP_DCA2DCB_R[ID];

    count_rcv[IDR] = Num_Rcv_Grid_DCA2DCB[IDR]*2;
    displ_rcv[IDR] = gp*2; 
  }

  count_snd[myid] = 0;
  count_rcv[myid] = 0;

  if (OFFT_measure_time==1) openfft_dtime(&time00);

  MPI_Alltoallv(array0, count_snd, displ_snd, MPI_DOUBLE,
		array1, count_rcv, displ_rcv, MPI_DOUBLE, mpi_comm_level1);

  if (OFFT_measure_time==1){
    openfft_dtime(&time11);
    timeT += time11 - time00;
    timeC += time11 - time00;
  }

  /* copy them to Rhor */

  for (ID=0; ID<NN_DCA2DCB_R; ID+=1){

    IDR = ID_NN_DCA2DCB_R[ID];
    gp = GP_DCA2DCB_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_DCA2DCB[IDR]; i++){
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = array1[gp+i].r;
	Rhor[BN_DCB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    else{
      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDR]; i++){
	BN_DCA = Index_Snd_Grid_DCA2DCB[IDR][i];
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = Rhok[BN_DCA].r;
	Rhor[BN_DCB].i = Rhok[BN_DCA].i;
	/*
	memcpy(&Rhor[BN_DCB],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      }
    }
  }

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += time1 - time0 - timeT;
  }

break;

} /* opt == 1 */

case 2: {
  
  blocksize = BLOCK_SIZE;
  if(blocksize > numprocs) blocksize = numprocs;
  
  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for(n=0; n<numprocs; n+=blocksize){
   
  NN_RS = 0;
  bsize = numprocs-n < blocksize ? numprocs-n : blocksize;

  /* MPI_Irecv */

  for (b=0; b<bsize; b++){

    IDR = (myid+b+n) % numprocs;
    gp = GP_DCA2DCB_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_DCA2DCB[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_arr[NN_RS]);
      NN_RS++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_DCA2DCB_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDS]; i++){
        BN_DCA = Index_Snd_Grid_DCA2DCB[IDS][i];
        array0[gp+i].r = Rhok[BN_DCA].r;
        array0[gp+i].i = Rhok[BN_DCA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_DCA2DCB[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_arr[NN_RS]);
      NN_RS++;  
    }
  }
  
  /* MPI_Waitall */

  MPI_Waitall(NN_RS,request_arr,stat_arr);

  } /* n */

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

  /* copy them to Rhor */

  for (ID=0; ID<NN_DCA2DCB_R; ID++){

    IDR = ID_NN_DCA2DCB_R[ID];
    gp = GP_DCA2DCB_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_DCA2DCB[IDR]; i++){
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = array1[gp+i].r;
	Rhor[BN_DCB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDR]; i++){
	BN_DCA = Index_Snd_Grid_DCA2DCB[IDR][i];
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = Rhok[BN_DCA].r;
	Rhor[BN_DCB].i = Rhok[BN_DCA].i;
	/*
	memcpy(&Rhor[BN_DCB],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt = 2*/

case 3: {

  NN_S = 0;
  NN_R = 0;
  recvCount = 0;

  if (OFFT_measure_time==1) openfft_dtime(&time0);
   
  /* MPI_Irecv */

  for (ID=0; ID<NN_DCA2DCB_R; ID++){

    IDR = ID_NN_DCA2DCB_R[ID];
    gp = GP_DCA2DCB_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_DCA2DCB[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_DCA2DCB_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_DCA2DCB_S; ID++){

    IDS = ID_NN_DCA2DCB_S[ID];
    gp = GP_DCA2DCB_S[ID];

    if (IDS!=myid){
      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDS]; i++){
        BN_DCA = Index_Snd_Grid_DCA2DCB[IDS][i];
        array0[gp+i].r = Rhok[BN_DCA].r;
        array0[gp+i].i = Rhok[BN_DCA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      }
      MPI_Isend( &array0[gp], Num_Snd_Grid_DCA2DCB[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_DCA2DCB_S[NN_S]);
      NN_S++;  
    }
  }

  /* MPI_Waitsome */ 

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,RQ_DCA2DCB_R,&recvNum,recvID,ST_DCA2DCB_R);
    for(j=0;j<recvNum;j++){
      IDR = ST_DCA2DCB_R[j].MPI_SOURCE;
      gp = GP_DCA2DCB_R_A[IDR];

      for (i=0; i<Num_Rcv_Grid_DCA2DCB[IDR]; i++){
        BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = array1[gp+i].r;
	Rhor[BN_DCB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  /* Local data */
 
  IDR = myid;
  for (i=0; i<Num_Snd_Grid_DCA2DCB[IDR]; i++){
    BN_DCA = Index_Snd_Grid_DCA2DCB[IDR][i];
    BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
    Rhor[BN_DCB].r = Rhok[BN_DCA].r;
    Rhor[BN_DCB].i = Rhok[BN_DCA].i;
    /*
    memcpy(&Rhor[BN_DCB],&Rhok[BN_DCA],SIZE_DCOMPLEX);
    */
  }

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_DCA2DCB_S,ST_DCA2DCB_S);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

break;

} /* opt == 3 */

case 4: {
  
  blocksize = BLOCK_SIZE;
  if(blocksize > numprocs) blocksize = numprocs;
  
  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for(n=0; n<numprocs; n+=blocksize){
   
  NN_S = 0;
  NN_R = 0;
  recvCount = 0;
  bsize = numprocs-n < blocksize ? numprocs-n : blocksize;

  /* MPI_Irecv */

  for (b=0; b<bsize; b++){

    IDR = (myid+b+n) % numprocs;
    gp = GP_DCA2DCB_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_DCA2DCB[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_DCA2DCB_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDS]; i++){
        BN_DCA = Index_Snd_Grid_DCA2DCB[IDS][i];
	array0[gp+i].r = Rhok[BN_DCA].r;
	array0[gp+i].i = Rhok[BN_DCA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_DCA2DCB[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;  
    }
  }
  
  /* MPI_Waitsome */

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,request_recv,&recvNum,recvID,stat_recv);
    for(j=0;j<recvNum;j++){
      IDR = stat_recv[j].MPI_SOURCE;
      gp = GP_DCA2DCB_R_A[IDR];

      for (i=0; i<Num_Rcv_Grid_DCA2DCB[IDR]; i++){
        BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = array1[gp+i].r;
	Rhor[BN_DCB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);

  } /* n */

  /* Local data */
 
  IDR = myid;
  for (i=0; i<Num_Snd_Grid_DCA2DCB[IDR]; i++){
    BN_DCA = Index_Snd_Grid_DCA2DCB[IDR][i];
    BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
    Rhor[BN_DCB].r = Rhok[BN_DCA].r;
    Rhor[BN_DCB].i = Rhok[BN_DCA].i;
    /*
    memcpy(&Rhor[BN_DCB],&Rhok[BN_DCA],SIZE_DCOMPLEX);
    */
  }

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

break;

} /* opt = 4*/


case 5: {

  if (OFFT_measure_time==1) openfft_dtime(&time0);
  
  for(n=0; n<numprocs; n++){
    IDR = (myid - n + numprocs) % numprocs;
    IDS = (myid + n) % numprocs;

    if(GP_DCA2DCB_R_A[IDR]>-1 && IDR != myid){
      recvbuf = &array1[GP_DCA2DCB_R_A[IDR]];
      recvcnt = Num_Rcv_Grid_DCA2DCB[IDR]*2;
    }
    else{
      IDR = MPI_PROC_NULL;
      recvcnt = 0;
      recvbuf = NULL;
    }

    if(GP_DCA2DCB_S_A[IDS]>-1 && IDS != myid){
      sendbuf = &array0[GP_DCA2DCB_S_A[IDS]];
      sendcnt = Num_Snd_Grid_DCA2DCB[IDS]*2;

      gp = GP_DCA2DCB_S_A[IDS];
      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDS]; i++){
        BN_DCA = Index_Snd_Grid_DCA2DCB[IDS][i];
	array0[gp+i].r = Rhok[BN_DCA].r;
	array0[gp+i].i = Rhok[BN_DCA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      } 
    }
    else{
      IDS = MPI_PROC_NULL;
      sendbuf = NULL;
      sendcnt = 0;
    }
    MPI_Sendrecv(sendbuf, sendcnt, MPI_DOUBLE, IDS, tag,
		 recvbuf, recvcnt, MPI_DOUBLE, IDR, tag, 
		 MPI_COMM_WORLD, &stat);
  }
   
  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

  /* copy them to Rhor */

  for (ID=0; ID<NN_DCA2DCB_R; ID++){

    IDR = ID_NN_DCA2DCB_R[ID];
    gp = GP_DCA2DCB_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_DCA2DCB[IDR]; i++){
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = array1[gp+i].r;
	Rhor[BN_DCB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDR]; i++){
	BN_DCA = Index_Snd_Grid_DCA2DCB[IDR][i];
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = Rhok[BN_DCA].r;
	Rhor[BN_DCB].i = Rhok[BN_DCA].i;
	/*
	memcpy(&Rhor[BN_DCB],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt == 5 */

case 60: {

  for (ID=0; ID<NN_DCA2DCB_S; ID++){

    IDS = ID_NN_DCA2DCB_S[ID];
    gp = IDS*Max_Num_Snd_Grid_DCA2DCB_A2A;

    if(IDS!=myid){
      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDS]; i++){
	BN_DCA = Index_Snd_Grid_DCA2DCB[IDS][i];
	array0[gp+i].r = Rhok[BN_DCA].r;
	array0[gp+i].i = Rhok[BN_DCA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      } 
    }
  }

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  MPI_Alltoall(array0, Max_Num_Snd_Grid_DCA2DCB_A2A, MPI_DOUBLE_COMPLEX,
	       array1, Max_Num_Rcv_Grid_DCA2DCB_A2A, MPI_DOUBLE_COMPLEX, 
	       mpi_comm_level1);
  /*
  MPI_Alltoall(array0, Max_Num_Snd_Grid_DCA2DCB_A2A*2, MPI_DOUBLE,
	       array1, Max_Num_Rcv_Grid_DCA2DCB_A2A*2, MPI_DOUBLE, 
	       mpi_comm_level1);
  */
  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeC += time1 - time0;
  }

  /* copy them to Rhor */

  for (ID=0; ID<NN_DCA2DCB_R; ID++){

    IDR = ID_NN_DCA2DCB_R[ID];
    gp = IDR*Max_Num_Rcv_Grid_DCA2DCB_A2A;

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_DCA2DCB[IDR]; i++){
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = array1[gp+i].r;
	Rhor[BN_DCB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    else{
      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDR]; i++){
	BN_DCA = Index_Snd_Grid_DCA2DCB[IDR][i];
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = Rhok[BN_DCA].r;
	Rhor[BN_DCB].i = Rhok[BN_DCA].i;
	/*
	memcpy(&Rhor[BN_DCB],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt == 60 */

default: {

  NN_S = 0;
  NN_R = 0;

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }
   
  /* MPI_Irecv */

  for (ID=0; ID<NN_DCA2DCB_R; ID++){

    IDR = ID_NN_DCA2DCB_R[ID];
    gp = GP_DCA2DCB_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_DCA2DCB[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_DCA2DCB_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */
  

  for (ID=0; ID<NN_DCA2DCB_S; ID++){

    IDS = ID_NN_DCA2DCB_S[ID];
    gp = GP_DCA2DCB_S[ID];

    if (IDS!=myid){

      if (OFFT_measure_time==1) openfft_dtime(&time00);

      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDS]; i++){
        BN_DCA = Index_Snd_Grid_DCA2DCB[IDS][i];
	array0[gp+i].r = Rhok[BN_DCA].r;
	array0[gp+i].i = Rhok[BN_DCA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      }     

      if (OFFT_measure_time==1){
	openfft_dtime(&time11);
	timeT += time11 - time00;
      }

      MPI_Isend( &array0[gp], Num_Snd_Grid_DCA2DCB[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_DCA2DCB_S[NN_S]);
      NN_S++;  
    }
  }
  
  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_DCA2DCB_S,ST_DCA2DCB_S);
  if (NN_R!=0) MPI_Waitall(NN_R,RQ_DCA2DCB_R,ST_DCA2DCB_R);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += timeT;
    timeC += time1 - time0 - timeT;
  }

  /* copy them to Rhor */

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for (ID=0; ID<NN_DCA2DCB_R; ID++){

    IDR = ID_NN_DCA2DCB_R[ID];
    gp = GP_DCA2DCB_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_DCA2DCB[IDR]; i++){
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = array1[gp+i].r;
	Rhor[BN_DCB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_DCB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_DCA2DCB[IDR]; i++){
	BN_DCA = Index_Snd_Grid_DCA2DCB[IDR][i];
	BN_DCB = Index_Rcv_Grid_DCA2DCB[IDR][i];
	Rhor[BN_DCB].r = Rhok[BN_DCA].r;
	Rhor[BN_DCB].i = Rhok[BN_DCA].i;
	/*
	memcpy(&Rhor[BN_DCB],&Rhok[BN_DCA],SIZE_DCOMPLEX);
	*/
      }
    }
  }

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += time1 - time0;
  }


break;

} /* opt else */

} /* switch */

  /*------------------ FFT along the A-axis in the DCB partition ------------------*/

 if (OFFT_measure_time==1){
   openfft_dtime(&time0);    
 }

  fftw_execute_dft(OFFT_p_mA,(fftw_complex*)Rhor,(fftw_complex*)Rhok);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeF += time1 - time0;

    openfft_dtime(&Etime_proc);
    time0 = Etime_proc -Stime_proc;

    timeD = timeDf + timeDc;
    OFFT_timers[0] = time0;
    OFFT_timers[1] = timeC;
    OFFT_timers[2] = timeF;
    OFFT_timers[3] = timeD;
    OFFT_timers[4] = time0 - timeC - timeF - timeD;
    OFFT_timers[5] = timeDf;
    OFFT_timers[6] = timeDc;
  }

  return time0;
}

