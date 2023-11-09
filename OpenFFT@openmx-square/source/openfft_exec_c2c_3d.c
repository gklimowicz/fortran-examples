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
  openfft_exec_c2c_3d performs the c2c transformation of the 3D data. 
  Must be called after openfft_init_c2c_3d.

     Input : Rhor (double complex).

     Output: Rhok (double complex).

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "openfft.h"

double openfft_exec_c2c_3d(dcomplex *Rhor, dcomplex *Rhok) 
{
  long long int i,j,BN_AB,BN_CB,BN_CA,gp,NN_S,NN_R,NN_RS;
  dcomplex *sendbuf,*recvbuf;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  int blocksize,n,bsize,b,recvCount,recvNum,sendcnt,recvcnt;
  double Stime_proc,Etime_proc;
  double time0,time1,time00,time11,timeC,timeF,timeDf,timeDc,timeD,timeT;
  MPI_Status stat;
  MPI_Comm mpi_comm_level1 = MPI_COMM_WORLD;

  time0 = 0.0;
  timeC = 0.0;
  timeF = 0.0;
  timeD = 0.0;
  timeDf = 0.0;
  timeDc = 0.0;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /*------------------ FFT along the C-axis in the AB partition ------------------*/

  if (OFFT_measure_time==1){
    openfft_dtime(&Stime_proc);
    openfft_dtime(&time0);
    timeT = 0.0;
  }
   
  fftw_execute_dft(OFFT_p_mC,(fftw_complex*)Rhor,(fftw_complex*)Rhok);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeF += time1 - time0 - timeT;
  }
  /*------------------ MPI: AB to CA partitions ------------------*/

switch (COMM_PATT){

case 1: {

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }

  for(i=0;i<numprocs;i++){
    count_snd[i] = 0;
    count_rcv[i] = 0;
    displ_snd[i] = 0;
    displ_rcv[i] = 0;
  }

  for (ID=0; ID<NN_AB2CA_S; ID++){

    IDS = ID_NN_AB2CA_S[ID];
    gp = GP_AB2CA_S[ID];

    count_snd[IDS] = Num_Snd_Grid_AB2CA[IDS]*2;
    displ_snd[IDS] = gp*2; 

    if(IDS!=myid){
      for (i=0; i<Num_Snd_Grid_AB2CA[IDS]; i++){
	BN_AB = Index_Snd_Grid_AB2CA[IDS][i];
	array0[gp+i].r = Rhok[BN_AB].r;
	array0[gp+i].i = Rhok[BN_AB].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_AB],SIZE_DCOMPLEX);
	*/
      } 
    }
  }
  count_snd[myid] = 0;

  for (ID=0; ID<NN_AB2CA_R; ID++){

    IDR = ID_NN_AB2CA_R[ID];
    gp = GP_AB2CA_R[ID];

    count_rcv[IDR] = Num_Rcv_Grid_AB2CA[IDR]*2;
    displ_rcv[IDR] = gp*2; 
  }
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

  for (ID=0; ID<NN_AB2CA_R; ID++){

    IDR = ID_NN_AB2CA_R[ID];
    gp = GP_AB2CA_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_AB2CA[IDR]; i++){
	BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = array1[gp+i].r;
	Rhor[BN_CA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    else{
      for (i=0; i<Num_Snd_Grid_AB2CA[IDR]; i++){
	BN_AB = Index_Snd_Grid_AB2CA[IDR][i];
	BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = Rhok[BN_AB].r;
	Rhor[BN_CA].i = Rhok[BN_AB].i;
	/*
	memcpy(&Rhor[BN_CA],&Rhok[BN_AB],SIZE_DCOMPLEX);
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
    gp = GP_AB2CA_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_AB2CA[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_arr[NN_RS]);
      NN_RS++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_AB2CA_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_AB2CA[IDS]; i++){
        BN_AB = Index_Snd_Grid_AB2CA[IDS][i];
        array0[gp+i].r = Rhok[BN_AB].r;
        array0[gp+i].i = Rhok[BN_AB].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_AB],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_AB2CA[IDS]*2, 
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

  for (ID=0; ID<NN_AB2CA_R; ID++){

    IDR = ID_NN_AB2CA_R[ID];
    gp = GP_AB2CA_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_AB2CA[IDR]; i++){
	BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = array1[gp+i].r;
	Rhor[BN_CA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_AB2CA[IDR]; i++){
	BN_AB = Index_Snd_Grid_AB2CA[IDR][i];
	BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = Rhok[BN_AB].r;
	Rhor[BN_CA].i = Rhok[BN_AB].i;
	/*
	memcpy(&Rhor[BN_CA],&Rhok[BN_AB],SIZE_DCOMPLEX);
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

  for (ID=0; ID<NN_AB2CA_R; ID++){

    IDR = ID_NN_AB2CA_R[ID];
    gp = GP_AB2CA_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_AB2CA[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_AB2CA_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_AB2CA_S; ID++){

    IDS = ID_NN_AB2CA_S[ID];
    gp = GP_AB2CA_S[ID];

    if (IDS!=myid){
      for (i=0; i<Num_Snd_Grid_AB2CA[IDS]; i++){
        BN_AB = Index_Snd_Grid_AB2CA[IDS][i];
        array0[gp+i].r = Rhok[BN_AB].r;
        array0[gp+i].i = Rhok[BN_AB].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_AB],SIZE_DCOMPLEX);
	*/
      }
      MPI_Isend( &array0[gp], Num_Snd_Grid_AB2CA[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_AB2CA_S[NN_S]);
      NN_S++;  
    }
  }

  /* MPI_Waitsome */ 

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,RQ_AB2CA_R,&recvNum,recvID,ST_AB2CA_R);
    for(j=0;j<recvNum;j++){
      IDR = ST_AB2CA_R[j].MPI_SOURCE;
      gp = GP_AB2CA_R_A[IDR];

      for (i=0; i<Num_Rcv_Grid_AB2CA[IDR]; i++){
        BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = array1[gp+i].r;
	Rhor[BN_CA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  /* Local data */
 
  IDR = myid;
  for (i=0; i<Num_Snd_Grid_AB2CA[IDR]; i++){
    BN_AB = Index_Snd_Grid_AB2CA[IDR][i];
    BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
    Rhor[BN_CA].r = Rhok[BN_AB].r;
    Rhor[BN_CA].i = Rhok[BN_AB].i;
    /*
    memcpy(&Rhor[BN_CA],&Rhok[BN_AB],SIZE_DCOMPLEX);
    */
  }

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_AB2CA_S,ST_AB2CA_S);

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
    gp = GP_AB2CA_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_AB2CA[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_AB2CA_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_AB2CA[IDS]; i++){
        BN_AB = Index_Snd_Grid_AB2CA[IDS][i];
	array0[gp+i].r = Rhok[BN_AB].r;
	array0[gp+i].i = Rhok[BN_AB].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_AB],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_AB2CA[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;  
    }
  }
  
  /* MPI_Waitsome */

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,request_recv,&recvNum,recvID,stat_recv);
    for(j=0;j<recvNum;j++){
      IDR = stat_recv[j].MPI_SOURCE;
      gp = GP_AB2CA_R_A[IDR];

      for (i=0; i<Num_Rcv_Grid_AB2CA[IDR]; i++){
        BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = array1[gp+i].r;
	Rhor[BN_CA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);

  } /* n */

  /* Local data */
 
  IDR = myid;
  for (i=0; i<Num_Snd_Grid_AB2CA[IDR]; i++){
    BN_AB = Index_Snd_Grid_AB2CA[IDR][i];
    BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
    Rhor[BN_CA].r = Rhok[BN_AB].r;
    Rhor[BN_CA].i = Rhok[BN_AB].i;
    /*
    memcpy(&Rhor[BN_CA],&Rhok[BN_AB],SIZE_DCOMPLEX);
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

    if(GP_AB2CA_R_A[IDR]>-1 && IDR != myid){
      recvbuf = &array1[GP_AB2CA_R_A[IDR]];
      recvcnt = Num_Rcv_Grid_AB2CA[IDR]*2;
    }
    else{
      IDR = MPI_PROC_NULL;
      recvcnt = 0;
      recvbuf = NULL;
    }

    if(GP_AB2CA_S_A[IDS]>-1 && IDS != myid){
      sendbuf = &array0[GP_AB2CA_S_A[IDS]];
      sendcnt = Num_Snd_Grid_AB2CA[IDS]*2;

      gp = GP_AB2CA_S_A[IDS];
      for (i=0; i<Num_Snd_Grid_AB2CA[IDS]; i++){
        BN_AB = Index_Snd_Grid_AB2CA[IDS][i];
	array0[gp+i].r = Rhok[BN_AB].r;
	array0[gp+i].i = Rhok[BN_AB].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_AB],SIZE_DCOMPLEX);
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

  for (ID=0; ID<NN_AB2CA_R; ID++){

    IDR = ID_NN_AB2CA_R[ID];
    gp = GP_AB2CA_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_AB2CA[IDR]; i++){
	BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = array1[gp+i].r;
	Rhor[BN_CA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_AB2CA[IDR]; i++){
	BN_AB = Index_Snd_Grid_AB2CA[IDR][i];
	BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = Rhok[BN_AB].r;
	Rhor[BN_CA].i = Rhok[BN_AB].i;
	/*
	memcpy(&Rhor[BN_CA],&Rhok[BN_AB],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt == 5 */

default: {

  NN_S = 0;
  NN_R = 0;

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }
   
  /* MPI_Irecv */

  for (ID=0; ID<NN_AB2CA_R; ID++){

    IDR = ID_NN_AB2CA_R[ID];
    gp = GP_AB2CA_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_AB2CA[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_AB2CA_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */
  

  for (ID=0; ID<NN_AB2CA_S; ID++){

    IDS = ID_NN_AB2CA_S[ID];
    gp = GP_AB2CA_S[ID];

    if (IDS!=myid){

      if (OFFT_measure_time==1) openfft_dtime(&time00);

      for (i=0; i<Num_Snd_Grid_AB2CA[IDS]; i++){
        BN_AB = Index_Snd_Grid_AB2CA[IDS][i];
	array0[gp+i].r = Rhok[BN_AB].r;
	array0[gp+i].i = Rhok[BN_AB].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_AB],SIZE_DCOMPLEX);
	*/
      }     

      if (OFFT_measure_time==1){
	openfft_dtime(&time11);
	timeT += time11 - time00;
      }

      MPI_Isend( &array0[gp], Num_Snd_Grid_AB2CA[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_AB2CA_S[NN_S]);
      NN_S++;  
    }
  }
  
  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_AB2CA_S,ST_AB2CA_S);
  if (NN_R!=0) MPI_Waitall(NN_R,RQ_AB2CA_R,ST_AB2CA_R);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += timeT;
    timeC += time1 - time0 - timeT;
  }

  /* copy them to Rhor */

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for (ID=0; ID<NN_AB2CA_R; ID++){

    IDR = ID_NN_AB2CA_R[ID];
    gp = GP_AB2CA_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_AB2CA[IDR]; i++){
	BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = array1[gp+i].r;
	Rhor[BN_CA].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CA],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_AB2CA[IDR]; i++){
	BN_AB = Index_Snd_Grid_AB2CA[IDR][i];
	BN_CA = Index_Rcv_Grid_AB2CA[IDR][i];
	Rhor[BN_CA].r = Rhok[BN_AB].r;
	Rhor[BN_CA].i = Rhok[BN_AB].i;
	/*
	memcpy(&Rhor[BN_CA],&Rhok[BN_AB],SIZE_DCOMPLEX);
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

  /*------------------ FFT along the B-axis in the CA partition ------------------*/

 if (OFFT_measure_time==1){
   openfft_dtime(&time0);    
   timeT = 0.0;
 }

 fftw_execute_dft(OFFT_p_mB,(fftw_complex*)Rhor,(fftw_complex*)Rhok);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeF += time1 - time0 - timeT;
  }

  /*------------------ MPI: CA to CB partitions ------------------*/

switch(COMM_PATT){

case 1:{

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }

  for(i=0;i<numprocs;i++){
    count_snd[i] = 0;
    count_rcv[i] = 0;
    displ_snd[i] = 0;
    displ_rcv[i] = 0;
  }

  for (ID=0; ID<NN_CA2CB_S; ID++){

    IDS = ID_NN_CA2CB_S[ID];
    gp = GP_CA2CB_S[ID];

    count_snd[IDS] = Num_Snd_Grid_CA2CB[IDS]*2;
    displ_snd[IDS] = gp*2; 

    if(IDS!=myid){
      for (i=0; i<Num_Snd_Grid_CA2CB[IDS]; i++){
	BN_CA = Index_Snd_Grid_CA2CB[IDS][i];
	array0[gp+i].r = Rhok[BN_CA].r;
	array0[gp+i].i = Rhok[BN_CA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_CA],SIZE_DCOMPLEX);
	*/
      } 
    }
  }
  count_snd[myid] = 0;

  for (ID=0; ID<NN_CA2CB_R; ID++){

    IDR = ID_NN_CA2CB_R[ID];
    gp = GP_CA2CB_R[ID];

    count_rcv[IDR] = Num_Rcv_Grid_CA2CB[IDR]*2;
    displ_rcv[IDR] = gp*2; 
  }
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

  for (ID=0; ID<NN_CA2CB_R; ID++){

    IDR = ID_NN_CA2CB_R[ID];
    gp = GP_CA2CB_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_CA2CB[IDR]; i++){
	BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
	Rhor[BN_CB].r = array1[gp+i].r;
	Rhor[BN_CB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    else{
      for (i=0; i<Num_Snd_Grid_CA2CB[IDR]; i++){
	BN_CA = Index_Snd_Grid_CA2CB[IDR][i];
	BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
	Rhor[BN_CB].r = Rhok[BN_CA].r;
	Rhor[BN_CB].i = Rhok[BN_CA].i;
	/*
	memcpy(&Rhor[BN_CB],&Rhok[BN_CA],SIZE_DCOMPLEX);
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
    gp = GP_CA2CB_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_CA2CB[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_arr[NN_RS]);
      NN_RS++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_CA2CB_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_CA2CB[IDS]; i++){
        BN_CA = Index_Snd_Grid_CA2CB[IDS][i];
	array0[gp+i].r = Rhok[BN_CA].r;
	array0[gp+i].i = Rhok[BN_CA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_CA],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_CA2CB[IDS]*2, 
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

  for (ID=0; ID<NN_CA2CB_R; ID++){

    IDR = ID_NN_CA2CB_R[ID];
    gp = GP_CA2CB_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_CA2CB[IDR]; i++){
        BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
	Rhor[BN_CB].r = array1[gp+i].r;
	Rhor[BN_CB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_CA2CB[IDR]; i++){
	BN_CA = Index_Snd_Grid_CA2CB[IDR][i];
	BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
	Rhor[BN_CB].r = Rhok[BN_CA].r;
	Rhor[BN_CB].i = Rhok[BN_CA].i;
	/*
	memcpy(&Rhor[BN_CB],&Rhok[BN_CA],SIZE_DCOMPLEX);
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

  for (ID=0; ID<NN_CA2CB_R; ID++){

    IDR = ID_NN_CA2CB_R[ID];
    gp = GP_CA2CB_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_CA2CB[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_CA2CB_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_CA2CB_S; ID++){

    IDS = ID_NN_CA2CB_S[ID];
    gp = GP_CA2CB_S[ID];

    if (IDS!=myid){
      for (i=0; i<Num_Snd_Grid_CA2CB[IDS]; i++){
        BN_CA = Index_Snd_Grid_CA2CB[IDS][i];
        array0[gp+i].r = Rhok[BN_CA].r;
        array0[gp+i].i = Rhok[BN_CA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_CA],SIZE_DCOMPLEX);
	*/
      }     
      MPI_Isend( &array0[gp], Num_Snd_Grid_CA2CB[IDS]*2, 
  	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_CA2CB_S[NN_S]);
      NN_S++;  
    }
  }

  /* MPI_Waitsome */

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,RQ_CA2CB_R,&recvNum,recvID,ST_CA2CB_R);
    for(j=0;j<recvNum;j++){
      IDR = ST_CA2CB_R[j].MPI_SOURCE;
      gp = GP_CA2CB_R_A[IDR];
      for (i=0; i<Num_Rcv_Grid_CA2CB[IDR]; i++){
        BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
        Rhor[BN_CB].r = array1[gp+i].r;
        Rhor[BN_CB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  IDR = myid;
  for (i=0; i<Num_Snd_Grid_CA2CB[IDR]; i++){
    BN_CA = Index_Snd_Grid_CA2CB[IDR][i];
    BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
    Rhor[BN_CB].r = Rhok[BN_CA].r;
    Rhor[BN_CB].i = Rhok[BN_CA].i;
    /*
    memcpy(&Rhor[BN_CB],&Rhok[BN_CA],SIZE_DCOMPLEX);
    */
  }

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_CA2CB_S,ST_CA2CB_S);

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
    gp = GP_CA2CB_R_A[IDR];
    if (gp>-1 && IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_CA2CB[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (b=0; b<bsize; b++){

    IDS = (myid-b-n+numprocs) % numprocs;
    gp = GP_CA2CB_S_A[IDS];
    if (gp>-1 && IDS!=myid){

      for (i=0; i<Num_Snd_Grid_CA2CB[IDS]; i++){
        BN_CA = Index_Snd_Grid_CA2CB[IDS][i];
	array0[gp+i].r = Rhok[BN_CA].r;
	array0[gp+i].i = Rhok[BN_CA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_CA],SIZE_DCOMPLEX);
	*/
      }     

      MPI_Isend( &array0[gp], Num_Snd_Grid_CA2CB[IDS]*2, 
    	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;  
    }
  }
  
  /* MPI_Waitsome */

  while(recvCount < NN_R){
    MPI_Waitsome(NN_R,request_recv,&recvNum,recvID,stat_recv);
    for(j=0;j<recvNum;j++){
      IDR = stat_recv[j].MPI_SOURCE;
      gp = GP_CA2CB_R_A[IDR];

      for (i=0; i<Num_Rcv_Grid_CA2CB[IDR]; i++){
        BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
	Rhor[BN_CB].r = array1[gp+i].r;
	Rhor[BN_CB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }
    recvCount += recvNum;
  }

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);

  } /* n */

  /* Local data */

  IDR = myid;
  for (i=0; i<Num_Snd_Grid_CA2CB[IDR]; i++){
    BN_CA = Index_Snd_Grid_CA2CB[IDR][i];
    BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
    Rhor[BN_CB].r = Rhok[BN_CA].r;
    Rhor[BN_CB].i = Rhok[BN_CA].i;
    /*
    memcpy(&Rhor[BN_CB],&Rhok[BN_CA],SIZE_DCOMPLEX);
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

    if(GP_CA2CB_R_A[IDR]>-1 && IDR != myid){
      recvbuf = &array1[GP_CA2CB_R_A[IDR]];
      recvcnt = Num_Rcv_Grid_CA2CB[IDR]*2;
    }
    else{
      IDR = MPI_PROC_NULL;
      recvbuf = NULL;
      recvcnt = 0;
    }

    if(GP_CA2CB_S_A[IDS]>-1 && IDS != myid){
      sendbuf = &array0[GP_CA2CB_S_A[IDS]];
      sendcnt = Num_Snd_Grid_CA2CB[IDS]*2;
      
      gp = GP_CA2CB_S_A[IDS];
      for (i=0; i<Num_Snd_Grid_CA2CB[IDS]; i++){
        BN_CA = Index_Snd_Grid_CA2CB[IDS][i];
	array0[gp+i].r = Rhok[BN_CA].r;
	array0[gp+i].i = Rhok[BN_CA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_CA],SIZE_DCOMPLEX);
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

  for (ID=0; ID<NN_CA2CB_R; ID++){

    IDR = ID_NN_CA2CB_R[ID];
    gp = GP_CA2CB_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_CA2CB[IDR]; i++){
        BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
	Rhor[BN_CB].r = array1[gp+i].r;
	Rhor[BN_CB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_CA2CB[IDR]; i++){
	BN_CA = Index_Snd_Grid_CA2CB[IDR][i];
	BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
	Rhor[BN_CB].r = Rhok[BN_CA].r;
	Rhor[BN_CB].i = Rhok[BN_CA].i;
	/*
	memcpy(&Rhor[BN_CB],&Rhok[BN_CA],SIZE_DCOMPLEX);
	*/
      }
    }
  }

break;

} /* opt == 5 */

default: {

  NN_S = 0;
  NN_R = 0;

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);
    timeT = 0.0;
  }

  /* MPI_Irecv */

  for (ID=0; ID<NN_CA2CB_R; ID++){

    IDR = ID_NN_CA2CB_R[ID];
    gp = GP_CA2CB_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &array1[gp], Num_Rcv_Grid_CA2CB[IDR]*2, 
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &RQ_CA2CB_R[NN_R]);
      NN_R++; 
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_CA2CB_S; ID++){

    IDS = ID_NN_CA2CB_S[ID];
    gp = GP_CA2CB_S[ID];

    if (IDS!=myid){

      if (OFFT_measure_time==1) openfft_dtime(&time00);

      for (i=0; i<Num_Snd_Grid_CA2CB[IDS]; i++){
        BN_CA = Index_Snd_Grid_CA2CB[IDS][i];
	array0[gp+i].r = Rhok[BN_CA].r;
	array0[gp+i].i = Rhok[BN_CA].i;
	/*
	memcpy(&array0[gp+i],&Rhok[BN_CA],SIZE_DCOMPLEX);
	*/
      }     

      if (OFFT_measure_time==1){
	openfft_dtime(&time11);
	timeT += time11 - time00;
      }

      MPI_Isend( &array0[gp], Num_Snd_Grid_CA2CB[IDS]*2, 
  	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &RQ_CA2CB_S[NN_S]);
      NN_S++;  
    }
  }

  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,RQ_CA2CB_S,ST_CA2CB_S);
  if (NN_R!=0) MPI_Waitall(NN_R,RQ_CA2CB_R,ST_CA2CB_R);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeDc += timeT;
    timeC += time1 - time0 - timeT;
  }

  /* copy them to Rhor */

  if (OFFT_measure_time==1) openfft_dtime(&time0);

  for (ID=0; ID<NN_CA2CB_R; ID++){

    IDR = ID_NN_CA2CB_R[ID];
    gp = GP_CA2CB_R[ID];

    if (IDR!=myid){
      for (i=0; i<Num_Rcv_Grid_CA2CB[IDR]; i++){
        BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
	Rhor[BN_CB].r = array1[gp+i].r;
	Rhor[BN_CB].i = array1[gp+i].i;
	/*
	memcpy(&Rhor[BN_CB],&array1[gp+i],SIZE_DCOMPLEX);
	*/
      }
    }

    else{
      for (i=0; i<Num_Snd_Grid_CA2CB[IDR]; i++){
	BN_CA = Index_Snd_Grid_CA2CB[IDR][i];
	BN_CB = Index_Rcv_Grid_CA2CB[IDR][i];
	Rhor[BN_CB].r = Rhok[BN_CA].r;
	Rhor[BN_CB].i = Rhok[BN_CA].i;
	/*
	memcpy(&Rhor[BN_CB],&Rhok[BN_CA],SIZE_DCOMPLEX);
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

  /*------------------ FFT along the A-axis in the CB partition ------------------*/

  if (OFFT_measure_time==1){
    openfft_dtime(&time0);    
    timeT = 0.0;
  }
   
  fftw_execute_dft(OFFT_p_mA,(fftw_complex*)Rhor,(fftw_complex*)Rhok);

  if (OFFT_measure_time==1){
    openfft_dtime(&time1);
    timeF += time1 - time0 - timeT;

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

