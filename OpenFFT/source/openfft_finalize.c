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
  openfft_finalize cleans up the calculations. 
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openfft.h"

double openfft_finalize()
{

  int numprocs,myid,ID;
  double Stime_proc,Etime_proc,time0;
  MPI_Comm mpi_comm_level1 = MPI_COMM_WORLD;

  time0 = 0.0;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (OFFT_measure_time==1) openfft_dtime(&Stime_proc);
  
  if(Num_Snd_Grid_AB2CA != NULL){

   free(Num_Snd_Grid_AB2CA);
   free(Num_Rcv_Grid_AB2CA);
   free(Num_Snd_Grid_CA2CB);
   free(Num_Rcv_Grid_CA2CB);

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

   free(RQ_AB2CA_S);
   free(RQ_AB2CA_R);
   free(ST_AB2CA_S);
   free(ST_AB2CA_R);
   free(RQ_CA2CB_S);
   free(RQ_CA2CB_R);
   free(ST_CA2CB_S);
   free(ST_CA2CB_R);
   
   if(Rhot != NULL){
    free(Rhot);
    fftw_destroy_plan(OFFT_p_A);  
    fftw_destroy_plan(OFFT_p_B);  
    fftw_destroy_plan(OFFT_p_C);  
    fftw_free(OFFT_in);
    fftw_free(OFFT_rin);
    fftw_free(OFFT_out);
   }
   else{
     fftw_destroy_plan(OFFT_p_mA);  
     fftw_destroy_plan(OFFT_p_mB);  
     fftw_destroy_plan(OFFT_p_mC);  
     fftw_free(OFFT_min);
     fftw_free(OFFT_mout);
   }
  }

  if(Num_Snd_Grid_ABC2ABD != NULL){

   free(Num_Snd_Grid_ABC2ABD);
   free(Num_Rcv_Grid_ABC2ABD);
   free(Num_Snd_Grid_ABD2DCA);
   free(Num_Rcv_Grid_ABD2DCA);
   free(Num_Snd_Grid_DCA2DCB);
   free(Num_Rcv_Grid_DCA2DCB);

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

   free(RQ_ABC2ABD_S);
   free(RQ_ABC2ABD_R);
   free(ST_ABC2ABD_S);
   free(ST_ABC2ABD_R);
   free(RQ_ABD2DCA_S);
   free(RQ_ABD2DCA_R);
   free(ST_ABD2DCA_S);
   free(ST_ABD2DCA_R);
   free(RQ_DCA2DCB_S);
   free(RQ_DCA2DCB_R);
   free(ST_DCA2DCB_S);
   free(ST_DCA2DCB_R);
   
   fftw_destroy_plan(OFFT_p_mA);  
   fftw_destroy_plan(OFFT_p_mB);  
   fftw_destroy_plan(OFFT_p_mC);  
   fftw_destroy_plan(OFFT_p_mD);  

   fftw_free(OFFT_min);
   fftw_free(OFFT_mout);
  }
  
  free(count_snd);
  free(count_rcv);
  free(displ_snd);
  free(displ_rcv); 

  free(request_arr);
  free(stat_arr);
  free(recvID);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  free(array0);
  free(array1);
  
  fftw_cleanup();

  if (OFFT_measure_time==1){
    openfft_dtime(&Etime_proc);
    time0 = Etime_proc -Stime_proc;
  }

  return time0;
}
