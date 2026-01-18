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
  openfft.h: C header file.
***********************************************************************/

#include "mpi.h"
#include "omp.h"

#ifndef ___dcomplex_definition___
typedef struct { double r,i; } dcomplex;
#define ___dcomplex_definition___ 
#endif

#define BLOCK_SIZE 256
#define COMM_NUM 6
#define COMM_DEFAULT 3

#define SIZE_DOUBLE 8
#define SIZE_DCOMPLEX 16

#define NumTest 10

#include <fftw3.h> 

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
        with communication-computation with overlap.

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
double openfft_init_c2c_3d(int N1, int N2, int N3,
			   int *My_Max_NumGrid, 
			   int *My_NumGrid_In, int *My_Index_In,
			   int *My_NumGrid_Out, int *My_Index_Out,
			   int offt_measure,
			   int measure_time, int print_memory);


/*******************************************************
  openfft_init_r2c_3d prepares data for the r2c transformation of 3D data.

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
        with communication-computation with overlap.

   + 5: MPI_Sendrecv.

   + 6: MPI_Isend and MPI_Irecv.

   + Others: default communication, which is 3.

 - measure_time and print_memory (0: disabled, 1: enabled).

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
double openfft_init_r2c_3d(int N1, int N2, int N3, 
			   int *My_Max_NumGrid, 
			   int *My_NumGrid_In, int *My_Index_In,
			   int *My_NumGrid_Out, int *My_Index_Out,
			   int offt_measure,
			   int measure_time, int print_memory);


/*******************************************************
  openfft_init_c2c_4d prepares data for the c2c transformation of 4D data.

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
        with communication-computation with overlap.

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
double openfft_init_c2c_4d(int N1, int N2, int N3, int N4,
			   int *My_Max_NumGrid, 
			   int *My_NumGrid_In, int *My_Index_In,
			   int *My_NumGrid_Out, int *My_Index_Out,
			   int offt_measure,
			   int measure_time, int print_memory);


/**********************************************************************
  openfft_exec_c2c_3d performs the c2c transformation of the 3D data. 
  Must be called after openfft_init_c2c_3d.

     Input : Rhor (double complex).

     Output: Rhok (double complex).

***********************************************************************/
double openfft_exec_c2c_3d(dcomplex *Rhor, dcomplex *Rhok); 

/**********************************************************************
  openfft_exec_r2c_3d performs the r2c transformation of the 3D data. 
  Must be called after openfft_init_r2c_3d.

     Input : Rhor (double).

     Output: Rhok (double complex).

***********************************************************************/
double openfft_exec_r2c_3d(double *Rhor, dcomplex *Rhok); 


/**********************************************************************
  openfft_exec_c2c_4d performs the c2c transformation of the 4D data. 
  Must be called after openfft_init_c2c_4d.

     Input : Rhor (double complex).

     Output: Rhok (double complex).

***********************************************************************/
double openfft_exec_c2c_4d(dcomplex *Rhor, dcomplex *Rhok); 


/*******************************************************
  openfft_finalize cleans up the calculations. 
******************************************************/
double openfft_finalize();  

/*******************************************************
  openfft_dtime measures elapsed time. 
******************************************************/
void openfft_dtime(double *time);
 

/*******************************************************
 int *Num_Snd_Grid_AB2CA

  Num_Snd_Grid_AB2CA gives the number of grid points
  sent from the AB partition to the CA partition.
  size: Num_Snd_Grid_AB2CA[numprocs]
*******************************************************/
int *Num_Snd_Grid_AB2CA;
int *Num_Snd_Grid_ABC2ABD;
 
/*******************************************************
 int *Num_Rcv_Grid_AB2CA

  Num_Rcv_Grid_AB2CA gives the number of grid points
  received in the CA partition and sent from the AB 
  partition.
  size: Num_Rcv_Grid_AB2CA[numprocs]
*******************************************************/
int *Num_Rcv_Grid_AB2CA;
int *Num_Rcv_Grid_ABC2ABD;

/*******************************************************
 int *Num_Snd_Grid_CA2CB

  Num_Snd_Grid_CA2CB gives the number of grid points
  sent from the CA partition to the CB partition.
  size: Num_Snd_Grid_CA2CB[numprocs]
*******************************************************/
int *Num_Snd_Grid_CA2CB;
int *Num_Snd_Grid_ABD2DCA;
int *Num_Snd_Grid_DCA2DCB;

/*******************************************************
 int *Num_Rcv_Grid_CA2CB

  Num_Rcv_Grid_CA2CB gives the number of grid points
  received in the CB partition and sent from the CA
  partition.
  size: Num_Rcv_Grid_CA2CB[numprocs]
*******************************************************/
int *Num_Rcv_Grid_CA2CB;
int *Num_Rcv_Grid_ABD2DCA;
int *Num_Rcv_Grid_DCA2DCB;

/*******************************************************
 int **Index_Snd_Grid_AB2CA

  Index_Snd_Grid_AB2CA gives index, BN_AB in the partition 
  AB associated with the grid points sent to ID.
  size: Index_Snd_Grid_AB2CA[numprocs][Num_Snd_Grid_AB2CA[ID]]
*******************************************************/
int **Index_Snd_Grid_AB2CA;
int **Index_Snd_Grid_ABC2ABD;

/*******************************************************
 int **Index_Rcv_Grid_AB2CA

  Index_Rcv_Grid_AB2CA gives index, BN_AB in the partition 
  AB associated with the grid point received from ID.
  size: Index_Rcv_Grid_AB2CA[numprocs][Num_Rcv_Grid_AB2CA[ID]]
*******************************************************/
int **Index_Rcv_Grid_AB2CA;
int **Index_Rcv_Grid_ABC2ABD;

/*******************************************************
 int **Index_Snd_Grid_CA2CB

  Index_Snd_Grid_CA2CB gives index, BN_CA in the partition 
  CA associated with the grid points sent to ID.
  size: Index_Snd_Grid_CA2CB[numprocs][Num_Snd_Grid_CA2CB[ID]]
*******************************************************/
int **Index_Snd_Grid_CA2CB;
int **Index_Snd_Grid_ABD2DCA;
int **Index_Snd_Grid_DCA2DCB;

/*******************************************************
 int **Index_Rcv_Grid_CA2CB

  Index_Rcv_Grid_CA2CB gives index, BN_CA in the partition 
  CA associated with the grid points received from ID.
  size: Index_Rcv_Grid_CA2CB[numprocs][Num_Rcv_Grid_CA2CB[ID]]
*******************************************************/
int **Index_Rcv_Grid_CA2CB;
int **Index_Rcv_Grid_ABD2DCA;
int **Index_Rcv_Grid_DCA2DCB;

/*******************************************************
  int *ID_NN_AB2CA_S;
  int *ID_NN_AB2CA_R;

  global process ID used for sending and receiving data
  in MPI commucation (AB to CA).
*******************************************************/
int *ID_NN_AB2CA_S;
int *ID_NN_AB2CA_R;
int *ID_NN_ABC2ABD_S;
int *ID_NN_ABC2ABD_R;

/*******************************************************
  int *GP_AB2CA_S;
  int *GP_AB2CA_R;

  starting index to data used for sending and receiving data
  in MPI commucation (AB to CA).

*******************************************************/
int *GP_AB2CA_S;
int *GP_AB2CA_R;
int *GP_AB2CA_S_A;
int *GP_AB2CA_R_A;
int *GP_ABC2ABD_S;
int *GP_ABC2ABD_R;
int *GP_ABC2ABD_S_A;
int *GP_ABC2ABD_R_A;

/*******************************************************
  int *ID_NN_CA2CB_S;
  int *ID_NN_CA2CB_R;

  global process ID used for sending and receiving data
  in MPI commucation (CA to CB).

*******************************************************/
int *ID_NN_CA2CB_S;
int *ID_NN_CA2CB_R;
int *ID_NN_ABD2DCA_S;
int *ID_NN_ABD2DCA_R;
int *ID_NN_DCA2DCB_S;
int *ID_NN_DCA2DCB_R;

/*******************************************************
  int *GP_CA2CB_S;
  int *GP_CA2CB_R;

  starting index to data used for sending and receiving data
  in MPI commucation (CA to CB).

*******************************************************/
int *GP_CA2CB_S;
int *GP_CA2CB_R;
int *GP_CA2CB_S_A;
int *GP_CA2CB_R_A;
int *GP_ABD2DCA_S;
int *GP_ABD2DCA_R;
int *GP_DCA2DCB_S;
int *GP_DCA2DCB_R;
int *GP_ABD2DCA_S_A;
int *GP_ABD2DCA_R_A;
int *GP_DCA2DCB_S_A;
int *GP_DCA2DCB_R_A;

int *count_snd,*count_rcv,*displ_snd,*displ_rcv; 
MPI_Status *stat_arr,*ST_AB2CA_S,*ST_AB2CA_R,*ST_CA2CB_S,*ST_CA2CB_R;
MPI_Status *ST_ABC2ABD_S,*ST_ABC2ABD_R,*ST_ABD2DCA_S,*ST_ABD2DCA_R;
MPI_Status *ST_DCA2DCB_S,*ST_DCA2DCB_R;
MPI_Request *request_arr,*RQ_AB2CA_S,*RQ_AB2CA_R,*RQ_CA2CB_S,*RQ_CA2CB_R;
MPI_Request *RQ_ABC2ABD_S,*RQ_ABC2ABD_R,*RQ_ABD2DCA_S,*RQ_ABD2DCA_R;
MPI_Request *RQ_DCA2DCB_S,*RQ_DCA2DCB_R;
MPI_Status *stat_send,*stat_a2a;
MPI_Status *stat_recv;
MPI_Request *request_send,*request_a2a;
MPI_Request *request_recv;

int *recvID; 

dcomplex *array0,*array1,*Rhot;

double OFFT_timers[10];

/*
  added by hagita@nda.ac.jp
*/
double *OFFT_rin;

fftw_complex *OFFT_in, *OFFT_out;
fftw_plan OFFT_p_A,OFFT_p_B,OFFT_p_C,OFFT_p_D;

fftw_complex *OFFT_fin_A, *OFFT_fout_A,*OFFT_fin_B, *OFFT_fout_B;
fftw_complex *OFFT_fin_C, *OFFT_fout_C,*OFFT_fin_D, *OFFT_fout_D;
fftw_plan OFFT_p_fA,OFFT_p_fB,OFFT_p_fC,OFFT_p_fD;

fftw_complex *OFFT_min_A, *OFFT_mout_A,*OFFT_min_B, *OFFT_mout_B;
fftw_complex *OFFT_min_C, *OFFT_mout_C,*OFFT_min_D, *OFFT_mout_D;
fftw_plan OFFT_p_mA,OFFT_p_mB,OFFT_p_mC,OFFT_p_mD;

fftw_complex **OFFT_tin,**OFFT_tout;
fftw_plan *OFFT_p_tA,*OFFT_p_tB,*OFFT_p_tC,*OFFT_p_tD;

fftw_complex *OFFT_min, *OFFT_mout;

/*
  added by hagita@nda.ac.jp
*/
int OFFT_Ngrid3f;

int OFFT_Ngrid1,OFFT_Ngrid2,OFFT_Ngrid3,OFFT_Ngrid4;
int OFFT_My_Max_NumGrid,Max_OneD_Grids;
int My_NumGrid_AB,My_NumGrid_CB,My_NumGrid_CA;
int My_NumGrid_ABC,My_NumGrid_ABD,My_NumGrid_DCA,My_NumGrid_DCB;
int OFFT_My_NumGrid_In,OFFT_My_NumGrid_Out;
int Max_Num_Snd_Grid_AB2CA,Max_Num_Rcv_Grid_AB2CA;
int Max_Num_Snd_Grid_CA2CB,Max_Num_Rcv_Grid_CA2CB;
int Max_Num_Snd_Grid_AB2CA_A2A,Max_Num_Rcv_Grid_AB2CA_A2A;
int Max_Num_Snd_Grid_CA2CB_A2A,Max_Num_Rcv_Grid_CA2CB_A2A;
int Max_Num_Snd_Grid_ABC2ABD,Max_Num_Rcv_Grid_ABC2ABD;
int Max_Num_Snd_Grid_ABD2DCA,Max_Num_Rcv_Grid_ABD2DCA;
int Max_Num_Snd_Grid_DCA2DCB,Max_Num_Rcv_Grid_DCA2DCB;
int Max_Num_Snd_Grid_ABC2ABD_A2A,Max_Num_Rcv_Grid_ABC2ABD_A2A;
int Max_Num_Snd_Grid_ABD2DCA_A2A,Max_Num_Rcv_Grid_ABD2DCA_A2A;
int Max_Num_Snd_Grid_DCA2DCB_A2A,Max_Num_Rcv_Grid_DCA2DCB_A2A;
int Max_Num_Snd_Grid_A2A,Max_Num_Rcv_Grid_A2A;
int NN_AB2CA_S,NN_AB2CA_R,NN_CA2CB_S,NN_CA2CB_R;
int NN_ABC2ABD_S,NN_ABC2ABD_R,NN_ABD2DCA_S,NN_ABD2DCA_R,NN_DCA2DCB_S,NN_DCA2DCB_R;
int OFFT_alloc_first,OFFT_alloc_first1,OFFT_alloc_first2;
int OFFT_measure_time,OFFT_print_memory;
int size_Index_Snd_Grid_AB2CA,size_Index_Rcv_Grid_AB2CA;
int size_Index_Snd_Grid_CA2CB,size_Index_Rcv_Grid_CA2CB;
int size_Index_Snd_Grid_ABC2ABD,size_Index_Rcv_Grid_ABC2ABD;
int size_Index_Snd_Grid_ABD2DCA,size_Index_Rcv_Grid_ABD2DCA;
int size_Index_Snd_Grid_DCA2DCB,size_Index_Rcv_Grid_DCA2DCB;

int OFFT_MEASURE,COMM_PATT;
