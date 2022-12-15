!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    This file is part of OpenFFT. 
!
!    OpenFFT is an open-source parallel package for 3-D FFTs
!    built on a communication-optimal domain decomposition method.
! 
!    Copyright (C) 2013-2015  Truong Vinh Truong Duy and Taisuke Ozaki, 
!                             The University of Tokyo.
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  timing_c2c_3d.f90:
!  
!   This program is used for benchmarking the performance of OpenFFT 
!   with timing and GLOPS results.
!   It can be executed with an arbitrary number of processes.
!   Time is measured by MPI_Wtime().
!   A numeric input parameter can be provided for specifying the size of 
!   the 3 dimensions. If no input parameter is provided, it will be executed 
!   with a default size of 128^3 data points.     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program timing_c2c_3d

  use MPI
  use, intrinsic :: iso_c_binding
  implicit none
  include 'openfft.fi'  
  
  integer(C_INT) :: N1, N2, N3  
  integer(C_INT), parameter :: NTEST = 10 
  
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: Rhor, Rhok, in1

  integer(C_INT), dimension(6) :: My_Index_In,My_Index_Out

  real(C_DOUBLE) :: nn,flops
  integer(C_INT) :: ierror, nproc, nrank
  integer(C_INT) :: i,m,offt_measure,measure_time,print_memory
  integer(C_INT) :: My_Max_NumGrid,My_NumGrid_In,My_NumGrid_Out 
  real(C_DOUBLE) :: t1, t2, t3 ,t4

  integer :: num_args, ix
  character(len=12), dimension(:), allocatable :: args
  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)

  num_args = command_argument_count()
  if(num_args.eq.0) then
     N1 = 128
     N2 = N1
     N3 = N1
  elseif(num_args.eq.1) then
     allocate(args(num_args))
     ix = 1
     call get_command_argument(ix,args(ix))
     read (args(ix),*) N1
     N2 = N1
     N3 = N1
  endif

! Select auto-tuning of communication 

  offt_measure = 0

! Set whether to use the timing and print memory functions of OpenFFT or not.
! Default=0 (not use) 

  measure_time = 0
  print_memory = 0
  
! Initialize OpenFFT  

  t3 = openfft_init_c2c_3d(%VAL(N1),%VAL(N2),%VAL(N3),&
       My_Max_NumGrid,My_NumGrid_In,My_Index_In,My_NumGrid_Out,My_Index_Out,&
       %VAL(offt_measure),%VAL(measure_time),%VAL(print_memory))

! Allocate local input and output arrays 

  allocate (Rhor(My_Max_NumGrid))
  allocate (Rhok(My_Max_NumGrid))
  allocate (in1(My_Max_NumGrid))

! Set local input 

  do i=1,My_NumGrid_In
     Rhor(i) = (0.1, 0.2)
  end do
  
  in1 = Rhor
  
  t2 = 0.0D0
  t4 = 0.0D0
  t3 = 0.0D0

! Repeat FFT transform for NTEST times 

  do m=1,NTEST
     t1 = MPI_WTIME()
     t3 = openfft_exec_c2c_3d(Rhor, Rhok)
     t2 = t2 + MPI_WTIME() - t1

! Re-set local input 

     Rhor = in1
  end do
  
! Average time of all processes */

  call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
       MPI_COMM_WORLD,ierror)

  t1 = t1 / real(NTEST)
  
! Calculate flops */

  nn = real((N1*N2))
  nn = real((nn*N3))
  nn = nn ** (1.0/3.0)
  flops = 5.0 * nn * log(nn) / log(2.0)
  flops = flops * 3.0 * nn**2  
  flops = flops / t1

! Finalize OpenFFT 

  t3 = openfft_finalize()

! Free arrays 

  deallocate(Rhor,Rhok,in1)

! Print results 

  if (nrank==0) then
     write(*,*) '=======OpenFFT complex2complex========'
     write(*,'(A,I4,A,I4,A,I4,A,I5)') 'N1=',N1,' N2=',N2,' N3=',N3,' numprocs=',nproc
     write(*,'(A,F10.6)') 'Time openfft_execute Total =  ', t1
     write(*,'(A,F10.6)') 'Gflops                     =  ', flops / real(1000**3)
  end if
  
  call MPI_FINALIZE(ierror)

end program timing_c2c_3d

