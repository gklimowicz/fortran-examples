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
!  check_c2c_3d.f90:
!  
!   This program transforms input data values to output data values.
!   It can be executed with an arbitrary number of processes.
!   Its input and output should match the corresponding values in 
!   check_c2c_3d.dat. 
!   This program does not require any input parameter.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program check_c2c_3d

  use MPI
  use, intrinsic :: iso_c_binding
  implicit none
  include 'openfft.fi'  
  
  integer(C_INT),parameter  :: N1=2, N2=3, N3=4  
  
  complex(C_DOUBLE_COMPLEX), dimension(N1,N2,N3) :: Input,Output
  complex(C_DOUBLE_COMPLEX), dimension(N1,N2,N3) :: Out,Output_ref
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: Rhor, Rhok
  integer(C_INT), dimension(6) :: My_Index_In,My_Index_Out

  real(C_DOUBLE) :: t3,factor 
  integer(C_INT) :: ierror, nproc, nrank
  integer(C_INT) :: i,j,k,l
  integer(C_INT) :: offt_measure,measure_time,print_memory
  integer(C_INT) :: My_Max_NumGrid,My_NumGrid_In,My_NumGrid_Out 

  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)

! Set global input 

  Input(1,1,1) = (1.000,  0.000)
  Input(1,1,2) = (0.999, -0.040)
  Input(1,1,3) = (0.987, -0.159)
  Input(1,1,4) = (0.936, -0.352)
  Input(1,2,1) = (0.994, -0.111)
  Input(1,2,2) = (0.989, -0.151)
  Input(1,2,3) = (0.963, -0.268)
  Input(1,2,4) = (0.891, -0.454)
  Input(1,3,1) = (0.903, -0.430)
  Input(1,3,2) = (0.885, -0.466)
  Input(1,3,3) = (0.823, -0.568)
  Input(1,3,4) = (0.694, -0.720)
  Input(2,1,1) = (0.500,  0.500)
  Input(2,1,2) = (0.499,  0.040)
  Input(2,1,3) = (0.487,  0.159)
  Input(2,1,4) = (0.436,  0.352)
  Input(2,2,1) = (0.494,  0.111)
  Input(2,2,2) = (0.489,  0.151)
  Input(2,2,3) = (0.463,  0.268)
  Input(2,2,4) = (0.391,  0.454)
  Input(2,3,1) = (0.403,  0.430)
  Input(2,3,2) = (0.385,  0.466)
  Input(2,3,3) = (0.323,  0.568)
  Input(2,3,4) = (0.194,  0.720)


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

  allocate(Rhor(My_Max_NumGrid),Rhok(My_Max_NumGrid)) 
  
! Set local input 

  call MPI_Barrier(MPI_COMM_WORLD,ierror)

  if(My_NumGrid_In > 0) then
     write (*,100) "myid=",nrank,": Input in the ABC(XYZ) order with ",&
          My_NumGrid_In," grid points from (A=",My_Index_In(1),",B=",&
          My_Index_In(2),",C=",My_Index_In(3)," to (A=",My_Index_In(4),&
          ",B=",My_Index_In(5),",C=",My_Index_In(6),")"
  else
     write(*,101) "myid=",nrank,": Input in the ABC(XYZ) order with ",&
          My_NumGrid_In," grid points"
  end if
100 format(A,I2,A,I2,A,I1,A,I1,A,I1,A,I1,A,I1,A,I1,A)
101 format(A,I2,A,I2,A)

  Rhor(:) = dcmplx(0.d0,0.d0)
  Rhok(:) = dcmplx(0.d0,0.d0)

  l=1
  if (My_NumGrid_In > 0) then 
  if (My_Index_In(1)+1.eq.My_Index_In(4)+1) then
     i=My_Index_In(1)+1
     do j=My_Index_In(2)+1,My_Index_In(5)+1
     do k=My_Index_In(3)+1,My_Index_In(6)+1
        Rhor(l) = Input(i,j,k)
        l=l+1
     enddo
     enddo
  elseif(My_Index_In(1)+1.lt.My_Index_In(4)+1) then
     do i=My_Index_In(1)+1,My_Index_In(4)+1
     if (i.eq.My_Index_In(1)+1) then
        do j=My_Index_In(2)+1,N2
        do k=My_Index_In(3)+1,My_Index_In(6)+1
           Rhor(l) = Input(i,j,k)
           l=l+1
        enddo
        enddo
     elseif(My_Index_In(1)+1.lt.i .and. i.lt.My_Index_In(4)+1) then
        do j=1,N2
        do k=My_Index_In(3)+1,My_Index_In(6)+1
           Rhor(l) = Input(i,j,k)
           l=l+1
        enddo
        enddo
     elseif(My_Index_In(4)+1.eq.i) then
        do j=1,My_Index_In(5)+1
        do k=My_Index_In(3)+1,My_Index_In(6)+1
           Rhor(l) = Input(i,j,k)
           l=l+1
        enddo
        enddo
     endif
     enddo
  endif
  endif

! Print global input
                                                    
  if (nrank==0) then
     write (*,*) 'Input values'
     do i=1,N1
        write(*,10) i
        write(*,*) ''
        do j=1,N2
           write(*,20) (real(Input(i,j,k)),k=1,N3)
           write(*,21) (aimag(Input(i,j,k)),k=1,N3)
           write(*,*) ''
        end do
     end do
  end if

10 format(1x,'Input(i,j,k) for i =', I6)
30 format(1x,'Output(i,j,k) for i =', I6)
20 format(1x,'Real ', 4F10.3)
21 format(1x,'Imag ', 4F10.3)


! FFT transform 

  t3 = openfft_exec_c2c_3d(Rhor, Rhok)

! Get local output 

  call MPI_Barrier(MPI_COMM_WORLD,ierror)

  if(My_NumGrid_Out > 0) then
     write(*,100) "myid=",nrank,": Output in the CBA(ZYX) order with ",&
          My_NumGrid_Out," grid points from (C=",My_Index_Out(1),",B=",&
          My_Index_Out(2),",A=",My_Index_Out(3)," to (C=",My_Index_Out(4),&
          ",B=",My_Index_Out(5),",A=",My_Index_Out(6),")"
  else
     write(*,101) "myid=",nrank,": Output in the CBC(ZYZ) order with ",&
          My_NumGrid_Out," grid points"
  end if
  
  factor = sqrt(dble(N1*N2*N3))

  Output(:,:,:) = dcmplx(0.d0,0.d0)
  Out(:,:,:) = dcmplx(0.d0,0.d0)

  l=1
  if (My_NumGrid_Out > 0) then 
  if (My_Index_Out(1)+1.eq.My_Index_Out(4)+1) then
     i=My_Index_Out(1)+1
     do j=My_Index_Out(2)+1,My_Index_Out(5)+1
     do k=My_Index_Out(3)+1,My_Index_Out(6)+1
        Out(k,j,i)=Rhok(l)/factor
        l=l+1
     enddo
     enddo
  elseif(My_Index_Out(1)+1.lt.My_Index_Out(4)+1) then
     do i=My_Index_Out(1)+1,My_Index_Out(4)+1
     if (i.eq.My_Index_Out(1)+1) then
        do j=My_Index_Out(2)+1,N2
        do k=My_Index_Out(3)+1,My_Index_Out(6)+1
           Out(k,j,i)=Rhok(l)/factor
           l=l+1
        enddo
        enddo
     elseif(My_Index_Out(1)+1.lt.i .and. i.lt.My_Index_Out(4)+1) then
        do j=1,N2
        do k=My_Index_Out(3)+1,My_Index_Out(6)+1
           Out(k,j,i)=Rhok(l)/factor
           l=l+1
        enddo
        enddo
     elseif(My_Index_Out(4)+1.eq.i) then
        do j=1,My_Index_Out(5)+1
        do k=My_Index_Out(3)+1,My_Index_Out(6)+1
           Out(k,j,i)=Rhok(l)/factor
           l=l+1
        enddo
        enddo
     endif
     enddo
  endif
  endif

! Gather results from all processes 

  call MPI_Allreduce(Out,Output,N1*N2*N3,MPI_DOUBLE_COMPLEX,MPI_SUM, &
       MPI_COMM_WORLD,ierror)
 
! Print global output                                                           
  if (nrank==0) then
     write (*,*) 'Output values'
     do i=1,N1
        write(*,30) i
        write(*,*) ''
        do j=1,N2
           write(*,20) (real(Output(i,j,k)),k=1,N3)
           write(*,21) (aimag(Output(i,j,k)),k=1,N3)
           write(*,*) ''
        end do
     end do
  end if

! Error check 

  Output_ref(1,1,1) = ( 3.292, 0.102)
  Output_ref(1,1,2) = ( 0.051,-0.042)
  Output_ref(1,1,3) = ( 0.113, 0.102)
  Output_ref(1,1,4) = ( 0.051, 0.246)
  Output_ref(1,2,1) = ( 0.143,-0.086)
  Output_ref(1,2,2) = ( 0.016, 0.153)
  Output_ref(1,2,3) = (-0.024, 0.127)
  Output_ref(1,2,4) = (-0.050, 0.086)
  Output_ref(1,3,1) = ( 0.143, 0.290)
  Output_ref(1,3,2) = (-0.050, 0.118)
  Output_ref(1,3,3) = (-0.024, 0.077)
  Output_ref(1,3,4) = ( 0.016, 0.051)
  Output_ref(2,1,1) = ( 1.225,-1.620)
  Output_ref(2,1,2) = ( 0.355, 0.083)
  Output_ref(2,1,3) = ( 0.000, 0.162)
  Output_ref(2,1,4) = (-0.355, 0.083)
  Output_ref(2,2,1) = ( 0.424, 0.320)
  Output_ref(2,2,2) = ( 0.020,-0.115)
  Output_ref(2,2,3) = ( 0.013,-0.091)
  Output_ref(2,2,4) = (-0.007,-0.080)
  Output_ref(2,3,1) = (-0.424, 0.320)
  Output_ref(2,3,2) = ( 0.007,-0.080)
  Output_ref(2,3,3) = (-0.013,-0.091)
  Output_ref(2,3,4) = (-0.020,-0.115)

  l = 0
  if (nrank==0) then
     do i=1,N1
        do j=1,N2
           do k=1,N3
              if(abs(real(Output(i,j,k)) - real(Output_ref(i,j,k))) &
                   .gt. 0.001 .or. (abs(aimag(Output(i,j,k)) - &
                   aimag(Output_ref(i,j,k))) .gt. 0.001)) then
                 l = 1
                 write(*,102) i,j,k,Output(i,j,k)
              end if
           end do
        end do
     end do
     if (l==0) then
        write(*,*) 'Check done. All output elements are correct.'
     else
        write(*,*) 'Check done. Some output elements are incorrect.'
     end if
  end if

102 format('ERROR Output(',I2,',',I2,',',I2,') (',F6.3,',',F6.3,')')

! Free arrays 

  deallocate(Rhor,Rhok) 

! Finalize OpenFFT 

  t3 = openfft_finalize()

  call MPI_FINALIZE(ierror)

end program check_c2c_3d

