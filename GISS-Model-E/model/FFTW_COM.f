      module fftw_com
!@sum This module provides variables and routines 
!@+   to compute FFTs using the external FFTW library
!@auth Denis Gueyffier

      implicit none
      include "fftw3.f"
      save
      integer*8 :: plan_forward,plan_backward
      integer :: Nsample
      real *8,allocatable :: infftw(:)
      double complex, allocatable ::  outfftw(:)

      end module fftw_com

      subroutine fft0(M)
      use fftw_com
      implicit none
      integer, intent(in) :: M

      Nsample=M
c***  Allocate input/output arrays
      allocate(infftw(Nsample),outfftw(Nsample/2+1))

c***  Create plan 
      call dfftw_plan_dft_r2c_1d(plan_forward,Nsample,infftw,
     &     outfftw,FFTW_ESTIMATE)

      call dfftw_plan_dft_c2r_1d(plan_backward,Nsample,outfftw,
     &     infftw,FFTW_ESTIMATE)



      end subroutine fft0
c*

      subroutine fftend
      use fftw_com
      implicit none

      call dfftw_destroy_plan(plan_forward)

      call dfftw_destroy_plan(plan_backward)

      deallocate(infftw,outfftw)

      end subroutine fftend
c*


      subroutine fft(F,A,B)
!@sum Computing FFT using FFTW
!@+   Input/Output is reorder to conform to G. Russell's convention
!@+   Gary Russell samples from 1 (->2 pi/ Nsample)  to Nsample(->2 pi )
!@+   FFTW uses a more widespread convention (same as Matlab for example): 
!@+   sampling from  0 (->0 )  to Nsample-1 (-> 2  pi (1 -1/ Nsample) ) 
      use fftw_com
      implicit none
      real*8, intent(in) :: F(1:Nsample)
      real*8, intent(out) :: A(0:Nsample/2),B(0:Nsample/2)          
      real*8 :: Atemp(1:Nsample/2+1),Btemp(1:Nsample/2+1)
      integer :: i, k

c***  Converting input from Gary Russell's convention to FFTW convention
      infftw(2:Nsample)=F(1:Nsample-1)/Nsample
      infftw(1)=F(Nsample)/Nsample

c***  Compute FFT
      call dfftw_execute_(plan_forward)

c      write(*,*) 'Afttw',real(outfftw)
c      write(*,*) 'Bfttw',imag(outfftw)

c***  Converting back output from FFTW convention to Gary Russell's convention 
      Atemp(:)=2.d0*real(outfftw)
      Atemp(1)=Atemp(1)/2.d0
      Atemp(Nsample/2+1)=Atemp(Nsample/2+1)/2.d0

      Btemp(:)=-2.d0*imag(outfftw)
      Btemp(1)=Btemp(1)/2.d0
      Btemp(Nsample/2+1)=Btemp(Nsample/2+1)/2.d0

      A(0:Nsample/2)=Atemp(1:Nsample/2+1)
      B(0:Nsample/2)=Btemp(1:Nsample/2+1)

      end subroutine fft
c*

      subroutine ffti(A,B,F)
!@sum inverse fft
!@auth Denis Gueyffier
      use fftw_com
      implicit none
      real*8, intent(out) :: F(1:Nsample)
      real*8, intent(in) :: A(0:Nsample/2),B(0:Nsample/2)  
      double complex :: cinput(0:Nsample/2)
      real*8 :: Atemp(1:Nsample/2+1),Btemp(1:Nsample/2+1)

c***  From Gary Russell's convention to FFTW convention
      Atemp(2:Nsample/2)=A(1:Nsample/2-1)/2.0
      Btemp(2:Nsample/2)=-B(1:Nsample/2-1)/2.0

      Atemp(1)=A(0)
      Btemp(1)=-B(0)

      Atemp(Nsample/2+1)=A(Nsample/2)
      Btemp(Nsample/2+1)=-B(Nsample/2)

      cinput=dcmplx(Atemp,Btemp)

      outfftw=cinput

c***  Compute InverseFFT
      call dfftw_execute_(plan_backward)
      
      F(Nsample)=infftw(1)
      F(1:Nsample-1)=infftw(2:Nsample)

      end subroutine ffti
c*


      subroutine FFTE(F,E)
!@sum  Spectral energy E from input gridpoint values F
!@auth Denis Gueyffier
      use fftw_com
      implicit none
      real*8, intent(in) :: F(1:Nsample)
      real*8, intent(out) :: E(0:Nsample/2)
      real*8 :: A(0:Nsample/2),B(0:Nsample/2)

      call fft(F,A,B)

      E(0)=0.5*(A(0)*A(0)+B(0)*B(0))*Nsample
      E(1:Nsample/2-1)=0.5*(A(1:Nsample/2-1)*A(1:Nsample/2-1)
     &   +B(1:Nsample/2-1)*B(1:Nsample/2-1))*Nsample/2.
      E(Nsample/2)=0.5*(A(Nsample/2)*A(Nsample/2)+B(Nsample/2)
     &   *B(Nsample/2))*Nsample


      end subroutine FFTE
c*

