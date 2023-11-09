      program r2cfft
!@sum  Testing FFTW and comparison with results from G. Russell's FFT package
!@+    Until Makefile is created, use /home/dgueyffi/fftwtests/comp to compile
!@auth Denis Gueyffier

      use fftw_com
      use FFT72
      implicit none
      real*8 :: dpi
      integer :: i,k
      
      real*8, allocatable :: A(:),B(:),Ar(:),Br(:),input(:),x(:),E(:)

      N=72

      allocate( A(0:N/2),B(0:N/2),Ar(0:N/2),Br(0:N/2),
     &     input(N),x(N),E(0:N/2))

      dpi = 4.0d0*atan(1.d0)


c***  Gary Russell samples from 1 (->2 pi/ N)  to N (->2 pi )
c***  FFTW uses more widespread convention : sampling 
c***  from  0 (->0 )  to N-1 (-> 2  pi (1 -1/ N) ) , e.g. same convention as Matlab

      do i = 1, N
          x(i) = 2.0d0*dpi*i/N
      enddo

      write(*,*) 'Sampling points, x='
      write(*,*) x

c***  Initialize FFTW
      call fft0(N)

c***  Define input function, sampled using Garry Russel convention
      input(:) = 0.d0

c***  Input : input = ( 1, 0,  ..., 0) dirac
c***  FFT :     (1 , 1, 1, ..., 1)
c      input(N)=N   !dirac function
 
c***  Input : input = cos(m*x)  one wavelength 
c***  FFT :     (0,..0, 1, 0..,0)
      input = sin( 10 * x )   ! mode 10

c***  Input : input = cste function = 1
c***  FFT :     (1,..0, 0, 0) dirac function 
c      input = 1.d0

c***  DFT transform using FFTW
      call fft(input,A,B)

c***  compare with Gary Russell's FFT functions
      call FFT0_G(N)
      call FFT_G(input,Ar,Br)

      write(*,*) "A=",A(0:N/2)
      write(*,*) "B=",B(0:N/2)

      write(*,*) "L2 norm error on A=",sqrt(sum((Ar(0:N/2)-A)
     &     *(Ar(0:N/2)-A),1))
      write(*,*) "L2 norm error on B=",sqrt(sum((Br(0:N/2)-B)
     &     *(Br(0:N/2)-B),1))


c***  compare Inverse

      call ffti(A,B,input)

      write(*,*) "FFTW inverse=",input

      call FFTI_G(Ar,Br,input)

      write(*,*) "Gar inverse=",input

c***  compare energy (Parseval) 
      call ffte(input,E)

      write(*,*) "energy stored in each mode FFTW=",E

      call ffte_G(input,E)

      write(*,*) "energy stored in each mode Gary=",E

      call fftend()

      end program r2cfft
c*
