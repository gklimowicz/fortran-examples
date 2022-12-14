      MODULE FFT144
!@sum  FFT144 calculates the Fast Fourier Transform
!@auth Gary Russell
      USE CONSTANT, only : twopi,rt2,rt3
      IMPLICIT NONE
      SAVE
      REAL*8, PARAMETER :: rt2p1=rt2+1., rt2m1=rt2-1., twort3=2.*rt3
!@param KM length of input array; to change KM => rewrite module !!!
      INTEGER, PARAMETER :: KM=144 !@param KM length of input array

      REAL*8, PARAMETER :: BYKM=1d0/KM   !@param BYKM  1/KM
      REAL*8, PARAMETER :: BYKMH=2d0/KM  !@param BYKMH 1/(KM/2)
      REAL*8, PARAMETER :: BYKM2=1d0/(2*KM)!@param BYKM2 1/(2*KM)
!@var C,S cos/sin evaluated on grid points
      REAL*8 :: C(0:KM),S(0:KM)
!@var CH,SH cos/sin evaluated on half points
      REAL*8 :: CH(KM/2-1),SH(KM/2-1)
C**** intermediate sums for FFT
      REAL*8 C48Q0(48),C48Q1(48),S48Q1(48),
     *  C16Q0(16),C16Q1(16),C16Q2(16),C16Q3(16),C16Q4(16),
     *            S16Q1(16),S16Q2(16),S16Q3(16),S16Q4(16),
     *  C8QN(8,0:9),S8QN(8,9),
     *  C40N(0:18),C41N(0:18),C42N(0:18),C43N(0:18),
     *  S40N(18),S41N(18),S42N(18),S43N(18),
     *  C20N(0:36),C21N(0:36),S20N(36),S21N(36)
      COMMON /FFTCOM/ C48Q0,C48Q1,S48Q1,C16Q0,C16Q1,C16Q2,C16Q3,C16Q4,
     *     S16Q1,S16Q2,S16Q3,S16Q4,C8QN,S8QN,C40N,C41N,C42N,C43N,S40N,
     *     S41N,S42N,S43N,C20N,C21N,S20N,S21N

      END MODULE FFT144
C****
      SUBROUTINE FFT0 (IM)
!@sum  FFT0 initializes sines and cosines used by FFT routines.
!@auth Gary Russell
      USE FFT144
      IMPLICIT NONE
      INTEGER IQ,N !@var IQ,N loop variables
      INTEGER, INTENT(IN) :: IM    !@var IM size of arrays (must=KM)
      IF(IM.NE.KM)  GO TO 100
      DO N=0,KM/4
         C(N) = COS(TWOPI*N/REAL(KM,KIND=8))
         S(KM/4   -N) =  C(N)
         C(KM/2   -N) = -C(N)
         S(KM/4   +N) =  C(N)
         C(KM/2   +N) = -C(N)
         S(3*KM/4 -N) = -C(N)
         C(KM     -N) =  C(N)
         S(3*KM/4 +N) = -C(N)
      END DO
      DO N=1,KM/2-1
         CH(N) = COS(TWOPI*N/REAL(2*KM,KIND=8))
         SH(N) = SIN(TWOPI*N/REAL(2*KM,KIND=8))
      END DO
      RETURN
  100 WRITE (6,*) ' This version of FFT is for ',KM,'. IM =',IM
      call stop_model('stopped in FFT144.f',255)
      END SUBROUTINE FFT0
C****
      SUBROUTINE FFT (F,A,B)
!@sum   FFT calculates fast fourier transform of an input array F
!@auth  Gary Russell
!@ver   1.0
!@calls DOCALC
C****
C**** FFT calculates a fast fourier transform of the input array F,
C**** producing the cosine and sine coefficients in the output
C**** arrays A and B.  F is dimensioned by 144 = KM, and A and B are
C**** dimensioned by 0:72.  Upon entering FFT, the total energy is:
C****   .5*sum(F(K)**2)
C**** with the sum being taken over all K from 1 to KM.  The Fourier
C**** coefficients are defined by:
C****   A(N)+i*B(N) = sum(F(K)*exp(-2*PI*i*N*K/KM))/KMH
C**** with the sum being taken over all K from 1 to KM.  KMH = KM
C**** when N = 0 or KM/2, and KMH = KM/2 otherwise.  In the
C**** program's notation, CPQN means:
C****   CPQN = sum(F(K)*cos(2*PI*N*K/KM))
C**** with the sum being taken over all K = Q mod(P).  SPQN has a
C**** similar definition, but cos is replaced by sin.  The notation
C**** A=10, B=11, etc. is used for P, Q and N.  The same total
C**** energy can be calculated from the spectral coefficients as:
C****   .5*sum(A(N)**2+B(N)**2)*KMH
C**** with the sum being taken over all wave numbers N from 0 to KM/2.
C****
      USE FFT144
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM)      !@var F input gridpoint array
      REAL*8, INTENT(OUT) :: A(0:KM/2) !@var A fourier coeffs. (cos)
      REAL*8, INTENT(OUT) :: B(0:KM/2) !@var B fourier coeffs. (sin)
      INTEGER IQ,N !@var IQ,N loop variables

      CALL DOCALC(F)

C**** Calculate final coefficients of fourier expansion
      A(0) = (C20N(0)+C21N(0))*BYKM
      B(0) = 0.
      DO N=1,KM/4
        A(N) = (C20N(N)+C21N(N))*BYKMH
        B(N) = (S20N(N)+S21N(N))*BYKMH
      END DO
      DO N=1,KM/4-1
        A(KM/2-N) = (C20N(N)-C21N(N))*BYKMH
        B(KM/2-N) = (S21N(N)-S20N(N))*BYKMH
      END DO
      A(KM/2) = (C20N(0)-C21N(0))*BYKM
      B(KM/2) = 0.
      RETURN
      END SUBROUTINE FFT
C****

      SUBROUTINE FFTI (A,B,F)
!@sum  FFTI performs an inverse fast fourier transform
!@auth Gary Russell
      USE FFT144
      IMPLICIT NONE
      REAL*8, INTENT(OUT) :: F(KM)     !@var F output gridpoint array
      REAL*8, INTENT(IN) :: A(0:KM/2)  !@var A fourier coeffs. (cos)
      REAL*8, INTENT(IN) :: B(0:KM/2)  !@var B fourier coeffs. (sin)
      INTEGER IQ,N !@var IQ,N loop variables
      INTEGER Q,Q3,Q9

C**** These equations have been divided by 36
      DO N=1,35
        S20N(N) =  B(N)-B(72-N)
        S21N(N) =  B(N)+B(72-N)
        C20N(N) =  A(N)+A(72-N)
        C21N(N) =  A(N)-A(72-N)
      END DO
      C20N(0) = (A(0)+A(72))*2.
      C21N(0) = (A(0)-A(72))*2.
      C21N(36) = 0.
      S20N(36) = 0.
      C20N(36) = A(36)*2.
      S21N(36) = B(36)*2.
C**** 1/36 * 2 = 1/18
      DO N=0,18
        C41N(N) =  C21N(N)+S21N(36-N)
        C43N(N) =  C21N(N)-S21N(36-N)
      END DO
      DO N=1,18
        S41N(N) =  S21N(N)+C21N(36-N)
        S43N(N) =  S21N(N)-C21N(36-N)
      END DO
      DO N=0,17
        C40N(N) =  C20N(N)+C20N(36-N)
        C42N(N) =  C20N(N)-C20N(36-N)
      END DO
      DO N=1,17
        S40N(N) =  S20N(N)-S20N(36-N)
        S42N(N) =  S20N(N)+S20N(36-N)
      END DO
      C40N(18) = C20N(18)*2.
      C42N(18) = 0.
      S42N(18) = S20N(18)*2.
      S40N(18) = 0.
C**** 1/18 * 2 = 1/9
      DO N=0,8
        C8QN(8,N) = C40N(N)+C40N(18-N)
        C8QN(4,N) = C40N(N)-C40N(18-N)
      END DO
      C8QN(8,9) = C40N(9)*2.
      C8QN(4,9) =  0.
      DO N=0,9
        C8QN(2,N) = C42N(N)+S42N(18-N)
        C8QN(6,N) = C42N(N)-S42N(18-N)
      END DO
      DO N=1,9
        S8QN(2,N) = S42N(N)+C42N(18-N)
        S8QN(6,N) = S42N(N)-C42N(18-N)
      END DO
      DO N=1,8
        S8QN(4,N) = S40N(N)+S40N(18-N)
        S8QN(8,N) = S40N(N)-S40N(18-N)
      END DO
      S8QN(4,9) = S40N(9)*2.
      S8QN(8,9) =  0.
      DO N=1,8
        S8QN(5,N) = (S41N(18-N)-C41N(18-N))*S(18)+S41N(N)
        S8QN(1,N) =-(S41N(18-N)-C41N(18-N))*S(18)+S41N(N)
      END DO
      DO N=0,8
        C8QN(1,N) = (S41N(18-N)+C41N(18-N))*S(18)+C41N(N)
        C8QN(5,N) =-(S41N(18-N)+C41N(18-N))*S(18)+C41N(N)
      END DO
      C8QN(1,0) =  S41N(18)*RT2+C41N(0)
      S8QN(5,9) = (S41N(9)*RT2P1-C41N(9))*S(18)
      S8QN(1,9) = (S41N(9)*RT2M1+C41N(9))*S(18)
      C8QN(1,9) = (C41N(9)*RT2P1+S41N(9))*S(18)
      C8QN(5,9) = (C41N(9)*RT2M1-S41N(9))*S(18)
      DO N=1,8
        S8QN(3,N) = (S43N(18-N)+C43N(18-N))*S(18)+S43N(N)
        S8QN(7,N) =-(S43N(18-N)+C43N(18-N))*S(18)+S43N(N)
      END DO
      DO N=0,8
        C8QN(3,N) = (S43N(18-N)-C43N(18-N))*S(18)+C43N(N)
        C8QN(7,N) =-(S43N(18-N)-C43N(18-N))*S(18)+C43N(N)
      END DO
      C8QN(3,0) =  S43N(18)*RT2 +C43N(0)
      C8QN(7,0) = -S43N(18)*RT2 +C43N(0)
      S8QN(3,9) = (S43N(9)*RT2P1+C43N(9))*S(18)
      S8QN(7,9) = (S43N(9)*RT2M1-C43N(9))*S(18)
      C8QN(3,9) = (C43N(9)*RT2M1+S43N(9))*S(18)
      C8QN(7,9) = (C43N(9)*RT2P1-S43N(9))*S(18)
C**** Multiply again by 2 ==> 2/9
      DO Q=1,7,2
        Q9 = Q*9
        C16Q0(Q)   = C8QN(Q,0)+S8QN(Q,9)*S(Q9)+C8QN(Q,9)*C(Q9)
        C16Q0(Q+8) = C8QN(Q,0)-S8QN(Q,9)*S(Q9)-C8QN(Q,9)*C(Q9)
        C16Q1(Q)   = C8QN(Q,1)+S8QN(Q,8)*S(Q9)+C8QN(Q,8)*C(Q9)
        C16Q1(Q+8) = C8QN(Q,1)-S8QN(Q,8)*S(Q9)-C8QN(Q,8)*C(Q9)
        S16Q1(Q)   = S8QN(Q,1)+C8QN(Q,8)*S(Q9)-S8QN(Q,8)*C(Q9)
        S16Q1(Q+8) = S8QN(Q,1)-C8QN(Q,8)*S(Q9)+S8QN(Q,8)*C(Q9)
        C16Q2(Q)   = C8QN(Q,2)+S8QN(Q,7)*S(Q9)+C8QN(Q,7)*C(Q9)
        C16Q2(Q+8) = C8QN(Q,2)-S8QN(Q,7)*S(Q9)-C8QN(Q,7)*C(Q9)
        S16Q2(Q)   = S8QN(Q,2)+C8QN(Q,7)*S(Q9)-S8QN(Q,7)*C(Q9)
        S16Q2(Q+8) = S8QN(Q,2)-C8QN(Q,7)*S(Q9)+S8QN(Q,7)*C(Q9)
        C16Q3(Q)   = C8QN(Q,3)+S8QN(Q,6)*S(Q9)+C8QN(Q,6)*C(Q9)
        C16Q3(Q+8) = C8QN(Q,3)-S8QN(Q,6)*S(Q9)-C8QN(Q,6)*C(Q9)
        S16Q3(Q)   = S8QN(Q,3)+C8QN(Q,6)*S(Q9)-S8QN(Q,6)*C(Q9)
        S16Q3(Q+8) = S8QN(Q,3)-C8QN(Q,6)*S(Q9)+S8QN(Q,6)*C(Q9)
        C16Q4(Q)   = C8QN(Q,4)+S8QN(Q,5)*S(Q9)+C8QN(Q,5)*C(Q9)
        C16Q4(Q+8) = C8QN(Q,4)-S8QN(Q,5)*S(Q9)-C8QN(Q,5)*C(Q9)
        S16Q4(Q)   = S8QN(Q,4)+C8QN(Q,5)*S(Q9)-S8QN(Q,5)*C(Q9)
        S16Q4(Q+8) = S8QN(Q,4)-C8QN(Q,5)*S(Q9)+S8QN(Q,5)*C(Q9)
      END DO
C     Q = 2
      C16Q0(2)  = C8QN(2,0)+(S8QN(2,9)+C8QN(2,9))*S(18)
      C16Q0(10) = C8QN(2,0)-(S8QN(2,9)+C8QN(2,9))*S(18)
      C16Q1(2)  = C8QN(2,1)+(S8QN(2,8)+C8QN(2,8))*S(18)
      C16Q1(10) = C8QN(2,1)-(S8QN(2,8)+C8QN(2,8))*S(18)
      S16Q1(2)  = S8QN(2,1)+(C8QN(2,8)-S8QN(2,8))*S(18)
      S16Q1(10) = S8QN(2,1)-(C8QN(2,8)-S8QN(2,8))*S(18)
      C16Q2(2)  = C8QN(2,2)+(S8QN(2,7)+C8QN(2,7))*S(18)
      C16Q2(10) = C8QN(2,2)-(S8QN(2,7)+C8QN(2,7))*S(18)
      S16Q2(2)  = S8QN(2,2)+(C8QN(2,7)-S8QN(2,7))*S(18)
      S16Q2(10) = S8QN(2,2)-(C8QN(2,7)-S8QN(2,7))*S(18)
      C16Q3(2)  = C8QN(2,3)+(S8QN(2,6)+C8QN(2,6))*S(18)
      C16Q3(10) = C8QN(2,3)-(S8QN(2,6)+C8QN(2,6))*S(18)
      S16Q3(2)  = S8QN(2,3)+(C8QN(2,6)-S8QN(2,6))*S(18)
      S16Q3(10) = S8QN(2,3)-(C8QN(2,6)-S8QN(2,6))*S(18)
      C16Q4(2)  = C8QN(2,4)+(S8QN(2,5)+C8QN(2,5))*S(18)
      C16Q4(10) = C8QN(2,4)-(S8QN(2,5)+C8QN(2,5))*S(18)
      S16Q4(2)  = S8QN(2,4)+(C8QN(2,5)-S8QN(2,5))*S(18)
      S16Q4(10) = S8QN(2,4)-(C8QN(2,5)-S8QN(2,5))*S(18)
C     Q = 6
      C16Q0(6)  = C8QN(6,0)+(S8QN(6,9)-C8QN(6,9))*S(18)
      C16Q0(14) = C8QN(6,0)-(S8QN(6,9)-C8QN(6,9))*S(18)
      C16Q1(6)  = C8QN(6,1)+(S8QN(6,8)-C8QN(6,8))*S(18)
      C16Q1(14) = C8QN(6,1)-(S8QN(6,8)-C8QN(6,8))*S(18)
      S16Q1(6)  = S8QN(6,1)+(C8QN(6,8)+S8QN(6,8))*S(18)
      S16Q1(14) = S8QN(6,1)-(C8QN(6,8)+S8QN(6,8))*S(18)
      C16Q2(6)  = C8QN(6,2)+(S8QN(6,7)-C8QN(6,7))*S(18)
      C16Q2(14) = C8QN(6,2)-(S8QN(6,7)-C8QN(6,7))*S(18)
      S16Q2(6)  = S8QN(6,2)+(C8QN(6,7)+S8QN(6,7))*S(18)
      S16Q2(14) = S8QN(6,2)-(C8QN(6,7)+S8QN(6,7))*S(18)
      C16Q3(6)  = C8QN(6,3)+(S8QN(6,6)-C8QN(6,6))*S(18)
      C16Q3(14) = C8QN(6,3)-(S8QN(6,6)-C8QN(6,6))*S(18)
      S16Q3(6)  = S8QN(6,3)+(C8QN(6,6)+S8QN(6,6))*S(18)
      S16Q3(14) = S8QN(6,3)-(C8QN(6,6)+S8QN(6,6))*S(18)
      C16Q4(6)  = C8QN(6,4)+(S8QN(6,5)-C8QN(6,5))*S(18)
      C16Q4(14) = C8QN(6,4)-(S8QN(6,5)-C8QN(6,5))*S(18)
      S16Q4(6)  = S8QN(6,4)+(C8QN(6,5)+S8QN(6,5))*S(18)
      S16Q4(14) = S8QN(6,4)-(C8QN(6,5)+S8QN(6,5))*S(18)
C     Q = 4
      C16Q0(4)  = C8QN(4,0)+S8QN(4,9)
      C16Q0(12) = C8QN(4,0)-S8QN(4,9)
      C16Q1(4)  = C8QN(4,1)+S8QN(4,8)
      C16Q1(12) = C8QN(4,1)-S8QN(4,8)
      S16Q1(4)  = S8QN(4,1)+C8QN(4,8)
      S16Q1(12) = S8QN(4,1)-C8QN(4,8)
      C16Q2(4)  = C8QN(4,2)+S8QN(4,7)
      C16Q2(12) = C8QN(4,2)-S8QN(4,7)
      S16Q2(4)  = S8QN(4,2)+C8QN(4,7)
      S16Q2(12) = S8QN(4,2)-C8QN(4,7)
      C16Q3(4)  = C8QN(4,3)+S8QN(4,6)
      C16Q3(12) = C8QN(4,3)-S8QN(4,6)
      S16Q3(4)  = S8QN(4,3)+C8QN(4,6)
      S16Q3(12) = S8QN(4,3)-C8QN(4,6)
      C16Q4(4)  = C8QN(4,4)+S8QN(4,5)
      C16Q4(12) = C8QN(4,4)-S8QN(4,5)
      S16Q4(4)  = S8QN(4,4)+C8QN(4,5)
      S16Q4(12) = S8QN(4,4)-C8QN(4,5)
C     Q = 8
      C16Q0(8)  = C8QN(8,0)-C8QN(8,9)
      C16Q0(16) = C8QN(8,0)+C8QN(8,9)
      C16Q1(8)  = C8QN(8,1)-C8QN(8,8)
      C16Q1(16) = C8QN(8,1)+C8QN(8,8)
      S16Q1(8)  = S8QN(8,1)+S8QN(8,8)
      S16Q1(16) = S8QN(8,1)-S8QN(8,8)
      C16Q2(8)  = C8QN(8,2)-C8QN(8,7)
      C16Q2(16) = C8QN(8,2)+C8QN(8,7)
      S16Q2(8)  = S8QN(8,2)+S8QN(8,7)
      S16Q2(16) = S8QN(8,2)-S8QN(8,7)
      C16Q3(8)  = C8QN(8,3)-C8QN(8,6)
      C16Q3(16) = C8QN(8,3)+C8QN(8,6)
      S16Q3(8)  = S8QN(8,3)+S8QN(8,6)
      S16Q3(16) = S8QN(8,3)-S8QN(8,6)
      C16Q4(8)  = C8QN(8,4)-C8QN(8,5)
      C16Q4(16) = C8QN(8,4)+C8QN(8,5)
      S16Q4(8)  = S8QN(8,4)+S8QN(8,5)
      S16Q4(16) = S8QN(8,4)-S8QN(8,5)
C**** Multiply by 3 ==> 2/3
      DO Q=1,16
        Q3 = Q*3
        C48Q0(Q+16) = (S(Q3+48)*S16Q3(Q)+C(Q3+48)*C16Q3(Q))*2.+C16Q0(Q)
        C48Q0(Q)    = (C(Q3+96)*S16Q3(Q)-S(Q3+96)*C16Q3(Q))*TWORT3
     *       + C48Q0(Q+16)
        C48Q0(Q+32) =(-C(Q3)*S16Q3(Q)+C16Q3(Q)*S(Q3))*TWORT3+C48Q0(Q+16)
        S48Q1(Q+16) = ( (C16Q2(Q)-C16Q4(Q))*S(Q3+48)
     *       -(S16Q2(Q)-S16Q4(Q))*C(Q3+48))+S16Q1(Q)
        S48Q1(Q)    = ( (C16Q2(Q)-C16Q4(Q))*C(Q3+96)
     *       +(S16Q2(Q)-S16Q4(Q))*S(Q3+96))*RT3+S48Q1(Q+16)
        S48Q1(Q+32) = (-(C16Q2(Q)-C16Q4(Q))*C(Q3)
     *       -(S16Q2(Q)-S16Q4(Q))*S(Q3))*RT3+S48Q1(Q+16)
        C48Q1(Q+16) = ( (C16Q2(Q)+C16Q4(Q))*C(Q3+48)
     *       +(S16Q2(Q)+S16Q4(Q))*S(Q3+48))+C16Q1(Q)
        C48Q1(Q)    = (-(C16Q2(Q)+C16Q4(Q))*S(Q3+96)
     *       +(S16Q2(Q)+S16Q4(Q))*C(Q3+96))*RT3+C48Q1(Q+16)
        C48Q1(Q+32) = ( (C16Q2(Q)+C16Q4(Q))*S(Q3)
     *       -(S16Q2(Q)+S16Q4(Q))*C(Q3))*RT3+C48Q1(Q+16)
      END DO
C**** Put back factor of 3/2
      DO Q=1,48
        F(Q+48) = (S48Q1(Q)*S(Q+48)+C48Q1(Q)*C(Q+48))+C48Q0(Q)*.5
        F(Q)    = (S48Q1(Q)*C(Q+96)-C48Q1(Q)*S(Q+96))*RT3+F(Q+48)
        F(Q+96) =-(S48Q1(Q)*C(Q)   -C48Q1(Q)*S(Q))   *RT3+F(Q+48)
      END DO
C****
      RETURN
      END SUBROUTINE FFTI

      SUBROUTINE FFTE (F,E)
!@sum   FFTE calcs. the spectral energy E from input gridpoint values F
!@auth  Gary Russell
!@ver   1.0
!@calls DOCALC
      USE FFT144
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM)       !@var F input gridpoint array
      REAL*8, INTENT(OUT) :: E(0:KM/2)  !@var E spectral energy
      INTEGER IQ,N !@var IQ,N loop variables

      CALL DOCALC(F)
C****
      E(0) =  (C20N(0)+C21N(0))*(C20N(0)+C21N(0))*BYKM2
      DO N=1,KM/4
        E(N) = ((C20N(N)+C21N(N))*(C20N(N)+C21N(N)) +
     *       (S20N(N)+S21N(N))*(S20N(N)+S21N(N)))*BYKM
      END DO
      DO N=1,KM/4-1
        E(KM/2-N) = ((C20N(N)-C21N(N))*(C20N(N)-C21N(N)) +
     *       (S21N(N)-S20N(N))*(S21N(N)-S20N(N)))*BYKM
      END DO
      E(KM/2) = (C20N(0)-C21N(0))*(C20N(N)-C21N(N))*BYKM2
C****
      RETURN
      END SUBROUTINE FFTE

      SUBROUTINE FFT2 (F,A,B)
!@sum   FFT2 calculates Fourier coefficients on secondary grid
!@auth  Gary Russell
!@ver   1.0
!@calls FFT
      USE FFT144
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM)      !@var input array
      REAL*8, INTENT(OUT) :: A(0:KM/2) !@var fourier coefficients (cos)
      REAL*8, INTENT(OUT) :: B(0:KM/2) !@var fourier coefficients (sin)
      REAL*8  ANEW  !@var  ANEW  dummy variable
      INTEGER N !@var N loop variable

      CALL FFT (F,A,B)
      DO N=1,KM/2-1
        ANEW  = A(N)*CH(N)-B(N)*SH(N)
        B(N)  = A(N)*SH(N)+B(N)*CH(N)
        A(N)  = ANEW
      END DO
      B(KM/2) = A(KM/2)
      A(KM/2) = 0.
      RETURN
      END SUBROUTINE FFT2
C****
      SUBROUTINE DOCALC(F)
!@sum  DOCALC calculate intermediate expressions for FFT
!@auth Gary Russell
      USE FFT144
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM) !@var F input grid point array
      INTEGER IQ,N !@var IQ,N loop variables
      INTEGER Q,Q3,Q9
      REAL*8 X,Y

C**** Calculate expressions summed by increments of 48
      DO Q=1,48
        C48Q0(Q) = F(Q)+F(Q+48)+F(Q+96)
        C48Q1(Q) = F(Q)*C(Q)+F(Q+48)*C(Q+48)+F(Q+96)*C(Q+96)
        S48Q1(Q) = F(Q)*S(Q)+F(Q+48)*S(Q+48)+F(Q+96)*S(Q+96)
      END DO
C**** Calculate expressions summed by increments of 16
      DO Q=1,16
        C16Q0(Q) = C48Q0(Q)+C48Q0(Q+16)+C48Q0(Q+32)
        C16Q1(Q) = C48Q1(Q)+C48Q1(Q+16)+C48Q1(Q+32)
        Q3 = Q*3
        X = C48Q1(Q)*C(Q3)+C48Q1(Q+16)*C(Q3+48)+C48Q1(Q+32)*C(Q3+96)
        Y = S48Q1(Q)*S(Q3)+S48Q1(Q+16)*S(Q3+48)+S48Q1(Q+32)*S(Q3+96)
        C16Q2(Q) = X+Y
        C16Q3(Q) =
     *       C48Q0(Q)*C(Q3)+C48Q0(Q+16)*C(Q3+48)+C48Q0(Q+32)*C(Q3+96)
        C16Q4(Q) = X-Y
        S16Q1(Q) = S48Q1(Q)+S48Q1(Q+16)+S48Q1(Q+32)
        X = C48Q1(Q)*S(Q3)+C48Q1(Q+16)*S(Q3+48)+C48Q1(Q+32)*S(Q3+96)
        Y = S48Q1(Q)*C(Q3)+S48Q1(Q+16)*C(Q3+48)+S48Q1(Q+32)*C(Q3+96)
        S16Q2(Q) = X-Y
        S16Q3(Q) =
     *       C48Q0(Q)*S(Q3)+C48Q0(Q+16)*S(Q3+48)+C48Q0(Q+32)*S(Q3+96)
        S16Q4(Q) = X+Y
      END DO
C**** Calculate expressions summed by increments of 8
      DO Q=1,8
        C8QN(Q,0) = C16Q0(Q)+C16Q0(Q+8)
        C8QN(Q,1) = C16Q1(Q)+C16Q1(Q+8)
        C8QN(Q,2) = C16Q2(Q)+C16Q2(Q+8)
        C8QN(Q,3) = C16Q3(Q)+C16Q3(Q+8)
        C8QN(Q,4) = C16Q4(Q)+C16Q4(Q+8)
        Q9 = Q*9
        C8QN(Q,5) = (C16Q4(Q)-C16Q4(Q+8))*C(Q9)
     *       +(S16Q4(Q)-S16Q4(Q+8))*S(Q9)
        C8QN(Q,6) = (C16Q3(Q)-C16Q3(Q+8))*C(Q9)
     *       +(S16Q3(Q)-S16Q3(Q+8))*S(Q9)
        C8QN(Q,7) = (C16Q2(Q)-C16Q2(Q+8))*C(Q9)
     *       +(S16Q2(Q)-S16Q2(Q+8))*S(Q9)
        C8QN(Q,8) = (C16Q1(Q)-C16Q1(Q+8))*C(Q9)
     *       +(S16Q1(Q)-S16Q1(Q+8))*S(Q9)
        C8QN(Q,9) = (C16Q0(Q)-C16Q0(Q+8))*C(Q9)
        S8QN(Q,1) =  S16Q1(Q)+S16Q1(Q+8)
        S8QN(Q,2) =  S16Q2(Q)+S16Q2(Q+8)
        S8QN(Q,3) =  S16Q3(Q)+S16Q3(Q+8)
        S8QN(Q,4) =  S16Q4(Q)+S16Q4(Q+8)
        S8QN(Q,5) = (C16Q4(Q)-C16Q4(Q+8))*S(Q9)
     *       -(S16Q4(Q)-S16Q4(Q+8))*C(Q9)
        S8QN(Q,6) = (C16Q3(Q)-C16Q3(Q+8))*S(Q9)
     *       -(S16Q3(Q)-S16Q3(Q+8))*C(Q9)
        S8QN(Q,7) = (C16Q2(Q)-C16Q2(Q+8))*S(Q9)
     *       -(S16Q2(Q)-S16Q2(Q+8))*C(Q9)
        S8QN(Q,8) = (C16Q1(Q)-C16Q1(Q+8))*S(Q9)
     *       -(S16Q1(Q)-S16Q1(Q+8))*C(Q9)
        S8QN(Q,9) = (C16Q0(Q)-C16Q0(Q+8))*S(Q9)
      END DO
C**** Calculate expressions summed by increments of 4
      DO N=0,9
        C40N(N) = C8QN(8,N)+C8QN(4,N)
        C41N(N) = C8QN(1,N)+C8QN(5,N)
        C42N(N) = C8QN(2,N)+C8QN(6,N)
        C43N(N) = C8QN(3,N)+C8QN(7,N)
      END DO
      DO N=1,9
        S40N(N) = S8QN(8,N)+S8QN(4,N)
        S41N(N) = S8QN(1,N)+S8QN(5,N)
        S42N(N) = S8QN(2,N)+S8QN(6,N)
        S43N(N) = S8QN(3,N)+S8QN(7,N)
      END DO
      DO N=1,8
        C40N(18-N) =  C8QN(8,N)-C8QN(4,N)
        C41N(18-N) = (C8QN(1,N)-C8QN(5,N)+S8QN(1,N)-S8QN(5,N))*C(18)
        C42N(18-N) =  S8QN(2,N)-S8QN(6,N)
        C43N(18-N) = (C8QN(7,N)-C8QN(3,N)-S8QN(7,N)+S8QN(3,N))*C(18)
        S40N(18-N) =  S8QN(4,N)-S8QN(8,N)
        S41N(18-N) = (C8QN(1,N)-C8QN(5,N)-S8QN(1,N)+S8QN(5,N))*C(18)
        S42N(18-N) =  C8QN(2,N)-C8QN(6,N)
        S43N(18-N) = (C8QN(3,N)-C8QN(7,N)+S8QN(3,N)-S8QN(7,N))*C(18)
      END DO
      C40N(18) =  C8QN(8,0)-C8QN(4,0)
      C41N(18) = (C8QN(1,0)-C8QN(5,0))*C(18)
      C42N(18) =  0.
      C43N(18) = (C8QN(7,0)-C8QN(3,0))*C(18)
      S40N(18) =  0.
      S41N(18) = (C8QN(1,0)-C8QN(5,0))*C(18)
      S42N(18) =  C8QN(2,0)-C8QN(6,0)
      S43N(18) = (C8QN(3,0)-C8QN(7,0))*C(18)
C**** Calculate expressions summed by increments of 2
      DO N=0,18
        C20N(N) = C40N(N)+C42N(N)
        C21N(N) = C41N(N)+C43N(N)
      END DO
      DO N=1,18
        S20N(N) = S40N(N)+S42N(N)
        S21N(N) = S41N(N)+S43N(N)
      END DO
      DO N=1,17
        C20N(36-N) = C40N(N)-C42N(N)
        C21N(36-N) = S41N(N)-S43N(N)
        S20N(36-N) = S42N(N)-S40N(N)
        S21N(36-N) = C41N(N)-C43N(N)
      END DO
      C20N(36) = C40N(0)-C42N(0)
      C21N(36) = 0.
      S20N(36) = 0.
      S21N(36) = C41N(0)-C43N(0)
C****
      RETURN
      END SUBROUTINE DOCALC
