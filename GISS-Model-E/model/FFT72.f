      MODULE FFT72
!@sum  FFT72 calculates the Fast Fourier Transform
!@auth Gary Russell
      USE CONSTANT, only : twopi,rt2,rt3
      IMPLICIT NONE
      SAVE
!@param KM length of input array; to change KM => rewrite module !!!
      INTEGER, PARAMETER :: KM=72

      REAL*8, PARAMETER :: BYKM=1d0/KM   !@param BYKM  1/KM
      REAL*8, PARAMETER :: BYKMH=2d0/KM  !@param BYKMH 1/(KM/2)
      REAL*8, PARAMETER :: BYKM2=1d0/(2*KM)!@param BYKM2 1/(2*KM)
!@var C,S cos/sin evaluated on grid points
      REAL*8 :: C(0:KM),S(0:KM)
!@var CH,SH cos/sin evaluated on half points
      REAL*8 :: CH(KM/2-1),SH(KM/2-1)

!@var C240,C241,S241,C8,S8  intermediate sums for FFT
!@var C41,C42,C43,C44,S41,S42,S43,S44 intermediate sums for FFT
!@var C21,C22,S21,S22  intermediate sums for FFT
      REAL*8 :: C240(24), C241(24), S241(24)
      REAL*8 :: C8(8,0:4),S8(8,4)
      REAL*8 :: C41(0:9),C42(0:9),C43(0:9),C44(0:9),
     *          S41(  9),S42(  9),S43(  9),S44(  9)
      REAL*8 :: C21(0:18),C22(0:18),S21(0:18),S22(0:18)
      COMMON /FFTCOM/ C240,C241,S241,C8,S8,C41,C42,C43,C44,
     *                S41,S42,S43,S44,C21,C22,S21,S22

      END MODULE FFT72
C****
      SUBROUTINE FFT0 (IM)
!@sum  FFT0 initializes sines and cosines used by FFT routines.
!@auth Gary Russell
      USE FFT72
      IMPLICIT NONE
      INTEGER N !@var N loop variables
      INTEGER, INTENT(IN) :: IM    !@var IM size of arrays (must=KM)
C
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
      call stop_model('stopped in FFT72.f',255)
      END SUBROUTINE FFT0
C****
      SUBROUTINE FFT (F,A,B)
!@sum   FFT calculates fast fourier transform of an input array F
!@auth  Gary Russell
!@ver   1.0
!@calls DOCALC
C****
C**** FFT calculates a fast fourier transform of the input array F,
C**** producing the cosine and sine coefficients in the output arrays
C**** A and B.  F is dimensioned by KM; A and B are dimensioned
C**** 0:KM/2.  Upon entering FFT, the total energy is:
C****   .5*sum(F(K)**2)
C**** with the sum being taken over all K from 1 to KM.  The
C**** Fourier coefficients are defined by:
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
      USE FFT72
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM)      !@var F input gridpoint array
      REAL*8, INTENT(OUT) :: A(0:KM/2) !@var A fourier coeffs. (cos)
      REAL*8, INTENT(OUT) :: B(0:KM/2) !@var B fourier coeffs. (sin)
      INTEGER N !@var N loop variables

      CALL DOCALC(F)

C**** Calculate final coefficients of fourier expansion
      A(0) = (C22(0)+C21(0))*BYKM
      B(0) = 0.
      DO N=1,KM/4
         A(N) = (C22(N)+C21(N))*BYKMH
         B(N) = (S22(N)+S21(N))*BYKMH
      END DO
      DO N=1,KM/4-1
         A(KM/2-N) = (C22(N)-C21(N))*BYKMH
         B(KM/2-N) = (S21(N)-S22(N))*BYKMH
      END DO
      A(KM/2) = (C22(0)-C21(0))*BYKM
      B(KM/2) = 0.
C****
      RETURN
      END SUBROUTINE FFT
C****
      SUBROUTINE FFTI (A,B,F)
!@sum  FFTI performs an inverse fast fourier transform
!@auth Gary Russell
      USE FFT72
      IMPLICIT NONE
      REAL*8, INTENT(OUT) :: F(KM)     !@var F output gridpoint array
      REAL*8, INTENT(IN) :: A(0:KM/2)  !@var A fourier coeffs. (cos)
      REAL*8, INTENT(IN) :: B(0:KM/2)  !@var B fourier coeffs. (sin)
      INTEGER IQ,N !@var IQ,N loop variables

C**** These have been divided by 18
      C22(0) = (A(0)+A(36))*2.
      C21(0) = (A(0)-A(36))*2.
      DO N=1,18
         C22(N) = A(N)+A(36-N)
         C21(N) = A(N)-A(36-N)
         S22(N) = B(N)-B(36-N)
         S21(N) = B(N)+B(36-N)
      END DO
C**** Now multiply by 2, so FACTOR = 2/18 = 9
      DO N=0,9
         C42(N)=C22(N)-C22(18-N)
         C44(N)=C22(N)+C22(18-N)
      END DO
      DO N=0,8
         C41(N)=C21(N)+S21(18-N)
         C43(N)=C21(N)-S21(18-N)
      END DO
      DO N=1,9
         S42(N)=S22(N)+S22(18-N)
         S44(N)=S22(N)-S22(18-N)
         S41(N)=S21(N)+C21(18-N)
         S43(N)=S21(N)-C21(18-N)
      END DO
C**** Multiply by 2 again, so FACTOR = 2*2/18 = 2/9
      DO N=0,4
         C8(8,N)=C44(N)+C44(9-N)
         C8(4,N)=C44(N)-C44(9-N)
         C8(2,N)=C42(N)+S42(9-N)
         C8(6,N)=C42(N)-S42(9-N)
      END DO
      DO N=1,4
         S8(2,N)=S42(N)+C42(9-N)
         S8(6,N)=S42(N)-C42(9-N)
         S8(8,N)=S44(N)-S44(9-N)
         S8(4,N)=S44(N)+S44(9-N)
         S8(5,N)=S41(N)+(S41(9-N)-C41(9-N))*S(9)
         S8(1,N)=S41(N)-(S41(9-N)-C41(9-N))*S(9)
         S8(3,N)=S43(N)+(S43(9-N)+C43(9-N))*S(9)
         S8(7,N)=S43(N)-(S43(9-N)+C43(9-N))*S(9)
         C8(5,N)=C41(N)-(S41(9-N)+C41(9-N))*S(9)
         C8(1,N)=C41(N)+(S41(9-N)+C41(9-N))*S(9)
         C8(3,N)=C43(N)+(S43(9-N)-C43(9-N))*S(9)
         C8(7,N)=C43(N)-(S43(9-N)-C43(9-N))*S(9)
      END DO
      C8(5,0)=C41(0)-S41(9)*RT2
      C8(1,0)=C41(0)+S41(9)*RT2
      C8(3,0)=C43(0)+S43(9)*RT2
      C8(7,0)=C43(0)-S43(9)*RT2
C**** Multiply by 3, so FACTOR = 3*2*2/18 = 3*2/9 = 2/3
      C240(24)= C8(8,0)+C8(8,3)*2.
      C240(16)= C8(8,0)-C8(8,3)-S8(8,3)*RT3
      C240(8) = C8(8,0)-C8(8,3)+S8(8,3)*RT3
      C240(18)= C8(2,0)-S8(2,3)*2.
      C240(2) = C8(2,0)+S8(2,3)+C8(2,3)*RT3
      C240(10)= C8(2,0)+S8(2,3)-C8(2,3)*RT3
      C240(12)= C8(4,0)-C8(4,3)*2.
      C240(4) = C8(4,0)+C8(4,3)+S8(4,3)*RT3
      C240(20)= C8(4,0)+C8(4,3)-S8(4,3)*RT3
      C240(6) = C8(6,0)+S8(6,3)*2.
      C240(14)= C8(6,0)-S8(6,3)-C8(6,3)*RT3
      C240(22)= C8(6,0)-S8(6,3)+C8(6,3)*RT3
      C240(9) = C8(1,0)-(C8(1,3)-S8(1,3))*RT2
      C240(1) = C8(1,0)+C8(1,3)*2.*S(15)+S8(1,3)*2.*S(3)
      C240(17)= C8(1,0)-C8(1,3)*2.*S(3)-S8(1,3)*2.*S(15)
      C240(3) = C8(3,0)+(C8(3,3)+S8(3,3))*RT2
      C240(11)= C8(3,0)-C8(3,3)*2.*S(15)+S8(3,3)*2.*S(3)
      C240(19)= C8(3,0)+C8(3,3)*2.*S(3)-S8(3,3)*2.*S(15)
      C240(21)= C8(5,0)+(C8(5,3)-S8(5,3))*RT2
      C240(13)= C8(5,0)-C8(5,3)*2.*S(15)-S8(5,3)*2.*S(3)
      C240(5) = C8(5,0)+C8(5,3)*2.*S(3)+S8(5,3)*2.*S(15)
      C240(15)= C8(7,0)-(C8(7,3)+S8(7,3))*RT2
      C240(23)= C8(7,0)+C8(7,3)*2.*S(15)-S8(7,3)*2.*S(3)
      C240(7) = C8(7,0)-C8(7,3)*2.*S(3)+S8(7,3)*2.*S(15)
      C241(24)=  C8(8,2)+C8(8,4)+C8(8,1)
      S241(24)= -S8(8,2)+S8(8,4)+S8(8,1)
      C241(8) =-(C8(8,2)+C8(8,4))*S(6)+(S8(8,2)+S8(8,4))*S(12)+C8(8,1)
      C241(16)=-(C8(8,2)+C8(8,4))*S(6)-(S8(8,2)+S8(8,4))*S(12)+C8(8,1)
      S241(8) = (C8(8,2)-C8(8,4))*S(12)+(S8(8,2)-S8(8,4))*S(6)+S8(8,1)
      S241(16)=-(C8(8,2)-C8(8,4))*S(12)+(S8(8,2)-S8(8,4))*S(6)+S8(8,1)
      C241(18)= -S8(2,2)-S8(2,4)+C8(2,1)
      S241(18)= -C8(2,2)+C8(2,4)+S8(2,1)
      C241(2) = (C8(2,2)+C8(2,4))*S(12)+(S8(2,2)+S8(2,4))*S(6)+C8(2,1)
      C241(10)=-(C8(2,2)+C8(2,4))*S(12)+(S8(2,2)+S8(2,4))*S(6)+C8(2,1)
      S241(2) = (C8(2,2)-C8(2,4))*S(6)-(S8(2,2)-S8(2,4))*S(12)+S8(2,1)
      S241(10)= (C8(2,2)-C8(2,4))*S(6)+(S8(2,2)-S8(2,4))*S(12)+S8(2,1)
      S241(12)=  S8(4,2)-S8(4,4)+S8(4,1)
      C241(12)= -C8(4,2)-C8(4,4)+C8(4,1)
      S241(4) = (C8(4,2)-C8(4,4))*S(12)-(S8(4,2)-S8(4,4))*S(6)+S8(4,1)
      S241(20)=-(C8(4,2)-C8(4,4))*S(12)-(S8(4,2)-S8(4,4))*S(6)+S8(4,1)
      C241(4) = (C8(4,2)+C8(4,4))*S(6)+(S8(4,2)+S8(4,4))*S(12)+C8(4,1)
      C241(20)= (C8(4,2)+C8(4,4))*S(6)-(S8(4,2)+S8(4,4))*S(12)+C8(4,1)
      C241(6) =  S8(6,2)+S8(6,4)+C8(6,1)
      S241(6) =  C8(6,2)-C8(6,4)+S8(6,1)
      C241(22)= (C8(6,2)+C8(6,4))*S(12)-(S8(6,2)+S8(6,4))*S(6)+C8(6,1)
      C241(14)=-(C8(6,2)+C8(6,4))*S(12)-(S8(6,2)+S8(6,4))*S(6)+C8(6,1)
      S241(22)=-(C8(6,2)-C8(6,4))*S(6)-(S8(6,2)-S8(6,4))*S(12)+S8(6,1)
      S241(14)=-(C8(6,2)-C8(6,4))*S(6)+(S8(6,2)-S8(6,4))*S(12)+S8(6,1)
      C241(9) =(-C8(1,2)-C8(1,4)+S8(1,2)+S8(1,4))*S(9)+C8(1,1)
      S241(9) = (C8(1,2)-C8(1,4)+S8(1,2)-S8(1,4))*S(9)+S8(1,1)
      C241(1) = (C8(1,2)+C8(1,4))*S(15)+(S8(1,2)+S8(1,4))*S(3)+C8(1,1)
      C241(17)=-(C8(1,2)+C8(1,4))*S(3)-(S8(1,2)+S8(1,4))*S(15)+C8(1,1)
      S241(1) = (C8(1,2)-C8(1,4))*S(3)-(S8(1,2)-S8(1,4))*S(15)+S8(1,1)
      S241(17)=-(C8(1,2)-C8(1,4))*S(15)+(S8(1,2)-S8(1,4))*S(3)+S8(1,1)
      C241(3) = (C8(3,2)+C8(3,4)+S8(3,2)+S8(3,4))*S(9)+C8(3,1)
      S241(3) = (C8(3,2)-C8(3,4)-S8(3,2)+S8(3,4))*S(9)+S8(3,1)
      C241(11)=-(C8(3,2)+C8(3,4))*S(15)+(S8(3,2)+S8(3,4))*S(3)+C8(3,1)
      S241(19)=-(C8(3,2)-C8(3,4))*S(15)-(S8(3,2)-S8(3,4))*S(3)+S8(3,1)
      C241(19)= (C8(3,2)+C8(3,4))*S(3)-(S8(3,2)+S8(3,4))*S(15)+C8(3,1)
      S241(11)= (C8(3,2)-C8(3,4))*S(3)+(S8(3,2)-S8(3,4))*S(15)+S8(3,1)
      C241(21)= (C8(5,2)+C8(5,4)-S8(5,2)-S8(5,4))*S(9)+C8(5,1)
      S241(21)=(-C8(5,2)+C8(5,4)-S8(5,2)+S8(5,4))*S(9)+S8(5,1)
      C241(13)=-(C8(5,2)+C8(5,4))*S(15)-(S8(5,2)+S8(5,4))*S(3)+C8(5,1)
      S241(13)=-(C8(5,2)-C8(5,4))*S(3)+(S8(5,2)-S8(5,4))*S(15)+S8(5,1)
      C241(5) = (C8(5,2)+C8(5,4))*S(3)+(S8(5,2)+S8(5,4))*S(15)+C8(5,1)
      S241(5) = (C8(5,2)-C8(5,4))*S(15)-(S8(5,2)-S8(5,4))*S(3)+S8(5,1)
      C241(15)=-(C8(7,2)+C8(7,4)+S8(7,2)+S8(7,4))*S(9)+C8(7,1)
      S241(15)=(-C8(7,2)+C8(7,4)+S8(7,2)-S8(7,4))*S(9)+S8(7,1)
      C241(23)= (C8(7,2)+C8(7,4))*S(15)-(S8(7,2)+S8(7,4))*S(3)+C8(7,1)
      C241(7) =-(C8(7,2)+C8(7,4))*S(3)+(S8(7,2)+S8(7,4))*S(15)+C8(7,1)
      S241(23)=-(C8(7,2)-C8(7,4))*S(3)-(S8(7,2)-S8(7,4))*S(15)+S8(7,1)
      S241(7) = (C8(7,2)-C8(7,4))*S(15)+(S8(7,2)-S8(7,4))*S(3)+S8(7,1)
C**** Multiply by FACTOR = 3/2
      DO IQ=1,24
         F(IQ+24) = (S241(IQ)*S(IQ+24)+C241(IQ)*C(IQ+24))+C240(IQ)*.5
         F(IQ)    = (S241(IQ)*C(IQ+48)-C241(IQ)*S(IQ+48))*RT3+F(IQ+24)
         F(IQ+48) =-(S241(IQ)*C(IQ)   -C241(IQ)*S(IQ))   *RT3+F(IQ+24)
      END DO
C****
      RETURN
      END SUBROUTINE FFTI
C****
      SUBROUTINE FFTE (F,E)
!@sum   FFTE calcs. the spectral energy E from input gridpoint values F
!@auth  Gary Russell
!@ver   1.0
!@calls DOCALC
      USE FFT72
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM)       !@var F input gridpoint array
      REAL*8, INTENT(OUT) :: E(0:KM/2)  !@var E spectral energy
      INTEGER N !@var N loop variables

      CALL DOCALC(F)
C****
      E(0) =  (C22(0)+C21(0))*(C22(0)+C21(0))*BYKM2
      DO N=1,KM/4
         E(N) = ((C22(N)+C21(N))*(C22(N)+C21(N)) +
     *        (S22(N)+S21(N))*(S22(N)+S21(N)))*BYKM
      END DO
      DO N=1,KM/4-1
         E(KM/2-N) = ((C22(N)-C21(N))*(C22(N)-C21(N)) +
     *        (S21(N)-S22(N))*(S21(N)-S22(N)))*BYKM
      END DO
      E(KM/2) = (C22(0)-C21(0))*(C22(0)-C21(0))*BYKM2
C****
      RETURN
      END SUBROUTINE FFTE
C****
      SUBROUTINE FFT2 (F,A,B)
!@sum   FFT2 calculates Fourier coefficients on secondary grid
!@auth  Gary Russell
!@ver   1.0
!@calls FFT
C****
C**** Values in F(K) are located a half space right of those in FFT.
C****
      USE FFT72
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM)      !@var input array
      REAL*8, INTENT(OUT) :: A(0:KM/2) !@var fourier coefficients (cos)
      REAL*8, INTENT(OUT) :: B(0:KM/2) !@var fourier coefficients (sin)
      REAL*8  ANEW  !@var  ANEW  dummy variable
      INTEGER N !@var N loop variable

      CALL FFT (F,A,B)
      DO N=1,KM/2-1
         ANEW = A(N)*CH(N)-B(N)*SH(N)
         B(N) = A(N)*SH(N)+B(N)*CH(N)
         A(N) = ANEW
      END DO
      B(KM/2) = A(KM/2)
      A(KM/2) = 0.
      RETURN
      END SUBROUTINE FFT2
C****
      SUBROUTINE DOCALC(F)
!@sum  DOCALC calculate intermediate expressions for FFT
!@auth Gary Russell
      USE FFT72
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM) !@var F input grid point array
      INTEGER IQ,N !@var IQ,N loop variables

C**** Calculate expressions summed by increments of 24
      DO IQ=1,24
         C240(IQ) = F(IQ) + F(IQ+24) + F(IQ+48)
         C241(IQ) = F(IQ)*C(IQ) + F(IQ+24)*C(IQ+24) + F(IQ+48)*C(IQ+48)
         S241(IQ) = F(IQ)*S(IQ) + F(IQ+24)*S(IQ+24) + F(IQ+48)*S(IQ+48)
      END DO
C**** Calculate expressions summed by increments of 8
      DO IQ = 1,8
         C8(IQ,0) = C240(IQ)+C240(IQ+8)+C240(IQ+16)
         C8(IQ,1) = C241(IQ)+C241(IQ+8)+C241(IQ+16)
         S8(IQ,1) = S241(IQ)+S241(IQ+8)+S241(IQ+16)
      END DO

      C8(1,2)= (C241(1)-S241(17))*S(15)+(S241(9)-C241(9))*S(9)+(S241(1)
     *     -C241(17))*S(3)
      C8(2,2)= (C241(2)-C241(10))*S(12)+(S241(2)+S241(10))*S(6)-S241(18)
      C8(3,2)= (C241(3)+S241(3))*S(9)-(C241(11)+S241(19))*S(15)
     *     +(C241(19)+S241(11))*S(3)
      C8(4,2)= (C241(4)+C241(20))*S(6)+(S241(4)-S241(20))*S(12)-C241(12)
      C8(5,2)= (C241(5)-S241(13))*S(3)+(C241(21)-S241(21))*S(9)
     *     -(C241(13)-S241(5))*S(15)
      C8(6,2)= (C241(22)-C241(14))*S(12)-(S241(22)+S241(14))*S(6)
     *     +S241(6)
      C8(7,2)= (C241(23)+S241(7))*S(15)-(S241(15)+C241(15))*S(9)
     *     -(S241(23)+C241(7))*S(3)
      C8(8,2)= C241(24)-(C241(8)+C241(16))*S(6)+(S241(8)-S241(16))*S(12)

      C8(1,3)= C240(1)*S(15)-C240(9)*S(9)-C240(17)*S(3)
      C8(2,3)= (C240(2)-C240(10))*S(12)
      C8(3,3)= C240(3)*S(9)-C240(11)*S(15)+C240(19)*S(3)
      C8(4,3)= (C240(4)+C240(20))*S(6)-C240(12)
      C8(5,3)= C240(5)*S(3)-C240(13)*S(15)+C240(21)*S(9)
      C8(6,3)= (C240(22)-C240(14))*S(12)
      C8(7,3)= C240(23)*S(15)-C240(15)*S(9)-C240(7)*S(3)
      C8(8,3)= C240(24)-(C240(8)+C240(16))*S(6)

      C8(1,4)= (C241(1)+S241(17))*S(15)-(S241(9)+C241(9))*S(9)-(S241(1)
     *     +C241(17))*S(3)
      C8(2,4)= (C241(2)-C241(10))*S(12)-(S241(2)+S241(10))*S(6)+S241(18)
      C8(3,4)= (C241(3)-S241(3))*S(9)-(C241(11)-S241(19))*S(15)
     *     +(C241(19)-S241(11))*S(3)
      C8(4,4)= (C241(4)+C241(20))*S(6)-(S241(4)-S241(20))*S(12)-C241(12)
      C8(5,4)= (C241(5)+S241(13))*S(3)+(C241(21)+S241(21))*S(9)
     *     -(C241(13)+S241(5))*S(15)
      C8(6,4)= (C241(22)-C241(14))*S(12)+(S241(22)+S241(14))*S(6)
     *     -S241(6)
      C8(7,4)= (C241(23)-S241(7))*S(15)+(S241(15)-C241(15))*S(9)
     *     +(S241(23)-C241(7))*S(3)
      C8(8,4)= C241(24)-(C241(8)+C241(16))*S(6)-(S241(8)-S241(16))*S(12)

      S8(1,2)= (C241(1)+S241(17))*S(3)+(C241(9)+S241(9))*S(9)-(C241(17)
     *     +S241(1))*S(15)
      S8(2,2)= (C241(10)+C241(2))*S(6)+(S241(10)-S241(2))*S(12)-C241(18)
      S8(3,2)= (C241(3)-S241(3))*S(9)+(S241(11)-C241(19))*S(15)
     *     +(C241(11)-S241(19))*S(3)
      S8(4,2)= (C241(4)-C241(20))*S(12)-(S241(4)+S241(20))*S(6)+S241(12)
      S8(5,2)= (C241(5)+S241(13))*S(15)-(S241(21)+C241(21))*S(9)
     *     -(C241(13)+S241(5))*S(3)
      S8(6,2)= C241(6)-(C241(14)+C241(22))*S(6)+(S241(14)-S241(22))
     *     *S(12)
      S8(7,2)= (C241(7)-S241(23))*S(15)+(S241(15)-C241(15))*S(9)
     *     -(C241(23)-S241(7))*S(3)
      S8(8,2)= (C241(8)-C241(16))*S(12)+(S241(8)+S241(16))*S(6)-S241(24)

      S8(1,3)= C240(1)*S(3)+C240(9)*S(9)-C240(17)*S(15)
      S8(2,3)= (C240(10)+C240(2))*S(6)-C240(18)
      S8(3,3)= C240(3)*S(9)-C240(19)*S(15)+C240(11)*S(3)
      S8(4,3)= (C240(4)-C240(20))*S(12)
      S8(5,3)= C240(5)*S(15)-C240(21)*S(9)-C240(13)*S(3)
      S8(6,3)= C240(6)-(C240(14)+C240(22))*S(6)
      S8(7,3)= C240(7)*S(15)-C240(15)*S(9)-C240(23)*S(3)
      S8(8,3)= (C240(8)-C240(16))*S(12)

      S8(1,4)= (C241(1)-S241(17))*S(3)+(C241(9)-S241(9))*S(9)-(C241(17)
     *     -S241(1))*S(15)
      S8(2,4)= (C241(10)+C241(2))*S(6)-(S241(10)-S241(2))*S(12)-C241(18)
      S8(3,4)= (C241(3)+S241(3))*S(9)-(S241(11)+C241(19))*S(15)
     *     +(C241(11)+S241(19))*S(3)
      S8(4,4)= (C241(4)-C241(20))*S(12)+(S241(4)+S241(20))*S(6)-S241(12)
      S8(5,4)= (C241(5)-S241(13))*S(15)+(S241(21)-C241(21))*S(9)
     *     -(C241(13)-S241(5))*S(3)
      S8(6,4)= C241(6)-(C241(14)+C241(22))*S(6)-(S241(14)-S241(22))
     *     *S(12)
      S8(7,4)= (C241(7)+S241(23))*S(15)-(S241(15)+C241(15))*S(9)
     *     -(C241(23)+S241(7))*S(3)
      S8(8,4)= (C241(8)-C241(16))*S(12)-(S241(8)+S241(16))*S(6)+S241(24)
C**** Calculate expressions summed by increments of 4
      DO N=0,4
         C41(N)=C8(1,N)+C8(5,N)
         C42(N)=C8(2,N)+C8(6,N)
         C43(N)=C8(3,N)+C8(7,N)
         C44(N)=C8(8,N)+C8(4,N)
      END DO
      DO N=1,4
         S41(N)=S8(1,N)+S8(5,N)
         S42(N)=S8(2,N)+S8(6,N)
         S43(N)=S8(3,N)+S8(7,N)
         S44(N)=S8(8,N)+S8(4,N)
         C41(9-N) =(C8(1,N)-C8(5,N)+S8(1,N)-S8(5,N))*C(9)
         C42(9-N) = S8(2,N)-S8(6,N)
         C43(9-N) =(C8(7,N)-C8(3,N)+S8(3,N)-S8(7,N))*C(9)
         C44(9-N) = C8(8,N)-C8(4,N)
         S41(9-N) =(C8(1,N)-C8(5,N)+S8(5,N)-S8(1,N))*C(9)
         S42(9-N) = C8(2,N)-C8(6,N)
         S43(9-N) =(C8(3,N)-C8(7,N)+S8(3,N)-S8(7,N))*C(9)
         S44(9-N) = S8(4,N)-S8(8,N)
      END DO
      C41(9) =(C8(1,0)-C8(5,0))*C(9)
      C42(9) = 0.
      C43(9) =(C8(7,0)-C8(3,0))*C(9)
      C44(9) = C8(8,0)-C8(4,0)
      S41(9) =(C8(1,0)-C8(5,0))*C(9)
      S42(9) = C8(2,0)-C8(6,0)
      S43(9) =(C8(3,0)-C8(7,0))*C(9)
      S44(9) = 0.
C**** Calculate expressions summed by increments of 2
      DO N=0,9
         C21(N) = C41(N)+C43(N)
         C22(N) = C44(N)+C42(N)
      END DO
      DO N=1,8
         S21(N) = S41(N)+S43(N)
         S22(N) = S44(N)+S42(N)
      END DO
      DO N=1,9
         C21(18-N) = S41(N)-S43(N)
         C22(18-N) = C44(N)-C42(N)
         S21(18-N) = C41(N)-C43(N)
         S22(18-N) = S42(N)-S44(N)
      END DO
      C21(18) = 0.
      C22(18) = C44(0)-C42(0)
      S21(18) = C41(0)-C43(0)
      S22(18) = 0.
C****
      RETURN
      END SUBROUTINE DOCALC

