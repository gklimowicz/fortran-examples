      MODULE FFT36
!@sum  FFT36 calculates the Fast Fourier Transform
!@auth Gary Russell
      USE CONSTANT, only : twopi,rt3
      IMPLICIT NONE
      SAVE
!@param KM length of input array; to change KM => rewrite module !!!
      INTEGER, PARAMETER :: KM=36

      REAL*8, PARAMETER :: BYKM=1d0/KM   !@param BYKM  1/KM
      REAL*8, PARAMETER :: BYKMH=2d0/KM  !@param BYKMH 1/(KM/2)
      REAL*8, PARAMETER :: BYKM2=1d0/(2*KM)!@param BYKM2 1/(2*KM)
!@var S sin evaluated on grid points
      REAL*8 :: S(KM)
!@var CH,SH cos/sin evaluated on half points
      REAL*8 :: CH(KM/2-1),SH(KM/2-1)

!@var C2,S2,C4,S4 etc. intermediate coefficients
      REAL*8 :: C200,C210,C201,C211,C202,C212,C203,C213,C204,C214,C205,
     *     C215,C206,C216,C207,C217,C208,C218,C209,C400,C410,C420,C430,
     *     S201,S211,S202,S212,S203,S213,S204,S214,S205,S215,S206,
     *     S216,S207,S217,S208,S218,S219,S401,S421,S402,S422,S403,S423,
     *     S404,S424,S411,S431,S412,S432,S413,S433,S414,S434,C401,C421,
     *     C402,C422,C403,C423,C404,C424,C411,C431,C412,C432,C413,C433,
     *     C414,C434,CC00,CC40,CC80,CC90,CC10,CC50,CC60,CC20,CCA0,CC30,
     *     CC70,CCB0,CC01,SC01,SC41,SC81,CC41,CC81,CC61,SC61,CC21,CCA1,
     *     SC21,SCA1,CC91,SC91,CC11,CC51,SC11,SC51,CC31,SC31,SC71,SCB1,
     *     CC71,CCB1
      COMMON /FFTCOM/ C200,C210,C201,C211,C202,C212,C203,C213,
     *     C204,C214,C205,C215,C206,C216,C207,C217,C208,C218,C209,C400,
     *     C410,C420,C430,S201,S211,S202,S212,S203,S213,S204,S214,S205,
     *     S215,S206,S216,S207,S217,S208,S218,S219,S401,S421,S402,S422,
     *     S403,S423,S404,S424,S411,S431,S412,S432,S413,S433,S414,S434,
     *     C401,C421,C402,C422,C403,C423,C404,C424,C411,C431,C412,C432,
     *     C413,C433,C414,C434,CC00,CC40,CC80,CC90,CC10,CC50,CC60,CC20,
     *     CCA0,CC30,CC70,CCB0,CC01,SC01,SC41,SC81,CC41,CC81,CC61,SC61,
     *     CC21,CCA1,SC21,SCA1,CC91,SC91,CC11,CC51,SC11,SC51,CC31,SC31,
     *     SC71,SCB1,CC71,CCB1

      END MODULE FFT36
C****
      SUBROUTINE FFT0 (IM)
!@sum  FFT0 initializes sines and cosines used by FFT routines.
!@auth Gary Russell
      USE FFT36
      IMPLICIT NONE
      INTEGER N !@var IQ,N loop variables
      INTEGER, INTENT(IN) :: IM    !@var IM size of arrays (must=KM)
C
      IF(IM.NE.KM)  GO TO 100
      DO N=1,KM/4
        S(N) = DSIN(TWOPI*N/REAL(KM,KIND=8))
      END DO
      DO N=1,KM/2-1
        CH(N) = COS(TWOPI*N/REAL(2*KM,KIND=8))
        SH(N) = SIN(TWOPI*N/REAL(2*KM,KIND=8))
      END DO
      RETURN
 100  WRITE (6,*) ' This version of FFT is for ',KM,'. IM =',IM
      call stop_model('stopped in FFT36.f',255)
      END SUBROUTINE FFT0

      SUBROUTINE FFT (F,A,B)
!@sum   FFT calculates fast fourier transform of an input array F
!@auth  Gary Russell
!@ver   1.0
!@calls DOCALC
C****
C**** FFT calculates a fast fourier transform of the input array
C**** F, producing the cosine and sine coefficients in the output
C**** arrays A and B.  F is dimensioned by 36 = KM; A and B are
C**** dimensioned 0:18.  Upon entering FFT, the total energy is:
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
      USE FFT36
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM)      !@var F input gridpoint array
      REAL*8, INTENT(OUT) :: A(0:KM/2) !@var A fourier coeffs. (cos)
      REAL*8, INTENT(OUT) :: B(0:KM/2) !@var B fourier coeffs. (sin)

      CALL DOCALC(F)

C**** Calculate final coefficients of Fourier expansion
      A(0)  = (C200+C210)*BYKM
      A(1)  = (C201+C211)*BYKMH
      A(2)  = (C202+C212)*BYKMH
      A(3)  = (C203+C213)*BYKMH
      A(4)  = (C204+C214)*BYKMH
      A(5)  = (C205+C215)*BYKMH
      A(6)  = (C206+C216)*BYKMH
      A(7)  = (C207+C217)*BYKMH
      A(8)  = (C208+C218)*BYKMH
      A(9)  = (C400-C420)*BYKMH
      A(10) = (C208-C218)*BYKMH
      A(11) = (C207-C217)*BYKMH
      A(12) = (C206-C216)*BYKMH
      A(13) = (C205-C215)*BYKMH
      A(14) = (C204-C214)*BYKMH
      A(15) = (C203-C213)*BYKMH
      A(16) = (C202-C212)*BYKMH
      A(17) = (C201-C211)*BYKMH
      A(18) = (C200-C210)*BYKM
      B(0)  = 0.
      B(1)  = (S201+S211)*BYKMH
      B(2)  = (S202+S212)*BYKMH
      B(3)  = (S203+S213)*BYKMH
      B(4)  = (S204+S214)*BYKMH
      B(5)  = (S205+S215)*BYKMH
      B(6)  = (S206+S216)*BYKMH
      B(7)  = (S207+S217)*BYKMH
      B(8)  = (S208+S218)*BYKMH
      B(9)  = (C410-C430)*BYKMH
      B(10) = (S218-S208)*BYKMH
      B(11) = (S217-S207)*BYKMH
      B(12) = (S216-S206)*BYKMH
      B(13) = (S215-S205)*BYKMH
      B(14) = (S214-S204)*BYKMH
      B(15) = (S213-S203)*BYKMH
      B(16) = (S212-S202)*BYKMH
      B(17) = (S211-S201)*BYKMH
      B(18) = 0.
C****
      RETURN
      END SUBROUTINE FFT
C****
      SUBROUTINE FFTI (A,B,F)
!@sum  FFTI performs an inverse fast fourier transform
!@auth Gary Russell
      USE FFT36
      IMPLICIT NONE
      REAL*8, INTENT(OUT) :: F(KM)     !@var F output gridpoint array
      REAL*8, INTENT(IN) :: A(0:KM/2)  !@var A fourier coeffs. (cos)
      REAL*8, INTENT(IN) :: B(0:KM/2)  !@var B fourier coeffs. (sin)

C     S200 = 0
      S201 = B(1)-B(17)
      S211 = B(1)+B(17)
      S202 = B(2)-B(16)
      S212 = B(2)+B(16)
      S203 = B(3)-B(15)
      S213 = B(3)+B(15)
      S204 = B(4)-B(14)
      S214 = B(4)+B(14)
      S205 = B(5)-B(13)
      S215 = B(5)+B(13)
      S206 = B(6)-B(12)
      S216 = B(6)+B(12)
      S207 = B(7)-B(11)
      S217 = B(7)+B(11)
      S208 = B(8)-B(10)
      S218 = B(8)+B(10)
C     S210 = 0
      S219 = B(9)*2.

      C200 =(A(0)+A(18))*2.
      C210 =(A(0)-A(18))*2.
      C201 = A(1)+A(17)
      C211 = A(1)-A(17)
      C202 = A(2)+A(16)
      C212 = A(2)-A(16)
      C203 = A(3)+A(15)
      C213 = A(3)-A(15)
      C204 = A(4)+A(14)
      C214 = A(4)-A(14)
      C205 = A(5)+A(13)
      C215 = A(5)-A(13)
      C206 = A(6)+A(12)
      C216 = A(6)-A(12)
      C207 = A(7)+A(11)
      C217 = A(7)-A(11)
      C208 = A(8)+A(10)
      C218 = A(8)-A(10)
      C209 = A(9)*2.

      S401 = S201-S208
      S421 = S201+S208
      S402 = S202-S207
      S422 = S202+S207
      S403 = S203-S206
      S423 = S203+S206
      S404 = S204-S205
      S424 = S204+S205
      S411 = S211+C218
      S431 = S211-C218
      S412 = S212+C217
      S432 = S212-C217
      S413 = S213+C216
      S433 = S213-C216
      S414 = S214+C215
      S434 = S214-C215

      C400 = C200+C209
      C420 = C200-C209
      C401 = C201+C208
      C421 = C201-C208
      C402 = C202+C207
      C422 = C202-C207
      C403 = C203+C206
      C423 = C203-C206
      C404 = C204+C205
      C424 = C204-C205
      C410 = C210+S219
      C430 = C210-S219
      C411 = C211+S218
      C431 = C211-S218
      C412 = C212+S217
      C432 = C212-S217
      C413 = C213+S216
      C433 = C213-S216
      C414 = C214+S215
      C434 = C214-S215
C Factor here = 9/2
      CC00 =  C400+C403*2.
      CC40 =  C400-C403+S403*rt3
      CC80 =  C400-C403-S403*rt3
      CC90 =  C410-S413*2.
      CC10 =  C410+S413+C413*rt3
      CC50 =  C410+S413-C413*rt3
      CC60 =  C420-C423*2.
      CC20 =  C420+C423+S423*rt3
      CCA0 =  C420+C423-S423*rt3
      CC30 =  C430+S433*2.
      CC70 =  C430-S433-C433*rt3
      CCB0 =  C430-S433+C433*rt3
      CC01 =   C402+C404+C401
      SC01 =  -S402+S404+S401
      SC41 =  (S402-S404)*S(3)+(C402-C404)*S(6)+S401
      SC81 =  (S402-S404)*S(3)-(C402-C404)*S(6)+S401
      CC41 = -(C402+C404)*S(3)+(S402+S404)*S(6)+C401
      CC81 = -(C402+C404)*S(3)-(S402+S404)*S(6)+C401
      CC61 = -(C424+C422)+C421
      SC61 = -(S424-S422)+S421
      CC21 =  (C424+C422)*S(3)+(S424+S422)*S(6)+C421
      CCA1 =  (C424+C422)*S(3)-(S424+S422)*S(6)+C421
      SC21 =  (S424-S422)*S(3)-(C424-C422)*S(6)+S421
      SCA1 =  (S424-S422)*S(3)+(C424-C422)*S(6)+S421
      CC91 =   C411-S412-S414
      SC91 =   S411-C412+C414
      CC11 =  (S412+S414)*S(3)+(C412+C414)*S(6)+C411
      CC51 =  (S412+S414)*S(3)-(C412+C414)*S(6)+C411
      SC11 =  (C412-C414)*S(3)-(S412-S414)*S(6)+S411
      SC51 =  (C412-C414)*S(3)+(S412-S414)*S(6)+S411
      CC31 =   S432+S434+C431
      SC31 =   C432-C434+S431
      SC71 = -(C432-C434)*S(3)+(S432-S434)*S(6)+S431
      SCB1 = -(C432-C434)*S(3)-(S432-S434)*S(6)+S431
      CC71 = -(S432+S434)*S(3)-(C432+C434)*S(6)+C431
      CCB1 = -(S432+S434)*S(3)+(C432+C434)*S(6)+C431

C FACTOR HERE =(1/3)*(9/2)=3/2 (KM/24)
      F(36) =  CC00*S(3)+CC01
      F(35) =  CCB0*S(3)-SCB1*S(1)+CCB1*S(8)
      F(34) =  CCA0*S(3)-SCA1*S(2)+CCA1*S(7)
      F(33) = (CC90-SC91)*S(3)+CC91*S(6)
      F(32) =  CC80*S(3)-SC81*S(4)+CC81*S(5)
      F(31) =  CC70*S(3)-SC71*S(5)+CC71*S(4)
      F(30) = (CC60+CC61)*S(3)-SC61*S(6)
      F(29) =  CC50*S(3)-SC51*S(7)+CC51*S(2)
      F(28) =  CC40*S(3)-SC41*S(8)+CC41*S(1)
      F(27) =  CC30*S(3)-SC31
      F(26) =  CC20*S(3)-SC21*S(8)-CC21*S(1)
      F(25) =  CC10*S(3)-SC11*S(7)-CC11*S(2)
      F(24) = (CC00-CC01)*S(3)-SC01*S(6)
      F(23) =  CCB0*S(3)-CCB1*S(4)-SCB1*S(5)
      F(22) =  CCA0*S(3)-CCA1*S(5)-SCA1*S(4)
      F(21) = (CC90-SC91)*S(3)-CC91*S(6)
      F(20) =  CC80*S(3)-CC81*S(7)-SC81*S(2)
      F(19) =  CC70*S(3)-CC71*S(8)-SC71*S(1)
      F(18) =  CC60*S(3)-CC61
      F(17) =  CC50*S(3)-CC51*S(8)+SC51*S(1)
      F(16) =  CC40*S(3)-CC41*S(7)+SC41*S(2)
      F(15) = (CC30+SC31)*S(3)-CC31*S(6)
      F(14) =  CC20*S(3)-CC21*S(5)+SC21*S(4)
      F(13) =  CC10*S(3)-CC11*S(4)+SC11*S(5)
      F(12) =  1.5*CC00-F(24)-F(36)
      F(11) =  1.5*CCB0-F(23)-F(35)
      F(10) =  1.5*CCA0-F(22)-F(34)
      F( 9) =  1.5*CC90-F(21)-F(33)
      F( 8) =  1.5*CC80-F(20)-F(32)
      F( 7) =  1.5*CC70-F(19)-F(31)
      F( 6) =  1.5*CC60-F(18)-F(30)
      F( 5) =  1.5*CC50-F(17)-F(29)
      F( 4) =  1.5*CC40-F(16)-F(28)
      F( 3) =  1.5*CC30-F(15)-F(27)
      F( 2) =  1.5*CC20-F(14)-F(26)
      F( 1) =  1.5*CC10-F(13)-F(25)
C****
      RETURN
      END SUBROUTINE FFTI
C****
      SUBROUTINE FFTE (F,E)
!@sum   FFTE calcs. the spectral energy E from input gridpoint values F
!@auth  Gary Russell
!@ver   1.0
!@calls DOCALC
      USE FFT36
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM)       !@var F input gridpoint array
      REAL*8, INTENT(OUT) :: E(0:KM/2)  !@var E spectral energy

      CALL DOCALC(F)
C****
      E(0)  =  (C200+C210)*(C200+C210)*BYKM2
      E(1)  = ((C201+C211)*(C201+C211)+(S201+S211)*(S201+S211))*BYKM
      E(2)  = ((C202+C212)*(C202+C212)+(S202+S212)*(S202+S212))*BYKM
      E(3)  = ((C203+C213)*(C203+C213)+(S203+S213)*(S203+S213))*BYKM
      E(4)  = ((C204+C214)*(C204+C214)+(S204+S214)*(S204+S214))*BYKM
      E(5)  = ((C205+C215)*(C205+C215)+(S205+S215)*(S205+S215))*BYKM
      E(6)  = ((C206+C216)*(C206+C216)+(S206+S216)*(S206+S216))*BYKM
      E(7)  = ((C207+C217)*(C207+C217)+(S207+S217)*(S207+S217))*BYKM
      E(8)  = ((C208+C218)*(C208+C218)+(S208+S218)*(S208+S218))*BYKM
      E(9)  = ((C400-C420)*(C400-C420)+(C410-C430)*(C410-C430))*BYKM
      E(10) = ((C208-C218)*(C208-C218)+(S218-S208)*(S218-S208))*BYKM
      E(11) = ((C207-C217)*(C207-C217)+(S217-S207)*(S217-S207))*BYKM
      E(12) = ((C206-C216)*(C206-C216)+(S216-S206)*(S216-S206))*BYKM
      E(13) = ((C205-C215)*(C205-C215)+(S215-S205)*(S215-S205))*BYKM
      E(14) = ((C204-C214)*(C204-C214)+(S214-S204)*(S214-S204))*BYKM
      E(15) = ((C203-C213)*(C203-C213)+(S213-S203)*(S213-S203))*BYKM
      E(16) = ((C202-C212)*(C202-C212)+(S212-S202)*(S212-S202))*BYKM
      E(17) = ((C201-C211)*(C201-C211)+(S211-S201)*(S211-S201))*BYKM
      E(18) =  (C200-C210)*(C200-C210)*BYKM2
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
      USE FFT36
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
      USE FFT36
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: F(KM) !@var F input grid point array

C**** Calculate expressions summed by increments of 12
      CC00 = F(36)+F(12)+F(24)
      CC10 = F(1)+F(13)+F(25)
      CC20 = F(2)+F(14)+F(26)
      CC30 = F(3)+F(15)+F(27)
      CC40 = F(4)+F(16)+F(28)
      CC50 = F(5)+F(17)+F(29)
      CC60 = F(6)+F(18)+F(30)
      CC70 = F(7)+F(19)+F(31)
      CC80 = F(8)+F(20)+F(32)
      CC90 = F(9)+F(21)+F(33)
      CCA0 = F(10)+F(22)+F(34)
      CCB0 = F(11)+F(23)+F(35)
      CC01 = F(36)-(F(12)+F(24))*S(3)
      CC11 = F(1)*S(8)-F(13)*S(4)-F(25)*S(2)
      CC21 = F(2)*S(7)-F(14)*S(5)-F(26)*S(1)
      CC31 = (F(3)-F(15))*S(6)
      CC41 = F(4)*S(5)-F(16)*S(7)+F(28)*S(1)
      CC51 = F(5)*S(4)-F(17)*S(8)+F(29)*S(2)
      CC61 = (F(6)+F(30))*S(3)-F(18)
      CC71 = F(7)*S(2)-F(19)*S(8)+F(31)*S(4)
      CC81 = F(8)*S(1)-F(20)*S(7)+F(32)*S(5)
      CC91 = (F(33)-F(21))*S(6)
      CCA1 = F(34)*S(7)-F(10)*S(1)-F(22)*S(5)
      CCB1 = F(35)*S(8)-F(11)*S(2)-F(23)*S(4)
      SC01 = (F(12)-F(24))*S(6)
      SC11 = F(1)*S(1)+F(13)*S(5)-F(25)*S(7)
      SC21 = F(2)*S(2)+F(14)*S(4)-F(26)*S(8)
      SC31 = (F(3)+F(15))*S(3)-F(27)
      SC41 = F(4)*S(4)+F(16)*S(2)-F(28)*S(8)
      SC51 = F(5)*S(5)+F(17)*S(1)-F(29)*S(7)
      SC61 = (F(6)-F(30))*S(6)
      SC71 = F(7)*S(7)-F(19)*S(1)-F(31)*S(5)
      SC81 = F(8)*S(8)-F(20)*S(2)-F(32)*S(4)
      SC91 = F(9)-(F(21)+F(33))*S(3)
      SCA1 = F(10)*S(8)-F(22)*S(4)-F(34)*S(2)
      SCB1 = F(11)*S(7)-F(23)*S(5)-F(35)*S(1)
C**** Calculate expressions summed by increments of 4
      C400 = CC00+CC40+CC80
      C410 = CC10+CC50+CC90
      C420 = CC20+CC60+CCA0
      C430 = CC30+CC70+CCB0
      C401 = CC01+CC41+CC81
      C411 = CC11+CC51+CC91
      C421 = CC21+CC61+CCA1
      C431 = CC31+CC71+CCB1
      C402 = (CC01-(CC41+CC81)*S(3))+(SC41-SC81)*S(6)
      C412 = (CC11-CC51)*S(6)+((SC11+SC51)*S(3)-SC91)
      C422 = ((CC21+CCA1)*S(3)-CC61)+(SC21-SCA1)*S(6)
      C432 = (CCB1-CC71)*S(6)+(SC31-(SC71+SCB1)*S(3))
      C403 = CC00-(CC40+CC80)*S(3)
      C413 = (CC10-CC50)*S(6)
      C423 = (CC20+CCA0)*S(3)-CC60
      C433 = (CCB0-CC70)*S(6)
      C404 = (CC01-(CC41+CC81)*S(3))-(SC41-SC81)*S(6)
      C414 = (CC11-CC51)*S(6)-((SC11+SC51)*S(3)-SC91)
      C424 = ((CC21+CCA1)*S(3)-CC61)-(SC21-SCA1)*S(6)
      C434 = (CCB1-CC71)*S(6)-(SC31-(SC71+SCB1)*S(3))
      S401 = SC01+SC41+SC81
      S411 = SC11+SC51+SC91
      S421 = SC21+SC61+SCA1
      S431 = SC31+SC71+SCB1
      S402 = ((CC41-CC81)*S(6))+((SC41+SC81)*S(3)-SC01)
      S412 = ((CC11+CC51)*S(3)-CC91)+((SC51-SC11)*S(6))
      S422 = ((CC21-CCA1)*S(6))+(SC61-(SC21+SCA1)*S(3))
      S432 = (CC31-(CC71+CCB1)*S(3))+((SC71-SCB1)*S(6))
      S403 = (CC40-CC80)*S(6)
      S413 = (CC10+CC50)*S(3)-CC90
      S423 = (CC20-CCA0)*S(6)
      S433 = CC30-(CC70+CCB0)*S(3)
      S404 = (CC41-CC81)*S(6)-((SC41+SC81)*S(3)-SC01)
      S414 = ((CC11+CC51)*S(3)-CC91)-(SC51-SC11)*S(6)
      S424 = (CC21-CCA1)*S(6)-(SC61-(SC21+SCA1)*S(3))
      S434 = (CC31-(CC71+CCB1)*S(3))-(SC71-SCB1)*S(6)
C**** Calculate expressions summed by increments of 2
      C200 = C400+C420
      C210 = C410+C430
      C201 = C401+C421
      C211 = C411+C431
      C202 = C402+C422
      C212 = C412+C432
      C203 = C403+C423
      C213 = C413+C433
      C204 = C404+C424
      C214 = C414+C434
      C205 = C404-C424
      C215 = S414-S434
      C206 = C403-C423
      C216 = S413-S433
      C207 = C402-C422
      C217 = S412-S432
      C208 = C401-C421
      C218 = S411-S431
C     C209 = C400-C420
C     C219 = 0.
C     S200 = 0.
C     S210 = 0.
      S201 = S401+S421
      S211 = S411+S431
      S202 = S402+S422
      S212 = S412+S432
      S203 = S403+S423
      S213 = S413+S433
      S204 = S404+S424
      S214 = S414+S434
      S205 = S424-S404
      S215 = C414-C434
      S206 = S423-S403
      S216 = C413-C433
      S207 = S422-S402
      S217 = C412-C432
      S208 = S421-S401
      S218 = C411-C431
C     S209 = 0.
C     S219 = C410-C430
C****
      RETURN
      END SUBROUTINE DOCALC
