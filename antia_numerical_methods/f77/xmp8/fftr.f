!     PROGRAM TO CALCULATE FAST FOURIER TRANSFORM OF REAL DATA

      PROGRAM FOUR
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      PARAMETER(PI=3.14159265358979324D0)
      DIMENSION CG(4100),G(8200)
      EQUIVALENCE(CG,G)

!     EXAMPLE 10.5

      F(Y)=SIN(2.*PI*A*Y)

51    FORMAT(2X,1P2D14.6)
52    FORMAT('    IER =',I4,5X,'N =',I7,5X,'A =',1PD14.6/
     1       5X,'FOURIER TRANSFORM')

100   PRINT *,'TYPE A,N=NO. OF PTS IN DATA SET    (QUITS WHEN N.LE.0)'
      READ *,A,N
      IF(N.LE.0) STOP

!     GENERATING THE INPUT DATA SET

      H=1.D0/N
      DO 1000 I=1,N
        X=(I-1)*H
        G(I)=F(X)
1000  CONTINUE

      IFLG=1
      CALL FFTR(N,CG,IFLG,IER)
      WRITE(6,52) IER,N,A
      WRITE(6,51) (CG(I),I=1,N/2)

      GO TO 100
      END

!     --------------------------------------

!	To calculate the discrete Fourier Transform using FFT algorithm
!
!	N : (input) Number of points, which must be a power of 2
!	CG : (input/output) Complex array of length N containing the data points
!		After execution it will contain the Fourier transform of CG
!	IFLG : (input) Flag to decide whether to calculate forward or inverse
!		transform. If IFLG.GE.0 then Fourier transform is calculated
!		IF IFLG<0 then inverse Fourier transform is calculated
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=611 implies that N<2, no calculations are done
!		IER=631 implies that N is not a power of 2, in this case
!			contents of CG will be destroyed but will not
!			contain the Fourier transform.
!
!	Required routines : None
 
      SUBROUTINE FFT(N,CG,IFLG,IER)
      IMPLICIT COMPLEX*16(C)
!	Following declarations may be retained even for COMPLEX*8 calculations
      COMPLEX*16 CWF,CWJ
      REAL*8 PI,TH
      PARAMETER(PI=3.14159265358979324D0)
      DIMENSION CG(N)

      IF(N.LT.2) THEN
        IER=611
        RETURN
      ENDIF
      CI=(0.D0,1.D0)

!	Bit reversal
      J=1
      DO 2000 I=1,N
        IF(J.GT.I) THEN
!	exchange CG(I) with CG(J)
          CT=CG(I)
          CG(I)=CG(J)
          CG(J)=CT
        ENDIF
        M=N/2
1800    IF(M.GE.1.AND.J.GT.M) THEN
          J=J-M
          M=M/2
          GO TO 1800
        ENDIF
!	J-1 is the bit reverse of I
        J=J+M
2000  CONTINUE

      IER=0
      J0=1
      K0=N/2
      TH=PI/K0
      IF(IFLG.GE.0) THEN
!	For DFT
        IW=1
      ELSE
!	For Inverse DFT
        IW=-1
      ENDIF
      CWF=-1

!	Main loop for FFT executed Log_2(N) times
3000  CWJ=1
!	Inner loop over all elements
      DO 3600 JR=1,J0
        DO 3400 I=JR,N,2*J0
          I1=I+J0
          CT=CG(I1)*CWJ
          CG(I1)=CG(I)-CT
          CG(I)=CG(I)+CT
3400    CONTINUE
        CWJ=CWJ*CWF
3600  CONTINUE

      J0=2*J0
      K0=K0/2
      IF(J0.EQ.N) RETURN
      IF(J0.GT.N.OR.K0.EQ.0) THEN
!	N is not a power of 2
        IER=631
        RETURN
      ENDIF

      CWF=EXP(IW*K0*TH*CI)
      GO TO 3000
      END

!     --------------------------------------

!	To calculate the discrete Fourier Transform of real data using FFT algorithm
!
!	N : (input) Number of points, which must be a power of 2
!	CG : (input/output) Complex array of length N/2 containing the data points
!		After execution it will contain the Fourier transform of CG
!		In the calling program this array may be treated as a
!		real array of length N, though the Fourier transform
!		will be complex.
!	IFLG : (input) Flag to decide whether to calculate forward or inverse
!		transform. If IFLG.GE.0 then Fourier transform is calculated
!		IF IFLG<0 then inverse Fourier transform is calculated
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=611 implies that N<4, no calculations are done
!		IER=631 implies that N is not a power of 2, in this case
!			contents of CG will be destroyed but will not
!			contain the Fourier transform.
!
!	Required routines : FFT
!
      SUBROUTINE FFTR(N,CG,IFLG,IER)
      IMPLICIT COMPLEX*16(C)
      COMPLEX*16 CW,CWF
      REAL*8 PI,TH,GR,GI
      PARAMETER(PI=3.14159265358979324D0)
      DIMENSION CG(N/2)

      NN=N/2
      TH=PI/NN
      CI=(0.D0,1.D0)
      IF(IFLG.GE.0) THEN
        CW=EXP(CI*TH)
        CF=(0.D0,-0.5D0)
!	Calculate the DFT of complex array of length N/2
        CALL FFT(NN,CG,IFLG,IER)
        IF(IER.GT.0) RETURN
      ELSE
        CF=(0.0D0,0.5D0)
        CW=EXP(-CI*TH)
      ENDIF

!	Rearranging the DFT
      CWF=CW
      DO 2000 I=2,NN/2+1
        I1=NN+2-I
        C1=0.5*(CG(I)+CONJG(CG(I1)))+CF*CWF*(CG(I)-CONJG(CG(I1)))
        CG(I1)=0.5*(CG(I1)+CONJG(CG(I)))-CF*(CG(I1)-CONJG(CG(I)))/CWF
        CG(I)=C1
        CWF=CWF*CW
2000  CONTINUE

!	The end points
      GR=CG(1)
      GI=IMAG(CG(1))
      IF(IFLG.GE.0) THEN
        CG(1)=GR+GI+CI*(GR-GI)
      ELSE
        CG(1)=0.5*(GR+GI)+CI*0.5*(GR-GI)
!	Calculate the inverse DFT
        CALL FFT(NN,CG,IFLG,IER)
      ENDIF
      END
