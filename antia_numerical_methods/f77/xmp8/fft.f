!     PROGRAM TO CALCULATE FAST FOURIER TRANSFORM OF A COMPLEX DATA SET
!     SIMPLE FOURIER TRANSFORM IS CALCULATED BY DFT, WHICH SHOULD
!     WORK EVEN WHEN N IS NOT A POWER OF 2

      PROGRAM FOUR
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      PARAMETER(PI=3.14159265358979324D0)
      DIMENSION CG(4100),CF(4100),CF1(4100),F1(4100)

!     EXAMPLE 10.5

      F(Y)=SIN(2.*PI*A*Y)

51    FORMAT(2X,1P2D14.6,5X,2D14.6)
52    FORMAT('    IER =',I4,5X,'N =',I7,6X,'A =',1PD14.6/
     1       2X,'FOURIER TRANSFORM (USING FFT)',5X,
     1       'FOURIER TRANSFORM (USING DFT)')
53    FORMAT(/2X,'MAXIMUM DIFFERENCE :   WITH FFT =',1PD11.3,4X,
     1       'WITH DFT =',D11.3)

!     A IS THE COEFFICIENT IN F(Y)

100   PRINT *,'TYPE A, N=NO. OF PTS IN DATA SET  (QUITS WHEN N.LE.0)'
      READ *,A,N
      IF(N.LE.0) STOP

!     GENERATING THE INPUT DATA SET

      H=1.D0/N
      DO 1000 I=1,N
        X=(I-1)*H
        CG(I)=F(X)
        CF(I)=CG(I)
        F1(I)=CG(I)
1000  CONTINUE

      IFLG=1
      CALL DFT(N,CG,CF,IFLG,IER)
      CALL FFT(N,CG,IFLG,IER)
      WRITE(6,52) IER,N,A
      WRITE(6,51) (CG(I),CF(I),I=1,N)

!      Take inverse transform and compare with input data
      IFLG=-1
      CALL FFT(N,CG,IFLG,IER)
      CALL DFT(N,CF,CF1,IFLG,IER)

      DIF1=0.0
      DIF2=0.0
      DO 2000 I=1,N
!	The inverse transform should be divided by N before comparing
        DIF1=MAX(DIF1,ABS(CG(I)/N-F1(I)))
        DIF2=MAX(DIF2,ABS(CF1(I)/N-F1(I)))
2000  CONTINUE
      WRITE(6,53) DIF1,DIF2
      GO TO 100
      END

!     --------------------------------------------------------

!	To calculate the discrete Fourier Transform
!	Normal evaluation of the sum, valid for any number of points
!	Should be used when the number of points is small
!
!	N : (input) Number of points
!	CG : (input) Complex array of length N containing the data points
!	CF : (output) Complex array of length N which will contain the
!		Fourier transform of CG
!	IFLG : (input) Flag to decide whether to calculate forward or inverse
!		transform. If IFLG.GE.0 then Fourier transform is calculated
!		IF IFLG<0 then inverse Fourier transform is calculated
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=611 implies that N<2 and no calculations are done
!
!	Required routines : None
 
      SUBROUTINE DFT(N,CG,CF,IFLG,IER)
      IMPLICIT COMPLEX*16(C)
!	The following declarations may be retained even for single precision
      COMPLEX*16 CWF,CWJ
      REAL*8 PI,TH
      PARAMETER(PI=3.14159265358979324D0)
      DIMENSION CG(N),CF(N)
 
      IF(N.LT.2) THEN
        IER=611
        RETURN
      ENDIF
      IER=0
      IF(IFLG.GE.0) THEN
        IW=1
      ELSE
        IW=-1
      ENDIF
      CI=(0.D0,1.D0)
 
!	Loop for Fourier transform
      DO 3000 I=1,N
        CW=2.*CI*IW*(I-1)*PI/N
        CS=0.0
        DO 2000 J=1,N
          CS=CS+CG(J)*EXP(CW*(J-1))
2000    CONTINUE
        CF(I)=CS
3000  CONTINUE
      END

!     --------------------------------------------------------

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
