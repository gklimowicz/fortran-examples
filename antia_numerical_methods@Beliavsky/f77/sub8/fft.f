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
