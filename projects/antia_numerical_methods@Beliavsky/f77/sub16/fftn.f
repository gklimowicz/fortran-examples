!	To calculate the discrete Fourier Transform using FFT algorithm in n dimensions
!
!	ND : (input) Number of dimensions
!	NN : (input) Integer array of length ND containing the number of points
!		along each dimension. NN(I) is the number of points along
!		the Ith dimension, which should be a power of 2
!	CG : (input/output) Complex array of length NN(1)*NN(2)*...*NN(ND)
!		containing the data points. The dimensions of CG in the calling
!		program must exactly match the number of points, e.g.
!		CG(NN(1),NN(2),...,NN(ND))
!		After execution it will contain the Fourier transform of CG
!	IFLG : (input) Flag to decide whether to calculate forward or inverse
!		transform. If IFLG.GE.0 then Fourier transform is calculated
!		IF IFLG<0 then inverse Fourier transform is calculated
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=631 implies that at-least one of NN(I) is not a power of 2,
!			in this case contents of CG will be destroyed but will not
!			contain the Fourier transform.
!
!	Required routines : None

      SUBROUTINE FFTN(ND,NN,CG,IFLG,IER)
      IMPLICIT COMPLEX*32(C)
      COMPLEX*32 CWF,CWJ
      REAL*16 PI,TH
      PARAMETER(PI=3.14159265358979323846264338327950288Q0)
      DIMENSION CG(*),NN(ND)

      NTOT=1
      DO 1000 I=1,ND
        NTOT=NTOT*NN(I)
1000  CONTINUE
      CI=(0.Q0,1.Q0)

      NPR1=1
      IF(IFLG.GE.0) THEN
!	Calculate the DFT
        IW=1
      ELSE
!	Calculate the inverse DFT
        IW=-1
      ENDIF

!	Loop over each dimension
      DO 5000 ID=1,ND
        N=NN(ID)
        NPR=NPR1
        NPR1=NPR*N

!	Loop for bit reversal
        J=1
        DO 2000 I=1,NPR1,NPR
          IF(J.GT.I) THEN
            DO 1600 I1=I,NTOT,NPR1
              DO 1600 I2=I1,I1+NPR-1
                J2=I2+J-I
                CT=CG(I2)
                CG(I2)=CG(J2)
                CG(J2)=CT
1600        CONTINUE

          ENDIF
          M=NPR1/2
1800      IF(M.GE.NPR.AND.J.GT.M) THEN
            J=J-M
            M=M/2
            GO TO 1800
          ENDIF
          J=J+M
2000    CONTINUE

        IER=0
        J0=1
        K0=N/2
        TH=PI/K0
        CWF=-1

!	Loop for FFT calculation
3000    CWJ=1
        DO 3600 JR=1,J0
          JR0=(JR-1)*NPR+1
          DO 3400 IR=JR0,NTOT,2*J0*NPR
            DO 3400 I=IR,IR+NPR-1
              I1=I+J0*NPR
              CT=CG(I1)*CWJ
              CG(I1)=CG(I)-CT
              CG(I)=CG(I)+CT
3400      CONTINUE
          CWJ=CWJ*CWF
3600    CONTINUE

        J0=2*J0
        K0=K0/2
        IF(J0.EQ.N) GO TO 5000
        IF(J0.GT.N.OR.K0.EQ.0) THEN
!	N is not a power of 2
          IER=631
          RETURN
        ENDIF

        CWF=EXP(IW*K0*TH*CI)
        GO TO 3000
5000  CONTINUE
      END
