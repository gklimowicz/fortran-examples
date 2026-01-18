!     PROGRAM TO CALCULATE FAST FOURIER TRANSFORM IN N DIMENSIONS

      PROGRAM FOUR
      IMPLICIT COMPLEX*32(C)
      IMPLICIT REAL*16(A,B,D-H,O-Z)
      PARAMETER(PI=3.14159265358979323846264338327950288Q0)
      DIMENSION CG(34000),NN(20),CF(34000)

!     FOURIER TRANSFORM IN THREE DIMENSIONS

      F(X,Y,Z)=SIN(2.*PI*A*Y)*COS(2.*PI*A*X)*COS(2.*PI*A*Z)

51    FORMAT(3I6,2X,1P4D14.6)
52    FORMAT('   IER =',I4,5X,'NO. OF PTS ALONG THE AXES =',(10I5))
53    FORMAT('   A =',1PD14.6/5X,'I',5X,'J',5X,'K',
     1       '   FOURIER TRANSFORM')
54    FORMAT(/'   MAXIMUM DIFFERENCE IN INVERSE TRANSFORM =',1PD10.2)

      ND=3

100   PRINT *,'TYPE NN(1)...NN(ND)=NO. OF PTS ALONG EACH AXES'
      PRINT *,'                 (QUITS WHEN NN(1).LE.0)'
      READ *,(NN(I),I=1,ND)
      IF(NN(1).LE.0) STOP
      PRINT *,'TYPE A = COEF. IN SIN AND COS TERM'
      READ *,A

!     GENERATING THE INPUT DATA SET USING F(X,Y,Z)

      DO 1000 I=1,NN(1)
        DO 1000 J=1,NN(2)
          DO 1000 K=1,NN(3)
            X=(I-1.Q0)/NN(1)
            Y=(J-1.Q0)/NN(2)
            Z=(K-1.Q0)/NN(3)

!     ORDER OF (I,J,K) ELEMENT IN LINEAR ARRAY
!     ALTERNATELY THE ARRAY CG CAN BE DIMENSIONED AS CG(NN(1),NN(2),NN(3))
!     AND SET CG(I,J,K)=F(X,Y,Z)

            I1=I+(J-1)*NN(1)+(K-1)*NN(1)*NN(2)
            CG(I1)=F(X,Y,Z)
            CF(I1)=CG(I1)
1000  CONTINUE

      IFLG=1
      CALL FFTN(ND,NN,CG,IFLG,IER)
      WRITE(6,52) IER,(NN(I),I=1,ND)

!     WRITE ONLY NONZERO COMPONENTS OF THE FOURIER TRANSFORM Gijk
!     THE FIRST THREE INTEGERS IN EACH LINE GIVE THE VALUES OF ijk

      WRITE(6,53) A
      DO 2000 I=1,NN(1)*NN(2)*NN(3)
        I1=MOD(I-1,NN(1))
        I3=(I-1)/(NN(1)*NN(2))
        I2=(I-1-I3*NN(1)*NN(2))/NN(1)
        IF(ABS(CG(I)).GT.1.Q-4) WRITE(6,51) I1,I2,I3,CG(I)
2000  CONTINUE

!	Take the inverse transform and compare with original function
      IFLG=-1
      CALL FFTN(ND,NN,CG,IFLG,IER)
      N=NN(1)*NN(2)*NN(3)
      DIF=0.0
      DO 3000 I=1,NN(1)*NN(2)*NN(3)
!	The inverse transform must be divided by N before comparing
        DIF=MAX(DIF,ABS(CG(I)/N-CF(I)))
3000  CONTINUE
      WRITE(6,54) DIF
      GO TO 100
      END

!     --------------------------------------

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
