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
      IMPLICIT COMPLEX*32(C)
!	The following declarations may be retained even for single precision
      COMPLEX*32 CWF,CWJ
      REAL*16 PI,TH
      PARAMETER(PI=3.14159265358979323846264338327950288Q0)
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
      CI=(0.Q0,1.Q0)
 
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
