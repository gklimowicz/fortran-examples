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
      IMPLICIT COMPLEX*32(C)
      COMPLEX*32 CW,CWF
      REAL*16 PI,TH,GR,GI
      PARAMETER(PI=3.14159265358979323846264338327950288Q0)
      DIMENSION CG(N/2)

      NN=N/2
      TH=PI/NN
      CI=(0.Q0,1.Q0)
      IF(IFLG.GE.0) THEN
        CW=EXP(CI*TH)
        CF=(0.Q0,-0.5Q0)
!	Calculate the DFT of complex array of length N/2
        CALL FFT(NN,CG,IFLG,IER)
        IF(IER.GT.0) RETURN
      ELSE
        CF=(0.0Q0,0.5Q0)
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
