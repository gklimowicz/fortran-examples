!	Multiple integration over a hyper-rectangle in n-dimensions
!	using product Gauss-Legendre formulas
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions
!	M : (input/output) Integer array of length N specifying the formula
!		to be used along each dimension. M(J)-point formula will
!		be used along Jth dimension, M(J) should be 2,4,8,16 or 32
!		otherwise it will be set to a default value of 2
!	IND : (input/output) Integer array of length N specifying the number
!		of subintervals to be used along each dimension. IND(J)>0
!		otherwise it will be set to a default value of 1
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	RINT : (output) The calculated value of the integral
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=39 implies specified accuracy was not achieved in
!			which case DIF will contain the estimated accuracy
!		IER=305 implies N<1 or N>20 and no calculations are done
!		IER=307 implies that number of points exceeded MAXPT in
!			first attempt and no approximation of RINT is calculated
!	NUM : (output) Number of function evaluations used by subroutine
!	MAXPT : (input/output) Maximum number of function evaluations permitted
!		If MAXPT <1 it is set to a default value of MAXPTS (=1100000)
!
!	FUNCTION F(N,X) must be supplied by the user
!	
!	Required routines : NGAUSS, F
!
      SUBROUTINE MULINT(A,B,N,M,IND,F,RINT,REPS,AEPS,DIF,IER,NUM,MAXPT)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL F
      PARAMETER(MAXPTS=1100000)
      DIMENSION A(N),B(N),M(N),IND(N)

!	Set M(I) and IND(I) to default values if they are unacceptable
      DO 1000 I=1,N
        IF(M(I).NE.2.AND.M(I).NE.4.AND.M(I).NE.8.AND.M(I).NE.16.
     1     AND.M(I).NE.32) M(I)=2
        IF(IND(I).LT.1) IND(I)=1
1000  CONTINUE
      IF(MAXPT.LE.0) MAXPT=MAXPTS
      NUM=0

!	Evaluate the integral
      CALL NGAUSS(A,B,N,M,IND,F,RINT,IER,NO,MAXPT)
      NUM=NUM+NO
      IF(IER.GT.100) RETURN
      IER=39

!	Iteration to check and improve the accuracy of integral
      DO 3000 I=1,10
        QC=.TRUE.
        DIF=0.0

!	Check for convergence along each dimension by doubling the
!	number of points at which function is evaluated
        DO 2000 J=1,N
          IF(NUM.GT.MAXPT) RETURN
          M1=M(J)
          I1=IND(J)
          IF(M(J).LT.32) THEN
!	If M(J)<32 then double the order of formula
            M(J)=2*M(J)
          ELSE
!	otherwise double the number of subintervals
            IND(J)=2*IND(J)
          ENDIF

          CALL NGAUSS(A,B,N,M,IND,F,RINT1,IER1,NO,MAXPT)
          NUM=NUM+NO
          DIF=DIF+ABS(RINT-RINT1)
          IF(IER1.GT.100) RETURN

          IF(ABS(RINT1-RINT).LT.MAX(AEPS,ABS(RINT1)*REPS)) THEN
!	If satisfactory accuracy is achieved then revert back to old
!	values of M(J) and IND(J)
            M(J)=M1
            IND(J)=I1
          ELSE
!	otherwise use new values
            RINT=RINT1
            QC=.FALSE.
          ENDIF
2000    CONTINUE

        IF(QC) THEN
!	If satisfactory accuracy is achieved for all dimensions then return
          IER=0
          RETURN
        ENDIF
!	If the number of function evaluations exceeds MAXPT then quit
        IF(NUM.GT.MAXPT) RETURN
3000  CONTINUE
      END
