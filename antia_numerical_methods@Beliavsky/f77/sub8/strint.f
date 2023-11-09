!	Multiple integration over a hyper-rectangle in n-dimensions
!	using compound monomial rules
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions
!	M : (input/output) Integer specifying the formula to be used
!		M can be 1, 3 or 5, otherwise it will be set to a default value of 3
!		M=1 selects 1-point formula of degree 1
!		M=3 selects 2N-point formula of degree 3 due to Stroud
!		M=5 selects (2N*N+1)-point formula of degree 5
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
!		IER=308 implies that number of points exceeded MAXPT in
!			first attempt and no approximation of RINT is calculated
!		IER=309 implies N<1 or N>50 and no calculations are done
!	NUM : (output) Number of function evaluations used by subroutine
!	MAXPT : (input/output) Maximum number of function evaluations permitted
!		If MAXPT <1 it is set to a default value of MAXPTS (=1000000)
!
!	FUNCTION F(N,X) must be supplied by the user
!	
!	Required routines : STROUD, F
!
      SUBROUTINE STRINT(A,B,N,M,IND,F,RINT,REPS,AEPS,DIF,IER,NUM,MAXPT)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL F
      PARAMETER(MAXPTS=1100000)
      DIMENSION A(N),B(N),IND(N)

!	set M and IND(I) to default values if they are unacceptable
      IF(M.NE.1.AND.M.NE.3.AND.M.NE.5) M=3
      DO 1000 I=1,N
        IF(IND(I).LT.1) IND(I)=1
1000  CONTINUE
      IF(MAXPT.LE.0) MAXPT=MAXPTS
      NUM=0

!	Evaluate the first approximation to the integral
      CALL STROUD(A,B,N,M,IND,F,RINT,IER,NO,MAXPT)
      NUM=NUM+NO
      IF(IER.GT.100) RETURN
      IER=39

!	Iteration to check and improve the accuracy of integral
      DO 3000 I=1,10
        QC=.TRUE.
        DIF=0.0

        DO 2000 J=1,N
          IF(NUM.GT.MAXPT) RETURN
          I1=IND(J)
!	Double the number of subintervals along the Jth axis
          IND(J)=2*IND(J)

          CALL STROUD(A,B,N,M,IND,F,RINT1,IER1,NO,MAXPT)
          NUM=NUM+NO
          DIF=DIF+ABS(RINT-RINT1)
          IF(IER1.GT.100) RETURN

          IF(ABS(RINT1-RINT).LT.MAX(AEPS,ABS(RINT1)*REPS)) THEN
!	If satisfactory accuracy is achieved then restore the old value of IND(J)
            IND(J)=I1
          ELSE
!	else retain the new value
            RINT=RINT1
            QC=.FALSE.
          ENDIF
2000    CONTINUE

        IF(QC) THEN
!	If satisfactory accuracy is achieved, then return
          IER=0
          RETURN
        ENDIF
        IF(NUM.GT.MAXPT) RETURN
3000  CONTINUE
      END
