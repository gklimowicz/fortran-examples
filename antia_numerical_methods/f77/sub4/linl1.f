!	To calculate coefficients of linear L_1 approximation in terms of
!	specified basis functions for a tabulated function
!
!	M : (input) Number of basis functions in the required approximation
!	N : (input) Number of points in the table of values. N> M
!	A : (output) Real array of length M+1 containing the coefficients
!		of approximation. A(I) is the coefficient of Phi_I(x) in the
!		approximation
!	F : (input) Real array of length N containing the function values
!	G : (input) Real array of length IG*N, containing the values of
!		basis functions at each point in table. G(I,J) is the
!		value of Phi_I(x) at Jth tabular point.
!	IG : (input) First dimension of array G as declared in the calling
!		program, IG .GE. M
!	EPS : (input) Estimate of roundoff error, passed on to SIMPL1
!	ESUM : (output) L_1 norm of the residual at the calculated approximation
!	IER : (output) error parameter, IER=0 implies successful execution
!		IER=616 implies that M<1  or N<M+1
!			in this case no calculations are done
!		Other values of IER may be set by SIMPL1
!	WK : Real array of length (N+2)*(M+3) used as scratch space
!	IWK : Integer array of length N+M+3 used as scratch space
!
!	Required routines : SIMPL1
!
 
      SUBROUTINE LINL1(M,N,A,F,G,IG,EPS,ESUM,IER,WK,IWK)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION A(M+2),F(N),IWK(N+M+3),WK(N+2,M+3),X(N),G(IG,N)
 
      IF(M.LE.0.OR.N.LT.M+1) THEN
        IER=616
        RETURN
      ENDIF
 
      IER=0
      LJ=N+2
      NV=M+1
      M3=N
 
!	Setting up the tableau for simplex algorithm
      DO 2100 J=1,M+2
2100  WK(1,J)=0.0
      DO 2300 I=1,N
        FI=F(I)
        SI=SIGN(1.0,FI)
        IWK(I+1)=I*SI
        WK(I+1,1)=FI*SI
        WK(1,1)=WK(1,1)-FI*SI
        S=0.0
        DO 2200 J=1,M
          TI=G(J,I)*SI
          S=S+TI
          WK(I+1,J+1)=TI
          WK(1,J+1)=WK(1,J+1)-TI
2200    CONTINUE
        WK(I+1,M+2)=-S
        WK(1,M+2)=WK(1,M+2)+S
2300  CONTINUE
 
      DO 2400 J=1,M+1
        IWK(N+1+J)=N+J
2400  CONTINUE
 
      NC=NV+M3
      CALL SIMPL1(WK,LJ,NC,M3,IWK,IWK(M3+1),IER,EPS)

!	L_1 norm of the residual
      ESUM=-WK(1,1)
!	Finding the coefficients from the tableau
      DO 3200 I=M3+2,NC+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=0.0
3200  CONTINUE
      DO 3400 I=2,M3+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=WK(I,1)
3400  CONTINUE
      DO 3600 I=1,M
3600  A(I)=A(I)-A(M+1)
      END
