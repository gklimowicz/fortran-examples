!	To find specified eigenvalues and eigenvectors of a real symmetric
!	tridiagonal matrix using Sturm sequence coupled with inverse iteration
!
!	E : (input) Real array of length N containing the off-diagonal elements
!		of the tridiagonal matrix, E(i+1)=A(i,i+1)=A(i+1,i)
!	D : (input) Real array of length N containing the diagonal elements
!		of the tridiagonal matrix, D(i)=A(i,i)
!	N : (input) Order of the matrix
!	M1 : (input) Serial number of lowest eigenvalue to be determined.
!		The eigenvalues are sorted in increasing order
!	M2 : (input) Serial number of highest eigenvalue to be determined.
!		All eigenvalues from M1 to M2 are determined
!	EI : (output) Real array of length M2-M1+1 containing the calculated
!		eigenvalues
!	EPS1 : (input) Relative accuracy to which eigenvalues are located by bisection
!		before using inverse iteration. If the inverse iteration
!		does not converge to nearest eigenvalues, EPS1 can be reduced.
!		A value of 0.1-0.01 times typical eigenvalue is generally sufficient.
!	REPS : (input) Desired relative accuracy in eigenvalues and eigenvectors
!	EV : (output) Real array of length IV*(M2-M1+1) containing the
!		eigenvectors. EV(i,j) should contain the ith component of
!		the jth eigenvector
!	IV : (input) The first dimension of array EV as declared in the
!		calling program
!	WK : Real array of length 7N used as scratch space
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=109 implies that N.LE.1 or N>IV or M1<1 or M2>N
!		Other values may be set by TINVIT, only the last value is
!		returned.
!
!	Required routines : STURM, TINVIT, RAN1
!
!
      SUBROUTINE TRIDIA(E,D,N,M1,M2,EI,EPS1,REPS,EV,IV,WK,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION E(N),D(N),WK(N,7),EI(M2-M1+1),EV(IV,M2-M1+1)

      IF(N.LE.1.OR.N.GT.IV.OR.M1.LT.1.OR.M2.GT.N) THEN
        IER=109
        RETURN
      ENDIF
      IER=0
      IF(M1.GT.M2) RETURN
!	Locate the eigenvalues
      CALL STURM(E,D,N,M1,M2,WK(1,2),WK(1,3),NUM,EPS1,WK,IER1)

!	Loop for finding individual eigenvalues and eigenvectors
      DO 2000 I=M1,M2
        IFLG=0
        I1=I-M1+1

        IF(I.GT.M1) THEN
          IF(ABS(WK(I,2)-WK(I-1,2)).LT.3.*REPS*ABS(WK(I,2))) THEN
!	Set the flag for close eigenvalues
            IFLG=1
            DO 1500 J=1,N
1500        EV(J,I1)=EV(J,I1-1)
          ENDIF
        ENDIF
        CALL TINVIT(E,D,N,WK(I,2),WK(I,3),EI(I1),EV(1,I1),REPS,
     1              IFLG,IER1,WK(1,1),WK(1,4),NUM)
        IF(IER1.GT.0) IER=IER1
2000  CONTINUE
      END
