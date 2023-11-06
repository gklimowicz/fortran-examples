!	To find eigenvalues and eigenvectors of a complex Hermitian matrix
!
!	ZA : (input) Complex array of length IA*N containing the matrix
!	N : (input) Order of the matrix
!	IA : (input) First dimension of array ZA as declared in the calling program
!	EI : (output) Real array of length N containing the eigenvalues
!	ZV : (output) Complex array of length IZ*N containing the eigenvectors
!		ZV(i,j) is the ith component of Jth eigenvector
!	IZ : (input) First dimension of array ZV as declared in the calling program
!	WK : Real array of length 2N*(2N+2) used as scratch space
!	REPS : (input) Required tolerance, should be of the order of machine
!		accuracy
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=111 implies that N.LE.1 or N>IA or N>IZ
!			in which case no calculations are done
!		Other values may be set by TRED2 or TQL2
!
!	Required routines : TRED2, TQL2
 
      SUBROUTINE HEREVP(ZA,N,IA,EI,ZV,IZ,WK,REPS,IER)
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      DIMENSION ZA(IA,N),WK(2*N,2*N+2),EI(N),ZV(IZ,N)
 
      IF(N.LE.1.OR.N.GT.IA.OR.N.GT.IZ) THEN
        IER=111
        RETURN
      ENDIF

!	Setup the 2N*2N real matrix
      DO 2000 I=1,N
        DO 2000 J=1,N
          WK(I,J)=ZA(I,J)
          WK(I+N,J+N)=WK(I,J)
          WK(I+N,J)=IMAG(ZA(I,J))
          WK(I,J+N)=-WK(I+N,J)
2000  CONTINUE
      N2=2*N
      IW=N2

!	To reduce the 2N*2N matrix to tridiagonal form
      CALL TRED2(WK,N2,IW,WK(1,N2+1),WK(1,N2+2),IER)
      IF(IER.GT.100) RETURN

!	Find eigenvalues and eigenvectors of tridiagonal matrix
      CALL TQL2(WK,N2,IW,WK(1,N2+1),WK(1,N2+2),REPS,IER)
      IF(IER.GT.100) RETURN
 
!	Since all eigenvalues are repeated and sorted in ascending order
!	pick alternate eigenvalues and eigenvectors
!	May create problem for multiple eigenvalues.
      ZI=(0.D0,1.D0)
      DO 3000 I=1,N
        I1=2*I-1
        EI(I)=WK(I1,N2+1)
        DO 2500 J=1,N
          ZV(J,I)=WK(J,I1)+ZI*WK(J+N,I1)
2500    CONTINUE
3000  CONTINUE
 
      END
