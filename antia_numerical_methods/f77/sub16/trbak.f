!	To perform back-transformation on eigenvectors of reduced tridiagonal matrix
!	to obtain eigenvectors of original real symmetric matrix reduced by TRED2
!	This routine is not required if eigenvectors are computed using TQL2
!
!	A : (input/output) Real array of length IA*N containing the transformation
!		matrix. The last column of this array is used as scratch
!		space and will be destroyed during execution
!	IA : (input) The first dimension of array A as declared in the
!		calling program (IA.GE.N)
!	N : (input) Order of the matrix
!	Z : (input/output) Real array of length IZ*NZ containing the
!		eigenvectors of tridiagonal matrix. After execution it will
!		be overwritten by eigenvectors of original matrix
!	IZ : (input) First dimension of array Z as declared in the calling
!		program (IZ.GE.N)
!	NZ : (input) Number of eigenvectors to be transformed
!
!	Required routines : None
!
      SUBROUTINE TRBAK(A,IA,N,Z,IZ,NZ)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(IA,N),Z(IZ,NZ)

!	Loop on eigenvectors
      DO 4000 I=1,NZ
!	The matrix multiplication
        DO 3000 J=1,N
          S=0.0
          DO 2200 K=1,N-1
2200      S=S+A(J,K)*Z(K,I)
!	To take care of the last column of A, which is overwritten
          IF(J.EQ.N) S=S+Z(N,I)
          A(J,N)=S
3000    CONTINUE

        DO 3200 J=1,N
3200    Z(J,I)=A(J,N)
4000  CONTINUE
      END
