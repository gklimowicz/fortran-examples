!	To find eigenvalues and eigenvectors of Z T Z^T using QL algorithm
!	where T is a symmetric tridiagonal matrix and Z is an orthogonal matrix.
!	If Z is the transformation matrix to reduce original matrix to
!	tridiagonal matrix, it will calculate the eigenvectors of original matrix
!
!	Z : (input/output) Real array of length IZ*N which should contain
!		the transformation matrix required to reduce original real
!		symmetric matrix to tridiagonal form. To find eigenvectors
!		of a symmetric tridiagonal matrix, set Z to unit matrix
!		After execution Z will contain the eigenvector of the original
!		matrix Z T Z^T. Z(i,j) should contain the ith component of
!		jth eigenvector
!	N : (input) Order of matrix
!	IZ : (input) The first dimension of array Z as declared in the
!		calling program. (IZ.GE.N)
!	D : (input/output) Real array of length N, containing the diagonal
!		elements of the tridiagonal matrix, D(i)=T(i,i).
!		After execution it will contain the eigenvalues of the matrix
!	E : (input/output) Real array of length N containing the off-diagonal
!		elements of the tridiagonal matrix, E(i)=T(i,i+1)=T(i+1,i)
!		It is used as scratch space and its contents will be destroyed
!		during execution.
!	REPS : (input) Required tolerance, it should be of order of machine
!		accuracy
!	IER : (output) Error parameter; IER=0 implies successful execution
!		IER=108 implies that N.LE.1 or N>IZ, in which case no
!			calculations are performed
!		IER=143 implies that the QL algorithm failed to converge
!			for some eigenvalue, the calculations are abandoned
!
!	Required routines : None
!	
      SUBROUTINE TQL2(Z,N,IZ,D,E,REPS,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=30)
      DIMENSION Z(IZ,N),D(N),E(N)

      IF(N.LE.0.OR.N.GT.IZ) THEN
        IER=108
        RETURN
      ENDIF
      IER=0
      DO 2000 I=2,N
2000  E(I-1)=E(I)
      E(N)=0
      B=0
      F=0

      DO 4000 L=1,N
        H=REPS*(ABS(D(L))+ABS(E(L)))
        IF(B.LT.H) B=H
!	Look for small off-diagonal elements
        DO 2200 M=L,N
          IF(ABS(E(M)).LE.B) GO TO 2300
2200    CONTINUE
        M=N
2300    IF(M.EQ.L) THEN
!	one eigenvalue is isolated
          D(L)=D(L)+F
          GO TO 4000
        ENDIF

!	Loop for QL transformation 
        DO 3600 IT=1,NIT
!	Find shift
          G=D(L)
          P=(D(L+1)-G)/(2*E(L))
          R=SQRT(P*P+1)
          IF(P.LT.0.0) R=-R
          D(L)=E(L)/(P+R)
          H=G-D(L)
          DO 2600 I=L+1,N
2600      D(I)=D(I)-H
          F=F+H

!	The QL transformation
          P=D(M)
          C=1
          S=0
!	Given's rotations
          DO 3400 I=M-1,L,-1
            G=C*E(I)
            H=C*P
            IF(ABS(P).GE.ABS(E(I))) THEN
              C=E(I)/P
              R=SQRT(C*C+1)
              E(I+1)=S*P*R
              S=C/R
              C=1/R
            ELSE
              C=P/E(I)
              R=SQRT(C*C+1)
              E(I+1)=S*E(I)*R
              S=1./R
              C=C/R
            ENDIF
            P=C*D(I)-S*G
            D(I+1)=H+S*(C*G+S*D(I))

!	Transforming the eigenvectors
            DO 3000 K=1,N
              H=Z(K,I+1)
              Z(K,I+1)=S*Z(K,I)+C*H
              Z(K,I)=C*Z(K,I)-S*H
3000        CONTINUE
3400      CONTINUE

          E(L)=S*P
          D(L)=C*P
          IF(ABS(E(L)).LE.B) THEN
!	One eigenvalue is isolated
            D(L)=D(L)+F
            GO TO 4000
          ENDIF
3600    CONTINUE
!	QL iteration fails to converge
        IER=143
        RETURN
4000  CONTINUE

!	Sort eigenvalues in ascending order by straight selection
      DO 5000 I=1,N-1
        K=I
        P=D(I)
        DO 4200 J=I+1,N
          IF(D(J).LT.P) THEN
            K=J
            P=D(J)
          ENDIF
4200    CONTINUE
        IF(K.NE.I) THEN
!	exchange the eigenvalues and eigenvectors
          D(K)=D(I)
          D(I)=P
          DO 4400 J=1,N
            P=Z(J,I)
            Z(J,I)=Z(J,K)
            Z(J,K)=P
4400      CONTINUE
        ENDIF
5000  CONTINUE
      END
