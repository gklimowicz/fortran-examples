!	Reduction of real symmetric matrix to tridiagonal form using
!	Householder's method
!
!	A : (input/output) Real array of length IA*N containing the matrix
!		elements. After execution it will be overwritten by the
!		transformation matrix
!	N : (input) Order of matrix
!	IA : (input) The first dimension of A as specified in the calling program
!	D : (output) Diagonal elements of the transformed tridiagonal matrix
!		D(I) would contain A(I,I)
!	E : (output) Off-diagonal elements of the transformed tridiagonal matrix
!		E(I+1) would contain A(I,I+1)=A(I+1,I)
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=107 implies that N.LE.1 or N>IA, in which case no
!			calculations are done
!
!	Required routines : None

      SUBROUTINE TRED2(A,N,IA,D,E,IER)
      IMPLICIT REAL*8(A-H,O-Z)
!	For REAL*4 use REPS=1.E-30
!      PARAMETER(REPS=1.E-30)
      PARAMETER(REPS=1.D-300)
      DIMENSION A(IA,N),D(N),E(N)

      IF(N.LE.1.OR.N.GT.IA) THEN
        IER=107
        RETURN
      ENDIF
      IER=0

      DO 4000 I=N,2,-1
        F=A(I,I-1)
        G=0
        DO 2000 K=1,I-2
2000    G=G+A(I,K)*A(I,K)
        H=G+F*F
        IF(G.LE.REPS) THEN
!	Skip the transformation
          E(I)=F
          H=0
        ELSE

          G=SQRT(H)
          IF(F.GT.0.0) G=-G
          E(I)=G
          H=H-F*G
          A(I,I-1)=F-G
          F=0
          DO 2600 J=1,I-1
!	Elements of u_i/H_i
            A(J,I)=A(I,J)/H
            G=0
!	Form elements of A_iu_i
            DO 2200 K=1,J
2200        G=G+A(J,K)*A(I,K)
            DO 2400 K=J+1,I-1
2400        G=G+A(K,J)*A(I,K)
!	Components of p_i
            E(J)=G/H
            F=F+G*A(J,I)
2600      CONTINUE

!	calculate u_i^Tp_i/2H_i
          HH=0.5*F/H
          DO 3400 J=1,I-1
            F=A(I,J)
            G=E(J)-HH*F
!	Elements of q_i
            E(J)=G
            DO 3000 K=1,J
3000        A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
3400      CONTINUE
        ENDIF
        D(I)=H
4000  CONTINUE

      D(1)=0
      E(1)=0
!	accumulation of transformation matrix Q
      DO 5000 I=1,N
        IF(D(I).NE.0.0) THEN
          DO 4600 J=1,I-1
            G=0
            DO 4200 K=1,I-1
4200        G=G+A(I,K)*A(K,J)
            DO 4400 K=1,I-1
4400        A(K,J)=A(K,J)-G*A(K,I)
4600      CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1
        DO 4800 J=1,I-1
          A(I,J)=0.0
4800    A(J,I)=0.0
5000  CONTINUE
      END
