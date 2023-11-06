!	To find all eigenvalues of an upper-Hessenberg matrix using QR algorithm
!
!	H : (input/output) Real array of length IH*NN containing the matrix
!		elements. The contents of array are destroyed during execution
!	NN : (input) Order of matrix
!	IH : (input) First dimension of array H as declared in the calling
!		program
!	ER : (output) Real array of length NN containing the real part of
!		the eigenvalues.
!	EI : (output) Real array of length NN containing the imaginary part of
!		the eigenvalues.
!	REPS : (input) Required relative accuracy in eigenvalues. It should
!		be of the order of machine accuracy
!	IER : (output) error parameter; IER=0 implies successful execution
!		IER=114 implies that N.LE.1 or N>IH, in which case no
!			calculations are done
!		IER=145 implies that QR iteration failed to converge at
!			some stage and calculations are abandoned.
!
!	Required routines : None

      SUBROUTINE HQR(H,NN,IH,ER,EI,REPS,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=30)
      DIMENSION H(IH,NN),ER(NN),EI(NN)

      T=0
      N=NN
      IF(N.LE.1.OR.N.GT.IH) THEN
        IER=114
        RETURN
      ENDIF
      IER=0
!	If all eigenvalues are found, then quit
2000  IF(N.EQ.0) RETURN

!	Loop for double QR iteration
      DO 5000 IT=1,NIT
!	Search for a small sub-diagonal element
        DO 2200 L=N,2,-1
          IF(ABS(H(L,L-1)).LE.REPS*(ABS(H(L-1,L-1))+ABS(H(L,L))))
     1       GO TO 2300
2200    CONTINUE
        L=1
2300    X=H(N,N)

        IF(L.EQ.N) THEN
!	One eigenvalue is isolated
          ER(N)=X+T
          EI(N)=0.0
          N=N-1
          GO TO 2000
        ENDIF

        Y=H(N-1,N-1)
        W=H(N,N-1)*H(N-1,N)
        IF(L.EQ.N-1) THEN
!	A pair of eigenvalues are isolated
          P=(Y-X)/2.
          D=P*P+W
          Y=SQRT(ABS(D))
          X=X+T
          IF(D.GT.0) THEN
!	Pair of real eigenvalues
            IF(P.LT.0) Y=-Y
            Y=P+Y
            ER(N-1)=X+Y
            ER(N)=X-W/Y
            EI(N-1)=0
            EI(N)=0
          ELSE
!	Pair of complex conjugate eigenvalues
            ER(N)=X+P
            ER(N-1)=ER(N)
            EI(N-1)=Y
            EI(N)=-Y
          ENDIF
          N=N-2
!	Find next eigenvalue
          GO TO 2000
        ENDIF

!	Apply special shifts
        IF(IT.EQ.10.OR.IT.EQ.20) THEN
          T=T+X
          DO 2500 I=1,N
2500      H(I,I)=H(I,I)-X
          S=ABS(H(N,N-1))+ABS(H(N-1,N-2))
          X=0.75D0*S
          Y=X
          W=-0.4375D0*S*S
        ENDIF

!	Search for two consecutive small sub-diagonal elements
        DO 2800 M=N-2,L,-1
          Z=H(M,M)
          R=X-Z
          S=Y-Z
          P=(R*S-W)/H(M+1,M)+H(M,M+1)
          U=H(M+1,M+1)-Z-R-S
          R=H(M+2,M+1)
          S=ABS(P)+ABS(U)+ABS(R)
          P=P/S
          U=U/S
          R=R/S
          IF(M.EQ.L) GO TO 3000
          IF(ABS(H(M,M-1))*(ABS(U)+ABS(R)).LE.REPS*ABS(P)*
     1       (ABS(H(M-1,M-1))+ABS(Z)+ABS(H(M+1,M+1)))) GO TO 3000
2800    CONTINUE
3000    DO 3200 I=M+3,N
          H(I,I-3)=0.0
3200    H(I,I-2)=0.0
        H(M+2,M)=0.0

!	Double QR transform for rows L to N and columns M to N
        DO 4000 K=M,N-1
          IF(K.NE.M) THEN
            P=H(K,K-1)
            U=H(K+1,K-1)
            R=0
            IF(K.NE.N-1) R=H(K+2,K-1)
            X=ABS(P)+ABS(U)+ABS(R)
            IF(X.EQ.0.0) GO TO 4000
            P=P/X
            U=U/X
            R=R/X
          ENDIF

          S=SQRT(P*P+U*U+R*R)
          IF(P.LT.0) S=-S
          IF(K.NE.M) THEN
            H(K,K-1)=-S*X
          ELSE IF(L.NE.M) THEN
            H(K,K-1)=-H(K,K-1)
          ENDIF
          P=P+S
          X=P/S
          Y=U/S
          Z=R/S
          U=U/P
          R=R/P

!	Row modification
          DO 3400 J=K,N
            P=H(K,J)+U*H(K+1,J)
            IF(K.NE.N-1) THEN
              P=P+R*H(K+2,J)
              H(K+2,J)=H(K+2,J)-P*Z
            ENDIF
            H(K+1,J)=H(K+1,J)-P*Y
            H(K,J)=H(K,J)-P*X
3400      CONTINUE

!	Column modification
          DO 3600 I=L,MIN(N,K+3)
            P=X*H(I,K)+Y*H(I,K+1)
            IF(K.NE.N-1) THEN
              P=P+Z*H(I,K+2)
              H(I,K+2)=H(I,K+2)-P*R
            ENDIF
            H(I,K+1)=H(I,K+1)-P*U
            H(I,K)=H(I,K)-P
3600      CONTINUE
4000    CONTINUE
5000  CONTINUE

!	QR iteration fails to converge
      IER=145
      END
