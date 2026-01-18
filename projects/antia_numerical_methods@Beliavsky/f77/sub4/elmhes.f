!	To reduce a general real matrix to Hessenberg form using stabilised
!	elementary transformations
!	It is advisable to balance the matrix before applying this transformations
!
!	A : (input/output) Real array of length IA*N containing the matrix
!		After execution the reduced matrix will be overwritten on
!		the same array
!	N : (input) Order of matrix
!	IA : (input) First dimension of array A as declared in the calling
!		program
!	LOW, IGH : (input) After balancing only rows and columns from
!		LOW to IGH need to be considered, since other rows or
!		columns contain isolated eigenvalues. If the matrix is
!		not balanced use LOW=1, IGH=N
!	INC : (output) Integer array of length N, which will contain information
!		about row and column interchanges during reduction
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=113 implies that N.LE.1 or N>IA
!
!	Required routines : None

      SUBROUTINE ELMHES(A,N,IA,LOW,IGH,INC,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(IA,N),INC(N)

      IF(N.LE.1.OR.N.GT.IA) THEN
        IER=113
        RETURN
      ENDIF
      IER=0
      IF(LOW.GT.IGH-2) RETURN
      DO I=1,N
        INC(I)=I
      ENDDO

      DO 5000 M=LOW+1,IGH-1
        I=M
!	Find the pivot
        AMAX=0
        DO 2000 J=M,IGH
          IF(ABS(A(J,M-1)).GT.ABS(AMAX)) THEN
            AMAX=A(J,M-1)
            I=J
          ENDIF
2000    CONTINUE
        INC(M)=I

        IF(I.NE.M) THEN
!	Interchange the corresponding rows and columns
          DO 2200 J=M-1,N
            T=A(I,J)
            A(I,J)=A(M,J)
            A(M,J)=T
2200      CONTINUE
          DO 2400 J=1,IGH
            T=A(J,I)
            A(J,I)=A(J,M)
            A(J,M)=T
2400      CONTINUE
        ENDIF

        IF(AMAX.NE.0.0) THEN
!	Perform Gaussian elimination
          DO 4000 I=M+1,IGH
            T=A(I,M-1)
            IF(T.NE.0.0) THEN
              T=T/AMAX
              A(I,M-1)=T
              DO 3000 J=M,N
3000          A(I,J)=A(I,J)-T*A(M,J)
              DO 3200 J=1,IGH
3200          A(J,M)=A(J,M)+T*A(J,I)
            ENDIF
4000      CONTINUE
        ENDIF
5000  CONTINUE
      END
