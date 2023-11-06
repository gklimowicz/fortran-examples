!	To balance a general real matrix by exact diagonal similarity transformation
!
!	A : (input/output) Real array of length IA*N containing the matrix
!		elements. After execution it will contain the balanced matrix
!	N : (input) Order of matrix
!	IA : (input) First dimension of array A as declared in the calling
!		program
!	B : (input) Base of the arithmetic to be used for balancing.
!		Should generally be 2.
!	LOW : (output) Index of column such that in balanced matrix A(i,j)=0
!		if i>j and j<LOW, the first LOW-1 eigenvalues are isolated
!	IGH : (output) Index of row such that in balanced matrix A(i,j)=0
!		if i>j and i>IGH, the last N-IGH eigenvalues are isolated
!		After balancing only the sub-matrix between rows and columns
!		LOW to IGH need to be considered
!	D : (output) Real array of length N containing information about
!		transformations used for balancing. D(LOW) to D(IGH)
!		contain elements used for balancing, while other elements
!		contain information about permutations used
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=112 implies that N.LE.1 or N>IA
!
!	Required routines : None

      SUBROUTINE BALANC(A,N,IA,B,LOW,IGH,D,IER)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(GAM=0.95D0)
      DIMENSION A(IA,N),D(N)

      IF(N.LE.1.OR.N.GT.IA) THEN
        IER=112
        RETURN
      ENDIF
      IER=0
      B2=B*B
      LOW=1
      IGH=N

!	Search for rows isolating an eigenvalue
2000  DO 3000 J=IGH,1,-1
        R=0
        DO 2100 I=1,J-1
2100    R=R+ABS(A(J,I))
        DO 2200 I=J+1,IGH
2200    R=R+ABS(A(J,I))

        IF(R.EQ.0.0) THEN
!	Push the row to bottom
          D(IGH)=J
          IF(J.NE.IGH) THEN
            DO 2400 I=1,IGH
              T=A(I,J)
              A(I,J)=A(I,IGH)
              A(I,IGH)=T
2400        CONTINUE
            DO 2600 I=LOW,N
              T=A(J,I)
              A(J,I)=A(IGH,I)
              A(IGH,I)=T
2600        CONTINUE
          ENDIF
          IGH=IGH-1
          GO TO 2000
        ENDIF
3000  CONTINUE

!	Search for columns isolating an eigenvalue
3200  DO 4000 J=LOW,IGH
        C=0
        DO 3300 I=LOW,J-1
3300    C=C+ABS(A(I,J))
        DO 3400 I=J+1,IGH
3400    C=C+ABS(A(I,J))

!	Move the column to the left end
        IF(C.EQ.0) THEN
          D(LOW)=J
          IF(J.NE.LOW) THEN
            DO 3600 I=1,IGH
              T=A(I,J)
              A(I,J)=A(I,LOW)
              A(I,LOW)=T
3600        CONTINUE
            DO 3800 I=LOW,N
              T=A(J,I)
              A(J,I)=A(LOW,I)
              A(LOW,I)=T
3800        CONTINUE
          ENDIF
          LOW=LOW+1
          GO TO 3200
        ENDIF
4000  CONTINUE

!	Balance the sub-matrix in rows LOW to IGH
      DO 4200 I=LOW,IGH
4200  D(I)=1
4300  QC=.FALSE.
      DO 5000 I=LOW,IGH
        C=0
        R=0
        DO 4400 J=LOW,IGH
          IF(J.NE.I) THEN
            C=C+ABS(A(J,I))
            R=R+ABS(A(I,J))
          ENDIF
4400    CONTINUE
        G=R/B
        F=1
        S=C+R
4500    IF(C.LT.G) THEN
          F=F*B
          C=C*B2
          GO TO 4500
        ENDIF
        G=R*B
4600    IF(C.GE.G) THEN
          F=F/B
          C=C/B2
          GO TO 4600
        ENDIF

!	Apply the transformation
        IF((C+R)/F.LT.GAM*S) THEN
          G=1./F
          D(I)=D(I)*F
          QC=.TRUE.
          DO 4800 J=LOW,N
4800      A(I,J)=A(I,J)*G
          DO 4900 J=1,IGH
4900      A(J,I)=A(J,I)*F
        ENDIF
5000  CONTINUE

!	perform one more cycle of balancing
      IF(QC) GO TO 4300
      END
