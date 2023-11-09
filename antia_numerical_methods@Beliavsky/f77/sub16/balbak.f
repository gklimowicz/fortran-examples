!	Perform back-transformation of a set of right eigenvectors from
!	those of balanced matrix to that for original matrix
!
!	N : (input) Order of matrix
!	LOW, IGH : (input) After balancing only rows and columns from
!		LOW to IGH need to be considered, since other rows or
!		columns contain isolated eigenvalues
!	CZ : (input/output) Complex array of length IZ*M containing the
!		eigenvectors of balanced matrix. After execution it will
!		be overwritten by eigenvectors of the original matrix.
!	M : (input) Number of eigenvectors to be balanced
!	IZ : (input) First dimension of array CZ as declared in the calling
!		program
!	D : (input) Real array of length N, containing information about
!		transformation used for balancing.
!
!	Required routines : None

      SUBROUTINE BALBAK(N,LOW,IGH,CZ,M,IZ,D)
      IMPLICIT REAL*16(A,B,D-H,O-Z)
      IMPLICIT COMPLEX*32(C)
      DIMENSION CZ(IZ,M),D(N)

      DO 3000 I=LOW,IGH
        S=D(I)
!	For back-transforming left eigenvectors use the following
!	statement instead of the preceding one
!       S=1/D(I)

        DO 2200 J=1,M
2200    CZ(I,J)=CZ(I,J)*S
3000  CONTINUE

      DO 3400 I=LOW-1,1,-1
        K=D(I)
        IF(K.NE.I) THEN
!	Exchange the corresponding rows
          DO 3200 J=1,M
            CS=CZ(I,J)
            CZ(I,J)=CZ(K,J)
3200      CZ(K,J)=CS
        ENDIF
3400  CONTINUE

      DO 4000 I=IGH+1,N
        K=D(I)
        IF(K.NE.I) THEN
!	Exchange the corresponding rows
          DO 3800 J=1,M
            CS=CZ(I,J)
            CZ(I,J)=CZ(K,J)
3800      CZ(K,J)=CS
        ENDIF
4000  CONTINUE
      END
