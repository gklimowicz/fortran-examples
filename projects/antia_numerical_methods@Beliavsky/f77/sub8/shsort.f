!	To sort an array of length N in ascending order, using shell sort algorithm
!
!	A : (input/output) Real array of length N, which is to be sorted
!               The sorted array will be overwritten on A
!	N : (input) No. of elements of array A to be sorted
!
!	Required routines : None

      SUBROUTINE SHSORT(A,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N)

      L=LOG(FLOAT(N))/LOG(2.D0)+1
      M=2**L
      IF(M.GE.2*N) M=M/2
      DO J=1,L
        M=M/2
        IF(M.EQ.0) RETURN
        DO I=M+1,N
          T=A(I)
          DO I1=I,M+1,-M
            IF(A(I1-M).GT.T) THEN
              A(I1)=A(I1-M)
            ELSE
              EXIT
            ENDIF
          ENDDO
          IF(I1.LE.0) I1=I1+M
          A(I1)=T
        ENDDO
      ENDDO
      END

