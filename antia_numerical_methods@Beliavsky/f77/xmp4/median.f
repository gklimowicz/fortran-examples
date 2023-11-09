!       To calculate median, mean and higher moments of a distribution
      PROGRAM MEDIAN
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(500)
      DATA (X(I),I=1,21)/49,73,46,45,67,47,28,49,61,71,49,57,52,79,
     1    40,41,55,69,41,71,40/

51    FORMAT(' Median =',F8.2,',  Mean =',F8.2,',  Standard deviation ='
     1       ,F8.2,/,'  Skewness =',F8.2,',  Kurtosis =',F8.2)

      N=21
      CALL SHSORT(X,N)
      AMED=X((N+1)/2)
      AMEAN=0.0
      DO I=1,N
        AMEAN=AMEAN+X(I)
      ENDDO
      AMEAN=AMEAN/N
      SIG=0.0
      DO I=1,N
        SIG=SIG+(X(I)-AMEAN)**2
      ENDDO
      SIG=SQRT(SIG/(N-1))
      SKEW=0.0
      AKER=0.0
      DO I=1,N
        SKEW=SKEW+((X(I)-AMEAN)/SIG)**3
        AKER=AKER+((X(I)-AMEAN)/SIG)**4
      ENDDO
      SKEW=SKEW/(N-1)
      AKER=AKER/(N-1)-3
      WRITE(6,51) AMED,AMEAN,SIG,SKEW,AKER
      END

!       -------------------------------------------------------

!	To sort an array of length N in ascending order, using shell sort algorithm
!
!	A : (input/output) Real array of length N, which is to be sorted
!               The sorted array will be overwritten on A
!	N : (input) No. of elements of array A to be sorted
!
!	Required routines : None

      SUBROUTINE SHSORT(A,N)
!      IMPLICIT REAL*8(A-H,O-Z)
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

