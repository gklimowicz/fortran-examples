!	To find cascade sum of a series provided as an array of terms
!
!	TERM : (input) Real array of length N containing the terms to be summed
!	N : (input) Number of terms to be summed
!
!	Required routines : None


      FUNCTION CASSUM_A(TERM,N)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(N2MAX=30)
      DIMENSION S(N2MAX),Q(N2MAX)
!	Comment this line if TERM is required to be a function
      DIMENSION TERM(N)

      CASSUM_A=0.0
      DO 1000 I=1,N2MAX
1000  Q(I)=.TRUE.

      DO 2000 I=1,N
        IF(Q(1)) THEN
          S(1)=TERM(I)
          Q(1)=.FALSE.

!	If a pair is formed add the sum to higher level in the binary tree
        ELSE
          S(1)=S(1)+TERM(I)
          Q(1)=.TRUE.
          DO 1500 J=2,N2MAX
            IF(Q(J)) THEN
              S(J)=S(J-1)
              Q(J)=.FALSE.
              GO TO 2000
            ELSE
              S(J)=S(J)+S(J-1)
              Q(J)=.TRUE.
            ENDIF
1500      CONTINUE

          CASSUM_A=CASSUM_A+S(N2MAX)
        ENDIF

2000  CONTINUE

!	Find the sum by adding up the incomplete pairs
      DO 3000 J=1,N2MAX-1
        IF(.NOT.Q(J)) THEN
          IF(Q(J+1)) THEN
            S(J+1)=S(J)
            Q(J+1)=.FALSE.
          ELSE
            S(J+1)=S(J+1)+S(J)
          ENDIF
        ENDIF
3000  CONTINUE
      IF(.NOT.Q(N2MAX)) CASSUM_A=CASSUM_A+S(N2MAX)
      END
