!	To find cascade sum of a series
!
!	TERM : (input) Name of function routine to calculate the Ith term
!	      It is possible to supply the terms in a real array TERM
!	      if the dimension statement is uncommented. (see CASSUM_A)
!	N : (input) Number of terms to be summed
!
!	FUNCTION TERM(I) to calculate the Ith term must be supplied by
!	the user
!
!	Required routines : TERM

      FUNCTION CASSUM(TERM,N)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(N2MAX=30)
      DIMENSION S(N2MAX),Q(N2MAX)
!	Uncomment this line if TERM is required to be an array 
!	DIMENSION TERM(N)

      CASSUM=0.0
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

          CASSUM=CASSUM+S(N2MAX)
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
      IF(.NOT.Q(N2MAX)) CASSUM=CASSUM+S(N2MAX)
      END
