!	To calculate inverse of a square matrix
!
!	N : (input) Order of matrix
!	IA : (input) The first dimension of arrays A and AI as specified
!		in the calling program
!	A : (input) Real array of length IA*N containing the matrix
!	AI : (output) Real array of length IA*N which will contain the
!		calculated inverse of A
!	IWK : (output) Integer array of length N used as scratch space
!	IER : (output) The error parameter, IER=0 implies successful execution
!		Nonzero values may be set by subroutine GAUELM
!	It is possible to use CROUT instead of GAUELM for calculating
!	the triangular decomposition, but in that case an extra real scratch
!	array of size N will be required.
!
!	Required routines : GAUELM
 
      SUBROUTINE MATINV(N,IA,A,AI,IWK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(IA,N),AI(IA,N),IWK(N)
 
      DO 1000 I=1,N
        DO 800 J=1,N
          AI(J,I)=0.0
800     CONTINUE
        AI(I,I)=1.0
1000  CONTINUE
 
      NUM=N
      IFLG=0
      CALL GAUELM(N,NUM,A,AI,DET,IWK,IA,IER,IFLG)
 
      END
