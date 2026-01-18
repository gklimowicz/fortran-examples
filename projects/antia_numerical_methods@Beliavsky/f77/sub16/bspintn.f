!	To calculate coefficients for B-spline interpolation in n-dimensions
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N giving the number of
!		tabular points in each dimension
!	X : (input) Array of length NXD*N containing the abscissas
!		X(I,J) is the Ith abscissa in Jth dimension
!	NXD : (input) First dimension of arrays X, XF, INTX as specified
!		in the calling program. NXD must be greater than or equal
!		to the maximum of NK(I)
!	F : (input) Array of length NK(1)*NK(2)* ... *NK(N) containing
!		the table of values. The dimensions of F in the calling
!		program must exactly match the size of table so that
!		there are no gaps in memory allocation.
!	K : (input) Order of B-spline required. K=4 gives cubic B-splines
!	AX : (output) Real array of length NXD*3*K*N containing the
!		triangular decomposition of equation matrix in band form
!		for each dimension.
!	C : (output) Real array of length NK(1)*NK(2)*... *NK(N)
!		containing the coefficients of expansion.
!	XF : (output) Real array of size NXD*N, containing
!		the knots used for B-spline calculations.
!		XF(I,J) is the Ith knot along Jth dimension.
!	MK : (output) Integer array of length N containing number of
!		knots for B-splines in each dimension.
!	INTX : (output) Integer array of length NXD*N containing information
!		about pivoting during solution of system of linear equations
!		for each dimension.
!	WK : Scratch array of length NK(1)*NK(2)*...*NK(N)+3*K
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values may be set by BSPINT, BSPLIN or GAUBND
!
!	Required routines : BSPINT, BSPLIN, GAUBND

 
      SUBROUTINE BSPINTN(N,NK,X,NXD,F,K,AX,C,XF,MK,INTX,WK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(NXD,N),NK(N),F(*),AX(NXD*3*K,N),XF(NXD,N),C(*),
     1        MK(N),WK(*),INTX(NXD,N)
 
!	Calculate the triangular decomposition of matrices for interpolation
!	along each dimension
      DO 2000 I=1,N
        IFLG1=1
        LJ=NK(I)
        CALL BSPINT(NK(I),X(1,I),F,K,AX(1,I),LJ,C,XF(1,I),MK(I),
     1         IFLG1,INTX(1,I),WK,IER)
        IF(IER.GT.100) RETURN
2000  CONTINUE
 
      KB=K-1
      NX=1
      DO 2200 I=1,N-1
        NX=NX*NK(I)
2200  CONTINUE
      N0=NX*NK(N)+1
      NY=1
      JU=N
 
!	If N is odd interpolate along last dimension outside the loop
      IF(MOD(N,2).EQ.1) THEN
        LJ=NK(N)
!	Set up the RHS
        DO 2500 I=1,NX
          DO 2500 J=1,NK(N)
            WK(J+(I-1)*LJ)=F(I+(J-1)*NX)
2500    CONTINUE
        NUM=NX
        IFLG1=2
        CALL GAUBND(NK(N),KB,NUM,AX(1,N),WK,DET,IDET,INTX(1,N),LJ,
     1         IER,IFLG1,WK(N0))
        IF(IER.GT.100) RETURN
 
!	Set up the RHS for next interpolation along N-1 th dimension
        NX1=NX/NK(N-1)
        NY=NK(N)
        LJ1=NK(N-1)
        DO 3000 I1=1,NK(N-1)
          DO 3000 I=1,NX1
            DO 3000 J=1,NY
              C(I1+(I-1)*LJ1+(J-1)*NX1*LJ1)=WK(J+(I-1)*LJ+(I1-1)*NX1*LJ)
3000    CONTINUE
        NX=NX1
        JU=N-1
      ELSE
!	Set up the RHS for interpolation along N th dimension
        LJ=NK(N)
        DO 3200 I=1,NX
          DO 3200 J=1,NK(N)
            C(J+(I-1)*LJ)=F(I+(J-1)*NX)
3200    CONTINUE
 
      ENDIF
 
!	Loop for interpolation in each dimension, each pass
!	interpolates along 2 dimensions
      DO 5000 J1=JU,1,-2
        NUM=NX*NY
        IFLG1=2
        LJ=NK(J1)
        CALL GAUBND(NK(J1),KB,NUM,AX(1,J1),C,DET,IDET,INTX(1,J1),LJ,
     1         IER,IFLG1,WK(N0))
        IF(IER.GT.100) RETURN
 
!	Set up the RHS for interpolation along the next dimension
        NX1=NX/NK(J1-1)
        NY1=NY*NK(J1)
        LJ1=NK(J1-1)
        DO 3500 I1=1,NY
          DO 3500 I2=1,NK(J1)
            DO 3500 I=1,NK(J1-1)
              DO 3500 J=1,NX1
                WK(I+(J-1)*LJ1+(I2-1)*NX+(I1-1)*NX*NK(J1))=
     1          C(I2+(J-1)*LJ+(I-1)*NX1*LJ+(I1-1)*NX*LJ)
3500    CONTINUE
        NX=NX1
        NY=NY1
        NUM=NY*NX
        IFLG1=2
        LJ=NK(J1-1)
        CALL GAUBND(NK(J1-1),KB,NUM,AX(1,J1-1),WK,DET,IDET,
     1         INTX(1,J1-1),LJ,IER,IFLG1,WK(N0))
        IF(IER.GT.100) RETURN
        IF(J1.EQ.2) THEN
!	Store the coefficients in array C
          DO 3800 I=1,NK(1)
            DO 3800 J=1,NY
              C(I+(J-1)*NK(1))=WK(I+(J-1)*LJ)
3800      CONTINUE
        ELSE
 
!	Set up the RHS for interpolation along the next dimension
          LJ1=NK(J1-2)
          NX1=NX/NK(J1-2)
          NY1=NY*NK(J1-1)
          DO 4000 I1=1,NK(J1-1)
            DO 4000 I2=1,NK(J1-2)
              DO 4000 I=1,NX1
                DO 4000 J=1,NY
                  C(I2+(I-1)*LJ1+(I1-1)*NX+(J-1)*NX*NK(J1-1))=
     1            WK(I1+(I-1)*LJ+(I2-1)*LJ*NX1+(J-1)*LJ*NX)
4000      CONTINUE
          NX=NX1
          NY=NY1
        ENDIF
5000  CONTINUE
      END
