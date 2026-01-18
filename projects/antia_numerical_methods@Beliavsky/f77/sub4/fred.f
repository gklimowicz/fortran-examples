!	To solve linear Fredholm equation using quadrature method
!
!	M : (input) Number of abscissas to be used in quadrature formula
!	A : (input) Lower limit of the integral
!	B : (input) Upper limit of the integral
!	WT : (input/output) Real array of length M containing the
!		weights used in quadrature formula.
!		If IQ is negative the weights must be supplied by the user
!		otherwise they are calculated by the routine
!	X : (input/output) Real array of length M containing the
!		abscissas used in quadrature formula.
!		If IQ is negative the abscissas must be supplied by the user
!		otherwise they are calculated by the routine
!	F : (output) Real array of length M containing the calculated solution
!		F(I) is the value of solution at X(I).
!	FC : (output) Real array of length M containing the calculated solution
!		after applying the deferred correction. This will be
!		relevant only if IQ=1 and IT=1,2. In other cases a dummy
!		array of any length may be supplied.
!	FG : (input) Name of the function routine used to calculate the right
!		hand side g(X). For IT=3 g(x) is not required, but in that
!		case this function is used to calculate an initial guess
!		for the eigenfunctions. In most cases the inverse iteration
!		converges from essentially arbitrary initial guess and it
!		is enough to set FG(X) to any nonzero value.
!	FKER : (input) Name of the function routine used to calculate the
!		kernel K(x,t)
!	EI : (input/output) Initial guess for the eigenvalue. After execution
!		it will contain the computed eigenvalue. Used only for IT=3
!	WK : Real array used as scratch space. Its length should be M*M for
!		IT=1,2, while for IT=3 it should be 2M*M+M.
!	IWK : Integer array of length M used as a scratch space
!	IQ : (input) Integer variable to specify the quadrature formula to be used.
!		If IQ=1, then trapezoidal rule is used and deferred correction
!			is calculated using Gregory's formula
!		If IQ=2, the Simpson's 1/3 rule is used
!		If IQ=4,8,16,32 a composite rule using IQ point
!			Gauss-Legendre formula is used.
!		In all these cases the weights and abscissas are calculated
!		If IQ is negative then it is assumed that weights and abscissas
!		are supplied in arrays WT and X.
!		Other values of IQ will cause an error return.
!	IT : (input) Integer variable to specify the type of integral equation
!		If IT=1 Fredholm equation of the first kind is solved
!		If IT=2 Fredholm equation of the second kind is solved
!		If IT=3 Fredholm equation of the third kind (eigenvalue
!			problem) is solved
!	REPS : (input) Required relative accuracy in calculating eigenvalue
!		and eigenvectors. It is not used for IT=1,2. It only specifies
!		the accuracy to which inverse iteration converges. It does not
!		control the truncation error.
!	IER : (output) The error parameter, IER=0 implies successful execution
!		IER=-11 implies that the calculations are actually performed
!			using a smaller number of abscissas than M
!		IER=706 implies M<3, IT>3, IT.LE.0 or M is not sufficient
!			to apply the required quadrature formula. No calculations
!			are done.
!		IER=707 implies that IQ is not acceptable and no calculations
!			are done
!		Other values of IER may be set by GAUELM and INVIT
!
!	FUNCTION FG(X) and FUNCTION FKER(X,T) must be supplied by the user.
!		FG is the right hand side function g(x) and FKER(X,T) is 
!		the kernel. The integral equation is specified in the form
!		given by Eq.(12.1)
!
!	Required routines : GAUELM, INVIT, FG, FKER
!	
      SUBROUTINE FRED(M,A,B,WT,X,F,FC,FG,FKER,EI,WK,IWK,IQ,IT,REPS,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WT(M),F(M),FC(*),X(M),WK(M,*),IWK(M)
      DIMENSION W(31),XA(31)

!	Weights and abscissas for Gauss-Legendre quadrature.
!	For N-point formula W(K)=W(N-K+1) and XA(K)=-XA(N-K+1)
!		For K=1,2,...,N/2. Hence only half points are tabulated.
!	For 2-point W(1); 4-point W(2), W(3); 8-point W(4),...,W(7);
!	16-point W(8),...,W(15); 32-point W(16),...,W(31) are the
!	weights corresponding to abscissas XA(I).

      DATA W/1.0D0,
     1       0.34785484513745385737D0, 0.65214515486254614263D0,
     2       0.10122853629037625915D0, 0.22238103445337447054D0,
     3       0.31370664587788728734D0, 0.36268378337836198297D0,
     4       0.02715245941175409485D0, 0.06225352393864789286D0,
     5       0.09515851168249278481D0, 0.12462897125553387205D0,
     6       0.14959598881657673208D0, 0.16915651939500253819D0,
     7       0.18260341504492358887D0, 0.18945061045506849629D0,
     8       0.00701861000947009660D0, 0.01627439473090567061D0,
     9       0.02539206530926205945D0, 0.03427386291302143310D0,
     1       0.04283589802222668066D0, 0.05099805926237617620D0,
     2       0.05868409347853554714D0, 0.06582222277636184684D0,
     3       0.07234579410884850623D0, 0.07819389578707030647D0,
     4       0.08331192422694675522D0, 0.08765209300440381114D0,
     4       0.09117387869576388471D0, 0.09384439908080456564D0,
     5       0.09563872007927485942D0, 0.09654008851472780057D0/

      DATA XA/0.57735026918962576451D0,
     1        0.86113631159405257522D0, 0.33998104358485626480D0,
     2        0.96028985649753623168D0, 0.79666647741362673959D0,
     3        0.52553240991632898582D0, 0.18343464249564980494D0,
     4        0.98940093499164993260D0, 0.94457502307323257608D0,
     5        0.86563120238783174388D0, 0.75540440835500303390D0,
     6        0.61787624440264374845D0, 0.45801677765722738634D0,
     7        0.28160355077925891323D0, 0.09501250983763744019D0,
     8        0.99726386184948156354D0, 0.98561151154526833540D0,
     9        0.96476225558750643077D0, 0.93490607593773968917D0,
     1        0.89632115576605212397D0, 0.84936761373256997013D0,
     2        0.79448379596794240696D0, 0.73218211874028968039D0,
     3        0.66304426693021520098D0, 0.58771575724076232904D0,
     4        0.50689990893222939002D0, 0.42135127613063534536D0,
     5        0.33186860228212764978D0, 0.23928736225213707454D0,
     6        0.14447196158279649349D0, 0.04830766568773831623D0/

      IER=0
      IF(M.LT.3.OR.IT.GT.3.OR.IT.LE.0) THEN
        IER=706
        RETURN
      ENDIF

!	M should not be changed since it is used in dimension statement
      N=M
      IF(IQ.EQ.1) THEN
!	Use the trapezoidal rule
        H=(B-A)/(N-1)
        WT(1)=0.5*H
        WT(N)=WT(1)
        X(1)=A
        X(N)=B
        DO 2000 I=2,N-1
          X(I)=A+(I-1)*H
2000    WT(I)=H

      ELSE IF(IQ.EQ.2) THEN
!	Use the Simpson's 1/3 rule, if N is even, then reduce it by 1
        N=2*((N-1)/2)+1
        H=(B-A)/(N-1)
        WT(1)=H/3.
        X(1)=A
        DO 2100 I=2,N-1,2
          X(I)=A+(I-1)*H
          X(I+1)=A+I*H
          WT(I)=4.*H/3.
2100    WT(I+1)=2.*H/3.
        WT(N)=WT(1)
        X(N)=B

      ELSE IF(IQ.GE.0) THEN
!	Try Gauss-Legendre formulas
        NO=-1
        IF(IQ.EQ.4) NO=1
        IF(IQ.EQ.8) NO=3
        IF(IQ.EQ.16) NO=7
        IF(IQ.EQ.32) NO=15
        IF(NO.LT.0) THEN
!	If IQ is not acceptable then quit
          IER=707
          RETURN
        ENDIF
        N=(N/IQ)*IQ
        IF(N.LT.IQ) THEN
!	If the number of points is not sufficient, then quit
          IER=706
          RETURN
        ENDIF

!	Setup the weights and abscissas for Gauss-Legendre formula
        H=IQ*(B-A)/N
        DO 2300 I=1,N,IQ
          A1=A+(I-1)*H/IQ+H/2.
          DO 2300 I1=1,IQ/2
            WT(I+I1-1)=W(NO+I1)*H/2.
            WT(I+IQ-I1)=WT(I+I1-1)
            X(I+I1-1)=A1-XA(NO+I1)*H/2.
            X(I+IQ-I1)=A1+XA(NO+I1)*H/2.
2300    CONTINUE
      ENDIF
      IF(M.NE.N) IER=-11

!	Setting up the equation matrix and the right hand side
      DO 3000 I=1,N
        F(I)=FG(X(I))
        DO 2800 J=1,N
          WK(J,I)=WT(I)*FKER(X(J),X(I))
2800    CONTINUE
        IF(IT.EQ.2) WK(I,I)=WK(I,I)-1.
3000  CONTINUE

      NUM=1
      LJ=M
      IF(IT.LE.2) THEN
!	Solve the system of linear equations
        IFLG=0
        CALL GAUELM(N,NUM,WK,F,DET,IWK,LJ,IER,IFLG)
        IF(IER.GT.0) RETURN
      ELSE
!	Solve the eigenvalue problem
        IFLG=0
        P=EI
        NIT=0
        CALL INVIT(WK,N,LJ,P,F,IFLG,EI,RC,REPS,WK(1,N+1),IWK,NIT,IER)
      ENDIF

      IF(IQ.EQ.1.AND.IT.NE.3) THEN
!	Apply the deferred correction
        DO 3200 I=1,N
          U2=FKER(X(I),X(2))*F(2)
          U2M=FKER(X(I),X(N-1))*F(N-1)
          T1=FKER(X(I),X(1))*F(1)-U2+FKER(X(I),X(N))*F(N)-U2M
          T2=T1+FKER(X(I),X(3))*F(3)-U2+FKER(X(I),X(N-2))*F(N-2)-U2M
          FC(I)=H*T1/12.+H*T2/24.
3200    CONTINUE
        CALL GAUELM(N,NUM,WK,FC,DET,IWK,LJ,IER,IFLG)

        DO 4000 I=1,N
4000    FC(I)=F(I)+FC(I)
      ENDIF
      END
