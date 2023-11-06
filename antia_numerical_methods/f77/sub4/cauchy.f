!	To evaluate the Cauchy principal value of an integral over a finite interval
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit (B > A)
!	C : (input) Location of the singularity (A < C < B)
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	F : (input) Name of the function routine to calculate the integrand
!	FUNP : (input) Name of the function routine to calculate F(C+x)+F(C-x)
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=304 implies A>B, A>C or C>B in which case no calculations
!			are done
!		Other values may be set by subroutine ADPINT which is called
!		twice. The returned value of IER is IER1+IER2*2, where IER1
!		and IER2 are the returned values of IER from two calls to ADPINT
!		In this case DIF will contain the estimated accuracy
!	NPT : (output) Number of function evaluations used by subroutine
!
!		FUNCTION F(X) must be supplied by the user.
!		FUNCTION FUNP(X) to calculate F(C+X)+F(C-X) should also be
!		supplied by the user. The value of C is available through
!		common block. Simplest version for FUNP may be
!
!		FUNCTION FUNP(X)
!		IMPLICIT REAL*4(A-H,O-Z)
!		COMMON/CAUFN/C
!		FUNP=F(C+X)+F(C-X)
!		END
!
!		If F(C+X)+F(C-X) can be combined roundoff error may be reduced.
!		There is no provision to pass the name F to FUNP, so it
!		will have to put explicitly.
!
!	Required routines : ADPINT, KRONRD, F, FUNP
 
      SUBROUTINE CAUCHY(RI,A,B,C,REPS,AEPS,DIF,F,FUNP,IER,NPT)
!      IMPLICIT REAL*8(A-H,O-Z)
!	To pass the value of C to FUNCTION F or FUNP
      EXTERNAL F,FUNP
      COMMON/CAUFN/CC
 
      IF(A.GT.B.OR.A.GT.C.OR.C.GT.B) THEN
        IER=304
        RETURN
      ENDIF
 
!     FIRST EVALUATE THE SINGULAR PART
 
      CC=C
      R=MIN(C-A,B-C)
      AA=0.0
      NMAX=0
      NPT1=0
      DIF1=0.0
      RI1=0.0
      CALL ADPINT(RI1,AA,R,REPS,AEPS,DIF1,FUNP,IER,NPT1,NMAX)
 
!     EVALUATE THE REMAINING PORTION
 
      IF(C-A.GT.B-C) THEN
        AA=A
        BB=C-R
      ELSE
        AA=C+R
        BB=B
      ENDIF
      NPT2=0
      DIF2=0.0
      RI2=0.0
      IER1=0
      IF(ABS(BB-AA).GT.AEPS)
     1  CALL ADPINT(RI2,AA,BB,REPS,AEPS,DIF2,F,IER1,NPT2,NMAX)
 
      RI=RI1+RI2
!     FUNCTION FUNP REQUIRES TWO EVALUATIONS OF FUN
      NPT=2*NPT1+NPT2
      DIF=ABS(DIF1)+ABS(DIF2)
      IER=IER+2*IER1
      END

!	---------------------------

!		FUNCTION FUNP(X)
!		IMPLICIT REAL*4(A-H,O-Z)
!		COMMON/CAUFN/C
!		FUNP=F(C+X)+F(C-X)
!		END
