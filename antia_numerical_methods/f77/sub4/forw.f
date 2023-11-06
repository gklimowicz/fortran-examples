!	To solve the forward problem corresponding to a linear inverse
!	problem. This routine may be used to generate artificial data
!	for testing inversion techniques.
!
!	NP : (input) Number of points used in defining the kernels.
!	NM : (input) Number of data points in the inversion problem
!		which should be same as the number of kernels that are
!		supplied in array RKER.
!	R : (input) Real array of length NP containing the coordinates
!		of points at which kernels are available.
!	RKER : (input) Real array of length IK*NP containing the kernels
!		for the inverse problem. RKER(I,J) should contain the
!		value at R(J) for the Ith kernel.
!	IK : (input) The first dimension of RKER as declared in the calling
!		program
!	DI : (output) Real array of length NM containing the calculated
!		data points using the kernel.
!	F : (input/output) Real array of length NP containing the function
!		value at points in R. F(I) should contain the function
!		value at R(I). If IFLG=0, the function values are
!		calculated using user supplied function routine FUN,
!		otherwise, these values must be supplied while calling
!		the routine.
!	FUN : (input) Name of function routine to calculate the given
!		function. This is used only if IFLG=0, otherwise the
!		function values are to be supplied in array F.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=711 implies that IK<NM and no calculations are done
!	IFLG : (input/output) Integer parameter used as a flag to decide
!		the type of computation required.
!		If IFLG=0, then the function values are calculated using
!			a user supplied routine FUN. These values are stored
!			in array F and IFLG is set to 1 so that next time
!			the values need not be calculated.
!		For other values of IFLG the function values must be
!			supplied in array F.
!
!	FUNCTION FUN(X) must be supplied by the user
!
!	Required routines : FUN

      SUBROUTINE FORW(NP,NM,R,RKER,IK,DI,F,FUN,IER,IFLG)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(NP),RKER(IK,NP),DI(NM),F(NP)
 
      IF(IK.LT.NM) THEN
        IER=711
        RETURN
      ENDIF
      IER=0

      IF(IFLG.EQ.0) THEN
!     Calculate the function value using supplied routine
        DO 2000 I=1,NP
          F(I)=FUN(R(I))
2000    CONTINUE
        IFLG=1
      ENDIF
 
!     Calculate the integrals
      DO 3000 I=1,NM
        S1=0.0
        H=(R(2)-R(1))/2.
        DO 2500 IR=1,NP
          S1=S1+H*F(IR)*RKER(I,IR)
          IF(IR.LT.NP-1) THEN
            H=(R(IR+2)-R(IR))/2.0
          ELSE IF(IR.EQ.NP-1) THEN
            H=(R(IR+1)-R(IR))/2.0
          ENDIF
2500    CONTINUE
        DI(I)=S1
3000  CONTINUE
 
 
      END
