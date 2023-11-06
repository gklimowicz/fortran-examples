!     PROGRAM FOR ROUNDING A GIVEN FLOATING POINT NUMBER TO SPECIFIED
!     NO. OF DIGITS USING THE SPECIFIED BASE. ALL NUMBERS ARE PRINTED
!     OUT IN DECIMAL SYSTEM AND HENCE IT MAY BE DIFFICULT TO VERIFY THE
!     RESULT UNLESS THE BASE B=10.

      PROGRAM FROUND
      IMPLICIT REAL*16(A-H,O-Z)

51    FORMAT(1PD15.6,'  ROUNDED TO',D15.6)

100   PRINT *,'TYPE X=NO. TO BE ROUNDED,   N=NO. OF DIGITS,   B=BASE'
      PRINT *,'        (QUITS WHEN N.LE.0)'
      READ *,X,N,B
      IF(N.LE.0) STOP
      XR=ROUND(X,N,B)
      WRITE(6,51) X,XR
      GO TO 100
      END
 
!     -----------------------------------------------------
 
!	To round a number to N digits using base B.
!
!	X : (input) The number to be rounded
!	N : (input) The number of digits required in rounded number
!	B : (input) Base of number system to be used for rounding
!	The rounded value will be returned by ROUND
!
!	Required routines : None

      FUNCTION ROUND(X,N,B)
      IMPLICIT REAL*16(A-H,O-Z)

      ROUND=0.0
      IF(X.EQ.0.0) RETURN
      XA=ABS(X)
      LGX=LOG(XA)/LOG(B)
      IF(XA.LT.1.0) LGX=LGX-1
      FX=B**(N-LGX-1)
      NX=XA*FX+0.5
      XA=NX/FX
      ROUND=SIGN(XA,X)
      END
