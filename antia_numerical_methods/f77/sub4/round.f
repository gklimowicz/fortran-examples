!	To round a number to N digits using base B.
!
!	X : (input) The number to be rounded
!	N : (input) The number of digits required in rounded number
!	B : (input) Base of number system to be used for rounding
!	The rounded value will be returned by ROUND
!
!	Required routines : None

      FUNCTION ROUND(X,N,B)
!      IMPLICIT REAL*8(A-H,O-Z)

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
