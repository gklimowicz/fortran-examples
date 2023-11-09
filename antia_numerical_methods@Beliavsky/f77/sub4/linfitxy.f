!       Least squares straight line fit when there is error in both x and y
!       To fit equation of form y=a+b*x
!
!       N : (input) Number of data points
!       X : (input) Real array of length N containing the x values
!       Y : (input) Real array of length N containing the y values
!       SIGX : (input) estimated error in x values, assumed to be the
!               same for all points
!       SIGY : (input) estimated error in y values, assumed to be the
!               same for all points
!       RHO : (input) estimated correlation between errors in x and y
!               assumed to be the same for all points
!       XI : (output) Reall array of length N which will give the fitted x values
!       YI : (output) Reall array of length N which will give the fitted y values
!       A : (output) fitted value of the intercept
!       B : (output) fitted value of the slope
!       CHI : (output) value of chi^2 at the minimum
!       IER : (output) Error parameter, IER=0 implies successful execution
!               IER=617 implies that discriminant of quadratic equation
!                       is negative and calculations are aborted
!
!       Required routines : none


      SUBROUTINE LINFITXY(N,X,Y,SIGX,SIGY,RHO,XI,YI,A,B,CHI,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),Y(N),XI(N),YI(N)

      R=RHO*SIGY/SIGX
      S=(SIGY/SIGX)**2
      SX=0.0
      SY=0.0
      SXY=0.0
      SXX=0.0
      SYY=0.0
      DO I=1,N
        SX=SX+X(I)
        SY=SY+Y(I)
        SXY=SXY+X(I)*Y(I)
        SXX=SXX+X(I)*X(I)
        SYY=SYY+Y(I)*Y(I)
      enddo
      SX=SX/N
      SY=SY/N
      SXY=SXY/N
      SXX=SXX/N
      SYY=SYY/N
      IER=0
      
!       Find the quadratic in the slope and solve it
      C2=SIGX*SIGX*(SXY-SX*SY)+RHO*SIGX*SIGY*(SX*SX-SXX)
      C1=SIGX*SIGX*(SY*SY-SYY)+SIGY*SIGY*(SXX-SX*SX)
      C0=-RHO*SIGX*SIGY*(SY*SY-SYY)-SIGY*SIGY*(SXY-SX*SY)
      DEL=C1*C1-4*C0*C2
      IF(DEL.LT.0.0) THEN
        IER=617
        RETURN
      ENDIF
      DEL=SQRT(DEL)
      IF(C1.GT.0.0) THEN
        B1=(-C1-DEL)/(2*C2)
      ELSE
        B1=(-C1+DEL)/(2*C2)
      ENDIF
      B=C0/(C2*B1)
      A1=SY-B1*SX
      A=SY-B*SX
      T1=(B1*R-S)/(B1-R)
      T=(B*R-S)/(B-R)

!       Choose the solution with smaller chi^2
      CHI1=0.0
      CHI=0.0
      DO I=1,N
        X1=(T1*X(I)-Y(I)+A1)/(T1-B1)
        X2=(T*X(I)-Y(I)+A)/(T-B)
        Y1=A1+B1*X1
        Y2=A+B*X2
        CHI1=CHI1+((X(I)-X1)/SIGX)**2-2*RHO*(X(I)-X1)*(Y(I)-Y1)
     1          /(SIGX*SIGY) + ((Y(I)-Y1)/SIGY)**2
        CHI=CHI+((X(I)-X2)/SIGX)**2-2*RHO*(X(I)-X2)*(Y(I)-Y2)
     1          /(SIGX*SIGY) + ((Y(I)-Y2)/SIGY)**2
      ENDDO
      
      IF(CHI1.LT.CHI) THEN
        CHI=CHI1
        A=A1
        B=B1
        T=T1
      ENDIF

      DO I=1,N
        XI(I)=(T*X(I)-Y(I)+A)/(T-B)
        YI(I)=A+B*XI(I)
      ENDDO
      END
