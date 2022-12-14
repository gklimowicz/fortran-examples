      MODULE lookup_emis60ms_com

      IMPLICIT NONE

      INTEGER,PARAMETER :: Im=234,Jm=234,Km=22

      INTEGER :: w,count,i,j,k,bb
      REAL(KIND=8) :: deltaw,sum,a,b,bessi0,zsum
      REAL(KIND=8) :: s,xx,bys2,val,array(im,jm)
      REAL(KIND=8) :: tke(Jm),wb(Im),wt(Km)

      END MODULE lookup_emis60ms_com

      PROGRAM lookup_emis60ms

c**** Calculates lookup table for logarithm of soil dust aerosol
c**** emission for a range of index values of the surface wind speed
c**** (wsgcm), standard deviation of the surface wind speed (sigma) due
c**** to subgrid variability, and surface wind speed threshold of dust
c**** emission (soilvtrsh), which depends on soil moisture. In ModelE,
c**** the logarithm of dust emission is obtained by interpolating the
c**** logarithm of the emission from the lookup table.

c**** The parameters Im, Jm, and Km must be the same as the ones in the
c**** model defined in TRDUST_COM.f. The current parameters in this
c**** program are set to calculate the logarithm of dust emission for a
c**** range of wsgcm, sigma, and soilvtrsh up to around 60 m/s, 60 m/s,
c**** and 14 m/s, respectively. To increase the range increase Im, Jm,
c**** and Km in the module above, re-calculate the lookup table, and
c**** adjust the parameters in the model (Lim, Ljm, Lkm in
c**** TRDUST_COM.f), accordingly.

      USE lookup_emis60ms_com

      IMPLICIT NONE

      REAL(KIND=8) :: func

      OPEN(10,FILE='log_emis60ms.new',FORM='unformatted')

      DO k=1,Km
        wt(k)=6.D0+0.5D0*k
      END DO

      zsum=0.D0
      DO j=1,Jm
        IF (j <= 30) THEN
          zsum=zsum+0.0001D0+FLOAT(j-1)*0.00008D0
          tke(j)=zsum
        ELSE IF (j > 30) THEN
          zsum=zsum-0.055254D0+0.005471D0*DBLE(j)-
     &         1.938365D-4*DBLE(j)**2.D0+
     &         3.109634D-6*DBLE(j)**3.D0-
     &         2.126684D-8*DBLE(j)**4.D0+
     &         5.128648D-11*DBLE(j)**5.D0
          tke(j)=zsum
        END IF
      END DO
      
      zsum=0.D0
      DO i=1,Im
        IF (i <= 30) THEN
          zsum=zsum+0.0001D0+DBLE(i-1)*0.00008D0
          wb(i)=zsum
        ELSE IF (i > 30) THEN
          zsum=zsum-0.055254D0+0.005471D0*DBLE(i)-
     &         1.938365D-4*DBLE(i)**2.D0+
     &         3.109634D-6*DBLE(i)**3.D0-
     &         2.126684D-8*DBLE(i)**4.D0+
     &         5.128648D-11*DBLE(i)**5.D0
          wb(i)=zsum
        END IF
      END DO
      
      DO k=1,Km
        DO j=1,Jm
          DO i=1,Im

            sum=0.D0
            bys2=(1.D0/tke(j))**2.D0
            deltaw=0.0001D0
            
            w=INT(wt(k)/deltaw)
            xx=deltaw*w+deltaw/2.0D0

            b=wb(i)+30.D0*tke(j)
            bb=INT(b/deltaw)
            DO WHILE(w <= bb)
               sum=sum+deltaw*func(i,j,k,xx)
               w=w+1
               xx=deltaw*w+deltaw/2.0D0
            END DO
            IF (sum == 0.D0) THEN
               array(i,j)=-1000.D0
            ELSE
               array(i,j)=LOG(sum)
            END IF
c            WRITE(*,*) 'i,j,k,wb(i),tke(j),wt(k),pdfint:',i,j,k,wb(i),
c     &           tke(j),wt(k),array(i,j)
         END DO
      END DO
      WRITE(10) array
      END DO
      
      END PROGRAM lookup_emis60ms

      REAL(KIND=8) FUNCTION func(i,j,k,w)

      USE lookup_emis60ms_com,ONLY : bys2,wb,wt

      IMPLICIT NONE

      REAL(KIND=8),PARAMETER :: Pi=4.D0*ATAN(1.D0)
      INTEGER,INTENT(IN) :: i,j,k
      REAL(KIND=8),INTENT(IN) :: w

      REAL(KIND=8) :: bessi0

      IF ((w*wb(i)*bys2) < 100.D0) THEN
         func=(bys2*w**3.D0)*(w-wt(k))
     &       *EXP((bys2/2.D0)*(-w**2.D0-wb(i)**2.D0))
     &       *bessi0(w*wb(i)*bys2)
      ELSE
         func=(bys2*w**3.D0)*(w-wt(k))
     &       *EXP(-(bys2/2.D0)*(w-wb(i))**2.D0)
     &       *(1.D0/SQRT(2.D0*Pi*w*wb(i)*bys2))
      END IF

      RETURN
      END FUNCTION func

      FUNCTION bessi0(x)
      REAL(KIND=8) :: bessi0,x
      REAL(KIND=8) :: ax
      REAL(KIND=8) :: p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229D0,3.0899424D0,
     &     1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228D0,0.1328592D-1
     &     ,0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1
     &     ,0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(x) < 3.75D0) THEN
        y=(x/3.75D0)**2.D0
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      ELSE
        ax=ABS(x)
        y=3.75D0/ax
        bessi0=(EXP(ax)/SQRT(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y
     &       *(q7+y*(q8+y*q9))))))))
      END IF
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.
