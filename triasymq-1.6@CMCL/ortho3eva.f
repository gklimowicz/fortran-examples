cc Copyright (C) 2009: Zydrunas Gimbutas and Hong Xiao
cc Contact: Zydrunas Gimbutas <gimbutas@cims.nyu.edu>
cc          Hong Xiao <hxiao@ucdavis.edu>
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       Orthogonal polynomials on the standard tetrahedron
c       
c       (-1,-1/Sqrt[3],-1/Sqrt[6]), (0,2/Sqrt[3],-1/Sqrt[6]), 
c       (1,-1/Sqrt[3],-1/Sqrt[6]),  and (0,0,3/Sqrt[6])
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       This file contains a set of subroutines for the handling 
c       of the evaluation of orthogonal polynomials on the standard 
c       tetrahedron.  It contains 2 user-callable subroutines.
c
c   ortho3eva - evaluates at the user-supplied point w
c       a collection of polynomials (of x,y,z) orthogonal on the
c       standard tetrahedron
c
c   ortho3eva4 - evaluates at the user-supplied point w a collection
c       of polynomials (of x,y,z) orthogonal on the standard tetrahedron,
c       together with their derivatives with respect to x, y, and z.
c
c
c     NOTE: experimental code, needs better memory management!!!
c
c
c
        subroutine ortho3eva(maxorder,xyz,fvals,ncount)
c
        implicit none
c
        real *8 xyz(3)
c
        real *8 fvals(1),derxs(1),derys(1),derzs(1),
     1     f2s(1400 00),f3s(1400 00)
c
        real *8 f1(1400 00),
     1     dersx1(1400 00),dersy1(1400 00),
     2     dxs(1400 00), dys(1400 00)
c
        real *8 dummy1(1400 00),dummy2(1400 00)

        integer maxorder,ncount,needle,mp,mmax,n,k,m
        real *8 zold,yold,xold,zero,done,s6,s3,x,y,z,u,v,w,
     1     par1,par2,par5,xa,yb,par7,fval,derx,par6,q1,q2,q3,dery,
     2     derz,rr,scale,derynew,derznew

c
c       This subroutine evaluates the Koornwinder's 
c       orthogonal polynomial up to order maxorder
c       on the standard tetrahedron with the  vertices
c
c          (-1,-1/Sqrt[3],-1/Sqrt[6]), (0,2/Sqrt[3],-1/Sqrt[6]), 
c          (1,-1/Sqrt[3],-1/Sqrt[6]),  and
c          (0,0,3/Sqrt[6])
c
c       at the user specified point xyz=(xold,yold,zold),  together 
c       with their derivatives with respect to xold, yold, and zold.  
c       The polynomials are ordered by their order, and in each order,
c       they are ordered lexicographically in terms of their indices
c       (m,n,k).   The total number of polynomials should be
c       (6+11m+6 m^2 + m^3)/6.
c       
c
c       This subroutine is based on the Koornwinder's representation
c       of the orthogonal polynomials on the right triangle
c
c       (-1,-1,-1), (-1,1,-1), (1,-1,-1),(-1,-1,1)                    (2)
c       which is given by the formula
c
c       K_mnk(x,y,z) = 
c            P_m ((2x+2+y+z)/(-y-z)) * ((-y-z)/2)^m *
C                  P_n^{2m+1,0}((2y+1+z)/(1-z)) * ((1-z)/2)^{n}
C                  P_k^{2m+2n+2,0} (z) 
c       with (x,y,z) transformed as
c       x = -1/2 + xold -yold/s3 - zold/s6
c       y = -1/2 + 2/s3 * yold - zold/s6
c       z = -1/2 + s6/2 * zold
c       and
c       s3=sqrt(3)
c       s6=sqrt(6). 
c       Futhermore, P_m are the Legendge polynomials of order m
c       and P_n^{2m+1,0} are the Jacobi polynomials of order n with
c       the parameters alpha=2*m+1 and beta=0, etc. 
c
C
c
c                   Input parameters:
c
c  mmax - the maximum order to which the polynomials are to be evaluated;
c  z - the location in R^3 where the polynomials are to be evaluated;
c       normally, expected to be inside (including boundary) the 
c       standard triangle (1) above.
c  w - the work array. must be at least mmax+1+(mmax+1)**2 long
c     
c                   Output parameters:
c
c  fval - the orthogonal polynomial evaluated at the point z 
c  derx - the derivative with respect to x of the polynomial
c  dery - the partial derivative with respect to y of the polynomial
c  derz - the derivative with respect to y of the polynomial
c

c       . . . allocate memory for the work arrays
c
c       
c
        xold=xyz(1)
        yold=xyz(2)
        zold=xyz(3)

c       
c       convert the coordinates onto the right tetrahedron
c       

ccc        call prinf("evas: just entering *",1,1)

        zero=0
        done=1
        s6=sqrt(6*done)
        s3=sqrt(3*done)
c
          x=-done/2 + xold - yold/s3 - zold/s6 
          y=-done/2 + 2*yold/s3 - zold/s6 
          z=-done/2 + zold*s6/2 
c
           u=(2*x+2+y+z)/(-y-z)
           v=(2*y+1+z)/(1-z)
           w=z
CCC           
        
        par1=(2*x+2+y+z)/2
        par2=(-y-z)/2

ccc        call prinf("evas: before klegeypols *",2,1)

c       
c       gives one extra function values, since the order starts at 0.
c       
        call klegeypols3(par1,par2,maxorder,f1,dersx1,dersy1)

c       
c       
c       get f2s, (maxorder+1)**2 of them
c
        par5=(2*y+1+z)/(1-z)

        needle=1

        do 3400 m=0,maxorder

        par6=2*m+1


        xa=(2*y+1+z)/2.0d0
        yb=(1-z)/2.0d0


        call kjacoypols3(xa,yb,par6,zero,(maxorder-m),f2s(needle),
     1     dxs(needle),dys(needle))

        needle=needle+maxorder+1

ccc        call prin2("f2=*",f2s,(maxorder+1)*(maxorder+1))

c
 3400   continue
c

ccc        call prinf("evas: after f2s=*",needle,1)

c
c       get f3s, (maxorder+1)**2 of them
c
        needle=1
        do 3600 mp=2,2*maxorder+2,2
c
c
        par7=mp

        call kjacoypols3(z,done,par7,zero,(maxorder+1-mp/2),
     1     f3s(needle),dummy1(needle),dummy2)

        needle=needle+maxorder+1

ccc        call prin2("f3=*",f3s,(maxorder+1)*(maxorder+1))

 3600   continue

ccc        call prinf("evas: after f3s=*",needle,1)
c
c
        ncount=0

        do 6000 mmax=0,maxorder
c
        do 5000 m=0,mmax
c
        do 4000 n=0,mmax-m
c
        k=mmax-m-n
        ncount=ncount+1


        fval=f1(m+1)*f2s(n+1 + m*(maxorder+1))*
     1     f3s(k+1 + (m+n)*(maxorder+1))

c
c       partial derivative for x
c
        derx = dersx1(m+1) * f2s(n+1 + m*(maxorder+1)) * 
     1     f3s(k+1 + (m+n)*(maxorder+1))


c
c       partial derivative for y
c


        q1=(0.5d0 * dersx1(m+1) -0.5d0*dersy1(m+1))*
     1     f2s(n+1 +  m*(maxorder+1)) * f3s(k+1  + (m+n)*(maxorder+1))

        q2=dxs(n+1  +  m*(maxorder+1))*f1(m+1) *
     1     f3s(k+1  + (m+n)*(maxorder+1))

        dery=q1+q2

c
c       partial derivative for z
c

        q1=(0.5d0 * dersx1(m+1) - 0.5d0*dersy1(m+1))*
     1     f2s(n+1 +  m*(maxorder+1) )* f3s(k+1  + (m+n)*(maxorder+1) )

        q2=(0.5d0*dxs(n+1  +  m*(maxorder+1)) - 0.5d0 *
     1     dys(n+1  +  m*(maxorder+1)) ) *f1(m+1)*
     1     f3s(k+1  + (m+n)*(maxorder+1))
        
        q3=f1(m+1)*f2s(n+1  +  m*(maxorder+1)) * 
     1     dummy1(k+1 + (m+n)*(maxorder+1))

        derz=q1+q2+q3


c
c       scaling 
c
	rr= 4.0d0/3
	
ccc        mmax=m+n+k

ccc	scale=sqrt(rr/(done*2*mmax/3+1)/((2.0d0*m+1)*(n+m+1)))
	scale=sqrt(rr/(done*2*mmax/3+1)/((2.0d0*m+1)*(n+m+1))
     1     /sqrt(2.0d0))

        fvals(ncount)=fval/scale

 4000   continue

 5000   continue

 6000   continue

ccc        call prinf("evas: after loops, ncount=*",ncount,1)

        return
        end
c       
c       
c
c
c
        subroutine ortho3eva4(maxorder,xyz,fvals,
     1     derxs,derys,derzs,ncount)

        implicit none
c
        real *8 xyz(3)
c
        real *8 fvals(1),derxs(1),derys(1),derzs(1),
     1     f2s(1400 00),f3s(1400 00)
c
        real *8 f1(1400 00),
     1     dersx1(1400 00),dersy1(1400 00),
     2     dxs(1400 00), dys(1400 00)
c
        real *8 dummy1(1400 00),dummy2(1400 00)

        integer maxorder,ncount,needle,mp,mmax,n,k,m
        real *8 zold,yold,xold,zero,done,s6,s3,x,y,z,u,v,w,
     1     par1,par2,par5,xa,yb,par7,fval,derx,par6,q1,q2,q3,dery,
     2     derz,rr,scale,derynew,derznew

c
c       This subroutine evaluates the Koornwinder's 
c       orthogonal polynomial up to order maxorder
c       on the standard tetrahedron with the  vertices
c
c          (-1,-1/Sqrt[3],-1/Sqrt[6]), (0,2/Sqrt[3],-1/Sqrt[6]), 
c          (1,-1/Sqrt[3],-1/Sqrt[6]),  and
c          (0,0,3/Sqrt[6])
c
c       at the user specified point xyz=(xold,yold,zold),  together 
c       with their derivatives with respect to xold, yold, and zold.  
c       The polynomials are ordered by their order, and in each order,
c       they are ordered lexicographically in terms of their indices
c       (m,n,k).   The total number of polynomials should be
c       (6+11m+6 m^2 + m^3)/6.
c       
c
c       This subroutine is based on the Koornwinder's representation
c       of the orthogonal polynomials on the right triangle
c
c       (-1,-1,-1), (-1,1,-1), (1,-1,-1),(-1,-1,1)                    (2)
c       which is given by the formula
c
c       K_mnk(x,y,z) = 
c            P_m ((2x+2+y+z)/(-y-z)) * ((-y-z)/2)^m *
C                  P_n^{2m+1,0}((2y+1+z)/(1-z)) * ((1-z)/2)^{n}
C                  P_k^{2m+2n+2,0} (z) 
c       with (x,y,z) transformed as
c       x = -1/2 + xold -yold/s3 - zold/s6
c       y = -1/2 + 2/s3 * yold - zold/s6
c       z = -1/2 + s6/2 * zold
c       and
c       s3=sqrt(3)
c       s6=sqrt(6). 
c       Futhermore, P_m are the Legendge polynomials of order m
c       and P_n^{2m+1,0} are the Jacobi polynomials of order n with
c       the parameters alpha=2*m+1 and beta=0, etc. 
c
C
c
c                   Input parameters:
c
c  mmax - the maximum order to which the polynomials are to be evaluated;
c  z - the location in R^3 where the polynomials are to be evaluated;
c       normally, expected to be inside (including boundary) the 
c       standard triangle (1) above.
c  w - the work array. must be at least mmax+1+(mmax+1)**2 long
c     
c                   Output parameters:
c
c  fval - the orthogonal polynomial evaluated at the point z 
c  derx - the derivative with respect to x of the polynomial
c  dery - the partial derivative with respect to y of the polynomial
c  derz - the derivative with respect to y of the polynomial
c

c       . . . allocate memory for the work arrays
c
c       
c
        xold=xyz(1)
        yold=xyz(2)
        zold=xyz(3)

c       
c       convert the coordinates onto the right tetrahedron
c       

ccc        call prinf("evas: just entering *",1,1)

        zero=0
        done=1
        s6=sqrt(6*done)
        s3=sqrt(3*done)
c
          x=-done/2 + xold - yold/s3 - zold/s6 
          y=-done/2 + 2*yold/s3 - zold/s6 
          z=-done/2 + zold*s6/2 
c
           u=(2*x+2+y+z)/(-y-z)
           v=(2*y+1+z)/(1-z)
           w=z
CCC           
        
        par1=(2*x+2+y+z)/2
        par2=(-y-z)/2

ccc        call prinf("evas: before klegeypols *",2,1)

c       
c       gives one extra function values, since the order starts at 0.
c       
        call klegeypols3(par1,par2,maxorder,f1,dersx1,dersy1)

c       
c       
c       get f2s, (maxorder+1)**2 of them
c
        par5=(2*y+1+z)/(1-z)

        needle=1

        do 3400 m=0,maxorder

        par6=2*m+1


        xa=(2*y+1+z)/2.0d0
        yb=(1-z)/2.0d0


        call kjacoypols3(xa,yb,par6,zero,(maxorder-m),f2s(needle),
     1     dxs(needle),dys(needle))

        needle=needle+maxorder+1

ccc        call prin2("f2=*",f2s,(maxorder+1)*(maxorder+1))

c
 3400   continue
c

ccc        call prinf("evas: after f2s=*",needle,1)

c
c       get f3s, (maxorder+1)**2 of them
c
        needle=1
        do 3600 mp=2,2*maxorder+2,2
c
c
        par7=mp

        call kjacoypols3(z,done,par7,zero,(maxorder+1-mp/2),
     1     f3s(needle),dummy1(needle),dummy2)

        needle=needle+maxorder+1

ccc        call prin2("f3=*",f3s,(maxorder+1)*(maxorder+1))

 3600   continue

ccc        call prinf("evas: after f3s=*",needle,1)
c
c
        ncount=0

        do 6000 mmax=0,maxorder
c
        do 5000 m=0,mmax
c
        do 4000 n=0,mmax-m
c
        k=mmax-m-n
        ncount=ncount+1


        fval=f1(m+1)*f2s(n+1 + m*(maxorder+1))*
     1     f3s(k+1 + (m+n)*(maxorder+1))

c
c       partial derivative for x
c
        derx = dersx1(m+1) * f2s(n+1 + m*(maxorder+1)) * 
     1     f3s(k+1 + (m+n)*(maxorder+1))


c
c       partial derivative for y
c


        q1=(0.5d0 * dersx1(m+1) -0.5d0*dersy1(m+1))*
     1     f2s(n+1 +  m*(maxorder+1)) * f3s(k+1  + (m+n)*(maxorder+1))

        q2=dxs(n+1  +  m*(maxorder+1))*f1(m+1) *
     1     f3s(k+1  + (m+n)*(maxorder+1))

        dery=q1+q2

c
c       partial derivative for z
c

        q1=(0.5d0 * dersx1(m+1) - 0.5d0*dersy1(m+1))*
     1     f2s(n+1 +  m*(maxorder+1) )* f3s(k+1  + (m+n)*(maxorder+1) )

        q2=(0.5d0*dxs(n+1  +  m*(maxorder+1)) - 0.5d0 *
     1     dys(n+1  +  m*(maxorder+1)) ) *f1(m+1)*
     1     f3s(k+1  + (m+n)*(maxorder+1))
        
        q3=f1(m+1)*f2s(n+1  +  m*(maxorder+1)) * 
     1     dummy1(k+1 + (m+n)*(maxorder+1))

        derz=q1+q2+q3


c
c       scaling 
c
	rr= 4.0d0/3
	
ccc        mmax=m+n+k

ccc	scale=sqrt(rr/(done*2*mmax/3+1)/((2.0d0*m+1)*(n+m+1)))
	scale=sqrt(rr/(done*2*mmax/3+1)/((2.0d0*m+1)*(n+m+1))
     1     /sqrt(2.0d0))


        fvals(ncount)=fval/scale
        derxs(ncount)=derx/scale
        derys(ncount)=dery/scale
        derzs(ncount)=derz/scale

c
c       the final scaling for the derivatives
c        
        derynew = -1.0d0/s3*derxs(ncount) + 2.0d0/s3*derys(ncount)
        derznew = -1.0d0/s6*derxs(ncount) - 1.0d0/s6*derys(ncount) 
     1     + s6/2.0d0*derzs(ncount)

        derys(ncount)=derynew
        derzs(ncount)=derznew

 4000   continue

 5000   continue

 6000   continue

ccc        call prinf("evas: after loops, ncount=*",ncount,1)

        return
        end
c       
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c       
c       
        subroutine klegeypols3(x,y,n,pols,dersx,dersy)
        implicit real *8 (a-h,o-z)
        dimension pols(1),dersx(1),dersy(1)
c
c     evaluate a sequence of scaled Legendre polynomials P_n(x/y) y^n,
c     with the parameter y \in [0..1], together with their derivatives
c     with respect to the parameters x and y.
c
c       ... if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        dersx(1)=0
        dersy(1)=0
        if(n .eq. 0) return
c
        pols(2)=x
        dersx(2)=1
        dersy(2)=0
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        pols(1)=1
        dersx(1)=0
        dersy(1)=0
        pols(2)=x
        dersx(2)=1
        dersy(2)=0
c
        pkm1=1
        pk=x
        dkm1=0
        dk=1
        ykm1=0
        yk=0
c
        pk=1
        pkp1=x
        dk=0
        dkp1=1
        yk=0
        ykp1=0
c
        do 2000 k=1,n-1
c
        pkm1=pk
        pk=pkp1
        dkm1=dk
        dk=dkp1
        ykm1=yk
        yk=ykp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1*y*y )/(k+1)
        dkp1=( (2*k+1)*(x*dk+pk)-k*dkm1*y*y )/(k+1)
        ykp1=( (2*k+1)*(x*yk)-k*(pkm1*2*y+ykm1*y*y) )/(k+1)
        pols(k+2)=pkp1
        dersx(k+2)=dkp1
        dersy(k+2)=ykp1
 2000 continue
c
        return
        end
c
c
c       
c       
c       
c       
c       
        subroutine kjacoypols3(x,y,a,b,n,pols,dersx,dersy)
        implicit real *8 (a-h,o-z)
        dimension pols(1),dersx(1),dersy(1)
c       
c       evaluates a bunch of Jacobi polynomials multiplied by
c       specific polynomials given by the formula
c
c       P_n^{(a,b)} (x/y) y^n
c
c       at the user-provided point x/y with alpha=a, beta=b.
c       n is the max order  of these Jacobi polynomials. 
c       Partial derivatives of x and y are also returned in
c       the arrays dersx and dersy.   
c

c
c       ... if n=0 or n=1 - exit
c       

        if(n .ge. 2) goto 1200
        pols(1)=1
        dersx(1)=0
        dersy(1)=0
        if(n .eq. 0) return
c       
        pols(2)=(a/2-b/2)*y+(1+a/2+b/2)*x
        dersx(2)=(1+a/2+b/2)
        dersy(2)=(a/2-b/2)
        return
 1200   continue
c       
        pols(1)=1
        dersx(1)=0
        dersy(1)=0

        pols(2)=(a/2-b/2)*y+(1+a/2+b/2)*x
        dersx(2)=(1+a/2+b/2)
        dersy(2)=(a/2-b/2)
c       
c       n is greater than 2. conduct recursion
c       
        pkm1=1
        dxkm1=0
        dykm1=0

        pk=(a/2-b/2)*y+(1+a/2+b/2)*x
        dxk=1+a/2+b/2
        dyk=a/2-b/2

        pk=1
        dxk=0
        dyk=0

        pkp1=(a/2-b/2)*y+(1+a/2+b/2)*x
        dxkp1=1+a/2+b/2
        dykp1=a/2-b/2

c       
        do 2000 k=2,n
c       
        pkm1=pk
        dxkm1=dxk
        dykm1=dyk
        
        pk=pkp1
        dxk=dxkp1
        dyk=dykp1

        alpha=(2*k+a+b-1)*(a**2-b**2)*y
     1     +(2*k+a+b-1)*(2*k+a+b-2)*(2*k+a+b)*x
        beta=2*(k+a-1)*(k+b-1)*(2*k+a+b) *y*y

        pkp1=(alpha*pk-beta*pkm1)/(2*k*(k+a+b)*(2*k+a+b-2))

        dxkp1=(alpha*dxk-beta*dxkm1)/(2*k*(k+a+b)*(2*k+a+b-2)) +
     1     (2*k+a+b-1)*(2*k+a+b-2)*(2*k+a+b) * pk 
     2     /(2*k*(k+a+b)*(2*k+a+b-2))

        dykp1=(alpha*dyk-beta*dykm1)/(2*k*(k+a+b)*(2*k+a+b-2)) +
     $     (2*k+a+b-1)*(a**2-b**2) *  pk /(2*k*(k+a+b)*(2*k+a+b-2)) -
     1     2*y *2*(k+a-1)*(k+b-1)*(2*k+a+b) * pkm1
     2     /(2*k*(k+a+b)*(2*k+a+b-2))

        pols(k+1)=pkp1
        dersx(k+1)=dxkp1
        dersy(k+1)=dykp1

 2000   continue
c       
        return
        end
       
        
        





