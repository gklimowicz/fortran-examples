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
c       Orthogonal polynomials on the standard cube [-1,1]^3
c       
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains a set of subroutines for the handling 
c       of the evaluation of othogonal polynomials on the standard 
c       cube.  It contains 2 user-callable subroutines.
c
c   lege3eva - evaluates at the user-supplied point z
c       a collection of polynomials (of x,y,z) orthogonal on the
c       standard cube
c
c   lege3eva3 - evaluates at the user-supplied point z a collection
c       of polynomials (of x,y,z) orthogonal on the standard cube,
c       together with their derivatives with respect to x, y and z.


c       
c
      subroutine lege3eva(mmax,z,pols,w)
c
c     This subroutine evaluates at the user-supplied point z
c     a collection of polynomials (of x,y,z) orthogonal on the
c     standard cube
c
c                   Input parameters:
c
c  mmax - the maximum order to which the polynomials are to be evaluated;
c  z - the location in R^3 where the polynomials are to be evaluated;
c
c  w - the work array. must be at least 3*(mmax+10) long
c     
c                   Output parameters:
c
c  pols - the orthogonal polynomials evaluated at the point z 
c       ( (mmax+1)*(mmax+2)*(mmax+3)/6 of them things)
c
c       . . . allocate memory for the work arrays
c
      implicit real *8 (a-h,o-z)
      dimension z(2), pols(1), w(1)
      iw1=1
      iw2=iw1+mmax+10
      iw3=iw2+mmax+10
      iw4=iw3+mmax+10
      call lege3eva0(mmax,z,pols,w(iw1),w(iw2),w(iw3))
      return
      end
c
c
c
c
c
      subroutine lege3eva0(mmax,z,pols,f1,f2,f3)
c
c     evaluate the orthonormal polynomials on the cube
c
      implicit real *8 (a-h,o-z)
      dimension z(3),pols(1),f1(1),f2(1),f3(1)
ccc      data init/1/
ccc      save
c
c
      call llegepols1(z(1),mmax,f1)
      call llegepols1(z(2),mmax,f2)
      call llegepols1(z(3),mmax,f3)
c
c
      kk=0
      do 2400 m=0,mmax
         do 2200 n2=0,m
         do 2000 n1=0,n2
         kk=kk+1
c
c     ... evaluate the polynomial (m-n2, n2-n1, n1)
c
         pols(kk)=f1(m-n2+1)*f2(n2-n1+1)*f3(n1+1)
c
c     ... and normalize it
c
         scale=1
         scale=scale*dsqrt((1.0d0+n1+n1)/2)
         scale=scale*dsqrt((1.0d0+n2-n1+n2-n1)/2)
         scale=scale*dsqrt((1.0d0+m-n2+m-n2)/2)
         pols(kk)=pols(kk)*scale
 2000	 continue
 2200 continue
 2400 continue
      return
      end
c
c
c
c
c
      subroutine lege3eva3(mmax,z,pols,dersx,dersy,dersz,w)
c
c     This subroutine evaluates at the user-supplied point z
c     a collection of polynomials (of x,y,z) orthogonal on the
c     standard cube, together with their derivatives with 
c     respect to x, y and z.
c
c                   Input parameters:
c
c  mmax - the maximum order to which the polynomials are to be evaluated;
c  z - the location in R^3 where the polynomials are to be evaluated;
c
c  w - the work array. must be at least 6*(mmax+10) long
c     
c                   Output parameters:
c
c  pols - the orthogonal polynomials evaluated at the point z 
c       ( (mmax+1)*(mmax+2)*(mmax+3)/6 of them things)
c  dersx - the derivatives with respect to x of the polynomials 
c       returned in array pols
c  dersy - the derivatives with respect to y of the polynomials 
c       returned in array pols
c  dersz - the derivatives with respect to z of the polynomials 
c       returned in array pols
c
c
c       . . . allocate memory for the work arrays
c
      implicit real *8 (a-h,o-z)
      dimension z(2), pols(1), w(1)
      iw1=1
      iw2=iw1+mmax+10
      iw3=iw2+mmax+10
      iw4=iw3+mmax+10
      iw5=iw4+mmax+10
      iw6=iw5+mmax+10
      iw7=iw6+mmax+10
      call lege3eva30(mmax,z,pols,dersx,dersy,dersz,
     $     w(iw1),w(iw2),w(iw3),w(iw4),w(iw5),w(iw6))
      return
      end
c
c
c
c
c
      subroutine lege3eva30(mmax,z,pols,dersx,dersy,dersz,
     $     f1,f2,f3,f4,f5,f6)
c
c     evaluate the orthonormal polynomials on the cube, 
c     together with their derivatives with respect to x, y and z.
c
      implicit real *8 (a-h,o-z)
      dimension z(3), pols(1), dersx(1), dersy(1), dersz(1)
      dimension f1(1), f2(1)
      dimension f3(1), f4(1)
      dimension f5(1), f6(1)
ccc      data init/1/
ccc      save
c
c
c
      call llegepols2(z(1),mmax,f1,f4)
      call llegepols2(z(2),mmax,f2,f5)
      call llegepols2(z(3),mmax,f3,f6)
c

      kk=0
      do 2400 m=0,mmax
         do 2200 n2=0,m
         do 2000 n1=0,n2
         kk=kk+1
c
c     ... evaluate the polynomial (m-n2, n2-n1, n1)
c
         pols(kk)=f1(m-n2+1)*f2(n2-n1+1)*f3(n1+1)
c
c     ... and the derivatives with respect in x and y
c
         dersx(kk)=f4(m-n2+1)*f2(n2-n1+1)*f3(n1+1)
         dersy(kk)=f1(m-n2+1)*f5(n2-n1+1)*f3(n1+1)
         dersz(kk)=f1(m-n2+1)*f2(n2-n1+1)*f6(n1+1)
c
c     ... and normalize it
c
         scale=1
         scale=scale*dsqrt((1.0d0+n1+n1)/2)
         scale=scale*dsqrt((1.0d0+n2-n1+n2-n1)/2)
         scale=scale*dsqrt((1.0d0+m-n2+m-n2)/2)
         pols(kk)=pols(kk)*scale
         dersx(kk)=dersx(kk)*scale
         dersy(kk)=dersy(kk)*scale
         dersz(kk)=dersz(kk)*scale
 2000	 continue
 2200 continue
 2400 continue

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
        subroutine llegepols1(x,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(1)
c
c       evaluates a bunch of Legendre polynomials 
c       at the user-provided point x
c
c       ... if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=x
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        pols(1)=1
        pols(2)=x
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
        do 2000 k=1,n-1
c
        pkm1=pk
        pk=pkp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine llegepols2(x,n,pols,ders)
        implicit real *8 (a-h,o-z)
        dimension pols(1),ders(1)
c
c       evaluates a bunch of Legendre polynomials (together
c       with their derivatives) at the user-provided point x
c
c       ... if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        ders(1)=0
        if(n .eq. 0) return
c
        pols(2)=x
        ders(2)=1
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        pols(1)=1
        ders(1)=0
        pols(2)=x
        ders(2)=1
c
        pkm1=1
        pk=x
        dkm1=0
        dk=1
c
        pk=1
        pkp1=x
        dk=0
        dkp1=1
c
        do 2000 k=1,n-1
c
        pkm1=pk
        pk=pkp1
        dkm1=dk
        dk=dkp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        dkp1=( (2*k+1)*(x*dk+pk)-k*dkm1 )/(k+1)
        pols(k+2)=pkp1
        ders(k+2)=dkp1
 2000 continue
c
        return
        end







c
        subroutine lege3gauc(n,vert1,vert2,rnodes,
     1      weights,ifinit,w)
        implicit real *8 (a-h,o-z)
        dimension w(1),vert1(1),vert2(1),
     1       rnodes(1),weights(1)
c
c
c       this subroutine constructs a tensor product Gaussian
c       quadrature formula on the cube in the space 
c       specified by its vertices. 
c

        iw1=1
        iw2=iw1+n+1
        iw3=iw2+n+1
        iw4=iw3+n+1
        iw5=iw4+n+1
        iw5=iw5+n+1
        iw6=iw6+n+1

        call lege3gauc0(n,vert1,vert2,rnodes,
     1      weights,ifinit,w(iw1),w(iw2),w(iw3),w(iw4),w(iw5),w(iw6))

        return
        end




        subroutine lege3gauc0(n,vert1,vert2,rnodes,
     1      weights,ifinit,tsx,wsx,tsy,wsy,tsz,wsz)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3)
        dimension rnodes(3,n,n,n),weights(n,n,n)
        dimension tsx(n),wsx(n),tsy(n),wsy(n),tsz(n),wsz(n)
c
c
c       this subroutine constructs a tensor product Gaussian
c       quadrature formula on the square in the plane 
c       specified by its vertices. 
c

        if( ifinit.eq.1) then
           ifinit = 1 
           call legewhts(n,tsx,wsx,ifinit)
           call legewhts(n,tsy,wsy,ifinit)
           call legewhts(n,tsz,wsz,ifinit)
        endif

        u1=(vert2(1)-vert1(1))/2
        v1=(vert2(1)+vert1(1))/2
        u2=(vert2(2)-vert1(2))/2
        v2=(vert2(2)+vert1(2))/2
        u3=(vert2(3)-vert1(3))/2
        v3=(vert2(3)+vert1(3))/2

ccc        call prin2('vert1=*',vert1,2)
ccc        call prin2('vert2=*',vert2,2)

        do i=1,n
        do j=1,n
        do k=1,n
           rnodes(1,i,j,k)=u1*tsx(i)+v1
           rnodes(2,i,j,k)=u2*tsy(j)+v2
           rnodes(3,i,j,k)=u3*tsz(k)+v3
           weights(i,j,k)=abs(u1)*abs(u2)*abs(u3)*wsx(i)*wsy(j)*wsz(k)
        enddo
        enddo
        enddo
        
        return
        end


