cc Copyright (C) 2009-2010: Zydrunas Gimbutas and Hong Xiao
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       Quadratures for smooth functions on tetrahedra
c       
c       Tensor product rules
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine tetragauc(n,vert1,vert2,vert3,vert4,
     1     rnodes,weights,numnodes,w)
        implicit real *8 (a-h,o-z)
c
c       this subroutine constructs a tensor product Gaussian
c       quadrature formula on the tetrahedron in the space
c       specified by its vertices. Actually, this 
c       is simply a memory manager for the subroutine tetragau1,
c       which does all the work.
c
c                              input parameters:
c
c  n - the number of nodes in each direction. note that the 
c       total number of nodes to be created is n**3, and the 
c       algebraic order of the quadrature is 2*n-1
c  vert1,vert2,vert3,vert4 - the vertices of the tetrahedron on which the 
c       quadrature rule is to be constructed
c
c                              output parameters:
c
c  rnodes - the array of n**3 nodes in the plane (all inside the
c       user-supplied tetrahedron)
c  weights - the quadrature weights corresponding to the nodes rnodes
c
c                              work arrays:
c
c  w - must be at least 2*n+2 real *8 locations for tetragau0
c
c  w - must be at least 2*n+2 + 3*n*n+2 real *8 locations for tetragau1
c
        dimension vert1(3),vert2(3),vert3(3),vert4(3)
        dimension rnodes(3,1),weights(1),w(1)
c
        itype=2
c
        if( itype .eq. 1 ) then
c
        its=1
        lts=n+1
c
        iws=its+lts
        lws=n+1
c        
        call tetragau0(n,vert1,vert2,vert3,vert4,
     1     rnodes,weights,numnodes,w(its),w(iws))        
c
        endif
c
c
        if( itype .eq. 2 ) then
c
        its=1
        lts=n+1
c
        iws=its+lts
        lws=n+1
c        
        irs2=iws+lws
        lrs2=2*n*n+1
c        
        iws2=irs2+lrs2
        lws2=n*n+1
c        
        call tetragau1(n,vert1,vert2,vert3,vert4,
     1     rnodes,weights,numnodes,w(its),w(iws),w(irs2),w(iws2))     
c
        endif

        return
        end
c
c
c
c
c
        subroutine tetragau0(n,vert1,vert2,vert3,vert4,
     1     rnodes,weights,numnodes,ts,ws)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3),vert4(3)
        dimension rnodes(3,n,n,n),weights(n,n,n),ts(n),ws(n)
        dimension ww(12),xy1(3),xy2(3),xy3(3),xy4(3),xy0(3)
        dimension endpts(2),b(10000),ts2(10000),ws2(10000)
c
        ifinit=1
        call legewhts(n,ts,ws,ifinit)
        do i=1,n
           ts(i)=(ts(i)+1)/2
           ws(i)=ws(i)/2
        enddo
c
        kind=5
        alpha=2.0d0
        beta=0.0d0
        kpts=0
        call gaussq(kind,n,alpha,beta,kpts,endpts,b,ts2,ws2)
        do i=1,n
           ts2(i)=(ts2(i)+1)/2
           ws2(i)=ws2(i)/2 
        enddo
        call prin2('ts2=*',ts2,n)
        call prin2('ws2=*',ws2,n)
c
        call tetraini(vert2,vert3,vert4,ww)
        call tetrafor(ww,vert1,xy1)
        call tetrafor(ww,vert2,xy2)
        call tetrafor(ww,vert3,xy3)
        call tetrafor(ww,vert4,xy4)
ccc        call prin2('xy1=*',xy1,3)
ccc        call prin2('xy2=*',xy2,3)
ccc        call prin2('xy3=*',xy3,3)
ccc        call prin2('xy4=*',xy4,3)
c
        do 1600 k=1,n
        do 1400 j=1,n
        do 1200 i=1,n
c
        rnodes(1,i,j,k)=ts(i)*(1-ts(j))*(1-ts(k))
        rnodes(2,i,j,k)=ts(j)*(1-ts(k))
        rnodes(3,i,j,k)=ts(k)
        weights(i,j,k)=ws(i)*ws(j)*ws(k)*(1-ts(j))*(1-ts(k))**2
c

        if( n.eq.1 ) then
        rnodes(1,i,j,k)=1/4.0d0
        rnodes(2,i,j,k)=1/4.0d0
        rnodes(3,i,j,k)=1/4.0d0
        weights(i,j,k)=1/6.0d0
        endif

c
        x=rnodes(1,i,j,k)
        y=rnodes(2,i,j,k)
        z=rnodes(3,i,j,k)
        xy0(1)=xy3(1)*x+xy4(1)*y+xy1(1)*z
        xy0(2)=         xy4(2)*y+xy1(2)*z
        xy0(3)=                  xy1(3)*z
c     
        scale=abs(xy3(1)*xy4(2)*xy1(3))
        call tetrabak(ww,xy0,rnodes(1,i,j,k))
        weights(i,j,k)=weights(i,j,k)*scale
c
ccc        write(18,*) rnodes(1,i,j,k),rnodes(2,i,j,k),rnodes(3,i,j,k)
c
 1200   continue
 1400   continue
 1600   continue
c
        numnodes=n**3
c
ccc        call prin2('rnodes=*',rnodes,3*numnodes)
ccc        call prin2('weights=*',weights,numnodes)
c
        return
        end
c
c
c
c
c
        subroutine tetragau1(n,vert1,vert2,vert3,vert4,
     1     rnodes,weights,numnodes,ts,ws,rs2,ws2)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3),vert4(3)
        dimension rnodes(3,1),weights(1),ts(n),ws(n)
        dimension ww(12),xy1(3),xy2(3),xy3(3),xy4(3),xy0(3)
        dimension v1(2),v2(2),v3(2),rs2(2,1),ws2(1)
c
        ifinit=1
        n1=n+1
        call legewhts(n1,ts,ws,ifinit)
        do i=1,n1
           ts(i)=(ts(i)+1)/2
           ws(i)=ws(i)/2
        enddo
c
c
        v1(1)=0
        v1(2)=0
        v2(1)=1
        v2(2)=0
        v3(1)=0
        v3(2)=1

        call triasymq(n*2-1,v1,v2,v3,rs2,ws2,n2)
c
        call tetraini(vert2,vert3,vert4,ww)
        call tetrafor(ww,vert1,xy1)
        call tetrafor(ww,vert2,xy2)
        call tetrafor(ww,vert3,xy3)
        call tetrafor(ww,vert4,xy4)
ccc        call prin2('xy1=*',xy1,3)
ccc        call prin2('xy2=*',xy2,3)
ccc        call prin2('xy3=*',xy3,3)
ccc        call prin2('xy4=*',xy4,3)
c
        kk=0
        do k=1,n1
        do i=1,n2
c
        kk=kk+1
        rnodes(1,kk)=rs2(1,i)*(1-ts(k))
        rnodes(2,kk)=rs2(2,i)*(1-ts(k))
        rnodes(3,kk)=ts(k)
        weights(kk)=ws2(i)*ws(k)*(1-ts(k))**2
c
        enddo
        enddo



        if( n.eq.1 ) then
        rnodes(1,1)=1/4.0d0
        rnodes(2,1)=1/4.0d0
        rnodes(3,1)=1/4.0d0
        weights(1)=1/6.0d0
        kk=1
        endif


        do 1600 k=1,kk
c
        x=rnodes(1,k)
        y=rnodes(2,k)
        z=rnodes(3,k)
        xy0(1)=xy3(1)*x+xy4(1)*y+xy1(1)*z
        xy0(2)=         xy4(2)*y+xy1(2)*z
        xy0(3)=                  xy1(3)*z
c     
        scale=abs(xy3(1)*xy4(2)*xy1(3))
        call tetrabak(ww,xy0,rnodes(1,k))
        weights(k)=weights(k)*scale
c
ccc        write(18,*) rnodes(1,k),rnodes(2,k),rnodes(3,k)
c
 1600   continue
c
        numnodes=kk
c
ccc        call prin2('rnodes=*',rnodes,3*numnodes)
ccc        call prin2('weights=*',weights,numnodes)
c
        return
        end
c
c
c
c
c
        subroutine tetraini(vert1,vert2,vert3,w)     
        implicit real *8 (a-h,o-z)
c
c     given a tetrahedron in R^3 with vertices (vert1,vert2,vert3,vert4)
c
c     this subroutine constructs an affine transformation putting the
c     triangle spoecified by the vertices T=(vert1,vert2,vert3) of the
c     user-specified tetrahedron on the x-y plane in such a manner that
c     one side of the triangle T is on the x-axis.
c
c     Note: the triangle (vert1,vert2,vert3) is the opposite to the
c     fouth vertex vert4 in R^3.
c
c     this subroutine also constructs the transformation inverse to the 
c     above one. The actual transformations are applied to particular
c     user-specified points by the entries tetrafor and tetrabak of this
c     subroutine
c                              input parameters:
c
c  vert1,vert2,vert3 - the vertices of the triangle to be put on x-y plane,
c                      (one of whose sides is to be put on the x axis)
c
c                              output parameters:
c
c  w - the real *8 array of length 12 containing the transformations;
c       it is to be used by the entries tetrafor, tetrabak (see below)
c

        dimension vert1(3),vert2(3),vert3(3),w(12)
c       
        call tetraortho3d(vert1,vert2,vert3,w(1),w(10)) 
        return
c
c
c
        entry tetrafor(w,zin,zout)
c
c       this entry applies to the user-specified point zin \in R^3
c       the first of the transformations constructed by the entry 
c       tetraini (see above), obtaining the point zout \in R^3.
c       
        call tetraprod3s(zin,w(1),w(10),zout)
        return
c
c
c
        entry tetrabak(w,zin,zout)
c
c       this entry applies to the user-specified point zin \in R^3
c       the second of the transformations constructed by the entry 
c       tetraini (see above), obtaining the point zout \in R^3.
c       
c
        call tetraprod3r(zin,w(1),w(10),zout)
c
        return
        end
c       
c
c
c
c       
        subroutine tetraprod3s(z,u,shift,y)     
        implicit real *8 (a-h,o-z)
        dimension u(3,3), shift(3), x(3), y(3), z(3)
        x(1)=z(1)+shift(1)
        x(2)=z(2)+shift(2)
        x(3)=z(3)+shift(3)
        y(1)=u(1,1)*x(1)+u(2,1)*x(2)+u(3,1)*x(3)
        y(2)=u(1,2)*x(1)+u(2,2)*x(2)+u(3,2)*x(3)
        y(3)=u(1,3)*x(1)+u(2,3)*x(2)+u(3,3)*x(3)
        return
        end
c
c
c
c
c
        subroutine tetraprod3r(z,u,shift,y)     
        implicit real *8 (a-h,o-z)
        dimension u(3,3), shift(3), x(3), y(3), z(3)
        x(1)=u(1,1)*z(1)+u(1,2)*z(2)+u(1,3)*z(3)
        x(2)=u(2,1)*z(1)+u(2,2)*z(2)+u(2,3)*z(3)
        x(3)=u(3,1)*z(1)+u(3,2)*z(2)+u(3,3)*z(3)
        y(1)=x(1)-shift(1)
        y(2)=x(2)-shift(2)
        y(3)=x(3)-shift(3)
        return
        end
c
c
c
c
c
        subroutine tetraortho3d(x,y,z,u,shift)     
        implicit real *8 (a-h,o-z)
        dimension x(3), y(3), z(3), shift(3), v(3,5), u(3,3)
        dimension rnorms(3)
        shift(1)=-x(1)
        shift(2)=-x(2)
        shift(3)=-x(3)
        u(1,1)=y(1)+shift(1)
        u(2,1)=y(2)+shift(2)
        u(3,1)=y(3)+shift(3)
        u(1,2)=z(1)+shift(1)
        u(2,2)=z(2)+shift(2)
        u(3,2)=z(3)+shift(3)
        eps=0
        ifpivot=0
        n=3
        m=2
        call tetragrmp2(u,n,m,rnorms,eps,ncols,ifpivot)
        call tetravprod3d(u(1,1), u(1,2), u(1,3))
        call tetranorm3d(u(1,3))
ccc      call prin2('inside orto3d, u=*', u, 3*2)
ccc      call prin2('inside orto3d, rnorms=*', rnorms, 2)
        return
ccc      call tetravprod3d(u(1,1), u(1,2), u(1,3))
ccc      call tetranorm3d(u(1,3))
ccc      call tetravprod3d(u(1,3), u(1,1), u(1,2))
ccc      call tetranorm3d(u(1,2))
ccc      call tetravprod3d(u(1,2), u(1,3), u(1,1))
ccc      call tetranorm3d(u(1,1))
        return
        end

c
c
c
c
c
        subroutine tetragrmp2(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        dimension b(n,m),rnorms(1)
c       
c       this subroutine applies a pivoted double gram-schmidt 
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is 
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c
c       NOTE: this version of a gram-schmidt procedure preserves
c       the orientation 
c
c                    input paramneters:
c
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy 
c  
c                     output parameters:
c
c  b - the matrix of gram-schmidt vectors of the matrix a. note 
c        that on exit from this subroutine only the first 
c        ncols columns of  b  are meaningful, the rest making no 
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c
c
c       . . . prepare the array of values of norms 
c             of columns
c
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)**2
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c
c       . . . conduct gram-schmidt iterations
c
        thresh=dtot*eps**2
        do 4000 i=1,m
c
c       find the pivot  
c
        if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c
c       put the column number ipivot in the i-th place
c
        write(*,*) 'pivot ', i, ipivot
        do 2600 j=1,n
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
c
c       preserve the orientation
c
        if( i.ne.ipivot ) b(j,ipivot)=-d
 2600 continue
c
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c
c       orthogonalize the i-th column to all preceeding ones
c
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c
        call tetrascap2(b(1,i),b(1,j),n,d)
c
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c
c       normalize the i-th column
c
        call tetrascap2(b(1,i),b(1,i),n,d)
c
c       if the norm of the remaining part of the matrix 
c       is sufficiently small - terminate the process
c
        if( (ifpivot .eq. 1) .and.
     1      (d .lt. thresh) ) return
        ncols=i
c
c
c       prevent division by zero
c
        if(d.ne.0) d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c
c        orthogonalize everything else to it 
c
        do 3200 j=i+1,m
c
        if(rnorms(j) .lt. thresh) goto 3200
c
        call tetrascap2(b(1,i),b(1,j),n,d)
c     
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue 
 3400 continue
 4000 continue

 4200 continue

        return
        end
c
c
c
        subroutine tetrascap2(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1)
c
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end

c
c
c
c
c
      subroutine tetravprod2d(x,y,z)     
      implicit real *8 (a-h,o-z)
c     
c     vector product in R^2
c
      dimension x(2), y(2)
      z= (x(1)*y(2)-y(2)*x(1))
      return
      end
c
c
c
c
c
      subroutine tetravprod3d(x,y,z)     
      implicit real *8 (a-h,o-z)
c     
c     vector product in R^3
c
      dimension x(3), y(3), z(3)
      z(1)= (x(2)*y(3)-y(2)*x(3))
      z(2)=-(x(1)*y(3)-y(1)*x(3))
      z(3)= (x(1)*y(2)-y(1)*x(2))
      return
      end
c
c
c
c
      subroutine tetranorm3d(z)     
      implicit real *8 (a-h,o-z)
c
c     normalize vector z in R^3
c
      dimension z(3)
      d = dsqrt(z(1)**2+z(2)**2+z(3)**2)
      z(1)=z(1)/d
      z(2)=z(2)/d
      z(3)=z(3)/d
      return
      end


c
c
c
      subroutine tetraprod3u(x,u,y)     
      implicit real *8 (a-h,o-z)
c
c        multiply matrices u' and x getting y
c
      dimension u(3,3), x(3), y(3)
      y(1)=u(1,1)*x(1)+u(2,1)*x(2)+u(3,1)*x(3)
      y(2)=u(1,2)*x(1)+u(2,2)*x(2)+u(3,2)*x(3)
      y(3)=u(1,3)*x(1)+u(2,3)*x(2)+u(3,3)*x(3)
      return
      end

c
c
c
c
c
      subroutine tetraplot(iw,vert1,vert2,vert3,vert4)     
      implicit real *8 (a-h,o-z)
      dimension vert1(3),vert2(3),vert3(3),vert4(3)
      write(iw,1000) vert1
      write(iw,1000) vert2
      write(iw,2000) 
      write(iw,2000) 
      write(iw,1000) vert2
      write(iw,1000) vert3
      write(iw,2000) 
      write(iw,2000) 
      write(iw,1000) vert3
      write(iw,1000) vert1
      write(iw,2000) 
      write(iw,2000) 
      write(iw,1000) vert1
      write(iw,1000) vert4
      write(iw,2000) 
      write(iw,2000) 
      write(iw,1000) vert2
      write(iw,1000) vert4
      write(iw,2000) 
      write(iw,2000) 
      write(iw,1000) vert3
      write(iw,1000) vert4
      write(iw,2000) 
      write(iw,2000) 
 1000 format(3(2x,e11.5))
 2000 format(80a1)
      return
      end
c
c
c
c
c
