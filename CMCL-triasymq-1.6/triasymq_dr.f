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
c
c       ... retrieve D_3 symmetric quadrature rule for triangles
c
c       Input:
c
c       n - the degree of the quadrature (must not exceed 50)
c       vert1, vert2, vert3 - the vertices of the triangle
c
c       Output:
c
c       fort.11 - quadrature nodes and weights
c       fort.12 - quadrature nodes in gnuplot compatible format
c       fort.14 - triangle as constructed in gnuplot compatible format
c       
c       in gnuplot:
c
c       plot 'fort.14' w l, 'fort.12'
c
c
        implicit real *8 (a-h,o-z)
        dimension rnodes(2,2000)
        dimension weights(2000)
        dimension vert1(2),vert2(2),vert3(2)
        dimension rints(2000),z(3),pols(100 000)
        dimension work(10 000 000)
c
c       SET ALL PARAMETERS
c       
        call prini(6,13)
c
        PRINT *, 'ENTER n (1..50)'
        READ *, n
        call prinf('n=*',n,1)
c
c
        itype=1
c
        if( itype .eq. 1 ) then
c
c       ... standard equilateral triangle
c        
c        (-1,-1/sqrt(3)), (1,-1/sqrt(3)), (0,2/sqrt(3))
c
        vert1(1)=-1
        vert1(2)=-1/sqrt(3.0d0)
        vert2(1)=+1
        vert2(2)=-1/sqrt(3.0d0)
        vert3(1)=0
        vert3(2)=2/sqrt(3.0d0)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... standard simplex
c        
        vert1(1)=0
        vert1(2)=0
        vert2(1)=1
        vert2(2)=0
        vert3(1)=0
        vert3(2)=1
c
        endif
c
c
c
c       ... retrieve D_3 symmetric quadrature rule
c
        call triasymq(n,vert1,vert2,vert3,rnodes,
     1     weights,numnodes)
c
        call prinf('nummodes=*',numnodes,1)
        call prin2('rnodes=*',rnodes,2*numnodes)       
        call prin2('weights=*',weights,numnodes)       
c
        d=0
        do i=1,numnodes
        d=d+weights(i)
        enddo
c
        call prin2('sum of weights=*',d,1)
c        
c       ... plot the triangle
c
        write(14,*) vert1(1),vert1(2)
        write(14,*) vert2(1),vert2(2)
        write(14,*) vert3(1),vert3(2)
        write(14,*) vert1(1),vert1(2)
c
c       ... dump the nodes into a file
c       
        write(11,*) numnodes
        do i=1,numnodes
        write(11,1000) rnodes(1,i),rnodes(2,i),weights(i)
        write(12,*) rnodes(1,i),rnodes(2,i)
        enddo
 1010   format(i2)
 1000   format(4(e21.15,2x))
c

c
        mmax=n
c
c       construct the matrix of values of the orthogonal polynomials
c       at the user-provided nodes        
c
        npols=(mmax+1)*(mmax+2)/2
c
c
        do 1200 j=1,npols
        rints(j)=0
 1200 continue
c
c
        do 2400 i=1,numnodes
c
        z(1)=rnodes(1,i)
        z(2)=rnodes(2,i)
c       
        call ortho2eva(mmax,z,pols,work)
c
        do 2200 j=1,npols
c
        rints(j)=rints(j)+weights(i)*pols(j)

 2200 continue

 2400 continue
c
c
        call prin2('rints=*',rints,npols)
c
        area=sqrt(3.0d0)
c
        d=0
        d=(rints(1)-sqrt(area))**2
        do i=2,npols
        d=d+rints(i)**2
        enddo
c
        d=sqrt(d)/npols
c
        call prin2('error=*',d,1)
c

        stop
        end
