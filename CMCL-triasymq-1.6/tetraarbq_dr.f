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
c       
c       Quadratures for smooth functions on a tetrahedron
c
c       All quadratures are accurate to 15 digits
c       All weights are positive and inside the standard tetrahedron
c
c       (-1,-1/Sqrt[3],-1/Sqrt[6]), (0,2/Sqrt[3],-1/Sqrt[6]), 
c       (1,-1/Sqrt[3],-1/Sqrt[6]),  and (0,0,3/Sqrt[6])
c
c       Input:
c
c       n - the degree of the quadrature (must not exceed 9)
c
c       Output:
c
c       fort.11 - quadrature nodes and weights
c       fort.12 - quadrature nodes in gnuplot compatible format
c       fort.14 - standard tetrahedron as constructed 
c                   in gnuplot compatible format
c       
c       in gnuplot:
c
c       splot 'fort.14' w l, 'fort.12'
c
c
c
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,1000),weights(1000)
        dimension rints(1000),z(3),pols(1000)
        dimension work(100000),vert(3,4)
c
c
        call prini(6,13)
c
c       SET ALL PARAMETERS
c
        PRINT *, 'ENTER mmax (1..15)'
        READ *, mmax
c
        call prinf('mmax=*',mmax,1)
c
c

c
c       ... retrieve the quadrature rule 
c
        call tetraarbq(mmax,rnodes,weights,numnodes)
c
        call prinf('nummodes=*',numnodes,1)
        call prin2('rnodes=*',rnodes,3*numnodes)       
        call prin2('weights=*',weights,numnodes)       
c
        d=0
        do i=1,numnodes
        d=d+weights(i)
        enddo
c
        call prin2('sum of weights=*',d,1)
c        
c       ... plot the standard tetrahedron
c

        vert(1,1)=-1
        vert(2,1)=-1/sqrt(3.0d0)
        vert(3,1)=-1/sqrt(6.0d0)
c
        vert(1,2)=0
        vert(2,2)=2/sqrt(3.0d0)
        vert(3,2)=-1/sqrt(6.0d0)
c 
        vert(1,3)=1
        vert(2,3)=-1/sqrt(3.0d0)
        vert(3,3)=-1/sqrt(6.0d0)
c
        vert(1,4)=0
        vert(2,4)=0
        vert(3,4)=3/sqrt(6.0d0)
c

c
        write(14,*) vert(1,1),vert(2,1),vert(3,1)
        write(14,*) vert(1,2),vert(2,2),vert(3,2)
        write(14,*) vert(1,3),vert(2,3),vert(3,3)
        write(14,*) vert(1,1),vert(2,1),vert(3,1)
        write(14,*)
c
        write(14,*) vert(1,1),vert(2,1),vert(3,1)
        write(14,*) vert(1,4),vert(2,4),vert(3,4)
        write(14,*)
        write(14,*) vert(1,2),vert(2,2),vert(3,2)
        write(14,*) vert(1,4),vert(2,4),vert(3,4)
        write(14,*)
        write(14,*) vert(1,3),vert(2,3),vert(3,3)
        write(14,*) vert(1,4),vert(2,4),vert(3,4)
        write(14,*)
c
c
c       ... dump the nodes into a file
c       
        write(11,*) numnodes
        do i=1,numnodes
        write(11,1020) rnodes(1,i),rnodes(2,i),rnodes(3,i),weights(i)
        write(12,*) rnodes(1,i),rnodes(2,i),rnodes(3,i)
        enddo
 1010   format(i2)
 1020   format(4(e21.15,2x))
c

c
c       construct the matrix of values of the orthogonal polynomials
c       at the user-provided nodes        
c
        npols=(mmax+1)*(mmax+2)*(mmax+3)/6
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
        z(3)=rnodes(3,i)
c       
        call symeval(itype,mmax,z(1),z(2),z(3),pols,dersx,dersy,dersz)
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
        vol = 2.0d0 * sqrt(2.0d0) / 3.0d0
c
        d=0
        d=(rints(1)-sqrt(vol))**2
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
c
c
c
c
c
        subroutine symeval(itype,mmax,x,y,z,polsout,dersx,dersy,dersz)
        implicit real *8 (a-h,o-z)
        dimension z1(3),polsout(1),dersx(1),dersy(1),dersz(1)
c
        dimension polsout1(10 000),
     $       dersx1(10 000),dersy1(10 000),dersz1(10 000)
        dimension work(100 000)
c
        dimension v(4 000 000), v0(4 000 000), rnorms(10000)
        data ifinit/0/
c
        save
c
c        evaluate all polynomials 
c
        z1(1)=x
        z1(2)=y
        z1(3)=z
c
        npols7=(mmax+1)*(mmax+2)*(mmax+3)/6
c
        scale=1
c
        if(itype .eq. 1)
     1    call ortho3eva4
     2       (mmax,z1,polsout1,dersx1,dersy1,dersz1,work)
c
        if(itype .eq. 0)
     1    call ortho3eva(mmax,z1,polsout1,work)
c   
        if(itype .eq. 0) goto 4400
c  
        do 4200 i=1,npols7
c
        polsout(i)=polsout1(i) *scale
        dersx(i)=dersx1(i) *scale
        dersy(i)=dersy1(i) *scale
        dersz(i)=dersz1(i) *scale
c
 4200 continue
      return
c
        return
c
 4400 continue
c
        do 4600 i=1,npols7
c
        polsout(i)=polsout1(i) *scale
c
 4600 continue
      return
c
      return
c
c
c
        entry symevalget(mmax,ncolsout)
c
        npolsout=(mmax+1)*(mmax+2)*(mmax+3)/6
        npols7=npolsout
c
        return
        end        
c
c
c
c
c
