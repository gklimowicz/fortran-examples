#include "hycom_mpi_hacks.h"
      subroutine bigrid(depth)
c
c --- set loop bounds for irregular basin in c-grid configuration
c --- q,u,v,p are vorticity, u-velocity, v-velocity, and mass points, resp.
c --- 'depth' = basin depth array, zero values indicate land
c
c --- this version works for both cyclic and noncyclic domains.
c --- land barrier at i=ii and/or j=jj signals closed-basin conditions
c
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT, pack_data
      USE HYCOM_DIM
      implicit none
      integer i,j,k,l,ia,ib,ja,jb
c
      real depth(idm,jdm)
      integer nfill,nzero,jsec,jfrst,jlast
      integer, allocatable :: ip_glob(:,:)
      character fmt*12,char2*2
      data fmt/'(i4,1x,75i1)'/
c
      do 17 j=J_0H,J_1H !1,jj
      do 17 i=1,ii
      ip(i,j)=0
      iq(i,j)=0
      iu(i,j)=0
 17   iv(i,j)=0
c
      go to 888    ! no change in topo after weights are done
c --- fill single-width inlets
 16   nfill=0
      do 15 j=J_0,J_1 !1,jj
      ! for now depth() is global array, so leave ja,jb as they are
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 15 i=1,ii
      ia=mod(i-2+ii,ii)+1
      ib=mod(i     ,ii)+1
      nzero=0
      if (depth(i,j).gt.0.) then
        if (depth(ia,j).le.0.) nzero=nzero+1
        if (depth(ib,j).le.0.) nzero=nzero+1
        if (depth(i,ja).le.0.) nzero=nzero+1
        if (depth(i,jb).le.0.) nzero=nzero+1
        if (nzero.ge.3) then
          write (*,'(a,i4,a,i4,a)') ' depth(',i,',',j,') set to zero'
          stop 'pre-process depth'
          depth(i,j)=0.
          nfill=nfill+1
        end if
      end if
 15   continue
      if (nfill.gt.0) go to 16
c
 888  continue
c --- mass points are defined where water depth is greater than zero
      do 2 j=J_0H,J_1H !1,jj
      do 2 i=1,ii
      if (depth(i,mod(j-1+jj,jj)+1).gt.0.) ip(i,j)=1
  2   continue
      !call halo_update(ogrid, ip)
c
c --- u,v points are located halfway between any 2 adjoining mass points
      !do 3 j=J_0,J_1 !1,jj  mkb
      do 3 j=J_0H,J_1H !1,jj
      do 3 i=1,ii
      ia=mod(i-2+ii,ii)+1
      ! if (ip(ia,j).gt.0.and.ip(i,j).gt.0) iu(i,j)=1 mkb
      if (depth(ia,mod(j-1+jj,jj)+1).gt.0.and.
     &    depth(i,mod(j-1+jj,jj)+1).gt.0) iu(i,j)=1
  3   continue
!     do 4 j=J_0,J_1 !1,jj   mkb
      do 4 j=J_0H,J_1H !1,jj
      ja=mod(j-2+jj,jj)+1
      do 4 i=1,ii
      !if (ip(i,ja).gt.0.and.ip(i,j).gt.0) iv(i,j)=1 mkb
      if (depth(i,mod(ja-1+jj,jj)+1).gt.0.and.
     &    depth(i,mod(j-1+jj,jj)+1).gt.0) iv(i,j)=1
  4   continue
c
c --- 'interior' q points require water on all 4 sides.
      do 5 j=J_0,J_1 !1,jj
      !ja=mod(j-2+jj,jj)+1
      !ja = j-1
      ja = PERIODIC_INDEX(j-1, jj)
      do 5 i=1,ii
      ia=mod(i-2+ii,ii)+1
      if (min(ip(i,j),ip(ia,j),ip(i,ja),ip(ia,ja)).gt.0) iq(i,j)=1
  5   continue
c
c --- 'promontory' q points require water on 3 (or at least 2 diametrically
c --- opposed) sides
      do 10 j=J_0,J_1 !1,jj
      !ja=mod(j-2+jj,jj)+1
      !ja = j-1
      ja = PERIODIC_INDEX(j-1, jj)
      do 10 i=1,ii
      ia=mod(i-2+ii,ii)+1
      if ((ip(i ,j).gt.0.and.ip(ia,ja).gt.0).or.
     .    (ip(ia,j).gt.0.and.ip(i ,ja).gt.0)) iq(i,j)=1
 10   continue
c
c --- determine loop bounds for vorticity points, including interior and
c --- promontory points
      call indxi(iq,ifq,ilq,isq)
      call indxj(iq,jfq,jlq,jsq)
c
c --- determine loop indices for mass and velocity points
      call indxi(ip,ifp,ilp,isp)
      call indxj(ip,jfp,jlp,jsp)
      call indxi(iu,ifu,ilu,isu)
      call indxj(iu,jfu,jlu,jsu)
      call indxi(iv,ifv,ilv,isv)
      call indxj(iv,jfv,jlv,jsv)
c
c --- write out  -ip-  array
c --- data are written in strips 75 points wide
      if (AM_I_ROOT()) then
        allocate( ip_glob(idm,jdm) )
      else
        allocate( ip_glob(1,1) )
      endif
      call pack_data(ogrid, ip, ip_glob)

      if (AM_I_ROOT()) then
        jsec=(jj-1)/75
        do 9 jfrst=0,75*jsec,75
          jlast=min(jj,jfrst+75)
          write (char2,'(i2)') jlast-jfrst
          fmt(8:9)=char2
          write (*,'(''ip array, cols'',i5,'' --'',i5)') jfrst+1,jlast
          write (*,fmt) (i,(10*ip_glob(i,j),j=jfrst+1,jlast),i=1,ii)
 9    continue
      endif
      deallocate( ip_glob )
c
      return
      end
c
c
      subroutine indxi(ipt,if,il,is)
c
c --- input array ipt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays if, il, is  where
c --- if(j,k) gives row index of first point in column j for k-th section
c --- il(j,k) gives row index of last point
c --- is(j) gives number of sections in column j (maximum: ms)
c
      USE HYCOM_DIM
      implicit none
      integer i,j,k,l
c
      integer ipt(idm,J_0H:J_1H),if(J_0H:J_1H,ms),il(J_0H:J_1H,ms),
     &     is(J_0H:J_1H)
      do 1 j=J_0,J_1 !1,jj
      is(j)=0
      do 4 k=1,ms
      if(j,k)=0
 4    il(j,k)=0
      i=1
      k=1
 3    if (ipt(i,j).ne.0) go to 2
      i=i+1
      if (i.le.ii) go to 3
      go to 1
 2    if (k.gt.ms) then
      write(0,*) "k,ms", k,ms
      write (*,'('' error in indxi - ms too small at i,j ='',2i5)') i,j
      write (*,'('' j-th line of ipt array:'',/(7(1x,10i1)))')
     .   (ipt(l,j),l=1,ii)
      stop '(indxi)'
      end if
      if(j,k)=i
 6    i=i+1
      if (i.le.ii) go to 5
      il(j,k)=ii
      is(j)=k
      go to 1
 5    if (ipt(i,j).ne.0) go to 6
      il(j,k)=i-1
      is(j)=k
      k=k+1
      go to 3
 1    continue
      return
      end
c
c
      subroutine indxj(jpt_loc,jf,jl,js)
c
c --- input array jpt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays jf, jl, js  where
c --- jf(i,k) gives column index of first point in row i for k-th section
c --- jl(i,k) gives column index of last point
c --- js(i) gives number of sections in row i (maximum: ms)
c
      USE DOMAIN_DECOMP_1D, only : pack_data, broadcast
      USE HYCOM_DIM
      implicit none
      integer i,j,k,l
c
      integer jpt_loc(idm,J_0H:J_1H),jf(idm,ms),jl(idm,ms),js(idm)
      integer jpt(idm,JDM)

      call pack_data(ogrid, jpt_loc, jpt)
      call broadcast(ogrid, jpt)
      do 1 i=1,ii
      js(i)=0
      do 4 k=1,ms
      jf(i,k)=0
 4    jl(i,k)=0
      j=1
      k=1
 3    if (jpt(i,j).ne.0) go to 2
      j=j+1
      if (j.le.jj) go to 3
      go to 1
 2    if (k.gt.ms) then
      write (*,'('' error in indxj - ms too small at i,j ='',2i5)') i,j
      write (*,'('' i-th line of jpt array:'',/(7(1x,10i1)))')
     .   (jpt(i,l),l=1,jj)
      stop '(indxj)'
      end if
      jf(i,k)=j
 6    j=j+1
      if (j.le.jj) go to 5
      jl(i,k)=jj
      js(i)=k
      go to 1
 5    if (jpt(i,j).ne.0) go to 6
      jl(i,k)=j-1
      js(i)=k
      k=k+1
      go to 3
 1    continue
      return
      end
c
c
c> Revision history:
c>
c> May  2000 - routine generalized to accomodate cyclic & noncyclic b.c.
c> Mar. 2006 - added bering strait exchange logic
