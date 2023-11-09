      subroutine bigrid(depth)
c
c --- set loop bounds for irregular basin in c-grid configuration
c --- q,u,v,p are vorticity, u-velocity, v-velocity, and mass points, resp.
c --- 'depth' = basin depth array, zero values indicate land
c
c --- this version works for both cyclic and noncyclic domains.
c --- land barrier at i=idm and/or j=jdm signals closed-basin conditions
c
      use hycom_dimen
      implicit none
      real, intent (in) :: depth(idm,jdm)
      integer ia,ib,ja,jb,nfill,nzero,jsec,jfrst,jlast,nwet,npts
c
      jchunk=23			! optimal for 8 threads
      nwet=0
      do 17 j=1,jdm
      do 17 i=1,idm
      if (depth(i,j).gt.0.) nwet=nwet+1
      ip(i,j)=0
      iq(i,j)=0
      iu(i,j)=0
 17   iv(i,j)=0
c
c --- mass points are defined where water depth is greater than zero
      do 2 j=1,jdm
      do 2 i=1,idm
      if (depth(i,j).gt.0.) ip(i,j)=1
  2   continue
c
c --- u,v points are located halfway between any 2 adjoining mass points
      do 3 j=1,jdm
      do 3 i=1,idm
      ia=mod(i-2+idm,idm)+1
      if (ip(ia,j).gt.0.and.ip(i,j).gt.0) iu(i,j)=1
  3   continue
      do 4 j=1,jdm
      ja=mod(j-2+jdm,jdm)+1
      do 4 i=1,idm
      if (ip(i,ja).gt.0.and.ip(i,j).gt.0) iv(i,j)=1
  4   continue
c
c --- 'interior' q points require water on all 4 sides.
      do 5 j=1,jdm
      ja=mod(j-2+jdm,jdm)+1
      do 5 i=1,idm
      ia=mod(i-2+idm,idm)+1
      if (min(ip(i,j),ip(ia,j),ip(i,ja),ip(ia,ja)).gt.0) iq(i,j)=1
  5   continue
c
c --- 'promontory' q points require water on 3 (or at least 2 diametrically 
c --- opposed) sides
      do 10 j=1,jdm
      ja=mod(j-2+jdm,jdm)+1
      do 10 i=1,idm
      ia=mod(i-2+idm,idm)+1
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
      npts=0
      do j=1,jdm
       do i=1,idm
        if (ip(i,j).eq.1) npts=npts+1
       end do
      end do
c     write (*,*) 'wet points in orig. depth array:',nwet
c     write (*,*) 'final number of grid points:    ',npts
c
c --- check load balance
c     call thrsiz(depth)
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
      use hycom_dimen
      implicit none
c
      integer ipt(idm,jdm),if(jdm,ms),il(jdm,ms),is(jdm)
      do 1 j=1,jdm
      is(j)=0
      do 4 k=1,ms
      if(j,k)=0
 4    il(j,k)=0
      i=1
      k=1
 3    if (ipt(i,j).ne.0) go to 2
      i=i+1
      if (i.le.idm) go to 3
      go to 1
 2    if (k.gt.ms) then
      write (*,'('' error in indxi - ms too small at i,j ='',2i5)') i,j
      write (*,'('' j-th line of ipt array:'',/(7(1x,10i1)))')
     .   (ipt(l,j),l=1,idm)
      stop '(indxi)'
      end if
      if(j,k)=i
 6    i=i+1
      if (i.le.idm) go to 5
      il(j,k)=idm
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
      subroutine indxj(jpt,jf,jl,js)
c
c --- input array jpt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays jf, jl, js  where
c --- jf(i,k) gives column index of first point in row i for k-th section
c --- jl(i,k) gives column index of last point
c --- js(i) gives number of sections in row i (maximum: ms)
c
      use hycom_dimen
      implicit none
c
      integer jpt(idm,jdm),jf(idm,ms),jl(idm,ms),js(idm)
      do 1 i=1,idm
      js(i)=0
      do 4 k=1,ms
      jf(i,k)=0
 4    jl(i,k)=0
      j=1
      k=1
 3    if (jpt(i,j).ne.0) go to 2
      j=j+1
      if (j.le.jdm) go to 3
      go to 1
 2    if (k.gt.ms) then
      write (*,'('' error in indxj - ms too small at i,j ='',2i5)') i,j
      write (*,'('' i-th line of jpt array:'',/(7(1x,10i1)))')
     .   (jpt(i,l),l=1,jdm)
      stop '(indxj)'
      end if
      jf(i,k)=j
 6    j=j+1
      if (j.le.jdm) go to 5
      jl(i,k)=jdm
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
