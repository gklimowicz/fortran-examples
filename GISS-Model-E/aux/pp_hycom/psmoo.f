      subroutine psmo1(alist,pbot)
c
c --- ragged-boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where ip > 0
c
c --- psmo1 is specially set up for interface smoothing.
c --- it only alters -alist- values that don't coincide with -pbot-.
c
      use hycom_dimen

      implicit none
      real, intent(INOUT) :: alist(idm,jdm)
      real, intent(IN)    :: pbot(idm,jdm)
      integer ia,ib,ja,jb
      real blist(idm,jdm),flxlo,flxhi
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ia)
      do 1 j=1,jdm
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
      ia=mod(i-2+idm,idm)+1
      if (ip(ia,j).gt.0) then
        flxhi= .25*(pbot(i ,j)-alist(i ,j))
        flxlo=-.25*(pbot(ia,j)-alist(ia,j))
        blist(i,j)=min(flxhi,max(flxlo,wgt*(alist(ia,j)-alist(i,j))))
      else
        blist(i,j)=0.
      end if
 1    continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ib)
      do 2 j=1,jdm
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
      ib=mod(i,idm)+1
      if (ip(ib,j).eq.0) blist(ib,j)=0.
      alist(i,j)=alist(i,j)-(blist(ib,j)-blist(i,j))
 2    continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ja)
      do 3 j=1,jdm
      do 3 l=1,isp(j)
      do 3 i=ifp(j,l),ilp(j,l)
      ja=mod(j-2+jdm,jdm)+1
      if (ip(i,ja).gt.0) then
        flxhi= .25*(pbot(i,j )-alist(i,j ))
        flxlo=-.25*(pbot(i,ja)-alist(i,ja))
        blist(i,j)=min(flxhi,max(flxlo,wgt*(alist(i,ja)-alist(i,j))))
      else
        blist(i,j)=0.
      end if
 3    continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(jb)
      do 4 j=1,jdm
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
      jb=mod(j,jdm)+1
      if (ip(i,jb).eq.0) blist(i,jb)=0.
      alist(i,j)=alist(i,j)-(blist(i,jb)-blist(i,j))
 4    continue
c$OMP END PARALLEL DO
c
      return
      end
c
c
      subroutine psmoo(alist)
c
c --- ragged-boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where ip > 0
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      use hycom_dimen

      implicit none
      real, intent(INOUT) :: alist(idm,jdm)
      integer ia,ib,ja,jb
      real blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb)
      do 1 j=1,jdm
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
      ja=mod(j-2+jdm,jdm)+1
      if (ip(i,ja).eq.0) ja=j
      jb=mod(j     ,jdm)+1
      if (ip(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib)
      do 2 j=1,jdm
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(idm,i+1)
      if (ip(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
      return
      end
c
c
      subroutine usmoo(alist)
c
c --- ragged-boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where iu > 0
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      use hycom_dimen

      implicit none
      real, intent(INOUT) :: alist(idm,jdm)
      integer ia,ib,ja,jb
      real blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb)
      do 1 j=1,jdm
      do 1 l=1,isu(j)
      do 1 i=ifu(j,l),ilu(j,l)
      ja=mod(j-2+jdm,jdm)+1
      if (iu(i,ja).eq.0) ja=j
      jb=mod(j     ,jdm)+1
      if (iu(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib)
      do 2 j=1,jdm
      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
      ia=max( 1,i-1)
      if (iu(ia,j).eq.0) ia=i
      ib=min(idm,i+1)
      if (iu(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
      return
      end
c
c
      subroutine vsmoo(alist)
c
c --- ragged-boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at poins where iv > 0
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      use hycom_dimen

      implicit none
      real, intent(INOUT) :: alist(idm,jdm)
      integer ia,ib,ja,jb
      real blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb)
      do 1 j=1,jdm
      do 1 l=1,isv(j)
      do 1 i=ifv(j,l),ilv(j,l)
      ja=mod(j-2+jdm,jdm)+1
      if (iv(i,ja).eq.0) ja=j
      jb=mod(j     ,jdm)+1
      if (iv(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib)
      do 2 j=1,jdm
      do 2 l=1,isv(j)
      do 2 i=ifv(j,l),ilv(j,l)
      ia=max( 1,i-1)
      if (iv(ia,j).eq.0) ia=i
      ib=min(idm,i+1)
      if (iv(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
      return
      end
c
c
      subroutine qsmoo(alist)
c
c --- ragged-boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at poins where iq > 0
c --- this routine is set up to smooth data carried at -q- points
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      use hycom_dimen

      implicit none
      real, intent(INOUT) :: alist(idm,jdm)
      integer ia,ib,ja,jb
      real blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb)
      do 1 j=1,jdm
      do 1 l=1,isq(j)
      do 1 i=ifq(j,l),ilq(j,l)
      ja=mod(j-2+jdm,jdm)+1
      if (iq(i,ja).eq.0) ja=j
      jb=mod(j     ,jdm)+1
      if (iq(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib)
      do 2 j=1,jdm
      do 2 l=1,isq(j)
      do 2 i=ifq(j,l),ilq(j,l)
      ia=max( 1,i-1)
      if (iq(ia,j).eq.0) ia=i
      ib=min(idm,i+1)
      if (iq(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
      return
      end
c
c
      subroutine psmoo2(alist)
c
c --- ragged-boundary version of 'heavy' 1-2-2-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where ip > 0
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      use hycom_dimen

      implicit none
      real,intent(INOUT) :: alist(idm,jdm)
      integer ia,ib,ja,jb
      real blist(idm,jdm)
c
c$OMP PARALLEL DO PRIVATE(ja,jb)
      do 1 j=1,jdm
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
      blist(i,j)=alist(i,j)
      ja=mod(j-2+jdm,jdm)+1
      if (ip(i,ja).eq.0) ja=j
      jb=mod(j     ,jdm)+1
      if (ip(i,jb).eq.0) jb=j
      blist(i,j)=.5*alist(i,j)+.25*(alist(i,ja)+alist(i,jb))
 1    continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib)
      do 2 j=1,jdm
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(idm,i+1)
      if (ip(ib,j).eq.0) ib=i
      alist(i,j)=.5*blist(i,j)+.25*(blist(ia,j)+blist(ib,j))
 2    continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ja,jb)
      do 3 j=1,jdm
      do 3 l=1,isp(j)
      do 3 i=ifp(j,l),ilp(j,l)
      blist(i,j)=alist(i,j)
      ja=mod(j-2+jdm,jdm)+1
      if (ip(i,ja).eq.0) ja=j
      jb=mod(j     ,jdm)+1
      if (ip(i,jb).eq.0) jb=j
      blist(i,j)=.5*(alist(i,ja)+alist(i,jb))
 3    continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib)
      do 4 j=1,jdm
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(idm,i+1)
      if (ip(ib,j).eq.0) ib=i
      alist(i,j)=.5*(blist(ia,j)+blist(ib,j))
 4    continue
c$OMP END PARALLEL DO
      return
      end

