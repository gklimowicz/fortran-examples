#include "hycom_mpi_hacks.h"
      subroutine reflux_pr(uflxo,vflxo,pold,
     .                     uflxn,vflxn,pnew,kold,knew)
c
c --- convert an idm x jdm array of mass fluxes in layers marked by
c --- interfaces -pold- to fluxes in isobaric layers
c
c --- input  variables: uflxo,vflxo,pold,kold,knew
c --- output variables: uflxn,vflxn,pnew
c
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : acurcy,onem,itest,jtest
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE,SOUTH
      implicit none
      integer i,j,k,l,ja
c
      integer kold,knew,ko
c
c-------------------------------------------------------------
      real uflxo(idm,J_0H:J_1H,kold),vflxo(idm,J_0H:J_1H,kold),
     .     pold(idm,J_0H:J_1H,kold+1)
c
      real uflxn(idm,J_0H:J_1H,knew),vflxn(idm,J_0H:J_1H,knew),
     .     pnew(idm,J_0H:J_1H,knew+1)
c
      real cloutu(idm),cloutv(idm),
     .     colinu(idm),colinv(idm),
     .     uold(idm,kdm),vold(idm,kdm),
     .     uinteg,vinteg,phi,plo,pa,pb,q,
     .     plft,prgt,pblft,pbrgt,delp,uvscal,pr1d

      logical abort,vrbos
      data abort/.false./
      data uvscal/1.e5/                        !  velocity x mesh size  --  SI
ccc   data uvscal/1.e9/                        !  velocity x mesh size  --  cgs
      pr1d(k)=5000.*float(k-1)/float(knew)* onem  ! isobaric depth levels

      vrbos(i,j)=i.eq.itest .and. j.eq.jtest
c
      do 1 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 1 l=1,isp(j)
c
      do 1 i=ifp(j,l),ilp(j,l)
      if (vrbos(i,j)) write (*,102) itest,jtest,
     .  ' (reflux_pr) old profile:    dpthu      u        dpthv      v',
     .  (k,.5*(pold(i,j,k+1)+pold(i-1,j,k+1))/onem,uflxo(i,j,k),
     .     .5*(pold(i,j,k+1)+pold(i,ja ,k+1))/onem,vflxo(i,j,k),
     .   k=1,kold)
 102  format (2i5,a/(i34,2(f10.2,es10.2)))
      pnew(i,j,     1)=pold(i,j,     1)
      pnew(i,j,knew+1)=pold(i,j,kold+1)
      do 2 k=2,knew
 2    pnew(i,j,k)=max(pnew(i,j,k-1),pr1d(k))
      do 3 k=knew,2,-1
 3    pnew(i,j,k)=min(pnew(i,j,k+1),pr1d(k))
 1    continue
c
       CALL HALO_UPDATE(ogrid,pold, FROM=SOUTH)
       CALL HALO_UPDATE(ogrid,pnew, FROM=SOUTH)
c
      do 21 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 24 k=1,kold
      do 24 i=1,ii
      uold(i,k)=uflxo(i,j,k)
 24   vold(i,k)=vflxo(i,j,k)
c
      do 22 l=1,isu(j)
c
c --- convert -uflx- to -u- and integrate over new depth intervals
c
      do 17 i=ifu(j,l),ilu(j,l)
 17   colinu(i)=0.
c
      do 18 k=kold,1,-1
      do 18 i=ifu(j,l),ilu(j,l)
      delp=min(depthu(i,j),.5*(pold(i,j,k+1)+pold(i-1,j,k+1)))-
     .     min(depthu(i,j),.5*(pold(i,j,k  )+pold(i-1,j,k  )))
      if (delp.gt.0.) then
        colinu(i)=colinu(i)+uold(i,k)
        uold(i,k)=uold(i,k)/delp
      else
        if (k.eq.1) then
          write (*,'(a,2i5)') 'reflux_pr error in loop 18 -- i,j =',i,j
          abort=.true.
        else
          uold(i,k-1)=uold(i,k-1)+uold(i,k)
          uold(i,k)=0.
        end if
      end if
 18   continue
c
      do 11 i=ifu(j,l),ilu(j,l)
 11   cloutu(i)=0.
c
      do 12 k=1,knew
      do 12 i=ifu(j,l),ilu(j,l)
      uinteg=0.
      phi=min(depthu(i,j),.5*(pnew(i,j,k+1)+pnew(i-1,j,k+1)))
      plo=min(depthu(i,j),.5*(pnew(i,j,k  )+pnew(i-1,j,k  )))
      pb=plo
      do 13 ko=1,kold
      q=min(depthu(i,j),.5*(pold(i,j,ko+1)+pold(i-1,j,ko+1)))
      if (q.le.plo) go to 13
      pa=pb
      pb=min(phi,q)
      uinteg=uinteg+uold(i,ko)*(pb-pa)
      if (pb.ge.phi) go to 26
 13   continue
 26   cloutu(i)=cloutu(i)+uinteg
 12   uflxn(i,j,k)=uinteg
c
      do 22 i=ifu(j,l),ilu(j,l)
      if (abs(cloutu(i)-colinu(i)).gt.acurcy*uvscal*depthu(i,j))
     .  write (*,100) i,j,'  reflux - bad u intgl.',colinu(i),
     .      cloutu(i),(cloutu(i)-colinu(i))/colinu(i)
 100  format (2i5,a,2es16.8,es9.1)
 22   continue
c
c --- convert -vflx- to -v- and integrate over new depth intervals
c
      do 23 l=1,isv(j)
c
      do 19 i=ifv(j,l),ilv(j,l)
 19   colinv(i)=0.
c
      do 20 k=kold,1,-1
      do 20 i=ifv(j,l),ilv(j,l)
      delp=min(depthv(i,j),.5*(pold(i,j,k+1)+pold(i,ja ,k+1)))-
     .     min(depthv(i,j),.5*(pold(i,j,k  )+pold(i,ja ,k  )))
      if (delp.gt.0.) then
        colinv(i)=colinv(i)+vold(i,k)
        vold(i,k)=vold(i,k)/delp
      else
        if (k.eq.1) then
          write (*,'(a,2i5)') 'reflux_pr error in loop 20 -- i,j =',i,j
          abort=.true.
        else
          vold(i,k-1)=vold(i,k-1)+vold(i,k)
          vold(i,k)=0.
        end if
      end if
 20   continue
c
      do 14 i=ifv(j,l),ilv(j,l)
 14   cloutv(i)=0.
c
      do 15 k=1,knew
      do 15 i=ifv(j,l),ilv(j,l)
      vinteg=0.
      phi=min(depthv(i,j),.5*(pnew(i,j,k+1)+pnew(i,ja ,k+1)))
      plo=min(depthv(i,j),.5*(pnew(i,j,k  )+pnew(i,ja ,k  )))
      pb=plo
      do 16 ko=1,kold
      q=min(depthv(i,j),.5*(pold(i,j,ko+1)+pold(i,ja ,ko+1)))
      if (q.le.plo) go to 16
      pa=pb
      pb=min(phi,q)
      vinteg=vinteg+vold(i,ko)*(pb-pa)
      if (pb.ge.phi) go to 27
 16   continue
 27   cloutv(i)=cloutv(i)+vinteg
 15   vflxn(i,j,k)=vinteg
c
      do 23 i=ifv(j,l),ilv(j,l)
      if (abs(cloutv(i)-colinv(i)).gt.acurcy*uvscal*depthv(i,j))
     .  write (*,100) i,j,'  reflux - bad v intgl.',colinv(i),
     .      cloutv(i),(cloutv(i)-colinv(i))/colinv(i)
 23   continue
c
      do 21 l=1,isp(j)
      do 21 i=ifp(j,l),ilp(j,l)
      if (vrbos(i,j)) write (*,102) itest,jtest,
     .  ' (reflux_pr) new profile:    dpthu      u        dpthv      v',
     .  (k,.5*(pnew(i,j,k+1)+pnew(i-1,j,k+1))/onem,uflxn(i,j,k),
     .     .5*(pnew(i,j,k+1)+pnew(i,ja ,k+1))/onem,vflxn(i,j,k),
     .   k=1,knew)
c
 21   continue
      if (abort) stop '(reflux_pr)'
      return
      end subroutine reflux_pr
